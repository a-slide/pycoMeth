# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Standard library imports
from collections import *
import csv
import datetime
from random import randint

# Third party imports
from tqdm import tqdm
import numpy as np

# Local imports
from pycoMeth.common import *
from pycoMeth import __version__ as package_version
from pycoMeth import __name__ as package_name

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MAIN CLASS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

class Freq_meth_calculate():

    def __init__ (self,
        input_fn:"str",
        fasta_index:"str"="",
        output_bed_fn:"str"="",
        output_tsv_fn:"str"="",
        min_depth:"int"=10,
        sample_id:"str"="",
        strand_specific:"bool"=False,
        min_llr:"float"=2,
        verbose:"bool"=False,
        quiet:"bool"=False):
        """
        Calculate methylation frequency at genomic CpG sites from the output of nanopolish call-methylation
        * input_fn
            Path to a nanopolish call_methylation tsv output file
        * fasta_index
            fasta index file obtained with samtools faidx. Required for coordinate sorting
        * output_bed_fn
            Path to write a summary result file in BED format
        * output_tsv_fn
            Path to write an more extensive result report in TSV format
        * min_depth
            Minimal number of reads covering a site to be reported
        * sample_id
            Sample ID to be used for the bed track header
        * strand_specific
            If True, output strand specific sites
        * min_llr
            Minimal log likelyhood ratio to consider a site significantly methylated or unmethylated
        * verbose
            Increase verbosity
        * quiet
            Reduce verbosity
        """
        # Save init options in dict for later
        kwargs = locals()

        # Define overall verbose level
        log = get_logger(name="Freq_meth_calculate", verbose=verbose, quiet=quiet)

        # Print option summary log
        log.debug ("## Options summary ##")
        log.debug ("\tpackage_name: {}".format(package_name))
        log.debug ("\tpackage_version: {}".format(package_version))
        log.debug ("\ttimestamp: {}".format(str(datetime.datetime.now())))
        log.debug (dict_to_str(kwargs, nsep=1, exclude_list=["self"]))

        # Verify parameters validity
        log.warning ("## Checking arguments ##")

        # Try to read input file if not a stream
        log.debug("\tTesting input file readability")
        if input_fn != 0 and not file_readable (input_fn):
            raise IOError ("Cannot read input file")

        # Verify that at least one output file is given:
        log.debug("\tCheck output file")
        if not output_bed_fn and not output_tsv_fn:
            raise pycoMethError("At least one output file should be given")
        if output_bed_fn:
            if os.path.dirname(output_bed_fn):
                mkdir (os.path.dirname(output_bed_fn), exist_ok=True)
            log.debug("\t\tOutput results in bed format")
        if output_tsv_fn:
            if os.path.dirname(output_tsv_fn):
                mkdir (os.path.dirname(output_tsv_fn), exist_ok=True)
            log.debug("\t\tOutput results in tsv format")

        # Create self variables
        counter = Counter()

        log.warning ("## Parsing methylation_calls file ##")
        # Init SGC class with fasta_index
        if fasta_index:
            SGC.set_chrom_list(fasta_index)

        # Create collection to store results
        site_dict = defaultdict(list)

        try:
            input_fp = open (input_fn, "r")

            log.info ("\tStarting to parse file Nanopolish methylation call file")
            header_line = input_fp.readline()
            byte_offset = len(header_line)
            lp = LineParser(header_line, sep="\t", cast_numeric_field=True)

            for line in tqdm(input_fp, desc="\t", unit=" lines", disable=log.level>=30):
                counter["Total read lines"]+=1
                byte_len = len(line)
                l = lp(line)

                if not l:
                    # Failsafe if line is malformed
                    counter["Invalid read line"]+=1
                else:
                    # Store byte offset corresponding to appropriate line
                    counter["Valid read lines"]+=1
                    if strand_specific:
                        coord = SGC(l.chromosome, l.start, l.strand)
                    else:
                        coord = SGC(l.chromosome, l.start)
                    site_dict[coord].append(byte_offset)
                    byte_offset += byte_len

            log.info ("\tFiltering out low coverage sites")
            filtered_site_dict = defaultdict(list)
            for k, offset_list in site_dict.items():
                counter["Total sites"]+=1

                # If low coverage unset list to release memory
                if len(offset_list) < min_depth:
                    counter["Low coverage sites"]+=1
                else:
                    counter["Valid sites"]+=1
                    filtered_site_dict[k]=offset_list
            del site_dict

            if not filtered_site_dict:
                raise pycoMethError ("No valid sites left after coverage filtering")

            if fasta_index:
                log.info ("\tSorting by coordinates")
                filtered_site_dict = OrderedDict(sorted(filtered_site_dict.items(), key=lambda t: t[0]))

            log.info ("\tProcessing valid sites found")
            Site.set_class_param(strand_specific=strand_specific, min_llr=min_llr)

            log.debug ("\t\tWrite output file header")
            if output_bed_fn:
                output_bed_fp = open (output_bed_fn, "w")
                output_bed_fp.write(Site.BED_header(sample_id)+"\n")
            if output_tsv_fn:
                output_tsv_fp = open (output_tsv_fn, "w")
                output_tsv_fp.write(Site.TSV_header()+"\n")

            for k, offset_list in tqdm(filtered_site_dict.items(), desc="\t", unit=" sites", disable=log.level>=30):
                # Get all read lines corresponding to current site
                ll = []
                for offset in offset_list:
                    input_fp.seek(offset, 0)
                    ll.append(lp(input_fp.readline()))

                # Parse list with helper class Site
                site = Site(ll=ll)
                counter["Valid sites"]+=1
                if output_bed_fn:
                    output_bed_fp.write(site.to_bed()+"\n")
                if output_tsv_fn:
                    output_tsv_fp.write(site.to_tsv()+"\n")

        finally:
            input_fp.close()
            if output_bed_fn:
                output_bed_fp.close()
            if output_tsv_fn:
                output_tsv_fp.close()

        log.info ("## Results summary ##")
        log.info (dict_to_str(counter, nsep=1))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~HELPER CLASS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Site():
    """Structure like class to store site information"""

    # Class variables and setter
    strand_specific = False
    min_llr = 2
    @classmethod
    def set_class_param (cls, strand_specific=False, min_llr=2):
        cls.strand_specific = strand_specific
        cls.min_llr = min_llr

    @classmethod
    def BED_header (cls, sample_id=""):
        return "track name=Methylation_{} itemRgb=On".format(sample_id)

    @classmethod
    def TSV_header (cls):
        return "\t".join([
            "chromosome",
            "start",
            "end",
            "strand",
            "site_id",
            "methylated_reads",
            "unmethylated_reads",
            "ambiguous_reads",
            "sequence",
            "num_motifs",
            "llr_list"])

    def __init__ (self, ll):
        """"""
        self.total = len(ll)
        self.methylated = 0
        self.unmethylated = 0
        self.ambiguous = 0
        self.id = "{}:{}-{}".format(ll[0].chromosome, ll[0].start, ll[0].end+1)
        self.sequence = ll[0].sequence
        self.num_motifs = ll[0].num_motifs
        self.chromosome = ll[0].chromosome
        self.start = ll[0].start
        self.end = ll[0].end+1
        self.strand = ll[0].strand if self.strand_specific else "."

        self.llr_list = []
        for l in ll:
            self.llr_list.append(float(l.log_lik_ratio))
            # Count read methylation call per site
            if l.log_lik_ratio >= self.min_llr:
                self.methylated+=1
            elif l.log_lik_ratio <= -self.min_llr:
                self.unmethylated+=1
            else:
                self.ambiguous+=1

        self.med_llr = np.mean(self.llr_list)
        if self.med_llr <= -self.min_llr:
            self.color = '8,121,207'
        elif self.med_llr < self.min_llr:
            self.color = '100,100,100'
        else:
            self.color = '235,5,79'

    def __repr__(self):
        return "{}:{}-{}({}) / id:{} / reads:{}".format(
            self.chromosome,
            self.start,
            self.end,
            self.strand,
            self.id,
            self.total)

    def to_bed (self):
        """"""
        return "{}\t{}\t{}\t{}\t{:.6f}\t{}\t{}\t{}\t'{}'".format(
            self.chromosome,
            self.start,
            self.end,
            self.id,
            self.med_llr,
            self.strand,
            self.start,
            self.end,
            self.color)

    def to_tsv (self):
        """"""
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
            self.chromosome,
            self.start,
            self.end,
            self.strand,
            self.id,
            self.methylated,
            self.unmethylated,
            self.ambiguous,
            self.sequence,
            self.num_motifs,
            ",".join([str(i) for i in self.llr_list]))

class SGC():
    """Sortable genomic coordinate object"""

    chr_list = OrderedDict()
    @classmethod
    def set_chrom_list (cls, index):
        if isinstance(index, dict):
            cls.chr_list = index
        else:
            with open(index) as fp:
                for i, line in enumerate(fp):
                    chrom = line.split()[0]
                    cls.chr_list[chrom]=i

    def __init__ (self, chrom, start=0, strand=""):
        # Verify and store chromosome name
        self.chrom = chrom
        self.start = int(start)
        self.strand = strand

    def __repr__ (self):
        return "{}:{} ({})".format(self.chrom, self.start, self.strand)

    def __lt__(self, other):
        if self.chrom == other.chrom:
            if self.start == other.start:
                # If everything is the same, just pick a random one
                if self.strand == other.strand:
                    return randint(0,1)
                else:
                    return self.strand < other.strand
            else:
                return self.start < other.start
        else:
            return self.chr_list.get(self.chrom, 0) < self.chr_list.get(other.chrom, 0)

    def __hash__(self):
        return hash((self.chrom, self.start, self.strand))

    def __eq__(self, other):
        return (self.chrom, self.start, self.strand) == (other.chrom, other.start, other.strand)

    def __ne__(self, other):
        return not(self == other)
