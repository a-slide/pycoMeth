# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Standard library imports
from collections import *
import datetime

# Third party imports
from tqdm import tqdm
import numpy as np
import pandas as pd

# Local imports
from pycoMeth.common import *
from pycoMeth.FileParser import FileParser
from pycoMeth import __version__ as package_version
from pycoMeth import __name__ as package_name

#~~~~~~~~~~~~~~~~~~~~~~~~Freq_meth_calculate MAIN CLASS~~~~~~~~~~~~~~~~~~~~~~~~#

class Freq_meth_calculate():

    def __init__ (self,
        input_fn:"str",
        fasta_index:"str",
        output_bed_fn:"str"="",
        output_tsv_fn:"str"="",
        min_depth:"int"=10,
        sample_id:"str"="",
        min_llr:"float"=2,
        verbose:"bool"=False,
        quiet:"bool"=False):
        """
        Calculate methylation frequency at genomic CpG sites from the output of nanopolish call-methylation
        * input_fn
            Path to a nanopolish call_methylation tsv output file
        * fasta_index
            fasta index file obtained with samtools faidx needed for coordinate sorting
        * output_bed_fn
            Path to write a summary result file in BED format
        * output_tsv_fn
            Path to write an more extensive result report in TSV format
        * min_depth
            Minimal number of reads covering a site to be reported
        * sample_id
            Sample ID to be used for the bed track header
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
        if output_bed_fn:
            log.debug("\t\tOutput results in BED format")
        if output_tsv_fn:
            log.debug("\t\tOutput results in TSV format")
        if not output_bed_fn and not output_tsv_fn:
            log.debug("\t\tInteractive mode. Save result in memory")

        log.warning ("## Parsing methylation_calls file ##")
        # Init SitesIndex object with fasta_index to aggregate data at genomic position level
        sites_index = SitesIndex(chr_index=fasta_index)

        # Open file parser
        with FileParser(
            fn=input_fn, sep="\t", first_line_header=True, include_byte_offset=True,
            dtypes={"start":int, "end":int, "log_lik_ratio":float, "num_motifs":int}) as np_call_fp:

            log.info ("\tStarting to parse file Nanopolish methylation call file")

            with tqdm (total=len(np_call_fp), desc="\t", unit=" bytes", unit_scale=True, disable=log.level>=30) as pbar:
                prev_byte_len = 0
                for lt in np_call_fp:
                    sites_index (lt.chromosome, lt.start, lt.byte_offset)
                    pbar.update(lt.byte_offset-prev_byte_len)
                    prev_byte_len = lt.byte_offset

            log.info ("\tFiltering out low coverage sites")
            sites_index.filter_low_count(min_depth)

            log.info ("\tSorting by coordinates")
            sites_index.sort()

            log.info ("\tProcessing valid sites found and write to file")
            if output_bed_fn: log.info ("\t\tStart writing BED output")
            if output_tsv_fn: log.info ("\t\tStart writing TSV output")

            with SitesWriter(
                bed_fn=output_bed_fn,
                tsv_fn=output_tsv_fn,
                sample_id=sample_id,
                min_llr=min_llr) as sites_writer:

                for chrom, pos, byte_offset_list in tqdm(sites_index, desc="\t", unit=" sites", unit_scale=True, disable=log.level>=30):
                    # Summarize and write all lines corresponding to a given genomic position
                    sites_writer (ll=np_call_fp.get_lines(byte_offset_list))

                log.info ("## Results summary ##")
                log.info (dict_to_str (np_call_fp.counter, nsep=1))
                log.info (dict_to_str (sites_index.counter, nsep=1))
                log.info (dict_to_str (sites_writer.counter, nsep=1))


                if not output_bed_fn and not output_tsv_fn:
                    self.df = sites_writer.df

#~~~~~~~~~~~~~~~~~~~~~~~~~~~SitesIndex HELPER CLASS~~~~~~~~~~~~~~~~~~~~~~~~~~~#

class SitesIndex():
    def __init__ (self, chr_index):
        """"""
        # Import chr list
        if isinstance(chr_index, list):
            chr_list = chr_index
        else:
            chr_list=[]
            with open(chr_index) as fp:
                for line in fp:
                    chr_list.append(line.split()[0])

        # Init dict with chromosomes names
        self.sites=OrderedDict()
        for c in chr_list:
            self.sites[c] = OrderedDict()

        # Init other self variables
        self.counter = Counter()

    #~~~~~~~~~~~~~~MAGIC AND PROPERTY METHODS~~~~~~~~~~~~~~#

    def __repr__(self):
        return dict_to_str(self.counter)

    def __len__(self):
        i=0
        for v in self.sites.values():
            i+=len(v)
        return i

    def __call__ (self, chrom, pos, byte_offset):
        """"""
        #assert type(pos) == int
        if chrom not in self.sites:
            self.counter ["Invalid chromosome lines"]+=1
        else:
            if not pos in self.sites[chrom]:
                self.counter ["Initial Sites"]+=1
                self.sites[chrom][pos]=[]
            self.counter ["Total Valid Lines"]+=1
            self.sites[chrom][pos].append(byte_offset)

    def __iter__(self):
        for chrom, chrom_d in self.sites.items():
            for pos, byte_offset_list in chrom_d.items():
                yield (chrom, pos, byte_offset_list)

    #~~~~~~~~~~~~~~PUBLIC METHODS~~~~~~~~~~~~~~#

    def filter_low_count (self, min_count=0):
        i = 0
        for k in self.sites.keys():
            filtered_d = OrderedDict()
            for i,j in self.sites[k].items():
                if len(j)<min_count:
                    self.counter ["Low Count Sites"]+=1
                else:
                    self.counter ["Valid Sites Found"]+=1
                    filtered_d[i]=j
                    i+=1
            self.sites[k] = filtered_d
        if i == 0:
            raise pycoMethError ("No valid sites left after coverage filtering")

    def sort (self):
        for k in self.sites.keys():
            self.sites[k] = OrderedDict(sorted(self.sites[k].items(), key=lambda t: t[0]))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~SitesWriter HELPER CLASS~~~~~~~~~~~~~~~~~~~~~~~~~~~#

class SitesWriter():
    """Extract data for valid sites and write to BED and/or TSV file"""

    def __init__ (self, bed_fn="", tsv_fn="", sample_id="", min_llr=2,):
        """"""
        self.sample_id = sample_id
        self.min_llr = min_llr
        self.counter = Counter()

        # Init output files handlers is needed
        self.bed_fp = self._init_bed (bed_fn)
        self.tsv_fp = self._init_tsv (tsv_fn)

        # Init list if running in non file mode
        self._ll = []
        self._l = namedtuple("l", [
            "chromosome",
            "start",
            "end",
            "strand",
            "methylated_reads",
            "unmethylated_reads",
            "ambiguous_reads",
            "sequence",
            "num_motifs",
            "median_llr",
            "llr_list"])

    #~~~~~~~~~~~~~~MAGIC AND PROPERTY METHODS~~~~~~~~~~~~~~#

    @property
    def df (self):
        return pd.DataFrame(self._ll)

    def __repr__ (self):
        return dict_to_str (self.counter)

    def __call__ (self, ll):
        """"""
        self.counter["Total Sites Written"]+=1

        # Shortcut for first line
        l0 = ll[0]

        # Collect llr values and count methylated, unmethylated and ambiguous reads
        llr_list = []
        reads_c = Counter()
        for lt in ll:
            reads_c["total sites"]+=1
            llr_list.append(lt.log_lik_ratio)
            # Count read methylation call per site
            if lt.log_lik_ratio >= self.min_llr:
                reads_c["methylated"]+=1
            elif lt.log_lik_ratio <= -self.min_llr:
                reads_c["unmethylated"]+=1
            else:
                reads_c["ambiguous"]+=1
        med_llr = np.median(llr_list)

        # Update counters and define color for bed file
        if med_llr >= self.min_llr:
            self.counter["Methylated sites"]+=1
            color = '235,5,79'
        elif med_llr <= -self.min_llr:
            self.counter["Unmethylated sites"]+=1
            color = '8,121,207'
        else:
            self.counter["Ambiguous sites"]+=1
            color = '100,100,100'

        # Write tsv line
        if self.tsv_fp:
            line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                l0.chromosome,
                l0.start,
                l0.end+1,
                ".",
                reads_c.get("methylated", 0),
                reads_c.get("unmethylated", 0),
                reads_c.get("ambiguous", 0),
                l0.sequence,
                l0.num_motifs,
                med_llr,
                ";".join([str(i) for i in llr_list]))
            self.tsv_fp.write (line)

        if self.bed_fp:
            line = "{}\t{}\t{}\t{}\t{:.2f}\t{}\t{}\t{}\t'{}'\n".format(
                l0.chromosome,
                l0.start,
                l0.end+1,
                ".",
                med_llr,
                ".",
                l0.start,
                l0.end+1,
                color)
            self.bed_fp.write (line)

        if not self.tsv_fp and not self.bed_fp:
            self._ll.append(
                self._l(
                    l0.chromosome,
                    l0.start,
                    l0.end+1,
                    ".",
                    reads_c.get("methylated", 0),
                    reads_c.get("unmethylated", 0),
                    reads_c.get("ambiguous", 0),
                    l0.sequence,
                    l0.num_motifs,
                    med_llr,
                    llr_list))

    def __enter__ (self):
        return self

    def __exit__(self, exception_type, exception_val, trace):
        for fp in (self.bed_fp, self.tsv_fp):
            if fp:
                try:
                    fp.close()
                except Exception as E:
                    print (E)

    #~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~#

    def _init_bed (self, fn):
        """Open BED file and write file header"""
        if fn:
            mkbasedir (fn, exist_ok=True)
            fp = open(fn, "w")
            if self.sample_id:
                fp.write("track name=Methylation_{} itemRgb=On\n".format(self.sample_id))
            else:
                fp.write("track name=Methylation itemRgb=On\n")
            return fp
        else:
            return False

    def _init_tsv (self, fn):
        """Open TSV file and write file header"""
        if fn:
            mkbasedir (fn, exist_ok=True)
            fp = open(fn, "w")
            fp.write("chromosome\tstart\tend\tstrand\tmethylated_reads\tunmethylated_reads\tambiguous_reads\tsequence\tnum_motifs\tmedian_llr\tllr_list\n")
            return fp
        else:
            return False
