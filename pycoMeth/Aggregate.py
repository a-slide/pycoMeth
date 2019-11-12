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

#~~~~~~~~~~~~~~~~~~~~~~~~Aggregate MAIN CLASS~~~~~~~~~~~~~~~~~~~~~~~~#

def Aggregate(
    input_fn:str,
    fasta_index:str,
    output_bed_fn:str="",
    output_tsv_fn:str="",
    min_depth:int=10,
    sample_id:str="",
    min_llr:float=2,
    **kwargs):
    """
    Calculate methylation frequency at genomic CpG sites from the output of nanopolish call-methylation
    * input_fn
        Path to a nanopolish call_methylation tsv output file
    * fasta_index
        fasta index file obtained with samtools faidx needed for coordinate sorting
    * output_bed_fn
        Path to write a summary result file in BED format (At least 1 output file is requires in CLI mode)
    * output_tsv_fn
        Path to write an more extensive result report in TSV format (At least 1 output file is requires in CLI mode)
    * min_depth
        Minimal number of reads covering a site to be reported
    * sample_id
        Sample ID to be used for the bed track header
    * min_llr
        Minimal log likelyhood ratio to consider a site significantly methylated or unmethylated
    * kwargs
        Allow to pass extra options such as verbose, quiet and progress
    """

    # Save init options in dict for later
    args = locals()

    verbose = kwargs.get("verbose", False)
    quiet = kwargs.get("quiet", False)
    progress = kwargs.get("progress", False)
    log = get_logger (name="pycoMeth_Aggregate", verbose=verbose, quiet=quiet)

    # Print option summary log
    log.debug ("## Options summary ##")
    log.debug ("\tpackage_name: {}".format(package_name))
    log.debug ("\tpackage_version: {}".format(package_version))
    log.debug ("\ttimestamp: {}".format(str(datetime.datetime.now())))
    log.debug (dict_to_str(args, nsep=1))

    # Verify parameters validity
    log.warning ("## Checking arguments ##")

    # At least one output file is required, otherwise it doesn't make any sense
    if not output_bed_fn and not output_tsv_fn:
        raise pycoMethError ("At least 1 output file is requires (-t or -b)")

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

    log.warning("## Parsing methylation_calls file ##")
    # Init SitesIndex object with fasta_index to aggregate data at genomic position level
    sites_index = SitesIndex(chr_index=fasta_index)

    # Open file parser
    # Possible fields chromosome	strand	start	end	read_name	log_lik_ratio	log_lik_methylated	log_lik_unmethylated	num_calling_strands	num_motifs	sequence
    dtypes = {"start":int, "end":int, "log_lik_ratio":float, "num_motifs":int}
    with FileParser(fn=input_fn, dtypes=dtypes, verbose=verbose, quiet=quiet, include_byte_offset=progress) as np_call_fp:

        log.info ("\tStarting to parse file Nanopolish methylation call file")
        with tqdm (total=len(np_call_fp), desc="\t", unit=" bytes", unit_scale=True, disable=not progress) as pbar:
            prev_byte_len = 0
            for lt in np_call_fp:
                sites_index.add(lt)

                # Update progress_bar
                if progress:
                    pbar.update(lt.byte_offset-prev_byte_len)
                    prev_byte_len = lt.byte_offset

        log.info ("\tFiltering out low coverage sites")
        sites_index.filter_low_count(min_depth)

        log.info ("\tSorting each chromosome by coordinates")
        sites_index.sort()

        log.info ("## Parsing summary ##")
        log.info (dict_to_str (np_call_fp.counter, nsep=1))
        log.info (dict_to_str (sites_index.counter, nsep=1))

    log.info ("\tProcessing valid sites found and write to file")
    with SitesWriter( bed_fn=output_bed_fn, tsv_fn=output_tsv_fn, sample_id=sample_id, min_llr=min_llr) as sites_writer:
        for coord, val_dict in tqdm(sites_index, desc="\t", unit=" sites", unit_scale=True, disable=not progress):
            sites_writer.write (coord, val_dict)

        log.info ("## Results summary ##")
        log.info (dict_to_str (sites_writer.counter, nsep=1))

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

    def __iter__(self):
        for chrom, chrom_d in self.sites.items():
            for (start, end), val_dict in chrom_d.items():
                yield ((chrom, start, end), val_dict)

    #~~~~~~~~~~~~~~PUBLIC METHODS~~~~~~~~~~~~~~#

    def add (self, lt):
        """"""
        chrom = lt.chromosome
        pos = (lt.start,lt.end+1)

        if chrom not in self.sites:
            self.counter ["Invalid chromosome lines"]+=1

        else:
            if not pos in self.sites[chrom]: ####################################################### Could also collect read_name, log_lik_methylated, log_lik_unmethylated
                self.counter ["Initial Sites"]+=1
                self.sites[chrom][pos]={
                    "sequence":lt.sequence,
                    "num_motifs":lt.num_motifs,
                    "llr":[],
                    "n_reads":0}

            self.counter ["Total Valid Lines"]+=1
            self.sites[chrom][pos]["llr"].append(lt.log_lik_ratio)
            self.sites[chrom][pos]["n_reads"]+=1

    def filter_low_count (self, min_count=0):
        """"""
        filtered_sites = OrderedDict()
        for chrom, chrom_d in self.sites.items():
            filtered_pos = OrderedDict()
            for pos, val_dict in chrom_d.items():
                if val_dict["n_reads"]<min_count:
                    self.counter ["Low Count Sites"]+=1
                else:
                    self.counter ["Valid Sites Found"]+=1
                    filtered_pos[pos] = val_dict
            if filtered_pos:
                filtered_sites[chrom] = filtered_pos

        self.sites = filtered_sites
        if len(self) == 0:
            raise pycoMethError ("No valid sites left after coverage filtering")

    def sort (self):
        """"""
        for k in self.sites.keys():
            self.sites[k] = OrderedDict(sorted(self.sites[k].items(), key=lambda t: t[0]))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~SitesWriter HELPER CLASS~~~~~~~~~~~~~~~~~~~~~~~~~~~#

class SitesWriter():
    """Extract data for valid sites and write to BED and/or TSV file"""

    def __init__ (self, bed_fn="", tsv_fn="", sample_id="", min_llr=2,):
        """"""
        self.min_llr = min_llr
        self.counter = Counter()
        self.bed = True if bed_fn else False
        self.tsv = True if tsv_fn else False

        # Init output files handlers is needed
        if self.bed:
            self.bed_fp = self._init_bed (bed_fn, sample_id)
        if self.tsv:
            self.tsv_fp = self._init_tsv (tsv_fn)

    #~~~~~~~~~~~~~~MAGIC AND PROPERTY METHODS~~~~~~~~~~~~~~#

    def __repr__ (self):
        return dict_to_str (self.counter)

    def __enter__ (self):
        return self

    def __exit__(self, exception_type, exception_val, trace):
        if self.bed:
            try:
                self.bed_fp.close()
            except Exception as E:
                pass
        if self.tsv:
            try:
                self.tsv_fp.close()
            except Exception as E:
                pass

    #~~~~~~~~~~~~~~PUBLIC METHODS~~~~~~~~~~~~~~#

    def write (self, coord, val_dict):
        """"""
        self.counter["Total Sites Written"]+=1

        # Collect llr values and count methylated, unmethylated and ambiguous reads
        llr_list = []
        reads_c = Counter()

        # Compute median llr and update counters
        med_llr = np.median(val_dict["llr"])
        if med_llr >= self.min_llr:
            self.counter["Methylated sites"]+=1
        elif med_llr <= -self.min_llr:
            self.counter["Unmethylated sites"]+=1
        else:
            self.counter["Ambiguous sites"]+=1

        # Write tsv and/or bed files
        if self.bed:
            self._write_bed (coord, med_llr)
        if self.tsv:
            self._write_tsv (coord, val_dict, med_llr)

    #~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~#

    def _init_bed (self, fn, sample_id):
        """Open BED file and write file header"""
        mkbasedir (fn, exist_ok=True)
        fp = open(fn, "w")
        if sample_id:
            fp.write("track name=Methylation_{} itemRgb=On\n".format(sample_id))
        else:
            fp.write("track name=Methylation itemRgb=On\n")
        return fp

    def _write_bed (self, coord, med_llr):
        """Write line to BED file"""
        # Define color for bed file
        if med_llr >= self.min_llr:
            color = '235,5,79'
        elif med_llr <= -self.min_llr:
            color = '8,121,207'
        else:
            color = '100,100,100'

        # Write line
        self.bed_fp.write ("{}\t{}\t{}\t.\t{:.3f}\t.\t{}\t{}\t'{}'\n".format(
            coord[0], coord[1], coord[2], med_llr, coord[1], coord[2], color))

    def _init_tsv (self, fn):
        """Open TSV file and write file header"""
        mkbasedir (fn, exist_ok=True)
        fp = open(fn, "w")
        fp.write("chromosome\tstart\tend\tsequence\tnum_motifs\tmedian_llr\tllr_list\n")
        return fp

    def _write_tsv (self, coord, val_dict, med_llr):
        """Write line to TSV file"""
        llr_str = ";".join([str(i) for i in val_dict["llr"]])
        self.tsv_fp.write ("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            coord[0], coord[1], coord[2], val_dict["sequence"], val_dict["num_motifs"], med_llr, llr_str))
