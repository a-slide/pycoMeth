# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Standard library imports
from collections import OrderedDict, namedtuple, Counter
import gzip

# Third party imports
from tqdm import tqdm
import numpy as np
from pyfaidx import Fasta

# Local imports
from pycoMeth.common import *
from pycoMeth.FileParser import FileParser

#~~~~~~~~~~~~~~~~~~~~~~~~CpG_Aggregate MAIN CLASS~~~~~~~~~~~~~~~~~~~~~~~~#
def CpG_Aggregate(
    nanopolish_fn:[str],
    ref_fasta_fn:str,
    output_bed_fn:str="",
    output_tsv_fn:str="",
    min_depth:int=10,
    sample_id:str="",
    min_llr:float=2,
    verbose:bool=False,
    quiet:bool=False,
    progress:bool=False,
    **kwargs):
    """
    Calculate methylation frequency at genomic CpG sites from the output of `nanopolish call-methylation`
    * nanopolish_fn
        Path to a nanopolish call_methylation tsv output file or a list of files or a regex matching several files (can be gzipped)
    * ref_fasta_fn
        Reference file used for alignment in Fasta format (ideally already indexed with samtools faidx)
    * output_bed_fn
        Path to write a summary result file in BED format (At least 1 output file is required) (can be gzipped)
    * output_tsv_fn
        Path to write a more extensive result report in TSV format (At least 1 output file is required) (can be gzipped)
    * min_depth
        Minimal number of reads covering a site to be reported
    * sample_id
        Sample ID to be used for the BED track header
    * min_llr
        Minimal log likelyhood ratio to consider a site significantly methylated or unmethylated in output BED file
    """

    # Init package
    opt_summary_dict = opt_summary(local_opt=locals())
    log = get_logger (name="pycoMeth_CpG_Aggregate", verbose=verbose, quiet=quiet)

    log.debug ("Options summary")
    for i,j in opt_summary_dict.items():
        log.debug ("\t{}: {}".format(i,j))

    log.warning("Checking options and input files")

    # At least one output file is required, otherwise it doesn't make any sense
    if not output_bed_fn and not output_tsv_fn:
        raise pycoMethError ("At least 1 output file is requires (-t or -b)")

    # Init SitesIndex object with ref_fasta_fn to aggregate data at genomic position level
    log.warning ("Parsing methylation_calls file")
    sites_index = SitesIndex(ref_fasta_fn=ref_fasta_fn)

    # Open file parser
    # Possible fields chromosome	strand	start	end	read_name	log_lik_ratio	log_lik_methylated	log_lik_unmethylated	num_calling_strands	num_motifs	sequence
    dtypes = {"start":int, "end":int, "log_lik_ratio":float, "num_motifs":int}
    with FileParser(fn=nanopolish_fn, dtypes=dtypes, verbose=verbose, quiet=quiet, include_byte_len=progress) as input_fp:

        log.info ("Starting to parse file Nanopolish methylation call file")
        with tqdm (total=len(input_fp), unit=" bytes", unit_scale=True, disable=not progress) as pbar:
            for lt in input_fp:
                sites_index.add(lt)
                # Update progress_bar
                if progress: pbar.update(lt.byte_len)

        log.info ("Filtering out low coverage sites")
        sites_index.filter_low_count(min_depth)

        log.info ("Sorting each chromosome by coordinates")
        sites_index.sort()

        log.info ("Parsing summary")
        for i,j in input_fp.counter.items():
            log.info ("\t{}: {:,}".format(i,j))
        for i,j in sites_index.counter.items():
            log.info ("\t{}: {:,}".format(i,j))

    log.warning("Processing valid sites found and write to file")

    with CpG_Writer(bed_fn=output_bed_fn, tsv_fn=output_tsv_fn, sample_id=sample_id, min_llr=min_llr, verbose=verbose) as writer:
        for coord, val_dict in tqdm(sites_index, unit=" sites", unit_scale=True, disable=not progress):
            writer.write (coord, val_dict)

        log.info ("Results summary")
        for i,j in writer.counter.items():
            log.info ("\t{}: {:,}".format(i,j))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~SitesIndex HELPER CLASS~~~~~~~~~~~~~~~~~~~~~~~~~~~#

class SitesIndex():
    def __init__ (self, ref_fasta_fn):
        """"""

        # Init dict with chromosomes names
        self.sites=OrderedDict()
        with Fasta(ref_fasta_fn) as fa:
            for ref in fa:
                self.sites[ref.name] = OrderedDict()

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
        ct = namedtuple("ct", ["chr_name", "start", "end"])
        for chrom, chrom_d in self.sites.items():
            for (start, end), val_dict in chrom_d.items():
                yield (ct(chrom, start, end), val_dict)

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~CpG_Writer HELPER CLASS~~~~~~~~~~~~~~~~~~~~~~~~~~~#

class CpG_Writer():
    """Extract data for valid sites and write to BED and/or TSV file"""

    def __init__ (self, bed_fn="", tsv_fn="", sample_id="", min_llr=2, verbose=True):
        """"""
        self.log = get_logger (name="pycoMeth_CpG_Writer", verbose=verbose,)
        self.min_llr = min_llr
        self.sample_id = sample_id
        self.counter = Counter()
        self.bed_fn = bed_fn
        self.tsv_fn = tsv_fn
        self.bed_fp = self._init_bed () if bed_fn else None
        self.tsv_fp = self._init_tsv () if tsv_fn else None

        # Color score tables
        self.pos_colors = OrderedDict()
        self.pos_colors[min_llr+4]='172,0,38'
        self.pos_colors[min_llr+3]='205,11,33'
        self.pos_colors[min_llr+2]='231,34,30'
        self.pos_colors[min_llr+1]='249,73,40'
        self.pos_colors[min_llr]='252,118,53'
        self.pos_colors[0]='230,230,230'

        self.neg_colors = OrderedDict()
        self.neg_colors[-min_llr-4]='28,45,131'
        self.neg_colors[-min_llr-3]='35,70,156'
        self.neg_colors[-min_llr-2]='33,102,171'
        self.neg_colors[-min_llr-1]='29,140,190'
        self.neg_colors[-min_llr]='52,168,194'
        self.neg_colors[0]='230,230,230'

    #~~~~~~~~~~~~~~MAGIC AND PROPERTY METHODS~~~~~~~~~~~~~~#
    def __repr__ (self):
        return dict_to_str (self.counter)

    def __enter__ (self):
        return self

    def __exit__(self, exception_type, exception_val, trace):
        for fp in (self.bed_fp, self.tsv_fp):
            try:
                fp.close()
            except:
                pass

    #~~~~~~~~~~~~~~PUBLIC METHODS~~~~~~~~~~~~~~#
    def write (self, coord, val_dict):
        """"""
        self.counter["Total Sites Writen"]+=1

        # Compute median llr and update counters
        med_llr = round(np.median(val_dict["llr"]), 3)
        if med_llr >= self.min_llr:
            self.counter["Methylated sites"]+=1
        elif med_llr <= -self.min_llr:
            self.counter["Unmethylated sites"]+=1
        else:
            self.counter["Ambiguous sites"]+=1

        # Write tsv and/or bed files
        if self.bed_fn:
            self._write_bed (coord, med_llr)
        if self.tsv_fn:
            self._write_tsv (coord, val_dict, med_llr)

    #~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~#

    def _init_bed (self):
        """Open BED file and write file header"""
        self.log.debug("Initialise output bed file")
        mkbasedir (self.bed_fn, exist_ok=True)
        fp = gzip.open(self.bed_fn, "wt") if self.bed_fn.endswith(".gz") else open(self.bed_fn, "w")
        fp.write("track name={}_CpG itemRgb=On\n".format(self.sample_id))
        return fp

    def _write_bed (self, coord, med_llr):
        """Write line to BED file"""

        # define track color dependign on med llr
        if med_llr>=0:
            for min_llr, color in self.pos_colors.items():
                if med_llr >= min_llr:
                    break
        else:
            for min_llr, color in self.neg_colors.items():
                if med_llr <= min_llr:
                    break
        # Write line
        self.bed_fp.write ("{}\t{}\t{}\t.\t{:.3f}\t.\t{}\t{}\t{}\n".format(
            coord.chr_name, coord.start, coord.end, med_llr, coord.start, coord.end, color))

    def _init_tsv (self):
        """Open TSV file and write file header"""
        self.log.debug("Initialise output tsv file")
        mkbasedir (self.tsv_fn, exist_ok=True)
        fp = gzip.open(self.tsv_fn, "wt") if self.tsv_fn.endswith(".gz") else open(self.tsv_fn, "w")
        fp.write("chromosome\tstart\tend\tsequence\tnum_motifs\tmedian_llr\tllr_list\n")
        return fp

    def _write_tsv (self, coord, val_dict, med_llr):
        """Write line to TSV file"""
        llr_str = list_to_str(val_dict["llr"])
        self.tsv_fp.write ("{}\t{}\t{}\t{}\t{}\t{:.3f}\t{}\n".format(
            coord.chr_name, coord.start, coord.end, val_dict["sequence"], val_dict["num_motifs"], med_llr, llr_str))
