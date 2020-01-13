# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Standard library imports
from collections import OrderedDict, namedtuple, Counter

# Third party imports
from tqdm import trange

# Third party imports
from pyfaidx import Fasta

# Local imports
from pycoMeth.common import *

#~~~~~~~~~~~~~~~~~~~~~~~~CpG_Comp MAIN CLASS~~~~~~~~~~~~~~~~~~~~~~~~#

def CGI_Finder (
    ref_fasta_fn:str,
    output_tsv_fn:str=None,
    output_bed_fn:str=None,
    merge_gap:int=0,
    min_win_len:int=200,
    min_CG_freq:float=0.5,
    min_obs_CG_ratio:float=0.6,
    verbose:bool=False,
    quiet:bool=False,
    progress:bool=False,
    **kwargs):
    """
    Simple method to find putative CpG islands in DNA sequences by using a sliding window and merging
    overlapping windows satisfying the CpG island definition.
    Results can be saved in bed and tsv format
    * ref_fasta_fn
        Reference file used for alignment in Fasta format (ideally already indexed with samtools faidx)
    * output_bed_fn
        Path to write a summary result file in BED format (At least 1 output file is required)
    * output_tsv_fn
        Path to write an more extensive result report in TSV format (At least 1 output file is required)
    * merge_gap
        Merge close CpG island within a given distance in bases
    * min_win_len
        Length of the minimal window containing CpG. Used as the sliding window length
    * min_CG_freq
        Minimal C+G frequency in a window to be counted as a valid CpG island
    * min_obs_CG_ratio
        Minimal Observed CG dinucleotidefrequency over expected distribution in a window to be counted as a valid CpG island
    """

    # Init method
    opt_summary_dict = opt_summary(local_opt=locals())
    log = get_logger (name="pycoMeth_CGI_Finder", verbose=verbose, quiet=quiet)

    log.warning("Checking options and input files")

    log.debug ("Options summary")
    for i,j in opt_summary_dict.items():
        log.debug ("\t{}: {}".format(i,j))

    # Init collections
    counter = Counter()

    # At least one output file is required, otherwise it doesn't make any sense
    log.debug ("Checking required output")
    if not output_bed_fn and not output_tsv_fn:
        raise pycoMethError ("At least 1 output file is requires (-t or -b)")

    log.warning("Parsing reference fasta file")
    try:
        with Fasta(ref_fasta_fn) as fasta_fp:
            with CGI_Writer(bed_fn=output_bed_fn, tsv_fn=output_tsv_fn, verbose=verbose) as writer:

                # Iterate over reference sequences in fasta file
                for seq in fasta_fp:

                    # Hard copy of seq and cast to lower case
                    seq_name = seq.name
                    seq = str(seq).lower()
                    log.info ("Parsing Reference sequence: {}".format(seq_name))
                    counter["Number of reference sequences"]+=1

                    # Loop control counters
                    valid_win_start = valid_win_end = 0
                    previous_valid = False

                    # Compute and evaluate first window
                    c_count = g_count = cg_count = 0
                    for i in range(min_win_len):
                        if seq[i] == "c":
                            c_count+=1
                        elif seq[i] == "g":
                            g_count+=1
                        if seq[i:i+2] =="cg":
                            cg_count+=1

                    if valid_window (c_count, g_count, cg_count, min_win_len, min_CG_freq, min_obs_CG_ratio):
                        counter["Valid minimal size windows"]+=1
                        valid_win_start = i
                        valid_win_end = i+min_win_len
                        previous_valid = True

                    for i in trange(1, len(seq)-min_win_len, unit=" bases", unit_scale=True, disable=not progress):

                        # Decrement counters based on changes at previous start
                        prev_start = i-1
                        if seq[prev_start] == "c":
                            c_count-=1
                        elif seq[prev_start] == "g":
                            g_count-=1
                        if seq[prev_start:prev_start+2] =="cg":
                            cg_count-=1

                        # Increment counters based on changes at new end
                        end = i+min_win_len
                        if seq[end] == "c":
                            c_count+=1
                        elif seq[end] == "g":
                            g_count+=1
                        if seq[end-1:end+1] =="cg":
                            cg_count+=1

                        # Evaluate windows
                        if valid_window (c_count, g_count, cg_count, min_win_len, min_CG_freq, min_obs_CG_ratio):
                            counter["Valid minimal size windows"]+=1

                            # Previous valid overlapping windows was valid => extend end
                            if previous_valid:
                                valid_win_end = i+min_win_len

                                # Special case where reaching end with valid window => save window
                                if i == len(seq)-min_win_len-1:
                                    counter["Valid merged windows"]+=1
                                    win_len, win_cg_count, win_cg_freq, win_obs_exp = compute_win (seq[valid_win_start:valid_win_end]) # Sometimes valid merged windows are not excatly matching the definition of a proper CpG_island
                                    writer.write(seq_name, valid_win_start, valid_win_end, win_len, win_cg_count, win_cg_freq, win_obs_exp)

                            # Not overlapping a previous valid window => Start new one
                            else:
                                valid_win_start = i
                                valid_win_end = i+min_win_len

                            previous_valid=True

                        # A previous overlapping windows was valid and getting out of the valid window or reaching end of sequence => save window
                        elif i > (valid_win_end+merge_gap) or i == len(seq)-min_win_len-1:
                            if previous_valid:
                                counter["Valid merged windows"]+=1
                                win_len, win_cg_count, win_cg_freq, win_obs_exp = compute_win (seq[valid_win_start:valid_win_end])
                                writer.write(seq_name, valid_win_start, valid_win_end, win_len, win_cg_count, win_cg_freq, win_obs_exp)

                            previous_valid=False


    finally:
        # Print counters
        log.info ("Results summary")
        for i,j in counter.items():
            log.info ("\t{}: {:,}".format(i,j))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~Comp_Writer HELPER CLASS~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class CGI_Writer():
    """Extract data for valid sites and write to BED and/or TSV file"""

    def __init__ (self, bed_fn=None, tsv_fn=None, verbose=True):
        """"""
        self.log = get_logger (name="Comp_Writer", verbose=verbose)
        self.bed_fn = bed_fn
        self.tsv_fn = tsv_fn
        self.bed_fp = self._init_bed () if bed_fn else None
        self.tsv_fp = self._init_tsv () if tsv_fn else None

    #~~~~~~~~~~~~~~PUBLIC METHODS~~~~~~~~~~~~~~#
    def write (self, chrom, start, end, length, n_cpg, cg_freq, obs_exp):
        """"""
        if self.bed_fn:
            self._write_bed(chrom, start, end)
        if self.tsv_fn:
            self._write_tsv(chrom, start, end, length, n_cpg, cg_freq, obs_exp)

    def __enter__ (self):
        self.log.debug("Opening Writer")
        return self

    def __exit__(self, exception_type, exception_val, trace):
        self.log.debug("Closing Writer")
        for fp in (self.bed_fp, self.tsv_fp):
            try:
                fp.close()
            except:
                pass

    #~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~#
    def _init_bed (self):
        """Open BED file and write file header"""
        self.log.debug("Initialise output bed file")
        mkbasedir (self.bed_fn, exist_ok=True)
        fp = open(self.bed_fn, "w")
        fp.write("track name=CpG_islands\n")
        return fp

    def _write_bed (self, chrom, start, end):
        """Write line to BED file"""
        self.bed_fp.write("{}\t{}\t{}\n".format(chrom, start, end))

    def _init_tsv (self):
        """Open TSV file and write file header"""
        self.log.debug("Initialise output tsv file")
        mkbasedir (self.tsv_fn, exist_ok=True)
        fp = open(self.tsv_fn, "w")
        fp.write("chromosome\tstart\tend\tlength\tnum_CpG\tCG_freq\tobs_exp_freq\n")
        return fp

    def _write_tsv (self, chrom, start, end, length, n_cpg, cg_freq, obs_exp):
        """Write line to TSV file"""
        self.tsv_fp.write ("{}\t{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}\n".format(chrom, start, end, length, n_cpg, cg_freq, obs_exp))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~Helper Functions~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def valid_window (c_count, g_count, cg_count, min_win_len=200, min_CG_freq=0.5, min_obs_CG_ratio=0.6):
    cg_freq = (c_count+g_count)/min_win_len
    if cg_freq < min_CG_freq:
        return False

    # Safely compute obs_exp in case of no Cs or Gs on windows
    try:
        obs_exp = cg_count/(c_count*g_count/min_win_len)
    except ZeroDivisionError:
        obs_exp = 0

    if obs_exp < min_obs_CG_ratio:
        return False

    return True

def compute_win (win_seq):
    win_len = len(win_seq)
    c_count = g_count = cg_count = 0

    for j in range(win_len):
        if win_seq[j] == "c":
            c_count+=1
        elif win_seq[j] == "g":
            g_count+=1
        if win_seq[j:j+2] =="cg":
            cg_count+=1

    cg_freq = (c_count+g_count)/win_len
    obs_exp = cg_count/(c_count*g_count/win_len)

    return (win_len, cg_count, cg_freq, obs_exp)
