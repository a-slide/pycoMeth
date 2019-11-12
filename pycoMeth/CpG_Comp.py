# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Standard library imports
from collections import *
import datetime

# Third party imports
from tqdm import tqdm ###################################################################################
import numpy as np
import pandas as pd
import scipy as sp
from statsmodels.stats.multitest import multipletests

# Local imports
from pycoMeth.common import *
from pycoMeth.FileParser import FileParser
from pycoMeth import __version__ as package_version
from pycoMeth import __name__ as package_name

#~~~~~~~~~~~~~~~~~~~~~~~~CpG_Comp MAIN CLASS~~~~~~~~~~~~~~~~~~~~~~~~#

def CpG_Comp (
    fn_list:list,
    fasta_index:str,
    output_tsv_fn:str="",
    max_missing:int=0,
    **kwargs):
    """

    * kwargs
        Allow to pass extra options such as verbose, quiet and progress
    """

    # Save init options in dict for later
    args = locals()

    verbose = kwargs.get("verbose", False)
    quiet = kwargs.get("quiet", False)
    progress = kwargs.get("progress", False)  ###################################### NOT USED
    log = get_logger (name="pycoMeth_CpG_Comp", verbose=verbose, quiet=quiet)

    # Print option summary log
    log.debug ("## Options summary ##")
    log.debug ("\tpackage_name: {}".format(package_name))
    log.debug ("\tpackage_version: {}".format(package_version))
    log.debug ("\ttimestamp: {}".format(str(datetime.datetime.now())))
    log.debug (dict_to_str(args, nsep=1))

    # Verify parameters validity
    log.warning ("## Checking arguments ##")

    # If CLI mode ask for a least one output file, otherwise it doesn't make any sense
    if "func" in kwargs:
        log.debug("\tRunning in CLI mode")
        if not output_tsv_fn:
            raise pycoMethError ("output file required in in CLI mode")
    else:
        log.debug("\tRunning in API mode")

    # Define collections and variables
    min_samples = len(fn_list)-max_missing
    counter = Counter()

    log.warning("## Parsing files ##")
    try:

        log.info("\tRead input files header and check consistancy between files")
        colnames = set()
        fp_list = []
        for fn in fn_list:
            dtypes={"start":int,"end":int,"median_llr":float}
            fp = FileParser(fn=fn, dtypes=dtypes, verbose=verbose, quiet=quiet)
            # Check colnames
            if not colnames:
                colnames = set(fp.colnames)
            elif not colnames == set(fp.colnames):
                raise ValueError (f"Invalid field {fp.colnames} in file {fn}")
            fp_list.append(fp)

        log.info("\tStart asynchronous file parsing")

        # Define coordinate and result collection
        ct = CoordTuple(fasta_index=fasta_index)
        rt = namedtuple("results", ["chromosome","start","end", "n_samples", "statistic", "pvalue"])
        res_list = []

        # Read first line from each file
        coord_d = defaultdict(list)
        for fp in fp_list:
            line = fp.next_line ()
            coord_d[ct(line.chromosome,line.start,line.end)].append(fp)

        # Read following lines
        done = 0
        while True:
            # Evaluate if lower coord has enough samples
            lower_coord = sorted(coord_d.keys())[0]
            coord_fp_list = coord_d[lower_coord]

            # If enough samples we have a winner
            if len(coord_fp_list) >= min_samples:
                counter["Sites with enough samples"]+=1

                # check with both positive and negative values
                coord_med_list = []
                for fp in coord_fp_list:
                    line = fp.current_line()
                    coord_med_list.append(line.median_llr)

                # Test only if both neg and pos vals
                if min(coord_med_list)<0 and max(coord_med_list)>0:
                    coord_llr_list = []
                    for fp in coord_fp_list:
                        line = fp.current_line()
                        llr_list = [float(i) for i in line.llr_list.split(";")]
                        coord_llr_list.append(llr_list)

                    # Run a KW test
                    res = sp.stats.kruskal(*coord_llr_list)
                    res_list.append(rt(
                        ct.id_to_name(lower_coord[0]),
                        lower_coord[1],
                        lower_coord[2],
                        len(coord_llr_list),
                        res.statistic,
                        res.pvalue))
            else:
                counter["Sites with insufficient samples"]+=1

            # Remove lower entry and move fp to next sequence
            del(coord_d[lower_coord])
            for fp in coord_fp_list:
                try:
                    line = fp.next_line()
                    coord = ct(line.chromosome,line.start,line.end)

                    # Check if file is sorted by coordinates
                    prev_line = fp.previous_line()
                    if prev_line:
                        p_coord = ct(prev_line.chromosome,prev_line.start,prev_line.end)
                        if p_coord >= coord:
                            raise ValueError("Unsorted coordinate found in file {}: {}:{}-{} found after {}:{}-{}".format(
                                fp.fn, prev_line.chromosome,prev_line.start,prev_line.end, line.chromosome,line.start,line.end))
                    coord_d[coord].append(fp)

                except StopIteration:
                    done+=1

            # Exit condition = all file are finished
            if done == len(fp_list):
                break

        # Convert results to dataframe and correct pvalues for multiple tests
        log.info("\tCollect results and adjust pvalues")
        df = pd.DataFrame(res_list)
        df["adj_pvalue"] = multipletests(df["pvalue"], method="fdr_bh")[1]
        counter["Sites with significant adj pvalue (< 0.1)"] = len(df.query("adj_pvalue<=0.1"))

        if output_tsv_fn:
            log.info("\tWrite output tsv file")
            df.to_csv(output_tsv_fn, sep="\t", index=False)
        else:
            return df

    finally:
        log.info ("## Results summary ##")
        log.info (dict_to_str (counter, nsep=1))

        for fp in fp_list:
            fp.close()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~CoordTuple HELPER CLASS~~~~~~~~~~~~~~~~~~~~~~~~~~~#

class CoordTuple():

    def __init__ (self, fasta_index):
        """"""
        # Import chr list
        self.chr_name_id = OrderedDict()
        self.chr_name_len = OrderedDict()
        self.chr_id_name = OrderedDict()

        with open(fasta_index) as fp:
            for i, line in enumerate(fp):
                c_id = line.split()[0]
                c_len = int(line.split()[1])
                self.chr_name_id[c_id]=i
                self.chr_name_len[c_id]=c_len
                self.chr_id_name[i]=c_id

    def __call__ (self, chr, start, end):
        """"""
        # Check Chr
        if not chr in self.chr_name_id:
            raise ValueError ("Invalid chromosome name: {}".format(chr))
        # Check Start
        try:
            start = int(start)
        except ValueError:
            raise ValueError ("Start coordinate is not a valid integer: {} ".format(start))
        if start < 0 or start > self.chr_name_len[chr]:
            raise ValueError ("Invalid value for start coordinate: {} [0:{}]".format(start, self.chr_name_len[chr]))
        # Check End
        try:
            end = int(end)
        except ValueError:
            raise ValueError ("End coordinate is not a valid integer: {} ".format(end))
        if end < start or end > self.chr_name_len[chr]:
            raise ValueError ("Invalid value for end coordinate: {} [{}:{}]".format(start, start, self.chr_name_len[chr]))

        return (self.chr_name_id[chr], start, end)

    def id_to_name (self, chr_id):
        return self.chr_id_name.get(chr_id, None)
