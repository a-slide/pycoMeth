#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import argparse
import sys

# Local imports
import pycoMeth as pkg
from pycoMeth.common import *
from pycoMeth.CpG_Aggregate import CpG_Aggregate
from pycoMeth.Interval_Aggregate import Interval_Aggregate
from pycoMeth.Meth_Comp import Meth_Comp
from pycoMeth.CGI_Finder import CGI_Finder

#~~~~~~~~~~~~~~TOP LEVEL ENTRY POINT~~~~~~~~~~~~~~#
def main(args=None):
    """ Main entry point for pycoMeth command line interface"""

    # Parser and subparsers for command
    parser = argparse.ArgumentParser (description=pkg.__description__)
    parser.add_argument("--version", action="version", version="{} v{}".format(pkg.__name__, pkg.__version__))
    subparsers = parser.add_subparsers (description="%(prog)s implements the following subcommands", dest="subcommands")
    subparsers.required = True

    # CpG_Aggregate subparser
    f = CpG_Aggregate
    sp_cpg = subparsers.add_parser("CpG_Aggregate", description=doc_func(f))
    sp_cpg.set_defaults(func=f)
    sp_cpg_io = sp_cpg.add_argument_group("Input/Output options")
    arg_from_docstr(sp_cpg_io, f, "nanopolish_fn", "i")
    arg_from_docstr(sp_cpg_io, f, "ref_fasta_fn", "f")
    arg_from_docstr(sp_cpg_io, f, "output_bed_fn", "b")
    arg_from_docstr(sp_cpg_io, f, "output_tsv_fn", "t")
    sp_cpg_ms = sp_cpg.add_argument_group("Misc options")
    arg_from_docstr(sp_cpg_ms, f, "min_depth", "d")
    arg_from_docstr(sp_cpg_ms, f, "sample_id", "s")
    arg_from_docstr(sp_cpg_ms, f, "min_llr", "l")

    # Interval_Aggregate subparser
    f = Interval_Aggregate
    sp_int = subparsers.add_parser("Interval_Aggregate", description=doc_func(f))
    sp_int.set_defaults(func=f)
    sp_int_io = sp_int.add_argument_group("Input/Output options")
    arg_from_docstr(sp_int_io, f, "cpg_aggregate_fn", "i")
    arg_from_docstr(sp_int_io, f, "ref_fasta_fn", "f")
    arg_from_docstr(sp_int_io, f, "interval_bed_fn", "a")
    arg_from_docstr(sp_int_io, f, "output_bed_fn", "b")
    arg_from_docstr(sp_int_io, f, "output_tsv_fn", "t")
    sp_int_ms = sp_int.add_argument_group("Misc options")
    arg_from_docstr(sp_int_ms, f, "interval_size", "n")
    arg_from_docstr(sp_int_ms, f, "min_cpg_per_interval", "m")
    arg_from_docstr(sp_int_ms, f, "sample_id", "s")
    arg_from_docstr(sp_int_ms, f, "min_llr", "l")

    # Meth_Comp subparser
    f = Meth_Comp
    sp_met = subparsers.add_parser("Meth_Comp", description=doc_func(f))
    sp_met.set_defaults(func=f)
    sp_met_io = sp_met.add_argument_group("Input/Output options")
    arg_from_docstr(sp_met_io, f, "aggregate_fn_list", "i")
    arg_from_docstr(sp_met_io, f, "ref_fasta_fn", "f")
    arg_from_docstr(sp_met_io, f, "output_bed_fn", "b")
    arg_from_docstr(sp_met_io, f, "output_tsv_fn", "t")
    sp_met_ms = sp_met.add_argument_group("Misc options")
    arg_from_docstr(sp_met_ms, f, "max_missing", "m")
    arg_from_docstr(sp_met_ms, f, "min_diff_llr", "l")
    arg_from_docstr(sp_met_ms, f, "sample_id_list", "s")
    arg_from_docstr(sp_met_ms, f, "pvalue_adj_method")
    arg_from_docstr(sp_met_ms, f, "pvalue_adj_alpha")

    # CGI_Finder subparser
    f = CGI_Finder
    sp_cgi = subparsers.add_parser("CGI_Finder", description=doc_func(f))
    sp_cgi.set_defaults(func=f)
    sp_cgi_io = sp_cgi.add_argument_group("Input/Output options")
    arg_from_docstr(sp_cgi_io, f, "ref_fasta_fn", "f")
    arg_from_docstr(sp_cgi_io, f, "output_bed_fn", "b")
    arg_from_docstr(sp_cgi_io, f, "output_tsv_fn", "t")
    sp_cgi_ms = sp_cgi.add_argument_group("Misc options")
    arg_from_docstr(sp_cgi_ms, f, "merge_gap", "m")
    arg_from_docstr(sp_cgi_ms, f, "min_win_len", "w")
    arg_from_docstr(sp_cgi_ms, f, "min_CG_freq", "c")
    arg_from_docstr(sp_cgi_ms, f, "min_obs_CG_ratio", "r")

    # Add common group parsers
    for sp in [sp_cpg, sp_int, sp_met, sp_cgi]:
        sp_vb = sp.add_argument_group("Verbosity options")
        sp_vb.add_argument("-v", "--verbose", action="store_true", default=False, help="Increase verbosity")
        sp_vb.add_argument("-q", "--quiet", action="store_true", default=False, help="Reduce verbosity")
        sp_vb.add_argument("-p", "--progress", action="store_true", default=False, help="Display a progress bar")

    # Parse args and call subfunction
    args = parser.parse_args()
    args.func(**vars(args))
