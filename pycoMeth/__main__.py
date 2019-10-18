#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import argparse
import sys

# Local imports
from pycoMeth import __version__ as package_version
from pycoMeth import __name__ as package_name
from pycoMeth import __description__ as package_description
from pycoMeth.common import *
from pycoMeth.Aggregate import Aggregate

#~~~~~~~~~~~~~~TOP LEVEL ENTRY POINT~~~~~~~~~~~~~~#
def main(args=None):
    """
    Main entry point for pycoMeth command line interface
    """

    # Parser and subparsers for command
    parser = argparse.ArgumentParser (description=package_description)
    parser.add_argument("--version", action="version", version="{} v{}".format(package_name, package_version))
    subparsers = parser.add_subparsers (description="%(prog)s implements the following subcommands", dest="subcommands")
    subparsers.required = True

    # Aggregate subparser
    f = Aggregate
    sp_agg = subparsers.add_parser("Aggregate", description=doc_func(f))
    sp_agg.set_defaults(func=f)
    sp_agg_io = sp_agg.add_argument_group("Input/Output options")
    arg_from_docstr(sp_agg_io, f, "input_fn", "i")
    arg_from_docstr(sp_agg_io, f, "fasta_index", "f")
    arg_from_docstr(sp_agg_io, f, "output_bed_fn", "b")
    arg_from_docstr(sp_agg_io, f, "output_tsv_fn", "t")
    sp_agg_ms = sp_agg.add_argument_group("Misc options")
    arg_from_docstr(sp_agg_ms, f, "min_depth", "d")
    arg_from_docstr(sp_agg_ms, f, "min_llr", "l")
    arg_from_docstr(sp_agg_ms, f, "sample_id", "s")

    # Add common group parsers
    for sp in [sp_agg]:
        sp_vb = sp.add_mutually_exclusive_group()
        sp_vb.add_argument("-v", "--verbose", action="store_true", default=False, help="Increase verbosity")
        sp_vb.add_argument("-q", "--quiet", action="store_true", default=False, help="Reduce verbosity")
        sp_vb.add_argument("-p", "--progress", action="store_true", default=False, help="Display a progress bar")

    # Parse args and call subfunction
    args = parser.parse_args()
    args.func(**vars(args))
