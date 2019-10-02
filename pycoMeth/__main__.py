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
from pycoMeth.Freq_meth_calculate import Freq_meth_calculate

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

    # Freq_meth_calculate subparser
    subparser_fm = subparsers.add_parser("Freq_meth_calculate", description="Calculate methylation frequency at genomic CpG sites from the output of nanopolish call-methylation")
    subparser_fm.set_defaults(func=Freq_meth_calculate_main)
    subparser_fm_io = subparser_fm.add_argument_group("Input/Output options")
    subparser_fm_io.add_argument("-i", "--input_fn", default=0, help="Path to a nanopolish call_methylation tsv output file. If not specified read from std input")
    subparser_fm_io.add_argument("-b", "--output_bed_fn", type=str, default="", help="Path to write a summary result file in BED format")
    subparser_fm_io.add_argument("-t", "--output_tsv_fn", type=str, default="", help="Path to write an more extensive result report in TSV format")
    subparser_fm_fo = subparser_fm.add_argument_group("Filtering options")
    subparser_fm_fo.add_argument("-d", "--min_depth", type=int, default=10, help="Minimal number of reads covering a site to be reported (default: %(default)s)")
    subparser_fm_other = subparser_fm.add_argument_group("Other options")
    subparser_fm_other.add_argument("-f", "--fasta_index", type=str, default="", help="fasta index file obtained with samtools faidx. Required for coordinate sorting")
    subparser_fm_other.add_argument("-s", "--sample_id", type=str, default="", help="Sample ID to be used for the bed track header (default: %(default)s)")
    subparser_fm_other.add_argument("--strand_specific", action="store_true", default=False, help="Output strand specific sites")
    subparser_fm_other.add_argument("--min_llr", type=float, default=2, help="Minimal log likelyhood ratio to consider a site significantly methylated or unmethylated (default: %(default)s)")

    # Add common group parsers
    for sp in [subparser_ec, subparser_fm]:
        sp_verbosity = sp.add_mutually_exclusive_group()
        sp_verbosity.add_argument("-v", "--verbose", action="store_true", default=False, help="Increase verbosity")
        sp_verbosity.add_argument("-q", "--quiet", action="store_true", default=False, help="Reduce verbosity")

    # Parse args and call subfunction
    args = parser.parse_args()
    args.func(args)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SUBPARSERS FUNCTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

def Freq_meth_calculate_main (args):
    """"""
    # Run corresponding class
    Freq_meth_calculate (
        input_fn = args.input_fn,
        fasta_index=args.fasta_index,
        output_bed_fn = args.output_bed_fn,
        output_tsv_fn = args.output_tsv_fn,
        min_depth = args.min_depth,
        sample_id = args.sample_id,
        strand_specific = args.strand_specific,
        min_llr = args.min_llr,
        verbose = args.verbose,
        quiet = args.quiet)
