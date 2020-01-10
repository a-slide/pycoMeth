# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import os
from collections import OrderedDict, namedtuple, Counter

# Third party imports
from pyfaidx import Fasta

# Local imports
from pycoMeth.common import *

#~~~~~~~~~~~~~~CLASS~~~~~~~~~~~~~~#
class CoordGen():
    def __init__ (self,
        ref_fasta_fn,
        verbose:bool=False,
        quiet:bool=False):
        """
        Initialise CoordTuple
        * ref_fasta_fn
            Reference file used for alignment in Fasta format (ideally already indexed with samtools faidx)
        """
        # Init logger
        self.log = get_logger (name="pycoMeth_CpG_Comp", verbose=verbose, quiet=quiet)

        # Import chr list
        self.chr_name_id = OrderedDict()
        self.chr_name_len = OrderedDict()

        self.log.debug ("Loading Fasta index:{}".format(ref_fasta_fn))
        with Fasta(ref_fasta_fn) as fa:
            for i, ref in enumerate(fa):
                self.chr_name_id[ref.name]=i
                self.chr_name_len[ref.name]=len(ref)

        self.log.debug ("Found {} reference sequences".format(len(self.chr_name_id)))

    def __iter__ (self):
        for name in self.chr_name_id.keys():
            yield("Name:{}\tID:{}\tLength:{}".format(name, self.chr_name_id[name], self.chr_name_len[name]))

    def __call__ (self, chr_name, start, end):
        """
        Check passed coordinates and generate a namedtuple containing a chromosome id, start end End
        The chromosome id is an integer corrersponding to the order of the reference in the index.
        The returned objects are thus easily sortable.
        * chr_name
            Name of the chromosome
        * start
            Start of the interval (has to be between 0 and chromosome length)
        * end
            End of the interval (has to be between start and chromosome length)
        """
        # Check Chr
        if not chr_name in self.chr_name_id:
            raise CoordTupleError ("Invalid chromosome name: {}".format(chr_name))

        # Extract chromosome len and id
        chr_len = self.chr_name_len[chr_name]
        chr_id = self.chr_name_id[chr_name]

        # Check Start
        try:
            start = int(start)
        except CoordTupleError:
            raise CoordTupleError ("Start coordinate is not a valid integer: {} ".format(start))
        if start < 0 or start > chr_len:
            raise CoordTupleError ("Invalid value for start coordinate: {} [0:{}]".format(start, chr_len))

        # Check End
        try:
            end = int(end)
        except CoordTupleError:
            raise CoordTupleError ("End coordinate is not a valid integer: {} ".format(end))
        if end < start or end > chr_len:
            raise CoordTupleError ("Invalid value for end coordinate: {} [{}:{}]".format(end, start, chr_len))

        return Coord(chr_id, chr_name, start, end)

class Coord():
    def __init__ (self, chr_id, chr_name, start, end):
        self.chr_id = chr_id
        self.chr_name = chr_name
        self.start = start
        self.end = end

    def __repr__ (self):
        return("chr_name:{}, start:{}, end:{}".format(self.chr_name, self.start, self.end))

    def __str__(self):
        return("{}:{}-{}".format(self.chr_name, self.start, self.end))

    def __hash__(self):
        return hash((self.chr_id, self.start, self.end))

    @property
    def center(self):
        return self.start+(self.end-self.start)/2

    # Comparison methods
    def __eq__(self, other):
        return (self.chr_id, self.start, self.end) == (other.chr_id, other.start, other.end)

    def __ne__(self, other):
        return (self.chr_id, self.start, self.end) != (other.chr_id, other.start, other.end)

    def __lt__(self, other):
        return (self.chr_id, self.start, self.end) < (other.chr_id, other.start, other.end)

    def __gt__(self, other):
        return (self.chr_id, self.start, self.end) > (other.chr_id, other.start, other.end)

    def __le__(self, other):
        return (self.chr_id, self.start, self.end) <= (other.chr_id, other.start, other.end)

    def __ge__(self, other):
        return (self.chr_id, self.start, self.end) >= (other.chr_id, other.start, other.end)

    def center_comp (self, other):
        # different chromosome
        if other.chr_id < self.chr_id:
            return "lower"
        elif other.chr_id > self.chr_id:
            return "greater"
        # same chromosome
        else:
            if other.center < self.start:
                return "lower"
            elif other.center > self.end:
                return "greater"
            else:
                return "inside"

class CoordTupleError (Exception):
    """ Basic exception class for FileParserError """
    pass
