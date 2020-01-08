# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import os
from collections import *
import gzip
from glob import iglob

# Local imports
from pycoMeth.common import *

#~~~~~~~~~~~~~~CLASS~~~~~~~~~~~~~~#

class FileParser ():
    def __init__ (self,
        fn,
        label="",
        colnames=False,
        first_line_header=True,
        sep="\t",
        comment="#",
        auto_numeric=False,
        include_byte_len=False,
        dtypes={},
        force_dtypes=False,
        force_col_len=True,
        verbose=False,
        quiet=False,
        **kwargs):
        """
        Open a parser ++ for field delimited file
        * fn
            Path to a field delimited file
        * label
            Label for the file of file group
        * colnames
            List of column names to use if not in first file line
        * sep
            field separator
        * comment
            skip any line starting with this string
        * auto_numeric
            Try to automatically cast fields values in int or float
        * include_byte_len
            Add byte len corresponding to each line
        * dtypes
            Dict corresponding to fields (based on colnames) to cast in a given python type
        * force_dtypes
            Raise an error if type casting fails
        * kwargs
            Allow to pass extra options such as verbose, quiet and progress
        """

        # Init logger and counter
        self.log = get_logger (name="pycoMeth_FileParser", verbose=verbose, quiet=quiet)
        self.counter = Counter()

        # Save self variables
        self.label = label
        self.sep = sep
        self.comment = comment
        self.first_line_header= first_line_header
        self.include_byte_len = include_byte_len
        self.auto_numeric = auto_numeric
        self.force_dtypes = force_dtypes
        self.force_col_len = force_col_len

        # Input file opening
        self.f_list = self._open_files (fn)

        # Init extra private variables
        self._previous_index = -1
        self._current_index = 0
        self._header_len = 0
        self._current = None
        self._previous = None

        # Define colname based on provided list of names
        if colnames and isinstance( colnames, (list, tuple)):
            self.colnames=colnames

        # Define colnames based on file header. Need to be the same for all the files
        elif first_line_header:
            self.colnames = []
            for fn, fp in self.f_list:
                if not self.colnames:
                    self.log.debug("Reading header from file: {}".format(fn))
                    self.colnames = self._get_first_line_header(fp)
                elif self.colnames != self._get_first_line_header(fp):
                    raise FileParserError ("Inconsistant headers between input files {}".format(fn))
            self.log.debug("Column names from header: '{}'".format(" / ".join(self.colnames)))
        else:
            raise ValueError("Invalid column name option")

        # Save initial number of columns
        self.ncols = len(self.colnames)

        # Define custom namedtuple to be returned as a line
        if include_byte_len:
            self.colnames.append("byte_len")
        self.lt = namedtuple("lt", self.colnames)

        # Set types to try to cast data in
        self.dtypes_index = self._set_types(dtypes)

    #~~~~~~~~~~~~~~MAGIC AND PROPERTY METHODS~~~~~~~~~~~~~~#

    def __len__ (self):
        size = 0
        for fn, fp in self.f_list:
            size+= int(os.path.getsize(fn))
        return size-self._header_len

    def __enter__ (self):
        return self

    def close (self):
        for fn, fp in self.f_list:
            try:
                self.log.debug ("Closing file:{}".format(fn))
                fp.close()
            except Exception as E:
                self.log.debug.warning (E)

    def __exit__(self, exception_type, exception_val, trace):
        self.close()

    def __iter__ (self):
        for i, (fn, fp) in enumerate(self.f_list):
            self.log.debug("Starting to parse file {}".format(fn))
            self._current_index = i
            for line in fp:
                self.counter["Lines Parsed"]+=1
                if line.startswith(self.comment):
                    self.counter["Comment lines skipped"]+=1
                    continue
                try:
                    line = self._parse_line(line)
                    self._previous = self._current
                    self._current = line
                    self.counter["Line successfully parsed"]+=1
                    yield line
                except (FileParserError, TypeError) as E:
                    self.counter["Malformed or Invalid Lines"]+=1
                    self.log.debug(E)
                    self.log.debug("File {}: Invalid line {}".format(fn, line))
            self.log.debug("End of file: {}".format(fn))
        self.log.debug("All files done")

    #~~~~~~~~~~~~~~PUBLIC METHODS~~~~~~~~~~~~~~#

    def current (self):
        return self._current

    def previous (self):
        return self._previous

    def next (self):
        # try to read a line
        while True:
            try:
                fn, fp = self.f_list[self._current_index]
                if self._current_index > self._previous_index:
                    self.log.debug("Starting to parse file {}".format(fn))
                    self._previous_index=self._current_index
                line = next(fp)
                if line.startswith(self.comment):
                    self.counter["Comment lines skipped"]+=1
                    continue
                line = self._parse_line(line)
                self._previous = self._current
                self._current = line
                self.counter["Line successfully parsed"]+=1
                return line

            # If one of the file is finished, start another one
            except StopIteration:
                self.log.debug("End of file: {}".format(fn))
                self._current_index+=1

            # End condition if all files where read
            except IndexError:
                self.log.debug("All files done")
                raise StopIteration

            except FileParserError:
                self.counter["Malformed or Invalid Lines"]+=1

    #~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~#

    def _get_first_line_header (self, fp):
        header_line = next(fp)
        self._header_len+=len(header_line)
        return header_line.rstrip().split(self.sep)

    def _parse_line (self, line):
        byte_len = len(line)
        line = line.rstrip().split(self.sep)

        #if the length of the line is inconsistent with the header
        if len(line) != self.ncols:
            # Raise error is required (default)
            if self.force_col_len:
                raise FileParserError("Invalid Number of fields found")
            # Else truncate extra fields
            else:
                line = line[:self.ncols]

        # Try to autocast in int or float
        if self.auto_numeric:
            for i in range(len(line)):
                val = line[i]
                try:
                    line[i] = int(val)
                except ValueError:
                    try:
                        line[i] = float(val)
                    except ValueError:
                        pass

        # Cast values according to provided types
        elif self.dtypes_index:
            for i, dtype in self.dtypes_index.items():
                try:
                    line[i] = dtype(line[i])
                except Exception:
                    if self.force_dtypes:
                        raise FileParserError("Cannot cast field in required type")

        # Add byte length if needed
        if self.include_byte_len:
            line.append(byte_len)

        # Return nametuple
        return self.lt(*line)

    def _set_types (self, dtypes):
        """"""
        dtypes_index = OrderedDict()
        if dtypes:
            for i, name in enumerate(self.colnames):
                if name in dtypes:
                    dtypes_index[i] = dtypes[name]
        return dtypes_index

    def _open_files (self, fn_list):
        """Transparently open files, lists, regex, gzipped or not"""
        f_list = []

        # Cast single file or regex to list
        if isinstance(fn_list, str):
            fn_list = [fn_list]

        if isinstance(fn_list, (list, tuple, set)):
            for fn_regex in fn_list:
                for fn in iglob(fn_regex):
                    self.counter["Input files"]+=1
                    if fn.endswith(".gz"):
                        self.log.debug("Opening file {} in gzip mode".format(fn))
                        fp = gzip.open(fn, "rt")
                    else:
                        self.log.debug("Opening file {} in normal mode".format(fn))
                        fp = open(fn, "r")
                    f_list.append((fn,fp))

            return f_list

        else:
            raise ValueError ("Invalid file type")

class FileParserError (Exception):
    """ Basic exception class for FileParserError """
    pass
