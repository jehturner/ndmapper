# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

# Define supported file formats, mapped to the corresponding file extension
# names. For each <key> in the dictionary below, there should also exist a
# module ndmapper.io._<key>.py, containing the back-end loader functions for
# that file format, matching the generic function defs in io.py.

# If more than one back end happens to support the same file format, the order
# of precedence is that written here:
formats = {'fits' : ('fit', 'fits', 'fits.gz', 'fits.bz2')}

