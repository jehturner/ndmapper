# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

import pyfits

def load_common_meta(filename):
    filename = str(filename)
    return pyfits.getheader(filename)

