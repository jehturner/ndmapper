# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

from ..gmos import *

# Default calibration dependence for GMOS spectroscopy:
CAL_DEPS = {'target' : ['specphot', 'flat', 'bias'],
            'specphot' : ['flat', 'bias'],
            'flat' : ['bias'],
            'bias' : []
           }

