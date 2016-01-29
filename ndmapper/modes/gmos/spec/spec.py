# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

from ..gmos import *

# Default calibration dependence for GMOS spectroscopy. The 'spectwilight'
# type is left out for now, as it's currently not used here.
CAL_DEPS = {'target' : ['specphot', 'flat', 'arc', 'bias'],
            'specphot' : ['flat', 'arc', 'bias'],
            'flat' : ['arc', 'bias'],
            'arc' : ['bias'],
            'bias' : []
           }

# For convenience & avoidance of visual noise, initialize these dictionaries
# that can be used for managing calibrations within user scripts:
biases, traces, arcs, flats, standards = {}, {}, {}, {}, {}

