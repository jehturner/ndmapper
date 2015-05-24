import os.path
import numpy as np
from astropy.utils.data import get_pkg_data_filename
from ..data import FileName, DataFile, DataFileList


def test_FileName_gemini_IRAF_1():

    fn = FileName('/home/james/rgS20120827S0066.fits')

    assert fn.dir == '/home/james' and \
           fn.prefix == 'rg' and \
           fn.base == 'S20120827S0066' and \
           fn.suffix == [] and \
           fn.ext == 'fits' and \
           fn.standard == True and \
           str(fn) == '/home/james/rgS20120827S0066.fits'

def test_FileName_nonconforming_1():

    fn = FileName('/home/fred/blah.fits')

    assert fn.dir == '/home/fred' and \
           fn.prefix == '' and \
           fn.base == 'blah' and \
           fn.suffix == [] and \
           fn.ext == 'fits' and \
           fn.standard == False and \
           str(fn) == '/home/fred/blah.fits'

def test_DataFile_gemini_IRAF_new_1():

    # Use a filename that can't exist, just to be sure it doesn't...
    df = DataFile(filename='rgS20120827S9999.fits')

    assert df.filename.dir == '' and \
           df.filename.prefix == 'rg' and \
           df.filename.base == 'S20120827S9999' and \
           df.filename.suffix == [] and \
           df.filename.ext == 'fits' and \
           df.filename.standard == True and \
           len(df) == 0

def test_DataFile_gemini_IRAF_existing():

    fn = get_pkg_data_filename('data/eqbprgS20120827S0069_flat.fits')
    df = DataFile(filename=fn)

    print df.filename.dir

    # Ignore the directory name because the above call translates it to the
    # absolute location of this package on the system.
    assert df.filename.prefix == 'eqbprg' and \
           df.filename.base == 'S20120827S0069' and \
           df.filename.suffix == ['_flat'] and \
           df.filename.ext == 'fits' and \
           df.filename.standard == True and \
           len(df) == 2 and \
           abs(np.mean(df[0]) - 1.0623) < 0.0002 and \
           abs(np.mean(df[1]) - 0.9461) < 0.0002


