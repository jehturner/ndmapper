import pytest
import os.path
import numpy as np
from astropy.utils import find_current_module
from astropy.nddata import NDDataArray
from ..io import NDMapIO
from ..data import FileName, DataFile, DataFileList, NDLater

# Determine path to data in common to some of the tests below:
# Using get_pkg_data_filename() from astropy here was causing problems with
# py.test expecting remote data.
#fn_mefnodq = get_pkg_data_filename('data/eqbprgS20120827S0069_flat.fits')
module = find_current_module(0, True)
module_path = os.path.dirname(module.__file__)
fn_mefnodq = os.path.join(module_path, 'data', 'eqbprgS20120827S0069_flat.fits')


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


def test_FileName_multiple_ext():

    fn = FileName('some.dir/rgS20120827S0066.fits.gz')

    assert fn.dir == 'some.dir' and \
           fn.prefix == 'rg' and \
           fn.base == 'S20120827S0066' and \
           fn.suffix == [] and \
           fn.ext == 'fits.gz' and \
           fn.standard == True and \
           str(fn) == 'some.dir/rgS20120827S0066.fits.gz'

def test_DataFile_gemini_IRAF_new_1():

    # Use a filename that can't exist, just to be sure it doesn't...
    df = DataFile('rgS20120827S9999.fits', mode='new')

    assert df.filename.dir == '' and \
           df.filename.prefix == 'rg' and \
           df.filename.base == 'S20120827S9999' and \
           df.filename.suffix == [] and \
           df.filename.ext == 'fits' and \
           df.filename.standard == True and \
           len(df) == 0


def test_DataFile_gemini_IRAF_existing_1():

    df = DataFile(filename=fn_mefnodq)

    # Ignore the directory name because the above call translates it to the
    # absolute location of this package on the system.
    assert df.filename.prefix == 'eqbprg' and \
           df.filename.base == 'S20120827S0069' and \
           df.filename.suffix == ['_flat'] and \
           df.filename.ext == 'fits' and \
           df.filename.standard == True and \
           df.meta['OBJECT'] == 'GCALflat' and \
           len(df) == 2 and \
           abs(np.mean(df[0]) - 1.0623) < 0.0001 and \
           abs(np.mean(df[1]) - 0.9461) < 0.0001


def test_DataFileList_gemini_IRAF_existing_1():

    dfl = DataFileList(fn_mefnodq)

    assert len(dfl) == 1 and \
           dfl[0].filename.prefix == 'eqbprg' and \
           dfl[0].filename.base == 'S20120827S0069' and \
           dfl[0].filename.suffix == ['_flat'] and \
           dfl[0].filename.ext == 'fits' and \
           dfl[0].filename.standard == True and \
           dfl[0].meta['OBJECT'] == 'GCALflat' and \
           len(dfl[0]) == 2 and \
           abs(np.mean(dfl[0][0]) - 1.0623) < 0.0001 and \
           abs(np.mean(dfl[0][1]) - 0.9461) < 0.0001


def test_DataFileList_replacing_data_1():

    # Replace the data & header from file with blank ones:
    dfl = DataFileList(filenames=fn_mefnodq, data=[], meta={})

    assert len(dfl) == 1 and len(dfl[0]) == 0 and len(dfl[0].meta) == 0


def test_DataFileList_len_mismatch_1():

    # Cannot have multiple data objects per filename:
    with pytest.raises(ValueError):
        dfl = DataFileList(filenames=fn_mefnodq, data=[DataFile(), DataFile()])


def test_DataFileList_broadcast_data_1():

    # But we can map the same dataset to multiple files:
    dfl = DataFileList([fn_mefnodq, 'test_name'],
                       data=DataFile(data=NDDataArray([1,2,3])),
                       mode='overwrite')

    # Produces 2 separate DataFiles referring to the same data array (as
    # long as the data aren't lazily-loaded separately).
    assert len(dfl) == 2 and dfl[0] is not dfl[1] \
        and len(dfl[0]) == len(dfl[1])\
        and [ndda.data is nddb.data for ndda, nddb in zip(dfl[0], dfl[1])]


def test_DataFileList_copy_self_1():

    dfl1 = DataFileList(fn_mefnodq)
    dfl2 = DataFileList(data=dfl1)

    assert dfl1 is not dfl2 and dfl1 == dfl2


def test_DataFileList_append_1():

    dfl = DataFileList(filenames=fn_mefnodq, mode='overwrite')
    dfl.append(DataFile('some_file', mode='new'))

    assert len(dfl) == 2 and str(dfl[0].filename) == fn_mefnodq \
        and str(dfl[1].filename) == 'some_file'


def test_DataFileList_nested_nddata_1():

    # Here the outer list maps to DataFiles and the inner one (if applicable)
    # to data extensions/groups within a file. The inner list can be omitted
    # where the file only contains a single NDData instance.
    dfl = DataFileList(filenames=['test_name_1', 'test_name_2'],
                       data=[[NDDataArray([1,2,3]), NDDataArray([4])], \
                             NDDataArray([5,6])], mode='new')

    assert len(dfl[0]) == 2 and len(dfl[1]) == 1


# Test our modified NDData init API & lazy loading behaviour:
def test_NDLater_1():

    ndd = NDLater(iomap=NDMapIO(fn_mefnodq, ident=2, data_idx=2))
    # For now, just check that the API works rather than trying to analyze
    # memory usage here (has been checked separately though).
    mean = np.mean(ndd.data)
    del ndd.data

    assert abs(mean - 0.9461) < 0.0001

