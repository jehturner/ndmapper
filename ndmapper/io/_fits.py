# Copyright(c) 2015 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

import os
import astropy.io.fits as pyfits

from .. import config
from .mapio import NDMapIO, TabMapIO


def load_common_meta(filename):

    return pyfits.getheader(filename)


def load_array_meta(filename, index):

    # This should probably have an optional header
    return pyfits.getheader(filename, index)


def load_array(filename, index):

    # Treat any int (flags) array as unsigned for the appropriate BZERO/BSCALE
    # (to avoid scaling int16 DQ to float32).
    return pyfits.getdata(filename, index, uint=True)


def save_array(filename, index, data, meta=None):

    # Convert the inputs to a PyFITS HDU:
    hdu = pyfits.ImageHDU(data=data, header=_convert_meta(meta), uint=True)

    # Open & parse the existing file:
    hdulist = pyfits.open(filename, mode='update', memmap=True, uint=True)

    # Overwrite or append our new HDU at the specified index, producing an
    # error if it doesn't already exist and isn't the next one.
    narr = len(hdulist)
    if index == narr:
        hdulist.append(hdu)
    elif index < narr:
        hdulist[index] = hdu
    else:
        raise IndexError('index %d is out of range for %s' % \
                         (index, filename))

    hdulist.close()


def load_table_meta(filename, index):
    return pyfits.getheader(filename, index)


def load_table(filename, index):
    # Treat any int (flags) array as unsigned for the appropriate BZERO/BSCALE
    # (to avoid scaling int16 DQ to float32). Not sure whether this applies to
    # binary tables but it should be harmless.
    return pyfits.getdata(filename, index, uint=True)


def save_list(filename, data, array_meta, identifiers, types, common_meta):

    narr = len(data)

    if array_meta is None:
        array_meta = [None for arr in data]
    if identifiers is None:
        identifiers = [(None, None) for arr in data]
    if types is None:
        types = ['image' for arr in data]
    else:
        if not all([arrtype in ['image', 'table'] for arrtype in types]):
            raise ValueError('types values must be \'image\' or \'table\'')

    if len(array_meta) != narr or len(identifiers) != narr \
                               or len(types) != narr:
        raise ValueError('lengths of array_meta, identifiers & types must '\
                         'match data')

    phu = pyfits.PrimaryHDU(header=_convert_meta(common_meta))
    phu.header['EXTEND'] = True  # required when adding MEF extensions

    exists = os.path.exists(filename)
    mode = 'update' if exists else 'append'

    hdulist = pyfits.open(filename, mode=mode, uint=True)

    oldlen = len(hdulist)-1

    if oldlen == -1:
        hdulist.append(phu)  # file either new or just empty/corrupt (no PHU)
    else:
        hdulist[0] = phu

    # Loop over the image extensions/inputs:
    n = 0
    for n, (arr, meta, (name, ver), arrtype) in \
        enumerate(zip(data, array_meta, identifiers, types), start=1):

        # Update only those extensions for which data and/or a header have
        # been provided; otherwise it's understood that the caller wants to
        # re-use whatever is already in the file (if applicable), to minimize
        # unnecessary writes (which io.fits does automatically).
        if arr is not None or meta is not None or n > oldlen:

            hduclass = pyfits.BinTableHDU if arrtype is 'table' \
                       else pyfits.ImageHDU
            hdu = hduclass(data=arr, header=_convert_meta(meta), name=name,
                           uint=True)
            # Set INHERIT = F here?
            hdu.ver = -1 if ver is None else ver

            if n <= oldlen:
                hdulist[n] = hdu
            else:
                hdulist.append(hdu)

    # Remove any unused extensions remaining at the end of an existing file
    # (unless explicitly given None values for those indices):
    del hdulist[n+1:]

    hdulist.close()


def _convert_meta(meta):
    if meta:
        if isinstance(meta, pyfits.Header):  # should use hasttr 'cards'?
            hdr = meta  # preserve comments
        else:
            hdr = pyfits.Header([pyfits.Card(key, val, verify='warn') \
                                 for key, val in meta.iteritems()])
        return hdr
    # (Return None by default if meta is None)


def map_file(filename, labels):

    hdulist = pyfits.open(filename, mode='readonly')

    # Use package default extname labels if not specified:
    if not labels:
        labels = config['labels']

    # A dict of empty lists to sort recognized extensions into:
    idx_dict = {'data' : [], 'uncertainty' : [], 'flags' : [], 'undef' : [],
                'tables' : []}

    # Sort any FITS image extensions by EXTNAME into SCI/VAR/DQ lists
    # (remembering the original MEF index for later I/O):
    have_names = False
    for idx, hdu in enumerate(hdulist):

        # Classify image extensions by extname (NB. any data in a FITS primary
        # header must be an image array according to the std & PyFITS docs):
        if (isinstance(hdu, pyfits.ImageHDU) or idx==0) and hdu.size > 0:

            # The name/ver attributes are semi-documented but seem to be
            # part of the public API now. The name seems to default to ''
            # when undefined but the following would also work for None:
            if hdu.name and idx > 0:  # ignore 'PRIMARY'
                have_names = True

            if hdu.name == labels['data']:
                idx_dict['data'].append(idx)

            elif hdu.name == labels['uncertainty']:
                idx_dict['uncertainty'].append(idx)

            elif hdu.name == labels['flags']:
                idx_dict['flags'].append(idx)

            elif not hdu.name or idx == 0:  # ignore 'PRIMARY'
                idx_dict['undef'].append(idx)

            # else:
            #     ignore any extensions with unrecognized names

            # print(idx_dict)

        # Also collect any table extensions for separate propagation. Any
        # more exotic extensions are currently ignored but should really also
        # be propagated in extras somehow.
        elif isinstance(hdu, (pyfits.BinTableHDU, pyfits.TableHDU)):
            idx_dict['tables'].append(idx)

    # If there are no named image extensions, treat the unnamed ones as our
    # "data" (SCI) arrays (otherwise ignore them):
    if not have_names:
        idx_dict['data'] = idx_dict['undef']

    # List existing data (SCI) array extension numbers for reference:
    extvers = [hdulist[idx].ver for idx in idx_dict['data']]

    # Create the NDMapIO instances for NDLater. Since the main data array is
    # mandatory, we just ignore any uncertainty (VAR) or flags (DQ) extensions
    # without a corresponding data (SCI) array and loop over the latter:
    lastver = 0
    maplist = []
    for data_idx in idx_dict['data']:

        data_hdu = hdulist[data_idx]

        # Give any unnumbered SCI extensions the next available EXTVER after
        # the last one used:
        if data_hdu.ver < 1:  # seems to default to -1; also works for None
            thisver = lastver + 1
            while thisver in extvers:
                thisver += 1
            data_hdu.ver = thisver
            uncert_idx = None
            flags_idx = None

        # Otherwise, if the EXTVER was defined to begin with, look for
        # associated uncertainty & flags (VAR/DQ) extensions:
        else:
            # Find uncertainty & flags HDUs matching this EXTVER:
            uncert_idx = None
            for idx in idx_dict['uncertainty']:
                hdu = hdulist[idx]
                if hdu.ver == data_hdu.ver:
                    uncert_idx = idx
                    break

            flags_idx = None
            for idx in idx_dict['flags']:
                hdu = hdulist[idx]
                if hdu.ver == data_hdu.ver:
                    flags_idx = idx
                    break

        lastver = data_hdu.ver

        # Instantiate the NDMapIO instance, recording the original FITS
        # extension indices and the group extver (== data extver).
        maplist.append(NDMapIO(filename, ident=data_hdu.ver,
            data_idx=data_idx, uncertainty_idx=uncert_idx,
            flags_idx=flags_idx))

    # Use another map list to propagate & access any extra extensions
    # (currently just tables) separately from the main array groups:
    tables = [TabMapIO(filename, tab_idx, label=hdulist[tab_idx].name,
              ident=hdulist[tab_idx].ver) for tab_idx in idx_dict['tables']]

    # We don't keep the file open continually, since it may get updated later
    # by IRAF or whatever (this means some trickery later on to keep io.fits
    # happy, since we haven't read in the lazy-loaded data arrays yet).
    hdulist.close()

    return maplist, tables

