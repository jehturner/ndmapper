NDMapper change log
===================

0.3 (unreleased)
----------------

- No changes yet. Planning to split out NDProcess module.


0.2 (2016-10-06)
----------------

- Add wrapper steps for ``gmos.spec.ifu`` data cube mosaicking and logarithmic
  re-binning, using my PyFU Python package.

- Add ``cosmetics.clean_cosmic_rays`` step, based on McCully's lacosmicx (fast
  Cython implementation of LA Cosmic algorithm; currently requires forked
  version at https://github.com/jehturner). Also supporting code for fitting.

- Add ``gemini.clean_pixels`` wrapper for gemfix (& work around IRAF corrupting
  'MDF' extensions in the latter).

- Add ``gmos.spec.normalize_QE`` wrapper for gqecorr (to correct differences in
  QE variation with wavelength between CCD chips).

- Add ``gmos.spec.ifu.subtract_bg`` wrapper for gfscatsub (scattered light
  subtraction) and associated ``gmos.spec.ifu.background_regions`` wrapper for
  gffindblocks (to identify gaps between blocks of fibres). Generalize the
  calibration association slightly to support the latter.

- Pass through additional options in ``gmos.spec.ifu.extract_spectra``, to help
  control the results.

- Allow associating MDFs as calibrations (so they can differ between inputs
  when the numbers of detectable fibres differ).

- Automatically map 'LTT' standard names to IRAF's onedstds 'l' convention,
  when finding the appropriate look-up table based on ``meta`` in ``gmos.spec``.

- Add a ``shift_spectra`` step for ``gmos`` that can apply a zero-point
  correction to the IRAF wavelength database(s) associated with a DataFileList
  (with the shift values also optionally associated as calibrations).

- New ``reload`` & ``save`` convenience methods in ``DataFileList``.

- Add initial arithmetic operators to ``NDLater`` and ``DataFile`` (with
  support for specifying an output ``filename`` when called as methods). In
  the process, add some missing NDDataArray attributes.

- Support assignment to a ``DataFile`` index (__setitem__).

- ``NDLater`` compatibility with AstroPy 1.2.1.

- Support assignment to ``TabMapIO.table`` (for updating MDFs; still missing a
  proper public API).

- Change filename prefix for ``cosmetics.add_bpm`` from 'b' to 'k', to avoid
  (unexpected) confusion with scattered light subtraction.


0.1.2 (2016-03-26)
------------------

- Fix ``libutils.splitext`` to ignore the path when splitting on the first dot,
  in case directory names contain dots (specifically temp directories created
  by py.test).

- Include FITS data in MANIFEST.in & setup.py so the tests can be run from an
  installed copy of the package (in particular one created by py.test).

- Configure py.test not to capture stdout (which causes importing pyraf to fail
  due to incompatible redirection) and to skip a doctest with no data, avoiding
  failures when using py.test instead of nose.

- Determine the path to FITS data in the tests without using astropy
  ``get_pkg_data_filename``, to avoid remote data errors from py.test.

- Enabled testing with Travis.

- Initial Python 3 compatibility changes, mainly for the `data` submodule.


0.1.1 (2016-03-12)
------------------

- Clarify that processing library ``outputs`` parameters are now intended only
  to specify filenames (rather than supporting shared data, at least for now)
  and rename them accordingly, to ``out_names``.

- Fixed a run_task failure when a DataFile modified in memory was passed as
  the input to multiple IRAF steps, which raised an ``IOError`` second time,
  due to sharing identically the same NDLater instances (& TabMapIO list) when
  instatiating one DataFile from another, such that run_task ended up
  modifying the NDMapIO instances of its inputs with a temporary filename.

- Include missing ``add_bpm`` in ``cosmetics`` name space.


0.1 (2016-03-01)
----------------

- PyPI test release.

