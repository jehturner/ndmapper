NDMapper change log
===================


0.1.2 (unreleased)
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

