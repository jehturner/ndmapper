NDMapper change log
===================


0.1.2 (unreleased)
------------------

- Nothing changed yet.


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

