NDMapper change log
===================

0.1.1 (unreleased)
------------------

- Nothing changed yet.


0.1 (2015-03-01)
----------------

- PyPI test release.

  Still pending a fix for a run_task failure when a DataFile modified in memory
  is passed as the input to multiple steps, which occurs the second time due to
  DataFiles instantiated from another instance sharing identically the same
  NDLater instances, such that run_task ends up modifying the iomap instances
  (for lazy loading) of its inputs with the temporary filename.

