.. _ndmapper_process:


Data processing steps
*********************

Initially, most of these steps wrap existing tasks in the Gemini IRAF package
(reimplementing their APIs) as if they were native Python functions; these and
the steps implemented in Python both operate on DataFileList objects and can be
mixed & matched relatively transparently (to the extent practical when using
IRAF databases etc.; see their individual docstrings).


Reference/API
=============

.. automodapi:: ndmapper.lib.cosmetics
    :no-inheritance-diagram:

.. automodapi:: ndmapper.lib.gemini
    :no-inheritance-diagram:

.. automodapi:: ndmapper.lib.gmos
    :no-inheritance-diagram:

.. automodapi:: ndmapper.lib.gmos.spec
    :no-inheritance-diagram:

.. automodapi:: ndmapper.lib.gmos.spec.ifu
    :no-inheritance-diagram:

Entries without links are inherited from their parent module, in which the full
documentation can be found.

