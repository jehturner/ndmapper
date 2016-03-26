========
NDMapper
========

A module for constructing intelligible astronomical data reduction processes,
based on a data abstraction that maps files on disk to collections of AstroPy
NDData objects.

Currently, this module also includes a small library of data processing steps
for GMOS data, a wrapper for integrating IRAF tasks into the processing scheme
and some simple functionality for managing calibration associations; these
things are intended to be moved to a separate "ndprocess" module as this
prototype nears completion.

NDMapper is built on the AstroPy affiliated package template but has no
official association with AstroPy at present.

.. image:: https://travis-ci.org/jehturner/ndmapper.svg
    :target: https://travis-ci.org/jehturner/ndmapper

