"""
I/O routines used internally by NDMapper. These have a public API for the
convenience of anyone who wants to use their file format abstraction directly
but users should nearly always work with the higher-level DataFile interface
instead (which is also a bit less liable to change).
"""

from .io import *

