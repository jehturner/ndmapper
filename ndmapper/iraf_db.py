# Copyright(c) 2016 Association of Universities for Research in Astronomy, Inc.
# by James E.H. Turner.

import re

__all__ = ['add_db_entry']

def add_db_entry(db_file, key, value):
    """
    Add or replace a parameter name and value for every aperture entry in the
    named IRAF database file (eg. 'shift', 10.0).
    """

    # Divide buffer into sections on "begin", preceded by optional comment(s):
    delim = re.compile('(?:^[ \t]*#.*?\n)*(?:^begin[ \t])', flags=re.M)

    with open(db_file, 'r+') as fileobj:

        # Read the whole file so we can easily match multi-line regexes:
        buff = fileobj.read()

        # Get a list of section start/end indices (slightly slow):
        indices = [m.start() for m in delim.finditer(buff)] + [len(buff)]

        # Construct modified buffer copy:
        outbuff = ''
        for start, end in zip(iter(indices), iter(indices[1:])):

            # Get entry with any existing instances of the key removed:
            entry = re.sub('^[ \t]*{0}[ \t][^\n]*\n?'.format(key), '\n',
                           buff[start:end], flags=re.M)

            # Strip any trailing space and append the new field:
            outbuff += entry.rstrip() + '\n\t{0}\t{1}\n\n'.format(key, value)

        # Replace file contents:
        fileobj.seek(0)
        fileobj.write(outbuff)
        fileobj.truncate()

