# Module to execute IRAF tasks with PyRAF

from pyraf import iraf

def run_task(name, input, prefix, params):
    """
    Input & output are file objects.
    """
    pass

    # Run the IRAF task, check it didn't barf & that the expected files
    # were created etc. and return a file object with the right prefix.

    return output

def run_imstat(input):
    iraf.images(_doprint=0)
    iraf.imstat(input)

