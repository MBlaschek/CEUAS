#!/usr/bin/env python3
__version__ = '0.1'
__author__ = 'MB'
__status__ = 'dev'
__date__ = 'Mon Jul 13 10:29:21 CEST 2020'
__institute__ = 'UNIVIE'
__github__ = 'git@github.com:MBlaschek/CEUAS.git'
__doc__ = """
This script is used in C3S_C311c_Lot2.2 Data Quality Assessment
QC v%s
Maintained by %s at %s
Github: %s [%s]
Updated: %s
""" % (__version__, __author__, __institute__, __github__, __status__, __date__)


def quality_report(filename: str):
    """ READ CDM data,
    0. Check Units and CDM definition (missing attributes, ...)
    1. Check datetime
    2. Check pressure data
    3. Check Variable limits
    4. Check Climatology + likelihood of outlier
    5. Check Departure Statistics

    Write results to a text file and return a summary of issues
    Args:
        filename:

    Returns:
        dict : Found problems
    """
    pass


if __name__ == '__main__':
    import sys

    for ifile in sys.argv[1:]:
        print(ifile)
        quality_report(ifile)
