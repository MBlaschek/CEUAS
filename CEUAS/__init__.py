
__version__ = '0.1'
__author__ = 'LH,MB,FA'
__status__ = 'dev'
__date__ = 'Di Feb 26 11:07:48 CET 2019'
__institute__ = 'Univie, IMGW'
__github__ = 'git@github.com:MBlaschek/CEUAS.git'
__doc__ = """
Radiosonde Tools Collection v%s
Maintained by %s at %s
Github: %s [%s]
Updated: %s
""" % (__version__, __author__, __institute__, __github__, __status__, __date__)

from . import data
from . import meta
from . import utils
from . import adjust
