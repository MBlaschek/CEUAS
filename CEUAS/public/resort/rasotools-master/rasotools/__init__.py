__version__ = '0.4'
__author__ = 'MB'
__status__ = 'dev'
__date__ = 'Mit Nov 20 08:36:56 CET 2019'
__institute__ = 'Univie, IMGW'
__github__ = 'git@github.com:MBlaschek/rasotools.git'
__doc__ = """
Radiosonde Tools Collection v%s
Maintained by %s at %s
Github: %s [%s]
Updated: %s
""" % (__version__, __author__, __institute__, __github__, __status__, __date__)

import os as _os

from . import fun
from . import cls
from . import net
from . import met
from . import grid
from . import bp
from . import plot


def _getlibs():
    """Get library version information for printing"""
    import numpy
    import pandas
    import xarray

    return __version__, numpy.__version__, pandas.__version__, xarray.__version__


stations = fun.station.read_igrav2_stationlist(None)

config = cls.Bunch()
config.homedir = _os.getenv("HOME")
config.wkdir = _os.getcwd()
config.igradir = ''
config.marsdir = ''
config.outdir = _os.getcwd() + '/results'
config.rasodir = _os.getcwd() + '/archive'
config.std_plevels = [1000., 2000., 3000., 5000., 7000., 10000., 15000., 20000., 25000., 30000., 40000., 50000., 70000.,
                      85000., 92500., 100000.]
config.era_plevels = [1000., 2000., 3000., 5000., 7000., 10000., 12500., 15000., 17500., 20000., 22500., 25000., 30000.,
                      35000., 40000., 45000., 50000., 55000., 60000., 65000., 70000., 75000., 77500., 80000., 82500.,
                      85000., 87500., 90000., 92500., 95000., 97500., 100000.]

config.rttov_profile_limits = None
config.libinfo = "RT(%s) NP(%s) PD(%s) XR(%s)" % _getlibs()


def load_config(filename='rasoconfig.json'):
    """ load config from JSON file

    Parameters
    ----------
    filename : str
        Filename of config file

    Example File
    ------------
    {
    "marsdir": "",
    "igradir": "",
    "rasodir": "/home/mblaschek/workspace/raso_archive"
    }
    """
    import os
    import json

    if os.path.isfile(filename):
        with open(filename, 'r') as file:
            variables = json.loads(file.read())

        for ikey, ival in variables.items():
            setattr(config, ikey, ival)
            print("[CONFIG] ", ikey, ":", repr(ival))
        print("[CONFIG] ", filename)


def dump_config(filename='rasoconfig.json'):
    """ Write config to file

    Parameters
    ----------
    filename : str
        Filename to save config
    """
    import json
    variables = vars(config)
    with open(filename, 'w') as file:
        file.write(json.dumps(variables, indent=0))

    print("Configuration written: ", filename)


def open_radiosonde(name, ident=None, filename=None, directory=None, **kwargs):
    """ Create a Radiosonde object from opening a dataset

    Args:
        name (str): used as filename and/or as dataset name
        ident (str): radiosonde wmo or igra id
        filename (str): filename to read from (netcdf)
        directory (str): directory of radiosonde store, default config rasodir

    Returns:
        Radiosonde : Radiosonde class object
    """
    from .cls import Radiosonde
    from .fun import get_data
    if name == 'example':
        ident = 'AUM00011035'
        name = 'IGRAv2'
        filename = get_data('AUM00011035_XARRAY.nc')

    if ident is None:
        ident = "Unknown"

    out = Radiosonde(ident)
    out.add(name, filename=filename, directory=directory, **kwargs)
    return out


def open_network(directory, **kwargs):
    """ Load pickle dump of a network
    
    Args:
        directory (str): directory or filename to pickle
    
    Returns:
        Network : Radiosonde network object
    """
    import pickle
    import os
    if os.path.isdir(directory):
        for ifile in sorted(os.listdir(directory)):
            if os.path.isfile(directory + '/' + ifile) and 'network_' in ifile:
                print("Restoring:", directory + '/' + ifile)
                return pickle.load(open(directory + '/' + ifile, 'rb'))
            
    elif os.path.isfile(directory):
        print("Restoring:", directory)
        return pickle.load(open(directory, 'rb'))
    else:
        pass
