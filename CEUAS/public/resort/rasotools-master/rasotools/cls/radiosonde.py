# -*- coding: utf-8 -*-
import os

import xarray as xr

from .bunch import Bunch
from ..fun import message, update_kw, dict_add, print_fixed, netcdf, levelup

__all__ = ["Radiosonde"]


class Radiosonde(object):
    def __init__(self, ident, data=None, directory=None, **kwargs):
        self.ident = ident
        if data is not None:
            self.data = Bunch(**data)
        else:
            self.data = Bunch()  # empty
        self.attrs = Bunch(**kwargs)
        self.directory = directory

    def __repr__(self):
        summary = u"Radiosonde (%s)\n" % self.ident
        if len(list(self.data)) > 0:
            summary += "Data: \n"
            summary += repr(self.data)
        if len(list(self.attrs)) > 0:
            summary += "\nGlobal Attributes: \n"
            summary += repr(self.attrs)
        if self.directory is not None:
            summary += "\nDirectory: " + self.directory
        return summary

    def __getitem__(self, item):
        if item in self.data:
            return getattr(self.data, item)
        return None

    def __setitem__(self, key, value):
        if key in self.data:
            print("Warning overwriting ...", key)
        setattr(self.data, key, value)

    def __delitem__(self, key):
        if key in self.data:
            delattr(self.data, key)

    def add(self, name, filename=None, directory=None, rename=None, cfunits=False, close=True, xwargs=None, **kwargs):
        """ Add dataset to radiosonde class object [container]

        Args:
            name (str): used as filename and/or as dataset name
            filename (str): filename to read from (netcdf)
            directory (str): directory of radiosonde store, default config rasodir
            rename (dict): rename variables
            cfunits (bool): apply cfunits
            xwargs (dict): xarray open_dataset keywords
        """
        if rename is None:
            rename = {}

        if xwargs is None:
            xwargs = {}

        from .. import config

        if self.directory is not None:
            directory = self.directory

        if 'engine' not in xwargs.keys():
            xwargs['engine'] = 'h5netcdf'

        if directory is None:
            if 'rasodir' in config and os.path.isdir(config.rasodir):
                directory = config.rasodir + '/' + str(self.ident) + '/'
            else:
                directory = './' + str(self.ident) + '/'
            message(directory, **kwargs)

        if filename is not None:
            if '.nc' not in filename:
                print("Warning can read only NetCDF / Xarray open_dataset formats")

            if '*' not in filename:
                if not os.path.isfile(filename):
                    if directory is not None:
                        if not os.path.isfile(directory + '/' + filename):
                            filename = directory + '/' + str(self.ident) + '/' + filename
                        else:
                            filename = directory + '/' + filename

        if filename is None:
            if os.path.isfile(directory + '/' + name + '.nc'):
                filename = directory + '/' + name + '.nc'
            else:
                message("Not found:", name, directory, **kwargs)
                return

        if '*' in filename:
            if 'combine' not in xwargs.keys():
                xwargs['combine'] = 'by_coords'

            try:
                message("Reading ...", filename, **kwargs)
                self.data[name] = xr.open_mfdataset(filename, **xwargs)
            except OSError:
                message("Reading ...", directory + '/' + filename, **kwargs)
                self.data[name] = xr.open_mfdataset(directory + '/' + filename, **xwargs)

        else:
            message("Reading ...", filename, **kwargs)
            self.data[name] = xr.open_dataset(filename, **xwargs)

        if close:
            # add this to make sure that a file is read and closed (some error with h5netcdf, h5py)
            if hasattr(self.data[name], 'close'):
                self.data[name].load()
                self.data[name].close()

        self.data[name] = self.data[name].rename(rename)

        if cfunits:
            import cf2cdm
            for i, j in self.data[name].coords.items():
                if str(j.dtype) == 'datetime64[ns]':
                    if 'standard_name' in j.attrs:
                        del j.attrs['standard_name']
            self.data[name] = cf2cdm.translate_coords(self.data[name], cf2cdm.CDS)

        # merge dictionaries and append common
        self.attrs.__dict__ = dict_add(vars(self.attrs), dict(self.data[name].attrs))
        if 'ident' in self.attrs:
            if self.ident != self.attrs['ident']:
                message("Warning different idents: ", self.ident, ">", self.attrs['ident'],
                        **update_kw('level', -1, **kwargs))
            if self.ident == "":
                self.ident = self.attrs['ident']

    def clear(self, exclude=None, **kwargs):
        if exclude is None:
            exclude = []
        else:
            if not isinstance(exclude, list):
                exclude = [exclude]

        for idata in list(self.data):
            if idata not in exclude:
                if hasattr(self.data[idata], 'close'):
                    self.data[idata].close()  # for xarray netcdf open files
                    message(idata, 'closed', **update_kw('mname', self.ident, **kwargs))
                del self.data[idata]
                message(idata, **update_kw('mname', self.ident, **kwargs))

        for iatt in list(self.attrs):
            if iatt not in exclude:
                del self.attrs[iatt]

    def rename(self, old_name, new_name, **kwargs):
        if old_name in self.data:
            self.data.__dict__[new_name] = self.data.__dict__.pop(old_name)
            message(old_name, " > ", new_name, **kwargs)

    def list_store(self, pattern=None, directory=None, varinfo=False, ncinfo=False):
        import time
        from numpy import ceil
        from .. import config

        if self.directory is not None:
            directory = self.directory

        if directory is None:
            directory = config.rasodir
            directory += "/%s/" % self.ident

        if os.path.isdir(directory):
            summe = 0
            print("Available Files (Datasets): ")
            print(directory)
            print('_' * 80)
            for ifile in sorted(os.listdir(directory)):
                if os.path.isdir(directory + ifile):
                    continue
                if pattern is not None:
                    if pattern not in ifile:
                        continue

                current = os.path.getsize(directory + ifile) / 1024 / 1024
                itime = time.ctime(os.path.getctime(directory + ifile))
                print("%-40s : %4d MB  : %s" % (ifile, ceil(current), itime))
                if '.nc' in ifile:
                    if varinfo:
                        with xr.open_dataset(directory + ifile) as f:
                            print("Variables:")
                            print(print_fixed(list(f.data_vars), ',', 80, offset=10))
                    elif ncinfo:
                        netcdf.view(directory + ifile)
                    else:
                        pass
                # print('_' * 80)
                summe += current

            print("\n%40s : %4d MB" % ('Total', summe))
        else:
            print("Store not found!")

    def in_store(self, name, directory=None, **kwargs):
        from .. import config

        if self.directory is not None:
            directory = self.directory

        if directory is None:
            directory = config.rasodir
            directory += "/%s/" % self.ident

        if os.path.isdir(directory):
            for ifile in os.listdir(directory):
                if os.path.isdir(directory + ifile):
                    continue

                if name in directory + ifile:
                    message("Found: ", name, directory + ifile, **kwargs)
                    return True
                message(name, "!=", ifile, **levelup(**kwargs))
        return False

    def to_netcdf(self, name, filename=None, directory=None, force=False, add_global=True, xwargs=None, **kwargs):
        """Write each dataset variable to NetCDF 4

        Args:
            name (str): dataset name
            filename (str): filename
            directory (str): directory
            force (bool): force new file
            add_global (bool):
            xwargs (dict): xarray to_netcdf Options

        Result:
            bool : Status
        """
        if xwargs is None:
            xwargs = {}

        from .. import config

        if not isinstance(name, list):
            name = [name]

        if 'engine' not in xwargs.keys():
            xwargs['engine'] = 'h5netcdf'
            xwargs['format'] = 'netcdf4'

        if directory is None:
            if self.directory is not None:
                directory = self.directory

            elif 'rasodir' in config:
                if os.path.isdir(config.rasodir):
                    directory = config.rasodir + '/' + str(self.ident) + '/'
            else:
                directory = './' + str(self.ident) + '/'

        message("to dir:", directory, **kwargs)
        attrs = vars(self.attrs)

        for iname in name:
            if iname in self.data:
                if isinstance(self.data[iname], (xr.DataArray, xr.Dataset)):

                    if filename is None:
                        ifilename = directory + '/' + iname + '.nc'
                    else:
                        ifilename = filename

                    # Create directory if necessary
                    if not os.path.isdir(os.path.dirname(ifilename)):
                        idir = os.path.dirname(ifilename)
                        message("makedir: ", idir, **kwargs)
                        os.makedirs(idir, exist_ok=True)  # no Errors

                    message("Writing", ifilename, **kwargs)

                    iobj = getattr(self.data, iname)
                    if len(iobj.attrs) or add_global:
                        iobj.attrs.update(attrs)  # add global attributes

                    for i in list(iobj.attrs.keys()):
                        if iobj.attrs[i] is None:
                            del iobj.attrs[i]

                    if 'engine' in xwargs.keys():
                        if xwargs['engine'] == 'h5netcdf':
                            xwargs['encoding'] = {i: {'compression': 'gzip', 'compression_opts': 9} for i in
                                                  list(iobj.data_vars)}

                    if force:
                        xwargs.update({'mode': 'w'})
                        force = False
                    else:
                        if os.path.isfile(ifilename):
                            xwargs.update({'mode': 'a'})
                        else:
                            xwargs.update({'mode': 'w'})
                    iobj.to_netcdf(ifilename, **xwargs)

                else:
                    message("Error unknown data type (xarray) ", iname, type(self.data[iname]), **kwargs)

    def map(self, **kwargs):
        from ..plot.map import station
        lon = None
        lat = None
        for iatt in self.attrs:
            if 'station_lon' in iatt:
                lon = self.attrs[iatt]
                if isinstance(lon, str):
                    for i in lon.split(' '):
                        try:
                            lon = float(i)
                        except:
                            pass

            if 'station_lat' in iatt:
                lat = self.attrs[iatt]
                if isinstance(lat, str):
                    for i in lat.split(' '):
                        try:
                            lat = float(i)
                        except:
                            pass
        print(lon, lat)
        if lon is None or lat is None:
            raise ValueError("No location information found ??? (station_lon, station_lat")
        return station(lon, lat, **kwargs)
