# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from . import Bunch

__all__ = ['Network']


class Network(object):

    def __init__(self, name, idents=None):
        if not isinstance(name, str):
            raise ValueError("Name needs to be s string")

        self.name = name
        if idents is not None:
            if not isinstance(idents, list):
                idents = [idents]
            self.idents = idents
        else:
            self.idents = list()

        self.active = np.array([True] * len(self.idents))
        self.lon = None
        self.lat = None
        self.distance = None
        self.data = Bunch()
        self.stations = pd.DataFrame()
        self.files = Bunch()

    def __repr__(self):
        head = "%s (%d / %d)" % (self.name, sum(self.active), len(self.idents))
        body = ""
        if self.stations is not None:
            body = "\n" + self.stations[self.active].to_string(line_width=None, max_rows=10) + "\n"

        foot = "Stations: {}\n".format(len(self.idents))
        foot += "Active: {} \ {} \n".format(sum(self.active), len(self.active))
        foot += "Dist: {} \n".format(self.shape if self.distance is not None else False)
        foot += "Data: {} \n".format(list(self.data))
        return head + body + foot

    def __iter__(self):
        # iter(self.__dict__)   # does the same
        return (x for x in self.idents[self.active])

    def __getitem__(self, item):
        if item in self.idents:
            return self.stations.loc[item]
        return None

    def __setattr__(self, key, value):
        if key in ['lon', 'lat', 'stations', 'idents']:
            if key in ['stations']:
                assert type(value) == pd.DataFrame, 'requires a pandas DataFrame [ident, lon, lat]'
            else:
                value = np.asarray(value)
            if hasattr(self, 'idents'):
                self.__dict__['active'] = np.array([True] * len(self.idents))
        self.__dict__[key] = value

    # def __delitem__(self, key):
    #     if key in self.__dict__.keys():
    #         delattr(self, key)

    def _calculate_dist_matrix(self):
        from ..fun.cal import distance
        if self.lon is not None and self.lat is not None:
            print("Calculating Distance Matrix ...")
            result = [distance(self.lon, self.lat, ilon, ilat) for ilon, ilat in zip(self.lon, self.lat)]
            self.distance = pd.DataFrame(result, index=self.idents, columns=self.idents)
        else:
            raise RuntimeError("Network is missing lon, lat information")

    def neighboring_stations(self, ident, distance):
        if self.distance is None:
            self._calculate_dist_matrix()

        if ident not in self.idents:
            raise ValueError("Ident not found", ident)

        result = self.distance.loc[ident].sort_values()
        return pd.Series(result[result <= distance], name='dist_%s_km' % ident)

    def add_station(self, ident, lon, lat, **kwargs):
        pass

    def set_station(self, ident, lon, lat, **kwargs):
        pass

    def del_station(self, ident):
        pass

    def find_data_files(self, name, pattern=None, directory=None, **kwargs):
        from ..fun import find_files, message, bool_verbose
        from .. import config

        if directory is None:
            directory = config.rasodir

        files = find_files(directory, pattern)
        message("Using:", directory, pattern, "Found:", len(files), **kwargs)
        status = {}
        for ifile in files[:]:
            inside = False
            for ident in self.idents:
                if ident in ifile:
                    inside = True
                    status[ident] = ifile
                    break
            if not inside:
                files.remove(ifile)
        self.files[name] = status
        message("Summary: ", len(self.files[name]), '/', len(self.idents), **kwargs)
        if not all(self.files[name].status):
            message("Failed:", **kwargs)
            if bool_verbose(**kwargs):
                print(self.files[name].iloc[~self.files[name].status])

    # load_data(filename, name from network directory)
    # rename to load_data_from_archive
    def load_data(self, name, files=None, pattern=None, directory=None, variables=None, invert_selection=False, npp=10,
                  combine=True, **kwargs):
        """ Read and load data for Radiosondes

        Args:
            name (str): store name in self.data[name]
            files (str, list): files to read
            pattern (str): search pattern
            directory (str): directory path
            variables (list): variables to select from data
            invert_selection (bool): invert variable selection
            npp (int): multiprocessing Pool size
            combine (bool): combine output to one large xarray ?
            **kwargs:

        Returns:
            self.data[name] : dict or Dataset
        """
        from ..fun.mp import make_process_list, execute_process
        from ..net.process import read_add_process, align_concat_fill
        from ..fun import levelup, message
        from .. import config

        if variables is not None:
            if not isinstance(variables, list):
                variables = [variables]

        if directory is None:
            directory = config.rasodir

        if files is None:
            self.find_data_files(name, pattern=pattern, directory=directory, **kwargs)
            # files = find_files(directory, pattern)
            # message("Using:", directory, pattern, "Found:", len(files), **kwargs)
            # for ifile in files[:]:
            #     inside = False
            #     for ident in self.idents[self.active]:
            #         if ident in ifile:
            #             inside = True
            #             break
            #     if not inside:
            #         files.remove(ifile)
            # message("Selecting only: ", len(files), "for current idents", **kwargs)
            files = list(self.files[name].values())

        if len(files) > 0:
            message("Setting up Process ...", **kwargs)
            stuff = make_process_list(files, read_add_process,
                                      args=(variables,),
                                      invert=invert_selection,
                                      **levelup(**kwargs))

            message("Running ...", **kwargs)
            errors, data = execute_process(stuff,
                                           npp=npp,
                                           return_value=True,
                                           loop='loop' in kwargs.keys())
            message("Errors:", errors, **kwargs)
            if combine:
                message("Combining ...", **kwargs)
                #
                # Store in class
                #
                self.data[name] = align_concat_fill(*data, dim='sonde')
                #
                # remove some individual attributes
                #
                for i in list(self.data[name].attrs.keys()):
                    if 'station' in i:
                        del self.data[name].attrs[i]
                self.data[name].attrs['title'] = 'RNetwork({}, {}) - {}'.format(self.name,
                                                                                len(self.idents[self.active]),
                                                                                files[0].replace('.nc', ''))
                self.data[name].attrs['source'] = 'e.g.: ' + files[0]
            else:
                # as dictionary
                self.data[name] = {idata.sonde.item(0): idata for idata in data}
            self.files[name] = files

    def save_network(self, directory=None):
        """ Save Network Class with pickle to file

        Args:
            directory (str): output directory

        Notes:
              Creates a network_[timestamp].pkl file in directory or
              config.outdir + [NETWORK NAME]

        """
        import pickle
        import os
        from ..fun import now
        from .. import config

        if directory is None:
            directory = config.outdir + '/' + self.name

        os.makedirs(directory, exist_ok=True)
        directory = directory + '/network_' + now(timespec='minutes') + '.pkl'
        pickle.dump(self, open(directory, 'wb'), protocol=-1)
        print('Saved to:', directory)

    def save_data(self, name=None, directory=None):
        import os
        import pickle
        from .. import config

        if directory is None:
            directory = config.outdir + '/' + self.name

        os.makedirs(directory, exist_ok=True)
        if name is None:
            name = list(self.data)

        for iname in name:
            if iname in self.data:
                filename = directory + '/' + iname
                if isinstance(self.data[iname], (xr.DataArray, xr.Dataset)):
                    self.data[iname].to_netcdf(filename + '.nc')
                    print("Saved [Xarray] netcdf")
                else:
                    # pickle
                    pickle.dump(self.data[name], open(filename + '.pkl', 'wb'), protocol=-1)
                    print("Saved [",type(self.data[iname]),"] pickle")
        print('Done')

    # what does it do?
    # run function on files or run function on data set?
    def apply_function(self, variables=None, npp=10, **kwargs):
        """ Apply a function

        Args:
            variables:
            npp:
            **kwargs:

        Returns:
            list

        Notes:
            Example what does it do?

        Examples:
            ???

        """
        from ..fun.mp import make_process_list, execute_process
        from ..net.process import process_function
        from ..fun import levelup
        # kwargs  : exec
        stuff = make_process_list(self.idents[self.active], process_function, args=(variables,), **levelup(**kwargs))
        errors, data = execute_process(stuff, npp=npp, return_value=True,
                                       loop='loop' in kwargs.keys())
        # add combine to Dataset ?
        return data

    def plot_map(self, filename=None, **kwargs):
        """ Create a Network Map with active stations

        Args:
            filename (str): output filename
            **kwargs:

        Returns:
            matplotlib.axes : plot axes
        """
        from ..plot import map
        kwargs.update({'title': self.name + kwargs.get('title','')})
        ax = map.points(self.lon[self.active], self.lat[self.active], **kwargs)

        if filename is not None:
            plt.savefig(filename, **kwargs)

        return ax

    def plot_neighboring_stations(self, ident, distance, filename=None, **kwargs):
        """ Plot neighboring stations within a certain distance

        Args:
            ident (str): radiosonde id
            distance (int, float): distance in km
            filename (str): output filename  
            **kwargs:

        Returns:

        """
        import matplotlib.patches as mpatches
        from ..plot import map

        neighbors = self.neighboring_stations(ident, distance * 2)
        neighbors.name = 'distance'

        idx = np.in1d(self.idents, neighbors.index.values)

        stations = pd.DataFrame({'lon': self.lon[idx], 'lat': self.lat[idx]}, index=np.array(self.idents)[idx])
        stations['distance'] = neighbors
        stations = stations.sort_values('distance')
        print(stations)
        stations['color'] = 'b'
        stations.loc[ident, 'color'] = 'r'
        stations.loc[stations.distance > distance, 'color'] = 'w'

        # todo add title with distance
        ax = map.points(stations.lon, stations.lat, labels=stations.index.values, color=stations.color, **kwargs)
        lon_radius = distance / (111.32 * np.cos(stations.lat[0] * np.pi / 180.))
        lat_radius = distance / 110.574
        print(lon_radius, lat_radius)
        ax.add_patch(
            mpatches.Ellipse(xy=[stations.lon[0], stations.lat[0]], width=lon_radius * 2., height=lat_radius * 2.,
                             color='red', alpha=0.2,
                             zorder=2, transform=ax.projection))

        if filename is not None:
            plt.savefig(filename, **kwargs)

        return ax

    def activate(self, logic=None, ident=None):
        """ activate a sonde

        Args:
            logic (list, array): boolean logic for activating/ deactivating sondes
            ident (list): active stations / others will be deactivated

        """
        if logic is not None:
            print(sum(self.active), sum(logic))
            self.active = np.asarray(logic)

        elif ident is not None:
            if not isinstance(ident, list):
                ident = [ident]

            print(sum(self.active), len(ident))
            for i in ident:
                if i in self.idents:
                    self.active[self.idents == i] = False if self.active[self.idents == i] else True
            print(sum(self.active), len(ident))
        else:
            pass

    def availability(self, filename=None, reference_year=None):
        """ Output availability based on stations as ascii graphic

        Args:
            filename (str): output filename
            reference_year (int): Reference Year to be placed inside Ascii graphic

        Returns:
            str
        """
        if hasattr(self, 'stations'):
            if sum(self.stations.columns.isin(['start', 'end'])) == 2:
                imin = self.stations[self.active]['start'].min()
                imax = self.stations[self.active]['end'].max()
                summary = "Availability of %d Stations between %d - %d \n" % (len(self.idents[self.active]), imin, imax)
                if reference_year is not None and imin < reference_year < imax:
                    line = np.array([' '] * (imax - imin + 1))
                    iref = reference_year - imin - 1
                else:
                    line = np.array([' '] * (imax - imin))
                    iref = -1

                for i in sorted(self.idents[self.active]):
                    istart = self.stations.loc[i, 'start']
                    iend = self.stations.loc[i, 'end']
                    line[istart - imin: iend - imin] = 'X'
                    if iref > 0:
                        line[iref] = '|'
                    try:
                        ii = "%06d" % int(i)
                    except:
                        ii = i
                    summary += "%12s [%s]\n" % (ii, "".join(line))
                    line[:] = ' '

                if filename is not None:
                    with open(filename, 'w') as f:
                        f.write(summary)
                    print(filename)
                else:
                    return summary

# def _add_id(ds, ivars, invert=False):
#     if 'ident' in ds.attrs:
#         ds.coords['sonde'] = ds.attrs['ident']
#     elif 'station_id' in ds.attrs:
#         ds.coords['sonde'] = ds.attrs['station_id']
#     else:
#         ds.coords['sonde'] = ds.encoding['source'].split('/')[-2]
#
#     if ivars is not None:
#         if invert:
#             ds = ds[[i for i in list(ds.data_vars) if i not in ivars]]
#         else:
#             ds = ds[[i for i in list(ds.data_vars) if i in ivars]]
#     return ds
