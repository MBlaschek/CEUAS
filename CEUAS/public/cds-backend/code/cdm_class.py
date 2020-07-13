#!/usr/bin/env python3
__version__ = '0.1'
__author__ = 'MB'
__status__ = 'dev'
__date__ = 'Mon Jul 13 10:29:21 CEST 2020'
__institute__ = 'UNIVIE'
__github__ = 'git@github.com:MBlaschek/CEUAS.git'
__doc__ = """
CDM Class to read ragged radiosonde profiles in Python
CDM_class v%s
Maintained by %s at %s
Github: %s [%s]
Updated: %s
""" % (__version__, __author__, __institute__, __github__, __status__, __date__)

"""

Create a class object based on Xarray that can read the CDM and use the information
or h5Py


tmp = cdm_dataset('..merged.nc')
tmp
 - groups
    - variables + metadata + shape
    - read information for datetime and z_coordinate
    - observed variables

read_group(variables, datetime_var) - convert to datetimes

write_group(name, data, attributes)

read_to_cube(group, variables)
 - group : observation_table
 - groupby:
    - observed_variable -> 85, ...
    - date_time
    - z_coordinate


"""
from datetime import datetime

import h5py
import numpy as np
import xarray as xr

odb_codes = {'t': 85, 'rh': 38, 'td': 36, 'dpd': 34, 'z': 117, 'dd': 106, 'ff': 107, 'u': 104, 'v': 105,
             'q': 39}
cdm_codes = {'z': 1, 't': 2, 'u': 3, 'v': 4, 'dd': 111, 'ff': 112, 'td': 59, 'dpd': 299, 'rh': 29, 'p': 999}

std_plevs = np.asarray([10, 20, 30, 50, 70, 100, 150, 200, 250, 300, 400, 500, 700, 850, 925, 1000])



@xr.register_dataset_accessor('cdm')
class CDMAccessor(object):
    def __init__(self, index, **kwargs):
        self.index = index
        self._names = []
        for ikey, ival in kwargs.items():
            setattr(self, ikey, ival)
            self._names.append(ikey)

    def keys(self):
        return self._names


class CDMVariable:
    def __init__(self, origin, name: str, **kwargs):
        self._name = name
        self._origin = origin
        self._names = []
        for ikey, ival in kwargs.items():
            setattr(self, ikey, ival)
            self._names.append(ikey)

    def __setitem__(self, key, value):
        self._names.append(key)
        self.__setattr__(key, value)

    def __getitem__(self, item):
        return self._origin[item]

    def __repr__(self):
        return " ".join(["{}".format(str(getattr(self, i))) for i in self._names])

    def keys(self):
        return self._names

    def attrs(self, name=None):
        if name is not None:
            if name in self._origin.attrs.keys():
                return self._origin.attrs[name]
            else:
                return None
        else:
            return self._origin.attrs.keys()


class CDMGroup(CDMVariable):
    def __repr__(self):
        return self._name + ":\n" + "\n".join(["{:_<50} : {}".format(i, getattr(self, i)) for i in self.keys()])


class CDMDataset:
    # __slots__ = ['filename', 'file', 'groups', 'data']

    def __init__(self, filename: str):
        self.filename = filename
        self.file = h5py.File(filename, 'r')
        print("[OPEN] %s" % self.filename)
        self.groups = []
        # Check everything
        self.inquire()

    def __getitem__(self, item):
        return self.__getattribute__(item)

    def __repr__(self):
        text = "Filename: " + self.filename
        text += "\nGroups: \n"
        text += "\n".join(["- {} - {} : {}".format('G' if isinstance(self[igroup], CDMGroup) else 'V', igroup,
                                                   getattr(self[igroup], 'shape')) for igroup in self.groups])
        return text

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.file.close()

    def inquire(self):
        try:
            for igroup in self.file.keys():
                self.groups.append(igroup)
                if isinstance(self.file[igroup], h5py.Group):
                    jgroup = CDMGroup(self.file[igroup], igroup)

                    for ivar in self.file[igroup].keys():
                        shape = self.file[igroup][ivar].shape
                        jgroup[ivar] = CDMVariable(self.file[igroup][ivar], ivar, shape=shape)

                    jgroup['shape'] = len(jgroup.keys())
                    setattr(self, igroup, jgroup)  # attach to class
                else:
                    setattr(self, igroup, CDMVariable(self.file[igroup], igroup, shape=self.file[igroup].shape))

        except Exception as e:
            print(repr(e))
            self.close()

    def close(self):
        self.file.close()
        print("[CLOSED] %s" % self.filename)

    def reopen(self, mode: str = 'r', **kwargs):
        self.file.close()
        print("[REOPEN] %s [%s]" % (self.filename, mode))
        self.file = h5py.File(self.filename, mode=mode, **kwargs)

    def get_observed_variable(self, variable: str or int, name: str = None, dates: list = None, plev: list = None,
                              observed_variable_name: str = 'observed_variable',
                              date_time_name: str = 'date_time',
                              date_time_in_seconds: bool = False,
                              z_coordinate_name: str = 'z_coordinate',
                              group: str = 'observations_table',
                              return_coordinates: bool = False
                              ):
        """ Get the indices of a variable in the CDM file

        Args:
            name:
            variable:
            dates:
            plev:
            observed_variable_name:
            date_time_name:
            date_time_in_seconds:
            z_coordinate_name:
            group:

        Returns:

        """

        if group not in self.groups:
            raise ValueError('Missing observations_table group?')

        if isinstance(variable, str):
            if variable in odb_codes.keys():
                variable = odb_codes[variable]

        if variable not in odb_codes.values():
            raise ValueError('Variable Code not in CDM Codes', variable, str(odb_codes))

        if observed_variable_name not in self[group].keys():
            raise ValueError('Observed variable not found:', observed_variable_name, self[group].keys())

        return_index = True
        if name is not None:
            if name not in self[group].keys():
                raise ValueError('Variable not found:', name, self[group].keys())
            return_index = False

        logic = (self[group][observed_variable_name][()] == variable)
        xdates = None
        if dates is not None:
            if date_time_name not in self[group].keys():
                raise ValueError('Datetime variable not found:', date_time_name, self[group].keys())

            if not isinstance(dates, list):
                dates = [dates]

            units = self[group][date_time_name].attrs('units')
            if units is None:
                units = 'seconds since 1900-01-01 00:00:00'
            else:
                units = units.decode()

            # todo add option if units are different
            # recordindex, recordtimestamp -> start end of profiles
            # profiles = np.split(times, recidx)  -> list of profiles []
            if not date_time_in_seconds:
                # to seconds since
                refdate = datetime(year=1900, month=1, day=1)
                if len(dates) == 1:
                    dates = [int((dates[0] - refdate).total_seconds())]
                else:
                    dates = [int((dates[0] - refdate).total_seconds()), int((dates[1] - refdate).total_seconds())]

            xdates = self[group][date_time_name][()]
            if len(dates) == 1:
                # searches in date_time
                # but could search in recordtimestamp (shorter) -
                # everything before false, everything after false
                logic = logic & (xdates == dates[0])
            else:
                logic = logic & ((xdates >= dates[0]) & (xdates <= dates[1]))

        xplevs = None
        if plev is not None:
            if z_coordinate_name not in self[group].keys():
                raise ValueError('Pressure variable not found:', z_coordinate_name, self[group].keys())
            if not isinstance(plev, list):
                plev = [plev]
            xplevs = self[group][z_coordinate_name][()]
            if len(plev) == 1:
                logic = logic & (xplevs == plev[0])
            else:
                logic = logic & (np.in1d(xplevs, plev))

        # logic = np.where(logic)[0]  # makes indices
        if return_coordinates:
            if xdates is None:
                xdates = self[group][date_time_name][()]
            if xplevs is None:
                xplevs = self[group][z_coordinate_name][()]
        if return_index:
            if return_coordinates:
                return logic, xdates[logic], xplevs[logic]
            return logic
        if return_coordinates:
            return self[group][name][()][logic], xdates[logic], xplevs[logic]
        return self[group][name][()][logic]

    def read_group(self, group: str, variables: str or list, decode_byte_array: bool = True, **kwargs):
        """ Return data from variables and concat string-byte arrays


        Args:
            group:
            variables:
            decode_byte_array:
            **kwargs:

        Returns:

        """
        if not isinstance(variables, list):
            variables = [variables]

        data = {}
        if group in self.groups:
            jgroup = getattr(self, group)
            for ivar in variables:
                if ivar in jgroup.keys():
                    if len(jgroup[ivar].shape) > 1 and decode_byte_array:
                        data[ivar] = jgroup[ivar][()].astype(object).sum(axis=1).astype(str)  # concat byte-char-array
                    else:
                        data[ivar] = jgroup[ivar][()]  # get the full numpy array
        return data

    def write_group(self, group: str, data: dict, index, **kwargs):
        # need to have an index for mapping
        if group in self.groups:
            # WOWOWOWOW
            pass
        else:
            if self.file.mode == 'r':
                self.reopen(mode='a')
                # igroup = self.file.create_group()

    def read_data_to_cube(self, variables: dict,
                          dates: list = None,
                          plevs: list = None,
                          plev_in_hPa: bool = True,
                          date_time_name: str = 'date_time',
                          **kwargs):
        """

        Args:
            variables: {variable : observed_variable}
            dates:
            plevs:
            plev_in_hPa:
            **kwargs:

        Returns:


        Notes:
            Converting from ragged arraay to Cube requires a index
            save the index for the group
            for later writing
            return a Xarray Dataset with cdm Accessor
        """
        for ivar in variables.keys():
            igroup = 'observations_table'
            if '/' in ivar:
                igroup, ivar = ivar.split('/')

            if ivar not in self[igroup].keys():
                raise ValueError("Variable not in Group", ivar, igroup)

        if plevs is not None:
            plevs = np.array(plevs)
            if plev_in_hPa:
                plevs *= 100  # to Pa
            if any((plevs > 110000) | (plevs < 500)):
                raise ValueError('Pressure levels outside range [5, 1100] hPa')
        else:
            plevs = std_plevs * 100  # in Pa

        pernum = {}
        for ikey, ival in variables.items():
            if ival not in pernum.keys():
                pernum[ival] = [ikey]
            else:
                pernum[ival].append(ikey)

        std_plevs_indices = np.zeros(1001, dtype=np.int32)
        for i, j in enumerate(std_plevs):
            std_plevs_indices[j] = i

        igroup = 'observations_table'
        time_unit = self[igroup][date_time_name].attrs('units')

        # Read the data and Coordinate Information & convert to Cube
        dataset = {}
        for ivar in pernum.keys():
            if len(pernum[ivar]) > 1:
                index, xdates, xplevs = self.get_observed_variable(ivar, return_coordinates=True, **kwargs)
                newdates = np.datetime64(" ".join(time_unit.split(' ')[-2:])) + xdates * np.timedelta64(1,
                                                                                                        time_unit.split(
                                                                                                            ' ')[
                                                                                                            0][0])
                for jvar in pernum[ivar]:
                    if '/' in jvar:
                        igroup, jvar = jvar.split('/')
                    data = self[igroup][jvar][()][index]

            else:
                igroup = 'observations_table'
                jvar = pernum[ivar][0]
                index, xdates, xplevs = self.get_observed_variable(ivar, return_coordinates=True, **kwargs)
                if '/' in pernum[ivar][0]:
                    igroup, jvar = jvar.split('/')
                data = self[igroup][jvar][()][index]
                itime, iobs = table_to_cube(xdates, std_plevs_indices[xplevs.astype(np.int32) // 100], data)
                newdates = np.datetime64(" ".join(time_unit.split(' ')[-2:])) + xdates * np.timedelta64(1,
                                                                                                        time_unit.split(
                                                                                                            ' ')[
                                                                                                            0][0])
                iname = "{}_{}".format(odb_codes)
                dataset[iname] = xr.DataArray(iobs, coords=(newdates, plevs), dims=('time', 'plev'), name=jvar)


def table_to_cube(time: np.array, plev: np.array, obs: np.array, nplev: int = 16):
    xtime, jtime, itime = np.unique(time, return_index=True, return_inverse=True)
    data = np.full((xtime.size, nplev), np.nan, dtype=obs.dtype)
    data[itime, plev] = obs
    return jtime, data


"""
Speed test of h5py fancy indexing:

tmp = CDMDataset('01001')

%%time  
    ...: idx = np.where(tmp.observations_table.observed_variable[()]==85) 
    ...: tmp.observations_table.observation_value[()][idx], tmp.observations_table.date_time[()][idx], tmp.observations_table.z_coordinate[()][idx] 
    ...:  
    ...:                                                                                                                                                                         
CPU times: user 18.8 s, sys: 5.17 s, total: 24 s
Wall time: 24 s


# insane indexing time

%%time  
    ...: idx = np.where(tmp.observations_table.observed_variable[()]==85) 
    ...: tmp.observations_table.observation_value[idx], tmp.observations_table.date_time[idx], tmp.observations_table.z_coordinate[idx] 
    ...:  
    ...:   

"""
