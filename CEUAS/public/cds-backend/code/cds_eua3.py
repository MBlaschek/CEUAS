#!/usr/bin/env python3
__version__ = '0.3'
__author__ = 'MB'
__status__ = 'dev'
__date__ = 'Thu Jul 16 09:25:53 UTC 2020'
__institute__ = 'UNIVIE'
__github__ = 'git@github.com:MBlaschek/CEUAS.git'
__doc__ = """
CDS_EUA Functions v%s
Maintained by %s at %s
Github: %s [%s]
Updated: %s
""" % (__version__, __author__, __institute__, __github__, __status__, __date__)

"""

Context
- This class CDMDataset with process_flat will replace the hug cds_eua2 function sets
- This class CDMDataset will be used in adjust and quality control 

Performance
- HDF5 Netcdf files should be chunked or sorted by variable and a variable lie=ke recordindex , e.g. varindex could
  give a slice per Variable to speed up reading performance.
- HDF5 files could be split by variable, would improve reading performance dramatically, especially in MP
"""

import logging
import time
from datetime import datetime, timedelta

import h5py
import numpy as np
import pandas as pd
import xarray as xr
from numba import njit

odb_codes = {'t': 85, 'rh': 38, 'td': 36, 'dpd': 34, 'z': 117, 'dd': 106, 'ff': 107, 'u': 104, 'v': 105,
             'q': 39}
cdm_codes = {'z': 1, 't': 2, 'u': 3, 'v': 4, 'dd': 111, 'ff': 112, 'td': 59, 'dpd': 299, 'rh': 29, 'p': 999}
std_plevs = np.asarray([10, 20, 30, 50, 70, 100, 150, 200, 250, 300, 400, 500, 700, 850, 925, 1000])
logger = logging.getLogger('upperair.cdm')

if not logger.hasHandlers():
    import sys

    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.DEBUG)  # respond only to Warnings and above
    # create formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s | %(funcName)s - %(levelname)s - %(message)s')
    # add formatter to ch
    ch.setFormatter(formatter)
    # add ch to logger
    logger.addHandler(ch)


@njit(cache=True)
def calc_trajindexfast(z, zidx, idx, trajectory_index):
    # zidx=numpy.zeros(z.shape[0],dtype=numpy.int32)
    z0 = zidx[0]
    j = 0
    l = 0
    i = 0
    for i in range(z.shape[0] - 1):
        jold = j
        while idx[j] >= z[i] and idx[j] < z[i + 1]:
            trajectory_index[j] = l
            j += 1
            if j == idx.shape[0]:
                break
        if j > jold:
            zidx[l] = z0 + i
            l += 1
        if j == idx.shape[0]:
            break

    if j < idx.shape[0]:

        if z.shape[0] > 1:
            i += 1
        jold = j
        while idx[j] >= z[i]:
            trajectory_index[j] = l
            j += 1
            if j == idx.shape[0]:
                break
        if j > jold:
            zidx[l] = z0 + i
            l += 1
    zidx = zidx[:l]

    return zidx


@njit(cache=True)
def andisin_t(mask, x, v):
    jsave = 0
    for i in range(mask.shape[0]):
        if mask[i]:
            if x[i] == v[jsave]:
                mask[i] = True
            else:
                found = False
                for j in range(jsave, v.shape[0]):
                    if x[i] == v[j]:
                        found = True
                        jsave = j
                        break
                mask[i] = found
    return


@njit(cache=True)
def andisin(mask, x, v):
    for i in range(mask.shape[0]):
        if mask[i]:
            found = False
            for j in range(v.shape[0]):
                if x[i] == v[j]:
                    found = True
                    break
            mask[i] = found


@njit(cache=True)
def tohourday(hhilf, hilf, ohilf, dshift):
    ohilfold = -1
    for i in range(ohilf.shape[0]):
        if ohilfold == ohilf[i]:
            hhilf[i] = hhilf[i - 1]
            hilf[i] = hilf[i - 1]
            dshift[i] = dshift[i - 1]
        else:
            ohilfold = ohilf[i]
            hhilf[i] = ((ohilf[i] + 3599) // 3600) % 24
            if hhilf[i] == 0:
                dshift[i] = ohilf[i] % 3600 != 0
            hilf[i] = ohilf[i] // 86400
            # if dshift[i]:
            # x=0
    return


def read_standardnames() -> dict:
    """ Read the Climate and Forcast Convention data

    Returns:
        dict : Standard Names and Units for variables
    Notes:
            https://cfconventions.org/

    """
    import urllib.request
    import xml.etree.ElementTree as ET

    url = 'http://cfconventions.org/Data/cf-standard-names/69/src/cf-standard-name-table.xml'
    response = urllib.request.urlopen(url).read()
    tree = ET.fromstring(response)
    snames = ['platform_id', 'platform_name', 'latitude', 'longitude', 'time', 'air_pressure',
              'air_temperature', 'dew_point_temperature', 'relative_humidity', 'specific_humidity',
              'eastward_wind', 'northward_wind', 'wind_speed', 'wind_direction', 'geopotential', 'trajectory_label',
              'obs_minus_bg', 'obs_minus_an', 'bias_estimate']
    cdmnames = ['header_table/primary_station_id', 'header_table/station_name', 'observations_table/latitude',
                'observations_table/longitude', 'observations_table/date_time', 'observations_table/z_coordinate']
    cdmnames += 9 * ['observations_table/observation_value']
    cdmnames += ['header_table/report_id', 'era5fb/fg_depar@body', 'era5fb/an_depar@body', 'era5fb/biascorr@body']
    cf = {}
    for c, cdm in zip(snames, cdmnames):
        cf[c] = {'cdmname': cdm, 'units': 'NA', 'shortname': c}
        if c not in 'latitude longitude time air_pressure':
            cf[c]['coordinates'] = 'lat lon time plev'  # short names
    l = 0
    for child in tree:
        # print(child.tag,child.attrib)
        try:
            c = child.attrib['id']
            if c in snames:
                i = snames.index(c)
                logger.debug(c)

                cf[c]['cdmname'] = cdmnames[i]
                if child[0].text is not None:
                    cf[c]['units'] = child[0].text
                if child[2].text is not None:
                    cf[c]['shortname'] = child[2].text
                cf[c]['standard_name'] = c
        except:
            pass
        l += 1
    cf['latitude']['shortname'] = 'lat'
    cf['longitude']['shortname'] = 'lon'
    cf['air_pressure']['shortname'] = 'plev'
    cf['time']['shortname'] = 'time'
    cf['bias_estimate']['cdsname'] = 'bias_estimate'
    cf['bias_estimate']['cdmcode'] = 0
    cf['bias_estimate']['odbcode'] = 0
    cf['obs_minus_bg']['cdsname'] = 'obs_minus_bg'
    cf['obs_minus_bg']['cdmcode'] = 0
    cf['obs_minus_bg']['odbcode'] = 0
    cf['obs_minus_an']['cdsname'] = 'obs_minus_an'
    cf['obs_minus_an']['cdmcode'] = 0
    cf['obs_minus_an']['odbcode'] = 0
    cf['air_temperature']['cdsname'] = 'temperature'
    cf['air_temperature']['cdmcode'] = 85
    cf['air_temperature']['odbcode'] = 2
    cf['eastward_wind']['cdsname'] = 'u_component_of_wind'
    cf['eastward_wind']['cdmcode'] = 104
    cf['eastward_wind']['odbcode'] = 3
    cf['northward_wind']['cdsname'] = 'v_component_of_wind'
    cf['northward_wind']['cdmcode'] = 105
    cf['northward_wind']['odbcode'] = 4
    cf['wind_speed']['cdsname'] = 'wind_speed'
    cf['wind_speed']['cdmcode'] = 107
    cf['wind_speed']['odbcode'] = 112
    cf['wind_direction']['cdsname'] = 'wind_direction'
    cf['wind_direction']['cdmcode'] = 106
    cf['wind_direction']['odbcode'] = 111
    cf['relative_humidity']['cdsname'] = 'relative_humidity'
    cf['relative_humidity']['cdmcode'] = 38
    cf['relative_humidity']['odbcode'] = 29
    cf['specific_humidity']['cdsname'] = 'specific_humidity'
    cf['specific_humidity']['cdmcode'] = 39
    cf['specific_humidity']['odbcode'] = 7
    cf['dew_point_temperature']['cdsname'] = 'dew_point_temperature'
    cf['dew_point_temperature']['cdmcode'] = 36
    cf['dew_point_temperature']['odbcode'] = 59
    cf['geopotential']['cdsname'] = 'geopotential'
    cf['geopotential']['cdmcode'] = -1
    cf['geopotential']['odbcode'] = 1

    # vdict={'111':'windDirection','112':'windSpeed','1':'geopotentialHeight',
    # '2':'airTemperature','59':'dewpointTemperature','29':'relativeHumidity'}
    return cf


def do_cfcopy(fout, fin, group, idx, cf, dim0, var_selection=[]):
    """ Copy H5PY variables and apply subsetting (idx)

    Args:
        fout: output file
        fin: input file
        group: group
        idx: selection (mask)
        cf: cdm mapping of names
        dim0: record dimension name
        var_selection: variables

    """
    # cuts vars and copies attributes of observation, feedback and header tables
    tt = time.time()
    if not var_selection:
        var_selection = fin[group].keys()

    vlist = []
    clist = []
    for i in cf.keys():
        if i not in ['platform_id', 'platform_name']:
            # , 'latitude', 'longitude', 'time', 'air_pressure']:
            clist.append(i)
            # todo this should not be defined here, but from cftable
            if i in ['air_temperature', 'dew_point_temperature', 'relative_humidity', 'specific_humidity',
                     'eastward_wind', 'northward_wind', 'wind_speed', 'wind_direction', 'geopotential']:
                for fb in ['obs_minus_bg', 'obs_minus_an', 'bias_estimate']:
                    try:
                        cf[fb]['units'] = cf[i]['units']
                        cf[fb]['standard_name'] = i
                        cf[fb]['long_name'] = group.split('fb')[0].upper() + ' reanalysis ' + fb
                    except:
                        pass

    for cfk, cfv in cf.items():
        for v in var_selection:
            if group + '/' + v == cfv['cdmname']:
                vlist.append(cfv['shortname'])
                try:
                    if fin[group][v].ndim == 1:
                        try:
                            fout.create_dataset_like(vlist[-1], fin[group][v],
                                                     shape=idx.shape,
                                                     chunks=True)
                            hilf = fin[group][v][idx[0]:idx[-1] + 1]
                            if 'time' in v:
                                # convert time units
                                us = fin[group][v].attrs['units']
                                if b'hours' in us:
                                    hilf = hilf * 3600  # hilf+=int(dh[0])
                                elif b'minutes' in us:
                                    hilf = hilf * 60  # +int(dh[0])
                                elif b'seconds' in us:
                                    hilf = hilf  # //60//60+int(dh[0])
                                elif b'days' in us:
                                    hilf *= 24 * 3600

                            fout[vlist[-1]][:] = hilf[idx - idx[0]]
                        except:
                            logger.warning('not found: %s %s', group, v)
                            pass
                    else:
                        s1 = fin[group][v].shape[1]
                        fout.create_dataset_like(vlist[-1], fin[group][v],
                                                 shape=(idx.shape[0], s1),
                                                 chunks=True)
                        sname = 'string{}'.format(s1)
                        if sname not in fout.keys():
                            fout.create_dataset(sname,
                                                data=np.zeros(s1, dtype='S1'),
                                                chunks=True)
                            fout[sname].attrs['NAME'] = np.string_(
                                'This is a netCDF dimension but not a netCDF variable.')
                            fout[sname].make_scale(sname)
                        hilf = fin[group][v][idx[0]:idx[-1] + 1, :]
                        if hilf.shape[0] == 0:
                            print('x')
                        fout[vlist[-1]][:] = hilf[idx - idx[0], :]
                except:
                    # fix for missing report_id SHOULD BE REMOVED
                    hilf = np.zeros(shape=(idx.shape[0]), dtype='S10')
                    for i in range(hilf.shape[0]):
                        hilf[i] = '{:0>10}'.format(i)
                    fout.create_dataset(vlist[-1],
                                        data=hilf.view('S1'),
                                        shape=(idx.shape[0], 10),
                                        chunks=True)

                try:
                    for a in fin[group][v].attrs.keys():
                        if a not in ['DIMENSION_LIST', 'CLASS', 'external_table']:
                            if type(fin[group][v].attrs[a]) is str:
                                fout[vlist[-1]].attrs[a] = np.string_(fin[group][v].attrs[a])
                            else:
                                fout[vlist[-1]].attrs[a] = fin[group][v].attrs[a]

                    for a in cfv.keys():
                        if a not in ['shortname', 'odbcode', 'cdmcode']:
                            fout[vlist[-1]].attrs[a] = np.string_(cfv[a])
                        if a == 'units' and cfv[a] == 'NA':
                            fout[vlist[-1]].attrs[a] = np.string_('')
                        if a == 'units' and vlist[-1] == 'time':
                            ahilf = np.bytes_(fin[group][v].attrs[a])
                            fout[vlist[-1]].attrs[a] = ahilf
                except:
                    # quick fix should be removed
                    logger.warning('%s/%s has no attributes', group, v)
                l = 0
                for d in fout[cfv['shortname']].dims:
                    if len(d) > 0:
                        if l == 0:
                            fout[vlist[-1]].dims[l].attach_scale(fout[dim0])
                        else:
                            fout[vlist[-1]].dims[l].attach_scale(fout[sname])
                    l += 1
    tt = time.time() - tt
    if tt > 0.4:
        logger.warning('slow copy: %s %f s', group, tt)


def totimes(tinput):
    if type(tinput[0]) is str:
        if '-' in tinput[0]:
            ilist = np.array(tinput[0].split('-'), dtype=np.int32)
            if ilist[0] <= ilist[1]:
                ilist = np.arange(ilist[0], ilist[1] + 1)
            else:
                ilist = np.array(list(range(ilist[0], 24)) + list(range(ilist[1] + 1)), dtype=np.int32)
            out = ilist
        else:
            out = np.array(tinput, dtype=np.int32)
    else:
        out = np.array(tinput, dtype=np.int32)

    if np.min(out) < 0 or np.max(out) > 23:
        raise ValueError
    return out


def to_seconds_since(time_elapsed: int, time_unit: str, reference=None):
    """ return seconds since a Reference date applying time unit

    Args:
        time_elapsed: seconds or other time unit
        time_unit: 'hours since 1900-01-01 00:00:00'
        reference: '1900-01-01 00:00:00'

    Returns:
        int : seconds since
    """
    if reference is None:
        reference = datetime(1900, 1, 1)

    fstart = datetime.strptime(" ".join(time_unit.split(' ')[-2:]), '%Y-%m-%d %H:%M:%S')
    offset = fstart - reference
    offsets = offset.days * 24 * 3600 + offset.seconds
    fak = 1
    if 'hours' in time_unit:
        fak = 3600
    elif 'minutes' in time_unit:
        fak = 60

    secs = time_elapsed * fak + offsets
    return secs


def seconds_to_datetime(seconds, ref='1900-01-01'):
    return pd.to_datetime(seconds, unit='s', origin=ref).values


# too slow
# @np.vectorize
# def seconds_to_datetime(seconds, ref=datetime(1900, 1, 1)):
#     return ref + timedelta(seconds=np.int(seconds))
# @np.vectorize
# def seconds_to_datetime(seconds, ref='1900-01-01T00:00:00'):
#     return np.datetime64(ref) + np.timedelta64(seconds, 's')


def last_day_of_month(any_day):
    next_month = any_day.replace(day=28) + timedelta(days=4)  # this will never fail
    return next_month - timedelta(days=next_month.day)


# @np.vectorize
def convert_datetime(idate, return_string=False, return_seconds=False, reference=None, format=None):
    if isinstance(idate, str):
        # 19000101 or 01-01-1900 00:00:00
        if '-' in idate:
            if ':' in idate:
                for iformat in ['%Y-%m-%d %H:%M', '%Y-%m-%d %H:%M:%S']:
                    try:
                        idate = datetime.strptime(idate, iformat)
                        break
                    except ValueError:
                        pass
            else:
                idate = datetime.strptime(idate, '%Y-%m-%d')

    if not isinstance(idate, (datetime, np.datetime64)):
        try:
            d = int(idate)
            idate = datetime(year=d // 10000, month=d % 10000 // 100, day=d % 100)
        except ValueError:
            if (d % 100) > 28:
                idate = last_day_of_month(datetime(year=d // 10000, month=d % 10000 // 100, day=1))

    if return_seconds:
        reference = reference if reference is not None else datetime(1900, 1, 1)
        return int((idate - reference).total_seconds())
    if return_string:
        if format is None:
            return idate.strftime('%Y%m%d')  # used by the request library
        else:
            return idate.strftime(format)
    return idate


def table_to_cube(time: np.ndarray, plev: np.ndarray, obs: np.ndarray, nplev: int = 16) -> tuple:
    """ Convert a ragged variable (table) to a data cube

    Args:
        time: datetime values
        plev: pressure index values
        obs: observation values
        nplev: number of pressure levels

    Returns:
        list : time indices, obs cube [time x plev]
    """
    xtime, jtime, itime = np.unique(time, return_index=True, return_inverse=True)
    data = np.full((xtime.size, nplev), np.nan, dtype=obs.dtype)
    data[itime, plev] = obs
    return jtime, data


def process_flat(outputdir: str, cftable: dict, datadir: str, request_variables: dict) -> tuple:
    """ Process a station file with the requested variables

    Args:
        outputdir: output directory
        cftable: CF convention definitions table
        datadir: data directory
        request_variables: request dictionary

    Returns:
        str : filename of results
        str : message or error
    """
    import os
    # mimicks process_flat from cds_eua2
    try:
        msg = ''  # Message or error
        filename = ''  # Filename
        statid = request_variables.pop('statid', None)
        if statid is None:
            logger.error('No station ID (statid) specified. %s', filename)
            raise ValueError('No station ID (statid) specified')
        if statid[:3] == '0-2':
            suffix = ['']
        else:
            suffix = ['0-20000-0-', '0-20001-0-']

        for ss in suffix:
            filename = os.path.expandvars(datadir + '/' + ss + statid + '_CEUAS_merged_v0.nc')
            if os.path.isfile(filename):
                break
        cdmnamedict = {}
        for igroup, v in cftable.items():
            if "odbcode" in v.keys():
                cdmnamedict[v['cdsname']] = igroup

        # todo this could be changed to the cf.keys() -> cdm names of the variables
        filename_out = outputdir + '/dest_' + statid + '_' + cdmnamedict[
            request_variables['variable']] + '.nc'

        with CDMDataset(filename=filename) as data:
            data.read_write_request(filename_out=filename_out,
                                    request=request_variables,
                                    cf_dict=cftable)

    except Exception as e:
        logger.error('Exception %s occurred while reading %s', repr(e), filename)
        return '', 'exception "{}" occurred while reading {}'.format(e, filename)

    return filename_out, msg


###############################################################################
#
# Xarray return Accessor for index attachment, might be useful for reading
# and rewriting to h5py CDM backend files
# as it is currently used in adjust to write back adjustments and interpolated
# results.
#
###############################################################################


@xr.register_dataarray_accessor('cdm')
class CDMAccessor:
    def __init__(self, data_array):
        self._data_array = data_array
        self.ragged_index = None  # np.array
        self.filename = None  # str
        self.info = None  # variable, plevs, date


def xarray_to_cdm_xarray(data: xr.DataArray, index: np.ndarray):
    new = CDMAccessor(data)
    new.cdm.ragged_index = index
    return new


###############################################################################
#
# CDM Backend / Frontend Classes
# 1. CDMVariable
#       * origin -> Link to H5PY file location
#       * Attributes, e.g.: shape : shape of H5PY Variable
# 2. CDMGroup
#       just a wrapper class for different printing
# 3. CDMDataset
#       Main class with H5PY File handlers
#       * file  : H5PY file handler (open)
#       * filename : string name of the file that is open
#       * groups : list of CDMGroups and CDMVariables
#       * hasgroups : True CDM Backend file / False CDM Frontend file
#                     depending on the file different variables are present
###############################################################################


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
        return self._name + " : " + " ".join(["{}".format(str(getattr(self, i))) for i in self._names])

    def keys(self):
        return self._names


class CDMGroup(CDMVariable):
    def __repr__(self):
        text = ''
        for i in self.keys():
            ivar = getattr(self, i)
            if 'HDF5' in str(getattr(ivar, '_origin', 'HDF5')):
                istatus = ' '
            else:
                istatus = 'L'
            text += "\n{:_<50} :{:1}: {}".format(i, istatus, getattr(ivar, 'shape', ''))
        return self._name + ":\n" + text

    def __getitem__(self, item):
        return self.__getattribute__(item)


class CDMDataset:
    # __slots__ = ['filename', 'file', 'groups', 'data']

    def __init__(self, filename: str):
        self.filename = filename
        self.file = h5py.File(filename, 'r')
        logger.debug("[OPEN] %s", self.filename)
        self.hasgroups = False
        self.groups = []
        # Check everything
        self.inquire()

    def __getitem__(self, item):
        return self.__getattribute__(item)

    def __repr__(self):
        text = "Filename: " + self.filename
        text += "\n(G)roups/(V)ariables: \n"
        for igroup in self.groups:
            ivar = getattr(self, igroup)
            if 'HDF5' in str(getattr(ivar, '_origin', 'HDF5')):
                istatus = ' '
            else:
                istatus = 'L'
            text += "\n - {} | {:_<45} :{:1}: {}".format('G' if isinstance(self[igroup], CDMGroup) else 'V',
                                                         igroup, istatus, getattr(ivar, 'shape'))
        # text += "\n".join(["- {} - {} : {}".format('G' if isinstance(self[igroup], CDMGroup) else 'V', igroup,
        #                                           getattr(self[igroup], 'shape')) for igroup in self.groups])
        return text

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

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
                    self.hasgroups = True
                else:
                    setattr(self, igroup, CDMVariable(self.file[igroup], igroup, shape=self.file[igroup].shape))

        except Exception as e:
            logger.debug(repr(e))
            self.close()

    def close(self):
        """
        Close H5py file
        """
        self.file.close()
        logger.debug("[CLOSED] %s", self.filename)

    def reopen(self, mode: str = 'r', **kwargs):
        """ Reopen the H5py file with different mode

        Args:
            mode: r, r+, w, w+, a
            **kwargs: other options to h5py.File( )
        """
        self.file.close()
        logger.debug("[REOPEN] %s [%s]", (self.filename, mode))
        self.file = h5py.File(self.filename, mode=mode, **kwargs)

    def read_attributes(self, name: str, group: str = None, subset: list = None):
        """ Return all attributes or a subset for a variable

        Args:
            name: variable
            group: h5py group
            subset: list of attributes

        Returns:
            dict : Attributes of Variable
        """
        if group is not None:
            if not self.hasgroups:
                raise ValueError('CDM Frontend files do not have groups. Remove group keyword!', group)
            if group not in self.groups:
                raise ValueError('Group not found', group)
            if name not in self[group].keys():
                raise ValueError('Variable not found', name)
            name = "{}/{}".format(group, name)  # HDF5 access

        if subset is not None:
            if not isinstance(subset, list):
                subset = [subset]

        attributes = {}
        # can raise a KeyError if name not in file
        for ikey, ival in self.file[name].attrs.items():
            if subset is not None and ikey in subset:
                attributes[ikey] = ival.decode()
            else:
                if ikey not in ['DIMENSION_LIST']:
                    attributes[ikey] = ival.decode()
        return attributes

    def read_write_request(self, filename_out: str, request: dict, cf_dict: dict):
        """ This is the basic request used in the cds_eua2 script

        Args:
            filename_out: request output filename, e.g. : /data/public//tmp//049532808458/dest_0-20000-0-10393_air_temperature.nc
            request: Request dictionary, minimum: {'variable' : 'temperature'}
            cf_dict: CF Convention for Names (results from read_standardnames() Function)

        """
        # version of process_flat with Leo's optimizations
        rname = self.filename.split('/')[-1]

        if 'variable' not in request.keys():
            logger.error('No variable specified %s %s', str(request.keys()), rname)
            raise ValueError('No variable specified')

        cdmdict = {}
        cdmnamedict = {}
        for igroup, v in cf_dict.items():
            if "odbcode" in v.keys():
                cdmdict[v['cdsname']] = v['cdmcode']
                cdmnamedict[v['cdsname']] = igroup
        time0 = time.time()
        # get time information from file
        # recordindex -> index of start position [start next[
        # recordtimestamp -> datetime seconds since for index position
        # Array [starttimestamp, index]
        dateindex = np.array((self.file['recordtimestamp'][:],
                              self.file['recordindex'][:self.file['recordtimestamp'].shape[0]]))
        if 'date' in request.keys():
            #
            # format: [DATE], [START, END], [DATE, DATE, DATE, ...] (min 3)
            #
            if len(request['date']) > 2:
                dsec = []
                for ievent in request['date']:
                    dsec.append(convert_datetime(ievent, return_seconds=True))
                dsec = np.asarray(dsec, dtype=np.int)
            else:
                dsec = np.arange(convert_datetime(request['date'][0], return_seconds=True),
                                 convert_datetime(request['date'][-1], return_seconds=True) + 1,
                                 86400, dtype=np.int)
            # if time interval e.g. 21h-3h is chosen, the previous day must be extracted as well.
            prevday = 0
            if 'time' in request.keys():
                tlist = totimes(request['time'])
                if tlist[0] > tlist[-1]:
                    prevday = 1

            dsec[-1] += 86399
            # request range [from - to]  (can be really large)
            # -86400 in order to get late launches of previous day
            didx = np.where(np.logical_and(dateindex[0, :] >= dsec[0] - 86400,
                                           dateindex[0, :] <= dsec[-1]
                                           )
                            )[0]
            if didx.shape[0] == 0:
                logger.warning('No data in time interval %s', rname)
                raise ValueError('No data in specified time interval')

            didx = [didx[0], didx[-1]]
        else:
            # read everything ? todo check if recordindex[-1] is really everything. I don't think so.
            didx = [0, dateindex.shape[1] - 1]
            request_date = ['']
            prevday = 0

        if didx[-1] + 1 == dateindex.shape[1]:
            trange = [dateindex[1, didx[0]],
                      self.file['observations_table']['observation_value'].shape[0]]  # Maximum
        else:
            trange = [dateindex[1, didx[0]], dateindex[1, didx[1] + 1]]  # Well within

        logger.debug('Datetime selection: %d - %d [%5.2f s] %s', trange[0], trange[1],
                     time.time() - time0, rname)
        mask = np.ones(trange[1] - trange[0], dtype=np.bool)
        criteria = {'variable': 'observations_table/observed_variable',
                    'pressure_level': 'observations_table/z_coordinate',
                    'time': 'observations_table/date_time'}
        # todo understand why that is not necessary ?
        # this causes the performance issue
        # Maybe that is not executed in cds_eua2
        #            'date': 'observations_table/date_time'}

        for ckey, cval in criteria.items():
            if ckey in request.keys():
                try:
                    if ckey == 'time':
                        us = np.string_(self.file[cval].attrs['units'])
                        dh = us.split(b' ')[-1].split(b':')
                        if b'seconds' not in us:
                            logger.warning('Units not given in seconds, %s %s', us, rname)
                            raise RuntimeError('Units not given in seconds')

                        ohilf = self.file[cval][trange[0]:trange[1]]
                        # add almost an hour (3600-1 sec) to account for difference between
                        # releasetime and nominal time. Time unit in CDS interface is hours.
                        hhilf = np.empty_like(ohilf, dtype=np.int32)
                        dshift = np.zeros_like(ohilf, dtype=np.int32)
                        hilf = np.empty_like(ohilf, dtype=np.int32)
                        tohourday(hhilf, hilf, ohilf, dshift)
                        # tohour(hhilf,ohilf,dshift)
                        # today(hilf,ohilf)
                        tlist = totimes(request[ckey])
                        if prevday == 1:
                            ttlist = tlist[tlist < tlist[0]]  # use only late hours of the day before
                        else:
                            andisin(mask, hhilf, tlist)

                    elif ckey == 'date':
                        us = np.string_(self.file[cval].attrs['units'])
                        dh = us.split(b' ')[-1].split(b':')
                        if b'seconds' not in us:
                            logger.warning('Units not given in seconds, %s %s', us, rname)
                            raise RuntimeError('Units not given in seconds')

                        if 'ohilf' not in locals():
                            ohilf = self.file[cval][trange[0]:trange[1]]
                            hhilf = np.empty_like(ohilf)
                            dshift = np.zeros_like(ohilf, dtype=np.int32)
                            hilf = np.empty_like(ohilf)
                            tohourday(hhilf, hilf, ohilf, dshift)
                            tlist = totimes(['0-23'])
                        # tohour(hhilf,ohilf,dshift)
                        # today(hilf,ohilf)
                        else:
                            tlist = totimes(request['time'])

                        dsec = np.array(dsec) // 86400
                        if prevday == 1:
                            logger.debug('selecting previous day %s', rname)
                            # imask=numpy.zeros_like(mask)
                            ttlist = tlist[tlist < tlist[0]]  # use only late hours of the day before
                            # imask[:]=numpy.logical_and(numpy.isin(hilf,dsec),hhilf<=ttlist[-1])
                            imask = hhilf <= ttlist[-1]
                            andisin_t(imask, hilf + dshift, dsec)
                            ttlist = tlist[tlist >= tlist[0]]  # use only early hours of last day
                            imask2 = hhilf >= ttlist[0]
                            # imask[:]=numpy.logical_or(imask,
                            # numpy.logical_and(numpy.isin(hilf,dsec-1),hhilf>=ttlist[0]))
                            andisin_t(imask2, hilf, dsec - 1)
                            # print('nach andisin', time.time() - t)
                            imask = np.logical_or(imask, imask2)
                            mask = np.logical_and(mask, imask)
                        else:
                            andisin_t(mask, hilf + dshift, dsec)

                    elif ckey == 'variable':
                        andisin(mask, self.file[cval][trange[0]:trange[1]],
                                np.int32(np.unique(cdmdict[request[ckey]])))

                    else:
                        andisin(mask, self.file[cval][trange[0]:trange[1]],
                                np.int32(np.unique(request[ckey])))
                    logger.debug('Finished %s [%5.2f s] %s', ckey, time.time() - time0, rname)
                except MemoryError as e:
                    logger.error('Error %s occurred while checking criteria', repr(e))
                    raise MemoryError('"{}" occurred while checking criteria'.format(e))

        logger.debug('Finished mask %s [%5.2f s]', rname, time.time() - time0)
        idx = np.where(mask)[0] + trange[0]  # make index for file subsetting
        if len(idx) == 0:
            logger.warning('No matching data found %s', rname)
            raise ValueError('No matching data found')

        logger.debug('Data found: %d %s', len(idx), rname)
        trajectory_index = np.zeros_like(idx, dtype=np.int32)
        recordindex = dateindex[1, :]  # recordindex
        zidx = np.where(np.logical_and(recordindex >= trange[0], recordindex < trange[1]))[0]
        recordindex = recordindex[zidx]
        zidx = calc_trajindexfast(recordindex, zidx, idx, trajectory_index)

        dims = {'obs': np.zeros(idx.shape[0], dtype=np.int32),
                'trajectory': np.zeros(zidx.shape[0], dtype=np.int32)}
        globatts = {'Conventions': "CF-1.7",
                    'source': "radiosonde",
                    'featureType': "trajectory"}

        snames = ['report_id', 'platform_id', 'platform_name', 'observation_value', 'latitude',
                  'longitude', 'time', 'air_pressure', 'trajectory_label']

        logger.debug('Request-keys: %s', str(request.keys()))

        if isinstance(request['variable'], list):
            snames.append(cdmnamedict[request['variable'][0]])
        else:
            snames.append(cdmnamedict[request['variable']])

        if 'fbstats' in request.keys():
            if isinstance(request['fbstats'], list):
                for c in request['fbstats']:
                    snames.append(cdmnamedict[c])
            else:
                snames.append(cdmnamedict[request['fbstats']])
        name_to_cdm = {}
        for ss in snames:
            try:
                name_to_cdm[ss] = cf_dict[ss]
            except:
                pass

        logger.debug('Writing: %s', filename_out)
        with h5py.File(filename_out, 'w') as fout:
            i = 0
            for d, v in dims.items():
                fout.create_dataset(d, data=v)
                fout[d].attrs['NAME'] = np.string_('This is a netCDF dimension but not a netCDF variable.')
                fout[d].make_scale(d)  # resolves phony_dim problem
                i += 1
            fout.create_dataset('trajectory_index', data=trajectory_index)
            fout['trajectory_index'].attrs['long_name'] = np.string_(
                "index of trajectory this obs belongs to")
            fout['trajectory_index'].attrs['instance_dimension'] = np.string_("trajectory")
            fout['trajectory_index'].attrs['coordinates'] = np.string_("lat lon time plev")
            for igroup in self.file.keys():
                if not isinstance(self.file[igroup], h5py.Group):
                    continue

                if igroup in ['observations_table']:
                    # only obs, feedback fitting criteria (idx) is copied
                    do_cfcopy(fout, self.file, igroup, idx, name_to_cdm, 'obs',
                              var_selection=['observation_id', 'latitude', 'longitude', 'z_coordinate',
                                             'observation_value', 'date_time'])
                    # 'observed_variable','units'
                    logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)

                elif igroup in ['era5fb']:
                    # only obs, feedback fitting criteria (idx) is copied
                    if 'fbstats' in request.keys():
                        try:
                            do_cfcopy(fout, self.file, igroup, idx, name_to_cdm, 'obs',
                                      var_selection=['fg_depar@body', 'an_depar@body',
                                                     'biascorr@body'])
                            # ['vertco_reference_1@body','obsvalue@body','fg_depar@body'])
                            logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)
                        except KeyError:
                            raise KeyError('no {} found in {}'.format(name_to_cdm[request['fbstats'][0]][
                                                                          'cdmname'], self.filename))

                elif igroup in ['header_table']:
                    # only records fitting criteria (zidx) are copied
                    do_cfcopy(fout, self.file, igroup, zidx, name_to_cdm, 'trajectory',
                              var_selection=['report_id'])
                    logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)
                    # ,'station_name','primary_station_id'])
                    # todo could be read from the observations_table

                elif 'station_configuration' in igroup:
                    # only records fitting criteria (zidx) are copied
                    try:
                        sh = self.file[igroup]['primary_id'].shape[1]
                        fout.attrs['primary_id'] = self.file[igroup]['primary_id'][0].view('S{}'.format(sh))[0]
                        sh = self.file[igroup]['station_name'].shape[1]
                        fout.attrs['station_name'] = self.file[igroup]['station_name'][0].view('S{}'.format(sh))[0]
                    except:
                        logger.warning('No primary_id in %s', filename_out)
                    logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)

                else:
                    # groups that are simply copied
                    pass
            fout['trajectory_label'].attrs['cf_role'] = np.string_('trajectory_id')
            fout['trajectory_label'].attrs['long_name'] = np.string_('Label of trajectory')
            for a, v in globatts.items():
                fout.attrs[a] = np.string_(v)

            fout.attrs['history'] = np.string_(
                'Created by Copernicus Early Upper Air Service Version 0, ' + datetime.now().strftime(
                    "%d-%b-%Y %H:%M:%S"))
            fout.attrs['license'] = np.string_('https://apps.ecmwf.int/datasets/licences/copernicus/')
        logger.debug('Finished %s [%5.2f s]', rname, time.time() - time0)

    def read_observed_variable(self, varnum: int,
                               variable: str = 'observation_value',
                               dates: list = None,
                               plevs: list = None,
                               observed_variable_name: str = 'observed_variable',
                               date_time_name: str = 'date_time',
                               date_time_in_seconds: bool = False,
                               z_coordinate_name: str = 'z_coordinate',
                               group: str = 'observations_table',
                               return_coordinates: bool = False,
                               return_index: bool = False,
                               use_odb_codes: bool = True,
                               return_xarray: bool = False
                               ):
        """ Read a variable from a CDM backend file
        Uses recordtimestamp and observed_variable for subsetting and z_coordinate as well

        Args:
            varnum: 85 or 2 for temperature
            variable: CDM variable name: observation_value
            dates: [start end]
            plevs: [plevs] in Pa
            observed_variable_name:
            date_time_name:
            date_time_in_seconds:
            z_coordinate_name:
            group:
            return_coordinates: dates and pressure levels
            return_index: subset and logic array
            use_odb_codes: ODB Codes or CDM Codes
            return_xarray: convert to xarray object

        Returns:
            values
            values, dates, pressure
            trange, index
            trange, index, dates, pressure
            DataArray

        Examples:
            Read all Temperatures at 500 hPa from 2000 to 2019
            >>> tmp.read_observed_variable(85, dates=['2000-01-01','2019-12-31'], plevs=[50000])
            2020-07-15 12:58:02,518 - cdm | read_observed_variable - INFO - [READ] recordtimestamp: slice(6350527, 137146806, None)
            2020-07-15 12:59:13,630 - cdm | read_observed_variable - INFO - [READ] Observed variable 85
            2020-07-15 12:59:14,748 - cdm | read_observed_variable - INFO - [READ] pressure levels [50000]
            array([246.1 , 238.1 , 239.7 , ..., 239.92, 238.77, 237.7 ], dtype=float32)

        """
        if not self.hasgroups:
            raise RuntimeError('This function only works with CDM Backend files')

        if group not in self.groups:
            raise ValueError('Missing observations_table group?')

        if not isinstance(variable, str):
            raise ValueError('(variable) Requires a string name, not ', str(variable))

        if not isinstance(varnum, int):
            raise ValueError('(varnum) Requires a integer number, not', str(varnum))

        if use_odb_codes:
            if varnum not in odb_codes.values():
                raise ValueError('(varnum) Code not in ODB Codes', variable, str(odb_codes))
        else:
            if varnum not in cdm_codes.values():
                raise ValueError('(varnum) Code not in CDM Codes', variable, str(cdm_codes))

        if observed_variable_name not in self[group].keys():
            raise ValueError('Observed variable not found:', observed_variable_name, self[group].keys())

        if dates is not None:
            if not isinstance(dates, list):
                dates = [dates]

        if plevs is not None:
            if not isinstance(plevs, (list, np.ndarray)):
                plevs = [plevs]

        xdates = None
        trange = slice(None, None)
        #
        # Select based on recordtimestamp (only for CDM Group files)
        #
        if 'recordtimestamp' in self.groups:
            timestamp_units = None
            if dates is not None:
                # loading variables allows reusing them in memory
                timestamp = self.load_variable('recordtimestamp', return_data=True)
                #
                # Make sure units are in seconds
                #
                timestamp_units = self['recordtimestamp'].attrs.get('units', None)
                if timestamp_units is None:
                    timestamp_units = 'seconds since 1900-01-01 00:00:00'
                if 'seconds' not in timestamp_units:
                    timestamp = to_seconds_since(timestamp, timestamp_units)
                #
                # Index
                #
                timeindex = self.load_variable('recordindex', return_data=True)
                if not date_time_in_seconds:
                    # to seconds since
                    if len(dates) == 1:
                        dates = [convert_datetime(dates[0], return_seconds=True)]
                    else:
                        dates = [convert_datetime(dates[0], return_seconds=True),
                                 convert_datetime(dates[-1], return_seconds=True)]
                logic = (timestamp >= dates[0]) & (timestamp <= dates[-1])
                timeindex = timeindex[logic]
                trange = slice(timeindex[0], timeindex[-1])
                if timestamp_units != self[group][date_time_name].attrs.get('units', '').decode():
                    raise ValueError('Timeunits missmatch?', timestamp_units,
                                     self[group][date_time_name].attrs.get('units', '').decode())
                xdates = self[group][date_time_name][trange]
            else:
                trange = slice(0, self[group][date_time_name].shape[0])
                xdates = self.load_variable(date_time_name, group=group, return_data=True)
            logger.info('[READ] recordtimestamp: %s', str(trange))
        #
        # Observed Code
        #
        if dates is None:
            self.load_variable(observed_variable_name, group=group)

        logic = (self[group][observed_variable_name][trange] == varnum)
        logger.info('[READ] Observed variable %s', varnum)
        #
        # Pressure levels
        #
        xplevs = None
        if plevs is not None:
            if z_coordinate_name not in self[group].keys():
                raise ValueError('Pressure variable not found:', z_coordinate_name, self[group].keys())
            p_attrs = self.read_attributes(z_coordinate_name, group=group)
            p_units = p_attrs.get('units', 'Pa')
            if p_units != 'Pa':
                RuntimeWarning('Pressure variable wrong unit [Pa], but', p_units)

            if dates is None:
                xplevs = self.load_variable(z_coordinate_name, group=group, return_data=True)
            else:
                xplevs = self[group][z_coordinate_name][trange]
            if len(plevs) == 1:
                logic = logic & (xplevs == plevs[0])
            else:
                logic = logic & (np.in1d(xplevs, plevs))
            logger.info('[READ] pressure levels %s', str(plevs))

        if return_xarray:
            return_coordinates = True

        if return_coordinates:
            if xdates is None:
                if dates is None:
                    xdates = self.load_variable(date_time_name, group=group, return_data=True)
                else:
                    xdates = self[group][date_time_name][trange]

            if xplevs is None:
                if plevs is None:
                    xplevs = self.load_variable(z_coordinate_name, group=group, return_data=True)
                else:
                    xplevs = self[group][z_coordinate_name][trange]

        if return_xarray:
            # todo add attributes from CF Table
            data = xr.DataArray(self[group][variable][trange][logic], dims=('obs'))
            logger.info('[READ] xarray ... %s', data.shape)
            if xdates is not None:
                data[date_time_name] = ('obs', seconds_to_datetime(xdates[logic]))
                logger.debug('[READ] datetime conversion')
            if xplevs is not None:
                data[z_coordinate_name] = ('obs', xplevs[logic])
                data[z_coordinate_name].attrs.update({'units': 'Pa', 'standard_name': 'air_pressure'})
            if len(data.coords) > 0:
                data = data.set_index(obs=list(data.coords))
            data.cdm.ragged_index = trange.start + np.where(logic)[0]
            return data

        if return_index:
            # logic = trange.start + np.where(logic)[0]
            if return_coordinates:
                return trange, logic, xdates[logic], xplevs[logic]
            return trange, logic
        if return_coordinates:
            return self[group][variable][trange][logic], xdates[logic], xplevs[logic]
        return self[group][variable][trange][logic]

    def read_variable(self, name,
                      dates: list = None,
                      plevs: list = None,
                      date_time_name: str = 'time',
                      date_time_in_seconds: bool = False,
                      z_coordinate_name: str = 'plev',
                      return_coordinates: bool = False,
                      return_index: bool = False,
                      return_xarray: bool = False):
        """ Read a variable from a CDM frontend file

        Args:
            name: name of the variable
            dates: [start end] datetime selection
            plevs: [pressure levels]
            date_time_name: Name of the datetime variable
            date_time_in_seconds: dates are in seconds since 1900-01-01 00:00:00
            z_coordinate_name: Name of the pressure level variable
            return_coordinates: add coordinate information
            return_index: no data only indices
            return_xarray: convert to xarray object including coordinates and attributes

        Returns:
            values
            values, dates, pressure
            trange, index
            trange, index, dates, pressure
            DataArray
        """

        if self.hasgroups:
            raise RuntimeError('Only for CDS frontend files')

        if not isinstance(name, str):
            raise ValueError('(variable) Requires a string name, not ', str(name))

        if name not in self.groups:
            raise ValueError('Variable not found', name)

        if dates is not None:
            if date_time_name not in self.groups:
                raise ValueError('date_time_name Variable not found', date_time_name)

            if not isinstance(dates, list):
                dates = [dates]

        if plevs is not None:
            if z_coordinate_name not in self.groups:
                raise ValueError('z_coordinate_name Variable not found', z_coordinate_name)

            if not isinstance(plevs, (list, np.ndarray)):
                plevs = [plevs]

        trange = slice(None)
        xdates = None
        logic = np.ones(self[name].shape, dtype=np.bool)  # all True
        if dates is not None:
            timestamp = self[date_time_name][()]
            d_attrs = self.read_attributes(date_time_name)
            time_units = d_attrs.get('units', 'seconds since 1900-01-01 00:00:00')
            if 'seconds' not in time_units:
                dates = to_seconds_since(timestamp, time_units)
            if not date_time_in_seconds:
                # to seconds since
                if len(dates) == 1:
                    dates = [convert_datetime(dates[0], return_seconds=True)]
                else:
                    dates = [convert_datetime(dates[0], return_seconds=True),
                             convert_datetime(dates[-1], return_seconds=True)]
            logic = (timestamp >= dates[0]) & (timestamp <= dates[-1])
            timeindex = np.where(logic)[0]
            trange = slice(timeindex[0], timeindex[-1])
            xdates = timestamp[trange]
            logger.info('[READ] %s : %s [%s]', date_time_name, str(trange), time_units)

        xplevs = None
        if plevs is not None:
            xplevs = self[z_coordinate_name][trange]
            p_attrs = self.read_attributes(z_coordinate_name)
            p_units = p_attrs.get('units', 'Pa')

            if len(plevs) == 1:
                logic = logic & (xplevs == plevs[0])
            else:
                logic = logic & (np.in1d(xplevs, plevs))
            logger.info('[READ] %s : %s [%s]', z_coordinate_name, str(plevs), p_units)

        if return_xarray:
            return_coordinates = True

        if return_coordinates:
            if xdates is None:
                try:
                    xdates = self[date_time_name][trange]
                except:
                    pass
            if xplevs is None:
                try:
                    xplevs = self[z_coordinate_name][trange]
                except:
                    pass

        if return_xarray:
            data = xr.DataArray(self[name][trange][logic], dims=('obs'), attrs=self.read_attributes(name))
            if xdates is not None:
                data[date_time_name] = ('obs', seconds_to_datetime(xdates[logic]))
            if xplevs is not None:
                data[z_coordinate_name] = ('obs', xplevs[logic])
                data[z_coordinate_name].attrs.update(self.read_attributes(z_coordinate_name))
            if len(data.coords) > 0:
                data = data.set_index(obs=list(data.coords))
            data.cdm.ragged_index = trange.start + np.where(logic)[0]
            return data

        if return_index:
            if return_coordinates:
                return trange, logic, xdates[logic], xplevs[logic]
            return trange, logic
        if return_coordinates:
            return self[name][trange][logic], xdates[logic], xplevs[logic]
        return self[name][trange][logic]

    def load_variable(self, name, group: str = None, return_data: bool = False):
        """ Allow to load a variable from a group

        Args:
            name: name of the variable
            group: group
            return_data: read data and return?

        Returns:
            data
        """
        if group is not None:
            if group not in self.groups:
                raise ValueError('Group not found', group)
            if name not in self[group].keys():
                raise ValueError('Variable not found', group, name)
            self.read_group(group, name)
            if return_data:
                return self[group][name][()]
        else:
            data = self[name][()]  # Read data
            setattr(self, name, CDMVariable(data, name, shape=data.shape))
            if return_data:
                return self[name][()]

    def read_group(self, group: str, variables: str or list, decode_byte_array: bool = True, **kwargs):
        """ Return data from variables and concat string-byte arrays


        Args:
            group:
            variables:
            decode_byte_array:
            **kwargs:

        Returns:

        Examples:
            >>> tmp = CDMDataset('/raid60/scratch/federico/JUNE_TEST_MERGING_ALL/0-20000-0-01001_CEUAS_merged_v0.nc')
            This will load two Variables from the observations_table into Memory
            >>> tmp.read_group('observations_table', variables=['observed_variable', 'z_coordinate'])
            2020-07-15 12:48:29,702 - cdm | read_group - INFO - Loading ... observations_table observed_variable (137205730,)
            2020-07-15 12:49:18,027 - cdm | read_group - INFO - Loading ... observations_table z_coordinate (137205730,)
            2020-07-15 12:49:56,094 - cdm | read_group - DEBUG - observations_table [loaded]

            Loaded Variables are indicated by an L
            >>> tmp.observations_table
            observations_table:

            adjustment_id_____________________________________ : : (137205730,)
            advanced_assimilation_feedback____________________ : : (137205730,)
            bbox_max_latitude_________________________________ : : (137205730,)
            bbox_max_longitude________________________________ : : (137205730,)
            bbox_min_latitude_________________________________ : : (137205730,)
            bbox_min_longitude________________________________ : : (137205730,)
            date_time_________________________________________ : : (137205730,)
            index_____________________________________________ : : (137205730,)
            latitude__________________________________________ : : (137205730,)
            location_precision________________________________ : : (137205730,)
            longitude_________________________________________ : : (137205730,)
            numerical_precision_______________________________ : : (137205730,)
            observation_height_above_station_surface__________ : : (137205730,)
            observation_id____________________________________ : : (137205730, 11)
            observation_value_________________________________ : : (137205730,)
            observed_variable_________________________________ :L: (137205730,)
            original_precision________________________________ : : (137205730,)
            original_value____________________________________ : : (137205730,)
            report_id_________________________________________ : : (137205730, 11)
            sensor_id_________________________________________ : : (137205730,)
            source_id_________________________________________ : : (137205730, 9)
            string11__________________________________________ : : (11,)
            string9___________________________________________ : : (9,)
            units_____________________________________________ : : (137205730,)
            z_coordinate______________________________________ :L: (137205730,)
            z_coordinate_type_________________________________ : : (137205730,)
            shape_____________________________________________ : :
        """
        if not self.hasgroups:
            raise RuntimeError('CDM File has no Groups, only for CDM Backend Files')

        if not isinstance(variables, list):
            variables = [variables]

        if group in self.groups:
            jgroup = getattr(self, group)
            for ivar in variables:
                if ivar in jgroup.keys():
                    logger.info('Loading ... %s %s %s', group, ivar, jgroup[ivar].shape)
                    if len(jgroup[ivar].shape) > 1 and decode_byte_array:
                        # concat byte-char-array
                        data = jgroup[ivar][()].astype(object).sum(axis=1).astype(str)
                    else:
                        data = jgroup[ivar][()]  # get the full numpy array
                    setattr(self[group], ivar, CDMVariable(data, ivar, shape=data.shape))
                else:
                    logger.info('Skipping ... %s', ivar)

    def profile_to_dataframe(self, groups, variables, date,
                             date_time_name: str = 'date_time',
                             group: str = 'observations_table',
                             date_is_index: bool = False,
                             **kwargs):
        if not self.hasgroups:
            raise RuntimeError('Only for CDM Backend Files')

        if isinstance(groups, str):
            groups = [groups]

        if isinstance(variables, str):
            variables = [variables]

        if date is not None:
            if date_is_index:
                date = slice(date, date + 1)
            else:
                if not isinstance(date, (str, int, datetime, np.datetime64)):
                    raise ValueError('Unknown date format, only str, int, datetime allowed', date)

                date = convert_datetime(date, return_seconds=True)
                if 'recordtimestamp' in self.groups:
                    timestamp = self.load_variable('recordtimestamp', return_data=True)
                    logic = np.where(timestamp == date)[0][0]
                    timeindex = self.load_variable('recordindex', return_data=True)
                    date = slice(timeindex[logic], timeindex[logic + 1])
                else:
                    timestamp = self.load_variable(date_time_name, group=group, return_data=True)
                    logic = (timestamp == date)
                    timeindex = np.where(logic)[0]
                    date = slice(timeindex[0], timeindex[-1])
        else:
            date = slice(None)

        logger.info("Reading Profile on %s", str(date))
        data = {}
        for igroup in groups:
            for ivar in variables:
                if ivar in self[igroup].keys():
                    data[ivar] = self[igroup][ivar][date]

        logger.debug('Read variables: %s', str(data.keys()))
        for ivar in data.keys():
            if len(data[ivar].shape) > 1:
                data[ivar] = data[ivar].astype(object).sum(1).astype(str)
        return pd.DataFrame(data)

    def read_data_to_cube(self, variables: dict,
                          dates: list = None,
                          plevs: list = None,
                          plev_in_hpa: bool = False,
                          **kwargs):
        """

        Args:
            variables: {varnum : [group/variable]}
            dates: [start, stop], str, int or datetime
            plevs: [list] in Pa or hPa
            plev_in_hpa: convert input pressure or not?
            **kwargs:

        Optional Keywords:
            date_time_name: Name of the datetime variable
            z_coordinate_name: Name of the pressure level variable

        Returns:


        Notes:
            Converting from ragged arraay to Cube requires a index
            save the index for the group
            for later writing
            return a Xarray Dataset with cdm Accessor
        """
        data = {}
        if plevs is not None:
            if plev_in_hpa:
                plevs = np.array(plevs) * 100  # need Pa
            else:
                plevs = np.asarray(plevs)

            if any((plevs > 110000) | (plevs < 500)):
                raise ValueError('Pressure levels outside range [5, 1100] hPa')
        else:
            plevs = std_plevs * 100  # in Pa

        std_plevs_indices = np.zeros(1001, dtype=np.int32)
        for i, j in enumerate(plevs // 100):
            std_plevs_indices[j] = i

        if self.hasgroups:
            #
            # check variables
            #
            varnum = variables.keys()
            var_names = dict({j: i for i, j in odb_codes.items()})  # reverse
            for ivarnum in varnum:
                if ivarnum not in odb_codes.values():
                    raise ValueError('(varnum) Code not in ODB Codes', ivarnum, str(odb_codes))

                if not isinstance(variables[ivarnum], list):
                    variables[ivarnum] = [variables[ivarnum]]

                igroup = 'observations_table'
                for ivar in variables[ivarnum]:
                    if '/' in ivar:
                        try:
                            self.file[ivar]
                        except KeyError:
                            raise ValueError('(variable) not found', ivar)
                    else:
                        if ivar not in self[igroup].keys():
                            raise ValueError('(variable) not found', ivar)
            #
            # really read stuff
            #
            for ivarnum in varnum:
                logger.info('Reading ... %d  %s', ivarnum, str(variables[ivarnum]))
                if len(variables[ivarnum]) > 1:
                    # TIME SLICE, INDEX, SECONDS ARRAY, PRESSURE LEVELS
                    trange, indices, secarray, pressure = self.read_observed_variable(ivarnum,
                                                                                      dates=dates,
                                                                                      plevs=plevs,
                                                                                      return_coordinates=True,
                                                                                      return_index=True,
                                                                                      **kwargs)
                    logger.info('[CUBE] Variable Group %d %s %s', ivarnum, str(trange), str(secarray.shape))
                    igroup = 'observations_table'
                    # todo multi process read -- this part could be multi processing
                    for ivar in variables[ivarnum]:
                        if '/' in ivar:
                            obs = self.file[ivar][trange][indices]  # Observations
                            _, ivar = ivar.split('/')
                        else:
                            obs = self[igroup][ivar][trange][indices]  # Observations
                        #
                        # to Cube
                        #
                        # requires hPa for indices
                        itime, iobs = table_to_cube(secarray,
                                                    std_plevs_indices[pressure.astype(np.int32) // 100],
                                                    obs,
                                                    nplev=plevs.size)
                        if ivar == 'observation_value':
                            ivar = var_names[ivarnum]
                        else:
                            ivar = var_names[ivarnum] + '_' + ivar
                        logger.info('[CUBE] %s %s', ivar, iobs.shape)
                        # Convert to Xarray [time x plev]
                        data[ivar] = xr.DataArray(iobs,
                                                  coords=(seconds_to_datetime(secarray[itime]), plevs),
                                                  dims=('time', 'plev'),
                                                  name=ivar)
                        # todo index does not fit array because of missing levels
                        data[ivar].cdm.ragged_index = trange.start + np.where(indices)[0]
                        # read attributes ?
                else:
                    ivar = variables[ivarnum][0]
                    igroup = 'observations_table'  # Default
                    # DATA, SECONDS ARRAY, PRESSURE LEVELS
                    if '/' in ivar:
                        igroup, ivar = ivar.split('/')
                    trange, indices, secarray, pressure = self.read_observed_variable(ivarnum, ivar,
                                                                                      dates=dates,
                                                                                      plevs=plevs,
                                                                                      return_coordinates=True,
                                                                                      return_index=True,
                                                                                      group=igroup,
                                                                                      **kwargs)

                    obs = self[igroup][ivar][trange][indices]  # Observations
                    logger.info('[CUBE] Variable %d %s', ivarnum, secarray.shape)
                    # requires hPa for indices
                    itime, iobs = table_to_cube(secarray,
                                                std_plevs_indices[pressure.astype(np.int32) // 100],
                                                obs)
                    if ivar == 'observation_value':
                        ivar = var_names[ivarnum]
                    else:
                        ivar = var_names[ivarnum] + '_' + ivar
                    logger.info('[CUBE] %s %s', ivar, iobs.shape)
                    # Convert to Xarray [time x plev]
                    data[ivar] = xr.DataArray(iobs,
                                              coords=(seconds_to_datetime(secarray[itime]), std_plevs),
                                              dims=('time', 'plev'),
                                              name=ivar)
                    data[ivar].cdm.ragged_index = trange.start + np.where(indices)[0]
                    data[ivar]['plev'].attrs.update({'units': 'hPa', 'standard_name': 'air_pressure'})
                    # read attributes ?
                    # -> standardnames
        else:
            # todo replace with just self.read_variable()
            date_time_name = kwargs.get('date_time_name', 'time')
            if date_time_name not in self.groups:
                raise ValueError('date_time_name not found, specify', date_time_name)
            z_coordinate_name = kwargs.get('z_coordinate_name', 'plev')
            if z_coordinate_name not in self.groups:
                raise ValueError('z_coordinate_name not found, specify', z_coordinate_name)

            d_attrs = self.read_attributes(date_time_name)
            # units of pressure ?
            p_attrs = self.read_attributes(z_coordinate_name)
            p_unit = p_attrs.get("units", 'Pa')
            if p_unit == 'Pa':
                pfactor = 100
            elif p_unit == 'hPa':
                pfactor = 1
            else:
                raise AttributeError('Unknown pressure unit', p_unit, z_coordinate_name)

            secarray = self[date_time_name][()]
            pressure = self[z_coordinate_name][()] // pfactor

            for ivar, iobs in variables.items():
                if iobs not in self.groups:
                    raise ValueError('Known variable', iobs)

                # Read Attributes Variable
                v_attrs = self.read_attributes(iobs)
                itime, iobs = table_to_cube(secarray,
                                            std_plevs_indices[pressure.astype(np.int32)],
                                            self[iobs][()])
                logger.info('[CUBE] %s %s', ivar, iobs.shape)
                # Convert to Xarray [time x plev]
                data[ivar] = xr.DataArray(iobs,
                                          coords=(seconds_to_datetime(secarray[itime]), std_plevs),
                                          dims=('time', 'plev'),
                                          name=ivar,
                                          attrs=v_attrs)
                data[ivar]['time'].attrs.update(d_attrs)
                data[ivar]['plev'].attrs.update(p_attrs)
                # data[ivar].cdm.ragged_index = trange.start + np.where(indices)[0]

        # when the dict is converted to a dataset the cdm Accessor is emptied ? Why?
        return data  # xr.Dataset(data)

    def trajectory_data(self):
        # todo finish this function
        # return trajectory data with lon, lat, label
        if self.hasgroups:
            # what to read here? station_configuration ?
            pass
        else:
            # read trajectory information
            pass

    def write_group(self, group: str, name: str, data: np.ndarray, index: np.ndarray, varnum: int = None, **kwargs):
        """ Write a new group to HDF5 file
        or add to group ?

        Args:
            group:
            name:
            data:
            index:
            varnum:
            **kwargs:

        Returns:

        """
        # need to have an index for mapping
        # if index is present use it
        # search file for different
        # new datetime index ?
        if group in self.groups:
            # WOWOWOWOW

            pass
        else:
            if self.file.mode == 'r':
                self.reopen(mode='a')
                # igroup = self.file.create_group()

    def report_quality(self, filename: str = None):
        """

        Args:
            filename:

        Returns:

        Plans:
            Station : ID
            Region: EUR
            Country: AUT
            Location : lon, lat (most of the time)
                      other lon, lat combinations

            Merged Stations : IDs
            Distance between Stations :
                1.2 km ID
            Temporal Coverage: [Start, End]
            before 1940: xxx
            before 1950: xxx
            before 1979: xxx
            after 1979: xxx
            Current Status: active / inactive ?
            Launch times: 00 (xxx), 12 (xxx), 6 (xxx)
            Observations: variables (xxx)
            Temperature Breakpoinst: from (adjust)
            Humidity Breakpoints:
            Wind Breakpoints:
            Temperature Outliers: (based on climatology)
            Humidity Outliers:
            Wind Outliers:
        """
        # run a quality control procedure to check if the CDM file is within the standard
        # run a duplicated controller
        # run a values consistency check
        # Calculate climatology and outlier statistics
        if not self.hasgroups:
            raise RuntimeError('Only available for CDM Backend files')

        report = {'Station': '', 'Region': '', 'Country': '', 'Location': [], 'Merged_Stations': []}

        igroup = 'header_table'
        if igroup in self.groups:
            pass

        igroup = 'observations_table'
        if igroup in self.groups:
            # Load groups (full arrays)
            self.read_group(igroup, ['observed_variable', 'units', 'observation_value',
                                     'z_coordinate', 'date_time'])
            #
            # observed_variables
            #
            varcodes = self[igroup]['observed_variable'][()]
            unique_varcodes = np.unique(varcodes)
            for ivar in unique_varcodes:
                _, idx = self.read_observed_variable(ivar, return_index=True)
                # Check units
                units = self[igroup]['units'][idx]
                # Check observation_value
                obs = self[igroup]['observation_value'][idx]
                # only on standard pressure levels ?

        if 'era5fb' in self.groups:
            pass

        # if 'station_'
        if filename is not None:
            pass

        return report

    def check_cdm_duplicates(self, groups, variables,
                             observed_variable_name: str = 'observed_variable',
                             date_time_name: str = 'date_time',
                             z_coordinate_name: str = 'z_coordinate',
                             **kwargs):

        if not self.hasgroups:
            raise RuntimeError('Only for CDM Backend files')

        data = self.profile_to_dataframe(groups, variables, None, **kwargs)
        logger.info("Evaluating for duplicates ... (pandas)")
        data = data[
            data.duplicated([date_time_name, z_coordinate_name, observed_variable_name], keep=False)].sort_values(
            [date_time_name, z_coordinate_name])
        return data


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
