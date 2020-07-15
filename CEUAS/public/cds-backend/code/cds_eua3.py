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
import logging
import sys
import time
from datetime import datetime, timedelta

import h5py
import numpy as np
import xarray as xr
from numba import njit

odb_codes = {'t': 85, 'rh': 38, 'td': 36, 'dpd': 34, 'z': 117, 'dd': 106, 'ff': 107, 'u': 104, 'v': 105,
             'q': 39}
cdm_codes = {'z': 1, 't': 2, 'u': 3, 'v': 4, 'dd': 111, 'ff': 112, 'td': 59, 'dpd': 299, 'rh': 29, 'p': 999}

std_plevs = np.asarray([10, 20, 30, 50, 70, 100, 150, 200, 250, 300, 400, 500, 700, 850, 925, 1000])

logger = logging.getLogger('upperair.cdm')
# logger.setLevel(logging.DEBUG)
# if logger.hasHandlers():
#     logger.handlers.clear()
#
# if not logger.hasHandlers():
#     ch = logging.StreamHandler(sys.stdout)
#     ch.setLevel(logging.DEBUG)  # respond only to Warnings and above
#     # create formatter
#     formatter = logging.Formatter('%(asctime)s - %(name)s | %(funcName)s - %(levelname)s - %(message)s')
#     # add formatter to ch
#     ch.setFormatter(formatter)
#     # add ch to logger
#     logger.addHandler(ch)


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


@np.vectorize
def seconds_to_datetime(seconds, ref=datetime(1900, 1, 1)):
    return ref + timedelta(seconds=np.int(seconds))


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


def table_to_cube(time: list, plev: list, obs: list, nplev: int = 16) -> list:
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


#
#
@xr.register_dataset_accessor('cdm')
class CDMAccessor(object):
    def __init__(self, **kwargs):
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
        return self._name + " : " + " ".join(["{}".format(str(getattr(self, i))) for i in self._names])

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
        text = ''
        for i in self.keys():
            ivar = getattr(self, i)
            if 'HDF5' in str(getattr(ivar, '_origin', 'HDF5')):
                istatus = ' '
            else:
                istatus = 'L'
            text += "\n{:_<50} :{:1}: {}".format(i, istatus, getattr(ivar, 'shape', ''))
        return self._name + ":\n" + text


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
        text += "\nGroups: \n"
        text += "\n".join(["- {} - {} : {}".format('G' if isinstance(self[igroup], CDMGroup) else 'V', igroup,
                                                   getattr(self[igroup], 'shape')) for igroup in self.groups])
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
        self.file.close()
        logger.debug("[CLOSED] %s", self.filename)

    def reopen(self, mode: str = 'r', **kwargs):
        self.file.close()
        logger.debug("[REOPEN] %s [%s]", (self.filename, mode))
        self.file = h5py.File(self.filename, mode=mode, **kwargs)

    def read_write_request(self, filename_out: str, request: dict, cdm_dict: dict, debug:bool=False):
        """

        Args:
            filename_out:
            request:
            cdm_dict:
            debug:

        Returns:

        Notes:
            This function runs about 10 sec slower than the exact copy in cds_eua2
            Not sure why that is ... need to find the performance tweak.
        """
        # version of process_flat with Leo's optimizations
        # in the new process_flat routine, an instance will be created and this will be executed
        rname = self.filename.split('/')[-1]
        cdmdict = {}
        cdmnamedict = {}
        for igroup, v in cdm_dict.items():
            if "odbcode" in v.keys():
                cdmdict[v['cdsname']] = v['cdmcode']
                cdmnamedict[v['cdsname']] = igroup
        if 'variable' not in request.keys():
            logger.error('No variable specified %s %s', str(request.keys()), rname)
            raise ValueError('No variable specified')

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
            # read everything
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
                name_to_cdm[ss] = cdm_dict[ss]
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
        return filename_out

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
        """ Get the indices of a variable in the CDM file

        Args:
            varnum: 85 or 2 for temperature
            variable: CDM variable name: observation_value
            dates:
            plevs:
            observed_variable_name:
            date_time_name:
            date_time_in_seconds:
            z_coordinate_name:
            group:
            use_odb_codes:
            return_xarray: convert to xarray object

        Returns:
            values
            values, dates, pressure
            trange, index
            trange, index, dates, pressure

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
                timestamp = self.recordtimestamp[()]
                #
                # Make sure units are in seconds
                #
                timestamp_units = self.recordtimestamp.attrs('units')
                if timestamp_units is None:
                    timestamp_units = 'seconds since 1900-01-01 00:00:00'
                if 'seconds' not in timestamp_units:
                    timestamp = to_seconds_since(timestamp, timestamp_units)
                #
                # Index
                #
                timeindex = self['recordindex'][()]  # index
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
            else:
                trange = slice(0, self['recordindex'][-1])
            logger.info('[READ] recordtimestamp: %s', str(trange))
            # Check time units !?
            if timestamp_units is not None and timestamp_units != self[group][date_time_name].attrs.get('units',
                                                                                                        '').decode():
                raise ValueError('Timeunits missmatch?', timestamp_units,
                                 self[group][date_time_name].attrs.get('units', '').decode())
            xdates = self[group][date_time_name][trange]
        #
        # Observed Code
        #
        logic = (self[group][observed_variable_name][trange] == varnum)
        logger.info('[READ] Observed variable %s', varnum)
        #
        # Pressure levels
        #
        xplevs = None
        if plevs is not None:
            if z_coordinate_name not in self[group].keys():
                raise ValueError('Pressure variable not found:', z_coordinate_name, self[group].keys())
            # todo add units here if really in Pa
            xplevs = self[group][z_coordinate_name][trange]
            if len(plevs) == 1:
                logic = logic & (xplevs == plevs[0])
            else:
                logic = logic & (np.in1d(xplevs, plevs))
            logger.info('[READ] pressure levels %s', str(plevs))

        if return_coordinates:
            if xdates is None:
                xdates = self[group][date_time_name][trange]

            if xplevs is None:
                xplevs = self[group][z_coordinate_name][trange]

        if return_index:
            # logic = trange.start + np.where(logic)[0]
            if return_coordinates:
                return trange, logic, xdates[logic], xplevs[logic]
            return trange, logic
        if return_coordinates:
            return self[group][variable][trange][logic], xdates[logic], xplevs[logic]
        return self[group][variable][trange][logic]

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

            logger.debug('%s [loaded]', group)

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

        Returns:


        Notes:
            Converting from ragged arraay to Cube requires a index
            save the index for the group
            for later writing
            return a Xarray Dataset with cdm Accessor
        """
        data = {}
        if self.hasgroups:
            if plevs is not None:
                if plev_in_hpa:
                    plevs = np.array(plevs) * 100  # need Pa

                if any((plevs > 110000) | (plevs < 500)):
                    raise ValueError('Pressure levels outside range [5, 1100] hPa')
            else:
                plevs = std_plevs * 100  # in Pa

            std_plevs_indices = np.zeros(1001, dtype=np.int32)
            for i, j in enumerate(std_plevs):
                std_plevs_indices[j] = i

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
                        data[ivar].cdm['trange'] = trange
                        data[ivar].cdm['index'] = indices

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
                    data[ivar].cdm['trange'] = trange
                    data[ivar].cdm['index'] = indices
                    # add to CDM accessor : indices
                    # data[ivar].cdm()
        else:
            # todo add function for CDM frontend file
            # time, obs, plev -> cube
            pass
        # todo read all attributes and attach
        return xr.Dataset(data)

    # def trajectory_data(self):
    #     # return trajectory data with lon, lat, label
    #     pass

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


def process_flat(wroot: str, randdir: str, cdmtable: dict, request_variables: dict) -> tuple:
    """ Process a station file with the requested variables

    Args:
        wroot: station file root
        randdir: random request directory for output
        cdmtable: CDM definitions table
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
            filename = os.path.expandvars('$RSCRATCH/era5/odbs/merged/' + ss + statid + '_CEUAS_merged_v0.nc')
            if os.path.isfile(filename):
                break
        cdmnamedict = {}
        for igroup, v in cdmtable.items():
            if "odbcode" in v.keys():
                cdmnamedict[v['cdsname']] = igroup

        filename_out = wroot + '/' + randdir + '/dest_' + statid + '_' + cdmnamedict[
            request_variables['variable']] + '.nc'

        with CDMDataset(filename=filename) as data:
            filename_out = data.read_write_request(filename_out=filename_out,
                                                   request=request_variables,
                                                   cdm_dict=cdmtable)
        # logger.info("Result %s", str(filename_out))

    except Exception as e:
        logger.error('Exception %s occurred while reading %s', repr(e), filename)
        return '', 'exception "{}" occurred while reading {}'.format(e, filename)

    return filename_out, msg


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
