#!/usr/bin/env python3
__version__ = '0.3'
__author__ = 'MB'
__status__ = 'dev'
__date__ = 'Mon Jul 13 10:29:21 CEST 2020'
__institute__ = 'UNIVIE'
__github__ = 'git@github.com:MBlaschek/CEUAS.git'
__doc__ = """
CDS_EUA Functions v%s
Maintained by %s at %s
Github: %s [%s]
Updated: %s
""" % (__version__, __author__, __institute__, __github__, __status__, __date__)

"""
New Version of cds_eua[] including much more documentation

Bugs:
 + Performance of make_selection is 6s slower than before ?
 
"""
import os
import copy
import glob
import time
import logging
import urllib.request
from datetime import datetime
from functools import partial
from multiprocessing import Pool

import h5py
import numpy
from numba import *

logger = logging.getLogger('upperair.eua')


def secsince(time_elapsed: int, time_unit: str, ref=datetime.strptime('1900-01-01 00:00:00', '%Y-%m-%d %H:%M:%S')):
    """ return seconds since a Reference date applying time unit

    Args:
        time_elapsed: seconds or other time unit
        time_unit: 'hours since 1900-01-01 00:00:00'
        ref: '1900-01-01 00:00:00'

    Returns:
        int : seconds since
    """
    fstart = datetime.strptime(time_unit[-19:], '%Y-%m-%d %H:%M:%S')
    offset = fstart - ref
    offsets = offset.days * 24 * 3600 + offset.seconds
    fak = 1
    if 'hours' in time_unit:
        fak = 3600
    elif 'minutes' in time_unit:
        fak = 60

    secs = time_elapsed * fak + offsets
    return secs


def read_standardnames(default_file='~/.tmp/cf.json') -> dict:
    """ Read CDM standard name of Variables

    Args:
        default_file: output filename, required by default.py

    Returns:
        dict : CDM definitons
    """
    import json
    import xml.etree.ElementTree as ET

    try:
        with open(os.path.expanduser(default_file)) as f:
            cf = json.load(f)
    except:
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
        with open(os.path.expanduser(default_file), 'w') as f:
            json.dump(cf, f)
    return cf


def do_cfcopy(fout: h5py.File, fin: h5py.File, group: str, idx: numpy.ndarray, cf: dict, dim0: str,
              var_selection: list = None):
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
                                                data=numpy.zeros(s1, dtype='S1'),
                                                chunks=True)
                            fout[sname].attrs['NAME'] = numpy.string_(
                                'This is a netCDF dimension but not a netCDF variable.')
                            fout[sname].make_scale(sname)
                        hilf = fin[group][v][idx[0]:idx[-1] + 1, :]
                        if hilf.shape[0] == 0:
                            print('x')
                        fout[vlist[-1]][:] = hilf[idx - idx[0], :]
                except:
                    # fix for missing report_id SHOULD BE REMOVED
                    hilf = numpy.zeros(shape=(idx.shape[0]), dtype='S10')
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
                                fout[vlist[-1]].attrs[a] = numpy.string_(fin[group][v].attrs[a])
                            else:
                                fout[vlist[-1]].attrs[a] = fin[group][v].attrs[a]

                    for a in cfv.keys():
                        if a not in ['shortname', 'odbcode', 'cdmcode']:
                            fout[vlist[-1]].attrs[a] = numpy.string_(cfv[a])
                        if a == 'units' and cfv[a] == 'NA':
                            fout[vlist[-1]].attrs[a] = numpy.string_('')
                        if a == 'units' and vlist[-1] == 'time':
                            ahilf = numpy.bytes_(fin[group][v].attrs[a])
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


@njit(cache=True)
def calc_trajindexfast(z, zidx, idx, trajectory_index):
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


def totimes(tinput):
    """ launch times to standard times ?

    Args:
        tinput: hours

    Returns:

    """
    if isinstance(tinput[0], str):
        if '-' in tinput[0]:
            ilist = numpy.array(tinput[0].split('-'), dtype=numpy.int32)
            if ilist[0] <= ilist[1]:
                ilist = numpy.arange(ilist[0], ilist[1] + 1)
            else:
                ilist = numpy.array(list(range(ilist[0], 24)) + list(range(ilist[1] + 1)), dtype=numpy.int32)
            out = ilist
        else:
            out = numpy.array(tinput, dtype=numpy.int32)
    else:
        out = numpy.array(tinput, dtype=numpy.int32)

    if numpy.min(out) < 0 or numpy.max(out) > 23:
        raise ValueError

    return out


@njit(cache=True)
def tohourday(hhilf, hilf, ohilf, dshift):
    """ Calculate what?

    Args:
        hhilf:
        hilf:
        ohilf: seconds array
        dshift:

    Returns:

    """
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
    return


@njit(cache=True)
def today(hilf, ohilf):
    ohilfold = -1
    for i in range(ohilf.shape[0]):
        if ohilf[i] == ohilfold:
            hilf[i] = hilf[i - 1]
        else:
            ohilfold = ohilf[i]
            hilf[i] = ohilf[i] // 86400
    return


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
def orisin(mask, x, v):
    for i in range(mask.shape[0]):
        if not mask[i]:
            found = False
            for j in range(v.shape[0]):
                if x[i] == v[j]:
                    found = True
                    break
            mask[i] = found


def make_selection(finput: h5py.File, request: dict, cdmdict: dict, rname,
                   observed_variable='observations_table/observed_variable',
                   z_coordinate='observations_table/z_coordinate',
                   date_time='observations_table/date_time',
                   debug=False):
    """

    Args:
        finput:
        request:
        cdmdict:
        rname:
        observed_variable:
        z_coordinate:
        date_time:
        debug:

    Returns:

    Time:
        Long request takes about 17 s
        +
    """
    time0 = time.time()
    # todo read reference date from unit string
    refdate = datetime(year=1900, month=1, day=1)
    # recordindex -> index of start position [start next[
    # recordtimestamp -> datetime seconds since for index position
    # Array [starttimestamp, index]
    dateindex = numpy.array((finput['recordtimestamp'][()],
                             finput['recordindex'][:finput['recordtimestamp'].shape[0]]))
    if request.get('date', None) is not None:
        #
        # format: [DATE], [START, END], [DATE, DATE, DATE, ...] (min 3)
        #
        if len(request['date']) > 2:
            dsec = []
            for ievent in request['date']:
                dsec.append(seconds_since_ref(int(ievent), refdate))
            dsec = numpy.asarray(dsec, dtype=numpy.int)
        else:
            dsec = numpy.arange(seconds_since_ref(int(request['date'][0]), refdate),
                                seconds_since_ref(int(request['date'][-1]), refdate) + 1,
                                86400, dtype=numpy.int)

        # if time interval e.g. 21h-3h is chosen, the previous day must be extracted as well.
        prevday = 0
        if request.get('time', None) is not None:
            tlist = totimes(request['time'])
            if tlist[0] > tlist[-1]:
                prevday = 1

        dsec[-1] += 86399
        # request range [from - to]  (can be really large)
        # -86400 in order to get late launches of previous day
        didx = numpy.where(numpy.logical_and(dateindex[0, :] >= dsec[0] - 86400,
                                             dateindex[0, :] <= dsec[-1]
                                             )
                           )[0]

        if didx.shape[0] == 0:
            logger.warning('No data in time interval')
            return '', 'No data in specified time interval'

        didx = [didx[0], didx[-1]]
    else:
        # all dates are considered
        didx = [0, dateindex.shape[1] - 1]  # todo not sure if that value is correct
        prevday = 0
        # dsec = numpy.arange(0, -1, 86400, dtype=numpy.int)  # todo wrong values here

    if didx[-1] + 1 == dateindex.shape[1]:
        trange = [dateindex[1, didx[0]], finput['observations_table']['observation_value'].shape[0]]  # Maximum
    else:
        trange = [dateindex[1, didx[0]], dateindex[1, didx[1] + 1]]  # Well within

    logger.debug('Mask recordtime: %d - %d [%5.2f s] %s', trange[0], trange[1], time.time() - time0, rname)
    mask = numpy.ones(trange[1] - trange[0], dtype=numpy.bool)
    try:
        if 'variable' in request.keys():

            if False:
                # almost no speed up
                mask = (finput[observed_variable][trange[0]:trange[1]] == numpy.int32(numpy.unique(cdmdict[request['variable']])))
            else:
                andisin(mask,
                        finput[observed_variable][slice(*trange)],
                        numpy.int32(numpy.unique(cdmdict[request['variable']]))
                        )
            logger.debug('Mask variable [%5.2f s] %s', time.time() - time0, rname)

        if 'pressure_level' in request.keys():
            # pressure levels
            andisin(mask, finput[z_coordinate][slice(*trange)],
                    numpy.int32(numpy.unique(request['pressure_level'])))

            logger.debug('Mask pressure_level [%5.2f s] %s', time.time() - time0, rname)

        ohilf = finput[date_time][slice(*trange)]

        if 'time' in request.keys():
            us = numpy.string_(finput[date_time].attrs['units'])
            # dh = us.split(b' ')[-1].split(b':')
            if b'seconds' not in us:
                logger.warning('Units not given in seconds, %s %s', us, rname)
                return '', 'Units not given in seconds'

            # ohilf = finput[date_time][trange[0]:trange[1]]
            # add almost an hour (3600-1 sec) to account for difference between
            # releasetime and nominal time. Time unit in CDS interface is hours.
            hhilf = numpy.empty_like(ohilf, dtype=numpy.int32)
            dshift = numpy.zeros_like(ohilf, dtype=numpy.int32)
            hilf = numpy.empty_like(ohilf, dtype=numpy.int32)
            tohourday(hhilf, hilf, ohilf, dshift)
            tlist = totimes(request['time'])
            if prevday == 1:
                # todo what is the meaning of this, not used anymore
                tlist = tlist[tlist < tlist[0]]  # use only late hours of the day before
            else:
                andisin(mask, hhilf, tlist)
            logger.debug('Mask time [%5.2f s] %s', time.time() - time0, rname)

        if 'date' in request.keys():
            us = numpy.string_(finput[date_time].attrs['units'])
            # dh = us.split(b' ')[-1].split(b':')
            if b'seconds' not in us:
                logger.warning('Units not given in seconds, %s', us)
                return '', 'Units not given in seconds'

            if 'time' not in request.keys():
                # ohilf = finput[date_time][trange[0]:trange[1]]
                hhilf = numpy.empty_like(ohilf)
                dshift = numpy.zeros_like(ohilf, dtype=numpy.int32)
                hilf = numpy.empty_like(ohilf)
                tohourday(hhilf, hilf, ohilf, dshift)
                tlist = totimes(['0-23'])   # All times allowed
            else:
                tlist = totimes(request['time'])  # Use only these times

            dsec = dsec // 86400  # request dsec for ???
            if prevday == 1:
                logger.debug('Mask selecting previous day %s', rname)
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
                imask = numpy.logical_or(imask, imask2)
                mask = numpy.logical_and(mask, imask)
            else:
                andisin_t(mask, hilf + dshift, dsec)
            logger.debug('Mask date [%5.2f s] %s', time.time() - time0, rname)

    except MemoryError as e:
        if debug:
            raise e
        logger.error('Error %s occurred while checking criteria', repr(e))
        return '"' + str(e) + '" occurred while checking criteria', ''
    return mask, trange


def seconds_since_ref(idate, refdate):
    return (datetime(year=idate // 10000, month=idate % 10000 // 100, day=idate % 100) - refdate).days * 86400


def process_flat(wroot: str, randdir: str, cdmtable: dict, request_variables: dict, debug: bool = False):
    """ Main Reading and writing Routine

    1. READ -> indices
    2. WRITE -> trajectory

    Args:
        wroot: path to public workdir
        randdir: path to random request dir
        cdmtable: CDM definitions table
        request_variables: {statid : ..., variable: ..., date: ..., time: ..., pressure_level: ...}
        debug: raise Errors or proceed?
    Returns:
        str, str : File path, Error message
    """
    vdict = {}
    cdmdict = {}
    cdmnamedict = {}
    for igroup, v in cdmtable.items():
        if "odbcode" in v.keys():
            vdict[v['cdsname']] = v['odbcode']
            cdmdict[v['cdsname']] = v['cdmcode']
            cdmnamedict[v['cdsname']] = igroup

    filename = ''
    request_variables = copy.copy(request_variables)  # copy the dictionary
    statid = request_variables.pop('statid', None)
    request_keys = request_variables.keys()

    # cost=calculate_cost(rvars) # estimate size of output file

    if statid is None:
        logger.error('No station ID (statid) specified. %s', filename)
        if debug:
            raise ValueError('No station ID (statid) specified')
        return filename, 'No station ID (statid) specified'

    if statid[:3] == '0-2':
        suffix = ['']
    else:
        suffix = ['0-20000-0-', '0-20001-0-']

    for ss in suffix:
        filename = os.path.expandvars('$RSCRATCH/era5/odbs/merged/' + ss + statid + '_CEUAS_merged_v0.nc')
        if os.path.isfile(filename):
            break

    rname = filename.split('/')[-1]
    logger.debug('Current: %s', filename)

    if len(request_keys) == 0:
        if debug:
            raise ValueError('Missing variable?', str(request_variables))
        return filename, ''

    try:
        with h5py.File(filename, 'r') as finput:
            time0 = time.time()
            # 1. Make Selection
            # 1.1. Datetime selection based on recordtimestamp
            # 1.2. Variable
            # 1.3. Pressure levels
            # 1.4. Datetime + time selection (final)
            # 1.5. Mask

            mask, trange = make_selection(finput,
                                          request_variables,
                                          cdmdict,
                                          rname,
                                          observed_variable='observations_table/observed_variable',
                                          z_coordinate='observations_table/z_coordinate',
                                          date_time='observations_table/date_time',
                                          debug=debug
                                          )
            if isinstance(mask, str):
                return '', mask  # Error (Memory Error)
            logger.debug('Make selection [%5.2f s] %s', time.time() - time0, rname)

            idx = numpy.where(mask)[0] + trange[0]  # make index for file subsetting
            if len(idx) == 0:
                logger.warning('No matching data found %s', rname)
                return '', 'No matching data found'

            logger.debug('Data found: %d %s', len(idx), rname)
            trajectory_index = numpy.zeros_like(idx, dtype=numpy.int32)
            recordindex = finput['recordindex'][()]  # recordindex
            zidx = numpy.where(numpy.logical_and(recordindex >= trange[0], recordindex < trange[1]))[0]
            recordindex = recordindex[zidx]
            zidx = calc_trajindexfast(recordindex, zidx, idx, trajectory_index)

            dims = {'obs': numpy.zeros(idx.shape[0], dtype=numpy.int32),
                    'trajectory': numpy.zeros(zidx.shape[0], dtype=numpy.int32)}
            globatts = {'Conventions': "CF-1.7",
                        'source': "radiosonde",
                        'featureType': "trajectory"}

            snames = ['report_id', 'platform_id', 'platform_name', 'observation_value', 'latitude',
                      'longitude', 'time', 'air_pressure', 'trajectory_label']

            logger.debug('Request-keys: %s', str(list(request_keys)))
            if 'variable' not in request_keys:
                logger.error('No variable specified %s %s', str(request_keys), rname)
                return '', 'No variable specified'

            logger.debug('Variable: %s', str(request_variables['variable']))
            if isinstance(request_variables['variable'], list):
                snames.append(cdmnamedict[request_variables['variable'][0]])
            else:
                snames.append(cdmnamedict[request_variables['variable']])

            if 'fbstats' in request_keys:
                if isinstance(request_variables['fbstats'], list):
                    for c in request_variables['fbstats']:
                        snames.append(cdmnamedict[c])
                else:
                    snames.append(cdmnamedict[request_variables['fbstats']])
            name_to_cdm = {}
            for ss in snames:
                try:
                    name_to_cdm[ss] = cdmtable[ss]
                except:
                    pass

            filename_out = wroot + '/' + randdir + '/dest_' + statid + '_' + cdmnamedict[
                request_variables['variable']] + '.nc'
            logger.debug('Writing: %s', filename_out)
            with h5py.File(filename_out, 'w') as fout:
                i = 0
                for d, v in dims.items():
                    fout.create_dataset(d, data=v)
                    fout[d].attrs['NAME'] = numpy.string_('This is a netCDF dimension but not a netCDF variable.')
                    fout[d].make_scale(d)  # resolves phony_dim problem
                    i += 1
                fout.create_dataset('trajectory_index', data=trajectory_index)
                fout['trajectory_index'].attrs['long_name'] = numpy.string_(
                    "index of trajectory this obs belongs to")
                fout['trajectory_index'].attrs['instance_dimension'] = numpy.string_("trajectory")
                fout['trajectory_index'].attrs['coordinates'] = numpy.string_("lat lon time plev")
                for igroup in finput.keys():
                    if not isinstance(finput[igroup], h5py.Group):
                        continue

                    if igroup in ['observations_table']:
                        # only obs, feedback fitting criteria (idx) is copied
                        do_cfcopy(fout, finput, igroup, idx, name_to_cdm, 'obs',
                                  var_selection=['observation_id', 'latitude', 'longitude', 'z_coordinate',
                                                 'observation_value', 'date_time'])
                        # 'observed_variable','units'
                        logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)
                    elif igroup in ['era5fb']:
                        # only obs, feedback fitting criteria (idx) is copied
                        if 'fbstats' in request_keys:
                            try:
                                do_cfcopy(fout, finput, igroup, idx, name_to_cdm, 'obs',
                                          var_selection=['fg_depar@body', 'an_depar@body',
                                                         'biascorr@body'])
                                # ['vertco_reference_1@body','obsvalue@body','fg_depar@body'])
                                logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)
                            except KeyError:
                                return '', 'no ' + name_to_cdm[request_variables['fbstats'][0]][
                                    'cdmname'] + ' found in ' + filename
                    elif igroup in ['header_table']:
                        # only records fitting criteria (zidx) are copied
                        do_cfcopy(fout, finput, igroup, zidx, name_to_cdm, 'trajectory',
                                  var_selection=['report_id'])
                        logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)
                        # ,'station_name','primary_station_id'])
                        # todo could be read from the observations_table
                    elif 'station_configuration' in igroup:
                        # only records fitting criteria (zidx) are copied
                        try:
                            sh = finput[igroup]['primary_id'].shape[1]
                            fout.attrs['primary_id'] = finput[igroup]['primary_id'][0].view('S{}'.format(sh))[0]
                            sh = finput[igroup]['station_name'].shape[1]
                            fout.attrs['station_name'] = finput[igroup]['station_name'][0].view('S{}'.format(sh))[0]
                        except:
                            logger.warning('No primary_id in %s', filename_out)
                        logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)

                    else:
                        pass
                fout['trajectory_label'].attrs['cf_role'] = numpy.string_('trajectory_id')
                fout['trajectory_label'].attrs['long_name'] = numpy.string_('Label of trajectory')
                for a, v in globatts.items():
                    fout.attrs[a] = numpy.string_(v)

                fout.attrs['history'] = numpy.string_(
                    'Created by Copernicus Early Upper Air Service Version 0, ' + datetime.now().strftime(
                        "%d-%b-%Y %H:%M:%S"))
                fout.attrs['license'] = numpy.string_('https://apps.ecmwf.int/datasets/licences/copernicus/')
        logger.debug('Finished %s [%5.2f s]', rname, time.time() - time0)
        return filename_out, ''
    except Exception as e:
        if debug:
            raise e
        logger.error('Exception %s occurred while reading %s', repr(e), filename)
        return '', 'exception "' + str(e) + '" occurred while reading ' + filename


if __name__ == '__main__':
    read_standardnames()
    os.chdir(os.path.expandvars('$HOME/python/web2py'))
    print(('called with ', sys.argv[1], sys.argv[2]))
    rvars = eval(sys.argv[2])
    wroot = os.path.expandvars('$RSCRATCH/tmp/')
    try:

        if '[0' in rvars['statid']:
            rvars['statid'] = "['" + rvars['statid'][1:7] + "']"
    except:
        pass

    # df=cdmexplainer(rvars)
    # print(df)

    os.chdir('/raid60/scratch/leo/scratch/era5/odbs/1')
    body = rvars
    bodies = []
    if type(body['statid']) is list:
        for s in body['statid']:
            bodies.append(dict(body))
            bodies[-1]['statid'] = s
    else:
        if body['statid'] == 'all':
            slist = glob.glob('/raid60/scratch/leo/scratch/era5/odbs/1/chera5.conv._?????.nc')
            statids = []
            for s in slist:
                bodies.append(dict(body))
                bodies[-1]['statid'] = s[-8:-3]
        else:
            bodies.append(dict(body))
    p = Pool(20)
    cf = read_standardnames()
    randdir = '{:012d}'.format(100000000000)
    try:
        os.mkdir(wroot + '/' + randdir)
    except:
        pass
    func = partial(process_flat, randdir, cf)
    t = time.time()
    results = list(map(func, bodies))
    # t=time.time()
    # rfile,error=process_flat(rvars)
    print(results)
    print(time.time() - t)
    print(results[0])
