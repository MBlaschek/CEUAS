#!/usr/bin/env python3
__version__ = '0.2'
__author__ = 'LH'
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

import os
import copy
import glob
import io
import subprocess
import time
import zipfile
from datetime import datetime, timedelta
from functools import partial
from multiprocessing import Pool

import h5py
import matplotlib.pylab as plt
import numpy
from numba import *

import logging
logger = logging.getLogger('upperair.eua')

#
# Plotting functions for testing
#

def readandplot(rfile, body):
    with zipfile.ZipFile(rfile, 'r') as a:
        for r in a.filelist:
            try:
                with h5py.File(io.BytesIO(a.read(r)), 'r') as hf:
                    print(r.filename, hf.keys())
                    qs = 'select date,time,vertco_reference_1,obsvalue,fg_depar@body,biascorr@body where date=20190101 and time>=40000 and time<50000 and varno=2'
                    # qs='select * where date=20190101 and time>=40000 and time<50000 and varno=2'
                    odbfile = os.path.expandvars('$RSCRATCH/era5/odbs/1/era5.conv.201901.10393')
                    rdata = subprocess.check_output(["odb", "sql", "-q", qs, "-i", odbfile, '--no_alignment'])
                    npos = rdata.index(b'\n') + 1
                    xx = numpy.fromstring(rdata[npos:], dtype='float', sep='\t')
                    xx = xx.reshape((xx.shape[0] // 6, 6))
                    header = rdata[:npos].split()
                    check = {}
                    for i in range(len(header)):
                        check[header[i]] = xx[:, i]
                    header = rdata[:npos].split()
                    datetime(1900, 1, 1) + timedelta(days=int(hf['time'][0] // 86400),
                                                     seconds=int(hf['time'][0] % 86400))
                    if 'fbstats' in body.keys():
                        f, (ax1, ax2) = plt.subplots(1, 2)
                        if type(body['fbstats']) is not list:
                            body['fbstats'] = [body['fbstats']]
                        for fbs in body['fbstats']:
                            p = ax2.semilogy(hf[fbs], hf['plev'][:] / 100, label=fbs)

                        ax2.set_ylim(1030, 5)
                        ax2.legend()
                        ax2.set_title(r.filename)
                    else:
                        f, ax1 = plt.subplots(1, 1)

                    p = ax1.semilogy(hf['ta'], hf['plev'][:] / 100, label='ta')
                    ax1.set_ylim(1030, 5)
                    ax1.legend()
                    ax1.set_title(r.filename)
                    f.savefig(os.path.expanduser('~/plots/' + r.filename + '.png'))
                    plt.close()


            except MemoryError:
                pass
    return


def readandplot_ts(rfile, body):
    with zipfile.ZipFile(rfile, 'r') as a:
        for r in a.filelist:
            try:
                with h5py.File(io.BytesIO(a.read(r)), 'r') as hf:
                    print(r.filename, hf.keys())
                    qs = 'select date,time,vertco_reference_1,obsvalue,fg_depar@body,biascorr@body where varno=2 and vertco_reference_1=10000'
                    # qs='select * where date=20190101 and time>=40000 and time<50000 and varno=2'
                    odbfile = os.path.expandvars('$RSCRATCH/era5/odbs/1/era5.conv.201901.10393')
                    rdata = subprocess.check_output(["odb", "sql", "-q", qs, "-i", odbfile, '--no_alignment'])
                    npos = rdata.index(b'\n') + 1
                    rds = b'NaN'.join(rdata[npos:].split(b'NULL'))
                    xx = numpy.fromstring(rds, dtype='float', sep='\t')
                    xx = xx.reshape((xx.shape[0] // 6, 6))
                    header = rdata[:npos].split()
                    check = {}
                    for i in range(len(header)):
                        check[header[i]] = xx[:, i]
                    header = rdata[:npos].split()
                    datetime(1900, 1, 1) + timedelta(days=int(hf['time'][0] // 86400),
                                                     seconds=int(hf['time'][0] % 86400))
                    if 'fbstats' in body.keys():
                        f, (ax1, ax2) = plt.subplots(2, 1)
                        if type(body['fbstats']) is not list:
                            body['fbstats'] = [body['fbstats']]
                        for fbs in body['fbstats']:
                            mask = numpy.logical_and(hf['time'][:] % 86400 >= 9 * 3600,
                                                     hf['time'][:] % 86400 < 15 * 3600)
                            p = ax2.plot(hf['time'][mask] / 86400 / 365.25, hf[fbs][mask], label=fbs)

                        # ax2.set_ylim(1030,5)
                        ax2.legend()
                        ax2.set_title(r.filename)
                    else:
                        f, ax1 = plt.subplots(1, 1)

                    p = ax1.plot(hf['time'][mask] / 86400 / 365.25, hf['dew_point_temperature'][mask], label='td')
                    # ax1.set_ylim(1030,5)
                    ax1.legend()
                    ax1.set_title(r.filename)
                    f.savefig(os.path.expanduser('~/plots/' + r.filename + '.png'))
                    plt.close()


            except MemoryError:
                pass
    return


def secsince(t, funits, ref=datetime.strptime('1900-01-01 00:00:00', '%Y-%m-%d %H:%M:%S')):
    fstart = datetime.strptime(funits[-19:], '%Y-%m-%d %H:%M:%S')
    offset = fstart - ref
    offsets = offset.days * 24 * 3600 + offset.seconds
    fak = 1
    if 'hours' in funits:
        fak = 3600
    elif 'minutes' in funits:
        fak = 60

    secs = t * fak + offsets

    return secs


# @njit(cache=True)
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


def find_dateindex(y, x):
    """ creates the indices list from the dates, for quick access 
        nb the benchmark script will not work with these files since the definition of the array size is swapped i.e. (x.shape[0], 3)"""

    # x=y#numpy.unique(y)
    z = numpy.zeros((3, x.shape[0]), dtype=numpy.int32)
    z -= 1
    j = 0
    for i in range(len(y)):
        m = y[i]
        if x[j] == y[i]:
            if z[1, j] == -1:
                z[1, j] = i
                # print(j,i)
            else:
                if z[2, j] < i:
                    z[2, j] = i
        elif x[j] < y[i]:
            j += 1
            if x[j] == y[i]:
                if z[1, j] == -1:
                    z[1, j] = i
                    # print(j,i)
                else:
                    if z[2, j] < i:
                        z[2, j] = i
            else:
                print('Error')
        else:
            j -= 1
            if x[j] == y[i]:
                if z[1, j] == -1:
                    z[1, j] = i
                    # print(j,i)
                else:
                    if z[2, j] < i:
                        z[2, j] = i
            else:
                print('Error')
    z[0, :] = x
    return z


@njit
def find_dateindex_cg(y):
    x = numpy.unique(y)
    z = numpy.zeros((x.shape[0], 3), dtype=numpy.int32)
    z -= 1
    j = 0
    for i in range(len(y)):
        m = y[i]
        if x[j] == y[i]:
            if z[j, 0] == -1:
                z[j, 0] = i
                # print(j,i)
            else:
                if z[j, 1] < i:
                    z[j, 1] = i
        elif x[j] < y[i]:
            j += 1
            if x[j] == y[i]:
                if z[j, 0] == -1:
                    z[j, 0] = i
                    # print(j,i)
                else:
                    if z[j, 1] < i:
                        z[j, 1] = i
            else:
                print('Error')
        else:
            j -= 1
            if x[j] == y[i]:
                if z[j, 0] == -1:
                    z[j, 0] = i
                    # print(j,i)
                else:
                    if z[j, 1] < i:
                        z[j, 1] = i
            else:
                print('Error')
    z[:, 2] = x
    return z


def do_copy(fd, f, k, idx, cut_dimension,
            var_selection=[]):  # cuts vars and copies attributes of observation, feedback and header tables
    if not var_selection:
        var_selection = f[k].keys()
    for v in var_selection:
        if f[k][v].ndim == 1:
            if f[k][v].dtype != 'S1':

                fd[k].create_dataset_like(v, f[k][v], shape=idx.shape, chunks=True)
                hilf = f[k][v][idx[0]:idx[-1] + 1]
                fd[k][v][:] = hilf[idx - idx[0]]
            else:
                if v in [cut_dimension]:

                    fd[k].create_dataset_like(v, f[k][v], shape=idx.shape, chunks=True)
                    hilf = f[k][v][idx[0]:idx[-1] + 1]
                    fd[k][v][:] = hilf[idx - idx[0]]
                else:
                    fd[k].create_dataset_like(v, f[k][v])
                    fd[k][v][:] = f[k][v][:]
                    pass
        else:
            fd[k].create_dataset_like(v, f[k][v], shape=(idx.shape[0], f[k][v].shape[1]), chunks=True)
            hilf = f[k][v][idx[0]:idx[-1] + 1, :]
            fd[k][v][:] = hilf[idx - idx[0], :]
        for a in f[k][v].attrs.keys():
            if a not in ['DIMENSION_LIST', 'CLASS', 'external_table']:
                fd[k][v].attrs[a] = f[k][v].attrs[a]
    for v in var_selection:
        l = 0
        for d in f[k][v].dims:
            if len(d) > 0:
                logger.debug('%s/%s : %s', k, v, f[k][v].dims[l][0].name)
                fd[k][v].dims[l].attach_scale(fd[k][f[k][v].dims[l][0].name])
            l += 1


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


''' Main routine for parsing the CDS request, writing into flat netCDF file
c.retrieve(
    'reanalysis-era5-pressure-levels',
    {
        'product_type':'reanalysis',
        'format':'netcdf',
        'variable':[
            'geopotential','relative_humidity','specific_humidity',
            'specific_snow_water_content','temperature','u_component_of_wind',
            'v_component_of_wind'
        ],
        'pressure_level':[
            '450','650'
        ],
        'year':'1999',
        'month':'03',
        'day':'09',
        'time':'10:00'
    },
    'download.nc')
'''


@njit(cache=True)
def calc_trajindex(hh, hilf):
    # print(type(hh),type(hilf))
    j = 0  # hh[0]-hh[0]
    hilf[j] = hh[j]
    for i in range(hh.shape[0] - 1):
        x = j
        # print(i,hh[i+1],hh[i])
        if hh[i + 1] > hh[i]:
            j += 1
            hilf[j] = hh[i + 1]
        hh[i] = x
    hh[-1] = j

    # zid=numpy.zeros(len(zidx))
    # for i in range(len(zidx)):
    # zid[i]=zidx[i]

    return hilf[:j + 1]


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


def totimes(tinput):
    if type(tinput[0]) is str:
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


def ecdt(seconds, ref=datetime(1900, 1, 1)):
    x = ref + timedelta(seconds=numpy.int(seconds))
    return str(x)  # .year*10000+x.month*100+x.day)+' '+str(x.hour)+':'+str(x.minute)+':'+str(x.second)


@njit(cache=True)
def tohour(hilf, ohilf, dshift):
    ohilfold = -1
    for i in range(ohilf.shape[0]):
        if ohilfold == ohilf[i]:
            hilf[i] = hilf[i - 1]
            dshift[i] = dshift[i - 1]
        else:
            ohilfold = ohilf[i]
            hilf[i] = ((ohilf[i] + 3599) // 3600) % 24
            if hilf[i] == 0:
                dshift[i] = int(ohilf[i] % 3600 != 0)
                # if dshift[i]:
                # x=0
    return


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


def seconds_since_ref(idate, refdate):
    return (datetime(year=idate // 10000, month=idate % 10000 // 100, day=idate % 100) - refdate).days * 86400


def process_flat(outputdir: str, cftable: dict, datadir: str, request_variables: dict):
    """ Main Reading and writing Routine

    Args:
        outputdir: output directory
        cftable: CDM definitions table
        datadir: data directory
        request_variables: request variables

    Returns:
        str, str : File path, Error message
    """
    time0 = time.time()
    vdict = {}
    cdmdict = {}
    cdmnamedict = {}
    for igroup, v in cftable.items():
        if "odbcode" in v.keys():
            vdict[v['cdsname']] = v['odbcode']
            cdmdict[v['cdsname']] = v['cdmcode']
            cdmnamedict[v['cdsname']] = igroup

    error = ''
    filename = ''
    request_variables = copy.copy(request_variables)  # copy the dictionary
    statid = request_variables.pop('statid', None)
    request_keys = request_variables.keys()

    # cost=calculate_cost(rvars) # estimate size of output file

    if statid is None:
        logger.error('No station ID (statid) specified. %s', filename)
        return filename, 'No station ID (statid) specified'
    if statid[:3] == '0-2':
        suffix = ['']
    else:
        suffix = ['0-20000-0-', '0-20001-0-']

    for ss in suffix:
        filename = os.path.expandvars(datadir + '/' + ss + statid + '_CEUAS_merged_v0.nc')
        if os.path.isfile(filename):
            break

    rname = filename.split('/')[-1]
    logger.debug('Current: %s', filename)

    if len(request_keys) == 0:
        return filename, ''

    # for igroup in request_keys:
    #     if not isinstance(request_variables[igroup], list):
    #         if request_variables[igroup] not in cdmdict.keys():
    #             request_variables[igroup] = [request_variables[igroup], request_variables[igroup]]
    #         else:
    #             request_variables[igroup] = [cdmdict[request_variables[igroup]], cdmdict[request_variables[igroup]]]

    refdate = datetime(year=1900, month=1, day=1)
    try:
        with h5py.File(filename, 'r') as finput:
            time0 = time.time()
            # get time information from file
            # recordindex -> index of start position [start next[
            # recordtimestamp -> datetime seconds since for index position
            # Array [starttimestamp, index]
            dateindex = numpy.array((finput['recordtimestamp'][:],
                                     finput['recordindex'][:finput['recordtimestamp'].shape[0]]))
            try:
                request_date = request_variables.pop('date')
                #
                # format: [DATE], [START, END], [DATE, DATE, DATE, ...] (min 3)
                #
                if len(request_date) > 2:
                    dsec = []
                    for ievent in request_date:
                        dsec.append(seconds_since_ref(int(ievent), refdate))
                    dsec = numpy.asarray(dsec, dtype=numpy.int)
                else:
                    dsec = numpy.arange(seconds_since_ref(int(request_date[0]), refdate),
                                        seconds_since_ref(int(request_date[-1]), refdate) + 1,
                                        86400, dtype=numpy.int)

                if False:
                    # old version with ranges and date list
                    for ievent in request_date:
                        if '-' in ievent:
                            ievent = ievent.split('-')  # Date range from - to
                            dsec.extend(range(seconds_since_ref(int(ievent[0]), refdate),
                                              seconds_since_ref(int(ievent[-1]), refdate) + 1,
                                              86400))
                        else:
                            dsec.append(seconds_since_ref(int(ievent), refdate))

                # if time interval e.g. 21h-3h is chosen, the previous day must be extracted as well.
                prevday = 0
                if 'time' in request_keys:
                    tlist = totimes(request_variables['time'])
                    if tlist[0] > tlist[-1]:
                        prevday = 1

                if False:
                    if '-' in request_date[0]:
                        rvdp = request_date[0].split('-')
                    else:
                        rvdp = request_date
                    sdate = numpy.array(rvdp, dtype='int')

                    # if time interval e.g. 21h-3h is chosen, the previous day must be extracted as well.
                    prevday = 0
                    if 'time' in request_variables:
                        tlist = totimes(request_variables['time'])
                        if tlist[0] > tlist[-1]:
                            prevday = 1

                    dsec = []
                    for d in sdate:
                        dsec.append(((datetime(year=d // 10000, month=d % 10000 // 100,
                                               day=d % 100) - refdate)).days * 86400)

                    if '-' in request_date[0]:
                        dsec = numpy.arange(dsec[0], dsec[-1] + 1, 86400, dtype=numpy.int)

                    if len(dsec) == 1:
                        dsec = numpy.concatenate((dsec, dsec))

                dsec[-1] += 86399
                # request range [from - to]  (can be really large)
                # -86400 in order to get late launches of previous day
                didx = numpy.where(numpy.logical_and(dateindex[0, :] >= dsec[0] - 86400,
                                                     dateindex[0, :] <= dsec[-1]
                                                     )
                                   )[0]

                if didx.shape[0] == 0:
                    logger.warning('No data in time interval %s', rname)
                    return '', 'No data in specified time interval'

                didx = [didx[0], didx[-1]]
            except KeyError:  # all dates are considered
                didx = [0, dateindex.shape[1] - 1]
                request_date = ['']
                prevday = 0

            if didx[-1] + 1 == dateindex.shape[1]:
                trange = [dateindex[1, didx[0]], finput['observations_table']['observation_value'].shape[0]]  # Maximum
            else:
                trange = [dateindex[1, didx[0]], dateindex[1, didx[1] + 1]]  # Well within

            logger.debug('Datetime selection: %d - %d [%5.2f s] %s', trange[0], trange[1],
                         time.time() - time0, rname)
            mask = numpy.ones(trange[1] - trange[0], dtype=numpy.bool)
            # criteria = {'variable': 'era5fb/varno@body', 'level': 'era5fb/vertco_reference_1@body'}
            criteria = {'variable': 'observations_table/observed_variable',
                        'pressure_level': 'observations_table/z_coordinate',
                        'time': 'observations_table/date_time'}
            if True or '-' not in request_date[0]:
                criteria['date'] = 'observations_table/date_time'

            # time0 = time.time()
            for ckey, cval in criteria.items():
                if ckey in request_keys:
                    # mask=numpy.logical_and(mask,f[cv][trange[0]:trange[1]]>=rvdict[ck][0])
                    # mask=numpy.logical_and(mask,f[cv][trange[0]:trange[1]]<=rvdict[ck][1])
                    try:
                        if ckey == 'time':
                            us = numpy.string_(finput[cval].attrs['units'])
                            dh = us.split(b' ')[-1].split(b':')
                            if b'seconds' not in us:
                                logger.warning('Units not given in seconds, %s %s', us, rname)
                                return '', 'Units not given in seconds'

                            ohilf = finput[cval][trange[0]:trange[1]]
                            # add almost an hour (3600-1 sec) to account for difference between
                            # releasetime and nominal time. Time unit in CDS interface is hours.
                            hhilf = numpy.empty_like(ohilf, dtype=numpy.int32)
                            dshift = numpy.zeros_like(ohilf, dtype=numpy.int32)
                            hilf = numpy.empty_like(ohilf, dtype=numpy.int32)
                            tohourday(hhilf, hilf, ohilf, dshift)
                            # tohour(hhilf,ohilf,dshift)
                            # today(hilf,ohilf)
                            tlist = totimes(request_variables[ckey])
                            if prevday == 1:
                                ttlist = tlist[tlist < tlist[0]]  # use only late hours of the day before
                            else:
                                andisin(mask, hhilf, tlist)
                        elif ckey == 'date':
                            us = numpy.string_(finput[cval].attrs['units'])
                            dh = us.split(b' ')[-1].split(b':')
                            if b'seconds' not in us:
                                logger.warning('Units not given in seconds, %s %s', us, rname)
                                return '', 'Units not given in seconds'

                            if 'ohilf' not in locals():
                                ohilf = finput[cval][trange[0]:trange[1]]
                                hhilf = numpy.empty_like(ohilf)
                                dshift = numpy.zeros_like(ohilf, dtype=numpy.int32)
                                hilf = numpy.empty_like(ohilf)
                                tohourday(hhilf, hilf, ohilf, dshift)
                                tlist = totimes(['0-23'])
                            # tohour(hhilf,ohilf,dshift)
                            # today(hilf,ohilf)
                            else:
                                tlist = totimes(request_variables['time'])

                            dsec = numpy.array(dsec) // 86400
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
                                imask = numpy.logical_or(imask, imask2)
                                mask = numpy.logical_and(mask, imask)
                            else:
                                andisin_t(mask, hilf + dshift, dsec)
                        elif ckey == 'variable':
                            andisin(mask, finput[cval][trange[0]:trange[1]], numpy.int32(numpy.unique(cdmdict[request_variables[ckey]])))
                        else:
                            andisin(mask, finput[cval][trange[0]:trange[1]], numpy.int32(numpy.unique(request_variables[ckey])))
                        logger.debug('Finished %s [%5.2f s] %s', ckey, time.time() - time0, rname)

                    except MemoryError as e:
                        logger.error('Error %s occurred while checking criteria', repr(e))
                        return '', '"' + str(e) + '" occurred while checking criteria'

            logger.debug('Finished mask %s [%5.2f s]', rname, time.time() - time0)
            idx = numpy.where(mask)[0] + trange[0]  # make index for file subsetting
            if len(idx) == 0:
                logger.warning('No matching data found %s', rname)
                return '', 'No matching data found'
            logger.debug('Data found: %d %s', len(idx), rname)
            if False:
                # for debugging:
                fcold = 0
                fc = finput[criteria['time']][idx[0]:idx[-1] + 1]
                for i in idx:
                    e = fc[i - idx[0]]
                    if fcold != e:
                        # print(ecdt(e))
                        fcold = e

            trajectory_index = numpy.zeros_like(idx, dtype=numpy.int32)
            recordindex = dateindex[1, :]  # recordindex
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

            logger.debug('Request-keys: %s', str(request_keys))
            if 'variable' not in request_keys:
                logger.error('No variable specified %s %s', str(request_keys), rname)
                return '', 'No variable specified'

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
                    name_to_cdm[ss] = cftable[ss]
                except:
                    pass

            filename_out = outputdir + '/dest_' + statid + '_' + cdmnamedict[request_variables['variable']] + '.nc'
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
                        # groups that are simply copied
                        if False:
                            # for now commented out
                            fout.create_group(igroup)
                            for v in finput[igroup].keys():
                                finput.copy(finput[igroup][v], fout[igroup], name=v, without_attrs=True)
                            for a in finput[igroup].attrs.keys():
                                fout[igroup].attrs[a] = finput[igroup].attrs[a]
                            for v in finput[igroup].keys():
                                # print(k,v)
                                # fd[k].create_dataset_like(v,f[k][v])
                                for a in finput[igroup][v].attrs.keys():
                                    print(igroup, v, a)
                                    if a not in ['DIMENSION_LIST', 'CLASS']:
                                        if type(finput[igroup][v].attrs[a]) is str:
                                            fout[igroup][v].attrs[a] = numpy.string_(finput[igroup][v].attrs[a])
                                        else:
                                            fout[igroup][v].attrs[a] = finput[igroup][v].attrs[a]
                                l = 0
                                for d in finput[igroup][v].dims:
                                    if len(d) > 0:
                                        fout[igroup][v].dims[l].attach_scale(fout[igroup][finput[igroup][v].dims[l][0].name])
                                    l += 1
                            print(igroup, 'copied', time.time() - time0)
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
