#!/usr/bin/env python
# -*- coding: utf-8 -*-

__doc__ = """

formerly known as numpyquery.py -> LH

modified by MB
Modifications:
    + names of variables
    + attributes, metadata
    + some functions
    + interpolation (still broken)
    + MAIN
    
Need:
    + rewrite functions and remove unecessary stuff
    
Example:
    python odb_netcdf_format.py [odb netzcdf file]

Works on SRVX 8,2

Last used: Fr Mär  1 18:25:14 CET 2019

"""

import gzip
import os.path
import time
from datetime import datetime
from multiprocessing import Pool

import numpy
import pandas
import xarray
import numba
from numba import *

_metadata = {'t': {'units': 'K', 'standard_name': 'air_temperature'},
             't_fg_dep': {'units': 'K', 'standard_name': 'air_temperature_first_guess_departure'},
             't_an_dep': {'units': 'K', 'standard_name': 'air_temperature_analysis_departure'},
             't_bias': {'units': 'K', 'standard_name': 'air_temperature_bias_adjustment'},
             'rh': {'units': '1', 'standard_name': 'relative_humidity'},
             'rh_fg_dep': {'units': '1', 'standard_name': 'relative_humidity_first_guess_departure'},
             'rh_an_dep': {'units': '1', 'standard_name': 'relative_humidity_analysis_departure'},
             'rh_bias': {'units': '1', 'standard_name': 'relative_humidity_bias_adjusstment'},
             'q': {'units': 'kg/kg', 'standard_name': 'specific_humidity'},
             'q_fg_dep': {'units': 'kg/kg', 'standard_name': 'specific_humidity_first_guess_departure'},
             'q_an_dep': {'units': 'kg/kg', 'standard_name': 'specific_humidity_analysis_departure'},
             'q_bias': {'units': 'kg/kg', 'standard_name': 'specific_humidity_bias_adjustment'},
             'td': {'units': 'K', 'standard_name': 'dew_point'},
             'td_fg_dep': {'units': 'K', 'standard_name': 'dew_point_first_guess_departure'},
             'td_an_dep': {'units': 'K', 'standard_name': 'dew_point_analysis_departure'},
             'td_bias': {'units': 'K', 'standard_name': 'dew_point_bias_adjustment'},
             'u': {'units': 'm/s', 'standard_name': 'eastward_wind'},
             'u_fg_dep': {'units': 'm/s', 'standard_name': 'eastward_wind_first_guess_departure'},
             'u_an_dep': {'units': 'm/s', 'standard_name': 'eastward_wind_analysis_departure'},
             'u_bias': {'units': 'm/s', 'standard_name': 'eastward_wind_bias_adjustment'},
             'v': {'units': 'm/s', 'standard_name': 'northward_wind'},
             'v_fg_dep': {'units': 'm/s', 'standard_name': 'northward_wind_first_guess_departure'},
             'v_an_dep': {'units': 'm/s', 'standard_name': 'northward_wind_analysis_departure'},
             'v_bias': {'units': 'm/s', 'standard_name': 'northward_wind_bias_adjustment'},
             'dpd': {'units': 'K', 'standard_name': 'dew_point_depression'},
             'z': {'units': 'gpm', 'standard_name': 'geopotential_height'},
             'z_fg_dep': {'units': 'gpm', 'standard_name': 'geopotential_heighte_first_guess_departure'},
             'z_an_dep': {'units': 'gpm', 'standard_name': 'geopotential_height_analysis_departure'},
             'z_bias': {'units': 'gpm', 'standard_name': 'geopotential_height_bias_adjustment'},
             }


@njit(cache=True)
def fselect(x, y, z):
    l = []
    for i in range(x.shape[0]):
        if y[i] == 50000:
            if z[i] == 2:
                l.append(x[i])

    return numpy.array(l)


@njit(cache=True)
def std_date_time(idate):
    itime = idate % 1e6
    idate = idate // 1e6

    if itime < 60000:
        # 0 - 6
        return numpy.array((idate, 0))
    elif itime < 180000:
        # 0 - 18
        return numpy.array((idate, 12))
    else:
        # 18 - 24 > next day
        idate += 1
        day = idate % 100
        if day >= 29:
            year = idate // 10000
            month = (idate % 10000) // 100
            if month == 2:
                if year % 4 != 0 or day == 30:
                    month = 3
                    day = 1
            elif month == 4 or month == 6 or month == 9 or month == 11:
                if day > 30:
                    month = month + 1
                    day = 1
            elif month == 12:
                if day > 31:
                    year = year + 1
                    month = 1
                    day = 1
            else:
                if day > 31:
                    month = month + 1
                    day = 1
            idate = year * 10000 + month * 100 + day
        return numpy.array((idate, 0))


@njit(cache=True)
def myintp(pnew, pobs, xobs):
    res = numpy.interp(numpy.log(pnew), numpy.log(pobs), xobs)
    pmin = numpy.nanmin(pobs)
    pmax = numpy.nanmax(pobs)
    return numpy.where((pnew <= pmax) & (pnew >= pmin), res, numpy.nan)


# @njit(cache=True)
def filldaily_int(idatetime, vco, varnos, obs, fg, bias, ps, varno, miss_val):
    n = idatetime.shape[0]
    dailyobs = numpy.empty((n, ps.shape[0]), dtype=numpy.float32)
    dailyobs.fill(numpy.nan)
    dailyfg = numpy.empty((n, ps.shape[0]), dtype=numpy.float32)
    dailyfg.fill(numpy.nan)
    dailybias = numpy.empty((n, ps.shape[0]), dtype=numpy.float32)
    dailybias.fill(numpy.nan)
    #
    # Loop over unique indices (datetime)
    #
    j = -1
    for i in range(idatetime.shape[0] - 1):
        #
        # Select by date (positions might be different)
        # (considering python numbering)
        #
        k = (varnos[j + 1:idatetime[i] + 1] == varno).sum()
        j = idatetime[i]  # for next
        if k < 2:
            continue

        iobs = numpy.empty((k), dtype=numpy.float32)
        iobs.fill(numpy.nan)
        ips = numpy.empty((k), dtype=numpy.float32)
        ips.fill(numpy.nan)
        ifg = numpy.empty((k), dtype=numpy.float32)
        ifg.fill(numpy.nan)
        ibias = numpy.empty((k), dtype=numpy.float32)
        ibias.fill(numpy.nan)
        ni = 0
        for k in range(j + 1, idatetime[i] + 1):
            if varnos[k] == varno:
                iobs[ni] = obs[k]
                ips[ni] = vco[k]
                ifg[ni] = fg[k]
                ibias[ni] = bias[k]
                ni += 1

        if numpy.isfinite(iobs).sum() < 2:
            continue
        #
        # Sort for interpolation
        #
        itx = numpy.argsort(ips)  # missing at the end
        ips = ips[itx]
        iobs = iobs[itx]
        ibias = ibias[itx]
        ifg = ifg[itx]

        for k in range(ni - 1):
            if ips[k] == ips[k + 1]:
                ips[k] = numpy.nan
                ni -= 1

        itx = numpy.argsort(ips)  # missing at the end
        ips = ips[itx]
        iobs = iobs[itx]
        ibias = ibias[itx]
        ifg = ifg[itx]

        #
        # Remove Duplicates
        #
        # nitx = (itx > 0)
        # for j in range(itx.shape[0] - 1):
        #    if vco[itx[j]] == vco[itx[j + 1]]:
        #        nitx[j] = ~nitx[j]
        # itx = itx[nitx]
        # print(vco[itx])
        # print(obs[itx])
        #
        # Interpolate with log(p)
        #
        dailyobs[i, :] = myintp(numpy.log(ps), numpy.log(ips[:ni]), iobs[:ni])
        dailyfg[i, :] = myintp(numpy.log(ps), numpy.log(ips[:ni]), ifg[:ni])
        dailybias[i, :] = myintp(numpy.log(ps), numpy.log(ips[:ni]), ibias[:ni])
        # dailyobs[i, :] = numpy.interp(numpy.log(ps), numpy.log(vco[itx]), obs[itx], left=numpy.nan, right=numpy.nan)
        # dailyfg[i, :] = numpy.interp(numpy.log(ps), numpy.log(vco[itx]), fg[itx], left=numpy.nan, right=numpy.nan)
        # dailybias[i, :] = numpy.interp(numpy.log(ps), numpy.log(vco[itx]), bias[itx], left=numpy.nan, right=numpy.nan)

    print('Total', n)
    return dailyobs, dailyfg, dailybias


@njit(cache=True)
def filldaily_012(ldate, ltime, vco, varnos, obs, fg, bias, ps, at, ati, varno, miss_val):
    n = numpy.sum(ltime[ati] >= 180000)
    atn = [at[0]]
    dailyobs = numpy.empty((2, ps.shape[0], ati.shape[0] + n), dtype=numpy.float32)
    dailyobs.fill(miss_val)
    dailyfg = numpy.empty((2, ps.shape[0], ati.shape[0] + n), dtype=numpy.float32)
    dailyfg.fill(miss_val)
    dailybias = numpy.empty((2, ps.shape[0], ati.shape[0] + n), dtype=numpy.float32)
    dailybias.fill(miss_val)
    dailyhours = numpy.empty((2, ati.shape[0] + n), dtype=numpy.int32)
    dailyhours.fill(-1)
    ocount = 0
    l = 0  #
    m = 0  # missing
    #
    # Loop over all indices
    #
    for i in range(ati.shape[0] - 2):
        # if l<i-1:
        # print(i,l)
        lfound = False
        #
        # Loop over all indices per day (including one extra day (because not really sorted))
        #
        for j in range(ati[i + 2] - ati[i]):
            atij = ati[i] + j  # index of section
            #
            # not out of bounds, only same day (ldate)
            #
            if atij < ltime.shape[0] and ldate[atij] == ldate[ati[i]]:
                #
                # select only that variable
                #
                if varnos[atij] == varno:
                    #
                    # Loop over pressure levels (find)
                    #
                    for ip in range(ps.shape[0] - 1, -1, -1):
                        #
                        # correct pressure level
                        #
                        if vco[atij] == ps[ip]:
                            #
                            # not missing
                            #
                            if obs[atij] != miss_val:
                                lfound = True
                                #
                                # add new date (??? < 18 Z)
                                #
                                if atn[-1] < at[i] and ltime[ati[i]] < 180000:
                                    atn.append(at[i])
                                    l += 1
                                #
                                # TODO: potential mixing of profiles (missing at one time, present at other)
                                # no soultion to this problem, due to iteration -> more sounding times
                                #
                                if ltime[atij] < 60000:
                                    if numpy.isnan(dailyobs[0, ip, l]):
                                        dailyobs[0, ip, l] = obs[atij]
                                        dailyfg[0, ip, l] = fg[atij]
                                        dailybias[0, ip, l] = bias[atij]
                                        dailyhours[0, l] = ltime[atij] // 10000
                                elif ltime[atij] < 180000:
                                    if numpy.isnan(dailyobs[1, ip, l]):
                                        dailyobs[1, ip, l] = obs[atij]
                                        dailyfg[1, ip, l] = fg[atij]
                                        dailybias[1, ip, l] = bias[atij]
                                        dailyhours[1, l] = ltime[atij] // 10000
                                else:
                                    if atn[-1] < at[i] + 1 and ltime[atij] >= 180000:
                                        lt = at[i] + 1
                                        day = lt % 100
                                        if day >= 29:
                                            year = lt // 10000
                                            month = (lt % 10000) // 100
                                            if month == 2:
                                                if year % 4 != 0 or day == 30:
                                                    month = 3
                                                    day = 1
                                            elif month == 4 or month == 6 or month == 9 or month == 11:
                                                if day > 30:
                                                    month = month + 1
                                                    day = 1
                                            elif month == 12:
                                                if day > 31:
                                                    year = year + 1
                                                    month = 1
                                                    day = 1
                                            else:
                                                if day > 31:
                                                    month = month + 1
                                                    day = 1
                                            lt = year * 10000 + month * 100 + day

                                        atn.append(lt)
                                        l += 1
                                    if numpy.isnan(dailyobs[0, ip, l]):
                                        dailyobs[0, ip, l] = obs[atij]
                                        dailyfg[0, ip, l] = fg[atij]
                                        dailybias[0, ip, l] = bias[atij]
                                        dailyhours[0, l] = ltime[atij] // 10000
                            # break pressure level loop
                            break

            if not lfound:
                m += 1
    print('l', l, 'm', m, ati.shape[0])
    return dailyobs[:, :, :l + 1], dailyfg[:, :, :l + 1], dailybias[:, :, :l + 1], numpy.array(atn), dailyhours[:,:l + 1]


@njit(cache=True)
def filldaily(ldate, ltime, vco, varnos, obs, fg, bias, ps, at, ati, varno, miss_val):
    n = numpy.sum(ltime[ati] >= 180000)
    atn = [at[0]]
    dailyobs = numpy.empty((4, ps.shape[0], ati.shape[0] + n), dtype=numpy.float32)
    dailyobs.fill(miss_val)
    dailyfg = numpy.empty((4, ps.shape[0], ati.shape[0] + n), dtype=numpy.float32)
    dailyfg.fill(miss_val)
    dailybias = numpy.empty((4, ps.shape[0], ati.shape[0] + n), dtype=numpy.float32)
    dailybias.fill(miss_val)
    dailyhours = numpy.empty((4, ati.shape[0] + n), dtype=numpy.float32)
    dailyhours.fill(miss_val)
    ocount = 0
    l = 0  #
    m = 0  # missing
    #
    # Loop over all indices
    #
    for i in range(ati.shape[0] - 2):
        # if l<i-1:
        # print(i,l)
        lfound = False
        #
        # Loop over all indices per day (including one extra day (because not really sorted))
        #
        for j in range(ati[i + 2] - ati[i]):
            atij = ati[i] + j   # index of section
            #
            # not out of bounds, only same day (ldate)
            #
            if atij < ltime.shape[0] and ldate[atij] == ldate[ati[i]]:
                #
                # select only that variable
                #
                if varnos[atij] == varno:
                    #
                    # Loop over pressure levels (find)
                    #
                    for ip in range(ps.shape[0] - 1, -1, -1):
                        #
                        # correct pressure level
                        #
                        if vco[atij] == ps[ip]:
                            #
                            # not missing
                            #
                            if obs[atij] != miss_val:
                                lfound = True
                                #
                                # add new date (??? < 18 Z)
                                #
                                if atn[-1] < at[i] and ltime[ati[i]] <= 210000:
                                    atn.append(at[i])
                                    l += 1
                                #
                                # 21 -> 03 == 0
                                # 03 -> 09 == 6
                                # 09 -> 15 == 12
                                # 15 -> 21 == 18
                                #
                                if ltime[atij] <= 30000:
                                    # 00 Z
                                    dailyobs[0, ip, l] = obs[atij]
                                    dailyfg[0, ip, l] = fg[atij]
                                    dailybias[0, ip, l] = bias[atij]
                                    dailyhours[0, l] = ltime[atij] // 10000
                                elif ltime[atij] <= 90000:
                                    # 06 Z
                                    dailyobs[1, ip, l] = obs[atij]
                                    dailyfg[1, ip, l] = fg[atij]
                                    dailybias[1, ip, l] = bias[atij]
                                    dailyhours[1, l] = ltime[atij] // 10000
                                elif ltime[atij] <= 150000:
                                    # 12 Z
                                    dailyobs[2, ip, l] = obs[atij]
                                    dailyfg[2, ip, l] = fg[atij]
                                    dailybias[2, ip, l] = bias[atij]
                                    dailyhours[2, l] = ltime[atij] // 10000
                                elif ltime[atij] <= 210000:
                                    # 18 Z
                                    dailyobs[1, ip, l] = obs[atij]
                                    dailyfg[1, ip, l] = fg[atij]
                                    dailybias[1, ip, l] = bias[atij]
                                    dailyhours[1, l] = ltime[atij] // 10000
                                else:
                                    # > 21 + one day
                                    if atn[-1] < at[i] + 1 and ltime[atij] > 210000:
                                        lt = at[i] + 1
                                        day = lt % 100
                                        if day >= 29:
                                            year = lt // 10000
                                            month = (lt % 10000) // 100
                                            if month == 2:
                                                if year % 4 != 0 or day == 30:
                                                    month = 3
                                                    day = 1
                                            elif month == 4 or month == 6 or month == 9 or month == 11:
                                                if day > 30:
                                                    month = month + 1
                                                    day = 1
                                            elif month == 12:
                                                if day > 31:
                                                    year = year + 1
                                                    month = 1
                                                    day = 1
                                            else:
                                                if day > 31:
                                                    month = month + 1
                                                    day = 1
                                            lt = year * 10000 + month * 100 + day

                                        atn.append(lt)
                                        l += 1
                                    dailyobs[0, ip, l] = obs[atij]
                                    dailyfg[0, ip, l] = fg[atij]
                                    dailybias[0, ip, l] = bias[atij]
                                    dailyhours[0, l] = ltime[atij] // 10000
                            # break pressure level loop
                            break

            if not lfound:
                m += 1
    print('ocount', ocount, 'l', l, 'm', m, ati.shape[0])
    return dailyobs[:, :, :l + 1], dailyfg[:, :, :l + 1], dailybias[:, :, :l + 1], numpy.array(atn), dailyhours[:,
                                                                                                     :l + 1]


@njit(cache=True)
def myunique(dates):
    #
    # Unique dates (at) and unique index of first date (ati)
    #
    at = numpy.zeros(45000, dtype=numpy.int32)
    ati = numpy.zeros(45000, dtype=numpy.int32)
    l = 0
    at[0] = dates[0]
    for i in range(1, dates.shape[0]):
        dat = dates[i]
        #
        # new date (later?)
        #
        if dat > at[l]:
            l += 1
            at[l] = dat
            ati[l] = i
        #
        # new date (older?)
        #
        elif dat < at[l]:
            #
            # is it new?
            #
            if dat > at[l - 1]:
                #
                # Switch places
                #
                at[l + 1] = at[l]
                ati[l + 1] = ati[l]
                at[l] = dat
                ati[l] = i
                l += 1
        else:
            pass

    return at[:l + 1], ati[:l + 1]


def make_xarray(fpv, template, dailyobs, dailyfgd, dailybias, atn, dailyhours, name, units):
    hpt = ['hour', 'pressure', 'time']

    spl = xarray.Dataset(data_vars={'lat': ('station', [fpv['lat@hdr'][0]]),
                                    'lon': ('station', [fpv['lon@hdr'][0]]),
                                    'alt': ('station', [fpv['stalt@hdr'][0]]),
                                    'press': ('pressure', template['press']),
                                    'datum': (('numdat', 'time'), [atn]),
                                    'hours': (('hour', 'time'), dailyhours),
                                    name: (hpt, dailyobs),
                                    'fg_dep': (hpt, dailyfgd),
                                    'bias': (hpt, dailybias),
                                    },
                         attrs=template.attrs)
    for v in [name, 'fg_dep', 'bias']:
        spl[v].attrs['units'] = units
    spl['datum'].attrs['units'] = 'days since 1900-01-01 0:0:0'
    spl = spl.assign_coords(lat=(spl.lat), lon=(spl.lon),
                            hour=(spl.hours[:, 0] - spl.hours[:, 0] + numpy.array([0, 12])), press=(spl.press),
                            datum=(spl.datum[0, :]))
    return spl


def process_obj(fpv, fn):
    t = time.time()
    ps = numpy.array(
        [1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 85000, 92500, 100000])

    ident = '0' + fn.split('.')[-2][1:]

    global_attrs = {'Conventions': 'CF-1.1',
                    'title': 'station daily temperature series',
                    'institution': 'University of Vienna',
                    'history': '19/03/20',
                    'source': 'radiosonde, ERA-5, ERA-Interim, ERA-40, RAOBCORE',
                    'references': 'www.univie.ac.at/theoret-met/research/raobcore',
                    'levels': 'plevs [%d -%d] #%d' % (min(ps), max(ps), len(ps)),
                    'libs': "NP(%s) XR(%s) NB(%s)" % (numpy.__version__, xarray.__version__, numba.__version__)}

    at, ati = myunique(fpv['date@hdr'].values)   # just dates
    alldates = []
    for a in at:
        alldates.append(datetime(a // 10000, (a % 10000) // 100, a % 100))
    alldates = numpy.array(alldates)

    td0 = time.time()
    #
    # too few data
    #
    if len(at) < 100:
        raise RuntimeError("Too few dates (<100)")

    variables = {2: 't',
                 3: 'u',
                 4: 'v',
                 59: 'td',
                 7: 'q',
                 29: 'rh',
                 # 1: 'z',
                 # 39: 't2m',
                 # 40: 'td2m',
                 # 41: 'u10m',
                 # 42: 'v10m',
                 # 58: 'rh2m',
                 # 111: 'dd',
                 # 112: 'ff'
                 }

    paras = list(variables.keys())
    names = list(variables.values())
    variables = list(fpv.data_vars)
    #
    # some specific variables present in ODB
    #
    default_drop = ["type", "expver", "class", "stream", "andate", "antime", "reportype", "numtsl@desc",
                    "timeslot@timeslot_index", "seqno@hdr", "source@hdr", "bufrtype@hdr", "subtype@hdr", "groupid@hdr",
                    "statid@hdr", "report_status@hdr", "report_event1@hdr", "report_rdbflag@hdr", "entryno@body",
                    "vertco_type@body", "ppcode@conv_body", "datum_anflag@body", "datum_status@body",
                    "datum_event1@body", "datum_rdbflag@body"]
    #
    # typical information that is only relevant for one sounding (repeated)
    #
    sel = [ivar for ivar in variables if '@hdr' in ivar or '@modsurf' in ivar or '@surfbody' in ivar or '@conv' in ivar]
    sel = [ivar for ivar in sel if ivar not in default_drop]
    #
    #
    #
    sel.remove('date@hdr')
    sel.remove('time@hdr')
    launches = None
    #
    # Make the directory for output
    #
    try:
        os.mkdir(ident)
    except:
        pass
    #
    # Iterate variables
    #
    for l, var in enumerate(paras):
        tdx = time.time()
        dailyobs, dailyfgd, dailybias, atn, dailyhours = filldaily(fpv['date@hdr'].values, fpv['time@hdr'].values,
                                                                   fpv[u'vertco_reference_1@body'].values,
                                                                   fpv[u'varno@body'].values,
                                                                   fpv[u'obsvalue@body'].values,
                                                                   fpv[u'fg_depar@body'].values,
                                                                   fpv[u'biascorr@body'].values, ps, at, ati, var,
                                                                   numpy.nan)
        index = []
        for a in atn:
            # Datetime as from variable -> will be aligned later
            index.append(datetime(a // 10000, (a % 10000) // 100, a % 100))

        #
        # Xarray
        #
        spl = xarray.Dataset()
        spl[names[l]] = (('hour', 'plev', 'time'), dailyobs)
        spl[names[l] + '_fg_dep'] = (('hour', 'plev', 'time'), dailyfgd)
        spl[names[l] + '_bias'] = (('hour', 'plev', 'time'), dailybias)

        spl['time'] = index
        spl['plev'] = ps
        spl['hour'] = numpy.array([0, 6, 12, 18])
        spl['plev'].attrs.update({'units': 'Pa', 'standard_name': 'air_pressure', 'axis': 'Z'})
        spl['time'].attrs.update({'axis': 'T'})
        spl['hour'].attrs.update({'units': 'h', 'standard_name': 'standard_launch_time'})
        #
        # Launch times
        #
        if l == 0 or launches is None:
            # once it's ok
            spl['launch'] = (('hour', 'time'), dailyhours)

        spl = spl.reindex(time=alldates)
        if l == 0 or launches is None:
            launches = spl['launch'].copy()
            spl = spl.drop('launch')
        #
        # Metadata
        #
        for v in [names[l], names[l] + '_fg_dep', names[l] + '_bias']:
            spl[v].attrs.update(_metadata[v])
        #
        # Attributes
        #
        spl.attrs.update(global_attrs)
        spl.attrs.update({'station_id': ident,
                          'station_lat': "%.2f N" % fpv['lat@hdr'].values[-1],
                          'station_lon': "%.2f E" % fpv['lon@hdr'].values[-1],
                          'station_alt': "%.1f m" % fpv['stalt@hdr'].values[-1]
                          })
        if 'hdrlen' in spl.coords:
            spl = spl.drop('hdrlen')
        #
        # IO
        #
        if '2402' in fn:
            fno = os.path.expandvars('./' + ident + '/ERAI_' + ident + '_' + names[l] + '.nc')
        else:
            fno = os.path.expandvars('./' + ident + '/ERA5_' + ident + '_' + names[l] + '.nc')

        spl.to_netcdf(fno)
        #
        #
        #
        td1 = time.time() - tdx
        print(ident, var, names[l], dailyobs.shape, td1)

    station = fpv[sel].isel(hdrlen=ati)
    station = station.rename({'hdrlen': 'time'})
    station = station.assign_coords(time=alldates)
    station = station.rename({i: i.replace('@hdr', '') for i in list(station.data_vars)})
    station = station.rename({i: i.replace('@modsurf', '_msurf') for i in list(station.data_vars)})
    station = station.rename({i: i.replace('@surfbody_feedback', '_surf_fb') for i in list(station.data_vars)})
    station = station.rename({i: i.split('@conv')[0] for i in list(station.data_vars)})

    if launches is not None:
        station['launch'] = launches

    station.attrs.update(global_attrs)
    # sort by index
    station = station.sortby('time')
    if '2402' in fn:
        station.to_netcdf('./' + ident + '/ERAI_' + ident + '_station.nc')
    else:
        station.to_netcdf('./' + ident + '/ERA5_' + ident + '_station.nc')

    td2 = time.time() - td0
    td3 = time.time() - t
    print(ident, "DATE: ", td0 - t, "VAR:", td2, "TOTAL:", td3, "COUNT:", fpv[u'obsvalue@body'].shape)


def doquery(fno):
    t = time.time()
    #
    #
    #
    if '.gz' in fno:
        fn = fno[:-3]  # remove .gz
    else:
        fn = fno
    ident = '0' + fn.split('.')[-2][1:]
    try:
        #
        # ASCII ODB Dump
        #
        if 'ascii' in fno:
            data = pandas.read_csv(sys.argv[1], sep='\t', error_bad_lines=False, engine='c',
                                quotechar="'", low_memory=False, skipinitialspace=True)
            data.index.name = 'hdrlen'  # compatibility
            data.columns = data.columns.str.strip()  # can happen
            g = data.to_xarray()

        elif '.gz' in fno:
            gz = gzip.open(fno)
        #
        # XArray (requires a netcdf3 - 64 bit version)
        #
            g = xarray.open_dataset(gz)
        else:
            g = xarray.open_dataset(fno, chunks={'hdrlen': 50000})
        #
        # Rearange
        #
        print(ident, time.time() - t, "Start processing")
        process_obj(g, fn)
        #
        #
        #
    except Exception as e:
        print(ident, "Error", repr(e))
    else:
        g.close()
        #
        #
        #
    print(ident, time.time() - t)


if __name__ == '__main__':
    import sys
    import os

    if len(sys.argv) == 2:
        print(sys.argv)
        doquery(sys.argv[1])
    elif len(sys.argv) > 2:
        # fns = '/raid8/srvx1/mblaschek/tmp/gzipped/era5.2402.conv.*.nc.gz'
        # fns = '/raid60/scratch/leo/scratch/era5/odbs/1/era5.conv.*.nc.gz'
        # files = glob.glob(fns)
        files = sys.argv[1:]

        for ifile in files[:]:
            if not os.path.isfile(ifile):
                files.remove(ifile)
                # print(ifile)

        print("ERA-5 Stations:", len(files))
        print("Starting a Pool ...")
        p = Pool(20)
        res = p.map(doquery, files)
        print("Done")
    else:
        print("Usage: numpyquery.py [files]")
