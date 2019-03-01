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
def standard_sounding_times(data, dim='date', times=(0, 12), span=6, freq='12h', return_indices=False, only_time=False,
                            **kwargs):
    """ Standardize datetime to times per date, try to fill gaps

    Args:
        data (xarray.DataArray): Input DataArray
        dim (str): datetime dimension
        times (tuple): sounding times
        span (int): plus minus sounding time for filling gaps
        freq (str): frequency of sounding times
        return_indices (bool): return indices for alignment
        only_time (bool): when multiple times can be set std times, use only time difference for selection

    Returns:
        xarray.DataArray : datetime standardized DataArray
    """

    kwargs['mname'] = kwargs.get('mname', 'std_hours')

    if not isinstance(data, xarray.DataArray):
        raise ValueError('Requires a DataArray', type(data))

    if dim not in data.dims:
        raise ValueError('Requires a datetime dimension', dim)

    dates = data[dim].values.copy()

    #  all possible dates
    alldates = pandas.date_range(pandas.Timestamp(dates.min()).replace(hour=np.min(times), minute=0, second=0),
                                 pandas.Timestamp(dates.max()).replace(hour=np.max(times), minute=0, second=0),
                                 freq=freq)

    new = data.reindex(**{dim: alldates})  # complete reindex to new dates (implicit copy)
    new['delay'] = (dim, np.zeros(alldates.size))  # new coordinate for delays
    # matching dates
    new_logic = np.isin(new[dim].values, dates)
    # not found directly
    jtx = np.where(~new_logic)[0]  # newdates not fitting dates (indices)
    new['delay'].values[jtx] = np.nan  # not fitting dates
    newindex = [slice(None)] * data.ndim
    oldindex = [slice(None)] * data.ndim
    axis = data.dims.index(dim)
    nn = 0
    indices = []
    # All times not yet filled
    # Is there some data that fits within the given time window
    for itime in new[dim].values[jtx]:
        diff = (itime - dates) / np.timedelta64(1, 'h')  # closest sounding
        n = np.sum(np.abs(diff) <= span)  # number of soundings within time window
        if n > 0:
            i = np.where(alldates == itime)[0]  # index for new array
            if n > 1 and not only_time:
                # many choices, count data
                k = np.where(np.abs(diff) <= span)[0]
                oldindex[axis] = k
                # count data of candidates
                # weight by time difference (assuming, that a sounding at the edge of the window is less accurate)
                distance = np.abs(diff[k])
                counts = np.sum(np.isfinite(data.values[tuple(oldindex)]), axis=1) / np.where(distance != 0, distance,
                                                                                              1)
                if np.any(counts > 0):
                    j = k[np.argmax(counts)]  # use the one with max data / min time diff
                else:
                    continue

            else:
                j = np.argmin(np.abs(diff))  # only one choice

            newindex[axis] = i
            oldindex[axis] = j
            counts = np.sum(np.isfinite(data.values[tuple(oldindex)]), axis=0)
            if counts > 0:
                new.values[tuple(newindex)] = data.values[tuple(oldindex)]  # update data array
                # new['delay'].values[i] = -1 * diff[j]  # pd.Timestamp(dates[j]).hour  # datetime of minimum
                # indices += [(i[0], j)]
                # nn += 1

    new.attrs['std_times'] = str(times)
    # new['delay'].attrs['updated'] = nn
    # new['delay'].attrs['missing'] = new['delay'].isnull().sum().values
    # new['delay'].attrs['times'] = str(times)
    # if return_indices:
    #     return new, np.array(indices)
    return new


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
def filldaily(ldate, ltime, vco, varnos, obs, fg, bias, ps, at, ati, varno, miss_val):
    n = numpy.sum(ltime[ati] >= 180000)
    atn = [at[0]]
    dailyobs = numpy.empty((2, ps.shape[0], ati.shape[0] + n), dtype=numpy.float32)
    dailyobs.fill(numpy.nan)
    dailyfg = numpy.empty((2, ps.shape[0], ati.shape[0] + n), dtype=numpy.float32)
    dailyfg.fill(numpy.nan)
    dailybias = numpy.empty((2, ps.shape[0], ati.shape[0] + n), dtype=numpy.float32)
    dailybias.fill(numpy.nan)
    dailyhours = numpy.empty((2, ati.shape[0] + n), dtype=numpy.int32)
    dailyhours.fill(numpy.nan)
    ocount = 0
    iold = 0
    l = 0  #
    m = 0  # missing
    for i in range(ati.shape[0] - 2):
        # if l<i-1:
        # print(i,l)
        lfound = False
        for j in range(ati[i + 2] - ati[i]):
            atij = ati[i] + j
            if atij < ltime.shape[0] and ldate[atij] == ldate[ati[i]]:
                if varnos[atij] == varno:
                    #
                    # Loop over pressure levels (find)
                    #
                    for ip in range(ps.shape[0] - 1, -1, -1):
                        if vco[atij] == ps[ip]:
                            if obs[atij] != miss_val:
                                lfound = True
                                #
                                if atn[-1] < at[i] and ltime[ati[i]] < 180000:
                                    atn.append(at[i])
                                    l += 1

                                if ltime[atij] < 60000:
                                    # if dailyspl[0,ip,i]!=dailyspl[0,ip,i]:
                                    # dailyspl[0,ip,i]=obs[atij]
                                    # else:
                                    dailyobs[0, ip, l] = obs[atij]
                                    dailyfg[0, ip, l] = fg[atij]
                                    dailybias[0, ip, l] = bias[atij]
                                    dailyhours[0, l] = ltime[atij] // 10000
                                # ocount+=1
                                elif ltime[atij] < 180000:
                                    # if dailyspl[1,ip,i]!=dailyspl[1,ip,i]:
                                    # dailyspl[1,ip,i]=obs[atij]
                                    # else:
                                    dailyobs[1, ip, l] = obs[atij]
                                    dailyfg[1, ip, l] = fg[atij]
                                    dailybias[1, ip, l] = bias[atij]
                                    dailyhours[1, l] = ltime[atij] // 10000
                                # ocount+=1
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
                                        # print(atn[-1])
                                        l += 1
                                    # !!
                                    # if dailyspl[0,ip,i+1]!=dailyspl[0,ip,i+1]:
                                    # dailyspl[0,ip,i+1]=obs[atij]
                                    # else:
                                    dailyobs[0, ip, l] = obs[atij]
                                    dailyfg[0, ip, l] = fg[atij]
                                    dailybias[0, ip, l] = bias[atij]
                                    dailyhours[0, l] = ltime[atij] // 10000
                                    # ocount+=1
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
        # if times[i]>=180000:
        # dat+=1
        # if dat==19870329.:
        # print(at[l],dates[i])
        if dat > at[l]:
            l += 1
            at[l] = dat
            ati[l] = i
        elif dat < at[l]:
            if dat > at[l - 1]:
                at[l + 1] = at[l]
                ati[l + 1] = ati[l]
                at[l] = dat
                ati[l] = i
                l += 1
                # print(at[l],dates[i])

        # if l==1897:
        # break

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

    # ps = numpy.array([1000., 2000., 3000., 5000., 7000., 10000., 12500., 15000., 17500., 20000., 22500., 25000., 30000.,
    #                   35000., 40000., 45000., 50000., 55000., 60000., 65000., 70000., 75000., 77500., 80000., 82500.,
    #                   85000., 87500., 90000., 92500., 95000., 97500., 100000.])

    ident = '0' + fn.split('.')[-2][1:]

    global_attrs = {'Conventions': 'CF-1.1',
                    'title': 'station daily temperature series',
                    'institution': 'University of Vienna',
                    'history': '19/01/05',
                    'source': 'radiosonde, ERA-5, ERA-Interim, ERA-40, RAOBCORE',
                    'references': 'www.univie.ac.at/theoret-met/research/raobcore',
                    'ident': ident,
                    'levels': 'plevs [%d -%d] #%d' % (min(ps), max(ps), len(ps)),
                    'libs': "NP(%s) XR(%s) NB(%s)" % (numpy.__version__, xarray.__version__, numba.__version__)}
    if False:
        # ,fpv['antime'].values)
        hdrdate = (fpv['date@hdr'].values * 100 + fpv['time@hdr'] / 1e4).astype(int)
        ati = numpy.flatnonzero(numpy.diff(hdrdate))  # end index of date groups
        ati = numpy.append(ati, hdrdate.size - 1)  # last
        at = hdrdate[ati]  # Dates
        alldates = []
        for i in range(at.shape[0]):
            a = int(at[i] // 100)
            b = at[i] % 1e2
            # Year, month, day, hour
            alldates.append(datetime(a // 10000, (a % 10000) // 100, a % 100, int(b)))
            # a,b = std_date_time(at[i])
            # stddates.append(datetime(int(a // 10000), int((a % 10000) // 100), int(a % 100), int(b)))

        alldates = numpy.array(alldates)
    else:
        at, ati = myunique(fpv['date@hdr'].values)   # just dates
        alldates = []
        for a in at:
            alldates.append(datetime(a // 10000, (a % 10000) // 100, a % 100))
        alldates = numpy.array(alldates)

    td0 = time.time()
    if len(at) < 100:
        raise RuntimeError("Too few dates (<100)")

    # names = ('temperatures', 'uwind', 'vwind', 'dewpoint')
    variables = {2: 't', 3: 'u', 4: 'v', 59: 'td', 7: 'q', 29: 'rh'}  # , 1: 'z',
    # {39: 't2m', 40: 'td2m', 41: 'u10m', 42: 'v10m', 58: 'rh2m', 111: 'dd', 112: 'ff'}
    # units = ('K', 'm/s', 'm/s', 'K')
    # paras = (2, 3, 4, 59)
    paras = list(variables.keys())
    names = list(variables.values())

    variables = list(fpv.data_vars)
    default_drop = ["type", "expver", "class", "stream", "andate", "antime", "reportype", "numtsl@desc",
                    "timeslot@timeslot_index", "seqno@hdr", "source@hdr", "bufrtype@hdr", "subtype@hdr", "groupid@hdr",
                    "statid@hdr", "report_status@hdr", "report_event1@hdr", "report_rdbflag@hdr", "entryno@body",
                    "vertco_type@body", "ppcode@conv_body", "datum_anflag@body", "datum_status@body",
                    "datum_event1@body", "datum_rdbflag@body"]
    # typical information that is only relevant for one sounding (repeated)
    sel = [ivar for ivar in variables if '@hdr' in ivar or '@modsurf' in ivar or '@surfbody' in ivar or '@conv' in ivar]
    sel = [ivar for ivar in sel if ivar not in default_drop]
    sel.remove('date@hdr')
    sel.remove('time@hdr')
    if False:
        station = fpv[sel].isel(hdrlen=ati)
        station = station.rename({'hdrlen': 'date'})
        station = station.assign_coords(date=alldates)
        station = station.rename({i: i.replace('@hdr', '') for i in list(station.data_vars)})
        station = station.rename({i: i.replace('@modsurf', '_msurf') for i in list(station.data_vars)})
        station = station.rename({i: i.replace('@surfbody_feedback', '_surf_fb') for i in list(station.data_vars)})
        station = station.rename({i: i.split('@conv')[0] for i in list(station.data_vars)})

        station.attrs.update(global_attrs)
        # sort by index
        station = station.sortby('date')

    try:
        os.mkdir(ident)
    except:
        pass

    if False:
        station.to_netcdf('./' + ident + '/ERA5_1_' + ident + '_station.nc')

    for l, var in enumerate(paras):
        tdx = time.time()
        dailyobs = numpy.zeros((1))
        if True:
            # refdate=date(1900,1,1)
            dailyobs = numpy.zeros(1)
            dailyfgd = numpy.zeros(1)
            # todo antime or time ? -> dailyhours is useless like this
            # todo this is not really useful?
            dailyobs, dailyfgd, dailybias, atn, dailyhours = filldaily(fpv['date@hdr'].values, fpv['antime'].values,
                                                                       fpv[u'vertco_reference_1@body'].values,
                                                                       fpv[u'varno@body'].values,
                                                                       fpv[u'obsvalue@body'].values,
                                                                       fpv[u'fg_depar@body'].values,
                                                                       fpv[u'biascorr@body'].values, ps, at, ati, var,
                                                                       numpy.nan)
            index = []
            for a in atn:
                # Year, month, day
                index.append(datetime(a // 10000, (a % 10000) // 100, a % 100))
                # index.append((date(a//10000,(a%10000)//100,a%100)-refdate).days)

            spl = xarray.Dataset()
            spl[names[l]] = (('hour', 'pres', 'date'), dailyobs)
            spl[names[l] + '_fg_dep'] = (('hour', 'pres', 'date'), dailyfgd)
            spl[names[l] + '_bias'] = (('hour', 'pres', 'date'), dailybias)
            spl['date'] = index
            spl['pres'] = ps
            spl['hour'] = numpy.array([0, 12])
            spl['pres'].attrs.update({'units': 'Pa', 'standard_name': 'air_pressure', 'axis': 'Z'})
            spl['date'].attrs.update({'axis': 'T'})
            spl['hour'].attrs.update({'units': 'h', 'standard_name': 'standard_launch_time'})
            if l == 0:
                spl['launch'] = (('hour', 'date'), dailyhours)
            spl = spl.reindex(date=alldates)
            if l == 0:
                launches = spl['launch'].copy()
                spl = spl.drop('launch')
        else:
            #
            # Interpolation to other levels
            # Problem with Speed, Time selection and indices
            # todo get it working
            dailyobs, dailyfgd, dailybias = filldaily_int(ati,
                                                          fpv[u'vertco_reference_1@body'].values,
                                                          fpv[u'varno@body'].values,
                                                          fpv[u'obsvalue@body'].values,
                                                          fpv[u'fg_depar@body'].values,
                                                          fpv[u'biascorr@body'].values,
                                                          ps,
                                                          var,
                                                          numpy.nan)

            spl = xarray.Dataset()
            spl[names[l]] = (('date', 'pres'), dailyobs)
            spl[names[l] + '_fg_dep'] = (('date', 'pres'), dailyfgd)
            spl[names[l] + '_bias'] = (('date', 'pres'), dailybias)
            if True:
                spl['date'] = alldates
                spl['pres'] = ps
                spl['pres'].attrs.update({'units': 'Pa', 'standard_name': 'air_pressure', 'axis': 'Z'})
                spl['date'].attrs.update({'axis': 'T'})

        for v in [names[l], names[l] + '_fg_dep', names[l] + '_bias']:
            spl[v].attrs.update(_metadata[v])

        spl['station_lat'] = fpv['lat@hdr'][-1]
        spl['station_lon'] = fpv['lon@hdr'][-1]
        spl['station_alt'] = fpv['stalt@hdr'][-1]
        spl['station_id'] = ident
        spl.attrs.update(global_attrs)
        #
        # IO
        #
        fno = os.path.expandvars('./' + ident + '/ERAI_' + ident + '_' + names[l] + '.nc')
        spl.to_netcdf(fno)
        td1 = time.time() - tdx
        print(ident, var, names[l], dailyobs.shape, td1)

    if True:
        station = fpv[sel].isel(hdrlen=ati)
        station = station.rename({'hdrlen': 'date'})
        station = station.assign_coords(date=alldates)
        station = station.rename({i: i.replace('@hdr', '') for i in list(station.data_vars)})
        station = station.rename({i: i.replace('@modsurf', '_msurf') for i in list(station.data_vars)})
        station = station.rename({i: i.replace('@surfbody_feedback', '_surf_fb') for i in list(station.data_vars)})
        station = station.rename({i: i.split('@conv')[0] for i in list(station.data_vars)})
        station['launch'] = launches  # (('hour','date'), dailyhours)
        station.attrs.update(global_attrs)
        # sort by index
        station = station.sortby('date')
        station.to_netcdf('./' + ident + '/ERAI_' + ident + '_station.nc')

    td2 = time.time() - td0
    td3 = time.time() - t
    print(ident, "DATE: ", td0 - t, "VAR:", td2, "TOTAL:", td3, "COUNT:", fpv[u'obsvalue@body'].shape)


def doquery(fno):
    t = time.time()
    #
    #
    #
    fn = fno[:-3]  # remove .gz
    ident = '0' + fn.split('.')[-2][1:]
    try:
        gz = gzip.open(fno)
        #
        # XArray (requires a netcdf3 - 64 bit version)
        #
        g = xarray.open_dataset(gz)
        #
        # Rearange
        #
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
    import glob
    import os

    if len(sys.argv) == 1:

        fns = '/raid8/srvx1/mblaschek/tmp/gzipped/era5.2402.conv.*.nc.gz'
        # fns = '/raid60/scratch/leo/scratch/era5/odbs/1/era5.conv.*.nc.gz'
        files = glob.glob(fns)

        for ifile in files[:]:
            if not os.path.isfile(ifile):
                files.remove(ifile)
                print(ifile)

        print("ERA-5 Stations:", len(files))
        print("Starting a Pool ...")
        p = Pool(20)
        res = p.map(doquery, files)
        print("Done")
    else:
        print(sys.argv)
        doquery(sys.argv[1])
