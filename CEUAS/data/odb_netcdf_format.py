#!/usr/bin/env python
# -*- coding: utf-8 -*-

__doc__ = """

formerly known as numpyquery.py -> LH

modified by MB
Modifications:
    + names of variables
    + attributes, metadata
    + some functions
    + interpolation (fixed)
    + datetime conversion and alignment into hours
    
Example:
    python odb_netcdf_format.py [odb file]
    
    File can be .nc, .ascii or .nc.gz

Works on SRVX 8

Last used: Don Okt 24 12:15:33 CEST 2019

"""
import gzip
import os.path
import time
from datetime import datetime
from multiprocessing import Pool

import numpy
import pandas
import xarray
from numba import __version__ as numbaversion

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

eraplevs = [1000., 2000., 3000., 5000., 7000., 10000., 12500., 15000., 17500., 20000., 22500., 25000., 30000., 35000.,
            40000., 45000., 50000., 55000., 60000., 65000., 70000., 75000., 77500., 80000., 82500., 85000., 87500.,
            90000., 92500., 95000., 97500., 100000.]
stdplevs = [1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 85000, 92500, 100000]


def fix_datetime(itime, span=6, debug=False):
    """ Fix datetime to standard datetime with hour precision

    Args:
        itime (datetime): Datetime
        span (int): allowed difference to standard datetime (0,6,12,18)

    Returns:
        datetime : standard datetime
    """
    itime = pandas.Timestamp(itime)
    # span=6 -> 0, 12
    # [18, 6[ , [6, 18[
    # span=3 -> 0, 6, 12, 18
    # [21, 3[, [3,9[, [9,15[, [15,21[
    for ihour in range(0, 24, span * 2):
        # 0 - 6 + 24 = 18
        lower = (ihour - span + 24) % 24
        # 0 + 6 + 24 = 6
        upper = (ihour + span + 24) % 24
        # 18 >= 18 or 18 < 6  > 00
        # 0 >= 18 or 0 < 6    > 00
        if debug:
            print("%d [%d] %d >= %d < %d" % (ihour, span, lower, itime.hour, upper))

        if (ihour - span) < 0:
            if itime.hour >= lower or itime.hour < upper:
                rx = itime.replace(hour=ihour, minute=0, second=0, microsecond=0)
                if itime.hour >= (24 - span):
                    rx = rx + pandas.DateOffset(days=1)
                return rx.to_datetime64()
        else:
            if lower <= itime.hour < upper:
                rx = itime.replace(hour=ihour, minute=0, second=0, microsecond=0)
                if itime.hour >= (24 - span):
                    rx = rx + pandas.DateOffset(days=1)
                return rx.to_datetime64()


def dataframe_to_array(data, dim='time', plev='plev', levels=None):
    if levels is None:
        levels = [1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 85000, 92500,
                  100000]

    # copy attributes
    attrs = data.attrs.copy()
    tatts = data[dim].attrs
    vatt = {i: data[i].attrs.copy() for i in data.data_vars}
    # dimensions for output
    varis = [dim, plev]
    # to pandas dataframe
    data = data.to_dataframe()
    # select only valid levels
    data = data[data[plev].isin(levels)]
    # convert to xarray
    data = data.reset_index().set_index(varis).to_xarray()  # 1D -> 2D
    # add attributes again
    for i, j in vatt.items():
        data[i].attrs.update(j)
    data.attrs.update(attrs)
    data[dim].attrs.update(tatts)
    return data


def dataframe(data, level_column, levels=None, variables=None, min_levels=3, keep_old_levels=False, **kwargs):
    """ Interpolate a database DataFrame according to pressure levels in level_column

    Interpolate:
    1. Select only levels with enough (min_levels) t and r values
    2. Interpolate each profile (date) vertically to levels

    Interpolation is only done at dates with enough Data

    Args:
        data (DataFrame):  Database with columns of non-uniform pressure levels
        level_column (str):  Database column with pressure levels
        levels (list, ndarray):  new pressure levels for interpolation
        variables (list): Variables to interpolate
        min_levels (int): minimum required levels per profile for interpolation
        keep_old_levels (bool) : keep old levels in database ?
        verbose (int): verbosness

    Returns:
    DataFrame : interpolated DataFrame with new pressure levels
    """
    import pandas as pd

    if not isinstance(data, pd.DataFrame):
        raise ValueError()

    if levels is None:
        levels = [1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 85000, 92500,
                  100000]

    data.index.name = 'time'

    if variables is not None:
        variables = list(set(variables + [level_column]))  # add p
        data = data.loc[:, data.columns.isin(variables)]
    # only use numeric columns
    data = data.select_dtypes(include=['number'])
    # Is there anything to work with?
    if len(data.columns.tolist()) < 2:
        raise ValueError("Requires at least 2 columns(%s,+) %s" % (level_column, ",".join(variables)))
    # Interpolate
    data = data.groupby(data.index).apply(table, level_column, levels, min_levels=min_levels,
                                          keep=keep_old_levels)
    if data.empty:
        raise RuntimeError("Not enough data for interpolation")
    # Change multi-index
    data = data.reset_index().drop('level_1', axis=1).sort_values(by=['time', level_column]).set_index('time',
                                                                                                       drop=True)
    return data


def table(data, level_column, levels, min_levels=3, keep=False):
    """ Wrapper Function for _np_profile to handle a DataFrame

    Args:
        data (DataFrame): Input DataFrame for a timestamp
        level_column (str): pressure level column
        levels (ndarray or list): new pressure levels
        min_levels (int): minimum required pressure levels
        keep (bool): keep old pressure levels

    Returns:
    DataFrame : new DataFrame with size of levels
    """
    import numpy as np
    import pandas as pd

    # dataset = dataset.iloc[np.unique(dataset[level_column], return_index=True)[1], :]   # subset
    data = data.sort_values(level_column)  # no subset
    pin = data[level_column].values
    data.drop(level_column, 1, inplace=True)
    if (data.count() > min_levels).sum() == 0:
        data = data.reset_index(drop=True)
        data[level_column] = pin
        if keep:
            data['orig'] = 0
            return data
        else:
            return data.iloc[np.in1d(pin, levels)]

    names = data.columns.tolist()
    if keep:
        alllevels = np.unique(np.sort(np.concatenate([pin, levels])))
        # 0 RAW, 1 Both, 2 NEW
        # orig = np.where(np.in1d(alllevels, levels), 2, 1) - np.where(np.in1d(alllevels, pin), 1, 0)
        orig = np.where(np.in1d(alllevels, pin), 0, 1)  # RAW=0, INT=1
        levels = alllevels

    data = np.apply_along_axis(profile, 0, data.values, pin, levels)
    data = pd.DataFrame(data, columns=names)
    data[level_column] = levels
    if keep:
        data['orig'] = orig
    return data


def profile(data, plevs, new_plevs):
    """ Modified numpy.interp Function for filtering nan

    Args:
        data (numpy.ndarray): Input dataset
        plevs (numpy.ndarray): Input pressure levels
        new_plevs (numpy.ndarray): Output pressure levels

    Returns:
        numpy.ndarray :  size of new_plevs
    """
    import numpy as np
    data = np.squeeze(data)  # remove 1-dims
    ix = np.isfinite(data)  # only finite values
    s = ix.sum()  # enough dataset left ?
    if s > 0:
        plevs, data = np.unique([plevs[ix], data[ix]], axis=1)
        data = np.interp(np.log(new_plevs), np.log(plevs), data, left=np.nan, right=np.nan)
        return data
    return np.full_like(new_plevs, np.nan)  # Nothing to do, but keep shape


def process_odb_obj(spl, fn, iname, ident, gattrs, t):
    ps = numpy.array(stdplevs)
    out = {}
    # gattrs = spl.attrs.copy()
    #
    # remove some variables incompatible with this transform /
    # unnecessary (could be usefull in the future)
    #
    for ivar in spl.data_vars:
        if ivar in [iname + '_orig', 'plev']:
            continue

        out[ivar] = dataframe_to_array(spl[[ivar, 'plev']], levels=ps)

    spl = xarray.merge(out.values())
    spl.attrs.update(gattrs)
    #
    # fix duplicates ?
    #
    _fix_datetime = numpy.vectorize(fix_datetime)
    newdates = _fix_datetime(spl.time.values, span=3)  # 3 > 6 > [0,6,12,18] UTC
    u, c = numpy.unique(newdates, return_counts=True)
    conflicts = u[c > 1]
    # todo print conflicts and resolution
    if conflicts.size > 0:
        for i in conflicts:
            indices = numpy.where(newdates == i)[0]
            a = (spl.time.values[indices[0]] - i) / numpy.timedelta64(1, 'h')
            b = (spl.time.values[indices[1]] - i) / numpy.timedelta64(1, 'h')
            if a < b:
                # change b back
                newdates[indices[1]] = spl.time.values[indices[1]]
            elif a > b:
                # change a back
                newdates[indices[0]] = spl.time.values[indices[0]]
            else:
                # a == b
                # really duplicated (check count)
                a = int(spl[iname].isel(time=indices[0]).count())
                b = int(spl[iname].isel(time=indices[1]).count())
                if a > b:
                    newdates[indices[1]] = newdates[indices[1]] + numpy.timedelta64(1, 'h')  # reset
                else:
                    newdates[indices[0]] = newdates[indices[0]] + numpy.timedelta64(1, 'h')
    td1 = time.time() - t
    print(ident, iname, "... Standard Time ...", numpy.size(conflicts), td1)
    #
    #
    #
    spl[iname + '_delay'] = ('time', ((spl.time.values - newdates) / numpy.timedelta64(1, 'h')).astype(int))
    spl = spl.assign_coords(time=newdates)
    spl[iname + '_delay'].attrs.update({'conflicts': conflicts.size})
    spl[iname + '_delay'].attrs['times'] = '0,6,12,18'
    #
    # Check for duplicates
    #
    u, c = numpy.unique(spl.time.values, return_counts=True)
    conflicts = u[c > 1]
    if conflicts.size > 0:
        for i in conflicts:
            indices = numpy.where(spl.time.values == i)[0]
            a = int(spl[iname].isel(time=indices[0]).count())
            b = int(spl[iname].isel(time=indices[1]).count())
            if a > b:
                spl.time.values[indices[1]] += numpy.timedelta64(1, 'h')  # reset
            else:
                spl.time.values[indices[0]] += numpy.timedelta64(1, 'h')
    td1 = time.time() - t
    print(ident, iname, "... Duplicates ...", numpy.size(conflicts), td1)
    u, c = numpy.unique(spl.time.values, return_counts=True)
    conflicts = u[c > 1]
    if conflicts.size > 0:
        for i in conflicts:
            indices = numpy.where(spl.time.values == i)[0]
            a = int(spl[iname].isel(time=indices[0]).count())
            b = int(spl[iname].isel(time=indices[1]).count())
            if a > b:
                spl.time.values[indices[1]] += numpy.timedelta64(1, 'h')  # reset
            else:
                spl.time.values[indices[0]] += numpy.timedelta64(1, 'h')
    td1 = time.time() - t
    print(ident, iname, "... Duplicates ...", numpy.size(conflicts), td1)

    #
    # convert to hour x plev x time Array
    #
    spl = dict(spl.sel(time=spl.time.dt.hour.isin([0, 6, 12, 18])).groupby('time.hour'))
    for ikey in spl.keys():
        spl[ikey] = spl[ikey].assign_coords(
            **{'time': spl[ikey]['time'].to_index().to_period('D').to_timestamp().values})

    spl = xarray.concat(spl.values(), dim=pandas.Index(spl.keys(), name='hour'))
    # make sure the shape is as promissed:
    spl = spl.reindex({'hour': [0, 6, 12, 18]})
    #
    # IO
    #
    if 'era5' in fn:
        fno = os.path.expandvars('./' + ident + '/ERA5_' + ident + '_' + iname + '.nc')
    else:
        fno = os.path.expandvars('./' + ident + '/ERAI_' + ident + '_' + iname + '.nc')

    encoding = {i: {'compression': 'gzip', 'compression_opts': 9} for i in list(spl.data_vars)}
    spl.to_netcdf(fno, mode='w', engine='h5netcdf', format='netcdf4', encoding=encoding)


def process_obj(fpv, fn, interpolate=True, convertarray=True):
    t = time.time()
    ps = numpy.array(stdplevs)

    ident = '0' + fn.split('.')[-2][1:]

    global_attrs = {'Conventions': 'CF-1.1',
                    'title': 'station daily temperature series',
                    'institution': 'University of Vienna',
                    'history': '19/03/20',
                    'source': 'radiosonde, ERA-5, ERA-Interim, ERA-40, RAOBCORE',
                    'references': 'www.univie.ac.at/theoret-met/research/raobcore',
                    'levels': 'plevs [%d -%d] #%d' % (min(ps), max(ps), len(ps)),
                    'libs': "NP(%s) XR(%s) NB(%s)" % (numpy.__version__, xarray.__version__, numbaversion)}

    td0 = time.time()
    #
    # too few dataset
    #
    if fpv['date@hdr'].size < 100:
        raise RuntimeError("Too few dates (<100)")

    global_attrs.update({'station_id': ident,
                         'station_lat': "%.2f N" % fpv['lat@hdr'].values[-1],
                         'station_lon': "%.2f E" % fpv['lon@hdr'].values[-1],
                         'station_alt': "%.1f m" % fpv['stalt@hdr'].values[-1]
                         })

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
                    "datum_event1@body", "datum_rdbflag@body", "datum_status@surfbody_feedback",
                    "datum_sfc_event@surfbody_feedback"]
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
        selection = [u'date@hdr', u'time@hdr', u'vertco_reference_1@body', u'obsvalue@body', u'fg_depar@body',
                     u'an_depar@body', u'biascorr@body']
        # only for these variables there are error statistics from the model
        if var in [2, 29]:
            if 'obs_error@errstat' in fpv.data_vars:
                selection += ['obs_error@errstat']
            if 'fg_error@errstat' in fpv.data_vars:
                selection += ['fg_error@errstat']
            if 'final_obs_error@errstat' in fpv.data_vars:
                selection += ['final_obs_error@errstat']
        #
        # Selecting variable
        #
        obs = fpv[selection].isel(hdrlen=fpv[u'varno@body'].values == var)
        if obs.hdrlen.size == 0:
            print(ident, var, names[l], "no data")
            continue

        dates = obs['date@hdr'].values.astype(int)
        times = obs['time@hdr'].values.astype(int)
        index = []
        for i in range(dates.shape[0]):
            # Datetime as from variable -> will be aligned later
            index.append(
                datetime(dates[i] // 10000, (dates[i] % 10000) // 100, dates[i] % 100, times[i] // 10000))

        #
        # Xarray
        #
        spl = xarray.Dataset()
        for ivar in obs.data_vars:
            if ivar in [u'date@hdr', u'time@hdr']:
                continue
            if ivar == u'vertco_reference_1@body':
                spl['plev'] = ('time', obs[ivar].values)
            if ivar == u'obsvalue@body':
                spl[names[l]] = ('time', obs[ivar].values)
            if ivar == u'fg_depar@body':
                spl[names[l] + '_fg_dep'] = ('time', obs[ivar].values)
            if ivar == u'an_depar@body':
                spl[names[l] + '_an_dep'] = ('time', obs[ivar].values)
            if ivar == u'biascorr@body':
                spl[names[l] + '_bias'] = ('time', obs[ivar].values)
            if ivar == 'obs_error@errstat':
                spl[names[l] + '_preerr'] = ('time', obs[ivar].values)
            if ivar == 'fg_error@errstat':
                spl[names[l] + '_fgerr'] = ('time', obs[ivar].values)
            if ivar == 'final_obs_error@errstat':
                spl[names[l] + '_finerr'] = ('time', obs[ivar].values)

        spl['time'] = index
        #
        # interpolation to standard pressure levels
        #
        if interpolate:
            td1 = time.time() - tdx
            print(ident, var, names[l], "... Interpolation ... <", numpy.size(ps),"> #", len(index), td1)
            try:
                data = dataframe(spl.to_dataframe(), 'plev', levels=ps, keep_old_levels=True)
            except RuntimeError:
                print(ident, var, names[l], "... Interpolation failed ... #", spl[names[l]].count().values, td1)
                continue

            spl = data.to_xarray()
            spl = spl.rename({'orig': names[l] + '_orig'})
            spl[names[l] + '_orig'].attrs.update({'standard_name': 'interpolation_flag',
                                                  'interpretation': '0: raw, 1: int'})

        spl['time'].attrs.update({'axis': 'T'})
        #
        # Metadata
        #
        for v in [names[l], names[l] + '_fg_dep', names[l] + '_an_dep', names[l] + '_bias']:
            spl[v].attrs.update(_metadata[v])
        #
        # Attributes
        #
        spl.attrs.update(global_attrs)
        # spl.attrs.update({'station_id': ident,
        #                   'station_lat': "%.2f N" % fpv['lat@hdr'].values[-1],
        #                   'station_lon': "%.2f E" % fpv['lon@hdr'].values[-1],
        #                   'station_alt': "%.1f m" % fpv['stalt@hdr'].values[-1]
        #                   })
        if 'hdrlen' in spl.coords:
            spl = spl.drop('hdrlen')
        #
        # IO
        #
        if '2402' in fn:
            fno = os.path.expandvars('./' + ident + '/odb_16plev_erai_' + ident + '_' + names[l] + '.nc')
        else:
            fno = os.path.expandvars('./' + ident + '/odb_16plev_era5_' + ident + '_' + names[l] + '.nc')

        encoding = {i: {'compression': 'gzip', 'compression_opts': 9} for i in list(spl.data_vars)}
        spl.to_netcdf(fno, engine='h5netcdf', format='netcdf4', encoding=encoding)
        #
        # convert Table to 2D (time x plev) Array
        #
        if convertarray:
            process_odb_obj(spl, fno, names[l], ident, global_attrs, tdx)
        #
        #
        #
        td1 = time.time() - tdx
        print(ident, var, names[l], len(index), td1)

    station = fpv[sel]
    dates = fpv['date@hdr'] * 1e6 + fpv['time@hdr']
    _, index = numpy.unique(dates, return_index=True)
    station = station.isel(hdrlen=index)
    alldates = []
    for a in dates[index]:
        b = int(a % 1e6)
        a = int(a // 1e6)
        alldates.append(datetime(a // 10000, (a % 10000) // 100, a % 100, b // 10000))

    alldates = numpy.array(alldates)
    station = station.assign_coords(hdrlen=alldates)
    station = station.rename({'hdrlen': 'time'})
    station = station.rename({i: i.replace('@hdr', '') for i in list(station.data_vars)})
    station = station.rename({i: i.replace('@modsurf', '_msurf') for i in list(station.data_vars)})
    station = station.rename({i: i.replace('@surfbody_feedback', '_surf_fb') for i in list(station.data_vars)})
    station = station.rename({i: i.split('@conv')[0] for i in list(station.data_vars)})
    if 'obstype' in station.data_vars:
        station['obstype'] = station['obstype'].fillna(-1).astype(int)
        station['obstype'].encoding['_FillValue'] = -1

    if 'codetype' in station.data_vars:
        station['codetype'] = station['codetype'].fillna(-1).astype(int)
        station['codetype'].encoding['_FillValue'] = -1

    if 'delay' in station.data_vars:
        station['delay'] = station['delay'].fillna(-1).astype(int)
        station['delay'].encoding['_FillValue'] = -1

    if 'sonde_type' in station.data_vars:
        station['sonde_type'] = station['sonde_type'].fillna(-1).astype(int)
        station['sonde_type'].encoding['_FillValue'] = -1

    if 'sensor' in station.data_vars:
        station['sensor'] = station['sensor'].fillna(-1).astype(int)
        station['sensor'].encoding['_FillValue'] = -1

    station.attrs.update(global_attrs)
    #
    # Index is not the same as variable netcdfs !!!
    # sort by index
    #
    station = station.sortby('time')
    if '2402' in fn:
        station.to_netcdf('./' + ident + '/' + ident + '_ERAI_station.nc', mode='w')
    else:
        station.to_netcdf('./' + ident + '/' + ident + '_ERA5_station.nc', mode='w')

    td2 = time.time() - td0
    td3 = time.time() - t
    print(ident, "DATE: ", td0 - t, "VAR:", td2, "TOTAL:", td3, "COUNT:", fpv[u'obsvalue@body'].shape)


def doquery(fno, debug=False, interpolate=True, convertarray=True):
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
            g = xarray.open_dataset(gz)   # ? engine='h5netcdf'  # ? faster
        else:
            g = xarray.open_dataset(fno)
        #
        # Rearange
        #
        print(ident, time.time() - t, "Start processing")
        if 'odb_16plev' in fn:
            pieces = fn.replace('.nc','').split('_')
            iname = pieces[-1]   # variable
            gattrs = g.attrs.copy()
            ident = pieces[-2]
            process_odb_obj(g, fn, iname, ident, gattrs, time.time())
            print(ident, iname, len(g.time), time.time() - t)
        else:
            process_obj(g, fn, interpolate=interpolate, convertarray=convertarray)
        #
        #
        #
    except Exception as e:
        print(ident, "Error", repr(e))
        if debug:
            raise e

    else:
        g.close()
        #
        #
        #
    print(ident, time.time() - t)


if __name__ == '__main__':
    import sys
    import os

    debug = True
    interpolate = True
    convertarray = True

    if len(sys.argv) == 2:
        print(sys.argv)
        doquery(sys.argv[1], debug=debug)

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
        p = Pool(5)
        res = p.map(doquery, files)
        print("Done")

    else:
        print("Usage: %s [files]" % sys.argv[0])
