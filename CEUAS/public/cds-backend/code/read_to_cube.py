#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__all__ = ['read_ragged_array_to_cube']

import xarray as xr
import logging

logger = logging.getLogger(__name__)
# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(30)  # respond only to Warnings and above
# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s | %(funcName)s - %(levelname)s - %(message)s')
# add formatter to ch
ch.setFormatter(formatter)
# add ch to logger
logger.addHandler(ch)


def read_ragged_cdm_to_array(filename: str, odb_codes: bool = True, daynight: bool = False,
                             variables: list = None, read_feedback: bool = True,
                             optional_variables: dict = None,
                             **kwargs) -> xr.Dataset:
    """
    Read netCDF 4 HDF5 group data into a data cube on standard pressure levels
    this meant only for the CDS backend, not public

    duplicates are overwritten

    Args:
        filename: Filename of netcdf
        odb_codes: Search for ODB or CDM variable codes
        daynight: Add hour dimension with 00Z and 12Z
        variables: List of variables as in the ODB Codes
        read_feedback: Look for ERA5 Feedback information
        optional_variables: Dictionary with group and list of variables to read

    Notes:
        ODB Codes:
         85 : t
         38 : rh
         36 : td
         34 : dpd
         117 : z
         106 : dd
         107 : ff
         104 : u
         105 : v
         39 : q
        CDM Codes:
         1 : z
         2 : t
         3 : u
         4 : v
         111 : dd
         112 : ff
         59 : td
         299 : dpd
         29 : rh
         999 : p
    Returns:
        Dataset (plev, time)
    """
    import os
    import h5py
    import numpy as np
    import pandas as pd
    import xarray as xr

    if not os.path.isfile(filename):
        raise IOError('File can not be found', filename)

    if odb_codes:
        var_names = {'t': 85, 'rh': 38, 'td': 36, 'dpd': 34, 'z': 117, 'dd': 106, 'ff': 107, 'u': 104, 'v': 105,
                     'q': 39}
    else:
        var_names = {'z': 1, 't': 2, 'u': 3, 'v': 4, 'dd': 111, 'ff': 112, 'td': 59, 'dpd': 299, 'rh': 29, 'p': 999}

    var_names = dict({j: i for i, j in var_names.items()})
    std_plevs = np.asarray([10, 20, 30, 50, 70, 100, 150, 200, 250, 300, 400, 500, 700, 850, 925, 1000])
    dim = 'date_time'
    lev = 'z_coordinate'
    # Pressure Array -> indices
    ip = np.zeros(1001, dtype=np.int32)
    for i, j in enumerate(std_plevs):
        ip[j] = i

    data = {}
    with h5py.File(filename, 'r') as f:
        #
        # Header Group ?
        #
        if 'header_table' in f.keys():
            igroup = 'header_table'
            logger.info("Header Table found: %s", igroup)
            #
            # What to read from here?
            # and how ?
        #
        # Observation Table
        #
        var_index = {}
        if 'observations_table' in f.keys():
            igroup = 'observations_table'
            logger.info("Observations Table found: %s", igroup)
            # get date_time
            time = f[igroup][dim][:]  # should be seconds since 1900 01 01

            if 'units' not in f[igroup][dim].attrs.keys():
                logger.info("Datetime dimension %s has no units!!! Assuming seconds", dim)

            else:
                if 'seconds' not in f[igroup][dim].attrs['units'].decode():
                    logger.error("Datetime dimension %s not in seconds since ????, but %s", dim,
                                 f[igroup][dim].attrs['units'].decode())
                    raise RuntimeError("Datetime dimension %s not in seconds since ????, but %s" % (
                        dim, f[igroup][dim].attrs['units'].decode()))

            if 'units' in f[igroup][dim].attrs.keys():
                time_unit = f[igroup][dim].attrs['units'].decode()
            else:
                time_unit = 'seconds since 1900-01-01 00:00:00'

            # convert to datetime64 from seconds
            try:
                date = np.datetime64(" ".join(time_unit.split(' ')[-2:])) + time * np.timedelta64(1,
                                                                                                  time_unit.split(' ')[
                                                                                                      0][0])
            except:
                date = np.datetime64(" ".join(time_unit.split(' ')[-2:])) + time * np.timedelta64(1,
                                                                                                  time_unit.split(' ')[
                                                                                                      0][0].upper())
            # get pressures
            plev = f[igroup][lev][:]
            try:
                logger.info("%s units [%s]", dim, f[igroup][lev].attrs['units'].decode())
            except:
                pass
            iplev = np.in1d(plev, std_plevs * 100)  # only std pressure (can have NaN)
            plev = plev.astype(np.int32) // 100  # needs to be integer for indexing, hPa as well
            # get observed_variables
            ovariables = f[igroup]['observed_variable'][:]
            unique_variables = np.unique(ovariables)
            if variables is None:
                variables = unique_variables

            # get observation_value
            for i in unique_variables:
                if i not in var_names.keys() or i not in variables:
                    logger.info("Skipping.... %s not in %s ", i, var_names.keys())
                    continue
                # only this variable + only standard pressure levels
                ivar = (ovariables == i) & iplev

                # check dates and pressure levels for duplicates
                timeplev = np.stack((time[ivar], plev[ivar],), axis=0)
                u, c = np.unique(timeplev, axis=1, return_counts=True)
                if c[c > 1].size > 0:
                    logger.warning("%s %s Duplicates: %d overwritten", igroup, i, c[c > 1].size)
                    if kwargs.get('debug', False): raise RuntimeError("Duplicates")

                var_index[i] = ivar
                logger.info("Converting to Cube .... %s %s", i, ivar.sum())

                itime, iobs = table_to_cube(time[ivar], ip[plev[ivar]], f[igroup]['observation_value'][ivar])
                logger.info("Converting to Xarray ... %s", iobs.shape)

                data[var_names[i]] = xr.DataArray(iobs, coords=(std_plevs, date[ivar][itime]),
                                                  dims=('plev', 'time'), name=var_names[i])
            # end for i in unique_variables
        #
        # Feedback Group
        #
        if read_feedback and 'era5fb' in f.keys():
            igroup = 'era5fb'
            logger.info("ERA5 Feedback  found: %s", igroup)
            #
            # Analysis
            #
            if 'an_depar@body' in f[igroup].keys():
                logger.info('ERA5 Feedback: an_depar@body')
                for i, ivar in var_index.items():
                    itime, iobs = table_to_cube(time[ivar], ip[plev[ivar]], f[igroup]['an_depar@body'][ivar],
                                                )

                    data[var_names[i] + '_an_dep'] = xr.DataArray(iobs, coords=(std_plevs, date[ivar][itime]),
                                                                  dims=('plev', 'time'),
                                                                  name=var_names[i] + '_an_dep')

            #
            # First Guess
            #
            if 'fg_depar@body' in f[igroup].keys():
                logger.info('ERA5 Feedback: fg_depar@body')
                for i, ivar in var_index.items():
                    itime, iobs = table_to_cube(time[ivar], ip[plev[ivar]], f[igroup]['fg_depar@body'][ivar],
                                                )

                    data[var_names[i] + '_fg_dep'] = xr.DataArray(iobs, coords=(std_plevs, date[ivar][itime]),
                                                                  dims=('plev', 'time'),
                                                                  name=var_names[i] + '_fg_dep')

            #
            # Biascorr
            #
            if 'biascorr@body' in f[igroup].keys():
                logger.info('ERA5 Feedback: biascorr@body')
                for i, ivar in var_index.items():
                    itime, iobs = table_to_cube(time[ivar], ip[plev[ivar]], f[igroup]['biascorr@body'][ivar],
                                                )

                    data[var_names[i] + '_biascorr'] = xr.DataArray(iobs, coords=(std_plevs, date[ivar][itime]),
                                                                    dims=('plev', 'time'),
                                                                    name=var_names[i] + '_biascorr')
        #
        # Optional Variables (mostly in observations_table)
        #
        if len(optional_variables.keys()) > 0:
            for igroup in optional_variables.keys():
                if igroup not in f.keys():
                    logger.warning('Optional Group %s not in file. skipping', igroup)
                    continue

                if igroup != 'observations_table':
                    logger.info('not obs_table %s', igroup)
                    continue

                for jvar in optional_variables[igroup]:
                    if jvar not in f[igroup].keys():
                        logger.warning('Optional Variable %s not in group %s, file. skipping', jvar, igroup)
                        continue

                    logger.info('Reading: %s/%s', igroup, jvar)
                    if len(f[igroup][jvar].shape) > 1:
                        xobs = f[igroup][jvar][:, :]
                        for i, ivar in var_index.items():
                            iobs = xobs[ivar, :].astype(object).sum(1).astype(str)
                            itime, iobs = table_to_cube(time[ivar], ip[plev[ivar]], iobs)
                            data[var_names[i] + '_' + jvar] = xr.DataArray(iobs, coords=(std_plevs, date[ivar][itime]),
                                                                           dims=('plev', 'time'),
                                                                           name=var_names[i] + '_' + jvar)
                    else:
                        for i, ivar in var_index.items():
                            itime, iobs = table_to_cube(time[ivar], ip[plev[ivar]], f[igroup][jvar][ivar])
                            data[var_names[i] + '_' + jvar] = xr.DataArray(iobs, coords=(std_plevs, date[ivar][itime]),
                                                                       dims=('plev', 'time'),
                                                                       name=var_names[i] + '_' + jvar)
        #
        # end
        #
    logger.info("Converting to Dataset....")
    data = xr.Dataset(data)
    data['plev'].attrs.update({'units': 'hPa', 'standard_name': 'air_pressure'})

    if daynight:
        # convert datetime to standard 00 and 12 times and add dimensions hour
        # ctmp = rt.met.std.to_hours(, times=(0, 12))
        data = align_datetime(data, times=(0, 12), span=kwargs.get('span', 6), freq='12h')
        dim = 'time'
        data = data.sel(**{dim: (data[dim].dt.hour.isin((0, 12)) & (data[dim].dt.minute == 0))})
        data = dict(data.groupby(dim + '.hour'))
        for ikey in data.keys():
            data[ikey] = data[ikey].assign_coords(
                {dim: data[ikey][dim].to_index().to_period('D').to_timestamp().values})

        data = xr.concat(data.values(), dim=pd.Index(data.keys(), name='hour'))
        data['flag_stdtime'] = data['flag_stdtime'].fillna(0)
        # make sure the shape is as promissed:
        data = data.reindex({'hour': list((0, 12))})
    return data


#
# Align launch times to standard sounding times
#

def align_datetime(data, dim='time', plev='plev', times=(0, 12), span=6, freq='12h', **kwargs):
    """ Standardize datetime to times per date, try to fill gaps

    Args:
        data (DataArray, Dataset): Input data
        dim (str): datetime dimension
        plev (str): pressure level dimension
        times (tuple): sounding times
        span (int): plus minus times (smaller than freq/2)
        freq (str): frequency of output times

    Returns:
        xarray.DataArray : datetime standardized DataArray

    """
    import numpy as np
    from pandas import DatetimeIndex
    from xarray import DataArray, Dataset

    if not isinstance(data, (DataArray, Dataset)):
        raise ValueError('Requires a DataArray or Dataset', type(data))

    if dim not in data.dims:
        raise ValueError('Requires a datetime dimension', dim)

    if int(24 / (len(times) * 2)) < span:
        raise ValueError("Times and span do not agree!?", times, span)

    if int(24 / int(freq[:-1])) != len(times):
        raise ValueError("Times and freq do not match:", times, freq)

    if span > int(freq[:-1]) // 2:
        raise ValueError("Frequency and Span need to be consistent (span < freq/2): ", freq, span)

    dates = data[dim].values.copy()
    #
    # Count levels per date
    #
    _fix_datetime = np.vectorize(fix_datetime)
    newdates = _fix_datetime(dates, span=span)  # (time: 33%)
    resolution = np.zeros(newdates.size)
    #
    # check for duplicates in standard launch times
    #
    u, c = np.unique(newdates, return_counts=True)
    conflicts = u[c > 1]
    if conflicts.size > 0:
        counts = _count_data(data, dim=dim, plev=plev)
        logger.warning("Conflicts: %d in %d", conflicts.size, newdates.size)
        for i in conflicts:
            indices = np.where(newdates == i)[0]  # numbers  (time: 45%)
            #
            # Count available data (DataArray or Dataset)
            #
            # slow
            # counts = data.isel(**{dim: indices}).count(plev).values
            # counts = _count_data(data.isel(**{dim: indices}), dim=dim, plev=plev)
            # slow end
            icounts = counts[indices]
            #
            # offsets to standard launch time
            #
            offset = np.abs((dates[indices] - i) / np.timedelta64(1, 'h'))
            j = np.argsort(offset)  # sort time offsets (first we want)
            jmax = np.argmax(icounts[j])  # check if counts from other time is larger
            if jmax != 0:
                #
                # there is a sounding with more level data (+/- 1 hour)
                #
                if (offset[j][0] + 1) <= offset[j][jmax]:
                    # ok close enough
                    jj = j.copy()
                    jj[j == 0] = jmax  # first pos is now at the position of the maximum
                    jj[j == jmax] = 0  # maximum is now first
                    j = jj
                #
                # there is a sounding with + 2 more levels
                #
                elif (icounts[j][0] + 2) <= icounts[j][jmax]:
                    # a lot more
                    jj = j.copy()
                    jj[j == 0] = jmax  # first pos is now at the position of the maximum
                    jj[j == jmax] = 0  # maximum is now first
                    j = jj
                else:
                    pass  # keep time sorting

            for m, k in enumerate(offset[j]):
                if m == 0:
                    continue  # this is the minimum

                # change back the others or add a delay to remove duplicates
                if k == 0:
                    newdates[indices[j][m]] += np.timedelta64(1, 'h')  # add offset
                    resolution[indices[j][m]] = 1  # add hour
                else:
                    newdates[indices[j][m]] = dates[indices[j][m]]  # revert back
                    resolution[indices[j][m]] = -1  # revert back
    #
    # recheck for standard times
    #
    idx_std = DatetimeIndex(newdates).hour.isin(times)
    u, c = np.unique(newdates[idx_std], return_counts=True)  # check only standard times
    conflicts = u[c > 1]
    if conflicts.size > 0:
        logger.warning("Conflicts remain: %d  Std: %d New: %d", conflicts.size, idx_std.sum(), newdates.size)
        if kwargs.get('debug', False): raise RuntimeError("Duplicates")
    #
    # new dates / new object
    #
    data = data.assign_coords({dim: newdates})
    #
    # delay
    #
    nn = (resolution > 0).sum()
    nx = (~idx_std).sum()
    data['hours'] = (dim, DatetimeIndex(dates).hour.astype(int))  # new coordinate for delays
    data.attrs['std_times'] = str(times)
    data['hours'].attrs.update({'long_name': 'Launch time', 'units': 'h', 'times': str(times)})
    data['flag_stdtime'] = (dim, resolution.astype(int))
    data['flag_stdtime'].attrs.update({'units': '1', 'standard_name': 'flag_standard_time_conflict_resolution',
                                       'info': '0: preferred, -1: lesser candidate, 1: duplicate, less data'})

    logger.info("Modified: %d No Standard: %d of %d", nn, nx, newdates.size)

    if not all(data[dim].values == np.sort(data[dim].values)):
        logger.info("Sorting by %s", dim)
        data = data.sortby(dim)
    return data


def _count_data(data, dim='time', plev='plev'):
    from xarray import DataArray
    #
    # Count data per pressure level (if it is a dimension)
    #
    if plev in data.dims:
        #
        # Array format
        #
        if isinstance(data, DataArray):
            return data.count(plev).values
        else:
            return data.count(plev).to_dataframe().sum(axis=1).values  # sum across variables

    elif data[dim].to_index().is_unique:
        #
        # has not pressure levels
        #
        if isinstance(data, DataArray):
            return data.count(dim).values
        else:
            return data.count(dim).to_dataframe().sum(axis=1).values
    else:
        #
        # Table format
        #
        return data.groupby(dim).count().to_dataframe().max(axis=1).values


def fix_datetime(itime, span=6, debug=False):
    """ Fix datetime to standard datetime with hour precision

    Args:
        itime (datetime): Datetime
        span (int): allowed difference to standard datetime (0,6,12,18)

    Returns:
        datetime : standard datetime
    """
    import pandas as pd
    itime = pd.Timestamp(itime)  # (time: 34%)
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
                    rx = rx + pd.DateOffset(days=1)
                return rx.to_datetime64()
        else:
            if lower <= itime.hour < upper:
                rx = itime.replace(hour=ihour, minute=0, second=0, microsecond=0)
                if itime.hour >= (24 - span):
                    rx = rx + pd.DateOffset(days=1)
                return rx.to_datetime64()


def table_to_cube(time, plev, obs, return_indexes=False):
    """ Helper Function to rearange the table to cube

    Args:
        time (array): time in seconds
        plev (array): pressure level indices (not pressure itself)
        obs (array): observation values
        hours (bool): flag to enable hour
        return_indexes (bool): do nothing, return only indices (deprecated)

    Returns:
        array, array : time indices, cube
    """
    import numpy as np

    xtime, jtime, itime = np.unique(time, return_index=True, return_inverse=True)
    if return_indexes:
        return jtime, xtime.size, plev, itime
    data = np.full((16, xtime.size), np.nan, dtype=obs.dtype)
    data[plev, itime] = obs
    return jtime, data


def read_ragged_cdm(filename, odb_codes=True, **kwargs):
    import os
    import h5py
    import numpy as np
    import pandas as pd
    import xarray as xr

    if not os.path.isfile(filename):
        raise IOError('File can not be found', filename)

    if odb_codes:
        var_names = {'t': 85, 'rh': 38, 'td': 36, 'dpd': 34, 'z': 117, 'dd': 106, 'ff': 107, 'u': 104, 'v': 105,
                     'q': 39}
    else:
        var_names = {'z': 1, 't': 2, 'u': 3, 'v': 4, 'dd': 111, 'ff': 112, 'td': 59, 'dpd': 299, 'rh': 29, 'p': 999}
    var_names = dict({j: i for i, j in var_names.items()})
    # standard pressure levels
    std_plevs = [10, 20, 30, 50, 70, 100, 150, 200, 250, 300, 400, 500, 700, 850, 925, 1000]
    groups = ['observations_table', 'era5fb']
    # Variables
    variables = ['observed_variable', 'observation_value', 'an_depar@body', 'fg_depar@body',
                 'biascorr@body']  # list(f.keys())

    dim = 'date_time'
    lev = 'z_coordinate'

    data = {}
    with h5py.File(filename, 'r') as f:
        # Groups
        for igroup in groups:
            if igroup not in f.keys():
                logger.info("Group not found: %s", igroup)
                continue

            if dim in f[igroup].keys():
                # time
                time = f[igroup][dim][:]  # should be seconds since 1900 01 01

                if 'units' not in f[igroup][dim].attrs.keys():
                    logger.info("Datetime dimension %s has no units!!! Assuming seconds", dim)
                else:
                    if 'seconds' not in f[igroup][dim].attrs['units'].decode():
                        logger.error("Datetime dimension %s not in seconds since ????, but %s", dim,
                                     f[igroup][dim].attrs['units'].decode())
                        raise RuntimeError("Datetime dimension %s not in seconds since ????, but %s" % (dim,
                                                                                                        f[igroup][
                                                                                                            dim].attrs[
                                                                                                            'units'].decode()))
                time_unit = f[igroup][dim].attrs['units'].decode() if 'units' in f[igroup][
                    dim].attrs.keys() else 'seconds since 1900-01-01 00:00:00'
                # convert to datetime64 from seconds
                time = np.datetime64(" ".join(time_unit.split(' ')[-2:])) + time * np.timedelta64(1, 's')
                data['time'] = xr.DataArray(time)

            if lev in f[igroup].keys():
                # pressure to hPa
                plev = f[igroup][lev][:].astype(np.int32) // 100  # requires hPa
                data['plev'] = xr.DataArray(plev, attrs={'units': 'hPa', 'standard_name': 'air_pressure'})
            # Read data
            for ivar in variables:
                if ivar not in f[igroup].keys():
                    continue

                # Create Xarray DataArray
                data[ivar] = xr.DataArray(f[igroup][ivar][:])
                # Copy Attributes
                data[ivar].attrs.update(
                    {i: j.decode() for i, j in f[igroup][ivar].attrs.items() if isinstance(j, bytes)})
        # Copy Global and Coordinate Attributes
        # global_attributes = {i: j.decode() for i, j in f.attrs.items()}
        # coord_attributes = {lev: {i: j.decode() for i, j in f[lev].attrs.items() if isinstance(j, bytes)}}

    data = xr.Dataset(data)  # .sortby(dim)
    data = data.swap_dims({'dim_0': 'time'})
    # split into separate variables
    # convert to xarray
    logger.info("Converting to 2D ...")
    logger.info(str(np.unique(data.observed_variable)))
    data = dict(data.groupby(data.observed_variable))
    new = []
    # Loop and rename
    for i, j in data.items():
        if i not in var_names.keys():
            logger.warning("Skipping.... %s not in %s", i, var_names.keys())
            continue
        iname = var_names[i]
        j = j.drop_vars('observed_variable').rename({'observation_value': iname})
        j = j.rename({i: iname + "_" + i.replace('@body', '') for i in j.data_vars if i not in ['plev', iname]})
        # Make 2D
        j = j.to_dataframe()
        j.index.name = 'time'
        #
        # select only valid levels
        #
        j = j[j['plev'].isin(std_plevs)]
        #
        # convert to xarray
        #
        j = j.reset_index().set_index(['time', 'plev'])
        if not j.index.is_unique:
            logger.warning("%s Non-unique index, removing duplicates... %d", iname, j.index.duplicated().sum())
            logger.info(str(j.loc[j.index.duplicated(keep=False)]))
            j = j.loc[~j.index.duplicated()]  # remove duplicated

        j = j.to_xarray()  # 1D -> 2D
        new.append(j)
    data = xr.merge(new)
    return data


def read_ragged_array_to_cube(filename, dim='time', lev='plev', variables=None, std_plevs=None, **kwargs):
    """
    The function reads a CF compliant upper air data file from the CDS backend and transforms the data into
    a three dimensional data cube with dimensions hour (2), plev (16), time (days).
    This representation is convenient for statistical breakpoint analysis.

    Args:
        filename (str): Input NetCDF file
        dim (str): Name of datetime dimension in file
        lev (str): Name of pressure dimension in file
        variables (list): List of variables to read
        std_plevs (list): List of pressure levels to use, default 16

    Returns:
        xarray.Dataset : variables and coordinates and metadata

    Note:
        At the moment this function is limited to 00 and 12 UTZ, intermediate times might overlap
        only valid for files that have std pressure levels
    """
    import os
    import h5py
    import numpy as np
    import pandas as pd
    import xarray as xr

    if not os.path.isfile(filename):
        raise IOError('File can not be found', filename)

    # standard pressure levels
    if std_plevs is None:
        std_plevs = [10, 20, 30, 50, 70, 100, 150, 200, 250, 300, 400, 500, 700, 850, 925, 1000]
    else:
        logger.warning("Standard pressure levels need to be in hPa")

    # Pressure Array -> indices
    ip = np.zeros(1001, dtype=np.int32)
    for i, j in enumerate(std_plevs):
        ip[j] = i

    data = {}
    with h5py.File(filename, 'r') as f:
        # time
        time = f[dim][:]  # should be seconds since 1900 01 01
        if 'units' not in f[dim].attrs.keys():
            logger.warning("Datetime dimension", dim, "has no units!!! Assuming seconds")
        else:
            if 'seconds' not in f[dim].attrs['units'].decode():
                logger.error("Datetime dimension", dim, "not in seconds since ????",
                             f[dim].attrs['units'].decode())
                raise RuntimeError("Datetime dimension", dim, "not in seconds since ????",
                                   f[dim].attrs['units'].decode())
        # pressure
        # todo causes a problem with non standard pressure levels see read_ragged_cdm_to_array fix
        plev = f[lev][:].astype(np.int32) // 100  # requires hPa
        if not np.all(np.unique(plev) == plev):
            logger.warning("non standard pressure levels problem")
        # Calculate Indices for cube positions
        date_index = np.array(time // 86400, dtype=np.int32)
        date_index -= date_index[0]  # index to 0
        hour_index = np.array(((time + 21600) % 86400) // 43200, dtype=np.int32)
        plev_index = ip[plev]
        # Convert time unit to datetime
        time_unit = f['time'].attrs['units'].decode() if 'units' in f[
            dim].attrs.keys() else 'seconds since 1900-01-01 00:00:00'
        time = pd.to_datetime(time // 86400, unit='D', origin=" ".join(time_unit.split(' ')[-2:]))
        new = np.full(date_index.max() + 1, '1850-01-01',
                      dtype='datetime64[ns]')  # missing date is 1850 -> remove later
        new[date_index] = time.values  # fill in dates
        time = new
        # Variables
        if variables is None:
            variables = list(f.keys())
        # remove dimension and trajectory variables
        variables = [i for i in variables if
                     i not in [dim, lev, 'trajectory', 'trajectory_index', 'trajectory_label'] and 'string' not in i]
        # Read data
        for ivar in variables:
            if ivar not in f.keys():
                continue
            # Create Xarray DataArray
            data[ivar] = xr.DataArray(np.full((2, 16, time.size), np.nan, dtype=np.float32),
                                      coords=([0, 12], std_plevs, time), dims=('hour', lev, dim), name=ivar)
            # Use Indices to fill in data
            data[ivar].values[hour_index, plev_index, date_index] = f[ivar][:]
            # Copy Attributes
            data[ivar].attrs.update({i: j.decode() for i, j in f[ivar].attrs.items() if isinstance(j, bytes)})
        # Copy Global and Coordinate Attributes
        global_attributes = {i: j.decode() for i, j in f.attrs.items()}
        coord_attributes = {lev: {i: j.decode() for i, j in f[lev].attrs.items() if isinstance(j, bytes)}}

    data = xr.Dataset(data).sortby(dim)
    # remove redundant information
    if 'lat' in data.data_vars:
        data['lat'] = data['lat'].reduce(np.nanmean, ('hour', lev))
    if 'lon' in data.data_vars:
        data['lon'] = data['lon'].reduce(np.nanmean, ('hour', lev))
    # remove missing times
    data = data.sel({dim: slice('1900', None)})
    # set Attributes
    data.attrs.update(global_attributes)
    for i, j in coord_attributes.items():
        data[i].attrs.update(j)
    return data


def read_index(ovars={}, offset=0, ip=None, sname=''):
    '''
    read_index reads a CF compliant upper air data file from the CDS backend such that
    its content can be easily transformed into a three dimensional data cube with dimensions
    hour (2), pressure (16), day (days since 19000101 + offset).
    This is quite convenient for statistical breakpoint analysis.

    The function expects
    - a dictionary ovars with variables to be read, e.g.
    ovars={'obs':'ta','obs_minus_bg'='obs_minus_bg','obs_minus_an'='obs_minus_an','bias_estimate'='bias_estimate'}.
    The values of this dictionary should be valid variable names of the file read.
    - an offset from 19000101 in days. This can be useful if one wants to produce a data cube which starts e.g. in 1950 instead
    of 1900, to save memory.
    - the name of the netCDF file to be read. The netCDF file contains CF compliant ragged arrays,
    but only on the above stated pressure levels. One can ensure with the backend request that only those levels are selected.
    A working example is:
    'curl -H "Content-Type: application/json" -X POST --digest --data '{"statid":"01001","date":[19390101,20191231],
    "pressure_level":[1000,2000,3000,5000,7000,10000,15000,20000,30000,40000,50000,70000,85000,92500,100000],
    "variable":"temperature"}'
    -o download.zip http://early-upper-air.copernicus-climate.eu'

    The function returns a dictionary d, which contains:
    sname = WMO number,
    lat = latitude,
    lon = longitude,
    hindex = index of hour dimension (0 or 1)
    pindex = index of pressure dimension (0 for 10 hPa, 1 for 20 hPa, ... 13 for 850 hPa, 14 for 925 hPa, 15 for 1000 hPa
    dindex = index of day dimension (days since 19000101 + offset in days
    datetime = the date and time expressed as seconds since 19000101
    optional output variables, as requested from the ovars dictionary.
    ovars={'obs':'ta'}, then the output dictionary will contain
    d['obs]=f['ta'], i.e. all air temperature values.

    As such the returned dictionary contains the most relevant content of the file in a quite compact format

    A variable can then be read into a data cube as:
    cube(d[hindex],d[pindex],d[dindex])=d[ovars[v]],
    which is very efficient.

    Leo Haimberger, 19 February 2020
    '''
    import time
    import h5py
    import numpy

    if sname == '' or ip is None:
        return dict()
    tt = time.time()
    with h5py.File(sname) as f:
        ft = f['time'][:]  # time must be in seconds since 19000101
        fp = numpy.array(f['plev'], dtype=numpy.int32) // 100  # f['plev'] must be in hPa!

        dindex = numpy.array(ft // 86400 - offset, dtype=numpy.int32)  # offset must be divisable by 86400!
        hindex = numpy.array(((ft + 21600) % 86400) // 43200, dtype=numpy.int32)
        pindex = ip[fp]

        opt = {}
        for k, v in ovars.items():
            try:
                opt[k] = f[v][:]
            except:
                print(sname, v, 'not available')

        d = dict(sname=sname.split('_')[1], lat=f['lat'][0], lon=f['lon'][0], hindex=hindex, pindex=pindex,
                 dindex=dindex,
                 datetime=ft, **opt)
        print(time.time() - tt)
    return d


def logger_level(level: int):
    """ logging level setter

    Args:
        level: 0, 10, 20, 30, 40, 50

    """
    global logger
    logger.setLevel(level)
    logger.handlers[0].setLevel(level)


if __name__ == '__main__':

    """
    LEO Code to run read_index on a file
    
    import os
    import h5py  # as h5py
    import numpy
    import xarray
    import pandas as pd
    import copy
    import time
    import glob
    from multiprocessing import Pool
    from functools import partial
    import psutil
    import matplotlib.pylab as plt

    os.chdir(os.path.expandvars('$RSCRATCH/era5/odbs/1/test'))
    os.chdir(os.path.expandvars('$HOME'))

    slist = glob.glob('dest_?????_air_temperature*')

    indexmax = 45000
    ip = numpy.zeros(1001, dtype=numpy.int32)
    ip[10], ip[20], ip[30], ip[50], ip[70], ip[100], ip[150], ip[200], ip[250], ip[300], ip[400], ip[500], ip[700], ip[
        850], \
    ip[925], ip[1000] = numpy.arange(16)

    talist = []
    bglist = []
    sdict = {}
    tt = time.time()
    tcube = numpy.empty((2, 16, indexmax), dtype=numpy.float32)
    tcube.fill(numpy.nan)
    bgcube = numpy.empty((2, 16, indexmax), dtype=numpy.float32)

    p = Pool(25)
    func = partial(read_index, {'obs': 'ta', 'obs_minus_bg': 'obs_minus_bg'}, 0, ip)
    sdicts = list(map(func, slist[:100]))

    # remove empty dicts
    for k in range(len(sdicts) - 1, -1, -1):
        if not sdicts[k]:
            del sdicts[k]

    # fill data cube with dictionary from first valid file
    tcube[sdicts[0]['hindex'], sdicts[0]['pindex'], sdicts[0]['dindex']] = sdicts[0]['obs']

    process = psutil.Process(os.getpid())
    print(process.memory_info().rss / 1024 / 1024)
    print(time.time() - tt)
    plt.plot(numpy.arange(45000) / 365.25, tcube[1, 5, :] - tcube[0, 5, :])
    """

    # 53845
    ifile = '/raid60/scratch/federico/20MARCH2020_SmallStations/0-20001-0-53845_CEUAS_merged_v0.nc'
    # newer version (4.2020)
    ifile = '/raid60/scratch/federico/0-20000-0-53845_CEUAS_merged_v0.nc'

    # 41883
    ifile = '/raid60/scratch/federico/20MARCH2020_SmallStations/0-20000-0-41883_CEUAS_merged_v0.nc'

    # 72764
    ifile = "/raid60/scratch/leo/scratch/era5/odbs/2/era5_2/" + "0-20000-0-STATIONID_era5_2_harvested_era5.conv._1:72764.nc"  # 1950 to 1978

    # Read with an older version of the function
    # converts to DataFrame and uses Xarray
    # either it has odb_codes or cdm_codes
    if False:
        xx = read_ragged_cdm(ifile, odb_codes=True)
        print(xx)

    # Read with an modern version
    #
    xx = read_ragged_cdm_to_array(ifile, odb_codes=True)
    print(xx)

    jfile = "/raid60/scratch/leo/scratch/era5/odbs/1/era5_1/" + '0-20000-0-72764_era5_1_harvested_era5.conv.72764.txt.gz.nc'  # 1979 to present
    xy = read_ragged_cdm_to_array(jfile, odb_codes=True)
    print(xy)
