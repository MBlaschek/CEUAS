#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__all__ = ['read_ragged_array_to_cube']


def read_ragged_cdm_to_array(filename, odb_codes=True, hours=False):
    """
    Read netCDF 4 HDF5 group data into a data cube on standard pressure levels
    this meant only for the CDS backend, not public

    duplicates are overwritten

    Args:
        filename (str): Filename of netcdf
        odb_codes (bool): Search for ODB or CDM variable codes
        hours (bool): enable hour split up (Caution does not work as expected), overwrites even more data

    Returns:
        xarray.Dataset : Dataset (plev, time)
    """
    import os
    import h5py
    import numpy as np
    import xarray as xr

    if not os.path.isfile(filename):
        raise IOError('File can not be found', filename)

    if odb_codes:
        var_names = {'t': 85, 'rh': 38, 'td': 36, 'dpd': 34, 'z': 117, 'dd': 106, 'ff': 107, 'u': 104, 'v': 105,
                     'q': 39}
    else:
        var_names = {'z': 1, 't': 2, 'u': 3, 'v': 4, 'dd': 111, 'ff': 112, 'td': 59, 'dpd': 299, 'rh': 29, 'p': 999}

    var_names = dict({j: i for i, j in var_names.items()})
    std_plevs = [10, 20, 30, 50, 70, 100, 150, 200, 250, 300, 400, 500, 700, 850, 925, 1000]
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
            print("Header Table found", igroup)
            #
            # What to read from here?
            #

        #
        # Observation Table
        #
        var_index = {}
        if 'observations_table' in f.keys():
            igroup = 'observations_table'
            print("Observations Table found", igroup)
            # get date_time
            time = f[igroup][dim][:]  # should be seconds since 1900 01 01

            if 'units' not in f[igroup][dim].attrs.keys():
                print("Datetime dimension", dim, "has no units!!! Assuming seconds")

            else:
                if 'seconds' not in f[igroup][dim].attrs['units'].decode():
                    raise RuntimeError("Datetime dimension", dim, "not in seconds since ????",
                                       f[igroup][dim].attrs['units'].decode())

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
            plev = f[igroup][lev][:].astype(np.int32) // 100  # requires hPa / because we only want std pressure levels
            iplev = np.in1d(plev, std_plevs)  # only std pressure (can have NaN)
            # get observed_variables
            variables = f[igroup]['observed_variable'][:]
            unique_variables = np.unique(variables)
            # get observation_value
            for i in unique_variables:
                if i not in var_names.keys():
                    print("Skipping.... ", i, "not in ", var_names.keys())
                    continue
                # only this variable + only standard pressure levels
                ivar = (variables == i) & (iplev)
                #
                # Experimental Code
                # Handle duplciates ?
                if False:
                    sec_var = None
                    # Can we use departures to find minimum ? -> but how to use the correct profile?
                    # this would potentially merge profiles
                    if 'era5fb' in f.keys() and 'an_depar@body' in f['era5fb'].keys():
                        sec_var = f['era5fb']['an_depar@body'][ivar][:]

                    ii = handle_duplicates(time[ivar], plev[ivar], f[igroup]['observation_value'][ivar][:],
                                           secondary_var=sec_var)
                    ivar = np.where(ivar)[0][ii]  # update only unique ivar
                else:
                    # check dates
                    timeplev = np.stack((time[ivar], plev[ivar],), axis=0)
                    # Duplicates will be overwritten if False
                    if False:
                        # Check for duplicates ? // take first one
                        u, ii, c = np.unique(timeplev, axis=1, return_counts=True, return_index=True)
                        ivar = np.where(ivar)[0][ii]  # update only unique ivar
                        print(igroup, i, "Duplicates", c[c > 1].size, "selected first")
                    else:
                        u, c = np.unique(timeplev, axis=1, return_counts=True)
                        print(igroup, i, "Duplicates", c[c > 1].size, "overwritten")

                var_index[i] = ivar
                itime, iobs = table_to_cube(time[ivar], ip[plev[ivar]], f[igroup]['observation_value'][ivar][:],
                                            hours=hours)
                if hours:
                    data[var_names[i]] = xr.DataArray(iobs, coords=([0, 12], std_plevs, date[ivar][itime]),
                                                      dims=('hour', 'plev', 'time'), name=var_names[i])
                else:
                    data[var_names[i]] = xr.DataArray(iobs, coords=(std_plevs, date[ivar][itime]),
                                                      dims=('plev', 'time'), name=var_names[i])
        #
        # Feedback Group
        #
        if 'era5fb' in f.keys():
            igroup = 'era5fb'
            print("ERA5 Feedback  found", igroup)
            #
            # Analysis
            #
            if 'an_depar@body' in f[igroup].keys():
                for i, ivar in var_index.items():
                    itime, iobs = table_to_cube(time[ivar], ip[plev[ivar]], f[igroup]['an_depar@body'][ivar][:],
                                                hours=hours)
                    if hours:
                        data[var_names[i] + '_an_dep'] = xr.DataArray(iobs,
                                                                      coords=([0, 12], std_plevs, date[ivar][itime]),
                                                                      dims=('hour', 'plev', 'time'),
                                                                      name=var_names[i] + '_an_dep')
                    else:
                        data[var_names[i] + '_an_dep'] = xr.DataArray(iobs, coords=(std_plevs, date[ivar][itime]),
                                                                      dims=('plev', 'time'),
                                                                      name=var_names[i] + '_an_dep')

            #
            # First Guess
            #
            if 'fg_depar@body' in f[igroup].keys():
                for i, ivar in var_index.items():
                    itime, iobs = table_to_cube(time[ivar], ip[plev[ivar]], f[igroup]['fg_depar@body'][ivar][:],
                                                hours=hours)
                    if hours:
                        data[var_names[i] + '_fg_dep'] = xr.DataArray(iobs,
                                                                      coords=([0, 12], std_plevs, date[ivar][itime]),
                                                                      dims=('hour', 'plev', 'time'),
                                                                      name=var_names[i] + '_fg_dep')
                    else:
                        data[var_names[i] + '_fg_dep'] = xr.DataArray(iobs, coords=(std_plevs, date[ivar][itime]),
                                                                      dims=('plev', 'time'),
                                                                      name=var_names[i] + '_fg_dep')

            #
            # Biascorr
            #
            if 'biascorr@body' in f[igroup].keys():
                for i, ivar in var_index.items():
                    itime, iobs = table_to_cube(time[ivar], ip[plev[ivar]], f[igroup]['biascorr@body'][ivar][:],
                                                hours=hours)
                    if hours:
                        data[var_names[i] + '_biascorr'] = xr.DataArray(iobs,
                                                                        coords=([0, 12], std_plevs, date[ivar][itime]),
                                                                        dims=('hour', 'plev', 'time'),
                                                                        name=var_names[i] + '_biascorr')
                    else:
                        data[var_names[i] + '_biascorr'] = xr.DataArray(iobs, coords=(std_plevs, date[ivar][itime]),
                                                                        dims=('plev', 'time'),
                                                                        name=var_names[i] + '_biascorr')

    data = xr.Dataset(data)
    data['plev'].attrs.update({'units': 'hPa', 'standard_name': 'air_pressure'})
    return data


def table_to_cube(time, plev, obs, hours=False, return_indexes=False):
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

    if hours:
        ihour = np.array(((time + 21600) % 86400) // 43200, dtype=np.int32)  # 0 or 12
        time = time // 86400

    xtime, jtime, itime = np.unique(time, return_index=True, return_inverse=True)
    if hours:
        ihour = np.array(((time + 21600) % 86400) // 43200, dtype=np.int32)  # 0 or 12
        if return_indexes:
            return jtime, xtime.size, ihour, plev, itime
        data = np.full((2, 16, xtime.size), np.nan, dtype=np.float32)
        data[ihour, plev, itime] = obs
    else:
        if return_indexes:
            return jtime, xtime.size, plev, itime
        data = np.full((16, xtime.size), np.nan, dtype=np.float32)
        data[plev, itime] = obs
    return jtime, data


def handle_duplicates(time, plev, obs, first=True, secondary_var=None):
    import numpy as np
    timeplev = np.stack((time, plev,), axis=0)
    u, ii, c = np.unique(timeplev, axis=1, return_counts=True, return_index=True)
    conflicts = np.where(c > 1)[0]
    print("Duplicates: ", conflicts.size)
    if conflicts.size > 0:
        for iconflict in conflicts:
            itx = np.where(timeplev == u[iconflict])[0]  # Find all instances
            if np.size(np.unique(obs[itx])) == 1:
                #
                # all the same value
                #
                ii[iconflict] = itx[0]
            else:
                #
                # different values
                #
                if first:
                    for i in itx:
                        if np.isfinite(obs[i]):
                            break
                elif secondary_var is not None:
                    i = np.argmin(secondary_var[itx])
                    i = itx[i]
                else:
                    i = itx[-1]  # LAST

                ii[iconflict] = i  # update indices to be applied
    return ii


def read_ragged_cdm(filename, odb_codes=True):
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
                print("Group not found: ", igroup)
                continue

            if dim in f[igroup].keys():
                # time
                time = f[igroup][dim][:]  # should be seconds since 1900 01 01

                if 'units' not in f[igroup][dim].attrs.keys():
                    print("Datetime dimension", dim, "has no units!!! Assuming seconds")
                else:
                    if 'seconds' not in f[igroup][dim].attrs['units'].decode():
                        raise RuntimeError("Datetime dimension", dim, "not in seconds since ????",
                                           f[igroup][dim].attrs['units'].decode())
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
    print("Converting to 2D ...")
    print(np.unique(data.observed_variable))
    data = dict(data.groupby(data.observed_variable))
    new = []
    # Loop and rename
    for i, j in data.items():
        if i not in var_names.keys():
            print("Skipping.... ", i, "not in ", var_names.keys())
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
            print(iname, "Non-unique index, removing duplicates...", j.index.duplicated().sum())
            print(j.loc[j.index.duplicated(keep=False)])
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
        print("Standard pressure levels need to be in hPa")

    # Pressure Array -> indices
    ip = np.zeros(1001, dtype=np.int32)
    for i, j in enumerate(std_plevs):
        ip[j] = i

    data = {}
    with h5py.File(filename, 'r') as f:
        # time
        time = f[dim][:]  # should be seconds since 1900 01 01
        if 'units' not in f[dim].attrs.keys():
            print("Datetime dimension", dim, "has no units!!! Assuming seconds")
        else:
            if 'seconds' not in f[dim].attrs['units'].decode():
                raise RuntimeError("Datetime dimension", dim, "not in seconds since ????",
                                   f[dim].attrs['units'].decode())
        # pressure
        plev = f[lev][:].astype(np.int32) // 100  # requires hPa
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
                     i not in ['string5', dim, lev, 'trajectory', 'trajectory_index', 'trajectory_label']]
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
