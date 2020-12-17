# -*- coding: utf-8 -*-

__all__ = ['read_ragged_array_to_cube', 'read_ragged_cdm_to_array']


def read_ragged_cdm_to_array(filename, odb_codes=True, hours=False, driver=None, add_index_array=False, verbose=0, debug=False, **kwargs):
    """
    Read netCDF 4 HDF5 group data into a data cube on standard pressure levels

    Args:
        filename (str): Filename of netcdf
        odb_codes (bool): Search for ODB or CDM variable codes
        hours (bool): split into 00 and 12 Z
        driver (str): h5py driver, default None, or core

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
    # standard pressure levels
    std_plevs = np.asarray([10, 20, 30, 50, 70, 100, 150, 200, 250, 300, 400, 500, 700, 850, 925, 1000])
    # Variables
    variables = ['observed_variable', 'observation_value', ]
    dim = 'date_time'
    lev = 'z_coordinate'
    # Pressure Array -> indices
    ip = np.zeros(1001, dtype=np.int32)
    for i, j in enumerate(std_plevs):
        ip[j] = i

    data = {}
    with h5py.File(filename, 'r', driver=driver) as f:
        #
        # Header Group ?
        #
        if 'header_table' in f.keys():
            igroup = 'header_table'
            if verbose >0: print("Header Table found", igroup)
            #
            # What to read from here?
            #
            
            # sources ? / meta data ?
            

        #
        # Observation Table
        #
        var_index = {}
        if 'observations_table' in f.keys():
            igroup = 'observations_table'
            
            if verbose >0: print("Observations Table found", igroup)
            # date time 
            time = f[igroup][dim][:]  # should be seconds since 1900 01 01
            # check time units
            if 'units' not in f[igroup][dim].attrs.keys():
                if verbose >0: print("Datetime dimension", dim, "has no units!!! Assuming seconds")

            else:
                if 'seconds' not in f[igroup][dim].attrs['units'].decode():
                    raise RuntimeError("Datetime dimension", dim, "not in seconds since ????",
                                       f[igroup][dim].attrs['units'].decode())
            # decode units or set to default
            if 'units' in f[igroup][dim].attrs.keys():
                time_unit = f[igroup][dim].attrs['units'].decode()
            else:
                time_unit = 'seconds since 1900-01-01 00:00:00'

            # convert to datetime64 from seconds using the units
            try:
                date = np.datetime64(" ".join(time_unit.split(' ')[-2:])) + time * np.timedelta64(1,
                                                                                                  time_unit.split(' ')[
                                                                                                      0][0])
            except:
                date = np.datetime64(" ".join(time_unit.split(' ')[-2:])) + time * np.timedelta64(1,
                                                                                                  time_unit.split(' ')[
                                                                                                      0][0].upper())

            # get pressures levels
            plev = f[igroup][lev][:] 
            iplev = np.in1d(plev, std_plevs*100)  # only std pressure (can have NaN)
            plev = plev.astype(np.int32) // 100  # needs to be integer for indexing, hPa as well 
            # get observed_variables -> variable codes
            # can be odb or cdm ?
            variables = f[igroup]['observed_variable'][:]
            unique_variables = np.unique(variables)
            #
            # loop variables and get observation_value's
            #
            for i in unique_variables:
                
                if i not in var_names.keys():
                    if verbose >0: print("Skipping.... ", i, "not in ", var_names.keys())
                    continue
                
                # only this variable + only standard pressure levels
                ivar = (variables == i) & (iplev)
                # combine datetime and pressure levels for unique record
                timeplev = np.stack((time[ivar], plev[ivar],), axis=0)
                # Duplicates will be overwritten if False
                if False:
                    # Check for duplicates ? // take first one
                    u, ii, c = np.unique(timeplev, axis=1, return_counts=True, return_index=True)
                    ivar = np.where(ivar)[0][ii]  # update only unique ivar
                    if verbose >0: print(igroup, i, "Duplicates", c[c > 1].size, "selected first")
                else:
                    u, c = np.unique(timeplev, axis=1, return_counts=True)
                    if verbose >0: print(igroup, i, "Duplicates", c[c > 1].size, "overwritten")
                
                if debug:
                    if c[c > 1].size > 0:
                    #if True:
                        raise RuntimeError()
                    
                # store variable indices for other groups
                var_index[i] = ivar
                # convert from table to cube
                itime, iobs = table_to_cube(time[ivar], ip[plev[ivar]], f[igroup]['observation_value'][ivar],
                                            hours=hours)
                # 
                if hours:
                    data[var_names[i]] = xr.DataArray(iobs, coords=([0, 12], std_plevs, date[ivar][itime]),
                                                      dims=('hour', 'plev', 'time'), name=var_names[i])
                else:
                    data[var_names[i]] = xr.DataArray(iobs, coords=(std_plevs, date[ivar][itime]),
                                                      dims=('plev', 'time'), name=var_names[i])
                if add_index_array:
                    itime, iivar = table_to_cube(time[ivar], ip[plev[ivar]], np.where(ivar)[0], hours=hours)
                    data[var_names[i]+'_index'] = xr.DataArray(iivar, coords=(std_plevs, date[ivar][itime]),
                                                      dims=('plev', 'time'), name=var_names[i]+'_index')
        #
        # Feedback Group
        #
        if 'era5fb' in f.keys():
            igroup = 'era5fb'
            if verbose >0: print("ERA5 Feedback  found", igroup)
            #
            # Analysis
            #
            if 'an_depar@body' in f[igroup].keys():
                for i, ivar in var_index.items():
                    itime, iobs = table_to_cube(time[ivar], ip[plev[ivar]], f[igroup]['an_depar@body'][ivar],
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
                    itime, iobs = table_to_cube(time[ivar], ip[plev[ivar]], f[igroup]['fg_depar@body'][ivar],
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
                    itime, iobs = table_to_cube(time[ivar], ip[plev[ivar]], f[igroup]['biascorr@body'][ivar],
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


def table_to_cube(time, plev, obs, hours=False, return_indexes=False):
    import numpy as np

    if hours:
        ihour = np.array(((time + 21600) % 86400) // 43200, dtype=np.int32) # 0 or 12
        time = time // 86400

    xtime, jtime, itime = np.unique(time, return_index=True, return_inverse=True)
    if hours:
        ihour = np.array(((time + 21600) % 86400) // 43200, dtype=np.int32) # 0 or 12
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