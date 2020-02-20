#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__all__ = ['read_ragged_array_to_cube']


def read_ragged_array_to_cube(filename, dim='time', lev='plev', variables=None, std_plevs=None, **kwargs):
    """
    The function reads a CF compliant upper air data file from the CDS backend and transforms the data into
    a three dimensional data cube with dimensions hour (2), plev (16), time (days).
    This representation is convienient for statistical breakpoint analysis.

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
