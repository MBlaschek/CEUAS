# -*- coding: utf-8 -*-

__all__ = ['read_eccodes', 'read_cfgrib']


def read_cfgrib(filename, **kwargs):
    """ Read grib files using cfgrib

    Args:
        filename:
        **kwargs:

    Returns:

    """
    import xarray as xr
    try:
        import cfgrib
        import cf2cdm

    except ImportError:
        print("requires eccodes library, and xarray")
        raise ImportError()

    kwargs['engine'] = 'cfgrib'
    kwargs['backend_kwargs'] = kwargs.get('backend_kwargs', {})
    if 'indexpath' not in kwargs['backend_kwargs']:
        kwargs['backend_kwargs'].update({'indexpath':''})  # no indexfiles

    data = xr.open_dataset(filename, **kwargs)
    return cf2cdm.translate_coords(data, cf2cdm.CDS)


def read_eccodes(filename, keys=None, **kwargs):
    """ Read a GRIB file into DataArray class

    Expects: var x date x level x lat x lon

    Args:
        filename (str):
        keys (dict):

    Returns:

    """
    import os
    from time import time
    import numpy as np
    import xarray as xr
    try:
        import eccodes

    except ImportError:
        print("This function requires eccodes from ECMWF installed")
        print("eccodes requires python 2, for python 3 use cfgrib")
        raise ImportError()

    from ..fun import message

    if not os.path.isfile(filename):
        raise IOError("File not found: %s" % filename)

    data = {}
    attrs = {}
    lat = None
    lon = None

    start_time = time()
    try:
        with eccodes.GribFile(filename) as grib:
            for msg in grib:
                # shortName
                ivar = msg.get('shortName')
                if lat is None:
                    # latlon = msg.get('latLonValues')  # same as grid Lat, Lon
                    lat = msg.get('distinctLatitudes')  # longitudes
                    lon = msg.get('distinctLongitudes')  # latitudes
                    # Nj,i
                    ni = msg.get('Ni')  # number of Columns
                    nj = msg.get('Nj')  # number of Rows ?

                ilevel = msg.get('level')  # vertical level ?
                values = msg.get('values')  # 1D array
                idate = "%04d-%02d-%02dT%02d:%02d:%02d" % (msg.get('year'), msg.get('month'), msg.get('day'),
                                                           msg.get('hour'), msg.get('minute'), msg.get('second'))

                if ivar not in data.keys():
                    data[ivar] = {}
                    attrs[ivar] = {}
                    attrs[ivar]['units'] = msg.get('units')
                    attrs[ivar]['long_name'] = msg.get('name')

                if idate not in data[ivar].keys():
                    data[ivar][idate] = {}

                data[ivar][idate][ilevel] = values.reshape(nj, ni)  # how to figure out which comes first ?
                message("%s [%s] [%s]" % (ivar, idate, str(ilevel)), **kwargs)
    except:
        raise

    # latlon = latlon.reshape(nj, ni)
    # lon = latlon[0, :]
    # lat = latlon[:, 0]
    # lon = lon.reshape(nj,ni)[0,:]  # Lons
    # lat = lat.reshape(nj,ni)[:,0]  # Lats

    for ivar in data.keys():
        dates = data[ivar].keys()
        levels = data[ivar][dates[0]].keys()
        order = ('date', 'level', 'latitude', 'longitude')
        tmp = []
        dim_attrs = {'level': {'units': 'hPa'}, 'longitude': {'units': 'degrees_east'},
                     'latitude': {'units': 'degrees_north'}}

        dims = {'date': np.array([np.datetime64(i) for i in dates]), 'level': np.array(levels).astype(float),
                'longitude': lon, 'latitude': lat}

        for idate in dates:
            tmp += [np.stack(data[ivar][idate].values())]

        message("DataArray ", ivar, **kwargs)
        data[ivar] = xr.DataArray(ivar, np.array(tmp), order, dims, attrs=attrs[ivar], dim_attrs=dim_attrs, axes=['T', 'Z', 'Y', 'X'])
        data[ivar].sort(verbose=0)  # funktioniert das?

    message("--- %s seconds ---" % (time() - start_time), **kwargs)
    # missing global attributes and history
    return data  # needs to be an xarray
