# -*- coding: utf-8 -*-
import xarray as xr

from ..fun import message


def update_shape(data, order):
    # xarray.DataArray -> xarray.DataArray
    iorder = ()
    for idim in order:
        if idim in data.dims:
            iorder += (idim,)

    if len(iorder) > 0:
        return data.transpose(*iorder)
    return data


def add_id(ds, ivars, invert=False, fuzzy=False):
    if 'ident' in ds.attrs:
        ds.coords['sonde'] = ds.attrs['ident']
    elif 'station_id' in ds.attrs:
        ds.coords['sonde'] = ds.attrs['station_id']
    else:
        ds.coords['sonde'] = ds.encoding['source'].split('/')[-2]

    if ivars is not None:
        # todo add fuzzy search
        if invert:
            ds = ds[[i for i in list(ds.data_vars) if i not in ivars]]
        else:
            ds = ds[[i for i in list(ds.data_vars) if i in ivars]]
    return ds


def align_concat_fill(*args, **kwargs):
    #
    # Put Datasets together
    #
    import numpy as np
    import xarray as xr

    xr.set_options(keep_attrs=True)
    #
    # align dimensions
    #
    args = xr.align(*args, join='outer', copy=False)
    message("Alignment completed", **kwargs)
    #
    # align datavars
    #
    dummy = xr.Dataset()
    for iset in args:
        ivars = list(iset.data_vars)
        for ivar in ivars:
            if ivar not in dummy.data_vars:
                dummy[ivar] = xr.full_like(iset[ivar], np.nan)  # copies attributes as well

    message("Initialization completed", dummy, **kwargs)
    #
    # add missing variables
    #
    for iset in args:
        for ivar in list(dummy.data_vars):
            if ivar not in iset.data_vars:
                iset[ivar] = dummy[ivar]
    #
    # concat
    #
    out = xr.concat(args, **kwargs)
    message("Concatenation completed", **kwargs)
    #
    # restore some Attributes, as they are removed by concat 'equal'
    #
    for ivar in out.data_vars:
        out[ivar].attrs.update(dict(dummy[ivar].attrs))
    return out


def read_add_process(ifile, ivars, invert=False, **kwargs):
    ds = xr.open_dataset(ifile, **kwargs.get('xarray', {'engine': 'h5netcdf'}))
    try:
        ds = add_id(ds, ivars, invert=invert)  # sonde coords
        ds.load()
        #
        # execute function
        #
        kwargs['mname'] = ifile
        if 'exec' in kwargs.keys():
            ds = process_function(ds, **kwargs)
    except Exception as e:
        print("Error: ", ifile)
        print(repr(e))

    finally:
        ds.close()
    return ds


def process_function(data, **kwargs):
    #
    # Copy attributes
    #
    gattrs = data.attrs.copy()
    #0
    # Module
    #
    iname = kwargs['exec'].split('.')
    imodule = __import__(iname[0])
    #
    # Function
    #
    ifunc = getattr(imodule, iname[1])
    message("Executing: ", kwargs['exec'], imodule.__file__, ifunc.__name__, **kwargs)
    #
    # Execution
    #
    data = ifunc(data, **kwargs)
    #
    # Attributes
    #
    data.attrs.update(gattrs)
    return data
