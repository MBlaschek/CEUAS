# -*- coding: utf-8 -*-

__all__ = ['standard_datetime_hours', 'saturation_water_vapor', 'total_precipitable_water',
           'vertical_interpolation', 'access_odb_table']


def standard_datetime_hours(data, dim='time', times=(0, 6, 12, 18), span=3, freq='6h', levels=None, **kwargs):
    from .std import align_datetime, to_hours
    from .. import config

    if levels is None:
        levels = config.std_plevels

    data = align_datetime(data.sel(plev=levels), dim=dim, times=times, span=span, freq=freq, **kwargs)
    return to_hours(data, dim=dim, times=times, **kwargs)


def access_odb_table(data, dim='time', pvar=None, sel=None):
    from numpy import unique
    import xarray  as xr

    data = data.sel(**{dim: sel}).set_coords(pvar)
    print(unique(data[dim].values))
    tmp = dict(data.groupby(dim))
    for ikey in tmp.keys():
        tmp[ikey] = tmp[ikey].swap_dims({dim: pvar}).assign_coords({dim: ikey})
    data = xr.concat(tmp.values(), dim=dim)
    return data


def saturation_water_vapor(temp, press=None, method='HylandWexler', precision=9, **kwargs):
    """ Calculate saturation water vapor

    Args:
        temp: temperatur
        method: method
        **kwargs:

    Returns:
        svp: satuartion water vapor pressure [Pa]
    """
    from numpy import around
    import xarray as xr
    from .convert import _conform
    from .esat import svp

    if not isinstance(temp, xr.DataArray):
        raise ValueError("Requires a DataArray", type(temp))

    evar = temp.copy()
    if press is not None:
        if isinstance(press, str):
            if press in evar.dims:
                press = evar[press].values
                press = _conform(press, evar.values.shape)
        elif isinstance(press, xr.DataArray):
            press = press.values
        else:
            pass

    evar.values = svp(temp.values, method=method, p=press, **kwargs)
    origin = 't' if press is None else 't,p'

    r_att = {'svp': method, 'standard_name': 'saturation_water_vapor_pressure',
             'long_name': 'saturation water vapor pressure', 'units': 'Pa',
             'origin': origin}
    if press is not None:
        r_att['enhancement_factor'] = "yes"

    r_att['precision'] = precision
    evar.attrs.update(r_att)
    evar.values = around(evar.values, decimals=precision)
    evar.name = 'td'
    return evar


def total_precipitable_water(data, dim='plev', levels=None, min_levels=8, fill_nan=True, **kwargs):
    """ Calculate total preciptable water by vertical integration of
    specific humidity profiles

     W = np.trapz(q, x=p) / 9.81  # kg/m2 == mm

    Args:
        data        (DataArray): Specific Humidity DataArray Class
        dim         (str): pressure dimension
        levels      (list): List of required pressure levels
        min_levels  (int):  minimum required levels for valid integration (exclusive with levels)
        fill_nan    (bool): convert missing numbers to 0

    Returns:
        DataArray : integrated TPW (vertical coordinate removed)

    Notes:
        Both Methods work fine
            W = np.trapz( q, x=p ) / 9.81
            W = np.sum( q * dp ) / 9.81   requires dp calculation dpres (NCL)
        The integral is with rho_water (assumed 1000) and neglected for conversion of m to mm
    """
    import xarray as xr
    from ..fun.xarray import xarray_function_wrapper
    from .tpw import tpw

    if not isinstance(data, xr.DataArray):
        raise ValueError("Requires a DataArray", type(data))

    if dim not in data.dims:
        raise ValueError("dim not found", dim)

    if 'standard_name' in data.attrs:
        if 'specific' not in data.attrs['standard_name']:
            raise RuntimeError("requires specific humidity, found:", data.attrs['standard_name'])
    else:
        RuntimeWarning("requires specfifc humidty, no standard_name present")

    data = data.copy()
    if levels is not None:
        data = data.sel(**{dim: levels})  # subset with only these pressure levels

    axis = data.get_axis_num(dim)

    if levels is not None:
        min_levels = len(levels)  # make sure we have exactly these levels

    counts = data.count(dim)

    data = xarray_function_wrapper(data, wfunc=tpw, dim=dim, axis=axis, min_levels=min_levels, fill_nan=fill_nan,
                                   pin=data[dim].values)
    # Update Name and Units
    data.name = 'tpw'
    r_att = {'standard_name': 'total_precipitable_water', 'min_levels': min_levels,
             'long_name': 'total precipitable water', 'units': 'mm',
             'cell_method': 'integral: specific_humidity/9.81'}

    data.attrs.update(r_att)
    data = data.to_dataset()
    data['counts'] = counts
    return data


def vertical_interpolation(data, dim='plev', levels=None, **kwargs):
    """ Apply a vertical log-pressure interpolation, no extrapolation

    todo: add table interpolation

    Args:
        data (DataArray): input data
        dim (str): vertical coordinate
        levels (list, ndarray): new vertical levels
        **kwargs:

    Returns:
        DataArray : interpolated Array
    """
    import xarray as xr
    import numpy as np
    from ..fun.interp import profile
    from .. import config

    if not isinstance(data, xr.DataArray):
        raise ValueError("Requires a DataArray", type(data))

    if dim not in data.dims:
        raise ValueError("Requires a valid dimension", dim, "of", data.dims)

    if levels is None:
        levels = config.std_plevels

    data = data.copy()
    axis = data.get_axis_num(dim)
    pin = data[dim].values
    values = np.apply_along_axis(profile, axis, data.values, pin, levels)
    data = data.reindex({dim: levels})  # can fail with duplicated values
    data.values = values
    cmethod = "%s: intp(%d > %d)" % (dim, len(pin), len(levels))
    if 'cell_method' in data.attrs:
        data.attrs['cell_method'] = cmethod + data.attrs['cell_method']
    else:
        data.attrs['cell_method'] = cmethod
    return data

#
# def adjust_dpd30(data, num_years=10, dim='time', subset=slice(None, '1995'), value=30,**kwargs):
#     """ Specifc Function to remove a certain value from the Histogram (DPD)
#
#     2. Criteria:
#         * Values at 30K are the end of the histogram
#         * before 1994
#         * clearly sticks out
#
#     Args:
#         data:
#         num_years:
#         dim:
#         subset (slice): time period, default before 1995 (mostly: US sondes)
#         value:
#         **kwargs:
#
#     Returns:
#
#     """
#     import numpy as np
#     from xarray import DataArray
#     from ..fun import message
#     from .manual import remove_spurious_values
#     if not isinstance(data, DataArray):
#         raise ValueError("Requires a DataArray", type(data))
#
#     if dim not in data.dims:
#         raise ValueError("DataArray class has no datetime dimension")
#
#     data = data.copy()
#     axis = data.dims.index(dim)
#     dates = data[dim].sel({dim: subset}).values
#     message("Using Subset %s" % str(subset), **kwargs)
#     values = data.sel({dim: subset}).values
#
#     if np.sum(np.isfinite(values)) == 0:
#         return data
#
#     count, mask = remove_spurious_values(dates, values, axis=axis, num_years=num_years, value=value, **kwargs)
#     dmask = DataArray(data=np.zeros_like(data.values), coords=data.coords, dims=data.dims,
#                       name=data.name + '_qc',
#                       attrs={'QC': count, 'DPD30':count,'flag_values': (0, 1, 2, 3),
#                              'flag_meanings': "no_qc good_data outside_range dpd30",
#                              'valid_range': (0, 3)})
#
#     dmask.loc[{dim: subset}] = mask*3   # DPD30
#     return dmask
