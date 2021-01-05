# -*- coding: utf-8 -*-

__all__ = ['polar_to_uv', 'uv_to_polar']


def polar_to_uv(speed, direction, suffix=''):
    """Calculate the U, V wind vector components from the speed and direction.

    Args:
        data (Dataset): Xarray Input data
        speed (str):  The wind speed (magnitude)
        direction (str): The wind direction, specified as the direction from which the wind is
                blowing (0-2 pi radians or 0-360 degrees), with 360 degrees being North.
        suffix(str): output name

    Returns:
        DataArray, DataArray : u,v wind components

    Notes: from MetPy (https://unidata.github.io/MetPy)
    """
    import numpy as np
    from xarray import DataArray
    if not isinstance(speed, DataArray):
        raise ValueError("Requires a DataArray")

    if not isinstance(direction, DataArray):
        raise ValueError("Requires a DataArray")

    if 'units' in direction.attrs.keys():
        if direction.attrs['units'] != 'rad':
            direction = direction.copy() * np.pi / 180
    else:
        print('Warning requires rad for Wind direction, no units found')

    u = -speed * np.sin(direction)
    u.attrs.update({'standard_name': 'eastward_wind', 'units': 'm s-1', "long_name": "U wind component"})
    u.name = 'u' + suffix
    v = -speed * np.cos(direction)
    v.attrs.update({'standard_name': 'northward_wind', 'units': 'm s-1', "long_name": "V wind component"})
    v.name = 'v' + suffix
    return u, v


def uv_to_polar(u, v, suffix='', convention='from'):
    """Compute the wind direction from u and v-components.

    Args:
        u(DataArray): Wind component in the X (East-West) direction
        v(DataArray): Wind component in the Y (North-South) direction
        suffix (str): output name
        convention (str): Convention to return direction. 'from' returns the direction the wind is coming from (meteorological convention). 'to' returns the direction the wind is going towards (oceanographic convention). Default is 'from'.
    Returns:
        DataArray, DataArray : ws, wd

    Notes:
        In the case of calm winds (where `u` and `v` are zero), this function returns a direction of 0.
        from MetPy (https://unidata.github.io/MetPy)
    """
    import numpy as np
    from xarray import DataArray
    if not isinstance(u, DataArray):
        raise ValueError("Requires a DataArray")
    if not isinstance(v, DataArray):
        raise ValueError("Requires a DataArray")

    ws = np.sqrt(u ** 2 + v ** 2)
    wd = 90 - np.arctan2(-v, -u) * 180 / np.pi
    # Handle oceanographic convection
    if convention == 'to':
        wd -= 180
    wd = wd.where(wd > 0, wd + 360)
    # avoid unintended modification of `pint.Quantity` by direct use of magnitude
    wd = wd.where((u != 0.) & (v != 0.), 0.)
    wd.attrs.update(
        {'convention': convention, "standard_name": "wind_" + convention + "_direction", "long_name": "Wind Direction",
         "units": "degree"})
    ws.attrs.update({"standard_name": "wind_speed", "long_name": "Wind Speed", "units": "m s-1", })
    wd.name = 'wd' + suffix
    ws.name = 'ws' + suffix
    return ws, wd
