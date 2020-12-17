#!/usr/bin/env python
# -*- coding: utf-8 -*-

_coords = {
    'X': {
        "name": "lon",
        "axis": "X",
        "direction": "increasing",
        "standard_name": "longitude",
        "long_name": "longitude",
        "valid min/max": "0/360",
        "units": "degrees_east",
        "positive": "east"
    },
    'Y': {
        "name": "lat",
        "axis": "Y",
        "direction": "increasing",
        "standard_name": "latitude",
        "long_name": "latitude",
        "valid min/max": "-90/90",
        "units": "degrees_north",
        "positive": "north"
    },
    'T': {
        "name": "time",
        "axis": "T",
        "standard_name": "time",
        "long_name": "time"
    },
    'Z': {
        "name": "plev",
        "axis": "Z",
        "direction": "decreasing",
        "standard_name": "air_pressure",
        "long_name": "pressure",
        "units": "Pa",
        "positive": "down"
    },
    'H': {
        "name": "height",
        "axis": "Z",
        "direction": "increasing",
        "standard_name": "height",
        "long_name": "height",
        "units": "m",
        "positive": "up"
    },
    'D': {
        "name": "depth",
        "axis": "Z",
        "direction": "increasing",
        "standard_name": "depth",
        "long_name": "depth",
        "units": "m",
        "positive": "down"
    },
    'R': {
        "name": "realization",
        "axis": "E",
        "direction": "increasing",
        "standard_name": "realization",
        "long_name": "realization",
        "units": "1"
    },
    'G': {
        "name": "gph",
        "axis": "Z",
        "standard_name": "geopotential_height",
        "long_name": "geopotential height",
        "units": "m",
        "positive": "up"
    }
}

_vars = {
    'tvar': {
        "name": "t",
        "standard_name": "air_temperature",
        "long_name": "air temperature",
        "short_name": "ta",
        "units": "K"
    },
    'tdvar': {
        "name": "td",
        "standard_name": "dew_point_temperature",
        "long_name": "dew point temperature",
        "short_name": "tda",
        "units": "K"
    },
    'dpdvar': {
        "name": "dpd",
        "standard_name": "dew_point_depression",
        "long_name": "dew point depression",
        "short_name": "dpd",
        "units": "K"
    },
    'gvar': {
        "name": "z",
        "standard_name": "geopotential",
        "long_name": "Geopotential",
        "short_name": "z",
        "units": "m2 s-2"
    },
    'qvar': {
        "name": "q",
        "long_name": "Specific humidity",
        "standard_name": "specific_humidity",
        "short_name": "hus",
        "units": "kg kg-1"
    },
    'rvar': {
        "name": "rh",
        "long_name": "Relative humidity",
        "standard_name": "relative_humidity",
        "short_name": "hur",
        "units": "1"
    },
    'uvar': {
        "name": "u",
        "long_name": "U wind component",
        "standard_name": "eastward_wind",
        "short_name": "ua",
        "units": "m s-1"
    },
    'vvar': {
        "name": "v",
        "long_name": "V wind component",
        "standard_name": "northward_wind",
        "short_name": "va",
        "units": "m s-1"
    },
    'ffvar': {
        "name": "ws",
        "standard_name": "wind_speed",
        "long_name": "Wind Speed",
        "units": "m s-1",
        "short_name": "ws"
    },
    'ddvar': {
        "name": "wd",
        "standard_name": "wind_to_direction",
        "long_name": "Wind Direction",
        "units": "degree",
        "short_name": "wd"
    }
}


def _inquire_variable(data, name, isdim=False):
    if isdim:
        dim = {'T': ['date', 'time', 'datetime', 'hours since', 'days since'],
               'Z': ['lev', 'pres', 'pressure', 'air_pressure', 'plev'],
               'X': ['lon', 'long', 'longitude', 'east'],
               'Y': ['lat', 'latt', 'latitude', 'north']}
        # R, D, H

        if 'axis' in data.attrs:
            if data.attrs['axis'] in dim.keys():
                return data.attrs['axis']  # T, Z ,X, Y

        if 'standard_name' in data.attrs:
            for idim, ival in dim.items():
                if data.attrs['standard_name'] in ival:
                    return idim

        for idim, ival in dim.items():
            for ikey in ival:
                for iatt in data.attrs.values():
                    if ikey in str(iatt):
                        return idim

        if data.name != '':
            for idim, ival in dim.items():
                if data.name in ival:
                    return idim

    else:
        varis = {'tvar': ['temp', 't', 'temperature', 'air_temperature'],
                 'rvar': ['rhumi', 'rh', 'relative_humidity', 'rel. humidity', 'relative humidity'],
                 'qvar': ['qhumi', 'q', 'specific_humidity', 'spec. humidity', 'kg kg**-1'],
                 'tdvar': ['tdew', 'td', 'dew_point', 'dew point', 'dewpoint'],
                 'dpdvar': ['dpd', 'dew_point_depression', 'dewpoint_depression'],
                 'uvar': ['uwind', 'u', 'horizontal_wind_speed'],
                 'vvar': ['vwind', 'v', 'meridonal_wind_speed'],
                 'ffvar': ['winds', 'wind_speed', 'ws', 'ff'],
                 'ddvar': ['windd', 'wind_direction', 'wd', 'dd', 'wind_to_direction']}

        if 'standard_name' in data.attrs:
            for idim, ival in varis.items():
                if data.attrs['standard_name'] in ival:
                    return idim

        for idim, ival in varis.items():
            for ikey in ival:
                if ikey in data.attrs.values():
                    return idim

        for idim, ival in varis.items():
            for ikey in ival:
                if len(ikey) <= 2:
                    continue

                for iatt in data.attrs.values():
                    if ikey in str(iatt).lower():
                        return idim

        if data.name != '':
            for idim, ival in varis.items():
                if data.name in ival:
                    return idim
    return name


def apply_cfunits(data, u_in, u_out, **kwargs):
    from rasotools.fun import message
    if u_in is not None and u_out is not None:
        if u_in != u_out:
            import pint
            ureg = pint.UnitRegistry()
            data.values = data.values * float(ureg(u_in).to(ureg(u_out)).magnitude)
            message("Converting ", u_in, " to ", u_out, " by ", float(ureg(u_in).to(ureg(u_out)).magnitude), **kwargs)
    else:
        message("Warning no units: ", data.name, u_in, u_out, **kwargs)
    return data


def apply_cflongitude(data, idim):
    import numpy as np
    data[idim].values = np.where(data[idim].values < 0, data[idim].values + 360., data[idim].values)


def apply_cfnames(data, name, src=None, method=None, common_coords=None, common_vars=None, **kwargs):
    import numpy as np
    from xarray import DataArray
    from rasotools.fun import dict_in_dict, dict2str, message

    if not isinstance(data, DataArray):
        raise ValueError('Requires a DataArray, ', type(data))

    if common_coords is None:
        common_coords = _coords

    if common_vars is None:
        common_vars = _vars

    data = data.copy()
    # check dims
    # check attributes
    for idim in data.coords:
        jdim = _inquire_variable(data[idim], idim, isdim=True)
        if idim == 'hour' or idim == 'hours':
            continue

        if jdim in common_coords.keys():
            # units
            if 'units' in common_coords[jdim].keys():
                data = data.assign_coords(**{idim: apply_cfunits(data[idim],
                                                                 data[idim].attrs.get('units'),
                                                                 common_coords[jdim]['units'],
                                                                 **kwargs)})
            # attributes
            _attrs = data[idim].attrs.copy()
            data[idim].attrs.update(common_coords[jdim])
            if 'direction' in data[idim].attrs:
                try:
                    if np.all(np.sort(data[idim].values) == data[idim].values):
                        data[idim].attrs['direction'] = 'increasing'
                    else:
                        data[idim].attrs['direction'] = 'decreasing'
                except:
                    pass
            if jdim == 'X':
                apply_cflongitude(data, idim)

            # old name
            data[idim].attrs['name'] = idim
            _attrsx = data[idim].attrs.copy()
            # rename Dimension
            data = data.rename({idim: common_coords[jdim]['name']})
            # Report
            message("Renaming ", idim, " to ", common_coords[jdim]['name'], **kwargs)
            message("Attributes: ", dict2str(dict_in_dict(_attrsx, _attrs)), **kwargs)

    iname = _inquire_variable(data, name, isdim=False)
    if iname in common_vars.keys():
        _attrs = data.attrs.copy()
        data = apply_cfunits(data, data.attrs.get('units'), common_vars[iname]['units'])
        data.attrs.update(common_vars[iname])
        data.attrs['name'] = name  # old name
        if src is not None:
            data.attrs['src'] = src
        if method is not None:
            data.attrs['method'] = method
        name = common_vars[iname]['name']
        # Report
        message("Renaming ", data.attrs['name'], " to ", name, **kwargs)
        message("Attributes: ", dict2str(dict_in_dict(data.attrs, _attrs)), **kwargs)

    return name, data


def cfprofile(data, dim='time', starttime='days since 1979-01-01 00:00:00'):
    from numpy import arange
    from rasotools.met import std

    for idim in data.dims:
        if idim in ['time', 'date']:
            dim = idim
    #
    # remove hour
    #
    if 'hour' in data.dims:
        data = std.from_hours(data, dim=dim)
        data = data.drop('hour')
    #
    # fix time dimension
    #
    data[dim].encoding['units'] = starttime
    data[dim].attrs['standard_name'] = 'time'
    data['profile'] = (dim, arange(data[dim].size, dtype='i4'))
    data = data.swap_dims({dim: 'profile'})
    data['profile'].attrs['cf_role'] = 'profile_id'
    if 'lon' not in data.data_vars and 'lon' not in data.dims:
        if 'station_lon' in data.attrs.keys():
            data['lon'] = ('profile', [data.attrs['station_lon']] * data[dim].size)
        elif 'lon' in data.attrs.keys():
            data['lon'] = ('profile', [data.attrs['lon']] * data[dim].size)
        else:
            pass
        if 'station_lat' in data.attrs.keys():
            data['lat'] = ('profile', [data.attrs['station_lat']] * data[dim].size)
        elif 'lat' in data.attrs.keys():
            data['lat'] = ('profile', [data.attrs['lat']] * data[dim].size)
        else:
            pass
        if 'lon' in data.data_vars and 'lat' in data.data_vars:
            data = data.set_coords(['lon', 'lat'])
    return data


def main():
    from rasotools import _getlibs
    import sys
    import xarray as xr
    option = {}
    if len(sys.argv) > 1:
        for ifile in sys.argv[1:]:
            if ifile[0:2] == '--':
                # option
                itmp = ifile.replace('--', '').split('=')
                option[itmp[0]] = itmp[1]
                continue

            print(ifile)
            #
            # Open
            #
            data = xr.open_dataset(ifile)
            out = {}
            # Apply cf convention

            for ivar in list(data.data_vars)[:]:
                src = None
                if '_fg_' in ivar:
                    src = 'first_guess'
                if '_an_' in ivar:
                    src = 'analysis'
                method = None
                if '_dep' in ivar:
                    method = 'time: departure'
                if '_bias' in ivar:
                    method = 'time: bias adjustment'
                _, out[ivar] = apply_cfnames(data[ivar], ivar, verbose=1, src=src, method=method, **option)
            #
            # Write with suffix
            #
            out = xr.Dataset(out)
            out.attrs.update(data.attrs)
            out.attrs['Conventions'] = 'CF-1.7'
            out.to_netcdf(ifile.replace('.nc', '_cf.nc'))
            print(ifile.replace('.nc', '_cf.nc'))
    else:
        print(sys.argv[0] + " [files] ")
        print("Apply CF-1.7 convention to NetCDF radiosonde data")
        print(_getlibs())


if __name__ == '__main__':
    main()
