#!/usr/bin/env python
# -*- coding: utf-8 -*-


def fix(unit, values):
    try:
        import pint
        ureg = pint.UnitRegistry()
        if unit != 'Pa':
            values = values * float(ureg(unit) / ureg('Pa'))

    except ImportError:
        if unit != 'Pa':
            values = values * 100.
    return values


def check_level_units(data, dim='level', want='Pa'):
    from xarray import DataArray, Dataset
    if not isinstance(data, (Dataset, DataArray)):
        raise ValueError("Requires a xarray Dataset or DataArray")

    print("[vars]", list(data.data_vars))
    converted = False
    if dim in data.coords:
        if 'units' in data[dim].attrs:
            i = data[dim].attrs['units']
            if i.lower() != want.lower():
                print("[converting]", dim, " ", i, " => ", want)
                data.load()  # copy to Memory
                data = data.assign_coords(**{dim: fix(i, data[dim].values)})
                data[dim].attrs.update({'units': want})
                if 'history' in data.attrs:
                    data.attrs[
                        'history'] = datetime.datetime.now().isoformat() + " level(" + i + " to " + want + ")\n" + \
                                     data.attrs['history']
                converted = True
            else:
                print("[doing nothing] " + want)
        else:
            print("[doing nothing] no units")
    else:
        print("[doing nothing] no level coord")
    return data, converted


def usage(name):
    print("""
%s [files]
Fix level dimension values to Pa
""" % name)


if __name__ == "__main__":
    import sys
    import xarray
    import datetime

    lev = 'level'

    if len(sys.argv) == 1:
        usage(sys.argv[0])
    else:
        for ifile in sys.argv[1:]:
            print("[file]", ifile)
            data = xarray.open_dataset(ifile, decode_times=False)
            data, status = check_level_units(data, lev, 'Pa')
            if status:
                data.to_netcdf(ifile)
            data.close()
        print("[Done]")
