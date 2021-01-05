# -*- coding: utf-8 -*-

__all__ = ['refraction']


def refraction(temp, pres, vp):
    import xarray as xr

    # T, p, vp ?
    # dry refractivity or wet ? from radiosonde
    c1 = 0
    c2 = 0
    N = c1 * pres / temp - c2 * vp / (temp*temp)
    # data['ref'] = data.eval("c1*p/t - c2*vp/t**2")
    return N
