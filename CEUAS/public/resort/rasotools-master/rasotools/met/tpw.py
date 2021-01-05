# -*- coding: utf-8 -*-
import numpy as np

__all__ = ['tpw']


def tpw(values, pin, axis=0, min_levels=8, fill_nan=True):
    """ Calculate total preciptable water by vertical integration of
    specific humidity profiles

     W = np.trapz(q, x=p) / 9.81 # kg/m2 == mm

    Args:
        values      (numpy.ndarray):    Specific Humidity [kg/kg]
        pin         (numpy.ndarray):    Pressure [Pa]
        axis        (int):      Axis for vertical integration
        min_levels  (int):      minimum required levels for valid integration (exclusive with levels)
        fill_nan    (bool):     convert missing numbers to 0

    Returns:
        numpy.ndarray : integrated TPW (vertical coordinate removed)

    Notes:
        Both Methods work fine
            W = np.trapz( q, x=p ) / 9.81

            W = np.sum( q*dp ) / 9.81   requires dp calculation dpres (NCL)
    """
    # Integration over pressure levels
    itx = np.isfinite(values)
    if fill_nan:
        values = np.where(itx, values, 0)  # Otherwise integration will be nan

    values = np.trapz(values, x=pin, axis=axis) / 9.81 # Integration
    # if not enough levels are present -> NAN
    values = np.where(np.sum(itx, axis=axis) >= min_levels, values, np.nan)
    return values


def iwv(values, pin, axis=0, min_levels=8, fill_nan=True):
    #
    # Pa * Pa
    # N / m2
    #
    return values

def dpres(p, sfc=101300., top=1000.):
    """Calculate delta p intervals for integration
    Copy from NCL function dpres

    Input:       p      [Array]
    Optional:
                 sfc    Surface Pressure
                 top    Column upper most pressure
    Output:      dp     delta p
    """
    p = np.asarray(p)
    i = np.argsort(p)  # sort ascending
    p = p[i][::-1]  # descending
    if np.min(p) < top:
        top = np.min(p)

    if np.max(p) > sfc:
        sfc = np.max(p)

    dp = np.zeros(len(p))
    # first
    dp[0] = sfc - p[0] + (p[0] - p[1]) / 2.
    # bottom to top
    for j in range(1, len(p) - 1):
        dp[j] = (p[j - 1] - p[j]) / 2. + (p[j] - p[j + 1]) / 2.
    # Letzte: kleinster p wert
    dp[-1] = (p[j] - p[j + 1]) / 2. + (p[j + 1] - top) / 2.
    # selbe ordnung wie p am anfang?
    return dp[i][::-1]
