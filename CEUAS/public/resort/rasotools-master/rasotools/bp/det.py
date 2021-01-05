# -*- coding: utf-8 -*-
import numpy as np
from numba import njit

from ..fun import message

__all__ = ['detector', 'test', 'detector_ensemble']


def detector_ensemble(data, axis=0, ndist=None, nthres=None, nlevels=None, **kwargs):
    """ Detect breakpoints using multiple threshold values

    Args:
        data:
        axis:
        ndist:
        nthres:
        nlevels:
        **kwargs:

    Returns:

    """
    if ndist is None and nthres is None and nlevels is None:
        raise ValueError("requires either ndist, nthres or nlevels")

    tmp = np.full(data.shape, 0)
    n = -1
    if ndist is not None:
        if isinstance(ndist, int):
            ndist = np.linspace(180, 1460, ndist)
        # todo use prange
        for i in ndist:
            tmp += detector(data, axis=axis, dist=i, **kwargs)

    if nthres is not None:
        if isinstance(nthres, int):
            nthres = np.linspace(5, 100, nthres)
        # n = len(nthres) if n == -1 else n + len(nthres)
        for i in nthres:
            tmp += detector(data, axis=axis, thres=i, **kwargs)

    if nlevels is not None:
        if isinstance(nlevels, int):
            iaxis = 1 if axis == 0 else 0
            nlevels = range(1, data.shape[iaxis])
        # n = len(nlevels) if n == -1 else n + len(nlevels)
        for i in nlevels:
            tmp += detector(data, axis=axis, min_levels=i, **kwargs)

    return tmp  # 2*(tmp / np.max(tmp))  # (1)


def detector(data, axis=0, dist=365, thres=50, min_levels=3, use_slopes=False, use_first=False, use_exp=False, **kwargs):
    """ Detect breakpoints given some parameters

    Args:
        data (np.ndarray): snht test data
        axis (int): datetime dimension
        dist (int): minimum distance between breakpoints
        thres (int): threshold for significant levels
        min_levels (int): minimum of significant levels
        use_slopes (bool): Use first derivative for detection
        use_first (bool): Use the beginning of a peak, rather than the peak
        **kwargs:

    Returns:
        np.ndarray : breakpoint data
    """
    if not isinstance(data, np.ndarray):
        raise ValueError("requires an array: %s" % str(type(data)))

    breaks = (data >= thres).astype(int)  # 0, 1

    if data.ndim > 1:
        # number of breaks per timeunit
        jbreak = np.sum(breaks, axis=1 if axis == 0 else 0)
        # combine #-significant and sum of test
        jbreak = np.where(jbreak >= min_levels, 1, 0)
        ibreak = jbreak * np.sum(data, axis=1 if axis == 0 else 0)

    else:
        ibreak = data * breaks  # truncate at threshold

    # find local maxima (within distance)
    if use_slopes:
        imax = local_maxima(np.diff(ibreak), dist=dist, exp=use_exp)

    elif use_first:
        # find the beginning of breakpoint
        imax = local_maxima(ibreak, dist=dist, exp=use_exp)  # Maximum
        ifirst = []
        for i in imax:
            ifirst += [i + np.argmin(ibreak[i:i + dist])]
        imax = ifirst

    else:
        imax = local_maxima(ibreak, dist=dist, exp=use_exp)

    if len(imax) > 0:
        imax = np.asarray(imax)
        message("Breaks: " + str(imax), **kwargs)
        for i in imax:
            breaks[idx2shp(i, axis, data.shape)] += 2  # maximum Breakpoint

    return breaks


def percentile_detector(data, reference, freq, min_levels, percentiles=None, weights=None, verbose=0):
    # use running percentiles counts as input for test
    if percentiles is None:
        percentiles = np.arange(10, 90, 10)

    # subset to
    # calc. percentiles from reference
    # Apply percentiles to full dataset and count occurances /
    # when doe the occurances change?
    # nanpercentile ->

    return


def test(x, window, missing):
    """Standard Normal Homogeneity Test (SNHT)
    over a running window

    Wrapper function for numba_snhtmov
    @Leo Haimberger

    window = 2 years
    missing = 1/2 year

    """

    snhtparas = np.asarray([window, missing, 10])
    tsa = np.zeros(x.shape[0])
    tmean = np.zeros(x.shape[0])
    tsquare = np.zeros(x.shape[0])
    count = np.zeros(x.shape[0], dtype=np.int32)

    numba_snhtmov(np.squeeze(np.asarray(x)),
                  tsa,
                  snhtparas,
                  count,
                  tmean,
                  tsquare)

    return tsa


@njit
def numba_snhtmov(t, tsa, snhtparas, count, tmean, tsquare):
    """Standard Normal Homogeneity Test Moving Window

    t         = np.random.randn(1000)
    snhtparas = np.asarray([100,50,10])
    tsa       = np.zeros(1000)
    tmean     = np.zeros(1000)
    tsquare   = np.zeros(1000)
    count     = np.zeros(1000,dtype=np.int32)

    snhtmov2(t,tsa,snhtparasmcount,tmean,tsquare)

    Output: tsa
    """
    n = snhtparas[0]
    max_miss = snhtparas[1]
    # ninc=snhtparas[2]

    ni = t.shape[0]
    good = 0
    tmean[0] = 0.
    tsquare[0] = 0.
    for j in range(ni):
        count[j] = 0
        # compare_lists if nan ?
        if t[j] == t[j]:
            if good > 0:
                tmean[good] = tmean[good - 1] + t[j]
                tsquare[good] = tsquare[good - 1] + t[j] * t[j]
            else:
                tmean[good] = t[j]
                tsquare[good] = t[j] * t[j]
            good += 1
        if good > 0:
            count[j] = good - 1

    if good > n - 2 * max_miss:
        rm = int(n / 2)  # needs to be an integer
        # k 1460/2=730 - 650=80, n-80
        for k in range(rm - max_miss, ni - (rm - max_miss)):
            xm = k - rm  # 80-730
            if xm < 0:
                xm = 0
            xp = k + rm
            if xp > ni - 1:
                xp = ni - 1
            if (count[k] - count[xm] > rm - max_miss) and (count[xp] - count[k] > rm - max_miss):
                x = (tmean[count[k]] - tmean[count[xm]]) / (count[k] - count[xm])  # Mittelwert 1 Periode
                y = (tmean[count[xp]] - tmean[count[k]]) / (count[xp] - count[k])  # Mittelwert 2 Periode
                xy = (tmean[count[xp]] - tmean[count[xm]]) / (count[xp] - count[xm])  # Mittelwert ganze Periode

                sig = (tsquare[count[xp]] - tsquare[count[xm]]) / (count[xp] - count[xm])  # t*t ganze Periode
                if sig > xy * xy:
                    sig = np.sqrt(sig - xy * xy)  # standard deviation of the whole window
                    # n1 * (m1-m)**2 + n2 * (m2-m)**2 / stddev
                    tsa[k] = ((count[k] - count[xm]) * (x - xy) * (x - xy) + (count[xp] - count[k]) * (y - xy) * (
                            y - xy)) / (sig * sig)
                else:
                    tsa[k] = 0.
    return


@njit
def local_maxima(x, dist=365, exp=False):
    maxima = []  # Leere Liste
    # Iteriere von 2 bis vorletzten Element
    if exp:
        weights = np.interp(np.arange(dist - 1), [0, dist // 2, dist - 2], [1, 1, 0.5])
    else:
        weights = np.ones(dist - 1)

    # ist iterator integer
    for i in range(dist, len(x) - dist):
        # Element davor und danach größer
        if np.all((x[i] > x[slice(i + 1, i + dist)] * weights)) and np.all(
                (weights[::-1] * x[slice(i - dist+1, i)] < x[i])):
            maxima.append(i)
    return maxima


@njit
def local_minima(x, dist=365):
    minima = []  # Leere Liste
    # Iteriere von 2 bis vorletzten Element
    # ist iterator integer
    for i in range(dist, len(x) - dist):
        # Element davor und danach größer / zwei Spitzen können aber durch viel nix getrennt sein
        if np.all((x[i] < x[slice(i + 1, i + dist)])) and np.all((x[slice(i - dist, i)] > x[i])):
            minima.append(i)
    return minima


def find_constant(x, v, n=50):
    for i in range(0, len(x) - n):
        if np.sum(np.abs(np.diff(x[slice(i, i + n + 1)])) <= v) == n:
            break
    return i


def idx2shp(idx, axis, shape):
    index = [slice(None)] * len(shape)
    index[axis] = idx
    return tuple(index)
