#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__version__ = '0.1'
__author__ = 'MB'
__status__ = 'dev'
__date__ = 'Mit Jän 15 15:30:59 CET 2020'
__institute__ = 'Univie, IMGW'
__github__ = 'git@github.com:MBlaschek/CEUAS.git'
__doc__ = """
Radiosonde Homogenisation Software v%s
Maintained by %s at %s
Github: %s [%s]
Updated: %s
""" % (__version__, __author__, __institute__, __github__, __status__, __date__)

import sys

import numpy as np
from numba import njit

#
# Suppress warnings
#
np.seterr(invalid='ignore')

# in Pa
std_plevels = np.array(
    [1000., 2000., 3000., 5000., 7000., 10000., 15000., 20000., 25000., 30000., 40000., 50000., 70000.,
     85000., 92500., 100000.])


def usage(name):
    return """
Run standardized radiosonde homogenisation software on CDM compliant file

{} -f [CDM file]

    """.format(name)


###############################################################################
#
# Helper functions
#
###############################################################################
def now(timespec='auto'):
    """ Datetime string

    Returns:
        str : datetime now
    """
    import datetime
    return datetime.datetime.now().isoformat(timespec=timespec)


def _print_string(*args, adddate=False, **kwargs):
    if adddate:
        return "[" + now() + "] " + " ".join([str(i) for i in args])
    else:
        return " ".join([str(i) for i in args])


def message(*args, mname=None, verbose=0, level=0, logfile=None, **kwargs):
    if logfile is not None:
        # with open(kwargs['filename'], 'a' if not kwargs.get('force', False) else 'w') as f:
        with open(logfile, 'a') as f:
            f.write(_print_string(*args, **kwargs) + "\n")

    elif verbose > level:
        text = _print_string(*args, **kwargs)
        if mname is not None:
            text = "[%s] " % mname + text

        print(text)
    else:
        pass


def conform(data, shape):
    """ Make numpy array conform to a certain shape

    Args:
        data:
        shape:

    Returns:

    """
    if not isinstance(data, np.ndarray):
        raise ValueError('Requires a numpy array')

    if not isinstance(shape, (tuple, list)):
        raise ValueError('Requires a tuple or list')

    data = data.copy()
    n = data.shape

    assert np.any([i in shape for i in n]), "Shapes do not allign?!"

    for i, j in enumerate(shape):
        if j not in n:
            data = np.expand_dims(data, axis=i)
    return data


def nancount(x, axis=0, keepdims=False):
    """

    Args:
        x (ndarray): input dataset
        axis (int): axis
        keepdims (bool): keep dimensions
    """
    return np.sum(np.isfinite(x), axis=axis, keepdims=keepdims)


def nanfunc(data, n=130, axis=0, nmax=1460, borders=0, ffunc=None, flip=False, fargs=(), **kwargs):
    """ Nan omitting function (numpy)

    Args:
        data (np.ndarray): dataset including NaN
        n (int): minimum sample size
        axis (int): datetime axis
        nmax (int): maximum sample size
        borders (int): border sample to ignore
        ffunc (callable): function to call
        flip (bool): reverse dataset before applying the function
        args (tuple): function arguments

    Returns:
        np.ndarray : func of values at axis, with sample size, borders and maximum
    """
    if ffunc is None:
        ffunc = np.nanmean
    return np.apply_along_axis(sample, axis, data, n, nmax, ffunc, borders=borders, flip=flip, fargs=fargs)


def sample(values, nmin, nmax, func, borders=0, flip=False, fargs=(), **kwargs):
    itx = np.isfinite(values)
    n = itx.sum()
    j = 0
    if n > nmax:
        if n > (nmax + borders):
            j = borders
        if flip:
            return func(np.flip(values[itx])[j:(nmax + j)], *fargs)  # reversed
        return func(values[itx][j:(nmax + j)], *fargs)  # normal

    elif n < nmin:
        # raises all nan warnings !!!
        return func(values, *fargs) * np.nan

    else:
        if n > (nmin * 2 + borders):
            j = borders
        if flip:
            return func(np.flip(values[j:]), *fargs)
        return func(values[j:], *fargs)


###############################################################################
#
# Detection of Breakpoints
#
###############################################################################
def detector(data, axis=0, dist=365, thres=50, min_levels=3, use_slopes=False, use_first=False, **kwargs):
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
        imax = local_maxima(np.diff(ibreak), dist=dist)

    elif use_first:
        # find the beginning of breakpoint
        imax = local_maxima(ibreak, dist=dist)  # Maximum
        ifirst = []
        for i in imax:
            ifirst += [i + np.argmin(ibreak[i:i + dist])]
        imax = ifirst

    else:
        imax = local_maxima(ibreak, dist=dist)

    if len(imax) > 0:
        imax = np.asarray(imax)
        message("Breaks: " + str(imax), **kwargs)
        for i in imax:
            breaks[idx2shp(i, axis, data.shape)] += 2  # maximum Breakpoint

    return breaks


def idx2shp(idx, axis, shape):
    index = [slice(None)] * len(shape)
    index[axis] = idx
    return tuple(index)


@njit
def local_maxima(x, dist=365):
    maxima = []  # Leere Liste
    # Iteriere von 2 bis vorletzten Element
    # ist iterator integer
    for i in range(dist, len(x) - dist):
        # Element davor und danach größer
        if np.all((x[i] > x[slice(i + 1, i + dist)])) and np.all((x[slice(i - dist, i)] < x[i])):
            maxima.append(i)
    return maxima


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


###############################################################################
#
# Adjustments
#
###############################################################################
def adj_mean(data, breaks, axis=0, sample_size=130, borders=30, max_sample=1460, recent=False, ratio=False,
             meanvar=False, **kwargs):
    """ Mean Adjustment of breakpoints

    Args:
        data (array): radiosonde data
        breaks (list): breakpoint indices
        axis (int): axis of datetime
        sample_size (int): minimum sample size
        borders (int): adjust breaks with borders
        max_sample (int): maximum sample size
        recent (bool): use full reference period
        ratio (bool): calculate ratio, instead of difference
        meanvar (bool): mean and var
        **kwargs:

    Returns:
        array : adjusted data
    """
    if not isinstance(data, np.ndarray):
        raise ValueError("requires a numpy array")
    if not isinstance(breaks, (np.ndarray, list)):
        raise ValueError('requires a numpy array')

    data = data.copy()
    dshape = data.shape  # Shape of data (date x levs)
    imax = dshape[axis]  # maximum index
    breaks = np.sort(np.asarray(breaks))  # sort
    breaks = np.append(np.insert(breaks, 0, 0), imax)  # 0 ... ibreaks ... Max
    breaks = breaks.astype(int)
    nb = breaks.size
    # print(breaks)

    for i in range(nb - 2, 0, -1):
        # Indices
        im = breaks[i - 1]  # earlier
        ib = breaks[i]  # current breakpoint
        if recent:
            ip = imax  # Max
        else:
            ip = breaks[i + 1]  # later

        # Slices all axes
        iref = slice(ib, ip)
        isample = slice(im, ib)
        isample = idx2shp(isample, axis, dshape)
        iref = idx2shp(iref, axis, dshape)
        # Before Adjustments
        before = np.nanmean(data[isample], axis=axis)
        # Apply Adjustments
        data[isample] = mean(data[iref], data[isample], axis=axis, sample_size=sample_size,
                             max_sample=max_sample,
                             borders=borders, ratio=ratio, **kwargs)
        #
        # Border zone (Break, Break + borders)
        # Linear interpolation
        #
        if borders > 0:
            zone = slice(ib, ib + borders)
            zone = idx2shp(zone, axis, dshape)
            linzone = conform(np.linspace(0, 1, data[zone].shape[axis]), data[zone].shape)
            if ratio:
                zsample = np.nanmean(data[isample], axis=axis) / np.nanmean(data[iref], axis=axis)
                data[zone] = data[zone] * linzone * zsample
            else:
                zsample = np.nanmean(data[isample], axis=axis) - np.nanmean(data[iref], axis=axis)
                data[zone] = data[zone] + linzone * zsample

    return data


def mean(sample1, sample2, axis=0, sample_size=130, borders=0, max_sample=1460, ratio=True,
         median=False, **kwargs):
    """ Adjustment method using mean differences or ratios

    ratio=False
    dataset[sampleout]  + (MEAN(dataset[sample1]) - MEAN(dataset[sample2]))

    ratio=True
    dataset[sampleout]  * (MEAN(dataset[sample1]) / MEAN(dataset[sample2]))

    Args:
        sample1 (np.ndarray): reference
        sample2 (np.ndarray): sample
        axis (int): date axis
        sample_size (int): minimum sample size
        ratio (bool): use ratio or difference?
        median (bool): use median instead of mean?
        borders (int): around breakpoint
        max_sample (int): maximum sample size

    Returns:
        np.ndarray : mean adjusted dataset
    """
    # minimum sample size, maximum sample size
    if median:
        s1 = nanfunc(sample1,
                     axis=axis,
                     n=sample_size,
                     nmax=max_sample,
                     ffunc=np.nanmedian,
                     borders=borders)
        s2 = nanfunc(sample2,
                     axis=axis,
                     n=sample_size,
                     nmax=max_sample,
                     ffunc=np.nanmedian,
                     borders=borders,
                     flip=True)
    else:
        s1 = nanfunc(sample1,
                     axis=axis,
                     n=sample_size,
                     nmax=max_sample,
                     ffunc=np.nanmean,
                     borders=borders)
        s2 = nanfunc(sample2,
                     axis=axis,
                     n=sample_size,
                     nmax=max_sample,
                     ffunc=np.nanmean,
                     borders=borders,
                     flip=True)

    if ratio:
        # Todo factor amplifies extreme values
        dep = s1 / s2
        dep = np.where(np.isfinite(dep), dep, 1.)  # replace NaN with 1
        sample2 *= dep
    else:
        dep = s1 - s2
        sample2 += dep
    return sample2


def percentile(sample1, sample2, percentiles, axis=0, sample_size=130, borders=0, max_sample=1460, ratio=True,
               apply=None, noise=False, **kwargs):
    """ Adjustment method using percentile differences or ratios

    ratio=False
    dataset[sample1] + ( percentiles(dataset[sample1]) - percentiles(dataset[sample2]) )

    ratio=True
    dataset[sample1] * ( percentiles(dataset[sample1]) / percentiles(dataset[sample2]) )

    Args:
        sample1 (np.ndarray): reference
        sample2 (np.ndarray): sample
        percentiles (list): percentiles to use
        axis (int): date axis
        sample_size (int): minimum sample size
        ratio (bool): use ratio or difference?
        borders (int): around breakpoint
        max_sample (int): maximum sample size

    Returns:
        np.ndarray : percentile adjusted dataset
    """
    # Add 0 and 100, and remove them
    percentiles = np.unique(np.concatenate([[0], percentiles, [100]]))
    percentiles = percentiles[1:-1]  # remove 0 and 100

    # Sample sizes are enough?
    # nsample1 = np.isfinite(dataset[sample1]).sum(axis=axis) > sample_size
    # nsample2 = np.isfinite(dataset[sample2]).sum(axis=axis) > sample_size

    # Percentiles of the samples
    # if special:
    # s1 = np.rollaxis(np.nanpercentile(sample1, percentiles, axis=axis),0, axis)
    # s2 = np.rollaxis(np.nanpercentile(sample2, percentiles, axis=axis),0, axis)
    # print(s1.shape)
    # print(s1[:, 0, 5])
    # print(s2[:, 0, 5])
    # else:
    #
    # Percentiles can be duplicated (because DPD might be integers)
    # limit calculations by sample_size, max_sample, borders
    #
    #     (part A)    |    (part B)
    #               break
    #             >borders<
    #  <max sample         max sample>
    #
    # (part B)
    s1 = nanfunc(sample1,
                 axis=axis,
                 n=sample_size,
                 nmax=max_sample,
                 ffunc=np.nanpercentile,
                 borders=borders,
                 fargs=(percentiles,))
    # flip means that beginning from the back (part A)
    s2 = nanfunc(sample2,
                 axis=axis,
                 n=sample_size,
                 nmax=max_sample,
                 ffunc=np.nanpercentile,
                 borders=borders,
                 fargs=(percentiles,),
                 flip=True)

    # print(s1.shape)
    # print(s1[0, :, 5])
    # print(s2[0, :, 5])

    if ratio:
        dep = np.divide(s1, s2, where=(s2 != 0), out=np.full(s2.shape, 1.))
        dep = np.where(np.isfinite(dep), dep, 1.)  # replace NaN
    else:
        dep = s1 - s2
        dep = np.where(np.isfinite(dep), dep, 0.)

    # Interpolate adjustments
    if apply is None:
        apply = sample2.copy()
    else:
        apply = apply.copy()

    dep = apply_percentile_adjustments(apply, s2, dep, axis=axis, noise=noise)

    if ratio:
        dep = np.where(np.isfinite(dep), dep, 1.)
        apply *= dep
    else:
        dep = np.where(np.isfinite(dep), dep, 0.)
        apply += dep

    return apply


def apply_percentile_adjustments(data, percentiles, adjustment, axis=0, noise=False):
    """ Helper Function for applying percentile adjustments

    Args:
        data (np.ndarray): dataset
        percentiles (np.ndarray): percentiles, points of adjustments
        adjustment (np.ndarray): adjustments to be interpolated
        axis (int): axis of datetime

    Returns:
        np.ndarray : interpolated adjustment, same shape as dataset
    """
    in_dims = list(range(data.ndim))
    # last dim == axis, Last dim should be time/date
    # print(data.shape)
    data = np.transpose(data, in_dims[:axis] + in_dims[axis + 1:] + [axis])
    # print(data.shape)
    percentiles = np.transpose(percentiles, in_dims[:axis] + in_dims[axis + 1:] + [axis])
    adjustment = np.transpose(adjustment, in_dims[:axis] + in_dims[axis + 1:] + [axis])
    # print(percentiles.shape, adjustment.shape)
    adjusts = np.zeros(data.shape)
    # Indices for iteration + expand
    inds = np.ndindex(data.shape[:-1])  # iterate all dimensions but last
    inds = (ind + (Ellipsis,) for ind in inds)  # add last as ':' == Ellipsis == all
    # k = 0
    for ind in inds:
        # INTERP -> Xnew, Xpoints, Fpoints
        iperc, idx = np.unique(percentiles[ind], return_index=True)
        iadj = adjustment[ind][idx]
        if noise:
            adjusts[ind] = np.interp(data[ind] + np.random.normal(size=data[ind].size, scale=0.5), iperc, iadj,
                                     left=np.nan, right=np.nan)
        else:
            adjusts[ind] = np.interp(data[ind], iperc, iadj, left=np.nan,
                                     right=np.nan)
        # if k == 5:
        #     print(ind, data[ind], adjusts[ind], iperc, iadj)
        # k+=1

    # Transform back to original shape
    return np.transpose(adjusts, in_dims[:axis] + in_dims[axis + 1:] + [axis])


"""
Steps:
1. read Netcdf
2. Select only std p levels
3. Make Xarray
4. Run SNHT
5. Detect Breaks
6. Adjust breaks
7. Return results
8. Interpolate to other levels
"""


def main(argv):
    if len(argv) == 0:
        usage(argv[0])
        exit(1)

    # Steps
    #
    # READ NETCDF -> select values + meta information
    # -> LEO routine to read cdm files quickly
    #

    # Use only std pressure levels
    # run test
    # run detector(data, axis=0, dist=365, thres=50, min_levels=3, use_slopes=False, use_first=False, **kwargs):
    # run adj_mean()
    # interpolate
    # write results / metadata (thresholds and so on)


if __name__ == "__main__":
    main(sys.argv)
