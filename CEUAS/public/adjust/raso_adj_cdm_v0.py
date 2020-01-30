#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# This script has been developed in the service contract for C3S
# Calculates radiosonde humidity adjustments based on CDM files
#
# (c) University of Vienna, M. Blaschek, Vienna, Austria
# Released under GNU Public License (GPL)
# email michael.blaschek (at) univie.ac.at
# -----------------------------------------------------------------------------

__version__ = '0.1'
__author__ = 'MB'
__status__ = 'delivered'
__date__ = 'Mit Jan 15 15:30:59 CET 2020'
__institute__ = 'Univie, IMGW'
__github__ = 'git@github.com:MBlaschek/CEUAS.git'
__doc__ = """
Radiosonde Homogenisation Software v%s
Maintained by %s at %s
Github: %s [%s]
License: C3S
Updated: %s
""" % (__version__, __author__, __institute__, __github__, __status__, __date__)

import os
import sys

import numpy as np
import pandas as pd
import xarray as xr
from numba import njit

np.seterr(invalid='ignore')

# in Pa
std_plevels = np.array(
    [1000., 2000., 3000., 5000., 7000., 10000., 15000., 20000., 25000., 30000., 40000., 50000., 70000.,
     85000., 92500., 100000.])


def usage():
    print("""
Run standardized radiosonde homogenisation software on CDM compliant file

{} -h -f [file] -o [name] 

Options:
    -h              Help
    --help      
    -f []           Input CDM compliant file
    --file []       
    -o []           Output name
    --output []
    
Optional Keyword Options:
    --thres []          Threshold value for SNHT, default: 50
    --window []         Moving Window for SNHT, default: 1470 (in days, 4 years)
    --missing []        Maximum allowed missing values in window, default: 600 (in days)
    --min_levels []     Minimum required levels for significant breakpoint, default: 3
    --dist []           Minimum distance between breakpoints, default: 730 (in days, 2 years)
    --sample_size []    Minimum sample size for statistics, default: 130 (in days)
    --borders []        Breakpoint zone, default: 90 (in days)
    --ratio []          Use ratio instead of differences, default: 0 (not)

    --logfile []        Write messages to a log file

Experimental Keyword Options:
    --donotwrite          Returns xarray Dataset
    --enable_ta_feature   Apply Temperature adjustments
    --interpolate_missing Interpolate Adjustments to non-standard times and pressure levels
    
    """.format(__file__.split('/')[-1]))


# -----------------------------------------------------------------------------
#
# Helper functions
#
# -----------------------------------------------------------------------------

def now(timespec='auto'):
    """ Datetime string
    Returns:
        str : datetime now
    """
    import datetime
    return datetime.datetime.now().isoformat(timespec=timespec)


def _print_string(*args, mname=None, adddate=False, **kwargs):
    text = " ".join([str(i) for i in args])
    if mname is not None:
        text = "[{}] {}".format(mname, text)
    if adddate:
        text = "[" + now() + "] " + text
    return text


def message(*args, verbose=0, level=0, logfile=None, **kwargs):
    """ Message function
    Args:
        *args:  text to be printed
        verbose (int): level of verbosness
        level (int): level of visibility (verbose > level: printed)
        logfile (str): logfile
        **kwargs:
    Returns:
        str : message
    """
    if logfile is not None:
        # with open(kwargs['filename'], 'a' if not kwargs.get('force', False) else 'w') as f:
        with open(logfile, 'a') as f:
            f.write(_print_string(*args, **kwargs) + "\n")

    elif verbose > level:
        text = _print_string(*args, **kwargs)
        print(text)
    else:
        pass


def update_kw(name, value, **kwargs):
    """ Update keyword dictionary on the fly
    """
    kwargs.update({name: value})
    return kwargs


def conform(data, shape):
    """ Make numpy array conform to a certain shape
    Args:
        data (np.ndarray): input data
        shape (tuple, list): desired shape

    Returns:
        np.ndarray : reshaped
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
    """ Count values excluding NaN
    Args:
        x (np.ndarray): input dataset
        axis (int): axis
        keepdims (bool): keep dimensions
    Returns:
        np.ndarray : sum of values
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
        fargs (tuple): function arguments

    Returns:
        np.ndarray : func of values at axis, with sample size, borders and maximum
    """
    if ffunc is None:
        ffunc = np.nanmean
    return np.apply_along_axis(sample, axis, data, n, nmax, ffunc, borders=borders, flip=flip, fargs=fargs)


def sample(values, nmin, nmax, func, borders=0, flip=False, fargs=(), **kwargs):
    """ Apply a function (func) to a sample of defined size

    Args:
        values (np.ndarray): input values
        nmin (int): minimum required values
        nmax (int): maximum number of values
        func (callable): function to execute
        borders (int): number of values to skip (start-borders , end-borders)
        flip (bool): reverse order of array
        fargs (tuple): arguments to function func
        **kwargs:

    Returns:
        np.ndarray : according to func
    """
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


def table_to_dataset(data, dim='time', plev='plev', levels=None, **kwargs):
    """ Convert pandas Dataframe to xarray Dataset

    Args:
        data (pd.DataFrame): input dataframe (columns are variables)
        dim (str): datetime dimension
        plev (str): pressure dimension
        levels (list): pressure levels to consider
        **kwargs:

    Returns:
        xr.Dataset : 2d (datetime x pressure levels) x variables
    """
    from xarray import Dataset
    if levels is None:
        levels = [1000., 2000., 3000., 5000., 7000., 10000., 15000., 20000., 25000., 30000., 40000., 50000., 70000.,
                  85000., 92500., 100000.]
    # dimensions for output
    varis = [dim, plev]
    attrs = None
    if isinstance(data, Dataset):
        # copy attributes
        attrs = data.attrs.copy()
        tatts = data[dim].attrs
        vatt = {i: data[i].attrs.copy() for i in data.data_vars}
        #
        # to pandas dataframe
        #
        data = data.to_dataframe()
        data.index.name = dim

    #
    # select only valid levels
    #
    message("Selecting only standard pressure levels", **kwargs)
    data = data[data[plev].isin(levels)]
    #
    # convert to xarray
    #
    data = data.reset_index().set_index(varis)
    if not data.index.is_unique:
        data = data.loc[~data.index.duplicated()]  # remove duplicated

    data = data.to_xarray()  # 1D -> 2D
    if attrs is not None:
        # add attributes again
        for i, j in vatt.items():
            data[i].attrs.update(j)
        data.attrs.update(attrs)
        data[dim].attrs.update(tatts)

    return data


@np.vectorize
def fix_datetime(itime, span=6, debug=False):
    """ Fix datetime to standard datetime with hour precision

    Args:
        itime (datetime): Datetime
        span (int): allowed difference to standard datetime (0,6,12,18)

    Returns:
        datetime : standard datetime
    """
    import pandas as pd
    itime = pd.Timestamp(itime)  # (time: 34%)
    # span=6 -> 0, 12
    # [18, 6[ , [6, 18[
    # span=3 -> 0, 6, 12, 18
    # [21, 3[, [3,9[, [9,15[, [15,21[
    for ihour in range(0, 24, span * 2):
        # 0 - 6 + 24 = 18
        lower = (ihour - span + 24) % 24
        # 0 + 6 + 24 = 6
        upper = (ihour + span + 24) % 24
        # 18 >= 18 or 18 < 6  > 00
        # 0 >= 18 or 0 < 6    > 00
        if debug:
            print("%d [%d] %d >= %d < %d" % (ihour, span, lower, itime.hour, upper))

        if (ihour - span) < 0:
            if itime.hour >= lower or itime.hour < upper:
                rx = itime.replace(hour=ihour, minute=0, second=0, microsecond=0)
                if itime.hour >= (24 - span):
                    rx = rx + pd.DateOffset(days=1)
                return rx.to_datetime64()
        else:
            if lower <= itime.hour < upper:
                rx = itime.replace(hour=ihour, minute=0, second=0, microsecond=0)
                if itime.hour >= (24 - span):
                    rx = rx + pd.DateOffset(days=1)
                return rx.to_datetime64()


# -----------------------------------------------------------------------------
#
# Detection of Breakpoints
#
# -----------------------------------------------------------------------------
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
    """ Generate indices

    Args:
        idx (int): index
        axis (int): axis
        shape (tuple): shape

    Returns:
        tuple : selection slice
    """
    index = [slice(None)] * len(shape)
    index[axis] = idx
    return tuple(index)


@njit
def local_maxima(x, dist=365):
    """ Find local maxima, using Numba

    Args:
        x (np.ndarray): input
        dist (int): minimum distance between peaks

    Returns:
        list : list of maxima
    """
    maxima = []  # Leere Liste
    # Iteriere von 2 bis vorletzten Element
    # ist iterator integer
    for i in range(dist, len(x) - dist):
        # Element davor und danach größer
        if np.all((x[i] > x[slice(i + 1, i + dist)])) and np.all((x[slice(i - dist, i)] < x[i])):
            maxima.append(i)
    return maxima


def test(x, window, missing):
    """Standard Normal Homogeneity Test (SNHT) with a running window
    Wrapper function for numba_snhtmov

    Args:
        x (np.ndarray) : input data
        window (int) : window size (in days)
        missing (int) : allowed missing values (in days)
    Returns:
        np.ndarray : SNHT
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


def get_breakpoints(data, value=2, dim='time', return_startstop=False, startstop_min=0, **kwargs):
    """ Return breakpoints from breakpoint data

    Args:

        data (xr.DataArray): input DataArray
        value (int): breakpoint indicator value
        dim (str): datetime dim
        startstop_min (int):
        return_startstop (bool):
        **kwargs:

    Returns:
        list : breakpoints
    """
    if not isinstance(data, xr.DataArray):
        raise ValueError("Require a DataArray / Dataset object", type(data))

    if dim not in data.dims:
        raise ValueError("Requires a datetime dimension", data.dims)

    if len(data.dims) > 2:
        RuntimeWarning("More than two dimensions found: ", str(data.dims))

    #
    # Dimension of time
    #
    axis = data.dims.index(dim)
    #
    # Search Threshold
    #
    tmp = np.where(data.values >= value)
    i = list(map(int, np.unique(tmp[axis])))
    dates = np.datetime_as_string(data[dim].values, unit='D')
    e = []
    s = []
    #
    # multi-dimension / combine to only time axis
    #
    if data.ndim > 1:
        summe = data.values.sum(axis=1 if axis == 0 else 0)
    else:
        summe = data.values

    for k in i:
        l = np.where(summe[:k][::-1] <= startstop_min)[0][0]
        m = np.where(summe[k:] <= startstop_min)[0][0]
        e += [k - l]
        s += [k + m]

    if len(i) > 0:
        message("Breakpoints for ", data.name, **kwargs)
        message("[%8s] [%8s] [%8s] [%8s] [ #]" % ('idx', 'end', 'peak', 'start'), **kwargs)
        for j, k, l in zip(i, s, e):
            message("[%8s] %s %s %s %4d" % (j, dates[l], dates[j], dates[k], k - l), **kwargs)

    if return_startstop:
        return i, e, s
    return i


# -----------------------------------------------------------------------------
#
# Adjustments
#
# -----------------------------------------------------------------------------

def adjustments(data, breaks, use_mean=True, axis=0, sample_size=130, borders=30, max_sample=1460, recent=False,
                ratio=False, **kwargs):
    """ Adjustment of breakpoints

    Args:
        data (np.ndarray): radiosonde data
        breaks (list): breakpoint indices
        use_mean (bool): mean or quantile adjustments
        axis (int): axis of datetime
        sample_size (int): minimum sample size
        borders (int): adjust breaks with borders
        max_sample (int): maximum sample size
        recent (bool): use full reference period
        ratio (bool): calculate ratio, instead of difference
        **kwargs:

    Returns:
        np.ndarray : adjusted data
    """
    if not isinstance(data, np.ndarray):
        raise ValueError("requires a numpy array")
    if not isinstance(breaks, (np.ndarray, list)):
        raise ValueError('requires a numpy array')

    data = data.copy()
    if not use_mean:
        percentilen = kwargs.get('percentiles', [0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100])
        nq = len(percentilen)
        #
        # Sample size per Percentile
        #
        sample_size = sample_size // nq
        if sample_size < 3:
            sample_size = 3
        message('Sample size:', sample_size, 'N-Q:', nq, **kwargs)

    dshape = data.shape  # Shape of data (date x levs)
    imax = dshape[axis]  # maximum index
    breaks = np.sort(np.asarray(breaks))  # sort
    breaks = np.append(np.insert(breaks, 0, 0), imax)  # 0 ... ibreaks ... Max
    breaks = breaks.astype(int)
    nb = breaks.size

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
        if use_mean:
            data[isample] = mean(data[iref], data[isample],
                                 axis=axis,
                                 sample_size=sample_size,
                                 max_sample=max_sample,
                                 borders=borders,
                                 ratio=ratio, **kwargs)
        else:
            data[isample] = percentile(data[iref], data[isample], percentilen,
                                       axis=axis,
                                       sample_size=sample_size,
                                       max_sample=max_sample,
                                       borders=borders,
                                       ratio=ratio, **kwargs)

        #
        # Border zone (Break, Break + borders)
        # Linear interpolation
        #
        if use_mean and borders > 0:
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
        apply (np.ndarray): apply adjustments to this array
        noise (bool): add random noise to adjustments

    Returns:
        np.ndarray : percentile adjusted dataset
    """
    # Add 0 and 100, and remove them
    percentiles = np.unique(np.concatenate([[0], percentiles, [100]]))
    percentiles = percentiles[1:-1]  # remove 0 and 100
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
        noise (bool): add random noise to adjustments

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


# -----------------------------------------------------------------------------
#
# Interpolation
#
# -----------------------------------------------------------------------------
def level_interpolation(idata, dim='time', method='linear', fill_value=None, extrapolate=True, **kwargs):
    """ Interpolate CDM variable per dimension (time, plev)

    Args:
        idata (xr.DataArray): Inputdata
        dim (str): coordinate of input, e.g. time, plev
        method (str): interpolation method, e.g. linear, log-linear for plev
        fill_value (any): Interpolation fill_value: 'extrapolate', None, np.nan
        extrapolate (bool): fill missing values [True]
        **kwargs:

    Returns:
        xr.DataArray : Interpolate data, same as input

    Examples:
        Interpolate adjustments per level backwards and forward in time
        >>> data['hur_q'].groupby('plev').apply(level_interpolation)

        Interpolate adjustments per profile up and down in pressure using a log-linear interpolation
        >>> data['hur_q'].groupby('time').apply(level_interpolation, dim='plev', method='log-linear')

        Interpolate relative humidity per profile up and down in pressure using a log-linear interpolation, but no
        extrapolation
        >>> data['hur'].groupby('time').apply(level_interpolation, dim='plev', method='log-linear', extrapolate=False)

    """
    if not isinstance(idata, xr.DataArray):
        raise ValueError("Requires a DataArray, not ", type(idata))
    #
    # maybe check if dimensions are present
    #
    if idata.isnull().all().item():
        return idata

    idim = idata.dims[0]
    obs = idata[idim].values.copy()
    #
    # swap to interpolation dimension
    #
    idata = idata.swap_dims({idim: dim})
    #
    # Find duplicated values /not allowed in interpolation
    #
    itx = idata[dim].to_index()
    ittx = itx.duplicated(keep=False) & ~np.isfinite(idata.values)  # BUG duplicates in dim & not a value -> skip
    ittx = np.where(~np.isfinite(itx), True, ittx)  # BUG missing values can be in the coordinate (e.g. plev)
    #
    # some duplicates might remain
    #
    if itx[~ittx].duplicated().any():
        ittx[~ittx] = itx[~ittx].duplicated()  # duplicates with same values / should not happen

    message(itx.size, ittx.sum(), mname='DUPLICATES', **update_kw('level', 1, **kwargs))
    if method == 'log-linear':
        #
        # Pressure levels to log
        #
        idata.values[~ittx] = idata[~ittx] \
            .assign_coords({'plev': np.log(idata['plev'].values[~ittx])}) \
            .interpolate_na(dim, method='linear', fill_value=fill_value) \
            .values
    else:
        idata[~ittx] = idata[~ittx].interpolate_na(dim, method=method, fill_value=fill_value)
    #
    # Extrapolate in dimension (backward and forward) or (up and down)
    #
    if extrapolate:
        idata = idata.bfill(dim).ffill(dim)
    #
    # revert to initial dimension
    #
    idata[idim] = (dim, obs)  # add coordinate / dimension
    return idata.swap_dims({dim: idim}).drop('obs')  # swap dimension back, remove obs coordinate


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def _cmd_arguments(args, longs):
    add = longs[:]  # copy
    names = []
    for i, iarg in enumerate(args):
        if iarg[0] == '-':
            if iarg[1] == '-':
                if iarg[2:] not in longs:
                    if len(args) > i + 1:
                        jarg = args[i + 1]
                        if jarg[0] != '-':
                            add.append(iarg[2:] + '=')
                        else:
                            add.append(iarg[2:])
                    else:
                        add.append(iarg[2:])
                    names.append(iarg[2:])

    return add, names


# -----------------------------------------------------------------------------
#
# Main function that does all the steps in adjusting
#
# -----------------------------------------------------------------------------

def main(ifile=None, ofile=None, ta_feature_enabled=False, interpolate_missing=False, return_cube=False, metadata=False,
         **kwargs):
    """ Main function for RASO_ADJ_CDM_v0
    This function executes the homogenization routines and deals with input and output

    Args:
        ifile (str): Filename or string pattern, e.g.: example_data/*[!_out].nc
        ofile (str): Output filename or None to use default output naming: [ifile]_out.nc
        ta_feature_enabled (bool): experimental flag for temperature adjustments
        interpolate_missing (bool): interpolate adjustments to non standard times and pressure levels from the input file?
        return_cube (bool): return data cube [hour x time x plev], that is used during the adjustment process
        metadata (bool): experimental flag, should point to the metadata for breakpoint identification
        thres (int): Threshold value for SNHT, default: 50
        window (int): Moving Window for SNHT, default: 1470 (in days, 4 years)
        missing (int): Maximum allowed missing values in window, default: 600 (in days)
        min_levels (int): Minimum required levels for significant breakpoint, default: 3
        dist (int): Minimum distance between breakpoints, default: 730 (in days, 2 years)
        sample_size (int): Minimum sample size for statistics, default: 130 (in days)
        borders (int): Breakpoint zone, default: 90 (in days)
        ratio (int): Use ratio instead of differences, default: 0 (not)
        logfile (str): Write messages to a log file
        donotwrite(bool): Returns xarray Dataset

    Returns:
        None : written to output file [ofile]
        Dataset :

    Examples:
        Use all input files in example_data directory:
        >>> main(ifile="example_data/*[!_out].nc")

        Return the output Dataset to the current session:
        >>> data = main(ifile="example_data/*[!_out].nc", donotwrite=True)
    """

    kwargs.update({'verbose': 1})

    if ifile is None:
        import getopt

        known = ["help", "file", "output"]
        ifile = None
        ofile = None

        try:
            known, knames = _cmd_arguments(sys.argv[1:], known)
            opts, args = getopt.getopt(sys.argv[1:], "f:ho:", known)

        except getopt.GetoptError as err:
            usage()
            message(str(err), mname='ERROR', verbose=1)
            return 2

        for opt, arg in opts:

            if opt in ("-f", "--file"):
                ifile = arg

            elif opt in ("-o", "--output"):
                ofile = arg

            elif opt in ("-h", "--help"):
                usage()
                return 0

            elif opt == "--enable_ta_feature":
                ta_feature_enabled = True

            elif opt == "--interpolate_missing":
                interpolate_missing = True

            elif any([opt[2:] in i for i in knames]):
                if arg != '':
                    kwargs[opt[2:]] = eval(arg)
                else:
                    kwargs[opt[2:]] = True

            else:
                assert False, "unhandled option"

    if ifile is None:
        usage()
        message("Missing input file", ifile, mname='ERROR', verbose=1)
        return 1

    variables = []
    #
    # 1. Step (Read data from NetCDF)
    #
    if "*" in ifile:
        import glob
        ifiles = glob.glob(ifile)
        message("Multiple input files: ", ifile, ' | #', len(ifiles), mname='INFO', **kwargs)
        data = []
        traj_data = []
        ifile = {}
        for jfile in ifiles:
            idata = xr.load_dataset(jfile)
            ivar = [k for k in list(idata.data_vars) if '_' not in k][0]
            message(jfile, ivar, mname='INPUT', **kwargs)
            idata = idata.rename({i: '{}_{}'.format(ivar, i) for i in list(idata.data_vars) if i != ivar})
            #
            # Trajectory data
            #
            to_be_dropped = []
            for jvar in idata.data_vars:
                if 'trajectory' in jvar:
                    to_be_dropped.append(jvar)

            if len(to_be_dropped) > 0:
                traj_data.append(idata[to_be_dropped])
                idata = idata.drop(to_be_dropped)

            data.append(idata)
            variables.append(ivar)
            ifile[ivar] = jfile
        #
        # to larger Dataset
        #
        data = xr.merge(data)
        traj_data = xr.merge(traj_data)
    else:
        if not os.path.isfile(ifile):
            usage()
            message("Missing input file", ifile, mname='ERROR', verbose=1)
            return 1

        traj_data = xr.Dataset()
        data = xr.load_dataset(ifile)
        ivar = [k for k in list(data.data_vars) if '_' not in k][0]
        message(ifile, ivar, mname='INPUT', **kwargs)
        data = data.rename({i: '{}_{}'.format(ivar, i) for i in list(data.data_vars) if i != ivar})
        #
        # Trajectory data
        #
        to_be_dropped = []
        for jvar in data.data_vars:
            if 'trajectory' in jvar:
                to_be_dropped.append(jvar)

        if len(to_be_dropped) > 0:
            traj_data = data[to_be_dropped].copy()
            data = data.drop(to_be_dropped)
        variables.append(ivar)
        ifile = {ivar: ifile}
    #
    # 2. Step (Convert to DataCube, rename dimension to time)
    #
    message("Converting to DataCube ...", mname='CONVERT', **kwargs)
    #
    # Make sure we have the original dataset
    #
    raw_data = data.copy()  # Backup for later
    dim = 'time'
    data = data.swap_dims({'obs': 'time'})
    #
    # Add reverse index to convert back to obs
    #
    data['obs'] = ('time', np.arange(0, data.time.size))
    message("Observations", ",".join(["{}: {}".format(i, j) for i, j in raw_data.dims.items()]), mname='BKP', **kwargs)
    message("Trajectories", ",".join(["{}: {}".format(i, j) for i, j in traj_data.dims.items()]), mname='BKP', **kwargs)
    message("Dates:", data[dim].min().dt.strftime("%Y-%m-%d %H:%M:%S").item(), ' - ',
            data[dim].max().dt.strftime("%Y-%m-%d %H:%M:%S").item(), "Duplicates:",
            data[dim].to_index().duplicated().sum(),
            mname=dim.upper(),
            **kwargs)
    #
    # Select only standard pressure data
    #
    data = table_to_dataset(data, mname='CONVERT', **kwargs)
    message("Done", mname='CONVERT', **kwargs)
    #
    # Add Standard time variable (needed for detection parameters)
    # todo add functionality to use hourly sondes, adjust detection parameters
    #   based on data availability (resample to hour freq, lots of missing values?)
    #
    # e.g. 00Z (21Z -1 day to 3Z same day)
    #
    dates = data[dim].to_index()
    newdates = fix_datetime(dates, span=3)  # use plus minus 3 hours
    #
    # find duplicates
    #
    u, c = np.unique(newdates, return_counts=True)
    conflicts = u[c > 1]
    for iconf in conflicts:
        # get duplicates
        indices = np.where(newdates == iconf)[0]
        #
        # check which is closer to standard time
        #
        offset = np.abs((dates[indices] - iconf) / np.timedelta64(1, 'h'))
        j = np.argsort(offset)  # sort time offsets (first we want)
        for m, k in enumerate(offset[j]):
            if m == 0:
                continue  # this is the minimum

            # change back the others or add a delay to remove duplicates
            newdates[indices[j][m]] = dates[indices[j][m]]  # revert back
    #
    # Add Standard time as variable
    #
    data['time_orig'] = ('time', dates)
    message("Standard time calculated, Duplicates resolved:", conflicts.size, mname='CONVERT', **kwargs)
    data = data.assign_coords({dim: newdates})  # overwrite dim with new dates
    #
    # Day-Night Split / selection
    #
    times = (0, 12)
    data = data.sel(**{dim: data[dim].dt.hour.isin(times) & (data[dim].dt.minute == 0)})
    data = dict(data.groupby(dim + '.hour'))
    for ikey in data.keys():
        data[ikey] = data[ikey].assign_coords(
            {dim: data[ikey][dim].to_index().to_period('D').to_timestamp().values})

    data = xr.concat(data.values(), dim=pd.Index(data.keys(), name='hour'))
    # make sure the shape is as promissed:
    data = data.reindex({'hour': list(times)})
    message("Converting to day-night Array [hour x time x pressure]", mname='CONVERT', **kwargs)
    #
    # What variables are present?
    # Feedback information is present ?
    #
    analysis = {}
    firstg = {}
    bias = {}
    for ivar in variables:
        if "{}_obs_minus_an".format(ivar) in data.data_vars:
            message("Departures found", ivar, "obs_minus_an", mname='CHECK', **kwargs)
            analysis[ivar] = "{}_obs_minus_an".format(ivar)

        if "{}_obs_minus_fg".format(ivar) in data.data_vars:
            message("Departures found", ivar, "obs_minus_fg", mname='CHECK', **kwargs)
            firstg[ivar] = "{}_obs_minus_fg".format(ivar)

        if "{}_bias_estimate".format(ivar) in data.data_vars:
            message("Departures found", ivar, "bias_estimate", mname='CHECK', **kwargs)
            bias[ivar] = "{}_bias_estimate".format(ivar)
    #
    # Apply Bias_estimates to Departures ?
    #
    # Default behaviour is to apply the bias to get absolute values
    #

    #
    # Analysis Departures are required
    #
    for ivar in variables:
        if ivar not in analysis.keys():
            message("Missing Analysis departures for", ivar, "\n Homogenization requires Analysis departures.",
                    mname='ERROR', verbose=1)
            return 1
    #
    # Some Variables for Breakpoint detection
    #
    window = kwargs.get('window', 1460)  # means 4 years on daily basis
    missing = kwargs.get('missing', 600)
    thres = kwargs.get('thres', 50)
    min_levels = kwargs.get('min_levels', 3)
    dist = kwargs.get('dist', 730)  # two years on daily basis
    #
    #
    #
    for ivar in variables:
        dim = 'time'
        jvar = analysis[ivar]
        dims = list(data[jvar].dims)
        axis = dims.index(dim)
        #
        # 3. SNHT
        #
        stest = np.apply_along_axis(test, axis, data[jvar].values, window, missing)
        #
        # Attributes
        #
        svar = '{}_snht'.format(jvar)
        data[svar] = (dims, stest)
        message("Test statistics calculted", svar, mname='SNHT', **kwargs)
        #
        # Day-Night Departures
        #
        stest = np.apply_along_axis(test, axis,
                                    data[ivar].sel(hour=12).values - data[ivar].sel(hour=0),
                                    window, missing)
        # Add Day-Night Departures
        data[svar] = data[svar] + stest
        attrs = {'units': '1', 'window': window, 'missing': missing, 'standard_name': svar, 'day-night': 'added'}
        data[svar].attrs.update(attrs)
        #
        # Check for Metadata in CDM
        #
        if metadata:
            #
            # In a future version Metadata can be read and used here.
            #
            pass

        #
        # 4. Detect Breakpoints
        #
        breaks = np.full_like(data[svar].values, 0)
        for i, ihour in enumerate(data.hour.values):
            breaks[i, ::] = detector(data[svar].values[i, ::], axis - 1, dist=dist, thres=thres, min_levels=min_levels,
                                     **kwargs)
        #
        bvar = '{}_breaks'.format(svar)
        attrs = {'units': '1', 'dist': dist, 'thres': thres, 'min_levels': min_levels, 'standard_name': bvar}
        data[bvar] = (dims, breaks)
        data[bvar].attrs.update(attrs)
        message("Test statistics calculated", bvar, mname='DETECT', **kwargs)
        attrs = data[ivar].attrs.copy()
        attrs.update({'sample_size': kwargs.get('sample_size', 130),
                      'borders': kwargs.get('borders', 90),
                      'ratio': kwargs.get('ratio', 0)})
        #
        # 5. Adjust Breakpoints
        #
        if ivar in ['ta'] and ta_feature_enabled:
            #
            # Temperature -> MEAN
            #
            avar = '{}_m'.format(ivar)
            data[avar] = xr.full_like(data[ivar], 0)
            for i, idata in data.groupby('hour'):
                breaks = get_breakpoints(idata[bvar], **kwargs)
                message("Breakpoints: ", len(breaks), mname='ADJUST', **kwargs)
                #
                adjv = adjustments(idata[jvar].values, breaks,
                                   axis=axis - 1,
                                   mname='ADJUST', **kwargs)
                # check limits [0 - 1]
                # obs + obs-an-adj - obs-an
                vadj = (idata[ivar].values + adjv - idata[jvar].values)
                adjv = np.where((vadj < 0) | (vadj > 1), idata[jvar].values, adjv)
                # new = obs-an-adj + obs - (obs-an)
                data[avar].loc[{'hour': i}] = (adjv - idata[jvar].values)
            data[avar].attrs.update(attrs)
            data[avar].attrs['standard_name'] += '_adjustments'
            data[avar].attrs['biascor'] = 'mean'

        if ivar in ['hur']:
            #
            # Relative Humidity -> QUANTILE
            #
            avar = '{}_q'.format(ivar)
            data[avar] = xr.full_like(data[ivar], 0)
            for i, idata in data.groupby('hour'):
                breaks = get_breakpoints(idata[bvar], **kwargs)
                message("Breakpoints: ", len(breaks), mname='ADJUST', **kwargs)
                #
                adjv = adjustments(idata[jvar].values, breaks,
                                   use_mean=False,
                                   axis=axis - 1,
                                   mname='ADJUST', **kwargs)
                # check limits [0 - 1]
                # obs + obs-an-adj - obs-an
                vadj = (idata[ivar].values + adjv - idata[jvar].values)
                adjv = np.where((vadj < 0) | (vadj > 1), idata[jvar].values, adjv)
                # new = obs-an-adj - (obs-an)
                data[avar].loc[{'hour': i}] = (adjv - idata[jvar].values)
            data[avar].attrs.update(attrs)
            data[avar].attrs['standard_name'] += '_adjustments'
            data[avar].attrs['biascor'] = 'quantile'

    if return_cube:
        return data

    #
    # Restore original time [hour x date] => datetime
    #
    message("hour x time -> datetime ", mname='ALIGN', **kwargs)
    data = dict(data.groupby('hour'))
    for ikey in data.keys():
        # use the backup of datetimes 'time_orig'
        data[ikey] = data[ikey].assign_coords({dim: data[ikey]['time_orig'].to_index()})

    data = xr.concat(data.values(), dim=dim).sortby(dim)
    data = data.reset_coords('hour').rename({'hour': 'standard_hour'})
    data = data.sel(time=np.isfinite(data.time))
    #
    # Merge Variables back into raw_data
    #
    message("Selecting variables", mname='REVERSE', **kwargs)
    outvars = ['obs']
    for ivar in variables:
        if ivar + '_m' in data.data_vars:
            outvars.append(ivar + '_m')
            outvars.append(ivar + '_obs_minus_an_snht')
            outvars.append(ivar + '_obs_minus_an_snht_breaks')
        if ivar + '_q' in data.data_vars:
            outvars.append(ivar + '_q')
            outvars.append(ivar + '_obs_minus_an_snht')
            outvars.append(ivar + '_obs_minus_an_snht_breaks')
    #
    # Convert back to obs
    #
    message("time to obs", mname='REVERSE', **kwargs)
    #
    # Processing to DataFrame and back to xarray
    #
    attrs = {i: data[i].attrs.copy() for i in outvars}
    data = data[outvars] \
        .to_dataframe() \
        .reset_index() \
        .dropna(subset=['obs']) \
        .set_index('obs') \
        .to_xarray()
    # put metadata back
    for i in outvars:
        data[i].attrs.update(attrs[i])
    # convert to integer
    data['obs'] = data.obs.astype(int)
    data = data.set_coords(['time', 'plev'])
    raw_data['obs'] = ('obs', raw_data.obs.values)
    raw_data = raw_data.merge(data, join='left')
    raw_data = raw_data.drop('obs')
    #
    # Fill SNHT, breakpoints with zeros
    #
    for ivar in list(raw_data.data_vars):
        if 'snht' in ivar or 'breaks' in ivar:
            raw_data[ivar] = raw_data[ivar].fillna(0)
    #
    # Interpolate adjustments back to intermediate times & significant levels
    #
    if interpolate_missing:
        for ivar in variables:
            if ivar + '_m' in raw_data.data_vars:
                message("Interpolating non-standard times ", ivar + '_m', mname='INTERP', **kwargs)
                #
                # Interpolate back/forward in time
                # extrapolate=True
                #
                subset = raw_data['plev'].isin(std_plevels)  # only standard pressure levels
                raw_data[ivar + '_m'][subset] = raw_data[ivar + '_m'][subset] \
                    .groupby('plev').apply(level_interpolation, dim=dim, **kwargs)

                message("Interpolating non-standard pressure ", ivar + '_m', mname='INTERP', **kwargs)
                #
                # Interpolate up/down levels
                #
                raw_data[ivar + '_m'] = raw_data[ivar + '_m'] \
                    .groupby('time').apply(level_interpolation, dim='plev',
                                           **update_kw('method', 'log-linear', **kwargs))

            if ivar + '_q' in raw_data.data_vars:
                message("Interpolating non-standard times ", ivar + '_q', mname='INTERP', **kwargs)
                #
                # Interpolate back/forward in time
                # extrapolate=True
                #
                subset = raw_data['plev'].isin(std_plevels)  # only standard pressure levels
                raw_data[ivar + '_q'][subset] = raw_data[ivar + '_q'][subset] \
                    .groupby('plev').apply(level_interpolation, dim=dim, **kwargs)
                message("Interpolating non-standard pressure ", ivar + '_q', mname='INTERP', **kwargs)
                #
                # Interpolate up/down levels
                #
                raw_data[ivar + '_q'] = raw_data[ivar + '_q'] \
                    .groupby('time').apply(level_interpolation, dim='plev',
                                           **update_kw('method', 'log-linear', **kwargs))

        message("Interpolation complete", mname='INTERP', **kwargs)
    #
    # Add Trajectory data again
    #
    message('Add trajectory information', mname='TRAJ', **kwargs)
    raw_data = raw_data.merge(traj_data, join='left')
    #
    # Store information in file
    #
    if kwargs.get('donotwrite', False):
        message("Returning data", mname='OUT', **kwargs)
        return raw_data
    #
    # Use Xarray to netcdf library
    #
    message("Writing ...", mname='OUT', **kwargs)
    if ofile is None:
        ofile = {ivar: ifile[ivar].replace('.nc', '_out.nc') for ivar in variables}

    for ivar in variables:
        if ivar == 'ta' and not ta_feature_enabled:
            continue

        selection = [jvar for jvar in list(raw_data.data_vars) if ivar + '_' in jvar] + [ivar]
        if len(selection) == 0:
            continue

        tmp = raw_data[selection]
        tmp = tmp.rename({jvar: jvar.replace(ivar + '_', '') for jvar in list(tmp.data_vars) if
                          jvar not in [ivar, ivar + '_q', ivar + '_m']})
        if False:
            # Alternative Package with improved performance and data compression
            # h5netcdf
            tmp.to_netcdf(ofile[ivar],
                          mode='w',
                          format='netCDF4',
                          engine=kwargs.get('engine', 'h5netcdf'),
                          encoding={i: {'compression': 'gzip', 'compression_opts': 9} for i in
                                    list(tmp.data_vars)})

        message(ofile[ivar], ", ".join(tmp.data_vars), mname='OUTPUT', **kwargs)
        tmp.to_netcdf(ofile[ivar],
                      mode='w',
                      format='netCDF4',
                      engine=kwargs.get('engine', 'netcdf4'),
                      encoding={i: {'zlib': True, 'complevel': 9} for i in
                                list(tmp.data_vars)})

    message("Finished", mname='MAIN', **kwargs)


# -----------------------------------------------------------------------------
#
# Script entry point
#
# -----------------------------------------------------------------------------


if __name__ == "__main__":
    message("Executing ...", mname='MAIN', verbose=1)
    sys.exit(main())
