#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# This script has been developed in the service contract for C3S
# Calculates radiosonde humidity adjustments based on CDM files
#
# (c) University of Vienna, M. Blaschek, Vienna, Austria
# Copernicus Climate Change Service, 2020
# https://apps.ecmwf.int/datasets/licences/copernicus/
# email michael.blaschek (at) univie.ac.at
# Created: Vienna, 26 August, 2019
# Last Modifed: 15 August, 2020
# Version: 0.1
# -----------------------------------------------------------------------------
__version__ = '0.1'
__author__ = 'MB'
__status__ = 'dev'
__date__ = 'Di 11 Aug 2020 20:27:45 CEST'
__institute__ = 'UNIVIE'
__github__ = 'git@github.com:MBlaschek/CEUAS.git'
__doc__ = """
CDS_adjust Functions v%s
Maintained by %s at %s
Github: %s [%s]
Updated: %s
""" % (__version__, __author__, __institute__, __github__, __status__, __date__)

import logging
import os
import sys
import warnings

import numpy as np
import pandas as pd
import xarray as xr
from numba import njit

try:
    sys.path.append('../../cds-backend/code/')
    import cds_eua3 as eua
except Exception as e:
    print('CDS_EUA3 Module is required. Add to path')
    raise e

warnings.simplefilter("ignore")
np.seterr(invalid='ignore')
logger = logging.getLogger('upperair.adjust')

if not logger.hasHandlers():
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.DEBUG)  # respond only to Warnings and above
    # create formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s | %(funcName)s - %(levelname)s - %(message)s')
    # add formatter to ch
    ch.setFormatter(formatter)
    # add ch to logger
    logger.addHandler(ch)


def logging_set_level(level: int):
    """ Set Logging Level, Default: 10 (DEBUG)"""
    for ihandle in logger.handlers:
        ihandle.setLevel(level)


# in Pa
std_plevs = np.asarray([10, 20, 30, 50, 70, 100, 150, 200, 250, 300, 400, 500, 700, 850, 925, 1000])


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


def write_adjustments(iofile, stest, breaks, adj_depar, variable, interpolate_missing=False):
    if variable is not None:
        # SNHT Test statistics
        iofile.write_observed_data(stest.name,
                                   varnum=eua.cdm_codes[variable],
                                   cube=stest,
                                   interpolate=interpolate_missing,
                                   interpolate_datetime=interpolate_missing,
                                   extrapolate_plevs=interpolate_missing)
        # Breakpoints
        iofile.write_observed_data(breaks.name,
                                   varnum=eua.cdm_codes[variable],
                                   cube=breaks,
                                   interpolate=interpolate_missing,
                                   interpolate_datetime=interpolate_missing,
                                   extrapolate_plevs=interpolate_missing)
        # Adjusted Values
        iofile.write_observed_data(adj_depar.name,
                                   varnum=eua.cdm_codes[variable],
                                   cube=adj_depar,
                                   interpolate=interpolate_missing,
                                   interpolate_datetime=interpolate_missing,
                                   extrapolate_plevs=interpolate_missing)
    else:
        # SNHT Test statistics
        iofile.write_variable(stest.name,
                              cube=stest,
                              interpolate=interpolate_missing,
                              interpolate_datetime=interpolate_missing,
                              extrapolate_plevs=interpolate_missing)
        # Breakpoints
        iofile.write_variable(breaks.name,
                              cube=breaks,
                              interpolate=interpolate_missing,
                              interpolate_datetime=interpolate_missing,
                              extrapolate_plevs=interpolate_missing)
        # Adjusted Values
        iofile.write_variable(adj_depar.name,
                              cube=adj_depar,
                              interpolate=interpolate_missing,
                              interpolate_datetime=interpolate_missing,
                              extrapolate_plevs=interpolate_missing)

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

def adjustment_procedure(obs: xr.DataArray, depar: xr.DataArray, dim: str = 'time', plev: str = 'plev',
                         metadata: bool = False, times: list = [0, 12], span: int = 3, freq: str = '12h',
                         return_dataset: bool = False, mean_adjustments:bool = True, quantile_adjustments:bool=False,
                         **kwargs):
    """

    Args:
        obs:
        depar:
        dim:
        plev:
        metadata:
        times:
        span:
        freq:
        return_dataset:
        mean_adjustments:
        quantile_adjustments:
        **kwargs:

    Returns:

    """
    # Main function for RASO_ADJ_CDM_v0
    # This function executes the homogenization routines and deals with input and output
    #
    # Args:
    #     ifile (str): Filename or string pattern, e.g.: example_data/*[!_out].nc
    #     ofile (str): Output filename or None to use default output naming: [ifile]_out.nc
    #     ta_feature_enabled (bool): experimental flag for temperature adjustments
    #     interpolate_missing (bool): interpolate adjustments to non standard times and pressure levels from the input file?
    #     return_cube (bool): return data cube [hour x time x plev], that is used during the adjustment process
    #     metadata (bool): experimental flag, should point to the metadata for breakpoint identification
    #     thres (int): Threshold value for SNHT, default: 50
    #     window (int): Moving Window for SNHT, default: 1470 (in days, 4 years)
    #     missing (int): Maximum allowed missing values in window, default: 600 (in days)
    #     min_levels (int): Minimum required levels for significant breakpoint, default: 3
    #     dist (int): Minimum distance between breakpoints, default: 730 (in days, 2 years)
    #     sample_size (int): Minimum sample size for statistics, default: 130 (in days)
    #     borders (int): Breakpoint zone, default: 90 (in days)
    #     ratio (int): Use ratio instead of differences, default: 0 (not)
    #     logfile (str): Write messages to a log file
    #     donotwrite(bool): Returns xarray Dataset
    #
    # Returns:
    #     xr.DataArray : adjusted values

    if not isinstance(obs, xr.DataArray):
        raise ValueError('Requires a xarray DataArray, not', type(obs))

    if not isinstance(depar, xr.DataArray):
        raise ValueError('Requires a xarray DataArray, not', type(depar))

    # todo check coordinates of obs and depar (need to be the same)
    assert obs.shape == depar.shape, 'Observation Array and Departures Array do not match in shape?'
    #
    # e.g. 00Z (21Z -1 day to 3Z same day)
    # standard_time
    #
    obs = eua.align_datetime(obs, times=times, span=span, freq=freq, dim=dim, plev=plev)
    icoord = 'standard_%s' % dim
    standard_index = np.where(obs[icoord + '_flag'].values == 1)[0]  # Standard datetime +  selection + times
    # orig_datetime = obs[dim].copy()
    obs = obs.assign_coords({dim: obs[icoord]})
    del obs[icoord]
    del obs[icoord + '_flag']
    # Check if result will be sorted
    reverse_sort = False
    if not obs.isel(**{dim: standard_index})[dim].to_index().is_monotonic:
        logger.warning('Datetime index is not monotonic %s', dim)
        idx = np.argsort(obs.isel(**{dim: standard_index})[dim].values)
        standard_index = standard_index[idx]
        reverse_sort = True
    #
    # Convert to day-night Cube
    #
    xobs = eua.stack_cube_by_hour(obs.isel(**{dim: standard_index}), dim=dim, times=times)
    xdepar = eua.stack_cube_by_hour(depar.isel(**{dim: standard_index}), dim=dim, times=times)
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
    axis = xobs.dims.index(dim)
    name = xobs.name
    sname = xobs.attrs.get('standard_name', name)
    #
    # 3. SNHT
    #
    stest = np.apply_along_axis(test, axis, xdepar.values, window, missing)
    #
    # Attributes
    #
    stest = xr.full_like(xdepar, stest, dtype=stest.dtype)
    stest.name = '{}_snht'.format(name)
    logger.info("SNHT from %s [%s]", xdepar.name, stest.name)
    #
    # Day-Night Departures
    #
    stest += np.apply_along_axis(test, axis,
                                 xobs.sel(hour=12).values - xobs.sel(hour=0),
                                 window, missing)
    stest.attrs.update({'units': '1', 'window': window,
                        'missing': missing, 'standard_name': '{}_snht'.format(sname),
                        'day-night': 'added'})
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
    breaks = xr.full_like(xobs, 0, dtype=np.int)
    breaks.name = '{}_breaks'.format(stest.name)
    breaks.attrs.update({'units': '1', 'dist': dist, 'thres': thres,
                         'min_levels': min_levels, 'standard_name': '{}_breaks'.format(sname)})
    attrs = xobs.attrs.copy()
    attrs.update({'sample_size': kwargs.get('sample_size', 130),
                  'borders': kwargs.get('borders', 90),
                  'ratio': kwargs.get('ratio', 0)})
    for i, ihour in enumerate(times):
        breaks[i, ::] = detector(xobs.values[i, ::], axis - 1, dist=dist, thres=thres, min_levels=min_levels,
                                 **kwargs)
        logger.info("Breakpoints detected %d/%d : %s", i, ihour, breaks.name)
    #
    # 5. Adjust Breakpoints
    #
    if mean_adjustments:
        #
        # Temperature -> MEAN
        #
        adj_depar = xr.full_like(xdepar, 0)
        adj_depar.name = '{}_m'.format(name)
        for i in range(len(times)):
            ibreaks = get_breakpoints(breaks[i, ::], **kwargs)
            logger.info("Breakpoints: %d", len(ibreaks))
            #
            adjv = adjustments(xdepar[i, ::].values, ibreaks,
                               axis=axis - 1,
                               **kwargs)
            # check limits [0 - 1]
            # obs + obs-an-adj - obs-an
            vadj = (xobs[i, ::].values + adjv - xdepar[i, ::].values)
            adjv = np.where((vadj < 0) | (vadj > 1), xdepar[i, ::].values, adjv)
            # new = obs-an-adj + obs - (obs-an)
            adj_depar[i, ::] = (adjv - xdepar[i, ::].values)
        adj_depar.attrs.update(attrs)
        adj_depar.attrs['standard_name'] += '_adjustments'
        adj_depar.attrs['biascor'] = 'mean'

    elif quantile_adjustments:
        #
        # Relative Humidity -> QUANTILE
        #
        adj_depar = xr.full_like(xdepar, 0)
        adj_depar.name = '{}_m'.format(name)
        for i in range(len(times)):
            breaks = get_breakpoints(breaks[i, ::], **kwargs)
            logger.info("Breakpoints: %d", len(breaks))
            #
            adjv = adjustments(xdepar[i, ::].values, breaks,
                                use_mean=False,
                                axis=axis - 1,
                                **kwargs)
            # check limits [0 - 1]
            # obs + obs-an-adj - obs-an
            vadj = (xobs[i, ::].values + adjv - xdepar[i, ::].values)
            adjv = np.where((vadj < 0) | (vadj > 1), xdepar[i, ::].values, adjv)
            # new = obs-an-adj - (obs-an)
            adj_depar[i, ::] = (adjv - xdepar[i, ::].values)
        adj_depar.attrs.update(attrs)
        adj_depar.attrs['standard_name'] += '_adjustments'
        adj_depar.attrs['biascor'] = 'quantile'
    else:
        raise RuntimeError('Either mean_adjustment or quantile_adjustment needs to be set')
    #
    # Convert back to time x plevs
    #
    xadj_depar = eua.unstack_cube_by_hour(adj_depar, dim=dim)
    # fill back (to the original input data)
    if reverse_sort:
        idx = idx.sort()
        standard_index = standard_index[idx]  # undo sorting

    adj_depar = xr.full_like(xdepar, 0)
    adj_depar[:] = xadj_depar
    # use _flag for
    #
    # Return results
    #
    if return_dataset:
        return xr.Dataset([obs, depar, stest, breaks, adj_depar])
    return obs, depar, stest, breaks, adj_depar


def adjustment_procedure_wind(obs_ws, obs_wd, dep_ws, dep_wd, dim:str='time', plev: str = 'plev',
                         metadata: bool = False, times: list = [0, 12], span: int = 3, freq: str = '12h',
                         return_dataset: bool = False, mean_adjustments:bool = True, quantile_adjustments:bool=False,
                         **kwargs):
    """

    Args:
        obs_ws:
        obs_wd:
        dep_ws:
        dep_wd:
        dim:
        plev:
        metadata:
        times:
        span:
        freq:
        return_dataset:
        mean_adjustments:
        quantile_adjustments:
        **kwargs:

    Returns:

    """
    if not isinstance(obs_ws, xr.DataArray):
        raise ValueError('Requires a xarray DataArray, not', type(obs_ws))

    if not isinstance(dep_ws, xr.DataArray):
        raise ValueError('Requires a xarray DataArray, not', type(dep_ws))

    if not isinstance(obs_wd, xr.DataArray):
        raise ValueError('Requires a xarray DataArray, not', type(obs_ws))

    if not isinstance(dep_wd, xr.DataArray):
        raise ValueError('Requires a xarray DataArray, not', type(dep_ws))

    # todo check coordinates of obs and depar (need to be the same)
    assert obs_ws.shape == dep_ws.shape, 'Observation Array and Departures Array do not match in shape?'
    assert obs_wd.shape == dep_wd.shape, 'Observation Array and Departures Array do not match in shape?'
    idim = obs_ws.dims.index(dim)
    #
    # Syncronize wind speed and direction arrays (can have different shapes)
    #
    if obs_ws.shape[idim] != obs_wd.shape[idim]:
        info = "{}<>{}".format(obs_ws.shape, obs_wd.shape)
        obs_ws, obs_wd = xr.broadcast(obs_ws, obs_wd)
        dep_ws, dep_wd = xr.broadcast(dep_ws, dep_wd)
        logger.info('Broadcasting dimensions: %s', info)
    #
    # e.g. 00Z (21Z -1 day to 3Z same day)
    # standard_time
    #
    obs_ws = eua.align_datetime(obs_ws, times=times, span=span, freq=freq, dim=dim, plev=plev)
    icoord = 'standard_%s' % dim
    standard_index = np.where(obs_ws[icoord + '_flag'].values == 1)[0]  # Standard datetime +  selection + times
    # orig_datetime = obs_ws[dim].copy()
    obs_ws = obs_ws.assign_coords({dim: obs_ws[icoord]})
    del obs_ws[icoord]
    del obs_ws[icoord + '_flag']
    # Check if result will be sorted
    reverse_sort = False
    if not obs_ws.isel(**{dim: standard_index})[dim].to_index().is_monotonic:
        logger.warning('Datetime index is not monotonic %s', dim)
        idx = np.argsort(obs_ws.isel(**{dim: standard_index})[dim].values)
        standard_index = standard_index[idx]
        reverse_sort = True
    #
    # Convert to day-night Cube
    #
    xobs_ws = eua.stack_cube_by_hour(obs_ws.isel(**{dim: standard_index}), dim=dim, times=times)
    xobs_wd = eua.stack_cube_by_hour(obs_wd.isel(**{dim: standard_index}), dim=dim, times=times)
    xdep_ws = eua.stack_cube_by_hour(dep_ws.isel(**{dim: standard_index}), dim=dim, times=times)
    xdep_wd = eua.stack_cube_by_hour(dep_wd.isel(**{dim: standard_index}), dim=dim, times=times)
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
    axis = xobs_ws.dims.index(dim)
    name = xobs_ws.name
    # sname = xobs_ws.attrs.get('standard_name', name)
    #
    if True:
        #
        # fix wind direction departures (-180, 180)
        # 5 - 355 = 350 -> 10
        dep_wd = dep_wd.where(dep_wd < 180, dep_wd )
        bg_wd = obs_wd - dep_wd  # Absolute background wind direction
        bg_wd = bg_wd.where(bg_wd < 180, bg_wd - 360.)
        bg_wd = bg_wd.where(bg_wd > -180, 360 + bg_wd)
        bg_ws = obs_ws - dep_ws  # Absolute background wind speed
    #
    # 3. SNHT
    #
    stest = np.apply_along_axis(test, axis, xdep_ws.values, window, missing)
    stest += np.apply_along_axis(test, axis, xdep_wd.values, window, missing)
    #
    # Attributes
    #
    stest = xr.full_like(xdep_ws, stest, dtype=stest.dtype)
    stest.name = '{}_snht'.format(name)
    logger.info("SNHT from %s [%s]", xdep_ws.name, stest.name)
    #
    # Day-Night Departures
    #
    if 0 in times and 12 in times:
        stest += np.apply_along_axis(test, axis,
                                     xobs_ws.sel(hour=12).values - xobs_ws.sel(hour=0),
                                     window, missing)
    stest.attrs.update({'units': '1', 'window': window,
                        'missing': missing, 'standard_name': '{}_snht'.format(sname),
                        'day-night': 'added'})


# -----------------------------------------------------------------------------
#
# Script entry point
#
# -----------------------------------------------------------------------------


if __name__ == "__main__":
    import sys
    import getopt
    import glob

    # handle arguments
    kwargs = {'verbose': 1}
    known = ["help", "backend", "frontend", "outdir", "temperature", "humidity", "winds", "dates", "plevs", "feedback", "fbgroup"]
    ifile = None
    odir = None
    dates = None
    plevs = None
    feedback = None
    feedback_group = 'era5fb'
    interpolate_missing = False
    do_temperature = False
    do_humidity = False
    do_winds = False

    try:
        known, knames = _cmd_arguments(sys.argv[1:], known)
        opts, args = getopt.getopt(sys.argv[1:], "b:f:ho:v:d:p:", known)

    except getopt.GetoptError as err:
        usage()
        raise err

    for opt, arg in opts:
        if opt in ("-f", "--frontend"):
            input_frontend = True
            ifile = arg
        
        elif opt in ("-b", "--backend"):
            input_frontend = False
            ifile = arg

        elif opt in ("-o", "--outdir"):
            odir = arg

        elif opt in ("temperature"):
            do_temperature = True

        elif opt in ("humidity"):
            do_humidity = True

        elif opt in ("winds"):    
            do_winds = True

        elif opt in ("-d", "--dates"):
            if ',' in arg:
                arg = arg.split(',')
            dates = arg  # --dates date,date

        elif opt in ("-p", "--plevs"):
            if ',' in arg:
                arg = arg.split(',')
            else:
                arg = [arg]
            plevs = arg

        elif opt in ("--feedback"):
            feedback = arg
        
        elif opt in ("--fbgroup"):
            feedback_group = arg  

        elif opt in ("-h", "--help"):
            usage()
            sys.exit(0)

        elif opt == "--interpolate_missing":
            interpolate_missing = True

        elif any([opt[2:] in i for i in knames]):
            if arg != '':
                kwargs[opt[2:]] = eval(arg)
            else:
                kwargs[opt[2:]] = True

        else:
            assert False, "unhandled option"
    #
    # Check input
    #
    if ifile is None:
        usage()
        raise RuntimeError("Missing input file (-f, -b)", ifile)
    
    if not (do_winds or do_temperature or do_humidity):
        raise RuntimeError('Please specify at least one option: --temperature, --humidty or --winds')

    if plevs is None:
        plevs = std_plevs * 100
    else:
        plevs = list(map(int, plevs))

    if input_frontend:
        #
        # FRONTEND File
        #
        if '*' in ifile:
            ifile = glob.glob(ifile)
        elif ',' in ifile:
            ifile = ifile.split(',')
        else:
            pass
        if not isinstance(ifile, list):
            ifile = [ifile]

        for i in ifile:
            assert os.path.isfile(i), i

        if feedback is None:
            feedback = 'obs_minus_an'

        if do_temperature or do_humidity:
            iofile = eua.CDMDataset(ifile)
            #
            # Temperature or Humdity adjustment
            #
            variable = 'ta' if do_temperature else 'hur'
            data = iofile.read_data_to_cube(variable,
                                            dates=dates,
                                            plevs=plevs)[variable]
            departures = iofile.read_data_to_cube(variable,
                                                  dates=dates,
                                                  plevs=plevs,
                                                  feedback=feedback)[variable]
            _, _, stest, breaks, adj_depar = adjustment_procedure(data, departures,
                                                                  metadata=False,
                                                                  times=[0, 12],
                                                                  dim='time',
                                                                  plev='plev',
                                                                  return_dataset=False)
            #
            # Write back adjusted (interpolation, extrapolation)
            #
            write_adjustments(iofile, stest, breaks, adj_depar, None,
                              interpolate_missing=interpolate_missing)
        else:
            #
            # Wind speed and wind direction adjustment
            #
            assert len(ifile) ==2, 'Inputfiles: wind_direction and wind_speed are required'
            filepool = {}
            for i in ifile:
                tmp = eua.CDMDataset(i)
                if 'wind_from_direction' in tmp.groups:
                    variable = 'wind_from_direction'
                    wd_data = tmp.read_data_to_cube(variable,
                                                    dates=dates,
                                                    plevs=plevs)[variable]
                    wd_departures = tmp.read_data_to_cube(variable,
                                                          dates=dates,
                                                          plevs=plevs,
                                                          feedback=feedback)[variable]
                    filepool[variable] = tmp
                else:
                    variable = 'wind_speed'
                    ws_data = tmp.read_data_to_cube(variable,
                                                    dates=dates,
                                                    plevs=plevs)[variable]
                    ws_departures = tmp.read_data_to_cube(variable,
                                                          dates=dates,
                                                          plevs=plevs,
                                                          feedback=feedback)[variable]
                    filepool[variable] = tmp
            #
            _, _, stest, breaks, adj_depar_wd, adj_depar_ws = adjustment_procedure_wind()
            # wind direction
            write_adjustments(filepool['wind_from_direction'], stest, breaks, adj_depar_wd, None,
                              interpolate_missing=interpolate_missing)
            # wind speed
            write_adjustments(filepool['wind_speed'], stest, breaks, adj_depar_ws, None,
                              interpolate_missing=interpolate_missing)
    else:
        #
        # BACKEND File
        #
        iofile = eua.CDMDataset(ifile)
        if feedback is None:
            feedback = 'an_depar@body'

        if not feedback_group in iofile.groups:
            raise ValueError('Feedback Group', feedback_group, 'not found')
        if not feedback in iofile[feedback_group].keys():
            raise ValueError('Feedback ',feedback,' not in Feedback Group ',feedback_group)
        
        if odir is not None:
            odir = "{}/{}".format(odir,os.path.basename(ifile))

        iofile.reopen(write_to_filename=odir, mode='r+')
        if do_temperature:
            # Code: 85
            variable = 'temperature'
            #
            # Write back adjusted (interpolation, extrapolation)
            #
            write_adjustments(iofile, stest, breaks, adj_depar, variable, interpolate_missing=interpolate_missing)

        if do_humidity:
            # Code 38 , relative humidity
            variable = 'relative_humidity'
            #
            # Write back adjusted (interpolation, extrapolation)
            #
            write_adjustments(iofile, stest, breaks, adj_depar, variable, interpolate_missing=interpolate_missing)
        
        if do_winds:
            # Code 106 (wind_direction), 107 (wind_speed)
            variable = 'wind_direction'

            #
            # Write back adjusted (interpolation, extrapolation)
            #
            write_adjustments(iofile, stest, breaks, adj_depar, variable, interpolate_missing=interpolate_missing)
            variable = 'wind_speed'
            write_adjustments(iofile, stest, breaks, adj_depar, variable, interpolate_missing=interpolate_missing)
 
    # FIN
