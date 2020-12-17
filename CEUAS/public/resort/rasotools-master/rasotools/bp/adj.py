# -*- coding: utf-8 -*-
import numpy as np

from . import dep
from ..fun import message, check_kw
from ..fun.cal import nancount
from ..met.convert import _conform as conform

np.seterr(invalid='ignore')


__all__ = ['mean', 'percentile', 'percentile_reference']


def mean(data, breaks, axis=0, sample_size=130, borders=30, max_sample=1460, recent=False, ratio=False,
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
        # print(i, isample, iref)

        # add seasonal sampling check ?
        # if the two samples have different seasonal sampling -> an adjustment bias might occur ?
        # if the seasonal cycle is badly represented this poses a problem

        # adjust data
        if meanvar:
            data[isample] = dep.meanvar(data[iref], data[isample], axis=axis, sample_size=sample_size,
                                        max_sample=max_sample,
                                        borders=borders, **kwargs)
        else:
            data[isample] = dep.mean(data[iref], data[isample], axis=axis, sample_size=sample_size,
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

        # Debug infos
        if check_kw('verbose', value=2, **kwargs):
            sdata = stats(data, iref, isample, axis=axis, a=before)
            sdata = np.array_str(sdata, precision=2, suppress_small=True)
            message(sdata, **kwargs)
    return data


def percentile(data, breaks, axis=0, percentilen=None, sample_size=130, borders=30, max_sample=1460,
               recent=False, ratio=False, **kwargs):
    """ Percentile Adjustment of breakpoints

    Args:
        data (array): radiosonde data
        breaks (list): breakpoint indices
        axis (int): axis of datetime dimension
        percentilen (list): percentile ranges
        sample_size (int): minimum sample size
        borders (int): adjust breaks with borders
        max_sample (int): maximum sample size
        recent (bool): use full reference period
        ratio (bool): calculate ratio, instead of difference
        **kwargs:

    Returns:
        array : adjusted data
    """
    if not isinstance(data, np.ndarray):
        raise ValueError("requires a numpy array")

    if not isinstance(breaks, (list, np.ndarray)):
        raise ValueError('requires a numpy array of list')

    data = data.copy()
    #
    # Percentiles
    #
    if percentilen is None:
        percentilen = np.arange(0, 101, 10)

    nq = len(percentilen)
    #
    # Sample size per Percentile
    #
    sample_size = sample_size // nq
    if sample_size < 3:
        sample_size = 3

    message('Sample size:', sample_size, 'N-Q:', nq, **kwargs)

    dshape = data.shape
    imax = dshape[axis]
    breaks = np.sort(np.asarray(breaks))  # sort
    breaks = np.append(np.insert(breaks, 0, 0), imax)  # 0 ... ibreaks ... None
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

        # slices all axes
        iref = slice(ib, ip)
        isample = slice(im, ib)
        isample = idx2shp(isample, axis, dshape)
        iref = idx2shp(iref, axis, dshape)

        before = np.nanmean(data[isample], axis=axis)
        # Apply Adjustments
        data[isample] = dep.percentile(data[iref], data[isample], percentilen, axis=axis, sample_size=sample_size,
                                       max_sample=max_sample, borders=borders, ratio=ratio)
        # Debug infos
        if check_kw('verbose', value=2, **kwargs):
            sdata = stats(data, iref, isample, axis=axis, a=before)
            sdata = np.array_str(sdata, precision=2, suppress_small=True)
            message(sdata, **kwargs)

    return data


def percentile_reference_period(xdata, ydata, breaks, axis=0, percentilen=None, sample_size=130, ratio=False,
                                ref_period=None, **kwargs):
    """ Adjust a dataset based on a reference period and percentile matching

    Args:
        xdata (np.ndarray): Reference dataset
        ydata (np.ndarray): Adjust dataset
        breaks (list): Breakpoints if no reference period is refered to
        axis (int): datetime dimension axis
        percentilen (list): percentiles
        sample_size (int): minimum sample size
        ratio (bool):
        ref_period (slice): indices for reference period
        **kwargs:

    Returns:
        np.ndarray : ydata adjusted
    """
    if not isinstance(xdata, np.ndarray):
        raise ValueError("requires a numpy array")
    if not isinstance(ydata, np.ndarray):
        raise ValueError("requires a numpy array")
    if not isinstance(breaks, (list, np.ndarray)):
        raise ValueError('requires a numpy array of list')

    xdata = xdata.copy()
    ydata = ydata.copy()
    if percentilen is None:
        percentilen = np.arange(0, 101, 10)

    nq = len(percentilen)
    sample_size = sample_size // nq
    if sample_size < 3:
        sample_size = 3

    message('Sample size:', sample_size, 'N-Q:', nq, **kwargs)

    dshape = xdata.shape
    imax = dshape[axis]
    breaks = np.sort(np.asarray(breaks))  # sort
    breaks = np.append(np.insert(breaks, 0, 0), imax)  # 0 ... ibreaks ... None
    breaks = breaks.astype(int)
    #
    # 1. Adjust Reference to match distribution of unbiased period
    #
    isample = slice(0, imax)  # Everything
    if ref_period is None:
        iref = slice(breaks[-2], imax)  # Reference before 1st Breakpoint
    else:
        iref = ref_period
    isample = idx2shp(isample, axis, dshape)
    iref = idx2shp(iref, axis, dshape)
    #
    # Apply Dist. from xdata[iref] to all ydata  (Match dists.)
    #
    ydata[isample] = dep.percentile(xdata[iref], ydata[iref], percentilen,
                                    axis=axis,
                                    sample_size=sample_size,
                                    max_sample=np.nan,
                                    ratio=ratio,
                                    apply=ydata[isample])

    return ydata


def percentile_reference(xdata, ydata, breaks, axis=0, percentilen=None, sample_size=130, borders=30, max_sample=1460,
                         recent=False, ratio=False, **kwargs):
    # xdata is RASO
    # ydata is Reference: ERA, CERA, JRA
    if not isinstance(xdata, np.ndarray):
        raise ValueError("requires a numpy array")
    if not isinstance(ydata, np.ndarray):
        raise ValueError("requires a numpy array")
    if not isinstance(breaks, (list, np.ndarray)):
        raise ValueError('requires a numpy array of list')

    xdata = xdata.copy()
    ydata = ydata.copy()
    if percentilen is None:
        percentilen = np.arange(0, 101, 10)

    nq = len(percentilen)
    sample_size = sample_size // nq
    if sample_size < 3:
        sample_size = 3

    message('Sample size:', sample_size, 'N-Q:', nq, **kwargs)

    dshape = xdata.shape
    imax = dshape[axis]
    breaks = np.sort(np.asarray(breaks))  # sort
    breaks = np.append(np.insert(breaks, 0, 0), imax)  # 0 ... ibreaks ... None
    breaks = breaks.astype(int)
    nb = breaks.size

    #
    # 1. Loop Breakpoints and adjust backwards using ydata as reference
    # An Absolute ?
    #
    for i in range(nb - 2, 0, -1):
        # Indices
        im = breaks[i - 1]  # earlier
        ib = breaks[i]  # current breakpoint
        if recent:
            ip = imax  # Max
        else:
            ip = breaks[i + 1]  # later

        # slices all axes
        iref = slice(ib, ip)  # B part
        isample = slice(im, ib)  # A part
        isample = idx2shp(isample, axis, dshape)
        iref = idx2shp(iref, axis, dshape)

        before = np.nanmean(xdata[isample], axis=axis)
        #
        # Use same sample for both data
        # compare xdata to ydata (REF)
        #
        xdata[isample] = dep.percentile(ydata[isample], xdata[isample], percentilen,
                                        axis=axis,
                                        sample_size=sample_size,
                                        borders=borders,
                                        max_sample=max_sample,
                                        ratio=ratio,
                                        **kwargs
                                        )
        # Debug infos
        if check_kw('verbose', value=2, **kwargs):
            sdata = stats(xdata, iref, isample, axis=axis, a=before)
            sdata = np.array_str(sdata, precision=2, suppress_small=True)
            if i == (nb - 2):
                sdata = '[    i     #S          B      dB       S      dS       R   #R   ]\n' + sdata
            message(sdata, **kwargs)

    return xdata


def idx2shp(idx, axis, shape):
    index = [slice(None)] * len(shape)
    index[axis] = idx
    return tuple(index)


def stats(data, ref, sample, axis=0, a=None):
    # print counts, means, dep
    sn = nancount(data[sample], axis=axis)
    s = np.nanmean(data[sample], axis=axis)
    rn = nancount(data[ref], axis=axis)
    r = np.nanmean(data[ref], axis=axis)
    n = np.squeeze(np.arange(0, np.size(sn)))
    if a is not None:
        # counts, before, sample, diff, ref, counts
        return np.array([n, sn, a, r - a, s, r - s, r, rn]).T
    return np.array([n, sn, s, r - s, r, rn]).T

# def index_samples(breaks, ibreak, axis, dshape, recent=False, sample_size=130, borders=30, max_sample=1460):
#     """ Apply Breakpoints to data shape, return index lists
#
#     Args:
#         breaks (list, np.ndarray): Breakpoints
#         ibreak (int): index of current Breakpoint
#         axis (int): datetime axis
#         dshape (tuple): shape of data array
#         recent (bool): always use the whole recent part
#         sample_size (int): minimum sample size
#         borders (int): Breakpoint borders/uncertainty
#         max_sample: maximum samples
#
#     Returns:
#         list, list, list : Biased , Sample, Reference
#     """
#     imax = dshape[axis]
#     biased, ref = sample_indices(breaks, ibreak, imax, recent=recent)
#     sample, ref = adjust_samples(biased, ref, sample_size, borders=borders, max_sample=max_sample)
#     return idx2shp(biased, axis, dshape), idx2shp(sample, axis, dshape), idx2shp(ref, axis, dshape)
# def sample_indices(breaks, ibreak, imax, recent=False):
#     """ Indices of Samples before and after a bp
#
#     Args:
#         breaks (list, np.ndarray) :        Breakpoints
#         ibreak (int) :          current Breakpoint
#         imax (int) :            maximum index
#         recent (bool) :         use all newest Data
#
#     Returns:
#         tuple: sample indices
#     """
#     n = len(breaks)
#
#     if ibreak > 0:
#         anfang = breaks[ibreak - 1]
#     else:
#         anfang = 0  # takes all the stuff? or only sometime after the break?
#
#     mitte = breaks[ibreak]  # Mittelpunkt ist der Bruchpunkt
#
#     if ibreak == (n - 1) or recent:
#         ende = imax  # most recent
#     else:
#         ende = breaks[ibreak + 1]  # bp before
#
#     sample1 = slice(anfang, mitte)  # erste Teil (indices niedriger)
#     sample2 = slice(mitte, ende)  # Zweite Teil (indices höher)
#     return sample1, sample2

#
# def adjust_samples(ibiased, iref, sample_size, borders, max_sample):
#     # start -> kleiner index (nächster Bruchpunkt, früher)
#     # stop -> grosser index (Bruchpunkt)
#     n = ibiased.stop - ibiased.start  # sample size
#     isample = slice(ibiased.start, ibiased.stop)  # isample == ibiased
#     if n - 2 * borders > sample_size:
#         isample = slice(ibiased.start + borders, ibiased.stop - borders)  # ohne Borders
#         if n - 2 * borders > max_sample:
#             isample = slice(ibiased.stop - borders - max_sample, ibiased.stop - borders)  # nur max_sample
#
#     n = iref.stop - iref.start
#     if n - 2 * borders > sample_size:
#         iref = slice(iref.start + borders, iref.stop - borders)
#         if n - 2 * borders > max_sample:
#             iref = slice(iref.start, iref.start + max_sample)
#
#     return isample, iref
#
# def breakpoint_zone(x, thres, axis=0, k=200, recent=False, target=None):
#     n = x.shape[axis]
#
#     if target is not None:
#         m1 = target
#     else:
#         m1 = np.nanmean(x, axis=axis)  # total median as backup
#
#     dep = x - m1
#     j = None
#
#     for i in range(n - k, 0, -k):
#         if not recent:
#             j = i + k if i + k < n else n
#         itx = idx2shp(slice(i, j), axis, x.shape)
#         s1 = np.where(np.isfinite(dep[itx]).sum(axis=axis) > 0, np.nanmean(dep[itx], axis=axis), m1)
#         l = i - k if i - k > 0 else 0
#         jtx = idx2shp(slice(l, i), axis, x.shape)
#         s2 = np.where(np.isfinite(dep[jtx]).sum(axis=axis) > 0, np.nanmean(dep[jtx], axis=axis), m1)
#         dep[jtx] += np.where(np.abs(s1 - s2) > thres, (s1 - s2), 0.)
#     return m1 + dep
