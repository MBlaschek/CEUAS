# -*- coding: utf-8 -*-
import numpy as np

from ..fun.cal import nanfunc

__all__ = ['mean', 'percentile']


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


def meanvar(sample1, sample2, axis=0, sample_size=130, borders=0, max_sample=1460, **kwargs):
    """ Adjustment method using mean differences or ratios

    dataset[sampleout]  + (MEAN(dataset[sample1]) - MEAN(dataset[sample2]))

    Args:
        sample1 (np.ndarray): reference
        sample2 (np.ndarray): sample
        axis (int): date axis
        sample_size (int): minimum sample size
        borders (int): around breakpoint
        max_sample (int): maximum sample size

    Returns:
        np.ndarray : mean adjusted dataset
    """
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
    s1v = nanfunc(sample1,
                  axis=axis,
                  n=sample_size,
                  nmax=max_sample,
                  ffunc=np.nanvar,
                  borders=borders)
    s2v = nanfunc(sample2,
                  axis=axis,
                  n=sample_size,
                  nmax=max_sample,
                  ffunc=np.nanvar,
                  borders=borders,
                  flip=True)

    # MEAN
    dep = s1 - s2
    # VAR
    fac = np.divide(s1v, s2v, out=np.ones(s2v.shape), where=s2v != 0)
    sample2 += (dep * fac)
    return sample2


def percentile(sample1, sample2, percentiles, axis=0, sample_size=130, borders=0, max_sample=1460, ratio=True,
               apply=None, noise=False,**kwargs):
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
    #s1 = np.rollaxis(np.nanpercentile(sample1, percentiles, axis=axis),0, axis)
    #s2 = np.rollaxis(np.nanpercentile(sample2, percentiles, axis=axis),0, axis)
    #print(s1.shape)
    #print(s1[:, 0, 5])
    #print(s2[:, 0, 5])
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


#
# Helper functions
#


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
            adjusts[ind] = np.interp(data[ind] + np.random.normal(size=data[ind].size, scale=0.5), iperc, iadj, left=np.nan, right=np.nan)
        else:
            adjusts[ind] = np.interp(data[ind], iperc, iadj, left=np.nan,
                                     right=np.nan)
        # if k == 5:
        #     print(ind, data[ind], adjusts[ind], iperc, iadj)
        # k+=1

    # Transform back to original shape
    return np.transpose(adjusts, in_dims[:axis] + in_dims[axis + 1:] + [axis])
