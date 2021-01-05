# -*- coding: utf-8 -*-
import numpy as np
from pandas import Index
from xarray import Dataset, DataArray, set_options, concat

from .. import fun as ff

__all__ = ['apply_threshold', 'snht',
           'adjust_mean', 'adjust_percentiles', 'adjust_percentiles_ref', 'adjust_reference_period',
           'breakpoint_statistics', 'get_breakpoints']


def snht(data, dim='time', var=None, dep=None, suffix=None, window=1460, missing=600, **kwargs):
    """ Calculate a Standard Normal Homogeinity Test

    Args:
        data (DataArray, Dataset):
        dim (str): datetime dimension
        var (str): variable if Dataset
        dep (str, DataArray): departure variable
        suffix (str): add to name of new variables
        window (int): running window (timesteps)
        missing (int): allowed missing values in window

    Returns:
        Dataset : test statistics
        or
        DataArray
    """
    from .det import test

    if not isinstance(data, (DataArray, Dataset)):
        raise ValueError('Requires an xarray DataArray or Dataset', type(data))

    if isinstance(data, DataArray):
        idata = data.copy()
        var = ff.suche123(idata.name, var, 'var')
    else:
        ivars = list(data.data_vars)
        if len(ivars) == 1:
            var = ivars[0]
        elif var is None:
            raise ValueError("Dataset requires a var")
        else:
            pass
        idata = data[var].copy()

    if dim not in idata.dims:
        raise ValueError('requires a datetime dimension', dim)

    if suffix is not None:
        if suffix[0] != '_':
            suffix = '_' + suffix
            raise Warning('suffix needs an _. Added:', suffix)
    else:
        suffix = ''

    axis = idata.dims.index(dim)
    attrs = idata.attrs.copy()
    dep_add = False

    if dep is not None:
        if isinstance(dep, str) and isinstance(data, Dataset):
            dep = data[dep]

        elif isinstance(dep, DataArray):
            dep = dep
            dep_add = True

        else:
            raise ValueError("dep var not present")

        #
        with set_options(keep_attrs=True):
            idata = (idata - dep.reindex_like(idata))

        ff.xarray.set_attrs(idata, 'standard_name', add='_departure', default='departure')

    stest = np.apply_along_axis(test, axis, idata.values, window, missing)
    attrs.update({'units': '1', 'window': window, 'missing': missing})

    if isinstance(data, DataArray):
        data = data.to_dataset(name=var)

    if dep is not None:
        data[var + '_dep' + suffix] = idata
        if dep_add:
            data[dep.name if dep.name is not None else 'dep'] = dep
        attrs['cell_method'] = 'departure ' + dep.name + attrs.get('cell_method', '')

    data[var + '_snht' + suffix] = (list(idata.dims), stest)
    ff.xarray.set_attrs(data[var + '_snht' + suffix], 'standard_name', add='_snht', default='snht')
    data[var + '_snht' + suffix].attrs.update(attrs)
    return data


def apply_threshold(data, dim='time', var=None, name='breaks', suffix=None, thres=50, dist=730, min_levels=3,
                    ensemble=False, **kwargs):
    """ Apply threshold on SNHT to detect breakpoints

    Args:
        data (DataArray, Dataset):
        dim (str): datetime dimension
        var (str): variable if Dataset
        name (str): name of new variable with above threshold (breaks)
        suffix (str): add to name of new variables
        thres (int, float): threshold value
        dist (int): distance between breaks
        min_levels (int): minimum significant levels for breaks
        ensemble (bool): run ensemble on thresholds, nthres=50,
    Returns:
        Dataset
    """
    from xarray import DataArray, Dataset
    from .det import detector, detector_ensemble

    if not isinstance(data, (DataArray, Dataset)):
        raise ValueError('Requires an xarray DataArray or Dataset', type(data))

    if suffix is not None:
        if suffix[0] != '_':
            suffix = '_' + suffix
            raise Warning('suffix needs an _. Added:', suffix)
    else:
        suffix = ''

    if isinstance(data, DataArray):
        idata = data.copy()  # link
        var = idata.name if idata.name is not None else 'var'
    else:
        if var is None or var not in list(data.data_vars):
            raise ValueError('Requires a variable name: var=', list(data.data_vars))

        idata = data[var]  # link

    if idata.ndim > 2:
        raise ValueError("Maximum of 2 dimensions: ", idata.shape)

    if dim not in idata.dims:
        raise ValueError('requires a datetime dimension', dim)

    axis = idata.dims.index(dim)
    params = {'units': '1', 'thres': thres, 'dist': dist, 'min_levels': min_levels, 'standard_name': 'break_flag',
              'flag_valus': [0, 1, 2, 3], 'valid_range': (0, 3),
              'flag_meanings': 'not_significant significant significant_at_other_level significant_at_level'}

    if ensemble:
        kwargs['nthres'] = kwargs.get('nthres', 50)
        breaks = detector_ensemble(idata.values, axis=axis, **kwargs)
        params['thres'] = 'ens%d' % kwargs.get('nthres')

    else:
        breaks = detector(idata.values, axis=axis, dist=dist, thres=thres, min_levels=min_levels,
                          **ff.levelup(**kwargs))

    name = var + '_' + name + suffix
    if isinstance(data, DataArray):
        data = idata.to_dataset(name=var)

    data[name] = (list(idata.dims), breaks)
    data[name].attrs.update(params)
    return data


def get_breakpoints(data, value=2, dim='time', return_startstop=False, startstop_min=0, **kwargs):
    """ Return breakpoints

    Args:

        data (DataArray): input dataset
        value (int): breakpoint indicator value
        dim (str): datetime dim
        startstop_min (int):
        return_startstop (bool):
        **kwargs:

    Returns:
        list : breakpoints
    """
    if not isinstance(data, DataArray):
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

    if kwargs.get('verbose', 0) > 0:
        if len(i) > 0:
            ff.message("Breakpoints for ", data.name, **kwargs)
            ff.message("[%8s] [%8s] [%8s] [%8s] [ #]" % ('idx', 'end', 'peak', 'start'), **kwargs)
            ff.message("\n".join(
                ["[%8s] %s %s %s %4d" % (j, dates[l], dates[j], dates[k], k - l) for j, k, l in zip(i, s, e)]),
                **kwargs)

    if return_startstop:
        return i, e, s
    return i


def adjust_mean(data, name, breakname, dim='time', suffix='_m', ratio=False, **kwargs):
    """Detect and Correct Radiosonde biases from Departure Statistics
Use ERA-Interim departures to detect breakpoints and
adjustments these with a mean  adjustment going back in time.


    Args:
        data (Dataset): Input Dataset with different variables
        name (str): Name of variable to adjust
        breakname (str): Name of variable with breakpoint information
        dim (str): datetime dimension
        suffix (str): add to name of new variables
        ratio (bool): differences or ratios?

    Optional Args:
        sample_size (int):  minimum Sample size [130]
        borders (int):  biased sample before and after a break [None]
        recent (bool):  Use all recent Data to adjustments
        ratio (bool):  Use ratio instead of differences

    Returns:
        Dataset
    """
    from . import adj
    if not isinstance(data, Dataset):
        raise ValueError("Requires a Dataset object", type(data))

    if not isinstance(name, str):
        raise ValueError("Requires a string name", type(name))

    if name not in data.data_vars:
        raise ValueError("dataset var not present")

    if breakname not in data.data_vars:
        raise ValueError("requires a breaks dataset var")
    #
    # Copy data
    #
    idata = data[name].copy()
    values = idata.values
    #
    # get Breakpoints
    #
    breaks = get_breakpoints(data[breakname], dim=dim, **ff.levelup(**kwargs))
    axis = idata.dims.index(dim)
    params = idata.attrs.copy()  # deprecated (xr-patch)
    ff.message(name, str(values.shape), 'A:', axis, **kwargs)
    params.update(
        {'sample_size': kwargs.get('sample_size', 130), 'borders': kwargs.get('borders', 90), 'ratio': int(ratio)})
    ff.message(ff.dict2str(params), **ff.levelup(**kwargs))
    stdn = data[name].attrs.get('standard_name', name)

    data[name + suffix] = (idata.dims, adj.mean(values, breaks, axis=axis, ratio=ratio, **kwargs))
    data[name + suffix].attrs.update(params)
    data[name + suffix].attrs['biascor'] = 'mean'
    if 'niter' in data[name + suffix].attrs:
        data[name + suffix].attrs['niter'] += 1
    else:
        data[name + suffix].attrs['niter'] = 1
        data[name + suffix].attrs['standard_name'] = stdn + '_mean_adj'

    return data


def adjust_percentiles(data, name, breakname, dim='time', dep_var=None, suffix='_q', percentilen=None, **kwargs):
    """Detect and Correct Radiosonde biases from Departure Statistics
Use ERA-Interim departures to detect breakpoints and
adjustments these with a mean and a percentile adjustment going back in time.


    Args:
        data (Dataset): Input Dataset with different variables
        name (str): Name of variable to adjust
        breakname (str): Name of variable with breakpoint information
        dim (str): datetime dimension
        dep_var (str): Name of variable to use as a departure
        suffix (str): add to name of new variables
        percentilen (list): percentiles for percentile_cor

    Optional Args:
        sample_size (int):  minimum Sample size [130]
        borders (int):  biased sample before and after a break [None]
        bounded (tuple):  limit correction to bounds
        recent (bool):  Use all recent Data to adjustments
        ratio (bool):  Use ratio instead of differences

    Returns:
        Dataset
    """
    from . import adj
    if not isinstance(data, Dataset):
        raise ValueError("Requires a Dataset object", type(data))

    if not isinstance(name, str):
        raise ValueError("Requires a string name", type(name))

    if name not in data.data_vars:
        raise ValueError("dataset var not present")

    if breakname not in data.data_vars:
        raise ValueError("requires a breaks dataset var")

    idata = data[name].copy()

    if dep_var is not None:
        if dep_var not in data.data_vars:
            raise ValueError("dep var not present", data.data_vars)

        with set_options(keep_attrs=True):
            idata = (idata - data[dep_var].reindex_like(idata))

    if percentilen is None:
        percentilen = np.arange(0, 101, 10)

    values = idata.values
    breaks = get_breakpoints(data[breakname], dim=dim, **ff.levelup(**kwargs))
    axis = idata.dims.index(dim)
    params = idata.attrs.copy()  # deprecated (xr-patch)

    ff.message(name, str(values.shape), 'A:', axis, 'Q:', np.size(percentilen), "Dep:", str(dep_var), **kwargs)

    params.update({'sample_size': kwargs.get('sample_size', 130), 'borders': kwargs.get('borders', 90)})

    ff.message(ff.dict2str(params), **ff.levelup(**kwargs))
    stdn = data[name].attrs.get('standard_name', name)

    data[name + suffix] = (
        idata.dims, adj.percentile(values, breaks, axis=axis, percentilen=percentilen, **kwargs))
    data[name + suffix].attrs.update(params)
    data[name + suffix].attrs['biascor'] = 'percentil'
    data[name + suffix].attrs['standard_name'] = stdn + '_percentil_adj'

    return data


def adjust_percentiles_ref(data, name, adjname, breakname, dim='time', suffix='_qa', percentilen=None,
                           adjust_reference=True, **kwargs):
    """Detect and Correct Radiosonde biases from Departure Statistics
Use ERA-Interim departures to detect breakpoints and
adjustments these with a mean and a percentile adjustment going back in time.


    Args:
        data (Dataset): Input Dataset with different variables
        name (str): Name of variable to adjust
        adjname (str): Name of adjust variable
        breakname (str): Name of variable with breakpoint information
        dim (str): datetime dimension
        suffix (str): add to name of new variables
        percentilen (list): percentilen
        adjust_reference (bool): return adjusted reference?

    Optional Args:
        sample_size (int):  minimum Sample size [130]
        borders (int):  biased sample before and after a break [None]
        recent (bool):  Use all recent Data to adjustments
        ratio (bool):  Use ratio instead of differences
        ref_period (slice): period to use for quantile matching of reference

    Returns:
        Dataset
    """
    from . import adj
    if not isinstance(data, Dataset):
        raise ValueError("Requires a Dataset object", type(data))

    if not isinstance(name, str):
        raise ValueError("Requires a string name", type(name))

    if name not in data.data_vars:
        raise ValueError("dataset var not present")

    if adjname not in data.data_vars:
        raise ValueError("dataset var not present")

    if breakname not in data.data_vars:
        raise ValueError("requires a breaks dataset var")

    if suffix is not None:
        if suffix[0] != '_':
            suffix = '_' + suffix
            Warning('suffix needs an _. Added:', suffix)
    else:
        suffix = ''

    if percentilen is None:
        percentilen = np.arange(0, 101, 10)

    values = data[name].values.copy()
    avalues = data[adjname].values.copy()
    breaks = get_breakpoints(data[breakname], dim=dim, **ff.levelup(**kwargs))
    axis = data[name].dims.index(dim)
    params = data[name].attrs.copy()  # deprecated (xr-patch)

    ff.message(name, str(values.shape), 'A:', axis, 'Q:', np.size(percentilen), "Adj:", adjname, **kwargs)

    params.update({'sample_size': kwargs.get('sample_size', 130), 'borders': kwargs.get('borders', 90)})

    ff.message(ff.dict2str(params), **ff.levelup(**kwargs))
    #
    # Adjust reference to a reference period?
    #
    if adjust_reference:
        avalues = adj.percentile_reference_period(values, avalues, breaks, axis=axis, percentilen=percentilen, **kwargs)
    #
    # Adjust according to reference dataset
    #
    if False:
        values = adj.percentile_reference(values, avalues, breaks, axis=axis, percentilen=percentilen, **kwargs)
    else:
        #
        # use adjusted era as reference and calculate departures -> adj departures
        #
        values = values - avalues
        values = adj.percentile(values, breaks, axis=axis, percentilen=percentilen, **kwargs)
        values = values + avalues
        ff.message(name, 'using QA-adj departures', **kwargs)
        # values = adj.percentile_reference(values, avalues, breaks, axis=axis, percentilen=percentilen, **kwargs)
        
    data[name + suffix] = (data[name].dims, values)
    data[name + suffix].attrs.update(params)
    data[name + suffix].attrs['biascor'] = 'percentil_ref'
    data[name + suffix].attrs['reference'] = adjname

    if adjust_reference:
        #
        # fix for no breakpoints
        #
        if len(breaks) > 0:
            ref_period = data[dim].values[breaks[-1]].astype('M8[M]').astype('str') + ' -'
        else:
            ref_period = '-'

        data[adjname + suffix] = (data[adjname].dims, avalues)
        data[adjname + suffix].attrs.update(params)
        data[adjname + suffix].attrs['ref_period'] = kwargs.get('ref_period', ref_period)
        data[adjname + suffix].attrs['reference'] = name
    return data


def adjust_reference_period(data, name, refname, breakname, dim='time', suffix='_qa', percentilen=None, **kwargs):
    from . import adj
    if not isinstance(data, Dataset):
        raise ValueError("Requires a Dataset object", type(data))

    if not isinstance(name, str):
        raise ValueError("Requires a string name", type(name))

    if name not in data.data_vars:
        raise ValueError("dataset var not present")

    if refname not in data.data_vars:
        raise ValueError("dataset var not present")

    if breakname not in data.data_vars:
        raise ValueError("requires a breaks dataset var")

    if suffix is not None:
        if suffix[0] != '_':
            suffix = '_' + suffix
            Warning('suffix needs an _. Added:', suffix)
    else:
        suffix = ''

    if percentilen is None:
        percentilen = np.arange(0, 101, 10)

    values = data[refname].values.copy()  # RASO
    avalues = data[name].values.copy()  # Reanalysis (ERA)

    breaks = get_breakpoints(data[breakname], dim=dim, **ff.levelup(**kwargs))
    axis = data[name].dims.index(dim)
    params = data[name].attrs.copy()  # deprecated (xr-patch)

    ff.message(name, str(values.shape), 'A:', axis, 'Q:', np.size(percentilen), "Adj:", refname, **kwargs)

    params.update({'sample_size': kwargs.get('sample_size', 130), 'borders': kwargs.get('borders', 90)})

    ff.message(ff.dict2str(params), **ff.levelup(**kwargs))
    stdn = data[name].attrs.get('standard_name', name)
    #
    # Adjust name with refname in reference period
    #
    avalues = adj.percentile_reference_period(values, avalues, breaks, axis=axis, percentilen=percentilen, **kwargs)
    data[name + suffix] = (data[name].dims, avalues)
    data[name + suffix].attrs.update(params)
    data[name + suffix].attrs['standard_name'] = stdn + '_percentil_adj'
    return data


# def apply_bounds(data, name, other, lower, upper):
#     "Apply bounds and replace"
#     logic = data[name].values < lower
#     n = np.sum(logic)
#     data[name].values = np.where(logic, data[other].values, data[name].values)
#     logic = data[name].values > upper
#     n += np.sum(logic)
#     data[name].values = np.where(logic, data[other].values, data[name].values)
#     data[name].attrs['bounds'] = "[%d , %d]" % (lower, upper)
#     print("Outside bounds [", lower, "|", upper, "] :", n)
#

#
# def correct_loop(dataset, dep_var=None, use_dep=False, mean_cor=False, percentile_cor=False, percentile_adj=None,
#                  percentilen=None, clim_ano=True, **kwargs):
#     funcid = "[DC] Loop "
#     if not isinstance(dataset, DataArray):
#         raise ValueError(funcid + "Requires a DataArray class object")
#
#     if not mean_cor and not percentile_cor and percentile_adj is None:
#         raise RuntimeError(funcid + "Requires a correction: mean_cor, percentile_cor or percentile_adj")
#
#     if np.array([mean_cor, percentile_cor, percentile_adj is not None]).sum() > 1:
#         raise RuntimeError(funcid + "Only one Method at a time is allowed!")
#
#     xdata = dataset.copy()
#
#     # Make Large Arrays with all iterations ?
#     dataset = dataset.copy()
#     dims = dataset.get_dimension_values()
#     dims['iter'] = [0]
#     order = dataset.dims.list[:] + ['iter']
#     dataset.update_values_dims_remove(np.expand_dims(dataset.values, axis=-1), order, dims)
#     # dataset.dims['iter'].set_attrs({''})  # ?
#     sdata = dataset.copy()
#     sdata.values[:] = 0.
#     sdata.name += '_snht'
#     bdata = dataset.copy()
#     bdata.values[:] = 0.
#     bdata.name += '_breaks'
#     status = True
#     i = 1
#     while status:
#         status, stest, breaks, xdata = adjustments(xdata, dep_var=dep_var, use_dep=use_dep, mean_cor=mean_cor,
#                                                    percentile_cor=percentile_cor, percentile_adj=percentile_adj,
#                                                    percentilen=percentilen, clim_ano=clim_ano,
#                                                    **kwargs)
#         # combine
#         dataset.values = np.concatenate((dataset.values, np.expand_dims(xdata.values, axis=-1)), axis=-1)
#         # dataset.update_values_dims()
#         bdata.values = np.concatenate((bdata.values, np.expand_dims(breaks.values, axis=-1)), axis=-1)
#         sdata.values = np.concatenate((sdata.values, np.expand_dims(stest.values, axis=-1)), axis=-1)
#
#         # Does the adjustments still change anything ?
#         test = np.abs(np.nansum(dataset.values[:, :, i - 1] - xdata.values))  # sum of differences
#         if test < 0.1:
#             break
#         message(funcid + "%02d Breaks: \n" % i, **kwargs)
#         i += 1
#     # SAVE
#     dataset.update_values_dims(dataset.values, {'iter': range(i + 1)})
#     sdata.update_values_dims(sdata.values, {'iter': range(i + 1)})
#     bdata.update_values_dims(bdata.values, {'iter': range(i + 1)})
#     sdata.attrs['iterations'] = i
#
#     params = {'sample_size': kwargs.get('sample_size', 730),
#               'borders': kwargs.get('borders', 180),
#               'bounded': str(kwargs.get('bounded', '')),
#               'recent': kwargs.get('recent', False),
#               'ratio': kwargs.get('ratio', True)}
#
#     message(funcid + "Breaks: \n", **kwargs)
#     # print_breaks(bdata.subset(dims={'iter': i - 1}), verbose)
#
#     if mean_cor:
#         dataset.name += '_m_iter'
#         dataset.attrs['biascor'] = 'mean'
#         dataset.attrs['standard_name'] += '_mean_adj'
#         dataset.attrs.set_items(params)
#
#     elif percentile_cor:
#         dataset.name += '_q_iter'
#         dataset.attrs['biascor'] = 'percentile'
#         dataset.attrs['standard_name'] += '_percentile_adj'
#         dataset.attrs.set_items(params)
#
#     elif percentile_adj is not None:
#         dataset.name += '_qe_iter'
#         dataset.attrs['biascor'] = 'percentile_era_adjusted'
#         dataset.attrs['standard_name'] += '_percentile_era_adj'
#         dataset.attrs.set_items(params)
#     else:
#         pass
#
#     return status, sdata, bdata, dataset

#
# def adjust_table(data, name, analysis, dim='time', **kwargs):
#     """
#     test
# Out[23]:
# {'dpd':       dataset
#  mean     2
#  rmse     3
#  var      2, 'era':       M  Q
#  mean -3 -3
#  rmse  2  2
#  var   4  4}
#
#     pd.concat(test, axis=1)
# Out[22]:
#       dpd era
#      dataset   M  Q
# mean    2  -3 -3
# rmse    3   2  2
# var     2   4  4
#
#     Args:
#         data:
#         name:
#         analysis:
#         dim:
#         **kwargs:
#
#     Returns:
#
#     """
#     import pandas as pd
#     from ..fun import rmse
#
#     axis = data[name].dims.index(dim)
#     # for all reanalysis
#     out = {}
#     out[name] = {'dataset': {'RMSE': rmse(data[name], np.nanmean(data[name], axis=axis)),
#                              'MEAN': np.nanmean(data[name]),
#                              'VAR': np.nanvar(data[name])}}
#     for i, iana in enumerate(analysis):
#         tmp = data[[name, iana]].copy()
#         # snht
#         tmp = snht(tmp, dim=dim, var=name, dep=iana, **kwargs)
#         # threshold
#         tmp = apply_threshold(tmp, var=name + '_snht', dim=dim)
#         out[iana] = {}
#         out[iana]['n'] = len(get_breakpoints(tmp, dim=dim, var=name + '_snht_breaks'))
#         out[iana] = {'dataset': {'RMSE': rmse(tmp[name], tmp[iana]),
#                                  'MEAN': np.nanmean(tmp[name] - tmp[iana]),
#                                  'VAR': np.nanvar(tmp[name] - tmp[iana])}}
#         # adjust Mean
#         tmp = adjust_mean(tmp, name, name + '_snht_breaks', dim=dim, **kwargs)
#         out[iana]['mdiff'] = {'RMSE': rmse(tmp[name + '_m'], tmp[iana]),
#                               'MEAN': np.nanmean(tmp[name + '_m'] - tmp[iana]),
#                               'VAR': np.nanvar(tmp[name + '_m'] - tmp[iana])}
#         # adjust Percentiles
#         tmp = adjust_percentiles(tmp, name, name + '_snht_breaks', dim=dim, **kwargs)
#         out[iana]['qdiff'] = {'RMSE': rmse(tmp[name + '_q'], tmp[iana]),
#                               'MEAN': np.nanmean(tmp[name + '_q'] - tmp[iana]),
#                               'VAR': np.nanvar(tmp[name + '_q'] - tmp[iana])}
#         # adjust Reference
#         tmp = adjust_reference_period(tmp, iana, name, name + '_snht_breaks', dim=dim, **kwargs)
#         out[iana]['qrdiff'] = {'RMSE': rmse(tmp[iana + '_qa'], tmp[iana]),
#                                'MEAN': np.nanmean(tmp[iana + '_qa'] - tmp[iana]),
#                                'VAR': np.nanvar(tmp[iana + '_qa'] - tmp[iana])}
#         # adjust Percentiles using a Reference
#         tmp = adjust_percentiles_ref(tmp, name, iana, name + '_snht_breaks', dim=dim, **kwargs)
#         out[iana]['qadiff'] = {'RMSE': rmse(tmp[name + '_qa'], tmp[iana]),
#                                'MEAN': np.nanmean(tmp[name + '_qa'] - tmp[iana]),
#                                'VAR': np.nanvar(tmp[name + '_qa'] - tmp[iana])}
#
#     for ikey, idata in out.items():
#         out[ikey] = pd.DataFrame(idata)
#
#     return pd.concat(out, axis=1)

#
# def correct_2var(xdata, ydata):
#     # Make a 3D (time, var1, var2) per level Test Statistics
#     # Use that to adjustments both variables at the same time
#     # ? water vapor transform -> how to unsplit vp to t,rh ?
#     # t, rh -> td (esatfunc) -> vp
#     # large errors -> temperature problem ?
#     # smaller errors -> humidity problem ?
#     # t, rh percentage of contribution to vp
#     # vp (esat_inv) -> td
#     pass


def breakpoint_statistics(data, breakname, dim='time', variables=None, borders=None, inbetween=True, nmax=None,
                          **kwargs):
    """

    Args:
        data (Dataset): experiment data
        breakname (str): SNHT break variable
        dim (str): datetime dimension
        variables (list): variables to use
        borders (int): breakpoint borders
        inbetween (bool): calculate bordered area
        nmax (int): maximum values to use
        **kwargs:

    Returns:
        statistics (Dataset) : default nanmean breakpoint statistics before (B, later) and after (A, earlier) a breakpoint
    """
    from .. import fun as ff

    if not isinstance(data, Dataset):
        raise ValueError("Requires a Dataset class object", type(data))

    if dim not in data.coords:
        raise ValueError("Requires a datetime dimension", data.coords)

    if breakname not in data.data_vars:
        raise ValueError("Variable breakname not present", breakname, data.data_vars)

    ibreaks = get_breakpoints(data[breakname], value=kwargs.pop('breakpoint_threshold', 2), dim=dim)
    nb = len(ibreaks)
    if nb == 0:
        ff.message("Warning no Breakpoints found", **ff.leveldown(**kwargs))  # Always print
        return
    #
    # Variables to use ?
    #
    if variables is None:
        variables = list(data.data_vars)
    else:
        variables = [i for i in variables if i in data.data_vars]

    ibreakdates = list(data[dim].values[ibreaks].astype('M8[D]').astype('str'))
    variables.remove(breakname)
    variables = [i for i in variables if 'snht' not in i]

    if borders is None:
        borders = 0

    if nmax is None:
        nmax = 100000

    data = data[variables].copy()
    gattrs = data.attrs.copy()
    axis = data[variables[0]].dims.index(dim)
    wfunc = kwargs.pop('wfunc', ff.cal.nanfunc)
    region = {}
    j = 0
    ibreaks = ibreaks + [data[dim].size - 1]
    for i, k in enumerate(ibreakdates):
        #
        # Region left of breakpoint (A)
        #
        m = ibreaks[i + 1]
        i = ibreaks[i]
        region['A' + k] = data.isel(**{dim: slice(j, i)}).apply(ff.xarray.xarray_function_wrapper,
                                                                wfunc=wfunc,
                                                                dim=dim,
                                                                axis=axis,
                                                                borders=borders,
                                                                nmax=nmax,
                                                                **kwargs)
        #
        # Region right of breakpoint (B)
        #
        region['B' + k] = data.isel(**{dim: slice(i, m)}).apply(ff.xarray.xarray_function_wrapper,
                                                                wfunc=wfunc,
                                                                dim=dim,
                                                                axis=axis,
                                                                borders=borders,
                                                                nmax=nmax,
                                                                **kwargs)
        #
        # Region between borders at breakpoint (I)
        #
        if borders > 0 and inbetween:
            # Area around breakpoint [bordered]
            region['I' + k] = data.isel(**{dim: slice(i - borders, i + borders)}).apply(
                ff.xarray.xarray_function_wrapper,
                wfunc=wfunc,
                dim=dim,
                axis=axis,
                borders=0,
                nmax=nmax,
                **kwargs)
        ff.message("Break", j, i, m, k, **kwargs)
        j = i + borders

    data = concat(region.values(), dim=Index(region.keys(), name='region'))

    if hasattr(wfunc, '__name__'):
        if wfunc.__name__ == 'nanfunc':
            gattrs['statistic'] = "nanfunc(" + kwargs.get('func', 'nanmean') + ")"
        else:
            gattrs['statistic'] = wfunc.__name__
    else:
        gattrs['statistic'] = str(wfunc)

    if nmax is not 100000:
        gattrs['max_sample'] = nmax

    if borders > 0:
        gattrs['borders'] = borders
        if inbetween:
            gattrs['inbetween'] = True

    data.attrs.update(gattrs)
    return data


def breakpoint_info(data, snhtname, breakname, dim='time', thres=50, **kwargs):
    if not isinstance(data, Dataset):
        raise ValueError('Requires an xarray Dataset', type(data))

    if snhtname not in data.data_vars:
        raise ValueError('Requires a variable name: snhtname', list(data.data_vars))

    if breakname not in data.data_vars:
        raise ValueError('Requires a variable name: breakname', list(data.data_vars))

    if dim not in data.dims:
        raise ValueError('requires a datetime dimension', dim)

    # get breakpoints
    breaks = get_breakpoints(data[breakname], dim=dim, **kwargs)
    # how many levels
    for ibreak in breaks:
        # 0: no significant, 1: significant, 2: significant at other level, 3; significant at level
        print(data[dim][ibreak], data[breakname].isel({dim: ibreak}).values,
              data[snhtname].isel({dim: ibreak}).values)
        # message()

    # look at snht and check how close
    #


# shape = list(dataset[name].values.shape)
#
# dep = {getattr(ifunc, '__name__'): [] for ifunc in functions}
#
# dep['counts'] = []
# dates = dataset.coords[dim].values
# jbreaks = sorted(ibreaks, reverse=True)
# jbreaks.append(0)
# idims = list(dataset[name].dims)
# jdims = idims.copy()
# jdims.pop(axis)
# func_kwargs.update({'axis': axis})
# #
# # iterate from now to past breakpoints
# #
# for i, ib in enumerate(break_iterator(ibreaks, axis, shape, borders=borders, max_sample=max_sample)):
#     period = vrange(dates[ib[axis]])
#     idate = dates[jbreaks[i]]
#     tmp = np.sum(np.isfinite(dataset[name][ib]), axis=axis)  # is an DataArray
#     tmp.coords[dim] = idate
#     tmp.coords['start'] = period[0]
#     tmp.coords['stop'] = period[1]
#     dep['counts'].append(tmp.copy())  # counts
#
#     for j, ifunc in enumerate(functions):
#         iname = getattr(ifunc, '__name__')
#         # Requires clear mapping of input and output dimensions
#         tmp = apply_ufunc(ifunc, dataset[name][ib],
#                           input_core_dims=[idims],
#                           output_core_dims=[jdims],
#                           kwargs=func_kwargs)
#         # tmp = ifunc(dataset[name][ib], axis=axis, **func_kwargs)
#         # only for functions with ufunc capability
#         tmp.coords[dim] = idate
#         tmp.coords['start'] = period[0]
#         tmp.coords['stop'] = period[1]
#         dep[iname].append(tmp.copy())
#
# for ifunc, ilist in dep.items():
#     dep[ifunc] = concat(ilist, dim=dim)
#
# dep = Dataset(dep)
# return dep


def reference_period(data, dim='time', dep_var=None, period=None, **kwargs):
    from ..met.time import anomaly

    if not isinstance(data, DataArray):
        raise ValueError("Requires a DataArray class object", type(data))

    if dim not in data.dims:
        raise ValueError("Requires a datetime dimension", data.dims)

    data = data.copy()
    # attrs ?

    if dep_var is not None:
        if not isinstance(dep_var, DataArray):
            raise ValueError("Requires a DataArray class object", type(dep_var))

        dep = data - dep_var  # Departures (do the units match?)
    else:
        dep, _ = anomaly(data, dim=dim, period=period)

    #
    # find best matching period (lowest differences)
    # run SNHT
    # Split into pieces
    # run stats on each piece (RMSE)
    # choose
    # return piece + index
    return None


def combine_metadata(data, dim='time', window=30,
                   lon=None, lat=None, distance_weight=1, distance_threshold=10,
                   read_igra=True, meta_ident=None, meta_weight=1,
                   stype=None, sonde_weight=1, **kwargs):
    from .meta import location_change, metadata, sondetype
    if not isinstance(data, Dataset):
        raise ValueError()

    # can be redundant


    if lon is not None and lat is not None:
        # distance in [km] of location changes
        dinfo = location_change(lon, lat, dim=dim, **kwargs)
        dinfo.values = np.where(dinfo.values > distance_threshold, distance_weight, 0)
        # triangle shape
        dinfo = dinfo.rolling(**{dim: window}, min_periods=1, center=True).sum().rolling(**{dim: window},
                                                                                               min_periods=1,
                                                                                               center=True).mean()
    if stype is not None:
        sinfo = sondetype(stype, dim, window=window, **kwargs)

    if read_igra:
        if meta_ident is None:
            raise ValueError('')
        minfo = metadata(meta_ident, data[dim].values, dim=dim, window=window, **kwargs)

    return data

def apply_biasadjustments(adjdata, data, isodb=False, **kwargs):
    # calculate the bias adjustmens and for sounding times
    # interpolate between sounding times
    # interpolate to table format?
    # todo finish that function
    # cal. adjustments adjdata - data = adj
    # interpolate?, unify across missing
    # quantile adjustments ? how to interpolate these across unknown
    pass
