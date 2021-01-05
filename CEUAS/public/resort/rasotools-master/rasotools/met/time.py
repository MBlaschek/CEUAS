# -*- coding: utf-8 -*-


__all__ = ['anomaly', 'trend', 'trend_per_month', 'trend_mon_percentile',
           'day_night_departures', 'correlate', 'covariance', 'resample_nanmean']
#
# global Variables
#
month_to_season = {1: 'DJF', 2: 'DJF', 3: 'MAM', 4: 'MAM', 5: 'MAM', 6: 'JJA', 7: 'JJA', 8: 'JJA', 9: 'SON',
                   10: 'SON', 11: 'SON', 12: 'DJF'}


#
# # Counting and Selecting
#


def select_period(data, dim='time', period=None):
    from ..fun.cal import nanrange
    from xarray import DataArray, Dataset

    if not isinstance(data, (DataArray, Dataset)):
        raise ValueError('requires an xarray DataArray, Dataset', type(data))

    if dim not in data.dims:
        raise ValueError("datetime dimension not found")

    if period is None:
        data.attrs['period'] = '%d-%d' % nanrange(data[dim].dt.year.values)
        return data
    else:
        iperiod = '%d-%d' % nanrange(data[dim].to_series()[period].index.year.values)
        data = data.sel(**{dim: period})
        data.attrs['period'] = iperiod
        return data


def estimate_sample_size(data, ratio=0.6, freq='12h'):
    """ Estimate the sample size from a timeseries, given freq and ratio

    Args:
        data (DataArray): Inputdata
        ratio (float): percentage of dataset as ratio
        freq (str): Pandas freq str

    Returns:
        int : Sample size according to freq and ratio
    """
    import numpy as np
    import pandas as pd
    from xarray import DataArray
    from ..fun import cal as fc

    if not isinstance(data, DataArray):
        raise ValueError("Requires a DataArray class object")

    date_dim = data.get_date_dimension()
    dates = pd.DatetimeIndex(data.dims[date_dim])
    axis = data.order.index(date_dim)
    dates = pd.date_range(dates.min().replace(hour=np.min(dates.hour.values)),
                          dates.max().replace(hour=np.max(dates.hour.values)), freq=freq)
    years = fc.nanrange(dates.year)
    print("Estimate Sample size (%d%%, F:%s): %d [%d : %d] %d %d" % (
        int(ratio * 100), freq, int(dates.size * ratio), years[0], years[1], np.diff(years) + 1, dates.size))
    print(100 * data.apply(fc.nancount, axis=axis).to_pandas() / float(dates.size))
    return int(dates.size * ratio)


#
# Climatology and Anomalies
#


def climatology(data, dim='time', period=None, keep_attrs=True):
    """

    Args:
        data (DataArray): Input Data
        dim (str): datetime dimension
        period (slice): datetime selection
        keep_attrs (bool) : xarray keep attributes
    Returns:
        DataArray : Climate Monthly Means
    """
    from xarray import DataArray, Dataset
    if not isinstance(data, (DataArray, Dataset)):
        raise ValueError("Requires a xarray DataArray, Dataset", type(data))

    if dim not in data.dims:
        raise ValueError("datetime dimension not found")

    data = select_period(data, dim=dim, period=period)  # adds metadata and selects if not None
    return data.groupby(dim + '.month').mean(dim, keep_attrs=keep_attrs)


def anomaly(data, dim='time', period=None, keep_attrs=True):
    """ Calculates the anomaly from the climatology per month of a time series

    Args:
        data (DataArray) : Inputdata
        dim (str) : datetime dimension
        period (slice, str) : Indices of Dates for calculation
        keep_attrs (bool) : xarray keep attributes
    Returns:
        DataArray : Anomalies
    """
    from xarray import DataArray, Dataset, set_options
    from ..fun.xarray import set_attrs
    if not isinstance(data, (DataArray, Dataset)):
        raise ValueError("Requires a xarray DataArray, Dataset", type(data))

    if dim not in data.dims:
        raise ValueError("datetime dimension not found")

    data = data.copy()
    # Calculate Climatology
    clim = climatology(data, dim=dim, period=period, keep_attrs=keep_attrs)
    # Calculate Anomaly
    with set_options(keep_attrs=keep_attrs):
        data = data.groupby(dim + '.month') - clim

    data = data.drop('month')
    if isinstance(data, Dataset):
        for ivar in data.data_vars:
            set_attrs(data[ivar], 'standard_name', add='_ano', default='anomaly')
            data[ivar].attrs['period'] = clim.attrs['period']
    else:
        set_attrs(data, 'standard_name', add='_ano', default='anomaly')
        data.attrs['period'] = clim.attrs['period']
    return data


#
# Trend Estimation
#


def trend(data, dim='time', use_anomalies=True, period=None, min_periods=3, method='theil_sen',
          alpha=0.95, keep_attrs=True, only_slopes=False, **kwargs):
    """ Calculate Trend estimates

    Args:
        data (DataArray): input dataset array
        dim (str): datetime dimension
        use_anomalies (bool): calc. trends from anomalies (climatology removed)
        period (slice): time period for climatology
        min_periods (int): minimum number of values for trend estimate
        method (str): polyfit, theil_sen, linregress, lsq
        alpha (float): get confidence levels for that p value
        keep_attrs (bool): keep DataArray Attributes?
        only_slopes (bool): return only slopes (e.g. for Datasets)

    Returns:
        DataArray : trends
    """
    import numpy as np
    from xarray import DataArray, Dataset
    from .. import fun as ff

    if not isinstance(data, (DataArray, Dataset)):
        raise ValueError("Requires a DataArray class object")

    if dim not in data.dims:
        raise ValueError("datetime dimension not found")

    if method not in ['polyfit', 'theil_sen', 'linregress', 'lsq']:
        raise ValueError("Requires either polyfit, theil_sen, linregress or lsq")

    if isinstance(data, Dataset):
        if use_anomalies:
            data = anomaly(data, dim=dim, period=period, keep_attrs=keep_attrs)

        out = {}
        for ivar in data.data_vars:
            out[ivar] = trend(data[ivar], dim=dim, use_anomalies=False, method=method, keep_attrs=keep_attrs,
                              only_slopes=True)

        out = Dataset(out)
        out.attrs.update(data.attrs.copy())
        return out

    data = data.copy()
    per = np.timedelta64(1, 'D') / np.timedelta64(1, 'ns')  # factor for trends
    axis = data.dims.index(dim)
    coords = {idim: data[idim].copy() for idim in data.dims if idim != dim}
    dimens = list(data.dims[:])
    dimens.remove(dim)
    attrs = data.attrs.copy()

    if use_anomalies:
        data = anomaly(data, dim=dim, period=period, keep_attrs=keep_attrs)
        attrs['period'] = data.attrs['period']  # copy

    # Convert to standard time axis
    idates = data[dim].values.astype('long')  # Nano Seconds
    idates -= idates[0]  # relative Times
    # Trends
    # k = [unit]/time
    params = ff.cal.linear_trend(data.values, idates, method=method, alpha=alpha, nmin=min_periods, axis=axis)
    # slope and intercept
    idx = [slice(None)] * params.ndim
    idx[axis] = 0  # slope
    slope = DataArray(params[tuple(idx)] * per, coords=coords, dims=dimens, name='slope', attrs=attrs)
    ff.xarray.set_attrs(slope, 'units', add='/day', default='1/day')
    ff.xarray.set_attrs(slope, 'standard_name', add='_trend', default='trend')
    slope.attrs['cell_method'] = 'daily trend of anomalies' if use_anomalies else 'daily trend'
    if only_slopes:
        return slope

    idx[axis] = 1  # slope
    interc = DataArray(params[tuple(idx)], coords=coords, dims=dimens, name='intercept', attrs=attrs)
    ff.xarray.set_attrs(interc, 'standard_name', add='_intercept', default='intercept')

    if params.shape[axis] > 2:
        if method == 'theil_sen':
            idx[axis] = 2  # slope lower
            aslope = DataArray(params[tuple(idx)] * per, coords=coords, dims=dimens, name='slope_min', attrs=attrs)
            ff.xarray.set_attrs(aslope, 'units', add='/day', default='1/day')
            ff.xarray.set_attrs(aslope, 'standard_name', add='_trend_min', default='trend_min')
            aslope.attrs['alpha'] = alpha
            aslope.attrs['cell_method'] = 'daily trend of anomalies' if use_anomalies else 'daily trend'

            idx[axis] = 3  # slope upper
            bslope = DataArray(params[tuple(idx)] * per, coords=coords, dims=dimens, name='slope_max', attrs=attrs)
            ff.xarray.set_attrs(bslope, 'units', add='/day', default='1/day')
            ff.xarray.set_attrs(bslope, 'standard_name', add='_trend_max', default='trend_max')
            bslope.attrs['alpha'] = alpha
            bslope.attrs['cell_method'] = 'daily trend of anomalies' if use_anomalies else 'daily trend'

            return Dataset({'slope': slope, 'intercept': interc, 'lower': aslope, 'upper': bslope})

        # r_value, p_value, std_err
        idx[axis] = 2  # R-value
        rslope = DataArray(params[tuple(idx)] ** 2, coords=coords, dims=dimens, name='r_squared', attrs=attrs)
        rslope.attrs['units'] = '1'
        ff.xarray.set_attrs(rslope, 'standard_name', add='_r_squared', default='r_squared')

        idx[axis] = 3  # p-value
        bslope = DataArray(params[tuple(idx)], coords=coords, dims=dimens, name='p_value', attrs=attrs)
        bslope.attrs['units'] = '1'
        ff.xarray.set_attrs(bslope, 'standard_name', add='_p_value', default='p_value')
        bslope.attrs['cell_method'] = 'p-value for null hypothesis(slope==0)'

        idx[axis] = 4  # std err
        sslope = DataArray(params[tuple(idx)], coords=coords, dims=dimens, name='std_err', attrs=attrs)
        ff.xarray.set_attrs(sslope, 'units', add='/day', default='1/day')
        ff.xarray.set_attrs(sslope, 'standard_name', add='_std_err', default='std_err')
        sslope.attrs['cell_method'] = 'standard error of slope'

        return Dataset({'slope': slope, 'intercept': interc, 'r_squared': rslope, 'p_value': bslope, 'std_err': sslope})

    return Dataset({'slope': slope, 'intercept': interc})


def trend_mon_percentile(data, dim='time', percentile=None, period=None,
                         min_periods=3, min_per_month=15, method='lsq', **kwargs):
    """ Monthly percentile trends

    Args:
        data (DataArray): input dataset
        dim (str): datetime dimension
        percentile (list): percentiles, int 1-99
        period (slice): datetime period for climatology
        min_periods (int): minimum values for trend
        min_per_month (int): minimum monthly count
        method (str): trend method
        **kwargs:

    Returns:
        Dataset : slope_perc_XX  for each percentile
    """
    import numpy as np
    from xarray import DataArray
    from .. import fun as ff
    if not isinstance(data, DataArray):
        raise ValueError("Requires a DataArray class object")

    if percentile is None:
        percentile = [25, 50, 75]  # Quartils
    else:
        if any([iq < 1 for iq in percentile]):
            raise ValueError('Percentiles need to be integers [1, 99]')

    data = data.copy()
    axis = data.dims.index(dim)
    #
    # Call wrapper for nanpercentile -> add as new dimension
    #
    tmp = data.resample(**{dim: 'M'}).apply(ff.xarray.xarray_function_wrapper,
                                            wfunc=ff.cal.sample_wrapper,
                                            add_dim='prc',
                                            dim=dim,
                                            axis=axis,
                                            ffunc=np.nanpercentile,
                                            nmin=min_per_month,
                                            q=percentile)
    #
    # Call trend with
    #
    trends = trend(tmp, dim=dim, period=period, use_anomalies=False,
                   min_periods=min_periods, method=method, only_slopes=True, **kwargs)
    #
    # Add metadata
    #
    ff.xarray.set_attrs(trends, 'standard_name', add='_perc', default='percentiles')
    trends.attrs['cell_method'] = 'daily trend of monthly percentiles'
    trends.attrs['min_per_month'] = min_per_month
    #
    # Assign Coordinate information
    #
    trends = trends.assign_coords({'prc': percentile})
    trends['prc'].attrs.update({'units': '%', 'standard_name': 'percentile'})
    return trends


def trend_per_month(data, dim='time', **kwargs):
    """ Trends per month

    Args:
        data (DataArray): input dataset
        dim (str): datetime dimension
        **kwargs:

    Returns:
        DataArray :
    """
    from xarray import DataArray, Dataset
    if not isinstance(data, (DataArray, Dataset)):
        raise ValueError("Requires a xarray DataArray, Dataset", type(data))

    if dim not in data.dims:
        raise ValueError("datetime dimension not found")

    trends = data.groupby(dim + '.month').apply(trend, dim=dim, use_anomalies=False, **kwargs)
    return trends


#
# Correlations
#


def correlate(x, y, dim='time', period=None, method='spearman', **kwargs):
    """ Correlation between Arrays

    Args:
        x (DataArray): input dataset
        y (DataArray): input dataset
        dim (str): datetime dimension
        period (slice): consider only that datetime period
        method (str): either spearman or pearson
        **kwargs:

    Returns:
        DataArray : correlation coefficients
    """
    from xarray import DataArray, align
    from .. import fun as ff

    if not isinstance(x, DataArray):
        raise ValueError("Requires a DataArray class object")

    if not isinstance(y, DataArray):
        raise ValueError("Requires a DataArray class object")

    if method not in ['spearman', 'pearson']:
        raise ValueError('Only spearman or pearson allowed')

    if dim not in x.dims or dim not in y.dims:
        raise ValueError('Dimension must be present in both Arrays')

    x = select_period(x, dim=dim, period=period)
    # Align
    x, y = align(x, y, join='left')
    axis = x.dims.index(dim)

    # def sp_corr(xx, yy, d, a):
    #     jdims = list(xx.dims)
    #     jdims.remove(d)
    #     return apply_ufunc(ff.cal.spearman_correlation, xx, yy,
    #                        input_core_dims=[xx.dims, yy.dims],
    #                        output_core_dims=[jdims],
    #                        output_dtypes=[float],
    #                        kwargs={'axis': a},
    #                        keep_attrs=True)
    #
    # def ps_corr(xx, yy, d, a):
    #     jdims = list(xx.dims)
    #     jdims.remove(d)
    #     return apply_ufunc(ff.cal.pearson_correlation, xx, yy,
    #                        input_core_dims=[xx.dims, yy.dims],
    #                        output_core_dims=[jdims],
    #                        output_dtypes=[float],
    #                        kwargs={'axis': a},
    #                        keep_attrs=True)

    if method == 'spearman':
        corr = ff.xarray.xarray_function_wrapper(x, wfunc=ff.cal.spearman_correlation, dim=dim, y=y, axis=axis)
        # corr = sp_corr(x, y, dim, axis)
    else:
        corr = ff.xarray.xarray_function_wrapper(x, wfunc=ff.cal.pearson_correlation, dim=dim, y=y, axis=axis)
        # corr = ps_corr(x, y, dim, axis)

    ff.xarray.set_attrs(corr, 'standard_name', add='_corr', default='correlation')
    corr.attrs['units'] = '1'
    corr.attrs['cell_method'] = '%s correlation with %s' % (method, y.name)
    return corr


def covariance(x, y, dim='time', period=None):
    """ Covariance

    Args:
        x:
        y:
        dim:
        period:

    Returns:

    """
    from xarray import DataArray, align
    from .. import fun as ff

    if not isinstance(x, DataArray):
        raise ValueError("Requires a DataArray class object")

    if not isinstance(y, DataArray):
        raise ValueError("Requires a DataArray class object")

    if dim not in x.dims or dim not in y.dims:
        raise ValueError('Dimension must be present in both Arrays')

    x = select_period(x, dim=dim, period=period)
    # Align
    x, y = align(x, y, join='left')
    axis = x.dims.index(dim)

    # def nancov(xx, yy, d, a):
    #     jdims = list(xx.dims)
    #     jdims.remove(d)
    #     return apply_ufunc(ff.cal.covariance, xx, yy,
    #                        input_core_dims=[xx.dims, yy.dims],
    #                        output_core_dims=[jdims],
    #                        output_dtypes=[float],
    #                        kwargs={'axis': a},
    #                        keep_attrs=True)
    #
    # corr = nancov(x, y, dim, axis)
    corr = ff.xarray.xarray_function_wrapper(x, wfunc=ff.cal.covariance, dim=dim, y=y, axis=axis)
    ff.xarray.set_attrs(corr, 'standard_name', add='_cov', default='covariance')
    ff.xarray.set_attrs(corr, 'units', add='2', default='2')
    ff.xarray.set_attrs(corr, 'cell_method', set='covariance with %s' % y.name)
    return corr


def day_night_departures(data, dim='time', day=12, night=0, **kwargs):
    """ Day-Night departures form dataset

    Args:
        data (DataArray): input dataset
        dim (str): datetime dimension
        day (int): hour of day: 12Z
        night (int): hour of night: 0Z
        **kwargs:

    Returns:

    """
    from ..fun.xarray import set_attrs
    from .std import to_hours
    from xarray import DataArray, Dataset, set_options

    if not isinstance(data, (DataArray, Dataset)):
        raise ValueError('Requires a DataArray, Dataset', type(data))

    if dim not in data.dims:
        raise ValueError('Requires a datetime dimension', dim)

    data = to_hours(data, dim=dim, times=[day, night], **kwargs)
    attrs = data.attrs.copy()
    with set_options(keep_attrs=True):
        data = data.sel(hour=day) - data.sel(hour=night)
    data.attrs.update(attrs)
    if isinstance(data, DataArray):
        # data.name = data.name + '_dep'
        set_attrs(data, 'standard_name', add='_day_night_dep', default='day_night_departure')
        data.attrs['cell_method'] = 'noon - night'
        data.attrs['info'] = 'Day(%dZ)-Night(%dZ)' % (day, night)
    else:
        # data = data.rename({iname: iname + '_dep' for iname in data.data_vars})
        for iname in data.data_vars:
            set_attrs(data[iname], 'standard_name', add='_day_night_dep', default='day_night_dep')
            data[iname].attrs.update({'cell_method': 'noon - night', 'info': 'Day(%dZ)-Night(%dZ)' % (day, night)})
    return data


def statistics(x, f='rmse', y=None, dim='time', period=None, **kwargs):
    from xarray import DataArray
    from .. import fun as ff

    if not isinstance(x, DataArray):
        raise ValueError("Requires a DataArray, not", type(x))

    if isinstance(f, str):
        try:
            f = getattr(ff.cal, f)
        except Exception as e:
            print('Function', f, 'not found in', ff.cal)
            raise e

    if period is not None:
        x = x.sel(**{dim: period})
        if y is not None:
            y = y.sel(**{dim: period})

    return ff.xarray.xarray_function_wrapper(x, wfunc=f, dim=dim, y=y, axis=x.dims.index(dim))


def resample_nanmean(data, dim='time', resample='M', nmin=15, **kwargs):
    """ Resample Dataset and apply a minimum for resampling mean

    Args:
        data (DataArray, Dataset): Input fields
        dim (str): datetime dimension
        resample (str): upsampling frequency
        nmin (int): minimum of samples per frequency
        **kwargs:

    Returns:
        DataArray, Dataset : means on freq
    """
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        dd = data.resample(**{dim: resample}).count(dim)
        return data.resample(**{dim: resample}).mean(dim).where(dd > nmin)

#
# def fill_missing_hours(data, dim='time', hour='hour', times=(0, 12), **kwargs):
#     import numpy as np
#     from pandas import Index
#     from xarray import DataArray, concat
#     from ..fun import message
#
#     if not isinstance(data, DataArray):
#         raise ValueError()
#
#     if dim not in data.coords.keys():
#         raise ValueError()
#     if hour not in data.coords.keys():
#         raise ValueError()
#
#     if not np.isin(np.array(times), data[hour].values).all():
#         raise ValueError()
#
#     #
#     # hours not 0, 12 -> fill in to
#     #
#     data = data.copy()
#     data = dict(data.groupby(hour))
#     for ikey, idata in data.items():
#         if ikey in times:
#             continue
#
#         if ikey >= 18 or ikey < 6:
#             jkey = 0
#         else:
#             jkey = 12
#         #
#         # from earlier times ?
#         #
#         logic = np.isfinite(data[jkey].values)
#         data[jkey].values = np.where(logic, data[jkey].values, idata.values)
#         message(ikey, " (%s)>(%s) %d" % (",".join(["%d" % i for i in data.keys()]),
#                                          ",".join(["%d" % i for i in times]),
#                                          sum(~logic & np.isfinite(idata.values)).sum()), **kwargs)
#     return concat([data[0], data[12]], dim=Index(times, name=hour, )).sortby(dim)
