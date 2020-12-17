

# try:
#     from xData import DataArray
#     from xData.fun import date_selection, date_groupby
#
# except ImportError:
#     pass

__all__ = ['percentile_count', 'percentiles']


def percentiles(data, percentiles=None, freq='M', **kwargs):
    """ Calculate percentiles per freq

    Args:
        data:
        percentiles:
        freq:
        **kwargs:

    Returns:

    """
    import numpy as np
    import pandas as pd
    from xarray import DataArray

    if not isinstance(data, DataArray):
        raise ValueError("Requires a DataArray class object")

    data = data.copy()
    date_dim = data.get_date_dimension()
    dates = data.dims[date_dim].values.copy()

    if percentiles is None:
        percentiles = range(10, 100, 10)
    else:
        percentiles = sorted(list(set(list(percentiles) + [0, 100]).difference(set([0, 100]))))  # remove 0 and 100

    print("Percentiles: ", percentiles)

    axis = data.order.index(date_dim)
    grouped = date_groupby(dates, axis, len(data.order), freq='M')  # monthly
    dates = pd.DatetimeIndex(dates).to_period('M').unique().to_timestamp().values
    values = []
    icounts = np.array([np.isfinite(data.values[ig]).sum(axis=axis) for ig in grouped])  # ratio per month
    for iq in percentiles:
        with np.errstate(invalid='ignore', over='ignore'):
            # time (freq) x other dimensions
            tmp = np.array([np.nanpercentile(data.values[ig], iq, axis=axis) for ig in grouped])  # > axis 0
            tmp = np.where(icounts > 30., tmp, np.nan)  # filter for too little dataset
            values.append(tmp)

    values = np.array(values)  # percentiles x time x other
    dims = data.get_dimension_values()
    order = list(data.order)
    order.insert(0, 'perc')
    dims['perc'] = percentiles  # add percentile dim
    dims[date_dim] = dates
    data.update_values_dims_remove(values, order, dims)
    data.name += '_perc'
    data.attrs['standard_name'] += '_%s_percentiles' % freq
    return data


def percentile_count(data, percentiles=None, freq='M', period=None):
    """ Calculate monthly percentile counts from a reference period

    Args:
        data:
        percentiles:
        freq:
        period:

    Returns:

    """
    if not isinstance(data, DataArray):
        raise ValueError("Requires a DataArray class object")

    data = data.copy()
    date_dim = data.get_date_dimension()
    order = list(data.order)  # copy?
    lev_dim = data.get_dimension_by_axis('Z')  # Vertical
    dates = data.dims[date_dim].values.copy()
    idate = order.index(date_dim)
    if percentiles is None:
        percentiles = range(10, 100, 10)
    else:
        percentiles = sorted(list(set(list(percentiles) + [0, 100]).difference(set([0, 100]))))  # remove 0 and 100

    nperc = len(percentiles) + 1  # +1 for residual class

    if period is None:
        period = slice('2010', None)

    itx = [slice(None, None)] * data.values.ndim
    itx[idate] = date_selection(dates, period)  # only the reference period

    ref_q = np.apply_along_axis(_unique_percentiles, idate, data.values[itx], percentiles)

    # old
    # nlev = dataset.values.shape[1]
    # ref_q = np.array([_unique_percentiles(dataset.values[itx, ilev], percentiles) for ilev in range(nlev)])  #

    grouped = date_groupby(dates, idate, len(data.order), freq=freq)
    dates = pd.DatetimeIndex(dates).to_period(freq).unique().to_timestamp().values

    shapes = list(data.values.shape)  # original shape of dataset
    oshapes = shapes[:]  # copy
    shapes.insert(idate+1, nperc)  # insert perc dim
    shapes[idate] = dates.size  # update date dim

    counts = np.zeros(tuple(shapes))
    oshapes.pop(idate)
    # (:,x,y,:) <- (x,y) <- (:,x,:) (:,x)
    for i in np.ndindex(tuple(oshapes)):
        idx = list(i)
        idx.insert(idate, slice(None))
        idy = idx[:]
        idy.insert(idate, slice(None))
        counts[idy] = np.array([_count(data.values[ig][idx], ref_q[idx]) for ig in grouped])
        # counts[idx] = _count(dataset.values[idx], ref_q[idx])

    # counts = np.array(counts).reshape(shapes)
    # counts = np.array([[_count(dataset.values[ig][:, ilev], ref_q[ilev, :]) for ilev in range(nlev)] for ig in grouped])
    tmp = data.copy()
    dims = tmp.get_dimension_values()
    dims.pop(date_dim)
    order[idate] = 'perc'
    order.insert(idate, date_dim)
    dims[date_dim] = dates  # months.to_datetime().shift(15, 'D').values
    dims['perc'] = percentiles + [100]
    tmp.update_values_dims_remove(counts, order, dims)
    tmp.dims[date_dim].set_attrs({'axis': 'T', 'freq': freq})
    tmp.axes[0] = 'T'
    tmp.dims['perc'].set_attrs({'units': '1', 'standard_name': 'percentile'})
    tmp.attrs['standard_name'] += '_percentile_counts'
    tmp.attrs['units'] = '1'
    tmp.attrs['name'] = 'pcounts'

    order = order[:]
    order.remove(date_dim)
    dims = {}
    dattrs = tmp.get_dimension_attrs()
    for idim in order:
        if idim == 'perc':
            dims[idim] = tmp.dims[idim].values[:-1]
        else:
            dims[idim] = tmp.dims[idim].values

    qs = DataArray('percentiles', ref_q, tuple(order), dims,
                   attrs={'standard_name': data.attrs['standard_name'] + '_percentile', 'units': data.attrs['units'],
                          'period': str(period)},
                   dim_attrs=dattrs)

    return tmp, qs


def _count(data, percentiles):
    n = len(percentiles) + 1
    counts = np.zeros(n)
    itx = np.squeeze(np.isfinite(data))
    if np.sum(itx) < 0.3 * data.size:
        return counts

    with np.errstate(invalid='ignore'):
        binned = np.digitize(data[itx], np.squeeze(percentiles))  # qs-1 <= x < qs
        counts[:] = np.bincount(binned, minlength=n)  # count values inside percentile bins
        counts /= np.sum(counts)
    return counts * 100.  # Percent


def _unique_percentiles(data, percentiles):
    qs = np.nanpercentile(data, percentiles)  # Monthly percentiles
    iq = np.round(qs, decimals=6)   # round to make
    if len(np.unique(iq)) != len(iq):
        for j in range(1, len(percentiles)):
            if iq[j - 1] == iq[j]:
                iq[j] = np.ceil(qs[j] * 1e6) / 1e6
    return iq
