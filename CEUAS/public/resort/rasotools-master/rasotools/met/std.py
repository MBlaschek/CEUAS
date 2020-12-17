# -*- coding: utf-8 -*-

__all__ = ['to_hours', 'reindex', 'sel_hours', 'align_datetime', 'datetime_range']


def sel_hours(data, dim='time', times=(0, 12), minutes=None, **kwargs):
    """ Select only given hours

    Args:
        data (xarray.DataArray): Input DataArray
        dim (str): datetime dimension
        times (tuple, list): hours to consider
        minutes (int): select only these minutes

    Returns:
        xarray.DataArray : DataArray only at selected hours
    """
    from ..fun import message
    from xarray import DataArray, Dataset

    if not isinstance(data, (DataArray, Dataset)):
        raise ValueError('Requires a xarray DataArray, Dataset', type(data))

    if dim not in data.dims:
        raise ValueError('Requires a datetime dimension', dim)

    message('Selecting times:', times, dim, data[dim].size, **kwargs)
    #
    # indexing with an array
    #
    if minutes is not None:
        return data.sel(**{dim: (data[dim].dt.hour.isin(times) & (data[dim].dt.minute == minutes))})
    return data.sel(**{dim: data[dim].dt.hour.isin(times)})  # copy


def reindex(data, dim='time', times=(0, 12), freq='12h', span=6, **kwargs):
    """ Reindex to standard times

    Args:
        data:
        dim:
        times:
        freq:
        span:

    Returns:

    """
    from ..fun import message

    if dim not in data.dims:
        raise ValueError("Dimension not present", dim, data.dims)

    dates = data[dim].values
    alldates = datetime_range(dates, times=times, freq=freq, span=span, **kwargs)
    message("New Index:", alldates[0], "to", alldates[-1], "with", freq, alldates.size, **kwargs)
    return data.reindex({dim: alldates})


def to_hours(data, dim='time', times=(0, 12), as_dataset=False, hour='hour', minutes=None, **kwargs):
    """ Split Array into separate Arrays by time

    Args:
        data (xarray.DataArray): Input dataset
        dim (str): datetime dimension
        times (tuple, list): std hours
        as_dataset (bool): return hour dim as variables
        hour (str): name of hour dimension
        minutes (int): use only these minutes

    Returns:
        xarray.DataArray : datetime dimension split to days and hours
    """
    from .. import fun as ff
    from pandas import Index
    from xarray import DataArray, Dataset, concat

    if not isinstance(data, (DataArray, Dataset)):
        raise ValueError('Requires a xarray DataArray or Dataset', type(data))

    if dim not in data.dims:
        raise ValueError('Requires a datetime dimension', dim)
    #
    # check if duplicated times will occur because of minutes
    if data[dim].to_index().to_period('h').duplicated().any():
        if minutes is None:
            minutes = 0
            ff.message("Warning minutes=0, datetime duplicates", **ff.leveldown(**kwargs))

    data = sel_hours(data, dim=dim, times=times, minutes=minutes, **ff.levelup(**kwargs))  # selection (copy)
    #
    data = dict(data.groupby(dim + '.hour'))
    for ikey in data.keys():
        if isinstance(data[ikey], Dataset):
            ff.message(ikey, ff.dict2str(data[ikey].dims), **kwargs)
        else:
            ff.message(ikey, data[ikey].shape, **kwargs)
        data[ikey] = data[ikey].assign_coords({dim: data[ikey][dim].to_index().to_period('D').to_timestamp().values})

    data = concat(data.values(), dim=Index(data.keys(), name=hour))
    # make sure the shape is as promissed:
    data = data.reindex({hour: list(times)})
    if as_dataset:
        return ff.xarray.array_to_dataset(data, hour, rename={i: 't%02d' % i for i in times})
    return data


def from_hours(data, dim='time', hour='hour', **kwargs):
    """ Combine separate times to one datetime axis

    Args:
        data (DataArray): Inputdata
        dim (str): datetime dimension
        hour (str): time dimension
        **kwargs:

    Returns:
        DataArray : combined datetime axis DataArray
    """
    import pandas as pd
    from xarray import DataArray, Dataset, concat

    if not isinstance(data, (DataArray, Dataset)):
        raise ValueError('Requires a xarray DataArray, Dataset', type(data))

    if hour not in data.dims:
        raise ValueError('Requires an hour dimension', hour)

    if dim not in data.dims:
        raise ValueError('Requires a datetime dimension', dim)

    data = data.copy()
    data = dict(data.groupby(hour))
    for ikey in data.keys():
        data[ikey] = data[ikey].assign_coords({dim: data[ikey][dim].to_index() + pd.DateOffset(hours=int(ikey))})

    return concat(data.values(), dim=dim).sortby(dim)


#
# Align launch times to standard sounding times
#

def align_datetime(data, dim='time', plev='plev', times=(0, 12), span=6, freq='12h', **kwargs):
    """ Standardize datetime to times per date, try to fill gaps

    Args:
        data (DataArray, Dataset): Input data
        dim (str): datetime dimension
        plev (str): pressure level dimension
        times (tuple): sounding times
        span (int): plus minus times (smaller than freq/2)
        freq (str): frequency of output times

    Returns:
        xarray.DataArray : datetime standardized DataArray

    """
    import numpy as np
    from pandas import DatetimeIndex
    from xarray import DataArray, Dataset
    from .. import fun as ff

    if not isinstance(data, (DataArray, Dataset)):
        raise ValueError('Requires a DataArray or Dataset', type(data))

    if dim not in data.dims:
        raise ValueError('Requires a datetime dimension', dim)

    if int(24 / (len(times) * 2)) < span:
        raise ValueError("Times and span do not agree!?", times, span)

    if int(24 / int(freq[:-1])) != len(times):
        raise ValueError("Times and freq do not match:", times, freq)

    if span > int(freq[:-1]) // 2:
        raise ValueError("Frequency and Span need to be consistent (span < freq/2): ", freq, span)

    dates = data[dim].values.copy()
    #
    # Count levels per date
    #
    _fix_datetime = np.vectorize(ff.cal.fix_datetime)
    newdates = _fix_datetime(dates, span=span)  # (time: 33%)
    resolution = np.zeros(newdates.size)
    #
    # check for duplicates in standard launch times
    #
    u, c = np.unique(newdates, return_counts=True)
    conflicts = u[c > 1]
    if conflicts.size > 0:
        counts = _count_data(data, dim=dim, plev=plev)
        ff.message("Conflicts", conflicts.size, newdates.size, **kwargs)
        for i in conflicts:
            indices = np.where(newdates == i)[0]  # numbers  (time: 45%)
            #
            # Count available data (DataArray or Dataset)
            #
            # slow
            # counts = data.isel(**{dim: indices}).count(plev).values
            # counts = _count_data(data.isel(**{dim: indices}), dim=dim, plev=plev)
            # slow end
            icounts = counts[indices]
            #
            # offsets to standard launch time
            #
            offset = np.abs((dates[indices] - i) / np.timedelta64(1, 'h'))
            j = np.argsort(offset)  # sort time offsets (first we want)
            jmax = np.argmax(icounts[j])  # check if counts from other time is larger
            if jmax != 0:
                #
                # there is a sounding with more level data (+/- 1 hour)
                #
                if (offset[j][0] + 1) <= offset[j][jmax]:
                    # ok close enough
                    jj = j.copy()
                    jj[j == 0] = jmax  # first pos is now at the position of the maximum
                    jj[j == jmax] = 0  # maximum is now first
                    j = jj
                #
                # there is a sounding with + 2 more levels
                #
                elif (icounts[j][0] + 2) <= icounts[j][jmax]:
                    # a lot more
                    jj = j.copy()
                    jj[j == 0] = jmax  # first pos is now at the position of the maximum
                    jj[j == jmax] = 0  # maximum is now first
                    j = jj
                else:
                    pass  # keep time sorting

            for m, k in enumerate(offset[j]):
                if m == 0:
                    continue  # this is the minimum

                # change back the others or add a delay to remove duplicates
                if k == 0:
                    newdates[indices[j][m]] += np.timedelta64(1, 'h')  # add offset
                    resolution[indices[j][m]] = 1  # add hour
                else:
                    newdates[indices[j][m]] = dates[indices[j][m]]  # revert back
                    resolution[indices[j][m]] = -1  # revert back
    #
    # recheck for standard times
    #
    idx_std = DatetimeIndex(newdates).hour.isin(times)
    u, c = np.unique(newdates[idx_std], return_counts=True)  # check only standard times
    conflicts = u[c > 1]
    ff.message("Conflicts remain:", conflicts.size, idx_std.sum(), newdates.size)
    #
    # new dates / new object
    #
    data = data.assign_coords({dim: newdates})
    #
    # delay
    #
    nn = (resolution > 0).sum()
    nx = (~idx_std).sum()
    data['delay'] = (dim, ((dates - newdates) / np.timedelta64(1, 'h')).astype(int))  # new coordinate for delays
    data.attrs['std_times'] = str(times)
    data['delay'].attrs['updated'] = nn
    data['delay'].attrs['missing'] = data['delay'].isnull().sum().values
    data['delay'].attrs['times'] = str(times)
    data['flag_stdtime'] = (dim, resolution.astype(int))
    data['flag_stdtime'].attrs.update({'units': '1', 'standard_name': 'flag_standard_time_conflict_resolution',
                                       'info': '0: preferred, -1: lesser candidate, 1: duplicate, less data'})
    ff.message('Updated [', nn, "] No Standard [", nx, "] [", newdates.size, "]", **kwargs)

    if not all(data[dim].values == np.sort(data[dim].values)):
        ff.message("Sorting by", dim, **kwargs)
        data = data.sortby(dim)
    return data


def _count_data(data, dim='time', plev='plev'):
    from xarray import DataArray
    #
    # Count data per pressure level (if it is a dimension)
    #
    if plev in data.dims:
        #
        # Array format
        #
        if isinstance(data, DataArray):
            return data.count(plev).values
        else:
            return data.count(plev).to_dataframe().sum(axis=1).values  # sum across variables

    elif data[dim].to_index().is_unique:
        #
        # has not pressure levels
        #
        if isinstance(data, DataArray):
            return data.count(dim).values
        else:
            return data.count(dim).to_dataframe().sum(axis=1).values
    else:
        #
        # Table format
        #
        return data.groupby(dim).count().to_dataframe().max(axis=1).values


def datetime_range(dates, times=(0, 12), freq='12h', span=6, **kwargs):
    import numpy as np
    from pandas import DatetimeIndex, Timestamp, date_range
    from .. import fun as ff
    #
    # minimum and maximum
    #
    min_date = Timestamp(dates.min()).replace(hour=np.min(times), minute=0, second=0)
    max_date = Timestamp(ff.cal.fix_datetime(dates.max(), span=span)).replace(hour=np.max(times), minute=0, second=0)
    ff.message("Dates: ", min_date, " to ", max_date, "with non-standard:",
               (~DatetimeIndex(dates).hour.isin(times)).sum(), "of",
               dates.size, **kwargs)
    #
    # full size datetime range
    #
    return date_range(min_date, max_date, freq=freq)


# def align_to_stdtimes_dataframe():
# is_table = False
# if plev not in data.dims:
#     is_table = True
#     dates, indices = np.unique(dates, return_inverse=True)
# if reindex:
#     alldates = datetime_range(dates, times=times, freq=freq, span=span, **kwargs)
#     message("New Index:", alldates[0], "to", alldates[1], "with", freq, alldates.size, **kwargs)

# resolution = np.zeros_like(newdates)
#
# for idate in conflicts:
#     #
#     # indices of conflicting new dates
#     #
#     indices = np.where(newdates == idate)[0]
#     #
#     # number of old dates ?
#     #
#     old, idx = np.unique(dates[indices], return_index=True)
#     if old.size > 1:
#         print("Conflict", idate, old.size)
#         #
#         # multiple dates
#         # counts per date
#         counts = data.isel({dim: indices}).groupby(dim).count().to_dataframe().sum(axis=1).values
#         #
#         # offsets per date
#         #
#         offset = (np.abs((dates[indices] - newdates[indices])) / np.timedelta64(1, 'h'))[idx]
#         j = np.argsort(offset)
#         jmax = np.argmax(counts[j])
#         if jmax != 0:
#             # there is a sounding with more level data
#             if (offset[j][0] + 1) <= offset[j][jmax]:
#                 # ok close enough
#                 jj = j.copy()
#                 jj[j == 0] = jmax  # first pos is now at the position of the maximum
#                 jj[j == jmax] = 0  # maximum is now first
#                 j = jj
#             elif (counts[j][0] + 2) <= counts[j][jmax]:
#                 # a lot more
#                 jj = j.copy()
#                 jj[j == 0] = jmax  # first pos is now at the position of the maximum
#                 jj[j == jmax] = 0  # maximum is now first
#                 j = jj
#             else:
#                 pass  # keep time sorting
#
#         for m, k in enumerate(offset[j]):
#             if m == 0:
#                 continue  # this is the minimum
#             #
#             # get index for newdates
#             #
#             ondices = np.where(dates == old[j][m])
#             # change back the others or add a delay to remove duplicates
#             if k == 0:
#                 newdates[ondices] += np.timedelta64(1, 'h')  # add offset
#                 resolution[ondices] = 1  # add hour
#             else:
#                 newdates[ondices] = dates[indices[j][m]]  # revert back
#                 resolution[ondices] = -1  # revert back


# def _standard_sounding_times_index(dates, numlev, span=6, **kwargs):
#     import numpy as np
#     from ..fun import fix_datetime
#     #
#     # Convert to Standard launch times
#     #
#     _fix_datetime = np.vectorize(fix_datetime)
#     newdates = _fix_datetime(dates, span=span)
#     resol = np.zeros_like(newdates)
#     #
#     # check for duplicates in standard launch times
#     #
#
#     u, c = np.unique(newdates, return_counts=True)
#     conflicts = u[c > 1]
#     if conflicts.size > 0:
#         print("Conflicts: ", conflicts.size)
#         newdates, resol = _solve_conflict_2d(conflicts, newdates, dates, data, dim=dim, plev=plev, **kwargs)
#     # else:
#     #     #
#     #     # combined unique
#     #     #
#     #     itx = (dates != newdates)  # something changed
#     #     conflicts = np.unique(newdates[itx])  # only these dates
#     #     #
#     #     # todo check if there are any conflicts?
#     #     #
#     #     newdates, resol = _solve_conflict_table(conflicts, newdates, dates, data, dim=dim, **kwargs)
#
#     return newdates, resol

#
# def _solve_table_duplicated(data, plev='plev'):
#     #
#     # can not disentangle profiles anymore / median
#     #
#     return data.groupby(plev).median()


#
# OLD VERSION
#

def _standard_sounding_times(data, dim='time', times=(0, 12), span=6, freq='12h', return_indices=False,
                             fillin=False, maxcount=False, **kwargs):
    """ Standardize datetime to times per date, try to fill gaps

    Args:
        data (xarray.DataArray): Input DataArray
        dim (str): datetime dimension
        times (tuple): sounding times
        span (int): plus minus times (smaller than freq/2)
        freq (str): frequency of output times
        return_indices (bool): return indices for alignment
        fillin (bool): use data to fill gaps
        maxcount (bool): use maximum information profiles for conflicts

    Returns:
        xarray.DataArray : datetime standardized DataArray
    """
    import numpy as np
    import pandas as pd
    from xarray import DataArray
    from .. import fun as ff

    if not isinstance(data, DataArray):
        raise ValueError('Requires a DataArray', type(data))

    if dim not in data.dims:
        raise ValueError('Requires a datetime dimension', dim)

    if int(24 / (len(times) * 2)) < span:
        raise ValueError("Times and span do not agree!?", times, span)

    if int(24 / int(freq[:-1])) != len(times):
        raise ValueError("Times and freq do not match:", times, freq)

    if span > int(freq[:-1]) // 2:
        raise ValueError("Frequency and Span need to be consistent (span < freq/2): ", freq, span)

    dates = data[dim].values.copy()
    #
    #  all possible dates
    min_date = pd.Timestamp(dates.min()).replace(hour=np.min(times), minute=0, second=0)
    max_date = pd.Timestamp(ff.cal.fix_datetime(dates.max(), span=span)).replace(hour=np.max(times), minute=0, second=0)
    ff.message("Dates: ", min_date, " to ", max_date, "with", freq, dates.size, **kwargs)
    alldates = pd.date_range(min_date, max_date, freq=freq)
    ff.message("New Index:", alldates[0], "to", alldates[1], "with", freq, alldates.size, **kwargs)
    #
    # use new dates and data that fits these dates:
    # complete reindex to new dates (implicit copy)
    #
    new = data.reindex(**{dim: alldates})
    new['delay'] = (dim, np.zeros(alldates.size))  # new coordinate for delays
    #
    # matching dates (new in old)
    #
    old_logic = np.isin(dates, new[dim].values)
    #
    # old dates not matched to new ones
    #
    jtx = np.where(~old_logic)[0]
    new['delay'].values[~np.isin(new[dim].values, dates)] = np.nan
    ff.message("Standard [", old_logic.sum(), "] Between [", (~old_logic).sum(), "]", **kwargs)
    axis = data.dims.index(dim)
    #
    # Mask
    #
    # mask = np.full_like(data.values, False)
    # mindex = idx2shp(old_logic, axis, data.values.shape)
    # mask[mindex] = True   # same dates
    # mask = (mask & np.isfinite(data.values))  # same dates & finite
    nn = 0
    nx = 0
    indices = []
    if jtx.sum() > 0:
        # _fix_datetime = np.vectorize(fix_datetime)
        #
        # Is there some dataset that fits within the given time window to new dates?
        #
        for m, i in enumerate(jtx):
            #
            # make standard date ->
            #
            # idate = _fix_datetime(dates[i], span=span)
            idate = ff.cal.fix_datetime(dates[i], span=span)
            # j = np.where(alldates == idate)[0]
            # if len(j) == 0:
            if idate not in alldates:
                nx += 1
                continue

            j = alldates.get_loc(idate)
            #
            # Indices of old and new arrays
            #
            # newindex = idx2shp(j[0], axis, new.values.shape)
            newindex = ff.xarray.idx2shp(j, axis, new.values.shape)
            oldindex = ff.xarray.idx2shp(i, axis, new.values.shape)
            status = False
            ny = np.isfinite(new.values[newindex]).sum()
            #
            # What to do with the data ?
            #
            if fillin:
                #
                # Fill in whenever missing values
                # should be the maximum possible values
                #
                logic = (~np.isfinite(new.values[newindex]) & np.isfinite(data.values[oldindex]))
                new.values[newindex] = np.where(logic, data.values[oldindex], new.values[newindex])
                # mask[oldindex] = logic
                # mask[] todo unset i in mask (remove values from standard time)
                status = True
            elif maxcount:
                #
                # Use Profile with more data
                # should be the most complete profile version
                a = np.isfinite(new.values[newindex]).sum()
                b = np.isfinite(data.values[oldindex]).sum()
                if a < b:
                    new.values[newindex] = data.values[oldindex]
                    # mask[oldindex] = np.isfinite(data.values[oldindex])
                    status = True

            elif np.isfinite(new.values[newindex]).sum() == 0:
                #
                # when all are missing (whole profile) in the current
                #
                new.values[newindex] = data.values[oldindex]
                status = True
            else:
                pass

            #
            # still missing (no report necessary)
            #
            if np.isfinite(new.values[newindex]).sum() == 0:
                continue
            #
            # keep original sounding times
            #
            if status:
                diff = (dates[i] - idate) / np.timedelta64(1, 'h')
                new['delay'].values[j] = -1 * diff  # datetime of minimum
                if kwargs.get('verbose', 0) > 2:
                    ff.message(m, dates[i], "+", -1 * diff, " = ", idate, " [",
                               np.sum(np.isfinite(data.values[oldindex]), axis=0), ">",
                               ny, "]", **kwargs)

                # indices += [(i, j[0])]
                indices += [(i, j)]
                nn += 1
            else:
                # nothing could be done
                if kwargs.get('verbose', 0) > 1:
                    ff.message(m, dates[i], "[", np.sum(np.isfinite(data.values[oldindex]), axis=0), "] Passed [",
                               np.sum(np.isfinite(new.values[newindex]), axis=0), "]", idate,
                               **kwargs)
    #
    # Add delay coordinate
    #
    ff.message('Updated [', nn, "] No Standard [", nx, "] [", jtx.size, "]", **kwargs)
    new.attrs['std_times'] = str(times)
    new['delay'].attrs['updated'] = nn
    new['delay'].attrs['missing'] = new['delay'].isnull().sum().values
    new['delay'].attrs['times'] = str(times)
    #
    # return indices for multi array alignment (datasets)
    #
    if return_indices:
        return new, np.array(indices)
    return new
