# -*- coding: utf-8 -*-

__all__ = ["xarray_function_wrapper", "array_to_dataset", "update_shape", 'combine_datasets', 'table_to_dataset']


def xarray_function_wrapper(x, wfunc=None, add_dim=None, y=None, **kwargs):
    """ Map a numpy function that is not currently in xarray to xarray with apply_ufunc

    Args:
        x (DataArray): Input dataset
        y (DataArray): second Input dataset
        wfunc (callable): function to call, e.g.: np.nanpercentile
        add_dim (str): new dimension for output
        **kwargs: all arguments to function

    Keyword Args:
        dim (str): Dimension to remove
        axis (int): axis for numpy functions
        debug (bool): show debug information of call

    Returns:
        xarray.DataArray : result of function call retains attrs

    Examples:
        Aggregate some data and remove dimensions
        >>> import numpy as np; import pandas as pd; import xarray as xr
        >>> def myfunc(x, **kwargs):
        >>>     return np.isfinite(x).sum(**kwargs)
        >>> dataset = xr.DataArray(np.random.randn(1000,2), dims=('time','lev'), coords=[pd.date_range('1-1-2019', periods=1000), [10, 12]])
        >>> xarray_function_wrapper(dataset, wfunc=myfunc, dim='time', axis=0)

        Calculate Percentiles and add a dimension (perc) for 9 percentiles
        >>> xarray_function_wrapper(dataset, wfunc=rasotools.fun.nanfunc, add_dim='perc', dim='time', ffunc=np.nanpercentile, fargs=((np.arange(5,95,10),)))
    """
    import xarray as xr
    if not isinstance(x, xr.DataArray):
        raise ValueError('requires a DataArray')

    jdims = list(x.dims)
    if 'dim' in kwargs.keys():
        jdims.remove(kwargs.pop('dim'))
        # if 'axis' not in kwargs.keys():
        #     raise RuntimeWarning('axis keyword not present')
        #     kwargs['axis'] = x.dims.index(kwargs['dim'])   # add axis keyword

    if kwargs.pop('debug', False):
        print(x.dims, x.shape, wfunc)

    if add_dim is not None:
        jdims = [add_dim] + jdims

    if y is not None:
        return xr.apply_ufunc(wfunc, x, y, kwargs=kwargs,
                              input_core_dims=[x.dims, y.dims], output_core_dims=[jdims],
                              keep_attrs=True)
    else:
        return xr.apply_ufunc(wfunc, x, kwargs=kwargs,
                              input_core_dims=[x.dims], output_core_dims=[jdims],
                              keep_attrs=True)


def set_attrs(data, name, set='', add='', default=''):
    from xarray import Dataset, DataArray
    if not isinstance(data, (Dataset, DataArray)):
        raise ValueError("Requires an XArray Dataset or Array")

    if name not in data.attrs:
        data.attrs[name] = default
    elif set != '':
        data.attrs[name] = set
    else:
        data.attrs[name] += add


def idx2shp(idx, axis, shape):
    index = [slice(None)] * len(shape)
    index[axis] = idx
    return tuple(index)


def limits(data, var, raw):
    import numpy as np
    from xarray import Dataset
    if not isinstance(data, Dataset):
        raise ValueError("Requires an Xarray Dataset", type(data))
    with np.errstate(invalid='ignore'):
        #
        # values below 0 or missing -> reset to original
        #
        logic = ((data[var].values < 0) | np.isnan(data[var]))
        data[var].values = np.where(logic,
                                    data[raw].values,
                                    data[var].values)


def update_shape(data, order):
    from xarray import DataArray
    if not isinstance(data, DataArray):
        raise ValueError("Requires an xarray DataArray", type(data))

    iorder = ()
    for idim in order:
        if idim in data.dims:
            iorder += (idim,)

    if len(iorder) > 0:
        return data.transpose(*iorder)
    return data


def array_to_dataset(data, dim, rename=None):
    """ Convert a DataArray to Dataset and consider dependent coordinates

    Args:
        data (DataArray): Input DataArray
        dim (str): Coordinate to use as variables
        rename (dict): renaming policy

    Returns:
        Dataset :
    """
    from xarray import DataArray

    if not isinstance(data, DataArray):
        raise ValueError('Requires a DataArray', type(data))

    if dim not in data.dims:
        raise ValueError('Requires a datetime dimension', dim)

    data = data.copy()
    tmp = {}

    for i, j in data.coords.items():
        if dim in j.dims and i != dim:
            tmp[i] = data[i]
            data = data.drop(i)

    data = data.to_dataset(dim=dim)
    for i, j in tmp.items():
        data[i] = j  # add as Coords / Variables

    if rename is not None:
        return data.rename(rename)

    return data


def table_to_dataset(data, dim='time', plev='plev', levels=None, return_rejected=False, **kwargs):
    """ Convert pandas Dataframe to xarray Dataset

    Args:
        data (dataframe): input dataframe (columns are variables)
        dim (str): datetime dimension
        plev (str): pressure dimension
        levels (list): pressure levels to consider
        **kwargs:

    Returns:
        Dataset : 2d (datetime x pressure levels) x variables
    """
    from .. import config
    from xarray import Dataset

    if levels is None:
        levels = config.std_plevels

    # dimensions for output
    varis = [dim, plev]
    attrs = None
    if isinstance(data, Dataset):
        # copy / backup attributes
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
    data = data[data[plev].isin(levels)]
    #
    # convert to xarray
    #
    data = data.reset_index().set_index(varis)
    if not data.index.is_unique:
        if kwargs.get('verbose', 0) > 0: print("Non-unique index, removing duplicates...")
        data = data.loc[~data.index.duplicated()]  # remove duplicated

    data = data.to_xarray()  # 1D -> 2D
    if attrs is not None:
        # add attributes again
        for i, j in vatt.items():
            data[i].attrs.update(j)
        data.attrs.update(attrs)
        data[dim].attrs.update(tatts)

    return data


def dataset_to_table(data, variables=None, dim='time'):
    if variables is None:
        variables = list(data.data_vars)
    
    attrs = {i: data[i].attrs.copy() for i in variables}
    
    dims = list(data.dims)
    dims.remove(dim)
    dims = [dim] + dims
    print(dims, variables)
    data = data[variables] \
        .to_dataframe() \
        .reset_index() \
        .dropna(subset=dims) \
        .sort_values(dims) \
        .reset_index(drop=True) \
        .to_xarray()
    # put metadata back
    for i in variables:
        data[i].attrs.update(attrs[i])
    return data


def combine_datasets(a, b, subset=None, suffix=None, only2a=True, profiles=True, plev='plev', add_flags=True, **kwargs):
    """ Combine two datasets

    Args:
        a (Dataset, DataArray): dataset 1
        b (Dataset, DataArray): dataset 2
        subset (list): list of variables
        suffix (str): renaming
        only2a (bool): fill only to dataset 1

        **kwargs:

    Returns:
        Dataset, Dataset : filled dataset, flags
        DataArray, DataArray : filled dataset, flags
        
    Note:
        similar to combine_nested(), however combined vars and can be profiles only
    """
    import numpy as np
    from xarray import Dataset, DataArray, broadcast, align
    from . import message

    if not isinstance(a, (Dataset, DataArray)):
        raise ValueError("Requires a xarray Dataset or DataArray")

    if not isinstance(b, (Dataset, DataArray)):
        raise ValueError("Requires a xarray Dataset or DataArray")

    avars = []
    bvars = []
    loop = False
    if isinstance(a, DataArray) and not isinstance(b, DataArray):
        if a.name not in b.data_vars:
            raise RuntimeError("No matching variables:", a.name, list(b.data_vars))
        avars = [a]
        bvars = [b[a.name]]
        loop = False

    elif isinstance(a, Dataset) and not isinstance(b, Dataset):
        if b.name not in a.data_vars:
            raise RuntimeError("No matching variables:", b.name, list(a.data_vars))

        avars = [a[b.name]]
        bvars = [b]
        loop = False

    elif isinstance(a, Dataset) and isinstance(b, Dataset):
        if not any(np.isin(a.data_vars, b.data_vars)):
            raise RuntimeError("No matching variables found:", list(a.data_vars), list(b.data_vars))
        if subset is None:
            subset = np.asarray(a.data_vars)[np.isin(a.data_vars, b.data_vars)]
        else:
            tmp = np.asarray(a.data_vars)[np.isin(a.data_vars, b.data_vars)]
            subset = subset[np.isin(subset, tmp)]
        avars = list(dict(a[subset]).values())
        bvars = list(dict(b[subset]).values())
        loop = True if len(avars) > 1 else False
    else:
        # DataArray + DataArray
        if a.name != b.name:
            message("Warning Names do not match:", a.name, b.name, **kwargs)

    if loop:
        # lots of variables to process (Datasets)
        res = {}
        for ivar, jvar in zip(avars, bvars):
            i, j = combine_datasets(ivar, jvar, suffix=suffix, only2a=only2a, profiles=profiles, **kwargs)
            res[i.name] = i
            if add_flags:
                res[j.name] = j
        return Dataset(res)
    # Only one Variable to merge (DataArray)
    # Check dimensions
    if tuple(a.dims) != tuple(b.dims):
        # something needs to be done
        if len(a.dims) != len(b.dims):
            raise ValueError("DataArrays do not have matching dimensions:", a.dims, b.dims)

        message("Updating shape: ", b.dims, "to", a.dims, **kwargs)
        b = update_shape(b, tuple(a.dims))
    #
    # Make same length
    #
    if only2a:
        message("Aligning: ", a.shape, b.shape, **kwargs)
        a, b = align(a, b, join='left')
    else:
        a, b = broadcast(a, b)
        message("Broadcasting: ", a.shape, b.shape, **kwargs)
    #
    # Combine data
    #
    a = a.copy()
    if profiles:
        # check if one profile has more data than the other one
        # todo update meaings
        # 0,1,2 [nan, A, B]
        flags = (a.count(plev) > b.count(plev))
        a = a.where(flags, b)
        flags = flags.astype(int)
    else:
        flags = a.copy()
        logic = np.isfinite(a.values)  # full array
        # only select whole profiles, not individual values
        a.values = np.where(logic, a.values, b.values)
        # where filled in
        flags.values = (~logic) & np.isfinite(b.values)
        flags = flags.astype(int)
        
    if suffix is None:
        suffix = ''
    flags.name = flags.name + suffix + '_flag'
    flags.attrs.update({'meaning': 'from_B from_A', 'range': [0,1]})
    a.name = a.name + suffix
    return a, flags

#
# def dataset_to_hours(data, ref='infer', variables=None, dim='time', lev='plev', hour='hour', suffix='',
#                      interpolate=False, levels=None, **kwargs):
#     """ convert all DataArrays in the Dataset to hours
#
#     Args:
#         data (Dataset): xarray Dataset
#         ref (str): name of DataArray to act as target for the rest, if None only times will be selected
#         variables (list): variables to consider
#         dim (str): datetime dimension
#         lev (str): pressure level dimension
#         suffix (str): suffix for new Dataset
#         interpolate (bool): interpolate along pressure levels
#         levels (list): interpolation pressure levels
#         **kwargs
#     Returns:
#         Dataset : rearanged and standardized Dataset
#     """
#     from .. import fun as ff
#     from . import time, std, convert
#     from xarray import Dataset, concat
#     from pandas import Index
#
#     data = data.copy()
#     new = Dataset()
#     new.attrs.update(data.attrs.copy())
#     #
#     # always choose a standard
#     #
#     if ref is 'infer':
#         counts = data.count().to_array().to_series()
#         ref = counts.idxmax()
#         ff.message("Infering largest variable: ", ref, **kwargs)
#
#     if ref not in data.data_vars:
#         ref = None
#
#     if ref is not None:
#         """ need to update that routine
#         runs with one variable to make a choice based on available data
#         select closest delay (<=1 h) / most data (+2 levels)
#         return an index for the time variable to select and new datetime values (std)
#         delays can be calculated alone
#         indices can be directly applied to a DataArray or Dataset alike
#
#         """
#         new[ref + suffix], idx = std.align_datetime(data[ref], return_indices=True, dim=dim,
#                                                              **ff.levelup(**kwargs))
#         ff.message("Using Standard: ", ref, data[ref].shape, "to", new[ref + suffix].shape, "with",
#                    kwargs.get('times', [0, 12]), **kwargs)
#         new.attrs.update({'std_times': new[ref + suffix].attrs['std_times']})
#     else:
#         #
#         # Count loses
#         #
#         counts = data[dim].groupby(dim + '.hour').count().to_series()
#         times = kwargs.get('times', [0, 12])
#         other = ",".join(["%d:%d" %(i,j) for i,j in dict(counts.drop(times)).items()])
#         ff.message("Using", dim, data[dim].size, " with ", str(times), "[",counts[times].sum(),"]:",
#                 list(counts[times].values)," Other [", counts.drop(times).sum(),"]:", other, **kwargs)
#
#     for i in list(data.data_vars):
#         if i == ref:
#             continue
#
#         if variables is not None and i not in variables:
#             continue
#
#         if interpolate and lev in data[i].dims:
#             tmp = convert.vertical_interpolation(data[i], lev, levels=levels, **ff.levelup(**kwargs))
#             # tmp is a copy
#             new[i + suffix] = time.sel_hours(tmp, dim=dim, **ff.levelup(**kwargs))
#         else:
#             # tmp is a copy
#             new[i + suffix] = time.sel_hours(data[i], dim=dim, **ff.levelup(**kwargs))
#
#         if ref is not None and idx.size > 0:
#             #
#             #
#             #
#             if dim in data[i].dims:
#                 inew = [slice(None)] * new[i + suffix].ndim
#                 iold = [slice(None)] * data[i].ndim
#                 inew[new[i + suffix].dims.index(dim)] = idx[:, 1]
#                 iold[data[i].dims.index(dim)] = idx[:, 0]
#                 new[i + suffix].values[tuple(inew)] = data[i].values[tuple(iold)]  # copy same
#
#         ff.message("Converted:", i + suffix, ref is not None, int(data[i].count()), int(new[i + suffix].count()), **kwargs)
#
#     data = dict(new.groupby(dim + '.hour'))
#     for ikey in data.keys():
#         data[ikey] = data[ikey].assign_coords({dim: data[ikey][dim].to_index().to_period('D').to_timestamp().values})
#
#     data = concat(data.values(), dim=Index(data.keys(), name=hour))
#     if 'delay' in data.coords:
#         data = data.reset_coords('delay').rename({'delay': 'delay' + suffix})
#
#     return data
