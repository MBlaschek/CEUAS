# -*- coding: utf-8 -*-

__all__ = ['var', 'winds', 'boxplot', 'bars']


def var(data, dim='plev', ax=None, logy=False, yticklabels=None, showna=False, **kwargs):
    import numpy as np
    from xarray import DataArray
    from ._helpers import line, set_labels, get_info
    from ..fun import message

    if not isinstance(data, DataArray):
        raise ValueError('Requires a DataArray', type(data))

    if dim not in data.dims:
        raise ValueError('Requires a datetime dimension', dim)

    if data.ndim > 1:
        raise ValueError('Too many dimensions', data.dims, data.shape)

    values = data.values
    levels = data[dim].values.copy()
    lev_units = data[dim].attrs.get('units', 'Pa')
    if lev_units == 'Pa':
        levels = levels.astype(float) / 100.
        message('Converting', lev_units, 'to', 'hPa', levels, **kwargs)
        lev_units = 'hPa'
        if yticklabels is not None:
            yticklabels = np.int_(np.asarray(yticklabels) / 100)

    itx = np.isfinite(values)

    kwargs.update({'marker': kwargs.get('marker', 'o')})

    set_labels(kwargs, xlabel=get_info(data),
               title=get_info(data), ylabel=dim + ' [%s]' %lev_units)
    ax = line(values[itx], levels[itx], ax=ax, **kwargs)
    if np.sum(itx) != np.size(levels) and showna:
        itmp = ax.get_xlim()
        ax.plot([itmp[1]]*np.sum(~itx), levels[~itx], marker=kwargs.get('marker','x'), c='red')

    if logy:
        ax.set_yscale('log')

    if np.diff(ax.get_ylim())[0] > 0:
        ax.invert_yaxis()

    ax.set_yticks(levels)
    if yticklabels is not None:
        yticklabels = np.asarray(yticklabels)  # can not calc on list
        ax.set_yticks(yticklabels)
        ax.set_yticklabels(np.int_(yticklabels))
    else:
        ax.set_yticks(levels[::2])
        ax.set_yticklabels(np.int_(levels[::2]))

    ax.set_ylim(*kwargs.get('ylim', (None, None)))  # fixed
    ax.set_xlim(*kwargs.get('xlim', (None, None)))
    return ax


def winds(data, u='u', v='v', dim='plev', barbs=True, ax=None, logy=False, yticklabels=None, showna=True, **kwargs):
    import numpy as np
    from xarray import Dataset
    import matplotlib.pyplot as plt
    from ._helpers import set_labels, line
    from ..fun import message

    if not isinstance(data, Dataset):
        raise ValueError('Requires a Dataset', type(data))

    if dim not in data.dims:
        raise ValueError('Requires a datetime dimension', dim)

    if u not in data.data_vars:
        raise ValueError('Requires a u-wind component', u)

    if v not in data.data_vars:
        raise ValueError('Requires a v-wind component', v)

    if data[u].ndim > 1 or data[v].ndim > 1:
        raise ValueError('Too many dimensions', data.dims, data.shape)

    uvalues = data[u].values
    vvalues = data[v].values
    levels = data[dim].values.copy()
    lev_units = data[dim].attrs.get('units', 'Pa')
    if lev_units == 'Pa':
        levels = levels.astype(float) / 100.
        message('Converting', lev_units, 'to', 'hPa', levels, **kwargs)
        lev_units = 'hPa'
        if yticklabels is not None:
            yticklabels = np.int_(np.asarray(yticklabels) / 100)

    itx = np.isfinite(uvalues) & np.isfinite(vvalues)

    set_labels(kwargs, xlabel='Winds ['+data[u].attrs.get('units','m/s')+']',
               title='Winds', ylabel=dim + ' [%s]' %lev_units)

    if barbs:
        if ax is None:
            f, ax = plt.subplots()  # 1D SNHT PLOT
        speed = np.sqrt(uvalues * uvalues + vvalues * vvalues)
        ax.barbs(np.zeros_like(levels[itx]), levels[itx], uvalues[itx], vvalues[itx], speed[itx],
                 alpha=kwargs.get('alpha', 1))
    else:
        ax = line(uvalues[itx], levels[itx], label='u-wind', ax=ax, **kwargs)
        ax = line(vvalues[itx], levels[itx], label='v-wind', ax=ax, **kwargs)
        ax.legend()

    ax.grid('gray', ls='--')
    ax.set_title(kwargs.get('title'))
    ax.set_ylabel(kwargs.get('ylabel'))
    ax.set_xlabel(kwargs.get('xlabel'))

    if logy:
        ax.set_yscale('log')

    if np.diff(ax.get_ylim())[0] > 0:
        ax.invert_yaxis()

    ax.set_yticks(levels, minor=True)
    if yticklabels is not None:
        yticklabels = np.asarray(yticklabels)  # can not calc on list
        ax.set_yticks(yticklabels)
        ax.set_yticklabels(np.int_(yticklabels))
    else:
        ax.set_yticks(levels[::2])
        ax.set_yticklabels(np.int_(levels[::2]))
    ax.set_ylim(*kwargs.get('ylim', (None, None)))  # fixed
    ax.set_xlim(*kwargs.get('xlim', (None, None)))
    return ax


def boxplot(data, dim='plev', ax=None, vline=None, yticklabels=None, logy=False, **kwargs):
    import numpy as np
    import pandas as pd
    from xarray import DataArray
    import matplotlib.pyplot as plt
    from ._helpers import set_labels, get_info

    if not isinstance(data, DataArray):
        raise ValueError('Requires a DataArray', type(data))

    if dim not in data.dims:
        raise ValueError('Requires a level dimension', dim)

    if data.ndim != 2:
        raise ValueError('Too many/few dimensions', data.dims, data.shape)

    if ax is None:
        fig, ax = plt.subplots()

    axis = data.dims.index(dim)
    dims = list(data.dims)
    dims.remove(dim)
    odim = dims[0]
    levels = data[dim].values.copy()
    lev_units = data[dim].attrs.get('units', 'Pa')
    if lev_units == 'Pa':
        levels = levels.astype(float) / 100.
        lev_units = 'hPa'
        if yticklabels is not None:
            yticklabels = np.int_(np.asarray(yticklabels) / 100)

    levels = levels.astype(int)
    if axis == 0:
        idata = pd.DataFrame(data.values.T, index=data[odim].values, columns=levels)
    else:
        idata = pd.DataFrame(data.values, index=data[odim].values, columns=levels)

    set_labels(kwargs, xlabel=get_info(data),
               title=get_info(data), ylabel=dim + ' [%s]' % lev_units)
    idata = idata.sort_index(axis=1, ascending=False)
    idata.boxplot(ax=ax, vert=False, return_type='axes', sym='+')
    if vline is not None:
        ax.axvline(x=vline, color='k', lw=1)

    ax.grid(ls='--')
    ax.set_title(kwargs.get('title'))
    ax.set_ylabel(kwargs.get('ylabel'))
    ax.set_xlabel(kwargs.get('xlabel'))
    if logy:
        ax.set_yscale('log')

    if yticklabels is not None:
        yticklabels = np.asarray(yticklabels)  # can not calc on list
        ax.set_yticks(yticklabels)
        ax.set_yticklabels(np.int_(yticklabels))
        # for label in ax.yaxis.get_ticklabels():
        #     if int(label.get_text()) not in yticklabels:
        #         label.set_visible(False)
    # else:
    #     for label in ax.yaxis.get_ticklabels()[::2]:
    #         label.set_visible(False)
    ax.set_ylim(*kwargs.get('ylim', (None, None)))  # fixed
    ax.set_xlim(*kwargs.get('xlim', (None, None)))
    return ax


def bars(data, dim='plev', ax=None, vline=None, yticklabels=None, logy=False, use_levels=False, bar_kwargs={}, **kwargs):
    import numpy as np
    from xarray import DataArray
    import matplotlib.pyplot as plt
    from ._helpers import set_labels, get_info

    if not isinstance(data, DataArray):
        raise ValueError('Requires a DataArray', type(data))

    if dim not in data.dims:
        raise ValueError('Requires a level dimension', dim)

    if data.ndim != 1:
        raise ValueError('Too many/few dimensions', data.dims, data.shape)

    if ax is None:
        fig, ax = plt.subplots()

    levels = data[dim].values.copy()
    lev_units = data[dim].attrs.get('units', 'Pa')
    if lev_units == 'Pa':
        levels = levels.astype(float) / 100.
        lev_units = 'hPa'

    levels = levels.astype(int)
    set_labels(kwargs, xlabel=get_info(data),
               title=get_info(data), ylabel=dim + ' [%s]' % lev_units)
    if use_levels:
        ax.barh(levels, data.values, align='center', **bar_kwargs)
        # ax.set_yticklabels([str(i) for i in levels])
    else:
        ax.barh(np.arange(1, levels.size+1), data.values, align='center', **bar_kwargs)
        ax.set_yticklabels([str(i) for i in levels])

    if logy:
        ax.set_yscale('log')

    if np.diff(levels)[0] > 0:
        ax.invert_yaxis()

    if vline is not None:
        ax.axvline(x=vline, color='k', lw=1)

    ax.grid(ls='--')
    ax.set_title(kwargs.get('title'))
    ax.set_ylabel(kwargs.get('ylabel'))
    ax.set_xlabel(kwargs.get('xlabel'))

    if yticklabels is not None:
        for label in ax.yaxis.get_ticklabels():
            if int(label.get_text()) not in yticklabels:
                label.set_visible(False)
    # else:
    #     for label in ax.yaxis.get_ticklabels()[::2]:
    #         label.set_visible(False)
    ax.set_ylim(*kwargs.get('ylim', (None, None)))  # fixed
    ax.set_xlim(*kwargs.get('xlim', (None, None)))
    return ax
