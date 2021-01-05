# -*- coding: utf-8 -*-


def plot_levels(imin, m=None, p=0, n=13):
    """plotlevels with logspace
    """
    import numpy as np

    if not isinstance(imin, (int, float)):
        raise ValueError("Require a int or float")

    if m is not None and not isinstance(m, (int, float)):
        raise ValueError("Require a int or float")

    if np.log(np.abs(imin)) < 2 and p == 0:
        p += 1

    # Centered around 0
    if m is None:
        values = (-1 * np.round(imin * np.logspace(0, 2, n/2) / 100., p)).tolist() + [0] + np.round(
            imin * np.logspace(0, 2, n/2) / 100., p).tolist()

    else:
        if imin > m:
            tmp = imin
            imin = m
            m = tmp

        if imin == 0:
            # only Positive
            values = np.round(m * np.logspace(0, 2, n) / 100., p).tolist()

        elif imin < 0 and m < 0:
            # only negative
            values = np.round(imin * np.logspace(0, 2, n) / 100., p).tolist()

        else:
            # positve and negative
            values = np.round(imin * np.logspace(0, 2, n/2) / 100., p).tolist() + np.round(
                m * np.logspace(0, 2, n/2) / 100., p).tolist()

    return np.unique(np.sort(np.asarray(values))).tolist()


def plot_arange(imin, m=None, p=0, n=7):
    import numpy as np
    if not isinstance(imin, (int, float)):
        raise ValueError("Require a int or float")

    if m is not None and not isinstance(m, (int, float)):
        raise ValueError("Require a int or float")

    if np.log(np.abs(imin)) < 2 and p == 0:
        p += 1

    if m is None:
        values = np.linspace(-1 * imin, imin, n)
    else:
        values = np.linspace(imin, m, n)

    values = np.round(values, p)
    return np.unique(np.sort(np.asarray(values))).tolist()

def get_info(x):
    return x.attrs.get('standard_name', x.name if x.name is not None else 'var') + ' [' + x.attrs.get('units', '1') + ']'


def set_labels(known, **kwargs):
    for ikey, ival in kwargs.items():
        known.update({ikey: known.get(ikey, ival)})


def line(dates, values, title='', ylabel='', xlabel='', xerr=None, yerr=None, filled=False, minmax=False, ax=None, **kwargs):
    """

    Args:
        dates (ndarray): datetime
        values (ndarray): values
        title (str): title
        ylabel (str): y-axis label
        xlabel (str): x-axis label
        xerr (ndarray): x error
        yerr (ndarray): y error
        filled (bool): fill between error lines
        ax (axes): matplotlib axis
        **kwargs: optional keyword arguments for plotting

    Returns:
        axes : matplotlib axis
    """
    import matplotlib.pyplot as plt
    if ax is None:
        f, ax = plt.subplots(figsize=kwargs.get('figsize', None))  # 1D SNHT PLOT

    if xerr is None and yerr is None:
        ax.plot(dates, values,
                ls=kwargs.get('ls', '-'),
                lw=kwargs.get('lw', 1),
                label=kwargs.get('label', None),
                marker=kwargs.get('marker', None),
                alpha=kwargs.get('alpha', 1),
                color=kwargs.get('color', None),
                zorder=kwargs.get('zorder', 1))  # Line Plot
    elif filled:
        ax.plot(dates, values,
                ls=kwargs.get('ls', '-'),
                lw=kwargs.get('lw', 1),
                label=kwargs.get('label', None),
                marker=kwargs.get('marker', None),
                # alpha=kwargs.get('alpha', 1),
                color=kwargs.get('color', None))  # Line Plot
        low, high = lowhigh(dates, values, xerr=xerr, yerr=yerr, minmax=minmax)
        if xerr is None:
            ax.fill_between(dates, low, high,
                            alpha=kwargs.get('alpha', 0.5),
                            color=ax.get_lines()[-1].get_color(),
                            hatch=kwargs.get('hatch', None),
                            zorder=-1)
        else:
            ax.fill_betweenx(values, low, high,
                             alpha=kwargs.get('alpha', 0.5),
                             color=ax.get_lines()[-1].get_color(),
                             hatch=kwargs.get('hatch', None),
                             zorder=-1)
    else:
        ax.errorbar(dates, values, xerr=xerr, yerr=yerr,
                    ls=kwargs.get('ls', '-'),
                    lw=kwargs.get('lw', 1),
                    label=kwargs.get('label', None),
                    marker=kwargs.get('marker', None),
                    alpha=kwargs.get('alpha', 1),
                    color=kwargs.get('color', None),
                    zorder=kwargs.get('zorder', 1))  # Line Plot

    ax.grid('gray', ls='--')
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    return ax


def lowhigh(dates, values, xerr=None, yerr=None, minmax=False):
    import numpy as np
    if xerr is None:
        if hasattr(yerr, '__iter__') and len(np.shape(yerr)) == 2:
            if minmax:
                low = yerr[0]
                high = yerr[1]
            else:
                low = values - yerr[0]
                high = values + yerr[1]
        else:
            low = values - yerr
            high = values + yerr
    else:
        if hasattr(xerr, '__iter__') and len(np.shape(xerr)) == 2:
            if minmax:
                low = xerr[0]
                high = xerr[1]
            else:
                low = dates - xerr[0]
                high = dates + xerr[1]
        else:
            low = dates - xerr
            high = dates + xerr
    return low, high


def contour(ax, dates, plevs, test, logy=False, colorlevels=None, yticklabels=None, legend=True,
            title='', xlabel='', ylabel='', clabel='', **kwargs):
    import numpy as np
    import matplotlib.pyplot as plt

    if ax is None:
        f, ax = plt.subplots(figsize=kwargs.get('figsize', None))  # 1D SNHT PLOT

    if kwargs.get('use_pcolormesh', False):
        from matplotlib.colors import BoundaryNorm
        cmap = plt.get_cmap(kwargs.pop('cmap', 'RdYlBu_r'))
        norm = BoundaryNorm(colorlevels, ncolors=cmap.N, clip=True)
        cs = ax.pcolormesh(dates, plevs, test.T, cmap=cmap, norm=kwargs.pop('norm', norm),
                           vmin=kwargs.pop('vmin', None),
                         vmax=kwargs.pop('vmax', None))
    else:
        cs = ax.contourf(dates, plevs, test.T, levels=colorlevels,
                         cmap=kwargs.pop('cmap', 'RdYlBu_r'),
                         extend=kwargs.get('extend', 'neither'),
                         vmin=kwargs.pop('vmin', None),
                         vmax=kwargs.pop('vmax', None),
                         norm=kwargs.pop('norm', None)
                         )  # hatches=kwargs.pop('hatches', [])

    if logy:
        ax.set_yscale('log')
    # xlim auto range
    tmp = np.isfinite(test).sum(-1)
    tmp = np.where(tmp > 0)[0]
    ax.set_xlim(dates[np.min(tmp)], dates[np.max(tmp)])

    ax.set_yticks(plevs)
    if yticklabels is not None:
        yticklabels = np.asarray(yticklabels)  # can not calc on list
        ax.set_yticks(yticklabels)
        ax.set_yticklabels(np.int_(yticklabels / kwargs.get('levfactor', 100)))
    else:
        ax.set_yticks(plevs[::2])
        ax.set_yticklabels(np.int_(plevs[::2] / kwargs.get('levfactor', 100)))

    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)

    if np.diff(ax.get_ylim()) > 0:
        ax.invert_yaxis()

    ax.grid('gray', ls='--')
    if legend:
        cbar = plt.colorbar(cs, ax=ax, orientation=kwargs.get('orientation', "vertical"),
                            fraction=kwargs.get('fraction', 0.01),
                            aspect=kwargs.get('aspect', 50),
                            extend=kwargs.get('extend', 'neither'),
                            shrink=kwargs.get('shrink', 0.8),
                            ticks=kwargs.get('legend_ticks', None))
        if kwargs.get('legend_ticklabels', None) is not None:
            cbar.set_ticklabels(kwargs.get('legend_ticklabels', None))
        cbar.set_label(clabel)
    else:
        return ax, cs
    return ax


def cost(lon, lat, values):
    """ Estimate Cost between Points

    Parameters
    ----------
    lon         array/list      Longitudes
    lat         array/list      Latitudes
    values      array/list      Values

    Returns
    -------
    float   Cost
    """
    import numpy as np
    n = lon.shape[0]
    cost = np.zeros((n))
    for i in range(n):
        # Distance of all points * difference of values
        #
        cost[i] = np.nansum((distance(lon[i], lat[i], lat, lon) * (values[i] - values)) ** 2)

    return cost  # np.nansum(cost)/np.sum(np.isfinite(values))


def distance(ilon, ilat, lats, lons):
    """ Calculate Distance between one point and others

    Parameters
    ----------
    ilon
    ilat
    lats
    lons

    Returns
    -------
    array   Distances
    """
    import numpy as np
    ix = np.cos(ilat * np.pi / 180.) * np.cos(ilon * np.pi / 180.)
    iy = np.cos(ilat * np.pi / 180.) * np.sin(ilon * np.pi / 180.)
    iz = np.sin(ilat * np.pi / 180.)
    x = np.cos(lats * np.pi / 180.) * np.cos(lons * np.pi / 180.)
    y = np.cos(lats * np.pi / 180.) * np.sin(lons * np.pi / 180.)
    z = np.sin(lats * np.pi / 180.)
    dists = ix * x + iy * y + iz * z
    return np.arccos(dists * 0.999999)


def stats(data, dim='time'):
    from ..met.time import statistics
    med = data.median().values
    std = data.std().values
    rmse = statistics(data, f='rmse', dim=dim).median().values
    return "R:{:.2f} M:{:.2f} S:{:.2f}".format(rmse, med, std)


def discrete_colormap(values, cmap='jet'):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    cmap = plt.get_cmap(cmap)  # define the colormap
    # extract all colors from the .jet map
    cmaplist = [cmap(i) for i in range(cmap.N)]
    # create the new map
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'Custom cmap', cmaplist, cmap.N)

    # define the bins and normalize
    norm = mpl.colors.BoundaryNorm(values, cmap.N)
    return cmap, norm
