# -*- coding: utf-8 -*-


def histogram(data, y=None, bins=36, ax=None, quantiles=None, bottom=2, **kwargs):
    import numpy as np
    import pandas as pd
    from xarray import DataArray
    from ._helpers import set_labels, get_info
    import matplotlib.pyplot as plt

    if not isinstance(data, DataArray):
        raise ValueError('Requires a DataArray', type(data))

    if data.ndim > 1:
        raise ValueError('Too many dimensions', data.dims, data.shape)

    if ax is None:
        f = plt.figure(figsize=kwargs.get('figsize', None))
        ax = plt.subplot(111,polar=True)  # 1D SNHT PLOT

    # width of each bin on the plot
    width = (2*np.pi) / bins
    values = data.values
    itx = np.isfinite(values)
    values = values[itx]
    # make the histogram that bined on 24 hour
    counts, bins = np.histogram(values, bins=bins)
    if y is not None:
        if quantiles is None:
            quantiles = [0.05, 0.25, 0.5, 0.75, 0.95]
        idx = np.digitize(values, bins, right=True)  # indices
        yv = y.values[itx]
        yv = pd.DataFrame({'y': yv, 'groups': idx}).groupby('groups').quantile(quantiles).unstack(-1)
        for i in np.sort(quantiles)[::-1]:
            ax.bar(bins * np.pi / 180., yv.xs(i, axis=1, level=1).y.values, width=width, bottom=bottom,
                   label='Q%02d' % int(i*100))
        ax.legend(loc='right', bbox_to_anchor=(1.4,0.5))
        set_labels(kwargs, xlabel=get_info(y), title=get_info(y))
    else:
        # make a polar plot
        bar = ax.bar(bins[:-1] * np.pi / 180.,counts, width=width, bottom=bottom)
        set_labels(kwargs, xlabel=get_info(data), title=get_info(data))

    ax.set_title(kwargs.get('title', None))
    ax.set_ylabel(kwargs.get('ylabel', None))
    # set the lable go clockwise and start from the top
    ax.set_theta_zero_location("N")
    # clockwise
    ax.set_theta_direction(-1)
    return ax
