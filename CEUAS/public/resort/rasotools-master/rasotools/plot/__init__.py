# -*- coding: utf-8 -*-
from ._helpers import plot_arange, plot_levels, discrete_colormap
from . import map
from . import profile
from . import polar
from .time import *
from .analysis import *


def init_fig_vertical(n=2, ratios=(1, 4), figsize=None, sharex='row', **kwargs):
    import matplotlib.pyplot as plt
    gridspec_kw = {i: kwargs[i] for i in ['height_ratios', 'hspace', 'wspace'] if i in kwargs.keys()}
    gridspec_kw.update({'height_ratios': ratios})
    return plt.subplots(n, 1, sharex=sharex, gridspec_kw=gridspec_kw, figsize=figsize)


def init_fig_horizontal(n=2, ratios=(1, 3), figsize=None, sharey='col', **kwargs):
    import matplotlib.pyplot as plt
    gridspec_kw = {i: kwargs[i] for i in ['height_ratios', 'hspace', 'wspace'] if i in kwargs.keys()}
    gridspec_kw.update({'width_ratios': ratios})
    return plt.subplots(1, n, sharey=sharey, gridspec_kw=gridspec_kw, figsize=figsize)
