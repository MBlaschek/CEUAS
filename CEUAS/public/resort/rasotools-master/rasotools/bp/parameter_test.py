# -*- coding: utf-8 -*-
import itertools

import numpy as np
import pandas as pd

#__all__ = ['parameter_test']

#
# def parameter_test(data, window_range=None, dist_range=None, thres_range=None, level_range=None, mp_process=True):
#     """ Run a Parameter test on all cases given some Data
#
#     Parameters
#     ----------
#     Data                DataFrame       2D time x plevels xData
#     window_range        list
#     dist_range          list
#     thres_range         list
#     level_range         list
#
#     Returns
#     -------
#     DataFrame of Cases
#     """
#     global testdata  # is this necessary ?
#
#     if window_range is None:
#         window_range = np.arange(365, 4*365+1, 365/4.)  # 1 to 4 years by 1/4
#
#     if dist_range is None:
#         dist_range = np.arange(182.5, 2*365+1, 365/4.)  # 1/1 to 2 years by 1/4
#
#     if thres_range is None:
#         thres_range = np.arange(25, 100, 5)
#
#     if level_range is None:
#         level_range = np.arange(1, 11)
#
#     dist_range = np.int_(dist_range)  # requires int
#     window_range = np.int_(window_range)
#     thres_range = np.int_(thres_range)
#     # w * d * t * l = 10 * 10 * 10 * 10 = 10000 simulations/ combinations
#     cases = list(itertools.product(window_range, dist_range, thres_range, level_range))
#     num = len(cases)
#     print("Running %d combinations of ranges: " % num)
#     print("Window (%d): %f to %f" % (len(window_range), np.min(window_range), np.max(window_range)))
#     print("Dist   (%d): %f to %f" % (len(dist_range), np.min(dist_range), np.max(dist_range)))
#     print("Thres  (%d): %f to %f" % (len(thres_range), np.min(thres_range), np.max(thres_range)))
#     print("Level  (%d): %f to %f" % (len(level_range), np.min(level_range), np.max(level_range)))
#
#     arg = input("Shall we run these marvelous cases? (y/n)")
#     if arg not in ['y', 'Y', 'yes']:
#         print("What a pity.")
#         return
#
#     print("You are the boss")
#     testdata = data.values.copy()
#     res = []
#     for icase in cases:
#         res.append(detection_wrapper(icase))
#
#     allcases = []
#     allbreaks = {}
#     i = 0
#     for iparams, ibreaks in res:
#         allcases.append(iparams)
#         allbreaks[i] = pd.Series(1, index=data.dims.date.values[ibreaks], name=i)
#         i += 1
#
#     allcases = pd.DataFrame(allcases)
#     allbreaks = pd.DataFrame(allbreaks)
#     allbreaks.index = allbreaks.index.to_datetime()  # fix index
#
#     return allcases, allbreaks

#
# def detection_wrapper(arg):
#     from .det import detector_2d
#     global testdata
#     window, dist, thres, min_levels = arg
#     params = {'window': window, 'dist': dist, 'thres': thres, 'level': min_levels, 'nbreaks': 0, 'vert': -1}
#     ibreaks = []
#     status, stest, breaks = detector_2d(testdata, 0, window, dist, thres, min_levels)
#     if status:
#         ibreaks = np.unique(np.where(breaks > 0)[0])
#         itl = np.median(np.where(breaks > 1)[1])  # above threshold / mean of levels -> vertical position
#         params.update({'nbreaks': len(ibreaks), 'vert': itl})
#     return params, ibreaks

    # stest = np.squeeze(np.apply_along_axis(test, 0, testdata.values, window, window / 4))  # per level
    # imax = np.asarray(local_maxima(np.sum(stest, 1), dist=dist))  # for all levels
    # if len(imax) == 0:
    #     return params, []  # no local maxima
    #
    # tstat = np.sum(stest[imax, :] > thres, 1)
    # if not np.any(tstat >= min_levels):
    #     return params, []  # not enough above threshold
    #
    # itx = imax[tstat >= min_levels]  # indices of breaks
    # ibreaks = testdata.index[itx].tolist()
    # itl = np.median(np.where(stest[itx] > thres)[1])  # above threshold / mean of levels -> vertical position
    # params.update({'nbreaks': len(itx), 'vert': itl})
    # return params, ibreaks
