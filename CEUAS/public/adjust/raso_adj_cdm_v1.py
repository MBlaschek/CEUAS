#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# This script has been developed in the service contract for C3S
# Calculates radiosonde humidity adjustments based on CDM files
#
# (c) University of Vienna, M. Blaschek, Vienna, Austria
# Copernicus Climate Change Service, 2020
# https://apps.ecmwf.int/datasets/licences/copernicus/
# email michael.blaschek (at) univie.ac.at
# Created: Vienna, 26 August, 2019
# Last Modifed:  9 March, 2021
# Version: 0.2
# -----------------------------------------------------------------------------
__version__ = '0.2'
__author__ = 'MB'
__status__ = 'dev'
__date__ = 'Di 09 Mär 2021 16:18:08 CET'
__institute__ = 'UNIVIE'
__github__ = 'git@github.com:MBlaschek/CEUAS.git'
__doc__ = """
CDS_adjust Functions v%s
Maintained by %s at %s
Github: %s [%s]
Updated: %s
""" % (__version__, __author__, __institute__, __github__, __status__, __date__)

import logging
import os
import sys
import warnings

import numpy as np
from numpy.core.fromnumeric import var
import pandas as pd
import xarray as xr
from xarray.core import dataarray
from numba import njit
import ray
import copy
import time
import h5py
os.chdir(os.path.expandvars('$HOME/CEUAS/CEUAS/public/adjust'))
sys.path.insert(0,os.getcwd()+'/../resort/')
from convert_numbamod import Sonntag, vap2sh, frostpoint_Sonntag, dewpoint_Sonntag
sys.path.insert(0,os.getcwd()+'/../merge/')
from harvest_convert_to_netCDF_yearSplit import write_dict_h5
import pyRAOBCORE_numbamod as pnm
from datetime import datetime, timedelta


sys.path.insert(0,os.getcwd()+'/../resort/rasotools-master/')
import rasotools

try:
    sys.path.append('../cds-backend/code/')
    import cds_eua4 as eua
except Exception as e:
    print('CDS_EUA4 Module is required. Add to path')
    raise e

warnings.simplefilter("ignore")
np.seterr(invalid='ignore')

logger = logging.getLogger()
logger.name = 'adjust'
formatter = logging.Formatter('%(asctime)s - %(name)s | %(funcName)s - %(levelname)s - %(message)s')

# in Pa
std_plevs = np.asarray([10, 20, 30, 50, 70, 100, 150, 200, 250, 300, 400, 500, 700, 850, 925, 1000])


# -----------------------------------------------------------------------------
#
# Helper functions
#
# -----------------------------------------------------------------------------
@njit
def rms(x):
    #x = np.array(x)
    mask = ~np.isnan(x)
    if np.any(mask):
        return np.sqrt(np.mean(x[mask]*x[mask]))
    else:
        return np.nan
    
@njit
def rrms(x):
    #x = np.array(x)
    mask = ~np.isnan(x)
    if np.any(mask):
        q75 = np.quantile(x[mask], 0.75)
        return q75
    else:
        return np.nan

def printstats(data, variable,bc, func, tstart=None, tstop=None, extended=False):
    
    sl = slice(0, data.time.shape[0])
    if type(tstart) ==str and type(tstop) == str:
        try:
            dstart = np.datetime64(tstart)
            dstop = np.datetime64(tstop)
            sl = slice(*np.searchsorted(data.time.values, (dstart, dstop)))
        except:
            pass
    
    print(f'{variable}, {bc}, {tstart}-{tstop}')
    
    for ip in range(16):
        if extended or ip == 9:  # 300 hPa
            odat = func(data[variable][sl, ip].values)
            bcdat = func(data[variable][sl, ip].values-data[bc][sl, ip].values)
            if odat > 0:           
                print(f'{int(data.plev[ip]): >6d}: {odat:.4f}, {bcdat:.4f}, {bcdat / odat:.4f}')
            else:
                print(f'{int(data.plev[ip]): >6d}: {odat:.4f}, {bcdat:.4f}')
            
    return

@njit
def rmean(t,tmean,index,runmean):

    tret=np.zeros(t.shape[0])
    ni=t.shape[0]
    good=runmean-runmean
    if runmean<2:
        for i in range(ni):
            tret[i]=t[i]
    else:

        for j in range(ni):
            tret[j]=np.nan
            if t[j]==t[j]:
                index[good]=j
                good+=1

        if good>runmean+2:
            i=runmean//2
            tmean[:]=np.nan
            if runmean%2==1:
                tmean[i]=0.
                for k in range(-runmean//2+1,runmean//2+1):
                    tmean[i]+=t[index[i+k]]
                tmean[i]/=runmean

                for i in range(runmean//2+1,good-runmean//2):
                    tmean[i]=(tmean[i-1]*runmean+t[index[i+runmean//2]])/(runmean+1)

            else:

                i=runmean//2
                tmean[i]=0.
                for k in range(-runmean//2,runmean//2):
                    tmean[i]+=t[index[i+k]]
                tmean[i]/=runmean

                for i in range(runmean//2+1,good-runmean//2-1):
                    tmean[i]=(tmean[i-1]*runmean+t[index[i+runmean//2-1]])/(runmean+1)

            for i in range(good):
                tret[index[i]]=tmean[i]
        else:
            for i in range(good):
                tret[index[i]]=t[index[i]]

    return tret

def rmeanw(t,runmean):
    tmean=t.copy()
    index=np.zeros(tmean.shape[0],dtype='int')
    
    tret=rmean(t,tmean,index,runmean)
    tret[:runmean]=np.nan
    tret[-runmean:]=np.nan
    return tret


def update_kw(name, value, **kwargs):
    """ Update keyword dictionary on the fly
    """
    kwargs.update({name: value})
    return kwargs


def conform(data, shape):
    """ Make numpy array conform to a certain shape
    Args:
        data (np.ndarray): input data
        shape (tuple, list): desired shape

    Returns:
        np.ndarray : reshaped
    """
    if not isinstance(data, np.ndarray):
        raise ValueError('Requires a numpy array')

    if not isinstance(shape, (tuple, list)):
        raise ValueError('Requires a tuple or list')

    data = data.copy()
    n = data.shape

    assert np.any([i in shape for i in n]), "Shapes do not allign?!"

    for i, j in enumerate(shape):
        if j not in n:
            data = np.expand_dims(data, axis=i)
    return data


def nancount(x, axis=0, keepdims=False):
    """ Count values excluding NaN
    Args:
        x (np.ndarray): input dataset
        axis (int): axis
        keepdims (bool): keep dimensions
    Returns:
        np.ndarray : sum of values
    """
    return np.sum(np.isfinite(x), axis=axis, keepdims=keepdims)


def nanfunc(data, n=130, axis=0, nmax=1460, borders=0, ffunc=None, flip=False, fargs=(), **kwargs):
    """ Nan omitting function (numpy)
    Args:
        data (np.ndarray): dataset including NaN
        n (int): minimum sample size
        axis (int): datetime axis
        nmax (int): maximum sample size
        borders (int): border sample to ignore
        ffunc (callable): function to call
        flip (bool): reverse dataset before applying the function
        fargs (tuple): function arguments

    Returns:
        np.ndarray : func of values at axis, with sample size, borders and maximum
    """
    if ffunc is None:
        ffunc = np.nanmean
    return np.apply_along_axis(sample, axis, data, n, nmax, ffunc, borders=borders, flip=flip, fargs=fargs)


def sample(values, nmin, nmax, func, borders=0, flip=False, fargs=(), **kwargs):
    """ Apply a function (func) to a sample of defined size

    Args:
        values (np.ndarray): input values
        nmin (int): minimum required values
        nmax (int): maximum number of values
        func (callable): function to execute
        borders (int): number of values to skip (start-borders , end-borders)
        flip (bool): reverse order of array
        fargs (tuple): arguments to function func
        **kwargs:

    Returns:
        np.ndarray : according to func
    """
    itx = np.isfinite(values)
    n = itx.sum()
    j = 0
    if n > nmax:
        if n > (nmax + borders):
            j = borders
        if flip:
            return func(np.flip(values[itx])[j:(nmax + j)], *fargs)  # reversed
        return func(values[itx][j:(nmax + j)], *fargs)  # normal

    elif n < nmin:
        # raises all nan warnings !!!
        return func(values, *fargs) * np.nan

    else:
        if n > (nmin * 2 + borders):
            j = borders
        if flip:
            return func(np.flip(values[j:]), *fargs)
        return func(values[j:], *fargs)


def table_to_dataset(data, dim='time', plev='plev', levels=None, **kwargs):
    """ Convert pandas Dataframe to xarray Dataset

    Args:
        data (pd.DataFrame): input dataframe (columns are variables)
        dim (str): datetime dimension
        plev (str): pressure dimension
        levels (list): pressure levels to consider
        **kwargs:

    Returns:
        xr.Dataset : 2d (datetime x pressure levels) x variables
    """
    from xarray import Dataset
    if levels is None:
        levels = [1000., 2000., 3000., 5000., 7000., 10000., 15000., 20000., 25000., 30000., 40000., 50000., 70000.,
                  85000., 92500., 100000.]
    # dimensions for output
    varis = [dim, plev]
    attrs = None
    if isinstance(data, Dataset):
        # copy attributes
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
    logger.info("Selecting only standard pressure levels")
    data = data[data[plev].isin(levels)]
    #
    # convert to xarray
    #
    data = data.reset_index().set_index(varis)
    if not data.index.is_unique:
        data = data.loc[~data.index.duplicated()]  # remove duplicated

    data = data.to_xarray()  # 1D -> 2D
    if attrs is not None:
        # add attributes again
        for i, j in vatt.items():
            data[i].attrs.update(j)
        data.attrs.update(attrs)
        data[dim].attrs.update(tatts)

    return data


@np.vectorize
def fix_datetime(itime, span=6, debug=False):
    """ Fix datetime to standard datetime with hour precision

    Args:
        itime (datetime): Datetime
        span (int): allowed difference to standard datetime (0,6,12,18)

    Returns:
        datetime : standard datetime
    """
    import pandas as pd
    itime = pd.Timestamp(itime)  # (time: 34%)
    # span=6 -> 0, 12
    # [18, 6[ , [6, 18[
    # span=3 -> 0, 6, 12, 18
    # [21, 3[, [3,9[, [9,15[, [15,21[
    for ihour in range(0, 24, span * 2):
        # 0 - 6 + 24 = 18
        lower = (ihour - span + 24) % 24
        # 0 + 6 + 24 = 6
        upper = (ihour + span + 24) % 24
        # 18 >= 18 or 18 < 6  > 00
        # 0 >= 18 or 0 < 6    > 00
        if debug:
            logger.debug("%d [%d] %d >= %d < %d" % (ihour, span, lower, itime.hour, upper))

        if (ihour - span) < 0:
            if itime.hour >= lower or itime.hour < upper:
                rx = itime.replace(hour=ihour, minute=0, second=0, microsecond=0)
                if itime.hour >= (24 - span):
                    rx = rx + pd.DateOffset(days=1)
                return rx.to_datetime64()
        else:
            if lower <= itime.hour < upper:
                rx = itime.replace(hour=ihour, minute=0, second=0, microsecond=0)
                if itime.hour >= (24 - span):
                    rx = rx + pd.DateOffset(days=1)
                return rx.to_datetime64()


def add_adj(fi,adjustments, variable, varnum, mode='r'):
    
    #adjustments={'humidity_bias':fi,}
    adjname={'humidity_bias':'humidity_bias',}
    #outfile = fi
    if False: #xyz is None:
        statid=fi['sid'][-6:] #adjfile[-8:-3]
        adjfile = RC['outdir'] + statid + '/feedbackglobbincorrsave' + statid + '.nc'
    
        outfile = RC['outdir']+ statid+'/' + fi['ifileorig'].split('/')[-1]
        if mode == 'r+':
            try:
                with eua.CDMDataset(fi['ifileorig'], mode=mode) as data:
                    outfile = fi['ifileorig']
            except:
                mode = 'r'
                print('could not open read-write, opening readonly', fi['ifileorig'])
                
        print('writing to:', outfile)
    
        try:
            with eua.CDMDataset(fi['ifileorig'], mode=mode) as data:
                
            
                xyz = data.read_observed_variable(eua.cdm_codes['temperature'],return_xarray=True,date_time_in_seconds=True)
            
                ref=np.datetime64(datetime.datetime(1900,1,1),'ns')
                xyzt=(xyz.date_time.values-ref).astype('long')//1000000000
        except:
    
            pass
        
    try:
            
#            for k,v in adjustments.items():
                try:
                    xyz = fi.read_observed_variable(varnum,return_xarray=True,date_time_in_seconds=True)
                    ref=np.datetime64(datetime(1900,1,1),'ns')
                    xyzt=(xyz.date_time.values-ref).astype('long')//1000000000
                    
                    #atime0 = xyzt.copy()
                    #adjustments=xr.open_dataset(v,decode_times=False)
                    if adjustments.time.ndim==2:
                        atime0=(adjustments.datum[0].values.astype(int)-1)*86400.
                    else:
                        #atime0=(adjustments.datum.values.astype(int)-1)*86400.
                        #ref=np.datetime64(datetime(1900,1,1),'ns')
                        atime0=(adjustments.time.values-ref).astype('long')//1000000000
                   
                    
                    #mask=xyz[adjname[k]].values==-999.
                    #adjustments[adjname[k]].values[mask]=np.nan
                    tt=time.time()
                    sh = adjustments.values.shape
                    av = np.reshape(adjustments.values.T, (1,sh[1],sh[0] )) 
                    adj=pnm.add_biasestimate(xyz.values,xyzt,xyz.z_coordinate.values,atime0,
                                         av,adjustments.plev.values)
                    print('add:',time.time()-tt)
                    xyz.values=adj
                    
                    #idx=np.where(xyz.z_coordinate.values==50000)
                    #plt.plot(xyzt[idx]/86400/365.25,adj[idx])
                    #plt.plot(atime0/86400/365.25,adjustments.rasocorr.values[0,11,:])
                    #plt.plot(atime0/86400/365.25,adjustments.rasocorr.values[1,11,:])
                    #plt.show()
                    # Daten schreiben neue Variable monkey in neuer gruppe adjust
                    if mode == 'r':
                        if not os.path.isfile(outfile):
                            #shutil.copyfile(fi['ifileorig'], outfile)
                            do_copy(fi['ifileorig'], outfile)
                        with eua.CDMDataset(outfile, mode='r+') as odata:
                            odata.write_observed_data(k.upper()+'_bias_estimate',
                                                     ragged=xyz,  # input data
                                                     varnum=eua.cdm_codes['temperature'],  # observed_variable to be aligned with
                                                     group='advanced_homogenisation',   # name of the new group
                                                     data_time='date_time',  # named datetime coordinate
                                                     data_plevs='z_coordinate',  # named pressure coordinate
                                                     attributes={'version':RC['version']}
                                                    )
                            print(k, time.time() - tt)
                    else:
                        #key = variable
                        #alldict = pd.DataFrame({key:xyz})
                        #del fi.file['advanced_homogenisation']['humidity_bias_estimate2']
                        fi.write_observed_data(variable,
                                                 ragged=xyz,  # input data
                                                 varnum=varnum,  # observed_variable to be aligned with
                                                 group='advanced_homogenisation',   # name of the new group
                                                 data_time='date_time',  # named datetime coordinate
                                                 data_plevs='z_coordinate',  # named pressure coordinate
                                                 attributes={'version':'1.0'}
                                                )
                    print('write:',time.time()-tt)
                except FileNotFoundError as e:
                    print('no adjustments found', k,e)
    except MemoryError as e:
        print('could not write '+outfile, e)

# -----------------------------------------------------------------------------
#
# Detection of Breakpoints
#
# -----------------------------------------------------------------------------
def detector(data, axis=0, dist=365, thres=50, min_levels=3, use_slopes=False, use_first=False, **kwargs):
    """ Detect breakpoints given some parameters

    Args:
        data (np.ndarray): snht test data
        axis (int): datetime dimension
        dist (int): minimum distance between breakpoints
        thres (int): threshold for significant levels
        min_levels (int): minimum of significant levels
        use_slopes (bool): Use first derivative for detection
        use_first (bool): Use the beginning of a peak, rather than the peak
        **kwargs:

    Returns:
        np.ndarray : breakpoint data
    """
    if not isinstance(data, np.ndarray):
        raise ValueError("requires an array: %s" % str(type(data)))

    breaks = (data >= thres).astype(int)  # 0, 1

    if data.ndim > 1:
        # number of breaks per timeunit
        jbreak = np.sum(breaks, axis=1 if axis == 0 else 0)
        # combine #-significant and sum of test
        jbreak = np.where(jbreak >= min_levels, 1, 0)
        ibreak = jbreak * np.sum(data, axis=1 if axis == 0 else 0)

    else:
        ibreak = data * breaks  # truncate at threshold

    # find local maxima (within distance)
    if use_slopes:
        imax = local_maxima(np.diff(ibreak), dist=dist)

    elif use_first:
        # find the beginning of breakpoint
        imax = local_maxima(ibreak, dist=dist)  # Maximum
        ifirst = []
        for i in imax:
            ifirst += [i + np.argmin(ibreak[i:i + dist])]
        imax = ifirst

    else:
        imax = local_maxima(ibreak, dist=dist)

    if len(imax) > 0:
        imax = np.asarray(imax)
        logger.info("Breaks: " + str(imax))
        for i in imax:
            breaks[idx2shp(i, axis, data.shape)] += 2  # maximum Breakpoint

    return breaks


def idx2shp(idx, axis, shape):
    """ Generate indices

    Args:
        idx (int): index
        axis (int): axis
        shape (tuple): shape

    Returns:
        tuple : selection slice
    """
    index = [slice(None)] * len(shape)
    index[axis] = idx
    return tuple(index)


@njit
def local_maxima(x, dist=365):
    """ Find local maxima, using Numba

    Args:
        x (np.ndarray): input
        dist (int): minimum distance between peaks

    Returns:
        list : list of maxima
    """
    maxima = []  # Leere Liste
    # Iteriere von 2 bis vorletzten Element
    # ist iterator integer
    for i in range(dist, len(x) - dist):
        # Element davor und danach größer
        if np.all((x[i] > x[slice(i + 1, i + dist)])) and np.all((x[slice(i - dist, i)] < x[i])):
            maxima.append(i)
    return maxima


def test(x, window, missing):
    """Standard Normal Homogeneity Test (SNHT) with a running window
    Wrapper function for numba_snhtmov

    Args:
        x (np.ndarray) : input data
        window (int) : window size (in days)
        missing (int) : allowed missing values (in days)
    Returns:
        np.ndarray : SNHT
    """

    snhtparas = np.asarray([window, missing, 10])
    tsa = np.zeros(x.shape[0])
    tmean = np.zeros(x.shape[0])
    tsquare = np.zeros(x.shape[0])
    count = np.zeros(x.shape[0], dtype=np.int32)

    if len(x) > 1:
        
        numba_snhtmov(np.squeeze(np.asarray(x)),
                      tsa,
                      snhtparas,
                      count,
                      tmean,
                      tsquare)
    else:
        pass
    
    return tsa


@njit(cache='False')
def numba_snhtmov(t, tsa, snhtparas, count, tmean, tsquare):
    """Standard Normal Homogeneity Test Moving Window

    t         = np.random.randn(1000)
    snhtparas = np.asarray([100,50,10])
    tsa       = np.zeros(1000)
    tmean     = np.zeros(1000)
    tsquare   = np.zeros(1000)
    count     = np.zeros(1000,dtype=np.int32)

    Output: tsa
    """
    n = snhtparas[0]
    max_miss = snhtparas[1]
    # ninc=snhtparas[2]

    ni = t.shape[0]
    good = 0
    tmean[0] = 0.
    tsquare[0] = 0.
    for j in range(ni):
        count[j] = 0
        # compare_lists if nan ?
        if t[j] == t[j]:
            if good > 0:
                tmean[good] = tmean[good - 1] + t[j]
                tsquare[good] = tsquare[good - 1] + t[j] * t[j]
            else:
                tmean[good] = t[j]
                tsquare[good] = t[j] * t[j]
            good += 1
        if good > 0:
            count[j] = good - 1

    if good > n - 2 * max_miss:
        rm = int(n / 2)  # needs to be an integer
        # k 1460/2=730 - 650=80, n-80
        for k in range(rm - max_miss, ni - (rm - max_miss)):
            xm = k - rm  # 80-730
            if xm < 0:
                xm = 0
            xp = k + rm
            if xp > ni - 1:
                xp = ni - 1
            if (count[k] - count[xm] > rm - max_miss) and (count[xp] - count[k] > rm - max_miss):
                x = (tmean[count[k]] - tmean[count[xm]]) / (count[k] - count[xm])  # Mittelwert 1 Periode
                y = (tmean[count[xp]] - tmean[count[k]]) / (count[xp] - count[k])  # Mittelwert 2 Periode
                xy = (tmean[count[xp]] - tmean[count[xm]]) / (count[xp] - count[xm])  # Mittelwert ganze Periode

                sig = (tsquare[count[xp]] - tsquare[count[xm]]) / (count[xp] - count[xm])  # t*t ganze Periode
                if sig > xy * xy:
                    sig = np.sqrt(sig - xy * xy)  # standard deviation of the whole window
                    # n1 * (m1-m)**2 + n2 * (m2-m)**2 / stddev
                    tsa[k] = ((count[k] - count[xm]) * (x - xy) * (x - xy) + (count[xp] - count[k]) * (y - xy) * (
                            y - xy)) / (sig * sig)
                else:
                    tsa[k] = 0.
    return


def get_breakpoints(data, value=2, dim='time', return_startstop=False, startstop_min=0, **kwargs):
    """ Return breakpoints from breakpoint data

    Args:

        data (xr.DataArray): input DataArray
        value (int): breakpoint indicator value
        dim (str): datetime dim
        startstop_min (int):
        return_startstop (bool):
        **kwargs:

    Returns:
        list : breakpoints
    """
    if not isinstance(data, xr.DataArray):
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

    if len(i) > 0:
        logger.debug("Breakpoints for %s" % data.name)
        logger.debug("[%8s] [%8s] [%8s] [%8s] [ #]" % ('idx', 'end', 'peak', 'start'))
        for j, k, l in zip(i, s, e):
            logger.debug("[%8s] %s %s %s %4d" % (j, dates[l], dates[j], dates[k], k - l))

    if return_startstop:
        return i, e, s
    return i


# -----------------------------------------------------------------------------
#
# Adjustments
#
# -----------------------------------------------------------------------------

def adjustments(data, breaks, use_mean=False, axis=0, sample_size=130, borders=30, max_sample=1460, recent=False,
                ratio=False,meta=None, **kwargs):
    """ Adjustment of breakpoints

    Args:
        data (np.ndarray): radiosonde data
        breaks (list): breakpoint indices
        use_mean (bool): mean or quantile adjustments
        axis (int): axis of datetime
        sample_size (int): minimum sample size
        borders (int): adjust breaks with borders
        max_sample (int): maximum sample size
        recent (bool): use full reference period
        ratio (bool): calculate ratio, instead of difference
        **kwargs:

    Returns:
        np.ndarray : adjusted data
    """
    if not isinstance(data, np.ndarray):
        raise ValueError("requires a numpy array")
    if not isinstance(breaks, (np.ndarray, list)):
        raise ValueError('requires a numpy array')

    data = data.copy()
    if not use_mean:
        percentilen = kwargs.get('percentiles', [0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100])
        nq = len(percentilen)
        #
        # Sample size per Percentile
        #
        sample_size = sample_size // nq
        if sample_size < 3:
            sample_size = 3
        logger.info('Sample size: %d N-Q: %d', sample_size, nq)

    dshape = data.shape  # Shape of data (date x levs)
    imax = dshape[axis]  # maximum index
    breaks = np.sort(np.asarray(breaks))  # sort
    breaks = np.append(np.insert(breaks, 0, 0), imax)  # 0 ... ibreaks ... Max
    breaks = breaks.astype(int)
    nb = breaks.size

    for i in range(nb - 2, 0, -1):
        # Indices
        im = breaks[i - 1]  # earlier
        ib = breaks[i]  # current breakpoint
        if recent:
            ip = imax  # Max
        else:
            ip = breaks[i + 1]  # later

        # Slices all axes
        iref = slice(ib, ip)
        isample = slice(im, ib)
        isample = idx2shp(isample, axis, dshape)
        iref = idx2shp(iref, axis, dshape)
        # Before Adjustments
        before = np.nanmean(data[isample], axis=axis)
        vais={'ref': {'tup':np.unique(meta[0,slice(ib,ip)],return_counts=True,return_inverse=True),'vcount':0,'nvcount':0,'vi':[],'nvi':[]},
            'sample': {'tup':np.unique(meta[0,slice(im,ib)],return_counts=True,return_inverse=True),'vcount':0,'nvcount':0,'vi':[],'nvi':[]},}
        sdict=dict(
            rs90 = np.array(['71', '72', '73', '74','78', ],dtype='S4'),
            rs92 = np.array(['70', '79', '80', '81'  ],dtype='S4'),
            rs80 = np.array(['37', '52', '60', '61', '62', '63' , '66', '67','VC_'  ],dtype='S4'),
            rs41 = np.array(['123', '124', '125', '141', '142'],dtype='S4'),
            graw = np.array(['17'],dtype='S4'),
            viz  = np.array(['49','ZKB'],dtype='S4'),
            vschroeder = np.array(['VN', 'VP'],dtype='S4'),  
            )
        debug=False
        if debug:
            import matplotlib.pyplot as plt
        # Apply Adjustments
        if use_mean:
            data[isample] = mean(data[iref], data[isample],
                                 axis=axis,
                                 sample_size=sample_size,
                                 max_sample=max_sample,
                                 borders=borders,
                                 ratio=ratio, **kwargs)
        else:
            for k in vais.keys():
                for it in range(vais[k]['tup'][0].shape[0]):
                    t=vais[k]['tup'][0][it]
                    if t in sdict['rs92'] or t in sdict['rs90'] or t in sdict['rs41'] or b'VN' in t or b'VP' in t:
                        vais[k]['vcount']+=vais[k]['tup'][2][it]
                    elif t!=b'nan': # and t not in sdict['rs80']:
                        vais[k]['nvcount']+=vais[k]['tup'][2][it]
                        vais[k]['nvi'].append(np.where(vais[k]['tup'][1]==it)[0])
                    elif t==b'nan':
                        vais[k]['nvi'].append(np.where(vais[k]['tup'][1]==it)[0])
                        
            if debug:
                plt.subplot(1,2,1)
                plt.hist(data[isample][:,7],label='sample')
                plt.hist(data[iref][:,7],label='ref')
                plt.legend()        
            print('min0:',np.nanmin(data[iref]),np.nanmin(data[isample]))
            if vais['ref']['vcount']>=vais['ref']['nvcount'] or vais['sample']['vcount']<vais['sample']['nvcount']:
                    
                vais['sample']['nvi']=np.concatenate(vais['sample']['nvi'])
                #if vais['sample']['nvi'].shape[0]>data[isample].shape[0]:
                #   vais['sample']['nvi']=np.arange(data[isample].shape[0],dtype=int) 
                d=data[isample]
                d[vais['sample']['nvi']] = percentile(data[iref], data[isample][vais['sample']['nvi']], percentilen,
                                       axis=axis,
                                       sample_size=sample_size,
                                       max_sample=max_sample,
                                       borders=borders,
                                       ratio=ratio, **kwargs)
                d[d<np.nanmin(data[iref])]=np.nanmin(data[iref])
                data[isample]=d[:]
            else:
                vais['ref']['nvi']=np.concatenate(vais['ref']['nvi'])
                d=data[iref]
                d[vais['ref']['nvi']] = percentile(data[isample], data[iref][vais['ref']['nvi']], percentilen,
                                       axis=axis,
                                       sample_size=sample_size,
                                       max_sample=max_sample,
                                       borders=borders,
                                       ratio=ratio, **kwargs)
                d[d<np.nanmin(data[isample])]=np.nanmin(data[isample])
                data[iref]=d[:]
            if debug:
                plt.subplot(1,2,2)
                plt.hist(data[isample][:,7],label='sample after')
                plt.hist(data[iref][:,7],label='ref after')
                plt.legend()
        print('min:',np.nanmin(data[iref]),np.nanmin(data[isample]))
        if debug:
            plt.show()

        #
        # Border zone (Break, Break + borders)
        # Linear interpolation
        #
        if use_mean and borders > 0:
            zone = slice(ib, ib + borders)
            zone = idx2shp(zone, axis, dshape)
            linzone = conform(np.linspace(0, 1, data[zone].shape[axis]), data[zone].shape)
            if ratio:
                zsample = np.nanmean(data[isample], axis=axis) / np.nanmean(data[iref], axis=axis)
                data[zone] = data[zone] * linzone * zsample
            else:
                zsample = np.nanmean(data[isample], axis=axis) - np.nanmean(data[iref], axis=axis)
                data[zone] = data[zone] + linzone * zsample

    return data


def mean(sample1, sample2, axis=0, sample_size=130, borders=0, max_sample=1460, ratio=True,
         median=False, **kwargs):
    """ Adjustment method using mean differences or ratios

    ratio=False
    dataset[sampleout]  + (MEAN(dataset[sample1]) - MEAN(dataset[sample2]))

    ratio=True
    dataset[sampleout]  * (MEAN(dataset[sample1]) / MEAN(dataset[sample2]))

    Args:
        sample1 (np.ndarray): reference
        sample2 (np.ndarray): sample
        axis (int): date axis
        sample_size (int): minimum sample size
        ratio (bool): use ratio or difference?
        median (bool): use median instead of mean?
        borders (int): around breakpoint
        max_sample (int): maximum sample size

    Returns:
        np.ndarray : mean adjusted dataset
    """
    # minimum sample size, maximum sample size
    if median:
        s1 = nanfunc(sample1,
                     axis=axis,
                     n=sample_size,
                     nmax=max_sample,
                     ffunc=np.nanmedian,
                     borders=borders)
        s2 = nanfunc(sample2,
                     axis=axis,
                     n=sample_size,
                     nmax=max_sample,
                     ffunc=np.nanmedian,
                     borders=borders,
                     flip=True)
    else:
        s1 = nanfunc(sample1,
                     axis=axis,
                     n=sample_size,
                     nmax=max_sample,
                     ffunc=np.nanmean,
                     borders=borders)
        s2 = nanfunc(sample2,
                     axis=axis,
                     n=sample_size,
                     nmax=max_sample,
                     ffunc=np.nanmean,
                     borders=borders,
                     flip=True)

    if ratio:
        # Todo factor amplifies extreme values
        dep = s1 / s2
        dep = np.where(np.isfinite(dep), dep, 1.)  # replace NaN with 1
        sample2 *= dep
    else:
        dep = s1 - s2
        sample2 += dep
    return sample2


def percentile(sample1, sample2, percentiles, axis=0, sample_size=130, borders=0, max_sample=1460, ratio=True,
               apply=None, noise=False, **kwargs):
    """ Adjustment method using percentile differences or ratios

    ratio=False
    dataset[sample1] + ( percentiles(dataset[sample1]) - percentiles(dataset[sample2]) )

    ratio=True
    dataset[sample1] * ( percentiles(dataset[sample1]) / percentiles(dataset[sample2]) )

    Args:
        sample1 (np.ndarray): reference
        sample2 (np.ndarray): sample
        percentiles (list): percentiles to use
        axis (int): date axis
        sample_size (int): minimum sample size
        ratio (bool): use ratio or difference?
        borders (int): around breakpoint
        max_sample (int): maximum sample size
        apply (np.ndarray): apply adjustments to this array
        noise (bool): add random noise to adjustments

    Returns:
        np.ndarray : percentile adjusted dataset
    """
    # Add 0 and 100, and remove them
    percentiles = np.unique(np.concatenate([[0], percentiles, [100]]))
    percentiles = percentiles[1:-1]  # remove 0 and 100
    #
    # Percentiles can be duplicated (because DPD might be integers)
    # limit calculations by sample_size, max_sample, borders
    #
    #     (part A)    |    (part B)
    #               break
    #             >borders<
    #  <max sample         max sample>
    #
    # (part B)
    s1 = nanfunc(sample1,
                 axis=axis,
                 n=sample_size,
                 nmax=max_sample,
                 ffunc=np.nanpercentile,
                 borders=borders,
                 fargs=(percentiles,))
    # flip means that beginning from the back (part A)
    s2 = nanfunc(sample2,
                 axis=axis,
                 n=sample_size,
                 nmax=max_sample,
                 ffunc=np.nanpercentile,
                 borders=borders,
                 fargs=(percentiles,),
                 flip=True)
    if ratio:
        dep = np.divide(s1, s2, where=(s2 != 0), out=np.full(s2.shape, 1.))
        dep = np.where(np.isfinite(dep), dep, 1.)  # replace NaN
    else:
        dep = s1 - s2
        dep = np.where(np.isfinite(dep), dep, 0.)

    # Interpolate adjustments
    if apply is None:
        apply = sample2.copy()
    else:
        apply = apply.copy()

    dep = apply_percentile_adjustments(apply, s2, dep, axis=axis, noise=noise)

    if ratio:
        dep = np.where(np.isfinite(dep), dep, 1.)
        apply *= dep
    else:
        dep = np.where(np.isfinite(dep), dep, 0.)
        apply += dep

    return apply


def apply_percentile_adjustments(data, percentiles, adjustment, axis=0, noise=False):
    """ Helper Function for applying percentile adjustments

    Args:
        data (np.ndarray): dataset
        percentiles (np.ndarray): percentiles, points of adjustments
        adjustment (np.ndarray): adjustments to be interpolated
        axis (int): axis of datetime
        noise (bool): add random noise to adjustments

    Returns:
        np.ndarray : interpolated adjustment, same shape as dataset
    """
    in_dims = list(range(data.ndim))
    # last dim == axis, Last dim should be time/date
    # print(data.shape)
    data = np.transpose(data, in_dims[:axis] + in_dims[axis + 1:] + [axis])
    # print(data.shape)
    percentiles = np.transpose(percentiles, in_dims[:axis] + in_dims[axis + 1:] + [axis])
    adjustment = np.transpose(adjustment, in_dims[:axis] + in_dims[axis + 1:] + [axis])
    # print(percentiles.shape, adjustment.shape)
    adjusts = np.zeros(data.shape)
    # Indices for iteration + expand
    inds = np.ndindex(data.shape[:-1])  # iterate all dimensions but last
    inds = (ind + (Ellipsis,) for ind in inds)  # add last as ':' == Ellipsis == all
    # k = 0
    for ind in inds:
        # INTERP -> Xnew, Xpoints, Fpoints
        iperc, idx = np.unique(percentiles[ind], return_index=True)
        iadj = adjustment[ind][idx]
        if noise:
            adjusts[ind] = np.interp(data[ind] + np.random.normal(size=data[ind].size, scale=0.5), iperc, iadj,
                                     left=np.nan, right=np.nan)
        else:
            adjusts[ind] = np.interp(data[ind], iperc, iadj, left=np.nan,
                                     right=np.nan)
        # if k == 5:
        #     print(ind, data[ind], adjusts[ind], iperc, iadj)
        # k+=1

    # Transform back to original shape
    return np.transpose(adjusts, in_dims[:axis] + in_dims[axis + 1:] + [axis])


# -----------------------------------------------------------------------------
#
# Main function that does all the steps in adjusting
#
# -----------------------------------------------------------------------------

def adjustment_procedure(data: xr.Dataset, dim: str = 'time', plev: str = 'plev', obs_name: str = 'obs', dep_name:str = 'dep',
                         metadata: bool = False, times: list = [0, 12], span: int = 3, freq: str = '12h',
                         mean_adjustments:bool = False, quantile_adjustments:bool=True,
                         **kwargs) -> xr.Dataset:
    """Main function for RASO_ADJ_CDM_v0
    This function executes the homogenization routines and deals with input and output

    Args:
        data (xr.Dataset): [description]
        dim (str, optional): [description]. Defaults to 'time'.
        plev (str, optional): [description]. Defaults to 'plev'.
        obs_name (str, optional): [description]. Defaults to 'obs'.
        dep_name (str, optional): [description]. Defaults to 'dep'.
        metadata (bool, optional): [description]. Defaults to False.
        times (list, optional): [description]. Defaults to [0, 12].
        span (int, optional): [description]. Defaults to 3.
        freq (str, optional): [description]. Defaults to '12h'.
        mean_adjustments (bool, optional): [description]. Defaults to True.
        quantile_adjustments (bool, optional): [description]. Defaults to False.

    Raises:
        ValueError: [description]
        ValueError: [description]
        ValueError: [description]
        RuntimeError: [description]

    Returns:
        xr.Dataset: [description]
    """
    #     thres (int): Threshold value for SNHT, default: 50
    #     window (int): Moving Window for SNHT, default: 1470 (in days, 4 years)
    #     missing (int): Maximum allowed missing values in window, default: 600 (in days)
    #     min_levels (int): Minimum required levels for significant breakpoint, default: 3
    #     dist (int): Minimum distance between breakpoints, default: 730 (in days, 2 years)
    #     sample_size (int): Minimum sample size for statistics, default: 130 (in days)
    #     borders (int): Breakpoint zone, default: 90 (in days)
    #     ratio (int): Use ratio instead of differences, default: 0 (not)

    # data [time x plev]  - > [hour x time x plev] at [times]
    if not isinstance(data, xr.Dataset):
        raise ValueError('Requires a xarray DataArray, not', type(data))
    if obs_name not in data.variables:
        raise ValueError('Required observation variable not present: ', obs_name)
    if dep_name not in data.variables:
        raise ValueError('Required departure variable not present: ', dep_name)
    bkptime = data[dim].copy()
    #
    # e.g. 00Z (21Z -1 day to 3Z same day)
    # standard_time
    #
    data[obs_name] = eua.align_datetime(data[obs_name], times=times, span=span, freq=freq, dim=dim, plev=plev)
    icoord = 'standard_%s' % dim
    standard_index = np.where(data[obs_name][icoord + '_flag'].values == 1)[0]
    # create a backup copy of the index/datetime/obs -> restore
    # overwrite dim with standard_dim
    data = data.assign_coords({dim: data[obs_name][icoord]})
    del data[icoord]
    del data[icoord + '_flag']
    # Check if result will be sorted
    #reverse_sort = False
    #if not data.isel(**{dim: standard_index})[dim].to_index().is_monotonic:
        #logger.warning('Datetime index is not monotonic %s', dim)
        ##idx = np.argsort(obs.isel(**{dim: standard_index})[dim].values)
        ##standard_index = standard_index[idx]
        #reverse_sort = True
    #
    # Convert to day-night Cube / drop duplicates by selection standard_index
    #
    # adds a stacking coordinate to unstack later
    data = eua.stack_cube_by_hour(data.isel(**{dim: standard_index}), dim=dim, times=times)
    
    # nur für test: tdata = eua.unstack_cube_by_hour(data, dim=dim)

    #
    # Some Variables for Breakpoint detection
    #
    window = kwargs.get('window', 1460)  # means 4 years on daily basis
    missing = kwargs.get('missing', 600)
    thres = kwargs.get('thres', 50)
    min_levels = kwargs.get('min_levels', 3)
    dist = kwargs.get('dist', 730)  # two years on daily basis
    #
    #
    #
    axis = data[dep_name].dims.index(dim)
    name = data[obs_name].name
    sname = data[obs_name].attrs.get('standard_name', name)
    #
    # 3. SNHT
    #
    
    if  len(data[dep_name].values) > 0:
        
        data['test'] = (data[dep_name].dims , np.apply_along_axis(test, axis, data[dep_name].values, window, missing))
    else:
        data['test'] = data[dep_name] + np.nan
    logger.info("Updated axis")
    #
    # Attributes
    #
    data['test'].name = '{}_snht'.format(name)
    logger.info("SNHT from %s [%s]", data[dep_name].name, data['test'].name)
    #
    # Day-Night Departures
    #
    data['test'] += np.apply_along_axis(test, axis - 1,
                                 data[obs_name].sel(hour=12).values - data[obs_name].sel(hour=0),
                                 window, missing)
    data['test'].attrs.update({'units': '1', 'window': window,
                        'missing': missing, 'standard_name': '{}_snht'.format(sname),
                        'day-night': 'added'})
    #
    # Check for Metadata in CDM
    #
    if metadata:
        print('using metadata',np.unique(data.meta.values[0].astype('S4')))
        #
        # In a future version Metadata can be read and used here.
        #
        pass
    else:
        data['meta']=None

    #
    # 4. Detect Breakpoints
    #
    breaks = xr.full_like(data[obs_name], 0, dtype=np.int32)
    breaks.name = '{}_breaks'.format(data['test'].name)
    breaks.attrs.update({'units': '1', 'dist': dist, 'thres': thres,
                         'min_levels': min_levels, 'standard_name': '{}_breaks'.format(sname)})
    attrs = data[obs_name].attrs.copy()
    attrs.update({'sample_size': kwargs.get('sample_size', 130),
                  'borders': kwargs.get('borders', 90),
                  'ratio': kwargs.get('ratio', 0)})
    for i, ihour in enumerate(times):
        breaks[i, ::] = detector(data['test'].values[i, ::], axis - 1, dist=dist, thres=thres, min_levels=min_levels)
        logger.info("Breakpoints detected %d/%d : %s", i, ihour, breaks.name)
    
    data['breakpoints'] = breaks
    #
    # 5. Adjust Breakpoints
    #
    if mean_adjustments:
        #
        # Temperature -> MEAN
        #
        data['adjustments'] = xr.full_like(data[dep_name], 0)
        data['adjustments'].name = '{}_m'.format(name)
        for i in range(len(times)):
            ibreaks = get_breakpoints(breaks[i, ::], **kwargs)
            logger.info("Breakpoints: %d", len(ibreaks))
            #
            adjv = adjustments(data[dep_name][i, ::].values, ibreaks,
                               axis=axis - 1,
                               **kwargs)
            # check limits [0 - 1]
            # obs + obs-an-adj - obs-an
            vadj = (data[obs_name][i, ::].values + adjv - data[dep_name][i, ::].values)
            adjv = np.where((vadj < 0) | (vadj > 1), data[dep_name][i, ::].values, adjv)
            # new = obs-an-adj + obs - (obs-an)
            data['adjustments'][i, ::] = (adjv - data[dep_name][i, ::].values)
        data['adjustments'].attrs.update(attrs)
        data['adjustments'].attrs['standard_name'] += '_adjustments'
        data['adjustments'].attrs['biascor'] = 'mean'
        

    elif quantile_adjustments:
        #
        # Relative Humidity -> QUANTILE
        #
        data['adjustments'] = xr.full_like(data[dep_name], 0)
        data['adjustments'].name = '{}_q'.format(name)
        for i in range(len(times)):
            ibreaks = get_breakpoints(breaks[i, ::], **kwargs)
            logger.info("Breakpoints: %d", len(ibreaks))
            #
            adjv = adjustments(data[dep_name][i, ::].values, ibreaks,
                                use_mean=False,
                                axis=axis - 1,
                                meta=data.meta.values.astype('S4'),
                                **kwargs)
            # check limits [0 - 1] of relative humidity
            # obs + obs-an-adj - obs-an
            vadj = (data[obs_name][i, ::].values + adjv - data[dep_name][i, ::].values)
            for j in range(len(data.plev)):
                nm = np.nanmin(data[obs_name][i, ::].values[:,j])
                idx = np.where((vadj[:,j] < nm))
                print(nm)               
                adjv[:,j] = np.where(vadj[:,j] < nm, data[dep_name][i, ::].values[:,j]+nm, adjv[:,j])
                adjv[:,j] = np.where(vadj[:,j] > 1, data[dep_name][i, ::].values[:,j], adjv[:,j])
            # new = obs-an-adj - (obs-an)
            data['adjustments'][i, ::] = (adjv - data[dep_name][i, ::].values)
        data['adjustments'].attrs.update(attrs)
        data['adjustments'].attrs['standard_name'] += '_adjustments'
        data['adjustments'].attrs['biascor'] = 'quantile'

        data['qadjustments'] = copy.deepcopy(data['adjustments'])
        data['qadjustments'].attrs['standard_name'] += 'specific_humidity_adjustments'
        data['dewpadjustments'] = copy.deepcopy(data['adjustments'])
        data['dewpadjustments'].attrs['standard_name'] += 'dewpoint_adjustments'
        
            
        obshom = data['relative_humidity'].values - data['adjustments'].values
        vpdata = data['relative_humidity'].values * Sonntag(data['temperature'].values )
        vpdatahom = obshom * Sonntag(data['temperature'].values )
        
        dp = dewpoint_Sonntag(vpdata)
        dphom = dewpoint_Sonntag(vpdatahom)
        #mask = vpdatahom < 10.
        #dp[mask] = frostpoint_Sonntag(vpdata[mask])
        #dphom[mask] = frostpoint_Sonntag(vpdatahom[mask])
        #if vpdatahom < 10.:
            #dp = frostpoint_Sonntag(vpdata)
            #dphom = frostpoint_Sonntag(vpdatahom)
        #else:
            #dp = dewpoint_Sonntag(vpdata)
            #dphom = dewpoint_Sonntag(vpdatahom)
        pl = data['plev'].values
        for i in range(data['plev'].shape[0]):
            adj = data['specific_humidity'].values[:, :, i] - vap2sh(vpdatahom[:, :, i], pl[i])
            data['qadjustments'].values[:, :, i] = -adj  # here the adjustment is needed, not the bias estimate
            dpadj = dp[:, :, i] - dphom[:, :, i]
            data['dewpadjustments'].values[:, :, i] = -dpadj  # here the adjustment is needed, not the bias estimate
        
        
        x = 0
    else:
        raise RuntimeError('Either mean_adjustment or quantile_adjustment needs to be set')
    # 
    # Convention (obs - bias_estimate = true_obs)
    #
    data['adjustments'] *= -1
    #
    # Convert back to time x plevs
    #
    try:
        
        data = eua.unstack_cube_by_hour(data, dim=dim)
    except ValueError:
        data = None
        return data
    # fill back (to the original input data)
    #if reverse_sort:
        ## idx = idx.sort()
        ## standard_index = standard_index[idx]  # undo sorting
        #logger.warning("Reverse sorting?")

    data = data.assign_coords({dim : bkptime[standard_index]})
    data = data.reindex({dim: bkptime})
    # Make sure we fill it up
    data['adjustments'] = data['adjustments'].fillna(0)
    data['qadjustments'] = data['qadjustments'].fillna(0)
    data['test'] = data['test'].fillna(0)
    data['breakpoints'] = data['breakpoints'].fillna(0).astype(int)
    return data


def run_frontend_file(args, **kwargs):
    if args.feedback is None:
        args.feedback = 'obs_minus_bg'  # Background: first guess, or 
                                        # obs_minus_an : Analysis departures

    if args.temperature or args.humidity:
        iofile = eua.CDMDataset(args.frontend)
        if args.outdir is not None:
            args.outdir = "{}/{}_adjusted.nc".format(args.outdir, os.path.basename(args.frontend).replace('.nc',''))
        else:
            args.outdir = "{}_adjusted.nc".format(args.frontend.replace('.nc',''))
        # Make sure we have this file open for writing
        iofile.reopen(write_to_filename=args.outdir, mode='r+', strict=True)
        #
        # Temperature or Humdity adjustment
        #
        variable = 'ta' if args.temperature else 'hur'
        # 
        data = iofile.read_data_to_cube([variable, args.feedback],
                                        dates=args.dates,
                                        plevs=args.plevs)
        # should contain variables
        # e.g. temperature, temperature_an_depar
        variable, depar = list(data.keys())
        data = xr.Dataset(data)
        # run adjustment procedure
        data = adjustment_procedure(data,
                                    obs_name=variable,
                                    dep_name=depar,
                                    metadata=False,
                                    times=[0, 12],
                                    dim='time',
                                    plev='plev',
                                    return_dataset=False,
                                    #quantile_adjustments=False if args.temperature else True
                                    )
        #
        # Write back adjusted (interpolation, extrapolation)
        #
        # create a copy of the file and write test, adjustments and breakpoint stats into it
        iofile.write_variable(variable + '_bias_estimate',
                              cube=data['adjustments'],
                              interpolate=args.interpolate_missing,
                              interpolate_datetime=args.interpolate_missing,
                              extrapolate_plevs=args.interpolate_missing
                              )
        iofile.write_variable(variable + '_snht_test',
                              cube=data['test'],
                              interpolate=args.interpolate_missing,
                              interpolate_datetime=args.interpolate_missing,
                              extrapolate_plevs=args.interpolate_missing
                              )
        iofile.write_variable(variable + '_breakpoints',
                              cube=data['breakpoints'],
                              interpolate=args.interpolate_missing,
                              interpolate_datetime=args.interpolate_missing,
                              extrapolate_plevs=args.interpolate_missing
                              )

    else:
        # Check if multiple files are given?
        if '*' in args.frontend:
            args.frontend = glob.glob(args.frontend)
        elif ',' in args.frontend:
            args.frontend = args.frontend.split(',')
        else:
            pass
        
        # convert to list and check files exist
        if not isinstance(args.frontend, list):
            args.frontend = [args.frontend]

        for i in args.frontend:
            assert os.path.isfile(i), i
            logger.info(i)
        #
        # Wind speed and wind direction adjustment
        #
        raise NotImplementedError()


def run_backend_file(args, **kwargs):
    # only one file
    if os.path.isfile(args.backend+'.xxx'):
        print(args.backend+' already processed')
        return

    try:
        hbi = False
        with h5py.File(args.backend) as f:
            try:
                if f['recordindices']['138'].shape[0] < 20:
                    print(args.backend+' humidity record too short')
                    return
            except:
                print(args.backend+' humidity record does not exist')
                return
            try:
                if f['advanced_homogenisation']['humidity_bias_estimate'].shape[0] > 0:
                    print(args.backend+' humidity bias exists')
                    hbi = True
                    #return
            except:
                pass
        
        if hbi:
            with h5py.File(args.backend, 'r+') as f:
                del f['advanced_homogenisation']['humidity_bias_estimate']

    except:
        raise ValueError(args.backend)
                
        
            
            
    iofile = eua.CDMDataset(args.backend)
    
    if not iofile.hasgroups:
        raise IOError("not a CDM backend file")
    
    if not args.feedback:
        args.feedback = 'fg_depar@offline'
        args.feedback_group = 'era5fb'
    
    if args.feedback_group not in iofile.groups:
        raise ValueError('Feedback Group', args.feedback_group, 'not found')
    
    if not args.feedback in iofile[args.feedback_group].keys():
        raise ValueError('Feedback ',args.feedback,' not in Feedback Group ',args.feedback_group)
    
    if args.outdir is not None:
        args.outdir = "{}/{}".format(args.outdir, os.path.basename(args.backend))

    # Make sure we have this file open for writing
    if args.temperature:
        raise NotImplementedError()
        # Code: 85
        variable = 'air_temperature'
        data = iofile.read_data_to_cube(variable,
                                        dates=args.dates,
                                        plevs=args.plevs,
                                        feedback=args.feedback,
                                        feedback_group=args.feedback_group,
                                        **kwargs)
        # should contain variables
        # e.g. temperature, temperature_an_depar
        variable, depar = list(data.keys())
        data = xr.Dataset(data)
        # run adjustment procedure
        data = adjustment_procedure(data,
                                    obs_name=variable,
                                    dep_name=depar,
                                    metadata=False,
                                    times=[0, 12],
                                    dim='time',
                                    plev='plev',
                                    return_dataset=False,
                                    )
        # TODO Convert adjustments to other variables?
        #
        #
        # Write back adjusted (interpolation, extrapolation)
        #
        iofile.write_observed_data('temperature_bias_estimate',
                                   varnum=eua.cdm_codes[variable],
                                   cube=data['adjustments'],
                                   group=args.homogenisation,
                                   interpolate=args.interpolate_missing,
                                   interpolate_datetime=args.interpolate_missing,
                                   extrapolate_plevs=args.interpolate_missing)

    if args.humidity:
        # Code 38, relative humidity, 
        variable = 'relative_humidity'
        # Code 34, dew point departure
        # variable = 'dew_point_departure'
        if not os.path.isfile(iofile.filename+'.xx'):
            
            if '138' not in iofile.recordindices.keys():
                return
        
            iofile.reopen(write_to_filename=args.outdir, mode='r+', strict=True)
            try:
                data = iofile.read_data_to_cube([variable,'specific_humidity', 'temperature'], #
                                                dates=args.dates,
                                                plevs=args.plevs,
                                                feedback=args.feedback,
                                                feedback_group=args.feedback_group,
                                                **kwargs)
                if args.dates is None:
                    trange = iofile.recordindices[f'{eua.cdm_codes[variable]}'][[0, -1]]
                    ri=iofile.recordindices[f'{eua.cdm_codes[variable]}'][:]
                    ref=np.datetime64(datetime(1900,1,1),'ns')
                    secs=(data['relative_humidity'].time.values-ref).astype('long')//1000000000
                    idx=np.searchsorted(iofile.file['observations_table']['date_time'][ri[0]:ri[-1]],secs)


                else:
                    trange = np.searchsorted(iofile.recordindices[f'{eua.cdm_codes[variable]}'][:], args.dates)[[0, -1]]
               

                if data['relative_humidity'].shape[0] <200:
                    print(iofile.filename, ' too short humidity record')
                    return

                try:
                    
                    data['meta'] = xr.DataArray(iofile.file['observations_table']['sensor_id'][ri[0]:ri[-1]][idx].view('S4').flatten(),
                                            coords={'time':data[variable].time}).astype('S4')
                except Exception as e:
                    raise IndexError(str(iofile.filename)+' '+str(e))
                # should contain variables
                # e.g. temperature, temperature_an_depar
                variables = list(data.keys())
            
                
                    
                data = xr.Dataset(data)
                # run adjustment procedure
                data = adjustment_procedure(data,
                                            obs_name=variables[0],
                                            dep_name=variables[1],
                                            metadata=True,
                                            times=[0, 12],
                                            dim='time',
                                            plev='plev',
                                            return_dataset=False,
                                            quantile_adjustments=True,
                                            )
                if data is None:
                    return
                # TODO Convert adjustments to other variables?
                #
                #
                # Write back adjusted (interpolation, extrapolation)
                #
            except ValueError as e:
                logger.warning('could not adjust '+ os.path.basename(iofile.filename) )
                print(e)
                return
            
            
            tt = time.time()
            try:
                xyz = add_adj(iofile, data['adjustments'],'humidity_bias_estimate', eua.cdm_codes[variable], 'r+')
                print(f'add_adj: {time.time() - tt:.3f}s')
                #iofile.write_observed_data('humidity_bias_estimate2',
                                           #varnum=eua.cdm_codes[variable],
                                           #ragged=xyz,
                                           #group=args.homogenisation,
                                           #interpolate=args.interpolate_missing,
                                           #interpolate_datetime=args.interpolate_missing,
                                           #extrapolate_plevs=args.interpolate_missing)
                
                #print(f'wrote xyz: {time.time() - tt:.3f}s')
                #iofile.write_observed_data('humidity_bias_estimate',
                                           #varnum=eua.cdm_codes[variable],
                                           #cube=data['adjustments'],
                                           #group=args.homogenisation,
                                           #interpolate=args.interpolate_missing,
                                           #interpolate_datetime=args.interpolate_missing,
                                           #extrapolate_plevs=args.interpolate_missing)
                
                print('wrote', variable, iofile.filename, args.interpolate_missing, time.time() - tt)
                printstats(data, 'relative_humidity_fg_depar', 'adjustments', np.nanstd, tstart=args.tstart, tstop=args.tstop)
                #for ip in range(16):
                    #print(f'{int(data.plev[ip]): >6d}: {np.nanstd(data.relative_humidity_fg_depar[:, ip]):.4f}, {np.nanstd(data.relative_humidity_fg_depar[:, ip]-data.adjustments[:, ip]):.4f}, {(np.nanstd(data.relative_humidity_fg_depar[:, ip]-data.adjustments[:, ip])) / np.nanstd(data.relative_humidity_fg_depar[:, ip]):.4f}')

            except MemoryError as e:
                print(iofile.filename, e, 'could not write humidity',time.time() - tt)
            
            # add also dewpoint homogenisation:
            # read data:
            try:
                
                rh = iofile.read_data_to_cube('relative_humidity',
                                                dates=args.dates,
                                                plevs=args.plevs,
                                                feedback='humidity_bias_estimate',
                                                feedback_group=args.homogenisation,
                                               )
    
                dp = iofile.read_data_to_cube('dew_point_temperature',
                                                dates=args.dates,
                                                plevs=args.plevs,
                                                feedback=None,
                                                feedback_group=None,
                                               )
                dp['dew_point_temperature'].values[dp['dew_point_temperature'].values<160] = np.nan
                ta = iofile.read_data_to_cube('temperature',
                                                dates=args.dates,
                                                plevs=args.plevs,
                                                feedback=None,
                                                feedback_group=None,
                                             )
                sh = iofile.read_data_to_cube('specific_humidity',
                                                dates=args.dates,
                                                plevs=args.plevs,
                                                feedback=None,
                                                feedback_group=None,
                                             )
                # combine cube
                comb = xr.combine_by_coords([ta['temperature'], dp['dew_point_temperature'], sh['specific_humidity']])
                comb = xr.combine_by_coords([comb, rh['relative_humidity']])
                comb = xr.combine_by_coords([comb, rh['relative_humidity_humidity_bias_estimate'].rename('relative_humidity_humidity_bias_estimate')])
                
                # convert adjusted relative humidity to adjusted dewpoint temperature
                comb['dew_point_depression_adjusted'] = rasotools.met.convert.to_dpd(temp=comb['temperature'], 
                                                                                      press=comb['plev'],
                                                                                      rel_humi=comb['relative_humidity']-comb['relative_humidity_humidity_bias_estimate'],
                                                                                      svp_method='Sonntag',
                                                                                     )
                comb['dew_point_depression'] = rasotools.met.convert.to_dpd(temp=comb['temperature'], 
                                                                                      press=comb['plev'],
                                                                                      rel_humi=comb['relative_humidity'],
                                                                                      svp_method='Sonntag',
                                                                                     )
                # create bias estimates
                comb['dew_point_temperature_humidity_bias_estimate'] = comb['dew_point_temperature'] - (ta['temperature']- comb['dew_point_depression_adjusted'])
                comb['dew_point_depression_humidity_bias_estimate'] = - comb['dew_point_temperature_humidity_bias_estimate']
                
                # write to file, into the same variable, at different lines
                xyz = add_adj(iofile, comb['dew_point_temperature_humidity_bias_estimate'],'humidity_bias_estimate', eua.cdm_codes['dew_point_temperature'], 'r+')
                #iofile.write_observed_data('humidity_bias_estimate',
                                           #varnum=eua.cdm_codes['dew_point_temperature'],
                                           #cube=comb['dew_point_temperature_humidity_bias_estimate'],
                                           #group=args.homogenisation,
                                           #interpolate=args.interpolate_missing,
                                           #interpolate_datetime=args.interpolate_missing,
                                           #extrapolate_plevs=args.interpolate_missing)
                print('wrote dp', iofile.filename, args.interpolate_missing, time.time()-tt)
                printstats(comb, 'dew_point_temperature', 'dew_point_temperature_humidity_bias_estimate', rms, tstart=args.tstart, tstop=args.tstop)
                # write to file, into the same variable, at different lines
                xyz = add_adj(iofile, comb['dew_point_depression_humidity_bias_estimate'],'humidity_bias_estimate', eua.cdm_codes['dew_point_depression'], 'r+')
                #iofile.write_observed_data('humidity_bias_estimate',
                                           #varnum=eua.cdm_codes['dew_point_depression'],
                                           #cube=comb['dew_point_depression_humidity_bias_estimate'],
                                           #group=args.homogenisation,
                                           #interpolate=args.interpolate_missing,
                                           #interpolate_datetime=args.interpolate_missing,
                                           #extrapolate_plevs=args.interpolate_missing)
                print('wrote dpd', iofile.filename, args.interpolate_missing, time.time()-tt)
                printstats(comb, 'dew_point_depression', 'dew_point_depression_humidity_bias_estimate', rms, tstart=args.tstart, tstop=args.tstop)
            except MemoryError as e:
                logger.warning('could not adjust dewpoint '+ os.path.basename(iofile.filename) )
            
            try:
                
                ## convert adjusted relative humidity to adjusted dewpoint temperature
                #comb['specific_humidity_adjusted'] = rasotools.met.convert.to_sh(temp=comb['temperature'], 
                                                                                      #press=comb['plev'],
                                                                                      #rel_humi=comb['relative_humidity']-comb['relative_humidity_humidity_bias_estimate'],
                                                                                      #svp_method='Sonntag',
                                                                                     #)
                #comb['specific_humidity2'] = rasotools.met.convert.to_sh(temp=comb['temperature'], 
                                                                                      #press=comb['plev'],
                                                                                      #rel_humi=comb['relative_humidity'],
                                                                                      #svp_method='Sonntag',
                                                                                     #)
                ## create bias estimates
                #comb['specific_humidity_bias_estimate'] = comb['specific_humidity'] - comb['specific_humidity_adjusted']
                #comb['specific_humidity_bias_estimate'].values[ ~np.isnan(comb['specific_humidity_bias_estimate'].values)] = 0.
                
                # write to file, into the same variable, at different lines
                xyz = add_adj(iofile, data.qadjustments,'humidity_bias_estimate', eua.cdm_codes['specific_humidity'], 'r+')
                #iofile.write_observed_data('humidity_bias_estimate',
                                           #varnum=eua.cdm_codes['specific_humidity'],
                                           #cube=data.qadjustments, #comb['specific_humidity_bias_estimate'],
                                           #group=args.homogenisation,
                                           #interpolate=args.interpolate_missing,
                                           #interpolate_datetime=args.interpolate_missing,
                                           #extrapolate_plevs=args.interpolate_missing)
                print('wrote sh', iofile.filename, args.interpolate_missing, time.time()-tt)
                printstats(data, 'specific_humidity_fg_depar', 'qadjustments', rms, tstart=args.tstart, tstop=args.tstop)
        
                #for ip in range(16):
                    #print(f'{int(data.plev[ip]): >6d}: {np.nanstd(data.specific_humidity_fg_depar[:, ip]):.4f}, {np.nanstd(data.specific_humidity_fg_depar[:, ip]-comb['specific_humidity_bias_estimate'][:, ip]):.4f},'+
                          #f'{(np.nanstd(data.specific_humidity_fg_depar[:, ip]-comb['specific_humidity_bias_estimate'][:, ip])) / np.nanstd(data.specific_humidity_fg_depar[:, ip]):.4f}')
            except MemoryError:
                logger.warning('could not adjust specific humidity '+ os.path.basename(iofile.filename) )
                
            with open(iofile.filename+'.x', 'w') as f:
                f.write('')
            print(iofile.filename, 'finished', time.time()-tt)
        else:
            print(iofile.filename, 'already processed')
            pass
        
        
        

    if args.winds:
        # Code 106 (wind_direction), 107 (wind_speed)
        variable = 'wind_direction'
        variable = 'wind_speed'
        raise NotImplementedError()

ray_run_backend_file = ray.remote(run_backend_file)

# -----------------------------------------------------------------------------
#
# Script entry point
#
# -----------------------------------------------------------------------------


if __name__ == "__main__":
    import sys
    import argparse
    import glob

    # handle arguments
    kwargs = {'verbose': 1}
    epilog = """
Keyword options for Standard Normal Homogeneity Test (SNHT):
    --thres []          Threshold value for SNHT, default: 50
    --window []         Moving Window for SNHT, default: 1470 (in days, 4 years)
    --missing []        Maximum allowed missing values in window, default: 600 (in days)
    --min_levels []     Minimum required levels for significant breakpoint, default: 3
    --dist []           Minimum distance between breakpoints, default: 730 (in days, 2 years)
    --sample_size []    Minimum sample size for statistics, default: 130 (in days)
    --borders []        Breakpoint zone, default: 90 (in days)
    --ratio []          Use ratio instead of differences, default: 0 (not)

Examples:
    This will run the humidity bias estimation routine and write back the bias estimates into a group called
    advanced_homogenisation.

    >>> raso_adj_cdm_v1.py -b 0-20000-0-67001_CEUAS_merged_v1.nc --humidity 

version: {}
author: {}
date: {}
---------------------------------------------
""".format(__version__, __author__, __date__)

    parser = argparse.ArgumentParser(description="Run standardized radiosonde homogenisation software on CDM compliant file",
    formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=epilog)
    parser.add_argument("-f", "--frontend", help="CDM compliant file")
    parser.add_argument("-b", "--backend", help="CDM raw file")
    parser.add_argument("-o", "--outdir", help="Output directory")
    parser.add_argument("-d","--dates", help="datetime selection, e.g. 2000-01-01,2020-12-31")
    parser.add_argument("-p","--plevs", help="pressure level selection, e.g. 100,200,300,500,700,850")
    parser.add_argument("--temperature", help="run adjustment on temperatures", action="store_true")
    parser.add_argument("--humidity", help="run adjustment on humidities", action="store_true")
    parser.add_argument("--winds", help="run adjustment on winds", action="store_true")
    parser.add_argument("--dewpoint", help="run adjustment on dewpoint", action="store_true")
    parser.add_argument("--feedback", help="feedback variables")
    parser.add_argument("--feedback_group", help="feedback group name (only backend files)", default='era5fb')
    parser.add_argument("--homogenisation", help="homogenisation group name (only backend files)", default='advanced_homogenisation')
    parser.add_argument("--interpolate_missing", help="interpolate missing values", action="store_true", default=True)
    parser.add_argument("--copypart", help="copy only partial backendfile", action="store_true")
    parser.add_argument("--debug", help="debug information", action="store_true")
    parser.add_argument("--verbose", help="show more information", action="store_true")
    parser.add_argument("--logfile", help="Logfile", default='adjustments.log')
    parser.add_argument("--tstart", help="start time for calculating statistics", default='1900-01-01')
    parser.add_argument("--tstop", help="stop time for calculating statistics", default='2100-01-01')
    
    kwargs = {}
    # Parse Arguments from definition
    args, unknown = parser.parse_known_args()
    # Check for Unknown arguments
    i = 0 
    n = len(unknown)
    for iarg in unknown:
        if '--' in iarg:
            if i+1 < n:
                if '--' not in unknown[i+1]:
                    kwargs[iarg[2:]] = eval(unknown[i+1])
                    i+=1
                else:
                    kwargs[iarg[2:]] = True
            else:
                kwargs[iarg[2:]] = True
        i+=1
    # Inputs
    if args.dates:
        args.dates = args.dates.split(',') if ',' in args.dates else args.dates
    if args.plevs:
        args.plevs = args.plevs.split(',') if ',' in args.plevs else args.plevs
    # LOGGING 
    # Clear any existing loggers
    if (logger.hasHandlers()):
        logger.handlers.clear()
    
    ch = logging.StreamHandler() # to std.err
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    if args.debug:
        logger.setLevel(logging.DEBUG)
        ch.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO if args.verbose else logging.WARNING)
        ch.setLevel(logging.INFO if args.verbose else logging.WARNING)

    # LOG TO FILE?
    if args.logfile:
        ch = logging.FileHandler(args.logfile)
        ch.setLevel(logging.DEBUG)
        ch.setFormatter(formatter)
        logger.addHandler(ch)
        logger.info("Logging to: %s", args.logfile)
    # END LOGGING
    #
    # Check file
    #
    if not (args.backend or args.frontend):
        parser.print_help()
        logger.error("Missing input file (-f, -b)")
        sys.exit(0)
    #
    # Check input variable
    #
    if not (args.temperature or args.humidity or args.winds or args.dewpoint):
        logger.error('Please specify at least one option: --temperature, --humidty or --winds')
        sys.exit(0)
    #
    # Pressure levels
    #
    if args.plevs:
        args.plevs = list(map(int, args.plevs))
    else:
        args.plevs = std_plevs * 100



    if args.frontend:
        #
        # FRONTEND File
        #
        run_frontend_file(args, **kwargs)
    else:
        #
        # BACKEND File
        #
        #run_backend_file(args, **kwargs)
        fns = glob.glob(os.path.expandvars('$RSCRATCH/converted_v29/long/*v3.nc'))
        #fns = glob.glob(os.path.expandvars('$RSCRATCH/converted_v11/long/*v1.nc'))
        tt = time.time()
        #run_backend_file(copy.copy(args), **kwargs)
        #print(time.time() - tt)
        #exit()

        #ray.init(num_cpus=60)
        futures = []
        #fns = [fns[fns.index('/mnt/users/scratch/leo/scratch/converted_v21/long/0-20000-0-72357_CEUAS_merged_v3.nc')]]
        #fns = [fns[fns.index('/mnt/users/scratch/leo/scratch/converted_v19/long/0-20666-0-1:25428_CEUAS_merged_v3.nc')],
               #fns[fns.index('/mnt/users/scratch/leo/scratch/converted_v19/long/0-20666-0-20353_CEUAS_merged_v3.nc')]]
               
        for fn in fns[:]:
            fargs = copy.copy(args)
            fargs.backend = fn
            #if '0-20000-0-17607' in fn:
            #run_backend_file(fargs, **kwargs)
            futures.append(ray_run_backend_file.remote(fargs, **kwargs))
        ray.get(futures)
        print('finished', time.time() - tt)
# FIN
