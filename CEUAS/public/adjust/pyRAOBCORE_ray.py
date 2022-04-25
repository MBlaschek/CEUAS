#!/usr/bin/env python
# coding: utf-8

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 
import numpy
import numpy as np
from numba import config
config.THREADING_LAYER = 'safe'

import math
import time
import datetime
import h5py
import netCDF4
import matplotlib.pylab as plt
import os,sys,glob,psutil,shutil
sys.path.append(os.getcwd()+'/../adjust/rasotools/')
sys.path.append(os.path.expanduser('~/python/Rasotools/rasotools/'))
sys.path.append(os.path.expanduser('~/python/Rasotools/'))
#sys.path.append(os.getcwd()+'/../adjust/rasotools/')
from utils import tdist,extract,rmeanw
# import ftest
from multiprocessing import Pool
#import odb
#from eccodes import *
from functools import partial
#from collections import OrderedDict
import subprocess
import json
import gzip
# from retrieve_fb_jra55 import add_feedback
import copy
import pickle
#import xarray as xr

import ray

import pandas as pd
import matplotlib
import matplotlib.pylab as plt
import matplotlib.pyplot as maplt
#matplotlib.rcParams.update({'font.size': 16})

#plt.rcParams['lines.linewidth'] = 3

#RAOBCORE constants
try:
    
    with open('RCx.json') as f:
        json.load(f,RC)
except:
    
    RC=dict(
        snht_maxlen = 1460,
        snht_maxmiss = 650,
        snht_increment = 30,
        miss_val=math.nan,
        stdplevs = (10.0, 20.0, 30.0, 50.0, 70.0, 100.0, 150.0, 200.0, 
                             250.0, 300.0, 400.0, 500.0, 700.0, 850.0, 925.0, 1000.0),
        pidx= (0,1,2,3,4,5,6,7,8,9,10,11,12,13),
        )
    
    RC['refdate']=(1900,1,1)
    refdate=datetime.datetime(*RC['refdate'])
    lastdate=datetime.datetime(2023,1,1)
    numdays=(lastdate-refdate).days
    dates = [refdate + datetime.timedelta(days=x) for x in range(numdays)]
    RC['years']= [ x.year for x in dates]
    RC['months']= [ x.month for x in dates]
    RC['days']= [ x.day for x in dates]
    RC['snht_thresh'] = 25
    RC['mon_thresh'] = 2
    RC['mean_maxlen'] = RC['snht_maxlen']*2
    RC['plot']=False
    RC['CPUs']=20
    #RC['transdict']={'montemp':'mtemperatures','rasocorrmon':'madjustments','goodmon':'goodmon'}
    RC['transdict']={'montemp':'madjusted_temperatures','rasocorrmon':'madjustments','goodmon':'goodmon'}
    
    
    with open('RC.json','w') as f:
        json.dump(RC,f)
for k in 'years','months','days','stdplevs','pidx':
    if type(k) is int:
        RC[k]=np.array(RC[k],dtype=np.int32)
    else:
        RC[k]=np.array(RC[k])
    if k in  ('months','days'):
        RC[k]-=1

RC['plot']=True
RC['savetsa']=True

#sys.path.append(os.getcwd()+'/../cds-backend/code/')
#import cds_eua3 as eua

#import zarr
#import dask
from timeit import default_timer as timer
from numba import njit,prange,typeof
import SNHT

# ## Translated hom-functions

# In[2]:


def testeq(x, window, missing):
    """Standard Normal Homogeneity Test (SNHT) with a running window
    Wrapper function for numba_snhtmov

    Args:
        x (np.ndarray) : input data
        window (int) : window size (in days)
        missing (int) : allowed missing values (in days)
    Returns:
        np.ndarray : SNHT
    """

    tsa = np.zeros_like(x)
    tmean = np.zeros_like(x,shape=(x.shape[0],12))
    tsquare = np.zeros_like(x,shape=(x.shape[0],12))
    count = np.zeros((x.shape[0],12), dtype=np.int32)
    rc_months=RC['months'][:x.shape[0]]

    bmax = SNHT.numba_snhteqmov(x,rc_months,
                       tsa,
                       RC['snht_maxlen'],
                       RC['snht_maxmiss'],
                       RC['miss_val'],
                       count,
                       tmean,
                       tsquare)
    return tsa


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

    tsa = SNHT.numba_snhtmov(np.squeeze(np.asarray(x)),
                       tsa,
                       snhtparas,
                       count,
                       tmean,
                       tsquare)
    return tsa

def prolog():
    
    a = 5
    nodec = np.concatenate(([1]*300, [np.nan]*65))
    allnodec = nodec
    for i in range(a-1):
        allnodec = np.concatenate((allnodec, nodec))
        
    
    before_break = np.cos(np.linspace(0,np.pi*2*a,365*a))*2 + np.random.randn(365*a) * 1 + 240
    after_break = np.cos(np.linspace(0,np.pi*2*a,365*a))*2 + np.random.randn(365*a) * 1 + 240.0
    ts = np.concatenate((before_break*allnodec, after_break), axis=0)
    
    
    
    # before_break = np.random.rand(365*a) * 15
    # after_break = np.random.rand(365*a) * 15 + 2
    # test = np.concatenate((before_break, after_break), axis=0)
    
    print(len(ts))
    fig, ax = maplt.subplots(1, figsize = (15,5))
    ax.plot(np.array(range(len(ts))), ts, color = 'blue', alpha = 0.3, label='test' )
    ax.legend()
    ax.grid()
    #maplt.show()
    maplt.close()
    
    ni = len(ts)
    
    ref = np.zeros(ni)
    
    maxlen = left_maxlen = right_maxlen = 4*365
    
    max_miss = 650 #(or more)
    istart = istartorig = 0 
    istop = iststoporig = len(ts)
    increment = 30
    miss_val = np.nan #(-999.99)
    
    
    pcount=np.zeros(ni,dtype=np.int32)
    mcount=np.zeros(ni,dtype=np.int32)
    tsa = np.zeros(ni,dtype=np.float32)
    plus = np.zeros(ni,dtype=np.float32)
    minus = np.zeros(ni,dtype=np.float32)
    prms = np.zeros(ni,dtype=np.float32)
    mrms = np.zeros(ni,dtype=np.float32)
    
    SNHT.snhteqsamp(np.array(ts),np.array(ref),ni,istart,istop,maxlen,increment,miss_val,max_miss,
               tsa,plus,minus,prms,mrms,pcount,mcount)
    
    tts=np.empty_like(ts,dtype=np.float32)
    tts[:]=ts[:]
    #tts[np.isnan(tts)]=-999.
    tsa2=testeq(tts, maxlen, max_miss)
    tsa3=test(ts, maxlen, max_miss)
    
    tt=time.time()
    pcount=np.zeros(ni,dtype=np.int32)
    mcount=np.zeros(ni,dtype=np.int32)
    tsa = np.zeros(ni,dtype=np.float32)
    plus = np.zeros(ni,dtype=np.float32)
    minus = np.zeros(ni,dtype=np.float32)
    prms = np.zeros(ni,dtype=np.float32)
    mrms = np.zeros(ni,dtype=np.float32)
    SNHT.snhteqsamp(np.array(ts),np.array(ref),ni,istart,istop,maxlen,increment,miss_val,max_miss,
               tsa,plus,minus,prms,mrms,pcount,mcount)
    print('classic',time.time()-tt)
    tt=time.time()
    rc_months=RC['months'][:tts.shape[0]]
    snhtparas = np.asarray([maxlen, max_miss, 10])
    tsa = np.zeros(tts.shape[0],dtype=np.float32)
    tmean = np.zeros((tts.shape[0],12),dtype=np.float32)
    tsquare = np.zeros((tts.shape[0],12),dtype=np.float32)
    count = np.zeros((tts.shape[0],12), dtype=np.int32)
    #rc_months=np.array(RC['months'])[:tts.shape[0]]-1
    bmax = SNHT.numba_snhteqmov(tts, rc_months,tsa2, RC['snht_maxlen'], RC['snht_maxmiss'], RC['miss_val'],
                                count,tmean,tsquare)
    
    #tsa2=testeq(tts, maxlen, max_miss)
    print('csumeq',time.time()-tt)
    tt=time.time()
    tsa3=test(ts, maxlen, max_miss)
    print('csum',time.time()-tt)
    
    
    fig=maplt.figure(figsize = (15,5))
    ax = maplt.subplot(2,1,1)
    ax.plot(np.array(range(len(tsa))),tsa, color = 'blue', alpha = 0.3, label='tsa' )
    ax.plot(np.array(range(len(tsa2))),tsa2, color = 'red', alpha = 0.3, label='tsa2' )
    ax.plot(np.array(range(len(tsa3))),tsa3, color = 'green', alpha = 0.3, label='tsa3' )
    ax.legend()
    ax = maplt.subplot(2,1,2)
    #ax.plot(np.array(range(len(plus))),plus, color = 'red', alpha = 0.3, label='plus' )
    #ax.plot(np.array(range(len(minus))),minus, color = 'green', alpha = 0.3, label='minus' )
    ax.plot(np.array(range(len(minus))),plus-minus, color = 'green', alpha = 0.3, label='minus' )
    
    ax.legend()
    #ax.grid()
    #maplt.show()
    maplt.close()

@njit
def do_test(fg_dep,rc_month,snht_maxlen,snht_maxmiss,miss_val):
    #rc_month= RC['months'][days]
    #RC['snht_maxlen'],RC['snht_maxmiss'],RC['miss_val']
    tmean = np.zeros((fg_dep.shape[2],12),dtype=fg_dep.dtype)
    tsquare = np.zeros((fg_dep.shape[2],12),dtype=fg_dep.dtype)
    count = np.zeros(shape=(fg_dep.shape[2],12),dtype=np.int32)
    tsas=np.zeros(fg_dep.shape,dtype=fg_dep.dtype)
    tsarad=np.zeros((1,6,fg_dep.shape[2]),dtype=fg_dep.dtype)
    
    for ih in range(fg_dep.shape[0]):
        for ip in range(fg_dep.shape[1]):
            SNHT.numba_snhteqmov(fg_dep[ih,ip,:].flatten(), rc_month,tsas[ih,ip,:], 
                                                 snht_maxlen,snht_maxmiss,miss_val,
                                                 count,tmean,tsquare)
            #print(ih,ip)
            
    for ip in range(6):
        SNHT.numba_snhteqmov((fg_dep[1,ip,:]-fg_dep[0,ip,:]).flatten(), rc_month, tsarad[0,ip,:], 
                             snht_maxlen,snht_maxmiss,miss_val,count,tmean,tsquare)
        #print(ip)
    
    return tsas,tsarad

@njit
def areg1(x,z,alpha):
    x[0]=z[0]
    for i in range(1,z.shape[0]):
        x[i]=alpha*x[i-1]+z[i]
    return

@njit(parallel=True)
def break_simulator(x):
    m=x.shape[0]
    p=x.shape[1]
    n=x.shape[2]
    
    tsas=np.zeros(x.shape,dtype=np.float32)
    imax=np.zeros((m,p),dtype=np.int32)
    rmax=np.zeros((m,p),dtype=np.int32)
   
    for im in range(m):
        for ip in prange(p):
            tsa = np.zeros((x.shape[2]),dtype=np.float32)
            tmean = np.zeros((x.shape[2],12),dtype=np.float32)
            tsquare = np.zeros((x.shape[2],12),dtype=np.float32)
            count = np.zeros((x.shape[2],12),dtype=np.int32)

            z=np.random.randn(n)
            areg1(x[im,ip,:],z,0.0)
            x[im,ip,n//2:]+=im*0.2
            tsas[im,ip,:]=SNHT.numba_snhteqmov(x[im,ip,:].flatten(), tsa, 
                                               RC['snht_maxlen'],RC['snht_maxmiss'],RC['miss_val'],
                                               count,tmean,tsquare)
            imax[im,ip]=np.argmax(tsas[im,ip,:])
            rmax[im,ip]=np.max(tsas[im,ip,:])
        print(im)
    
    return imax,rmax,tsas  
    
def break_analyze(finfo):
    
    breaklist=np.zeros(50,dtype=np.int32)
    #plt.plot(finfo['days']/365.25,np.sum(finfo['tsas'][0,:,:],axis=0))
    #plt.plot(finfo['days']/365.25,np.sum(finfo['tsas'][1,:,:],axis=0))
    #plt.plot(finfo['days']/365.25,np.sum(finfo['tsarad'][0,:,:],axis=0))
    
    totalo=np.sum(finfo['tsarad'][0,:,:],axis=0)+np.sum(finfo['tsas'][0,:,:],axis=0)+np.sum(finfo['tsas'][1,:,:],axis=0)
    if RC['plot']:
        plt.plot(finfo['days']/365.25,totalo)
    total=np.sum(finfo['tsarad'][0,:,:],axis=0)+np.sum(finfo['tsas'][0,:,:],axis=0)+np.sum(finfo['tsas'][1,:,:],axis=0)
    l=0

    while np.max(total)>RC['snht_thresh']*5:
        breaklist[l]=np.argmax(total)
        lb=max(breaklist[l]-RC['snht_maxlen']//2,0)
        rb=min(total.shape[0],breaklist[l]+RC['snht_maxlen']//2)
        total[lb:rb]=0.
        l+=1
        #print(l,np.max(total))
        #plt.plot(finfo['days']/365.25,total)
        #plt.plot(finfo['days'][breaklist[:l]]/365.25,total[breaklist[:l]],'k*')
        #plt.show()
    if RC['plot']:
        plt.plot(finfo['days'][breaklist[:l]]/365.25,totalo[breaklist[:l]],'k*')
        plt.savefig(finfo['sid']+'_breakanalysis.png')
        plt.close()
    if RC['savetsa']:
        with open(finfo['sid']+'_breakanalysis.json','w') as f:
            json.dump(dict(days_since_1900=[int(d) for d in finfo['days']],
                           tsasum=[float(d) for d in totalo],
                           breaks=[int(d) for d in breaklist[:l]]),f)
        with open(finfo['sid']+'_breakanalysis.json') as f:
            d=json.load(f)
    
    return breaklist[:l]
    
    
@ray.remote
def RAOB_findbreaks(method,fn):
    tt=time.time()
    finfo={}

    with h5py.File(fn,'r') as f:
        if f['datum'].shape[0]<RC['snht_maxlen']:
            return finfo
        finfo['days']=f['datum'][:]
        print(fn.split('/')[-1],finfo['days'].shape)
        finfo['fg_dep']=f['era5_fgdep'][:,np.array(RC['pidx']),:]
        finfo['temperatures']=f['temperatures'][:,np.array(RC['pidx']),:]
        finfo['lon']=f['lon'][0]
        finfo['lat']=f['lat'][0]
        finfo['sid']=fn.split('/')[-1][:-3]
    
    print('before snht: {:6.4f}'.format(time.time()-tt)) 
    #print(typeof(finfo['fg_dep']),typeof(RC['months'][finfo['days']]),
          #typeof(RC['snht_maxlen']),typeof(RC['snht_maxmiss']),typeof(RC['miss_val']))
    finfo['tsas'],finfo['tsarad']=do_test(finfo['fg_dep'],RC['months'][finfo['days']],
                                          RC['snht_maxlen'],RC['snht_maxmiss'],RC['miss_val'])  
    
    print('before analyze: {:6.4f}'.format(time.time()-tt))  
    finfo['breaklist']=np.sort(break_analyze(finfo))
    
    print(finfo['sid']+' findbreaks: {:6.4f}'.format(time.time()-tt))  
    obj_ref=ray.put(finfo)
    return obj_ref

def RAOB_adjust(finfo):
    
    sys.path.append(os.path.expanduser('~/python/Rasotools/rasotools/'))
    sys.path.append(os.path.expanduser('~/python/Rasotools/'))
    from utils import tdist,extract,rmeanw
    tt=time.time()
    fg_dep=np.empty_like(finfo['fg_dep'])
    fg_dep[:]=finfo['fg_dep'][:]
    mask=np.isnan(fg_dep)
    tmean = np.zeros_like(fg_dep,shape=(fg_dep.shape[2],12))
    tsquare = np.zeros_like(fg_dep,shape=(fg_dep.shape[2],12))
    count = np.zeros(shape=(fg_dep.shape[2],12),dtype=np.int32)
    tsa=np.zeros_like(fg_dep,shape=(fg_dep.shape[2]))

    finfo['adjustments']=np.zeros_like(fg_dep)
    finfo['adjusted_temperatures']=np.empty_like(fg_dep)
    finfo['adjusted_temperatures'][:]=finfo['temperatures'][:]
    
    sh=finfo['adjustments'].shape
    fb=finfo['breaklist']
    RC_months=np.array(RC['months'])[finfo['days']]-1
    if RC['plot']:
        plt.plot(finfo['days']/365.25,rmeanw(finfo['fg_dep'][0,5,:],30))
    
    break_profiles=np.zeros_like(fg_dep,shape=(len(fb),sh[0],sh[1]))
    for ib in range(len(fb)-1,-1,-1):
        if ib>0:          
            istart=np.max((fb[ib]-RC['mean_maxlen'],fb[ib-1]))
        else:
            istart=np.max((fb[ib]-RC['mean_maxlen'],0))
        istop=np.min((fb[ib]+RC['mean_maxlen'],sh[2]))

        for ih in range(sh[0]):
            for ip in range(sh[1]):
                break_profiles[ib,ih,ip]=SNHT.numba_meaneqmov(fg_dep[ih,ip,istart:istop], 
                                                          RC_months[istart:istop], tsa, 
                                                          RC['mean_maxlen'],RC['mean_maxlen']//2-80,RC['miss_val'], 
                                                          count, tmean, tsquare, kref=fb[ib]-istart)
                
                if ~np.isnan(break_profiles[ib,ih,ip]):
                    finfo['adjustments'][ih,ip,:fb[ib]]-=break_profiles[ib,ih,ip]
                    fg_dep[ih,ip,:fb[ib]]+=break_profiles[ib,ih,ip]
                    finfo['adjusted_temperatures'][ih,ip,:fb[ib]]-=break_profiles[ib,ih,ip]
        #print(ib,fb[ib],break_profiles[ib,:,:])
    
    print('adjust: {:6.4f}'.format(time.time()-tt))  
    finfo['adjusted_temperatures'][mask]=np.nan
    finfo['adjustments'][mask]=np.nan
    if RC['plot']:
        plt.plot(finfo['days']/365.25,rmeanw(fg_dep[0,5,:],30))
        plt.plot(finfo['days']/365.25,finfo['adjustments'][0,5,:])
        plt.savefig(finfo['sid']+'_adjustments.png')
        #plt.show()
        plt.close()
        
    return ray.put((finfo['adjustments'],finfo['adjusted_temperatures']))
        
        
    
@ray.remote    
def RAOB_adjustbreaks(method,finfo):
    
    adjustments=RAOB_adjust(finfo)
    
    return adjustments
    
def fcopy(fn,exper='exp02'):
    fo=exper.join(fn.split('exp03'))
    if not os.path.isdir(os.path.dirname(fo)):
        os.mkdir(os.path.dirname(fo))
    
    try:
        
        shutil.copy(fn,fo)
    except:
        print(fn,'could not be copied to',fo)
    
    return

@njit(fastmath={'nsz','arcp','contract','afn','reassoc'},cache=True)
def goodmon(var,idx):

    l=0
    for i in range(len(idx)-1):
        if idx[i+1]>idx[i]:
            l+=1
    if l==0:
        print('idx not valid')
    out=np.empty((var.shape[0],var.shape[1],l),dtype=np.int32)
    for ih in range(var.shape[0]):
        for ip in range(var.shape[1]):
            l=0
            for i in range(len(idx)-1):
                if idx[i+1]>idx[i]:
                    mc=0
                    for j in range(idx[i],idx[i+1]):
                        v=var[ih,ip,j]
                        if not np.isnan(var[ih,ip,j]):
                            mc+=1
                    out[ih,ip,l]=mc
                    l+=1
                
    return out
    
@njit(fastmath={'nsz','arcp','contract','afn','reassoc'},cache=True)
def monmean(var,idx,thresh=15,goodmon=goodmon):
    l=0
    for i in range(len(idx)-1):
        if idx[i+1]>idx[i]:
            l+=1
    if l==0:
        print('idx not valid')
    out=np.empty((var.shape[0],var.shape[1],l),dtype=var.dtype)
    for ih in range(var.shape[0]):
        for ip in range(var.shape[1]):
            l=0
            for i in range(len(idx)-1):
                if idx[i+1]>idx[i]:
                    mc=0
                    mm=np.float32(0.)
                    for j in range(idx[i],idx[i+1]):
                        v=var[ih,ip,j]
                        if v==v:
                            mm+=v
                            mc+=1
                    if mc>thresh:
                        out[ih,ip,l]=mm/mc
                    else:
                        out[ih,ip,l]=np.nan
                    l+=1
                
    return out

#@njit
def do_monmean(vars,tidx,fi):
    
    mdays=np.zeros_like(tidx)
    dum=np.zeros_like(fi['temperatures'],shape=tidx.shape)
    idx=np.searchsorted(fi['days'],tidx)
    mdays=idx[1:]-idx[:-1]
    midx=np.where(mdays>0)[0]
    
    fi['months']=tidx[midx]
            
    for v in vars:
        fi['m'+v]=np.empty_like(fi[v],shape=(fi[v].shape[0],fi[v].shape[1],mdays.shape[0]))
        fi['m'+v]=monmean(fi[v],idx,RC['mon_thresh'])
        
    fi['goodmon']=goodmon(fi['temperatures'],idx)
            
                
    return

def write_monmean(ipath,opath,prefix,fi):
    
    if not os.path.isdir(os.path.dirname(opath)):
        os.mkdir(os.path.dirname(opath))
    
    fno=opath+'/'+prefix+fi['sid'][-6:]+'.nc'
    fni=ipath+'/'+prefix+fi['sid'][-6:]+'.nc'
    try:
        with netCDF4.Dataset(fni) as src, netCDF4.Dataset(fno, "w") as dst:
            # copy global attributes all at once via dictionary
            new_globals=src.__dict__
            new_globals['history']=datetime.datetime.today().strftime("%m/%d/%Y")
            new_globals['source']='ERA5, IGRA2, NCAR UADB'
            new_globals['references']='Copernicus Early Upper Air Dataset'
            dst.setncatts(new_globals)
            # copy dimensions
            for name, dimension in src.dimensions.items():
                if name =='time':
                    dst.createDimension(
                        name, (fi['mtemperatures'].shape[2] if not dimension.isunlimited() else None))
                else:
                    dst.createDimension(
                        name, (len(dimension) if not dimension.isunlimited() else None))
            # copy all file data except for the excluded
            for name, variable in src.variables.items():
                if name == 'datum':
                    x = dst.createVariable(name, variable.datatype, variable.dimensions)
                    dst[name][:] = fi['months'][:]
                    # copy variable attributes all at once via dictionary
                    dst[name].setncatts(src[name].__dict__)
                elif name in ['montemp','rasocorrmon','goodmon']: #,'goodmon','rasocorrmon','eracorrmon']:
                    x = dst.createVariable(name, variable.datatype, variable.dimensions)
                    if variable.datatype in [np.dtype('float32'),np.dtype('float64')]:
                        
                        dst[name][:]=np.nan
                    else:
                        dst[name][:]=0
                        
                    dst[name][:,:RC['pidx'].shape[0],:] = fi[RC['transdict'][name]][:]
                    # copy variable attributes all at once via dictionary
                    dst[name].setncatts(src[name].__dict__)
                    
                elif name in ['eracorrmon']:
                    pass
                else:
                    x = dst.createVariable(name, variable.datatype, variable.dimensions)
                    dst[name][:] = src[name][:]
                    # copy variable attributes all at once via dictionary
                    dst[name].setncatts(src[name].__dict__)
    except FileNotFoundError:
        print(fni,'could not be copied to',fno)
    
    return

@ray.remote    
def save_monmean(tfile,fi,res_fi):
    
    #fnames= fi['sid']
    fi['adjustments']=res_fi[0]
    fi['adjusted_temperatures']=res_fi[1]
    tidx=np.where(RC['days']==1)[0]
    exper='exp02'
    prefix='feedbackglobbincorrmon'
    ipath=fi['sid'][-6:].join(os.path.dirname(tfile).split(files[0][-9:-3]))
    do_monmean(['temperatures','adjusted_temperatures','fg_dep','adjustments'],tidx,fi)
    write_monmean(ipath,exper.join(ipath.split('exp03')),prefix,fi)
    
    
def write_era5bc(opath,fi):
    
    flag=False
    fn='../Temperature_adjustment/'+st+'/feedbackmerged'+st+'.nc'
    f = netCDF4.Dataset(fn,"r")
    fno='../Temperature_adjustment/'+st+'/ERA5bc_RAOBCORE_v'+version+'_'+st+'.nc'
    fo = netCDF4.Dataset(fno,"w", format='NETCDF4_CLASSIC')

    for i in f.ncattrs():
        if i=='history':
            setattr(fo,i,datetime.date.today().strftime("%Y/%m/%d"))
        elif i=='source':
            setattr(fo,i,'RAOBCORE/RICH v'+version+' solar elevation dependency (from 197901 onward)' )
        elif i=='title':
            setattr(fo,i,'Station daily temperature series with ERA5/20VRv3 background departure statistics and RISE bias estimates' )
        else:
            setattr(fo,i,getattr(f,i))
    for i in list(f.dimensions.keys()):
        if i=='time':
            fo.createDimension(i,s["mdatum"].shape[0])
        else:
            try:
                fo.createDimension(i,len(f.dimensions[i]))
            except:
                flag=True
                continue
    nalias=8    
    fo.createDimension('nalias',8)
    fo.createDimension('nchar',8)
    if flag:
        return
    #nogos=['flags',u's_type', u'eijra_fgdep', u'jra55_fgdep', u'jra55_andep', u'e20c_andep', u'n20c_andep', u'ce20c_andep', u'erapresat_andep']
    tobecopied=['datum','hours','lat','lon','alt','press','temperatures','an_dep',fgdepvar,'source','mergedstats']
    for i in list(f.variables.keys()):
        var=f.variables[i]
        if i=='datum':
            fo.createVariable(i,var.dtype,var.dimensions)
            fo.variables[i][:]=s["mdatum"][:]
        elif i=='hours':
            fo.createVariable(i,var.dtype,var.dimensions)
            fo.variables[i][:]=s["mhours"][:]
        elif i=='temperatures':
            fo.createVariable(i,var.dtype,var.dimensions)
            s["mtemperatures"][numpy.isnan(s["mtemperatures"])]=-999.
            fo.variables[i][:]=s["mtemperatures"][:]
        elif i==fgdepvar:
            fo.createVariable(i,var.dtype,var.dimensions)
            s["m"+fgdepname][numpy.isnan(s["m"+fgdepname])]=-999.
            fo.variables[i][:]=s["m"+fgdepname][:]
        #elif i=='an_dep':
            s["newbias"][numpy.isnan(s["newbias"])]=-999.
            s["newrichbias"][numpy.isnan(s["newrichbias"])]=-999.
            try:
                fo.createVariable('bias',var.dtype,var.dimensions)
            except:
                pass
            fo.variables['bias'][:]=s["newbias"][:]
            fo.createVariable('richbias',var.dtype,var.dimensions)
            fo.variables['richbias'][:]=s["newrichbias"][:]
        elif i in ('lon','lat','alt','press'):
            fo.createVariable(i,var.dtype,var.dimensions)
            fo.variables[i][:]=var[:]
        elif i in ('source'):
            #str_out = netCDF4.stringtochar(numpy.array(['test'], 'S4'))
            fo.createVariable(i,'S1',('time','nchar'))
            x=numpy.empty((var.shape[0]),dtype='S8')
            x.fill('        ')
            try:
                svar=var[:]
            except:
                print('could not read source, supplying BUFRDATA')
                svar=x
            str_out = netCDF4.stringtochar(numpy.asarray(svar,'S8'))
            fo.variables[i][:]=str_out
        elif i in ('mergedstats'):
            fo.createVariable(i,'S1',('time','nalias','nchar'))
            tt=time.time()
            x=numpy.empty((var.shape[0],nalias),dtype='S8')
            x.fill('        ')
            try:
                svar=var[:]
                for k in range(svar.shape[0]):
                    l=svar[k].split(',')
                    #if len(l)>1:
                        #print 'l>1'
                    for m in range(len(l)):
                        x[k,m]=l[m]+' '*(8-len(l[m]))
            except:
                print('could not read mergedstats, filling with WMO number')
                x.fill(st[1:]+'   ')

            str_out = netCDF4.stringtochar(x)
            fo.variables[i][:]=str_out
            print(('mergevar:',time.time()-tt))
        else:
            if i in tobecopied:
                print((i,'some unknown variable'))
                fo.createVariable(i,var.dtype,var.dimensions)

        for j in var.ncattrs():
            if j!='_FillValue' and j!='scale_factor' and j!='add_offset':
                if i in tobecopied:
                    if i=='an_dep':
                        setattr(fo.variables['bias'],j,getattr(var,j))
                        setattr(fo.variables['richbias'],j,getattr(var,j))
                    elif i=='datum' and j=='units':
                        setattr(fo.variables[i],j,'days since 1900-01-01 00:00:00')
                    else:
                        setattr(fo.variables[i],j,getattr(var,j))

    setattr(fo.variables[fgdepvar],'infos','obs-20CRv3 up to 1949; obs-ERA5 1950 onwards')
    setattr(fo.variables['bias'],'long_name','RAOBCORE v'+version+' RISE bias estimate')
    setattr(fo.variables['richbias'],'long_name','RICH v'+version+' RICH bias estimate')
    #setattr(fo.variables['source'],'info','Preferred source used for merged temperature record')
    #setattr(fo.variables['source'],'valid_entries','BUFRDATA: (ECMWF holdings), NCARUA(20,21,24): sources from NCAR, CHUAN2.1: ERA-CLIM(2) digitized data')
    #setattr(fo.variables['mergedstats'],'info','ODB statIDs that matched during merge - bias adjustments should be applied to those')


    fo.close()
    f.close()
    
def stdists(finfo):
    
            
    dists=np.zeros(len(finfo)*(len(finfo)+1)//2)
    lats=np.array([i['lat'] for i in finfo ])
    lons=np.array([i['lon'] for i in finfo ])
    dum=tdist(dists,lats.astype(np.float64),lons.astype(np.float64),0)
    for i in range(len(finfo)):  
        finfo[i]['stdists']=extract(dists,i,lats,lons)
    
    return ray.put([finfo[i]['stdists'] for i in range(len(finfo))])

    
if __name__ == '__main__':
    
    
    #prolog()
    process = psutil.Process(os.getpid())
    files = glob.glob('/users/staff/leo/fastscratch/rise/1.0/exp03/*/feedbackmerged*')[:]
    tt=time.time()
    # read data and analyze breakpoints
    #func=partial(RAOB_findbreaks,'SNHT')
    #with Pool(20) as P:        
        #finfo=list(P.map(func,files))
        
    #finfo=[i for i in finfo if i]
    
    ray.init(num_cpus=RC['CPUs'])
    #func=partial(RAOB_findbreaks,'SNHT')
    #with Pool(20) as P:        
        #obj_ref=list(P.map(func,files))
    

    futures = []
    for fn in files:
        futures.append(RAOB_findbreaks.remote('SNHT',fn))
    obj_ref  = ray.get(futures)
    obj_ref = [i for i in obj_ref if i]
    
    break_sum=np.sum([len(l['breaklist']) for l in ray.get(obj_ref)])
    print('Total number of breaks',break_sum)
    #obj_ref=[]
    #for li in finfo:
        #obj_ref.append(ray.put(li))

    
    sdist_ref=stdists(ray.get(obj_ref))
    
    #func=partial(RAOB_adjustbreaks,'SNHT')
    #with Pool(20) as P:        
        #adjustments=list(P.map(func,finfo))
    #for i in range(len(finfo)): 
        #finfo[i]['adjustments']=adjustments[i][0]
        #finfo[i]['adjusted_temperatures']=adjustments[i][1]
        

    
    ## copy feedbackmerged files
    ## do not forget to copy radcorpar, mergedstations.t
    #fnames= [fi['sid'] for fi in finfo]    
    futures=[]
    results_ref=[]
    for l in range(len(obj_ref)):
        futures.append(RAOB_adjustbreaks.remote('SNHT',obj_ref[l]))
    results_ref = ray.get(futures)
        
    
    

    #for fn in fnames:
        #for fi in files:
            #if fn in fi:
                #ipath=os.path.dirname(fi)
                #break
        #fcopy(ipath+'/'+fn+'.nc')
    
    #func=partial(save_monmean,files[0])    
    #with Pool(20) as P:
        #l=list(P.map(func,finfo))
    del futures ; futures=[]
    for l in range(len(obj_ref)):
        futures.append(save_monmean.remote(files[0],obj_ref[l],results_ref[l]))
    
    res = ray.get(futures)
        

        
    print('Mem [MiB]',process.memory_info().rss//1024//1024)
    print(time.time()-tt)
    ray.shutdown()
    print('end')
     
def epilog():
    

    window = 1460  # means 4 years on daily basis
    missing = 600
    
    df_night = df[dict(hour=0)]
    df_day = df[dict(hour=1)]
    
    fig, ax = maplt.subplots(len(stdplevs)*2, 1, figsize = (15,80))
    for i in range(len(stdplevs)):
        df_day_snht = df_day[dict(pressure=i)]
        df_night_snht = df_night[dict(pressure=i)]
        df1 = df_day_snht.to_dataframe()
        
        for yr in [2010, 2011]:
            for mn in [6, 7, 8]: # range(1,12,1): #
                df1 = df1[df1.datum < str(yr)+"-"+str(mn)].append(df1[df1.datum >= str(yr)+"-"+str(mn+1)])
    
    #     for yr in [2015, 2016, 2017, 2018]:
    #         df1 = df1[df1.datum < str(yr)].append(df1[df1.datum >= str(1+yr)])
    
        snht_day,nc= test(np.array(df1.era5_fgdep), window, missing)
        
    #     snht_night = test(np.array(df_night_snht.era5_fgdep), window, missing)
    #     snht_diff_dn = test(np.array(df_day_snht.era5_fgdep)-np.array(df_night_snht.era5_fgdep), window, missing)
    #     snht_diff_dn_t = test(np.array(df_day_snht.temperatures)-np.array(df_night_snht.temperatures), window, missing)
        
        ax[i*2].plot(np.array(df1.datum),snht_day,color = 'blue', alpha = 0.6, label='DAY- ' + str(stdplevs[i]) + ' hPa', )
    #     ax[i*2].plot(np.array(df1.datum),np.array(nc)*50,color = 'blue', alpha = 0.2, label='DAY- ' + str(stdplevs[i]) + ' hPa', )
        
        snht_day_corr = snht_day
        snht_day_corr[np.array(nc) == 1] = 0
        ax[i*2].plot(np.array(df1.datum),snht_day_corr,color = 'red', alpha = 0.6, label='DAY nan - ' + str(stdplevs[i]) + ' hPa', )
        ax[i*2].scatter(np.array(df1.datum),df1.era5_fgdep,color = 'green', alpha = 0.6, label='DAY era5_fgdep - ' + str(stdplevs[i]) + ' hPa', )
    #     ax[i].plot(np.array(df_night_snht.datum),snht_night,color = 'orange', alpha = 0.6, label='NIGHT- ' + str(stdplevs[i]) + ' Pa', )
    #     ax[i].plot(np.array(df_day_snht.datum),snht_diff_dn,color = 'purple', alpha = 0.6, label='DIFF - ' + str(stdplevs[i]) + ' Pa', )
    #     ax[i].plot(np.array(df_day_snht.datum),snht_diff_dn_t,color = 'red', alpha = 0.4, label='DIFF T- ' + str(stdplevs[i]) + ' Pa', )
    #     mean_snht = (np.array(snht_day) + np.array(snht_night) + np.array(snht_diff_dn) + np.array(snht_diff_dn_t))/4.
    #     ax[i].plot(np.array(df_day_snht.datum),mean_snht,color = 'black', alpha = 0.3, linewidth = 12, label='MEAN - ' + str(stdplevs[i]) + ' Pa', )
        ax[i*2].set_ylabel('SNHT')
        ax[i*2].set_xlabel('time')
        ax[i*2].legend(loc='center left')
        ax[i*2].grid()
        
        snht_day_o,nc= test(np.array(df_day_snht.era5_fgdep), window, missing)
        ax[i*2+1].plot(np.array(df_day_snht.datum),snht_day_o,color = 'blue', alpha = 0.6, label='DAY- ' + str(stdplevs[i]) + ' hPa', )
        ax[i*2+1].set_ylabel('SNHT')
        ax[i*2+1].set_xlabel('time')
        ax[i*2+1].legend(loc='center left')
        ax[i*2+1].grid()
        
    maplt.show()
    maplt.close()
    
    files = glob.glob('/users/staff/leo/fastscratch/rise/1.0/exp03/011035/feedbackmerged*')
    func=partial(RAOB_findbreaks,'SNHT')
    breaks=list(map(func,files))
    
    
    exit()
    # ## ACMANT4 Input
    
    # In[ ]:
    
    
    # creating input files for ACMANT4
    
    testdf = df_day.to_dataframe()
    snumb = 1
    
    from datetime import date, timedelta
    
    def daterange(start_date, end_date):
        for n in range(int((end_date - start_date).days)):
            yield start_date + timedelta(n)
            
    def last_day_of_month(any_day):
        next_month = any_day.replace(day=28) + datetime.timedelta(days=4)
        return next_month - datetime.timedelta(days=next_month.day)
    
    for i in testdf.press.drop_duplicates():
        pdf = testdf[testdf.press == i]
        lines = []
        lines.append(str(df.unique_source_identifier) + '\n')
        for yr in range(1950,2022,1):
            print(yr)
            for mn in range(1,13,1):
                line = str(yr) + '\t' + str(mn)
                start_date = date(yr, mn, 1)
                end_date = last_day_of_month(start_date)
                for single_date in daterange(start_date, end_date + datetime.timedelta(days=1)):
                    try:
                        temp = float(pdf[pdf.datum == single_date.strftime("%Y-%m-%d")].temperatures)
                        if temp == temp:
                            line = line + '\t' + str(temp-273.15)[:6]
                        else:
                            line = line + '\t' + '-999.9'
                    except:
                        line = line + '\t' + '-999.9'
                lines.append(line + '\n')
        with open('S'+str(snumb).rjust(4, '0')+str(int(i)).rjust(4, '0')+"t.txt","w+") as f:
            f.writelines(lines)
    
    
    # In[ ]:
    
    
    with open("text.txt","w+") as f:
        f.writelines(lines)
    
    
    # In[13]:
    
    
    str(int(50.0)).rjust(5, '0')
    
     
    
    with open('nearest_stations.p', "rb") as input_file:
        nearest_stations = pickle.load(input_file)
    
    
    # In[184]:
    
    
    with open('nearest_stations_with_snht.p', "rb") as input_file:
        nearest_stations_with_snht = pickle.load(input_file)
    
    
    # In[186]:
    
    
    vie = glob.glob('/users/staff/leo/fastscratch/rise/1.0/exp03/*11035*/feedbackmerged*')[0]
    print(vie)
    print()
    display(nearest_stations[vie])
    print()
    display(nearest_stations_with_snht[vie])
    
    

    
    stdplevs = [10.0, 20.0, 30.0, 50.0, 70.0, 100.0, 150.0, 200.0, 250.0, 300.0, 400.0, 500.0, 700.0, 850.0, 925.0, 1000.0]
    window = 1460  # means 4 years on daily basis
    missing = 600
    
    LOG_EVERY_N = 50
    
    files = glob.glob('/users/staff/leo/fastscratch/rise/1.0/exp03/*/feedbackmerged*')
    night_save = {}
    night_times = {}
    day_save = {}
    day_times = {}
    for i in stdplevs:
        night_save[i]=[]
        night_times[i]=[]
        day_save[i]=[]
        day_times[i]=[]
    t0 = time.time()
    for j in range(len(files)):
        if (j % LOG_EVERY_N) == 0: print(j)
        try:
            df = xr.open_dataset(files[j])
            for i in range(len(stdplevs)):
                df_new = df[dict(pressure=i)]
                if len(df_new.time) >= window:
                    df_snht = df_new[dict(hour=0)]
                    snht = test(np.array(df_snht.temperatures), window, missing)
                    night_save[stdplevs[i]].append(np.array(snht))
                    night_times[stdplevs[i]].append(np.array(df_snht.datum))
    
                    df_snht = df_new[dict(hour=1)]
                    snht = test(np.array(df_snht.temperatures), window, missing)
                    day_save[stdplevs[i]].append(np.array(snht))
                    day_times[stdplevs[i]].append(np.array(df_snht.datum))
        except:
            print('xxxxxxxxx', files[j], j)
            break
    print(time.time()-t0)
    
    
    # In[41]:
    
    
    stdplevs = [10.0, 20.0, 30.0, 50.0, 70.0, 100.0, 150.0, 200.0, 250.0, 300.0, 400.0, 500.0, 700.0, 850.0, 925.0, 1000.0]
    
    nlen = 0
    max_time = 0
    for i in day_times[500]:
        if len(i) > nlen:
            max_time = i
    
    for j in stdplevs:
        mean_snht = []
        print(j)
    #     counter = 0
    #     print(len(max_time))
        ic = []
        for i in max_time:
    #         if (counter % 100) == 0: print(counter)
    #         counter += 1
            summerizer = []
            inputcounter = 0
            for h in range(len(day_save[j][:100])):
                if i in day_times[j][h]:
                    summerizer.append(day_save[j][h][day_times[j][h] == i][0])
                    inputcounter += 1
            if divider != 0:
                mean_snht.append(np.nanmean(summerizer))
            else:
                mean_snht.append(0)
            ic.append(inputcounter)
                    
        fig, ax = maplt.subplots(1, figsize = (15,5))
        ax.plot(np.array(max_time),np.array(mean_snht),color = 'blue', alpha = 0.6, label='SNHT - ' + str(j) + ' hPa', )
        ax.plot(np.array(max_time),np.array(ic),color = 'green', alpha = 0.6, label='number of input stations', )
        ax.set_ylabel('SNHT')
        ax.set_xlabel('time')
        ax.legend(loc='upper right')
        ax.grid()
        maplt.show()
        maplt.close()
    
    
    # In[10]:
    
    
    len(day_save[10])
    
    
