#!/usr/bin/env python
# coding: utf-8

import warnings
warnings.filterwarnings("ignore")#, category=DeprecationWarning) 
import numpy
import numpy as np
from numba import config,version_info
config.THREADING_LAYER = 'tbb'

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
from anomaly import ndnanmean2,nnanmean,danomaly
# import ftest
#from multiprocessing import Pool
#import odb
#from eccodes import *
from functools import partial
#from collections import OrderedDict
#import subprocess
import json
#import gzip
# from retrieve_fb_jra55 import add_feedback
import copy
import pickle
#import xarray as xr
import ruptures as rpt
from pympler.asizeof import asizeof
import ray

#plt.rcParams['lines.linewidth'] = 3

#RAOBCORE constants
try:
    
    with open('RCx.json') as f:
        json.load(f,RC)
except:
    
    RC=dict(
        snht_maxlen = 1460,
        snht_min_sampsize = 80,
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
    RC['mon_thresh'] = 0
    RC['mean_maxlen'] = RC['snht_maxlen']*2
    RC['plot']=False
    RC['CPUs']=60
    RC['min_neighbours']=[[3,10],[10,30]]
    RC['weight_distance']=[3000,5000]
#    RC['ri_min_sampsize']=[330,30]
    RC['ri_min_sampsize']=[330,80]
    RC['initial_adjustments']='era5neighbours' #'era5' #'era5bc' #'era5neighbours'
    #RC['transdict']={'montemp':'mtemperatures','rasocorrmon':'madjustments','goodmon':'goodmon'}
    RC['transdict']={'montemp':'mini_adjusted_temperatures','rasocorrmon':'madjustments','goodmon':'goodmon'}
    RC['richfuture']=True 
    RC['findfuture']=True
    RC['goodsondes']=[141,142,123,124,125,113,114,79,80,81,70,152,177,183] # from WMO manual on codes, 2019 edition
    RC['apriori_prob_default']=0.01
    
    with open('RC.json','w') as f:
        json.dump(RC,f)
for k in 'years','months','days','stdplevs','pidx':
    if type(k) is int:
        RC[k]=np.array(RC[k],dtype=np.int32)
    else:
        RC[k]=np.array(RC[k])
    if k in  ('months','days'):
        RC[k]-=1

#RC['plot']=True
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


def ifray(dec, condition):
    def decorator(func):
        if not condition:
            # Return the function unchanged, not decorated.
            return func
        return dec(func)
    return decorator

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

    
def read_hadCRUT5(path,prefix,refdate):

    t1=time.time()
    fn=path+prefix+".anomalies.ensemble_mean.nc"
    fvar='tas_mean'
    fevar='tas'
    try:
        with netCDF4.Dataset(fn,"r") as f:
            f.set_auto_mask(False)
            hadmed=f.variables[fvar][:]
            hadmed[np.abs(hadmed)>1.e29]=np.nan
            sdate=f.variables['time'].getncattr('units')
            index=refdate[0]-int(sdate.split('-')[0].split()[-1])*12
    except:
        print((fn+' not found'))
    print('hadmed',time.time()-t1)

    return hadmed[index:,:,:] #,hadtem,hadens

#@njit
def do_test(fg_dep,temperature,rc_month,snht_maxlen,snht_min_sampsize,miss_val):

    tmean = np.zeros((fg_dep.shape[2]),dtype=fg_dep.dtype)
    tsquare = np.zeros((fg_dep.shape[2]),dtype=fg_dep.dtype)
    count = np.zeros(shape=(fg_dep.shape[2]),dtype=np.int32)
    tsas=np.zeros(fg_dep.shape,dtype=fg_dep.dtype)
    tsarad=np.zeros((1,6,fg_dep.shape[2]),dtype=fg_dep.dtype)
    
    parr=np.array((snht_maxlen,snht_min_sampsize))
    for ih in range(fg_dep.shape[0]):
        for ip in range(fg_dep.shape[1]):
            #SNHT.numba_snhteqmov(fg_dep[ih,ip,:].flatten(), rc_month,tsas[ih,ip,:], 
                                                 #snht_maxlen,snht_maxmiss,miss_val,
                                                 #count,tmean,tsquare)
            SNHT.numba_snhtmov_njit(fg_dep[ih,ip,:].flatten(), tsas[ih,ip,:], 
                                                 parr,
                                                 count,tmean,tsquare)
            #print(ih,ip)
            
    for ip in range(6):
        #SNHT.numba_snhteqmov((temperature[1,ip,:]-temperature[0,ip,:]).flatten(), rc_month, tsarad[0,ip,:], 
                             #snht_maxlen,snht_maxmiss,miss_val,count,tmean,tsquare)
        SNHT.numba_snhtmov_njit((temperature[1,ip,:]-temperature[0,ip,:]), tsarad[0,ip,:], 
                                             parr,
                                             count,tmean,tsquare)
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
                                               RC['snht_maxlen'],RC['snht_min_sampsize'],RC['miss_val'],
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
    
    totalo=np.concatenate((finfo['tsarad'][0,:,:],finfo['tsas'][0,:,:],finfo['tsas'][1,:,:]),axis=0)
    weight=np.array([2.0]*6+[1.0]*6+[0.8]*3+[0.5]*5+[1.0]*6+[0.8]*3+[0.5]*5)
    weight/=np.mean(weight)
    total=np.zeros(totalo.shape[1])
    for i in range(len(weight)):
        total+=totalo[i,:]*weight[i]
    print('totalo',totalo.shape)
    if RC['plot']:
        plt.figure(figsize=(6,12))
        lmax=6#totalo.shape[0]
        for l in range(lmax):    
            plt.subplot(lmax+1,1,l+1)
            plt.plot(finfo['days'][::5]/365.25,totalo[l][::5],label=str(l))
            plt.legend()

        plt.subplot(lmax+1,1,lmax+1)
        plt.plot(finfo['days'][::5]/365.25,total[::5],label='total')
        plt.legend()
        plt.title(finfo['sid'])
        plt.tight_layout()
    
    totalo=np.zeros_like(total)
    totalo[:]=total[:]
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
        with open(finfo['wigos']+'_breakanalysis.json','w') as f:
            json.dump(dict(days_since_1900=[int(d) for d in finfo['days']],
                           tsasum=[float(d) for d in totalo],
                           breaks=[int(d) for d in breaklist[:l]]),f)
        #with open(finfo['sid']+'_breakanalysis.json') as f:
            #d=json.load(f)
    
    return breaklist[:l]
    

#@njit
def metadata_analysis(sonde_type):
    schroedlist=np.array(['VN','Vn'],dtype='S4')
    apriori_probs=(np.empty_like(sonde_type,dtype=np.float32)).fill(RC['apriori_prob_default'])
    rcg=np.concatenate((np.array(RC['goodsondes'],dtype='S4'),schroedlist))
    
    bmatch=0
    gmatch=0
    count=0
    lbound=np.max([sonde_type.shape[0]-RC['mean_maxlen'],-1])
    i=sonde_type.shape[0]-1
    while i>lbound:
        found=False
        for r in rcg:            
            while r==sonde_type[i,0][:len(r)] and i>lbound:
                gmatch+=1
                i-=1
                found=True
        #if found:
            #break
        if not found and sonde_type[i,0] not in [b'nan',b'nnnn']:
                bmatch+=1
        count +=1   
        i-=1
    if gmatch>bmatch and gmatch >= RC['snht_min_sampsize']:
        lastsonde=1
    else:
        lastsonde=0
    
    return apriori_probs,lastsonde
#@njit
def bilin(hadmed,temperatures,rcy,rcm,days,lat,lon):
    ilat=int(np.floor((lat+90)/5))
    ilon=int(np.floor((180+lon)/5))
    if ilon>=hadmed.shape[2]:
        ilon-=hadmed.shape[2]
    t850=np.zeros(len(days))
    t850.fill(np.nan)
    tm=0
    mold=(rcy[days[0]]-1900)*12+rcm[days[0]]
    mc=0
    iref=0
    hadmedsmooth=rmeanw(hadmed[(rcy[days[:]]-1900)*12+rcm[days[:]],ilat,ilon],30)
    for i in range(len(days)):
        #mc+=1
        #if temperatures[0,13,i]==temperatures[0,13,i] and temperatures[1,13,i]==temperatures[1,13,i]:
            #tm+=temperatures[0,13,i]+temperatures[1,13,i]
        #elif temperatures[0,13,i]==temperatures[0,13,i]:
            #tm+=temperatures[0,13,i]
        #elif temperatures[0,13,i]==temperatures[0,13,i]:
            #tm+=temperatures[1,13,i]
        #else:
            #mc-=1
            
        #if (rcy[days[i]]-1900)*12+rcm[days[i]] != mold:
            #if mc>2: 
                #t850[iref:i]=tm/mc-hadmed[(rcy[days[iref]]-1900)*12+rcm[days[iref]],ilat,ilon]
            #else:
                #t850[iref:i]=np.nan
                  
            #mc=0
            #tm=0.
            #mold=(rcy[days[i]]-1900)*12+rcm[days[i]]
            #iref=i
        
        #if mc>2:
            
            #t850[iref:i]=tm/mc-hadmed[(rcy[days[iref]]-1900)*12+rcm[days[iref]],ilat,ilon]
        #else:
            #t850[iref:i]=np.nan
            
        t850[i]=np.nanmean(temperatures[:,13,i])-hadmedsmooth[i] #rmeanw(hadmed[(rcy[days[i]]-1900)*12+rcm[days[i]],ilat,ilon],30))
        
    return t850
    
    
#@ifray(ray.remote,RC['findfuture'])
def RAOB_findbreaks(method,fn):

    #import SNHT
    sys.path.append(os.path.expanduser('~/python/Rasotools/rasotools/'))
    sys.path.append(os.path.expanduser('~/python/Rasotools/'))
    from utils import rmeanw
    tt=time.time()
    finfo={}

    try:
        
        with h5py.File(fn,'r') as f:
            if f['datum'].shape[0]<RC['snht_maxlen']:
                return finfo
            finfo['days']=f['datum'][:]
            #print(fn.split('/')[-1],finfo['days'].shape)
            finfo['fg_dep']=-f['era5_fgdep'][:,np.array(RC['pidx']),:]
            finfo['era5bc']=-f['bias_estimate'][:,np.array(RC['pidx']),:] # ERA5 bias correction - to be used for adjusting most recent part of some series
            finfo['temperatures']=f['temperatures'][:,np.array(RC['pidx']),:]
            #finfo['fg']=finfo['temperatures']-finfo['fg_dep']
            finfo['lon']=f['lon'][0]
            finfo['lat']=f['lat'][0]
            finfo['sid']=fn.split('/')[-1][:-3]
            finfo['wigos']=f.attrs['unique_source_identifier'].decode("utf-8")
            if 'sonde_type' in f.keys():       
                #print(f['sonde_type'].shape)
                tt=time.time()
                finfo['apriori_probs'],finfo['lastsonde']=metadata_analysis(f['sonde_type'][0,0,:].view('S{}'.format(f['sonde_type'].shape[3])))
                #print(time.time()-tt)
            else:
                finfo['apriori_probs']=np.zeros_like(finfo['days'],dtype=np.float32)
                finfo['lastsonde']=0
        
        bins=RC['months'][finfo['days']]
        finfo['anomalies']=np.empty_like(finfo['temperatures'])
        finfo['anomalies'].fill(np.nan)
        s=finfo['temperatures'].shape
        finfo['climatologies']=np.zeros_like(finfo['temperatures'],shape=(s[0],s[1],12))
        ccount=np.zeros(12)
        good=danomaly(finfo['temperatures'],bins,ccount,finfo['anomalies'],finfo['climatologies'])
        finfo['hadcrut']=bilin(hadmed,finfo['anomalies'],RC['years'],RC['months'],finfo['days'],finfo['lat'],finfo['lon'])
        
        #print('before snht: {:6.4f}'.format(time.time()-tt)) 
        #print(typeof(finfo['fg_dep']),typeof(RC['months'][finfo['days']]),
              #typeof(RC['snht_maxlen']),typeof(RC['snht_maxmiss']),typeof(RC['miss_val']))
        if method=='SNHT':
            
            finfo['tsas'],finfo['tsarad']=do_test(finfo['fg_dep'],finfo['temperatures'],RC['months'][finfo['days']],
                                                  RC['snht_maxlen'],RC['snht_min_sampsize'],RC['miss_val'])  
            
            #print('before analyze: {:6.4f}'.format(time.time()-tt))  
            #print(version_info)
            finfo['breaklist']=np.sort(break_analyze(finfo))
        elif method=='Binseg':
            testset=np.concatenate((finfo['fg_dep'][0,2:8,:],finfo['fg_dep'][1,2:8,:],finfo['fg_dep'][1,2:8,:]-finfo['fg_dep'][0,2:8,:],finfo['hadcrut'].reshape(1,finfo['hadcrut'].shape[0])),axis=0)
            algo = rpt.Binseg(model="l2").fit(testset.T)
            finfo['breaklist'] = algo.predict(pen=1.5*np.log(testset.shape[1]) * testset.shape[0] * np.nanvar(testset))
            finfo['breaklist'] = np.sort(finfo['breaklist'])[:-1]
            del algo
            del testset
        
        
        print(finfo['sid']+' findbreaks: {:6.4f}'.format(time.time()-tt))  
        obj_ref= ray.put(finfo)
    except MemoryError as e:
        print(e)
        return None
    
    return obj_ref

ray_RAOB_findbreaks=ray.remote(RAOB_findbreaks)

def RAOB_adjustbs(finfo):
    
    sys.path.append(os.path.expanduser('~/python/Rasotools/rasotools/'))
    sys.path.append(os.path.expanduser('~/python/Rasotools/'))
    from utils import tdist,extract,rmeanw
    from anomaly import ndnanmean2,nnanmean,ndnanvar2
    tt=time.time()
    if type(finfo) is not dict:
        finfo=ray.get(finfo)
    finfo['adjusted_fg_dep']=np.empty_like(finfo['fg_dep'])
    finfo['adjusted_fg_dep'][:]=finfo['fg_dep'][:]
    mask=np.isnan(finfo['adjusted_fg_dep'])
    s=finfo['adjusted_fg_dep'].shape
    tmean = np.zeros_like(finfo['adjusted_fg_dep'],shape=(s[2],12))
    tsquare = np.zeros_like(finfo['adjusted_fg_dep'],shape=(s[2],12))
    count = np.zeros(shape=(s[2],12),dtype=np.int32)
    tsa=np.zeros_like(finfo['adjusted_fg_dep'],shape=(s[2]))

    fb=finfo['breaklist']
    finfo['adjustments']=np.zeros_like(finfo['adjusted_fg_dep'],shape=(s[0],s[1],fb.shape[0]))
    finfo['adjusted_temperatures']=np.empty_like(finfo['adjusted_fg_dep'])
    finfo['adjusted_temperatures'][:]=finfo['temperatures'][:]
    finfo['adjusted_anomalies']=np.empty_like(finfo['adjusted_fg_dep'])
    finfo['adjusted_anomalies'][:]=finfo['anomalies'][:]
    
    sh=finfo['adjustments'].shape
    RC_months=np.array(RC['months'])[finfo['days']]-1
    ipl=2
    if RC['plot']:
        plt.subplot(2,1,1)
        plt.plot(finfo['days']/365.25,rmeanw(finfo['fg_dep'][0,ipl,:],30))
        plt.subplot(2,1,2)
        plt.plot(finfo['days']/365.25,rmeanw(finfo['fg_dep'][1,ipl,:],30))
    
    break_profiles=np.zeros_like(finfo['adjusted_fg_dep'],shape=(len(fb),sh[0],sh[1]))
    break_confidence=np.zeros_like(finfo['adjusted_fg_dep'],shape=(len(fb),sh[0],sh[1]))

    goodbreaklist=[]

    for ib in range(len(fb)-1,-1,-1):
        if ib>0:
            if fb[ib]-fb[ib-1]>RC['snht_min_sampsize']+30:
                safety=0 #60
            else:
                safety=0
            istart=np.max((fb[ib]-RC['mean_maxlen']+safety,fb[ib-1]+safety))
        else:
            istart=np.max((fb[ib]-RC['mean_maxlen']+0,0))
        istop=np.min((fb[ib]+RC['mean_maxlen'],s[2]))

        sum1=np.sum(~np.isnan(finfo['fg_dep'][:,11,fb[ib]:istop]),axis=1)
        ri_min_sampsize=RC['ri_min_sampsize'][1]
        if np.any(sum1<3*ri_min_sampsize):                           
            sum2=np.sum(~np.isnan(finfo['fg_dep'][:,11,fb[ib]:]),axis=1)
            if any(sum2<3*ri_min_sampsize):
                pass
            else:
                x=np.cumsum(~np.isnan(finfo['fg_dep'][:,11,fb[ib]:]),axis=1)
                y00=np.searchsorted(x[0,:],4*ri_min_sampsize)
                y12=np.searchsorted(x[1,:],4*ri_min_sampsize)
                mml=np.max([y00,y12])-1
                    
                istop=fb[ib]+mml # probably a long data gap
                
        else:
            pass
        
        ini=np.zeros(finfo['adjusted_fg_dep'].shape[:2])
        ini2=np.zeros(finfo['adjusted_fg_dep'].shape[:2])
        delta=0
        break_profiles[ib,:,:]=-ndnanmean2(finfo['adjusted_fg_dep'][:,:,istart:fb[ib]-delta],ini,ri_min_sampsize)+\
                               ndnanmean2(finfo['adjusted_fg_dep'][:,:,fb[ib]+delta:istop],ini2,3*ri_min_sampsize)
        ini=np.zeros(finfo['adjusted_fg_dep'].shape[:2])
        ini2=np.zeros(finfo['adjusted_fg_dep'].shape[:2])
        break_confidence[ib,:,:]=np.sqrt(0.5*(ndnanvar2(finfo['adjusted_fg_dep'][:,:,istart:fb[ib]-delta],ini,ri_min_sampsize)+\
                               ndnanvar2(finfo['adjusted_fg_dep'][:,:,fb[ib]+delta:istop],ini2,3*ri_min_sampsize)))

        sig=0
        nonsig=0
        #print('break_confidence',break_confidence)
        for ih in range(sh[0]):
            for ip in range(sh[1]):
                if break_profiles[ib,ih,ip]==break_profiles[ib,ih,ip]:
                    if np.abs(break_profiles[ib,ih,ip])>2*1.96*break_confidence[ib,ih,ip]: # factor 2 corresponds to an autocorrelation of ~0.3
                        sig+=1
                    else:
                        nonsig+=1
        if sig<4:
            break_profiles[ib,:,:]=np.nan
        else:
            goodbreaklist.insert(0,ib)
    
        for ih in range(sh[0]):
            for ip in range(sh[1]):
                    
                if ~np.isnan(break_profiles[ib,ih,ip]):
                    finfo['adjustments'][ih,ip,:ib]-=break_profiles[ib,ih,ip]
                    finfo['adjusted_fg_dep'][ih,ip,:fb[ib]]+=break_profiles[ib,ih,ip]
                    finfo['adjusted_temperatures'][ih,ip,:fb[ib]]+=break_profiles[ib,ih,ip]
                    finfo['adjusted_anomalies'][ih,ip,:fb[ib]]+=break_profiles[ib,ih,ip]
        
        
        #print(ib,fb[ib],break_profiles[ib,:,:])
    finfo['goodbreaklist']=list(fb[np.array(goodbreaklist,dtype=int)]) 
    
    
    
    print('adjust: {:6.4f}'.format(time.time()-tt))  
    finfo['adjusted_temperatures'][mask]=np.nan
    if RC['plot']:
        plt.subplot(2,1,1)
        plt.plot(finfo['days']/365.25,rmeanw(finfo['adjusted_fg_dep'][0,ipl,:],30))
        plt.plot(finfo['days']/365.25,finfo['adjusted_temperatures'][0,ipl,:]-finfo['temperatures'][0,ipl,:])
        plt.plot(finfo['days'][::finfo['days'].shape[0]-1]/365.25,[0,0],'k')
        plt.title(finfo['sid'][-6:])
        plt.subplot(2,1,2)
        plt.plot(finfo['days']/365.25,rmeanw(finfo['adjusted_fg_dep'][1,ipl,:],30))
        plt.plot(finfo['days']/365.25,finfo['adjusted_temperatures'][1,ipl,:]-finfo['temperatures'][1,ipl,:])
        plt.plot(finfo['days'][::finfo['days'].shape[0]-1]/365.25,[0,0],'k')
        plt.title(finfo['sid'][-6:])
        plt.tight_layout()
        plt.savefig(finfo['sid']+'_adjustments.png')
        #plt.show()
        plt.close()

    return {'goodbreaklist':np.array(finfo['goodbreaklist']),'break_profiles':break_profiles[goodbreaklist,:,:],
                    'adjustments':finfo['adjustments'][:,:,goodbreaklist],
                    'adjusted_temperatures':finfo['adjusted_temperatures'],
                    'adjusted_anomalies':finfo['adjusted_anomalies'],
                    'adjusted_fg_dep':finfo['adjusted_fg_dep'],}
        
        

def RAOB_adjust(finfo):
    
    sys.path.append(os.path.expanduser('~/python/Rasotools/rasotools/'))
    sys.path.append(os.path.expanduser('~/python/Rasotools/'))
    from utils import tdist,extract,rmeanw
    from anomaly import ndnanmean2,nnanmean,ndnanvar2
    tt=time.time()
    if type(finfo) is not dict:
        finfo=ray.get(finfo)
    finfo['adjusted_fg_dep']=np.empty_like(finfo['fg_dep'])
    finfo['adjusted_fg_dep'][:]=finfo['fg_dep'][:]
    mask=np.isnan(finfo['adjusted_fg_dep'])
    s=finfo['adjusted_fg_dep'].shape
    tmean = np.zeros_like(finfo['adjusted_fg_dep'],shape=(s[2],12))
    tsquare = np.zeros_like(finfo['adjusted_fg_dep'],shape=(s[2],12))
    count = np.zeros(shape=(s[2],12),dtype=np.int32)
    tsa=np.zeros_like(finfo['adjusted_fg_dep'],shape=(s[2]))

    fb=finfo['breaklist']
    finfo['adjustments']=np.zeros_like(finfo['adjusted_fg_dep'],shape=(s[0],s[1],fb.shape[0]))
    finfo['adjusted_temperatures']=np.empty_like(finfo['adjusted_fg_dep'])
    finfo['adjusted_temperatures'][:]=finfo['temperatures'][:]
    finfo['adjusted_anomalies']=np.empty_like(finfo['adjusted_fg_dep'])
    finfo['adjusted_anomalies'][:]=finfo['anomalies'][:]
    
    sh=finfo['adjustments'].shape
    RC_months=np.array(RC['months'])[finfo['days']]-1
    ipl=2
    if RC['plot']:
        plt.subplot(2,1,1)
        plt.plot(finfo['days']/365.25,rmeanw(finfo['fg_dep'][0,ipl,:],30))
        plt.subplot(2,1,2)
        plt.plot(finfo['days']/365.25,rmeanw(finfo['fg_dep'][1,ipl,:],30))
    
    break_profiles=np.zeros_like(finfo['adjusted_fg_dep'],shape=(len(fb),sh[0],sh[1]))
    break_confidence=np.zeros_like(finfo['adjusted_fg_dep'],shape=(len(fb),sh[0],sh[1]))

    goodbreaklist=[]

    for ib in range(len(fb)-1,-1,-1):
        if ib>0:
            if fb[ib]-fb[ib-1]>RC['snht_min_sampsize']+30:
                safety=60
            else:
                safety=0
            istart=np.max((fb[ib]-RC['mean_maxlen']+safety,fb[ib-1]+safety))
        else:
            istart=np.max((fb[ib]-RC['mean_maxlen']+0,0))
        istop=np.min((fb[ib]+RC['mean_maxlen'],s[2]))
        count=np.sum(~np.isnan(finfo['fg_dep'][:,11,fb[ib]:istop]),axis=1)
        stops=np.zeros(2,dtype=np.int32)+istop
        rcss=np.int(3*RC['snht_min_sampsize'])
        for j in range(2):
            if count[j]<rcss and np.sum(~np.isnan(finfo['fg_dep'][j,11,fb[ib]:]))>rcss:
                i=istop
                while i<s[2] and count[j]<rcss:
                    if finfo['fg_dep'][j,11,i]==finfo['fg_dep'][j,11,i]:
                        count[j]+=1
                    i+=1
                stops[j]=i
        if np.any(stops!=istop):
            print(finfo['sid'][-6:],ib,istop,stops)
        istop=np.max(stops)
        ini=np.zeros(finfo['adjusted_fg_dep'].shape[:2])
        ini2=np.zeros(finfo['adjusted_fg_dep'].shape[:2])
        break_profiles[ib,:,:]=-ndnanmean2(finfo['adjusted_fg_dep'][:,:,istart:fb[ib]-30],ini,RC['snht_min_sampsize'])+\
                               ndnanmean2(finfo['adjusted_fg_dep'][:,:,fb[ib]+30:istop],ini2,RC['snht_min_sampsize'])
        ini=np.zeros(finfo['adjusted_fg_dep'].shape[:2])
        ini2=np.zeros(finfo['adjusted_fg_dep'].shape[:2])
        break_confidence[ib,:,:]=np.sqrt(0.5*(ndnanvar2(finfo['adjusted_fg_dep'][:,:,istart:fb[ib]-30],ini,RC['snht_min_sampsize'])+\
                               ndnanvar2(finfo['adjusted_fg_dep'][:,:,fb[ib]+30:istop],ini2,RC['snht_min_sampsize'])))

        sig=0
        nonsig=0
        #print('break_confidence',break_confidence)
        for ih in range(sh[0]):
            for ip in range(sh[1]):
                if break_profiles[ib,ih,ip]==break_profiles[ib,ih,ip]:
                    if np.abs(break_profiles[ib,ih,ip])>2*1.96*break_confidence[ib,ih,ip]: # factor 2 corresponds to an autocorrelation of ~0.3
                        sig+=1
                    else:
                        nonsig+=1
        if sig<4:
            break_profiles[ib,:,:]=np.nan
        else:
            goodbreaklist.insert(0,ib)
    
        for ih in range(sh[0]):
            for ip in range(sh[1]):
                #break_profiles[ib,ih,ip]=SNHT.numba_meaneqmov(fg_dep[ih,ip,istart:istop], 
                                                          #RC_months[istart:istop], tsa, 
                                                          #RC['mean_maxlen'],RC['snht_maxmiss'],RC['miss_val'], 
                                                          #count, tmean, tsquare, kref=fb[ib]-istart)
                
                
                #if np.sum(~np.isnan(fg_dep[ih,ip,istart:fb[ib]]))>RC['mean_maxlen']//4-RC['snht_maxmiss'] and np.sum(~np.isnan(fg_dep[ih,ip,fb[ib]:istop]))>RC['mean_maxlen']//4-RC['snht_maxmiss']:
                    
                    #break_profiles[ib,ih,ip]=-np.nanmean(fg_dep[ih,ip,istart:fb[ib]])+np.nanmean(fg_dep[ih,ip,fb[ib]:istop])
                #else:
                    #break_profiles[ib,ih,ip]=np.nan
                    
                if ~np.isnan(break_profiles[ib,ih,ip]):
                    finfo['adjustments'][ih,ip,:ib]-=break_profiles[ib,ih,ip]
                    finfo['adjusted_fg_dep'][ih,ip,:fb[ib]]+=break_profiles[ib,ih,ip]
                    finfo['adjusted_temperatures'][ih,ip,:fb[ib]]+=break_profiles[ib,ih,ip]
                    finfo['adjusted_anomalies'][ih,ip,:fb[ib]]+=break_profiles[ib,ih,ip]
        
        
        #print(ib,fb[ib],break_profiles[ib,:,:])
    finfo['goodbreaklist']=list(fb[np.array(goodbreaklist,dtype=int)]) 
    
    
    
    print('adjust: {:6.4f}'.format(time.time()-tt))  
    finfo['adjusted_temperatures'][mask]=np.nan
    if RC['plot']:
        plt.subplot(2,1,1)
        plt.plot(finfo['days']/365.25,rmeanw(finfo['adjusted_fg_dep'][0,ipl,:],30))
        plt.plot(finfo['days']/365.25,finfo['adjusted_temperatures'][0,ipl,:]-finfo['temperatures'][0,ipl,:])
        plt.plot(finfo['days'][::finfo['days'].shape[0]-1]/365.25,[0,0],'k')
        plt.title(finfo['sid'][-6:])
        plt.subplot(2,1,2)
        plt.plot(finfo['days']/365.25,rmeanw(finfo['adjusted_fg_dep'][1,ipl,:],30))
        plt.plot(finfo['days']/365.25,finfo['adjusted_temperatures'][1,ipl,:]-finfo['temperatures'][1,ipl,:])
        plt.plot(finfo['days'][::finfo['days'].shape[0]-1]/365.25,[0,0],'k')
        plt.title(finfo['sid'][-6:])
        plt.tight_layout()
        plt.savefig(finfo['sid']+'_adjustments.png')
        #plt.show()
        plt.close()

    return {'goodbreaklist':np.array(finfo['goodbreaklist']),'break_profiles':break_profiles[goodbreaklist,:,:],
                    'adjustments':finfo['adjustments'][:,:,goodbreaklist],
                    'adjusted_temperatures':finfo['adjusted_temperatures'],
                    'adjusted_anomalies':finfo['adjusted_anomalies'],
                    'adjusted_fg_dep':finfo['adjusted_fg_dep'],}
    
#@ifray(ray.remote,RC['findfuture'])    
def RAOB_adjustbreaks(method,dists,finfo,i):
    
    adjustments=copy.copy(finfo)
    adjustments['sdists']=extract(dists,i)

    adjustments.update(RAOB_adjust(finfo))
    
    return ray.put(adjustments)

ray_RAOB_adjustbreaks=ray.remote(RAOB_adjustbreaks)

    
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
def goodmon(var,idx,out):

    l=0
    for i in range(len(idx)-1):
        if idx[i+1]>idx[i]:
            l+=1
    if l==0:
        print('idx not valid')
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
def monmean(var,idx,out,thresh=0,goodmon=goodmon):
    l=0
    for i in range(len(idx)-1):
        if idx[i+1]>idx[i]:
            l+=1
    if l==0:
        print('idx not valid')
    out[:]=np.nan
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
                
    return out[:,:,:l]

#@njit
def do_monmean(vars,tidx,fi):
    
    idx=np.searchsorted(fi['days'],tidx)
    mdays=idx[1:]-idx[:-1]
    midx=np.where(mdays>0)[0]
    
    fi['months']=tidx[midx]
            
    for v in vars:
        out=np.empty((fi[v].shape[0],fi[v].shape[1],fi['months'].shape[0]),dtype=fi[v].dtype)
        fi['m'+v]=monmean(fi[v],idx,out,RC['mon_thresh'])
        
    out=np.empty((fi[v].shape[0],fi[v].shape[1],fi['months'].shape[0]),dtype=np.int32)
    fi['goodmon']=goodmon(fi['temperatures'],idx,out)
            
                
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
                elif name in ['montemp','goodmon','rasocorrmon']: #,'goodmon','rasocorrmon','eracorrmon']:
                    x = dst.createVariable(name, variable.datatype, variable.dimensions)
                    if variable.datatype in [np.dtype('float32'),np.dtype('float64')]:
                        
                        dst[name][:]=np.nan
                    else:
                        dst[name][:]=0
                    
                    if name in ['montemp','goodmon']:
                        
                        dst[name][:,:RC['pidx'].shape[0],:] = fi[RC['transdict'][name]][:]
                    elif name=='rasocorrmon':
                        dst[name][:,:RC['pidx'].shape[0],:] = fi['mtemperatures']-fi['mini_adjusted_temperatures']

                       
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
def write_bgmonmean(ipath,opath,prefix,fi):
    
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
                elif name in ['goodmon','montemp']: #,'goodmon','rasocorrmon','eracorrmon']:
                    x = dst.createVariable(name, variable.datatype, variable.dimensions)
                    if variable.datatype in [np.dtype('float32'),np.dtype('float64')]:
                        
                        dst[name][:]=np.nan
                    else:
                        dst[name][:]=0
                    
                    if name in ['goodmon']:
                        
                        dst[name][:,:RC['pidx'].shape[0],:] = fi[RC['transdict'][name]][:]
                    elif name=='montemp':
                        dst[name][:,:RC['pidx'].shape[0],:] = fi['mfg']
                    else:
                        pass
                       
                    # copy variable attributes all at once via dictionary
                    dst[name].setncatts(src[name].__dict__)
                    
                elif name in ['eracorrmon','rasocorrmon']:
                    pass
                else:
                    x = dst.createVariable(name, variable.datatype, variable.dimensions)
                    dst[name][:] = src[name][:]
                    # copy variable attributes all at once via dictionary
                    dst[name].setncatts(src[name].__dict__)
    except FileNotFoundError:
        print(fni,'could not be copied to',fno)
    
    return

#@njit
def addini(adjorig,iniad):
    
    ad=adjorig.shape
    res_with_ini=np.concatenate((adjorig[:,:,0:1],adjorig[:],np.zeros((ad[0],ad[1],1))),axis=2)
    for ih in range(ad[0]):
        for ip in range(ad[1]):
            res_with_ini[ih,ip,:-1]+=iniad[ih,ip]
            
    return res_with_ini

def write_adjustment(ipath,opath,prefix,fi,rich_ref0=None,rich_ref1=None,initial_adjust_RAOB=None,initial_adjust_RICH=None):
    
    if len(fi['goodbreaklist'])==0:
        print(fi['sid']+': no breaks for this station')
        return
    if 'ri' in prefix and not rich_ref0 and not rich_ref1:
        return
    
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
                    dst.createDimension(name, (len(fi['goodbreaklist'])+2 if not dimension.isunlimited() else None))
                else:
                    dst.createDimension(
                        name, (len(dimension) if not dimension.isunlimited() else None))
            # copy all file data except for the excluded
            for name, variable in src.variables.items():
                if name == 'datum':
                    x = dst.createVariable(name, variable.datatype, variable.dimensions)
                    dst[name][:] = np.concatenate((np.array([1]),fi['days'][fi['goodbreaklist'][:]],np.array([44998])),axis=0)
                    # copy variable attributes all at once via dictionary
                    dst[name].setncatts(src[name].__dict__)
                elif name in ['rasocorr','rasobreak']: #,'goodmon','rasocorrmon','eracorrmon']:
                    
                    if rich_ref0 and 'ri' in prefix:
                        suff=''
                        if rich_ref1:
                            suff='0'
                        x = dst.createVariable(name+suff, variable.datatype, variable.dimensions)
                        if variable.datatype in [np.dtype('float32'),np.dtype('float64')]:
                            
                            dst[name+suff][:]=np.nan
                        else:
                            dst[name+suff][:]=0
                        if 'rio' in prefix and rich_ref0:
                            ad=rich_ref0['obs']['rich_adjustments'].shape
                            dst[name+suff][:,:RC['pidx'].shape[0],:] = np.concatenate((rich_ref0['obs']['rich_adjustments'][:,:,0:1],
                                                                                  rich_ref0['obs']['rich_adjustments'][:],np.zeros((ad[0],ad[1],1))),axis=2)
                        elif 'rit' in prefix and rich_ref0:
                            ad=rich_ref0['tau']['rich_adjustments'].shape
                            dst[name+suff][:,:RC['pidx'].shape[0],:] = np.concatenate((rich_ref0['tau']['rich_adjustments'][:,:,0:1],
                                                                                  rich_ref0['tau']['rich_adjustments'][:],np.zeros((ad[0],ad[1],1))),axis=2)

                    if rich_ref1:
                        x = dst.createVariable(name, variable.datatype, variable.dimensions)
                        if variable.datatype in [np.dtype('float32'),np.dtype('float64')]:
                            
                            dst[name][:]=np.nan
                        else:
                                dst[name][:]=0
                        if 'rio' in prefix and rich_ref1:
                            ad=rich_ref1['obs']['rich_adjustments'].shape
                            #dst[name][:,:RC['pidx'].shape[0],:] = np.concatenate((rich_ref1['obs']['rich_adjustments'][:,:,0:1],
                            if initial_adjust_RICH:
                                if 'initial_adjustments' in initial_adjust_RICH.keys():
                                    
                                    dst[name][:,:RC['pidx'].shape[0],:] = addini(rich_ref1['obs']['rich_adjustments'],initial_adjust_RICH['initial_adjustments'])
                                else:
                                    dst[name][:,:RC['pidx'].shape[0],:] = addini(rich_ref1['obs']['rich_adjustments'],np.zeros(fi['adjustments'].shape[:2]))
                            else:
                                dst[name][:,:RC['pidx'].shape[0],:] = addini(rich_ref1['obs']['rich_adjustments'],np.zeros(fi['adjustments'].shape[:2]))
                            #rich_ref1['obs']['rich_adjustments'][:],np.zeros((ad[0],ad[1],1))),axis=2)
                        elif 'rit' in prefix and rich_ref1:
                            ad=rich_ref1['tau']['rich_adjustments'].shape
                            #dst[name][:,:RC['pidx'].shape[0],:] = np.concatenate((rich_ref1['tau']['rich_adjustments'][:,:,0:1],
                                                                                  #rich_ref1['tau']['rich_adjustments'][:],np.zeros((ad[0],ad[1],1))),axis=2)
                            if initial_adjust_RICH:
                                if 'initial_adjustments' in initial_adjust_RICH.keys():
                                    dst[name][:,:RC['pidx'].shape[0],:] = addini(rich_ref1['tau']['rich_adjustments'],initial_adjust_RICH['initial_adjustments'])
                                else:
                                    dst[name][:,:RC['pidx'].shape[0],:] = addini(rich_ref1['tau']['rich_adjustments'],np.zeros(fi['adjustments'].shape[1:]))
                            else:
                                dst[name][:,:RC['pidx'].shape[0],:] = addini(rich_ref1['tau']['rich_adjustments'],np.zeros(fi['adjustments'].shape[1:]))
                    else:
                        x = dst.createVariable(name, variable.datatype, variable.dimensions)
                        if variable.datatype in [np.dtype('float32'),np.dtype('float64')]:
                            
                            dst[name][:]=np.nan
                        else:
                                dst[name][:]=0
                        
                        #print(fi['sid'][-6:],initial_adjust_RAOB.keys(),fi['adjustments'].shape)
                        if initial_adjust_RAOB:
                            if 'initial_adjustments' in initial_adjust_RAOB.keys():
                                
                                dst[name][:,:RC['pidx'].shape[0],:] = addini(fi['adjustments'],initial_adjust_RAOB['initial_adjustments'])
                            else:
                                dst[name][:,:RC['pidx'].shape[0],:] = addini(fi['adjustments'],np.zeros(fi['adjustments'].shape[:2]))
                        else:
                            dst[name][:,:RC['pidx'].shape[0],:] = addini(fi['adjustments'],np.zeros(fi['adjustments'].shape[:2]))
                            
                    # copy variable attributes all at once via dictionary
                    dst[name].setncatts(src[name].__dict__)
                    
                else:
                    x = dst.createVariable(name, variable.datatype, variable.dimensions)
                    dst[name][:] = src[name][:]
                    # copy variable attributes all at once via dictionary
                    dst[name].setncatts(src[name].__dict__)
    except FileNotFoundError:
        print(fni,'could not be copied to',fno)
    
    return

#@njit
def make_adjusted_series(orig,adjustments,breaklist):
    
    adjusted_series=orig.copy()
    
    for ib in range(len(breaklist)-1,-1,-1):
        for ih in range(adjusted_series.shape[0]):
            for ip in range(adjusted_series.shape[1]):  
                if ib>0:
                    adjusted_series[ih,ip,breaklist[ib-1]:breaklist[ib]]-=adjustments[ih,ip,ib]
                else:
                    adjusted_series[ih,ip,:breaklist[ib]]-=adjustments[ih,ip,ib]
    
    return adjusted_series

#@ifray(ray.remote,RC['richfuture'])
def save_monmean(tfile,l,fi,rich_ref0=None,rich_ref1=None,initial_adjust_RAOB=None,initial_adjust_rich=None):
    
    #return ray.put({'break_profiles':break_profiles,'adjusted_temperatures':finfo['adjusted_temperatures']})
    #fnames= fi['sid']
    if type(fi[l]) is not dict:
        fi=ray.get(fi)[l]
    else:
        fi=fi[l]
    #fi['adjustments']=-fi['break_profiles']
    fi['ini_adjusted_temperatures']=fi['adjusted_temperatures'].copy()
    fi['fg']=fi['temperatures']-fi['fg_dep']
    fi['obs']={}
    fi['tau']={}
    
    tidx=np.where(RC['days']==1)[0]
    exper='exp02'
    ipath=fi['sid'][-6:].join(os.path.dirname(tfile).split(files[0][-9:-3]))
    plist=['temperatures','ini_adjusted_temperatures','fg']
    if rich_ref0:
        if type(rich_ref0[l]) is not dict:
            rich_ref0=ray.get(rich_ref0)[l]
        else:
            rich_ref0=rich_ref0[l]
            
        #if 'rich_adjustments' in rich_ref0[24]['obs'].keys():
            ##fi['obs']['rich_adjusted_temperatures_0']=make_adjusted_series(fi['temperatures'],rich_ref0[24]['obs']['rich_adjustments'],rich_ref0['goodbreaklist'])
            #plist=plist+['rich_adjusted_temperatures_0']
    if rich_ref1:
        if type(rich_ref1[l]) is not dict:
            rich_ref1=ray.get(rich_ref1)[l]
        else:
            rich_ref1=rich_ref1[l]
            
    if initial_adjust_RAOB:
        if type(initial_adjust_RAOB[l]) is not dict:
            initial_adjust_RAOB=ray.get(initial_adjust_RAOB)[l]
        else:
            initial_adjust_RAOB=initial_adjust_RAOB[l]

        if 'initial_adjustments' in initial_adjust_RAOB.keys():
            sh=fi['adjusted_temperatures'].shape
            for ih in range(sh[0]):
                for ip in range(sh[1]):
                    fi['ini_adjusted_temperatures'][ih,ip,:]-=initial_adjust_RAOB['initial_adjustments'][ih,ip]
    if initial_adjust_rich:
        if type(initial_adjust_rich[l]) is not dict:
            initial_adjust_rich=ray.get(initial_adjust_rich)[l]
        else:
            initial_adjust_rich=initial_adjust_rich[l]
        #if 'rich_adjustments' in rich_ref1[24]['obs'].keys():
            ##fi['obs']['rich_adjusted_temperatures_0']=make_adjusted_series(fi['temperatures'],rich_ref1[24]['obs']['rich_adjustments'],rich_ref1['goodbreaklist'])
            #plist=plist+['rich_adjusted_temperatures_1']
    do_monmean(plist,tidx,fi)
    prefix='feedbackglobbincorrmon'
    write_monmean(ipath,exper.join(ipath.split('exp03')),prefix,fi)
    prefix='feedbackglobbgmon'
    write_bgmonmean(ipath,exper.join(ipath.split('exp03')),prefix,fi) 
    prefix='feedbackglobbincorrsave'
    write_adjustment(ipath,exper.join(ipath.split('exp03')),prefix,fi,
                     initial_adjust_RAOB=initial_adjust_RAOB)#,initial_adjust_rich=initial_adjust_rich)
    if rich_ref0 or rich_ref1:
        for iens in [24]:#range(len(rich_ref0)):
            sens='{:02}'.format(iens)
            pfdict={'tau':'rit','obs':'rio'}
            for iens in [24]:#range(len(rich_ref0)):
                sh=fi['adjusted_temperatures'].shape
                if rich_ref0:
                    riref0=copy.copy(rich_ref0[iens])
                else:
                    riref0=rich_ref0
                riref1=copy.copy(rich_ref1[iens])
                for method in 'tau','obs':
                    prefix = 'feedbackglobbincorrsave_'+pfdict[method]+sens+'_' 
                    iar=None
                    if initial_adjust_rich:  
                        if len(initial_adjust_rich)>1:
                            
                            if 'initial_adjustments' in initial_adjust_rich[1][iens][method].keys():
                                iar=initial_adjust_rich[1][iens][method]

                    write_adjustment(ipath,exper.join(ipath.split('exp03')),prefix,fi,
                                     rich_ref0=riref0,
                                     rich_ref1=riref1,
                                     initial_adjust_RAOB=initial_adjust_RAOB,
                                     initial_adjust_RICH=iar)

ray_save_monmean=ray.remote(save_monmean)    
    
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
    
            
    tt=time.time()
    dists=np.zeros(len(finfo)*(len(finfo)+1)//2)
    lats=np.array([i['lat'] for i in finfo ])
    lons=np.array([i['lon'] for i in finfo ])
    print('latlon',time.time()-tt)
    dum=tdist(dists,lats.astype(np.float64),lons.astype(np.float64),0)
    print('tdist',time.time()-tt)
    #for i in range(len(finfo)):  
        #finfo[i]['stdists']=extract(dists,i,lats,lons)
    #print('extract',time.time()-tt)
    
    return ray.put(dists) #[finfo[i]['stdists'] for i in range(len(finfo))])

ray_stdists=ray.remote(stdists)


@njit(fastmath={'nsz','arcp','contract','afn','reassoc'},cache=True)
def lagmean(ldiff,sdiff,ldays,sdays,llmaxlen,slmaxlen,thresh,lprofl,sprofl,lprofr,sprofr):
    
    idx=np.searchsorted(sdays,ldays)
    for ih in range(ldiff.shape[0]):
        for ip in range(ldiff.shape[1]-1,-1,-1):
            cl=0
            for i in range(llmaxlen):
                j=idx[i]
                if j==sdays.shape[0]:
                    break
                if ldiff[ih,ip,i]==ldiff[ih,ip,i]:
                    if sdays[j]==ldays[i] and sdiff[ih,ip,j]==sdiff[ih,ip,j]:
                        lprofl[ih,ip]+=ldiff[ih,ip,i]
                        sprofl[ih,ip]+=sdiff[ih,ip,j]
                        cl+=1
            if cl>thresh:
                lprofl[ih,ip]/=cl
                sprofl[ih,ip]/=cl
            else:
                if ip<=6:
                    lprofl[ih,ip]=lprofl[ih,ip+1] # take value from next pressure level downward, assuming constant bias in the vertical
                    sprofl[ih,ip]=sprofl[ih,ip+1]
                else:    
                    lprofl[ih,ip]=np.nan
                    sprofl[ih,ip]=np.nan
                
            cr=0
            for i in range(llmaxlen,ldiff.shape[2]):
                j=idx[i]
                if j==sdays.shape[0]:
                    break
                if ldiff[ih,ip,i]==ldiff[ih,ip,i]:
                    if sdays[j]==ldays[i] and sdiff[ih,ip,j]==sdiff[ih,ip,j]:
                        lprofr[ih,ip]+=ldiff[ih,ip,i]
                        sprofr[ih,ip]+=sdiff[ih,ip,j]
                        cr+=1
            if cr>thresh:
                lprofr[ih,ip]/=cr
                sprofr[ih,ip]/=cr
            else:
                if ip<=6:
                    lprofr[ih,ip]=lprofr[ih,ip+1] # take value from next pressure level downward, assuming constant bias in the vertical
                    sprofr[ih,ip]=sprofr[ih,ip+1]
                else:
                    lprofr[ih,ip]=np.nan
                    sprofr[ih,ip]=np.nan
                
    return (lprofr-lprofl)-(sprofr-sprofl)
    

def rich_calcadjprofile(ldiff,sdiff,ldays,sdays,tbin,sbin,llmaxlen,slmaxlen):

    sprof=np.zeros(sdiff.shape[:2],dtype=ldiff.dtype)
    lprofl=np.zeros(ldiff.shape[:2],dtype=ldiff.dtype)
    sprofl=np.zeros(ldiff.shape[:2],dtype=ldiff.dtype)
    lprofr=np.zeros(ldiff.shape[:2],dtype=ldiff.dtype)
    sprofr=np.zeros(ldiff.shape[:2],dtype=ldiff.dtype)

    sprof=lagmean(ldiff,sdiff,ldays,sdays,llmaxlen,slmaxlen,RC['snht_min_sampsize'],lprofl,sprofl,lprofr,sprofr)
    return sprof

@njit
def calc_breakprofiles(adjlist,weight):
    breakprofile=np.zeros(adjlist.shape[1:],dtype=adjlist.dtype)
    
    for i in range(breakprofile.shape[0]):
        for j in range(breakprofile.shape[1]):
            bsum=0.
            wsum=0.
            count=0
            for k in range(adjlist.shape[0]):
                if adjlist[k,i,j]==adjlist[k,i,j]:
                    bsum+=adjlist[k,i,j]*weight[k]
                    wsum+=weight[k]
                    count+=1
            if count>1:
                breakprofile[i,j]=bsum/wsum
    return breakprofile
            
def add_breakprofile(ref,l,rich,richtemp,method,ib):
    
    bp=rich[method]['break_profiles']
    bm=bp.shape
    #if ib==16:
        #print(ib)
    sign={'tau':1.0,'obs':1.0}
    for ih in range(bm[1]):
        for ip in range(bm[2]):
            
            if ~np.isnan(bp[ib,ih,ip]):
                rich[method]['rich_adjustments'][ih,ip,:ib]-=sign[method]*bp[ib,ih,ip]
                if method=='tau':
                    
                    richtemp[method]['rich_adjusted_fg_dep'][ih,ip,:ref[l]['goodbreaklist'][ib]]+=sign[method]*bp[ib,ih,ip]
                else:
                    richtemp[method]['rich_adjusted_anomalies'][ih,ip,:ref[l]['goodbreaklist'][ib]]+=sign[method]*bp[ib,ih,ip]
                #rich[method]['rich_adjusted_temperatures'][ih,ip,:ref[l]['breaklist'][ib]]+=sign[method]*bp[ib,ih,ip]
    return


#@ifray(ray.remote,RC['richfuture'])
def rich_adjust(l,obj_ref,rich_iter,obj_rich_ref=None):
    sys.path.append(os.path.expanduser('~/python/Rasotools/rasotools/'))
    sys.path.append(os.path.expanduser('~/python/Rasotools/'))
    import warnings
    warnings.filterwarnings("ignore")
    tt=time.time()

    if type(obj_ref) is list: # if called by ray, the first dereferencing is done automatically
        # instead of deferencing the whole list of object references at once, we dereference only those needed -
        # Note that for each stations only at most 500 buddies (out of 1300) are tested.
        ref=[None for i in range(len(obj_ref))]
        ref[l]=ray.get(obj_ref[l])
        o_ref=obj_ref
        #ref=ray.get(obj_ref)
        if rich_iter!=0:
            rich_ref=[None for i in range(len(obj_ref))]
            rich_ref[l]=ray.get(obj_rich_ref[l])
            o_rich_ref=obj_rich_ref
    else: # if called directly from python, the first dereferencing must be done manually
        o_ref=ray.get(obj_ref) # first dereference - gets list ob object references
        
        if rich_iter!=0:
            o_rich_ref=ray.get(obj_rich_ref)
        ref=[None for i in range(len(o_ref))]
        ref[l]=ray.get(o_ref[l]) # second dereference - gets actual dictionaries, here only for tested station
        
        if rich_iter!=0:
            rich_ref=[None for i in range(len(o_ref))]
            rich_ref[l]=ray.get(o_rich_ref[l]) # second dereference - gets actual dictionaries, here only for tested station



    #print('rich_adjust:',l,type(ref),len(ref),type(ref),type(sdist))
    sdidx=np.argsort(ref[l]['sdists'])

    bll=len(ref[l]['goodbreaklist'])
    if bll>0:       
        if ref[l]['goodbreaklist'][-1]==44998: # initial adjustment was applied if this is the case
            bll=len(ref[l]['goodbreaklist'])-1
    
    richens=[]
    richcount=-1
    for daynight_method in ['normal','maxdata']:       
        for weight_distance in RC['weight_distance']: # 3000,5000
            weight=np.exp(-0.01-ref[l]['sdists'][sdidx]*6370./weight_distance)
            for eqsamp in [False, True]:
                for min_neighbours in RC['min_neighbours']:
                    for ri_min_sampsize in RC['ri_min_sampsize']:
                        richcount+=1
                        richens.append(dict(params=dict(daynight_method=daynight_method,weight_distance=weight_distance,eqsamp=eqsamp,min_neighbours=min_neighbours,ri_min_sampsize=ri_min_sampsize),
                                            found=dict(break_dates=['' for i in range(bll)],
                                             dists=['' for i in range(bll)],
                                             weight=weight.copy(),
                                             ref=[list() for i in range(bll)],
                                             refsid=[list() for i in range(bll)]),
                                  tau={},obs={}))
                        bm=(bll,ref[l]['fg_dep'].shape[0],ref[l]['fg_dep'].shape[1])
                        for method in 'tau','obs':
                            richens[-1][method]['buddies']=[list() for i in range(bll)]
                            richens[-1][method]['break_profiles']=np.zeros_like(ref[l]['fg_dep'],shape=bm)
                            richens[-1][method]['rich_adjustments']=np.zeros_like(richens[-1][method]['break_profiles'],shape=(bm[1],bm[2],bm[0]))
    

    richcount=24-1
    for daynight_method in ['normal']:       
        for weight_distance in [RC['weight_distance'][1]]: # 3000,5000
            weight=np.exp(-0.01-ref[l]['sdists'][sdidx]*6370./weight_distance)
            for eqsamp in [False]:
                for min_neighbours in [RC['min_neighbours'][1]]:
                    for ri_min_sampsize in [RC['ri_min_sampsize'][1]]:
                        richcount+=1
                        richens[richcount]=dict(params=dict(daynight_method=daynight_method,weight_distance=weight_distance,eqsamp=eqsamp,min_neighbours=min_neighbours,ri_min_sampsize=ri_min_sampsize),
                                            found=dict(break_dates=['' for i in range(bll)],
                                             dists=['' for i in range(bll)],
                                             weight=weight.copy(),
                                             ref=[list() for i in range(bll)],
                                             refsid=[list() for i in range(bll)]),
                                  tau={},obs={})
                        bm=(bll,ref[l]['fg_dep'].shape[0],ref[l]['fg_dep'].shape[1])
                        for method in 'tau','obs':
                            richens[richcount][method]['buddies']=[list() for i in range(bll)]
                            richens[richcount][method]['break_profiles']=np.zeros_like(ref[l]['fg_dep'],shape=bm)
                            richens[richcount][method]['rich_adjustments']=np.zeros_like(richens[-1][method]['break_profiles'],shape=(bm[1],bm[2],bm[0]))
    

    # Here would be the right place to generate the ensemble. 
    richcount=24-1
    for rich in [richens[richcount+1]]:
        richcount+=1
        richtemp={'tau':{},'obs':{}} # these arrays are for convenience but are too large to be stored in the object store
        richtemp['tau']['rich_adjusted_fg_dep']=ref[l]['fg_dep'].copy()
        richtemp['obs']['rich_adjusted_anomalies']=ref[l]['anomalies'].copy()
        print('richens',richcount,len(richens))
        bp=richens[richcount]['tau']['break_profiles']

        # main loop over breakpoints. 
        for ib in range(bll-1,-1,-1): # test station breaks       
            tib=ref[l]['days'][ref[l]['goodbreaklist'][ib]] # we need days since 1900 here! breaklist values are indices of time series, not days since 1900
            s=0
            refcount=0
            lblen=0
            si=0
            try:
                
                for si in sdidx[1:]: # try ref stations beginning with nearest:
                    if si>=len(ref):
                        continue
                    
                    
                    refcount +=1
                    if ref[l]['sdists'][si]*6370. > weight_distance and refcount>500: # second condition is for test stations at remote places
                        print('not enough buddies found','{:5.1f}'.format(ref[l]['sdists'][si]*6370.),refcount,lblen)
                        break
                    if ref[si] is None:
                        ref[si]=ray.get(o_ref[si])
                        if rich_iter!=0:
                            rich_ref[si]=ray.get(o_rich_ref[si])
                            
                    if  ref[si]['sid'][-5:-3] in ['42','43','48']: # Don't use Indian stations as buddies
                        continue
                    if ref[l]['sid'][-5]=='5' and (ref[si]['sid'][-5]=='5' or ref[si]['sid'][-5:-3] in ['42','43','48']):  # avoid Chinese sondes as buddies for Chinese sondes
                        continue
                    
        
                    sni=ref[si]['days'].shape[0]
                    #print(ref[si]['sid'],'breaklistref',[0],ref[si]['goodbreaklist'],[sni-1])
                    if len(ref[si]['goodbreaklist'])>0:
                        
                        breaklistref=ref[si]['days'][np.concatenate(([0],ref[si]['goodbreaklist'],[sni-1]))] # days since 1900 of breaks of reference station. First and last values added for convenience
                    else:
                        breaklistref=ref[si]['days'][np.concatenate(([0],[sni-1]))] # days since 1900 of breaks of reference station. First and last values added for convenience
                    
                    rbm=np.searchsorted(breaklistref,tib)-1
                    if rbm<breaklistref.shape[0]-1:
                        if rich_iter==0:
                            
                            rmaxlen=breaklistref[rbm+1]-tib 
                        else: 
                            gib=ref[l]['goodbreaklist'][ib] # here we need the index in the time series, not the day
                            sum1=np.sum(~np.isnan(ref[l]['fg_dep'][:,11,gib:gib+RC['mean_maxlen']]),axis=1)
                            if np.any(sum1<3*ri_min_sampsize):                           
                                sum2=np.sum(~np.isnan(ref[l]['fg_dep'][:,11,gib:]),axis=1)
                                if any(sum2<3*ri_min_sampsize):
                                    rmaxlen=np.min([breaklistref[-1]-tib,RC['mean_maxlen']])
                                else:
            
                                    x=np.cumsum(~np.isnan(ref[l]['fg_dep'][:,11,gib:]),axis=1)
                                    y00=np.searchsorted(x[0,:],4*ri_min_sampsize)
                                    y12=np.searchsorted(x[1,:],4*ri_min_sampsize)
                                    mml=np.max([y00,y12])-1
                                        
                                    rmaxlen=np.min([breaklistref[-1]-tib,mml]) # probably a long data gap
                            else:
                                rmaxlen=np.min([breaklistref[-1]-tib,RC['mean_maxlen']])
                            #if ib==11:
                                #print(rmaxlen)
                                
                    else:
                        rmaxlen=np.min([breaklistref[-1]-tib,RC['mean_maxlen']])
                    if rmaxlen<3*ri_min_sampsize:
                        continue
                    if rbm<breaklistref.shape[0]:
                        
                        if ib>0:
                            sibm1=np.searchsorted(ref[si]['days'],ref[l]['goodbreaklist'][ib-1])
                            tibm1=ref[l]['days'][ref[l]['goodbreaklist'][ib-1]] # we need days since 1900 here! breaklist values are indices of time series, not days since 1900
                            #lmaxlen=np.min((tib-breaklistref[rbm],tib-tibm1,RC['mean_maxlen'])) if rich_iter==0 else np.min([tib-tibm1,RC['mean_maxlen']])
                            lmaxlen=np.min((tib-breaklistref[rbm],tib-tibm1)) if rich_iter==0 else np.min([tib-tibm1])
                        else:
                            #lmaxlen=np.min((tib-breaklistref[rbm],RC['mean_maxlen'])) if rich_iter==0 else np.min([tib,RC['mean_maxlen']])            
                            lmaxlen=np.min([tib-breaklistref[rbm]]) if rich_iter==0 else np.min([tib])            
        
                        if lmaxlen<ri_min_sampsize:
                            continue
                    #ldiff=ref[l]['fg_dep'][tib-lmaxlen:tib]-ref[si]['fg_dep'][tib-lmaxlen:tib]
                    #rdiff=ref[l]['fg_dep'][tib:tib+rmaxlen]-ref[si][tib:tib+rmaxlen]
                    # need to do reverse indexing here!
                    lidx=np.searchsorted(ref[l]['days'], (tib-lmaxlen,tib,tib+rmaxlen))
                    sidx=np.searchsorted(ref[si]['days'], (tib-lmaxlen,tib,tib+rmaxlen))
                    
                    if np.any(sidx[1:]-sidx[:-1]<ri_min_sampsize) or np.any(lidx[1:]-lidx[:-1]<ri_min_sampsize):
                        continue
                    
                    sbin=RC['months'][ref[si]['days'][sidx[0]:sidx[2]]]
                    tbin=RC['months'][ref[l]['days'][lidx[0]:lidx[2]]]
                    
                    # define actual series that are compared. In first iteration the test series is compared with 
                    # background departures or anomalies (see definition of rich['tau']['rich_adjusted_fg_dep'] above), 
                    # in the second iteration with the adjusted departures or anomalies from the first iteration
                    # Note also that the tested series starts always as raw series but is updated at each break. The
                    # updated series is in fg_dep (not ref[l]['fg_dep']) or anomalies (not ref[l]['anomalies'])
                    if rich_iter==0:
                        comp={'tau':{'test':richtemp['tau']['rich_adjusted_fg_dep'],'ref':ref[si]['fg_dep']},
                                   'obs':{'test':richtemp['obs']['rich_adjusted_anomalies'],'ref':ref[si]['anomalies']}}
                    else:
                        
                        if not RC['richfuture']:
                            #!!!!! this is the debug setting change to rich_ref['tau']['rich_adjusted_anomalies'] or .. [obs] !!!!!!!!!!!
                            comp={'tau':{'test':richtemp['tau']['rich_adjusted_fg_dep'],
                                         'ref':make_adjusted_series(ref[si]['fg_dep'],ref[si]['adjustments'],ref[si]['goodbreaklist'])},
                                    'obs':{'test':richtemp['obs']['rich_adjusted_anomalies'],
                                            'ref':make_adjusted_series(ref[si]['anomalies'],ref[si]['adjustments'],ref[si]['goodbreaklist'])}}
                        else:
                            #!!!!! this is the full  setting, requires complete first RICH iteration
                            try:
                                
                                #comp={'tau':{'test':richtemp['tau']['rich_adjusted_fg_dep'],'ref':rich_ref[si]['tau']['rich_adjusted_fg_dep']},
                                           #'obs':{'test':richtemp['obs']['rich_adjusted_anomalies'],'ref':rich_ref[si]['obs']['rich_adjusted_anomalies']}}
                                comp={'tau':{'test':richtemp['tau']['rich_adjusted_fg_dep'],
                                             'ref':make_adjusted_series(ref[si]['fg_dep'],rich_ref[si][richcount]['tau']['rich_adjustments'],ref[si]['goodbreaklist'])},
                                           'obs':{'test':richtemp['obs']['rich_adjusted_anomalies'],
                                                  'ref':make_adjusted_series(ref[si]['anomalies'],rich_ref[si][richcount]['obs']['rich_adjustments'],ref[si]['goodbreaklist'])}}
                            except KeyError:
                                print(ref[l]['sid']+': reference '+ref[si]['sid']+' has no rich_adjusted values')
                                continue
                        
                    for method in ['tau','obs']:
                                            
                        ldiff=comp[method]['test'][:,:,lidx[0]:lidx[2]]
                        sdiff=comp[method]['ref'][:,:,sidx[0]:sidx[2]]
                        aprof=rich_calcadjprofile(ldiff,sdiff,
                                                  ref[l]['days'][lidx[0]:lidx[2]],ref[si]['days'][sidx[0]:sidx[2]],
                                                  tbin,sbin,
                                                  lidx[1]-lidx[0],sidx[1]-sidx[0])
                        if np.sum(~np.isnan(aprof[0,:]))>9 or np.sum(~np.isnan(aprof[1,:]))>9:
                            rich[method]['buddies'][ib].append(aprof)
                        if method=='tau':   
                            rich['found']['ref'][ib].append(ref[si]['sid'][-6:])
                            rich['found']['refsid'][ib].append('{:5.0f} km'.format(ref[l]['sdists'][si]*6370.))
                                    
                    lblen=len(rich['tau']['buddies'][ib])        
                    if lblen>=min_neighbours[rich_iter]:
                        for method in 'tau','obs':
                            try:
                                if len(rich[method]['buddies'][ib])>0:
                                    buddies=np.array(rich[method]['buddies'][ib])
                                    rich[method]['break_profiles'][ib,:,:]=calc_breakprofiles(buddies,
                                                                                          rich['found']['weight'][:buddies.shape[0]])
                                else:
                                    rich[method]['break_profiles'][ib,:,:]=np.nan
                            except MemoryError as e:
                                print(e,lblen,ib,rich[method]['break_profiles'].shape,rich['found']['weight'].shape,
                                np.array(rich[method]['buddies'][ib]).shape)
                                rich[method]['break_profiles'][ib,:,:]=np.nan
                                
                            if ref[l]['sid'][-6:]=='012425' and ib==11 and method=='tau':
                                plt.plot(ref[l]['days'][lidx[0]:lidx[2]],comp['tau']['test'][:,:,lidx[0]:lidx[2]][0,3,:])
                                plt.plot(ref[si]['days'][sidx[0]:sidx[2]],comp['tau']['ref'][:,:,sidx[0]:sidx[2]][0,3,:])
                                plt.plot(ref[si]['days'][sidx[1]],0,'ro')
                                plt.plot(ref[l]['days'][lidx[1]],0,'go')
                                
                                plt.plot(ref[l]['days'][lidx[0]:lidx[2]],comp['obs']['test'][:,:,lidx[0]:lidx[2]][0,3,:])
                                plt.plot(ref[si]['days'][sidx[0]:sidx[2]],comp['obs']['ref'][:,:,sidx[0]:sidx[2]][0,3,:])
        
                                plt.savefig('rich_{}_{}_{}_ts.png'.format(richcount,ref[l]['sid'][-6:],ib))
                                plt.semilogy(rich['tau']['buddies'][ib][0][0,:],RC['stdplevs'][RC['pidx']])
                                plt.semilogy(rich['tau']['buddies'][ib][0][1,:],RC['stdplevs'][RC['pidx']])
                                plt.ylim(1000,10)
                                plt.semilogy(rich['obs']['buddies'][ib][0][0,:],RC['stdplevs'][RC['pidx']])
                                plt.semilogy(rich['obs']['buddies'][ib][0][1,:],RC['stdplevs'][RC['pidx']])
                                plt.savefig('rich_{}_{}_{}_prof.png'.format(richcount,ref[l]['sid'][-6:],ib))
                                print('ib')
                            add_breakprofile(ref,l,rich,richtemp,method,ib)
                            #bstart= 0 if ib==0 else ref[l]['goodbreaklist'][ib-1]
                            #for ih in range(bm[1]):
                                #for ip in range(bm[2]):
                                    ##comp[method]['test'][ih,ip,bstart:ref[l]['goodbreaklist'][ib]]-=rich[method]['rich_adjustments'][ih,ip,ib]
                                    #if method=='tau':  
                                        #richtemp[method]['rich_adjusted_fg_dep'][ih,ip,bstart:ref[l]['goodbreaklist'][ib]]+=rich[method]['rich_adjustments'][ih,ip,ib]
                                    #else:
                                        #richtemp[method]['rich_adjusted_anomalies'][ih,ip,bstart:ref[l]['goodbreaklist'][ib]]+=rich[method]['rich_adjustments'][ih,ip,ib]
                        #print(ref[l]['sid'],tib, '{:5.0f} km'.format(sdists[l][si]*6370.))
                        break
                rich['found']['break_dates'][ib]=str(datetime.date(*RC['refdate'])+datetime.timedelta(days=int(tib)))
                rich['found']['dists'][ib]='{:5.0f} km'.format(ref[l]['sdists'][si]*6370.)
                if lblen<min_neighbours[rich_iter]:
                    print(ref[l]['sid'],tib,'ensemble too small, try supplying RAOBCORE estimate',lblen)
                    try:
                        for method in 'tau','obs':
                            #print(rich[method].keys(),ref[l].keys())
                            rich[method]['break_profiles'][ib,:,:]=ref[l]['break_profiles'][ib,:,:]
                            add_breakprofile(ref,l,rich,richtemp,method,ib)
                            #bstart= 0 if ib==0 else ref[l]['goodbreaklist'][ib-1]
                            #for ih in range(bm[1]):
                                #for ip in range(bm[2]):
                                    ##comp[method]['test'][ih,ip,bstart:ref[l]['goodbreaklist'][ib]]-=rich[method]['rich_adjustments'][ih,ip,ib]
                                    #if method=='tau':  
                                        #richtemp[method]['rich_adjusted_fg_dep'][ih,ip,bstart:ref[l]['goodbreaklist'][ib]]+=rich[method]['rich_adjustments'][ih,ip,ib]
                                    #else:
                                        #richtemp[method]['rich_adjusted_anomalies'][ih,ip,bstart:ref[l]['goodbreaklist'][ib]]+=rich[method]['rich_adjustments'][ih,ip,ib]
                                
                    except:
                        print(ref[l]['sid'],tib,'no RAOBCORE estimate could be supplied')
            except MemoryError as e:
                print(e,'adjustment failed for breakpoint at ',tib)
        
        
        print(richcount,l,len(ref),ref[l]['sid'],time.time()-tt)
        print([(rich['found']['break_dates'][i],len(rich['tau']['buddies'][i]),
                                  rich['found']['dists'][i]) for i in range(len(rich['tau']['buddies']))])
    
        
        
    return ray.put(richens)

ray_rich_adjust=ray.remote(rich_adjust)
        
def breakplot(l,obj_ref,obj_rich_ref0=None,obj_rich_ref1=None):
  
    #print(type(obj_ref),type(obj_ref[0]),type(results_ref),type(sdist_ref),type(rich_ref0))
    if type(obj_ref[0]) is ray._raylet.ObjectRef:
        ref=ray.get(obj_ref)
    else:
        ref=obj_ref
        
    rd={}
    if type(obj_rich_ref0[0]) is ray._raylet.ObjectRef: #obj_rich_ref0 is not dict:
        rd['RICH0']=ray.get(obj_rich_ref0)
    else:
        rd['RICH0']=obj_rich_ref0
        
    if obj_rich_ref1: 
        if type(obj_rich_ref0[0]) is ray._raylet.ObjectRef:
            rd['RICH1']=ray.get(obj_rich_ref1)
        else:
            rd['RICH1']=obj_rich_ref1
        
    nbreaks=len(ref[l]['goodbreaklist'])
    col=['b','r']
    spagcol={'tau':'brown','obs':'grey'}
    
    for k in rd.keys():
        for method in ['obs','tau']:
            
            plt.figure(figsize=(10,2+nbreaks//3*6))
            if l>=len(rd['RICH0']):
                rl=0
            else:
                rl=l
            buddies=rd[k][rl][24][method]['buddies']
            for ib in range(nbreaks):
                tib=ref[l]['days'][ref[l]['goodbreaklist'][ib]]
                plt.subplot(nbreaks//3+1,3,ib+1)
                #plt.subplot(1,1,1)
                xlim=np.nanmin(ref[l]['break_profiles'][ib,:,:])-1,np.nanmax(ref[l]['break_profiles'][ib,:,:]+1)
                for ih in 0,1:
                    #print(ref[l]['sid'][-6:],np.array(rd[k][rl][24][method]['buddies'][ib][:3]).shape)
                    ba=np.array(buddies[ib][:3])
                    if ba.ndim==3:
                        
                        plt.semilogy(np.array(buddies[ib][:3])[:,ih,:].T,RC['stdplevs'][RC['pidx']],color=spagcol[method],alpha=0.7)

                    plt.semilogy(ref[l]['break_profiles'][ib,ih,:],RC['stdplevs'][RC['pidx']],label='RAOBCORE',color=col[ih],linewidth=3,alpha=0.6)
                    
                    plt.semilogy(rd[k][rl][24][method]['break_profiles'][ib,ih,:],RC['stdplevs'][RC['pidx']],label=k+'-'+method,color=col[ih],linestyle='--')
                    tmin=np.min((3,len(rd[k][0][24]['found']['refsid'][ib])))
                    for i in range(tmin):
                        plt.text(xlim[0]+0.3,500+100*i, rd[k][0][24]['found']['ref'][ib][i]+','+rd[k][0][24]['found']['refsid'][ib][i],fontsize='x-small')
                
                plt.ylim(1000.,10.)
                plt.xlim(*xlim)
                plt.title(ref[l]['sid'][-6:]+' '+str(ib)+' '+str(datetime.date(*RC['refdate'])+datetime.timedelta(days=int(tib))))
                plt.xlabel('K')
                plt.ylabel('hPa')
                plt.legend()
                
            plt.tight_layout()
            plt.savefig(ref[l]['sid'][-6:]+'_'+k+'-'+method+'_breakprofiles.png')
            plt.close()
    
    return

ray_breakplot=ray.remote(breakplot)

@njit
def findstart(fg_dep,minlen):
    i=fg_dep.shape[2]-1
    igood=np.zeros(2,dtype=np.int32)
    isave=np.zeros(2,dtype=np.int32)
    while i>0 and (igood[0]<minlen or igood[1]<minlen):
        for j in range(fg_dep.shape[0]):
            if fg_dep[j,11,i]==fg_dep[j,11,i]:
                igood[j]+=1
                if igood[j]==minlen:
                    isave[j]=i
        i-=1
        
    istart=i
    #if igood[0]>=minlen and igood[1]<minlen/8:
        #istart=isave[0]
    #if igood[1]>=minlen and igood[0]<minlen/8:
        #istart=isave[1]
    return istart,igood

def initial_adjust(l,obj_ref,obj_rich_ref=None,initial_composite_ref=None,rich_initial_composite_ref=None):
    sys.path.append(os.path.expanduser('~/python/Rasotools/rasotools/'))
    sys.path.append(os.path.expanduser('~/python/Rasotools/'))
    
    tt=time.time()

    if type(obj_ref) is list: # if called by ray, the first dereferencing is done automatically
        # instead of deferencing the whole list of object references at once, we dereference only those needed -
        # Note that for each stations only at most 500 buddies (out of 1300) are tested.
        ref=[None for i in range(len(obj_ref))]
        ref[l]=ray.get(obj_ref[l])
        #ref=ray.get(obj_ref)
        if obj_rich_ref:
            rich_ref=[None for i in range(len(obj_ref))]
            rich_ref[l]=ray.get(obj_rich_ref[l])
        if initial_composite_ref:
            initial_composite=initial_composite_ref
        else:
            initial_composite=None
            
        if rich_initial_composite_ref:
            rich_initial_composite=rich_initial_composite_ref
        else:
            rich_initial_composite=None
    else: # if called directly from python, the first dereferencing must be done manually
        ref=ray.get(obj_ref) # first dereference - gets list ob object references
        ref=ray.get(ref) # second dereference - gets actual dictionaries
        if obj_rich_ref:
            rich_ref=ray.get(obj_rich_ref) # first dereference - gets list ob object references
            rich_ref=ray.get(rich_ref) # second dereference - gets actual dictionaries
            
            #rich_ref=[None for i in range(len(obj_ref))]
            #rich_ref[l]=ray.get(obj_rich_ref[l])
        
        if initial_composite_ref:
            initial_composite=initial_composite_ref
            pass
            #initial_composite=ray.get(initial_composite_ref)
        else:
            initial_composite=None
        if rich_initial_composite_ref:
            rich_initial_composite=rich_initial_composite_ref #ray.get(rich_initial_composite_ref)
        else:
            rich_initial_composite=None

        #print('rich_adjust:',l,type(ref),len(ref),type(ref),type(sdist))
    sdidx=np.argsort(ref[l]['sdists'])

    bll=len(ref[l]['goodbreaklist'])
    
    sh=ref[l]['fg_dep'].shape[:2]+(bll,)
    RAOBini={}
    RAOBini['buddies']=[list() for i in range(bll)]
    
    if obj_rich_ref:
        richensini=[]
        richcount=-1
        for i in range(32):
            richcount+=1
            richensini.append(dict(tau={},obs={}))
            for method in 'tau','obs':
                richensini[-1][method]['inibuddies']=[]
                #richensini[-1][method]['initial_adjustments']=np.zeros_like(ref[l]['fg_dep'],shape=ref[l]['fg_dep'].shape[:2])
    


    # add initial adjustments here
    if RC['initial_adjustments']=='era5bc' and int(ref[l]['days'][-1]/365.25)>116:
        start=np.max((ref[l]['days'].shape[0]-RC['mean_maxlen'],np.searchsorted(ref[l]['days'],int(115.5*365.25))))
        if ref[l]['days'].shape[0]-start>RC['snht_min_sampsize']:
            
            ini=np.zeros(ref[l]['era5bc'].shape[:2])
            RAOBini['initial_adjustments']= -ndnanmean2(ref[l]['era5bc'][:,:,start:],ini,RC['snht_min_sampsize'])
            RAOBini['initial_adjustments'][np.isnan(RAOBini['initial_adjustments'])]=0.

    elif RC['initial_adjustments']=='era5':
        start=np.max((0,ref[l]['days'].shape[0]-RC['mean_maxlen']))
        if ref[l]['days'].shape[0]-start>RC['snht_min_sampsize']:
            
            ini=np.zeros(ref[l]['adjusted_fg_dep'].shape[:2])
            RAOBini['initial_adjustments']= ndnanmean2(ref[l]['adjusted_fg_dep'][:,:,start:],ini,RC['snht_min_sampsize'])
            RAOBini['initial_adjustments'][np.isnan(RAOBini['initial_adjustments'])]=0.

                        
    #elif RC['initial_adjustments']=='era5neighbours':
    if True: 
        if ref[l]['lastsonde']: # and RC['initial_adjustments'] not in ['era5','era5bc']:
            #start,igood=findstart(ref[l]['fg_dep'],RC['mean_maxlen'])
            #ini=np.zeros(ref[l]['adjusted_fg_dep'].shape[:2])
            #RAOBini['initial_adjustments']= ndnanmean2(ref[l]['adjusted_fg_dep'][:,:,start:],ini,RC['snht_min_sampsize'])
            #RAOBini['initial_adjustments'][np.isnan(RAOBini['initial_adjustments'])]=0.
            RAOBini['initial_adjustments']=np.zeros_like(ref[l]['fg_dep'],shape=ref[l]['fg_dep'].shape[:2])
            if obj_rich_ref:
                for i in range(32):
                    for method in 'tau','obs':
                        richensini[i][method]['initial_adjustments']=np.zeros_like(ref[l]['fg_dep'],shape=ref[l]['fg_dep'].shape[:2])
            pass
        else:
            #start=np.max((0,ref[l]['days'].shape[0]-RC['mean_maxlen']))
            start,igood=findstart(ref[l]['fg_dep'],RC['mean_maxlen'])
            spagini=[]
            if any(igood>RC['snht_min_sampsize']):
                refcount=0
                for si in sdidx[1:]: # try ref stations beginning with nearest:
                    if si>=len(ref):
                        continue
                        
                    refcount +=1
                    if ref[l]['sdists'][si]*6370. > RC['weight_distance'][0] and refcount>500: # second condition is for test stations at remote places
                        print('not enough buddies found','{:5.1f}'.format(ref[l]['sdists'][si]*6370.),refcount)
                        break

                    if ref[si] is None:
                        ref[si]=ray.get(obj_ref[si])

                    try:

                        sistart=np.searchsorted(ref[si]['days'],ref[l]['days'][start])
                        if ref[si]['days'].shape[0]-sistart>RC['snht_min_sampsize']:
                            ini=np.zeros(ref[si]['adjusted_fg_dep'].shape[:2])
                            #spagest=ndnanmean2(ref[si]['adjusted_fg_dep'][:,:,sistart:],ini,RC['snht_min_sampsize'])
                            if initial_composite:
                                if initial_composite[si]:
                                    ic=ray.get(initial_composite[si])
                                    
                                    if 'initial_adjustments' in ic.keys():
                                        spagest=ndnanmean2(ref[si]['adjusted_fg_dep'][:,:,sistart:],ini,RC['snht_min_sampsize'])
                                        for ih in range(spagest.shape[0]):
                                            for ip in range(spagest.shape[1]):
                                                ia=ic['initial_adjustments'][ih,ip]
                                                if ia==ia:
                                                    spagest[ih,ip]-=ia
                                        spagini.append(spagest)
                            
                    
                    except MemoryError as e:
                        print(e)
                        pass
                    if len(spagini)>RC['min_neighbours'][0][1]:
                        spaginimean=np.nanmean(np.array(spagini),axis=0)
                        break
                        
                if len(spagini)>RC['min_neighbours'][0][1]:
                    ini=np.zeros(ref[l]['adjusted_fg_dep'].shape[:2])
                    RAOBini['initial_adjustments']= ndnanmean2(ref[l]['adjusted_fg_dep'][:,:,start:],ini,RC['snht_min_sampsize'])-spaginimean
                    RAOBini['initial_adjustments'][np.isnan(RAOBini['initial_adjustments'])]=0.
                else:
                    print(ref[l]['sid'][-6:]+': no initial composite found, substituting mean departure')
                    ini=np.zeros(ref[l]['adjusted_fg_dep'].shape[:2])
                    RAOBini['initial_adjustments']= ndnanmean2(ref[l]['adjusted_fg_dep'][:,:,start:],ini,RC['snht_min_sampsize'])
                    RAOBini['initial_adjustments'][np.isnan(RAOBini['initial_adjustments'])]=0.
                    
                # now the RICH initial adjustments
                if obj_rich_ref:
                    refcount=0
                    firstbuddy=True                                            
                    ri_test_adjusted_fg_dep=[{'tau':None,'obs':None} for i in range(len(richensini))]
                    
                    for si in sdidx[1:]: # try ref stations beginning with nearest:
                        if si>=len(ref):
                            continue
                            
                        refcount +=1
                        if ref[l]['sdists'][si]*6370. > RC['weight_distance'][0] and refcount>500: # second condition is for test stations at remote places
                            print('not enough buddies found','{:5.1f}'.format(ref[l]['sdists'][si]*6370.),refcount)
                            break
        
                        if rich_ref[si] is None:
                            rich_ref[si]=ray.get(obj_rich_ref[si])
                        if ref[si] is None:
                            ref[si]=ray.get(obj_ref[si])
                        if rich_initial_composite:
                            if rich_initial_composite[si]:
                                ric=ray.get(rich_initial_composite[si])
        
                        try:
        
                            sistart=np.searchsorted(ref[si]['days'],ref[l]['days'][start])
                            if ref[si]['days'].shape[0]-sistart>RC['snht_min_sampsize']:
                                    
                                ri_adjusted_fg_dep=[{'tau':None,'obs':None} for i in range(len(richensini))]
                                for iens in [24]: #range(len(richensini)):
                                    for method in 'tau','obs':
                                        # TO BE CHANGED - currently uses same initial estimate as RAOBCORE
                                        if rich_initial_composite:
                                            if rich_initial_composite[si]:
                                                if 'initial_adjustments' in ric[1][iens][method].keys():
                                                    
                                                    if firstbuddy:    # calculate rich adjusted test time series                               
                                                        ini=np.zeros(ref[l]['adjusted_fg_dep'].shape[:2])
                                                        ri_test_adjusted_fg_dep[iens][method]=make_adjusted_series(ref[l]['fg_dep'],rich_ref[l][iens][method]['rich_adjustments'],
                                                                            ref[l]['goodbreaklist'])
                                                        if iens==len(richensini)-1 and method=='obs':
                                                            firstbuddy=False

                                                    # calculate rich adjusted reference time series
                                                    ini=np.zeros(ref[si]['adjusted_fg_dep'].shape[:2])
                                                    ri_adjusted_fg_dep[iens][method]=make_adjusted_series(ref[si]['fg_dep'],rich_ref[si][iens][method]['rich_adjustments'],
                                                                                            ref[si]['goodbreaklist'])
                                                    rspagest=ndnanmean2(ri_adjusted_fg_dep[iens][method][:,:,start:],ini,RC['snht_min_sampsize'])
                                                    for ih in range(rspagest.shape[0]):
                                                        for ip in range(rspagest.shape[1]):
                                                            ia=ric[1][iens][method]['initial_adjustments'][ih,ip]
                                                            if ia==ia:
                                                                rspagest[ih,ip]-=ia
                                                    try:
                                                        if np.sum(~np.isnan(rspagest))>5:
                                                            
                                                            richensini[iens][method]['inibuddies'].append(rspagest)
                                                    except:
                                                        pass
                                        #if iens==24:
                                            #print('inibuddies')
                        
                        except MemoryError as e:
                            print(e)
                            pass
                        if len(richensini[24]['obs']['inibuddies'])>RC['min_neighbours'][0][1]:
                            for iens in [24]: #range(len(richensini)):
                                for method in 'tau','obs':
                                    # TO BE CHANGED - currently uses same initial estimate as RAOBCORE
                                    #ini=np.zeros(ref[si]['adjusted_fg_dep'].shape[:2])
                                    #rspagest=ndnanmean2(rich_ref[si]['adjusted_fg_dep'][:,:,start:],ini,RC['snht_min_sampsize'])
                                    richensini[iens][method]['spaginimean']=np.nanmean(np.array(richensini[iens][method]['inibuddies']),axis=0)
                            break
                            
                    if len(richensini[24]['obs']['inibuddies'])>RC['min_neighbours'][0][1]:
                        for iens in [24]: # range(len(richensini)):
                            for method in 'tau','obs':
                                # needs to be changed, currently takes RAOBCORE mean
                                ini=np.zeros(ref[l]['adjusted_fg_dep'].shape[:2])
                                richensini[iens][method]['initial_adjustments']=ndnanmean2(ri_test_adjusted_fg_dep[iens][method][:,:,start:],ini,RC['snht_min_sampsize'])-richensini[iens][method]['spaginimean']
                                richensini[iens][method]['initial_adjustments'][np.isnan(RAOBini['initial_adjustments'])]=0.
                    else:
                        print(ref[l]['sid'][-6:]+': no initial RICH composite found, substituting mean departure')
                        for iens in range(len(richensini)):
                            for method in 'tau','obs':
                                ini=np.zeros(ref[l]['adjusted_fg_dep'].shape[:2])
                                richensini[iens][method]['initial_adjustments']=ndnanmean2(ref[l]['adjusted_fg_dep'][:,:,start:],ini,RC['snht_min_sampsize'])
                                richensini[iens][method]['initial_adjustments'][np.isnan(RAOBini['initial_adjustments'])]=0.
                        
    elif RC['initial_adjustments']=='gpsro':
        start=np.max((0,ref[l]['days'].shape[0]-RC['mean_maxlen']))
        if ref[l]['days'].shape[0]-start>RC['snht_min_sampsize']:
            ini=np.zeros(ref[l]['adjusted_fg_dep'].shape[:2])
            RAOBini['initial_adjustments']= ndnanmean2(ref[l]['adjusted_fg_dep'][:,:,start:],ini,RC['snht_min_sampsize'])
            RAOBini['initial_adjustments'][np.isnan(RAOBini['initial_adjustments'])]=0.


    elif RC['initial_adjustments']=='rharm':
        start=np.max((0,ref[l]['days'].shape[0]-RC['mean_maxlen']))
        if ref[l]['days'].shape[0]-start>RC['snht_min_sampsize']:
            
            ini=np.zeros(ref[l]['adjusted_fg_dep'].shape[:2])
            RAOBini['initial_adjustments']= ndnanmean2(ref[l]['adjusted_fg_dep'][:,:,start:],ini,RC['snht_min_sampsize'])
            RAOBini['initial_adjustments'][np.isnan(RAOBini['initial_adjustments'])]=0.

    RAOBini['initial_adjustment_method']=RC['initial_adjustments']
    print(ref[l]['sid'][-6:],time.time()-tt)
    if obj_rich_ref:
        return ray.put((RAOBini,richensini))
    else:
        return ray.put(RAOBini)
        

ray_initial_adjust=ray.remote(initial_adjust)

'''
RAOBCORE/RICH python - main script
Leopold Haimberger, 6.6. 2022

In order to have a more flexible code base, RAOBCORE/RICH (Haimberger et al. 2012)
has been completely rewritten in python3 (the original version is in Fortran90)

It is not a mere translation, there are also significant differences.
Most notably, data gaps are handled differently. In the FORTRAN version, segments of series with
gaps larger than a certain size (2 years) were treated separately whereas now, data gaps are bridged.


'''

if __name__ == '__main__':
    
    
    #prolog()
    if len(sys.argv)>1:
        pattern=sys.argv[1]
        pkl='RAOBx'
    else:
        pattern=''
        pkl='RAOBx'
    process = psutil.Process(os.getpid())
    files = glob.glob('/users/staff/leo/fastscratch/rise/1.0/exp02/[01]*/feedbackmerged'+pattern+'*.nc')[:]
    tt=time.time()
    hadmed= read_hadCRUT5('/users/staff/leo/fastscratch/rise/1.0/common/','HadCRUT.5.0.1.0',RC['refdate'])
    # read data and analyze breakpoints
    #func=partial(RAOB_findbreaks,'SNHT')
    #with Pool(20) as P:        
        #finfo=list(P.map(func,files))
        
    #finfo=[i for i in finfo if i]
    
    ray.init(num_cpus=RC['CPUs'])
    #ray.init(address='131.130.157.11:49705', _redis_password='5241590000000000')
    #func=partial(RAOB_findbreaks,'SNHT')
    #with Pool(20) as P:        
        #obj_ref=list(P.map(func,files))
    
    if not RC['findfuture']:  
        func=partial(RAOB_findbreaks,'SNHT') # Binseg
        obj_ref=list(map(func,files))
        obj_ref = [i for i in obj_ref if i]
        finfo=ray.get(obj_ref)
    else:
        futures = []
        for fn in files:
            futures.append(ray_RAOB_findbreaks.remote('SNHT',fn))
        obj_ref  = ray.get(futures)
    obj_ref = [i for i in obj_ref if i]
    
    break_sum=np.sum([len(l['breaklist']) for l in ray.get(obj_ref)])
    if RC['plot']:
        
        breaks=np.zeros(45000,dtype=np.int)
        for l in ray.get(obj_ref):
            breaks[np.array(l['days'][l['breaklist']-1])]+=1
            if len(l['breaklist'])>30:
                print(l['sid'],len(l['breaklist']))
            
        plt.plot(np.arange(45000)[::30]/365.25,[np.sum(breaks[i:i+30]) for i in range(0,45000,30)])
        plt.title('Detected breaks: {}'.format(break_sum))
        plt.savefig('numberofbreaks.png')
        
    print('Total number of breaks',break_sum)
    
    print(time.time()-tt)
    # calculate distance between stations
    sdist_ref=stdists(ray.get(obj_ref))
    print(time.time()-tt)

    
    if not RC['findfuture']:  
        func=partial(RAOB_adjustbreaks,'SNHT')
        results=list(map(func,obj_ref))
        results_ref=results
    else:
        futures=[]
        results_ref=[]
        for l in range(len(obj_ref)):
            futures.append(ray_RAOB_adjustbreaks.remote('SNHT',sdist_ref,obj_ref[l],l))
        obj_ref = ray.get(futures)
    
    del futures
    
    break_sum=np.sum([len(l['goodbreaklist']) for l in ray.get(obj_ref)])
    print('Total number of "good" breaks',break_sum)

    
    sids=[f['sid'] for f in ray.get(obj_ref)] 
    if pattern=='':
        pattern='012425'
    lx= sids.index('feedbackmerged'+pattern)
    
    single_obj_ref=ray.put(obj_ref) 
    #single_results_ref=ray.put(results_ref)

    print('before initial_adjust',time.time()-tt)
    if True:
        lgood=[]
        lrecent=[]
        lneedscomposite=[]
        initial_composite=[{} for i in range(len(obj_ref))]
        for l in range(len(obj_ref)):
            d=ray.get(obj_ref[l])
            if d['lastsonde']==1:
                lgood.append(l)
            elif d['days'][-1]>int(110*365.25):
                lrecent.append(l)
            else:
                lneedscomposite.append(l)
            
        ifile='exp03'.join(files[0].split('exp02'))
        #for lx in lrecent:
            #initial_adjust_RAOB=initial_adjust(lx,single_obj_ref,single_results_ref,sdist_ref)  
            #save_monmean(ifile,0,[obj_ref[lx]],[results_ref[lx]],
                         #initial_adjust_RAOB=[initial_adjust_RAOB])
        #for lx in lneedscomposite:
            #initial_adjust_RAOB=initial_adjust(lx,single_obj_ref,single_results_ref,sdist_ref)  
            ##breakplot(lx,obj_ref,results_ref,sdist_ref,obj_rich_ref0=[rich_ref0[lx]],obj_rich_ref1=[rich_ref1])
            #save_monmean(ifile,0,[obj_ref[lx]],[results_ref[lx]],
                         #initial_adjust_RAOB=[initial_adjust_RAOB])
                         
# for some values of RC['initial_adjust'] even the "good" sondes are initial-adjusted
    #for l in lgood:
        #initial_adjust_RAOB=initial_adjust(l,single_obj_ref,single_results_ref,sdist_ref)  
    futures=[]
    for l in lgood: #range(len(obj_ref)):
        futures.append(ray_initial_adjust.remote(l,single_obj_ref))
    iadjustlist = ray.get(futures)
    for i in range(len(lgood)): 
        initial_composite[lgood[i]]=copy.copy(iadjustlist[i])
    #initial_composite_ref=ray.put(initial_composite)
    initial_composite_ref=initial_composite

# next those with records up to recent times but susptected biases are initial-adjusted
    for l in lrecent:
        if ray.get(obj_ref[l])['sid'][-6:] in ('022522'):
            
            initial_adjust_RAOB=initial_adjust(l,single_obj_ref,initial_composite_ref=initial_composite_ref)  
    futures=[]
    for l in lrecent: #range(len(obj_ref)):
        futures.append(ray_initial_adjust.remote(l,single_obj_ref,initial_composite_ref=initial_composite_ref))
    iadjustlist = ray.get(futures)
    for i in range(len(lrecent)): 
        initial_composite[lrecent[i]]=copy.copy(iadjustlist[i])
    #initial_composite_ref=ray.put(initial_composite)
    initial_composite_ref=initial_composite
    
    print('after initial_adjust',time.time()-tt)

    for l in lneedscomposite:
        if ray.get(obj_ref[l])['sid'][-6:] in ('022522'):
            
            initial_adjust_RAOB=initial_adjust(l,single_obj_ref,single_results_ref,sdist_ref,initial_composite_ref=initial_composite_ref)  
            ###breakplot(lx,obj_ref,results_ref,sdist_ref,obj_rich_ref0=[rich_ref0[lx]],obj_rich_ref1=[rich_ref1])

# finally those records that end early (before 2010) are initial-adjusted
    futures=[]
    for l in lneedscomposite: #range(len(obj_ref)):
        futures.append(ray_initial_adjust.remote(l,single_obj_ref,initial_composite_ref=initial_composite_ref))
    iadjustlist = ray.get(futures)

    for i in range(len(lneedscomposite)): 
        initial_composite[lneedscomposite[i]]=copy.copy(iadjustlist[i])
    #initial_adjust_ref=[ray.put(initial_composite[l]) for l in range(len(initial_composite))]
    initial_adjust_ref=[initial_composite[l] for l in range(len(initial_composite))]
    
    save_monmean(ifile,0,[obj_ref[lx]],initial_adjust_RAOB=[initial_adjust_ref[lx]])
    print('before write',time.time()-tt)
    futures=[]
    for l in range(len(obj_ref)):
        futures.append(ray_save_monmean.remote(ifile,0,[obj_ref[l]],initial_adjust_RAOB=[initial_adjust_ref[l]]))
    x = ray.get(futures)
    print('after write',time.time()-tt)
    
    func=partial(save_monmean,files[0])
    #res=list(map(func, [obj_ref[l]], [results_ref[l]]))
    #res=list(map(func, obj_ref, results_ref))

    # it is more efficient to have the whole list of object_references also as object_reference. 
    # Instead of copying the whole list ob object references for each station when calling rich_adjust, only the
    # object reference to the list of object references is copied.

    

    #single_sdist_ref=ray.put(sdist_ref)
    print('before first RICH iteration',time.time()-tt)
    futures=[]
    if True:
        for l in range(len(obj_ref)):
            futures.append(ray_rich_adjust.remote(l,single_obj_ref,0))
        
        rich_ref0=ray.get(futures)   
    else:
        rich_ref0=[]
        for l in range(len(obj_ref)):
            rich_ref0.append(rich_adjust(l,single_obj_ref,0))
        
    
    single_rich_ref0=ray.put(rich_ref0) 
    
    #if not RC['richfuture']:    
        #RC['richfuture']=True
    if False:
        rich_ref0=rich_adjust(lx,single_obj_ref,0)
        breakplot(lx,obj_ref,obj_rich_ref0=[rich_ref0])
        single_rich_ref0=ray.put([rich_ref0])
    elif True:
        rich_ref1=rich_adjust(lx,single_obj_ref,1,obj_rich_ref=single_rich_ref0)  
        breakplot(lx,obj_ref,obj_rich_ref0=[rich_ref0[lx]],obj_rich_ref1=[rich_ref1])
        save_monmean(ifile,0,[obj_ref[lx]],rich_ref0=[rich_ref0[lx]],rich_ref1=[rich_ref1])

    print('before second RICH iteration',time.time()-tt)
    #rich_ref1=[]
    #for l in range(len(obj_ref)):
        #rich_ref1.append(rich_adjust(l,single_obj_ref,single_results_ref,sdist_ref,1,obj_rich_ref=single_rich_ref0))
    futures=[]
    if True:
        for l in range(len(obj_ref)):
            futures.append(ray_rich_adjust.remote(l,single_obj_ref,1,obj_rich_ref=single_rich_ref0))
        rich_ref1=ray.get(futures)
    else:
        rich_ref1=[]
        for l in range(len(obj_ref)):
            rich_ref1.append(rich_adjust(l,single_obj_ref,1,obj_rich_ref=single_rich_ref0))
        
    single_rich_ref1=ray.put(rich_ref1) 
    
# For debugging:
    #futures=[]
    #for l in range(len(obj_ref)):
        #futures.append(ray_breakplot.remote(l,obj_ref,results_ref,sdist_ref,obj_rich_ref0=[rich_ref0[l]],obj_rich_ref1=[rich_ref1[l]]))
    #dum=ray.get(futures)

# initial adjustment of RICH composite        
    rich_initial_composite=[{} for i in range(len(obj_ref))]
    
    futures=[]
    for l in lgood: #range(len(obj_ref)):
        futures.append(ray_initial_adjust.remote(l,single_obj_ref,obj_rich_ref=single_rich_ref1))
    irichadjustlist = ray.get(futures)
    for i in range(len(lgood)): 
        rich_initial_composite[lgood[i]]=copy.copy(irichadjustlist[i])
    #rich_initial_composite_ref=ray.put(rich_initial_composite)
    rich_initial_composite_ref=rich_initial_composite

    for l in lrecent:
        if ray.get(obj_ref[l])['sid'][-6:] in ('034247'):
            initial_adjust_x=initial_adjust(l,single_obj_ref,
                                        obj_rich_ref=single_rich_ref1,initial_composite_ref=initial_composite_ref,
                                        rich_initial_composite_ref=rich_initial_composite_ref)  
    futures=[]
    for l in lrecent: #range(len(obj_ref)):
        futures.append(ray_initial_adjust.remote(l,single_obj_ref,obj_rich_ref=single_rich_ref1,initial_composite_ref=initial_composite_ref,
                                        rich_initial_composite_ref=rich_initial_composite_ref))
    irichadjustlist = ray.get(futures)
    for i in range(len(lrecent)): 
        rich_initial_composite[lrecent[i]]=copy.copy(irichadjustlist[i])
    #rich_initial_composite_ref=ray.put(rich_initial_composite)
    rich_initial_composite_ref=rich_initial_composite
    
    print('after initial_adjust recent',time.time()-tt)

    # test rich initial composite
    #for l in lneedscomposite:
        #initial_adjust_x=initial_adjust(l,single_obj_ref,single_results_ref,sdist_ref,
                                        #obj_rich_ref=single_rich_ref1,initial_composite_ref=initial_composite_ref,
                                        #rich_initial_composite_ref=rich_initial_composite_ref)  
            ###breakplot(lx,obj_ref,results_ref,sdist_ref,obj_rich_ref0=[rich_ref0[lx]],obj_rich_ref1=[rich_ref1])

    futures=[]
    for l in lneedscomposite: #range(len(obj_ref)):
        futures.append(ray_initial_adjust.remote(l,single_obj_ref,
                                        obj_rich_ref=single_rich_ref1,initial_composite_ref=initial_composite_ref,
                                        rich_initial_composite_ref=rich_initial_composite_ref))
    irichadjustlist = ray.get(futures)

    for i in range(len(lneedscomposite)): 
        rich_initial_composite[lneedscomposite[i]]=copy.copy(irichadjustlist[i])

# this object list is used for final writing.
    #rich_initial_adjust_ref=[ray.put(rich_initial_composite[l]) for l in range(len(rich_initial_composite))]
    rich_initial_adjust_ref=[rich_initial_composite[l] for l in range(len(rich_initial_composite))]

    print('after initial_adjust needscomposite',time.time()-tt)

    save_monmean(ifile,0,[obj_ref[lx]],rich_ref1=[rich_ref1[lx]],
                 initial_adjust_RAOB=[initial_adjust_ref[lx]],
                 initial_adjust_rich=[rich_initial_adjust_ref[lx]])

    print('before write',time.time()-tt)
    futures=[]
    for l in range(len(obj_ref)):
        #futures.append(ray_save_monmean.remote(files[0],l,single_obj_ref,single_results_ref,rich_ref0=single_rich_ref0,rich_ref1=single_rich_ref1))
        futures.append(ray_save_monmean.remote(ifile,0,[obj_ref[l]], rich_ref1=[rich_ref1[l]],
                                               initial_adjust_RAOB=[initial_adjust_ref[l]],
                                               initial_adjust_rich=[rich_initial_adjust_ref[l]]))
    
    res = ray.get(futures)
        

        
    print('Mem [MiB]',process.memory_info().rss//1024//1024)
    print(time.time()-tt)
    ray.shutdown()
    print('end')
     
