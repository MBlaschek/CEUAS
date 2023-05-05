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
import xarray as xr
#sys.path.append(os.getcwd()+'/../common/rasotools/')
#sys.path.append(os.path.expanduser('../common/Rasotools/rasotools/'))
#sys.path.append(os.path.expanduser('../common/Rasotools/'))
#sys.path.append(os.getcwd()+'/../adjust/rasotools/')
from utils import tdist,extract,rmeanw
from anomaly import ndnanmean2,nnanmean,danomaly
import pyRAOBCORE_numbamod as pnm  
from harvest_convert_to_netCDF import write_dict_h5
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

from add_solarangle_adjustments import add_solelev
ray_add_solelev=ray.remote(add_solelev)

#sys.path.append(os.getcwd()+'/../cds-backend/code/')
#import cds_eua3 as eua

#import zarr
#import dask
from timeit import default_timer as timer
from numba import njit,prange,typeof
import cds_eua4 as eua
import SNHT

#plt.rcParams['lines.linewidth'] = 3

#RAOBCORE constants
def RC_ini(RCfile):
    
    try:
        
        with open(os.path.expandvars(RCfile)) as f:
            RC = json.load(f)
    except Exception as e:
        
        RC=dict(
            snht_maxlen = 1460,
            snht_min_sampsize = 80,
            snht_increment = 30,
            miss_val=math.nan,
            stdplevs = (10.0, 20.0, 30.0, 50.0, 70.0, 100.0, 150.0, 200.0, 
                                 250.0, 300.0, 400.0, 500.0, 700.0, 850.0, 925.0, 1000.0),
            pidx= (0,1,2,3,4,5,6,7,8,9,10,11,12,13, 14, 15),
            )
        
        RC['refdate']=(1900,1,1)
        RC['lastdate']=(2023,1,1)
        RC['snht_thresh'] = 25
        RC['mon_thresh'] = 0
        RC['mean_maxlen'] = RC['snht_maxlen']*3
        RC['plot']=False
        RC['CPUs']=25
        RC['min_neighbours']=[[3,10],[10,30]]
        RC['weight_distance']=[3000,5000]
    #    RC['ri_min_sampsize']=[330,30]
        RC['ri_min_sampsize']=[330,80]
        RC['initial_adjustments']='era5neighbours' #'rharm' #'era5neighbours' #'era5' #'era5bc' #'era5neighbours'
        #RC['transdict']={'montemp':'mtemperatures','rasocorrmon':'madjustments','goodmon':'goodmon'}
        RC['transdict']={'montemp':'mini_adjusted_temperatures','rasocorrmon':'madjustments','goodmon':'goodmon'}
        RC['rharm'] = 'rcuon'
        RC['rich'] = True
        RC['richfuture']=True 
        RC['findfuture']=True
        RC['goodsondes']=[141,142,123,124,125,113,114,79,80,81,70,152,177,183] # from WMO manual on codes, 2019 edition
        RC['apriori_prob_default']=0.01
        RC['exp']='exp00'
        RC['CUON'] = 'converted_v11/'
        RC['fgdepname'] = 'era5_fgdep'
        RC['version'] = '1.8' #1.9.0'
        RC['add_solelev'] = False
        RC['write_to_backend'] = False
        
        with open(os.path.expandvars(RCfile),'w') as f:
            json.dump(RC,f, indent=4, sort_keys=True)
            
    refdate=datetime.datetime(*RC['refdate'])
    lastdate=datetime.datetime(*RC['lastdate'])
    numdays=(lastdate-refdate).days
    dates = [refdate + datetime.timedelta(days=x) for x in range(numdays)]
    RC['years']= [ x.year for x in dates]
    RC['months']= [ x.month for x in dates]
    RC['days']= [ x.day for x in dates]
    for k in 'years','months','days','stdplevs','pidx':
        if type(k) is int:
            RC[k]=np.array(RC[k],dtype=np.int32)
        else:
            RC[k]=np.array(RC[k])
        if k in  ('months','days'):
            RC[k]-=1
    
    RC['tidx'] = np.where(RC['days']==0)[0]
    #RC['plot']=True
    RC['savetsa']=True
    
    return RC



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
    fn=path+prefix+".analysis.anomalies.ensemble_mean.nc.1"
    fvar='tas_mean'
    fevar='tas'
    try:
        with netCDF4.Dataset(fn,"r") as f:
            f.set_auto_mask(False)
            sdate=f.variables['time'].getncattr('units')
            index=(refdate[0]-int(sdate.split('-')[0].split()[-1]))*12
            hadmed=f.variables[fvar][index:,:,:]
            hadmed[np.abs(hadmed)>1.e29]=np.nan
            ds1900=f.variables['time'][index:]-f.variables['time'][index]
            dslon=f.variables['longitude'][:]
            dslat=f.variables['latitude'][:]
    except:
        print((fn+' not found'))
    print('hadmed',time.time()-t1)

    return {'hadmed':hadmed,'time':ds1900,'lat':dslat,'lon':dslon} #,hadtem,hadens

#@njit
def do_test(fg_dep,temperature,anomaly,hadcrut_dep,rc_month,snht_maxlen,snht_min_sampsize,miss_val):

    tmean = np.zeros((fg_dep.shape[2]),dtype=fg_dep.dtype)
    tsquare = np.zeros((fg_dep.shape[2]),dtype=fg_dep.dtype)
    count = np.zeros(shape=(fg_dep.shape[2]),dtype=np.int32)
    tsas=np.zeros(fg_dep.shape,dtype=fg_dep.dtype)
    tsarad=np.zeros((1,6,fg_dep.shape[2]),dtype=fg_dep.dtype)
    tsahad=np.zeros((1,2,fg_dep.shape[2]),dtype=fg_dep.dtype)
    
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
    for ih in range(2):
        #SNHT.numba_snhteqmov((temperature[1,ip,:]-temperature[0,ip,:]).flatten(), rc_month, tsarad[0,ip,:], 
                             #snht_maxlen,snht_maxmiss,miss_val,count,tmean,tsquare)
        SNHT.numba_snhtmov_njit(hadcrut_dep[ih,:], tsahad[0,ih,:], 
                                             parr,
                                             count,tmean,tsquare)
            
        #print(ip)
    
    return tsas,tsarad,tsahad


def break_analyze(finfo):
    
    breaklist=np.zeros(50,dtype=np.int32)
    #plt.plot(finfo['days']/365.25,np.sum(finfo['tsas'][0,:,:],axis=0))
    #plt.plot(finfo['days']/365.25,np.sum(finfo['tsas'][1,:,:],axis=0))
    #plt.plot(finfo['days']/365.25,np.sum(finfo['tsarad'][0,:,:],axis=0))
    
    totalo=np.concatenate((finfo['tsarad'][0,:,:],finfo['tsas'][0,:13,:],finfo['tsas'][1,:13,:],finfo['tsahad'][0,:,:]),axis=0)
    weight=np.array([2.0]*6+[1.0]*6+[0.8]*3+[0.5]*4+[2.0]+[1.0]*6+[0.8]*3+[0.5]*4+[2.0])
    weight/=np.mean(weight)
    total=np.zeros(totalo.shape[1])
    for i in range(len(weight)):
        total+=totalo[i,:]*weight[i]
    #totalo[:6,:]*=3
    #total=np.max(totalo,axis=0)
    #print('totalo',totalo.shape)
    if RC['plot']:
        plt.figure(figsize=(6,12))
        lmax=totalo.shape[0]
        for l in range(0,lmax,2):    
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
        tbl=total[breaklist[l]]
        total[lb:rb]=0.
        #for i in range(lb,rb):
            #total[i]-=tbl*(((RC['snht_maxlen']/2-np.abs(breaklist[l]-i))/(RC['snht_maxlen']/2))**2)
            #if total[i]<0:
                #total[i]=0
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
def bilin(hadmeddict,temperatures,rcy,rcm,days,lat,lon):
    ilat=int(np.floor((lat+90)/5))
    ilon=int(np.floor((180+lon)/5))
    if ilon>=hadmeddict['hadmed'].shape[2]:
        ilon-=hadmeddict['hadmed'].shape[2]
    t850=np.zeros((2,len(days)))
    t850.fill(np.nan)
    tm=0
    mold=(rcy[days[0]]-1900)*12+rcm[days[0]]
    mc=0
    iref=0
    
    hadmedsmooth=rmeanw(hadmeddict['hadmed'][(rcy[days[:]]-1900)*12+rcm[days[:]],ilat,ilon],30)
    
    for ih in range(2):
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
            
        t850[ih,:]=temperatures[ih,13,:]-hadmedsmooth #rmeanw(hadmed[(rcy[days[i]]-1900)*12+rcm[days[i]],ilat,ilon],30))
        
    return t850
    
    
#@ifray(ray.remote,RC['findfuture'])
def RAOB_findbreaks(hadmeddict, method,fnorig):

    tt=time.time()
    from utils import rmeanw
    finfo={}
    
    variable = 'temperature'
    
    fidx = RC['filesorig'].index(fnorig)
    fn = RC['outdir'] + RC['statids'][fidx] + '/feedbackmerged' +RC['statids'][fidx] + '.nc'
    cache = False
    if os.path.isfile(fn):
        if os.path.getmtime(fn) > os.path.getmtime(fnorig):
            cache = True
            
    if not cache:
    
        with h5py.File(fnorig) as iofile:
            if '126' not in iofile['recordindices'].keys():
                    return None
            if np.unique(iofile['recordindices']['126'][:]).shape[0] < 365:
                    return None
            
            try:
                
                x = iofile['station_configuration/station_name'][:]
                sn = x.view('S'+str(x.shape[1]))[0][0].decode()
            except:
                sn = 'missing'

        with eua.CDMDataset(fnorig) as iofile:
            dq=iofile.read_data_to_3dcube(['temperature'], 
                    dates=None,
                    plevs=(RC['stdplevs']*100.).astype(int),
                    feedback=['fg_depar@body', 'biascorr@body'],
                    feedback_group='era5fb')
        #data = iofile.read_data_to_cube(variable,
                #dates=None,
                #plevs=(RC['stdplevs'][RC['pidx']]*100.).astype(int),
                #feedback=['fg_depar', 'bias'],
                #feedback_group='era5fb')#,**kwargs)

            if not dq:
                print(fnorig.split('/')[-1], 'failed', iofile['recordindices']['126'].shape[0])
                return None
        dql={}
        for k in dq.keys():
        ###dq[k].rename_dims({'plev':'pressure'})
            dql[RC['cdict'][k]]=dq[k].rename(RC['cdict'][k])
        xrdq=xr.Dataset(dql).rename_dims({'press':'pressure', 'datum': 'time'})
        xrdq.attrs['unique_source_identifier']=RC['wigos'][fidx] #fnf
        xrdq.attrs['station_name'] = sn
        for ih in range(xrdq['era5_fgdep'].shape[0]):
            for ip in range(xrdq['era5_fgdep'].shape[1]):
                v=xrdq['era5_fgdep'].values[ih,ip,:]
                o=xrdq['temperatures'].values[ih,ip,:]
                qs=numpy.nanquantile(v,[0.005,0.995])
                idx=numpy.where(numpy.logical_or(v<qs[0],v>qs[1]))
                #print(qs,idx)
                v[idx]=numpy.nan
                o[idx]=numpy.nan
                hilf=xrdq['bias_estimate'].values[ih,ip,:]
                hilf[numpy.isnan(hilf)]=0.
                xrdq['era5_fgdep'].values[ih,ip,:]=-v-hilf
                # ignore missing fgdep, bias_estimate
                #idx=numpy.where(numpy.logical_and(numpy.isnan(v),~numpy.isnan(o)))
                ##print(len(idx[0]))
                #xrdq['era5_fgdep'].values[ih,ip,idx]=0.


        try:
            os.mkdir(RC['outdir']+RC['statids'][fidx])
        except:
            pass
    
        fo = RC['outdir']+RC['statids'][fidx] + '/feedbackmerged' + RC['statids'][fidx] + '.nc'
        xrdq.to_netcdf(path=fo, format='NETCDF4_CLASSIC')
        print('wrote '+fo)

    try:
        
        #print('before read', time.time() - tt)
        with h5py.File(fn,'r') as f:
            if f['datum'].shape[0]<RC['snht_maxlen']:
                return finfo
            finfo['ifile'] = fn
            finfo['ifileorig'] = fnorig
            finfo['days']=f['datum'][:] - 1
            ndays = finfo['days'].shape[0]
            while  finfo['days'][ndays - 1] >= RC['months'].shape[0] and ndays > 0:
                ndays -= 1
            finfo['days'] = finfo['days'][:ndays]
            #print(fn.split('/')[-1],finfo['days'].shape)
            finfo['fg_dep']=-f['era5_fgdep'][:,np.array(RC['pidx']),:ndays]
            finfo['era5bc']=-f['bias_estimate'][:,np.array(RC['pidx']),:ndays] # ERA5 bias correction - to be used for adjusting most recent part of some series
            finfo['temperatures']=f['temperatures'][:,np.array(RC['pidx']),:ndays]
            #finfo['fg']=finfo['temperatures']-finfo['fg_dep']
            finfo['lon']=f['lon'][0]
            finfo['lat']=f['lat'][0]
            finfo['sid']=fn.split('/')[-1][:-3]
            finfo['station_name'] = f.attrs['station_name'].decode("utf-8")
            finfo['wigos']=f.attrs['unique_source_identifier'].decode("utf-8")
            if 'sonde_type' in f.keys():       
                #print(f['sonde_type'].shape)
                finfo['apriori_probs'],finfo['lastsonde']=metadata_analysis(f['sonde_type'][0,0,:ndays].view('S{}'.format(f['sonde_type'].shape[3])))
                #print(time.time()-tt)
            else:
                finfo['apriori_probs']=np.zeros_like(finfo['days'],dtype=np.float32)
                finfo['lastsonde']=0
        
        finfo['rharm'] =False
        if RC['rharm'] in ('rharmc', 'rcuon'):
            
            try:
                
                fnrharm_h = '/'.join(fn.split('/')[:-3]) + '/exp01/' + finfo['sid'][-6:] + '/rharm_h_' + finfo['sid'][-6:] +'.nc'
                #fnrharm = '/'.join(fn.split('/')[:-3]) + '/exp01/' + finfo['sid'][-6:] + '/rharm_' + finfo['sid'][-6:] +'.nc'
                with h5py.File(fnrharm_h,'r') as f:
                
                    radjust = f['ta'][:]
                    rdays = np.int64(f['datum'][0, :]) // 86400
                    finfo['rharmbc'] = np.empty_like(finfo['era5bc'])
                    finfo['rharmbc'].fill(np.nan)
                    l = 0
                    for i in range(finfo['days'].shape[0]):
                        while rdays[l] < finfo['days'][i]:
                            l += 1
                            if l == rdays.shape[0]:
                                l -= 1
                                break
                        if rdays[l] == finfo['days'][i]:
                            finfo['rharmbc'][:, :, i] =radjust[:, RC['pidx'], l] - finfo['temperatures'][:, :, i]
                            if RC['rharm'] == 'rharmc':
                                
                                finfo['temperatures'][:, :, i] += finfo['rharmbc'][:, :, i]
                                finfo['rharmbc'][:, :, i] = 0.
                finfo['rharm'] = True
            except FileNotFoundError:
                finfo['rharm'] = False
            
        
        bins=RC['months'][finfo['days']]
        finfo['anomalies']=np.empty_like(finfo['temperatures'])
        finfo['anomalies'].fill(np.nan)
        s=finfo['temperatures'].shape
        finfo['climatologies']=np.zeros_like(finfo['temperatures'],shape=(s[0],s[1],12))
        ccount=np.zeros(12)
        good=danomaly(finfo['temperatures'],bins,ccount,finfo['anomalies'],finfo['climatologies'])
        # finfo['hadcrut'] is difference between station anomalies at 850 hPa and (smoothed) hadcrut anomalies
        finfo['hadcrut_dep']=bilin(hadmeddict,finfo['anomalies'],RC['years'],RC['months'],finfo['days'],finfo['lat'],finfo['lon'])
        
        #print('before snht: {:6.4f}'.format(time.time()-tt)) 
        #print(typeof(finfo['fg_dep']),typeof(RC['months'][finfo['days']]),
              #typeof(RC['snht_maxlen']),typeof(RC['snht_maxmiss']),typeof(RC['miss_val']))
        if method=='SNHT':
            
            finfo['tsas'],finfo['tsarad'],finfo['tsahad']=do_test(finfo['fg_dep'],finfo['temperatures'],finfo['anomalies'],finfo['hadcrut_dep'],RC['months'][finfo['days']],
                                                  RC['snht_maxlen'],RC['snht_min_sampsize'],RC['miss_val'])  
            
            print('before analyze: {:6.4f}'.format(time.time()-tt))  
            #print(version_info)
            finfo['breaklist']=np.sort(break_analyze(finfo))
        elif method=='Binseg':
            testset=np.concatenate((finfo['fg_dep'][0,2:8,:],
                                    finfo['fg_dep'][1,2:8,:],finfo['fg_dep'][1,2:8,:]-finfo['fg_dep'][0,2:8,:],
                                    finfo['hadcrut'].reshape(1,finfo['hadcrut'].shape[0])),axis=0)
            algo = rpt.Binseg(model="l2").fit(testset.T)
            finfo['breaklist'] = algo.predict(pen=1.5*np.log(testset.shape[1]) * testset.shape[0] * np.nanvar(testset))
            finfo['breaklist'] = np.sort(finfo['breaklist'])[:-1]
            del algo
            del testset
        
        
        print(finfo['sid']+' findbreaks: found {:d}, {:6.4f}'.format(len(finfo['breaklist']), time.time()-tt))  
        obj_ref= ray.put(finfo)
        #print('after put', time.time()-tt)
    except MemoryError as e:
        print(fnorig, e)
        return None
    
    return obj_ref

ray_RAOB_findbreaks=ray.remote(RAOB_findbreaks)

def RAOB_adjustbs(finfo):
    
    sys.path.append(os.path.expanduser('../common/Rasotools/rasotools/'))
    sys.path.append(os.path.expanduser('../common/Rasotools/'))
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
    
    #sys.path.append(os.path.expanduser('../common/Rasotools/rasotools/'))
    #sys.path.append(os.path.expanduser('../common/Rasotools/'))
    from utils import tdist,extract,rmeanw
    from anomaly import ndnanmean2,nnanmean,ndnanvar2
    tt=time.time()
    if type(finfo) is not dict:
        finfo=ray.get(finfo)
    finfo['adjusted_fg_dep']=finfo['fg_dep'].copy()
    mask=np.isnan(finfo['adjusted_fg_dep'])
    s=finfo['adjusted_fg_dep'].shape
    tmean = np.zeros_like(finfo['adjusted_fg_dep'],shape=(s[2],12))
    tsquare = np.zeros_like(finfo['adjusted_fg_dep'],shape=(s[2],12))
    count = np.zeros(shape=(s[2],12),dtype=np.int32)
    tsa=np.zeros_like(finfo['adjusted_fg_dep'],shape=(s[2]))

    fb=finfo['breaklist']
    finfo['adjustments']=np.zeros_like(finfo['adjusted_fg_dep'],shape=(s[0],s[1],fb.shape[0]))
    finfo['adjusted_temperatures']=finfo['temperatures'].copy()
    finfo['adjusted_anomalies']=finfo['anomalies'].copy()
    finfo['adjusted_hadcrut_dep']=finfo['hadcrut_dep'].copy()
    
    sh=finfo['adjustments'].shape
    #RC_months=np.array(RC['months'])[finfo['days']]-1
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
        rcss=int(3*RC['snht_min_sampsize'])
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
        for ih in range(break_profiles.shape[1]):
            
            if np.sum(~np.isnan(finfo['adjusted_anomalies'][ih,13,istart:fb[ib]-30]))>RC['snht_min_sampsize'] and \
               np.sum(~np.isnan(finfo['adjusted_anomalies'][ih,13,fb[ib]+30:istop]))>RC['snht_min_sampsize']:
                y=-np.nanmean(finfo['adjusted_hadcrut_dep'][ih,istart:fb[ib]-30])+\
                    np.nanmean(finfo['adjusted_hadcrut_dep'][ih,fb[ib]+30:istop])
                break_profiles[ib,ih,13]=y
            else:
                break_profiles[ib,ih,13]=np.nan
                

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
                    
                if ~np.isnan(break_profiles[ib,ih,ip]) and (finfo['rharm'] == False or RC['rharm'] != 'rharmc'):
                    finfo['adjustments'][ih,ip,:ib]-=break_profiles[ib,ih,ip]
                    finfo['adjusted_fg_dep'][ih,ip,:fb[ib]]+=break_profiles[ib,ih,ip]
                    finfo['adjusted_temperatures'][ih,ip,:fb[ib]]+=break_profiles[ib,ih,ip]
                    finfo['adjusted_anomalies'][ih,ip,:fb[ib]]+=break_profiles[ib,ih,ip]
                    if ip==13:
                        
                        finfo['adjusted_hadcrut_dep'][ih,:fb[ib]]+=break_profiles[ib,ih,ip]
        
        
        #print(ib,fb[ib],break_profiles[ib,:,:])
    finfo['goodbreaklist']=list(fb[np.array(goodbreaklist,dtype=int)]) 
    
    
    
    print('adjust: {:d}, {:6.4f}'.format(len(finfo['goodbreaklist']), time.time()-tt))  
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
def RAOB_adjustbreaks(method,dists,active,finfo,i):
    
    tt=time.time()
    adjustments=copy.copy(finfo)
    adjustments['sdists']=extract(dists,i)
    adjustments['active']=active

    adjustments.update(RAOB_adjust(finfo))
    
    print(time.time()-tt)
    return ray.put(adjustments)

ray_RAOB_adjustbreaks=ray.remote(RAOB_adjustbreaks)

    
#@njit
def do_monmean(vars,tidx,fi):
    
    idx=np.searchsorted(fi['days'],tidx)
    mdays=idx[1:]-idx[:-1]
    midx=np.where(mdays>0)[0]
    
    fi['months']=tidx[midx]
            
    for v in vars:
        out=np.empty((fi[v].shape[0],fi[v].shape[1],fi['months'].shape[0]),dtype=fi[v].dtype)
        fi['m'+v]=pnm.monmean(fi[v],idx,out,RC['mon_thresh'])
        
    out=np.empty((fi[v].shape[0],fi[v].shape[1],fi['months'].shape[0]),dtype=np.int32)
    fi['goodmon']=pnm.goodmon(fi['temperatures'],idx,out)
            
                
    return

def write_monmean(opath,prefix,fi):
        
    if not os.path.isdir(os.path.dirname(opath+'/'+fi['sid'][-6:])):
        os.mkdir(os.path.dirname(opath+'/'+fi['sid'][-6:]))
    
    mtnames = {'feedbackglobbincorrmon': 'monthly adjusted temperature','feedbackglobbgmon': 'monthly background temperature',}
    
    fno=opath+'/'+fi['sid'][-6:]+'/'+prefix+fi['sid'][-6:]+'.nc'
    try:
        with netCDF4.Dataset(fno, "w") as dst:
            # copy global attributes all at once via dictionary
            new_globals={'Conventions': 'CF-1.1',
                         'title': 'Monthly radiosonde temperatures and -adjustments',
                         'institution': 'Institute for Meteorology and Geophysics, University of Vienna',
                         'Stationname': fi['station_name'], 
                         'history': datetime.datetime.today().strftime("%m/%d/%Y"), 
                         'source':'ERA5, IGRA2, NCAR UADB', 
                         'references':'Copernicus Early Upper Air Dataset', 
                         'url': 'early-upper-air.copernicus-climate.eu'
                         }
            dst.setncatts(new_globals)
            # copy dimensions
            dims = {'station': 1,'numdat': 4,'time': 1,'pressure': RC['stdplevs'].shape[0],'hour': fi['temperatures'].shape[0]}
            for name, dimension in dims.items(): #src.dimensions.items():
                if name =='time':
                    dst.createDimension(
                        name, (fi['mtemperatures'].shape[2] ))
                else:
                    dst.createDimension(
                        name, (dimension))
            # copy all file data except for the excluded
            dt = {'lat': np.dtype('float32'), 'lon': np.dtype('float32'), 'press': np.dtype('float32'), 'datum': np.dtype('int32'),
                  'montemp': np.dtype('float32'), 'goodmon': np.dtype('float32'), 'rasocorrmon': np.dtype('float32')}
            dd = {'lat': ('station',), 'lon': ('station',), 'press': ('pressure',), 'datum': ('numdat', 'time'),
                  'montemp': ('hour', 'pressure', 'time'), 'goodmon': ('hour', 'pressure', 'time'), 'rasocorrmon': ('hour', 'pressure', 'time')}

            d = {'lat': {'long_name': 'station latitude', 'units': 'degrees_north', 'axis': 'Y', 'valid_range': np.array([-90.,  90.]), 'missing_value': -999.0},
                 'lon': {'long_name': 'station longitude', 'units': 'degrees_east', 'axis': 'X', 'valid_range': np.array([-180.,  180.]), 'missing_value': -999.0},
                 'press': {'long_name': 'pressure levels', 'units': 'hPa', 'axis': 'Z', 'valid_range': np.array([   0., 1100.]), 'missing_value': -999.0},
                 'datum': {'long_name': 'datum', 'units': 'days since 1900-01-01 0:0:0', 'axis': 'T', 'calendar': 'gregorian', 'missing_value': -999.0},
                 'montemp': {'long_name': mtnames[prefix], 'units': 'K', 'missing_value': -999.0, 'valid_range': np.array([  0., 400.], dtype=np.float32), 'cell_methods': 'time: mean over months'},
                 'goodmon': {'long_name': 'number_of_values_in_month', 'missing_value': -999.0, 'valid_range': np.array([ 0., 31.], dtype=np.float32)},
                 'rasocorrmon': {'long_name': 'monthly_raso_correction', 'units': 'K', 'missing_value': -999.0, 'valid_range': np.array([-20.,  20.], dtype=np.float32), 'cell_methods': 'time: mean over months'}
                 }
            fi['press'] = RC['stdplevs']

            for name, atts in d.items():
                if name == 'datum':
                    x = dst.createVariable(name, dt[name], dd[name]) # type and dimensions
                    dst[name][:] = fi['months'][:]
                    # copy variable attributes all at once via dictionary
                    dst[name].setncatts(atts)
                elif name in ['montemp','goodmon','rasocorrmon']: #,'goodmon','rasocorrmon','eracorrmon']:
                    x = dst.createVariable(name, dt[name], dd[name])
                    if dt[name] in [np.dtype('float32'),np.dtype('float64')]:
                        
                        dst[name][:]=np.nan
                    else:
                        dst[name][:]=0
                    
                    if name == 'goodmon':
                        dst[name][:,:RC['pidx'].shape[0],:] = fi[RC['transdict'][name]][:]
                        
                    if name in ['montemp']:
                        if prefix == 'feedbackglobbincorrmon':
                            
                            dst[name][:,:RC['pidx'].shape[0],:] = fi[RC['transdict'][name]][:]
                        else:
                            dst[name][:,:RC['pidx'].shape[0],:] = fi['mfg']

                    elif name=='rasocorrmon':
                        dst[name][:,:RC['pidx'].shape[0],:] = fi['mtemperatures']-fi['mini_adjusted_temperatures']

                       
                    # copy variable attributes all at once via dictionary
                    dst[name].setncatts(atts)
                    
                else: # lat,lon,pressure levels
                    x = dst.createVariable(name,dt[name], dd[name] )
                    if fi[name] is np.ndarray:
                        
                        dst[name][:] = fi[name][:]
                    else:
                        dst[name][:] = np.array((fi[name],) )
                        
                    # copy variable attributes all at once via dictionary
                    dst[name].setncatts(atts)
    except Exception as e:
        print('could not write',fno, e)
    
    return

def write_adjustment(opath,prefix,fi,rich_ref0=None,rich_ref1=None,initial_adjust_RAOB=None,initial_adjust_RICH=None):
    
    if len(fi['goodbreaklist'])==0:
        print(fi['sid']+': no breaks for this station')
        return
    if 'ri' in prefix and not rich_ref0 and not rich_ref1:
        return
    
    if not os.path.isdir(os.path.dirname(opath)):
        os.mkdir(os.path.dirname(opath))
    
    fno=opath+'/'+fi['sid'][-6:]+'/'+prefix+fi['sid'][-6:]+'.nc'
    #fni=ipath+'/'+prefix+fi['sid'][-6:]+'.nc'
    try:
        with netCDF4.Dataset(fno, "w") as dst:
            # copy global attributes all at once via dictionary
            new_globals={'Conventions': 'CF-1.1',
                         'title': 'Monthly radiosonde temperatures and -adjustments',
                         'institution': 'Institute for Meteorology and Geophysics, University of Vienna',
                         'Stationname': fi['station_name'], 
                         'history': datetime.datetime.today().strftime("%m/%d/%Y"), 
                         'source':'ERA5, IGRA2, NCAR UADB', 
                         'references':'Copernicus Early Upper Air Dataset', 
                         'url': 'early-upper-air.copernicus-climate.eu'
                         }
            dst.setncatts(new_globals)
            # copy dimensions
            dims = {'station': 1,'numdat': 4,'time': 1,'pressure': RC['stdplevs'].shape[0],'hour': fi['temperatures'].shape[0]}
            for name, dimension in dims.items(): #src.dimensions.items():
                if name =='time':
                    dst.createDimension(
                        name, (len(fi['goodbreaklist'])+2 ))
                else:
                    dst.createDimension(
                        name, (dimension))
            # copy all file data except for the excluded
            dt = {'lat': np.dtype('float32'), 'lon': np.dtype('float32'), 'press': np.dtype('float32'), 'datum': np.dtype('int32'),
                  'rasocorr': np.dtype('float32'), 'rasobreak': np.dtype('float32')}
            dd = {'lat': ('station',), 'lon': ('station',), 'press': ('pressure',), 'datum': ('numdat', 'time'),
                  'rasocorr': ('hour', 'pressure', 'time'),'rasobreak': ('hour', 'pressure', 'time') }

            d = {'lat': {'long_name': 'station latitude', 'units': 'degrees_north', 'axis': 'Y', 'valid_range': np.array([-90.,  90.]), 'missing_value': -999.0},
                 'lon': {'long_name': 'station longitude', 'units': 'degrees_east', 'axis': 'X', 'valid_range': np.array([-180.,  180.]), 'missing_value': -999.0},
                 'press': {'long_name': 'pressure levels', 'units': 'hPa', 'axis': 'Z', 'valid_range': np.array([   0., 1100.]), 'missing_value': -999.0},
                 'datum': {'long_name': 'datum', 'units': 'days since 1900-01-01 0:0:0', 'axis': 'T', 'calendar': 'gregorian', 'missing_value': -999.0},
                 'rasocorr': {'long_name': 'raso_correct', 'units': 'K', 'missing_value': -999.0, 'valid_range': np.array([-20.,  20.], dtype=np.float32)}, 
                 'rasobreak': {'long_name': 'raso_correct', 'units': 'K', 'missing_value': -999.0, 'valid_range': np.array([-20.,  20.], dtype=np.float32)}, 
                 }
            fi['press'] = RC['stdplevs']

            # copy all file data except for the excluded
            for name, attrs in d.items():
                if name == 'datum':
                    x = dst.createVariable(name, dt[name], dd[name])
                    dst[name][:] = np.concatenate((np.array([1]),fi['days'][fi['goodbreaklist']],np.array([44998])),axis=0)
                    # copy variable attributes all at once via dictionary
                    dst[name].setncatts(attrs)
                elif name in ['rasocorr','rasobreak']: #,'goodmon','rasocorrmon','eracorrmon']:
                    
                    if rich_ref0 and 'ri' in prefix:
                        suff=''
                        if rich_ref1:
                            suff='0'
                        x = dst.createVariable(name+suff, dt[name], dd[name])
                        if dt[name] in [np.dtype('float32'),np.dtype('float64')]:
                            
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
                        x = dst.createVariable(name, dt[name], dd[name])
                        if dt[name] in [np.dtype('float32'),np.dtype('float64')]:
                            
                            dst[name][:]=np.nan
                        else:
                                dst[name][:]=0
                        if 'rio' in prefix and rich_ref1:
                            ad=rich_ref1['obs']['rich_adjustments'].shape
                            #dst[name][:,:RC['pidx'].shape[0],:] = np.concatenate((rich_ref1['obs']['rich_adjustments'][:,:,0:1],
                            if initial_adjust_RICH:
                                if 'initial_adjustments' in initial_adjust_RICH.keys():
                                    
                                    dst[name][:,:RC['pidx'].shape[0],:] = pnm.addini(rich_ref1['obs']['rich_adjustments'],initial_adjust_RICH['initial_adjustments'])
                                else:
                                    dst[name][:,:RC['pidx'].shape[0],:] = pnm.addini(rich_ref1['obs']['rich_adjustments'],np.zeros(fi['adjustments'].shape[:2]))
                            else:
                                dst[name][:,:RC['pidx'].shape[0],:] = pnm.addini(rich_ref1['obs']['rich_adjustments'],np.zeros(fi['adjustments'].shape[:2]))
                            #rich_ref1['obs']['rich_adjustments'][:],np.zeros((ad[0],ad[1],1))),axis=2)
                        elif 'rit' in prefix and rich_ref1:
                            ad=rich_ref1['tau']['rich_adjustments'].shape
                            #dst[name][:,:RC['pidx'].shape[0],:] = np.concatenate((rich_ref1['tau']['rich_adjustments'][:,:,0:1],
                                                                                  #rich_ref1['tau']['rich_adjustments'][:],np.zeros((ad[0],ad[1],1))),axis=2)
                            if initial_adjust_RICH:
                                if 'initial_adjustments' in initial_adjust_RICH.keys():
                                    dst[name][:,:RC['pidx'].shape[0],:] = pnm.addini(rich_ref1['tau']['rich_adjustments'],initial_adjust_RICH['initial_adjustments'])
                                else:
                                    dst[name][:,:RC['pidx'].shape[0],:] = pnm.addini(rich_ref1['tau']['rich_adjustments'],np.zeros(fi['adjustments'].shape[:2]))
                            else:
                                dst[name][:,:RC['pidx'].shape[0],:] = pnm.addini(rich_ref1['tau']['rich_adjustments'],np.zeros(fi['adjustments'].shape[:2]))
                    else:
                        x = dst.createVariable(name, dt[name], dd[name])
                        if dt[name] in [np.dtype('float32'),np.dtype('float64')]:
                            
                            dst[name][:]=np.nan
                        else:
                                dst[name][:]=0
                        
                        #print(fi['sid'][-6:],initial_adjust_RAOB.keys(),fi['adjustments'].shape)
                        if initial_adjust_RAOB:
                            if 'initial_adjustments' in initial_adjust_RAOB.keys():
                                
                                dst[name][:,:RC['pidx'].shape[0],:] = pnm.addini(fi['adjustments'],initial_adjust_RAOB['initial_adjustments'])
                            else:
                                dst[name][:,:RC['pidx'].shape[0],:] = pnm.addini(fi['adjustments'],np.zeros(fi['adjustments'].shape[:2]))
                        else:
                            dst[name][:,:RC['pidx'].shape[0],:] = pnm.addini(fi['adjustments'],np.zeros(fi['adjustments'].shape[:2]))
                            
                    # copy variable attributes all at once via dictionary
                    dst[name].setncatts(attrs)
                    
                else:
                    x = dst.createVariable(name, dt[name], dd[name])
                    if fi[name] is np.ndarray:
                        
                        dst[name][:] = fi[name][:]
                    else:
                        dst[name][:] = np.array((fi[name],) )
                    # copy variable attributes all at once via dictionary
                    dst[name].setncatts(attrs)
    except Exception as e:
        print(fno,'could not written',fno, e)
    
    return


#@ifray(ray.remote,RC['richfuture'])
def save_monmean(l,exper='exp02',fi={},rich_ref0=None,rich_ref1=None,initial_adjust_RAOB=None,initial_adjust_rich=None):
    
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
    
    tidx=np.where(RC['days']==0)[0] + 1
    #ipath=fi['sid'][-6:].join(os.path.dirname(tfile).split(tfile[-9:-3]))
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
    write_monmean(RC['outdir'],prefix,fi)
    prefix='feedbackglobbgmon'
    write_monmean(RC['outdir'],prefix,fi) 
    prefix='feedbackglobbincorrsave'
    write_adjustment(RC['outdir'],prefix,fi,
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

                    write_adjustment(RC['outdir'],prefix,fi,
                                     rich_ref0=riref0,
                                     rich_ref1=riref1,
                                     initial_adjust_RAOB=initial_adjust_RAOB,
                                     initial_adjust_RICH=iar)

ray_save_monmean=ray.remote(save_monmean)    
    
    
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


def rich_calcadjprofile(ldiff,sdiff,ldays,sdays,tbin,sbin,llmaxlen,slmaxlen):

    sprof=np.zeros(sdiff.shape[:2],dtype=ldiff.dtype)
    lprofl=np.zeros(ldiff.shape[:2],dtype=ldiff.dtype)
    sprofl=np.zeros(ldiff.shape[:2],dtype=ldiff.dtype)
    lprofr=np.zeros(ldiff.shape[:2],dtype=ldiff.dtype)
    sprofr=np.zeros(ldiff.shape[:2],dtype=ldiff.dtype)

    sprof=pnm.lagmean(ldiff,sdiff,ldays,sdays,llmaxlen,slmaxlen,RC['snht_min_sampsize'],lprofl,sprofl,lprofr,sprofr)
    return sprof

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
    sys.path.append(os.path.expanduser('../common/Rasotools/rasotools/'))
    sys.path.append(os.path.expanduser('../common/Rasotools/'))
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
        #print('richens',richcount,len(richens))
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
                    if ref[l]['active'][si][-1]<tib+RC['snht_min_sampsize'] or ref[l]['active'][si][0]>tib-RC['snht_min_sampsize']:
                        #print('reference ends early')
                        continue
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
                            #sibm1=np.searchsorted(ref[si]['days'],ref[l]['goodbreaklist'][ib-1])
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
                                         'ref':pnm.make_adjusted_series(ref[si]['fg_dep'],ref[si]['adjustments'],ref[si]['goodbreaklist'])},
                                    'obs':{'test':richtemp['obs']['rich_adjusted_anomalies'],
                                            'ref':pnm.make_adjusted_series(ref[si]['anomalies'],ref[si]['adjustments'],ref[si]['goodbreaklist'])}}
                        else:
                            #!!!!! this is the full  setting, requires complete first RICH iteration
                            try:
                                
                                #comp={'tau':{'test':richtemp['tau']['rich_adjusted_fg_dep'],'ref':rich_ref[si]['tau']['rich_adjusted_fg_dep']},
                                           #'obs':{'test':richtemp['obs']['rich_adjusted_anomalies'],'ref':rich_ref[si]['obs']['rich_adjusted_anomalies']}}
                                comp={'tau':{'test':richtemp['tau']['rich_adjusted_fg_dep'],
                                             'ref':pnm.make_adjusted_series(ref[si]['fg_dep'],rich_ref[si][richcount]['tau']['rich_adjustments'],ref[si]['goodbreaklist'])},
                                           'obs':{'test':richtemp['obs']['rich_adjusted_anomalies'],
                                                  'ref':pnm.make_adjusted_series(ref[si]['anomalies'],rich_ref[si][richcount]['obs']['rich_adjustments'],ref[si]['goodbreaklist'])}}
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
                            if ib==4 and ref[si]['sid'][-6:]=='045004':
                                print('bad buddy')
                            if method=='tau':   
                                rich['found']['ref'][ib].append(ref[si]['sid'][-6:])
                                rich['found']['refsid'][ib].append('{:5.0f} km'.format(ref[l]['sdists'][si]*6370.))
                                    
                    lblen=len(rich['tau']['buddies'][ib])        
                    if lblen>=min_neighbours[rich_iter]:
                        for method in 'tau','obs':
                            try:
                                if len(rich[method]['buddies'][ib])>0:
                                    buddies=np.array(rich[method]['buddies'][ib])
                                    rich[method]['break_profiles'][ib,:,:]=pnm.calc_breakprofiles(buddies,
                                                                                          rich['found']['weight'][:buddies.shape[0]])
                                else:
                                    rich[method]['break_profiles'][ib,:,:]=np.nan
                            except MemoryError as e:
                                print(e,lblen,ib,rich[method]['break_profiles'].shape,rich['found']['weight'].shape,
                                np.array(rich[method]['buddies'][ib]).shape)
                                rich[method]['break_profiles'][ib,:,:]=np.nan
                                
                            if ref[l]['sid'][-6:]=='159431' and ib==4 and method=='obs':
                                plt.plot(ref[l]['days'][lidx[0]:lidx[2]],comp['tau']['test'][:,:,lidx[0]:lidx[2]][0,3,:])
                                plt.plot(ref[si]['days'][sidx[0]:sidx[2]],comp['tau']['ref'][:,:,sidx[0]:sidx[2]][0,3,:])
                                plt.plot(ref[si]['days'][sidx[1]],0,'ro')
                                plt.plot(ref[l]['days'][lidx[1]],0,'go')
                                
                                plt.plot(ref[l]['days'][lidx[0]:lidx[2]],comp['obs']['test'][:,:,lidx[0]:lidx[2]][0,3,:])
                                plt.plot(ref[si]['days'][sidx[0]:sidx[2]],comp['obs']['ref'][:,:,sidx[0]:sidx[2]][0,3,:])
        
                                plt.savefig('rich_{}_{}_{}_ts.png'.format(richcount,ref[l]['sid'][-6:],ib))
                                plt.close()
                                taubuddies=np.array(rich['tau']['buddies'][ib])
                                obsbuddies=np.array(rich['obs']['buddies'][ib])
                                for m in 'tau','obs':
                                    
                                    for i in range(taubuddies.shape[0]):
                                    
                                        plt.semilogy(rich[m]['buddies'][ib][i][0,:],RC['stdplevs'][RC['pidx']],'b',label=rich['found']['ref'][ib][i]+'00')
                                        plt.semilogy(rich[m]['buddies'][ib][i][1,:],RC['stdplevs'][RC['pidx']],'r',label=rich['found']['ref'][ib][i]+'12')
                                    plt.ylim(1000,10)
                                    #plt.legend()
                                    plt.savefig('rich'+m+'_{}_{}_{}_prof.png'.format(richcount,ref[l]['sid'][-6:],ib))
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

def initial_adjust(l,obj_ref,obj_rich_ref=None,initial_composite_ref=None,rich_initial_composite_ref=None):
    sys.path.append(os.path.expanduser('../common/Rasotools/rasotools/'))
    sys.path.append(os.path.expanduser('../common/Rasotools/'))
    
    tt=time.time()

    if type(obj_ref) is list: # if called by ray, the first dereferencing is done automatically
        # instead of deferencing the whole list of object references at once, we dereference only those needed -
        # Note that for each stations only at most 500 buddies (out of 1300) are tested.
        ref=[None for i in range(len(obj_ref))]
        ref[l]=ray.get(obj_ref[l])
        oref=obj_ref
        #ref=ray.get(obj_ref)
        if obj_rich_ref:
            orref=obj_rich_ref
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
        oref=ray.get(obj_ref) # first dereference - gets list ob object references
        ref=[None for i in range(len(oref))]
        ref[l]=ray.get(oref[l]) # second dereference - gets actual dictionaries
        if obj_rich_ref:
            orref=ray.get(obj_rich_ref) # first dereference - gets list ob object references
            rich_ref=[None for i in range(len(orref))]
            rich_ref[l]=ray.get(orref[l]) # first dereference - gets list ob object references
            
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

                        
    elif RC['initial_adjustments']=='era5neighbours':
    #if True: 
        if ref[l]['lastsonde'] or ref[l]['rharm']: # and RC['initial_adjustments'] not in ['era5','era5bc']:
            #start,igood=findstart(ref[l]['fg_dep'],RC['mean_maxlen'])
            #ini=np.zeros(ref[l]['adjusted_fg_dep'].shape[:2])
            #RAOBini['initial_adjustments']= ndnanmean2(ref[l]['adjusted_fg_dep'][:,:,start:],ini,RC['snht_min_sampsize'])
            #RAOBini['initial_adjustments'][np.isnan(RAOBini['initial_adjustments'])]=0.
            RAOBini['initial_adjustments']=np.zeros_like(ref[l]['fg_dep'],shape=ref[l]['fg_dep'].shape[:2])
            if ref[l]['rharm']:
                start=np.max((0,ref[l]['days'].shape[0]-RC['mean_maxlen']))
                if ref[l]['days'].shape[0]-start>RC['snht_min_sampsize']:
                    
                    ini=np.zeros(ref[l]['adjusted_fg_dep'].shape[:2])
                    rharmini = np.zeros_like(ini)
                    rawini = np.zeros_like(ini)
                    rharmini = ndnanmean2(ref[l]['rharmbc'][:,:,start:],ini,RC['snht_min_sampsize'])[:]
                    rawini = ndnanmean2(ref[l]['adjusted_fg_dep'][:,:,start:],ini.copy(),RC['snht_min_sampsize'])[:]
                    RAOBini['initial_adjustments']= rawini - rharmini
                    RAOBini['initial_adjustments'][np.isnan(RAOBini['initial_adjustments'])]=0.
                
            if obj_rich_ref:
                for i in range(32):
                    for method in 'tau','obs':
                        richensini[i][method]['initial_adjustments']=RAOBini['initial_adjustments'].copy()
            pass
        
        else:
            #start=np.max((0,ref[l]['days'].shape[0]-RC['mean_maxlen']))
            start,igood=pnm.findstart(ref[l]['fg_dep'],RC['mean_maxlen'])
            spagini=[]
            if any(igood>RC['snht_min_sampsize']):
                refcount=0
                for si in sdidx[1:]: # try ref stations beginning with nearest:
                    if si>=len(ref):
                        continue
                    if ref[l]['active'][si][-1]<ref[l]['days'][-1]-RC['snht_min_sampsize']:
                        #print('reference ends early')
                        continue
                        
                    refcount +=1
                    if ref[l]['sdists'][si]*6370. > RC['weight_distance'][0] and refcount>500: # second condition is for test stations at remote places
                        print('not enough buddies found','{:5.1f}'.format(ref[l]['sdists'][si]*6370.),refcount)
                        break

                    if ref[si] is None:
                        ref[si]=ray.get(oref[si])
                        #print(ref[si]['sid'[-6:]],ref[si]['days'][-1],'reference ok')

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
                        if ref[l]['active'][si][-1]<ref[l]['days'][-1]-RC['snht_min_sampsize']:
                            #print('reference ends early')
                            continue
                        if ref[l]['sdists'][si]*6370. > RC['weight_distance'][0] and refcount>500: # second condition is for test stations at remote places
                            print('not enough buddies found','{:5.1f}'.format(ref[l]['sdists'][si]*6370.),refcount)
                            break
        
                        if rich_ref[si] is None:
                            rich_ref[si]=ray.get(orref[si])
                        if ref[si] is None:
                            ref[si]=ray.get(oref[si])
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
                                                        ri_test_adjusted_fg_dep[iens][method]=pnm.make_adjusted_series(ref[l]['fg_dep'],rich_ref[l][iens][method]['rich_adjustments'],
                                                                            ref[l]['goodbreaklist'])
                                                        if iens==len(richensini)-1 and method=='obs':
                                                            firstbuddy=False

                                                    # calculate rich adjusted reference time series
                                                    ini=np.zeros(ref[si]['adjusted_fg_dep'].shape[:2])
                                                    ri_adjusted_fg_dep[iens][method]=pnm.make_adjusted_series(ref[si]['fg_dep'],rich_ref[si][iens][method]['rich_adjustments'],
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
            ini2 = np.zeros_like(ini)
            rharmini = np.zeros_like(ini)
            rawini = np.zeros_like(ini)
            rharmini = ndnanmean2(ref[l]['rharmbc'][:,:,start:],ini,RC['snht_min_sampsize'])[:]
            rawini = ndnanmean2(ref[l]['adjusted_fg_dep'][:,:,start:],ini.copy(),RC['snht_min_sampsize'])[:]
            RAOBini['initial_adjustments']= rawini - rharmini
            RAOBini['initial_adjustments'][np.isnan(RAOBini['initial_adjustments'])]=0.

    RAOBini['initial_adjustment_method']=RC['initial_adjustments']
    print(ref[l]['sid'][-6:],time.time()-tt)
    if obj_rich_ref:
        return ray.put((RAOBini,richensini))
    else:
        return ray.put(RAOBini)
        

ray_initial_adjust=ray.remote(initial_adjust)

def do_copy(fni, fno, grvdict={'recordindices': [], 'observations_table': ['index', 'date_time', 'z_coordinate', 'observed_variable', 'observation_value']}, mode='r+'):
    
    if(not os.path.isfile(fno)):
        mode = 'w'
        
    with h5py.File(fni,'r') as f:
        try:
            
            with h5py.File(fno,mode) as g:
                
                for k, vv in grvdict.items():
                           
                    if k not in g.keys():
                        
                        g.create_group(k)
                    
                    if not vv:
                        vv = list(f[k].keys())
                        
                    for v in vv:
                        
                        try:
                            try:
                                del g[k][v]
                            except Exception as e:
                                #print(e)
                                pass
                            
                            g[k].create_dataset_like(v,f[k][v],compression='gzip')
                            g[k][v][:]=f[k][v][:]
                        except Exception as e:
                            print(e)
                            continue
                        #ll+=1
                        for a,av in f[k][v].attrs.items():
                            if a not in ('CLASS','NAME','REFERENCE_LIST','DIMENSION_LIST'):
                                #print(a,av)
                                g[k][v].attrs[a]=av
                    
                    for v in g[k].keys(): #var_selection:
                        l=0            
                        try:
                            fvv=g[k][v]
                            if 'string' not in v and v!='index':                    
                                g[k][v].dims[l].attach_scale(g[k]['index'])
                                #print(v,fvv.ndim,type(fvv[0]))
                                if fvv.ndim==2 or type(fvv[0]) in [str,bytes,numpy.bytes_]:
                                    slen=sdict[v]
                                    #slen=10
                                    g[k][v].dims[1].attach_scale(g[k]['string{}'.format(slen)])
                        except Exception as e:
                            print(fn.split('/')[-1],e)
                            pass
                
        except Exception as e:
            print(fn.split('/')[-1],e)
            pass
        
        print(fn.split('/')[-1]+' copied to '+fno)

def add_adj(fi, mode='r'):
    
    statid=fi['sid'][-6:] #adjfile[-8:-3]
    adjfile = RC['outdir'] + statid + '/feedbackglobbincorrsave' + statid + '.nc'
    adjustments={'raobcore':os.path.expandvars(adjfile),
                 'rich':os.path.expandvars('corrsave_rio24_'.join(adjfile.split('corrsave'))),
                 'rase':os.path.expandvars(('ERA5bc_RAOBCORE_v'+RC['version']+'_').join(adjfile.split('feedbackglobbincorrsave'))),
                 'rise':os.path.expandvars(('ERA5bc_RAOBCORE_v'+RC['version']+'_').join(adjfile.split('feedbackglobbincorrsave')))}
    adjname={'raobcore':'rasocorr',
                 'rich':'rasocorr',
                 'rase':'bias',
                 'rise':'richbias'}

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
    

            for k,v in adjustments.items():
                try:
                    
                    adjustments=xr.open_dataset(v,decode_times=False)
                    if adjustments.datum.ndim==2:
                        atime0=(adjustments.datum[0].values.astype(int)-1)*86400.
                    else:
                        atime0=(adjustments.datum.values.astype(int)-1)*86400.
                   
                    
                    mask=adjustments[adjname[k]].values==-999.
                    adjustments[adjname[k]].values[mask]=np.nan
                    tt=time.time()
                    adj=pnm.add_biasestimate(xyz.values,xyzt,xyz.z_coordinate.values,atime0,
                                         adjustments[adjname[k]].values,adjustments.press.values*100)
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
                        key = k.upper()+'_bias_estimate'
                        alldict = pandas.DataFrame({key:xyz})
                        write_dict_h5(outfile, alldict, 'advanced_homogenisation', {key: { 'compression': 'gzip' } }, [key])  
                        #data.write_observed_data(k.upper()+'_bias_estimate',
                                                 #ragged=xyz,  # input data
                                                 #varnum=eua.cdm_codes['temperature'],  # observed_variable to be aligned with
                                                 #group='advanced_homogenisation',   # name of the new group
                                                 #data_time='date_time',  # named datetime coordinate
                                                 #data_plevs='z_coordinate',  # named pressure coordinate
                                                 #attributes={'version':RC['version']}
                                                #)
                    print('write:',time.time()-tt)
                except MemoryError as e:
                    print('could not write', k,e)
    except MemoryError as e:
        print('could not write '+outfile, e)

ray_add_adj=ray.remote(add_adj)


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
    if len(sys.argv)>2:
        pattern=sys.argv[2]
        RCdir = sys.argv[1]
        pkl='RAOBx'
    elif len(sys.argv)>1:
        pattern=''
        RCdir = sys.argv[1]
        pkl='RAOBx'
    else:
        pattern='v1'
        RCdir = 'exp00'
        pkl='RAOBx'
    tt=time.time()
    process = psutil.Process(os.getpid())
    #files = glob.glob('/users/staff/leo/fastscratch/rise/1.0/'+RC['exp']+'/[01]*/feedbackmerged'+pattern+'*.nc')[:]
    RCfile = os.path.expandvars('/users/staff/leo/fastscratch/rise/1.0/')+RCdir+'/RC.json'

#    ray.init(num_cpus=RC['CPUs'], object_store_memory=1024*1024*1024*50)
    ray.init(num_cpus=25, object_store_memory=1024*1024*1024*50)

    RC = RC_ini(RCfile)
                                
    filesorig = glob.glob('/mnt/users/scratch/leo/scratch/'+RC['CUON']+'/*'+pattern+'*.nc')[:]
    filesorig.sort()
    lats = [];lons = []; wigos = []; wshort = [];statids = []
    goodfilesorig = []
    lold = 0
    llold = 0
    for fn in filesorig:
        if '00000' in fn:
            continue

        #lats.append(f['observations_table']['latitude'][-1])
        #lons.append(f['observations_table']['longitude'][-1])
        if 'orphan' not in fn:
            wigos.append(fn.split('_CEUAS')[0].split('/')[-1]) 
            wshort.append(fn.split('_CEUAS')[0].split('-')[-1][:5])
        else:
            wigos.append(fn.split('_CEUAS')[0].split('/')[-1]) 
            wshort.append('orpha')
            print(wigos[-1], wshort[-1])
                         
        if wshort[-1] == 'orpha':
            wshort[-1] = 'or'
            l = lold
            while '{:0>4}'.format(l) + wshort[-1] in statids:
                l += 1
            wshort[-1] = '{:0>4}'.format(l) + wshort[-1]
            lold = l
            statids.append(wshort[-1])
        else:
            if 'data' in wshort[-1]:
                print(wshort[-1], wigos[-1])
                
            while len(wshort[-1]) < 5:
                print(wshort[-1], wigos[-1])
                wshort[-1] = '0' + wshort[-1]
                print('')
            ll = llold
            while str(ll) + wshort[-1] in statids:
                ll += 1
            statids.append(str(ll)+wshort[-1])
            llold = ll
        goodfilesorig.append(fn)
        print(statids[-1])
    RC['filesorig'] = goodfilesorig
    RC['statids'] = statids
    RC['wigos'] = wigos
    RC['outdir'] = '/users/staff/leo/fastscratch/rise/1.0/'+RC['exp']+'/'
    RC['cdict'] ={'temperature':'temperatures','fg_depar@body':'era5_fgdep','biascorr@body':'bias_estimate','lat':'lat','lon':'lon','hours':'hours'}
    print(time.time() -tt)
    hadmeddict= read_hadCRUT5('/users/staff/leo/fastscratch/rise/1.0/common/','HadCRUT.5.0.1.0',RC['refdate'])
    hadmeddict_ref = ray.put(hadmeddict)

    
    #ray.init(address='131.130.157.11:49705', _redis_password='5241590000000000')

    
    if not RC['findfuture']:  
        func=partial(RAOB_findbreaks,hadmeddict, 'SNHT') # Binseg
        obj_ref=list(map(func,RC['filesorig']))
        obj_ref = [i for i in obj_ref if i]
        finfo=ray.get(obj_ref)
        #add_adj(finfo[0], mode='r')
    else:
        futures = []
        for fn in RC['filesorig']:
            futures.append(ray_RAOB_findbreaks.remote(hadmeddict_ref, 'SNHT', fn))
        obj_ref  = ray.get(futures)
        #futures = []
        #for o in obj_ref:
            #futures.append(ray_add_adj.remote(o, mode='r'))
        #ray.get(futures)
    
    for key in ('filesorig', 'statids', 'wigos'):
        l = 0
        m = 0
        for i in obj_ref:
            if i:
                RC[key][l] = RC[key][m]
                l += 1
            m += 1
        RC[key] = RC[key][:l]
        
    obj_ref = [i for i in obj_ref if i]
    
    if False and RC['write_to_backend']:
        
        futures = []
        for o in obj_ref:
            futures.append(ray_add_adj.remote(o, mode='r'))
        ray.get(futures)
        print('wrote to backend files',time.time()-tt)
    
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
    obj=ray.get(obj_ref)
    sids=[f['sid'] for f in obj] 
    sdist_ref=stdists(obj)
    active=[ob['days'][[0,-1]] for ob in obj]
    del obj
    active_ref=ray.put(active)
    print(time.time()-tt)

    
    futures=[]
    if not RC['findfuture']:  
        func=partial(RAOB_adjustbreaks,'SNHT',ray.get(sdist_ref), active,ray.get(obj_ref[0]))
        results=list(map(func,[0]))
        #results=list(map(func,obj_ref))
        obj_ref=results
        #sdist=ray.get(sdist_ref)
        #for l in range(len(obj_ref)):
            #futures.append(RAOB_adjustbreaks('SNHT',sdist,active,ray.get(obj_ref[l]),l))
        #obj_ref=ray.put(futures)
    else:
        for l in range(len(obj_ref)):
            futures.append(ray_RAOB_adjustbreaks.remote('SNHT',sdist_ref,active_ref,obj_ref[l],l))
        obj_ref = ray.get(futures)
    
    del futures
    
    #break_sum=np.sum([len(l['goodbreaklist']) for l in ray.get(obj_ref)])
    #print('Total number of "good" breaks',break_sum)

    
    if pattern=='':
        pattern='01001'
        
    lx= sids.index('feedbackmerged0'+pattern)
    
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
            if d['lastsonde']==1 or (d['rharm'] == 1 and RC['rharm'] in ('rharmc', 'rcuon')):
                lgood.append(l)
            elif d['days'][-1]>int(110*365.25):
                lrecent.append(l)
            else:
                lneedscomposite.append(l)
            
        ifile=ray.get(obj_ref[0])['ifile'] #'exp03'.join(files[0].split(RC['exp']))
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
    if(not RC['findfuture']):
        
        iadjustlist =[]    
        for l in lgood:
            iadjustlist.append(initial_adjust(l,single_obj_ref)) #,single_results_ref,sdist_ref)
    else:
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
    
    if lx >= len(obj_ref):
        lx = 0
    save_monmean(0,exper=RC['exp'],fi=[obj_ref[lx]],initial_adjust_RAOB=[initial_adjust_ref[lx]])
    print('before write',time.time()-tt)
    futures=[]
    for l in range(len(obj_ref)):
        futures.append(ray_save_monmean.remote(0,exper=RC['exp'],fi=[obj_ref[l]],initial_adjust_RAOB=[initial_adjust_ref[l]]))
    x = ray.get(futures)
    print('after write',time.time()-tt)
    
    if RC['add_solelev']:  # zum Testen
        print('before addsolelev',time.time()-tt)
        add_solelev(ray.get(obj_ref[lx]), RC)
        futures=[]
        for l in range(len(obj_ref)):
            futures.append(ray_add_solelev.remote(obj_ref[l], RC))
            #add_solelev(ray.get(obj_ref[l]), RC)
        x = ray.get(futures)
        add_adj(ray.get(obj_ref[0]), mode='r')
        print('after addsolelev',time.time()-tt)

    #if RC['write_to_backend']:
        
        #futures = []
        #for o in obj_ref:
            #futures.append(ray_add_adj.remote(o, mode='r'))
        #ray.get(futures)
        #print('wrote to backend files',time.time()-tt)

    #func=partial(save_monmean,files[0])
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
    sodr = ray.get(ray.get(single_rich_ref0))
    l = -1
    for s in sodr:
        l += 1
        sob = s[24]['obs']['break_profiles']
        for ib in range(len(sob)):
            for ih in range(2):
                if np.any(np.abs(sob[ib][ih]) >1.e4):
                    print(l,ib,ih, ray.get(obj_ref[l])['ifile'],sob[ib][ih] )
                    
    if True:
        rich_ref00=rich_adjust(lx,single_obj_ref,0)
        breakplot(lx,obj_ref,obj_rich_ref0=[rich_ref00])
        #save_monmean(ifile,0,exper=RC['exp'],fi=[obj_ref[lx]],rich_ref1=[rich_ref00])
    if True:
        rich_ref2=rich_adjust(lx,single_obj_ref,1,obj_rich_ref=single_rich_ref0)  
        breakplot(lx,obj_ref,obj_rich_ref0=[rich_ref0[lx]],obj_rich_ref1=[rich_ref2])
        #save_monmean(ifile,0,exper=RC['exp'],fi=[obj_ref[lx]],rich_ref0=[rich_ref0[lx]],rich_ref1=[rich_ref2])

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

    save_monmean(0,exper=RC['exp'],fi=[obj_ref[lx]],rich_ref1=[rich_ref1[lx]],
                 initial_adjust_RAOB=[initial_adjust_ref[lx]],
                 initial_adjust_rich=[rich_initial_adjust_ref[lx]])

    print('before write',time.time()-tt)
    futures=[]
    for l in range(len(obj_ref)):
        #futures.append(ray_save_monmean.remote(files[0],l,single_obj_ref,single_results_ref,rich_ref0=single_rich_ref0,rich_ref1=single_rich_ref1))
        futures.append(ray_save_monmean.remote(0,exper=RC['exp'],fi=[obj_ref[l]], rich_ref1=[rich_ref1[l]],
                                               initial_adjust_RAOB=[initial_adjust_ref[l]],
                                               initial_adjust_rich=[rich_initial_adjust_ref[l]]))
    
    res = ray.get(futures)
    if RC['add_solelev']:
        print('rich before addsolelev',time.time()-tt)
        add_solelev(ray.get(obj_ref[lx]), RC)
        futures=[]
        for l in range(len(obj_ref)):
            futures.append(ray_add_solelev.remote(obj_ref[l], RC))
        x = ray.get(futures)
        print('after addsolelev',time.time()-tt)
    
    if RC['write_to_backend']:
        
        futures = []
        for o in obj_ref:
            futures.append(ray_add_adj.remote(o, mode='r'))
        ray.get(futures)
        print('wrote to backend files',time.time()-tt)

        
    print('Mem [MiB]',process.memory_info().rss//1024//1024)
    print(time.time()-tt)
    ray.shutdown()
    print('end')
     
