#!/usr/bin/env
# coding: utf-8

from numba import njit
import numpy as np
import sys,glob
import zipfile, os, time
import urllib3
from datetime import datetime, timedelta
import glob, shutil
import hdf5plugin
import h5py
sys.path.append(os.getcwd()+'/../cds-backend/code/')
sys.path.append(os.getcwd()+'/../harvest/code/')
from harvest_convert_to_netCDF import write_dict_h5
from convert_numbamod import * 
import gc
import cds_eua3 as eua
#eua.logging_set_level(30)
#import xarray as xr
import cdsapi, zipfile, os, time
import copy
from shutil import copyfile
from multiprocessing import Pool
sys.path.insert(0,os.getcwd()+'/../resort/rasotools-master/')
import rasotools
dir(rasotools)
import warnings
from functools import partial
import matplotlib.pylab as plt
warnings.filterwarnings('ignore')
import pickle
import ray
from assemble_longrecord import h5concatenate
ray_h5concatenate = ray.remote(h5concatenate)
import pandas as pd

@ray.remote
class h5p:
    def __init__(self, fn):
        self.fn = fn
        self.fp = h5py.File(fn, 'r')
        
    
    def content(self):
        return list(self.fp.keys())
    
    def name(self):
        return self.fn
    
    def shape(self, v):
        return self.fp[v].shape
    
    def slice4d(self, v, tup):            
        return self.fp[v][tup[0], tup[1], tup[2], tup[3]]
    def slice3d(self, v, tup):            
        return self.fp[v][tup[0], tup[1], tup[2]]
    def slice2d(self, v, tup):            
        return self.fp[v][tup[0], tup[1]]
    def slice(self, v, tup):            
        return self.fp[v][tup]
    def full(self, v):            
        return self.fp[v][:]
    def vattr(self, v, attr):
        return self.fp[v].attrs[attr]
    def vattrs(self, v):
        d = {}
        for k, v in  self.fp[v].attrs.items():
            if type(v) is str or type(v) is int or type(v) is float:
                d[k] = v
                
        return d
        
    
    
#fn = '/mnt/users/scratch/leo/scratch/era5/gridded/era5t.196507.131.nc'
#ray.init(num_cpus=10)

#tt = time.time()
#x = h5p.remote(fn)
#print(time.time() - tt)
#for k in range(14):
    #y = ray.get(x.slice4d.remote('u', (slice(k, k+2), slice(k, k+2), slice(k, k+2), slice(k, k+2))))
    #print(time.time() - tt)

#del x
#print(ray.get(x.content.remote()))
#print(ray.get(x.shape.remote('u')))
#print(ray.get(x.slice4d.remote('u', (slice(0, 2), slice(0, 2), slice(0, 2), slice(0, 2)))))
#print(type(x))


@njit(cache=True, boundscheck=False)
def add_fb(loaded_obstab,loaded_feedback,ref20CR,refera5an,refera5fc):
    i20=0
    iera=0
    for i in range(loaded_obstab['date_time'].shape[0]):
        if loaded_feedback['fg_depar@body'][i]!=loaded_feedback['fg_depar@body'][i]:
            if loaded_obstab['observation_value'][i]==loaded_obstab['observation_value'][i]:
                if i < refera5fc.shape[0] and  i < ref20CR.shape[0]:
                    
                    if refera5fc[i]==refera5fc[i]:
                        loaded_feedback['fg_depar@body'][i]=loaded_obstab['observation_value'][i]-refera5fc[i]
                        loaded_feedback['an_depar@body'][i]=loaded_obstab['observation_value'][i]-refera5an[i]
                        loaded_feedback['biascorr@body'][i]=0.
                        iera+=1
                    elif ref20CR[i]==ref20CR[i]:
                        loaded_feedback['fg_depar@body'][i]=loaded_obstab['observation_value'][i]-ref20CR[i]
                        loaded_feedback['an_depar@body'][i]=loaded_obstab['observation_value'][i]-ref20CR[i]
                        loaded_feedback['biascorr@body'][i]=0.
                        i20+=1
        #if i%1000000==0:
            #print(i,i20,iera)

from numba.typed import List


@njit(cache=True, boundscheck=False)
def augment2(obstab, a_obstab, loaded_feedback, a_loaded_feedback,
            ri, ts, ps, humvar, wvar, zvar, fn):


    #recordindex=np.empty(ri.shape[0],obstab['date_time'].dtype)
    #recordtimestamp=np.empty(ts.shape[0],obstab['date_time'].dtype)
    rim = np.max(ri[1:]-ri[:-1])
    zlist = np.empty(rim, dtype=np.float32)
    ztype = np.empty(rim, dtype=np.int32)
    #allvar = np.concatenate(((0, 117, 126), humvar, wvar))
    allvar = np.zeros(11, dtype=np.int32)
    allvar[1:3] = np.array((117, 126))
    allvar[3:7] = humvar
    allvar[7:] = wvar
    alli = np.zeros(np.max(allvar) + 1, dtype=np.int32)
    a_obstab['observation_value'][:] = np.nan
    
    
    for i in range(len(allvar)):
        alli[allvar[i]] = i

    j=-1 # augmented index
    jsave=0
    p= obstab['z_coordinate'][0]-1 # z_coordinate[0]-1
    rts= obstab['date_time'][0] # date_time[0]

    j=-1
    porig = -1
    #good=True
    spurious_z = 0
    z_coordinate = obstab['z_coordinate']
    g = 9.80665
    oovs = obstab['observation_value'].shape[0]
    oidx = np.arange(oovs, dtype=np.int64)
    orig = np.empty(a_obstab.shape[0], dtype=np.int64)
    orig.fill(-1)
    ## adjust for relative humidity in %
    idx = np.where((obstab['observed_variable']==138)&(obstab['observation_value']>2.0))
    if len(idx) > 0:
        obstab['observation_value'][idx] /=100.
    idx = np.where((obstab['observed_variable']==138)&(obstab['observation_value']<0.))
    if len(idx) > 0:
        obstab['observation_value'][idx] = np.nan
       
    k = 0
    ak = 0
    while k < oovs:
        i = k
        j = 0
        plist = np.empty(zlist.shape[0]+16, dtype=np.float32)
        isave = np.empty(zlist.shape[0]+16, dtype=np.int32)
        isave[:] = -1
        while i < oovs and obstab['date_time'][i] ==rts:

            
            if z_coordinate[i]!=porig or i == k:
                porig = z_coordinate[i]
                zct = obstab['z_coordinate_type'][i]
                oov = obstab['observed_variable'][i]
                if zct == 1:
                    
                    plist[j]=z_coordinate[i]
                    ztype[j] = zct
                elif oov == 0:
                    plist[j]=z_coordinate[i]
                    ztype[j] = zct
                elif zct == 0: 
                    zlist[j] = z_coordinate[i]*g
                    plist[j] = z_to_p_ifs(zlist[j])
                    z_coordinate[i] = plist[j]
                    #if obstab['observed_variable'][i] == 126:
                        #print('126', obstab['observation_value'][i] )
                    ztype[j] = 1
                elif zct == 2: 
                    zlist[j] = z_coordinate[i]
                    plist[j] = z_to_p_ifs(zlist[j])
                    z_coordinate[i] = plist[j]
                    ztype[j] = 1
                    #if obstab['observed_variable'][i] == 126:
                        #print('126', obstab['observation_value'][i] )
                else:
                    spurious_z += 1
                    plist[j] = z_coordinate[i]
                    zlist[j] = np.nan
                    ztype[j] = zct
                isave[j] = i
                j += 1
            else:
                z_coordinate[i] = z_coordinate[i-1]
            i += 1
        
        
        ix = np.where((obstab['z_coordinate_type'][k:i]!=1) & (obstab['z_coordinate_type'][k:i]==126))[0]
        if len(ix) > 0:
            
            obstab['observation_value'][k + ix] = np.nan
        
        if j == 0:
            k = i
            rts = obstab['date_time'][i]
            continue
        plj = plist[:j]
        pmax = np.max(plj)
        pmin = np.min(plj)
        
        
        for ip in range(ps.shape[0]):
            if ps[ip] > pmin and ps[ip] < pmax:
                if ps[ip] not in plj:
                    plist[j] = ps[ip]
                    j += 1
        
        idx = np.argsort(plist[:j])
        plist = plist[idx]
        isave = isave[idx]

        fobs = np.empty((isave.shape[0], allvar.shape[0]), dtype=np.float32) # +15 to accomodate for additional standard pressure levels
        fidx = np.empty((isave.shape[0], allvar.shape[0]), dtype=np.int64) # +15 to accomodate for additional standard pressure levels
        fobsvar = np.empty((isave.shape[0], allvar.shape[0]), dtype=np.int64) # +15 to accomodate for additional standard pressure levels
        fidx.fill(-1)
        fobsvar.fill(-1)
        ffg_depar = np.empty_like(fobs) # +15 to accomodate for additional standard pressure levels
        fan_depar = np.empty_like(fobs)
        fbc = np.empty_like(fobs)
        fbc_fg = np.empty_like(fobs)
        fcf = np.zeros((isave.shape[0], allvar.shape[0]), dtype=np.int32)
        fcm = np.empty_like(fcf)
        int_min = -2147483647 - 1
        fcf[:] = int_min
        
        for f in fobs, ffg_depar, fan_depar, fbc: 
            f.fill(np.nan)
        
        i = k
        j = -1
        porig = -1
        while i < oovs and obstab['date_time'][i] ==rts:
          
            if z_coordinate[i]!=porig:
                j += 1
                porig = z_coordinate[i]
                ii = np.searchsorted(plist, porig)
                if ii == plist.shape[0]:
                    ii -= 1
                    ##print('x', plist[-1], porig)
            vi = alli[obstab['observed_variable'][i]]
            fobs[ii, vi] = obstab['observation_value'][i]
            fobsvar[ii, vi] = obstab['observed_variable'][i]
            fidx[ii, vi] = oidx[i]
            ffg_depar[ii, vi] = loaded_feedback['fg_depar@body'][i]
            fan_depar[ii, vi] = loaded_feedback['an_depar@body'][i]
            fbc[ii, vi] = loaded_feedback['biascorr@body'][i]
            fbc_fg[ii, vi] = loaded_feedback['biascorr_fg@body'][i]
            i += 1
            
        foo = fobs.copy()
        
        #fl = [np.sum(~np.isnan(fobs[:, k])) for k in range(7, 11)]
        #if fl[2] != fl[3]:
            #print(fobs[:, 7:11])
            #x = 0
        ## convert ws,wd to u,v before vertical interpolation        
        mask =   np.isnan(fobs[:, 7]+fobs[:, 8])
        if np.any(mask):
            fobs[mask, 7], fobs[mask, 8], fobs[mask, 9], fobs[mask, 10] = touv(fobs[mask, 9] , fobs[mask, 10])
            fcm[mask] = 1
            fcf[mask] = 0

        for f in fobs, ffg_depar, fan_depar, fbc:
            lvint(plist, f[:isave.shape[0]], np.array((0, 1,2,3, 4, 5, 7, 8))) # interpolate vertically all variables except specific humidity, wind speed and wind direction
        
        ####if np.any(~np.isnan(fobs[:, 2])):
            
            ####plt.semilogy(fobs[:, 2], plist/100., 'o')
            ####plt.semilogy(foo[:, 2], plist/100., 'x')
            ####plt.ylim([1020., 500.])
            ####plt.show()
            #### convert to missing humidity variables
            
        tmask = ~np.isnan(fobs[:, 2])
        if np.sum(tmask) > 0:
            
            # dpd present, rh missing
            mask =   tmask & (~np.isnan(fobs[:, 3])) & (np.isnan(fobs[:, 5]))
            if np.any(mask):
                dp = fobs[mask, 2] - fobs[mask, 3]
                fobs[mask, 5] = Sonntag(dp) / Sonntag(fobs[mask, 2]) 
                fcm[mask, 5] = 2
                fcf[mask, 5] = 0
    
            # dpd present, dp missing
            mask =  tmask & ( ~np.isnan(fobs[:, 3])) & (np.isnan(fobs[:, 4]))
            if np.any(mask):
                fobs[mask, 4] = fobs[mask, 2] - fobs[mask, 3]
                fcm[mask, 4] = 2
                fcf[mask, 4] = 0

            # td present, rh missing
            mask =  tmask & ( ~np.isnan(fobs[:, 4])) & (np.isnan(fobs[:, 5]))
            if np.any(mask):
                fobs[mask, 5] = Sonntag(fobs[mask, 4]) / Sonntag(fobs[mask, 2]) 
                fcm[mask, 5] = 2
                fcf[mask, 5] = 0
    
            # sh present, rh missing
            mask = tmask & (  ~np.isnan(fobs[:, 6])) & (np.isnan(fobs[:, 5])) & (plist > 0.)
            if np.any(mask):
                vpdata = sh2vap(fobs[mask, 6], plist[mask])
                fobs[mask, 5] = vpdata / svp(fobs[mask, 2], p=plist[mask])
                fcm[mask, 5] = 4
                fcf[mask, 5] = 0
                
            #rh present, dp missing
            mask =  tmask & ( ~np.isnan(fobs[:, 5])) & (np.isnan(fobs[:, 4]))
            if np.any(mask):
                vpdata = fobs[mask, 5] * np.exp(liquid(fobs[mask, 2]))  #Sonntag(fobs[mask, 2])
                fobs[mask, 4] = dewpoint_Sonntag(vpdata)
                fcm[mask, 4] = 3
                fcf[mask, 4] = 0
    
            #rh present, td, dpd missing
            mask =  tmask & ( ~np.isnan(fobs[:, 5])) & (np.isnan(fobs[:, 3]))
            if np.any(mask):
                vpdata = fobs[mask, 5] * np.exp(liquid(fobs[mask, 2]))  #Sonntag(fobs[mask, 2])
                fobs[mask, 3] = fobs[mask, 2] - dewpoint_Sonntag(vpdata)
                fcm[mask, 3] = 3
                fcf[mask, 3] = 0
    
            # rh present, q missing
            mask =  tmask & ( ~np.isnan(fobs[:, 5])) & (np.isnan(fobs[:, 6])) & (plist > 0.)
            if np.any(mask):
                vpdata = fobs[mask, 5] * Sonntag(fobs[mask, 2])
                if np.any(vpdata>0.2*plist[mask]):
                    #print('vapor pressure too large')
                    fobs[mask, 6] = np.nan
                else:
                    fobs[mask, 6] = vap2sh(vpdata, plist[mask])
                if np.any(fobs[mask, 6]<0):
                    raise ValueError('negative specific humidity')
                fcm[mask, 6] = 3
                fcf[mask, 6] = 0
    
    
            
        # convert vertically interpolated u,v into ws, wd
        mask =   ~np.isnan(fobs[:, 7]+   fobs[:, 8]) & np.isnan(foo[:, 9]+foo[:, 10])
        if np.any(mask):
            fobs[mask, 10], fobs[mask, 9] = wswd(fobs[mask, 7] , fobs[mask, 8])
            fcm[mask, 9:] = 2
            fcf[mask, 9:] = 0
            
            
        #fl = [np.sum(~np.isnan(fobs[:, k])) for k in range(7, 11)]
        #if fl[2] != fl[3]:
            #print(fobs[:, 7:11])
            #x = 0
            
        ##idz = np.where(~np.isnan(foo))
        ##if len(idz) > 0:
            ##if np.sum(foo[idz]!=fobs[idz]) > 0:
                ##raise ValueError
    
            ### make sure dpd is consistent with temperature, dewpoint
            ##if np.any(fobs[:, 3] != fobs[:, 2] - fobs[:, 4]):
                
                ##fobs[:, 3] = fobs[:, 2] - fobs[:, 4]
        #gcount = 0
        mask = ~np.isnan(foo)
        if np.all(np.isnan(fobs)):
            for j in range(len(plist)):
                if np.any(mask[j, 7:]):
                    mask[j, 7:] = True
            
        for j in range(len(plist)):
            for jj in range(11):
                #if foo[j, jj] == foo[j, jj]:
                    #a_obstab['observation_value'][ak] = fobs[j, jj]
                    
                    #a_obstab['z_coordinate'][ak] = plist[j]
                    #a_obstab['z_coordinate_type'][ak] = 1
                    #a_obstab['observed_variable'][ak] = allvar[jj]        
                    #a_obstab['date_time'][ak] = rts
                    #a_obstab['conversion_flag'][ak] = fcf[j, jj]
                    #a_obstab['conversion_method'][ak] = fcm[j, jj]
                    #a_loaded_feedback['fg_depar@body'][ak] = ffg_depar[j, jj]
                    #a_loaded_feedback['an_depar@body'][ak] = fan_depar[j, jj]
                    #a_loaded_feedback['biascorr@body'][ak] = fbc[j, jj]
                    #a_loaded_feedback['biascorr_fg@body'][ak] = fbc_fg[j, jj]
                    #orig[ak] = fidx[j, jj]
                    #if a_obstab['observation_value'][ak] != obstab['observation_value'][orig[ak]]:
                        #if jj == 8 and ~np.isnan(foo[j, 7]) or jj == 7 and ~np.isnan(foo[j, 8]):                            
                            #raise ValueError
                    #ak += 1
            #for jj in range(11):
                #if foo[j, jj] != foo[j, jj] and fobs[j, jj] == fobs[j, jj]:
                if mask[j, jj] or (fobs[j, jj] == fobs[j, jj]):
                    a_obstab['observation_value'][ak] = fobs[j, jj]
                    
                    a_obstab['z_coordinate'][ak] = plist[j]
                    a_obstab['z_coordinate_type'][ak] = 1
                    a_obstab['observed_variable'][ak] = allvar[jj]        
                    a_obstab['date_time'][ak] = rts
                    a_obstab['conversion_flag'][ak] = fcf[j, jj]
                    a_obstab['conversion_method'][ak] = fcm[j, jj]
                    a_loaded_feedback['fg_depar@body'][ak] = ffg_depar[j, jj]
                    a_loaded_feedback['an_depar@body'][ak] = fan_depar[j, jj]
                    a_loaded_feedback['biascorr@body'][ak] = fbc[j, jj]
                    a_loaded_feedback['biascorr_fg@body'][ak] = fbc_fg[j, jj]
                    #orig[ak] = fidx[j, jj]
                    ak += 1
                    #gcount += 1
        #if gcount == 0:
            #print('record empty after deletions')
                

            
        k = i
        if i < oovs:
            
            rts=obstab['date_time'][i]

    #test code, do not delete
    #k = 0
    #for jj in range(a_obstab['observation_value'].shape[0]):#[:20]:
        
        #if orig[j] != -1:
            #if a_obstab['observation_value'][j] != obstab['observation_value'][orig[j]]:
                #print(j, a_obstab['date_time'][j], obstab['date_time'][orig[j]],
                    #a_obstab['z_coordinate'][j], obstab['z_coordinate'][orig[j]], 
                    #a_obstab['observation_value'][j], obstab['observation_value'][orig[j]])
            #else:
                #k += 1
    #if k != obstab['observation_value'].shape[0]:
        #print('too few values')

    
        
    return a_obstab[:ak], a_loaded_feedback[:ak], ak, orig[:ak] #addedvar[:l+1]



def offline_fb3(fpattern,fns, p,pn, fdict,ts, obstype, latorig,lonorig,z, refs,ans):
    #from scipy.interpolate import RectBivariateSpline
    tt=time.time()
    
    if ans[3] == ans[2]: # no data
        return
    lons = lonorig[ans[2]:ans[3]]
    lats = latorig[ans[2]:ans[3]]
    tss = ts[ans[2]:ans[3]]
    obss = obstype[ans[2]:ans[3]]
    zs = z[ans[2]:ans[3]] / 100.
    try:

        fn = fpattern.format(ans[0],ans[1])
        with h5py.File(fn,'r') as g:
            
            if not fdict:
                fdict['level'] = g['level'][:]
                fdict['pidx'] = np.searchsorted(fdict['level'],refs['level'])
                
                fdict['longitude']= g['longitude'][:]
                fdict['latitude']= g['latitude'][:]
                
    
                tunits=g['time'].attrs['units']
                try:
                    tunits=tunits.decode('latin1')
                except:
                    pass
                fdict['tunits']=tunits.split()
    
            fdict['attribs'] =  dict(g[p].attrs)
            try:
                
                del fdict['attribs']['DIMENSION_LIST']
            except:
                pass
            for k,v in g.attrs.items():
                fdict['attribs'][k]=v
            fdict['attribs']['source_filename']=fn
                
            pres = fdict['level']
            pidx = fdict['pidx']
            sh = g[p].shape


            lonmin = np.nanmin(lons)
            lonmax = np.nanmax(lons)
            latmin = np.nanmin(lats)
            latmax = np.nanmax(lats)
            dx=360/sh[1]
            dy=(fdict['latitude'][0] - fdict['latitude'][-1])/(sh[0]-sh[0] % 2)
                
            if lonmin < 0:
                lonmin += 360.
                lons += 360.
            if lonmax < 0:
                lonmax += 360.
            try:
                ixrefmin=int(np.floor(lonmin/dx))
                ixrefmax=int(np.floor(lonmax/dx))
                if latmin== -90.:
                    latmin = fdict['latitude'][-1]
                if latmax == 90.:
                    latmax = fdict['latitude'][0]
                iyrefmin= int(np.floor((fdict['latitude'][0]-latmax)/dy))
                iyrefmax = int(np.floor((fdict['latitude'][0]-latmin)/dy))
            except Exception as e:
                print('spurious lat, lon', e)
                return
                
            if ixrefmax >= ixrefmin:
                
                xsten = np.arange(ixrefmin-1, ixrefmax+3)#[ix - 1, ix, ix + 1, ix + 2])
            else:
                ixrefmax = int(np.floor(lonmax+360./dx))
                xsten = np.arange(ixrefmin-1, ixrefmax+3) % fdict['longitude'].shape[0]
            ysten = np.arange(iyrefmin-2, iyrefmax+3)#[iy - 1, iy, iy + 1, iy + 2])
            ysten[ysten>(sh[0]-1)] = sh[0] - 1
            ysten[ysten<0] = 0
            xsten[xsten>(sh[1]-1)] -= sh[1]
            xsten[xsten<0] += sh[1]

            if xsten[0] == xsten[-1] - xsten.shape[0] + 1:
                tera5=g[p][ysten[0]:ysten[-1] + 1,xsten[0]:xsten[-1] + 1,:,:]
                tlon = fdict['longitude'][:][xsten[0]:xsten[-1] + 1]
            else:
                tera5=g[p][ysten[0]:ysten[-1] + 1, :, :, :][:,xsten,:,:]
                tlon = fdict['longitude'][:][xsten]
            tlat = fdict['latitude'][ysten[0]:ysten[-1] + 1]
            
            try:
                mv=np.where(tera5==fdict['attribs']['missing_value'])
                tera5 = tera5 * np.float32(fdict['attribs']['scale_factor'])+np.float32(fdict['attribs']['add_offset'])
                tera5[mv]=np.nan
            except:
                pass

            if '20CRv3' in fpattern:
                for it in range(tera5.shape[2]):                       
                    for ip in range(len(pidx)):
                        if latmax > 90 - 2 * dy:
                            tera5[:,:,it,pidx[ip]]+=refs[p][ans[1] - 1, ysten[0]:ysten[-1] + 1,:, it % 8, ip][:, xsten]
                        elif latmin < -90 + 2 * dy:
                            tera5[:,:,it,pidx[ip]]+=refs[p][ans[1] - 1, ysten[0]:ysten[-1] + 1,:,it % 8, ip][:, xsten]
                        else:    
                            if xsten[0] == xsten[-1]- xsten.shape[0] + 1:
                                tera5[:,:,it,pidx[ip]]+=refs[p][ans[1] - 1, ysten[0]:ysten[-1] + 1,xsten[0]:xsten[-1] + 1,it % 8, ip]
                            else:
                                tera5[:,:,it,pidx[ip]]+=refs[p][ans[1] - 1, ysten[0]:ysten[-1] + 1,xsten,it % 8, ip].T
            #attribs['missing_value'] = np.nan
            #for k,v in g[p].attrs.items():
                #if k not in ['scale_factor','add_offset','DIMENSION_LIST','_FillValue','name','long_name','standard_name','units']:
                    #attribs[k]=v
                    #if k=='missing_value':
                        #attribs[k]=np.nan
            tunits = fdict['tunits']
            if tunits[0]=='hours':
                secs=np.int64(g['time'][:])*3600
            if tunits[-2]!='1900-01-01':
                if tunits[-2] =='1800-01-01':
                    offset=datetime.strptime(tunits[-2],'%Y-%m-%d')-datetime(1900,1,1) #+ timedelta(days=1)
                else:
                    offset=datetime.strptime(tunits[-2],'%Y-%m-%d')-datetime(1900,1,1)
                    
                secs=secs+int(offset.total_seconds())
                delta = secs[1] - secs[0]
                if datetime(1900, 1, 1) + timedelta(seconds=int(secs[0])) != datetime(ans[0], ans[1], 1):
                    secs = (datetime(ans[0], ans[1], 1) - datetime(1900, 1, 1) ).total_seconds() + np.arange(secs.shape[0]) *delta
                    print('modified timestamp of '+fn)
                

#        print(val[0, 0],'\n', ans[0, 0],'\n', time.time()-tt)
        
        an1 = [ans[0], ans[1] + 1]
        if an1[1] == 13:
            an1 = [ans[0] + 1, 1]
        try:
            
            with h5py.File(fpattern.format(an1[0],an1[1]),'r') as g:
                
                if xsten[0] == xsten[-1] - xsten.shape[0] + 1:
                    tera51=g[p][ysten[0]:ysten[-1] + 1,xsten[0]:xsten[-1] + 1,0:1,:]
                else:
                    tera51=g[p][ysten[0]:ysten[-1] + 1, :, :, :][:,xsten,0:1,:]
                
                try:
                    mv=np.where(tera51==fdict['attribs']['missing_value'])
                    tera51 =tera51 * np.float32(g[p].attrs['scale_factor'])+ np.float32(g[p].attrs['add_offset'])
                    tera51[mv]=np.nan
                except:
                    pass
    
                if '20CRv3' in fpattern:
                    for it in 0, :                       
                        for ip in range(len(pidx)):
                            if latmax > 90 - 2 * dy:
                                tera51[:,:,it,pidx[ip]]+=refs[p][an1[1] - 1, ysten[0]:ysten[-1] + 1,:, it % 8, ip][:, xsten]
                            elif latmin < -90 + 2 * dy:
                                tera51[:,:,it,pidx[ip]]+=refs[p][an1[1] - 1, ysten[0]:ysten[-1] + 1,:,it % 8, ip][:, xsten]
                            else:    
                                if xsten[0] == xsten[-1]- xsten.shape[0] + 1:
                                    tera51[:,:,it,pidx[ip]]+=refs[p][an1[1] - 1, ysten[0]:ysten[-1] + 1,xsten[0]:xsten[-1] + 1,it % 8, ip]
                                else:
                                    tera51[:,:,it,pidx[ip]]+=refs[p][an1[1] - 1, ysten[0]:ysten[-1] + 1,xsten,it % 8, ip].T
            tera5 = np.concatenate((tera5, tera51), axis=2)
            secs = np.concatenate((secs, secs[[-1]]+secs[-1]-secs[-2]))

        except FileNotFoundError as e:
            print(fpattern.format(an1[0],an1[1]), ' extrapolating last hours of month')
            tera5 = np.concatenate((tera5, tera5[:, :, -1:, :]), axis=2)
            secs = np.concatenate((secs, secs[[-1]]+secs[-1]-secs[-2]))
            pass

        #print(time.time() - tt)
        if False:
            #X, Y = np.meshgrid(g['longitude'], g['latitude'])
            val = np.empty(sh[2:])
            off = 2
            #print('lin', time.time() - tt)
            #tt = time.time()
            gt = g[p][iy +1 - off:iy + 1 + off, ix + 1- off:ix + 1 + off, :, :][::-1]
            glat = g['latitude'][iy + 1- off:iy + 1+ off][::-1].copy()
            glon = g['longitude'][ix + 1- off:ix + 1+ off].copy()

            for it in range(gt.shape[2]):

                for ip in range(gt.shape[3]):

                    interp_spline = RectBivariateSpline( glat,glon, gt[:, :, it, ip].reshape(gt.shape[:2]))
                    val[it, ip] = interp_spline(lat, lon)

            val =val*np.float32(g[p].attrs['scale_factor'])+np.float32(g[p].attrs['add_offset'])

            #print(time.time() - tt)
        
        reatab = np.zeros_like(lons)
        #params = {'t': 126}
        tpi = tera5[:, :, :, pidx]
        pri = pres[pidx]
        idx = triint(tpi, reatab, tlon, tlat, lons, lats, secs, tss,obss,pri,zs, pn, dy, fns)

        print(os.path.basename(fns), ans,p, reatab[0], time.time()-tt)
    except FileNotFoundError as e:
        print(fn, 'not available, continuing')
        return 
    return reatab[:len(idx)], ans[2] + idx, secs,pres,fdict['attribs']

def offline_fb4(fpattern,p,pn, fdict,ts, obstype, latorig,lonorig,z, refs,gactor, ans):
    #from scipy.interpolate import RectBivariateSpline
    tt=time.time()
    
    if ans[3] == ans[2]: # no data
        return
    lons = lonorig[ans[2]:ans[3]]
    lats = latorig[ans[2]:ans[3]]
    tss = ts[ans[2]:ans[3]]
    obss = obstype[ans[2]:ans[3]]
    zs = z[ans[2]:ans[3]] / 100.
    try:

        fn = fpattern.format(ans[0],ans[1])
        if type(gactor) == dict:
            gact = gactor
            #print(gact[os.path.basename(fn)].content.remote())
        else:
            refs = ray.get(refs)
            gact = ray.get(gactor)
            #print(gact[os.path.basename(fn)].content.remote())
        #print(ray.get(gact[os.path.basename(fn)].content.remote()))
        g = gact[os.path.basename(fn)]
        if True: #with h5py.File(fn,'r') as g:
            
            if not fdict:
                fdict['level'] = ray.get(g.full.remote('level')) #g['level'])
                fdict['pidx'] = np.searchsorted(fdict['level'],refs['level'])
                
                fdict['longitude']= ray.get(g.full.remote('longitude'))#g['longitude'][:]
                fdict['latitude']= ray.get(g.full.remote('latitude'))#g['latitude'][:]
                
    
                tunits=ray.get(g.vattr.remote('time', 'units')) #g['time'].attrs['units']
                try:
                    tunits=tunits.decode('latin1')
                except:
                    pass
                fdict['tunits']=tunits.split()
    
            fdict['attribs'] =  ray.get(g.vattrs.remote(p)) #dict(g[p].attrs)
            try:
                
                del fdict['attribs']['DIMENSION_LIST']
            except:
                pass
            #for k,v in g.attrs.items():
                #fdict['attribs'][k]=v
            fdict['attribs']['source_filename']=fn
                
            pres = fdict['level']
            pidx = fdict['pidx']
            sh = ray.get(g.shape.remote(p)) #g[p].shape


            lonmin = np.min(lons)
            lonmax = np.max(lons)
            latmin = np.min(lats)
            latmax = np.max(lats)
            dx=360/sh[1]
            dy=(fdict['latitude'][0] - fdict['latitude'][-1])/(sh[0]-sh[0] % 2)
                
            if lonmin < 0:
                lonmin += 360.
                lons += 360.
            if lonmax < 0:
                lonmax += 360.
            ixrefmin=int(np.floor(lonmin/dx))
            ixrefmax=int(np.floor(lonmax/dx))
            if latmin== -90.:
                latmin = fdict['latitude'][-1]
            if latmax == 90.:
                latmax = fdict['latitude'][0]
            iyrefmin= int(np.floor((fdict['latitude'][0]-latmax)/dy))
            iyrefmax = int(np.floor((fdict['latitude'][0]-latmin)/dy))
            if ixrefmax >= ixrefmin:
                
                xsten = np.arange(ixrefmin-1, ixrefmax+3)#[ix - 1, ix, ix + 1, ix + 2])
            else:
                ixrefmax = int(np.floor(lonmax+360./dx))
                xsten = np.arange(ixrefmin-1, ixrefmax+3) % fdict['longitude'].shape[0]
            ysten = np.arange(iyrefmin-2, iyrefmax+3)#[iy - 1, iy, iy + 1, iy + 2])
            ysten[ysten>(sh[0]-1)] = sh[0] - 1
            ysten[ysten<0] = 0
            xsten[xsten>(sh[1]-1)] -= sh[1]
            xsten[xsten<0] += sh[1]

            if xsten[0] == xsten[-1] - xsten.shape[0] + 1:
                tera5=ray.get(g.slice2d.remote(p, (slice(ysten[0], ysten[-1] + 1), slice(xsten[0], xsten[-1] + 1)))).copy() #g[p][ysten[0]:ysten[-1] + 1,xsten[0]:xsten[-1] + 1,:,:]
                #tera5=g[p][ysten[0]:ysten[-1] + 1,xsten[0]:xsten[-1] + 1,:,:]
                tlon = fdict['longitude'][:][xsten[0]:xsten[-1] + 1]
            else:
                tera5=ray.get(g.slice1d.remote(p, slice(ysten[0], ysten[-1] + 1)))[:,xsten,:,:] #g[p][ysten[0]:ysten[-1] + 1, :, :, :][:,xsten,:,:]
                tlon = fdict['longitude'][:][xsten]
            tlat = fdict['latitude'][ysten[0]:ysten[-1] + 1]
            
            try:
                mv=np.where(tera5==fdict['attribs']['missing_value'])
                tera5 = tera5 * np.float32(fdict['attribs']['scale_factor'])+np.float32(fdict['attribs']['add_offset'])
                tera5[mv]=np.nan
            except:
                pass

            if '20CRv3' in fpattern:
                for it in range(tera5.shape[2]):                       
                    for ip in range(len(pidx)):
                        if latmax > 90 - 2 * dy:
                            tera5[:,:,it,pidx[ip]]+=refs[p][ans[1] - 1, ysten[0]:ysten[-1] + 1,:, it % 8, ip][:, xsten]
                        elif latmin < -90 + 2 * dy:
                            tera5[:,:,it,pidx[ip]]+=refs[p][ans[1] - 1, ysten[0]:ysten[-1] + 1,:,it % 8, ip][:, xsten]
                        else:    
                            if xsten[0] == xsten[-1]- xsten.shape[0] + 1:
                                tera5[:,:,it,pidx[ip]]+=refs[p][ans[1] - 1, ysten[0]:ysten[-1] + 1,xsten[0]:xsten[-1] + 1,it % 8, ip]
                            else:
                                tera5[:,:,it,pidx[ip]]+=refs[p][ans[1] - 1, ysten[0]:ysten[-1] + 1,xsten,it % 8, ip].T
            #attribs['missing_value'] = np.nan
            #for k,v in g[p].attrs.items():
                #if k not in ['scale_factor','add_offset','DIMENSION_LIST','_FillValue','name','long_name','standard_name','units']:
                    #attribs[k]=v
                    #if k=='missing_value':
                        #attribs[k]=np.nan
            tunits = fdict['tunits']
            if tunits[0]=='hours':
                secs=np.int64(ray.get(g.full.remote('time')))*3600
                #secs=np.int64(g['time'][:])*3600
            if tunits[-2]!='1900-01-01':
                if tunits[-2] =='1800-01-01':
                    offset=datetime.strptime(tunits[-2],'%Y-%m-%d')-datetime(1900,1,1) #+ timedelta(days=1)
                else:
                    offset=datetime.strptime(tunits[-2],'%Y-%m-%d')-datetime(1900,1,1)
                    
                secs=secs+int(offset.total_seconds())
                delta = secs[1] - secs[0]
                if datetime(1900, 1, 1) + timedelta(seconds=int(secs[0])) != datetime(ans[0], ans[1], 1):
                    secs = (datetime(ans[0], ans[1], 1) - datetime(1900, 1, 1) ).total_seconds() + np.arange(secs.shape[0]) *delta
                    print('modified timestamp of '+fn)
                

#        print(val[0, 0],'\n', ans[0, 0],'\n', time.time()-tt)
        
        an1 = [ans[0], ans[1] + 1]
        if an1[1] == 13:
            an1 = [ans[0] + 1, 1]
        try:
            
            with h5py.File(fpattern.format(an1[0],an1[1]),'r') as g:
                
                if xsten[0] == xsten[-1] - xsten.shape[0] + 1:
                    tera51=g[p][ysten[0]:ysten[-1] + 1,xsten[0]:xsten[-1] + 1,0:1,:]
                else:
                    tera51=g[p][ysten[0]:ysten[-1] + 1, :, :, :][:,xsten,0:1,:]
                
                try:
                    mv=np.where(tera51==fdict['attribs']['missing_value'])
                    tera51 =tera51 * np.float32(g[p].attrs['scale_factor'])+ np.float32(g[p].attrs['add_offset'])
                    tera51[mv]=np.nan
                except:
                    pass
    
                if '20CRv3' in fpattern:
                    for it in 0, :                       
                        for ip in range(len(pidx)):
                            if latmax > 90 - 2 * dy:
                                tera51[:,:,it,pidx[ip]]+=refs[p][an1[1] - 1, ysten[0]:ysten[-1] + 1,:, it % 8, ip][:, xsten]
                            elif latmin < -90 + 2 * dy:
                                tera51[:,:,it,pidx[ip]]+=refs[p][an1[1] - 1, ysten[0]:ysten[-1] + 1,:,it % 8, ip][:, xsten]
                            else:    
                                if xsten[0] == xsten[-1]- xsten.shape[0] + 1:
                                    tera51[:,:,it,pidx[ip]]+=refs[p][an1[1] - 1, ysten[0]:ysten[-1] + 1,xsten[0]:xsten[-1] + 1,it % 8, ip]
                                else:
                                    tera51[:,:,it,pidx[ip]]+=refs[p][an1[1] - 1, ysten[0]:ysten[-1] + 1,xsten,it % 8, ip].T
            tera5 = np.concatenate((tera5, tera51), axis=2)
            secs = np.concatenate((secs, secs[[-1]]+secs[-1]-secs[-2]))

        except FileNotFoundError as e:
            print(fpattern.format(an1[0],an1[1]), ' extrapolating last hours of month')
            tera5 = np.concatenate((tera5, tera5[:, :, -1:, :]), axis=2)
            secs = np.concatenate((secs, secs[[-1]]+secs[-1]-secs[-2]))
            pass

        #print(time.time() - tt)
        if False:
            #X, Y = np.meshgrid(g['longitude'], g['latitude'])
            val = np.empty(sh[2:])
            off = 2
            #print('lin', time.time() - tt)
            #tt = time.time()
            gt = g[p][iy +1 - off:iy + 1 + off, ix + 1- off:ix + 1 + off, :, :][::-1]
            glat = g['latitude'][iy + 1- off:iy + 1+ off][::-1].copy()
            glon = g['longitude'][ix + 1- off:ix + 1+ off].copy()

            for it in range(gt.shape[2]):

                for ip in range(gt.shape[3]):

                    interp_spline = RectBivariateSpline( glat,glon, gt[:, :, it, ip].reshape(gt.shape[:2]))
                    val[it, ip] = interp_spline(lat, lon)

            val =val*np.float32(g[p].attrs['scale_factor'])+np.float32(g[p].attrs['add_offset'])

            #print(time.time() - tt)
        
        reatab = np.zeros_like(lons)
        #params = {'t': 126}
        idx = triint(tera5[:, :, :, pidx], reatab, tlon, tlat, lons, lats, secs, tss,obss,pres[pidx],zs, pn, dy)

        print(ans,p, reatab[0], time.time()-tt)
    except FileNotFoundError as e:
        print(fn, 'not available, continuing')
        return 
    return reatab[:len(idx)], ans[2] + idx, secs,pres,fdict['attribs']


def load_20CRoffset(readict):
       
    from scipy.interpolate import RectBivariateSpline
    start = 1940
    try:
        with open(os.path.expandvars(wpath+'/rea/refs{}x.pkl'.format(start)),'rb') as f:
            refs=pickle.load(f)
    except:

        refs={}
        raw_20CR={}
        tt=time.time()
        l=0
        for iy in range(1940,1942):
            l+=1
            for im in range(1,13):
                k='era5'
                v=readict['era5']
                k20='20CRv3'
                v20=readict['20CRv3']
                for k,v in readict.items():
                    if k=='20CRv3':
                        continue
                    for ftype in v['ftype']:  
                        if 'fc' in ftype:
                            dtype='fc'
                            continue
                        else:
                            dtype='an'
                        if dtype not in readict[k].keys():
                            readict[k][dtype]={}

                        for p,pn in v['param'].items():
                            if k!='JRA55':

                                fpattern=v['path']+v['prefix']+ftype[:-1]+v['glue']+'{}{:0>2}'+v['glue']+pn+v['suffix']+'.nc'
                                if p in ['q', 'z']:
                                    continue
                            else:
                                fpattern=v['path']+v['prefix']+ftype+v['glue']+pn+'.reg_tl319.{}{:0>2}'+v['glue']+v['suffix']+'.nc'

                            try:
                                with h5py.File(fpattern.format(iy,im),'r') as f:
                                    print(f.keys())
                                    fpattern20=v20['path']+v20['prefix']+v20['glue']+'{}{:0>2}'+v20['glue']+v20['param'][p]+v20['suffix']+'.nc'
                                    with h5py.File(fpattern20.format(iy,im),'r') as g:
                                        print(g.keys())
                                        tfac = f['time'].shape[0] // g['time'].shape[0]
                                        tt = time.time()
                                        ifield=f[p][::tfac, :, :, :]
                                        print(time.time() - tt)
                                        ifields = np.zeros((ifield.shape[2], ifield.shape[3], 8, ifield.shape[1]), dtype=np.float32)
                                        for i in range(ifields.shape[2]):
                                            for ip in range(ifields.shape[3]):
                                                ifields[:, :, i, ip] = np.mean(ifield[i::8, ip, :, :], axis=0)*np.float32(f[p].attrs['scale_factor'])+np.float32(f[p].attrs['add_offset'])
                                        print(time.time() - tt)
                                        ilon=f['longitude'][:]
                                        ilat=f['latitude'][:]
                                        iplev=f['level'][:]

                                        oplev=g['level'][:]
                                        pindex=np.searchsorted(oplev,iplev)
                                        oshape=g[p].shape
                                        if p not in refs.keys():
                                            refs[p]=np.zeros((12, oshape[0],oshape[1],8, ifields.shape[3]),dtype=ifields.dtype)
                                            raw_20CR[p]=np.zeros(refs[p].shape, dtype=ifields.dtype)
                                            refs['level']=iplev[:]
                                        print(time.time() - tt)
                                        gmem = g[p][:]
                                        print(time.time() - tt)
                                        for i in range(ifields.shape[2]):
                                            for ip in range(ifields.shape[3]):

                                                raw_20CR[p][im - 1, :, :, i, ip]+=np.mean(gmem[:,:,i::8, pindex[ip]], axis=2)

                                        if 'olon' not in locals():
                                            glon = g['longitude'][:]
                                            glat = g['latitude'][:]
                                            olon,olat=np.meshgrid(g['longitude'][:],g['latitude'][:])

                                        print(time.time() - tt)


                                ttt = time.time()
                                for it in range(ifields.shape[2]):
                                    for ip in range(ifields.shape[3]):
                                        interp_spline = RectBivariateSpline( -ilat, ilon, ifields[:, :, it, ip].reshape(ifields.shape[:2]))
                                        refs[p][im - 1, :,:,it, ip] += interp_spline(-glat, glon)
                                    print(iy,im,it, p,np.std(refs[p][im - 1, :, :, it, :]),np.std(ifields[:, :, it, :]),time.time()-tt)
                                print('x', time.time() - ttt)

                            except MemoryError as e:
                                print(e,'could not read')

        for r in ['t','u','v']:
            refs[r]/=l
            raw_20CR[r]/=l
            refs[r]-=raw_20CR[r]

        try:
            os.mkdir(wpath+'/rea')
        except:
            pass
        with open(os.path.expandvars(wpath+'/rea/refs{}x.pkl'.format(start)),'wb') as f:
            pickle.dump(refs,f)

    return refs    

ray_load_20CRoffset = ray.remote(load_20CRoffset)

def readictplot(readict, tasks, plevs, figprefix, marker='*'):
    
        
    for plev in plevs:
        for sfunc in tasks:
    
            for p in readict['obstypes'].keys():
            

            #plev = plevs[0]

                idx=np.where(np.logical_and(readict['obs']['obstype']==readict['obstypes'][p],readict['obs']['z_coordinate']==plev))[0]
                if p == 'u':
                    idxv = np.where(np.logical_and(readict['obs']['obstype']==readict['obstypes']['v'],readict['obs']['z_coordinate']==plev))[0]
                if len(idx) < 50:
                    continue
    
                for k,v in readict.items():
                    if 'obs' in k:
                        continue
                    if k =='era5fb':
                        i = 0
                    
                    if 'ftype' in v.keys():                    
                        iterable = v['ftype']
                    else:
                        iterable = v
                        
                    for ftype in iterable:  
                        if 'fc' in ftype:
                            dtype='fc'
                        else:
                            dtype='an'
    
                        if dtype in readict[k].keys():
    
                            rms=[]
                            rmsv = []
                            years=[]
                            ref = datetime(1900, 1, 1)
                            tsy = ref.year + np.floor(readict['obs']['date_time'][idx]/365.25/86400)
                            
                            q = np.nanquantile(readict['obs']['obs'][idx]-readict[k][dtype]['refvalues'][idx], (0.01, 0.99))
                            print(p, k, q)
                            ystart = (ref + timedelta(seconds=int(readict['obs']['date_time'][idx[0]]))).year
                            ystop = (ref + timedelta(seconds=int(readict['obs']['date_time'][idx[-1]]))).year + 1
                            smarker = marker
                            if ystop - ystart == 1:
                                marker = '*'
                            
                            for iy in range(ystart, ystop):
                                #print(iy)
    
                                idy=np.where(tsy==iy)[0]
                                if p != 'u':
                                    
                                    rv = np.clip(readict['obs']['obs'][idx[idy]] - readict[k][dtype]['refvalues'][idx[idy]], q[0], q[1])
                                else:
                                    ou = readict['obs']['obs'][idx[idy]]
                                    fu = readict[k][dtype]['refvalues'][idx[idy]]
                                    fv = readict[k][dtype]['refvalues'][idxv[idy]]
                                    ov = readict['obs']['obs'][idxv[idy]]
                                    ff=np.sqrt(ou**2+ov**2)
                                    o2=270-np.arctan2(ov,ou)*180/np.pi
                                    fd=270-np.arctan2(fv,fu)*180/np.pi
                                    v=o2-fd
                                    idxx = np.where(~np.isnan(v))[0]
                                    if len(idxx) > 0:
                                        idyy = np.where(v[idxx]>180)[0]
                                        if len(idyy) > 0:
                                            v[idxx[idyy]] -= 360
                                        idyy = np.where(v[idxx]<-180)[0]
                                        if len(idyy) > 0:
                                            v[idxx[idyy]] += 360
                                        idyy = np.where(np.abs(v[idxx]) > 90)[0]
                                        if len(idyy) > 0:
                                            v[idxx[idyy]] =np.nan # larger deviations than 90 degrees are highly unlikely to come from wrong north alignment, are therefore discarded
                                    
                                    idxx = np.where(~np.isnan(ff))[0]
                                    if len(idxx) > 0:
                                        idyy = np.where(ff[idxx]<2.0)[0]
                                        if len(idyy) > 0:        
                                            v[idxx[idyy]]=np.nan # wind directions for small wind speeds too uncertain
                                            fd[idxx[idyy]]=np.nan # wind directions for small wind speeds too uncertain
                                            o2[idxx[idyy]]=np.nan # wind directions for small wind speeds too uncertain
                                    rv = v
                                if len(idy)>50 and np.sum(~np.isnan(rv)) > 50:
                                    if sfunc == 'rms':
    
                                        rms.append(np.sqrt(np.nanmean((rv)**2)))
                                    elif sfunc == 'std':
                                        rms.append(np.nanstd(rv))
                                    else:
                                        rms.append(np.nanmean(rv))
    
                                else:
                                    rms.append(np.nan)
                                years.append(iy)
    
                            plt.plot(np.array(years),np.array(rms),'-'+marker, 
                                     label='obs -'+k+'_'+dtype+', '+sfunc+'= {:5.3f}'.format(np.sqrt(np.nanmean(np.array(rms)**2))))
    
                plt.title('Monthly '+sfunc+' reanalysis '+p+' departures, '+str(int(plev/100)) + ' hPa, '+figprefix.split('/')[-1].split('_')[0])
                plt.ylabel(sfunc+' ['+readict['obsunits'][p]+']')
                plt.legend()
                plt.tight_layout()
                #plt.xlim(1935, 1965)
                plt.savefig(figprefix+'_'+p+'_'+str(int(plev/100))+'_'+sfunc+'stats.png')
                plt.close()

def retrieve_anfg(f, out, out_fb, readict, ts, tslice, refs, gactor, out_name,path_to_gridded):
    #from scipy.interpolate import RectBivariateSpline


    #era5=refdiff('/raid60/scratch/leo/scratch/ERA5/gridded/','era5fc.{0}{1:02}.'+par['gid'],1979,1983,tidx,fieldsperday=2,
                    #fcstep=12,mrange=list(range(1,13)),tgrid=None,tempname=tempname)

    #jra55_58_c[im,ipar,ip,:,:]=Interp2d(jra55_58_m[im,ipar,ip,:,:], jlon, jlat, eelon, eelat, order=1)

    #readict={'era5':{'ftype':('t','fct'),'param':{'t':'130','u':'131','v':'132','q':'133','z':'129'},
                     #'path':os.path.expandvars('$RSCRATCH/era5/gridded/'),'prefix':'era5','suffix':'','glue':'.'},
             ##'CERA20C':{'ftype':('t',),'param':{'t':'130','u':'131','v':'132','q':'133','z':'129'},
               ##'path':os.path.expandvars('$RSCRATCH/CERA20C/'),'prefix':'CERA20C','suffix':'','glue':'.'},
                          ##'JRA55':{'ftype':('fcst_mdl',),'param':{'t':'011_tmp','u':'033_ugrd','v':'034_vgrd','q':'051_spfh'},
               ##'path':os.path.expandvars('$RSCRATCH/JRA55/split/'),'prefix':'test.','suffix':'grb','glue':'.'},
                        #'20CRv3':{'ftype':('',),'param':{'t':'TMP','u':'UGRD','v':'VGRD'},#,'q':'SPFH'},#,'z':'HGT'},
                         #'path':os.path.expandvars('$RSCRATCH/20CRv3/'),'prefix':'anl_meant','suffix':'_pres','glue':'_'},
               #}

    tt=time.time()
    cyear = int(out_name.split('/')[-2])
    outdir = os.path.dirname(out_name) + '/'
    fn = f.filename

    if True:
        try:

            
            otsu, ori=np.unique(out['date_time'], return_index=True)
            tsu, ri=np.unique(ts, return_index=True)
            
            ilat=f['observations_table']['latitude'][tslice]
            ilon=f['observations_table']['longitude'][tslice]
            ilon[ilon>360.] -= 360.
            ilon[ilon<0] += 360.

            lat = np.empty_like(out['observation_value'])
            lon = np.empty_like(lat)
            for l in range(len(otsu)-1):
                lat[ori[l]:ori[l+1]] = ilat[ri[l]]
                lon[ori[l]:ori[l+1]] = ilon[ri[l]]
            lat[ori[-1]:] = ilat[-1]
            lon[ori[-1]:] = ilon[-1]
                
                
            obstype=out['observed_variable']
            obs=out['observation_value']
            z=out['z_coordinate']
            latdis = lat + (100000. - z) * 3. / 100000.
            londis = lon + (100000. - z) * 3. / 100000.
            zt=out['z_coordinate_type']
            ofb=True
        except Exception as e:
            print(fn, cyear, e)
            return None
        try:

            bc = out_fb['biascorr@body']
            bc[np.isnan(bc)] = 0.
            o_minus_bg=out_fb['fg_depar@body'] + bc
            o_minus_an=out_fb['an_depar@body'] + bc
        except:
            ofb=False

    ##fix geopotential
    #tsu, ri=np.unique(ts, return_index=True)
    #for i in range(len(ri) - 1):
        #o = obstype[ri[i]:ri[i+1]]
        #idx = np.where(np.logical_and(o==117, zt[ri[i]:ri[i+1]]==1))[0]
        #if len(idx) > 0:
            #m = np.argmax(obs[ri[i]+idx])
            ##print(287 * 25 * np.log10(100000./z[ri[i]+idx[m]]) / obs[ri[i]+idx[m]])
            #if 287 * 25 * np.log10(100000./z[ri[i]+idx[m]]) / obs[ri[i]+idx[m]] > 0.1:
                #obs[ri[i]+idx] *=  9.80665
                #if np.any(obs[ri[i]+idx]>500000):
                    #print(i, len(ri), obs[ri[i]+idx[:]])
            ##print(i, len(ri))
                    #print('spurious')
        ##if i % 1000 == 0:
            ##print(i, len(ri))



    ref=datetime(1900,1,1)
    yms=[]
    oldyear=0
    oldmonth=0
    l = 0
    x = [ref+timedelta(seconds=int(k)) for k in otsu]
    while l < otsu.shape[0]:
        if (x[l].year!=oldyear or x[l].month!=oldmonth): #and x[l].year >1977 and x[l].year < 1983:
            m = 0
            
            oldmonth=x[l].month
            oldyear=x[l].year
            while x[l+m].year ==oldyear and x[l+m].month == oldmonth and l + m < otsu.shape[0] :
                m += 1
                if l + m ==len(x):
                    break
                
            if l +m == otsu.shape[0]:                
                yms.append((oldyear, oldmonth, ori[l], out['date_time'].shape[0]))
            else:
                yms.append((oldyear, oldmonth, ori[l], ori[l+m]))
            l += m - 1
        l += 1

    #try:
        #os.remove(out_name)
    #except:
        #pass
        
    obstypes={'t':ipar[85],'u':ipar[104],'v':ipar[105],'q':ipar[39],'z':ipar[117]}    
    obsunits={'t':'K','u':'m/s','v':'m/s','q':'g/kg','z':'m^2/s^2'}
    for k,v in readict.items():
        #if k!='20CRv3':
            #continue
        for ftype in v['ftype']:  
            if 'fc' in ftype:
                dtype='fc'
            else:
                dtype='an'
            if dtype not in readict[k].keys():
                readict[k][dtype]={}
                readict[k][dtype]['refvalues']=np.empty(obs.shape,dtype=np.float32)
                readict[k][dtype]['refvalues'].fill(np.nan)

            for p,pn in v['param'].items():
                #if p == 't':
                    #continue
                if k!='JRA55':

                    fpattern=v['path']+v['prefix']+ftype+v['glue']+'{}{:0>2}'+v['glue']+pn+v['suffix']+'.nc'
                else:
                    fpattern=v['path']+v['prefix']+ftype+v['glue']+pn+'.reg_tl319.{}{:0>2}'+v['glue']+v['suffix']+'.nc'
                #found=len(glob.glob(fpattern.format(*yms[0])))
                #print(fpattern.format(*yms[0]),found)
                #if found:
                fdict = {}
                if gactor is None:
                    func=partial(offline_fb3,fpattern,fn, p,obstypes[p], fdict,out['date_time'],obstype, lat,lon, z, refs)
                else:
                    func=partial(offline_fb4,fpattern,p,obstypes[p], fdict,out['date_time'],obstype, lat,lon, z, refs, gactor)
#                func=partial(offline_fb3,fpattern,p,fdict, lat+latdis,lon+londis,ts,refs)
                #tups=list(P.map(func,yms[:]))
                tups=list(map(func,yms[:])) # no multiprocessing because function is already mapped
                ntups=[]
                for t in tups:
                    if t is not None:
                        ntups.append(t)
                if len(ntups)>0:

                    press=ntups[0][3]
                    if np.max(press)>1000:
                        press/=100.

                    readict[k][dtype][p]={}
                    readict[k][dtype][p]['time']=np.concatenate([ntups[i][2] for i in range(len(ntups))])
                    #readict[k][dtype][p]['values']=np.concatenate([ntups[i][0] for i in range(len(ntups))])
                    
                    for i in range(len(ntups)):
                        
                        readict[k][dtype]['refvalues'][ntups[i][1]] = ntups[i][0]
                        
                    readict[k][dtype][p]['attribs']=ntups[0][4]

            #tr[1]=117  # should change
            #tr[2]=85
            #tr[3]=104
            #tr[4]=105
            #tr[7]=39 #spec hum
            #tr[29]=38 #relative hum
            #tr[59]=36 # dew point
            #tr[111]=106 #dd
            #tr[112]=107  #ff

            if dtype in readict[k].keys():  
                for p in readict[k][dtype].keys():
                    if p =='refvalues':
                        continue
                    #try:

                        #ttt = time.time()
                        #interp(readict[k][dtype]['refvalues'],obstype,z,ts,press,obstypes[p],
                               #readict[k][dtype][p]['time'],readict[k][dtype][p]['values'])
                        #print(k,dtype,p,time.time()-tt)
                    #except MemoryError as e:
                        #print(e)
                    ts = out['date_time']
                    plev = 85000
                    idx=np.where(np.logical_and(obstype==obstypes[p],z==plev))
                    if len(idx[0]) > 0:

                        try:

                            plt.subplot(2,1,1)
                            plt.plot(1900+ts[idx]/365.25/86400,obs[idx], label='obs')
                            rv = readict[k][dtype]['refvalues'][idx]
                            plt.plot(1900+ts[idx]/365.25/86400,rv, label=k)
                            tmin = int(1900 + np.floor(ts[idx[0][0]]/365.25/86400))
                            tmax = int(1900 + np.floor(ts[idx[0][-1]]/365.25/86400) + 1)
                            plt.xlim(1900+ts[idx[0][0]]/365.25/86400, 1900+ts[idx[0][-1]]/365.25/86400)
                            ax = plt.gca()
                            ax.set_xticks(np.arange(tmin, tmax+1))
                            plt.title(p+', '+fn.split('/')[-1].split('_')[0]+' '+str(plev//100) +' hPa')
                            plt.legend()
                            #plt.xlim(1935, 1965)
                            plt.subplot(2,1,2)
                            plt.plot(1900+ts[idx]/365.25/86400,obs[idx]-rv,
                                     label='offline mean:{:5.3f} rms= {:5.3f}'.format(np.nanmean(obs[idx]-rv), np.sqrt(np.nanmean((obs[idx]-rv)**2))))
                            plt.plot(1900+ts[idx]/365.25/86400,o_minus_an[idx]+rv-rv,
                                     label='online mean:{:5.3f} rms:{:5.3f}'.format(np.nanmean(o_minus_an[idx]+rv-rv), np.sqrt(np.nanmean((o_minus_an[idx]+rv-rv)**2))))
                            plt.title(p+', obs -'+k )
                            yscale = 10. ** np.floor(np.log10(np.max(np.abs(obs[idx]-rv)))) + 1
                            plt.ylim(-yscale, yscale)
                            ax = plt.gca()
                            ax.set_xticks(np.arange(tmin, tmax+1))
                            plt.legend()
                            plt.tight_layout()
                            fnp=fn.split('/')[-1].split('CEUAS_merged_v1.nc')[0]
                            plt.savefig(outdir +fnp+p+'_'+k+dtype+'.png')
                            plt.close()

                        except Exception as e:

                            print('plotting reference', e)


    readict['tslice'] = tslice
    print(os.path.basename(fn),cyear, 'feedback calculated', time.time() - tt)
    if not ofb:
        return readict

    readict['era5fb']={'ftype':['an','fc']}                  
    readict['era5fb']['an']={'refvalues':obs-o_minus_an }                  
    readict['era5fb']['fc']={'refvalues':obs-o_minus_bg }
    readict['obs'] = {'z_coordinate': z,'obstype': obstype,'date_time': out['date_time'], 'obs': obs,}
    readict['obstypes'] = obstypes
    readict['obsunits'] = obsunits

    try:
        os.mkdir(outdir)
    except:
        pass
    
    #readictplot(readict, ('mean', 'std','rms'), (30000, ), outdir+os.path.basename(fn).split('_')[0])

    #P.close()
    #P.join()
    #del P

    return readict

def dim_attach(g, k):
    for v in g[k].keys(): #var_selection:
        l=0            
        try:
            fvv=g[k][v]
            if 'string' not in v and v!='index':                    
                g[k][v].dims[l].attach_scale(g[k]['index'])
                #print(v,fvv.ndim,type(fvv[0]))
                if fvv.ndim==2 or type(fvv[0]) in [str,bytes,np.bytes_]:
                    slen=fvv.shape[1] #sdict[v]
                    #slen=10
                    g[k][v].dims[1].attach_scale(g[k]['string{}'.format(slen)])
        except MemoryError as e:
            print(g.filename.split('/')[-1],k, e)
            pass

def read_station_height(fn, jj, inventory_list):
    heights = np.zeros(1, dtype=np.float32)
    wigos = fn.split('/')[-1].split('_')[0]
    pd_list = []
    for fns in inventory_list:
        try:
            
            pd_list.append(pd.read_csv(fns, delimiter='\t'))
            print(fns.split('/')[-1], len(pd_list[-1]))
            matches = pd_list[-1].index[pd_list[-1]['primary_id']==wigos].tolist()
            #matches = pd_list[-1].index[pd_list[-1]['StationId']==wigos].tolist()
            for l in matches:
                #heights[:] = pd_list[-1]['Hha'][l]
                #print(heights[0])
                print(pd_list[-1]['start_date'][l], pd_list[-1]['elevation'][l])
                heights[:] = pd_list[-1]['elevation'][l]
        except Exception as e:
            print(fns, e)
    return heights
    #fpattern=path_to_gridded+'/era5fct.{}{:0>2}.130.nc'
    #func=partial(offline_fb,fpattern,lat,lon)
    #tfgs=list(map(func,yms))

def convert_missing(refs, gactor, wpath,cyear, fn):

    tt=time.time()
    print(fn.split('/')[-1], cyear,end=' ' )
    wpathy = wpath+'/' + str(cyear) + '/'
    targetfile = wpathy + fn.split('/')[-1]
    logfile = wpathy+'/log/'+fn.split('/')[-1]+".txt"
    tfile = wpath + '2022/0-20000-0-96207_CEUAS_merged_v1.nc'
    if os.path.isfile('x'+logfile):
        print('already processed')
        wtime = os.path.getmtime(logfile) 
        rtime = os.path.getmtime(tfile) #fn
        if rtime>wtime:
            print(logfile, 'older than some input, processing')
        else:
            return targetfile
    else:
        print('processing...')
    sys.path.insert(0,os.getcwd()+'/../resort/rasotools-master/')
    
    readict={'era5':{'ftype':('t','fct'),'param':{'t':'130','u':'131','v':'132','q':'133','z':'129'},
                     'path':os.path.expandvars('$RSCRATCH/era5/gridded/'),'prefix':'era5','suffix':'','glue':'.'},
             #'CERA20C':{'ftype':('t',),'param':{'t':'130','u':'131','v':'132','q':'133','z':'129'},
               #'path':os.path.expandvars('$RSCRATCH/CERA20C/'),'prefix':'CERA20C','suffix':'','glue':'.'},
                          #'JRA55':{'ftype':('fcst_mdl',),'param':{'t':'011_tmp','u':'033_ugrd','v':'034_vgrd','q':'051_spfh'},
               #'path':os.path.expandvars('$RSCRATCH/JRA55/split/'),'prefix':'test.','suffix':'grb','glue':'.'},
                        '20CRv3':{'ftype':('',),'param':{'t':'TMP','u':'UGRD','v':'VGRD'},#,'q':'SPFH'},#,'z':'HGT'},
                         'path':os.path.expandvars('$RSCRATCH/20CRv3/'),'prefix':'anl_meant','suffix':'_pres','glue':'_'},
               }
    
    #inventory_list = glob.glob('../merge/example_notebooks/new_stat_conf/*.dat')
    #inventory_list = glob.glob('../../meta/inventory_comparison_2/data/tables/vola_legacy_report.txt')

    nanlist = [float('nan'), np.nan, 0, -2147483648]

### code for checking if z_coordinate_type=2 exists in particular year
    if False:
        with h5py.File(fn,'r') as f:
            try:
                if f['recordtimestamp'].shape[0] == 1:                
                    ts=f['recordtimestamp'][:]
                else:
                    ts=f['recordtimestamp'][[0, -1]]
                lon = f['observations_table']['longitude'][:]
                #if np.max(lon) < -2. or np.min(lon) > 2.:
                    #return
                
                if True:                               
                #if f['observations_table']['date_time'].shape[0] == 1:
                    
                    #ts=f['observations_table']['date_time'][:]
                #else:
                    #ts=f['observations_table']['date_time'][[0, -1]]
                    ref = datetime(1900, 1, 1)
                    tstart = int((datetime(cyear, 1, 1) - ref).total_seconds())
                    tstop = int((datetime(cyear+1, 1, 1) - ref-timedelta(seconds=1)).total_seconds())
                    if tstop < ts[0] or tstart > ts[-1]:
        #                print(fn, cyear, 'year missing in obs records')
                        return
        
                    ts=f['recordtimestamp'][:]
                    ri = f['recordindex'][:]
                    
                    sta, sto = np.searchsorted(ts, (tstart, tstop))
    
                    if sto == ri.shape[0]:
                        tslice =slice(ri[sta], f['observations_table']['date_time'].shape[0])
                    else:
                        tslice =slice(ri[sta], ri[sto])
    
                    #tsliceold = slice(*np.searchsorted(f['observations_table']['date_time'][:], (tstart, tstop)))
                    #print (tslice, tsliceold)
                    #assert tslice == tsliceold
                    #tslice =slice(*ri[np.searchsorted(ts, (tstart, tstop))])
                    
                    #ts=f['observations_table']['date_time'][:]
                    #tslice =slice(*np.searchsorted(ts, (tstart, tstop)))
        
                    if tslice.stop ==tslice.start:
                        return
                    
                    n = np.sum(f['observations_table']['z_coordinate_type'][tslice]==0)
                    if n > 0:
                        print(fn, n, 'with z as coordinate')
                    else:
                        return
            except Exception as e:
                print(fn, e)
                with open(logfile.split('.nc')[0]+'.err', 'w') as fe:
                    fe.write(str(fn) +'\n'+ str(e) +'\n')
                    
                return

### code for checking if z_coordinate_type=2 exists in particular year
    rscratch='/mnt/users/scratch/leo/scratch/'
    try:
        os.mkdir(wpathy)
        time.sleep(0.1)
    except:
        pass

    
    inventory_list = glob.glob('../../meta/inventory_comparison_2/code/station_configuration/CUON_station_configuration_extended.csv')
    station_elevation = read_station_height(fn,100, inventory_list)
    print('read inventory', time.time()-tt)

     # wpath+fn.split('/')[-1]
    if os.path.isfile(targetfile):
        try:
            os.remove(targetfile)
        except:
            print('file could not be removed - overwriting will lead to errors')

    #print('writing reanalysis reference series')
    #for k in readict.keys():
        #if k in ('tslice', ):
            #continue
        #for dtype in 'an', 'fc':
            
            #try:
                #os.mkdir(os.path.basename(targetfile))
            #except:
                #pass
            #try:
        
                #df = {dtype:readict[k][dtype]['refvalues']}  # making a 1 column dataframe
                #write_dict_h5(targetfile, df, k, {dtype: {'compression': 32015, 'compression_opts':(3,)}}, 
                              #var_selection=[], mode='a', attrs = {dtype:readict[k][dtype]['t']['attribs']} )  
            #except Exception as e:
                #print(e, 'no values from ',k,'for station ',os.path.basename(fn))


    out_name = wpathy + fn.split('/')[-1]  
    cyear = int(out_name.split('/')[-2])
    outdir = os.path.dirname(out_name) + '/'
    with h5py.File(fn,'r') as f:
        try:
            if f['recordtimestamp'].shape[0] == 1:                
                ts=f['recordtimestamp'][:]
            else:
                ts=f['recordtimestamp'][[0, -1]]
            #if f['observations_table']['date_time'].shape[0] == 1:
                
                #ts=f['observations_table']['date_time'][:]
            #else:
                #ts=f['observations_table']['date_time'][[0, -1]]
            ref = datetime(1900, 1, 1)
            tstart = int((datetime(cyear, 1, 1) - ref).total_seconds())
            tstop = int((datetime(cyear+1, 1, 1) - ref).total_seconds())
            if tstop < ts[0] or tstart > ts[-1]:
#                print(fn, cyear, 'year missing in obs records')
                return

            ts=f['recordtimestamp'][:]
            ri = f['recordindex'][:]
            sta, sto = np.searchsorted(ts, (tstart, tstop))
            i = 0
            if sto == ri.shape[0]:
                i = 1
                tslice =slice(ri[sta], f['observations_table']['date_time'].shape[0])
            else:
                tslice =slice(ri[sta], ri[sto])
            #tsliceold = slice(*np.searchsorted(f['observations_table']['date_time'][:], (tstart, tstop)))
            #print (tslice, tsliceold)
            #assert tslice == tsliceold
            #tslice =slice(*ri[np.searchsorted(ts, (tstart, tstop))])
            
            #ts=f['observations_table']['date_time'][:]
            #tslice =slice(*np.searchsorted(ts, (tstart, tstop)))

            if tslice.stop ==tslice.start:
                return
            else:
                ts = f['observations_table']['date_time'][tslice]
                if sto == ri.shape[0]:
                    ri = np.concatenate((f['recordindex'][sta:sto], [f['observations_table']['date_time'].shape[0]]))
                else:
                    ri = f['recordindex'][sta:sto + 1]
                if len(ri) == 1:
                    print('x')
        except Exception as e:
            print(fn, cyear, e)
            return 
    #tslice = readict['tslice']
    with eua.CDMDataset(fn) as data:
        keys = data.observations_table.keys()
        keys = [x for x in keys if not x.startswith('string')]
        #keys = [x for x in keys if x in ['conversion_flag','conversion_method','report_id','date_time','observation_value','observed_variable','observation_id','z_coordinate','z_coordinate_type']]
        if 'index' in keys: 
            keys.remove('index')
        if 'shape' in keys:
            keys.remove('shape')
        obskeys = keys

        ofb=True
        try:           
            keys = data.era5fb.keys()
            #keys = [x for x in keys if x in ['fg_depar@body','an_depar@body','biascorr@body','biascorr_fg@body']]
            keys = [x for x in keys if not x.startswith('string')]
            keys.remove('index')
            if 'shape' in keys:
                keys.remove('shape')
            fbkeys = keys
        except Exception as e:
            fbkeys = ['an_depar@body', 'fg_depar@body','biascorr@body','biascorr_fg@body']
            print(e)
            ofb=False


        # loading data:
        loaded_data=[]
        a_loaded_data=[]
        loaded_type = {'names':[],'formats':[]}
        ld=[]
        addmem=1000000
        for o in obskeys:
            index = data.observations_table.index.shape[0]
            if o in ['observed_variable','observation_value','observation_id','z_coordinate','z_coordinate_type','date_time','conversion_flag','conversion_method']:  
                if len(data.observations_table[o].shape)==1:
                    loaded_data.append(data.observations_table[o][tslice])
                    a_loaded_data.append(np.empty_like(loaded_data[-1],shape=2*len(loaded_data[-1])+addmem))
                else:
                    #loaded_data.append(data.observations_table[o][:index, :].view('S{}'.format(data.observations_table[o].shape[1])).flatten())
                    loaded_data.append(data.observations_table[o][tslice, :].view('S{}'.format(data.observations_table[o].shape[1])).flatten())
#                     a_loaded_data.append(np.empty((2*len(loaded_data[-1]),len(loaded_data[-1][0])), dtype=loaded_data[-1][0].dtype))
                    a_loaded_data.append(np.empty_like(loaded_data[-1],shape=2*len(loaded_data[-1])+addmem))
                    a_loaded_data[-1].fill(b' '*data.observations_table[o].shape[1])
                loaded_type['names'].append(o)
                loaded_type['formats'].append(loaded_data[-1].dtype)
                ld.append((o,loaded_data[-1].dtype))
        try:

            loaded_obstab = np.rec.fromarrays(loaded_data, dtype=ld)
        except Exception as e:

            for l in loaded_data:
                print(l.shape)
            print(fn, e)
            return None

        del loaded_data
        a_loaded_obstab = np.rec.fromarrays(a_loaded_data, dtype=ld)
        del a_loaded_data

        loaded_fb=[]
        a_loaded_fb=[]
        loaded_type = {'names':[],'formats':[]}
        lf=[]
        for o in fbkeys:
            if o in ['fg_depar@body','an_depar@body','biascorr@body','biascorr_fg@body']:
                if ofb:
                    loaded_fb.append((data.era5fb[o][tslice]))
                else:
                    try:
                        
                        loaded_fb.append(np.full_like(data.observations_table['observation_value'][tslice], np.nan))
                    except Exception as e:
                        print(os.path.basename(fn), e)
                        return
                    
                a_loaded_fb.append(np.zeros_like(loaded_fb[-1],shape=2*len(loaded_fb[-1])+addmem))
                if o in  ['fg_depar@body','an_depar@body']:
                    a_loaded_fb[-1].fill(np.nan)
                loaded_type['names'].append(o)
                loaded_type['formats'].append(loaded_fb[-1].dtype)
                lf.append((o,loaded_fb[-1].dtype))
        loaded_feedback = np.rec.fromarrays(loaded_fb, dtype=lf)
        del loaded_fb
        a_loaded_feedback = np.rec.fromarrays(a_loaded_fb, dtype=lf)
        del a_loaded_fb
            

        humvar=np.array((ipar[34],ipar[36],ipar[38],ipar[39])) #dpd,dp,rh,sh
        wvar=np.array((ipar[104],ipar[105],ipar[106],ipar[107])) #dpd,dp,rh,sh
        hvar = np.array([117]) # geopotential (if pressure is not present)
        ps = np.array((10., 20., 30., 50., 70., 100, 150, 200, 250, 300, 400, 500, 700, 850, 925, 1000)) *100.
    
        print('before augment', time.time()-tt)
        tt = time.time()
        
        tup =augment2(loaded_obstab, a_loaded_obstab, loaded_feedback, a_loaded_feedback,
                                                  ri,ts, ps, humvar,wvar, hvar, fn)
        
        print(tup[2])#print(time.time()-tt)
        #tup =augment1(loaded_obstab, a_loaded_obstab, loaded_feedback, a_loaded_feedback,
                                                  #ri,ts, ps, humvar,wvar, hvar, fn)
        print('after augment',time.time()-tt)
    
        out, fb_out, jj,addedvar = tup
        print(np.sum((out['observed_variable']==106)))
        
        allvar = np.concatenate((np.array([117, 126]), humvar , wvar))
        #plt.figure(figsize=(10, 12))
        #l = 0
        #for v in allvar:
            #idx = np.where((out['observed_variable']==v) & (out['z_coordinate']==70000))
            #if len(idx[0]) > 0:
                
                #l += 1
                #plt.subplot(5, 2, l)
                #plt.plot(out['date_time'][idx]/86400/365.25, out['observation_value'][idx], label=str(v))
                #idy = np.where((loaded_obstab['observed_variable']==v) & (loaded_obstab['observation_value']!=-999.) & (loaded_obstab['z_coordinate']==70000))
                #if len(idy[0]) > 0:
                    #plt.plot(loaded_obstab['date_time'][idy]/86400/365.25, loaded_obstab['observation_value'][idy])
                
                #plt.legend()
        #plt.show()
        #plt.figure(figsize=(10, 12))
        #l = 0
        #for v in allvar:
            #idx = np.where((out['observed_variable']==v) & (out['z_coordinate']==70000))
            #if len(idx[0]) > 0:
                
                #l += 1
                #plt.subplot(5, 2, l)
                #plt.plot(out['date_time'][idx]/86400/365.25, fb_out['fg_depar@body'][idx], label=str(v))
                #idy = np.where((loaded_obstab['observed_variable']==v) & (loaded_obstab['observation_value']!=-999.) & (loaded_obstab['z_coordinate']==70000))
                #if len(idy[0]) > 0:
                    #plt.plot(loaded_obstab['date_time'][idy]/86400/365.25, loaded_feedback['fg_depar@body'][idy])
                
                #plt.legend()
        #plt.show()
        
        try:
            reaname = os.path.expandvars(wpathy+fn.split('/')[-1].split('_CEUAS_merged_v1.nc')[0]+'_'+str(cyear) + '.pkl')
            wtime = os.path.getmtime(reaname)
            rtime = os.path.getmtime(fn)
            if rtime>wtime:
                print(reaname, 'older than some input')
                raise ValueError
            with open(reaname,'rb') as f:
                readict=pickle.load(f)
                if 'tslice' not in readict.keys():
                    print(reaname, 'something is wrong with readict, recreating...')
                    raise ValueError
        except:
    
            out_name = wpathy + fn.split('/')[-1]  
    
            path_to_gridded=os.path.expandvars(rscratch+'/era5/gridded/')
            readict=retrieve_anfg(data.file, out,fb_out, readict,ts, tslice, refs, gactor, out_name,path_to_gridded)
    
            if readict is None:
                print(fn.split('/')[-1], 'no data for year', cyear, time.time() - tt)
                return
            with open(reaname,'wb') as f:
                pickle.dump(readict,f)
        print ('read fb',time.time()-tt)
    
        if(ofb):
            if(np.any(loaded_feedback['fg_depar@body']>1.e26) or np.any(readict['era5']['fc']['refvalues']>1.e26)):
                print(cyear)
            add_fb(out,fb_out,readict['20CRv3']['an']['refvalues'],
                   readict['era5']['an']['refvalues'],readict['era5']['fc']['refvalues'])
    
        #plt.figure(figsize=(10, 12))
        #l = 0
        #for v in allvar:
            #idx = np.where((out['observed_variable']==v) & (out['z_coordinate']==70000))
            #if len(idx[0]) > 0:
                
                #l += 1
                #plt.subplot(5, 2, l)
                #plt.plot(out['date_time'][idx]/86400/365.25, fb_out['fg_depar@body'][idx], label=str(v))
                #idy = np.where((loaded_obstab['observed_variable']==v) & (loaded_obstab['observation_value']!=-999.) & (loaded_obstab['z_coordinate']==70000))
                #if len(idy[0]) > 0:
                    #plt.plot(loaded_obstab['date_time'][idy]/86400/365.25, loaded_feedback['fg_depar@body'][idy])
                
                #plt.legend()
        #plt.show()

        del readict
        try:
            
            dr = data.recordindex[:]
            rslice = slice(*np.searchsorted(dr, (tslice.start, tslice.stop)))
            recordindex = dr[rslice]
            del dr
        except:
            print('something is wrong with recordindex in', fn)
            return
    
            # --->


    reduced_obskeys=List(loaded_obstab.dtype.fields.keys())
    reduced_fbkeys=List(loaded_feedback.dtype.fields.keys())

    for ii in range(104,108):
        #idx=np.where(loaded_obstab['observed_variable']==ii)[0]
        idy=np.where(a_loaded_obstab['observed_variable'][:jj]==ipar[ii])[0]    
        print('wind check',ipar[ii],len(idy))


    # sorting:
    print('start sorting', time.time()-tt)

    #try:      
        #os.remove(targetfile)
    #except:
        #pass
    
    # recordtimestamps are only necessary once
    rt, ri = np.unique(out['date_time'], return_index=True)
    recordtimestamps = rt
    #rt = ts[0], ts[-1]
    
    try:
        with h5py.File(fn, 'r') as file:
            with h5py.File(targetfile, 'w') as newfile:
    
                headerslice = slice(*np.searchsorted(file['header_table']['record_timestamp'][:], (rt[0], rt[-1]+1)))
                rt = file['header_table']['record_timestamp'][headerslice]
    
                groups = []
                for i in file.keys():
                    if type(file[i]) == h5py._hl.group.Group:
                        if i not in ('observations_table','era5fb'):
                            groups.append(i)
                            newfile.create_group(i)
                            
                    elif i in ('recordindex', 'recordtimestamp', 'dateindex'):
                        pass
                    else:
                        newfile.create_dataset(i, data=file[i][:], compression= 'gzip')
                for i in groups:
                    if(i == 'recordindices' or i == 'observations_table' or i == 'era5fb'):
                        pass
                    elif i in ('header_table', 'source_configuration'):
                        if 'index' not in file[i].keys():
                            newfile[i].create_dataset('index', data=np.empty(headerslice.stop-headerslice.start, dtype='S1'), compression= 'gzip')                       
                        for j in file[i].keys():
                            #print(i, j)
                            if 'string' in j:
                                xdata = file[i][j][:]
                            else:
                                xdata = file[i][j][headerslice][:]
                            newfile[i].create_dataset(j, data=xdata, compression= 'gzip')
                        if i == 'header_table':
                            newfile[i]['height_of_station_above_sea_level'][:] = station_elevation
                            
                    else:
                        if 'index' not in file[i].keys():
                            if len(file[i].keys()) == 0:
                                print('empty group', i)
                                newfile[i].create_dataset('index', data=np.empty(1, dtype='S1'), compression= 'gzip')                       
                                continue
                            sh = file[i][list(file[i].keys())[0]].shape[0]
                            newfile[i].create_dataset('index', data=np.empty(sh, dtype='S1'), compression= 'gzip')                       
                        for j in file[i].keys():
                            newfile[i].create_dataset(j, data=file[i][j][:], compression= 'gzip')
                    try:
                        
                        dim_attach(newfile, i)
                    except:
                        print(i, 'dim_attach failed')
    except MemoryError as e:
        print('could not write to resorted file')
        return

    obsv = out['observed_variable'][:jj]
    
    is_sorted = lambda a: np.all(a[:-1] <= a[1:])

                
    allvars = np.sort(np.unique(obsv))
    #
    #
    # resorting the data
    #
    @njit(cache=True, boundscheck=False)
    def make_vrindex(vridx,ridx): # this function is similar to np.unique with return_index=True, but it expands the index to the  
        # original array dimensions
        l=0
        for i in range(len(ridx)): # to set the recordindices
            if i == 0:
                l +=1
            else:
                if ridx[i]>ridx[i-1]:
                    vridx[ridx[i-1] + 1 : ridx[i] + 1]=l # next record after l
                    l += 1
                else:
                    l += 1
        vridx[ridx[i] + 1 :]=l #len(idx) # next record for the last element is the len of the data


    ridxall=np.zeros(obsv.shape[0],dtype=np.int64) # reverse index - index of the record index
    j=-1
    for j in range(len(ri)-1):
        ridxall[ri[j]:ri[j+1]]=j
    j+=1
    ridxall[ri[j]:]=j # for the last elemenet
    #dta=a_loaded_obstab['date_time'][:jj]
    ridx=[]
    vridx=[]
    absidx=[]
    absidx=np.empty_like(obsv)
    abscount=0
    idx = []
    for j in range(len(allvars)):
        idx.append(np.where(obsv==allvars[j])[0]) # index of all elements form certain variable j
#         print(j,len(idx),',',end='')
        vridx.append(np.zeros(ri.shape[0]+1,dtype=np.int64)) # all zeros in lenght of record index
        ridx=ridxall[idx[-1]] # ridxall where variable is j
        make_vrindex(vridx[-1],ridx)
        vridx[-1]+=abscount # abscount for stacking the recordindex

        idy = np.lexsort((out['z_coordinate'][idx[-1]], out['date_time'][idx[-1]]) )
        absidx[abscount:abscount+len(idx[-1])]=idx[-1][idy] # why copy? - to make sure it's not just the ref. - maybe ok without the cp
        abscount+=len(idx[-1])
        vridx[-1][-1]=abscount
    #absidx=np.concatenate(absidx)


    # check integrity of indices
    ref = datetime(1900, 1, 1)
    od = out['date_time'][:jj][absidx]
    oov = out['observed_variable'][:jj][absidx]
    ooval = out['observation_value'][:jj][absidx]
    zo = out['z_coordinate'][:jj][absidx]
    zot = out['z_coordinate_type'][:jj][absidx]
    #pres = out['z_coordinate'][:jj][absidx]
    abscount=0
    for j in range(len(allvars)):
        for i in range(1, rt.shape[0]):
            
            abscount+=len(zo[vridx[j][i]:vridx[j][i + 1]])
            if allvars[j] != 0 and np.any(zot[vridx[j][i]:vridx[j][i + 1]]!=1):
                mask = ~np.isnan(zo[vridx[j][i]:vridx[j][i + 1]])
                if not is_sorted(zo[vridx[j][i]:vridx[j][i + 1]][mask]):
                    print(fn, list(zip(zo[vridx[j][i]:vridx[j][i + 1]],
                                   zot[vridx[j][i]:vridx[j][i + 1]],
                                   oov[vridx[j][i]:vridx[j][i + 1]])))
                    raise ValueError
                
            if any(od[vridx[j][i]:vridx[j][i + 1]]!=rt[i]):
                print('spurious',allvars[j], i, ref+timedelta(seconds=int(od[vridx[j][i - 1]])),
                      ref+timedelta(seconds=int(od[vridx[j][i + 1] - 1])), ref+timedelta(seconds=int(rt[ridxall[i]])) )
                print('spurious',allvars[j], i, oov[vridx[j][i]],oov[vridx[j][i + 1] - 1])
                print('spurious',allvars[j], i, ooval[vridx[j][i]:vridx[j][i + 1]])
            #else:
                #if vridx[j][i] == vridx[j][i + 1]:
                    ##print('no values', allvars[j], i)
                    #pass
                #else:
                    #print('ok',allvars[j], i, ref+timedelta(seconds=int(od[vridx[j][i]])),
                        #ref+timedelta(seconds=int(od[vridx[j][i + 1] - 1])), ref+timedelta(seconds=int(rt[i])) )
                    #print('ok', allvars[j], i, oov[vridx[j][i]],oov[vridx[j][i + 1] - 1])
                    #print('ok',allvars[j], i, ooval[vridx[j][i]:vridx[j][i + 1]])

    del ridxall

    print('elapsed converting: ',time.time()-tt)

    for i in obskeys:

        if i == 'observation_id':
            ov_vars = np.empty(jj,dtype=out[i].dtype)
            ov_vars[:]=out[i][:jj]
            dim1 = int(str(out[i].dtype).split('S')[-1])
            ov_vars=fill_obsid(ov_vars.view('S1').reshape((len(ov_vars),dim1)),out['conversion_flag'][:jj])

        elif i in reduced_obskeys:
            ov_vars = out[i][:jj]
        else:
            continue

        ov_vars = ov_vars[absidx]
        if i == 'index':
            pass
        elif i in ['observation_id', 'report_id', 'sensor_id', 'source_id']:
            alldict = {i:np.asarray(ov_vars, dtype='S1')}
            write_dict_h5(targetfile, alldict, 'observations_table', {i: { 'compression': 32015, 'compression_opts':(3,) } }, [i])
        else:
            alldict = {i:ov_vars}
            write_dict_h5(targetfile, alldict, 'observations_table', {i: { 'compression': 32015, 'compression_opts':(3,) } }, [i])  
        #print(i, time.time()-tt)
    del ov_vars
    del obsv
    #del dta
    outsh = out.shape[0]
    del out
    del a_loaded_obstab
    gc.collect()

    if True:
        for i in fbkeys:
            if i in reduced_fbkeys:
                ov_vars = fb_out[i][:jj]
            else: 
                continue

            ov_vars = ov_vars[absidx]

            if i == 'index':
                pass
            elif i in ['expver', 'source@hdr', 'source_id', 'statid@hdr']:
                alldict = {i:np.asarray(ov_vars, dtype='S1')}
                write_dict_h5(targetfile, alldict, 'era5fb', {i: { 'compression': 32015, 'compression_opts':(3,) } }, [i])
            else:
                #alldict = pandas.DataFrame({i:ov_vars})
                alldict = {i:ov_vars}
                
                write_dict_h5(targetfile, alldict, 'era5fb', {i: { 'compression': 32015, 'compression_opts':(3,) } }, [i]) 

            #print(i, time.time()-tt)
        del fb_out
        del a_loaded_feedback
        gc.collect()

        obskeys.append('station_elevation') # fix for adding station elevation 
        #if True:
        with eua.CDMDataset(fn) as data:
            for i in obskeys:
                if i in reduced_obskeys or i == 'observation_id' :
                    continue
                #if i not in ('latitude','longitude','report_id'):
                    #continue
                if i != 'station_elevation':                   
                        #i='z_coordinate'
                    rest_data = data.observations_table[i][tslice]
                    if rest_data.ndim==2: #i in ['observation_id', 'report_id', 'sensor_id', 'source_id']:
                        ov_vars = np.empty((addedvar.shape[0],len(rest_data[0])), dtype=rest_data[0].dtype)
                        ov_vars.fill('n')
                    else:
                        ov_vars = np.empty(addedvar.shape[0], dtype=rest_data.dtype)
                        try:
                            
                            ov_vars.fill(np.nan)
                        except:
                            ov_vars.fill(-2147483648)
                            
                    
    
                #print('vor fill_restdata', fn)

                    fill_restdata(ov_vars, rest_data, addedvar) #,out['z_coordinate'][:jj])
                    ov_vars = ov_vars[absidx]

                else:
                    ov_vars = np.full(jj, station_elevation,dtype=np.float32) # add station elevation
                    
                if i == 'index':
                    pass
                elif i in ['observation_id', 'report_id', 'sensor_id', 'source_id']:
                    alldict = {i:np.asarray(ov_vars, dtype='S1')}
                    write_dict_h5(targetfile, alldict, 'observations_table', {i: { 'compression': 32015, 'compression_opts':(3,) } }, [i])
                else:
                    #alldict = pandas.DataFrame({i:ov_vars})
                    alldict = {i:ov_vars}
                    write_dict_h5(targetfile, alldict, 'observations_table', {i: { 'compression': 32015, 'compression_opts':(3,) } }, [i])  

                #print(i, time.time()-tt)
        #if True:
        with eua.CDMDataset(fn) as data:
            for i in fbkeys:

                if i in reduced_fbkeys:
                    continue
                else: 
                    rest_data = data.era5fb[i][tslice]
                    if i in ['expver', 'source@hdr', 'source_id', 'statid@hdr']:
                        ov_vars = np.empty((addedvar.shape[0],len(rest_data[0])), dtype=rest_data[0].dtype)
                        ov_vars.fill('n')
                    else:
                        ov_vars = np.empty(addedvar.shape[0], dtype=rest_data[0].dtype)
                        try:
                            
                            ov_vars.fill(np.nan)
                        except:
                            ov_vars.fill(-2147483648)
                            
                ov_vars = fill_restdata(ov_vars, rest_data, addedvar)
                ov_vars = ov_vars[absidx]

                if i == 'index':
                    pass
                elif i in ['expver', 'source@hdr', 'source_id', 'statid@hdr']:
                    alldict = {i:np.asarray(ov_vars, dtype='S1')}
                    write_dict_h5(targetfile, alldict, 'era5fb', {i: { 'compression': 32015, 'compression_opts':(3,) } }, [i])
                else:
                    alldict = {i:ov_vars}
                    write_dict_h5(targetfile, alldict, 'era5fb', {i: { 'compression': 32015, 'compression_opts':(3,) } }, [i]) 
                #print(i, time.time()-tt)
        #
        # writing the recordindices and recordtimestamp.
        #       
    recordindices=vridx
    for i in range(len(recordindices)):
        testvar = {str(allvars[i]):recordindices[i]}
        write_dict_h5(targetfile, testvar, 'recordindices', {str(allvars[i]): { 'compression': None } }, [str(allvars[i])]) 

    write_dict_h5(targetfile, {'recordtimestamp':recordtimestamps}, 'recordindices', {'recordtimestamp': { 'compression': None } }, ['recordtimestamp'])

    print('elapsed writing '+targetfile+':',time.time()-tt)
    try:
        os.mkdir(wpathy+'/log')
    except:
        pass
    f= open(wpathy+'/log/'+fn.split('/')[-1]+".txt","w")
    f.write("done") 
    f.close()
    return targetfile

ray_convert_missing = ray.remote(convert_missing)

def rmsplot(files):
    
    readicts = []
    year = os.path.dirname(files[0]).split('/')[-1]
    pkname = files[0] +'_'+year + '.pkl'
    fno = wpath+'/plots/'+os.path.basename(pkname).split('.')[0].split('_')[0]
    if os.path.exists(fno+'.txtx'):
        print(fno+' already plotted')
        return
    
    print('rmsplot', files)
    for file in files:
        year = os.path.dirname(file).split('/')[-1]
        pkname = file.split('_CEUAS')[0] +'_'+year + '.pkl'
        
        print(pkname)
        try:
            
            with open(pkname,'rb') as f:
                readicts.append(pickle.load(f))
        except Exception as e:
            print (pkname + ' not found', e)
    if len(readicts) == 0:
        print(pkname+' no feedback information found')
        return
    
    
    readict = {'obs' : {}}
    print(readicts[0].keys())
    for k in readicts[0]['obs'].keys():
        readict['obs'][k] = np.concatenate([i['obs'][k] for i in readicts])
    for p in readicts[0]['obstypes'].keys():

        for k,v in readicts[0].items():
            if k in ('obs', 'tslice'):
                continue
            if 'obs' in k:
                readict[k] = readicts[0][k]
                continue
            if k =='era5fb':
                i = 0
            for ftype in v['ftype']:  
                if 'fc' in ftype:
                    dtype='fc'
                else:
                    dtype='an'
                if k not in readict.keys():
                    readict[k] = {}
                if dtype not in readict[k].keys():
                    readict[k][dtype] = {}
                readict[k][dtype]['refvalues'] = np.concatenate([i[k][dtype]['refvalues'] for i in readicts])

    try:
        os.mkdir(wpath+'/plots')
    except:
        pass
    
    readictplot(readict, ('mean', 'std','rms'), (10000, 30000,70000,85000,92500 ), fno, marker='')
    with open(fno+'.txt', 'w') as f:
        f.write('plots done')

    return

ray_rmsplot = ray.remote(rmsplot)

def do_station(refs_ref, file, years):
    
    futures = [ray_convert_missing.remote(refs_ref, wpath, year, file ) for year in years]
    

    return ray.get(ray_rmsplot.remote(futures))

ray_do_station = ray.remote(do_station)

def readgridded(readict, refs_ref, year, month):
    fdict = {}
    for k,v in readict.items():
        fdict[k] = {}
        if k!='20CRv3':
            continue
        for ftype in v['ftype']:  
            if 'fc' in ftype:
                dtype='fc'
            else:
                dtype='an'
            if dtype not in readict[k].keys():
                readict[k][dtype]={}
                #readict[k][dtype]['refvalues']=np.empty(obs.shape,dtype=np.float32)
                #readict[k][dtype]['refvalues'].fill(np.nan)
            fdict[k][dtype] = {}
            for p,pn in v['param'].items():
                if k!='JRA55':
        
                    fpattern=v['path']+v['prefix']+ftype+v['glue']+'{}{:0>2}'+v['glue']+pn+v['suffix']+'.nc'
                else:
                    fpattern=v['path']+v['prefix']+ftype+v['glue']+pn+'.reg_tl319.{}{:0>2}'+v['glue']+v['suffix']+'.nc'
            
                if k != '20CRv3':
                    fn = fpattern.format(year,month)
                    zname = 'level'
                else:
                    fn =  v['path']+v['prefix']+ftype+v['glue']+'{}'.format(year)+v['glue']+pn+v['suffix']+'.nc'
                    zname = 'isobaricInhPa'
                    
                
                with h5py.File(fn,'r') as g:
                    
                    if p not in fdict[k][dtype].keys():
                        if k != '20CRv3':
                            fdict['level'] = g['level'][:]
                            fdict['pidx'] = np.searchsorted(fdict['level'],refs['level'])
                        else:
                            fdict['level'] = g['isobaricInhPa'][:]                                  
                            fdict['pidx'] = np.searchsorted(-fdict['level'], -refs['level'])
                        
                        fdict['longitude']= g['longitude'][:]
                        fdict['latitude']= g['latitude'][:]
                        
                        fdict['attribs'] =  dict(g[p].attrs)
                        try:
                            
                            del fdict['attribs']['DIMENSION_LIST']
                        except:
                            pass
                        for kk,vv in g.attrs.items():
                            fdict['attribs'][kk]=vv
                        fdict['attribs']['source_filename']=fn
                        fdict[k][dtype][p] = g[p][:][:, fdict['pidx'], :, :]* np.float32(fdict['attribs']['scale_factor'])+np.float32(fdict['attribs']['add_offset'])
                        
                        print(k, dtype, p)
                        tunits=g['time'].attrs['units']
                        try:
                            tunits=tunits.decode('latin1')
                        except:
                            pass
                        fdict['tunits']=tunits.split()
            
        
    return fdict
        #found=len(glob.glob(fpattern.format(*yms[0])))
        #print(fpattern.format(*yms[0]),found)
        #if found:
    fdict = {}
    func=partial(offline_fb3,fpattern,p,obstypes[p], fdict,ts,obstype, lat,lon, z, refs)

ray_readgridded =ray.remote(readgridded)   

def plot_contents(wpath,cyear, fn):

    pardict = {'0': ['unknown',''], 
               '34': ['dewpoint depression','K'],
               '39': ['specific humidity','kg/kg'],
               '139': ['u-component of wind','m/s'],
               '140': ['v-component of wind','m/s'],
               '106': ['wind from direction','deg'],
               '107': ['wind speed','m/s'],
               '117': ['geopotential','J/kg'],
               '126': ['temperature','K'],
               '137': ['dewpoint','K'],
               '138': ['relative humidity','']
               }
    print(fn.split('/')[-1], cyear,end=' ' )
    wpathy = wpath+'/' + str(cyear) + '/'
    targetfile = wpathy + fn.split('/')[-1]
    #if os.path.isfile(wpathy+'/log/'+fn.split('/')[-1]+".txt"):
        #print('already processed')
        #return targetfile
    #else:
        #print('processing...')
    sys.path.insert(0,os.getcwd()+'/../resort/rasotools-master/')
    
    with h5py.File(targetfile, 'r') as f:
        rk = list(f['recordindices'].keys())[:-2]
        for pl in [70000.]:
            plt.figure(figsize=(10, (len(rk) +1) *1.5))
            l = 1
            for v in rk:
                tslice = slice(*f['recordindices'][v][[0, -1]])
                ts = f['observations_table']['date_time'][:][tslice]
                obs = f['observations_table']['observation_value'][tslice]
                pidx = np.where(f['observations_table']['z_coordinate'][tslice]==pl)[0]
                if len(pidx) < 2:
                    continue
                plt.subplot( len(rk),1, l)
                plt.plot(np.asarray(1900+ts[pidx]/86400/365.25, dtype='int'), rmeanw(obs[pidx], 30))
                #plt.plot(np.asarray(1900+thin2(ts[pidx], 10)/86400/365.25, dtype='int'), thin2(rmeanw(obs[pidx], 30), 10))
                plt.title(os.path.basename(targetfile).split('_')[0]+','+str(np.int(pl/100.)) + 'hPa, '+ pardict[v][0])
                plt.ylabel(pardict[v][1])
                l += 1
        plt.tight_layout()
        plt.savefig(targetfile[:-3]+'.png')
        plt.close()
                
                
        
################################################################
        
if __name__ == '__main__':


    wpath= os.path.expandvars('$RSCRATCH/converted_v13/') #'./'
    opath=wpath
    wlpath=wpath+'log/'
    for p in wpath,wlpath, wpath + '/long', wpath + '/rea':      
        try:
            os.mkdir(p)
        except:
            pass

    ipar=np.zeros(140,dtype=np.int32)-2100000000 # 
    ipar[0]=0
    ipar[34]=34
    ipar[39]=39
    ipar[85]=126
    ipar[106]=106
    ipar[107]=107
    ipar[117]=117
    #ipar[]=136
    ipar[36]=137 #dp
    ipar[38]=138 #rh
    ipar[104]=139
    ipar[105]=140

    readict={'era5':{'ftype':('t','fct'),'param':{'t':'130','u':'131','v':'132','q':'133','z':'129'},
                     'path':os.path.expandvars('$RSCRATCH/era5/gridded/'),'prefix':'era5','suffix':'','glue':'.'},
             #'CERA20C':{'ftype':('t',),'param':{'t':'130','u':'131','v':'132','q':'133','z':'129'},
               #'path':os.path.expandvars('$RSCRATCH/CERA20C/'),'prefix':'CERA20C','suffix':'','glue':'.'},
                          #'JRA55':{'ftype':('fcst_mdl',),'param':{'t':'011_tmp','u':'033_ugrd','v':'034_vgrd','q':'051_spfh'},
               #'path':os.path.expandvars('$RSCRATCH/JRA55/split/'),'prefix':'test.','suffix':'grb','glue':'.'},
                        '20CRv3':{'ftype':('',),'param':{'t':'TMP','u':'UGRD','v':'VGRD'},#,'q':'SPFH'},#,'z':'HGT'},
                         'path':os.path.expandvars('$RSCRATCH/20CRv3/'),'prefix':'anl_meant','suffix':'_pres','glue':'_'},
               }
    readict={'era5':{'ftype':('','fc'),'param':{'t':'130','u':'131','v':'132','q':'133','z':'129'},
                     'path':os.path.expandvars('$RSCRATCH/era5/gridded/'),'prefix':'era5','suffix':'','glue':'.'},
             #'CERA20C':{'ftype':('t',),'param':{'t':'130','u':'131','v':'132','q':'133','z':'129'},
               #'path':os.path.expandvars('$RSCRATCH/CERA20C/'),'prefix':'CERA20C','suffix':'','glue':'.'},
                          #'JRA55':{'ftype':('fcst_mdl',),'param':{'t':'011_tmp','u':'033_ugrd','v':'034_vgrd','q':'051_spfh'},
               #'path':os.path.expandvars('$RSCRATCH/JRA55/split/'),'prefix':'test.','suffix':'grb','glue':'.'},
                        '20CRv3':{'ftype':('',),'param':{'t':'TMP','u':'UGRD','v':'VGRD'},#,'q':'SPFH'},#,'z':'HGT'},
                         'path':os.path.expandvars('$RSCRATCH/20CRv3/'),'prefix':'anl_meant','suffix':'_pres','glue':'_'},
               }
    
    #files = glob.glob('/mnt/scratch/scratch/federico/YEAR_SPLIT_MERGING/[12]???/0-*v2.nc')
    
    tt = time.time()
    files = glob.glob('/mnt/scratch/scratch/federico/MERGED_FEB2023/*v1.nc')[:] #+ \
    files += glob.glob('/mnt/scratch/scratch/federico/MERGED_FEB2023_MOBILE_mobile/*v1.nc')[:] + \
        glob.glob('/mnt/scratch/scratch/federico/MERGED_FEB2023_coordinates_coordinate/*v1.nc')[:] #MERGED_FEB2023_big
    #files= glob.glob('/scratch/das/federico/MERGED_MOBILE_FIX_STATION_CONFIGURATION/*v1.nc')[:]
    #files +=glob.glob('/mnt/scratch/scratch/federico/MERGED_02APR2023_arctic_giub_NPSOUND_SHIPSOUND_coordinate/*v1.nc')[:] #MERGED_FEB2023_big
    #files +=glob.glob('/mnt/scratch/scratch/federico/MERGED_FEB2023_IGRAship_20300-0-99/*v1.nc')[:]
    #files += glob.glob('/mnt/scratch/scratch/federico/MERGED_02APR2023_arctic_giub/*v1.nc')[:]

    #files = ['/mnt/scratch/scratch/federico/MERGED_02APR2023_arctic_giub/0-20000-0-20891_CEUAS_merged_v1.nc']
    #files = ['/mnt/scratch/scratch/federico/MERGED_02APR2023_arctic_giub/0-20000-0-01107_CEUAS_merged_v1.nc']
    #files = ['/scratch/das/federico/COP2_HARVEST_MAR2023_giub_arctic_01APR2023//giub/0-20000-0-20891_giub_harvested_4777.txt_converted_csv.csv.nc']
    #files = ['/mnt/scratch/scratch/federico/MERGED_FEB2023/0-20000-0-01107_CEUAS_merged_v1.nc']
    #files = ['/mnt/scratch/scratch/federico/MERGED_FEB2023/0-124-0-73033_CEUAS_merged_v1.nc']
#    files = ['/mnt/scratch/scratch/federico/MERGED_FEB2023_MOBILE_mobile/0-20999-0-ELML_CEUAS_merged_v1.nc']
    

        
    print(len(files), time.time() - tt)

    print(files[0], files[-1])

    #already_done = glob.glob(wlpath+'*_yearly.txt')

    #files_to_convert = []#files #[]
    #for i in files:
        #if not wlpath+i.split('/')[-1]+'_yearly.txt' in already_done:
            #files_to_convert.append(i)
    files_to_convert = files
    #files_to_convert.sort()
    tt=time.time()

    
    
    refs = load_20CRoffset(readict)
    readicts = []
    #for i in range(1950, 1939, -1):
        
        #func=partial(convert_missing, refs, wpath, i)
        #result_list = list(map(func, files_to_convert[:]))
        
    ## simple plotting
    #for i in  ('long', ):#(1995, 1994, -1):
        
        #func=partial(plot_contents, wpath, i)
        #result_list = list(map(func, files_to_convert[:]))

    #exit()
               

    
    #fp = {}
    #for fn in glob.glob(os.path.expandvars('$RSCRATCH/era5/gridded/era5t.*')):
        #fp[os.path.basename(fn)] = h5py.File(fn, 'r')
        
    #ray_fp = ray.put(fp)
        
    ray.init(num_cpus=60, _temp_dir=os.path.expanduser('~/ray'))
    refs_ref = ray.put(refs)

    futures = []
    for year in range(2020, 1904, -1):
        """        glist = glob.glob(os.path.expandvars('$RSCRATCH/era5/gridded/era5t.{}??.*.nc'.format(year))) + \
            glob.glob(os.path.expandvars('$RSCRATCH/era5/gridded/era5fct.{}*.nc'.format(year))) + \
            glob.glob(os.path.expandvars('$RSCRATCH/20CRv3/anl_meant_{}*.nc'.format(year)))
        gactor = {}
        gp = {}
        futures = []
        for g in glist:
            gactor[os.path.basename(g)] = h5p.remote(g)
            #gp[os.path.basename(g)] = h5py.File(g, 'r')
            #print(os.path.basename(g))
            futures.append(gactor[os.path.basename(g)].name.remote())
        x = ray.get(futures)
        myactors =ray.put(gactor)
        #convert_missing(refs_ref, myactors, wpath, year, files_to_convert[0]) 
        futures = [ ray_convert_missing.remote(refs_ref, myactors, wpath, year, file )  for file in files_to_convert[:]]
        """
        for file in files_to_convert:
            if '0-11010' in file: #or '0-07602_CEUAS_merged_v1' in file:
                #convert_missing(refs, None, wpath, year, file)
                futures.append(ray_convert_missing.remote(refs_ref, None, wpath, year, file ))
        
        #futures = futures + [ ray_convert_missing.remote(refs_ref, None, wpath, year, file )  for file in files_to_convert]
    obj_result_list = ray.get(futures)
    exit()
        
    """
        del myactors
        for g in glist:
            del gactor[os.path.basename(g)]
    """

    print('')
    #exit()
    if False:
        pwd = os.getcwd()
        obj_result_list = glob.glob(wpath+'/[12]???/*v1.nc')
        os.chdir(pwd)
        shlist = np.unique([os.path.basename(ll) for ll in obj_result_list if ll is not None])
        futures = []
        vienna = False
        for file in shlist:
            if '11035' in file:
                vienna = True
            sublist = sorted([s for s in obj_result_list if s is not None and file in s])
            #rmsplot(sublist)
            #if True or vienna:
            futures.append(ray_rmsplot.remote(sublist))
        ray.get(futures)

    pwd = os.getcwd()
    os.chdir(wpath)
    longlist = glob.glob('[12]???/0-20999*v1.nc') + glob.glob('[12]???/*43185*v1.nc')
    os.chdir(pwd)
    shlist = []
    shlist = [os.path.basename(ll) for ll in longlist]
    fkeys = np.unique(shlist)
    
    futures =[]   
    for fks in fkeys:
        #if '0-62053' not in fks:
            #continue
        fkey = wpath+'/[12]???/' + fks
                
        #h5concatenate(fkey)
        futures .append(ray_h5concatenate.remote(fkey))
        
    ray.get(futures)
    
    print('total:',time.time()-tt)
    print("Don't for get to run orphan script")
    
    os.chdir(wpath+'/long/')
    fns=glob.glob('*orphan*.nc')
    fns.sort()
    i=0
    try:
        os.mkdir ('orph')
    except:
        pass
    for fn in fns:
        i+=1
        wigos=f'0-23001-2-orph{i:04d}'
        fno='orph/'+wigos+'_CEUAS_merged_v1.nc'
        wigos=np.string_(wigos)
        shutil.copyfile(fn,fno)
        f=h5py.File(fno,'r+')
        #print(f['station_configuration/primary_id'][:])
        print(f['station_configuration/primary_id'].shape)
        l=f['station_configuration/index'].shape[0]
        del f['station_configuration/primary_id']
        f['station_configuration'].create_dataset('primary_id',data=(np.vstack([wigos]*l)).view('S1'))
        
        l=f['header_table/index'].shape[0]
        del f['header_table/primary_station_id']
        f['header_table'].create_dataset('primary_station_id',data=(np.vstack([wigos]*l)).view('S1'))
        print(fn,fno)
     
    
    
    l = 0

