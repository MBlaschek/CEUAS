#!/usr/bin/env
# coding: utf-8

import pandas as pd
#import hdf5plugin
from numba import njit
import numpy as np

import sys,glob
import zipfile, os, time
import urllib3
from datetime import datetime, timedelta
os.chdir(os.path.expandvars("$HOME/python/CEUAS/CEUAS/public/resort/"))
x=pd.read_csv('../../meta/inventory_comparison_2/code/station_configuration/CUON_station_configuration_extended.csv',delimiter='\t')
import glob, shutil
import h5py
sys.path.append(os.getcwd()+'/../cds-backend/code/')
sys.path.append(os.getcwd()+'/../harvest/code/')
from harvest_convert_to_netCDF import write_dict_h5
from convert_numbamod import augment2,  p_to_z_ifs, triint,trisimple, pfill,pfilla, \
     dewpoint_Sonntag, frostpoint_Sonntag, sh2vap, svp , fill_obsid, extract_tuv, \
     fill_restdata, fill_restdata_abs, fill_restdata_addedvar, fill_restdata_addedvar_abs, fill_restdata_obsid_abs, \
     liquid, bisection, func, rms, nancheck, check_spurious, is_surf, is_surf2, is_surf3
import gc
#import cds_eua3 as eua
import cds_eua4 as eua
#eua.logging_set_level(30)
#import xarray as xr
import cdsapi, zipfile, os, time
import copy
from shutil import copyfile
#from multiprocessing import Pool
sys.path.insert(0,os.getcwd()+'/../resort/rasotools-master/')
import rasotools
dir(rasotools)
import warnings
from functools import partial
import matplotlib.pylab as plt
warnings.filterwarnings('ignore')
import pickle
import ray
import rs_drift as rsd
from tqdm import tqdm


#for fn in glob.glob('/mnt/scratch/scratch/federico/MERGED_YEARLY_22NOV_FULL_checkDimensions/*/*v3.nc'):
        #os.makedirs(os.path.dirname('/mnt/users/scratch/leo/scratch/test'), exist_ok=True)
        #target = '/mnt/users/scratch/leo/scratch/test/' + '/'.join(fn.split('/')[-2:])
        #os.makedirs(os.path.dirname(target), exist_ok=True)
        #shutil.copyfile(fn, target)
        #with h5py.File(target, 'r+') as f:
                #for gr in 'observations_table', 'era5fb', 'header_table':
                        #for vs in f[gr].keys():
                                #if vs != 'index' and 'string' not in vs:
                                        #var = f[gr][vs][:]
                                        #del f[gr][vs]
                                        #write_dict_h5(f, {vs: var}, gr, {vs: { 'compression': 32015, 'compression_opts':(3,) } }, mode='a', chunksize=100000)
                                        #print(target, gr, vs)
                        

def ptime(c, reftime, debug=True, end=None):
    if debug:
        print(f'{c}: {time.time() - reftime:.3f}', end=end)

def  copy_attrs(file, new_file, i):
    for v in file[i].keys():
        for a in  file[i][v].attrs.keys():
            if a not in ('_Netcdf4Dimid', 'CLASS', 'NAME', 'REFERENCE_LIST', 'DIMENSION_LIST'):
                
                new_file[i][v].attrs[a] = file[i][v].attrs[a]
    
    return
            
new = True
for wid in  []:#'0-20001-0-10393', : #0-20001-0-10393
    for vs in  'observation_id', :
        
        for chunks in 1000, 10000, 100000, 1000000:
            for filters in (32015,(3,)),('gzip', 3): #, 32004:
                if new:
                    fn = '/mnt/scratch/scratch/federico/MERGED_YEARLY_22NOV_FULL_checkDimensions/'+wid+'/'+wid+'_2022_CEUAS_merged_v3.nc'
                    with h5py.File(fn) as f:
                        tt = time.time()
                        var = f['observations_table'][vs][:]
                        print(time.time()-tt)
                        new = False
                    target = '/mnt/users/scratch/leo/scratch/test/' + '/'.join(fn.split('/')[-2:])
                    with h5py.File(target) as f:
                        tt = time.time()
                        var = f['observations_table'][vs][:]
                        print(time.time()-tt)
                        new = False
                               
                with h5py.File('/mnt/users/scratch/leo/scratch/converted_v15/2022/'+wid+'_CEUAS_merged_v3.nc') as f:
                    tt = time.time()
                    var = f['observations_table'][vs][:]
                    print(time.time()-tt)
                    tt = time.time()
                    target = f.filename.split('.nc')[0]+'_'+vs+'.nc'
                    write_dict_h5(target, {vs: var}, 'observations_table', {vs: { 'compression': filters[0], 'compression_opts':filters[1] } }, mode='w', chunksize=chunks)
                    print(time.time()-tt)
                with h5py.File('/mnt/users/scratch/leo/scratch/converted_v15/2022/'+wid+'_CEUAS_merged_v3_'+vs+'.nc') as f:
                    tt = time.time()
                    var = f['observations_table'][vs][:]
                    #target = f.filename.split('.nc')[0]+'_'+vs+'.nc'
                    #print(time.time()-tt)
                    #tt = time.time()
                    #write_dict_h5(target, {vs: var}, 'observations_table', {vs: { 'compression': filters, 'compression_opts':(3,) } }, mode='w', chunksize=chunks)
                    print(time.time()-tt, wid, vs, chunks, filters, os.path.getsize(target))


from assemble_longrecord import h5concatenate
ray_h5concatenate = ray.remote(h5concatenate)


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

@njit #(boundscheck=True)
def memest3(ts, z):
    
    #if len(z) != len(ts):raise ValueError('z!=ts')
    l = 0
    for i in range(1, len(ts)):
        if ts[i] != ts[i-1] or z[i] != z[i-1]:
            l += 1
    #tups, idx= np.unique(tuple(zip(ts, z)), return_index=True)
    
    return(l *10)
    
    
def memest2(ov, ri, z):
    
    levs = []
    for i in range(ri.shape[0]-1):
        levs.append(np.unique(z[ri[i]-ri[0]:ri[i+1]-ri[0]]).shape[0])
        
    if np.any(np.isin(ov, (ipar[34], ipar[36], ipar[38], ipar[39]))):
        fak = 12
    else:
        fak = 8
    
    return int(fak*np.mean(levs))
    
    
def memest(ov, z, ipar):
    #humvar=np.array((ipar[34],ipar[36],ipar[38],ipar[39])) #dpd,dp,rh,sh
    #wvar=np.array((ipar[104],ipar[105],ipar[106],ipar[107])) #dpd,dp,rh,sh
    #hvar = np.array([117]) # geopotential (if pressure is not present)
    ptup = (117,126, ipar[104],ipar[105],ipar[106],ipar[107], ipar[34],ipar[36],ipar[38],ipar[39])
    sdict = dict(zip(ptup,ptup ))
    
    for k in sdict.keys():
        sdict[k] = np.sum(ov==k)
    
    wmax = np.sum((sdict[ipar[104]], sdict[ipar[106]]))
    qmax = np.sum((sdict[ipar[34]], sdict[ipar[36]], sdict[ipar[38]], sdict[ipar[39]]))
    
    needed = sdict[117] + sdict[126] + 4 * wmax + 4 * qmax
    
    return needed
    
        
        


@njit(cache=False, boundscheck=True)
def add_fb(loaded_obstab,loaded_feedback,ref20CR,refera5an,refera5fc, refzidx):
    i20=0
    iera=0
    
    i = 0
    #while i < loaded_obstab['date_time'].shape[0]:
    while i < refzidx.shape[0]:
        m = refzidx[i]
        k = 0
        pold = loaded_obstab['z_coordinate'][m]
        dold = loaded_obstab['date_time'][m]
        while m + k < loaded_obstab['date_time'].shape[0] and \
              loaded_obstab['z_coordinate'][m + k] == pold and \
              loaded_obstab['date_time'][m + k] == dold:
            k += 1
        if k == 0:
            i = i + 1
            continue
    
        #imax = m + k
        #imin = m
        q = np.nan
        u = np.nan
        t = np.nan
        for j in range(k):
            if loaded_obstab['observed_variable'][m + j] == 39:
                if refera5fc[i+j] == refera5fc[i+j]:
                    q = refera5fc[i+j]
                else:
                    q = ref20CR[i+j]                    
            if loaded_obstab['observed_variable'][m + j] == 126:
                if refera5fc[i+j] == refera5fc[i+j]:
                    t = refera5fc[i+j]
                else:
                    t = ref20CR[i]+j
            if loaded_obstab['observed_variable'][m + j] == 139:
                if refera5fc[i] == refera5fc[i]:
                    u = refera5fc[i+j]
                else:
                    u = ref20CR[i+j]
            if loaded_obstab['observed_variable'][m + j] == 140:
                if refera5fc[i+j] == refera5fc[i+j]:
                    v = refera5fc[i+j]
                else:
                    v = ref20CR[i+j]
                
    
    #for j in range(imin,imax):
        for j in range(k):
            if loaded_obstab['observed_variable'][m +j] == 138:
                if q == q:                
                    p = loaded_obstab['z_coordinate'][m +j]
                    vpdata = sh2vap(q, p)
                    x = svp(np.array((t,) ), np.array((p,) ))
                    if x[0] > 0:                       
                        refera5fc[i+j] = vpdata / x[0]
                    else:
                        refera5fc[i+j] = np.nan
                
            if loaded_obstab['observed_variable'][m +j] == 137:
                if q == q:                
                    p = loaded_obstab['z_coordinate'][m +j]
                    vpdata = sh2vap(q, p)
                    dp = bisection(func, 150, 350, vpdata)
                    
                    refera5fc[i+j] = dp # dewpoint by inverting Sonntag's formula. Don't use dewpoint_Sonntag(vpdata)
                #if vpdata < 610.:
                    #refera5fc[i] = frostpoint_Sonntag(vpdata)
                #else:
                    #refera5fc[i] = dewpoint_Sonntag(vpdata)
                
                #refera5fc[j] = np.nan
    
            if loaded_obstab['observed_variable'][m +j] == 34:
                if q == q and t == t:         
                    p = loaded_obstab['z_coordinate'][m +j]
                    vpdata = sh2vap(q, p)
                    #if vpdata < 610.:
                        #dp = frostpoint_Sonntag(vpdata)
                    #else:
                        #dp = dewpoint_Sonntag(vpdata)
                    dp = bisection(func, 150, 350, vpdata)
                    refera5fc[i+j] = t - dp
            
            if loaded_obstab['observed_variable'][m +j] == 107:
                if u == u:
                    
                    refera5fc[i+j] = np.sqrt(u*u+v*v)
                    
            if loaded_obstab['observed_variable'][m +j] == 106:
                if u == u:
                    wd = - np.arctan2(v, u) * 180 / np.pi - 90.
                    if wd < 0:
                        wd+= 360.
                    refera5fc[i+j] = wd
    
                
                        
        #add offline calculated feedback
        #for j in range(k):
            if loaded_obstab['observation_value'][m + j]==loaded_obstab['observation_value'][m + j]:
                #if i < refera5fc.shape[0] and  i < ref20CR.shape[0]:
                    if refera5fc[i+j]==refera5fc[i+j]:
                        loaded_feedback['fg_depar@offline'][m + j]=loaded_obstab['observation_value'][m + j]-refera5fc[i+j]
                        #loaded_feedback['biascorr@offline'][i]=0.
                        iera+=1
                    elif ref20CR[i+j]==ref20CR[i+j]:
                        loaded_feedback['fg_depar@offline'][m + j]=loaded_obstab['observation_value'][m + j]-ref20CR[i+j]
                        #loaded_feedback['biascorr@offline'][i]=0.
                        i20+=1
                    
        i = i + k
    rawdir = loaded_feedback['fg_depar@offline'][loaded_obstab['observed_variable']==106]
    rawdir[rawdir>180.] -= 360.
    rawdir[rawdir<-180] += 360.
    loaded_feedback['fg_depar@offline'][loaded_obstab['observed_variable']==106] = rawdir
      
    return
        
        #if i%1000000==0:
            #print(i,i20,iera)

#from numba.typed import List
def make_fdict(fps, p, cyear, refs):
    
    l = 0
    g = fps[l]
    
    fdict = {}
    fdict['level'] = g['level'][:]
    fdict['pidx'] = np.searchsorted(fdict['level'],refs['level'])
    
    fdict['longitude']= g['longitude'][:]
    fdict['latitude']= g['latitude'][:]
    
    #if 'era5fc.0.25t' in fn and '133' in fn:
        #x = 0

    fdict['attribs'] =  dict(g[p].attrs)
    try:
        
        del fdict['attribs']['DIMENSION_LIST']
    except:
        pass
    #for k,v in g.attrs.items():
        #fdict['attribs'][k]=v
    fdict['attribs']['source_filename']=os.path.basename(g.filename)
    
    
    
    yplus = int((datetime(cyear+1, 1, 2) - datetime(1900, 1, 1) ).total_seconds()) # include one day of year after
        
    yminus = int((datetime(cyear, 1, 1) - datetime(1900, 1, 2) ).total_seconds()) # include one day of year before
    atup = 'fi','mti','sf', 'ao', 'secs'
    for s in atup:        
        fdict[s] = []

    l = 0
    for g in fps:
        if g is None:
            continue
        tunits=g['time'].attrs['units']
        try:
            tunits=tunits.decode('latin1')
        except:
            pass
        fdict['tunits']=tunits.split()
        offset = 0
        fak = 1
        if fdict['tunits'][0]=='hours':
            if tunits[-2]!='1900-01-01':
                if tunits[-2] =='1800-01-01':
                    offset=(datetime.strptime(fdict['tunits'][-2],'%Y-%m-%d')-datetime(1900,1,1)).total_seconds() #+ timedelta(days=1)
                else:
                    offset=(datetime.strptime(fdict['tunits'][-2],'%Y-%m-%d')-datetime(1900,1,1)).total_seconds()
            fak = 3600
            
        if offset != 0:
            
            fdict['offset'] = offset
                    
        a = 0.
        for sc in  'add_offset', 'scale_factor' :
            try:
                
                fdict['attribs'][sc] =np.float32(g[p].attrs[sc])
            except:
                fdict['attribs'][sc] = np.float32(a)
            a += 1

        fdict['secs'].append( np.int64(g['time'][:])*fak+int(offset))
        fdict['fi'].append(np.full_like(fdict['secs'][-1], l, dtype=np.int32))
        fdict['sf'].append((fdict['attribs']['scale_factor'], ))
        fdict['ao'].append((fdict['attribs']['add_offset'], ))
        fdict['mti'].append(np.arange(fdict['secs'][-1].shape[0], dtype=np.int32))
        l += 1
    
    for s in atup:   
        fdict[s] = np.concatenate(fdict[s])
    mask =(fdict['secs'] >= yminus) &   (fdict['secs'] <= yplus )
    for s in atup:
        #if s in ('fi', 'mti'):
            #fdict[s] = np.roll(fdict[s], 3)
        if s not in ('sf','ao'):            
            fdict[s] = fdict[s][mask]
    
    if 'offset' not in fdict.keys():
        fdict['offset'] =0        
    #fdict['sh'] = g[p].shape
    
    
    #fdict['ystenmin'] = np.min(lats) - np.abs(fps[1].latitude[1] - fps[1].latitude[0])
    #fdict['ystenmax'] = np.max(lats) + np.abs(fps[1].latitude[1] - fps[1].latitude[0])
    #fdict['xstenmin'] = np.min(lons) - np.abs(fps[1].longitude[1] - fps[1].longitude[0])
    #fdict['xstenmax'] = np.max(lons) + np.abs(fps[1].longitude[1] - fps[1].longitude[0])
    
    
    return fdict

def make_tc_sten(lats, lons,tss,lodt, fdict):
    
    sh = fdict['sh']

    dx=360/sh[1]
    dy=(fdict['latitude'][0] - fdict['latitude'][-1])/(sh[0]-sh[0] % 2)
        
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
    
    try:
        
        if ixrefmax >= ixrefmin:
            
            xsten = np.arange(ixrefmin-1, ixrefmax+3)#[ix - 1, ix, ix + 1, ix + 2])
        else:
            ixrefmax = int(np.floor((lonmax+360.)/dx))
            xsten = np.arange(ixrefmin-1, ixrefmax+3) % fdict['longitude'].shape[0]
        ysten = np.arange(iyrefmin-2, iyrefmax+3)#[iy - 1, iy, iy + 1, iy + 2])
        ysten[ysten>(sh[0]-1)] = sh[0] - 1
        ysten[ysten<0] = 0
        xsten[xsten>(sh[1]-1)] -= sh[1]
        xsten[xsten<0] += sh[1]
    except Exception as e:
        with open('errors.log', 'a') as f:
            f.write(fns+' size \n')
        if '0-20000-0-12330' in fns or '0-20000-0-26406' in fns or '0-20000-0-34247' in fns:
            return
        else:
            return
            #raise ValueError(fns)
    
    if len(ysten) *len(xsten) > 200000 and ('0-999' in fns or '20999-0-SHIP' in fns):
        print(fns, ': Too much latitude movement')
        return

    if xsten[0] == xsten[-1] - xsten.shape[0] + 1:
        #tt = time.time()
        #tera5=g[p][ysten[0]:ysten[-1] + 1,xsten[0]:xsten[-1] + 1,:,:]
        tlon = fdict['longitude'][:][xsten[0]:xsten[-1] + 1]
        ##print('5', time.time()-tt)
    else:
        ##tera5=g[p][ysten[0]:ysten[-1] + 1, :, :, :][:,xsten,:,:]
        #piv = np.where(xsten==0)[0][0]
        #tt = time.time()
        #tera5=np.concatenate((g[p][ysten[0]:ysten[-1] + 1, xsten[0]:xsten[piv-1] + 1, :, :],
                                #g[p][ysten[0]:ysten[-1] + 1, xsten[piv]:xsten[-1] + 1, :, :]),
                               #axis=1)
        ##print('51', time.time()-tt)
        tlon = fdict['longitude'][:][xsten]
    tlat = fdict['latitude'][ysten[0]:ysten[-1] + 1]
    
    fdict['tlon'] = tlat
    fdict['tlat'] = tlon
    fdict['xsten'] = xsten
    fdict['ysten'] = ysten
    
    return #tlon, tlat, xsten, ysten

def getchunk(fhash,fps, fis, p, ysten,xsten,i, idp,mti, sf,ao, tss, zs):
    
        
    if fps[fis[i]] is None:
        return np.full((4, 4, 2, 1), np.nan, dtype=np.float32)
    
    key = []
    for f in fis[i], fis[i] + 1:
        
        key.append((f, ysten.start, xsten[0]))
        if key[-1] not in fhash.keys() and f < len(fps) and fps[f] is not None: #
            
            if xsten[0] == xsten[-1] - xsten.shape[0] + 1:
                fhash[key[-1]]=fps[f][p][ysten, xsten[0]:xsten[0]+4,:,:] * sf[f] + ao[f]
                while fhash[key[-1]].shape[0] != 4:
                    fhash[key[-1]] = np.concatenate((fhash[key[-1]], fhash[key[-1]][-1:, :, :, :]), axis=0)
            else:
                piv = np.where(xsten==0)[0][0]
                fhash[key[-1]] = np.concatenate((fps[f][p][ysten, xsten[:piv] , :, :],
                                            fps[f][p][ysten, xsten[piv:], :, :]), axis=1)* sf[f] + ao[f]
                while fhash[key[-1]].shape[0] != 4:
                    fhash[key[-1]] = np.concatenate((fhash[key[-1]], fhash[key[-1]][-1:, :, :, :]), axis=0)
                    
        if mti < fhash[key[0]].shape[2] - 1:
            break
        #else:
            #if f < len(fps) and np.sum(fhash[key[-1]][:, :, mti:mti + 2, idp:idp + 1]- (fps[f][p][ysten, xsten[0]:(xsten[-1] + 1),mti:mti + 2, idp:idp + 1] * sf[f] + ao[f])) != 0:
                
                #print(fhash[key[-1]][:, :, mti:mti + 2, idp:idp + 1]- (fps[f][p][ysten, xsten[0]:(xsten[-1] + 1),mti:mti + 2, idp:idp + 1] * sf[f] + ao[f]))
                #x = 0
    if mti >= fhash[key[0]].shape[2] - 1:
        if fis[i] + 1 < len(fps) and fps[fis[i]+1] is not None:
            return np.concatenate((fhash[key[0]][:, :,mti:mti+ 1,idp:idp + 1], 
                                   fhash[key[1]][:, :,0:1,idp:idp + 1]), axis=2)
        else:
            return np.concatenate((fhash[key[0]][:, :,mti:mti+ 1,idp:idp + 1], 
                                   fhash[key[0]][:, :,mti:mti+ 1,idp:idp + 1]), axis=2)
    else:
        #if idp == 5:
            
            #print(key, fps[fis[i]], mti[i], fps[key[0][0]]['level'][idp], tss-np.int64(fps[key[0][0]]['time'][mti[i]])*3600, zs)
        return fhash[key[0]][:, :, mti:mti + 2, idp:idp + 1]

@njit(cache=True, boundscheck=True)
def getchunk3(fps, fis, ysten,xsten,i, idp,mti, sf,ao, tss, zs):
    
        
    #key = []
    #for f in fis[i], fis[i] + 1:
        
        #key.append((f, ysten.start, xsten[0]))
        #if key[-1] not in fhash.keys() and f < len(fps) and fps[f] is not None: #
            
            #if xsten[0] == xsten[-1] - xsten.shape[0] + 1:
                #fhash[key[-1]]=fps[f][p][ysten, xsten[0]:xsten[0]+4,:,:] * sf[f] + ao[f]
                #while fhash[key[-1]].shape[0] != 4:
                    #fhash[key[-1]] = np.concatenate((fhash[key[-1]], fhash[key[-1]][-1:, :, :, :]), axis=0)
            #else:
                #piv = np.where(xsten==0)[0][0]
                #fhash[key[-1]] = np.concatenate((fps[f][p][ysten, xsten[:piv] , :, :],
                                            #fps[f][p][ysten, xsten[piv:], :, :]), axis=1)* sf[f] + ao[f]
                #while fhash[key[-1]].shape[0] != 4:
                    #fhash[key[-1]] = np.concatenate((fhash[key[-1]], fhash[key[-1]][-1:, :, :, :]), axis=0)
                    
        #if mti < fhash[key[0]].shape[2] - 1:
            #break
        ##else:
            ##if f < len(fps) and np.sum(fhash[key[-1]][:, :, mti:mti + 2, idp:idp + 1]- (fps[f][p][ysten, xsten[0]:(xsten[-1] + 1),mti:mti + 2, idp:idp + 1] * sf[f] + ao[f])) != 0:
                
                ##print(fhash[key[-1]][:, :, mti:mti + 2, idp:idp + 1]- (fps[f][p][ysten, xsten[0]:(xsten[-1] + 1),mti:mti + 2, idp:idp + 1] * sf[f] + ao[f]))
                ##x = 0
    f =fis[i]        
    if mti >= fps[f].shape[2] - 1:
        if f + 1 < len(fps) and fps[f+1] is not None:
            return np.concatenate((fps[f][ysten[0]:ysten[0] + 4, xsten[0]:xsten[0] + 4,mti:mti+ 1,idp:idp + 1]* sf[f] + ao[f], 
                                   fps[f+1][ysten[0]:ysten[0] + 4, xsten[0]:xsten[0] + 4,0:1,idp:idp + 1]* sf[f+1] + ao[f+1]), axis=2)
        else:
            return np.concatenate((fps[f][ysten[0]:ysten[0] + 4, xsten[0]:xsten[0] + 4,mti:mti+ 1,idp:idp + 1], 
                                   fps[f][ysten[0]:ysten[0] + 4, xsten[0]:xsten[0] + 4,mti:mti+ 1,idp:idp + 1]), axis=2)* sf[f] + ao[f]
    else:
        #if idp == 5:
            
            #print(key, fps[f], mti[i], fps[key[0][0]]['level'][idp], tss-np.int64(fps[key[0][0]]['time'][mti[i]])*3600, zs)
        return fps[f][ysten[0]:ysten[0] + 4, xsten[0]:xsten[0] + 4, mti:mti + 2, idp:idp + 1]* sf[f] + ao[f]


def getchunk2(fhash,fps, fis, p, ysten,xsten,mti, idp):
    
        
    key = []
    for f in fis, fis+ 1:
        
        key.append((f, ysten.start, xsten[0]))
        if key[-1] not in fhash.keys() and f < len(fps) and fps[f] is not None: #
            
            if xsten[0] == xsten[-1] - xsten.shape[0] + 1:
                fhash[key[-1]]=fps[f][p][ysten, xsten[0]:xsten[0]+4,:,:] #* sf[f] + ao[f]
                while fhash[key[-1]].shape[0] != 4:
                    fhash[key[-1]] = np.concatenate((fhash[key[-1]], fhash[key[-1]][-1:, :, :, :]), axis=0)
            else:
                piv = np.where(xsten==0)[0][0]
                fhash[key[-1]] = np.concatenate((fps[f][p][ysten, xsten[0]:xsten[piv] , :, :],
                                            fps[f][p][ysten, xsten[piv]:, :, :]), axis=1)#* sf[f] + ao[f]
                while fhash[key[-1]].shape[0] != 4:
                    fhash[key[-1]] = np.concatenate((fhash[key[-1]], fhash[key[-1]][-1:, :, :, :]), axis=0)
        else:
            if np.sum(fhash[key[-1]][:, :, mti:mti + 2, idp:idp + 1]- (fps[f][p][ysten, xsten[0]:(xsten[-1] + 1),mti:mti + 2, idp:idp + 1])):# * sf[f] + ao[f])) != 0:
                
                print(fhash[key[-1]][:, :, mti:mti + 2, idp:idp + 1]- (fps[f][p][ysten, xsten[0]:(xsten[-1] + 1),mti:mti + 2, idp:idp + 1]))# * sf[f] + ao[f]))
                x = 0
            
    if mti >= fhash[key[0]].shape[2] - 1:
        if fis + 1 < len(fps) and fps[fis+1] is not None:
            return np.concatenate((fhash[key[0]][:, :,mti:mti+ 1,idp:idp + 1], 
                                   fhash[key[1]][:, :,0:1,idp:idp + 1]), axis=2)
        else:
            return np.concatenate((fhash[key[0]][:, :,mti:mti+ 1,idp:idp + 1], 
                                   fhash[key[0]][:, :,mti:mti+ 1,idp:idp + 1]), axis=2)
    else:
        #if idp == 5:
            
            #print(key, fps[fis[i]], mti[i], fps[key[0][0]]['level'][idp], tss-np.int64(fps[key[0][0]]['time'][mti[i]])*3600, zs)
        return fhash[key[0]][:, :, mti:mti + 2, idp:idp + 1]

def lin3d_new(fps,fn, p, fdict,tss,z,lons, lats, yms, out=None, out_fb=None ):
    #from scipy.interpolate import RectBivariateSpline
    tt=time.time()
    
    zs = z / 100.
    
    pres = fdict['level']
    #tera5 = np.array((z.shape[0], 2, 4, 4), dtype=np.float32) #dx * dy * 2tlevs * #tlevs
    
    its = np.searchsorted(fdict['secs'], tss) - 1
    its[its==-1] = 0
    fis = fdict['fi'][its]
    mti = fdict['mti'][its]
    iyref = np.searchsorted(-fdict['latitude'], -lats) - 2
    hlon = fdict['longitude']
    ixref = np.searchsorted(hlon, lons)
    ixref[ixref==hlon.shape[0]] = -1
    ixref -= 1
    
    ixrefmin = np.min(ixref)
    iyrefmin = np.min(iyref)
    iyrefmax = np.max(iyref+4)
    ixrefmax = np.max(ixref+4)
    
            
    mdict = {'xsten': [],'ysten': [],'secsten': [], 'crude': [],'crude2': [],}
    fhash = {}
    
    idps = np.searchsorted(pres, zs)
    
    ref =datetime(1900, 1, 1)
    
    reatab = np.zeros_like(lons, dtype=np.float32)

    xsten = np.empty((len(tss), 4), dtype=np.int32)
    ysten = np.empty((len(tss), 4), dtype=np.int32)
    #print(fps[1])
    for i in range(len(tss)):                        
            xsten[i, :] = np.arange(ixref[i], ixref[i] + 4, dtype=np.int32) % fdict['longitude'].shape[0]
            ysten[i, :] = np.clip(np.arange(iyref[i], iyref[i] + 4, dtype=np.int32),0,fdict['latitude'].shape[0]-1)
            
            #mdict['xsten'].append(xsten)
            #mdict['ysten'].append(ysten)
            
    fpm = []
    if False and ixrefmin > 0 and ixrefmax < fdict['longitude'].shape[0] - 4 and (ixrefmax - ixrefmin) * (iyrefmax - iyrefmin) < 400: # try if reading array into memory upfront to minimize I/O requests. Is not faster, thus deactivated.
        for i in range(len(fps)):
            if fps[i] is not None:
                fpm.append(fps[i][p][iyrefmin: iyrefmax, ixrefmin: ixrefmax,: , : ][:])
            else:
                fpm.append(None)
            #ixref -= ixrefmin
            #iyref -= iyrefmin
        fpm = tuple(fpm)
    

    cubes = np.empty((4, 4, 2*len(tss), 1), dtype=np.float32)
    
    
    for i in range(len(tss)):    
        try:
    
            #mdict['secsten'].append(its[i])
            #mdict['ps'].append(ps)
            #if len(mdict['idt']) > 0 and idt < mdict['idt'][-1]:
                #x = 0
            #mdict['idt'].append(idt)
            if False:
                offset = fdict['offset']
                deltat = tss[i] - fps[fis[i]]['time'][mti[i]]*3600 - offset
                if deltat < 0 or deltat > fdict['secs'][its[i] + 1] - fdict['secs'][its[i]]:
                    print(deltat,ref+timedelta(seconds=int(tss[i])), fps[fis[i]],
                          ref+timedelta(seconds=int(fps[fis[i]]['time'][mti[i]]*3600-offset)),
                          ref+timedelta(seconds=int(fdict['secs'][its[i] + 1] )))
                    x = 0
            if i == 4068:
                x = 0
            if len(fpm) == 0:
                
                cubes[:, :, 2*i:2*i+2, 0:1] =  getchunk(fhash,fps, fis, p, slice(ysten[i, 0], ysten[i, 0] + 4),xsten[i, :],i, idps[i],mti[i], fdict['sf'], fdict['ao'], tss[i], zs[i])
            else:
                cubes[:, :, 2*i:2*i+2, 0:1] =  getchunk3(fpm, fis, ysten[i, :]-iyrefmin,
                                   xsten[i, :]-ixrefmin,i, idps[i],mti[i], fdict['sf'], fdict['ao'], tss[i], zs[i])
                
            
            #try:
                
                #mask = tera5==fdict['attribs']['missing_value']
            #except:
                #mask = np.zeros_like(tera5, dtype=bool)
            #cubes.append(tera5) #*fdict['attribs']['scale_factor']+fdict['attribs']['add_offset'])
            
            if False:
                mdict['crude'].append(cubes[-1][1, 1, 0, 0] - out['observation_value'][i])
                fpsh = fps[fis[i]][p].shape[2]
                if mti[i] + 3 < fpsh:
                    
                    mdict['crude2'].append((fps[fis[i]][p][slice(ysten[0], ysten[0]+4), xsten[0]: xsten[0] + 4, mti[i] + 3:mti[i] + 4, 14:15] * fdict['sf'][fis[i]] + fdict['ao'][fis[i]])[1, 1, 0, 0] - out['observation_value'][i])
                else:
                    mdict['crude2'].append((fps[fis[i]+1][p][slice(ysten[0], ysten[0]+4), xsten[0] : xsten[0] + 4, (mti[i] + 3) % fpsh:(mti[i] + 3) % fpsh + 1, 14:15] * fdict['sf'][fis[i] + 1] + fdict['ao'][fis[i] + 1])[1, 1, 0, 0] - out['observation_value'][i])
                if out is not None and p == 't' and zs[i] == 925. and np.abs(cubes[-1][1, 1, 0, 0] - (out['observation_value'][i]-out_fb['fg_depar@body'][i])) > 10:
                    print(fps[fis[i]], zs[i], pres[idps[i]], cubes[-1][1, 1, 0, 0],
                          out['observation_value'][i], out['observation_value'][i]-out_fb['fg_depar@body'][i], out['observation_value'][i] - cubes[-1][1, 1, 0, 0], out_fb['fg_depar@body'][i])
                    x = 0
                        #cubes[-1][mask] = np.nan
                        
                        #if tss[k]-fdict['secs'][idt] < 0:
                            
                        #print(k, ps, tss[k]-fdict['secs'][idt],idt )
                    
    
        except MemoryError as e:
            print(fpattern.format(an1[0],an1[1]), ' extrapolating last hours of month')
            tera5 = np.concatenate((tera5, tera5[:, :, -1:, :], tera5[:, :, -1:, :]), axis=2)
            secs = np.concatenate((secs, secs[[-1]]+secs[-1]-secs[-2], secs[[-1]]+2*(secs[-1]-secs[-2])))
            pass
        #fi += 1
        #ptime(os.path.basename(fn), tt)
    
    #if len(idks) == 0:
        #reatab[:] = np.nan
        #return reatab, None, tss,zs,fdict['attribs']
        
    #udks = np.concatenate(idks)   
    #xsten = np.vstack(mdict['xsten'])
    #ysten = np.vstack(mdict['ysten'])
    #secsten = np.array(mdict['secsten'])
    #ps = np.vstack(mdict['ps'])[:, 0]
    #tera5 = np.concatenate(cubes, axis=2)
    
    #reatabraw = np.zeros_like(lons, dtype=np.float32)
    #if '0.25t' in fn:
        #x = 0
    #ptime(os.path.basename(fn)+ f' {len(fhash)} stencils, before trisimple', tt, end='')
    #ttt = time.time()
    trisimple(cubes, lons, lats, tss, xsten, ysten,its, 
              fdict['longitude'],
              fdict['latitude'],
              fdict['secs'], reatab)
    
    if False:
        idyy = np.where(zs==925.)[0]
        for m in range(12):        
            idm = np.where((tss[idyy]>=tss[idyy][0]+m*30.5*86400 )&(tss[idyy]<tss[idyy][0]+(m+1)*30.5*86400))[0]
            print(rms(reatab[idyy[idm]]-out['observation_value'][idyy[idm]]), rms(np.array(mdict['crude'])[idyy[idm]]), rms(np.array(mdict['crude2'])[idyy[idm]]), rms(out_fb['fg_depar@body'][idyy[idm]]), out_fb['fg_depar@body'][idyy[idm]].shape)            
        print(rms(reatab[idyy]-out['observation_value'][idyy]), rms(np.array(mdict['crude'])[idyy]), rms(np.array(mdict['crude2'])[idyy]), rms(out_fb['fg_depar@body'][idyy]), rms(reatab[idyy]-out['observation_value'][idyy]-out_fb['fg_depar@body'][idyy]), out_fb['fg_depar@body'][idyy].shape)            

    #if udks.shape[0] == reatab.shape[0]:
        #reatab[udks] = reatabraw
    #else:
        #print('WARNING: data gap in gridded files', os.path.basename(fns), os.path.basename(fn),reatab.shape[0],udks.shape[0] )
        #reatab[:] = np.nan
        #reatab[udks] = reatabraw[:udks.shape[0]]
        ##raise ValueError(fns)
        
    
    #ptime(os.path.basename(fn).split('_C')[0]+ f' {len(fhash)} stencils, before trisimple {ttt - tt:.4f} after trisimple', tt)


    if False and out is not None:
        idz = np.where(zs==100.)[0]
        plt.subplot(1, 2, 2)    
        plt.plot(tss[idz]/365.25/86400, out['observation_value'][idz]-reatab[idz], label=f'offline {rms(out["observation_value"][idz]-reatab[idz])[0]:.4f}')
        plt.plot(tss[idz]/365.25/86400, out_fb['fg_depar@body'][idz], label=f'online {rms(out_fb["fg_depar@body"][idz])[0]:.4f}')
        plt.legend()
        plt.subplot(1, 2, 1)
        plt.plot(tss[idz]/365.25/86400, out['observation_value'][idz], label='obs')
        plt.plot(tss[idz]/365.25/86400, out['observation_value'][idz]-out_fb['fg_depar@body'][idz], label='fg online')
        #plt.plot(tss[idz]/365.25/86400, tera5[1, 1, ::2, 0][idz])
        plt.plot(tss[idz]/365.25/86400, reatab[idz], label='fg offline')
        plt.legend()
        plt.show()
    if out is not None:
        pl = 100.
        idz = np.where(zs==pl)[0]
        print(f'{os.path.basename(fn)} {p} {int(pl)} hPa offline {rms(out[idz]-reatab[idz])[0]:.4f} online {rms(out_fb[idz])[0]:.4f}')

    return reatab, None, tss,zs,fdict['attribs']

def lin3d(fps,fns, p, fdict,tss,z,lons, lats, yms, out=None, out_fb=None ):
    #from scipy.interpolate import RectBivariateSpline
    tt=time.time()
    
    zs = z / 100.
    
    pres = fdict['level']
    pidx = fdict['pidx']
    tera5 = np.array((z.shape[0], 2, 4, 4), dtype=np.float32) #dx * dy * 2tlevs * #tlevs
    
    its = np.searchsorted(fdict['secs'], tss) - 1
    its[its==-1] = 0
    iyref = np.searchsorted(-fdict['latitude'], -lats) - 2
    hlon = fdict['longitude']
    ixref = np.searchsorted(hlon, lons)
    ixref[ixref==hlon.shape[0]] = -1
    ixref -= 1
    
    cubes = []
    mdict = {'xsten': [],'ysten': [],'idt': [], 'ps': [],}
    k = 0
    fi =0
    k0 = 0
    idks = []
    allcount = 0
    fhash = {}
    for fi in range(len(fps)):    
        try:
    
            #fn = fpattern.format(ans[0],ans[1])
            #with h5py.File(fn,'r') as g:
                
                if fps[fi] is None:
                    print('skipping file', fi)
                    if fi == len(fps) -1:
                        continue
                    elif fps[fi+1] is None:
                        continue
                else:
                    fn = os.path.basename(fps[fi].filename)
                    
                
                a = 0.
                for sc in  'add_offset', 'scale_factor' :
                    try:
                        
                        fdict['attribs'][sc] =fps[fi][p].attrs[sc]
                    except:
                        fdict['attribs'][sc] = a
                    a += 1
                    
                #print(fps[fi], fdict['attribs']['scale_factor'],fdict['attribs']['add_offset'] )
                    
               
            #sh = fps[fi][p].shape
                if fps[fi] is None:
                    pass
                    fsecs=np.array((fsecs[-1],)) #fsecs[-1]+(fsecs[-1]-fsecs[-2])))
                else:
                    fsecs=np.int64(fdict['offset']) + np.int64(fps[fi]['time'][:])*3600
                    if np.any((fsecs[1:]-fsecs[:-1])!=(fsecs[1]-fsecs[0])):
                        print('WARNING: inconsistent timestamps', fn)
                        continue
                    glevel = fps[fi]['level'][:]
                
                fsi = np.max((np.searchsorted(fsecs, np.min(tss)) - 1, 0)) # +1-(fsecs[1]-fsecs[0]))
                
                
                if fsi == fsecs.shape[0]:
                    if (np.min(tss)+1 - fsecs[-1]) < (fsecs[1]-fsecs[0]):
                        fsi = fsi - 1                        
                    else:
                        print(fps[fi], f'timedelta {np.min(tss)+1 - fsecs[-1]} too large, try next month')
                        fi += 1
                        continue
                    
                        #continue
                idt0 = np.searchsorted(fdict['secs'], fsecs[fsi])  # index where to put monthly data in large array
                
                if len(tss) > 1:
                    
                    maxleap = int(np.max((tss[1:]-tss[:-1])/(fsecs[1]-fsecs[0])+1))
                else:
                    maxleap = 1
                
                if '_pres' in fn:
                    x = 0

                l = 0
                for i in range(fsi, fsecs.shape[0]):
                                   
                    idt = idt0 + l #np.searchsorted(fdict['secs'], fsecs[i])  # index where to put monthly data in large array
                    l += 1
                    
                    if True:
                        idk = np.where(its[k0:k+maxleap+32]==idt)[0]
                        idkm1 = k0
                        if len(idk) == 0:
                            continue
                        
                        idk += k0
    
                        k0 = np.max((idkm1, idk[0] - 32))
                    else:
                        idk = np.where(its==idt)[0]
                    #if len(idk) == len(idk1):
                        #if np.any((idk-idk1)!=0):
                            #x = 0
                    #else:
                        #x = 0
                    
                    
                    #print(i, zs[idk], k0)
                    if 5881 in idk:
                        x =0
                    idks.append(idk)
                    j = 0
                    for k in idk:
                        
                        ps = zs[k]
                        idp = np.searchsorted(glevel, ps)
                        
                    
                        xsten = np.arange(ixref[k], ixref[k] + 4, dtype=np.int32) % fdict['longitude'].shape[0]
                        ysten = np.clip(np.arange(iyref[k], iyref[k] + 4, dtype=np.int32),0,fdict['latitude'].shape[0]-1)
                        
                        mdict['xsten'].append(xsten)
                        mdict['ysten'].append(ysten)
                        mdict['ps'].append(ps)
                        if len(mdict['idt']) > 0 and idt < mdict['idt'][-1]:
                            x = 0
                        mdict['idt'].append(idt)
                        
                                    
                        tera5 =  getchunk2(fhash,fps, fi, p, slice(ysten[0], ysten[0] + 4),xsten,i, idp)  
                        
                        try:
                            
                            mask = tera5==fdict['attribs']['missing_value']
                        except:
                            mask = np.zeros_like(tera5, dtype=bool)
                        cubes.append(tera5*fdict['attribs']['scale_factor']+fdict['attribs']['add_offset'])
                        
                        if out is not None and p == 't' and np.abs(cubes[-1][1, 1, 0, 0] - (out['observation_value'][k]-out_fb['fg_depar@body'][k])) > 30:
                            print(zs[k], glevel[idp], cubes[-1][1, 1, 0, 0],
                                  out['observation_value'][k], out['observation_value'][k]-out_fb['fg_depar@body'][k])
                            x = 0
                        cubes[-1][mask] = np.nan
                        
                        #if tss[k]-fdict['secs'][idt] < 0:
                            
                        #print(k, ps, tss[k]-fdict['secs'][idt],idt )
                    
    
        except MemoryError as e:
            print(fpattern.format(an1[0],an1[1]), ' extrapolating last hours of month')
            tera5 = np.concatenate((tera5, tera5[:, :, -1:, :], tera5[:, :, -1:, :]), axis=2)
            secs = np.concatenate((secs, secs[[-1]]+secs[-1]-secs[-2], secs[[-1]]+2*(secs[-1]-secs[-2])))
            pass
        fi += 1
        #ptime(os.path.basename(fn), tt)
    
    reatab = np.zeros_like(lons, dtype=np.float32)
    if len(idks) == 0:
        reatab[:] = np.nan
        return reatab, None, tss,zs,fdict['attribs']
        
    udks = np.concatenate(idks)   
    xsten = np.vstack(mdict['xsten'])
    ysten = np.vstack(mdict['ysten'])
    secsten = np.vstack(mdict['idt'])[:, 0]
    ps = np.vstack(mdict['ps'])[:, 0]
    tera5 = np.concatenate(cubes, axis=2)
    
                
    reatabraw = np.zeros_like(lons, dtype=np.float32)
    if '0.25t' in fn:
        x = 0
    trisimple(tera5, lons[udks], lats[udks], tss[udks], xsten, ysten,secsten, 
              fdict['longitude'],
              fdict['latitude'],
              fdict['secs'], reatabraw)
    
    if udks.shape[0] == reatab.shape[0]:
        reatab[udks] = reatabraw
    else:
        print('WARNING: data gap in gridded files', os.path.basename(fns), os.path.basename(fn),reatab.shape[0],udks.shape[0] )
        reatab[:] = np.nan
        reatab[udks] = reatabraw[:udks.shape[0]]
        #raise ValueError(fns)
        
    
    
    ptime(os.path.basename(fn)+ ' after trisimple', tt)

    if False and out is not None:
        idz = np.where(zs==100.)[0]
        plt.subplot(1, 2, 2)    
        plt.plot(tss[idz]/365.25/86400, out['observation_value'][idz]-reatab[idz], label=f'offline {rms(out["observation_value"][idz]-reatab[idz]):.4f}')
        plt.plot(tss[idz]/365.25/86400, out_fb['fg_depar@body'][idz], label=f'online {rms(out_fb["fg_depar@body"][idz]):.4f}')
        plt.legend()
        plt.subplot(1, 2, 1)
        plt.plot(tss[idz]/365.25/86400, out['observation_value'][idz])
        plt.plot(tss[idz]/365.25/86400, out['observation_value'][idz]-out_fb['fg_depar@body'][idz])
        #plt.plot(tss[idz]/365.25/86400, tera5[1, 1, ::2, 0][idz])
        plt.plot(tss[idz]/365.25/86400, reatab[idz])
        plt.show()
    if out is not None:
        idz = np.where(zs==100.)[0]
        print(f'{os.path.basename(fn)} 100 hPa offline {rms(out["observation_value"][idz]-reatab[idz])[0]:.4f} online {rms(out_fb["fg_depar@body"][idz])[0]:.4f}')

    return reatab, None, tss,zs,fdict['attribs']


def lin4d(fpattern,fns, out, p,pn, fdict,ts, ori, obstype, latorig,lonorig,z, refs,ans, drift=False, station_elevation=0., driftonly=False):
    #from scipy.interpolate import RectBivariateSpline
    tt=time.time()
    
    if ans[3] == ans[2]: # no data
        return
    if not driftonly:   # check if file exists before calculating drift, which can be costly. Skip if only drift is needed.
        fn = fpattern.format(ans[0],ans[1])
        if not os.path.isfile(fn):
            fn = ''.join(fn.split('.0.25'))       
            if not os.path.isfile(fn):
                return
    
    sl = slice(ans[2], ans[3])
    lons = lonorig[sl]
    lons[lons>360.] -= 360.
    lats = latorig[sl]
    tss = ts[sl]
    obss = obstype[sl]
    zs = z[sl] / 100.
    ovs = out['observation_value'][sl]
    
    fn = fpattern.format(ans[0],ans[1])
    if not os.path.isfile(fn):
        fn = ''.join(fn.split('.0.25'))
        
    try:

        #fn = fpattern.format(ans[0],ans[1])
        with h5py.File(fn,'r') as g:
            
            for sc in  'scale_factor', 'add_offset':
                
                fdict['attribs'][sc] =g[p].attrs[sc]    
                
            pres = fdict['level']
            pidx = fdict['pidx']
            sh = g[p].shape

            if False:
                if lats[0]<89.0 and lats[0] > -89.0:  # no drift if balloon started near pole
                    
                    lats = lats + out['latd'][sl]
                    lons = lons + out['lond'][sl]
                    
                lons[lons>=360.] -= 360
                lonmin = np.nanmin(lons)
                lonmax = np.nanmax(lons)
                atcycle = lonmax - lonmin > 300.
                if atcycle:
                    lons[lons>=180.] -= 360
                else:
                    lons[lons>=360.] -= 360
                    
                #lons[lons<0] += 360.
                tss = tss + out['timed'][sl]
                            
                
                lonmin = np.nanmin(lons)
                lonmax = np.nanmax(lons)
                latmin = np.nanmin(lats)
                latmax = np.nanmax(lats)
                dx=360/sh[1]
                dy=(fdict['latitude'][0] - fdict['latitude'][-1])/(sh[0]-sh[0] % 2)
                    
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
                
                try:
                    
                    if ixrefmax >= ixrefmin:
                        
                        xsten = np.arange(ixrefmin-1, ixrefmax+3)#[ix - 1, ix, ix + 1, ix + 2])
                    else:
                        ixrefmax = int(np.floor((lonmax+360.)/dx))
                        xsten = np.arange(ixrefmin-1, ixrefmax+3) % fdict['longitude'].shape[0]
                    ysten = np.arange(iyrefmin-2, iyrefmax+3)#[iy - 1, iy, iy + 1, iy + 2])
                    ysten[ysten>(sh[0]-1)] = sh[0] - 1
                    ysten[ysten<0] = 0
                    xsten[xsten>(sh[1]-1)] -= sh[1]
                    xsten[xsten<0] += sh[1]
                except Exception as e:
                    with open('errors.log', 'a') as f:
                        f.write(fns+' size \n')
                    if '0-20000-0-12330' in fns or '0-20000-0-26406' in fns or '0-20000-0-34247' in fns:
                        return
                    else:
                        return
                        #raise ValueError(fns)
                
                if len(ysten) *len(xsten) > 200000 and ('0-999' in fns or '20999-0-SHIP' in fns):
                    print(fns, ': Too much latitude movement')
                    return
    
                if xsten[0] == xsten[-1] - xsten.shape[0] + 1:
                    tt = time.time()
                    tera5=g[p][ysten[0]:ysten[-1] + 1,xsten[0]:xsten[-1] + 1,:,:]
                    tlon = fdict['longitude'][:][xsten[0]:xsten[-1] + 1]
                    #print('5', time.time()-tt)
                else:
                    #tera5=g[p][ysten[0]:ysten[-1] + 1, :, :, :][:,xsten,:,:]
                    piv = np.where(xsten==0)[0][0]
                    tt = time.time()
                    tera5=np.concatenate((g[p][ysten[0]:ysten[-1] + 1, xsten[0]:xsten[piv-1] + 1, :, :],
                                            g[p][ysten[0]:ysten[-1] + 1, xsten[piv]:xsten[-1] + 1, :, :]),
                                           axis=1)
                    #print('51', time.time()-tt)
                    tlon = fdict['longitude'][:][xsten]
                tlat = fdict['latitude'][ysten[0]:ysten[-1] + 1]
                
            else:
                if xsten[0] == xsten[-1] - xsten.shape[0] + 1:
                    tera5=g[p][ysten[0]:ysten[-1] + 1,xsten[0]:xsten[-1] + 1,:,:]
                    #print('5', time.time()-tt)
                else:
                    #tera5=g[p][ysten[0]:ysten[-1] + 1, :, :, :][:,xsten,:,:]
                    piv = np.where(xsten==0)[0][0]
                    tera5=np.concatenate((g[p][ysten[0]:ysten[-1] + 1, xsten[0]:xsten[piv-1] + 1, :, :],
                                            g[p][ysten[0]:ysten[-1] + 1, xsten[piv]:xsten[-1] + 1, :, :]),
                                           axis=1)
                
            try:
                if 'scale_factor' in  fdict['attribs']:
                    
                    if np.any(tera5==fdict['attribs']['missing_value']):
                        
                        mv=np.where(tera5==fdict['attribs']['missing_value'])
                        tera5 = tera5 * np.float32(fdict['attribs']['scale_factor'])+np.float32(fdict['attribs']['add_offset'])
                        tera5[mv]=np.nan
                    else:
                        tera5 = tera5 * np.float32(fdict['attribs']['scale_factor'])+np.float32(fdict['attribs']['add_offset'])
                else:
                    x = 0
                    
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
            if False:
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
            
        fn1 = fpattern.format(an1[0],an1[1])
        if not os.path.isfile(fn1):
            fn1 = ''.join(fn1.split('.0.25'))
        
        
        try:
            if '0.25' in fn1 and '0.25' not in fn or '0.25' not in fn1 and '0.25' in fn:
                raise FileNotFoundError
            with h5py.File(fn1,'r') as g:
                
                if xsten[0] == xsten[-1] - xsten.shape[0] + 1:
                    tt = time.time()
                    tera51=g[p][ysten[0]:ysten[-1] + 1,xsten[0]:xsten[-1] + 1,0:2,:]
                    #print('52', time.time()-tt)
                else:
                    #tera51=g[p][ysten[0]:ysten[-1] + 1, :, :, :][:,xsten,0:2,:]
                    
                    tt = time.time()
                    tera51=np.concatenate((g[p][ysten[0]:ysten[-1] + 1, xsten[0]:xsten[piv-1] + 1, 0:2, :],
                                            g[p][ysten[0]:ysten[-1] + 1, xsten[piv]:xsten[-1] + 1, 0:2, :]),
                                           axis=1)
                    #print('53', time.time()-tt)
                try:
                    if 'scale_factor' in  fdict['attribs'].keys():
                        
                        if np.any(tera51==fdict['attribs']['missing_value']):
                            mv=np.where(tera51==fdict['attribs']['missing_value'])
                            tera51 =tera51 * np.float32(g[p].attrs['scale_factor'])+ np.float32(g[p].attrs['add_offset'])
                            tera51[mv]=np.nan
                        else:
                            tera51 =tera51 * np.float32(g[p].attrs['scale_factor'])+ np.float32(g[p].attrs['add_offset'])
                except:
                    pass
    
                if '20CRv3' in fpattern:
                    for it in 0,1 :                       
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
            secs = np.concatenate((secs, secs[[-1]]+secs[-1]-secs[-2], secs[[-1]]+2*(secs[-1]-secs[-2])))

        except FileNotFoundError as e:
            print(fpattern.format(an1[0],an1[1]), ' extrapolating last hours of month')
            tera5 = np.concatenate((tera5, tera5[:, :, -1:, :], tera5[:, :, -1:, :]), axis=2)
            secs = np.concatenate((secs, secs[[-1]]+secs[-1]-secs[-2], secs[[-1]]+2*(secs[-1]-secs[-2])))
            pass

        #print(time.time() - tt)
        #if np.any(np.isnan(tera5)):
            #x = 0
        
        reatab = np.zeros_like(lons)
        reatab.fill(np.nan)
        #params = {'t': 126}
        if pidx.shape[0] == tera5.shape[3]:
            tpi = tera5
        else:
            tpi = tera5[:, :, :, pidx]
        pri = pres[pidx]
        #print(os.path.basename(fns).split('_')[0], os.path.basename(fn), tpi.shape, time.time()-tt)
        if ans[1] == 7:
            
            ptime(f'before {os.path.basename(fns).split('_')[0]} {os.path.basename(fn)}', tt)
        idx = triint(tpi, reatab, tlon, tlat, lons, lats, secs, tss,ori, obss,pri,zs, pn, dy, fns, 5)
        
        if ans[1] == 7:
            
            ptime(f'after {os.path.basename(fns).split('_')[0]} {os.path.basename(fn)}', tt)
            #ptime(f'after {os.path.basename(fns).split('_')[0]} {os.path.basename(fn)}, {np.abs(np.nanmax(reatab[:idx.shape[0]])):.4f} {np.nanmin(reatab[:idx.shape[0]]):.4f}', tt)
    except FileNotFoundError as e:
        print(fn, 'not available, continuing')
        return 
    return reatab[:len(idx)], ans[2] + idx, secs,pres,fdict['attribs']

# lin5d works with small lat/lon chunks in contrast to lin4d - suitable especially for moving platforms
def lin5d(fpattern,fns, out, p,pn, fdict,ts, obstype, latorig,lonorig,z, refs,ans, drift=False, station_elevation=0., driftonly=False):
    #from scipy.interpolate import RectBivariateSpline
    tt=time.time()
    
    if ans[3] == ans[2]: # no data
        return
    sl = slice(ans[2], ans[3])
    lons = lonorig[sl]
    lons[lons>360.] -= 360.
    lats = latorig[sl]
    tss = ts[sl]
    obss = obstype[sl]
    zs = z[sl] / 100.
    ovs = out['observation_value'][sl]
    
    extract_tuv(out, lats, lons, obss, ovs, zs, tss, sl, station_elevation)
    
    
        #latd,lond,timed = rsd.numba_drift.trajectory(lats[ii],lons[ii], np.array(df.u), np.array(df.v),
                                                 #zs[ii:iip], np.array(df.airTemperature))            
    idx = np.where((zs==500) & (obss==140))[0] 
    if np.any(out['timed'][sl][idx]<100):
        print(fns, 'X')
        print(out['timed'][sl][idx])
        if np.sum(out['timed'][sl][idx]<100) > 1 and len(np.unique(zs)) > 1:
            if '0-20999-0-0K1X' not in fns:
                raise ValueError(fns)
    idx = np.where(obss==140)[0]
    if len(idx) > 0 and np.any(np.max(np.abs(out['latd'][sl][idx]))>10.):
        print(fns, 'latd too big',np.max(np.abs(out['latd'][sl][idx]) ))
        #raise ValueError(fns)
        
    #if driftonly:
        #return
    fn = fpattern.format(ans[0],ans[1])
    if not os.path.isfile(fn):
        fn = ''.join(fn.split('.0.25'))
        
    try:

        #fn = fpattern.format(ans[0],ans[1])
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
            
            dx=360/sh[1]
            dy=(fdict['latitude'][0] - fdict['latitude'][-1])/(sh[0]-sh[0] % 2)

            lons = lons + out['lond'][sl]
            lons[lons>=360.] -= 360
            lonmin = np.nanmin(lons)
            lonmax = np.nanmax(lons)
            atcycle = lonmax - lonmin > 300.
            if atcycle:
                lons[lons>=180.] -= 360
            else:
                lons[lons>=360.] -= 360
                
            #lons[lons<0] += 360.
            lats = lats + out['latd'][sl]
            tss = tss + out['timed'][sl]
                        
            loni = np.searchsorted(g['longitude'][:], lons)
            lati = np.searchsorted(-g['latitude'][:], -lats)
            gtime = np.int64(g['time'][:]) * 3600
            ti = np.searchsorted(gtime, tss)
            if not np.any(ti<gtime.shape[0]):
                return

            tu, ri = np.unique(ts[sl], return_index=True)
            
            lri = len(ri)
            
            sh = g[p].shape
            ix = 0
            iy = 0
            for i in range(lri-1):
                if i == lri - 1:
                    rimax = len(loni)
                else:
                    rimax = ri[i+1]
                ix = np.max((ix, np.unique(loni[ri[i]:rimax]).shape[0]))
                iy = np.max((iy, np.unique(lati[ri[i]:rimax]).shape[0]))
                #print(i, np.max(loni[ri[i]:rimax])-np.min(loni[ri[i]:rimax]),
                      #np.max(lati[ri[i]:rimax])-np.min(lati[ri[i]:rimax]),
                      #)
                #print(ix, iy)
            
            tera5s = []
            tlons = []
            tlats = []
            tt = time.time()
            for i in range(lri):
                if i == lri - 1:
                    rimax = len(loni)
                    timax = len(loni)
                else:
                    rimax = ri[i+1]
                    timax = ti[ri[i+1]]
                xs0 = np.min(loni[ri[i]:rimax]) - 1
                if xs0 < 0:
                    xs0 += sh[1]
                xs = slice(xs0,xs0+ix+3)
                ys0 = np.min(lati[ri[i]:rimax]) - 1
                ys = slice(ys0,ys0+iy+3)
                it = np.min((sh[2]-1,ti[ri[i]] ))
                if xs.stop > sh[1]:
                    xs = np.arange(xs0, xs0+ix+3)
                    
                    xs[xs>=sh[1]] -= sh[1]
                    tera5s.append(np.concatenate((g[p][ys, xs[0]:, it:it + 1, :], g[p][ys, :xs[-1] + 1, it:it + 1, :]), axis=1))
                    tlons.append(np.concatenate((fdict['longitude'][xs[0]:],fdict['longitude'][:xs[-1] + 1])))
                else:
                    tera5s.append(g[p][ys, xs, it:it + 1, :])
                    tlons.append(fdict['longitude'][:][xs])
                tlats.append(fdict['latitude'][:][ys])

            tera5 = np.concatenate(tera5s, axis=2)
            tlon = np.concatenate(tlons)
            tlat = np.concatenate(tlats)
            print(time.time()-tt)
            
                    

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
            
        fn1 = fpattern.format(an1[0],an1[1])
        if not os.path.isfile(fn1):
            fn1 = ''.join(fn1.split('.0.25'))
        
        
        try:
            if '0.25' in fn1 and '0.25' not in fn or '0.25' not in fn1 and '0.25' in fn:
                raise FileNotFoundError
            with h5py.File(fn1,'r') as g:
                
                gtime = np.int64(g['time'][:]) * 3600
                ti = np.searchsorted(gtime, tss)
                tera5s = []
                tt = time.time()
                for i in range(lri):
                    if i == lri - 1:
                        rimax = len(loni)
                        timax = len(loni)
                    else:
                        rimax = ri[i+1]
                        timax = ti[ri[i+1]]
                    xs0 = np.min(loni[ri[i]:rimax]) - 1
                    if xs0 < 0:
                        xs0 += sh[1]
                    xs = slice(xs0,xs0+ix+3)
                    ys0 = np.min(lati[ri[i]:rimax]) - 1
                    ys = slice(ys0,ys0+iy+3)
                    if ti[ri[i]] > np.max(tss):
                        continue
                    it = ti[ri[i]]
                    if xs.stop > sh[1]:
                        xs = np.arange(xs0, xs0+ix+3)
                        
                        xs[xs>=sh[1]] -= sh[1]
                        tera5s.append(np.concatenate((g[p][ys, xs[0]:, it:it + 1, :], g[p][ys, :xs[-1] + 1, it:it + 1, :]), axis=1))
                    else:
                        tera5s.append(g[p][ys, xs, it:it + 1, :])
    
                tera51 = np.concatenate(tera5s, axis=2)
    
                if '20CRv3' in fpattern:
                    for it in 0,1 :                       
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
            secs = np.concatenate((secs, secs[[-1]]+secs[-1]-secs[-2], secs[[-1]]+2*(secs[-1]-secs[-2])))

            if np.any(tera5==fdict['attribs']['missing_value']):
                
                mv=np.where(tera5==fdict['attribs']['missing_value'])
                tera5 = tera5 * np.float32(fdict['attribs']['scale_factor'])+np.float32(fdict['attribs']['add_offset'])
                tera5[mv]=np.nan
            else:
                tera5 = tera5 * np.float32(fdict['attribs']['scale_factor'])+np.float32(fdict['attribs']['add_offset'])

        except FileNotFoundError as e:
            print(fpattern.format(an1[0],an1[1]), ' extrapolating last hours of month')
            tera5 = np.concatenate((tera5, tera5[:, :, -1:, :], tera5[:, :, -1:, :]), axis=2)
            secs = np.concatenate((secs, secs[[-1]]+secs[-1]-secs[-2], secs[[-1]]+2*(secs[-1]-secs[-2])))
            pass

        #print(time.time() - tt)
        #if np.any(np.isnan(tera5)):
            #x = 0
        
        reatab = np.zeros_like(lons)
        reatab.fill(np.nan)
        #params = {'t': 126}
        if pidx.shape[0] == tera5.shape[3]:
            tpi = tera5
        else:
            tpi = tera5[:, :, :, pidx]
        pri = pres[pidx]
        #print(os.path.basename(fns).split('_')[0], os.path.basename(fn), tpi.shape, time.time()-tt)
        #idx = triint(tpi, reatab, tlon, tlat, lons, lats, secs, tss,obss,pri,zs, pn, dy, fns)
        
        idx = np.array([0])
        if idx.shape[0] > 0 and ans[1] == 1:
            
            print(os.path.basename(fns).split('_')[0], os.path.basename(fn), np.abs(np.nanmax(reatab[:idx.shape[0]])), np.nanmin(reatab[:idx.shape[0]]), time.time()-tt)
    except FileNotFoundError as e:
        print(fn, 'not available, continuing')
        return 
    return reatab[:len(idx)], ans[2] + idx, secs,pres,fdict['attribs']



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
                if np.any(tera5==fdict['attribs']['missing_value']):
                    
                    mv=np.where(tera5==fdict['attribs']['missing_value'])
                    tera5 = tera5 * np.float32(fdict['attribs']['scale_factor'])+np.float32(fdict['attribs']['add_offset'])
                    tera5[mv]=np.nan
                else:
                    tera5 = tera5 * np.float32(fdict['attribs']['scale_factor'])+np.float32(fdict['attribs']['add_offset'])
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
                    if np.any(tera51==fdict['attribs']['missing_value']):
                        
                        mv=np.where(tera51==fdict['attribs']['missing_value'])
                        tera51 =tera51 * np.float32(g[p].attrs['scale_factor'])+ np.float32(g[p].attrs['add_offset'])
                        tera51[mv]=np.nan
                    else:
                        tera51 =tera51 * np.float32(g[p].attrs['scale_factor'])+ np.float32(g[p].attrs['add_offset'])
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
        if np.any(np.isnan(tera5)):
            x = 0
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
        
        if idx.shape[0] > 0:
            
            print(os.path.basename(fns), ans,p, np.abs(np.nanmax(reatab[:idx.shape[0]])), np.nanmin(reatab[:idx.shape[0]]), time.time()-tt)
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
       
    #from scipy.interpolate import RectBivariateSpline
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

def readictplot(readict, tasks, plevs, figprefix, marker='*', minvals=50, do_plot=True):
    
    res ={'fn' : figprefix}    
    for plev in plevs:
        iplev = np.int32(plev/100)
        res[iplev] = {}
        ida = np.where(readict['obs']['z_coordinate']==plev)[0]
        for p in readict['obstypes'].keys():
            res[iplev][p] = {}
        

        #plev = plevs[0]

            idx=ida[np.where(readict['obs']['obstype'][ida]==readict['obstypes'][p])[0]]
            if p == 'u':
                idxv = ida[np.where(readict['obs']['obstype'][ida]==readict['obstypes']['v'])[0]]
            if len(idx) < minvals:
                continue

            for sfunc in tasks:
                good = 0
                res[iplev][p][sfunc] = {}
                for k,v in readict.items():
                    if 'obs' in k or '20CR' in k:
                        continue
                    if k =='era5fb':
                        i = 0
                    
                    if 'ftype' in v.keys():                    
                        iterable = v['ftype']
                    else:
                        iterable = v
                        
                    res[iplev][p][sfunc][k] = {}
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
                            print(p, k, q, plev)
                            ystart = (ref + timedelta(seconds=int(readict['obs']['date_time'][idx[0]]))).year
                            ystop = (ref + timedelta(seconds=int(readict['obs']['date_time'][idx[-1]]))).year + 1
                            smarker = marker
                            if ystop - ystart == 1:
                                marker = '*'
                            
                            roo = readict['obs']['obs'][idx]
                            rref = readict[k][dtype]['refvalues'][idx]
                            
                            for iy in range(ystart, ystop):
                                #print(iy)
    
                                idy=np.where(tsy==iy)[0]
                                if p != 'u':
                                    
                                    rv = np.clip(roo[idy] - rref[idy], q[0], q[1])
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
                                if len(idy)>minvals and np.sum(~np.isnan(rv)) > minvals:
                                    if sfunc == 'rms':
    
                                        rms.append(np.sqrt(np.nanmean((rv)**2)))
                                    elif sfunc == 'std':
                                        rms.append(np.nanstd(rv))
                                    else:
                                        rms.append(np.nanmean(rv))
    
                                else:
                                    rms.append(np.nan)
                                years.append(iy)
                            #x = np.array(years),np.array(rms)
                            nrms = np.array(rms)
                            ngood = np.sum(~np.isnan(nrms))
                            res[iplev][p][sfunc][k][dtype] = (np.array(years), nrms)
                            if ngood > 0:
                                    
                                if ngood <= 2:
                                    marker = '*'
                                if do_plot:
                                    plt.plot(np.array(years),nrms,'-'+marker, 
                                         label='obs -'+k+'_'+dtype+', '+sfunc+'= {:5.3f}'.format(np.sqrt(np.nanmean(nrms)**2)))
                                good += 1
                if good == 0:                    
                    continue
                
                if p == 'u':
                    pt = 'winddir'
                    obsu = 'deg'
                else:
                    pt = p
                    obsu = readict['obsunits'][p]
                if do_plot:
                    plt.title('Monthly '+sfunc+' reanalysis '+pt+' departures, '+str(int(plev/100)) + ' hPa, '+figprefix.split('/')[-1].split('_')[0])
                    plt.ylabel(sfunc+' ['+obsu+']')
                    plt.legend()
                    plt.tight_layout()
                    #plt.xlim(1935, 1965)
                    plt.savefig(figprefix+'_'+p+'_'+str(int(plev/100))+'_'+str(minvals)+'_'+sfunc+'stats.png')
                    plt.close()
    return res

def retrieve_anfg(f, out, out_fb, readict, ts, tsu, tslice, refs, gactor, out_name,path_to_gridded,
                  drift=False, station_elevation=0., driftonly=False, ri=None, ori=None):

    
    tt=time.time()
    cyear = int(out_name.split('/')[-2])
    outdir = os.path.dirname(out_name) + '/'
    fn = f.filename

    if True:
        try:

            
            #otsu, ori=np.unique(out['date_time'], return_index=True)
            #_, ri=np.unique(ts, return_index=True)
            
            #otsu = ts[ri[:-1]]
            
            lon = np.empty(out['date_time'].shape[0], dtype=np.float32)
            lat = np.empty(out['date_time'].shape[0], dtype=np.float32)
            fill_restdata(lat, f['observations_table']['latitude'][tslice], ori, ri[:-1], '','')
            fill_restdata(lon, f['observations_table']['longitude'][tslice], ori, ri[:-1], '','')

                
            lon[lon>360.] -= 360.
            lon[lon<0] += 360.
                
            obstype=out['observed_variable']
            obs=out['observation_value']
            z=out['z_coordinate']
            #latdis = lat + (100000. - z) * 3. / 100000.
            #londis = lon + (100000. - z) * 3. / 100000.
            #zt=out['z_coordinate_type']
            ofb=True
        except MemoryError as e:
            print(fn, cyear, e)
            return None

    if drift:
        try:
            
            lodt = out['date_time'].copy()
            for i in range(ori.shape[0]-1):
                lodt[ori[i]:ori[i+1]] = tsu[i]
            i = ori.shape[0]-1
            lodt[ori[i]:] = tsu[-1]
        except MemoryError as e:
            raise ValueError(fn+str(e))
        
#        extract_tuv(out, lats, lons, obss, ovs, zs, tss, sl, station_elevation)  # calculates balloon drift
        sl = slice(0, lon.shape[0])
        extract_tuv(out, lat, lon, out['observed_variable'], out['observation_value'], z/100., lodt, sl, station_elevation)  # calculates balloon drift
    
        out['date_time'][:] = lodt
    
        #latd,lond,timed = rsd.numba_drift.trajectory(lats[ii],lons[ii], np.array(df.u), np.array(df.v),
                                                 #zs[ii:iip], np.array(df.airTemperature
        #some checks
        idx = np.where(obstype==140)[0]
        idy = np.where(z[idx]==50000)[0]
        if np.any(out['timed'][sl][idx[idy]]<100):
            print(os.path.basename(fn), 'X', np.sum(out['timed'][sl][idx[idy]]<100), 'ascent values<100')
            #print(out['timed'][sl][idx])

        #idx = np.where(obstype==140)[0]
        if len(idx) > 0 and np.any(np.max(np.abs(out['latd'][sl][idx]))>10.):
            
            idz = np.argmax(np.abs(out['latd'][sl][idx]) )
            print(fn, 'latd too big',out['latd'][sl][idx][idz],out['lond'][sl][idx][idz], out['timed'][sl][idx][idz])
            
            #raise ValueError(fn, 'latd too big',np.max(np.abs(out['latd'][sl][idx]) ))
            out['latd'][sl][idx] = np.clip(out['latd'][sl][idx], -1.0, 1.0)
            #raise ValueError(fn)
        
        ptime('after drift', tt)        
        
    if driftonly:
        return    
        
    ptime('before read gridded', tt)
    obstypes={'t':ipar[85],'u':ipar[104],'v':ipar[105],'q':ipar[39],'z':ipar[117]}    
    obsunits={'t':'K','u':'m/s','v':'m/s','q':'g/kg','z':'m^2/s^2'}
    plevs = np.array([10, 20, 30, 50, 70, 100, 150, 200, 250, 300, 400, 500, 700, 850, 925, 1000])
    idz = np.where(np.isin(z, plevs*100))[0]
    for k,v in readict.items():
        #if k!='era5fc':
            #continue
        for ftype in v['ftype']:  
            if 'fc' in ftype:
                dtype='fc'
            else:
                dtype='an'
            if dtype not in readict[k].keys():
                readict[k][dtype]={}
                readict[k][dtype]['refvalues']=np.empty(obs[idz].shape,dtype=np.float32)
                readict[k][dtype]['refvalues'].fill(np.nan)            

            vpk = list(v['param'].keys())
            for p,pn in v['param'].items():
                
                #if p == 'q':
                    #continue
                if k!='JRA55':

                    fpattern=v['path']+v['prefix'][ftype]+v['glue']+'{}{:0>2}'+v['glue']+pn+v['suffix']+'.nc'
                else:
                    fpattern=v['path']+v['prefix']+ftype+v['glue']+pn+'.reg_tl319.{}{:0>2}'+v['glue']+v['suffix']+'.nc'
                #found=len(glob.glob(fpattern.format(*yms[0])))
                #print(fpattern.format(*yms[0]),found)
                #if found:
                ref = datetime(1900, 1, 1)
                
                idw = np.where(obstype[idz]==obstypes[p])[0]
                idy = idz[idw]
                
                if len(idy) == 0:
                    continue
                
                tss = lodt[idy] + out['timed'][idy]
                
                tstart = ref + timedelta(seconds=int(np.min(tss))-43200)
                tstop = ref + timedelta(seconds=int(np.max(tss)-1)+43200)
                ylist = sorted(list(set((tstart.year, cyear, tstop.year))))
                fdict = {}
                yms = []
                fps = []
                goodfiles = 0
                for y in ylist:
                    mstart = 1
                    if y < cyear:
                        mstart = 12
                    for m in range(mstart, 13):
                        
                        yms.append((y, m))

                        fng = fpattern.format(*yms[-1])
                        if not os.path.isfile(fng):
                            fng = ''.join(fng.split('.0.25'))
                        try:
                            
                            fps.append(h5py.File(fng))
                            goodfiles += 1
                        except:
                            pass
                        #fps.append(None)
                        
                        if y > cyear:
                            break
                        
                if goodfiles <= 1: # do not process if only dec of past year or jan of next year is available
                    break
                
                   
                
                lats = lat[idy] # copy because drift will be added to lats, but not lat
                lons = lon[idy]
                if True:
                    if lats[0]<89.0 and lats[0] > -89.0:  # no drift if balloon started near pole
                        
                        lats = lats + out['latd'][idy]
                        lons = lons + out['lond'][idy]                
                        
                    lons[lons>=360.] -= 360
                    lons[lons<0.] += 360
                
                #if '0.25t' not in fpattern:
                    #continue
                fdict = make_fdict(fps, p, cyear, refs)
                debug = True
                #if debug:
                #tups = lin3d(fps,fn, p, fdict,tss,z[idy],lons,lats, yms, out[idy], out_fb[idy])
                #else:
                tups = lin3d_new(fps,fn, p, fdict,tss,z[idy],lons,lats, yms)#, out['observation_value'][idy], out_fb['fg_depar@body'][idy])
                                
                #idw = np.where(obstype[idz]==obstypes[p])[0]
                #idy = idz[idw]
                # readict contains only standard pressure levels. Thus needs idw as index, not idyi
                readict[k][dtype]['refvalues'][idw] = tups[0]
                if False:
                    i = 7
                    idp = np.where(z[idy]==plevs[i]*100.)[0]
                    while i < 15 and len(idp) < 2:
                        i += 1
                        idp = np.where(z[idy]==plevs[i]*100.)[0]
                        
                    r1,n1 = rms(out["observation_value"][idy[idp]]-readict[k][dtype]['refvalues'][idw[idp]])
                    if 'fc' in dtype:
                        r2,n2 = rms(out_fb["fg_depar@body"][idy[idp]])
                    else:
                        r2,n2 = rms(out_fb["an_depar@body"][idy[idp]])
                                        
                    try:
                        
                        om = int(np.log10(r1))
                    except:
                        om = 0
                    fak = 1
                    if om < -3:
                        fak *= 1000
                    if r2 == 0.:
                        r2 = np.nan
                    print(f'{os.path.basename(fn).split('_C')[0]} {fps[1].filename.split('/')[-1]} {plevs[i]} hPa, ratio {r1/r2:.4f} {n1},{n2}; offline '+\
                          f'{r1 * fak:.4f} online {r2 * fak:.4f}')
                
                
                readict[k][dtype][p] = {}
                readict[k][dtype][p]['attribs']=tups[4]

    try:
        #idz = np.where(np.isin(z/100, fdict['level'][fdict['pidx']]))[0]
        
        bc = out_fb['biascorr@body'][idz]
        bc[np.isnan(bc)] = 0.
        o_minus_bg=out_fb['fg_depar@body'][idz] + bc
        o_minus_an=out_fb['an_depar@body'][idz] + bc

        readict['era5fb']={'ftype':['an','fc']}                  
        readict['era5fb']['an']={'refvalues':obs[idz]-o_minus_an }                  
        readict['era5fb']['fc']={'refvalues':obs[idz]-o_minus_bg }
        readict['obs'] = {'z_coordinate': z[idz],'obstype': obstype[idz],'date_time': out['date_time'][idz], 'obs': obs[idz],}
        readict['obstypes'] = obstypes
        readict['obsunits'] = obsunits
        readict['idz'] = idz # index of standard pressure levels
    except:
        return readict
    
    if False:
        for k,v in readict.items():
            if 'obs' in k or k == 'idz':
                continue
            for ftype in v['ftype']:  
                if 'fc' in ftype:
                    dtype='fc'
                else:
                    dtype='an'
    
                if dtype in readict[k].keys():  
                    for p in readict[k][dtype].keys():
                        if p == 'refvalues':
                            continue
                        ts = out['date_time']
                        plevs = 10000, 92500
                        for plev in plevs:
                            idx=np.where(np.logical_and(obstype[idz]==obstypes[p],z[idz]==plev))
                            if len(idx[0]) > 0:
        
                                try:
        
                                    plt.subplot(2,1,1)
                                    plt.plot(1900+ts[idz][idx]/365.25/86400,obs[idz][idx],
                                             label='obs mean:{:5.3f} rms= {:5.3f}'.format(np.nanmean(obs[idz][idx]), np.sqrt(np.nanmean(obs[idz][idx]**2))))
                                    rv = readict[k][dtype]['refvalues'][idx]
                                    plt.plot(1900+ts[idz][idx]/365.25/86400,rv,
                                             label=k+' mean:{:5.3f} rms= {:5.3f}'.format(np.nanmean(rv), np.sqrt(np.nanmean(rv**2))))
                                    
                                    tmin = int(1900 + np.floor(ts[idz][idx[0][0]]/365.25/86400)+0.1)
                                    tmax = int(1900 + np.floor(ts[idz][idx[0][-1]]/365.25/86400) + 1)
                                    plt.xlim(1900+ts[idz][idx[0][0]]/365.25/86400-0.1, 1900+ts[idz][idx[0][-1]]/365.25/86400+0.1)
                                    ax = plt.gca()
                                    #ax.set_xticks(np.arange(tmin, tmax+1))
                                    plt.title(p+', '+fn.split('/')[-1].split('_')[0]+' '+str(plev//100) +' hPa')
                                    plt.legend()
                                    #plt.xlim(1935, 1965)
                                    plt.subplot(2,1,2)
                                    plt.plot(1900+ts[idz][idx]/365.25/86400,obs[idz][idx]-rv,
                                             label='offline mean:{:5.3f} rms= {:5.3f}'.format(np.nanmean(obs[idz][idx]-rv), np.sqrt(np.nanmean((obs[idz][idx]-rv)**2))))
                                    if np.any(~np.isnan(o_minus_an[idx]+rv)):
                                        plt.plot(1900+ts[idz][idx]/365.25/86400,o_minus_an[idx]+rv-rv,
                                                 label='online mean:{:5.3f} rms:{:5.3f}'.format(np.nanmean(o_minus_an[idx]+rv-rv), np.sqrt(np.nanmean((o_minus_an[idx]+rv-rv)**2))))
                                    plt.title(p+', obs -'+k )
                                    yscale = 10. ** np.floor(np.log10(np.nanmax(np.abs(obs[idz][idx]-rv)))) + 1
                                    if np.isfinite(yscale):
                                        
                                        plt.ylim(-yscale, yscale)
                                    
                                    plt.xlim(1900+ts[idz][idx[0][0]]/365.25/86400-0.1, 1900+ts[idz][idx[0][-1]]/365.25/86400+0.1)
                                    
                                    ax = plt.gca()
                                    #ax.set_xticks(np.arange(tmin, tmax+1))
                                    plt.legend()
                                    plt.tight_layout()
                                    fnp=fn.split('/')[-1].split('CEUAS_merged')[0]
                                    plt.savefig(outdir +fnp+p+'_'+k+dtype+'_'+str(plev//100)+'.png')
                                    plt.close()
        
                                except MemoryError as e:
        
                                    print('plotting reference', e)


    readict['tslice'] = tslice
    ptime(f'{os.path.basename(fn)},{cyear}, feedback calculated', tt)


    return readict

def dim_attach(g, k):
    for v in g[k].keys(): #var_selection:
        l=0            
        try:
            fvv=g[k][v]
            if 'string' not in v and v!='index':                    
                g[k][v].dims[l].attach_scale(g[k]['index'])
                #print(v,fvv.ndim,type(fvv[0]))
                if fvv.ndim==2 : #or type(fvv[0]) in [str,bytes,np.bytes_]:
                    #print(type(fvv[0]), fvv[0])
                    slen=fvv.shape[1] #sdict[v]
                    #slen=10
                    g[k][v].dims[1].attach_scale(g[k]['string{}'.format(slen)])
        except MemoryError as e:
            print(g.filename.split('/')[-1],k, e)
            pass

def read_station_height(fn, jj, inventory_list, pd_list=[]):
    heights = np.zeros(1, dtype=np.float32)
    wigos = fn.split('/')[-1].split('_')[0]
    
    if not pd_list:
        for fns in inventory_list:
            try:
                
                pd_list.append(pd.read_csv(fns, delimiter='\t'))
                print(fns.split('/')[-1], len(pd_list[-1]))
                matches = pd_list[-1].index[pd_list[-1]['primary_id']==wigos].tolist()
                #matches = pd_list[-1].index[pd_list[-1]['StationId']==wigos].tolist()
                for l in matches:
                    #heights[:] = pd_list[-1]['Hha'][l]
                    #print(heights[0])
                    #print(pd_list[-1]['start_date'][l], pd_list[-1]['elevation'][l])
                    print(f'{wigos} height: {pd_list[-1]['elevation'][l]}')
                    heights[:] = pd_list[-1]['elevation'][l]
            except Exception as e:
                print(fns, e)
    else:
        for pd_i in pd_list:
            try:
                
                #pd_list.append(pd.read_csv(fns, delimiter='\t'))
                #print(fns.split('/')[-1], len(pd_list[-1]))
                matches = pd_i.index[pd_i['primary_id']==wigos].tolist()
                #matches = pd_list[-1].index[pd_list[-1]['StationId']==wigos].tolist()
                for l in matches:
                    #heights[:] = pd_list[-1]['Hha'][l]
                    #print(heights[0])
                    #print(pd_i['start_date'][l], pd_i['elevation'][l])
                    print(f'{wigos} height: {pd_list[-1]['elevation'][l]}')
                    heights[:] = pd_i['elevation'][l]
            except Exception as e:
                print(fns, e)
    return pd_list, heights
    #fpattern=path_to_gridded+'/era5fct.{}{:0>2}.130.nc'
    #func=partial(offline_fb,fpattern,lat,lon)
    #tfgs=list(map(func,yms))

def convert_missing(refs, pd_list, gactor, wpath,cyear, fn,
                    compfilter={'compression': 32015, 'compression_opts': (3,)},
                    #compfilter={'compression': 'gzip', 'compression_opts': 3},
                    record_timestamp=False, drift=False):

    if type(refs) ==  ray._raylet.ObjectRef:
        refs = ray.get(refs)
        pd_list = ray.get(pd_list)

    tt=time.time()
    if type(compfilter) is not dict:
        compfilter = {}
    print(fn.split('/')[-1], cyear,'filter:', compfilter, record_timestamp, end=' ' )
    wpathy = wpath+'/' + str(cyear) + '/'
    targetfile = wpathy + fn.split('/')[-1].split('_')[0] +'_'+ '_'.join(fn.split('/')[-1].split('_')[2:])
    logfile = wpathy+'/log/'+fn.split('/')[-1]+".txt"
    #os.remove(logfile)
    tfile = wpath + '2022/0-20000-0-96207_CEUAS_merged_v1.nc'
    debug = False
    if os.path.isfile(logfile):
        print('already processed')
        wtime = os.path.getmtime(logfile) 
        rtime = os.path.getmtime(fn) #fn
        if rtime>wtime:
            print(logfile, 'older than some input, processing')
        else:
            return targetfile
    else:
        print('processing...')
    sys.path.insert(0,os.getcwd()+'/../resort/rasotools-master/')
    
    readictorig={'era5':{'ftype':('t','fct'),'param':{'t':'130','u':'131','v':'132','q':'133','z':'129'},
                     'path':os.path.expandvars('/mnt/scratch/scratch/leo/scratch/era5/gridded/'),'prefix':{'t':'era5t','fct': 'era5fc.0.25t'},'suffix':'','glue':'.'},
             #'CERA20C':{'ftype':('t',),'param':{'t':'130','u':'131','v':'132','q':'133','z':'129'},
               #'path':os.path.expandvars('$RSCRATCH/CERA20C/'),'prefix':'CERA20C','suffix':'','glue':'.'},
                          #'JRA55':{'ftype':('fcst_mdl',),'param':{'t':'011_tmp','u':'033_ugrd','v':'034_vgrd','q':'051_spfh'},
               #'path':os.path.expandvars('$RSCRATCH/JRA55/split/'),'prefix':'test.','suffix':'grb','glue':'.'},
                        '20CRv3':{'ftype':('t',),'param':{'t':'TMP','u':'UGRD','v':'VGRD'},#,'q':'SPFH'},#,'z':'HGT'},
                         'path':os.path.expandvars('$RSCRATCH/20CRv3/'),'prefix':{'t': 'anl_meant' },'suffix':'_pres','glue':'_'},
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
                    ref = datetime(1900, 1, 1)
                    tstart = int((datetime(cyear, 1, 1) - ref-timedelta(seconds=7200)).total_seconds())
                    tstop = int((datetime(cyear+1, 1, 1) - ref-timedelta(seconds=10799)).total_seconds())  # records between 31 12 21:00:01 and 31 12 22:00 will get lost. 
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

    
    inventory_list = ['../../meta/inventory_comparison_2/code/station_configuration/CUON_station_configuration_extended.csv']
    #inventory_list = glob.glob('../../meta/inventory_comparison_2/code/station_configuration/CUON_station_configuration_extended.csv')
    dum, station_elevation = read_station_height(fn,100, inventory_list, pd_list=pd_list)
    ptime('read inventory', tt, debug=debug)

     # wpath+fn.split('/')[-1]
    try:
#        os.remove(targetfile)
        pass
    except:
        pass

    try:
#        os.remove(logfile)
        pass
    except:
        pass
    
    try:
        
        out_name = wpathy + fn.split('/')[-1]  
        cyear = int(out_name.split('/')[-2])
        outdir = os.path.dirname(out_name) + '/'
        wigos = os.path.basename(fn).split('_')[0]
        with h5py.File(fn,'r') as f:
            try:
                if f['recordtimestamp'].shape[0] == 1:                
                    ts=f['recordtimestamp'][:]
                else:
                    ts=f['recordtimestamp'][[0, -1]]

                ref = datetime(1900, 1, 1)
                tstart = int((datetime(cyear, 1, 1) - ref).total_seconds()-3*3600)         ###
                tstop = int((datetime(cyear+1, 1, 1) - ref).total_seconds()-1)      ### records on 31 12 between 21:00 and 22:00 will get lost
                try:
                    with h5py.File(f'_{cyear + 1}_'.join(fn.split(f'_{cyear}_'))) as g:
                        tsn = g['recordtimestamp'][0]
                except:
                    tsn = tstop + 7200
                if(tsn < tstop + 7200):
                    tstop -= 10800
                if tstop < ts[0] or tstart > ts[-1]:
    #                print(fn, cyear, 'year missing in obs records')
                    try:
                        os.remove(targetfile)
                        print('removed', targetfile)
                    except:
                        pass
                    try:
                        dfile = targetfile.split('CEUAS')[0]+'.pkl'
                        os.remove(dfile)
                        print('removed', dfile)
                    except:
                        pass
                    return
                if(tstart - ts[0] > 86400 * 1):
                    print('wrong time, copying to corrupt', ref+timedelta(seconds=int(tstart)), ref+timedelta(seconds=int(ts[0])), fn)
                    os.makedirs(os.path.dirname(fn)+'/corrupt', exist_ok=True)
                    shutil.copyfile(fn, os.path.dirname(fn)+'/corrupt/'+os.path.basename(fn))
                    return
                               
                tsu=f['recordtimestamp'][:]
                ri = f['recordindex'][:]
                
                sta, sto = np.searchsorted(tsu, (tstart, tstop))
                rslice = slice(sta, sto)
                
                i = 0
                if sto == ri.shape[0]:
                    i = 1
                    tslice =slice(ri[sta], f['observations_table']['date_time'].shape[0])
                else:
                    tslice =slice(ri[sta], ri[sto])
        
                if tslice.stop ==tslice.start:
                    try:
                        os.remove(targetfile)
                        print('removed', targetfile)
                    except:
                        pass
                    try:
                        dfile = targetfile.split('CEUAS')[0]+'.pkl'
                        os.remove(dfile)
                        print('removed', dfile)
                    except:
                        pass
                    return                
                
                header_record = f['header_table']['record_timestamp'][rslice]
                #header_report = f['header_table']['report_timestamp'][rslice]
                so = f['source_configuration']['source_file'][rslice].view('S200')
                imask = np.zeros_like(header_record, dtype=bool)
                idx = np.flatnonzero(np.core.defchararray.find(so, b'igra' )!=-1)
                imask[idx] = True


                ts = f['observations_table']['date_time'][tslice]
                if sto == ri.shape[0]:
                    ri = np.concatenate((ri, [f['observations_table']['date_time'].shape[0]]))
                    tsu = tsu[sta:]
                else:
                    ri = f['recordindex'][sta:sto + 1]                
                    tsu = tsu[sta:sto]
                    
            except Exception as e:
                print(fn, cyear, e)
                return 
        #tslice = readict['tslice']
        fmatch = glob.glob(os.path.expandvars('$RSCRATCH/UH/CUON_HARVEST4/harvest_regular/igra2/')+\
                       wigos+'/'+wigos+'_'+str(cyear)+'_igra2*.nc')
        #fmatch = glob.glob('/mnt/scratch/scratch/federico/HARVEST_YEARLY_10NOV2023_Vienna/igra2/'+\
                       #wigos+'/'+wigos+'_'+str(cyear)+'_igra2*.nc')
        best_estimate = 1800 # WMO recommendation for launch time before nominal time
        tsunew = np.copy(tsu)
        print(cyear, wigos, fmatch)
        if len(fmatch) > 0:
            
            with h5py.File(fmatch[0],'r') as g:
                try:
                    its = g['recordtimestamp'][:]  # launch time in IGRA
                    hts = g['header_table']['record_timestamp'][:]  # nominal time in IGRA
                except:
                    
                    print('x')
                    pass
            # find matching report_timestamps
            imatch =np.searchsorted(hts, tsu)
            imatch[imatch==its.shape[0]] -= 1
            mask =tsu == hts[imatch]
            mask1 = tsu + 3600 == hts[imatch]  # sometimes IGRA and ERA5 records are shifted by 1 hour
            mask2 = tsu + 7200 == hts[imatch]  # sometimes IGRA and ERA5 records are shifted by 1 hour
            maskm1 = tsu - 3600 == hts[imatch-1]  # sometimes IGRA and ERA5 records are shifted by 1 hour
            
            print(os.path.basename(fn), 'IGRAshift', np.sum(mask), np.sum(mask1), np.sum(mask2), np.sum(maskm1))
            #itm = its[imatch][mask]
            #htm = hts[imatch][mask]
            # if IGRA and merged timestamps agree, estimate difference to release_time, if it is >0 or less than 3 hours
            dtm = np.zeros_like(tsu)
            dtm[mask] = hts[imatch][mask] - its[imatch][mask]
            dtm[mask1] = hts[imatch][mask1] - its[imatch][mask1] - 3600
            mask = mask | mask1 
            
            dmask = (dtm > -10800) & (dtm < 10800)
            if not np.any(dtm!=0):               
                dtm[:] = best_estimate

            mask = dmask & mask & ~imask
            if np.sum(mask) > 5 and np.any(hts-its) != 0 :
                
                best_estimate = np.int64(np.median(dtm[mask]))

            tsunew[mask] -= dtm[mask] # very important variable, used to modify record_timestamp and date_time
                        
            tsunew[ ~mask & (header_record==tsunew) & (tsunew%3600==0) & (ri[1:]-ri[:-1]<5000)] -= best_estimate 
        else:
                
            tsunew[(tsunew%3600==0) & (ri[1:]-ri[:-1]<5000)] -= best_estimate # very important variable, used to modify record_timestamp and date_time

        #launch time adjustment could introduce nonuniqueness or non-sortedness of timestamp. This loop corrects this
        idz =[0]    
        while len(idz) > 0:               
            tsdiff = tsunew[1:]-tsunew[:-1]
            idz = np.where(tsdiff<=0)[0]
            if len(idz) > 0:
                print(f'{wigos} fixed {len(idz)} timestamps')
                tsunew[idz+1] += -tsdiff[idz] +60 
    
        with eua.CDMDataset(fn) as data:
            keys = data.observations_table.keys()
            latitude = data.observations_table.latitude[0]
            longitude = data.observations_table.longitude[0]
            statid = fn.split('/')[-2]
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
                fbkeys = ['an_depar@body', 'fg_depar@body','biascorr@body']
                print(e)
                ofb=False
    
    
           # loading data:
            loaded_data=[]
            a_loaded_data=[]
            loaded_type = {'names':[],'formats':[]}
            ld=[]
            ptime('before memest', tt,debug=debug)
            a_mem = 500000
            rmin = 20

            try:
                
                a_mem = a_mem + memest3(ts, data.observations_table.z_coordinate[tslice])
            except Exception as e:
                raise ValueError(fn, e)
                
            ptime('after memest', tt,debug=debug)    
            #addmem=10000000
            okmin = ['observed_variable','observation_value','observation_id','z_coordinate','z_coordinate_type','date_time','conversion_flag','conversion_method', 'units']
            oktype = [np.int32,np.float32,data.observations_table['observation_id'][0].view('S{}'.format(data.observations_table['observation_id'].shape[1])).dtype,
                      np.float32,np.int32,np.int64,np.int32,np.int32, np.int32]
            okmin = ['observed_variable','observation_value','z_coordinate','z_coordinate_type','date_time','conversion_flag','conversion_method', 'units']
            oktype = [np.int32,np.float32, np.float32,np.int32,np.int64,np.int32,np.int32, np.int32]
            for o in okmin:
                if o not in obskeys:
                    raise ValueError(fn, 'obskey',o,' missing')
            
            struc =True    
            if struc:
                ld = list(zip(okmin, oktype))
                loaded_obstab = np.empty(ts.shape[0], dtype=ld)
            
            try:
                ptime('before load obstable', tt,debug=debug)
                for o in obskeys:
                    if o in okmin:  
                        if len(data.observations_table[o].shape)==1:
                            with np.errstate(invalid='raise'):
                                osl = data.observations_table[o][0]
                                #print(o, osl, oktype[okmin.index(o)])
                                if oktype[okmin.index(o)] in (np.int32,np.int64) and osl.dtype not in (np.int32, np.int64) :#and np.any(np.isnan(osl)):
                                    nancheck(loaded_obstab[o], data.observations_table[o][tslice], 5, nanlist[3])
                                                                
                                else:
                                    pfilla(loaded_obstab[o], data.observations_table[o][tslice], 5)
                            ptime(o+' a', tt,debug=debug)
                        else:
                            pfilla(loaded_obstab[o],data.observations_table[o][tslice].view('S{}'.format(data.observations_table[o].shape[1]))[:, 0] , 5)
                            ptime(o, tt,debug=debug)
            except MemoryError as e:
                print(fn, e)
                #return None
                raise ValueError(fn)
            
            # fix wrong pressure scaling for some hara files - max be removed when harvest is fixed
            
            struc = True
            ptime('after load obstable', tt,debug=debug)
            try:

                if record_timestamp:
                    
                    mustsort = np.all(tsunew[:-1] < tsunew[1:]) # sorting of all variables needed!
                    
                    #for i in range(tsunew.shape[0]):
                        #loaded_obstab['date_time'][ri[i]:ri[i+1]] = tsunew[i]
                    #ts[:] = loaded_obstab['date_time'][:]
                    if mustsort:
                        print(time.time()-tt, fn.split('/')[-1],'MUSTSORT true - date_time modified!')
                        x = 0
                        
            except Exception as e:
    
                #for l in loaded_data:
                    #print(l.shape)
                print(fn, e)
                return None            
            
            if ts.shape[0] < 10000:
                sl = ts.shape[0]
            else:
                sl = ts.shape[0] // 20
            if np.nanmean(loaded_obstab['z_coordinate'][:sl][loaded_obstab['z_coordinate_type'][:sl] == 1]) > 110000:
                if wigos in ('0-20000-0-37717', '0-20000-0-43599'):
                    print(f'WARNING {wigos} {cyear} spurious z')
                    return targetfile
                else:
                    print(f'WARNING {wigos} {cyear} spurious z')
                    return targetfile
                    #raise ValueError(fn, 'hara file fixed')
                loaded_obstab['z_coordinate'][loaded_obstab['z_coordinate_type'] == 1] /= 100
            ptime('after hara', tt,debug=debug)
            
                
            if struc:
                dr = []
                if drift:
                    dr =  [('latd',np.float32),('lond',np.float32),('timed',np.int32)]

                a_loaded_obstab = np.empty(a_mem, dtype=ld+dr)
                for n in dr:
                    pfill(a_loaded_obstab[n[0]], 0, 5)

            ptime('after drift arrays', tt,debug=debug)
    
            loaded_fb=[]
            a_loaded_fb=[]
            loaded_type = {'names':[],'formats':[]}
            lf=[]
            fbmin = ['fg_depar@body','an_depar@body','biascorr@body']
            for o in fbmin:
                if o not in fbkeys:
                    raise ValueError(fn, 'fbkey',o,' missing')
            
            if struc:
                lfb = list(zip(fbmin+['fg_depar@offline'], [np.float32]*5))
                loaded_feedback = np.empty(ts.shape[0], dtype=lfb)
                a_loaded_feedback = np.empty(a_mem, dtype=lfb)

            for o in fbkeys + ['fg_depar@offline']:
                if o in fbmin + ['fg_depar@offline']:
                    if ofb:
                        if o == 'fg_depar@offline':
                            #loaded_fb.append((data.era5fb['fg_depar@body'][tslice]))
                            #loaded_fb.append(np.empty_like(data.era5fb['fg_depar@body'][tslice], dtype=np.float32))
                            ptime(f'after {o} feedback', tt,debug=debug)
                            pfill(loaded_feedback[o], np.nan, 5)
                            #loaded_feedback[o].fill(np.nan)
                        else:    
                            pfilla(loaded_feedback[o] , data.era5fb[o][tslice] , 5)#).astype(np.float32))
                            loaded_feedback[o][loaded_feedback[o] == -2147483600.0] = np.nan
                    else:
                        try:
                            
                            pfill(loaded_feedback[o] , np.nan, 5)
                        except Exception as e:
                            print(os.path.basename(fn), e)
                            return
                                
                    ptime(f'before {o} a_feedback', tt,debug=debug)
                    pfill(a_loaded_feedback[o], np.nan,5)
                    ptime(f'after {o} a_feedback', tt,debug=debug)

            ptime('after recarray feedback', tt,debug=debug)
            
            ov = loaded_obstab['observation_value']
            ot = loaded_obstab['observed_variable']
            try:
                
                elsum,rhsum = check_spurious(ov, ot)
                if elsum > 0:               
                    print(f'{statid},{cyear}, {elsum} spurious obs values eliminated')
            except Exception as e:
                raise ValueError(fn, e)
            if rhsum > 5:
                print(f'WARNING: {rhsum} RH values above 1 have been found and corrected, {fn}')
            
                
                
            
            ptime('after spurious value elimination', tt,debug=debug)
        
            humvar=np.array((ipar[34],ipar[36],ipar[38],ipar[39])) #dpd,dp,rh,sh
            wvar=np.array((ipar[104],ipar[105],ipar[106],ipar[107])) #dpd,dp,rh,sh
            hvar = np.array([117]) # geopotential (if pressure is not present)
            ps = np.array((10., 20., 30., 50., 70., 100, 150, 200, 250, 300, 400, 500, 700, 850, 925, 1000)) *100.
        
            ptime('before augment', tt, True)
            #a_loaded_obstab['units'] = np.zeros_like(a_loaded_obstab['observed_variable'])
            
            try:
                
                slist = [ np.where(d==b'/')[0] for d in data.source_configuration.source_file]
                dlist = []
                m = 0
                for s in slist:
                    dlist.append(data.source_configuration.source_file[m, s[-3]+1:s[-2]])
                    dlist[-1] = dlist[-1].view(f'S{dlist[-1].shape[0]}')
                    m += 1
                dlist =tuple(dlist)    
                tup =augment2(loaded_obstab, a_loaded_obstab, loaded_feedback, a_loaded_feedback,
                                                      ri,ts, ps, humvar,wvar, hvar, fn)
            except Exception as e:
                print('ERROR', fn, e)
                raise ValueError(fn, loaded_obstab.shape[0], a_loaded_obstab.shape[0])
                    
            out, fb_out, jj,addedvar,tddmax = tup
            ptime(f'after augment {loaded_obstab.shape[0]}, {a_loaded_obstab.shape[0]}, {jj}', tt)

            try:
                
                elsum,rhsum = check_spurious( out['observation_value'], out['observed_variable'],wind=True)
                if elsum > 0:               
                    print(f'{statid},{cyear}, {elsum} spurious obs values eliminated')
            except MemoryError as e:
                raise ValueError(fn, e)
            if rhsum > 5:
                print(f'WARNING: {rhsum} RH values above 1 have been found and corrected, {fn}')

            #print(np.sum((out['observed_variable']==106)), f'tddmax: {tddmax:.3f}')
            if out.shape[0] == 0:
                print('WARNING', fn, 'no data after augment')
                return targetfile
                        
            _, ori = np.unique(out['date_time'], return_index=True)
            
            for a, b, c in zip(np.arange(ori.shape[0]), ts[ri[:-3]],out['date_time'][ori] ):
                if b != c:
                    print(a, b, c)
                    raise ValueError(fn, a, b, c)
            
            #ori = np.concatenate((ori, [out['date_time'].shape[0]]))
            
            if False:
                ari = np.where(out['date_time'][1:]-out['date_time'][:-1]>0 )
                ari = np.concatenate(([0], ari[0]))
                for i in range(len(ari)-1):
                    sl = slice(ari[i], ari[i+1])
                    tdidx = np.where((out['observed_variable'][sl]==137))
                    if len(tdidx[0]) == 0:
                        continue
                    tidx = np.where((out['observed_variable'][sl]==126))
                    qidx = np.where((out['observed_variable'][sl]==39))
                    rhidx = np.where((out['observed_variable'][sl]==138))
                    gpl = out['z_coordinate'][sl][rhidx]
                    tpl = out['z_coordinate'][sl][tidx]
                    g = np.searchsorted(tpl, gpl)
                    q = out['observation_value'][sl][qidx]
                    qc = out['conversion_flag'][sl][qidx]
                    qcidx = np.where(qc==0)
                    rh = out['observation_value'][sl][rhidx]
                    rhc = out['conversion_flag'][sl][rhidx]
                    rhcidx = np.where(rhc==0)
                    td = out['observation_value'][sl][tdidx]
                    tdc = out['conversion_flag'][sl][tdidx]
                    tdcidx = np.where(tdc==0)
                    t = out['observation_value'][sl][tidx][g]
                    
                    #print(f'{np.corrcoef(t-td, rh)[0][1]:.4f}, {np.corrcoef(rh[rhcidx], (t-td)[rhcidx])[0][1]:.4f}', len(tdcidx[0]), len(rhcidx[0]))
                    #plt.plot(rh, t-td, 'o')
                    #plt.plot(rh[rhcidx], (t-td)[rhcidx], 'o')
                    #plt.show()
                plev = 20000
                td = []
                rh = []
                t = []
                tdc = []
                rhc = []
                tc = []
                dt = []
                q = []
                qc = []
                for i in range(len(ari)-1):
                    sl = slice(ari[i], ari[i+1])
                    tdidx = np.where((out['observed_variable'][sl]==137)&(out['z_coordinate'][sl]==plev))
                    if len(tdidx[0]) == 0:
                        continue
                    rhidx = np.where((out['observed_variable'][sl]==138)&(out['z_coordinate'][sl]==plev))
                    tidx = np.where((out['observed_variable'][sl]==126)&(out['z_coordinate'][sl]==plev))
                    qidx = np.where((out['observed_variable'][sl]==39)&(out['z_coordinate'][sl]==plev))
                    dt.append(out['date_time'][sl][tdidx])
                    td.append(out['observation_value'][sl][tdidx])
                    t.append(out['observation_value'][sl][tidx])
                    q.append(out['observation_value'][sl][qidx])
                    qc.append( out['conversion_flag'][sl][qidx])
                    rh.append(out['observation_value'][sl][rhidx])
                    rhc.append( out['conversion_flag'][sl][rhidx])
                    tc.append( out['conversion_flag'][sl][tidx])
                    tdc.append( out['conversion_flag'][sl][tdidx])
                
                rh =np.array(rh).flatten()
                rhc = np.array(rhc).flatten()
                q =np.array(q).flatten()
                qc = np.array(qc).flatten()
                t = np.array(t).flatten()
                td = np.array(td).flatten()
                tdc = np.array(tdc).flatten()
                tc = np.array(tc).flatten()
                dt = np.array(dt).flatten()
                print(f'{np.corrcoef(t-td, rh)[0][1]:.4f}, {np.corrcoef(rh[rhc==0], (t-td)[rhc==0])[0][1]:.4f}', np.sum(tdc==0), np.sum(rhc==0))
            
            #plt.subplot(1, 2, 1)
            #plt.plot(rh, t-td, 'o')
            #plt.plot(rh[rhc==0], (t-td)[rhc==0], 'o')
            #plt.subplot(1, 2, 2)
            #plt.plot(rh, q, 'o')
            #plt.plot(rh[rhc==0], q[rhc==0], 'o')
            #plt.show()

                
            allvar = np.concatenate((np.array([117, 126]), humvar , wvar))
            
            os.makedirs(wpathy, exist_ok=True)
            path_to_gridded=os.path.expandvars(rscratch+'/era5/gridded/')
            try:
                reaname = os.path.expandvars(wpathy+fn.split('/')[-1].split('_CEUAS_merged')[0] + '.pkl')
                wtime = os.path.getmtime(reaname)
                rtime = os.path.getmtime(fn) #fn
                if rtime>wtime:
                    print(reaname, 'older than some input')
                    raise ValueError
                with open(reaname+'x','rb') as f:
                    
                    readict=retrieve_anfg(data.file, out,fb_out, copy.copy(readictorig),ts, tsunew, tslice, refs, gactor, out_name,path_to_gridded,
                                          drift=True, station_elevation=station_elevation, ri=ri, ori=ori, driftonly=True)
                    readict=pickle.load(f)
                    readict['lat'] = latitude
                    readict['lon'] = longitude
                    readict['statid'] = statid
                    with open(reaname,'wb') as f:
                        pickle.dump(readict,f)
                    if 'tslice' not in readict.keys():
                        print(reaname, 'something is wrong with readict, recreating...')
                        raise ValueError
                    ptime (f'{wigos} read fb',tt)
            except:
        
                out_name = wpathy + fn.split('/')[-1]
        
                #print(np.sum(out['date_time']))
                readict=retrieve_anfg(data.file, out,fb_out, copy.copy(readictorig),ts, tsunew, tslice, refs, gactor,
                                      out_name,path_to_gridded, drift=True, station_elevation=station_elevation, ri=ri, ori=ori)
                #print(np.sum(out['date_time']))
        
                if readict is None:
                    print(fn.split('/')[-1], 'no data for year', cyear, time.time() - tt)
                    return
                with open(reaname,'wb') as f:
                    pickle.dump(readict,f)
                    
                ptime (f'{wigos} created and wrote fb',tt)
            
        
            if(ofb):
                if(np.any(loaded_feedback['fg_depar@body']>1.e26) or np.any(readict['era5']['fc']['refvalues']>1.e26)):
                    print(cyear)
                ttt = time.time()
                try:
                    
                    add_fb(out,fb_out,readict['20CRv3']['an']['refvalues'],
                       readict['era5']['an']['refvalues'],readict['era5']['fc']['refvalues'], readict['idz'])
                except MemoryError as e:
                    print(e)
                    raise ValueError(fn)
                print('addfb', time.time()-ttt)
                
                if False:
                    ari = np.where(out['date_time'][1:]-out['date_time'][:-1]>0 )
                    ari = np.concatenate(([0], ari[0]))
                    dcalcl = {};dobsl = {};dfgdepl = {};drefl = {}
                    for v in 39, 137, 138, 126:
                        dcalcl[v] = np.full((len(ari), ps.shape[0]), np.nan)
                        dobsl[v] = dcalcl[v].copy() 
                        drefl[v] = dcalcl[v].copy() 
                        dfgdepl[v] = dcalcl[v].copy()
                    
                    pplot1 = False
                    pplot = True
                    nempty = 0
                    ttt = time.time()
                    for i in range(len(ari)-1):
                        sl = slice(ari[i], ari[i+1])
                        obsv = out['observed_variable'][sl]
                        if not np.any(obsv==138):
                            nempty += 1
                            continue
                        obs = out['observation_value'][sl]
                        z = out['z_coordinate'][sl]
                        fgdep = fb_out['fg_depar@body'][sl]
                        ref = readict['era5']['fc']['refvalues'][sl]
                        
                        didx = {}; dobs = {}; dfgdep = {};dref = {}; dcalc = {}; dz = {}
                        for v in 39, 137, 138, 126:
                            didx[v] = np.where((obsv==v))[0]
                            dz[v] = z[didx[v]]
                            if v == 126:
                                idy = np.searchsorted(dz[v], dz[138])
                                didx[v] = didx[v][idy]
                                dz[v] = dz[v][idy]
                            dobs[v] = obs[didx[v]]
                            dfgdep[v] = fgdep[didx[v]]
                            dref[v] = ref[didx[v]]
                        
                        #vpdata = dobs[138]* np.exp(liquid(dobs[126])) #* Sonntag(fobs[mask, 2])# * np.exp(liquid(fobs[mask, 2]))  #Sonntag(fobs[mask, 2])
                    #vpdata2 = fobs[mask, 5]* np.exp(ice(fobs[mask, 2])) #* Sonntag(fobs[mask, 2])# * np.exp(liquid(fobs[mask, 2]))  #Sonntag(fobs[mask, 2])
                    #tm = fobs[mask, 2] < 213.15
                    #vpdata[tm] = vpdata2[tm]
                        dcalc[137] = np.empty_like(dobs[137])
                        #l = 0
                        #for il in range(len(dcalc[137])):
                            #if il < dobs[126].shape[0]:
                                
                                #dcalc[137][il] = bisection(func, dobs[126][il]-70., dobs[126][il]+1., vpdata[il])
                            #else:
                                #print(i,dobs[126].shape[0],dobs[137].shape[0] )
                        
                        for j in range(len(ps)):
                            
                            try:
                                
                                si= np.searchsorted(dz[137], ps[j])
                                if si < dz[137].shape[0]:                          
                                    if dz[137][si] == ps[j]:
                                        for v in 39, 137, 138, 126:
                                            dobsl[v][i, j]=dobs[v][si]
                                            drefl[v][i, j]=dref[v][si]
                                            dfgdepl[v][i, j]=dfgdep[v][si]
                                            if v == 137:                                       
                                                vpdata = dobs[138][si]* np.exp(liquid(dobs[126][si])) #* Sonntag(fobs[mask, 2])# * np.exp(liquid(fobs[mask, 2]))  #Sonntag(fobs[mask, 2])
                                                dcalc[137][si] = bisection(func, dobs[126][si]-70., dobs[126][si]+1., vpdata)
                                                dcalcl[v][i, j]=dcalc[v][si]
                                                
                                                if dlist[i] == b'igra2':
                                                    print(ps[j], dfgdepl[v][i, j],dcalcl[v][i, j],dobsl[v][i, j] )
                                                    dcalcl[v][i, j] = np.nan
                                                if np.abs(dcalcl[v][i, j]-dobsl[v][i, j]) > 10:
                                                    x = 0
                                                
                            except:
                                pass
                                            
                        if pplot1:
                            z100 = np.searchsorted(dz[137], 10000)
                            if np.abs(dfgdep[137][z100]) < 2 and  np.abs(dobs[137][z100]-dcalc[137][z100]) > 10.:
                                x = 0
                            idx =np.where(~np.isnan(dfgdep[137]))[0]
                            plt.subplot(1, 3, 1)    
                            plt.semilogy(dobs[137][idx]-dcalc[137][idx],dz[137][idx]/100.,'.' )
                            plt.subplot(1, 3, 2)    
                            plt.semilogy(dfgdep[137][idx],dz[137][idx]/100.,'.' )
                            plt.subplot(1, 3, 3)    
                            plt.semilogy(dcalc[137][idx]-dref[137][idx],dz[137][idx]/100.,'.' )
                    if pplot1:
                        for sp in range(3):
                            
                            plt.subplot(1, 3, sp+1)    
                            plt.ylim(1000, 10)
                            plt.ylabel('p [hPa')
                            plt.xlabel('difference[K]')
                            plt.xlim(-10, 10)
                        plt.title(f'reported-calculated dewpoint (from rh,T),{statid},{cyear}')
                        
                        plt.show()
                    print(f'calctime,{time.time()-ttt:.4f}')
                    m = 0
                    plt.figure(figsize=(15, 10))
                    for pl in 5, 6, 7, 9:
                        m += 1
                        plt.subplot(4, 1, m)
                        plt.plot(out['date_time'][ari]/86400/365.25, dcalcl[137][:, pl]-drefl[137][:, pl],'.' , label=f'calc-ref,{rms(dcalcl[137][:, pl]-drefl[137][:, pl]):.2f}, #{np.sum(~np.isnan(dcalcl[137][:, pl]))}')
                        plt.title(f'Dewpoint, {ps[pl] / 100} hPa, {statid}')
                        plt.ylabel('K')
                        plt.plot(out['date_time'][ari]/86400/365.25, dobsl[137][:, pl]-dcalcl[137][:, pl],'.' , label=f'obs-calc,{rms(dobsl[137][:, pl]-dcalcl[137][:, pl]):.2f},{np.nanmean(dobsl[137][:, pl]-dcalcl[137][:, pl]):.2f}')
                        plt.ylabel('K')
                        plt.plot(out['date_time'][ari]/86400/365.25, dfgdepl[137][:, pl], '.', label=f'fgdep,{rms(dfgdepl[137][:, pl]):.2f}')
                        plt.ylabel('K')
                        plt.legend()
                    plt.tight_layout()
                    plt.savefig(f'{statid}_{cyear}_td.png')
                    plt.close()
    
    
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
    
    
        reduced_obskeys=list(loaded_obstab.dtype.fields.keys())
        reduced_fbkeys=list(loaded_feedback.dtype.fields.keys())
    
    
    
        # sorting:
        print('start sorting', time.time()-tt)
        rt = out['date_time'][ori]
        #rt, ori = np.unique(out['date_time'], return_index=True)
        #ori = np.concatenate((ori, [out['date_time'].shape[0]]))
        recordtimestamps = rt
    
        
        try:
            with h5py.File(fn, 'r') as file:
                with h5py.File(targetfile, 'w') as newfile:
        
                    groups = []
                    for i in file.keys():
                        if type(file[i]) == h5py._hl.group.Group:
                            if i not in ('observations_table','era5fb'):
                                groups.append(i)
                                newfile.create_group(i)
                                
                        elif i in ('recordindex', 'recordtimestamp', 'dateindex'):
                            pass
                        else:
                            newfile.create_dataset(i, data=file[i][:], **compfilter) #compression= 32015, compression_opts=(3,))
                    for i in groups:
                        if(i == 'recordindices' or i == 'observations_table' or i == 'era5fb'):
                            pass
                        elif i in ('header_table', 'source_configuration', 'station_configuration'):
                            if 'index' not in file[i].keys():
                                newfile[i].create_dataset('index', data=np.empty(sto-sta, dtype='S1'), **compfilter) #compression= 32015, compression_opts=(3,))                       
                            for j in file[i].keys():
                                #print(i, j)
                                if 'string' in j:
                                    xdata = file[i][j][:]
                                else:
                                    xdata = file[i][j][sta:sto]
                                    #print(xdata.shape, xdata.dtype)
                                newfile[i].create_dataset(j, data=xdata, **compfilter) #compression= 32015, compression_opts=(3,))
                            if i == 'header_table':
                                newfile[i]['height_of_station_above_sea_level'][:] = station_elevation
                                if record_timestamp:
                                    try:
                                        
                                        newfile[i]['report_timestamp'][:] = tsunew  # must be identical with date_time at surface
                                        #newfile[i]['record_timestamp'][:] = tsunew
                                    except MemoryError as e:
                                        print(e, '##ERROR##',fn, headerslice, file['header_table']['record_timestamp'][:][[0, -1]])
                                        #print(e, '##ERROR##',fn, headerslice, rt[0], rt[-1], file['header_table']['record_timestamp'][:][[0, -1]])
                                        return
                                        raise ValueError(fn)
                                                                            
                        else:
                            if 'index' not in file[i].keys():
                                if len(file[i].keys()) == 0:
                                    print('empty group', i)
                                    newfile[i].create_dataset('index', data=np.empty(1, dtype='S1'), **compfilter) #compression= 32015, compression_opts=(3,))                       
                                    continue
                                sh = file[i][list(file[i].keys())[0]].shape[0]
                                newfile[i].create_dataset('index', data=np.empty(sh, dtype='S1'), **compfilter) #compression= 32015, compression_opts=(3,))                       
                            for j in file[i].keys():
                                #var = file[i][j][:]
                                #print(var.shape, var.dtype)
                                newfile[i].create_dataset(j, data=file[i][j][rslice], **compfilter) #compression= 32015, compression_opts=(3,))
                        try:
                            
                            dim_attach(newfile, i)
                        except Exception:
                            print(i, 'dim_attach failed')
        except MemoryError as e:
            print('could not write to resorted file', e)
            raise ValueError(fn)
            return
        
        obsv = out['observed_variable'][:jj]
        allvars = np.sort(np.unique(obsv))
        absidx2 = np.lexsort((out['z_coordinate'][:jj], out['date_time'][:jj], out['observed_variable'][:jj]) )
        od2 = out['date_time'][:jj][absidx2]
        ov2 = out['observed_variable'][:jj][absidx2]
        _, idv = np.unique(ov2, return_index=True)
        idv = np.concatenate((idv, np.array([ov2.shape[0]])))
        vridx2 = []
        for j in range(len(allvars)):
            vridx2.append(np.searchsorted(od2[idv[j]:idv[j+1]], rt, side='left'))
            vridx2[-1] += idv[j]
            vridx2[-1] = np.concatenate((vridx2[-1], np.array([idv[j+1]])))
                
    
        vridx = vridx2
        absidx = absidx2
        od = od2
        ov = ov2
        
        if False:  # old method
            is_sorted = lambda a: np.all(a[:-1] <= a[1:])
            # resorting the data
            #
            #@njit(cache=True, boundscheck=False)
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
            
            j = -1
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
        
            ats = out['date_time'][:jj]
            # check integrity of indices
            od = out['date_time'][:jj][absidx]
            oov = out['observed_variable'][:jj][absidx]
            
        ref = datetime(1900, 1, 1)
        if False:  # just for debug
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
                            raise ValueError(fn)
                        
                    if any(od[vridx[j][i]:vridx[j][i + 1]]!=rt[i]):
                        try:
                            
                            print('spurious',fn, allvars[j], i, ref+timedelta(seconds=int(od[vridx[j][i - 1]])),
                                  ref+timedelta(seconds=int(od[vridx[j][i + 1] - 1])) )
                            #print('spurious',allvars[j], i, oov[vridx[j][i]],oov[vridx[j][i + 1] - 1])
                            print('spurious',fn, allvars[j], i, ooval[vridx[j][i]:vridx[j][i + 1]])
                            raise ValueError(fn)
    
                        except Exception as e:
                            raise ValueError(fn, e)
            
            print('vridx consistency:', [np.std(vridx[i]-vridx2[i]) for i in range(len(vridx))])
    
        print('elapsed converting: ',time.time()-tt)
        
        @njit(boundscheck=False)
        def inplace(a, b, idx):
            for i in range(a.shape[0]):
                a[i] = b[idx[i]]
                
            return a
        
        try:
            #ov_vars_int32 = np.empty_like(out['observed_variable'])
            ov_vars = np.empty(out['date_time'].shape[0], dtype='S20')
            ov_vars[:] = ''
            #ov_vars_float32 = np.empty_like(out['observation_value'])
            
    
            #ttx = time.time()
            with h5py.File(targetfile, 'a') as fd:
                    
                for i in obskeys + ['lond', 'latd', 'timed']:
            
                    if i == 'observation_id':
                        #ov_vars = np.empty(jj,dtype=out[i].dtype)
                        #ov_vars[:]=out[i]
                        #dim1 = int(str(out[i].dtype).split('S')[-1])
                        #ov_vars=fill_obsid(ov_vars.view('S1').reshape((len(ov_vars),dim1)),out['conversion_flag'])
                        continue
            
                    elif i in reduced_obskeys or drift and i in ('lond', 'latd', 'timed'):
                        pass
                        #ov_vars = out[i]
                    else:
                        continue
            
                    #ov_vars = ov_vars[absidx]
                    if i == 'index':
                        pass
                    elif i in ['report_id', 'sensor_id', 'source_id']:
                        #alldict = {i:ov_vars[absidx].view(dtype='S1')}
                        alldict = {i:inplace(np.frombuffer(ov_vars, out[i].dtype, out[i].shape[0]), out[i], absidx).view('S1')}
                        write_dict_h5(fd, alldict, 'observations_table', {i: compfilter  }, [i])  #{ 'compression': 32015, 'compression_opts':(3,) } 
                    else:
                        if out[i].dtype in (np.float32, np.float64):
                            dt = np.float32
                        elif i == 'date_time':
                            dt = out[i].dtype
                        else:
                            dt = np.int32
                        
                        alldict = {i:inplace(np.frombuffer(ov_vars, dt, out[i].shape[0]), out[i], absidx)}
                                       
                        write_dict_h5(fd, alldict, 'observations_table', {i: compfilter } , [i])  #{ 'compression': 32015, 'compression_opts':(3,) } 
                    #print(i, time.time()-tt)
                #del ov_vars
                del obsv
                #del dta
                outsh = out.shape[0]
           
                if True:
                    for i in fbkeys + ['fg_depar@offline']:
                        if i in reduced_fbkeys:
                            
                            #ov_vars = fb_out[i][:jj]
                            pass
                        else: 
                            continue
            
                        #ov_vars = ov_vars[absidx]
            
                        if i == 'index':
                            pass
                        elif i in ['expver', 'source@hdr', 'source_id', 'statid@hdr']:
                            #alldict = {i:np.asarray(ov_vars[absidx], dtype='S1')}
                            alldict = {i:inplace(np.frombuffer(ov_vars, fb_out[i].dtype, fb_out[i].shape[0]), fb_out[i], absidx).view('S1')}
                            write_dict_h5(fd, alldict, 'era5fb', {i: compfilter  }, [i]) #{ 'compression': 32015, 'compression_opts':(3,) }
                        else:
                            #alldict = pandas.DataFrame({i:ov_vars})
                            #alldict = {i:ov_vars[absidx]}
                            if fb_out[i].dtype in (np.float32, np.float64):
                                dt = np.float32
                            else:
                                dt = np.int32
                            alldict = {i:inplace(np.frombuffer(ov_vars, dt, fb_out[i].shape[0]), fb_out[i], absidx)}
                            
                          
                            write_dict_h5(fd, alldict, 'era5fb', {i: compfilter}, [i]) #{ 'compression': 32015, 'compression_opts':(3,) } 
            
            
                    obskeys.append('station_elevation') # fix for adding station elevation 
                    #if True:
                    with eua.CDMDataset(fn) as data:
                        for i in obskeys:
                            if i in reduced_obskeys:
                                continue
                            #if i not in ('latitude','longitude','report_id'):
                                #continue
                            if i == 'station_elevation':                   
                                    #i='z_coordinate'
                                #ov_vars = np.full(jj, station_elevation,dtype=np.float32) # add station elevation
                                write_dict_h5(fd, {i: np.full(jj, station_elevation,dtype=np.float32)}, 'observations_table', {i: compfilter  }, [i]) #{ 'compression': 32015, 'compression_opts':(3,) }
                                continue
                            elif i == 'observation_height_above_station_surface':
                                try:
                                    
                                    z = out['z_coordinate']
                                    evarno = np.full_like(out['observed_variable'], -1)
                                    evarno = fill_restdata_addedvar(evarno, data.era5fb['varno@body'][tslice], addedvar)
                                    oh = np.full_like(out['z_coordinate'], np.nan, dtype=np.float32)
                                    oh = fill_restdata_addedvar(oh, data.observations_table[i][tslice], addedvar)
                                    oh[oh!=1] = np.nan
                                    
                                    if cyear > 2019: # ncar surface information is not encoded properly after 2019, thus needs o be removed.
                                        
                                        so = data.source_configuration.source_file[rslice].view('S200')
                                        idf = np.flatnonzero(np.core.defchararray.find(so, b'ncar' )!=-1)
                                        for f in idf:
                                            if f == ori.shape[0] - 1:
                                                oh[ori[f]:] = np.nan
                                            else:    
                                                oh[ori[f]:ori[f+1]] = np.nan

                                    fmatch = glob.glob(os.path.expandvars('$RSCRATCH/UH/CUON_HARVEST4/harvest_regular/igra2/')+\
                                                   wigos+'/'+wigos+'_'+str(cyear)+'_igra2*.nc')
                                    #fmatch = glob.glob(os.path.expandvars('$RSCRATCH/UH/CUON_HARVEST4/harvest_regular/ncar/')+\
                                                   #wigos+'/'+wigos+'_'+str(cyear)+'_ncar*.nc')
                                    #print(cyear, wigos, fmatch)
                                    if len(fmatch) > 0: # search for IGRA surface information even if record was not used
                                        
                                        with h5py.File(fmatch[0],'r') as g:
                                            #try:
                                                iri = g['recordindex'][:]
                                                its = g['recordtimestamp'][:]  # launch time in IGRA
                                                hts = g['header_table']['record_timestamp'][:]  # nominal time in IGRA
                                                ohi = g['observations_table']['observation_height_above_station_surface'][:]
                                                zi = g['observations_table']['z_coordinate'][:]
                                                #idi = np.where(ohi==1.0)[0]
                                                #ti = g['observations_table']['date_time'][:]
                                            #except:
                                                #pass
                                    
                                    
                                        is_surf3(oh, z, ori, mask, maskm1, imatch, iri, zi, ohi)
                                                                                          
                                            
                                    #dt = out['date_time']
                                    zmed = np.median(z[oh==1])
                                    #oh[z<zmed-5000.] = np.nan
                                    
                                    #zdiff = is_surf2(out['z_coordinate'], out['observed_variable'], out['observation_value'], evarno, oh, ori, ps, station_elevation[0])
                                    idx = np.where((np.isin(evarno , (40, 39, 58, 281, 41, 42))))[0]
                                    #idy = np.where(zdiff==zdiff)[0]
                                    #idz = np.where((zdiff==zdiff)&(np.isin(evarno , (40, 39, 58, 281, 41, 42))))[0]
                                    oh[idx] = 1.0
                                    oh[(z<zmed-5000.) | (z>zmed+5000.)] = np.nan
                                    #ida = np.where(oh==1.0)[0]
                                    #plt.plot(idy, z[idy], '.', label=f'sv {len(idy)}')
                                    #plt.plot(dt[idx]/86400/365.25, z[idx], '.', label=f'sv {len(idx)}')
                                    #if len(fmatch) > 0: plt.plot(ti[idi]/86400/365.25, zi[idi], '.', label=f'igra {len(idi)}')
                                    #plt.plot(dt[ida]/86400/365.25, z[ida], '.', label=f'sv&igrancar {len(ida)}')
                                    #plt.legend()
                                    #plt.title(f'{wigos}_{cyear}')
                                    #plt.savefig(os.path.expandvars('$HOME/tmp/')+f'{wigos}_{cyear}_sv.png')
                                    #plt.close()
                                    
                                    #ov_vars = oh[absidx]
                                    #dt = ov_vars.dtype
                                    write_dict_h5(fd, {i: oh[absidx]}, 'observations_table', {i: compfilter  }, [i]) #{ 'compression': 32015, 'compression_opts':(3,) }
                                    continue

                                    ##p_to_z_ifs(z[ix]) /g > station_height + 100. or 
                                except MemoryError as e:
                                    raise ValueError(fn, e)
                            else:
                                try:
                                    
                                    rest_data = data.observations_table[i][tslice]
                                except Exception as e:
                                    raise ValueError(fn, e)
                                    
                                if rest_data.dtype in (np.float32, np.float64):
                                    dt = np.float32
                                elif rest_data.dtype in (np.int32, np.int64):
                                    dt = np.int32
                                else:
                                    dt = rest_data.dtype
                                    
                                pass
                                #rest_data = data.observations_table[i][tslice]
                                #if rest_data.ndim==2: #i in ['observation_id', 'report_id', 'sensor_id', 'source_id']:
                                    #ov_vars = np.frombuffer(ov_vars, out[i].dtype, out[i].shape[0]), out[i], absidx)
                                    ##ov_vars = np.empty((addedvar.shape[0],len(rest_data[0])), dtype=rest_data[0].dtype)
                                    ##ov_vars.fill('n')
                                #else:
                                    #if rest_data.dtype in (np.float64, np.float32):
                                        #ov_vars = np.empty(addedvar.shape[0], dtype=np.float32)
                                        #miss_val = np.nan
                                    #elif  rest_data.dtype == np.int64:
                                        #ov_vars = np.empty(addedvar.shape[0], dtype=np.int32)
                                        #miss_val = -2147483648
                                    #elif  rest_data.dtype == np.int32:
                                        #ov_vars = np.empty(addedvar.shape[0], dtype=np.int32)
                                        #miss_val = -2147483648
                                    #else:
                                        #raise ValueError('unknown dtype', rest_data.dtype)

                                #fill_restdata_abs(ov_vars, rest_data, ori, recordindex-tslice.start, absidx) #,out['z_coordinate'][:jj])
            
                                
                            if i == 'index' :
                                pass
                            elif i in ['observation_id']:
                                dlen = np.sum(addedvar==-1)
                                dum = np.array(np.arange(cyear*10000000000000000+9900000000000000, cyear*10000000000000000+9900000000000000+dlen), dtype = 'S20')
                                sh = dum.shape
                                dum = dum.view('S1').reshape(sh[0], 20)
                                sh = (out['date_time'].shape[0], rest_data.shape[1])
                                try:
                                    
                                    alldict = {i:fill_restdata_obsid_abs(np.frombuffer(ov_vars, dt, sh[0]*sh[1]).reshape(sh),
                                           rest_data, addedvar, absidx, dum)}
                                except:
                                    raise ValueError(fn)
                                write_dict_h5(fd, alldict, 'observations_table', {i: compfilter  }, [i]) #{ 'compression': 32015, 'compression_opts':(3,) }
                            elif i in [ 'report_id', 'sensor_id', 'source_id']:
                                try:
                                    if i == 'sensor_id' and rest_data.ndim != 2:
                                        rest_data = np.array(rest_data, dtype='S4').view('S1').reshape(rest_data.shape[0], 4)
                                        dt = np.dtype('S1')
                                        print('WARNING: no sensor_id', fn)
                                        
                                    sh = (out['date_time'].shape[0], rest_data.shape[1])
                                except:
                                    raise ValueError(fn+' '+i) 
                                #ov = fill_restdata_addedvar_abs(np.frombuffer(ov_vars, dt, sh[0]*sh[1]).reshape(sh), rest_data, addedvar, absidx, 'n')

                                alldict = {i:fill_restdata_abs(np.frombuffer(ov_vars, dt, sh[0]*sh[1]).reshape(sh),
                                           rest_data, ori, recordindex-tslice.start, addedvar, absidx)}
                                write_dict_h5(fd, alldict, 'observations_table', {i: compfilter  }, [i]) #{ 'compression': 32015, 'compression_opts':(3,) }
                            else:
                                alldict = {i:fill_restdata_abs(np.frombuffer(ov_vars, dt, out['date_time'].shape[0]), rest_data, ori, recordindex-tslice.start, addedvar, absidx)} #ov_vars
                                write_dict_h5(fd, alldict, 'observations_table', {i: compfilter }, [i])  #{ 'compression': 32015, 'compression_opts':(3,) }
                            #ptime(f'{wigos} {i}  written ', tt)
                                
            
                        ptime(f'{wigos} obstable written ', tt)
                    #if True:
            #        with eua.CDMDataset(fn) as data:
                        #tt = time.time()
                        for i in fbkeys + ['fg_depar@offline']:
            
                            if i in reduced_fbkeys:
                                continue
                            else: 
                                rest_data = data.era5fb[i][tslice]
                                #ptime(f'{wigos} feedback written ', tt)
                                #if rest_data.ndim == 1:
                                     #print(np.max(rest_data[:]))
                                #if i in ['expver', 'source@hdr', 'source_id', 'statid@hdr']:
                                if rest_data.ndim > 1: #shape[1]i in ['expver', 'source@hdr', 'source_id', 'statid@hdr']:
                                    if len(rest_data[0])>32:
                                        print ('WARNING, truncated', i, 'to 32', data.filename)
                                        rest_data = rest_data[:, :32]
                                    #ov = np.empty((addedvar.shape[0],len(rest_data[0])), dtype=rest_data.dtype)
                                    #ov_vars.fill('n')
                                    miss_val = 'n'
                                    dt = rest_data.dtype
                                else:
                                    if rest_data.dtype in (np.float32, np.float64):
                                        dt = np.float32
                                        miss_val = np.nan
                                    else:
                                        dt = np.int32
                                        miss_val = -2147483648
                                    pass
                                    #if rest_data.dtype in (np.float64, np.float32):
                                        #ov_vars = np.empty(addedvar.shape[0], dtype=np.float32)
                                        #miss_val = np.nan
                                    #elif  rest_data.dtype == np.int64:
                                        #ov_vars = np.empty(addedvar.shape[0], dtype=np.int32)
                                        #miss_val = -2147483648
                                    #elif  rest_data.dtype == np.int32:
                                        #ov_vars = np.empty(addedvar.shape[0], dtype=np.int32)
                                        #miss_val = -2147483648
                                    #else:
                                        #raise ValueError('unknown dtype', rest_data.dtype)

                                if i in ('an_depar@surfbody_feedback','codetype@hdr', 'collection_identifier@conv','datum_anflag@body', 'datum_event1@body',
                                         'datum_sfc_event@surfbody_feedback', 'datum_status@body', 'datum_status@surfbody_feedback', 'eda_spread@errstat',
                                         'entryno@body', 'fg_depar@surfbody_feedback', 'fg_error@errstat', 'final_obs_error@errstat', 'obs_error@errstat',
                                         'obstype@hdr', 'obsvalue@body', 'orography@modsurf', 'qc_pge@body', 'report_event1@hdr', 'report_status@hdr',
                                         'seqno@hdr', 'station_type@conv', 'subtype@hdr', 'unique_identifier@conv','varbc_ix@body', 'varno@body',
                                         'vertco_reference_1@body','vertco_type@body', 'reportype', 'timeseries_index@conv', 'datum_rdbflag@body'): 
                                    ov = fill_restdata_addedvar_abs(np.frombuffer(ov_vars, dt, out['date_time'].shape[0]), rest_data, addedvar, absidx, miss_val)
                                elif rest_data.ndim == 1:
                                    ov = fill_restdata_abs(np.frombuffer(ov_vars, dt, out['date_time'].shape[0]), rest_data, ori, recordindex-tslice.start,addedvar, absidx)
                                else:
                                    sh = (out['date_time'].shape[0], rest_data.shape[1])
                                    ov = fill_restdata_abs(np.frombuffer(ov_vars, dt, sh[0]*sh[1]).reshape(sh), rest_data, ori, recordindex-tslice.start,addedvar, absidx)
                                    
                
                                if i == 'index':
                                    pass
                                elif i in ['expver', 'source@hdr', 'source_id', 'statid@hdr']:
                                    alldict = {i:ov.view('S1')}
                                    write_dict_h5(fd, alldict, 'era5fb', {i: compfilter  }, [i])
                                else:
                                    alldict = {i:ov}
                                    write_dict_h5(fd, alldict, 'era5fb', {i: compfilter  }, [i])
                                #ptime(f'{wigos} {i}  written ', tt)
                                
                        ptime(f'{wigos} feedback written ', tt)
                            #print(i, time.time()-tt)
                    #
                    # writing the recordindices and recordtimestamp.
                    #
                #print('ttx', fn, time.time()-ttx)
                recordindices=vridx
                cn=list(eua.cdm_codes.keys())
                cv =list(eua.cdm_codes.values())

                for i in range(len(recordindices)):
                    testvar = {str(allvars[i]):recordindices[i]}
                    try:            
                        idx = cv.index(allvars[i])
                    except:
                        cv.append(allvars[i])
                        if allvars[i] == 139:
                            cn.append('eastward_wind_speed')
                        elif allvars[i] == 140:
                            cn.append('northward_wind_speed')
                        else:
                            cn.append('unknown')
                        idx = cv.index(allvars[i])
                    write_dict_h5(fd, testvar, 'recordindices', {str(allvars[i]): { 'compression': None , 'compression_opts': None } }, [str(allvars[i])],
                                  attrs={str(allvars[i]):{'content': 'indices of first rows of data records for given variable in observations_table and era5fb tables','variable_name': cn[idx]}}) 
            
                write_dict_h5(fd, {'recordtimestamp':recordtimestamps}, 'recordindices', {'recordtimestamp': { 'compression': None , 'compression_opts': None} }, ['recordtimestamp'],
                              attrs={'recordtimestamp':{'units': fd['observations_table']['date_time'].attrs['units']}})
                with h5py.File(fn) as f:
                    for g in f.keys():
                        if g not in ('recordindex', 'recordtimestamp'):
                            
                            copy_attrs(f, fd, g)        
        
            fd.close()

            ptime('elapsed writing '+targetfile+':',tt, debug=True)
            try:
                os.mkdir(wpathy+'/log')
            except:
                pass
            f= open(wpathy+'/log/'+fn.split('/')[-1]+".txt","w")
            f.write("done") 
            f.close()
        except MemoryError as e:
            print('failed writing '+targetfile, e)
            raise ValueError(fn)
    except MemoryError as e:
        print('ERROR: failed reading '+targetfile, e)
        raise ValueError(fn)
        
    return targetfile

ray_convert_missing = ray.remote(convert_missing)

def rmsplot(files, do_plot=True, adjust=True):
    
    tt = time.time()
    readicts = []
    year = os.path.dirname(files[0]).split('/')[-1]
    pkname = os.path.basename(files[0]).split('_CEUAS')[0] + '.pkl'
    fno = wpath+'/plots/'+os.path.basename(pkname).split('.')[0].split('_')[0]
    if os.path.exists(fno+'.txtx'):
        print(fno+' already plotted')
        return
    
    print('rmsplot', files)
    for file in files:
        year = os.path.dirname(file).split('/')[-1]
        ss = 'CEUAS'
        fs = file.split('CEUAS')
        pkname = fs[0] + year+'.pkl'
        
        print(pkname)
        try:
            
            with open(pkname,'rb') as f:
                readicts.append(pickle.load(f))
        except Exception as e:
            print (pkname + ' not found', e)
    if adjust:
        lpath = files[0].split('/')
        lpath[-2] = 'long'
        lfile = '/'.join(lpath)
        #with h5py.File(lfile) as f:
            #ladj = f['advanced_homogenisation']['wind_bias_estimate'][:]
            #lobstype = f['observations_table']['observed_variable'][:]
            #ldate_time = f['observations_table']['date_time'][:]
            #lz_coordinate = f['observations_table']['z_coordinate'][:]
            #lobs = f['observations_table']['observation_value'][:]
    if len(readicts) == 0:
        print(pkname+' no feedback information found')
        return
    
    
    readict = {'obs' : {}}
    plevs = np.array((10000., 30000,50000,70000,85000, 92500 ))
    plevs = np.array((10000., 92500 ))
    print(readicts[0].keys())
    #print(time.time()-tt)

    mask = []
    for i in readicts:
        mask.append(np.isin(i['obs']['z_coordinate'], plevs))
    
    for k in readicts[0]['obs'].keys():
        readict['obs'][k] = np.concatenate([readicts[i]['obs'][k][mask[i]] for i in range(len(readicts))])
    
    #readict[lobs] = {}
    #for k in readicts[0]['obs'].keys():
        #readict['lobs']['l' + k] = readict['obs'][k][:] - 100000.
        #for p in plevs:
            #readict['lobs']['l' + k]
    #print(time.time()-tt)
    for p in readicts[0]['obstypes'].keys():

        for k,v in readicts[0].items():
            if k in ('obs', 'tslice'):
                continue
            if 'obs' in k:
                readict[k] = readicts[0][k]
                continue
            if k =='era5fb':
                i = 0
            try:
                
                for ftype in v['ftype']:  
                    if 'fc' in ftype:
                        dtype='fc'
                    else:
                        dtype='an'
                    if k not in readict.keys():
                        readict[k] = {}
                    if dtype not in readict[k].keys():
                        readict[k][dtype] = {}
                    readict[k][dtype]['refvalues'] = np.concatenate([readicts[i][k][dtype]['refvalues'][mask[i]]  for i in range(len(readicts))])
            except:
                
                pass
                #print(time.time()-tt)


    try:
        os.mkdir(wpath+'/plots')
    except:
        pass
    tt1 = time.time()-tt
    
#    res = readictplot(readict, ('mean', 'std','rms'), plevs, fno, marker='', do_plot=do_plot)
    res = readictplot(readict, ('rms', ), plevs, fno, marker='', do_plot=do_plot)
    if do_plot:
        res = readictplot(readict, ('mean', 'std','rms'), plevs, fno, marker='', minvals=2)
    for k in 'lat', 'lon', 'statid':
        try:
            
            res[k] = readicts[-1][k]
        except:
            print(res['fn'], 'no lat', readicts[-1].keys())
            res[k] = np.nan
    print(fno, tt1, time.time()-tt, 'secs')
    with open(fno+'.txt', 'w') as f:
        f.write('plots done')

    return res

ray_rmsplot = ray.remote(rmsplot)

def adjplot(fn, pars, plevs, pdict, do_plot=True):
    
    tt = time.time()
            #137: {'ap': 'humidity_bias_estimate',},}    
    sid = os.path.basename(fn).split('_CEUAS')[0]
    x = 0
    with h5py.File(fn) as f:
        lat = f['observations_table']['latitude'][-1]
        lon = f['observations_table']['longitude'][-1]
        nrdict = {}
        ydict = {}
        dsdict = {}
        dmdict = {}
        for par in pars:
            try:
                
                sl = slice(f['recordindices'][str(par)][0], f['recordindices'][str(par)][-1])
                if sl.stop - sl.start < 30:
                    continue
                nrdict[par] = {}
                nrdict[par]['fobs'] = f['observations_table']['observation_value'][sl]
                #nrdict[par][plev]['fbc'] = f['era5fb']['biascorr@body'][sl]
                nrdict[par]['ffgdep'] = f['era5fb']['fg_depar@body'][sl] + f['era5fb']['biascorr@body'][sl]
                nrdict[par]['fnfgdep'] = f['era5fb']['fg_depar@offline'][sl]
                nrdict[par]['fdatetime'] = f['observations_table']['date_time'][sl]
                nrdict[par]['fah'] = f['advanced_homogenisation'][pdict[par]['ap']][sl]
                nrdict[par]['fz'] = f['observations_table']['z_coordinate'][sl]
                ydict[par] = {}
                dsdict[par] = {}
                dmdict[par] = {}
                for plev in plevs:

                    pi = np.where(nrdict[par]['fz']==plev)

                    nrdict[par][plev] = {}
                    #nrdict[par][plev]['obstype'] = f['observations_table']['observed_variable'][sl][pi]
                    nrdict[par][plev]['obs'] = nrdict[par]['fobs'][pi]
                    nrdict[par][plev]['fgdep'] = nrdict[par]['ffgdep'][pi]
                    nrdict[par][plev]['nfgdep'] = nrdict[par]['fnfgdep'][pi]
                    nrdict[par][plev]['datetime'] = nrdict[par]['fdatetime'][pi]
                    nrdict[par][plev]['ah'] = nrdict[par]['fah'][pi]
                    ydict[par][plev] = {}
                    dsdict[par][plev] = {}
                    dmdict[par][plev] = {}

                    
                    if np.sum(~np.isnan(nrdict[par][plev]['nfgdep'])) == 0:
                        nrdict[par][plev]['nfgdep'] = nrdict[par][plev]['fgdep'][:]
                        
                    if par == 39:
                        idx = np.searchsorted(nrdict[138][plev]['datetime'], nrdict[par][plev]['datetime'])
                        idx[idx==nrdict[138][plev]['datetime'].shape[0]] = nrdict[138][plev]['datetime'].shape[0] - 1
                        obshom = nrdict[138][plev]['obs'][idx] - nrdict[138][plev]['ah'][idx]
                        
                        idx = np.searchsorted(nrdict[126][plev]['datetime'], nrdict[par][plev]['datetime'])
                        idx[idx==nrdict[126][plev]['datetime'].shape[0]] = nrdict[126][plev]['datetime'].shape[0] - 1
                        vpdata = obshom * Sonntag(nrdict[126][plev]['obs'][idx])

                        adj = nrdict[par][plev]['obs'] - vap2sh(vpdata, plev)

                        nrdict[par][plev]['adj'] = nrdict[par][plev]['nfgdep'] - adj
                    else:
                        nrdict[par][plev]['adj'] = nrdict[par][plev]['nfgdep'] - nrdict[par]['fah'][pi]
                    
                    q = np.nanquantile(np.concatenate((nrdict[par][plev]['nfgdep'], nrdict[par][plev]['fgdep'])), (0.001, 0.999))
                    #print(par, plev, q)
                    bi = (nrdict[par][plev]['nfgdep'] < q[0]) | (nrdict[par][plev]['nfgdep'] > q[1])
                    for quant in 'obs', 'fgdep', 'nfgdep', 'adj':
                        nrdict[par][plev][quant][bi] = np.nan
                    
                    
                    if do_plot : plt.subplot(2, 1, 1)
                    #for quant in 'fgdep', 'nfgdep', 'adj':
                        #if do_plot :plt.plot(nrdict[par][plev]['datetime']/86400/365.25+1900, nrdict[par][plev][quant]*pdict[par]['fak'], label=f'{quant}: {np.nanstd(nrdict[par][plev][quant]*pdict[par]['fak']):5.3f}')
                    if do_plot: plt.title(f'ERA5 departures {os.path.basename(fn).split("_CEUAS")[0]}, {pdict[par]["short"]}, {np.int32(plev/100)} hPa')
                    #plt.legend()
                    #plt.ylabel(pdict[par]['units'])
                    #plt.subplot(2, 1, 1)
                    
                    for quant in 'fgdep', 'nfgdep', 'adj':
                        ylist = np.arange(124)
                        ydict[par][plev][quant] = ylist[:]
                        dlist = []
                        for iy in ylist:
                            idx = np.searchsorted(nrdict[par][plev]['datetime'], (iy*365.25*86400,(iy+1)*365.25*86400))
                            if idx[1] -idx[0] > 30:
                                dlist.append(np.nanstd(nrdict[par][plev][quant][idx[0]:idx[1]]))
                                x = 1
                            else:
                                dlist.append(np.nan)
                        dsdict[par][plev][quant] = np.array(dlist)*pdict[par]['fak']
                                
                        if do_plot: plt.plot(ylist+1900.5,dsdict[par][plev][quant], label=f'{quant}: {np.nanmean(dsdict[par][plev][quant]):5.3f}')
                    
                    if do_plot:
                        plt.legend()
                        plt.ylabel('Std '+ pdict[par]['units'])
                        plt.subplot(2, 1, 2)
                    
                    for quant in 'fgdep', 'nfgdep', 'adj':
                        dlist = []
                        for iy in ylist:
                            idx = np.searchsorted(nrdict[par][plev]['datetime'], (iy*365.25*86400,(iy+1)*365.25*86400))
                            if idx[1] -idx[0] > 30:
                                dlist.append(np.nanmean(nrdict[par][plev][quant][idx[0]:idx[1]]))
                                x = 1
                            else:
                                dlist.append(np.nan)
                        dmdict[par][plev][quant] = np.array(dlist) *pdict[par]['fak']
                                
                        if do_plot: plt.plot(ylist+1900.5,dmdict[par][plev][quant], label=f'{quant}: {np.nanstd(dmdict[par][plev][quant]):5.3f}')
                    
                    mask = ylist != ylist
                    for quant in 'fgdep', 'nfgdep', 'adj':
                        if np.any(~np.isnan(dmdict[par][plev][quant])):
                            xx = 0
                        mask = mask | ( np.isnan(dmdict[par][plev][quant]))
                    for quant in 'fgdep', 'nfgdep', 'adj':
                        dmdict[par][plev][quant][mask] = np.nan
                        dsdict[par][plev][quant][mask] = np.nan
                    if do_plot:
                        plt.legend()
                        plt.ylabel('Mean '+pdict[par]['units'])
                        plt.tight_layout()
                        pout = '/'.join(fn.split('/')[:-2])+f'/plots/{sid}_{pdict[par]["short"]}_{np.int32(plev/100)}.png'
                        #print(os.path.basename(pout))
                        plt.savefig(pout)
                        plt.close()
            except:
                pass #print(sid+' failed')
    if x == 1:
        print(sid, x, time.time()-tt)
        return ydict, dsdict, dmdict, sid, lat, lon
    else:
        return None

ray_adjplot = ray.remote(adjplot)

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
    #opath=wpath
    #wlpath=wpath+'log/'
    #for p in wpath,wlpath, wpath + '/long', wpath + '/rea':      
        #try:
            #os.mkdir(p)
        #except:
            #pass


    #import numba_drift
    
    df = rsd.numba_drift.load_test_data()
    df['u'] = - np.abs(df.windSpeed) * np.sin(np.radians(df.windDirection))
    df['v'] = - np.abs(df.windSpeed) * np.cos(np.radians(df.windDirection))
    latd,lond,timed = rsd.numba_drift.trajectory(df.lat[0], df.lon[0], np.array(df.u), np.array(df.v), np.array(df.pressure), np.array(df.airTemperature))

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
    #ofiles = glob.glob('/mnt/scratch/scratch/federico/MERGED_FEB2023/*v1.nc')[:] #+ \
    dsave = os.getcwd()
    #os.chdir('/mnt/scratch/scratch/federico/REMERGE_GIUB_12072023/')
    #remerged = glob.glob('*v1.nc')[:]
    #os.chdir('/mnt/scratch/scratch/federico/MERGED_14082023_npsound_coordinate/')
    remerged = glob.glob('*v1.nc')[:]
    os.chdir(dsave)
    files = []
    #for o in ofiles:
        #if os.path.basename(o) not in remerged:
            #files.append(o)
        #else:
            #continue
    
    #files += glob.glob('/mnt/scratch/scratch/federico/REMERGE_GIUB_12072023/*v1.nc')[:] #+ \
    #files += glob.glob('/mnt/scratch/scratch/federico/MERGED_14082023_npsound_coordinate/*v1.nc')[:] #+ \
    
    #files = glob.glob('/mnt/scratch/scratch/federico/MERGED_YEARLY_06DEC_FULLDATABASE/*/*v3.nc')[:]
    #files = glob.glob('/mnt/scratch/scratch/federico/MERGED_YEARLY_14FEB2024/*/*v3.nc')[:]
    #files = glob.glob('/mnt/scratch/scratch/federico/MERGED_YEARLY_14FEB2024*/*/*v3.nc')[:]
    #files = glob.glob('/mnt/users/scratch/leo/scratch/raymerge/MERGED_YEARLY_27FEB2024/*11035*/*v3.nc')[:]
    #files = glob.glob('/mnt/scratch/scratch/federico/MERGED_YEARLY_0MAR2024_REGULAR/*/*v3.nc')[:]
    try:
        
        with open('slist.pkl', 'rb') as f:
            files = pickle.load(f)
    except:
        
        #files = glob.glob('/mnt/users/scratch/leo/scratch/FH/MERGED_YEARLY_0MAR2024_REGULAR/*58457/*v3.nc')[:]
        #files = glob.glob('/mnt/users/scratch/leo/scratch/UH/MERGED_YEARLY_20NOV2024_REGULAR/*/*v3.nc')[:]
        #files += glob.glob('/mnt/users/scratch/leo/scratch/UH/MERGED_YEARLY_20NOV2024_REGULAR_mobile/*/*v3.nc')[:]
        #files += glob.glob('/mnt/users/scratch/leo/scratch/UH/MERGED_YEARLY_20NOV2024_REGULAR_orphan/*/*v3.nc')[:]
        #files = glob.glob('/mnt/users/scratch/leo/scratch/UH/MERGED_YEARLY_25JAN2025_REGULAR/*/*v3.nc')[:]
        #files += glob.glob('/mnt/users/scratch/leo/scratch/UH/MERGED_YEARLY_25JAN2025_REGULAR_mobile/*/*v3.nc')[:]
        #files += glob.glob('/mnt/users/scratch/leo/scratch/UH/MERGED_YEARLY_25JAN2025_REGULAR_orphan/*/*v3.nc')[:]
        files = glob.glob('/mnt/scratch/scratch/leo/scratch/UH/MERGED_YEARLY_19FEB25/*/*v3.nc')[:]
        files += glob.glob('/mnt/scratch/scratch/leo/scratch/UH/MERGED_YEARLY_19FEB25_mobile/*/*v3.nc')[:]
        files += glob.glob('/mnt/scratch/scratch/leo/scratch/UH/MERGED_YEARLY_19FEB25_orphan/*/*v3.nc')[:]
        #files = glob.glob('/mnt/users/scratch/leo/scratch/FH/MERGED_YEARLY_0MAR2024_REGULAR/*94975/*v3.nc')[:]

        #dirs = glob.glob('/mnt/users/scratch/leo/scratch/FH/MERGED_YEARLY_0MAY2024_REGULAR/*')[:]
        #for d in dirs:
            #dd = d.split('/')[-1]
            #rlist = glob.glob('/mnt/users/scratch/leo/scratch/converted_v20/*/*'+dd+'*')
            #for fl in rlist:
                #os.remove(fl)
                
        #files += glob.glob('/mnt/users/scratch/leo/scratch/FH/MERGED_YEARLY_0MAR2024_REGULAR_orphan/*/*v3.nc')[:]
        #files += glob.glob('/mnt/users/scratch/leo/scratch/FH/MERGED_YEARLY_0MAR2024_REGULAR_mobile/*/*v3.nc')[:]
        
        with open('slist.pkl', 'wb') as f:
            d = pickle.dump(files, f)
    
    #files = glob.glob('/mnt/users/scratch/leo/scratch/FH/MERGED_YEARLY_0MAR2024_REGULAR_orphan/*/*v3.nc')[:]
#    files += glob.glob('/mnt/users/scratch/leo/scratch/premerge/*/*merged*v3.nc')[:] #+ \
            
    bfiles = [os.path.basename(f) for f in files]

    
    inter = False
    if inter:

        files = glob.glob('/mnt/scratch/scratch/federico/INTERCOMPARISON/MAURITIUS/*.nc')
        fndict = {}
        newfiles = []
        i = 0
        for fn in files:
            if('vaisala') in fn and 'digi' not in fn:
                i += 1
                key = fn.split('/')[-1].split('_')[-1][:-3]
                fndict[key] = f'0-20100-0-019{i:0>2d}'
                newfiles.append('/tmp/'+fndict[key]+'_CEUAS_merged_v1.nc')
                os.system('cp '+fn+' '+newfiles[-1])
                #with h5py.File(newfiles[-1], 'r+') as f:
                    #f['observations_table/latitude'][:] = f['observations_table/latitude'][0]
                    #f['observations_table/longitude'][:] = f['observations_table/longitude'][0]
    
        df =pd.read_csv('/mnt/scratch/scratch/federico/databases_service2/MAURITIUS/vaisala_ascents.csv',sep=',')
        print(df.columns)
        h = df.height.values
        pc = df.pressure.values
        t = df.temperature.values
        dt, ri= np.unique(df.date_time.values, return_index=True)
        del df
        with h5py.File(newfiles[-1], 'r+') as g:
            sh = g['observations_table']['date_time'].shape
            ts,ris = np.unique(g['observations_table']['date_time'][:], return_index=True)
            hx = []
            pcx = []
            tx = []
            wdx = []
            wsx = []
            rhx = []
            for ii in range(ris.shape[0]):
                if ii < ris.shape[0] - 1:
                    
                    dslice = slice(ris[ii], ris[ii+1])
                    hslice = slice(ri[ii], ri[ii+1])
                else:
                    dslice = slice(ris[ii], sh[0])
                    hslice = slice(ri[ii], h.shape[0])
                ov = g['observations_table']['observed_variable'][dslice]
                vals = g['observations_table']['observation_value'][dslice]
                p = g['observations_table']['z_coordinate'][dslice][ov == 126]
                pcy,rip = np.unique(pc[hslice]* 100,return_index=True)
                
                hy = h[hslice][rip]
                mask = hy!=-999.
                hy = hy[mask]
                pcy = pcy[mask]
                rip = rip[mask]
                hx.append(hy)
                pcx.append(pcy)
                tx.append(vals[ov == 126][mask])
                wdx.append(vals[ov == 106][mask])
                wsx.append(vals[ov == 107][mask])
                rhx.append(vals[ov == 138][mask])
            
            g.attrs['oname'] = fn
            
                    
            for fn in files:
                
                if 'digi' in fn:
                    i += 1
                    print(i, fn)
                    key = fn.split('/')[-1].split('_')[-1][:-3]
                    fndict[key] = f'0-20100-0-019{i:0>2d}'
                    newfiles.append('/tmp/'+fndict[key]+'_CEUAS_merged_v1.nc')
                    os.system('cp '+fn+' '+newfiles[-1])
                    with h5py.File(newfiles[-1], 'r+') as f:
                        
                        sh = f['observations_table']['date_time'].shape
                        ts2,rid = np.unique(f['observations_table']['date_time'][:], return_index=True)
                        zs = []
                        ps = []
                        rhs = []
                        temps = []
                        wds = []
                        wss = []
                        oval = []
                        ov = []
                        rids = []
                        dts = []
                        
                        for ii in range(len(ts2)):
                            idx = np.searchsorted(ts, ts2[ii])
                            if idx < ts.shape[0] and ts[idx] == ts2[ii]:
                                if ii < ts2.shape[0] - 1:
                                    fslice = slice(rid[ii], rid[ii+1])
                                else:
                                    fslice = slice(rid[ii], sh[0])
                            else:
                                continue
                            fov = f['observations_table']['observed_variable'][fslice]
                            pdict = {126: tx,138: rhx}
                            rlen = 0
                            for par in 126, 138:
                                
                                fz = f['observations_table']['z_coordinate'][fslice][fov == par]
                                zidx,uidx = np.unique(np.searchsorted(fz, hx[idx][::-1]), return_index=True)
                                mask = zidx<fz.shape[0]
                                zidx = zidx[mask]
                                uidx = uidx[mask]
                                oval.append(f['observations_table']['observation_value'][fslice][fov == par][zidx][::-1])                        
                                ov.append(np.repeat(par, uidx.shape[0]))                        
                                ps.append(pcx[idx][::-1][uidx][::-1])
                                rlen += len(ps[-1])
    
                            for v in 117, 106, 107:
                                ov.append(np.repeat(v, uidx.shape[0]))
                            m = 0
                            for v in hx[idx][::-1], wdx[idx][::-1], wsx[idx][::-1]:
                                if m == 0: 
                                    oval.append(v[uidx][::-1]) *9.80665
                                else:
                                    oval.append(v[uidx][::-1])
                                    
                                ps.append(pcx[idx][::-1][uidx][::-1])
                                rlen += len(ps[-1])
    
                            rids .append(np.repeat([f['observations_table']['report_id'][fslice.start]], rlen, axis=0))
                            dts .append(np.repeat(f['observations_table']['date_time'][fslice.start], rlen))
                                
                        #print(f['observations_table']['z_coordinate'].shape)
                        sid = f['observations_table']['sensor_id'][0]
                        nlist = 'sensor_id','z_coordinate_type','z_coordinate', 'observation_value', 'observed_variable', 'report_id', 'observation_id', 'date_time'
                        for d in  nlist:                   
                            del f['observations_table'][d]
                            
                        for d in f['observations_table'].keys():
                            if np.any((f['observations_table'][d][0] != f['observations_table'][d][-1])):
                                print(d, f['observations_table'][d][0],f['observations_table'][d][-1] )
    
                        f['observations_table'].create_dataset('z_coordinate', data=np.concatenate(ps))
                        index = f['observations_table']['z_coordinate'].shape[0]
                        f['observations_table'].create_dataset('z_coordinate_type', data=np.repeat(1,index ))
                        f['observations_table'].create_dataset('observed_variable', data=np.concatenate(ov))
                        f['observations_table'].create_dataset('observation_value', data=np.concatenate(oval))
                        f['observations_table'].create_dataset('sensor_id', data=np.repeat([sid], index, axis=0))
                        f['observations_table'].create_dataset('report_id', data=np.concatenate(rids))
                        f['observations_table'].create_dataset('date_time', data=np.concatenate(dts))
                        obsid = np.arange(index)
                        sobsid = np.array([f'{i:0{10}}' for i in range(index)]).astype('S10')
                        f['observations_table'].create_dataset('observation_id', data=sobsid.view('S1').reshape(index, 10))
                        
                        for d in f['observations_table'].keys():
                            if d not in nlist and d != 'index' and 'string' not in d:
                                v = f['observations_table'][d][0]
                                del f['observations_table'][d]
                                if v.ndim > 0:           
                                    f['observations_table'][d] = np.repeat([v], index, axis=0)
                                else:
                                    f['observations_table'][d] = np.repeat(v, index)
                        
                        dim_attach(f,'observations_table')
                            
                        #print(f['observations_table']['z_coordinate'].shape)
                        f.attrs['oname'] = fn
                    x = 0
                                
                                #fig = plt.figure()
                                #plt.subplot(1, 4, 1)
                                #plt.semilogy( oval[-4], ps[-4]/100.)
                                #plt.semilogy( pdict[par][idx][::-1][uidx], ps[-4]/100.)
                                #plt.ylim(1000., 5)
                                #plt.subplot(1, 4, 2)
                                #plt.semilogy( pdict[par][idx][::-1][uidx]-oval[-4], ps[-4]/100.)
                                #plt.ylim(1000., 5)
                                
                                #plt.subplot(1, 4, 3)
                                #plt.semilogy( fz[zidx], ps[-4]/100.)
                                #plt.semilogy( hx[idx][::-1][uidx], ps[-4]/100.)
                                #plt.ylim(1000., 5)
                                #plt.subplot(1, 4, 4)
                                #plt.semilogy( hx[idx][::-1][uidx]-fz[zidx], ps[-4]/100.)
                                #plt.ylim(1000., 5)
                                #x = 0
                                #plt.close(fig)
                    
                    #f['observations_table/latitude'][:] = f['observations_table/latitude'][0]
                    #f['observations_table/longitude'][:] = f['observations_table/longitude'][0]
        
        
    
        files = glob.glob('/mnt/scratch/scratch/federico/INTERCOMPARISON/YANGJIANG/*.nc')
        i = 0
        for fn in files:
            i += 1
            key = fn.split('/')[-1].split('_')[-1][:-3]
            fndict[key] = f'0-20100-0-020{i:0>2d}'
            newfiles.append('/tmp/'+fndict[key]+'_CEUAS_merged_v1.nc')
            os.system('cp '+fn+' '+newfiles[-1])
        files =newfiles    
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

    #ray.init( num_cpus=64, _temp_dir=os.path.expandvars('/srvfs/fastscratch/scratch/leo/ray'))
#    ray.init(address="localhost:6379")
    
    
    refs = load_20CRoffset(readict)
    
    inventory_list = glob.glob('../../meta/inventory_comparison_2/code/station_configuration/CUON_station_configuration_extended.csv')
    pd_list, station_elevation = read_station_height(files_to_convert[0],100, inventory_list)
    print('read inventory', time.time()-tt)

    readicts = []

    
    refs_ref = ray.put(refs)
    pd_list_ref = ray.put(pd_list)
    del refs
    del pd_list

    finput = 'yearly'
    fdict = {}
    dat = datetime.now()
    ref = datetime(1970, 1, 1)
    age = (dat - ref).total_seconds()
    for fn in files:
        year = int(fn.split('_CEUAS')[0][-4:])
        if age > os.path.getmtime(fn) + 300:
            try:           
                fdict[year].append(fn)
            except:
                fdict[year] = [fn]
        
    if finput == 'yearly':
        wpath = '/mnt/users/scratch/leo/scratch/converted_v29/'
        #wpath = '/run/user/73643/'
        #wpath = '/mnt/ssdraid/scratch/leo/converted_v18/'
        rts = True
        drift = True
        futures = []
        for iyear in range(1972, 1900, -1):
            if iyear not in fdict:
                continue
            for file in tqdm(fdict[iyear]): #glob.glob('/mnt/scratch/scratch/federico/MERGED_YEARLY_06DEC_FULLDATABASE/*/*_'+str(iyear)+'_*v3.nc'):
                #if (dat - ref).total_seconds() > os.path.getmtime(file) + 300:
                    
                
                    #if  'np_05sound' in file or 'np_04sound' in file in file: #or '0-07602_CEUAS_merged_v1' in file:
                    #fn = '/mnt/users/scratch/leo/scratch/test/' + '/'.join(file.split('/')[-2:])
                    fn = file
                    wigosid = file.split('/')[-2]
                    try:
                        
                        #with h5py.File(f'{wpath}/{iyear}/{wigosid}_CEUAS_merged_v3.nc') as f:
                            #lons = f['observations_table']['longitude'][:]
                            #lond = f['observations_table']['lond'][:]
                        #if np.any((lons>=0)&((lons+lond)<0)):
                                
    #                    blist = ['0-20999-0-ASUK02', '0-20999-0-A', '0-20999-0-DBLK', '0-20999-0-WDK38HS', '0-20999-0-SMLQ']
                        if  True or '-0-11952' in file: #True or '0-20000-0-89009' in file: #or '0-20000-0-07354' in file or '0-20001-0-11035' in file: #or '0-20000-0-70219' in file or '0-20001-0-10393' in file or  '0-20000-0-03882' in file or '0-20000-0-07510' in file or '0-20001-0-11035' in file:
    #                    if True or wigosid in blist: #'47600' in wigosid or '94975' in wigosid or '94120' in wigosid or '11035' in wigosid  #True or '0-20001-0-11035' in wigosid: #'0-20000-0-94461' in wigosid: #wigos in ['0-20000-0-72357']: #'0-20001-0-10393']: #'0-20000-0-72357','0-20001-0-10393']: #'0-20000-0-72357', '0-20001-0-11035']: #, '0-20001-0-11035', '0-20001-0-10393']: #, '0-20000-0-72357' , '0-20000-0-70219' ,'0-20000-0-91413' , '0-20001-0-10393'] :
                            #if (iyear == 2023 and 'premerge' in file) or iyear < 2023:
                            convert_missing(refs_ref, pd_list_ref, None, wpath, iyear, file, record_timestamp=rts, drift=drift)
                            #
                            #futures.append(ray_convert_missing.remote(refs_ref, pd_list_ref,None,wpath, iyear, fn,record_timestamp=rts,drift=drift ))
                            #print('ZERO', iyear, fn)
                    except MemoryError:
                        pass
            
            #futures = futures + [ ray_convert_missing.remote(refs_ref, None, wpath, iyear, file )  for file in files_to_convert]
        obj_result_list = ray.get(futures)
    else:
        futures = []
        for year in range(-2005, 2004, -1):
    
            for file in files_to_convert[5:]:
                #if  'np_05sound' in file or 'np_04sound' in file in file: #or '0-07602_CEUAS_merged_v1' in file:
                if True or '11035' in file:
                    #convert_missing(refs, pd_list, None, wpath, year, file)
                    futures.append(ray_convert_missing.remote(refs_ref, pd_list_ref,None, wpath, year, file ))
            
            #futures = futures + [ ray_convert_missing.remote(refs_ref, None, wpath, year, file )  for file in files_to_convert]
        obj_result_list = ray.get(futures)
        
    
    #for fn in [f'/mnt/users/scratch/leo/scratch/converted_v24/{year}/0-20001-0-11035_CEUAS_merged_v3.nc' for year in range(2000, 2007)]:
        #with h5py.File(fn) as f:
            #mask = np.where((f['observations_table']['observed_variable'][:]==34)&(f['observations_table']['z_coordinate'][:]==10000))[0]
            #print(np.nanstd(f['era5fb']['fg_depar@offline'][:][mask]))

    print('yearly files finished')
    #exit(0)
    if False:
        flist = glob.glob(wpath+'/long/*v3.nc')
        plevs = np.array([3000., 5000., 10000., 30000., 50000., 70000., 85000., 92500.])
        pars = [140, 126, 138, 39, 139]
        nreadict = copy.copy(readict)
        fgstats =[]
        pdict ={140: {'ap': 'wind_bias_estimate','short': 'v','units': 'm/s','fak': 1,},
                139: {'ap': 'wind_bias_estimate','short': 'u','units': 'm/s','fak': 1,},
                126: {'ap': 'RISE_bias_estimate','short': 't','units': 'K','fak': 1, },
                138: {'ap': 'humidity_bias_estimate','short': 'rh','units': '','fak': 1, },    
                39: {'ap': 'humidity_bias_estimate','short': 'q','units': 'g/kg','fak': 1000}, }   
        
        for fn in flist[:]:
            futures.append(ray_adjplot.remote(fn, pars, plevs, pdict, True))
            #fgstats.append(adjplot(fn, pars, plevs, pdict))

        fgstats = ray.get(futures)
        fgstats = [f for f in fgstats if f is not None]
        
        #ylist = fgstats[0][0][pars[0]][plevs[0]]['fgdep']
        #dslist = fgstats[0][1][pars[0]][plevs[0]]['fgdep']
        ylist = np.arange(1900, 2024) + 0.5
        sid = 'all'
        for par in pars: #pars[1:2]:
            for plev in plevs: #plevs[6:7]: #plevs:
                plt.subplot(2, 1, 1)
                for quant in  fgstats[0][2][pars[0]][plevs[0]].keys():
                    dmlist = []
                    dslist = []
                    for f in fgstats:
                        try:
                            
                            dmlist.append(f[2][par][plev][quant])
                            dslist.append(f[1][par][plev][quant]**2)
                            print(f[3], np.nanmax(dslist[-1]))
                        except:
                            pass
                    ds = np.sqrt(np.nanmean(np.array(dslist), axis=0))
                    #plt.plot(ylist+1900.5,dm, label=f'{quant}: {np.nanstd(dm):5.3f}')
                    plt.plot(ylist,ds, label=f'{quant}: {np.nanmean(ds):5.3f}')
                    
                plt.title(f'Global average ERA5 departures {pdict[par]["short"]}, {np.int32(plev/100)} hPa')
                plt.ylabel('Std '+pdict[par]['units'])
                plt.legend()
                plt.subplot(2, 1, 2)
                for quant in  fgstats[0][2][pars[0]][plevs[0]].keys():
                    dmlist = []
                    dslist = []
                    for f in fgstats:
                        try:
                            
                            dmlist.append(f[2][par][plev][quant])
                            print(f[3], np.nanmax(dslist[-1]))
                        except:
                            pass
                    dm = np.nanmean(np.array(dmlist), axis=0)
                    #plt.plot(ylist+1900.5,dm, label=f'{quant}: {np.nanstd(dm):5.3f}')
                    plt.plot(ylist,dm, label=f'{quant}: {np.nanstd(dm):5.3f}')
                plt.ylabel('Mean '+pdict[par]['units'])
                plt.legend()
                plt.tight_layout()
                pout = '/'.join(fn.split('/')[:-2])+f'/plots/{sid}_{pdict[par]["short"]}_{np.int32(plev/100)}.png'
                #print(os.path.basename(pout))
                plt.savefig(pout)
                plt.close()
        
        
                
    if False:
        try:
            
            with open('xalldata.pkl', 'rb') as f:
                dlist = pickle.load(f)
        except:
        
            pwd = os.getcwd()
            obj_result_list = glob.glob(wpath+'/[12]???/*v3.nc')
            os.chdir(pwd)
            shlist = np.unique([os.path.basename(ll) for ll in obj_result_list if ll is not None])
            futures = []
            vienna = False
            dlist = []
            for file in shlist:
                if  '11035' in file:
                    vienna = True
                sublist = sorted([s for s in obj_result_list if s is not None and file in s])
                #dlist.append(rmsplot(sublist, do_plot=True))
                #if True or vienna:
                futures.append(ray_rmsplot.remote(sublist, do_plot=True))
            dlist = ray.get(futures)
            with open('alldata.pkl', 'wb') as f:
                pickle.dump(dlist, f)
        
        exit(0)        
        dlist = [d for d in dlist if d is not None]

        times = np.arange(2024-1900, dtype=np.int32)
        ttt = time.time()
        cond = {'gl': {'min': -90, 'max': 90, 'era5fb': True,},
                'nh': {'min': 20, 'max': 80, 'era5fb': True,},
                'sh': {'min': -80, 'max': -20, 'era5fb': True,},
                'tr': {'min': -20, 'max': 20, 'era5fb': True,}
                }
        for ck, cv in cond.items():
            
            for stat in 'rms', :
                for pl in 925, 850, 700, 500, 300, 100:
                    units = {'t':'K','q':'kg/kg','u': 'deg','v': 'm/s','z': 'J/kg',}
                    for par in 't','q','u','v','z' :
                        
                        for k in 'era5','era5fb':
                            for ftype in 'fc', 'an':
                                
                            
                        
                                res = []
                                for d in dlist:
                                    try:
                                        res.append({'data': d[pl][par][stat][k][ftype], 'lat': d['lat'],'lon': d['lon'],'statid':d['statid']})
                                        #print(d['statid'], pl, par, stat, k, ftype)
                                        x = 0
                                    except Exception:
                                        if d is not None:
                                            res.append({'data': None, 'lat': d['lat'],'lon': d['lon'],'statid':d['statid']})
                                    
                                            
                                        #print('no', pl, par, stat, k, ftype)
                                        pass
                                rmss = np.zeros(times.shape[0])
                                n = np.zeros(times.shape[0])
                                maxs = np.zeros(times.shape[0])
                                maxi = np.zeros(times.shape[0], dtype=np.int32) - 1
                                si = 0
                                for r in res:
                                    if r['lat'] > cv['min'] and r['lat'] < cv['max'] and r['data']:
                                        idx = r['data'][0] - 1900
                                        for i in range(len(idx)):
                                            if r['data'][1][i] == r['data'][1][i]:
                                                if r['data'][1][i] > maxs[idx[i]]:
                                                    
                                                    maxs[idx[i]] = r['data'][1][i]
                                                    maxi[idx[i]] = si
                                    si += 1
                                
                                amax = np.argmax(maxs)
                                u =np.unique(maxi)
                                if k == 'era5' and ftype == 'fc':
                                    
                                    print('amax', ck, pl, par, stat, k, ftype, f'{np.max(maxs):5.3f}',res[maxi[amax]]['statid'], 1900+amax)
                                #for us in u:
                                    #print(res[us]['statid'])
                                si = 0
                                for r in res:
                                    if r['lat'] > cv['min'] and r['lat'] < cv['max']   and r['data']:
                                        
                                        fb = dlist[si][pl][par][stat]['era5fb'][ftype]
                                        
                                        if np.all(fb[0]==r['data'][0]):
                                            
                                            idx = r['data'][0] - 1900
                                            mask = (~np.isnan(r['data'][1])) &  (~np.isnan(fb[1])) &  (r['data'][1] < maxs[idx])
                                            n[idx[mask]] += 1
                                            rmss[idx[mask]] += r['data'][1][mask] * r['data'][1][mask]
                                    si += 1
                                mask = n > 5
                                rms = rmss + np.nan
                                rms[mask] = np.sqrt(rmss[mask]/n[mask])
                                plt.plot(1900+times, rms, label=k+','+ftype)
                        plt.title(f'{par}_{pl}_{stat}_{ck}')
                        plt.ylabel(units[par])
                        plt.legend()
                        plt.savefig(f'all{ck}_{par}_{pl}_{stat}.png')
                        plt.close()
    
        print(time.time()-ttt)
        
    pwd = os.getcwd()
    os.chdir(wpath)
    longlist = glob.glob('[12]???/*v3.nc') 
    os.chdir(pwd)
    shlist = []
    shlist = [os.path.basename(ll) for ll in longlist]
    fkeys = np.unique(shlist)
    
    futures =[]   
    for fks in fkeys:
        #if '0-94998' not in fks : #and '0-11035' not in fks: #and '0-70219' not in fks:
            #continue
        fkey = wpath+'/[12]???/' + fks # .split('_')[0] + '_????_' + '_'.join(fks.split('_')[2:])
                
        if  True or '0-20001-0-11035' in fkey:
            #h5concatenate(fkey)
            futures .append(ray_h5concatenate.remote(fkey))
        
    ray.get(futures)
    
    print('total:',time.time()-tt)
    print("Don't forget to run orphan script")
    exit()
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

