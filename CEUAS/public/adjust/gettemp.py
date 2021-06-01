#!/usr/bin/env python
import traceback
import sys,glob
import os.path

import numpy
import datetime
from datetime import date
import netCDF4
import time
from numba import njit
sys.path.append('./')
from utils import *
from anomaly import *
import matplotlib.pylab as plt
import scipy.stats
import f90nml
import xarray as xr
sys.path.append('../cds-backend/code/')
import cds_eua3 as eua
import cdsapi
import urllib3
import json
import h5py
import zipfile, os

# with open(os.path.expanduser('~leo/python/hug2/config/active.json')) as f:
#     active=json.load(f)
# ids=list(active.keys())
# lats=numpy.asarray([active[x][2] for x in active.keys()])   
# lons=numpy.asarray([active[x][3] for x in active.keys()])
# starts=numpy.asarray([active[x][0] for x in active.keys()])   
# stops=numpy.asarray([active[x][1] for x in active.keys()])
# l=0
# for i in range(lats.shape[0],lats.shape[0]):
#     idx=numpy.where(numpy.logical_and(numpy.abs(lats[i]-lats)<0.1,numpy.abs(lons[i]-lons)<0.1))[0]
#     if len(idx)>1:
#         fak=86400*365.25
#         print('duplicate {:s},{:s},{:4.0f},{:4.0f},{:4.0f},{:4.0f}'.format(ids[idx[0]],ids[idx[1]],
#                                                                            starts[idx[0]]/fak,starts[idx[1]]/fak,stops[idx[0]]/fak,stops[idx[1]]/fak))
#         try:
            
#             with h5py.File('/raid60/scratch/leo/scratch/converted_v5/'+ids[idx[0]]+'_CEUAS_merged_v1.nc','r') as f:
#                 with h5py.File('/raid60/scratch/leo/scratch/converted_v5/'+ids[idx[1]]+'_CEUAS_merged_v1.nc','r') as g:
#                     try:
#                         print(f['observations_table']['latitude'][0],f['observations_table']['longitude'][0],
#                               g['observations_table']['latitude'][1],g['observations_table']['longitude'][1])
#                         l+=1
#                     except:
                        
#                         print('table read error')
#         except:
#             print('file open error')
            
# print(l,' duplicates')
        

http = urllib3.PoolManager()
r = http.request('GET', 'http://early-upper-air.copernicus-climate.eu/statlist/?mindate=1900-01-01&enddate=2020-12-31')
fns=r.data.split(b'\n')
for i in range(len(fns)):
    fns[i]=fns[i].split(b',')[0].decode()
# opath=os.path.expandvars('/raid60/raid/home/srvx7/lehre/users/a1400070/adjust/Temperature_adjustment/files')
opath=os.path.expandvars('Temperature_adjustment/files2')

print(opath)
os.chdir(opath)
#fns=glob.glob('0?????/')
#fns=[fns[fns.index('0-20000-0-26781')]]
cdict={'ta':'temperatures','obs_minus_bg':'era5_fgdep','bias_estimate':'bias_estimate','lat':'lat','lon':'lon','hours':'hours'}
fnu=[]
fnd=[]
for fnf in fns:
    if fnf == fns[0]:
        continue
    fn=fnf[-5:]
    prefix='0'
    if fn in fnu:
        print('duplicate '+fnf+', incrementing leading zero to 1')
        prefix='1'
    fnu.append(fn)
    fo='feedbackmerged'+prefix+fn+'.nc'
    print(fo)
    try:
        
        mt=os.path.getmtime(fo)
        ts=time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(mt))
        #if ts.split(' ')[0][-2:]=='11':
            #continue
    except:
        pass
    
    try:        
        c = cdsapi.Client()
        r = c.retrieve('insitu-comprehensive-upper-air-observation-network',
                       {'variable': 'temperature',
                        'optional':['obs_minus_bg','bias_estimate'],
                        'statid': fn,
                        'newdl':'000111000',
                        'pressure_level':[10,20,30,50,70,100,150,200,250,300,400,500,700,850,925,1000]
                       }
                      )
        r.download(target='download.zip')
        assert os.stat('download.zip').st_size == r.content_length, "Downloaded file is incomplete"
        z = zipfile.ZipFile('download.zip')
        z.extractall(path='./downloaded/downloaded_'+fn)
        z.close()
        files = glob.glob('./downloaded/downloaded_'+ fn +'/*.nc')
        data=eua.CDMDataset(files[0])

            
#         data=eua.vm_request_wrapper({'variable': 'temperature', 'optional':['obs_minus_bg','bias_estimate'],'statid': fn, 'pressure_level':[1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000,92500,100000]}, 
#                                     overwrite=True,vm_url='http://srvx8.img.univie.ac.at:8002')
    except Exception as e:
        print(e)
        continue
    dq=data.read_data_to_3dcube(list(cdict.keys()))
    dql={}
    for k in dq.keys():
    ###dq[k].rename_dims({'plev':'pressure'})
        dql[cdict[k]]=dq[k].rename(cdict[k])
    xrdq=xr.Dataset(dql).rename_dims({'press':'pressure','datum':'time'})
    xrdq.attrs['unique_source_identifier']=fnf
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

#     os.mkdir(opath+prefix+fn)

#     try:
#         os.mkdir(opath+prefix+fn)
#         print('dir was made')
#     except:
#         pass

    xrdq.to_netcdf(path=fo, format='NETCDF4_CLASSIC')
    print('wrote '+fo)