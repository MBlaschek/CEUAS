#!/usr/bin/env python
import traceback
import sys,glob
import os.path

import numpy
from datetime import date
import netCDF4
import time
from numba import njit
# from rasotools.utils import *
# from rasotools.anomaly import *
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


http = urllib3.PoolManager()
r = http.request('GET', 'http://early-upper-air.copernicus-climate.eu/statlist/?mindate=1900-01-01&enddate=2020-12-31')
fns=r.data.split(b'\n')
for i in range(len(fns)):
    fns[i]=fns[i].split(b',')[0].decode()
opath=os.path.expandvars('Wind_adjustment')

print(opath)
os.chdir(opath)
fnu=[]
fnd=[]
for fn in ['35229', '68994', '11035', '10393', '91413', '70219']:
# for fnf in fns:
#     if fnf == fns[0]:
#         continue
#     fn=fnf[-5:]
    try:        
        c = cdsapi.Client()
        r = c.retrieve('insitu-comprehensive-upper-air-observation-network',
                       {'variable': ['eastward_wind_speed', 'northward_wind_speed','wind_from_direction', 'wind_speed'],
                        'optional':['obs_minus_bg','bias_estimate'],
                        'statid': fn,
                        'skip':'7888',
                        'pressure_level':[10,20,30,50,70,100,150,200,250,300,400,500,700,850,925,1000]
                       }
                      )
        r.download(target='download.zip')
        assert os.stat('download.zip').st_size == r.content_length, "Downloaded file is incomplete"
        z = zipfile.ZipFile('download.zip')
        z.extractall(path='./downloaded/wind_downloaded_'+fn)
        z.close()
#         files = glob.glob('./downloaded/downloaded_'+ fn +'/*.nc')
#         data=eua.CDMDataset(files[0])
    except:
        pass

# a    