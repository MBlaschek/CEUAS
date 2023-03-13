

import numpy
import numpy as np
import time
import datetime
import netCDF4
import matplotlib.pylab as plt
import os,sys,glob
from multiprocessing import Pool
#import odb
from eccodes import *
from functools import partial
from collections import OrderedDict
import subprocess
import json
import gzip
import copy
import pickle
import xarray as xr
import pickle
import pandas as pd
import matplotlib
import matplotlib.pylab as plt
import matplotlib.pyplot as maplt
matplotlib.rcParams.update({'font.size': 20})
import hdf5plugin
from tqdm import tqdm

plt.rcParams['lines.linewidth'] = 3

import warnings
warnings.filterwarnings('ignore')

sys.path.append(os.getcwd()+'/../cds-backend/code/')
import cds_eua4 as eua

import h5py
import ray
import pyproj
# ray.init(num_cpus=7)

import urllib
import pycountry

import ray
from tqdm import tqdm

def datetime_to_seconds(dates, ref='1900-01-01T00:00:00'):
    """ from datetime64 to seconds since 1900-01-01 00:00:00"""
    return ((dates - np.datetime64(ref)) / np.timedelta64(1, 's')).astype(np.int64)

def seconds_to_datetime(seconds, ref='1900-01-01'):
    """ from seconds to datetime64 """
    seconds = np.asarray(seconds)
    return pd.to_datetime(seconds, unit='s', origin=ref)


def trend_station(i):
    sys.path.append(os.getcwd()+'/../resort/rasotools-master/')
    import rasotools
    df_dict = {}
    sout = []
    aout = []
    stats = []
    lats = []
    lons = []
    start = 1973
    end = 2003
    intervall = end - start
    dt_from = datetime_to_seconds(np.datetime64('1973-01-01'))
    dt_to = datetime_to_seconds(np.datetime64('2002-12-31'))
    
    try:
        # multiprocessing here
        with h5py.File(i, 'r') as file:
            rts = file['recordindices']['recordtimestamp'][:]
            idx = np.where(np.logical_and((rts >= dt_from), (rts <= dt_to)))[0]
            plevs = [70000]

            idx_d = {}
            # var_d = {'air_temperature':'126', 'relative_humidty':'138', 'geopotential':'117', 'eastward_wind_speed':'139', 'northward_wind_speed':'140', 'dew_point': '137', 'specific_humidity':'39'}
            var_d = {'relative_humidty':'138'}
            for j in var_d:
                idx_d[j] = file['recordindices'][var_d[j]][idx]

            masks = {}
            for j in idx_d:
                masks[j] = file['observations_table']['z_coordinate'][idx_d[j][0]:idx_d[j][-1]]
                masks[j] = np.isin(masks[j],plevs)
                # masks[i] = np.isfinite(masks[i])

            mask = masks['relative_humidty']
            t_idx = idx_d['relative_humidty']
            df_dict['observation_value'] = list(file['observations_table']['observation_value'][t_idx[0]:t_idx[-1]][mask])
            df_dict['z_coordinate'] = list(file['observations_table']['z_coordinate'][t_idx[0]:t_idx[-1]][mask])
            df_dict['date_time'] = seconds_to_datetime(list(file['observations_table']['date_time'][t_idx[0]:t_idx[-1]][mask]))
            df_dict['latitude'] = list(file['observations_table']['latitude'][t_idx[0]:t_idx[-1]][mask])
            df_dict['longitude'] = list(file['observations_table']['longitude'][t_idx[0]:t_idx[-1]][mask])
            df_dict['humidity_bias_estimate'] = list(file['advanced_homogenisation']['humidity_bias_estimate'][t_idx[0]:t_idx[-1]][mask])

            temp = pd.DataFrame.from_dict(df_dict)
            # display(temp)

            if len(temp) > 0:
                # temp.sort_values('time')
                temp['time'] = pd.to_datetime(temp['date_time'])
                temp['lat'] = numpy.array([temp.latitude.iloc[-1]]*len(temp))
                temp['lon'] = numpy.array([temp.longitude.iloc[-1]]*len(temp))
                temp['adjusted'] = temp['observation_value'] - temp['humidity_bias_estimate']
                temptime = temp.time
                if len(temp) >= 19*365 and len(numpy.unique(temptime.dt.year)) > 19 :
                    # print('enough data')
                    xa = temp.set_index(['lat', 'lon', 'time']).to_xarray()
                    # and do it twice for the adjusted values too!
                    out = rasotools.met.time.trend(xa.observation_value,only_slopes=True).to_dataframe(name='out')
                    sout.append(float(out.iloc[-1] *3650))
                    out_adj = rasotools.met.time.trend(xa.adjusted,only_slopes=True).to_dataframe(name='out_adj')
                    aout.append(float(out_adj.iloc[-1] *3650))
                else:
                    sout=np.nan
                    aout=np.nan
            else:
                sout=np.nan
                aout=np.nan
            try:
                lats=temp.latitude.iloc[-1]
                lons=temp.longitude.iloc[-1]
                stats=i
            except:
                lats=np.nan
                lons=np.nan
                stats=i
    except:
        lats=np.nan
        lons=np.nan
        stats=i
        sout=np.nan
        aout=np.nan
    return [stats, lats, lons, sout, aout]




def to_iterator(obj_ids):
    while obj_ids:
        done, obj_ids = ray.wait(obj_ids)
        yield ray.get(done[0])
        

files =  glob.glob('/mnt/users/scratch/leo/scratch/converted_v11/long/*.nc')
test_r = ray.remote(trend_station)
ray.init(num_cpus=20)
results = []
obj_ids = [test_r.remote(i) for i in files]
for x in tqdm(to_iterator(obj_ids), total=len(obj_ids)):
    results.append(x)

pickle.dump( results, open( "trends_700hPa_1973_2003_Trend.p", "wb" ))
ray.shutdown()
