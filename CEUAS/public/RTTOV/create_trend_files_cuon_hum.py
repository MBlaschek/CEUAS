import numpy
import numpy as np
import pandas
import pandas as pd
from numba import njit
import sys,glob
import zipfile, os, time
import urllib3
from datetime import datetime, timedelta
import glob
import h5py
import netCDF4 as nc
sys.path.append(os.getcwd()+'/../cds-backend/code/')
sys.path.append(os.getcwd()+'/../harvest/code/')
import rasotools
# from harvest_convert_to_netCDF_newfixes import write_dict_h5
import cds_eua4 as eua
# eua.logging_set_level(30)
import xarray as xr

import cdsapi, zipfile, os, time
#import schedule
import copy
from shutil import copyfile
import multiprocessing
import pickle

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.pylab as pylab
import warnings
warnings.filterwarnings('ignore')

from inspect import getmembers, isfunction


from IPython.display import Image
from IPython.core.display import HTML 
import rasotools


import matplotlib.pylab as plt
import matplotlib.colors as mcolors
import numpy

import ray
ray.init(num_cpus=15)

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


def datetime_to_seconds(dates, ref="1900-01-01T00:00:00"):
    """from datetime64 to seconds since 1900-01-01 00:00:00"""
    return ((dates - np.datetime64(ref)) / np.timedelta64(1, "s")).astype(np.int64)


def seconds_to_datetime(seconds, ref="1900-01-01"):
    """from seconds to datetime64"""
    seconds = np.asarray(seconds)
    return pd.to_datetime(seconds, unit="s", origin=ref, errors="coerce")

# create trend files from v13 rh only
@ray.remote
def calc_trend(file):
    trends_adj = {}
    trends_uadj = {}   
    start = datetime_to_seconds(pd.to_datetime(['1999-01-01']))[0] 
    end = datetime_to_seconds(pd.to_datetime(['2011-01-01']))[0]
    level = 30000

    with h5py.File(file, 'r') as file:
        if '138' in list(file['recordindices'].keys()):                                                                                 # check if there is humidity data for this station
            input = {}                                                                                                                  # setup input var for trend calculation
            lat_lon = str(file['observations_table']['latitude'][-1]) + '_' + str(file['observations_table']['longitude'][-1])          # create station identifier out of lat lon
            ts = file['recordindices']['recordtimestamp'][:]                                                                            # whole timestamp 
            idx = file['recordindices']['138'][:-1]                                                                                     # 138 for relative humidity
            idx = idx[np.logical_and(ts >= start, ts < end)]                                                                            # reduce to data between start and end
            if len(idx) < 1000:
                return trends_adj, trends_uadj
            z_coord = file['observations_table']['z_coordinate'][idx[0]:idx[-1]]                                                        # get pressure data
            zidx = idx[0] + np.where(z_coord == level)[0]                                                                               # create index, where on 50 000 Pa
            input['hums'] = file['observations_table']['observation_value'][zidx]                                                       # select relative humidity data
            try:
                input['adj_hums'] = input['hums'] - file['advanced_homogenisation']['humidity_bias_estimate'][zidx]                     # subtract adjustments for adjusted relative humidity
            except: return trends_adj, trends_uadj
            input['time'] = seconds_to_datetime(file['observations_table']['date_time'][zidx])                                          # convert time
        else:
            return trends_adj, trends_uadj                                                                                              # skip if no humidity data available
    df = pd.DataFrame.from_dict(input)                                                                                                  # create df from input
    df = df.dropna()
    if len(df) < 5*12*30:
        return trends_adj, trends_uadj
    df = df.resample('M', on='time').mean()                                                                                             # convert to monthly means
    if len(df) > ((10*12)-24):
        xdf = df.to_xarray()                                                                                                            # convert to xr and start trend calculation
        trends_adj[lat_lon] = float(rasotools.met.time.trend(xdf.hums,only_slopes=True)*3650)
        trends_uadj[lat_lon] = float(rasotools.met.time.trend(xdf.adj_hums,only_slopes=True)*3650)
    return trends_adj, trends_uadj

level = 30000
files = glob.glob('/mnt/users/scratch/leo/scratch/converted_v13/long/*v1.nc')

result_ids = []
for i in files[:]:
    result_ids.append(calc_trend.remote(i))
results = ray.get(result_ids)
ray.shutdown()

trends_adj = {}
trends_uadj = {}
for i in results:
    for j in i[0]:
        trends_adj[j] = i[0][j] 
        trends_uadj[j] = i[1][j] 

pickle.dump( trends_adj, open( "/users/staff/uvoggenberger/scratch/RTTOV_output/hum_worldmap/adj/hum_adj_"+str(level)+"cuonv13.p", "wb" ) )
pickle.dump( trends_uadj, open( "/users/staff/uvoggenberger/scratch/RTTOV_output/hum_worldmap/adj/hum_uadj_"+str(level)+"cuonv13.p", "wb" ) )