# System information
import numpy as np
import os, sys, glob
import xarray
print(sys.executable)
print(sys.version_info)
import pandas as pd
import pandas
pd.set_option('display.max_columns', None)
# pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)
sys.path.append(os.getcwd()+'/../cds-backend/code/')
import cds_eua3 as eua
import pickle
import multiprocessing
from functools import partial

def haversine_np(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)

    All args must be of equal length.    

    """
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = np.sin(dlat/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2

    c = 2 * np.arcsin(np.sqrt(a))
    km = 6367 * c
    return km


def calc_station(date, search_lat, search_lon, channel):
    yr = date[0]
    mn = date[1]
    print(yr, mn)
    folder = glob.glob('/users/staff/a1400070/CEUAS/CEUAS/public/testing/fiduceo/n16/dap.ceda.ac.uk/neodc/fiduceo/data/fcdr/microwave/v4.1/amsub/noaa16/'+str(yr)+'/'+"%02d" % (mn)+'/*/')
    counter = 0
    mean = 0
    for i in folder:
        try:
            files = glob.glob(i+'*.nc')
            df_night = None
            df_day = None
            for j in files:
                start = (j.split('_')[-5][-6:-4])
                ds = xarray.open_dataset(j)
                df = ds[channel].to_dataframe().dropna()
                df_day = df.append(df_day)
            dist_day = haversine_np(len(df_day)*[search_lon], len(df_day)*[search_lat], df_day.longitude, df_day.latitude)
            mean = mean + df_day[dist_day == dist_day.min()][channel].values
            counter = counter + 1
        except: pass
    try:
        mean = (mean/counter)
    except:
        mean = np.nan
    pickle.dump( mean, open( "./fiduceo/fiduceo_noaa16_"+str(chan)+"_"+str(yr)+"_"+"%02d" % (mn)+"_"+str(lat)+"_"+str(lon)+"_brightness_temperature.p", "wb" ) )


if __name__ == '__main__': 
    lat = 48.2
    lon = 16.37
    chan = 'Ch20_BT'
    try:
        os.makedirs("./fiduceo/")
    except:
        pass
    
    dates = []
    for yr in range(1999, 2012, 1):
        for mn in range(1, 13, 1):
            dates.append([yr, mn])

    pool = multiprocessing.Pool(processes=30)
    func=partial(calc_station, search_lat = lat, search_lon = lon, channel = chan)
    result_list = list(pool.map(func, dates))
    print(result_list)