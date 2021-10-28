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


def calc_station(date, search_lat, search_lon, channel, stat):
    yr = date[0]
    mn = date[1]
    print(yr, mn)
    folder = glob.glob('/users/staff/a1400070/CEUAS/CEUAS/public/testing/fiduceo/dap.ceda.ac.uk/neodc/fiduceo/data/fcdr/microwave/v4.1/mhs/noaa19/'+str(yr)+'/'+"%02d" % (mn)+'/*/')
#     folder = glob.glob('/users/staff/a1400070/CEUAS/CEUAS/public/testing/fiduceo/dap.ceda.ac.uk/neodc/fiduceo/data/fcdr/microwave/v4.1/amsub/noaa15/'+str(yr)+'/'+"%02d" % (mn)+'/*/')
#     folder = glob.glob('/users/staff/a1400070/CEUAS/CEUAS/public/testing/fiduceo/n16/dap.ceda.ac.uk/neodc/fiduceo/data/fcdr/microwave/v4.1/amsub/noaa16/'+str(yr)+'/'+"%02d" % (mn)+'/*/')
    counter = [0]*len(stat)
    mean = [0]*len(stat)
    times = [[]]*len(stat)
    dists = [[]]*len(stat)
    for i in folder: # for each day of selected month
        try:
            files = glob.glob(i+'*.nc')
            df_night = None
            df_day = None

    #         ds = xarray.open_mfdataset(files, combine='nested') 
    #         # ValueError: arguments without labels along dimension 'y' cannot be aligned because they have different dimension sizes: {2304, 2303}

            for j in files:
                ds = xarray.open_dataset(j)
                df = ds[[channel,'Time']].to_dataframe().dropna()
                df_day = df.append(df_day)
                
            for k in range(len(stat)):
                dist_day = haversine_np(len(df_day)*[search_lon[k]], len(df_day)*[search_lat[k]], df_day.longitude, df_day.latitude)
                mindist = dist_day.min()
                selection = df_day[dist_day == mindist]
                mean[k] = mean[k] + selection[channel].values
                times[k].append(selection.Time.values)
                dists[k].append(mindist)
                counter[k] = counter[k] + 1
            
        except: pass
    for k in range(len(stat)):
        try:
            out_mean = (mean[k]/counter[k])
            out_times = times[k]
            out_dists = dists[k]
        except:
            out_mean = np.nan
            out_times = [np.nan]
            out_dists = [np.nan]
            
        try:
            os.makedirs("./fiduceo/out/all_noaa19/"+stat[k])
        except:
            pass
        
        pickle.dump( [out_mean, out_times, out_dists], open( "./fiduceo/out/all_noaa19/"+stat[k]+"/fiduceo_mhs_noaa19_"+str(channel)+"_"+str(yr)+"_"+"%02d" % (mn)+"_"+str(search_lat[k])+"_"+str(search_lon[k])+"_brightness_temperature.p", "wb" ) )
        
    
    
#     pickle.dump( mean, open( "./fiduceo/fiduceo_amsub_noaa15_"+str(chan)+"_"+str(yr)+"_"+"%02d" % (mn)+"_"+str(lat)+"_"+str(lon)+"_brightness_temperature.p", "wb" ) )
#     pickle.dump( mean, open( "./fiduceo/fiduceo_amsub_noaa16_"+str(chan)+"_"+str(yr)+"_"+"%02d" % (mn)+"_"+str(lat)+"_"+str(lon)+"_brightness_temperature.p", "wb" ) )


if __name__ == '__main__': 
#     selstats = glob.glob('/users/staff/a1400070/CEUAS/CEUAS/public/testing/rttov_out_hum_noaa_19_mhs/*')
#     stats = []
#     for i in selstats:
#         stats.append(glob.glob('/mnt/users/scratch/leo/scratch/converted_v7/*'+ i.split('/')[-1] +'*_CEUAS_merged_v1.nc')[0])
        
#     lats = []
#     lons = []
#     names = []
#     for i in stats:
#         names.append(i.split('-')[-1][:5])
#         df = eua.CDMDataset(filename = i).to_dataframe(groups=['observations_table'], variables=['latitude', 'longitude'])
#         lats.append(df.latitude.iloc[0])
#         lons.append(df.longitude.iloc[0])
#     print('lat, lon collected')
    
    try:
        os.makedirs("./fiduceo/")
    except:
        pass
    
    dates = []
    for yr in range(2009, 2018, 1):
        for mn in range(1, 13, 1):
            dates.append([yr, mn])

    lats, lons, names = pickle.load(open("./stationlist.p", "rb" ))
    chan = 'Ch3_BT'
    pool = multiprocessing.Pool(processes=40)
    func=partial(calc_station, search_lat = lats, search_lon = lons, channel = chan, stat = names)
    result_list = list(pool.map(func, dates))
    print(result_list)
