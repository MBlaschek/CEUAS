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
import time
from numba import njit
import h5py

@njit(cache=True)
def haversine_np(lon1, lat1, lon2, lat2): #lon1,lat1 are scalars; lon2, lat2 are arrays
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)

    All args must be of equal length.    

    """
    abslat=np.abs(lat2-lat1)
    idx=np.where(abslat<0.5)[0]
        
    rlon1=np.radians(lon1)
    rlat1=np.radians(lat1)
    rlon2=np.radians(lon2[idx])
    rlat2=np.radians(lat2[idx])

    dlon = rlon2 - rlon1
    dlat = rlat2 - rlat1

    a = np.sin(dlat/2.0)**2 + np.cos(rlat1) * np.cos(rlat2) * np.sin(dlon/2.0)**2

    c = 2 * np.arcsin(np.sqrt(a))
    ic=np.argmin(c)
    #if abs(rlat2[ic]-rlat1)*180/np.pi>0.1:
        
        #print(abs(rlat2[ic]-rlat1)*180/np.pi,c[ic],np.min(abslat),rlat2[ic],rlat1)

    return idx[ic],6367*c[ic]


def calc_station(date, search_lat, search_lon, channel, stat):
    tt=time.time()
    yr = date[0]
    mn = date[1]
    print(yr, mn)
    folder = glob.glob('/users/staff/a1400070/CEUAS/CEUAS/public/testing/fiduceo/dap.ceda.ac.uk/neodc/fiduceo/data/fcdr/microwave/v4.1/mhs/noaa19/'+str(yr)+'/'+"%02d" % (mn)+'/*/')
#     folder = glob.glob('/users/staff/a1400070/CEUAS/CEUAS/public/testing/fiduceo/dap.ceda.ac.uk/neodc/fiduceo/data/fcdr/microwave/v4.1/amsub/noaa15/'+str(yr)+'/'+"%02d" % (mn)+'/*/')
#     folder = glob.glob('/users/staff/a1400070/CEUAS/CEUAS/public/testing/fiduceo/n16/dap.ceda.ac.uk/neodc/fiduceo/data/fcdr/microwave/v4.1/amsub/noaa16/'+str(yr)+'/'+"%02d" % (mn)+'/*/')
    counter = [0]*len(stat)
    mean3 = [0]*len(stat)
    mean4 = [0]*len(stat)
    mean5 = [0]*len(stat)
    times = [[]]*len(stat)
    dists = [[]]*len(stat)
    for i in folder: # for each day of selected month
        try:
            files = glob.glob(i+'*.nc')
            df_night = None
            df_day = None
            print(i.split('/')[-3:],time.time()-tt)

    #         ds = xarray.open_mfdataset(files, combine='nested') 
    #         # ValueError: arguments without labels along dimension 'y' cannot be aligned because they have different dimension sizes: {2304, 2303}

            df_day={'Ch3_BT':[], 'Ch4_BT':[], 'Ch5_BT':[], 'Time':[],'longitude':[],'latitude':[]}
            for j in files:
                #ds = xarray.open_dataset(j)
                #df = ds[['Ch3_BT', 'Ch4_BT', 'Ch5_BT', 'Time']].to_dataframe().dropna()
                #df_day = df.append(df_day)
                with h5py.File(j,'r') as f:
                    mask=None
                    for fvar in df_day.keys():
                        if 'scale_factor' in f[fvar].attrs.keys():
                            v=f[fvar][:]
                            if mask is None:
                                mask=v!=f[fvar].attrs['_FillValue']
                            vals=v[mask]*f[fvar].attrs['scale_factor']
                            
                            df_day[fvar].append(vals)
                        else:
                            df_day[fvar].append(np.asarray([f[fvar][:]]*90).T[mask])
            if not df_day['Ch3_BT']:
                continue
            for fvar in df_day.keys():
                df_day[fvar]=np.concatenate(df_day[fvar])
                
    
            print('nach read',time.time()-tt)
            
            for k in range(len(stat)):
                #dist_day = haversine_np(len(df_day)*[search_lon[k]], len(df_day)*[search_lat[k]], df_day.longitude, df_day.latitude)
                imin,dmin = haversine_np(search_lon[k],search_lat[k], df_day['longitude'], df_day['latitude'])
                if dmin>300.:
                    print(k,stat[k],dmin)
                    continue
                #print(k,'dist',time.time()-tt)
                #mindist = dist_day.min()
                #selection = df_day[dist_day == mindist]
                #mean3[k] = mean3[k] + selection['Ch3_BT'].values
                #mean4[k] = mean4[k] + selection['Ch4_BT'].values
                #mean5[k] = mean5[k] + selection['Ch5_BT'].values
                #times[k].append(selection.Time.values)
                #imin=np.argmin(dist_day)
                mean3[k] = mean3[k] + df_day['Ch3_BT'][imin]
                mean4[k] = mean4[k] + df_day['Ch4_BT'][imin]
                mean5[k] = mean5[k] + df_day['Ch5_BT'][imin]
                times[k].append(df_day['Time'][imin])
                dists[k].append(dmin)
                counter[k] = counter[k] + 1
                #print(k,time.time()-tt)
            print('nach stat',time.time()-tt)
            print('')
            
            
        except MemoryError: pass
    for k in range(len(stat)):
        try:
            out_mean = [mean3[k]/counter[k], mean4[k]/counter[k], mean5[k]/counter[k]]
            out_times = times[k]
            out_dists = dists[k]
        except:
            out_mean = [np.nan, np.nan, np.nan]
            out_times = [np.nan]
            out_dists = [np.nan]
            
        try:
            os.makedirs("./fiduceo/out/all_noaa19/"+stat[k])
        except:
            pass
        
        pickle.dump( [out_mean, out_times, out_dists], open( "./fiduceo/out/all_noaa19/"+stat[k]+"/fiduceo_mhs_noaa19_"+str(yr)+"_"+"%02d" % (mn)+"_"+str(search_lat[k])+"_"+str(search_lon[k])+"_brightness_temperature.p", "wb" ) )

    dt=time.time()-tt
    if dt>1:
        
        print(time.time()-tt)
        print('x')
    
    
#     pickle.dump( mean, open( "./fiduceo/fiduceo_amsub_noaa15_"+str(chan)+"_"+str(yr)+"_"+"%02d" % (mn)+"_"+str(lat)+"_"+str(lon)+"_brightness_temperature.p", "wb" ) )
#     pickle.dump( mean, open( "./fiduceo/fiduceo_amsub_noaa16_"+str(chan)+"_"+str(yr)+"_"+"%02d" % (mn)+"_"+str(lat)+"_"+str(lon)+"_brightness_temperature.p", "wb" ) )


if __name__ == '__main__': 
    #selstats = glob.glob('/users/staff/a1400070/CEUAS/CEUAS/public/testing/rttov_out_hum_noaa_19_mhs/*')
    #stats = []
    #for i in selstats:
        #stats.append(glob.glob('/mnt/users/scratch/leo/scratch/converted_v7/*'+ i.split('/')[-1] +'*_CEUAS_merged_v1.nc')[0])
        
    #lats = []
    #lons = []
    #names = []
    #for i in stats:
        #names.append(i.split('-')[-1][:5])
        #df = eua.CDMDataset(filename = i).to_dataframe(groups=['observations_table'], variables=['latitude', 'longitude'])
        #lats.append(df.latitude.iloc[0])
        #lons.append(df.longitude.iloc[0])
    #print('lat, lon collected')
    
    #pickle.dump([lats,lons,names],open("./stationlist.p", "wb" ))
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
    pool = multiprocessing.Pool(processes=20)
    func=partial(calc_station, search_lat = lats, search_lon = lons, channel = chan, stat = names)
    result_list = list(map(func, dates))
    print(result_list)
