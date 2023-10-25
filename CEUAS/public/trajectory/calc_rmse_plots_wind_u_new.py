#!/usr/bin/env
# coding: utf-8

import h5py
import trajectory as trj
import numpy as np
import pandas as pd
import os,sys,glob
import copy
import time
import xarray as xr
import pickle

import matplotlib
import matplotlib.pylab as plt
import matplotlib.pyplot as maplt
matplotlib.rcParams.update({'font.size': 20})

plt.rcParams['lines.linewidth'] = 3

import warnings
warnings.filterwarnings('ignore')

import ray


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def datetime_to_seconds(dates, ref='1900-01-01T00:00:00'):
    """ from datetime64 to seconds since 1900-01-01 00:00:00"""
    return ((dates - np.datetime64(ref)) / np.timedelta64(1, 's')).astype(np.int64)

def seconds_to_datetime(seconds, ref='1900-01-01'):
    """ from seconds to datetime64 """
    seconds = np.asarray(seconds)
    return pd.to_datetime(seconds, unit='s', origin=ref)


@ray.remote
def calc_station(sid, year, era_input_files):
    show_date = False
    diff = True
    stat = sid
    compare_to = 'fc' # fc
    
    maxtimediff = pd.Timedelta(hours=2)
    
    dt_from = datetime_to_seconds(np.datetime64(str(year)+'-01-01'))
    dt_to = datetime_to_seconds(np.datetime64(str(year)+'-12-31'))
    
    conv_file = glob.glob('/mnt/users/scratch/leo/scratch/converted_v9/newindex/*' + stat + '*_CEUAS_merged_v1.nc')[0]
#     print(conv_file)
#     print(dt_from)
#     print(dt_to)
    df_dict = {}
    h_df_dict = {}

    stdplevs = [1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000,92500]
    stat_marker = {}
    rmse_sum_shbase_sonde={}
    rmse_sum_shbase_adjsonde={}
    rmse_sum_shdisp_sonde={}
    rmse_sum_shdisp_adjsonde={}
    rms_sum_shbase={}
    rms_sum_adjsonde={}
    rms_sum_sonde={}
    rms_sum_shdisp={}
    rms_sum_dispminusbase={}

    for i in stdplevs:
        stat_marker[i] = []
        rmse_sum_shbase_sonde[i] = []
        rmse_sum_shbase_adjsonde[i] = []
        rmse_sum_shdisp_sonde[i] = []
        rmse_sum_shdisp_adjsonde[i] = []
        rms_sum_shbase[i] = []
        rms_sum_adjsonde[i] = []
        rms_sum_sonde[i] = []
        rms_sum_shdisp[i] = []
        rms_sum_dispminusbase[i] = []
    try:
        with h5py.File(conv_file, 'r') as file:
            rts = file['recordindices']['recordtimestamp'][:]
            idx = np.where(np.logical_and((rts >= dt_from), (rts <= dt_to)))[0]
            if len(idx) == 0:
                print('(1) NO DATA FOUND IN CONVERTED_V9: ', sid)
                return [rmse_sum_shbase_sonde, rmse_sum_shbase_adjsonde, rmse_sum_shdisp_sonde,
                        rmse_sum_shdisp_adjsonde, rms_sum_shbase, rms_sum_adjsonde, rms_sum_sonde,
                        rms_sum_shdisp, rms_sum_dispminusbase, stat_marker]

            h_idx = file['recordindices']['139'][idx]
            t_idx = file['recordindices']['126'][idx]
            plevs = [1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000,92500,100000]


            mask = file['observations_table']['z_coordinate'][t_idx[0]:t_idx[-1]]
            mask = np.isin(mask,plevs)

            h_mask = file['observations_table']['z_coordinate'][h_idx[0]:h_idx[-1]]
            h_mask = np.isin(h_mask,plevs)

            t_len = len(mask[mask == True])

            df_dict['z_coordinate'] = list(file['observations_table']['z_coordinate'][t_idx[0]:t_idx[-1]][mask])
            df_dict['date_time'] = list(file['observations_table']['date_time'][t_idx[0]:t_idx[-1]][mask])
    #             df_dict['observation_value'] = list(file['observations_table']['observation_value'][t_idx[0]:t_idx[-1]][mask])
            df_dict['latitude'] = list(file['observations_table']['latitude'][t_idx[0]:t_idx[-1]][mask])
            df_dict['longitude'] = list(file['observations_table']['longitude'][t_idx[0]:t_idx[-1]][mask])
            repid = np.asarray(file['observations_table']['report_id'][t_idx[0]:t_idx[-1]][mask])
            df_dict['report_id'] = list(repid.view('|S{}'.format(repid.shape[1])).flatten().astype(str))

            df_dict['RASE_bias_estimate'] = list(file['advanced_homogenisation']['RASE_bias_estimate'][t_idx[0]:t_idx[-1]][mask])
            df_dict['latitude_displacement'] = list(file['advanced_homogenisation']['latitude_displacement'][t_idx[0]:t_idx[-1]][mask])
            df_dict['longitude_displacement'] = list(file['advanced_homogenisation']['longitude_displacement'][t_idx[0]:t_idx[-1]][mask])
            df_dict['time_since_launch'] = list(file['advanced_homogenisation']['time_since_launch'][t_idx[0]:t_idx[-1]][mask])

            h_df_dict['z_coordinate'] = list(file['observations_table']['z_coordinate'][h_idx[0]:h_idx[-1]][h_mask])
            h_df_dict['date_time'] = list(file['observations_table']['date_time'][h_idx[0]:h_idx[-1]][h_mask])
            h_df_dict['observation_value'] = list(file['observations_table']['observation_value'][h_idx[0]:h_idx[-1]][h_mask])

            df_dict['variable'] = ['specific_humidity']*t_len

            df_dict['date_time'] = seconds_to_datetime(df_dict['date_time'])
            df = pd.DataFrame.from_dict(df_dict)

            h_df_dict['date_time'] = seconds_to_datetime(h_df_dict['date_time'])
            h_df = pd.DataFrame.from_dict(h_df_dict)

            # put dfs together:
            df = df.merge(h_df, how='inner', on=['date_time','z_coordinate'])
#         print(df)
    except:
        print('(0) NO DATA FOUND IN CONVERTED_V9: ', sid)
        return [rmse_sum_shbase_sonde, rmse_sum_shbase_adjsonde, rmse_sum_shdisp_sonde,
                rmse_sum_shdisp_adjsonde, rms_sum_shbase, rms_sum_adjsonde, rms_sum_sonde,
                rms_sum_shdisp, rms_sum_dispminusbase, stat_marker]
    df = df.dropna()
    if len(df) == 0:
        print('(2) NO DATA FOUND IN CONVERTED_V9: ', sid)
        return [rmse_sum_shbase_sonde, rmse_sum_shbase_adjsonde, rmse_sum_shdisp_sonde,
                rmse_sum_shdisp_adjsonde, rms_sum_shbase, rms_sum_adjsonde, rms_sum_sonde,
                rms_sum_shdisp, rms_sum_dispminusbase, stat_marker]

    for mon in [1,2,3,4,5,6,7,8,9,10,11,12]:

        df_mon = df[df.date_time.dt.month == mon]
    #     display(df_mon)
        t0 = time.time()
        ds_fc = era_input_files[mon]
        t0 = time.time()

        t0 = time.time()       
        for day in df_mon.date_time.drop_duplicates()[:]:

            input_data = df_mon[df_mon.date_time == day]
            ds_fc_time = ds_fc.sel(time=day, method='nearest')
#             print(day, pd.Timestamp(ds_fc_time.time.values))
            if (pd.Timestamp(ds_fc_time.time.values) - day) > maxtimediff:
                continue
            t_list = []
            for i in np.array(ds_fc_time.level): #10,20,...,1000
                step = find_nearest(input_data.z_coordinate, i*100)
                input_data_step = input_data[input_data.z_coordinate == step]
                station_lat = input_data.latitude.iloc[0] + np.array(input_data_step.latitude_displacement)[0]
                station_lon = input_data.longitude.iloc[0] + np.array(input_data_step.longitude_displacement)[0]
                lon = station_lon
                if lon < 0:
                    lon = 360.+lon
                ds_now = ds_fc_time.interp(latitude=[station_lat], longitude=[lon], method="linear")
                t = ds_now.u.sel(level = i)
                t_list.append(float(t))

            p_ml = [1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000,92500,100000]
            lon = input_data.longitude.iloc[0]
            if lon < 0:
                lon = 360.+lon
            base_t = np.array(ds_fc_time.interp(latitude=[input_data.latitude.iloc[0]], longitude=[lon], method="linear").u)
            for i in range(len(stdplevs)):
                if np.abs(stdplevs[i] - find_nearest(input_data.z_coordinate,stdplevs[i])) > 500:
                    stat_marker[stdplevs[i]].append(np.nan)
                    rmse_sum_shbase_sonde[stdplevs[i]].append(np.nan)
                    rmse_sum_shbase_adjsonde[stdplevs[i]].append(np.nan)
                    rmse_sum_shdisp_sonde[stdplevs[i]].append(np.nan)
                    rmse_sum_shdisp_adjsonde[stdplevs[i]].append(np.nan)
                    rms_sum_shbase[stdplevs[i]].append(np.nan)
                    rms_sum_adjsonde[stdplevs[i]].append(np.nan)
                    rms_sum_sonde[stdplevs[i]].append(np.nan)
                    rms_sum_shdisp[stdplevs[i]].append(np.nan)
                    rms_sum_dispminusbase[stdplevs[i]].append(np.nan)

                else:
                    sq_t = np.squeeze(base_t)
                    t_base = float(sq_t[p_ml == find_nearest(p_ml,stdplevs[i])])
                    t_disp = float(np.array(t_list)[p_ml == find_nearest(p_ml,stdplevs[i])])
                    input_data_step = input_data[input_data.z_coordinate == find_nearest(input_data.z_coordinate, stdplevs[i])]
                    t_sonde = float(input_data_step.observation_value)
        #             t_adjsonde = float(input_data_step.temperature) - float(input_data_step['fg_depar@body']) - float(input_data_step['an_depar@body'])
                    t_adjsonde = float(input_data_step.observation_value - input_data_step.RASE_bias_estimate)# - float(input_data_step['an_depar@body'])

                    stat_marker[stdplevs[i]].append(stat)
                    rmse_sum_shbase_sonde[stdplevs[i]].append(t_base - t_sonde)
                    rmse_sum_shbase_adjsonde[stdplevs[i]].append(t_base - t_adjsonde)
                    rmse_sum_shdisp_sonde[stdplevs[i]].append(t_disp - t_sonde)
                    rmse_sum_shdisp_adjsonde[stdplevs[i]].append(t_disp - t_adjsonde)
                    rms_sum_shbase[stdplevs[i]].append(t_base)
                    rms_sum_adjsonde[stdplevs[i]].append(t_adjsonde)
                    rms_sum_sonde[stdplevs[i]].append(t_sonde)
                    rms_sum_shdisp[stdplevs[i]].append(t_disp)
                    rms_sum_dispminusbase[stdplevs[i]].append(t_disp-t_base)

    with open('rmse_h_out/era5_h_fc_'+str(year)+'_'+str(sid)+'_rmse_data.p', 'wb') as file:
        pickle.dump(
            [rmse_sum_shbase_sonde, rmse_sum_shbase_adjsonde, rmse_sum_shdisp_sonde,
            rmse_sum_shdisp_adjsonde, rms_sum_shbase, rms_sum_adjsonde, rms_sum_sonde,
            rms_sum_shdisp, rms_sum_dispminusbase, stat_marker],
            file
        )
        print('calculating mon: ',mon, time.time() - t0)
    return [rmse_sum_shbase_sonde, rmse_sum_shbase_adjsonde, rmse_sum_shdisp_sonde,
            rmse_sum_shdisp_adjsonde, rms_sum_shbase, rms_sum_adjsonde, rms_sum_sonde,
            rms_sum_shdisp, rms_sum_dispminusbase, stat_marker]
 
    

if __name__ == '__main__':
#     year = 1990
    stdplevs = [1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000,92500]
    diff = True
    show_date = False
    for year in range(2010,2011,1): #[2000]: #  [1960, 1970, 1980, 1990, 2000, 2010, 2020]:
        t0 =  time.time()
        ray.init(num_cpus=20)
        ds_fcs = {}
        for mon in [1,2,3,4,5,6,7,8,9,10,11,12]:
            # if compare_to == 'fc':
            files = glob.glob('/mnt/scratch/scratch/leo/scratch/era5/gridded/era5fct.'+str(year)+str(mon).zfill(2)+'*.131.nc')[0]
            ds_fcs[mon] = xr.load_dataset(files)
        era_input_files = ray.put(ds_fcs)
        print('loading era5 data: ', time.time() - t0)
        
        file_list = []
        for i in glob.glob('/mnt/users/scratch/leo/scratch/converted_v9/newindex/*.nc')[:]:
            print(i)
            sid = i.split('_CEUAS_merged_v1.nc')[0][-5:]
            file_list.append(calc_station.remote(sid,year,era_input_files))
        results = ray.get(file_list)
                
        with open('era5_u_fc_'+str(year)+'_rmse_data.p', 'wb') as file:
            pickle.dump(results, file)
        print('final time: ', time.time() - t0)
        ray.shutdown()
            
