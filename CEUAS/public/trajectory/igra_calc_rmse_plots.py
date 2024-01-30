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

sys.path.insert(0, "/users/staff/uvoggenberger/uvpy/rasotools-master/")
import rasotools

plt.rcParams['lines.linewidth'] = 3

import warnings
warnings.filterwarnings('ignore')

import ray
ray.init(num_cpus=60)

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
def calc_station(sid, year, var, ds_fc_store, selected_mons = None, input_file = 'ERA5', bg_depar_check = True):
    import rs_drift.drift as trj
    sys.path.insert(0, "/users/staff/uvoggenberger/uvpy/rasotools-master/")
    import rasotools

    print(sid)
    
    show_date = False
    diff = True
    stat = sid

    varsel_dict = {'eastward windspeed':'u', 'northward windspeed':'v', 'air temperature':'t', 'specific humidity': 'q'}
    era_dict = {'t':'130','u':'131','v':'132','q':'133'}

    varsel = varsel_dict[var]
    varsel_era = era_dict[varsel]

    if selected_mons == None:
        selected_mons = [1,2,3,4,5,6,7,8,9,10,11,12]

    maxtimediff = pd.Timedelta(hours=2)

    dt_from = datetime_to_seconds(np.datetime64(str(year)+'-01-01'))
    dt_to = datetime_to_seconds(np.datetime64(str(year)+'-12-31'))

    if input_file == 'ERA5':
        if year > 1979:
            conv_file_igra = glob.glob('/scratch/das/federico/COP2_HARVEST_JAN2023/era5_1/*'+ stat +'.txt.gz.nc')
        else:
            conv_file_igra = glob.glob('/scratch/das/federico/COP2_HARVEST_JAN2023/era5_2/*'+ stat +'.txt.gz.nc')
    else:
        conv_file_igra = glob.glob('/mnt/users/scratch/leo/scratch/converted_v13/long/*'+ stat +'*CEUAS_merged_v1.nc')
    
    df_dict = {}
    df_dict_w = {}
    df_dict_h = {}

    stdplevs = [1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000,92500,100000]
    rmse_sum_shbase_sonde={}
    rmse_sum_shdisp_sonde={}
    rms_sum_shbase={}
    rms_sum_sonde={}
    rms_sum_shdisp={}
    rms_sum_dispminusbase={}

    for i in stdplevs:
        rmse_sum_shbase_sonde[i] = []
        rmse_sum_shdisp_sonde[i] = []
        rms_sum_shbase[i] = []
        rms_sum_sonde[i] = []
        rms_sum_shdisp[i] = []
        rms_sum_dispminusbase[i] = []


    try:
        
        with h5py.File(conv_file_igra[0], 'r') as file:
            if input_file == 'ERA5':
                rts = file['recordtimestamp'][:]
                idx = np.where(np.logical_and((rts >= dt_from), (rts <= dt_to)))[0]
                if len(idx) == 0:
                    print('NO DATA FOUND IN IGRA: ', sid)
                    return [rmse_sum_shbase_sonde, rmse_sum_shdisp_sonde,
                            rms_sum_shbase, rms_sum_sonde,
                            rms_sum_shdisp, rms_sum_dispminusbase]
                t_idx = file['recordindex'][idx]
                plevs = [1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000,92500,100000]

                p_mask = file['observations_table']['z_coordinate'][t_idx[0]:t_idx[-1]]
                v_mask = file['observations_table']['observed_variable'][t_idx[0]:t_idx[-1]]
                mask_t = np.logical_and(np.isin(p_mask,plevs), np.isin(v_mask, [126]))
                mask_wd = np.logical_and(np.isin(p_mask,plevs), np.isin(v_mask, [106]))
                mask_ws = np.logical_and(np.isin(p_mask,plevs), np.isin(v_mask, [107]))
                mask_rh = np.logical_and(np.isin(p_mask,plevs), np.isin(v_mask, [138]))
                t_len = len(mask_t[mask_t == True])

                # wind data
                df_dict_w['z_coordinate'] = list(file['observations_table']['z_coordinate'][t_idx[0]:t_idx[-1]][mask_wd])
                df_dict_w['date_time'] = list(file['observations_table']['date_time'][t_idx[0]:t_idx[-1]][mask_wd])
                df_dict_w['date_time'] = seconds_to_datetime(df_dict_w['date_time'])
                df_dict_w['wd'] = list(file['observations_table']['observation_value'][t_idx[0]:t_idx[-1]][mask_wd])
                df_dict_w['ws'] = list(file['observations_table']['observation_value'][t_idx[0]:t_idx[-1]][mask_ws])

                # #####
                df_dict_w['fg_depar_ws'] = list(file['era5fb']['fg_depar@body'][t_idx[0]:t_idx[-1]][mask_ws])
                # #####

                df_dict_w['u'] = - np.abs(df_dict_w['ws']) * np.sin(np.radians(df_dict_w['wd']))
                df_dict_w['v'] = - np.abs(df_dict_w['ws']) * np.cos(np.radians(df_dict_w['wd']))
                
                # humidity data
                df_dict_h['z_coordinate'] = list(file['observations_table']['z_coordinate'][t_idx[0]:t_idx[-1]][mask_rh])
                df_dict_h['date_time'] = list(file['observations_table']['date_time'][t_idx[0]:t_idx[-1]][mask_rh])
                df_dict_h['date_time'] = seconds_to_datetime(df_dict_h['date_time'])
                df_dict_h['rh'] = list(file['observations_table']['observation_value'][t_idx[0]:t_idx[-1]][mask_rh])

                #temperature data
                df_dict['z_coordinate'] = list(file['observations_table']['z_coordinate'][t_idx[0]:t_idx[-1]][mask_t])
                df_dict['date_time'] = list(file['observations_table']['date_time'][t_idx[0]:t_idx[-1]][mask_t])
                df_dict['date_time'] = seconds_to_datetime(df_dict['date_time'])
                df_dict['t'] = list(file['observations_table']['observation_value'][t_idx[0]:t_idx[-1]][mask_t])

                # #####
                df_dict['fg_depar_t'] = list(file['era5fb']['fg_depar@body'][t_idx[0]:t_idx[-1]][mask_t])
                # #####

            else:
                rts = file['recordindices']['recordtimestamp'][:]
                idx = np.where(np.logical_and((rts >= dt_from), (rts <= dt_to)))[0]
                if len(idx) == 0:
                    print('NO DATA FOUND IN IGRA: ', sid)
                    return [rmse_sum_shbase_sonde, rmse_sum_shdisp_sonde,
                            rms_sum_shbase, rms_sum_sonde,
                            rms_sum_shdisp, rms_sum_dispminusbase]

                t_idx = file['recordindices']['126'][idx]
                p_mask = file['observations_table']['z_coordinate'][t_idx[0]:t_idx[-1]]
                mask_t = np.isin(p_mask,plevs)

                wd_idx = file['recordindices']['106'][idx]
                p_mask = file['observations_table']['z_coordinate'][wd_idx[0]:wd_idx[-1]]
                mask_wd = np.isin(p_mask,plevs)

                ws_idx = file['recordindices']['107'][idx]
                p_mask = file['observations_table']['z_coordinate'][ws_idx[0]:ws_idx[-1]]
                mask_ws = np.isin(p_mask,plevs)

                rh_idx = file['recordindices']['138'][idx]
                p_mask = file['observations_table']['z_coordinate'][rh_idx[0]:rh_idx[-1]]
                mask_rh = np.isin(p_mask,plevs)

                t_len = len(mask_t[mask_t == True])
                plevs = [1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000,92500,100000]

                # wind data
                df_dict_w['z_coordinate'] = list(file['observations_table']['z_coordinate'][wd_idx[0]:wd_idx[-1]][mask_wd])
                df_dict_w['date_time'] = list(file['observations_table']['date_time'][wd_idx[0]:wd_idx[-1]][mask_wd])
                df_dict_w['date_time'] = seconds_to_datetime(df_dict_w['date_time'])
                df_dict_w['wd'] = list(file['observations_table']['observation_value'][wd_idx[0]:wd_idx[-1]][mask_wd])
                df_dict_w['ws'] = list(file['observations_table']['observation_value'][ws_idx[0]:ws_idx[-1]][mask_ws])

                # #####
                df_dict_w['fg_depar_ws'] = list(file['era5fb']['fg_depar@body'][ws_idx[0]:ws_idx[-1]][mask_ws])
                # #####

                df_dict_w['u'] = - np.abs(df_dict_w['ws']) * np.sin(np.radians(df_dict_w['wd']))
                df_dict_w['v'] = - np.abs(df_dict_w['ws']) * np.cos(np.radians(df_dict_w['wd']))
                
                # humidity data
                df_dict_h['z_coordinate'] = list(file['observations_table']['z_coordinate'][rh_idx[0]:rh_idx[-1]][mask_rh])
                df_dict_h['date_time'] = list(file['observations_table']['date_time'][rh_idx[0]:rh_idx[-1]][mask_rh])
                df_dict_h['date_time'] = seconds_to_datetime(df_dict_h['date_time'])
                df_dict_h['rh'] = list(file['observations_table']['observation_value'][rh_idx[0]:rh_idx[-1]][mask_rh])

                #temperature data
                df_dict['z_coordinate'] = list(file['observations_table']['z_coordinate'][t_idx[0]:t_idx[-1]][mask_t])
                df_dict['date_time'] = list(file['observations_table']['date_time'][t_idx[0]:t_idx[-1]][mask_t])
                df_dict['date_time'] = seconds_to_datetime(df_dict['date_time'])
                df_dict['t'] = list(file['observations_table']['observation_value'][t_idx[0]:t_idx[-1]][mask_t])

                # #####
                df_dict['fg_depar_t'] = list(file['era5fb']['fg_depar@body'][t_idx[0]:t_idx[-1]][mask_t])
                # #####

                



            #meta data
            df_dict['latitude'] = list(file['observations_table']['latitude'][t_idx[0]:t_idx[-1]][mask_t])
            df_dict['longitude'] = list(file['observations_table']['longitude'][t_idx[0]:t_idx[-1]][mask_t])
            repid = np.asarray(file['observations_table']['report_id'][t_idx[0]:t_idx[-1]][mask_t])
            df_dict['report_id'] = list(repid.view('|S{}'.format(repid.shape[1])).flatten().astype(str))

            df_t = pd.DataFrame.from_dict(df_dict)
            df_w = pd.DataFrame.from_dict(df_dict_w)
            df_h = pd.DataFrame.from_dict(df_dict_h)

            #####
            if bg_depar_check:
                df = pd.merge(df_t, df_w[['z_coordinate', 'date_time', 'u', 'v', 'fg_depar_ws']], on=['z_coordinate', 'date_time'], how='inner')
            else:
                df = pd.merge(df_t, df_w[['z_coordinate', 'date_time', 'u', 'v']], on=['z_coordinate', 'date_time'], how='inner')
            #####

            df = pd.merge(df, df_h[['z_coordinate', 'date_time', 'rh']], on=['z_coordinate', 'date_time'], how='inner')
            df['q'] = rasotools.met.humidity.vap2sh(rasotools.met.humidity.rh2vap(df.rh, df.t), df.z_coordinate)
            
            #####
            if bg_depar_check:
                t_pc01 = {}
                t_pc99 = {}
                ws_pc01 = {}
                ws_pc99 = {}

                for i in stdplevs:
                    t_pc01[i] = np.nanpercentile(df[df.z_coordinate == i].fg_depar_t, 1)
                    t_pc99[i] = np.nanpercentile(df[df.z_coordinate == i].fg_depar_t, 99)
                    ws_pc01[i] = np.nanpercentile(df[df.z_coordinate == i].fg_depar_ws, 1)
                    ws_pc99[i] = np.nanpercentile(df[df.z_coordinate == i].fg_depar_ws, 99)
            #####

            lat_disp, lon_disp, sec_disp = np.array([np.nan]*len(df)),np.array([np.nan]*len(df)),np.array([np.nan]*len(df))
            
            for rid in df.report_id.drop_duplicates():
                df_j = df[df.report_id == rid].copy()

                #####
                if bg_depar_check:
                    # remove this if no bg threshold check should be done
                    filter_j = []
                    for i in df_j.z_coordinate:
                        i = int(i)
                        df_j_i = df_j[df_j.z_coordinate == i]
                        if df_j_i.fg_depar_ws.iloc[0] > ws_pc99[i] or df_j_i.fg_depar_ws.iloc[0] < ws_pc01[i]:
                            filter_j.append(False)
                        elif df_j_i.fg_depar_t.iloc[0] > t_pc99[i] or df_j_i.fg_depar_t.iloc[0] < t_pc01[i]:
                            filter_j.append(False)
                        else:
                            filter_j.append(True)
                    
                    # remove all the outliers
                    df_j = df_j[filter_j]        
                #####

                df_j_cleanded = df_j.sort_values(by='z_coordinate', ascending=False).dropna(subset=['t', 'u', 'v'])
                if len(df_j_cleanded) > 3:

                    idx =  df_j_cleanded.index.values
                    lat_i, lon_i, sec_i = trj.trajectory(df_j_cleanded.latitude.iloc[0], 
                                                            df_j_cleanded.longitude.iloc[0], 
                                                            df_j_cleanded.u.values, 
                                                            df_j_cleanded.v.values, 
                                                            df_j_cleanded.z_coordinate.values, 
                                                            df_j_cleanded.t.values
                                                        )
                    lat_disp[idx] = lat_i
                    lon_disp[idx] = lon_i
                    sec_disp[idx] = sec_i
            df['latitude_displacement'] = lat_disp
            df['longitude_displacement'] = lon_disp
            df['time_displacement'] = sec_disp

    except:
        print('NO DATA FOUND IN IGRA: ', sid)
        return [rmse_sum_shbase_sonde, rmse_sum_shdisp_sonde,
                rms_sum_shbase, rms_sum_sonde,
                rms_sum_shdisp, rms_sum_dispminusbase]
        
    df = df.dropna()
    if len(df) == 0:
        print('NO DATA FOUND IN IGRA: ', sid)
        return [rmse_sum_shbase_sonde, rmse_sum_shdisp_sonde,
                rms_sum_shbase, rms_sum_sonde,
                rms_sum_shdisp, rms_sum_dispminusbase]

    for mon in selected_mons:

        df_mon = df[df.date_time.dt.month == mon]
    #     display(df_mon)
        t0 = time.time()


        ds_fc = ds_fc_store[mon] 
        print('calculating month: ', mon,  time.time() - t0)
        t0 = time.time()

        t0 = time.time()     
        ##
        ## change to report id
        ##
        for day in df_mon.date_time.drop_duplicates()[:]:

            input_data = df_mon[df_mon.date_time == day]
            ds_fc_time = ds_fc.sel(time=day, method='nearest')
    #             print(day, pd.Timestamp(ds_fc_time.time.values))
            if (pd.Timestamp(ds_fc_time.time.values) - day) > maxtimediff:
                continue
            t_list = []
            for i in np.array(ds_fc_time.level): #10,20,...,1000
                step = find_nearest(input_data.z_coordinate, i*100)
                ###
                input_data_step = input_data[input_data.z_coordinate == step]
                station_lat = input_data.latitude.iloc[0] + np.array(input_data_step.latitude_displacement)[0]
                station_lon = input_data.longitude.iloc[0] + np.array(input_data_step.longitude_displacement)[0]
                lon = station_lon
                if lon < 0:
                    lon = 360.+lon
                ds_now = ds_fc_time.interp(latitude=[station_lat], longitude=[lon], method="linear")
                t = ds_now[varsel].sel(level = i)
                t_list.append(float(t))

            p_ml = [1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000,92500,100000]
            lon = input_data.longitude.iloc[0]
            if lon < 0:
                lon = 360.+lon
            base_t = np.array(ds_fc_time.interp(latitude=[input_data.latitude.iloc[0]], longitude=[lon], method="linear")[varsel])
            for i in range(len(stdplevs)):
                if np.abs(stdplevs[i] - find_nearest(input_data.z_coordinate,stdplevs[i])) > 500:
                    rmse_sum_shbase_sonde[stdplevs[i]].append(np.nan)
                    rmse_sum_shdisp_sonde[stdplevs[i]].append(np.nan)
                    rms_sum_shbase[stdplevs[i]].append(np.nan)
                    rms_sum_sonde[stdplevs[i]].append(np.nan)
                    rms_sum_shdisp[stdplevs[i]].append(np.nan)
                    rms_sum_dispminusbase[stdplevs[i]].append(np.nan)

                else:
                    sq_t = np.squeeze(base_t)
                    t_base = float(sq_t[p_ml == find_nearest(p_ml,stdplevs[i])])
                    t_disp = float(np.array(t_list)[p_ml == find_nearest(p_ml,stdplevs[i])])
                    input_data_step = input_data[input_data.z_coordinate == find_nearest(input_data.z_coordinate, stdplevs[i])]
                    try:
                        t_sonde = float(input_data_step[varsel])
                        rmse_sum_shbase_sonde[stdplevs[i]].append(t_base - t_sonde)
                        rmse_sum_shdisp_sonde[stdplevs[i]].append(t_disp - t_sonde)
                        rms_sum_shbase[stdplevs[i]].append(t_base)
                        rms_sum_sonde[stdplevs[i]].append(t_sonde)
                        rms_sum_shdisp[stdplevs[i]].append(t_disp)
                        rms_sum_dispminusbase[stdplevs[i]].append(t_disp-t_base)
                    except:
                        rmse_sum_shbase_sonde[stdplevs[i]].append(np.nan)
                        rmse_sum_shdisp_sonde[stdplevs[i]].append(np.nan)
                        rms_sum_shbase[stdplevs[i]].append(np.nan)
                        rms_sum_sonde[stdplevs[i]].append(np.nan)
                        rms_sum_shdisp[stdplevs[i]].append(np.nan)
                        rms_sum_dispminusbase[stdplevs[i]].append(np.nan)

                    

        print('calculating mon: ',mon, time.time() - t0)

    return [rmse_sum_shbase_sonde, rmse_sum_shdisp_sonde,
            rms_sum_shbase, rms_sum_sonde,
            rms_sum_shdisp, rms_sum_dispminusbase]
 
    

if __name__ == '__main__':
#     year = 1990
    save_dict = {'eastward windspeed':'u', 'northward windspeed':'v', 'air temperature':'temperature', 'specific humidity':'q'}
    units_dict = {'eastward windspeed':'m/s', 'northward windspeed':'m/s', 'air temperature':'K', 'specific humidity':'1'}
    stdplevs = [1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000,92500]
    diff = True
    show_date = False
    varsel_dict = {'eastward windspeed':'u', 'northward windspeed':'v', 'air temperature':'t', 'specific humidity': 'q'}
    era_dict = {'t':'130','u':'131','v':'132','q':'133'}
    
    
    for var in ['air temperature', 'eastward windspeed', 'northward windspeed', 'specific humidity']: # ['eastward windspeed', 'northward windspeed', 'air temperature', 'specific humidity']: # ['eastward windspeed', 'northward windspeed', 'air temperature', 'specific humidity']
        for year in [2000]: # [1960, 1970, 1980, 1990, 2000, 2010, 2020]:
            ds_fc_store = {}
            for mon in [1,2,3,4,5,6,7,8,9,10,11,12]:
                varsel = varsel_dict[var]
                varsel_era = era_dict[varsel]
                files = glob.glob('/mnt/scratch/scratch/leo/scratch/era5/gridded/era5fct.' + str(year) + str(mon).zfill(2) + '*.' + varsel_era + '.nc')[0]
                ds_fc_store[mon] = xr.load_dataset(files)
            ds_fc_id = ray.put(ds_fc_store) 

            ###

            file_list = []
            for i in glob.glob('/scratch/das/federico/COP2_HARVEST_JAN2023/igra2/*.nc')[:]:
                sid = i.split('-data')[-2][-5:]
                print(sid)
                # results = calc_station(sid,year,var, ds_fc_id)
                # break
                file_list.append(calc_station.remote(sid,year,var, ds_fc_id))
            results = ray.get(file_list)
            with open('/users/staff/uvoggenberger/CEUAS/CEUAS/public/trajectory/era5source_newbgcheck_era5_' + save_dict[var] + '_fc_'+str(year)+'_rmse_data.p', 'wb') as file:
                pickle.dump(results, file)

            ###

            # with open('/users/staff/uvoggenberger/scratch/displacement_data/igra/world/era5_source_bgcheck/era5source_bgcheck_era5_' + save_dict[var] + '_fc_'+str(year)+'_rmse_data.p', 'rb') as file: # 'qc_era5_temperature_fc_2000_rmse_data.p' 
            #     results = pickle.load(file)

            ###

            # rmse_sum_shbase_sonde, rmse_sum_shdisp_sonde, rms_sum_shbase, rms_sum_sonde, rms_sum_shdisp, rms_sum_dispminusbase = copy.deepcopy(results[0])
            # for i in results[1:]:
            #     for k in [1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000,92500]:
            #         rmse_sum_shbase_sonde[k]  = rmse_sum_shbase_sonde[k] + i[0][k]
            #         rmse_sum_shdisp_sonde[k] = rmse_sum_shdisp_sonde[k] + i[1][k]
            #         rms_sum_shbase[k] = rms_sum_shbase[k] + i[2][k]
            #         rms_sum_sonde[k] = rms_sum_sonde[k] + i[3][k]
            #         rms_sum_shdisp[k] =  rms_sum_shdisp[k] + i[4][k]
            #         rms_sum_dispminusbase[k] = rms_sum_dispminusbase[k] + i[5][k]

            # print('valid ascents: ', len(rms_sum_shdisp[50000]))
            # t0 = time.time()
            # rmse_shbase_sonde=[]
            # rmse_shdisp_sonde=[]

            # rms_shbase=[]
            # rms_sonde=[]
            # rms_shdisp=[]
            # rms_dispmbase=[]

            # for i in range(len(stdplevs)):
            #     print(i, ' Pa')
            #     testdata = np.array(rms_sum_shdisp[stdplevs[i]])
            #     print('valid ascents: ', len(testdata[~ np.isnan(testdata)]))
            #     print()

            #     rmse_shbase_sonde.append(np.sqrt(np.nanmean((np.array(rmse_sum_shbase_sonde[stdplevs[i]])**2))))
            #     if show_date:    
            #         print('rmse_shbase_sonde - plev: ', stdplevs[i], ' RMSE: ', rmse_shbase_sonde[-1])
            #     rmse_shdisp_sonde.append(np.sqrt(np.nanmean((np.array(rmse_sum_shdisp_sonde[stdplevs[i]])**2))))
            #     if show_date:    
            #         print('rmse_shdisp_sonde - plev: ', stdplevs[i], ' RMSE: ', rmse_shdisp_sonde[-1])
            #     rms_shbase.append(np.sqrt(np.nanmean((np.array(rms_sum_shbase[stdplevs[i]])**2))))
            #     if show_date:    
            #         print('rms_shbase - plev: ', stdplevs[i], ' RMS: ', rms_shbase[-1])
            #     rms_sonde.append(np.sqrt(np.nanmean((np.array(rms_sum_sonde[stdplevs[i]])**2))))
            #     if show_date:    
            #         print('rms_sonde - plev: ', stdplevs[i], ' RMS: ', rms_sonde[-1])
            #     rms_shdisp.append(np.sqrt(np.nanmean((np.array(rms_sum_shdisp[stdplevs[i]])**2))))
            #     if show_date:
            #         print('rms_shdisp - plev: ', stdplevs[i], ' RMS: ', rms_shdisp[-1])
            #     rms_dispmbase.append(np.sqrt(np.nanmean((np.array(rms_sum_dispminusbase[stdplevs[i]])**2))))
            #     if show_date:
            #         print('rms_dispmbase - plev: ', stdplevs[i], ' RMS: ', rms_shdisp[-1])


            # print('')

            # fig, ax = maplt.subplots(1, 2, gridspec_kw={'width_ratios': [4, 1]}, figsize = (15,10))
            # ax1 = ax[0]
            # ax2 = ax[1] 
            # if var != 'specific humidity':
            #     ax1.set_yscale('log')
            #     ax2.set_yscale('log')
            # ax2.sharey(ax1)
            # if var == 'specific humidity':
            #     ax1.plot(100000*np.array(rmse_shbase_sonde),stdplevs,color='orange', label=r'RMSE vertical $\times 10^{-5}$')
            #     ax1.plot(100000*np.array(rmse_shdisp_sonde),stdplevs, color='red', label=r'RMSE slanted $\times 10^{-5}$')
            # else:
            #     ax1.plot(np.array(rmse_shbase_sonde),stdplevs,color='orange', label='RMSE vertical')
            #     ax1.plot(np.array(rmse_shdisp_sonde),stdplevs, color='red', label='RMSE slanted')

            # ax1_4 = ax1.twiny()
            # ax1_4.axvline(x=0, color='black', alpha=0.8, ls='--', lw=0.5)
            # if var == 'specific humidity':
            #     if diff:
            #         ax1_4.plot(100000*(np.array(rmse_shbase_sonde)-np.array(rmse_shdisp_sonde)),stdplevs,color='purple', label=r'RMSE difference vertical - slanted $\times 10^{-5}$')
            #         ax1_4.plot(np.array(rms_dispmbase)*100000,stdplevs, color='green', alpha=0.3, ls='--', label=r'RMS vertical - slanted $\times 10^{-5}$')
            # else:
            #     if diff:
            #         ax1_4.plot(np.array(rmse_shbase_sonde)-np.array(rmse_shdisp_sonde),stdplevs,color='purple', label='RMSE difference slanted - vertical')
            #         ax1_4.plot(np.array(rms_dispmbase),stdplevs, color='green', alpha=0.3, ls='--', label='RMS slanted - vertical')

            # ax1_4.legend(loc='upper right', prop={'size':14})
            # ax1.set_ylim(ax1.get_ylim()[::-1])
            # ax1.set_ylabel('pressure (Pa)')
            # ax1.set_xlabel(var+' RMSE (' +str(units_dict[var]) + ')')
            # ax1.legend(loc='lower left', prop={'size':14})
            # ax1.grid()

            # value_nr = []
            # for i in rmse_sum_shbase_sonde:
            #     value_nr.append(len(np.asarray(rmse_sum_shbase_sonde[i])[~np.isnan(rmse_sum_shbase_sonde[i])]))
            # if var == 'specific humidity':
            #     ax2.barh(stdplevs, value_nr, 3000, color='g', alpha = 0.4, align='center')
            # else:
            #     ax2.barh(stdplevs, value_nr, np.array(stdplevs)/7, color='g', alpha = 0.4, align='center')
            # ax2.set_xlabel('Observations')
            # ax2.tick_params(labelleft=False)
            # ax2.grid()

            # #         maplt.title(str(year)+' Temperature RMSE \n' + str(len(results)) + ' stations    ' +str(len(rms_sum_shdisp[50000])) +' valid ascents')
            # maplt.title(str(year)+' '+var+' RMSE \n' +str(len(rms_sum_shdisp[50000])) +' valid ascents')
            # maplt.savefig(str(year)+'_'+save_dict[var]+'_era5_fc_world_rmse_plot_era5_log_bg_check.png')
            # maplt.close()
            # print('RMSE calculation: ', time.time()-t0)