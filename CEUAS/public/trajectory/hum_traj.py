#!/usr/bin/env
# coding: utf-8

import numpy as np
import pandas as pd
import sys, os, glob
import h5py

sys.path.append(os.getcwd()+'/../cds-backend/code/')
import cds_eua4 as eua

import ray
ray.init(num_cpus=20)

def seconds_to_datetime(seconds, ref='1900-01-01'):
    """ from seconds to datetime64 """
    seconds = np.asarray(seconds)
    return pd.to_datetime(seconds, unit='s', origin=ref)

@ray.remote
def do_station(sid):
    print(sid)
    
    if len(glob.glob('/users/staff/a1400070/scratch/humidity_trajectories_new/humidity_trajectory_out_'+str(sid)+'.nc')) != 0:
        return 0
    if int(sid) in list(pd.read_csv('bad_requests.csv').bad_request):
        return 0
    
#     sys.path.append(os.getcwd()+'/../cds-backend/code/')
#     import cds_eua4 as eua
#     try:
#         request={
#             'variable': ['air_temperature', 'relative_humidity'],
#             'optional':['RASE_bias_estimate', 'humidity_bias_estimate', 'latitude_displacement', 'longitude_displacement', 'time_since_launch'],
#             'statid': sid,
#             'pressure_level':[1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,35000,40000,45000,50000,55000,60000,65000,70000,75000,80000,85000,90000, 92500, 95000, 100000],
#             'date': ['19940101', '20230102'],
#             'format':'csv',
#         } 
#         df = eua.vm_request_wrapper(request, overwrite=True, vm_url='http://127.0.0.1:8007')
#     except:
#         with open('bad_requests.csv','a') as fd:
#             fd.write(sid+'\n')
#         return '1 - '+sid

    df_dict = {}
    conv_file = glob.glob('/mnt/users/scratch/leo/scratch/converted_v9/*'+ str(sid)+'*.nc')[0]
    try:
        with h5py.File(conv_file, 'r') as file:
            rts = file['recordindices']['recordtimestamp'][:]
            idx = np.where(rts >= 2966371200)[0]
            if len(idx) == 0:
                return 1
            t_idx = file['recordindices']['126'][idx]
            h_idx = file['recordindices']['138'][idx]

            plevs = [1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,35000,40000,45000,50000,55000,60000,65000,70000,75000,80000,85000,90000, 92500, 95000, 100000]


            mask = file['observations_table']['z_coordinate'][t_idx[0]:t_idx[-1]]
            mask = np.isin(mask,plevs)
            t_len = len(mask[mask == True])

            df_dict['z_coordinate'] = list(file['observations_table']['z_coordinate'][t_idx[0]:t_idx[-1]][mask])
            df_dict['date_time'] = list(file['observations_table']['date_time'][t_idx[0]:t_idx[-1]][mask])
            df_dict['observation_value'] = list(file['observations_table']['observation_value'][t_idx[0]:t_idx[-1]][mask])
            df_dict['latitude'] = list(file['observations_table']['latitude'][t_idx[0]:t_idx[-1]][mask])
            df_dict['longitude'] = list(file['observations_table']['longitude'][t_idx[0]:t_idx[-1]][mask])
            repid = np.asarray(file['observations_table']['report_id'][t_idx[0]:t_idx[-1]][mask])
            df_dict['report_id'] = list(repid.view('|S{}'.format(repid.shape[1])).flatten().astype(str))

            df_dict['RASE_bias_estimate'] = list(file['advanced_homogenisation']['RASE_bias_estimate'][t_idx[0]:t_idx[-1]][mask])
            df_dict['latitude_displacement'] = list(file['advanced_homogenisation']['latitude_displacement'][t_idx[0]:t_idx[-1]][mask])
            df_dict['longitude_displacement'] = list(file['advanced_homogenisation']['longitude_displacement'][t_idx[0]:t_idx[-1]][mask])
            df_dict['time_since_launch'] = list(file['advanced_homogenisation']['time_since_launch'][t_idx[0]:t_idx[-1]][mask])

            df_dict['variable'] = ['air_temperature']*t_len
            df_dict['humidity_bias_estimate'] = [np.nan]*t_len



            h_mask = file['observations_table']['z_coordinate'][h_idx[0]:h_idx[-1]]
            h_mask = np.isin(h_mask,plevs)
            h_len = len(h_mask[h_mask == True])

            df_dict['z_coordinate'].extend(list(file['observations_table']['z_coordinate'][h_idx[0]:h_idx[-1]][h_mask]))
            df_dict['date_time'].extend(list(file['observations_table']['date_time'][h_idx[0]:h_idx[-1]][h_mask]))
            df_dict['observation_value'].extend(list(file['observations_table']['observation_value'][h_idx[0]:h_idx[-1]][h_mask]))
            df_dict['latitude'].extend(list(file['observations_table']['latitude'][h_idx[0]:h_idx[-1]][h_mask]))
            df_dict['longitude'].extend(list(file['observations_table']['longitude'][h_idx[0]:h_idx[-1]][h_mask]))
            repid = np.asarray(file['observations_table']['report_id'][h_idx[0]:h_idx[-1]][h_mask])
            df_dict['report_id'].extend(list(repid.view('|S{}'.format(repid.shape[1])).flatten().astype(str)))

            df_dict['humidity_bias_estimate'].extend(list(file['advanced_homogenisation']['humidity_bias_estimate'][h_idx[0]:h_idx[-1]][h_mask]))

            df_dict['variable'].extend(['relative_humidity']*h_len)
            df_dict['RASE_bias_estimate'].extend([np.nan]*h_len)
            df_dict['latitude_displacement'].extend([np.nan]*h_len)
            df_dict['longitude_displacement'].extend([np.nan]*h_len)
            df_dict['time_since_launch'].extend([np.nan]*h_len)    

            df_dict['date_time'] = seconds_to_datetime(df_dict['date_time'])

            df = pd.DataFrame.from_dict(df_dict)

    except:
        with open('bad_requests.csv','a') as fd:
            fd.write(sid+'\n')
        return '1 - '+sid
    
    if len(df) == 0:
        return '1 - '+sid
    
    try:
        df['z_coordinate'] = np.log(df['z_coordinate'])
        df.date_time = pd.to_datetime(df.date_time)
        df.report_id = df.report_id.astype('int')
        print('data loaded')

        target_pressure = np.log([1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,35000,40000,45000,50000,55000,60000,65000,70000,75000,80000,85000,90000, 92500, 95000, 100000])


        out = {}
        for i in ['humidity_bias_estimate','relative_humidity', 'RASE_bias_estimate', 'air_temperature', 'date_time', 'latitude', 'longitude', 
                  'report_id', 'z_coordinate', 'latitude_displacement', 'longitude_displacement', 'seconds_since_release']: # 'primary_id'
            out[i] = []
        for i in df.report_id.drop_duplicates():
            try:
#             print('report_id', i)
                i_df = df[df.report_id == i]
                hum_df = i_df[i_df.variable == 'relative_humidity']
                temp_df = i_df[i_df.variable == 'air_temperature']
                if (np.abs(len(hum_df) - len(temp_df)) > 3) or (len(hum_df) == 0) or (len(temp_df) == 0):
                    continue

                out['latitude_displacement'].extend(list(np.interp(target_pressure, temp_df['z_coordinate'], temp_df['latitude_displacement'])))
                out['longitude_displacement'].extend(list(np.interp(target_pressure, temp_df['z_coordinate'], temp_df['longitude_displacement'])))
                out['seconds_since_release'].extend(list(np.interp(target_pressure, temp_df['z_coordinate'], temp_df['time_since_launch'])))

                out['humidity_bias_estimate'].extend(list(np.interp(target_pressure, hum_df['z_coordinate'], hum_df['humidity_bias_estimate'])))
                out['relative_humidity'].extend(list(np.interp(target_pressure, hum_df['z_coordinate'], hum_df['observation_value'])))
                out['RASE_bias_estimate'].extend(list(np.interp(target_pressure, temp_df['z_coordinate'], temp_df['RASE_bias_estimate'])))
                out['air_temperature'].extend(list(np.interp(target_pressure, temp_df['z_coordinate'], temp_df['observation_value'])))
                out['date_time'].extend(list([i_df['date_time'].iloc[0]]*len(target_pressure)))
                out['latitude'].extend(list([i_df['latitude'].iloc[0]]*len(target_pressure)))
                out['longitude'].extend(list([i_df['longitude'].iloc[0]]*len(target_pressure)))
                out['report_id'].extend(list([i_df['report_id'].iloc[0]]*len(target_pressure)))
        #             out['primary_id'].extend(list([i_df['primary_id'].iloc[0]]*len(target_pressure)))
                out['z_coordinate'].extend([1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,35000,40000,45000,50000,55000,60000,65000,70000,75000,80000,85000,90000, 92500, 95000, 100000])
            except:
                pass
        out_df = pd.DataFrame.from_dict(out)
        out_df
        print('everything calculated')
        ds = out_df.set_index('date_time').to_xarray()
        ds.to_netcdf('/users/staff/a1400070/scratch/humidity_trajectories_new/humidity_trajectory_out_'+str(sid)+'.nc')
        print('done writing to file')
        try:
            with h5py.File('/users/staff/a1400070/scratch/humidity_trajectories_new/humidity_trajectory_out_'+str(sid)+'.nc', 'r+') as file:
                file['date_time'].attrs['units'] = np.string_('nano seconds since 1970')
                file['date_time'].attrs['description'] = np.string_('release time of the sonde')
                file['air_temperature'].attrs['units'] = np.string_('Kelvin')
                file['air_temperature'].attrs['description'] = np.string_('observed air temperature')
                file['RASE_bias_estimate'].attrs['units'] = np.string_('Kelvin')
                file['RASE_bias_estimate'].attrs['description'] = np.string_('bias estimates: air_temperature - RASE_bias_estimate = adjusted temperature')
                file['relative_humidity'].attrs['units'] = np.string_('1')
                file['relative_humidity'].attrs['description'] = np.string_('observed relative humidity')
                file['humidity_bias_estimate'].attrs['units'] = np.string_('1')
                file['humidity_bias_estimate'].attrs['description'] = np.string_('bias estimates: relative_humidity - humidity_bias_estimate = adjusted relative humidity')
                file['latitude'].attrs['units'] = np.string_('degree north')
                file['latitude'].attrs['description'] = np.string_('station position')
                file['longitude'].attrs['units'] = np.string_('degree east')
                file['longitude'].attrs['description'] = np.string_('station position')
                file['latitude_displacement'].attrs['units'] = np.string_('degree north')
                file['latitude_displacement'].attrs['description'] = np.string_('latitude - latitude_displacement = actual position latitude')
                file['longitude_displacement'].attrs['units'] = np.string_('degree east')
                file['longitude_displacement'].attrs['description'] = np.string_('longitude - longitude_displacement = actual position ')
                file['seconds_since_release'].attrs['units'] = np.string_('seconds')
                file['seconds_since_release'].attrs['description'] = np.string_('date_time + seconds_since_release = actual time')
                file['z_coordinate'].attrs['units'] = np.string_('Pascal')
                file['z_coordinate'].attrs['description'] = np.string_('height coordinate')
    #             file['primary_id'].attrs['description'] = np.string_('station identifier')
                file['report_id'].attrs['description'] = np.string_('ascent identifier')
                file.attrs['file'] = np.string_('Humidity Bias Estimates on Trajectories')
                file.attrs['info'] = np.string_('Created by IMGW - University of Vienna Service Version 0, ' + datetime.datetime.now().strftime("%d-%b-%Y %H:%M:%S"))
        except:
            print('writing_error: ',i)
            with open('failed_requests.csv','a') as fd:
                fd.write(sid+'\n')

        print('done with attrs')

        return '0 - '+sid
    except:
        with open('failed_requests.csv','a') as fd:
            fd.write(sid+'\n')
        return '2 - '+sid
    
if __name__ == '__main__':
    with open('bad_requests.csv','w') as fd:
        fd.write('bad_request'+'\n')
    with open('failed_requests.csv','a') as fd:
        fd.write('failed_requests'+'\n')
    file_list = []
    for i in glob.glob('/mnt/users/scratch/leo/scratch/converted_v9/*2000*.nc')[:]:
        print(i)
        sid = i.split('_CEUAS_merged_v1.nc')[0][-5:]
#         print(sid)
#         print(do_station(sid))
        file_list.append(do_station.remote(sid))
    results = ray.get(file_list)
    print(results)