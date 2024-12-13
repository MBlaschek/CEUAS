import pandas as pd 
import numpy as np 
import h5py
import hdf5plugin
import glob
import ray 
import sys, os
sys.path.append(os.getcwd() + "/../cds-backend/code/")
import cds_eua4 as eua



def datetime_to_seconds(dates, ref='1900-01-01T00:00:00'):
    """ from datetime64 to seconds since 1900-01-01 00:00:00"""
    return ((dates - np.datetime64(ref)) / np.timedelta64(1, 's')).astype(np.int64)

"""
'wigos_statid', 'report_id', 'pressure', 'time', 'lat', 'lon',
       'station_id', 'airTemperature', 'airTempBiasCorr', 'sonde_type',
       'latitudeDisplacement', 'longitudeDisplacement', 'timePeriod',
       'source_id', 'height', 'platform_type', 'dewpointTemperature',
       'dewpointTempBiasCorr', 'windSpeed', 'windSpeedBiasCorr',
       'windDirection', 'windDirectionBiasCorr', 'geopotential',
       'nonCoordinateGeopotentialHeight',
       'extendedVerticalSoundingSignificance', 'lat05', 'lon05'],

       'insitu-comprehensive-upper-air-observation-network', {'variable': ['dew_point_temperature'], 
       'format': 'nc', 'dummy': 'avoid_cach_2024-04-29 13:24:58.463435', 
       'cdm': ['observations_table/source_id', 'station_configuration/platform_type'],
       'optional': ['humidity_bias_estimate', 'station_elevation'], 
       'date': ['19410803']}, '/ec/res4/scratch/erc/CUON_DISPLACEMENT_FINAL/CUON_dew_point_temperature_19410803.zip'

"""

biasadj = {'dew_point_temperature':['station_elevation', 'humidity_bias_estimate',],
            'air_temperature':['station_elevation', 'RISE_bias_estimate', 'sonde_type', 'latitude_displacement', 'longitude_displacement', 'time_since_launch'],
              'wind_speed':['station_elevation', 'wind_bias_estimate'],
                'wind_direction':['station_elevation', 'wind_bias_estimate'],
                  'geopotential':['station_elevation']
                  }


for year in [2008]: # range(2007, 2024): # , 1980, 2020, 2021, 2022, 2023]:
    try:
        os.mkdir('/mnt/users/scratch/uvoggenberger/to_bufr_1/'+str(year))
    except:
        pass
    for month in range(12,13): # range(1,2): # !!! check if really range(1,13)
        for day in range(1,32): # range(1,3): # !!! check if really range(1,32)
            for var in ['air_temperature','dew_point_temperature',  'wind_speed', 'wind_direction', 'geopotential']:
                dt = [str(year) + str(month).zfill(2) + str(day).zfill(2)]
                rq = {
                    # "statid": "11035",
                    "date": dt,
                    "variable": [var],
                    "optional": biasadj[var],
                    'cdm': ['observations_table/source_id', 'station_configuration/platform_type', "observations_table/observation_height_above_station_surface"],
                    "format": "nc",
                }
                if var == "wind_direction":
                    df_v23 = eua.vm_request_wrapper(rq, overwrite=True, request_filename = "/mnt/users/scratch/uvoggenberger/to_bufr_1/"+str(year)+"/CUON_"+"wind_from_direction"+"_"+str(dt[0])+".zip" ,vm_url="http://127.0.0.1:8007", download_only=True)
                else:
                    df_v23 = eua.vm_request_wrapper(rq, overwrite=True, request_filename = "/mnt/users/scratch/uvoggenberger/to_bufr_1/"+str(year)+"/CUON_"+var+"_"+str(dt[0])+".zip" ,vm_url="http://127.0.0.1:8007", download_only=True)
                    # df_v23 = eua.vm_request_wrapper(rq, overwrite=True, request_filename = "/mnt/users/scratch/uvoggenberger/to_bufr_1/CUON_"+var+"_"+str(dt[0])+".zip" ,vm_url="http://127.0.0.1:8007", download_only=True)
        #         break
        #     break
        # break


# data_dir = "/mnt/users/scratch/leo/scratch/converted_v23/long/"

# target_date = np.datetime64("2000-01-01")
# target_date = datetime_to_seconds(target_date) + (12*60*60)

# files = glob.glob(data_dir + "*10393*.nc")
# print(files)

# # @ray.remote
# def read_data(file, dt):
#     with h5py.File(file, 'r') as fl:
#         rts = fl['recordindices']['recordtimestamp'][:]

#         idx = np.searchsorted(rts, dt)
#         if (idx > 0) and (idx < (len(rts) - 1)):
#             result = np.min(rts[idx - 1:idx + 1])
#         elif idx == 0:
#             result = np.min(rts[idx:idx + 1])
#         else:
#             result = np.min(rts[idx - 1:idx])

#         if np.abs(result - dt) <= 43200: # half a day in seconds
#             timestamps = rts[np.logical_and(rts >= dt , rts < dt + 86400)]
#             vars = fl["recordindices"].keys()

#         else:
#             return None

# read_data(files[0], target_date)
# # ray.init(num_cpus=20)

# # ray.shutdown()