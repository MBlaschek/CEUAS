""" Create a single netCDF file contaning gridded data.
      structure of file:
      - dimensions
          lat
          lon
          time
          pleve
"""

import os,sys

import xarray as xr
import numpy as np
import pandas as pd


gridded_files_dir = '/raid60/scratch/federico/GRIDDED_FILES_FEB2021/ta'
files =  os.listdir(gridded_files_dir)

lat = list(set( [ float(f.split('_')[5]) for f in os.listdir(gridded_files_dir) ] ) ) 
lat.sort()
lon = list(set( [ float(f.split('_')[6]) for f in os.listdir(gridded_files_dir) ] ) ) 
lon.sort()

Lat = np.array(lat)
Lon = np.array(lon)
Plev = np.array([1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 85000, 92500, 100000])
Hour = np.array([0,12])


Time = xr.open_dataset( gridded_files_dir + '/' + files[0] , engine = 'h5netcdf', decode_times = True )['time'].values





""" Creating hour array """
df = xr.open_dataset( gridded_files_dir + '/' + files[0] , engine = 'h5netcdf', decode_times = True ).to_dataframe()
q = pd.arrays.DatetimeArray( df['time'][:] )
hours = q.hour
df['hour'] = hours

Time = df.loc [ (df['hour'] == 12) & (df['plev']==100000) ]['time'].values  # just need one plev per  hour, i.e. this is the list of distinct  time stamps

""" Smaller example """
Lat = Lat[:3]
Lon = Lon[:5]
res = np.empty([len(Lat) , len(Lon), len(Hour) , len(Plev),  len(Time) ] )  # 2*16 is 2 hours x pressure levels 


"""
Lat = np.array(lat)
Lon = np.array(lon)
Plev = np.array([1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 85000, 92500, 100000])
Hour = np.array([0,12])


Time = xr.open_dataset( gridded_files_dir + '/' + files[0] , engine = 'h5netcdf', decode_times = True )['time'].values

"""
Lat = Lat[:3]
Lon = Lon[:5]

for lat in range(len(Lat)):
    for lon in range(len(Lon)):
        box_file = [f for f in files if str(Lat[lat])+'_'+str(Lon[lon]) in f ][0]
        #df = xr.open_dataset( gridded_files_dir + '/' + box_file , engine = 'h5netcdf', decode_times = True ).to_dataframe()
        #q = pd.arrays.DatetimeArray( df['time'][:] )
        #hh = q.hour
        #df['hour'] = hours  # it is the same for each dataframe !!!        
        for p in range(len(Plev)):
            a = 0 # read here input data 
            for h in range(len(Hour)):
                print(lat,lon,p,h)
                

                df_red = df.loc[ (df['plev'] == Plev[p] )   & ( df['hour'] == Hour[h] )]
                temp = df_red['ta_average'].values
                #temp_b= df_red['ta_average_bias'].values
                #ano = df_red['ta_anomaly'].values
                #ano_b = df_red['ta_anomaly_bias'].values
                
                res[lat,lon,h,p,:] = temp
                0

da = xr.DataArray (data = res, 
                                dims = ["lat","lon","hour","pressure","time"],
                                coords = dict( lat = Lat ,
                                                         lon =  Lon,
                                                         hour = Hour,
                                                         pressure = Plev,
                                                         time = Time,                  
                                                         ),
                                   
                                   attrs = dict ( 
                                                 title = 'CEUAS gridded data for temperatures and anomalies',
                                                 institution = 'University of Vienna',
                                                 source = 'Institut fuer Meteorologie und Geophysik, leopold.haimberger@univie.ac.at',
                                                 history = '2021-02-19 09:54:0') 
                                   )
                                                         
                                                 
                                           
dummy = da.to_netcdf('Try_temp_gridded.nc' , mode = 'w')                          
0
