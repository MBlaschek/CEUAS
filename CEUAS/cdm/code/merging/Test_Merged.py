""" Utility to extract basic information from the netCDF file of the merged database """


import os,sys
import netCDF4 as nc
import pandas as pd
import numpy as np
import random 
from datetime import datetime, timedelta
import xarray as xr 


print('Running the merged file checker **** ')

""" Loading """
a = '/raid8/srvx1/federico/GitHub/DEVELOP_LATEST_NOVEMBER/CEUAS/CEUAS/cdm/code/merging/FAST_INDEX_merged_chuadb_windc_82930.txt.nc' 


''' Print all '''
pd.set_option('display.max_rows', 50)
pd.set_option('display.max_columns', 20)
pd.set_option('display.width', 300)


""" Code mapping """
observation_ids_merged = {  'igra2':1 , 'ncar':2 , 'bufr':3,  'era5_1':4 , 'era5_1759' :5 , 'era5_1761':6 ,  'era5_3188' :7}  # values used to convert original record_id to the merged record_id, see merge_all_data 



def make_datetime(rts , ri):
    """ Make date_time """
    
    time_offset            = rts.units
    time_offset_value  = time_offset.split('since') [1]                            
    try:
        time_offset_value  = datetime.strptime(time_offset_value, '%Y-%m-%d %H:%M:%S')
    except:
        time_offset_value = time_offset_value.replace(' 19', '19')
        time_offset_value  = datetime.strptime(time_offset_value, '%Y-%m-%d %H:%M:%S')
        
        
    if 'minutes' in  time_offset:
        date_time_delta = [ timedelta(minutes = float(i) ) + time_offset_value for i in rts ]
    elif 'hours' in time_offset:
        date_time_delta = [ timedelta(hours = float(i) )  + time_offset_value  for i in rts ]    
    
    date_time_delta, ri = list(date_time_delta), list(ri) 
    
    dictionary = dict(zip(date_time_delta, ri)) 
    return dictionary  
    
    
            

File = nc.Dataset(a)
variables = list(File.variables )
groups     = list(File.groups )
                 
print('The file contains the following: ')
print('--- Variables: ' , variables )
print('--- Groups: ' , groups )

ri = File.variables['recordindex']
ri_values = list( File.variables['recordindex'][:] )

rts = File.variables['recordtimestamps']
rts_values=  list( File.variables['recordtimestamps'][:] )

#print('**** recordindex ****')
#print('Length: ' , len(ri))
#print(ri[:100])

#print('**** recordtimestamps ****')
#print('Length: ' , len(rts))
#print(rts[:100])



print(' I convert to proper date_time')
rts_datetime = make_datetime(rts , ri_values)

File.close()

""" Picking one random  time_stamp """
random = random.choice( list(rts_datetime.keys() ) )


""" Reading the tables """
obs_t =  xr.open_dataset( a , engine = 'h5netcdf' , group = 'observations_table' )
obs_df = obs_t.to_dataframe()   

head_t =  xr.open_dataset( a , engine = 'h5netcdf' , group = 'header_table' )
head_t = head_t.to_dataframe()   
head_df = head_t[ ['report_id', 'duplicates'] ] 

source_conf =  xr.open_dataset( a , engine = 'h5netcdf' , group = 'era5_1_source_configuration' )
print(' The source_configuration for era5_1 is : ' , source_conf )   
station_conf =  xr.open_dataset( a , engine = 'h5netcdf' , group = 'era5_1_station_configuration' ) 
print(' The station_configuration for era5_1 is : ' , station_conf )   




""" Loading a subset of variables to be put into a pd dataframe """
obs_df = obs_df [ ['date_time', 'z_coordinate' , 'z_coordinate_type', 'observed_variable' , 'observation_value' , 'report_id' , 'observation_id' , 'latitude' , 'longitude' , 'source_id', 'advanced_assimilation_feedback' , 'units']]
index_low , index_up = rts_datetime[random] ,   ri_values [    ri_values.index(rts_datetime[random]) + 1  ] 

print('\n\n\n I will analyze one random time_Stamp from the merged file:  ' ,   random , '  with time_record index: ' ,   rts_datetime[random])
print('\n\n\n The indices corresponding to the DF are: ' ,  index_low , index_up   )
print('\n\n\n I print one record below and one above the selected [index_low , index_up] range to check the correct time_stamp'   )




df = obs_df.iloc [index_low - 1 : index_up + 1 ]
head_df = head_df.iloc [index_low - 1 : index_up + 1 ]



print( '\n\n ***** Observations_table: ' ,  df )
print( '\n\n ***** Header_table: ' ,  head_df )




a = 0 