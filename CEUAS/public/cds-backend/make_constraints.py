""" Make a pandas dataframe then a csv file to create the table of constraints for the CDS interface.
    OLD version, September 2020"""


import xarray as xr
import pandas as pd
import numpy as np
import os,sys
import path


from multiprocessing import Pool
from functools import partial
import datetime

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)

from tqdm import tqdm





#merged_files = [] 
"""
Eastward wind component	m s-1	Atmospheric eastward wind component at given pressure level
Geopotential height	m	Geopotential altitude from corrected pressure product
Northward wind component	m s-1	Atmospheric northward wind component at given pressure level
Relative humidity	%	Ratio of actual water vapor pressure to saturation water vapor pressure
Specific humidity	Kg/kg	Mass fraction of water vapor in atmospheric volume
Temperature	K	Atmospheric temperature at given pressure level
Wind from direction	Deg (0=N, 90=E)	Atmospheric wind direction at given pressure or height level, following the metorological convention, i.e. wind direction is wind from direction, it increases clockwise and starts with 0 deg if wind comes from North
Wind speed
"""

 
# 1) create a list of all the unique date_time, with a list of the files where this date is available 



def make_df_perstation(file):
    
    #std_plevs    = [1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 85000, 92500, 100000]
    
    station = file.split('_CEUAS_')[0].split('/')[-1]
    
    processed = open('processed_stations.txt' , 'a')
    empty       = open('empty_stations.txt' , 'a')
        
    #columns = [ 'date_time' , 'z_coordinate' , 'Eastward_wind_component' , 'Geopotential_height' , 'Northward_wind_component',  'Relative_humidity' , 'Specific_humidity' , 'Temperature' ,  'Wind_from_direction' , 'Wind_speed' ]     
    #df = pd.DataFrame(columns = columns)  # only available data 
    
    to_drop = ['adjustment_id', 'advanced_assimilation_feedback', 'advanced_homogenisation', 'advanced_qc', 'advanced_uncertainty', 
                      'bbox_max_latitude', 'bbox_max_longitude', 'bbox_min_latitude', 'bbox_min_longitude', 'code_table', 'conversion_flag', 'conversion_method', 
                      'crs', 'data_policy_licence', 'date_time_meaning', 'exposure_of_sensor', 'latitude', 'location_method', 'location_precision', 'longitude', 'numerical_precision', 
                      'observation_duration', 'observation_height_above_station_surface', 'observation_id', 'original_code_table', 'original_precision', 'original_units', 'original_value', 
                      'processing_level', 'quality_flag', 'report_id', 'secondary_value', 'secondary_variable', 'sensor_automation_status', 'source_id', 'spatial_representativeness', 
                      'traceability', 'units', 'value_significance', 'z_coordinate_method']

    #variables=[ 'date_time', 'z_coordinate', 'observed_variable' , 'observation_value',  'z_coordinate_type', ]   
    
    ot = xr.open_dataset(file, engine = 'h5netcdf' , group = "observations_table", decode_times = True, drop_variables= to_drop ).to_dataframe()  #(slow, requires lot of memory)
    # reducing the observations table to standard pressure levels (slow, requires lot of memory)
    # - standard pressure levels
    # - valid observation values
    # - z_coordinate_type as pressure
    # - variables: see dict. (dp=dew point, dpd = dew point depression)
    
    var = { 'temp': 85   , 'hrel': 38      , 'hspecific': 39 ,  'dp': 36 , 'dpd': 34 ,
                'uwind': 104 , 'vwind': 105 ,  'speed': 107  ,  'direction': 106      ,
                'gph'   : 117   ,  }
    
    
    
    red = ot.loc [(ot['z_coordinate_type'] == 1) &
                  
                         ( (ot['z_coordinate'] == 1000)   | (ot['z_coordinate'] == 2000)   | (ot['z_coordinate'] == 3000)   | (ot['z_coordinate'] == 5000) | (ot['z_coordinate'] == 7000)   |
                         (ot['z_coordinate'] == 10000) | (ot['z_coordinate'] == 15000) | (ot['z_coordinate'] == 20000) | (ot['z_coordinate'] == 25000) | (ot['z_coordinate'] == 30000) |
                         (ot['z_coordinate'] == 40000) | (ot['z_coordinate'] == 50000) | (ot['z_coordinate'] == 70000) | (ot['z_coordinate'] == 85000) | (ot['z_coordinate'] == 92500) | (ot['z_coordinate'] == 100000) ) &
                         
                         ( (ot['observed_variable'] == 85)  | (ot['observed_variable'] == 38)   | (ot['observed_variable'] == 39)  | (ot['observed_variable'] == 36)   | (ot['observed_variable'] == 34)   |
                         ( ot['observed_variable'] == 104) | (ot['observed_variable'] == 105) | (ot['observed_variable'] == 107) | (ot['observed_variable'] == 106)  |
                         ( ot['observed_variable'] == 117)  ) &
                         
                         (ot['observation_value'] != np.nan )  ]
    
    del ot 
    
    if not red.empty:
        red['year']    = red['date_time'].dt.year 
        red['month'] = red['date_time'].dt.month 
        red['day']     = red['date_time'].dt.day 
        
        red.to_pickle( out_dir + '/' + station +  '_df_pickled.pkl' )
        #processed.write(station + '\n' )
        
    else:
        print('skipping EMPTY' , file )
        red.to_pickle( out_dir + '/' + station +  '_df_pickled_EMPTY.pkl' )
        
        #empty.write(station + '\n' )
    print('done *** ' , file )
    return 0
    
    
    
    
out_dir = '/raid60/scratch/federico/correct_pickled_df/' 
os.system('mkdir ' + out_dir )

merged_directory = '/raid60/scratch/federico/DATABASE_V1_11SEPTEMBER2020' 
merged_files = os.listdir(merged_directory)







""" Checking already processed files """
def get_processed(out_dir):
    processed_stations =  [ s.split('_df_pickled')[0] for s in os.listdir(out_dir) if '62131' not in s ]
    return processed_stations
processed = get_processed(out_dir)


""" Creating list of missing files """
to_do = []
for s in merged_files:
    station = s.split('_CEUAS_')[0].split('/')[-1]

    if station not in processed:
        to_do.append(s)
to_do_files = [ merged_directory + '/' + s for s in to_do ]


pickle_files = [ f for f in os.listdir( '/raid60/scratch/federico/correct_pickled_df' ) if 'EMPTY' not in f ]

def make_pickles(merged_files):
    p      = Pool(25)
    func = partial(make_df_perstation)    
    out   = p.map(func, merged_files)
    0    

def merge_pickles(pickle_files, read_pickle = True ) :
    
    if read_pickle:
        first = pd.read_pickle( out_dir + '/' + pickle_files[0])
        first = first[ ['year' , 'month', 'day', 'z_coordinate' , 'observed_variable'] ]
        first = first.drop_duplicates( keep='first')
        partial = first 
    else:
        partial = pickle_files[0]
        
        
    #for p, num in zip( tqdm(pickle_files) ):
    for p,num  in  zip(pickle_files, range(len(pickle_files) ) ) :
    
        try:
            
            if read_pickle: # if read=True, I extract the df from the pickle file, otherwise I am already using a df 
                p = out_dir + '/' + p
                temp = pd.read_pickle( p )
            else:
                print('Merging the partial dataframes +++ ')
                temp = p 
            
            temp = temp[ ['year' , 'month', 'day', 'z_coordinate' , 'observed_variable' ] ]
            
            temp = temp.drop_duplicates( keep='first')
            
            a = pd.concat([partial, temp]) .drop_duplicates(keep='first')
            partial = a 
            print(' Length after concatenating: ' , len(partial) , ' ' , num  )
        except:
            pass
    
    partial = partial.sort_values( by=['year', 'month' , 'day' , 'z_coordinate' , 'observed_variable'] )
    #partial.to_pickle( 'df_constraints_pickled.pkl' )
    ### total length of th erecords should be 5372800 i.e. 8 variables * 16 plevels * 365 days * 115 years 
    print('DONE +++ ' )
    return partial 

def split_list(lista, length):
    """ Split the list 'lista' of files into chunks of maximum length 'length' """
    chunks = [lista[x:x+length] for x in range(0, len(lista), length)]
    return chunks
    

if __name__ == '__main__':
    
    #a = make_pickles(to_do_files)
    
    #pickle_list = [ pickle_files[0:500]       , pickle_files[1000:1500] ,  pickle_files[2000:2500] , pickle_files[3000:3500] , pickle_files[4000:], 
    #                      pickle_files[500:1000] , pickle_files[1500:2000] ,  pickle_files[2500:3000] , pickle_files[3500:4000]  ]
    
    
    pickle_list = split_list(pickle_files[:10], 11)
    
    pools          = len(pickle_list)
    p                = Pool( pools )
    func           = partial(merge_pickles)    
    partial_out = p.map(func, pickle_list)    
    output_df = merge_pickles(partial_out, read_pickle = False )
    output_df.to_csv('constraints.csv')
    
    #output_df.to_pickle('df_pickled_constraints_prova.pkl' )
    
    print(0)    
    

    


#for f in merged_files[:100]:
    #dummy = make_df_perstation( f )
 
#np.save('timestamps_per_station', timestamps_per_station, allow_pickle = True)

        




'''
def check_sensor_id(file):
    
    ot = xr.open_dataset(file, engine = 'h5netcdf' , group = "observations_table", decode_times = False )
    if 'sensor_id' not in ot.keys():
        print('MISSING SENSOR ID ::: ' , file , os.path.getsize(file )/1000000 )
    else:
        print('OK! ========================================= ' , file )

        
for f in merged_files:
    check_sensor_id(f)    
'''
