""" Make a pandas dataframe then a csv file to create the table of constraints for the CDS interface """


import xarray as xr
import pandas as pd
import numpy as np
import os,sys
import path
import h5py

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
    
    try:
        
        #std_plevs    = [1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 85000, 92500, 100000]
        
        station = file.split('_CEUAS_')[0].split('/')[-1]
        
        processed = open('processed_stations.txt' , 'a')
        empty       = open('empty_stations.txt' , 'a')
        
        """
        to_drop = ['adjustment_id', 'advanced_assimilation_feedback', 'advanced_homogenisation', 'advanced_qc', 'advanced_uncertainty', 
                          'bbox_max_latitude', 'bbox_max_longitude', 'bbox_min_latitude', 'bbox_min_longitude', 'code_table', 'conversion_flag', 'conversion_method', 
                          'crs', 'data_policy_licence', 'date_time_meaning', 'exposure_of_sensor', 'latitude', 'location_method', 'location_precision', 'longitude', 'numerical_precision', 
                          'observation_duration', 'observation_height_above_station_surface', 'observation_id', 'original_code_table', 'original_precision', 'original_units', 'original_value', 
                          'processing_level', 'quality_flag', 'report_id', 'secondary_value', 'secondary_variable', 'sensor_automation_status', 'source_id', 'spatial_representativeness', 
                          'traceability', 'units', 'value_significance', 'z_coordinate_method']
        """
        
        # reading and converting date_time with xarray 
        ot = xr.open_dataset(file, engine = 'h5netcdf' , group = "observations_table", decode_times = True ) #(slow, requires lot of memory)
        dt = ot['date_time'].values
    
        
        ot.close()
        
        # reading the data with h5py
        f =  h5py.File(file, 'r')
        
        obs_vars = ['observed_variable', 'z_coordinate' ]
        
        # get std plevels
        # z_coord = np.array(f["observations_table"]['z_coordinate'][:])
        obs_val = f["observations_table"]['observation_value'][:]
    
        # get lat and lon
        lat, lon =  f["observations_table"]['latitude'][0] , f["observations_table"]['longitude'][0]
        
        
        #plevel_filter = np.where( (z_coord == 1000) | (z_coord == 2000) | ( z_coord == 3000) |  (z_coord == 5000 )|  (z_coord == 7000) | (z_coord == 10000) | 
        #                                              (z_coord == 15000) |  (z_coord == 20000) |  (z_coord == 25000) |  (z_coord == 30000) |  (z_coord == 40000) |  (z_coord == 50000) | 
        #                                              (z_coord == 70000) |  (z_coord == 85000) |  (z_coord == 92500) |  (z_coord == 100000) ) 
        
        
        nan_filter = np.where (obs_val != np.nan )
        
        # applying nan mask and plevels mask 
        #total_mask = plevel_filter[nan_filter]
        
        #total_mask = np.intersect1d(plevel_filter, nan_filter)
        #total_mask =total_mask[0]
        
        #total_mask = plevel_filter & nan_filter  # total dat to be extracted 
    
        total_mask = nan_filter
        
        res = {}
    
        for c in obs_vars:
            res[c] = f['observations_table'][c][total_mask]
        
        res['date_time'] = dt[total_mask]
    
        red = pd.DataFrame.from_dict(res)
    
        del res 
        
        if not red.empty:
            red['year']  = red['date_time'].dt.year 
            red['month'] = red['date_time'].dt.month 
            red['day']   = red['date_time'].dt.day 
            
            #red = red.drop(columns ='date_time')
            #red = red.drop_duplicates()
            
            #red.to_pickle( out_dir + '/' + station +  '_df_pickled.pkl' )
            #red.to_csv( out_dir + '/' + station + '.csv', sep = '\t' )
            
            ydm_df = red[['year', 'month' , 'day', 'observed_variable']]
            
            for var in np.unique( ydm_df['observed_variable']):
                ind = np.where ( ydm_df['observed_variable'] == var )
                ydm_df_var = ydm_df.iloc[ ind ]
                ydm_df_var = ydm_df_var.drop_duplicates()
                
                try:
                    ydm_df_var.to_csv( out_dir + '/' + str(var) + '/' + str(var) + '_' + str(lat) + '_' + str(lon) + '_' + station + '.csv', sep = '\t' , index = False)
                except:
                    os.mkdir(out_dir + '/' + str(var) )
                    ydm_df_var.to_csv( out_dir + '/' + str(var) + '/' + str(var) + '_' + str(lat) + '_' + str(lon) + '_' + station + '.csv', sep = '\t' , index = False )
                    
            processed.write(file + '\n' )
            
        else:
            print('skipping EMPTY' , file )
    
            empty.write(station +'\n')
        print('done *** ' , file )
        return 0
    
    except:
        pass

#a = make_constraints(merged_files[0])


def make_constraints(merged_files):
    p      = Pool(25)
    func = partial(make_df_perstation)    
    out   = p.map(func, merged_files)
    0    





""" Main run """


out_dir = '/raid60/scratch/federico/PROVAAAAAA/' 

merged_directory = '/raid60/scratch/leo/scratch/converted_v5/' 
merged_files = [ merged_directory + f for f in os.listdir(merged_directory) if '.nc' in f ][:30]


make_pickles = False 
single = False


    
    
if make_pickles:
    
    os.system('mkdir ' + out_dir )

    if single:   
        for f in merged_files:
            make_df_perstation(f)

    else:
        dummy = make_constraints(merged_files)
    
else:
    
    out_dir_2 = out_dir + '/by_date/' 
    os.system('mkdir ' + out_dir_2 )
    
    for var in [v for v in os.listdir( out_dir) if v != '0' and 'by' not in v ]:   
        
        out_dir_2_var = out_dir_2 + '/' + var
        os.system('mkdir ' + out_dir_2_var )
        print("*** Doing variable: " , var )
        files = os.listdir(out_dir + '/' + var)
                    
        res = {} # keeping all results 
        
        for f in tqdm(files):
            lat, lon, station = f.split('_')[1] , f.split('_')[2] , f.split('_')[3].replace('.csv','')
            file = out_dir + '/' + var + '/' + f
            df = pd.read_csv(file, sep = '\t')
            
            years, months, days = df['year'] , df['month'], df['day']
            for y,m,d in zip( years, months, days ):
                
                dic = var + '_' + str(y) + '_' + str(m) + '_' + str(d)
                
                try:
                    res[dic]['lat'].append(lat)
                    res[dic]['lon'].append(lon)
                    res[dic]['station'].append(station)
                    
                except:
                
                    res[dic] = {}
                    

                    res[dic]['station'] = []
                    res[dic]['lon'] = [] 
                    res[dic]['lat'] = [] 
                    
                    res[dic]['lat'].append(lat)
                    res[dic]['lon'].append(lon)
                    res[dic]['station'].append(station)                    
                    

                #out_file = var + '_' + y + '_' + m + '_' + d + '.csv'
                #file = out_dir_2 + out_file
                #out = open(file, 'a+' )
                #out.write(lat + ',' + lon + ',' + station + '\n')
                
        for key in res.keys():
            df = pd.DataFrame.from_dict(res[key])
            out_file = out_dir_2_var  + '/'  +  key + '.csv'
            df.to_csv(out_file, sep = ',', index=False)





