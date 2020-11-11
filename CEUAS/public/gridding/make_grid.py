""" Module to create a grid in lat-long boxes.
      1. Read the station_configuration from all the merged files
      2. Create a CSV (\t separated) file with all the station configurations
      3. Create a grid in lat-lon coordinates 
      4. Extract a list of stations included in each box  """



import pandas as pd
import xarray as xr
import os,sys
from tqdm import tqdm
import numpy as np

from multiprocessing import Pool
from functools  import partial

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)


merged_database = '/raid60/scratch/federico/MERGED_DATABASE_OCTOBER2020_sensor/'
station_configuration = []

file_list = os.listdir(merged_database)

def get_sc(file):
    ''' Extract the station_configuration form the merged file '''
    f = merged_database + '/' + file
    sc = xr.open_dataset (f, engine='h5netcdf', group = 'station_configuration', decode_times = True ).to_dataframe()    
    first_sc = sc[0:1]
    first_sc['file']= f 
    print('Done ***' , f )
    return first_sc
    

if not os.path.isfile('all_merged_station_configurations.csv'):
    print('*** Will create the station_configuration file out of the merged files ***')
    p = Pool(20)
    func = partial(get_sc)
    out = p.map(func, file_list)

    concatenated = pd.concat(out)
    os.system('rm all_merged_station_configurations.csv')
    concatenated.to_csv('all_merged_station_configurations.csv' , sep = '\t')
else:
    print('*** Will read the station_configuration file ***')
    sc = pd.read_csv('all_merged_station_configurations.csv', delimiter = '\t')
    
    """ Fixing wrong convention for values > 180 degrees """
    lat, lon = sc['latitude'].values , sc['longitude'].values
    lat_fix, lon_fix = [ -360 + l if l > 180 else l for l in lat ] ,  [ -360 + lg if lg > 180 else lg for lg in lon ] 
    
    sc['latitude']   = lat_fix
    sc['longitude'] = lon_fix
    
    print(0)
print('*** DONE extracting all the station_configuration tables ')

    
"""
for f in tqdm( os.listdir(merged_database) ):
    #print(f)
    f = merged_database + '/' + f
    sc = xr.open_dataset (f, engine='h5netcdf', group = 'station_configuration', decode_times = True ).to_dataframe()
    
    first_sc = sc[0:1]
    first_sc['file']= f 
    print(first_sc) 
    print('Done with ', f  )
    station_configuration.append(first_sc)

concatenated = pd.concat(station_configuration)
concatenated.to_csv('all_merged_station_configurations.csv' , sep = '\t')
"""



def make_grid(size=10):
    """ Divide the globe into squared boxes of the given size (default 10 degrees) """
    
    boxes = {}
    for i, n in zip ( range(-180,180, size) , range(1,100000) ):     # longitude
        for j,m in zip( range(-180,180, size) , range (1,100000) ): # latitude 

            center_lat = i + size/2
            center_lon = j + size/2
            
            boxes[str(n)+'_'+str(m)] = {}
            boxes[str(n)+'_'+str(m)]['lat'] = []
            boxes[str(n)+'_'+str(m)]['lon'] = []            
            
            lat_min , lon_min  = center_lat - size/2 ,  center_lon - size/2
            lat_max , lon_max = center_lat + size/2  ,  center_lon + size/2 
            
            boxes[str(n)+'_'+str(m)]['lat'].append( float( "{:.2f}".format(lat_min) ) ) 
            boxes[str(n)+'_'+str(m)]['lat'].append( float( "{:.2f}".format(lat_max) ) ) 
            boxes[str(n)+'_'+str(m)]['lon'].append(float( "{:.2f}".format (lon_min) ) )
            boxes[str(n)+'_'+str(m)]['lon'].append(float( "{:.2f}".format (lon_max) ) )
    
    return boxes


""" Creating the gridded boxes """
boxes = make_grid(size=10)

def get_box(sc= '', boxes = ''):
    """ For each file in the station_configuration, it associates the correct gridded box
          by comparing latitude and longitude """
    
    found_boxes = {}
    
    for c,r in sc.iterrows():
        print('Done row::: ' , c)
        lat = r.latitude
        lon = r.longitude
        file = r.file
        primary = r.primary_id
        
        for k in boxes.keys():
            #print(' Checking box ::::' , k )
            if k not in found_boxes.keys():     
                found_boxes[k] = {}
                found_boxes[k]['lat'] = []
                found_boxes[k]['lon'] = []
                found_boxes[k]['files'] = []
            try:
                if  lat < boxes[k]['lat'][0] or lat > boxes[k]['lat'][1]:
                    continue
                elif  lon < boxes[k]['lon'][0] or lon > boxes[k]['lon'][1]:
                    continue
                else:
                    found_boxes[k]['lat'].append(boxes[k]['lat'][0] ) , found_boxes[k]['lat'].append(boxes[k]['lat'][1] )
                    found_boxes[k]['lon'].append(boxes[k]['lon'][0] ), found_boxes[k]['lon'].append(boxes[k]['lon'][1] )
                    found_boxes[k]['files'].append( file ) 
            except:
                continue
            
    np.save('gridded_stations' , found_boxes )
    


#""" Chekc that all stations are asisgned a box """
#num = 0
#for k in assigned_boxes.keys():
#    num += len(assigned_boxes[k]['files'])
    
if os.path.isfile('gridded_stations.npy'):
    print('*** Loading the gridded stations :::')
    assigned_boxes = np.load('gridded_stations.npy' , allow_pickle= True ).item()
else:
    get_box(sc= sc, boxes = boxes)
    





def calculate_anomaly(file, years_window=20, min_days=10, min_months=10, variables= [85] ):
    """ main utiliy to calculate the monthly anomalies.
          Read the monthly file (on std pressure levels) and calculate the anomaly.
          input :: 
                         file = path to the file
                         years = number of years to consider when calculating the averages
                         min_days = minimum allowed number of days when calculating monthly averages
                         min_months =  minimum allowed number of months when calculating anomalies
          """

    std_plevs    = [1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 85000, 92500, 100000]
    
    to_drop = ['adjustment_id', 'advanced_assimilation_feedback', 'advanced_homogenisation', 'advanced_qc', 'advanced_uncertainty', 
                          'bbox_max_latitude', 'bbox_max_longitude', 'bbox_min_latitude', 'bbox_min_longitude', 'code_table', 'conversion_flag', 'conversion_method', 
                          'crs', 'data_policy_licence', 'date_time_meaning', 'exposure_of_sensor', 'latitude', 'location_method', 'location_precision', 'longitude', 'numerical_precision', 
                          'observation_duration', 'observation_height_above_station_surface', 'observation_id', 'original_code_table', 'original_units', 'original_value', 
                          'processing_level', 'quality_flag', 'report_id', 'secondary_variable', 'sensor_automation_status', 'source_id', 'spatial_representativeness', 
                          'traceability', 'units', 'z_coordinate_method',  'sensor_id',  'value_significance' ]
    
    obs_tab = xr.open_dataset(file , engine = 'h5netcdf' , group = 'observations_table' , decode_times = True, drop_variables = to_drop ).to_dataframe()
    
    print('Loading the observations_table as a pandas dataframe')
    #dates = np.unique(obs_tab['date_time'].values)
    #dates = pd.arrays.DatetimeArray(dates)
    
    """ Exracting the first date of the station """
    obs_tab['year']     = pd.arrays.DatetimeArray (obs_tab['date_time'].values[:] ).year
    obs_tab['month'] =pd.arrays.DatetimeArray (obs_tab['date_time'].values[:] ).month
    obs_tab['hour']    = pd.arrays.DatetimeArray (obs_tab['date_time'].values[:] ).hour
    
    first_date =  pd.to_datetime (obs_tab['date_time'].values[0] )
    first_date_year, first_date_month = first_date.year , first_date.month

    start_year = first_date_year + years_window # since I calclate anomalies I need to start from 20 years after the first available date 

    """
    # observation_value_reanalysis -> monthly averages calculated based on reanalyses deviations
    # secondary_value_reanalysis -> number of days used to calculate the monthly average. Must be > min_days value 
    # original_precision_reanalysis -> std of the the monthly average
    """
    anomaly, date_time, year, month, hour, average, num_months = [],[],[],[],[],[],[] 
    
    """ - loop over the variable,
          - loop over the time [0,12] 
          - loop over the available dates 
          - loop over the pressure levels """
    
    for v in variables:
        for t in [0,12]:  # loop over the time values i.e. 00 and 12
            
            df_red = obs_tab.loc [ (obs_tab['observed_variable']==v) 
                                   &  (obs_tab['hour']==t)  
                                   &  (obs_tab['secondary_value_reanalysis']> min_days) ] 
            
            dates = np.unique(df_red['date_time'].values)
            dates = pd.arrays.DatetimeArray(dates)
            
            for d in tqdm(dates):
                
                y, m = d.year , d.month 
                
                if y < start_year:
                    continue

                else:           
                    for m in range(1,13):  # loop over months
                        
                        """ First entry available for Lindenberg: 206/562, @Timestamp('1974-01-15 00:00:00'), 10 mindays, 10 min years  """
                        for p in std_plevs:
                            df = df_red.loc [ (df_red['z_coordinate']==p) 
                                                      & (df_red['month'] == m)
                                                      & (df_red['year'] > y - years_window )
                                                      & (df_red['year']  <= y - 1) ]
                        
                        """ This reduced df contains only the monthly averages for the particular variable, pressure, time, 
                              that are calculated using more than the minimum required number of days """
                        
                        if len(df) < 2:
                            continue
                        
                        if len(df) > min_months:
                            current_average =  df_red.loc [ (df_red['z_coordinate']==p) 
                                                      & (df_red['month'] == m)
                                                      & (df_red['year'] == y )  ]['observation_value_reanalysis'].values
                            
                            reference_average = np.mean(df['observation_value_reanalysis'].values)
                            
                            anomaly = current_average - reference_average
                            
                            
                            
                            print('Found!' )
                            
                        
                        
    
f = '/raid60/scratch/federico/MONTHLY_MEAN_NOVEMBER2020/0-20000-0-10393_monthly_averages_VARIABLE.nc'

calculate_anomaly(file=f, years_window=20, min_days=10, min_months=10, variables= [85] )

