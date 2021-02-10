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


def reduce_dataframe(obs_tab, variable=85, low_year = 1980 , max_year=1990, min_days=10, hour = 12):
    """ Reduces the total initial dataframe removing variables and values that are not needed """
    
    if not low_year:
        low_year = 1900
    if not max_year:
        max_year = 2020
        
    obs_tab['year']     = pd.arrays.DatetimeArray (obs_tab['date_time'].values[:] ).year
    obs_tab['month'] = pd.arrays.DatetimeArray (obs_tab['date_time'].values[:] ).month
    obs_tab['hour']    = pd.arrays.DatetimeArray (obs_tab['date_time'].values[:] ).hour
    
    df_red = obs_tab.loc [ (obs_tab['observed_variable']==v) 
                           &  (obs_tab['secondary_value_reanalysis']> min_days) 
                           &  (obs_tab['year'] >= low_year)  
                           &  (obs_tab['year'] <= max_year) 
                           &  (obs_tab['hour'] == hour ) ] 
    return df_red

    
def calculate_anomaly(obs_tab, years_window=20, min_months=10, month ='', year ='', hour ='' ):
    """ main utility to calculate the monthly anomalies.
          Read the monthly file (on std pressure levels) and calculate the anomaly.
          input :: 
                         file = path to the file
                         years = number of years to consider when calculating the averages
                         min_days = minimum allowed number of days when calculating monthly averages
                         min_months =  minimum allowed number of months when calculating anomalies
                         low_year, max_year = years of interest (of which I will calculate the anomaly)
          """

    #print('Loading the observations_table as a pandas dataframe')

    """ Exracting the first date of the station """
    #first_date =  pd.to_datetime (obs_tab['date_time'].values[0] )
    #first_date_year, first_date_month = first_date.year , first_date.month

    """
    # observation_value_reanalysis -> monthly averages calculated based on reanalyses deviations
    # secondary_value_reanalysis -> number of days used to calculate the monthly average. Must be > min_days value 
    # original_precision_reanalysis -> std of the the monthly average
    """

    """ First entry available for Lindenberg: 206/562, @Timestamp('1974-01-15 00:00:00'), 10 mindays, 10 min years  """
    
    averages, averages_bias, anomalies, anomalies_bias, plevels = [],[],[],[],[]
    
    for p in std_plevs:
            # df of the previous 20 years of which I calculate the climatology 
            climatology_df = obs_tab.loc [ (obs_tab['z_coordinate']==p) 
                               & (obs_tab['month'] == month)
                               & (obs_tab['hour'] == hour)
                               & (obs_tab['year'] > year - years_window )
                               & (obs_tab['year']  <= year - 1) ]
                        
            """ This reduced df contains only the monthly averages for the particular variable, pressure, time, 
                  that are calculated using more than the minimum required number of days """

            if len(climatology_df) > min_months:
                # dataframe of the year-month of which I want to calculate the anomaly  
                current_df = obs_tab.loc [ (obs_tab['z_coordinate']==p) 
                                   & (obs_tab['month'] == month)
                                   & (obs_tab['hour'] == hour)
                                   & (obs_tab['year'] == year) ]

                average = current_df['observation_value_reanalysis'].values
                average_bias =  current_df['observation_value_bias'].values
                
                climatology_average = np.mean(climatology_df['observation_value_reanalysis'].values)
                climatology_average_bias = np.mean(climatology_df['observation_value_bias'].values)
                
                anomaly = average - climatology_average
                anomaly_bias = average_bias - climatology_average_bias

            else:                
                average, average_bias, anomaly, anomaly_bias = np.nan, np.nan, np.nan, np.nan 
                
            averages.append(average) # no bias correction average
            averages_bias.append(average_bias) # bias corrected average
            anomalies.append(anomaly) # no bias correction anomaly 
            anomalies_bias.append(anomaly_bias) # bias corrected anomaly
            plevels.append(p)
                
    return averages, averages_bias, anomalies, anomalies_bias, plevels 
                                                       
    
f = '/raid60/scratch/federico/MONTHLY_MEAN_NOVEMBER2020/0-20000-0-10393_monthly_averages_VARIABLE.nc'

monthly_file_dir = '/raid60/scratch/federico/MONTHLY_MEAN_NOVEMBER2020/'

std_plevs    = [1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 85000, 92500, 100000]

def make_all_datetime(start_year=1905, end_year=2020):
    """ Create all date_time, on the day 15 of each month, at hours 00 and 12, 
          from start_year to end_year """
    date_time = []
    for y in range (start_year, end_year):
        for m in ['01' , '02' , '03' , '04' , '05' , '06' , '07' , '08' , '09', '10', '11', '12' ]:
            day = '15'
            for h in ['00','12']:
                #DT = np.datetime64('1978-01-02T06:00:00.000000000')
                timestamp = str(y) + '-' + m + '-' + day + 'T' + str(h) + ':00'
                ts = np.datetime64(timestamp)
                date_time.append(ts)
    return date_time
    
date_time = make_all_datetime(start_year=1905, end_year=2021)

def make_empty(results, date_time, std_plevs ):
    ''' empty dataframe, will work also as a fixed place holder '''
    results['res']['empty']['averages']         = ([np.nan] * len(date_time) * 16 )  
    results['res']['empty']['averages_bias'] = ([np.nan] * len(date_time) * 16)
    results['res']['empty']['anomalies']        = ([np.nan] * len(date_time) * 16)
    results['res']['empty']['anomalies_bias']= ([np.nan] * len(date_time) * 16)
    results['res']['empty']['plevels']             = ( std_plevs * len(date_time) )   

    results['res']['empty']['date_time'] = []
    for dt in date_time:
        results['res']['empty']['date_time'].extend( [dt] * 16 )  
        
    return results
    
for box in assigned_boxes.keys():   # loop over each grid box
    results = {}
    results['res'] = {}
    results['res']['empty'] = {}
    
    results = make_empty(results, date_time, std_plevs ) # make a default empty result 
    
    bbox = assigned_boxes[box]
    if not bbox['files']:
        continue

    for f in bbox['files']:              # loop over each station file (if any)
        station = f.split('sensor//')[1].split('_CEUAS')[0]
        file = monthly_file_dir + '/' + station + '_monthly_averages_VARIABLE.nc'
        if not os.path.isfile(file):
            #results = make_empty(results, date_time)
            print('!!! NOT FOUND ::: the file ',  file , '  could not be found')
            # remember that not for all the files it was possible to extract the standard_p_levels or the averages,
            # while the list of station in the boxes was extracted using only coordinates lat and lon.
            # So it makes sense that a large number of files is missing 
        else:
            for v in [85]: # looping over variables 
                for h in [0,12]:    
                    try:   
                        obs_tab = xr.open_dataset(file , engine = 'h5netcdf' , group = 'observations_table' , decode_times = True, drop_variables = to_drop ).to_dataframe()
                    except:
                        print('---> Wrong dataframe, skipping!' , file )
                        continue 
                    red_df = reduce_dataframe(obs_tab, variable= v, low_year = 1905 , max_year=2021, min_days=10, hour = h)
                    
                    if red_df.empty : 
                        continue
                        #results = make_empty(results, date_time)                
                    results['res'][station]  = {}
                    
                    results['res'][station]['averages']  , results['res'][station] ['averages_bias']  = [], []
                    results['res'][station]['anomalies'] , results['res'][station]['anomalies_bias'] = [], []
                    
                    results['res'][station]['plevels']      = []
                    results['res'][station]['date_time'] = []
                    
                    date_time_df = list(red_df['date_time'])
                    for dt in date_time: # all possible date_time starting from 1905 until 2020
                        if dt in date_time_df: # available date_time for this station 
                            
                            dt = pd.to_datetime(dt)
                            m = dt.month
                            y = dt.year
                            
                            average, average_bias, anomaly, anomaly_bias, plevels = calculate_anomaly(red_df, years_window=20, min_months=10, month =dt.month , year = dt.year , hour = h )
                       
                            for p, a, ab, an, anb in zip (plevels, average, average_bias, anomaly, anomaly_bias ):
                                results['res'][station]['averages'].append(a)  
                                results['res'][station]['averages_bias'].append(ab)
                                results['res'][station]['anomalies'].append(an)
                                results['res'][station]['anomalies_bias'].append(anb)
                                results['res'][station]['plevels'].append(p)
                                                                
                        else:
                            results['res'][station]['averages']          .extend ( [np.nan]  * 16 )  
                            results['res'][station]['averages_bias']  .extend ( [np.nan] * 16 )
                            results['res'][station]['anomalies']        .extend ( [np.nan] * 16 )
                            results['res'][station]['anomalies_bias'].extend ( [np.nan] * 16 )
                            results['res'][station]['plevels']             .extend ( std_plevs )   
                            
                            results['res'][station]['date_time'].extend ( [dt] * 16 )  
                            
    def make_write_box_df(results, box ):
        """ Combine data from differnet stations and save netcdf file output """
        output_df_columns =   ['average' , 'anomaly' , 'average_bias' , 'anomaly_bias' ]
        names_res              =   ['averages' , 'anomalies' , 'averages_bias' , 'anomalies_bias' ]
        output_df = pd.DataFrame(columns = ['z_coordinate' , 'date_time' , 'average' , 'anomaly' , 'average_bias' , 'anomaly_bias' ] )        
        stations = [ s for s in results['res'].keys() if 'empty' not in s ]

        """ Here, loop over the calculated  results.
              If no stations are available, use dummy nan results. """
        for column, name in zip ( output_df_columns , names_res ):
            # need to add a sum!!! 
            means = []
            if stations:             
                for v in range( len(results['res']['empty']['date_time'] )  ) :
                    values = [ results['res'][s][name][v] for s in stations ]

                    average = np.nanmean(values)
                    if not average:
                        average = np.nan 
                    try:
                        means.append(average[0])
                    except:
                        means.append(average)
                        
                #val = np.array( [ len(results['res']['empty']['date_time']) * [np.nan] ] ).reshape(len(results['res']['empty']['date_time']) , 1)
                #values =   [ val + np.array(results['res'][s][name]).reshape( len(results['res']['empty']['date_time']) , 1 )  for s in stations ][0]
                
            else:
                means = results['res']['empty'][name]               
            
            output_df[column] = means 
        
        output_df['z_coordinate'] = results['res']['empty']['plevels']
        output_df['date_time']     = results['res']['empty']['date_time']
        lat, lon = boxes[box]['lat'][0] , boxes[box]['lon'][0]
        size      = abs( boxes[box]['lat'][0] - boxes[box]['lat'][1] )
        
        Lat , Lon = lat + size/2 , lon + size/2
        output_df['latitude']    = ( [Lat] * len(results['res']['empty']['date_time'] ) )
        output_df['longitude'] = ( [Lon] * len(results['res']['empty']['date_time'] ) )
        output_df['grid_size']  = ( [size] * len(results['res']['empty']['date_time'] ) )

        """ Saving array with stations """
        if stations:
            stat = '_'.join(stations)
            0
        else:
            stat = 'None'
        output_df['primary_ids']      = ( [stat] * len(results['res']['empty']['date_time'] ) )

        box_name = 'lat_lon_size_' + str(Lat) + '_' + str(Lon) + '_' + str(size) + '_' + stat + '_stations.nc'
        xarray = output_df.to_xarray()
        out = xarray.to_netcdf ( 'output_monthly/' + box_name + '.nc' , mode = 'w' )
        # xarray =output_df[43500:43600]
        # check::: 43536  1000         2018-05-15 12:00:00  []  
     
        
    done = make_write_box_df(results, box )     
    print('+++ File done ' , file , 'for box ' , box )
        
    '''
    try:
        done = make_write_box_df(results, box )     
        print('+++ File done ' , file , 'for box ' , box )
        
    except:
        print('+++ File FAILED ' , file , 'for box ' , box )        
        pass 
        # Loop over all possible dates from 1905 till 2020 ???
    '''    
        
        

                
 
    