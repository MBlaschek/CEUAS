import os
import sys
import netCDF4 as nc
import pandas as pd
pd.set_option('display.max_rows', 50)
pd.set_option('display.max_columns', 20)
pd.set_option('display.width', 300)
from pathlib import Path

import numpy as np
#import argparse
from datetime import datetime, timedelta  
import numpy.ma as ma
#import math
import h5py as h5py
import xarray as xr 
#from numba import njit
import psutil
import copy
from numba import njit
import code
import urllib.request
from functools import reduce
from tqdm import tqdm

from multiprocessing import Pool
from functools  import partial

from numba import njit
from numba.typed import Dict
from numba.core import types
# http://numba.pydata.org/numba-doc/latest/reference/pysupported.html#typed-list    --> to use dictionaries 

from collections import Counter

import cProfile
#cProfile.run('__main__')

import time as T 
t=T.time()


sys.path.append('../harvest/code')
#from harvest_convert_to_netCDF_newfixes import load_cdm_tables , write_dict_h5 

#from Desrozier import * 


# nan int = -2147483648 


import warnings
warnings.simplefilter(action='ignore', category=FutureWarning) # deactivates Pandas warnings 

import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')  # up to INFO level, DEBUG statements will not be printed 

import argparse


"""
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)
"""

process = psutil.Process(os.getpid())

# https://stackoverflow.com/questions/287871/how-to-print-colored-text-in-python
cend   = '\033[0m'
blue   = '\033[34m'
red     = '\33[91m'
green = '\33[92m'


def datetime_toseconds(date_time):
    """ Converts a generic date_time array to seconds since '1900-01-01 00:00:00' """
    offset = np.datetime64('1900-01-01 00:00:00')       
    
    ### I cant put nans as int32 with the largest number toherwise it will still convert it to a date_time which makes no sense 
    to_seconds = [] 
    for dt in date_time:
        if dt != 0:
            try:
                delta = np.datetime64(dt) - offset 
                to_seconds.append(  delta.item().total_seconds()  )                
            except:
                to_seconds.append(0)
        else:
            to_seconds.append(0)
    a = np.array(to_seconds).astype(np.int64)
    return a # replacing with seconds from 1900-01-01 00:00:00     

def remove_outliers(data= '', min_p= 25, max_p= 75, cut= 1, skewed= False):
    """ Finds outliers, and replace them with np.nan (to keep vector of same length)                                                                                                                                                                                              

         input ::       data = list of values 
                           min_p , max_p = minimum and maximum values of the percentile to consider to determine outliers 
                           skewed = use True to consider a skewed (not symmetricla Gaussian) distribution 
                           cut = factor to allow slight deviation from given percentiles 
         returns ::   cleaned   = list of values without outliers                                                                                                                                                                    
                           outliers   = list of outlier values                                                                                                                                                                                                               
                           lower,upper, median = outliers delimiter and median values """

    q_min, q_max = np.nanpercentile(data, min_p), np.nanpercentile(data, max_p)
    cut_off = (q_max - q_min) * cut
    lower, upper = q_min-cut_off, q_max+cut_off

    if skewed==True:
        q50 = np.nanpercentile(data, 50)
        lower , upper = q_min-(q50-q_min)*cut ,  q_max+(q_max-q50)*cut  # the higher the cut, the more relaxed the contition for exclusion 

    median = np.nanmedian(data)
    cleaned, outliers = [],[]

    for d in np.asarray(data):
        if d >= lower and d <= upper:
            cleaned.append(d)

        else: # only storing non nans values 
            if not np.isnan(d):
                outliers.append(d)
                
    return cleaned, outliers, lower, upper, median




class Monthly_Average(object):
    
    def __init__(self, out_dir = '', file='', variable = '' ):
        """ Initializing the file
              var = name of the variable inside the file, e.g. temperature = 'ta' """
        station_id =  [ f.split("_")[1] for f in file.split('/') if 'dest' in f ][0]
        
        self.out_dir = out_dir 
        
        if not os.path.isdir(out_dir):
            print('+++ Creating output directory: ' , out_dir )
            os.mkdir(out_dir)
            
        self.station_id = station_id
        self.file = file 
        #self.summary_file =  station_id + '_' + datetime.now().strftime("%Y-%m-%d_%H-%M-%S") + '_summary.txt'
        self.std_plevs    = [1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 85000, 92500, 100000]
        self.variable = variable
        self.retrieved_variables = ['bias_estimate', 'lat', 'lon', 'obs_minus_an', 'obs_minus_bg', 'plev', 'sonde_type', 'time', self.variable]
        self.monthly_data_columns =  [self.variable , 'time' , 'plev' , 'sonde_type']
        self.monthly_data_columns_full = self.monthly_data_columns + [ 'lat', 'lon']
        
        # list of columns that will be used in the output files 
        
        
    def load(self):
        """ Read the input data and store in a dictionary of numpy arrays """
        data = {}

        file_data = xr.open_dataset(self.file )

        for v in self.retrieved_variables:
            if v == 'sonde_type':
                temp =  file_data[v].values.view('S3')
                values = [d[0] for d in temp]
            else:
                values = file_data[v].values
                
            data[v] = list(values)
            
        self.data = data 
        file_data.close()
        

    def get_dates_indices(self):
        """ This function analyzes the time of the observations, and extracts the corresponding indices."""
        time = self.data['time']
        timestamps, indices = np.unique(time, return_index=True) # unique timestamps and their starting indices 

        indices_plus = [i-1 for i in indices[1:]]# getting the upper index of each record  
        indices_plus.append (len(self.data['time']) )  # adding the last index 
        
        dates_indices = {}
        limit = 0
                
        for rts, index, index_p  in zip ( tqdm(timestamps[limit:]) , indices[limit:], indices_plus[limit:] ) :

            pd_rts = pd.DatetimeIndex([rts]) # using pandas to extract day, month, year and hour 
            
            this_hour = pd_rts[0].hour # if the observation is after 10 pm we report it to the following day. Will handle automatically case of February 
            
            if not ( this_hour >= 22 or ( this_hour >= 10 and this_hour <= 14 ) or this_hour == 0 or ( this_hour <= 2 ) ):
                #print(' Time is too far from limits, cannot use this hour ', this_hour)
                continue      # the distance of the hour from 00 or 12 is larger than 2 hrs to be included in the averages. Must skip the record 
            #else:
            #    pass
            #    #print('Will keep this hour! ')
                
            if this_hour >= 22 or this_hour == 0 or this_hour <= 2 :
                pd_rts = pd_rts + timedelta(minutes=120)
                std_time = '00'
            else:
                std_time = '12'              
                         
            date_string = str(pd_rts[0].year) + '_' + str(pd_rts[0].month) 
            
            if date_string not in dates_indices.keys():
                dates_indices[date_string] = {}
                
                dates_indices[date_string]['00'] = {}
                dates_indices[date_string]['12'] = {}
                
                dates_indices[date_string]['00']['indices'] = []
                dates_indices[date_string]['12']['indices'] = []
                
                dates_indices[date_string][std_time]['indices'].append( [index,index_p] )
            
            else:
                dates_indices[date_string][std_time]['indices'].append( [index,index_p] )
            
        return dates_indices


    def write_monthly_file(self):
        out_file = self.out_dir + '/' + self.station_id + '_monthly_averages_' + self.variable + '.nc'
        
        os.system('rm ' +  out_file )
        index = np.array ( range(0,len(self.monthly_table['time'])) , dtype='S1')       
        
        # converting to appropriate variable types
        var_type = {}
        for var in self.monthly_table.keys(): 
            var_type[var] = ''
            if var == 'time':
                kind = np.int64
            else: 
                if var != 'sonde_type':
                    kind = np.float32      
                else:
                    kind = np.dtype('|S3')
            var_type[var] = kind
        
        '''
            data['observations_table'].create_dataset('sensor_id', data = s_id_vec.view('S1').reshape( len(s_id_vec), 3 )  , compression = 'gzip' ,  chunks=True)                  
                data['observations_table']['sensor_id'].dims[0].attach_scale( data['observations_table']['index'] )  
                data['observations_table'].create_dataset( s ,  data=stringa[:slen]  )                
                data['observations_table']['sensor_id'].dims[1].attach_scale(data['observations_table']['string{}'.format(slen)])
                data['observations_table']['string{}'.format(slen)].attrs['NAME']=np.bytes_('This is a netCDF dimension but not a netCDF variable.')             
                print(' *** Done with the attributes of the dimension *** ')
        '''
        
        with h5py.File(out_file, 'a' ) as out_netcdf:
            out_netcdf.create_dataset('index', index.shape, index.dtype, compression= 'gzip', chunks=True)
            for var in self.monthly_table.keys(): 
                print('--- Written var ' , var )
                data = np.array(self.monthly_table[var]).astype(var_type[var])
                out_netcdf.create_dataset(var, data = data, compression= 'gzip', chunks=True)
        
                try:
                    out_netcdf[var].dims[0].attach_scale(out_netcdf['index'])
                except:
                    print('Failed dimensions +++++++++++++++++= ')        
                    
            out_netcdf['time'].attrs['units'] = np.bytes_('seconds since 1900-01-01 00:00:00')            
                    
        #if variable in self.MergedFile.attr_dic.keys():
        #    out_netcdf['observations_table'][variable].attrs['description']       = np.bytes_(self.MergedFile.attr_dic[variable]['description'])
        #    out_netcdf['observations_table'][variable].attrs['external_table' ] = np.bytes_(self.MergedFile.attr_dic[variable]['external_table'])
                        
        out_netcdf.close() 
        print('Written output file in ', out_file )
        
    def make_monthly_data(self):
        """ Main utility to produce monthly averages observations_table. 
              We only consider one variable at a time and write separated netCDF files """
        
        ''' # not in use atm
        def get_mostcommon_value(indices, column, index_min='', index_max=''):
            """ Extracts the most common values (in a given moth) for the variable considered. 
                  E.g. the sensor_id of the monthly average would be the most common sensor id found in the records for that specific  month. 
                  Indices is the list of indices, of the chucnked obs_tab[index_min:index_max] where I have the observation values.
                  So index_min+indices[0] is the first available value of the desired variable. """
            
            # HAD TO SIMPLY to make it faster: returns the first value of the list 
            # I simply pick the first value available for that chunk, for these three columns. Otherwise I set it to nan 
            if column in ['sensor_id', 'latitude', 'longitude']:
                most_common =self.data['obs_tab'][column][index_min ].values
            else:
                most_common = np.nan 
            return most_common
        '''
        
        def get_observations(indices = ''):
            """ Loop over the indices of the records and extract the specific indices of the observations, and the values of the observations,
                  and the departures,
                  indices = list of two lists, i.e. [ind_min, ind_max ]"""

            #  ['bias_estimate', 'lat', 'lon', 'obs', 'obs_minus_an', 'obs_minus_bg', 'plev', 'sonde_type', 'time', self.variable]
            results = {}  
            results['is_there_data'] = False 
            columns = ['bias_estimate', 'obs_minus_an', self.variable, 'observation_indices' ] # columns that will be stored in the results. Time will be added later on. 
            
            for i in range(len(indices[0]) ): # Looping over the entire list of indices 
                
                ind_min, ind_max = indices[0][i] , indices[1][i]+1  # you have to increase by 1 unit otherwise you miss the last element when slicing 
                plevels = np.array(self.data['plev'])[ind_min: ind_max]  # available pressure levels in the record between ind_min:ind_max 
                #time = self.data['time'][ind_min: ind_max] # is it needed? Time of the observation for that particular record between ind_min:ind_max 
                
                for p in self.std_plevs:  # Loop over the stanard pressure level list             
                    if p not in results.keys(): # check if the pressure level has already been checked  
                        results[p] = {}      
                        for v in columns:
                            results[p][v] = []   
                        
                    if p not in list(plevels): # there is NO data for the selected pressure level. Filling with np.nan values 
                        for v in columns :
                                results[p][v].append(np.nan)
                           
                    else: # case where the standard p level IS  availabe in the data pressure levels  
                        unique_ind = np.where(plevels == p)[0] # selecting the specific observation for the selected pressure level 
                        # ['bias_estimate', 'obs', 'obs_minus_an', self.variable ] 
                        if len(unique_ind)>=1: # this if should be redundant, I must always have >=1
                            val = self.data[self.variable][ind_min: ind_max][unique_ind[0]]
                            
                            results[p][self.variable].append(val)
                            results[p]['observation_indices'].append(unique_ind[0])
                            
                            results['is_there_data'] = True 
                            
                            dep         = self.data['obs_minus_an'][ind_min: ind_max][unique_ind[0]] # departures from era5fb
                            era5_b    = self.data['bias_estimate'][ind_min: ind_max][unique_ind[0]]
                                
                            results[p]['obs_minus_an']             .append(dep)        
                            results[p]['bias_estimate'].append(era5_b)

                        else:
                            for v in columns:   
                                results[p][v].append(np.nan)      
                 
            return results , ind_min, ind_max # returns the dictionary of the results and the correpsonding new indices in the results 
                        
        def fill_columns(monthly_table , index = '', dt = '' , data = True ):
            """ Fill the dictionary of the observations_table. 
            When there is no data, all the columns except time and pressure levels are set to nans.
            Otherwise, they will be filled with the proper values taken with the indices from the data lists. """

            # these are to be always filled
            tot_observations = len(self.std_plevs) # 16 std pressure levels 
            monthly_table['plev'] .extend ( self.std_plevs )
            monthly_table['time'] .extend ( [dt] * tot_observations ) 
            
            # nan values for data missing for specific dates
            if not data:
                all_columns = self.averages_columns
                for c in all_columns:
                    monthly_table[c].extend  ( [np.nan] * tot_observations )
                monthly_table['sonde_type'].extend  ( [b'NA '] * tot_observations )                    
                monthly_table['lat'].extend  ( [self.data['lat'][0]] * tot_observations )
                monthly_table['lon'].extend ( [self.data['lon'][0]] * tot_observations )
                
            else:
                monthly_table['lat'].extend  ( [self.data['lat'][index]] * tot_observations )
                monthly_table['lon'].extend ( [self.data['lon'][index]] * tot_observations )
                monthly_table['sonde_type'].extend ( [self.data['sonde_type'][index]] * tot_observations )
                
            return 
        
        
        def get_all_obs_per_month( ):
            """ Extract all the values for all the months, pressure levels and time.
                  Will be used to calculate the global statistics. ERA5 biases are not needed in this step. 
                  Return: a dictionary, for each month and time, containing the list of all values """
            
            # ['bias_estimate', 'lat', 'lon', 'obs', 'obs_minus_an', 'obs_minus_bg', 'plev', 'sonde_type', 'time' ]
            monthly_dates = self.monthly_dates
            
            #z_coord          = self.data['plev']
            #obs_val          = self.data[self.variable]
            #std_plevels    = self.std_plevs  

            results = {}                
            for i in range(1,13): # lop over the months of the year 1,...,12 
                results[i] = {}                
                    
                for m in monthly_dates.keys(): # loop over all the months_years in the data e.g. 1970_11 , 1970_12 etc. 
                    month = int(m.split('_')[1])
                    #year    = int(m.split('_')[0]) # not needed ?
                    if month == i:
                        for time in ['00','12'] :
                            if time not in results[i].keys():
                                results[i][time] = {}  
                            
                            indices = monthly_dates[m][time]['indices'] 

                            if indices:  
                                ind = np.array ( ([ i[0] for i in indices ] , [ i[1] for i in indices ]  ) , dtype = np.int32 )
                                extracted_data , index_min, index_max = get_observations(indices = ind )
                                # extracted_data[1000].keys() -> dict_keys(['bias_estimate', 'obs', 'obs_minus_an', 'ta', 'observation_indices'])
                                for p in self.std_plevs:
                                    if p not in results[i][time].keys():
                                        results[i][time][p] = {}
                                        results[i][time][p]['all_observations'] = []
                                        results[i][time][p]['all_departures']    = []
                                    for v,d in zip(extracted_data[p][self.variable], extracted_data[p]['obs_minus_an'] ) :
                                            results[i][time][p]['all_observations'].append(v)
                                            results[i][time][p]['all_departures'].append(d)             
                                            
            return results                
        
        def global_statistics(all_obs):
            """ Calculate the average, std_dev, quartiles for the complete set of observations and departures (all available data) """
            
            stat = {}
            for i in range(1,13) :
                stat[i] = {}                
                for time in ['00','12'] :
                    stat[i][time] = {}                
                    for p in self.std_plevs:
                        stat[i][time][p] = {}                
                        
                        try:
                            data = all_obs[i][time][p]['all_observations']  # if the month_year is not available at all, the data dict does not have this key 
                        except:
                            data = []
                            
                        """ calculating statistics of observed values """
                        if data:                              
                            cleaned_obs, outliers_obs , lower_obs, upper_obs, median_obs = remove_outliers(data= np.array(data), min_p=25, max_p=75, cut= 2, skewed=False)
                            mean, std = np.nanmean(cleaned_obs)  , np.std(cleaned_obs) 
                        else:
                            mean, std , lower_obs, upper_obs = np.nan, np.nan, np.nan, np.nan 
                            
                        stat[i][time][p]['obs_average'] = mean       
                        stat[i][time][p]['obs_std_dev']  = std    
                        stat[i][time][p]['obs_lower']      = lower_obs    # lower threshold for valid data (i.e. data < lower_obs are considered outliers ) 
                        stat[i][time][p]['obs_upper']     = upper_obs    # upper threshold for valid data (i.e. data > upper_obs are considered outliers )                       
                        
                        """ calculating the statistics of departures  """
                        try:
                            data_dep = all_obs[i][time][p]['all_departures']      
                        except:
                            data_dep = []
                            
                        if data and data_dep: # calculating statistics of departures  
                            cleaned_dep, outliers_dep, lower_dep, upper_dep, median_dep = remove_outliers(data= np.array(data_dep), min_p=25, max_p=75, cut= 2, skewed=False)
                            mean_dep, std_dep = np.nanmean(cleaned_dep)  , np.std(cleaned_dep)
                        else:
                            mean_dep, std_dep, lower_dep, upper_dep =  np.nan, np.nan, np.nan, np.nan
                            
                        stat[i][time][p]['dep_average'] = mean_dep        
                        stat[i][time][p]['dep_std_dev']  = std_dep  
                        stat[i][time][p]['dep_lower']     = lower_dep     # lower threshold for valid data (i.e. data < lower_obs are considered outliers ) 
                        stat[i][time][p]['dep_upper']     = upper_dep    # upper threshold for valid data (i.e. data > upper_obs are considered outliers )                  
                          
            print(red + '*** Calculated global statistics ***' + cend)
            return stat            
                   

        monthly_dates = self.get_dates_indices() # Create the monthly dates and time (every 15th of the month, 12 and 00 time )
        # monthly_dates = {'1966_11':{'00': {'indices': [[0, 5], [6, 10], [23, 28], [41, 46], [53, 58], [74, 78]]}, 
        #                                                  '12': {'indices': [[11, 16], [17, 22], [29, 34], [35, 40], [47, 52], [59, 62], [63, 67], [68, 73]]}}, ...
        # So inside the original data list, I find the right data located at these indices 
        self.monthly_dates = monthly_dates 
        
        #z_coord          = self.data['plev']
        #obs_val          = self.data[self.variable]
        #std_plevels    = self.std_plevs                 
        #era5fb_departures = self.data['obs_minus_an']
        #obs_era5       = self.data['obs']  # the 'obs' variable is not used right now 
        #era5_bias       = self.data['bias_estimate']
        
        print(blue + "*** Done Loading Data ::: "  + str(T.time()) + cend )
        
        """ Calculating the global statistics i.e. for all the avaiable months """
        print ('*** Extracting the monthly data\n')
        all_obs = get_all_obs_per_month()  # extracts all the observations and departures for each plevel
        self.global_statistics = global_statistics(all_obs)  # calculate all the global statistics
        
        
        monthly_table = {} # dictionary to store final averages results
        for c in self.monthly_data_columns_full:
            monthly_table[c] = []
            
        # names of columns containing the averages in the output
        averages_columns = [self.variable               , self.variable + '_std_dev'          , self.variable + '_num_obs' , 
                                           self.variable + '_glob' , self.variable + '_std_dev_glob' , self.variable + '_num_obs_glob' , 
                                           self.variable + '_dep'  , self.variable + '_std_dev_dep'  , self.variable + '_num_obs_dep'  , 
                                           self.variable + '_bias' , self.variable + '_std_dev_bias'  , self.variable + '_num_obs_bias' ]
        
        self.averages_columns = averages_columns
        
        for c in averages_columns : 
            monthly_table[c] = []

          
        # TODO chage the loop so it is fixed between 19000101 and 20210101 
        
        for m in tqdm(monthly_dates.keys() ):  # monthly_dates.keys() = ['1966_11' , '1966_12' , ... ]
            value = monthly_dates[m] # loop over the months    
            month = int(m.split('_')[1] )
            print("*** Processing: " , m , '   ', T.time() )
            for time in ['00','12'] : # loop over the times 

                dt = datetime.strptime(m + '_15 ' + time + ':00:00', '%Y_%m_%d %H:%M:%S')   # create the date time for the 15th day of that particular month at 12 or 00               
                indices = value[time]['indices']         

                if not indices: # if I have no indices, I fill every column with nans except the time,plevels,lat,lon
                    dummy = fill_columns(monthly_table  , dt = dt, data = False ) 
                    continue
                
                elif indices:  # check if I have records within this period of time  
                                            
                    ind = np.array ( ([ i[0] for i in indices ] , [ i[1] for i in indices ]  ) , dtype = np.int32 )  # can simplify this, it is a leftover from trying numba  

                    extracted_data , index_min, index_max = get_observations(indices = ind )
                    
                    if extracted_data['is_there_data']: # check if I have at least one valid observation, otherwise do nothing 
                        
                        dummy = fill_columns(monthly_table, index = index_min , dt = dt , data = True ) 

                        """ loop over the available std pressure levels """
                        for p in self.std_plevs: 
                            #print("month, level, hour " , m , '  ' , p , '  ' , time , '  ' , T.time())

                            values = [f for f in extracted_data[p][self.variable] if not np.isnan(f)  ]
                            
                            if values:
                                """ Method 1: calculate average and outliers for each month individually - might suffer from low statistics """
                                values_cleaned, outliers, lower, upper, median = remove_outliers(data= np.array(values), min_p=25, max_p=75, cut= 2, skewed=False)
                                mean, std, Len = np.nanmean(values_cleaned), np.std(values_cleaned), len(values_cleaned)

                                """ Method 2: calculate average and outliers for all the available months together """
                                lower_global , upper_global = self.global_statistics[month][time][p]['obs_lower'] , self.global_statistics[month][time][p]['obs_upper'] 
                                values_global_stat = [ v for v in values if v > lower_global and v < upper_global ] # remove form current data the outliers 
                                
                                """
                                if values_global_stat:
                                    mean_glob, std_glob, len_glob = np.nanmean(values_global_stat) , np.std(values_global_stat) , len(values_global_stat) 
                                else:
                                    mean_glob, std_glob, len_glob = np.nan, np.nan, np.nan 
                                """
                                mean_glob, std_glob, len_glob = np.nanmean(values_global_stat) , np.std(values_global_stat) , len(values_global_stat) 
                                
                                """ Method 3: if available, check departure statistics, and remove observation data if the corresponding departure is an outlier """
                                observations = extracted_data[p][self.variable]
                                departures, biases =  extracted_data[p]['obs_minus_an'], extracted_data[p]['bias_estimate']
                                
                                values_departures = [] 
                                values_departures_unbiased = [] 
                                for o,d,b in zip(observations, departures, biases ):
                                    if not np.isnan(o) and not np.isnan(d):
                                        # check if the departures values are not outliers, so that I keep the observations
                                        if d > self.global_statistics[month][time][p]['dep_lower'] and d < self.global_statistics[month][time][p]['dep_upper']:
                                            values_departures.append(o)
                                        if not np.isnan(b):
                                            unbiased = o + b 
                                            values_departures_unbiased.append(unbiased) # TODOhere
                                                
                                if values_departures:
                                    mean_dep, std_dep, len_dep = np.nanmean(values_departures) , np.std(values_departures) , len(values_departures)
                                else:
                                    mean_dep, std_dep, len_dep = np.nan , np.nan , np.nan 
                                    
                                if values_departures_unbiased:
                                    mean_bias, std_bias, len_bias = np.nanmean(values_departures_unbiased) , np.std(values_departures_unbiased) , len(values_departures_unbiased)
                                else:
                                    mean_bias, std_bias, len_bias =  np.nan , np.nan , np.nan 
                                            
                            else:
                                mean, std, Len                         = np.nan , np.nan , np.nan 
                                mean_glob, std_glob, len_glob = np.nan , np.nan , np.nan 
                                mean_dep, std_dep, len_dep    = np.nan , np.nan , np.nan 
                                mean_bias, std_bias, len_bias  = np.nan , np.nan , np.nan
                                
                                #for column in other_columns: # other_columns: columns for which we want to keep the most freque value, e.g. latitude or sensor_id 
                                #    monthly_observations_table[column].append(np.nan)
                                
                            """ Filling the values """    
                            # example columns ['ta' , 'ta_std_dev' , 'ta_num_obs'] , ['ta_glob' , 'ta_std_dev_glob' , 'ta_num_obs_glob'] , 
                            # ['ta_dep' , 'ta_std_dev_dep' , 'ta_num_obs_dep'] , ['ta_bias' , 'ta_std_dev_bias' , 'ta_num_obs_bias'] 
                            
                            monthly_table[self.variable].append( mean )
                            monthly_table[self.variable + '_std_dev'] .append( std ) 
                            monthly_table[self.variable + '_num_obs']  .append( Len )
                            
                            monthly_table[self.variable + '_glob'].append(mean_glob)
                            monthly_table[self.variable + '_std_dev_glob'] .append(std_glob ) 
                            monthly_table[self.variable + '_num_obs_glob']  .append( len_glob )            
                            
                            monthly_table[self.variable + '_dep'].append(mean_dep)
                            monthly_table[self.variable + '_std_dev_dep'] .append(std_dep ) 
                            monthly_table[self.variable + '_num_obs_dep']  .append( len_dep )  
                            
                            monthly_table[self.variable + '_bias'].append(mean_bias)
                            monthly_table[self.variable + '_std_dev_bias'] .append(std_bias ) 
                            monthly_table[self.variable + '_num_obs_bias']  .append( len_bias )  
                            
                    else:
                        dummy = fill_columns(monthly_table,  index = index_min , dt = dt , data = False ) 
                        

        monthly_table['time'] =  datetime_toseconds(monthly_table['time'] ) 
        #or c in ['source_id' ]:
        #    monthly_table[c].extend([np.nan] * len(monthly_table['date_time']) )
            
        self.monthly_table = monthly_table

        for c in monthly_table.keys():
            print(c , '     ' , len(monthly_table[c] ) )
        """ Finally writing the output file """ 
        print(0)
        dummy = self.write_monthly_file()
            

        
  


files_directory = '/raid60/scratch/federico/CDS_DATABASE_01FEB2021/temperature'
out_dir = '/raid60/scratch/federico/MONTHLY_FEB2021'
 
if __name__ == '__main__':
        
            parser = argparse.ArgumentParser(description="Utility to create Monthly Averages Files")

            parser.add_argument('--force_run' , '-f', 
                                  help="Force running the file(s)"  ,
                                  type = str,
                                  default = 'False' )
            
            parser.add_argument('--monthly_averages' , '-m', 
                                  help="Calculate monthly averages"  ,
                                  type = str,
                                  default = 'False'  )    
            
            args = parser.parse_args()
            
            force_run                = args.force_run
            monthly_averages  = args.monthly_averages

            stations_list = [f for f in os.listdir(files_directory) if '.nc' in f ]
            stations_list = [f for f in stations_list if '82930' in f ]
            
            POOL = False  
            
            if POOL:
                files = [ files_directory + '/' + s for s in stations_list ]
                def run(out_dir, variable, file ):
                    try:
                        Average = Monthly_Average( out_dir = out_dir, file= file , variable = variable) 
                        load = Average.load()   
                        monthly = Average.make_monthly_data()
                        
                    except:
                        print('Failed file: ' , file )
                        a = open('Failed_monthly_extraction.txt' , 'a+')
                        a.write(file + '\n')
                
                p = Pool(30)
                func = partial(run, out_dir, 'ta')
                out = p.map(func, files)                      
        
        
            else:
                for s in stations_list:
                    file =  files_directory + '/' + s 
                    
                    Average = Monthly_Average( out_dir = out_dir , file= file , variable = 'ta' ) 
                    load = Average.load()   
                    monthly = Average.make_monthly_data()
                            
          
            
                        
