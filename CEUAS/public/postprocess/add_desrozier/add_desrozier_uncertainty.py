""" This postprocessing script reads newly-sorted netCDF files,
      calculates Desrozier's statistics for the standard pressure level and 
      writes them in a new group called advanced_uncertainty """


import os
import sys
import netCDF4 as nc
import pandas as pd
from pathlib import Path

import numpy as np
from datetime import datetime, timedelta  
import numpy.ma as ma
import h5py as h5py
import xarray as xr 
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

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)

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

data = {}


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

 
def remove_outliers(data= '', min_p= 25, max_p= 75, cut= 1, skewed= False, only_clean = True):
    """ Finds outliers, and replace them with np.nan (to keep vector of same length)                                                                                                                                                                                              

         input ::       data = list of values 
                           min_p , max_p = minimum and maximum values of the percentile to consider to determine outliers 
                           skewed = use True to consider a skewed (not symmetricla Gaussian) distribution 
                           cut = factor to allow slight deviation from given percentiles 
         returns ::   cleaned   = list of values without outliers                                                                                                                                                                    
                           outliers   = list of outlier values                                                                                                                                                                                                               
                           lower,upper, median = outliers delimiter and median values.
                           if only_clean == True, return only list of cleaned values. """

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
                
    if only_clean:
        return cleaned
    else:
        
        return cleaned, outliers, lower, upper, median




class Desroziers():
    
    def __init__(self, out_dir = '' , station_id = '' , file='' ):
        self.out_dir = out_dir 
        
        if not os.path.isdir(out_dir):
            print('+++ Creating output directory: ' , out_dir )
            os.mkdir(out_dir)
            
        self.station_id = station_id
        self.file = file 
        #self.summary_file =  station_id + '_' + datetime.now().strftime("%Y-%m-%d_%H-%M-%S") + '_summary_desrozier.txt'
        self.std_plevs = [1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 85000, 92500, 100000]
        self.variables = [85,104,105,107] # 104,105 u and v wind, 106 wind from direction, 107 wind speed
        #self.dic_type_attributes = np.load('../merge/dic_type_attributes.npy', allow_pickle= True).item()
        #self.encodings = np.load('../merge/groups_encodings.npy' , allow_pickle = True ).item()
        
    def load_data(self):
        
        """ Global container of the data.
              From observations_table: ['date_time', 'z_coordinate', 'z_coordinate_type' , 'observed_variable', 'observation_value'] 
              From era5fb: ['an_depar@body' , 'fg_depar@body', 'biascorr@body']  
              
        """
        data = {}
        obs_var = ['date_time', 'z_coordinate', 'z_coordinate_type' , 'observed_variable', 'observation_value'] 
        era5fb_var = ['an_depar@body' , 'fg_depar@body', 'biascorr@body']
        
        
        h5py_file = h5py.File(self.file, 'r+')
        
        self.h5py_file = h5py_file 
        
        """ Reading only needed input data and store it inside an xarray """
        #data = xr.Dataset()
        
        for v in obs_var:
            data[v] = h5py_file['observations_table'][v]
        
        for c in era5fb_var:
            data[c] = h5py_file['era5fb'][c]
            
        self.recordindices = h5py_file['recordindices']
        
        #df = pd.DataFrame.from_dict(data)
        
        self.data = data 
      
        
    def calculate_Desroziers_errors(self):
        
        """ Errors xarray """
        #errors_out = xr.Dataset()
        errors_out = {}
        
        error_day_window= ['30','60','90','180']
        #error_day_window= ['30']
                
        """ Assuming a day window around each observation to calculate the time-average
             e.g. 30 days around 30 january = from 15 Jan to 14 Feb. """
        
        """ Creating the output place-holder filled with np.nan """
        for d in error_day_window: #  ['30','60','90','180']
            error, num = 'desroziers_' + d , 'num_' +d
            errors_out[error]     = np.full( (len(self.data['date_time'] ) ) , np.nan ) # create place holder for uncertainties filled with nans 
            errors_out[num]      = np.full( (len(self.data['date_time'] ) ) , -2147483647 ) 
            #errors_out[std_dev] = np.full( (len(self.data['date_time'] ) ) , np.nan ) # not in use

        for var in self.variables:
            ind = self.recordindices[str(var)][:]
            imin, imax = min(ind) , max(ind) + 1000 # I do not know the maximum, so I take max+1000
                       
            red_dic = {}
            for col in self.data.keys():
                red_dic[col] = self.data[col][imin:imax]
            
            """ Building a data frame usin gonly data between the correct indices """
            df = pd.DataFrame(red_dic)
            
            start_date =  np.datetime64('1900-01-01T00:00') 
            
            df = df.loc[df['observed_variable'] == var ] # since I do not know a priori when the indices of the next variable start, check the observed_var number
            
            delta = pd.to_timedelta(df['date_time'] , unit = 's')   
            stamps = start_date + delta 
            df['timestamp'] = stamps
            df['date'] = df['timestamp'].dt.date
            df['hour'] = df['timestamp'].dt.hour
            
            """ Reduce dataframe to sandard plevels 
            #         self.std_plevs    = [1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 85000, 92500, 100000]
            These indices will be used to fill the output array containing the Desroziers errors """
            
            #indices = np.where ( (df['z_coordinate'] == 1000) | (df['z_coordinate'] == 10000) |  (df['z_coordinate'] == 100000 ) )
            
            """ Converting all the hours, that can have random values, to standardized 0,6,12,18 if they are within 2 hours delta """
            df = df.replace({'hour': {1:0 , 2:0, 3: np.nan,  4:6, 5:6, 7:6, 8:6 , 9: np.nan , 
                                      10:12, 11:12, 13:12, 14:12, 15: np.nan, 16:18,17:18,18:18,20:18, 21:np.nan, 22:0, 23:0 }} )
            
            #for plev in [10000,] : #[50000, 92500, 100000]
            for plev in self.std_plevs : 
                
                plev_indices = np.where ( (df['z_coordinate'] == plev ))[0] # first indices relative to pressure level 
                
                
                #index_offest_m = min(plev_indices)
                df_plev = df.iloc[ plev_indices ]
                #df_plev['index'] = np.array(range(len(df_plev['date_time'])))
                
                #total_index_offest += index_offest_m
                #for d in error_day_window: #  ['30','60','90','180']
                for d in error_day_window: #  ['30','60','90','180']
                
                    indices_to_insert = []
                    errors = []
                    num_day = []
                    #stddev = []
                    
                    day_delta = pd.Timedelta( int(d)/2 , unit = 'd' )
                    #delta_plus = df['timestamp'] + day_delta
                    #delta_minus =  df['timestamp'] - day_delta
                    
                    df_plev['day_delta_plus']    = df_plev['timestamp'] + day_delta  # creating columns with the time window given by error_day_window around the observation timestamp
                    df_plev['day_delta_minus'] = df_plev['timestamp'] - day_delta
                    
                    """ Here I select the data from the dataframe which contains only the plevel selected.
                    Form here on I will use only numpy arrays for better speed. """
                    
                    timestamp = np.array(df_plev['timestamp'])
                    df_hour = np.array(df_plev['hour'])
                    day_delta_plus = np.array(df_plev['day_delta_plus'])
                    day_delta_minus = np.array(df_plev['day_delta_minus'])
                    
                    fg_depar = np.array(df_plev['fg_depar@body'])
                    an_depar = np.array(df_plev['an_depar@body'])
                    bias = np.array(df_plev['biascorr@body'])
                                                                          
                    for date,hour, plus, minus, fg_dep, an_dep, bias, num in zip( tqdm(timestamp), 
                                                                                  df_hour, 
                                                                                  day_delta_plus, 
                                                                                  day_delta_minus, 
                                                                                  fg_depar, 
                                                                                  an_depar, 
                                                                                  bias, 
                                                                                  range(len(timestamp)) ) :

                        #if date < np.datetime64('1950-01-01'):
                        #    continue 
                        

                        
                        #selected_indices_t = np.where( (timestamp[num-1000:num+1000] >= minus) &  (timestamp[num-1000:num+1000] <= plus) )[0]
                        #selected_indices_h = np.where ( df_hour[num-1000:num+1000] == hour ) [0] 
                        #selected_indices =  np.intersect1d(selected_indices_t, selected_indices_h)
                        
                        #selected_indices_t = np.where( (timestamp[num-lim:num+lim] >= minus) &  (timestamp[num-lim:num+lim] <= plus) )[0]
                        #selected_indices_h = np.where ( df_hour[num-lim:num+lim] == hour ) [0] 
                        #selected_indices =  np.intersect1d(selected_indices_t, selected_indices_h)      
                        
                        selected_indices_t = np.where( (timestamp>= minus) &  (timestamp <= plus) )[0]
                        selected_indices_h = np.where ( df_hour == hour ) [0] 
                        selected_indices =  np.intersect1d(selected_indices_t, selected_indices_h)                         
                        
                        
                        """ to try 
                        lim = 4* int(d) + 1 # since it is a running wndow, I onyl want days up to date-lim, where lim = 4 (maximum hours available)* number of days of the window 
                        selected_indices_t = np.where( (timestamp[num-lim:num+lim] >= minus) &  (timestamp[num-lim:num+lim] <= plus) )[0]
                        
                        time_index_offest = min(selected_indices_t)
                        total_index_offest += time_index_offest
                        
                        selected_indices_h = np.where ( df_hour.take[selected_indices_t] == hour ) [0] 
                        hour_index_offest = min(selected_indices_h)
                        total_index_offest += hour_index_offest
                        """ 
                        
                        
                        #selected_indices_t = np.searchsorted( (timestamp[num-1000:num+1000] >= minus) &  (timestamp[num-1000:num+1000] <= plus) )[0]
                        #0
                                                
                        # indices correpsonding to selected hour, plevel, aorund a window of days 
                        # at these indices the desrozier values will replace the dummy values in the output array
                        #window_df = df.loc[ selected_indices ]
                        
                        if len(selected_indices)>0:
                            """ Here is the core of the Desroziers error.
                                  Extract the data around a time window,
                                  get the analysis and background departures,
                                  calculate product, sum over the values and extract the root square """
                            # product = df['an_depar@body'][selected_indices] * df['fg_depar@body'][selected_indices]
                            
                            indices_to_insert.append(num)
                            
                            product = an_depar[selected_indices] * fg_depar[selected_indices]
                            valid_data = [p for p in product if not np.isnan(p)]
                            
                            valid_data = remove_outliers(data= valid_data, min_p= 25, max_p= 75, cut= 1, skewed= False, only_clean = True)
                            
                            len_valid = len(valid_data)
                            
                            #if len_valid > int(d)/2:
                            #    print(0)
                            #    0
                                
                            """ If not enough data available (i.e. less than 50% of the dates), will only write the number of data avilable and nan for all the other variable.
                                 Otherwise, procede with the full error estimation.
                                 NB we check if the variable is temperature = 85.
                                 In this case, if bias is not available, we skip the caculation.
                                 For wind variables (speed and direction) we ignore the absence of the bias. """
                            
                            if var == 104 or var == 105 or var == 107:
                                bias = 1 # dummy value to make the following check = False for wind speed, which has no bias
                                
                            if len_valid < int(d)/2 or np.isnan(bias) or np.isnan(fg_dep) or np.isnan(an_dep): # require at least 50% of valid data for the error calculation 
                                desroziers_error = np.nan 
                            
                            else:
                                desroziers_error = np.sqrt(abs(sum(valid_data))/len_valid)
                                #std = np.std(valid_data)
                                
                            errors.append(desroziers_error)
                            num_day.append(len_valid )    
                            #stddev.append(std)
                                
                        else: # no error could be calculated 
                            #errors.append(np.nan)
                            #num_day.append(np.nan )                        
                            #stddev.append(np.nan)
                            """ Everything is nan so nothing to do """
                            indices_to_insert.append(num)      
                            errors.append(np.nan)
                            num_day.append(-2147483647  )                                    
                            pass # nothing to do since the vectors are already initialized with nans 
                        
                    # now we have a list of calculated des. error.
                    # these values will replace the values of the dummy values in the vectors  error, num, std_dev,
                    # at the position indicated by the indices selected_indices
                    
                    if len(indices_to_insert) == 0:
                        continue
                    
                    error, nums = 'desroziers_' + d , 'num_' +d
                    new_indices = imin + indices_to_insert + plev_indices
                    # need to add the starting indices for this variable AND for the selected timestamp
                    np.put(errors_out[error]     , new_indices, errors      ) # replace in the errors_out vectors the values of the vector errors at indices=selected_indices
                    np.put(errors_out[nums]      , new_indices, num_day ) 
                    
                    """
                    valid = np.where( np.array(num_day) > int(d)/2 )
                    if len (valid[0]) > 0:
                        ind = new_indices[valid[0]]
                        values = errors[list(valid[0])]
                        print(0)
                    """
                    
                        
                    #np.put(errors_out[std_dev] , new_indices, stddev   ) 
                    
                    print('Finished: ' , d , ' plevel: ' , plev )
                    
            #print(0)
        
        self.errors_out = errors_out 
        
        
        
    def write_output(self):
        out_data = self.errors_out
        
        encodings = {}
        #val1, val2, val3, val4 = 
        
        errors = xr.Dataset()
        #val30, val90 = out_data['30_desrozier'] , out_data['90_desrozier']
        
        for var in self.errors_out.keys():
            
            if 'desroziers' in var:
                attr = 'Desroziers uncertainty v 1.0 - xx days window'.replace('xx', var.split('_')[-1] )
                encodings[var] = {'dtype': np.float32 , 'compression': 'gzip'}                
            elif 'num' in var:
                attr = 'Number of records - xx days window'.replace('xx', var.split('_')[-1] )
                encodings[var] = {'dtype': int , 'compression': 'gzip'}
                
                
            errors[var] = out_data[var]
            errors[var].attrs['description']= attr
          
        """ Write the new group using h5py """
        group = 'advanced_uncertainty' 
        self.h5py_file.create_group(group)
        #index = np.array( range(len(data['date_time'] ) )  )
        
        
        """ Writing  the uncertainties in the advanced_uncertainty group """
        index = self.h5py_file['observations_table']['index']        
        self.h5py_file[group].create_dataset('index', data=index, compression = 'gzip')
        
        for v in out_data.keys():
            print('*** Writing to output file: ***' , v )
            data =  out_data[v]
            if 'num' in v:
                data = data.astype(int)
            else:
                data = data.astype(float)
                
                
            self.h5py_file[group].create_dataset(v, self.h5py_file['observations_table']['date_time'].shape , encodings[v]['dtype'] , 
                                                 compression= encodings[v]['compression'] , chunks = True , data = data  )
            
            #self.h5py_file[group][v][:] = out_data[v]
            self.h5py_file[group][v].attrs['description'] = errors[v].attrs['description']
            self.h5py_file[group][v].dims[0].attach_scale( self.h5py_file['advanced_uncertainty']['index'] )  
            
            
        """ Adding some extra variables  """    
        
        """
        for v in [ 'date_time' , 'z_coordinate' , 'observation_value']:
            
            self.h5py_file[group].create_dataset(v, self.data[v].shape,  chunks = True, data = self.h5py_file['observations_table'][v] )
            #self.h5py_file[group][v][:] = self.h5py_file['observations_table']['date_time']
            self.h5py_file[group][v].dims[0].attach_scale( self.h5py_file['observations_table']['index'] )  
            
        self.h5py_file[group]['date_time'].attrs['description'] = np.bytes_('seconds since 1900-01-01 00:00:00')
        
        self.h5py_file[group].create_dataset('biascorr@body' , self.data['date_time'].shape,  chunks = True, data =  self.h5py_file['era5fb']['biascorr@body'] )
        #self.h5py_file[group]['biascorr@body'][:] = self.h5py_file['era5fb']['biascorr@body']
        self.h5py_file[group]['biascorr@body'].dims[0].attach_scale( self.h5py_file['advanced_uncertainty']['index'] )          
        """
        
        self.h5py_file.close()
        
        os.system('mv  ' + self.file + '    ' + self.out_dir )

        print(" *** Done writing Desroziers error to output file --->  " , self.file  )
 
if __name__ == '__main__':
        
            parser = argparse.ArgumentParser(description="Utility for Desroziers statistics")

            parser.add_argument('--force_run' , '-f', 
                                  help="Force running the file(s)"  ,
                                  type = str,
                                  default = 'False' )
            
            parser.add_argument('--multi_processing' , '-p', 
                                  help="Running multiprocesses"  ,
                                  type = str,
                                  default = 'False' ) 
            
            args = parser.parse_args()
            
            force_run                = args.force_run
            multiproc                = args.multi_processing
            
            """ Wrapper for running """
            def run(directory, out_dir, force_run, station):
                
                file = directory + '/' + station
              
                station_id = file.split("_CEUAS")[0].split(merged_directory)[1].replace('/','')
                
                """ Initialize classes """
                DS = Desroziers(out_dir = out_dir , station_id = station_id , file = file  )
                """ Read data """
                DS.load_data()

                if force_run in ['yes', 'y', 'YES', 'Y', 'True', 'true']:   # still to implement 
                    print("    === Running in force mode ===     ")
                    
                    try:                        
                        DS.calculate_Desroziers_errors()
                        DS.write_output()      
                        
                    except:
                        print(" The file " + file + ' hase failed! MUST REDO! ')
                        a = open('Failed_Desrozier.txt', 'a')
                        a.write(file + '\n')
                    
                else:
                    DS.calculate_Desroziers_errors()
                    DS.write_output()
            
            """ File source direcotry """
            #merged_directory = '/raid60/scratch/uli/converted_v5/'
            merged_directory = '/raid60/scratch/federico/TO_PROCESS_DESROZIERS'

            """ Moving postprocessed files to new directory """
            out_dir = '/raid60/scratch/federico/DESROZIERS_25MARCH2021'
            
            #os.system('rm -r  ' + out_dir)     
            os.system('mkdir ' + out_dir)

            """
            stations_list = [ s for s in os.listdir(merged_directory) if 'TEST0'  in s ]   
            for s in stations_list:
                os.system('cp ' + merged_directory + '/'+s   + '  '  +  merged_directory + '/'+s.replace('_TEST0.nc','.nc')   )
            
            stations_list = [ s for s in os.listdir(merged_directory) if 'TEST.nc'  in s ]   
            for s in stations_list:
                    os.system('cp ' + merged_directory + '/'+s   + '  '  +  merged_directory + '/'+s.replace('_TEST.nc','.nc')   )
            """        
            stations_list = [ s for s in os.listdir(merged_directory)  ]   
            
            

            cleaned_list = []
            if os.path.isdir(out_dir):
                try:
                    
                    processed = [ s.split('_')[0] for s in os.listdir(out_dir) ]  # skipping processed files 
                except:
                    processed = []
                for file in stations_list:
                        station_id = file.split("_CEUAS")[0]
                        
                        if station_id in processed:
                            print('Already processed:::' , file )
                        else:
                            cleaned_list.append(file)

                
            #cleaned_list = ['0-20000-0-94463_CEUAS_merged_v0.nc', ]
            
            #for s in cleaned_list:
            #    a = run(merged_directory, out_dir, force_run, s)
                
            
            p = Pool(10)
            func = partial(run, merged_directory, out_dir, force_run)
            out = p.map(func, cleaned_list)        

