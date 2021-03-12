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




class Desrozier():
    
    def __init__(self, out_dir = '' , station_id = '' , file='' ):
        self.out_dir = out_dir 
        
        if not os.path.isdir(out_dir):
            print('+++ Creating output directory: ' , out_dir )
            os.mkdir(out_dir)
            
        self.station_id = station_id
        self.file = file 
        #self.summary_file =  station_id + '_' + datetime.now().strftime("%Y-%m-%d_%H-%M-%S") + '_summary_desrozier.txt'
        self.std_plevs    = [1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 85000, 92500, 100000]
        self.variables = [85]
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
        
        os.system('rm 0-20000-0-82930_CEUAS_merged_v0.nc')        
        os.system('cp  0-20000-0-82930_CEUAS_merged_v0_TEST.nc 0-20000-0-82930_CEUAS_merged_v0.nc') # to remove
        
        
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
      
        
    def calculate_Desrozier_errors(self):
        
        """ Errors xarray """
        #errors_out = xr.Dataset()
        errors_out = {}
        
        error_day_window= ['30','60','90','180']
        
        """ Assuming a day window around each observation to calculate the time-average
             e.g. 30 days around 30 january = from 15 Jan to 14 Feb. """
        
        """ Creating the output place-holder filled with np.nan """
        for d in error_day_window: #  ['30','60','90','180']
            error, num, std_dev = d+'_desrozier' , d+'_num' , d+'_std_dev'
            errors_out[error]     = np.full( (len(self.data['date_time'] ) ) , np.nan ) # create place holder for uncertainties filled with nans 
            errors_out[num]      = np.full( (len(self.data['date_time'] ) ) , np.nan ) # create place holder for uncertainties filled with nans 
            errors_out[std_dev] = np.full( (len(self.data['date_time'] ) ) , np.nan ) # create place holder for uncertainties filled with nans 

        for var in self.variables:
            ind = self.recordindices[str(var)][:]
            imin, imax = min(ind) , max(ind) + 1000 # I do not know the maximum, so I take max+1000
                       
            red_dic = {}
            for col in self.data.keys():
                red_dic[col] = self.data[col][imin:imax]
            
            #date_time = start_date + dt
            #delta = np.timedelta64(2244283200, 's')
            #n = start_date + delta 
            
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
            These indices will be used to fill the output array containing the desroziers errors """
            
            #indices = np.where ( (df['z_coordinate'] == 1000) | (df['z_coordinate'] == 10000) |  (df['z_coordinate'] == 100000 ) )
            
            """ Converting all the hours, that can have random values, to standardized 0,6,12,18 if they are within 2 hours delta """
            df = df.replace({'hour': {1:0 , 2:0, 3: np.nan,  4:6, 5:6, 7:6, 8:6 , 9: np.nan , 10:12, 11:12, 13:12, 14:12, 15: np.nan, 16:18,17:18,18:18,20:18, 21:np.nan, 22:0, 23:0 }} )
            
            for plev in [50000, 92500, 100000]: #[1000,50000]
                #plev_indices = np.where ( (df['z_coordinate'] == plev )) # first indices relative to pressure level 
                #df = df.iloc[ plev_indices ]
                
                #for d in error_day_window: #  ['30','60','90','180']
                for d in ['30','90']: #  ['30','60','90','180']
                
                    errors = []
                    num_day = []
                    stddev = []
                    
                    day_delta = pd.Timedelta( int(d)/2 , unit = 'd' )
                    delta_plus = df['timestamp'] + day_delta
                    delta_minus =  df['timestamp'] - day_delta
                    
                    df['day_delta_plus']    = delta_plus
                    df['day_delta_minus'] = delta_minus
                    
                    for date,hour, plus, minus, fg_dep, an_dep in zip( tqdm(df['timestamp'][25000:]), df['hour'],
                                                       df['day_delta_plus'], df['day_delta_minus'], df['fg_depar@body'], df['an_depar@body'] ) :
                        if hour not in [0,6,12,18]:
                            continue 
                        if np.isnan(an_dep) or np.isnan(fg_dep):
                            continue 
                        selected_indices = np.where( (df['timestamp'] >= minus) &  (df['timestamp'] <= plus) 
                                                    & (df['z_coordinate'] == plev )
                                                    & (df['hour']== hour ) )[0] # indices correpsonding to selected hour, plevel, aorund a window of days 
                        # at these indices the desroyier values will replace the dummy values in the output array
                        #window_df = df.loc[ selected_indices ]
                        
                        if len(selected_indices)>0:
                            """ Here is the core of the Desrozier error.
                                  Extract the data around a time window,
                                  get the analysis and background departures,
                                  calculate product, sum over the values and extract the root square """
                            product = df['an_depar@body'][selected_indices] * df['fg_depar@body'][selected_indices]
                            valid_data = [p for p in product if not np.isnan(p)]
                            try:
                                desrozier_error = np.sqrt(abs(sum(valid_data))/len(valid_data)  )
                                errors.append(desrozier_error)
                                num_day.append(len(valid_data) )    
                                stddev.append(np.std(valid_data))
                            except:
                                pass
                        else: # no error could be calculated 
                            #errors.append(np.nan)
                            #num_day.append(np.nan )                        
                            #std_dev.append(np.nan)
                            pass # nothing to do since the vectors are already initialized with nans 
                        
                    # now we have a list of calculated des. error.
                    # these values will replace the values of the dummy values in the vectors  error, num, std_dev,
                    # at the position indicated by the indices selected_indices
                    error, num, std_dev = d+'_desrozier' , d+'_num' , d+'_std_dev'
                    np.put(errors_out[error]     , selected_indices, errors      ) # replace in the errors_out vectors the values of the vector errors at indices=selected_indices
                    np.put(errors_out[num]      , selected_indices, num_day ) 
                    np.put(errors_out[std_dev] , selected_indices, stddev   ) 
                    print('Finished: ' , d , ' plevel: ' , plev )
                    
            print(0)
        
        self.errors_out = errors_out 
        
        
        
    def write_output(self):
        out_data = self.errors_out
        
        encodings = {}
        #val1, val2, val3, val4 = 
        
        errors = xr.Dataset()
        #val30, val90 = out_data['30_desrozier'] , out_data['90_desrozier']
        
        for var,attr in zip(['30_desrozier','90_desrozier', '30_std_dev', '30_num'] , 
                            ['Desrozier uncertainity-30 days average' , 'Desrozier uncertainity- 90 days average',
                             'Standard deviation over 30 days window', 'Number of records over 30 days window'] ):
            errors[var] = out_data[var]
            errors[var].attrs['description']= attr
            encodings[var] = {'dtype': np.float32 , 'compression': 'gzip'}
          
        """ Write the new group using h5py """
        group = 'advanced_uncertainty'
        self.h5py_file.create_group(group)
        #index = np.array( range(len(data['date_time'] ) )  )
        

        for v in ['30_desrozier', '90_desrozier' , '30_std_dev', '30_num']:
            print('*** Writing to output file: ***' , v )
            self.h5py_file[group].create_dataset(v , errors[v].shape, encodings[v]['dtype'] , compression= encodings[v]['compression'] , chunks = True )
            self.h5py_file[group][v][:] = out_data[v]
            self.h5py_file[group][v].attrs['description'] = errors[v].attrs['description']
            
        self.h5py_file[group].create_dataset('date_time' , self.data['date_time'].shape,  chunks = True )
        
        self.h5py_file.close()

        print(" *** Done writing desrozier's error to output file --->  " , self.file  )

 
if __name__ == '__main__':
        
            parser = argparse.ArgumentParser(description="Utility for Desrozier statistics")

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
                DS = Desrozier(out_dir = out_dir , station_id = station_id , file = file  )
                """ Read data """
                DS.load_data()

                if force_run in ['yes', 'y', 'YES', 'Y', 'True', 'true']:   # still to implement 
                    print("    === Running in force mode ===     ")
                    
                    try:                        
                        DS.calculate_Desrozier_errors()
                        DS.write_output()      
                        
                    except:
                        print(" The file " + file + ' hase failed! MUST REDO! ')
                        a = open('Failed_Desrozier.txt', 'a')
                        a.write(file + '\n')
                    
                else:
                    DS.calculate_Desrozier_errors()
                    DS.write_output()
            
            """ File source direcotry """
            #merged_directory = '/raid60/scratch/uli/converted_v5/'
            merged_directory = os.getcwd()

            """ Moving postprocessed files to new directory """
            out_dir = '/raid60/scratch/federico/DESROZIERS_MARCH2021'
            os.system('mkdir ' + out_dir)

            stations_list = [ s for s in os.listdir(merged_directory) if 'empty'  not in s ]           
            processed = [ s.split('_')[0] for s in os.listdir(out_dir) ]  # skipping processed files 

            cleaned_list = []

            stations_list = ['0-20000-0-82930_CEUAS_merged_v0.nc']
            for file in stations_list:
                station_id = file.split("_CEUAS")[0]
                
                if station_id in processed:
                    print('Already processed:::' , file )
                else:
                    cleaned_list.append(file)
            
            
            #cleaned_list = ['0-20000-0-94463_CEUAS_merged_v0.nc', ]
            
            for s in cleaned_list:
                a = run(merged_directory, out_dir, force_run, s)
                
            
            #p = Pool(30)
            #func = partial(run, merged_directory, postprocessed_new, force_run)
            #out = p.map(func, cleaned_list)        

