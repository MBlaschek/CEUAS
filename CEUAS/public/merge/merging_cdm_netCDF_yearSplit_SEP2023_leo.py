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
import time 
from datetime import datetime
import shutil

from tqdm import tqdm

import matplotlib.pyplot as plt

#from numba import njit
import psutil
import copy
from numba import njit
from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
import warnings

from multiprocessing import Pool
from functools import partial

import ray

warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)

import code

sys.path.append('../harvest/code_cop2/')
sys.path.append('../postprocess/add_sensor/')

import plotly.graph_objects as go
from plotly.subplots import make_subplots


from add_sensor_to_merged_OCT2023 import Sensor, datetime_toseconds, wrapper, MergedFile
from harvest_convert_to_netCDF_yearSplit import  clean_station_configuration , write_dict_h5_old
#from harvest_convert_to_netCDF_yearSplit import   write_dict_h5

# nan int = -2147483648 
#from harvest_convert_to_netCDF import datetime_toseconds   # importing the function to write files with h5py 


"""
@njit
def replace_global(replaced_indices, replaced_vector, replacing_vector):
    #FAILED Attempt to make vector item replacement faster 
    for insert_index, replacing_value in zip(replaced_indices, replacing_vector):
        replaced_vector[insert_index] = replacing_value 
"""    

        
"""
PRIORITY
- ERA5 over all IF LEN >=
- GIUB over ERA5 pre 1950
- IGRA over ERA5 if max(plevel_igra) > max(plevel_era5)
- longest record if ERA5 and IGRA not available
- what to do with NCAR => set as lowest priority, only select if only one available
"""








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

if not os.path.isdir('logs'):
    os.mkdir('logs')
    
def now(time):
    """ For easily print a readable current time stamp """
    a = datetime.fromtimestamp(time).strftime('%Y-%m-%d %H:%M:%S')
    return  a



    
    
    
data = {}
class Merger():
    """ Main class for the merging of the data from different netCDF files. """

    def __init__(self, 
                 add_sensor = '',  
                 out_dir = 'output' , 
                 min_year='1880' , 
                 max_year='2023',
                 processed_stats =''):
        
        """ Define the attributes (some will be defined in other parts of the code) . 
        Attr :: 
                self.data : read the dictionary mapping the dataset and the netCDF cdm file for each observation station  
                self.datasets : store the
                self.datasets_keys : ''
                self.datasets_all : hard coded list of all the potentially available datasets. Nothe that the ncar is split according to wind (_w) and temperature (_t) variables as in the original dataset source
                self.observation_ids_merged  : dictionary used to calculate the new merdeg record_id and report_id 
                self.unique_dates : dictionary containing the lists of observation date_times and indices for each dataset and of the merged dataset             
        """

        self.data = {}                              # will contain the data for each different dataset 
        self.datasets = ''                          # will contain the input datasets (original dictionary)
        self.datasets_keys = ''                 # will contain the input datasets names only (i.e. keys of the datasets dictionary)
        #self.datasets_all = ['era5_2_2']    # all possibly available datasets                          

        self.unique_dates = {}            
        self.attributes = {} # will keep the original attributes from the CDM tables, read from the netCDF files 
        self.id_string_length = 14 # fixed length for record_id and observation_id values 
        self.out_dir = out_dir 
        
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        if not os.path.exists(out_dir+'/logs'):
                os.makedirs(out_dir+'/logs')            
        
        self.variable_types = {}
        self.processed_stats = processed_stats 
        '''
        self.observation_ids_merged  = {  'igra2':b'3' , 
                                          'ncar':b'4', 
                                          'bufr':b'5',  
                                          'era5_1':b'1' , 
                                          'era5_1_mobile':b'1' , 
                                          
                                          'era5_2':b'2', 
                                          'era5_2_mobile':b'2', 
                                          
                                          'era5_1759' :b'6' , 
                                          'era5_1761':b'7' ,  
                                          'era5_3188' :b'8' ,
                                          'amma': b'9' }
        '''
        
        # values used to convert original record_id to the merged record_id, see method merge_all_data 

        logging.info('*** Initialising the Merging procedure ***' )   
        #self.era5b_columns = []  # stores the columns of the era5fb 
        self.standard_cdm = [ 'crs' , 'observed_variable', 'units' , 'z_coordinate_type' , 'station_type', 'station_configuration_codes'] 

        self.hour_time_delta = 2 # decide up to which time shift in HOURS separate records are considered identical  

        self.only_std_plevels = False  # set to True to store only standard pressure level data 
        self.std_plevs    = [1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 85000, 92500, 100000]

        self.add_sensor = add_sensor
            
        self.copy = True # make a copy of the merged file before adding the sensor. Argument for the add_sensor wrapper function  


    def initialize_data(self , datasets, station ):
        """ Initialize dataset; store relevant data as attributes.
                   Args ::     dic{}  datasets (dataframe with year, dataset and file path 
                                   e.g. tateno = { 'ncar'    : 'example_stations/ncar/chuadb_trhc_47646.txt.nc'   ,
                                                           'igra2'   : 'example_stations/igra2/chJAM00047646-data.txt.nc'  } 
        """       
        self.datasets          = datasets
        self.datasets_keys = datasets.keys()
        self.station             = station
        
        self.observed_ts_values= {} ### will store observation values for temp and wind for each current timestamp 
        
        # list of meaningful variables in the observations_table, all others are empty/nans filled
        self.observations_table_vars = ['date_time', 'z_coordinate' , 'z_coordinate_type', 'observed_variable' , 
                                                           'observation_value' , 'report_id' , 'observation_id' , 'latitude' , 'longitude',
                                                           'units', 'source_id', 'data_policy_licence', 'observation_duration', 'value_significance',
                                                           'advanced_assimilation_feedback' , 'advanced_uncertainty']
                
        """ Loading the econding of the tables created from the harvester script and to be applied again """
        
        shutil.copy2('../harvest/code_cop2/encodings.txt', 'encodings.txt')
        df = pd.read_csv('encodings.txt' , sep = '\t' , names = ['variable' , 'table'  , 'type'])
        self.df_variable_type = df 
        
        #shutil.copy2('../harvest/code_cop2/groups_encodings.npy', 'groups_encodings.npy')
        self.encodings = np.load('groups_encodings.npy' , allow_pickle = True ).item()
        #shutil.copy2('../harvest/code_cop2/era5fb_encodings_all.npy', 'era5fb_encodings_all.npy')        
        self.encodings['era5fb'] = np.load('era5fb_encodings_all.npy' , allow_pickle = True ).item()
        for k in self.encodings.keys():
            for kk in self.encodings[k].keys():
                
                if 'dtype' in self.encodings[k][kk].keys():
                    del self.encodings[k][kk]['dtype']
                if 'chunksizes' in self.encodings[k][kk].keys():
                    del self.encodings[k][kk]['chunksizes']
                if 'compression' in self.encodings[k][kk].keys():
                    self.encodings[k][kk]['compression'] = 32015
                    self.encodings[k][kk]['compression_opts'] = (3, )

        
        self.dic_type_attributes = np.load('dic_type_attributes.npy',allow_pickle= True).item()

        self.era5fb_columns = self.dic_type_attributes['era5fb'].keys()
        self.header_columns = self.dic_type_attributes['header_table'].keys()

        self.empty_cdm_var = [ v for v in self.dic_type_attributes['observations_table'].keys() if v not in self.observations_table_vars ]  # variables to be filled with nans   with proper data type         


        self.obstab_nans_filled = False      
        self.obs_in_header = {} # dic placeholder for data from obs_tab to be kept for header_table
        self.fill_cdm()
        
        self.last_timestamp = ''  # to be filled progressively
        
    def fill_cdm(self):
            
        self.data['cdm_tables'] = {}              

        """ Loop over all the datasets                                                                                                                                     
                k: name of the dataset                                                                                                                    
                v: list of file paths,    eg 'era5_1':[filepath_1, filepath_2 ] """      
        
        #######                                           
        # STATION CONFIGURATION
        #######   
        # Read the station primary_id, and find the entry in the CUON dataset (or global dataset)
        stat_conf = pd.read_csv("CUON_station_configuration_extended.csv", sep = '\t')
        self.stat_conf_CUON = stat_conf 
        s = stat_conf.loc[ stat_conf.primary_id == self.station ]   
        self.stat_conf = s
        
        if s.empty:
            a = open('logs/failed_stat_conf.txt' , 'a+')
            a.write('Cannot_find_station_configuration_for_' + self.station + '\n')
        else:
            for k in s.columns:
                if pd.isnull(s[k].values[0]):
                    try:
                        v = s[k].astype(np.float64)            
                    except:
                        v = b'NA'           
                    s[k] = v 
                    
        if '20999' in self.station:
            sc = pd.DataFrame( 0, index=[0], columns = stat_conf.columns)  # empty dataframe 
        else:
            #sc = clean_station_configuration(s)  ### FIX THIS, dont knwo why it happens 
            sc = s 
        
        # dropping false columns
        for c in sc.columns:
            if 'Unnamed' in c:
                sc = sc.drop(columns=c)
            
        self.data['station_configuration'] = sc
        self.statconf_columns = sc.keys()

        """ Reading the CDM tables that do not depend on specific stations or observations (fixed values), for the first file only """
        for t in self.standard_cdm: # [ 'crs' , 'observed_variable', 'units' , 'z_coordinate_type' , 'station_type', 'station_configuration_codes']   
            try:
                if t not in self.data['cdm_tables'].keys():
                    #data['cdm_tables'][t] = ''
                    cdm = xr.open_dataset('example_standard_cdm_tables.nc'  , engine = 'h5netcdf' , group = t )
                    self.data['cdm_tables'][t] = cdm 
            except:
                print("FAIL standard cdm group " , t )
                pass
                
        
    def extract_file_per_year(self):
        """ Extract each unique year, and files for each dataset """
        years = np.unique( self.datasets.year)
        
        all_data = {}
        
        df = self.datasets
        for y in years:
            all_data[y] = {}
            
            df_r = df.loc[df.year == y ]
            datasets = np.unique(df_r.dataset)
            for ds in datasets:
                df_r_ds = df_r.loc[df_r.dataset == ds ]
                all_data[y][ds] = list(df_r_ds.files)
            
        return all_data
        
        
    def get_null(self, tipo = ''):
        """ Simply returns the proper format for ''null' value """        
        if tipo in  [np.int32, np.int64]  :
            void = -2147483648
        elif tipo == np.float32 :
            void = np.nan
        elif tipo == np.bytes_ :
            void = b'nan'
        else:
            return np.nan 
        return void
        
        
    def open_data(self, data):
        """ Open each file with hyp5, store in dictionary
        
        data=
        {'era5_1759': ['/scratch/das/federico/HARVEST_YEARLY_17AUG2023/era5_1759/0-20000-0-71879_1962_era5_1759_harvested_era5.1759.conv._1:71879.gz.nc'],
        'ncar': ['scratch/das/federico/HARVEST_YEARLY_17AUG2023/ncar/0-20000-0-71879_1962_ncar_harvested_uadb_trhc_71879.txt.nc'],
        }
        """
        
        dic = {}
        for ds in data.keys():
            dic[ds]={}
            for file in data[ds]:
                
                files = file.split('/')[-1].split(',')  # there might be two different files for the same dataset and year, for example windc and trhc NCAR 
                if len(files) > 1:
                    for f in files:
                        fpath = file.split(self.station)[0] + '/' + self.station + '/' + files[0]   # e.g. '/scratch/das/federico/HARVEST_YEARLY_18SEP2023//ncar//0-20001-0-11035/0-20001-0-11035_1951_ncar_harvested_uadb_trhc_11035.txt.nc'
 
                        dic[ds][fpath]  = h5py.File(fpath, 'r')
                else:
                    
                    dic[ds][file]  = h5py.File(file, 'r')
                    
        self.dic_h5py = dic
        
        
        
    def make_unique_datetime(self):
        """ Building the global set of date_times and indices from the various datasets. 
              The datetimeindex is read from the original netCDF file. 
              Will compare the unique date_time of each dataset and extract the global unique date_times
              
              Return a dictionary with the dataset as key, then file and 
              {1966028400 : { 'era5_1759': {'/scratch/das/federico/HARVEST_YEARLY_17AUG2023/era5_1759/0-20000-0-71879/0-20000-0-71879_1962_era5_1759_harvested_era5.1759.conv._1:71879.gz.nc': [1966028400, 0, 90]}, 
              'ncar': {'/scratch/das/federico/HARVEST_YEARLY_17AUG2023/ncar/0-20000-0-71879/0-20000-0-71879_1962_ncar_harvested_uadb_trhc_71879.txt.nc': [1966028400, 598, 235]}} }
        """

        logging.info('\n *** Running make_all_datetime ' )

        """ Loop over all the datasets 
                k: name of the dataset
                v: list of file paths,    eg 'era5_1':[filepath_1, filepath_2 ]"""

        all_timestamps_dic = {}
        
        # containing all files and h5py open file (per year)
        dic= self.dic_h5py
        for ds,files in dic.items() :
            k,v = ds, files #rename
            self.unique_dates[k] = {}  # unique dates for each file per dataset 

            for F in v:                 
                #print('FILE ' , F )
                data = dic[k][F]  # h5py open file 
                self.unique_dates[k][F] = {}
                self.unique_dates[k][F]['datetimes'] = {}
                
                timestamps = data['recordtimestamp'][:]
                #indices = data['recordindex'][:]
                indices_inf = data['recordindex'][:] # starting index of the record
                indices_sup = np.concatenate((data['recordindex'][:][1:], [data['observations_table']['index'].shape[0]]))
                
                #data.[ indices_inf[1:][i] for i in  range(len(indices_inf[:-1])) ] # ending index of the record
                #indices_sup.append(len(h5py.File(F, 'r')['observations_table']['date_time'] ) ) # dummy high value, to be used with last record of the list 
                
                ### quick check
                # df_check = pd.DataFrame( {'ts':timestamps , 'inf': indices_inf, 'sup':indices_sup } ) 
                
                for ts,inf,sup, index in zip(timestamps, indices_inf, indices_sup, range(len(indices_inf)) ) :
                    if ts not in all_timestamps_dic.keys():
                        all_timestamps_dic[ts] = {}
                    if ds not in all_timestamps_dic[ts].keys():
                        all_timestamps_dic[ts][ds] = {}   
                              
                    all_timestamps_dic[ts][ds][F] = [inf, sup, index]  # store the records limits in the obesrvation_table and the record_index in the headet_table

                #timestamps_conv = pd.to_datetime( timestamps, unit='s',  origin=pd.Timestamp('1900-01-01') )
                 
        #unique_timestamps = list(np.unique(list(all_timestamps_dic.keys() ))).sort()
        #if hasattr(self, 'all_timestamps_dic') :
            
            #self.all_timestamps_dic = self.all_timestamps_dic | all_timestamps_dic
        #else:
        self.all_timestamps_dic = all_timestamps_dic
            
                    
    
    def reduce_timestamps(self):
        """ Simplify/reduce all timestamps by flagging possible duplicates """
        
        #1 check duplicates 
        time_delta = self.hour_time_delta * 60*60 # timestamps are still in seconds, so must convert hours in self.hour_time_delta to seconds 
        unique_timestamps = list(np.unique(list(self.all_timestamps_dic.keys() )))
        unique_timestamps.sort()
        
        #d = 1988064000 # 658335600
        
        if self.last_timestamp: # this is the last timestamp from previous year i.e. previous file 
            if unique_timestamps[0] - self.last_timestamp < time_delta:
                # remove first entry from new timestamps -> easy way to remove duplicate 
                unique_timestamps = unique_timestamps[1:]
        else:
            pass
        
        duplicated_ts = []
        unique_ts = [] 
        
        if len(unique_timestamps) ==1: #only one timestamp 
            unique_ts.append(unique_timestamps[0] )
            
        else:
            for i in range(1, len(unique_timestamps)):
                current = unique_timestamps[i] 
                previous = unique_timestamps[i-1]
                
                a= current - previous < time_delta 
                #print('index: ' , i, 'previous: ' , previous, 'current: ' , current  ,  '  ' , str(a) + ' is duplicate')
                
                if (current - previous) < time_delta :  #here: found duplicated 
                    #print('Must check duplicated timestamps', current , previous, current - previous )
                    if len(duplicated_ts) > 0: # check for triple duplicates
                        if previous in duplicated_ts[-1]:
                            duplicated_ts[-1].append(current)
                        else:
                            duplicated_ts.append( [previous, current] )
                    else:
                        duplicated_ts.append( [previous, current] )

                else:
                    if i==len(unique_timestamps)-1:
                         # last record, will not be looped over again
                        unique_ts.append(previous)
                        unique_ts.append(current)
                        
                    else:
                        if len(duplicated_ts) > 0:
                            if previous not in duplicated_ts[-1]:
                                unique_ts.append(previous)
                        else:
                            unique_ts.append(previous)                        

        return unique_ts, duplicated_ts    
    
    def reduce_timestamps_u(self):
        """ Simplify/reduce all timestamps by flagging possible duplicates """
        
        #1 check duplicates 
        time_delta = self.hour_time_delta * 60*60 # timestamps are still in seconds, so must convert hours in self.hour_time_delta to seconds 
        unique_timestamps = np.unique(list(self.all_timestamps_dic.keys() ))
        
        duplicated_ts = []
        for i in range(0, len(unique_timestamps)):
            duplicated_ts.append([])
            for k in range(-5, 6):
                if i + k >= 0 and i + k < unique_timestamps.shape[0]:
                    if np.abs(unique_timestamps[i+k] - unique_timestamps[i]) < time_delta:
                        duplicated_ts[-1].append(unique_timestamps[i+k])
                
        
        #d = 1988064000 # 658335600
        
        #if self.last_timestamp: # this is the last timestamp from previous year i.e. previous file 
            #if unique_timestamps[0] - self.last_timestamp < time_delta:
                ## remove first entry from new timestamps -> easy way to remove duplicate 
                #duplicated_ts[0].append(self.last_timestamp)
        #else:
            #pass
        
        unique_ts = []

        return unique_ts, duplicated_ts, list(unique_timestamps)   
    

    def reduce_timestamps_u_old(self):
        """ Simplify/reduce all timestamps by flagging possible duplicates """
        
        #1 check duplicates 
        time_delta = self.hour_time_delta * 60*60 # timestamps are still in seconds, so must convert hours in self.hour_time_delta to seconds 
        unique_timestamps = list(np.unique(list(self.all_timestamps_dic.keys() )))
        unique_timestamps.sort()
        
        #d = 1988064000 # 658335600
        
        #if self.last_timestamp: # this is the last timestamp from previous year i.e. previous file 
            #if unique_timestamps[0] - self.last_timestamp < time_delta:
                ## remove first entry from new timestamps -> easy way to remove duplicate 
                #unique_timestamps = unique_timestamps[1:]
        #else:
            #pass
        
        duplicated_ts = []
        unique_ts = [] 
        
        if len(unique_timestamps) ==1: #only one timestamp 
            unique_ts.append(unique_timestamps[0] )
            
        else:
            for i in range(1, len(unique_timestamps)):
                current = unique_timestamps[i] 
                previous = unique_timestamps[i-1]
                
                a= current - previous < time_delta 
                #print('index: ' , i, 'previous: ' , previous, 'current: ' , current  ,  '  ' , str(a) + ' is duplicate')
                
                if (current - previous) < time_delta :  #here: found duplicated 
                    #print('Must check duplicated timestamps', current , previous, current - previous )
                    if len(duplicated_ts) > 0: # check for triple duplicates
                        if current - duplicated_ts[-1][0] < time_delta:
                            if previous in duplicated_ts[-1]:
                                duplicated_ts[-1].append(current)
                            else:
                                duplicated_ts.append( [previous, current] )
                        else:
                            duplicated_ts.append( [previous, current] )
                    else:
                        duplicated_ts.append( [previous, current] )

                else:
                    if i==len(unique_timestamps)-1:
                         # last record, will not be looped over again
                        unique_ts.append(previous)
                        unique_ts.append(current)
                        
                    else:
                        if len(duplicated_ts) > 0:
                            if previous not in duplicated_ts[-1]:
                                unique_ts.append(previous)
                        else:
                            unique_ts.append(previous)                        

        return unique_ts, duplicated_ts, unique_timestamps    
    
    
    def extract_record_data(self, dt, ds, file ):
        
        """ Extracting the length of the temp and wind valid observations, and the maximum height (or min pressure) 
        TODO must check what happens with z_coordinate type in height or gph """
        
        ###TO DO CAN LOAD OBS TAB adnd HEAD TAB in memory so you only slice and not read everytime 
        
        #if dt == 2114769600: #2115331200
        #    a = 0 
            
        h5_file = self.dic_h5py[ds][file]
        ind_min, ind_max = self.all_timestamps_dic[dt][ds][file][0] ,  self.all_timestamps_dic[dt][ds][file][1] 
        
        ot = h5_file['observations_table']
        
        otzt = ot['z_coordinate_type'][ind_min:ind_max]
        ### Data on Pressure 
        try:
            #pressure_ind = np.where(otzt == 1 )[0] #
            pim = otzt == 1
        except:
            print("CHEEEECK FILE " , file  )
        ### Data on Height 
        #height_ind = np.where(otzt == 0 )[0] #
        him = otzt == 0
        
        #if len(pressure_ind) >0:
        if np.any(pim):
            #z_ind = pressure_ind
            zim = pim
            z = 'pressure'
        else:
            #z_ind = height_ind
            zim = him
            z= 'height'


        # subset of indices with pressure / height
        otov = ot['observed_variable'][ind_min:ind_max]
        #temp_ind = np.where(otov == 126)[0] # temperature
        #temp_ind = [i for i in temp_ind if i in z_ind ]
        temp_ind = np.where((otov==126)&zim)[0]
        
        #wind_ind = np.where(otov  == 107)[0] # wind speed 
        #wind_ind = [i for i in wind_ind if i in z_ind ]
        wind_ind = np.where((otov==107)&zim)[0]
        
        #gph_ind = np.where(otov  == 117)[0] # geopotential 
        #gph_ind = [i for i in gph_ind if i in z_ind ]
        gph_ind = np.where((otov==117)&zim)[0]

        
        # length of valid data for temp and wind
        otoval = ot['observation_value'][ind_min:ind_max]
        temp_values = otoval[temp_ind]
        num_valid_temp = len(np.unique(temp_values[~np.isnan(temp_values)]))
        
        wind_values = otoval[wind_ind]
        num_valid_wind = len(np.unique(wind_values[~np.isnan(wind_values)]))
        
        otz = ot['z_coordinate'][ind_min:ind_max]
        values_dict = {'temp' :  temp_values , 
                                'temp_z' : otz[temp_ind] ,
                                'z_type' : z ,
                                'wspeed' : wind_values , 
                                'wind_z':  otz[wind_ind] } 
        
        self.observed_ts_values[file] = values_dict
        
        #checking pressure or height or gph 
        #press_temp_ind = ot['z_coordinate'][ind_min:ind_max][temp_ind]
        #press_wind_ind = ot['z_coordinate'][ind_min:ind_max][wind_ind]
        
            
        if len(temp_ind) >0:
            min_temp = min( otz[temp_ind] )  
            max_temp =  max( otz[temp_ind] ) 
        else:
            min_temp, max_temp = 999999, -999999 
            
        if len(wind_ind) >0:
            min_wind =  min( otz[wind_ind] ) 
            max_wind =  max( otz[wind_ind] )
            
        else:
            min_wind, max_wind = 999999, -999999        


        ### min, MAX pressure
        #if len(pressure_ind) >0:
        if np.any(pim):
            
            if len(temp_ind) > 0:
                
                max_pressure_t = max( otz[temp_ind] )
                min_pressure_t = min( otz[temp_ind] )
            else:
                max_pressure_t = 0
                min_pressure_t = 999999
                
            if len(wind_ind) > 0:
                max_pressure_w = max( otz[wind_ind] )
                min_pressure_w = min( otz[wind_ind] )
            else:
                min_pressure_w = 999999
                max_pressure_w = 0                
                
            max_pressure = max( max_pressure_t, max_pressure_w )
            min_pressure = min( min_pressure_t, min_pressure_w )
            
            
        else:
            max_pressure = -999999
            min_pressure = 999999
            
        ### min, MAX height
        #if len(height_ind) >0:
        if np.any(him):
            max_height = max( otz[him] )
            min_height = min( otz[him] )
            
        else:
            max_height = -999999
            min_height = 999999
            
        if len(gph_ind) >0:
            max_gph = max( otoval[gph_ind] )
            min_gph = min( otoval[gph_ind] )
            
        else:
            max_gph = -999999                    
            min_gph = 999999
            
        return [num_valid_temp, num_valid_wind], [min_temp, min_wind, max_temp, max_wind], [min_pressure, min_gph, min_height], [max_pressure, max_gph, max_height]
    
    
        
    
        
 
    def plot_profile_extended(self, all_times, real_time = '', best_file=''):
        """ Create a simple plot of the temperature and wind profile for a specific timestamp (including the duplicated possibilities).
        Save a png file with matplolib and an interactive HTML file from plotly """
        
        # output directory 
        plot_dir = self.out_dir + '/' + self.station + '/merging_plots/' 
        if not os.path.isdir(plot_dir):
            os.makedirs(plot_dir) 
            
        ### MATPLOTLIB
        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(15,8))
        
        for dt in all_times:
            dic_dataset = self.all_timestamps_dic[dt] 
            ds = [d for d in dic_dataset if 'data' not in d and 'record' not in d ]
            
            all_files = []
            for d in ds:
                all_files.extend ( list(dic_dataset[d].keys()) ) 
                
            date = pd.to_datetime( dt, unit='s',  origin=pd.Timestamp('1900-01-01') )
            data= self.observed_ts_values
                
            for file in list(all_files):
                label = str(date) + ' / ' + file.split('/')[-1]
                if file == best_file:
                    s = 4
                    label = 'BEST ' + label
                else:
                    s = 2
                try:
                    ax[0].plot(  data[file]['temp']     , data[file]['temp_z'] , label = label + '  [' + str(min(data[file]['temp_z'] ) ) + ']' , lw=s ) #row=0, col=0
                except:
                    pass
                try:
                    ax[1].plot(  data[file]['wspeed'] , data[file]['wind_z'] , label = label + '  [' + str(min(data[file]['wind_z'] ) ) + ']' , lw=s ) #row=0, col=0
                except:
                    pass
        ax[0].grid(color='lightgray' , ls=':')
        ax[1].grid(color='lightgray' , ls=':')        
        if not ax[0].get_legend_handles_labels() == ([], []):
            ax[0].legend(fontsize=8)
        if not ax[1].get_legend_handles_labels() == ([], []):
            ax[1].legend(fontsize=8)
    
        #plt.show()

        plt.savefig(plot_dir + str(date) + '_profiles.png' , dpi=150)
        plt.close()
                
        ### PLOTLY

        fig = make_subplots(rows=1, cols=3,
                            column_widths=[0.20, 0.20, 0.60], 
                            specs=[[{"type": "scatter"},{"type": "scatter"},{"type": "table"}]],
                            subplot_titles=("Temperature [K]","Wind Speed [m/s]", "Summary")
                            )
        
        dic_data = {'Date' : [],
                            'File' : [], 
                            'Dataset' : [],
                            'Min P Temp' : [] ,
                            'Max P Temp' : [] ,
                            'Min P Wind' : [] ,
                            'Max P Wind' : [] ,
                            'Best' : []
                             }
                
        all_datasets =[]
            
        for dt in all_times:
            dic_dataset = self.all_timestamps_dic[dt] 
            ds = [d for d in dic_dataset if 'data' not in d and 'record' not in d ]
            
            all_files = []
            for d in ds:
                all_files.extend ( list(dic_dataset[d].keys()) ) 
                all_datasets.extend( [d] * len( list(dic_dataset[d].keys()) )  )
                
            date = pd.to_datetime( dt, unit='s',  origin=pd.Timestamp('1900-01-01') )
            
            data= self.observed_ts_values
                
            for file in list(all_files):
                isBest = '-'
                label = file.split('/')[-1].replace('harvested_','').replace('.txt.nc','').replace(self.station, '').replace('_'+self.current_year , self.current_year )
                
                if file == best_file:
                    label_v = 'BEST ' + label
                    isBest = 'BEST'
                else:
                    label_v = label
                    
                fig.add_trace(go.Scatter(x=data[file]['temp'] , y=data[file]['temp_z'],
                                    mode='lines',
                                    name=label_v, ),
                                row=1, col=1,
                                    )
                
                fig.add_trace(go.Scatter(x=data[file]['wspeed'] , y=data[file]['wind_z'],
                                    mode='lines',
                                    name=label_v ),
                              row=1, col=2,
                              )
                
                ### preparing dataframe for table 
                dic_data['Date'].append( date )
                dic_data['File'].append( label )
                dic_data['Best'].append( isBest )
                
                try:
                    dic_data['Min P Temp'].append(  '%.2f'%(min(data[file]['temp_z'] )  ) ) 
                    dic_data['Max P Temp'].append( '%.2f'%(max(data[file]['temp_z'] ) ) )
                except:
                    dic_data['Min P Temp'].append( 'NA' ) 
                    dic_data['Max P Temp'].append( 'NA' )
                try:
                    dic_data['Min P Wind'].append( '%.2f'%(min(data[file]['wind_z']) ) )
                    dic_data['Max P Wind'].append('%.2f'%(max(data[file]['wind_z']) ) )
                except:
                    dic_data['Min P Wind'].append( 'NA' )
                    dic_data['Max P Wind'].append( 'NA' )                    
                
            dic_data['Dataset'] = list(np.array(all_datasets).flatten() )
            df = pd.DataFrame.from_dict(dic_data )
              
            fig.add_trace( go.Table(
                
                columnwidth = [250,480,100,150,150,150,150,100],
                
                header=dict(values=list(df.columns),
                            fill_color='paleturquoise',
                            align='left'),
                cells=dict(values = [ df[c] for c in df.columns ],
                           fill_color='lavender',
                           align='left')),
                row=1, col=3,
                           )

        fig.update_layout( title= self.station + '  -  Profile for Record @ ' + str(date) + ' ' + str(dt),
                                 width=2200,
                                 height=1100,                             
                                 )
        fig.update_yaxes( autorange="reversed",
                    row=1, col=1
                             )
        fig.update_yaxes( autorange="reversed",
                    row=1, col=2
                             )        

        fig.update_yaxes(tickformat="f" , row=1, col=1 )
        fig.update_yaxes(tickformat="f" , row=1, col=2 )
        
        fig.update_layout(legend=dict(
            yanchor="bottom",
            y=-0.25,
            xanchor="left",
            x=0.01
        ))
        
        fig.write_html(self.out_dir + '/' + self.station + '/merging_plots/' + str(date).replace(' ','--')  + '_profiles.html')
        
        plt.close()
        
    
        
    def find_best_record(self, ts):
        """ Effective merging procedure, applying hierarchical selection rules 
                   ts: timestamp in seconds after 1900 
                   datasets: list of datasets per record e.g. ['era5_1', 'igra2' , 'ncar', ... ] """
        
        ### Variable placeholders
        best_ds = False
        best_file = False

        # alla available datasets
        datasets = [ l for l in list(self.all_timestamps_dic[ts].keys() ) if 'extract_record_data' not in l and 'best_record' not in l ]   # might remove the datasets as a variable since it is in an attirbute of the class already. Think about it....
        
        ### Setting data policy (might be redundant)
        ### see https://github.com/glamod/common_data_model/blob/master/tables/data_policy_licence.dat        
        if 'ncar' in datasets or  'igra2' in datasets or  'era5_1759' in datasets or 'era5_1761' in datasets or 'giub' in datasets:
            data_policy = 0
        else:
            data_policy = 4        
        
        # list of all era5 datasets to be preferred except special cases 
        era5_ds = [ f for f in datasets if f in ['era5_1', 'era5_1_mobile' , 'era5_2' , 'era5_2_mobile'] ]

        '''
        ###  1) giub over era5 if year < 1950
        ###  1950-01-01-00:00:00 = 1577836800 seconds
        ###  pd.to_datetime( a[:], unit='s',  origin=pd.Timestamp('1900-01-01') )
        '''
        
        ### placeholder for best data 
        extracted_best_record = {}
        

        if ts == 2492985600:
            a = 0
            
        ### only one dataset available 
        if len(datasets) == 1:
            best_ds = datasets[0]
            files =  list(self.all_timestamps_dic[ts][best_ds].keys() )
            if len(files) ==1:  ### must only loop over all files   => only one dataset, only one file available easiest case 
                best_file = files[0]
                # [num_valid_temp, num_valid_wind], [min_press_temp, min_press_wind, max_press_temp, max_press_wind], [max_gph, max_pressure, max_height]
            for f in files:
                    valid_t_w, min_t_w_max_t_w, min_p_h_gph, max_p_h_gph  = self.extract_record_data(ts, best_ds, f)   
                    extracted_best_record[f] = [ best_ds, valid_t_w, min_t_w_max_t_w, min_p_h_gph, max_p_h_gph] 

        else:  # hierarchical choices, according to datasets and files 
            ### first find best datasets, then check length of files 
            if ts < 1577836800 and  'giub' in datasets :
                    best_ds = 'giub'
                        
            else:
                if len(era5_ds) >0 and 'igra2' not in datasets:
                    best_ds = era5_ds[0]
                elif len(era5_ds) ==0 and 'igra2' in datasets:
                    best_ds = 'igra2'
                else:
                    best_ds = False #to be determined  
                
            for d in datasets:
                files =  list(self.all_timestamps_dic[ts][d].keys() )
                for f in files:
                    
                    if f == '/scratch/das/federico/HARVEST_YEARLY_10NOV2023_Vienna_Lin_Bethel//ncar/0-20001-0-11035/0-20001-0-11035_1979_ncar_harvested_uadb_trhc_11035.txt.nc':
                        a = 0

                    valid_t_w, min_t_w_max_t_w, min_p_h_gph, max_p_h_gph  = self.extract_record_data(ts, d, f)                     
                    extracted_best_record[f] = [d, valid_t_w, min_t_w_max_t_w, min_p_h_gph, max_p_h_gph] 

        ### HERE: analyze extracted_best_record dictionary with data 
                        
        if not best_file:
            if len(extracted_best_record.keys()) ==1:
                best_file = list(extracted_best_record.keys() )[0]
                best_ds = extracted_best_record[best_file][0]
                valid_t_w, min_t_w_max_t_w, min_p_h_gph, max_p_h_gph  = self.extract_record_data(ts, best_ds, best_file)                       
                extracted_best_record[best_file] = [best_ds, valid_t_w, min_t_w_max_t_w, min_p_h_gph, max_p_h_gph] 

                
            else:  # loop and determine best record 
                
                max_height = -999999
                min_height = 999999
                min_pressure = 999999
                max_pressure = -999999
                best_file = ''
                best_ds = ''
                
                for file in list(extracted_best_record.keys() ) :                    
                    # [num_valid_temp, num_valid_wind], [min_press_temp, min_press_wind, max_press_temp, max_press_wind], [min_pressure, min_gph, min_height], [max_pressure, max_gph, max_height]
                    # where [dataset, valid_t_w, min_t_w_max_t_w, min_p_h_gph, max_p_h_gph]
                    #current_sum = extracted_best_record[file][1][0] + extracted_best_record[file][1][1]  # [d, valid_t_w, min_t_w_max_t_w, min_p_h_gph, max_p_h_gph] 
                    current_min_height = extracted_best_record[file][3][2]                    
                    current_max_height = extracted_best_record[file][4][2]
                    
                    current_min_pressure = extracted_best_record[file][3][0]
                    current_max_pressure = extracted_best_record[file][4][0]
                    
                    current_ds =  extracted_best_record[file][0]
                    current_file = file 
                    
                    
                    more_data = bool(current_min_pressure < min_pressure or current_max_pressure> max_pressure or  current_max_height > max_height or current_min_height < min_height) # boolean variable for more available data 
                    
                    if  more_data :  # current pressure is lower than previous OR max pressure higher than previous
                        
                        max_dist_level = max ( abs(current_min_pressure-min_pressure) , abs(current_max_pressure-max_pressure) , abs(current_max_height-max_height) , abs(current_min_height- min_height) )
    
                        if max_dist_level < 10: # very little difference (pascal or meters) 
                            
                            if best_ds in era5_ds:  # preference to era5 
                                continue
                            
                            else:
                                if current_ds in era5_ds or current_ds == 'igra2':  # preference to era5 or igra2 
                                    best_file = current_file
                                    best_ds = current_ds
                                    min_pressure = current_min_pressure
                                    max_pressure = current_max_pressure        
                                    min_height = current_min_height
                                    max_height = current_max_height
                                    
                                else:
                                    if current_ds == 'ncar' and best_ds:  # discard ncar as much as possible 
                                        continue
                                    
                                    else:
                                        best_file = current_file
                                        best_ds = current_ds
                                        min_pressure = current_min_pressure
                                        max_pressure = current_max_pressure                                           
                                        min_height = current_min_height
                                        max_height = current_max_height                                            
                                                
                        elif max_dist_level >= 10: # significantly more data
                            best_file = current_file
                            best_ds = current_ds
                            min_pressure = current_min_pressure
                            max_pressure = current_max_pressure   
                            min_height = current_min_height
                            max_height = current_max_height                            
                            
                    else:
                        if current_ds in era5_ds:  # preference to era5  
                            best_file = current_file
                            best_ds = current_ds
                            min_pressure = current_min_pressure
                            max_pressure = current_max_pressure    
                            min_height = current_min_height
                            max_height = current_max_height
                            
                        elif current_ds == 'igra2' and best_ds in era5_ds: # preference to era5
                            continue
                        
                        elif current_ds == 'igra2' and best_ds not in era5_ds:  # preference to igra if era5 not available
                            best_file = current_file
                            best_ds = current_ds
                            min_pressure = current_min_pressure
                            max_pressure = current_max_pressure         
                            min_height = current_min_height
                            max_height = current_max_height

                        else:
                            if current_ds == 'ncar' and best_ds and best_ds != 'ncar':  # discard ncar as much as possible 
                                        continue
                            else:
                                best_file = current_file
                                best_ds = current_ds
                                min_pressure = current_min_pressure
                                max_pressure = current_max_pressure   
                                min_height = current_min_height
                                max_height = current_max_height                                
                     
                        
                        
                    #if  current_max_height >= max_height and max_height >0 and  not (current_min_pressure <= min_pressure): # current pressure higher than previous (pressure and height should behave oppositely)
                    #    a = 0
                    #    print('CHECK WRONG')
                    #    # must be something wrong with pressure and height ??? 

        if not best_file:
            best_file = list(extracted_best_record.keys() )[0]
            best_ds = extracted_best_record[best_file][0]
            
        #if best_ds == 'bufr':
        #    a = 0 
            
        self.all_timestamps_dic[ts]['extract_record_data'] = extracted_best_record[best_file]
        self.all_timestamps_dic[ts]['extract_record_data_all'] = extracted_best_record
        self.all_timestamps_dic[ts]['best_record'] = [ best_ds , best_file ]
        
        #dummy = self.plot_profile(ts)
        
        return best_ds, best_file, data_policy 
        
        
    def merge_timestamp_u(self):
        """ Apply merging selection criteria to each reduced timestamp (no more duplicated timestamp)"""
        
        print('=== Merging all time stamps ')
        # unique ts without duplicates, duplicated timestamps
        #unique_ts2, duplicated_ts2 = self.reduce_timestamps()
        unique_ts, duplicated_ts,all_timestamps = self.reduce_timestamps_u()
        # all available ts (unique or possibly duplicated)
        #all_timestamps = list(self.all_timestamps_dic.keys())
        
        #all_timestamps = all_timestamps[:300]  # speed it up TO DO CHANGE !!!!! 
        #all_timestamps.sort()

        # container for processed ts 
        processed_timestamps = []
        all_combined_timestamps ={}        
            
        ## dictionary holding all the data for a given timestamp 
        ## keys = [ 'timestamp' , 'policy' , 'all_duplicated_dt', 'all_duplicated_datasets', 'all_duplicated_files', 'best_ds' , 'best_file', 'ind_inf', 'ind_sup'] 
        ## all_combined_timestamp_data = {}
        ## for k in keys:
        ##    all_combined_timestamp_data[k] = []
            
        ### TO DO
        all_timestamps = [ f for f in all_timestamps ] #if f != self.last_timestamp ] # TODO it has to do with the duplicate check from previous year! see Lindenberg dt=1988060400 (1962-1963)
        all_era5 = ['era5_1', 'era5_1_mobile' , 'era5_2' , 'era5_2_mobile'] 

        i = -1
        duplicates_dic={}
        for dt,index in zip( tqdm(all_timestamps), range(len(all_timestamps)) ) :  #TODO TO DO WRONG CHANGE HERE !!!! 
            
            i += 1
            #if dt == 2492985600:
            #    a = 0
            self.observed_ts_values = {} # restore empty dic 
            
            #if dt == 3339824820 or index == 746:
                #x = 0
            #if dt in processed_timestamps[-5:]:   # all timestamps are sorted so should not check the entire list 
                #continue 

            dtsave = dt

            x = 1
            if len(duplicated_ts[i]) == 1: # no time duplicate detected, apply standard merging procedure 
                real_time = dt 
                all_times = [dt] # used for plotting 
                
                best_ds, best_file, policy = self.find_best_record(dt)
                
                processed_timestamps.append(dt)
                
                sel = {'unique': True}
                duplicates_dic[real_time] = [best_ds, best_file, 0, 0, False]
#if dt == 3339824820:
                    #x = 0
                   
            else:
                if dt == 3689434800:
                    x = 0
                possible_duplicates = duplicated_ts[i] #[ p for p in duplicated_ts if dt in p ]
                
                #if len(possible_duplicates) >0:
##                    possible_duplicates2= [p for p in possible_duplicates[0] if np.abs(p-dtsave) <= 10800]
                    #possible_duplicates2= [p for p in possible_duplicates[0]]
                    #possible_duplicates= possible_duplicates[0] #[p for p in possible_duplicates[0] if np.abs(p-dtsave) <= 7200]
                #else:
                    #print('Check inconsistency with dt , might be due to merging of different years ')
                    #continue
                ## apply hierarchical selection 
                
                #if len(possible_duplicates) >=5: ### (should always be true...)
                    #x = 0
                for t in possible_duplicates:
                    
                    #if np.abs(t-dt) > 10800:    
                        #continue
                    #duplicate_data = {}  # becomign unecessary
                    
                    
                    if len(possible_duplicates) >=2: ### (should always be true...)
                     
                        max_length = 0
                        for dt2 in possible_duplicates:
                            if dt2 not in duplicates_dic.keys():
                                
                                best_ds_low, best_file_low, policy    = self.find_best_record(dt2)   
                                total_length = self.all_timestamps_dic[dt2]['extract_record_data'][1][0] + self.all_timestamps_dic[dt2]['extract_record_data'][1][0]  # wind plus temperature observations
                                if total_length > max_length:
                                    max_length = total_length
                                duplicates_dic[dt2] = [best_ds_low, best_file_low, policy, total_length,False ]  # must extract data also for non selected dt to make plots
                            else:
                                if duplicates_dic[dt2][3] > max_length:
                                    max_length = duplicates_dic[dt2][3]
                                

                            
                        ### skip ncar if multiple available...
                        
                        high_ts = [t for t in possible_duplicates if  duplicates_dic[t][3] == max_length and not duplicates_dic[t][4]]
                        if len(high_ts) > 1:
                            valid = [ np.sum(self.all_timestamps_dic[t]['extract_record_data'][1]) for t in high_ts ]
                            high_ts = high_ts[np.argmax(valid)]
                        else:
                            if len(high_ts) > 0:
                                high_ts = high_ts[0]
                            else:
                                continue
                            
                        
                        all_datasets = np.unique ( [ duplicates_dic[k][0] for k in possible_duplicates ] )
                        
                        if len(all_datasets) >2:  # skipping ncar as much as possible 
                            
                            try:
                                low_ts = min( [i for i in possible_duplicates if duplicates_dic[i][0] != 'ncar'  and i != high_ts and not duplicates_dic[i][4]]  )
                            except:
                                low_ts = high_ts
                            #high_ts = max( [i for i in possible_duplicates if duplicates_dic[i][0] != 'ncar' ])   
                            
                        else:  # worst case, maximum is equal to minimum
                            try:
                                
                                low_ts = min( [i for i in possible_duplicates if i != high_ts and not duplicates_dic[i][4]]  )
                            except:
                                low_ts = high_ts


                            
                        #datasets_low = list(self.all_timestamps_dic[low_ts].keys() ) 
                        best_ds_low, best_file_low, policy    =   duplicates_dic[low_ts][0] , duplicates_dic[low_ts][1]  , duplicates_dic[low_ts][2] 
                        best_ds_high, best_file_high, policy =   duplicates_dic[high_ts][0] , duplicates_dic[high_ts][1]  , duplicates_dic[high_ts][2] 
               
                        ### [ best_ds, valid_t_w, min_t_w_max_t_w, min_p_h_gph, max_p_h_gph ]
                        ### example: 
                        ### self.all_timestamps_dic[low_ts]['extract_record_data'] = ['igra2', [0, 2], [1961.33, 1961.33, 3432.3274, 3432.3274], [999999, 1961.33, 999999], [0, 3432.3274, 0]]
                        
                        #1. check total temp+wind record length
                        l_low =  self.all_timestamps_dic[low_ts]['extract_record_data'][1][0] + self.all_timestamps_dic[low_ts]['extract_record_data'][1][0]
                        l_high =  self.all_timestamps_dic[high_ts]['extract_record_data'][1][0] + self.all_timestamps_dic[high_ts]['extract_record_data'][1][0]
                        
                        min_p_low =  self.all_timestamps_dic[low_ts]['extract_record_data'][3][0] # minimum pressure
                        min_p_high = self.all_timestamps_dic[high_ts]['extract_record_data'][3][0]
                        
                        max_p_low =  self.all_timestamps_dic[low_ts]['extract_record_data'][4][0] # maximum pressure
                        max_p_high = self.all_timestamps_dic[high_ts]['extract_record_data'][4][0]
                        
                        max_h_low =  self.all_timestamps_dic[low_ts]['extract_record_data'][4][1] # maximum height
                        max_h_high = self.all_timestamps_dic[high_ts]['extract_record_data'][4][1]

                        sel = {'unique': False,}
                        if max ( abs(min_p_low - min_p_high) ,  abs(max_p_low - max_p_high) ) < 10:  # small negligible difference, get ERA5 over IGRA over rest
                            
                            if best_ds_low in all_era5 and best_ds_high not in all_era5:
                                high = False
                                
                            elif best_ds_high in all_era5 and best_ds_low not in all_era5:
                                high = True
                                
                            elif best_ds_high == 'igra2':
                                high = True
                            
                            elif best_ds_low == 'igra2':
                                high = False
                                
                            elif best_ds_low in all_era5 and best_ds_high in all_era5: # take longest records within era5
                                
                                if ( l_low >= l_high ):
                                    high = False
                                else:
                                    high = True
                                    
                            elif best_ds_low == best_ds_high:
                                high = True
                            else:
                                
                                if ( l_low >= l_high ):
                                    high = False
                                else:
                                    high = True
                            
                        else:
                            # selecting best ds by number of records
                            if ( l_low >= l_high ):
                                high = False
                            else:
                                high = True
                            # selecting best ds by lowest pressure          
                            if abs(min_p_low-min_p_high) <10:
                                pass
                            else:
                                if min_p_low < min_p_high:
                                    high = False
                                else:
                                    high = True
                            # selecting best ds by highest height   
                            if abs(max_h_low-max_h_high) <10:
                                pass
                            else:
                                if  max_h_low > max_h_high: # current pressure higher than previous (pressure and height should go together)
                                    high = False
                                else:
                                    high = True
                                        
                
                        try:
                            
                            if high:
                                best_ds = best_ds_high
                                best_file = best_file_high
                                real_time = high_ts                                
                                sel['l'] = l_high
                                sel['min_p'] = min_p_high
                                sel['max_h'] = max_h_high
                            else:
                                best_ds = best_ds_low
                                best_file = best_file_low
                                real_time = low_ts
                                sel['l'] = l_low
                                sel['min_p'] = min_p_low
                                sel['max_h'] = max_h_low
                        except Exception as e:
                            raise ValueError(self.dic_h5py)
                            
                if 'real_time' not in locals():
                    continue
                duplicated_time = [t for t in possible_duplicates if t != real_time ]
                for d in duplicated_time:
                    duplicates_dic[d][4] = True
                all_times = possible_duplicates
                    
                    
            processed_timestamps.extend(all_times)

            ### Producing check plots
            try:
                
                if int(self.current_year) > 2004 and  int(self.current_year) < 2011 and  index % 50 == 0 :
                    dummy = self.plot_profile_extended(all_times, real_time = real_time, best_file= best_file )
                    
                if int(self.current_year) < 1950 and  index % 50 == 0:
                    dummy = self.plot_profile_extended(all_times, real_time = real_time, best_file= best_file )
                else:
                    if index % 500 == 0: #or x == 0:
                        #dummy = self.plot_profile(all_times, real_time = real_time, best_file= best_file )
                        dummy = 0 #self.plot_profile_extended(all_times, real_time = real_time, best_file= best_file )
            except:
                
                print('could not plot')
                
            ### Saving the extracted data
            replace = True
            if real_time in  all_combined_timestamps.keys() :
                #if all_combined_timestamps[real_time]['best_ds'] != best_ds:                   
                #print(f'DUPLICATE {real_time} !')
                replace = False
                if sel['min_p'] < all_combined_timestamps[real_time]['sel']['min_p']:
                    replace = True
                elif  sel['min_p'] == all_combined_timestamps[real_time]['sel']['min_p']:                   
                    if sel['l'] > all_combined_timestamps[real_time]['sel']['l']:
                        replace = True
                if sel['max_h'] > all_combined_timestamps[real_time]['sel']['max_h']:
                    replace = True
                elif  sel['max_h'] == all_combined_timestamps[real_time]['sel']['max_h']:                   
                    if sel['l'] > all_combined_timestamps[real_time]['sel']['l']:
                        replace = True
                        
            if replace:
                all_combined_timestamps[real_time]= {}
                all_combined_timestamps[real_time]['sel'] = copy.deepcopy(sel)
                if not duplicates_dic[real_time][4]:
                    all_combined_timestamps[real_time]['policy'] = policy
        
                    all_combined_timestamps[real_time]['best_ds'] = best_ds
                    all_combined_timestamps[real_time]['best_file'] = best_file
                    #all_combined_timestamps[real_time]['sel'] = copy.deepcopy(sel)
            #print(real_time, all_combined_timestamps[real_time]['sel'])
            
            # duplicated sources
            if len(duplicates_dic.keys()) >1 :
                all_ds = np.unique([duplicates_dic[k][0] for k in duplicates_dic.keys() ])            
                duplicated = ','.join(all_ds) 
                all_combined_timestamps[real_time]['all_duplicated_records'] = duplicated
                status=0
            else:
                duplicated = best_ds
                status =1
                all_combined_timestamps[real_time]['all_duplicated_records'] = duplicated

            # duplicated status
            all_combined_timestamps[real_time]['duplicated_status'] = status
            
            all_combined_timestamps[real_time]['real_time'] = real_time
            all_combined_timestamps[real_time]['duplicated_time'] = all_times
            
        if 'dt' in locals():        
            #self.last_timestamp = dt # last timestamp updated in the loop
            self.merged_timestamp = {}
            for k, v in  all_combined_timestamps.items():
                if not duplicates_dic[k][4]:
                    self.merged_timestamp[k] = v
                    self.last_timestamp = k # last timestamp updated in the loop
            #self.merged_timestamp = [all_combined_timestamps[t] for t in  all_combined_timestamps.keys() if not duplicates_dic[t][4]] 
        else:
            print('ERROR, no dt exists')
            raise ValueError(self.dic_h5py)
        
        print("DONE Merging all timestamps ")
        
        
    def merge_timestamp(self):
        print("Don't use me")
        return
        """ Apply merging selection criteria to each reduced timestamp (no more duplicated timestamp)"""
        
        print('=== Merging all time stamps ')
        # unique ts without duplicates, duplicated timestamps
        unique_ts, duplicated_ts = self.reduce_timestamps()
        # all available ts (unique or possibly duplicated)
        all_timestamps = list(self.all_timestamps_dic.keys())
        
        #all_timestamps = all_timestamps[:300]  # speed it up TO DO CHANGE !!!!! 
        all_timestamps.sort()

        # container for processed ts 
        processed_timestamps = []
        all_combined_timestamps ={}        
            
        ## dictionary holding all the data for a given timestamp 
        ## keys = [ 'timestamp' , 'policy' , 'all_duplicated_dt', 'all_duplicated_datasets', 'all_duplicated_files', 'best_ds' , 'best_file', 'ind_inf', 'ind_sup'] 
        ## all_combined_timestamp_data = {}
        ## for k in keys:
        ##    all_combined_timestamp_data[k] = []
            
        ### TO DO
        all_timestamps = [ f for f in all_timestamps if f != self.last_timestamp ] # TODO it has to do with the duplicate check from previous year! see Lindenberg dt=1988060400 (1962-1963)
        all_era5 = ['era5_1', 'era5_1_mobile' , 'era5_2' , 'era5_2_mobile'] 

        if len(all_timestamps) > 25000:
            print('too many timestamps, skipping..')
            dt = 0
        else:

            for dt,index in zip( tqdm(all_timestamps), range(len(all_timestamps)) ) :  #TODO TO DO WRONG CHANGE HERE !!!! 
                
                #if dt == 2492985600:
                #    a = 0
                self.observed_ts_values = {} # restore empty dic 
                
                if dt in processed_timestamps[-5:]:   # all timestamps are sorted so should not check the entire list 
                    continue 
    
                duplicates_dic={}
    
                if dt in unique_ts: # no time duplicate detected, apply standard merging procedure 
                    real_time = dt 
                    all_times = [dt] # used for plotting 
                    
                    best_ds, best_file, policy = self.find_best_record_u(dt)
                    
                    processed_timestamps.append(dt)
                    #print(dt)
                        
                else:
                    possible_duplicates = [ p for p in duplicated_ts if dt in p ]
                    
                    if len(possible_duplicates) >0:
                        possible_duplicates= possible_duplicates[0]
                        #timedelta = self.hour_time_delta * 60 * 60
                        #possible_duplicates= [p for p in possible_duplicates[0] if np.abs(dt - p) <= timedelta]
                        #if  len(possible_duplicates) >20:
                            #pass
                        #print('too many possible duplicates, limiting to 20')
                        #possible_duplicates= possible_duplicates[:20]
    
                        #continue
                    else:
                        print('Check inconsistency with dt , might be due to merging of different years ')
                        continue
                    # apply hierarchical selection 
                    
                    i = 0
                    for t in possible_duplicates:
                        
                        i += 1
                        #print(i, t, dt)
                        #duplicate_data = {}  # becomign unecessary 
                        
                        if len(possible_duplicates) >=2: ### (should always be true...)
                         
                            max_length = 0
                            for pdt in possible_duplicates:
                                
                                best_ds_low, best_file_low, policy    = self.find_best_record_u(pdt)   
                                total_length = self.all_timestamps_dic[pdt]['extract_record_data'][1][0] + self.all_timestamps_dic[pdt]['extract_record_data'][1][0]  # wind plus temperature observations
                                if total_length > max_length:
                                    max_length = total_length
                                duplicates_dic[pdt] = [best_ds_low, best_file_low, policy, total_length ]  # must extract data also for non selected dt to make plots 
    
                                
                            ### skip ncar if multiple available...
                            
                            high_ts = [t for t in possible_duplicates if  duplicates_dic[t][3] == max_length ][0]
                            
                            all_datasets = np.unique ( [ duplicates_dic[k][0] for k in possible_duplicates ] )
                            
                            if len(all_datasets) >2:  # skipping ncar as much as possible 
                                
                                low_ts = min( [i for i in possible_duplicates if duplicates_dic[i][0] != 'ncar'  and i != high_ts ]  )
                                #high_ts = max( [i for i in possible_duplicates if duplicates_dic[i][0] != 'ncar' ])   
                                
                            else:  # worst case, maximum is equal to minimum 
                                low_ts = min( [i for i in possible_duplicates if i != high_ts ]  )
    
    
                                
                            #datasets_low = list(self.all_timestamps_dic[low_ts].keys() ) 
                            best_ds_low, best_file_low, policy    =   duplicates_dic[low_ts][0] , duplicates_dic[low_ts][1]  , duplicates_dic[low_ts][2] 
                            best_ds_high, best_file_high, policy =   duplicates_dic[high_ts][0] , duplicates_dic[high_ts][1]  , duplicates_dic[high_ts][2] 
                   
                            ### [ best_ds, valid_t_w, min_t_w_max_t_w, min_p_h_gph, max_p_h_gph ]
                            ### example: 
                            ### self.all_timestamps_dic[low_ts]['extract_record_data'] = ['igra2', [0, 2], [1961.33, 1961.33, 3432.3274, 3432.3274], [999999, 1961.33, 999999], [0, 3432.3274, 0]]
                            
                            #1. check total temp+wind record length
                            l_low =  self.all_timestamps_dic[low_ts]['extract_record_data'][1][0] + self.all_timestamps_dic[low_ts]['extract_record_data'][1][0]
                            l_high =  self.all_timestamps_dic[high_ts]['extract_record_data'][1][0] + self.all_timestamps_dic[high_ts]['extract_record_data'][1][0]
                            
                            min_p_low =  self.all_timestamps_dic[low_ts]['extract_record_data'][3][0] # minimum pressure
                            min_p_high = self.all_timestamps_dic[high_ts]['extract_record_data'][3][0]
                            
                            max_p_low =  self.all_timestamps_dic[low_ts]['extract_record_data'][4][0] # maximum pressure
                            max_p_high = self.all_timestamps_dic[high_ts]['extract_record_data'][4][0]
                            
                            max_h_low =  self.all_timestamps_dic[low_ts]['extract_record_data'][4][1] # maximum height
                            max_h_high = self.all_timestamps_dic[high_ts]['extract_record_data'][4][1]
    
                            if max ( abs(min_p_low - min_p_high) ,  abs(max_p_low - max_p_high) ) < 10:  # small negligible difference, get ERA5 over IGRA over rest
                                
                                if best_ds_low in all_era5 and best_ds_high not in all_era5:
                                    best_ds = best_ds_low
                                    best_file = best_file_low
                                    real_time = low_ts
                                    
                                elif best_ds_high in all_era5 and best_ds_low not in all_era5:
                                    best_ds = best_ds_high
                                    best_file = best_file_high                                  
                                    real_time = high_ts
                                    
                                elif best_ds_high == 'igra2':
                                    best_ds = best_ds_high
                                    best_file = best_file_high                                  
                                    real_time = high_ts            
                                
                                elif best_ds_low == 'igra2':
                                    best_ds = best_ds_low
                                    best_file = best_file_low
                                    real_time = low_ts      
                                    
                                elif best_ds_low in all_era5 and best_ds_high in all_era5: # take longest records within era5
                                    
                                    if ( l_low >= l_high ):
                                        best_ds = best_ds_low
                                        best_file = best_file_low
                                        real_time = low_ts
                                    else:
                                        best_ds = best_ds_high
                                        best_file = best_file_high                        
                                        real_time = high_ts           
                                        
                                elif best_ds_low == best_ds_high:
                                    best_ds = best_ds_high
                                    best_file = best_file_high
                                    real_time = high_ts                                
                                
                            else:
                                # selecting best ds by number of records 
                                if ( l_low >= l_high ):
                                    best_ds = best_ds_low
                                    best_file = best_file_low
                                    real_time = low_ts
                                else:
                                    best_ds = best_ds_high
                                    best_file = best_file_high                        
                                    real_time = high_ts
                                # selecting best ds by lowest pressure          
                                if abs(min_p_low-min_p_high) <10:
                                    continue
                                else:
                                    if min_p_low < min_p_high:
                                        best_ds = best_ds_low
                                        best_file = best_file_low   
                                        real_time = low_ts
                                    else:
                                        best_ds = best_ds_high
                                        best_file = best_file_high                                  
                                        real_time = high_ts
                                # selecting best ds by highest height   
                                if abs(max_h_low-max_h_high) <10:
                                    continue
                                else:
                                    if  max_h_low > max_h_high: # current pressure higher than previous (pressure and height should go together)
                                        best_ds = best_ds_low
                                        best_file = best_file_low    
                                        real_time = low_ts
                                    else:
                                        best_ds = best_ds_high
                                        best_file = best_file_high                           
                                        real_time = high_ts
                                            
                    if real_time not in locals():
                        real_time = possible_duplicates[0] +1   # dirty fix
                    duplicated_time = [t for t in possible_duplicates if t != real_time ][0]
                    all_times = possible_duplicates
                        
                        
                processed_timestamps.extend(all_times)
    
                ### Producing check plots 
                if int(self.current_year) > 2004 and  int(self.current_year) < 2011 and  index % 50 == 0 :
                    dummy = self.plot_profile_extended(all_times, real_time = real_time, best_file= best_file )
                    
                if int(self.current_year) < 1950 and  index % 50 == 0:
                    dummy = self.plot_profile_extended(all_times, real_time = real_time, best_file= best_file )
                else:
                    if index % 500 == 0:
                        #dummy = self.plot_profile(all_times, real_time = real_time, best_file= best_file )
                        dummy = self.plot_profile_extended(all_times, real_time = real_time, best_file= best_file )
    
                ### Saving the extracted data             
                all_combined_timestamps[real_time]= {}            
                all_combined_timestamps[real_time]['policy'] = policy
    
                all_combined_timestamps[real_time]['best_ds'] = best_ds
                all_combined_timestamps[real_time]['best_file'] = best_file
                
                # duplicated sources
                if len(duplicates_dic.keys()) >1 :
                    all_ds = np.unique([duplicates_dic[k][0] for k in duplicates_dic.keys() ])            
                    duplicated = ','.join(all_ds) 
                    all_combined_timestamps[real_time]['all_duplicated_records'] = duplicated
                    status=0
                else:
                    duplicated = best_ds
                    status =1
                    all_combined_timestamps[real_time]['all_duplicated_records'] = duplicated
    
                # duplicated status
                all_combined_timestamps[real_time]['duplicated_status'] = status
                
                all_combined_timestamps[real_time]['real_time'] = real_time
                all_combined_timestamps[real_time]['duplicated_time'] = all_times
        
        if 'dt' in locals():        
            self.last_timestamp = dt # last timestamp updated in the loop 
            self.merged_timestamp = all_combined_timestamps
        else:
            print('ERROR, no dt exists')
        
        print("DONE Merging all timestamps ")
        
        
    def get_null(self, tipo = ''):
        ''' Simply returns the proper format for ''null' value '''        
        if tipo == np.int32 or tipo == np.int64:
            void = -2147483648
        elif tipo == np.float32 or tipo == np.float64 :
            void = np.nan
        elif tipo == np.bytes_ :
            void = b'NA '
        else:
            return -2147483648
        return void
    

    def merge_all_data(self):       
        """ Construct a dictionary with all the dataframes of each dataset, either reading it from saved pickles or reading it from memory """

        logging.info('***** Starting the merging process merge_all_data')

        # avoidable loops, should not matter much  # self.all_years:  
        #for this_year in [1981, 1982, 1983]:  # loop over all available years in the date_times 
        #for this_year in self.all_years:   # loop over all available years in the date_times 
        
        ### here: extract all available years
        
        #for this_year in self.all_years:   # loop over all available years in the date_times
        data_all_years = self.extract_file_per_year()
        all_years = list(data_all_years.keys() )
        
        ### filtering years 
        if year_range:
            all_years = [y for y in all_years if int(y) >= year_range[0] and int(y) <= year_range[1] ]  ### TO DO TODO HERE
        
        if self.station in self.processed_stats.keys():
            print("Filtering out already processed year::: " )
            all_years = [y for y in all_years if y not in self.processed_stats[self.station] ]
            if 'completed' in self.processed_stats[self.station]:
                all_years = []
        
        if len(all_years) ==0:
            print('+++ WARNING: no file checked the selection criteria, check if correct! ')
            return
        else:
            
            for this_year in all_years:   # loop over all available years in the date_times 
                print('=== Running year ::: ' , this_year )
                self.current_year = this_year 
                
                data_this_year = data_all_years[this_year]
                
                #following_year = str(int(this_year) + 1) 
                #if following_year in all_years:
                #    data_following_year = data_all_years[following_year]
                #else:
                #    data_following_year = ''
                
                dummy = self.open_data(data_this_year)
                """ Making all date_time for the selected year  """
    
                dummy = self.make_unique_datetime()            
                
                # principal merging algoritm, applying selection rules
                    
                merged = self.merge_timestamp_u() # selection of the best data
                
                dummy = self.initialize_out_file()
                
                if len(self.merged_timestamp) == 0:
                    continue
                    
                print('=== Making and writing standard fixed CDM tables')                                    
                dummy = self.make_standard_cdm_table()

                print('=== Making and writing observations_table')
                dummy = self.make_merged_observation_table()
                
                print('=== Making and writing recordindex')
                dummy = self.make_merged_recordindex()
                
                print('=== Making and writing feedback_table')                        
                dummy = self.make_merged_feedback_table()  ### TODO TO DO HERE 
                
                print('=== Making and writing source_conf table')                                    
                dummy = self.make_sourceconf_table()
                
                print('=== Making and writing header_table')                                    
                dummy = self.make_merged_header_table()
                
                
                a = 'HERE'  # TO DO HERE
            
                ### Adding sensor_id to observations_table 
                if self.add_sensor:
                    print('*** Adding sensor *** ')        
                    add_sensor = wrapper(out_dir = self.out_dir , station_id = self.station.split('-')[-1] , file = self.out_file , copy = self.copy )
                    print('*** DONE Adding sensor *** ')
                else:
                    pass
                
                ### Writing log of completed year
                a = open(self.out_dir+'/'+ self.station  + '/completed.txt' , 'a+')
                if str(this_year)+'\n' not in a.readlines():
                    a.write(this_year + '\n')
                    
        try:
            a = open(self.out_dir+'/'+ self.station  + '/completed.txt' , 'a+')
            lines = a.readlines()
            if 'completed\n' not in lines:        
                a.write('completed\n')
        except:
            pass

    def initialize_out_file(self):
        """ Create output directory, initialize out file """
        
        out_dir = self.out_dir  + '/' + self.station 
        out_file = out_dir + '/'  + self.station + '_' + self.current_year + '_CEUAS_merged_v3.nc'  
        
        if not os.path.isdir(out_dir):
            Path(out_dir).mkdir(parents=True, exist_ok=True)  
            
        self.out_file = out_file 
        
        ### TO DO REMOVE
        if os.path.isfile(out_file):
            os.remove(out_file)        
            
    def make_standard_cdm_table(self, out_file='' ):
        """ Create standard CDM table """
            
        self.write_merged_new(var='', table = 'cdm_tables', data='')
        
    def make_sourceconf_table(self, out_file='' ):
        """ Create source_configuration CDM table """
        
        ### HEADER TABLE
        res = self.merged_timestamp
        all_ts = list(res.keys() )
        files = []
        for ts in all_ts :
            best_file = res[ts]['best_file']
            if len(best_file) < 200:
                best_file = best_file + ' '*(200-len(best_file))
            files.append(np.bytes_(best_file))
        if len(files) == 0:
            files = [b'NA ']
        
        self.write_merged_new(var='source_file', table = 'source_configuration', data=np.array(files))
        
        
    def make_merged_feedback_table(self, out_file='' ):
        """ Create merged era5fb table """
        ### FEEDBACK TABLE
        tt = time.time()
        res = self.merged_timestamp
        
        # all date_time for this year 
        all_ts = list(res.keys() )
        
        # variables to read from the era5fb. Note that the effective variables contained in each era5fb table is not fixed, it depends on the exact year
        variables = self.dic_type_attributes['era5fb'].keys()
                
        # here, only extract files that are actually selected as best_file (will not pre-load the others) ( I dont remember why this, check! )
        files = list( np.unique( [ res[dt]['best_file'] for dt in all_ts ] ) ) 
        all_files_to_preload = list( [ f for f in files if f.split('/')[-3] in ['era5_1',  'era5_1_mobile' , 'era5_2' , 'era5_2_mobile'] ]  )
        
        #era5_vars =  list( h5py.File(all_files[0], 'r')['era5fb'].keys() )
        if 'datum_sfc_event@surfbody_feedback' not in variables:
            variables.append('datum_sfc_event@surfbody_feedback')
        if 'source_id' not in variables:
            variables.append('source_id')
        
        fps =[]    
        for f in all_files_to_preload:
            fps.append(h5py.File(f, 'r'))
        for v in variables:
            
            #if v in ['datum_anflag@body','datum_event1@body','datum_rdbflag@body','report_event1@hdr','report_rdbflag@hdr', 'varbc_ix@body' ]:
            #    a=0
            
            if v == 'source_id':  # why?
                continue
            
            '''
            if v in ['collection_identifier@conv', 'expver', 'source@hdr', 'source_id' , 'statid@hdr']:
                continue
            if v in ['expver' , 'source@hdr', ]:
                continue
            if v in [ 'expver', 'source@hdr', 'source_id' , 'statid@hdr']:
                    a=0                
            '''
            # keep full data in memory, for better efficiency
            #print("        Variable " , v )
            load_full_data = {}
            for fp in fps:
                fpe = fp['era5fb']
                file_era5fb_variables = list(fpe.keys() ) # this list changes with the year 
            
                if v in file_era5fb_variables:
                    data_v = fpe[v][:]
                else:
                    n = self.get_null( self.encodings['era5fb'][v]['type']  )
                    data_v = np.full( len( fp['observations_table']['date_time']), n)                          

                load_full_data[fp.filename] = data_v     
            
            data = []
        
            for ts, size in zip(all_ts, self.record_size):
                
                best_file = res[ts]['best_file']
                best_ds = res[ts]['best_ds']
                
                ind_min = self.all_timestamps_dic[ts][best_ds][best_file][0]
                ind_max = self.all_timestamps_dic[ts][best_ds][best_file][1]
            
                if best_ds in ['era5_1' , 'era5_1_mobile' , 'era5_2' , 'era5_2_mobile']: #only for these datasets the feedback exists
                    sliced_data = load_full_data[best_file][ind_min:ind_max]
                    if v in ['collection_identifier@conv', 'expver', 'source@hdr' , 'statid@hdr', 'source_id']:
                        if type(sliced_data[0]) in [np.int32 , np.int64, np.float64]:
                            pass
                        else:
                            #sliced_data = np.array( [ b''. join(d) if isinstance(d, np.ndarray) else np.bytes_(d) for d in sliced_data ] )
                            if sliced_data.ndim == 2:
                                
                                sliced_data = sliced_data.view(f'S{sliced_data.shape[1]}' ).flatten()
                            else:
                                sliced_data = sliced_data.view(f'S{sliced_data.shape[0]}' )
                                
                        
                else: # must write a dummz feedback
                    a = self.get_null(tipo= self.dic_type_attributes['era5fb'][v]['type'])
                    try:
                        
                        sliced_data = np.full( (ind_max-ind_min ) , a )
                    except:
                        #return
                        raise ValueError(f'{files[0]},{ind_max},{ind_min}' )
                    #print('************' , v , '  ' , self.dic_type_attributes['era5fb'][v]['type'] , '  ' ,   sliced_data[:2] , '  ' , a )
                    
                if ts in self.no_duplicated_indices.keys():
                    sliced_data = np.array([sliced_data[i] for i in self.no_duplicated_indices[ts] ])
                        
                #if len(sliced_data) != size:
                #    print(a)
                    
                data.append(sliced_data)
                #print(size, len(sliced_data))
                
            d = np.concatenate(data)
            
            '''
            if v in ['collection_identifier@conv', 'expver', 'source@hdr' , 'statid@hdr']:
                if type(d[0]) in [np.int32 , np.int64]:
                    continue
                
                d = np.array( [ b''. join(d) if isinstance(d, np.ndarray) else np.bytes_(d) for d in data ] )
            '''
            #print(v,d.dtype,d.shape, time.time()-tt)
            dummy_write = self.write_merged_new( var=v, table = 'era5fb', data=d)
                                                 
        for v in ['datum_anflag@body','datum_event1@body','datum_rdbflag@body','report_event1@hdr','report_rdbflag@hdr', 'varbc_ix@body' ]:
            if v not in variables:
                data_v = np.full( len( h5py.File(f, 'r')['observations_table']['date_time']), np.nan )  
                dummy_write = self.write_merged_new( var=v, table = 'era5fb', data=data_v)
        print(time.time()-tt)


    def make_merged_header_table(self, out_file='' ):
        """ Create merged header_table """ 
        #a = self.header_table_rep_id
      
        ### HEADER TABLE
        res = self.merged_timestamp
        
        # all date_time for this year 
        all_ts = list(res.keys() )
        
        # variables to read from the observations_table of each file 
        other_variables = ['report_id' , 'station_record_number' , 'source_id' , 'duplicate_status' , 'duplicates' ,]
        
        variables = [ v for v in self.dic_type_attributes['header_table'].keys() if v not in other_variables  and v != 'primary_station_id' ]
                
        # here, only extract files that are actually selected as best_file (will not pre-load the others)
        all_files = list( np.unique( [ res[dt]['best_file'] for dt in all_ts ] ) )

        lats, lons = [],[]
        
        fps = []
        for f in all_files:
            fps.append(h5py.File(f, 'r'))
            
        for v in variables:
            #if v == 'report_meaning_of_timestamp':
            #    a=0
            # keep full data in memory, for better efficiency 
            load_full_data = {}
            for fp in fps:
                fpn = fp.filename
                load_full_data[fpn] = {}
                if v in fp['header_table'].keys():
                    load_full_data[fpn][v] = fp['header_table'][v][:]
                else:
                    n = self.get_null( self.encodings['header_table'][v]['type']  )
                    load_full_data[fpn][v] = np.full( len( fp['header_table']['report_timestamp']), n)      
            
            data = []
            source_files = []
            source_ids, duplicates, duplicate_status = [],[],[]
            
            for ts in all_ts :
                
                best_file = res[ts]['best_file']
                source_files.append(best_file)

                best_ds = res[ts]['best_ds']
                source_ids.append(best_ds)

                # duplicates 
                duplicates.append(res[ts]['all_duplicated_records'])
                duplicate_status.append(res[ts]['duplicated_status'])
                
                
                #ind_min = self.all_timestamps_dic[ts][best_ds][best_file][0]
                #ind_max = self.all_timestamps_dic[ts][best_ds][best_file][1] # not useful in this case 
                index =  self.all_timestamps_dic[ts][best_ds][best_file][2]
                value = load_full_data[best_file][v][index]
                    
                data.append(value)
                    
            d = np.array(data)
            dummy_write = self.write_merged_new( var=v, table = 'header_table', data=d)
              
            if v == 'latitude':
                lats = d
            if v == 'longitude':
                lons = d 
                
        # station_id
        stat_ids = np.full( (len(data)), np.bytes_(self.station) )
        dummy_write = self.write_merged_new( var='primary_station_id', table = 'header_table', data= stat_ids )
        
        # report_id
        dummy_write = self.write_merged_new( var='report_id', table = 'header_table', data= self.header_table_rep_id )
        
        #source_id, duplicates, duplicate status
        
        dummy_write = self.write_merged_new( var='source_id', table = 'header_table',  data= np.array(source_ids).astype('|S10') ) 
        dummy_write = self.write_merged_new( var='duplicates', table = 'header_table', data= np.array(source_ids).astype('|S30') )
        dummy_write = self.write_merged_new( var='duplicate_status', table = 'header_table', data= duplicate_status )
        
        source_ids, duplicates, duplicate_status = [],[],[]

        
        #self source_configuration_files = source_files
        #### source_configuration file 
        #dummy_write = self.write_merged_new( var='source_files', table = 'source_configuration', data=np.array(source_files ) )
        
        ### make station_configuration
        ### will replicate what I read from CUON and replace the file lat and lon (not the station lat and lon)
        sc = self.stat_conf
        
        '''
        variables = ['primary_id' , 'secondary_id' , 'city' , 'station_type' ,'primary_id_scheme' , 'secondary_id_scheme' ,
            'station_crs' , 'station_type' , 'platform_type' , 'platform_sub_type' , 'operating_institute' , 
            'operating_territory' , 'observed_variables', 'metadata_contact' , 'metadata_contact_role' ]
        for v in self.stat_conf.columns:
            if v in variables:
                sc[v] = np.bytes_( np.array(self.stat_conf[v].values[0]).astype("|S20") )  # TO DO TODO HERE TO DO CHANGE
        '''
        
        sc_all = [ sc for i in range(len(lons)) ] 
        sc_all_df = pd.concat(sc_all)
        sc_all_df['latitude'] = lats
        sc_all_df['longitude'] = lons
        
        variables_str = ['primary_id' , 'secondary_id' , 'city' , 'operating_institute' , 'operating_territory',
            'observed_variables', 'metadata_contact' , 'metadata_contact_role', 'station_name' , ]
        
        variables_int = [ 'station_type' ,'primary_id_scheme' , 'secondary_id_scheme',
            'station_crs' , 'station_type' , 'platform_type' , 'platform_sub_type'  ,  ]      
        
        for v in sc_all_df.columns:
            if "Unnamed" in v: continue 
            
            data = np.array(sc_all_df[v].values )
            
            if v in variables_str:
                if sc_all_df[v].values.dtype in [np.float64, np.float32, np.dtype('O')] and type(sc_all_df[v].values[0] is not str):
                    data = [ str(val) for val in sc_all_df[v].values ]
                else:    
                    data = [ val.encode('utf8') for val in sc_all_df[v].values ]
                try:
                    data = np.array(data).astype('|S30')
                except:
                    try:
                        
                        data = [ val.encode('utf8') for val in sc_all_df[v].values ]
                        data = np.array(data).astype('|S30')
                    except Exception as e:
                        raise ValueError(self.station)
                    
            elif v in variables_int:
                try:
                    data = []
                    for val in sc_all_df[v].values:
                        if val == val:
                            
                            data.append(int(val))
                        else:
                            data.append(-2147483648)
                except:
                    data = [-2147483648 for val in sc_all_df[v].values ]
                data = np.array(data)
                a=0
                
            else:
                data = np.array(sc_all_df[v].values )
                
            dummy_write = self.write_merged_new( var=v, table = 'station_configuration', data=data)
        
        a = 0
        
              
    def make_merged_recordindex(self):
        """ Create merged recordindex """
        
        data = self.date_time
        di=xr.Dataset()
        datetimes, recordindex = np.unique( data, return_index=True )
        di['recordtimestamp'] = ( {'recordtimestamp' : datetimes.shape } , datetimes )
        di['recordtimestamp'].attrs['units']='seconds since 1900-01-01 00:00:00'
        di['recordindex']          = ( {'recordindex' : recordindex.shape } ,  recordindex )
        self.write_merged_new(var='record_index', table='record_index', data=di)   
        return
        
    def  make_merged_observation_table(self, out_file='' ):
        """ Create merged observations,header tables """
        
        tt = time.time()
        res = self.merged_timestamp
        
        # all date_time for this year 
        all_ts = list(res.keys() )
        
        # variables to read from the observations_table of each file 
        all_variables = []
        
        ### list of variables in hte header/observations table that are taken from the harvested files (and not empty variables)
        variables = ['date_time' , 'z_coordinate_type', 'z_coordinate' , 'observed_variable', 'observation_value' ,'longitude', 'latitude' , 'original_units'  ]
        all_variables.extend(variables)
        
        # here, only extract files that are actually selected as best_file (will not pre-load the others)
        all_files = list( np.unique( [ res[dt]['best_file'] for dt in all_ts ] ) ) 
        
        # pre-load once for all
        temp_data = { 'observed_variable':'', 
                    'z_coordinate_type':'' }
        
        fps = []
        for f in all_files:
            fps.append(h5py.File(f, 'r'))

        for fp in fps:
            temp_data[fp.filename] = { 'observed_variable':'', 
                    'z_coordinate_type':'' }
            #with  h5py.File(f, 'r') as fp:
                
            temp_data[fp.filename]['observed_variable']= fp['observations_table']['observed_variable'][:]
            temp_data[fp.filename]['z_coordinate']= fp['observations_table']['z_coordinate'][:]
            
        # pre-loop to check if there are duplicated pressure values within one record
        
        ### for each timestamp, keep the subset of indices to be kept when building the observation and feedback table
        no_duplicated_indices = {}
        
        
        for ts in all_ts :  ### TO DO might reuce the double loops over ths
            
            #if ts == 3760556520:
            #    a = 0 
                
            best_file = res[ts]['best_file']
            best_ds = res[ts]['best_ds']
            
            ind_min = self.all_timestamps_dic[ts][best_ds][best_file][0]
            ind_max = self.all_timestamps_dic[ts][best_ds][best_file][1]
        
            #obs_val = h5py.File(best_file, 'r')['observations_table']['observation_value'][ind_min:ind_max]
            obs_var = temp_data[best_file]['observed_variable'][ind_min:ind_max]
            z_coord = temp_data[best_file]['z_coordinate'][ind_min:ind_max]
            
            ### check if there are internal duplicated entries i.e. duplicated observations 
            unique_var = np.unique(obs_var)
            is_duplicate = False
            
            for v in unique_var: # here I check if, for at least one variable, I find some duplicates i.e. the length of the distinct pressure values is > length of the data for that variable 
                obs_var_ind = np.where( obs_var == v )[0]
                if len( np.unique(z_coord[obs_var_ind])) == len(obs_var_ind):
                    pass
                else:
                    is_duplicate = True
                    print(np.where(z_coord[obs_var_ind][1:]-z_coord[obs_var_ind][:-1]==0))
                    
            if is_duplicate:  # If I find at least one, then I check with pandas and remove the duplicates 
                dic = {}
                ### in this case, there are internal dupliacted entries.
                ### only int his case I create pandas df to easily check duplicated,
                ### then extract the indices to be kept when building the merged observations and feedback tables
                with  h5py.File(best_file, 'r') as fp:
                    
                    for v in ['observed_variable' , 'date_time' , 'z_coordinate_type' , 'z_coordinate']:
                        dic[v] = fp['observations_table'][v][ind_min:ind_max]
                df = pd.DataFrame.from_dict(dic)
                no_dupl = df.drop_duplicates()
                no_duplicated_indices[ts] = list(no_dupl.index)
                #b = 0 ### here must check how to remove duplicates 

        self.no_duplicated_indices = no_duplicated_indices 
        
        
        sizes = []
        
        for v in variables:
            #print("      Variable " , v )  
            # keep full data in memory, for better efficiency 
            load_full_data = {}
            #for f in all_files:
                #load_full_data[f] = h5py.File(f, 'r')['observations_table'][v][:]
            for f in fps:
                load_full_data[f.filename] = f['observations_table'][v][:]
            
            data = []
        
            for ts in all_ts :
                #if ts == 1893283200:
                #    a = 0 
                    
                best_file = res[ts]['best_file']
                best_ds = res[ts]['best_ds']
                
                ind_min = self.all_timestamps_dic[ts][best_ds][best_file][0]
                ind_max = self.all_timestamps_dic[ts][best_ds][best_file][1]
            
                sliced_data = load_full_data[best_file][ind_min:ind_max]
                if ts in no_duplicated_indices.keys():  # here: check if timestamp has internal pressure duplicate to remove 
                    #sliced_data = sliced_data[ no_duplicated_indices[ts] ] # removing duplicated entries 
                    sliced_data = np.array([sliced_data[i] for i in no_duplicated_indices[ts] ])
                    
                data.append(sliced_data)
                
                ### building a vector with the length of the records 
                if v == 'date_time':
                    sizes.append( len(sliced_data) )
                    
                    #if sliced_data[0] == 2115331200:
                    #    a = 0
            d = np.concatenate(data)
            #data = []
        
            #for ts in all_ts :
                ##if ts == 1893283200:
                ##    a = 0 
                    
                #best_file = res[ts]['best_file']
                #best_ds = res[ts]['best_ds']
                
                #ind_min = self.all_timestamps_dic[ts][best_ds][best_file][0]
                #ind_max = self.all_timestamps_dic[ts][best_ds][best_file][1]
            
                #sliced_data = list(load_full_data[best_file][ind_min:ind_max])
                #if ts in no_duplicated_indices.keys():  # here: check if timestamp has internal pressure duplicate to remove 
                    ##sliced_data = sliced_data[ no_duplicated_indices[ts] ] # removing duplicated entries 
                    #sliced_data = [sliced_data[i] for i in no_duplicated_indices[ts] ] 
                    
                #data.extend(sliced_data)
                
                #### building a vector with the length of the records 
                #if v == 'date_time':
                    #sizes.append( len(sliced_data) )
                    
                    ##if sliced_data[0] == 2115331200:
                    ##    a = 0
                        
            #d = np.array(data)
            #print(time.time()-tt)            
            dummy_write = self.write_merged_new( var=v, table = 'observations_table', data=d)
            
            ### writing RECORD INDEX 
            if v == 'date_time':
                self.date_time = d

        ### build and write observations_id
        fixed_size = 20  # empty spaced will be filled with zeroes so that they will all have same length TO DO 
        
        #obs_id = range(len(data))
        year = self.current_year
        ref = int(year) * 10000000000000000
        obs_id_resized = np.array(np.arange(len(d), dtype=np.int64) + ref, dtype='S')
        #obs_id_resized = []
        #tt = time.time()
        #for i in obs_id:
            #zeroes = fixed_size - len(year) - len(str(i))
            #obs_id = year+ '0'*zeroes +str(i)
            #obs_id_resized.append(np.bytes_(obs_id))
        #print(time.time()-tt)
        
        dummy_write = self.write_merged_new( var='observation_id', table = 'observations_table', data=obs_id_resized)
        all_variables.extend(['observation_id'])
        
        ### build report_id
        rep_id_obs=[] # must be as long as observations_table i.e. one per observation
        rep_id_head=[] # must be as long as header table i.e. one per record
        source_ids = [] # storing source_id 
        policies = []
        
        for i,ts in enumerate(all_ts):
            
            best_file = res[ts]['best_file']
            best_ds = res[ts]['best_ds']
            
            ind_min = self.all_timestamps_dic[ts][best_ds][best_file][0]
            ind_max = self.all_timestamps_dic[ts][best_ds][best_file][1]    
            
            record_size = sizes[i]
            zeroes = fixed_size - len(year) - len(str(i))
            rep_id =  year+ '0'*zeroes +str(i)
            rep_id_head.append(rep_id)            
            
            rep_id_obs.extend( record_size*[rep_id] )
            source_ids.extend( record_size*[np.bytes_(best_ds)]  )
            policies.extend( record_size* [res[ts]['policy'] ] )
            
        self.header_table_rep_id =  np.array( rep_id_head).astype(np.bytes_)
        
        dummy_write = self.write_merged_new( var='report_id', table = 'observations_table', data=np.array(rep_id_obs).astype(np.bytes_))
        all_variables.extend(['report_id'])
        del rep_id_obs
        
        ### build source_id
        dummy_write = self.write_merged_new( var='source_id', table = 'observations_table', data=np.array(source_ids).astype(np.bytes_))
        all_variables.extend(['source_id'])
        
        ### build data policy license 
        dummy_write = self.write_merged_new( var='data_policy_licence', table = 'observations_table', data=np.array(policies))
        all_variables.extend(['data_policy_licence'])
        del policies 
        
        # ### build data value_significance         
        ### 12	Instantaneous value of observed parameter
        ### https://github.com/glamod/common_data_model/blob/master/tables/observation_value_significance.dat
        
        # all_variables.extend(['value_significance'])
        all_variables.extend(['value_significance'])        
        sig = np.full (  len(d), 12 ) 
        dummy_write = self.write_merged_new( var='value_significance', table = 'observations_table', data=np.array(sig))
        
        ### build conversion_flag, method  ## TO DO TODO can be implemented in harvester already 
        variables.extend(['conversion_flag' , 'conversion_method'])
        conv = np.full (  len(d), 2 ) 
        dummy_write = self.write_merged_new( var='conversion_flag', table = 'observations_table', data=np.array(conv))
        meth = np.full (len(d), np.nan) 
        dummy_write = self.write_merged_new( var='conversion_method', table = 'observations_table', data=np.array(meth))
        
        all_variables.extend(['conversion_flag','conversion_method'])
        del conv, meth
        
        ### advanced_assimilation_feedback         
        ass = [ 1 if s in ['era5_1' , 'era5_1_mobile' , 'era5_2' , 'era5_2_mobile'] else 0  for s in source_ids  ] 
        dummy_write = self.write_merged_new( var='advanced_assimilation_feedback', table = 'observations_table', data=np.array(ass))
        del ass, source_ids 
        all_variables.extend(['advanced_assimilation_feedback' ])
        
        ### Writing remainig of the untouched i.e. void variables 
        missing_cdm_var = [ v for v in self.dic_type_attributes['observations_table'].keys() if v not in all_variables]
        for v in missing_cdm_var:
            if v in all_variables:
                continue
            try:
                n = self.get_null( self.encodings['observations_table'][v]['type']  )
            except:
                n = np.nan 
            da = np.full( len(d), n) 
            dummy_write = self.write_merged_new( var=v, table = 'observations_table', data=np.array(da))                    
        
        self.record_size = sizes
        print(time.time()-tt)
        
    def write_merged_new(self, var='', table = '', data=''):
        """ Module to write the output file as netCDF.
        NB a table is made of only one variable"""

        out_name = self.out_file
        df_variable_type = self.df_variable_type  ### NB variable types are stored as strings 
        
        
        ### write 'observations_table','header_table','era5fb', 'station_configuration' tables 
        if table in ['observations_table','header_table','era5fb', 'station_configuration']:
                try:
                    descr   = bytes( self.dic_type_attributes[table][var]['description']    , 'utf-8' )
                except:
                    descr    = bytes( 'missing'    , 'utf-8' )
                    #print(' FFF FAILING WITH DESCRIPTION: ', var , ' ' ,  self.dic_type_attributes[table][var]['description']) 
                try:
                    ext_tab = bytes( self.dic_type_attributes[table][var]['external_table'] , 'utf-8' )
                except:
                    ext_tab = bytes( 'missing' , 'utf-8' )
                    #print(' FFF FAILING WITH EXTERNAL TABLE : ', var )                                                         

                
                ''' trying to convert the variable types to the correct types stored as attribute, read from the numpy dic file '''
                try:
                    var_type = df_variable_type[df_variable_type.variable == var]['type'].values[0] 
                    if  str(  type(data[0]) ) == var_type:
                        pass
                    else:
                        if type(data[0]) in [str, bytes, np.bytes_, np.str_] :
                            pass
                except:
                    pass
                    #print ('FAILED converting column ' , var, ' type ', type(data[0]) , )
                    
        
                if 'collection' in var:
                    a=0
                dic = {var: np.array(data)}  # making a 1 colum dictionary to write 
                if type(data) == np.dtype('O'):
                    if type(data[0]) in [float, np.float64, np.float32]:
                        
                        dic = {var: np.array(data, dtype=np.float64)}  # making a 1 colum dictionary to write
                    elif  type(data[0]) in [str]:
                        dic = {var: np.array(data, dtype=np.dtype('S1'))}  # making a 1 colum dictionary to write
                            
                    #print(var, type(data[0]))
                #print('SHAPE IS FFF ', table[k].shape )
                try:
                    #write_dict_h5(out_name, dic , table, self.encodings[table], var_selection=[], mode='a', attrs = {'description':descr, 'external_table':ext_tab}  )  # TO DO HERE TODO CHANGE write_dict_h5_old or write_dict_h5
                    write_dict_h5_old(out_name, dic , table, self.encodings[table], var_selection=[], mode='a', attrs = {var: {'description':descr, 'external_table':ext_tab}}  )  
                    
                except:
                    print("+++ FAILED table " ,  table , '  ' , var )
                    print(var, type(data[0]))
                    a=0
                    
        elif table == 'cdm_tables':
            for k in self.data['cdm_tables'].keys():
                table = self.data['cdm_tables'][k]
                if len(table) > 0:
                    
                    for t in table.keys():
                        table[t]['description'] = np.string_(table[t].description)
                        try:
                            
                            if table[t].external_table == ' ':
                            
                                del table[t].attrs['external_table']
                        except:
                            pass

                    table.to_netcdf(out_name, format='netCDF4', engine='h5netcdf', mode='a', group = k)
                #a=0
                #logging.info('Writing the cdm table %s to the netCDF output ', k)

        #elif table == 'source_configuration':  
        #    data.to_netcdf(out_name, format='netCDF4', engine='h5netcdf', mode='a', group = table)
        #    #logging.info('Writing the source_configuration table to the netCDF output ')
        
        elif table == 'source_configuration':  
            dic = {var: data} 
            write_dict_h5_old(out_name, dic , table, self.encodings[table], var_selection=[], mode='a', attrs = {'description':'Filename for data from source', 'external_table':''}  )  

        elif table == 'record_index':  
            data.to_netcdf(out_name, format='netCDF4', engine='h5netcdf', mode='a')
            

    def merge(self ,  datasets, station ):        

        """ Call to the merge_all_data() and write_merged_file() methods """
        
        mode='test'
        if mode == "test":      
            a = self.initialize_data( datasets, station ) # reading the input files 
            dummy = self.merge_all_data() 
            #dummy = self.create_merged_tables()
            
            #logging.info('*** Finished merging, now writing the output netCDF file ***' )         
            #a = self.write_merged_file()
            #logging.info('*** Done writing the output ! ***')

            self.write_merging_summary()
            return True

        else:
            try:
                a = self.initialize_data( station = station, datasets = datasets ) # reading the input files 
                dummy = self.merge_all_data()            
                #logging.info('*** Finished merging, now writing the output netCDF file ***' )         
                #a = self.write_merged_file()
                #logging.info('*** Done writing the output ! ***')
                self.write_merging_summary()                        
                return True          
            
            except MemoryError:
                o = open("FAILED_MERGING_LIST.txt", 'a+')                          
                print('Failed: ' , station )
                o.write(station + '\n' )
                self.write_merging_summary()                        
                return False 

    def write_merging_summary(self ):
        
        return 0

# 

def make_merge_list_old(harvested_base_dir, merged_out_dir, kind=''):
    """ Utility to summarize the content of the harvested directory.
    Kind is a flag for the four types of stations:
    - 'regular' i.e. normal stations 
    - 'mobile' i.e. known moving stations , WIGOS  0-20999-0-xxxxx
    - 'NoCoordMatchNoIDMatch' i.e. completely unidentified stations , WIGOS  0-20888-0-xxxxx
    - 'noDistMatch' i.e. it exixts a valid WIGOS id but it is not compatible with the coordinates, WIGOS 0-20777-0-xxxxx  
    - 'inconsistentCoordinates' i.e. stations with problems in lat and lon but not moving stations , WIGOS  0-20666-0-xxxxx
    """
    
    print('+++ Creating merged list from the harvested directory +++ ')
    dic = { 'station' : [] , 
                 'dataset' : [] ,
                 'file' : [],
                 'years' : [],
                 'MAXsize(MB)' : [] }
    
    
    for db in tqdm(os.listdir(harvested_base_dir) ) :
        
        if kind== 'mobile':
            if db not in ['era5_1_mobile' , 'era5_2_mobile' , 'igra2']:
                continue 
        
        if 'logs' in db:
            continue
        
        stations = os.listdir( harvested_base_dir + '/' + db )
        
        for station in stations:
            
            ### filter station type 
            
            if kind == 'regular':
                if '0-20666-' in station or  '0-20777-'  in station or   '0-20888-'  in station:
                    continue
            
            elif kind == 'mobile':
                if '0-20999-' not in station:
                    continue                
            
            elif kind == 'orphan':
                if '0-20666-' not in station or  '0-20777-'  not in station or   '0-20888-'  not in station:
                        continue            
            # extract the file with the correctly processed years per file 
            if '.dat' in station:
                continue
            
            files = [f for f in os.listdir(harvested_base_dir + '/' + db  + '/' + station ) if 'harvested' in f and 'processed' not in f ]
            
            years=[]
            sizes = []
            
            for f in files:

                station = f.split('_')[0]
                year = f.split('_')[1].replace('_','') 
                years.append(year)

                file_size = os.stat( harvested_base_dir + '/' + db  + '/' + station +'/' + f).st_size / (1000*1000)  ### bytes to MB
                sizes.append(file_size)
                
            dic['station'].append(station)
            dic['years'].append( ','.join(years) )
            dic['MAXsize(MB)'].append(max(sizes))
            dic['dataset'].append(db)
            dic['file'].append(files[0])
            
    df = pd.DataFrame.from_dict(dic)
    df = df.sort_values(by=['station' , 'dataset' ])
    df.to_csv(merged_out_dir + '/logs/merging_list_from_harvested_dir.csv' , sep = '\t')
    
    return list( np.unique( df.station.values ) ) 



def make_merge_list(data_directories, merged_out_dir, kind=''):
    """ Utility to summarize the content of the harvested directory.
    Kind is a flag for the four types of stations:
    - 'regular' i.e. normal stations 
    - 'mobile' i.e. known moving stations , WIGOS  0-20999-0-xxxxx
    - 'NoCoordMatchNoIDMatch' i.e. completely unidentified stations , WIGOS  0-20888-0-xxxxx
    - 'noDistMatch' i.e. it exixts a valid WIGOS id but it is not compatible with the coordinates, WIGOS 0-20777-0-xxxxx  
    - 'inconsistentCoordinates' i.e. stations with problems in lat and lon but not moving stations , WIGOS  0-20666-0-xxxxx
    """
    
    print('+++ Creating merged list from the harvested directory +++ ')
    dic = { 'station' : [] , 
                 'dataset' : [] ,
                 'file' : [],
                 'years' : [],
                 'MAXsize(MB)' : [] }
    
    
    if kind == 'mobile':
        datasets = [ f for f in data_directories.keys() if f in ['era5_1_mobile' , 'era5_2_mobile' , 'npsound', 'igra2', 'shipsound'] ]
    else:
        datasets = [ f for f in data_directories.keys() if f not in ['era5_1_mobile' , 'era5_2_mobile' , 'npsound', 'igra2', 'shipsound'] ]
        #datasets = [  f for f in data_directories.keys() if f in ['era5_1' , 'era5_2' ,
                                                                  #'era5_1759' , 'era5_1761',
                                                                  #'igra2',
                                                                  #'ncar',
                                                                  #'hara' , 'giub' , 'amma' , ]  ]
        
    for db in tqdm( datasets ) :
        if not os.path.isdir( data_directories[db] ):
            continue # not all datasets e.g. GIUB have orhpahs 
        stations = os.listdir( data_directories[db] )
        for station in stations:
            ### filter station type 
            if kind == 'regular':
                if '0-20666-' in station or  '0-20777-'  in station or   '0-20888-'  in station or '0-20999' in station:
                    continue
            
            elif kind == 'mobile':
                if '0-20999-' not in station:
                    continue                
                
            elif kind in ['orphan','orphans']:
                if not ( '0-20666-' in station or  '0-20777-'  in station or   '0-20888-'  in station):
                        continue        
            # extract the file with the correctly processed years per file 
            if '.dat' in station:
                continue
            
            files = [f for f in os.listdir( data_directories[db] + '/' + station ) if 'harvested' in f and 'processed' not in f ]
            if len(files)==0:
                continue
            
            years=[]
            sizes = []
            
            for f in files:

                #station = f.split(db)[0][-4:]
                if 'restored' in db:
                    station = f.split('_')[0]
                else:
                    station = f.split(db)[0][:-6]
                
                year = f.split('_')[1].replace('_','') 
                years.append(year)

                try:
                    
                    file_size = os.stat( data_directories[db]  + '/' + station +'/' + f).st_size / (1000*1000)  ### bytes to MB
                    sizes.append(file_size)
                except:
                    print (data_directories[db]  + '/' + station +'/' + f, 'not found')
                    sizes.append(0)
                
            dic['station'].append(station)
            dic['years'].append( ','.join(years) )
            dic['MAXsize(MB)'].append(max(sizes))
            dic['dataset'].append(db)
            dic['file'].append(files[0])
            
    df = pd.DataFrame.from_dict(dic)
    df = df.sort_values(by=['station' , 'dataset' ])
    if not os.path.isdir(merged_out_dir + '/logs'):
        os.makedirs(merged_out_dir + '/logs')
    df.to_csv(merged_out_dir + '/logs/' + kind + '_merging_list_from_harvested_dir.csv' , sep = '\t')
    
    return list( np.unique( df.station.values ) ) 
    
def create_stat_summary(data_directories, stat_id):
    """ Looks in the database for files matching with the given station id.
    Returns a dictionary with the files per dataset
    """
    station_dic = { 'year':[], 'dataset':[] , 'files': [] }
    
    ### loop over possible years
    for y in range(1880, 2030):
        y = str(y)    
        
        for d,i in data_directories.items():
            print(d, i)
            if not os.path.isdir(i):
                continue
            stations = os.listdir(i) 
            if stat_id not in stations:
                continue
            else:
                try:
                    
                    harvested_files = [f for f in os.listdir(i+'/'+stat_id) if '.nc' in f and y in f.split(stat_id+'_')[1].split('_')[0]  ]
                except:
                    harvested_files = []
                
                
                ### filter giub wrong files

                    
                
                if len(harvested_files) >0:
                    
                    
                    if d == 'giub':
                        
                        black_list_giub = open( 'giub_problematic_file.txt' , 'r').readlines()
                        black_list_giub = [ l.replace('\n','').split('/')[-1] for l in black_list_giub ] 
                        harvested_files = [f for f in harvested_files if f not in black_list_giub ]
                        a=0
                    #if y not in station_dic.keys():
                    #    station_dic = {'year':[], 'dataset':[] , 'files': [] }
                    
                    ''' # old 
                    station_dic['files'].append(i+'/' +stat_id + '/' + ','.join(harvested_files)  ) 
                    station_dic['year'].append(y) 
                    station_dic['dataset'].append(d) 
                    '''
                    for f in harvested_files:
                        station_dic['files'].append(i+'/' +stat_id + '/' + f  ) 
                        station_dic['year'].append(y) 
                        station_dic['dataset'].append(d)       
                        
    df = pd.DataFrame.from_dict(station_dic)
    
    return df

                        
       

    
    
    
    


#from merging_yearly_parameters import harvested_base_dir, merged_out_dir, data_directories, run_exception 
from merging_yearly_parameters import  merged_out_dir,  data_directories, run_exception, POOL, pool_number, check_missing_stations
from merging_yearly_parameters import   add_sensor, stations, station_kind, create_merging_list

run_mode = ''

def run_wrapper(data_directories, run_exception, station):
    
    station_df = create_stat_summary(data_directories, station)

    if run_exception:  # try to run whatever working stations 
        try:
            a = Merging.merge(station_df, station) # single station dictionary
        except:
            a = open('FAILED_files.txt' , 'a+')
            a.write(station + '\n')
            a.close()
            pass
        
    else:
        try:
            
            a = Merging.merge(station_df, station ) # single station dictionary
        except MemoryError as e:
            print(e)
            raise ValueError(station)
        
ray_run_wrapper = ray.remote(run_wrapper)

def check_missing(merged_dir, skip_completed=False):
    """ Creates a table of the already merged files """
    
    dic = {}
    stations = [s for s in os.listdir(merged_dir) if 'logs' not in s ]
    for s in stations:
        completed = merged_dir + '/' +s + '/completed.txt' 
        if os.path.isfile( completed ):
            lines = open(completed).readlines()
            
            dic[s] = [ y.replace('\n', '') for y in lines ]
        else:
            pass
        
    fully_completed=[]
    if skip_completed:
        for s in dic.keys():
            if 'completed' in dic[s]:
                fully_completed.append(s)
            
    return dic, fully_completed
            
    

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Make CDM compliant netCDFs")
    parser.add_argument('--min_year' , '-min_y', 
                    help="Earliest year to merge"  ,
                    type = str,
                    default = '1880')

    parser.add_argument('--max_year' , '-max_y', 
                    help="Latest year to merge"  ,
                    type = str,
                    default = '2023')
    
    args = parser.parse_args()
    min_year = args.min_year 
    max_year = args.max_year
    
    try:
        min_year = int(min_year)
        max_year = int(max_year)
    except:
        min_year = '1880'
        max_year = '2023'
        print('Please check the range of the years to merge - Setting to default minimum and maximum values')

    year_range = [ int(min_year) , int(max_year) ]
        
        
    """ Initialize the merging class """
    ### TO BE FIXED 
    if station_kind in ['mobile', 'orphans','orphan']  :
        merged_out_dir = merged_out_dir + '_' + station_kind 
    
    
    #import glob
    #xstations = glob.glob('/users/staff/uvoggenberger/scratch/harvest_20240403/hara/*')
    #stations = [x.split('/')[-1] for x in xstations]
    #for s in stations:
        #try:
            
            #os.remove('/mnt/users/scratch/leo/scratch/FH/MERGED_YEARLY_0MAR2024_REGULAR/'+s+'/completed.txt')
            #print(s, 'removed')
        #except:
            #print(s, 'not found')

    
    ### Saving list of file to merge, create list of stations 
    if create_merging_list:
        all_stat = make_merge_list(data_directories, merged_out_dir, kind=station_kind )
    else:
        try:
            df =  pd.read_csv( merged_out_dir + '/logs/' + station_kind + '_merging_list_from_harvested_dir.csv', sep = '\t')
            all_stat =  list( np.unique( df.station.values ) )              
        except:
            print('Merging list file not found, will create it +++ ')
            all_stat = make_merge_list(data_directories, merged_out_dir, kind=station_kind )
    
    #stations = ['0-20000-0-82930']

    # select from parameters file or extracted list
    if not stations:
        stations = all_stat
        
    pool_number = 50
    #stations = stations[0:500]
 
    skip_completed = False ### set to FALSE to rerun the station, will check each year afterwards
    
    if check_missing_stations:
        processed_stats, fully_completed = check_missing(merged_out_dir, skip_completed=skip_completed)
        stations = [p for p in all_stat if p not in fully_completed ]
    else:
        processed_stats = {}
            
        
    #stations = stations[4501:]
    Merging = Merger(add_sensor=add_sensor, 
                     out_dir = merged_out_dir,
                     min_year= min_year,
                     max_year=max_year,
                     processed_stats=processed_stats)

    run_exception = False
    print('run exception is ', run_exception )
        
    ### Here: must create a list of stations 
    # stations = ['0-20001-0-11035', '0-20001-0-10393' , '0-20000-0-70219' , '0-20000-0-06610']  # 0-20000-0-71879 , 0-20000-0-82900                                                                                                                 
    # stations = ['0-20001-0-10393']
    #stations = [s for s in stations if s in os.listdir('/scratch/das/federico/HARVEST_YEARLY_22FEB2024_amma/amma')]
    #stations = ['0-20000-0-45004']
    
    ray.init(num_cpus=20) # runtime_env={"working_dir": "/users/staff/uvoggenberger/CEUAS/CEUAS/public/merge"}, 

    POOL = True
    #stations = stations[0:1]
    #for bad in ['0-20000-0-78970', '0-20000-0-93891']:
       #del stations[stations.index(bad)]
    
    #stations = ['0-20001-0-11035']
    #stations = ['0-20666-0-20353',  '0-20666-0-25428']
    
    #import glob
    #xstations = glob.glob('/users/staff/uvoggenberger/scratch/harvest_20240403/hara/*')
    #stations = [x.split('/')[-1] for x in xstations]
    #for s in stations:
        #os.remove('/mnt/users/scratch/leo/scratch/FH/MERGED_YEARLY_0MAR2024_REGULAR/s/completed.txt')
    
    #stations = ['0-20000-0-71938']
    if False:
        import glob
        rstations = glob.glob('/mnt/users/scratch/leo/scratch/FH_orphan/*/restored/*')
        stations = [os.path.basename(r) for r in rstations]
        
    if False:
        import glob
        if 'fixedgiub' in data_directories['giub']:
            fixeddirs = glob.glob(data_directories['giub']+'/*/*.nc')
            if len(fixeddirs) != 0:
                import reharvest_giub
                #reharvest_giub.fix(data_directories['igra2'].split('igra2')[0]+'giub')
                fixeddirs = glob.glob(data_directories['igra2'].split('igra2')[0]+'giub'+'/*/*.nc')
                print('fixed', len(fixeddirs), 'giub files')
            for d in fixeddirs:
                sid = d.split('/')[-2]
                yr = d.split(sid)[-1][1:5]

                x = glob.glob(merged_out_dir+'/'+sid+'/*.txt')
                if len(x) == 1:     
                    os.remove(x[0])
                    print(sid, yr)
                else:
                    print(yr)
                    #print('old giub', sid, 'not deleted')
    #i = stations.index('0-200data_directories['giub']00-0-80413')
    #del stations[i]
    stations = ['0-20000-0-60769']
    try :
        
        import glob
        os.remove(glob.glob(merged_out_dir+'/'+stations[0]+'/*.txt')[0])
    except:
        pass
    POOL = True
    if len(stations)== 1:
        POOL = False
        
    if POOL:
        #p=Pool(pool_number)
        #func=partial(run_wrapper, data_directories,  run_exception )
        #dummy=p.map(func,stations)
        rp_data_directories = ray.put(data_directories)
        result_ids = []
        for station_i in stations:
            result_ids.append(ray_run_wrapper.remote(rp_data_directories,  run_exception, station_i))
        results = ray.get(result_ids)
        ray.shutdown()
                
    else:
        for station in stations[:1]:
            dummy=run_wrapper(data_directories, run_exception, station)
    '''
    else:
        for station in stations:
            station_df = create_stat_summary(station, data_directories)
    
            if run_exception:
                try:
                    a = Merging.merge(station,  datasets = station_df , mode = run_mode ) # single station dictionary
                except:
                    print(' ERROR, SKIPPED !!!')
                    out = open('Failed_merged.txt' , 'a+')                    
                    out.write( str(station) + '\n')
            else:
                a = Merging.merge(station,  datasets = station_df , mode = run_mode ) # single station dictionary                
    
            print(' --------------------------------> finished station ', station )
    '''



### MAURITIUS /users/staff/uvoggenberger/scratch/mauritius_data/temp 
