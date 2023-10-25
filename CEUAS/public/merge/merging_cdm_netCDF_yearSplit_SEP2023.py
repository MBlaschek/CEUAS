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


#from numba import njit
import psutil
import copy
from numba import njit
from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
import warnings

warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)

import code

sys.path.append('../harvest/code_cop2')
sys.path.append('../postprocess/add_sensor')

from add_sensor_to_merged import Sensor, datetime_toseconds, wrapper, MergedFile
from harvest_convert_to_netCDF_yearSplit import write_dict_h5, clean_station_configuration 

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

    def __init__(self, out_dir = 'output'  ):
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
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
        self.variable_types = {}
        
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
        self.index_offset = 0 # will be replaced when running 
        self.hour_time_delta = 2 # decide up to which time shift in HOURS separate records are considered identical  

        self.only_std_plevels = False  # set to True to store only standard pressure level data 
        self.std_plevs    = [1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 85000, 92500, 100000]

        self.add_sensor = True
            
        self.copy = True # make a copy of the merged file before adding the sensor. Argument for the add_sensor wrapper function  


    def initialize_data(self , station = '', datasets = {} ):
        """ Initialize dataset; store relevant data as attributes.
                   Args ::     dic{}  datasets (dataframe with year, dataset and file path 
                                   e.g. tateno = { 'ncar'    : 'example_stations/ncar/chuadb_trhc_47646.txt.nc'   ,
                                                           'igra2'   : 'example_stations/igra2/chJAM00047646-data.txt.nc'  } 
        """       
        self.datasets          = datasets
        self.datasets_keys = datasets.keys()
        self.station             = station
        
        self.counting_record = 1 # will be used to build record_id in header table
        self.counting_observations = 1
        
        # list of meaningful variables in the observations_table, all others are empty/nans filled
        
        self.observations_table_vars = ['date_time', 'z_coordinate' , 'z_coordinate_type', 'observed_variable' , 
                                                           'observation_value' , 'report_id' , 'observation_id' , 'latitude' , 'longitude',
                                                           'units', 'source_id', 'data_policy_licence', 'observation_duration', 'value_significance',
                                                           'advanced_assimilation_feedback' , 'advanced_uncertainty']
                
        """ Loading the econding of the tables created from the harvester script and to be applied again """
        self.encodings = np.load('groups_encodings.npy' , allow_pickle = True ).item()
        self.encodings['era5fb'] = np.load('era5fb_encodings_all.npy' , allow_pickle = True ).item()            
        self.dic_type_attributes = np.load('dic_type_attributes.npy',allow_pickle= True).item()

        self.era5fb_columns = self.dic_type_attributes['era5fb'].keys()
        self.header_columns = self.dic_type_attributes['header_table'].keys()
        #self.statconf_columns = self.dic_type_attributes['station_configuration'].keys()

        self.empty_cdm_var = [ v for v in self.dic_type_attributes['observations_table'].keys() if v not in self.observations_table_vars ]  # variables to be filled with nans   with proper data type         


        self.obstab_nans_filled = False      
        self.obs_in_header = {} # dic placeholder for data from obs_tab to be kept for header_table
        self.fill_cdm()
        
    def fill_cdm(self):
            
        self.data['cdm_tables'] = {}              

        """ Loop over all the datasets                                                                                                                                     
                k: name of the dataset                                                                                                                    
                v: list of file paths,    eg 'era5_1':[filepath_1, filepath_2 ] """      
        
        #######                                           
        # STATION CONFIGURATION
        #######   
        # Read the station primary_id, and find the entry in the CUON dataset (or global dataset)
        stat_conf = pd.read_csv("CUON_station_configuration.csv", sep = '\t')
        self.stat_conf_CUON = stat_conf 
        s = stat_conf.loc[ stat_conf.primary_id == self.station ]   
        
        if s.empty:
            a = open('logs/failed_stat_conf.txt' , 'a+')
            a.write(self.station + '\n')
        else:
            for k in s.columns:
                if pd.isnull(s[k].values[0]):
                    try:
                        v = s[k].astype(np.float64)            
                    except:
                        v = b'NA'           
                    s[k] = v 
                    
        if '20999' in station:
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
                    cdm = xr.open_dataset(self.datasets.files.values[0] , engine = 'h5netcdf' , group = t )
                    self.data['cdm_tables'][t] = cdm 
            except:
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
        
        
        
    '''
    def read_obs_data(self):
        """ Reading observations and header table data """
            
        ### Observations_table
        for k,v in self.datasets.items() :
                data[k] = {}
                for F in v:
    
                    logging.info(' Dataset ::: *** %s %s ' , k , F   )                  
    
                    data[k][F] = {}
    
                    h5py_file = h5py.File(F, 'r')
                    data[k][F]['h5py_file'] = h5py_file 
    
                    a = h5py_file['recordtimestamp']
    
                    data[k][F]['recordtimestamp'] = pd.to_datetime( a[:], unit='s',  origin=pd.Timestamp('1900-01-01') )
                    data[k][F]['recordindex']         = h5py_file['recordindex']
                    data[k][F]['dateindex']            = h5py_file['dateindex']
                    
                    data[k][F]['observations_table'] = h5py_file['observations_table']
                    data[k][F]['header_table']           = h5py_file['header_table']
                    
                    a = h5py_file['recordtimestamp']
                    data[k][F]['max_date'] = max(a)
                    data[k][F]['min_date']  = min(a)
    
                    data[k][F]['counter'] = 0


                #######
                # HEADER TABLE
                #######
                head_tab = h5py_file['header_table']
                logging.info('*** header_table')
                data[k][F]['header_table'] = {}
                for var in head_tab.keys():
                    if ('string' in var or 'hdrlen' in var or 'duplicates' in var ): continue
                    try:                                  
                        data[k][F]['header_table'][var] = (np.array(head_tab[var][:])).astype(self.dic_type_attributes['header_table'][var]['type'] )
                    except:
                        print('failed convertion type header' , k , ' ' , F , ' ' ,  var )


                #######                                                                                    
                # SOURCE CONFIGURATION                                             
                #######      
                d = xr.open_dataset(F , engine = 'h5netcdf' , group = 'source_configuration' , decode_times = False )
                data[k][F]['source_configuration'] = d
                logging.debug('Done with %s source_configuration' , str(k) )
                d.close()
        self.data = data
        
    '''

    '''
    def make_quick_df(self, ds, file, ind_low, ind_sup, indices=None ):
        """ Quickly create pandas dataframe (debug purpose) """
        
        data=self.data[ds][file]['observations_table']
        variables = ['date_time', 'z_coordinate', 'observed_variable', 'observation_value']
        data_dic = {}
        if not indices:
            indices = range(ind_low, ind_sup)
        for v in variables:
            d = data[v][indices]
            if v == 'date_time':
                d = pd.to_datetime ( pd.to_datetime( d[:], unit='s',  origin=pd.Timestamp('1900-01-01') ) )
            data_dic[v] = d 
        
        df = pd.DataFrame.from_dict( data_dic )
        
        return df
    '''
        
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
            
            if ds == 'era5_1759':
                a = 0
            #last_index = 0
            for F in v:                 
                print('FILE ' , F )
                data = dic[k][F]  # h5py open file 
                self.unique_dates[k][F] = {}
                self.unique_dates[k][F]['datetimes'] = {}
                
                timestamps = data['recordtimestamp'][:]
                #indices = data['recordindex'][:]
                indices_inf = data['recordindex'][:] # starting index of the record
                indices_sup = [ indices_inf[1:][i] for i in  range(len(indices_inf[:-1])) ] # ending index of the record
                indices_sup.append(999999) # dummy high value, to be used with last record of the list 
                
                ### quick check
                # df_check = pd.DataFrame( {'ts':timestamps , 'inf': indices_inf, 'sup':indices_sup } ) 
                
                for ts,inf,sup, index in zip(timestamps, indices_inf, indices_sup, range(len(timestamps)) ) :
                    if ts == 2114769600:
                        a = 0
                    if ts not in all_timestamps_dic.keys():
                        all_timestamps_dic[ts] = {}
                    if ds not in all_timestamps_dic[ts].keys():
                        all_timestamps_dic[ts][ds] = {}   
                              
                    all_timestamps_dic[ts][ds][F] = [inf, sup, index]  # store the records limits in the obesrvation_table and the record_index in the headet_table

                #timestamps_conv = pd.to_datetime( timestamps, unit='s',  origin=pd.Timestamp('1900-01-01') )
                 
        #unique_timestamps = list(np.unique(list(all_timestamps_dic.keys() ))).sort()
        self.all_timestamps_dic = all_timestamps_dic
                    
    
    def reduce_timestamps(self):
        """ Simplify/reduce all timestamps by flagging possible duplicates """
        
        #1 check duplicates 
        time_delta = self.hour_time_delta * 60*60 # timestamps are still in seconds, so must convert hours in self.hour_time_delta to seconds 
        unique_timestamps = list(np.unique(list(self.all_timestamps_dic.keys() )))
        unique_timestamps.sort()
        
        duplicated_ts = []
        unique_ts = [] 
        for i in range(1, len(unique_timestamps)):

            if (unique_timestamps[i] - unique_timestamps[i-1]) < time_delta :
                print('Must check duplicated timestamps',  unique_timestamps[i] )
                duplicated_ts.append( [unique_timestamps[i-1] ,  unique_timestamps[i] ] )
                #unique_ts.pop()
                #duplicated_ts.append( unique_timestamps[i-1] )                                      
                #duplicated_ts.append( unique_timestamps[i]    )
            else:
                if i ==1:  # since loop starts from 2nd item, must add the 1st timestamp if the 2nd is not duplicated 
                    unique_ts.append( unique_timestamps[i-1]  )
                unique_ts.append( unique_timestamps[i]  )
        
        #if len(duplicated_ts) >0:
        #    pass

        return unique_ts, duplicated_ts    
    
        
    def extract_record_data(self, dt, ds, file ):
        
        """ Extracting the length of the temp and wind valid observations, and the maximum height (or min pressure) 
        TODO must check what happens with z_coordinate type in height or gph """
        
        ###TO DO CAN LOAD OBS TAB adnd HEAD TAB in memory so you only slice and not read everytime 
        
        if dt == 2114769600: #2115331200
            a = 0 
            
        h5_file = self.dic_h5py[ds][file]
        ind_min, ind_max = self.all_timestamps_dic[dt][ds][file][0] ,  self.all_timestamps_dic[dt][ds][file][1] 
        
        ot = h5_file['observations_table']
        
        temp_ind = np.where(ot['observed_variable'][ind_min:ind_max] == 126)[0] # temperature
        wind_ind = np.where(ot['observed_variable'][ind_min:ind_max]  == 107)[0] # wind speed 
        gph_ind = np.where(ot['observed_variable'][ind_min:ind_max]  == 117)[0] # geopotential 
        
        # length of valid data for temp and wind 
        temp_values = ot['observation_value'][ind_min:ind_max][temp_ind]
        num_valid_temp = len(np.unique(temp_values[~np.isnan(temp_values)]))
        
        wind_values = ot['observation_value'][ind_min:ind_max][wind_ind]
        num_valid_wind = len(np.unique(wind_values[~np.isnan(wind_values)]))
        
        #checking pressure or height or gph 
        press_temp_ind = ot['z_coordinate'][ind_min:ind_max][temp_ind]
        press_wind_ind = ot['z_coordinate'][ind_min:ind_max][wind_ind]
        
        if len(gph_ind) > 0:
            
            max_gph = max( ot['observation_value'][ind_min:ind_max][gph_ind] )
        else:
            max_gph = -999999
            
        if len(temp_ind) >0:
            min_press_temp = min( ot['z_coordinate'][ind_min:ind_max][temp_ind] )  
            max_press_temp =  max( ot['z_coordinate'][ind_min:ind_max][temp_ind] ) 
        else:
            min_press_temp, max_press_temp = 999999, -999999 
            
        if len(wind_ind) >0:
            min_press_wind =  min( ot['z_coordinate'][ind_min:ind_max][wind_ind] ) 
            max_press_wind =  max( ot['z_coordinate'][ind_min:ind_max][wind_ind] )
            
        else:
            min_press_wind, max_press_wind = 999999, -999999        

        
        pressure_ind = np.where(ot['z_coordinate_type'][ind_min:ind_max] == 1 )[0] #
        
        if len(pressure_ind) >0:
            max_pressure = max( ot['z_coordinate'][ind_min:ind_max][pressure_ind] )
            min_pressure = min( ot['z_coordinate'][ind_min:ind_max][pressure_ind] )
            
        else:
            max_pressure = -999999
            min_pressure = 999999
            
        height_ind = np.where(ot['z_coordinate_type'][ind_min:ind_max] == 0 )[0] #        
        if len(height_ind) >0:
            max_height = max( ot['z_coordinate'][ind_min:ind_max][height_ind] )
            min_height = min( ot['z_coordinate'][ind_min:ind_max][height_ind] )
            
        else:
            max_height = -999999
            min_height = 999999
            
        if len(gph_ind) >0:
            max_gph = max( ot['observation_value'][ind_min:ind_max][gph_ind] )
            min_gph = min( ot['observation_value'][ind_min:ind_max][gph_ind] )
            
        else:
            max_gph = -999999                    
            min_gph = 999999
            
        return [num_valid_temp, num_valid_wind], [min_press_temp, min_press_wind, max_press_temp, max_press_wind], [min_pressure, min_gph, min_height], [max_pressure, max_gph, max_height]
    
    
    
    
    
    
    def find_best_record(self, ts, extract_record_data =True):
        """ Effective merging procedure, applying hierarchical selection rules 
                   ts: timestamp in seconds after 1900 
                   datasets: list of datasets per record e.g. ['era5_1', 'igra2' , 'ncar', ... ] """
        
        ### Variable placeholders
        best_ds = False
        best_file = False
 
        if ts == 2115331200:
            a = 0
        if ts == 2114769600:
            a = 0
            

        # alla available datasets
        datasets = [ l for l in list(self.all_timestamps_dic[ts].keys() ) if 'extract_record_data' not in l ]   # might remove the datasets as a variable since it is in an attirbute of the class already. Think about it....
        
        ### Setting data policy (might be redundant)
        ### see https://github.com/glamod/common_data_model/blob/master/tables/data_policy_licence.dat        
        if 'ncar' in datasets or  'igra2' in datasets or  'era5_1759' in datasets or 'era5_1761' in datasets or 'giub' in datasets:
            data_policy = 0
        else:
            data_policy = 4        
        
        # list of all era5 datasets to be preferred except special cases 
        era5_ds = [ f for f in datasets if f in ['era5_1', 'era5_1_mobile' , 'era5_2' , 'era5_2_mobile'] ]

        
        '''
        ### 1) giub over era5 if year < 1950
        ###  1950-01-01-00:00:00 = 1577836800 seconds
        ###  pd.to_datetime( a[:], unit='s',  origin=pd.Timestamp('1900-01-01') )
        '''
        
        ### placeholder for best data 
        extracted_best_record = {}
        
        ### only one dataset available 
        if len(datasets) == 1:
            best_ds = datasets[0]
            files =  list(self.all_timestamps_dic[ts][best_ds].keys() )
            if len(files) ==1:  ### must only loop over all files 
                best_file = files[0]
                if extract_record_data:
                    # [num_valid_temp, num_valid_wind], [min_press_temp, min_press_wind, max_press_temp, max_press_wind], [max_gph, max_pressure, max_height]
                    valid_t_w, min_t_w_max_t_w, min_p_h_gph, max_p_h_gph  = self.extract_record_data(ts, best_ds, best_file)  # TO REMOVE, not needed 
                    extracted_best_record[best_file] = [best_ds, valid_t_w, min_t_w_max_t_w, min_p_h_gph, max_p_h_gph]
                    
            else:
                for f in files:
                    valid_t_w, min_t_w_max_t_w, min_p_h_gph, max_p_h_gph  = self.extract_record_data(ts, best_ds, f)   
                    extracted_best_record[f] = [ best_ds, valid_t_w, min_t_w_max_t_w, min_p_h_gph, max_p_h_gph] 
                
            
        else:  # hierarchical choices, according to datasets and files 
            extract_record_data = True
            ### first find best datasets, then check length of files 
            if ts < 1577836800:
                if 'giub' in datasets:
                    best_ds = 'giub'
                
            else:
                if len(era5_ds) >0 and 'igra2' not in datasets:
                    best_ds = era5_ds[0]
                elif len(era5_ds) ==0 and 'igra2' in datasets:
                    best_ds = 'igra2'
                else:
                    best_ds = False #to be determined  
                
                #else:  # no era5, no igra => must check length of files EXCEPT ncar to be taken as last resort due to wind problem (some data in knots, some other in m/s without clear rule )
                #     valid_t_w, min_t_w_max_t_w, min_p_h_gph, max_p_h_gph = self.extract_record_data(ts, best_ds, f)
                    
                if best_ds:
                    files =  list(self.all_timestamps_dic[ts][best_ds].keys() )
                    for f in files:
                        valid_t_w, min_t_w_max_t_w, min_p_h_gph, max_p_h_gph  = self.extract_record_data(ts, best_ds, f)  # TO REMOVE, not needed                     
                        extracted_best_record[f] = [best_ds, valid_t_w, min_t_w_max_t_w, min_p_h_gph, max_p_h_gph] 
                    
                else:
                    for d in datasets:
                        files =  list(self.all_timestamps_dic[ts][d].keys() )
                        for f in files:
                            valid_t_w, min_t_w_max_t_w, min_p_h_gph, max_p_h_gph  = self.extract_record_data(ts, d, f)  # TO REMOVE, not needed                     
                            extracted_best_record[f] = [d, valid_t_w, min_t_w_max_t_w, min_p_h_gph, max_p_h_gph] 

        #return best_ds, best_file, data_policy, highest_level, total_temp_records, total_wind_records 
                        
                        
        ### HERE: analyze extracted_best_record dictionary with data 
                        
        if not best_file:
            if len(extracted_best_record.keys()) ==1:
                best_file = list(extracted_best_record.keys() )[0]
                best_ds = extracted_best_record[best_file][0]
            else:  # loop and determine best record 
                
                max_height = -1
                min_pressure = 999999
                best_file = ''
                best_ds = ''
                
                for file in list(extracted_best_record.keys() ) :                    
                    # [num_valid_temp, num_valid_wind], [min_press_temp, min_press_wind, max_press_temp, max_press_wind], [min_pressure, min_gph, min_height], [max_pressure, max_gph, max_height]
                    # where [dataset, valid_t_w, min_t_w_max_t_w, min_p_h_gph, max_p_h_gph]
                    
                    #current_sum = extracted_best_record[file][1][0] + extracted_best_record[file][1][1]  # [d, valid_t_w, min_t_w_max_t_w, min_p_h_gph, max_p_h_gph] 
                    current_max_height = extracted_best_record[file][4][1]
                    current_min_pressure = extracted_best_record[file][3][0]
                    
                    if  current_min_pressure < min_pressure:  # current pressure is lower than previous
                        best_file = file
                        best_ds = extracted_best_record[file][0]
                        min_pressure = current_min_pressure
                        #max_height = current_max_height 
                        
                    if  current_max_height >= max_height and max_height >0 and  not (current_min_pressure <= min_pressure): # current pressure higher than previous (pressure and height should behave oppositely)
                        a = 0
                        print('CHECK WRONG')
                        # must be something wrong with pressure and height ??? 

        self.all_timestamps_dic[ts]['extract_record_data'] =  extracted_best_record[best_file]
        

            
        return best_ds, best_file, data_policy 
        
        
        
    def merge_timestamp(self):
        """ Apply merging selection criteria to each reduced timestamp (no more duplicated timestamp)"""
        
        # unique ts without duplicates, duplicated timestamps
        unique_ts, duplicated_ts = self.reduce_timestamps()
        # all available ts (unique or possibly duplicated)
        all_timestamps = list(self.all_timestamps_dic.keys())
        
        #all_timestamps = all_timestamps[:300]  # speed it up TO DO CHANGE !!!!! 
        all_timestamps.sort()

        # container for processed ts 
        processed_timestamps = []
        all_combined_timestamps ={}        
            
        # dictionary holding all the data for a given timestamp 
        keys = [ 'timestamp' , 'policy' , 'all_duplicated_dt', 'all_duplicated_datasets', 'all_duplicated_files', 'best_ds' , 'best_file', 'ind_inf', 'ind_sup'] 
        all_combined_timestamp_data = {}
        for k in keys:
            all_combined_timestamp_data[k] = []
            
        for dt,index in zip(all_timestamps, range(len(all_timestamps)) ) :
            
            if dt == 2115331200: ### TODO for debugging  2114769600
                a = 0
                
            if dt == 2114769600: ### TODO for debugging  2114769600
                    a = 0                
            # already processed timestamp 
            if dt in processed_timestamps[-5:]:   # all timestamps are sorted so should not check the entire list 
                continue 
            # temporary dictionary holding info on duplicates 

            # all_combined_timestamps[dt]['indices'] = []
            # NO all_combined_timestamps[dt]['index_header'] = []

            # if the time delta is less than the specified limit in self.hour_time_delta
            # here: remove duplicated timestamps 
            #datasets = [ l for l in list(self.all_timestamps_dic[dt].keys() ) if 'extract_record_data' not in l ]   # might remove the datasets as a variable since it is in an attirbute of the class already. Think about it....
            
            if dt in unique_ts: # no time duplicate detected, apply standard merging procedure 
                
                real_time = dt 
                duplicated_time = dt 
                
                best_ds, best_file, policy = self.find_best_record(dt)
                
                processed_timestamps.append(dt)
                
            else:
                possible_duplicates = [ p for p in duplicated_ts if p[0] == dt or p[1] == dt ][0]
                # apply hierarchical selection 
                for t in possible_duplicates:
                    
                    #two possible timestamps, analyze each then compare best ds and best file of each 
                    ### earlier timestamp
                    low_ts = possible_duplicates[0]
                    datasets_low = list(self.all_timestamps_dic[low_ts].keys() ) 
                    best_ds_low, best_file_low, policy = self.find_best_record(low_ts)
                    
                    ### later timestamp
                    high_ts = possible_duplicates[1]
                    datasets_high = list(self.all_timestamps_dic[high_ts].keys() ) 
                    best_ds_high, best_file_high, policy = self.find_best_record(high_ts)
                    
                                        
                    ### valid_t_w, min_t_w_max_t_w, min_p_h_gph, max_p_h_gph
                    ### example: 
                    ### self.all_timestamps_dic[low_ts]['extract_record_data'] = ['igra2', [0, 2], [1961.33, 1961.33, 3432.3274, 3432.3274], [999999, 1961.33, 999999], [0, 3432.3274, 0]]
                    
                    #1. check total temp+wind record length
                    l_low =  self.all_timestamps_dic[low_ts]['extract_record_data'][1][0] + self.all_timestamps_dic[low_ts]['extract_record_data'][1][0]
                    l_high =  self.all_timestamps_dic[high_ts]['extract_record_data'][1][0] + self.all_timestamps_dic[high_ts]['extract_record_data'][1][0]
                    
                    min_p_low =  self.all_timestamps_dic[low_ts]['extract_record_data'][3][0] # minimum pressure
                    min_p_high = self.all_timestamps_dic[high_ts]['extract_record_data'][3][0]
                    
                    max_h_low =  self.all_timestamps_dic[low_ts]['extract_record_data'][4][1] # maximum height
                    max_h_high = self.all_timestamps_dic[high_ts]['extract_record_data'][4][1]
                    
                    # selecting best ds by number of records 
                    if ( l_low > l_high ):
                        best_ds = best_ds_low
                        best_file = best_file_low
                        real_time = low_ts
                    else:
                        best_ds = best_ds_high
                        best_file = best_file_high                        
                        real_time = high_ts
                    
                    # selecting best ds by lowest pressure                      
                    if min_p_low < min_p_high:
                        best_ds = best_ds_low
                        best_file = best_file_low   
                        real_time = low_ts
                        
                    else:
                        best_ds = best_ds_high
                        best_file = best_file_high                                  
                        real_time = high_ts

                    # selecting best ds by highest height            
                    if  max_h_low > max_h_high: # current pressure higher than previous (pressure and height should go together)
                        best_ds = best_ds_low
                        best_file = best_file_low    
                        real_time = low_ts
                    
                    else:
                        best_ds = best_ds_high
                        best_file = best_file_high                           
                        real_time = high_ts
                                        
                    duplicated_time = [t for t in possible_duplicates if t != real_time ][0]
                    
                    processed_timestamps.append(low_ts)
                    processed_timestamps.append(high_ts)


            ### Saving the extracted data             
            all_combined_timestamps[real_time]= {}            
            all_combined_timestamps[real_time]['policy'] = policy
            all_combined_timestamps[real_time]['best_ds'] = best_ds
            all_combined_timestamps[real_time]['best_file'] = best_file
            all_combined_timestamps[real_time]['all_duplicated_files'] = 0
            all_combined_timestamps[real_time]['all_duplicated_records'] = 0
            all_combined_timestamps[real_time]['real_time'] = real_time
            all_combined_timestamps[real_time]['duplicated_time'] = duplicated_time
            
        self.merged_timestamp = all_combined_timestamps
        
        

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
    

    def make_header_var(self, tab='', var='' ):
        """  """
        combined_records = self.best_record  
        # extract datasets
        datasets = [ d for d in combined_records.keys() if 'length' not in d and 'policy' not in d]
        
        tipo = self.dic_type_attributes[tab][var]['type']                   
        void = self.get_null(tipo)
        
        if var in ['latitude' , 'longitude', 'report_id']:
            vector = self.obs_in_header[var][self.record_index]
            return vector 
            
        if var in ['duplicates', 'source_id']:
            vector = np.full( (combined_records['total_length_headertab']) , b'NA                                                          ' ) # filling a vector with dataset id i.e. source_id field
        else:
            vector = np.full( (combined_records['total_length_headertab']) , void ) # filling a vector with dataset id i.e. source_id field

        for ds in datasets:
            for file in combined_records[ds].keys():
                
                if var in ['duplicates']:
                    duplicates = combined_records[ds][file]['duplicates']
                    np.put (vector , combined_records[ds][file]['all_inserting_indices_header_tab'], duplicates)                        
                else:
                    reading_indices = combined_records[ds][file]['all_reading_indices_header_tab']                    
                    sliced_data = self.data[ds][file]['h5py_file'][tab][var][:][reading_indices]  
                    if var == 'source_id':
                        sliced_data = [b''.join(f) for f in sliced_data ] 
                    np.put (vector , combined_records[ds][file]['all_inserting_indices_header_tab'], sliced_data)                        
                    
        return vector
    

    def merge_all_data(self):       
        """ Construct a dictionary with all the dataframes of each dataset, either reading it from saved pickles or reading it from memory """

        logging.info('***** Starting the merging process merge_all_data')

        # avoidable loops, should not matter much  # self.all_years:  
        #for this_year in [1981, 1982, 1983]:  # loop over all available years in the date_times 
        #for this_year in self.all_years:   # loop over all available years in the date_times 
        
        ### here: extract all available years
        
        #for this_year in self.all_years:   # loop over all available years in the date_times
        a = 0
        data_all_years = self.extract_file_per_year()
        all_years = list(data_all_years.keys() )
        
        #a = self.initialize_all_data(data_per_year)
        
        #all_years = [y for y in all_years if int(y) > 1978 ]
        #all_years = ['1967']
        for this_year in all_years:   # loop over all available years in the date_times 
            print('=== Running year ::: ' , this_year )
            self.current_year = this_year 
            
            data_this_year = data_all_years[this_year]
            
            following_year = str(int(this_year) + 1) 
            if following_year in all_years:
                data_following_year = data_all_years[following_year]
            else:
                data_following_year = ''
            
            
            dummy = self.open_data(data_this_year)
            """ Making all date_time for the selected year  """

            dummy = self.make_unique_datetime()            

            
            # principal merging algoritm, applying selection rules 
            merge = self.merge_timestamp() # selection of the best data 
            
            dummy = self.initialize_out_file()
            dummy = self.make_merged_observation_table()
            dummy = self.make_merged_feedback_table()
            dummy = self.make_sourceconf_table()
            
            dummy = self.make_merged_header_sourceconf_table()
            dummy = self.make_standard_cdm_table()
            
            a = 'HERE'  # TO DO HERE
        
            """
            ### WRITING RECORD TIMESTAMPS AND INDEX
            di=xr.Dataset()
            datetimes, recordindex = np.unique( date_time, return_index=True )
            di['recordtimestamp'] = ( {'recordtimestamp' : datetimes.shape } , datetimes )
            di['recordtimestamp'].attrs['units']='seconds since 1900-01-01 00:00:00'
            di['recordindex']          = ( {'recordindex' : recordindex.shape } ,  recordindex )

            self.write_merged_new(var='record_index', table='record_index', data=di)                      
            
            ### WRITING STATION CONFIGURATION TABLE 
            #for d in self.statconf_columns:
            for d in self.statconf_columns:
                try:
                    data = self.make_statconf_var(tab='station_configuration', var=d)
                    self.write_merged_new(var = d, table = 'station_configuration', data=data )    
                    #print('done with station_configuration variable ' , d )
                except:
                    pass
                    #print('wrong variable ' , d )
            
            ### Adding sensor_id to observations_table 
            if self.add_sensor:
                print('*** Adding sensor *** ')            
                add_sensor = wrapper(out_dir = self.out_dir , station_id = self.station.split('-')[-1] , file = self.out_name , copy = self.copy )
                print('*** Added sensor *** ')
            else:
                pass
            

        """

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
            files.append(np.bytes_(best_file))
        
        self.write_merged_new(var='source_file', table = 'source_configuration', data=np.array(files))
        
        
    def make_merged_feedback_table(self, out_file='' ):
        """ Create merged era5fb table """
        ### FEEDBACK TABLE
        res = self.merged_timestamp
        
        # all date_time for this year 
        all_ts = list(res.keys() )
        
        # variables to read from the observations_table of each file 
        variables = self.dic_type_attributes['era5fb'].keys()
                
        # here, only extract files that are actually selected as best_file (will not pre-load the others)
        all_files = [ f for f in list( np.unique( [ res[dt]['best_file'] for dt in all_ts ] ) ) if 'era5_1_' in f or 'era5_2_' in f ]

        for v in variables:
            # keep full data in memory, for better efficiency 
            load_full_data = {}
            for f in all_files:
                if v in h5py.File(f, 'r')['era5fb'].keys():
                    load_full_data[f] = h5py.File(f, 'r')['era5fb'][v][:]
                else:
                    n = self.get_null( self.encodings['era5fb'][v]['type']  )
                    load_full_data[f] = np.full( len( h5py.File(f, 'r')['observations_table']['date_time']), n)      
            
            data = []
        
            for ts in all_ts :
                
                best_file = res[ts]['best_file']
                best_ds = res[ts]['best_ds']
                
                ind_min = self.all_timestamps_dic[ts][best_ds][best_file][0]
                ind_max = self.all_timestamps_dic[ts][best_ds][best_file][1]
            
                if 'era5_1_' in best_ds or 'era5_2_' in best_ds:
                    sliced_data = load_full_data[f][ind_min:ind_max]
                else:
                    a = self.get_null(tipo= self.dic_type_attributes['era5fb'][v])
                    sliced_data = np.full( (ind_max-ind_min ) , a )
                    
                data.extend(sliced_data)

            d = np.array(data)
            dummy_write = self.write_merged_new( var=v, table = 'era5fb', data=d)
                                                 
                                                 
    def make_merged_header_sourceconf_table(self, out_file='' ):
        """ Create merged era5fb table """ 
        a = self.header_table_rep_id
      
        ### HEADER TABLE
        res = self.merged_timestamp
        
        # all date_time for this year 
        all_ts = list(res.keys() )
        
        # variables to read from the observations_table of each file 
        other_variables = ['report_id' , 'station_record_number' , 'source_id' , 'duplicate_status' , 'duplicates' ,]
        
        variables = [ v for v in self.dic_type_attributes['header_table'].keys() if v not in other_variables  and v != 'primary_station_id' ]
                
        # here, only extract files that are actually selected as best_file (will not pre-load the others)
        all_files = list( np.unique( [ res[dt]['best_file'] for dt in all_ts ] ) )

        for v in variables:
            # keep full data in memory, for better efficiency 
            load_full_data = {}
            for f in all_files:
                load_full_data[f] = {}
                if v in h5py.File(f, 'r')['header_table'].keys():
                    load_full_data[f][v] = h5py.File(f, 'r')['header_table'][v][:]
                else:
                    n = self.get_null( self.encodings['header_table'][v]['type']  )
                    load_full_data[f][v] = np.full( len( h5py.File(f, 'r')['header_table']['report_timestamp']), n)      
            
            data = []
            source_files = []
            
            for ts in all_ts :
                
                best_file = res[ts]['best_file']
                source_files.append(best_file)
                
                best_ds = res[ts]['best_ds']
                
                #ind_min = self.all_timestamps_dic[ts][best_ds][best_file][0]
                #ind_max = self.all_timestamps_dic[ts][best_ds][best_file][1] # not useful in this case 
                index =  self.all_timestamps_dic[ts][best_ds][best_file][2]
                value = load_full_data[best_file][v][index]
                    
                data.append(value)
                    
            d = np.array(data)
            dummy_write = self.write_merged_new( var=v, table = 'header_table', data=d)
              
        # station_id
        stat_ids = np.full( (len(data)), np.bytes_(self.station) )
        dummy_write = self.write_merged_new( var='station_primary_id', table = 'header_table', data= stat_ids )
        
        # report_id
        dummy_write = self.write_merged_new( var='report_id', table = 'header_table', data=np.array(self.header_table_rep_id ) )
        
        #self source_configuration_files = source_files
        #### source_configuration file 
        #dummy_write = self.write_merged_new( var='source_files', table = 'source_configuration', data=np.array(source_files ) )
        
              
    def  make_merged_observation_table(self, out_file='' ):
        """ Create merged observations,header tables """
        
        res = self.merged_timestamp
        
        # all date_time for this year 
        all_ts = list(res.keys() )
        
        # variables to read from the observations_table of each file 
        all_variables = []
        
        variables = ['date_time' , 'z_coordinate_type', 'z_coordinate' , 'observed_variable', 'observation_value' ,'longitude', 'latitude' , 'original_units']
        all_variables.extend(variables)
        
        # here, only extract files that are actually selected as best_file (will not pre-load the others)
        all_files = list( np.unique( [ res[dt]['best_file'] for dt in all_ts ] ) ) 
        
        # loop over all variables, will extract and write one by one 
        
        ### save to write later
        
        # pre-loop to check if there are duplicated pressure values within one record
        for ts in all_ts :
            
            best_file = res[ts]['best_file']
            best_ds = res[ts]['best_ds']
            
            ind_min = self.all_timestamps_dic[ts][best_ds][best_file][0]
            ind_max = self.all_timestamps_dic[ts][best_ds][best_file][1]
        
            obs_val = h5py.File(best_file, 'r')['observations_table']['observation_value'][ind_min:ind_max]
            obs_var = h5py.File(best_file, 'r')['observations_table']['observed_variable'][ind_min:ind_max]
            z_coord =  h5py.File(best_file, 'r')['observations_table']['z_coordinate'][ind_min:ind_max]
            
            unique_var = np.unique(obs_var)
            for v in unique_var:
                obs_var_ind = np.where( obs_var == v )[0]
                if len( np.unique(z_coord[obs_var_ind])) == len(obs_var_ind):
                    pass
                else:
                    dic = {}
                    for v in ['observation_value' , 'observed_variable' , 'date_time' , 'z_coordinate_type' , 'z_coordinate']:
                        dic[v] = h5py.File(best_file, 'r')['observations_table'][v][ind_min:ind_max]
                    df = pd.DataFrame.from_dict(dic)
                    b = 0 ### here must check how to remove duplicates 

        sizes = []
        
        for v in variables:
            
            # keep full data in memory, for better efficiency 
            load_full_data = {}
            for f in all_files:
                load_full_data[f] = h5py.File(f, 'r')['observations_table'][v][:]
            
            data = []
        
            for ts in all_ts :
                
                best_file = res[ts]['best_file']
                best_ds = res[ts]['best_ds']
                
                ind_min = self.all_timestamps_dic[ts][best_ds][best_file][0]
                ind_max = self.all_timestamps_dic[ts][best_ds][best_file][1]
            
                sliced_data = list(load_full_data[best_file][ind_min:ind_max])
                data.extend(sliced_data)
                
                if v == 'date_time':
                    sizes.append( len(sliced_data) )
                    
                    if sliced_data[0] == 2115331200:
                        a = 0
                        
            d = np.array(data)
            dummy_write = self.write_merged_new( var=v, table = 'observations_table', data=d)

        ### build and write observations_id
        fixed_size = 20  # empty spaced will be filled with zeroes so that they will all have same length TO DO 
        
        obs_id = range(len(data))
        year = self.current_year
        obs_id_resized = []
        for i in obs_id:
            zeroes = fixed_size - len(year) - len(str(i))
            obs_id = year+ '0'*zeroes +str(i)
            obs_id_resized.append(obs_id)
        
        dummy_write = self.write_merged_new( var='observation_id', table = 'observations_table', data=np.array(obs_id_resized))
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
            
        self.header_table_rep_id =  rep_id_head 
        
        dummy_write = self.write_merged_new( var='report_id', table = 'observations_table', data=np.array(rep_id_obs))
        all_variables.extend(['report_id'])
        del rep_id_obs
        
        ### build source_id
        dummy_write = self.write_merged_new( var='source_id', table = 'observations_table', data=np.array(source_ids))
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
        sig = np.full (  len(data), 12 ) 
        dummy_write = self.write_merged_new( var='value_significance', table = 'observations_table', data=np.array(sig))
        
        ### build conversion_flag, method  ### TODO can be implemented in harvester already 
        variables.extend(['conversion_flag' , 'convertion_method'])
        conv = np.full (  len(data), 2 ) 
        dummy_write = self.write_merged_new( var='conversion_flag', table = 'observations_table', data=np.array(conv))
        meth = np.full (len(data), np.nan) 
        dummy_write = self.write_merged_new( var='convertion_method', table = 'observations_table', data=np.array(meth))
        
        all_variables.extend(['conversion_flag','convertion_method'])
        del conv, meth
        
        ### advanced_assimilation_feedback         
        ass = [ 1 if s in ['era5_1' , 'era5_1_mobile' , 'era5_2' , 'era5_2_mobile'] else 0  for s in source_ids  ] 
        dummy_write = self.write_merged_new( var='advanced_assimilation_feedback', table = 'observations_table', data=np.array(ass))
        del ass, source_ids 
        all_variables.extend(['advanced_assimilation_feedback' ])

        self.header_table_rep_id =  rep_id_head # to be inserted in header_table
        
        ### Writing remainig of the untouched i.e. void variables 
        missing_cdm_var = [ v for v in self.dic_type_attributes['observations_table'].keys() if v not in all_variables]
        for v in missing_cdm_var:
            if v in all_variables:
                continue
            try:
                n = self.get_null( self.encodings['observations_table'][v]['type']  )
            except:
                n = np.nan 
            da = np.full( len(data), n) 
            dummy_write = self.write_merged_new( var=v, table = 'observations_table', data=np.array(da))                    
        
        
    def write_merged_new(self, var='', table = '', data=''):
        """ Module to write the output file as netCDF.
        NB a table is made of only one variable"""

        out_name = self.out_file

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
                    
                try:
                    var_type = self.dic_type_attributes[table][var]['type']
                except:
                    var_type = np.int32

                ''' trying to convert the variable types to the correct types stored as attribute, read from the numpy dic file '''
                if type(table[0]) != var_type:
                    if type(table[0]) == np.bytes_:
                        pass
                    else:
                        try:
                            data = data.astype( var_type ) 
                        except:
                            print ('FAILED converting column ' , var, ' type ', type(data[0]) , ' to type ', var_type )
    
                dic = {var: data}  # making a 1 colum dictionary to write 
                #print('SHAPE IS FFF ', table[k].shape )
                write_dict_h5(out_name, dic , table, self.encodings[table], var_selection=[], mode='a', attrs = {'description':descr, 'external_table':ext_tab}  )  

        elif table == 'cdm_tables':
            for k in self.data['cdm_tables'].keys():
                table = self.data['cdm_tables'][k]
                table.to_netcdf(out_name, format='netCDF4', engine='h5netcdf', mode='a', group = k)
                #logging.info('Writing the cdm table %s to the netCDF output ', k)

        #elif table == 'source_configuration':  
        #    data.to_netcdf(out_name, format='netCDF4', engine='h5netcdf', mode='a', group = table)
        #    #logging.info('Writing the source_configuration table to the netCDF output ')
        
        elif table == 'source_configuration':  
            dic = {var: data} 
            write_dict_h5(out_name, dic , table, self.encodings[table], var_selection=[], mode='a', attrs = {'description':'Filename for data from source', 'external_table':''}  ) 

        elif table == 'record_index':  
            data.to_netcdf(out_name, format='netCDF4', engine='h5netcdf', mode='a')
            




            
    def merge(self ,  station = '' , datasets = '' , mode = 'test'):        

        """ Call to the merge_all_data() and write_merged_file() methods """
        if mode == "test":      
            a = self.initialize_data( station = station, datasets = datasets ) # reading the input files 
            dummy = self.merge_all_data() 
            dummy = self.create_merged_tables()
            
            #logging.info('*** Finished merging, now writing the output netCDF file ***' )         
            #a = self.write_merged_file()
            #logging.info('*** Done writing the output ! ***')
            #self.write_merging_summary()
            return True

        else:
            o = open("FAILED_MERGING_LIST.txt", 'a+')          
            try:
                a = self.initialize_data( station = station, datasets = datasets ) # reading the input files 
                dummy = self.merge_all_data()            
                #logging.info('*** Finished merging, now writing the output netCDF file ***' )         
                #a = self.write_merged_file()
                #logging.info('*** Done writing the output ! ***')
                self.write_merging_summary()                        
                return True          
            except MemoryError:
                print('Failed: ' , station )
                o.write(station + '\n' )
                self.write_merging_summary()                        
                return False 

    def write_merging_summary(self ):
        
        return 0
        '''
        out = open(self.out_dir + '/' + self.station + '/merging_summary.txt', 'a+')
        
        out.write(self.out_name + '\n')
        for d in self.datasets.keys():
            out.write(d + '\t')
            for f in self.datasets[d]:
                out.write(f + '\t')
            out.write('\n')
        '''

'''
# can be changed singularly if needed 
data_directories   = { 'era5_1'       : harvested_base_dir + '/era5_1' ,
                                   'era5_1_mobile'       : harvested_base_dir + '/era5_1_mobile' ,
                                   'era5_2'       : harvested_base_dir + '/era5_2' ,
                                   'era5_3188' : harvested_base_dir + '/era5_3188' ,
                                   'era5_1759' : harvested_base_dir + '/era5_1759' ,
                                   'era5_1761' : harvested_base_dir + '/era5_1761' ,
                                   'ncar'           : harvested_base_dir + '/ncar' ,
                                   'igra2'          : harvested_base_dir + '/igra2' ,
                                   'bufr'            : harvested_base_dir + '/bufr' , 
                                   'amma'        : harvested_base_dir + '/amma' ,
                                   }



'''


def create_stat_summary(stat_id, data_directories):
    """ Looks in the database for files matching with the given station id.
    Returns a dictionary with the files per dataset
    """
    station_dic = { 'year':[], 'dataset':[] , 'files': [] }
    
    ### loop over possible years
    for y in range(1880, 2030):
        y = str(y)    
        
        for d,i in data_directories.items():
            if not os.path.isdir(i):
                continue
            stations = os.listdir(i) 
            if stat_id not in stations:
                continue
            else:
                harvested_files = [f for f in os.listdir(i+'/'+stat_id) if '.nc' in f and y in f.split(stat_id)[1]  ]
                if len(harvested_files) >0:
                    
                    #if y not in station_dic.keys():
                    #    station_dic = {'year':[], 'dataset':[] , 'files': [] }
                    
                    station_dic['files'].append(i+'/' +stat_id + '/' + ','.join(harvested_files)  ) 
                    station_dic['year'].append(y) 
                    station_dic['dataset'].append(d) 
                
    df = pd.DataFrame.from_dict(station_dic)
    
    return df

                        
       

    
    
    
    


from merging_yearly_parameters import harvested_base_dir, merged_out_dir, data_directories, run_exception 

kind = ''
run_mode = ''


#b = 0
#a = create_stat_summary('0-20000-0-63661', data_directories) 




if __name__ == '__main__':

    """ Initialize the merging class """
    ### TO BE FIXED 
    if kind == 'mobile':
        merged_out_dir = merged_out_dir + '_mobile'
        
    Merging = Merger(merged_out_dir)

    run_exception = False
    print('run exception is ', run_exception )

    
    # here: must create stations list 
    stations = ['0-20001-0-11035']  # 0-20000-0-71879 , 0-20000-0-82900
    """
     for g in os.listdir('.'):
    ...:     files = os.listdir(g)
    ...:     for f in files:
    ...:         if f not in a.keys():
    ...:             a[f] = []
    ...:         a[f].append(g)

    find stations in multiple datasets in harvested directory
    """
    for station in stations:
        #print(station)

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
              try:
                    print(' I am running a total compressed size of ::: ', (size / 10**9) , '  GB   since minimum allowed:' , min_size , '  and maximum allowed: ', max_size , '   station: ' , station ) 
                    a = Merging.merge(station,  datasets = station_dic , mode = run_mode ) # single station dictionary
              except:
                    print('*** Failed :::  ', station , ' ::: must re-do ***')
                    out.write(station + '\n' )
              '''


# run: 
# -s 0-20000-0-82930 -l 100 -m 1000000000000 
# 0-20000-0-82930,0-20000-0-03976 
# -s 0-20000-0-97760 -l 100 -m 1000000000000 


# 0-20999-0-99007 
### 34737


### MAURITIUS /users/staff/uvoggenberger/scratch/mauritius_data/temp 
