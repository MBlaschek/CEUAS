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
from harvest_convert_to_netCDF import write_dict_h5, clean_station_configuration 

# nan int = -2147483648 
#from harvest_convert_to_netCDF import datetime_toseconds   # importing the function to write files with h5py 


"""
@njit
def replace_global(replaced_indices, replaced_vector, replacing_vector):
    #FAILED Attempt to make vector item replacement faster 
    for insert_index, replacing_value in zip(replaced_indices, replacing_vector):
        replaced_vector[insert_index] = replacing_value 
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
                   Args ::     dic{}  datasets (dictionary where keys are the dataset names e.g. bufr, igra2 etc. , and the value is the path to the corresponding netCDF file 
                                   e.g. tateno = { 'ncar'    : 'example_stations/ncar/chuadb_trhc_47646.txt.nc'   ,
                                                           'igra2'   : 'example_stations/igra2/chJAM00047646-data.txt.nc'  } 
        """       
        self.datasets          = datasets
        self.datasets_keys = datasets.keys()
        self.station             = station
        
        self.counting_record = 1 # will be used to build record_id in header table
        self.counting_observations = 1
        
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
        
        data['cdm_tables'] = {}              

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
            
        for c in sc.columns:
            if 'Unnamed' in c:
                sc = sc.drop(columns=c)
            
        data['station_configuration'] = sc
                        
        self.statconf_columns = sc.keys()
        
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

            data['cdm_tables'] = {}
            """ Reading the CDM tables that do not depend on specific stations or observations (fixed values), for the first file only """
            for t in self.standard_cdm: # [ 'crs' , 'observed_variable', 'units' , 'z_coordinate_type' , 'station_type', 'station_configuration_codes']              
                if t not in data['cdm_tables'].keys():
                    #data['cdm_tables'][t] = ''
                    cdm = xr.open_dataset(F , engine = 'h5netcdf' , group = t )
                    data['cdm_tables'][t] = cdm 

        print(blue + 'Memory used after reading data: ', process.memory_info().rss/1000000000 , cend)

        self.data = data

        """ Making all date_times  """
        self.make_unique_datetime()

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
        
    def get_null(self, tipo = ''):
        ''' Simply returns the proper format for ''null' value '''        
        if tipo in  [np.int32, np.int64]  :
            void = -2147483648
        elif tipo == np.float32 :
            void = np.nan
        elif tipo == np.bytes_ :
            void = b'nan'
        else:
            return np.nan 
        return void
        
        
    def make_unique_datetime(self):
        """ Building the global set of date_times and indices from the various datasets. 
              The datetimeindex is read from the original netCDF file. 
              Will compare the unique date_time of each dataset and extract the global unique date_times
              """

        logging.info('\n *** Running make_all_datetime ' )

        """ Loop over all the datasets 
                k: name of the dataset
                v: list of file paths,    eg 'era5_1':[filepath_1, filepath_2 ]"""

        all_timestamps = []
        
        for ds,files in self.datasets.items() :
            k,v = ds, files #rename
            self.unique_dates[k] = {}
            
            #last_index = 0
            for F in v: 
                print('FILE ' , F )
                data = self.data[k][F]
                self.unique_dates[k][F] = {}
                self.unique_dates[k][F]['datetimes'] = {}
                
                timestamps = data['recordtimestamp']
                ind_inf = data['recordindex']
                
                # inf and sup indices used when using list slices. It will not work if you try the access the last index as in "vector[last_index]" since it is out of boud. 
                # This will work: "vector[last_index-1:last_index]"
                len_obstab = len(data['observations_table']['date_time'][:])
                len_records =  len(data['recordtimestamp'])
                
                ind_sup = [ data['recordindex'][i+1] if i < ( len_records -1 ) else  len_obstab for i in range( len_records )  ]                
                
                for dt,inf,sup in zip(timestamps, ind_inf, ind_sup):
                    self.unique_dates[k][F]['datetimes'][dt] = {}
                    indices = [inf, sup, sup-inf] 
                    
                    ### filter data 
                    # remove internal duplicates i.e. data for the same records but multiple pressure levels
                    # remove also non finite plevel data (i.e. where plve == nan ) 
                    #ind_filtered = np.arange(inf, sup)
                    
                    plev =  self.data[k][F]['observations_table']['z_coordinate'][inf:sup]
                    plev_finite_ind = np.where(np.isfinite(plev) )[0]
                    obs_var = self.data[k][F]['observations_table']['observed_variable'][inf:sup][plev_finite_ind]
                    plev_finite = plev[plev_finite_ind]
                    
                    processed_p, indices_temp = [],[]
                    for value, press,index in zip(plev_finite, obs_var, plev_finite_ind):
                        if press in processed_p:
                            continue
                        if np.isfinite(value):
                            indices_temp.append(inf + index) # appending the lower index 
                            processed_p.append(press)
                    
                    '''
                    # TO DO IMPROVE SPEED, too slow 
                    df = pd.DataFrame.from_dict ( {'plev': plev , 
                                                                       'observed_variable':  self.data[k][F]['observations_table']['observed_variable'][plev_finite] , 
                                                                       'indices' : plev_finite    } )
                    
                    dfr =df.drop_duplicates (subset = ['plev', 'observed_variable'] , ignore_index=True )
                    # add_offset 
                    '''
                    
                    indices.append( indices_temp ) # append list of real indices and true length
                    indices.append( len(indices_temp) )
                    
                    ### TRY 
                    #inserting_indices =  list(range(len(indices_temp)))  # the new vector where the values are replaced is created from scratch, index starts from 0 each time 
                    inserting_indices = np.arange(len(indices_temp))
                    
                    indices.append( inserting_indices) # where to put the read data in the empty template vector 
                    # INDICES DESCRIPTION ::: indices = [inf, sup, total length, indices_to_read, length reduced indices, indices_to_insert ]
                    
                    ### end filtering data 
                    #self.unique_dates[k][F]['datetimes'][dt][0]
                    self.unique_dates[k][F]['datetimes'][dt]['indices'] = indices

                    #unfiltered = self.make_quick_df(  k, F, inf, sup, indices= list(range(sup-inf)))
                    #filtered = self.make_quick_df(  k, F, inf, sup, indices=None )
                    #last_index = last_index + len(inserting_indices) 
                    all_timestamps.append(dt)

        logging.info('\n *** Done calculating reduced indices ' )
                    
        # Storing all the info per timestamp, dataset, file 
        all_timestamps_dic = {} 
        all_timestamps = np.unique( all_timestamps )  
        all_timestamps.sort()
        
        ### build intermediate dictionary with dates 
        logging.info('\n *** Build intermediate dictionary with dates ' )

        # CHECK TIMESTAMP Timestamp('1979-06-20 12:00:00') 
        for dt in all_timestamps: # <- loop over all possible timestamps 
            all_timestamps_dic[dt] = {}
            for  ds,files in self.datasets.items():  # v is a list of files from the dataset
                k,v = ds, files #rename
                for F in v:
                    if dt in self.unique_dates[k][F]['datetimes'].keys():
                        
                        if k not in all_timestamps_dic[dt].keys():
                            all_timestamps_dic[dt][k]= {}
                        all_timestamps_dic[dt][k][F] = self.unique_dates[k][F]['datetimes'][dt]['indices']
                             
        logging.info('\n *** Done build intermediate dictionary with dates ' )
        
        ### flag for public data, extract duplicates
        processed_timestamps = []
        all_combined_timestamps ={}        
        
        for dt,index in zip(all_timestamps, range(len(all_timestamps)) ) :
            
            # will keep lowest timestamp of all possible duplicates i.e. only dt 
            
            #if dt == pd.Timestamp('1979-06-03 12:00:00'):
            #    print(0)

            if dt in processed_timestamps:
                continue 

            # temporary dictionary
            all_combined_timestamps[dt]= {}
            all_combined_timestamps[dt]['real_time']= ''
            
            all_combined_timestamps[dt]['policy'] = 0
            all_combined_timestamps[dt]['all_duplicated_dt'] = []
            all_combined_timestamps[dt]['all_duplicated_datasets'] = []
            all_combined_timestamps[dt]['best_ds'] = ''  
            all_combined_timestamps[dt]['best_file'] = ''  
            all_combined_timestamps[dt]['all_duplicated_files'] = []
            all_combined_timestamps[dt]['all_duplicated_records'] = []
            all_combined_timestamps[dt]['indices'] = []
            all_combined_timestamps[dt]['index_header'] = []
            
            # if the time delta is less than the specified limit in self.hour_time_delta
            possible_duplicate_timestamps = [d for d in all_timestamps[index+1:index+3] if d >=dt and abs(d-dt)< pd.Timedelta(hours= self.hour_time_delta)  ] 
            for t in possible_duplicate_timestamps:
                processed_timestamps.append(t)
            all_combined_timestamps[dt]['all_duplicated_dt'].extend(possible_duplicate_timestamps)
            
            # list with all possible duplicated date_times 
            possible_duplicate_timestamps.append(dt)
            all_datasets =  [ list(all_timestamps_dic[g].keys())[0] for g in possible_duplicate_timestamps ]

            all_combined_timestamps[dt]['all_duplicated_datasets'].append( np.bytes_([''.join(d) for d in all_datasets ][0]))
            
            if 'ncar' in all_datasets or  'igra2' in all_datasets or  'era5_1759' in all_datasets or 'era5_1761' in all_datasets:
                    all_combined_timestamps[dt]['policy'] = 4
            
            era5 = [ f for f in all_datasets if f == 'era5_1' or f== 'era5_1_mobile' or f=='era5_2' or f== 'era5_2_mobile' ]
            
            if len(era5) >0:
                best_ds = era5[0]  # only one case can happen, not more than 1 since era5 datasets are mutually excluding 
            elif 'igra2' in all_datasets:
                best_ds = 'igra2'
            elif 'ncar' in all_datasets:
                best_ds = 'ncar'
            else:
                best_ds = False # must look for longest records 
                

            for time in possible_duplicate_timestamps:
                best_dt, best_file, best_record = time,'',0        
                
                if best_ds and best_ds not in all_timestamps_dic[time].keys():
                    continue
                
                elif  best_ds and best_ds in all_timestamps_dic[time].keys():
                    
                    files = list(all_timestamps_dic[time][best_ds].keys())
                    all_combined_timestamps[dt]['all_duplicated_files'].extend(files)  
                    
                    for file in files:
                        records = all_timestamps_dic[time][best_ds][file][4]
                        if records >= best_record:
                            best_dt = time
                            best_file = file
                            best_record = records 
                            all_combined_timestamps[dt]['all_duplicated_records'].append(records)
                            all_combined_timestamps[dt]['all_duplicated_files'].append(file)  
                                
                elif not best_ds:
                    datasets = all_timestamps_dic[time].keys()
                    for ds in datasets:
                        files = all_timestamps_dic[time][ds]
                        all_combined_timestamps[dt]['all_duplicated_files'].extend(files)  
                        
                        for file in files:
                            records = all_timestamps_dic[time][ds][file][4]
                            if records >= best_record:
                                best_dt = time
                                best_file = file
                                best_record = records 
                                best_ds = ds 
                                all_combined_timestamps[dt]['all_duplicated_records'].append(records)
                                all_combined_timestamps[dt]['all_duplicated_files'].append(file)  

                all_combined_timestamps[dt]['best_ds'] = best_ds  
                all_combined_timestamps[dt]['best_file'] = best_file  
                all_combined_timestamps[dt]['real_time'] = best_dt # might not coincide 
                all_combined_timestamps[dt]['indices'] = all_timestamps_dic[best_dt][best_ds][best_file] # might not coincide 
                #if dt == pd.Timestamp('1979-07-27 11:00:00') or best_dt == pd.Timestamp('1979-07-27 11:00:00'):
                #    print(0)
                all_combined_timestamps[dt]['index_header'] = list(self.unique_dates[best_ds][best_file]['datetimes'].keys()).index(best_dt) # index of this merged time in the original file, for constructing the header table 
   
        # reloop to select candidates, duplicates 
        all_unique_dates = np.unique( list(all_combined_timestamps.keys()) )
        self.all_unique_dates = all_unique_dates
        self.best_duplicates_dic = all_combined_timestamps
        
        self.all_combined_timestamps = all_combined_timestamps
        
        all_years =  np.unique( pd.to_datetime( all_unique_dates ).year )
        self.all_years = all_years
        
        logging.debug('*** make_all_datetime finished ')         
        
        
    def get_null(self, tipo = ''):
        ''' Simply returns the proper format for ''null' value '''        
        if tipo == np.int32 or tipo == np.int64:
            void = -2147483648
        elif tipo == np.float32 or tipo == np.float64 :
            void = np.nan
        elif tipo == np.bytes_ :
            void = b'NA '
        return void


    def merge_all_timestamps(self, selected_timestamps_this_year ):
        """ Loop over the timestamps, takes the best record, store the indices to be read at the end of the loop """
        
        #all_timestamps = self.all_combined_timestamps
        #best = self.best_duplicates_dic         

        combined = {}
        
        # counter for the total length of the observations table 
        total_length_obstab,  total_length_headertab = 0,0

        # looping over the best file only 
        offset_obs, offset_head = 0,0 
        
        records_dic = {}
        keys = ['date_time', 'source_id' , 'file' , 'original_record_index' , 'policy' , 'duplicates' ]
        for k in keys:
            records_dic[k] = []

        
        
        
        for dt in selected_timestamps_this_year:
            best_ds =  self.best_duplicates_dic[dt]['best_ds'] 
            
            if not best_ds: # skip empty
                continue 
            
            best_file = self.best_duplicates_dic[dt]['best_file']
            policy = self.best_duplicates_dic[dt]['policy']
            
            if best_ds not in list(combined.keys() ):
                combined[best_ds] = {}
                
            if best_file not in combined[best_ds] .keys():
                combined[best_ds][best_file]= {}
                combined[best_ds][best_file]['all_indices_obs_tab'] = []
                combined[best_ds][best_file]['all_indices_header'] = []                
                combined[best_ds][best_file]['all_indices_header_in_obstab'] = []
                combined[best_ds][best_file]['data_policy_licence'] = []
                combined[best_ds][best_file]['report_id'] = []
                # TODO
                combined[best_ds][best_file]['duplicates'] = []
                combined[best_ds][best_file]['files'] = []
                combined[best_ds][best_file]['source_id'] = []
                
                
                
                combined[best_ds][best_file]['all_reading_indices_header_tab'] = []
                combined[best_ds][best_file]['all_inserting_indices_header_tab'] = []
                
                combined[best_ds][best_file]['all_reading_indices_obs_tab'] = []
                combined[best_ds][best_file]['all_inserting_indices_obs_tab'] = []

            # reading selected indices
            # indices =  self.best_duplicates_dic[dt]['indices'][5]  
            
            # storing indices for the observations table and the header table
            
            combined[best_ds][best_file]['all_inserting_indices_obs_tab'].extend( [ offset_obs + i for i in self.best_duplicates_dic[dt]['indices'][5] ] ) # must be sequential, since the best file has already been selected 
            combined[best_ds][best_file]['all_reading_indices_obs_tab'].extend( self.best_duplicates_dic[dt]['indices'][3] )
            
            combined[best_ds][best_file]['all_indices_header_in_obstab'].append( self.best_duplicates_dic[dt]['indices'][0] ) # indices of the header wrt the observations table (i.e. first occurrence of each date_time) or record_index
            index_header = list( self.unique_dates[best_ds][best_file]['datetimes'].keys() ).index(self.best_duplicates_dic[dt]['real_time'])      # self.best_duplicates_dic[dt]['real_time']
            
            combined[best_ds][best_file]['all_reading_indices_header_tab'].append(index_header) # indices of the header wrt the observations table (i.e. first occurrence of each date_time)
            combined[best_ds][best_file]['all_inserting_indices_header_tab'].append(offset_head) 

            combined[best_ds][best_file]['data_policy_licence'].extend([policy]* len(combined[best_ds][best_file]['all_inserting_indices_obs_tab']) )        
            combined[best_ds][best_file]['report_id'].extend( np.full( len(combined[best_ds][best_file]['all_inserting_indices_obs_tab']), offset_head ) )
            combined[best_ds][best_file]['files'].append(best_file)
            combined[best_ds][best_file]['source_id'].append(best_ds)
            
            combined[best_ds][best_file]['duplicates'].append(  self.best_duplicates_dic[dt]['all_duplicated_datasets'])
            
            
            # saving a small header table to simplify checks 
            records_dic['date_time'].append(dt)
            records_dic['source_id'].append(best_ds)
            records_dic['file'].append(best_file)
            records_dic['original_record_index'].append(best_file)
            records_dic['policy'].append(policy)
            records_dic['duplicates'].append(self.best_duplicates_dic[dt]['all_duplicated_datasets'])
            
            #if self.best_duplicates_dic[dt]['indices'][4] == 0:
            #    print('check ')
            # storing length of tables
            total_length_obstab = total_length_obstab +  self.best_duplicates_dic[dt]['indices'][4]  # length of selected observations
            total_length_headertab = total_length_headertab + 1

            offset_head = offset_head +   1            
            offset_obs = offset_obs +   self.best_duplicates_dic[dt]['indices'][4] 

        combined['total_length_obstab'] = total_length_obstab
        combined['total_length_headertab'] = total_length_headertab
        
        self.best_record = combined
        record_df = pd.DataFrame.from_dict(records_dic)
        self.record_df = record_df 
        
        return 0
    
    
    
    
    
    def make_obstab_var(self, var='date_time' ):
        """ Creates a vector for the given variable, reading the h5py data files at specific indices.
        Loop first over all possible datasets, then over all possible files to extract only once the indices for the particular file,
        when selected as best_record """
        
        combined_records = self.best_record  
        self.record_var = {} # storing a dictionary with data to be read for the station_configuration, header etc. e.g. latitude, longitude 
        
        # extract datasets
        datasets = [ d for d in combined_records.keys() if 'length' not in d and 'policy' not in d]
            
        # dummy observation_id
        if var == 'observation_id':
            vector = np.array(list(range(combined_records['total_length_obstab'])) )
            return vector # to eb filled later 
        
        elif var in self.empty_cdm_var:
            # special case for empty variables
            var_type = self.dic_type_attributes['observations_table'][var]['type']
            if var_type == np.int32 :
                nan = np.int32(-2147483648)
            else:
                nan = np.float32(np.nan)       
            vector = np.empty( (combined_records['total_length_obstab']) , dtype=np.dtype(nan) ) 
                
        elif var == 'source_id': 
            vector = np.full( (combined_records['total_length_obstab']) , b'dummy_source' ) # filling a vector with dataset id i.e. source_id field
        else: # other regular variable
            vector = np.empty( (combined_records['total_length_obstab'])  )                
        
        total_l = 0
        for ds in datasets:
            for file in combined_records[ds].keys():
                indices_insert = combined_records[ds][file]['all_inserting_indices_obs_tab']
                indices_read = combined_records[ds][file]['all_reading_indices_obs_tab']
                
                if var == 'source_id': 
                    np.put (vector , indices_insert, ds)
                    
                elif var  in ['sensor_id', 'adjustment_id']: # this will be replaced by the dedicated sensor script 
                    vector = np.empty( (combined_records['total_length_obstab']) )
                    
                elif var in ['advanced_assimilation_feedback' , 'advanced_uncertainty']:
                    if ds in ['era5_1' , 'era5_2']:
                        np.put (vector , indices_insert, 1)
                    else:
                        np.put (vector , indices_insert, 0)
                
                elif var in ['data_policy_licence', 'report_id']:
                    data = combined_records[ds][file][var]
                    np.put (vector , indices_insert, data)
                    
                else:
                    # accessing h5py data and slicing 
                    data = self.data[ds][file]['h5py_file']['observations_table'][var][:][indices_read]

                    #sliced_data = self.data[ds][file]['h5py_file']['era5fb'][var][indices_reading]

                    if var == 'observed_variable':  # check if the convention for CDM variable is correctly used
                        old_to_new_cdm = { 
                                        85:126, 
                                        104:139,
                                        105:140,
                                        38:138,
                                        36:137
                                    }
                        data = np.array( [ old_to_new_cdm[v]  if v in old_to_new_cdm.keys() else v for v in list(data) ] )
                    try: # replaced_indices, replaced_vector, replacing_vector)
                        #dummy = replace_global(indices_insert, vector, data)
                        np.put(vector , indices_insert, data)
                    except:
                        pass
                            
                    # updating the vector 
                    #np.put ( vector, indices_insert, data )                        
                    
        # obtain record index i.e. start of each record 
        if var == 'date_time':
            record_index = np.unique(vector, return_index=True)[1]
            self.record_index = record_index 
        if var in ['latitude' , 'longitude', 'report_id']:
            self.obs_in_header[var] = vector
            
        return vector
                 
            
    def make_era5fb_var(self, var='' ):
        """ Creates a vector for the given variable, reading the h5py data files at specific indices.
        Loop first over all possible datasets, then over all possible files to extract only once the indices for the particular file,
        when selected as best_record """
        
        combined_records = self.best_record  
        
        # extract datasets
        datasets = [ d for d in combined_records.keys() if 'length' not in d and 'policy' not in d]

        tipo = self.dic_type_attributes['era5fb'][var]['type']                   
        void = self.get_null(tipo)
        #vector = np.full( (combined_records['total_length_obstab']) , void ) # filling a vector with dataset id i.e. source_id field
        
        if var in [ 'collection_identifier@conv' , 'timeseries_index@conv', 'unique_identifier@conv', 'source@hdr', 'statid@hdr', 'expver']:
            vector = np.array( [b'NA'] * combined_records['total_length_obstab']  ) # filling a vector with dataset id i.e. source_id field
        #elif var in ['expver']:
        #    vector = np.empty( (combined_records['total_length_obstab'])  , dtype=int ) # filling a vector with dataset id i.e. source_id field
            
        else:
            vector = np.empty( (combined_records['total_length_obstab'])  ) # filling a vector with dataset id i.e. source_id field
        
        for ds in datasets:
            if 'era5' in ds:  # in this case we read the fb from the table 
                for file in combined_records[ds].keys():
                    if var not in  self.data[ds][file]['h5py_file']['era5fb'].keys(): # some files have more fb variables 
                        pass
                    else:
                        indices_reading = combined_records[ds][file]['all_reading_indices_obs_tab']
                        indices_writing  = combined_records[ds][file]['all_inserting_indices_obs_tab']
                        
                        data = self.data[ds][file]['h5py_file']['era5fb'][var][:][indices_reading]

                        try: # replaced_indices, replaced_vector, replacing_vector)
                            #dummy = replace_global(indices_writing, vector, data)
                            np.put(vector , indices_writing, data)
                            
                            #sliced_data = h5py_data[indices_reading]
                            #np.put (vector , indices_writing, sliced_data)                        
                        except:  
                            data = [ b''.join(i) for i in data ]
                            #dummy = replace_global(indices_writing, vector, data)
                            np.put(vector , indices_writing, data)
                            print(0)
                            #pass
                            ''' 
                            data = [ b''.join(i) for i in d ]
                            vector = np.empty( (combined_records['total_length_obstab'])  , dtype=np.bytes_ ) # filling a vector with dataset id i.e. source_id field                            
                            dummy = replace_global(indices_writing, indices_reading, vector, data)
                            
                            # collection_identifier@conv , timeseries_index@conv, unique_identifier@conv, stat_id ???,  does not exist in all feedback
                            # expver -> bytes
                            # source@hdr
                            # source_id
                            print("CHECK NUMBA " , var )
                            '''
            else: # here we leave empty vectors 
                pass                  
                        
        return vector

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
    
    def make_sourceconf_var(self, tab='', var='' ):
        """  """
        combined_records = self.best_record  
        
        # extract datasets
        datasets = [ d for d in combined_records.keys() if 'length' not in d and 'policy' not in d]
        
        vector = np.full( (combined_records['total_length_headertab']) , 'NA ' ).astype('|S100') # filling a vector with dataset id i.e. source_id field

        for ds in datasets:
            for file in combined_records[ds].keys():
                    sliced_data = combined_records[ds][file]['files']
                    sliced_data = [f.split('/')[-1] for f in sliced_data ] 
                    np.put (vector , combined_records[ds][file]['all_inserting_indices_header_tab'] , sliced_data)                        

        return vector
    
    def make_statconf_var(self, tab='', var='' ):
        """  """
        
        combined_records = self.best_record  

        if var in ['latitude' , 'longitude']:
            vector = 0
            return vector 
        
        values = self.data['station_configuration'][var].values[0]
        vector = np.full( (combined_records['total_length_headertab']) , values ) # filling a vector with dataset id i.e. source_id field

        return vector
    
    
        
    def merge_all_data(self):       
        """ Construct a dictionary with all the dataframes of each dataset, either reading it from saved pickles or reading it from memory """

        logging.info('***** Starting the merging process merge_all_data')

        # avoidable loops, should not matter much  # self.all_years:  
        #for this_year in [1981, 1982, 1983]:  # loop over all available years in the date_times 
        #for this_year in self.all_years:   # loop over all available years in the date_times 
        obs_offset = 0
        head_offset = 0
        
        #for this_year in self.all_years:   # loop over all available years in the date_times 
        for this_year in [1980,1981]:   # loop over all available years in the date_times 
            
            print('=== Running year ::: ' , this_year )
            self.current_year = this_year 
            selected_timestamps_this_year = [f for f in self.all_unique_dates if f.year == this_year ]
            
            # extract all indices per file, create aggregated vectors with proper lenght
            a = self.merge_all_timestamps(selected_timestamps_this_year)
            
            # extract data from harvested files for observations_table variables, except date_time
            
            all_obst_tab_vars =  self.observations_table_vars + self.empty_cdm_var
            all_obst_tab_vars.remove('date_time')
            ### WRITING OBSERVATIONS TABLE 
            print('Writing OBSERVATIONS ' ,  datetime.now()  )         
            
            for v in ['date_time']:
                date_time = self.make_obstab_var(var=v)
                self.write_merged_new(var = v, table = 'observations_table', data=date_time )
            
            for v in all_obst_tab_vars:
                ### Observations table
                #print('Extracting data ' , v , '   ' ,  datetime.now()  )      
                if v == 'observation_id':
                    data = np.array ( [  obs_offset +i for i in list(range(self.best_record['total_length_obstab']) )  ]  )
                    obs_offset = obs_offset + data[-1] # adding last entry
                    
                elif v == 'report_id':
                    dt, length = np.unique(date_time, return_counts=True )
                    report_id_obs, report_id_head = [],[]
                    for i,c,index in zip(dt, length, list(range(len(dt))) ):
                        report_id_obs.extend( [index]*c )
                        report_id_head.append(head_offset + index)
                    head_offset = head_offset + report_id_head[-1]
                    self.write_merged_new(var = v, table = 'report_id', data=report_id_obs )
                elif v == 'source_id':
                    print(0)
                    data = self.make_obstab_var(var=v)
                    
                else:
                    data = self.make_obstab_var(var=v)
                #print('Writing data ' ,  v, '   ', datetime.now()  )                                        
                self.write_merged_new(var = v, table = 'observations_table', data=data )
            print('Written OBSERVATIONS ',  datetime.now()  )       
            
            
            ### WRITING FEEDBACK TABLE 
            print('Writing FEEDBACK ' ,  datetime.now()   )            
            for d in self.era5fb_columns:
                if d in ['source_id', 'stat_id']:
                    continue 
                #print('Extracting data ' , d , '   ' ,  datetime.now()  )                                                        
                data = self.make_era5fb_var(var=d)
                #print('Writing data ' ,  d, '   ', datetime.now()  )                                                        
                self.write_merged_new(var = d, table = 'era5fb', data=data )                
            print('Written FEEDBACK ' ,  datetime.now()   )            
                
                
            ### WRITING HEADER TABLE 
            print('Writing HEADER ')
            for d in self.header_columns:
                if d == 'report_id':
                    data = np.array(report_id_head)
                    
                else:
                    data = self.make_header_var(tab='header_table', var=d)
                self.write_merged_new(var = d, table = 'header_table', data=data )             
            print('Written HEADER ')
            
            
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
                         
            ### WRITING SOURCE CONFIGURATION
            source_files =  self.make_sourceconf_var(tab='source_configuration', var='source_file' )                        
            self.write_merged_new(var = 'source_file', table = 'source_configuration', data=source_files )              
            # OLD,, slow
            #print('Writing source_configuration ' )                        
            #source_conf=xr.Dataset()
            #source_files =  self.make_sourceconf_var(tab='source_configuration', var='source_file' )            
            #source_files = np.array(source_files).astype(dtype='|S70')
            #source_conf['source_file'] = ( {'source_file' : source_files.shape } , source_files )
            #self.write_merged_new(var='source_file', table = 'source_configuration', data= source_conf )

            ### WRITING CDM METADATA TABLES
            self.write_merged_new(table = 'cdm_tables', data= self.data , var='')
            
            ### Adding sensor_id to observations_table 
            if self.add_sensor:
                print('*** Adding sensor *** ')            
                add_sensor = wrapper(out_dir = self.out_dir , station_id = self.station.split('-')[-1] , file = self.out_name , copy = self.copy )
                print('*** Added sensor *** ')
            else:
                pass
            
            for a in [1,2]:
                print(a)
                ### era5fb table 
                
                # write the obstab variable immediately 
                ### write the data immediately 

                '''
                ### Station configuration table 
                if '20999' in self.station:
                    primary, name = self.station, np.bytes_('UNKNOWN')
                else:
                    primary, name = self.data['station_configuration']['primary_id'].values[0] , self.data['station_configuration']['station_name'].values[0] 
                
                combined_head_tab['primary_station_id'] = np.array( [primary] )
                combined_head_tab['station_name']         = np.array( [name] )
                   
                """ New merged recordindex and recordtimestamps indices """
                combined_indices.append(len(combined_obs_tab['date_time']))                 
                combined_date_time.append(dt)
    
                del cleaned_df_container 
            
                #print(blue + 'Memory used after deleting the cleaned_df_container: ', process.memory_info().rss/1000000000 , cend)
        
            ### STATION CONFIGURATION 
            
            if '20999' not in self.station:
                stat_conf = pd.DataFrame(np.repeat(self.data['station_configuration'].values , len(all_combined_head) , axis=0), columns=self.data['station_configuration'].columns)                  
        
            else:
                #stat_conf = self.stat_conf_CUON 
                for v in self.data['station_configuration'].columns:
                        a = self.get_null( tipo = self.dic_type_attributes['station_configuration'][v]['type'])
                        print(v,  '   ' ,  self.dic_type_attributes['station_configuration'][v]['type'] , '   ' , a )
                        print(0)
                        
                stat_conf = pd.DataFrame(np.repeat(self.data['station_configuration'].values , len(all_combined_head) , axis=0), columns=self.data['station_configuration'].columns)  
                stat_conf['primary_id'] = np.bytes_(self.station)
                stat_conf['primary_id'] = stat_conf['primary_id'].astype('|S15')
                stat_conf['platform_type'] = 2
                stat_conf['platform_type'] = 2
                metadata_contact = 'L.Haimberger'
                metadata_contact_role = '0'
        
            # check other variables !!! 
            stat_conf.latitude = [i['latitude'][0] for i in all_combined_head ]
            stat_conf.longitude = [i['longitude'][0] for i in all_combined_head ]
            stat_conf.record_number = np.array ( range(len(all_combined_head) ) )
            
            stat_conf = clean_station_configuration(stat_conf)        
            
            self.data['station_configuration'] = stat_conf 
                
            for k in self.data['station_configuration'].columns: # try writing only one entry

                a =np.array( self.data['station_configuration'][k])
                self.write_merged(content = 'station_configuration', table= {k:a})
                logging.info('*** Written station_configuration %s             ', k)
                        
            '''




    def write_merged_new(self, var='', table = '', data=''):
        """ Module to write the output file as netCDF.
        NB a table is made of only one variable"""

        #logging.debug('Writing the variable ' , var, ' for the table: ', table )
        out_dir = self.out_dir + '/' + str(self.current_year)
        if not os.path.isdir(out_dir):
            Path(out_dir).mkdir(parents=True, exist_ok=True)                  
        out_name = out_dir + '/' + self.station + '_CEUAS_merged_v2.nc'  
        self.out_name = out_name 
        
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

                var_type = self.dic_type_attributes[table][var]['type']
    
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
            for k in data['cdm_tables'].keys():
                table = data['cdm_tables'][k]
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
            

        """
        elif table == 'station_configuration':
            var_type = self.dic_type_attributes[table][k]['type']

            ''' trying to convert the variable types to the correct types stored as attribute, read from the numpy dic file '''
            if type(table[k][0]) != var_type:
                try:
                    table[k] = table[k].astype( var_type ) 
                    #print('Done station_conf' , k )
                except:
                    if k == 'secondary_id':
                        table[k] = table[k].astype( bytes ) 
                    else:
                        table[k] = table[k].astype( np.float32 )                             
                        print ('FAILED converting column ' , k, ' type ', type(table[k][0]) , ' to type ', var_type , '  so used np.float64')

            dic = {k:table[k]}  

            write_dict_h5(out_name, dic , table, self.encodings[table], var_selection=[], mode='a', attrs = attrs_dic  )
        """
            
    def merge(self ,  station = '' , datasets = '' , mode = 'test'):        

        """ Call to the merge_all_data() and write_merged_file() methods """
        if mode == "test":      
            a = self.initialize_data( station = station, datasets = datasets ) # reading the input files 
            dummy = self.merge_all_data()            
            #logging.info('*** Finished merging, now writing the output netCDF file ***' )         
            #a = self.write_merged_file()
            #logging.info('*** Done writing the output ! ***')
            self.write_merging_summary()
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
        out = open(self.out_dir + '/merging_summary.txt', 'a+')
        out.write(self.out_name + '\n')
        for d in self.datasets.keys():
            out.write(d + '\t')
            for f in self.datasets[d]:
                out.write(f + '\t')
            out.write('\n')





'''
def create_stat_dic(stat_id, data_directories, kind):
    """ Looks in the database for files matching with the given station id.
    Returns a dictionary with the files per dataset,
    the total size of the files t be processed,
    the station_id to be used for the output file. """
    station_dic = {}
    total_dim = []
    found_20001 = False

    for d,i in data_directories.items():
        #if kind == 'mobile':
        #    if d not in ['era5_1_mobile' , 'era5_2_mobile']:
        #       continue
        files = os.listdir(i)
        for f in files: # normal check for any file
                Id = f.split('_'+d)[0]
                if Id == stat_id:
                    if d not in station_dic.keys():
                        station_dic[d] = []                            
                    station_dic[d].append(i + '/' + f)
    
                    total_dim. append( os.path.getsize (i + '/' + f) )
                if '-20001-' in f:  # here, look for possible alternatives identified as e.g. 0-20001-0-10393 to be merged together
                    stat = stat_id.split('-')[-1]
                    stat_id_2 = '0-20001-0-' + stat 
                    if Id == stat_id_2:
                        if d not in station_dic.keys():
                            station_dic[d] = []                            
                        station_dic[d].append(i + '/' + f)
                        found_20001 = True
    
        size = sum(total_dim)        
        if found_20001:
            stat_id = stat_id.replace('-20000-','-20001-')
    return station_dic, size, stat_id
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

                    
                    station_dic['files'].append(i+'/' + ','.join(harvested_files)  ) 
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
    stations = ['0-20000-0-71879']  # 0-20000-0-71879 , 0-20000-0-82900
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
