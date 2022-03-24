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
#from numba import njit
import psutil
import copy
from numba import njit
import code

sys.path.append('../harvest/code')

sys.path.append('../postprocess/add_sensor')

from add_sensor_to_merged import Sensor, datetime_toseconds, wrapper, MergedFile
from harvest_convert_to_netCDF import write_dict_h5 

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
        self.observation_ids_merged  = {  'igra2':b'3' , 
                                                               'ncar':b'4', 
                                                                   'bufr':b'5',  
                                                                   'era5_1':b'1' , 
                                                                   'era5_2':b'2', 
                                                                   'era5_1759' :b'6' , 
                                                                   'era5_1761':b'7' ,  
                                                                   'era5_3188' :b'8' }  # values used to convert original record_id to the merged record_id, see method merge_all_data 

        logging.info('*** Initialising the Merging procedure ***' )   
        #self.era5b_columns = []  # stores the columns of the era5fb 
        self.standard_cdm = [ 'crs' , 'observed_variable', 'units' , 'z_coordinate_type' , 'station_type', 'station_configuration_codes'] 
        self.slice_size = 3000 # size of the chunk of the subset of data 
        self.index_offset = 0 # will be replaced when running 
        self.hour_time_delta = 60 * 60 * 2 # decide up to which time shift records are considered identical  [unit = seconds]

        self.only_std_plevels = False  # set to True to store only standard pressure level data 
        self.std_plevs    = [1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 85000, 92500, 100000]

        self.add_sensor = True
        self.copy = True # make a copy of the merged file before adding the sensor. Argument for the add_sensor wrapper function  
        
    def load_obstab_feedback_sliced(self, dataset='' , file ='' , datetime='' ):
        """ Specific function to load the observations_table and feedback_table for the selected dataset. 
             If index==True, it will check the self.slice_size and look for the correpsonding recordtimestamp and recordinidex in the era5_1 data.
             It will load only a numer or records equal to self.slice_size + 1 , not to use massive memory.
        """
        k = dataset 
        F = file 
        dt = datetime

        if dt != self.unique_dates[k][F]['up_to_dt_slice']:
            print("Error! the dit does not correspond to the dt I calculated in the previous loading! ")
            return 0

        logging.debug(" === (Re)Load data for %s      file %s      counter %s" , dataset, file, data[k][F]["counter"])
        print(blue + 'Memory used before reading data: ', process.memory_info().rss/1000000000 , cend)

        slice_size = self.slice_size

        file     = data[k][F]['h5py_file']
        rts, ri = data[k][F]["recordtimestamp"][:] , data[k][F]["recordindex"][:]

        index_min  = self.unique_dates[k][F]['indices'][dt]['low']  # here no offset since I am reading the original data 
        ind             = np.where(rts==dt)[0][0]                               # index of specific dt , I need the extremes indices of the next date_time after slicing 

        try:            
            up_to_dt_slice = rts[ind + slice_size  ]  # 
            index_max = self.unique_dates[k][F]['indices'][up_to_dt_slice]['low']  # maximum index in the array of date_time to slice on
            update_index = True
        except:
            """ If the dt is too large, I take the whole array """
            index_max = 1000000000000000
            update_index = False 


        ####################
        # OBSERVATIONS TABLE
        ####################                  
        logging.debug ('*** Loading observations_table' )
        obs_tab = file['observations_table']                                                                                                      

        #print('CHECKING THE INDICES:::: ' , k , ' index_min ', index_min , ' index_max ', index_max )
        obs_dic= {} 
        for ov in self.observations_table_vars:
            v = copy.deepcopy( obs_tab[ov][index_min:index_max ] )
            obs_dic[ov] = v 
        data[k][F]['observations_table']= obs_dic 

        ###########
        # ERA5FB
        ###########
        if k == 'era5_1' or k == 'era5_2':
            logging.debug('*** Loading era5fb ' )
            era5fb_tab = file['era5fb']
            fb_dic = {} 
            for ov in self.era5fb_columns:
                try:
                    v = copy.deepcopy( era5fb_tab[ov][index_min:index_max ] )
                    fb_dic[ov] = v 
                except:
                    continue
                    #print("CANNOT FIND  ", ov ) 

            data[k][F]['era5fb_tab']= fb_dic

        print(blue + 'Memory used after reading data: ', process.memory_info().rss/1000000000 , cend)

        """ Updating the indices """ 
        self.unique_dates[k][F]['index_offset']          = copy.deepcopy( self.unique_dates[k][F]['index_offset_next'] )    

        if update_index:    
            self.unique_dates[k][F]['index_offset_next'] = index_max 
            self.unique_dates[k][F]['up_to_dt_slice']       = up_to_dt_slice

        return 0 


    def initialize_data(self , station = '', datasets = {} ):
        """ Initialize dataset; store relevant data as attributes.
                   Args ::     dic{}  datasets (dictionary where keys are the dataset names e.g. bufr, igra2 etc. , and the value is the path to the corresponding netCDF file 
                                   e.g. tateno = { 'ncar'    : 'example_stations/ncar/chuadb_trhc_47646.txt.nc'   ,
                                                           'igra2'   : 'example_stations/igra2/chJAM00047646-data.txt.nc'  } 
        """       
        
        self.datasets          = datasets
        self.datasets_keys = datasets.keys()
        self.station             = station

        self.observations_table_vars = ['date_time', 'z_coordinate' , 'z_coordinate_type', 'observed_variable' , 
                                                           'observation_value' , 'report_id' , 'observation_id' , 'latitude' , 'longitude',
                                                           'units', 'source_id', 'data_policy_licence' ]

        """ Loading the econding of the tables created from the harvester script and to be applied again """
        self.encodings = np.load('groups_encodings.npy' , allow_pickle = True ).item()
        self.encodings['era5fb'] = np.load('era5fb_encodings_all.npy' , allow_pickle = True ).item()            
        self.dic_type_attributes = np.load('dic_type_attributes.npy',allow_pickle= True).item()

        self.era5fb_columns = self.dic_type_attributes['era5fb'].keys()

        self.obstab_nans_filled = False      

        data['cdm_tables'] = {}              

        """ Loop over all the datasets                                                                                                                                     
                k: name of the dataset                                                                                                                    
                v: list of file paths,    eg 'era5_1':[filepath_1, filepath_2 ] """    
        
        # preliminary loop, quick check for consistency of station_configuration

        stat_conf = pd.read_csv("CUON_station_configuration.csv", sep = '\t')
        s = stat_conf.loc[ stat_conf.primary_id == self.station ]
        if s.empty:
            a = open('logs/failed_stat_conf.txt' , 'a+')
            a.write(self.station + '\n')

        data['station_configuration'] = s

        for k,v in self.datasets.items() :
            data[k] = {}
            for F in v:

                logging.info(' Dataset ::: *** %s %s ' , k , F   )                  

                data[k][F] = {}

                h5py_file = h5py.File(F, 'r')
                data[k][F]['h5py_file'] = h5py_file 

                a = h5py_file['recordtimestamp']

                data[k][F]['recordtimestamp'] = a
                data[k][F]['recordindex']         = h5py_file['recordindex']
                data[k][F]['dateindex']            = h5py_file['dateindex']
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
                    if ('string' in var or 'hdrlen' in var): continue
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
        self.make_all_datetime()
        
        return True 


    def delete_ds(self, dt):
        """ Delete the dataset from the memory once the maximum date of data availability has been reached;
             load the era5_1 once the date_time is in the valid range """

        for k in self.datasets_keys:
            for F in self.datasets[k]:
                if F not in data[k].keys():
                    continue  
                max_date = data[k][F]['max_date'] 
                """ Deleting unecessary ds """
                if dt > max_date : # check max date and check if data is still loaded
                    print(blue + 'Memory used before deleting : ' , process.memory_info().rss/1000000000 , cend)                               
                    del data[k][F] 
                    print("*** Erasing dataset: " , k , ' ' , F ) 
                    print(blue + 'Memory used after deleting : '    , process.memory_info().rss/1000000000 , cend)                               

                else:
                    continue


    def make_all_datetime(self):
        """ Building the global set of date_times and indices from the various datasets. 
              The datetimeindex is read from the original netCDF file. 
              Will compare the unique date_time of each dataset and extract the global unique date_times
              """

        logging.info('\n *** Running make_all_datetime ' )

        all_uniques     = []  # storing a list with all the unique date_times            
        which_k_in_dt = {}  # list of avilable dataset for each unique date_time, so that when looping over the distinct date_times, only the proper dataset will be read and compared 

        """ Loop over all the datasets 
                k: name of the dataset
                v: list of file paths,    eg 'era5_1':[filepath_1, filepath_2 ]"""

        for k,v in self.datasets.items() :
            self.unique_dates[k] = {}
            for F in v: 
                self.unique_dates[k][F] = {}

                self.unique_dates[k][F]['indices'] = {}                             
                self.unique_dates[k][F]['index_offset_next'] = 0 # to be replaced later when slicing 
                self.unique_dates[k][F]['index_offset'] = 0 # to be replaced later when slicing 

                unique_dt = list(data[k][F]['recordtimestamp'])

                indices   = list(data[k][F]['recordindex'])
                all_uniques += unique_dt   # adding to the total unique date_times 

                """ Loop over all the date_times of each dataset """
                for dt, index_low, count in zip (unique_dt,  indices, range(len(unique_dt))  ):

                    if dt not in which_k_in_dt.keys():
                        which_k_in_dt[dt] = {}
                    if k not in which_k_in_dt[dt].keys():
                        which_k_in_dt[dt][k] = []                             
                    if F not in which_k_in_dt[dt][k]:
                        which_k_in_dt[dt][k].append(F)
                    # at this point I have e.g.  which_k_in_dt= {1990-01-01-12-00: {era5_1:[file1,file2] , ncar:[file3] } }

                    self.unique_dates[k][F]['indices'][dt] = {}
                    self.unique_dates[k][F]['indices'][dt]['low'] = index_low                       
                    try:
                        index_up =  indices[ count + 1 ]  # works until the last available recordindex
                    except:                            
                        index_up = max(indices)+1000000 # dummy large number    

                    self.unique_dates[k][F]['indices'][dt]['up'] = index_up
                    self.unique_dates[k][F]['up_to_dt_slice']     = data[k][F]['min_date'] 


        self.dataset_per_dt = which_k_in_dt             
        self.merged_unique_dates = np.unique(np.array(all_uniques) )  # storing the set of *ALL* distinct dt values of all datasets and all files            
        logging.debug('*** make_all_datetime finished ')         
    print(blue + 'Memory used after makind all date_times : ', process.memory_info().rss/1000000000 , cend)

    def get_header_table(self , dt, ds = '' , File = ''):
        """ Extracting the header_table """

        index= np.searchsorted(data[ds][File]['recordtimestamp'], dt)
        hd = {}
        for v in data[ds][File]['header_table'].keys():
            hd[v]    = np.array( [data[ds][File]['header_table'][v][index] ])

        return hd

    def make_obstab_era5fb_dic(self, dataset = '' , date_time = '', File = ''):            
        """ Create obs_tab and feedback tables """
        index_offset = self.unique_dates[dataset][File]['index_offset']

        # Removing the index_offset, which is defined only if any slicing was done 
        index      = self.unique_dates[dataset][File]['indices'][date_time]['low'] - index_offset
        index_up = self.unique_dates[dataset][File]['indices'][date_time]['up'] - index_offset

        """ Extract the current chunks from the dict """
        obs_dic = {}                                    
        for v in self.observations_table_vars:
            obs_dic[v]    = data[dataset][File]['observations_table'][v][index:index_up]
            #print('v is : ',  v )

        """ Loop over the obs_tab to find duplicates.
                I fill a dictionary for each distinct pressure level, and I put inside
                the observed_variable number.
                If the list lready contains the combination pressure level - observed variable,
                then the record is skipped """

        indices = [] # these are the only non-duplicates to be kept

        already_selected = { }
        for p in list(set(obs_dic['z_coordinate'])):
            already_selected[p] = [] # create a list for every distinct plevel                  

        #print('starting the loop: ' , date_time, ' '  , dataset, ' ', index, ' ' , index_up)
        for p,var,val,ind in zip ( obs_dic['z_coordinate'] , obs_dic['observed_variable'], obs_dic['observation_value'] ,range(len(obs_dic['z_coordinate'])) ):
            #print(p,var,val,ind)
            #if date_time > 2354300000:
            #      print('looping :::', var, ' ' , val, ' ' , ind , ' ', dataset, ' ' , index_up, ' ' , index, ' ', File)

            if self.only_std_plevels:
                if p not in self.std_plevs:
                    continue 

            # create a list for each available plevel 
            if p not in already_selected.keys():
                already_selected[p] = []
            """ to test
                  a = {'z': obs_dic['z_coordinate'] , 'var': obs_dic['observed_variable'] , 'val':  obs_dic['observation_value'], 'time': obs_dic['date_time']  }
                  b = pd.DataFrame.from_dict(a)
                  """
            if np.isfinite(val): # check what it does here exactly ??? 
                if var not in already_selected[p]:
                    already_selected[p].append(var)
                    indices.append(ind) # record to be kept
                else:
                    #print('There is a duplicate *** ')
                    pass
            else: # skipping nans
                pass

        #print('done with the loop')
        red_obs_dic = {} # dictionary for the reduced (removed duplicates) obs_tab
        for v in self.observations_table_vars:
            red_obs_dic[v] = obs_dic[v][indices]

        ''' Simply returns the proper format for ''null' value '''
        def get_null( tipo = ''):
            if tipo == np.int32 :
                void = -2147483648
            elif tipo == np.float32 :
                void = np.nan
            elif tipo == np.bytes_ :
                void = b'nan'
            return void

        ''' Filling the feedback table. Only feedback for era5_1 and era5_2 are currently available. 
                Reads the total number of possible columns from the dic_type_attributes dictionary.
                Era5_1 and era5_2 fb have different columns.
                If data for a variable is not available, it fills with the appropriate null value '''

        #print('making the era5fb  ', date_time, ' ' , dataset)
        red_era5fb_dic = {}
        for v in self.era5fb_columns:
            tipo = self.dic_type_attributes['era5fb'][v]['type']                   
            if dataset == 'era5_1' or dataset == 'era5_2':
                if v in data[dataset][File]['era5fb_tab'].keys():                            
                    red_era5fb_dic[v] = data[dataset][File]['era5fb_tab'][v][index:index_up][indices]
                else:
                    void = get_null(tipo = tipo)
                    red_era5fb_dic[v]= np.full(len(indices), void)                              
            else:       # no feedback for non era%-1 or era5_2 datasets 
                void = get_null(tipo = tipo)
                red_era5fb_dic[v]= np.full(len(indices), void)

        #print('done making_obstab_era5fb')
        """
            try:
                  if  len(red_obs_dic['date_time']) > 2:
                        print('yes')
                  else:
                        print('check')                        
            except:
                  print('check')
            """      
        return red_obs_dic , red_era5fb_dic    

    def is_same_record(self, time_shift = '' , dt = ''):

        if (dt - self.processed_dt[-1] ) < time_shift :
            is_same_record = True
        else:
            is_same_record = False

        return is_same_record 


    def merge_all_data(self):       
        """ Construct a dictionary with all the dataframes of each dataset, either reading it from saved pickles or reading it from memory """

        logging.info('***** Starting the merging process merge_all_data')

        """ All possible unique_dates to loop on """
        date_times = self.merged_unique_dates
        date_times.sort()
        date_times = np.array(date_times) 

        """ List storing the indices of the date_index of the merged dataset """
        all_combined_obs ,  all_combined_head, all_combined_era5fb , combined_indices , combined_date_time,  = [] , [] , [] , [] , []
        best_ds_list = [] 
        source_files = []

        """ The items contained in the lists in the list below can be removed from the list when the record that was previously stored is removed. """
        all_list = [all_combined_obs ,  all_combined_head, all_combined_era5fb , combined_indices , combined_date_time,  best_ds_list, source_files ]  
        # list holder of all the above lists

        #all_list_name = ['all_combined_obs' ,  'all_combined_head', 'all_combined_era5fb' , 'combined_indices' , 'combined_date_time'  , 'best_ds_list', 'source_files' ] 
        #removed_record, kept_record = [], []

        """ Dictionary that will contain the merged file. """            
        # rand = datetime.strptime('1981-01-03 12:00:00', '%Y-%m-%d %H:%M:%S')  
        #dt_bestds_dic = {} # store the selected best dataset for each dt     
        #date_times=date_times[0:30000]
        tot = len(date_times)
        tt=time.time()
        print('*** Merging ' , tot, '  records ***')

        #early_datasets = True

        self.processed_dt = [] 

        for dt, c in zip(date_times, range(tot) ): # loop over all the possible date_times 

            if (c+1)%1000==0:
                print('Analize : ', str(c+1) , '/',  str(tot)  , ' ', str(dt/(365.25*3600*24)) , ' ',
                              now(time.time()),'{:5.3f}'.format(time.time()-tt ))

            delete = self.delete_ds(dt) # check if there is a dataset to delete 

            """ Finding if this record is the same as the previous one analyzed, according to the given time_shift """
            if c == 0:
                is_same_record = False
            else:
                is_same_record = self.is_same_record( time_shift = self.hour_time_delta , dt = dt)

            """ Updating list of processed datetimes """
            self.processed_dt.append(dt) # cannot put it before the check_timeshift or it will check itself 


            cleaned_df_container = {}          
            all_len = [] # will hold the length of all the obs_tabs 

            for k in self.dataset_per_dt[dt].keys() :  # checking the list of available datasets  
                ''' {'era5_2': ['example_stations/0-20000-0-82930_era5_2_harvested_era5.conv._1:82930.gz.nc', 
                                        'example_stations/0-20000-0-82930_era5_2_harvested_era5.conv._82930.gz.nc']}
                        '''                        
                for F in self.dataset_per_dt[dt][k]: # checking the list of available files for the dataset

                    if data[k][F]["counter"] %self.slice_size==0 or data[k][F]["counter"]  == 0:  # loading the data only at specific slices 
                        load = self.load_obstab_feedback_sliced(datetime=dt, dataset=k, file = F)

                    data[k][F]["counter"] =  data[k][F]["counter"]  + 1 

                    obs_tab, era5fb_tab = self.make_obstab_era5fb_dic(dataset = k , date_time = dt, File = F )

                    if len(obs_tab['date_time'][:])==0: # go to next file if obs_tab is empty 
                        continue                              

                    all_len.append( len(obs_tab['date_time'][:] ) )

                    if k not in cleaned_df_container.keys():
                        cleaned_df_container[k] = {}

                    cleaned_df_container[k][F] = {}
                    cleaned_df_container[k][F]['obs_tab']      = obs_tab         # cleaned dataframe 
                    cleaned_df_container[k][F]['era5fb_tab'] = era5fb_tab     # cleaned dataframe  

            """ Merging the different records found in the different sources """
            
            if len(all_len)>0: # skipping empty container dictionary. At this point I certainly have one valid record 
                best_ds, combined_obs_tab, combined_era5fb_tab, combined_head_tab, selected_file, best_file = self.combine_record(dt, container = cleaned_df_container)

                #if is_same_record and  dt/(3600*24*365)  + 1900 > 1961:
                if is_same_record:

                    remove_previous = False 
                    
                    #print('Found a time-shifted record, current' + best_ds + '- previous ' + best_ds_list[-1] + ' = ' , float( (combined_obs_tab['date_time'][0] - temporary_previous['date_time'][0])/3600)  ) # decide what to keep in case of same record

                    """ Check if the licence needs to be updated (at least one is free so we set it as free for all) """
                    current_licence = combined_obs_tab['data_policy_licence'][0] 
                    previous_licence = all_combined_obs[-1] ['data_policy_licence'][0] 
                    
                    if current_licence != previous_licence: # if the policy are not the same it means that at least one is free so set to free 
                        if current_licence != 0:
                            combined_obs_tab['data_policy_licence'] = np.full ( len (combined_obs_tab['date_time']) , 0 , dtype = int)
                        if previous_licence !=0:
                            all_combined_obs[-1]['data_policy_licence'] = np.full ( len (all_combined_obs[-1]['date_time']) , 0 , dtype = int)
                            
                    temporary_previous = all_combined_obs[-1] # keep the temporary previous record 
                    
                    if best_ds in ['era5_1','era5_2']:  # best_ds from era5
                        if  best_ds_list[-1] not in ['era5_1','era5_2']: # remove previous non era5_1 or era5_2 record 
                            remove_previous = True
                            #removed_record.append(temporary_previous)
                            #kept_record.append(combined_obs_tab)                                                  

                        elif best_ds_list[-1] in ['era5_1','era5_2']:
                            if len(combined_obs_tab['date_time']) <= len(temporary_previous['date_time'] ):
                                #kept_record.append(temporary_previous)  
                                #removed_record.append(combined_obs_tab)
                                all_combined_head[-1]['duplicates'] = np.array( [ all_combined_head[-1]['duplicates'][0] + b',' + combined_head_tab['duplicates'][0] ], dtype='|S70') 
                                continue  # nothing to do, will keep the previous records -> go to next dt 

                            else: # case where both the current and previous are from era5_1 and era5_2, but the previous has fewer data 
                                remove_previous = True
                                #removed_record.append(temporary_previous)
                                #kept_record.append(combined_obs_tab)                                                        

                    else:  # best_ds not from era5
                        if best_ds_list[-1] in ['era5_1','era5_2']:
                            #print('This best ds is ' , best_ds , '  but I will keep ' ,  best_ds_list[-1] )
                            #kept_record.append(temporary_previous)  
                            #removed_record.append(combined_obs_tab)   
                            all_combined_head[-1]['duplicates'] = np.array( [ all_combined_head[-1]['duplicates'][0] + b',' + combined_head_tab['duplicates'][0] ], dtype='|S70') 
                            continue 

                        else:
                            if len(combined_obs_tab['date_time']) < len(temporary_previous['date_time'] ):
                                #kept_record.append(temporary_previous)  
                                #removed_record.append(combined_obs_tab)
                                all_combined_head[-1]['duplicates'] = np.array( [ all_combined_head[-1]['duplicates'][0] + b',' + combined_head_tab['duplicates'][0] ], dtype='|S70') 
                                continue  # nothing to do, will keep the previous records -> go to next dt 

                            elif len(combined_obs_tab['date_time']) > len(temporary_previous['date_time'] ): # remove previous, keep current 
                                remove_previous = True                                
                                #kept_record.append(combined_obs_tab)  
                                #removed_record.append(temporary_previous)

                            elif len(combined_obs_tab['date_time']) == len(temporary_previous['date_time'] ): # prefer igra2, otherwise
                                if best_ds == 'igra2':
                                    remove_previous = True
                                    #removed_record.append(temporary_previous)
                                    #kept_record.append(combined_obs_tab)  

                                else: # case where data source is not important, I keep the previous and do nothing 
                                    #kept_record.append(temporary_previous)  
                                    #removed_record.append(combined_obs_tab)
                                    all_combined_head[-1]['duplicates'] = np.array( [ all_combined_head[-1]['duplicates'][0] + b',' + combined_head_tab['duplicates'][0] ], dtype='|S70') 
                                    continue    
                                
                    if remove_previous:
                        # remove from lists the previous data
                        # copy the duplicated from before into new header 
                        combined_head_tab['duplicates'] = np.array( [ all_combined_head[-1]['duplicates'][0] + b',' + combined_head_tab['duplicates'][0] ], dtype='|S70') 
                        for lista in all_list:
                            lista.pop()   
                            
                else: # not the same record, nothing special to do, keep both previous and current 
                    pass                     
            else:
                print(' Found an empty record')
                continue


            """ Fill the best_ds list """
            best_ds_list.append(best_ds)

            """ Storing the selected file for the source_configuration """
            source_files.append(selected_file)

            """ Storing the combined era5fb, header and observations tables"""
            all_combined_era5fb.append(combined_era5fb_tab)
            all_combined_obs   .append(combined_obs_tab)

            primary, name = self.data['station_configuration']['primary_id'].values[0] , self.data['station_configuration']['station_name'].values[0]

            combined_head_tab['primary_station_id'] = np.array( [primary] )
            combined_head_tab['station_name']         = np.array( [name] )

            all_combined_head  .append(combined_head_tab)

            """ Dictionary to fill the best_ds for duplicates """
            #dt_bestds_dic[dt] = {}
            #dt_bestds_dic[dt]['best_ds'] = best_ds
            #dt_bestds_dic[dt]['len'] = len(combined_obs_tab['date_time'])

            """ New merged recordindex and recordtimestamps indices """
            combined_indices.append(len(combined_obs_tab['date_time']))                 
            combined_date_time.append(dt)

            del cleaned_df_container 


        """ Removing remaining loaded df """
        for k in self.datasets_keys:
            for F in self.datasets[k]:
                try:
                    del data[k][F]['era5fb_tab']
                    print('=== removed era5fb ' , k , F )
                except:
                    pass
                try:
                    del data[k][F]['observations_table']
                    print('=== removed obstab ' , k , F )   
                except:
                    pass


        """ Saving a numpy dictionary """
        print(" === Saving the numpy dictionary of removed and kept records +++ ")
        #dic_records = { 'kept' : kept_record , 'removed': removed_record }
        #np.save(self.station + '_time_shift_removed_kept.npy',dic_records )


        """ Storing the merged date_time values and indices """
        di=xr.Dataset()
        combined_date_time = np.array(combined_date_time)
        di['recordtimestamp'] = ( {'recordtimestamp' : combined_date_time.shape } , combined_date_time )
        di['recordtimestamp'].attrs['units']='seconds since 1900-01-01 00:00:00'

        """ Creating the merged indices mi """
        mi = [] 
        mi.append(0)
        for i in range(len(combined_indices)):
            mi.append( combined_indices[i] + mi[-1] )
        mi.pop()
        pop = np.array(mi) # removing last unecessary index  
        di['recordindex']          = ( {'recordindex' : pop.shape } , pop )


        """ Creating the combined data """
        logging.debug('*** Concatenating the observations_table ' )      
        combined_obs = {}
        ####  Writing combined observations_table dic
        logging.info(' ***** Writing the observations_table to the netCDF output ***** ' )             
        for k in all_combined_obs[0].keys():        
            a = np.concatenate([all_combined_obs[i][k][:] for i in range(len(all_combined_obs))])
            if k == 'date_time':
                combined_obs[k]= a  
                self.tot_records = len(combined_obs[k])
            if k =='observed_variable':
                # check if new numberidng convention is really satisfied 
                old_to_new_cdm = { 
                                   85:126, 
                                   104:139,
                                   105:140,
                                   38:138,
                                   36:137
                                   }
                
                a = np.array( [  old_to_new_cdm[v]  if v in old_to_new_cdm.keys() else v for v in list(a) ] )
                
            self.write_merged(content = 'observations_table', table= {k:a})
            #logging.info('*** Written observations table %s: ', k)


        #self.tot_records = len(combined_obs['date_time'])
        del all_combined_obs
        print(blue + 'Memory used after deleting all_combined_obs dic: ', process.memory_info().rss/1000000000 , cend )

        dateindex = combined_obs['date_time']//86400                                                                                                                                  
        date_times, indices, counts = np.unique(dateindex, return_counts = True, return_index= True)                                                                                  
        di['dateindex'] = ( {'dateindex' : indices.shape } , indices )  # considers the day only                                                                                      
        del combined_obs

        combined_era5fb = {}
        ####  Writing combined era5fb_table dic                                                                                                                      
        for k in all_combined_era5fb[0].keys():
            try:
                #combined_era5fb[k]=np.concatenate([all_combined_era5fb[i][k][:] for i in range(len(all_combined_era5fb))])
                #self.write_merged(content = 'era5fb', table= {k:combined_era5fb[k]})
                """ try replacing , remove combined_era5fb = {} """
                a = np.concatenate([all_combined_era5fb[i][k][:] for i in range(len(all_combined_era5fb))])
                self.write_merged(content = 'era5fb', table= {k:a})
                logging.debug('*** Written era5fb %s:  ', k)
            except:
                #print("Failed feedback variable " , k)
                pass
        del all_combined_era5fb
        print(blue + 'Memory used after deleting era5fb_tab dic: ', process.memory_info().rss/1000000000 , cend)


        ####  Writing combined header_table dic                                                                                   
        for k in all_combined_head[0].keys():
            #print('head variable is', k )
            if ( k == 'comments' or k == 'history'):
                continue
            try:
                tab=np.concatenate([all_combined_head[i][k][:] for i in range(len(all_combined_head))])
                self.write_merged(content = 'header_table', table= {k: tab})  # { key: np.array([])}
                logging.info('*** Written header table %s: ', k)
            except:
                #print('Failed variable in header table', k )
                pass
        del all_combined_head
        print(blue + 'Memory used after deleting all_merged head_tab dic: ', process.memory_info().rss/1000000000 , cend)

        self.write_merged(content = 'recordindex', table = di)                      
        self.write_merged(content = 'cdm_tables', table= '')


        source_conf=xr.Dataset()
        source_files = np.array(source_files).astype(dtype='|S70')
        source_conf['source_file'] = ( {'source_file' : source_files.shape } , source_files )
        self.write_merged(content = 'source_configuration', table= source_conf )


        ### STATION CONFIGURATION                                                                                                                                                                                                                                                         

        for k in self.data['station_configuration'].columns: # try writing only one entry                                                                                                                                                                                                 
            if "Unnamed" in k:
                continue
            a =np.array( self.data['station_configuration'][k])
            self.write_merged(content = 'station_configuration', table= {k:a})
            logging.info('*** Written station_configuration %s:  ', k)

        
        """ Adding sensor_id to observations_table """
        if self.add_sensor:
            add_sensor = wrapper(out_dir = self.out_dir , station_id = self.station.split('-')[-1] , file = self.out_name , copy = self.copy )
            print('*** Added sensor *** ')
        else:
            return 0
                
                
        return 0      


    def combine_record(self, dt, container = ''):
        """ This is the main function that analize each record (i.e. separate ascent) and decides which one to keep as merged. 
            Extracs the observations_table, header_table and era5fb accordingly """

        record_dataset_legth ={}     
        other_ds   = []

        ''' I fill the dic e.g. record_dataset_legth{100:['era5_1','ncar'], 80:['bufr','igra2'] }
                i.e. the keys are the lengths, the entries are the lists of datasets '''

        duplicates = []

        found_era5 = False
        if 'era5_1' in container.keys() or 'era5_2' in container.keys() :
            found_era5 = True

        for k in container.keys(): # loop over the dataset
            if k not in other_ds:
                other_ds.append(k)
            for f in container[k]: # loop over the file per dataset
                num_rec = len(container[k][f]['obs_tab']["date_time"])

                """ Storing all the reports id with the proper prefix (different for each different dataset) """
                rep_id = b''.join(container[k][f]["obs_tab"]['report_id'][0]) 
                rep_id = self.observation_ids_merged[k] + rep_id 
                duplicates.append( rep_id )  

                if num_rec not in record_dataset_legth.keys():
                    record_dataset_legth[num_rec] = {}
                    record_dataset_legth[num_rec]['best_ds'] = []
                    record_dataset_legth[num_rec]['file'] = []

                record_dataset_legth[num_rec]['best_ds'].append(k)
                record_dataset_legth[num_rec]['file'].append(f)
                
        entries = list(record_dataset_legth.keys())
        entries.sort(reverse= True)

        if found_era5:                       
            for e in entries:
                best_datasets = record_dataset_legth[e]['best_ds']
                if 'era5_2' in best_datasets:  # era5_1 and era5_2 should never be both present anyway...
                    best_ds = 'era5_2'      
                    break
                elif 'era5_1' in best_datasets:
                    best_ds = 'era5_1'
                    break  # will pick either era5_1 or era5_2 if available

        else:
            for e in entries:
                if 'igra2' in record_dataset_legth[entries[0]]['best_ds']:
                    best_ds = 'igra2' # will pick igra2 if available
                    break
                elif 'igra2' not in record_dataset_legth[entries[0]]['best_ds']:
                    best_ds = record_dataset_legth[entries[0]]['best_ds'][0] # will pick anything else if available
                    break

        best_file = record_dataset_legth[e]['file'][record_dataset_legth[e]['best_ds'].index(best_ds)]

        ''' If more file are available for the same best_ds, pick the first one from the list '''
        selected_obstab, selected_era5fb = container[best_ds][best_file]['obs_tab'] , container[best_ds][best_file]['era5fb_tab']

        ''' Update the data policy licence '''
        if 'ncar' in container.keys() or 'igra2' in container.keys() or 'era5_1759' in container.keys() or 'era5_1761' in container.keys():
            dp = 0 # free data
        else:
            dp = 4 # restricted data
        licence = np.full ( len (selected_obstab['date_time']) , dp , dtype = int)
        selected_obstab['data_policy_licence'] = licence 
        #print('check')
            
        
        ''' Creating the correct observations and record ids. 
                All the bytes variable are shrunk to a long |S1 byte variable type, otherwise 
                writing in h5py will not work. '''

        for var in ['observation_id']:
            if type (selected_obstab[var] ) == np.ndarray and type (selected_obstab[var][0] ) == np.bytes_:
                selected_obstab[var] = np.array ([self.observation_ids_merged[best_ds] + b''.join(l) for l in selected_obstab[var] ] )
            elif type (selected_obstab[var] ) == np.ndarray and type (selected_obstab[var][0] ) == np.ndarray:
                selected_obstab[var] = np.array ([self.observation_ids_merged[best_ds] + b''.join(l) for l in selected_obstab[var][:] ] )

        for var in ['report_id']:
            val = selected_obstab[var][0]
            if type (selected_obstab[var] ) == np.ndarray and type (val) == np.bytes_:
                value = self.observation_ids_merged[best_ds] + b''.join(val)  # it is the same for each row in the table
            elif  type (selected_obstab[var] ) == np.ndarray and type (val) == np.ndarray:
                value = self.observation_ids_merged[best_ds] + b''.join(val)                  
                arr = np.full( (1, len( selected_obstab['date_time']) ) , value )[0] # np.full returns a list of lists

            selected_obstab[var] = arr


        for var in selected_era5fb.keys():
            if type (selected_era5fb[var]) == np.ndarray and type (selected_era5fb[var][0] ) == np.ndarray:
                try:
                    selected_era5fb[var] = np.array( [b''.join(l) for l in selected_era5fb[var][:] ] )
                except:
                    value = [b''.join(l) for l in selected_era5fb[var][0] ][0]
                    selected_era5fb[var] = np.array( (1, len( selected_obstab[var]) ) ).fill(value)

        """ Extracting the header """
        selected_head = self.get_header_table(dt, ds = best_ds, File = best_file )
        for var in selected_head.keys():
            if type (selected_head[var] ) == np.ndarray and type (selected_head[var][0] ) == np.bytes_:
                selected_head[var] = np.array( [b''.join(l) for l in selected_head[var][:] ] )

        if  'best_ds' == 'era5_1' or best_ds == 'era5_2' :
            selected_obstab['advanced_assimilation_feedback'] = np.array([1]*len(selected_obstab['date_time']) )
            selected_obstab['advanced_uncertainty'] = np.array([1]*len(selected_obstab['date_time']) )
            
        else:
            selected_obstab['advanced_assimilation_feedback'] = np.array([0]*len(selected_obstab['date_time']) )
            selected_obstab['advanced_uncertainty'] = np.array([0]*len(selected_obstab['date_time']) )

        #best_ds_byte = np.bytes_(best_ds, ndtype = '|S10') # converting to bytes object
        best_ds_byte = np.bytes_(best_ds) # converting to bytes object            
        arr = np.full( (1, len( selected_obstab['date_time']) ) , best_ds_byte )[0]
        selected_obstab['source_id'] = arr

        duplicate = b','.join(duplicates)
        #selected_head['duplicates'] = np.array(duplicate)

        duplicate = np.array(duplicate).astype(dtype='|S70')
        selected_head['duplicates']               = np.array([duplicate])
        selected_head['report_id']                = np.array([selected_obstab['report_id'][0]])
        selected_head['source_id']                = np.array([selected_obstab['source_id'][0]])
        selected_head['record_timestamp']  = np.array([selected_obstab['date_time'][0]])

        selected_file = np.bytes_(best_file.split('/')[-1])

             
        #print (best_ds , best_file.split("/")[-1], selected_head['source_id'][0].decode("utf-8")  ,  selected_obstab['date_time'][0] )
        return  best_ds, selected_obstab, selected_era5fb, selected_head, selected_file, best_file


    def retrieve_attr_dic(self, content, var):
        attrs_dic = {}
        attrs_dic[var]= {}
        
        try:
            attrs_dic[var]['description']    = bytes( self.dic_type_attributes[content][var]['description']    , 'utf-8' )
        except:
            attrs_dic[var]['description']    = bytes( 'missing'    , 'utf-8' )
            print(' FFF FAILING WITH DESCRIPTION: ', var , ' ' ,  self.dic_type_attributes[content][var]['description']) # todo CHECK WHY SOME ARE FAILING

        try:
            attrs_dic[var]['external_table'] = bytes( self.dic_type_attributes[content][var]['external_table'] , 'utf-8' )
        except:
            attrs_dic[var]['external_table'] = bytes( 'missing' , 'utf-8' )
        
        return attrs_dic
    
    def write_merged(self, content = '', table=''):
        """ Module to write the output file as netCDF """

        if not os.path.isdir(self.out_dir):
            Path(self.out_dir).mkdir(parents=True, exist_ok=True)                  
        out_name = self.out_dir + '/' + self.station + '_CEUAS_merged_v1.nc'  
        self.out_name = out_name 
        
        #self.out_name =  self.station + '_CEUAS_merged_v1.nc' 
        '''
            if os.path.isfile('dic_obstab_attributes.npy'):
                  attrs_dic = np.load('dic_obstab_attributes.npy' , allow_pickle = True).item()
            else:
                  attrs_dic = {}
            '''
        #attrs_dic = {}


        
                
        if content in ['observations_table','header_table','era5fb', 'station_configuration']:
            for var in table.keys():
                if var == 'comments':
                    continue 
                
                attrs_dic = self.retrieve_attr_dic(content, var)
                
        '''  
        """ Retrieving the attributes """
        if content in ['observations_table','header_table','era5fb', 'station_configuration']:
            for var in table.keys():
                if var == 'comments':
                    continue 

                attrs_dic[var] = {}
                try:
                    attrs_dic[var]['description']    = bytes( self.dic_type_attributes[content][var]['description']    , 'utf-8' )
                except:
                    attrs_dic[var]['description']    = bytes( 'missing'    , 'utf-8' )
                    print(' FFF FAILING WITH DESCRIPTION: ', var , ' ' ,  self.dic_type_attributes[content][var]['description']) # todo CHECK WHY SOME ARE FAILING

                try:
                    attrs_dic[var]['external_table'] = bytes( self.dic_type_attributes[content][var]['external_table'] , 'utf-8' )
                except:
                    attrs_dic[var]['external_table'] = bytes( 'missing' , 'utf-8' )
                    #print(' FFF FAILING WITH EXTERNAL TABLE : ', var ) # FFF CHECK WHY SOME ARE FAILING                                                          
        '''

        if content == 'recordindex':  # writing the recordindex, recordtimestamp, dateindex
            #logging.info('Writing the merged record indices to the netCDF output ')
            table.to_netcdf(out_name, format='netCDF4', engine='h5netcdf', mode='a')

        elif content == 'cdm_tables':
            for k in data['cdm_tables'].keys():
                table = data['cdm_tables'][k]
                table.to_netcdf(out_name, format='netCDF4', engine='h5netcdf', mode='a', group = k)
                #logging.info('Writing the cdm table %s to the netCDF output ', k)

        elif content == 'source_configuration':                 
            table.to_netcdf(out_name, format='netCDF4', engine='h5netcdf', mode='a', group = content)
            #logging.info('Writing the source_configuration table to the netCDF output ')

        elif content == 'station_configuration':
            for k in table.keys(): 

                    
                var_type = self.dic_type_attributes[content][k]['type']

                ''' trying to convert the variable types to the correct types stored as attribute, read from the numpy dic file '''
                if type(table[k][0]) != var_type:
                    try:
                        table[k] = table[k].astype( var_type ) 
                        #print('Done station_conf' , k )
                    except:
                        if k == 'secondary_id':
                            table[k] = table[k].astype( bytes ) 

                        #print ('FAILED converting column ' , k, ' type ', type(table[k][0]) , ' to type ', var_type )

                dic = {k:table[k]}  
                write_dict_h5(out_name, dic , content, self.encodings[content], var_selection=[], mode='a', attrs = attrs_dic  )


        # Writing the observations_table, header_table, era5fb 
        elif content in ['observations_table', 'era5fb', 'header_table']: 

            shape = ''
            for k in table.keys(): 
                if k == 'index' or k == 'hdrlen' or 'string' in k :
                    continue


                var_type = self.dic_type_attributes[content][k]['type']

                ''' trying to convert the variable types to the correct types stored as attribute, read from the numpy dic file '''
                if type(table[k][0]) != var_type:

                    if k == 'hdrlen': 
                        continue
                    try:
                        #table[k] = table[k].astype( bytes ) 
                        table[k] = table[k].astype( var_type ) 

                    except:
                        #print ('FAILED converting column ' , k, ' type ', type(table[k][0]) , ' to type ', var_type )
                        pass
                dic = {k:table[k]}  # making a 1 colum dictionary
                shape = table[k].shape
                #print('SHAPE IS FFF ', table[k].shape )
            
                write_dict_h5(out_name, dic , content, self.encodings[content], var_selection=[], mode='a', attrs = attrs_dic  )

            if content == 'observations_table' and not self.obstab_nans_filled :
                missing_cdm_var = [ v for v in self.dic_type_attributes[content].keys() if v not in self.observations_table_vars]  # variables to be filled with nans            
                for k in missing_cdm_var:
                    if k not in ['advanced_assimilation_feedback', 'advanced_uncertainty']:
                        var_type = self.dic_type_attributes[content][k]['type']
                        if var_type == np.int32 :
                            nan = np.int32(-2147483648)
                        else:
                            nan = np.float32(np.nan)       
                        logging.debug('Adding missing cdm colum with empty values: %s ' , k )
                        dic={k:np.empty(shape,dtype=np.dtype(nan))}
                        dic[k].fill(nan)
                        d = self.retrieve_attr_dic(content, k)
                        write_dict_h5(out_name, dic, 'observations_table', self.encodings['observations_table'], var_selection=[], 
                                      mode='a', attrs = d  ) ### TO DO
                self.obstab_nans_filled = True

            elif content == 'observations_table' and self.obstab_nans_filled:
                return

    def merge(self ,  station = '' , datasets = '' , mode = 'test'):        

        """ Call to the merge_all_data() and write_merged_file() methods """

        if mode == "test":      
            a = self.initialize_data( station = station, datasets = datasets ) # reading the input files 
            if not a:
                print("Quitting - nothing I can do")
                return False
            
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
                if not a:
                    print("Quitting - nothing I can do")
                    return False
                
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

       


base_dir = "/scratch/das/federico/COP2_HARVEST_FEB2022/"

data_directories   = { 'era5_1'       : base_dir + '/era5_1'     ,
                                   'era5_2'       : base_dir + '/era5_2'     ,
                                   'era5_3188' : base_dir + '/era5_3188'  ,
                                   'era5_1759' : base_dir + '/era5_1759'  ,
                                   'era5_1761' : base_dir + '/era5_1761'  ,
                                   'ncar'           : base_dir + '/ncar'       ,
                                   'igra2'          : base_dir + '/igra2'      ,
                                   'bufr'            : base_dir + '/bufr'       ,                
                                   }




out_dir = '/scratch/das/federico/MERGED_25FEB2022/'

#out_dir = 'PROVA'

run_mode = 'dummy'

#os.system('rm  ' + out_dir +  '/0-20000-0-82930_CEUAS_merged_v1_beforeSensor.nc' )
#os.system('rm  ' + out_dir +  '/0-20000-0-82930_CEUAS_merged_v1.nc' )


def create_stat_dic(stat_id, data_directories):
    """ Looks in the database for files matching with the given station id.
    Returns a dictionary with the files per dataset,
    the total size of the files t be processed,
    the station_id to be used for the output file. """
    
    station_dic = {}
    total_dim = []
    found_20001 = False
    for d,i in data_directories.items():
        files = os.listdir(i)
        for f in files: # normal check for any file
            Id = f.split('_'+d)[0]
            if 'igra' in Id:
                Id = Id.split('-igra2')[0]
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
                    total_dim. append( os.path.getsize (i + '/' + f) )

    size = sum(total_dim)        
    if found_20001:
        stat_id = stat_id.replace('-20000-','-20001-')
    return station_dic, size, stat_id





if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Make CDM compliant netCDFs")
    parser.add_argument('--stations' , '-s', 
                          help="List of station ids to merge"  ,
                      type = str)

    parser.add_argument('--maximum' , '-m', 
                          help="Maximum size of a single station file "  ,
                      type = int ) # 1 GB

    parser.add_argument('--lightest' , '-l',
                          help="Minimum size of a single station file "  ,
                          type = int ) # 0.1 MB

    parser.add_argument('--run_exception' , '-e',
                          help="Run the station without raising errors"  ,
                          type = bool,
                          default = False ) # 0.1 MB 

    args = parser.parse_args()
    stations = args.stations
    max_size = args.maximum
    min_size = args.lightest
    run_exception = args.run_exception

    """ Initialize the merging class """
    Merging = Merger(out_dir)

    run_exception = False
    print('run exception is ', run_exception )

    out = open('Failed_merged.txt' , 'a+')
    out.write('# Failed stations \n')
    for station in stations.split(','):
        if '-20600-' in station or '-20999-' in station:
            continue 
        station_dic, size, station = create_stat_dic(station, data_directories)

        if size < 1:
            print("No station found, please check the station ID !!!")
            
        elif size < max_size and size > min_size:
            if run_exception:
                try:
                    a = Merging.merge(station,  datasets = station_dic , mode = run_mode ) # single station dictionary
                except:
                    print(' ERROR, SKIPPED !!!')
                    out.write( str(station) + '\n')
            else:
                a = Merging.merge(station,  datasets = station_dic , mode = run_mode ) # single station dictionary                

            print(' --------------------------------> finished station ', station )

            '''
                  try:
                        print(' I am running a total compressed size of ::: ', (size / 10**9) , '  GB   since minimum allowed:' , min_size , '  and maximum allowed: ', max_size , '   station: ' , station ) 
                        a = Merging.merge(station,  datasets = station_dic , mode = run_mode ) # single station dictionary
                  except:
                        print('*** Failed :::  ', station , ' ::: must re-do ***')
                        out.write(station + '\n' )
                  '''
        else:
            print(' The file size ' , size/10**9 , '   is not compatible with minimum allowed:' , min_size , '  and maximum allowed: ', max_size ) 


# run: -s 0-20000-0-82930 -l 100 -m 1000000000000 
# 0-20000-0-82930,0-20000-0-03976 
# -s 0-20000-0-97760 -l 100 -m 1000000000000 
# -s 0-20000-0-72274 -l 100 -m 1000000000000 
