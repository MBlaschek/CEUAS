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
#import h5py as h5py
import xarray as xr 
import time 
#from numba import njit

sys.path.append('../harvest/code')
from harvest_convert_to_netCDF_newfixes import write_dict_h5 

# nan int = -2147483648 
#from harvest_convert_to_netCDF import datetime_toseconds   # importing the function to write files with h5py 


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

def now(time):
      """ For easily print a readable current time stamp """
      a = datetime.fromtimestamp(time).strftime('%Y-%m-%d %H:%M:%S')
      return  a


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
            self.datasets_all = ['igra2' , 'era5_1' , 'ncar_w' ,  'ncar_t', 'bufr' , 'era5_1759' , 'era5_1761' , 'era5_3188']    # all possibly available datasets                          
            #self.observation_ids_merged  = {  'igra2':1 , 'ncar_t':2 , 'ncar_w':2, 'bufr':3,  'era5_1':4 , 'era5_1759' :5 , 'era5_1761':6 ,  'era5_3188' :7}  # values used to convert original record_id to the merged record_id, see method merge_all_data 
            
            self.observation_ids_merged  = {  'igra2':1 , 'ncar':2, 'bufr':3,  'era5_1':4 , 'era5_1759' :5 , 'era5_1761':6 ,  'era5_3188' :7}  # values used to convert original record_id to the merged record_id, see method merge_all_data 
            
            self.unique_dates = {}            
            self.attributes = {} # will keep the original attributes from the CDM tables, read from the netCDF files 
            self.id_string_length = 14 # fixed length for record_id and observation_id values 
            self.out_dir = out_dir 
            self.variable_types = {}
            logging.info('*** Initialising the Merging procedure ***' )   
            
      def initialize_data(self , station = '', datasets = {} ):
            """ Initialize dataset; store relevant data as attributes.
                       Args ::     dic{}  datasets (dictionary where keys are the dataset names e.g. bufr, igra2 etc. , and the value is the path to the corresponding netCDF file 
                                       e.g. tateno = { 'ncar'    : 'example_stations/ncar/chuadb_trhc_47646.txt.nc'   ,
                                                               'igra2'   : 'example_stations/igra2/chJAM00047646-data.txt.nc'  } 
            """       
            self.datasets = datasets
            self.datasets_keys = datasets.keys()
            self.station = station
            self.out_name = self.out_dir + '/' + self.station + '_CEUAS_merged_v0.nc'

            self.observations_table_vars = ['date_time', 'z_coordinate' , 'z_coordinate_type', 'observed_variable' , 'observation_value' , 'report_id' , 'observation_id' , 'latitude' , 'longitude', 'units']

            """ Loading the econding of variables created from the harvester script """
            self.encodings = np.load('groups_encodings.npy' , allow_pickle = True ).item()

            """ Loading the econding of era5fb variables created from the harvester script """
            self.era5fb_encodings = np.load('era5fb_encodings.npy' , allow_pickle = True ).item()

            data = {}   # container for the data of each dataset
            #source_configuration = {}  # container for the source_configuration of each dataset
            
     
            """ Looping over the datasets """
            logging.info('*** Reading and Initializing the data from the netCDF files for the station %s ' , station )
                            
            for k,v in datasets.items() :
                  logging.info(' Initialising the dataset: *** %s ' , k  )
                  data[k] = {} 
                  
            data['cdm_tables'] = {}             
            
            """ Reading the header_table, station_configuration, source_configuration """
            for k,v in datasets.items() :  
                  
                  try:                      
                        d = xr.open_dataset(v , engine = 'h5netcdf' , decode_times = False )    # leaving units as seconds from 1900-01-01 00:00:00 
                        data[k]['recordtimestamp'] = d['recordtimestamp']                  
                        data[k]['recordindex']         = d['recordindex']
                        data[k]['dateindex']            = d['dateindex']
                        
                        logging.debug('Done with %s observations_table' , str(k) )
                  except:
                        print('Failed opening indices: ', v )
                  
                  ########################################################## station_configuration, source_configuration             

                  try:                       
                        d = xr.open_dataset(v , engine = 'h5netcdf' , group = 'observations_table'  , decode_times = False )    
                        
                        """
                        if 'observations_table' not in list( self.attributes.keys() ):  # saving the attributes to be re-applied at the end
                              self.attributes['observations_table'] = {}
                              for var in d.variables  :
                                    if d == 'index':
                                          continue
                                    self.attributes['observations_table'][var] = {}
                                    self.attributes['observations_table'][var]['description'] = d[var].description
                                    try:
                                          self.attributes['observations_table'][var]['external_table'] = d[var].external_table
                                    except:
                                          self.attributes['observations_table'][var]['external_table'] = ''
                        """
                        self.obs_table_columns = d.variables # storing all the variables from the CDM table 
                        
                        if k == list(datasets.keys())[0]:
                              self.get_variable_type(k, d, 'observations_table' )
                        
                        to_drop = [ n for n in d.variables if n not in self.observations_table_vars ] # removing still not implemented variables 
                        d = d.drop( to_drop ) 
                        df = d.to_dataframe()
                        
                        #for c in df.columns:
                        #      print(c , ' ' , type(df[c][:][0]))    
                              
                        data[k]['observations_table'] = df
                        logging.debug('Done with %s observations_table' , str(k) )
                        d.close() 
                  except:
                        print('Failed opening observations_table: ', v )
                        
                        
                  ###   station_configuration
                  try:
                        d = xr.open_dataset(v , engine = 'h5netcdf' , group = 'station_configuration' , decode_times = False )                 
                        #data[k]['station_configuration'] = d.to_dataframe()  
                        data[k]['station_configuration'] = d                    
                        logging.debug('Done with %s station_configuration' , str(k) )
                        d.close()                   
                  except:
                        pass
                  
                  ###   source_configuration
                  try:                        
                        d = xr.open_dataset(v , engine = 'h5netcdf' , group = 'source_configuration' , decode_times = False )                 
                        if len(d) == 0:
                              source_file = v.split('/')[-1]
                              df = pd.DataFrame({'source_file': [ source_file    ]    }).to_xarray() 
                              data[k]['source_configuration'] = df    
                        else:
                              data[k]['source_configuration'] = d                            
                        logging.debug('Done with %s source_configuration' , str(k) )
                        d.close() 
                  except:
                        pass
                  
                  
                  d = xr.open_dataset(v , engine = 'h5netcdf' , group = 'header_table' ,  decode_times = False  )  
                  logging.debug('Loading the header_table')            
                  
                  
                  if 'header_table' not in list( self.attributes.keys() ):  # saving the attributes to be re-applied at the end
                        self.attributes['header_table'] = {}
                        for var in d.variables:
                              try:
                                    self.attributes['header_table'][var] = {}
                                    self.attributes['header_table'][var]['description'] = d[var].description
                                    self.attributes['header_table'][var]['external_table'] = d[var].external_table           
                                    #print(var, ' ' , d[var].description , ' ' , type(d[var].values[0]) , ' ' )
                              except:
                                    print('Cannot extract description or external table for', var )
            
                  if k == list(datasets.keys())[0]:
                        self.get_variable_type(k, d, 'header_table' )                    
                  data[k]['header_table'] = d.to_dataframe()   
                  logging.debug('Done with %s ' , k )
  
                        
                  logging.info("*** Loading the observations_table (might take time) %s" , k )          
                  d.close() 

                  if k == 'era5_1': # reading the whole era5_1 feedback (including reanalysis)
                        d = xr.open_dataset(v , engine = 'h5netcdf' , group = 'era5fb')     
                        data[k]['era5fb'] = d.to_dataframe()                          
                        self.get_variable_type(k, d, 'era5fb' )                        
                        logging.debug('Done with %s era5 feedback ', k )
            
                  """ Reading the CDM tables that do not depend on specific stations or observations (fixed values), for the first file only """           
                  if list(datasets.keys()).index(k) == 0  : # reading the tables from the first dataset file only
                        for t in [ 'crs' , 'observed_variable', 'units' , 'z_coordinate_type' , 'station_type', 'station_configuration_codes']:                                                            
                              d = xr.open_dataset(v , engine = 'h5netcdf' , group = t)                 
                              #data['cdm_tables'][t] = d.to_dataframe()   ### USELESS ?
                              data['cdm_tables'][t]  =    d                                      
                  d.close()    
            
            
            self.data = data     

            
            """ Making all date_times  """   
            self.make_all_datetime()            
            

            
      def get_variable_type(self, dataset, df, table ):
            """ Reads the variable type for all the tables and variable as stored in the netCDF files. """            
            var_type_dic = {}            
            array = df                  
            try:                       
                  for v in array.variables:
                        tipo = type( array[v].values[0] ) 
                        if tipo == np.float64 :  # there should be no float64 types 
                              tipo = np.float32                       
                        var_type_dic[v] = tipo   
            except:                       
                  for v in array.columns:
                        tipo = type( array[v].values[0] )
                        if tipo == np.float64 :
                              tipo = np.float32                                                     
                        var_type_dic[v] = tipo            
            self.variable_types[table] = var_type_dic
            
            
      def make_all_datetime(self):
            """ Building the global set of date_times and indices from the various datasets. 
                  The datetimeindex is read from the original netCDF file. 
                  Will compare the unique date_time of each dataset and extract the global unique date_times
                  """
            
            logging.info('\n *** Running make_all_datetime ' )
            
            all_uniques = []  # storing a list with all the unique date_tmes            
            which_k_in_dt = {}  # list of avilable dataset for each unique date_time, so that when looping over the distinct date_times, only the proper dataset will be read and compared 

   
            for k,v in self.datasets.items() :                  
                  self.unique_dates[k] = {}
                  self.unique_dates[k]['indices'] = {}                             
                                   
                  unique_dt = list(self.data[k]['recordtimestamp'].data.astype(int))
                  indices = list(self.data[k]['recordindex'].data)
                     
                  all_uniques += unique_dt   # adding to the total unique date_times 
            
                  
                  """ Loop over all the date_times of each dataset """
                  for dt, index_low, count in zip (unique_dt,  indices, range(len(unique_dt))  ):
                        
                        try:                          
                              which_k_in_dt[dt].append(k)
                        except:
                              which_k_in_dt[dt] = []
                              which_k_in_dt[dt].append(k)                             
                        
                        self.unique_dates[k]['indices'][dt] = {}
                        self.unique_dates[k]['indices'][dt]['low'] = index_low                       
                        try:
                              index_up =  indices[ count + 1 ]  # works until the last available recordindex
                        except:                             
                              #index_up = len(indices-1)    
                              index_up = len(indices)-1    
                              
                        self.unique_dates[k]['indices'][dt]['up'] = index_up
                             
                  #self.unique_dates[k]['indices'].append(index) 
                  #self.unique_dates[k]['indices_up'].append(index_up) 
                                    
            self.dataset_per_dt = which_k_in_dt             
            self.merged_unique_dates = np.unique(np.array(all_uniques) )  # storing the set of all distinct dt values            
            logging.debug('make_all_datetime finished ')         
      
 
      
      def clean_dataframe(self, df_in , what = ''):  # NOT IN USE 
            """ Remove empty (nan) or wrong values from the original data """         
            
            if what == 'era5fb':  # cleaning the era5 feedback only 
                  df = df_in[np.isfinite(df_in['obsvalue@body'])]
                  try:                        
                        df = df.loc[ df['vertco_type@body'] != 2 ]   
                  except:
                        pass
                  df = df.reindex()
                  df = df[np.isfinite(df['vertco_reference_1@body'])]
                  #print('check lengths: ' , len(df_in) , len(df) )
                  new_ind = np.array ( range(len(df))) 
                  df['index'] =new_ind
                  df = df.set_index('index')
                  
            else:             
                  ### check if can be optimized ???
                  df =  df_in.loc[ df_in['z_coordinate_type'] != 2 ]  # case where the levels are given in terms of geopotential only (pressure not available)
                  
                  df = df.loc[ (df['observation_value'] != -99999.0) 
                                       & (df['observation_value'] != -999.0) 
                                       & (df['observation_value'] != -9999)                                        
                                       & (df['observation_value'] != -9999.0) 
                                       & (df['observation_value'] != -999.9) 
                                       & (df['observation_value'] != -8888 )
                                       & (df['observation_value'] != -8888.0 )
                                       
                                       #& (df['z_coordinate_type'] != 2)  
                                       & (df['z_coordinate'] != -99999.0) 
                                       & (df['z_coordinate'] != -9999.0 )
                                       & (df['z_coordinate'] != 999 )
                                       & (df['z_coordinate'] != 999.0 )
                                          
                                          
                                       ] #cleaning the values                        
                  #clean = clean.loc[ (clean['z_coordinate_type'] != 2)] #cleaning the values
                  #clean = clean.loc[ (clean['z_coordinate'] != -99999.0 )] #cleaning the values
                  df = df[np.isfinite(df['observation_value'])]  # excluding nan values 
                  df = df[np.isfinite(df['z_coordinate'])]                  
            return df 
      
   

      def get_reanalysis_feedback(self, dt_bestds_dic, reanalysis = ''):
            """ Extracts the renanalysis feedback from the dt_bestds dictionary.
                  Currently only era5fb available (era5_1 dataset). """
            if 'era5_1' in self.datasets_keys:

                  fb = xr.open_dataset(self.datasets['era5_1'] , engine = 'h5netcdf' , group = 'era5fb')
                  self.get_variable_type('era5_1', fb, 'era5fb' )
                  era5fb = fb.to_dataframe()
                  fb_columns = list(era5fb.columns)
                  fb.close()
                  logging.debug('Loaded era5 feedback ')

            else:
                  fb_columns = ['empty']


            feedbacks = [] # container for each feedback
            for dt, obstab in dt_bestds_dic.items():
                  dataset = obstab['best_ds']
                  len_obstab = obstab['len']
                  if dataset == 'era5_1' and reanalysis == 'era5fb':   # reading the feedback for era5_1 
                        index_low = self.unique_dates[dataset]['indices'][dt]['low']
                        index_up  = self.unique_dates[dataset]['indices'][dt]['up']
                        fb = era5fb.iloc[index_low:index_up]                 
                        feedbacks.append(fb)
                        #input('check FB', fb)
                  else:  # empty feedback otherwise 
                        """ Creating empty feedback tables. 
                        If era5_1 is never used as best_ds, the era5fb tables are never used and will store an empty, single column DF as feedback. 
                        If era5_1 is used somewehre and its own fb is used, we must use the format of that df (from the odb files) also for the other empty fb """
                  
                        #len_odf = len(merged_observations_table)
                        empty_fb = pd.DataFrame(np.nan, columns= fb_columns ,  index = range(len_obstab) )                  
                        feedbacks.append(empty_fb)
                        
            merged_fb = pd.concat(feedbacks)
            return merged_fb
                                               
      def get_header_table(self , dt, ds = '' , all_ds = '', length = '', report_id = '' , obs_tab = '' ):
            """ Extracting the header_table, and replacing the "duplicates" entries with the list of alternative available datasets """
            #index_low = self.unique_dates[ds]['indices'][dt]['low']
            #index_up  = self.unique_dates[best_ds]['indices'][dt]['up']            
            
            index = np.where( self.data[ds]['recordtimestamp'][:] == dt)[0][0]
            #data[k]['recordtimestamp'] = d['recordtimestamp']                  
            #data[k]['recordindex']         = d['recordindex']
            #data[k]['dateindex']            = d['dateindex']
            
            hd = self.data[ds]['header_table'].iloc[index]  
            #hd = hd.replace(-2147483648 , np.nan )
            
            frame = {}
            for c in self.data[ds]['header_table'].columns:
                  frame[c] = hd[c]
                  
                  
            frame['duplicates'] = all_ds 
            
            """ Copying variables from observations table """
            frame['report_id'] = obs_tab['report_id']
            source_id = ds.rjust(10) 
            frame['source_id'] = source_id 
            frame['report_timestamp'] = obs_tab['date_time'] 
            
            hd = pd.DataFrame(frame)
            
            hd['source_id'] = (hd['source_id']).astype('S10') 
            
            #hd = hd.replace(-2147483648 , np.nan )
            
            #for v in hd.columns:
            #      if v in self.data[ds]['station_configuration'].variables:
            #                  hd[v] = self.data[ds]['station_configuration'][v].values[0].astype( self.variable_types['station_configuration'][v])
            return hd


      def merge_all_data(self):       
            """ Construct a dictionary with all the dataframes of each dataset, either reading it from saved pickles or reading it from memory """
            
            logging.info('***** Starting the merging process ')

            """ All possible unqiue_dates to loop on """
            date_times = self.merged_unique_dates
            date_times.sort()
            date_times = np.array(date_times) 
                                        
            """ List storing the indices of the date_index of the merged dataset """
            all_merged_obs ,  all_merged_head, all_merged_fb , merged_indices , merged_date_time, merged_indices_out , mi = [] , [] , [] , [] , [], [], [] 
                     
            """ Dictionary that will contain the merged file. """            
            # rand = datetime.strptime('1981-01-03 12:00:00', '%Y-%m-%d %H:%M:%S')  
            #for dt in date_times[3008:3100]: # loop over all the possible date_times 
            dt_bestds_dic = {} # store the selected best dataset for each dt             
            tot = len(date_times)
            for dt, c in zip(date_times, range(tot) ): # loop over all the possible date_times 
                  print('Analize : ', str(c+1) , '/',  str(tot)  , ' ', dt , ' ', now(time.time()) )
                  
                  #logging.info('Analize : %s %s /', str(c+1) ,  str(tot)  )
            
                  cleaned_df_container = {}                  
                  
                  for k in self.dataset_per_dt[dt] :  # checking the list of available datasets  
                                                
                        index, index_up = self.unique_dates[k]['indices'][dt]['low'] , self.unique_dates[k]['indices'][dt]['up']  # extracting the exact chunk of the dataframe where the data of this are stored                           
                        chunk = self.data[k]['observations_table'].iloc[index:index_up]
                        #chunk = chunk.replace(-2147483648 , np.nan )     

                        if len(chunk)==0:
                              continue
                        
                        cleaned_df_container[k] = {}                                                
                        cleaned_df_container[k]['df']                   = chunk         # cleaned dataframe 
                  
                  if all(value == 0 for value in cleaned_df_container.values()):
                        logging.debug('No data were found! ')
                        continue
                  
                  merged_observations_table, best_ds, duplicates, header = self.merge_record(dt, container = cleaned_df_container)
                  dt_bestds_dic[dt] = {}
                  dt_bestds_dic[dt]['best_ds'] = best_ds
                  dt_bestds_dic[dt]['len'] = len(merged_observations_table)

                  """ Adding source_id to the observations_table equal to the best_ds chosen """
                  merged_observations_table['source_id'] = best_ds.rjust(10)                             
                  merged_observations_table['source_id'] = merged_observations_table['source_id'].astype('S10')
                  all_merged_obs.append(merged_observations_table)

                  #OLD VERSION """ Extracting the merged feedback, flagging the advanced_observations_feedback flag = 1 """
                  #feedback, merged_obs = self.get_reanalysis_feedback( dt, merged_observations_table , reanalysis='era5fb', best_ds= best_ds)
                  #all_merged_fb.append(feedback)                  
                  #all_merged_obs.append(merged_obs)
                  
                  """ Setting the correct report_id in the header table """
                  merged_report_id = merged_observations_table['report_id'].values[0]  # same report_id as calculated in the observation_table 
                  header['report_id'] = merged_report_id 
                  all_merged_head.append(header)
                                    
                  """ New merged recordindex and recordtimestamps indices """
                  merged_indices.append(len(merged_observations_table))                 
                  merged_date_time.append(dt)


            """ Storing the merged date_time values and indices """
            di=xr.Dataset()
            merged_date_time = np.array(merged_date_time)
            di['recordtimestamp'] = ( {'recordtimestamp' : merged_date_time.shape } , merged_date_time )
            di['recordtimestamp'].attrs['units']='seconds since 1900-01-01 00:00:00'
                     
                     
            """ Creating the merged indices """
            mi.append(0)
            for i in range(len(merged_indices)):
                  mi.append( merged_indices[i] + mi[-1] )
            del mi[-1]       
            mi = np.array(mi) 
            di['recordindex']          = ( {'recordindex' : mi.shape } , mi )

                                                      
            """ Creating the merged dataframes """
            logging.debug('*** Concatenating the observations_table dataframes' )      
            merged_obs = pd.concat(all_merged_obs)


            del all_merged_obs

            logging.debug('*** Finished concatenating theobservations_table  dataframes' )
            self.write_merged(content = 'observations_table', table= merged_obs)
            #del merged_obs 
        
            dateindex = np.array( [ i // 86400 for i in merged_obs['date_time'] ] )            
            date_times, indices, counts = np.unique(dateindex, return_counts = True, return_index= True)
            di['dateindex'] = ( {'dateindex' : indices.shape } , indices )  # considers the day only
            
            self.write_merged(content = 'recordindex', table = di)                                                                                      
            del di , merged_obs

            logging.debug('*** Concatenating the header_table dataframes' )      
            merged_hd = pd.concat (all_merged_head)
            #merged_hd = merged_hd.replace( -2147483648 , np.nan )           
            self.write_merged(content = 'header_table', table= merged_hd)                          
            logging.debug('*** Finished concatenating the header_table dataframes'  )  
            del merged_hd

            """ Make and write feedback tables """
            fb_tab = self.get_reanalysis_feedback(dt_bestds_dic, reanalysis = 'era5fb')            
            del dt_bestds_dic
            self.write_merged(content = 'era5fb', table= fb_tab)
            del fb_tab

            return 0      
      
      
      def add_cdm_missing_columns(self, all_merged_obs , table = 'observations_table'):
            """ Add the CDM observations_table columns for which no data are available at the end of the merging """
            #cdm_keys = self.obs_table_columns 

            for k in self.obs_table_columns:
                  if self.variable_types[table][k] == np.int32 :
                        nan = -2147483648
                  else:
                        nan = np.float32(np.nan) 

                  if k not in list(all_merged_obs.columns ):
                        logging.debug('Adding missing cdm colum with empty values: %s' , k )
                        all_merged_obs[k] = nan 
                        
            return all_merged_obs
         
            
      def merge_record(self, dt, container = ''):
            """ This is the main function that analize each record (i.e. separate ascent) and decides which one to keep as merged. """            
            record_dataset_legth ={}     
            
            
            """ Combining the ncar_t and ncar_w files.
                  If both are present, select the ncar_t data and rename it as 'ncar'. 
                  If only one is present, simply rename it as 'ncar'. 
            """           
            if ('ncar_t' in list(container.keys())  ):
                  container['ncar'] = {}                                    
                  container['ncar']['df']  = container['ncar_t']['df'] 
                  
            elif ( 'ncar_w' in list(container.keys()) and 'ncar_t' not in list(container.keys())  ) :
                  container['ncar'] = {}                  
                  container['ncar']['df']  =  container['ncar_w']['df'] 

            
            for k in container.keys():
                  if k == 'ncar_t' or k == 'ncar_w': 
                        continue 
                  record_dataset_legth[k] = len(container[k]['df'] )
                
                
            """ For now, choosing the dataset with more records of all or igra2>ncar>rest  data if available and with same number of records """
            best_ds, all_ds , best_datasets, all_ds_reports = 'dummy' , [] , [], [] # total number of records, name of the chosen dataset , list of other possible dataset with available data  
            
            most_records = max( [ v for v in  record_dataset_legth.values() ] )  # maximum number of records per date_time           
            
            for k, v in record_dataset_legth.items():                 
                  if v == 0:
                        continue
                  if v == most_records:
                        best_datasets.append(k)                                  
                  if v > 0:
                        all_ds.append(k) # all other datasets with smaller number of records than the maximum found              
                        try: 
                              all_ds_reports.append( self.observation_ids_merged[k] * 1000000000  + int(container[k]['df']['report_id'].values[0].decode('utf-8')) )  # converting the original report id using the same convention as for observation_id
                        except:
                              #all_ds_reports.append( self.observation_ids_merged[k] * 1000000000  + int( ( container[k]['df']['report_id'].values[0]).decode('latin1')   ) )  # converting the original report id using the same convention as for observation_id
                              #input('check wrong!')
                              all_ds_reports.append( np.nan)  # converting the original report id using the same convention as for observation_id
                              
                              
                              #all_ds_reports.append(np.nan)
                              #print ( type(container[k]['df']['report_id'].values) )
                              #all_ds_reports.append( self.observation_ids_merged[k] * 1000000000  + float(container[k]['df']['report_id'].values[0].decode('latin1') ))
                              
            if len(best_datasets) ==0:
                  print('wrong??? please check')
                  return 0,0,0,0        
   
            if 'igra2' in best_datasets:
                  best_ds = 'igra2'
            elif 'ncar' in best_datasets:
                  best_ds = 'ncar'
            elif 'era5_1' in best_datasets:
                  best_ds = 'era5_1'                  
            else:
                  best_ds = best_datasets[0]
            
            """ Extract container """            
            selected_df = container[best_ds]['df'].copy(deep = True)  # might take extra time, dont know how to get rid of this 

            try:
                  merged_report = self.observation_ids_merged[best_ds] * 10**(self.id_string_length-1)  + int(container[k]['df']['report_id'].values[0].decode('utf-8') )  
            except:
                  input('check what is wrong with the report_id')                  
                  merged_report = np.nan 

            """ Calculate new unique observation id """
            try: 
                  obs_ids_merged =  [   self.observation_ids_merged[best_ds] * 10**(self.id_string_length-1)  + int( i.decode('utf-8')  ) for i in selected_df['observation_id'].values   ]
            except:
                  input('check what is wrong with the observation_id')
                  obs_ids_merged =  [ np.nan for i in selected_df['observation_id'] ]
                  
            obs_ids_merged =  np.bytes_(obs_ids_merged)
            selected_df['observation_id'] = obs_ids_merged
            #selected_df['observation_id'] = np.char.zfill ( selected_df['observation_id'].astype('S' + str(self.id_string_length-1) ) , self.id_string_length )
            
            

            """ Calculate new unique report id """            
            selected_df['report_id'] = merged_report
            selected_df['report_id'] = selected_df['report_id'].astype('S' + str(self.id_string_length) )

            """ Returning a string with the alternative available datasets data """
            if len(all_ds_reports) > 1: 
                  duplicates =   ",".join( [ str(i) for i in all_ds_reports] )
            else:
                  duplicates = str(all_ds_reports[0])
 
                  
            """ Extracting the merged header_table.
                  Again, must consider the special case where best_ds == ncar. 
                  Note that the header table *should* be identical for ncar_w or ncar_t """            
            if best_ds != 'ncar':
                  header = self.get_header_table(dt, ds= best_ds,  all_ds = duplicates , length= len(selected_df) , obs_tab = selected_df  )
                  
            elif ( best_ds == 'ncar' and  'ncar_t' in list(container.keys()) ) :
                  header = self.get_header_table(dt, ds = 'ncar_t', all_ds = duplicates, length= len(selected_df) , obs_tab = selected_df )
                  
            elif ( best_ds == 'ncar' and 'ncar_t' not in list(container.keys())  ) :
                  header = self.get_header_table(dt, ds = 'ncar_w', all_ds = duplicates, length= len(selected_df) , obs_tab = selected_df )            
                  
            logging.debug('I use %s  record since it has more entries: %s but other available datasets are : %s' , best_ds , str(most_records) ,  all_ds ) 

            """ Setting the flag for the presence of feedback information """
            if best_ds == 'era5_1':  
                  selected_df['advanced_assimilation_feedback'] = 1
            else:
                  selected_df['advanced_assimilation_feedback'] = 0

            #print ('duplicates are: ', duplicates)
            #best_ds = best_ds.rjust(10)   
            return  selected_df, best_ds , duplicates, header
      

      def write_merged(self, content = '', table=''):
            """ Module to write the output file as netCDF """
            
            #out_name = os.getcwd() + '/FAST_INDEX_merged_' + [ x for x in self.datasets[ list(self.datasets_keys)[0]].split('/') if '.nc' in x   ] [0]             
            #""" Loading the econding of variables created from the harvester script """
            #encodings = np.load('groups_encodings.npy' , allow_pickle = True ).item()
            
            #""" Loading the econding of era5fb variables created from the harvester script """
            #era5fb_encodings = np.load('era5fb_encodings.npy' , allow_pickle = True ).item()                                    
            if not os.path.isdir(self.out_dir):
                  Path(self.out_dir).mkdir(parents=True, exist_ok=True)                  
            out_name = self.out_dir + '/' + self.station + '_CEUAS_merged_v0.nc'  

            if os.path.isfile('dic_obstab_attributes.npy'):
                  attrs_dic = np.load('dic_obstab_attributes.npy' , allow_pickle = True).item()
            else:
                  attrs_dic = {}
            
            if content == 'recordindex':  # writing the recordindex, recordtimestamp, dateindex
                  logging.info('Writing the merged record indices to the netCDF output ')
                  di = table
                  di.to_netcdf(out_name, format='netCDF4', engine='h5netcdf', mode='a')


            elif content == 'standard_cdm': # writing the station_configuration and source_configuration
                  logging.info('Writing the station_configuration and source_configurations tables to the netCDF output via xarray ')                             
                  for k in self.data.keys():                                                                                                                       
                        if k == 'cdm_tables':                                                                                                                      
                              continue                                                                                                                           
 
                        group_name = k + '_station_configuration'                                                                                                 
                        sc = self.data[k]['station_configuration']                                                                                              
                        sc.to_netcdf(out_name, format='netCDF4', encoding= self.encodings['station_configuration'] , engine='h5netcdf', mode='a' , group = group_name ) 

                        group_name = k + '_source_configuration'                                                                                                
                        sc = self.data[k]['source_configuration']                                                                                               
                        sc.to_netcdf(out_name, format='netCDF4', encoding=self.encodings['source_configuration'] , engine='h5netcdf', mode='a' , group = group_name )  
                                                                                                                                                            
                        del self.data[k]['source_configuration'] , self.data[k]['station_configuration']                                
                  
            
            
            #elif contet == 'fixed_cdm':
            #     logging.info('Writing the standard cdm tables to the netCDF output ')                  
            #     for t in self.data['cdm_tables'].keys():   
                  #if t in ['units']:
                  #      continue   # TODO FIX !!!!
            #          d = self.data['cdm_tables'][t]
                  
            #          try: 
            #              d.to_netcdf(out_name, format='netCDF4', engine='h5netcdf', mode='a'  , group = t )
            #          except:
            #              print('exception: ' , t)
                        #pass
                        #d.to_netcdf(out_name, format='netCDF4', engine='h5netcdf', mode='a'  , group = t )
            #     del self.data['cdm_tables']
            
            elif content == 'observations_table': # header table
                  logging.info('Writing the observations_tables to the netCDF output via xarray to_netcdf() ')
                  obs_tab = table             
                  obs_tab = self.add_cdm_missing_columns(obs_tab)     
            
                  #if os.path.isfile('dic_obstab_attributes.npy'):
                  #      attrs_dic = np.load('dic_obstab_attributes.npy' , allow_pickle = True).item()
                  #else:
                  #      attrs_dic = {}
                  
            
                  for k in obs_tab.columns:
                        if k == 'index' :
                              continue
                        if not type(obs_tab[k].values[0]) == self.variable_types['observations_table'][k]:
                              #print ('Converting column' , k , '   type  ', type(obs_tab[k].values[0]) , ' to type   ' ,  self.variable_types['observations_table'][k]  )
                              obs_tab[k] =  obs_tab[k].astype(self.variable_types['observations_table'][k] ) 
                        #else:      
                             #print ('ALREADY OK ' , k , '   type  ', type(obs_tab[k].values[0]) , ' to type   ' ,  self.variable_types['observations_table'][k]  )                              
                        
                  #print('Writing the observation table using h5py new method for the variable: ' , k )
                        df = obs_tab[ [k] ]  # making a 1 column dataframe  
                        write_dict_h5(out_name, df, 'observations_table', self.encodings['observations_table'], var_selection=[], mode='a', attrs = attrs_dic  ) ### TO DO !!! 
                            
            
            
            elif content == 'header_table':
                  logging.info('Writing the header_table to the netCDF output via xarray ')
                  head_tab = table
                  for v in head_tab.columns:   
                        if v == "index" or v == "hdrlen" or v == "string80":
                              continue            
                  
                        if not type(head_tab[v].values[0]) == self.variable_types['header_table'][v]:
                              try:                              
                                    head_tab[v] =  head_tab[v].astype(self.variable_types['header_table'][v] ) 
                              except:
                                    print ('FAILED converting column' , v , '   type  ', type(head_tab[v].values[0]) , ' to type   ' ,  self.variable_types['header_table'][v]  )                              
                              
                        
                        if 'time' not in v and type(head_tab[v].values[0]) == np.int64 :
                              head_tab[v] = head_tab[v].astype(np.int32)
                        if type(head_tab[v].values[0]) == np.float64:
                              head_tab[v] = head_tab[v].astype(np.float32)
                        a = head_tab[ [v] ]  # making a 1 column dataframe  
                        try: 
                              write_dict_h5(out_name, a, 'header_table', self.encodings['header_table'], var_selection=[], mode='a', attrs = attrs_dic  ) ### TO DO !!!             
                        except:
                              print('ERROR with attributes codec')
                              write_dict_h5(out_name, a, 'header_table', self.encodings['header_table'], var_selection=[], mode='a', attrs = {} )

            elif content == 'era5fb':
                  logging.info('Writing the merged feedback to the netCDF output ')      
                  group_name = 'era5fb'        
                  era5fb = table
                  era5fb = era5fb.reset_index(drop = True)
             
                  
                  for k in era5fb.columns:
                        if k == 'source@hdr':
                              continue 
                        if k == 'index' :
                              continue 
                        #print('Writing the era5fb table using h5py new method for the variable: ' , k )
                        df = era5fb[ [k] ]  # making a 1 column dataframe  
                        try:
                              write_dict_h5(out_name, df, 'era5fb', self.era5fb_encodings,  var_selection=[], mode='a' ) ### TO DO !!! 
                        except:
                              print('Could not read the ' , k , ' column of the era5fb! ')
                        
                        #del self.MergedFeedback , era5fb
            


  
            
             
      def merge(self ,  station = '' , datasets = '' , mode = 'test'):        
            
            """ Call to the merge_all_data() and write_merged_file() methods """
            
            if mode == "test":      
                  a = self.initialize_data( station = station, datasets = datasets ) # reading the input files 
                  dummy = self.merge_all_data()            
                  #logging.info('*** Finished merging, now writing the output netCDF file ***' )         
                  #a = self.write_merged_file()
                  #logging.info('*** Done writing the output ! ***')
                  return True
            
            else:
                  o = open("FAILED_MERGING_LIST.txt", 'a+')          
                  try:
                        a = self.initialize_data( station = station, datasets = datasets ) # reading the input files 
                        dummy = self.merge_all_data()            
                        #logging.info('*** Finished merging, now writing the output netCDF file ***' )         
                        #a = self.write_merged_file()
                        #logging.info('*** Done writing the output ! ***')
                        return True          
                  except:
                        print('Failed: ' , station )
                        o.write(station + '\n' )
                        return False 




""" Example input dictionary """


small = {   'ncar_w'           : '/raid8/srvx1/federico/GitHub/August_develop/CEUAS/CEUAS/cdm/code/23JAN2020_TEST/ncar/chuadb_windc_82930.txt.nc'       ,
                  'ncar_t'            : '/raid8/srvx1/federico/GitHub/August_develop/CEUAS/CEUAS/cdm/code/23JAN2020_TEST/ncar/chuadb_trhc_82930.txt.nc'       ,
                  'igra2'              : '/raid8/srvx1/federico/GitHub/August_develop/CEUAS/CEUAS/cdm/code/23JAN2020_TEST/igra2/chBRM00082930-data.txt.nc'  ,
                  'era5_1'            :  '/raid8/srvx1/federico/GitHub/August_develop/CEUAS/CEUAS/cdm/code/23JAN2020_TEST/era5_1/chera5.conv._82930.nc',
                  'era5_1759'      : '/raid8/srvx1/federico/GitHub/August_develop/CEUAS/CEUAS/cdm/code/23JAN2020_TEST/era5_1759/chera5.1759.conv.1:82930.nc',
                  'bufr'                : '/raid8/srvx1/federico/GitHub/August_develop/CEUAS/CEUAS/cdm/code/23JAN2020_TEST/bufr/chera5.82930.bfr.nc',                          }









#out_dir = '/raid60/scratch/federico/MERGED'

""" Select the output directory """
out_dir = '/raid60/scratch/federico/20MARCH2020_SmallStations'



#out_dir = 'output_test'
#os.system('rm -r output_test' )

if not os.path.isdir(out_dir):
      os.mkdir(out_dir)
""" Dictionary containing the data directories of each dataset """


"""
data_directories = { 'bufr'  : 'example_stations/bufr/' ,
                     'era5_1': 'example_stations/era5_1/', } 
"""

data_directories = {  'era5_1'       :  '/raid60/scratch/federico/netCDF_converted_Jan2020/ready_for_merging/era5_1'       ,
                                  'era5_1759' : '/raid60/scratch/federico/netCDF_converted_Jan2020/ready_for_merging/era5_1759'   ,
                                  'era5_1761' : '/raid60/scratch/federico/netCDF_converted_Jan2020/ready_for_merging/era5_1761'   ,
                                  'era5_3188' : '/raid60/scratch/federico/netCDF_converted_Jan2020/ready_for_merging/era5_3188'   ,     
                                  
                                  
                                  'ncar'          : '/raid60/scratch/federico/12MARCH2020_harvested/ncar'    ,
                                  'igra2'         : '/raid60/scratch/federico/12MARCH2020_harvested//igra2'  , 
                                  'bufr'           : '/raid60/scratch/federico/12MARCH2020_harvested/bufr'         }







def create_station_dics(data_directories):
      """ Create a list of dictionaries, one for each station to be merged, from the database directories. """
      
      files_all = {} 
      for k,v in data_directories.items() :
            files = os.listdir(v)
                              
            for f in files:
                  station = f.split('_')[0] 
                  if 'uknown' in station:
                        continue
                  if station not in files_all.keys():
                        files_all[station] = {}
                  
                  if k == 'ncar':  # separating ncar temperature and wind files 
                        if 'trhc' in f:
                              d = 'ncar_t'
                        elif 'windc' in f:
                              d = 'ncar_w'                          
                        files_all[station][d] = v + '/' + f   # complete path to the netCDF file 
                  else:
                        files_all[station][k] = v + '/' + f   # complete path to the netCDF file 
                        
                  #print('check')     
      return files_all


def get_processed_stations(out_dir):
      """ Retrieve a list of already processed stations in the chosen output directory """
      lista = [ f.split('_')[0] for f in os.listdir(out_dir) if '.nc' in f  ]
      #print('Skipping these stations: ' , lista )
      return lista 
      
      
""" main block """
if __name__ == '__main__':

      """ Parsing the string containing the list of stations to be merged """ 
      parser = argparse.ArgumentParser(description="Make CDM compliant netCDFs")
      parser.add_argument('--stations' , '-s', 
                    help="List of stations id separated by commas, e.g. -f 0-20000-0-48649,0-20000-0-98642"  ,
                    type = str)
      
      parser.add_argument('--dim' , '-d', 
                    help="Maximum size of the data to be merged. Will skip merging if the total size of files exceeds this amount. "  ,
                    default = 5000000,
                    type = float)
      
      parser.add_argument('--mode' , '-m', 
                    help="Switch on the test running mode."  ,
                    default = 'test',
                    type = str)
      
      
      args = parser.parse_args()
      stations            = args.stations
      maximum_size = args.dim
      run_mode        = args.mode
      
      stations = stations.split(',')


      if run_mode == 'test':            
            out_dir = 'output_test'
            os.system('rm -r  ' + out_dir )
            os.system('mkdir output_test')

            Merging = Merger(out_dir)
            
            base_dir = '/raid8/srvx1/federico/GitHub/CEUAS_master_FEB202/CEUAS/CEUAS/public/merge/example_stations'
            data_example = {  'ncar'            : base_dir + '/ncar'              ,
                                              'igra2'           : base_dir + '/igra2'          ,
                                              'era5_1'        : base_dir + '/era5_1'        ,
                                              'era5_1759'  : base_dir + '/era5_1759'  ,
                                              'bufr'            : base_dir + '/bufr'                     }
            
            all_stations_dics = create_station_dics(data_example)
            
      else:
            Merging = Merger(out_dir)
            
            all_stations_dics = create_station_dics(data_directories)
      
      
  
      
      
      """ Running looping over the stations in the list provided  """      
      for stat in stations:
            if len(stat) < 2:
                  print('Only one station, no need for merging, skipping! ' , stat )
                  continue
            processed_stations = get_processed_stations(out_dir)
            
            if stat in processed_stations:
                  logging.info('Skipping already processed station: %s. Please select a different output directory. ' , stat )
                  continue
            
            stat_dic  = all_stations_dics[stat]
            #datasets = list(stat_dic.keys())
           
            if len (stat_dic.keys()) <2: # in case there is only one station, and nothing to be merged 
                  continue             

            """ Calculating the total dimension of the input files to be merged """
            total_dim = []                  
            for s,p in stat_dic.items():
                  total_dim. append( os.path.getsize (p) )
                  total = sum(total_dim)
            
            if total < maximum_size:                  
                  logging.info('The total dimension of the loaded files is: %s [Gb] ' , total/10**9. )    
                  a = Merging.merge(stat,  datasets = stat_dic , mode = run_mode ) # single station dictionary 
                  
            else:                  
                  logging.info('The total size: %s [Gb] is greater than %s , so skipping the station %s *** ' , total/10**9. , maximum_size/10**9., stat )          
                  
            logging.info('*** Convertion completed! ***')


""" To test a single example station: 
# -s 0-20000-0-82930 -d 1000000000000 
""" 
