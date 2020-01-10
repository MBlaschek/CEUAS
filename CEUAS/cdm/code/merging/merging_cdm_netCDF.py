""" Merging the station configuration files """

import os
import sys
import netCDF4 as nc
import pandas as pd
pd.set_option('display.max_rows', 50)
pd.set_option('display.max_columns', 20)
pd.set_option('display.width', 300)

import numpy as np
#import argparse
from datetime import datetime, timedelta
import numpy.ma as ma
#import math
import h5py as h5py
import xarray as xr 
import time 
#from numba import njit

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning) # deactivates Pandas warnings 

import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')  # up to INFO level, DEBUG statements will not be printed 





def now(time):
      """ For easily print a readable current time stamp """
      a = datetime.fromtimestamp(time).strftime('%Y-%m-%d %H:%M:%S')
      return  a


class Merger():
      """ Main class for the merging of the data from different netCDF files. """
      
      def __init__(self ):
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
              
      def initialize_data(self, datasets = {} ):
            """ Initialize dataset; store relevant data as attributes.
                       Args ::     dic{}  datasets (dictionary where keys are the dataset names e.g. bufr, igra2 etc. , and the value is the path to the corresponding netCDF file 
                                       e.g. tateno = { 'ncar'    : 'example_stations/ncar/chuadb_trhc_47646.txt.nc'   ,
                                                               'igra2'   : 'example_stations/igra2/chJAM00047646-data.txt.nc'  } 
            """
        
            data = {}   # container for the data of each dataset
            source_configuration = {}  # container for the source_configuration of each dataset
            
            self.datasets = datasets
            self.datasets_keys = datasets.keys()
            
            """ Looping over the datasets """
            logging.info('*** Reading and Initializing the data from the netCDF files ')
                           
            for k,v in datasets.items() :
                  logging.info('Initialising the dataset: *** ' , k  )
                  data[k] = {} 
                  data['cdm_tables'] = {} 
                  
                  ### alternative with xarray                  
                  #ds =  xr.load_dataset(v)   
                  #observations_table =  xr.open_dataset(v , engine = 'h5netcdf' , group = 'observations_table')                                       
                  
                  ### alternative with netCDF4
                  #ds =  nc.Dataset(v)                   
                  #data[k]['dateindex'] = ds.variables['dateindex'][0,:]  # storing the dateindex 
                  
                  ###for h5py but cant extract date time units !!!
                  ds =  h5py.File(v , driver="core" , )   
                  data[k]['df'] = ds # storing the entire file                
                  data[k]['source_file']           = ds['source_configuration']['source_file'][0]
                  #data[k]['product_code']       = ds['source_configuration']['product_code'][0]  
                  data[k]['recordtimestamp'] = ds['recordtimestamp'].value
                  data[k]['recordindex']         = ds['recordindex'].value                                    
                  #ds.close()                 
                  logging.debug('Reading the file with xarray ')
             
                  
            # add here appending datasets for the case of   ncar_w   and   ncar_t   
            
            
            self.data = data
            self.make_dataframe()
            ds.close()                 

            """ Reading the header_table, station_configuration, source_configuration """
            for k,v in datasets.items() :   
                  d = xr.open_dataset(v , engine = 'h5netcdf' , group = 'station_configuration')                 
                  data[k]['station_configuration'] = d.to_dataframe()   
                  logging.debug('Done with ', k , ' station_configuration')
                  
                  d = xr.open_dataset(v , engine = 'h5netcdf' , group = 'header_table')                 
                  data[k]['header_table'] = d.to_dataframe()   
                  logging.debug('Done with ', k , ' header_table')
                  
                  d = xr.open_dataset(v , engine = 'h5netcdf' , group = 'source_configuration')
                  d = d.isel(hdrlen=[0])
                  data[k]['source_configuration'] = d.to_dataframe()   
                  logging.debug('Done with ', k , ' source_configuration')
                                                     
                  if k == 'era5_1': # reading the whole era5_1 feedback (including reanalysis)
                        d = xr.open_dataset(v , engine = 'h5netcdf' , group = 'era5fb')                 
                        data[k]['era5fb'] = d.to_dataframe()   
                        logging.debug('Done with ', k , ' era5 feedback')
            
                  """ Reading the CDM tables that do not depend on specific stations or observations (fixed values), for the first file only """           
                  if list(datasets.keys()).index(k) == 0  :
                        #for t in ['units' , 'z_coordinate_type' , 'crs' , 'observed_variable']:  
                        for t in [ 'crs' , 'observed_variable']:                              
                              
                              d = xr.open_dataset(v , engine = 'h5netcdf' , group = t)                 
                              data['cdm_tables'][t] = d.to_dataframe()   
                                    
                                   
                  d.close() 
                  ds.close()

                  """ Reading the name of the original source file """
                  source_configuration[k] = {} 
                  source_configuration[k]['source_file'] = [ c for c in v.split('/') if '.nc' in c][0]

                  """ cant fix this right now... Dont know how to extract only the first entries """
                  #d = xr.open_dataset(v , engine = 'h5netcdf' , group = 'source_configuration')                  
                  #data[k]['source_configuration'] = d.to_dataframe()   
                  #print('Done with source conf.')                                    

                                    
            """ Storing the station configurations  """   
            self.source_configuration =  source_configuration      
            
            """ Making all date_times  """   
            self.make_all_datetime()
            
            """ feedback columns """
            if 'era5_1' in list (self.data.keys() ):
                  self.fb_columns = list(self.data['era5_1']['era5fb'].columns ) 
            else:
                  self.fb_columns = ['empty']


      def make_all_datetime(self):
            """ Building the global set of date_times and indices from the various datasets. 
                  The datetimeindex is read from the original netCDF file. 
                  Will compare the unique date_time of each dataset and extract the global unique date_times
                  """
            
            logging.info('\n *** Running make_all_datetime ' )
            
            all_uniques = []  # storing a list with all the unique date_tmes            
            which_k_in_dt = {}  # list of avilable dataset for each unique date_time, so that when looping over the distinct date_times, only the proper dataset will be read and compared 

            def add_time_delta(time_offset_value, date_time):
                  """ Converting to proper date_time adding the time_delta.  
                        Removes minutes rounding to closest integer hour. """ 
                  if 'minutes' in  time_offset:
                        date_time_delta = [ timedelta(minutes = float(i) ) + time_offset_value for i in date_time ]
                  elif 'hours' in time_offset:
                        date_time_delta = [ timedelta(hours = float(i) )  + time_offset_value  for i in date_time ]    
                  else:
                        print('check if time is wrong !!!! (should never happen)')
                        sys.exit()                                               
                  #unique_dt = [i for i in  [  time_offset_value +  j for j in delta  ] ]    
                  #unique_dt = [ i +0 ]
                  date_time_delta = [ i.replace(minute=0, second=0) for i in date_time_delta ]                 
                  return date_time_delta                 


            for k,v in self.datasets.items() :                  
                  self.unique_dates[k] = {}
                  
                  self.unique_dates[k]['indices'] = {}                   
                  #self.unique_dates[k]['indices_low'] = {}     
                  #self.unique_dates[k]['index_up'] = {}                               
                  
                  """ recordtimestamp from the input file """
                  unique = self.data[k]['recordtimestamp']
                  
                  """ Convert to proper date_time using the add_time_delta funtion """
                  logging.debug(' Calculating the time_delta for : ', k )
                  
                  time_offset            = nc.Dataset(self.datasets[k])   
                  time_offset            = time_offset.groups['observations_table']['date_time'].units
                  time_offset_value  = time_offset.split('since ') [1]      
                  time_offset_value  = datetime.strptime(time_offset_value, '%Y-%m-%d %H:%M:%S')
                 
                  unique_dt = add_time_delta (time_offset_value, unique) 
                  
                  all_uniques += unique_dt   # adding to the total unique date_times 
            
                  """ Extracting the recordindex low and up from the input file """
                  indices = self.data[k]['recordindex']
                  
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
                              index_up = len(indices-1)    
                              
                        self.unique_dates[k]['indices'][dt]['up'] = index_up
                             
                  #self.unique_dates[k]['indices'].append(index) 
                  #self.unique_dates[k]['indices_up'].append(index_up) 
                                    
            self.dataset_per_dt = which_k_in_dt             
            self.merged_unique_dates = np.unique(np.array(all_uniques) )  # storing the set of all distinct dt values            
            logging.debug('make_all_datetime finished ')         
      
 
      def clean_dataframe(self, df_in , what = ''):
            """ Remove empty (nan) or wrong values from the original data """         
            
            if what == 'era5fb':  # cleaning the era5 feedback only 
                  df = df_in[np.isfinite(df_in['obsvalue@body'])]
                  df = df.loc[ df['vertco_type@body'] != 2 ]   
                  df = df[np.isfinite(df_in['vertco_reference_1@body'])]
                  #print('check lengths: ' , len(df_in) , len(df) )
                  
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
            
            
      def make_dataframe(self):
            """ Convert netCDF files into panda dataframes. No manipulation of data here; only the CDM columns with real data are included """
            logging.info('*** Creating the dataframes from the source files ' )
            
            for k in self.datasets_keys:
            #for k in ['igra2' , 'ncar']:
            
                  logging.debug('*** Creating the dataframe for the dataset:  ' , k )        
                                    
                  p_levels               = self.data[k]['df']['observations_table']['z_coordinate'][:]
                  logging.debug('     Loaded the  z_coordinate')
                  
                  z_type                 = self.data[k]['df']['observations_table']['z_coordinate_type'][:]
                  logging.debug('     Loaded the  z_coordinate_type')
                  
                  obs_variable        = self.data[k]['df']['observations_table']['observed_variable'][:]
                  logging.debug('     Loaded the  observed_variable')
                  
                  obs_values          = self.data[k]['df']['observations_table']['observation_value'][:]
                  logging.debug('     Loaded the  observation_value')
                  
                  observation_id    = self.data[k]['df']['observations_table']['observation_id'][:]
                  logging.debug('     Loaded the  observation_id')
                  
                  units    = self.data[k]['df']['observations_table']['units'][:].astype(int)
                  logging.debug('     Loaded the  units')                  
                  
                  report_id             = self.data[k]['df']['observations_table']['report_id'][:]                  
                  logging.debug('     Loaded the  report_id')
                  
                  date_time           = self.data[k]['df']['observations_table']['date_time'][:]
                  logging.debug('     Loaded the  date_time (deltas)')
                  
                  lat , lon = self.data[k]['df']['observations_table']['latitude'][:] , self.data[k]['df']['observations_table']['longitude'][:]
                  logging.debug('     Loaded the lat,lon ')
                  
                  
                  self.obs_table_columns = list(self.data[k]['df']['observations_table'].keys() )
                  
                  self.data[k]['df'].close()
 
                  """ Creating a dataframe """
                  columns = ['date_time', 'z_coordinate' , 'z_coordinate_type', 'observed_variable' , 'observation_value' , 'report_id' , 'observation_id' , 'latitude' , 'longitude', 'units']
                  df = pd.DataFrame( list(zip( date_time, p_levels, z_type, obs_variable , obs_values, report_id,  observation_id , lat , lon, units ) ) , columns = columns )       
                  
                        
                  """ Storing the dataframe """      ### try using xarrays ??? 
                  logging.debug('Storing the DF ' )                  
                  self.data[k]['dataframe'] = df
                  
                  logging.debug('  PD dataframe created !!! ')
   

      def get_reanalysis_feedback(self, dt, merged_observations_table, reanalysis = '' , best_ds = ''):
            """ Extracts the renanalysis feedback from the dataset used in the merged file.
                  For now the only available is from era5_1. 
                  Return an empty df otherwise. """
       
            if best_ds == 'era5_1' and reanalysis == 'era5fb':   # reading the feedback for era5_1 
                  index_low = self.unique_dates[best_ds]['indices'][dt]['low']
                  index_up  = self.unique_dates[best_ds]['indices'][dt]['up']
                  
                  fb = self.data[best_ds][reanalysis][index_low:index_up]
                  fb = self.clean_dataframe(fb, what= reanalysis ) # I have to clean the fb exactly the same way I clean the obs_table otherwise they will not match anymore with the indices 
                  
                  merged_observations_table['advanced_assimilation_feedback'] = 1  # filling the flag for the presence of advanced assimilation feedback   
                  
                  return fb , merged_observations_table 
            
            else:  # empty feedback otherwise 
                  """ Creating empty feedback tables. 
                       If era5_1 is never used as best_ds, the era5fb tables are never used and will store an empty, single column DF as feedback. 
                       If era5_1 is used somewehre and its own fb is used, we must use the format of that df (from the odb files) also for the other empty fb """
                  
                  len_odf = len(merged_observations_table)
                  
                  empty_feedback= pd.DataFrame(np.nan, columns= self.fb_columns ,  index = range(len_odf) ) 
                  merged_observations_table['advanced_assimilation_feedback'] = 0              
                  
                  return empty_feedback , merged_observations_table   
                            
                                               
      def get_header_table(self , dt, ds = '' , all_ds = '', length = ''):
            """ Extracting the header_table, and replacing the "duplicates" entries with the list of alternative available datasets """
            index_low = self.unique_dates[ds]['indices'][dt]['low']
            #index_up  = self.unique_dates[best_ds]['indices'][dt]['up']            
            hd = self.data[ds]['header_table'][index_low:index_low+length]                                              
            hd['duplicates'] = all_ds 
            
            return hd


      def merge_all_data(self):       
            """ Construct a dictionary with all the dataframes of each dataset, either reading it from saved pickles or reading it from memory """
            
            logging.info('***** Starting the merging process ')

            
            """ All possible unqiue_dates to loop on """
            date_times = self.merged_unique_dates
            date_times.sort()
               
            date_times = np.array(date_times) 
                                        
            """ List storing the indices of the date_index of the merged dataset """
            all_merged_obs ,  all_merged_head, all_merged_fb , merged_indices , merged_date_time, mi= [] , [] , [] , [] , [], []
           
            
            """ Dictionary that will contain the merged file. """            
            #Merged = {}             
            #for dt in date_times[0:4]: # loop over all the possible date_times 
            
            #chunk = self.data['ncar']['dataframe'] [100:150]

            
            #for dt in date_times[3008:3100]: # loop over all the possible date_times 
            tot = len(date_times)
            for dt, c in zip(date_times[2000:3780] , range(tot)): # loop over all the possible date_times 
            #for dt, c in zip(date_times, range(tot)): # loop over all the possible date_times 
                 
                  print('Analize : ', str(c) , '/',  str(tot)  , ' ', dt , ' ', now(time.time()) )
            
                  cleaned_df_container = {}                  
                  chunk = ''
                  
                  for k in self.dataset_per_dt[dt] :  # checking the list of available datasets  
                                                
                        index, index_up = self.unique_dates[k]['indices'][dt]['low'] , self.unique_dates[k]['indices'][dt]['up']  # extracting the exact chunk of the dataframe where the data of this are stored   
                        
                        chunk = self.data[k]['dataframe'].iloc[index:index_up]
                        
                        #chunk['date_time'].replace( {self.unique_dates[k]['dt_mapping'][dt] : dt} )
                        chunk['date_time'] = dt
                        
                        #self.unique_dates[k]['dt_mapping']  # converting to proper date_time using the self.unique_dates[k]['dt_mapping'] dictionary 

                        chunk = self.clean_dataframe(chunk) # cleaning from wrong values 
 
                        if len(chunk)==0:
                              continue
                        
                        cleaned_df_container[k] = {}                                                
                        cleaned_df_container[k]['df']                   = chunk         # cleaned dataframe 
                        #cleaned_df_container[k]['tot_len']           = len(chunk)   # number of records 

                  
                  if all(value == 0 for value in cleaned_df_container.values()):
                        #print('No data were found! ')
                        continue
                  #print(cleaned_df_container.keys() )
                  merged_observations_table, best_ds, duplicates, header = self.merge_record(dt, container = cleaned_df_container)
                  
                  #if best_ds == 'igra2':
                  #      print('check')

                  merged_observations_table['source_id'] = best_ds   # adding extra columns i.e. chosen dataset, other dataset with data, number of pressure levels 
                  merged_observations_table['z_coordinate_type']  = 1   # only pressure inn [Pa] available at the moment. Check z_coordinate_type table for the correpsonding code 
                                   
                                    
                  """ Extracting the merged feedback, flagging the advanced_observations_feedback flag = 1"""
                  feedback, merged_obs = self.get_reanalysis_feedback( dt, merged_observations_table , reanalysis='era5fb', best_ds= best_ds)
                  all_merged_fb.append(feedback)                  
                  all_merged_obs.append(merged_obs)
                  
                  #""" Extracting the merged header_table """
                  #len_obs = len(merged_observations_table)                  
                  #header = self.get_header_table(dt, best_ds= best_ds,  all_ds = duplicates , length= len_obs)
                  
                  """ Setting the correct report_id in the header table """
                  merged_report_id = merged_obs['report_id'].values[0]  # same report_id as calculated in the observation_table 
                  header['report_id'] = merged_report_id 
                  all_merged_head.append(header)
                                    
                  #if  len(merged_observations_table) !=   len(header):                       
                  #print('lengths check best ds: ', best_ds , '         obs_merged: ' , len(merged_observations_table), '       feedback:' , len(feedback)  , '   header: ' , len(header) )
                  #print( len(merged_observations_table), '       ' , len(feedback)  )

                  """ New merged recordindex and recordtimestamps indices """
                  merged_indices.append(len(merged_observations_table))                 
                  merged_date_time.append(dt)


            """ Storing the merged date_time values and indices """
            di=xr.Dataset()
            merged_date_time = np.array(merged_date_time)
            di['recordtimestamps'] = ( {'recordtimestamps' : merged_date_time.shape } , merged_date_time )
                     
                     
            """ Creating the merged indices """
            mi.append(0)
            for i,ind  in zip(merged_indices[0:], range(len(merged_indices[0:]) ) ) :
                  mi.append(mi[ind] + i  )
            mi = np.array(mi) 
            di['recordindex']          = ( {'recordindex' : mi.shape } , mi )
            self.MergedRecordIndex = di 
                              
                              
            """ Creating the merged dataframes """
            logging.debug('*** Concatenating the observations_table dataframes' )      
            merged_obs = pd.concat (all_merged_obs)
            self.MergedObs = merged_obs                   
            logging.debug('*** Finished concatenating theobservations_table  dataframes' )             
            
            logging.debug('*** Concatenating the header_table dataframes' )      
            merged_hd = pd.concat (all_merged_head)
            self.MergedHead = merged_hd                   
            logging.debug('*** Finished concatenating the header_table dataframes'  )  
            
            logging.debug('*** Concatenating the feedback dataframes' )      
            merged_fb = pd.concat (all_merged_fb)
            self.MergedFeedback = merged_fb                   
            logging.debug('*** Finished concatenating the feedback dataframes' )              

            return 0      
      
      
      def add_cdm_missing_columns(self, all_merged_obs):
            """ Add the CDM observations_table columns for which no data are available at the end of the merging """
            cdm_keys = self.obs_table_columns 
            nan_array = np.empty( all_merged_obs['observed_variable'].shape )
            nan_array[:] = np.nan
            for k in self.obs_table_columns:
                  if k not in list(all_merged_obs.columns ):
                        print('Adding missing cdm colum with empty values: ' , k )
                        all_merged_obs[k] = ( nan_array )
                        
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
                        all_ds_reports.append( self.observation_ids_merged[k] * 1000000000  + container[k]['df']['report_id'].values[0]  )  # converting the original report id using the same convention as for observation_id
                        
            if len(best_datasets) ==0:
                  print('wrong?check')
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
            
            merged_report = self.observation_ids_merged[best_ds] * 1000000000 + selected_df['report_id'].values[0] 
            
            """ Calculate new unique observation id """            
            obs_ids_merged =  [ self.observation_ids_merged[best_ds] * 1000000000 + i for i in selected_df['observation_id'] ]
            selected_df['observation_id'] = obs_ids_merged
            
            """ Calculate new unique report id """            
            selected_df['report_id'] = merged_report

            """ Returning a string with the alternative available datasets data """
            if len(all_ds_reports) > 1: 
                  duplicates =   ",".join( [ str(i) for i in all_ds_reports] )
            else:
                  duplicates = str(all_ds_reports[0])
                  
            
                  
            """ Extracting the merged header_table.
                  Again, must consider the special case where best_ds == ncar. 
                  Note that the header table *should* be identical for ncar_w or ncar_t """            
            if best_ds != 'ncar':
                  header = self.get_header_table(dt, ds= best_ds,  all_ds = duplicates , length= len(selected_df) )
                  
            elif ( best_ds == 'ncar' and  'ncar_t' in list(container.keys()) ) :
                  header = self.get_header_table(dt, ds = 'ncar_t', all_ds = duplicates, length= len(selected_df))
                  
            elif ( best_ds == 'ncar' and 'ncar_t' not in list(container.keys())  ) :
                  header = self.get_header_table(dt, ds = 'ncar_w', all_ds = duplicates, length= len(selected_df) )            
                  
            logging.debug('I use ' , best_ds , '   record since it has more entries: ', most_records , ' but other available datasets are : ' , all_ds ) 
            
            print ('duplicates are: ', duplicates)
            return  selected_df, best_ds , duplicates, header
      

      def write_merged_file(self):
            """ Module to write the output file as netCDF """
            
            out_name = os.getcwd() + '/FAST_INDEX_merged_' + [ x for x in self.datasets[ list(self.datasets_keys)[0]].split('/') if '.nc' in x   ] [0] 
            
            logging.info('Writing the observations_tables to the netCDF output via xarray ')
            #obs_tab = self.MergedObs[ ['date_time' , 'latitude', 'longitude' ,  'observation_value' , 'observed_variable' , 'source_id' , 'observation_id',  'z_coordinate' ]     ] # including only some columns 
            obs_tab = self.MergedObs   # including only some columns             
            obs_tab = self.add_cdm_missing_columns(obs_tab)            
            obs_tab = obs_tab.to_xarray() 
            obs_tab.to_netcdf(out_name, format='netCDF4', engine='h5netcdf', mode='w' , group = 'observations_table')  # writing the merged observations_table 
    
            logging.info('Writing the header_table to the netCDF output via xarray ')
            head_tab = self.MergedHead.to_xarray()
            head_tab.to_netcdf(out_name, format='netCDF4', engine='h5netcdf', mode='a' , group = 'header_table')  # writing the merged observations_table 
            
            logging.info('Writing the station_configuration and source_configurations tables to the netCDF output via xarray ')         
            for k in self.data.keys():
                  if k == 'cdm_tables':
                        continue                  
                  group_name = k + '_station_configuration'
                  sc = self.data[k]['station_configuration'].to_xarray()
                  sc.to_netcdf(out_name, format='netCDF4', engine='h5netcdf', mode='a' , group = group_name )
                  
                  group_name = k + '_source_configuration'
                  sc = self.data[k]['source_configuration'].to_xarray()
                  sc.to_netcdf(out_name, format='netCDF4', engine='h5netcdf', mode='a' , group = group_name )
                  
                  """ To be fixed ! """
                  #group_name = k + '_source_configuration'
                  #sc = self.data[k]['source_configuration'][:1].to_xarray()
                  #sc.to_netcdf(out_name, format='netCDF4', engine='h5netcdf', mode='a' , group = group_name )            
            
            logging.info('Writing the merged record indices to the netCDF output ')      
            di = self.MergedRecordIndex
            di.to_netcdf(out_name, format='netCDF4', engine='h5netcdf', mode='a')
             
            logging.info('Writing the merged feedback to the netCDF output ')      
            group_name = 'era5fb'        
            di = self.MergedFeedback
            di = di.to_xarray()
            di.to_netcdf(out_name, format='netCDF4', engine='h5netcdf', mode='a'  , group = group_name )
                        
            logging.info('Writing the standard cdm tables to the netCDF output ')                  
            for t in self.data['cdm_tables'].keys():                
                  d = self.data['cdm_tables'][t]
                  d = d.to_xarray()
                  d.to_netcdf(out_name, format='netCDF4', engine='h5netcdf', mode='a'  , group = t )
                              
            logging.info('*** Done writing the output netCDF file ')       
            
            
      def merge(self, limit = False , dataframeIsPickled = False):                                           
            """ Call to the merge_all_data() and write_merged_file() methods """                         
            dummy = self.merge_all_data()            
            logging.info('*** Finished merging, now writing the output netCDF file' )         
            a = self.write_merged_file()
            logging.info('*** Done writing output !!! ')




""" 
# keys of the station_configuration table
['primary_id', 'primary_id_scheme', 'record_number', 'secondary_id', 'secondary_id_scheme', 'station_name', 'station_abbreviation', 'alternative_name', 'station_crs', 'longitude', 'latitude',
'local_gravity', 'start_date', 'end_date', 'station_type', 'platform_type', 'platform_sub_type', 'operating_institute', 'operating_territory', 'city', 'contact', 'role', 'observing_frequency',
'reporting_time', 'telecommunication_method', 'station_automation', 'measuring_system_model', 'measuring_system_id', 'observed_variables', 'comment', 'optional_data', 'bbox_min_longitude',
'bbox_max_longitude', 'bbox_min_latitude', 'bbox_max_latitude', 'metadata_contact', 'metadata_contact_role']


['primary_id',  'secondary_id', 'station_name', 'alternative_name', 'station_crs', 'longitude', 'latitude',
'start_date', 'end_date',  'city', 
'observed_variables', ]
"""



""" Example input dictionary """
tateno = {'ncar'    : 'example_stations_big/ncar/chuadb_trhc_47646.txt.nc'   ,
               'igra2'   : 'example_stations_big/igra2/chJAM00047646-data.txt.nc'  ,  
               'bufr'     : 'example_stations_big/bufr/chera5.47646.bfr.nc'  ,
               'era5_1' : 'example_stations_big/era5_1/chera5.conv._47646.nc' , 
               'era5_1759' : 'example_stations_big/era5_1759/chera5.1759.conv.1:47646.nc' , 
               'era5_1761' : 'example_stations_big/era5_1761/chera5.1761.conv.1:47646.nc' , 
               'era5_3188' : 'example_stations_big/era5_3188/chera5.3188.conv.C:5357.nc' , 
               }


small = {   'ncar_w'           : 'example_stations/ncar/chuadb_windc_82930.txt.nc'       ,
                  'ncar_t'            : 'example_stations/ncar/chuadb_trhc_82930.txt.nc'       ,
                  'igra2'              : 'example_stations/igra2/chBRM00082930-data.txt.nc'  ,
                  'era5_1'            :  'example_stations/era5_1/chera5.conv._82930.nc',
                  'era5_1759'      : 'example_stations/era5_1759/chera5.1759.conv.1:82930.nc',
                  'bufr'                : 'example_stations/bufr/chera5.82930.bfr.nc',                          
}



""" main block """
if __name__ == '__main__':
      """ Initialize the Merger class """
      Merging = Merger()
      logging.info('*** Initialising the data ***' )      
      #Merging.initialize_data( datasets = small_other ) #  Read each dataset netCDF file, initialize the dataframes, calculated proper date_time arrays       
      Merging.initialize_data( datasets = small ) #  Read each dataset netCDF file, initialize the dataframes, calculated proper date_time arrays  
      
      """ Merging the data, writing the merged output file """
      Merging.merge(limit = '', dataframeIsPickled = False)
           
      logging.info( '     *** Done  ***     ' ) 





