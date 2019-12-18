""" Merging the station configuration files """

import os,sys
import netCDF4 as nc
import pandas as pd
import numpy as np
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pylab as plt
#import argparse
from datetime import datetime, timedelta
import numpy.ma as ma
import math
import h5py as h5py
import xarray as xr 
import time 
from numba import njit

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning) # deactivates Pandas warnings 


pd.set_option('display.max_rows', 50)
pd.set_option('display.max_columns', 20)
pd.set_option('display.width', 300)


def now(time):
      a = datetime.fromtimestamp( time  ).strftime('%Y-%m-%d %H:%M:%S')
      return  a

"""
For v.0 of the merged dataset

- skip all the data with no pressure information (might have only gepotential; to use these, we will rely on reanalyses data. To Do for next versions)
- igra removed data which did not pass the quality check. If the record is selected from other datasets for some reasons, we keep the data even if removed
- we create a unique record id for the merged dataset in the following way:
     create a nu,bering scheme for each dataset, see dictionary  dataset_id_scheme = {} 
     multiply this number by 1 billion
     add the original dataset record_id to the number
This way, the merged observation_id is unique and it also sotres the original dataset observation_id 
"""
      
class Merger():
      """ Main class for the merging of the data from different netCDF files """
      
      def __init__(self ):
            self.data = {}
            self.datasets = ''
            self.datasets_keys = ''
            self.datasets_all                     = ['igra2' , 'era5_1' , 'ncar' , 'bufr' , 'era5_1759' , 'era5_1761' , 'era5_3188']                          
            self.observation_ids_merged = {  'igra2':1 , 'ncar':2 , 'bufr':3,  'era5_1':4 , 'era5_1759' :5 , 'era5_1761':6 ,  'era5_3188' :7}  # values used to convert original record_id to the merged record_id, see merge_all_data 
            self.unique_dates = {}
            
             
      def initialize_data(self, datasets = {} , fast = False ):
            """ Initialize dataset; store relevant data as attributes.
                       Args ::     dic{}  datasets (dictionary where keys are the dataset names e.g. bufr, igra2 etc. , and the value is the path to the corresponding netCDF file 
                                       e.g. tateno = { 'ncar'    : 'example_stations/ncar/chuadb_trhc_47646.txt.nc'   ,
                                                            'igra2'   : 'example_stations/igra2/chJAM00047646-data.txt.nc'  } """
        
            data = {}
            source_configuration = {}
            
            self.datasets = datasets
            self.datasets_keys = datasets.keys()
            
            """ Looping over the avilable datasets """
            print('*** Reading and Initializing the data from the netCDF files ')
            for k,v in datasets.items() :
                  print('Initialising the dataset: *** ' , k )
                  data[k] = {} 
                  
                  ### xarray                  
                  #ds =  xr.load_dataset(v)   
                  #observations_table =  xr.open_dataset(v , engine = 'h5netcdf' , group = 'observations_table')   
                  #print('xarray')
                  #data[k]['df'] = ds # storing the entire file                
                  #data[k]['dateindex']       = ds['dateindex'][0,:]  # storing the dateindex 
                  #data[k]['source_file']      = ds['source_configuration']['source_file'][0]
                  #data[k]['product_code'] = ds['source_configuration']['product_code'][0]                                       
                  
                  ###for h5py but cant extract date time units !!!
                  #print('Reading the file with h5py ' , now (time.time() ) )
                  ds =  h5py.File(v , driver="core" )   
                  data[k]['df'] = ds # storing the entire file                
                  #data[k]['dateindex']       = ds['dateindex'][0,:]  # storing the dateindex 
                  data[k]['source_file']           = ds['source_configuration']['source_file'][0]
                  data[k]['product_code']       = ds['source_configuration']['product_code'][0]  
                  data[k]['recordtimestamp'] = ds['recordtimestamp'].value
                  data[k]['recordindex']         = ds['recordindex'].value                                    
                  #ds.close()                 
                  print('Reading the file with xarray ' , now (time.time() ) )

            self.data = data
            self.make_dataframe()
            ds.close()                 

            """ Reading the header_table, station_configuration, source_configuration """
            for k,v in datasets.items() :
                  
                  #d = xr.open_dataset(v , engine = 'h5netcdf' , group = 'station_configuration')
                  
                  d = xr.open_dataset(v , engine = 'h5netcdf' , group = 'station_configuration')                 
                  data[k]['station_configuration'] = d.to_dataframe()   
                  #print('Done with ', k , ' station_configuration')
                  
                  d = xr.open_dataset(v , engine = 'h5netcdf' , group = 'header_table')                 
                  data[k]['header_table'] = d.to_dataframe()   
                  #print('Done with ', k , ' header_table')
                  
                  d = xr.open_dataset(v , engine = 'h5netcdf' , group = 'source_configuration')
                  d = d.isel(hdrlen=[0])
                  data[k]['source_configuration'] = d.to_dataframe()   
                  #print('Done with ', k , ' source_configuration')
                  
                  
                  
                  if k == 'era5_1': # reading the whole era5_1 feedback (including reanalysis)
                        d = xr.open_dataset(v , engine = 'h5netcdf' , group = 'era5fb')                 
                        data[k]['era5fb'] = d.to_dataframe()   
                        #print('Done with ', k , ' era5 feedback')
                        
                  d.close() # close ?
                  ds.close()

                  """ Reading the name of the original source file """
                  source_configuration[k] = {} 
                  source_configuration[k]['source_file'] = [ c for c in v.split('/') if '.nc' in c][0]


                  """ cant fix this right now... Dont know how to extract only the first entries """
                  #d = xr.open_dataset(v , engine = 'h5netcdf' , group = 'source_configuration')                  
                  #data[k]['source_configuration'] = d.to_dataframe()   
                  #print('Done with source conf.')
                                    
                  ### for netCDF4
                  #ds =  nc.Dataset(v)                   
                  #data[k]['dateindex'] = ds.variables['dateindex'][0,:]  # storing the dateindex                
                  #data[k]['source_file']      = ds.groups['source_configuration']['source_file'][0]
                  #data[k]['product_code']  = ds.groups['source_configuration']['product_code'][0]
               
            """ Storing the station configurations  """   
            self.source_configuration =  source_configuration      
            
            """ Making all date_times  """   
            self.make_all_datetime()
            
            """ fb columns """
            if 'era5_1' in list (self.data.keys() ):
                  self.fb_columns = list(self.data['era5_1']['era5fb'].columns ) 
            else:
                  self.fb_columns = ['empty']

            

      '''  
      def MakeCDMOutputTables(self):
            """ Store the CDM tables for the output netCDF file """
            print('MakeCDMOutputTables ' , now(time.time()) )
            
            """
            tables = { 'observations_table' : { 'observation_value'   : [] ,
                                                                    'observed_variable'  : [] , 
                                                                    'z_coordinate_type'  : [] ,
                                                                    'date_time'               : [] ,
                                                                    'longitude'                : [] ,
                                                                    'latitude'                   : [] , 
                                                                    'observation_id'        : [] ,     
                                                                    'source_id'                : [] } ,
                       
                             'source_configuration' : {'source_id' : [] ,
                                                                     'source_file' : [] }                      
                       }
            """
            
            source_configuration = {}
            
            for k,v in self.datasets.items():
                  source_configuration[k] = {} 
                  source_configuration[k]['source_file'] = [ c for c in v.split('/') if '.nc' in c][0]
             
            self.source_configuration =  source_configuration
      '''

      def make_all_datetime(self):
            """ Building the global set of date_times and indices from the various datasets. 
                  The datetimeindex is read from the original netCDF file """
            
            print('\n *** Running make_all_datetime ' , now (time.time()) )
            
            all_uniques = []  # storing a list with all the unique date_tmes
            
            """ Storing a list of avilable dataset for each unique date_time, so that when looping over the distinct date_times, only the proper dataset will be read and compared """                                               
            which_k_in_dt = {} 

            def add_time_delta(time_offset_value, date_time):
                  """ Converting to proper date_time adding the time_delta.  """ ### check There should be only hours in the input files !!!
                  if 'minutes' in  time_offset:
                        date_time_delta = [ timedelta(minutes = float(i) ) + time_offset_value for i in date_time ]
                  elif 'hours' in time_offset:
                        date_time_delta = [ timedelta(hours = float(i) )  + time_offset_value  for i in date_time ]    
                        
                  #unique_dt = [i for i in  [  time_offset_value +  j for j in delta  ] ]    
                  #unique_dt = [ i +0 ]
                  return date_time_delta                 


            for k,v in self.datasets.items() :                  
                  self.unique_dates[k] = {}
                  
                  self.unique_dates[k]['indices']      = {}                   
                  self.unique_dates[k]['indices_low']  = {}     # for each dt, I store a the lower and upper index of the dataframe 
                  self.unique_dates[k]['index_up']    = {}                               
                  #self.unique_dates[k]['indices']        = {}  
                  #self.unique_dates[k]['indices_up']  = []
                  
                  """ Extracting the recordtimestamp from the input file """
                  unique = self.data[k]['recordtimestamp']
                  
                  """ Convert to proper date_time usig the time_delta """
                  time_offset            = nc.Dataset(self.datasets[k])   
                  time_offset            = time_offset.groups['observations_table']['date_time'].units
                  time_offset_value  = time_offset.split('since') [1]                            
                  try:
                        time_offset_value  = datetime.strptime(time_offset_value, '%Y-%m-%d %H:%M:%S')
                  except:
                        time_offset_value = time_offset_value.replace(' 19', '19')
                        time_offset_value  = datetime.strptime(time_offset_value, '%Y-%m-%d %H:%M:%S')
                  
                  #print(' Calculating the time_delta for : ', k )
                  
                  unique_dt = add_time_delta (time_offset_value, unique) 
                  
                  #date_time_mapping = dict(zip(unique_dt , unique))
                  #self.unique_dates[k]['dt_mapping'] = date_time_mapping # will be used to replace the values in the dataframe 
                  
                  all_uniques += unique_dt   
            
                  """ Extracting the recordindex low and up from the input file """
                  indices = self.data[k]['recordindex']
                  for dt, low, count in zip (unique_dt,  indices, range(len(unique_dt))  ):
                        
                        try:                          
                              which_k_in_dt[dt].append(k)
                        except:
                              which_k_in_dt[dt] = []
                              which_k_in_dt[dt].append(k)                             
                        
                        self.unique_dates[k]['indices'][dt] = {}
                        self.unique_dates[k]['indices'][dt]['low'] = low                       
                        try:
                              index_up =  indices[ count + 1 ]
                        except:                             
                              index_up = len(indices-1)   # should not be a problem if the index exceeds the length when you slice pandas df  
                              
                        self.unique_dates[k]['indices'][dt]['up'] = index_up
                             
                  #self.unique_dates[k]['indices'].append(index) 
                  #self.unique_dates[k]['indices_up'].append(index_up) 
                  
                  
            self.dataset_per_dt = which_k_in_dt             
            self.merged_unique_dates = np.unique(np.array(all_uniques) )  # storing the set of all distinct dt values            
            #print('make_all_datetime finished ' , now (time.time()) )         
      
      
      
            
            
      '''     
      def MakeDateTime_build(self, date_time = '', time_offset = '' , k = '', shortener = False ):   # only consider a small number of entries           
            """  OLD version, building the date_time from the netCDF , not reading the indices 
            Calculating the actual date_time from the time offset and the time deltas stored in ['observations_table']['date_time'] """
            
            print('MakeDateTime ' , now (time.time()) )
                                 
            unique =   np.unique(date_time) # Getting the unique date_time 
                        
            time_offset_value           = time_offset.split('since ')[1]                            
            time_offset_value           = datetime.strptime(time_offset_value, '%Y-%m-%d %H:%M:%S')
           
            """ Converting to proper date_time """
            if 'minutes' in  time_offset:
                  delta = [ timedelta(minutes = float(i) ) for i in unique ]
            elif 'hours' in time_offset:
                  delta = [ timedelta(hours = float(i) )    for i in unique ]                  
            unique_dt = [i for i in  [  time_offset_value + i  for i in delta  ] ]     
            
            print('MakeDateTime making dic  ' , now (time.time()) )                          
            
            """ Mapping the original date_time list to the new one """
            zipp = zip(unique, unique_dt)
            dic   = dict(zipp)
            
            """ Initialising the containers for the indices """
            self.unique_dates[k] = {}
            self.unique_dates[k]['unique_dates'] = unique_dt 
            self.unique_dates[k]['index'] = {} 
            self.unique_dates[k]['index_up'] = {} 
            self.unique_dates[k]['indices'] = [] 
            self.unique_dates[k]['indices_up'] = []
            
            print('Getting the indices ' , now (time.time()) )                
            for dt in unique:
                  self.unique_dates[k]['index'][dt] = ''
                  matching = np.where( date_time == dt)                       
                  index = min(matching[0] )
                  index_up = max(matching[0] ) + 1
                  
                  self.unique_dates[k]['index'][dic[dt]] = index
                  self.unique_dates[k]['index_up'][dic[dt]] = index_up                         
                  self.unique_dates[k]['indices'].append(index) 
                  self.unique_dates[k]['indices_up'].append(index_up) 
   
            print('MakeDateTime finished ' , now (time.time()) )         
            
            return dic , unique 
        
        
      '''  
      def clean_dataframe(self, df_in , what = ''):
            """ Remove empty or wrong values from a dataframe """         
            
            if what == 'era5fb':  # cleaning the era5 feedback only 
                  df = df_in[np.isfinite(df_in['obsvalue@body'])]
                  df = df.loc[ df['vertco_type@body'] != 2 ]   
                  df = df[np.isfinite(df_in['vertco_reference_1@body'])]
                  #print('check lengths: ' , len(df_in) , len(df) )
                  
            else:
                  
                  df =  df_in.loc[ df_in['z_coordinate_type'] != 2 ]  
                  df = df.loc[ (df['observation_value'] != -99999.0) 
                                       & (df['observation_value'] != -999.0) 
                                       & (df['observation_value'] != -9999)                                        
                                       & (df['observation_value'] != -9999.0) 
                                       & (df['observation_value'] != -999.9) 
                                       #& (df['z_coordinate_type'] != 2)  
                                       & (df['z_coordinate'] != -99999.0) 
                                       & (df['z_coordinate'] != -9999.0) 
                                       
                                       ] #cleaning the values           
                  #clean = clean.loc[ (clean['z_coordinate_type'] != 2)] #cleaning the values
                  #clean = clean.loc[ (clean['z_coordinate'] != -99999.0 )] #cleaning the values
                  df = df[np.isfinite(df['observation_value'])]
                  df = df[np.isfinite(df['z_coordinate'])]
                  
            return df 
            
            
      def make_dataframe(self):
            """ Convert netCDF files into panda dataframes. No manipulation of data here; only the CDM columns with real data are included """
            print('*** Creating the dataframes ' )
            
            for k in self.datasets_keys:
            #for k in ['igra2' , 'ncar']:
            
                  #print('*** Creating the dataframe for the dataset:  ' , k )                  
                  p_levels               = self.data[k]['df']['observations_table']['z_coordinate'][:]
                  #print('     Loading the  z_coordinate')
                  z_type                 = self.data[k]['df']['observations_table']['z_coordinate_type'][:]
                  #print('     Loading the  z_coordinate_type')
                  obs_variable        = self.data[k]['df']['observations_table']['observed_variable'][:]
                  #print('     Loading the  observed_variable')
                  obs_values          = self.data[k]['df']['observations_table']['observation_value'][:]
                  #print('     Loading the  observation_value')
                  observation_id    = self.data[k]['df']['observations_table']['observation_id'][:]
                  #print('     Loading the  observation_id')
                  units    = self.data[k]['df']['observations_table']['units'][:].astype(int)
                  #print('     Loading the  units')                  
                  #report_id             = self.data[k]['df']['header_table']['report_id'][:]
                  report_id             = self.data[k]['df']['observations_table']['report_id'][:]                  
                  #print('     Loading the  report_id')
                  date_time           = self.data[k]['df']['observations_table']['date_time'][:]
                  #print('     Loading the  date_time (deltas)')
                  lat , lon = self.data[k]['df']['observations_table']['latitude'][:] , self.data[k]['df']['observations_table']['longitude'][:]
                  
                  
                  self.obs_table_columns = list(self.data[k]['df']['observations_table'].keys() )
                  
                  self.data[k]['df'].close()
                  
                  '''
                  """ Converting to proper date_time using the time offset i.e. earliest date_time of observation """ 
                  
                  time_offset = nc.Dataset(self.datasets[k])
                  time_offset = time_offset.groups['observations_table']['date_time'].units # cannot use h5py                  
                  dic_dt, unique = self.MakeDateTime(date_time= date_time, time_offset = time_offset  , k = k ) # dictionary mapping odl dt to correct date_time calculated with deltas  
                  
                  self.data[k]['unique_dt']  = unique 
                                 
                  dt_correct = [dic_dt[i] for i in date_time  ]     
                  
                  '''
                  #print('Creating the DF ' , now (time.time() ) )

                  columns = ['date_time', 'z_coordinate' , 'z_coordinate_type', 'observed_variable' , 'observation_value' , 'report_id' , 'observation_id' , 'latitude' , 'longitude', 'units']
                  df = pd.DataFrame( list(zip( date_time, p_levels, z_type, obs_variable , obs_values, report_id,  observation_id , lat , lon, units ) ) , columns = columns )       
                  
                  '''
                  print('Mapping the date_time' , now (time.time() ) )                  
                  df['date_time'] = dt_correct
                  
                  print('Cleaning the DF ' , now (time.time() ) )                  
                  df = df.loc[ (df['observation_value'] > -800 )  & (df['observation_value'] != np.nan)  ]   
                  df = df.loc[ (df['z_coordinate_type'] != 2)] #cleaning the values
                  df = df.loc[ (df['z_coordinate'] != -99999.0 )] #cleaning the values
                  '''

    
                  """ Controllare se serve """
                  #print('Mapping the date_time' , now (time.time() ) )
                  #df = df.replace( {'date_time' : dic_dt } )               
                  #df.sort_values(by = ['date_time' ] )    # should be already sorted 
                
                  """ Converting variables to a specific data type """
                  #variables_types = { 'z_coordinate' : 'float' , 'z_coordinate_type': 'int', 'observation_value':'float' , 'observed_variable':'int' , 'report_id': 'int' , 'observation_id':'int'} 
                  #if 'era5' in k:
                  #      lista = df.index[:]
                  #      df['report_id'] = lista
                  #      df['observation_id'] = lista
                        
                  #for p,v in variables_types.items(): 
                  #      df[p].astype(v) 
                        
                  """ Storing the dataframe """      ### try using xarrays ??? 
                  #print('Storing the DF ' , now (time.time() ) )                  
                  self.data[k]['dataframe'] = df
                  
                  #print('  PD dataframe created !!! ')
   
      '''           
      def MakeAllData(self , pickle = False):            
            """ Creates a global dataframe from the dataframe of each dataset.
                 """
            
            #if pickle:
            #      print('Loading pickle file')
            #      a = np.load('all_data.npy' , allow_pickle = True)
            #      print('Pickle Loaded')
            #      return a 
            
            def vectorize_panda_numpy(dt='', all_data = '', k='' , index = '' , index_up = ''):
                  dataframe = self.data[k]['dataframe'][index:index_up]
                  observed_variable = dataframe['observed_variable'].values
                  press    = dataframe['z_coordinate'].values
                  z_type    = dataframe['z_coordinate_type'].values                  
                  for v,p,z in zip (observed_variable, press, z_type):
                        if p == -99999.0 :
                              continue
                         #if z != 2 :                    # i.e. not equal to geopotential           
                        try:
                              all_data[dt][k][v].append( p )
                        except:
                              all_data[dt][k][v] = []
                              all_data[dt][k][v].append( p )
                        
                        #TODO deal with geopotential 
                        elif z==2:
                              try:
                                    all_data[dt][k][v].append( np.nan )
                              except:
                                    all_data[dt][k][v] = []
                                    all_data[dt][k][v].append( np.nan )                              
                        
      
            all_data = { }
            print('   Making all data')
            for k in self.datasets_keys:
                  print('--- ', k )
                  unique = self.unique_dates[k]['unique_dates'] 
                  indices, indices_up = self.unique_dates[k]['indices'] , self.unique_dates[k]['indices_up']

                  print('Starting the loop over the dataframe at ',  datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') )
                  
                  for dt, indices, indices_up in zip(unique, indices, indices_up):                                            
                        try:
                              all_data[dt][k] = {}                                      
                        except:
                              all_data[dt] = {}                           
                              all_data[dt][k] = {}           
                         
                        """ For each differnet dt I extract the different pressures I have for each different variable v :: all_data[dt][k][v].append(p) """      
                        a = vectorize_panda_numpy(dt=dt, all_data=all_data, k=k, index=indices, index_up=indices_up)

                        #loc = smaller_df.loc[smaller_df['date_time'] == dt]  # extract the items for this particular date_time 
         
                  print('Finished the loop over the dataframe at ',  datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))
            return all_data
      '''
    
      '''
       def MakeAllDateTime_NEW(self): 
             """ Creates a dictionary of unique date_times for each dataset, and a global list of date_times from all the different dataset. 
                   For each date_time, it also stores the index where it is found in each dataset """
             print('MakeAllDateTime ', now (time.time()) )
       
             all_dates = []
             
             for k in self.datasets_keys:                  
                   time_offset = nc.Dataset(self.datasets[k])  # calculating the proper date_time from the netCDF time deltas 
                   time_offset = time_offset.groups['observations_table']['date_time'].units
                   
                   
                   all_dates = all_dates + list( self.unique_dates[k]['unique_dates'] ) 
 
             all_dates = list(set (all_dates) )
             all_dates.sort()
             
             print('Done with the date_times ***')
             
             self.all_dates = all_dates
             
       '''



      def get_reanalysis_feedback(self, dt, merged_observations_table, reanalysis = '' , best_ds = ''):
            """ Extracts the renanalysis feedback from the dataset used in the merged file.
                  For now the only available is from era5_1. 
                  Return an empty df otherwise. """
            
            index_low = self.unique_dates[best_ds]['indices'][dt]['low']
            index_up  = self.unique_dates[best_ds]['indices'][dt]['up']
            
            if best_ds == 'era5_1' and reanalysis == 'era5fb':   # reading the feedback for era5_1 
                  
                  
                  fb = self.data[best_ds][reanalysis][index_low:index_up]
                  fb = self.clean_dataframe(fb, what= reanalysis ) # I have to clean the fb exactly the same way I clean the obs_table otherwise they will not match anymore with the indices 
                  
                  merged_observations_table['advanced_assimilation_feedback'] = 1  # filling the flag for the presence of advanced assimilation feedback                   
                  return fb , merged_observations_table 
            
            else:  # empty feedback otherwise 
                  """ Creating empty feedback tables. 
                       If era5_1 is never used as best_ds, the era5fb tables are never used hence I will store an empty, single column DF as feedback. 
                       If era5_1 is used somewehre and its own fb is used, we must use the format of that df (from the odb files) also for the other empty fb """
                  
                  len_odf = len(merged_observations_table)
                  
                  empty_feedback= pd.DataFrame(np.nan, columns= self.fb_columns ,  index = range(len_odf) ) 
                  merged_observations_table['advanced_assimilation_feedback'] = 0              
                  
                  return empty_feedback , merged_observations_table   
                  
                     
                    
                  """  
                  try:   # case where the era5_1 columns are defined so even the empty feedbacks must maintain the same structure. Example:
                        empty_feedback= pd.DataFrame(np.nan, columns= self.data[k]['era5fb'].columns,  index = range(len_odf) )  # only the first one to copy the structure 
                        #empty_feedback.append([empty_feedback]*len_odf, ignore_index=True)
                        merged_observations_table['advanced_assimilation_feedback'] = 0 
                  
                  except:  # case where the era5_1 columns are not defined and the feedback is everywhere empty
                        #empty_feedback= pd.DataFrame(np.nan, columns= ['empty'] , index=[1] )
                        empty_feedback= pd.DataFrame(np.nan, columns= ['empty'] , index = range(len_odf) )                        
                        #empty_feedback.append([empty_feedback]*len_odf, ignore_index=True)                        
                        merged_observations_table['advanced_assimilation_feedback'] = 0     
                        
                  return empty_feedback , merged_observations_table   


                  """ 
                  
                  
                                   
      def get_header_table(self , dt, best_ds = '' , all_ds = '', length = ''):
            """ Extracting the header_table, and replacing the "duplicates" entries with the list of alternative available datasets """
            index_low = self.unique_dates[best_ds]['indices'][dt]['low']
            #index_up  = self.unique_dates[best_ds]['indices'][dt]['up']            
            hd = self.data[best_ds]['header_table'][index_low:index_low+length]                                              
            hd['duplicates'] = all_ds 
            
            return hd


      def merge_all_data(self):       
            """ Construct a dictionary with all the dataframes of each dataset, either reading it from saved pickles or reading it from memory """
            
            print('***** Starting the merging process ')

            
            """ All possible unqiue_dates to loop on """
            date_times = self.merged_unique_dates
            date_times.sort()
               
            date_times = np.array(date_times) 
                        
    
            #all_var = [38, 85, 106, 107, 117, -7777 ]   # all possible variables (hardcoded) , -7777 is dummy
            
            """ List storing the indices of the date_index of the merged dataset """
            all_merged_obs ,  all_merged_head, all_merged_fb , merged_indices , merged_date_time, mi= [] , [] , [] , [] , [], []
           
            #merged_df = ''   
            
            """ Dictionary that will contain the merged file. """            
            #Merged = {}             
            #for dt in date_times[0:4]: # loop over all the possible date_times 
            
            #chunk = self.data['ncar']['dataframe'] [100:150]

            
            #for dt in date_times[3008:3100]: # loop over all the possible date_times 
            for dt in date_times: # loop over all the possible date_times 
                 
                  print('Analize the date_time number ', str(np.where(date_times == dt)) , ' ' , dt ,  ' ', now(time.time()) )
            
                  cleaned_df_container = {}                  
                  chunk = ''
                  
                  for k in self.dataset_per_dt[dt] :  # checking the list of available datasets  
                                                
                        index, index_up = self.unique_dates[k]['indices'][dt]['low'] , self.unique_dates[k]['indices'][dt]['up']  # extracting the exact chunk of the dataframe where the data of this are stored   
                        
                        chunk = self.data[k]['dataframe'].iloc[index:index_up]
                        
                        #chunk['date_time'].replace( {self.unique_dates[k]['dt_mapping'][dt] : dt} )
                        chunk['date_time'] = dt
                        
                        #self.unique_dates[k]['dt_mapping']  # converting to proper date_time using the self.unique_dates[k]['dt_mapping'] dictionary 

                        chunk = self.clean_dataframe(chunk) # cleaning from wrong values 
                        
                        #print('check')
                        #cleaned_df_container[k] = {}
                        #cleaned_df_container[k]['df']                   = clean
                        #cleaned_df_container[k]['tot_len']           = len(clean)
                        
                        ''' # FAST VERSION
                        clean = small_df
                        clean = small_df.loc[ (small_df['observation_value'] != -7777)  &   (small_df['observation_value'] != -99999.0) & (small_df['observation_value'] != -999.0) & (small_df['observation_value'] != np.nan)  ] #cleaning the values
                        clean = clean.loc[ (clean['z_coordinate_type'] != 2)] #cleaning the values
                        clean = clean.loc[ (clean['z_coordinate'] != -99999.0 )] #cleaning the values
                        
                        if len(clean) == 0:
                              continue 
                        '''
                        if len(chunk)==0:
                              continue
                        
                        cleaned_df_container[k] = {}                                                
                        cleaned_df_container[k]['df']                   = chunk         # cleaned dataframe 
                        cleaned_df_container[k]['tot_len']           = len(chunk)   # number of records 
                        
                        #cleaned_df_container[k]['var']                 = np.unique( clean['observed_variable'].values ) 
                        #cleaned_df_container[k]['p_levels']          = np.unique( clean['z_coordinate'].values ) 
                        #cleaned_df_container[k]['num_var']         = len( cleaned_df_container[k]['var'] )
                        #cleaned_df_container[k]['num_p_levels']  = len( cleaned_df_container[k]['p_levels'] )

                  #if len(cleaned_df_container.keys() ) == 0: # no useful data (not possible anymore I think since I clean the df at the beginning now)
                  #      print('11111111111111111111   No dataset with usable value: I continue ')                        
                  #      continue # next dt 
                  '''
                  merged , best_ds , all_ds = self.AnalizeDecideMergeRecord_new( container = cleaned_df_container )
                  
                  if best_ds == 0:
                        print('2222222222222222   No dataset with usable value: I continue ')                        
                        continue # should never happen 
                  '''
                  
                  if all(value == 0 for value in cleaned_df_container.values()):
                        #print('No data were found! ')
                        continue
                  
                  merged_observations_table, best_ds, duplicates = self.merge_record(container = cleaned_df_container)
                  
                  #if best_ds == 'igra2':
                  #      print('check')

                  merged_observations_table['source_id']        = best_ds   # adding extra columns i.e. chosen dataset, other dataset with data, number of pressure levels 
                  #merged_df['other_dataset'] = all_ds  # TODO dont know where to store this information !!!                
                                   
                  """ Calculate new unique observation id """
                  
                  ''' 
                  merged_id = []                 
                  
                  
                  for ids,d in zip( merged_observations_table['source_id'], merged_observations_table['observation_id' ] ):
                        if d < 999999999:    # HP: not more than 1 billion records for each file ? 
                              merged_id.append( self.observation_ids_merged[ids] * 1000000000 + d ) # converting original observation_id to merged 
                        else: 
                              print('Must increase the unique id !!!')
                              #raise ValueError('Fix the new unique_id !!!!!')                
                              break
                  '''
                                    
                  """ Extracting the merged feedback, flagging the advanced_observations_feedback flag = 1"""
                  feedback, merged_obs = self.get_reanalysis_feedback( dt, merged_observations_table , reanalysis='era5fb', best_ds= best_ds)
                  all_merged_fb.append(feedback)
                  all_merged_obs.append(merged_obs)
                  
                  
                  """ Extracting the merged header_table """
                  len_obs = len(merged_observations_table)                  
                  header = self.get_header_table(dt, best_ds= best_ds,  all_ds = duplicates , length= len_obs)
                  all_merged_head.append(header)
                                    
                  #if  len(merged_observations_table) !=   len(header):                       
                  #print('lengths check best ds: ', best_ds , '         obs_merged: ' , len(merged_observations_table), '       feedback:' , len(feedback)  , '   header: ' , len(header) )
                  #print( len(merged_observations_table), '       ' , len(feedback)  )

                  
                  '''                 
                  try: 
                        Merg = pd.concat( [Merg, Merged] , ignore_index = True )
                  except:
                        Merg = Merged
                  '''
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
            print('*** Concatenating the observations_table dataframes  ' ,  now(time.time()) )      
            merged_obs = pd.concat (all_merged_obs)
            self.MergedObs = merged_obs                   
            print('*** Finished concatenating theobservations_table  dataframes  ' , now(time.time()) )             
            
            print('*** Concatenating the header_table dataframes  ' ,  now(time.time()) )      
            merged_hd = pd.concat (all_merged_head)
            self.MergedHead = merged_hd                   
            print('*** Finished concatenating the header_table dataframes  ' , now(time.time()) )  
            
            print('*** Concatenating the feedback dataframes  ' ,  now(time.time()) )      
            merged_fb = pd.concat (all_merged_fb)
            self.MergedFeedback = merged_fb                   
            print('*** Finished concatenating the feedback dataframes  ' , now(time.time()) )              

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
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
      def merge_record(self, container = ''):
            """ This is the main function that analize each record (i.e. separate ascent) and decides which one to keep as merged.
                 Args ::
                                 summedRecord         , dataframe including all datasets 
                 Return ::
                                 mergedRecord.index , indices of the record to keep
                                 ds                              , name of the dataset to keep
                                 num_record               , number of records (plevels * variables)
                                 all_ds [string]             , other datasets with available data but not selected 
               
               """            
            record_dataset_legth ={}     
            
            for k in container.keys():                 
                  record_dataset_legth[k] = container[k]['tot_len']
                
                
            """ For now, choosing the dataset with more records of all or igra2>ncar>rest  data if available and with same number of records """
            best_ds, all_ds , best_datasets, all_ds_reports = 'dummy' , [] , [], [] # total number of records, name of the chosen dataset , list of other possible dataset with available data                      
            most_records = max( [ container[v]['tot_len'] for v in  container.keys() ] )  # maximum number of records per date_time           
            
            for k, v in record_dataset_legth.items():                 
                  if v == 0:
                        continue
                  if v == most_records:
                        best_datasets.append(k)                                  
                  if v > 0:
                        all_ds.append(k) # all other datasets with smaller number of records than the maximum found 
                        all_ds_reports.append( self.observation_ids_merged[k] * 1000000000  + container[k]['df']['report_id'].values[0]  )  # converting the original report id using the same convention as for observation_id
                        
            if len(best_datasets) ==0:
                  return 0,0,0         
   
            if 'igra2' in best_datasets:
                  best_ds = 'igra2'
            elif 'ncar' in best_datasets:
                  best_ds = 'ncar'
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
            try:     
                  duplicates = duplicates.replace( str(merged_report)+',' , '').replace(str(merged_report), '') 
                  duplicates =   ",".join( [ str(i) for i in all_ds_reports] )
            except:
                  duplicates = ''
                  
            
            #print ('I use ' , best_ds , '   record since it has more entries: ', most_records , ' but other available datasets are : ' , all_ds ) 
            return  selected_df, best_ds , duplicates
      

      def write_merged_file(self):
            """ Module to write the output file as netCDF """
            
            out_name = os.getcwd() + '/FAST_INDEX_merged_' + [ x for x in self.datasets[ list(self.datasets_keys)[0]].split('/') if '.nc' in x   ] [0] 
            
            print('Writing the observations_tables to the netCDF output via xarray ')
            #obs_tab = self.MergedObs[ ['date_time' , 'latitude', 'longitude' ,  'observation_value' , 'observed_variable' , 'source_id' , 'observation_id',  'z_coordinate' ]     ] # including only some columns 
            obs_tab = self.MergedObs   # including only some columns             
            obs_tab = self.add_cdm_missing_columns(obs_tab)            
            obs_tab = obs_tab.to_xarray() 
            obs_tab.to_netcdf(out_name, format='netCDF4', engine='h5netcdf', mode='w' , group = 'observations_table')  # writing the merged observations_table 
    
    
            print('Writing the header_table to the netCDF output via xarray ')
            head_tab = self.MergedHead.to_xarray()
            head_tab.to_netcdf(out_name, format='netCDF4', engine='h5netcdf', mode='a' , group = 'header_table')  # writing the merged observations_table 
            
            print('Writing the station_configuration and source_configurations tables to the netCDF output via xarray ')         
            for k in self.data.keys():
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
            
            print('Writing the merged record indices to the netCDF output ')      
            di = self.MergedRecordIndex
            di.to_netcdf(out_name, format='netCDF4', engine='h5netcdf', mode='a')
             
            print('Writing the merged feedback to the netCDF output ')      
            group_name = 'era5fb'        
            di = self.MergedFeedback
            di = di.to_xarray()
            di.to_netcdf(out_name, format='netCDF4', engine='h5netcdf', mode='a'  , group = group_name )
            
            print('***** Done writing the output netCDF file !')       
            
            
      def merge(self, limit = False , dataframeIsPickled = False):                                           
            #a = self.MakeAllDateTime_NEW()     # creates a dictionary with the unique observation date_time for each dataset , plus a list with all the dates from any datase
            #dictionary_data = self.MakeAllData(pickle = dataframeIsPickled)    
            
            # if dataframeIsPickled: 
            #       dataframe_data = self.Read_Pickle()                             
            dummy = self.merge_all_data()            
            print('Finished merging !!! \n *** NOW: witing the output file' )         
            a = self.write_merged_file()
            print('Done writing output !!! ')
      '''      
      def Read_Pickle(self):
            """ Reads the panda dataframe stored in a pickle file (faster in testing) """            
            print(' Loading the pickled dataframes, storing in self.[d]["dataframe"] *** ')
            for d in self.datasets_keys:
                  print(' Loading -- ' + d )
                  self.data[d]['dataframe']= pd.read_pickle(d+'_dic')

      def Make
      '''


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



tateno = { 'ncar'    : 'example_stations/ncar/chuadb_trhc_47646.txt.nc'   ,
               'igra2'   : 'example_stations/igra2/chJAM00047646-data.txt.nc'  ,
               'bufr'     : 'example_stations/bufr/chera5.47646.bfr.nc'  ,
               # 'era5_1' : 'example_stations/era5_1/chera5.conv._47646.nc' , 
               'era5_1759' : 'example_stations/era5_1759/chera5.1759.conv.1:47646.nc' , 
               'era5_1761' : 'example_stations/era5_1761/chera5.1761.conv.1:47646.nc' , 
               'era5_3188' : 'example_stations/era5_3188/chera5.3188.conv.C:5357.nc' , 
               }



very_small = { 'era5_1759'    : 'example_stations/era5_1759/chera5.1759.conv.1:82930.nc'       ,
                         'bufr'    : 'example_stations/bufr/chera5.82930.bfr.nc'  , }


small = {  'ncar'           : 'example_stations/ncar/chuadb_windc_82930.txt.nc'       ,
                  'igra2'          : 'example_stations/igra2/chBRM00082930-data.txt.nc'  ,
                           'era5_1'       :  'example_stations/era5_1/chera5.conv._82930.nc',
                           'era5_1759' : 'example_stations/era5_1759/chera5.1759.conv.1:82930.nc',
                           'bufr'           : 'example_stations/bufr/chera5.82930.bfr.nc',                          
}


'''
small_other = {  'ncar'           : 'example_stations/ncar/chuadb_windc_82930.txt.nc'       ,

                           'bufr'           : 'example_stations/bufr/chera5.82930.bfr.nc',                          
}
'''

'''
if __name__ == '__main__':

      """ Initialize the Merger class """
      Merging = Merger()
      
      fast  = False # se to True if the dataframe for the station have been stored and can be loaded # TODO check, it might not work anymore 
      
      if fast == False:
            print('*** Initialising the data ***' , now(time.time()) )      
            Merging.InitializeData( datasets = small_other ) #  Read each dataset netCDF file, initialize the dataframes, calculated proper date_time arrays 
            save = Merging.Save_DataFrame()
            print('Dataframes have been saved! ')
            Merging.Merge_new(limit = '', dataframeIsPickled = False)
           
      elif fast == True:
            Merging.InitializeData( datasets = full_data , fast = True )       
            Merging.Merge_new(limit ='', pickled = True) # Merging procedure 
'''



if __name__ == '__main__':
      """ Initialize the Merger class """
      Merging = Merger()
      print('*** Initialising the data ***' , now(time.time()) )      
      #Merging.initialize_data( datasets = small_other ) #  Read each dataset netCDF file, initialize the dataframes, calculated proper date_time arrays       
      Merging.initialize_data( datasets = small ) #  Read each dataset netCDF file, initialize the dataframes, calculated proper date_time arrays       
      Merging.merge(limit = '', dataframeIsPickled = False)
           

            


print('Done ALL ' , now(time.time() ) ) 


