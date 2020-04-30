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

process = psutil.Process(os.getpid())

cend   = '\033[0m'
blue   = '\033[34m'

def now(time):
      """ For easily print a readable current time stamp """
      a = datetime.fromtimestamp(time).strftime('%Y-%m-%d %H:%M:%S')
      return  a

#@njit(cache=True)
def nobsfill(y,x,bds):
      y[:,0]=bds[0]
      y[:,y.shape[1]-x.shape[1]:]=x[:y.shape[0],:]
      return

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
            self.datasets_all = ['era5_2_2']    # all possibly available datasets                          

            self.unique_dates = {}            
            self.attributes = {} # will keep the original attributes from the CDM tables, read from the netCDF files 
            self.id_string_length = 14 # fixed length for record_id and observation_id values 
            self.out_dir = out_dir 
            self.variable_types = {}
            logging.info('*** Initialising the Merging procedure ***' )   
            self.era5b_columns = []  # stores the columns of the era5fb 

      def initialize_data(self , station = '', datasets = {} ):
            """ Initialize dataset; store relevant data as attributes.
                       Args ::     dic{}  datasets (dictionary where keys are the dataset names e.g. bufr, igra2 etc. , and the value is the path to the corresponding netCDF file 
                                       e.g. tateno = { 'ncar'    : 'example_stations/ncar/chuadb_trhc_47646.txt.nc'   ,
                                                               'igra2'   : 'example_stations/igra2/chJAM00047646-data.txt.nc'  } 
            """       
            self.datasets = datasets
            self.datasets_keys = datasets.keys()
            self.station = station
            self.out_name = self.out_dir + '/' + self.station + '_CEUAS_premerged_v0.nc'

            self.observations_table_vars = ['date_time', 'z_coordinate' , 'z_coordinate_type', 'observed_variable' , 'observation_value' , 'report_id' , 'observation_id' , 'latitude' , 'longitude', 'units']

            """ Loading the econding of the tables created from the harvester script and to be applied again """
            self.encodings = np.load('groups_encodings.npy' , allow_pickle = True ).item()
            self.encodings['era5fb'] = np.load('era5fb_encodings.npy' , allow_pickle = True ).item()

            self.dic_type_attributes = np.load('dic_type_attributes.npy',allow_pickle= True).item()
            
            self.obstab_nans_filled = False      

            data = {}   # container for the data of each dataset
            #source_configuration = {}  # container for the source_configuration of each dataset


            """ Looping over the datasets """
            logging.info('*** Reading and Initializing the data from the netCDF files for the station %s ' , station )

            for k,v in datasets.items() :
                  logging.info(' Initialising the dataset: *** %s ' , k  )
                  data[k] = {} 

            data['cdm_tables'] = {}             

            print(blue + 'Memory used before reading data: ', process.memory_info().rss/1000000000 , cend)





            """ Reading the header_table, station_configuration, source_configuration """
            for k,v in datasets.items() :  

                  try:                      
                        d = xr.open_dataset(v , engine = 'h5netcdf' , decode_times = False )    # leaving units as seconds from 1900-01-01 00:00:00 
                        data[k]['recordtimestamp'] = d['recordtimestamp']                  
                        data[k]['recordindex']     = d['recordindex']
                        data[k]['dateindex']       = d['dateindex']

                        logging.debug('Done with %s observations_table' , str(k) )
                  except:
                        print('Failed opening indices: ', v )


                  h5py_file = h5py.File(v, 'r')
                                                                                                                               
                  ####################
                  # OBSERVATIONS TABLE
                  ####################

                  logging.info ('*** Reading the observations_table for %s', k )
                  obs_tab = h5py_file['observations_table']                                                                                                      
                  '''
                  self.variable_types['observations_table'] = {}
                  for i in obs_tab.items():
                        var = i[0]
                        #print(var)
                        tipo = type(obs_tab[var][0]) 
                        if tipo == np.float64 :  # there should be no float64 types                                                 
                              tipo = np.float32
                        self.variable_types['observations_table'][var] = tipo

                  self.variable_types['observations_table']['processing_code'] = np.int32
                  '''

                  data[k]['observations_table']={}
                  for ov in self.observations_table_vars:
                        data[k]['observations_table'][ov] = obs_tab[ov][:]                                                                                                                       

                  logging.debug('Done with %s observations_table with h5py' , str(k) )                                                                                            
                  ###########
                  # ERA5FB
                  ###########

                  logging.info('*** Initializing era5fb for %s', k )
                  era5fb_tab = h5py_file['era5fb']
                  #self.variable_types['era5fb_table'] = {}

                  era5fb_columns = [c for c in era5fb_tab.keys() if 'string' not in c and c!='index' ]
                  self.era5fb_columns = era5fb_columns

                  '''
                  for i in era5fb_tab.items():
                        var = i[0]
                        #print(var)                                                                                                                                                                             tipo = type(era5fb_tab[var][0])
                        if tipo == np.float64 :  # there should be no float64 types                                                                                           
                              tipo = np.float32
                        #self.variable_types['era5fb_tab'][var] = tipo   #### FIX THIS !!!
                  '''
                  data[k]['era5fb_tab']={}
                  for ov in era5fb_columns:
                        data[k]['era5fb_tab'][ov] = era5fb_tab[ov][:]                                                                                                    


                  logging.debug('Done with %s era5fb with h5py' , str(k) ) 

                  #######
                  # HEADER TABLE
                  #######
                  head_tab = h5py_file['header_table']
                  logging.info('*** Loading the header_table')
                  data[k]['header_table'] = {}
                  for var in head_tab.keys():
                        if ('string' in var or 'hdrlen' in var): continue
                        #print('reading var', var , ' ' , len(head_tab[var][:] ) )                        
                        #data[k]['header_table'][var] = np.array(head_tab[var][:])
                        #print('FFFF ', self.dic_type_attributes['header_table'][var]['type'] )
                        data[k]['header_table'][var] = (np.array(head_tab[var][:])).astype(self.dic_type_attributes['header_table'][var]['type'] )

                  #data[k]['header_table'] = (d.to_dataframe().astype({'duplicates':np.dtype('S44')})).to_records(index=False)
                  ###   station_configuration


                  # FFF TODO REMOVE TRY
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
                              df = pd.DataFrame({'source_file': [ source_file ] }).to_xarray() 
                              data[k]['source_configuration'] = df    
                        else:
                              data[k]['source_configuration'] = d                            
                        logging.debug('Done with %s source_configuration' , str(k) )
                        d.close() 
                  except:
                        pass


                  d.close() 


                  """ Reading the CDM tables that do not depend on specific stations or observations (fixed values), for the first file only """           
                  if list(datasets.keys()).index(k) == 0  : # reading the tables from the first dataset file only
                        for t in [ 'crs' , 'observed_variable', 'units' , 'z_coordinate_type' , 'station_type', 'station_configuration_codes']:                                                            
                              d = xr.open_dataset(v , engine = 'h5netcdf' , group = t)                 
                              data['cdm_tables'][t]  =    d                                      
                  d.close()    


            print(blue + 'Memory used after reading data (ons_tables, stat_conf, header_tab): ', process.memory_info().rss/1000000000 , cend)

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
            if 'duplicates' in var_type_dic.keys():
                  var_type_dic['duplicates']=np.bytes_
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
      print(blue + 'Memory used after makind all date_times : ', process.memory_info().rss/1000000000 , cend)



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

      def get_header_table(self , dt, ds = '' ):
            """ Extracting the header_table """
            
            index= np.searchsorted(self.data[ds]['recordtimestamp'].values, dt)
            hd = {}
            for v in self.data[ds]['header_table'].keys():
                  hd[v]    = np.array( [self.data[ds]['header_table'][v][index] ])

            return hd
      '''
      def make_obstab_era5fb_dic(self, dataset = '' , date_time = ''):
            """ Create obs_tab and feedback tables, using only the data for which report_status@hdr==1 """
            index, index_up = self.unique_dates[dataset]['indices'][date_time]['low'] , self.unique_dates[dataset]['indices'][date_time]['up']
            obs_dic , era5fb_dic = {} , {}

            # I select only report_status@hdr==1 in the era5fb, and extract the indices.
            #print(self.data[dataset]['era5fb_tab']['report_status@hdr'][index:index_up])
            indexes_rep_status_1 = np.where(self.data[dataset]['era5fb_tab']['report_status@hdr'][index:index_up] == 1)[0] 
            # I exclude the data which do not have report_status@hdr==1
            for v in self.observations_table_vars:
                  obs_dic[v]    = np.take(self.data[dataset]['observations_table'][v][index:index_up], indexes_rep_status_1)

            for v in self.era5fb_columns:
                  era5fb_dic[v] = np.take(self.data[dataset]['era5fb_tab'][v][index:index_up], indexes_rep_status_1)


            return obs_dic , era5fb_dic
      '''

      def make_obstab_era5fb_dic(self, dataset = '' , date_time = ''):                                                                                                                      
            """ Create obs_tab and feedback tables, using only the data for which report_status@hdr==1 """                                                                                  
            index, index_up = self.unique_dates[dataset]['indices'][date_time]['low'] , self.unique_dates[dataset]['indices'][date_time]['up']                                              
            obs_dic , era5fb_dic = {} , {}                                                                                                                                                  
                                                                                                                                                                                            
            # I exclude the data which do not have report_status@hdr==1                                                                                                                     
            for v in self.observations_table_vars:                                                                                                                                          
                  obs_dic[v]    = self.data[dataset]['observations_table'][v][index:index_up]
                                                                                                                                                                                            
            for v in self.era5fb_columns:                                                                                                                                                   
                  era5fb_dic[v] = self.data[dataset]['era5fb_tab'][v][index:index_up]                                                                                                                                                                                                                    
            return obs_dic , era5fb_dic    


      def merge_all_data(self):       
            """ Construct a dictionary with all the dataframes of each dataset, either reading it from saved pickles or reading it from memory """

            logging.info('***** Starting the merging process ')

            """ All possible unqiue_dates to loop on """
            date_times = self.merged_unique_dates
            date_times.sort()
            date_times = np.array(date_times) 

            """ List storing the indices of the date_index of the merged dataset """
            all_combined_obs ,  all_combined_head, all_combined_era5fb , combined_indices , combined_date_time, combined_indices_out , mi = [] , [] , [] , [] , [], [], [] 

            """ Dictionary that will contain the merged file. """            
            # rand = datetime.strptime('1981-01-03 12:00:00', '%Y-%m-%d %H:%M:%S')  
            #for dt in date_times[3008:3100]: # loop over all the possible date_times 
            dt_bestds_dic = {} # store the selected best dataset for each dt     
            date_times=date_times[20000:22000]
            tot = len(date_times)
            tt=time.time()
            for dt, c in zip(date_times, range(tot) ): # loop over all the possible date_times 
                  if (c+1)%100==0:
                        print('Analize : ', str(c+1) , '/',  str(tot)  , ' ', dt , ' ',
                              now(time.time()),'{:5.3f}'.format(time.time()-tt ))

                  cleaned_df_container = {}                  

                  for k in self.dataset_per_dt[dt] :  # checking the list of available datasets  

                        obs_tab, era5fb_tab = self.make_obstab_era5fb_dic(dataset = k , date_time = dt)

                        if len(obs_tab['date_time'][:])==0: 
                              #print('skipping empty observations')
                              continue

                        cleaned_df_container[k] = {}                                                
                        cleaned_df_container[k]['obs_tab']    = obs_tab         # cleaned dataframe 
                        cleaned_df_container[k]['era5fb_tab'] = era5fb_tab       # cleaned dataframe                                                                                                                     
                  #if len(cleaned_df_container) == 0:
                  #      #print('Skipping empty container FFFF ')
                  #      continue 

                  best_ds, combined_obs_tab, combined_era5fb_tab = self.combine_record(dt, container = cleaned_df_container)

                  if 'advanced_assimilation_feedback' not in combined_obs_tab.keys() and 'era5' in best_ds:
                        combined_obs_tab['advanced_assimilation_feedback'] = np.array([1]*len(combined_obs_tab['date_time']) )
                        #code.interact(local = locals())

                  """ Storing the combined era5fb """
                  all_combined_era5fb.append(combined_era5fb_tab)

                  """ Extracting and storing the header """
                  combined_head = self.get_header_table(dt, ds = best_ds )
                  all_combined_head.append(combined_head)


                  dt_bestds_dic[dt] = {}
                  dt_bestds_dic[dt]['best_ds'] = best_ds
                  dt_bestds_dic[dt]['len'] = len(combined_obs_tab['date_time'])

                  all_combined_obs.append(combined_obs_tab)

                  """ New merged recordindex and recordtimestamps indices """
                  combined_indices.append(len(combined_obs_tab['date_time']))                 
                  combined_date_time.append(dt)

                  del cleaned_df_container 

            #print(blue + 'Memory used after deleting the cleaned_df_container: ', process.memory_info().rss/1000000000 , cend)


            """ Storing the merged date_time values and indices """
            di=xr.Dataset()
            combined_date_time = np.array(combined_date_time)
            di['recordtimestamp'] = ( {'recordtimestamp' : combined_date_time.shape } , combined_date_time )
            di['recordtimestamp'].attrs['units']='seconds since 1900-01-01 00:00:00'


            """ Creating the merged indices """
            mi.append(0)
            for i in range(len(combined_indices)):
                  mi.append( combined_indices[i] + mi[-1] )
            mi = np.array(mi) 
            di['recordindex']          = ( {'recordindex' : mi.shape } , mi )


            """ Creating the combined data """
            logging.debug('*** Concatenating the observations_table ' )      
            combined_obs = {}
            ####  Writing combined observations_table dic
            for k in all_combined_obs[0].keys():                        
                  combined_obs[k]=np.concatenate([all_combined_obs[i][k][:] for i in range(len(all_combined_obs))])
                  self.tot_records = len(combined_obs[k])

                  self.write_merged(content = 'observations_table', table= {k:combined_obs[k]})
                  if k not in ['date_time']:
                        del combined_obs[k]

            #self.tot_records = len(combined_obs['date_time'])
            del all_combined_obs
            print(blue + 'Memory used after deleting all_merged obs_tab dic: ', process.memory_info().rss/1000000000 , cend )
 
            dateindex = combined_obs['date_time']//86400                                                                                                                                  
            date_times, indices, counts = np.unique(dateindex, return_counts = True, return_index= True)                                                                                  
            di['dateindex'] = ( {'dateindex' : indices.shape } , indices )  # considers the day only                                                                                      
            self.write_merged(content = 'recordindex', table = di)                                                                                                                        
            del di , combined_obs
            
            combined_era5fb = {}
            ####  Writing combined era5fb_table dic                                                                                                                      
            for k in all_combined_era5fb[0].keys():
                  combined_era5fb[k]=np.concatenate([all_combined_era5fb[i][k][:] for i in range(len(all_combined_era5fb))])
                  self.write_merged(content = 'era5fb', table= {k:combined_era5fb[k]})

            del all_combined_era5fb
            print(blue + 'Memory used after deleting era5fb_tab dic: ', process.memory_info().rss/1000000000 , cend)

            combined_head = {}
            ####  Writing combined header_table dic                                                                                                                   
   
            #code.interact(local = locals()) 
            for k in all_combined_head[0].keys():
                  print('head variable is', k )
                  if ( k == 'comments' or k == 'history'):
                        continue
                  try:
                        tab=np.concatenate([all_combined_head[i][k][:] for i in range(len(all_combined_head))])
                        self.write_merged(content = 'header_table', table= {k: tab})
                  except:
                        print('FFF FAILED variable in header table', k )


            del all_combined_head
            print(blue + 'Memory used after deleting all_merged head_tab dic: ', process.memory_info().rss/1000000000 , cend)
            



            '''
            dateindex = combined_obs['date_time']//86400           
            date_times, indices, counts = np.unique(dateindex, return_counts = True, return_index= True)
            di['dateindex'] = ( {'dateindex' : indices.shape } , indices )  # considers the day only
            self.write_merged(content = 'recordindex', table = di)                                                                                      
            del di , combined_obs
            '''

            #logging.debug('*** Concatenating the header_table dataframes' )      
            #merged_hd = {}

            #x=np.concatenate([all_combined_head])
            #l=0
            #for k in self.data['igra2']['header_table'].dtype.names:
            #      combined_hd[k]=np.array([x[i][k] for i in range(x.shape[0])])# np.concatenate([all_merged_head])
            #      l+=1
            #for k in all_merged_head[0].keys():
                  #merged_hd[k]=np.concatenate([all_merged_head[i][k][:] for i in range(len(all_merged_head))])
#            merged_hd = pd.concat (all_merged_head)
            #merged_hd = merged_hd.replace( -2147483648 , np.nan )           
            #self.write_merged(content = 'header_table', table= combined_hd)                          
            #logging.debug('*** Finished concatenating the header_table dataframes'  )  
            #del combined_hd , all_combined_head,x


            return 0      


      def add_cdm_missing_columns(self, all_merged_obs , table = 'observations_table'):
            """ Add the CDM observations_table columns for which no data are available at the end of the merging """
            #cdm_keys = self.obs_table_columns 
            print(self.encodings['observations_table'].keys() )
            for k in self.encodings['observations_table'].keys():
                  if self.variable_types[table][k] == np.int32 :
                        nan = np.int32(-2147483648)
                  else:
                        nan = np.float32(np.nan) 

                  if k not in list(all_merged_obs.keys() ):
                        logging.debug('Adding missing cdm colum with empty values: %s' , k )
                        if k=='sensor_id':

                              print(k)
                        all_merged_obs[k] = np.empty_like(all_merged_obs['date_time'],dtype=np.dtype(nan))
                        all_merged_obs[k].fill(nan)

            return all_merged_obs


      def combine_record(self, dt, container = ''):
            """ This is the main function that analize each record (i.e. separate ascent) and decides which one to keep as merged. """            
            record_dataset_legth ={}     

            for k in container.keys():
                  record_dataset_legth[k] = len(container[k]['obs_tab'] )

            if len(container) == 1:
                  return list(container.keys())[0], container[k]['obs_tab'] , container[k]['era5fb_tab']

            """ For now, choosing the dataset with more records """
            best_ds, best_datasets, all_ds = 'dummy' , [] , [] # total number of records, name of the chosen dataset , list of other possible dataset with available data  

            most_records = max( [ v for v in  record_dataset_legth.values() ] )  # maximum number of records per date_time           

            for k, v in record_dataset_legth.items():                 
                  if v == 0:
                        continue
                  if v == most_records:
                        best_datasets.append(k)                                  
                  if v > 0:
                        all_ds.append(k) # all other datasets with smaller number of records than the maximum found              

            if len(best_datasets) ==0:
                  print('wrong??? please check')
                  return 0,0,0,0        

            elif 'era5_2_1' in best_datasets:
                  best_ds = 'era5_2_1'                  
            else:
                  best_ds = 'era5_2_2'


            selected_obstab, selected_era5fb = container[best_ds]['obs_tab'] , container[best_ds]['era5fb_tab']

            selected_obstab['advanced_assimilation_feedback'] = np.array([1]*len(selected_obstab['date_time']) )

            #code.interact(local = locals()) 
            return  best_ds, selected_obstab, selected_era5fb


      def write_merged(self, content = '', table=''):
            """ Module to write the output file as netCDF """

            if not os.path.isdir(self.out_dir):
                  Path(self.out_dir).mkdir(parents=True, exist_ok=True)                  
            out_name = self.out_dir + '/' + self.station + '_CEUAS_merged_v0.nc'  

            '''
            if os.path.isfile('dic_obstab_attributes.npy'):
                  attrs_dic = np.load('dic_obstab_attributes.npy' , allow_pickle = True).item()
            else:
                  attrs_dic = {}
            '''
            attrs_dic = {}

            if content in ['observations_table','header_table','era5fb']:
                  for var in self.dic_type_attributes[content].keys():
                        if var == 'comments':
                              continue 

                        attrs_dic[var] = {}
                        try:
                              attrs_dic[var]['description']    = bytes( self.dic_type_attributes[content][var]['description']    , 'utf-8' )
                        except:
                              attrs_dic[var]['description']    = bytes( 'missing'    , 'utf-8' )
                              print(' FFF FAILING WITH DESCRIPTION: ', var ) # FFF CHECK WHY SOME ARE FAILING

                        try:
                              attrs_dic[var]['external_table'] = bytes( self.dic_type_attributes[content][var]['external_table'] , 'utf-8' )
                        except:
                              attrs_dic[var]['external_table'] = bytes( 'missing' , 'utf-8' )
                              print(' FFF FAILING WITH EXTERNAL TABLE : ', var ) # FFF CHECK WHY SOME ARE FAILING                                                          


            if content == 'recordindex':  # writing the recordindex, recordtimestamp, dateindex
                  logging.info('Writing the merged record indices to the netCDF output ')
                  table.to_netcdf(out_name, format='netCDF4', engine='h5netcdf', mode='a')

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

            # Writing the observations_table, header_table, era5fb 
            elif content == 'observations_table' or content =='era5fb' or content == 'header_table': 
                  logging.info('Writing the observations_tables to the netCDF output ') 

                  shape = ''
                  for k in table.keys(): 

                        if k == 'index' or k == 'hdrlen' or 'string' in k :
                              continue

                        var_type = self.dic_type_attributes[content][k]['type']

                        ''' trying to convert the variable types to the correct types stored as attribute, read from the numpy dic file '''
                        if type(table[k][0]) != var_type:

                              if k == 'hdrlen': continue
                              #print ('Converting column' , k , '   type  ', type(table[k][0]) , ' to type   ' ,  var_type  )
                              #try: ### FIX REPORT_ID
                              #      if k == 'report_id':
                              #            table[k] =  table[k].astype(self.variable_types[content][k] ) 
                              #except:
                              #      pass
                              try:
                                    table[k] = table[k].astype( var_type ) 
                              except:
                                    print ('FAILED converting column ' , k, ' type ', type(table[v][0]) , ' to type ', var_type )

                        print('*** Writing the table ', content, ' variable ',  k)

                        dic = {k:table[k]}  # making a 1 colum dictionary
                        #shape = table[k].shape
                        #print('SHAPE IS FFF ', table[k].shape )
                        write_dict_h5(out_name, dic , content, self.encodings[content], var_selection=[], mode='a', attrs = attrs_dic  )


                  if content == 'observations_table' and self.obstab_nans_filled == False:
                        missing_cdm_var = [ v for v in self.dic_type_attributes[content].keys() if v not in self.observations_table_vars]  # variables to be filled with nans            
                        for k in missing_cdm_var:
                              if k == 'advanced_assimilation_feedback':
                                    continue
                              var_type = self.dic_type_attributes[content][k]['type']
                              if var_type == np.int32 :
                                    nan = np.int32(-2147483648)
                              else:
                                    nan = np.float32(np.nan) 
      
                              print('fff Adding missing cdm colum with empty values: %s' , k, ' ' , shape )
                              
                              dic={k:np.empty(shape,dtype=np.dtype(nan))}

                              dic[k].fill(nan)

                              write_dict_h5(out_name, dic, 'observations_table', self.encodings['observations_table'],var_selection=[], mode='a', attrs = attrs_dic  ) ### TO DO
                        self.obstab_nans_filled = True

                  elif content == 'observations_table' and self.obstab_nans_filled == True:
                        return

            # Writing the header_table
            #elif content == 'header_table':
            #      logging.info('Writing the header_table to the netCDF output via xarray ')
            #      head_tab = table
            #      for v in head_tab.keys():   
            #            if v == "index" or v == "hdrlen" or v == "string80":
            #                  continue            

            #            if not type(head_tab[v][0]) == self.variable_types['header_table'][v]:
            #                  try:                              
            #                        head_tab[v] =  head_tab[v].astype(self.variable_types['header_table'][v] ) 
            #                  except:
            #                        print ('FAILED converting column' , v , '   type  ', type(head_tab[v][0]) , ' to type   ' ,  self.variable_types['header_table'][v]  )                              
            #            if 'time' not in v and type(head_tab[v][0]) == np.int64 :
            #                  head_tab[v] = head_tab[v].astype(np.int32)
            #            if type(head_tab[v][0]) == np.float64:
            #                  head_tab[v] = head_tab[v].astype(np.float32)
            #            a = {v:head_tab[v]}  # making a 1 column dataframe  
            #            try: 
            #                  write_dict_h5(out_name, a, 'header_table', self.encodings['header_table'], var_selection=[], mode='a', attrs = attrs_dic  ) ### TO DO !!!             
            #            except:
            #                  print('ERROR with attributes codec')
            #                  write_dict_h5(out_name, a, 'header_table', self.encodings['header_table'], var_selection=[], mode='a', attrs = {} )






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
                  except MemoryError:
                        print('Failed: ' , station )
                        o.write(station + '\n' )
                        return False 






f1 = '/raid60/scratch/federico/BISMARK/0-20000-0-72764_era5_2_harvested_era5.conv._1:72764.gz.nc'
f2 = '/raid60/scratch/federico/BISMARK/0-20000-0-72764_era5_2_harvested_era5.conv._72764.gz.nc'

out_dir = '/raid60/scratch/federico/BISMARK/combined/'

os.system('rm /raid60/scratch/federico/BISMARK/combined/0-20000-0-72764_CEUAS_merged_v0.nc')

if not os.path.isdir(out_dir):
      os.mkdir(out_dir)



""" main block """
if __name__ == '__main__':

      Merging = Merger(out_dir)

      stat = '0-20000-0-72764'
      run_mode = 'dummy'
      stat_dic  = {'era5_2_1': f1, 'era5_2_2': f2 }

      a = Merging.merge(stat,  datasets = stat_dic , mode = run_mode ) # single station dictionary

      logging.info('*** Convertion completed! ***')


