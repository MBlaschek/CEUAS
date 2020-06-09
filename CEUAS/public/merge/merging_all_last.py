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
            #self.datasets_all = ['era5_2_2']    # all possibly available datasets                          

            self.unique_dates = {}            
            self.attributes = {} # will keep the original attributes from the CDM tables, read from the netCDF files 
            self.id_string_length = 14 # fixed length for record_id and observation_id values 
            self.out_dir = out_dir 
            self.variable_types = {}
            self.observation_ids_merged  = {  'igra2':b'3' , 'ncar':b'4', 'bufr':b'5',  'era5_1':b'1' , 'era5_2':b'2', 'era5_1759' :b'6' , 'era5_1761':b'7' ,  'era5_3188' :b'8' }  # values used to convert original record_id to the merged record_id, see method merge_all_data 

            logging.info('*** Initialising the Merging procedure ***' )   
            #self.era5b_columns = []  # stores the columns of the era5fb 
            self.standard_cdm = [ 'crs' , 'observed_variable', 'units' , 'z_coordinate_type' , 'station_type', 'station_configuration_codes'] 

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

            self.observations_table_vars = ['date_time', 'z_coordinate' , 'z_coordinate_type', 'observed_variable' , 'observation_value' , 'report_id' , 'observation_id' , 'latitude' , 'longitude', 'units', 'source_id']

            """ Loading the econding of the tables created from the harvester script and to be applied again """
            self.encodings = np.load('groups_encodings.npy' , allow_pickle = True ).item()
            #self.encodings['era5fb'] = np.load('era5fb_encodings.npy' , allow_pickle = True ).item() # old incomplete
            self.encodings['era5fb'] = np.load('era5fb_encodings_all.npy' , allow_pickle = True ).item()            
            self.dic_type_attributes = np.load('dic_type_attributes.npy',allow_pickle= True).item()
            
            self.era5fb_columns = self.dic_type_attributes['era5fb'].keys()

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


            """ Loop over all the datasets                                                                                                                                     
                k: name of the dataset                                                                                                                    
                v: list of file paths,    eg 'era5_1':[filepath_1, filepath_2 ] """

            for k,v in self.datasets.items() :
                  data[k] = {}
                  for F in v:
                        data[k][F] = {}

                        h5py_file = h5py.File(F, 'r')
             
                        
                        data[k][F]['recordtimestamp'] = h5py_file['recordtimestamp']
                        data[k][F]['recordindex']         = h5py_file['recordindex']
                        data[k][F]['dateindex']            = h5py_file['dateindex']
                        
                                                                                                                  
                        ####################
                        # OBSERVATIONS TABLE
                        ####################

                        logging.info ('*** Reading the observations_table for %s', k )
                        obs_tab = h5py_file['observations_table']                                                                                                      


                        data[k][F]['observations_table']={}
                        for ov in self.observations_table_vars:
                              data[k][F]['observations_table'][ov] = obs_tab[ov][:]

                        logging.debug('Done with %s observations_table with h5py' , str(k) )                                                                                            

                        ###########
                        # ERA5FB
                        ###########
                        if k == 'era5_1' or k == 'era5_2':
                              logging.info('*** Initializing era5fb for %s', k )
                              era5fb_tab = h5py_file['era5fb']

                              #era5fb_columns = [c for c in era5fb_tab.keys() if 'string' not in c and c!='index' and c != 'source_id']
                              #self.era5fb_columns = era5fb_columns

                              data[k][F]['era5fb_tab']={}

                              for ov in self.era5fb_columns:
                                    try:
                                          data[k][F]['era5fb_tab'][ov] = era5fb_tab[ov][:]
                                    except:
                                          print("CANNOT FIND  ", ov ) 
                              logging.debug('Done with %s era5fb with h5py' , str(k) ) 

                        #######
                        # HEADER TABLE
                        #######
                        head_tab = h5py_file['header_table']
                        logging.info('*** Loading the header_table')
                        data[k][F]['header_table'] = {}
                        for var in head_tab.keys():
                              if ('string' in var or 'hdrlen' in var): continue
                              try:                                  
                                    data[k][F]['header_table'][var] = (np.array(head_tab[var][:])).astype(self.dic_type_attributes['header_table'][var]['type'] )
                              except:
                                    print('failed convertion type header' , k , ' ' , F , ' ' ,  var )
                        #######                                           
                        # STATION CONFIGURATION
                        #######   
                        d = xr.open_dataset(F , engine = 'h5netcdf' , group = 'station_configuration' , decode_times = False )
                        data[k][F]['station_configuration'] = d
                        logging.debug('Done with %s station_configuration' , str(k) )
                        d.close()

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

            #code.interact(local = locals())
            print(blue + 'Memory used after reading data: ', process.memory_info().rss/1000000000 , cend)

            self.data = data

            """ Making all date_times  """
            self.make_all_datetime()





                  # FFF TODO REMOVE TRY
                  #try:
                  #      d = xr.open_dataset(v , engine = 'h5netcdf' , group = 'station_configuration' , decode_times = False )                 
                  #      #data[k]['station_configuration'] = d.to_dataframe()  
                  #      data[k]['station_configuration'] = d                    
                  #      logging.debug('Done with %s station_configuration' , str(k) )
                  #      d.close()                   
                  #except:
                  #      pass

                  ###   source_configuration
                  #try:                        
                  #      d = xr.open_dataset(v , engine = 'h5netcdf' , group = 'source_configuration' , decode_times = False )                 
                  #      if len(d) == 0:
                  #            source_file = v.split('/')[-1]
                  #            df = pd.DataFrame({'source_file': [ source_file ] }).to_xarray() 
                  #            data[k]['source_configuration'] = df    
                  #      else:
                  #            data[k]['source_configuration'] = d                            
                  #      logging.debug('Done with %s source_configuration' , str(k) )
                  #      d.close() 
                  #except:
                  #       pass

                  #d.close() 
                  
                  

                  #""" Reading the CDM tables that do not depend on specific stations or observations (fixed values), for the first file only """           
                  #for t in self.standard_cdm: # [ 'crs' , 'observed_variable', 'units' , 'z_coordinate_type' , 'station_type', 'station_configuration_codes']
                  #      if t not in data.keys():
                  #            data[t] = {}
                  #            for v in h5py_file[t].keys():
                  #                  if 'string' in v:
                  #                        continue
                  #                  data[t][v] = np.array(h5py_file[t][v])

            #print(blue + 'Memory used after reading data (ons_tables, stat_conf, header_tab): ', process.memory_info().rss/1000000000 , cend)
            #self.data = data     
            #""" Making all date_times  """   
            #self.make_all_datetime()            


      def make_all_datetime(self):
            """ Building the global set of date_times and indices from the various datasets. 
                  The datetimeindex is read from the original netCDF file. 
                  Will compare the unique date_time of each dataset and extract the global unique date_times
                  """

            logging.info('\n *** Running make_all_datetime ' )

            all_uniques = []  # storing a list with all the unique date_times            
            which_k_in_dt = {}  # list of avilable dataset for each unique date_time, so that when looping over the distinct date_times, only the proper dataset will be read and compared 

            """ Loop over all the datasets 
                k: name of the dataset
                v: list of file paths,    eg 'era5_1':[filepath_1, filepath_2 ]"""

            for k,v in self.datasets.items() :
                  self.unique_dates[k] = {}
                  for F in v: 
                        self.unique_dates[k][F] = {}
                  
                        self.unique_dates[k][F]['indices'] = {}                             

                        unique_dt = list(self.data[k][F]['recordtimestamp'])
                  
                        indices   = list(self.data[k][F]['recordindex'])
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
                                    index_up = max(indices)+10000    

                              self.unique_dates[k][F]['indices'][dt]['up'] = index_up


            self.dataset_per_dt = which_k_in_dt             
            self.merged_unique_dates = np.unique(np.array(all_uniques) )  # storing the set of *ALL* distinct dt values of all datasets and all files            
            logging.debug('*** make_all_datetime finished ')         
      print(blue + 'Memory used after makind all date_times : ', process.memory_info().rss/1000000000 , cend)


      '''
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
      '''


      def get_header_table(self , dt, ds = '' , File = ''):
            """ Extracting the header_table """
            
            index= np.searchsorted(self.data[ds][File]['recordtimestamp'], dt)
            hd = {}
            for v in self.data[ds][File]['header_table'].keys():
                  hd[v]    = np.array( [self.data[ds][File]['header_table'][v][index] ])

            return hd


      def make_obstab_era5fb_dic(self, dataset = '' , date_time = '', File = ''):                                                                                                                      
            """ Create obs_tab and feedback tables """                                                                                  
            index, index_up = self.unique_dates[dataset][File]['indices'][date_time]['low'] , self.unique_dates[dataset][File]['indices'][date_time]['up']                                             
            # I exclude the data which do not have report_status@hdr==1 NO: this approach was abandoned                                                                                 
            obs_dic = {}                                    
            for v in self.observations_table_vars:   
                  obs_dic[v]    = self.data[dataset][File]['observations_table'][v][index:index_up]

            """ Loop over the obs_tab to find duplicates.
                I fill a dictionary for each distinct pressure level, and I put inside
                the observed_variable number.
                If the list lready contains the combination pressure level - observed variable,
                then the record is skipped """

            indices = [] # these are the only non-duplicates to be kept

            already_selected = { }
            
            for p,var,val,ind in zip ( obs_dic['z_coordinate'] , obs_dic['observed_variable'],obs_dic['observation_value'] ,range(len(obs_dic['z_coordinate'])) ):
                  if p not in already_selected.keys():
                        already_selected[p] = []
                  
                  if np.isfinite(val):
                        if var not in already_selected[p]:
                              already_selected[p].append(var)
                              indices.append(ind) # record to be kept
                        else:
                              pass
                  else: # skipping nans
                        pass

            red_obs_dic = {} # dictionary for the reduced (removed duplicates) obs_tab
            for v in self.observations_table_vars:
                  red_obs_dic[v] = obs_dic[v][indices]

            ''' Simply returns the proper format for ''null' value '''
            def get_null( tipo = ''):
                  if tipo == np.int32 :
                        void = 0
                  elif tipo == np.float32 :
                        void = 0.0
                  elif tipo == np.bytes_ :
                        void = b'nan'
                  return void
                  
            ''' Filling the feedback table. Only feednack for era5_1 and era5_2 are currently available. 
                Reads the total number of possible columns from the dic_type_attributes dictionary.
                Era5_1 and era5_2 fb have different columns.
                If data for a variable is not available, it fills with the appropriate null value '''
            
            red_era5fb_dic = {}
            for v in self.era5fb_columns:
                  tipo = self.dic_type_attributes['era5fb'][v]['type']                   
                  if dataset == 'era5_1' or dataset == 'era5_2':
                        if v in self.data[dataset][File]['era5fb_tab'].keys():                            
                              red_era5fb_dic[v] = self.data[dataset][File]['era5fb_tab'][v][index:index_up][indices]
                        else:
                              void = get_null(tipo = tipo)
                              red_era5fb_dic[v]= np.full(len(indices), void)                              
                  else:       # no feedback for non era%-1 or era5_2 datasets 
                        void = get_null(tipo = tipo)
                        red_era5fb_dic[v]= np.full(len(indices), void)
                          
            return red_obs_dic , red_era5fb_dic    


      def merge_all_data(self):       
            """ Construct a dictionary with all the dataframes of each dataset, either reading it from saved pickles or reading it from memory """

            logging.info('***** Starting the merging process ')

            """ All possible unique_dates to loop on """
            date_times = self.merged_unique_dates
            date_times.sort()
            date_times = np.array(date_times) 

            """ List storing the indices of the date_index of the merged dataset """
            all_combined_obs ,  all_combined_head, all_combined_era5fb , combined_indices , combined_date_time, combined_indices_out , mi = [] , [] , [] , [] , [], [], [] 
            source_files = []
            
            """ Dictionary that will contain the merged file. """            
            # rand = datetime.strptime('1981-01-03 12:00:00', '%Y-%m-%d %H:%M:%S')  
            dt_bestds_dic = {} # store the selected best dataset for each dt     
            #date_times=date_times[0:2000]
            tot = len(date_times)
            tt=time.time()
            for dt, c in zip(date_times, range(tot) ): # loop over all the possible date_times 
                  if (c+1)%100==0:
                        print('Analize : ', str(c+1) , '/',  str(tot)  , ' ', dt , ' ',
                              now(time.time()),'{:5.3f}'.format(time.time()-tt ))

                  cleaned_df_container = {}                  

                  for k in self.dataset_per_dt[dt].keys() :  # checking the list of available datasets  
                        ''' {'era5_2': ['example_stations/0-20000-0-82930_era5_2_harvested_era5.conv._1:82930.gz.nc', 
                                        'example_stations/0-20000-0-82930_era5_2_harvested_era5.conv._82930.gz.nc']}
                        '''
                        cleaned_df_container[k] = {}

                        all_len = []
                        for F in self.dataset_per_dt[dt][k]:
                              cleaned_df_container[k][F] = {}
                              obs_tab, era5fb_tab = self.make_obstab_era5fb_dic(dataset = k , date_time = dt, File = F)
                              all_len.append( len(obs_tab['date_time'][:] ) )
                              if len(obs_tab['date_time'][:])==0: 
                                    print('skipping empty observations')
                                    del cleaned_df_container[k]
                                    continue

                              cleaned_df_container[k][F]['obs_tab']    = obs_tab         # cleaned dataframe 
                              cleaned_df_container[k][F]['era5fb_tab'] = era5fb_tab       # cleaned dataframe                                                                                                                     
                  if max(all_len) >  0: # skipping empty container dictionary
                        best_ds, combined_obs_tab, combined_era5fb_tab, combined_head_tab, selected_file = self.combine_record(dt, container = cleaned_df_container)
                  else:
                        print('Found an empty record ')
                        continue
      
                  #if 'advanced_assimilation_feedback' not in combined_obs_tab.keys() and 'era5' in best_ds:
                  #      combined_obs_tab['advanced_assimilation_feedback'] = np.array([1]*len(combined_obs_tab['date_time']) )

                  """ Storing the combined era5fb, header and observations tables"""
                  all_combined_era5fb.append(combined_era5fb_tab)
                  all_combined_obs   .append(combined_obs_tab)
                  all_combined_head  .append(combined_head_tab)


                  """ Dictionary to fill the best_ds for duplicates """
                  dt_bestds_dic[dt] = {}
                  dt_bestds_dic[dt]['best_ds'] = best_ds
                  dt_bestds_dic[dt]['len'] = len(combined_obs_tab['date_time'])


                  """ New merged recordindex and recordtimestamps indices """
                  combined_indices.append(len(combined_obs_tab['date_time']))                 
                  combined_date_time.append(dt)

                  del cleaned_df_container 
                  
                  """ Storing the selected file for the source_configuration """
                  source_files.append(selected_file)
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
            mi.pop()
            pop = np.array(mi) # removing last unecessary index  
            di['recordindex']          = ( {'recordindex' : pop.shape } , pop )

            """ Creating the combined data """
            logging.debug('*** Concatenating the observations_table ' )      
            combined_obs = {}
            ####  Writing combined observations_table dic
            for k in all_combined_obs[0].keys():                        
                  combined_obs[k]=np.concatenate([all_combined_obs[i][k][:] for i in range(len(all_combined_obs))])
                  self.tot_records = len(combined_obs[k])
                  self.write_merged(content = 'observations_table', table= {k:combined_obs[k]})
                  logging.info('*** Written observations table %s: ', k)

                  if k not in ['date_time']:
                        del combined_obs[k]

            #self.tot_records = len(combined_obs['date_time'])
            del all_combined_obs
            print(blue + 'Memory used after deleting all_merged obs_tab dic: ', process.memory_info().rss/1000000000 , cend )
 
            dateindex = combined_obs['date_time']//86400                                                                                                                                  
            date_times, indices, counts = np.unique(dateindex, return_counts = True, return_index= True)                                                                                  
            di['dateindex'] = ( {'dateindex' : indices.shape } , indices )  # considers the day only                                                                                      
            del combined_obs
            
            combined_era5fb = {}
            ####  Writing combined era5fb_table dic                                                                                                                      
            for k in all_combined_era5fb[0].keys():
                  try:
                        combined_era5fb[k]=np.concatenate([all_combined_era5fb[i][k][:] for i in range(len(all_combined_era5fb))])
                        self.write_merged(content = 'era5fb', table= {k:combined_era5fb[k]})
                        logging.info('*** Written era5fb %s: ', k)
                  except:
                        print("FAILED feedback variable " , k)

            del all_combined_era5fb
            print(blue + 'Memory used after deleting era5fb_tab dic: ', process.memory_info().rss/1000000000 , cend)

            combined_head = {}

            ####  Writing combined header_table dic                                                                                   
            for k in all_combined_head[0].keys():
                  print('head variable is', k )
                  if ( k == 'comments' or k == 'history'):
                        continue
                  try:
                        tab=np.concatenate([all_combined_head[i][k][:] for i in range(len(all_combined_head))])
                        self.write_merged(content = 'header_table', table= {k: tab})
                        logging.info('*** Written header table %s: ', k)
                  except:
                        print('FFF FAILED variable in header table', k )

            del all_combined_head
            print(blue + 'Memory used after deleting all_merged head_tab dic: ', process.memory_info().rss/1000000000 , cend)
            
            self.write_merged(content = 'recordindex', table = di)                      
            self.write_merged(content = 'cdm_tables', table= '')


            source_conf=xr.Dataset()
            source_files = np.array(source_files).astype(dtype='|S70')
            source_conf['source_file'] = ( {'source_file' : source_files.shape } , source_files )
            self.write_merged(content = 'source_configuration', table= source_conf )

            print(0)

            '''
            dateindex = combined_obs['date_time']//86400           
            date_times, indices, counts = np.unique(dateindex, return_counts = True, return_index= True)
            di['dateindex'] = ( {'dateindex' : indices.shape } , indices )  # considers the day only
            self.write_merged(content = 'recordindex', table = di)                                                                                      
            del di , combined_obs
            '''

            return 0      


            #self.tot_records = len(combined_obs[k])
      '''
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
      '''

      def combine_record(self, dt, container = ''):
            """ This is the main function that analize each record (i.e. separate ascent) and decides which one to keep as merged. 
                Extracs the observations_table, header_table and era5fb accordingly """
            
            record_dataset_legth ={}     

            other_ds   = []

            ''' I fill the dic e.g. record_dataset_legth{100:['era5_1','ncar'], 80:['bufr','igra2'] }
                i.e. the keys are the lengths, the entries are the lists of datasets '''

            duplicates = []


            for k in container.keys(): # loop over the dataset
                  if k not in other_ds:
                        other_ds.append(k)
                  for f in container[k]: # loop over the file per dataset
                          num_rec = len(container[k][f]['obs_tab']["date_time"])
                          
                          """ Storing all the reports id with the proper prefix (for each different dataset) """
                          rep_id = b''.join(container[k][f]["obs_tab"]['report_id'][0]) 
                          rep_id = self.observation_ids_merged[k] + rep_id 
                          duplicates.append( rep_id )  
                          
                          
                          #print("This " , f , ' has ' , num_rec , ' record ')
                          if num_rec not in record_dataset_legth.keys():
                                record_dataset_legth[num_rec] = {}
                                record_dataset_legth[num_rec]['best_ds'] = []
                                record_dataset_legth[num_rec]['file'] = []

                          record_dataset_legth[num_rec]['best_ds'].append(k)
                          record_dataset_legth[num_rec]['file'].append(f)

            max_entries = max(record_dataset_legth.keys())

            
            ''' best_ds is the list of longest datasets, best_datasets the list of all the datasets available including best_ds '''
            best_datasets = record_dataset_legth[max_entries]



            #if len(best_datasets) ==0: # should never happen
            #      print('wrong??? please check')
            #      return 0,0,0,0        

            """ Choosing the priority of the datasets:
                - if era5_1 or era5_2 are present, pick them (they cant be both present for the same date_time)
                - else, if igra2 is present, pick it
                - else, one of the remaining ones """

            if 'era5_2' in best_datasets and 'era5_1' not in best_datasets:  # era5_1 and era5_2 should never be both present anyway...
                  best_ds = 'era5_2'                   
            elif 'era5_1' in best_datasets and 'era5_2' not in best_datasets:
                  best_ds = 'era5_1'
            elif 'era5_1' not in best_datasets and 'era5_2' not in best_datasets and 'igra2' in best_datasets:
                  best_ds = 'igra2'
            elif 'era5_1' not in best_datasets and 'era5_2' not in best_datasets and 'igra2' not in best_datasets:
                  best_ds =  record_dataset_legth[max_entries]['best_ds'][0]  # pick the first of the list 

            best_file = record_dataset_legth[max_entries]['file'][0]

            ''' If more file are available for the same best_ds, pick the first one from the list '''
            selected_obstab, selected_era5fb = container[best_ds][best_file]['obs_tab'] , container[best_ds][best_file]['era5fb_tab']

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
                              #print('MANAGED FFF',  var)
                        except:
                              value = [b''.join(l) for l in selected_era5fb[var][0] ][0]
                              #print('VALUE IS FFF', value)
                              selected_era5fb[var] = np.array( (1, len( selected_obstab[var]) ) ).fill(value)

            """ Extracting the header """
            selected_head = self.get_header_table(dt, ds = best_ds, File = best_file )
            for var in selected_head.keys():
                  if type (selected_head[var] ) == np.ndarray and type (selected_head[var][0] ) == np.bytes_:
                        selected_head[var] = np.array( [b''.join(l) for l in selected_head[var][:] ] )

            if  'best_ds' == 'era5_1' or best_ds == 'era5_2' :
                  selected_obstab['advanced_assimilation_feedback'] = np.array([1]*len(selected_obstab['date_time']) )
            else:
                  selected_obstab['advanced_assimilation_feedback'] = np.array([0]*len(selected_obstab['date_time']) )

            #code.interact(local = locals() )

            best_ds_byte = np.bytes_(best_ds, ndtype = '|S10') # converting to bytes ojbect
            arr = np.full( (1, len( selected_obstab['date_time']) ) , best_ds_byte )[0]
            selected_obstab['source_id'] = arr
            #code.interact(local = locals() )                                                                                                                 


            duplicate = b','.join(duplicates)
            selected_head['duplicates'] = np.array(duplicate)
            
            selected_file = np.bytes_(best_file.split('/')[-1])
            
            return  best_ds, selected_obstab, selected_era5fb, selected_head, selected_file


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

            """ Retrieving the attributes """
            if content in ['observations_table','header_table','era5fb']:
                  for var in table.keys():
                        if var == 'comments':
                              continue 

                        attrs_dic[var] = {}
                        try:
                              attrs_dic[var]['description']    = bytes( self.dic_type_attributes[content][var]['description']    , 'utf-8' )
                        except:
                              attrs_dic[var]['description']    = bytes( 'missing'    , 'utf-8' )
                              #print(' FFF FAILING WITH DESCRIPTION: ', var , ' ' ,  self.dic_type_attributes[content][var]['description']) # FFF CHECK WHY SOME ARE FAILING

                        try:
                              attrs_dic[var]['external_table'] = bytes( self.dic_type_attributes[content][var]['external_table'] , 'utf-8' )
                        except:
                              attrs_dic[var]['external_table'] = bytes( 'missing' , 'utf-8' )
                              #print(' FFF FAILING WITH EXTERNAL TABLE : ', var ) # FFF CHECK WHY SOME ARE FAILING                                                          


            if content == 'recordindex':  # writing the recordindex, recordtimestamp, dateindex
                  logging.info('Writing the merged record indices to the netCDF output ')
                  table.to_netcdf(out_name, format='netCDF4', engine='h5netcdf', mode='a')

            elif content == 'cdm_tables':
                  for k in self.data['cdm_tables'].keys():
                        table = self.data['cdm_tables'][k]
                        table.to_netcdf(out_name, format='netCDF4', engine='h5netcdf', mode='a', group = k)
                        logging.info('Writing the cdm table %s to the netCDF output ', k)
                  
            elif content == 'source_configuration':                 
                  table.to_netcdf(out_name, format='netCDF4', engine='h5netcdf', mode='a', group = content)
                  logging.info('Writing the source_configuration table to the netCDF output ')
            
            
            #elif content == 'standard_cdm': # writing the station_configuration and source_configuration
            #      logging.info('Writing the station_configuration and source_configurations tables to the netCDF output via xarray ')                             
            #      for k in self.standard_cdm:
            #            
            #            for v in self.data[k].keys():
            #                  dic = {v: self.data[k][v]}
            #                  try:
            #                        write_dict_h5(out_name, dic , k, self.encodings[k], var_selection=[], mode='a', attrs = ''  )
            #                        logging.info('*** Wrote the cdm table: %s', k)
            #                  except:
            #                        print('FFF FAILED: ', k , ' ' , v , ' ' )
            #
 
            # Writing the observations_table, header_table, era5fb 
            elif content in ['observations_table', 'era5fb', 'header_table']: 
                  #logging.info('Writing the observations_tables to the netCDF output ') 

                  shape = ''
                  for k in table.keys(): 
                        if k == 'index' or k == 'hdrlen' or 'string' in k :
                              continue
                        var_type = self.dic_type_attributes[content][k]['type']

                        ''' trying to convert the variable types to the correct types stored as attribute, read from the numpy dic file '''
                        if type(table[k][0]) != var_type:

                              if k == 'hdrlen': continue
                              try:
                                    table[k] = table[k].astype( var_type ) 
                              except:
                                    print ('FAILED converting column ' , k, ' type ', type(table[v][0]) , ' to type ', var_type )

                        print('*** Writing the table ', content, ' variable ',  k)

                        dic = {k:table[k]}  # making a 1 colum dictionary
                        shape = table[k].shape
                        #print('SHAPE IS FFF ', table[k].shape )
                        write_dict_h5(out_name, dic , content, self.encodings[content], var_selection=[], mode='a', attrs = attrs_dic  )



                  if content == 'observations_table' and self.obstab_nans_filled == False:
                        missing_cdm_var = [ v for v in self.dic_type_attributes[content].keys() if v not in self.observations_table_vars]  # variables to be filled with nans            
                        for k in missing_cdm_var:
                              if k not in ['advanced_assimilation_feedback']:
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




'''
# BISMARK
f1 = '/raid60/scratch/federico/BISMARK/0-20000-0-72764_era5_2_harvested_era5.conv._1:72764.gz.nc'
f2 = '/raid60/scratch/federico/BISMARK/0-20000-0-72764_era5_2_harvested_era5.conv._72764.gz.nc'
f3 = '/raid60/scratch/federico/BISMARK/0-20000-0-72764_era5_1_harvested_era5.conv.72764.txt.gz.nc'

bismark = '0-20000-0-72764'
bismark_dic = { 'era5_2'  : [f1, f2],
                          'era5_1'  : [f3] }


out_dir = '/raid60/scratch/federico/BISMARK/combined/'

out_dir = 'TEST_MERGEALL'

#os.system('rm /raid60/scratch/federico/BISMARK/combined/0-20000-0-72764_CEUAS_merged_v0.nc')

os.system('rm -r TEST_MERGEALL')


bufr = 'example_stations/0-20000-0-82930_bufr_harvested_era5.82930.bfr.nc'
era5_2_1 = 'example_stations/0-20000-0-82930_era5_2_harvested_era5.conv._1:82930.gz.nc'
era5_2_2 = 'example_stations/0-20000-0-82930_era5_2_harvested_era5.conv._82930.gz.nc'
era5_1 = 'example_stations/0-20000-0-82930_era5_1_harvested_era5.conv._82930.gz.nc'

if not os.path.isdir(out_dir):
      os.mkdir(out_dir)

'''

""" main block 
if __name__ == '__main__':

      Merging = Merger(out_dir)

      stat = '0-20000-0-82930'
      run_mode = 'dummy'
      stat_dic  = {'era5_2'  : [era5_2_1, era5_2_2],
                   'era5_1'  : [era5_1],
                   'bufr'    : [bufr]
 }



      stat = bismark
      stat_dic = bismark_dic

      #stat_dic  = {'era5_2_1': f1,
       #            'era5_2_2': f2}

      a = Merging.merge(stat,  datasets = stat_dic , mode = run_mode ) # single station dictionary

      logging.info('*** Convertion completed! ***')
"""

data_directories   = {     'era5_1'       :  '/raid60/scratch/federico/TEST_MAY_ERA5/era5_1'          ,
                                       'era5_2'       :  '/raid60/scratch/federico/TEST_MAY_ERA5/era5_2'          ,
                                       'era5_3188' :  '/raid60/scratch/federico/TEST_MAY_ERA5/era5_3188'    ,
                                       'era5_1759' :  '/raid60/scratch/federico/TEST_MAY_ERA5/era5_1759'   ,
                                       'era5_1761' :  '/raid60/scratch/federico/TEST_MAY_ERA5/era5_1761'   ,
                                       'ncar'           :  '/raid60/scratch/federico/TEST_MAY_ERA5/ncar'    ,
                                       'igra2'          :  '/raid60/scratch/federico/TEST_MAY_ERA5/igra2'           ,
                                       'bufr'            :  '/raid60/scratch/federico/TEST_MAY_ERA5/bufr'            ,                
                                       }

def create_stat_dic(stat_id, data_directories):
      """ Looks in the database dirdctory for files matching with the given station id """
      station_dic = {}
      total_dim = []
      for d,i in data_directories.items():
            files = os.listdir(i)
            for f in files:
                  Id = f.split('_'+d)[0]
                  if Id == stat_id:
                        if d not in station_dic.keys():
                              station_dic[d] = []                            
                        station_dic[d].append(i + '/' + f)
                        
                        total_dim. append( os.path.getsize (i + '/' + f) )
                        
                        print('FOUND!' , d , '   ' , f )
                        
      size = sum(total_dim)                 
      return station_dic, size 



out_dir = '/raid60/scratch/federico/JUNE_TEST_MERGING_PARALLEL_THREE/'
run_mode = 'dummy'

if __name__ == '__main__':

      parser = argparse.ArgumentParser(description="Make CDM compliant netCDFs")
      parser.add_argument('--stations' , '-s', 
                      help="List of station ids to merge"  ,
                      type = str)

      parser.add_argument('--maximum' , '-m', 
                      help="Maximum size of a single station file "  ,
                      type = int,
                      default = 1000000000) # 1 GB
      
      args = parser.parse_args()
      stations = args.stations
      max_size = args.maximum


      """ Initialize the merging class """
      Merging = Merger(out_dir)

      
      for station in stations.split(','):
            print(station)
            station_dic, size = create_stat_dic(station, data_directories)
            
            if size > max_size and size < 50*10**6:
                  print('For now skipping large files ')
                  continue 
 
            try:            
                  a = Merging.merge(station,  datasets = station_dic , mode = run_mode ) # single station dictionary
            except:
                  print('Failed: ', station )
            
            
            
            
            
# 0-20000-0-82930,0-20000-0-03976 
