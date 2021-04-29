import os
import sys
import netCDF4 as nc
import pandas as pd



from pathlib import Path

import numpy as np
#import argparse
from datetime import datetime, timedelta  
import numpy.ma as ma
#import math
import h5py as h5py
import xarray as xr 
#from numba import njit
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


sys.path.append('../../harvest/code')
from harvest_convert_to_netCDF import load_cdm_tables , write_dict_h5 

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


#os.system(' rm 0-20000-0-82930_CEUAS_merged_v0.nc ')
#os.system(' cp 0-20000-0-82930_CEUAS_merged_v0_KEEP_NORMALMERGED.nc 0-20000-0-82930_CEUAS_merged_v0.nc ')


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


class MergedFile(object):
    
    def __init__(self, out_dir = '' , station_id = '' , file='', copy = ''):
        self.out_dir = out_dir 
        
        if not os.path.isdir(out_dir):
            print('+++ Creating output directory: ' , out_dir )
            os.mkdir(out_dir)
            
        self.station_id = station_id
        self.file = file 
        self.summary_file =  station_id + '_' + datetime.now().strftime("%Y-%m-%d_%H-%M-%S") + '_summary.txt'
        self.std_plevs    = [1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 85000, 92500, 100000]
        
        #self.dic_type_attributes = np.load('../../merge/dic_type_attributes.npy', allow_pickle= True).item()
        #self.encodings = np.load('../../merge/groups_encodings.npy' , allow_pickle = True ).item()
        
        self.dic_type_attributes = np.load('dic_type_attributes.npy', allow_pickle= True).item()
        self.encodings = np.load('groups_encodings.npy' , allow_pickle = True ).item()        
        
        self.copy = copy 
        
    def load(self):
        
        """ Global container of the data """
        data = {}
        data['cdm_tables'] = {}
        
        
        #cdm_tables = ['crs', 'observed_variable', 'sensor_configuration', 'station_configuration_codes', 'station_type', 'units', 'z_coordinate_type']
        cdm_tables = ['crs', 'observed_variable', 'station_configuration_codes', 'station_type', 'units', 'z_coordinate_type']
        
        h5py_file = h5py.File(self.file, 'r+')

        data['h5py_file'] = h5py_file 
        data['station_id'] = self.station_id 
        
        data['recordtimestamp']               = xr.open_dataset (self.file, engine = 'h5netcdf' , decode_times = False )['recordtimestamp']
        data['recordtimestampdecoded'] = xr.open_dataset (self.file, engine = 'h5netcdf' , decode_times = True )['recordtimestamp'].values
        
        data['recordindex']         = xr.open_dataset (self.file, engine = 'h5netcdf' , decode_times = False )['recordindex'].values
        data['dateindex']            = h5py_file['dateindex']
        
        #data['sensor_id'] = h5py_file['observations_table']['sensor_id']
        data['length_max'] = len(h5py_file['observations_table']['date_time'] )

        h5py_file.close() # to check
        
        for t in cdm_tables:
            table = xr.open_dataset (self.file, engine = 'h5netcdf' , group = t)
            data['cdm_tables'][t] = table
            table.close()
            
        self.data = data 
    
    def load_obstab_era5fb(self, era5fb=False, obs_tab=True):
        """ Loading era5fb and observations_table if necessary """

        if obs_tab:
            obs_tab = xr.open_dataset (self.file, engine = 'h5netcdf'      , group = 'observations_table', decode_times = True )
            #obs_tab_vars = ['z_coordinate_type' , 'z_coordinate' , 'observation_value' , 'observed_variable', 'date_time', 'sensor_id']
            self.data['obs_tab'] = obs_tab
            attr_dic = self.retrieve_attributes(obs_tab)
            self.attr_dic = attr_dic 
            
            self.fixed_columns            = ['observed_variable' , 'observation_value', 'z_coordinate' , 'z_coordinate_type', 'secondary_value' , 'value_significance' , 'original_precision', 'date_time'] # to be defined explicitly                     
            self.unavailable_columns = ['observation_id' , 'report_id' , 'source_id'] 
            self.other_columns           = [c for c in self.data['obs_tab'].keys() if c not in self.fixed_columns and c not in self.unavailable_columns] # to be calculated by averaging 
            
        if era5fb:
            #era5fb_tab = xr.open_dataset (file, engine = 'h5netcdf' , group = 'era5fb', decode_times = True ,  drop_variables=None ) # biascorr@body ,  biascorr_fg@body 
            era5fb_tab = xr.open_dataset (self.file, engine = 'h5netcdf' , group = 'era5fb', decode_times = True )
            
            self.data['era5fb'] = era5fb_tab
        
    def write_summary(self, what = '', done = False):
        """ Write report summary """
        a = open(self.summary_file , 'a')
        a.write(self.file + '_' + what + '_' + str(done) )      
        
    def retrieve_attributes(self, obs_tab ):
        """ Retrieving the attributes for the observations table columns, stored in the dic_type_attributes.npy file """
        attrs_dic = {}
        for var in obs_tab:
            attrs_dic[var] = {}
            try:
                attrs_dic[var]['description']    = bytes( self.dic_type_attributes['observations_table'][var]['description']    , 'utf-8' )
            except:
                attrs_dic[var]['description']    = bytes( 'missing'    , 'utf-8' )    
            try:
                attrs_dic[var]['external_table'] = bytes( self.dic_type_attributes['observations_table'][var]['external_table'] , 'utf-8' )
            except:
                attrs_dic[var]['external_table'] = bytes( 'missing' , 'utf-8' )
        return attrs_dic
                
                
                
                
class Sensor(MergedFile):
    
    def __init__(self, MF = '', copy = ''  ):
        
        #MF = MergedFile.__init__(self, out_dir = out_dir, station_id = station_id, file= file )
        self.MergedFile = MF
        self.data = MF.data 
        self.copy = copy 
        
    def load_cdm_tables(self):
        """ Download the cdm tables definitions from the glamod GitHub. Taken from the harvester script. """
        cdmpath='https://raw.githubusercontent.com/glamod/common_data_model/master/tables/'             
        
        cdmtabledeflist=['id_scheme', 'crs', 'station_type', 'observed_variable', 'station_configuration', 'station_configuration_codes', 'observations_table', 
                                     'header_table', 'source_configuration', 'sensor_configuration','units' , 'z_coordinate_type']  
        cdm_tabdef = dict()
        for key in cdmtabledeflist:
            url='table_definitions'.join(cdmpath.split('tables'))+key+'.csv' # https://github.com/glamod/common_data_model/tree/master/table_definitions/ + ..._.dat 
            f=urllib.request.urlopen(url)
            col_names=pd.read_csv(f, delimiter='\t',quoting=3,nrows=0,comment='#')
            f=urllib.request.urlopen(url)
            tdict={col: str for col in col_names}
            cdm_tabdef[key]=pd.read_csv(f,delimiter='\t',quoting=3,dtype=tdict,na_filter=False,comment='#')
            
        self.cdm_tabdef = cdm_tabdef 
        
    def load_Schroeder_tables(self):
        """ Load the Schroeder's tables """
        
        sch_file = 'data/vapor.library.2' 
        
        cdm = {} # dictionary holding the CDM tables and Schr. table (that will become the sensor_configuration table)
        
        dtypes={ 'station_id': np.int32, 'latitude':str, 'longitude':str, 'altitude':str,
                        'rstype':'S4','datetime':np.int,'date_flag':'S2','Station Name':'S60'}
        
        names=list(dtypes.keys())
        
        # Metadata Schroeder -> columns: Index(['station_id', 'rstype', 'datetime', 'date_flag', 'Station Name', 'latitude', 'longitude', 'altitude'], dtype='object')
        cdm['metadata_schroeder']=pd.read_csv(sch_file,  sep=':', header=0, dtype=dtypes, names=names )        
        for l in 'latitude','longitude','altitude':
            cdm['metadata_schroeder'][l]=pd.to_numeric(cdm['metadata_schroeder'].pop(l), errors='coerce')
            
        # Sensor_configuration CDM table 
        cdms=pd.read_csv('data/vapor.instruments.all' , sep=':', names=('sensor_id','comments') )
        
        # CDM tables definitions, loaded with load_cdm_tables
        cdmd = self.cdm_tabdef
        
        cdm['sensor_configuration']=pd.DataFrame(columns=cdmd['sensor_configuration'].element_name)
        for c in cdm['sensor_configuration'].columns:
            if c not in ('sensor_id','comments'):   
                cdm['sensor_configuration'][c]=cdm['sensor_configuration'].pop(c).astype('int64')
        cdm['sensor_configuration'].sensor_id=cdms['sensor_id']
        cdm['sensor_configuration'].comments=cdms['comments']
    
        for k in range(len(cdm['sensor_configuration'].sensor_id)):
            cdm['sensor_configuration'].sensor_id.values[k]=cdm['sensor_configuration'].sensor_id.values[k].strip()
        cdm['sensor_configuration']['sensor_id']=cdm['sensor_configuration'].pop('sensor_id').astype('|S4')
            
        cdm['sensor_configuration']['comments']=cdm['sensor_configuration'].pop('comments').astype('|S200')
        
        
        """ Adding the table from WMO gruan """
        wmo = pd.read_csv('data/table_BUFR_radiosonde.csv' , sep=',' , header=1 , names = ['date', 'table_1', 'sensor_id', 'comments'] )
        wmo = wmo[['sensor_id', 'comments']]        
        
        for c in [f for f in  cdm['sensor_configuration'].columns ]:
            if c in ['sensor_id', 'comments']:
                continue
            col = np.full( (len(wmo)) , np.nan , dtype=type(cdm['sensor_configuration'][c][1] ) )
            wmo[c] = col 
            
        sensor_conf = pd.concat( [cdm['sensor_configuration'] , wmo ]) # concatenating the two tables 
        sensor_conf['comments'] = sensor_conf['comments'].str.encode('utf-8')
        cdm['sensor_configuration'] = sensor_conf
        self.cdm = cdm # saving the tables
        
    def write_sensorconfig(self):
        
        cdmd = self.cdm_tabdef
        cdm = self.cdm
        
        groups={}
        groupencodings={}
        
        for k in ['sensor_configuration']:  
            groupencodings[k]={} 
            
            groups[k]=xr.Dataset() 

            for i in range(len(cdmd[k])): 
                d=cdmd[k].iloc[i] 
                
                print(' Analyzing ', d.element_name)
                
                groups[k][d.element_name]=({k+'_len':len(cdm[k])}, cdm[k][d.element_name].values) # element_name is the netcdf variable name, which is the column name of the cdm table k 
                
                try:
                    groups[k][d.element_name].attrs['external_table']  = d .external_table # defining variable attributes that point to other tables (3rd and 4th columns)
                    groups[k][d.element_name].attrs['description']       = d.description  # it fails when trying with the observations_table 
                    print('*** Setting the attributes to the sensor_configuration table ' , d.element_name )
                    
                except KeyError:
                    print('Failed --- ')
                    pass
                
                groupencodings[k][d.element_name]={'compression': 'gzip'}

        # self.data['crs'].to_netcdf(self.file, format='netCDF4', engine='h5netcdf',group='ciao', mode='a') #
                
        for k in groups.keys():  
            try:           
                groups[k].to_netcdf(self.MergedFile.file, format='netCDF4', engine='h5netcdf', encoding=groupencodings[k], group=k, mode='a') #
                print('+++ Written  group: ' , k  )
            except  KeyError:
                print('--- Passing variable ' )
                    
        return 'done'           
                    
                    
    def extract_sensor_id(self):
        """ Extract the sensor id from the Schroeder's data,
              map the sensor id to the datetime list from the station data """
        
        try:
            
            station_id = int(self.data['station_id'])
        except:
            station_id = 99999999 # dummy station id whcih does not exist in Schroeder's table
        # select the data from the Schroeder's DF for this station_id
        
        sch_df = self.cdm['metadata_schroeder']
        located_df = sch_df.loc[ ( sch_df['station_id'] ==  station_id )]
        
        sensor_datetime = {} # I fill a dic with the datetime (closest found in the station file) and the sensor id (extracted from Schroed. data)
        """ The idea is the following:
              - I need to map the date_time contained in the merged file with the timestamps I read from Schroeder file.
              For this, I check which date_time is the closest to the time stamp, since they will hardly ever be identical. This is why I minimize the time distance to match them.
              - this works in most cases, however there are weird time stamps in the Schroeder's table, e.g. 194611000000 that has not day.
                 For now I have to skip these, since I have no clear solution. However, I find that the sensor type is unidentified in such cases, most of the times.
                 So the information is somewhat irrelevant.       
              """
        def nearest(datetimes_list, dt):
            return min(datetimes_list, key=lambda x: abs(x - dt))
        
        for i,j in located_df.iterrows():
            dt = str(j.datetime)
            sensor = j.rstype
            
            if dt == '0':
                ddt =  np.datetime64( datetime(1900, 1, 1, 0 , 0 ) ) 
                sensor_datetime[ddt] = {'sensor' : sensor, 'min_index' : 0}
                
            else:
                year, month, day, hour, minute =  int(dt[0:4]) , int(dt[4:6]) , int(dt[6:8]) , int(dt[8:10])  , int(dt[10:12] ) 
                try:
                    dt = np.datetime64( datetime ( year, month, day, hour, minute ) )
                    near = nearest(self.data['recordtimestampdecoded'], dt)
                    sensor_datetime[near] = {'sensor' : sensor } 
                    
                except: # this does not work because I cannot match the date_time anymore   
                    pass
                

        """ If I have only one entry, this instrument will be applied to the entire list of source_id  """
        lista = list(sensor_datetime.keys())
        lista.sort() 
        if len ( lista ) == 1:
            sensor_datetime[ lista[0] ] ['max_index'] = -1
            sensor_datetime[ lista[0] ] ['min_index']  = -1
            
        else: # if I have more entries, I have different instruments for different period
            
            for dt, num in zip(lista, range(len(lista)) ) :
                if num == 0 :
                    index_dt = np.where( self.data['recordtimestampdecoded'] == lista[num+1] )[0][0] # getting the following datetime 
                    index_max = self.data['recordindex'][index_dt]
                    
                    sensor_datetime[dt]['max_index'] = index_max # temporary arbitrary large number (to be reduced to the actual size of the observations_table)
                    sensor_datetime[dt]['min_index'] = 0
                    
                elif num > 0 and dt != lista[-1]: # all the values until the last one in the list (most recent datetime available in Schroeder's data)
                    index_dt = np.where( self.data['recordtimestampdecoded'] == lista[num+1] )[0][0]  # getting the following datetime 
                    index_max = self.data['recordindex'][index_dt]
    
                    sensor_datetime[dt]['max_index'] = index_max 
                    sensor_datetime[dt]['min_index'] = sensor_datetime[lista[num-1]]['max_index']
                    
                else:
                    sensor_datetime[dt]['max_index'] = self.data['length_max']                      
                    sensor_datetime[dt]['min_index'] = sensor_datetime[lista[num-1]]['max_index']

        self.sensor_datetime = sensor_datetime

    def replace_sensor_id(self):
        
        self.data['h5py_file'] = h5py.File(self.MergedFile.file, 'r+')
        
        print(' *** Replacing the sensor_id *** ')
        #replace = self.data['h5py_file']['observations_table']['sensor_id'][:]
        try:    
            del self.data['h5py_file']['observations_table']['sensor_id']
        except:
            pass
        
        sensor_datetime = self.sensor_datetime 
        lista = list(sensor_datetime.keys())
        lista.sort()
        
        """ Looping over the datetime and extracting the sensor id with the indices in the observations_table """
        temp_sensor_list = [] 
        for dt in lista:
            sensor_id = sensor_datetime[dt]['sensor']
            while len(sensor_id) < 4:
                sensor_id = sensor_id + b' ' 
                #print(sensor_id)
                
            """ If I have only one sensor, it will be applied to the whole length of the observations_table """    
            
            if  sensor_datetime[dt]['max_index'] == -1:
                length = self.data['length_max']
            else:
                length = sensor_datetime[dt]['max_index'] - sensor_datetime[dt]['min_index']
            sensor_id_array = np.full( ( length )  , sensor_id ).astype(  np.dtype('|S4')  ) 
            temp_sensor_list.append(sensor_id_array)
        
        #if not temp_sensor_list:
        #    print('No sensor found ::: ')
        #    return 'NoSensor'  # stop here if no sensor is found 
        
        """ Checking if index dimension variable exists, if not create it """
        if 'index' not in self.data['h5py_file']['observations_table'].keys():
            index = np.zeros (  self.data['h5py_file']['observations_table']['date_time'] .shape[0], dtype= int)           
            self.data['h5py_file']['observations_table'].create_dataset('index', data=index)
            
        #temp_sensor_list = []  # TO DO CHANGE !!!!!  ONLY FOR DEVELOPMENT / TESTING 
        
        
        # must replace sensor_ids starting from 2013-01-01, find the correct index in the obstable / era5fb table
        index_replace = np.searchsorted( self.data['recordtimestampdecoded'], np.datetime64('2013-01-01') )
        
        add_era5 = True
        if index_replace == len(self.data['recordtimestampdecoded'] ):
            add_era5 = False
            
        index_max = self.data['recordindex'][index_replace-1]
        
        if temp_sensor_list: # case where I found some sensor_ids inside Schroeder's table 
            
            sensor_list = np.concatenate(temp_sensor_list)
            #slen=len(sensor_list[0]) # =3
            
            slen= 4
            
            # combining the two sensor lists 
            if add_era5:
                sensor_list_combined = sensor_list[:index_max] # extracting sensors only before 2013 
                era5_sensor = self.data['h5py_file']['era5fb']['sonde_type@conv'][index_max:]
                sensor_list_combined = np.append(sensor_list_combined,era5_sensor).astype('|S4')
            else:
                sensor_list_combined = sensor_list
                
            self.data['h5py_file']['observations_table'].create_dataset('sensor_id', data = sensor_list_combined.view('S1').reshape(sensor_list_combined.shape[0], slen ), 
                                                                        compression = 'gzip' ,  chunks=True)                
            
            s = 'string{}'.format(slen)
            stringa=np.zeros(slen, dtype='S1')         
            #try: # TO DO check what this is for 
            #    del  self.data['h5py_file']['observations_table'][s]
            #    del self.data['observations_table']['index']
            #except:
            #    pass
            
            """ Adding missing dimensions """
            try:            
                self.data['h5py_file']['observations_table'].create_dataset( 'string{}'.format(slen) ,  data=stringa[:slen]  )                
                self.data['h5py_file']['observations_table'][ 'string{}'.format(slen) ].attrs['NAME']=np.bytes_('This is a netCDF dimension but not a netCDF variable.')         
                #self.data['h5py_file']['observations_table']['sensor_id'].dims[0].attach_scale( self.data['h5py_file']['observations_table']['index'] )
                #self.data['h5py_file']['observations_table']['sensor_id'].dims[1].attach_scale( self.data['h5py_file']['observations_table'][ 'string{}'.format(slen)  ] )
                print(' *** Done with the attributes of the dimension *** ')
            except ValueError:
                print('Dimension already exist, passing ')
            
        else:
            sensor_id = b'NA '
            slen = len(sensor_id)
            s = 'string{}'.format(slen)
            stringa=np.zeros(slen,dtype='S1')
            sensor_list = np.full( (len(self.data['h5py_file']['observations_table']['index']) ), sensor_id).astype(  np.dtype('|S4')  )
            
            # combining the two sensor lists 
            if add_era5:
                sensor_list_combined = sensor_list[:index_max]
                era5_sensor = self.data['h5py_file']['era5fb']['sonde_type@conv'][index_max:]
                sensor_list_combined = np.append(sensor_list_combined, era5_sensor).astype('|S3')
            else:
                sensor_list_combined = sensor_list
                
            try:
                self.data['h5py_file']['observations_table'].create_dataset('sensor_id', data = sensor_list_combined.view('S1').reshape( len(sensor_list_combined), slen ) , 
                                                                            compression = 'gzip' ,  chunks=True)       
                
                self.data['h5py_file']['observations_table']['sensor_id'].dims[0].attach_scale(  self.data['h5py_file']['observations_table']['index'] )  
                self.data['h5py_file']['observations_table'].create_dataset( s ,  data=stringa[:slen]  )                
                self.data['h5py_file']['observations_table']['string{}'.format(slen)].attrs['NAME']=np.bytes_('This is a netCDF dimension but not a netCDF variable.')                             
                self.data['h5py_file']['observations_table']['sensor_id'].dims[1].attach_scale( self.data['h5py_file']['observations_table'][ s ])
                
                print(' *** Done with the attributes of the dimension *** ')
            except ValueError:
                print('Dimension already exist, passing ')
                
        self.data['h5py_file'].close()
            
    def run(self):

        if self.copy: # if Ttue, then I create a copy of the file before adding the sensor id to avoid possible corruptions of the merged file
            os.system('cp  ' + self.MergedFile.file + '   ' +  self.MergedFile.file.replace('.nc', '_beforeSensor.nc') )
            
        cdm_tables = self.load_cdm_tables()      
        load_Schroeder_table = self.load_Schroeder_tables()
        dummy = self.extract_sensor_id()
        status = self.write_sensorconfig()
        status = self.replace_sensor_id()

        if not self.copy:
            os.system('mv  ' + self.MergedFile.file + '   ' +  self.MergedFile.out_dir )
        
        print(' --- Done writing the output file ' + self.MergedFile.file + '  ! ---  ' )            
            
            



def wrapper(out_dir = '' , station_id = '' , file = '', copy = copy ):
    """ To use to run a single file, during merging """

    MF = MergedFile(out_dir = out_dir , station_id = station_id , file = file, copy = copy  )
    data_dummy = MF.load()
    sensor = Sensor(MF = MF, copy = copy)  # loading sensor class 
    run_dummy = sensor.run() # running module     



if __name__ == '__main__':
        
            parser = argparse.ArgumentParser(description="Postprocessing Utility for the addition if instrumentation metadata")
            
            
            parser.add_argument('--force_run' , '-f', 
                                  help="Force running the file(s)"  ,
                                  type = str,
                                  default = 'False' )
 
            
            args = parser.parse_args()
            
            force_run                = args.force_run
            

            
            def run(merged_directory, postprocessed_new, force_run, station):
                
                file = merged_directory + '/' + station
              
                #station_id = file.split("_CEUAS")[0].split(merged_directory)[1].replace('/','').split('-')[-1]
                station_id = file.split("_CEUAS")[0].split(merged_directory)[1].replace('/','').split('-')[-1]
                '''
                processed = [ s.split('_')[0] for s in os.listdir(postprocessed_new) ]  # skipping processed files 

                if station_id in processed:
                    print('Already processed:::' , file )
                    return
                print (' I will process the file ::: ' , file , ' ::: station_id ::: ' , station_id )  
                '''

                """ Initialize classes """
                MF = MergedFile(out_dir = postprocessed_new , station_id = station_id , file = file, copy = True )
                """ Read data """
                data = MF.load()
            
                if force_run in ['yes', 'y', 'YES', 'Y', 'True', 'true']:   # still to implement 
                    print("    === Running in force mode ===     ")
                
                    try:                        
                        """ Running sensor module """
                        sensor = Sensor( MF = MF , copy = True )  # loading sensor class 
                        run = sensor.run() # running module 
                            
                        print(' *** Done Writing Sensor Information for file ' , file )    
                        
                    except:
                        print(" The file " + file + ' hase failed! MUST REDO! ')
                        a = open('Failed.txt', 'a')
                        a.write(file + '\n')
                    
                else:
                    print('   +++ Running in normal mode +++')
                    """ Running sensor module """
                    sensor = Sensor( MF = MF )  # loading sensor class 
                    run = sensor.run() # running module 

                    print(' *** Done Writing Sensor Information for file ' , file )    


            #station = stations_list[0]
            #a = run(merged_directory, postprocessed_new, force_run, station)
            

            """ File source directory """
            merged_directory = '../../merge/PROVA/'
            
            """ Moving postprocessed files to new directory """
            #postprocessed_new = '/raid60/scratch/federico/DATABASE_MARCH2021_sensor'

            postprocessed_new = '/raid60/scratch/federico/PROVA_sensor_newsensors'
            os.system('rm -r  /raid60/scratch/federico/PROVA_sensor_newsensors/')

            os.system('mkdir ' + postprocessed_new)
            
        
            stations_list = [ s for s in os.listdir(merged_directory) if 'empty'  not in s ]           
            #stations_list = [ s for s in stations_list if 'Sensor'   in s ]           
            
            processed = [ s.split('_')[0] for s in os.listdir(postprocessed_new) ]  # skipping processed files 
            processed = []
            
            cleaned_list = []

            for file in stations_list:
                station_id = file.split("_CEUAS")[0]
                
                if station_id in processed:
                    print('Already processed:::' , file )
                else:
                    cleaned_list.append(file)
            


            print(cleaned_list)
            for s in cleaned_list:
                a = run(merged_directory, postprocessed_new, force_run, s)
                
            
            #p = Pool(30)
            #func = partial(run, merged_directory, postprocessed_new, force_run)
            #out = p.map(func, cleaned_list)        
            
