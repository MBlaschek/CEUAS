""" Utility to postprocess the merged files from v1 of the CEAUS database
    - Add instrument type from Schroeder's list
    - ...
"""


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
import urllib.request

sys.path.append('../harvest/code')
from harvest_convert_to_netCDF_newfixes import load_cdm_tables 


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


os.system(' rm 0-20000-0-82930_CEUAS_merged_v0.nc ')
os.system(' cp 0-20000-0-82930_CEUAS_merged_v0_KEEP_NORMALMERGED.nc 0-20000-0-82930_CEUAS_merged_v0.nc ')


class MergedFile(object):
    
    def __init__(self, out_dir = '' , station_id = '' , file='' ):
        self.out_dir = out_dir 
        self.station_id = station_id
        self.file = file 
    
    def load(self, file=''):
        
        data = {}
        
        h5py_file = h5py.File(file, 'r+')
        
        data['h5py_file'] = h5py_file 
        data['station_id'] = self.station_id 
        
        data['recordtimestamp']               = xr.open_dataset (file, engine = 'h5netcdf' , decode_times = False )['recordtimestamp']
        data['recordtimestampdecoded'] = xr.open_dataset (file, engine = 'h5netcdf' , decode_times = True )['recordtimestamp']
        data['crs']                                     = xr.open_dataset (file, engine = 'h5netcdf' , decode_times = True, group = 'crs')
        
        data['recordindex']         = h5py_file['recordindex']
        data['dateindex']            = h5py_file['dateindex']
        
        #data['sensor_id'] = h5py_file['observations_table']['sensor_id']
        data['length_max'] = len(h5py_file['observations_table']['date_time'] )
        
        print(0)
        
        return data 


        

class Sensor(MergedFile):
    """ Main class to extract the sensor type from Schroeder's list, 
        and write it back to the observations_table
        Moreover, it creates the sensor_configuration table """

    def __init__(self, data = '', file='' , out_dir = '' , station_id = '' ):
        
        self.data = data 
        
        MergedFile.__init__(self, out_dir = out_dir, station_id = station_id, file= file )


    def load_cdm_tables(self):
        
        cdmpath='https://raw.githubusercontent.com/glamod/common_data_model/master/tables/' # cdm tables            
        
        """ Selecting the list of table definitions. Some of the entires do not have the corresponding implemented tables """
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
        
        cdm = {} # dictionary holding the CDM tables and the Schr table 
        
        dtypes={ 'station_id': np.int32, 'latitude':str, 'longitude':str, 'altitude':str,
                        'rstype':'S4','datetime':np.int,'date_flag':'S2','Station Name':'S60'}
        
        names=list(dtypes.keys())
        
        # Metadata Schroeder -> columns: Index(['station_id', 'rstype', 'datetime', 'date_flag', 'Station Name', 'latitude', 'longitude', 'altitude'], dtype='object')
        cdm['metadata_schroeder']=pd.read_csv(sch_file,  sep=':', header=0, dtype=dtypes, names=names )        
        for l in 'latitude','longitude','altitude':
            cdm['metadata_schroeder'][l]=pd.to_numeric(cdm['metadata_schroeder'].pop(l), errors='coerce')
            
        # Sensor_configuration CDM table 
        cdms=pd.read_csv('data/vapor.instruments.all' , sep=':', names=('sensor_id','comments') )
        
        # CDM tables definitions 
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
        
        self.cdm = cdm 
        
        
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
                    groups[k][d.element_name].attrs['external_table']=d.external_table # defining variable attributes that point to other tables (3rd and 4th columns)
                    groups[k][d.element_name].attrs['description']=d.description  # it fails when trying with the observations_table 
                    print('*** Setting the attributes to the sensor_configuration table ' , d.element_name )
                    
                except KeyError:
                    print('Failed --- ')
                    pass
                
                  
                groupencodings[k][d.element_name]={'compression': 'gzip'}

        # self.data['crs'].to_netcdf(self.file, format='netCDF4', engine='h5netcdf',group='ciao', mode='a') #
                
        for k in groups.keys():  
            try:           
                groups[k].to_netcdf(self.file, format='netCDF4', engine='h5netcdf', encoding=groupencodings[k], group=k, mode='a') #
                print('+++ Written sensor_configuration ' )
            except  KeyError:
                print('--- Passing variable ' )
                    
                    
                    
                    
    def extract_sensor_id(self):
        """ Extract the sensor id from the Schroeder's data,
              map the sensor id to the datetime list from the station data """
        
        station_id = int(self.data['station_id'])

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
                sensor_datetime[ddt] = {'sensor': sensor, 'min_index': 0}
                
            else:
                year, month, day, hour, minute =  int(dt[0:4]) , int(dt[4:6]) , int(dt[6:8]) , int(dt[8:10])  , int(dt[10:12] ) 
                try:
                    dt = np.datetime64( datetime ( year, month, day, hour, minute ) )
                    near = nearest(self.data['recordtimestampdecoded'].values, dt)
                    sensor_datetime[near] = {'sensor': sensor } 
                    
                except: # this does not work because I cannot match the date_time anymore   
                    pass
                    """ 
                    if month <= 0:
                        month = 1
                    if day <=0:
                        day = 1
    
                    dt = np.datetime64( datetime ( year, month, day, hour, minute ) )
                    """

            
            print (0)
            

        """ If I have only one entry, this instrument will be applied to the entire list of observations """
        lista = list(sensor_datetime.keys())
        lista.sort() 
        if len ( lista ) == 1:
            return sensor_datetime
        
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
            length = sensor_datetime[dt]['max_index'] - sensor_datetime[dt]['min_index']
            sensor_id_array = np.full( ( length )  , sensor_id ).astype(  np.dtype('|S3')  ) 
            temp_sensor_list.append(sensor_id_array)
            print(2)
        
        sensor_list = np.concatenate(temp_sensor_list)
        
        data = self.data['h5py_file']
        slen=len(sensor_list[0]) # =3
        
        # making an array of 1-byte elements
        data['observations_table'].create_dataset('sensor_id', data = sensor_list.view('S1').reshape(sensor_list.shape[0], slen ), compression = 'gzip' ,  chunks=True)
        
        s = 'string{}'.format(slen)
        stringa=np.zeros(slen,dtype='S1')
        
        try:    
            del  self.data['h5py_file']['observations_table'][s]
            del self.data['observations_table']['index']
            
        except:
            pass
        
        """ Checking if index dimension variable exists, if not create it """
        if 'index' not in self.data['h5py_file']['observations_table'].keys():
            index = np.zeros (  self.data['h5py_file']['observations_table']['date_time'] .shape[0], dtype='S1')           
            data['observations_table'].create_dataset('index', data=index)

        """ Adding missing dimensions """
        try:            
            data['observations_table']['sensor_id'].dims[0].attach_scale( data['observations_table']['index'] )
            data['observations_table'].create_dataset( s ,  data=stringa[:slen]  )
            data['observations_table']['sensor_id'].dims[1].attach_scale(data['observations_table']['string{}'.format(slen)])

            data['observations_table']['string{}'.format(slen)].attrs['NAME']=np.bytes_('This is a netCDF dimension but not a netCDF variable.')
            print(' *** Done with the attributes of the dimension *** ')
            
        except ValueError:
            print('Dimension already exist, passing ')
            
        print(0)

        self.data['h5py_file'].close()
        

if __name__ == '__main__':
        
            parser = argparse.ArgumentParser(description="Postprocessing Utility")
            
            parser.add_argument('--stations' , '-s', 
                                  help="Station to postprocess"  ,
                                  type = str,
                                  default = 'a')
        
            parser.add_argument('--instrument' , '-i', 
                                  help="Add instrument type"  ,
                                  type = bool,
                                  default = True )
            
            args = parser.parse_args()
            
            stations       = args.stations
            get_instrument = args.instrument

            # 0-20000-0-70316_CEUAS_merged_v0.nc
            
            """ for now, the file and stations are here hardcoded """

            file = '0-20000-0-70316_CEUAS_merged_v0.nc' 
            station_id = 70316

            """ Initialize classes """
            MF = MergedFile(out_dir = '' , station_id = station_id , file = file  )
            
            """ Read data """
            data = MF.load(file = file )
            
            """ Load Schroeder """
            Sensor = Sensor(data= data, out_dir = '' , station_id = station_id , file = file )
            cdm_tables = Sensor.load_cdm_tables()      
            
            load_Schroeder_table = Sensor.load_Schroeder_tables()
            
            dummy = Sensor.extract_sensor_id()
            dummy = Sensor.replace_sensor_id()
            dummy = Sensor.write_sensorconfig()
            
            print('*** Done ***')




""" To run simply type 
                         python postprocess.py 
    (you might want to modify the file path and the station id variables above) """
