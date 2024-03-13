import os
import sys
import netCDF4 as nc
import pandas as pd
pd.options.mode.chained_assignment = None


import numpy as np
#import argparse
from datetime import datetime, timedelta  
import numpy.ma as ma
#import math
import h5py as h5py

import urllib.request

import copy
from numba import njit
import code
from functools import reduce
from tqdm import tqdm

from multiprocessing import Pool
from functools  import partial

import xarray as xr

from collections import Counter

import cProfile
#cProfile.run('__main__')

import time as T 
t=T.time()

sys.path.append('../../harvest/code_cop2')

# nan int = -2147483648 


import warnings
warnings.simplefilter(action='ignore', category=FutureWarning) # deactivates Pandas warnings 

import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')  # up to INFO level, DEBUG statements will not be printed 

import argparse

from harvest_convert_to_netCDF_yearSplit import  clean_station_configuration , write_dict_h5_old

"""
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)
"""

#process = psutil.Process(os.getpid())

# https://stackoverflow.com/questions/287871/how-to-print-colored-text-in-python
cend   = '\033[0m'
blue   = '\033[34m'
red     = '\33[91m'
green = '\33[92m'

data = {}


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
     
        
    def load(self):
        
        """ Global container of the data """
        data = {}
        data['cdm_tables'] = {}
        
        h5py_file = h5py.File(self.file, 'r+')

        data['h5py_file'] = h5py_file 
        data['station_id'] = self.station_id 

        
        recordtimestamp = h5py_file['recordtimestamp'][:]
        indices = h5py_file['recordindex'][:]
                            
        data['recordtimestamp'] = recordtimestamp
        data['recordindex']  = indices
        data['recordtimestampdecoded'] = pd.to_datetime( recordtimestamp ,  unit='s',  origin=pd.Timestamp('1900-01-01') ) 

        # storing the length of each record
        record_length = [indices[i] - indices[i-1] for i in range(1, len(indices) ) ] 
        lmax = len(h5py_file['observations_table']['date_time'][:] )
        data['length_max'] = lmax
        last = lmax - indices[-1]
        record_length.append(last)
        # find length of records per timestamp
        data['record_length'] = record_length        
            
        self.data = data 



    
class Sensor(MergedFile):
    
    def __init__(self, MF = '', copy = ''  ):
        #MF = MergedFile.__init__(self, out_dir = out_dir, station_id = station_id, file= file )
        self.MergedFile = MF
        self.data = MF.data 
        self.copy = copy 
        
        
    def load_Schroeder_tables(self):
        """ Load the Schroeder's tables and the WMO codes table. 
        Create a combined table to be merged into the CDM sensor_configuration table """
        
        ### Schroeder 
        sch_file = 'data/vapor.library.2'
        cdm = {} # dictionary holding the CDM tables and Schr. table (that will become the sensor_configuration table)
        dtypes={ 'station_id': np.int32, 'latitude':str, 'longitude':str, 'altitude':str,
                        'rstype':'S4','datetime':np.int64,'date_flag':'S2','Station Name':'S60'}
        
        names=list(dtypes.keys())
        
        # Metadata Schroeder -> columns: Index(['station_id', 'rstype', 'datetime', 'date_flag', 'Station Name', 'latitude', 'longitude', 'altitude'], dtype='object')
        schr = pd.read_csv(sch_file,  sep=':', header=0, dtype=dtypes, names=names )   
        stations = np.array( [ '0' + str(s) if len(str(s)) ==4 else str(s) for s in schr['station_id'] ] )
        schr['station_id'] = stations
        
        cdm['metadata_schroeder']= schr          
        
        for l in 'latitude','longitude','altitude':
            cdm['metadata_schroeder'][l]=pd.to_numeric(cdm['metadata_schroeder'].pop(l), errors='coerce')
            
        cdms=pd.read_csv('data/vapor.instruments.all' , sep=':', names=('sensor_id','comments') )
        
        ### Sensor_configuration CDM table         
        cdmpath='https://raw.githubusercontent.com/glamod/common_data_model/master/tables/'             

        cdmtabledeflist=['sensor_configuration']  
        cdmd = dict()
        for key in cdmtabledeflist:
            url='table_definitions'.join(cdmpath.split('tables'))+key+'.csv' # https://github.com/glamod/common_data_model/tree/master/table_definitions/ + ..._.dat 
            f=urllib.request.urlopen(url)
            col_names=pd.read_csv(f, delimiter='\t',quoting=3,nrows=0,comment='#')
            f=urllib.request.urlopen(url)
            tdict={col: str for col in col_names}
            cdmd[key]=pd.read_csv(f,delimiter='\t',quoting=3,dtype=tdict,na_filter=False,comment='#')

        self.cdm_tabdef = cdmd 
        
        
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
        
        ### Adding the table from WMO gruan 
        wmo = pd.read_csv('data/table_BUFR_radiosonde.csv' , sep=',' , header=1 , names = ['date', 'table_1', 'sensor_id', 'comments'] )

            
        def make_wmo_dates(wmo):
            """ Build dates for the wmo df.
            Corresponds to the starting date of the validity of the sensor id from era5 fb. """  ### TO DO TODO boh 
            
            dates_n = []
            for d in wmo.date.values:
                if type(d) == str and '/' in d :
                    s = d.split('/')
                    d,m,y = s[0] , s[1], s[2]
                    if len(y) <3:
                        y = '20' + y 
                    dt = pd.Timestamp( '-'.join([y,m,d]) ).date()
                    dates_n.append(dt)
                        
                else:
                    dates_n.append(np.nan)
            
            wmo["date"] = dates_n
            return wmo 
                
        wmo = make_wmo_dates(wmo)
        w = wmo.dropna(subset = ['date'])
        self.wmo_sensor_update = w 
        
        #df = pd.DataFrame( {'date':dates_n , 'old_id': old_sensor , 'new_id': new_sensor}  )        
        wmo['date_start'] = wmo['date']
        wmo = wmo[['sensor_id', 'comments', 'date_start']]        
        
        for c in [f for f in  cdm['sensor_configuration'].columns ]:
            if c in ['sensor_id', 'comments' , 'date_end']:
                continue
            col = np.full( (len(wmo)) , np.nan , dtype=type(cdm['sensor_configuration'][c][1] ) )
            wmo[c] = col 
            
        sensor_conf = pd.concat( [cdm['sensor_configuration'] , wmo ]) # concatenating the two tables 
        #sensor_conf['comments'] = sensor_conf['comments'].str.encode('utf-8') TODO, this was needed in old srvx2 / srvx8
        sensor_conf['comments'] = sensor_conf['comments']      
        sensor_conf.to_csv('sensor_configuration_all.csv', sep = '\t')
        self.cdm = cdm # saving the tables
        self.wmo = wmo # (converted to CDM convention )
        
        
  
    def extract_Schroeder_id(self):
        """ Extract the useful starting and ending dates from Schroeder table """
        stat = self.data['station_id'].split('_')[0]
        schr = self.cdm['metadata_schroeder']
        schr_red = schr.loc[schr['station_id'] == stat ]

        datetimes = schr_red.datetime
        converted_timestamps = [] 
        for dt in datetimes:
            try:
                dt = str(dt)
                try:
                    year, month, day, hour, minute =  int(dt[0:4]) , int(dt[4:6]) , int(dt[6:8]) , int(dt[8:10])  , int(dt[10:12] ) 
                except:
                    year, month, day, hour, minute =  int(dt[0:4]) , int(dt[4:6]) , int(dt[6:8]) , 0, 0 

                timestamp = pd.Timestamp( year=year, month=month, day=day, hour=hour) 
            except:
                timestamp = np.nan 
            converted_timestamps.append(timestamp )

        schr_red['datetimes'] = converted_timestamps
        schr_red = schr_red.dropna(subset=['datetimes'])
        schr_red = schr_red.reset_index()

        ### check validity (min_date, max_date) of each sensor
        ### Each sensor will be valid until min(next_sensor, 2 year)
        ### i.e. if the gap is larger than 5 years, the sensor will lose its applicability 

        start, end = [], [] 
        allowed_delta = 5 * 365 # in days 
        for i, dt in  enumerate( schr_red['datetimes']):
            start.append(dt)
            if i != ( len(schr_red['datetimes'])-1 ):
                next_dt =  schr_red['datetimes'][i+1]
                if (next_dt - dt) < pd.Timedelta(days=allowed_delta ):
                    end.append(next_dt  )
                else:
                    end.append( dt + pd.Timedelta ( days=allowed_delta )  )
            else:
                end.append( dt + pd.Timedelta (days=allowed_delta )  )

        schr_red['start_date'] = start
        schr_red['end_date'] = end 
        #for dt

        self.extracted_schoreder = schr_red

        a =0

    
    def write_sensorconfiguration_table(self):
        """ Write the sensor configuration table (standard CDM)  """
        cdmd = self.cdm_tabdef
        cdm = self.cdm
        
        groups={}
        groupencodings={}
        
        for k in ['sensor_configuration']:  
            groupencodings[k]={} 
            
            groups[k]=xr.Dataset() 

            for i in range(len(cdmd[k])): 
                d=cdmd[k].iloc[i] 
                #print(' Analyzing ', d.element_name)
                groups[k][d.element_name]=({k+'_len':len(cdm[k])}, cdm[k][d.element_name].values) # element_name is the netcdf variable name, which is the column name of the cdm table k 
                
                try:
                    groups[k][d.element_name].attrs['external_table']  = d .external_table # defining variable attributes that point to other tables (3rd and 4th columns)
                    groups[k][d.element_name].attrs['description']       = d.description  # it fails when trying with the observations_table 
                    print('*** Setting the attributes to the sensor_configuration table ' , d.element_name )
                    
                except KeyError:
                    print('Failed --- ' , k , ' ', i )
                    pass
                groupencodings[k][d.element_name]={'compression': 'gzip'}
        # self.data['crs'].to_netcdf(self.file, format='netCDF4', engine='h5netcdf',group='ciao', mode='a') #
                
        for k in groups.keys():
            try:           
#                groups[k].to_netcdf(self.MergedFile.file, format='netCDF4', engine='h5netcdf', encoding=groupencodings[k], group=k, mode='a') #
                dic = {}
                for v in groups[k].keys():
                    dic[v] = groups[k][v].values
                write_dict_h5_old(self.MergedFile.file, dic , k, self.encodings[k], var_selection=[], mode='a',
                                  attrs = {'description':'Schroeder and WMO Sensor codes'}  )  
                
                #print('+++ Written  group: ' , k  )
            except  KeyError:
                print('--- Passing group ' , k )
        return 'done'           
                    
                    
                    
    def update_wmo(self, era5_sensors):
        
        updated_wmo = self.wmo_sensor_update
        file = self.MergedFile.file
        
        """ Must update the sensor id due to incosistent notation in the WMO guideline.
            See table page A-398
            https://library.wmo.int/doc_num.php?explnum_id=10235"""
        
        limit_date_era5 =  np.datetime64('2013-01-01')
        
        updated_wmo = updated_wmo[['date', 'table_1', 'sensor_id']]
        # if date in updated_wmo is < 2013: automatically convert old entries to new one
        updated_wmo_after2013    = updated_wmo.loc[updated_wmo['date'] >= limit_date_era5 ]
        updated_wmo_after2013['table_1'] = updated_wmo_after2013['table_1'].astype(int)
        
        updated_wmo_before2013 = updated_wmo.loc[updated_wmo['date'] < limit_date_era5 ]
        updated_wmo_before2013['table_1'] = updated_wmo_before2013['table_1'].astype(int)
        
        ids = list(np.unique( era5_sensors ) )
        to_replace_after = [ i for i in  ids if i in updated_wmo_after2013.table_1.astype(int).values ] 
        to_replace_before = [ i for i in  ids if i in updated_wmo_before2013.table_1.astype(int).values ] 
        
        
        conv_dic = dict(zip(updated_wmo.table_1.values , updated_wmo.sensor_id.values ))
        
        if len(to_replace_after) + len(to_replace_before)==0:
            return  era5_sensors
        
        if len(to_replace_before):
                era5_sensors = [ conv_dic[i] if i in conv_dic.keys() else  i for i in era5_sensors   ]
        
        if len(to_replace_after):
            dt = pd.to_datetime( h5py.File(self.MergedFile.file, 'r')['observations_table']['date_time'], unit='s',  origin=pd.Timestamp('1900-01-01') )
            for v in to_replace_after:
                # here I need to check the exact datetime of the sensor before replacing,
                # and compare with the exact time of the change in the id numbers.
                date = pd.to_datetime (updated_wmo_after2013.loc[updated_wmo_after2013.table_1 == v ].date.values[0] )
                try:
                    index = np.where(dt >  date)[0][0]
                except:
                    continue 
                # only replacing the value for the chunk after the index 
                to_update = era5_sensors[index:]                        
                to_update = [ conv_dic[i] if i in conv_dic.keys() else  i for i in to_update   ] 
                
                era5_sensors = list(era5_sensors[:index]) + to_update 
                        
        return era5_sensors
    
    
    def extract_sensor_id(self):
        """ Extract the sensor id from the Schroeder's data,
              map the sensor id to the datetime list from the station data """

        ### Procedure 
        # 1. load era5fb codes, use them as standard
        # 2. search for the WIGOS id in the Schroeder's table (only  of the form [0-20000-0- or 0-20001-0-] )
        # 3. replace era5fb when Possible with Schroeder 
        
        try:
            if '_' in self.data['station_id']:
                station_id = self.data['station_id'].split('_')[0] 
            else:
                station_id = self.data['station_id']
                
        except:
            station_id = 99999999 # dummy station id that does not exist in Schroeder's table

        ### Extracting all WMO possible sensor_ids
        
        wmo_sensor_id = self.data['h5py_file']['era5fb']['sonde_type@conv'][:]
        wmo_sensor_id = self.update_wmo(wmo_sensor_id)
        
        ### update the wmo codes accoridng to the changing date 
        self.data['wmo_sensor_id'] =  wmo_sensor_id 
            

        sch_df = self.extracted_schoreder
        located_df = pd.DataFrame()
        
        ### Checking lat/lon compatibility  of the Schroeder data with the station 
        if not sch_df.empty:    
            stat_lat, stat_lon = self.data['h5py_file']['observations_table']['latitude'][:][0] , self.data['h5py_file']['observations_table']['longitude'][:][0]
        
            if '0-20000-0-' in self.MergedFile.file or '0-20001-0-' in self.MergedFile.file: #only WMO codes from OSCAR must be compared with Schroeder's table
                located_df = sch_df.loc[ ( sch_df['station_id'] ==  station_id )]
                located_df = located_df.reset_index()
                
                if abs(located_df.latitude.values[0] - stat_lat) < 0.5 and  abs(located_df.longitude.values[0] - stat_lon) < 0.5 :
                    pass
                else:
                    located_df = sch_df.loc[ ( sch_df['station_id'] ==  'DUMMY' )] # Schroeder's data does not satisfy lat/lon compatibility 
            else:
                located_df = sch_df.loc[ ( sch_df['station_id'] ==  'DUMMY' )] # dummy emtpy         
        
        
        sensor_ids_all = []        
        
        if not located_df.empty: ### Extract Schroeder's data for this station
            
            datetimes = self.data['recordtimestampdecoded']
            record_lenghts = self.data['record_length']
            indices = self.data['recordindex']

            for dt, length, ind in zip(datetimes, record_lenghts, indices ) :
                find_sch = located_df.loc[ (located_df.start_date >= dt  ) &  ( located_df.end_date < dt ) ]
                
                if not find_sch.empty: # found Schroeder compatible data
                    try:
                        
                        sensor = find_sch.sensor_id.values[0]
                    except:
                        sensor = 'NA  '
                        print('Schroeder Sensor not found')
                    lista = [sensor] * length
                else:
                    lista = wmo_sensor_id[ind: (ind+length)]
                    
                sensor_ids_all.extend( lista )
                
        else:
            sensor_ids_all.extend( wmo_sensor_id )
            
        
        self.sensor_ids_all = sensor_ids_all 
        
        
            
        '''    
        def nearest(datetimes_list, dt):
            inf = np.searchsorted(datetimes_list, dt)
            sup = inf+1
            
            if abs(dt - datetimes_list[inf]) < abs(dt - datetimes_list[sup]):
                return datetimes_list[inf]
            else:
                return datetimes_list[sup] 
                
        """ If I have only one entry, this instrument will be applied to the entire list of source_id  """  ### TO DO TODO NO!!!!! 
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
        '''
        
        
    def apply_sensor_id(self):
        self.data['h5py_file'] = h5py.File(self.MergedFile.file, 'r+')
        print(' *** Replacing the sensor_id *** ')
        #replace = self.data['h5py_file']['observations_table']['sensor_id'][:]

        try:    
            del self.data['h5py_file']['observations_table']['sensor_id']
        except:
            pass  
        
        
        ### Formatting the sensor id 
        sensors = np.array(self.sensor_ids_all).astype(str)
        sensors = np.array( [ s.replace('-2147483600.0','NA ').replace('.0','') for s in sensors ] )
        sensors = sensors.astype('|S4')
        
        slen= 4 
        self.data['h5py_file']['observations_table'].create_dataset('sensor_id', data = sensors.view('S1').reshape(sensors.shape[0], slen ), 
                                                                        compression = 'gzip' ,  chunks=True)          
        
        
        # writing dimensions for the strings
        
        stringa=np.zeros(slen,dtype='S1')
        if 'string{}'.format(slen) not in self.data['h5py_file']['observations_table'].keys():          
            self.data['h5py_file']['observations_table'].create_dataset( 'string{}'.format(slen) ,  data=stringa[:slen]  )                
            self.data['h5py_file']['observations_table'][ 'string{}'.format(slen) ].attrs['NAME']=np.bytes_('This is a netCDF dimension but not a netCDF variable.')         
        #self.data['h5py_file']['observations_table']['sensor_id'].dims[0].attach_scale( self.data['h5py_file']['observations_table']['index'] )
        #self.data['h5py_file']['observations_table']['sensor_id'].dims[1].attach_scale( self.data['h5py_file']['observations_table'][ 'string{}'.format(slen)  ] )
        
        

        
        
        '''
        sensor_datetime = self.sensor_datetime 
        lista = list(sensor_datetime.keys())
        lista.sort()
        
        """ Looping over the datetime and extracting the sensor id with the indices in the observations_table """
        temp_sensor_list = [] 
        for dt in lista:
            sensor_id = sensor_datetime[dt]['sensor']
            while len(sensor_id) < 4:
                sensor_id = sensor_id + b' ' 
                
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
        # index_add_era5 == index from which one must use ERA5 
        if '20999' in self.MergedFile.file:
            index_add_era5 = 0
            
            #elif '0-20000-0' in self.MergedFile.file or '0-20001-0' in self.MergedFile.file:
        else:
            index_timestamp = np.searchsorted( self.data['recordtimestampdecoded'], np.datetime64('2013-01-01') )
            if index_timestamp ==  len(self.data['recordindex']): # the index_timestamp where era5 should be replaced must be less than the total length, or the date 01-01-2013 is not included in the data 
                index_add_era5 = -1
            elif index_timestamp < len(self.data['recordindex']):
                index_add_era5 = self.data['recordindex'][index_timestamp]
                
            
        # if index_replace == 0, the file has only ERA5 data after 2013
        # if index_replace == len(self.data['recordtimestampdecoded'] ), the file has no ERA5 data after 2013
        # if in between, part of sensor from Schroeder, part from ERA5 
        
        slen= 4 # define strings length for the sensor id
        stringa=np.zeros(slen,dtype='S1')



        def make_empty_vec(l, slen= 4):
            """ Creates an empty vector i.e. filled with b'NA ' (standard for non available sensor_id) of given length """
            sensor_id = b'NA '
            slen = len(sensor_id)
            s = 'string{}'.format(slen)
            lista = np.full( (l) , sensor_id).astype(  np.dtype('|S4')  )
            return lista 
        
      
        
        ### only use ERA5:
        # dump the whole content as sensor_id from ERA5, before correcitng the ids from the WMO table witht he correct date update
        if index_add_era5 == 0:
            sensor_list_combined = update_wmo(self.wmo_sensor_update, 
                                              self.data['h5py_file']['era5fb']['sonde_type@conv'][:].astype(int).astype('|S4'),
                                              self.MergedFile.file)

        else:
            #if index_add_era5 == len(self.data['recordtimestampdecoded'] ):  # case I only use Schroeder data if available
            if index_add_era5 == -1:  # case I only use Schroeder data if available
                if temp_sensor_list: # case where I found some sensor_ids inside Schroeder's table 
                    sensor_list_sch = np.concatenate(temp_sensor_list)
                    sensor_list_combined = sensor_list_sch
                    #sensor_list_combined = sensor_list_sch[:(index_add_era5-1)] 
                else: # empty Schroeder and no ERA5 -> all NANS
                    sensor_list_combined = make_empty_vec(len(self.data['h5py_file']['observations_table']['index']),  slen= 4)

            else: # here I have to check the dates and see if Schr is available. There is for sure data both before and after 2013 
                if temp_sensor_list:
                    sensor_list_sch = np.concatenate(temp_sensor_list)
                    sensor_list_sch = sensor_list_sch[ :index_add_era5 ]      
                    #sensor_list_sch = sensor_list_sch[:(index_add_era5-1)]      
                    
                    sensor_list_era5 = self.data['h5py_file']['era5fb']['sonde_type@conv'][index_add_era5:].astype(int)
                    sensor_list_era5 = update_wmo(self.wmo_sensor_update, sensor_list_era5, self.MergedFile.file)
                    sensor_list_combined = np.append(sensor_list_sch, sensor_list_era5).astype('|S4')                    
                else:
                    # try this fix 
                    #sensor_list_combined = make_empty_vec( index_add_era5-1,  slen= 4) # until 2013 here empty Schroeder
                    sensor_list_combined = make_empty_vec( index_add_era5,  slen= 4) # until 2013 here empty Schroeder                    
                    #sensor_list_era5 = self.data['h5py_file']['era5fb']['sonde_type@conv'][index_add_era5-1:].astype(int)
                    sensor_list_era5 = self.data['h5py_file']['era5fb']['sonde_type@conv'][index_add_era5:].astype(int)
                    
                    sensor_list_era5 = update_wmo(self.wmo_sensor_update, sensor_list_era5, self.MergedFile.file)                    
                    sensor_list_combined = np.append(sensor_list_sch,sensor_list_era5).astype('|S4')                    

        slen= 4 # define strings length for the sensor id
        stringa=np.zeros(slen,dtype='S1') 
        
        try:
            self.data['h5py_file']['observations_table'].create_dataset('sensor_id', data = sensor_list_combined.view('S1').reshape(sensor_list_combined.shape[0], slen ), 
                                                                            compression = 'gzip' ,  chunks=True)                
        except:
            sensor_list_combined = make_empty_vec( len(self.data['h5py_file']['observations_table']['index']) ,  slen= 4)
            self.data['h5py_file']['observations_table'].create_dataset('sensor_id', data =sensor_list_combined, compression = 'gzip' ,  chunks=True )     
             
             
             
        """ Adding missing dimensions ??? """
        try:            
            self.data['h5py_file']['observations_table'].create_dataset( 'string{}'.format(slen) ,  data=stringa[:slen]  )                
            self.data['h5py_file']['observations_table'][ 'string{}'.format(slen) ].attrs['NAME']=np.bytes_('This is a netCDF dimension but not a netCDF variable.')         
            #self.data['h5py_file']['observations_table']['sensor_id'].dims[0].attach_scale( self.data['h5py_file']['observations_table']['index'] )
            #self.data['h5py_file']['observations_table']['sensor_id'].dims[1].attach_scale( self.data['h5py_file']['observations_table'][ 'string{}'.format(slen)  ] )
            print(' *** Done with the attributes of the dimension *** ')
        except:
            print('WRONG adding missing dimension, might not be necessary'  )

        self.data['h5py_file'].close()
            
        '''
        
        
    def run(self):
        
        if self.copy: # if Ttue, then I create a copy of the file before adding the sensor id to avoid possible corruptions of the merged file
            os.system('cp  ' + self.MergedFile.file + '   ' +  self.MergedFile.file.replace('.nc', '_beforeSensor.nc') )
            
        load_Schroeder_table = self.load_Schroeder_tables()
        dummy = self.extract_Schroeder_id()
        
        status = self.write_sensorconfiguration_table() # write sensorconfig table to netcdf
        
        status = self.extract_Schroeder_id() # extract possible schroeder data 
        
        dummy = self.extract_sensor_id() # find best available sensor  
        
        #status = self.make_sensor_id()
        
        status = self.apply_sensor_id()
        
        a = 0
        
        if not self.copy:
            os.system('mv  ' + self.MergedFile.file + '   ' +  self.MergedFile.out_dir.replace('_beforeSensor','') )
        
        print(' --- Done writing the output file ' + self.MergedFile.file + '  ! ---  ' )            
            
            
            
def wrapper(out_dir = '' , station_id = '' , file = '', copy = copy ):
    """ To use to run a single file, to be used during merging """

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
                station_id = file.split("_CEUAS")[0].split(merged_directory)[1].replace('/','').split('-')[-1]

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
            merged_directory = '/mnt/users/scratch/leo/scratch/raymerge/MERGED_YEARLY_27FEB2024'
            
            """ Moving postprocessed files to new directory """
            postprocessed_new = '/mnt/users/scratch/leo/scratch/sense/' 
                        
            #source = ['/scratch/das/federico/MERGED_YEARLY_22NOV_FULL_checkDimensions/0-20001-0-11035/0-20001-0-11035_2014_CEUAS_merged_v3.nc',
                      #'/scratch/das/federico/MERGED_YEARLY_22NOV_FULL_checkDimensions/0-20001-0-11035/0-20001-0-11035_2015_CEUAS_merged_v3.nc' ]
            
            ## /scratch/das/federico/MERGED_YEARLY_22NOV_FULL_checkDimensions/0-20001-0-11035
            
            #for s in source:
                #os.system('cp ' + s + '   ' + merged_directory  )
                
            os.system('rm -r  ' + postprocessed_new )
            os.system('mkdir ' + postprocessed_new )
        
            stations_list = [ s for s in os.listdir(merged_directory) ]           
            #stations_list = [ s for s in stations_list if  '11035' in s ]           

            processed = []
            cleaned_list = []

            for file in stations_list:
                station_id = file.split("_CEUAS")[0]
                if station_id in processed:
                    print('Already processed:::' , file )
                else:
                    cleaned_list.append(file)
            
            POOL = False
            
            if POOL:
                p = Pool(30)
                func = partial(run, merged_directory, postprocessed_new, force_run)
                out = p.map(func, cleaned_list)     
            else:
                for s in cleaned_list:
                    a = run(merged_directory, postprocessed_new, force_run, s)
  