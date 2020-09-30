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
#from numba import njit
import psutil
import copy
from numba import njit
import code
import urllib.request
from functools import reduce
from tqdm import tqdm

from numba import njit
from numba.typed import Dict
from numba.core import types
# http://numba.pydata.org/numba-doc/latest/reference/pysupported.html#typed-list    --> to use dictionaries 

from collections import Counter

import cProfile
#cProfile.run('__main__')

import time as T 
t=T.time()


sys.path.append('../harvest/code')
from harvest_convert_to_netCDF_newfixes import load_cdm_tables , write_dict_h5 

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


os.system(' rm 0-20000-0-82930_CEUAS_merged_v0.nc ')
os.system(' cp 0-20000-0-82930_CEUAS_merged_v0_KEEP_NORMALMERGED.nc 0-20000-0-82930_CEUAS_merged_v0.nc ')


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

def remove_outliers(data= '', min_p= 25, max_p= 75, cut= 1, skewed= False):
    """ Finds outliers, and replace them with np.nan (to keep vector of same length)                                                                                                                                                                                              

         input ::       data = list of values 
                           min_p , max_p = minimum and maximum values of the percentile to consider to determine outliers 
                           skewed = use True to consider a skewed (not symmetricla Gaussian) distribution 
                           cut = factor to allow slight deviation from given percentiles 
         returns ::   cleaned   = list of values without outliers (filled with nans replacing outliers)                                                                                                                                                                   
                           outliers   = list of outlier values                                                                                                                                                                                                               
                           lower,upper, median = outliers delimiter and median values """

    q_min, q_max = np.nanpercentile(data, min_p), np.nanpercentile(data, max_p)
    cut_off = (q_max - q_min) * cut
    lower, upper = q_min-cut_off, q_max+cut_off

    if skewed==True:
        q50 = np.nanpercentile(data, 50)
        lower , upper = q_min-(q50-q_min)*cut ,  q_max+(q_max-q50)*cut

    median = np.nanmedian(data)
    cleaned, outliers = [],[]

    for d in np.asarray(data):
        if d >= lower and d <= upper:
            cleaned.append(d)

        else: # only storing non nans values 
            if not np.isnan(d):
                outliers.append(d)
                
    return cleaned, outliers, lower, upper, median




class MergedFile(object):
    
    def __init__(self, out_dir = '' , station_id = '' , file='' ):
        self.out_dir = out_dir 
        
        if not os.path.isdir(out_dir):
            print('+++ Creating output directory: ' , out_dir )
            os.mkdir(out_dir)
            
        self.station_id = station_id
        self.file = file 
        self.summary_file =  station_id + '_' + datetime.now().strftime("%Y-%m-%d_%H-%M-%S") + '_summary.txt'
        self.std_plevs    = [1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 85000, 92500, 100000]
        self.dic_type_attributes = np.load('../merge/dic_type_attributes.npy', allow_pickle= True).item()
        self.encodings = np.load('../merge/groups_encodings.npy' , allow_pickle = True ).item()
        
    def load(self):
        
        """ Global container of the data """
        data = {}
        data['cdm_tables'] = {}
        
        h5py_file = h5py.File(self.file, 'r+')
        
        cdm_tables = ['crs', 'observed_variable', 'sensor_configuration', 'station_configuration_codes', 'station_type', 'units', 'z_coordinate_type']
        cdm_tables = ['crs', 'observed_variable', 'station_configuration_codes', 'station_type', 'units', 'z_coordinate_type']
        
        for t in cdm_tables:
            table = xr.open_dataset (file, engine = 'h5netcdf' , group = t)
            data['cdm_tables'][t] = table
        
        data['h5py_file'] = h5py_file 
        data['station_id'] = self.station_id 
        
        data['recordtimestamp']               = xr.open_dataset (file, engine = 'h5netcdf' , decode_times = False )['recordtimestamp']
        data['recordtimestampdecoded'] = xr.open_dataset (file, engine = 'h5netcdf' , decode_times = True )['recordtimestamp']
        
        data['recordindex']         = xr.open_dataset (file, engine = 'h5netcdf' , decode_times = False )['recordindex']
        data['dateindex']            = h5py_file['dateindex']
        
        #data['sensor_id'] = h5py_file['observations_table']['sensor_id']
        data['length_max'] = len(h5py_file['observations_table']['date_time'] )

        self.data = data 
    
    
    def load_obstab_era5fb(self, era5fb=False, obs_tab=True):
        """ Loading era5fb and observations_table if necessary """

        if obs_tab:
            obs_tab = xr.open_dataset (file, engine = 'h5netcdf'      , group = 'observations_table', decode_times = True )
            #obs_tab_vars = ['z_coordinate_type' , 'z_coordinate' , 'observation_value' , 'observed_variable', 'date_time', 'sensor_id']
            self.data['obs_tab'] = obs_tab
            attr_dic = self.retrieve_attributes(obs_tab)
            self.attr_dic = attr_dic 
            
            self.fixed_columns            = ['observed_variable' , 'observation_value', 'z_coordinate' , 'z_coordinate_type', 'secondary_value' , 'value_significance' , 'original_precision', 'date_time'] # to be defined explicitly                     
            self.unavailable_columns = ['observation_id' , 'report_id' , 'source_id'] 
            self.other_columns           = [c for c in self.data['obs_tab'].keys() if c not in self.fixed_columns and c not in self.unavailable_columns] # to be calculated by averaging 
            
        if era5fb:
            era5fb_tab = xr.open_dataset (file, engine = 'h5netcdf' , group = 'era5fb', decode_times = True )
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
    
    def __init__(self, MF = ''  ):
        
        #MF = MergedFile.__init__(self, out_dir = out_dir, station_id = station_id, file= file )
        self.MergedFile = MF
        self.data = MF.data 
        
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
                print('+++ Written sensor_configuration ' )
            except  KeyError:
                print('--- Passing variable ' )
                    
        return 'done'           
                    
                    
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
                sensor_datetime[ddt] = {'sensor' : sensor, 'min_index' : 0}
                
            else:
                year, month, day, hour, minute =  int(dt[0:4]) , int(dt[4:6]) , int(dt[6:8]) , int(dt[8:10])  , int(dt[10:12] ) 
                try:
                    dt = np.datetime64( datetime ( year, month, day, hour, minute ) )
                    near = nearest(self.data['recordtimestampdecoded'].values, dt)
                    sensor_datetime[near] = {'sensor' : sensor } 
                    
                except: # this does not work because I cannot match the date_time anymore   
                    pass
                
            print (0)

        """ If I have only one entry, this instrument will be applied to the entire list of observations """
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
            while len(sensor_id) < 3:
                sensor_id = sensor_id + b' ' 
                print(sensor_id)
                
            """ If I have only one sensor, it will be applied to the whole length of the observations_table """    
            
            if  sensor_datetime[dt]['max_index'] == -1:
                length = self.data['length_max']
            else:
                length = sensor_datetime[dt]['max_index'] - sensor_datetime[dt]['min_index']
            sensor_id_array = np.full( ( length )  , sensor_id ).astype(  np.dtype('|S3')  ) 
            temp_sensor_list.append(sensor_id_array)
        
        if not temp_sensor_list:
            print('No sensor found ::: ')
            return 'NoSensor'  # stop here if no sensor is found 
        
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

        self.data['h5py_file'].close()
        
    def run(self):

        cdm_tables = self.load_cdm_tables()      
        load_Schroeder_table = self.load_Schroeder_tables()
        dummy = self.extract_sensor_id()
        status = self.replace_sensor_id()
        
        if status != 'NoSensor':
            status = sensor.write_sensorconfig()

        write = self.MergedFile.write_summary( what = 'sensor', done = True)
        
        os.system('mv  ' + self.MergedFile.file + '   ' +  self.MergedFile.out_dir )
        print(' --- Done writing the output file ' + self.MergedFile.file + '  ! ---  ' )
        return 0 

        
class MonthlyAverages(MergedFile):
    """ Class holding the utilties to produce monthly average files out of merged files """
    
    def __init__(self, MF = MergedFile ):
        self.MergedFile = MF
        self.data           = MF.data     
        
    def write_monthly_file(self , monthly_tab):
        """
        self.out_dir = out_dir 
        self.station_id = station_id
        self.file = file """
        
        out_file = self.MergedFile.out_dir + '/' + self.MergedFile.station_id + '_monthly_averages_VARIABLE.nc' 
        
        os.system('rm ' + out_file )
        """ Writing standard CDM tables """
        for t in self.data['cdm_tables'].keys():
            table = self.data['cdm_tables'][t]
            table.to_netcdf(out_file, format='netCDF4', engine='h5netcdf', mode='a', group = t)
            print('written ==== ' , t )

        """ Initialize h5py output netcdf """
        with h5py.File(out_file, 'a' ) as out_netcdf:
            out_netcdf.create_group('observations_table')
            index = np.zeros (  len(monthly_tab['date_time']) , dtype='S1')                       
            out_netcdf['observations_table'].create_dataset('index', data=index)
                        
            """ Converitng to proper data type """                   
            for variable in monthly_tab.keys(): 
                var_type = self.MergedFile.dic_type_attributes['observations_table'][variable]['type']
        
                '''
                """ trying to convert the variable types to the correct types stored as attribute, read from the numpy dic file """
                if type(self.data['monthly_obs_tab'][variable][0]) != var_type:       
                    try:
                        self.data['monthly_obs_tab'][variable] = self.data['monthly_obs_tab'][variable].astype( var_type ) 
                    except:
                        print ('FAILED converting column ' , variable, ' type ', type(self.data['monthly_obs_tab'][variable][0]) , ' to type ', var_type )
                '''
                
                try:
                    data = np.array( monthly_tab[variable] ).astype(var_type)       
                except:
                    data = np.array( monthly_tab[variable] )      
                    
                out_netcdf['observations_table'].create_dataset(variable, data.shape, data.dtype, compression=self.MergedFile.encodings['observations_table'][variable]['compression'], chunks=True)
                out_netcdf['observations_table'][variable][:]=data[:]
                
                """ Setting the attributes """
                if variable == 'date_time':
                    out_netcdf['observations_table'][variable].attrs['units'] = np.bytes_('seconds since 1900-01-01 00:00:00')            
                    
                if variable in self.MergedFile.attr_dic.keys():
                        out_netcdf['observations_table'][variable].attrs['description']       = np.bytes_(self.MergedFile.attr_dic[variable]['description'])
                        out_netcdf['observations_table'][variable].attrs['external_table' ] = np.bytes_(self.MergedFile.attr_dic[variable]['external_table'])
                        
                print('Done writing ::: ' , variable )

                """ Adding missing dimensions """
                try:
                    out_netcdf['observations_table'][variable].dims[0].attach_scale(out_netcdf['observations_table']['index'])
                except:
                    print('Failed dimensions +++++++++++++++++= ')

        """ ### Checkign if index dimension variable exists, if not create it 
        if 'index' not in self.data['h5py_file']['observations_table'].keys():
            index = np.zeros (  self.data['h5py_file']['observations_table']['date_time'] .shape[0], dtype='S1')           
            data['observations_table'].create_dataset('index', data=index)

        try:            
            data['observations_table']['sensor_id'].dims[0].attach_scale( data['observations_table']['index'] )
            data['observations_table'].create_dataset( s ,  data=stringa[:slen]  )
            data['observations_table']['sensor_id'].dims[1].attach_scale(data['observations_table']['string{}'.format(slen)])
            data['observations_table']['string{}'.format(slen)].attrs['NAME']=np.bytes_('This is a netCDF dimension but not a netCDF variable.')
            print(' *** Done with the attributes of the dimension *** ')
            
        except ValueError:
            print('Dimension already exist, passing ')
                
            """

    def get_dates_indices(self):
        """ This function analyzes the recordtimestamp and extracts the indices for the observations_table """
        timestamps = self.data['recordtimestampdecoded'].values
        indices         = self.data['recordindex'].values
        indices_plus = [i-1 for i in indices[1:]]# getting the upper index of each record  
        indices_plus.append (len(self.data['obs_tab']['date_time']) )  
        
        dates_indices = {}
        limit = 0
                
        for rts, index, index_p  in zip ( tqdm(timestamps[limit:]) , indices[limit:], indices_plus[limit:] ) :

            pd_rts = pd.DatetimeIndex([rts]) # using pandas to extract day, month, year and hour 
            
            this_hour = pd_rts[0].hour # if the observation is after 10 pm we report it to the following day. Will handle automatically case of February 
            
            if not ( this_hour >= 22 or ( this_hour >= 10 and this_hour <= 14 ) or this_hour == 0 or ( this_hour <= 2 ) ):
                #print(' Not usable hour ', this_hour)
                continue      # the distance of the hour from 00 or 12 is too large to be included in the averages 
            else:
                pass
                #print('Will keep this hour! ')
                
            if this_hour >= 22 or this_hour == 0 or this_hour <= 2 :
                pd_rts = pd_rts + timedelta(minutes=120)
                std_time = '00'
            else:
                std_time = '12'              
                         
            this_month , this_year = pd_rts[0].month , pd_rts[0].year 

            date_string = str(this_year) + '_' + str(this_month)
            if date_string not in dates_indices.keys():
                dates_indices[date_string] = {}
                
                dates_indices[date_string]['00'] = {}
                dates_indices[date_string]['12'] = {}
                
                dates_indices[date_string]['00']['indices'] = []
                dates_indices[date_string]['12']['indices'] = []
                
                dates_indices[date_string][std_time]['indices'].append( [index,index_p] )
            
            else:
                dates_indices[date_string][std_time]['indices'].append( [index,index_p] )
            
        return dates_indices


    def make_monthly_observations_table(self, variables = ''):
        """ Main utility to produce monthly averages observations_table. 
              We only consider one variable at a time and write separated netCDF files """
        
        print("LOADING variables::: "  , T.time() )
        monthly_dates = self.get_dates_indices() # independent on the variables, depends only on the records time_stamps 
        obs_var          = self.data['obs_tab']['observed_variable'].values
        z_coord_type = self.data['obs_tab']['z_coordinate_type'].values
        z_coord          = self.data['obs_tab']['z_coordinate'].values
        obs_val          = self.data['obs_tab']['observation_value'].values
        std_plevels    = np.array(self.MergedFile.std_plevs)           
        print("FONE LOADING variables::: "  , T.time() )
        
        
        """ Auxilliary functions """
        def get_mostcommon_value(indices, column, index_min='', index_max=''):
            """ Extracts the most common values (in a given moth) for the variable considered. 
                  E.g. the sensor_id of the monthly average would be the most common sensor id found in the records for that specific  month. 
                  Indices is the list of indices, of the chucnked obs_tab[index_min:index_max] where I have the observation values.
                  So index_min+indices[0] is the first available value of the desired variable. """
            
            #print("EXTRACTING values get_mostcommon_value " , T.time() )
            #values = self.data['obs_tab'][column][index_min:index_max][indices]
            
            #most_common = self.data['obs_tab'][column][index_min:index_max][indices[0]]

            #most_common = self.data['obs_tab'][column][index_min:index_max][indices[0]]
            
            ### HAD TO SIMPLY: returns the first value of the list 
            #most_common =self.data['obs_tab'][column][index_min + indices[0] ].values
            
            most_common =self.data['obs_tab'][column][index_min].values # try this. pick the fiurst value of the obs table 
            

            return most_common
            
            """
            try:
                (values,counts) = np.unique(values,return_counts=True)
                ind=np.argmax(counts)
                most_common =  values[ind]
                
                #occurence_count = Counter( list(values.values) )
                #most_common=occurence_count.most_common(1)[0][0]
            except:
                most_common = np.nan 
                
            #print(" DONE return get_mostcommon_value " , T.time() )
            
            return most_common 
            """

        def get_observation_indices(v, obs_var, z_coord_type, z_coord, obs_val, std_plevels, indices = '' ):
            """ Loop over the indices of the records and extract the specific indices of the observations, and the values of the observations """
            #results_dic = Dict.empty( key_type=types.unicode_type ,  value_type=types.float32[:],  ) ## for numba but it does not work

            results = {} # firsty solum,nc will hold variable number, second the value 
            results['is_there_data'] = False 
            for i in range(len(indices[0]) ):
                ind_min, ind_max = indices[0][i] , indices[1][i]
                                
                plevels = z_coord[ind_min: ind_max]  # available pressure levels in the records 
                                
                for p in std_plevels:
                    
                    results[p] = {}
                    results[p]['observation_values']  = []
                    results[p]['observation_indices'] = []
                    
                    if p not in list(plevels):
                        results[p]['observation_values'].append(np.nan)
                        results[p]['observation_indices'].append(np.nan )
                        
                    else:
                        unique_ind = np.where ( ( obs_var[ind_min: ind_max] == v ) & (z_coord_type[ind_min: ind_max] == 1.)   &   (z_coord[ind_min: ind_max] == p)  ) [0]         
                        if bool(unique_ind):
                            val = obs_val[ind_min: ind_max][unique_ind][0] 
                            if val < -1*10**(-18) and val > -1*10**(-15):
                                continue
                            if val > 99999999:
                                continue
    
                            results[p]['observation_values'].append(val)
                            results[p]['observation_indices'].append(unique_ind[0] )
                            results['is_there_data'] = True 
                            
                        else:
                            results[p]['observation_values'].append(np.nan)
                            results[p]['observation_indices'].append(np.nan )
                 
            return results ,  ind_min, ind_max 
                        
        def fill_columns(monthly_observations_table , tot_observations , unavailable_columns = '' , other_columns = '', index_min = '', indices = '' ):
            """ Fill the dictionary of the observations_table. when there is no observation data available (dummy empty record) """

            monthly_observations_table['z_coordinate']         .extend ( self.MergedFile.std_plevs )
            monthly_observations_table['date_time']             .extend ( [dt] * tot_observations )  # must be converted to seconds after 1900-01-01
            monthly_observations_table['z_coordinate_type'].extend ( [1] * tot_observations )
            monthly_observations_table['value_significance'].extend ( [2] * tot_observations )
            monthly_observations_table['observed_variable'].extend ( [v] * tot_observations ) 
            
            # setting nans for all the columns with unavailable data e.g. observations_id 
            for column in other_columns:
                monthly_observations_table[column].extend( [np.nan] * tot_observations )   
                    
            for column in other_columns:
                if index_min: # case where I am filling  data for an existing record 
                    value = get_mostcommon_value(indices, column, index_min= index_min , index_max = '')
                else:
                    value = np.nan
                monthly_observations_table[column].extend( [value] * tot_observations )

            return 
        
        for v in variables:
            """ Start the loop over the records """
            monthly_observations_table = {k: [] for k in self.data['obs_tab'].keys() }  # container for the averages observations_table 
            
            unavailable_columns = self.MergedFile.unavailable_columns 
            other_columns           = self.MergedFile.other_columns 
            tot_observations        = len(self.MergedFile.std_plevs)
            
            print('*** Extracting and calculating averages *** ')
            
            for m in tqdm(monthly_dates.keys() ): 
                value = monthly_dates[m] # loop over the months    
                print("\n \n MONTH IS::: " , m , '   ', T.time() )
                for time in ['00','12'] :
                    dt = datetime.strptime(m + '_15 ' + time + ':00:00', '%Y_%m_%d %H:%M:%S')   # create the date time for the 15th day of that particular month at 12 o'clock               
                    indices = value[time]['indices']         

                    if not indices: # if I have no indices, I fill with nans
                        dummy = fill_columns(monthly_observations_table , tot_observations , unavailable_columns = unavailable_columns , other_columns = other_columns  ) 

                    else:  # check if I have records within this period of time  
                                                
                        ind = np.array ( ([ i[0] for i in indices ] , [ i[1] for i in indices ]  ) , dtype = np.int32 )  # can simplify this, it is a leftover from trying numba  
                        #print('DO get_observation_indices +++ ' ,  T.time()  )                                   
                        extracted_data, index_min, index_max = get_observation_indices(v, obs_var, z_coord_type, z_coord, obs_val, std_plevels, indices = ind )
                        #print('DONE get_observation_indices +++ ' ,  T.time()  )           
                        
                        if extracted_data['is_there_data']: # check if I have at least one valid observation, otherwise do nothing 
                            
                            dummy = fill_columns(monthly_observations_table , tot_observations , unavailable_columns = unavailable_columns , other_columns = other_columns , index_min = index_min, indices = '' ) 

                            """ loop over the available std pressure levels """
                            for p in self.MergedFile.std_plevs:  # (self.MergedFile.std_plevs to keep same order ### extracted_data.keys()
                                #print("month, level, hour " , m , '  ' , p , '  ' , time , '  ' , T.time())
                                
                                values = [f for f in extracted_data[p]['observation_values'] if not np.isnan(f)  ]
                                
                                if values:
                                    values_cleaned, outliers, lower, upper, median = remove_outliers(data= np.array(values), min_p=25, max_p=75, cut= 2, skewed=False)
                                    #if m == '1979_3' and p == 7000:

                                    #print("CALCULATING ::: " , T.time() )
                                    monthly_observations_table['observation_value'].append(np.nanmean(values_cleaned) )
                                    monthly_observations_table['original_precision'] .append(np.std(values_cleaned) ) 
                                    monthly_observations_table['secondary_value']  .append( len(values_cleaned) )
                                    #print("DONE CALCULATING ::: " , T.time() )
   
                                else:
                                    monthly_observations_table['observation_value'].append( np.nan )
                                    monthly_observations_table['original_precision'] .append( np.nan ) 
                                    monthly_observations_table['secondary_value']  .append( np.nan )
                                    for column in other_columns: # other_columns: columns for which we want to keep the most freque value, e.g. latitude or sensor_id 
                                        monthly_observations_table[column].append(np.nan)


            monthly_observations_table['date_time'] =  datetime_toseconds( monthly_observations_table['date_time'] ) 
            self.data['monthly_obs_tab'] = monthly_observations_table

            
            """ Finally writing the output file """
            a = self.write_monthly_file(monthly_observations_table)
            

        
  

    
""" File source direcotry """
#merged_directory = '/raid60/scratch/federico/do/'

merged_directory = '/raid8/srvx1/federico/GitHub/CEUAS_master_MAY/CEUAS/CEUAS/public/merge/PROVA_stdplevels_only'

""" Moving postprocessed files to new directory """
postprocessed_new = '/raid60/scratch/federico/CIAONE'
os.system('mkdir ' + postprocessed_new)

""" Set if running is enforced (False: will crash at errors) """

if __name__ == '__main__':
        
            parser = argparse.ArgumentParser(description="Postprocessing Utility")
            
            parser.add_argument('--stations' , '-s', 
                                  help="List of Station Ids to Postprocess"  ,
                                  type = str,
                                  default = 'a')
        
            parser.add_argument('--instrument' , '-i', 
                                  help="Add Instrument Type"  ,
                                  type = str,
                                  default = 'False'  )
            
            parser.add_argument('--force_run' , '-f', 
                                  help="Force running the file(s)"  ,
                                  type = str,
                                  default = 'False' )
            
            parser.add_argument('--desrozier' , '-d', 
                                  help="Calculate Desroziers' statistics"  ,
                                  type = str,
                                  default = 'False'  )
            
            parser.add_argument('--monthly_averages' , '-m', 
                                  help="Calculate monthly averages"  ,
                                  type = str,
                                  default = 'False'  )    
            
            args = parser.parse_args()
            
            stations                   = args.stations
            get_instrument       = args.instrument
            force_run                = args.force_run
            get_desrozier          = args.desrozier
            monthly_averages  = args.monthly_averages
            
            # 0-20000-0-70316_CEUAS_merged_v0.nc   ### check if this file works now 
            #file = '0-20000-0-70316_CEUAS_merged_v0.nc' 
            
            stations_list = os.listdir(merged_directory)                        

            for s in stations_list:
                file = merged_directory + '/' + s
                station_id = file.split("_CEUAS")[0].split(merged_directory)[1].replace('/','').split('-')[-1]

                print (' I will process the file ::: ' , file , ' ::: station_id ::: ' , station_id )  
                
                """ Initialize classes """
                MF = MergedFile(out_dir = postprocessed_new , station_id = station_id , file = file  )
                """ Read data """
                data = MF.load()
                
                if force_run in ['yes', 'y', 'YES', 'Y', 'True', 'true']: 
                    print("    === Running in force mode ===     ")
                    
                    try:                        
                        """ Running sensor module """
                        if get_instrument in  ['yes', 'y', 'YES', 'Y', 'True', 'true'] :                                       
                            sensor = Sensor( MF = MF  )  # loading sensor class 
                            run = sensor.run() # running module 
                            
                        if get_desrozier in ['yes', 'y', 'YES', 'Y', 'True', 'true'] :                       
                            tabs = MF.load_obstab_era5fb(obs_tab=True, era5fb=True)
                            print('Running Desroziers statistics')
                            
                        if monthly_averages in ['yes', 'y', 'YES', 'Y', 'True', 'true'] :                       
                            tabs = MF.load_obstab_era5fb(obs_tab=True, era5fb=False)
                            print('Extracting Monthly Averages')                          
                            
                    except:
                        print(" The file " + file + ' hase failed! MUST REDO! ')
                    
                else:
                    print('   +++ Running in normal mode +++')
                    """ Running sensor module """
                    if get_instrument in  ['yes', 'y', 'YES', 'Y', 'True', 'true']:
                        sensor = Sensor( MF = MF )  # loading sensor class 
                        run = sensor.run() # running module 
                        
                    if get_desrozier in ['yes', 'y', 'YES', 'Y', 'True', 'true'] :
                        tabs = MF.load_obstab_era5fb()
                        print('Running Desroziers statistics')     
                            
                    if monthly_averages in ['yes', 'y', 'YES', 'Y', 'True', 'true'] :                       
                        tabs = MF.load_obstab_era5fb(obs_tab=True, era5fb=False)
                        monthly = MonthlyAverages(MF)
                        a = monthly.make_monthly_observations_table( variables = [85])
                        
                        print(' *** Done Monthly Averages for file ' , file )    
