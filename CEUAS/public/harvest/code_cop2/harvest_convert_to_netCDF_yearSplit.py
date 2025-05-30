#!/usr/bin/env python
import sys
import os.path 
import glob 
import psutil
import subprocess
import urllib.request
import xarray as xr
import h5py
from datetime import date, datetime,timedelta
import time
from multiprocessing import Pool
#from netCDF4 import Dataset
import gzip
import pandas as pd   
pd.options.mode.chained_assignment = None

import zipfile
from functools import partial
from numba import njit
import argparse
from io import StringIO
import numpy as np
import warnings
import numpy
from numba import njit
from eccodes import *
from tqdm import tqdm 
# from harvester_yearsplit_parameters import *

pv=sys.version.split('.')
# if pv[1]<'8':
    # from eccodes import *

import ssl
ssl._create_default_https_context = ssl._create_unverified_context


warnings.simplefilter(action='ignore', category=FutureWarning) # deactivates Pandas warnings 
from warnings import filterwarnings
filterwarnings(action='ignore', category=DeprecationWarning, message='`np.int` is a deprecated alias')

# pd.set_option('display.max_columns', None)
# pd.set_option('display.max_rows', None)
# pd.set_option('display.width', None)
# pd.set_option('display.max_colwidth', -1)

debug = False

try:
    os.mkdir('logs')
except:
    pass


""" Some colors for pretty printout """ 
red    = '\033[91m' 
cend   = '\033[0m'
blue   = '\033[34m'
green  = '\033[92m'
yellow = '\033[33m'

long_string_len = 100  # maximum allowed length of strings in header and observations_table
fixed_string_len = 20  # maximum allowed length of strings in header and observations_table
id_string_length = 10 # maximum allowed length of strings for observation_id and report_id in header and observations_table


### nan for integer  number 
int_void = -2147483648


""" Possible variable types as listed in the CMD tables """
okinds={'varchar (pk)':np.dtype('|S' + str(fixed_string_len) ),
        'varchar':np.dtype('|S' + str(fixed_string_len) ), 
               'numeric':np.float32, 
               'int':np.int32,
               #'timestamp with timezone':np.datetime64,
               'timestamp with timezone':np.int64,
               'int[]*':list,
               'int[]':list,
               'varchar[]*':list,
               'varchar[]':list}

""" Variable types to be used in the compressed netCDF files """
kinds={'varchar (pk)':str,
       'varchar(pk)':str,
             'varchar ':str,
             'varchar':str,
             ' varchar':str,             
             'numeric':np.float32,
             'numeric ':np.float32,             
             'int':np.int32,
             'int ':np.int32,             
             'int(pk)' : np.int32,
             'timestamp with timezone':np.float32,             
             'int[]*':list,
             'int[]':list,
             'varchar[]*':list,
             'varchar[]':list}

gkinds={'varchar (pk)':numpy.dtype('|S'+ str(fixed_string_len)  ),
        'varchar':numpy.dtype('|S'+ str(fixed_string_len)  ),
               'numeric':numpy.float32,

               'int':numpy.int32,
               'timestamp with timezone':numpy.datetime64,
               'int[]*':numpy.int32,
               'int[]':numpy.int32,
               'varchar[]*':numpy.dtype('|S'+ str(long_string_len)  ),
               'varchar[]':numpy.dtype('|S'+ str(long_string_len)  )}

@njit
def dfdateinput(dvari,tvari,dhash,dsecs,dt):
    j=0
    for i in range(dvari.shape[0]):
        while dhash[j]<dvari[i]:
            j+=1
        dt[i]=dsecs[j]+(tvari[i]//10000)*3600+((tvari[i]%10000)//100)*60+tvari[i]%100

    return

def make_datetime(dvar,tvar):
    """ Converts into date-time """
    datelist = pd.date_range(pd.Timestamp(year=1900,month=1,day=1), periods=50000) #### !!!
    dhash=numpy.array(datelist.year)*10000+numpy.array(datelist.month)*100+numpy.array(datelist.day)
    dt=numpy.empty_like(dvar,dtype=numpy.int64)
    dsecs=((datelist.values-datelist.values[0])//1000000000).astype(numpy.int64)

    #df=pd.DataFrame({'year':dvari//10000,'month':(dvari%10000)//100,'day':dvari%100,
                        #'hour':tvari//10000,'minute':(tvari%10000)//100,'second':tvari%100})

    dfdateinput(dvar.values,tvar.values,dhash,dsecs,dt)

    return dt

@njit(cache=True)
def narange(x,y):
#    y=numpy.arange(48,58,dtype=numpy.int32)
    #y=numpy.array('0123456789').view('S1')
    l=0
    i=numpy.arange(10,dtype=numpy.int32)
    for i[9] in range(10):
        for i[8] in range(10):
            for i[7] in range(10):
                for i[6] in range(10):
                    for i[5] in range(10):
                        for i[4] in range(10):
                            for i[3] in range(10):
                                for i[2] in range(10):
                                    for i[1] in range(10):
                                        for i[0] in range(10):
                                            for j in range(10):
                                                if i[j]>0:
                                                    x[l][9-j]=y[i[j]]

                                            l+=1
                                            if l==x.shape[0]:
                                                return

    return

def make_obsid(ivar):
    x=numpy.zeros(ivar.shape[0],dtype='S'+ str(id_string_length))
    x.fill(b'0000000000')
    y=numpy.array(['0','1','2','3','4','5','6','7','8','9'],dtype='S1')

    narange(x.view(dtype='S1').reshape(x.shape[0],10),y)
    #x=numpy.char.zfill(numpy.arange(ivar.shape[0]).astype('S'+ str(id_string_length)), id_string_length  )    #slow old version
    return x

def make_recid(ivar):
    x=numpy.char.zfill(numpy.arange(ivar.values.shape[0]).astype('S'+ str(id_string_length)) , id_string_length )
    return x

@njit(cache=True)
def nrecfill(y,x,ivv):
    for i in range(x.shape[0]-1):
        for j in range(ivv[i],ivv[i+1]):
            y[j]=x[i]
    if x.shape[0]>1:
        for j in range(ivv[i+1],y.shape[0]):
            y[j]=x[i+1]
    return

def make_obsrecid(fbvar,ivar):
    x=numpy.char.zfill(numpy.arange(ivar.values.shape[0]).astype('S' + str(id_string_length )) , id_string_length )
    y=numpy.zeros(fbvar.shape[0]).astype('S' + str(id_string_length ) )
    y.fill(b'0000000000')
    nrecfill(y,x,ivar.values)
    #for i in range(ivar.values.shape[0]-1):
        #y[ivar.values[i]:ivar.values[i+1]]=x[i]
    #if ivar.values.shape[0]>1:
        #y[ivar.values[i+1]:]=x[i+1]
    return y


def make_vars(ivar):
    tr=numpy.zeros(300,dtype=int) 
    """ translates odb variables number to Lot3 numbering convention 
    CDM
    57: pressure of air column at a specified height
    117: geopotential height
    106,107: wind speed and wind direction
    139,140: eastward and northward wind speed from profile measurement
    126: air temperature from profile measurement
    137,138: air dew point and air relative humidity from profile measurement
    39: specific humidity
    """

    tr[1]=117  # should change
    tr[2]=126
    tr[3]=139
    tr[4]=140
    tr[7]=39 #spec hum
    tr[29]=138 #relative hum
    tr[59]=137 # dew point
    tr[111]=106 #dd
    tr[112]=107  #ff
    #
    tr[39]= 126 # 2m T according to proposed CDM standard  # bug fixed on Feb 08 2022 (136->126)
    tr[40]= 137 # 2m Td according to proposed CDM standard
    tr[41]= 139 #10m U according to proposed CDM standard
    tr[42]= 140  #10m V according to proposed CDM standard
    tr[58]= 138 # 2m rel hum according to proposed CDM standard

    x=tr[ivar.astype(int)] # reads the varno from the odb feedback and writes it into the variable id of the cdm
    return x


def make_units(ivar):
    tr=numpy.zeros(113,dtype='int32')
    tr[1]=631  # geopotential gpm
    tr[2]=5   #temperature K
    tr[3]=731  #uwind m/s
    tr[4]=731  #vwind
    tr[7]=622 #spec hum kg/kg
    tr[29]=0 #relative hum e/es Pa/Pa
    tr[59]=5 # dew point
    tr[111]=110 #dd
    tr[112]=731  #ff

    tr[39]= 5 # 2m T
    tr[40]= 5 # 2m Td
    tr[41]= 731 #10m U
    tr[42]= 731  #10m V
    tr[58]= 0 # 2m rel hum

    x=tr[ivar.astype(int)] # reads the varno from the odb feedback and writes it into the variable id of the cdm
    return x


""" Translates the odb variables name  into cdm var """
cdmfb={'observation_value':'obsvalue@body',
       'observed_variable':[make_vars,'varno@body'],
       'observation_id':[make_obsid,'date@hdr'],
       'observations_table.report_id':[make_obsrecid,'date@hdr','recordindex'],
       'header_table.report_id':[make_recid,'recordindex'],
       'z_coordinate_type':'vertco_type@body',
       'z_coordinate':'vertco_reference_1@body',
       'date_time':[make_datetime,'date@hdr','time@hdr'],
       'report_timestamp':[make_datetime,'date@hdr','time@hdr'],
       'record_timestamp':[make_datetime,'date@hdr','time@hdr'],
       'longitude':'lon@hdr',
       'latitude':'lat@hdr',
       'source_id': 'source_id',       
       'units':[make_units,'varno@body'],
       'primary_station_id':'statid@hdr',       
       }


# dic containig all the variables which are present in the header or obervations table, but not in the era5fb table 

cdmfb_noodb={'observation_value':'obsvalue@body',
             'observed_variable':'varno@body',
                          'z_coordinate_type':'vertco_type@body',
                          'z_coordinate':'vertco_reference_1@body',
                          'date_time': 'report_timestamp', 
                          'report_timestamp': 'report_timestamp' , # only available for igra2 
                          'longitude':'lon@hdr',
                          'latitude':'lat@hdr',
                          'observation_id':'observation_id' ,
                          'source_file':'source_file',
                          #'product_code': 'product_code' ,
                          'report_meaning_of_timestamp': 'report_meaning_of_timestamp',
                          'report_id':'report_id' ,
                          'number_of_pressure_levels' : 'num_lev',  # this entry does not exist in CDM header_table 
                          'units' : 'units',
                          'source_id': 'source_id',
                          'primary_id' : 'statid@hdr' ,                  # station_configuration
                          'primary_station_id':'statid@hdr',
                          'sensor_id':'sensor_id',
                          'observation_height_above_station_surface':'observation_height_above_station_surface'}  # header_table 


## NB CONVENTIONS FOR TIMESTAMPS
### chaged and fixed in November 2023 
# report_timestamp = exact release time, which will be copied to all observations (if available)
# will be flagged with the report_meaning_of_timestamp = 1 if this is given, otherwise leave as nan if not given 
# and the report_timestamp is simply the copy of the record_timestamp
# record_timestamp = nominal time 

#NB::: only IGRA@ and NPSOUND data have the release_time information that might differ from nominal time
# Some intercomparison eg Mauritius and Mauritius digitized have the release time and not nominal time

### NB CONVENTION FOR z_coordinate
# according to the official glamod table, height is in meters with a code 0 
# https://github.com/glamod/common_data_model/blob/master/tables/z_coordinate_type.dat
# We use:
# 0 -> height in meters
# 1 -> Pressure [Pa]
# 2 -> geopotential i.e. height * 9.8...


def check_read_file(file='', read= False):

    if not os.path.isfile(file):
        raise IOError("File not Found! %s" % file)

    if read:
        if '.zip' in file:
            archive = zipfile.ZipFile(file, 'r')
            inside = archive.namelist()
            tmp = archive.open(inside[0])
            tmp = io.TextIOWrapper(tmp, encoding='utf-8')
            tmp = tmp.read()
            archive.close()
            data = tmp.splitlines()  # Memory (faster) 

        elif '.gz' in file:
            with gzip.open(file, 'rt',  encoding='utf-8') as infile:
                tmp = infile.read()  # alternative readlines (slower)                                                                                                                                                                                                            
                data = tmp.splitlines()  # Memory (faster)                                                                                                                                                                                                                              
        else:
            with open(file, 'rt') as infile:
                tmp = infile.read()  # alternative readlines (slower)                                                                                                                                                                                                                   
                data = tmp.splitlines()  # Memory (faster)  

        return data




# https://github.com/glamod/common_data_model/blob/master/tables/observed_variable.dat
# https://github.com/glamod/common_data_model/blob/master/tables/units.dat
# https://apps.ecmwf.int/odbgov/varno/

""" Dictionary mapping names, odb_codes and cdm_codes . """ 
cdmvar_dic = {'temperature'          : { 'odb_var': 2      , 'cdm_unit': 5        , 'cdm_var': 126}     ,        # K, Air temperature (from profile measurement) 
                'wind_direction'      : { 'odb_var': 111   , 'cdm_unit': 110    , 'cdm_var': 106}   ,            # degree (angle) direction from which the wind is blowing Lot 1 uses dd  - WMO abbrev.
                'wind_speed'           : { 'odb_var': 112  , 'cdm_unit': 731     , 'cdm_var': 107 } ,            # m/s , Speed is the magnitude of velocity.
                'uwind'                    : { 'odb_var': 3       , 'cdm_unit': 731     , 'cdm_var': 139}   ,    # m/s, Eastward wind speed (from profile measurement)
                'vwind'                    : { 'odb_var': 4       , 'cdm_unit': 731     , 'cdm_var': 140}    ,   # m/s, Northward wind speed (from profile measurement)
                'dew_point'             : { 'odb_var': 59             , 'cdm_unit': 5  , 'cdm_var': 137}     ,   # K, Dewpoint measurement (from profile measurement) 
                'dew_point_depression' : { 'odb_var': 299  , 'cdm_unit': 5     , 'cdm_var': 34}   ,              # K fake number, does not exhist in ODB file 
                'relative_humidity'  : { 'odb_var': 29     , 'cdm_unit': 0    , 'cdm_var': 138}     ,            # Relative humidity (from profile measurement)
                'gph'                       : { 'odb_var': 1       , 'cdm_unit': 631    , 'cdm_var': 117} , # geopotential , it iw wrong, it shoudl be m2/s2 and not the code 631
                'pressure'               : { 'odb_var': 999    , 'cdm_unit': 32       , 'cdm_var': 57}      ,    # Pa  (it goes into z_coordinate type, does not exist in odb files )
            }

""" Common columns of the read dataframes (not from ODB files, which have their own fixed column names definitions) """
column_names = [ 'source_id', 'report_id',  'observation_id', 'record_timestamp' , 'iday', 'statid@hdr', 'lat@hdr', 'lon@hdr', 'vertco_reference_1@body',
                 'obsvalue@body', 'varno@body' , 'units',  'number_of_pressure_levels', 'vertco_type@body']




def bufr_to_dataframe(file=''):
    """ Read a bufr file and convert to a Pandas DataFrame 
        Variables used inside the DataFrame are already CDM compliant                                                                                                                                                                                                               
        Args:                                                                                                                                                                                                                                                                       
             file (str): path to the bufr  file  

        Returns:                                                                                                                                                                                                                                                                    
             Pandas DataFrame with cdm compliant column names        
    """

    if debug:
        print("Running bufr_to_dataframe for: ", file)

    check_read_file (file = file, read= False)
    f = open(file)
    #source_file = [l for l in file.split('/') if '.bfr' in l][0]
    read_data = []

    """ Name of the columns as they will appear in the pandas dataframe (not necessarily CDM compliant) """
    #column_names = ['report_timestamp' , 'iday',  'station_id', 'latitude', 'longitude', 'pressure', 'value','varno@body']

    lat, lon, alt, blockNumber, stationNumber, statid = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

    obs_id, report_id = -1, 0 # progressive observation id
    stations_id = [] 

    while 1:
        #lista = [] # temporary list
        bufr = codes_bufr_new_from_file(f)

        if bufr is None:
            break

        codes_set(bufr, 'unpack', 1) # eCcodes must expand all the descriptors and unpack the data section

        if 'AMMA' in file or 'MAESTRO':
            date = codes_get_array(bufr, "typicalDate")[0]
        else:
            date = '19'+codes_get_array(bufr, "typicalDate")[0][2:]  # must fix this cause when reading from e.g. era5/odbs/ai_bfr/era5.10106.bfr , I obtain codes_get_array(bufr, "typicalDate") -> 20710627 i.e. it starts from 2000 and not 1900 !
        timePeriod = codes_get_array(bufr, "typicalTime")[0]   

        year, month, day =  date[0:4], date[4:6] , date[6:8]
        hour, minutes = timePeriod[0:2] , timePeriod[2:4]

        idate =  datetime.strptime(year + month + day + hour + minutes, '%Y%m%d%H%M')
        iday = int(year + month + day )
        try:

            pressure          = codes_get_array(bufr, "pressure") 
        except:
            continue # next iteration if pressure not available


        wind_direction = codes_get_array(bufr, "windDirection")
        wind_speed     = codes_get_array(bufr, "windSpeed")

        try:  # not all the bufr files have the dewpoint 
            temperature    = codes_get_array(bufr, "airTemperature")           
        except:
            temperature= np.empty((len(pressure)))
            temperature[:] = np.nan

        try:  # not all the bufr files have the dewpoint 
            dew_point          = codes_get_array(bufr, "dewpointTemperature")
        except:
            #dew_point= np.empty((1, len(pressure)))  # I do not understnad why this works i.e. (1, len()) gives an array of array
            dew_point= np.empty( len(pressure))
            dew_point[:] = np.nan

        num_lev             = len(pressure) # number of  distinct pressure levels 

        try:
            geopotential   = codes_get_array(bufr, "nonCoordinateGeopotentialHeight")         
        except:
            geopotential = np.full( (1,len(temperature)) , np.nan )[0,:]

        if report_id == 0:
            ''' Check again but these values should remain the same for all cnt, so it makes no sense to read them every time '''
            try:
                lat                     = codes_get(bufr, "latitude")
                lon                    = codes_get(bufr, "longitude")
            except:
                continue
            try:
                alt                     = float(codes_get(bufr, "heightOfStation"))
            except:
                alt = 99999
            try:
                blockNumber    = codes_get(bufr, "blockNumber")
                stationNumber = codes_get(bufr, "stationNumber")
                #statid                = str(blockNumber*1000+stationNumber) # changed to int instead of str
                statid                = blockNumber*1000+stationNumber
                if statid not in     stations_id:
                    stations_id.append(statid) 
            except:
                statid = '99999'
                if statid not in     stations_id:
                    stations_id.append(statid)                 

        codes_release(bufr)

        miss_value = -1.e100     

        for i in range(len(temperature)):
            obs_id = obs_id + 1 
            airT         = temperature[i]
            winds      = wind_speed[i]
            windd      = wind_direction[i]
            press       = pressure[i]
            gph         =  geopotential[i]
            dp = dew_point[i]
            if press == miss_value:
                press = np.nan 
            if dp == miss_value:
                dp = np.nan
            if airT == miss_value :    # replacing none values with numpy nans
                airT = np.nan 
            if winds == miss_value:
                winds = np.nan
            if gph == miss_value:
                gph = np.nan                
            if windd == 2147483647 or windd == -2147483647:
                windd = np.nan 

            for value,var in zip( [gph, airT, winds, windd, dp],  ['gph', 'temperature', 'wind_speed', 'wind_direction', 'dew_point'] ):
                obs_id = obs_id + 1 
                if not np.isnan(press):     # when pressure is available, z_coord== pressure and z_type==1
                    z_type = 1                      
                    read_data.append( ( 'BUFR'.rjust(10), report_id, int(obs_id),  idate, iday, statid, lat, lon, press, value, cdmvar_dic[var]['cdm_var'] , int(cdmvar_dic[var]['cdm_unit']) , num_lev , z_type  ) ) 
                elif  (np.isnan(press) and  not np.isnan(gph) ) :  # when pressure is not available, z_coord== gph and z_type==2 
                    z_type = 2              
                    read_data.append( ( 'BUFR'.rjust(10), report_id, int(obs_id),  idate, iday, statid, lat, lon, gph, value, cdmvar_dic[var]['cdm_var'] , int(cdmvar_dic[var]['cdm_unit']) , num_lev , z_type  ) ) 
                else:
                    z_type = -2147483648              
                    read_data.append( ( 'BUFR'.rjust(10), report_id, int(obs_id),  idate, iday, statid, lat, lon, press, value, cdmvar_dic[var]['cdm_var'] , int(cdmvar_dic[var]['cdm_unit']) , num_lev , z_type  ) ) 

        report_id += 1

    df = pd.DataFrame(data= read_data, columns= column_names)       

    df['observation_id']  = np.chararray.zfill( (df['observation_id'].astype(int)) .astype('S'+str(id_string_length ) ), id_string_length  )  #converting to fixed length bite objects 
    df['report_id']           = np.chararray.zfill( (df['report_id'].astype(int)).astype ('S'+str(id_string_length ) ), id_string_length  )

    col = [c for c in df.columns if '_id' not in c]
    df_new = df[col].replace([-999.9, -9999, -999, -999.0, -99999.0, -99999.9, 99999.0,  -99999.00 ], np.nan)    
    for c in ['observation_id', 'report_id']:
        df_new[c] = df[c]

    df_new['report_timestamp'] = df_new['record_timestamp'] 

    #df_new = df_new.sort_values(by = ['record_timestamp', 'vertco_reference_1@body' ] )    

    return df_new, stations_id



def uadb_ascii_to_dataframe(file=''):
    """ Read an uadb stationfile in ASCII format and convert to a Pandas DataFrame.                                                                                                                                                                                          
        Adapted from https://github.com/MBlaschek/CEUAS/tree/master/CEUAS/data/igra/read.py                                                                                                                                                                                         
        Variables used inside the DataFrame are already CDM compliant   
        Documentation available at: https://rda.ucar.edu/datasets/ds370.1/docs/uadb-format-ascii.pdf

        Args:
             file (str): path to the uadb station file

        Returns:
             Pandas DataFrame with cdm compliant column names            
    """     

    if debug:
        print("Running uadb_ascii_to_dataframe for: ", file)    

    data = check_read_file(file=file, read=True) 

    nmiss = 0
    search_h = False   
    read_data = []

    usi,idate, usi, lat, lon, lat, stype, press, gph, temp, rh, wdir, wspd = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

    #usi,idate, usi, lat, lon, lat, stype, press, gph, temp, rh, wdir, wspd, iday, ident, numlev= 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    obs_id = 0
    stations_id = [] 

    ### Complication
    ### NCAR dataset is not well documented and the quality of data is quite low.
    ### In particular they provide values for pressure which are large outliers, and they sometimes express values of wind speed in knots instead of m/s for standard levels.
    ### See dedicated notebook for the example of the station VIENNA, uadb_trhc_11035
    ### record: 1987  3  25  600
    ### Rule:
    ### - will keep only record with ltyp = 4,3 or 5
    ### - if the report type flag contains the values [3,4,5] : keep only data with flag 3 or 5 
    ### - if the report type flag contains only the value 4 and not 3,5: keep only data with flag 4 
    
    
    ### 
    temp_read_data = []
    ltyp_flag_all=[]
    
    
    for i, line in enumerate(data):
        if line[0] == 'H':
            
            if len(temp_read_data) >0:
                ltyp_flags = np.unique(ltyp_flag_all)
                if 4 in ltyp_flags and 3 not in ltyp_flags and 5 not in ltyp_flags: # keep 2,4
                    read_data.extend(temp_read_data) # keep 2,3,5

                elif 4 in ltyp_flags and (3 in ltyp_flags or 5 in ltyp_flags): # keep 4,5 
                    ind_to_keep = list( np.where (np.array(ltyp_flag_all) !=4)[0] )
                    for ind in ind_to_keep:
                        read_data.append(temp_read_data[ind])
                elif 4 not in ltyp_flags:
                    read_data.extend(temp_read_data) # keep 2,3,5
                    
                temp_read_data = []
                ltyp_flag_all=[]                
                
            else:
                pass
            
            try:
                # Header
                usi      = int(line[2:14])  # unique station identifier

                ident    = int(line[15:21].replace(' ',''))# WMO
                if ident not in stations_id:
                    stations_id.append(ident)

                #if len(ident) == 4:
                #    ident = '0' + ident 
                #idflag   = int(line[22:24])  # id flag
                #d_src    = int(line[25:28])  # source dataset
                #version  = float(line[29:34])  # version
                #dateflag = int(line[35:37])  # date flag
                year     = line[38:42]  # year
                month    = "%02d" % int(line[43:45])
                day      = "%02d"  % int(line[46:48])
                hour     = line[49:53]
                #locflag  = int(line[54:56])  # Location Flag
                lat      = float(line[57:67])

                lon      = float(line[68:78])
                if lon > 180:
                    lon = - (360-lon)

                #ele      = float(line[79:85])
                #stype    = int(line[86:88])
                numlev   = int(line[89:93])
                #pvers    = line[94:102]

                if '99' in hour:
                    hour = hour.replace('99', '00')

                if '99' in day:
                    search_h = True
                    continue

                minutes = int(hour) % 100                
                hour = "%02d" % (int(hour) // 100)
                if minutes > 60 or minutes < 0:
                    minutes = 0
                minutes = "%02d" % minutes
                idate = datetime.strptime(year + month + day + hour + minutes, '%Y%m%d%H%M')
                iday = int(year + month + day)
                #pday = int(day)
                search_h = False

            except Exception as e:
                search_h = True

        elif search_h:
            nmiss += 1
            continue  # Skipping block

        else:
            ltyp      = int(line[0:4])
            
            if ltyp not in [1,2,3,4,5]: # now includes "1" for the surface observation
                continue 
            
            ltyp_flag_all.append(ltyp)
            p   = float(line[5:13])

            if p != -99999.0 and p != 9999.9: 
                press   = float(line[5:13])*100  # converting to Pa, since P is given in mb (1 mb = 1 hPa) 
                if press > 110000:  # remove clearly wrong outliers 
                    press = np.nan 
            else:
                press = np.nan                 

            gph = float(line[14:22]) # gph [m]

            if gph == -999.0 or gph == -99999.00 or gph >= 99999.0:
                gph = np.nan
            else:
                gph = gph * 9.80665  # gph is expressed in meters 

            temp = float(line[23:29])
            if temp == -999.0:
                temp = np.nan 
            else:
                temp = temp + 273.15

            rh  = float(line[30:36])  # %
            if rh == -999.0:
                rh = np.nan
            else:
                rh = rh / 100.  # convert to absolute ratio  

            wdir    = float(line[37:43])            
            if wdir == -999.0 or wdir == -999 :
                wdir = np.nan

            wspd   = float(line[44:50])  # [m/s], module of the velocity
            if wspd <0 :
                wspd = np.nan     

            if ltyp ==2: # correcting wrong speed in knots, remove also wind direction
                wspd, wdir = np.nan, np.nan 
                
                
            try:

                for value,var in zip([ gph, temp, wspd, wdir, rh],  [ 'gph', 'temperature', 'wind_speed', 'wind_direction', 'relative_humidity'] ):
                    obs_id = obs_id +1
                    if not np.isnan(press):     # when pressure is available, z_coord== pressure and z_type==1
                        z_type = 1                    
                        temp_read_data.append( ( 'NCAR'.rjust(10), int(usi), int(obs_id), idate, iday, ident, lat, lon, press, value, cdmvar_dic[var]['cdm_var'] , int(cdmvar_dic[var]['cdm_unit']), numlev , z_type, ltyp) )
                    elif  (np.isnan(press) and  not np.isnan(gph) ) :  # when pressure is not available, z_coord== gph and z_type==2 
                        z_type = 2    # geopotential = 2          
                        temp_read_data.append( ( 'NCAR'.rjust(10), int(usi), int(obs_id), idate, iday, ident, lat, lon, gph, value, cdmvar_dic[var]['cdm_var'] , int(cdmvar_dic[var]['cdm_unit']), numlev , z_type, ltyp) )
                    else:
                        z_type = -2147483648             
                        temp_read_data.append( ( 'NCAR'.rjust(10), int(usi), int(obs_id), idate, iday, ident, lat, lon, press, value, cdmvar_dic[var]['cdm_var'] , int(cdmvar_dic[var]['cdm_unit']), numlev , z_type, ltyp) )

                '''
                ltyp_flags = np.unique(ltyp_flag_all)
                if 4 in ltyp_flags and 3 not in ltyp_flags and 5 not in ltyp_flags: # keep 2,4
                    read_data.extend(temp_read_data) # keep 2,3,5

                elif 4 in ltyp_flags and 3 in ltyp_flags or 5 in ltyp_flags: # keep 4,5 
                    ind_to_keep = np.where(ltyp_flag_all !=4)
                    for ind in ind_to_keep:
                        read_data.append(temp_read_data[ind])
                elif 4 not in ltyp_flags:
                    read_data.extend(temp_read_data) # keep 2,3,5
                '''
                
            except:
                print('WRONG NCAR DATA LINE **** ')
                0

    column_names = [ 'source_id', 'report_id',  'observation_id', 'record_timestamp' , 'iday', 'statid@hdr', 'lat@hdr', 'lon@hdr', 'vertco_reference_1@body',
                 'obsvalue@body', 'varno@body' , 'units',  'number_of_pressure_levels', 'vertco_type@body', 'observation_height_above_station_surface']

    #column_names = ['source_file', 'product_code', 'report_id', 'observation_id', 'report_timestamp' , 'iday', 'station_id', 'lat@hdr', 'lon@hdr', 'vertco_reference_1@body', 'obsvalue@body', 'varno@body' ,  'units',  'number_of_pressure_levels' ]

    df = pd.DataFrame(data= read_data, columns= column_names)   

    # save ohass only as 1 for surface, or nan
    df.observation_height_above_station_surface = np.where(df.observation_height_above_station_surface != 1, np.nan, df.observation_height_above_station_surface)

    df['observation_id']  = np.chararray.zfill( (df['observation_id'].astype(int)) .astype('S'+str(id_string_length ) ), id_string_length  )  #converting to fixed length bite objects 
    df['report_id']           = np.chararray.zfill( (df['report_id'].astype(int)).astype ('S'+str(id_string_length ) ), id_string_length  )

    col = [c for c in df.columns if '_id' not in c]
    df_new = df[col].replace([-999.9, -9999, -999, -999.0, -99999.0, -99999.9, 99999.0,  -99999.00 ], np.nan)    
    for c in ['observation_id', 'report_id']:
        df_new[c] = df[c]


    df_new['report_timestamp'] = df_new['record_timestamp'] 

    print('Done reading DF')
    return df_new , stations_id 


def read_giub(file=''):
    """ Read giub 2.1 station files saved in csv format.

        Args:
             file (str): path to the giub station file

        Returns:
             Pandas DataFrame with cdm compliant column names
    """    

    # column names standardized:  ['source_id', 'report_id', 'observation_id', 'record_timestamp', 'iday', 'statid@hdr', 'lat@hdr', 'lon@hdr', 
    # 'vertco_reference_1@body', 'obsvalue@body', 'varno@body', 'units', 'number_of_pressure_levels', 'vertco_type@body']

    #col = ['index', 'statid', 'latitude', 'longitude', 'pressure',
    #       'airTemperature', 'nonCoordinateGeopotentialHeight', 'typicalDate',
    #       'typicalTime', 'windSpeed', 'windDirection', 'dewpointTemperature',
    #       'heightOfStation']    

    df = pd.read_csv( file, sep = '\t').astype(str)
    df['rh'] = df['rh'].astype(float).apply(lambda x: x / 100 if x != -999 else x).astype(str)

    names = ["index"  , "date"   , "time"   , "temp"  ,  "gph"   , 'wspeed'  ,"wdir" ,   'uwind'   ,"vwind"  , "rh"    ,  "sh"    ,  'z_coordinate'  ,  'z_coordinate_type'       ,'days'  ,  'months'  ,"geopotential" ]

    statid = file.split('/')[-1].split('_')[0].replace('.txt','')

    all_standard_variables =  ['source_id', 'report_id', 'observation_id', 'record_timestamp', 'iday', 'statid@hdr', 'lat@hdr', 'lon@hdr', 
                               'vertco_reference_1@body', 'obsvalue@body', 'varno@body', 'units', 'number_of_pressure_levels', 'vertco_type@body'] 

    all_data = {}
    for v in all_standard_variables:
        all_data[v] = []

    ''' # just a reminder 
    cdmvar_dic = {'temperature'          : { 'odb_var': 2      , 'cdm_unit': 5        , 'cdm_var': 126}     ,  # K, Air temperature (from profile measurement) 
                             'wind_direction'      : { 'odb_var': 111   , 'cdm_unit': 110    , 'cdm_var': 106}   ,  # degree (angle) direction from which the wind is blowing Lot 1 uses dd  - WMO abbrev.
                             'wind_speed'           : { 'odb_var': 112  , 'cdm_unit': 731     , 'cdm_var': 107 } ,  # m/s , Speed is the magnitude of velocity.
                             'uwind'                    : { 'odb_var': 3       , 'cdm_unit': 731     , 'cdm_var': 139}   ,  # m/s, Eastward wind speed (from profile measurement)
                             'vwind'                    : { 'odb_var': 4       , 'cdm_unit': 731     , 'cdm_var': 140}    ,  # m/s, Northward wind speed (from profile measurement)
                             'dew_point'             : { 'odb_var': 59             , 'cdm_unit': 5  , 'cdm_var': 137}     ,  # K, Dewpoint measurement (from profile measurement) 
                             'dew_point_depression' : { 'odb_var': 299  , 'cdm_unit': 5     , 'cdm_var': 34}   ,  # K fake number, does not exhist in ODB file
                             'relative_humidity'  : { 'odb_var': 29     , 'cdm_unit': 0    , 'cdm_var': 138}     ,  # Relative humidity (from profile measurement)
                             'gph'                       : { 'odb_var': 1       , 'cdm_unit': 631         , 'cdm_var': 117}    ,  # geopotential height
                             'pressure'               : { 'odb_var': 999    , 'cdm_unit': 32       , 'cdm_var': 57}      , # Pa  (it goes into z_coordinate type, does not exist in odb files )
                              }
    '''

    # data placeholder  
    read_data = []
    proc_dates = []

    obs_id = 0
    report_id = 0 

    # must find lat and lon from station configuration
    stat_inv = pd.read_csv('/srvfs/home/uvoggenberger/CEUAS/CEUAS/public/harvest/data/station_configurations/giub_station_configuration_extended.csv' , sep = '\t')
    if "A" in statid or 'B' in statid:
        statid = statid[:4]
    s = stat_inv.loc[stat_inv['secondary_id'].astype(str) == statid] 
    try:
        lat, lon = s.latitude.values[0] , s.longitude.values[0]
    except:
        lat, lon = -999, -999 

    for i in range(len(df)):

        temp_v = df['temp'].values[i]
        if 'temp' in temp_v:
            continue

        if '999' in temp_v:
            temp_v = np.nan 
        else:
            temp_v = float(temp_v) + 273.15

        date_v = str(df['date'].values[i])
        time_v = str(df['time'].values[i])
        lat_v = lat
        lon_v = lon

        gp_v = df['geopotential'][i]  # converted from geopotential height
        if float(gp_v) <0:
            gp_v = np.nan 

        gph_v = df['gph'][i] # geopotential height
        if float(gph_v) <0:
            gph_v = np.nan 

        rh_v = df['rh'].values[i]
        if '999' in rh_v or '-' in rh_v or float(rh_v) < 0:
            rh_v = np.nan 

        wind_sp_v = df['wspeed'].values[i]
        if '999' in wind_sp_v or '-' in wind_sp_v:
            wind_sp_v = np.nan 

        wind_dir_v = df['wdir'].values[i]
        if '999' in wind_dir_v or '-' in wind_dir_v:
            wind_dir_v = np.nan 

        wind_u_v = df['vwind'].values[i]
        if '999' in wind_u_v or '-' in wind_u_v:
            wind_u_v = np.nan 

        wind_v_v = df['uwind'].values[i]   
        if '999' in wind_v_v or '-' in wind_v_v:
            wind_v_v = np.nan                

        z_coord = df['z_coordinate'].values[i]
        if '999' in z_coord or '-' in z_coord:
            z_coord = np.nan 

        z_type = df['z_coordinate_type'].values[i] # CHECK what is 1 and 2 for this convention 

        # drop lines with all empty variables


        #timestamp = pd.Timestamp(date_v[0:4] + '-'  + date_v[4:6] + '-' + date_v[6:8] + '-' + time_v )
        #if time_v != '0':
        #    time_v = time_v.replace('0','')
        # sometimes issues with hours and minutes 
        try:
            h = time_v.split(':')[0]
            m =  time_v.split(':')[1]
            if h =='00':
                h='0'
            if m == '00':
                m = '0'
            m = int(m.replace('.0',''))
            h = int(h.replace('.0',''))            
            timestamp = pd.Timestamp(year=int(date_v.split('-')[0]), month=int(date_v.split('-')[1]), day=int(date_v.split('-')[2]) , hour=h, minute=m ) 
        except:
            continue
        if timestamp not in proc_dates:
            proc_dates.append(timestamp)
            report_id = report_id + 1


        ### here we need to solve the issue for missing values ofr geopotential or proessure
        ### sometimes the gph is given as a variable but not as a z_coordinate, which might appear to be nan
        ### in this case we copy the variable gph inside the z_coordinate, and change the z_coordinate type
        ### NB first we need to check if at least one variable e.g. temp or wind is not nan, otherwise the data is useless IF only z_coordinate is given
        #missing = ['-999', '-999.0' , -999, -999.0]
        if pd.isnull(z_coord):

            if  pd.isnull(gph_v): #nothing to do, no valid z coordinate at all 
                continue
            else:  # taking gph as valid coordinate 
                z_coord = gph_v
                z_type = 0  # usign gph HEIGHT in m as coordinate 

        if z_type ==1 or z_type =='1' :
            z_coord = float(z_coord) * 100 # converting  from [hPa] to [Pa] NB it is a guess throguh comparison with merged data for same station 

        # should not be necessary anymore 
        #if pd.isnull(gph_v) and z_type in [2 , '2']:
        #    gph_v = float(z_coord) * 9.80665 

        for value,var in zip([ gp_v, temp_v, wind_sp_v, wind_dir_v, wind_u_v, wind_v_v, rh_v ],  [ 'gph', 'temperature', 'wind_speed', 'wind_direction', 'uwind' , 'vwind' , 'relative_humidity'] ):
            if pd.isnull(temp_v) and  pd.isnull(wind_dir_v)  and  pd.isnull(wind_u_v) and  pd.isnull(rh_v)  and  pd.isnull(wind_v_v)  and  pd.isnull(wind_sp_v):  # not a single value available
                continue
            else:
                obs_id = obs_id +1          
                read_data.append( ( 'GIUB'.rjust(10), int(obs_id), report_id,  timestamp, date_v, statid, lat_v, lon_v, z_coord, value, cdmvar_dic[var]['cdm_var'] , int(cdmvar_dic[var]['cdm_unit']), z_type) )

    column_names = [ 'source_id', 'observation_id', 'report_id', 'record_timestamp' , 'iday', 'station_id', 'lat@hdr', 'lon@hdr', 'vertco_reference_1@body', 'obsvalue@body', 'varno@body' ,  'units', 'vertco_type@body'  ]    

    df = pd.DataFrame(data= read_data, columns=column_names )       

    df['observation_id']  = np.chararray.zfill( (df['observation_id'].astype(int)) .astype('S'+str(id_string_length ) ), id_string_length  )  #converting to fixed length bite objects 
    df['report_id']           = np.chararray.zfill( (df['report_id'].astype(int)).astype ('S'+str(id_string_length ) ), id_string_length  )

    df['report_timestamp'] = df['record_timestamp'] 
    #df = df.sort_values(by = ['report_timestamp', 'vertco_reference_1@body' ] ) 

    print('Done reading DF')
    return df , statid 




def read_amma_csv(file=''):
    """ Read AMMA station files, which were previously created via jupyter notebbok,
    by splitting BUFR files. 

        Args:
             file (str): path to the AMMA station file

        Returns:
             Pandas DataFrame with cdm compliant column names
    """    

    # column names standardized:  ['source_id', 'report_id', 'observation_id', 'record_timestamp', 'iday', 'statid@hdr', 'lat@hdr', 'lon@hdr', 
    # 'vertco_reference_1@body', 'obsvalue@body', 'varno@body', 'units', 'number_of_pressure_levels', 'vertco_type@body']

    #col = ['index', 'statid', 'latitude', 'longitude', 'pressure',
    #       'airTemperature', 'nonCoordinateGeopotentialHeight', 'typicalDate',
    #       'typicalTime', 'windSpeed', 'windDirection', 'dewpointTemperature',
    #       'heightOfStation']    

    df = pd.read_csv( file, sep = '\t')
    statid = file.split('/')[-1].split('_')[0]

    # dictionary placeholder 
    all_standard_variables =  ['source_id', 'report_id', 'observation_id', 'record_timestamp', 'iday', 'statid@hdr', 'lat@hdr', 'lon@hdr', 
                               'vertco_reference_1@body', 'obsvalue@body', 'varno@body', 'units', 'number_of_pressure_levels', 'vertco_type@body'] 
    all_data = {}
    for v in all_standard_variables:
        all_data[v] = []

    ''' # just a reminder 
    cdmvar_dic = {'temperature'          : { 'odb_var': 2      , 'cdm_unit': 5        , 'cdm_var': 126}     ,  # K, Air temperature (from profile measurement) 
                             'wind_direction'      : { 'odb_var': 111   , 'cdm_unit': 110    , 'cdm_var': 106}   ,  # degree (angle) direction from which the wind is blowing Lot 1 uses dd  - WMO abbrev.
                             'wind_speed'           : { 'odb_var': 112  , 'cdm_unit': 731     , 'cdm_var': 107 } ,  # m/s , Speed is the magnitude of velocity.
                             'uwind'                    : { 'odb_var': 3       , 'cdm_unit': 731     , 'cdm_var': 139}   ,  # m/s, Eastward wind speed (from profile measurement)
                             'vwind'                    : { 'odb_var': 4       , 'cdm_unit': 731     , 'cdm_var': 140}    ,  # m/s, Northward wind speed (from profile measurement)
                             'dew_point'             : { 'odb_var': 59             , 'cdm_unit': 5  , 'cdm_var': 137}     ,  # K, Dewpoint measurement (from profile measurement) 
                             'dew_point_depression' : { 'odb_var': 299  , 'cdm_unit': 5     , 'cdm_var': 34}   ,  # K fake number, does not exhist in ODB file
                             'relative_humidity'  : { 'odb_var': 29     , 'cdm_unit': 0    , 'cdm_var': 138}     ,  # Relative humidity (from profile measurement)
                             'gph'                       : { 'odb_var': 1       , 'cdm_unit': 631         , 'cdm_var': 117}    ,  # geopotential height
                             'pressure'               : { 'odb_var': 999    , 'cdm_unit': 32       , 'cdm_var': 57}      , # Pa  (it goes into z_coordinate type, does not exist in odb files )
                              }
    '''
    # data placeholder  
    read_data = []
    proc_dates = []

    obs_id = 0
    report_id = 0 

    for i in range(len(df)):
        temp_v = df['airTemperature'].values[i]
        press_v = df['pressure'].values[i]
        date_v = str(df['typicalDate'].values[i])
        time_v = str(df['typicalTime'].values[i])
        lat_v = df['latitude'].values[i]
        lon_v = df['longitude'].values[i]
        gph_v = df['nonCoordinateGeopotentialHeight'][i]
        dp_v = df['dewpointTemperature'].values[i]
        wind_sp_v = df['windSpeed'].values[i]
        wind_dir_v = df['windDirection'].values[i]

        #timestamp = pd.Timestamp(date_v[0:4] + '-'  + date_v[4:6] + '-' + date_v[6:8] + '-' + time_v )
        if time_v != '0':
            time_v = time_v.replace('0','')
        timestamp = pd.Timestamp(year=int(date_v[0:4]), month=int(date_v[4:6]), day=int(date_v[6:8]) , hour=int( time_v ) ) 

        if timestamp not in proc_dates:
            proc_dates.append(timestamp)
            report_id = report_id + 1

        for value,var in zip([ gph_v, temp_v, wind_sp_v, wind_dir_v, dp_v ],  [ 'gph', 'temperature', 'wind_speed', 'wind_direction', 'dew_point'] ):
            obs_id = obs_id +1
            z_type = 1                              

            read_data.append( ( 'AMMA'.rjust(10), int(obs_id), report_id,  timestamp, date_v, statid, lat_v, lon_v, press_v, value, cdmvar_dic[var]['cdm_var'] , int(cdmvar_dic[var]['cdm_unit']), z_type) )


    #column_names = ['source_file', 'product_code', 'report_id', 'observation_id', 'report_timestamp' , 'iday', 'station_id', 'lat@hdr', 'lon@hdr', 'vertco_reference_1@body', 'obsvalue@body', 'varno@body' ,  'units',  'number_of_pressure_levels' ]    

    column_names = [ 'product_code', 'observation_id', 'report_id', 'record_timestamp' , 'iday', 'station_id', 'lat@hdr', 'lon@hdr', 'vertco_reference_1@body', 'obsvalue@body', 'varno@body' ,  'units', 'vertco_type@body'  ]    

    df = pd.DataFrame(data= read_data, columns=column_names )       

    df['observation_id']  = np.chararray.zfill( (df['observation_id'].astype(int)) .astype('S'+str(id_string_length ) ), id_string_length  )  #converting to fixed length bite objects 
    df['report_id']           = np.chararray.zfill( (df['report_id'].astype(int)).astype ('S'+str(id_string_length ) ), id_string_length  )

    df['report_timestamp'] = df['record_timestamp'] 
    #df = df.sort_values(by = ['record_timestamp', 'vertco_reference_1@body' ] ) 

    print('Done reading DF')
    return df , statid 

def read_bufr_cnr_csv(file=''):
    """ Read converted station files from the CNR BUFR data, which were created by splitting BUFR files. 

        Args:
             file (str): path to the BURF_CNR station file

        Returns:
             Pandas DataFrame with cdm compliant column names
    """    

    df = pd.read_csv( file, sep = '\t')
    statid = df.statid.iloc[-1]

    # dictionary placeholder 
    all_standard_variables =  ['source_id', 'report_id', 'observation_id', 'record_timestamp', 'iday', 'statid@hdr', 'lat@hdr', 'lon@hdr', 'radiosondeType',
                               'vertco_reference_1@body', 'obsvalue@body', 'varno@body', 'units', 'number_of_pressure_levels', 'vertco_type@body'] 
    all_data = {}
    for v in all_standard_variables:
        all_data[v] = []

    ''' # just a reminder 
    cdmvar_dic = {'temperature'          : { 'odb_var': 2      , 'cdm_unit': 5        , 'cdm_var': 126}     ,  # K, Air temperature (from profile measurement) 
                             'wind_direction'      : { 'odb_var': 111   , 'cdm_unit': 110    , 'cdm_var': 106}   ,  # degree (angle) direction from which the wind is blowing Lot 1 uses dd  - WMO abbrev.
                             'wind_speed'           : { 'odb_var': 112  , 'cdm_unit': 731     , 'cdm_var': 107 } ,  # m/s , Speed is the magnitude of velocity.
                             'uwind'                    : { 'odb_var': 3       , 'cdm_unit': 731     , 'cdm_var': 139}   ,  # m/s, Eastward wind speed (from profile measurement)
                             'vwind'                    : { 'odb_var': 4       , 'cdm_unit': 731     , 'cdm_var': 140}    ,  # m/s, Northward wind speed (from profile measurement)
                             'dew_point'             : { 'odb_var': 59             , 'cdm_unit': 5  , 'cdm_var': 137}     ,  # K, Dewpoint measurement (from profile measurement) 
                             'dew_point_depression' : { 'odb_var': 299  , 'cdm_unit': 5     , 'cdm_var': 34}   ,  # K fake number, does not exhist in ODB file
                             'relative_humidity'  : { 'odb_var': 29     , 'cdm_unit': 0    , 'cdm_var': 138}     ,  # Relative humidity (from profile measurement)
                             'gph'                       : { 'odb_var': 1       , 'cdm_unit': 631         , 'cdm_var': 117}    ,  # geopotential height
                             'pressure'               : { 'odb_var': 999    , 'cdm_unit': 32       , 'cdm_var': 57}      , # Pa  (it goes into z_coordinate type, does not exist in odb files )
                              }
    '''
    # data placeholder  
    read_data = []
    proc_dates = []

    obs_id = 0
    report_id = 0 

    for i in range(len(df)):  # len(df)
        temp_v = df['airTemperature'].values[i]
        press_v = df['pressure'].values[i]
        date_v = str(df['date'].values[i])
        time_v = str(df['time'].values[i])
        lat_v = df['latitude'].values[i]
        lon_v = df['longitude'].values[i]
        gph_v = df['nonCoordinateGeopotentialHeight'][i]
        dp_v = df['dewpointTemperature'].values[i]
        wind_sp_v = df['windSpeed'].values[i]
        wind_dir_v = df['windDirection'].values[i]
        lat_d = df['latitudeDisplacement'].values[i]
        lon_d = df['longitudeDisplacement'].values[i]
        time_d = df['timePeriod'].values[i]
        sonde_type = df['radiosondeType'].values[i]

        #timestamp = pd.Timestamp(date_v[0:4] + '-'  + date_v[4:6] + '-' + date_v[6:8] + '-' + time_v )
        original_time = time_v
        if "2147483647" in time_v:
            time_v = time_v.replace("2147483647", "00")
        time_v = time_v.zfill(6)
        try:
            timestamp = pd.Timestamp(year=int(date_v[0:4]), month=int(date_v[4:6]), day=int(date_v[6:8]) , hour=int(time_v[0:2]), minute=int(time_v[2:4]), second=int(time_v[4:6]) ) 
        except:
            with open('errors.txt', 'a') as f:
                f.write('-----------------------\n')
                f.write('ERROR: ' + str(original_time) + ' --> ' + str(time_v))
                f.write('\nFile: ' + str(file))
            print('-----------------------')
            print('ERROR: ', original_time, ' --> ', time_v)
            print('File: ', file)

        if timestamp not in proc_dates:
            proc_dates.append(timestamp)
            report_id = report_id + 1

        for value,var in zip([ gph_v, temp_v, wind_sp_v, wind_dir_v, dp_v ],  [ 'gph', 'temperature', 'wind_speed', 'wind_direction', 'dew_point'] ):
            obs_id = obs_id +1
            z_type = 1                              

            read_data.append( ( 'BUFR_CNR'.rjust(10), int(obs_id), report_id,  timestamp, date_v, statid, lat_v, lon_v, press_v, value, cdmvar_dic[var]['cdm_var'] , int(cdmvar_dic[var]['cdm_unit']), z_type, lat_d, lon_d, time_d, sonde_type))


    column_names = [ 'product_code', 'observation_id', 'report_id', 'record_timestamp' , 'iday', 'station_id', 'lat@hdr', 'lon@hdr', 'vertco_reference_1@body', 'obsvalue@body', 'varno@body' ,  'units', 'vertco_type@body', 'observed_latitude_displacement', 'observed_longitude_displacement', 'observed_time_displacement', 'sonde_type@conv']    

    df = pd.DataFrame(data= read_data, columns=column_names )       

    df['observation_id']  = np.chararray.zfill( (df['observation_id'].astype(int)) .astype('S'+str(id_string_length ) ), id_string_length  )  #converting to fixed length bite objects 
    df['report_id']           = np.chararray.zfill( (df['report_id'].astype(int)).astype ('S'+str(id_string_length ) ), id_string_length  )

    df['report_timestamp'] = df['record_timestamp'] 
    #df = df.sort_values(by = ['record_timestamp', 'vertco_reference_1@body' ] ) 

    print('Done reading DF')
    return df , statid 


def read_woudc_csv(file=''):
    """ Read converted station files from the WOUDC data, which were created by splitting BUFR files. 

        Args:
             file (str): path to the WOUDC station file

        Returns:
             Pandas DataFrame with cdm compliant column names
    """    

    df = pd.read_csv( file, sep = ',')
    statid = df.primary_station_id.iloc[-1][2:-1]
    # dictionary placeholder 
    all_standard_variables =  ['source_id', 'report_id', 'observation_id', 'record_timestamp', 'iday', 'statid@hdr', 'lat@hdr', 'lon@hdr', 'radiosondeType',
                               'vertco_reference_1@body', 'obsvalue@body', 'varno@body', 'units', 'number_of_pressure_levels', 'vertco_type@body'] 
    all_data = {}
    for v in all_standard_variables:
        all_data[v] = []

    # data placeholder  
    read_data = []
    proc_dates = []

    obs_id = 0
    report_id = 0 

    dfs = {category: sub_df for category, sub_df in df.groupby('observed_variable')}

    for i in range(len(dfs[126.0])):  # len(df)
        temp_v = dfs[126.0]['observation_value'].values[i]
        press_v = dfs[126.0]['z_coordinate'].values[i]
        date_v = str(dfs[126.0]['report_timestamp'].values[i])[:10]
        time_v = str(dfs[126.0]['report_timestamp'].values[i])[11:]
        lat_v = dfs[126.0]['latitude|header_table'].values[i]
        lon_v = dfs[126.0]['longitude|header_table'].values[i]
        gph_v = dfs[117.0]['observation_value'].values[i]
        rh_v = dfs[138.0]['observation_value'].values[i]
        wind_sp_v = dfs[107.0]['observation_value'].values[i]
        wind_dir_v = dfs[106.0]['observation_value'].values[i]
        sonde_type = dfs[126.0]['sensor_id'].values[i][2:-1]

        #timestamp = pd.Timestamp(date_v[0:4] + '-'  + date_v[4:6] + '-' + date_v[6:8] + '-' + time_v )
        original_time = time_v
        try:
            timestamp = pd.Timestamp(year=int(date_v[0:4]), month=int(date_v[5:7]), day=int(date_v[8:10]) , hour=int(time_v[0:2]), minute=int(time_v[3:5]), second=int(time_v[6:8]) ) 
        except:
            with open('errors.txt', 'a') as f:
                f.write('-----------------------\n')
                f.write('ERROR: ' + str(original_time) + ' --> ' + str(time_v))
                f.write('\nFile: ' + str(file))
            print('-----------------------')
            print('ERROR: ', original_time, ' --> ', time_v)
            print('File: ', file)

        if timestamp not in proc_dates:
            proc_dates.append(timestamp)
            report_id = report_id + 1

        for value,var in zip([ gph_v, temp_v, wind_sp_v, wind_dir_v, rh_v ],  [ 'gph', 'temperature', 'wind_speed', 'wind_direction', 'relative_humidity'] ):
            obs_id = obs_id +1
            z_type = 1                              

            read_data.append( ( 'WOUDC'.rjust(10), int(obs_id), report_id,  timestamp, date_v, statid, lat_v, lon_v, press_v, value, cdmvar_dic[var]['cdm_var'] , int(cdmvar_dic[var]['cdm_unit']), z_type, sonde_type))


    column_names = [ 'product_code', 'observation_id', 'report_id', 'record_timestamp' , 'iday', 'station_id', 'lat@hdr', 'lon@hdr', 'vertco_reference_1@body', 'obsvalue@body', 'varno@body' ,  'units', 'vertco_type@body', 'sonde_type@conv']    

    df = pd.DataFrame(data= read_data, columns=column_names )       

    df['observation_id']  = np.chararray.zfill( (df['observation_id'].astype(int)) .astype('S'+str(id_string_length ) ), id_string_length  )  #converting to fixed length bite objects 
    df['report_id']           = np.chararray.zfill( (df['report_id'].astype(int)).astype ('S'+str(id_string_length ) ), id_string_length  )

    df['report_timestamp'] = df['record_timestamp'] 
    #df = df.sort_values(by = ['record_timestamp', 'vertco_reference_1@body' ] ) 

    print('Done reading DF')
    return df , statid 



def read_yangjiang_csv(file='', metadata=''):
    """ Read the Yangjiang intercomparison """

    # '/scratch/das/federico/databases_service2/YANGJIANG/intercomparison_yangjiang_2010//Graw'

    files = glob.glob(file + '/*')
    files.sort()

    # extracting names of columns 
    df = pd.read_csv(files[0], nrows=1)
    cols = [ c for c in list(df.columns)[0].split('\t') if c ]    


    all_df = []

    ### Extracting the sensor from the file name 
    ### ['3therm', 'CFH', 'Changf', 'Daqiao', 'Graw', 'Huayun', 'IntermSA', 'Jinyang', 'LM', 'ML', 'ML_SW', 'Meis_REF', 'Meisei', 'Modem', 'Vais_ref', 'Vaisala']


    sensors = [ file.split('/')[-1] ]     

    ### List of available variable 
    variables = ['Humidity' , 'Direction', 'Velocity' , 'Temperature']

    yang_sensor_map = { '3therm' : '3therm_y',
                        'CFH': 'CFH_y',
                    'Changf' : 'HCB' ,
                    'Daqiao' : '143',
                    'Graw' : '117',
                    'Huayun' : 'HY1', 
                    'IntermSA' : 'BfA',
                    'Jinyang' : '112',
                    'LM' : '111' , # lock head martin sippican 
                    'ML' : 'SGR' ,  # meteolabor
                    'ML_SW' : 'SRQ' , #meteolabor snow white 
                    'Meisei': '130',
                    'Modem': 'FGD',
                    'Vaisala' : 'VNH',                       

                    'Vais_ref' : 'Vais_ref_y',
                    'Meis_REF' : 'Meis_REF_y',
                    }


    def create_df(files, cols= '', variables='',  variable='Humidity', sensor=''):
        """ Extract data from each single file and combine """

        dic_map = { 'Humidity':'relative_humidity' , 
                    'Temperature' : 'temperature' , 
                              'Direction' : 'wind_direction', 
                              'Velocity': 'wind_speed',
                              'Dew_point': 'dew_point',
                              }

        temp_df = []

        for record,f in enumerate(files):  # keeping humidity data 
            record = int( f.split('.')[1] )

            timestamp_record = metadata.iloc[record-1].Date_Time 
            ts = pd.to_datetime(timestamp_record)
            # entries are spaced by empty characters ' ' with variable length, so cannot really use pandas 
            df = pd.read_csv(f, skiprows=1, names= cols, header=None ,  delim_whitespace=True )        


            df = df.loc[df[variable] != 'NODATA']        

            if variable == 'Humidity':
                df['Humidity'] = df['Humidity'].astype(float)/ 100

            if variable == 'Dew_point':
                df['Dew_point'] = df['Dew_point'].astype(float) + 273.15

            if variable == 'Temperature':
                df['Temperature'] = df['Temperature'].astype(float) + 273.15

            df['observed_variable'] = cdmvar_dic[dic_map[variable]]['cdm_var'] 
            df['unit'] = int(cdmvar_dic[dic_map[variable]]['cdm_unit'])
            df['varno@body'] = cdmvar_dic[dic_map[variable]]['cdm_var']         
            df['obsvalue@body'] = df[variable]

            if 'Lat' in cols:
                df['lat@hdr'] = df['Lat']
                df['lon@hdr'] = df['Long']

                df = df.loc[df['lon@hdr'] != 'NODATA']        
                df = df.loc[df['lat@hdr'] != 'NODATA']        


            df['iday'] =  str(ts.date()).replace('-','').replace('-','')


            if sensor != 'Meis_REF':
                df = df.loc[df['Pressure'] != 'NODATA']        
                df['Pressure'] =  [float(p) * 100 for p in df['Pressure'].values ]  # converting from [hPa] to [Pa]                
                df['vertco_reference_1@body'] = df['Pressure']
                df['vertco_type@body'] = 1
            else:
                df['vertco_reference_1@body'] = df['Height']
                df['vertco_type@body'] = 0

            sd = yang_sensor_map[s]  #mapping the sensor id to Schroeder/WMO entries 

            df['sensor_id']  = np.full(  len(df), sd.rjust(10)).astype('S'+str(id_string_length )  ) 
            df['source_id']  = np.full(  len(df), 'YANGJIANG').astype('S'+str(id_string_length )  )             
            df['report_timestamp'] = ts 


            df = df.loc[df['vertco_reference_1@body'] != 'NODATA']
            df = df.loc[df[variable] != 'NODATA']

            to_drop = [v for v in variables if v != variable ]
            to_drop.extend(cols)

            to_drop = [c for c in  to_drop if c in df.columns ]
            df = df.drop(columns=to_drop )
            if not df.empty:
                temp_df.append(df)

        try:

            df = pd.concat(temp_df)
        except:
            for i in range(len(temp_df)):
                try:

                    a = pd.concat(temp_df[:i])
                    print(i, '  works ')
                except:
                    print(i , ' does not work ')
        return df 

    variables = ['Humidity' , 'Direction', 'Velocity' , 'Temperature']
    if 'Dew_point' in cols:
        variables.append('Dew_point')
    for s in sensors:
        for v in variables:
            if s=='Meis_REF' and v != 'Temperature':
                continue
            df_var =  create_df(files, variables= variables, cols=cols, variable=v, sensor=s)
            all_df.append(df_var)

    df_res = pd.concat(all_df)

    #read_data.append( ( product.rjust(20), int(obs_id), report_id,  timestamp, int(date_v), statid, lat_v, lon_v, z_coordinate_v, value, cdmvar_dic[var]['cdm_var'] , int(cdmvar_dic[var]['cdm_unit']), z_type) )

    ### Adding coordinates from station configuration 
    if 'Lat' not in cols:
        df_res['lat@hdr'] = cdm_tab['station_configuration'].latitude.values[0]
        df_res['lon@hdr'] = cdm_tab['station_configuration'].longitude.values[0] 

    df_res['product_code'] = np.full(  len(df_res), 'YANGJIANG').astype('S'+str(id_string_length )  ) 

    ### sorting by date-time 
    df_res = df_res.sort_values(by=['report_timestamp', 'vertco_reference_1@body'])

    ### creating report_id i.e. same for
    ts,counts = np.unique(df_res.report_timestamp , return_counts=True )

    report_id = []
    for index,c in enumerate(counts):
        report_id.extend( [index]*c)  # report id 

    df_res['report_id'] = report_id
    df_res['report_id']           = np.chararray.zfill( (df_res['report_id'].astype(int)).astype ('S'+str(id_string_length ) ), id_string_length  )    
    df_res['observation_id'] = list(range(len(df_res)))
    df_res['observation_id']  = np.chararray.zfill( (df_res['observation_id'].astype(int)) .astype('S'+str(id_string_length ) ), id_string_length  )  #converting to fixed length bite objects 

    df_res['record_timestamp'] = df_res['report_timestamp'] 

    statid = 'YANGJIANG'

    print('Done reading DF')
    return df_res , statid 


def read_mauritius_csv(file=''):
    """ Read the Mauritius intercomparison data for Vaisala and Maisei sondes 
        Args:
             file (str): path to the intercomparison file
        Returns:
             Pandas DataFrame with cdm compliant column names
    """    
    # column names standardized:  ['source_id', 'report_id', 'observation_id', 'record_timestamp', 'iday', 'statid@hdr', 'lat@hdr', 'lon@hdr', 
    # 'vertco_reference_1@body', 'obsvalue@body', 'varno@body', 'units', 'number_of_pressure_levels', 'vertco_type@body']

    #col =Index(['date_time', 'gph', 'temp', 'dew', 'windsp', 'winddir', 'lat', 'lon',
    #   'pressure', 'month', 'day', 'year', 'hour', 'elev'], 

    df = pd.read_csv( file, sep = ',' ).astype(str)

    statid = file.split('/')[-1].split('_')[0] + '_mauritius'
    """
    MEISEI
    Index(['Unnamed: 0', 'Temperature 0', 'Observation Time', 'Wind Direction',
       'Wind Speed', 'Humidity 0', 'Barometric Pressure 0',
       'Positioning latitude', 'Positioning Longitude'],
      dtype='object')

    VAISALA
    Index(['Unnamed: 0', 'number', 'pressure', 'temperature', 'relative_humidity',
       'geopotential', 'wdir', 'wspeed', 'height', 'date_time'],
      dtype='object')

    """


    sensor_map = { 
        'Meisei' : 'J0A',
                   'Vaisala' : 'VN2' ,
                   'Vaisala-GPS' : 'VN2' ,
    }

    if 'meisei' in file:
        lat, lon = 'Positioning latitude' , 'Positioning Longitude' 
        relhum, temp = 'Humidity 0', 'Temperature 0' 
        wspeed, wdir = 'Wind Speed' , 'Wind Direction' 
        datetime ='Observation Time'
        product = 'MEISEI'
        z_coordinate =  'Barometric Pressure 0' 
        date = 'Date'
        z_type=1
        sensor = sensor_map['Meisei']

    elif 'vaisala' in file:
        lat, lon = 'number' , 'number' 
        relhum, temp = 'relative_humidity', 'temperature' 
        wspeed, wdir = 'wspeed' , 'wdir' 
        press, datetime = 'pressure' , 'date_time'
        product = 'VAISALA'
        z_coordinate =  'pressure' 
        z_type=1
        sensor = sensor_map['Vaisala']

    # dictionary placeholder 
    all_standard_variables =  ['source_id', 'report_id', 'observation_id', 'report_timestamp', 'iday', 'statid@hdr', 'lat@hdr', 'lon@hdr', 
                               'vertco_reference_1@body', 'obsvalue@body', 'varno@body', 'units', 'number_of_pressure_levels', 'vertco_type@body'] 

    all_data = {}
    for v in all_standard_variables:
        all_data[v] = []
    # data placeholder  
    read_data = []
    proc_dates = []

    obs_id = 0
    report_id = 0 

    for i in tqdm(range(len(df)),  miniters=int(len(df)/10000) ):

        date_time_v = df[datetime].values[i]

        # z_coordinate
        z_coordinate_v = df[z_coordinate].values[i]
        try:    
            z_coordinate_v = float(z_coordinate_v)  *100 # P in hPa
        except:
            continue

        # dates, timestamps 
        if 'meisei' in file:    
            date_v = df[date].values[i]
            time = date_time_v.split(':')

            year =  int(date_v[0:4]) 
            month = int(date_v[4:6])
            day = int(date_v[6:8]) 

            hour = int(time[0])
            minute = int(time[1]) 
            sec = int(time[2])

            timestamp = pd.Timestamp(year=year, month=month, day=day , hour=hour, minute=minute, second=sec  ) 

        elif 'vaisala' in file:
            date_v = df[datetime].values[i]

            date = date_v.split('-')            

            day = int(date[0])
            month = int(date[1])
            year = int(date[2][0:4])

            hour = int(date[2][5:7])
            minute = int(date[2][8:10])
            sec = int(date[2][11:13])
            timestamp = pd.Timestamp(year=year, month= month, day= day , hour= hour, minute= minute , second= sec  ) 
            date_v = str(timestamp.date()).replace('-','')

        # lat, lon 
        try:
            lat_v = float(df[lat].values[i])
            lon_v = float(df[lon].values[i])
        except:
            lat_v = -20.2972  # using values from MAISEI +> NEED To FIX THIS for VAISALA 
            lon_v = 57.49692
        if lon_v > 180:
            lon_v = -180.0 + (lon_v - 180)

        # temp, wind speed, wind dir 
        temp_v = df[temp].values[i] # 

        try:
            temp_v = float(temp_v) + 273.15 # temp is given in degree Celsius -> convert to Kelvin
        except:
            temp_v = np.nan

        relhum_v = df[relhum].values[i] # 
        try:
            relhum_v = float(relhum_v)/100.0 # normalize to [0,1] range 
        except:
            relhum_v = np.nan

        wind_sp_v = df[wspeed].values[i]  # 
        try:
            wind_sp_v = float(wind_sp_v)
        except:
            wind_sp_v = np.nan

        wind_dir_v = df[wdir].values[i]  
        try:
            wind_dir_v = float(wind_dir_v)
        except:
            wind_dir_v = np.nan

        # to remove duplicates
        if timestamp not in proc_dates:
            proc_dates.append(timestamp)
            report_id = report_id + 1

        for value,var in zip([ temp_v, relhum_v, wind_sp_v, wind_dir_v],  [ 'temperature', 'relative_humidity', 'wind_speed', 'wind_direction' ] ):
            obs_id = obs_id +1
            read_data.append( ( product.rjust(20), int(obs_id), report_id,  timestamp, int(date_v), statid, lat_v, lon_v, z_coordinate_v, value, cdmvar_dic[var]['cdm_var'] , int(cdmvar_dic[var]['cdm_unit']), z_type) )


    column_names = [ 'product_code', 'observation_id', 'report_id', 'report_timestamp' , 'iday', 'station_id', 'lat@hdr', 'lon@hdr', 'vertco_reference_1@body', 'obsvalue@body', 'varno@body' ,  'units', 'vertco_type@body'  ]    


    df = pd.DataFrame(data= read_data, columns=column_names )       

    df['observation_id']  = np.chararray.zfill( (df['observation_id'].astype(int)) .astype('S'+str(id_string_length ) ), id_string_length  )  #converting to fixed length bite objects 
    df['report_id']           = np.chararray.zfill( (df['report_id'].astype(int)).astype ('S'+str(id_string_length ) ), id_string_length  )

    df['record_timestamp'] = df['report_timestamp'] 
    df['vertco_type@body'] = df['vertco_type@body'].astype(float)

    #df['vertco_reference_1@body'] = df['vertco_reference_1@body'].astype(float)
    #df = df.sort_values(by = ['report_timestamp', 'vertco_reference_1@body' ] ) 

    print('Done reading DF')


    # need to extract record timestamps for MEISEI case
    print(' === Creating record_timestamps === ')
    if 'meisei' in file:
        record_ts = []
        report_id = []

        for t in  tqdm(range(len(df)),  miniters=int(len(df)/100) ): 
            ts = df.report_timestamp.values[t] 

            previous = df.report_timestamp.values[t-1]

            if t ==0:
                record_ts.append(ts)
                report_id.append(0)
                continue
            timedelta = pd.Timedelta(4, unit='h')

            if (ts - previous) <  timedelta:
                record_ts.append(record_ts[-1])
                report_id.append(report_id[-1])

            else:
                record_ts.append(ts)
                report_id.append(report_id[-1]+1)

        df['report_id'] = np.chararray.zfill( ( np.array(report_id).astype(int)).astype ('S'+str(id_string_length ) ), id_string_length  )
        df['record_timestamp'] = record_ts


    # Adding sensor id 
    df['sensor_id'] = np.full(  len(df), sensor.rjust(10)).astype('S'+str(id_string_length )  ) 

    return df , statid 


def read_mauritius_csv_digitized(direc=''):
    """ Read the Mauritius intercomparison data for Vaisala and Maisei sondes 
        Args:
             file (str): path to the intercomparison file
             sensors list: list of sensor to process. 
        Returns:
             Pandas DataFrame with cdm compliant column names
    """    
    sensors = list( np.unique( [ s.split('_')[1].replace('.csv','') for s in os.listdir(direc +'/temp') ] ) )     
    def make_dt(dt):
        """ Internal utility to create pd datetime objects """

        date = dt.split(' ')[0]
        time = dt.split(' ')[1]

        day = int(date.split('-')[0])
        month =  int(date.split('-')[1])
        year =  int(date.split('-')[2])

        hour = int(time.split(':')[0])
        minute =  int(time.split(':')[1])
        sec =  int(time.split(':')[2].split('.')[0])

        timestamp = pd.Timestamp(year=year, month= month, day= day , hour= hour, minute= minute , second= sec  ) 
        date_v = str(timestamp.date()).replace('-','')

        return timestamp, date_v

    all_df = []

    sensor_map = { 
        'Graw' : 'DGL',
                   'Graw-GPS' : 'DGL',

                   'Meisei' : 'J0A',

                   'Vaisala' : 'VN2' ,
                   'Vaisala-GPS' : 'VN2' ,

                   'Modem' : 'FG2',

                   'Sip' : '110' ,  # Sippican LMS5
                   'MKII' : 'ZSm' ,  # Sippican multithermistor, see official report

                   'SRS' : 'SRd' , #Meteolabor
    }

    for s in sensors: 
        print('Harvesting Sensor::: ' , s)
        files = [direc+'/hum/'+f for f in os.listdir(direc+'/hum') if s == f.split('/')[-1].split('_')[1].replace('.csv','')]

        files.extend( [direc+'/temp/' + f for f in os.listdir(direc+'/temp') if s==f.split('/')[-1].split('_')[1].replace('.csv','') ])

        for f in tqdm(files, miniters=10):

            df = pd.read_csv( f, sep = ',' ).astype(str)

            if 'hum' in f:
                df['hum'] = df['hum'].astype(float)/ 100
                df['observed_variable'] = cdmvar_dic['relative_humidity']['cdm_var'] 

                df['unit'] = int(cdmvar_dic['relative_humidity']['cdm_unit'])
                df['varno@body'] = cdmvar_dic['relative_humidity']['cdm_var']              
                df['obsvalue@body'] = df['hum']

            elif 'temp' in f:
                df['temp'] = df['temp'].astype(float)

                df['unit'] = int(cdmvar_dic['temperature']['cdm_unit'])
                df['varno@body'] = cdmvar_dic['temperature']['cdm_var']  # HERE 
                df['obsvalue@body'] = df['temp']

            # creating timestamps 
            dt = df['datetime'][0]
            timestamp,date = make_dt(dt)

            df['iday'] = date 
            df['observation_id'] = list(range(len(df)))
            df['vertco_reference_1@body'] = df['press']
            df['vertco_type@body'] = 1
            df['report_timestamp'] = timestamp

            #df['sensor_id'] = np.bytes_(s)

            ss = sensor_map[s]
            df['sensor_id']  = np.full(  len(df), ss.rjust(10)).astype('S'+str(id_string_length )  ) 

            df['source_id']  = np.full(  len(df), 'MAURITIUS_DIGITIZED').astype('S'+str(id_string_length )  ) 

            #converting to fixed length bite objects 

            all_df.append(df)


    df_res = pd.concat(all_df)
    print(np.unique(df_res.sensor_id))
    cols = [c for c in df.columns if c not in ['datetime', 'press', 'hum', 'temp']]
    df_res = df_res[cols]

    ts = np.unique( df_res.iday)
    dt_map = {}
    for i,t in enumerate(ts):
        dt_map[int(t)] = i

    df_res['iday'] = df_res['iday'].astype(int)
    df_res['report_id'] = df_res['iday'].astype(int)

    lat_v = -20.2972  # using values from MAISEI  
    lon_v = 57.49692

    df_res['lat@hdr'] = lat_v
    df_res['lon@hdr'] = lon_v

    df_res = df_res.replace({"report_id": dt_map}) 

    #rep_ids = [dt_map[i] for i in df_res.report_id ]    

    df_res['product_code'] = np.full(  len(df_res), 'MAURITIUS_DIGITIZED').astype('S'+str(id_string_length )  ) 

    df_res['observation_id']  = np.chararray.zfill( (df_res['observation_id'].astype(int)) .astype('S'+str(id_string_length ) ), id_string_length  )  #converting to fixed length bite objects 
    df_res['report_id']           = np.chararray.zfill( (df_res['report_id'].astype(int)).astype ('S'+str(id_string_length ) ), id_string_length  )

    df_res['record_timestamp'] = df_res['report_timestamp'] 

    statid = 'MAURITIUS_DIGITIZED'
    #print(len(df_res))

    return df_res , statid 




def read_hara_csv(file=''):
    """ Read HARA station files, which were previously created via jupyter notebook 
        Args:
             file (str): path to the HARA station file
        Returns:
             Pandas DataFrame with cdm compliant column names

        See https://nsidc.org/sites/default/files/nsidc-0008-v001-userguide.pdf


        Geopotential conversion: see  https://confluence.ecmwf.int/pages/viewpage.action?pageId=111155328
    """    
    # column names standardized:  ['source_id', 'report_id', 'observation_id', 'record_timestamp', 'iday', 'statid@hdr', 'lat@hdr', 'lon@hdr', 
    # 'vertco_reference_1@body', 'obsvalue@body', 'varno@body', 'units', 'number_of_pressure_levels', 'vertco_type@body']

    #col =Index(['date_time', 'gph', 'temp', 'dew', 'windsp', 'winddir', 'lat', 'lon',
    #   'pressure', 'month', 'day', 'year', 'hour', 'elev'], 

    df = pd.read_csv( file, sep = '\t').astype(str)
    df = df.replace( { '999.0': np.nan , '999' : np.nan , '9999.0': np.nan , '9999' : np.nan , '99999': np.nan , '99999.0': np.nan , '99999 ': np.nan , '-999.0': np.nan } ) 
    
    statid = file.split('/')[-1].split('_')[0]

    # check that this eliminates all problems
    
    for v in ['date_time', 'gph', 'temp', 'dew', 'windsp', 'winddir', 'lat', 'lon',
       'pressure', 'month', 'day', 'year', 'hour', 'elev']:
        
        values = [ val for val in df[v] if not pd.isna(val) and '999' in val  ]
        # print(v , values)
    

    # data placeholder  
    read_data = []
    proc_dates = []

    obs_id = 0
    report_id = 0 

    for i in range(len(df)): 
        dt = df['date_time'].values[i]
        h = df['hour'][i] 

        press_v = df['pressure'].values[i]

        try:    
                press_v = float(press_v) * 10 # Pressure is tenth of millibar, 1 millibar = 100 Pa = 1 hPa hence I multiply by 10
        except:
                continue

        lat_v = float(df['lat'].values[i])

        lon_v = float(df['lon'].values[i])
        if lon_v > 180:
            lon_v = -180.0 + (lon_v - 180)

        temp_v = df['temp'].values[i] # temp is given in tenth of degree Celsius 


        temp_v = float(temp_v) / 10 + 273.15 # temp is given in tenth of degree Celsius 

        gph_v = df['gph'][i]      

        gph_v = float(gph_v)  #* 9.80665 # geopotential height in meters 

        if not np.isnan(press_v):
            z_type = 1
            z_coordinate = press_v
        elif not np.isnan(gph_v) and np.isnan(press_v):
            z_type = 2
            z_coordinate = gph_v            
        elif np.isnan(gph_v) and np.isnan(press_v):  # no valid z cooridnate in this case
            continue


        dp_v = df['dew'].values[i]      
        dp_v = float(dp_v) / 10 + 273.15

        wind_sp_v = df['windsp'].values[i]  # meter per second 

        wind_dir_v = df['winddir'].values[i]  # Wind direction, 0 to 360 degrees, measured clockwise from north (e.g., 90 degrees is east).

        #timestamp = pd.Timestamp(date_v[0:4] + '-'  + date_v[4:6] + '-' + date_v[6:8] + '-' + time_v )
        date_v = dt.replace('-','')
        timestamp = pd.Timestamp(year=int(date_v[0:4]), month=int(date_v[4:6]), day=int(date_v[6:8]) , hour=int( h )  ) 

        # to remove duplicates
        if timestamp not in proc_dates:
            proc_dates.append(timestamp)
            report_id = report_id + 1

        for value,var in zip([ gph_v, temp_v, wind_sp_v, wind_dir_v, dp_v],  [ 'gph', 'temperature', 'wind_speed', 'wind_direction' , 'dew_point'] ):
            obs_id = obs_id +1
            read_data.append( ( 'HARA'.rjust(10), int(obs_id), report_id,  timestamp, int(date_v), int(statid), lat_v, lon_v, z_coordinate, value, cdmvar_dic[var]['cdm_var'] , int(cdmvar_dic[var]['cdm_unit']), z_type) )

    column_names = [ 'product_code', 'observation_id', 'report_id', 'report_timestamp' , 'iday', 'station_id', 'lat@hdr', 'lon@hdr', 'vertco_reference_1@body', 'obsvalue@body', 'varno@body' ,  'units', 'vertco_type@body'  ]    

    df = pd.DataFrame(data= read_data, columns=column_names )       

    df['observation_id']  = np.chararray.zfill( (df['observation_id'].astype(int)) .astype('S'+str(id_string_length ) ), id_string_length  )  #converting to fixed length bite objects 
    df['report_id']           = np.chararray.zfill( (df['report_id'].astype(int)).astype ('S'+str(id_string_length ) ), id_string_length  )

    df['record_timestamp'] = df['report_timestamp'] 
    df['vertco_type@body'] = df['vertco_type@body'].astype(float)

    #df['vertco_reference_1@body'] = df['vertco_reference_1@body'].astype(float)
    #df = df.sort_values(by = ['report_timestamp', 'vertco_reference_1@body' ] ) 

    print('Done reading DF')
    return df , statid 

def read_shipsound_csv(file=''):
    """ Read SHIPSOUND station files, which were previously created via jupyter notebook 
        Args:
             file (str): path to the SHIPSOUND station file
        Returns:
             Pandas DataFrame with cdm compliant column names

        See https://nsidc.org/sites/default/files/nsidc0054_pdf.pdf
        file = '/scratch/das/federico/databases_service2/SHIPSOUND-NSIDC-0054/shipsound7696.out

        Geopotential conversion: see  https://confluence.ecmwf.int/pages/viewpage.action?pageId=111155328
    """    

    df = pd.read_csv( file, sep = '\t').astype(str)

    """
    'date_time', 'height', 'press', 'temp', 'dpd', 'wspeed',
       'wdir', 'latitude', 'longitude', 'elevation', 'ship_type'
    """
    # dictionary placeholder 
    all_standard_variables =  ['source_id', 'report_id', 'observation_id', 'record_timestamp', 'iday', 'statid@hdr', 'lat@hdr', 'lon@hdr', 
                               'vertco_reference_1@body', 'obsvalue@body', 'varno@body', 'units', 'number_of_pressure_levels', 'vertco_type@body'] 
    all_data = {}
    for v in all_standard_variables:
        all_data[v] = []
    # data placeholder  
    read_data = []
    proc_dates = []

    obs_id = 0
    report_id = 0 

    for i in range(len(df)):
        timestamp = df['date_time'].values[i] # already as a timestamp 
        date_v = pd.Timestamp(timestamp)

        press_v = df['press'].values[i]
        try:    
            press_v = float(press_v) * 100  # The information about pressure is missing since the doc is cut. But other pressure are given in millibar , 1 millibar = 100 Pa = 1 hPa so I assume it is the same. hence I convert to [Pa] 
        except:
            press_v = np.nan 

        lat_v = df['latitude'].values[i]
        lon_v = float(df['longitude'].values[i])

        if lon_v > 180:
            lon_v = -180.0 + (lon_v - 180)

        temp_v = df['temp'].values[i] # temp is given in Celsius 
        if '9999' in temp_v:
            temp_v = np.nan
        else:
            try:
                temp_v = float(temp_v) +  273.15 # temp is given in Celsius 
            except:
                print('wrong TEMP at index ' , i )
                temp_v = float(np.nan)
        gph_v = df['height'][i]      # given in meters 
        if '9999' in gph_v or '      ' in gph_v:
            gph_v = np.nan

        else:
            try:
                gph_v = float(gph_v)  * 9.80665 # geopotential height in meters 
            except:
                print('wrong GPH at index ' , i)
                gph_v = float(np.nan)
                pass 

        if not np.isnan(press_v):
            z_type = 1
            z_coordinate = press_v
        elif not np.isnan(gph_v) and np.isnan(press_v):
            z_type = 2
            z_coordinate = gph_v            
        elif np.isnan(gph_v) and np.isnan(press_v):  # no valid z coordinate in this case
            continue

        dp_v = df['dpd'].values[i]      
        if '9999' in dp_v:
            dp_v=np.nan
       
        else:
            try:
                dp_v = float(dp_v) * 10 
            except:
                print('wrong DP at index ' , i)
                dp_v = float(np.nan)
                pass 

        wind_sp_v = df['wspeed'].values[i]  # meter per second 
        if '9999' in wind_sp_v:
            wind_sp_v = np.nan

        wind_dir_v = df['wdir'].values[i]  # Wind direction, 0 to 360 degrees
        if '9999' in wind_dir_v:
            wind_dir_v = np.nan        

        # to remove duplicates
        if timestamp not in proc_dates:
            proc_dates.append(timestamp)
            report_id = report_id + 1

        for value,var in zip([ gph_v, temp_v, wind_sp_v, wind_dir_v, dp_v],  [ 'gph', 'temperature', 'wind_speed', 'wind_direction' , 'dew_point_depression'] ):
            obs_id = obs_id +1
            read_data.append( ( 'SHIPSOUND'.rjust(10), int(obs_id), report_id,  date_v, str(date_v.date()), 99999, lat_v, lon_v, z_coordinate, value, cdmvar_dic[var]['cdm_var'] , int(cdmvar_dic[var]['cdm_unit']), z_type) )

    column_names = [ 'product_code', 'observation_id', 'report_id', 'record_timestamp' , 'iday', 'station_id', 'lat@hdr', 'lon@hdr', 'vertco_reference_1@body', 'obsvalue@body', 'varno@body' ,  'units', 'vertco_type@body'  ]    

    df = pd.DataFrame(data= read_data, columns=column_names )       

    df['obsvalue@body'] = df['obsvalue@body'].astype(float)
    
    df['observation_id']  = np.chararray.zfill( (df['observation_id'].astype(int)) .astype('S'+str(id_string_length ) ), id_string_length  )  #converting to fixed length bite objects 
    df['report_id']           = np.chararray.zfill( (df['report_id'].astype(int)).astype ('S'+str(id_string_length ) ), id_string_length  )

    df['report_timestamp'] = df['record_timestamp']

    print('Done reading DF')
    return df , '20999' 



def read_npsound_csv(file=''):
    """ Read NPSOUND station files, which were previously created via jupyter notebook 
        Args:
             file (str): path to the NPSOUND station file
        Returns:
             Pandas DataFrame with cdm compliant column names

        See https://nsidc.org/sites/default/files/nsidc0054_pdf.pdf
        Geopotential conversion: see  https://confluence.ecmwf.int/pages/viewpage.action?pageId=111155328
    """    

    df = pd.read_csv( file, sep = '\t').astype(str)
    statid = file.split('/')[-1].split('.dat')[0].replace('_','')

    # dictionary placeholder 
    """ Not used ? 
    all_standard_variables =  ['source_id', 'report_id', 'observation_id', 'record_timestamp', 'iday', 'statid@hdr', 'lat@hdr', 'lon@hdr', 
    'vertco_reference_1@body', 'obsvalue@body', 'varno@body', 'units', 'number_of_pressure_levels', 'vertco_type@body'] 
    all_data = {}
    for v in all_standard_variables:
        all_data[v] = []
    """


    # data placeholder  
    read_data = []
    proc_dates = []

    obs_id = 0
    report_id = 0 

    for i in range(len(df)):

        date_time = df['date_time'].values[i]
        date_time = pd.Timestamp(date_time)

        release_time = df['release_time'].values[i]
        release_time = pd.Timestamp(release_time)

        press_v = df['pressure'].values[i]
        if  '99999' in press_v:
            press_v = np.nan 
        else:
            try:    
                press_v = float(press_v) * 100  # Original Pressure is in [hPa]
            except:
                continue

        lat_v = df['latitude'].values[i]

        lon_v = float(df['longitude'].values[i])
        if lon_v > 180:
            lon_v = -180.0 + (lon_v - 180)

        temp_v = df['temp'].values[i] # AIR TEMPERATURE in degrees Celsius
        if '9999' in temp_v:
            temp_v = np.nan
        else:
            temp_v = float(temp_v) + 273.15  # converting to Kelvin  

        gph_v = df['height'][i]      
        if '9999' in gph_v:
            gph_v = np.nan
        else:
            gph_v = float(gph_v)  * 9.80665 # height in meters to geopotential height 

        if not np.isnan(press_v):
            z_type = 1
            z_coordinate = press_v
        elif not np.isnan(gph_v) and np.isnan(press_v):
            z_type = 2
            z_coordinate = gph_v            
        elif np.isnan(gph_v) and np.isnan(press_v):  # no valid z coordinate in this case
            continue

        dp_v = df['dpdp'].values[i]    #  DEW POINT DEPRESSION at the current level in degrees Celsius
        if '9999' in dp_v:
            dp_v=np.nan
        else:
            try:
                dp_v = float(dp_v) + 273.15 # converting to Kelvin 
            except:
                dp_v = np.nan 

        wind_sp_v = df['wspeed'].values[i]  # meter per second 
        if '9999' in wind_sp_v:
            wind_sp_v = np.nan

        wind_dir_v = df['wdir'].values[i]  # Wind direction, 0 to 360 degrees, measured clockwise from north (e.g., 90 degrees is east).
        if '9999' in wind_dir_v:
            wind_dir_v = np.nan        

        rh_v = df['rh'].values[i]  # Wind direction, 0 to 360 degrees, measured clockwise from north (e.g., 90 degrees is east).
        if '-99.99' in rh_v:
            rh_v = np.nan     
        else:
            rh_v = float(rh_v)/100.0 # RELATIVE HUMIDITY ranging from 0 to 100 percent, reset to range [0-1]

        # to remove duplicates
        if date_time not in proc_dates:
            proc_dates.append(date_time)
            report_id = report_id + 1

        for value,var in zip([ gph_v, temp_v, wind_sp_v, wind_dir_v, dp_v, rh_v],  [ 'gph', 'temperature', 'wind_speed', 'wind_direction' , 'dew_point_depression', 'relative_humidity'] ):
            obs_id = obs_id +1
            read_data.append( ( 'NPSOUND'.rjust(10), int(obs_id), report_id,   release_time, int(str(release_time.date()).replace('-','')), statid, lat_v, lon_v, 
                                z_coordinate, value, cdmvar_dic[var]['cdm_var'] , int(cdmvar_dic[var]['cdm_unit']), z_type, date_time , 1 ) )

    column_names = [ 'product_code', 'observation_id', 'report_id', 'report_timestamp' , 'iday', 'station_id', 'lat@hdr', 'lon@hdr', 
                     'vertco_reference_1@body', 'obsvalue@body', 'varno@body' ,  'units', 'vertco_type@body' , 'record_timestamp' , 'report_meaning_of_timestamp'   ]    

    df = pd.DataFrame(data= read_data, columns=column_names )       

    df['observation_id']  = np.chararray.zfill( (df['observation_id'].astype(int)) .astype('S'+str(id_string_length ) ), id_string_length  )  #converting to fixed length bite objects 
    df['report_id']           = np.chararray.zfill( (df['report_id'].astype(int)).astype ('S'+str(id_string_length ) ), id_string_length  )

    print('Done reading DF')
    return df , [statid] 



def igra2_ascii_to_dataframe(file=''):
    """ Read an igra2 stationfile in ASCII format and convert to a Pandas DataFrame. 
        Adapted from https://github.com/MBlaschek/CEUAS/tree/master/CEUAS/data/igra/read.py 
        Variables used inside the DataFrame are already CDM compliant
        df.types

        Args:
             file (str): path to the igra2 station file

        Returns:
             Pandas DataFrame with cdm compliant column names
    """
    if debug:
        print("Running igra2_ascii_to_dataframe for: ", file)    

    data = check_read_file(file=file, read=True)
    #source_file = [l for l in file.split('/') if '.txt' in l][0]
    read_data = [] #  Lists containing the raw data from the ascii file, and the observation dates
    """ Data to be extracted and stored from the igra2 station files 
        Some info is contained in the header of each ascent, some in the following data """

    """ Initialize the variables that can be read from the igra2 files """
    ident,year,month,day,hour,reltime,p_src,np_src,lat, lon = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan 
    lvltyp1,lvltyp2,etime,press,pflag,gph,zflag,temp,tflag,rh,dpdep,wdir,wspd = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan # initialize to zeros
    stations_id = []
    idate = np.nan
    count = 0
    head_count = 0

    obs_id = 0

    def make_release_time(date_time,  hour, release):
        """ build a sonde release time 
              ex 2019 02 20 00 2349 
              ex 2019 01 10 00 0011 
              They round the release time to the closest hour. 
              It can be the same day or the following !!!
              date_time = date_time python object, 
              date, time, release = original strings 
              """
        release_h  = int(release[:2])
        release_m = int(release[2:4])

        if release_h == 99:
            release_date_time = date_time 
            timestamp_flag = int_void #largest integer number int 64   

        else:
            timestamp_flag = 1 

            if release_m == 99:
                release_m = 0
            release_date_time = date_time.replace(hour= release_h, minute= release_m) 

            """ Here, I have to subtract one day to the release time stamp if the hour of the time stamp is in the night,
                  but the nominal time is reported at midnight hence in the following day. For example  2019 02 20 00 2349 from file VMM00048820 """
            if hour == '00':
                if release_h > 20:
                    release_date_time = release_date_time - timedelta(days=1)
                else:
                    pass

        return release_date_time, timestamp_flag 

    #data = data[3706329:]  # TO DO HERE TO DO CHANGE!!!
    #for i, line in enumerate(data):
    #      print(i , '  ' ,  line[0] )

    for i, line in enumerate(data):
        if line[0] == '#':
            head_count = head_count +1 
            # Info from the Header line of each ascent                                                                                                                                                                                                                   
            ident     = line[1:12]               # station identifier
            #ident     = ident[6:12]
            if ident not in stations_id:
                stations_id.append(ident)

            year      = line[13:17]               # year, months, day, hour of the observation
            month   = line[18:20]
            day       = line[21:23]
            hour      = line[24:26]               
            reltime  = line[27:31]            # release time of the sounding.
            numlev  = int(line[32:36])        # number of levels in the sounding == number of data recorded in the ascent
            p_src     = line[37:45]              # data source code for the pressure levels 
            np_src   = line[46:54]             # data source code for non-pressure levels
            lat         = int(line[55:62]) / 10000.  # latitude and longitude
            lon        = int(line[63:71]) / 10000.
            #observation_id = i
            if int(hour) == 99:
                time = reltime + '00'
            else:
                time = hour + '0000'

            if '99' in time:
                time = time.replace('99', '00')

            idate = datetime.strptime(year + month + day + time, '%Y%m%d%H%M%S') # constructed according to CDM

            ### making release time and its flag
            release_time , report_timeflag = make_release_time(idate, hour, reltime) # making the release time 

            ### report_meaning_of_timestamp	int	meaning_of_time_stamp:meaning	Report time - beginning, middle or end of reporting period
            ### 1	beginning	Date / time specified indicates the start of the period over which the observation was made.
            ### end	Date / time specified indicates the end of the period over which the observation was made.
            ### middle	Date / time specified indicates the middle of the period over which the observation was made.

            iday =  int(year + month + day)
            count = count + 1

        else:
            # Data of each ascent
            lvltyp1 = int(line[0])            # 1-  1   integer major level type indicator
            lvltyp2 = int(line[1])            # 2-  2   integer minor level type indicator
            etime   = int(line[3:8])          # 4-  8   integer elapsed time since launch
            press   = int(line[9:15])         # 10- 15  integer reported pressure

            if press == -9999:
                press = np.nan

            pflag   = line[15]                # 16- 16  character pressure processing flag

            gph     = int(line[16:21])        # 17- 21  integer geopotential height  [m]

            if gph == -9999 or gph == -8888:   # reading the values andh check if they are missing or removed as -9999 or -8888 before dividing by 10 as the instructions say 
                gph = np.nan # 23- 27  integer temperature, [Celsius to Kelvin ] 
            else:
                gph = gph * 9.80665 

            zflag   = line[21]                # 22- 22  character gph processing flag, 

            temp    = int(line[22:27])              
            if temp != -9999 and temp != -8888:   # reading the values andh check if they are missing or removed as -9999 or -8888 before dividing by 10 as the instructions say 
                temp = temp / 10.   + 273.15 # 23- 27  integer temperature, [Celsius to Kelvin ]    
            else:
                temp = np.nan 

            tflag   = line[27]    # 28- 28  character temperature processing flag

            rh      = int(line[28:33])  # 30- 34  integer relative humidity [%]           
            if rh != -8888 and rh != -9999:
                rh = rh / 1000.  # converting from percentage to absolute ratio 
            else:
                rh = np.nan

            dpdp    = int(line[34:39]) 
            if dpdp != -9999 and dpdp !=-8888:                
                dpdp    = dpdp / 10.  # 36- 40  integer dew point depression (degrees to tenth e.g. 11=1.1 C)    
            else:
                dpdp = np.nan 

            wdir    = int(line[40:45])        # 41- 45  integer wind direction (degrees from north, 90 = east)
            if wdir == -8888 or wdir == -9999 :
                wdir = np.nan         

            wspd    = int(line[46:51])   # 47- 51  integer wind speed (meters per second to tenths, e.g. 11 = 1.1 m/s  [m/s]
            if wspd != -8888 and wspd != -9999 :
                wspd = wspd / 10.  
            else:
                wspd = np.nan                  
            if reltime == 9999.0:
                reltime = np.nan 

            z_type = np.nan
            if not (np.isnan(press)):
                z_type = 1
            elif  (np.isnan(press) and  not np.isnan(gph) ) :
                z_type = 2     

            for value,var in zip([gph, temp, wspd, wdir, rh, dpdp],  ['gph', 'temperature', 'wind_speed', 'wind_direction', 'relative_humidity' , 'dew_point_depression'] ):
                obs_id = obs_id +1 
                if not np.isnan(press):     # when pressure is available, z_coord== pressure and z_type==1  
                    z_type = 1            
                    z_value = press
                    #read_data.append ( ( 'IGRA2'.rjust(10), head_count,  int(obs_id),  idate, iday, ident, lat, lon, press, value, cdmvar_dic[var]['cdm_var'], int(cdmvar_dic[var]['cdm_unit']), numlev, z_type, release_time ) )
                elif  (np.isnan(press) and  not np.isnan(gph) ) :  # when pressure is not available, z_coord== gph and z_type==2 
                    z_type = 2       
                    z_value = gph
                    #read_data.append ( ( 'IGRA2'.rjust(10), head_count,  int(obs_id),  idate, iday, ident, lat, lon, gph, value, cdmvar_dic[var]['cdm_var'], int(cdmvar_dic[var]['cdm_unit']), numlev, z_type, release_time ) )
                else:
                    z_type = -2147483648     
                    z_value = press
                    #read_data.append ( ( 'IGRA2'.rjust(10), head_count,  int(obs_id),  idate, iday, ident, lat, lon, press, value, cdmvar_dic[var]['cdm_var'], int(cdmvar_dic[var]['cdm_unit']), numlev, z_type, release_time ) )
                read_data.append ( ( 'IGRA2'.rjust(10), head_count,  int(obs_id),  idate, iday, ident, lat, lon, z_value, value, cdmvar_dic[var]['cdm_var'], int(cdmvar_dic[var]['cdm_unit']), 
                                     numlev, z_type, release_time, report_timeflag, lvltyp2 ) )


    column_names_igra2 = [ 'source_id', 'report_id',  'observation_id', 'record_timestamp' , 'iday', 'statid@hdr', 'lat@hdr', 'lon@hdr', 'vertco_reference_1@body',
                           'obsvalue@body', 'varno@body' , 'units',  'number_of_pressure_levels', 'vertco_type@body', 'report_timestamp', 'report_meaning_of_timestamp',
                           'observation_height_above_station_surface']

    df = pd.DataFrame(data= read_data, columns= column_names_igra2)

    # set obs height a.s.s. to nan, where it's unknown. 
    # df.observation_height_above_station_surface.replace(0, np.nan, inplace=True)
    df.loc[df.observation_height_above_station_surface != 1, 'observation_height_above_station_surface'] = np.nan


    df['observation_id']  = np.chararray.zfill( (df['observation_id'].astype(int)) .astype('S'+str(id_string_length ) ), id_string_length  )  #converting to fixed length bite objects 
    df['report_id']           = np.chararray.zfill( (df['report_id'].astype(int)).astype ('S'+str(id_string_length ) ), id_string_length  )

    col = [c for c in df.columns if '_id' not in c]
    df_new = df[col].replace([-999.9, -9999, -999, -999.0, -99999.0, -99999.9, 99999.0,  -99999.00 ], np.nan)    
    for c in ['observation_id', 'report_id']:
        df_new[c] = df[c]

    #df_new = df.sort_values(by = ['record_timestamp', 'vertco_reference_1@body' ] )    # FF check here !!!! 

    print(np.unique(df_new['report_meaning_of_timestamp'] ) )

    #df_new = df_new[:5000]  # TO DO HERE TO DO CHANGE!!!
    return df_new, stations_id



def make_odb_header(odbfile, dataset):
    """ Create the header from the odb file, if not found in the 'headers/' directory.
          Headers contain the columsn names and their respective variable types.
          Only for ODB files. """

    if 'mobile' in dataset:
        dataset = dataset.replace('_mobile', '')

    #if dataset in ('era5_1', 'era5_2'):
    #    year = odbfile.split('.conv.')[1][0:4]
            
            
    header = '/srvfs/home/uvoggenberger/CEUAS/CEUAS/public/harvest/code_cop2/headers/' + dataset + '_header.dat' ## changed path to hard
    
    if dataset in ('era5_1'):
        year = odbfile.split('.conv.')[1][0:4]
        header = header.replace('.dat' , '_'+year+'.dat')
        
        
    if not os.path.isfile ( header ):
        print(' Creating the header file for the dataset: ', dataset )
        if dataset in ('era5_1'):

            odbfile = odbfile.replace('.gz','').replace('.txt','')
        else:
            odbfile = odbfile.replace('.gz','').replace('.conv._','.conv.')
        
        if "new" in odbfile:
            odbfile = odbfile.replace('new/', '')

        rdata=subprocess.check_output(["odc","header",  odbfile ])

        with open( header , 'wb' ) as f:
            f.write(rdata) 

        f = open(header , 'rb')
        rdata=f.read()
        rdata=rdata.decode('utf-8').split('\n')   

    else:
        f = open(header , 'rb')
        rdata=f.read()
        rdata=rdata.decode('utf-8').split('\n')
        #print(' Done reading the existing header file for the dataset: ', dataset )

    columns, kinds, tdict =[] , [] , {} 

    for r in rdata[2:-2]:
        try:

            if r[:6]=='Header':
                break
            else:    
                columns.append(r.split('name: ')[1].split(',')[0])
                kinds.append(r.split('type: ')[1].split(',')[0])
                if kinds[-1]=='REAL':
                    tdict[columns[-1]]=numpy.float32
                elif 'INTEGER' in kinds[-1] or 'BITFIELD' in kinds[-1]:
                    #print(columns[-1])
                    if columns[-1]=='sonde_type@conv' or columns[-1]=='station_type@conv':
                        tdict[columns[-1]]=numpy.float32
                    else: 
                        #tdict[columns[-1]]=numpy.int32  # 
                        tdict[columns[-1]]='Int64' # needed when pandas read a csv that has to fill with NANs
                else:
                    tdict[columns[-1]]=numpy.dtype('S') # dict containng column name and type

        except IndexError:
            pass       

    """ This is done otherwise for the era5 databases (1759,1761,3188) the tdict has different length than the columns list.
          So the following call alldict=pd.read_csv(f,delimiter='\t', usecols=columns, quoting=3,comment='#', skipinitialspace=True, dtype=tdict) breaks  """    
    for t in tdict.keys():
        if t not in columns:
            del tdict[t]

    """ These values must be removed rom the fb, since they have NULL values and it creates problem with 
    alldict=pd.read_csv(f,delimiter='\t', usecols=columns, quoting=3,comment='#', skipinitialspace=True, dtype=tdict) """                  

    if dataset in ["era5_1759", "era5_1761", "era5_3188"]:
        remove = ['sonde_type@conv' , "eda_spread@errstat", "bias_volatility@body" , "timeseries_index@conv"]
        for c in remove:
            #print("Removing wrong fb column: " , c)
            try:
                columns.remove(c)
                del tdict[c]
            except:
                pass
    return columns, kinds, tdict


def read_all_odbsql_stn_withfeedback(dataset, odbfile):
    """ Read the text gzip files, created from the original ODB files. """
    columns, kinds, tdict = make_odb_header(odbfile, dataset)      

    try:            
        t=time.time()                 
        try:
            # f=gzip.open(odbfile) 
            f=odbfile                
        except:
            print(odbfile, 'The zipped ODB file was not found !')
            return

        tdict['sensor@hdr']=numpy.float32
        tdict['ppcode@conv_body']=numpy.float32
        tdict['expver']=numpy.dtype('S')

        # restrict feedback to certain columns        
        #for c in columns:
        #    if c not in d:
        #        del tdict[c]

        if 'mobile' not in dataset:
            # mobile dataset have no column names, so you have to assign them manually
            alldict=pd.read_csv(f,delimiter='\t', usecols=columns, quoting=3,comment='#', skipinitialspace=True, dtype=tdict,compression='gzip' ) #nrows=1000000) # 
            alldict.to_csv('prova_stat_74005.csv', sep = '\t')
            
        else:
            alldict=pd.read_csv(f,delimiter='\t', names=columns, quoting=3,comment='#', skipinitialspace=True, skiprows=1, index_col=False, compression='gzip') #nrows=1000000) #  dtype=tdict, 

            
        """ Case where erafb is not available """
        if 'fg_depar@body' not in columns:
            alldict['fg_depar@body']=numpy.float32(numpy.NaN)
            alldict['an_depar@body']=numpy.float32(numpy.NaN)
            alldict['biascorr@body']=numpy.float32(numpy.NaN)
            alldict['sonde_type@conv']=numpy.int32(-2147483648)
            alldict['reportype']=numpy.int32(-2147483648)

        #print(time.time()-t,sys.getsizeof(alldict)//1024//1024)
        idx=numpy.where(numpy.logical_or(alldict.reportype.values==16045,alldict.reportype.values==16068))[0]

        if len(idx)>0:

        #alldict.drop(index=alldict.index[idx],inplace=True)
            y=numpy.int64(alldict['date@hdr'].values)*1000000+alldict['time@hdr'].values
            x=numpy.unique(y)
            dropindex=[]
            for i in range(1,x.shape[0]):
                if x[i]-x[i-1]<60:
                    idx=numpy.where(y==x[i-1])[0]
                    if idx.shape[0]>0:
                        dropindex.append(idx)
                    else:
                        print('empty index')
            if dropindex:          
                dropindex = numpy.concatenate(dropindex).ravel()
                alldict.drop(index=alldict.index[dropindex],inplace=True)

        alldict['source_id'] = dataset.rjust(10)

        for c in alldict.columns:                
            if type(alldict[c].iloc[0]) in [str,bytes]:
                l=alldict[c].shape[0]
                slen=len(alldict[c].values[0])
                alldict[c]=numpy.array(alldict.pop(c).values,dtype='S{}'.format(slen))
                #alldict[c]=np.bytes_(alldict[c])

            if type(alldict[c].iloc[0]) is numpy.int64:
                try:
                    alldict[c]=numpy.int32(alldict[c]) 
                except:
                    pass

            if type(alldict[c].iloc[0]) is numpy.float64:
                try:
                    alldict[c]=numpy.float32(alldict[c])
                except:
                    pass
                
        #print('after odb:',time.time()-t)

    #except MemoryError:
    except:
        print('Reading ODB failed !  ' + odbfile)
        return None

    #print(odbfile,time.time()-t)#, sys.getsizeof(alldict))
    return alldict


def fromfb_l(fbv,di, cdmfb,cdmkind):
    x=0
    if type(cdmfb) is list:
        if len(cdmfb)==3:
            if cdmfb[2] in fbv.keys():    
                x=cdmfb[0](fbv[cdmfb[1]],fbv[cdmfb[2]])
            else:
                x=cdmfb[0](fbv[cdmfb[1]],di[cdmfb[2]])
        else:
            x=cdmfb[0](fbv[cdmfb[1]])

    else:    
        x=fbv[cdmfb]

    return x


def hdrfromfb(fbv,di,cdmfb,cdmkind):
    x=0
    if type(cdmfb) is list:
        if len(cdmfb)==3:
            x=cdmfb[0](fbv[cdmfb[1]],fbv[cdmfb[2]])
            if di[list(di.keys())[0]].shape[0]<fbv.shape[0]:
                x=numpy.unique(x)

        else:
            if cdmfb[1] in fbv.keys():    
                x=cdmfb[0](fbv[cdmfb[1]])
            else:
                x=cdmfb[0](di[cdmfb[1]])

    else:
        x=np.array(fbv[cdmfb])[di['recordindex'].values]

    return x


def ttrans(cdmtype, kinds=kinds):
    """ convert the cdm types to numpy types """    
    nptype=numpy.float32
    try:
        nptype=kinds[cdmtype.strip()]
    except:
        pass
        #print(cdmtype,'not found, using numpy.float32')   
    return nptype


def make_datetime_indices(date_time):
    """ Extracts the list of observation dates, and store the indices of the first and last observations 
          Args:: 
                  list of *sorted* observation date_times
          Return:: numpy array (3,:) where the first list is the list of unique obs date_times, the second are the indices of first obs  date_time, the third the list of last obs date_time. 
          """            

    date_times, indices, counts = numpy.unique(date_time, return_counts = True, return_index= True)

    return indices , date_times,  counts  


@njit(cache=True)
def find_dateindex_l(y,x):
    #x=y#numpy.unique(y)
    z=numpy.zeros((3,x.shape[0]),dtype=numpy.int32)
    z-=1
    j=0
    for i in range(len(y)):
        m=y[i]
        if x[j]==y[i]:
            if z[1,j]==-1:
                z[1,j]=i
                #print(j,i)
            else:
                if z[2,j]<i:
                    z[2,j]=i
        elif x[j]<y[i]:
            j+=1
            if x[j]==y[i]:
                if z[1,j]==-1:
                    z[1,j]=i
                    #print(j,i)
                else:
                    if z[2,j]<i:
                        z[2,j]=i
            else:
                print('Error')
        else:
            j-=1
            if x[j]==y[i]:
                if z[1,j]==-1:
                    z[1,j]=i
                    #print(j,i)
                else:
                    if z[2,j]<i:
                        z[2,j]=i
            else:
                print('Error')
    z[0,:]=x
    return z

@njit(cache=True)
def find_recordindex_l(y,x):
    #x=y#numpy.unique(y)
    z=numpy.zeros((3,x.shape[0]),dtype=numpy.int32)
    z-=1
    j=0
    for i in range(len(y)):
        m=y[i]
        if x[j]==y[i]:
            if z[1,j]==-1:
                z[1,j]=i
                #print(j,i)
            else:
                if z[2,j]<i:
                    z[2,j]=i
        elif x[j]<y[i]:
            j+=1
            if x[j]==y[i]:
                if z[1,j]==-1:
                    z[1,j]=i
                    #print(j,i)
                else:
                    if z[2,j]<i:
                        z[2,j]=i
            else:
                print('Error')
        else:
            j-=1
            if x[j]==y[i]:
                if z[1,j]==-1:
                    z[1,j]=i
                    #print(j,i)
                else:
                    if z[2,j]<i:
                        z[2,j]=i
            else:
                print('Error')
    #z[0,:]=x
    return z


# write_dict_h5(fno, groups[k], k, groupencodings[k], var_selection=[],mode='a', attrs={'date_time':('units','seconds since 1900-01-01 00:00:00')})

def write_dict_h5(dfile, f, k, fbencodings, var_selection=[], mode='a', attrs={}, chunksize=100000): 
    """ Writes each separate variable from the observation or feedback tables inot netcdf using h5py.
          f is a pandas dataframe with one column, one for each variable
          k is either 'era5fb' or 'observations_table'
          fbencodings is the encodings of variable types, e.g. {'observations_id': { 'compression': 'gzip' } ,'compression_opts': 4 }} or
          {'observations_id': { 'compression': 32015 } ,'compression_opts': 3 }}

          I THINK IT IS WRONG...SHOULD BE          
          {'observations_id': { 'compression': 32015 ,'compression_opts': 3 }}
          ...


          attrs to set variable attributes
          mode can be 'a' or 'w'
          chunksize is set by default to 100000. Auto is not a good choice especially for character variables. Those have chunksize (chunksize,strlen)
    """


    """
    TO DO 
    write_dict_h5(dfile, f, k, fbencodings={'observations_id': { 'compression': 32015 } ,'compression_opts': 3 }}...)
    """


    if os.path.isfile( 'encodings.txt'  ):
        lines =  open('encodings.txt' , 'r').readlines()
        all_vars = [ l.split('\t')[0] for l in lines]
        all_vars = list (np.unique(all_vars) )
    else:
        all_vars = []


    if isinstance(dfile, h5py._hl.files.File):
        fd = dfile
    else:
        os.makedirs(os.path.dirname(dfile), exist_ok=True)
        fd = h5py.File(dfile,mode)

        if(fbencodings):
            for v in fbencodings.values():
                if 'compression' in v.keys():
                    comp = v['compression']
                else:
                    comp = None
                if 'compression_opts' in v.keys():
                    compopt = v['compression_opts']
                else:
                    compopt = None
        else:
            comp = None
            compopt = None

        if fd:
            if not var_selection:
                var_selection=list(f.keys())

            try:
                # variable 'index' is needed only to attach a scale to it (needed for netCDF compatibility)
                fd.create_group(k)
                idl = f[list(f.keys())[0]].shape[0]
                index=numpy.zeros (idl, dtype=np.int8)
                fd[k].create_dataset('index', data=index, compression=comp, compression_opts=compopt, chunks=True)
            except:
                pass

            sdict={}
            slist=[]

            for v in var_selection:

                if v in ['source_file', 'observation_id', 'report_id']:
                    string_length = numpy.zeros(200,dtype='S1')

                else:
                    string_length =numpy.zeros(fixed_string_len,dtype='S1')  ### fixed_string_len=10

                #if v == 'sensor_id':
                #    a = 0
                if v in [ 'report_event1@hdr' , 'report_rdbflag@hdr' , 'datum_anflag@body', 'datum_event1@body', 'datum_rdbflag@body']:
                    continue 
                if v in [ 'index']:
                    continue                 


                if type(f[v]) == pd.core.series.Series:
                    fvv=f[v].values
                else:
                    fvv=f[v]
                    
                if v not in all_vars:
                    out_en = open('encodings.txt' , 'a+')
                    out_en.write(v + '\t' + k + '\t' +  str(type(fvv[0]) ) + '\n')
                    out_en.close()
                    all_vars.append(v)
                        
                try:
                    if fvv.dtype ==pd.Int64Dtype(): ### TO DO 
                        continue
                except:
                    pass

                if type(fvv[0]) not in [str,bytes,numpy.bytes_]:  ### HORRIBLE HANDLING of types, dtypes, strings, bytes...          
                    if fvv.dtype !='S1':
                        try:
                            if fvv.shape[0] > chunksize:
                                fd[k].create_dataset(v,data=fvv,compression=comp,compression_opts=compopt,
                                                     chunks=(chunksize, ))                            
                            else:
                                fd[k].create_dataset(v,data=fvv)
                        except:
                            #print('except',dfile, k, v, fd[k].keys())
                            fd[k].create_dataset(v,data=fvv, chunks=True)

                        '''
                        try:
                            fd[k][v][:]=fvv[:]
                        except:
                            fd[k][v][:] = np.empty( (len( fvv)) )
                        '''    
                        if attrs:   
                            if v in attrs.keys():
                                for kk,vv in attrs[v].items():
                                    if type(vv) is str:  
                                        fd[k][v].attrs[kk]=numpy.bytes_(vv)
                                    else:
                                        fd[k][v].attrs[kk]=vv

                        if v in ['date_time','report_timestamp','record_timestamp']:
                            fd[k][v].attrs['units']=numpy.bytes_('seconds since 1900-01-01 00:00:00')                            

                    else:
                        if fvv.shape[0] > chunksize:                       
                            fd[k].create_dataset(v,data=fvv,compression=comp,compression_opts=compopt,
                                                 chunks=(chunksize, fvv.shape[1]))
                        else:
                            fd[k].create_dataset(v,data=fvv)

                        slen=fvv.shape[1]
                        sdict[v]=slen
                        if slen not in slist:
                            slist.append(slen)
                            try:
                                fd[k].create_dataset( 'string{}'.format(slen),  data=string_length[:slen]  )
                            except:
                                pass               
                        if v in attrs.keys():
                            for kk,vv in attrs[v].items():
                                if type(vv) is str:  
                                    fd[k][v].attrs[kk]=numpy.bytes_(vv)
                                else:
                                    fd[k][v].attrs[kk]=vv
                            # fd[k][v].attrs['description']=numpy.bytes_(attrs[v]['description'])
                            # fd[k][v].attrs['external_table']=numpy.bytes_(attrs[v]['external_table'])

                else:
                    sleno=len(fvv[0])
                    slen=sleno
                    try:
                        slen=int(fvv.dtype.descr[0][1].split('S')[1])
                    except:  
                        slen=15

                    sdict[v]=slen
                    if slen not in slist:
                        slist.append(slen)

                        try:
                            fd[k].create_dataset( 'string{}'.format(slen),  data=string_length[:slen]  )
                        except:
                            pass          

                    try:
                        if fvv.shape[0] > chunksize:
                            fd[k].create_dataset(v,data=fvv.view('S1').reshape(fvv.shape[0],slen),compression=comp,
                                                 compression_opts=compopt, chunks=(chunksize, slen))
                        else:
                            fd[k].create_dataset(v,data=fvv)

                    except:
                        #fd[k].create_dataset(v,data=np.bytes_(fvv).view('S1').reshape(fvv.shape[0],slen),compression=fbencodings[v]['compression'],chunks=True)                    
                        pass
                    if v in attrs.keys():
                        for kk,vv in attrs[v].items():
                            if type(vv) is str:  
                                fd[k][v].attrs[kk]=numpy.bytes_(vv)
                            else:
                                fd[k][v].attrs[kk]=vv
                        # fd[k][v].attrs['description']     =numpy.bytes_(attrs[v]['description'])
                        # fd[k][v].attrs['external_table']=numpy.bytes_(attrs[v]['external_table'])                


                for v in fd[k].keys(): #var_selection:
                    if 'string'  in v or v== 'index' :                    
                        continue 
                    if v not in  f.keys():
                        continue

                    try:
                        if type(f[v]) == pd.core.series.Series:
                            fvv=f[v].values
                        else:
                            fvv=f[v]

                        fd[k][v].dims[0].attach_scale(fd[k]['index'])

                        if fvv.ndim==2 or type(fvv[0]) in [str,bytes,numpy.bytes_]:
                            slen=sdict[v]
                            #slen=10
                            fd[k][v].dims[1].attach_scale(fd[k]['string{}'.format(slen)])

                    except:
                        pass




        for v in slist:
            s='string{}'.format(v)
            for a in ['NAME']:
                fd[k][s].attrs[a]=numpy.bytes_('This is a netCDF dimension but not a netCDF variable.')


    return


def write_dict_h5_old_leo(dfile, f, k, fbencodings, var_selection=[], mode='a', attrs={}): 
    """ Writes each separate variable from the observation or feedback tables inot netcdf using h5py.
          f is a pandas dataframe with one column, one for each variable
          k is either 'era5fb' or 'observations_table'
          fbencodings is the encodings of variable types, e.g. {'observations_id': { 'compression': 'gzip' } }
    """

    #attrs=  {'date_time':('units','seconds since 1900-01-01 00:00:00')}
    #attrs = {'observation_id': ('description', 'unique ID for observation'), 'report_id': ('description', 'Link to header information') , 'date_time':('units','seconds since 1900-01-01 00:00:00') }

    with h5py.File(dfile,mode) as fd:
        try:
            fd.create_group(k)
            index=numpy.zeros (f[list(f.keys())[0]].shape[0], dtype='S1')
            fd[k].create_dataset('index', data=index)
        except:
            pass
        if not var_selection:
            var_selection=list(f.keys())

        string10=numpy.zeros(fixed_string_len,dtype='S1')
        sdict={}
        slist=[]

        for v in var_selection:
            if v == 'source_file':
                a = 0
            if v in [ 'report_event1@hdr' , 'report_rdbflag@hdr' , 'datum_anflag@body', 'datum_event1@body', 'datum_rdbflag@body', 'index' , 'varbc_ix@body']:
                a=0
                #continue 

            if type(f[v]) == pd.core.series.Series:
                fvv=f[v].values
            else:
                fvv=f[v]

            try:
                if fvv.dtype ==pd.Int64Dtype(): ### TO DO 
                    continue
            except:
                pass
            
            if v == 'source_id':
                x = 0


            if type(fvv[0]) not in [str,bytes,numpy.bytes_]:  ### HORRIBLE HANDLING of types, dtypes, strings, bytes... 
                #print(v, '  ', type(fvv[0]) , '  ' , fvv.dtype )

                if fvv.dtype !='S1':
                    #if fvv.dtype == "Int64":
                    #    0
                    #vtype = np.int32
                    #else:
                    #    vtype = fvv.dtype

                    try:
                        fd[k].create_dataset(v,fvv.shape,fvv.dtype, **fbencodings[v])
                    except:
                        #fd[k].create_dataset(v,fvv.shape,'int32',compression=fbencodings[v]['compression'], chunks=True)  
                        fd[k].create_dataset(v,fvv.shape,fvv.dtype,compression='gzip', chunks=True)

                    try:
                        fd[k][v][:]=fvv[:]
                    except:
                        fd[k][v][:] = np.empty( (len( fvv)) )

                    if attrs:    #  attrs={'date_time':('units','seconds since 1900-01-01 00:00:00')}
                        if v in attrs.keys():
                            for kk,vv in attrs[v].items():
                                if type(vv) is str:  
                                    fd[k][v].attrs[kk]=numpy.bytes_(vv)
                                else:
                                    fd[k][v].attrs[kk]=vv

                    if v in ['date_time','report_timestamp','record_timestamp']:
                        fd[k][v].attrs['units']=numpy.bytes_('seconds since 1900-01-01 00:00:00')                            #print (  fk, ' ' , v , ' ' ,   ) 

                else:
                    fd[k].create_dataset(v,fvv.shape,fvv.dtype, **fbencodings[v])
                    fd[k][v][:]=fvv[:]
                    slen=fvv.shape[1]
                    sdict[v]=slen
                    if slen not in slist:
                        slist.append(slen)
                        try:
                            fd[k].create_dataset( 'string{}'.format(slen),  data=string10[:slen]  )
                        except:
                            pass               
                    if v in attrs.keys():
                        fd[k][v].attrs['description']=numpy.bytes_(attrs[v]['description'])
                        fd[k][v].attrs['external_table']=numpy.bytes_(attrs[v]['external_table'])

            else:
                sleno=len(fvv[0])
                slen=sleno
                try:
                    slen=int(fvv.dtype.descr[0][1].split('S')[1])
                except:  
                    slen=15

                sdict[v]=slen
                if slen not in slist:
                    slist.append(slen)
                    try:
                        fd[k].create_dataset( 'string{}'.format(slen),  data=string10[:slen]  )
                    except:
                        pass               
                try:

                    fd[k].create_dataset(v,data=fvv.view('S1').reshape(fvv.shape[0],slen),**fbencodings[v])
                except KeyError:
                    fd[k].create_dataset(v,data=fvv.view('S1').reshape(fvv.shape[0],slen),compression='gzip',chunks=True)

                if v in attrs.keys():
                    fd[k][v].attrs['description']     =numpy.bytes_(attrs[v]['description'])
                    fd[k][v].attrs['external_table']=numpy.bytes_(attrs[v]['external_table'])                


            #variables_dic[v] = f[v].values.dtype

        for v in fd[k].keys(): #var_selection:
            l=0      
            if 'string'  in v or v== 'index' :                    
                continue 
            try:
                if type(f[v]) == pd.core.series.Series:
                    fvv=f[v].values
                else:
                    fvv=f[v]
                fd[k][v].dims[l].attach_scale(fd[k]['index'])
                #print(v,fvv.ndim,type(fvv[0]))
                if fvv.ndim==2 or type(fvv[0]) in [str,bytes,numpy.bytes_]:
                    slen=sdict[v]
                    #slen=10
                    fd[k][v].dims[1].attach_scale(fd[k]['string{}'.format(slen)])
            except:
                pass

        #i=4        
        for v in slist:
            s='string{}'.format(v)
            for a in ['NAME']:
                fd[k][s].attrs[a]=numpy.bytes_('This is a netCDF dimension but not a netCDF variable.')
            #i+=1

    return
    
    
def write_dict_h5_new(dfile, f, k, fbencodings, var_selection=[], mode='a', attrs={}, chunksize=100000): 
    """ Writes each separate variable from the observation or feedback tables inot netcdf using h5py.
          f is a pandas dataframe with one column, one for each variable
          k is either 'era5fb' or 'observations_table'
          fbencodings is the encodings of variable types, e.g. {'observations_id': { 'compression': 'gzip' } ,'compression_opts': 4 }} or
          {'observations_id': { 'compression': 32015 } ,'compression_opts': 3 }}
          attrs to set variable attributes
          mode can be 'a' or 'w'
          chunksize is set by default to 100000. Auto is not a good choice especially for character variables. Those have chunksize (chunksize,strlen)
    """

    #attrs=  {'date_time':('units','seconds since 1900-01-01 00:00:00')}
    #attrs = {'observation_id': ('description', 'unique ID for observation'), 'report_id': ('description', 'Link to header information') , 'date_time':('units','seconds since 1900-01-01 00:00:00') }
    
    #tt = time.time()
    if isinstance(dfile, h5py._hl.files.File):
        fd = dfile
    else:
        os.makedirs(os.path.dirname(dfile), exist_ok=True)
        fd = h5py.File(dfile,mode)
            
            
            
#    with h5py.File(dfile,mode) as fd:
    if(fbencodings):
        for v in fbencodings.values():
            if 'compression' in v.keys():
                comp = v['compression']
            else:
                comp = None
            if 'compression_opts' in v.keys():
                compopt = v['compression_opts']
            else:
                compopt = None
    else:
        comp = None
        compopt = None
        
    if fd:
        if not var_selection:
            var_selection=list(f.keys())

        try:
            # variable 'index' is needed only to attach a scale to it (needed for netCDF compatibility)
            fd.create_group(k)
            idl = f[list(f.keys())[0]].shape[0]
            index=numpy.zeros (idl, dtype=np.int8)
            fd[k].create_dataset('index', data=index, compression=comp,compression_opts=compopt, chunks=True)
        except:
            pass
        
        string10=numpy.zeros(fixed_string_len,dtype='S1')
        sdict={}
        slist=[]

        #groupencodings     
        
        #print('start', time.time() - tt)
        for v in var_selection:          
            #variables_dic[v] = ''
            if type(f[v]) == pd.core.series.Series:
                fvv=f[v].values
            else:
                fvv=f[v]
                
            if type(fvv[0]) not in [str,bytes,numpy.bytes_]:

                if fvv.dtype !='S1':
                    try:
                        
                        #fd[k].create_dataset(v,fvv.shape,fvv.dtype,compression=fbencodings[v]['compression'], chunks=True)
                        if True or fvv.shape[0] > chunksize:                       
                            fd[k].create_dataset(v,data=fvv,compression=comp,compression_opts=compopt,
                                                 chunks=(np.min([chunksize, fvv.shape[0]]), ))
                                                 #chunks=(np.int32(np.sqrt(fvv.shape[0]))*10, ))
                        else:    
                            fd[k].create_dataset(v,data=fvv)
                    except:
                        print('except',dfile, k, v, fd[k].keys())
                        try:
                            fd[k].create_dataset(v,data=fvv, chunks=True)
                        except:
                            pass
                        
                    #fd[k][v][:]=fvv[:]
                    if attrs:    #  attrs={'date_time':('units','seconds since 1900-01-01 00:00:00')}
                        if v in attrs.keys():
                            for kk,vv in attrs[v].items():
                                if type(vv) is str:  
                                    fd[k][v].attrs[kk]=numpy.bytes_(vv)
                                else:
                                    fd[k][v].attrs[kk]=vv
                                                                
                    if v in ['date_time','report_timestamp','record_timestamp']:
                        fd[k][v].attrs['units']=numpy.bytes_('seconds since 1900-01-01 00:00:00')                            #print (  fk, ' ' , v , ' ' ,   ) 
                                
                else:
                    #fd[k].create_dataset(v,fvv.shape,fvv.dtype,compression=fbencodings[v]['compression'], chunks=True)
                    #fd[k][v][:]=fvv[:]
                    if True or fvv.shape[0] > chunksize:                       
                        fd[k].create_dataset(v,data=fvv,compression=comp,compression_opts=compopt,
                                             chunks=(np.min([chunksize, fvv.shape[0]]), fvv.shape[1]))
                    else:
                        fd[k].create_dataset(v,data=fvv)
                        
                    #fd[k][v][:]=fvv[:]
                    slen=fvv.shape[1]
                    sdict[v]=slen
                    if slen not in slist:
                        slist.append(slen)
                        try:
                            fd[k].create_dataset( 'string{}'.format(slen),  data=string10[:slen]  )
                        except:
                            pass               
                    if v in attrs.keys():
                        fd[k][v].attrs['description']=numpy.bytes_(attrs[v]['description'])
                        fd[k][v].attrs['external_table']=numpy.bytes_(attrs[v]['external_table'])
                        
            else:
                sleno=len(fvv[0])
                slen=sleno
                try:
                    slen=int(fvv.dtype.descr[0][1].split('S')[1])
                except:  
                    pass

                sdict[v]=slen
                if slen not in slist:
                    slist.append(slen)
                    try:
                        fd[k].create_dataset( 'string{}'.format(slen),  data=string10[:slen]  )
                    except:
                        pass               
                try:
                    if v == 'observation_id':
                        x = 0
                    
                    fd[k].create_dataset(v,data=fvv.view('S1').reshape(fvv.shape[0],slen),compression=comp,
                                         compression_opts=compopt, chunks=(np.min([chunksize, fvv.shape[0]]), slen))
                except:
                    #fd[k].create_dataset(v,data=np.bytes_(fvv).view('S1').reshape(fvv.shape[0],slen),compression=fbencodings[v]['compression'],chunks=True)                    
                    pass
                if v in attrs.keys():
                    fd[k][v].attrs['description']     =numpy.bytes_(attrs[v]['description'])
                    fd[k][v].attrs['external_table']=numpy.bytes_(attrs[v]['external_table'])                

                        
            #variables_dic[v] = f[v].values.dtype
             
            #print('v', time.time() - tt)
        for v in fd[k].keys(): #var_selection:
            l=0      
        
            try:
                if v in f.keys():
                    
                    if type(f[v]) == pd.core.series.Series:
                        fvv=f[v].values
                    else:
                        fvv=f[v]
                else:
                    continue
                if 'string' not in v and v!='index':                    
                    fd[k][v].dims[l].attach_scale(fd[k]['index'])
                    #print(v,fvv.ndim,type(fvv[0]))
                    if fvv.ndim==2 or type(fvv[0]) in [str,bytes,numpy.bytes_]:
                        slen=sdict[v]
                        #slen=10
                        fd[k][v].dims[1].attach_scale(fd[k]['string{}'.format(slen)])
            except:
                pass
            
            
            
        i=4        
        for v in slist:
            s='string{}'.format(v)
            for a in ['NAME']:
                fd[k][s].attrs[a]=numpy.bytes_('This is a netCDF dimension but not a netCDF variable.')
            
            i+=1
        #print('el', time.time() - tt)
    return


def write_dict_h5_old(dfile, f, k, fbencodings, var_selection=[], mode='a', attrs={}): 
    """ Writes each separate variable from the observation or feedback tables inot netcdf using h5py.
          f is a pandas dataframe with one column, one for each variable
          k is either 'era5fb' or 'observations_table'
          fbencodings is the encodings of variable types, e.g. {'observations_id': { 'compression': 'gzip' } }
    """

    #attrs=  {'date_time':('units','seconds since 1900-01-01 00:00:00')}
    #attrs = {'observation_id': ('description', 'unique ID for observation'), 'report_id': ('description', 'Link to header information') , 'date_time':('units','seconds since 1900-01-01 00:00:00') }

    with h5py.File(dfile,mode) as fd:
        try:
            fd.create_group(k)
            index=numpy.zeros (f[list(f.keys())[0]].shape[0], dtype='S1')
            fd[k].create_dataset('index', data=index)
        except:
            pass
        if not var_selection:
            var_selection=list(f.keys())

        string10=numpy.zeros(fixed_string_len,dtype='S1')
        sdict={}
        slist=[]

        for v in var_selection:
            #if v == 'source_file':
            #    a = 0
            #if v in [ 'report_event1@hdr' , 'report_rdbflag@hdr' , 'datum_anflag@body', 'datum_event1@body', 'datum_rdbflag@body', 'index' , 'varbc_ix@body']:
            #    a=0
                #continue 
            if v in [ 'index']:
                    continue
                    
            if type(f[v]) == pd.core.series.Series:
                fvv=f[v].values
            else:
                fvv=f[v]

            try:
                if fvv.dtype ==pd.Int64Dtype(): ### TO DO 
                    continue
            except:
                pass

            # this is needed to remove the <NA> types from pandas in the array or it crashes, due to type Int64 NAN in pandas 
            if v in [ 'report_event1@hdr' , 'report_rdbflag@hdr' , 'datum_anflag@body', 'datum_event1@body', 'datum_rdbflag@body', 'index' , 'varbc_ix@body']:
                fvv = np.array( [int_void if pd.isna(i) else i for i in fvv   ] )
            
            if type(fvv[0]) not in [str,bytes,numpy.bytes_]:  ### HORRIBLE HANDLING of types, dtypes, strings, bytes... 
                #print(v, '  ', type(fvv[0]) , '  ' , fvv.dtype )

                if fvv.dtype !='S1':
                    #if fvv.dtype == "Int64":
                    #    0
                    #vtype = np.int32
                    #else:
                    #    vtype = fvv.dtype

                    try:
                        fd[k].create_dataset(v,fvv.shape,fvv.dtype,compression=fbencodings[v]['compression'], chunks=True)
                    except:
                        #fd[k].create_dataset(v,fvv.shape,'int32',compression=fbencodings[v]['compression'], chunks=True)  
                        fd[k].create_dataset(v,fvv.shape,fvv.dtype,compression='gzip', chunks=True)

                    try:
                        fd[k][v][:]=fvv[:]
                    except:
                        fd[k][v][:] = np.empty( (len( fvv)) )

                    if attrs:    #  attrs={'date_time':('units','seconds since 1900-01-01 00:00:00')}
                        if v in attrs.keys():
                            for kk,vv in attrs[v].items():
                                if type(vv) is str:  
                                    fd[k][v].attrs[kk]=numpy.bytes_(vv)
                                else:
                                    fd[k][v].attrs[kk]=vv

                    if v in ['date_time','report_timestamp','record_timestamp']:
                        fd[k][v].attrs['units']=numpy.bytes_('seconds since 1900-01-01 00:00:00')                            #print (  fk, ' ' , v , ' ' ,   ) 

                else:
                    fd[k].create_dataset(v,fvv.shape,fvv.dtype,compression=fbencodings[v]['compression'], chunks=True)
                    fd[k][v][:]=fvv[:]
                    slen=fvv.shape[1]
                    sdict[v]=slen
                    if slen not in slist:
                        slist.append(slen)
                        try:
                            fd[k].create_dataset( 'string{}'.format(slen),  data=string10[:slen]  )
                        except:
                            pass               
                    if v in attrs.keys():
                        fd[k][v].attrs['description']=numpy.bytes_(attrs[v]['description'])
                        fd[k][v].attrs['external_table']=numpy.bytes_(attrs[v]['external_table'])

            else:
                sleno=len(fvv[0])
                slen=sleno
                try:
                    slen=int(fvv.dtype.descr[0][1].split('S')[1])
                except:  
                    slen=15

                sdict[v]=slen
                if slen not in slist:
                    slist.append(slen)
                    try:
                        fd[k].create_dataset( 'string{}'.format(slen),  data=string10[:slen]  )
                    except:
                        pass               
                try:

                    fd[k].create_dataset(v,data=fvv.view('S1').reshape(fvv.shape[0],slen),compression=fbencodings[v]['compression'],chunks=True)
                except KeyError:
                    fd[k].create_dataset(v,data=fvv.view('S1').reshape(fvv.shape[0],slen),compression='gzip',chunks=True)

                if v in attrs.keys():
                    fd[k][v].attrs['description']     =numpy.bytes_(attrs[v]['description'])
                    fd[k][v].attrs['external_table']=numpy.bytes_(attrs[v]['external_table'])                


            #variables_dic[v] = f[v].values.dtype

        for v in fd[k].keys(): #var_selection:
            l=0      
            if 'string'  in v or v== 'index' :                    
                continue 
            try:
                if type(f[v]) == pd.core.series.Series:
                    fvv=f[v].values
                else:
                    fvv=f[v]
                fd[k][v].dims[l].attach_scale(fd[k]['index'])
                #print(v,fvv.ndim,type(fvv[0]))
                if fvv.ndim==2 or type(fvv[0]) in [str,bytes,numpy.bytes_]:
                    slen=sdict[v]
                    #slen=10
                    fd[k][v].dims[1].attach_scale(fd[k]['string{}'.format(slen)])
            except:
                pass

        #i=4        
        for v in slist:
            s='string{}'.format(v)
            for a in ['NAME']:
                fd[k][s].attrs[a]=numpy.bytes_('This is a netCDF dimension but not a netCDF variable.')
            #i+=1

    return

def initialize_output(fn, output_dir, station_id, dataset, year):
    """ Simple initializer for writing the output netCDF file """
    if not os.path.isdir(output_dir + '/' + station_id):
        os.mkdir(output_dir + '/' + station_id)


    source_file = fn # fn.split('/')[-1]
    output_file = output_dir + '/' + station_id + '/' + station_id +  '_' + str(year) + '_' +  dataset + '_harvested_' + fn.split('/')[-1] + '.nc'  # creating an output file name e.g. chera5.conv._10393.nc  , try 01009 faster

    if len(source_file) < 200:
        source_file = source_file + ' '*(200-len(source_file))
    return output_file , source_file 

def convert_variable_type_n(df):
    """ Converts the variable type of the input DataFrame  (igra2,ncar,bufr) """
        # available columns
    """
    'source_file', 'source_id', 'report_id', 'observation_id',
       'record_timestamp', 'iday', 'station_id', 'lat@hdr', 'lon@hdr',
       'vertco_reference_1@body', 'obsvalue@body', 'varno@body', 'units',
       'number_of_pressure_levels'
    """
    dic_var_type = { 'int32'    : ['varno@body', 'number_of_pressure_levels' , 'units', 'z_coordinate_type' , 'vertco_type@body' ]  ,
                     'float32' : ['lat@hdr', 'lon@hdr' ,  'vertco_reference_1@body', 'obsvalue@body',  'iday'  ]  ,
                                'string'   : ['source_id' , 'station_id' ,  'source_file' , 'report_id', 'observation_id',   ] ,
                                'int64'    : ['report_timestamp' , 'date_time', 'record_timestamp'] }    

    convert = { 'int32'    : np.int32    , 
                'string'   : np.bytes_  ,
                        'float32' : np.float32 ,
                        'float64' : np.float64

                        }
    # creating a dictionary variable - nptype 
    mapping = {}
    for k in dic_var_type.keys():
        for l in dic_var_type[k]:
            mapping[l] = k 

    for c in df.columns:
        try:
            #print('converting ' , c , ' to type ' , mapping[c] )
            df[c] =  df[c].astype( convert[mapping[c]] )
            #print('converted: ', c )

        except:
            #print('could not convert type column ' , c )
            pass 

    return df


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
    # old
    #deltas =  [ dt - offset for dt in date_time ]
    #date_times_seconds  =  [ i.total_seconds() for i in deltas ] 
    #return np.array(date_times_seconds).astype(np.int64) # replacing with seconds from 1900-01-01 00:00:00 

def read_df_to_cdm(cdm, dataset, fn, metadata='' ):
    # era5 analysis feedback is read from compressed netcdf files era5.conv._?????.nc.gz in $RSCRATCH/era5/odbs/1
    """ Reading the odb and convert to xarray """  
    if 'bufr_cnr' in dataset:
        df, stations_id= read_bufr_cnr_csv(fn)
    elif 'woudc' in dataset:
        df, stations_id= read_woudc_csv(fn)
    elif  'bufr' in dataset :
        df, stations_id= bufr_to_dataframe(fn) # fdbs: the xarray converted from the pandas dataframe 
    elif  'maestro' in dataset :
        df, stations_id= bufr_to_dataframe(fn) # fdbs: the xarray converted from the pandas dataframe
    elif 'amma' in dataset:
        df, stations_id= read_amma_csv(fn)
    elif  'ncar' in dataset:
        df, stations_id = uadb_ascii_to_dataframe(fn)
    elif  'igra2' in dataset:
        df, stations_id = igra2_ascii_to_dataframe(fn)
    elif 'giub' in dataset:
        df, stations_id = read_giub(fn)
    elif 'hara' in dataset:
        df, stations_id = read_hara_csv(fn)
    elif 'shipsound' in dataset:
        df, stations_id = read_shipsound_csv(fn)
    elif 'npsound' in dataset:
        df, stations_id = read_npsound_csv(fn)   
    elif dataset == 'mauritius' in dataset :
        df, stations_id = read_mauritius_csv(fn)   
    elif dataset == 'mauritius_digitized':
        df, stations_id = read_mauritius_csv_digitized(fn)   # here fn is the name of the directory containing the files for each hum/temp sensor 
    elif dataset =='yangjiang':
        df, stations_id = read_yangjiang_csv(fn, metadata=metadata)   # here fn is the name of the directory containing the files for each hum/temp sensor 
    else:
        #print('Unidentified file is: ', fn)
        raise ValueError('Cannot identify the type of file to be analized!!! ')

    if df.empty:
        return None, None, None, None, None 


    # save = ['correct', 'wrong']                                                                                                                                                                              

    df = df.reset_index()
    df['vertco_reference_1@body'] = df['vertco_reference_1@body'].astype(float)



    station_configuration_retrieved = get_station_configuration_cuon(stations_id=stations_id, 
                                                                     station_configuration = cdm['station_configuration'],
                                                                     lat=df['lat@hdr'].values[0],  lon=df['lon@hdr'].values[0], 
                                                                      fn=fn, 
                                                                      db=dataset,
                                                                      change_lat=False )              
    # try:
    primary_id = station_configuration_retrieved['primary_id'].values[0].decode('utf-8')      
    # except:
    #     pass

    correct_data, df, most_freq_lat, most_freq_lon = check_lat_lon(df, fn, save='correct')

    try:
        sc_lat, sc_lon = station_configuration_retrieved.latitude , station_configuration_retrieved.longitude 
    except:
        a = open('logs/' + dataset + '_failed.txt' , 'a+')
        a.write(fn + '\n')
        sc_lat, sc_lon = 999 , 999
        #return None                
    # checking that coordinates in the retrieved stat_conf and file are compatible 
    stat_conf_check = True          
    
    if dataset not in ['npsound' , 'shipsound']:

        if isinstance(station_configuration_retrieved, pd.DataFrame):
            if ('20666' in primary_id) or ('20777' in primary_id) or ('20888' in primary_id) or ('20999' in primary_id):        # stations with known coordinates problems ## edit: added 20777,... also to be skipped because of orphans!
                stat_conf_check = True
            else:
                if not station_configuration_retrieved.empty:
                    sc_lat, sc_lon = station_configuration_retrieved.latitude.values[0] , station_configuration_retrieved.longitude.values[0]
    
                    if abs(sc_lat - most_freq_lat ) < 1.0 and  abs(sc_lon - most_freq_lon ) < 1.0 :
                        stat_conf_check = True
                    else:
                        stat_conf_check = False

    # try:
    primary_id = station_configuration_retrieved['primary_id'].values[0].decode('utf-8')             
    # except:
    #     pass

    '''
    try:
        primary_id = station_configuration_retrieved['primary_id'].values[0].decode('utf-8')                
    except:
        out =open('logs/' + dataset + "_wrong_ids.txt" , 'a+')
        out.write(fn + '\n')
        primary_id = '0-20999-0-' + str(stations_id[0]) + '-coordinateIssue-'
    '''
    
    if dataset in ['shipsound', 'npsound']:
        correct_data = True 
    if dataset in ['mauritius', 'mauritius_digitized' , 'yangjiang' ]:
        correct_data = True 
        stat_conf_check = True
    if '20666' in primary_id: #stations with know problems in lat and lo
        correct_data = True 
        stat_conf_check = True        


    """ Casting the original variable types to appropriate numpy types """     
    df = convert_variable_type_n(df)
    years = [str(d)[:4] for d in df['iday'] ]
    df['year'] = years

    """ Dropping columns """
    dcols=[]
    for d in df.columns:
        if d not in ['date@hdr','time@hdr','statid@hdr','vertco_reference_1@body','varno@body', 'lon@hdr','lat@hdr','seqno@hdr',
                         'obsvalue@body','fg_depar@body','an_depar@body','biascorr@body','sonde_type@conv',  'record_timestamp' , 'report_timestamp',
                               'observation_id', 'report_id' , 'units' , 'vertco_type@body' , 'year' , 'iday',  'sensor_id', 'report_meaning_of_timestamp', 
                               'observation_height_above_station_surface'] :

            dcols.append(d)

    df.drop(columns=dcols,inplace=True)

    # removing duplicated observations within one single record (same z coordinate, variable number)
    df = df.drop_duplicates(subset=['report_timestamp', 'vertco_reference_1@body', 'varno@body'] )

    # dropping NANS
    df = df.dropna(subset=['obsvalue@body'] )    

    # sorting values 
    df = df.sort_values(by = ['report_timestamp', 'vertco_reference_1@body' ] )    

    if dataset == 'mauritius_digitized':
        df.sensors = sensors 
    return df, stat_conf_check, station_configuration_retrieved, primary_id, min(years)



def write_df_to_cdm(df, stat_conf_check, station_configuration_retrieved, cdm, cdmd, out_dir, dataset, dic_obstab_attributes, fn, year, sensor=False):
    """ Convert the  pandas dataframe from the file into cdm compliant netCDF files. Use with bufr, igra2 and ncar databases.
        According to each file type, it will use the appropriate reading function to extract a Pandas DataFrame
        input:
              fn       :: odb file name (e.g. era5.conv._10393)
              cdm   :: cdm tables (read with pandas)
              cdmd :: cdm tables definitions ("")  """    

    #if dataset == 'mauritius_digitized':
    #    sensor = df.sensor 

    # chunking year 
    year = str(year)
    df =df.loc[df['year'] == year ]

    primary_id = station_configuration_retrieved.primary_id.values[0].decode('utf-8')

    fno,  source_file = initialize_output(fn, output_dir, primary_id, dataset, year)   
    log_name = fno.replace('.nc', '_correctly_processed_year.txt') # fno.replace('_'+str(year),'').replace('.nc', '_correctly_processed_year.txt')

    # splitting sensor
    #if sensor:
    #   ss = sensor.rjust(10)
    #    df = df.loc[df.sensor_id == np.bytes_(ss )]

    if df.empty:
        print('No data for year === ' , year )
        return log_name

    #fno,  source_file = initialize_output(fn, output_dir, primary_id, dataset, year)   
    #log_name = fno.replace('_'+str(year),'').replace('.nc', '_correctly_processed_year.txt')


    if dataset in ['mauritius_digitized', 'yangjiang'] and sensor:
        fno = fno.replace('intercomparison' , 'intercomparison_' + sensor )
    """ Extract the unique indices of each date observation, one for only dates, one for date_time (i.e. record index). 
        Converts the time variables in seconds since 1900-01-01 00:00:00 """  
    di=xr.Dataset() 

    df['report_timestamp'] = datetime_toseconds( df['report_timestamp'] )  #  
    df['record_timestamp'] =  datetime_toseconds( df['record_timestamp'] ) # replacing with seconds from 1900-01-01 00:00:00 

    indices, day, counts = make_datetime_indices( df['iday'].values )   #only date information
    di['dateindex']  = ( { 'dateindex' :  day.shape } , indices )          


    indices, date_times , counts  = make_datetime_indices( df['report_timestamp'].values ) #date_time plus indices           
    di['recordindex']          = ( {'recordindex' : indices.shape }, indices )
    di['recordtimestamp']  = ( {'recordtimestamp' : date_times.shape }, date_times  )

    di['recordtimestamp'].attrs['units']='seconds since 1900-01-01 00:00:00'

    di.to_netcdf( fno, format='netCDF4', engine='h5netcdf', mode='w' )


    """ Storing the variable encodings """
    '''
    fbencodings={}    # useless for non -odb files 
    for d,v in df.items():              
        if v.dtype==numpy.dtype('float64'):
                fbencodings[d]={'dtype':numpy.dtype('float32'), 'compression': 'gzip'}  

        elif v.dtype==numpy.dtype('float32'):
            fbencodings[d]={'dtype':numpy.dtype('float32'), 'compression': 'gzip'}      

        elif v.dtype==numpy.dtype('int32'):
            fbencodings[d]={'dtype':numpy.dtype('int32'), 'compression': 'gzip'} 

        elif v.dtype==numpy.dtype('int64'):
            fbencodings[d]={'dtype':numpy.dtype('int64'), 'compression': 'gzip'}                     

        elif type(v.values[0])==bytes:
                fbencodings[d]={'compression': 'gzip', 'chunksizes': ( min( [10000,v.shape[0] ] ), 10 ) }#,'chunksizes':(10000,10)
        else:
                fbencodings[d]={'compression': 'gzip'}
    '''

    # {'observations_id': { 'compression': 'gzip' } ,'compression_opts': 4 }} or
    #      {'observations_id': { 'compression': 32015 } ,'compression_opts': 3 }}

    '''
    for d,v in df.items():              
        fbencodings[d] = { d : { 'compression': 32015 } ,'compression_opts': 3 } 

    fbencodings['index']={'index' : { 'compression': 'gzip' } ,'compression_opts': 4 } 
    '''

    groups={}
    groupencodings={}

    #for k in ['observations_table']:

    encodings_variables = {} 
    for k in cdmd.keys(): # loop over all the table definitions 

        if k not in ('observations_table'):
            groups[k]=xr.Dataset() # create an  xarray
        groupencodings[k]={} # create a dict of group econding
        encodings_variables[k] = {}

        for i in range(len(cdmd[k])):
            d=cdmd[k].iloc[i]
            if d.element_name == 'observation_value':
                print()
            try:

                d_type = numpy.dtype(ttrans(d.kind,kinds=okinds) ) 
                #a = open('encodings.txt')
            except:
                print(d.element_name ,  '   WRONG dtype '     )


            encodings_variables[k][d.element_name] = {}
            encodings_variables[k][d.element_name]['type'] = d_type 
            """ Filling the observations_table """
            if k in ('observations_table'):
                groups[k]=pd.DataFrame()  # creating dataframes that will be written to netcdf via h5py methods by the write_dict_h5() method 
                if d.element_name == 'observation_height_above_station_surface':
                    print()
                try:         
                    groups[k][d.element_name]= fromfb_l(df, di._variables, cdmfb_noodb[d.element_name], ttrans(d.kind,kinds=okinds))                                                       
                except KeyError:

                    print('Missing::: ' ,  d.element_name , ' in the  obsTab' , '  ' ,  d_type   )

                    if d_type == np.dtype('O'):
                        print(d.element_name , '  ' , d_type )
                        d_type = np.int32
                    if d_type == np.int32:
                        x=numpy.zeros( df['report_timestamp'].shape[0], dtype= d_type)
                        x.fill(int_void)
                    else:
                        x=numpy.zeros( df['report_timestamp'].shape[0], dtype= d_type  )
                        x.fill(numpy.nan)
                    groups[k][d.element_name]=x


            elif k in ('header_table'):
                if 'latitude' in d.element_name: #'report_synoptic_time' in d.element_name:
                    a = 0
                if d.element_name in  ['record_timestamp', 'report_timestamp']:
                    groups[k][d.element_name]= ( {'hdrlen':di['recordindex'].shape[0]} , df[d.element_name].values[di['recordindex'].values]) 
                    groups[k][d.element_name].attrs['units'] = 'seconds since 1900-01-01 00:00:00'    

                elif d.element_name  == 'report_meaning_of_timestamp' and d.element_name in df.columns:
                    groups[k][d.element_name]= ( {'hdrlen':di['recordindex'].shape[0]} , list(df[d.element_name].iloc[di['recordindex']]))                     
                else:
                    try:  # WHY IS THIS (SO MESSY)
                        # basically no variable in the header table is truly set with the exception og the ones above OR with the ones already present in the station_configuration
                        # so you have to fill them as empty
                        # it is very badly coded... 

                        if d.element_name not in station_configuration_retrieved.columns: # variables might be from the df (input file) or from the retrieved station configuration  
                            try:          
                            #  groups[k][d.element_name] = ({'hdrlen':di['recordindex'].shape[0]}, hdrfromfb(df,di._variables, cdmfb[d.element_name],ttrans(d.kind,kinds=gkinds) ) )
                                groups[k][d.element_name]= (di['recordindex'].shape[0], hdrfromfb(df, di._variables, cdmfb_noodb[d.element_name], ttrans(d.kind,kinds=gkinds) ) )
                            except:  
                                if d.element_name == 'report_synoptic_time':
                                    x=numpy.zeros(di['recordindex'].shape[0], dtype=numpy.int32)
                                    x.fill(int_void)
                                else:
                                    x=numpy.zeros(di['recordindex'].shape[0], dtype=numpy.dtype(ttrans(d.kind,kinds=gkinds)))
                                    x.fill(numpy.nan)
                                groups[k][d.element_name]=({'hdrlen':di['recordindex'].shape[0]},x)                                   
                        else:            
                            x=numpy.zeros(di['recordindex'].shape[0], dtype=numpy.dtype(ttrans(d.kind,kinds=gkinds)))
                            x.fill( station_configuration_retrieved[d.element_name].values[0] )
                            groups[k][d.element_name]= ({'hdrlen':di['recordindex'].shape[0]},x) 
                    except: # in case I cannot retrieve the station configuration file 
                        try:                        
                            groups[k][d.element_name] = (di['recordindex'].shape[0], hdrfromfb(df, di._variables, cdmfb_noodb[d.element_name], ttrans(d.kind,kinds=gkinds) ) ) # hdrfromfb(df, di._variables, cdmfb_noodb[d.element_name], ttrans(d.kind,kinds=gkinds)) # edit (di['recordindex'].shape[0], hdrfromfb(df, di._variables, cdmfb_noodb[d.element_name], ttrans(d.kind,kinds=gkinds) ) )
                        except:
                            if d.element_name in ['latitude','longitude']:
                                groups[k][d.element_name]=hdrfromfb(df, di._variables, cdmfb_noodb[d.element_name], ttrans(d.kind,kinds=gkinds))
                            else:
                                x=numpy.zeros(di['recordindex'].shape[0], dtype=numpy.dtype(ttrans(d.kind,kinds=gkinds)))
                                if numpy.dtype(ttrans(d.kind,kinds=gkinds)) == np.int32:
                                    x.fill(int_void)
                                else:
                                    x.fill(numpy.nan)
                                groups[k][d.element_name] = ({'hdrlen':di['recordindex'].shape[0]},x) 

            elif k in ('station_configuration'): # station_configurationt contains info of all the stations, so this extracts only the one line for the wanted station with the numpy.where
                try: # case when the station+conf cannot be retrieved 
                    if d.element_name in station_configuration_retrieved.columns:                            
                        groups[k][d.element_name]=({'hdrlen': 1}, np.full( 1 , station_configuration_retrieved[d.element_name].values[0] ) )
                except:
                    pass

            elif k in ('source_configuration'): # storing the source configuration info, e.g. original file name, 
                if d.element_name=='source_file':
                    
                    #  groups[k][d.element_name] = ( {'hdrlen':fbds.variables['date@hdr'].shape[0] } ,  np.full( fbds.variables['date@hdr'].shape[0] , source_file  ) ) 
                    groups[k][d.element_name] = ( {'hdrlen': 1 },   np.full( 1 , source_file) )
                else:
                    try:   
                        groups[k][d.element_name] = (  {'hdrlen': 1 }, np.full( 1 , np.nan) ) # element_name is the netcdf variable name, which is the column name of the cdm table k 
                    except KeyError:
                        pass

            else : # this is the case where the cdm tables DO exist
                try:   
                    groups[k][d.element_name]=({k+'_len':len(cdm[k])}, cdm[k][d.element_name].values) # element_name is the netcdf variable name, which is the column name of the cdm table k 
                except KeyError:
                    pass

            """ Tryin to add attributes, e.g. description and external tables """     
            try:
                groups[k][d.element_name].attrs['external_table']=d.external_table # defining variable attributes that point to other tables (3rd and 4th columns)
                groups[k][d.element_name].attrs['description']=d.description  # it faisl when trying with the observations_table 
            except:
                pass
            try:  
                if k == 'sensor_configuration':
                    print('')
                # TO DO TODO the new writ_dic does not work, I dont have time to debug it now and I dont need it right now, do it later 
                OLD = True
                if OLD:
                    # OLD nonsense
                    if type(groups[k][d.element_name].values[0])==str:
                        s=groups[k][d.element_name].values.shape
                        groupencodings[k][d.element_name]={'dtype':numpy.dtype('S80'),'compression': 'gzip','chunksizes':(min(100000,s[0]),80)}  ### OLD VERSION WORKING with the old write_dict   
                    else:
                        groupencodings[k][d.element_name]={'compression': 'gzip'}

                    if k in ('observations_table'):     # edit ???      
                        if d.element_name =='source_id':
                            a=0
                        write_dict_h5_new(fno, groups[k], k, groupencodings[k], var_selection=[d.element_name],mode='a', attrs= dic_obstab_attributes ) # old
                else:
                    ### new version of wrtoe_dict, does not work 
                    if d.element_name != 'index':
                        groupencodings[k][d.element_name]= { 'compression': 32015  ,'compression_opts': 3 } 
                    else:
                        groupencodings[k]['index']= { 'compression': 'gzip' ,'compression_opts': 4 } 

                    # {'observations_id': { 'compression': 'gzip' } ,'compression_opts': 4 }} or
                    #      {'observations_id': { 'compression': 32015 } ,'compression_opts': 3 }} 

                    if k in ('observations_table'):                
                        write_dict_h5(fno, groups[k], k, groupencodings[k], mode='a', attrs= dic_obstab_attributes ) 

            except:
                print('bad:',k,d.element_name)
                pass

    """
                        if type(groups[k]) is dict:
                        gkev=groups[k][d.element_name]
                    else:
                        gkev=groups[k][d.element_name].values
                    if type(gkev[0])==str:
                        s=gkev.shape
                        groupencodings[k][d.element_name]={'dtype':numpy.dtype('S80'),'compression': 'gzip','chunksizes':(min(100000, s[0] ) , 80 ) }
                    else: 
                        groupencodings[k][d.element_name]={'compression': 'gzip'}

                    if k in ('observations_table'):
                        #print(k,d.element_name,time.time()-tt,' mem:',process.memory_info().rss//1024//1024)
                        write_dict_h5(fno, groups[k], k, groupencodings[k], var_selection=[],mode='a', attrs= dic_obstab_attributes  )
                except:
                    #print('bad:',k,d.element_name)
                    pass
                #if k=='observations_table':
                #    print(k,d.element_name,time.time()-tt,' mem:',process.memory_info().rss//1024//1024)

    for k in groups.keys():            
            ##this appends group by group to the netcdf file
            if k not in ['observations_table'] :   
                try:
                    groups[k].to_netcdf(fno,format='netCDF4',engine='h5netcdf',encoding=groupencodings[k],group=k,mode='a') #
                except:
                    pass

    """
    for k in groups.keys():            
        if k not in ('observations_table') :   
            if k in ['id_scheme', 'observed_variable', 'crs']:
                continue      
            #write_dict_h5(fno, groups[k], k, groupencodings[k], var_selection=[],mode='a', attrs= dic_obstab_attributes )
            groups[k].to_netcdf(fno,format='netCDF4',engine='h5netcdf',encoding=groupencodings[k],group=k,mode='a') 
            #old groups[k].to_netcdf(fno,format='netCDF4',engine='h5netcdf',encoding=groupencodings[k],group=k,mode='a') 

    del df

    #print(0)
    log_name = fno.replace('.nc', '_correctly_processed_year.txt') # fno.replace('_'+str(year),'').replace('.nc', '_correctly_processed_year.txt')
    a = open(log_name, 'a+')
    a.write(year+ '\n')
    a.close()

    print(' Finished writing ' , fno )
    return log_name


def get_station_configuration_cuon(stations_id='', station_configuration='', lat='', lon='', fn = '', db='', change_lat=False):

    """ Gets the primary station_id from the station_configuration table. 
         station_id is the id taken from the input file.
         First it checks if a primary_id in th estation_conf file matches the station_id, 
         otherwise it looks for an element in the list of secondary ids.      
         New version using CUON combined station_configuration files
         """

    f = fn.split("/")[-1]
    if db in  ["era5_1", 'era5_1_mobile']:
        fname = str.encode(  "era5.conv._" + stations_id[0] )
    else:
        if 'era5_1' in db or 'era5_3188' in db :
            fname = str.encode( fn.split("/")[-1].replace(".gz","").replace("_","") )
        elif "era5_2" in db:
            fname = str.encode( fn.split("/")[-1].replace(".gz","") ) 
        elif db == "ncar":
            fname = str.encode(fn.split("/")[-1] )
        elif db == "igra2":
            fname = str.encode(stations_id[0] + "-data.txt")
        elif db == "igra2_mobile":
            fname = str.encode(stations_id[0] + "-data.txt")
        elif db in ["bufr", 'amma', 'giub', 'hara', 'maestro']:
            fname = str.encode(fn.split('/')[-1] )
        elif db in ['bufr_cnr']:
            fname = str.encode(fn.split('/')[-1])
        elif db in ['woudc']:
            fname = str.encode(fn.split('/')[-1])
        elif db in ['shipsound', 'npsound']:
            fname = str.encode(fn.split('/')[-1]) # .split('.')[0].replace('_', '')
        elif 'mauritius' in db or 'yangjiang' in db :
            return station_configuration

    d = station_configuration.loc[station_configuration['file'] == fname ]
    if d.empty:
        # d = station_configuration.loc[[stations_id[0].encode() in s for s in station_configuration['file']]]
        d = station_configuration.loc[[stations_id[0] in s.decode().split('.')[-1] for s in station_configuration['file']]]
    if d.empty:
        d = station_configuration.loc[[stations_id[0] in s.decode().split('.')[-2] for s in station_configuration['file']]]
    #print(0) # there must always be a matching stat conf since we are checking the file name now
    if  d.empty:
        a = open( 'logs/' + db + "_wrong_stat_conf.dat" , "a+")
        a.write(fn+'\n')
        #print('must check wrong coordinates' )
        return None
    return d





def check_lat_lon(fbds, fn, save= 'correct'):
    """ Check if the data frame contains consisten latitude and longitude.
          If more than 99% of the data has inconsistent values for lat/lon,
          the entire file will be flagged wih '-20999-' code. 

          If less than 1%, the majority of the df with conisstent lat/lon will be kept,
          while the rest will be discarded and flagged with  '-20999-' 

          NB should not happen in cop2, since incosistent files were already flagged and removed from the 
          new station_configuration files.

          NB
          The check for consistency is made based on valued of distance in km [30km] and not coordinates in the inventory,
          so the two methods might give slightly different results.
          """

    import operator

    df_lat_lon = fbds[['lat@hdr', 'lon@hdr']].astype(float)
    df_lat_lon = df_lat_lon.drop_duplicates()

    fbds['lat@hdr'] = fbds['lat@hdr'] .astype(float)
    fbds['lon@hdr'] = fbds['lon@hdr'].astype(float)

    lats, lats_indices, lats_counts = np.unique(fbds['lat@hdr']  , return_index=True, return_counts=True )
    lons, lons_indices, lons_counts  = np.unique(fbds['lon@hdr'] , return_index=True, return_counts=True )

    lat_to_counts_dic = dict(zip(lats, lats_counts))
    lon_to_counts_dic = dict(zip(lons, lons_counts))

    most_freq_lat = max(lat_to_counts_dic.items(), key=operator.itemgetter(1))[0]
    most_freq_lon = max(lon_to_counts_dic.items(), key=operator.itemgetter(1))[0]

    len_data = len(fbds)

    # find the rows with lat and lon compatible with the most frequent value
    good_lats = np.where ( abs(fbds['lat@hdr'] - most_freq_lat) < 0.5 ) [0]
    good_lons = np.where ( abs(fbds['lon@hdr']  - most_freq_lon) < 0.5 ) [0]

    # if all good data, return the df as it is 
    if len(good_lats) == len_data and len_data == len(good_lons):
        return True, fbds, most_freq_lat, most_freq_lon


    # indices of the majority of the 
    combined_indices = np.intersect1d(good_lats, good_lons)

    # case where I have > 99% consistent data
    # will save either the good data or the bad one
    tot_good_data = len(combined_indices)/len_data


    if tot_good_data  >= 0.99:
        if save=='correct':
            good_df = fbds.iloc[combined_indices]
            a = open('log_lat-lon_checks.log', 'a')
            a.write(fn+'_good_'+ str(tot_good_data) + '_bad_' + str(1-tot_good_data) + '\n' )

            return True, good_df, most_freq_lat, most_freq_lon

        else:
            all_ind = list(range(len_data))
            wrong_ind = [ f for f in all_ind if f not in combined_indices ]
            bad_df = fbds.loc[wrong_ind]
            return False, bad_df, most_freq_lat, most_freq_lon

    # will save the bad data
    else:
        return False, fbds, most_freq_lat, most_freq_lon

    #primary_id = primary_id.replace('-1_', '0-20999-').replace('-20000-','-20999-').replace('-20300-','-20999-')  
    #primary_id = primary_id.replace('-20400-', '0-20999-').replace('-20500-','-20999-').replace('-20600-','-20999-')  

    return 0




def write_odb_to_cdm(fbds, cdm, cdmd, output_dir,  dataset, dic_obstab_attributes, fn, fns, change_lat, change_lon, year):
    """ Write the data to file 
    If coming from era5_1/ era5_2 + mobile, no need for year selection, otherwise must split """
    tt=time.time()

    """ Read the station_id, getting the station_configuration from the table list, extracting primary_id """      
    station_id =  [ fbds['statid@hdr'].values[0][1:-1].decode('utf-8') ]

    # TO DO verify it still works with igra, ncar etc. 
    #station_configuration_retrieved = stations_id='', station_configuration='', lat='', lon='', fn = '', db='', change_lat=False          
    station_configuration_retrieved = get_station_configuration_cuon(stations_id=station_id, 
                                                                    station_configuration = cdm['station_configuration'],
                                                                    lat=fbds['lat@hdr'].values[0],  lon=fbds['lon@hdr'].values[0], 
                                                                    fn=fn, 
                                                                    db=dataset,
                                                                    change_lat=change_lat  )       

    primary_id = station_configuration_retrieved.primary_id.values[0].decode('utf-8')

    fno,  source_file = initialize_output(fn, output_dir, primary_id, dataset, year)         

    log_name = fno.replace('.nc', '_correctly_processed_year.txt') # fno.replace('_'+str(year),'').replace('.nc', '_correctly_processed_year.txt')
    if dataset not in ['era5_1' , 'era5_1_mobile']:
        fbds = fbds.loc[fbds['year'].astype(int) == year ]
        if fbds.empty:
            return log_name



    # checking for missing minus sign from era5 1759
    if change_lat and 'era5_1759' in dataset:
        print('Changing latitude - WBAN missing sign ')
        fbds['lat@hdr'] = - fbds['lat@hdr']
    elif change_lat != False:
        print('Changing latitude - WBAN missing sign ')
        fbds['lat@hdr'][:] = change_lat
    elif change_lon != False:
        print('Changing longitude - WBAN missing sign ')
        fbds['lon@hdr'][:] = change_lon


    # check consistent lat and long throughout the file
    # save = ['correct', 'wrong']

    if 'mobile' not in dataset:
        correct_data, fbds, most_freq_lat, most_freq_lon = check_lat_lon(fbds, fn, save='correct')

        fbds = fbds.reset_index()
        fbds = fbds.drop(columns = ["index", "level_0"])

        # check if retireved station inventory lat and lon are compatible with file 
        stat_conf_check = True

        ### TODO TO DO if this really needed ???
        if isinstance(station_configuration_retrieved, pd.DataFrame):
            if '20666' in primary_id: # these are stations with known lat/lon problems or inconsistencies 
                stat_conf_check = True 
            else:
                if not station_configuration_retrieved.empty:
                    sc_lat, sc_lon = station_configuration_retrieved.latitude.values[0] , station_configuration_retrieved.longitude.values[0]
                    try:
                        
                        if type(sc_lat) == bytes:
                                sc_lat = float(sc_lat.decode('utf-8'))
                                sc_lon = float(sc_lon.decode('utf-8'))
                    except:
                        print("CHECK MISTAKE IN DECODING LAT LON" , fn, '  -  ' , fns  )
                        
                        
                    try:
                        if abs(sc_lat - most_freq_lat ) < 1.0 and  abs(sc_lon - most_freq_lon ) < 1.0:
                            stat_conf_check = True
                    except:
                        print("CHECK MISTAKE IN DECODING LAT LON" , fn, '  -  ' , fns  )
                        
                    else:
                        stat_conf_check = False

        if dataset == "era5_1759":
            try: # the stat_cof_retr might be a None, so it wont work 
                if (station_configuration_retrieved['latitude'].values[0] < 0 and fbds['lat@hdr'][0] > 0) :
                    print('*** Correcting latitude missing minus sign from WBAN archive ') # files from era5_1759 might have missing minus sign from the WBAN inventory 
                    fbds['lat@hdr'] =- fbds['lat@hdr']
            except:
                pass

        if dataset == "era5_2":
            try: # the stat_cof_retr might be a None, so it wont work 
                if (station_configuration_retrieved['latitude'].values[0] < 0 and fbds['lat@hdr'][0] > 0) :
                    print('*** Correcting latitude missing minus sign from WBAN archive ') # files from era5_2 might have missing minus sign from the WBAN inventory 
                    fbds['lat@hdr'] =- fbds['lat@hdr']
            except:
                pass

        '''
        try:
            primary_id = station_configuration_retrieved['primary_id'].values[0].decode('utf-8')     
        except:
            primary_id = '0-209990-0-' + str(station_id[0]) + '-coordinateIssue-'
            if 'mobile' in dataset:
                primary_id = primary_id + 'mobileStation-'

            primary_id = primary_id.replace('-1_', '0-20999-').replace('-20000-','-20999-').replace('-20300-','-20999-') .replace("-20001-","-20999-")
            primary_id = primary_id.replace('-20400-', '0-20999-').replace('-20500-','-20999-').replace('-20600-','-20999-')  
        if not stat_conf_check and '999' not in primary_id:
            primary_id = 'stat_conf_inconsistent_' + primary_id
        '''



    fno='.'.join(fno.split('.')[:2]+fno.split('.')[2:])

    if fbds is None:
        return

    y=numpy.int64(fbds['date@hdr'].values)*1000000+fbds['time@hdr'].values
    tt=time.time()
    idx=numpy.lexsort((fbds['vertco_reference_1@body'].values,y))  # sorting by z_coordinate and dat_time 
    y=y[idx]
    #print(time.time()-tt)
    for fb in fbds.keys():
        fbds[fb].values[:]=fbds[fb].values[idx]
    #print(time.time()-tt)
    x=numpy.unique(y)
    #print(time.time()-tt)
    z=find_recordindex_l(y,x)
    #print(time.time()-tt)
    di=xr.Dataset() 
    di['recordindex']=({'record':z.shape[1]},z[1])
    #x=make_datetime(x//1000000,x%1000000)
    df=pd.DataFrame({'year':x//10000000000,'month':(x%10000000000)//100000000,'day':x%100000000//1000000,
                     'hour':x%1000000//10000,'minute':(x%10000)//100,'second':x%100})
    dt=pd.to_datetime(df).values
    di['recordtimestamp']=({'record':z.shape[1]}, numpy.array(dt-numpy.datetime64('1900-01-01'),dtype=int)//1000000000)
    di['recordtimestamp'].attrs['units']='seconds since 1900-01-01 00:00:00'
    del dt,df
    y=fbds['date@hdr'].values
    x=numpy.unique(y)
    z=find_dateindex_l(y,x)
    di['dateindex']=({'days':z.shape[1],'drange':z.shape[0]},z) # date, index of the first occurrance, index of the last
    del y

    di.to_netcdf(fno,format='netCDF4',engine='h5netcdf',mode='w')

    ### storing all econdings for variables in both the era5fb and other CDM tables 
    encodings_variables = {} 

    # setting encodings 
    fbencodings={}
    for d,v in fbds.items():
        if v.dtype==numpy.dtype('float64'):
            if d!='date@hdr':             
                fbencodings[d]={'dtype':numpy.dtype('float32'),'compression': 'gzip'} # probably dtype not neccessary, but compression must be there
            else:
                fbencodings[d]={'dtype':numpy.dtype('int32'),'compression': 'gzip'}               
        else:
            if type(v.values[0])==bytes:
                fbencodings[d]={'compression': 'gzip','chunksizes':(min([10000,v.shape[0]]),10)}#,'chunksizes':(10000,10)
            else:
                fbencodings[d]={'compression': 'gzip'}
    fbencodings['index']={'compression': 'gzip'}  

    fbds = fbds.drop(columns=['year'])
    
    #fbds = fbds[:2000]  # TO DO HERE TO DO CHANGE!!!
    write_dict_h5_new(fno, fbds, 'era5fb', fbencodings, var_selection=[],mode='a') #old


    dcols=[]
    for d in fbds.columns:
        if d not in ['date@hdr','time@hdr','statid@hdr','vertco_reference_1@body','varno@body', 'lon@hdr','lat@hdr','seqno@hdr',
                     'obsvalue@body','fg_depar@body','an_depar@body','biascorr@body','sonde_type@conv', 'vertco_type@body' , 'source_id']:
            dcols.append(d)
    fbds.drop(columns=dcols,inplace=True)

    #if dataset not in ['era5_1' , 'era5_1_mobile' , 'era5_2', 'era5_2_mobile']:
    #    fbds['year'] = 0



    groups={}
    groupencodings={}
    sc = {}         
    
    if os.path.isfile( 'encodings.txt'  ):
        lines =  open('encodings.txt' , 'r').readlines()
        all_vars = [ l.split('\t')[0] for l in lines]
        all_vars = list (np.unique(all_vars) )
    else:
        all_vars = []
        
    for k in cdmd.keys(): # loop over all the table definitions 
        if k in ('observations_table'):
            pass 
        else:
            groups[k]=xr.Dataset() # create an  xarray
        groupencodings[k]={} # create a dict of group econdingsep=..., end=..., file=..., flush=...

        for i in range(len(cdmd[k])): # in the cdm table definitions you always have the element(column) name, the type, the external table and the description 
            d=cdmd[k].iloc[i] # so here loop over all the rows of the table definition . iloc is just the index of the element
            # two possibilities: 1. the corresponding  table already exists in the cdm (case in the final else)
            #                    2. the corr table needs to be created from the local data sources (e.g. the feedback or IGRA or smt else). 
            # These are the observation_tables, the header_tables and the station_configuration.
            # These tables are contained in the CEUAS GitHub but not in the cdm GitHub
            if k in ('observations_table'):
                groups[k]=dict() #pd.DataFrame()
                try:
                    # fbds is an xarray dataset , fbds._variables is a dict of the variables 
                    if d.element_name=='report_id':                            
                        groups[k][d.element_name]=fromfb_l(fbds,di._variables,cdmfb[k+'.'+d.element_name], ttrans(d.kind,kinds=okinds))

                    else:
                        groups[k][d.element_name]=fromfb_l(fbds,di._variables,cdmfb[d.element_name], ttrans(d.kind,kinds=okinds))
                        if d.element_name=='date_time':
                            fbds["date_time"] = groups['observations_table']['date_time']

                except KeyError:
                    x=numpy.zeros(  fbds['date@hdr'].shape[0],dtype=numpy.dtype(ttrans(d.kind,kinds=okinds) ) )
                    try:
                        x.fill(numpy.nan)
                    except:
                        x.fill(-2147483648)
                    groups[k][d.element_name]=x


            elif k in ('header_table'):
                try:
                    if d.element_name=='report_id':
                        groups[k][d.element_name]=( {'hdrlen':di['recordindex'].shape[0]}, hdrfromfb(fbds,di._variables, cdmfb[k+'.'+d.element_name],ttrans(d.kind,kinds=gkinds) ) )
                    else:
                        groups[k][d.element_name]=( {'hdrlen':di['recordindex'].shape[0]}, hdrfromfb(fbds,di._variables, cdmfb[d.element_name],ttrans(d.kind,kinds=gkinds) ) )
                    j=0

                except KeyError:   
                    pass
                    #print ('FFF ', d.element_name , ' ' , numpy.dtype(ttrans(d.kind,kinds=gkinds)), time.time()-tt )
                    #print ('FFF ', d.element_name , ' ' , numpy.dtype(ttrans(d.kind,kinds=gkinds)) )

                    if d.element_name in cdm['station_configuration'].columns:
                        x=numpy.zeros(di['recordindex'].shape[0], dtype=numpy.dtype(ttrans(d.kind,kinds=gkinds)))
                        try:

                            idx=numpy.where('0-20000-0-'+fnl[-1].split('_')[-1] == cdm['station_configuration']['primary_id'])[0][0]
                            groups[k][d.element_name]=x.fill(cdm['station_configuration'][d.element_name][idx])
                        except:
                            groups[k][d.element_name]= ({'hdrlen':di['recordindex'].shape[0]},x)
                    else:
                        x=numpy.zeros(di['recordindex'].shape[0], dtype=numpy.dtype(ttrans(d.kind,kinds=okinds)))
                        try:    
                            x.fill(numpy.nan)
                        except:
                            x.fill(-2147483648)
                    groups[k][d.element_name]=({'hdrlen':di['recordindex'].shape[0]},x)

            elif k in ('station_configuration'): # station_configurationt contains info of all the stations, so this extracts only the one line for the wanted station with the numpy.where
                try:                        
                    groups[k][d.element_name]=({'hdrlen': 1}, np.full( 1 , station_configuration_retrieved[d.element_name].values[0] ) )

                except:
                    #print("Failing station_conf" , k )
                    pass

            elif k in ('source_configuration'): # storing the source configuration info, e.g. original file name, 
                if d.element_name=='source_file':
                    #  groups[k][d.element_name] = ( {'hdrlen':fbds.variables['date@hdr'].shape[0] } ,  np.full( fbds.variables['date@hdr'].shape[0] , source_file  ) ) 
                    groups[k][d.element_name]=({'hdrlen': 1 },   np.full( 1 , source_file) )
                else:
                    try:   
                        groups[k][d.element_name]=({'hdrlen': 1 },   np.full (1, hdrfromfb(fbds,di._variables, cdmfb[k+'.'+d.element_name],ttrans(d.kind,kinds=gkinds) )  ) ) # element_name is the netcdf variable name, which is the column name of the cdm table k 
                    except KeyError:
                        pass

            else : # this is the case where the cdm tables DO exist
                try:   
                    groups[k][d.element_name]=({k+'_len':len(cdm[k])}, cdm[k][d.element_name].values) # element_name is the netcdf variable name, which is the column name of the cdm table k 
                    print('OK group' , k )
                    
                except KeyError:
                    print("FAILING " , k )
                    pass

            """ Tryin to add attributes, e.g. description and external tables """ 
            if k not in ('observations_table'):               
                try:
                    groups[k][d.element_name].attrs['external_table']= d.external_table # defining variable attributes that point to other tables (3rd and 4th columns)
                    groups[k][d.element_name].attrs['description']     = d.description  # it fails when trying with the observations_table 

                    sc[d.element_name] = {}
                    sc[d.element_name]['description'] = d.description
                    sc[d.element_name]['type'] = type(groups[k][d.element_name].values[0] )
                    sc[d.element_name]['external_table'] = d.external_table

                except:
                    pass

            try:
                if type(groups[k]) is dict:
                    gkev=groups[k][d.element_name]
                else:
                    gkev=groups[k][d.element_name].values
                    
                gkev = np.array(gkev)
                if k in ('observations_table') and d.element_name == 'source_id':
                    a= 0
                    
                if k in ('source_configuration') and d.element_name == 'source_id':
                        a= 0
                        
                if k in ('observations_table') and d.element_name == 'processing_code':
                    gkev = gkev.astype('float')
                    groups[k][d.element_name] = gkev 
                     
                if type(gkev[0])==str: # TO DO HERE TODO
                    s=gkev.shape
                    groupencodings[k][d.element_name]={'dtype':numpy.dtype('S80'),'compression': 'gzip','chunksizes':(min(100000, s[0] ) , 80 ) }
                else: 
                    groupencodings[k][d.element_name]={'compression': 'gzip'}

                if k in ('observations_table'):
                    #if d.element_name == 'z_coordinate_type':
                    #    a=0
                    write_dict_h5_new(fno, groups[k], k, groupencodings[k], var_selection=[],mode='a', attrs= dic_obstab_attributes  ) # old
            except:
                #print('bad:',k, d.element_name)
                pass
                #if k=='observations_table':
                #    print(k,d.element_name,time.time()-tt,' mem:',process.memory_info().rss//1024//1024)

                ### updating the data type encoding 
            if d.element_name not in all_vars and k != 'observations_table':
                #if k =='header_table':
                #    a=0
                try:
                    out_en = open('encodings.txt' , 'a+')
                    out_en.write(d.element_name + '\t' + k + '\t' +  str(  type(groups[k][d.element_name].values[0]) )  + '\n')
                    out_en.close()
                    all_vars.append(d.element_name)   
                except:
                    pass
                        
    for k in groups.keys():            
        ##this appends group by group to the netcdf file
        if k not in ['observations_table'] :   
            try:
                groups[k].to_netcdf(fno,format='netCDF4',engine='h5netcdf',encoding=groupencodings[k],group=k,mode='a') #
            except:
                pass
    ##print('sizes: in: {:6.2f} out: {:6.2f}'.format(os.path.getsize(fn+'.gz')/1024/1024, os.path.getsize(fno)/1024/1024))
    del fbds

    """ Storing the group encodings in a numpy dictionary to be reused by the merging script """
    np.save('groups_encodings',  groupencodings)
    np.save('station_configuration_encodings',  sc)
    np.save('era5fb_encodings',  fbencodings)

    #print(0)
    log_name = fno.replace('.nc', '_correctly_processed_year.txt') #  fno.replace('_'+str(year),'').replace('.nc', '_correctly_processed_year.txt')
    a = open(log_name, 'a+')
    a.write(str(year)+ '\n')
    a.close()

    return log_name




def read_odb_to_cdm(output_dir, dataset, dic_obstab_attributes, fn, fns):
    """ Convert the  file into cdm compliant netCDF files. 
        According to each file type, it will use the appropriate reading function to extract a Pandas DataFrame
        input:
              fn       :: odb file name (e.g. era5.conv._10393)
              cdm   :: cdm tables (read with pandas)
              cdmd :: cdm tables definitions ("")  """
    """ Index(['type', 'class', 'stream', 'andate', 'antime', 'reportype',
       'numtsl@desc', 'timeslot@timeslot_index', 'seqno@hdr', 'bufrtype@hdr',
       'subtype@hdr', 'groupid@hdr', 'obstype@hdr', 'codetype@hdr',
       'sensor@hdr', 'date@hdr', 'time@hdr', 'report_status@hdr',
       'report_event1@hdr', 'report_rdbflag@hdr', 'lat@hdr', 'lon@hdr',
       'lsm@modsurf', 'orography@modsurf', 'windspeed10m@modsurf',
       'tsfc@modsurf', 'albedo@modsurf', 'seaice@modsurf',
       'snow_depth@modsurf', 'sonde_type@conv', 'station_type@conv',
       'collection_identifier@conv', 'timeseries_index@conv',
       'unique_identifier@conv', 'entryno@body', 'obsvalue@body', 'varno@body',
       'vertco_type@body', 'vertco_reference_1@body',
       'vertco_reference_2@body', 'ppcode@conv_body', 'datum_anflag@body',
       'datum_status@body', 'datum_event1@body', 'datum_rdbflag@body',
       'biascorr@body', 'biascorr_fg@body', 'varbc_ix@body', 'qc_pge@body',
       'an_depar@body', 'an_sens_obs@body', 'fg_depar@body',
       'an_depar@surfbody_feedback', 'fg_depar@surfbody_feedback',
       'snow_depth@surfbody_feedback', 'snow_density@surfbody_feedback',
       'datum_status@surfbody_feedback', 'datum_sfc_event@surfbody_feedback',
       'lsm@surfbody_feedback', 'stalt@hdr', 'obs_error@errstat',
       'final_obs_error@errstat', 'fg_error@errstat', 'eda_spread@errstat',
       'expver', 'source@hdr', 'statid@hdr', 'source_id'], """

    process = psutil.Process(os.getpid())
    t=time.time()

    ### only era5_1 and era5_2 (mobile) have year splitting 
    
    # if 'mobile' in dataset:
    #     fbds = []
    #     for fnsi in fns:
    #         fbds.append(read_all_odbsql_stn_withfeedback(dataset, fnsi))
    # else:

    p=Pool(3)
    #func=partial(read_all_odbsql_stn_withfeedback,dataset)
    ### uncomment to debug one single file 
    #for f in fns:
    #    d = read_all_odbsql_stn_withfeedback(dataset, f )

    #fns =[f for f in fns if '202212' in f ]

    if len(fns)==0:
        #a = open('missing_gzipped_files_era5_2.txt' , 'a')
        #a.write(fn + '\n' )
        #a.close()
        #print('COULD NOT FIND file ++++ ' , fn )
        return 

    func=partial(read_all_odbsql_stn_withfeedback,dataset)

    print('Reading ' + str(len(fns)) + ' ODB files ')

    if len(fns)==1:       
        fbds= read_all_odbsql_stn_withfeedback(dataset,fns[0]) # list(map(func,fns))
        p.close()
        del p
    else:
        fbds=list(p.map(func,fns))
        p.close()
        del p

        try:
            fbds=pd.concat(fbds,axis=0,ignore_index=True)
        except:
            return None
    fbds = fbds.reset_index()
    fbds['vertco_reference_1@body'] = fbds['vertco_reference_1@body'].astype(float)
    fbds = fbds.sort_values(by = ['date@hdr', 'time@hdr', 'vertco_reference_1@body' ] )    

    years = [str(s)[:4] for s in fbds['date@hdr'] ]
    fbds['year'] = years
    #fbds = fbds.replace( -2147483648 , np.nan ) 

    # dropping NANS
    fbds = fbds.dropna(subset=['obsvalue@body'] )   
    if len(fbds) == 0:
        print('#############################')
        print('NO DATA IN: ', 'obsvalue@body')
        print('#############################')
        return None


    # removing duplicated observations within one single record (same z coordinate, variable number)
    fbds = fbds.drop_duplicates(subset=['date@hdr', 'time@hdr', 'vertco_reference_1@body', 'varno@body'] )

    # sorting values 
    fbds = fbds.sort_values(by = ['date@hdr', 'time@hdr', 'vertco_reference_1@body' ] )   

    print('+++ Created pandas dataframe ')
    
    """
    #fbds = fbds[:10]  # TO DO HERE TO DO CHANGE WRONG  !!!!!!!!!
    fbds = fbds[:20000]    
    fbds['date_time'] = fbds['date@hdr'].astype(str) + fbds['time@hdr'].astype(str) 
    print(np.unique(fbds['date_time']))
    
    ## I obtain:
    ### '20160120163500' '20160120170000' '20160120203400' 
    """

 
    return fbds, min(years)



def process_nasa(fdict,cdmd,cdm,fn,data):
    if(not fdict['observations_table']):

        for k in cdmd['observations_table'].element_name.values:
            fdict['observations_table'][k]=[]
        for k in cdmd['header_table'].element_name.values:
            fdict['header_table'][k]=[]
    refdate=datetime(year=1900,month=1,day=1)
    header=True
    aux={}
    l=0
    if 'phase4' in fn:
        l-=1
        with open('/'.join(fn.split('/')[:-1])+'/readme') as g:
            gdata=g.read().split('\r\n')
        idx=np.where(cdm['station_configuration'].primary_id==key)[0][0]
        for g in gdata:
            if 'WMO' in d:
                if 'measurement_programme' not in fdict['header_table'].keys():
                    fdict['header_table']['measurement_programme']=[]
                fdict['header_table']['measurement_programme'].append(d.strip())
            if 'nbo' in g:
                glist=g.split()
                relset=datetime.strptime(glist[2]+' '+glist[3],'%d/%b/%y %H:%M:%S')
                fdict['header_table']['report_timestamp'].append(int((relset-refdate).total_seconds()))
                fdict['header_table']['record_timestamp'].append(int((relset-refdate).total_seconds()))
                fdict['header_table']['primary_station_id'].append(cdm['station_configuration'].primary_id.values[idx])
                fdict['header_table']['report_id'].append(fdict['header_table']['primary_station_id'][-1]+b'-'+key+b'-'+glist[1].split('-')[-1].strip())
    else:
        for d in data:
            if b'WMO' in d:
                if 'measurement_programme' not in fdict['header_table'].keys():
                    fdict['header_table']['measurement_programme']=[]
                fdict['header_table']['measurement_programme'].append(d.strip())
            if b'ASCENSION NUMBER' in d:
                dlist=d.split(b'=')
                key=np.bytes_('0-20100-0-0000'+fn[-8])
                idx=np.where(cdm['station_configuration'].primary_id==key)[0][0]
                fdict['header_table']['primary_station_id'].append(cdm['station_configuration'].primary_id.values[idx])
                fdict['header_table']['report_id'].append(fdict['header_table']['primary_station_id'][-1]+b'-'+key+b'-'+dlist[1].split()[0].strip())
                relset=datetime.strptime(dlist[2].split()[0].decode()+' '+dlist[3].split()[0].decode(),'%d-%b-%y %H:%M:%S')
                fdict['header_table']['report_timestamp'].append(int((relset-refdate).total_seconds()))
                fdict['header_table']['record_timestamp'].append(int((relset-refdate).total_seconds()))
            if b'PARTICIPANT' in d:
                aux['sensor_id']=d.split(b'=')[2].strip()
                try:
                    idx=np.where(cdm['iri_rstype_map'].riname==np.bytes_(fn[-9:-4]))[0][0]
                    aux['sensor_id']=cdm['iri_rstype_map'].vapor_name[idx]
                except:
                    #print('could not find rs type key '+fn[-9:-4]+'!')
                    pass
                aux['report_id']=fdict['header_table']['report_id'][-1]
            if b'MMM' in d:
                header=False
                break
            l+=1

    m=0
    for d in data[l+1:]:
        dlist=d.split()
        if dlist:
            fdict['observations_table']['date_time']+=5*[fdict['header_table']['report_timestamp'][-1]+int(dlist[0])*60]
            fdict['observations_table']['z_coordinate']+=5*[float(dlist[1])]
            fdict['observations_table']['observed_variable']+=[85,138,117,106,107]
            fdict['observations_table']['observation_value']+=[float(i) for i in dlist[2:6]]
            fdict['observations_table']['observation_id']+=['{:0>8}'.format(m+i+1) for i in range(5)]
            m+=5
        else:
            break

    for k,v in aux.items():
        fdict['observations_table'][k]+=[v]*m

    return



def process_ubern(fdict,cdmd,cdm,fn,proflist,campaigns,data):
    """ Analyze data and write info into observations_table and header_table (contained in fdict)
    cdmd is cdm table definitions
    """
    good=False
    for d in data:
        if 'Min' in d:
            good=True
        if 'ri_name' in d:
            key=d.split()[-1]
            if key not in fdict.keys():
                fdict[key]={'observations_table':{},'header_table':{}}

    if not good:
        #print(fn,'could not be processed')
        return
    if(not fdict[key]['observations_table']):

        for k in cdmd['observations_table'].element_name.values:
            fdict[key]['observations_table'][k]=[]
        for k in cdmd['header_table'].element_name.values:
            fdict[key]['header_table'][k]=[]
    refdate=datetime(year=1900,month=1,day=1)
    header=True
    aux={}
    l=0
    comp=fn[:7]
    for d in data:
        if 'Min' in d:
            head=d.split()
            break
        l+=1
    m=0
    for d in data[l+1:]:
        dlist=d.split()
        if dlist:
            if dlist[0]=='NA':
                continue
                #fdict[key]['observations_table']['date_time']+=ddlen*[fdict[key]['header_table']['report_timestamp'][-1]]

            try:
                x=float(dlist[1])*100
            except ValueError:
                try:    
                    if 'PPPcorr' in head:
                        x=float(dlist[2])*100
                    else:
                        continue
                except ValueError:
                    continue
            m+=1

    if m==0:
        #print(fn,'no valid obs')
        return

    try:
        idx=np.where(cdm['iri_rstype_map'].ri_name==key.split('_')[0])[0][0]
        fdict[key]['sensor_id']=np.bytes_(cdm['iri_rstype_map'].vapor_name[idx]) #.split(b',')[0])
    except:
        fdict[key]['sensor_id']=np.bytes_('NA')

    l=0
    for d in data:
        if 'ri_name' in d:
            key=d.split()[1]
            if 'measurement_programme' not in fdict[key]['header_table'].keys():
                fdict[key]['header_table']['measurement_programme']=[]
            idx=campaigns.index[campaigns.ID==comp]
            fdict[key]['header_table']['measurement_programme'].append(campaigns.Name[idx].values[0])        

            id_ascent=fn.split('/')[-1][:11]
            #key=np.bytes_('0-20100-0-0000'+fn[-8])
            sidx=np.where(cdm['station_configuration'].measuring_system_id==np.bytes_(key))[0][0]
            pidx=np.where(proflist['id_ascent']==id_ascent)[0][0]
            fdict[key]['header_table']['primary_station_id'].append(cdm['station_configuration'].primary_id.values[sidx])
            fdict[key]['header_table']['report_id'].append(fdict[key]['header_table']['primary_station_id'][-1]+
                                                           b'-'+np.bytes_(fn.split('_')[-2]))
            fdict[key]['header_table']['latitude'].append(cdm['station_configuration'].latitude.values[sidx])
            fdict[key]['header_table']['longitude'].append(cdm['station_configuration'].longitude.values[sidx])
            fdict[key]['header_table']['height_of_station_above_sea_level'].append(
                proflist['alt'].values[pidx])

            for k in 'latitude','longitude':
                aux[k]=fdict[key]['header_table'][k][-1]

        if 'date ' in d:
            dat=d.split()[-1]
        if 'hour_utc' in d:
            hour=d.split()[-1]
            hh=hour.split(':')
            if len(hh)==2:
                hour=hour+':00'
            relset=datetime.strptime(dat+' '+hour,'%Y-%m-%d %H:%M:%S')
            fdict[key]['header_table']['report_timestamp'].append(int((relset-refdate).total_seconds()))
            fdict[key]['header_table']['record_timestamp'].append(int((relset-refdate).total_seconds()))
        if 'Min' in d:
            head=d.split()
            break
        l+=1

    m=0
    dlen=len(head)-2
    ddlen=dlen
    dstart=2
    if 'Notes' in head:
        dlen-=1
        ddlen-=1
    if 'PPPcorr' in head:
        dstart+=1
        ddlen-=2
    ovar={'TTT':85,'TTTcorr':85,'UU':38,'H':117,'DEWPT':36,'DIR':106,'WSPEED':107}
    offset={'TTT':273.15,'TTTcorr':273.15,'UU':0,'H':0,'DEWPT':273.15,'DIR':0,'WSPEED':0}
    scale={'TTT':1,'TTTcorr':1,'UU':1,'H':9.80655,'DEWPT':1,'DIR':1,'WSPEED':1}
    units={'TTT':5,'TTTcorr':5,'UU':0,'H':806,'DEWPT':5,'DIR':110,'WSPEED':731}
    for d in data[l+1:]:
        dlist=d.split()
        if dlist:
            if dlist[0]=='NA':
                continue
                #fdict[key]['observations_table']['date_time']+=ddlen*[fdict[key]['header_table']['report_timestamp'][-1]]

            try:
                fdict[key]['observations_table']['z_coordinate']+=ddlen*[float(dlist[1])*100]
            except ValueError:
                try:    
                    if 'PPPcorr' in head:
                        fdict[key]['observations_table']['z_coordinate']+=ddlen*[float(dlist[2])*100]
                    else:
                        continue
                except ValueError:
                    continue
                    #fdict[key]['observations_table']['z_coordinate']+=ddlen*[numpy.nan]

            fdict[key]['observations_table']['date_time']+=ddlen*[fdict[key]['header_table']['report_timestamp'][-1]+int(float(dlist[0])*60)]
            if ddlen<len(head[dstart:dlen+2]):
                fdict[key]['observations_table']['observed_variable']+=[ovar[i] for i in head[dstart+1:dlen+2]]
                fdict[key]['observations_table']['units']+=[units[i] for i in head[dstart+1:dlen+2]]
            else:
                fdict[key]['observations_table']['observed_variable']+=[ovar[i] for i in head[dstart:dlen+2]]
                fdict[key]['observations_table']['units']+=[units[i] for i in head[dstart:dlen+2]]
            for i in range(dstart,dlen+2):

                if head[i]!='TTTcorr':                        
                    try:
                        fdict[key]['observations_table']['observation_value'].append(offset[head[i]]+scale[head[i]]*float(dlist[i]))
                    except ValueError:
                        fdict[key]['observations_table']['observation_value'].append(numpy.nan)
                else:
                    try:
                        fdict[key]['observations_table']['observation_value'][-1]=offset[head[i]]+scale[head[i]]*float(dlist[i])
                    except ValueError:
                        pass

            fdict[key]['observations_table']['observation_id']+=[np.bytes_('{:0>8}'.format(m+i+1)) for i in range(ddlen)]
            m+=ddlen
        else:
            break

    aux['sensor_id']='NA'
    try:
        idx=np.where(cdm['iri_rstype_map'].ri_name==key)[0][0]
        aux['sensor_id']=cdm['iri_rstype_map'].vapor_name[idx]
    except:
        #print('could not find rs type key for '+key+'!')
        pass
    aux['report_id']=fdict[key]['header_table']['report_id'][-1]

    for k,v in aux.items():
        fdict[key]['observations_table'][k]+=[v]*m

    return

def process_ubern_aggregated(fdict,cdmd,cdm,fn,proflist,campaigns,data):

    comp=fn[:7]

    refdate=datetime(year=1900,month=1,day=1)
    header=True
    aux={}

    ovar={'temperature':85,'pressure':142,'rel.humidity':38,'relhumidity':38,'geopotential':117,'DEWPT':36,'winddirection':106,'windspeed':107,'radar':117}
    #offset={'TTT':273.15,'TTTcorr':273.15,'UU':0,'H':0,'DEWPT':273.15,'DIR':0,'WSPEED':0}
    scale={'temperature':1,'pressure':100,'rel.humidity':1,'relhumidity':1,'geopotential':9.80655,'DEWPT':1,'winddirection':1,'windspeed':1,'radar':9.80655}
    units={'temperature':5,'pressure':32,'rel.humidity':0,'relhumidity':0,'geopotential':806,'DEWPT':5,'winddirection':110,'windspeed':731,'radar':806}

    cols=list(data.iloc[9])
    for i in range(len(cols)):
        data.rename(columns={data.columns[i]:cols[i]},inplace=True)

    delta=0
    if 'day' in fn:
        delta=12
    elif '1800gmt' in fn:
        delta=18

    mo='mo'
    day='day'
    if proflist is not None:
        relsetf=datetime(proflist['yr'][0],proflist[mo][0],proflist[day][0],delta)
        lpl=len(proflist['yr'])-1
        relsetl=datetime(proflist['yr'][lpl],proflist[mo][lpl],proflist[day][lpl],delta)
    else:
        yrs=data[data.columns[1]][6].split('-')
        if '-' in data[data.columns[1]][6]:
            relsetf=datetime(int(yrs[0]),1,1,delta)
            relsetl=datetime(int(yrs[1]),12,31,delta)
        else:
            relsetf=datetime(int(yrs[0]),1,1,delta)
            relsetl=datetime(int(yrs[0]),12,31,delta)

    m=0
    fdc=[]
    for f in fdict.keys():
        if '_' in f:
            fdc.append(f)
    for i in range(10,len(data)):
        key=data['ri_name'][i]+'_'+data['ri_name_ref'][i]
        try:
            fdi=fdc.index(key)
        except ValueError:
            fdict[key]={'observations_table':{},'header_table':{}}

            for k in cdmd['observations_table'].element_name.values:
                fdict[key]['observations_table'][k]=[]
            for k in cdmd['header_table'].element_name.values:
                fdict[key]['header_table'][k]=[]

            fdc.append(key)
            fdi=len(fdc)-1

            try:
                idx=np.where(cdm['iri_rstype_map'].ri_name==key.split('_')[0])[0][0]
                fdict[key]['sensor_id']=np.bytes_(cdm['iri_rstype_map'].vapor_name[idx]) #.split(b',')[0])
            except:
                fdict[key]['sensor_id']=np.bytes_('NA')
            try:    
                idx=np.where(cdm['iri_rstype_map'].ri_name==key.split('_')[1])[0][0]
                fdict[key]['reference_sensor_id']=np.bytes_(cdm['iri_rstype_map'].vapor_name[idx]) #.split(b',')[0])
            except:
                fdict[key]['reference_sensor_id']=np.bytes_('NA')

        fdict[key]['observations_table']['sensor_id'].append(fdict[key]['sensor_id'])
        fdict[key]['observations_table']['reference_sensor_id'].append(fdict[key]['reference_sensor_id'])

        try:
            fdict[key]['observations_table']['observation_value'].append(numpy.float(data['mean_diff'][i]))
        except ValueError:
            fdict[key]['observations_table']['observation_value'].append(numpy.nan)
        try:
            fdict[key]['observations_table']['original_precision'].append(numpy.float(data['sd'][i]))
        except ValueError:
            fdict[key]['observations_table']['original_precision'].append(numpy.nan)
        try:
            fdict[key]['observations_table']['secondary_value'].append(numpy.float(data['n'][i]))
        except ValueError:
            fdict[key]['observations_table']['secondary_value'].append(numpy.nan)

        fdict[key]['observations_table']['value_significance'].append(102) # mean value obs-ref
        #fdict[key]['observations_table']['sensor_id'].append(np.bytes_(data['ri_name'][i])) # mean value obs-ref
        #fdict[key]['observations_table']['reference_sensor_id'].append(np.bytes_(data['ri_name'][i])) # mean value obs-ref

        if 'pl2' in data.columns:
            fdict[key]['observations_table']['z_coordinate'].append(numpy.float(data['pl2'][i]))
            fdict[key]['observations_table']['reference_z_coordinate'].append(numpy.float(data['pl2'][i]))
        else:
            fdict[key]['observations_table']['z_coordinate'].append(numpy.float(data['pl'][i])*100) #Pa
            fdict[key]['observations_table']['reference_z_coordinate'].append(numpy.float(data['pl'][i])*100) #Pa

        otl=len(fdict[key]['observations_table']['observation_id'])      
        fdict[key]['observations_table']['observation_id'].append(np.bytes_('{:0>8}'.format(otl+1)))

        try:
            fdict[key]['observations_table']['latitude'].append(numpy.float(data[data.columns[1]][5]))
            fdict[key]['observations_table']['longitude'].append(numpy.float(data[data.columns[1]][4]))
        except ValueError: # lat lon range instead of point
            if 'to' in data[data.columns[1]][5]:
                s=data[data.columns[1]][5].split('to')
                fdict[key]['observations_table']['latitude'].append(numpy.mean(numpy.array(s,dtype='float')))
                s=data[data.columns[1]][4].split('to')
                fdict[key]['observations_table']['longitude'].append(numpy.mean(numpy.array(s,dtype='float')))
            else:
                fdict[key]['observations_table']['latitude'].append(numpy.nan)
                fdict[key]['observations_table']['longitude'].append(numpy.nan)

        fdict[key]['observations_table']['date_time'].append((relsetf-refdate).total_seconds())
        fdict[key]['observations_table']['date_time_meaning'].append(1)
        fdict[key]['observations_table']['observation_duration'].append((relsetl-relsetf).total_seconds())
        fdict[key]['observations_table']['observed_variable'].append(ovar[data[data.columns[1]][7]])
        fdict[key]['observations_table']['units'].append(units[data[data.columns[1]][7]])

        id_ascent=key.split('_')[0]
        fdcshort=[]
        for fd in fdc:
            if comp.lower() in fd:
                fdcshort.append(fd)
        fdwi=fdcshort.index(key)
        wigos=np.bytes_('0-20200-0-'+comp[-3:]+'{:0>2}'.format(fdwi+1))

        try:
            fdhi=fdict[key]['header_table']['report_timestamp'].index(int((relsetf-refdate).total_seconds()))
        except:

            fdict[key]['header_table']['report_timestamp'].append(int((relsetf-refdate).total_seconds()))
            fdict[key]['header_table']['record_timestamp'].append(int((relsetf-refdate).total_seconds()))
            htl=np.bytes_(str(len(fdict[key]['header_table']['report_timestamp'])))
            try:

                sidx=np.where(cdm['station_configuration'].measuring_system_id==np.bytes_(id_ascent))[0][0]
                #pidx=np.where(proflist['id_ascent']==id_ascent)[0][0]
                #fdict[key]['header_table']['primary_station_id'].append(cdm['station_configuration'].primary_id.values[sidx])
                fdict[key]['header_table']['primary_station_id'].append(np.bytes_(wigos))
                fdict[key]['header_table']['report_id'].append(fdict[key]['header_table']['primary_station_id'][-1]+b'-'+htl)

                fdict[key]['header_table']['latitude'].append(cdm['station_configuration'].latitude.values[sidx])
                fdict[key]['header_table']['longitude'].append(cdm['station_configuration'].longitude.values[sidx])
            except:
                fdict[key]['header_table']['primary_station_id'].append(np.bytes_(wigos))
                fdict[key]['header_table']['report_id'].append(fdict[key]['header_table']['primary_station_id'][-1]+b'-'+htl)

                fdict[key]['header_table']['latitude'].append(fdict[key]['observations_table']['latitude'][0])
                fdict[key]['header_table']['longitude'].append(fdict[key]['observations_table']['longitude'][0])
            fdhi=-1

        fdict[key]['observations_table']['report_id'].append(fdict[key]['header_table']['report_id'][fdhi])

    #fdict[key]['header_table']['height_of_station_above_sea_level'].append(
        #proflist['alt'].values[pidx])
        m+=1



    return



def ir_to_cdm(cdm, cdmd, output_dir, dataset, dic_obstab_attributes, fn):
    process = psutil.Process(os.getpid())
    t=time.time()
    #fno = initialize_convertion(fn, output_dir) 
    #station_id = ''    
    cdm['iri_rstype_map']=pd.read_csv(os.path.expanduser('../data/tables/iri_rstype_map.dat'),sep='\t',
                                      dtype={'riname':np.dtype('S10'),'vapor_name':'S4','wmo_C2011_code':np.float32,'ri_long_name':'S80'})
    cdm['iri_rstype_map']['wmo_C2011_code']=numpy.int32(cdm['iri_rstype_map']['wmo_C2011_code'])
    cdmd['iri_rstype_map']=pd.DataFrame(data=None, columns=cdmd['observations_table'].columns)
    l=0
    for k in cdm['iri_rstype_map'].columns:
        cdmd['iri_rstype_map']=cdmd['iri_rstype_map'].append({'element_name':k,
                                                              'kind':type(cdm['iri_rstype_map'][k].values[0]),
                                       'description':k},ignore_index=True)
        l+=1

    # to do remove if wrong 
    cdms=pd.read_csv(os.path.expanduser('../data/tables/vapor.instruments.all'),sep=':',names=('sensor_id','comments'))
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

    fdict={'observations_table':{},'header_table':{}}
    flist=glob.glob(fn)
    for fn in flist:

        with zipfile.ZipFile(fn,'r') as a:
            for r in a.filelist:
                if r.filename[0]=='m' or 'min' in r.filename:
                    try:
                        sp=b'\n'
                        if 'min' in r.filename:
                            sp=b'\r\n'
                        data=a.read(r).split(sp)
                        process_nasa(fdict,cdmd,cdm,fn,data)
                    except MemoryError:
                        pass

                #if type(v.values[0])==bytes:
                    #fbencodings[d]={'compression': 'gzip','chunksizes':(min([10000,v.shape[0]]),10)}#,'chunksizes':(10000,10)
                #else:
                    #fbencodings[d]={'compression': 'gzip'}
        #fbencodings['index']={'compression': 'gzip'}

        x=np.array(fdict['observations_table']['date_time'][:])
        idx=[]
        l=0
        iold=b''
        for i in fdict['observations_table']['report_id']:
            if iold!=i:
                idx.append(l)
                iold=i
            l+=1
        #idx.append(len(fdict['observations_table']['report_id']))
        idx=np.array(idx)
        di=xr.Dataset() 
        di['recordindex']=({'record':idx.shape[0]},idx)
        di['recordtimestamp']=({'record':idx.shape[0]}, 
                               numpy.array(fdict['header_table']['record_timestamp']))
        di['recordtimestamp'].attrs['units']='seconds since 1900-01-01 00:00:00'

        l=0
        iold=0
        refdate=datetime(1900,1,1)
        z=[]
        for i in fdict['header_table']['record_timestamp']:
            if i//86400!=iold:
                dat=refdate+timedelta(seconds=i)
                if l==di['recordindex'].shape[0]-1:
                    z.append((dat.year*10000+dat.month*100+dat.day,di['recordindex'].values[l],len(fdict['observations_table']['date_time'])))    
                else:
                    z.append((dat.year*10000+dat.month*100+dat.day,di['recordindex'].values[l],di['recordindex'].values[l+1]))
                iold=i//86400
            l+=1

        z=np.array(z)
        di['dateindex']=({'days':z.shape[1],'drange':z.shape[0]},z) # date, index of the first occurrance, index of the last

        wigos='-'.join(fdict['header_table']['report_id'][0].decode().split('-')[:-2])
        fno='/'.join(fn.split('/')[:-2])+'/nasa/'+wigos+'.nc'
        di.to_netcdf(fno,format='netCDF4',engine='h5netcdf',mode='w')

        print(time.time()-t)
        tt=time.time()

        # each cdm table is written into an hdf group, groups is the dict of all the groups
        # to write the group to the disk, you need the group encoding dict
        groups={}
        groupencodings={}
        for k in cdmd.keys(): # loop over all the table definitions 
            groupencodings[k]={} # create a dict of group econding

            for i in range(len(cdmd[k])): # in the cdm table definitions you always have the element(column) name, the type, the external table and the description 
                d=cdmd[k].iloc[i] 
                if k in ('observations_table'):
                    try:
                        groups[k]={d.element_name:fdict['observations_table'][d.element_name]}
                    except KeyError:
                        x=numpy.zeros(fdict['observations_table']['date_time'].shape[0],
                                      dtype=numpy.dtype(ttrans(d.kind,kinds=okinds) ) )
                        x.fill(numpy.nan)
                        groups[k][d.element_name]=x

                elif k in ('header_table'):
                    try:
                        groups[k]={d.element_name:fdict['header_table'][d.element_name]} 

                    except KeyError:   
                        x=numpy.zeros(di['recordindex'].shape[0], dtype=numpy.dtype(ttrans(d.kind,kinds=okinds)))
                        x.fill(numpy.nan)
                        groups[k][d.element_name]=({'hdrlen':di['recordindex'].shape[0]},x)

                elif k in ('station_configuration'): # station_configurationt contains info of all the stations, so this extracts only the one line for the wanted station with the numpy.where
                    try:                        
                        groups[k]={d.element_name:np.full( 1 , station_configuration_retrieved[d.element_name].values[0] ) }
                    except:
                        pass

                elif k in ('source_configuration'): # storing the source configuration info, e.g. original file name, 
                    try:   
                        groups[k]={d.element_name:np.full (1, fn,dtype='S{}'.format(len(fn)))}
                    except KeyError:
                        pass


                else : # this is the case where the cdm tables DO exist
                    groups[k]={d.element_name:cdm[k][d.element_name].values}
                    #try:   
                        #groups[k][d.element_name]=({k+'_len':len(cdm[k])}, cdm[k][d.element_name].values) # element_name is the netcdf variable name, which is the column name of the cdm table k 
                    #except KeyError:
                        #pass

                """ Tryin to add attributes, e.g. description and external tables """ 
                if k not in ('observations_table'):               
                    try:
                        groups[k][d.element_name].attrs['external_table']=d.external_table # defining variable attributes that point to other tables (3rd and 4th columns)
                        groups[k][d.element_name].attrs['description']=d.description  # it fails when trying with the observations_table 
                    except:
                        pass

                try:

                    if type(groups[k]) is dict:
                        idx=np.where(cdmd[k].element_name==d.element_name)[0][0]
                        if len(groups[k][d.element_name])>0:   
                            if type(groups[k][d.element_name][0]) in (str,bytes):
                                groups[k][d.element_name]=np.array(groups[k][d.element_name],dtype='S')
                            else:
                                groups[k][d.element_name]=np.array(groups[k][d.element_name])
                        else:
                            for m in fdict[k].values():
                                if len(m)>0:
                                    break
                            dt=gkinds[cdmd[k].kind[idx]]
                            if dt==numpy.datetime64:
                                dt=numpy.int64
                            groups[k][d.element_name]=np.empty(len(m),dtype=dt)
                            groups[k][d.element_name].fill(np.nan)
                        gkev=groups[k][d.element_name]
                    else:
                        gkev=groups[k][d.element_name].values
                    if type(gkev[0])==str:
                        s=gkev.shape
                        groupencodings[k][d.element_name]={'dtype':numpy.dtype('S80'),'compression': 'gzip','chunksizes':(min(100000, s[0] ) , 80 ) }
                    else:
                        groupencodings[k][d.element_name]={'compression': 'gzip'}

                    write_dict_h5(fno, groups[k], k, groupencodings[k], var_selection=[],mode='a', attrs= dic_obstab_attributes  )
                except KeyError:
                    #print('bad:',k,d.element_name)
                    pass

                print(k,d.element_name,time.time()-tt,' mem:',process.memory_info().rss//1024//1024)

    """ Storing the group encodings in a numpy dictionary to be reused by the merging script """
    np.save('groups_encodings',  groupencodings)

    return 0

def ubern_to_cdm(cdm, cdmd, output_dir, dataset, dic_obstab_attributes, fn):

    process = psutil.Process(os.getpid())
    t=time.time()

    dtypes={'station_id':numpy.int32,'latitude':str,'longitude':str,
            'altitude':str,
                                                 'rstype':'S4','datetime':numpy.int32,'date_flag':'S2','Station Name':'S60'}
    names=list(dtypes.keys())
    cdm['metadata_schroeder']=pd.read_csv(os.path.expanduser('../data/tables/vapor.library.2'),sep=':',header=0,
                                          dtype=dtypes,names=names)
    for l in 'latitude','longitude','altitude':
        cdm['metadata_schroeder'][l]=pd.to_numeric(cdm['metadata_schroeder'].pop(l), errors='coerce')
    cdm['iri_rstype_map']=pd.read_csv('../data/tables/iri_rstype_map_unibe.dat',sep='\t',
                                      dtype={'riname':np.dtype('S10'),'vapor_name':'S20','wmo_C2011_code':np.float32,'ri_long_name':'S80'})
    cdm['iri_rstype_map']['wmo_C2011_code']=numpy.int32(cdm['iri_rstype_map']['wmo_C2011_code'])
    cdmd['iri_rstype_map']=pd.DataFrame(data=None, columns=cdmd['observations_table'].columns)
    l=0
    for k in cdm['iri_rstype_map'].columns:
        cdmd['iri_rstype_map']=cdmd['iri_rstype_map'].append({'element_name':k,
                                                              'kind':type(cdm['iri_rstype_map'][k].values[0]),
                                       'description':k},ignore_index=True)
        l+=1

    cdms=pd.read_csv(os.path.expanduser('../data/tables/vapor.instruments.all'),sep=':',names=('sensor_id','comments'))
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


    #flist=glob.glob(fn+'/COMP004*')  # cant understand this piece
    flist=glob.glob(fn+'/COMP/COMP*') # fixed by F 
    campaigns=pd.read_excel(fn+'/Table1_Campaigns.xls')
    fdicts={}
    for fn in flist:

        os.chdir(fn)
        comp=fn[-7:]
        try:
            proflist=pd.read_excel(comp+'_ListProfiles.xls')
        except:
            try:
                proflist=pd.read_excel(comp+'_ListProfiles.xlsx')
            except:
                proflist=None

        try:
            soundlist=pd.read_excel(comp+'_ListSoundings.xls')
        except:
            soundlist=None

        #mon='mo'
        #if proflist is None:
            #print('No proflist')
            #proflist=soundlist
            #mon='mon'

        ll= sorted(glob.glob(comp+'_Aggregated/*_s_*.tsv'))

        for r in ll:
            if 'COMP002' in r:
                continue
            if 'xlsx' in r:
                data=pd.read_excel(r)
                continue
            else:
                try: 
                    data=pd.read_csv(r,sep='\t')
                except:
                    #print(r+' could not be read')
                    continue

            process_ubern_aggregated(fdicts,cdmd,cdm,r,proflist,campaigns,data)

        ll= sorted(glob.glob(comp+'_SoundingData/*m.txt'))
        if not ll:
            ll= sorted(glob.glob(comp+'_SoundingData/*s.txt'))

        for r in ll:
            try:
                sp='\n'
                with open(r,'rb') as a:
                    data=a.read().decode('latin1').split(sp)
                crfound=False
                for i in range(len(data)):
                    if '\r' in data[i]:
                        data[i]=data[i][:-1]
                    crfound=True
                if crfound:
                    x=0#print('CR found')
                else:
                    print('CR not found')

                process_ubern(fdicts,cdmd,cdm,r,proflist,campaigns,data)
            except MemoryError:
                pass

    fnos=[]
    for key,fdict in fdicts.items():

        idx=numpy.lexsort((fdict['observations_table']['z_coordinate'],fdict['observations_table']['report_id']))
        for i in fdict['observations_table'].keys():
            try:
                fdict['observations_table'][i]=numpy.array(fdict['observations_table'][i])[idx]
                #print(i,'success')
            except:
                print(i,'fail')

        x=np.array(fdict['observations_table']['date_time'][:])
        idx=[]
        l=0
        iold=b''
        for i in fdict['observations_table']['report_id']:
            if iold!=i:
                idx.append(l)
                iold=i
            l+=1
        #idx.append(len(fdict['observations_table']['report_id']))
        idx=np.array(idx)
        di=xr.Dataset() 
        di['recordindex']=({'record':idx.shape[0]},idx)
        di['recordtimestamp']=({'record':idx.shape[0]}, 
                               numpy.array(fdict['header_table']['record_timestamp']))
        di['recordtimestamp'].attrs['units']='seconds since 1900-01-01 00:00:00'

        l=0
        iold=0
        refdate=datetime(1900,1,1)
        z=[]
        for i in fdict['header_table']['record_timestamp']:
            if i//86400!=iold:
                dat=refdate+timedelta(seconds=i)
                if l==di['recordindex'].shape[0]-1:
                    z.append((dat.year*10000+dat.month*100+dat.day,di['recordindex'].values[l],len(fdict['observations_table']['date_time'])))    
                else:
                    z.append((dat.year*10000+dat.month*100+dat.day,di['recordindex'].values[l],di['recordindex'].values[l+1]))
                iold=i//86400
            l+=1

        z=np.array(z)
        di['dateindex']=({'days':z.shape[1],'drange':z.shape[0]},z) # date, index of the first occurrance, index of the last

        wigos='-'.join(fdict['header_table']['report_id'][0].decode().split('-')[:-1])
        wigos=fdict['header_table']['primary_station_id'][0].decode()
        fno='/'.join(fn.split('/')[:-2])+'/nc/'+wigos+'.nc'
        fnos.append(fno)
        os.remove(fno)
        di.to_netcdf(fno,format='netCDF4',engine='h5netcdf',mode='w')

        print(time.time()-t)
        tt=time.time()

        groups={}
        groupencodings={}
        for k in cdmd.keys(): # loop over all the table definitions 
            groupencodings[k]={} # create a dict of group econding

            for i in range(len(cdmd[k])): 
                d=cdmd[k].iloc[i] 
                if k in ('observations_table'):
                    try:
                        groups[k]={d.element_name:fdict['observations_table'][d.element_name]}
                    except KeyError:
                        x=numpy.zeros(fdict['observations_table']['date_time'].shape[0],
                                      dtype=numpy.dtype(ttrans(d.kind,kinds=okinds) ) )
                        x.fill(numpy.nan)
                        groups[k][d.element_name]=x


                elif k in ('header_table'):
                    try:
                        groups[k]={d.element_name:fdict['header_table'][d.element_name]} 

                    except KeyError:   
                        x=numpy.zeros(di['recordindex'].shape[0], dtype=numpy.dtype(ttrans(d.kind,kinds=okinds)))
                        x.fill(numpy.nan)
                        groups[k][d.element_name]=({'hdrlen':di['recordindex'].shape[0]},x)

                elif k in ('station_configuration'): # station_configurationt contains info of all the stations, so this extracts only the one line for the wanted station with the numpy.where
                    try:                        
                        groups[k]={d.element_name:np.full( 1 , station_configuration_retrieved[d.element_name].values[0] ) }
                    except:
                        pass

                elif k in ('source_configuration'): # storing the source configuration info, e.g. original file name, 
                    try:   
                        groups[k]={d.element_name:np.full (1, fn,dtype='S{}'.format(len(fn)))}
                    except KeyError:
                        pass


                else : # this is the case where the cdm tables DO exist
                    groups[k]={d.element_name:cdm[k][d.element_name].values}
                    #try:   
                        #groups[k][d.element_name]=({k+'_len':len(cdm[k])}, cdm[k][d.element_name].values) # element_name is the netcdf variable name, which is the column name of the cdm table k 
                    #except KeyError:
                        #pass

                """ Tryin to add attributes, e.g. description and external tables """ 
                if k not in ('observations_table'):               
                    try:
                        groups[k][d.element_name].attrs['external_table']=d.external_table # defining variable attributes that point to other tables (3rd and 4th columns)
                        groups[k][d.element_name].attrs['description']=d.description  # it fails when trying with the observations_table 
                    except:
                        pass

                try:
                    if type(groups[k]) is dict:
                        idx=np.where(cdmd[k].element_name==d.element_name)[0][0]
                        if len(groups[k][d.element_name])>0:   
                            if type(groups[k][d.element_name][0]) in (str,bytes):
                                groups[k][d.element_name]=np.array(groups[k][d.element_name],dtype='S')
                            else:
                                groups[k][d.element_name]=np.array(groups[k][d.element_name])
                        else:
                            for m in fdict[k].values():
                                if len(m)>0:
                                    break
                            dt=gkinds[cdmd[k].kind[idx]]
                            if dt==numpy.datetime64:
                                dt=numpy.int64
                            groups[k][d.element_name]=np.empty(len(m),dtype=dt)
                            groups[k][d.element_name].fill(np.nan)
                        gkev=groups[k][d.element_name]
                    else:
                        gkev=groups[k][d.element_name].values
                    if type(gkev[0])==str:
                        s=gkev.shape
                        groupencodings[k][d.element_name]={'dtype':numpy.dtype('S80'),'compression': 'gzip','chunksizes':(min(100000, s[0] ) , 80 ) }
                    else:
                        groupencodings[k][d.element_name]={'compression': 'gzip'}

                    write_dict_h5(fno, groups[k], k, groupencodings[k], var_selection=[],mode='a', attrs= dic_obstab_attributes  )
                except KeyError:
                    #print('bad:',k,d.element_name)
                    pass


                #print(k,d.element_name,time.time()-tt,' mem:',process.memory_info().rss//1024//1024)

    fnu=[]
    for f in fnos:
        if f not in fnu:
            fnu.append(f)

    for f in fnu:
        counter=0
        for g in fnos:
            if g==f:
                counter+=1
        if counter>1:
            print(f,counter)

    print(fnos)
    with open(fn+'/../../nc/mapping.tsv','w') as fo:
        data=['Comparison\tFile\tSensor\tReferenceSensor']
        i=0
        for k,v in fdicts.items():
            if '_' in k:
                s=k+'\t'+fnos[i].split('/')[-1]+'\t'+v['sensor_id'].decode()+'\t'+v['reference_sensor_id'].decode()
            else:
                s=k+'\t'+fnos[i].split('/')[-1]+'\t'+v['sensor_id'].decode()+'\t'
            data.append(s)
            i+=1
        data.sort()
        fo.write('\n'.join(data)+'\n')

    """ Storing the group encodings in a numpy dictionary to be reused by the merging script """
    np.save('groups_encodings',  groupencodings)

    return 0


def load_cdm_tables():
    """ Load the cdm tables into Panda DataFrames, reading the tables from the cdm GitHub page FF To do 

    # Uncomment to get the list of all the .csv files present at the url specified
    # url = 'https://github.com/glamod/common_data_model/tree/master/table_definitions'
    # cdmtabledeflist = csvListFromUrls(url)
    """
    tpath = '/srvfs/home/uvoggenberger/CEUAS/CEUAS/public/harvest/data' # os.getcwd() + '/../data' ## overwritten by hardcoded path
    cdmpath='https://raw.githubusercontent.com/glamod/common_data_model/master/tables/' # cdm tables            

    """ Selecting the list of table definitions. Some of the entires do not have the corresponding implemented tables """
    cdmtabledeflist=['id_scheme', 'crs', 'station_type', 'observed_variable', 'station_configuration', 'station_configuration_codes', 'observations_table', 
                     'header_table', 'source_configuration', 'sensor_configuration', 'units' , 'z_coordinate_type']  
    cdm_tabdef = dict()
    for key in cdmtabledeflist:
        url='table_definitions'.join(cdmpath.split('tables'))+key+'.csv' # https://github.com/glamod/common_data_model/tree/master/table_definitions/ + ..._.dat 
        f=urllib.request.urlopen(url)
        col_names=pd.read_csv(f, delimiter='\t',quoting=3,nrows=0,comment='#')
        f=urllib.request.urlopen(url)
        tdict={col: str for col in col_names}
        cdm_tabdef[key]=pd.read_csv(f,delimiter='\t',quoting=3,dtype=tdict,na_filter=False,comment='#')


    """ Selecting the list of tables. 'station_configuration_codes','observations_table','header_table' are not implemented in the CDM GitHub"""        
    cdmtablelist=['id_scheme', 'crs', 'station_type', 'observed_variable', 'station_configuration_codes','units']        
    cdm_tab=dict() # dictionary where each key is the name of the cdm table, and the value is read from the .dat file    
    for key in cdmtablelist:
        f=urllib.request.urlopen(cdmpath+key+'.dat')
        col_names=pd.read_csv(f,delimiter='\t',quoting=3,nrows=0)
        f=urllib.request.urlopen(cdmpath+key+'.dat')
        tdict={col: str for col in col_names}
        if key != 'observed_variable':
            cdm_tab[key]=pd.read_csv(f,delimiter='\t',quoting=3,dtype=tdict,na_filter=False)
        else:
            cdm_tab[key]=pd.read_csv(f,delimiter='\t',quoting=3,dtype=tdict,na_filter=False, nrows=150)  # TO DO FIX, latest update broke compatibility, cannot read last lines

    """ Adding the  tables that currently only have the definitions but not the implementation in the CDM, OR    need extensions """  
    cdm_tabdef['header_table']          = pd.read_csv(tpath+'/table_definitions/header_table.csv',delimiter='\t',quoting=3,comment='#')
    #cdm_tabdef['observations_table'] = pd.read_csv(tpath+'/table_definitions/observations_table.csv',delimiter='\t',quoting=3,comment='#')

    id_scheme={ cdm_tabdef['id_scheme'].element_name.values[0]:[0,1,2,3,4,5,6],
                cdm_tabdef['id_scheme'].element_name.values[1]:['WMO Identifier','Volunteer Observing Ships network code',
                                                                         'WBAN Identifier','ICAO call sign','CHUAN Identifier',
                                                             'WIGOS Identifier','Specially constructed Identifier']}

    cdm_tab['id_scheme']=pd.DataFrame(id_scheme)
    cdm_tab['crs']=pd.DataFrame({'crs':[0],'description':['wgs84']})

    """ Here we add missing entries, e.g. in the z_coordinate_type for the pressure levels in Pascal (the available CDM table in the glamod GitHub rep.  contains onle the altitude in [meter] """
    cdm_tab['station_type']=pd.DataFrame({'type':[0,1],'description':['Radiosonde','Pilot']}) 
    cdm_tab['z_coordinate_type']=pd.DataFrame({'type':[0,1],'description':['height (m) above sea level','pressure (Pa)']})  # only the m above sea level is available currently in the GitHub cdm table, added pressure 


    """ Make dictionary of variables and attributes for the observations table """ 
    dic_obstab_attributes = {}
    for index, row in cdm_tabdef['observations_table'].iterrows():
        dic_obstab_attributes[row['element_name'] ] = {}
        dic_obstab_attributes[row['element_name'] ]['description'] = row.description 
        dic_obstab_attributes[row['element_name'] ]['external_table'] = row.external_table 

    #dic_obs['date_time'] = ['units',  'seconds since 1900-01-01 00:00:00' ]

    if not os.path.isfile('dic_obstab_attributes.npy'): 
        np.save( 'dic_obstab_attributes' , dic_obstab_attributes )

    return cdm_tabdef, cdm_tab, tdict , dic_obstab_attributes 


def csvListFromUrls(url=''):
    """ Return a list of csv files, as fond in the url on the cdm GitHub """   
    urlpath = urlopen(url)
    string = urlpath.read().decode('utf-8')
    split = string.split(' ')
    csv_files_list = [m.replace('"','') for m in [n.split('title="')[1] for n in split if '.csv' in n and "title" in n] ] 
    return csv_files_list


def filelist_cleaner(lista, dataset=''):
    """ Removes unwanted files that might be present in the dataset directories """
    if dataset == 'ncar':
        cleaned = [ l for l in lista if '.nc' not in l ]
    if dataset == 'bufr':
        cleaned = [ l for l in lista if '.bfr' in l ]
    if dataset == 'maestro':
        cleaned = [ l for l in lista if '.bufr' in l ]
    if dataset =='amma':
        cleaned = [ l for l in lista if '.' not in l ]
    if 'era5' in dataset:
        cleaned = [ l for l in lista if '.nc' not in l and '.conv.' in l ]
    else:
        cleaned = lista
    return cleaned


def clean_station_configuration(cdm_tab ):
    """ Replace wrong characters from the station configuration tables """
    subs={'o':[240,242,243,244,245,246,248],'O':[210,211,212,213,214,216],
          'a':[224,225,226,227,228,229,230],'A':[192,193,194,195,196,197,198],
          'u':[249,250,251,252,253],'U':[217,218,219,220],
          'i':[236,237,238,239],'I':[204,205,206,207,304],
          'S':[350],'n':[241],'c':[231],'C':[199],'e':[232,233,234,235],'E':[200,201,202,203]}
    for k in cdm_tab['station_configuration'].columns:
        if type(cdm_tab['station_configuration'][k][0]) is str:
            try:               
                cdm_tab['station_configuration'][k].values[:]=cdm_tab['station_configuration'][k].values[:].astype('S')
            except:
                for l in range(cdm_tab['station_configuration'][k].values.shape[0]):
                    try:                       
                        cdm_tab['station_configuration'][k].values[l]=np.bytes_(cdm_tab['station_configuration'][k].values[l])
                    except:
                        for m in range(len(cdm_tab['station_configuration'][k].values[l])):
                            mychar=cdm_tab['station_configuration'][k].values[l][m]
                            if ord(mychar)>128:
                                for n,v in subs.items():
                                    if ord(mychar) in v:
                                        cdm_tab['station_configuration'][k].values[l]=n.join(cdm_tab['station_configuration'][k].values[l].split(mychar))

                        cdm_tab['station_configuration'][k].values[l]=np.bytes_( (cdm_tab['station_configuration'][k].values[l] ).encode('utf-8') ) 
                cdm_tab['station_configuration'][k]=np.bytes_(cdm_tab['station_configuration'][k])

    print('Cleaned station_configuration')





# on srvx1, srvx8 
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Make CDM compliant netCDFs")
    parser.add_argument('--dataset' , '-d', 
                        help="Select the dataset to convert. Available options: all, era5_1, era5_1759, era5_1761, bufr, maestro, igra2, ncar, test. If not selected or equal to 'test', the script will run the example files in the /examples directory."  ,
                    type = str)

    parser.add_argument('--output' , '-o', 
                        help="Select the output directory. If not selected, converted files will be stored in the 'converted_files' directory."  ,
                    default = 'converted_files',
                    type = str)    

    parser.add_argument('--files' , '-f',
                        help = "File to be processed."  ,
                    default = '',
                        type = str)    

    parser.add_argument('--run_missing' , '-r',
                        help = "Only run missing year"  ,
                    default = 'False',
                        type = str)    

    parser.add_argument('--max_year' , '-maxy',
                        help = "Most recent correctly harvested year"  ,
                    default = '2025',
                        type = str)    

    parser.add_argument('--primary_id' , '-p',
                        help = "Primary id of the station"  ,
                    default = '',
                        type = str)    

    parser.add_argument('--kind' , '-k',
                        help = "Type of station [regular,mobile,orphan]"  ,
                        default = 'regular',
                        type = str)    

    parser.add_argument('--parameter_file' , '-pf', 
                    help="Chose parameter file"  ,
                    type = str,
                    default = '.' )

    
    args = parser.parse_args()
    dataset = args.dataset 
    out_dir = args.output
    Files = args.files
    kind = args.kind
    parameter_file = args.parameter_file

    if parameter_file != '.': 
        sys.path.append(parameter_file.split('modded')[0])
        from modded_harvester_yearsplit_parameters import *
    else:
        from harvester_yearsplit_parameters import *  


    vlist= [
        
        'era5_1', 'era5_2', 'era5_3188', 'era5_1759', 'era5_1761', 
            
            'bufr_cnr', 'bufr', 'igra2', 'ncar', 
            'maestro',
            'woudc',
            'amma',
            
            'era5_1_mobile',
            'era5_2_mobile',
            'igra2_mobile',
            
            'giub' ,
            
            'hara',
            'npsound',
            'shipsound' ,
            
            # intercomparisons
            'mauritius',
            'mauritius_digitized',
            'yangjiang',

            'ubern',
            'nasa']

    ### to fix inconsistency in station kind and dataset, if ever happens
    
    if dataset in ['era5_1_mobile','era5_2_mobile','igra2_mobile']:
        kind = 'mobile'
        
    if dataset not in vlist:
        print('wrong dataset', dataset)
        raise ValueError(" The selected dataset is not valid. Please choose from ['era5_1', 'era5_1759', 'era5_1761', 'era5_3188', 'bufr', 'maestro' 'igra2', 'ncar', 'amma', 'giub', 'hara', 'npsound', 'shipsound' , 'era5_1_mobile', 'era5_2_mobile' , 'yangjiang']  ")    

    """ Loading the CDM tables into pandas dataframes """
    cdm_tabdef  , cdm_tab, tdict , dic_obstab_attributes= load_cdm_tables()

    """ Paths to the output directory """    
    if not os.path.isdir(out_dir):
        os.system('mkdir ' + out_dir )       

    output_dir = out_dir + '/' + dataset           
    if not os.path.isdir(output_dir):
        os.system('mkdir ' + output_dir )

    stat_conf_path = '/srvfs/home/uvoggenberger/CEUAS/CEUAS/public/harvest/data/station_configurations/'     # overwritten path

    if kind in ['orphan']:
        stat_conf_file = stat_conf_path +   '/' + dataset+ '_orphans_station_configuration_extended.csv' 

    elif kind in ['mobile']:
        if dataset not in ['igra2_mobile']:
            stat_conf_file = stat_conf_path +   '/' + dataset+ '_station_configuration_extended.csv'  
        else:
            stat_conf_file = stat_conf_path +   '/igra2_mobile_station_configuration_extended.csv'  
            
    #if 'era5_1_mobile' in dataset:
    #    stat_conf_file = stat_conf_path +   '/era5_1_mobile_station_configuration_extended.csv'    

    elif 'mauritius' in dataset:
        stat_conf_file = stat_conf_path +   '/station_configuration_mauritius.dat'    

    elif 'yangjiang' in dataset:
        if not os.path.isfile( stat_conf_path +   '/station_configuration_yangjiang.csv' ):
            df = pd.read_csv(  stat_conf_path + '/CUON_station_configuration_extended.csv', sep = '\t') 
            stat_conf_file = df.loc[df['station_name'] == 'YANGJIANG']
            stat_conf_file.to_csv( stat_conf_path+'/station_configuration_yangjiang.csv' , sep = '\t')
        else:
            stat_conf_file =  stat_conf_path + '/station_configuration_yangjiang.csv'

    else:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
        stat_conf_file = stat_conf_path +   '/' + dataset + '_station_configuration_extended.csv'    

    tdict['primary_id'] = np.bytes_
    tdict['secondary_id'] = np.bytes_    

    '''
    if dataset == 'shipsound' or dataset == 'npsound': # there is no station_configuration for this one 
        cdm_tab['station_configuration'] = []
    else:
        cdm_tab['station_configuration']=pd.read_csv(stat_conf_file,  delimiter='\t', quoting=3, dtype=tdict, na_filter=False, comment='#')
        clean_station_configuration(cdm_tab)              
    '''
    cdm_tab['station_configuration']=pd.read_csv(stat_conf_file,  delimiter='\t', quoting=3, dtype=tdict, na_filter=False, comment='#')
    clean_station_configuration(cdm_tab)       
    a=0

    ### splitting files from input string 
    Files = Files.split(',')

    if not os.path.isdir(out_dir):
        os.system('mkdir ' + out_dir ) 

    output_dir = out_dir + '/' + dataset      
    if not os.path.isdir(output_dir):
        os.system('mkdir ' + output_dir )

    for File in Files:

        primary_id = File.split('/')[0]

        File = File.split('/')[-1] 
        original_file_name = File

        if dataset in [ 'era5_1759' , 'era5_1761']:
            File = File.replace('.conv.', '.conv._')+'.gz'

        elif dataset in ['era5_1', 'era5_1_mobile']:
            print(File)
            File = File.replace(File.split('.')[-2],'??????')+'.gz' ## changed txt.gz to .gz
            original_file_name = original_file_name.replace(original_file_name.split('.')[-2],'??????') # original_file_name.replace('.conv._','.conv.??????.')

        # elif dataset in ['era5_1_mobile']:
        #     File = File.replace(File.split('.')[-3],'??????') ## changed txt.gz to .gz # File.replace('.conv.','.conv.??????.')+'.gz'  # +'.txt.gz' -> '.gz'
        #     original_file_name = original_file_name.replace(original_file_name.split('.')[-3],'??????') # original_file_name.replace('.conv._','.conv.??????.')

        elif dataset in ['igra2']:
            File = File # + '-data.txt'

        elif dataset in ['igra2', 'igra2_mobile']:
            File = File # + '-data.txt'

        elif dataset in ['era5_2', 'era5_2_mobile']:
            File = File + '.gz'
            #File = File

        File = datasets_path[dataset] + '/' + File 
        tt=time.time()              

        if parameter_file != '.': 
            sys.path.append(parameter_file.split('modded')[0])
            from modded_harvester_yearsplit_parameters import run_only_missing_stations   # will run only missing station and year; if False, will rerun the whole station 
            from modded_harvester_yearsplit_parameters import max_year_to_process, min_year_to_process # min and max year to harvest (independent of the station data)

        else:
            from harvester_yearsplit_parameters import run_only_missing_stations   # will run only missing station and year; if False, will rerun the whole station 
            from harvester_yearsplit_parameters import max_year_to_process, min_year_to_process # min and max year to harvest (independent of the station data)

        ### here: check if a file exists for a specific year and skip rerunning
        if run_only_missing_stations:
            out_d = out_dir + '/' + dataset + '/' + primary_id 
            if os.path.isdir(out_d):

                file = [out_d+'/' + f for f in os.listdir(out_d) if 'correctly' in f and original_file_name in f ] ### check if file with the correctly processed years exists 
                if len(file) > 0:
                    print('::: Found existing data for station/file ')
                    try:
                        years = open(file[0]).readlines()
                        last_year = max([ int(y.replace('\n','')) for y in years ])
                        min_year_to_process = last_year
                    except:
                        pass 
            else:
                pass

        ### ERA5 BLOCK 
        print(' Starting harvesting from the YEAR::: ' , min_year_to_process )
        out_dir = ''
        if 'era5' in dataset:   
            change_lat = False
            change_lon = False
            if '1759' in dataset:  # oading the files list to apply the WBAN latitude correction (missing minus sign)
                lat_mismatch = stat_conf_path + '/era5_1759_WBAN_latitude_mismatch.dat'                
                wban_lat_mismatch=pd.read_csv( lat_mismatch, delimiter='\t' ) 
                files_mismatch = list(wban_lat_mismatch['file'].values) 

                f = File.split('/')[-1].replace('.gz','').replace('._','.') # -f /raid60/scratch/leo/scratch/era5/odbs/1759/era5.1759.conv.2:82405.gz -d era5_1759 -o OUTPUT

                # decide if positive latitute should be changed to negative 
                if f in files_mismatch:
                    change_lat = True
                else:
                    pass

            if 'era5_2' in dataset:  # oading the files list to apply the WBAN latitude correction (missing minus sign)
                lat_mismatch = stat_conf_path + '/era5_2_lat_lon_mismatch.csv'                
                wban_mismatch=pd.read_csv(lat_mismatch) 
                files_mismatch = list(wban_mismatch['file'].values) 

                f = File.split('/')[-1]

                # decide if positive latitute should be changed to negative 
                if f in files_mismatch:
                    if wban_mismatch[wban_mismatch.file == f]['lat/lon'].values[0] == 'latitude':
                        change_lat = wban_mismatch[wban_mismatch.file == f]['value'].values[0]
                    else:
                        change_lon = wban_mismatch[wban_mismatch.file == f]['value'].values[0]
                else:
                    pass

            # extracting once for all the list of yearly files for this station 
            #fnl=File.split('/')   #TODO remove, useless ???
            #fnl[-1]='ch'+fnl[-1]
            fns=sorted(glob.glob(File))

            if len(fns) == 0 and dataset in ['era5_1_mobile']:
                f = File.replace('1_mobile','1').replace("??????._","??????.")  
                fns=sorted(glob.glob( f ) )
                
                
            # problem: era5_1 can contain also era5_1_mobile
            # so must check both 1/ and 1_mobile directory 
            
            ### example '/mnt/users/scratch/leo/scratch/era5/odbs/1/era5.conv.??????.38475.txt.gz'
            ###  fns=sorted(glob.glob(/mnt/users/scratch/leo/scratch/era5/odbs/1/era5.conv.??????.38475.txt.gz) - > ['/mnt/users/scratch/leo/scratch/era5/odbs/1/era5.conv.197911.38475.txt.gz', '/mnt/users/scratch/leo/scratch/era5/odbs/1/era5.conv.198601.38475.txt.gz']

            #if dataset in ['era5_1', 'era5_1_mobile', 'era5_2']:  # these are already split by year in the original files 
            if dataset in ['era5_1', 'era5_1_mobile']:  # these are already split by year in the original files 

                for year in range(min_year_to_process, max_year_to_process):
                # for year in range(2022, 2023):       ### !!!
                    dummy_writing = None
                    if run_only_missing_stations and year <= int(min_year_to_process):
                        print('Skipping already processed year: ' , year)
                    else:
                        year = str(year)
                        ff = [f for f in fns if year in  f.split('/')[-1].split('.')[2][0:4] ]  #here: pre-selecting the files already split by year 

                        force_this_run = True # set False to stop at errors
                        
                        if len(ff) > 0:
                                fbds, min_year_data  = read_odb_to_cdm(output_dir, dataset, dic_obstab_attributes, File, ff)
                                if fbds.empty:
                                    print("No valid data in source file ::: " , File)
                                    a = open('logs/' + dataset + '_empty_files.txt' , 'a+')
                                    a.write(File + '\t' + dataset + '\n')
                                    a.close()
                                    continue
                                
                                if force_this_run:
                                    try:
                                        dummy_writing = write_odb_to_cdm( fbds, cdm_tab, cdm_tabdef, output_dir, dataset,  dic_obstab_attributes, File, ff, change_lat, change_lon, year)
                                        print('DONE --- ' , year )
                                    except:
                                        print('FAILING ::: ' ,  File, ' ', ff )
                                        a = open("logs/" + dataset + '_wrong_files.txt' , 'a+')
                                        a.write(File + '\t' + dataset + '\n')
                                        a.close()
                                        
                                        #sys.exit() ## TO DO HERE TODO CHANGE                                         
                                else:
                                    dummy_writing = write_odb_to_cdm( fbds, cdm_tab, cdm_tabdef, output_dir, dataset,  dic_obstab_attributes, File, ff, change_lat, change_lon, year)
                                    print('DONE --- ' , year )

                        else:
                            print("No ERA5 1 files for year " , year )   
                            dummy_writing = None
                if dummy_writing != None:
                    a=open(dummy_writing, 'a+')
                    a.write('completed\n')
                    a.close()

            else:
                fbds, min_year_data  = read_odb_to_cdm(output_dir, dataset, dic_obstab_attributes, File, fns)
                for year in range(min_year_to_process, max_year_to_process):                
                    if year < min_year_to_process:
                        print('Year ' + year + ' but no data before ' + str(min_year_data) )
                        continue
                    if run_only_missing_stations and year <= int(min_year_to_process):
                        print('Skipping already processed year: ' , year )   

                    dummy_writing = write_odb_to_cdm( fbds, cdm_tab, cdm_tabdef, output_dir, dataset,  dic_obstab_attributes, File, fns, change_lat, change_lon, year)
                a=open(dummy_writing, 'a+')
                a.write('completed\n')
                a.close()                

        elif 'nasa' in dataset:   
            ir_to_cdm( cdm_tab, cdm_tabdef, output_dir, dataset,  dic_obstab_attributes, File)

        elif 'ubern' in dataset:   
            ubern_to_cdm( cdm_tab, cdm_tabdef, output_dir, dataset,  dic_obstab_attributes, File)      

        elif 'yangjiang' in dataset:
            sensors = list( np.unique( [ s for s in os.listdir(datasets_path['yangjiang'] ) if '.' not in s  ] ) )
            timestamps_df = pd.read_csv(datasets_path['yangjiang'] + '/intercomparison_2010_meta_data.csv')
            timestamps_df = timestamps_df[['FlightNo.' , 'Date_Time']]

            #dummy_writing = write_df_to_cdm(df, stat_conf_check, station_configuration_retrieved, cdm_tab, cdm_tabdef, output_dir, dataset,  dic_obstab_attributes, 'mauritius_2005_intercomparison', '2005', sensor=False)             

            for s in sensors:
                #if s != 'Modem':
                #    continue 
                df, stat_conf_check, station_configuration_retrieved, primary_id, min_year = read_df_to_cdm(cdm_tab, dataset, datasets_path['yangjiang']+'/' + s , metadata= timestamps_df )                    
                dummy_writing = write_df_to_cdm(df, stat_conf_check, station_configuration_retrieved, cdm_tab, cdm_tabdef, output_dir, dataset,  dic_obstab_attributes, 'yangjiang_2010_intercomparison', '2010', sensor=s) 

        elif 'mauritius_digitized' in dataset:
            sensors = list( np.unique( [ s.split('_')[1].replace('.csv','') for s in os.listdir(datasets_path['mauritius_digitized'] +'/temp') ] ) )
            df, stat_conf_check, station_configuration_retrieved, primary_id, min_year_data = read_df_to_cdm(cdm_tab, dataset, datasets_path['mauritius_digitized'])

            #write combined file 
            dummy_writing = write_df_to_cdm(df, stat_conf_check, station_configuration_retrieved, cdm_tab, cdm_tabdef, output_dir, dataset,  dic_obstab_attributes, 'mauritius_2005_intercomparison', '2005', sensor=False) 

            # write single sensor file
            for s in sensors:
                dummy_writing = write_df_to_cdm(df, stat_conf_check, station_configuration_retrieved, cdm_tab, cdm_tabdef, output_dir, dataset,  dic_obstab_attributes, 'mauritius_2005_intercomparison', '2005', sensor=s) 

        else:
            #try:
            # include the sensor_type from bufr_cnr: 
            if 'bufr_cnr' in dataset:
                cdmfb_noodb['sensor_id'] = 'sonde_type@conv'

            df, stat_conf_check, station_configuration_retrieved, primary_id, min_year_data= read_df_to_cdm(cdm_tab, dataset, File)

            if not isinstance(df, pd.DataFrame):
                print("=== Found an empty file " + File + "     ==>  nothing to harvest ")
                a=open(output_dir + '/' + dataset + '_empty_files.dat', 'a+')
                a.write(File + '\n')
                a.close()
                continue

            for year in range(min_year_to_process, max_year_to_process):         
                if year < min_year_to_process:
                    print('Year ' + year + ' but no data before ' + str(min_year_data) )
                if run_only_missing_stations and year <= int(min_year_to_process):
                    print('Skipping already processed year: ' , year )       
                dummy_writing = write_df_to_cdm(df, stat_conf_check, station_configuration_retrieved, cdm_tab, cdm_tabdef, output_dir, dataset,  dic_obstab_attributes, File, year) 
                print(' ***** Convertion of  ' , File ,  '    ' , dummy_writing , '    year     ' ,    year ,  '         completed after {:5.2f} seconds! ***** '.format(time.time()-tt))   
                
            #except:
            a=open(dummy_writing, 'a+')
            a.write('completed\n')
            a.close()            
        print(' ***** Convertion of  ' , File ,  '    ' , dummy_writing , '    completed after {:5.2f} seconds! ***** '.format(time.time()-tt))   

""" Examples for running 

# use monthly input files, can be read in parallel
small file

-f '/mnt/users/scratch/leo/scratch/era5/odbs/1/era5.conv.??????.11035.txt.gz'  -d era5_1 -o test_float
-f '/mnt/users/scratch/leo/scratch/era5/odbs/1_mobile/era5.conv.202108._4DC8UUK.txt.gz'  -d era5_1_mobile  -o COP2
-f /mnt/users/scratch/leo/scratch/era5/odbs/1759/era5.1759.conv._6:99024.gz -d era5_1759 -o COP2

-f '/mnt/users/scratch/leo/scratch/era5/odbs/1_mobile/era5.conv.??????.__FNPH.txt.gz' -d era5_1_mobile  -o COP2

-f '/mnt/users/scratch/leo/scratch/era5/odbs/2/era5.conv._UMEY.txt.gz' -d era5_2 -o COP2
-f /mnt/users/scratch/leo/scratch/era5/odbs/2_mobile/era5.conv.195212 -d era5_2_mobile -o COP2

# /mnt/users/scratch/leo/scratch/era5/odbs/2/era5.conv._00856.gz

# era5.conv._57707 
# era5.conv._82930 

# era5.conv._1:48379
-f  /mnt/users/scratch/leo/scratch/era5/odbs/ai_bfr/era5.10106.bfr -d bufr -o COP2

NCAR
-f 0-20000-0-82930/uadb_windc_82930.txt -d ncar -o COP3 
-f 0-20001-0-11035/uadb_windc_11035.txt -d ncar -o COP3 


ERA5 1
-f  0-20000-0-38475/era5.conv._38475 -d era5_1  -o COP3 
-f  0-20000-0-82930/era5.conv._82930 -d era5_1  -o COP3

-f  0-20000-0-82930/era5.conv._74005 -d era5_1  -o COP3


IGRA2
-f  0-20000-0-82930/BRM00082930  -d igra2 -o COP3
-f 0-20000-0-72230/USM00072230 -d igra2 -o COP3
ERA5 1 MOBILE
-f  0-20999-0-DBDK/era5.conv._DBDK  -d era5_1_mobile  

ERA5 2 MOBILE  *** this 
-d era5_2_mobile   -f  0-20999-0-00630/era5.conv._00630  -r False   -o COP3

ERA5 2
-f 0-20000-0-72575/era5.conv._2:24126 -o COP2 -d era5_2

AMMA
-f /scratch/das/federico/databases_service2/AMMA_BUFR/HRT2006071606 OLD , read_BUFR

GIUB
# /scratch/das/federico/COP2_HARVEST_JAN2023/GIUB_22062023_reduced//5382A.txt_converted_csv.csv
# -f 0-20000-0-01025/5628.txt_converted_csv.csv  -d giub -o COP2 
# 0-20000-0-38987/5194A.txt_converted_csv.csv


HARA
-f 0-20666-0-25428/25428_stationfile.csv -d hara -o COP3

NPSOUND
-f dummy/np_11sound.dat_converted.csv -d npsound -o COP3

SHIPSOUND (only one file)
-f dummy/shipsound7696.csv -d shipsound -o COP2

MAURITIUS original excel files from vendors 0-20000-0-61995
-f 0-20000-0-61995/vaisala_ascents.csv  -d mauritius -o COP2 
-f 0-20000-0-61995/meisei_ascents.csv  -d mauritius -o COP2

MAURITIUS digitized by Ulrich 0-20000-0-61995 
-f 0-20000-0-61995/dummy  -d mauritius_digitized -o COP2 

YANGJIANG # 0-20000-0-59663 
/srvfs/home/uvoggenberger/scratch/intercomparison_2010
-f dummy  -d yangjiang -o COP2 

-f 0-20000-0-59663/dummy -d yangjiang -o COP2
"""


# MOBILE ERA5 1 weird character
# -f 0-20999-0-E§ML7/era5.conv._E§ML7 

