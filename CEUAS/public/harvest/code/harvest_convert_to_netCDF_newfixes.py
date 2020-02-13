#!/usr/bin/env python
import sys
import os.path 
import glob 
import psutil
import subprocess
import urllib.request
import xarray as xr
#import numpy
import h5py
from datetime import date, datetime,timedelta
import time
from multiprocessing import Pool
from netCDF4 import Dataset
import gzip
import pandas as pd    
from functools import partial
#from rasotools.utils import *
from numba import *
import cartopy.crs as ccrs
import argparse
#import copy
from io import StringIO
#import h5netcdf
import numpy as np
from eccodes import *
import warnings
import numpy
from numba import njit

warnings.simplefilter(action='ignore', category=FutureWarning) # deactivates Pandas warnings 


debug = False


""" Some colors for pretty printout """ 
red    = '\033[91m' 
cend   = '\033[0m'
blue   = '\033[34m'
green  = '\033[92m'
yellow = '\033[33m'

""" Possible variable types as listed int he CMD tables S80 is a fixed 80 char. string. Compression does not work well with normal strings. """
okinds={'varchar (pk)':np.dtype('|S80'),
               'varchar':np.dtype('|S80'), 
               'numeric':np.float32, 
               'int':np.int32,
               #'timestamp with timezone':np.datetime64,
               'timestamp with timezone':np.float32,
               
               'int[]*':list,
               'int[]':list,
               'varchar[]*':list,
               'varchar[]':list}

""" Variable types to be used in the compressed netCDF files """
kinds={'varchar (pk)':str,
             'varchar':str,
             'numeric':np.float32,
             'int':np.int32,
             'int(pk)' : np.int32,
             #'timestamp with timezone':np.datetime64,  # gives errore when wiritng out the netCDF file with the np.datetime64 
             'timestamp with timezone':np.float32,             
             'int[]*':list,
             'int[]':list,
             'varchar[]*':list,
             'varchar[]':list}

gkinds={'varchar (pk)':numpy.dtype('|S80'),
               'varchar':numpy.dtype('|S80'),
               'numeric':numpy.float32,
               'int':numpy.int32,
               'timestamp with timezone':numpy.datetime64,
               'int[]*':list,
               'int[]':list,
               'varchar[]*':list,
               'varchar[]':list}



def make_datetime(dvar,tvar):
    """ Converts into dat-time """
    dvari=dvar.astype(numpy.int)
    tvari=tvar.astype(numpy.int)
    df=pd.DataFrame({'year':dvar//10000,'month':(dvar%10000)//100,'day':dvar%100,
                        'hour':tvar//10000,'minute':(tvar%10000)//100,'second':tvar%100})
    dt=pd.to_datetime(df).values
    return numpy.array(dt-numpy.datetime64('1900-01-01'),dtype=int)//1000000000

def make_obsid(ivar):
    x=numpy.char.zfill(numpy.arange(ivar.shape[0]).astype('S8'), 10)
    return x

def make_recid(ivar):
    x=numpy.char.zfill(numpy.arange(ivar.values.shape[0]).astype('S5'), 10)
    return x

def make_obsrecid(fbvar,ivar):
    x=numpy.char.zfill(numpy.arange(ivar.values.shape[0]).astype('S5'), 10)
    y=numpy.zeros(fbvar.shape[0]).astype('S5')
    for i in range(ivar.values.shape[0]-1):
        y[ivar.values[i]:ivar.values[i+1]]=x[i]
    if ivar.values.shape[0]>1:
        y[ivar.values[i+1]:]=x[i+1]
    return y

def make_vars(ivar):
    tr=numpy.zeros(113,dtype=int) 
    """ translates odb variables number to Lot3 numbering convention """
    tr[1]=117  # should change
    tr[2]=85
    tr[3]=104
    tr[4]=105
    tr[7]=39 #spec hum
    tr[29]=38 #relative hum
    tr[59]=36 # dew point
    tr[111]=106 #dd
    tr[112]=107  #ff
    #
    tr[39]= 136 # 2m T according to proposed CDM standard
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
    #
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
       'units':[make_units,'varno@body'],
       'primary_station_id':'statid@hdr',       
       }


cdmfb_noodb={'observation_value':'obsvalue@body',
                          'observed_variable':'varno@body',
                          'z_coordinate_type':'vertco_type@body',
                          'z_coordinate':'vertco_reference_1@body',
                          'date_time':'record_timestamp', 
                          'release_time': 'report_timestamp' , # only available for igra2 
                          'longitude':'lon@hdr',
                          'latitude':'lat@hdr',
                          'observation_id':'observation_id' ,
                          'source_file':'source_file',
                          #'product_code': 'product_code' ,
                          'report_id':'report_id' ,
                          'number_of_pressure_levels' : 'num_lev',  # this entry does not exist in CDM header_table 
                          'units' : 'units',
                          'source_id': 'source_id',
                          'primary_id' : 'statid@hdr' ,                  # station_configuration
                          'primary_station_id':'statid@hdr'}        # header_table 



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
cdmvar_dic = {'temperature'          : { 'odb_var': 2      , 'cdm_unit': 5        , 'cdm_var': 85}     ,  # K
              
                         'wind_direction'      : { 'odb_var': 111   , 'cdm_unit': 110    , 'cdm_var': 106}   ,  # degree (angle)
                         'wind_speed'           : { 'odb_var': 112  , 'cdm_unit': 731     , 'cdm_var': 107 } ,  # m/s 
                         'uwind'                    : { 'odb_var': 3       , 'cdm_unit': 731     , 'cdm_var': 104}   ,  # m/s
                         'vwind'                    : { 'odb_var': 4       , 'cdm_unit': 731     , 'cdm_var': 105}    ,  # m/s
                         
                         'dew_point'             : { 'odb_var': 59             , 'cdm_unit': 5  , 'cdm_var': 36}     ,  # K
                         'dew_point_depression' : { 'odb_var': 299  , 'cdm_unit': 5     , 'cdm_var': 34}   ,  # K fake number, does not exhist in ODB file 
                         
                         'relative_humidity'  : { 'odb_var': 29     , 'cdm_unit': 0    , 'cdm_var': 38}     ,  # absolute value (ratio) 
                         'gph'                       : { 'odb_var': 1       , 'cdm_unit': 631         , 'cdm_var': 117}    ,  # need to verify from original data source

                         'pressure'               : { 'odb_var': 999    , 'cdm_unit': 32       , 'cdm_var': 57}      , # Pa  (it goes into z_coordinate type)
                          }


""" CDM variable codes for the corresponding ODB variables """
cdm_odb_var_dic = { 1    : 117    , # geopotential
                                   2    : 85        , # temperature K
                                   
                                   3    : 104    , # uwind m/s , upper air u component 
                                   4    : 105    ,  # vwind m/s
                                  111 : 106    , # degree (angle) , wind from direction 
                                  112 : 107    , # m/s , wind force 
                                  
                                  29   : 38      , # relative humidity in %
                                 59    : 36      , # dew point (available in ODB files )
                                 
                                 999  : 57      , # Pa  (NOTE: it goes into z_coordinate type, not in the observed variables)                                 
                                 99    : 34   , # dew_point depression (non existing in ODB files )
                          }


""" Common columns of the read dataframes (not from ODB files, which have their own fixed column names definitions) """
column_names = [ 'source_id', 'report_id',  'observation_id', 'record_timestamp' , 'iday', 'statid@hdr', 'lat@hdr', 'lon@hdr', 'vertco_reference_1@body', 'obsvalue@body', 'varno@body' , 'units',  'number_of_pressure_levels' ]
column_names_igra2 = [ 'source_id', 'report_id',  'observation_id', 'record_timestamp' , 'iday', 'statid@hdr', 'lat@hdr', 'lon@hdr', 'vertco_reference_1@body', 'obsvalue@body', 'varno@body' , 'units',  'number_of_pressure_levels', 'report_timestamp' ] 
# record_timestamp: same as date_time in the observations_table
# report_timestamp: only for igra2, representes the release time of the sonde 


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
    bufr_values = []
    
    """ Name of the columns as they will appear in the pandas dataframe (not necessarily CDM compliant) """
    #column_names = ['report_timestamp' , 'iday',  'station_id', 'latitude', 'longitude', 'pressure', 'value','varno@body']
       
    lat, lon, alt, blockNumber, stationNumber, statid = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    
    obs_id, report_id = -1, 0 # progressive observation id
    
    while 1:
        #lista = [] # temporary list
        bufr = codes_bufr_new_from_file(f)
   
        if bufr is None:
            break
   
        codes_set(bufr, 'unpack', 1) # eCcodes must expand all the descriptors and unpack the data section
    
        date = '19'+codes_get_array(bufr, "typicalDate")[0][2:]
        timePeriod = codes_get_array(bufr, "typicalTime")[0]   
        
        year, month, day =  date[0:4], date[4:6] , date[6:8]
        hour, minutes = timePeriod[0:2] , timePeriod[2:4]
                       
        idate =  datetime.strptime(year + month + day + hour + minutes, '%Y%m%d%H%M')
        iday = int(year + month + day )

        pressure             = codes_get_array(bufr, "pressure") 
        temperature       = codes_get_array(bufr, "airTemperature")           
        wind_direction    = codes_get_array(bufr, "windDirection")
        wind_speed        = codes_get_array(bufr, "windSpeed")
        
        try:  # not all the bufr files have the dewpoint 
            dew_point          = codes_get_array(bufr, "dewpointTemperature")
        except:
            dew_point= np.empty((1, len(temperature)))
            dew_point[:] = np.nan
            
        num_lev             = len(pressure) # number of  distinct pressure levels 
        
        try:
            geopotential   = codes_get_array(bufr, "nonCoordinateGeopotentialHeight")         
        except:
            geopotential = np.full( (1,len(temperature)) , np.nan )[0,:]
                
        if report_id == 0:
            ''' Check again but these values should remain the same for all cnt, so it makes no sense to read them every time '''
            lat                     = codes_get(bufr, "latitude")
            lon                    = codes_get(bufr, "longitude")
            alt                     = float(codes_get(bufr, "heightOfStation"))
            blockNumber    = codes_get(bufr, "blockNumber")
            stationNumber = codes_get(bufr, "stationNumber")
            statid                = str(blockNumber*1000+stationNumber)
            
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
            if dp == miss_value:
                dp = np.nan
            if airT == miss_value :    # replacing none values with numpy nans
                airT = np.nan 
            if winds == miss_value:
                winds = np.nan
            if windd == 2147483647:
                windd = np.nan 
                
            for value,var in zip( [gph, airT, winds, windd, dp],  ['gph', 'temperature', 'wind_speed', 'wind_direction', 'dew_point'] ):
                obs_id = obs_id + 1 
                bufr_values.append( ( 'BUFR', report_id, obs_id,  idate, iday, statid, lat, lon, press, value, cdmvar_dic[var]['odb_var'] , cdmvar_dic[var]['cdm_unit'], num_lev  ) ) 
        
        report_id += 1
            
    #column_names = ['source_file', 'product_code', 'report_id',  'observation_id', 'report_timestamp' , 'iday', 'station_id', 'lat@hdr', 'lon@hdr', 'vertco_reference_1@body', 'obsvalue@body', 'varno@body' , 'units',  'number_of_pressure_levels' ]
    df = pd.DataFrame(data= bufr_values, columns= column_names)
    df['vertco_type@body'] = 1        
    df.sort_values(by = ['record_timestamp', 'vertco_reference_1@body' ] )    
    df['report_id'] = numpy.int64 (df['report_id'] ) 
    df['observation_id'] = numpy.int64 (df['observation_id'] )     
    return df 
    #return df.to_xarray()
    
    

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
    
    #source_file = [l for l in file.split('/') if '.txt' in l][0]

    nmiss = 0
    search_h = False   
    read_data = []
    
    usi,idate, usi, lat, lon, lat, stype, press, gph, temp, rh, wdir, wspd = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

    obs_id = 0
    
    for i, line in enumerate(data):
        if line[0] == 'H':
            try:
                # Header
                usi      = int(line[2:14])  # unique station identifier
                ident    = line[15:21].replace(' ','')# WMO
                if len(ident) == 4:
                    ident = '0' + ident 
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
                print("Error: ", i, line, repr(e), "Skipping Block:")
                search_h = True
                #iprev = i

        elif search_h:
            nmiss += 1
            continue  # Skipping block

        else:
            # Data
            #ltyp      = int(line[0:4])
            p   = float(line[5:13])
            
            if p != -99999.0 and p != 9999.9: 
                press   = float(line[5:13])*100  # converting to Pa, since P is given in mb (1 mb = 1 hPa) 
            else:
                press = np.nan                 
                
            gph = float(line[14:22]) # gph [m]
            
            ###if gph == -99999.0 or gph == -99999.00 or gph >= 99999.0:
            ###    gph = np.nan
         
            temp = float(line[23:29])
            if temp == -999.0:
                temp = np.nan 
            else:
                temp = temp + 273.15
                
            rh  = float(line[30:36])  # %
            if rh == -999.0:
                rh = np.nan
            else:
                rh = rh / 100.  # convert to absolute ratio  TODO

            wdir    = float(line[37:43])            
            if wdir == -999.0 or wdir == -999 :
                wdir = np.nan
            
            wspd   = float(line[44:50])  # [m/s], module of the velocity
            if wspd <0 :
                wspd = np.nan             
                
            for value,var in zip([gph, temp, wspd, wdir, rh],  ['gph', 'temperature', 'wind_speed', 'wind_direction', 'relative_humidity'] ):
                obs_id = obs_id +1
                read_data.append( ( 'NCAR', usi, obs_id, idate, iday, ident, lat, lon, press, value, cdmvar_dic[var]['odb_var'] , cdmvar_dic[var]['cdm_unit'], numlev) )
              
              
    
    #column_names = ['source_file', 'product_code', 'report_id', 'observation_id', 'report_timestamp' , 'iday', 'station_id', 'lat@hdr', 'lon@hdr', 'vertco_reference_1@body', 'obsvalue@body', 'varno@body' ,  'units',  'number_of_pressure_levels' ]
    
    df = pd.DataFrame(data= read_data, columns= column_names)       
    df = df.replace([-999.9, -9999, -999, -999.0, -99999.0, -99999.9, 99999.0,  -99999.00 ], np.nan)
    
        
    df['vertco_type@body'] = 1
    df.sort_values(by = ['record_timestamp', 'vertco_reference_1@body' ] )    
    df['report_id'] = numpy.int64 (df['report_id'] ) 
    df['observation_id'] = numpy.int64 (df['observation_id'] ) 
    
    #return df.to_xarray()
    return df


def igra2_ascii_to_dataframe(file=''):
    """ Read an igra2 stationfile in ASCII format and convert to a Pandas DataFrame. 
        Adapted from https://github.com/MBlaschek/CEUAS/tree/master/CEUAS/data/igra/read.py 
        Variables used inside the DataFrame are already CDM compliant
        
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
    idate = np.nan
    count = 0
    head_count = 0
    
    obs_id = 0
    for i, line in enumerate(data):
        if line[0] == '#':
            head_count = head_count +1 
            # Info from the Header line of each ascent                                                                                                                                                                                                                   
            ident     = line[1:12]               # station identifier
            ident     = ident[6:12]
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
            iday =  int(year + month + day)
            count = count + 1
        else:
           # Data of each ascent
            lvltyp1 = int(line[0])            # 1-  1   integer major level type indicator
            lvltyp2 = int(line[1])            # 2-  2   integer minor level type indicator
            etime   = int(line[3:8])          # 4-  8   integer elapsed time since launch
            press   = int(line[9:15])         # 10- 15  integer reported pressure
            pflag   = line[15]                # 16- 16  character pressure processing flag
            
            gph     = int(line[16:21])        # 17- 21  integer geopotential height  [m]
            
            if gph == -9999 or gph == -8888:   # reading the values andh check if they are missing or removed as -9999 or -8888 before dividing by 10 as the instructions say 
                gph = np.nan # 23- 27  integer temperature, [Celsius to Kelvin ]    
                
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
                
            for value,var in zip([gph, temp, wspd, wdir, rh, dpdp],  ['gph', 'temperature', 'wind_speed', 'wind_direction', 'relative_humidity' , 'dew_point_depression'] ):
                obs_id = obs_id +1 
                #read_data.append( (source_file, 'NCAR', usi, obs_id,                idate, iday, ident, lat, lon, press, value, cdmvar_dic[var]['odb_var'] , cdmvar_dic[var]['cdm_unit'], numlev) )
                read_data.append ( ( 'IGRA', head_count,  obs_id,  idate, iday, ident, lat, lon, press, value, cdmvar_dic[var]['odb_var'], cdmvar_dic[var]['cdm_unit'], numlev, reltime ) )
                #print('check', value, cdmvar_dic[var] , var )                   
            #column_names = ['product_code', 'report_id', 'observation_id', 'report_timestamp' , 'iday', 'station_id', 'lat@hdr', 'lon@hdr', 'vertco_reference_1@body', 'obsvalue@body', 'varno@body', 'number_of_pressure_levels' , 'units']
            
    df = pd.DataFrame(data= read_data, columns= column_names_igra2)
    df['vertco_type@body'] = 1    
    df.sort_values(by = ['record_timestamp', 'vertco_reference_1@body' ] )    
    df['report_id'] = numpy.int64 (df['report_id'] ) 
    df['observation_id'] = numpy.int64 (df['observation_id'] ) 
    
    return df


def read_all_odbsql_stn_withfeedback(odbfile):

    #alldata=''

    alldict={} #xr.Dataset()
    t=time.time()
    #sonde_type=True
    #obstype=True
    ol=odbfile.split('.')
    if ol[-1][0]!='_':
        odbfile='.'.join(ol[:-1])+'._'+ol[-1]
    if os.path.getsize(odbfile+'.gz')>0:
        """ Read first the odbb header to extract the column names and type """
        try:
            try:                
                with open(os.path.dirname(odbfile)+'/odbheader','rb') as f:
                    rdata=f.read()
                    print('read the odb file' , odbfile)
            except:
#                rdata=subprocess.check_output(["odb","header",'/3645/'.join(odbfile.split('/1/'))])
                #rdata=subprocess.check_output(["odb","header",'.'.join(odbfile.split('._'))])
                rdata=subprocess.check_output(["odb","header",'.'.join(odbfile.split('._'))])
                
                with open('odbheader','wb') as f:
                    f.write(rdata)
            rdata=rdata.decode('latin-1').split('\n')
            columns=[]
            kinds=[]
            tdict={}
            for r in rdata[2:-2]:
                try:
                    #print(r[:6])
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
                                tdict[columns[-1]]=numpy.int32
                        else:
                            tdict[columns[-1]]=numpy.dtype('S') # dict containng column name and type
                                     
                            
                except IndexError:
                    pass
        except KeyError:
            print('could not read odbfile '+odbfile)
            return alldict
        try:
            # returns a byte string
            t=time.time()
            
            try:
                f=gzip.open(odbfile+'.gz')                 
            except:
                print(odbfile+'.gz','not found')
                return
                
            # access the string like a file, return a file pointer to read the string with pandas
            # nb  if you have null values, reading of integer fails and are read as floats
            # to improve, you can convert the columns with nabs in the alldicts (pandas data frame) into int(np.nans)
            # date values are large so float 32 precision is not sufficient 

            #d=['date@hdr','time@hdr','statid@hdr','vertco_reference_1@body','varno@body','reportype','andate','antime',
            #                 'obsvalue@body','fg_depar@body','an_depar@body','biascorr@body','sonde_type@conv','collection_identifier@conv','source@hdr']
            
            # had to remove 'collection_identifier@conv'
            d=['date@hdr','time@hdr','statid@hdr','vertco_reference_1@body','varno@body','lon@hdr','lat@hdr','seqno@hdr',
                             'obsvalue@body','source@hdr' , 'vertco_type@body']
            if 'fg_depar@body' in columns:
                d=d+['fg_depar@body','an_depar@body','biascorr@body','sonde_type@conv','reportype','andate','antime']
            
            
            for c in columns:
                if c not in d:
                    del tdict[c]
                    
            columns=d.copy()
            
            alldict=pd.read_csv(f,delimiter='\t',usecols=columns,quoting=3,comment='#', skipinitialspace=True,dtype=tdict)#,nrows=1000000)
            
            """ Case where erafb is not available """
            if 'fg_depar@body' not in columns:
                alldict['fg_depar@body']=numpy.float32(numpy.NaN)
                alldict['an_depar@body']=numpy.float32(numpy.NaN)
                alldict['biascorr@body']=numpy.float32(numpy.NaN)
                alldict['sondetype@conv']=numpy.int32(-2147483648)
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
    
                print(time.time()-t,sys.getsizeof(alldict)//1024//1024)
                
                idx=numpy.where(alldict.reportype.values==16045)[0]
                if idx.shape[0]>0:
                    idy=numpy.where(numpy.logical_and(alldict.reportype.values!=16045,alldict.reportype.values!=16068))[0]
                    if idy.shape[0]>0:
                        idz=numpy.isin(alldict.andate.values[idy],alldict.andate.values[idx])
                        if numpy.sum(idz)>0:
                            alldict.drop(index=alldict.index[idy[idz]],inplace=True)
                           
                idx=numpy.where(alldict.reportype.values==16068)[0]
                if idx.shape[0]>0:
                    idy=numpy.where(numpy.logical_and(alldict.reportype.values!=16045,alldict.reportype.values!=16068))[0]
                    if idy.shape[0]>0:
                        idz=numpy.isin(alldict.andate.values[idy],alldict.andate.values[idx])
                        if numpy.sum(idz)>0:
                            alldict.drop(index=alldict.index[idy[idz]],inplace=True)
                          

            print(time.time()-t,sys.getsizeof(alldict)//1024//1024)

            for c in alldict.columns:
                
                if type(alldict[c].iloc[0]) in [str,bytes]:
                    l=alldict[c].shape[0]
                    slen=len(alldict[c].values[0])
                    #alldict[c]=numpy.array(alldict.pop(c).values,dtype='S{}'.format(slen))
                    alldict[c]=numpy.string_(alldict[c])
                    
                if type(alldict[c].iloc[0]) is numpy.int64:
                    alldict[c]=numpy.int32(alldict[c])
                    
                if type(alldict[c].iloc[0]) is numpy.float64:
                    alldict[c]=numpy.float32(alldict[c])

            print('after odb:',time.time()-t)
   
        except subprocess.CalledProcessError as e:
            print('odb failed!:'+' '+odbfile)
            return alldict

    print(odbfile,time.time()-t)
    print(odbfile,time.time()-t,sys.getsizeof(alldict))

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
        x=fbv[cdmfb][di['recordindex'].values]
        
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


def make_datetime_seconds_indices(date_time):
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

def write_dict_h5(dfile, f, k, fbencodings, var_selection=[], mode='a', attrs={}): # cuts vars and copies attributes of observation, feedback and header tables
    """ Writes each separate variable from the observation or feedback tables inot netcdf using h5py.
          f is a pandas dataframe with one column, one for each variable
          k is either 'era5fb' or 'observations_table'
          fbencodings is the encodings of variable types, e.g. {'observations_id': { 'compression': 'gzip' } }
    """

    print('writing write_dict_h5 for ' , k  )
    #attrs=  {'date_time':('units','seconds since 1900-01-01 00:00:00')}
    attrs = {'observation_id': ('description', 'unique ID for observation'), 'report_id': ('description', 'Link to header information') , 'date_time':('units','seconds since 1900-01-01 00:00:00') }
    
    with h5py.File(dfile,mode) as fd:
        try:
            fd.create_group(k)
            index=numpy.zeros(f[f.columns[0]].shape[0],dtype='S1')
            fd[k].create_dataset('index', data=index)
        except:
            pass
        if not var_selection:
            var_selection=list(f.keys())
        
        string10=numpy.zeros(80,dtype='S1')
        sdict={}
        slist=[]

         
        for v in var_selection:
            if type(f[v].values[0]) not in [str,bytes,numpy.bytes_]:
                if f[v].values.dtype !='S1':
                    
                    fd[k].create_dataset(v,f[v].values.shape,f[v].values.dtype,compression=fbencodings[v]['compression'],chunks=True)
                    fd[k][v][:]=f[v]
                    if attrs:    #  attrs={'date_time':('units','seconds since 1900-01-01 00:00:00')}
                        if v in attrs.keys():
                            fd[k][v].attrs[attrs[v][0]]=numpy.bytes_(attrs[v][1])
                            print (  fk, ' ' , v , ' ' ,   ) 
                else:
                    fd[k].create_dataset(v,f[v].values.shape,f[v].values.dtype,compression=fbencodings[v]['compression'],chunks=True)
                    fd[k][v][:]=f[v][:]
            
            else:
                sleno=len(f[v].values[0])
                slen=sleno
                x=numpy.array(f[v].values,dtype='S').view('S1')
                slen=x.shape[0]//f[v].values.shape[0]
                sdict[v]=slen
                if slen not in slist:
                    slist.append(slen)
                    try:
                        fd[k].create_dataset('string{}'.format(slen),data=string10[:slen])
                    except:
                        pass
                #if attrs:    #  attrs={'date_time':('units','seconds since 1900-01-01 00:00:00')}
                #    if v in attrs.keys():
                #        print  ('v is: !!! ' , v ) 
                #        fd[k][v].attrs[attrs[v][0]]=numpy.bytes_(attrs[v][1])
                            
                x=x.reshape(f[v].values.shape[0],slen)
                fd[k].create_dataset(v,data=x,compression=fbencodings[v]['compression'],chunks=True)


                
        for v in fd[k].keys(): #var_selection:
            l=0
            
            try:
                if 'string' not in v and v!='index':
                    
                    fd[k][v].dims[l].attach_scale(fd[k]['index'])
                    if type(f[v].values[0]) in [str,bytes,numpy.bytes_]:
                        slen=sdict[v]
                        #slen=10
                        fd[k][v].dims[1].attach_scale(fd[k]['string{}'.format(slen)])
                        #print('')
            except:
                pass
            
            
        i=4        
        for v in slist:
            s='string{}'.format(v)
            for a in ['NAME']:
                fd[k][s].attrs[a]=numpy.bytes_('This is a netCDF dimension but not a netCDF variable.')
            
            i+=1
    return

        


                
def initialize_output(fn, output_dir, station_id, dataset):
    """ Simple initializer for writing the output netCDF file """
    
    source_file =fn.split('/')[-1]
    output_file = output_dir + '/' + station_id + '_' + dataset + '_harvested_' + source_file + '.nc'  # creating an output file name e.g. chera5.conv._10393.nc  , try 01009 faster
    return output_file , source_file 



def convert_variable_type(df):
    """ Converts the variable type of the input DataFrame  (igra2,ncar,bufr) """
    ''' # available columns
    'source_file', 'source_id', 'report_id', 'observation_id',
       'record_timestamp', 'iday', 'station_id', 'lat@hdr', 'lon@hdr',
       'vertco_reference_1@body', 'obsvalue@body', 'varno@body', 'units',
       'number_of_pressure_levels'
    '''  
    dic_var_type = { 'int32'     : ['varno@body', 'number_of_pressure_levels']  ,
                                'float32' :  ['lat@hdr', 'lon@hdr' ,  'vertco_reference_1@body', 'obsvalue@body',  'iday' ]  ,
                                'string'   :  ['source_id' , 'station_id' ,  'source_file' , 'report_id', 'observation_id',   'units',  ]  }
    
    for c in df.columns:
        if c == 'record_timestamp':
            continue 
        if c == 'observation_id' or c == 'report_id':
            df[c] = numpy.string_ (df[c] ) 
            continue 
        if c in dic_var_type['int32']:
            df[c] = numpy.int32 ( df[c] )
        elif c in dic_var_type['string']:
            df[c] = numpy.string_( df[c] )
        else:
            df[c] = numpy.float32( df[c] )
            
    return df

 
def datetime_toseconds(date_time):
    """ Converts a generic date_time array to seconds since '1900-01-01 00:00:00' """
    offset = np.datetime64('1900-01-01 00:00:00')            
    deltas =  [ dt - offset for dt in date_time ]
    date_times_seconds  =  [ i.total_seconds() for i in deltas ] 
    return date_times_seconds # replacing with seconds from 1900-01-01 00:00:00 
    
    
    
    
    
def df_to_cdm(cdm, cdmd, out_dir, dataset, fn):
    """ Convert the  pandas dataframe from the file into cdm compliant netCDF files. Use with bufr, igra2 and ncar databases.
        According to each file type, it will use the appropriate reading function to extract a Pandas DataFrame
        input:
              fn       :: odb file name (e.g. era5.conv._10393)
              cdm   :: cdm tables (read with pandas)
              cdmd :: cdm tables definitions ("")  """    
        
    #station_id_fails = open('station_id_fail.log' , 'a') 
    #station_id_ok = open('station_id_correct.log' , 'a')
    
    t=time.time()
    if not False:        
            # era5 analysis feedback is read from compressed netcdf files era5.conv._?????.nc.gz in $RSCRATCH/era5/odbs/1
            """ Reading the odb and convert to xarray """  
            if  'bufr' in dataset :
                df= bufr_to_dataframe(fn) # fdbs: the xarray converted from the pandas dataframe 
                
            elif  'ncar' in dataset:
                df = uadb_ascii_to_dataframe(fn)
                
            elif  'igra2' in dataset:
                df = igra2_ascii_to_dataframe(fn)
                
            else:
                print('Unidentified file is: ', fn)
                raise ValueError('Cannot identify the type of file to be analized!!! ')
              
              
            station_id = str( df['statid@hdr'].values[0].replace(' ','') )                          
            station_configuration_retrieved = get_station_configuration( station_id, cdm['station_configuration'] )    
            primary_id = station_configuration_retrieved['primary_id'].values[0].decode('latin1')      
            fno,  source_file = initialize_output(fn, output_dir, primary_id, dataset)   
            
                             
            """ Casting the original variable types to appropriate numpy types """     
            df = convert_variable_type(df)  # ->think this is overly complicated. 
                          
            """ Storing the variable encodings """
            fbencodings={}
            for d,v in df.items():
                
                if v.dtype==numpy.dtype('float64'):
                        fbencodings[d]={'dtype':numpy.dtype('float32'),'compression': 'gzip'}           
                        
                elif v.dtype==numpy.dtype('int32'):
                    fbencodings[d]={'dtype':numpy.dtype('int32'),'compression': 'gzip'} 
                    
                elif type(v.values[0])==bytes:
                        fbencodings[d]={'compression': 'gzip', 'chunksizes': ( min( [10000,v.shape[0] ] ), 10 ) }#,'chunksizes':(10000,10)
                else:
                        fbencodings[d]={'compression': 'gzip'}
            fbencodings['index']={'compression': 'gzip'}

    
            """ Extract the unique indices of each date observation, one for only dates, one for date_time (i.e. record index). 
                Converts the time variables in seconds since 1900-01-01 00:00:00 """  
            di=xr.Dataset() 
            
            

            record_timestamp_seconds = datetime_toseconds( df['record_timestamp'] )  # will fill the header_table 
            df['record_timestamp'] = record_timestamp_seconds # replacing with seconds from 1900-01-01 00:00:00 

            try:               
                report_timestamp_seconds = datetime_toseconds( df['report_timestamp'] )  # will fill the header_table 
                df['report_timestamp'] = report_timestamp_seconds # replacing with seconds from 1900-01-01 00:00:00 
            except KeyError:
               print('Skipping report_timestamp, only available for IGRA2')
               
            date_time_seconds = datetime_toseconds( df['date_time'] ) # will fill the observatiosn_table 
            df['date_time'] = date_time_seconds # replacing with seconds from 1900-01-01 00:00:00 
           
           

           
            indices, date_times, counts = make_datetime_seconds_indices( df['iday'].values )   #only date information
            di['dateindex']  = ( { 'dateindex' :  date_times.shape } , date_times )          
            
            indices, date_times , counts  = make_datetime_seconds_indices( df['record_timestamp'].values ) #date_time plus indices           
            di['recordindex']          = ( {'recordindex' : indices.shape }, indices )
            di['recordtimestamp']  = ( {'recordtimestamp' : date_times.shape }, date_times  )
            
            di['recordtimestamp'].attrs['units']='seconds since 1900-01-01 00:00:00'

            di.to_netcdf( fno, format='netCDF4', engine='h5netcdf', mode='w' )
         
            #write_dict_h5(fno, df, 'era5fb', fbencodings, var_selection=[],mode='a')
            dcols=[]
            for d in df.columns:
                    if d not in ['date@hdr','time@hdr','statid@hdr','vertco_reference_1@body','varno@body', 'lon@hdr','lat@hdr','seqno@hdr',
                                       'obsvalue@body','fg_depar@body','an_depar@body','biascorr@body','sonde_type@conv',  'record_timestamp' ,
                                       'observation_id', 'report_id' ] :
                         
                        dcols.append(d)
            df.drop(columns=dcols,inplace=True)
          
            groups={}
            groupencodings={}
            for k in cdmd.keys(): # loop over all the table definitions 
                if k in ('observations_table'):
                    pass #groups[k]=pd.DataFrame()
                else:
                    groups[k]=xr.Dataset() # create an  xarray
                groupencodings[k]={} # create a dict of group econding
          
                for i in range(len(cdmd[k])):
                    d=cdmd[k].iloc[i] 
                    
                    if k in ('observations_table'):
                        groups[k]=pd.DataFrame()
                        try:

                            if d.element_name=='report_id':                            
                                groups[k][d.element_name]=fromfb_l(df, di._variables, cdmfb[k+'.'+d.element_name], ttrans(d.kind, kinds=okinds))
                            else:
                                groups[k][d.element_name]=fromfb_l(df,di._variables,cdmfb[d.element_name], ttrans(d.kind,kinds=okinds) )
                            print('no key error')
                            #groups[k][d.element_name]=fromfb_l(df,di._variables,cdmfb_noodb[d.element_name], ttrans(d.kind, kinds=okinds) )
                            
                        except KeyError:
                            print('key error!! )')
                            x=numpy.zeros( df['record_timestamp'].shape[0], dtype=numpy.dtype(ttrans(d.kind,kinds=okinds) ) )
                            x.fill(numpy.nan)
                            groups[k][d.element_name]=x
                                  
                    elif k in ('header_table'):
                              # if the element_name is found in the cdmfb dict, then it copies the data from the odb into the header_table
                        try:
                            
                            if d.element_name=='report_id':
                                groups[k][d.element_name]=({'hdrlen':di['recordindex'].shape[0]}, hdrfromfb(df, di._variables, cdmfb[k+'.'+d.element_name], ttrans(d.kind,kinds=gkinds) ) )
                            else:
                                groups[k][d.element_name]=({'hdrlen':di['recordindex'].shape[0]}, hdrfromfb(df, di._variables, cdmfb[d.element_name], ttrans(d.kind,kinds=gkinds) ) )
                            j=0
                            
                            
                        except KeyError:                            
                            if d.element_name in cdm['station_configuration'].columns:
                                x=numpy.zeros(di['recordindex'].shape[0],dtype=numpy.dtype(ttrans(d.kind,kinds=gkinds)))
                                try:
                                    idx=numpy.where('0-20000-0-'+ station_id == cdm['station_configuration']['primary_id'])[0][0]
                                    groups[k][d.element_name]=x.fill(cdm['station_configuration'][d.element_name][idx])
                                except:
                                    groups[k][d.element_name]=x
                            else:
                                x=numpy.zeros(di['recordindex'].shape[0],dtype=numpy.dtype(ttrans(d.kind,kinds=okinds)))
                                x.fill(numpy.nan)
                            groups[k][d.element_name]=({'hdrlen':di['recordindex'].shape[0]},x)
                            
 
                    elif k in ('station_configuration'): # station_configurationt contains info of all the stations, so this extracts only the one line for the wanted station with the numpy.where
                        try:                        
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
                        if type(groups[k][d.element_name].values[0])==str:
                            s=groups[k][d.element_name].values.shape
                            groupencodings[k][d.element_name]={'dtype':numpy.dtype('S80'),'compression': 'gzip','chunksizes':(min(100000,s[0]),80)}
                        else:
                            groupencodings[k][d.element_name]={'compression': 'gzip'}
                        
                        if k in ('observations_table'):
                            write_dict_h5(fno, groups[k], k, groupencodings[k], var_selection=[],mode='a', attrs={'date_time':('units','seconds since 1900-01-01 00:00:00')})
                            print('written observations_table' , k )
                            
                    except:
                        #print('bad:',k,d.element_name)
                        pass
                
            for k in groups.keys():            
                if k not in ('observations_table') :   
                    print('Doing: ', k )
                    groups[k].to_netcdf(fno,format='netCDF4',engine='h5netcdf',encoding=groupencodings[k],group=k,mode='a') #
            del df
          
    return 0
          
          
          

def get_station_configuration(station_id, station_configuration):
    """ Gets the primary station_id from the station_configuration table. 
          station_id is the id taken from the input file.
          First it checks if a primary_id in th estation_conf file matches the station_id, 
          otherwise it looks for an element in the list of secondary ids.         
          """
    si = station_id.decode('latin1')
    
    if ':' in si:  # primary ids can have either the 20000 or 20001 numerical flag 
        station_id_primary = numpy.string_( '0-20000-0-' + si.split(':')[1] )   # remove the prefix to the station id 
        station_id_primary_alternative = numpy.string_( '0-20001-0-' + si.split(':')[1] )
    else:
        station_id_primary = numpy.string_('0-20000-0-' + si )
        station_id_primary_alternative = numpy.string_( '0-20001-0-' + si )
        
            
    #station_id = numpy.string_( '1:68351' )
    
    """ First, check for matching primary_id. 
          If not found, check for secondary id. Note that secondary is a list, so must loop over the entry to find a matching one """
    
    matching_primary = station_configuration.loc[station_configuration['primary_id'] == station_id_primary ]
    matching_primary_alt= station_configuration.loc[station_configuration['primary_id'] == station_id_primary_alternative ]
    
    if len(matching_primary) > 0:
        return matching_primary 
    
    elif   len(matching_primary_alt) > 0  :
        return matching_primary_alt 

    else:
        secondary = station_configuration['secondary_id'] 
        for s in secondary:
            if s == station_id:
                print('Found element in secondary ids')
                sc = station_configuration.loc[station_configuration['secondary_id'] == s ]
                return sc 
        return 0        
    

def odb_to_cdm(cdm, cdmd, output_dir, dataset, obs_tab_attrs, fn):
    """ Convert the  file into cdm compliant netCDF files. 
        According to each file type, it will use the appropriate reading function to extract a Pandas DataFrame
        input:
              fn       :: odb file name (e.g. era5.conv._10393)
              cdm   :: cdm tables (read with pandas)
              cdmd :: cdm tables definitions ("")  """

    process = psutil.Process(os.getpid())
    t=time.time()
    #fno = initialize_convertion(fn, output_dir) 
    #station_id = ''    
    
    # restored 
    fnl=fn.split('/')
    fnl[-1]='ch'+fnl[-1]
        
    if not False:
        
        # era5 analysis feedback is read from compressed netcdf files era5.conv._?????.nc.gz in $RSCRATCH/era5/odbs/1
        
        fbds=read_all_odbsql_stn_withfeedback(fn) 
        
        """ Read the station_id, getting the station_configuration from the table list, extracting primary_id """
       
        #station_id = numpy.string_(fbds['statid@hdr'][0][1:-1].decode('latin1') )    
        station_id = fbds['statid@hdr'][0][1:-1]    
        station_configuration_retrieved = get_station_configuration( station_id, cdm['station_configuration'] )            
        primary_id = station_configuration_retrieved['primary_id'].values[0].decode('latin1')  
        fno,  source_file = initialize_output(fn, output_dir, primary_id, dataset)        
        
        #fbds['source_file'] = source_file
        if fbds is None:
            return
        #fbds=xr.open_dataset(f)
        #print(time.time()-t) # to check the reading of the odb
        # the fbencodings dictionary specifies how fbds will be written to disk in the CDM compliant netCDF file.
        # float64 is often not necessary. Also int64 is often not necessary. 
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
        
        y=numpy.int64(fbds['date@hdr'].values)*1000000+fbds['time@hdr'].values
        tt=time.time()
        idx=numpy.lexsort((fbds['vertco_reference_1@body'].values,y))
        y=y[idx]
        for fb in fbds.keys():
            fbds[fb].values[:]=fbds[fb].values[idx]
        #print(time.time()-tt)
        x=numpy.unique(y)
        z=find_recordindex_l(y,x)
        di=xr.Dataset() 
        di['recordindex']=({'record':z.shape[1]},z[1])
        #x=make_datetime(x//1000000,x%1000000)
        df=pd.DataFrame({'year':x//10000000000,'month':(x%10000000000)//100000000,'day':x%100000000//1000000,
                            'hour':x%1000000//10000,'minute':(x%10000)//100,'second':x%100})
        dt=pd.to_datetime(df).values
        di['recordtimestamp']=({'record':z.shape[1]},numpy.array(dt-numpy.datetime64('1900-01-01'),dtype=int)//1000000000)
        di['recordtimestamp'].attrs['units']='seconds since 1900-01-01 00:00:00'
        del dt,df
        y=fbds['date@hdr'].values
        x=numpy.unique(y)
        z=find_dateindex_l(y,x)
        di['dateindex']=({'days':z.shape[1],'drange':z.shape[0]},z) # date, index of the first occurrance, index of the last
        del y

        #this writes the dateindex to the netcdf file. For faster access it is written into the root group
        di.to_netcdf(fno,format='netCDF4',engine='h5netcdf',mode='w')
          
        write_dict_h5(fno, fbds, 'era5fb', fbencodings, var_selection=[],mode='a')
        dcols=[]
        for d in fbds.columns:
            if d not in ['date@hdr','time@hdr','statid@hdr','vertco_reference_1@body','varno@body', 'lon@hdr','lat@hdr','seqno@hdr',
                         'obsvalue@body','fg_depar@body','an_depar@body','biascorr@body','sonde_type@conv', 'vertco_type@body']:
                dcols.append(d)
        fbds.drop(columns=dcols,inplace=True)

        #print(sys.getsizeof(fbds)//1024//1024,process.memory_info().rss//1024//1024)        
        #print(time.time()-t)

        # each cdm table is written into an hdf group, groups is the dict of all the groups
        # to write the group to the disk, you need the group encoding dict
        groups={}
        groupencodings={}
        for k in cdmd.keys(): # loop over all the table definitions 
            if k in ('observations_table'):
                pass #groups[k]=pd.DataFrame()
            else:
                groups[k]=xr.Dataset() # create an  xarray
            groupencodings[k]={} # create a dict of group econding

            for i in range(len(cdmd[k])): # in the cdm table definitions you always have the element(column) name, the type, the external table and the description 
                d=cdmd[k].iloc[i] # so here loop over all the rows of the table definition . iloc is just the index of the element
                # two possibilities: 1. the corresponding  table already exists in the cdm (case in the final else)
                #                    2. the corr table needs to be created from the local data sources (e.g. the feedback or IGRA or smt else). 
                # These are the observation_tables, the header_tables and the station_configuration.
                # These tables are contained in the CEUAS GitHub but not in the cdm GitHub
                if k in ('observations_table'):
                    groups[k]=pd.DataFrame()
                    try:
                        # fbds is an xarray dataset , fbds._variables is a dict of the variables 
                        if d.element_name=='report_id':                            
                            groups[k][d.element_name]=fromfb_l(fbds,di._variables,cdmfb[k+'.'+d.element_name], ttrans(d.kind,kinds=okinds))
                        else:
                            groups[k][d.element_name]=fromfb_l(fbds,di._variables,cdmfb[d.element_name], ttrans(d.kind,kinds=okinds))
                    except KeyError:
                        x=numpy.zeros(  fbds['date@hdr'].shape[0],dtype=numpy.dtype(ttrans(d.kind,kinds=okinds) ) )
                        x.fill(numpy.nan)
                        groups[k][d.element_name]=x
                        
                elif k in ('header_table'):
                    print(d.element_name)
                    try:
                        if d.element_name=='report_id':
                            groups[k][d.element_name]=( {'hdrlen':di['recordindex'].shape[0]}, hdrfromfb(fbds,di._variables, cdmfb[k+'.'+d.element_name],ttrans(d.kind,kinds=gkinds) ) )
                        else:
                            groups[k][d.element_name]=( {'hdrlen':di['recordindex'].shape[0]}, hdrfromfb(fbds,di._variables, cdmfb[d.element_name],ttrans(d.kind,kinds=gkinds) ) )
                        j=0
                            
                    except KeyError:                            
                        if d.element_name in cdm['station_configuration'].columns:
                            x=numpy.zeros(di['recordindex'].shape[0],dtype=numpy.dtype(ttrans(d.kind,kinds=gkinds)))
                            try:
                                
                                idx=numpy.where('0-20000-0-'+fnl[-1].split('_')[-1] == cdm['station_configuration']['primary_id'])[0][0]
                                groups[k][d.element_name]=x.fill(cdm['station_configuration'][d.element_name][idx])
                            except:
                                groups[k][d.element_name]= ({'hdrlen':di['recordindex'].shape[0]},x)
                        else:
                            x=numpy.zeros(di['recordindex'].shape[0],dtype=numpy.dtype(ttrans(d.kind,kinds=okinds)))
                            x.fill(numpy.nan)
                        groups[k][d.element_name]=({'hdrlen':di['recordindex'].shape[0]},x)
               
                        
                #elif k in ('station_configuration'): # station_configuration contains info of all the stations, so this extracts only the one line for the wanted station with the numpy.where
                    
                    #try:   
                        #if 'sci' not in locals(): 
                            #sci=numpy.where( cdm[k]['primary_id']== numpy.string_('0-20000-0-'+fbds['statid@hdr'][0][1:-1].decode('latin1')  )   )[0]
                        #if len(sci)>0:
                            #groups[k][d.element_name]=({k+'_len':1}, cdm[k][d.element_name].values[sci])
                            #fno_new = fno.replace('STATIONID' , station_id)
                            
                            
                            
                        #else:                                                                               
                            #if dataset != 'era5_3188':
                                #sci=numpy.where(cdm[k]['secondary_id']== numpy.string_(fbds['statid@hdr'][0][1:-1].decode('latin1') )   )[0]                                
                            #elif dataset == 'era5_3188':
                                ##if len ( str (numpy.string_(fbds['statid@hdr'][0].decode('latin1')) ) ):                                    
                                #sci=numpy.where( cdm[k]['secondary_id']== numpy.string_(fbds['statid@hdr'][0][3:-1].decode('latin1') )   )[0]
                                
                            #if len(sci)>0:
                                #groups[k][d.element_name]=({k+'_len':1}, cdm[k][d.element_name].values[sci])
                                ##fno_new = fno.replace('STATIONID' , cdm[k]['secondary_id'].split('-')[-1] )
                            
                            ##print('statconf:',k,groups[k][d.element_name])
                            
                    #except KeyError:
                        ##print('x')
                        #pass
                
                
                
                elif k in ('station_configuration'): # station_configurationt contains info of all the stations, so this extracts only the one line for the wanted station with the numpy.where
                    try:                        
                        groups[k][d.element_name]=({'hdrlen': 1}, np.full( 1 , station_configuration_retrieved[d.element_name].values[0] ) )
                        #print('working sc')
                    except:
                        #print('not working sc')
                        pass


                elif k in ('source_configuration'): # storing the source configuration info, e.g. original file name, 
                    if d.element_name=='source_file':
                        #  groups[k][d.element_name] = ( {'hdrlen':fbds.variables['date@hdr'].shape[0] } ,  np.full( fbds.variables['date@hdr'].shape[0] , source_file  ) ) 
                        groups[k][d.element_name]=({'hdrlen': 1 },   np.full( 1 , source_file) )
                    else:
                        try:   
                            groups[k][d.element_name]=({{'hdrlen': 1 },   np.full (1, cdm[k][d.element_name].values[0] ) # element_name is the netcdf variable name, which is the column name of the cdm table k 
                        except KeyError:
                            pass
                       
   
                else : # this is the case where the cdm tables DO exist
                    try:   
                        groups[k][d.element_name]=({k+'_len':len(cdm[k])}, cdm[k][d.element_name].values) # element_name is the netcdf variable name, which is the column name of the cdm table k 
                    except KeyError:
                        pass
                    
                """ Tryin to add attributes, e.g. description and external tables """ 
                if k not in ('observations_table'):               
                    try:
                        groups[k][d.element_name].attrs['external_table']=d.external_table # defining variable attributes that point to other tables (3rd and 4th columns)
                        groups[k][d.element_name].attrs['description']=d.description  # it faisl when trying with the observations_table 
                    except:
                        pass
                
                try:
                    
                    if type(groups[k][d.element_name].values[0])==str:
                        s=groups[k][d.element_name].values.shape
                        groupencodings[k][d.element_name]={'dtype':numpy.dtype('S80'),'compression': 'gzip','chunksizes':(min(100000, s[0] ) , 80 ) }
                    else:
                        groupencodings[k][d.element_name]={'compression': 'gzip'}
                    
                    if k in ('observations_table'):
                       write_dict_h5(fno, groups[k], k, groupencodings[k], var_selection=[],mode='a', attrs= obs_tab_attrs  )
                       print ('writing observations_table')
                except:
                    #print('bad:',k,d.element_name)
                    pass
        
        
        for k in groups.keys():            
            ##this appends group by group to the netcdf file
            if k not in ('observations_table') :           
                print('ERA5 doing k: ' , k )
                groups[k].to_netcdf(fno,format='netCDF4',engine='h5netcdf',encoding=groupencodings[k],group=k,mode='a') #
        ##print('sizes: in: {:6.2f} out: {:6.2f}'.format(os.path.getsize(fn+'.gz')/1024/1024, os.path.getsize(fno)/1024/1024))
        del fbds
    
    """ Storing the group encodings in a numpy dictionary to be reused by the merging script """
    np.save('groups_encodings',  groupencodings)
    #print(fno,time.time()-t)

    return 0


            
            
            

def load_cdm_tables():
    """ Load the cdm tables into Panda DataFrames, reading the tables from the cdm GitHub page FF To do 
    
          Return:
                      dictionary with a Panda DataFrame for each table """
    
    
    
    """ # Uncomment to get the list of all the .csv files present at the url specified
    url = 'https://github.com/glamod/common_data_model/tree/master/table_definitions'
    cdmtabledeflist = csvListFromUrls(url)
    """
    tpath = os.getcwd() + '/../data'
    cdmpath='https://raw.githubusercontent.com/glamod/common_data_model/master/tables/' # cdm tables            
    
    """ Selecting the list of table definitions. Some of the entires do not have the corresponding implemented tables """
    cdmtabledeflist=['id_scheme', 'crs', 'station_type', 'observed_variable', 'station_configuration', 'station_configuration_codes', 'observations_table', 'header_table', 'source_configuration', 'units' , 'z_coordinate_type']  
    cdm_tabdef = dict()
    for key in cdmtabledeflist:
        url='table_definitions'.join(cdmpath.split('tables'))+key+'.csv' # https://github.com/glamod/common_data_model/tree/master/table_definitions/ + ..._.dat 
        f=urllib.request.urlopen(url)
        col_names=pd.read_csv(f,delimiter='\t',quoting=3,nrows=0,comment='#')
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
        cdm_tab[key]=pd.read_csv(f,delimiter='\t',quoting=3,dtype=tdict,na_filter=False)


    """ Adding the  tables that currently only have the definitions but not the implementation in the CDM, OR    need extensions """  
    cdm_tabdef['header_table']          = pd.read_csv(tpath+'/table_definitions/header_table.csv',delimiter='\t',quoting=3,comment='#')
    cdm_tabdef['observations_table'] = pd.read_csv(tpath+'/table_definitions/observations_table.csv',delimiter='\t',quoting=3,comment='#')

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
    dic_obs = {}
    for index, row in cdm_tabdef['observations_table'].iterrows():
        dic_obs[row['element_name'] ] = []
        dic_obs[row['element_name'] ] = ['description' ,  row.description ]     
    dic_obs['date_time'] = ['units',  'seconds since 1900-01-01 00:00:00' ]
    

    return cdm_tabdef  , cdm_tab, tdict , dic_obs 



def csvListFromUrls(url=''):
    """ Return a list of csv files, as fond in the url on the cdm GitHub """   
    urlpath = urlopen(url)
    string = urlpath.read().decode('utf-8')
    split = string.split(' ')
    csv_files_list = [m.replace('"','') for m in [n.split('title="')[1] for n in split if '.csv' in n and "title" in n] ] 
    return csv_files_list



def filelist_cleaner(lista, dataset=''):
    """ Removes unwanted files that might be present in the database directories """
    if dataset == 'ncar':
        cleaned = [ l for l in lista if '.nc' not in l ]
    if dataset == 'bufr':
        cleaned = [ l for l in lista if '.bfr' in l ]
    if 'era5' in dataset:
        cleaned = [ l for l in lista if '.nc' not in l and '.conv.' in l ]
    else:
        cleaned = lista
    
    return cleaned
   



def clean_station_configuration(cdm_tab ):
    """ Replace wrong characters from the station cofniguration tables """
    subs={'o':[240,242,243,244,245,246,248],'O':[210,211,212,213,214,216],
          'a':[224,225,226,227,228,229,230],'A':[192,193,194,195,196,197,198],
          'u':[249,250,251,252,253],'U':[217,218,219,220],
          'i':[236,237,238,239],'I':[204,205,206,207,304],
          'S':[350],'n':[241],'c':[231],'C':[199],'e':[232,233,234,235],'E':[200,201,202,203]}
    for k in cdm_tab['station_configuration'].columns:
        #print(k,type(cdm['station_configuration'][k][0]))
        if type(cdm_tab['station_configuration'][k][0]) is str:
            try:               
                cdm_tab['station_configuration'][k].values[:]=cdm_tab['station_configuration'][k].values[:].astype('S')
            except:
                for l in range(cdm_tab['station_configuration'][k].values.shape[0]):
                    try:                       
                        cdm_tab['station_configuration'][k].values[l]=numpy.string_(cdm_tab['station_configuration'][k].values[l])
                    except:
                        for m in range(len(cdm_tab['station_configuration'][k].values[l])):
                            mychar=cdm_tab['station_configuration'][k].values[l][m]
                            if ord(mychar)>128:
                                for n,v in subs.items():
                                        if ord(mychar) in v:
                                            cdm_tab['station_configuration'][k].values[l]=n.join(cdm_tab['station_configuration'][k].values[l].split(mychar))
                        
                        cdm_tab['station_configuration'][k].values[l]=numpy.string_(cdm_tab['station_configuration'][k].values[l])
                cdm_tab['station_configuration'][k]=numpy.string_(cdm_tab['station_configuration'][k])

    print('Cleaned station_configuration')


    

db   = { 'era5_1'       : { 'dbpath' : '/raid60/scratch/leo/scratch/era5/odbs/1'               } ,
                  'era5_3188' : { 'dbpath' : '/raid60/scratch/leo/scratch/era5/odbs/3188'     } ,
                  'era5_1759' : { 'dbpath' : '/raid60/scratch/leo/scratch/era5/odbs/1759'     } ,
                  'era5_1761' : { 'dbpath' : '/raid60/scratch/leo/scratch/era5/odbs/1761'     } ,
                  'ncar'           : { 'dbpath' : '/raid60/scratch/federico/databases/UADB'        } ,
                  'ncar'           : { 'dbpath' : '/raid60/scratch/federico/databases/UADB'        } ,
                  'igra2'          : { 'dbpath' : '/raid60/scratch/federico/databases/IGRAv2'      } ,
                  'bufr'            : { 'dbpath' : '/raid60/scratch/leo/scratch/era5/odbs/ai_bfr'    }   }
        

    

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Make CDM compliant netCDFs")
    parser.add_argument('--dataset' , '-d', 
                    help="Select the dataset to convert. Available options: all, era5_1, era5_1759, era5_1761, bufr, igra2, ncar, test. If not selected or equal to 'test', the script will run the example files in the /examples directory."  ,
                    default = 'test',
                    type = str)
    
    parser.add_argument('--output' , '-o', 
                    help="Select the output directory. If not selected, converted files will be stored in the 'converted_files' directory."  ,
                    default = 'converted_files',
                    type = str)    
    
    parser.add_argument('--files' , '-f',
                    help = "File to be processed."  ,
                        default = '',
                        type = str)    
         
    args = parser.parse_args()
    dataset = args.dataset 
    out_dir = args.output
    Files = args.files

    if dataset not in ['era5_1', 'era5_2', 'era5_3188', 'era5_1759', 'era5_1761', 'bufr', 'igra2', 'ncar']:
        raise ValueError(" The selected dataset is not valid. Please choose from ['era5_1', 'era5_1759', 'era5_1761', 'era5_3188', 'bufr', 'igra2', 'ncar' ]  ")    
                        
    """ Loading the CDM tables into pandas dataframes """
    cdm_tabdef  , cdm_tab , tdict, obs_tab_attrs = load_cdm_tables()
     
    """ Paths to the output directory """    
    if not os.path.isdir(out_dir):
            os.system('mkdir ' + out_dir )       
            
    output_dir = out_dir + '/' + dataset           
    if not os.path.isdir(output_dir):
            os.system('mkdir ' + output_dir )
    
    stat_conf_path = '../data/station_configurations/'     
    stat_conf_file = stat_conf_path +   '/station_configuration_' + dataset + '.dat'    
    cdm_tab['station_configuration']=pd.read_csv(stat_conf_file,  delimiter='\t', quoting=3, dtype=tdict, na_filter=False, comment='#')
    clean_station_configuration(cdm_tab)              
            
            
    """ Leo run         
    Files=glob.glob(Files)
    print( blue + '*** Processing the database ' + dataset + ' ***  \n \n *** file: ' + Files[0] + '\n'  + cend)
    stat_conf_path = '../data/station_configurations/'     
    stat_conf_file = stat_conf_path +   '/station_configuration_' + dataset + '.dat'
    # adding the station configuration to the cdm tables      
    cdm_tab['station_configuration']=pd.read_csv(stat_conf_file,  delimiter='\t', quoting=3, dtype=tdict, na_filter=False, comment='#')
    clean_station_configuration(cdm_tab)             
    p=Pool(12)
    if 'era5' in dataset and 'bufr' not in dataset:   
        func=partial(odb_to_cdm,cdm_tab, cdm_tabdef, output_dir, dataset)
        out=list(map(func,Files))
    else:
        func=partial(df_to_cdm,cdm_tab, cdm_tabdef, output_dir, dataset)
        out=list(map(func, Files))
    print('*** CONVERTED: ' , Files[-1] )      
    print(' ***** Convertion of  ' , Files,  '  completed ! ***** ')
    """
    
    
    """ Federico run """
    Files = Files.split(',')
    for File in Files:
             
        if not os.path.isdir(out_dir):
            os.system('mkdir ' + out_dir ) 
            
        output_dir = out_dir + '/' + dataset      
        if not os.path.isdir(output_dir):
            os.system('mkdir ' + output_dir )
                    
        print( blue + '*** Processing the database ' + dataset + ' ***  \n \n *** file: ' + File + '\n'  + cend)
                                    
        if 'era5' in dataset and 'bufr' not in dataset:   
            odb_to_cdm( cdm_tab, cdm_tabdef, output_dir, dataset,  obs_tab_attrs, File)
        else:
            df_to_cdm( cdm_tab, cdm_tabdef, output_dir, dataset, File)

          
        print(' ***** Convertion of  ' , File ,  '  completed ! ***** ')   
        
    

      







""" Examples for running 
-f /raid60/scratch/leo/scratch/era5/odbs/1/era5.conv._82930 -d era5_1 -o OUTPUT
-f /raid60/scratch/leo/scratch/era5/odbs/1759/era5.1759.conv.1:82930.gz -d era5_1759 -o OUTPUT
-f /raid60/scratch/federico/databases/IGRAv2/BRM00082930-data.txt -d igra2 -o OUTPUT
-f /raid60/scratch/federico/databases/UADB//uadb_windc_82930.txt -d ncar -o OUTPUT
-f /raid60/scratch/leo/scratch/era5/odbs/ai_bfr/era5.82930.bfr -d bufr -o OUTPUT
-f /raid8/srvx1/federico/GitHub/CEUAS_master_FEB202/CEUAS/CEUAS/public/harvest/data/example_stations/era5_3188/era5.3188.conv.C:4629  -d era5_3188 -o OUTPUT

-f /raid60/scratch/leo/scratch/era5/odbs/3188/era5.3188.conv.C:8022  -d era5_3188 -o OUTPUT

-f /raid60/scratch/federico/databases/UADB//uadb_windc_27962.txt  -d  ncar  -o OUTPUT  (1553915536 vs 29374566 )

"""