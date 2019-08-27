#!/usr/bin/env python
import sys, os.path, glob 
import subprocess
import urllib.request
import xarray as xr
import numpy
#import h5pickle as h5py
import h5py
from datetime import date, datetime
import time
from multiprocessing import Pool
from netCDF4 import Dataset
import gzip
import pandas as pd    
from functools import partial
#from rasotools.utils import *
from numba import *
import matplotlib.pylab as plt
import cartopy.crs as ccrs
import argparse
import copy
from io import StringIO
import h5netcdf
import numpy as np
from eccodes import *

debug = False


""" Some colors for pretty printout """ 
red    = '\033[91m' 
cend   = '\033[0m'
blue   = '\033[34m'
green  = '\033[92m'
yellow = '\033[33m'

""" Possible variable types as listed int he CMD tables S80 is a fixed 80 char. string. Compression does not worw well with normal strings. """
okinds={'varchar (pk)':numpy.dtype('|S80'),'varchar':numpy.dtype('|S80'),'numeric':numpy.float32,'int':numpy.int32,
       'timestamp with timezone':numpy.datetime64,
       'int[]*':list,'int[]':list,'varchar[]*':list,'varchar[]':list}

""" Variable types to be used in the compressed netCDF files """
kinds={'varchar (pk)':str,'varchar':str,'numeric':numpy.float32,'int':numpy.int32,
       'timestamp with timezone':numpy.datetime64,
       'int[]*':list,'int[]':list,'varchar[]*':list,'varchar[]':list}


def make_datetime(dvar,tvar):
    """ Converts into date-time standard format """
    dvari=dvar.values.astype(numpy.int)
    tvari=tvar.values.astype(numpy.int)
    df=pd.DataFrame({'year':dvar//10000,'month':(dvar%10000)//100,'day':dvar%100,
                        'hour':tvar//10000,'minute':(tvar%10000)//100,'second':tvar%100})
    dt=pd.to_datetime(df).values
    return dt


""" Translates some of the odb variables name  into cdm var """
cdmfb={'observation_value':'obsvalue@body',
       'observed_variable':'varno@body',
       'z_coordinate_type':'vertco_type@body',
       'z_coordinate':'vertco_reference_1@body',
       'date_time':[make_datetime,'date@hdr','time@hdr'],
       'longitude':'lon@hdr',
       'latitude':'lat@hdr'}

def check_read_file(file='', read= False):
    """ Simple utility to check if file exists and uncompress it, then optionally read the lines (in case of text files e.g. igra2 and UADB)
        and store them as entry of a list. Return the list.
        Used to prepare the igra2 and UADB files. 
        Adapted from https://github.com/MBlaschek/igra/blob/master/igra/read.py 
        
        Args:
             file (str): path to the file
            read (bool): set to True to return the list of lines in case of txt files 
        Returns:
             lines (list) read from the file """

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


""" Dictionary mapping generic names of the variables to the numbering scheme in the CDM """
cdmvar_dic = {'temperature'         : 85, 
                         'wind_direction'     : 106,  
                         'wind_speed'          : 107, 
                         'dew_point'            : 36, 
                         'relative_humidity' : 38 ,
                         'pressure'               : 117 }



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
    cnt = 0
    
    
    """ Dicts mapping the names of the variables in the bufr files (keys) to the cdm names (values) """
    cdm_header_dic =  {'stationNumber'  : 'station_name' ,
                       'blockNumber'    : ''             ,
                       'heightOfStation': ''             ,
                       'latitude'       : 'latitude'     ,
                       'longitute'      : 'longitude'   ,
                       'typicalDate': '', 'typicalTime':'' 
                       }

    #cdmvar_dic = { 'airTemperature' : {'units':'C_to_tenths'       , 'cdm_name': 'temperature'         },
    #              'windDirection'  : {'units':'ms_to_tenths'      , 'cdm_name': 'wind_speed'          },
    #                'windSpeed'      : {'units':'degree'            , 'cdm_name': 'wind_direction'      },
    #                'dewpointTemperature' : {'units': ' '           , 'cdm_name': 'dew_point'           },
    #                'pressure'       : {'units': ' '                , 'cdm_name': 'pressure'            } }
    
    
    
    bufr_values = []
    
    """ Name of the columns as they will appear in the pandas dataframe (not necessarily CDM compliant) """
    column_names = ['report_timestamp' , 'iday',  'station_id', 'latitude', 'longitude', 'pressure', 'varno@body', 'obsvalue@body']
    
    while 1:
        lista = [] # temporary list
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
        
        pressure            = codes_get_array(bufr, "pressure") 
        temperature      = codes_get_array(bufr, "airTemperature")           
        wind_direction    = codes_get_array(bufr, "windDirection")
        wind_speed        = codes_get_array(bufr, "windSpeed")
        
        lat, lon, alt, blockNumber, stationNumber, statid = 0,0,0,0,0,0
        
        if cnt ==0:

            lat                     = codes_get(bufr, "latitude")
            lon                    = codes_get(bufr, "longitude")
            alt                     = float(codes_get(bufr, "heightOfStation"))
            blockNumber    = codes_get(bufr, "blockNumber")
            stationNumber = codes_get(bufr, "stationNumber")
            statid                = str(blockNumber*1000+stationNumber)
        
        codes_release(bufr)
   
        miss_value = -1.e100     
        
        for i in range(len(temperature)): 
            airT         = temperature[i]
            winds      = wind_speed[i]
            windd      = wind_direction[i]
            press       = pressure[i]
            
            if airT == miss_value :    # replacing none values with numpy nans
                airT = np.nan 
            if winds == miss_value:
                winds = np.nan
            if windd == 2147483647:
                windd = np.nan 
                
            # column_names = ['report_timestamp' , 'iday',  'station_id', 'latitude', 'longitude', 'varno@body', 'obsvalue@body']    
            for value,var in zip([airT, winds, windd],  ['temperature', 'wind_speed', 'wind_direction'] ):
                   bufr_values.append( (idate, iday, statid, lat, lon, press, cdmvar_dic[var] , value) )
        
            cnt += 1
               
    df = pd.DataFrame(data= bufr_values, columns= column_names)
    
    return df.to_xarray()
    
    



def uadb_ascii_to_dataframe(file=''):
    """ Read an uadb stationfile in ASCII format and convert to a Pandas DataFrame.                                                                                                                                                                                          
        Adapted from https://github.com/MBlaschek/CEUAS/tree/master/CEUAS/data/igra/read.py                                                                                                                                                                                         
        Variables used inside the DataFrame are already CDM compliant   
        
        Args:
             file (str): path to the uadb station file

        Returns:
             Pandas DataFrame with cdm compliant column names
    """     
    if debug:
         print("Running uadb_ascii_to_dataframe for: ", file)    
         
    data = check_read_file(file=file, read=True)

    raw = []
    headers = []
    dates = []
    nmiss = 0
    iprev = 0
    search_h = False
    i = 0
     
    read_data = []

    idate, usi, lat, lon, lat, stype, press, gph, temp, rh, wdir, wspd = 0,0,0,0,0,0,0,0,0,0,0,0

    for i, line in enumerate(data):
        if line[0] == 'H':
            try:
                # Header
                usi      = int(line[2:14])  # unique station identifier
                ident    = line[15:21]  # WMO
                idflag   = int(line[22:24])  # id flag
                d_src    = int(line[25:28])  # source dataset
                version  = float(line[29:34])  # version
                dateflag = int(line[35:37])  # date flag
                year     = line[38:42]  # year
                month    = "%02d" % int(line[43:45])
                day      = "%02d"  % int(line[46:48])
                hour     = line[49:53]
                locflag  = int(line[54:56])  # Location Flag
                lat      = float(line[57:67])
                lon      = float(line[68:78])
                ele      = float(line[79:85])
                stype    = int(line[86:88])
                numlev   = int(line[89:93])
                pvers    = line[94:102]

                if '99' in hour:
                    hour = hour.replace('99', '00')
        
                if '99' in day:
                    search_h = True
                    continue

                hour = "%02d" % (int(hour) // 100)
                minutes = int(hour) % 100
                if minutes > 60 or minutes < 0:
                    minutes = 0
                minutes = "%02d" % minutes
                idate = datetime.strptime(year + month + day + hour + minutes, '%Y%m%d%H%M')
                iday = int(year + month + day)
                #headers.append((idate, usi, numlev, lat, lon, ele, stype))
                pday = int(day)
                search_h = False

            except Exception as e:
                print("Error: ", i, line, repr(e), "Skipping Block:")
                #if kwargs.get('debug', False):
                #    raise e

                search_h = True
                iprev = i

        elif search_h:
            nmiss += 1
            continue  # Skipping block

        else:
            # Data
            ltyp   = int(line[0:4])
            press  = float(line[5:13])   # hPa
            gph    = float(line[14:22])
            temp   = float(line[23:29])  # degree
            rh     = float(line[30:36])  # %
            wdir   = float(line[37:43])
            wspd   = float(line[44:50])  # m/s

        for value,var in zip([temp, wspd, wdir, rh],  ['temperature', 'wind_speed', 'wind_direction', 'relative_humidity'] ):
               read_data.append( (idate, iday, ident, lat, lon, press, cdmvar_dic[var] , value) )
                            
    column_names = ['report_timestamp' , 'iday',  'station_id', 'latitude', 'longitude', 'pressure', 'varno@body', 'obsvalue@body']

    df = pd.DataFrame(data= read_data, columns= column_names)
    #pdf = pdf.replace([-999.9, -9999, -999, -999.0, -99999.0, -99999.9], np.nan)
    
    df['pressure'] *= 100.  # need Pa
    df.sort_values(by = ['iday' , 'report_timestamp'] )    
    
    return df.to_xarray()


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
                 
    read_data = [] #  Lists containing the raw data from the ascii file, and the observation dates

    """ Data to be extracted and stored from the igra2 station files 
        Some info is contained in the header of each ascent, some in the following data """
    #columns=['ident','year','month','day','hour','reltime','p_src','np_src','lat','lon',
     #        'lvltyp1', 'lvltyp2', 'etime', 'press', 'pflag', 'gph', 'zflag', 'temp', 'tflag', 'rh' ,'dpdep', 'wdir','wspd']
 
    column_names = ['report_timestamp' , 'iday',  'station_id', 'latitude', 'longitude', 'pressure', 'varno@body', 'obsvalue@body']
    
 
    """ Initialize the variables that can be read from the igra2 files """
    ident,year,month,day,hour,reltime,p_src,np_src,lat, lon = 0,0,0,0,0,0,0,0,0,0
    lvltyp1,lvltyp2,etime,press,pflag,gph,zflag,temp,tflag,rh,dpdep,wdir,wspd = 0,0,0,0,0,0,0,0,0,0,0,0,0 # initialize to zeros
    idate = 0
    
    for i, line in enumerate(data):
        if line[0] == '#':
            # Info from the Header line of each ascent                                                                                                                                                                                                                   
            ident     = line[1:12]               # station identifier
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

            if int(hour) == 99:
                time = reltime + '00'
            else:
                time = hour + '0000'
       
            if '99' in time:
                time = time.replace('99', '00')

            idate = datetime.strptime(year + month + day + time, '%Y%m%d%H%M%S') # constructed according to CDM
            iday =  int(year + month + day)
        else:
           # Data of each ascent
            lvltyp1 = int(line[0])            # 1-  1   integer major level type indicator
            lvltyp2 = int(line[1])            # 2-  2   integer minor level type indicator
            etime   = int(line[3:8])          # 4-  8   integer elapsed time since launch
            press   = int(line[9:15])         # 10- 15  integer reported pressure
            pflag   = line[15]                # 16- 16  character pressure processing flag
            gph     = int(line[16:21])        # 17- 21  integer geopotential height
            zflag   = line[21]                # 22- 22  character gph processing flag
            temp    = int(line[22:27]) / 10.  # 23- 27  nteger temperature
            tflag   = line[27]                # 28- 28  character temperature processing flag
            rh      = int(line[28:33]) / 10.  # 30- 34  integer relative humidity
            dpdp    = int(line[34:39]) / 10.  # 36- 40  integer dew point depression (degrees to tenth e.g. 11=1.1 C)
            wdir    = int(line[40:45])        # 41- 45  integer wind direction (degrees from north, 90 = east)
            wspd    = int(line[46:51]) / 10.  # 47- 51  integer wind speed (meters per second to tenths, e.g. 11 = 1.1 m/s 

        
            for value,var in zip([temp, wspd, wdir, rh, dpdp],  ['temperature', 'wind_speed', 'wind_direction', 'relative_humidity', 'dew_point'] ):
                   read_data.append( (idate, iday, ident, lat, lon, press, cdmvar_dic[var] , value) )


    df = pd.DataFrame(data= read_data, columns= column_names)

    return df.to_xarray()

 


    
def read_all_odbsql_stn_withfeedback(odbfile):
    
    if debug: 
        print("Running read_all_odbsql_stn_withfeedback for: ", odbfile)
        
    alldata=''

    alldict=xr.Dataset()
    t=time.time()
    sonde_type=True
    obstype=True
    if os.path.getsize(odbfile)>0:
        """ Read first the odb header to extract the column names and type """
        try:
            rdata=subprocess.check_output(["odb","header",odbfile])
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
                        elif kinds[-1] in ('INTEGER','BITFIELD'):
                            if columns[-1]=='date@hdr':
                                tdict[columns[-1]]=numpy.int32
                            else: 
                                tdict[columns[-1]]=numpy.float32
                        else:
                            tdict[columns[-1]]=numpy.dtype('S8') # dict containing column name and type
                                     
                            
                except IndexError:
                    pass
        except:
            print('could not read odbfile '+odbfile)
            return alldict
        try:
            rdata=subprocess.check_output(["odb","sql","-q","select *","-i",odbfile,'--no_alignment']) # after reading the header it does the query
            # returns a byte string
            #print('after odb:',time.time()-t)
            rdata=''.join(rdata.decode('latin-1').split("'")) # decoding the string into a unicode
            f=StringIO(rdata) # access the string like a file, return a file pointer to read the string with pandas
            # nb  if you have null values, reading of integer fails and are read as floats
            # to improve, you can convert the columns with nans in the alldicts (pandas data frame) into int(np.nans)
            # date values are large so float 32 precision is not sufficient  FF TODO
            alldict=pd.read_csv(f,delimiter='\t',quoting=3,comment='#',dtype=tdict)
            del f,rdata

 
            """ alternative method to read the odb           
            if False:
                
                rdata='nan'.join(rdata.decode('latin-1').split('NULL'))
                rdata=''.join(rdata.split("'"))
                rdata='\t'.join(rdata.split('\n')[1:-1])
                rdata=tuple(rdata.split('\t'))
                
                print('after odb:',time.time()-t)
                #xdata=numpy.fromstring(rdata,sep='\t')
                cl=len(columns)
                rl=len(rdata)//cl
                #for k in range(cl):
                    #if kinds[k]=='REAL':
                        #alldict[columns[k]]=({'obslen':rl},numpy.empty(rl,dtype=numpy.float32))
                    #elif kinds[k] in ('INTEGER','BITFIELD'):
                        #alldict[columns[k]]=({'obslen':rl},numpy.empty(rl,dtype=numpy.int32))
                    #else:
                        #alldict[columns[k]]=({'obslen':rl},numpy.empty(rl,dtype='|S8'))
                #print(odbfile,time.time()-t)
                for k in range(cl):
                    if kinds[k]=='REAL':
                        alldict[columns[k]]=({'obslen':rl},numpy.float32(rdata[k::cl]))
                    elif kinds[k] in ('INTEGER','BITFIELD'):
                        rds=rdata[k::cl]
                        if 'nan' in rds: 
                            alldict[columns[k]]=({'obslen':rl},numpy.asarray(rds,dtype=float).astype(numpy.int32))
                        else:
                            alldict[columns[k]]=({'obslen':rl},numpy.asarray(rds,dtype=numpy.int32))
                    else:
                        alldict[columns[k]]=({'obslen':rl},numpy.asarray(rdata[k::cl],dtype='|S8'))
             """
       
            #input('FF check alldict')
            # alldict.head() # pandas method to print the columns name . alldict is now a panda dataframe with e.g. 22000*63 columns where the 63 columns are read from the odb file 
            
        except subprocess.CalledProcessError as e:
            print('odb failed!:'+' '+odbfile)
            return alldict

    #print(odbfile,time.time()-t)
    #idy=numpy.lexsort((alldict['varno@body'],
                       #-alldict['vertco_reference_1@body'],
                       #alldict['time@hdr'],
                       #alldict['date@hdr']))
    #for k in alldict.columns:
        #alldict[k]=alldict[k][idy]

    #print(odbfile,time.time()-t)

    """ may not be necessary to convert into x_array sicne you can write a pandas df into an HDF file """

    #print('FF alldict_to_xarray is: ', alldict.to_xarray )
   # input('verifica alldicts.to_xarray')

    return alldict.to_xarray()


            
def fromfb(fbv, cdmfb):
    """ input: 
               fbv    : xarray converted from the original input file 
               cdmfb  :  
    """
    x=0
    # checks if the type of the variable is a list, so that it uses the function to extract the date time 
    if type(cdmfb) is list:
        x=cdmfb[0](fbv[cdmfb[1]], fbv[cdmfb[2]])
    else:
        if cdmfb=='varno@body':
            tr=numpy.zeros(113,dtype=int) # fixed length of 113 since the highest var number is 112 
            """ translate odb variables (left) number to CDM numbering convention (right) """
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
            tr[39]= 85 # 2m T
            tr[40]= 36 # 2m Td
            tr[41]= 104 #10m U
            tr[42]= 105  #10m V
            tr[58]=38 # 2m rel hum
            
            x=tr[fbv[cdmfb].values.astype(int)] # reads the varno from the odb feedback and writes it into the variable id of the cdm 
        else:    
            x=fbv[cdmfb].values
        
    return x

def ttrans(cdmtype, kinds=kinds):
    """ convert the cdm types to numpy types """    
    nptype=numpy.float32
    try:
        nptype=kinds[cdmtype.strip()]
    except:
        print(cdmtype,'not found, using numpy.float32')   
    return nptype

#@njit
def find_dateindex(y):
    """ creates the indices list from the dates, for quick access 
        nb the benchmark script will not work with these files since the definition of the array size is swapped i.e. (x.shape[0], 3)"""        

    
    if debug: 
        print("Running find_dateindex for: ", y)
        
    x=numpy.unique(y)
    z=numpy.zeros((3, x.shape[0]), dtype=numpy.int32 )
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

def find_date_indices(iday):
    """ Extracts the list of observation dates, and store the indices of the first and last observations 
          Args: 
                  list of observation dates
          Return: numpy array (3,:) where the first list is the list of unique obs dates, the second are the indices of first obs  dates, the second the list of last obs dates. 
          """
    if debug: 
        print("Running find_dateindex")
        
    days, counts = numpy.unique(iday, return_counts = True)
    z=numpy.zeros((3, days.shape[0]), dtype=numpy.int32 )    

    first, last = [], []    
    
    for day,count in zip(days, counts):
          f = int( min(np.where(iday == day )[0] ) )
          
          first.append(f)
          last.append(f + count - 1) 

    z[0,:] = days
    z[1,:] = np.array(first)
    z[2,:] = np.array(last)
    
    return z 
          
def readvariables_encodings(fbds):
    """ Extracts the list of variables from the xarray, read from different sources, and assign a variable type for the encodings 
        
        Args:
               pandas dataframe 
        Return:
               dictionary where each key is a cdm variable, and the values are the variable type and compression methods to be used in the hdf files """
    
    fbencodings={}
    for d in fbds._variables.keys():
        if fbds.variables[d].dtype==numpy.dtype('float64'):
            if d!='date@hdr' or d!='report_timestamp':             
                fbencodings[d]={'dtype':numpy.dtype('float32'),'compression': 'gzip'} # !! probably dtype not necessary, but compression must be there
            else:
                print('var is timestamp!!', d)
                fbencodings[d]={'dtype':numpy.dtype('int32'),'compression': 'gzip'}               
        else:
                fbencodings[d]={'compression': 'gzip'}
                
    return  fbencodings          
                
def initialize_convertion(fn, output_dir):
     """ Simple initializer for writing the output netCDF file """
     
     fnl=fn.split('/')
     fnl[-1]='ch'+fnl[-1]
     fno=output_dir + '/' + fnl[-1] + '.nc' # creating an output file name e.g. chera5.conv._10393.nc  , try 01009 faster
     return fno 
    
#def get_encodings(cdm, cdmd,fdbs):
#    """ Extract the correct encodings for the variables in the tables """
    
    

def df_to_cdm(cdm, cdmd, out_dir, fn):
    """ Convert the  pandas dataframe from the file into cdm compliant netCDF files. Use with bufr, igra2 and ncar databases.
        According to each file type, it will use the appropriate reading function to extract a Pandas DataFrame
        input:
              fn       :: odb file name (e.g. era5.conv._10393)
              cdm   :: cdm tables (read with pandas)
              cdmd :: cdm tables definitions ("")  """    
    
    if debug:
        print("Running df_to_cdm for: ", fn)
        
    station_id_fails = open('station_id_fail.log' , 'a') 
    
    t=time.time()
    fno = initialize_convertion(fn, out_dir ) 
    log_file = open(out_dir +  '/' + 'logger.log', 'a') 
    if not False:
        
            # era5 analysis feedback is read from compressed netcdf files era5.conv._?????.nc.gz in $RSCRATCH/era5/odbs/1
            """ Reading the odb and convert to xarray """  
            if  '.bfr' in fn:
                fbds= bufr_to_dataframe(fn) # fdbs: the xarray converted from the pandas dataframe 
            elif  'uadb' in fn:
                fbds = uadb_ascii_to_dataframe(fn)
            elif  'data' in fn:
                #print(" analysing an igra dataset")
                fbds = igra2_ascii_to_dataframe(fn)
            else:
                print('Unidentified file is: ', fn)
                raise ValueError('Cannot identify the type of file to be analized!!! ')
                            
            """ Extract the unique indices of each date observation. Create an xarray to be saved as a separate variable in the netCDF """                      
            y=fbds['iday'].values # variables iday contains only dates and no time information
            #zz=find_dateindex(y)
            z = find_date_indices(y)
            di=xr.Dataset() 
            di['dateindex']= ( {'days' : z.shape[1] , 'drange' : z.shape[0]},  z)   # di xarray will be saved to netCDF         
                                 
            """ Extracting the variables type encodings """
            fbencodings = readvariables_encodings(fbds)    
    
            """ Here we extract the name of the variables from the xarray, and we check the type. 
                  According to it, we set the type of the variable in the hdf compressed variable and store the encoding in the fbencodings dictionary
                  Looks like fbencodings = { 'obsvalue@body': {'compression': 'gzip'}, 'sonde_type@conv': {'compression': 'gzip'}, ... } 
 
                  Empty dictionaries for hdf groups (one for each table) and their  hdf encoding description.
                  Each entry in the group dict is an xarray mapping to a single table_definition """
            
            groups={} 
            groupencodings={}
        
            """" Loop over every table_definition in the cdmd dictionary of pandas DF 
                   Two possibilities: 1. the corresponding  table already exists in the cdm (case in the final else)
                                              2. the corr table needs to be created from the local data sources (e.g. the feedback or IGRA or smt else). 
                                                  e.g. the observation_tables, the header_tables and the station_configuration.
                                                  These tables are contained in the CEUAS GitHub but not in the cdm GitHub """
            
            for k in cdmd.keys(): # loop over all the table definitions e.g. ['id_scheme', 'crs', 'station_type', 'observed_variable', 'station_configuration', 'station_configuration_codes', 'observations_table', 'header_table']
                groups[k]=xr.Dataset() # create an xarray for each table
                groupencodings[k]={} # create a dict of group econding for each table
                
                """ Loop over each row in the definition table.
                     In the CDM, the def. tables have the 4 columns ['element_name', 'kind', 'external_table', 'description' ] 
                     d is each element in the DF, so e.g. d.external_table will return the correpsonding value """    
                
                for i in range(len(cdmd[k])): 
                    d=cdmd[k].iloc[i] # extract the element at the i-th row, e.g.   element_name           z_coordinate_type, kind                                 int, external_table    z_coordinate_type:type, description         Type of z coordinate
                                                                             
                    if k in ('observations_table'):
                        try:                                                     
                            groups[k][d.element_name]=( {'hdrlen':fbds.variables['report_timestamp'].shape[0] }, fromfb(fbds._variables, cdmfb[d.element_name] ) )
                            
                        except KeyError:
                            x=numpy.zeros(fbds.variables['report_timestamp'].values.shape[0], dtype= numpy.dtype(ttrans(d.kind, kinds=okinds)))
                            x.fill(numpy.nan)
                            groups[k][d.element_name]= ({'hdrlen':fbds.variables['report_timestamp'].shape[0]}, x)
                            
                    elif k in ('header_table'):
                        # if the element_name is found in the cdmfb dict, then it copies the data from the odb into the header_table
                        try:
                            groups[k][d.element_name]= ({'hdrlen':fbds.variables['report_timestamp'].shape[0] }, fromfb(fbds._variables, cdmfb[d.element_name] ) )
                        except KeyError:
                            # if not found, it fills the columns with nans of the specified kind. Same for the observation_tables 
                            x= numpy.zeros(fbds.variables['report_timestamp'].values.shape[0],dtype=numpy.dtype(ttrans(d.kind,kinds=okinds) ) )
                            x.fill(numpy.nan)
                            groups[k][d.element_name]= ({'hdrlen':fbds.variables['report_timestamp'].shape[0]}, x)
                            
                    elif k in ('station_configuration'): # station_configurationt contains info of all the stations, so this extracts only the one line for the wanted station with the numpy.where
                        try:   
                            if 'sci' not in locals(): 
                                #sci=numpy.where(cdm[k]['primary_id']=='0-20000-0-'+fbds['station_id'].values[0].decode('latin1'))[0]
                                station_id = str( fbds['station_id'].values[0].replace(' ','') )
                                sci=numpy.where(cdm[k]['primary_id']=='0-20000-0-'+ station_id) [0]
                                
                            if len(sci)>0:
                                groups[k][d.element_name]=({k+'_len':1},  cdm[k][d.element_name].values[sci] )
                                
                            else:  
                                log_file.write('stationid_' + str(station_id_) + fn + '\n')
                                print( red + ' The station id ', station_id, ' for the file: ', fn, 'could not be found in the station_configuration! Please check! ' + cedn )
                                #input('check station id')
     
                        except KeyError:
                            log_file.write('station_id_keyerror_' + fn + '\n' )
                            pass
                            
                    else : # this is the case where the cdm tables DO exist in th CDM GitHub 
                        try:   
                            groups[k][d.element_name]=({k+'_len':len(cdm[k] ) }, cdm[k][d.element_name].values)  # element_name is the netcdf variable name, which is the column name of the cdm table k 
                        except KeyError:
                            log_file.write('table_error_notexisting_' + fn + '\n')
                            pass
                    try:
                        groups[k][d.element_name].attrs['external_table'] = d.external_table # defining variable attributes that point to other tables (3rd and 4th columns)
                        groups[k][d.element_name].attrs['description']    = d.description
                        #print('good element in cdm table: ' , k, d.element_name ) 
                        groupencodings[k][d.element_name] = {'compression': 'gzip'}
                    except KeyError:
                        log_file.write('k_d.element_name_error_'  + fn + '\n') 
                        #print('bad:', k, d.element_name)
                        pass
    
            #this writes the dateindex to the netcdf file. For faster access it is written into the root group
            """ Wiriting the di (date index) as a separate xarray. """
            di.to_netcdf(fno, format='netCDF4', engine='h5netcdf', mode='w')
   
            """ Writing the content of the original odb file to the netCDF. For each variable, use the proper type encoding."""
            fbds.to_netcdf(fno, format= 'netCDF4', engine= 'h5netcdf', encoding= fbencodings, group= 'era5fb', mode= 'a')
                        
            """ Writing each separate CDM table to the netCDF file """
            for k in groups.keys():
                groups[k].to_netcdf(fno, format='netCDF4', engine='h5netcdf', encoding=groupencodings[k], group=k, mode='a') 
                
            print('sizes: in: {:6.2f} out: {:6.2f}'.format( os.path.getsize( fn) /1024/1024, os.path.getsize( fno )/1024/1024) )
            log_file.close()
            del fbds
        
    #print(fno,time.time()-t)
    
    return 0



def odb_to_cdm(cdm, cdmd, output_dir, fn):
    """ Convert the  file into cdm compliant netCDF files. 
        According to each file type, it will use the appropriate reading function to extract a Pandas DataFrame
        input:
              fn   :: odb file name (e.g. era5.conv._10393)
              cdm  :: cdm tables (read with pandas)
              cdmd :: cdm tables definitions ("")  """

    if debug:
        print("Running odb_to_cdm for: ", fn)
        
    t=time.time()
    fno = initialize_convertion(fn, output_dir) 
    
    if not False:
        
            # era5 analysis feedback is read from compressed netcdf files era5.conv._?????.nc.gz in $RSCRATCH/era5/odbs/1
            """ Reading the odb and convert to xarray """                      
            fbds=read_all_odbsql_stn_withfeedback(fn) # fdbs: the xarray converted from the pandas dataframe 
        
        
            """ Extract the unique indices of each date observation. Create an xarray to be saved as a separate variable in the netCDF """                      
            y=fbds['date@hdr'].values
            z=find_dateindex(y)
            di=xr.Dataset() 
            di['dateindex']= ( {'days' : z.shape[1] , 'drange' : z.shape[0]},  z)   # di xarray will be saved to netCDF         
                        

            
            """ Extracting the variables type encodings """
            fbencodings = readvariables_encodings(fbds)
            
            """ Here we extract the name of the variables from the xarray, and we check the type. 
                  According to it, we set the type of the variable in the hdf compressed variable and store the encoding in the fbencodings dictionary
                  Looks like fbencodings = { 'obsvalue@body': {'compression': 'gzip'}, 'sonde_type@conv': {'compression': 'gzip'}, ... } 
 
                  Empty dictionaries for hdf groups (one for each table) and their  hdf encoding description.
                  Each entry in the group dict is an xarray mapping to a single table_definition """
            
            groups={} 
            groupencodings={}
        
            """" Loop over every table_definition in the cdmd dictionary of pandas DF 
                   Two possibilities: 1. the corresponding  table already exists in the cdm (case in the final else)
                                              2. the corr table needs to be created from the local data sources (e.g. the feedback or IGRA or smt else). 
                                                  e.g. the observation_tables, the header_tables and the station_configuration.
                                                  These tables are contained in the CEUAS GitHub but not in the cdm GitHub """
            
            for k in cdmd.keys(): # loop over all the table definitions e.g. ['id_scheme', 'crs', 'station_type', 'observed_variable', 'station_configuration', 'station_configuration_codes', 'observations_table', 'header_table']
                groups[k]=xr.Dataset() # create an xarray for each table
                groupencodings[k]={} # create a dict of group econding for each table
                
                """ Loop over each row in the definition table.
                     In the CDM, the def. tables have the 4 columns ['element_name', 'kind', 'external_table', 'description' ] 
                     d is each element in the DF, so e.g. d.external_table will return the correpsonding value """    
                
                for i in range(len(cdmd[k])): 
                    d=cdmd[k].iloc[i] # extract the element at the i-th row, e.g.   element_name           z_coordinate_type, kind                                 int, external_table    z_coordinate_type:type, description         Type of z coordinate
                                                                             
                    if k in ('observations_table'):
                        try:                                                     
                            groups[k][d.element_name]=( {'hdrlen':fbds.variables['date@hdr'].shape[0] }, fromfb(fbds._variables, cdmfb[d.element_name] ) )
                        except KeyError:
                            x=numpy.zeros(fbds.variables['date@hdr'].values.shape[0], dtype= numpy.dtype(ttrans(d.kind, kinds=okinds)))
                            x.fill(numpy.nan)
                            groups[k][d.element_name]= ({'hdrlen':fbds.variables['date@hdr'].shape[0]}, x)
                            
                    elif k in ('header_table'):
                        # if the element_name is found in the cdmfb dict, then it copies the data from the odb into the header_table
                        try:
                            groups[k][d.element_name]= ({'hdrlen':fbds.variables['date@hdr'].shape[0] },
                                        fromfb(fbds._variables, cdmfb[d.element_name] ) )
                        except KeyError:
                            # if not found, it fills the columns with nans of the specified kind. Same for the observation_tables 
                            x= numpy.zeros(fbds.variables['date@hdr'].values.shape[0],dtype=numpy.dtype(ttrans(d.kind,kinds=okinds)))
                            x.fill(numpy.nan)
                            groups[k][d.element_name]= ({'hdrlen':fbds.variables['date@hdr'].shape[0]}, x)
                            
                    elif k in ('station_configuration'): # station_configurationt contains info of all the stations, so this extracts only the one line for the wanted station with the numpy.where
                        try:   
                            if 'sci' not in locals(): 
                                sci=numpy.where(cdm[k]['primary_id']=='0-20000-0-'+fbds['statid@hdr'].values[0].decode('latin1'))[0]
                            if len(sci)>0:
                                groups[k][d.element_name]=({k+'_len':1},
                                        cdm[k][d.element_name].values[sci])
                        except KeyError:
                            pass
                            
                    else : # this is the case where the cdm tables DO exist
                        try:   
                            groups[k][d.element_name]=({k+'_len':len(cdm[k])},
                                        cdm[k][d.element_name].values) # element_name is the netcdf variable name, which is the column name of the cdm table k 
                        except KeyError:
                            pass
                    try:
                        groups[k][d.element_name].attrs['external_table'] = d.external_table # defining variable attributes that point to other tables (3rd and 4th columns)
                        groups[k][d.element_name].attrs['description']      = d.description
                        #print('good element in cdm table: ',k,d.element_name)
                        groupencodings[k][d.element_name]={'compression': 'gzip'}
                    except KeyError:
                        print('bad:', k, d.element_name)
                        pass
    
            #this writes the dateindex to the netcdf file. For faster access it is written into the root group
            """ Wiriting the di (date index) as a separate xarray. """
            di.to_netcdf(fno, format='netCDF4', engine='h5netcdf', mode='w')
   
            """ Writing the content of the original odb file to the netCDF. For each variable, use the proper type encoding."""
            fbds.to_netcdf(fno, format='netCDF4', engine='h5netcdf', encoding=fbencodings,group='era5fb',mode='a')
                        
            """ Writing each separate CDM table to the netCDF file """
            for k in groups.keys():
                groups[k].to_netcdf(fno, format='netCDF4', engine='h5netcdf', encoding=groupencodings[k], group=k, mode='a') 
                
            print('sizes: in: {:6.2f} out: {:6.2f}'.format( os.path.getsize( fn) /1024/1024, os.path.getsize( fno )/1024/1024) )
            del fbds
        
    print(fno,time.time()-t)
    
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
    cdmtabledeflist=['id_scheme', 'crs', 'station_type', 'observed_variable', 'station_configuration', 'station_configuration_codes', 'observations_table', 'header_table']  
    cdm_tabdef = dict()
    for key in cdmtabledeflist:
        url='table_definitions'.join(cdmpath.split('tables'))+key+'.csv' # https://github.com/glamod/common_data_model/tree/master/table_definitions/ + ..._.dat 
        f=urllib.request.urlopen(url)
        col_names=pd.read_csv(f,delimiter='\t',quoting=3,nrows=0,comment='#')
        f=urllib.request.urlopen(url)
        tdict={col: str for col in col_names}
        cdm_tabdef[key]=pd.read_csv(f,delimiter='\t',quoting=3,dtype=tdict,na_filter=False,comment='#')
        
    
    """ Selecting the list of tables. 'station_configuration_codes','observations_table','header_table' are not implemented in the CDM GitHub"""        
    cdmtablelist=['id_scheme', 'crs', 'station_type', 'observed_variable', 'station_configuration_codes']        
    cdm_tab=dict() # dictionary where each key is the name of the cdm table, and the value is read from the .dat file    
    for key in cdmtablelist:
        f=urllib.request.urlopen(cdmpath+key+'.dat')
        col_names=pd.read_csv(f,delimiter='\t',quoting=3,nrows=0)
        f=urllib.request.urlopen(cdmpath+key+'.dat')
        tdict={col: str for col in col_names}
        cdm_tab[key]=pd.read_csv(f,delimiter='\t',quoting=3,dtype=tdict,na_filter=False)


    """ Adding the  tables that currently only have the definitions but not the implementation in the CDM, or need extensions """  
    cdm_tabdef['header_table']          = pd.read_csv(tpath+'/table_definitions/header_table.csv',delimiter='\t',quoting=3,comment='#')
    cdm_tabdef['observations_table'] = pd.read_csv(tpath+'/table_definitions/observations_table.csv',delimiter='\t',quoting=3,comment='#')

    id_scheme={cdm_tabdef['id_scheme'].element_name.values[0]:[0,1,2,3,4,5,6],
               cdm_tabdef['id_scheme'].element_name.values[1]:['WMO Identifier','Volunteer Observing Ships network code',
                                                             'WBAN Identifier','ICAO call sign','CHUAN Identifier',
                                                             'WIGOS Identifier','Specially constructed Identifier']}

    cdm_tab['id_scheme']=pd.DataFrame(id_scheme)
    #cdm['id_scheme'].to_csv(tpath+'/id_scheme_ua.dat')
    cdm_tab['crs']=pd.DataFrame({'crs':[0],'description':['wgs84']})
    #cdm['crs'].to_csv(tpath+'/crs_ua.dat')
    cdm_tab['station_type']=pd.DataFrame({'type':[0,1],'description':['Radiosonde','Pilot']})
    #cdm['station_type'].to_csv(tpath+'/station_type_ua.dat')
    #cdm['observed_variable']=pd.read_csv(tpath+'/observed_variable.dat',delimiter='\t',quoting=3,dtype=tdict,na_filter=False,comment='#')   
    
    
    return cdm_tabdef  , cdm_tab, tdict


def csvListFromUrls(url=''):
    """ Return a list of csv files, as fond in the url on the cdm GitHub """   
    urlpath = urlopen(url)
    string = urlpath.read().decode('utf-8')
    split = string.split(' ')
    csv_files_list = [m.replace('"','') for m in [n.split('title="')[1] for n in split if '.csv' in n and "title" in n] ] 
    return csv_files_list



def filelist_cleaner(lista, dataset=''):
       """ Removes unwanted files that might be present in the database directories """
       print('Cleaning the list of files to be converted')
       if dataset == 'ncar':
          cleaned = [ l for l in lista if '.nc' not in l ]
       if dataset == 'bufr':
          cleaned = [ l for l in lista if '.bfr' in l ]
       if 'era5' in dataset:
           cleaned = [ l for l in lista if '.nc' not in l and '.conv.' in l ]
       else:
           cleaned = lista
       
       return cleaned
   
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
    
    
    
    """
    parser.add_argument('--auxtables_dir' , '-a', 
                    help="Optional: path to the auxiliary tables directory. If not given, will use the files in the data/tables directory" ,
                    default = '../../cdm/data/',
                    type = str)
    parser.add_argument('--odbtables_dir' , '-odb', 
                    help="Optional: path to the odb tables directory. If not given, will use the files in the data/tables directory" ,
                    #default = '../../cdm/data/tables/',
                    default = '/raid60/scratch/leo/scratch/era5/odbs/',
                    type = str)
    parser.add_argument('--output_dir' , '-out',
                    help="Optional: path to the netcdf output directory" ,
                    default = '../../cdm/data/tables',
                    type = str)
    

    args = parser.parse_args()
    dpath = args.database_dir
    tpath = args.auxtables_dir
    output_dir = args.output_dir
    odbpath = args.odbtables_dir
    """
    
    args = parser.parse_args()
    dataset = args.dataset 
    out_dir = args.output

    if dataset not in ['era5_1', 'era5_3188', 'era5_1759', 'era5_1761', 'bufr', 'igra2', 'ncar', 'test', 'all' ]:
        raise ValueError(" The selected dataset is not valid. Please choose from ['era5_1', 'era5_1759', 'era5_1761', 'era5_3188', 'bufr', 'igra2', 'ncar', 'test', 'all' ]  ")
    
   


    """ Loading the CDM tables into pandas dataframes """
    cdm_tabdef  , cdm_tab , tdict = load_cdm_tables()
    
 
    """ Paths to the databases """
    
    examples_dir = os.getcwd() + '/examples'
    stat_conf_dir = os.getcwd() + '/stations_configurations/'   
    
    
 
    """ Sources of the files """
    db                    = { 'era5_1'       : { 'dbpath' : '/raid60/scratch/leo/scratch/era5/odbs/1'            , 'stat_conf' : 'station_configuration_era5_1.dat'       , 'example': 'era5.conv._01009'             } ,
                                   'era5_3188' : { 'dbpath' : '/raid60/scratch/leo/scratch/era5/odbs/3188'      , 'stat_conf' : 'station_configuration_era5_3188.dat'  , 'example': 'era5.3188.conv.C:6072'   } ,
                                   'era5_1759' : { 'dbpath' : '/raid60/scratch/leo/scratch/era5/odbs/1759'      , 'stat_conf' : 'station_configuration_era5_1759.dat'  , 'example': 'era5.1759.conv.6:99041' } ,
                                   'era5_1761' : { 'dbpath' : '/raid60/scratch/leo/scratch/era5/odbs/1761'      , 'stat_conf' : 'station_configuration_era5_1761.dat'  , 'example': 'era5.1761.conv.9:967'     } ,
                                   'ncar'           : { 'dbpath' : '/raid60/scratch/federico/databases/UADB'         , 'stat_conf' : 'station_configuration_ncar.dat'            , 'example': 'uadb_trhc_81405.txt'       } ,
                                   'igra2'          : { 'dbpath' : '/raid60/scratch/federico/databases/IGRAv2'      , 'stat_conf' : 'station_configuration_igra2.dat'           , 'example': 'BRM00082571-data.txt'   } ,
                                   'bufr'            : { 'dbpath' : '/raid60/scratch/leo/scratch/era5/odbs/ai_bfr'    , 'stat_conf' : 'station_configuration_bufr.dat'             , 'example': 'era5.94998.bfr'                }    }
    
    
    
    
    
    

    if dataset == 'test':
        output_dir = out_dir + '/tests'
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
            
        """  To run one example file included in the examples/ directory """
        
        print( blue + '*** Running the example files stored in ' + examples_dir + ' ***  \n \n ' + cend)


        #for s in db.keys():
        for s in ['ncar'] :
            stat_conf_file = stat_conf_dir +  db[s]['stat_conf']            
            f = examples_dir + '/' + db[s]['example']    
            
            print('Analyzing the file: *** ', f  )            
            cdm_tab['station_configuration']=pd.read_csv(stat_conf_file,  delimiter='\t', quoting=3, dtype=tdict, na_filter=False, comment='#')
            
            if 'era5' in s and 'bufr' not in s:                             
                odb_to_cdm(cdm_tab, cdm_tabdef, output_dir, f)    

            else:    
                df_to_cdm(cdm_tab, cdm_tabdef, output_dir, f)
        
        print('****** \n \n \n Finished processing the file : ', f)
        
    else:
        if dataset == 'all':
            datasets = db.keys()
        elif dataset in ['era5_1', 'era5_1759', 'era5_3188', 'era5_1761', 'bufr', 'igra2', 'ncar'] :
            datasets = [dataset]
      
        os.environ["OMP_NUM_THREADS"] = "1" # to avoid multithreading complainings

        npool = 30          
        p=Pool(npool)   
        
        print(red + '*** Running on ' + str(npool) + ' cores \n' + cend )
        for d in datasets:
            output_dir = out_dir + '/' + d
            
            if not os.path.isdir(output_dir):
                os.makedirs(output_dir)
                
            print( blue + '*** Processing the database ' + d + ' ***  \n \n ' + cend)
            
            stat_conf_file = stat_conf_dir +  db[d]['stat_conf']            
            cdm_tab['station_configuration']=pd.read_csv(stat_conf_file,  delimiter='\t', quoting=3, dtype=tdict, na_filter=False, comment='#')
            
            files_list = [ db[d]['dbpath'] + '/' + f for f in os.listdir(db[d]['dbpath']) if os.path.isfile(db[d]['dbpath'] + '/' + f)] # extracting the files list stores in the database path 
            
            files_list = [ f for f in files_list if os.path.getsize(f) > 10 ] # cleaning the list of files in the original database directories
                        
            files_list = filelist_cleaner(files_list, d) 
            
            print(green + '*** The total number of file in the dataset ' + d +  'to be processed is: ' + str(len(files_list)) + ' ***' + cend)
            print('*** Check files list: ', files_list)
             
            if 'era5' in d and 'bufr' not in d:   
                
                func= partial(odb_to_cdm, cdm_tab, cdm_tabdef, output_dir)
                transunified= p.map(func, files_list ) 
                
            else:    
                func= partial(df_to_cdm, cdm_tab, cdm_tabdef, output_dir)
                transunified= p.map(func, files_list )     
      
    print('Convertion completed !')
