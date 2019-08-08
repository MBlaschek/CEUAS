#!/usr/bin/env python
import sys
import os.path
import glob

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
#from eccodes import *
from numba import *
import matplotlib.pylab as plt
import cartopy.crs as ccrs
import argparse
import copy
from io import StringIO
import h5netcdf

"""  using fixed lenght strings (80 chars) 
Compression does not work well with normal strings 

"""
okinds={'varchar (pk)':numpy.dtype('|S80'),'varchar':numpy.dtype('|S80'),'numeric':numpy.float32,'int':numpy.int32,
       'timestamp with timezone':numpy.datetime64,
       'int[]*':list,'int[]':list,'varchar[]*':list,'varchar[]':list}


kinds={'varchar (pk)':str,'varchar':str,'numeric':numpy.float32,'int':numpy.int32,
       'timestamp with timezone':numpy.datetime64,
       'int[]*':list,'int[]':list,'varchar[]*':list,'varchar[]':list}



def make_datetime(dvar,tvar):
    """ Converts into date-time """
    dvari=dvar.values.astype(numpy.int)
    tvari=tvar.values.astype(numpy.int)
    df=pd.DataFrame({'year':dvar//10000,'month':(dvar%10000)//100,'day':dvar%100,
                        'hour':tvar//10000,'minute':(tvar%10000)//100,'second':tvar%100})
    dt=pd.to_datetime(df).values
    return dt


""" Translates the odb variables name  into cdm var """
cdmfb={'observation_value':'obsvalue@body',
       'observed_variable':'varno@body',
       'z_coordinate_type':'vertco_type@body',
       'z_coordinate':'vertco_reference_1@body',
       'date_time':[make_datetime,'date@hdr','time@hdr'],
       'longitude':'lon@hdr',
       'latitude':'lat@hdr'}

igra_cdm = { 'longitude': 'lon'}



def igra2_ascii_to_dataframe(file=''):
    """ Read the igra2 stations in ASCII format and convert to a Pandas DataFrame. 
        Adapted from https://github.com/MBlaschek/CEUAS/tree/master/CEUAS/data/igra/read.py 
       
        Args:
             file (str): path to the igra2 station file

        Returns:
             Pandas DataFrame
    """


    if not os.path.isfile(file):
        raise IOError("File not Found! %s" % file)

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

    
    read_data, dates = [], [] #  Lists containing the raw data from the ascii file, and the observation dates

    """ Data to be extracted and stored from the igra2 station files 
        Some info is contained in the header of each ascent, some in the following data """

    columns=['ident','year','month','day','hour','reltime','p_src','np_src','lat','longi',
             'lvltyp1', 'lvltyp2', 'etime', 'press', 'pflag', 'gph', 'zflag', 'temp', 'tflag', 'rh' ,'dpdep', 'wdir','wspd']
 
    ident,year,month,day,hour,reltime,p_src,np_src,lat, longi = 0,0,0,0,0,0,0,0,0,0
    lvltyp1,lvltyp2,etime,press,pflag,gph,zflag,temp,tflag,rh,dpdep,wdir,wspd = 0,0,0,0,0,0,0,0,0,0,0,0,0 # initialize to zeros

    idate = 0
    for i, line in enumerate(data):
        if line[0] == '#':
            # Info from the Header line of each ascent                                                                                                                                                                                                                   
            ident = line[1:12]               # station identifier
            year = line[13:17]               # year, months, day, hour of the observation
            month = line[18:20]
            day = line[21:23]
            hour = line[24:26]               
            reltime = line[27:31]            # release time of the sounding.
            numlev = int(line[32:36])        # number of levels in the sounding == number of data recorded in the ascent
            p_src = line[37:45]              # data source code for the pressure levels 
            np_src = line[46:54]             # data source code for non-pressure levels
            lat = int(line[55:62]) / 10000.  # latitude and longitude
            lon = int(line[63:71]) / 10000.

            """ FF TODO: check if this format is cdm compliant """
            if int(hour) == 99:
                time = reltime + '00'
            else:
                time = hour + '0000'
       
            if '99' in time:
                time = time.replace('99', '00')

            
            idate = datetime.strptime(year + month + day + time, '%Y%m%d%H%M%S')
      
        else:
           # Data of each ascent
            lvltyp1 = int(line[0])  # 1-  1   integer
            lvltyp2 = int(line[1])  # 2-  2   integer
            etime = int(line[3:8])  # 4-  8   integer
            press = int(line[9:15])  # 10- 15   integer
            pflag = line[15]  # 16- 16   character
            gph = int(line[16:21])  # 17- 21   integer
            zflag = line[21]  # 22- 22   character
            temp = int(line[22:27]) / 10.  # 23- 27   integer
            tflag = line[27]  # 28- 28   character
            rh = int(line[28:33]) / 10.  # 30- 34   integer
            dpdp = int(line[34:39]) / 10.  # 36- 40   integer
            wdir = int(line[40:45])  # 41- 45   integer
            wspd = int(line[46:51]) / 10.  # 47- 51   integer

            # raw.append((lvltyp1, lvltyp2, etime, press, pflag, gph, zflag, temp, tflag, rh, dpdp, wdir, wspd))
            #raw.append((press, gph, temp, rh, dpdp, wdir, wspd))
            #dates.append(idate)
        
        read_data.extend( (ident,year,month,day,hour,reltime,p_src,np_src,lat, longi) )
        read_data.extend( (lvltyp1,lvltyp2,etime,press,pflag,gph,zflag,temp,tflag,rh,dpdep,wdir,wspd) )
        dates.append(idate)
        
    #print(read_data, dates)
    #input('check the line')

    #df = pd.DataFrame(data= read_data, index=dates, columns= columns)
    df = pd.DataFrame(data= read_data, columns= columns)
    return df

 


    
def read_all_odbsql_stn_withfeedback(odbfile):

    #countlist=glob.glob(opath+'/*.count')
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
                    print(r[:6])
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
            print('after odb:',time.time()-t)
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

    print(odbfile,time.time()-t)
    #idy=numpy.lexsort((alldict['varno@body'],
                       #-alldict['vertco_reference_1@body'],
                       #alldict['time@hdr'],
                       #alldict['date@hdr']))
    #for k in alldict.columns:
        #alldict[k]=alldict[k][idy]

    print(odbfile,time.time()-t)

    """ may not be necessary to convert into x_array sicne you can write a pandas df into an HDF file """

    #print('FF alldict_to_xarray is: ', alldict.to_xarray )
   # input('verifica alldicts.to_xarray')

    return alldict.to_xarray()

"""
def par_read_bufr_stn_nofeedback(varno,bufrfile):

    '''       [  1.95301010e+07,   1.30000000e+04,   5.00000000e+00,                                                                                                                                                                                                                      
          3.50000000e+01,   5.59300000e+01,   3.75200000e+01,                                                                                                                                                                                                                             
          1.87000000e+02,   1.00000000e+00,   3.34000000e+04,                                                                                                                                                                                                                             
          2.00000000e+00,   2.22550003e+02]])                                                                                                                                                                                                                                             
alldict['header']                                                                                                                                                                                                                                                                         
['date', 'time', 'obstype', 'codetype', 'lat', 'lon', 'stalt', 'vertco_type', 'vertco_reference_1', 'varno', 'obsvalue']                                                                                                                                                                  
len(alldict['header'])                                                                                                                                                                                                                                                                    
11                                                                                                                                                                                                                                                                                        
    '''

    alldata=''
    alldict=dict()

    bufrlist=[]
    tx=time.time()
    try:
        f = open(bufrfile)
        cnt = 0
        # loop over the messages in the file                                                                                                                                                                                                                                              
        while 1:
            # get handle for message                                                                                                                                                                                                                                                      
            bufr = codes_bufr_new_from_file(f)
            if bufr is None:
                break
            # we need to instruct ecCodes to expand all the descriptors                                                                                                                                                                                                                   
            # i.e. unpack the data section                                                                                                                                                                                                                                                
            codes_set(bufr, 'unpack', 1)
            # get all the timePeriods                                                                                                                                                                                                                                                     
            #iterid = codes_bufr_keys_iterator_new(bufr)                                                                                                                                                                                                                                  

            # loop over the keys                                                                                                                                                                                                                                                           
            #if codes_get_array(bufr,'dataSubCategory')[0]!=101:                
            
            ##            print codes_get_array(bufr,'dataSubCategory')[0]                                                                                                                                                                                                                        
                        #codes_release(bufr)                                                                                                                                                                                                                                                      
                        #continue                                                                                                                                                                                                                                                                 
                    #while codes_bufr_keys_iterator_next(iterid):                                                                                                                                                                                                                                 
        
                        ## print key name                                                                                                                                                                                                                                                         
                        #keyname = codes_bufr_keys_iterator_get_name(iterid)                                                                                                                                                                                                                      
                        #print keyname,codes_get_array(bufr,keyname)                                                                                                                                                                                                                              
        
                    ## delete the key iterator                                                                                                                                                                                                                                                    
                    #codes_bufr_keys_iterator_delete(iterid)                                                                                                                                                                                                                                      
        
                    
            datum = float('19'+codes_get_array(bufr, "typicalDate")[0][2:])
            timePeriod = float(codes_get_array(bufr, "typicalTime")[0])
            pressure = codes_get_array(bufr, "pressure")
            #        extendedVerticalSoundingSignificance = codes_get_array(bufr, "extendedVerticalSoundingSignificance")                                                                                                                                                                         
            #        geopotentialHeight = codes_get_array(bufr, "nonCoordinateGeopotentialHeight")                                                                                                                                                                                                
            #        latitudeDisplacement = codes_get_array(bufr, "latitudeDisplacement")                                                                                                                                                                                                         
            #        longitudeDisplacement = codes_get_array(bufr, "longitudeDisplacement")                                                                                                                                                                                                       
            if varno==2:
                airTemperature = codes_get_array(bufr, "airTemperature")
            elif varno==111:
                    #dewpointTemperature = codes_get_array(bufr, "dewpointTemperature")                                                                                                                                                                                                           
                windDirection = codes_get_array(bufr, "windDirection")
                windSpeed = codes_get_array(bufr, "windSpeed")
            else:
                print('unimplemented varno',varno)
                return alldict
            if cnt==0:
                lat = codes_get(bufr, "latitude")
                lon = codes_get(bufr, "longitude")
                alt = float(codes_get(bufr, "heightOfStation"))
                blockNumber = codes_get(bufr, "blockNumber")
                stationNumber = codes_get(bufr, "stationNumber")
        
            codes_release(bufr)
            #print 'station %d%d' % (blockNumber,stationNumber)                                                                                                                                                                                                                           
        
            #print 'timePeriod pressure geopotentialHeight latitudeDisplacement longitudeDisplacement airTemperature windDirection windSpeed significance'                                                                                                                                
            miss_val=-1.e100

            if varno==2:
                for i in range(0,len(airTemperature)):
                    if airTemperature[i]!=miss_val:
                        line=numpy.asarray((datum,timePeriod,5.,35.,lat,lon,alt,1.,pressure[i],2.0,airTemperature[i]))
                        bufrlist.append(line)#print pressure[i],airTemperature[i],windDirection[i],windSpeed[i]                                                                                                                                                                           
                        cnt += 1
            else:
                miss_val=-1.e100
                for i in range(0,len(windDirection)):
                    if windSpeed[i]!=miss_val and windDirection[i]!=2147483647:
                        line=numpy.asarray((datum,timePeriod,5.,35.,lat,lon,alt,1.,pressure[i],111.0,windDirection[i]))
                        bufrlist.append(line)#print pressure[i],airTemperature[i],windDirection[i],windSpeed[i]                                                                                                                                                                           
                        line=numpy.asarray((datum,timePeriod,5.,35.,lat,lon,alt,1.,pressure[i],112.0,windSpeed[i]))
                        bufrlist.append(line)#print pressure[i],airTemperature[i],windDirection[i],windSpeed[i]                                                                                                                                                                           
                        cnt += 1
                
        f.close()
                #        print '/'.join(bufrfile.split('/')[-1:]),cnt,"messages",time.time()-tx                                                                                                                                                                                                           
    except:
       
        try:
            codes_release(bufr)
        except:
            pass
        try:
            f.close()
        except:
            pass
        return alldict
                
    if len(bufrlist)==0:
        return alldict
    alldict['header']=['date', 'time', 'obstype', 'codetype', 'lat', 'lon', 'stalt', 'vertco_type', 'vertco_reference_1',
                                       'varno', 'obsvalue']
                
    ad=numpy.asarray(bufrlist)
    h=alldict['header']
    statid=str(blockNumber*1000+stationNumber)
    alldict[statid]=dict()
    alldict[statid]=dict()
    idy=numpy.lexsort((ad[:,h.index('varno')],
                       ad[:,h.index('vertco_reference_1')],
                       ad[:,h.index('time')],
                       ad[:,h.index('date')]))
    alldict[statid]['data']=ad[idy,:]
    alldict[statid]['source']=['BUFRDATA']
    alldict[statid]['odbstatid']=[statid]
    alldict[statid]['odbfile']=bufrfile

    print('/'.join(bufrfile.split('/')[-1:]),statid,cnt,"messages",time.time()-tx)

    return alldict
    


def par_read_igra_stn_nofeedback(varno,odbfile):

    alldata=''
    alldict=dict()


    t=time.time()
    sonde_type=True
    if os.path.getsize(odbfile)>0:
        try:
            with open(odbfile) as f:
                rdata=f.read().split('\n')
                hdr=rdata[0].split()
                l=0
                idx=[]
                h=[]
                rdata.pop()
                for r in rdata:
                    if r[0]=='#':
                        idx.append(l)
                        h.append(numpy.fromstring(r[6:36]+r[54:71],sep=' ',dtype=numpy.int))
                    l+=1
            for l in idx:
                rdata[l]=' '*50
            rdata='\n'.join(rdata)
            rdata=' '.join(rdata.split('A'))
            rdata=' '.join(rdata.split('B'))
            xdata=numpy.fromstring(rdata,sep=' ',dtype=numpy.float64)
            print(xdata.shape)


        except KeyError:
            print('igra read failed!:'+' '+odbfile)
            return alldict
    else:
        return alldict
    tx=time.time()


    cols=9
    xdata=numpy.reshape(xdata,(xdata.shape[0]/cols,cols))
    ydata=numpy.empty((xdata.shape[0]*2,11),numpy.float64)
    k=fill_ydata(xdata,ydata,numpy.asarray(idx),numpy.asarray(h),varno)

    #print time.time()-tx
    alldict['header']=['date', 'time', 'obstype', 'codetype', 'lat', 'lon', 'stalt', 'vertco_type', 'vertco_reference_1', 'varno', 'obsvalue']

#    statid=stationrow[0]
    statid=hdr[0][6:12]
    alldict[statid]=dict()

    alldict[statid]['data']=ydata[:k,:]
    alldict[statid]['source']=[hdr[7]] 
    alldict[statid]['odbstatid']=[statid] 
    alldict[statid]['odbfile']=odbfile

    print(odbfile,time.time()-tx)

    return alldict
"""


            
            
def fromfb(fbv,cdmfb,cdmkind):
    """ input: 
               fbv    : feedback variable (cdm compliant)
               cdmfb  :  
               cdmkind: data type of the cdmfb
    """
    x=0
    # checks if the type of the variable is a list, so that it uses the function to extract the date time 
    if type(cdmfb) is list:
        x=cdmfb[0](fbv[cdmfb[1]],fbv[cdmfb[2]])
    else:
        if cdmfb=='varno@body':
            tr=numpy.zeros(113,dtype=int) 
            """ translate odb variables number to Lot3 numbering convention """
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

def ttrans(cdmtype,kinds=kinds):
    """ convert the cdm types to numpy types """    
    nptype=numpy.float32
    try:
        nptype=kinds[cdmtype.strip()]
    except:
        print(cdmtype,'not found, using numpy.float32')
        
    
    return nptype

@njit
def find_dateindex(y):
    """ creates the indices list from the dates, for quick access 
        nb the benchmark script will not work with these files since the definition of the array size is swapped i.e. (x.shape[0], 3)"""        


    x=numpy.unique(y)
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

def odb_to_cdm(cdm, cdmd, fn):
    """ Convert the odb file into cdm compliant netCDF files
        input:
              fn   :: odb file name (e.g. era5.conv._10393)
              cdm  :: cdm tables (read with pandas)
              cdmd :: cdm tables definitions ("") """

    recl=0    
    t=time.time()
    #f=gzip.open(fn)
    #fn=fn[:-3]
    fnl=fn.split('/')
    fnl[-1]='ch'+fnl[-1]
    fno=output_dir + '/' + fnl[-1] + '.nc' # creating an output file name e.g. chera5.conv._10393.nc  , try 01009 faster
    if not False:
        
        # era5 analysis feedback is read from compressed netcdf files era5.conv._?????.nc.gz in $RSCRATCH/era5/odbs/1
        fbds=read_all_odbsql_stn_withfeedback(fn) # this is the xarray converted from the pandas dataframe 
        #fbds=xr.open_dataset(f)
        print(time.time()-t) # to check the reading of the odb
        # the fbencodings dictionary specifies how fbds will be written to disk in the CDM compliant netCDF file.
        # float64 is often not necessary. Also int64 is often not necessary. 
        fbencodings={}
        for d in fbds._variables.keys():
            if fbds.variables[d].dtype==numpy.dtype('float64'):
                if d!='date@hdr':             
                    fbencodings[d]={'dtype':numpy.dtype('float32'),'compression': 'gzip'} # probably dtype not neccessary, but compression must be there
                else:
                    fbencodings[d]={'dtype':numpy.dtype('int32'),'compression': 'gzip'}               
            else:
                fbencodings[d]={'compression': 'gzip'}
                
        y=fbds['date@hdr'].values
        z=find_dateindex(y)
        di=xr.Dataset() 
        di['dateindex']=({'days':z.shape[1],'drange':z.shape[0]},z) # date, index of the first occurrance, index of the last 
    

        # odb is read into xarray. now we must encode the cdm into several xarray datasets
        # each cdm table is written into an hdf group, groups is the dict of all the groups
        # to write the group to the disk, you need the group encoding dict
        groups={}
        groupencodings={}
        for k in cdmd.keys(): # loop over all the table definitions 
            groups[k]=xr.Dataset() # create an  xarray
            groupencodings[k]={} # create a dict of group econding

            for i in range(len(cdmd[k])): # in the cdm table definitions you always have the element(column) name, the type, the external table and the description 
                d=cdmd[k].iloc[i] # so here loop over all the rows of the table definition . iloc is just the index of the element
                # two possibilities: 1. the corresponding  table already exists in the cdm (case in the final else)
                #                    2. the corr table needs to be created from the local data sources (e.g. the feedback or IGRA or smt else). 
                # These are the observation_tables, the header_tables and the station_configuration.
                # These tables are contained in the CEUAS GitHub but not in the cdm GitHub
                if k in ('observations_table'):
                    try:
                        # fbds is an xarray dataset , fbds._variables is a dict of the variables 
                        groups[k][d.element_name]=({'hdrlen':fbds.variables['date@hdr'].shape[0]},
                                    fromfb(fbds._variables,cdmfb[d.element_name],ttrans(d.kind,kinds=okinds)))
                    except KeyError:
                        x=numpy.zeros(fbds.variables['date@hdr'].values.shape[0],dtype=numpy.dtype(ttrans(d.kind,kinds=okinds)))
                        x.fill(numpy.nan)
                        groups[k][d.element_name]=({'hdrlen':fbds.variables['date@hdr'].shape[0]},x)
                        
                elif k in ('header_table'):
                    # if the element_name is found in the cdmfb dict, then it copies the data from the odb into the header_table
                    try:
                        groups[k][d.element_name]=({'hdrlen':fbds.variables['date@hdr'].shape[0]},
                                    fromfb(fbds._variables,cdmfb[d.element_name],ttrans(d.kind,kinds=okinds)))
                    except KeyError:
                        # if not found, it fills the columns with nans of the specified kind. Same for the observation_tables 
                        x=numpy.zeros(fbds.variables['date@hdr'].values.shape[0],dtype=numpy.dtype(ttrans(d.kind,kinds=okinds)))
                        x.fill(numpy.nan)
                        groups[k][d.element_name]=({'hdrlen':fbds.variables['date@hdr'].shape[0]},x)
                        
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
                    groups[k][d.element_name].attrs['external_table']=d.external_table # defining variable attributes that point to other tables (3rd and 4th columns)
                    groups[k][d.element_name].attrs['description']=d.description
                    #print('good element in cdm table: ',k,d.element_name)
                    groupencodings[k][d.element_name]={'compression': 'gzip'}
                except KeyError:
                    print('bad:',k,d.element_name)
                    pass

        #this writes the dateindex to the netcdf file. For faster access it is written into the root group
        di.to_netcdf(fno,format='netCDF4',engine='h5netcdf',mode='w')
        # add era5 feedback to the netcdf file. 
        # fbds is the 60 col. xarray
        fbds.to_netcdf(fno,format='netCDF4',engine='h5netcdf',encoding=fbencodings,group='era5fb',mode='a')
        
        for k in groups.keys():
            #this code deletes some variables to check how well they are compressed. Variable strings are badly compressed
            #gk=list(groups[k].keys())
            #for l in gk:
                #if groups[k][l].dtype==numpy.dtype('<U1'):
                    #del groups[k][l]
                    #del groupencodings[k][l]
            
            #this appends group by group to the netcdf file
            groups[k].to_netcdf(fno,format='netCDF4',engine='h5netcdf',encoding=groupencodings[k],group=k,mode='a') #
        print('sizes: in: {:6.2f} out: {:6.2f}'.format(os.path.getsize(fn)/1024/1024,
                                              os.path.getsize(fno)/1024/1024))
        del fbds
    

    # speed test for accessing dateindex
    #for k in range(10):
        
        #t=time.time()
        #with h5py.File(fno,'r') as f:
            #di=f['dateindex'][:]
            #idx=numpy.where(di[0,:]==19950101)[0]
            #print(numpy.nanmean(f['observations_table']['observation_value'][di[1,idx[0]]:di[2,idx[0]]+1]))
        #print('ch',time.time()-t)
        #fng='cgera5'.join(fno.split('chera5'))
        #t=time.time()
        #with h5py.File(fng,'r') as g:
            #do=g['dateindex'][:]
            #idx=numpy.where(do[:,2]==19950101)[0]
            #print(numpy.nanmean(g['obsvalue@body'][do[idx[0],0]:do[idx[0],1]+1]))
        #print('cg',time.time()-t)
        
        
    # first implementation of inner join  - resolving numeric variable code 
    #with h5py.File(fno,'r') as f:
        #ext=f['observations_table']['observed_variable'].attrs['external_table']
        #lext=ext.split(':')
        
        #l=[]
        #lidx=[]
        #llen=len(f['observations_table']['observed_variable'][:])
        #for i in range(llen):
            #obv=f['observations_table']['observed_variable'][i]
            #if obv not in l:
                #idx=numpy.where(f[lext[0]][lext[1]][:]==obv)[0][0]
                #l.append(obv)
                #lidx.append(idx)
            #else:
                #idx=lidx[l.index(obv)]
            #for k in f[lext[0]].keys():
                #print(k,':',f[lext[0]][k][idx],)
            #print('value:',f['observations_table']['observation_value'][i])
        #print(f)

    
    print(fno,time.time()-t)

    return recl


def csvListFromUrls(url=''):
    """ Return a list of csv files, as fond in the url on the cdm GitHub """   
    urlpath = urlopen(url)
    string = urlpath.read().decode('utf-8')
    split = string.split(' ')
    csv_files_list = [m.replace('"','') for m in [n.split('title="')[1] for n in split if '.csv' in n and "title" in n] ] 
    return csv_files_list







if __name__ == '__main__':


    parser = argparse.ArgumentParser(description="Make CDM compliant netCDFs")
    parser.add_argument('--database_dir' , '-d', 
                    help="Optional: path to the database directory. If not given, will use the files in the data directory" ,
                    default = '../../cdm/data/',
                    type = str)
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
    print ('THE DPATH IS', dpath)
    if not dpath:
        dpath = '../cdm/code/data/'
   
    if not tpath:
        tpath = dpath+'/tables/'
   
    """
    print ('Analysing the databases: ')
    print ('dpath is ' , dpath)
    print ('tpath is ' , tpath )


    cdmpath='https://raw.githubusercontent.com/glamod/common_data_model/master/tables/' # cdm tables            
    

    # TODO: get the list of files in the tables
    cdmtablelist=['id_scheme','crs','station_type','observed_variable','station_configuration_codes']        
    cdm_tab=dict() # dictionary where each key is the name of the cdm table, and the value is read from the .dat file

    """ e.g. 
    'id_scheme':    scheme                            description
    0       0                               WIGOS ID
    1       1                               GRUAN ID
    2       2                             IMO Number
    3       3                            National ID
    4       4              WMO buoy / station number
    5       5               Ship / platform callsign
    6       6       Generic ID (e.g. SHIP, PLAT etc)
    7       7                           Station name
    8       8                           ICOADS other
    9       9                         ICOADS unknown
    10     10                       ICOADS composite
    11     11  Oceangraphic platform / cruise number
    12     12          Other buoy number (e.g. Argo)
    13     13                C3S 311a Lot 2 Internal,
    """ 
    for key in cdmtablelist:
        f=urllib.request.urlopen(cdmpath+key+'.dat')
        col_names=pd.read_csv(f,delimiter='\t',quoting=3,nrows=0)
        f=urllib.request.urlopen(cdmpath+key+'.dat')
        tdict={col: str for col in col_names}
        cdm_tab[key]=pd.read_csv(f,delimiter='\t',quoting=3,dtype=tdict,na_filter=False)


    #print('the cdm dictionary is:', cdm )
    #input('FF vedi')
    cdm_tabdef=dict()


    """ # Uncomment to get the list of all the .csv files present at the url specified
    url = 'https://github.com/glamod/common_data_model/tree/master/table_definitions'
    cdmtabledeflist = csvListFromUrls(url)
    """

    # see that there are more entries than in the rpevious list, since e.g. station_configuration does not exist int he cdm GitHub but was created in the CEUAS GitHub 
    cdmtabledeflist=['id_scheme','crs','station_type','observed_variable','station_configuration','station_configuration_codes','observations_table','header_table']        
    for key in cdmtabledeflist:
        url='table_definitions'.join(cdmpath.split('tables'))+key+'.csv' # https://github.com/glamod/common_data_model/tree/master/table_definitions/ + ..._.dat 
        f=urllib.request.urlopen(url)
        col_names=pd.read_csv(f,delimiter='\t',quoting=3,nrows=0,comment='#')
        f=urllib.request.urlopen(url)
        tdict={col: str for col in col_names}
        cdm_tabdef[key]=pd.read_csv(f,delimiter='\t',quoting=3,dtype=tdict,na_filter=False,comment='#')


    # up to here: only information read from the public cdm github
    # header table is instead read from CEUAS github, cdm directory 
    #cdmd['header_table']=pd.read_csv(tpath+'../table_definitions/header_table.csv',delimiter='\t',quoting=3,comment='#')


    #print('FF cdmd dictionary', cdm_tabdef)
   # input('check')
    cdm_tabdef['header_table']=pd.read_csv(tpath+'/table_definitions/header_table.csv',delimiter='\t',quoting=3,comment='#')
    cdm_tabdef['observations_table']=pd.read_csv(tpath+'/table_definitions/observations_table.csv',delimiter='\t',quoting=3,comment='#')

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



    igra = 'igra2_example/BRM00082571-data.txt'
    a = igra2_ascii_to_dataframe(file= igra)


    """ Up to here, we have two different dictionaries.
    cdm_tab: contains the dictionary for the cdm tables
    cdm_tabdef: contains the dictionary for the cdm tables definitions
    """





    ''' 

    #func=partial(odb_to_cdm, cdm=cdm_tab, cdmd=cdm_tabdef)
    func=partial(odb_to_cdm, cdm_tab, cdm_tabdef)

    #bfunc=partial(read_bufr_stn_meta,2)
    #rfunc=partial(read_rda_meta) 
    tu=dict()
    p=Pool(25)
    
    dbs=['1759'] #,'igra2','ai_bfr','rda','3188','1759','1761']
    
    debug = True
    
    for odir in dbs: 
        print('*** Processing a sample file from the odb dataset: ', odir)
        #cdm_tab['station_configuration']=pd.read_csv(odbpath+'/'+odir+'/station_configuration.dat',delimiter='\t',quoting=3,dtype=tdict,na_filter=False,comment='#') #fix this for the cases of 3188,1759,1761 . This file does not exist for the other db
        if debug:
            cdm_tab['station_configuration']= pd.read_csv('/raid60/scratch/leo/scratch/era5/odbs/1/station_configuration.dat',delimiter='\t',quoting=3,dtype=tdict,na_filter=False,comment='#')
        
        """
        if 'ai' in odir:
            pass
        elif 'rda' in odir:
            pass
        elif 'igra2' in odir:
            pass

        else:
            #flist=glob.glob(odbpath+odir+'/'+'era5.conv.*01009.nc.gz')
            #flist=glob.glob(odbpath+odir+'/'+'era5.conv._01009')
            flist = glob.glob(odbpath+odir+'/'+'era5.3188.conv.C:6072')
            transunified=list(map(func,flist))
        """
        if odir == '1':
            flist=glob.glob(odbpath+odir+'/'+'era5.conv.*01009.nc.gz')
        elif odir =='3188':
            flist=glob.glob(odbpath+odir+'/'+'era5.3188.conv.C:6072')
        elif odir =='1759':
            flist=glob.glob(odbpath+odir+'/'+'era5.1759.conv.6:99041')
        elif odir =='1761':
            flist=glob.glob(odbpath+odir+'/'+'era5.1761.conv.9:967')
            
        transunified=list(map(func,flist))
        
    '''
