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
import numpy as np
 
from eccodes import * 


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

            datum = float('19'+codes_get_array(bufr, "typicalDate")[0][2:])
            timePeriod = float(codes_get_array(bufr, "typicalTime")[0])
            pressure = codes_get_array(bufr, "pressure")

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



bufr_example = 'bufr_example/era5.94998.bfr'

varno = 2
alldic = par_read_bufr_stn_nofeedback(varno, bufr_example)

print('all dict is', alldic)
