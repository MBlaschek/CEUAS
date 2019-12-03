#!/usr/bin/env python
import sys
import os.path
import glob
import psutil

import subprocess
import urllib.request
import xarray as xr
import numpy
#import h5pickle as h5py
import h5py
import netCDF4
from datetime import date, datetime,timedelta
import time
from multiprocessing import Pool
#from netCDF4 import Dataset
import gzip
import pandas as pd
from functools import partial
#from rasotools.utils import *
#from eccodes import *
from numba import njit
#import matplotlib.pylab as plt
#import cartopy.crs as ccrs
import argparse
import copy
import matplotlib.pylab as plt
from collections import OrderedDict
#from io import StringIO
#import h5netcdf

"""  using fixed lenght strings (80 chars) 
Compression does not work well with normal strings 

"""
okinds={'varchar (pk)':numpy.dtype('|S8'),'varchar':numpy.dtype('|S8'),'numeric':numpy.float32,'int':numpy.int32,
        'timestamp with timezone':numpy.datetime64,
        'int[]*':list,'int[]':list,'varchar[]*':list,'varchar[]':list}

gkinds={'varchar (pk)':numpy.dtype('|S80'),'varchar':numpy.dtype('|S80'),'numeric':numpy.float32,'int':numpy.int32,
        'timestamp with timezone':numpy.datetime64,
        'int[]*':list,'int[]':list,'varchar[]*':list,'varchar[]':list}


kinds={'varchar (pk)':str,'varchar':str,'numeric':numpy.float32,'int':numpy.int32,
       'timestamp with timezone':numpy.datetime64,
       'int[]*':list,'int[]':list,'varchar[]*':list,'varchar[]':list}



def substitute_highascii(cdm):

    subs={'o':[240,242,243,244,245,246,248],'O':[210,211,212,213,214,216],
          'a':[224,225,226,227,228,229,230],'A':[192,193,194,195,196,197,198],
              'u':[249,250,251,252,253],'U':[217,218,219,220],
              'i':[236,237,238,239],'I':[204,205,206,207,304],
              'S':[350],'n':[241],'c':[231],'C':[199],'e':[232,233,234,235],'E':[200,201,202,203]}
    for k in cdm['station_configuration'].columns:
        #print(k,type(cdm['station_configuration'][k][0]))
        if type(cdm['station_configuration'][k][0]) is str:
            #for l in range(cdm['station_configuration'][k].values.shape[0]):

                #cdm['station_configuration'][k].values[l]=numpy.string_(cdm['station_configuration'][k].values[l])
            try:

                cdm['station_configuration'][k].values[:]=cdm['station_configuration'][k].values[:].astype('S')
            except:
                for l in range(cdm['station_configuration'][k].values.shape[0]):
                    try:

                        cdm['station_configuration'][k].values[l]=numpy.string_(cdm['station_configuration'][k].values[l])
                    except:
                        for m in range(len(cdm['station_configuration'][k].values[l])):
                            mychar=cdm['station_configuration'][k].values[l][m]
                            if ord(mychar)>128:
                                for n,v in subs.items():
                                    #print(n,v,ord(mychar))
                                    if ord(mychar) in v:
                                        cdm['station_configuration'][k].values[l]=n.join(cdm['station_configuration'][k].values[l].split(mychar))

                        cdm['station_configuration'][k].values[l]=numpy.string_(cdm['station_configuration'][k].values[l])
                cdm['station_configuration'][k]=numpy.string_(cdm['station_configuration'][k])



    #print('converted station_configuration')

def make_datetime(dvar,tvar):
    """ Converts into dat-time """
    dvari=dvar.astype(numpy.int)
    tvari=tvar.astype(numpy.int)
    df=pd.DataFrame({'year':dvar//10000,'month':(dvar%10000)//100,'day':dvar%100,
                     'hour':tvar//10000,'minute':(tvar%10000)//100,'second':tvar%100})
    dt=pd.to_datetime(df).values
    return dt

def make_obsid(ivar):
    x=numpy.char.zfill(numpy.arange(ivar.shape[0]).astype('S8'),8)
    return x

def make_recid(ivar):
    x=numpy.char.zfill(numpy.arange(ivar.values.shape[0]).astype('S5'),5)
    return x

def make_obsrecid(fbvar,ivar):
    x=numpy.char.zfill(numpy.arange(ivar.values.shape[0]).astype('S5'),5)
    y=numpy.zeros(fbvar.shape[0]).astype('S5')
    for i in range(ivar.values.shape[0]-1):
        y[ivar.values[i]:ivar.values[i+1]]=x[i]
    if ivar.values.shape[0]>1:
        y[ivar.values[i+1]:]=x[i+1]
    return y

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
       'primary_station_id':'statid@hdr'}

def read_all_odbsql_stn_withfeedback(odbfile):

    #countlist=glob.glob(opath+'/*.count')
    alldata=''

    alldict={} #xr.Dataset()
    t=time.time()
    sonde_type=True
    obstype=True
    if os.path.getsize(odbfile)>0:
        """ Read first the odbb header to extract the column names and type """
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
                        elif 'INTEGER' in kinds[-1] or 'BITFIELD' in kinds[-1]:
                            print(columns[-1])
                            if columns[-1]=='sonde_type@conv':
                                tdict[columns[-1]]=numpy.float32
                            else: 
                                tdict[columns[-1]]=numpy.int32
                        else:
                            tdict[columns[-1]]=numpy.dtype('S') # dict containng column name and type


                except IndexError:
                    pass
        except:
            print('could not read odbfile '+odbfile)
            return alldict
        try:
            #rdata=subprocess.check_output(["odb","sql","-q","select *","-i",odbfile,'--no_alignment']) # after reading th eheader it does the query
            # returns a byte string
            #print('after odb:',time.time()-t)
            t=time.time()
            #rdata=''.join(rdata.decode('latin-1').split("'")) # decoding the stirng into a unicode
            #f=StringIO(rdata) # access the string like a file, return a file pointer to read the string with pandas

            try:
                f=gzip.open(odbfile+'.gz') 
            except:
                print(odbfile+'.gz','not found')
                return

            # access the string like a file, return a file pointer to read the string with pandas
            # nb  if you have null values, reading of integer fails and are read as floats
            # to improve, you can convert the columns with nabs in the alldicts (pandas data frame) into int(np.nans)
            # date values are large so float 32 precision is not sufficient 
            print('vor decode',odbfile+'.gz')
            #try:
            #try:
            for c in ['sensor@hdr','vertco_reference_2@body','ppcode@conv_body','timeslot@timeslot_index']:
                del tdict[c]
                columns.remove(c)

            alldict=pd.read_csv(f,delimiter='\t',usecols=columns,quoting=3,comment='#',
                                skipinitialspace=True,dtype=tdict)#,nrows=100000)
            #except:
                ##os.remove(odbfile+'.gz')
                ##print('removed',odbfile+'.gz')
                #return
            #alldict.rename(str.strip,axis='columns',inplace=True)
            for c in alldict.columns:
                if type(alldict[c][0]) in [str,bytes]:
                    l=alldict[c].shape[0]
                    slen=len(alldict[c].values[0])
                    #alldict[c]=numpy.array(alldict.pop(c).values,dtype='S{}'.format(slen))
                    alldict[c]=numpy.string_(alldict[c])
                if type(alldict[c][0]) is numpy.int64:
                    alldict[c]=numpy.int32(alldict[c])
                if type(alldict[c][0]) is numpy.float64:
                    alldict[c]=numpy.float32(alldict[c])
            #for k in tdict.keys():
                #alldict[k]=[]
            #for line in f:
                #l=line.split(b'\t')
                #i=0
                #for k in tdict.keys():
                    #alldict[k].append(l[i])
                    #i+=1
                #l=line.split()
            #rdata=f.read().split()
            print('after odb:',time.time()-t)
            #xdata=numpy.fromstring(rdata,sep='\t')
            #cl=len(tdict.keys())
            #rl=len(rdata)//cl

            #i=0
            #for k,v in tdict.items():
                #alldict[k]=numpy.array(rdata[i+cl:-cl:cl],dtype=v)
                #i+=1
                ##print(i)
            #print(time.time()-t)
            #print('after read')
            #except:
                #pass

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

    return alldict

def fromfb(fbv,di,cdmfb,cdmkind):
    """ input: 
               fbv    : feedback variable (cdm compliant)
               di     : record index variable (cdm compliant)
               cdmfb  :  
               cdmkind: data type of the cdmfb
    """
    x=0
    # checks if the type of the variable is a list, so that it uses the function to extract the date time 
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
    """ input: 
               fbv    : feedback variable (cdm compliant)
               cdmfb  :  
               cdmkind: data type of the cdmfb
    """
    x=0
    # checks if the type of the variable is a list, so that it uses the function to extract the date time 
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

def ttrans(cdmtype,kinds=kinds):
    """ convert the cdm types to numpy types """    
    nptype=numpy.float32
    try:
        nptype=kinds[cdmtype.strip()]
    except:
        print(cdmtype,'not found, using numpy.float32')


    return nptype

@njit(cache=True)
def find_dateindex(y,x):
    """ creates the indices list from the dates, for quick access 
        nb the benchmark script will not work with these files since the definition of the array size is swapped i.e. (x.shape[0], 3)"""        


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
def find_recordindex(y,x):
    """ creates the indices list from the dates, for quick access 
        nb the benchmark script will not work with these files since the definition of the array size is swapped i.e. (x.shape[0], 3)"""        


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

def write_dict_h5(dfile,f,k,fbencodings,var_selection=[],mode='a'): # cuts vars and copies attributes of observation, feedback and header tables

    #g=h5py.File('9'.join(dfile.split('7')))

    with h5py.File(dfile,mode) as fd:
        try:
            fd.create_group(k)
        except:
            pass
        #f[fbencodings[v]['dims']]=numpy.arange(f[v].shape[0],dtype=numpy.int)
        if not var_selection:
            var_selection=list(f.keys())

        index=numpy.arange(f['obsvalue@body'].shape[0],dtype=numpy.int)
        string10=numpy.zeros(80,dtype='S1')
        fd[k].create_dataset('index',data=index)
        sdict={}
        slist=[]
        #var_selection.append('index')
        #var_selection.append('string10')
        for v in var_selection:
            print(v,f[v].values.ndim,f[v].values.dtype)#,f[v].values[0].dtype)
            if type(f[v].values[0]) not in [str,bytes,numpy.bytes_]:
                if f[v].values.dtype!='S1':

                    fd[k].create_dataset(v,f[v].values.shape,f[v].values.dtype,compression=fbencodings[v]['compression'],chunks=True)
                    fd[k][v][:]=f[v]
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
                    fd[k].create_dataset('string{}'.format(slen),data=string10[:slen])

                x=x.reshape(f[v].values.shape[0],slen)
                #x[:,sleno:]=' '
                fd[k].create_dataset(v,data=x,compression=fbencodings[v]['compression'],chunks=True)
                #fd[k].create_dataset(v,data=x,compression=fbencodings[v]['compression'],chunks=True)
                #fd[k].create_dataset(v,x.shape,x.dtype,compression=fbencodings[v]['compression'],chunks=True)
                #fd[k][v][:]=x[:]

        #for v in fd[k].keys(): #var_selection:
            #for a in g[k][v].attrs.keys():
                #if a not in ['DIMENSION_LIST','CLASS']:
                    #fd[k][v].attrs[a]=g[k][v].attrs[a]
                    #if 'string' in v:
                    #print (v,a,type(g[k][v].attrs[a]))
                    #x=0
        #for v in var_selection:
            #l=0
            #for d in g[k][v].dims:
                #if len(d)>0:
                    #print(k,v,g[k][v].dims[l][0].name)
                    #fd[k][v].dims[l].attach_scale(fd[k][g[k][v].dims[l][0].name])
                #l+=1

        for v in fd[k].keys(): #var_selection:
            l=0
            try:
                if 'string' not in v and v!='index':

                    fd[k][v].dims[l].attach_scale(fd[k]['index'])
                    if type(f[v].values[0]) in [str,bytes,numpy.bytes_]:
                        slen=sdict[v]
                        #slen=10
                        fd[k][v].dims[1].attach_scale(fd[k]['string{}'.format(slen)])
                        print('')
            except:
                pass

            #for d in f[v]['dims']:
                #if len(d)>0:
                    #print(k,v,f[v]['dims'][l]['name'])
                    #fd[k][v].dims[l].attach_scale(fd[k][f[v]['dims'][l]['name']])
                #l+=1

        i=4        
        for v in slist:
            s='string{}'.format(v)
            for a in ['NAME']:
                fd[k][s].attrs[a]=numpy.bytes_('This is a netCDF dimension but not a netCDF variable.')

            for a in ['_Netcdf4Dimid']:           
                fd[k][s].attrs[a]=numpy.int64(i)
            i+=1

    #g.close()
    return
        #for v in fd[k].keys():
        #for a in fd[k][v].attrs.keys():
            #print (v,a,fd[k][v].attrs[a])
        #print('ready')

def odb_to_cdm(cdm,cdmd,fn):
    """ input:
              fn: odb file name (e.g. era5.conv._10393)
              cdm: cdm tables (read with pandas)
              cdmd: cdm tables definitions ("") """


    recl=0

    process = psutil.Process(os.getpid())
    print(process.memory_info().rss/1024/1024)        
    t=time.time()
    #f=gzip.open(fn)
    #fn=fn[:-3]
    fnl=fn.split('/')
    fnl[-1]='ch'+fnl[-1]
    fno=output_dir + '/' + fnl[-1] + '.nc' # creating an output file name e.g. chera5.conv._10393.nc  , try 01009 faster
    if not False:

        # era5 analysis feedback is read from compressed netcdf files era5.conv._?????.nc.gz in $RSCRATCH/era5/odbs/1
        fbds=read_all_odbsql_stn_withfeedback(fn) # i.e. the xarray 
        if fbds is None:
            return
        #fbds=xr.open_dataset(f)
        print(time.time()-t) # to check the reading of the odb
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
                print(d,type(v[0]),v.size)
                if type(v[0])==bytes:
                    fbencodings[d]={'compression': 'gzip','chunksizes':(min([10000,v.shape[0]]),10)}#,'chunksizes':(10000,10)
                else:
                    fbencodings[d]={'compression': 'gzip'}
        fbencodings['index']={'compression': 'gzip'}
        y=fbds['date@hdr'].values*100+fbds['time@hdr'].values//10000
        #dt=pd.DataFrame({'year': y//10000,
                #'month': (y%10000)//100,
                #'day': y%100,
                #'hour': fbds['time@hdr']//10000})

        #y=pd.to_datetime(dt).values
        tt=time.time()
        idx=numpy.lexsort((fbds['vertco_reference_1@body'].values,y))
        y=y[idx]
        for fb in fbds.keys():
            fbds[fb].values[:]=fbds[fb].values[idx]
        print(time.time()-tt)
        x=numpy.unique(y)
        z=find_recordindex(y,x)
        di=xr.Dataset() 
        di['recordindex']=({'record':z.shape[1]},z[1])
        di['recordtimestamp']=({'record':z.shape[1]},x)

        y=fbds['date@hdr'].values
        x=numpy.unique(y)
        z=find_dateindex(y,x)
        di['dateindex']=({'days':z.shape[1],'drange':z.shape[0]},z) # date, index of the first occurrance, index of the last


        #this writes the dateindex to the netcdf file. For faster access it is written into the root group
        di.to_netcdf(fno,format='netCDF4',engine='h5netcdf',mode='w')
        #sid=fbds['statid@hdr'][0].split(b"'")[1]
        #sid=sid.decode('latin1')
        #for k in cdm['station_configuration']['primary_id'].values:
            #if sid in k[-5:]:
            #print(k)

        write_dict_h5(fno, fbds, 'era5fb', fbencodings, var_selection=[],mode='a')
        # add era5 feedback to the netcdf file. 
        # fbds is the 60 col. xarray
        #for k in fbds.keys():
            #if fbds[k].values.dtype=='object':
            #print(k,len(fbds[k].values[0]))
        #fbds.to_netcdf(fno,format='netCDF4',engine='h5netcdf',encoding=fbencodings,group='era5fb',mode='a')
        print(process.memory_info().rss/1024/1024)        
        print(time.time()-t)

        #import gc
        #for obj in gc.get_objects():   # Browse through ALL objects
            #if isinstance(obj, h5py.File):   # Just HDF5 files
            #try:
                #obj.close()
            #except:
                #pass        # odb is read into xarray. now we must encode the cdm into several xarray datasets
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
                        if d.element_name=='report_id':

                            groups[k][d.element_name]=({'obslen':fbds['date@hdr'].shape[0]},
                                                       fromfb(fbds,di._variables,cdmfb[k+'.'+d.element_name],ttrans(d.kind,kinds=okinds)))
                        else:
                            groups[k][d.element_name]=({'obslen':fbds['date@hdr'].shape[0]},
                                                       fromfb(fbds,di._variables,cdmfb[d.element_name],ttrans(d.kind,kinds=okinds)))
                    except KeyError:
                        x=numpy.zeros(fbds['date@hdr'].shape[0],dtype=numpy.dtype(ttrans(d.kind,kinds=okinds)))
                        x.fill(numpy.nan)
                        groups[k][d.element_name]=({'obslen':fbds['date@hdr'].shape[0]},x)

                elif k in ('header_table'):
                    # if the element_name is found in the cdmfb dict, then it copies the data from the odb into the header_table
                    try:
                        if d.element_name=='report_id':
                            groups[k][d.element_name]=({'hdrlen':di['recordindex'].shape[0]},
                                                       hdrfromfb(fbds,di._variables,cdmfb[k+'.'+d.element_name],ttrans(d.kind,kinds=gkinds)))
                        else:
                            groups[k][d.element_name]=({'hdrlen':di['recordindex'].shape[0]},
                                                       hdrfromfb(fbds,di._variables,cdmfb[d.element_name],ttrans(d.kind,kinds=gkinds)))
                        print('got',k,d.element_name)
                        j=0

                    except KeyError:

                        # if not found, it fills the columns with nans of the specified kind. Same for the observation_tables 
                        print('fail',k,d.element_name)
                        if d.element_name in cdm['station_configuration'].columns:
                            x=numpy.zeros(di['recordindex'].shape[0],dtype=numpy.dtype(ttrans(d.kind,kinds=gkinds)))
                            try:

                                idx=numpy.where('0-20000-0-'+fnl[-1].split('_')[-1] == cdm['station_configuration']['primary_id'])[0][0]
                                groups[k][d.element_name]=x.fill(cdm['station_configuration'][d.element_name][idx])
                            except:
                                groups[k][d.element_name]=x
                        else:
                            x=numpy.zeros(di['recordindex'].shape[0],dtype=numpy.dtype(ttrans(d.kind,kinds=okinds)))
                            x.fill(numpy.nan)
                        groups[k][d.element_name]=({'hdrlen':di['recordindex'].shape[0]},x)

                elif k in ('station_configuration'): # station_configurationt contains info of all the stations, so this extracts only the one line for the wanted station with the numpy.where
                    try:   
                        if 'sci' not in locals(): 
                            sci=numpy.where(cdm[k]['primary_id']==numpy.string_('0-20000-0-'+fbds['statid@hdr'][0][1:-1].decode('latin1')))[0]
                        if len(sci)>0:
                            groups[k][d.element_name]=({k+'_len':1},
                                                       cdm[k][d.element_name].values[sci])
                            print('statconf:',k,groups[k][d.element_name])
                    except KeyError:
                        print('x')
                        pass

                else : # this is the case where the cdm tables DO exist
                    if 'expver' in d.element_name:
                        print('x')
                    try:   
                        groups[k][d.element_name]=({k+'_len':len(cdm[k])},
                                                   cdm[k][d.element_name].values) # element_name is the netcdf variable name, which is the column name of the cdm table k 
                    except KeyError:
                        pass
                try:
                    groups[k][d.element_name].attrs['external_table']=d.external_table # defining variable attributes that point to other tables (3rd and 4th columns)
                    groups[k][d.element_name].attrs['description']=d.description
                    print('good:',k,d.element_name)
                    if type(groups[k][d.element_name].values[0])==str:
                        s=groups[k][d.element_name].values.shape
                        groupencodings[k][d.element_name]={'dtype':numpy.dtype('S80'),'compression': 'gzip','chunksizes':(min(100000,s[0]),80)}
                    else:
                        groupencodings[k][d.element_name]={'compression': 'gzip'}
                except:
                    print('bad:',k,d.element_name)
                    pass


        print(process.memory_info().rss/1024/1024)        
        for k in groups.keys():
            #this code deletes some variables to check how well they are compressed. Variable strings are badly compressed
            #gk=list(groups[k].keys())
            #for l in gk:
                #if groups[k][l].dtype==numpy.dtype('<U1'):
                #del groups[k][l]
                #del groupencodings[k][l]

            #this appends group by group to the netcdf file
            print('write group',k)
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

# simple selector: choose longest record
def select_record(sdict,foundrecs):
    selectedrec=foundrecs['ilist'][0]
    for f in foundrecs['ilist']:
        #t=sdict['ifpointers'][f[0]]['observations_table']['date_time']
        #dates = netCDF4.num2date(t[f[2]], t.attrs['units'], t.attrs['calendar'])
        #print('selected:',sdict['ifiles'][f[0]], dates)
        
        if f[3]-f[2]>selectedrec[3]-selectedrec[2]:
            selectedrec=f
            
    #t=sdict['ifpointers'][selectedrec[0]]['observations_table']['date_time']
    #dates = netCDF4.num2date(t[selectedrec[2]], t.attrs['units'], t.attrs['calendar'])
    #print('selected:',sdict['ifiles'][selectedrec[0]], dates)
            
    return f


def do_merge(fd,f,sdict,mergedts,k,idx,cut_dimension,var_selection=[]): # cuts vars and copies attributes of observation, feedback and header tables
    if not var_selection:
        var_selection=f[k].keys()
    findices=[{'tindex':[],'sindex':[]} for i in range(len(sdict['ifpointers']))]
    l=0
    # append one element per value for header table
    # append many elements per value for observations and era5fb table
    if idx[0]==1:
        
        for v in mergedts.values():
            findices[v['selected'][0]]['tindex'].append(l)
            findices[v['selected'][0]]['sindex'].append(v['selected'][1])
            l+=1
    else:
        for v in mergedts.values():
            dv=v['selected'][idx[1]]-v['selected'][idx[0]]
            findices[v['selected'][0]]['tindex']+=list(numpy.arange(l,l+dv))
            findices[v['selected'][0]]['sindex']+=list(numpy.arange(v['selected'][idx[0]],v['selected'][idx[1]]))
            
            l+=dv
    shp0=l        
    for v in var_selection:
        print(k,v,f[k][v].shape)
        if f[k][v].ndim==1:
            if f[k][v].dtype=='S1':
                hilf=numpy.zeros(shape=f[k][v].shape,dtype=sdict['ifpointers'][1][k][v][0].dtype)
            else:
                hilf=numpy.zeros(shape=(shp0,),dtype=sdict['ifpointers'][1][k][v][0].dtype)
        else:
            hilf=numpy.zeros(shape=(shp0,sdict['ifpointers'][1][k][v].shape[1]),
                                 dtype=sdict['ifpointers'][1][k][v][0].dtype)
                
        l=0
        for ff in sdict['ifpointers']:
            hilf2=ff[k][v][:]
            if 'calendar' in ff[k][v].attrs.keys():
                try:
                    
                    start=datetime.strptime(ff[k][v].attrs['units'][11:].lstrip(), '%Y-%m-%d %H:%M:%S')
                    startdate=list(mergedts.keys())[0]
                    mstart=datetime(startdate//1000000,(startdate%1000000)//10000,(startdate%10000)//100,hour=startdate%100)
                    mstarts=mstart.strftime('hours since %Y-%m-%d %H:%M:%S')
                    startdelta=start-mstart
                    deltahours=startdelta.days*24+startdelta.seconds//3600
                    hilf2=hilf2+deltahours
                    if not numpy.isnan(hilf2[0]):
                        print(k,v,'start:',start,mstart,hilf2[0])
                    #unitstring='hours since {}-{:02d}-{:02d} 00:00:00'.format(startdate//1000000,(startdate%1000000)//10000,(startdate%10000)//100)
                    #hilf2 = netCDF4.date2num(hilf2, unitstring, ff[k][v].attrs['calendar']).astype(numpy.int32)
                except:
                    print('date',k,v,'could not be converted')
                    pass
                
            if hilf.dtype=='S1' and hilf2.dtype!='S1': # BUFR files have this problem
                hilf2=hilf2.astype('S{}'.format(hilf.shape[1])).view('S1').reshape(hilf2.shape[0],hilf.shape[1])
            if ff[k][v].ndim==1 and ff[k][v].dtype=='S1':
                hilf=hilf2
            else:
                if v=='source_id':
                    hilf3=hilf2[findices[l]['sindex'][:]].view('S80')
                    hilf3[:,0]=numpy.string_(sdict['ifiles'][l])
                    hilf2[findices[l]['sindex'][:]]=hilf3.view('S1')
                    hilf[findices[l]['tindex'][:]]=hilf2[findices[l]['sindex'][:]]
                else:
                    hilf[findices[l]['tindex'][:]]=hilf2[findices[l]['sindex'][:]]
            l+=1
        
        if f[k][v].ndim==1:
            if f[k][v].dtype!='S1': 
                
                fd[k].create_dataset_like(v,f[k][v],shape=(shp0,),chunks=True)
                #hilf=f[k][v][idx[0]:idx[-1]+1]
                fd[k][v][:]=hilf[:] #[idx-idx[0]]
            else:
                if v in [cut_dimension]:
                    
                    fd[k].create_dataset_like(v,f[k][v],shape=(shp0,),chunks=True)
                    #hilf=f[k][v][idx[0]:idx[-1]+1]                                                
                    fd[k][v][:]=numpy.zeros(shp0,dtype='S1') #hilf[:] #[idx-idx[0]]
                else:
                    fd[k].create_dataset_like(v,f[k][v])
                    fd[k][v][:]=f[k][v][:]
                    pass
        else:
            fd[k].create_dataset_like(v,f[k][v],shape=(shp0,f[k][v].shape[1]),chunks=True)
            #hilf=f[k][v][idx[0]:idx[-1]+1,:]
            fd[k][v][:]=hilf[:] #[idx-idx[0],:]
        for a in f[k][v].attrs.keys():
            if a not in ['DIMENSION_LIST','CLASS','external_table']:
                fd[k][v].attrs[a]=f[k][v].attrs[a]
                if v=='date_time' and a=='units':
                    fd[k][v].attrs[a]=mstarts
    for v in var_selection:
        l=0
        for d in f[k][v].dims:
            if len(d)>0:
                print(k,v,f[k][v].dims[l][0].name)
                fd[k][v].dims[l].attach_scale(fd[k][f[k][v].dims[l][0].name])
            l+=1



def write_merged(wigid,mergedts,sdict):
    
    f=sdict['ifpointers'][1] # era5_1
    with h5py.File('merged_'+wigid+'.nc','w') as fd:
        for k in f.keys():
            if isinstance(f[k],h5py.Group):
                t=time.time()
                fd.create_group(k)
                if k in ['era5fb']: # only obs, feedback fitting criteria (idx) is copied       
                    #do_merge(fd,f,sdict,mergedts,k,idx,'obslen')  #,var_selection=['index','string10','vertco_reference_1@body','obsvalue@body','fg_depar@body','statid@hdr'])
                    print('x')
                elif k in ['observations_table']:
                    do_merge(fd,f,sdict,mergedts,k,[2,3],'hdrlen')#,var_selection=['obslen','string8','observation_id','latitude','longitude','z_coordinate','observed_variable','observation_value','units','date_time'])
                    print('x')
                elif k in ['header_table']:  # only records fitting criteria (zidx) are copied
                    do_merge(fd,f,sdict,mergedts,k,[1],'hdrlen')#,var_selection=['hdrlen','string8','station_name'])

                else: # groups that are simply copied
                    #print(k)
                    for v in f[k].keys():
                        f.copy(f[k][v],fd[k],name=v,without_attrs=True)
                    for a in f[k].attrs.keys():
                        fd[k].attrs[a]=f[k].attrs[a]
                    for v in f[k].keys():
                        #print(k,v)
                        #fd[k].create_dataset_like(v,f[k][v])
                        for a in f[k][v].attrs.keys():
                            print(k,v,a)
                            if a not in ['DIMENSION_LIST','CLASS']:
                                fd[k][v].attrs[a]=f[k][v].attrs[a]
                        l=0
                        for d in f[k][v].dims:
                            if len(d)>0:
                                fd[k][v].dims[l].attach_scale(fd[k][f[k][v].dims[l][0].name])
                            l+=1
                print('do_copy',k,time.time()-t)
                
            else:  # add vars in root group - recordindex and recordtimestamp need to be recreated
                #print(k)
                fd.create_dataset_like(k,f[k])
                fd[k][:]=f[k][:]
                for a in f[k].attrs.keys():
                        if a not in ['DIMENSION_LIST','CLASS']:
                            fd[k].attrs[a]=f[k].attrs[a]
        for k in f.keys():  # attach dims in root group
            if not isinstance(f[k],h5py.Group):
                l=0
                for d in f[k].dims:
                    if len(d)>0:
                        fd[k].dims[l].attach_scale(fd[f[k].dims[l][0].name])
                    l+=1
        
        print(wigid)
        
        
    

if __name__ == '__main__':


    parser = argparse.ArgumentParser(description="Make CDM compliant netCDFs")
    parser.add_argument('--database_dir' , '-d', 
                        help="Optional: path to the database directory. If not given, will use the files in the data directory" ,
                        default = '../../../cdm/data/',
                    type = str)
    parser.add_argument('--auxtables_dir' , '-a', 
                        help="Optional: path to the auxiliary tables directory. If not given, will use the files in the data/tables directory" ,
                        default = '../../../cdm/data/tables/',
                    type = str)
    parser.add_argument('--odbtables_dir' , '-odb', 
                        help="Optional: path to the odb tables directory. If not given, will use the files in the data/tables directory" ,
                        default = '../data/tables/',
                    type = str)

    parser.add_argument('--output_dir' , '-out',
                        help="Optional: path to the netcdf output directory" ,
                        default = '../data/tables',
                    type = str)


    args = parser.parse_args()
    dpath = args.database_dir
    tpath = args.auxtables_dir
    output_dir = args.output_dir
    odbpath = args.odbtables_dir

    print ('THE DPATH IS', dpath)
    if not dpath:
        dpath = '../../cdm/code/data/'

    if not tpath:
        tpath = dpath+'/tables/'

    print ('Analysing the databases stored in ', dpath)
    cdmpath='https://raw.githubusercontent.com/glamod/common_data_model/master/tables/'                                                                                                                                                                                               


    # TODO: get the list of files in the tables directory of the cdm 
    cdmtablelist=['id_scheme','crs','station_type','observed_variable','station_configuration_codes','units']        
    cdm=dict()
    for key in cdmtablelist:
        f=urllib.request.urlopen(cdmpath+key+'.dat')
        col_names=pd.read_csv(f,delimiter='\t',quoting=3,nrows=0)
        f=urllib.request.urlopen(cdmpath+key+'.dat')
        tdict={col: str for col in col_names}
        cdm[key]=pd.read_csv(f,delimiter='\t',quoting=3,dtype=tdict,na_filter=False)

    cdmd=dict()

    # TODO: get the list of files in the tables_definition directory
    # you need the list of file in the directory, there must be some urllib function for that
    # see that there are more entries than in the rpevious list, since e.g. station_configuration does not exist int he cdm GitHub but was created in the CEUAS GitHub 
    cdmtabledeflist=['id_scheme','crs','station_type','observed_variable','station_configuration','station_configuration_codes',
                     'observations_table','header_table','units']        
    for key in cdmtabledeflist:
        url='table_definitions'.join(cdmpath.split('tables'))+key+'.csv'
        f=urllib.request.urlopen(url)
        col_names=pd.read_csv(f,delimiter='\t',quoting=3,nrows=0,comment='#')
        f=urllib.request.urlopen(url)
        tdict={col: str for col in col_names}
        cdmd[key]=pd.read_csv(f,delimiter='\t',quoting=3,dtype=tdict,na_filter=False,comment='#')

    # up to here: only information read from the public cdm github


    # header table is instead read from CEUAS github, cdm directory 
    cdmd['header_table']=pd.read_csv(tpath+'../table_definitions/header_table.csv',delimiter='\t',quoting=3,comment='#')
    cdmd['observations_table']=pd.read_csv(tpath+'../table_definitions/observations_table.csv',delimiter='\t',quoting=3,comment='#')
    id_scheme={cdmd['id_scheme'].element_name.values[0]:[0,1,2,3,4,5,6],
               cdmd['id_scheme'].element_name.values[1]:['WMO Identifier','Volunteer Observing Ships network code',
                                                         'WBAN Identifier','ICAO call sign','CHUAN Identifier',
                                                         'WIGOS Identifier','Specially constructed Identifier']}

    cdm['id_scheme']=pd.DataFrame(id_scheme)
    #cdm['id_scheme'].to_csv(tpath+'/id_scheme_ua.dat')
    cdm['crs']=pd.DataFrame({'crs':[0],'description':['wgs84']})
    #cdm['crs'].to_csv(tpath+'/crs_ua.dat')
    cdm['station_type']=pd.DataFrame({'type':[0,1],'description':['Radiosonde','Pilot']})

    #cdm['station_configuration']=pd.read_csv(os.path.expandvars('$RSCRATCH/era5/odbs/1/meta.csv'),delimiter='\t',quoting=3,comment='#')
    #cdm['station_type'].to_csv(tpath+'/station_type_ua.dat')
    #cdm['observed_variable']=pd.read_csv(tpath+'/observed_variable.dat',delimiter='\t',quoting=3,dtype=tdict,na_filter=False,comment='#')


    func=partial(odb_to_cdm,cdm,cdmd)
    #bfunc=partial(read_bufr_stn_meta,2)
    #rfunc=partial(read_rda_meta) 
    tu=dict()
    p=Pool(25)
    dbs=['1','igra2','ai_bfr','rda','3188','1759','1761']
    sids=set()
    for odir in dbs:
        cdm[odir]={}
        cdm[odir]['station_configuration']=pd.read_csv(odbpath+'/'+odir+'/meta.csv',delimiter='\t',quoting=3,dtype=tdict,na_filter=False,comment='#')
        #cdm[odir]['station_configuration']=pd.read_csv(odbpath+'/'+odir+'/station_configuration.dat',delimiter='\t',quoting=3,dtype=tdict,na_filter=False,comment='#')
        substitute_highascii(cdm[odir])
        sids=sids|set(cdm[odir]['station_configuration']['primary_id'])

    lids=[]
    for s in sids:
        if b'20000' in s:
            lids.append(s)
    lids.sort()
    
    tt=time.time()
    os.chdir('example_stations')
    slist=['82930']
    #idx=lids.index('0-20000-0-'+slist[-1])
    sdict={}
    for s in slist:
        wigid='0-20000-0-'+s
        sdict[wigid]={}
        sdict[wigid]['ifiles']=glob.glob('*/*'+s+'*.nc')
        sdict[wigid]['ifpointers']=[]
        sdict[wigid]['recordindex']=[]
        sdict[wigid]['recordtimestamp']=[]
        for ifile in sdict[wigid]['ifiles']:
            sdict[wigid]['ifpointers'].append(h5py.File(ifile))
            cfp=sdict[wigid]['ifpointers'][-1]
            print(ifile,cfp.keys())
            if 'recordindex' in cfp.keys():
                sdict[wigid]['recordindex'].append(cfp['recordindex'][:])
                sdict[wigid]['recordtimestamp'].append(cfp['recordtimestamp'][:])
            else:
        # this is slow, can be replaced with h5py read once recordindex is in the source files
                t=sdict[wigid]['ifpointers'][-1]['observations_table']['date_time']
                dates = netCDF4.num2date(t[:], t.attrs['units'], t.attrs['calendar'])
                y=[]
                for d in dates:
                    y.append(d.year*1000000+d.month*10000+d.day*100+d.hour)
                y=numpy.array(y,dtype=numpy.int32)
                # this should not be necessary and is incomplete here. All columns would need to be sorted. 
                idx=numpy.lexsort((cfp['observations_table']['z_coordinate'][:],y))
                y=y[idx]
                
                x=numpy.unique(y)
                z=find_recordindex(y,x)
    
                sdict[wigid]['recordindex'].append(z[1])
                sdict[wigid]['recordtimestamp'].append(x)

        mergedts={}
        l=0
        for ifile in sdict[wigid]['ifiles']:
            m=0
            rlen=len(sdict[wigid]['recordtimestamp'][l])
            for rec in sdict[wigid]['recordtimestamp'][l]:
                if m+1<rlen:
                    endrec=sdict[wigid]['recordindex'][l][m+1]
                else:
                    endrec=len(sdict[wigid]['ifpointers'][l]['observations_table']['date_time'])
                if rec in mergedts.keys():
                    mergedts[rec]['ilist'].append((l,m,sdict[wigid]['recordindex'][l][m],endrec))
                else:
                    mergedts[rec]={'ilist':[(l,m,sdict[wigid]['recordindex'][l][m],endrec)]} # file index, recordindex index, recordindex 
                m+=1
            l+=1
        
        mergedts=dict(sorted(mergedts.items()))
        for rec,v in mergedts.items():
            mergedts[rec]['selected']=select_record(sdict[wigid],v)
            print(rec, mergedts[rec]['selected'][0])
            
        write_merged(wigid,mergedts,sdict[wigid])
            
        with h5py.File('merged_'+wigid+'.nc') as fd:
            x=fd['observations_table']['observed_variable'][:]
            y=fd['observations_table']['observation_value'][:]
            z=fd['observations_table']['z_coordinate'][:]
            y[y<-800]=numpy.nan
            mask=numpy.ones_like(x,dtype=numpy.bool)
            mask=numpy.logical_and(mask,x==107)
            mask=numpy.logical_and(mask,z==50000)
            t=fd['observations_table']['date_time'][:]
            plt.plot(t[mask],y[mask])
            
            l=0
            for check in [b'bufr',b'era5_1/',b'era5_1759',b'igra',b'ncar']:
                
                hilf=fd['observations_table']['source_id'][:].view('S80')[:,0]
                idxb=numpy.flatnonzero(numpy.core.defchararray.find(hilf,check)!=-1)
                fdstarts=fd['observations_table']['date_time'].attrs['units']
                bstarts=sdict[wigid]['ifpointers'][l]['observations_table']['date_time'].attrs['units']
                fdstart=datetime.strptime(fdstarts[11:].lstrip(), '%Y-%m-%d %H:%M:%S')
                bstart=datetime.strptime(bstarts[11:].lstrip(), '%Y-%m-%d %H:%M:%S')
                startdelta=bstart-fdstart
                deltahours=startdelta.days*24+startdelta.seconds//3600
                tb=sdict[wigid]['ifpointers'][l]['observations_table']['date_time'][:]+deltahours
                idb=numpy.where(tb==t[idxb[0]])
                idd=numpy.where(t==t[idxb[0]])
                
                print ('source '+check.decode()+':',sdict[wigid]['ifpointers'][l]['observations_table']['observation_value'][:][idb])
                print ('merged  early   :',fd['observations_table']['observation_value'][:][idd])
                idb=numpy.where(tb==t[idxb[-1]])
                idd=numpy.where(t==t[idxb[-1]])
                
                print ('source '+check.decode()+':',sdict[wigid]['ifpointers'][l]['observations_table']['observation_value'][:][idb])
                print ('merged late    :',fd['observations_table']['observation_value'][:][idd])
                l+=1
            
            process = psutil.Process(os.getpid())
            print(process.memory_info().rss/1024/1024, time.time()-tt)        
            print('x')
            
        for f in sdict[wigid]['ifpointers']:
            f.close()
            
            
            
            
        

    #transunified=list(map(func,hlist))
    print(time.time()-tt)
    print('finished')

