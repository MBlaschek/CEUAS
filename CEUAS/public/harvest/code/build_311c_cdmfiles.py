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
from rasotools.utils import *
from eccodes import *
import matplotlib.pylab as plt
import cartopy.crs as ccrs
import argparse

okinds={'varchar (pk)':numpy.dtype('|S80'),'varchar':numpy.dtype('|S80'),'numeric':numpy.float32,'int':numpy.int32,
       'timestamp with timezone':numpy.datetime64,
       'int[]*':list,'int[]':list,'varchar[]*':list,'varchar[]':list}
kinds={'varchar (pk)':str,'varchar':str,'numeric':numpy.float32,'int':numpy.int32,
       'timestamp with timezone':numpy.datetime64,
       'int[]*':list,'int[]':list,'varchar[]*':list,'varchar[]':list}

def make_datetime(dvar,tvar):
    dvari=dvar.values.astype(numpy.int)
    tvari=tvar.values.astype(numpy.int)
    df=pd.DataFrame({'year':dvar//10000,'month':(dvar%10000)//100,'day':dvar%100,
                        'hour':tvar//10000,'minute':(tvar%10000)//100,'second':tvar%100})
    dt=pd.to_datetime(df).values
    return dt

cdmfb={'observation_value':'obsvalue@body',
       'observed_variable':'varno@body',
       'z_coordinate_type':'vertco_type@body',
       'z_coordinate':'vertco_reference_1@body',
       'date_time':[make_datetime,'date@hdr','time@hdr'],
       'longitude':'lon@hdr',
       'latitude':'lat@hdr'}
    
def fromfb(fbv,cdmfb,cdmkind):
    x=0
    if type(cdmfb) is list:
        x=cdmfb[0](fbv[cdmfb[1]],fbv[cdmfb[2]])
    else:
        if cdmfb=='varno@body':
            tr=numpy.zeros(113,dtype=int)
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
            
            x=tr[fbv[cdmfb].values.astype(int)]
        else:    
            x=fbv[cdmfb].values
        
    return x
def ttrans(cdmtype,kinds=kinds):
    
    nptype=numpy.float32
    try:
        nptype=kinds[cdmtype.strip()]
    except:
        print(cdmtype,'not found, using numpy.float32')
        
    
    return nptype

@njit
def find_dateindex(y):
    
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

def odb_to_cdm(cdm,cdmd,fn):
    recl=0
    
    t=time.time()
    f=gzip.open(fn)
    fn=fn[:-3]
    fnl=fn.split('/')
    fnl[-1]='ch'+fnl[-1]
    fno='/'.join(fnl)
    if not False:
        
        fbds=xr.open_dataset(f)
        #ds=xr.open_dataset(fn,engine='h5netcdf')
        print(time.time()-t)
        fbencodings={}
        for d in fbds._variables.keys():
            if fbds.variables[d].dtype==numpy.dtype('float64'):
                if d!='date@hdr':             
                    fbencodings[d]={'dtype':numpy.dtype('float32'),'compression': 'gzip'}
                else:
                    fbencodings[d]={'dtype':numpy.dtype('int32'),'compression': 'gzip'}               
            else:
                fbencodings[d]={'compression': 'gzip'}
                
        y=fbds['date@hdr'].values
        z=find_dateindex(y)
        di=xr.Dataset()
        di['dateindex']=({'days':z.shape[1],'drange':z.shape[0]},z)
        #obds=xr.Dataset()
        #obencodings={}
        #for i in range(len(cdmd['observations_table'])):
            #d=cdmd['observations_table'].iloc[i]
            #obds[d.element_name]=({'hdrlen':fbds.variables['date@hdr'].shape[0]},
                                #numpy.zeros_like(fbds.variables['date@hdr'].values,dtype=numpy.dtype(ttrans(d.kind))))
            #obds[d.element_name].attrs['external_table']=d.external_table
            #obds[d.element_name].attrs['description']=d.description
            #obencodings[d.element_name]={'compression': 'gzip'}
    
        groups={}
        groupencodings={}
        for k in cdmd.keys():
            groups[k]=xr.Dataset()
            groupencodings[k]={}
            for i in range(len(cdmd[k])):
                d=cdmd[k].iloc[i]
                if k in ('observations_table'):
                    try:
                        groups[k][d.element_name]=({'hdrlen':fbds.variables['date@hdr'].shape[0]},
                                    fromfb(fbds._variables,cdmfb[d.element_name],ttrans(d.kind,kinds=okinds)))
                    except KeyError:
                        x=numpy.zeros(fbds.variables['date@hdr'].values.shape[0],dtype=numpy.dtype(ttrans(d.kind,kinds=okinds)))
                        x.fill(numpy.nan)
                        groups[k][d.element_name]=({'hdrlen':fbds.variables['date@hdr'].shape[0]},x)
                        
                elif k in ('header_table'):
                    try:
                        groups[k][d.element_name]=({'hdrlen':fbds.variables['date@hdr'].shape[0]},
                                    fromfb(fbds._variables,cdmfb[d.element_name],ttrans(d.kind,kinds=okinds)))
                    except KeyError:
                        x=numpy.zeros(fbds.variables['date@hdr'].values.shape[0],dtype=numpy.dtype(ttrans(d.kind,kinds=okinds)))
                        x.fill(numpy.nan)
                        groups[k][d.element_name]=({'hdrlen':fbds.variables['date@hdr'].shape[0]},x)
                        
                elif k in ('station_configuration'):
                    try:   
                        groups[k][d.element_name]=({k+'_len':len(cdm[k])},
                                    cdm[k][d.element_name].values)#,dtype=numpy.dtype(ttrans(d.kind)))
                    except KeyError:
                        pass
                        
                else:
                    try:   
                        groups[k][d.element_name]=({k+'_len':len(cdm[k])},
                                    cdm[k][d.element_name].values)#,dtype=numpy.dtype(ttrans(d.kind)))
                    except KeyError:
                        pass
                try:
                    groups[k][d.element_name].attrs['external_table']=d.external_table
                    groups[k][d.element_name].attrs['description']=d.description
                    print('good:',k,d.element_name)
                    groupencodings[k][d.element_name]={'compression': 'gzip'}
                except KeyError:
                    print('bad:',k,d.element_name)
                    pass
    
        di.to_netcdf(fno,format='netCDF4',engine='h5netcdf',mode='w')
        fbds.to_netcdf(fno,format='netCDF4',engine='h5netcdf',encoding=fbencodings,group='era5fb',mode='a')
        
        for k in groups.keys():
            #gk=list(groups[k].keys())
            #for l in gk:
                #if groups[k][l].dtype==numpy.dtype('<U1'):
                    #del groups[k][l]
                    #del groupencodings[k][l]
            
            groups[k].to_netcdf(fno,format='netCDF4',engine='h5netcdf',encoding=groupencodings[k],group=k,mode='a') #
        print('sizes: in: {:6.2f} out: {:6.2f}'.format(os.path.getsize(fn+'.gz')/1024/1024,
                                              os.path.getsize(fno)/1024/1024))
        del fbds
    
    #del obds
    
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

if __name__ == '__main__':


    parser = argparse.ArgumentParser(description="Make CDM compliant netCDFs")
    parser.add_argument('--database_dir' , '-d', 
                    help="Optional: path to the database directory. If not given, will use the files in the data directory" ,
                    default = '../data/',
                    type = str)
    parser.add_argument('--auxtables_dir' , '-a', 
                    help="Optional: path to the auxiliary tables directory. If not given, will use the files in the data/tables directory" ,
                    default = '../data/tables/',
                    type = str)
    args = parser.parse_args()
    dpath = args.database_dir
    tpath = args.auxtables_dir

    print ('THE DPATH IS', dpath)
    if not dpath:
        dpath = '../data/'
   
    if not tpath:
        tpath = '../data/tables/'
   
    print ('Analysing the databases stored in ', dpath)
    cdmpath='https://raw.githubusercontent.com/glamod/common_data_model/master/tables/'                                                                                                                                                                                               



    cdmtablelist=['id_scheme','crs','station_type','observed_variable','station_configuration_codes']        
    cdm=dict()
    for key in cdmtablelist:
        f=urllib.request.urlopen(cdmpath+key+'.dat')
        col_names=pd.read_csv(f,delimiter='\t',quoting=3,nrows=0)
        f=urllib.request.urlopen(cdmpath+key+'.dat')
        tdict={col: str for col in col_names}
        cdm[key]=pd.read_csv(f,delimiter='\t',quoting=3,dtype=tdict,na_filter=False)

    cdmd=dict()
    cdmtabledeflist=['id_scheme','crs','station_type','observed_variable','station_configuration','station_configuration_codes','observations_table','header_table']        
    for key in cdmtabledeflist:
        url='table_definitions'.join(cdmpath.split('tables'))+key+'.csv'
        f=urllib.request.urlopen(url)
        col_names=pd.read_csv(f,delimiter='\t',quoting=3,nrows=0,comment='#')
        f=urllib.request.urlopen(url)
        tdict={col: str for col in col_names}
        cdmd[key]=pd.read_csv(f,delimiter='\t',quoting=3,dtype=tdict,na_filter=False,comment='#')
    cdmd['header_table']=pd.read_csv(tpath+'header_table.csv',delimiter='\t',quoting=3,comment='#')
    cdmd['observations_table']=pd.read_csv(tpath+'observations_table.csv',delimiter='\t',quoting=3,comment='#')
    id_scheme={cdmd['id_scheme'].element_name.values[0]:[0,1,2,3,4,5,6],
               cdmd['id_scheme'].element_name.values[1]:['WMO Identifier','Volunteer Observing Ships network code',
                                                             'WBAN Identifier','ICAO call sign','CHUAN Identifier',
                                                             'WIGOS Identifier','Specially constructed Identifier']}

    cdm['id_scheme']=pd.DataFrame(id_scheme)
    cdm['id_scheme'].to_csv(tpath+'/id_scheme_ua.dat')
    cdm['crs']=pd.DataFrame({'crs':[0],'description':['wgs84']})
    cdm['crs'].to_csv(tpath+'/crs_ua.dat')
    cdm['station_type']=pd.DataFrame({'type':[0,1],'description':['Radiosonde','Pilot']})
    cdm['station_type'].to_csv(tpath+'/station_type_ua.dat')
    cdm['observed_variable']=pd.read_csv(tpath+'/observed_variable.dat',delimiter='\t',quoting=3,dtype=tdict,na_filter=False,comment='#')


    func=partial(odb_to_cdm,cdm,cdmd)
    #bfunc=partial(read_bufr_stn_meta,2)
    #rfunc=partial(read_rda_meta) 
    tu=dict()
    p=Pool(25)
    dbs=['1','igra2','ai_bfr','rda','3188','1759','1761']
    for odir in dbs: 
        cdm['station_configuration']=pd.read_csv(dpath+'/'+odir+'/station_configuration.dat',delimiter='\t',quoting=3,dtype=tdict,na_filter=False,comment='#')
        if 'ai' in odir:
            pass
        elif 'rda' in odir:
            pass
        elif 'igra2' in odir:
            pass

        else:
            flist=glob.glob(dpath+odir+'/'+'era5.conv.*01009.nc.gz')
            transunified=list(map(func,flist))

