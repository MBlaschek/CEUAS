import os,sys
import h5py #ickle as h5py
import numpy
from numba import *
#from numba.typed import List
import xarray
import pandas as pd
import copy
import time
import glob
from multiprocessing import Pool
from functools import partial

import xml.etree.ElementTree as ET
import urllib.request
import json
import datetime

#@njit(cache=True)
def read_standardnames():

    try:
        with open('cf.json') as f:
            cf=json.load(f)
            
    except:
        
        url = 'http://cfconventions.org/Data/cf-standard-names/69/src/cf-standard-name-table.xml'
        response = urllib.request.urlopen(url).read()
        tree = ET.fromstring(response)
        snames=['platform_id','platform_name','latitude','longitude','time','air_pressure',
                'air_temperature','dew_point_temperature','relative_humidity','specific_humidity',
                'eastward_wind','northward_wind','wind_speed','wind_direction','geopotential','trajectory_label',
                'obs_minus_bg','obs_minus_an','bias_estimate']
        cdmnames=['header_table/primary_station_id','header_table/station_name','observations_table/latitude',
                  'observations_table/longitude','observations_table/date_time','observations_table/z_coordinate']
        cdmnames+=9*['observations_table/observation_value']
        cdmnames+=['header_table/report_id','era5fb/fg_depar@body','era5fb/an_depar@body','era5fb/biascorr@body']
        cf={}
        for c,cdm in zip(snames,cdmnames):
            cf[c]={'cdmname':cdm,'units':'NA','shortname':c}
            if c not in 'latitude longitude time air_pressure':
                cf[c]['coordinates']='lat lon time plev' # short names
        l=0
        for child in tree:
            #print(child.tag,child.attrib)
            try:
                c=child.attrib['id']
                if c in snames:
                    i=snames.index(c)
                    print(c)
                    
                    cf[c]['cdmname']=cdmnames[i]
                    if child[0].text is not None:
                        cf[c]['units']=child[0].text
                    if child[2].text is not None:
                        cf[c]['shortname']=child[2].text
                    cf[c]['standard_name']=c
            except:
                pass
            l+=1
        cf['latitude']['shortname']='lat'
        cf['longitude']['shortname']='lon'
        cf['air_pressure']['shortname']='plev'
        cf['time']['shortname']='time'
        cf['bias_estimate']['cdsname']='bias_estimate'
        cf['bias_estimate']['cdmcode']=0
        cf['bias_estimate']['odbcode']=0
        cf['obs_minus_bg']['cdsname']='obs_minus_bg'
        cf['obs_minus_bg']['cdmcode']=0
        cf['obs_minus_bg']['odbcode']=0
        cf['obs_minus_an']['cdsname']='obs_minus_an'
        cf['obs_minus_an']['cdmcode']=0
        cf['obs_minus_an']['odbcode']=0
        cf['air_temperature']['cdsname']='temperature'
        cf['air_temperature']['cdmcode']=85
        cf['air_temperature']['odbcode']=2
        cf['eastward_wind']['cdsname']='u_component_of_wind'
        cf['eastward_wind']['cdmcode']=104
        cf['eastward_wind']['odbcode']=3
        cf['northward_wind']['cdsname']='v_component_of_wind'
        cf['northward_wind']['cdmcode']=105
        cf['northward_wind']['odbcode']=4
        cf['wind_speed']['cdsname']='wind_speed'
        cf['wind_speed']['cdmcode']=107
        cf['wind_speed']['odbcode']=112
        cf['wind_direction']['cdsname']='wind_direction'
        cf['wind_direction']['cdmcode']=106
        cf['wind_direction']['odbcode']=111
        cf['relative_humidity']['cdsname']='relative_humidity'
        cf['relative_humidity']['cdmcode']=38
        cf['relative_humidity']['odbcode']=29
        cf['specific_humidity']['cdsname']='specific_humidity'
        cf['specific_humidity']['cdmcode']=39
        cf['specific_humidity']['odbcode']=7
        cf['dew_point_temperature']['cdsname']='dew_point_temperature'
        cf['dew_point_temperature']['cdmcode']=36
        cf['dew_point_temperature']['odbcode']=59
        cf['geopotential']['cdsname']='geopotential'
        cf['geopotential']['cdmcode']=-1
        cf['geopotential']['odbcode']=1
       
     #vdict={'111':'windDirection','112':'windSpeed','1':'geopotentialHeight',
           #'2':'airTemperature','59':'dewpointTemperature','29':'relativeHumidity'}

        with open('cf.json','w') as f:
            
            json.dump(cf,f)
            
    return cf
        
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

@njit
def find_dateindex_cg(y):
    
    x=numpy.unique(y)
    z=numpy.zeros((x.shape[0],3),dtype=numpy.int32)
    z-=1
    j=0
    for i in range(len(y)):
        m=y[i]
        if x[j]==y[i]:
            if z[j,0]==-1:
                z[j,0]=i
                #print(j,i)
            else:
                if z[j,1]<i:
                    z[j,1]=i
        elif x[j]<y[i]:
            j+=1
            if x[j]==y[i]:
                if z[j,0]==-1:
                    z[j,0]=i
                    #print(j,i)
                else:
                    if z[j,1]<i:
                        z[j,1]=i
            else:
                print('Error')
        else:
            j-=1
            if x[j]==y[i]:
                if z[j,0]==-1:
                    z[j,0]=i
                    #print(j,i)
                else:
                    if z[j,1]<i:
                        z[j,1]=i
            else:
                print('Error')
    z[:,2]=x
    return z

def do_copy(fd,f,k,idx,cut_dimension,var_selection=[]): # cuts vars and copies attributes of observation, feedback and header tables
    if not var_selection:
        var_selection=f[k].keys()
    for v in var_selection:
        if f[k][v].ndim==1:
            if f[k][v].dtype!='S1':
                
                fd[k].create_dataset_like(v,f[k][v],shape=idx.shape,chunks=True)
                hilf=f[k][v][idx[0]:idx[-1]+1]
                fd[k][v][:]=hilf[idx-idx[0]]
            else:
                if v in [cut_dimension]:
                    
                    fd[k].create_dataset_like(v,f[k][v],shape=idx.shape,chunks=True)
                    hilf=f[k][v][idx[0]:idx[-1]+1]                                                
                    fd[k][v][:]=hilf[idx-idx[0]]
                else:
                    fd[k].create_dataset_like(v,f[k][v])
                    fd[k][v][:]=f[k][v][:]
                    pass
        else:
            fd[k].create_dataset_like(v,f[k][v],shape=(idx.shape[0],f[k][v].shape[1]),chunks=True)
            hilf=f[k][v][idx[0]:idx[-1]+1,:]
            fd[k][v][:]=hilf[idx-idx[0],:]
        for a in f[k][v].attrs.keys():
            if a not in ['DIMENSION_LIST','CLASS','external_table']:
                fd[k][v].attrs[a]=f[k][v].attrs[a]
    for v in var_selection:
        l=0
        for d in f[k][v].dims:
            if len(d)>0:
                print(k,v,f[k][v].dims[l][0].name)
                fd[k][v].dims[l].attach_scale(fd[k][f[k][v].dims[l][0].name])
            l+=1

def do_cfcopy(fd,f,k,idx,cf,dim0,var_selection=[]): # cuts vars and copies attributes of observation, feedback and header tables
    
    tt=time.time()
    if not var_selection:
        var_selection=f[k].keys()
    
    vlist=[]    
    clist=[]
    for i in cf.keys():
        if i not in ['platform_id', 'platform_name']: #, 'latitude', 'longitude', 'time', 'air_pressure']:
            clist.append(i)
            if i in ['air_temperature','dew_point_temperature','relative_humidity','specific_humidity',
                'eastward_wind','northward_wind','wind_speed','wind_direction','geopotential']:
                for fb in ['obs_minus_bg','obs_minus_an','bias_estimate']:
                    try:
                        cf[fb]['units']=cf[i]['units']
                        cf[fb]['standard_name']=i
                        cf[fb]['long_name']=k.split('fb')[0].upper()+ ' reanalysis '+fb
                    except:
                        pass
                    
    for cfk,cfv in cf.items():
        for v in var_selection:
            #print(k+'/'+v,cfv['cdmname'])
            if k+'/'+v==cfv['cdmname']:
                vlist.append(cfv['shortname'])
            
                if f[k][v].ndim==1:
                    try:               
                        fd.create_dataset_like(vlist[-1],f[k][v],shape=idx.shape,chunks=True)
                        hilf=f[k][v][idx[0]:idx[-1]+1]
                        if 'time' in v: # convert time units 
                            us=f[k][v].attrs['units']
                            #dh=us.split(' ')[-1].split(':')
                            
                            if b'hours' in us:
                                hilf=hilf*3600 # hilf+=int(dh[0])
                            elif b'minutes' in us:
                                hilf=hilf*60  #+int(dh[0])
                            elif b'seconds' in us:
                                hilf=hilf #//60//60+int(dh[0])
                            elif b'days' in us:
                                hilf*=24*3600

                        fd[vlist[-1]][:]=hilf[idx-idx[0]]
                    except:
                        print(k,v)
                        pass
                else:
                    s1=f[k][v].shape[1]
                    fd.create_dataset_like(vlist[-1],f[k][v],shape=(idx.shape[0],s1),chunks=True)
                    #if k=='header_table':
                        #print(k,v,time.time()-tt, 'nach create')       
                    #sname=f[k][v].dims[1][0].name.split('/')[-1]
                    sname='string{}'.format(s1)
                    #if k=='header_table':
                        #print(k,v,time.time()-tt)       
                    if sname not in fd.keys():
                        #if k=='header_table':
                            #print(k,v,time.time()-tt)       
                        fd.create_dataset(sname,data=numpy.zeros(s1,dtype='S1'),chunks=True)
                        fd[sname].attrs['NAME']=numpy.string_('This is a netCDF dimension but not a netCDF variable.')
        
                    #if k=='header_table':
                        #print(k,v,time.time()-tt)       
                    hilf=f[k][v][idx[0]:idx[-1]+1,:]
                    if hilf.shape[0]==0:
                        print('x')
                    fd[vlist[-1]][:]=hilf[idx-idx[0],:]
                    
                #if k=='header_table':
                    #print(k,v,time.time()-tt)       
                for a in f[k][v].attrs.keys():
                    if a not in ['DIMENSION_LIST','CLASS','external_table']:
                        if type(f[k][v].attrs[a]) is str: 
                            fd[vlist[-1]].attrs[a]=numpy.string_(f[k][v].attrs[a])
                        else:
                            fd[vlist[-1]].attrs[a]=f[k][v].attrs[a]
                
                for a in cfv.keys():
                    if a not in ['shortname','odbcode','cdmcode']:
                        fd[vlist[-1]].attrs[a]=numpy.string_(cfv[a]) 
                    if a=='units' and cfv[a]=='NA':
                        fd[vlist[-1]].attrs[a]=numpy.string_('')
                    if a=='units' and vlist[-1]=='time':
                        ahilf=numpy.bytes_(f[k][v].attrs[a])
                        fd[vlist[-1]].attrs[a]=ahilf
                        #if b'seconds' not in ahilf:
                            #aa=ahilf.split()
                            #fd[vlist[-1]].attrs[a]=b'seconds since '+aa[2]+b' '+aa[3].split(b':')[0]+b':00:00'
                        #else:
                            #fd[vlist[-1]].attrs[a]=ahilf
                    
                
            
                
                #print(k,time.time()-tt)       
                l=0
                for d in f[k][v].dims:
                    if len(d)>0:
                        #print(k,v,f[k][v].dims[l][0].name)
                        if l==0:
                            
                            fd[vlist[-1]].dims[l].attach_scale(fd[dim0])
                        else:
                            #fd[vlist[-1]].dims[l].attach_scale(fd[f[k][v].dims[l][0].name.split('/')[-1]])
                            fd[vlist[-1]].dims[l].attach_scale(fd[sname])
                    l+=1
    
    tt=time.time()-tt                
    if tt >0.4:
        print('slow:',k,tt)
    
    

        
''' Main routine for parsing the CDS request, writing into flat netCDF file
c.retrieve(
    'reanalysis-era5-pressure-levels',
    {
        'product_type':'reanalysis',
        'format':'netcdf',
        'variable':[
            'geopotential','relative_humidity','specific_humidity',
            'specific_snow_water_content','temperature','u_component_of_wind',
            'v_component_of_wind'
        ],
        'pressure_level':[
            '450','650'
        ],
        'year':'1999',
        'month':'03',
        'day':'09',
        'time':'10:00'
    },
    'download.nc')
'''
@njit(cache=True)
def calc_trajindex(hh,hilf):
    
    #print(type(hh),type(hilf))
    j=0#hh[0]-hh[0]
    hilf[j]=hh[j]
    for i in range(hh.shape[0]-1):
        x=j
        #print(i,hh[i+1],hh[i])
        if hh[i+1]>hh[i]:
            j+=1
            hilf[j]=hh[i+1]
        hh[i]=x
    hh[-1]=j
    
    #zid=numpy.zeros(len(zidx))
    #for i in range(len(zidx)):
        #zid[i]=zidx[i]
    
    return hilf[:j+1]
                    
@njit(cache=True)
def calc_trajindexfast(z,zidx,idx,trajectory_index):
    
    #zidx=numpy.zeros(z.shape[0],dtype=numpy.int32)
    z0=zidx[0]
    j=0
    l=0
    i=0
    for i in range(z.shape[0]-1):
        jold=j
        while idx[j]>=z[i] and idx[j]<z[i+1]:
            trajectory_index[j]=l
            j+=1
            if j==idx.shape[0]:
                break
        if j>jold:
            zidx[l]=z0+i
            l+=1
        if j==idx.shape[0]:
            break
            
    if j<idx.shape[0]:
        
        if z.shape[0]>1:
            i+=1
        jold=j
        while idx[j]>=z[i]:
            trajectory_index[j]=l
            j+=1
            if j==idx.shape[0]:
                break
        if j>jold:
            zidx[l]=z0+i
            l+=1
    zidx=zidx[:l]
    
    return zidx
            
def process_flat(randdir,cf,rvars):

    t=time.time()
    vdict={}
    cdmdict={}
    cdmnamedict={}
    for k,v in cf.items():
        if "odbcode" in v.keys():
            vdict[v['cdsname']]=v['odbcode']
            cdmdict[v['cdsname']]=v['cdmcode']
            cdmnamedict[v['cdsname']]=k
        
    error=''
    rfile=''
    rvkeys=rvars.keys()
    statid=rvars['statid']
    #cost=calculate_cost(rvars) # estimate size of output file
    if 'statid' in rvkeys:
        #rfile='/fio/srvx7/leo/python/CEUAS/CEUAS/public/harvest/data/tables/'+'chera5.conv._'+rvars['statid']+'.nc'
        rfile=os.path.expandvars('$EUA_ROOT/subdaily/v0.1/source/ERA5_1/obs/0-20000-0-'+rvars['statid']+
                                 '/eua_subdaily_v0.1_source_ERA5_1_obs_0-20000-0-'+rvars['statid']+'_t.nc')
        rfile=os.path.expandvars('$RSCRATCH/era5/odbs/1/'+'chera5.conv._'+rvars['statid']+'.nc')
        print(rfile)
        if len(rvkeys)>0:
            rvdict=copy.copy(rvars)
            del rvdict['statid']
            rvdk=rvdict.keys()
            #for k in rvdk:
                #try:
                    #rvdict[k]=eval(rvdict[k])
                    #print(k,type(rvdict[k]),rvdict[k])
                #except:
                    #print('could not evaluate '+rvdict[k])
            for k in rvdk:
                if type(rvdict[k]) is not list:
                    if rvdict[k] not in cdmdict.keys():
                        rvdict[k]=[rvdict[k],rvdict[k]]
                    else:
                        rvdict[k]=[cdmdict[rvdict[k]],cdmdict[rvdict[k]]]
                                    
            #print(rfile)
            try:
                
                with h5py.File(rfile,'r') as f:
                    t=time.time()
                    if False:
                        
                        di=f['dateindex']
        
                        try:
                            sdate=rvdict.pop('date')
                            didx=numpy.where(numpy.logical_and(di[0,:]>=sdate[0],di[0,:]<=sdate[1]))[0]
                            if didx.shape[0]==0:
                                return '','Error while checking criteria'
                            didx=[didx[0],didx[-1]]
                        except KeyError:
                            didx=[0,di.shape[1]-1]
                        
                        trange=[di[1,didx[0]],di[2,didx[1]]+1]    
                    else:
                        di=numpy.array((f['recordtimestamp'][:],f['recordindex'][:]))
                        try:
                            sdate=numpy.array(rvdict.pop('date'))
                            dsec=[]
                            for d in sdate:
                                dsec.append(((datetime.datetime(year=d//10000,month=d%10000//100,day=d%100)-datetime.datetime(year=1900,month=1,day=1))).days*86400)
                            if len(dsec)==1:
                                dsec=numpy.concatenate((dsec,dsec))
                            dsec[-1]+=86399
                            didx=numpy.where(numpy.logical_and(di[0,:]>=dsec[0],di[0,:]<=dsec[-1]))[0]
                            if didx.shape[0]==0:
                                return '','Error while checking criteria'
                            didx=[didx[0],didx[-1]]
                        except KeyError:
                            didx=[0,di.shape[1]-1]
                        if didx[-1]+1==di.shape[1]:
                            trange=[di[1,didx[0]],f['observations_table']['observation_value'].shape[0]]
                        else:
                            trange=[di[1,didx[0]],di[1,didx[1]+1]]    
                    print('didx:',time.time()-t)
                    mask=numpy.ones(trange[1]-trange[0],dtype=numpy.bool)
                    criteria={'variable':'era5fb/varno@body','level':'era5fb/vertco_reference_1@body'}
                    criteria={'variable':'observations_table/observed_variable',
                              'pressure_level':'observations_table/z_coordinate',
                              'time':'observations_table/date_time'}
                    ck=criteria.keys()
                    t=time.time()
                    for ck,cv in criteria.items():
                        if ck in rvdk:
                            #mask=numpy.logical_and(mask,f[cv][trange[0]:trange[1]]>=rvdict[ck][0])
                            #mask=numpy.logical_and(mask,f[cv][trange[0]:trange[1]]<=rvdict[ck][1])
                            try:
                                if ck=='time':
                                    us=f[cv].attrs['units']
                                    dh=us.split(' ')[-1].split(':')
                                    if 'hours since' in f[cv].attrs['units']:
                                        hilf=(f[cv][trange[0]:trange[1]]+int(dh[0]))%24
                                        mask=numpy.logical_and(mask,numpy.isin(hilf,rvdict[ck]))
                                    elif 'minutes since'  in f[cv].attrs['units']:
                                        hilf=f[cv][trange[0]:trange[1]]
                                        hilf=(hilf//60)+int(dh[0])%24
                                        mask=numpy.logical_and(mask,numpy.isin(hilf,rvdict[ck]))
                                    if 'seconds' in f[k][v].attrs['units']:
                                        hilf=f[cv][trange[0]:trange[1]]
                                        hilf=(hilf//60//60+int(dh[0]))%24
                                else:
                                    mask=numpy.logical_and(mask,numpy.isin(f[cv][trange[0]:trange[1]],rvdict[ck]))
                                #print('sum',cv,rvdict[ck],numpy.sum(mask))
                            except:
                                return '','Error while checking criteria'
                            
                    print('mask:',time.time()-t)
                    #t=time.time()
                    #y=f['era5fb/date@hdr'][mask]
                    idx=numpy.where(mask)[0]+trange[0]
                    try:
                        if len(idx)>0:            
                            print(len(idx),'values found')
                        else:
                            return '','No data found'
                           
                    except:
                        print ('no data')
                        return '','No data found'
                    ## slow version ##
                    #h=numpy.asarray(f['observations_table/report_id'][trange[0]:trange[1],:].view('S5'),dtype=numpy.int32)[:,0]
                    #trajectory_index=h[idx-trange[0]]
                    ## medium version ##
                    #h=numpy.asarray(f['observations_table/report_id'][trange[0]:trange[1],:].view('S5')[idx-trange[0]],dtype=numpy.int32)[:,0]
                    #trajectory_index=h[:]
                    ##hilf=numpy.zeros(idx.shape[0],dtype=numpy.int32)
                    #print('recordindex:',time.time()-t)
                    #zidx=calc_trajindex(trajectory_index,h)
    
                    #trajectory_index_orig=trajectory_index[:]
                    ##zidx=numpy.unique(h[idx-trange[0]])
                    #z=f['recordindex'][:]
                    z=di[1,:]
                    trajectory_index=numpy.zeros_like(idx,dtype=numpy.int32)
                    zidx=numpy.where(numpy.logical_and(z>=trange[0],z<trange[1]))[0]
                    z=z[zidx]
                    
                    zidx=calc_trajindexfast(z,zidx,idx,trajectory_index)
                   
                    dims={'obs':numpy.zeros(idx.shape[0],dtype=numpy.int32),
                          'trajectory':numpy.zeros(zidx.shape[0],dtype=numpy.int32)}
                    globatts={'Conventions':"CF-1.7" ,
                              'source':"radiosonde",
                              'featureType' : "trajectory"}
                    
                    snames=['report_id','platform_id','platform_name','observation_value','latitude','longitude','time','air_pressure','trajectory_label']
                    print('rvdk',rvdk)
                    if 'variable' in rvdk:
                        if type(rvars['variable']) is list:
                            #for r in rvars['variable']:
                                #snames.append(cdmnamedict[r])
                            snames.append(cdmnamedict[rvars['variable'][0]])
                        else:
                            snames.append(cdmnamedict[rvars['variable']])
                    else:
                        return '','No variable specified'
                        
                    if 'fbstats' in rvdk:  
                        if type(rvars['fbstats']) is list:
                            for c in rvars['fbstats']:
                                snames.append(cdmnamedict[c])
                        else:
                            snames.append(cdmnamedict[rvars['fbstats']])
                    #else:
                        #rvars['fbstats']=['observation_value']
                        
    
                    ccf={}
                    for s in snames:
                        try:       
                            ccf[s]=cf[s]
                        except:
                            pass
                    
                    
                    print('recordindex:',time.time()-t)
                    dfile=randdir+'/dest_'+statid+'_'+cdmnamedict[rvars['variable']]+'.nc'
                    print(os.getcwd()+'/'+dfile)
                    with h5py.File(dfile,'w') as fd:
                        i=0
                        for d,v in dims.items():
                            fd.create_dataset(d,data=v)
                            fd[d].attrs['NAME']=numpy.string_('This is a netCDF dimension but not a netCDF variable.')
                            #fd[d].attrs['_Netcdf4Dimid']=numpy.int64(i)
                            i+=1
                        fd.create_dataset('trajectory_index',data=trajectory_index)
                        fd['trajectory_index'].attrs['long_name'] = numpy.string_("index of trajectory this obs belongs to")
                        fd['trajectory_index'].attrs['instance_dimension'] = numpy.string_("trajectory") 
                        fd['trajectory_index'].attrs['coordinates'] = numpy.string_("lat lon time plev") 
                        print('trajectory stuff:',time.time()-t)
                        for k in f.keys():
                            if isinstance(f[k],h5py.Group):
                                #t=time.time()
                                if k in ['observations_table']: # only obs, feedback fitting criteria (idx) is copied
                                    do_cfcopy(fd,f,k,idx,ccf,'obs',var_selection=['observation_id','latitude','longitude','z_coordinate',
                                                                                 'observation_value','date_time'])#'observed_variable','units'
                                    print(k,'copied',time.time()-t)
                                elif k in ['era5fb']: # only obs, feedback fitting criteria (idx) is copied
                                    do_cfcopy(fd,f,k,idx,ccf,'obs',var_selection=['fg_depar@body','an_depar@body','biascorr@body'])#['vertco_reference_1@body','obsvalue@body','fg_depar@body'])
                                    print(k,'copied',time.time()-t)
                                elif k in ['header_table']:  # only records fitting criteria (zidx) are copied
                                    do_cfcopy(fd,f,k,zidx,ccf,'trajectory',var_selection=['report_id']) #,'station_name','primary_station_id'])
                                    print(k,'copied',time.time()-t)
                                elif k in ['station_configuration']:  # only records fitting criteria (zidx) are copied
                                    #sh=f['header_table']['primary_station_id'][0].shape[0]
                                    #sid=f['header_table']['primary_station_id'][0].view('S{}'.format(sh))[0].split(b"'")[1]
                                    #sid=sid.decode('latin1')
                                    #for k in f['station_configuration']['primary_id'].values:
                                        #if sid in k[-5:]:
                                    try:
                                        
                                        sh=f[k]['primary_id'].shape[1]                                
                                        fd.attrs['primary_id']=f[k]['primary_id'][0].view('S{}'.format(sh))[0]
                                        sh=f[k]['station_name'].shape[1]                                
                                        fd.attrs['station_name']=f[k]['station_name'][0].view('S{}'.format(sh))[0]
                                        print(k,'copied',time.time()-t)
                                    except:
                                        
                                        print('no primary_id:',dfile)
                                    #fd['primary_id'][:]=f[k]['primary_id'][:]
                                    #fd['station_name'][:]=f[k]['station_name'][:]
                                    
                                else: # groups that are simply copied
                                    #print(k)
                                    if False: ## for now commented out
                                        fd.create_group(k)
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
                                                    if type(f[k][v].attrs[a]) is str:
                                                        fd[k][v].attrs[a]=numpy.string_(f[k][v].attrs[a])
                                                    else:
                                                        fd[k][v].attrs[a]=f[k][v].attrs[a]
                                            l=0
                                            for d in f[k][v].dims:
                                                if len(d)>0:
                                                    fd[k][v].dims[l].attach_scale(fd[k][f[k][v].dims[l][0].name])
                                                l+=1
                                        print(k,'copied',time.time()-t)
                        
                        fd['trajectory_label'].attrs['cf_role']=numpy.string_('trajectory_id')
                        fd['trajectory_label'].attrs['long_name']=numpy.string_('Label of trajectory')
                        
                        for a,v in globatts.items():
                            fd.attrs[a]=numpy.string_(v)
                        fd.attrs['history']=numpy.string_('Created by Copernicus Early Upper Air Service Version 0, '+ datetime.datetime.now().strftime("%d-%b-%Y %H:%M:%S"))
    
                        print(os.getcwd()+'/'+dfile)
            except MemoryError:
                return '','error while reading '+rfile

            print(time.time()-t)

        return dfile,''
    else:
        return rfile,'No station ID (statid) specified'
    
    
if __name__ == '__main__':


    read_standardnames()

    os.chdir(os.path.expandvars('$HOME/python/web2py'))
    print(('called with ',sys.argv[1],sys.argv[2]))
    rvars=eval(sys.argv[2])
    try:
        
        if '[0' in rvars['statid']:
            rvars['statid']="['"+rvars['statid'][1:7]+"']"
    except:
        
        pass
    
    #df=cdmexplainer(rvars)   
    #print(df)

    os.chdir('/raid60/scratch/leo/scratch/era5/odbs/1')
    body=rvars
    bodies=[]
    if type(body['statid']) is list:
        for s in  body['statid']:
            bodies.append(dict(body))
            bodies[-1]['statid']=s
    else:
        if body['statid'] == 'all':
            slist=glob.glob('/raid60/scratch/leo/scratch/era5/odbs/1/chera5.conv._?????.nc')
            statids=[]
            for s in slist:
                bodies.append(dict(body))
                bodies[-1]['statid']=s[-8:-3]
        else:
            bodies.append(dict(body))
    p=Pool(20)
    cf=read_standardnames()
    randdir='{:012d}'.format(100000000000)
    try:
        os.mkdir(randdir)
    except:
        pass
    func=partial(process_flat,randdir,cf)
    t=time.time()
    results=list(map(func,bodies))
    #t=time.time()
    #rfile,error=process_flat(rvars)
    print(results)
    print(time.time()-t)
    print(results[0])
    
    
