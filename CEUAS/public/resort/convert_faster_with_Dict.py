#!/usr/bin/env
# coding: utf-8

import numpy
import numpy as np
import pandas
import pandas as pd
import numba
from numba import njit
import sys,glob
import zipfile, os, time
import urllib3
from datetime import datetime, timedelta
import glob
import h5py
sys.path.append(os.getcwd()+'/../cds-backend/code/')
sys.path.append(os.getcwd()+'/../harvest/code/')
from harvest_convert_to_netCDF_newfixes import write_dict_h5
import cds_eua3 as eua
eua.logging_set_level(30)
import xarray as xr
import cdsapi, zipfile, os, time
import copy
from shutil import copyfile
import multiprocessing
sys.path.insert(0,os.getcwd()+'/../resort/rasotools-master/')
import rasotools
import warnings
warnings.filterwarnings('ignore')

opath='/raid60/scratch/uli/converted_test/'
# opath='/raid60/scratch/leo/scratch/converted_v2/'
# if there are nan values in the pressure level - we will just sort without any converting!
def do_resort(fn):
    targetfile = opath+fn.split('/')[-1]  
    
    with h5py.File(fn, 'r') as file:
        with h5py.File(targetfile, 'w') as newfile:
            groups = []
            for i in file.keys():
                if type(file[i]) == h5py._hl.group.Group:
                    newfile.create_group(i)
                    groups.append(i)
                elif i == 'recordindex' or i == 'recordtimestamp':
                    pass
                else:
                    newfile.create_dataset(i, data=file[i][:])
            for i in groups:
                if(i == 'recordindices' or i == 'observations_table' or i == 'era5fb'):
                    pass
                else:
                    for j in file[i].keys():
                        newfile[i].create_dataset(j, data=file[i][j][:])
    
    data =  eua.CDMDataset(fn)
    allvars = data.observations_table.observed_variable[()]
    allvars.sort()
    allvars = numpy.unique(allvars)
    #
    ri = data.recordindex[()]
#     print('recordindex: ', len(ri))
    rt = data.recordtimestamp[()]
    keys = data.observations_table.keys()[:-1]
    fbkeys = data.era5fb.keys()[:-1]
    # dropping all keys, where dimensions won't work - just help variabels for dimensions
    pops = []
    for i in range(len(keys)):
        if 'string' in keys[i]:
            pops.append(keys[i])
    for i in pops: keys.remove(i)
    pops = []
    for i in range(len(fbkeys)):
        if 'string' in fbkeys[i]:
            pops.append(fbkeys[i])
    for i in pops: fbkeys.remove(i)

    recordindices = [[] for i in range(len(allvars))]
    recordtimestamps = [[] for i in range(len(allvars))]

    # output variables (from observations_table)
    ov = []
    for o in keys:
        ov.append([[] for i in range(len(allvars))])
    fb = []
    for o in fbkeys:
        fb.append([[] for i in range(len(allvars))])
    #
    # loading the observed_variables
    #
    obsv = data.observations_table.observed_variable[:]
    #
    # resorting the data
    #
#     print('resort:start')
    @njit
    def make_vrindex(vridx,ridx,idx):
        l=0
        for i in range(1,len(idx)): # to set the recordindices
            if ridx[i]>ridx[i-1]:
                vridx[ridx[i-1]:ridx[i]]=l # next record after l
                l=i
        vridx[ridx[i]:]=len(idx) # next record for the last element is the len of the data


    tt=time.time()

    ridxall=np.zeros(obsv.shape[0],dtype=np.int64) # reverse index - index of the record index
    j=-1
    for j in range(len(ri)-1):
        ridxall[ri[j]:ri[j+1]]=j
    j+=1
    ridxall[ri[j]:]=j # for the last elemenet
    ridx=[]
    vridx=[]
    absidx=[]
    abscount=0
    for j in range(len(allvars)):
        idx=np.where(obsv==allvars[j])[0] # index of all elements form certain variable j
#         print(j,len(idx),',',end='')
        vridx.append(np.zeros(ri.shape[0],dtype=np.int64)) # all zeros in lenght of record index
        ridx=ridxall[idx] # ridxall where variable is j
        make_vrindex(vridx[-1],ridx,idx)
        vridx[-1]+=abscount # abscount for stacking the recordindex

        absidx.append(copy.copy(idx)) # why copy? - to make sure it's not just the ref. - maybe ok without the cp
        abscount+=len(idx)

#     print('')
    #
    # finishing the sorting 
    #
    absidx=np.concatenate(absidx)
    #
    # recordtimestamps are only necessary once
    #
    recordtimestamps = recordtimestamps[0]
    #
    # targetfile has to be a copy of the original file
    #
#     targetfile = '/raid60/scratch/uli/converted_v2/'+fn.split('/')[-1]# 0-20000-0-63894_CEUAS_merged_v0.nc'
    if os.path.isfile(targetfile):
        mode='r+'
    else:
        mode='w'
#     print()
#     print('writing '+targetfile)

    for i in range(len(keys)):
        ov_vars = data.observations_table[keys[i]][:]
        ov_vars = ov_vars[absidx]
        if keys[i] == 'index':
            pass
        elif keys[i] == 'observation_id' or keys[i] == 'report_id' or keys[i] == 'sensor_id' or keys[i] == 'source_id':
            alldict = {keys[i]:np.asarray(ov_vars, dtype='S1')}
            write_dict_h5(targetfile, alldict, 'observations_table', {keys[i]: { 'compression': 'gzip' } }, [keys[i]])
        else:
            alldict = pandas.DataFrame({keys[i]:ov_vars})
            write_dict_h5(targetfile, alldict, 'observations_table', {keys[i]: { 'compression': 'gzip' } }, [keys[i]])  

    for i in range(len(fbkeys)):
        fb_vars = data.era5fb[fbkeys[i]][:]
        fb_vars = fb_vars[absidx]
        if fbkeys[i] == 'index' or fbkeys[i] == 'string6' or fbkeys[i] == 'string7' or fbkeys[i] == 'string10':
            pass
        elif fbkeys[i] == 'expver' or fbkeys[i] == 'source@hdr' or fbkeys[i] == 'source_id' or fbkeys[i] == 'statid@hdr':
            alldict = {fbkeys[i]:np.asarray(fb_vars, dtype='S1')}
            write_dict_h5(targetfile, alldict, 'era5fb', {fbkeys[i]: { 'compression': 'gzip' } }, [fbkeys[i]])
        else:
            alldict = pandas.DataFrame({fbkeys[i]:fb_vars})
            write_dict_h5(targetfile, alldict, 'era5fb', {fbkeys[i]: { 'compression': 'gzip' } }, [fbkeys[i]]) 
    #
    # writing the recordindices and recordtimestamp.
    #       
    recordindices=vridx
    for i in range(len(recordindices)):
        testvar = pandas.DataFrame({str(allvars[i]):recordindices[i]})
        write_dict_h5(targetfile, testvar, 'recordindices', {str(allvars[i]): { 'compression': None } }, [str(allvars[i])]) 

    write_dict_h5(targetfile, {'recordtimestamp':rt}, 'recordindices', {'recordtimestamp': { 'compression': None } }, ['recordtimestamp'])

    print('elapsed:',time.time()-tt)

from numba.typed import List

@njit(boundscheck=True)
def ipl(observed_variable,observation_value,z_coordinate,z_coordinate_type,recordtimestamp):
    jdx=numpy.where(numpy.logical_and(observed_variable==85,~numpy.isnan(observation_value)))[0]
    idx=numpy.empty(len(jdx),dtype=numpy.int64)
    press=numpy.empty(len(jdx),dtype=z_coordinate.dtype)
    relhum=numpy.empty(len(jdx),dtype=observation_value.dtype)
    dewpoint=numpy.empty(len(jdx),dtype=observation_value.dtype)
    dpd=numpy.empty(len(jdx),dtype=observation_value.dtype)
    spechum=numpy.empty(len(jdx),dtype=observation_value.dtype)
    
    uwind=numpy.empty(len(jdx),dtype=observation_value.dtype)
    vwind=numpy.empty(len(jdx),dtype=observation_value.dtype)
    ws=numpy.empty(len(jdx),dtype=observation_value.dtype)
    wd=numpy.empty(len(jdx),dtype=observation_value.dtype)

    for v in relhum,dewpoint,dpd,spechum:
        v.fill(numpy.nan)
    temp=numpy.empty(len(jdx),dtype=observation_value.dtype)
    p=z_coordinate[0]-1.
    rts=recordtimestamp[0]-1
    j=0
    good=True
    for i in range(observed_variable.shape[0]):
        if z_coordinate[i]==z_coordinate[i] and (z_coordinate[i]!=p or recordtimestamp[i]!=rts):
            if not good:
                j-=1
            if j<temp.shape[0]:
                idx[j]=i
                press[j]=z_coordinate[i]
            good=False
            p=z_coordinate[i]
            rts=recordtimestamp[i]
            j+=1
            if j==2446:
                print(j)

        if observed_variable[i]==38 and observation_value[i]==observation_value[i]:
            relhum[j-1]=observation_value[i]
        if observed_variable[i]==39 and observation_value[i]==observation_value[i]:
            if j<=spechum.shape[0]:
                spechum[j-1]=observation_value[i]
        if observed_variable[i]==34 and observation_value[i]==observation_value[i]:
            dpd[j-1]=observation_value[i]
        if observed_variable[i]==36 and observation_value[i]==observation_value[i] and z_coordinate[i]==z_coordinate[i]:
            dewpoint[j-1]=observation_value[i]
        if observed_variable[i]==85 and observation_value[i]==observation_value[i]:
            temp[j-1]=observation_value[i]
            good=True

    return idx,press,temp,relhum,spechum,dpd,dewpoint

@njit(boundscheck=True)
def ipl2(observed_variable,observation_value,z_coordinate,z_coordinate_type,recordtimestamp):
    
    jdx=0
    dpress=-1.
    dts=-1
    for i in range(observation_value.shape[0]):
        if z_coordinate[i]!=dpress or recordtimestamp[i]!=dts:
            dpress=z_coordinate[i]
            dts=recordtimestamp[i]
            jdx+=1
                
                
    #jdx=len(numpy.where(numpy.logical_and(observed_variable==85,~numpy.isnan(observation_value)))[0]
    idx=numpy.empty(jdx,dtype=numpy.int64)
    press=numpy.empty(jdx,dtype=z_coordinate.dtype) # either height or pressure, depending on coordinate type
    relhum=numpy.empty(jdx,dtype=observation_value.dtype)
    dewpoint=numpy.empty(jdx,dtype=observation_value.dtype)
    dpd=numpy.empty(jdx,dtype=observation_value.dtype)
    spechum=numpy.empty(jdx,dtype=observation_value.dtype)
    
    uwind=numpy.empty(jdx,dtype=observation_value.dtype)
    vwind=numpy.empty(jdx,dtype=observation_value.dtype)
    ws=numpy.empty(jdx,dtype=observation_value.dtype)
    wd=numpy.empty(jdx,dtype=observation_value.dtype)

    temp=numpy.empty(jdx,dtype=observation_value.dtype)
    
    for v in temp,relhum,dewpoint,dpd,spechum,uwind,vwind,wd,ws:
        v.fill(numpy.nan)
    p=z_coordinate[0]-1.
    rts=recordtimestamp[0]-1
    j=-1
    good=True
    for i in range(observation_value.shape[0]):
        if z_coordinate[i]!=p or recordtimestamp[i]!=rts:
            p=z_coordinate[i]
            rts=recordtimestamp[i]
            j+=1
            press[j]=p
            idx[j]=i
    
        if observed_variable[i]==38:
            relhum[j]=observation_value[i]
        elif observed_variable[i]==39:
            spechum[j]=observation_value[i]
        elif observed_variable[i]==34:
            dpd[j]=observation_value[i]
        elif observed_variable[i]==36:
            dewpoint[j]=observation_value[i]
        elif observed_variable[i]==85:
            temp[j]=observation_value[i]
        elif observed_variable[i]==104:
            uwind[j]=observation_value[i]
        elif observed_variable[i]==105:
            vwind[j]=observation_value[i]
        elif observed_variable[i]==106:
            wd[j]=observation_value[i]
        elif observed_variable[i]==107:
            ws[j]=observation_value[i]
        else:
            pass
                
            
    print(j,jdx)
    return idx,press,temp,relhum,spechum,dpd,dewpoint,uwind,vwind,wd,ws

@njit(boundscheck=True)
def qconvert(j,k,h,a_observation_value,a_conversion_flag,a_conversion_method,temp,cdpddp,cdpdrh,crhdpd,cshrh,cshdpd):
    if h==34:
        if cdpddp[k]==cdpddp[k]:
            a_observation_value[j]=cdpddp[k]
            if numpy.abs(cdpddp[k])>50:
                #print(k,cdpddp[k],cdpdrh[k],temp[k],press[k],dewpoint[k],i-1)
                a_observation_value[j]=numpy.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=2
        else:
            a_observation_value[j]=cdpdrh[k]
            a_conversion_flag[j]=0
            a_conversion_method[j]=3
            
    elif h==36:
        a_observation_value[j]=temp[k]-cdpdrh[k]
        a_conversion_flag[j]=0
        a_conversion_method[j]=3
    elif h==38:
        a_observation_value[j]=crhdpd[k]
        a_conversion_flag[j]=0
        a_conversion_method[j]=2
    elif h==39:
        if cshrh[k]==cshrh[k]:
            a_observation_value[j]=cshrh[k]
            a_conversion_flag[j]=0
            a_conversion_method[j]=3
        else:
            a_observation_value[j]=cshdpd[k]
            a_conversion_flag[j]=0
            a_conversion_method[j]=2
    else:
        print('qconvert called with wrong variable')
        
    return

@njit(boundscheck=True)
def wconvert(j,k,h,a_observation_value,a_conversion_flag,a_conversion_method,cuwind,cvwind,cwd,cws):
    if h==104:
        if cuwind[k]==cuwind[k]:
            a_observation_value[j]=cuwind[k]
            a_conversion_flag[j]=0
            a_conversion_method[j]=1
    elif h==105:
        if cvwind[k]==cvwind[k]:
            a_observation_value[j]=cvwind[k]
            a_conversion_flag[j]=0
            a_conversion_method[j]=1
    elif h==106:
        if cwd[k]==cwd[k]:
            a_observation_value[j]=cwd[k]
            a_conversion_flag[j]=0
            a_conversion_method[j]=2
    elif h==107:
        if cws[k]==cws[k]:
            a_observation_value[j]=cws[k]
            a_conversion_flag[j]=0
            a_conversion_method[j]=2
    else:
        print('wconvert called with wrong variable')
    

@njit(boundscheck=True)          
def augment(obstab_f, f_keys, obstab_i32, i32_keys, obstab_i64, i64_keys, observation_id, report_id, sensor_id, source_id,
            # observed_variable,observation_value,z_coordinate,z_coordinate_type,date_time,conversion_flag,conversion_method,
            idx,temp,press,relhum,spechum,dpd,dewpoint,uwind,vwind,wd,ws,
            cdpddp,cdpdrh,cshrh,cshdpd,crhdpd,cuwind,cvwind,cwd,cws,humvar,wvar):
    
    a_obstab_i32 = numba.typed.Dict()
    for o in i32_keys:
#         a_obstab_i32[o] = numpy.empty((obstab_i32[o].shape[0]*3,obstab_i32[o].shape[1]),obstab_i32[o].dtype)
        a_obstab_i32[o] = numpy.empty(obstab_i32[o].shape[0]*3,obstab_i32[o].dtype)
    a_obstab_i64 = numba.typed.Dict()
    for o in i64_keys:
        a_obstab_i64[o] = numpy.empty(obstab_i64[o].shape[0]*3,obstab_i64[o].dtype)
    a_obstab_f = numba.typed.Dict()
    for o in f_keys:
        a_obstab_f[o] = numpy.empty(obstab_f[o].shape[0]*3,obstab_f[o].dtype)

    a_observation_id=numpy.empty((observation_id.shape[0]*3,observation_id.shape[1]),observation_id.dtype)
    a_report_id=numpy.empty((report_id.shape[0]*3,report_id.shape[1]),report_id.dtype)
    a_sensor_id=numpy.empty((sensor_id.shape[0]*3,sensor_id.shape[1]),sensor_id.dtype)
    a_source_id=numpy.empty((source_id.shape[0]*3,source_id.shape[1]),source_id.dtype)

    recordindex=numpy.empty(idx.shape[0],obstab_i64['date_time'].dtype)
    recordtimestamp=numpy.empty(idx.shape[0],obstab_i64['date_time'].dtype)
    
    j=-1 # augmented index
    p= obstab_f['z_coordinate'][0]-1 # z_coordinate[0]-1
    rts= obstab_i64['date_time'][0] # date_time[0]
#     from numpy import recarray
    humlist=List()
    wlist=List()
    recordindex[0]=0
    recordtimestamp[0]=obstab_i64['date_time'][0] # date_time[0]
    ri=1
    for k in range(idx.shape[0]):  # no new levels are created, but variables per level will increase
        if k==idx.shape[0]-2:
            print('vor Schluss')
        if k==idx.shape[0]-1:
            idxu= obstab_f['observation_value'].shape[0] # observation_value.shape[0]
        else:
            idxu=idx[k+1]
        for i in range(idx[k],idxu):
            j+=1
            if obstab_f['observation_value'][i]==obstab_f['observation_value'][i]: # observation_value[i]==observation_value[i]:
                for o in f_keys:
                    a_obstab_f[o][j]=obstab_f[o][i]
                for o in i32_keys:
                    a_obstab_i32[o][j]=obstab_i32[o][i]
                for o in i64_keys:
                    a_obstab_i64[o][j]=obstab_i64[o][i]
                a_observation_id[j]=observation_id[i]
                a_report_id[j]=report_id[i]
                a_sensor_id[j]=sensor_id[i]
                a_source_id[j]=source_id[i]
                if int(obstab_i64['observed_variable'][i]) in humvar:
                    humlist.append(int(obstab_i64['observed_variable'][i]))
                elif int(obstab_i64['observed_variable'][i]) in wvar:
                    wlist.append(int(obstab_i64['observed_variable'][i]))
            else:
                # writing all obstab vars and overwriting them if needed:
                for o in f_keys:
                    a_obstab_f[o][j]=obstab_f[o][i]
                for o in i32_keys:
                    a_obstab_i32[o][j]=obstab_i32[o][i]
                for o in i64_keys:
                    a_obstab_i64[o][j]=obstab_i64[o][i]
                a_observation_id[j]=observation_id[i]
                a_report_id[j]=report_id[i]
                a_sensor_id[j]=sensor_id[i]
                a_source_id[j]=source_id[i]
                if int(obstab_i64['observed_variable'][i]) in humvar:
                    qconvert(j,k,obstab_i64['observed_variable'][i],a_obstab_f['observation_value'],a_obstab_i32['conversion_flag'],a_obstab_i32['conversion_method'],temp,cdpddp,cdpdrh,crhdpd,cshrh,cshdpd)
                    # set rest of obstable too
                elif int(obstab_i64['observed_variable'][i]) in wvar:
                    wconvert(j,k,obstab_i64['observed_variable'][i],a_obstab_f['observation_value'],a_obstab_i32['conversion_flag'],a_obstab_i32['conversion_method'],cuwind,cvwind,cwd,cws)
                    # set rest of obstable too
#                 else:
#                     a_observation_value[j]=observation_value[i]
#                     a_conversion_flag[j]=conversion_flag[i]
#                     a_conversion_method[j]=a_conversion_method[i]
        if humlist:
            for h in humvar:
                if h not in humlist:
                    j+=1
                    # writing all obstab vars and overwriting them if needed:
                    for o in f_keys:
                        a_obstab_f[o][j]=obstab_f[o][i]
                    for o in i32_keys:
                        a_obstab_i32[o][j]=obstab_i32[o][i]
                    for o in i64_keys:
                        a_obstab_i64[o][j]=obstab_i64[o][i]
                    a_observation_id[j]=observation_id[i]
                    a_report_id[j]=report_id[i]
                    a_sensor_id[j]=sensor_id[i]
                    a_source_id[j]=source_id[i]
                    
                    a_obstab_i64['observed_variable'][j]=h
                    qconvert(j,k,h,a_obstab_f['observation_value'],a_obstab_i32['conversion_flag'],a_obstab_i32['conversion_method'],temp,cdpddp,cdpdrh,crhdpd,cshrh,cshdpd)

                    if a_obstab_f['observation_value'][j]!=a_obstab_f['observation_value'][j]:
                        j-=1
            humlist.clear()
        if wlist:
            for h in wvar:
                if h not in wlist:
                    j+=1
                    # writing all obstab vars and overwriting them if needed:
                    for o in f_keys:
                        a_obstab_f[o][j]=obstab_f[o][i]
                    for o in i32_keys:
                        a_obstab_i32[o][j]=obstab_i32[o][i]
                    for o in i64_keys:
                        a_obstab_i64[o][j]=obstab_i64[o][i]
                    a_observation_id[j]=observation_id[i]
                    a_report_id[j]=report_id[i]
                    a_sensor_id[j]=sensor_id[i]
                    a_source_id[j]=source_id[i]
                    
                    a_obstab_i64['observed_variable'][j]=h
                    wconvert(j,k,h,a_obstab_f['observation_value'],a_obstab_i32['conversion_flag'],a_obstab_i32['conversion_method'],cuwind,cvwind,cwd,cws)

                    if a_obstab_f['observation_value'][j]!=a_obstab_f['observation_value'][j]:
                        j-=1
            wlist.clear()
        if idxu!=obstab_f['observation_value'].shape[0]:        
            if obstab_i64['date_time'][idxu] != obstab_i64['date_time'][idx[k]]:
                recordindex[ri]=j+1
                recordtimestamp[ri]=obstab_i64['date_time'][idxu]
                ri+=1

    
        if k%100000==0:
            print(k,idx.shape[0])
#     print(k,j,i,observed_variable.shape[0],a_observed_variable.shape[0])
    j=j+1
    out = {}
    for o in f_keys:
        out[o]=obstab_f[o][:j]
    for o in i32_keys:
        out[o]=obstab_i32[o][:j]
    for o in i64_keys:
        out[o]=obstab_i64[o][:j]
    out['observation_id']=observation_id[:j]
    out['report_id']=report_id[:j]
    out['sensor_id']=sensor_id[:j]
    out['source_id']=source_id[:j]
    out['ri']=recordindex[:ri]
    out['rt']=recordtimestamp[:ri]
    return out #a_observed_variable[:j],a_observation_value[:j],a_z_coordinate[:j],a_z_coordinate_type[:j],a_date_time[:j],a_conversion_flag[:j],a_conversion_method[:j],recordindex[:ri],recordtimestamp[:ri]
    

def convert_missing(fn, destination: str = opath):
    tt=time.time()
    nanlist = [float('nan'), np.nan, 0, -2147483648]
    with eua.CDMDataset(fn) as data:
        arrayconverter = data.to_dataframe(groups='observations_table', variables=['observed_variable'])
        arrayconverter = arrayconverter.observed_variable.head(1).to_xarray()
        rto = data.recordtimestamp[:]
        rio = data.recordindex[:]
        keys = data.observations_table.keys()
#         keys = [x for x in keys if x in ['observation_id','date_time','observed_variable','z_coordinate','z_coordinate_type','observation_value','conversion_method','conversion_flag']]
        keys = [x for x in keys if not x.startswith('string')]
        keys.remove('index')
        keys.remove('shape')
        obskeys = keys
        obstab_writetofile = [[] for i in range(len(obskeys))]
        
        keys = data.era5fb.keys()
        keys = [x for x in keys if x in ['fg_depar@body','an_depar@body','biascorr@body','biascorr_fg@body']]
        #keys = [x for x in keys if not x.startswith('string')]
        #keys.remove('index')
        #keys.remove('shape')
        fg_depar = keys.index('fg_depar@body')
        depar = keys.index('an_depar@body')
        biascorr = keys.index('biascorr@body')
        fg_biascorr = keys.index('biascorr_fg@body')
        fbkeys = keys
        fb_writetofile = [[] for i in range(len(fbkeys))]
        
        recidxlen = len(data.recordindex[:])
        
        addtorecordindex = [] # will be filled with the count of how many variables have been added before !! addtorecordindex[0] has to be added to recordindex [1] !!
        addedvarscount = 0 # will grow with every variable added
        onlyone = False
        if recidxlen == 1:
            recidxlen = 2
            onlyone = True
        
        # loading data:
        loaded_data_f = numba.typed.Dict() # {}
        loaded_data_i32 = numba.typed.Dict() # {}
        loaded_data_i64 = numba.typed.Dict() # {}
        a_keys = copy.copy(obskeys)
        a_keys.remove('observation_id')
        a_keys.remove('report_id')
        a_keys.remove('sensor_id')
        a_keys.remove('source_id')
        i32_keys = []
        i64_keys = []
        f_keys = []
        for o in a_keys:
            if type(data.observations_table[o][0]) == numpy.float32 :
                f_keys.append(o)
                loaded_data_f[o] = numpy.asarray(data.observations_table[o][:])                 
            elif type(data.observations_table[o][0]) == numpy.int32:
                i32_keys.append(o)
                loaded_data_i32[o] = numpy.asarray(data.observations_table[o][:])
            elif type(data.observations_table[o][0]) == numpy.int64:
                i64_keys.append(o)
                loaded_data_i64[o] = numpy.asarray(data.observations_table[o][:])
        loaded_observation_id = numpy.asarray(data.observations_table['observation_id'][:])
        loaded_report_id= numpy.asarray(data.observations_table['report_id'][:])
        loaded_sensor_id = numpy.asarray(data.observations_table['sensor_id'][:])
        loaded_source_id = numpy.asarray(data.observations_table['source_id'][:])
        loaded_fb = {}
        
        for o in fbkeys:
            loaded_fb[o] = data.era5fb[o][:]
            
        recordindex = data.recordindex[:]
        # --->

    print(time.time()-tt)
    idx,press,temp,relhum,spechum,dpd,dewpoint,uwind,vwind,wd,ws=ipl2(loaded_data_i64['observed_variable'],loaded_data_f['observation_value'],loaded_data_f['z_coordinate'],
                                  loaded_data_i32['z_coordinate_type'],loaded_data_i64['date_time'])
#    idx=numpy.array(result[0])
    xtemp=xr.DataArray(temp)
    xpress=xr.DataArray(press)
    xrelhum=xr.DataArray(relhum)
    xspechum=xr.DataArray(spechum)
    xdpd=xr.DataArray(dpd)
    xdewpoint=xr.DataArray(dewpoint)
    cdpddp=temp-dewpoint
    cdpdrh=rasotools.met.convert.to_dpd(temp=xtemp,press=xpress,rel_humi=xrelhum).values
#    cdpdsh=rasotools.met.convert.to_dpd(temp=xtemp,press=xpress,spec_humi=xspechum).values
    cshrh = rasotools.met.convert.to_sh(temp=xtemp, press=xpress, rel_humi=xrelhum).values
    cshdpd = rasotools.met.convert.to_sh(dpd=xtemp-xdewpoint, press=xpress, temp=xtemp).values
#    crhsh = rasotools.met.convert.to_rh(temp=xtemp, spec_humi=xspechum, press=xpress).values
    crhdpd = rasotools.met.convert.to_rh(temp=xtemp,dpd=xtemp-xdewpoint).values

    idy=numpy.where(loaded_data_i32['z_coordinate_type'][idx]==2) # do not convert humidity if data are not on pressure coordinates
    for c in cdpdrh,cshrh,cshdpd,crhdpd:
        c[idy]=numpy.nan
    
    cuwind = ws * np.cos(np.radians(wd))
    cvwind = ws * np.sin(np.radians(wd))
    cws = np.sqrt(uwind ** 2 + vwind ** 2)
    cwd = 90 - np.arctan2(-vwind, -uwind) * 180 / np.pi - 180.
    cwd = np.where(cwd > 0., cwd, 360.+cwd)

    humvar=numpy.array([34,36,38,39]) #dpd,dp,rh,sh
    wvar=numpy.array([104,105,106,107]) #dpd,dp,rh,sh
#     alist=augment(loaded_data, obskeys, 
#                   # loaded_data['observed_variable'],loaded_data['observation_value'],loaded_data['z_coordinate'], 
#                   # loaded_data['z_coordinate_type'],loaded_data['date_time'],loaded_data['conversion_flag'],loaded_data['conversion_method'],
#                   idx,temp,press,relhum,spechum,dpd,dewpoint,uwind,vwind,wd,ws,
#                   cdpddp,cdpdrh,cshrh,cshdpd,crhdpd,cuwind,cvwind,cwd,cws,humvar,wvar)
    avars=augment(loaded_data_f, f_keys, loaded_data_i32, i32_keys, loaded_data_i64, i64_keys, loaded_observation_id, loaded_report_id, loaded_sensor_id, loaded_source_id,
                  # loaded_data['observed_variable'],loaded_data['observation_value'],loaded_data['z_coordinate'], 
                  # loaded_data['z_coordinate_type'],loaded_data['date_time'],loaded_data['conversion_flag'],loaded_data['conversion_method'],
                  idx,temp,press,relhum,spechum,dpd,dewpoint,uwind,vwind,wd,ws,
                  cdpddp,cdpdrh,cshrh,cshdpd,crhdpd,cuwind,cvwind,cwd,cws,humvar,wvar)
    
#     avarkeys='observed_variable','observation_value','z_coordinate','z_coordinate_type','date_time','conversion_flag','conversion_method','ri','rt'
    avarkeys = copy.copy(obskeys)
    avarkeys.extend(['ri','rt'])
#     avars=dict()
#     for i in range(len(avarkeys)):
#         avars[avarkeys[i]]=alist[i]
    
    print(time.time()-tt)
    
    import matplotlib.pylab as plt
#    plt.plot(loaded_data['date_time'][idx],crhdpd*100)
    idz=numpy.where(numpy.logical_and(cdpdrh==cdpdrh ,press==50000))[0]
    plt.plot(loaded_data_i64['date_time'][idx[idz]],cdpdrh[idz],linewidth=3)
    idy=np.where(avars['observed_variable']==34)[0]
    idzz=numpy.where(numpy.logical_and(~numpy.isnan(avars['observation_value'][idy]),avars['z_coordinate'][idy]==50000))[0]
    plt.plot(avars['date_time'][idy[idzz]],avars['observation_value'][idy[idzz]])
    #plt.show()
    
          
    # sorting:
    print('start sorting')
    targetfile = destination+fn.split('/')[-1]
    if os.path.isfile(targetfile):
        try:
            os.remove(targetfile)
        except:
            print('file could not be removed - overwriting will lead to errors')
    
    with h5py.File(fn, 'r') as file:
        with h5py.File(targetfile, 'w') as newfile:
            groups = []
            for i in file.keys():
                if type(file[i]) == h5py._hl.group.Group:
                    if i not in ('observations_table','era5fb'):
                        newfile.create_group(i)
                        groups.append(i)
                elif i == 'recordindex' or i == 'recordtimestamp':
                    pass
                else:
                    newfile.create_dataset(i, data=file[i][:])
            for i in groups:
                if(i == 'recordindices' or i == 'observations_table' or i == 'era5fb'):
                    pass
                else:
                    for j in file[i].keys():
                        newfile[i].create_dataset(j, data=file[i][j][:])
    
#     data =  eua.CDMDataset(fn)
    allvars = copy.copy(avars['observed_variable'])
    allvars.sort()
    allvars = numpy.unique(allvars)
    #
    obsv = avars['observed_variable']
    #
    # resorting the data
    #
#     print('resort:start')
    @njit
    def make_vrindex(vridx,ridx,idx):
        l=0
        for i in range(1,len(idx)): # to set the recordindices
            if ridx[i]>ridx[i-1]:
                vridx[ridx[i-1]:ridx[i]]=l # next record after l
                l=i
        vridx[ridx[i]:]=len(idx) # next record for the last element is the len of the data

    ri=avars['ri']
    ridxall=np.zeros(obsv.shape[0],dtype=np.int64) # reverse index - index of the record index
    j=-1
    for j in range(len(ri)-1):
        ridxall[ri[j]:ri[j+1]]=j
    j+=1
    ridxall[ri[j]:]=j # for the last elemenet
    ridx=[]
    vridx=[]
    absidx=[]
    abscount=0
    for j in range(len(allvars)):
        idx=np.where(obsv==allvars[j])[0] # index of all elements form certain variable j
#         print(j,len(idx),',',end='')
        vridx.append(np.zeros(ri.shape[0],dtype=np.int64)) # all zeros in lenght of record index
        ridx=ridxall[idx] # ridxall where variable is j
        make_vrindex(vridx[-1],ridx,idx)
        vridx[-1]+=abscount # abscount for stacking the recordindex

        absidx.append(copy.copy(idx)) # why copy? - to make sure it's not just the ref. - maybe ok without the cp
        abscount+=len(idx)

    # finishing the sorting 
    #
    absidx=np.concatenate(absidx)
#     try:
#         absidx=np.concatenate(absidx)
#     except:
#         absidx = absidx[0]
    #
    # recordtimestamps are only necessary once
    #
    recordtimestamps = avars['rt']
    #
    # targetfile has to be a copy of the original file
    #
    print('elapsed converting: ',time.time()-tt)
    tt=time.time()
    if os.path.isfile(targetfile):
        mode='r+'
    else:
        mode='w'
#     print()
#     print('writing '+targetfile)
    
    for i in range(len(obskeys)):
        ov_vars = avars[obskeys[i]]
        ov_vars = ov_vars[absidx]
        if obskeys[i] == 'index':
            pass
        elif obskeys[i] == 'observation_id' or obskeys[i] == 'report_id' or obskeys[i] == 'sensor_id' or obskeys[i] == 'source_id':
            alldict = {obskeys[i]:np.asarray(ov_vars, dtype='S1')}
            write_dict_h5(targetfile, alldict, 'observations_table', {obskeys[i]: { 'compression': 'gzip' } }, [obskeys[i]])
        else:
            alldict = pandas.DataFrame({obskeys[i]:ov_vars})
            write_dict_h5(targetfile, alldict, 'observations_table', {obskeys[i]: { 'compression': 'gzip' } }, [obskeys[i]])  

# geht noch nicht
    #for i in range(len(fbkeys)):
        #fb_vars = np.asarray(fb_writetofile[i]) # data.era5fb[fbkeys[i]][:]
        #fb_vars = fb_vars[absidx]
        #if fbkeys[i] == 'index' or fbkeys[i] == 'string6' or fbkeys[i] == 'string7' or fbkeys[i] == 'string10':
            #pass
        #elif fbkeys[i] == 'expver' or fbkeys[i] == 'source@hdr' or fbkeys[i] == 'source_id' or fbkeys[i] == 'statid@hdr':
            #alldict = {fbkeys[i]:np.asarray(fb_vars, dtype='S1')}
            #write_dict_h5(targetfile, alldict, 'era5fb', {fbkeys[i]: { 'compression': 'gzip' } }, [fbkeys[i]])
        #else:
            #alldict = pandas.DataFrame({fbkeys[i]:fb_vars})
            #write_dict_h5(targetfile, alldict, 'era5fb', {fbkeys[i]: { 'compression': 'gzip' } }, [fbkeys[i]]) 
    #
    # writing the recordindices and recordtimestamp.
    #       
    recordindices=vridx
    for i in range(len(recordindices)):
        testvar = pandas.DataFrame({str(allvars[i]):recordindices[i]})
        write_dict_h5(targetfile, testvar, 'recordindices', {str(allvars[i]): { 'compression': None } }, [str(allvars[i])]) 

    write_dict_h5(targetfile, {'recordtimestamp':recordtimestamps}, 'recordindices', {'recordtimestamp': { 'compression': None } }, ['recordtimestamp'])

    print('elapsed writing:',time.time()-tt)
    return
    

files = glob.glob('/raid60/scratch/federico/MERGED_DATABASE_OCTOBER2020_sensor/0-20000-0-01*.nc')
# print(files[:10])

# convert_missing(files[6020])
# convert_missing('/raid60/scratch/federico/MERGED_DATABASE_OCTOBER2020_sensor/0-20000-0-03414_CEUAS_merged_v0.nc')

if __name__ == '__main__':
#    pool = multiprocessing.Pool(processes=20)
#    result_list = pool.map(convert_missing, files[100:1000])
    idx=0
    for f in files:
        if '20000-0-01384' in f:
            print(idx)
            break
        idx+=1
    result_list = list(map(convert_missing, [files[idx]]))
    print(result_list)