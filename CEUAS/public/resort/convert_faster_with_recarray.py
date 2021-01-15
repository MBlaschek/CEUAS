#!/usr/bin/env
# coding: utf-8

import numpy
import numpy as np
import pandas
import pandas as pd
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
def ipl2(lobs):
    
    observed_variable=lobs['observed_variable']
    observation_value=lobs['observation_value']
    z_coordinate=lobs['z_coordinate']
    z_coordinate_type=lobs['z_coordinate_type']
    recordtimestamp=lobs['date_time']
    
    jdx=0
    dpress=-1.
    dts=-1
    for i in range(observation_value.shape[0]):
        if z_coordinate[i]!=dpress or recordtimestamp[i]!=dts:
            dpress=z_coordinate[i]
            dts=recordtimestamp[i]
            jdx+=1
                
                
    ##jdx=len(numpy.where(numpy.logical_and(observed_variable==85,~numpy.isnan(observation_value)))[0]
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
   #good=True
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
                
    #return        
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

@njit
def do_copy(a_obstab,obstab,j,i):
    # all written here - * will be overwritten if it's a converted variable
    a_obstab['date_time'][j]=obstab['date_time'][i]
    a_obstab['observation_id'][j]=obstab['observation_id'][i] # *
    a_obstab['observation_value'][j]=obstab['observation_value'][i] # *
    a_obstab['observed_variable'][j]=obstab['observed_variable'][i] # *
    a_obstab['z_coordinate'][j]=obstab['z_coordinate'][i]
    a_obstab['z_coordinate_type'][j]=obstab['z_coordinate_type'][i]
    a_obstab['conversion_flag'][j]=obstab['conversion_flag'][i] # *
    a_obstab['conversion_method'][j]=obstab['conversion_method'][i] # *
    return

@njit(boundscheck=True)          
def augment(obstab, a_obstab,obskeys,
             idx,temp,press,relhum,spechum,dpd,dewpoint,uwind,vwind,wd,ws,
             cdpddp,cdpdrh,cshrh,cshdpd,crhdpd,cuwind,cvwind,cwd,cws,humvar,wvar):
    
    print(obskeys,humvar)
    
    recordindex=numpy.empty(idx.shape[0],obstab['date_time'].dtype)
    recordtimestamp=numpy.empty(idx.shape[0],obstab['date_time'].dtype)
    
    j=-1 # augmented index
    p= obstab['z_coordinate'][0]-1 # z_coordinate[0]-1
    rts= obstab['date_time'][0] # date_time[0]
    
    addedvar=List()
    humlist=List()
    wlist=List()
    recordindex[0]=0
    recordtimestamp[0]=obstab['date_time'][0] # date_time[0]
    ri=1
    for k in range(idx.shape[0]):  # no new levels are created, but variables per level will increase
        if k==idx.shape[0]-2:
            print('vor Schluss')
        if k==idx.shape[0]-1:
            idxu= obstab['observation_value'].shape[0] # observation_value.shape[0]
        else:
            idxu=idx[k+1]
        for i in range(idx[k],idxu):
            j+=1
            if obstab['observation_value'][i]==obstab['observation_value'][i]: 
                do_copy(a_obstab,obstab,j,i)
                if obstab['observed_variable'][i] in humvar:
                    humlist.append(obstab['observed_variable'][i])
                elif obstab['observed_variable'][i] in wvar:
                    wlist.append(obstab['observed_variable'][i])
            else:
                do_copy(a_obstab,obstab,j,i)
                if obstab['observed_variable'][i] in humvar:
                    qconvert(j,k,obstab['observed_variable'][i],a_obstab['observation_value'],a_obstab['conversion_flag'],a_obstab['conversion_method'],temp,cdpddp,cdpdrh,crhdpd,cshrh,cshdpd)
                elif obstab['observed_variable'][i] in wvar:
                    wconvert(j,k,obstab['observed_variable'][i],a_obstab['observation_value'],a_obstab['conversion_flag'],a_obstab['conversion_method'],cuwind,cvwind,cwd,cws)
                    
        if humlist:
            for h in humvar:
                if h not in humlist:
                    j+=1
                    do_copy(a_obstab,obstab,j,i)
                    a_obstab['observed_variable'][j]=h
#                     a_obstab['observation_id'][j] = b'99' + obstab['observation_id'][j][2:]
                    qconvert(j,k,h,a_obstab['observation_value'],a_obstab['conversion_flag'],a_obstab['conversion_method'],temp,cdpddp,cdpdrh,crhdpd,cshrh,cshdpd)
                    if a_obstab['observation_value'][j]!=a_obstab['observation_value'][j]:
                        j-=1
            humlist.clear()
        if wlist:
            for h in wvar:
                if h not in wlist:
                    j+=1
                    do_copy(a_obstab,obstab,j,i)
                    a_obstab['observed_variable'][j]=h
#                     a_obstab['observation_id'][j] = b'99' + obstab['observation_id'][j][2:]
                    wconvert(j,k,h,a_obstab['observation_value'],a_obstab['conversion_flag'],a_obstab['conversion_method'],cuwind,cvwind,cwd,cws)
                    if a_obstab['observation_value'][j]!=a_obstab['observation_value'][j]:
                        j-=1
            wlist.clear()
        if idxu!=obstab['observation_value'].shape[0]:        
            if obstab['date_time'][idxu] != obstab['date_time'][idx[k]]:
                recordindex[ri]=j+1
                recordtimestamp[ri]=obstab['date_time'][idxu]
                ri+=1

        addedvar.append([i, j])
        if k%100000==0:
            print(k,idx.shape[0])
    j=j+1
    addedvar.append([i, j])
    return a_obstab, recordindex[:ri], recordtimestamp[:ri], j, addedvar
    
@njit
def fill_obsid(avar,conversion_flag):
    for o in range(avar.shape[0]):
        if conversion_flag[o] == 0:
            for i in range(2):
                avar[o,i]=b'9'
    return avar

@njit
def fill_restdata(final, rest_data, addedvar, j):
    fidx=0
    for o in addedvar:
        cc = o[0]
        final[fidx]=rest_data[cc]
        fidx+=1
        while fidx < o[1]:
            final[fidx]=rest_data[cc]
            fidx+=1
    while len(final) < j:
        final[fidx]=rest_data[cc]
        fidx+=1
    return final

def split(x): 
    return [i.encode() for i in x.decode()]

def convert_missing(fn, destination: str = opath):
    tt=time.time()
    nanlist = [float('nan'), np.nan, 0, -2147483648]
    with eua.CDMDataset(fn) as data:
        keys = data.observations_table.keys()
        keys = [x for x in keys if not x.startswith('string')]
        keys.remove('index')
        keys.remove('shape')
        obskeys = keys
        obstab_writetofile = [[] for i in range(len(obskeys))]
        
        keys = data.era5fb.keys()
#         keys = [x for x in keys if x in ['fg_depar@body','an_depar@body','biascorr@body','biascorr_fg@body']]
        keys = [x for x in keys if not x.startswith('string')]
        keys.remove('index')
        keys.remove('shape')
        fg_depar = keys.index('fg_depar@body')
        depar = keys.index('an_depar@body')
        biascorr = keys.index('biascorr@body')
        fg_biascorr = keys.index('biascorr_fg@body')
        fbkeys = keys
        fb_writetofile = [[] for i in range(len(fbkeys))]
        
        
        # loading data:
        loaded_data=[]
        a_loaded_data=[]
        loaded_type = {'names':[],'formats':[]}
        ld=[]
        for o in obskeys:
            if o in ['observed_variable','observation_value','z_coordinate','z_coordinate_type','date_time','observation_id','conversion_flag','conversion_method']:  
                if len(data.observations_table[o].shape)==1:
                    loaded_data.append((data.observations_table[o][:]))
                    a_loaded_data.append(numpy.empty_like(loaded_data[-1],shape=3*len(loaded_data[-1])))
                else:
                    loaded_data.append(data.observations_table[o][:].view('S{}'.format(data.observations_table[o].shape[1])).flatten())   
#                     a_loaded_data.append(numpy.empty((3*len(loaded_data[-1]),len(loaded_data[-1][0])), dtype=loaded_data[-1][0].dtype))
                    a_loaded_data.append(numpy.empty_like(loaded_data[-1],shape=3*len(loaded_data[-1])))
                    a_loaded_data[-1].fill(b' '*data.observations_table[o].shape[1])
                loaded_type['names'].append(o)
                loaded_type['formats'].append(loaded_data[-1].dtype)
                ld.append((o,loaded_data[-1].dtype))
        loaded_obstab = numpy.rec.fromarrays(loaded_data, dtype=ld)
        a_loaded_obstab = numpy.rec.fromarrays(a_loaded_data, dtype=ld)
        
        loaded_fb=[]
        a_loaded_fb=[]
        loaded_type = {'names':[],'formats':[]}
        lf=[]
        for o in fbkeys:
            if o in ['fg_depar@body','an_depar@body','biascorr@body','biascorr_fg@body']:  
                loaded_fb.append((data.era5fb[o][:]))
                a_loaded_fb.append(numpy.empty_like(loaded_fb[-1],shape=3*len(loaded_fb[-1])))
                loaded_type['names'].append(o)
                loaded_type['formats'].append(loaded_fb[-1].dtype)
                lf.append((o,loaded_fb[-1].dtype))
        loaded_feedback = numpy.rec.fromarrays(loaded_fb, dtype=lf)
        a_loaded_feedback = numpy.rec.fromarrays(a_loaded_fb, dtype=lf)
            
        recordindex = data.recordindex[:]
        # --->

    print(time.time()-tt)
    
    idx,press,temp,relhum,spechum,dpd,dewpoint,uwind,vwind,wd,ws=ipl2(loaded_obstab)
    
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

    idy=numpy.where(loaded_obstab['z_coordinate_type'][idx]==2) # do not convert humidity if data are not on pressure coordinates
    for c in cdpdrh,cshrh,cshdpd,crhdpd:
        c[idy]=numpy.nan
    
    cuwind = ws * np.cos(np.radians(wd))
    cvwind = ws * np.sin(np.radians(wd))
    cws = np.sqrt(uwind ** 2 + vwind ** 2)
    cwd = 90 - np.arctan2(-vwind, -uwind) * 180 / np.pi - 180.
    cwd = np.where(cwd > 0., cwd, 360.+cwd)

    humvar=numpy.array((34,36,38,39)) #dpd,dp,rh,sh
    wvar=numpy.array((104,105,106,107)) #dpd,dp,rh,sh
                       
    reduced_obskeys=List(loaded_obstab.dtype.fields.keys())
    reduced_fbkeys=List(loaded_feedback.dtype.fields.keys())
    out, ri, rt, j, addedvar=augment(loaded_obstab, a_loaded_obstab, reduced_obskeys,
                  idx,temp,press,relhum,spechum,dpd,dewpoint,uwind,vwind,wd,ws,
                  cdpddp,cdpdrh,cshrh,cshdpd,crhdpd,cuwind,cvwind,cwd,cws,humvar,wvar)
    avars = {}
    for i in reduced_obskeys:
        avars[i] = out[i][:j]                
        
    print(time.time()-tt)    
          
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
    
    allvars = copy.copy(avars['observed_variable'])
    allvars.sort()
    allvars = numpy.unique(allvars)
    #
    obsv = avars['observed_variable']
    #
    # resorting the data
    #
    @njit
    def make_vrindex(vridx,ridx,idx):
        l=0
        for i in range(1,len(idx)): # to set the recordindices
            if ridx[i]>ridx[i-1]:
                vridx[ridx[i-1]:ridx[i]]=l # next record after l
                l=i
        vridx[ridx[i]:]=len(idx) # next record for the last element is the len of the data

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
    absidx=np.concatenate(absidx)
                       
    # recordtimestamps are only necessary once
    recordtimestamps = rt 
    
    print('elapsed converting: ',time.time()-tt)
    tt=time.time()
    if os.path.isfile(targetfile):
        mode='r+'
    else:
        mode='w'
    
    for i in obskeys:
        print(i)
        print(time.time()-tt)
        
        if i == 'observation_id':
            ov_vars = copy.copy(avars[i])
#             ov_vars = numpy.array([b'99'+x[:2] for x in ov_vars[avars['conversion_flag'] == 0]])
#             for o in range(len(ov_vars)):
#                 if avars['conversion_flag'][o] == 0:
#                     ov_vars[o]=b'99'+ov_vars[o][:2]
#             ov_vars = ov_vars.view('S1').reshape((len(avars[i]),11))

            ov_vars=fill_obsid(ov_vars.view('S1').reshape((len(ov_vars),11)),avars['conversion_flag'])

        elif i in reduced_obskeys:
            ov_vars = avars[i]
            
        else: 
            with eua.CDMDataset(fn) as data:
                rest_data = data.observations_table[i][:]
            if i in ['observation_id', 'report_id', 'sensor_id', 'source_id']:
                final = numpy.empty((addedvar[-1][1],len(rest_data[0])), dtype=rest_data[0].dtype)
            else:
                final = numpy.empty(addedvar[-1][1], dtype=rest_data[0].dtype)
            ov_vars = fill_restdata(final, rest_data, addedvar, j)
        
        ov_vars = ov_vars[absidx]
        if i == 'index':
            pass
        elif i == 'observation_id' or i == 'report_id' or i == 'sensor_id' or i == 'source_id':
            alldict = {i:np.asarray(ov_vars, dtype='S1')}
            write_dict_h5(targetfile, alldict, 'observations_table', {i: { 'compression': 'gzip' } }, [i])
        else:
            alldict = pandas.DataFrame({i:ov_vars})
            write_dict_h5(targetfile, alldict, 'observations_table', {i: { 'compression': 'gzip' } }, [i])  

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

    print('elapsed writing '+targetfile+':',time.time()-tt)
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