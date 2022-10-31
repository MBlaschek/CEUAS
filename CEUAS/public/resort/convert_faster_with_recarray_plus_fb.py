#!/usr/bin/env
# coding: utf-8

from numba import njit
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
from harvest_convert_to_netCDF import write_dict_h5
import gc
import cds_eua3 as eua
#eua.logging_set_level(30)
import xarray as xr
import cdsapi, zipfile, os, time
import copy
from shutil import copyfile
import multiprocessing
sys.path.insert(0,os.getcwd()+'/../resort/rasotools-master/')
import rasotools
dir(rasotools)
import warnings
from functools import partial
import matplotlib.pylab as plt
warnings.filterwarnings('ignore')
import pickle

#opath='/raid60/scratch/uli/converted_v5/'
#opath='/raid60/scratch/leo/scratch/converted_v5/'
# if there are nan values in the pressure level - we will just sort without any converting!
def do_resort(fn):
    targetfile = opath+fn.split('/')[-1] 
    
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
    del obsv
    del dt
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

        if observed_variable[i]==ipar[38] and observation_value[i]==observation_value[i]:
            relhum[j-1]=observation_value[i]
        if observed_variable[i]==ipar[39] and observation_value[i]==observation_value[i]:
            if j<=spechum.shape[0]:
                spechum[j-1]=observation_value[i]
        if observed_variable[i]==ipar[34] and observation_value[i]==observation_value[i]:
            dpd[j-1]=observation_value[i]
        if observed_variable[i]==ipar[36] and observation_value[i]==observation_value[i] and z_coordinate[i]==z_coordinate[i]:
            dewpoint[j-1]=observation_value[i]
        if observed_variable[i]==85 and observation_value[i]==observation_value[i]:
            temp[j-1]=observation_value[i]
            good=True

    return idx,press,temp,relhum,spechum,dpd,dewpoint

@njit(boundscheck=True)
def ipl2(lobs, fb):
    
    observed_variable=lobs['observed_variable']
    observation_value=lobs['observation_value']
    z_coordinate=lobs['z_coordinate']
    z_coordinate_type=lobs['z_coordinate_type']
    recordtimestamp=lobs['date_time']
    
    departure=fb['an_depar@body']
    fg_departure=fb['fg_depar@body']
    
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
    
    d_temp=numpy.empty(jdx,dtype=observation_value.dtype)
    d_relhum=numpy.empty(jdx,dtype=observation_value.dtype)
    d_spechum=numpy.empty(jdx,dtype=observation_value.dtype)
    d_dpd=numpy.empty(jdx,dtype=observation_value.dtype)
    d_dewpoint=numpy.empty(jdx,dtype=observation_value.dtype)
    d_uwind=numpy.empty(jdx,dtype=observation_value.dtype)
    d_vwind=numpy.empty(jdx,dtype=observation_value.dtype)
    d_wd=numpy.empty(jdx,dtype=observation_value.dtype)
    d_ws=numpy.empty(jdx,dtype=observation_value.dtype)
    
    fgd_temp=numpy.empty(jdx,dtype=observation_value.dtype)
    fgd_relhum=numpy.empty(jdx,dtype=observation_value.dtype)
    fgd_spechum=numpy.empty(jdx,dtype=observation_value.dtype)
    fgd_dpd=numpy.empty(jdx,dtype=observation_value.dtype)
    fgd_dewpoint=numpy.empty(jdx,dtype=observation_value.dtype)
    fgd_uwind=numpy.empty(jdx,dtype=observation_value.dtype)
    fgd_vwind=numpy.empty(jdx,dtype=observation_value.dtype)
    fgd_wd=numpy.empty(jdx,dtype=observation_value.dtype)
    fgd_ws=numpy.empty(jdx,dtype=observation_value.dtype)
    
    for v in temp,relhum,dewpoint,dpd,spechum,uwind,vwind,wd,ws,d_temp,d_relhum,d_spechum,d_dpd,d_dewpoint,d_uwind,d_vwind,d_wd,d_ws,fgd_temp,fgd_relhum,fgd_spechum,fgd_dpd,fgd_dewpoint,fgd_uwind,fgd_vwind,fgd_wd,fgd_ws:
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
    
        if observed_variable[i]==ipar[38]:
            relhum[j]=observation_value[i]
            d_relhum[j]=observation_value[i]-departure[i]
            fgd_relhum[j]=observation_value[i]-fg_departure[i]
        elif observed_variable[i]==ipar[39]:
            spechum[j]=observation_value[i]
            d_spechum[j]=observation_value[i]-departure[i]
            fgd_spechum[j]=observation_value[i]-fg_departure[i]
        elif observed_variable[i]==ipar[34]:
            dpd[j]=observation_value[i]
            d_dpd[j]=observation_value[i]-departure[i]
            fgd_dpd[j]=observation_value[i]-fg_departure[i]
        elif observed_variable[i]==ipar[36]:
            dewpoint[j]=observation_value[i]
            d_dewpoint[j]=observation_value[i]-departure[i]
            fgd_dewpoint[j]=observation_value[i]-fg_departure[i]
        elif observed_variable[i]==85:
            temp[j]=observation_value[i]
            d_temp[j]=observation_value[i]-departure[i]
            fgd_temp[j]=observation_value[i]-fg_departure[i]
        elif observed_variable[i]==ipar[104]:
            uwind[j]=observation_value[i]
            d_uwind[j]=observation_value[i]-departure[i]
            fgd_uwind[j]=observation_value[i]-fg_departure[i]
        elif observed_variable[i]==ipar[105]:
            vwind[j]=observation_value[i]
            d_vwind[j]=observation_value[i]-departure[i]
            fgd_vwind[j]=observation_value[i]-fg_departure[i]
        elif observed_variable[i]==ipar[106]:
            wd[j]=observation_value[i]
            d_wd[j]=observation_value[i]-departure[i]
            fgd_wd[j]=observation_value[i]-fg_departure[i]
        elif observed_variable[i]==ipar[107]:
            ws[j]=observation_value[i]
            d_ws[j]=observation_value[i]-departure[i]
            fgd_ws[j]=observation_value[i]-fg_departure[i]
        else:
            pass
                
    #return        
    print(j,jdx)
    return idx,press,temp,relhum,spechum,dpd,dewpoint,uwind,vwind,wd,ws,d_temp,d_relhum,d_spechum,d_dpd,d_dewpoint,d_uwind,d_vwind,d_wd,d_ws,fgd_temp,fgd_relhum,fgd_spechum,fgd_dpd,fgd_dewpoint,fgd_uwind,fgd_vwind,fgd_wd,fgd_ws

@njit(boundscheck=True)
def qconvert(j,k,h,a_observation_value,a_conversion_flag,a_conversion_method,
             a_an_depar,a_fg_depar,
             temp,cdpddp,cdpdrh,crhdpd,cshrh,cshdpd,crhsh,cdpdsh,
             d_cdpddp,d_cdpdrh,d_cdpdsh,d_cshrh,d_cshdpd,d_crhsh,d_crhdpd,
             fgd_cdpddp,fgd_cdpdrh,fgd_cdpdsh,fgd_cshrh,fgd_cshdpd,fgd_crhsh,fgd_crhdpd):
    if h==ipar[34]:
        if cdpddp[k]==cdpddp[k]:
            a_observation_value[j]=cdpddp[k]
            a_an_depar[j]=cdpddp[k]-d_cdpddp[k]
            a_fg_depar[j]=cdpddp[k]-fgd_cdpddp[k]
            if (numpy.abs(cdpddp[k])>80) or (cdpddp[k]<0.01):
                a_observation_value[j]=numpy.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=2
        elif cdpdrh[k]==cdpdrh[k]:
            a_observation_value[j]=cdpdrh[k]
            a_an_depar[j]=cdpdrh[k]-d_cdpdrh[k]
            a_fg_depar[j]=cdpdrh[k]-fgd_cdpdrh[k]
            if (numpy.abs(cdpdrh[k])>80) or (cdpdrh[k]<0.01):
                a_observation_value[j]=numpy.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=3
        else:
            a_observation_value[j]=cdpdsh[k]
            a_an_depar[j]=cdpdsh[k]-d_cdpdsh[k]
            a_fg_depar[j]=cdpdsh[k]-fgd_cdpdsh[k]
            if (numpy.abs(cdpdsh[k])>80) or (cdpdsh[k]<0.01):
                a_observation_value[j]=numpy.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=4
            
    elif h==ipar[36]:
        if cdpdrh[k]==cdpdrh[k]:
            a_observation_value[j]=temp[k]-cdpdrh[k]
            a_an_depar[j]=(temp[k]-cdpdrh[k])-(temp[k]-d_cdpdrh[k])
            a_fg_depar[j]=(temp[k]-cdpdrh[k])-(temp[k]-fgd_cdpdrh[k])
            if (numpy.abs(cdpdrh[k])>80) or (cdpdrh[k]<0.01):
                a_observation_value[j]=numpy.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=3
        elif cdpdrh[k]==cdpdrh[k]:
            a_observation_value[j]=temp[k]-cdpddp[k]
            a_an_depar[j]=(temp[k]-cdpddp[k])-(temp[k]-d_cdpddp[k])
            a_fg_depar[j]=(temp[k]-cdpddp[k])-(temp[k]-fgd_cdpddp[k])
            if (numpy.abs(cdpddp[k])>80) or (cdpddp[k]<0.01):
                a_observation_value[j]=numpy.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=2
        else:
            a_observation_value[j]=temp[k]-cdpdsh[k]
            a_an_depar[j]=(temp[k]-cdpdsh[k])-(temp[k]-d_cdpdsh[k])
            a_fg_depar[j]=(temp[k]-cdpdsh[k])-(temp[k]-fgd_cdpdsh[k])
            if (numpy.abs(cdpdsh[k])>80) or (cdpdsh[k]<0.01):
                a_observation_value[j]=numpy.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=4
        
    elif h==ipar[38]:
        if crhdpd[k]==crhdpd[k]:
            a_observation_value[j]=crhdpd[k]
            a_an_depar[j]=crhdpd[k]-d_crhdpd[k]
            a_fg_depar[j]=crhdpd[k]-fgd_crhdpd[k]
            if (crhdpd[k]<0.) or (crhdpd[k]>1.03):
                a_observation_value[j]=numpy.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=2
        else: 
            a_observation_value[j]=crhsh[k]
            a_an_depar[j]=crhsh[k]-d_crhsh[k]
            a_fg_depar[j]=crhsh[k]-fgd_crhsh[k]
            if (crhsh[k]<0.) or (crhsh[k]>1.03):
                a_observation_value[j]=numpy.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=4
            
    elif h==ipar[39]:
        if cshrh[k]==cshrh[k]:
            a_observation_value[j]=cshrh[k]
            a_an_depar[j]=cshrh[k]-d_cshrh[k]
            a_fg_depar[j]=cshrh[k]-fgd_cshrh[k]
            if (cshrh[k]<0.) or (cshrh[k]>50.):
                a_observation_value[j]=numpy.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=3
        else:
            a_observation_value[j]=cshdpd[k]
            a_an_depar[j]=cshdpd[k]-d_cshdpd[k]
            a_fg_depar[j]=cshdpd[k]-fgd_cshdpd[k]
            if (cshdpd[k]<0.) or (cshdpd[k]>50.):
                a_observation_value[j]=numpy.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=2
    else:
        print('qconvert called with wrong variable')
        
    return

@njit(boundscheck=True)
def wconvert(j,k,h,a_observation_value,a_conversion_flag,a_conversion_method,
             a_an_depar,a_fg_depar,
             cuwind,cvwind,cwd,cws,
             d_cwd,d_cws,
             fgd_cwd,fgd_cws):
    if h==ipar[104]:
        if cuwind[k]==cuwind[k]:
            a_observation_value[j]=cuwind[k]
            a_an_depar[j]=numpy.nan
            a_fg_depar[j]=numpy.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=1
    elif h==ipar[105]:
        if cvwind[k]==cvwind[k]:
            a_observation_value[j]=cvwind[k]
            a_an_depar[j]=numpy.nan
            a_fg_depar[j]=numpy.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=1
    elif h==ipar[106]:
        if cwd[k]==cwd[k]:
            a_observation_value[j]=cwd[k]
            a_an_depar[j]=d_cwd[k]
            a_fg_depar[j]=fgd_cwd[k]
            if (cwd[k]<0.) or (cwd[k]>360.):
                a_observation_value[j]=numpy.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=2
    elif h==ipar[107]:
        if cws[k]==cws[k]:
            a_observation_value[j]=cws[k]
            a_an_depar[j]=d_cws[k]
            a_fg_depar[j]=fgd_cws[k]
            a_conversion_flag[j]=0
            a_conversion_method[j]=2
    else:
        print('wconvert called with wrong variable')

@njit
def do_copy(a_obstab,obstab,j,i):
    # all written here - * will be overwritten if it's a converted variable
    a_obstab['date_time'][j]=obstab['date_time'][i]
    a_obstab['observation_id'][j]=obstab['observation_id'][i] # *
    #a_obstab['report_id'][j]=obstab['report_id'][i] 
    a_obstab['observation_value'][j]=obstab['observation_value'][i] # *
    a_obstab['observed_variable'][j]=obstab['observed_variable'][i] # *
    a_obstab['z_coordinate'][j]=obstab['z_coordinate'][i]
    a_obstab['z_coordinate_type'][j]=obstab['z_coordinate_type'][i]
    a_obstab['conversion_flag'][j]=obstab['conversion_flag'][i] # *
    a_obstab['conversion_method'][j]=obstab['conversion_method'][i] # *
    return

@njit
def do_fb_copy(a_loaded_feedback,loaded_feedback,j,i):
    # all written here - * will be overwritten if it's a converted variable
    a_loaded_feedback['fg_depar@body'][j]=loaded_feedback['fg_depar@body'][i] # *
    a_loaded_feedback['an_depar@body'][j]=loaded_feedback['an_depar@body'][i] # *
    a_loaded_feedback['biascorr@body'][j]=loaded_feedback['biascorr@body'][i]
    a_loaded_feedback['biascorr_fg@body'][j]=loaded_feedback['biascorr_fg@body'][i]
    return

@njit(boundscheck=True)          
def augment(obstab, a_obstab, loaded_feedback, a_loaded_feedback,
            idx,temp,press,relhum,spechum,dpd,dewpoint,uwind,vwind,wd,ws,
            cdpddp,cdpdrh,cshrh,cshdpd,crhdpd,crhsh,cdpdsh,cuwind,cvwind,cwd,cws,
            d_cdpddp,d_cdpdrh,d_cdpdsh,d_cshrh,d_cshdpd,d_crhsh,d_crhdpd,
            fgd_cdpddp,fgd_cdpdrh,fgd_cdpdsh,fgd_cshrh,fgd_cshdpd,fgd_crhsh,fgd_crhdpd,
            d_cwd,d_cws,fgd_cwd,fgd_cws,
            humvar,wvar):
    
    
    recordindex=numpy.empty(idx.shape[0],obstab['date_time'].dtype)
    recordtimestamp=numpy.empty(idx.shape[0],obstab['date_time'].dtype)
    
    j=-1 # augmented index
    jsave=0
    p= obstab['z_coordinate'][0]-1 # z_coordinate[0]-1
    rts= obstab['date_time'][0] # date_time[0]
    
    addedvar=numpy.zeros((a_obstab['observation_value'].shape[0],2),dtype=numpy.int32)
    humlist=List()
    wlist=List()
    recordindex[0]=0
    recordtimestamp[0]=obstab['date_time'][0] # date_time[0]
    ri=1
    l=1
    for k in range(idx.shape[0]):  # no new levels are created, but variables per level will increase
        jsave=j
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
                do_fb_copy(a_loaded_feedback,loaded_feedback,j,i)
                if obstab['observed_variable'][i] in humvar:
                    humlist.append(obstab['observed_variable'][i])
                elif obstab['observed_variable'][i] in wvar:
                    wlist.append(obstab['observed_variable'][i])
            else:
                do_copy(a_obstab,obstab,j,i)
                if obstab['observed_variable'][i] in humvar:
                    a_loaded_feedback['biascorr@body'][j]=numpy.nan
                    a_loaded_feedback['biascorr_fg@body'][j]=numpy.nan
                    qconvert(j,k,obstab['observed_variable'][i],a_obstab['observation_value'],a_obstab['conversion_flag'],a_obstab['conversion_method'],
                             a_loaded_feedback['an_depar@body'],a_loaded_feedback['fg_depar@body'],
                             temp,cdpddp,cdpdrh,crhdpd,cshrh,cshdpd,crhsh, cdpdsh, 
                             d_cdpddp,d_cdpdrh,d_cdpdsh,d_cshrh,d_cshdpd,d_crhsh,d_crhdpd,
                             fgd_cdpddp,fgd_cdpdrh,fgd_cdpdsh,fgd_cshrh,fgd_cshdpd,fgd_crhsh,fgd_crhdpd)
                elif obstab['observed_variable'][i] in wvar:
                    a_loaded_feedback['biascorr@body'][j]=numpy.nan
                    a_loaded_feedback['biascorr_fg@body'][j]=numpy.nan
                    wconvert(j,k,obstab['observed_variable'][i],a_obstab['observation_value'],a_obstab['conversion_flag'],a_obstab['conversion_method'],
                             a_loaded_feedback['an_depar@body'],a_loaded_feedback['fg_depar@body'],
                             cuwind,cvwind,cwd,cws,
                             d_cwd,d_cws,
                             fgd_cwd,fgd_cws)
                    #lens=np.zeros(4)
                    #for ii in range(104,108):
                        #idy=np.where(a_obstab['observed_variable'][:j+1]==ii)[0]
                        #lens[ii-104]=len(idy)
                        #print(ii,lens[ii-104])
                        #if ii>104 and lens[ii-104]!=lens[ii-105]:
                            #print('inconsistent')
                    
        if humlist:
            for h in humvar:
                if h not in humlist:
                    j+=1
                    do_copy(a_obstab,obstab,j,i)
                    a_loaded_feedback['biascorr@body'][j]=numpy.nan
                    a_loaded_feedback['biascorr_fg@body'][j]=numpy.nan
                    a_obstab['observed_variable'][j]=h
                    qconvert(j,k,h,a_obstab['observation_value'],a_obstab['conversion_flag'],a_obstab['conversion_method'],
                             a_loaded_feedback['an_depar@body'],a_loaded_feedback['fg_depar@body'],
                             temp,cdpddp,cdpdrh,crhdpd,cshrh,cshdpd,crhsh, cdpdsh, 
                             d_cdpddp,d_cdpdrh,d_cdpdsh,d_cshrh,d_cshdpd,d_crhsh,d_crhdpd,
                             fgd_cdpddp,fgd_cdpdrh,fgd_cdpdsh,fgd_cshrh,fgd_cshdpd,fgd_crhsh,fgd_crhdpd)
                    if a_obstab['observation_value'][j]!=a_obstab['observation_value'][j]:
                        j-=1
            humlist.clear()
        if len(wlist)>0:
            for h in wvar:
                if h not in wlist:
                    j+=1
                    do_copy(a_obstab,obstab,j,i)
                    a_loaded_feedback['biascorr@body'][j]=numpy.nan
                    a_loaded_feedback['biascorr_fg@body'][j]=numpy.nan
                    a_obstab['observed_variable'][j]=h
                    wconvert(j,k,h,a_obstab['observation_value'],a_obstab['conversion_flag'],a_obstab['conversion_method'],
                             a_loaded_feedback['an_depar@body'],a_loaded_feedback['fg_depar@body'],
                             cuwind,cvwind,cwd,cws,
                             d_cwd,d_cws,
                             fgd_cwd,fgd_cws)
                    if a_obstab['observation_value'][j]!=a_obstab['observation_value'][j]:
                        j-=1
            #lens=np.zeros(4)
            #for ii in range(104,108):
                #idy=np.where(a_obstab['observed_variable'][jsave:j+1]==ii)[0]
                #lens[ii-104]=len(idy)
                #print(ii,lens[ii-104])
                #if ii>104 and lens[ii-104]!=lens[ii-105]:
                    #print('inconsistent')
            #jsave=j
            wlist.clear()
        if idxu!=obstab['observation_value'].shape[0]:        
            if obstab['date_time'][idxu] != obstab['date_time'][idx[k]]:
                recordindex[ri]=j+1
                recordtimestamp[ri]=obstab['date_time'][idxu]
                ri+=1
        else:
            print('spurious idxu')

        addedvar[l,0]=i
        addedvar[l,1]=j#.append([i, j])
        l+=1
        if k%100000==0:
            print(k,idx.shape[0])
    j=j+1
    addedvar[l,0]=i+1
    addedvar[l,1]=j#.append([i, j])
    l+=1
    #addedvar[i]=j#.append([i, j])
    #addv=numpy.empty((len(addedvar),2),dtype=numpy.int64)
    #for k in range(len(addedvar)):
        #for l in range(2):
            #addv[k,l]=addedvar[k][l]
    return a_obstab, a_loaded_feedback, recordindex[:ri], recordtimestamp[:ri], j, addedvar[:l]

    
@njit
def fill_obsid(avar,conversion_flag):
    for o in range(avar.shape[0]):
        if conversion_flag[o] == 0:
            for i in range(2):
                avar[o,i]=b'9'
    return avar

@njit()
def fill_restdata(final, rest_data, addedvar, j): #,test):

    for l in range(addedvar.shape[0]-1):
        diff=(addedvar[l+1,1]-addedvar[l,1])-(addedvar[l+1,0]-addedvar[l,0])
        final[addedvar[l,1]:addedvar[l+1,1]-diff]=rest_data[addedvar[l,0]:addedvar[l+1,0]]
        if diff>0:
            for i in range(diff):
                final[addedvar[l+1,1]-diff+i]=rest_data[addedvar[l+1,0]-1]
        #x=numpy.nansum(final[addedvar[l,1]:addedvar[l+1,1]]-test[addedvar[l,1]:addedvar[l+1,1]])
        #if ~numpy.isnan(x) and x!=0:
                     #print(addedvar[l:l+2])                             
                     #print('x')
        
    return final

def split(x): 
    return [i.encode() for i in x.decode()]

def offline_fb(fpattern,p,lat,lon,refs,ans):
    tt=time.time()
    try:

        with h5py.File(fpattern.format(ans[0],ans[1]),'r') as g:
            pres=g['level'][:]
            pidx=numpy.searchsorted(pres,refs['level'])
            dx=360/g[p].shape[1]
            dy=180/(g[p].shape[0]-1)
            ix=int(np.floor(lon/dx))  
            if lat==90.:
                lat=89.9999
            iy=int(np.floor((90-lat)/dy))
            if ix==g[p].shape[1]-1:
                tera5=np.empty_like(g['t'][iy:iy+2,ix:ix+2,:,:])
                tera5[0,:,:,:]=g['t'][iy:iy+2,ix,:,:]
                tera5[1,:,:,:]=g['t'][iy:iy+2,0,:,:]
                if '20CRv3' in fpattern:
                    for it in range(tera5.shape[2]):                       
                        for ip in range(len(pidx)):
                            tera5[0,:,it,pidx[ip]]+=refs[p][iy:iy+2,ix,ip]
                            tera5[1,:,it,pidx[ip]]+=refs[p][iy:iy+2,0,ip]
            else:
                tera5=g[p][iy:iy+2,ix:ix+2,:,:]
                if '20CRv3' in fpattern:
                    for it in range(tera5.shape[2]):                       
                        for ip in range(len(pidx)):
                            tera5[:,:,it,pidx[ip]]+=refs[p][iy:iy+2,ix:ix+2,ip]
            try:
                mv=np.where(tera5==g[p].attrs['missing_value'])
                tera5=tera5*np.float32(g[p].attrs['scale_factor'])+np.float32(g[p].attrs['add_offset'])
                tera5[mv]=np.nan
            except:
                pass
            attribs={}
            for k,v in g[p].attrs.items():
                if k not in ['scale_factor','add_offset','DIMENSION_LIST','_FillValue','name','long_name','standard_name','units']:
                    attribs[k]=v
                    if k=='missing_value':
                        attribs[k]=np.nan
            for k,v in g.attrs.items():
                attribs[k]=v
            attribs['source_filename']=fpattern

            tunits=g['time'].attrs['units']
            try:
                tunits=tunits.decode('latin1')
            except:
                pass
            tunits=tunits.split()

            if tunits[0]=='hours':
                secs=np.int64(g['time'][:])*3600
            if tunits[-2]!='1900-01-01':
                offset=datetime.strptime(tunits[-2],'%Y-%m-%d')-datetime(1900,1,1)
                secs=secs+int(offset.total_seconds())

        if '20CRv3' in fpattern:
            tera5
        w=np.zeros(4)            
        w[0]=lon-ix*dx
        w[1]=1.-w[0]
        w[3]=lat-(90-iy*dy)
        w[2]=1.-w[3]
        if tera5.shape[0]==1:     
            tan=w[0]*tera5[0,1,:,:]+w[1]*tera5[0,0,:,:]
        else:
            tan=w[0]*w[2]*tera5[1,1,:,:]+w[1]*w[2]*tera5[1,0,:,:]+w[0]*w[3]*tera5[0,1,:,:]+w[1]*w[3]*tera5[0,0,:,:]

        #print(ans,time.time()-tt)
    except:
        return 
    return tan,secs,pres,attribs

def Interp2d(datain,xin,yin,xout,yout,order=0):

    """
              Interpolates a 2D array onto a new grid (only works for linear grids), 
              with the Lat/Lon inputs of the old and new grid. Can perfom nearest
              neighbour interpolation or bilinear interpolation (of order 1)'

       This is an extract from the basemap module (truncated)
    """

    if order==0:
           interpolation='NearestNeighbour' 
    else:
           interpolation = 'Bilinear'
           
    # Mesh Coordinates so that they are both 2D arrays
 #   xout,yout = numpy.meshgrid(xout,yout)

   # compute grid coordinates of output grid.
    delx = xin[1:]-xin[0:-1]
    dely = yin[1:]-yin[0:-1]

    xcoords = (len(xin)-1)*(xout-xin[0])/(xin[-1]-xin[0])
    ycoords = (len(yin)-1)*(yout-yin[0])/(yin[-1]-yin[0])


    xcoords = numpy.clip(xcoords,0,len(xin)-1)
    ycoords = numpy.clip(ycoords,0,len(yin)-1)

    ## Interpolate to output grid using nearest neighbour
    if interpolation == 'NearestNeighbour':
        pass
        #xcoordsi = numpy.around(xcoords).astype(numpy.int32)
        #ycoordsi = numpy.around(ycoords).astype(numpy.int32)
        #dataout = datain[ycoordsi,xcoordsi]

    # Interpolate to output grid using bilinear interpolation.
    elif interpolation == 'Bilinear':
        xi = xcoords.astype(numpy.int32)
        yi = ycoords.astype(numpy.int32)
        xip1 = xi+1
        yip1 = yi+1
        xip1 = numpy.clip(xip1,0,len(xin)-1)
        yip1 = numpy.clip(yip1,0,len(yin)-1)
        delx = xcoords-xi.astype(numpy.float32)
        dely = ycoords-yi.astype(numpy.float32)
        dataout = (1.-delx)*(1.-dely)*datain[yi,xi] + \
                  delx*dely*datain[yip1,xip1] + \
                  (1.-delx)*dely*datain[yip1,xi] + \
                  delx*(1.-dely)*datain[yi,xip1]

    return dataout

@njit(cache=True)
def interp(reatab,obstype,z,ts,press,iobstype,times,values):
    print('x')
    ipress=0
    itimesold=0
    tsold=0
    noval=0
    dt=times[1]-times[0]
    for i in range(obstype.shape[0]):
        if obstype[i]==iobstype:
            if ts[i]!=tsold:
                ipress=0
                tsold=ts[i]
            zz=z[i]/100
            while press[ipress]<zz:
                ipressold=ipress
                if ipress==press.shape[0]-1:
                    break
                ipress+=1
            if press[ipress]==zz:
                itimes=itimesold
                while times[itimes]<ts[i]:
                    itimesold=itimes
                    if itimes==times.shape[0]-1:
                        break
                    itimes+=1
                if itimes>0 and times[itimes]>=ts[i]:
                    w0=(times[itimes]-ts[i])/dt
                    w1=1.-w0
                    reatab[i]=values[itimes-1,ipress]*w0+values[itimes,ipress]*w1
                    if obstype[i]==85 and z[i]==10000 and reatab[i]>400.:
                        print(i,reatab[i],itimes,w0,w1)
                else:
                    noval+=1
            else:
                #if ipress==press.shape[0]-1:
                ipress=0

    if noval>0:
        print('no reanalysis data for',noval,'obsvalues' )
    return

def retrieve_anfg(fn,out_name,path_to_gridded):


    #era5=refdiff('/raid60/scratch/leo/scratch/ERA5/gridded/','era5fc.{0}{1:02}.'+par['gid'],1979,1983,tidx,fieldsperday=2,
                 #fcstep=12,mrange=list(range(1,13)),tgrid=None,tempname=tempname)
    
    #jra55_58_c[im,ipar,ip,:,:]=Interp2d(jra55_58_m[im,ipar,ip,:,:], jlon, jlat, eelon, eelat, order=1)
    
    readict={'era5':{'ftype':('t','fct'),'param':{'t':'130','u':'131','v':'132','q':'133','z':'129'},
                       'path':os.path.expandvars('$RSCRATCH/era5/gridded/'),'prefix':'era5','suffix':'','glue':'.'},
               #'CERA20C':{'ftype':('t',),'param':{'t':'130','u':'131','v':'132','q':'133','z':'129'},
                          #'path':os.path.expandvars('$RSCRATCH/CERA20C/'),'prefix':'CERA20C','suffix':'','glue':'.'},
               #'JRA55':{'ftype':('fcst_mdl',),'param':{'t':'011_tmp','u':'033_ugrd','v':'034_vgrd','q':'051_spfh'},
                        #'path':os.path.expandvars('$RSCRATCH/JRA55/split/'),'prefix':'test.','suffix':'grb','glue':'.'},
               '20CRv3':{'ftype':('',),'param':{'t':'TMP','u':'UGRD','v':'VGRD'},#,'q':'SPFH'},#,'z':'HGT'},
                         'path':os.path.expandvars('$RSCRATCH/20CRv3/'),'prefix':'anl_meant','suffix':'_pres','glue':'_'},
               }
    try:
        with open(os.path.expandvars(wpath+'/rea/refs.pkl'),'rb') as f:
            refs=pickle.load(f)
    except:
        
        refs={}
        raw_20CR={}
        tt=time.time()
        l=0
        for iy in range(1950,1952):
            for im in range(1,13):
                l+=1
                k='era5'
                v=readict['era5']
                k20='20CRv3'
                v20=readict['20CRv3']
                for k,v in readict.items():
                    if k=='20CRv3':
                        continue
                    for ftype in v['ftype']:  
                        if 'fc' in ftype:
                            dtype='fc'
                            continue
                        else:
                            dtype='an'
                        if dtype not in readict[k].keys():
                            readict[k][dtype]={}
            
                        for p,pn in v['param'].items():
                            if k!='JRA55':
            
                                fpattern=v['path']+v['prefix']+ftype+v['glue']+'{}{:0>2}'+v['glue']+pn+v['suffix']+'.nc'
                            else:
                                fpattern=v['path']+v['prefix']+ftype+v['glue']+pn+'.reg_tl319.{}{:0>2}'+v['glue']+v['suffix']+'.nc'
                    
                            try:
                                with h5py.File(fpattern.format(iy,im),'r') as f:
                                    print(f.keys())
                                    ifield=f[p][:]*np.float32(f[p].attrs['scale_factor'])+np.float32(f[p].attrs['add_offset'])
                                    ilon=f['longitude'][:]
                                    ilat=f['latitude'][:]
                                    iplev=f['level'][:]

                                fpattern20=v20['path']+v20['prefix']+v20['glue']+'{}{:0>2}'+v20['glue']+v20['param'][p]+v20['suffix']+'.nc'
                                with h5py.File(fpattern20.format(iy,im),'r') as f:
                                    print(f.keys())
                                    oplev=f['level'][:]
                                    pindex=np.searchsorted(oplev,iplev)
                                    oshape=f[p].shape
                                    if p not in refs.keys():
                                        refs[p]=np.zeros((oshape[0],oshape[1],ifield.shape[3]),dtype=ifield.dtype)
                                        raw_20CR[p]=np.zeros_like(refs[p])
                                        refs['level']=iplev[:]
                                    raw_20CR[p]+=np.mean(f[p][:],axis=2)[:,:,pindex]
                                    olon,olat=numpy.meshgrid(f['longitude'][:],f['latitude'][:])
                                
                                    
                                for ip in range(ifield.shape[3]):
                                    
                                    refs[p][:,:,ip]+=Interp2d(np.mean(ifield[:,:,:,ip],axis=2),ilon,ilat,olon,olat, order=1)
                                print(iy,im,p,np.std(refs[p]),np.std(ifield),time.time()-tt)
                                         
                            except Exception as e:
                                print(e,'could not read')
    
        for r in ['t','u','v']:
            refs[r]/=l
            raw_20CR[r]/=l
            refs[r]-=raw_20CR[r]
        
        try:
            os.mkdir(wpath+'/rea')
        except:
            pass
        with open(os.path.expandvars(wpath+'/rea/refs.pkl'),'wb') as f:
            pickle.dump(refs,f)
            
    #P=multiprocessing.Pool(12)
    with h5py.File(fn,'r') as f:   
        lat=f['header_table']['latitude'][0]
        lon=f['header_table']['longitude'][0]
        obstype=f['observations_table']['observed_variable'][:]
        obs=f['observations_table']['observation_value'][:]
        z=f['observations_table']['z_coordinate'][:]
        ts=f['observations_table']['date_time'][:]
        ofb=True
        try:
            
            o_minus_bg=f['era5fb']['fg_depar@body'][:]
            o_minus_an=f['era5fb']['an_depar@body'][:]
        except:
            ofb=False
    tsu=np.unique(ts)

    ref=datetime(1900,1,1)
    yms=[]
    oldyear=0
    oldmonth=0
    for k in tsu:
        x=ref+timedelta(seconds=int(k))
        if x.year!=oldyear or x.month!=oldmonth:
            yms.append((x.year,x.month))
            oldmonth=x.month
            oldyear=x.year

    for k,v in readict.items():
        #if k!='20CRv3':
            #continue
        for ftype in v['ftype']:  
            if 'fc' in ftype:
                dtype='fc'
            else:
                dtype='an'
            if dtype not in readict[k].keys():
                readict[k][dtype]={}
                readict[k][dtype]['refvalues']=np.empty(obs.shape,dtype=np.float32)
                readict[k][dtype]['refvalues'].fill(np.nan)

            for p,pn in v['param'].items():
                if k!='JRA55':

                    fpattern=v['path']+v['prefix']+ftype+v['glue']+'{}{:0>2}'+v['glue']+pn+v['suffix']+'.nc'
                else:
                    fpattern=v['path']+v['prefix']+ftype+v['glue']+pn+'.reg_tl319.{}{:0>2}'+v['glue']+v['suffix']+'.nc'
                #found=len(glob.glob(fpattern.format(*yms[0])))
                #print(fpattern.format(*yms[0]),found)
                #if found:
                func=partial(offline_fb,fpattern,p,lat,lon,refs)
                #tups=list(P.map(func,yms[:]))
                tups=list(map(func,yms[:])) # no multiprocessing because function is already mapped
                ntups=[]
                for t in tups:
                    if t is not None:
                        ntups.append(t)
                if len(ntups)>0:

                    press=ntups[0][2]
                    if np.max(press)>1000:
                        press/=100.

                    readict[k][dtype][p]={}
                    readict[k][dtype][p]['time']=np.concatenate([ntups[i][1] for i in range(len(ntups))])
                    readict[k][dtype][p]['values']=np.concatenate([ntups[i][0] for i in range(len(ntups))])
                    readict[k][dtype][p]['attribs']=ntups[0][3]
            obstypes={'t':ipar[85],'u':ipar[104],'v':ipar[105],'q':ipar[39],'z':ipar[0]}    
            #tr[1]=117  # should change
            #tr[2]=85
            #tr[3]=104
            #tr[4]=105
            #tr[7]=39 #spec hum
            #tr[29]=38 #relative hum
            #tr[59]=36 # dew point
            #tr[111]=106 #dd
            #tr[112]=107  #ff

            if dtype in readict[k].keys():  
                tt=time.time()
                for p in readict[k][dtype].keys():
                    if p!='refvalues':
                        try:

                            interp(readict[k][dtype]['refvalues'],obstype,z,ts,press,obstypes[p],
                                                 readict[k][dtype][p]['time'],readict[k][dtype][p]['values'])
                            print(k,dtype,p,time.time()-tt)
                        except Exception as e:
                            print(e)
                idx=np.where(np.logical_and(obstype==85,z==10000))
                plt.subplot(2,1,1)
                plt.plot(1900+ts[idx]/365.25/86400,obs[idx])
                plt.plot(1900+ts[idx]/365.25/86400,readict[k][dtype]['refvalues'][idx])
                plt.title(k)
                plt.subplot(2,1,2)
                plt.plot(1900+ts[idx]/365.25/86400,obs[idx]-readict[k][dtype]['refvalues'][idx])
                plt.title('obs -'+k+', rms= {:5.3f}'.format(np.sqrt(np.nanmean((obs[idx]-readict[k][dtype]['refvalues'][idx])**2))))
                plt.tight_layout()
                fnp=fn.split('/')[-1].split('CEUAS_merged_v1.nc')[0]
                plt.savefig(os.path.expanduser('~/tmp/'+fnp+k+dtype+'.png'))
                plt.close()

                df = {dtype:readict[k][dtype]['refvalues']}  # making a 1 column dataframe  
                print('writing',k)
                try:

                    write_dict_h5(out_name, df, k, {dtype: {'compression': 'gzip'}}, 
                                            var_selection=[], mode='a', attrs = {dtype:readict[k][dtype]['t']['attribs']} )  
                except Exception as e:
                    print(e, 'no values from ',k,'for station ',os.path.basename(fn))


    if not ofb:
        return readict
    
    readict['era5fb']={'ftype':['an','fc']}                  
    readict['era5fb']['an']={'refvalues':obs-o_minus_an }                  
    readict['era5fb']['fc']={'refvalues':obs-o_minus_bg }                  
    for k,v in readict.items():
        #if k!='era5fb':
            #continue
        for ftype in v['ftype']:  
            if 'fc' in ftype:
                dtype='fc'
            else:
                dtype='an'

            if dtype in readict[k].keys():

                idx=np.where(np.logical_and(obstype==85,z==20000))[0]
                rms=[]
                years=[]
                for iy in range(120):
                    #print(iy)

                    idy=np.where(np.logical_and(np.floor(ts[idx]/365.25/86400)==iy,
                                                np.abs(obs[idx]-readict[k][dtype]['refvalues'][idx])<12))[0]
                    if len(idy)>10:
                        rms.append(np.sqrt(np.nanmean((obs[idx[idy]]-readict[k][dtype]['refvalues'][idx[idy]])**2)))
                    else:
                        rms.append(np.nan)
                    years.append(iy)

                plt.plot(1900+np.array(years),np.array(rms),
                                 label='obs -'+k+'_'+dtype+', rms= {:5.3f}'.format(np.sqrt(np.nanmean(np.array(rms)**2))))

    plt.title('Monthly rms reanalysis t departures, 200 hPa, '+fn.split('/')[-1].split('_')[0])
    plt.ylabel('rms / K')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.expanduser('~/tmp/'+fnp+'rmsstats.png'))
    plt.close()

    #P.close()
    #P.join()
    #del P

    return readict


    #fpattern=path_to_gridded+'/era5fct.{}{:0>2}.130.nc'
    #func=partial(offline_fb,fpattern,lat,lon)
    #tfgs=list(map(func,yms))

def convert_missing(wpath,fn):
    tt=time.time()
    nanlist = [float('nan'), np.nan, 0, -2147483648]
    
    
    rscratch='/mnt/users/scratch/leo/scratch/'
    try:
        
        with open(os.path.expandvars(wpath+'/rea/'+fn.split('/')[-1].split('_CEUAS_merged_v1.nc')[0]+'pkl'),'rb') as f:
            readict=pickle.load(f)
    except:
        
        out_name = wpath+fn.split('/')[-1]  
        
        path_to_gridded=os.path.expandvars(rscratch+'/era5/gridded/')
        readict=retrieve_anfg(fn,out_name,path_to_gridded)
        with open(os.path.expandvars(wpath+'/rea/'+fn.split('/')[-1].split('_CEUAS_merged_v1.nc')[0]+'pkl'),'wb') as f:
            pickle.dump(readict,f)
    print (time.time()-tt)
    with eua.CDMDataset(fn) as data:
        keys = data.observations_table.keys()
        keys = [x for x in keys if not x.startswith('string')]
        #keys = [x for x in keys if x in ['conversion_flag','conversion_method','report_id','date_time','observation_value','observed_variable','observation_id','z_coordinate','z_coordinate_type']]
        if 'index' in keys: 
            keys.remove('index')
        if 'shape' in keys:
            keys.remove('shape')
        obskeys = keys

        ofb=True
        try:           
            keys = data.era5fb.keys()
            #keys = [x for x in keys if x in ['fg_depar@body','an_depar@body','biascorr@body','biascorr_fg@body']]
            keys = [x for x in keys if not x.startswith('string')]
            keys.remove('index')
            if 'shape' in keys:
                keys.remove('shape')
            fbkeys = keys
        except Exception as e:
            print(e)
            ofb=False


        # loading data:
        loaded_data=[]
        a_loaded_data=[]
        loaded_type = {'names':[],'formats':[]}
        ld=[]
        addmem=1000000
        for o in obskeys:
            if o in ['observed_variable','observation_value','observation_id','z_coordinate','z_coordinate_type','date_time','conversion_flag','conversion_method']:  
                if len(data.observations_table[o].shape)==1:
                    loaded_data.append((data.observations_table[o][:]))
                    a_loaded_data.append(numpy.empty_like(loaded_data[-1],shape=2*len(loaded_data[-1])+addmem))
                else:
                    loaded_data.append(data.observations_table[o][:].view('S{}'.format(data.observations_table[o].shape[1])).flatten())   
#                     a_loaded_data.append(numpy.empty((2*len(loaded_data[-1]),len(loaded_data[-1][0])), dtype=loaded_data[-1][0].dtype))
                    a_loaded_data.append(numpy.empty_like(loaded_data[-1],shape=2*len(loaded_data[-1])+addmem))
                    a_loaded_data[-1].fill(b' '*data.observations_table[o].shape[1])
                loaded_type['names'].append(o)
                loaded_type['formats'].append(loaded_data[-1].dtype)
                ld.append((o,loaded_data[-1].dtype))
        loaded_obstab = numpy.rec.fromarrays(loaded_data, dtype=ld)
        del loaded_data
        a_loaded_obstab = numpy.rec.fromarrays(a_loaded_data, dtype=ld)
        del a_loaded_data

        if ofb:
            loaded_fb=[]
            a_loaded_fb=[]
            loaded_type = {'names':[],'formats':[]}
            lf=[]
            for o in fbkeys:
                if o in ['fg_depar@body','an_depar@body','biascorr@body','biascorr_fg@body']:  
                    loaded_fb.append((data.era5fb[o][:]))
                    a_loaded_fb.append(numpy.empty_like(loaded_fb[-1],shape=2*len(loaded_fb[-1])+addmem))
                    loaded_type['names'].append(o)
                    loaded_type['formats'].append(loaded_fb[-1].dtype)
                    lf.append((o,loaded_fb[-1].dtype))
            loaded_feedback = numpy.rec.fromarrays(loaded_fb, dtype=lf)
            del loaded_fb
            a_loaded_feedback = numpy.rec.fromarrays(a_loaded_fb, dtype=lf)
            del a_loaded_fb

        @njit
        def add_fb(loaded_obstab,loaded_feedback,ref20CR,refera5an,refera5fc):
            i20=0
            iera=0
            for i in range(loaded_obstab['date_time'].shape[0]):
                if loaded_feedback['fg_depar@body'][i]!=loaded_feedback['fg_depar@body'][i]:
                    if loaded_obstab['observation_value'][i]==loaded_obstab['observation_value'][i]:
                        if refera5fc[i]==refera5fc[i]:
                            loaded_feedback['fg_depar@body'][i]=loaded_obstab['observation_value'][i]-refera5fc[i]
                            loaded_feedback['an_depar@body'][i]=loaded_obstab['observation_value'][i]-refera5an[i]
                            loaded_feedback['biascorr@body'][i]=0.
                            iera+=1
                        elif ref20CR[i]==ref20CR[i]:
                            loaded_feedback['fg_depar@body'][i]=loaded_obstab['observation_value'][i]-ref20CR[i]
                            loaded_feedback['an_depar@body'][i]=loaded_obstab['observation_value'][i]-ref20CR[i]
                            loaded_feedback['biascorr@body'][i]=0.
                            i20+=1
                if i%1000000==0:
                    print(i,i20,iera)
        if(ofb):      
            add_fb(loaded_obstab,loaded_feedback,readict['20CRv3']['an']['refvalues'],
                   readict['era5']['an']['refvalues'],readict['era5']['fc']['refvalues'])
        
        del readict    
        recordindex = data.recordindex[:]

        # --->

    print(time.time()-tt)
    idx,press,temp,relhum,spechum,dpd,dewpoint,uwind,vwind,wd,ws,d_temp,d_relhum,d_spechum,d_dpd,d_dewpoint,d_uwind,d_vwind,d_wd,d_ws,fgd_temp,fgd_relhum,fgd_spechum,fgd_dpd,fgd_dewpoint,fgd_uwind,fgd_vwind,fgd_wd,fgd_ws=ipl2(loaded_obstab, loaded_feedback)

    xtemp=xr.DataArray(temp)
    xpress=xr.DataArray(press)
    xrelhum=xr.DataArray(relhum)
    xspechum=xr.DataArray(spechum)
    xdpd=xr.DataArray(dpd)
    xdewpoint=xr.DataArray(dewpoint)
    cdpddp=temp-dewpoint
    cdpdrh=rasotools.met.convert.to_dpd(temp=xtemp,press=xpress,rel_humi=xrelhum).values
    cdpdsh=rasotools.met.convert.to_dpd(temp=xtemp,press=xpress,spec_humi=xspechum).values
    cshrh = rasotools.met.convert.to_sh(temp=xtemp, press=xpress, rel_humi=xrelhum).values
    cshdpd = rasotools.met.convert.to_sh(dpd=xtemp-xdewpoint, press=xpress, temp=xtemp).values
    crhsh = rasotools.met.convert.to_rh(temp=xtemp, spec_humi=xspechum, press=xpress).values
    crhdpd = rasotools.met.convert.to_rh(temp=xtemp,dpd=xtemp-xdewpoint).values

    d_xtemp=xr.DataArray(d_temp)
    d_xrelhum=xr.DataArray(d_relhum)
    d_xspechum=xr.DataArray(d_spechum)
    d_xdpd=xr.DataArray(d_dpd)
    d_xdewpoint=xr.DataArray(d_dewpoint)
    d_cdpddp=d_temp-d_dewpoint
    d_cdpdrh=rasotools.met.convert.to_dpd(temp=d_xtemp,press=xpress,rel_humi=d_xrelhum).values
    d_cdpdsh=rasotools.met.convert.to_dpd(temp=d_xtemp,press=xpress,spec_humi=d_xspechum).values
    d_cshrh = rasotools.met.convert.to_sh(temp=d_xtemp, press=xpress, rel_humi=d_xrelhum).values
    d_cshdpd = rasotools.met.convert.to_sh(dpd=d_xtemp-d_xdewpoint, press=xpress, temp=d_xtemp).values
    d_crhsh = rasotools.met.convert.to_rh(temp=d_xtemp, spec_humi=d_xspechum, press=xpress).values
    d_crhdpd = rasotools.met.convert.to_rh(temp=d_xtemp,dpd=d_xtemp-d_xdewpoint).values

    fgd_xtemp=xr.DataArray(fgd_temp)
    fgd_xrelhum=xr.DataArray(fgd_relhum)
    fgd_xspechum=xr.DataArray(fgd_spechum)
    fgd_xdpd=xr.DataArray(fgd_dpd)
    fgd_xdewpoint=xr.DataArray(fgd_dewpoint)
    fgd_cdpddp=fgd_temp-fgd_dewpoint
    fgd_cdpdrh=rasotools.met.convert.to_dpd(temp=fgd_xtemp,press=xpress,rel_humi=fgd_xrelhum).values
    fgd_cdpdsh=rasotools.met.convert.to_dpd(temp=fgd_xtemp,press=xpress,spec_humi=fgd_xspechum).values
    fgd_cshrh = rasotools.met.convert.to_sh(temp=fgd_xtemp, press=xpress, rel_humi=fgd_xrelhum).values
    fgd_cshdpd = rasotools.met.convert.to_sh(dpd=fgd_xtemp-fgd_xdewpoint, press=xpress, temp=fgd_xtemp).values
    fgd_crhsh = rasotools.met.convert.to_rh(temp=fgd_xtemp, spec_humi=fgd_xspechum, press=xpress).values
    fgd_crhdpd = rasotools.met.convert.to_rh(temp=fgd_xtemp,dpd=fgd_xtemp-fgd_xdewpoint).values

    idy=numpy.where(loaded_obstab['z_coordinate_type'][idx]==2) # do not convert humidity if data are not on pressure coordinates
    for c in cdpdrh,cshrh,cshdpd,crhdpd:
        c[idy]=numpy.nan

    cuwind = ws * np.cos(np.radians(270.-wd))
    cvwind = ws * np.sin(np.radians(270.-wd))
    cws = np.sqrt(uwind ** 2 + vwind ** 2)
    d_ws = np.sqrt(d_uwind ** 2 + d_vwind ** 2)
    fgd_ws = np.sqrt(fgd_uwind ** 2 + fgd_vwind ** 2)
    cwd = 90 - np.arctan2(-vwind, -uwind) * 180 / np.pi - 180.
    cwd = np.where(cwd > 0., cwd, 360.+cwd)
    d_cwd = 90 - np.arctan2(-d_vwind, -d_uwind) * 180 / np.pi - 180.
    d_cwd = np.where(cwd > 0., cwd, 360.+cwd)
    fgd_cwd = 90 - np.arctan2(-fgd_vwind, -fgd_uwind) * 180 / np.pi - 180.
    fgd_cwd = np.where(cwd > 0., cwd, 360.+cwd)

    humvar=numpy.array((ipar[34],ipar[36],ipar[38],ipar[39])) #dpd,dp,rh,sh
    wvar=numpy.array((ipar[104],ipar[105],ipar[106],ipar[107])) #dpd,dp,rh,sh

    reduced_obskeys=List(loaded_obstab.dtype.fields.keys())
    reduced_fbkeys=List(loaded_feedback.dtype.fields.keys())
    out, fb_out, ri, rt, jj, addedvar=augment(loaded_obstab, a_loaded_obstab, loaded_feedback, a_loaded_feedback,
                                             idx,temp,press,relhum,spechum,dpd,dewpoint,uwind,vwind,wd,ws,
                                             cdpddp,cdpdrh,cshrh,cshdpd,crhdpd,crhsh,cdpdsh,cuwind,cvwind,cwd,cws,
                                             d_cdpddp,d_cdpdrh,d_cdpdsh,d_cshrh,d_cshdpd,d_crhsh,d_crhdpd,
                                             fgd_cdpddp,fgd_cdpdrh,fgd_cdpdsh,fgd_cshrh,fgd_cshdpd,fgd_crhsh,fgd_crhdpd,
                                             d_cwd,d_ws,fgd_cwd,fgd_ws,
                                             humvar,wvar)
    
    for ii in range(104,108):
        #idx=np.where(loaded_obstab['observed_variable']==ii)[0]
        idy=np.where(a_loaded_obstab['observed_variable'][:jj]==ipar[ii])[0]    
        print('wind check',ipar[ii],len(idy))
    
    del temp,press,relhum,spechum,dpd,dewpoint,uwind,vwind,wd,ws,\
                                             cdpddp,cdpdrh,cshrh,cshdpd,crhdpd,crhsh,cdpdsh,cuwind,cvwind,cwd,cws,\
                                             d_cdpddp,d_cdpdrh,d_cdpdsh,d_cshrh,d_cshdpd,d_crhsh,d_crhdpd,\
                                             fgd_cdpddp,fgd_cdpdrh,fgd_cdpdsh,fgd_cshrh,fgd_cshdpd,fgd_crhsh,fgd_crhdpd,\
                                             d_cwd,d_ws,fgd_cwd,fgd_ws,xtemp,xpress,xrelhum,xspechum,xdpd,xdewpoint,\
                                             d_xtemp,d_xrelhum,d_xspechum,d_xdpd,d_xdewpoint,\
                                             fgd_xtemp,fgd_xrelhum,fgd_xspechum,fgd_xdpd,fgd_xdewpoint

    #avars = {}
    #fb_avars = {}
    #for i in reduced_obskeys:
        #avars[i] = out[i][:jj]  

    #for i in reduced_fbkeys:
        #fb_avars[i] = fb_out[i][:jj]

    print(time.time()-tt)    

    # sorting:
    print('start sorting')
    targetfile = rscratch+'converted_v8/'+fn.split('/')[-1] # wpath+fn.split('/')[-1]
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

#    allvars = copy.copy(avars['observed_variable'])
#    allvars.sort()
    obsv = out['observed_variable'][:jj]
    allvars = numpy.sort(numpy.unique(obsv))
    #
    #
    # resorting the data
    #
    @njit
    def make_vrindex(vridx,ridx): # this function is similar to np.unique with return_index=True, but it expands the index to the  
        # original array dimensions
        l=0
        for i in range(len(ridx)): # to set the recordindices
            if i == 0:
                l +=1
            else:
                if ridx[i]>ridx[i-1]:
                    vridx[ridx[i-1] + 1 : ridx[i] + 1]=l # next record after l
                    l += 1
                else:
                    l += 1
        vridx[ridx[i] + 1 :]=l #len(idx) # next record for the last element is the len of the data

    @njit
    def make_vrindexold(vridx,ridx,idx):
        l=0
        for i in range(1,len(idx)): # to set the recordindices
            if ridx[i]>ridx[i-1]:
                vridx[ridx[i-1]:ridx[i]]=l # next record after l
                l=i
        vridx[ridx[i]:]=l #len(idx) # next record for the last element is the len of the data

    @njit
    def make_vrindexn(vridx,ridx,idx,dta):
        #lold=0
        l=0
        for i in range(1,len(idx)): # to set the recordindices
            if dta[idx[i]]>dta[idx[i-1]]:
                vridx[ridx[l]:ridx[i]]=l # next record after l
                l=i
        vridx[ridx[i]:]=l
#        vridx[-1]=len(idx)
#                l+=1
#        vridx[ridx[i]:]=len(idx) # next record for the last element is the len of the data

    @njit
    def make_vrindex2(vridx,ridx,idx,dta):
        l=0
        for i in range(1,len(idx)): # to set the recordindices
            if ridx[i]>ridx[i-1] or dta[idx[i]]>dta[idx[i-1]]:
                vridx[ridx[i-1]:ridx[i]]=l # next record after l
                l=i
        vridx[ridx[i]:]=len(idx) # next record for the last element is the len of the data

    ridxall=np.zeros(obsv.shape[0],dtype=np.int64) # reverse index - index of the record index
    j=-1
    for j in range(len(ri)-1):
        ridxall[ri[j]:ri[j+1]]=j
    j+=1
    ridxall[ri[j]:]=j # for the last elemenet
    #dta=a_loaded_obstab['date_time'][:jj]
    ridx=[]
    vridx=[]
    absidx=[]
    absidx=np.empty_like(obsv)
    abscount=0
    idx = []
    for j in range(len(allvars)):
        idx.append(np.where(obsv==allvars[j])[0]) # index of all elements form certain variable j
#         print(j,len(idx),',',end='')
        vridx.append(np.zeros(ri.shape[0]+1,dtype=np.int64)) # all zeros in lenght of record index
        ridx=ridxall[idx[-1]] # ridxall where variable is j
        make_vrindex(vridx[-1],ridx)
        ##begin debugcopy
        #l=0
        ##print(ridx[0])
        
        ##vridx[-1][:ridx[0]]=l
        ##l += 1
        ##vridx[-1][ridx[0]] = l
        #for i in range(len(ridx)): # to set the recordindices
            ##print(ridx[i])
            #if i == 0:
                #l +=1
            #else:
                #if ridx[i]>ridx[i-1]:
                    #vridx[-1][ridx[i-1] + 1:ridx[i] + 1]=l # next record after l
                    #l += 1
                #else:
                    #l += 1
            ##print(i, vridx[-1][:20])
        #vridx[-1][ridx[i] + 1:]=l #len(idx) # next record for the last element is the len of the data
        ##end debugcopy
        vridx[-1]+=abscount # abscount for stacking the recordindex

        absidx[abscount:abscount+len(idx[-1])]=idx[-1] # why copy? - to make sure it's not just the ref. - maybe ok without the cp
        abscount+=len(idx[-1])
        vridx[-1][-1]=abscount
    #absidx=np.concatenate(absidx)

    # recordtimestamps are only necessary once
    recordtimestamps = rt
    
    # check integrity of indices
    ref = datetime(1900, 1, 1)
    od = out['date_time'][:jj][absidx]
    oov = out['observed_variable'][:jj][absidx]
    ooval = out['observation_value'][:jj][absidx]
    pres = out['z_coordinate'][:jj][absidx]
    for j in range(len(allvars)):
        for i in range(1, rt.shape[0]):
            
            if any(od[vridx[j][i]:vridx[j][i + 1]]!=rt[i]):
                print('spurious',allvars[j], i, ref+timedelta(seconds=int(od[vridx[j][i - 1]])),
                      ref+timedelta(seconds=int(od[vridx[j][i + 1] - 1])), ref+timedelta(seconds=int(rt[ridxall[i]])) )
                print('spurious',allvars[j], i, oov[vridx[j][i]],oov[vridx[j][i + 1] - 1])
                print('spurious',allvars[j], i, ooval[vridx[j][i]:vridx[j][i + 1]])
            #else:
                #if vridx[j][i] == vridx[j][i + 1]:
                    ##print('no values', allvars[j], i)
                    #pass
                #else:
                    #print('ok',allvars[j], i, ref+timedelta(seconds=int(od[vridx[j][i]])),
                      #ref+timedelta(seconds=int(od[vridx[j][i + 1] - 1])), ref+timedelta(seconds=int(rt[i])) )
                    #print('ok', allvars[j], i, oov[vridx[j][i]],oov[vridx[j][i + 1] - 1])
                    #print('ok',allvars[j], i, ooval[vridx[j][i]:vridx[j][i + 1]])
                
    del ridxall

    print('elapsed converting: ',time.time()-tt)
    
    #out_name = fn  
    #path_to_gridded=os.path.expandvars('$RSCRATCH/era5/gridded/')
    #readict=retrieve_anfg(out_name,path_to_gridded)

    tt=time.time()
    if os.path.isfile(targetfile):
        mode='r+'
    else:
        mode='w'

    for i in obskeys:
        print(i)
        print(time.time()-tt)

        if i == 'observation_id':
            ov_vars = numpy.empty(jj,dtype='S11')
            ov_vars[:]=out[i][:jj]
            ov_vars=fill_obsid(ov_vars.view('S1').reshape((len(ov_vars),11)),out['conversion_flag'][:jj])

        elif i in reduced_obskeys:
            ov_vars = out[i][:jj]
        else:
            continue

        ov_vars = ov_vars[absidx]
        if i == 'index':
            pass
        elif i in ['observation_id', 'report_id', 'sensor_id', 'source_id']:
            alldict = {i:np.asarray(ov_vars, dtype='S1')}
            write_dict_h5(targetfile, alldict, 'observations_table', {i: { 'compression': 'gzip' } }, [i])
        else:
            alldict = pandas.DataFrame({i:ov_vars})
            write_dict_h5(targetfile, alldict, 'observations_table', {i: { 'compression': 'gzip' } }, [i])  
    del ov_vars
    del obsv
    #del dta
    del out
    del a_loaded_obstab
    gc.collect()
    
    if True:
        for i in fbkeys:
            print(i)
            print(time.time()-tt)
    
            if i in reduced_fbkeys:
                ov_vars = fb_out[i][:jj]
            else: 
                continue
    
            ov_vars = ov_vars[absidx]
            
            if i == 'index':
                pass
            elif i in ['expver', 'source@hdr', 'source_id', 'statid@hdr']:
                alldict = {i:np.asarray(ov_vars, dtype='S1')}
                write_dict_h5(targetfile, alldict, 'era5fb', {i: { 'compression': 'gzip' } }, [i])
            else:
                alldict = pandas.DataFrame({i:ov_vars})
                write_dict_h5(targetfile, alldict, 'era5fb', {i: { 'compression': 'gzip' } }, [i]) 
    
        del fb_out
        del a_loaded_feedback
        gc.collect()
        
        if True:
            for i in obskeys:
                if i in reduced_obskeys or i == 'observation_id' :
                    continue
                #if i not in ('latitude','longitude','report_id'):
                    #continue
                    
                print(i)
                print(time.time()-tt)
                with eua.CDMDataset(fn) as data:
                    #i='z_coordinate'
                    rest_data = data.observations_table[i][:]
                if rest_data.ndim==2: #i in ['observation_id', 'report_id', 'sensor_id', 'source_id']:
                    final = numpy.empty((addedvar[-1][1],len(rest_data[0])), dtype=rest_data[0].dtype)
                else:
                    final = numpy.empty(addedvar[-1][1], dtype=rest_data[0].dtype)
                ov_vars = fill_restdata(final, rest_data, addedvar, jj) #,out['z_coordinate'][:jj])
        
                ov_vars = ov_vars[absidx]
                if i == 'index':
                    pass
                elif i in ['observation_id', 'report_id', 'sensor_id', 'source_id']:
                    alldict = {i:np.asarray(ov_vars, dtype='S1')}
                    write_dict_h5(targetfile, alldict, 'observations_table', {i: { 'compression': 'gzip' } }, [i])
                else:
                    alldict = pandas.DataFrame({i:ov_vars})
                    write_dict_h5(targetfile, alldict, 'observations_table', {i: { 'compression': 'gzip' } }, [i])  
        
        if True:
            for i in fbkeys:
                print(i)
                print(time.time()-tt)
        
                if i in reduced_fbkeys:
                    continue
                else: 
                    with eua.CDMDataset(fn) as data:
                        rest_data = data.era5fb[i][:]
                    if i in ['expver', 'source@hdr', 'source_id', 'statid@hdr']:
                        final = numpy.empty((addedvar[-1][1],len(rest_data[0])), dtype=rest_data[0].dtype)
                    else:
                        final = numpy.empty(addedvar[-1][1], dtype=rest_data[0].dtype)
                    ov_vars = fill_restdata(final, rest_data, addedvar, jj)
        
                ov_vars = ov_vars[absidx]
                
                if i == 'index':
                    pass
                elif i in ['expver', 'source@hdr', 'source_id', 'statid@hdr']:
                    alldict = {i:np.asarray(ov_vars, dtype='S1')}
                    write_dict_h5(targetfile, alldict, 'era5fb', {i: { 'compression': 'gzip' } }, [i])
                else:
                    alldict = pandas.DataFrame({i:ov_vars})
                    write_dict_h5(targetfile, alldict, 'era5fb', {i: { 'compression': 'gzip' } }, [i]) 
        #
        # writing the recordindices and recordtimestamp.
        #       
    recordindices=vridx
    for i in range(len(recordindices)):
        testvar = pandas.DataFrame({str(allvars[i]):recordindices[i]})
        write_dict_h5(targetfile, testvar, 'recordindices', {str(allvars[i]): { 'compression': None } }, [str(allvars[i])]) 

    write_dict_h5(targetfile, {'recordtimestamp':recordtimestamps}, 'recordindices', {'recordtimestamp': { 'compression': None } }, ['recordtimestamp'])

    print('elapsed writing '+targetfile+':',time.time()-tt)
    f= open(wlpath+fn.split('/')[-1]+".txt","w+")
    f.write("done") 
    f.close()
    return
    

# files = glob.glob('/raid60/scratch/federico/MERGED_DATABASE_OCTOBER2020_sensor/0-20000-0-01*.nc')
# files = glob.glob('/raid60/scratch/federico/DATABASE_JANUARY2021_sensor/0-20000-0-97690*.nc')
# files = glob.glob('/raid60/scratch/federico/DATABASE_JANUARY2021_sensor/0-20500-0-93954*.nc')
# files = glob.glob('/raid60/scratch/federico/DATABASE_JANUARY2021_sensor/0-20400-0-04665*.nc')
# files = glob.glob('/raid60/scratch/federico/DATABASE_JANUARY2021_sensor/*.nc')

# print(files[:10])

# convert_missing(files[6020])
# convert_missing('/raid60/scratch/federico/MERGED_DATABASE_OCTOBER2020_sensor/0-20000-0-03414_CEUAS_merged_v0.nc')

if __name__ == '__main__':
    
    no_height = ['/raid60/scratch/federico/DATABASE_JANUARY2021_FIXED_sensor/0-20000-0-41915_CEUAS_merged_v0.nc',
                 '/raid60/scratch/federico/DATABASE_JANUARY2021_FIXED_sensor/0-20000-0-94231_CEUAS_merged_v0.nc',
                 '/raid60/scratch/federico/DATABASE_JANUARY2021_FIXED_sensor/0-20000-0-43009_CEUAS_merged_v0.nc',
                 '/raid60/scratch/federico/DATABASE_JANUARY2021_FIXED_sensor/0-20000-0-40951_CEUAS_merged_v0.nc',
                 '/raid60/scratch/federico/DATABASE_JANUARY2021_FIXED_sensor/0-20000-0-43042_CEUAS_merged_v0.nc',
                 '/raid60/scratch/federico/DATABASE_JANUARY2021_FIXED_sensor/0-20000-0-94323_CEUAS_merged_v0.nc',
                 '/raid60/scratch/federico/DATABASE_JANUARY2021_FIXED_sensor/0-20000-0-41565_CEUAS_merged_v0.nc',
                 '/raid60/scratch/federico/DATABASE_JANUARY2021_FIXED_sensor/0-20000-0-41738_CEUAS_merged_v0.nc',
                 '/raid60/scratch/federico/DATABASE_JANUARY2021_FIXED_sensor/0-20000-0-42484_CEUAS_merged_v0.nc',
                 '/raid60/scratch/federico/DATABASE_JANUARY2021_FIXED_sensor/0-20000-0-42354_CEUAS_merged_v0.nc',
                 '/raid60/scratch/federico/DATABASE_JANUARY2021_FIXED_sensor/0-20000-0-43259_CEUAS_merged_v0.nc',
                 '/raid60/scratch/federico/DATABASE_JANUARY2021_FIXED_sensor/0-20000-0-62381_CEUAS_merged_v0.nc',
                 '/raid60/scratch/federico/DATABASE_JANUARY2021_FIXED_sensor/0-20500-0-80417_CEUAS_merged_v0.nc',
                 '/raid60/scratch/federico/DATABASE_JANUARY2021_FIXED_sensor/0-20000-0-95721_CEUAS_merged_v0.nc',
                 '/raid60/scratch/federico/DATABASE_JANUARY2021_FIXED_sensor/0-20000-0-78088_CEUAS_merged_v0.nc',
                 '/raid60/scratch/federico/DATABASE_JANUARY2021_FIXED_sensor/0-20000-0-87544_CEUAS_merged_v0.nc',
                 '/raid60/scratch/federico/DATABASE_JANUARY2021_FIXED_sensor/0-20000-0-38696_CEUAS_merged_v0.nc',
                 '/raid60/scratch/federico/DATABASE_JANUARY2021_FIXED_sensor/0-20000-0-04085_CEUAS_merged_v0.nc',]    

    wpath= os.path.expandvars('$RSCRATCH/converted_v8/') #'./'
    opath=wpath
    wlpath=wpath+'log/'
    for p in wpath,wlpath:      
        try:
            os.mkdir(p)
        except:
            pass
    
    #for i in no_height:
        #do_resort(i)
        #f= open(wlpath+i.split('/')[-1]+".txt","w+")
        #f.write("done") 
        #f.close()

    files = glob.glob('/mnt/scratch/scratch/federico/MERGED_15JUNE2022/*64600*v1.nc')
    files_to_convert = glob.glob('/raid60/scratch/federico/MERGED_JUNE2021/*72357*v1.nc')
    #files = glob.glob('/scratch/das/federico/TRY_MERGED_FEB2022/*89564*v1.nc')
    files_to_convert = files #glob.glob('/scratch/das/federico/TRY_MERGED_JAN2022/*72357*v1.nc')
    
    ipar=np.zeros(140,dtype=np.int32)-2100000000 # 
    ipar[0]=0
    ipar[34]=34
    ipar[39]=39
    ipar[85]=126
    ipar[106]=106
    ipar[107]=107
    ipar[117]=117
    #ipar[]=136
    ipar[36]=137 #dp
    ipar[38]=138 #rh
    ipar[104]=139
    ipar[105]=140

    print(files)
    #nfiles=glob.glob(wpath+'*v1.nc')
    #nfile=os.path.dirname(nfiles[0])+'/'+files[0].split('/')[-1]
    #ofiles=['/mnt/users/scratch/leo/scratch/converted_v7/0-20000-0-89564_CEUAS_merged_v1.nc']
    #nfiles=['/mnt/users/scratch/leo/scratch/converted_v8/0-20000-0-89564_CEUAS_merged_v1.nc']
    #with h5py.File(ofiles[0],'r') as f:
        #vcode=f['observations_table']['observed_variable'][:]
        #vcodes=np.unique(vcode)
        #print(vcodes)
        #print(ofiles[0])
        #for v in vcodes:
            #print(v,np.sum(v==vcode))
        #with h5py.File(nfiles[0],'r') as g:
            #vcode=g['observations_table']['observed_variable'][:]
            #vcodes=np.unique(vcode)
            #print(vcodes)
            #print(nfiles[0])
            #for v in vcodes:
                #print(v,np.sum(v==vcode))
            
            #nstart=g['recordindices']['126'][0]
            #nstop=g['recordindices']['126'][-1]
            #plt.plot(g['observations_table']['date_time'][nstart:nstop]/86400/365.25,g['observations_table']['observation_value'][nstart:nstop])
            #ostart=f['recordindices']['85'][0]
            #ostop=f['recordindices']['85'][-1]
            #plt.plot(f['observations_table']['date_time'][ostart:ostop]/86400/365.25,f['observations_table']['observation_value'][ostart:ostop])
            #print(nstop-nstart,ostop-ostart)
    
    already_done = glob.glob(wlpath+'*.txt')

#     files_to_convert = []
#     for i in files:
#         if not wlpath+i.split('/')[-1]+'.txt' in already_done:
#             files_to_convert.append(i)
#     files_to_convert.sort()
    tt=time.time()
    
    #for i in files_to_convert:
        #print(i)
        #convert_missing(wpath, i)
        

    pool = multiprocessing.Pool(processes=20)
#     #result_list = pool.map(convert_missing, files_to_convert)
    func=partial(convert_missing,wpath)
    result_list = list(map(func, files_to_convert))
#     print(result_list)
    print('total:',time.time()-tt)

