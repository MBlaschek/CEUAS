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
import gc
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

#opath='/raid60/scratch/uli/converted_v5/'
opath='/raid60/scratch/leo/scratch/converted_v5/'
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
    
        if observed_variable[i]==38:
            relhum[j]=observation_value[i]
            d_relhum[j]=observation_value[i]-departure[i]
            fgd_relhum[j]=observation_value[i]-fg_departure[i]
        elif observed_variable[i]==39:
            spechum[j]=observation_value[i]
            d_spechum[j]=observation_value[i]-departure[i]
            fgd_spechum[j]=observation_value[i]-fg_departure[i]
        elif observed_variable[i]==34:
            dpd[j]=observation_value[i]
            d_dpd[j]=observation_value[i]-departure[i]
            fgd_dpd[j]=observation_value[i]-fg_departure[i]
        elif observed_variable[i]==36:
            dewpoint[j]=observation_value[i]
            d_dewpoint[j]=observation_value[i]-departure[i]
            fgd_dewpoint[j]=observation_value[i]-fg_departure[i]
        elif observed_variable[i]==85:
            temp[j]=observation_value[i]
            d_temp[j]=observation_value[i]-departure[i]
            fgd_temp[j]=observation_value[i]-fg_departure[i]
        elif observed_variable[i]==104:
            uwind[j]=observation_value[i]
            d_uwind[j]=observation_value[i]-departure[i]
            fgd_uwind[j]=observation_value[i]-fg_departure[i]
        elif observed_variable[i]==105:
            vwind[j]=observation_value[i]
            d_vwind[j]=observation_value[i]-departure[i]
            fgd_vwind[j]=observation_value[i]-fg_departure[i]
        elif observed_variable[i]==106:
            wd[j]=observation_value[i]
            d_wd[j]=observation_value[i]-departure[i]
            fgd_wd[j]=observation_value[i]-fg_departure[i]
        elif observed_variable[i]==107:
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
    if h==34:
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
            
    elif h==36:
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
        
    elif h==38:
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
            
    elif h==39:
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
    if h==104:
        if cuwind[k]==cuwind[k]:
            a_observation_value[j]=cuwind[k]
            a_an_depar[j]=numpy.nan
            a_fg_depar[j]=numpy.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=1
    elif h==105:
        if cvwind[k]==cvwind[k]:
            a_observation_value[j]=cvwind[k]
            a_an_depar[j]=numpy.nan
            a_fg_depar[j]=numpy.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=1
    elif h==106:
        if cwd[k]==cwd[k]:
            a_observation_value[j]=cwd[k]
            a_an_depar[j]=d_cwd[k]
            a_fg_depar[j]=fgd_cwd[k]
            if (cwd[k]<0.) or (cwd[k]>360.):
                a_observation_value[j]=numpy.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=2
    elif h==107:
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
        if len(wlist)>1:
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
            wlist.clear()
        if idxu!=obstab['observation_value'].shape[0]:        
            if obstab['date_time'][idxu] != obstab['date_time'][idx[k]]:
                recordindex[ri]=j+1
                recordtimestamp[ri]=obstab['date_time'][idxu]
                ri+=1

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

def convert_missing(fn, destination: str = opath):
    tt=time.time()
    nanlist = [float('nan'), np.nan, 0, -2147483648]
    with eua.CDMDataset(fn) as data:
        keys = data.observations_table.keys()
        keys = [x for x in keys if not x.startswith('string')]
        #keys = [x for x in keys if x in ['conversion_flag','conversion_method','report_id','date_time','observation_value','observed_variable','observation_id','z_coordinate','z_coordinate_type']]
        keys.remove('index')
        keys.remove('shape')
        obskeys = keys

        keys = data.era5fb.keys()
        #keys = [x for x in keys if x in ['fg_depar@body','an_depar@body','biascorr@body','biascorr_fg@body']]
        keys = [x for x in keys if not x.startswith('string')]
        keys.remove('index')
        keys.remove('shape')
        fbkeys = keys


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

    cuwind = ws * np.cos(np.radians(wd))
    cvwind = ws * np.sin(np.radians(wd))
    cws = np.sqrt(uwind ** 2 + vwind ** 2)
    d_ws = np.sqrt(d_uwind ** 2 + d_vwind ** 2)
    fgd_ws = np.sqrt(fgd_uwind ** 2 + fgd_vwind ** 2)
    cwd = 90 - np.arctan2(-vwind, -uwind) * 180 / np.pi - 180.
    cwd = np.where(cwd > 0., cwd, 360.+cwd)
    d_cwd = 90 - np.arctan2(-d_vwind, -d_uwind) * 180 / np.pi - 180.
    d_cwd = np.where(cwd > 0., cwd, 360.+cwd)
    fgd_cwd = 90 - np.arctan2(-fgd_vwind, -fgd_uwind) * 180 / np.pi - 180.
    fgd_cwd = np.where(cwd > 0., cwd, 360.+cwd)

    humvar=numpy.array((34,36,38,39)) #dpd,dp,rh,sh
    wvar=numpy.array((104,105,106,107)) #dpd,dp,rh,sh

    reduced_obskeys=List(loaded_obstab.dtype.fields.keys())
    reduced_fbkeys=List(loaded_feedback.dtype.fields.keys())
    out, fb_out, ri, rt, jj, addedvar=augment(loaded_obstab, a_loaded_obstab, loaded_feedback, a_loaded_feedback,
                                             idx,temp,press,relhum,spechum,dpd,dewpoint,uwind,vwind,wd,ws,
                                             cdpddp,cdpdrh,cshrh,cshdpd,crhdpd,crhsh,cdpdsh,cuwind,cvwind,cwd,cws,
                                             d_cdpddp,d_cdpdrh,d_cdpdsh,d_cshrh,d_cshdpd,d_crhsh,d_crhdpd,
                                             fgd_cdpddp,fgd_cdpdrh,fgd_cdpdsh,fgd_cshrh,fgd_cshdpd,fgd_crhsh,fgd_crhdpd,
                                             d_cwd,d_ws,fgd_cwd,fgd_ws,
                                             humvar,wvar)
    
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

#    allvars = copy.copy(avars['observed_variable'])
#    allvars.sort()
    obsv = out['observed_variable'][:jj]
    allvars = numpy.sort(numpy.unique(obsv))
    #
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
    for j in range(len(allvars)):
        idx=np.where(obsv==allvars[j])[0] # index of all elements form certain variable j
#         print(j,len(idx),',',end='')
        vridx.append(np.zeros(ri.shape[0]+1,dtype=np.int64)) # all zeros in lenght of record index
        ridx=ridxall[idx] # ridxall where variable is j
        make_vrindex(vridx[-1],ridx,idx)
        #make_vrindex(vridx[-1],ridx,idx)
        vridx[-1]+=abscount # abscount for stacking the recordindex

        absidx[abscount:abscount+len(idx)]=idx # why copy? - to make sure it's not just the ref. - maybe ok without the cp
        abscount+=len(idx)
        vridx[-1][-1]=abscount
    #absidx=np.concatenate(absidx)

    # recordtimestamps are only necessary once
    recordtimestamps = rt 
    del ridxall

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
        
        for i in obskeys:
            if i in reduced_obskeys or i == 'observation_id':
                continue
            print(i)
            print(time.time()-tt)
            with eua.CDMDataset(fn) as data:
                #i='z_coordinate'
                rest_data = data.observations_table[i][:]
            if i in ['observation_id', 'report_id', 'sensor_id', 'source_id']:
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

    wpath='/raid60/scratch/leo/scratch/converted_v5/'
    wlpath=wpath+'log/'
    try:
        os.mkdir(wpath)
        os.mkdir(wlpath)
    except:
        pass
    
    #for i in no_height:
        #do_resort(i)
        #f= open(wlpath+i.split('/')[-1]+".txt","w+")
        #f.write("done") 
        #f.close()

    files = glob.glob('/raid60/scratch/federico/DATABASE_JANUARY2021_FIXED_sensor/*6610*.nc')
    already_done = glob.glob(wlpath+'*.txt')

    files_to_convert = []
    for i in files:
        if not wlpath+i.split('/')[-1]+'.txt' in already_done:
            files_to_convert.append(i)
    
#     for i in files_to_convert:
#         print(i)
#         convert_missing(i)

    pool = multiprocessing.Pool(processes=10)
    #result_list = pool.map(convert_missing, files_to_convert)
    result_list = list(map(convert_missing, files_to_convert))
    print(result_list)

