#!/usr/bin/env
# coding: utf-8

from numba import njit
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
from multiprocessing import Pool
sys.path.insert(0,os.getcwd()+'/../resort/rasotools-master/')
import rasotools
dir(rasotools)
import warnings
from functools import partial
import matplotlib.pylab as plt
warnings.filterwarnings('ignore')
import pickle
import ray

#opath='/raid60/scratch/uli/converted_v5/'
#opath='/raid60/scratch/leo/scratch/converted_v5/'
# if there are nan values in the pressure level - we will just sort without any converting!
def rmeanw(t,runmean):
    tmean=t.copy()
    index=np.zeros(tmean.shape[0],dtype='int')
    
    tret=rmean(t,tmean,index,runmean)
    tret[:runmean]=np.nan
    tret[-runmean:]=np.nan
    return tret

@njit(cache=True)
def rmean(t,tmean,index,runmean):

    tret=np.zeros(t.shape[0])
    ni=t.shape[0]
    good=runmean-runmean
    if runmean<2:
        for i in range(ni):
            tret[i]=t[i]
    else:

        for j in range(ni):
            tret[j]=np.nan
            if t[j]==t[j]:
                index[good]=j
                good+=1

        if good>runmean+2:
            i=runmean//2
            tmean[:]=np.nan
            if runmean%2==1:
                tmean[i]=0.
                for k in range(-runmean//2+1,runmean//2+1):
                    tmean[i]+=t[index[i+k]]
                tmean[i]/=runmean

                for i in range(runmean//2+1,good-runmean//2):
                    tmean[i]=(tmean[i-1]*runmean+t[index[i+runmean//2]])/(runmean+1)

            else:

                i=runmean//2
                tmean[i]=0.
                for k in range(-runmean//2,runmean//2):
                    tmean[i]+=t[index[i+k]]
                tmean[i]/=runmean

                for i in range(runmean//2+1,good-runmean//2-1):
                    tmean[i]=(tmean[i-1]*runmean+t[index[i+runmean//2-1]])/(runmean+1)

            for i in range(good):
                tret[index[i]]=tmean[i]
        else:
            for i in range(good):
                tret[index[i]]=t[index[i]]

    return tret

@njit(cache=True)
def thin2(t,n):

    ni=t.shape[0]
    index=np.zeros(ni//n,dtype=np.int32)-1
    if n<2:
        for i in range(t.shape[0]):
            index[i]=i
    else:
        ni=t.shape[0]//n
        for i in range(ni):
            index[i]=i*n
            for j in range(n):
                idx=i*n+j
                if t[idx]==t[idx]:
                    index[i]=idx
                    break

    return index

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
    allvars = np.unique(allvars)
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
    jdx=np.where(np.logical_and(observed_variable==85,~np.isnan(observation_value)))[0]
    idx=np.empty(len(jdx),dtype=np.int64)
    press=np.empty(len(jdx),dtype=z_coordinate.dtype)
    relhum=np.empty(len(jdx),dtype=observation_value.dtype)
    dewpoint=np.empty(len(jdx),dtype=observation_value.dtype)
    dpd=np.empty(len(jdx),dtype=observation_value.dtype)
    spechum=np.empty(len(jdx),dtype=observation_value.dtype)

    uwind=np.empty(len(jdx),dtype=observation_value.dtype)
    vwind=np.empty(len(jdx),dtype=observation_value.dtype)
    ws=np.empty(len(jdx),dtype=observation_value.dtype)
    wd=np.empty(len(jdx),dtype=observation_value.dtype)

    for v in relhum,dewpoint,dpd,spechum:
        v.fill(np.nan)
    temp=np.empty(len(jdx),dtype=observation_value.dtype)
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


    ##jdx=len(np.where(np.logical_and(observed_variable==85,~np.isnan(observation_value)))[0]
    idx=np.empty(jdx,dtype=np.int64)
    press=np.empty(jdx,dtype=z_coordinate.dtype) # either height or pressure, depending on coordinate type
    relhum=np.empty(jdx,dtype=observation_value.dtype)
    dewpoint=np.empty(jdx,dtype=observation_value.dtype)
    dpd=np.empty(jdx,dtype=observation_value.dtype)
    spechum=np.empty(jdx,dtype=observation_value.dtype)
    uwind=np.empty(jdx,dtype=observation_value.dtype)
    vwind=np.empty(jdx,dtype=observation_value.dtype)
    ws=np.empty(jdx,dtype=observation_value.dtype)
    wd=np.empty(jdx,dtype=observation_value.dtype)
    temp=np.empty(jdx,dtype=observation_value.dtype)

    d_temp=np.empty(jdx,dtype=observation_value.dtype)
    d_relhum=np.empty(jdx,dtype=observation_value.dtype)
    d_spechum=np.empty(jdx,dtype=observation_value.dtype)
    d_dpd=np.empty(jdx,dtype=observation_value.dtype)
    d_dewpoint=np.empty(jdx,dtype=observation_value.dtype)
    d_uwind=np.empty(jdx,dtype=observation_value.dtype)
    d_vwind=np.empty(jdx,dtype=observation_value.dtype)
    d_wd=np.empty(jdx,dtype=observation_value.dtype)
    d_ws=np.empty(jdx,dtype=observation_value.dtype)

    fgd_temp=np.empty(jdx,dtype=observation_value.dtype)
    fgd_relhum=np.empty(jdx,dtype=observation_value.dtype)
    fgd_spechum=np.empty(jdx,dtype=observation_value.dtype)
    fgd_dpd=np.empty(jdx,dtype=observation_value.dtype)
    fgd_dewpoint=np.empty(jdx,dtype=observation_value.dtype)
    fgd_uwind=np.empty(jdx,dtype=observation_value.dtype)
    fgd_vwind=np.empty(jdx,dtype=observation_value.dtype)
    fgd_wd=np.empty(jdx,dtype=observation_value.dtype)
    fgd_ws=np.empty(jdx,dtype=observation_value.dtype)

    for v in temp,relhum,dewpoint,dpd,spechum,uwind,vwind,wd,ws,d_temp,d_relhum,d_spechum,d_dpd,d_dewpoint,d_uwind,d_vwind,d_wd,d_ws,fgd_temp,fgd_relhum,fgd_spechum,fgd_dpd,fgd_dewpoint,fgd_uwind,fgd_vwind,fgd_wd,fgd_ws:
        v.fill(np.nan)
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
            if (np.abs(cdpddp[k])>80) or (cdpddp[k]<0.01):
                a_observation_value[j]=np.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=2
        elif cdpdrh[k]==cdpdrh[k]:
            a_observation_value[j]=cdpdrh[k]
            a_an_depar[j]=cdpdrh[k]-d_cdpdrh[k]
            a_fg_depar[j]=cdpdrh[k]-fgd_cdpdrh[k]
            if (np.abs(cdpdrh[k])>80) or (cdpdrh[k]<0.01):
                a_observation_value[j]=np.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=3
        else:
            a_observation_value[j]=cdpdsh[k]
            a_an_depar[j]=cdpdsh[k]-d_cdpdsh[k]
            a_fg_depar[j]=cdpdsh[k]-fgd_cdpdsh[k]
            if (np.abs(cdpdsh[k])>80) or (cdpdsh[k]<0.01):
                a_observation_value[j]=np.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=4

    elif h==ipar[36]:
        if cdpdrh[k]==cdpdrh[k]:
            a_observation_value[j]=temp[k]-cdpdrh[k]
            a_an_depar[j]=(temp[k]-cdpdrh[k])-(temp[k]-d_cdpdrh[k])
            a_fg_depar[j]=(temp[k]-cdpdrh[k])-(temp[k]-fgd_cdpdrh[k])
            if (np.abs(cdpdrh[k])>80) or (cdpdrh[k]<0.01):
                a_observation_value[j]=np.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=3
        elif cdpdrh[k]==cdpdrh[k]:
            a_observation_value[j]=temp[k]-cdpddp[k]
            a_an_depar[j]=(temp[k]-cdpddp[k])-(temp[k]-d_cdpddp[k])
            a_fg_depar[j]=(temp[k]-cdpddp[k])-(temp[k]-fgd_cdpddp[k])
            if (np.abs(cdpddp[k])>80) or (cdpddp[k]<0.01):
                a_observation_value[j]=np.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=2
        else:
            a_observation_value[j]=temp[k]-cdpdsh[k]
            a_an_depar[j]=(temp[k]-cdpdsh[k])-(temp[k]-d_cdpdsh[k])
            a_fg_depar[j]=(temp[k]-cdpdsh[k])-(temp[k]-fgd_cdpdsh[k])
            if (np.abs(cdpdsh[k])>80) or (cdpdsh[k]<0.01):
                a_observation_value[j]=np.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=4

    elif h==ipar[38]:
        if crhdpd[k]==crhdpd[k]:
            a_observation_value[j]=crhdpd[k]
            a_an_depar[j]=crhdpd[k]-d_crhdpd[k]
            a_fg_depar[j]=crhdpd[k]-fgd_crhdpd[k]
            if (crhdpd[k]<0.) or (crhdpd[k]>1.03):
                a_observation_value[j]=np.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=2
        else: 
            a_observation_value[j]=crhsh[k]
            a_an_depar[j]=crhsh[k]-d_crhsh[k]
            a_fg_depar[j]=crhsh[k]-fgd_crhsh[k]
            if (crhsh[k]<0.) or (crhsh[k]>1.03):
                a_observation_value[j]=np.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=4

    elif h==ipar[39]:
        if cshrh[k]==cshrh[k]:
            a_observation_value[j]=cshrh[k]
            a_an_depar[j]=cshrh[k]-d_cshrh[k]
            a_fg_depar[j]=cshrh[k]-fgd_cshrh[k]
            if (cshrh[k]<0.) or (cshrh[k]>50.):
                a_observation_value[j]=np.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=3
        else:
            a_observation_value[j]=cshdpd[k]
            a_an_depar[j]=cshdpd[k]-d_cshdpd[k]
            a_fg_depar[j]=cshdpd[k]-fgd_cshdpd[k]
            if (cshdpd[k]<0.) or (cshdpd[k]>50.):
                a_observation_value[j]=np.nan
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
            a_an_depar[j]=np.nan
            a_fg_depar[j]=np.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=1
    elif h==ipar[105]:
        if cvwind[k]==cvwind[k]:
            a_observation_value[j]=cvwind[k]
            a_an_depar[j]=np.nan
            a_fg_depar[j]=np.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=1
    elif h==ipar[106]:
        if cwd[k]==cwd[k]:
            a_observation_value[j]=cwd[k]
            a_an_depar[j]=d_cwd[k]
            a_fg_depar[j]=fgd_cwd[k]
            if (cwd[k]<0.) or (cwd[k]>360.):
                a_observation_value[j]=np.nan
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

#@njit(cache=True, boundscheck=True)          
def augment(obstab, a_obstab, loaded_feedback, a_loaded_feedback,
            idx,temp,press,relhum,spechum,dpd,dewpoint,uwind,vwind,wd,ws,
            cdpddp,cdpdrh,cshrh,cshdpd,crhdpd,crhsh,cdpdsh,cuwind,cvwind,cwd,cws,
            d_cdpddp,d_cdpdrh,d_cdpdsh,d_cshrh,d_cshdpd,d_crhsh,d_crhdpd,
            fgd_cdpddp,fgd_cdpdrh,fgd_cdpdsh,fgd_cshrh,fgd_cshdpd,fgd_crhsh,fgd_crhdpd,
            d_cwd,d_cws,fgd_cwd,fgd_cws,
            humvar,wvar):


    recordindex=np.empty(idx.shape[0],obstab['date_time'].dtype)
    recordtimestamp=np.empty(idx.shape[0],obstab['date_time'].dtype)

    j=-1 # augmented index
    jsave=0
    p= obstab['z_coordinate'][0]-1 # z_coordinate[0]-1
    rts= obstab['date_time'][0] # date_time[0]

    addedvar=np.zeros((a_obstab['observation_value'].shape[0],2),dtype=np.int32)
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
                    a_loaded_feedback['biascorr@body'][j]=np.nan
                    a_loaded_feedback['biascorr_fg@body'][j]=np.nan
                    qconvert(j,k,obstab['observed_variable'][i],a_obstab['observation_value'],a_obstab['conversion_flag'],a_obstab['conversion_method'],
                             a_loaded_feedback['an_depar@body'],a_loaded_feedback['fg_depar@body'],
                             temp,cdpddp,cdpdrh,crhdpd,cshrh,cshdpd,crhsh, cdpdsh, 
                             d_cdpddp,d_cdpdrh,d_cdpdsh,d_cshrh,d_cshdpd,d_crhsh,d_crhdpd,
                             fgd_cdpddp,fgd_cdpdrh,fgd_cdpdsh,fgd_cshrh,fgd_cshdpd,fgd_crhsh,fgd_crhdpd)
                elif obstab['observed_variable'][i] in wvar:
                    a_loaded_feedback['biascorr@body'][j]=np.nan
                    a_loaded_feedback['biascorr_fg@body'][j]=np.nan
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
                    a_loaded_feedback['biascorr@body'][j]=np.nan
                    a_loaded_feedback['biascorr_fg@body'][j]=np.nan
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
                    a_loaded_feedback['biascorr@body'][j]=np.nan
                    a_loaded_feedback['biascorr_fg@body'][j]=np.nan
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
    #addv=np.empty((len(addedvar),2),dtype=np.int64)
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
        #x=np.nansum(final[addedvar[l,1]:addedvar[l+1,1]]-test[addedvar[l,1]:addedvar[l+1,1]])
        #if ~np.isnan(x) and x!=0:
                        #print(addedvar[l:l+2])                             
                        #print('x')

    return final

def split(x): 
    return [i.encode() for i in x.decode()]

@njit(cache=True)
def polinl(XA,YA,X):


#    REAL XA(4),YA(NFELD,4),X,Y(NFELD)
#    REAL XX1,XX2,XX3,XX4,X1X2

    Y = np.zeros(YA.shape[0])
    XX1=X-XA[0]
    XX2=X-XA[1]
    XX3=X-XA[2]
    XX4=X-XA[3]
    X1X2=XA[0]-XA[1]
    X1X2=X1X2*X1X2*X1X2

    for I in range(Y.shape[0]):
        Y[I]=XX2*XX3*XX4/3.0*YA[I,0]-XX1*XX3*XX4*YA[I,1]+XX1*XX2*XX4*YA[I,2]-XX1*XX2*XX3/3.0*YA[I,3]
        Y[I]=Y[I]/X1X2/2.0

    #for I in range(Y.shape[0]):
        #if Y[I]>YA[I,1] and Y[I]> YA[I,2]:
            #Y[I]=np.max((YA[I,1],YA[I,2]))
        #if Y[I] < YA[I,1] and Y[I] < YA[I,2]:
            #Y[I]=np.min((YA[I,1],YA[I,2]))

    return Y

@njit(cache=True)
def triint(tera5, reatab, tlon, tlat, lons, lats, secs, ts,obss,press, zs, p, dy):
    sh = tera5.shape
    pval = np.zeros((sh[2],sh[3]))
    yval = np.zeros((1, 4))
    
    idx = np.zeros(len(reatab), dtype=np.int32)
    oldlats = -100.
    oldlons = -400.
    oldts = -1
    im = 0
    for il in range(len(reatab)):
        if obss[il] == p:
            found = False
            for ip in range(sh[3]):
                if press[ip] == zs[il]:
                    found = True
                    break
            if not found:
                continue
            if ts[il] != oldts:
                
                its =np.searchsorted(secs, ts[il]) - 1
                if its == -1:
                    its = 0
                oldts = ts[il]
            
            #if its + 1 >= tera5.shape[2]:
                #continue
            w0=(secs[its+1]-ts[il])/(secs[its+1] - secs[its])
            w1=1.-w0
            if lats[il] != oldlats: 
                iyref = np.searchsorted(-tlat, -lats[il]) - 2
                if iyref == 2:
                    print(iyref)
                oldlats = lats[il]
            if lons[il] != oldlons:
                if tlon[-1] - tlon[0] > 0:
                    ixref = np.searchsorted(tlon, lons[il]) - 2
                    hlon = tlon
                else:
                    hlon = tlon.copy()
                    if lons[il] > hlon[0]:
                        
                        hlon[hlon<hlon[0]] += 360.
                    else:
                        hlon[hlon>hlon[-1]] -= 360.
                    ixref = np.searchsorted(hlon, lons[il]) - 2
                    #ixref = np.where(np.abs(hlon-lons[il]) <=0.5)[0][0] - 2
                if ixref <0:
                    ixref = 0
                oldlons = lons[il]
            
            for it in its, its + 1:
                if lats[il] > 90 - 2 * dy:
                    pval[it, ip] = tera5[0, 0, it, ip]          
                elif lats[il] < -90 + 2 * dy:
                    pval[it, ip] = tera5[-1, 0, it, ip]          
                else:
                    i = 0
                    for iy in range(4):
                        #print(yval.dtype, tlon.dtype, tera5[iy, :, it, ip:ip+1].T.dtype, lons[il])
                        #print(yval.shape, tlon.shape, tera5[iy, :, it, ip:ip+1].T.shape)
                        if ixref == 2:
                            print(ixref)
                        xx = polinl(hlon[ixref:ixref+4],tera5[iyref+iy, ixref:ixref+4, it, ip:ip+1].T,lons[il])
                        #print(xx.shape, xx.dtype)
                        yval[i, iy] = xx[0] #polinl(tlon,tera5[iy, :, it, ip:ip+1].T,lons[il])
                    xx = polinl(tlat[iyref:iyref+4], yval, lats[il])
                    pval[it, ip] = xx[0] #polinl(tlat, yval, lats[il])          
    
            reatab[im]=pval[its,ip]*w0+pval[its+1,ip]*w1
            idx[im] = il
            im += 1
    
    return idx[:im]

def offline_fb3(fpattern,p,pn, fdict,ts, obstype, latorig,lonorig,z, refs,ans):
    #from scipy.interpolate import RectBivariateSpline
    tt=time.time()
    
    if ans[3] == ans[2]: # no data
        return
    lons = lonorig[ans[2]:ans[3]]
    lats = latorig[ans[2]:ans[3]]
    tss = ts[ans[2]:ans[3]]
    obss = obstype[ans[2]:ans[3]]
    zs = z[ans[2]:ans[3]] / 100.
    try:

        fn = fpattern.format(ans[0],ans[1])
        with h5py.File(fn,'r') as g:
            
            if not fdict:
                fdict['level'] = g['level'][:]
                fdict['pidx'] = np.searchsorted(fdict['level'],refs['level'])
                
                fdict['longitude']= g['longitude'][:]
                fdict['latitude']= g['latitude'][:]
                
    
                tunits=g['time'].attrs['units']
                try:
                    tunits=tunits.decode('latin1')
                except:
                    pass
                fdict['tunits']=tunits.split()
    
            fdict['attribs'] =  dict(g[p].attrs)
            try:
                
                del fdict['attribs']['DIMENSION_LIST']
            except:
                pass
            for k,v in g.attrs.items():
                fdict['attribs'][k]=v
            fdict['attribs']['source_filename']=fn
                
            pres = fdict['level']
            pidx = fdict['pidx']
            sh = g[p].shape


            lonmin = np.min(lons)
            lonmax = np.max(lons)
            latmin = np.min(lats)
            latmax = np.max(lats)
            dx=360/sh[1]
            dy=(fdict['latitude'][0] - fdict['latitude'][-1])/(sh[0]-sh[0] % 2)
                
            if lonmin < 0:
                lonmin += 360.
                lons += 360.
            if lonmax < 0:
                lonmax += 360.
            ixrefmin=int(np.floor(lonmin/dx))
            ixrefmax=int(np.floor(lonmax/dx))
            if latmin== -90.:
                latmin = fdict['latitude'][-1]
            if latmax == 90.:
                latmax = fdict['latitude'][0]
            iyrefmin= int(np.floor((fdict['latitude'][0]-latmax)/dy))
            iyrefmax = int(np.floor((fdict['latitude'][0]-latmin)/dy))
            if ixrefmax >= ixrefmin:
                
                xsten = np.arange(ixrefmin-1, ixrefmax+3)#[ix - 1, ix, ix + 1, ix + 2])
            else:
                ixrefmax = int(np.floor(lonmax+360./dx))
                xsten = np.arange(ixrefmin-1, ixrefmax+3) % fdict['longitude'].shape[0]
            ysten = np.arange(iyrefmin-2, iyrefmax+3)#[iy - 1, iy, iy + 1, iy + 2])
            ysten[ysten>(sh[0]-1)] = sh[0] - 1
            ysten[ysten<0] = 0
            xsten[xsten>(sh[1]-1)] -= sh[1]
            xsten[xsten<0] += sh[1]

            if xsten[0] == xsten[-1] - xsten.shape[0] + 1:
                tera5=g[p][ysten[0]:ysten[-1] + 1,xsten[0]:xsten[-1] + 1,:,:]
                tlon = fdict['longitude'][:][xsten[0]:xsten[-1] + 1]
            else:
                tera5=g[p][ysten[0]:ysten[-1] + 1, :, :, :][:,xsten,:,:]
                tlon = fdict['longitude'][:][xsten]
            tlat = fdict['latitude'][ysten[0]:ysten[-1] + 1]
            
            try:
                mv=np.where(tera5==fdict['attribs']['missing_value'])
                tera5 = tera5 * np.float32(fdict['attribs']['scale_factor'])+np.float32(fdict['attribs']['add_offset'])
                tera5[mv]=np.nan
            except:
                pass

            if '20CRv3' in fpattern:
                for it in range(tera5.shape[2]):                       
                    for ip in range(len(pidx)):
                        if latmax > 90 - 2 * dy:
                            tera5[:,:,it,pidx[ip]]+=refs[p][ans[1] - 1, ysten[0]:ysten[-1] + 1,:, it % 8, ip][:, xsten]
                        elif latmin < -90 + 2 * dy:
                            tera5[:,:,it,pidx[ip]]+=refs[p][ans[1] - 1, ysten[0]:ysten[-1] + 1,:,it % 8, ip][:, xsten]
                        else:    
                            if xsten[0] == xsten[-1]- xsten.shape[0] + 1:
                                tera5[:,:,it,pidx[ip]]+=refs[p][ans[1] - 1, ysten[0]:ysten[-1] + 1,xsten[0]:xsten[-1] + 1,it % 8, ip]
                            else:
                                tera5[:,:,it,pidx[ip]]+=refs[p][ans[1] - 1, ysten[0]:ysten[-1] + 1,xsten,it % 8, ip].T
            #attribs['missing_value'] = np.nan
            #for k,v in g[p].attrs.items():
                #if k not in ['scale_factor','add_offset','DIMENSION_LIST','_FillValue','name','long_name','standard_name','units']:
                    #attribs[k]=v
                    #if k=='missing_value':
                        #attribs[k]=np.nan
            tunits = fdict['tunits']
            if tunits[0]=='hours':
                secs=np.int64(g['time'][:])*3600
            if tunits[-2]!='1900-01-01':
                if tunits[-2] =='1800-01-01':
                    offset=datetime.strptime(tunits[-2],'%Y-%m-%d')-datetime(1900,1,1) #+ timedelta(days=1)
                else:
                    offset=datetime.strptime(tunits[-2],'%Y-%m-%d')-datetime(1900,1,1)
                    
                secs=secs+int(offset.total_seconds())
                delta = secs[1] - secs[0]
                if datetime(1900, 1, 1) + timedelta(seconds=int(secs[0])) != datetime(ans[0], ans[1], 1):
                    secs = (datetime(ans[0], ans[1], 1) - datetime(1900, 1, 1) ).total_seconds() + np.arange(secs.shape[0]) *delta
                    print('modified timestamp of '+fn)
                

#        print(val[0, 0],'\n', ans[0, 0],'\n', time.time()-tt)
        
        an1 = [ans[0], ans[1] + 1]
        if an1[1] == 13:
            an1 = [ans[0] + 1, 1]
        try:
            
            with h5py.File(fpattern.format(an1[0],an1[1]),'r') as g:
                
                if xsten[0] == xsten[-1] - xsten.shape[0] + 1:
                    tera51=g[p][ysten[0]:ysten[-1] + 1,xsten[0]:xsten[-1] + 1,0:1,:]
                else:
                    tera51=g[p][ysten[0]:ysten[-1] + 1, :, :, :][:,xsten,0:1,:]
                
                try:
                    mv=np.where(tera51==fdict['attribs']['missing_value'])
                    tera51 =tera51 * np.float32(g[p].attrs['scale_factor'])+ np.float32(g[p].attrs['add_offset'])
                    tera51[mv]=np.nan
                except:
                    pass
    
                if '20CRv3' in fpattern:
                    for it in 0, :                       
                        for ip in range(len(pidx)):
                            if latmax > 90 - 2 * dy:
                                tera51[:,:,it,pidx[ip]]+=refs[p][an1[1] - 1, ysten[0]:ysten[-1] + 1,:, it % 8, ip][:, xsten]
                            elif latmin < -90 + 2 * dy:
                                tera51[:,:,it,pidx[ip]]+=refs[p][an1[1] - 1, ysten[0]:ysten[-1] + 1,:,it % 8, ip][:, xsten]
                            else:    
                                if xsten[0] == xsten[-1]- xsten.shape[0] + 1:
                                    tera51[:,:,it,pidx[ip]]+=refs[p][an1[1] - 1, ysten[0]:ysten[-1] + 1,xsten[0]:xsten[-1] + 1,it % 8, ip]
                                else:
                                    tera51[:,:,it,pidx[ip]]+=refs[p][an1[1] - 1, ysten[0]:ysten[-1] + 1,xsten,it % 8, ip].T
            tera5 = np.concatenate((tera5, tera51), axis=2)
            secs = np.concatenate((secs, secs[[-1]]+secs[-1]-secs[-2]))

        except FileNotFoundError as e:
            print(fpattern.format(an1[0],an1[1]), ' extrapolating last hours of month')
            tera5 = np.concatenate((tera5, tera5[:, :, -1:, :]), axis=2)
            secs = np.concatenate((secs, secs[[-1]]+secs[-1]-secs[-2]))
            pass

        #print(time.time() - tt)
        if False:
            #X, Y = np.meshgrid(g['longitude'], g['latitude'])
            val = np.empty(sh[2:])
            off = 2
            #print('lin', time.time() - tt)
            #tt = time.time()
            gt = g[p][iy +1 - off:iy + 1 + off, ix + 1- off:ix + 1 + off, :, :][::-1]
            glat = g['latitude'][iy + 1- off:iy + 1+ off][::-1].copy()
            glon = g['longitude'][ix + 1- off:ix + 1+ off].copy()

            for it in range(gt.shape[2]):

                for ip in range(gt.shape[3]):

                    interp_spline = RectBivariateSpline( glat,glon, gt[:, :, it, ip].reshape(gt.shape[:2]))
                    val[it, ip] = interp_spline(lat, lon)

            val =val*np.float32(g[p].attrs['scale_factor'])+np.float32(g[p].attrs['add_offset'])

            #print(time.time() - tt)
        
        reatab = np.zeros_like(lons)
        #params = {'t': 126}
        idx = triint(tera5[:, :, :, pidx], reatab, tlon, tlat, lons, lats, secs, tss,obss,pres[pidx],zs, pn, dy)

        print(ans,p, reatab[0], time.time()-tt)
    except FileNotFoundError as e:
        print(fn, 'not available, continuing')
        return 
    return reatab[:len(idx)], ans[2] + idx, secs,pres,fdict['attribs']

def offline_fb(fpattern,p,lat,lon,ts, refs,ans):
    from scipy.interpolate import RectBivariateSpline
    tt=time.time()
    try:

        with h5py.File(fpattern.format(ans[0],ans[1]),'r') as g:
            pres=g['level'][:]
            pidx=np.searchsorted(pres,refs['level'])
            dx=360/g[p].shape[1]
            dy=180/(g[p].shape[0]-1)
            ix=int(np.floor(lon/dx))  
            if lat==90.:
                lat=89.9999
            iy=int(np.floor((90-lat)/dy))
            #if ix==g[p].shape[1]-1:
                #tera5=np.empty_like(g['t'][iy:iy+2,ix:ix+2,:,:])
                #tera5[0,:,:,:]=g['t'][iy:iy+2,ix,:,:]
                #tera5[1,:,:,:]=g['t'][iy:iy+2,0,:,:]
                #if '20CRv3' in fpattern:
                    #for it in range(tera5.shape[2]):                       
                        #for ip in range(len(pidx)):
                            #tera5[0,:,it,pidx[ip]]+=refs[p][iy:iy+2,ix,it % 8, ip]
                            #tera5[1,:,it,pidx[ip]]+=refs[p][iy:iy+2,0,it % 8, ip]
            #else:
                #tera5=g[p][iy:iy+2,ix:ix+2,:,:]
                #if '20CRv3' in fpattern:
                    #for it in range(tera5.shape[2]):                       
                        #for ip in range(len(pidx)):
                            #tera5[:,:,it,pidx[ip]]+=refs[p][iy:iy+2,ix:ix+2,it % 8, ip]
            #try:
                #mv=np.where(tera5==g[p].attrs['missing_value'])
                #tera5=tera5*np.float32(g[p].attrs['scale_factor'])+np.float32(g[p].attrs['add_offset'])
                #tera5[mv]=np.nan
            #except:
                #pass
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

            #X, Y = np.meshgrid(g['longitude'], g['latitude'])
            val = np.empty(g[p].shape[2:])
            off = 2
            #print('lin', time.time() - tt)
            #tt = time.time()
            gt = g[p][iy - off:iy + off, ix - off:ix + off, :, :][::-1]
            glat = g['latitude'][iy - off:iy + off][::-1].copy()
            glon = g['longitude'][ix - off:ix + off].copy()
            print(time.time() - tt)
            #bounds = np.searchsorted(ts, secs)
            #tsi = ts[bounds[0]:bounds[1]+1]
            for it in range(gt.shape[2]):

                for ip in range(gt.shape[3]):

                    interp_spline = RectBivariateSpline( glat,glon, gt[:, :, it, ip].reshape(gt.shape[:2]))
                    val[it, ip] = interp_spline(lat, lon+360)

            val =val*np.float32(g[p].attrs['scale_factor'])+np.float32(g[p].attrs['add_offset'])    
            #print(ans, p, it, time.time() - tt)


#        if '20CRv3' in fpattern:
#            tera5
        #w=np.zeros(4)            
        #w[1]=lon-ix*dx
        #w[0]=1.-w[1]
        #w[2]=lat-(90-iy*dy)
        #w[3]=1.-w[2]
        #if tera5.shape[0]==1:     
            #tan=w[0]*tera5[0,1,:,:]+w[1]*tera5[0,0,:,:]
        #else:
            #tan=w[0]*w[2]*tera5[1,1,:,:]+w[1]*w[2]*tera5[1,0,:,:]+w[0]*w[3]*tera5[0,1,:,:]+w[1]*w[3]*tera5[0,0,:,:]

        print(ans,time.time()-tt)
    except Exception as e:
        print(e)
        return 
    return val,secs,pres,attribs

#@njit
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
    #   xout,yout = np.meshgrid(xout,yout)

    # compute grid coordinates of output grid.
    delx = xin[1:]-xin[0:-1]
    dely = yin[1:]-yin[0:-1]

    xcoords = (len(xin)-1)*(xout-xin[0])/(xin[-1]-xin[0])
    ycoords = (len(yin)-1)*(yout-yin[0])/(yin[-1]-yin[0])


    xcoords = np.clip(xcoords,0,len(xin)-1)
    ycoords = np.clip(ycoords,0,len(yin)-1)

    ## Interpolate to output grid using nearest neighbour
    if interpolation == 'NearestNeighbour':
        pass
        #xcoordsi = np.around(xcoords).astype(np.int32)
        #ycoordsi = np.around(ycoords).astype(np.int32)
        #dataout = datain[ycoordsi,xcoordsi]

    # Interpolate to output grid using bilinear interpolation.
    elif interpolation == 'Bilinear':
        
        xi = np.zeros_like(xcoords, dtype=int)
        yi = np.zeros_like(ycoords, dtype=int)
        yi = np.asarray(ycoords, dtype=int)
        xip1 = xi+1
        yip1 = yi+1
        xip1 = np.clip(xip1,0,len(xin)-1)
        yip1 = np.clip(yip1,0,len(yin)-1)
        delx = xcoords-xi.astype(np.float32)
        dely = ycoords-yi.astype(np.float32)
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
    reatab2 = np.zeros_like(reatab)
    eps = 1.e-5
    #dt=times[1]-times[0]
    for i in range(obstype.shape[0]):
        if obstype[i]==iobstype:
            if ts[i]!=tsold:
                ipress=0
                tsold=ts[i]
            zz=z[i]/100
            while press[ipress]<zz - eps:
                #ipressold=ipress
                if ipress==press.shape[0]-1:
                    break
                ipress+=1
            if np.abs(press[ipress]-zz) < eps:
                itimes=itimesold
                while times[itimes]<ts[i]:
                    itimesold=itimes
                    if itimes==times.shape[0]-1:
                        break
                    itimes+=1
                if itimes>0 and times[itimes]>=ts[i]:

                    w0=(times[itimes]-ts[i])/(times[itimes] - times[itimes-1])
                    w1=1.-w0
                    reatab[i]=values[itimes-1,ipress]*w0+values[itimes,ipress]*w1
                    if w0*w1 != 0:
                        if itimes>0 and itimes < times.shape[0] - 2:
                            reatab2[i] = polinl(times[itimes-2:itimes+2],values[itimes-2:itimes+2,ipress:ipress+1].T,ts[i])[0]
                        else:
                            reatab2[i] = reatab[i]

                    #if obstype[i]==126 and z[i]==30000:# and reatab[i]>400.:
                        #print(i,reatab[i],reatab2[i], ipress, itimes,w0,w1)

                else:
                    noval+=1
            else:
                #if ipress==press.shape[0]-1:
                ipress=0

    if noval>0:
        print('no reanalysis data for',noval,'obsvalues' )
    return

def Interp2dx(tup):
    return Interp2d(*tup)

def load_20CRoffset(readict):
       
    start = 1940
    try:
        with open(os.path.expandvars(wpath+'/rea/refs{}x.pkl'.format(start)),'rb') as f:
            refs=pickle.load(f)
    except:

        refs={}
        raw_20CR={}
        tt=time.time()
        l=0
        for iy in range(1940,1942):
            l+=1
            for im in range(1,13):
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

                                fpattern=v['path']+v['prefix']+ftype[:-1]+v['glue']+'{}{:0>2}'+v['glue']+pn+v['suffix']+'.nc'
                                if p in ['q', 'z']:
                                    continue
                            else:
                                fpattern=v['path']+v['prefix']+ftype+v['glue']+pn+'.reg_tl319.{}{:0>2}'+v['glue']+v['suffix']+'.nc'

                            try:
                                with h5py.File(fpattern.format(iy,im),'r') as f:
                                    print(f.keys())
                                    fpattern20=v20['path']+v20['prefix']+v20['glue']+'{}{:0>2}'+v20['glue']+v20['param'][p]+v20['suffix']+'.nc'
                                    with h5py.File(fpattern20.format(iy,im),'r') as g:
                                        print(g.keys())
                                        tfac = f['time'].shape[0] // g['time'].shape[0]
                                        tt = time.time()
                                        ifield=f[p][::tfac, :, :, :]
                                        print(time.time() - tt)
                                        ifields = np.zeros((ifield.shape[2], ifield.shape[3], 8, ifield.shape[1]), dtype=np.float32)
                                        for i in range(ifields.shape[2]):
                                            for ip in range(ifields.shape[3]):
                                                ifields[:, :, i, ip] = np.mean(ifield[i::8, ip, :, :], axis=0)*np.float32(f[p].attrs['scale_factor'])+np.float32(f[p].attrs['add_offset'])
                                        print(time.time() - tt)
                                        ilon=f['longitude'][:]
                                        ilat=f['latitude'][:]
                                        iplev=f['level'][:]

                                        oplev=g['level'][:]
                                        pindex=np.searchsorted(oplev,iplev)
                                        oshape=g[p].shape
                                        if p not in refs.keys():
                                            refs[p]=np.zeros((12, oshape[0],oshape[1],8, ifields.shape[3]),dtype=ifields.dtype)
                                            raw_20CR[p]=np.zeros(refs[p].shape, dtype=ifields.dtype)
                                            refs['level']=iplev[:]
                                        print(time.time() - tt)
                                        gmem = g[p][:]
                                        print(time.time() - tt)
                                        for i in range(ifields.shape[2]):
                                            for ip in range(ifields.shape[3]):

                                                raw_20CR[p][im - 1, :, :, i, ip]+=np.mean(gmem[:,:,i::8, pindex[ip]], axis=2)

                                        if 'olon' not in locals():
                                            glon = g['longitude'][:]
                                            glat = g['latitude'][:]
                                            olon,olat=np.meshgrid(g['longitude'][:],g['latitude'][:])

                                        print(time.time() - tt)



                                #tups = []
                                #for it in range(ifield.shape[2]):
                                    #for ip in range(ifield.shape[3]):
                                        #tups.append((ifield[:,:,it,ip], ilon,ilat,olon,olat,1))

                                #with multiprocessing.Pool(20) as pool:
                                    #res = pool.map(Interp2dx, tups)

                                ttt = time.time()
                                for it in range(ifields.shape[2]):
                                    for ip in range(ifields.shape[3]):
                                        interp_spline = RectBivariateSpline( -ilat, ilon, ifields[:, :, it, ip].reshape(ifields.shape[:2]))
                                        refs[p][im - 1, :,:,it, ip] += interp_spline(-glat, glon)
                                    print(iy,im,it, p,np.std(refs[p][im - 1, :, :, it, :]),np.std(ifields[:, :, it, :]),time.time()-tt)
                                print('x', time.time() - ttt)

                            except MemoryError as e:
                                print(e,'could not read')

        for r in ['t','u','v']:
            refs[r]/=l
            raw_20CR[r]/=l
            refs[r]-=raw_20CR[r]

        try:
            os.mkdir(wpath+'/rea')
        except:
            pass
        with open(os.path.expandvars(wpath+'/rea/refs{}x.pkl'.format(start)),'wb') as f:
            pickle.dump(refs,f)

    return refs    

ray_load_20CRoffset = ray.remote(load_20CRoffset)

def readictplot(readict, tasks, plevs, figprefix, marker='*'):
    
        
    for sfunc in tasks:

        for p in readict['obstypes'].keys():

            plev = plevs[0]

            idx=np.where(np.logical_and(readict['obs']['obstype']==readict['obstypes'][p],readict['obs']['z_coordinate']==plev))[0]
            if len(idx) < 50:
                continue

            for k,v in readict.items():
                if 'obs' in k:
                    continue
                if k =='era5fb':
                    i = 0
                
                if 'ftype' in v.keys():                    
                    iterable = v['ftype']
                else:
                    iterable = v
                    
                for ftype in iterable:  
                    if 'fc' in ftype:
                        dtype='fc'
                    else:
                        dtype='an'

                    if dtype in readict[k].keys():

                        rms=[]
                        years=[]
                        ref = datetime(1900, 1, 1)
                        tsy = ref.year + np.floor(readict['obs']['date_time'][idx]/365.25/86400)
                        
                        q = np.nanquantile(readict['obs']['obs'][idx]-readict[k][dtype]['refvalues'][idx], (0.01, 0.99))
                        print(p, k, q)
                        ystart = (ref + timedelta(seconds=int(readict['obs']['date_time'][idx[0]]))).year
                        ystop = (ref + timedelta(seconds=int(readict['obs']['date_time'][idx[-1]]))).year + 1
                        smarker = marker
                        if ystop - ystart == 1:
                            marker = '*'
                        
                        for iy in range(ystart, ystop):
                            #print(iy)

                            idy=np.where(tsy==iy)[0]
                            rv = np.clip(readict['obs']['obs'][idx[idy]] - readict[k][dtype]['refvalues'][idx[idy]], q[0], q[1])
                            if len(idy)>50 and np.sum(~np.isnan(rv)) > 50:
                                if sfunc == 'rms':

                                    rms.append(np.sqrt(np.nanmean((rv)**2)))
                                elif sfunc == 'std':
                                    rms.append(np.nanstd(rv))
                                else:
                                    rms.append(np.nanmean(rv))

                            else:
                                rms.append(np.nan)
                            years.append(iy)

                        plt.plot(np.array(years),np.array(rms),'-'+marker, 
                                 label='obs -'+k+'_'+dtype+', '+sfunc+'= {:5.3f}'.format(np.sqrt(np.nanmean(np.array(rms)**2))))

            plt.title('Monthly '+sfunc+' reanalysis '+p+' departures, '+str(int(plev/100)) + ' hPa, '+figprefix.split('/')[-1].split('_')[0])
            plt.ylabel(sfunc+' ['+readict['obsunits'][p]+']')
            plt.legend()
            plt.tight_layout()
            #plt.xlim(1935, 1965)
            plt.savefig(figprefix+'_'+p+'_'+sfunc+'stats.png')
            plt.close()

def retrieve_anfg(fn, readict, refs, out_name,path_to_gridded):
    from scipy.interpolate import RectBivariateSpline


    #era5=refdiff('/raid60/scratch/leo/scratch/ERA5/gridded/','era5fc.{0}{1:02}.'+par['gid'],1979,1983,tidx,fieldsperday=2,
                    #fcstep=12,mrange=list(range(1,13)),tgrid=None,tempname=tempname)

    #jra55_58_c[im,ipar,ip,:,:]=Interp2d(jra55_58_m[im,ipar,ip,:,:], jlon, jlat, eelon, eelat, order=1)

    #readict={'era5':{'ftype':('t','fct'),'param':{'t':'130','u':'131','v':'132','q':'133','z':'129'},
                     #'path':os.path.expandvars('$RSCRATCH/era5/gridded/'),'prefix':'era5','suffix':'','glue':'.'},
             ##'CERA20C':{'ftype':('t',),'param':{'t':'130','u':'131','v':'132','q':'133','z':'129'},
               ##'path':os.path.expandvars('$RSCRATCH/CERA20C/'),'prefix':'CERA20C','suffix':'','glue':'.'},
                          ##'JRA55':{'ftype':('fcst_mdl',),'param':{'t':'011_tmp','u':'033_ugrd','v':'034_vgrd','q':'051_spfh'},
               ##'path':os.path.expandvars('$RSCRATCH/JRA55/split/'),'prefix':'test.','suffix':'grb','glue':'.'},
                        #'20CRv3':{'ftype':('',),'param':{'t':'TMP','u':'UGRD','v':'VGRD'},#,'q':'SPFH'},#,'z':'HGT'},
                         #'path':os.path.expandvars('$RSCRATCH/20CRv3/'),'prefix':'anl_meant','suffix':'_pres','glue':'_'},
               #}

    cyear = int(out_name.split('/')[-2])
    outdir = os.path.dirname(out_name) + '/'
    with h5py.File(fn,'r') as f:
        try:
            if f['observations_table']['date_time'].shape[0] == 1:
                
                ts=f['observations_table']['date_time'][:]
            else:
                ts=f['observations_table']['date_time'][[0, -1]]
            ref = datetime(1900, 1, 1)
            tstart = int((datetime(cyear, 1, 1) - ref).total_seconds())
            tstop = int((datetime(cyear+1, 1, 1) - ref).total_seconds())
            if tstop < ts[0] or tstart > ts[-1]:
#                print(fn, cyear, 'year missing in obs records')
                return

            ts=f['observations_table']['date_time'][:]
            tslice =slice(*np.searchsorted(ts, (tstart, tstop)))
            if ts[tslice].shape[0] == 0:
                return
            else:
                ts = ts[tslice]
        except Exception as e:
            print(fn, cyear, e)
            return 

    #P=multiprocessing.Pool(12)
    with h5py.File(fn,'r') as f:
        try:

            #ts=f['observations_table']['date_time'][:]
            #ref = datetime(1900, 1, 1)
            #tstart = int((datetime(cyear, 1, 1) - ref).total_seconds())
            #tstop = int((datetime(cyear+1, 1, 1) - ref).total_seconds())
            #tslice =slice(*np.searchsorted(ts, (tstart, tstop)))
            #if ts[tslice].shape[0] == 0:
                #return
            #else:
                #ts = ts[tslice]
            
            lat=f['observations_table']['latitude'][tslice]
            lon=f['observations_table']['longitude'][tslice]
            lon[lon>360.] -= 360.
            lon[lon<0] += 360.
            obstype=f['observations_table']['observed_variable'][tslice]
            obs=f['observations_table']['observation_value'][tslice]
            z=f['observations_table']['z_coordinate'][tslice]
            latdis = lat + (100000. - z) * 3. / 100000.
            londis = lon + (100000. - z) * 3. / 100000.
            zt=f['observations_table']['z_coordinate_type'][tslice]
            if np.any(zt!=1):
                print('no pressure in', np.sum(zt!=1), 'observations')
                print('x')
            ofb=True
        except Exception as e:
            print(fn, cyear, e)
            return None
        try:

            bc = f['era5fb']['biascorr@body'][tslice]
            bc[np.isnan(bc)] = 0.
            o_minus_bg=f['era5fb']['fg_depar@body'][tslice] + bc
            o_minus_an=f['era5fb']['an_depar@body'][tslice] + bc
        except:
            ofb=False

    #fix geopotential
    tsu, ri=np.unique(ts, return_index=True)
    for i in range(len(ri) - 1):
        o = obstype[ri[i]:ri[i+1]]
        idx = np.where(np.logical_and(o==117, zt[ri[i]:ri[i+1]]==1))[0]
        if len(idx) > 0:
            m = np.argmax(obs[ri[i]+idx])
            #print(287 * 25 * np.log10(100000./z[ri[i]+idx[m]]) / obs[ri[i]+idx[m]])
            if 287 * 25 * np.log10(100000./z[ri[i]+idx[m]]) / obs[ri[i]+idx[m]] > 0.1:
                obs[ri[i]+idx] *=  9.80665
                if np.any(obs[ri[i]+idx]>500000):
                    print(i, len(ri), obs[ri[i]+idx[:]])
            #print(i, len(ri))
                    print('spurious')
        #if i % 1000 == 0:
            #print(i, len(ri))



    ref=datetime(1900,1,1)
    yms=[]
    oldyear=0
    oldmonth=0
    l = 0
    x = [ref+timedelta(seconds=int(k)) for k in tsu]
    while l < tsu.shape[0]:
        if (x[l].year!=oldyear or x[l].month!=oldmonth): #and x[l].year >1977 and x[l].year < 1983:
            m = 0
            
            oldmonth=x[l].month
            oldyear=x[l].year
            while x[l+m].year ==oldyear and x[l+m].month == oldmonth and l + m < tsu.shape[0] :
                m += 1
                if l + m ==len(x):
                    break
                
            if l +m == tsu.shape[0]:                
                yms.append((oldyear, oldmonth, ri[l], ts.shape[0]))
            else:
                yms.append((oldyear, oldmonth, ri[l], ri[l+m]))
            l += m - 1
        l += 1

    #try:
        #os.remove(out_name)
    #except:
        #pass
        
    obstypes={'t':ipar[85],'u':ipar[104],'v':ipar[105],'q':ipar[39],'z':ipar[117]}    
    obsunits={'t':'K','u':'m/s','v':'m/s','q':'g/kg','z':'m^2/s^2'}
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
                fdict = {}
                func=partial(offline_fb3,fpattern,p,obstypes[p], fdict,ts,obstype, lat,lon, z, refs)
#                func=partial(offline_fb3,fpattern,p,fdict, lat+latdis,lon+londis,ts,refs)
                #tups=list(P.map(func,yms[:]))
                tups=list(map(func,yms[:])) # no multiprocessing because function is already mapped
                ntups=[]
                for t in tups:
                    if t is not None:
                        ntups.append(t)
                if len(ntups)>0:

                    press=ntups[0][3]
                    if np.max(press)>1000:
                        press/=100.

                    readict[k][dtype][p]={}
                    readict[k][dtype][p]['time']=np.concatenate([ntups[i][2] for i in range(len(ntups))])
                    #readict[k][dtype][p]['values']=np.concatenate([ntups[i][0] for i in range(len(ntups))])
                    
                    for i in range(len(ntups)):
                        
                        readict[k][dtype]['refvalues'][ntups[i][1]] = ntups[i][0]
                        
                    readict[k][dtype][p]['attribs']=ntups[0][4]

            #tr[1]=117  # should change
            #tr[2]=85
            #tr[3]=104
            #tr[4]=105
            #tr[7]=39 #spec hum
            #tr[29]=38 #relative hum
            #tr[59]=36 # dew point
            #tr[111]=106 #dd
            #tr[112]=107  #ff

            if False and dtype in readict[k].keys():  
                tt=time.time()
                for p in readict[k][dtype].keys():
                    if p =='refvalues':
                        continue
                    #try:

                        #ttt = time.time()
                        #interp(readict[k][dtype]['refvalues'],obstype,z,ts,press,obstypes[p],
                               #readict[k][dtype][p]['time'],readict[k][dtype][p]['values'])
                        #print(k,dtype,p,time.time()-tt)
                    #except MemoryError as e:
                        #print(e)
                    plev = 30000
                    idx=np.where(np.logical_and(obstype==obstypes[p],z==plev))
                    if len(idx[0]) > 0:

                        try:

                            plt.subplot(2,1,1)
                            plt.plot(1900+ts[idx]/365.25/86400,obs[idx], label='obs')
                            rv = readict[k][dtype]['refvalues'][idx]
                            plt.plot(1900+ts[idx]/365.25/86400,rv, label=k)
                            plt.title(p+', '+fn.split('/')[-1].split('_')[0]+' '+str(plev//100) +' hPa')
                            plt.legend()
                            #plt.xlim(1935, 1965)
                            plt.subplot(2,1,2)
                            plt.plot(np.int32(1900+ts[idx]/365.25/86400),obs[idx]-rv,
                                     label='offline mean:{:5.3f} rms= {:5.3f}'.format(np.nanmean(obs[idx]-rv), np.sqrt(np.nanmean((obs[idx]-rv)**2))))
                            plt.plot(np.int32(1900+ts[idx]/365.25/86400),o_minus_an[idx]+rv-rv,
                                     label='online mean:{:5.3f} rms:{:5.3f}'.format(np.nanmean(o_minus_an[idx]+rv-rv), np.sqrt(np.nanmean((o_minus_an[idx]+rv-rv)**2))))
                            plt.title(p+', obs -'+k )
                            yscale = 10. ** np.floor(np.log10(np.max(np.abs(obs[idx]-rv)))) + 1
                            plt.ylim(-yscale, yscale)
                            plt.legend()
                            plt.tight_layout()
                            fnp=fn.split('/')[-1].split('CEUAS_merged_v1.nc')[0]
                            plt.savefig(outdir +fnp+p+'_'+k+dtype+'.png')
                            plt.close()

                        except Exception as e:

                            print('plotting reference', e)



    if not ofb:
        return readict

    readict['era5fb']={'ftype':['an','fc']}                  
    readict['era5fb']['an']={'refvalues':obs-o_minus_an }                  
    readict['era5fb']['fc']={'refvalues':obs-o_minus_bg }
    readict['obs'] = {'z_coordinate': z,'obstype': obstype,'date_time': ts, 'obs': obs,}
    readict['obstypes'] = obstypes
    readict['obsunits'] = obsunits

    readictplot(readict, ('mean', 'std','rms'), (30000, ), outdir+os.path.basename(fn).split('_')[0])

    #P.close()
    #P.join()
    #del P

    readict['tslice'] = tslice
    return readict

def dim_attach(g, k):
    for v in g[k].keys(): #var_selection:
        l=0            
        try:
            fvv=g[k][v]
            if 'string' not in v and v!='index':                    
                g[k][v].dims[l].attach_scale(g[k]['index'])
                #print(v,fvv.ndim,type(fvv[0]))
                if fvv.ndim==2 or type(fvv[0]) in [str,bytes,np.bytes_]:
                    slen=fvv.shape[1] #sdict[v]
                    #slen=10
                    g[k][v].dims[1].attach_scale(g[k]['string{}'.format(slen)])
        except MemoryError as e:
            print(g.filename.split('/')[-1],k, e)
            pass


    #fpattern=path_to_gridded+'/era5fct.{}{:0>2}.130.nc'
    #func=partial(offline_fb,fpattern,lat,lon)
    #tfgs=list(map(func,yms))

def convert_missing(refs, wpath,cyear, fn):

    print(fn.split('/')[-1], cyear,end=' ' )
    wpathy = wpath+'/' + str(cyear) + '/'
    targetfile = wpathy + fn.split('/')[-1]
    if os.path.isfile(wpathy+'/log/'+fn.split('/')[-1]+".txt"):
        print('already processed')
        return targetfile
    else:
        print('processing...')
    sys.path.insert(0,os.getcwd()+'/../resort/rasotools-master/')
    
    readict={'era5':{'ftype':('t','fct'),'param':{'t':'130','u':'131','v':'132','q':'133','z':'129'},
                     'path':os.path.expandvars('$RSCRATCH/era5/gridded/'),'prefix':'era5','suffix':'','glue':'.'},
             #'CERA20C':{'ftype':('t',),'param':{'t':'130','u':'131','v':'132','q':'133','z':'129'},
               #'path':os.path.expandvars('$RSCRATCH/CERA20C/'),'prefix':'CERA20C','suffix':'','glue':'.'},
                          #'JRA55':{'ftype':('fcst_mdl',),'param':{'t':'011_tmp','u':'033_ugrd','v':'034_vgrd','q':'051_spfh'},
               #'path':os.path.expandvars('$RSCRATCH/JRA55/split/'),'prefix':'test.','suffix':'grb','glue':'.'},
                        '20CRv3':{'ftype':('',),'param':{'t':'TMP','u':'UGRD','v':'VGRD'},#,'q':'SPFH'},#,'z':'HGT'},
                         'path':os.path.expandvars('$RSCRATCH/20CRv3/'),'prefix':'anl_meant','suffix':'_pres','glue':'_'},
               }
    

    tt=time.time()
    nanlist = [float('nan'), np.nan, 0, -2147483648]


    rscratch='/mnt/users/scratch/leo/scratch/'
    try:
        reaname = os.path.expandvars(wpathy+fn.split('/')[-1].split('x_CEUAS_merged_v1.nc')[0]+'_'+str(cyear) + '.pkl')
        with open(reaname,'rb') as f:
            readict=pickle.load(f)
    except:

        out_name = wpathy + fn.split('/')[-1]  

        path_to_gridded=os.path.expandvars(rscratch+'/era5/gridded/')
        readict=retrieve_anfg(fn,readict, refs, out_name,path_to_gridded)
        if readict is None:
            print(fn.split('/')[-1], 'no data for year', cyear)
            return
        with open(reaname,'wb') as f:
            pickle.dump(readict,f)
    print (time.time()-tt)
    
     # wpath+fn.split('/')[-1]
    try:
        os.mkdir(wpathy)
    except:
        pass
    if os.path.isfile(targetfile):
        try:
            os.remove(targetfile)
        except:
            print('file could not be removed - overwriting will lead to errors')

    #print('writing reanalysis reference series')
    #for k in readict.keys():
        #if k in ('tslice', ):
            #continue
        #for dtype in 'an', 'fc':
            
            #try:
                #os.mkdir(os.path.basename(targetfile))
            #except:
                #pass
            #try:
        
                #df = {dtype:readict[k][dtype]['refvalues']}  # making a 1 column dataframe
                #write_dict_h5(targetfile, df, k, {dtype: {'compression': 'gzip'}}, 
                              #var_selection=[], mode='a', attrs = {dtype:readict[k][dtype]['t']['attribs']} )  
            #except Exception as e:
                #print(e, 'no values from ',k,'for station ',os.path.basename(fn))


    tslice = readict['tslice']
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
            index = data.observations_table.index.shape[0]
            if o in ['observed_variable','observation_value','observation_id','z_coordinate','z_coordinate_type','date_time','conversion_flag','conversion_method']:  
                if len(data.observations_table[o].shape)==1:
                    loaded_data.append((data.observations_table[o][tslice]))
                    a_loaded_data.append(np.empty_like(loaded_data[-1],shape=2*len(loaded_data[-1])+addmem))
                else:
                    #loaded_data.append(data.observations_table[o][:index, :].view('S{}'.format(data.observations_table[o].shape[1])).flatten())
                    loaded_data.append(np.zeros(tslice.stop-tslice.start, dtype='S21'))
#                     a_loaded_data.append(np.empty((2*len(loaded_data[-1]),len(loaded_data[-1][0])), dtype=loaded_data[-1][0].dtype))
                    a_loaded_data.append(np.empty_like(loaded_data[-1],shape=2*len(loaded_data[-1])+addmem))
                    a_loaded_data[-1].fill(b' '*data.observations_table[o].shape[1])
                loaded_type['names'].append(o)
                loaded_type['formats'].append(loaded_data[-1].dtype)
                ld.append((o,loaded_data[-1].dtype))
        try:

            loaded_obstab = np.rec.fromarrays(loaded_data, dtype=ld)
        except Exception as e:

            for l in loaded_data:
                print(l.shape)
            print(fn, e)
            return None

        del loaded_data
        a_loaded_obstab = np.rec.fromarrays(a_loaded_data, dtype=ld)
        del a_loaded_data

        if ofb:
            loaded_fb=[]
            a_loaded_fb=[]
            loaded_type = {'names':[],'formats':[]}
            lf=[]
            for o in fbkeys:
                if o in ['fg_depar@body','an_depar@body','biascorr@body','biascorr_fg@body']:  
                    loaded_fb.append((data.era5fb[o][tslice]))
                    a_loaded_fb.append(np.zeros_like(loaded_fb[-1],shape=2*len(loaded_fb[-1])+addmem))
                    if o in  ['fg_depar@body','an_depar@body']:
                        a_loaded_fb[-1].fill(np.nan)
                    loaded_type['names'].append(o)
                    loaded_type['formats'].append(loaded_fb[-1].dtype)
                    lf.append((o,loaded_fb[-1].dtype))
            loaded_feedback = np.rec.fromarrays(loaded_fb, dtype=lf)
            del loaded_fb
            a_loaded_feedback = np.rec.fromarrays(a_loaded_fb, dtype=lf)
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
            if(np.any(loaded_feedback['fg_depar@body']>1.e26) or np.any(readict['era5']['fc']['refvalues']>1.e26)):
               print(cyear)
            add_fb(loaded_obstab,loaded_feedback,readict['20CRv3']['an']['refvalues'],
                   readict['era5']['an']['refvalues'],readict['era5']['fc']['refvalues'])

        del readict
        dr = data.recordindex[:]
        rslice = slice(*np.searchsorted(dr, (tslice.start, tslice.stop)))
        recordindex = dr[rslice]
        del dr

        # --->

    print(time.time()-tt)
    idx,press,temp,relhum,spechum,dpd,dewpoint,uwind,vwind,wd,ws,\
    d_temp,d_relhum,d_spechum,d_dpd,d_dewpoint,d_uwind,d_vwind,d_wd,d_ws,\
    fgd_temp,fgd_relhum,fgd_spechum,fgd_dpd,fgd_dewpoint,fgd_uwind,fgd_vwind,fgd_wd,fgd_ws=ipl2(loaded_obstab, loaded_feedback)

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

    idy=np.where(loaded_obstab['z_coordinate_type'][idx]==2) # do not convert humidity if data are not on pressure coordinates
    for c in cdpdrh,cshrh,cshdpd,crhdpd:
        c[idy]=np.nan

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

    humvar=np.array((ipar[34],ipar[36],ipar[38],ipar[39])) #dpd,dp,rh,sh
    wvar=np.array((ipar[104],ipar[105],ipar[106],ipar[107])) #dpd,dp,rh,sh

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

    with h5py.File(fn, 'r') as file:
        with h5py.File(targetfile, 'a') as newfile:

            headerslice = slice(*np.searchsorted(file['header_table']['record_timestamp'][:], (rt[0], rt[-1]+1)))

            groups = []
            for i in file.keys():
                if type(file[i]) == h5py._hl.group.Group:
                    if i not in ('observations_table','era5fb'):
                        newfile.create_group(i)
                        groups.append(i)
                elif i in ('recordindex', 'recordtimestamp', 'dateindex'):
                    pass
                else:
                    newfile.create_dataset(i, data=file[i][:])
            for i in groups:
                if(i == 'recordindices' or i == 'observations_table' or i == 'era5fb'):
                    pass
                elif i in ('header_table', 'source_configuration'):
                    if 'index' not in file[i].keys():
                        newfile[i].create_dataset('index', data=np.empty(headerslice.stop-headerslice.start, dtype='S1'))                       
                    for j in file[i].keys():
                        print(i, j)
                        if 'string' in j:
                            xdata = file[i][j][:]
                        else:
                            xdata = file[i][j][headerslice][:]
                        newfile[i].create_dataset(j, data=xdata)
                else:
                    if 'index' not in file[i].keys():
                        sh = file[i][list(file[i].keys())[0]].shape[0]
                        newfile[i].create_dataset('index', data=np.empty(sh, dtype='S1'))                       
                    for j in file[i].keys():
                        newfile[i].create_dataset(j, data=file[i][j][:])
                    
                dim_attach(newfile, i)

    obsv = out['observed_variable'][:jj]
    allvars = np.sort(np.unique(obsv))
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
    #pres = out['z_coordinate'][:jj][absidx]
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
            ov_vars = np.empty(jj,dtype='S11')
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

                print(i, time.time()-tt)
                with eua.CDMDataset(fn) as data:
                    #i='z_coordinate'
                    rest_data = data.observations_table[i][tslice]
                if rest_data.ndim==2: #i in ['observation_id', 'report_id', 'sensor_id', 'source_id']:
                    final = np.empty((addedvar[-1][1],len(rest_data[0])), dtype=rest_data[0].dtype)
                else:
                    final = np.empty(addedvar[-1][1], dtype=rest_data.dtype)

                #print('vor fill_restdata', fn)
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
                        rest_data = data.era5fb[i][tslice]
                    if i in ['expver', 'source@hdr', 'source_id', 'statid@hdr']:
                        final = np.empty((addedvar[-1][1],len(rest_data[0])), dtype=rest_data[0].dtype)
                    else:
                        final = np.empty(addedvar[-1][1], dtype=rest_data[0].dtype)
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
    try:
        os.mkdir(wpathy+'/log')
    except:
        pass
    f= open(wpathy+'/log/'+fn.split('/')[-1]+".txt","w")
    f.write("done") 
    f.close()
    return targetfile

ray_convert_missing = ray.remote(convert_missing)

def rmsplot(files):
    
    readicts = []
    print('rmsplot', files)
    for file in files:
        year = os.path.dirname(file).split('/')[-1]
        pkname = file +'_'+year + '.pkl'
        
        print(pkname)
        try:
            
            with open(pkname,'rb') as f:
                readicts.append(pickle.load(f))
        except:
            print (pkname + ' not found')
    if len(readicts) == 0:
        print(pkname+' no feedback information found')
        return
    
    readict = {'obs' : {}}
    print(readicts[0].keys())
    for k in readicts[0]['obs'].keys():
        readict['obs'][k] = np.concatenate([i['obs'][k] for i in readicts])
    for p in readicts[0]['obstypes'].keys():

        for k,v in readicts[0].items():
            if k in ('obs', 'tslice'):
                continue
            if 'obs' in k:
                readict[k] = readicts[0][k]
                continue
            if k =='era5fb':
                i = 0
            for ftype in v['ftype']:  
                if 'fc' in ftype:
                    dtype='fc'
                else:
                    dtype='an'
                if k not in readict.keys():
                    readict[k] = {}
                if dtype not in readict[k].keys():
                    readict[k][dtype] = {}
                readict[k][dtype]['refvalues'] = np.concatenate([i[k][dtype]['refvalues'] for i in readicts])


    readictplot(readict, ('mean', 'std','rms'), (30000, ), wpath+'/plots/'+os.path.basename(pkname).split('.')[0].split('_')[0], marker='')
    return

ray_rmsplot = ray.remote(rmsplot)

def do_station(refs_ref, file, years):
    
    futures = [ray_convert_missing.remote(refs_ref, wpath, year, file ) for year in years]
    

    return ray.get(ray_rmsplot.remote(futures))

ray_do_station = ray.remote(do_station)

def readgridded(readict, refs_ref, year, month):
    fdict = {}
    for k,v in readict.items():
        fdict[k] = {}
        if k!='20CRv3':
            continue
        for ftype in v['ftype']:  
            if 'fc' in ftype:
                dtype='fc'
            else:
                dtype='an'
            if dtype not in readict[k].keys():
                readict[k][dtype]={}
                #readict[k][dtype]['refvalues']=np.empty(obs.shape,dtype=np.float32)
                #readict[k][dtype]['refvalues'].fill(np.nan)
            fdict[k][dtype] = {}
            for p,pn in v['param'].items():
                if k!='JRA55':
        
                    fpattern=v['path']+v['prefix']+ftype+v['glue']+'{}{:0>2}'+v['glue']+pn+v['suffix']+'.nc'
                else:
                    fpattern=v['path']+v['prefix']+ftype+v['glue']+pn+'.reg_tl319.{}{:0>2}'+v['glue']+v['suffix']+'.nc'
            
                if k != '20CRv3':
                    fn = fpattern.format(year,month)
                    zname = 'level'
                else:
                    fn =  v['path']+v['prefix']+ftype+v['glue']+'{}'.format(year)+v['glue']+pn+v['suffix']+'.nc'
                    zname = 'isobaricInhPa'
                    
                
                with h5py.File(fn,'r') as g:
                    
                    if p not in fdict[k][dtype].keys():
                        if k != '20CRv3':
                            fdict['level'] = g['level'][:]
                            fdict['pidx'] = np.searchsorted(fdict['level'],refs['level'])
                        else:
                            fdict['level'] = g['isobaricInhPa'][:]                                  
                            fdict['pidx'] = np.searchsorted(-fdict['level'], -refs['level'])
                        
                        fdict['longitude']= g['longitude'][:]
                        fdict['latitude']= g['latitude'][:]
                        
                        fdict['attribs'] =  dict(g[p].attrs)
                        try:
                            
                            del fdict['attribs']['DIMENSION_LIST']
                        except:
                            pass
                        for kk,vv in g.attrs.items():
                            fdict['attribs'][kk]=vv
                        fdict['attribs']['source_filename']=fn
                        fdict[k][dtype][p] = g[p][:][:, fdict['pidx'], :, :]* np.float32(fdict['attribs']['scale_factor'])+np.float32(fdict['attribs']['add_offset'])
                        
                        print(k, dtype, p)
                        tunits=g['time'].attrs['units']
                        try:
                            tunits=tunits.decode('latin1')
                        except:
                            pass
                        fdict['tunits']=tunits.split()
            
        
    return fdict
        #found=len(glob.glob(fpattern.format(*yms[0])))
        #print(fpattern.format(*yms[0]),found)
        #if found:
    fdict = {}
    func=partial(offline_fb3,fpattern,p,obstypes[p], fdict,ts,obstype, lat,lon, z, refs)

ray_readgridded =ray.remote(readgridded)   

def plot_contents(wpath,cyear, fn):

    pardict = {'0': ['unknown',''], 
               '34': ['dewpoint depression','K'],
               '39': ['specific humidity','kg/kg'],
               '139': ['u-component of wind','m/s'],
               '140': ['v-component of wind','m/s'],
               '106': ['wind from direction','deg'],
               '107': ['wind speed','m/s'],
               '117': ['geopotential','J/kg'],
               '126': ['temperature','K'],
               '137': ['dewpoint','K'],
               '138': ['relative humidity','']
               }
    print(fn.split('/')[-1], cyear,end=' ' )
    wpathy = wpath+'/' + str(cyear) + '/'
    targetfile = wpathy + fn.split('/')[-1]
    #if os.path.isfile(wpathy+'/log/'+fn.split('/')[-1]+".txt"):
        #print('already processed')
        #return targetfile
    #else:
        #print('processing...')
    sys.path.insert(0,os.getcwd()+'/../resort/rasotools-master/')
    
    with h5py.File(targetfile, 'r') as f:
        rk = list(f['recordindices'].keys())[:-2]
        for pl in [70000.]:
            plt.figure(figsize=(10, (len(rk) +1) *1.5))
            l = 1
            for v in rk:
                tslice = slice(*f['recordindices'][v][[0, -1]])
                ts = f['observations_table']['date_time'][:][tslice]
                obs = f['observations_table']['observation_value'][tslice]
                pidx = np.where(f['observations_table']['z_coordinate'][tslice]==pl)[0]
                if len(pidx) < 2:
                    continue
                plt.subplot( len(rk),1, l)
                plt.plot(np.asarray(1900+ts[pidx]/86400/365.25, dtype='int'), rmeanw(obs[pidx], 30))
                #plt.plot(np.asarray(1900+thin2(ts[pidx], 10)/86400/365.25, dtype='int'), thin2(rmeanw(obs[pidx], 30), 10))
                plt.title(os.path.basename(targetfile).split('_')[0]+','+str(np.int(pl/100.)) + 'hPa, '+ pardict[v][0])
                plt.ylabel(pardict[v][1])
                l += 1
        plt.tight_layout()
        plt.savefig(targetfile[:-3]+'.png')
        plt.close()
                
                
        

# files = glob.glob('/raid60/scratch/federico/MERGED_DATABASE_OCTOBER2020_sensor/0-20000-0-01*.nc')
# files = glob.glob('/raid60/scratch/federico/DATABASE_JANUARY2021_sensor/0-20000-0-97690*.nc')
# files = glob.glob('/raid60/scratch/federico/DATABASE_JANUARY2021_sensor/0-20500-0-93954*.nc')
# files = glob.glob('/raid60/scratch/federico/DATABASE_JANUARY2021_sensor/0-20400-0-04665*.nc')
# files = glob.glob('/raid60/scratch/federico/DATABASE_JANUARY2021_sensor/*.nc')

# print(files[:10])

# convert_missing(files[6020])
# convert_missing('/raid60/scratch/federico/MERGED_DATABASE_OCTOBER2020_sensor/0-20000-0-03414_CEUAS_merged_v0.nc')

if __name__ == '__main__':


    wpath= os.path.expandvars('$RSCRATCH/converted_v11/') #'./'
    opath=wpath
    wlpath=wpath+'log/'
    for p in wpath,wlpath, wpath + '/long', wpath + '/rea':      
        try:
            os.mkdir(p)
        except:
            pass

    #for i in no_height:
        #do_resort(i)
        #f= open(wlpath+i.split('/')[-1]+".txt","w+")
        #f.write("done") 
        #f.close()

#    files = glob.glob('/mnt/scratch/scratch/federico/MERGED_15JUNE2022/*02365*v1.nc')
#    files = glob.glob('/mnt/scratch/scratch/federico/MERGING_DEC2022_FIXED_1/*27612*v1.nc')
#    files = glob.glob('/mnt/scratch/scratch/federico/MERGING_DEC2022_FIXED_1/*-20???-0-[12]*v1.nc')
    #files = glob.glob('/mnt/scratch/scratch/federico/MERGING_DEC2022_FIXED_1/*-0-*v1.nc')
    #files = glob.glob('/mnt/scratch/scratch/federico/MERGING_DEC2022_FIXED_1/*-20???-0-82930*v1.nc')
    #files = glob.glob('/mnt/scratch/scratch/federico/VIENNA_SENSORFIX_JAN2023/*-20???-0-11035*v1.nc')
    files = glob.glob('/mnt/scratch/scratch/federico/YEAR_SPLIT_MERGING/[12]???/0-*v2.nc')
#    files = glob.glob('/mnt/scratch/scratch/federico/COP2_HARVEST_NOVEMBER2022/era5_1_mobile/*ASEU03*.nc')
    #files_to_convert = glob.glob('/raid60/scratch/federico/MERGED_JUNE2021/*72357*v1.nc')
    #files = glob.glob('/scratch/das/federico/TRY_MERGED_FEB2022/*89564*v1.nc')
    #files_to_convert = files #glob.glob('/scratch/das/federico/TRY_MERGED_JAN2022/*72357*v1.nc')

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

    print(files[0], files[-1])

    #already_done = glob.glob(wlpath+'*_yearly.txt')

    #files_to_convert = []#files #[]
    #for i in files:
        #if not wlpath+i.split('/')[-1]+'_yearly.txt' in already_done:
            #files_to_convert.append(i)
    files_to_convert = files
    #files_to_convert.sort()
    tt=time.time()

    
    readict={'era5':{'ftype':('t','fct'),'param':{'t':'130','u':'131','v':'132','q':'133','z':'129'},
                     'path':os.path.expandvars('$RSCRATCH/era5/gridded/'),'prefix':'era5','suffix':'','glue':'.'},
             #'CERA20C':{'ftype':('t',),'param':{'t':'130','u':'131','v':'132','q':'133','z':'129'},
               #'path':os.path.expandvars('$RSCRATCH/CERA20C/'),'prefix':'CERA20C','suffix':'','glue':'.'},
                          #'JRA55':{'ftype':('fcst_mdl',),'param':{'t':'011_tmp','u':'033_ugrd','v':'034_vgrd','q':'051_spfh'},
               #'path':os.path.expandvars('$RSCRATCH/JRA55/split/'),'prefix':'test.','suffix':'grb','glue':'.'},
                        '20CRv3':{'ftype':('',),'param':{'t':'TMP','u':'UGRD','v':'VGRD'},#,'q':'SPFH'},#,'z':'HGT'},
                         'path':os.path.expandvars('$RSCRATCH/20CRv3/'),'prefix':'anl_meant','suffix':'_pres','glue':'_'},
               }
    readict={'era5':{'ftype':('','fc'),'param':{'t':'130','u':'131','v':'132','q':'133','z':'129'},
                     'path':os.path.expandvars('$RSCRATCH/era5/gridded/'),'prefix':'era5','suffix':'','glue':'.'},
             #'CERA20C':{'ftype':('t',),'param':{'t':'130','u':'131','v':'132','q':'133','z':'129'},
               #'path':os.path.expandvars('$RSCRATCH/CERA20C/'),'prefix':'CERA20C','suffix':'','glue':'.'},
                          #'JRA55':{'ftype':('fcst_mdl',),'param':{'t':'011_tmp','u':'033_ugrd','v':'034_vgrd','q':'051_spfh'},
               #'path':os.path.expandvars('$RSCRATCH/JRA55/split/'),'prefix':'test.','suffix':'grb','glue':'.'},
                        '20CRv3':{'ftype':('',),'param':{'t':'TMP','u':'UGRD','v':'VGRD'},#,'q':'SPFH'},#,'z':'HGT'},
                         'path':os.path.expandvars('$RSCRATCH/20CRv3/'),'prefix':'anl_meant','suffix':'_pres','glue':'_'},
               }
    
    
    refs = load_20CRoffset(readict)
    readicts = []
    for i in range(2022, 1957, -1):
        
        func=partial(convert_missing, refs, wpath, i)
        result_list = list(map(func, files_to_convert[:]))
        
    # simple plotting
    for i in  ('long', ):#(1995, 1994, -1):
        
        func=partial(plot_contents, wpath, i)
        result_list = list(map(func, files_to_convert[:]))

    exit()
               

    ray.init(num_cpus=40, _temp_dir=os.path.expanduser('~/ray'))
    #refs_ref = ray.put(refs)

    #for year in range(2015, 1904, -1):
            
        #futures = [ ray_convert_missing.remote(refs_ref, wpath, year, file )  for file in files_to_convert]
    #obj_result_list = ray.get(futures)

    pwd = os.getcwd()
    wpath = os.path.expandvars('$RSCRATCH/converted_v10/')
#    os.chdir(wpath)
    obj_result_list = glob.glob(wpath+'/[12]???/*v1.nc')
    os.chdir(pwd)
    shlist = np.unique([os.path.basename(ll) for ll in obj_result_list if ll is not None])
    futures = []
    vienna = False
    for file in shlist:
        if '11035' in file:
            vienna = True
        sublist = sorted([s for s in obj_result_list if s is not None and file in s])
        #rmsplot(sublist)
        if vienna:
            futures.append(ray_rmsplot.remote(sublist))
    ray.get(futures)

    #futures = []
    #for file in files_to_convert:
        #futures .append(ray_do_station.remote(refs_ref, file, range(2022, 1904, -1)))
#     print(result_list)
    print('total:',time.time()-tt)
    l = 0

