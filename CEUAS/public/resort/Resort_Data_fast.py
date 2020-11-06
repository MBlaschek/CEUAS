#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas
import numpy as np
from numba import njit
import sys
import zipfile, os, time
import urllib3
from datetime import datetime, timedelta
import glob
import h5py
import plotly.express as px
sys.path.append(os.getcwd()+'/../cds-backend/code/')
sys.path.append(os.getcwd()+'/../harvest/code/')
from  harvest_convert_to_netCDF_newfixes import write_dict_h5
import cds_eua3 as eua
eua.logging_set_level(30)
import xarray as xr

import cdsapi, zipfile, os, time
import schedule
import copy


# In[2]:


import numpy
import pandas as pd


# In[3]:


files = glob.glob('/raid60/scratch/federico/MERGED_DATABASE_OCTOBER2020_sensor/0-20000-0-11035*.nc')
#files = glob.glob('/raid60/scratch/uli/0-20000-0-96*.nc')
files[0]


# In[4]:


data = eua.CDMDataset(files[0])
data


# In[5]:


allvars = data.observations_table.observed_variable[()]
np.unique(allvars)


# In[ ]:


allvars = data.observations_table.observed_variable[()]
allvars.sort()
allvars = np.unique(allvars)
#
ri = data.recordindex[()]
print('recordindex: ', len(ri))
rt = data.recordtimestamp[()]
keys = data.observations_table.keys()[:-1]
# dropping all keys, where dimensions won't work - just help variabels for dimensions
keys.pop(43)
keys.pop(42)
keys.pop(41)
recordindices = [[] for i in range(len(allvars))]
recordtimestamps = [[] for i in range(len(allvars))]

# output variables (from observations_table)
ov = []
for o in keys:
    ov.append([[] for i in range(len(np.unique(allvars)))])

#
# load data into memory and decode byte arrays
#
print('loading data')
obsv = data.observations_table.observed_variable[:]
ov_vars = []
for o in range(len(keys)):
    ov_vars.append(data.observations_table[keys[o]][:])

#
# resorting the data
#
print('resort:start')
@njit
def make_vrindex(vridx,ridx,idx):
    l=0
    for i in range(1,len(idx)):
        if ridx[i]>ridx[i-1]:
            vridx[ridx[i-1]:ridx[i]]=l
            l=i
    vridx[ridx[i]:]=len(idx)       
    

tt=time.time()
ridxall=np.zeros(ov_vars[0].shape[0],dtype=np.int64)
for j in range(len(ri)-1):
    ridxall[ri[j]:ri[j+1]]=j
j+=1
ridxall[ri[j]:]=j
ridx=[]
vridx=[]
absidx=[]
abscount=0
for j in range(len(allvars)):
    idx=np.where(ov_vars[keys.index('observed_variable')]==allvars[j])[0]
    print(j,len(idx))
    vridx.append(np.zeros(ri.shape[0],dtype=np.int64))
    ridx=ridxall[idx]
    make_vrindex(vridx[-1],ridx,idx)
    vridx[-1]+=abscount
    #ridx.append(np.unique(ridxall[idx]))
    #vridx.append(np.zeros(ri.shape[0],dtype=np.int64))
    #vridx[-1][ridx[-1]]=abscount+np.arange(ridx[-1].shape[0],dtype=np.int64)
    #for i in range(len(vridx[-1])):
        #if vridx[-1][i]==-999:
            #vridx[-1][i]=vridx[-1][i-1]
    #for i in range(len(vridx[-1])-1,-1,-1):
        #if vridx[-1][i]==-999:
            #vridx[-1][i]=vridx[-1][i+1]
    
    absidx.append(copy.copy(idx))
    abscount+=len(idx)

absidx=np.concatenate(absidx)
for o in range(len(keys)):
    ov_vars[o]=ov_vars[o][absidx]
    print(o)
#
# recordtimestamps are only necessary once
#
recordtimestamps = recordtimestamps[0]

# 
# stacking all output variables
#
print('stacking output variables')


# targetfile has to be a copy of the original file
targetfile = '/raid60/scratch/leo/scratch/'+files[0].split('/')[-1]
    

#
# writing data into observations_table
#
for i in 0,14,25,26,27,44:#range(1): #len(keys)):
    try:
        with h5py.File(targetfile, 'r+') as file:
            del file['observations_table'][keys[i]]
    except:
        pass
    if keys[i] == 'index':
        pass
    elif keys[i] == 'observation_id' or keys[i] == 'report_id' or keys[i] == 'sensor_id' or keys[i] == 'source_id':
        slen = len(ov_vars[i][0])
        alldict = {keys[i]:np.asarray(ov_vars[i], dtype='S1')}
        write_dict_h5(targetfile, alldict, 'observations_table', {keys[i]: { 'compression': 'gzip' } }, [keys[i]])
    else:
        alldict = pandas.DataFrame({keys[i]:ov_vars[i]})
        write_dict_h5(targetfile, alldict, 'observations_table', {keys[i]: { 'compression': 'gzip' } }, [keys[i]])  
    print(i,keys[i])
#
# writing the recordindices and recordtimestamp.
#       
recordindices=vridx
for i in range(len(recordindices)):
    try:
        with h5py.File(targetfile, 'r+') as file:
            del file['recordindices'][str(allvars[i])]
    except:
        pass
    testvar = pandas.DataFrame({str(allvars[i]):recordindices[i]})
    write_dict_h5(targetfile, testvar, 'recordindices', {str(allvars[i]): { 'compression': None } }, [str(allvars[i])]) 

try:
    with h5py.File(targetfile, 'r+') as file:
        del file['recordindices']['recordtimestamp']
except:
    pass
testvar = pandas.DataFrame({'recordtimestamp':rt})
write_dict_h5(targetfile, testvar, 'recordindices', {'recordtimestamp': { 'compression': None } }, ['recordtimestamp']) 

try:
    with h5py.File(targetfile, 'r+') as file:
        del file['recordindex']
except:
    pass


# In[9]:


#targetfile = '/raid60/scratch/uli/testdata/0-20000-0-96035_CEUAS_merged_v0.nc'
data_test = eua.CDMDataset(targetfile)
data_test


# In[10]:


data_test.observations_table


# In[11]:


data_test.recordindices


# In[12]:


data_test.observations_table.observed_variable[0:-1:2000]


# In[13]:

print('read source')
data_test2 = eua.CDMDataset(files[0])
iidx=numpy.where(data_test2.observations_table.observed_variable[:]==85)
iidy=numpy.where(data_test2.observations_table.z_coordinate[:][iidx[0]]==50000)
iidx=iidx[0]
iidy=iidy[0]
iidz=iidx[iidy]
ix=data_test2.observations_table.observation_value[:][iidz]

import matplotlib.pylab as plt
print('read target')
idx=numpy.where(data_test.observations_table.observed_variable[:]==85)
idy=numpy.where(data_test.observations_table.z_coordinate[:][idx[0]]==50000)
idx=idx[0]
idy=idy[0]
idz=idx[idy]
print(idz.shape)
x=data_test.observations_table.observation_value[:][idz]

iri=data_test2.recordindex[:]
ri=data_test.recordindices['85'][:]
plt.plot(x,label='target')
plt.plot(ix,label='source')
plt.show()
keys[-3]


# In[14]:


data_test.recordindices['85'][:50]


# In[15]:


data_test.observations_table.z_coordinate[:36]


# ---

