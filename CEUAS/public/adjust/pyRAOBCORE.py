#!/usr/bin/env python
# coding: utf-8

# In[5]:


import numpy
import numpy as np
import time
import datetime
import h5py
import matplotlib.pylab as plt
import os,sys,glob
sys.path.append(os.getcwd()+'/../adjust/rasotools/')
sys.path.append(os.path.expanduser('~/python/Rasotools/rasotools/'))
sys.path.append(os.path.expanduser('~/python/Rasotools/'))
#sys.path.append(os.getcwd()+'/../adjust/rasotools/')
from utils import *
# import ftest
from multiprocessing import Pool
#import odb
#from eccodes import *
from functools import partial
#from collections import OrderedDict
import subprocess
import json
import gzip
# from retrieve_fb_jra55 import add_feedback
import copy
import pickle
#import xarray as xr

import pandas as pd
import matplotlib
import matplotlib.pylab as plt
import matplotlib.pyplot as maplt
matplotlib.rcParams.update({'font.size': 16})

plt.rcParams['lines.linewidth'] = 3

import warnings
warnings.filterwarnings('ignore')

#sys.path.append(os.getcwd()+'/../cds-backend/code/')
#import cds_eua3 as eua

#import zarr
#import dask
from timeit import default_timer as timer
from numba import njit
import SNHT

# ## Translated hom-functions

# In[2]:


def testeq(x, window, missing):
    """Standard Normal Homogeneity Test (SNHT) with a running window
    Wrapper function for numba_snhtmov

    Args:
        x (np.ndarray) : input data
        window (int) : window size (in days)
        missing (int) : allowed missing values (in days)
    Returns:
        np.ndarray : SNHT
    """

    snhtparas = np.asarray([window, missing, 10])
    tsa = np.zeros_like(x)
    tmean = np.zeros_like(x,shape=(x.shape[0],12))
    tsquare = np.zeros_like(x,shape=(x.shape[0],12))
    count = np.zeros((x.shape[0],12), dtype=np.int32)

    tsa = SNHT.numba_snhteqmov(x,
                       tsa,
                       snhtparas,
                       -999.,
                       count,
                       tmean,
                       tsquare)
    return tsa


def test(x, window, missing):
    """Standard Normal Homogeneity Test (SNHT) with a running window
    Wrapper function for numba_snhtmov

    Args:
        x (np.ndarray) : input data
        window (int) : window size (in days)
        missing (int) : allowed missing values (in days)
    Returns:
        np.ndarray : SNHT
    """

    snhtparas = np.asarray([window, missing, 10])
    tsa = np.zeros(x.shape[0])
    tmean = np.zeros(x.shape[0])
    tsquare = np.zeros(x.shape[0])
    count = np.zeros(x.shape[0], dtype=np.int32)

    tsa = SNHT.numba_snhtmov(np.squeeze(np.asarray(x)),
                       tsa,
                       snhtparas,
                       count,
                       tmean,
                       tsquare)
    return tsa


a = 5
nodec = np.concatenate(([1]*300, [np.nan]*65))
allnodec = nodec
for i in range(a-1):
    allnodec = np.concatenate((allnodec, nodec))
    

before_break = np.cos(np.linspace(0,np.pi*2*a,365*a))*2 + np.random.randn(365*a) * 1 + 240
after_break = np.cos(np.linspace(0,np.pi*2*a,365*a))*2 + np.random.randn(365*a) * 1 + 240.0
ts = np.concatenate((before_break*allnodec, after_break), axis=0)



# before_break = np.random.rand(365*a) * 15
# after_break = np.random.rand(365*a) * 15 + 2
# test = np.concatenate((before_break, after_break), axis=0)

print(len(ts))
fig, ax = maplt.subplots(1, figsize = (15,5))
ax.plot(np.array(range(len(ts))), ts, color = 'blue', alpha = 0.3, label='test' )
ax.legend()
ax.grid()
#maplt.show()
maplt.close()

ni = len(ts)

ref = np.zeros(ni)

maxlen = left_maxlen = right_maxlen = 4*365

max_miss = 650 #(or more)
istart = istartorig = 0 
istop = iststoporig = len(ts)
increment = 30
miss_val = np.nan #(-999.99)


pcount=np.zeros(ni,dtype=np.int32)
mcount=np.zeros(ni,dtype=np.int32)
tsa = np.zeros(ni,dtype=np.float32)
plus = np.zeros(ni,dtype=np.float32)
minus = np.zeros(ni,dtype=np.float32)
prms = np.zeros(ni,dtype=np.float32)
mrms = np.zeros(ni,dtype=np.float32)

SNHT.snhteqsamp(np.array(ts),np.array(ref),ni,istart,istop,maxlen,increment,miss_val,max_miss,
           tsa,plus,minus,prms,mrms,pcount,mcount)

tts=np.empty_like(ts,dtype=np.float32)
tts[:]=ts[:]
tts[np.isnan(tts)]=-999.
tsa2=testeq(tts, maxlen, max_miss)
tsa3=test(ts, maxlen, max_miss)

tt=time.time()
pcount=np.zeros(ni,dtype=np.int32)
mcount=np.zeros(ni,dtype=np.int32)
tsa = np.zeros(ni,dtype=np.float32)
plus = np.zeros(ni,dtype=np.float32)
minus = np.zeros(ni,dtype=np.float32)
prms = np.zeros(ni,dtype=np.float32)
mrms = np.zeros(ni,dtype=np.float32)
SNHT.snhteqsamp(np.array(ts),np.array(ref),ni,istart,istop,maxlen,increment,miss_val,max_miss,
           tsa,plus,minus,prms,mrms,pcount,mcount)
print('classic',time.time()-tt)
tt=time.time()
snhtparas = np.asarray([maxlen, max_miss, 10])
tsa = np.zeros(tts.shape[0],dtype=np.float32)
tmean = np.zeros((tts.shape[0],12),dtype=np.float32)
tsquare = np.zeros((tts.shape[0],12),dtype=np.float32)
count = np.zeros((tts.shape[0],12), dtype=np.int32)
tsa2 = SNHT.numba_snhteqmov(tts, tsa, snhtparas,-999.,count,tmean,tsquare)

#tsa2=testeq(tts, maxlen, max_miss)
print('csumeq',time.time()-tt)
tt=time.time()
tsa3=test(ts, maxlen, max_miss)
print('csum',time.time()-tt)


fig=maplt.figure(figsize = (15,5))
ax = maplt.subplot(2,1,1)
ax.plot(np.array(range(len(tsa))),tsa, color = 'blue', alpha = 0.3, label='tsa' )
ax.plot(np.array(range(len(tsa2))),tsa2, color = 'red', alpha = 0.3, label='tsa2' )
ax.plot(np.array(range(len(tsa3))),tsa3, color = 'green', alpha = 0.3, label='tsa3' )
ax.legend()
ax = maplt.subplot(2,1,2)
#ax.plot(np.array(range(len(plus))),plus, color = 'red', alpha = 0.3, label='plus' )
#ax.plot(np.array(range(len(minus))),minus, color = 'green', alpha = 0.3, label='minus' )
ax.plot(np.array(range(len(minus))),plus-minus, color = 'green', alpha = 0.3, label='minus' )

ax.legend()
#ax.grid()
#maplt.show()
maplt.close()

stdplevs = np.array([10.0, 20.0, 30.0, 50.0, 70.0, 100.0, 150.0, 200.0, 
                     250.0, 300.0, 400.0, 500.0, 700.0, 850.0, 925.0, 1000.0])
pidx=np.array([0,1,2,3,4,5,6,7,8,9,10,11,12,13])

def do_test(fg_dep,days,snhtparas):
    tsa = np.zeros_like(fg_dep,shape=(fg_dep.shape[2]))
    tmean = np.zeros_like(fg_dep,shape=(fg_dep.shape[2],12))
    tsquare = np.zeros_like(fg_dep,shape=(fg_dep.shape[2],12))
    count = np.zeros(shape=(fg_dep.shape[2],12),dtype=np.int32)
    tsas=np.zeros_like(fg_dep)
    tsarad=np.zeros_like(fg_dep,shape=(1,5,fg_dep.shape[2]))
   
    #fg_dep[np.isnan(fg_dep)]=-999.
    for ih in range(fg_dep.shape[0]):
        for ip in range(fg_dep.shape[1]):
            tsas[ih,ip,:] = SNHT.numba_snhteqmov(fg_dep[ih,ip,:].flatten(), tsa, snhtparas[:3],snhtparas[3],
                                                 count,tmean,tsquare)
            #print(ih,ip)
            
    for ip in range(5):
        tsarad[0,ip,:] = SNHT.numba_snhteqmov((fg_dep[1,ip,:]-fg_dep[0,ip,:]).flatten(), tsa, snhtparas[:3],snhtparas[3],
                                                     count,tmean,tsquare)
        #print(ip)
    
    return tsas,tsarad


def RAOB_findbreaks(method,fn):
    with h5py.File(fn,'r') as f:
        fg_dep=f['era5_fgdep'][:,pidx,:]
        days=f['datum'][:]
        print(fn.split('/')[-1],days.shape)
    
    snhtparas = (maxlen, max_miss, 10,-999.)
    tsas,tsarad=do_test(fg_dep,days,snhtparas)
    
    breaks=tsas #np.array([0])
    
    return breaks
    
files = glob.glob('/users/staff/leo/fastscratch/rise/1.0/exp03/*/feedbackmerged*')
func=partial(RAOB_findbreaks,'SNHT')
tt=time.time()
breaks=list(map(func,files))
print(time.time()-tt)

window = 1460  # means 4 years on daily basis
missing = 600

df_night = df[dict(hour=0)]
df_day = df[dict(hour=1)]

fig, ax = maplt.subplots(len(stdplevs)*2, 1, figsize = (15,80))
for i in range(len(stdplevs)):
    df_day_snht = df_day[dict(pressure=i)]
    df_night_snht = df_night[dict(pressure=i)]
    df1 = df_day_snht.to_dataframe()
    
    for yr in [2010, 2011]:
        for mn in [6, 7, 8]: # range(1,12,1): #
            df1 = df1[df1.datum < str(yr)+"-"+str(mn)].append(df1[df1.datum >= str(yr)+"-"+str(mn+1)])

#     for yr in [2015, 2016, 2017, 2018]:
#         df1 = df1[df1.datum < str(yr)].append(df1[df1.datum >= str(1+yr)])

    snht_day,nc= test(np.array(df1.era5_fgdep), window, missing)
    
#     snht_night = test(np.array(df_night_snht.era5_fgdep), window, missing)
#     snht_diff_dn = test(np.array(df_day_snht.era5_fgdep)-np.array(df_night_snht.era5_fgdep), window, missing)
#     snht_diff_dn_t = test(np.array(df_day_snht.temperatures)-np.array(df_night_snht.temperatures), window, missing)
    
    ax[i*2].plot(np.array(df1.datum),snht_day,color = 'blue', alpha = 0.6, label='DAY- ' + str(stdplevs[i]) + ' hPa', )
#     ax[i*2].plot(np.array(df1.datum),np.array(nc)*50,color = 'blue', alpha = 0.2, label='DAY- ' + str(stdplevs[i]) + ' hPa', )
    
    snht_day_corr = snht_day
    snht_day_corr[np.array(nc) == 1] = 0
    ax[i*2].plot(np.array(df1.datum),snht_day_corr,color = 'red', alpha = 0.6, label='DAY nan - ' + str(stdplevs[i]) + ' hPa', )
    ax[i*2].scatter(np.array(df1.datum),df1.era5_fgdep,color = 'green', alpha = 0.6, label='DAY era5_fgdep - ' + str(stdplevs[i]) + ' hPa', )
#     ax[i].plot(np.array(df_night_snht.datum),snht_night,color = 'orange', alpha = 0.6, label='NIGHT- ' + str(stdplevs[i]) + ' Pa', )
#     ax[i].plot(np.array(df_day_snht.datum),snht_diff_dn,color = 'purple', alpha = 0.6, label='DIFF - ' + str(stdplevs[i]) + ' Pa', )
#     ax[i].plot(np.array(df_day_snht.datum),snht_diff_dn_t,color = 'red', alpha = 0.4, label='DIFF T- ' + str(stdplevs[i]) + ' Pa', )
#     mean_snht = (np.array(snht_day) + np.array(snht_night) + np.array(snht_diff_dn) + np.array(snht_diff_dn_t))/4.
#     ax[i].plot(np.array(df_day_snht.datum),mean_snht,color = 'black', alpha = 0.3, linewidth = 12, label='MEAN - ' + str(stdplevs[i]) + ' Pa', )
    ax[i*2].set_ylabel('SNHT')
    ax[i*2].set_xlabel('time')
    ax[i*2].legend(loc='center left')
    ax[i*2].grid()
    
    snht_day_o,nc= test(np.array(df_day_snht.era5_fgdep), window, missing)
    ax[i*2+1].plot(np.array(df_day_snht.datum),snht_day_o,color = 'blue', alpha = 0.6, label='DAY- ' + str(stdplevs[i]) + ' hPa', )
    ax[i*2+1].set_ylabel('SNHT')
    ax[i*2+1].set_xlabel('time')
    ax[i*2+1].legend(loc='center left')
    ax[i*2+1].grid()
    
maplt.show()
maplt.close()

files = glob.glob('/users/staff/leo/fastscratch/rise/1.0/exp03/011035/feedbackmerged*')
func=partial(RAOB_findbreaks,'SNHT')
breaks=list(map(func,files))


exit()
# ## ACMANT4 Input

# In[ ]:


# creating input files for ACMANT4

testdf = df_day.to_dataframe()
snumb = 1

from datetime import date, timedelta

def daterange(start_date, end_date):
    for n in range(int((end_date - start_date).days)):
        yield start_date + timedelta(n)
        
def last_day_of_month(any_day):
    next_month = any_day.replace(day=28) + datetime.timedelta(days=4)
    return next_month - datetime.timedelta(days=next_month.day)

for i in testdf.press.drop_duplicates():
    pdf = testdf[testdf.press == i]
    lines = []
    lines.append(str(df.unique_source_identifier) + '\n')
    for yr in range(1950,2022,1):
        print(yr)
        for mn in range(1,13,1):
            line = str(yr) + '\t' + str(mn)
            start_date = date(yr, mn, 1)
            end_date = last_day_of_month(start_date)
            for single_date in daterange(start_date, end_date + datetime.timedelta(days=1)):
                try:
                    temp = float(pdf[pdf.datum == single_date.strftime("%Y-%m-%d")].temperatures)
                    if temp == temp:
                        line = line + '\t' + str(temp-273.15)[:6]
                    else:
                        line = line + '\t' + '-999.9'
                except:
                    line = line + '\t' + '-999.9'
            lines.append(line + '\n')
    with open('S'+str(snumb).rjust(4, '0')+str(int(i)).rjust(4, '0')+"t.txt","w+") as f:
        f.writelines(lines)


# In[ ]:


with open("text.txt","w+") as f:
    f.writelines(lines)


# In[13]:


str(int(50.0)).rjust(5, '0')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[49]:


# ftest.snhteqsamp2y(
#     test= ,
#     ref= ,
#     ni= ,
#     istartorig= 0, #gets overwritten in the function
#     istoporig= 0, #gets overwritten in the function
#     maxlen= ,
#     increment= ,
#     miss_val= np.nan,
#     max_miss= ,
#     critical_dates= ,
#     ncritical= ,
#     tsa=, 
#     plus= ,
#     minus= ,
#     prms= ,
#     mrms= ,
#     pcount= ,
#     mcount= 
#     )


# In[179]:


# files = glob.glob('/users/staff/leo/fastscratch/rise/1.0/exp03/*/feedbackmerged*')
# latlon = {}
# for i in files:
#     df = xr.open_dataset(i)
#     latlon[i]=[float(df.lat[0]), float(df.lon[0])]
# with open("stat_latlon.p", "wb") as output_file:
#     pickle.dump(latlon, output_file)
    
# with open('stat_latlon.p', "rb") as input_file:
#     latlon = pickle.load(input_file)
# distdict = {}
# for i in latlon:
#     lat = latlon[i][0]
#     lon = latlon[i][1]
#     dists = []
#     for j in latlon:
#         dists.append([j, haversine(lat1=lat, lon1=lon, lat2=latlon[j][0], lon2=latlon[j][1])])
#     distdict[i]=dists
    
# ct = 0
# out_to_file = {}
# for i in distdict:
#     print(i)
#     save_to_dict = []
#     stats, dists = np.array(distdict[i]).transpose()
#     dists = np.array(dists).astype(float)
#     out = [x for x in sorted(zip(dists, stats))]
#     for j in range(len(out))[:10]: # only save the 10 nearest stations includin self
#         save_to_dict.append([out[j][1], out[j][0], latlon[out[j][1]]])
#     out_to_file[i] = save_to_dict
# #     ct+=1
# #     if ct >5:break
# with open("nearest_stations.p", "wb") as output_file:
#     pickle.dump(out_to_file, output_file)


# In[182]:


# files = glob.glob('/users/staff/leo/fastscratch/rise/1.0/exp03/*/feedbackmerged*')
# window = 1460  # means 4 years on daily basis
# latlon = {}
# for i in files:
#     df = xr.open_dataset(i)
#     if len(df.time) >= window:
#         latlon[i]=[float(df.lat[0]), float(df.lon[0])]
# # with open("stat_latlon.p", "wb") as output_file:
# #     pickle.dump(latlon, output_file)
    
# # with open('stat_latlon.p', "rb") as input_file:
# #     latlon = pickle.load(input_file)
# distdict = {}
# for i in latlon:
#     lat = latlon[i][0]
#     lon = latlon[i][1]
#     dists = []
#     for j in latlon:
#         dists.append([j, haversine(lat1=lat, lon1=lon, lat2=latlon[j][0], lon2=latlon[j][1])])
#     distdict[i]=dists
    
# ct = 0
# out_to_file = {}
# for i in distdict:
#     print(i)
#     save_to_dict = []
#     stats, dists = np.array(distdict[i]).transpose()
#     dists = np.array(dists).astype(float)
#     out = [x for x in sorted(zip(dists, stats))]
#     for j in range(len(out))[:10]: # only save the 10 nearest stations includin self
#         save_to_dict.append([out[j][1], out[j][0], latlon[out[j][1]]])
#     out_to_file[i] = save_to_dict
# #     ct+=1
# #     if ct >5:break
# with open("nearest_stations_with_snht.p", "wb") as output_file:
#     pickle.dump(out_to_file, output_file)


# In[183]:


with open('nearest_stations.p', "rb") as input_file:
    nearest_stations = pickle.load(input_file)


# In[184]:


with open('nearest_stations_with_snht.p', "rb") as input_file:
    nearest_stations_with_snht = pickle.load(input_file)


# In[186]:


vie = glob.glob('/users/staff/leo/fastscratch/rise/1.0/exp03/*11035*/feedbackmerged*')[0]
print(vie)
print()
display(nearest_stations[vie])
print()
display(nearest_stations_with_snht[vie])


# In[ ]:





# In[ ]:





# In[ ]:





# In[8]:


stdplevs = [10.0, 20.0, 30.0, 50.0, 70.0, 100.0, 150.0, 200.0, 250.0, 300.0, 400.0, 500.0, 700.0, 850.0, 925.0, 1000.0]
window = 1460  # means 4 years on daily basis
missing = 600

LOG_EVERY_N = 50

files = glob.glob('/users/staff/leo/fastscratch/rise/1.0/exp03/*/feedbackmerged*')
night_save = {}
night_times = {}
day_save = {}
day_times = {}
for i in stdplevs:
    night_save[i]=[]
    night_times[i]=[]
    day_save[i]=[]
    day_times[i]=[]
t0 = time.time()
for j in range(len(files)):
    if (j % LOG_EVERY_N) == 0: print(j)
    try:
        df = xr.open_dataset(files[j])
        for i in range(len(stdplevs)):
            df_new = df[dict(pressure=i)]
            if len(df_new.time) >= window:
                df_snht = df_new[dict(hour=0)]
                snht = test(np.array(df_snht.temperatures), window, missing)
                night_save[stdplevs[i]].append(np.array(snht))
                night_times[stdplevs[i]].append(np.array(df_snht.datum))

                df_snht = df_new[dict(hour=1)]
                snht = test(np.array(df_snht.temperatures), window, missing)
                day_save[stdplevs[i]].append(np.array(snht))
                day_times[stdplevs[i]].append(np.array(df_snht.datum))
    except:
        print('xxxxxxxxx', files[j], j)
        break
print(time.time()-t0)


# In[41]:


stdplevs = [10.0, 20.0, 30.0, 50.0, 70.0, 100.0, 150.0, 200.0, 250.0, 300.0, 400.0, 500.0, 700.0, 850.0, 925.0, 1000.0]

nlen = 0
max_time = 0
for i in day_times[500]:
    if len(i) > nlen:
        max_time = i

for j in stdplevs:
    mean_snht = []
    print(j)
#     counter = 0
#     print(len(max_time))
    ic = []
    for i in max_time:
#         if (counter % 100) == 0: print(counter)
#         counter += 1
        summerizer = []
        inputcounter = 0
        for h in range(len(day_save[j][:100])):
            if i in day_times[j][h]:
                summerizer.append(day_save[j][h][day_times[j][h] == i][0])
                inputcounter += 1
        if divider != 0:
            mean_snht.append(np.nanmean(summerizer))
        else:
            mean_snht.append(0)
        ic.append(inputcounter)
                
    fig, ax = maplt.subplots(1, figsize = (15,5))
    ax.plot(np.array(max_time),np.array(mean_snht),color = 'blue', alpha = 0.6, label='SNHT - ' + str(j) + ' hPa', )
    ax.plot(np.array(max_time),np.array(ic),color = 'green', alpha = 0.6, label='number of input stations', )
    ax.set_ylabel('SNHT')
    ax.set_xlabel('time')
    ax.legend(loc='upper right')
    ax.grid()
    maplt.show()
    maplt.close()


# In[10]:


len(day_save[10])


# In[105]:


# window = 1460  # means 4 years on daily basis
# missing = 600

# files = glob.glob('/users/staff/leo/fastscratch/rise/1.0/exp03/011035/feedbackmerged*')
# df = xr.open_dataset(files[0]).to_dataframe()
# plevs = np.array(df.press.drop_duplicates())
# fig, ax = maplt.subplots(len(plevs), 1, figsize = (15,40))
# for i in range(len(plevs)):
#     df_new = df[df.press == plevs[i]]
#     resort = df_new.sort_values([('time')], ascending=True)
#     snht = test(np.array(resort.temperatures), window, missing)
    
#     ax[i].plot(np.array(resort.datum),snht,color = 'blue', label='SNHT - ' + str(plevs[i]) + ' Pa', )
#     ax[i].set_ylabel('SNHT')
#     ax[i].set_xlabel('time')
#     ax[i].legend(loc='upper right')
#     ax[i].grid()

# # maplt.title('SNHT')
# maplt.show()
# maplt.close()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




