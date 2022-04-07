#!/usr/bin/env python
# coding: utf-8

# In[5]:


import numpy
import numpy as np
import time
import datetime
import netCDF4
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
from collections import OrderedDict
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

sys.path.append(os.getcwd()+'/../cds-backend/code/')
import cds_eua3 as eua

#import zarr
#import dask
from timeit import default_timer as timer
from numba import njit


# ## Translated hom-functions

# In[2]:


# @njit
def numba_meaneqsamp(test,ref,ni,istart,left_maxlen,right_maxlen,increment,miss_val,max_miss):
    tsa = []
    plus = []
    minus = []
    prms = []
    mrms = []
    pcount = []
    mcount = []
    
    ni = len(t)
    istart=0
    istop=ni
    
    dbin = []
    for i in range(ni):
        dbin.append(i%365.)
        dbin[i]=np.round(dbin[i]/365.*12+1)

    if ref == None:
        diff = test
    else:
        diff = test - ref
    
    ini= True 
    
    for j in [istart]: # istart, istart, increment):
        jstart=j
        if(jstart  <  1):
            jstart=1
            
        jstartlen=j+left_maxlen-jstart
        jend=j+left_maxlen+right_maxlen-1
        if(jend  >  ni) :
            jend=ni
            jendlen=ni-j
        else:
            jendlen=jend-j-left_maxlen+1

        if(jstart  ==  1  or  jend  ==  ni  or  ini) :
            ip= False 
            ipc=0
            im= False 
            imc=0

            ip[:jstartlen]= np.not_equal(diff[jstart:j+left_maxlen-1],  miss_val)
            ipc=numpy.count_nonzero(ip[:jstartlen])
            im[:jendlen]= np.not_equal(diff[j+left_maxlen:jend], miss_val)
            imc=numpy.count_nonzero(im[:jendlen])
        
        minlen=min(left_maxlen,right_maxlen)
        if(minlen-ipc  <  max_miss)  and  (minlen-imc  <  max_miss) :
            for ibin in range(0, 12,1):
                iph[:jstartlen]=np.logical_and(ip[:startlen], dbin[jstart:j+left_maxlen] == ibin)
                pmon[ibin]=numpy.count_nonzero(iph[:jstartlen])
                imh[:jstartlen]=np.logical_and(im[:jendlen], dbin[j+left_maxlen:jend] == ibin)
                mmon[ibin]=numpy.count_nonzero(imh[:jendlen])
                mmonhilf=mmon[ibin]*left_maxlen/right_maxlen
                pmonhilf=pmon[ibin]*right_maxlen/left_maxlen

                if(pmonhilf  <  mmon[ibin]) :
                    i=0
                    k=0

                    while(i  <  (mmon[ibin]-pmonhilf)):
                        while(not imh[jendlen-k] and (k  <  jendlen-1)):
                            k=k+1
                        imh[jendlen-k]= False 
                        im[jendlen-k]= False 
                        i=i+1
                if(mmonhilf  <  pmon[ibin]) :
                    i=0
                    k=0
                    while(i  <  pmon[ibin]-mmonhilf):
                        while( not iph[1+k] and (k  <  jstartlen-1)):
                            k=k+1
                        iph[1+k]= False # remove 1+?
                        ip[1+k]= False # remove 1+?
                        i+=1
            ipc=numpy.count_nonzero(ip[:jstartlen])
            imc=numpy.count_nonzero(im[:jendlen])
            if((minlen-ipc < max_miss) and (inlen-imc < max_miss)):
                for ibin in range(0, 12):
                    pmon[ibin]=np.logical_and(numpy.count_nonzero(ip[:jstartlen]), dbin[jstart:j+left_maxlen-1] == ibin)
                    mmon[ibin]=np.logical_and(numpy.count_nonzero(im[:jendlen]), dbin[j+left_maxlen:jend] == ibin)
                if((jstart == 1) or (jend == ni) or ini) :
                    plusshort[:jstartlen][ip[:jstartlen]] = diff[jstart:j+left_maxlen-1][ip[:jstartlen]]
                    minusshort[:jendlen][im[:jendlen]] = diff[j+left_maxlen:jend][im[:jendlen]]
                    tp=np.sum(plusshort[:jstartlen][ip[:jstartlen]])
                    tm=np.sum(minusshort[:jendlen][im[:jendlen]])
                    tpp=np.sum(plusshort[:jstartlen]*plusshort[:jstartlen][ip[:jstartlen]])
                    tmm=np.sum(minusshort[:jendlen]*minusshort[:jendlen][im[:jendlen]])
                    ini= False 

                qquer=(tp+tm)/(ipc+imc)
                rms=(tpp+tmm)/(ipc+imc)
                sigq=np.sqrt(rms-qquer*qquer)  

                if(sigq  >  0.001) :
                    z1=[tp/ipc-qquer]/sigq
                    z2=[tm/imc-qquer]/sigq
                    tsa=ipc*z1*z1+imc*z2*z2
                    plus=tp/ipc
                    minus=tm/imc
                    prms=tpp/ipc
                    mrms=tmm/imc
                else:
                    ini= True 
        else:
            ini= True 

        ini= True 
        pcount=ipc
        mcount=imc


# In[12]:


@njit
def snhteqsamp(test,ref,ni,istart,istop,maxlen,increment,miss_val,max_miss,tsa,plus,minus,prms,mrms,pcount,mcount):
    
    #dbin = []
    #for i in range(ni):
        #dbin.append(i%365.)
        #dbin[i]=np.round(dbin[i]/365.*12+1)
    
    
    dbin=np.zeros(ni,dtype=np.int32)
    mon=np.array((31,28,31,30,31,30,31,31,30,31,30,31),dtype=np.int32)

    msum=0
    i=0
    while i<ni:
        for m in range(12):
            if msum+mon[m]>ni:
                dbin[msum:ni]=m
            else:
                dbin[msum:msum+mon[m]]=m
            i+=mon[m]
            msum+=mon[m]
        
    
    istart=0
    istop=ni-1 
    if(istart < 0 or istop > ni):
        print('istart,istop ',istart,istop,' must be in range ',1,ni)
        return

    m2=maxlen//2 # is index -> needs to be int
    
    # no need to prepare tsa, gets written, once done
#     if(istart  ==  0  and  istop  ==  ni) :
#         tsa=miss_val
    
    # no need to prepare diff, gets just filled
#     diff=miss_val
#     where(test  !=  miss_val  and  ref  !=  miss_val) 
    if miss_val==miss_val:  # to make sure missing values in diff are encoded as nan
        test[test==miss_val]=np.nan
        ref[ref==miss_val]=np.nan
    diff=test-ref # later add a np.dropna?
    
    
    
    while((diff[istart] != diff[istart]) and (istart < istop)):
        istart=istart+1
        print("skipped missing value on the start")
    while((diff[istop] != diff[istop]) and (istop > istart)):
        istop=istop-1 
        print("skipped missing value on the end")
    if(istart == ni): 
        print('no valid values')
        return
    # no need to prepare output arrays, are and stay empty
#     if(istart  ==  0  and  istop  ==  ni) :
#         plus=miss_val ; minus=miss_val ; prms=miss_val ; mrms=miss_val
#         pcount=0 ; mcount=0
#     else:
#         plus=miss_val ; minus=miss_val ; prms=miss_val ; mrms=miss_val
#         pcount=0 ; mcount=0
        
    ini= True 

    for j in range(istart-max_miss, istop-maxlen+max_miss, increment): # remove +1 of end, bc start is 0


        jstart=j
        if(jstart  <  0): # creates windows with not equal left and right sides
            jstart=0
        jstartlen=j+m2-jstart
        jend=j+maxlen-1
        if(jend  >  ni-1): # =^=
            jend=ni-1
        jendlen=m2+jend-(j+maxlen)

        if((jstart == 0) or (jend == ni-1) or ini) :
#             ip= False 
#             ipc=0

#             im= False 
#             imc=0

#             ip[:jstartlen]=np.not_equal(diff[jstart:j+m2-1], miss_val)

#             print('jstart', jstart)
#             print('jend',jend)
#             print('j',j)
#             print('m2',m2)
            ip=~np.isnan(diff[jstart:j+m2-1])
            ipc=np.sum(ip)
#             im[:jendlen]=np.not_equal(diff[j+m2:jend], miss_val)
            im=~np.isnan(diff[j+m2:jend])
            imc=np.sum(im[:jendlen])
            ipcall=ipc
            imcall=imc


        if((maxlen//2-ipc  <  max_miss)  and  (maxlen//2-imc  <  max_miss)):
            for ibin in range(0, 12,1):
                iph=np.logical_and(ip[:jstartlen], dbin[jstart:j+m2-1] == ibin)
                pmon=numpy.sum(iph) # ??? not sure why it would be pmon(ibin), when only this iteration is used in the whole loop and it needs to get written also once per iteration
                imh=np.logical_and(im[:jendlen], dbin[j+m2:jend] == ibin)
                mmon=numpy.sum(imh) #??? ==^==
                
                if(pmon < mmon) :
                    i=0
                    k=0
                    while(i  <  (mmon-pmon)):
                        while(not imh[jendlen-k-1] and (k  <  jendlen-1)):
                            k=k+1
                        imh[jendlen-k-1]= False 
                        im[jendlen-k-1]= False 
                        i=i+1
                if(mmon  <  pmon) :
                    i=0
                    k=0
                    while(i  <  pmon-mmon):
                        while( not iph[k] and (k  <  jstartlen-1)): # remove 1+? is just an index
                            k=k+1
                        iph[k]= False # ==^==
                        ip[k]= False # ==^==
                        i+=1
            ipc=numpy.sum(ip[:jstartlen])
            imc=numpy.sum(im[:jendlen])

            if((jstart  ==  0)  or  (jend  ==  ni)  or  ini) :
                plusshort=diff[jstart:j+m2-1][ip[:jstartlen]] # correct translation of where?
                minusshort=diff[j+m2:jend][im[:jendlen]] # ==^==
                
                tp=np.sum(plusshort)
                tm=np.sum(minusshort)
                
                tpp=np.sum(plusshort*plusshort)
                tmm=np.sum(minusshort*minusshort)
                ini= False 
    
    
            qquer=(tp+tm)/(ipc+imc)
            rms=(tpp+tmm)/(ipc+imc)
            sigq=np.sqrt(rms-qquer*qquer)  
    
            if(sigq  >  0.001) :
                z1=(tp/ipc-qquer)/sigq
                z2=(tm/imc-qquer)/sigq
                tsa[j+m2:j+m2+increment-1]=ipc*z1*z1+imc*z2*z2
                plus[j+m2:j+m2+increment-1]=tp/ipc
                minus[j+m2:j+m2+increment-1]=tm/imc
                prms[j+m2:j+m2+increment-1]=tpp/ipc
                mrms[j+m2:j+m2+increment-1]=tmm/imc
            # endif
            else:
                ini= True 
            ini= True 
            pcount[j+m2:j+m2+increment-1]=ipcall
            mcount[j+m2:j+m2+increment-1]=imcall

    return
        


# ### dec = nan

# In[14]:


a = 5
nodec = np.concatenate(([1]*334, [np.nan]*31))
allnodec = nodec
for i in range(a-1):
    allnodec = np.concatenate((allnodec, nodec))
    

before_break = np.sin(np.linspace(0,np.pi*2*a,365*a))*1 + np.random.randn(365*a) * 1 + 240
after_break = np.sin(np.linspace(0,np.pi*2*a,365*a))*1 + np.random.randn(365*a) * 1 + 240.0
test = np.concatenate((before_break*allnodec, after_break), axis=0)



# before_break = np.random.rand(365*a) * 15
# after_break = np.random.rand(365*a) * 15 + 2
# test = np.concatenate((before_break, after_break), axis=0)

print(len(test))
fig, ax = maplt.subplots(1, figsize = (15,5))
ax.plot(np.array(range(len(test))), test, color = 'blue', alpha = 0.3, label='test' )
ax.legend()
ax.grid()
maplt.show()
maplt.close()

ni = len(test)

ref = np.zeros(ni)

maxlen = left_maxlen = right_maxlen = 2*365

max_miss = 200 #(or more)
istart = istartorig = 0 
istop = iststoporig = len(test)
increment = 30
miss_val = np.nan #(-999.99)


pcount=np.zeros(ni,dtype=np.int32)
mcount=np.zeros(ni,dtype=np.int32)
tsa = np.zeros(ni,dtype=np.float32)
plus = np.zeros(ni,dtype=np.float32)
minus = np.zeros(ni,dtype=np.float32)
prms = np.zeros(ni,dtype=np.float32)
mrms = np.zeros(ni,dtype=np.float32)

snhteqsamp(np.array(test),np.array(ref),ni,istart,istop,maxlen,increment,miss_val,max_miss,
           tsa,plus,minus,prms,mrms,pcount,mcount)
fig=maplt.figure(figsize = (15,5))
ax = maplt.subplot(2,1,1)
ax.plot(np.array(range(len(tsa))),tsa, color = 'blue', alpha = 0.3, label='tsa' )
ax.legend()
ax = maplt.subplot(2,1,2)
#ax.plot(np.array(range(len(plus))),plus, color = 'red', alpha = 0.3, label='plus' )
#ax.plot(np.array(range(len(minus))),minus, color = 'green', alpha = 0.3, label='minus' )
ax.plot(np.array(range(len(minus))),plus-minus, color = 'green', alpha = 0.3, label='minus' )

ax.legend()
#ax.grid()
maplt.show()
maplt.close()


# ### dec is removed

# In[29]:


a = 5
nodec = np.concatenate(([1]*334, [np.nan]*31))
allnodec = nodec
for i in range(a-1):
    allnodec = np.concatenate((allnodec, nodec))
    

before_break = np.sin(np.linspace(0,np.pi*2*a,365*a))*0 + np.random.rand(365*a) * 15 + 240
after_break = np.sin(np.linspace(0,np.pi*2*a,365*a))*0 + np.random.rand(365*a) * 15 + 245
test = np.concatenate((before_break*allnodec, after_break), axis=0)
test = test[~numpy.isnan(test)]



# before_break = np.random.rand(365*a) * 15
# after_break = np.random.rand(365*a) * 15 + 2
# test = np.concatenate((before_break, after_break), axis=0)

print(len(test))
fig, ax = maplt.subplots(1, figsize = (15,5))
ax.plot(np.array(range(len(test))), test, color = 'blue', alpha = 0.3, label='test' )
ax.legend()
ax.grid()
maplt.show()
maplt.close()

ni = len(test)

ref = [0]*ni

maxlen = left_maxlen = right_maxlen = 2*365

max_miss = 200 #(or more)
istart = istartorig = 0 
istop = iststoporig = len(test)
increment = 30
miss_val = np.nan #(-999.99)


tsa = plus = minus = prms = mrms = pcount = mcount = np.array([0]*ni)

snhteqsamp(np.array(test),np.array(ref),ni,istart,istop,maxlen,increment,miss_val,max_miss,tsa,plus,minus,prms,mrms,pcount,mcount)
fig, ax = maplt.subplots(1, figsize = (15,5))
ax.plot(np.array(range(len(tsa))),tsa, color = 'blue', alpha = 0.3, label='tsa' )
ax.plot(np.array(range(len(plus))),plus, color = 'red', alpha = 0.3, label='plus' )
ax.plot(np.array(range(len(minus))),minus, color = 'green', alpha = 0.3, label='minus' )

ax.legend()
ax.grid()
maplt.show()
maplt.close()


# ### no break, dec is nan

# In[31]:


a = 5
nodec = np.concatenate(([1]*334, [np.nan]*31))
allnodec = nodec
for i in range(a-1):
    allnodec = np.concatenate((allnodec, nodec))
    

before_break = np.sin(np.linspace(0,np.pi*2*a,365*a))*15 + np.random.rand(365*a) * 15 + 240
after_break = np.sin(np.linspace(0,np.pi*2*a,365*a))*15 + np.random.rand(365*a) * 15 + 240
test = np.concatenate((before_break*allnodec, after_break), axis=0)



# before_break = np.random.rand(365*a) * 15
# after_break = np.random.rand(365*a) * 15 + 2
# test = np.concatenate((before_break, after_break), axis=0)

print(len(test))
fig, ax = maplt.subplots(1, figsize = (15,5))
ax.plot(np.array(range(len(test))), test, color = 'blue', alpha = 0.3, label='test' )
ax.legend()
ax.grid()
maplt.show()
maplt.close()

ni = len(test)

ref = [0]*ni

maxlen = left_maxlen = right_maxlen = 2*365

max_miss = 200 #(or more)
istart = istartorig = 0 
istop = iststoporig = len(test)
increment = 30
miss_val = np.nan #(-999.99)


tsa = plus = minus = prms = mrms = pcount = mcount = np.array([0]*ni)

snhteqsamp(np.array(test),np.array(ref),ni,istart,istop,maxlen,increment,miss_val,max_miss,tsa,plus,minus,prms,mrms,pcount,mcount)
fig, ax = maplt.subplots(1, figsize = (15,5))
ax.plot(np.array(range(len(tsa))),tsa, color = 'blue', alpha = 0.3, label='tsa' )
ax.plot(np.array(range(len(plus))),plus, color = 'red', alpha = 0.3, label='plus' )
ax.plot(np.array(range(len(minus))),minus, color = 'green', alpha = 0.3, label='minus' )

ax.legend()
ax.grid()
maplt.show()
maplt.close()


# In[23]:


allnodec = nodec
for i in range(a-1):
    allnodec = np.concatenate((allnodec, nodec))


# In[27]:


len(allnodec)/13


# In[ ]:





# In[ ]:





# In[ ]:





# ## Michael's Function

# In[104]:


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

    nc = numba_snhtmov(np.squeeze(np.asarray(x)),
                       tsa,
                       snhtparas,
                       count,
                       tmean,
                       tsquare)
    return tsa, nc


@njit
def numba_snhtmov(t, tsa, snhtparas, count, tmean, tsquare):
    """Standard Normal Homogeneity Test Moving Window

    t         = np.random.randn(1000)
    snhtparas = np.asarray([100,50,10])
    tsa       = np.zeros(1000)
    tmean     = np.zeros(1000)
    tsquare   = np.zeros(1000)
    count     = np.zeros(1000,dtype=np.int32)

    Output: tsa
    """
    n = snhtparas[0]
    max_miss = snhtparas[1]
    # ninc=snhtparas[2]

    ni = t.shape[0]
    good = 0
    tmean[0] = 0.
    tsquare[0] = 0.
    nancut = []
    for j in range(ni):
        count[j] = 0
        # compare_lists if nan ?
        if t[j] == t[j]:
            if good > 0:
                tmean[good] = tmean[good - 1] + t[j]
                tsquare[good] = tsquare[good - 1] + t[j] * t[j]
            else:
                tmean[good] = t[j]
                tsquare[good] = t[j] * t[j]
            good += 1
            nancut.append(0)
        else:
            nancut.append(1)
            if good > 0:
                tmean[good] = tmean[good - 1] + tmean[good - 2]
                tsquare[good] = tsquare[good - 1] + tsquare[good - 2]
            
        if good > 0:
            count[j] = good - 1

    if good > n - 2 * max_miss:
        rm = int(n / 2)  # needs to be an integer
        # k 1460/2=730 - 650=80, n-80
        for k in range(rm - max_miss, ni - (rm - max_miss)):
            xm = k - rm  # 80-730
            if xm < 0:
                xm = 0
            xp = k + rm
            if xp > ni - 1:
                xp = ni - 1
            if (count[k] - count[xm] > rm - max_miss) and (count[xp] - count[k] > rm - max_miss):
                x = (tmean[count[k]] - tmean[count[xm]]) / (count[k] - count[xm])  # Mittelwert 1 Periode
                y = (tmean[count[xp]] - tmean[count[k]]) / (count[xp] - count[k])  # Mittelwert 2 Periode
                xy = (tmean[count[xp]] - tmean[count[xm]]) / (count[xp] - count[xm])  # Mittelwert ganze Periode

                sig = (tsquare[count[xp]] - tsquare[count[xm]]) / (count[xp] - count[xm])  # t*t ganze Periode
                if sig > xy * xy:
                    sig = np.sqrt(sig - xy * xy)  # standard deviation of the whole window
                    # n1 * (m1-m)**2 + n2 * (m2-m)**2 / stddev
                    tsa[k] = ((count[k] - count[xm]) * (x - xy) * (x - xy) + (count[xp] - count[k]) * (y - xy) * (
                            y - xy)) / (sig * sig)
                else:
                    tsa[k] = 0.
    return nancut


# In[105]:


def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance in kilometers between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1 = numpy.radians(lon1)
    lat1 = numpy.radians(lat1)
    lon2 = numpy.radians(lon2)
    lat2 = numpy.radians(lat2)

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a)) 
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles. Determines return value units.
    return c * r


# In[106]:


stdplevs = [10.0, 20.0, 30.0, 50.0, 70.0, 100.0, 150.0, 200.0, 250.0, 300.0, 400.0, 500.0, 700.0, 850.0, 925.0, 1000.0]
files = glob.glob('/users/staff/leo/fastscratch/rise/1.0/exp03/011035/feedbackmerged*')
df = xr.open_dataset(files[0])

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


# day_snht_700_200 = test(np.array(df_day[dict(pressure=12)].temperatures) - np.array(df_day[dict(pressure=7)].temperatures), window, missing)
# day_snht_200_100 = test(np.array(df_day[dict(pressure=7)].temperatures) - np.array(df_day[dict(pressure=5)].temperatures), window, missing) 
# day_snht_200_50 = test(np.array(df_day[dict(pressure=7)].temperatures) - np.array(df_day[dict(pressure=3)].temperatures), window, missing) 
# day_snht_200_30 = test(np.array(df_day[dict(pressure=7)].temperatures) - np.array(df_day[dict(pressure=2)].temperatures), window, missing) 

# fig, ax = maplt.subplots(1, figsize = (15,5))
# ax.plot(np.array(df_day.datum),day_snht_700_200, color = 'blue', alpha = 0.3, label='700 - 200 hPa', )
# ax.plot(np.array(df_day.datum),day_snht_200_100, color = 'green', alpha = 0.3, label='200 - 100 hPa', )
# ax.plot(np.array(df_day.datum),day_snht_200_50, color = 'red', alpha = 0.3, label='200 - 50 hPa', )
# ax.plot(np.array(df_day.datum),day_snht_200_30, color = 'black', alpha = 0.3, label='200 - 30 hPa', )
# mean_day = (np.array(day_snht_700_200) + np.array(day_snht_200_100) + np.array(day_snht_200_50) + np.array(day_snht_200_30))/4.
# ax.plot(np.array(df_day.datum),mean_day, color = 'black', alpha = 0.8, linewidth = 3, label='MEAN', )
# ax.set_ylabel('SNHT - DAY')
# ax.set_xlabel('time')
# ax.legend()
# ax.grid()
# maplt.show()
# maplt.close()

# night_snht_700_200 = test(np.array(df_night[dict(pressure=12)].temperatures) - np.array(df_night[dict(pressure=7)].temperatures), window, missing)
# night_snht_200_100 = test(np.array(df_night[dict(pressure=7)].temperatures) - np.array(df_night[dict(pressure=5)].temperatures), window, missing) 
# night_snht_200_50 = test(np.array(df_night[dict(pressure=7)].temperatures) - np.array(df_night[dict(pressure=3)].temperatures), window, missing) 
# night_snht_200_30 = test(np.array(df_night[dict(pressure=7)].temperatures) - np.array(df_night[dict(pressure=2)].temperatures), window, missing) 

# fig, ax = maplt.subplots(1, figsize = (15,5))
# ax.plot(np.array(df_night.datum),night_snht_700_200, color = 'blue', alpha = 0.3, label='700 - 200 hPa', )
# ax.plot(np.array(df_night.datum),night_snht_200_100, color = 'green', alpha = 0.3, label='200 - 100 hPa', )
# ax.plot(np.array(df_night.datum),night_snht_200_50, color = 'red', alpha = 0.3, label='200 - 50 hPa', )
# ax.plot(np.array(df_night.datum),night_snht_200_30, color = 'black', alpha = 0.3, label='200 - 30 hPa', )
# mean_night = (np.array(night_snht_700_200) + np.array(night_snht_200_100) + np.array(night_snht_200_50) + np.array(night_snht_200_30))/4.
# ax.plot(np.array(df_night.datum),mean_night, color = 'black', alpha = 0.8, linewidth = 3, label='MEAN', )

# ax.set_ylabel('SNHT - NIGHT')
# ax.set_xlabel('time')
# ax.legend()
# ax.grid()
# maplt.show()
# maplt.close()


# In[107]:


stdplevs = [10.0, 20.0, 30.0, 50.0, 70.0, 100.0, 150.0, 200.0, 250.0, 300.0, 400.0, 500.0, 700.0, 850.0, 925.0, 1000.0]
files = glob.glob('/users/staff/leo/fastscratch/rise/1.0/exp03/011035/feedbackmerged*')
df = xr.open_dataset(files[0])

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

    snht_day,nc= test(np.array(df1.temperatures), window, missing)
    
#     snht_night = test(np.array(df_night_snht.era5_fgdep), window, missing)
#     snht_diff_dn = test(np.array(df_day_snht.era5_fgdep)-np.array(df_night_snht.era5_fgdep), window, missing)
#     snht_diff_dn_t = test(np.array(df_day_snht.temperatures)-np.array(df_night_snht.temperatures), window, missing)
    
    ax[i*2].plot(np.array(df1.datum),snht_day,color = 'blue', alpha = 0.6, label='DAY- ' + str(stdplevs[i]) + ' hPa', )
#     ax[i*2].plot(np.array(df1.datum),np.array(nc)*50,color = 'blue', alpha = 0.2, label='DAY- ' + str(stdplevs[i]) + ' hPa', )
    
    snht_day_corr = snht_day
    snht_day_corr[np.array(nc) == 1] = 0
#     ax[i*2].plot(np.array(df1.datum),snht_day_corr,color = 'red', alpha = 0.6, label='DAY nan - ' + str(stdplevs[i]) + ' hPa', )
    ax[i*2].scatter(np.array(df1.datum),df1.temperatures,color = 'green', alpha = 0.6, label='DAY era5_fgdep - ' + str(stdplevs[i]) + ' hPa', )
#     ax[i].plot(np.array(df_night_snht.datum),snht_night,color = 'orange', alpha = 0.6, label='NIGHT- ' + str(stdplevs[i]) + ' Pa', )
#     ax[i].plot(np.array(df_day_snht.datum),snht_diff_dn,color = 'purple', alpha = 0.6, label='DIFF - ' + str(stdplevs[i]) + ' Pa', )
#     ax[i].plot(np.array(df_day_snht.datum),snht_diff_dn_t,color = 'red', alpha = 0.4, label='DIFF T- ' + str(stdplevs[i]) + ' Pa', )
#     mean_snht = (np.array(snht_day) + np.array(snht_night) + np.array(snht_diff_dn) + np.array(snht_diff_dn_t))/4.
#     ax[i].plot(np.array(df_day_snht.datum),mean_snht,color = 'black', alpha = 0.3, linewidth = 12, label='MEAN - ' + str(stdplevs[i]) + ' Pa', )
    ax[i*2].set_ylabel('SNHT')
    ax[i*2].set_xlabel('time')
    ax[i*2].set_xlim([datetime.date(2008, 1, 1), datetime.date(2013, 12, 31)])
    ax[i*2].legend(loc='center left')
    ax[i*2].grid()
    
    snht_day_o,nc= test(np.array(df_day_snht.temperatures), window, missing)
    ax[i*2+1].plot(np.array(df_day_snht.datum),snht_day_o,color = 'blue', alpha = 0.6, label='DAY- ' + str(stdplevs[i]) + ' hPa', )
    ax[i*2+1].scatter(np.array(df_day_snht.datum),df_day_snht.temperatures,color = 'green', alpha = 0.6, label='DAY era5_fgdep - ' + str(stdplevs[i]) + ' hPa', )
    ax[i*2+1].set_ylabel('SNHT')
    ax[i*2+1].set_xlabel('time')
    ax[i*2+1].set_xlim([datetime.date(2008, 1, 1), datetime.date(2013, 12, 31)])
    ax[i*2+1].legend(loc='center left')
    ax[i*2+1].grid()
    
maplt.show()
maplt.close()


# day_snht_700_200 = test(np.array(df_day[dict(pressure=12)].temperatures) - np.array(df_day[dict(pressure=7)].temperatures), window, missing)
# day_snht_200_100 = test(np.array(df_day[dict(pressure=7)].temperatures) - np.array(df_day[dict(pressure=5)].temperatures), window, missing) 
# day_snht_200_50 = test(np.array(df_day[dict(pressure=7)].temperatures) - np.array(df_day[dict(pressure=3)].temperatures), window, missing) 
# day_snht_200_30 = test(np.array(df_day[dict(pressure=7)].temperatures) - np.array(df_day[dict(pressure=2)].temperatures), window, missing) 

# fig, ax = maplt.subplots(1, figsize = (15,5))
# ax.plot(np.array(df_day.datum),day_snht_700_200, color = 'blue', alpha = 0.3, label='700 - 200 hPa', )
# ax.plot(np.array(df_day.datum),day_snht_200_100, color = 'green', alpha = 0.3, label='200 - 100 hPa', )
# ax.plot(np.array(df_day.datum),day_snht_200_50, color = 'red', alpha = 0.3, label='200 - 50 hPa', )
# ax.plot(np.array(df_day.datum),day_snht_200_30, color = 'black', alpha = 0.3, label='200 - 30 hPa', )
# mean_day = (np.array(day_snht_700_200) + np.array(day_snht_200_100) + np.array(day_snht_200_50) + np.array(day_snht_200_30))/4.
# ax.plot(np.array(df_day.datum),mean_day, color = 'black', alpha = 0.8, linewidth = 3, label='MEAN', )
# ax.set_ylabel('SNHT - DAY')
# ax.set_xlabel('time')
# ax.legend()
# ax.grid()
# maplt.show()
# maplt.close()

# night_snht_700_200 = test(np.array(df_night[dict(pressure=12)].temperatures) - np.array(df_night[dict(pressure=7)].temperatures), window, missing)
# night_snht_200_100 = test(np.array(df_night[dict(pressure=7)].temperatures) - np.array(df_night[dict(pressure=5)].temperatures), window, missing) 
# night_snht_200_50 = test(np.array(df_night[dict(pressure=7)].temperatures) - np.array(df_night[dict(pressure=3)].temperatures), window, missing) 
# night_snht_200_30 = test(np.array(df_night[dict(pressure=7)].temperatures) - np.array(df_night[dict(pressure=2)].temperatures), window, missing) 

# fig, ax = maplt.subplots(1, figsize = (15,5))
# ax.plot(np.array(df_night.datum),night_snht_700_200, color = 'blue', alpha = 0.3, label='700 - 200 hPa', )
# ax.plot(np.array(df_night.datum),night_snht_200_100, color = 'green', alpha = 0.3, label='200 - 100 hPa', )
# ax.plot(np.array(df_night.datum),night_snht_200_50, color = 'red', alpha = 0.3, label='200 - 50 hPa', )
# ax.plot(np.array(df_night.datum),night_snht_200_30, color = 'black', alpha = 0.3, label='200 - 30 hPa', )
# mean_night = (np.array(night_snht_700_200) + np.array(night_snht_200_100) + np.array(night_snht_200_50) + np.array(night_snht_200_30))/4.
# ax.plot(np.array(df_night.datum),mean_night, color = 'black', alpha = 0.8, linewidth = 3, label='MEAN', )

# ax.set_ylabel('SNHT - NIGHT')
# ax.set_xlabel('time')
# ax.legend()
# ax.grid()
# maplt.show()
# maplt.close()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





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




