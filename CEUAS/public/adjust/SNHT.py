#!/usr/bin/env python
# coding: utf-8

# In[5]:

import matplotlib.pylab as plt
import numpy as np
import time

from numba import njit


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

    tsa = numba_snhteqmov(x,
                       tsa,
                       snhtparas,
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

    tsa = numba_snhtmov(np.squeeze(np.asarray(x)),
                       tsa,
                       snhtparas,
                       count,
                       tmean,
                       tsquare)
    return tsa

#@njit(fastmath={'nsz','arcp','contract','afn','reassoc'},cache=True)
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
    min_sampsize = snhtparas[1]
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

    if good > min_sampsize:
        rm = int(n / 2)  # needs to be an integer
        # k 1460/2=730 - 650=80, n-80
        for k in range(min_sampsize, ni - min_sampsize):
            xm = k - rm  # 80-730
            if xm < 0:
                xm = 0
            xp = k + rm
            if xp > ni - 1:
                xp = ni - 1
            if (count[k] - count[xm] > min_sampsize) and (count[xp] - count[k] > min_sampsize):
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
    return tsa

numba_snhtmov_njit =njit(numba_snhtmov,fastmath={'nsz','arcp','contract','afn','reassoc'},debug=True,cache=False)

# Michis function with equal sampling
# t is time series
# dbin=RC['months'][days since reftime] number of month (0-11)  - should be indexed with days since reftime !!
# tsa=test statistic (array with lenght of t, initialized in calling function)
# snhtparas: 0=window length (default 2x2 years), 1=allowed missing values in one window half (default 650), 2 not used
#            3= missing value - NaN does not work with fastmath=True
# count (integer),tmean (float),tsquare (float are auxiliary arrays) initialized in calling function (to avoid reallocation)
# kref - if nonzero, the SNHT is calculated only at time series position kref 
#     (typically used to actually calculate the means of the intervals prior and after kref)

@njit(fastmath={'nsz','arcp','contract','afn','reassoc'},cache=True)
def numba_snhteqmov(t, dbin, tsa, n, min_sampsize,miss_val, count, tmean, tsquare,kref=0):
#def numba_snhtmov(t, tsa, snhtparas, count, tmean, tsquare):
    """Standard Normal Homogeneity Test Moving Window

    t         = np.random.randn(1000)
    snhtparas = np.asarray([100,50,10])
    tsa       = np.zeros(1000)
    tmean     = np.zeros(1000)
    tsquare   = np.zeros(1000)
    count     = np.zeros(1000,dtype=np.int32)

    Output: tsa
    """

    ni = t.shape[0]
    #dbin=month[days]
    #mon=np.array((31,28,31,30,31,30,31,31,30,31,30,31),dtype=np.int32)

    #msum=0
    #i=0
    #while i<ni:
        #for m in range(12):
            #if msum+mon[m]>ni:
                #dbin[msum:ni]=m
            #else:
                #dbin[msum:msum+mon[m]]=m
            #i+=mon[m]
            #msum+=mon[m]

    tmean[0,:] = 0.
    tsquare[0,:] = 0.
    count[0,:]=0
    #print(t[0],miss_val,np.isnan(t[0]),np.isnan(miss_val))
    if ~np.isnan(t[0]):
        tmean[0,0]=t[0]
        tsquare[0,0]=t[0]
        count[0,0]=1

    for j in range(1,ni):
        # compare_lists if nan ?
        tmean[j,:]=tmean[j-1,:]
        tsquare[j,:]=tsquare[j-1,:]
        count[j,:]=count[j-1,:]
        dj=dbin[j]
#        if t[j] == t[j]:  
        if ~np.isnan(t[j]):  
            tmean[j,dj]+= t[j]
            tsquare[j,dj]+= t[j]*t[j]
            count[j,dj]+=1
    
    bmax=tsa[0]
    tsamax=tsa[0]
    kmax=0
    if np.sum(count[-1,:]) > min_sampsize:
        rm = n//2  # needs to be an integer
        # k 1460/2=730 - 650=80, n-80
        xm=0
        k=min_sampsize-1
        xp=k+rm
        if xp>ni-1:
            xp=ni-1 

        spsave=np.sum(count[xp,:] - count[k,:])
        smsave=np.sum(count[k,:] - count[xm,:])
        k=min_sampsize
        ll=0
        while k<ni - min_sampsize:
            xm = k - rm  # 80-730
            if xm < 0:
                xm = 0
            xp = k + rm
            if xp > ni - 1:
                xp = ni - 1
            ck=count[k,dbin[k]]-count[k-1,dbin[k]]
            if xm>0:
                smsave += ck-(count[xm,dbin[xm]]-count[xm-1,dbin[xm]]) #if True: #(np.sum(count[k,:] - count[xm,:]) > rm - max_miss) and (np.sum(count[xp,:] - count[k,:]) > rm - max_miss):
            else:
                smsave += ck #if True: #(np.sum(count[k,:] - count[xm,:]) > rm - max_miss) and (np.sum(count[xp,:] - count[k,:]) > rm - max_miss):
            spsave += (count[xp,dbin[xp]]-count[xp-1,dbin[xp]])-ck #if True: #(np.sum(count[k,:] - count[xm,:]) > rm - max_miss) and (np.sum(count[xp,:] - count[k,:]) > rm - max_miss):
            #print(k,smsave,spsave)
            # skip tsa calculation if previous tsa value was small
            #if ( tsa[k-1]<10 or tsa[k-11]>tsa[k-1] ) and ll<10 and tsa[k-1]!=0:
                #tsa[k]=tsa[k-1]
                #k+=1
                #ll+=1
                #continue
            ll=0
            if smsave > min_sampsize and spsave > min_sampsize:
                cp=0; cm=0; c=0; x=np.float32(0); y=x; xx=x; xy=x
                sm=smsave
                sp=spsave
                mcbin=0
                for ibin in range(12):
                    kp=k
                    km=k
                    cdiff=count[k,ibin] - count[xm,ibin]-(count[xp,ibin] - count[k,ibin])
                    #r=cdiff - count[k,ibin]
                    while(kp<xp and cdiff<0 and sp>min_sampsize):
                        dbdiff=(ibin-dbin[kp]+12)%12
                        if dbdiff>1:
                            kp+=(dbdiff-1)*30
                            if kp>=xp:
                                kp=xp-1
                                break
                        kp+=1
                        while kp<xp and count[kp,ibin]>count[kp-1,ibin] and sp>min_sampsize:
                            cdiff+=1 #=r+ count[kp,ibin]
                            kp+=1
                            #print('p',k,kp,ibin,dbin[kp],cdiff)
                            sp-=1
                    #r=count[xm,ibin]+(count[xp,ibin] - count[k,ibin])
                    while(km>xm and cdiff>0 and sm>min_sampsize):
                        dbdiff=(dbin[km]-ibin+12)%12
                        if dbdiff>1:
                            km-=(dbdiff-1)*30
                            if km<xm:
                                km=xm
                                break
                        km-=1
                        while count[km+1,ibin]>count[km,ibin] and sm>min_sampsize:
                            cdiff-=1 #=count[km,ibin] - r
                            #print('m',k,km,ibin,dbin[km],cdiff)
                            km-=1   
                            sm-=1
                    
                    if count[km,ibin] - count[xm,ibin]>0 and count[xp,ibin] - count[kp,ibin]>0:
                        
                        x+= (tmean[km,ibin] - tmean[xm,ibin]) / (count[km,ibin] - count[xm,ibin])  # Mittelwert 1 Periode
                        mcbin+=1
                        
                        y+= (tmean[xp,ibin] - tmean[kp,ibin]) / (count[xp,ibin] - count[kp,ibin])  # Mittelwert 2 Periode
                        xy+= (tmean[xp,ibin] - tmean[xm,ibin]-(tmean[kp,ibin]-tmean[km,ibin])) / (count[xp,ibin] - count[xm,ibin]-(count[kp,ibin]-count[km,ibin]))  # Mittelwert ganze Periode
                        xx+= (tsquare[xp,ibin] - tsquare[xm,ibin]-(tsquare[kp,ibin]-tsquare[km,ibin])) / (count[xp,ibin] - count[xm,ibin]-(count[kp,ibin]-count[km,ibin]))  # Mittelwert ganze Periode
                        
                    #c+= (count[xp,ibin] - count[xm,ibin]-(count[kp,ibin]-count[km,ibin]))
                    cp+=count[xp,ibin] - count[kp,ibin]
                    cm+=count[km,ibin] - count[xm,ibin]
                    #xy+= (tmean[xp,ibin] - tmean[xm,ibin]-(tmean[kp,ibin]-tmean[km,ibin])) #/ (count[xp,ibin] - count[xm,ibin]-(count[kp,ibin]-count[km,ibin]))  # Mittelwert ganze Periode
                    #xx+= (tsquare[xp,ibin] - tsquare[xm,ibin]-(tsquare[kp,ibin]-tsquare[km,ibin])) #/ (count[xp,ibin] - count[xm,ibin]-(count[kp,ibin]-count[km,ibin]))  # Mittelwert ganze Periode
                
                if cm+1>=min_sampsize and cp+1>=min_sampsize and mcbin>3:
                    
                    sig = xx / np.float32(mcbin)  # t*t ganze Periode
                    xym=xy/ np.float32(mcbin)
                    if sig > xym * xym:
                        sig = sig - xym * xym  # standard deviation of the whole window
                        # n1 * (m1-m)**2 + n2 * (m2-m)**2 / stddev
                        tsa[k] = 30.*(np.float32(mcbin)* (x/np.float32(mcbin) - xym)*(x/np.float32(mcbin) - xym)  + np.float32(mcbin) * (y/np.float32(mcbin) - xym)*(y/np.float32(mcbin) - xym)) / sig
                        if tsa[k]>tsamax:
                            bmax=y/np.float32(mcbin)-x/np.float32(mcbin)
                            tsamax=tsa[k]
                            kmax=k
                    else:
                        tsa[k] = 0.
                else:
                    #print(k,cm,cp,'too small')
                    tsa[k]=0.
            else:
                tsa[k]=0.
            k+=1
            
    if kref:
        return y/np.float32(mcbin)-x/np.float32(mcbin)  # break estimate at kref
    else:
        
        return bmax # break estimate at position with maximum tsa

# Averaging with equal sampling
@njit(fastmath={'nsz','arcp','contract','afn','reassoc'},cache=True)
def numba_meaneqmov(t, dbin,tsa, n,min_sampsize, miss_val, count, tmean, tsquare,kref=0):

    ni = t.shape[0]

    tmean[0,:] = 0.
    tsquare[0,:] = 0.
    count[0,:]=0
    if ~np.isnan(t[0]):
        tmean[0,0]=t[0]
        tsquare[0,0]=t[0]
        count[0,0]=1

    for j in range(1,ni):
        # compare_lists if nan ?
        tmean[j,:]=tmean[j-1,:]
        tsquare[j,:]=tsquare[j-1,:]
        count[j,:]=count[j-1,:]
        dj=dbin[j]
#        if t[j] == t[j]:  
        if ~np.isnan(t[j]) :  
            tmean[j,dj]+= t[j]
            tsquare[j,dj]+= t[j]*t[j]
            count[j,dj]+=1
    
    bmax=miss_val
    tsamax=tsa[0]
    kmax=0
    
    if np.sum(count[ni-1,:]) > min_sampsize:
        rm = n//4  # needs to be an integer
        # k 1460/2=730 - 650=80, n-80
        xm=0
        k=min_sampsize-1
        xp=k+rm
        if xp>ni-1:
            xp=ni-1 

        spsave=np.sum(count[xp,:] - count[k,:])
        smsave=np.sum(count[k,:] - count[xm,:])
        if kref:
            k=kref
            #gen=range(kref,kref+1)
        #else:
            #gen=range(rm - max_miss, ni - (rm - max_miss))
        #for k in gen:
            xm = k - rm  # 80-730
            if xm < 0:
                xm = 0
            xp = k + rm
            if xp > ni - 1:
                xp = ni - 1

            smsave=np.sum(count[k,:] - count[xm,:])
            spsave=np.sum(count[xp,:] - count[k,:]) 
            #print(k,smsave,spsave)
            if smsave > min_sampsize and spsave > min_sampsize:
                cp=0; cm=0; c=0; x=np.float32(0); y=x; xx=x; xy=x
                sm=smsave
                sp=spsave
                pcbin=0
                mcbin=0
                for ibin in range(12):
                    kp=k
                    km=k
                    cdiff=count[k,ibin] - count[xm,ibin]-(count[xp,ibin] - count[k,ibin])
                    #r=cdiff - count[k,ibin]
                    while(kp<xp and cdiff<0 and sp>min_sampsize):
                        dbdiff=(ibin-dbin[kp]+12)%12
                        if dbdiff>1:
                            kp+=(dbdiff-1)*30
                            if kp>xp-1:
                                kp=xp-1
                                break
                        kp+=1
                        while kp<xp and count[kp,ibin]>count[kp-1,ibin] and sp>min_sampsize:
                            cdiff+=1 #=r+ count[kp,ibin]
                            kp+=1
                            #print('p',k,kp,ibin,dbin[kp],cdiff)
                            sp-=1
                    #r=count[xm,ibin]+(count[xp,ibin] - count[k,ibin])
                    while(km>xm and cdiff>0 and sm>min_sampsize):
                        dbdiff=(dbin[km]-ibin+12)%12
                        if dbdiff>1:
                            km-=(dbdiff-1)*30
                            if km<xm:
                                km=xm+1
                                break
                        km-=1
                        while km>0 and count[km+1,ibin]>count[km,ibin] and sm>min_sampsize:
                            cdiff-=1 #=count[km,ibin] - r
                            #print('m',k,km,ibin,dbin[km],cdiff)
                            km-=1
                            sm-=1

                    if (count[km,ibin] - count[xm,ibin]>5) and (count[xp,ibin] - count[kp,ibin])>5:
                        
                        x+= (tmean[km,ibin] - tmean[xm,ibin]) / (count[km,ibin] - count[xm,ibin])  # Mittelwert 1 Periode
                        y+= (tmean[xp,ibin] - tmean[kp,ibin]) / (count[xp,ibin] - count[kp,ibin])  # Mittelwert 2 Periode
                        mcbin+=1
                        
                    #c+= (count[xp,ibin] - count[xm,ibin]+(count[kp,ibin]-count[km,ibin]))
                    cp+=count[xp,ibin] - count[kp,ibin]
                    cm+=count[km,ibin] - count[xm,ibin]
                    #xy+= (tmean[xp,ibin] - tmean[xm,ibin]-(tmean[kp,ibin]-tmean[km,ibin])) #/ (count[xp,ibin] - count[xm,ibin]-(count[kp,ibin]-count[km,ibin]))  # Mittelwert ganze Periode
                    #xx+= (tsquare[xp,ibin] - tsquare[xm,ibin]-(tsquare[kp,ibin]-tsquare[km,ibin])) #/ (count[xp,ibin] - count[xm,ibin]-(count[kp,ibin]-count[km,ibin]))  # Mittelwert ganze Periode
                
                if cm>=min_sampsize and cp>=min_sampsize and mcbin>3:
                    
                    bmax=y/np.float32(mcbin)-x/np.float32(mcbin)
                else:
                    #print(k,cm,cp,'too small')
                    tsa[k]=0.
            else:
                tsa[k]=0.
    #if np.abs(bmax)>5:
        #print(bmax)
    return bmax
    
@njit(cache=True,fastmath={'nsz','arcp','contract','afn','reassoc'})
def meqm_ini(t,tmean,tsquare,count,dbin):
    
    ni = t.shape[0]
    tmean[0,:] = 0.
    tsquare[0,:] = 0.
    count[0,:]=0
    if ~np.isnan(t[0]):
        tmean[0,0]=t[0]
        tsquare[0,0]=t[0]
        count[0,0]=1

    for j in range(1,ni):
        # compare_lists if nan ?
        tmean[j,:]=tmean[j-1,:]
        tsquare[j,:]=tsquare[j-1,:]
        count[j,:]=count[j-1,:]
        dj=dbin[j]
        if ~np.isnan(t[j]) :  
            tmean[j,dj]+= t[j]
            tsquare[j,dj]+= t[j]*t[j]
            count[j,dj]+=1
    return
    
def meaneqmov(t, dbin,tsa, n,min_sampsize, miss_val, count, tmean, tsquare,kref=0):

    ni = t.shape[0]

    meqm_ini(t,tmean,tsquare,count,dbin)
    
    bmax=miss_val
    tsamax=tsa[0]
    kmax=0
    
    if np.sum(count[ni-1,:]) > min_sampsize:
        rm = n//4  # needs to be an integer
        # k 1460/2=730 - 650=80, n-80
        xm=0
        k=min_sampsize-1
        xp=k+rm
        if xp>ni-1:
            xp=ni-1 

        spsave=np.sum(count[xp,:] - count[k,:])
        smsave=np.sum(count[k,:] - count[xm,:])
        if kref:
            gen=range(kref,kref+1)
        else:
            gen=range(min_sampsize, min_sampsize)
        for k in gen:
            xm = k - rm  # 80-730
            if xm < 0:
                xm = 0
            xp = k + rm
            if xp > ni - 1:
                xp = ni - 1

            smsave=np.sum(count[k,:] - count[xm,:])
            spsave=np.sum(count[xp,:] - count[k,:]) 
            #print(k,smsave,spsave)
            if smsave > min_sampsize and spsave > min_sampsize:
                cp=0; cm=0; c=0; x=np.float32(0); y=x; xx=x; xy=x
                sm=smsave
                sp=spsave
                for ibin in range(12):
                    kp=k
                    km=k
                    cdiff=count[k,ibin] - count[xm,ibin]-(count[xp,ibin] - count[k,ibin])
                    #r=cdiff - count[k,ibin]
                    while(kp<xp and cdiff<0 and sp>rmin_sampsize):
                        dbdiff=(ibin-dbin[kp]+12)%12
                        if dbdiff>1:
                            kp+=(dbdiff-1)*30
                            if kp>xp-1:
                                kp=xp-1
                                break
                        kp+=1
                        while kp<xp and count[kp,ibin]>count[kp-1,ibin] and sp>min_sampsize:
                            cdiff+=1 #=r+ count[kp,ibin]
                            kp+=1
                            #print('p',k,kp,ibin,dbin[kp],cdiff)
                            sp-=1
                    #r=count[xm,ibin]+(count[xp,ibin] - count[k,ibin])
                    while(km>xm and cdiff>0 and sm>min_sampsize):
                        dbdiff=(dbin[km]-ibin+12)%12
                        if dbdiff>1:
                            km-=(dbdiff-1)*30
                            if km<xm:
                                km=xm+1
                                break
                        km-=1
                        while km>0 and count[km+1,ibin]>count[km,ibin] and sm>min_sampsize:
                            cdiff-=1 #=count[km,ibin] - r
                            #print('m',k,km,ibin,dbin[km],cdiff)
                            km-=1
                            sm-=1

                    x+= (tmean[km,ibin] - tmean[xm,ibin]) #/ (count[km,ibin] - count[xm,ibin])  # Mittelwert 1 Periode
                    y+= (tmean[xp,ibin] - tmean[kp,ibin]) #/ (count[xp,ibin] - count[kp,ibin])  # Mittelwert 2 Periode
                    c+= (count[xp,ibin] - count[xm,ibin]-(count[kp,ibin]-count[km,ibin]))
                    cp+=count[xp,ibin] - count[kp,ibin]
                    cm+=count[km,ibin] - count[xm,ibin]
                    xy+= (tmean[xp,ibin] - tmean[xm,ibin]-(tmean[kp,ibin]-tmean[km,ibin])) #/ (count[xp,ibin] - count[xm,ibin]-(count[kp,ibin]-count[km,ibin]))  # Mittelwert ganze Periode
                    xx+= (tsquare[xp,ibin] - tsquare[xm,ibin]-(tsquare[kp,ibin]-tsquare[km,ibin])) #/ (count[xp,ibin] - count[xm,ibin]-(count[kp,ibin]-count[km,ibin]))  # Mittelwert ganze Periode
                
                if cm>=min_sampsize and cp>=min_sampsize:
                    
                    if c==0 or cm==0 or cp==0:
                        print(min_sampsize,c,cm,cp,'zero!')
                    sig = xx / np.float32(c)  # t*t ganze Periode
                    xym=xy/ np.float32(c)
                    if sig > xym * xym:
                        sig = sig - xym * xym  # standard deviation of the whole window
                        # n1 * (m1-m)**2 + n2 * (m2-m)**2 / stddev
                        if sig==0:
                            print(sig,'zero!')
                            bmax=0.
                            kmax=k
                        else:
                            tsa[k] = (np.float32(cm)* (x/np.float32(cm) - xym)*(x/np.float32(cm) - xym)  + np.float32(cp) * (y/np.float32(cp) - xym)*(y/np.float32(cp) - xym)) / sig
                            bmax=y/np.float32(cp)-x/np.float32(cm)
                            tsamax=tsa[k]
                            kmax=k
                    else:
                        tsa[k] = 0.
                else:
                    #print(k,cm,cp,'too small')
                    tsa[k]=0.
            else:
                tsa[k]=0.
    #if np.abs(bmax)>5:
        #print(bmax)
    return bmax
    

# @njit
def numba_meaneqsamp(test,ref,ni,istart,left_maxlen,right_maxlen,increment,miss_val,min_sampsize):
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
            ipc=np.count_nonzero(ip[:jstartlen])
            im[:jendlen]= np.not_equal(diff[j+left_maxlen:jend], miss_val)
            imc=np.count_nonzero(im[:jendlen])
        
        minlen=min(left_maxlen,right_maxlen)
        if(minlen-ipc  >  min_sampsize)  and  (minlen-imc  >  min_sampsize) :
            for ibin in range(0, 12,1):
                iph[:jstartlen]=np.logical_and(ip[:startlen], dbin[jstart:j+left_maxlen] == ibin)
                pmon[ibin]=np.count_nonzero(iph[:jstartlen])
                imh[:jstartlen]=np.logical_and(im[:jendlen], dbin[j+left_maxlen:jend] == ibin)
                mmon[ibin]=np.count_nonzero(imh[:jendlen])
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
            ipc=np.count_nonzero(ip[:jstartlen])
            imc=np.count_nonzero(im[:jendlen])
            if((minlen-ipc > min_sampsize) and (inlen-imc > min_sampsize)):
                for ibin in range(0, 12):
                    pmon[ibin]=np.logical_and(np.count_nonzero(ip[:jstartlen]), dbin[jstart:j+left_maxlen-1] == ibin)
                    mmon[ibin]=np.logical_and(np.count_nonzero(im[:jendlen]), dbin[j+left_maxlen:jend] == ibin)
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
def snhteqsamp(test,ref,ni,istart,istop,maxlen,increment,miss_val,min_sampsize,tsa,plus,minus,prms,mrms,pcount,mcount):
    
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

    for j in range(istart-min_sampsize, istop-min_sampsize, increment): # remove +1 of end, bc start is 0


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


        if((maxlen//2-ipc  >min_sampsize)  and  (maxlen//2-imc  >  min_sampsize)):
            for ibin in range(0, 12,1):
                iph=np.logical_and(ip[:jstartlen], dbin[jstart:j+m2-1] == ibin)
                pmon=np.sum(iph) # ??? not sure why it would be pmon(ibin), when only this iteration is used in the whole loop and it needs to get written also once per iteration
                imh=np.logical_and(im[:jendlen], dbin[j+m2:jend] == ibin)
                mmon=np.sum(imh) #??? ==^==
                
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
            ipc=np.sum(ip[:jstartlen])
            imc=np.sum(im[:jendlen])

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
        

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance in kilometers between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1 = np.radians(lon1)
    lat1 = np.radians(lat1)
    lon2 = np.radians(lon2)
    lat2 = np.radians(lat2)

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a)) 
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles. Determines return value units.
    return c * r


