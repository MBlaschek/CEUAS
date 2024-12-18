#!/usr/bin/env python
# coding: utf-8

import warnings
warnings.filterwarnings("ignore")#, category=DeprecationWarning) 
import numpy
import numpy as np
from numba import config,version_info, njit, prange
config.THREADING_LAYER = 'tbb'


@njit
def areg1(x,z,alpha):
    x[0]=z[0]
    for i in range(1,z.shape[0]):
        x[i]=alpha*x[i-1]+z[i]
    return

#@njit(parallel=True)
def break_simulator(x):
    m=x.shape[0]
    p=x.shape[1]
    n=x.shape[2]
    
    tsas=np.zeros(x.shape,dtype=np.float32)
    imax=np.zeros((m,p),dtype=np.int32)
    rmax=np.zeros((m,p),dtype=np.int32)
   
    for im in range(m):
        for ip in prange(p):
            tsa = np.zeros((x.shape[2]),dtype=np.float32)
            tmean = np.zeros((x.shape[2],12),dtype=np.float32)
            tsquare = np.zeros((x.shape[2],12),dtype=np.float32)
            count = np.zeros((x.shape[2],12),dtype=np.int32)

            z=np.random.randn(n)
            areg1(x[im,ip,:],z,0.0)
            x[im,ip,n//2:]+=im*0.2
            tsas[im,ip,:]=SNHT.numba_snhteqmov(x[im,ip,:].flatten(), tsa, 
                                               RC['snht_maxlen'],RC['snht_min_sampsize'],RC['miss_val'],
                                               count,tmean,tsquare)
            imax[im,ip]=np.argmax(tsas[im,ip,:])
            rmax[im,ip]=np.max(tsas[im,ip,:])
        #print(im)
    
    return imax,rmax,tsas  
    

@njit(fastmath={'nsz','arcp','contract','afn','reassoc'},cache=True)
def goodmon(var,idx,out):

    l=0
    for i in range(len(idx)-1):
        if idx[i+1]>idx[i]:
            l+=1
    #if l==0:
        #print('idx not valid')
    for ih in range(var.shape[0]):
        for ip in range(var.shape[1]):
            l=0
            for i in range(len(idx)-1):
                if idx[i+1]>idx[i]:
                    mc=0
                    for j in range(idx[i],idx[i+1]):
                        v=var[ih,ip,j]
                        if not np.isnan(var[ih,ip,j]):
                            mc+=1
                    out[ih,ip,l]=mc
                    l+=1
                
    return out
    
@njit(fastmath={'nsz','arcp','contract','afn','reassoc'},cache=True)
def monmean(var,idx,out,thresh=0,goodmon=goodmon):
    l=0
    for i in range(len(idx)-1):
        if idx[i+1]>idx[i]:
            l+=1
    #if l==0:
        #print('idx not valid')
    out[:]=np.nan
    for ih in range(var.shape[0]):
        for ip in range(var.shape[1]):
            l=0
            for i in range(len(idx)-1):
                if idx[i+1]>idx[i]:
                    mc=0
                    mm=np.float32(0.)
                    for j in range(idx[i],idx[i+1]):
                        v=var[ih,ip,j]
                        if v==v:
                            mm+=v
                            mc+=1
                    if mc>thresh:
                        out[ih,ip,l]=mm/mc
                    else:
                        out[ih,ip,l]=np.nan
                    l+=1
                
    return out[:,:,:l]

#@njit
def do_monmean(vars,tidx,fi):
    
    idx=np.searchsorted(fi['days'],tidx)
    mdays=idx[1:]-idx[:-1]
    midx=np.where(mdays>0)[0]
    
    fi['months']=tidx[midx]
            
    for v in vars:
        out=np.empty((fi[v].shape[0],fi[v].shape[1],fi['months'].shape[0]),dtype=fi[v].dtype)
        fi['m'+v]=monmean(fi[v],idx,out,RC['mon_thresh'])
        
    out=np.empty((fi[v].shape[0],fi[v].shape[1],fi['months'].shape[0]),dtype=np.int32)
    fi['goodmon']=goodmon(fi['temperatures'],idx,out)
            
                
    return


@njit
def addini(adjorig,iniad):
    
    ad=adjorig.shape
    res_with_ini=np.concatenate((adjorig[:,:,0:1],adjorig[:],np.zeros((ad[0],ad[1],1))),axis=2)
    for ih in range(ad[0]):
        for ip in range(ad[1]):
            res_with_ini[ih,ip,:-1]+=iniad[ih,ip]
            
    return res_with_ini


@njit
def make_adjusted_series(orig,adjustments,breaklist):
    
    adjusted_series=orig.copy()
    
    for ib in range(len(breaklist)-1,-1,-1):
        for ih in range(adjusted_series.shape[0]):
            for ip in range(adjusted_series.shape[1]):  
                if ib>0:
                    adjusted_series[ih,ip,breaklist[ib-1]:breaklist[ib]]-=adjustments[ih,ip,ib-1]
                else:
                    adjusted_series[ih,ip,:breaklist[ib]]-=adjustments[ih,ip,ib]
    
    return adjusted_series



@njit(fastmath={'nsz','arcp','contract','afn','reassoc'},cache=True)
def lagmean(ldiff,sdiff,ldays,sdays,llmaxlen,slmaxlen,thresh,lprofl,sprofl,lprofr,sprofr):
    
    idx=np.searchsorted(sdays,ldays)
    for ih in range(ldiff.shape[0]):
        for ip in range(ldiff.shape[1]-1,-1,-1):
            cl=0
            for i in range(llmaxlen):
                j=idx[i]
                if j==sdays.shape[0]:
                    break
                if ldiff[ih,ip,i]==ldiff[ih,ip,i]:
                    if sdays[j]==ldays[i] and sdiff[ih,ip,j]==sdiff[ih,ip,j]:
                        lprofl[ih,ip]+=ldiff[ih,ip,i]
                        sprofl[ih,ip]+=sdiff[ih,ip,j]
                        cl+=1
            if cl>thresh:
                lprofl[ih,ip]/=cl
                sprofl[ih,ip]/=cl
            else:
                if False: #ip<=6:
                    lprofl[ih,ip]=lprofl[ih,ip+1] # take value from next pressure level downward, assuming constant bias in the vertical
                    sprofl[ih,ip]=sprofl[ih,ip+1]
                else:    
                    lprofl[ih,ip]=np.nan
                    sprofl[ih,ip]=np.nan
                
            cr=0
            for i in range(llmaxlen,ldiff.shape[2]):
                j=idx[i]
                if j==sdays.shape[0]:
                    break
                if ldiff[ih,ip,i]==ldiff[ih,ip,i]:
                    if sdays[j]==ldays[i] and sdiff[ih,ip,j]==sdiff[ih,ip,j]:
                        lprofr[ih,ip]+=ldiff[ih,ip,i]
                        sprofr[ih,ip]+=sdiff[ih,ip,j]
                        cr+=1
            if cr>thresh:
                lprofr[ih,ip]/=cr
                sprofr[ih,ip]/=cr
            else:
                if ip<=6:
                    lprofr[ih,ip]=lprofr[ih,ip+1] # take value from next pressure level downward, assuming constant bias in the vertical
                    sprofr[ih,ip]=sprofr[ih,ip+1]
                else:
                    lprofr[ih,ip]=np.nan
                    sprofr[ih,ip]=np.nan
                
    return (lprofr-lprofl)-(sprofr-sprofl)
    


@njit
def calc_breakprofiles(adjlist,weight):
    breakprofile=np.zeros(adjlist.shape[1:],dtype=adjlist.dtype)
    
    for i in range(breakprofile.shape[0]):
        for j in range(breakprofile.shape[1]):
            bsum=0.
            wsum=0.
            count=0
            for k in range(adjlist.shape[0]):
                if adjlist[k,i,j]==adjlist[k,i,j]:
                    bsum+=adjlist[k,i,j]*weight[k]
                    wsum+=weight[k]
                    count+=1
            if count>1:
                breakprofile[i,j]=bsum/wsum
    return breakprofile
            

@njit
def findstart(fg_dep,minlen):
    i=fg_dep.shape[2]-1
    igood=np.zeros(2,dtype=np.int32)
    isave=np.zeros(2,dtype=np.int32)
    while i>0 and (igood[0]<minlen or igood[1]<minlen):
        for j in range(fg_dep.shape[0]):
            if fg_dep[j,11,i]==fg_dep[j,11,i]:
                igood[j]+=1
                if igood[j]==minlen:
                    isave[j]=i
        i-=1
        
    istart=i
    #if igood[0]>=minlen and igood[1]<minlen/8:
        #istart=isave[0]
    #if igood[1]>=minlen and igood[0]<minlen/8:
        #istart=isave[1]
    return istart,igood

@njit
def add_biasestimate(xyzv,xyzt,press,atime0,adj,adjpress):
    adjv=np.empty_like(xyzv)
    adjv.fill(np.nan)
    isec18=64800
    isec06=21600
    idtold=0
    for it in range(atime0.shape[0]):
        if it==atime0.shape[0]-1:
            idt=press.shape[0]
        else:
            idt=np.searchsorted(xyzt,atime0[it+1])
        for idx in range(idtold,idt):
            ih=xyzt[idx]%86400
            if adj.shape[0] == 2:
                
                if ih>=isec18 or ih<isec06:
                    ihi = 0
                else:
                    ihi = 1
            else:
                ihi = 0
            ip = 0
            while ip<adjpress.shape[0]:
                if press[idx] < adjpress[ip]: 
                    break  # no extrapolation upwards
                if press[idx]==adjpress[ip]:
                    adjv[idx]=adj[ihi,ip,it]
                    break
                else:
                    ipn = ip + 1
                    while ipn < adjpress.shape[0] and press[idx] > adjpress[ipn]:
                        ipn += 1
                    if ipn == adjpress.shape[0]:
                        break # no extrapolation downwards
                    w0 = (adjpress[ipn] - press[idx]) / (adjpress[ipn] - adjpress[ipn-1])
                    w1 = 1. - w0
                    adjv[idx]=adj[ihi,ipn-1,it] * w0 + adj[ihi,ipn,it] * w1
                    break
                    
                        
                    
        idtold=idt


    return adjv

