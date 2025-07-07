from numba import njit
import numpy as np
#import matplotlib.pyplot as plt
from scipy import stats
#import utils
import math
import os, sys, inspect
import pytz
import astral
import datetime

    # realpath() will make your script run, even if you symlink it :)
cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)

    # use this if you want to include modules from a subfolder
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"Rasotools/rasotools")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from utils import *

miss_val=np.nan

# pythran export manomaly(float[],float[],float[])
def manomaly( series,anomaly,climatology):
#   import np
#   from scipy.stats.stats import nanmean
    interval=np.int32([1980,2009])
    startyear=np.int32(1957)
    n=series.shape[0]
    for j in range(12):
        hanomaly=series[j:n:12]

        climatology[j]=sum(hanomaly[interval[0]-startyear:interval[1]-startyear])/(n/12)
        anomaly[j:n:12]=hanomaly-climatology[j]
    return

## pythran export anomalies_and_slopes(float[][][][],int,int[],int,int,int,float[],float[],float[][][][],float[],int[][][],float[][][][])
@njit(parallel=True,cache=True)
def anomalies_and_slopes(series,startyear,interval,tolerance,iens,itime,anomalies,climatologies,good,slopes):
    sshape=series.shape #
    for si in prange(sshape[0]):
        if sshape[2] > 5:
            gg = False
            for ipar in range(sshape[1]):
                g = 0
                for j in range((interval[0]-startyear)*12, (interval[1]-startyear+1)*12):
                    if ~np.isnan(series[si,ipar,10,j]):
                        g += 1
                if g > (interval[1] - interval[0]) * 6:
                    gg = True
            if gg:
                panomalies_and_slopes(si,series,startyear,interval,tolerance,iens,itime,anomalies,climatologies,good,slopes)
            else:
                anomalies[si, :, :, :] = np.nan
                climatologies[si, :, :, :] = np.nan
        else:
            
            panomalies_and_slopes(si,series,startyear,interval,tolerance,iens,itime,anomalies,climatologies,good,slopes) 
    return sshape[0]
@njit(cache=True,fastmath={'nsz','arcp','contract','afn','reassoc'})
def panomalies_and_slopes(si,series,startyear,interval,tolerance,iens,itime,anomalies,climatologies,good,slopes):

    sshape=series.shape
    climatology=np.zeros(12)
    anomaly=np.empty(sshape[3])
    orig=np.empty(sshape[3])
    stop=interval[1]
    start=interval[0]
#    for si in prange(sshape[0]):
    for ipar in range(sshape[1]):
        for ip in range(sshape[2]):
            g = 0
            for j in range(sshape[3]):
                orig[j]=series[si,ipar,ip,j]
                if orig[j] == orig[j]:
                    g += 1
            if g > (interval[1] - interval[0]) * 6:
                
                good[si,ipar,ip]=janomaly(orig,startyear,interval,anomaly,climatology)
            else:
                good[si,ipar,ip] = g
            slopes[si,iens,ipar,ip]  = miss_val
#            print(good[si,ipar,ip],(stop-start+1-tolerance)*12)
            if good[si,ipar,ip]>(stop-start+1-tolerance)*12:
                for j in range(sshape[3]):
                    anomalies[si,ipar,ip,j]=anomaly[j]
                for j in range(12):
                    climatologies[si,ipar,ip,j]=climatology[j]
                tgood=0
                for j in range(itime.shape[0]):
                    orig[j]=anomaly[itime[j]]     
                    if orig[j]==orig[j]:
                        tgood+=1
                if tgood>(stop-start+1-tolerance)*12:
                    val  = fastlinregress(itime,orig)
                #if si==1171 and ipar==0 and ip==0:
                    #print val*120
                    slopes[si,iens,ipar,ip]  = val*120.
            else:
                for j in range(sshape[3]):
                    anomalies[si,ipar,ip,j]=miss_val
                for j in range(12):
                    climatologies[si,ipar,ip,j]=miss_val




    return

# this routine expects two series as input, the series and the reference series.
# it calculates the mean difference, given the interval and the tolerance.
# only the difference is returned.
@njit(cache=True)
def calc_diffs(series,refseries,startyear,interval,tolerance,iens,good,slopes):

    sshape=series.shape
    stop=interval[1]
    start=interval[0]
    while (stop-startyear+1)*12>series.shape[3]:
        stop-=1
    
    for si in range(sshape[0]):
        for ipar in range(sshape[1]):
            for ip in range(sshape[2]):
                b=0
                g=0
                s=0.
                for i in range((start-startyear)*12,(stop-startyear+1)*12):
                    if series[si,ipar,ip,i]==series[si,ipar,ip,i] and refseries[si,ipar,ip,i]==refseries[si,ipar,ip,i]:
                        s+=series[si,ipar,ip,i]-refseries[si,ipar,ip,i]
                        g+=1
                    else:
                        b+=1
                good[si,ipar,ip]=g
                if b>=tolerance or g==0:
                    slopes[si,iens,ipar,ip]=np.nan
                else:
                    slopes[si,iens,ipar,ip]=s/g

    return

#@njit(cache=True)
def guanomalies_and_slopes(series,startyear,interval,tolerance,iens,itime,orig,anomaly,anomalies,climatology,good,slopes):

    sshape=series.shape
    stop=interval[1]-startyear+1
    start=interval[0]-startyear
    jansimple( series,start,stop,anomalies,good)
    gulinregress(np.asarray(itime,dtype='float64'),anomalies[:,:,:,itime],good,np.asarray([tolerance]),slopes[:,0,:2,:anomalies.shape[2]])
    slopes*=120
    return
    for si in range(sshape[0]):
        for ipar in range(sshape[1]):
            for ip in range(sshape[2]):
                for j in range(sshape[3]):
                    orig[j]=series[si,ipar,ip,j]
#                good[si,ipar,ip]=janomaly(orig,startyear,interval,anomaly,climatology)

                slopes[si,iens,ipar,ip]  = miss_val
                if good[si,ipar,ip]>(stop-start+1-tolerance)*12:
                    for j in range(sshape[3]):
                        anomalies[si,ipar,ip,j]=anomaly[j]
                    tgood=0
                    for j in range(itime.shape[0]):
                        orig[j]=anomaly[itime[j]]     
                        if orig[j]==orig[j]:
                            tgood+=1
                    if tgood>(stop-start+1-tolerance)*12:
                        val  = fastlinregress(itime,orig)
                    #if si==1171 and ipar==0 and ip==0:
                        #print val*120
                        slopes[si,iens,ipar,ip]  = val*120.
                else:
                    for j in range(sshape[3]):
                        anomalies[si,ipar,ip,j]=miss_val



    return


## pythran export janomaly(float[],int,int[],float[],float[])
@njit(cache=True,fastmath={'nsz','arcp','contract','afn','reassoc'})
def janomaly( series,startyear,interval,anomaly,climatology):
    n=series.shape[0]
    start=interval[0]-startyear
    stop=interval[1]-startyear+1
    good=0
#   print range(start,stop)
    for j in range(12):
        l=0
        climatology[j]=0.
        for k in range(start,stop):
            s=series[k*12+j]
            if s==s:
                climatology[j]+=s
                l=l+1
        if l > (stop-start+1)/2:
            climatology[j]/=l
            for k in range(n//12):
                anomaly[k*12+j]=series[k*12+j]-climatology[j]
                if anomaly[k*12+j]==anomaly[k*12+j]: 
                    good+=1
        else:
            climatology[j]=miss_val
            for k in range(n//12):
                anomaly[k*12+j]=miss_val


    return good

#@guvectorize([(float64[:,:],float64[:,:])],'(m,n)->(m,n)', target='parallel')
def vanomalies(vseries,vanomalies):
    good=0
    for i in range(vseries.shape[0]):
        vanomalies[i,:]=vseries[i,:]
        #jandummy(vseries[i,:],vanomalies[i,:])
        #        jansimple(vseries[i,:],vanomalies[i,:])
        
#@njit(void(float64[:],float64[:]) ,cache=True)
#@guvectorize([(float64[:],float64[:])],'(n)->(n)', target='parallel')
def jandummy(series,anomaly):
    anomaly[:]=series[:]

#@njit(cache=True)
#@guvectorize([(float64[:],int64,int64,float64[:],int64[:])],'(n),(),()->(n),()', target='parallel',nopython=True)
def jansimple( series,start,stop,anomaly,good):
    climatology=np.zeros(12)
    n=series.shape[0]
    good[0]=0
#   print range(start,stop)
    for j in range(12):
        l=0
        climatology[j]=0.
        for k in range(start,stop):
            s=series[k*12+j]
            if s==s:
                climatology[j]=climatology[j]+s
                l=l+1
        if l > (stop-start+1)/2:
            climatology[j]/=l
            for k in range(n/12):
                anomaly[k*12+j]=series[k*12+j]-climatology[j]
                if anomaly[k*12+j]==anomaly[k*12+j]: 
                    if k>=start and k<stop: 
                        good[0]+=1
        else:
            climatology[j]=miss_val
            for k in range(n/12):
                anomaly[k*12+j]=miss_val


    return

# pythran export nnanstd(float[])
@njit(cache=True)
def nnanstd(field):

    s=0.
    ss=0.
    c=0
    for i in range(field.shape[0]):
        if field[i]==field[i]:
            s+=field[i]
            ss+=field[i]*field[i]
            c=c+1
    return math.sqrt(ss/c-s*s/c/c)

# pythran export nnanmean(float[])
@njit(cache=True,fastmath={'nsz','arcp','contract','afn','reassoc'})
def nnanmean(field,thresh=0):

    s=0.
    c=0
    for i in range(field.shape[0]):
        if field[i]==field[i]:
            s+=field[i]
            c=c+1
    if c>0:
        s/=c
    else:
        s=miss_val

    return s

@njit(cache=True,fastmath={'nsz','arcp','contract','afn','reassoc'})
def ndnanmean(field,thresh):

    fieldm=np.zeros((field.shape[0],field.shape[3],field.shape[4]))
    for l in range(field.shape[0]):
        for k in range(field.shape[3]):
            for j in range(field.shape[4]):
                s=0.
                c=0
                for m in range(field.shape[1]):
                    for i in range(field.shape[2]):
                        f=field[l,m,i,k,j]
                        if f==f:
                            s+=f
                            c=c+1
                if c>thresh:
                    fieldm[l,k,j]=s/c
                else:
                    fieldm[l,k,j]=miss_val

    return fieldm

@njit(cache=True,fastmath={'nsz','arcp','contract','afn','reassoc'})
def ndnanmean2(field,fieldm,thresh):

    #fieldm=np.zeros(field.shape[:2]) #,dtype=field.dtype)
    for l in range(field.shape[0]):
        for k in range(field.shape[1]):
            s=0.
            c=0
            for i in range(field.shape[2]):
                f=field[l,k,i]
                if f==f:
                    s+=f
                    c=c+1
            if c>thresh:
                fieldm[l,k]=s/c
            else:
                fieldm[l,k]=miss_val

    return fieldm

# calculate variance of mean of time series (which is in index 2)
@njit(cache=True,fastmath={'nsz','arcp','contract','afn','reassoc'})
def ndnanvar2(field,fieldvar,thresh):

    #fieldvar=np.zeros(field.shape[:2],dtype=field.dtype)
    for l in range(field.shape[0]):
        for k in range(field.shape[1]):
            s=0.
            sq=0.
            c=0
            for i in range(field.shape[2]):
                f=field[l,k,i]
                if f==f:
                    s+=f
                    sq+=f*f
                    c=c+1
            if c>thresh:
                fieldvar[l,k]=(sq/c-(s*s/c/c))/c # divide by c for variance of mean
            else:
                fieldvar[l,k]=miss_val

    return fieldvar

# pythran export nnanmeanlist0(float[][][][],float[][][][],int[][][],int list)
@njit(cache=True)
def nnanmeanlist0(field,fieldm,counts,index):

    for i in index:
        for j in range(field.shape[1]):
            for k in range(field.shape[2]):
                for l in range(field.shape[3]):
                    f=field[i,j,k,l]
                    if f==f:
                        fieldm[0,j,k,l]+=f
                        counts[j,k,l]+=1
    for j in range(field.shape[1]):
        for k in range(field.shape[2]):
            for l in range(field.shape[3]):
                if counts[j,k,l]>0:
                    fieldm[0,j,k,l]/=counts[j,k,l]
                else:
                    fieldm[0,j,k,l]=miss_val

    return

## pythran export xyanomalies(float[][][][],float[][][][],float[][][][],int,int[],float[],float[],int)
@njit(cache=True)
def xyanomalies(field,fieldplus,fieldminus,syear,interval,anomaly,climatology,hilf,mlistlen):

    for i in range(field.shape[2]):
        for j in range(field.shape[3]):
            for k in range(field.shape[0]):
                hilf[k]=field[k,0,i,j]
            n=janomaly(hilf,syear,interval,anomaly,climatology)
            if n>1:
                alpha=alphaest(anomaly)
                if alpha<1:
                    delta=1.64*nnanstd(anomaly)/math.sqrt(n*(1-alpha)/(1+alpha)/12*mlistlen)
                else:
                    delta=0.
                fieldminus[0,0,i,j]-=delta
                fieldplus[0,0,i,j]+=delta
            else:
                fieldminus[0,0,i,j]=miss_val
                fieldplus[0,0,i,j]=miss_val

## pythran export xyanomalies(float[][][][],float[][][][],float[][][][],int,int[],float[],float[],int[])
#@njit(cache=True)
def xyanomaliesmbb(field,fieldplus,fieldminus,syear,interval,anomaly,climatology,hilf,season):

    n=40
    ind=np.random.randint(field.shape[0],size=[field.shape[0]*n])
    Xpos=np.zeros([n,field.shape[0]],dtype=int)
    Xres=np.zeros([field.shape[0]])
    est=np.zeros([n])
    hilf[:]=np.nan
    try:
        xyanomaliesmbbcore(field,fieldplus,fieldminus,syear,interval,anomaly,climatology,n,ind,Xpos,Xres,est,hilf,season)
    except ZeroDivisionError:
        fieldplus[:]=np.nanmean(field,axis=0)
        fieldminus[:]=fieldplus[:]

    return

## pythran export xyanomalies(float[][][][],float[][][][],float[][][][],int,int[],float[],float[],int,int,int,int,int,float[],int[])
@njit(cache=True)
def xyanomaliesmbbcore(field,fieldplus,fieldminus,syear,interval,anomaly,climatology,n,ind,Xpos,Xres,est,hilf,season):

    for i in range(field.shape[2]):
#       print i
        for j in range(field.shape[3]):
            for k in range(field.shape[0]/12):
                for m in range(season.shape[0]):
                    hilf[k*12+season[m]]=field[k*12+season[m],0,i,j]
            l=-1
            janomaly(hilf,syear,interval,anomaly,climatology)
            mbbfast(anomaly,Xres,n,l,ind,Xpos,est)
#           est,Xpos=utils.mbb(anomaly,n,func=nnanmean,ind=ind,Xpos=Xpos)
#           sest=sorted(est)
#            bubblesort(est)
            np.sort(est)
            fieldminus[0,0,i,j]+=est[int(np.floor(0.05*n))]
            fieldplus[0,0,i,j]+=est[int(np.floor(0.95*n))]

    return

def dailyanomalydriver(d,ens,plotinterval):
    h=d["ddata"]
    sh=h.shape
    anomaly=np.zeros(sh[4])
    climatology=np.zeros(365)
    for ist in range(sh[0]):
        for ipar in range(sh[2]):
            for ip in range(sh[3]):
                good=dailyanomaly(h[ist,ens,ipar,ip,:],np.int(d["startdate"]/10000),
                                  np.array(plotinterval),anomaly,climatology)
                h[ist,ens,ipar,ip,:]=anomaly

    return h[:,ens,:,:,:]

@njit(cache=True)
def danomaly(series,bins,ccount,anomaly,climatology):
    n=series.shape[2]
    s=0
    #   print range(start,stop)
    for ih in range(series.shape[0]):
        for ip in range(series.shape[1]):
            ccount[:]=0
            for i in range(n):
                if series[ih,ip,i]==series[ih,ip,i]:
                    j=bins[i]
                    climatology[ih,ip,j]+=series[ih,ip,i]
                    ccount[j]+=1
                    s+=1
            for j in range(12):
                if ccount[j]>10:
                    climatology[ih,ip,j]/=ccount[j]
                else:
                    climatology[ih,ip,j]=np.nan
            for i in range(n):
                if series[ih,ip,i]==series[ih,ip,i]:
                    j=bins[i]
                    anomaly[ih,ip,i]=series[ih,ip,i]-climatology[ih,ip,j]
    
    return s

##pythran export dailyanomaly(float[],int,int[],float[],float[])
#@njit(cache=True)
def dailyanomaly( series,startyear,interval,anomaly,climatology):
    n=series.shape[0]
    start=int(interval[0]-startyear)
    stop=int(interval[1]-startyear+1)
    good=0
#   print range(start,stop)
    for j in range(365):
        l=0
        climatology[j]=0.
        for k in range(start,stop):
            ik=int(k*365.25)+j
            s=series[ik]
            if s==s:
                climatology[j]=climatology[j]+s
                l=l+1
        if l > 5: #(stop-start+1)/2:
            climatology[j]=climatology[j]/l
            for k in range(n/365):
                ik=int(k*365.25)+j
                anomaly[ik]=series[ik]-climatology[j]
                if anomaly[ik]==anomaly[ik]: 
                    good+=1
        else:
            climatology[j]=miss_val
            for k in range(n/365):
                anomaly[ik]=miss_val

        if good>0:
            pass


    return good

@njit(cache=False,parallel=True)
def ncat(a,b,c):
    for i in prange(a.shape[0]):
        pncat(a,b,c,i)
        
@njit(cache=True)
def pncat(a,b,c,i):
#    for i in prange(a.shape[0]):
        for j in range(a.shape[1]):
            for k in range(a.shape[2]):
                for l in range(a.shape[3]):
                    c[i,j,k,l]=a[i,j,k,l]
        ash=a.shape[2]
    #for i in prange(a.shape[0]):
        for j in range(a.shape[1]):
            for k in range(b.shape[2]):
                for l in range(a.shape[3]):
                    c[i,j,k+ash,l]=b[i,j,k,l]
                    
        return 

@njit(parallel=True)
def grid_anomalies(anomaly,good,tolerance,gstatindex,ganomalies,gslopes,start,stop,itime):
    nj=gstatindex.shape[0]
    for j in prange(nj):
        pgrid_anomalies(anomaly,good,tolerance,gstatindex,ganomalies,gslopes,start,stop,itime,j)

@njit(cache=True)
def pgrid_anomalies(anomaly,good,tolerance,gstatindex,ganomalies,gslopes,start,stop,itime,j):
    nj=gstatindex.shape[0]
    ni=gstatindex.shape[1]
    orig=np.empty(itime.shape[0])
    gst=np.empty(anomaly.shape[3],dtype=np.int32)

#    for j in range(nj):
    for i in range(ni):
        nh=gstatindex[j,i,0]
        for ipar in range(2):
            for ip in range(anomaly.shape[2]):
                if ip<gslopes.shape[3]:
                    gslopes[j,i,ipar,ip]  = miss_val
                    if nh>0:
                        for m in range(anomaly.shape[3]):
                            gst[m]=0
                            ganomalies[j,i,ipar,ip,m]=0.
                        for h in range(nh):
                            si=gstatindex[j,i,h+1]
                            g=0
                            for m in range(itime.shape[0]):
                                if itime[m]<anomaly.shape[3]:
                                    if anomaly[si,ipar,ip,itime[m]]==anomaly[si,ipar,ip,itime[m]]:
                                        g+=1
                            if g>=(stop-start+1-tolerance)*12:
                                for m in range(anomaly.shape[3]):
                                    a=anomaly[si,ipar,ip,m]
                                    if a==a:
                                        ganomalies[j,i,ipar,ip,m]+=a
                                        gst[m]+=1
    
                        for m in range(anomaly.shape[3]):
                            if gst[m]>0:
                                ganomalies[j,i,ipar,ip,m]/=gst[m]
                            else:
                                ganomalies[j,i,ipar,ip,m]=miss_val
    
                        tgood=0
                        for k in range(itime.shape[0]):
                            if itime[k] < ganomalies.shape[4]:
                                orig[k]=ganomalies[j,i,ipar,ip,itime[k]]     
                                if orig[k]==orig[k]:
                                    tgood+=1
                        if tgood>(stop-start+1-tolerance)*12 and anomaly.shape[3]>12:
                            val  = fastlinregress(itime,orig)
                            gslopes[j,i,ipar,ip]  = val*120.
    
                    else:
                        for m in range(anomaly.shape[3]):
                            ganomalies[j,i,ipar,ip,m]=miss_val

    return

@njit(cache=True,parallel=True)
def grid_diffs(diff,good,tolerance,gstatindex,gdiffs,iens):
    nj=gstatindex.shape[0]
    for j in prange(nj):
        pgrid_diffs(diff,good,tolerance,gstatindex,gdiffs,iens,j)

@njit(cache=True,parallel=False)
def pgrid_diffs(diff,good,tolerance,gstatindex,gdiffs,iens,j):
    nj=gstatindex.shape[0]
    ni=gstatindex.shape[1]

#    for j in range(nj):
    for i in range(ni):
        nh=gstatindex[j,i,0]
        for ipar in range(diff.shape[2]):
            for ip in range(gdiffs.shape[3]):
                gdiffs[j,i,ipar,ip]  = miss_val
                if nh>0:
                    gst=0
                    gdiffs[j,i,ipar,ip]=0.
                    for h in range(nh):
                        si=gstatindex[j,i,h+1]
                        if diff[si,iens,ipar,ip]==diff[si,iens,ipar,ip]:
                            gdiffs[j,i,ipar,ip]+=diff[si,iens,ipar,ip]
                            gst+=1

                    if gst>0:
                        gdiffs[j,i,ipar,ip]/=gst
                    else:
                        gdiffs[j,i,ipar,ip]=miss_val

                else:
                    gdiffs[j,i,ipar,ip]=miss_val

    return

@njit(cache=False,parallel=True)
def allhad(hadens,hadtem,hadmed,hadslopes,hadtime,first,last,hadstop):
    for iens in prange(hadens.shape[0]):
        pallhad(hadens,hadtem,hadmed,hadslopes,hadtime,first,last,hadstop,iens)

@njit(cache=True)
def pallhad(hadens,hadtem,hadmed,hadslopes,hadtime,first,last,hadstop,iens):
    for ilat in range(hadmed.shape[2]):
        for ilon in range(hadmed.shape[3]):
            hadslopes[iens,ilat,ilon]=fastlinregress(hadtime[0:hadstop],hadens[iens,first:last,ilat,ilon]
                                                     )*120.

@njit(fastmath={'nsz','arcp','contract','afn','reassoc'},cache=True)
def fastlinregress(x,y):
#    print typeof(x),typeof(y)
    m=x.shape[0]
    
    xmean=0.
    ymean=0.
    xcov=0.
    xsq=0.
    n=0
    for k in range(m):
        if x[k]==x[k] and y[k]==y[k]:
            xmean+=x[k]
            ymean+=y[k]
            n+=1
    if n>2:   
        xmean/=n
        ymean/=n
        for k in range(m):
            if x[k]==x[k] and y[k]==y[k]:
                xcov+=x[k]*y[k]
                xsq+=x[k]*x[k]
        slope=(xcov-n*xmean*ymean)/(xsq-n*xmean*xmean)
    else:
        slope=miss_val

    return slope

#@guvectorize([(float64[:],float64[:],int64[:],int64,float64[:])],'(n),(n),(),()->()',target='parallel',nopython=True)
def gulinregress(x,y,good,tolerance,slope):
#    print typeof(x),typeof(y)
    m=x.shape[0]
    if good[0]<=m-12*tolerance:
        slope[0]=np.nan
        return
    xmean=0.
    ymean=0.
    xcov=0.
    xsq=0.
    n=0
    for k in range(m):
        if x[k]==x[k] and y[k]==y[k]:
            xmean+=x[k]
            ymean+=y[k]
            n+=1
    if n>2:   
        xmean/=n
        ymean/=n
        for k in range(m):
            if x[k]==x[k] and y[k]==y[k]:
                xcov+=x[k]*y[k]
                xsq+=x[k]*x[k]
        slope[0]=(xcov-n*xmean*ymean)/(xsq-n*xmean*xmean)
    else:
        slope[0]=np.nan

    return

@njit(cache=True)
def fastlinregressnonanerr(x,y):
    m=x.shape[0]

    xmean=0.
    ymean=0.
    xcov=0.
    xsq=0.
    ysq=0.
    if m>2:   
        for k in range(m):
            xmean+=x[k]
            ymean+=y[k]
        xmean/=m
        ymean/=m
        for k in range(m):
            xcov+=x[k]*y[k]
            xsq+=x[k]*x[k]
            ysq+=y[k]*y[k]
        slope=(xcov-m*xmean*ymean)/(xsq-m*xmean*xmean)
        intercept=ymean-slope*xmean
        sig=(ysq-m*ymean*ymean-slope*slope*(xsq-m*xmean*xmean))/(m-2)
        serr=sig/(xsq-m*xmean*xmean)
    else:
        slope=miss_val
        serr=miss_val

    return slope,np.sqrt(serr)

@njit(cache=True)
def fastlinregressnonan(x,y):
    m=x.shape[0]

    xmean=0.
    ymean=0.
    xcov=0.
    xsq=0.
    if m>2:   
        for k in range(m):
            xmean+=x[k]
            ymean+=y[k]
        xmean/=m
        ymean/=m
        for k in range(m):
            xcov+=x[k]*y[k]
            xsq+=x[k]*x[k]
        slope=(xcov-m*xmean*ymean)/(xsq-m*xmean*xmean)
    else:
        slope=miss_val

    return slope




def mbb(X,n,ind=None,l=None,func=np.mean,timearg=None,Xpos=None):
    """Moving Block Bootstrap as described in Mudelsee, 2011"""

#    t1=time.time()
    if timearg is None:
        Xresid=X
    else:
        slope=func(timearg,X)
        Xresid=X-slope*timearg

    if l==None:
        alpha=alphaest(Xresid)
        alpha=np.abs(alpha)
        l=int((math.sqrt(6.)*alpha/(1-alpha*alpha))**(2./3.)*X.shape[0]**(1./3.))
        if l<1:
            l=1

#    print alpha,l    
    if(ind==None):
        ind=np.random.randint(X.shape[0]-l+1,size=[X.shape[0]*n/l])

    est=np.zeros(n)
    pos=0
    #Xres=np.zeros([n,X.shape[0]])
    #for i in range(n):
        #pos=block_resample(Xres[i,:],Xresid,ind,l,pos)
        #if timearg is not None:
            #est[i]=slope+func(timearg,Xres)
        #else:
            #est[i]=func(Xres[i,:])
    Xres=np.zeros([X.shape[0]])
    if Xpos is not None:
        if Xpos[0,0]==0 and Xpos[n-1,X.shape[0]-1]==0:
            Xpos=np.zeros([n,X.shape[0]],dtype=int)
    else:
        Xpos=np.zeros([n,X.shape[0]],dtype=int)

    pos=fastblock_resample_mean(est,Xres,Xpos,Xresid,ind,l,pos)
#    print 'mbb: ',time.time()-t1

#    if np.nanstd(Xresid)>0:
#        print est
    return est,Xpos

@njit(cache=True)
def mbbfast(X,Xres,n,l,ind,Xpos,est):
    """Moving Block Bootstrap as described in Mudelsee, 2011"""

    if l==-1:
        alpha=alphaest(X)
        alpha=np.abs(alpha)
        l=int((math.sqrt(6.)*alpha/(1-alpha*alpha))**(2./3.)*X.shape[0]**(1./3.))
        if l<1:
            l=1
    pos=0

    pos=fastblock_resample_mean(est,Xres,Xpos,X,ind,l,pos)

    return

@njit(cache=True)
def areg1(X,Z,alpha):

    n=Z.shape[0]
    X[0]=Z[0]
    for i in range(n-1):
        X[i+1]=alpha*X[i]+Z[i]

    return

@njit(cache=True)
def alphaest(X):
# X must have zero mean
    rho=0.
    var=0.
    n=X.shape[0]
    m=0
    for i in range(n-1):
        if X[i]==X[i] and X[i+1]==X[i+1]:
            rho+=X[i]*X[i+1]
            var+=X[i]*X[i]
            m=m+1
    if var>0:
        if m>4:
            alpha=rho/var
        else:
            alpha=0.
    else:
        alpha=0.

    if abs(alpha)>1:
        alpha=0.

    return alpha

@njit(cache=True)
def bootstrap_resample(A,X,ind,pos):

    n=X.shape[0]
    for i in range(n):
        for j in range(n):
            A[i,j]=X[ind[pos+j]]
        pos=pos+n
    return pos

@njit(cache=True)
def block_resample(A,X,ind,l,pos):

    n=X.shape[0]
    jpos=0
    for j in range((n-l)/l):
        for i in range(l):
            A[jpos+i]=X[ind[pos+j]+i]
        jpos+=l
    pos=pos+(n-l)/l
    return pos

@njit(cache=True)
def fastblock_resample_mean(est,Y,A,X,ind,l,pos):

    n=A.shape[1]
    if A[0,0]==0 and A[A.shape[0]-1,n-1]==0:
        for k in range(A.shape[0]):
            jpos=0
            for j in range((n-l)/l):
                for i in range(l):
                    A[k,jpos+i]=ind[pos+j]+i
                jpos+=l
            pos=pos+(n-l)/l
    for k in range(A.shape[0]):
        for j in range(n):
            Y[j]=X[A[k,j]]
        est[k]=nnanmean(Y)


    return pos

def alb(vflux):

    s=vflux.shape
    zmflux=vflux[:]
    zmflux[vflux<0]=0.
    zmflux[vflux>1]=1.

    return zmflux

def globalflux(vflux):

    s=vflux.shape
    zmflux=globalfluxub(vflux,spread=False)
    for il in range(1,s[2]):
        zmflux[:,0,il,0]-=zmflux[:,0,s[2]-1,0]*il/s[2]
    for i in range(1,s[3]):
        zmflux[:,:,:,i]=zmflux[:,:,:,0]

    return zmflux

def tend(av):

    t=np.arange(144)/120.
    s=np.empty((12,av.shape[2],av.shape[3]))
    for im in range(12):
            for iy in range(av.shape[2]):
                for ix in range(av.shape[3]):
        #            s[iy,ix]=fastlinregress(t,av[:,0,iy,ix])
                    s[im,iy,ix]=(av[-12+im,0,iy,ix]-av[im,0,iy,ix])
    
    for it in range(0,av.shape[0]/12,12):
        for im in range(12):
            av[it*12+im,0,:,:]=s[im,:,:]
       
    return av

def globalfluxdreh(vflux):

    s=vflux.shape
    zmflux=globalfluxubdreh(vflux,spread=False)
    for il in range(1,s[2]):
        zmflux[:,0,il,0]-=zmflux[:,0,s[2]-1,0]*il/s[2]
    for i in range(1,s[3]):
        zmflux[:,:,:,i]=zmflux[:,:,:,0]

    return zmflux

def globalfluxcorr(vflux):

    s=vflux.shape
    zmflux=globalfluxub(vflux,spread=False)
    for il in range(1,s[2]):
        zmflux[:,0,il,0]=-zmflux[:,0,s[2]-1,0]*il/s[2]
    for i in range(1,s[3]):
        zmflux[:,:,:,i]=zmflux[:,:,:,0]

    return zmflux

def globalfluxub(vflux,spread=True):

    s=vflux.shape
    erad=6370000.
    coslats=np.sin((-90.+np.arange(s[2]+1)*180./s[2])*math.pi/180.)
    deltacos=coslats[1:s[2]+1]-coslats[0:s[2]]
    zmflux=np.zeros(s)
    zmflux[:,:,:,0]=np.nanmean(vflux,axis=3)
    for it in range(s[0]):
        zmflux[it,0,:,0]=zmflux[it,0,:,0]*deltacos*2*math.pi*erad*erad
    zmflux[:,0,0,0]=zmflux[:,0,1,0]/2.
    for il in range(1,s[2]):
        zmflux[:,0,il,0]=zmflux[:,0,il-1,0]+zmflux[:,0,il,0]
    if spread:
        for i in range(1,s[3]):
            zmflux[:,:,:,i]=zmflux[:,:,:,0]

    return zmflux

def globalfluxubdreh(vflux,spread=True):

    s=vflux.shape
    erad=6370000.
    coslats=np.sin((-90.+np.arange(s[2]+1)*180./s[2])*math.pi/180.)
    deltacos=coslats[1:s[2]+1]-coslats[0:s[2]]
    coslats=np.cos((-90.+np.arange(s[2]+1)*180./s[2])*math.pi/180.)
    coslats2=(coslats[1:s[2]+1]+coslats[0:s[2]])/2.
    zmflux=np.zeros(s)
    zmflux[:,:,:,0]=np.nanmean(vflux,axis=3)
    for it in range(s[0]):
        zmflux[it,0,:,0]=zmflux[it,0,:,0]*coslats2*deltacos*2*math.pi*erad*erad*erad
    zmflux[:,0,0,0]=zmflux[:,0,1,0]/2.
    for il in range(1,s[2]):
        zmflux[:,0,il,0]=zmflux[:,0,il-1,0]+zmflux[:,0,il,0]
    if spread:
        for i in range(1,s[3]):
            zmflux[:,:,:,i]=zmflux[:,:,:,0]

    return zmflux


def meanflux(vflux):

    s=vflux.shape
    zmflux=meanfluxub(vflux,spread=False)
    for il in range(1,s[2]):
        zmflux[:,0,il,0]-=zmflux[:,0,s[2]-1,0]*il/s[2]
    for i in range(1,s[3]):
        zmflux[:,:,:,i]=zmflux[:,:,:,0]

    return zmflux

def meanfluxub(vflux,spread=True):

    s=vflux.shape
    erad=6370000.
    coslats=np.cos((-90.+np.arange(s[2]+1)*180./s[2])*math.pi/180.)
    mcoslats=0.5*(coslats[1:s[2]+1]+coslats[0:s[2]])
    deltacos=coslats[1:s[2]+1]-coslats[0:s[2]]
    zmflux=np.zeros(s)
    zmflux[:,:,:,0]=np.nanmean(vflux,axis=3)
    for it in range(s[0]):
        zmflux[it,0,:,0]=zmflux[it,0,:,0]*2*mcoslats*math.pi*erad#/6.37
#    zmflux[:,0,0,0]=zmflux[:,0,1,0]/2.
#    for il in range(1,s[2]):
#        zmflux[:,0,il,0]=zmflux[:,0,il-1,0]+zmflux[:,0,il,0]
    if spread:
        for i in range(1,s[3]):
            zmflux[:,:,:,i]=zmflux[:,:,:,0]

    return zmflux
def meanfluxubdreh(vflux,spread=True):

    s=vflux.shape
    erad=6370000.
    coslats=np.cos((-90.+np.arange(s[2]+1)*180./s[2])*math.pi/180.)
    mcoslats=0.5*(coslats[1:s[2]+1]+coslats[0:s[2]])
    deltacos=coslats[1:s[2]+1]-coslats[0:s[2]]
    zmflux=np.zeros(s)
    zmflux[:,:,:,0]=np.nanmean(vflux,axis=3)
    for it in range(s[0]):
        zmflux[it,0,:,0]=zmflux[it,0,:,0]*2*mcoslats*mcoslats*math.pi*erad*erad*erad/6.37
#    zmflux[:,0,0,0]=zmflux[:,0,1,0]/2.
#    for il in range(1,s[2]):
#        zmflux[:,0,il,0]=zmflux[:,0,il-1,0]+zmflux[:,0,il,0]
    if spread:
        for i in range(1,s[3]):
            zmflux[:,:,:,i]=zmflux[:,:,:,0]

    return zmflux

def todreh(vflux):

    s=vflux.shape
    erad=6370000.
    coslats=np.cos((-90.+np.arange(s[2]+1)*180./s[2])*math.pi/180.)
    mcoslats=0.5*(coslats[1:s[2]+1]+coslats[0:s[2]])
    deltacos=coslats[1:s[2]+1]-coslats[0:s[2]]
    zmflux=np.zeros(s)
    for il in range(s[2]):
        zmflux[:,:,il,:]=vflux[:,:,il,:]*mcoslats[il]*erad

    return zmflux-np.mean(zmflux)/3.

def tanfluxub(vflux,spread=True):

    s=vflux.shape
    erad=6370000.
    coslats=np.cos((-90.+np.arange(s[2]+1)*180./s[2])*math.pi/180.)
    sinlats=np.sin((-90.+np.arange(s[2]+1)*180./s[2])*math.pi/180.)
    mcoslats=0.5*(coslats[1:s[2]+1]+coslats[0:s[2]])
    msinlats=0.5*(sinlats[1:s[2]+1]+sinlats[0:s[2]])
    tanlats=msinlats/mcoslats
    deltacos=coslats[1:s[2]+1]-coslats[0:s[2]]
    zmflux=np.zeros(s)
    zmflux[:,:,:,0]=np.nanmean(vflux,axis=3)
    for it in range(s[0]):
        zmflux[it,0,:,0]=zmflux[it,0,:,0]*tanlats/erad
#    zmflux[:,0,0,0]=zmflux[:,0,1,0]/2.
#    for il in range(1,s[2]):
#        zmflux[:,0,il,0]=zmflux[:,0,il-1,0]+zmflux[:,0,il,0]
    if spread:
        for i in range(1,s[3]):
            zmflux[:,:,:,i]=zmflux[:,:,:,0]

    return zmflux


def divflux(vflux,spread=True):

    s=vflux.shape
    erad=6370000.
    coslats=np.cos((-90.+np.arange(s[2]+1)*180./s[2])*math.pi/180.)
    mcoslats=0.5*(coslats[1:s[2]+1]+coslats[0:s[2]])
    deltacos=coslats[1:s[2]+1]-coslats[0:s[2]]
    zmflux=np.zeros(s)
    zmflux[:,:,:,0]=np.nanmean(vflux,axis=3)
    dzmflux=np.zeros(s)
    for it in range(s[0]):
        zmflux[it,0,:,0]=zmflux[it,0,:,0]*mcoslats
    for il in range(1,s[2]):
        dzmflux[:,0,il,0]=(zmflux[:,0,il,0]-zmflux[:,0,il-1,0])/erad/math.pi*s[2]/coslats[il]
    if spread:
        for i in range(1,s[3]):
            dzmflux[:,:,:,i]=dzmflux[:,:,:,0]

    return dzmflux


def annmean(vflux):

    tmflux=np.nanmean(vflux,axis=0)
    for i in range(vflux.shape[0]):
        vflux[i,:,:,:]=tmflux

    return vflux

class TZ(datetime.tzinfo):
    def utcoffset(self, dt): return datetime.timedelta(minutes=0)     


def calc_elevangles(bgdep,lats,lons):

#    tz=pytz.timezone('Europe/London')
    day=datetime.timedelta(days=1)
    h12=datetime.timedelta(hours=12)
    elev=np.zeros([lats.shape[0],2,45000])
    for l in range(lats.shape[0]):
        loc=astral.Location(info=('London','Europe',lats[l],lons[l],'Europe/London'))
        for it in range(366):
            ltime=datetime.datetime(1900,1,1,1,0)+it*day
#            ltime=ltime+day
#               loc.lat=lats[l]
#               loc.lon=lons[l]
            elev[l,0,it]=loc.solar_elevation(ltime)
            elev[l,1,it]=loc.solar_elevation(ltime+h12)
            
#            print ltime,elev[l,:,it]

    return elev

def calc_biases(bgdep,elev):
    X=np.zeros(bgdep.shape)
#     for it in range(1,bgdep.shape[3]):

    A=np.ones([2,2*bgdep.shape[0]])
    A[0,bgdep.shape]
    c=np.linalg.lstsq(A,bgdep[0,:,0,:])

    return c



