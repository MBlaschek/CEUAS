#from numba import *
import numpy
import matplotlib.pyplot as plt
from scipy import stats
#import utils
import math
import os, sys, inspect
import astral
import datetime
import random

from .putils import *


# pythran export anomalies_and_slopes(float[][][][],int,int[],int,int,int[],float[],float[],float[][][][],float[],int[][][],float[][][][])

def anomalies_and_slopes(series,startyear,interval,tolerance,iens,itime,orig,anomaly,anomalies,climatology,good,slopes):

    sshape=series.shape
    stop=interval[1]
    start=interval[0]
    ##omp parallel for 
    for si in range(sshape[0]):
        for ipar in range(sshape[1]):
            for ip in range(sshape[2]):
                clim=numpy.zeros(12)
                good[si,ipar,ip]=janomaly(series[si,ipar,ip,:],startyear,interval,anomalies[si,ipar,ip,:],clim)

                slopes[si,iens,ipar,ip]  = numpy.nan
                if good[si,ipar,ip]>(stop-start+1-tolerance)*12:
                    tgood=0
                    for j in range(itime.shape[0]):
                        if anomalies[si,ipar,ip,itime[j]]==anomalies[si,ipar,ip,itime[j]]:
                            tgood+=1
                    if tgood<(stop-start+1-tolerance)*12:
                        anomalies[si,ipar,ip,:]=numpy.nan
                else:
                    anomalies[si,ipar,ip,:]=numpy.nan
    #omp parallel for 
    for si in range(sshape[0]):
        for ipar in range(sshape[1]):
            for ip in range(sshape[2]):
                ftime=numpy.asarray(itime,dtype=float)
                anom=numpy.zeros(itime.shape[0])
                for j in range(itime.shape[0]):
                    anom[j]=anomalies[si,ipar,ip,itime[j]]
                val = fastlinregress(ftime,anom)*120.
                slopes[si,iens,ipar,ip]=val


    return

# pythran export janomaly(float[],int,int[],float[],float[])

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
                climatology[j]=climatology[j]+s
                l=l+1
        if l > (stop-start+1)/2:
            climatology[j]/=l
            for k in range(n/12):
                anomaly[k*12+j]=series[k*12+j]-climatology[j]
                if anomaly[k*12+j]==anomaly[k*12+j]: 
                    good+=1
        else:
            climatology[j]=numpy.nan
            for k in range(n/12):
                anomaly[k*12+j]=numpy.nan


    return good

# pythran export nnanstd(float[])

def nnanstd(field):

    s=0.
    ss=0.
    c=0
    for i in range(field.shape[0]):
        if field[i]==field[i]:
            s = s + field[i]
            ss = ss + field[i]*field[i]
            c=c+1
    return math.sqrt(ss/c-s*s/c/c)

# pythran export nnanmean(float[])

def nnanmean(field):

    s=0.
    c=0
    for i in range(field.shape[0]):
        if field[i]==field[i]:
            s+=field[i]
            c=c+1
    if c>0:
        s/=c
    else:
        s=numpy.nan

    return s

# pythran export nnanmeanlist0(float[][][][],float[][][][],int[][][],int list)

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
                    fieldm[0,j,k,l]=numpy.nan

    return

# pythran export xyanomalies(float[][][][],float[][][][],float[][][][],int,int[],float[],float[],float[],int)

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
                fieldminus[0,0,i,j]=numpy.nan
                fieldplus[0,0,i,j]=numpy.nan


## pythran export xyanomaliesmbb(float[][][][],float[][][][],float[][][][],int,int[],float[],float[],float[],int[])

#def xyanomaliesmbb(field,fieldplus,fieldminus,syear,interval,anomaly,climatology,hilf,season):

    #n=40
    #ind=numpy.empty(field.shape[0]*n,dtype='int')
    #for i in range(field.shape[0]*n):
        #ind[i]=int(numpy.floor(random.random()*field.shape[0]))
##    ind=numpy.random.randint(field.shape[0],size=[field.shape[0]*n])
    #Xpos=numpy.zeros([n,field.shape[0]],dtype='int')
    #Xres=numpy.zeros(field.shape[0])
    #est=numpy.zeros(n)
    #hilf[:]=numpy.nan
    #xyanomaliesmbbcore(field,fieldplus,fieldminus,syear,interval,anomaly,climatology,n,ind,Xpos,Xres,est,hilf,season)

    #return

# pythran export bubblesort(float[])
def bubblesort(X):
    for i in range(X.shape[0]):
        for j in range(1,X.shape[0]-i):
            if X[j]<X[j-1]:
                swap=X[j]
                X[j]=X[j-1]
                X[j-1]=swap

    return X

# pythran export xyanomaliesmbbcore(float[][][][],float[][][][],float[][][][],int,int[],float[],float[],int,int,int,int,float[],float[],int[])
def xyanomaliesmbbcore(field,fieldplus,fieldminus,syear,interval,anomaly,climatology,n,ind,Xpos,Xres,est,hilf,season):

    for i in range(field.shape[2]):
#       print i
        for j in range(field.shape[3]):
            for k in range(field.shape[0]/12):
                for m in range(season.shape[0]):
                    hilf[k*12+season[m]]=field[k*12+season[m],0,i,j]
            l=-1
            janomaly(hilf,syear,interval,anomaly,climatology)
#            mbbfast(anomaly,Xres,n,l,ind,Xpos,est)
#           est,Xpos=utils.mbb(anomaly,n,func=nnanmean,ind=ind,Xpos=Xpos)
#           sest=sorted(est)
            bubblesort(est)
            fieldminus[0,0,i,j]+=est[int(numpy.floor(0.05*n))]
            fieldplus[0,0,i,j]+=est[int(numpy.floor(0.95*n))]

    return

def dailyanomalydriver(d,ens,plotproperties):
    h=d["ddata"]
    sh=h.shape
    anomaly=numpy.zeros(sh[4])
    climatology=numpy.zeros(365)
    for ist in range(sh[0]):
        for ipar in range(sh[2]):
            for ip in range(sh[3]):
                good=dailyanomaly(h[ist,ens,ipar,ip,:],int(d["startdate"]/10000),
                                  numpy.array(plotproperties["plotinterval"]),anomaly,climatology)
                h[ist,ens,ipar,ip,:]=anomaly

    return h[:,ens,:,:,:]

# pythran export dailyanomaly(float[],int,int[],float[],float[])

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
            climatology[j]=numpy.nan
            for k in range(n/365):
                anomaly[ik]=numpy.nan

        if good>0:
            pass


    return good

# pythran export grid_anomalies(float[][][][],int[][][],int,int[][][],float[][][][][],float[][][][],float[],int,int,int[],float[])

def grid_anomalies(anomaly,good,tolerance,gstatindex,ganomalies,gslopes,gst,start,stop,itime,orig):
    nj=gstatindex.shape[0]
    ni=gstatindex.shape[1]

    for j in range(nj):
        for i in range(ni):
            nh=gstatindex[j,i,0]
            for ipar in range(2):
                for ip in range(anomaly.shape[2]):
#                    if j==16 and i==10 and ip==9:
#                        print nh
                    gslopes[j,i,ipar,ip]  = numpy.nan
                    if nh>0:
                        for m in range(anomaly.shape[3]):
                            gst[m]=0
                            ganomalies[j,i,ipar,ip,m]=0.
                        for h in range(nh):
                            si=gstatindex[j,i,h+1]
#                            if good[si,ipar,ip]>(stop-start+1-tolerance)*12:
                            g=0
                            for m in range(itime.shape[0]):
                                if anomaly[si,ipar,ip,itime[m]]==anomaly[si,ipar,ip,itime[m]]:
                                    g+=1
                            if g>(stop-start+1-tolerance)*12:
                                for m in range(anomaly.shape[3]):
                                    a=anomaly[si,ipar,ip,m]
                                    if a==a:
                                        ganomalies[j,i,ipar,ip,m]+=a
                                        gst[m]+=1

                        for m in range(anomaly.shape[3]):
                            if gst[m]>0:
                                ganomalies[j,i,ipar,ip,m]/=gst[m]
                            else:
                                ganomalies[j,i,ipar,ip,m]=numpy.nan

                        tgood=0
                        for k in range(itime.shape[0]):
                            orig[k]=ganomalies[j,i,ipar,ip,itime[k]]     
                            if orig[k]==orig[k]:
                                tgood+=1
                        if tgood>(stop-start+1-tolerance)*12:
                            val  = fastlinregress(itime,orig)
#                            if val==0.0:
#                                 print j,i,ipar,ip,tgood,
                            gslopes[j,i,ipar,ip]  = val*120.

                    else:
                        for m in range(anomaly.shape[3]):
                            ganomalies[j,i,ipar,ip,m]=numpy.nan

    return

# pythran export fastlinregress(float[],float[])

def fastlinregress(x,y):
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
        slope=numpy.nan

    return float(slope)

# pythran export fastlinregressnonan(float[],float[])

def fastlinregressnonan(x,y):
    m=x.shape[0]

    xmean=0.
    ymean=0.
    xcov=0.
    xsq=0.
    if m>2:   
        for k in range(m):
            xmean = xmean+x[k]
            ymean = ymean + y[k]
        xmean = xmean/m
        ymean = ymean/m
        for k in range(m):
            xcov=xcov+x[k]*y[k]
            xsq=xsq+x[k]*x[k]
        slope=(xcov-m*xmean*ymean)/(xsq-m*xmean*xmean)
    else:
        slope=numpy.nan

    return slope


#def mbb(X,n,ind=None,l=None,func=numpy.mean,timearg=None,Xpos=None):
    #"""Moving Block Bootstrap as described in Mudelsee, 2011"""
##
###    t1=time.time()
    #if timearg is None:
        #Xresid=X
    #else:
        #slope=func(timearg,X)
        #Xresid=X-slope*timearg
##
    #if l==None:
        #alpha=alphaest(Xresid)
        #alpha=numpy.abs(alpha)
        #l=int((math.sqrt(6.)*alpha/(1-alpha*alpha))**(2./3.)*X.shape[0]**(1./3.))
        #if l<1:
            #l=1
##
###    print alpha,l    
    #if(ind==None):
        #ind=numpy.empty(X.shape[0]*n/l,dtype=int)
        #for i in range(X.shape[0]*n/l):
            #ind[i]=int(random.random()*(X.shape[0]-l+1))
###        ind=numpy.random.randint(X.shape[0]-l+1,size=[X.shape[0]*n/l])

    #est=numpy.zeros(n)
    #pos=0
    ##Xres=numpy.zeros([n,X.shape[0]])
    ##for i in range(n):
        ##pos=block_resample(Xres[i,:],Xresid,ind,l,pos)
        ##if timearg is not None:
            ##est[i]=slope+func(timearg,Xres)
        ##else:
            ##est[i]=func(Xres[i,:])
    #Xres=numpy.zeros([X.shape[0]])
    #if Xpos is not None:
        #if Xpos[0,0]==0 and Xpos[n-1,X.shape[0]-1]==0:
            #Xpos=numpy.zeros([n,X.shape[0]],dtype=int)
    #else:
        #Xpos=numpy.zeros([n,X.shape[0]],dtype=int)

    #pos=fastblock_resample_mean(est,Xres,Xpos,Xresid,ind,l,pos)
###    print 'mbb: ',time.time()-t1

###    if numpy.nanstd(Xresid)>0:
###        print est
    #return est,Xpos

# pythran export mbbfast(float[],float[],int,int,int[],int[][],float[])

def mbbfast(X,Xres,n,l,ind,Xpos,est):
        """Moving Block Bootstrap as described in Mudelsee, 2011"""

        if l==-1:
            alpha=alphaest(X)
            alpha=numpy.abs(alpha)
            l=int((math.sqrt(6.)*alpha/(1-alpha*alpha))**(2./3.)*X.shape[0]**(1./3.))
            if l<1:
                l=1
        pos=0

        pos=fastblock_resample_mean(est,Xres,Xpos,X,ind,l,pos)

        return

## pythran export areg1(float[],float[],float)

def areg1(X,Z,alpha):

        n=Z.shape[0]
        X[0]=Z[0]
        for i in range(n-1):
            X[i+1]=alpha*X[i]+Z[i]

        return

# pythran export alphaest(float[])

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

# pythran export bootstrap_resample(float[][],float[],int[],int)

def bootstrap_resample(A,X,ind,pos):

    n=X.shape[0]
    for i in range(n):
        for j in range(n):
            A[i,j]=X[ind[pos+j]]
        pos=pos+n
    return pos

# pythran export block_resample(float[],float[],int[],int,int)

def block_resample(A,X,ind,l,pos):

    n=X.shape[0]
    jpos=0
    for j in range((n-l)/l):
        for i in range(l):
            A[jpos+i]=X[ind[pos+j]+i]
        jpos+=l
    pos=pos+(n-l)/l
    return pos

# pythran export fastblock_resample_mean(float[],float[],int[][],float[],int[],int,int)

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

## pythran export globalflux(float[])
#def globalflux(vflux):

    #s=vflux.shape
    #zmflux=globalfluxub(vflux,spread=False)
    #for il in range(1,s[2]):
        #zmflux[:,0,il,0]-=zmflux[:,0,s[2]-1,0]*il/s[2]
    #for i in range(1,s[3]):
        #zmflux[:,:,:,i]=zmflux[:,:,:,0]
##
    #return zmflux

## pythran export globalfluxub(float[],bool)
#def globalfluxub(vflux,spread):

    #s=vflux.shape
    #erad=6370000.
    #coslats=numpy.sin((-90.+numpy.arange(s[2]+1)*180./s[2])*math.pi/180.)
    #deltacos=coslats[1:s[2]+1]-coslats[0:s[2]]
    #zmflux=numpy.zeros(s)
    #for i in range(vflux.shape[0]):
        #for j in range(vflux.shape[1]):
            #for k in range(vflux.shape[2]):
                #zmflux[i,j,k,0]=nnanmean(vflux[i,j,k,:])
###    zmflux[:,:,:,0]=numpy.nanmean(vflux,axis=3)
    #for it in range(s[0]):
        #zmflux[it,0,:,0]=zmflux[it,0,:,0]*deltacos*2*math.pi*erad*erad
    #zmflux[:,0,0,0]=zmflux[:,0,1,0]/2.
    #for il in range(1,s[2]):
        #zmflux[:,0,il,0]=zmflux[:,0,il-1,0]+zmflux[:,0,il,0]
    #if spread:
        #for i in range(1,s[3]):
            #zmflux[:,:,:,i]=zmflux[:,:,:,0]
##
    #return zmflux

# pythran export annmean(float[][][][])
def annmean(vflux):

    for k in range(vflux.shape[1]):
        for l in range(vflux.shape[2]):
            for m in range(vflux.shape[3]):
                c=0
                s=0.
                for j in range(vflux.shape[0]):
                    if vflux(j,k,l,m)==vflux(j,k,l,m):
                        s=s+vflux(j,k,l,m)
                        c=c+1
                if c>0:
                    vflux[:,k,l,m]=s/c
                else:
                    vflux[:,k,l,m]=numpy.nan


    return vflux

# pythran export ndnanmean(float[][][][][])

def ndnanmean(field):

    fieldm=numpy.zeros((field.shape[0],field.shape[3],field.shape[4]))
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
                if c>0:
                    fieldm[l,k,j]=s/c
                else:
                    fieldm[l,k,j]=numpy.nan

    return fieldm



