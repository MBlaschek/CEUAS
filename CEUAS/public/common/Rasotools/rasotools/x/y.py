import numpy
import math
import matplotlib.pyplot as plt
from datetime import date
from rasotools.anomaly import *
import datetime
import time

# pythran export anomalies_and_slopes(float[][][][],int,int[],int,int,int[],float[],float[],float[][][][],float[],int[][][],float[][][][])

def anomalies_and_slopes(series,startyear,interval,tolerance,iens,itime,orig,anomaly,anomalies,climatology,good,slopes):

    sshape=series.shape
    stop=interval[1]
    start=interval[0]
    for si in range(sshape[0]):
        for ipar in range(sshape[1]):
            for ip in range(sshape[2]):
                for j in range(sshape[3]):
                    orig[j]=series[si,ipar,ip,j]
                good[si,ipar,ip]=janomaly(orig,startyear,interval,anomaly,climatology)

                slopes[si,iens,ipar,ip]  = numpy.nan
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
                        anomalies[si,ipar,ip,j]=numpy.nan



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

    return slope


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


# pythran export bubblesort(float[])
def bubblesort(X):
    for i in range(X.shape[0]):
        for j in range(1,X.shape[0]-i):
            if X[j]<X[j-1]:
                swap=X[j]
                X[j]=X[j-1]
                X[j-1]=swap
                
    return X

