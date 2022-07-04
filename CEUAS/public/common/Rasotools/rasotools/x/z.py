import numpy
import math
import matplotlib.pyplot as plt
from datetime import date
from rasotools.anomaly import *
import datetime
import time

def alphaest(X):
# X must have zero mean
    rho=0.
    var=0.
    n=X.shape[0]
    m=0
    alpha=0.
    #for i in range(n-1):
        #if X[i]==X[i] and X[i+1]==X[i+1]:
            #rho+=X[i]*X[i+1]
            #var+=X[i]*X[i]
            #m=m+1
    #if var>0:
        #if m>4:
            #alpha=rho/var
        #else:
            #alpha=0.
    #else:
        #alpha=0.
    
    #if abs(alpha)>1.:
        #alpha=0.
    
    return alpha

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

## pythran export mbbfast(float[],float[],int,int,int[],float[],float[])
#def mbbfast(X,Xres,n,l,ind,Xpos,est):
##pos=fastblock_resample_mean(est,Xres,Xpos,X,ind,l,pos)
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

