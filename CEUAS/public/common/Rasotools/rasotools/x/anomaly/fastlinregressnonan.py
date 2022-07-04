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
from putils import *
# pythran export fastlinregressnonan(float[],float[])

def fastlinregressnonan(x,y):

    miss_val=-1.e30

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
        slope=numpy.nan

    return slope




#def mbb(X,n,ind=None,l=None,func=numpy.mean,timearg=None,Xpos=None):
#    """Moving Block Bootstrap as described in Mudelsee, 2011"""
#
##    t1=time.time()
#    if timearg is None:
#        Xresid=X
#    else:
#        slope=func(timearg,X)
#        Xresid=X-slope*timearg
#
#    if l==None:
#        alpha=alphaest(Xresid)
#        alpha=numpy.abs(alpha)
#        l=int((math.sqrt(6.)*alpha/(1-alpha*alpha))**(2./3.)*X.shape[0]**(1./3.))
 #       if l<1:
#            l=1
#
##    print alpha,l    
#    if(ind==None):
#        ind=numpy.empty(X.shape[0]*n/l,dtype=int)
#        for i in range(X.shape[0]*n/l):
#            ind[i]=int(random.random()*(X.shape[0]-l+1))
##        ind=numpy.random.randint(X.shape[0]-l+1,size=[X.shape[0]*n/l])

#    est=numpy.zeros(n)
#    pos=0
#    #Xres=numpy.zeros([n,X.shape[0]])
#    #for i in range(n):
#        #pos=block_resample(Xres[i,:],Xresid,ind,l,pos)
#        #if timearg is not None:
#            #est[i]=slope+func(timearg,Xres)
#        #else:
#            #est[i]=func(Xres[i,:])
#    Xres=numpy.zeros([X.shape[0]])
#    if Xpos is not None:
#        if Xpos[0,0]==0 and Xpos[n-1,X.shape[0]-1]==0:
#            Xpos=numpy.zeros([n,X.shape[0]],dtype=int)
#    else:
#        Xpos=numpy.zeros([n,X.shape[0]],dtype=int)

#    pos=fastblock_resample_mean(est,Xres,Xpos,Xresid,ind,l,pos)
##    print 'mbb: ',time.time()-t1

##    if numpy.nanstd(Xresid)>0:
##        print est
#    return est,Xpos

