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
# pythran export mbbfast(float[],float[],int,int,int[],float[],float[])

def mbbfast(X,Xres,n,l,ind,Xpos,est):

    miss_val=-1.e30

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

