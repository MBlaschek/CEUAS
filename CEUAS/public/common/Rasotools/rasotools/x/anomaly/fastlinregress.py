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
# pythran export fastlinregress(float[],float[])

def fastlinregress(x,y):

    miss_val=-1.e30

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

