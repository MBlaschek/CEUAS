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
## pythran export xyanomalies(float[][][][],float[][][][],float[][][][],int,int[],float[],float[],float[],int)
#
#def xyanomalies(field,fieldplus,fieldminus,syear,interval,anomaly,climatology,hilf,mlistlen):

    miss_val=-1.e30

#
#    for i in range(field.shape[2]):
#        for j in range(field.shape[3]):
#            for k in range(field.shape[0]):
#                hilf[k]=field[k,0,i,j]
#            n=janomaly(hilf,syear,interval,anomaly,climatology)
#            if n>1:
#                alpha=alphaest(anomaly)
#                if alpha<1:
#                    delta=1.64*nnanstd(anomaly)/math.sqrt(n*(1-alpha)/(1+alpha)/12*mlistlen)
#                else:
 #                   delta=0.
 #               fieldminus[0,0,i,j]-=delta
#                fieldplus[0,0,i,j]+=delta
#            else:
#                fieldminus[0,0,i,j]=numpy.nan
#                fieldplus[0,0,i,j]=numpy.nan


