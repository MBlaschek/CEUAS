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
# pythran export janomaly(float[],int,int[],float[],float[])

def janomaly( series,startyear,interval,anomaly,climatology):

    miss_val=-1.e30

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

