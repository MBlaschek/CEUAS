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
# pythran export dailyanomaly(float[],int,int[],float[],float[])

def dailyanomaly( series,startyear,interval,anomaly,climatology):

    miss_val=-1.e30

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

