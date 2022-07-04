from numba import *
import numpy
import matplotlib.pyplot as plt
from scipy import stats
#import utils
import math
import os, sys, inspect
import astral
import datetime
from utils import *
#   import numpy
#   from scipy.stats.stats import nanmean
# pythran export manomaly(float[],float[],float[])
def manomaly( series,anomaly,climatology):

    miss_val=-1.e30

    interval=numpy.int32([1980,2009])
    startyear=numpy.int32(1957)
    n=series.shape[0]
    for j in range(12):
        hanomaly=series[j:n:12]

        climatology[j]=sum(hanomaly[interval[0]-startyear:interval[1]-startyear])/(n/12)
        anomaly[j:n:12]=hanomaly-climatology[j]
    return

