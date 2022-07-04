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
# pythran export anomalies_and_slopes(float[][][][],int,int[],int,int,int,float[],float[],float[][][][],float[],int[][][],float[][][][])

def anomalies_and_slopes(series,startyear,interval,tolerance,iens,itime,orig,anomaly,anomalies,climatology,good,slopes):

    miss_val=-1.e30


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

