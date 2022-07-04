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
# pythran export grid_anomalies(float[][][][],int[][][],int,int[][][],float[][][][][],float[][][][],float[],int,int,int[],float[])

def grid_anomalies(anomaly,good,tolerance,gstatindex,ganomalies,gslopes,gst,start,stop,itime,orig):

    miss_val=-1.e30

    nj=gstatindex.shape[0]
    ni=gstatindex.shape[1]

    for j in range(nj):
        for i in range(ni):
            nh=gstatindex[j,i,0]
            for ipar in range(2):
                for ip in range(anomaly.shape[2]):
#                    if j==16 and i==10 and ip==9:
#                        print nh
                    gslopes[j,i,ipar,ip]  = numpy.nan
                    if nh>0:
                        for m in range(anomaly.shape[3]):
                            gst[m]=0
                            ganomalies[j,i,ipar,ip,m]=0.
                        for h in range(nh):
                            si=gstatindex[j,i,h+1]
#                            if good[si,ipar,ip]>(stop-start+1-tolerance)*12:
                            g=0
                            for m in range(itime.shape[0]):
                                if anomaly[si,ipar,ip,itime[m]]==anomaly[si,ipar,ip,itime[m]]:
                                    g+=1
                            if g>(stop-start+1-tolerance)*12:
                                for m in range(anomaly.shape[3]):
                                    a=anomaly[si,ipar,ip,m]
                                    if a==a:
                                        ganomalies[j,i,ipar,ip,m]+=a
                                        gst[m]+=1

                        for m in range(anomaly.shape[3]):
                            if gst[m]>0:
                                ganomalies[j,i,ipar,ip,m]/=gst[m]
                            else:
                                ganomalies[j,i,ipar,ip,m]=numpy.nan

                        tgood=0
                        for k in range(itime.shape[0]):
                            orig[k]=ganomalies[j,i,ipar,ip,itime[k]]     
                            if orig[k]==orig[k]:
                                tgood+=1
                        if tgood>(stop-start+1-tolerance)*12:
                            val  = fastlinregress(itime,orig)
#                            if val==0.0:
#                                 print j,i,ipar,ip,tgood,
                            gslopes[j,i,ipar,ip]  = val*120.

                    else:
                        for m in range(anomaly.shape[3]):
                            ganomalies[j,i,ipar,ip,m]=numpy.nan

    return
