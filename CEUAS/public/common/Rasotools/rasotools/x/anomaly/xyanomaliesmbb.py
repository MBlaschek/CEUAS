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

# pythran export xyanomaliesmbb(float[][][][],float[][][][],float[][][][],int,int[],float[],float[],float[],int[])

def xyanomaliesmbb(field,fieldplus,fieldminus,syear,interval,anomaly,climatology,hilf,season):

    miss_val=-1.e30


    n=40
    ind=numpy.zeros(field.shape[0]*n,'int')
    for i in range(field.shape[0]*n):
        ind[i]=int(random.random()*field.shape[0])
#    ind=numpy.random.randint(field.shape[0],size=[field.shape[0]*n])
    Xpos=numpy.zeros([n,field.shape[0]],'int')
    Xres=numpy.zeros(field.shape[0])
    est=numpy.zeros(n)
    hilf[:]=numpy.nan
#    xyanomaliesmbbcore(field,fieldplus,fieldminus,syear,interval,anomaly,climatology,n,ind,Xpos,Xres,est,hilf,season)

    return

