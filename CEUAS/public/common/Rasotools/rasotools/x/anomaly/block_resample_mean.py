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
# pythran export block_resample_mean(float[],float[],float[][],float[][],int[],int,int)

def fastblock_resample_mean(est,Y,A,X,ind,l,pos):

    miss_val=-1.e30


    n=A.shape[1]
    if A[0,0]==0 and A[A.shape[0]-1,n-1]==0:
        for k in range(A.shape[0]):
            jpos=0
            for j in range((n-l)/l):
                for i in range(l):
                    A[k,jpos+i]=ind[pos+j]+i
                jpos+=l
            pos=pos+(n-l)/l
    for k in range(A.shape[0]):
        for j in range(n):
            Y[j]=X[A[k,j]]
        est[k]=nnanmean(Y)


    return pos
