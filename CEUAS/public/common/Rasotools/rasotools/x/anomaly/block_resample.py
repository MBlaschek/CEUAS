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
# pythran export block_resample(float[],float[],int[],int,int)

def block_resample(A,X,ind,l,pos):

    miss_val=-1.e30


    n=X.shape[0]
    jpos=0
    for j in range((n-l)/l):
        for i in range(l):
            A[jpos+i]=X[ind[pos+j]+i]
        jpos+=l
    pos=pos+(n-l)/l
    return pos

