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
# pythran export bootstrap_resample(float[][],float[],int[],int)

def bootstrap_resample(A,X,ind,pos):

    miss_val=-1.e30


    n=X.shape[0]
    for i in range(n):
        for j in range(n):
            A[i,j]=X[ind[pos+j]]
        pos=pos+n
    return pos

