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
# pythran export areg1(float[],float[],float)

def areg1(X,Z,alpha):

    miss_val=-1.e30


    n=Z.shape[0]
    X[0]=Z[0]
    for i in range(n-1):
        X[i+1]=alpha*X[i]+Z[i]

    return

