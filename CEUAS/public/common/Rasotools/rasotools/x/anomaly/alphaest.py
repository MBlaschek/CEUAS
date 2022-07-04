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
# pythran export alphaest(float[])

def alphaest(X):

    miss_val=-1.e30

# X must have zero mean
    rho=0.
    var=0.
    n=X.shape[0]
    m=0
    for i in range(n-1):
        if X[i]==X[i] and X[i+1]==X[i+1]:
            rho+=X[i]*X[i+1]
            var+=X[i]*X[i]
            m=m+1
    if var>0:
        if m>4:
            alpha=rho/var
        else:
            alpha=0.
    else:
        alpha=0.

    if abs(alpha)>1:
        alpha=0.

    return alpha

