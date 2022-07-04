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
# pythran export bubblesort(float[])
def bubblesort(X):

    miss_val=-1.e30

    for i in range(X.shape[0]):
        for j in range(1,X.shape[0]-i):
            if X[j]<X[j-1]:
                swap=X[j]
                X[j]=X[j-1]
                X[j-1]=swap
                
    return X

