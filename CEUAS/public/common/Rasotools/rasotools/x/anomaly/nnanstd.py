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
# pythran export nnanstd(float[])

def nnanstd(field):

    miss_val=-1.e30


    s=0.
    ss=0.
    c=0
    for i in range(field.shape[0]):
        if field[i]==field[i]:
            s+=field[i]
            ss+=field[i]*field[i]
            c=c+1
    return math.sqrt(ss/c-s*s/c/c)

