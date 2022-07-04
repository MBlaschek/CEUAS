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
# pythran export nnanmeanlist0(float[][][][],float[][][][],int[][][],int list)

def nnanmeanlist0(field,fieldm,counts,index):

    miss_val=-1.e30


    for i in index:
        for j in range(field.shape[1]):
            for k in range(field.shape[2]):
                for l in range(field.shape[3]):
                    f=field[i,j,k,l]
                    if f==f:
                        fieldm[0,j,k,l]+=f
                        counts[j,k,l]+=1
    for j in range(field.shape[1]):
        for k in range(field.shape[2]):
            for l in range(field.shape[3]):
                if counts[j,k,l]>0:
                    fieldm[0,j,k,l]/=counts[j,k,l]
                else:
                    fieldm[0,j,k,l]=numpy.nan

    return

