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
# pythran export annmean(float[][][][])
def annmean(vflux):

    miss_val=-1.e30


    for k in range(vflux.shape[1]):
       for l in range(vflux.shape[2]):
          for m in range(vflux.shape[3]):
              c=0
              s=0.
              for j in range(vflux.shape[0]):
                 if vflux(j,k,l,m)==vflux(j,k,l,m):
                     s+=vflux(j,k,l,m)
                     c+=1
              if c>0:
                  vflux[:,k,l,m]=s/c
              else:
                  vflux[:,k,l,m]=numpy.nan
                         

    return vflux




