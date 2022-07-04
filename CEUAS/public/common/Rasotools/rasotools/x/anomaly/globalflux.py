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
## pythran export globalflux(float[])
#def globalflux(vflux):

    miss_val=-1.e30


#    s=vflux.shape
#    zmflux=globalfluxub(vflux,spread=False)
#    for il in range(1,s[2]):
#        zmflux[:,0,il,0]-=zmflux[:,0,s[2]-1,0]*il/s[2]
#    for i in range(1,s[3]):
#        zmflux[:,:,:,i]=zmflux[:,:,:,0]
#
#    return zmflux

