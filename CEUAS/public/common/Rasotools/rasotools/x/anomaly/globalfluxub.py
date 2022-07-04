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
## pythran export globalfluxub(float[],bool)
#def globalfluxub(vflux,spread):

    miss_val=-1.e30


#    s=vflux.shape
#    erad=6370000.
#    coslats=numpy.sin((-90.+numpy.arange(s[2]+1)*180./s[2])*math.pi/180.)
#    deltacos=coslats[1:s[2]+1]-coslats[0:s[2]]
#    zmflux=numpy.zeros(s)
#    for i in range(vflux.shape[0]):
#       for j in range(vflux.shape[1]):
#           for k in range(vflux.shape[2]):
#              zmflux[i,j,k,0]=nnanmean(vflux[i,j,k,:])
##    zmflux[:,:,:,0]=numpy.nanmean(vflux,axis=3)
#    for it in range(s[0]):
#        zmflux[it,0,:,0]=zmflux[it,0,:,0]*deltacos*2*math.pi*erad*erad
#    zmflux[:,0,0,0]=zmflux[:,0,1,0]/2.
#    for il in range(1,s[2]):
#        zmflux[:,0,il,0]=zmflux[:,0,il-1,0]+zmflux[:,0,il,0]
#    if spread:
#        for i in range(1,s[3]):
#            zmflux[:,:,:,i]=zmflux[:,:,:,0]
#
#    return zmflux

