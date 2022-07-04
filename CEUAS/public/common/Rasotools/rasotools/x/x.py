import numpy
import math
import matplotlib.pyplot as plt
from datetime import date
from rasotools.anomaly import *
import datetime
import time


# pythran export bubblesort(float[])
def bubblesort(X):
    for i in range(X.shape[0]):
        for j in range(1,X.shape[0]-i):
            if X[j]<X[j-1]:
                swap=X[j]
                X[j]=X[j-1]
                X[j-1]=swap
                
    return X

# pythran export plus(float[],float[])

def plus(a,b):


    c=numpy.nan

    return a+b

# pythran export tdist(float[], float[], float[], float)


def tdist(dists, lats, lons, weight):
    # copy strides from up to 8 dimensional arrays

    if dists[0] != dists[0]:
        x = numpy.cos(lats*math.pi/180.)*numpy.cos(lons*math.pi/180.)
        y = numpy.cos(lats*math.pi/180.)*numpy.sin(lons*math.pi/180.)
        z = numpy.sin(lats*math.pi/180.)

        id = 0
        for l in range(lats.shape[0]):
            for k in range(l+1, lats.shape[0]):
                dists[id] = x[l]*x[k]+y[l]*y[k]+z[l]*z[k]
                id += 1
        sdist(dists,  x,  y,  z)
        dists[:] = numpy.arccos(dists*0.999999)
        if weight != 0:
            dists[:] = numpy.exp(-dists*40./2/math.pi)*numpy.sqrt(800.)

    return

# pythran export expand(float[][][][][], int[][], int[], float[][][][][][])


def expand(b, index, pindex, a):
    # copy strides from up to 8 dimensional arrays

    for l in range(b.shape[0]):
        for k in range(b.shape[2]):
            for m in range(pindex.shape[0]):
                for j in range(b.shape[4] - 1):
                    if index[l, j + 1] == 0:
                        break
                    a[l, 0, k, m, index[l, j], index[
                        l, j + 1]] = b[l, 0, k, pindex[m], j]

    return

# pythran export sdist(float[], float[], float[], float[])


def sdist(dists, x, y, z):

    id = 0
    for l in range(x.shape[0]):
        for k in range(l, x.shape[0]):
            dists[id] = x[l]*x[k]+y[l]*y[k]+z[l]*z[k]
            id += 1

    return

