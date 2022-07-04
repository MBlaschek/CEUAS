import numpy
import math
import matplotlib.pyplot as plt
from datetime import date
from rasotools.anomaly import *
import datetime
import time


# pythran export hmav13(float[],float[])

def hmav13(x,y):

    hmweights=numpy.asarray([
        [-0.01935,-0.02787,0,0.06549,0.14736,0.21434,0.24006,0.21434,0.14736,0.06549,0,-0.02787,-0.01935],
        [0,-0.01643,-0.02577,0.00127,0.06594,0.14698,0.21314,0.23803,0.21149,0.14368,0.06099,-0.00533,-0.03401],
        [0,0,-0.011,-0.02204,0.0033,0.06626,0.14559,0.21004,0.23324,0.20498,0.13547,0.05108,-0.01695],
        [0,0,0,-0.00814,-0.0202,0.00413,0.06608,0.14441,0.20784,0.23002,0.20076,0.13024,0.04483],
        [0,0,0,0,-0.01604,-0.02487,0.00267,0.06784,0.14939,0.21605,0.24144,0.2154,0.1481],
        [0,0,0,0,0,-0.04271,-0.03864,0.00182,0.0799,0.17436,0.25392,0.29223,0.2791],
        [0,0,0,0,0,0,-0.09187,-0.05812,0.01202,0.11977,0.2439,0.35315,0.42113]
    ])    

    for i in range(6):
        if x[i]==x[i]:
            y[i]=0.
            wsum=0.
            for j in range(5-i,6):
                if x[i+j]==x[i+j]:
                    y[i]+=x[i+j]*hmweights[6-i,j+6]
                    wsum+=hmweights[6-i,j+6]
            if wsum>0:
                y[i]/=wsum
            else:
                y[i]=numpy.nan
        else:
            y[i]=numpy.nan 

    last=x.shape[0]
    for i in range(last-1,6,-1):
        if x[i]!=x[i]:
            y[i]=numpy.nan
            last-=1
        else:
            break

    for i in range(6,last-6):
        y[i]=0.
        wsum=0.
        if x[i]==x[i]:
            for j in range(-6,6):
                if x[i+j]==x[i+j]:
                    y[i]+=x[i+j]*hmweights[0,j+6]
                    wsum+=hmweights[0,j+6]
            if wsum>0:
                y[i]/=wsum
            else:
                y[i]=numpy.nan
        else:
            y[i]=numpy.nan 

    for i in range(last-6,last):
        if x[i]==x[i]:
            wsum=0.
            y[i]=0.
            for j in range(-6,last-i):
                if x[i+j]==x[i+j]:
                    y[i]+=x[i+j]*hmweights[i-last+7,6-j]
                    wsum+=hmweights[i-last+7,6-j]
            if wsum>0:
                y[i]/=wsum
            else:
                y[i]=numpy.nan
        else:
            y[i]=numpy.nan 

    return y



# pythran export stationaverage(float[][][][],float[][][],int[][][], int)
def stationaverage(currentdata,currentdataav,avcount,thresh):

    s=currentdata.shape
    for istat in range(s[0]):
        for ipar in range(s[1]):
            for ip in range(s[2]):
                for it in range(s[3]):
                    if currentdata[istat,ipar,ip,it]==currentdata[istat,ipar,ip,it]:
                        currentdataav[ipar,ip,it]=currentdataav[ipar,ip,it]+currentdata[istat,ipar,ip,it]
                        avcount[ipar,ip,it]=avcount[ipar,ip,it]+1

    for ipar in range(s[1]):
        for ip in range(s[2]):
            for it in range(s[3]):
                if avcount[ipar,ip,it]>=thresh:
                    currentdataav[ipar,ip,it]=currentdataav[ipar,ip,it]/avcount[ipar,ip,it]
                else:
                    currentdataav[ipar,ip,it]=numpy.nan

# pythran export statcore(float[], float[], float[], int)


def statcore(feld, arr, weights, dim):
    s = 0.
    sq = 0.
    count = 0.

    if weights.shape[0] != 0:

        stride = 1
        for i in range(feld.ndim - 1, dim, -1):
            stride *= feld.shape[i]

        k = 0
        ki = 0
        for i in range(arr.shape[0]):
            if arr[i] == arr[i]:
                s += arr[i] * weights[k]
                sq += arr[i] * arr[i] * weights[k]
                count += weights[k]
            ki += 1
            if ki == stride:
                k += 1
                ki = 0
                if k == weights.shape[0]:
                    k = 0

    else:
        for i in range(arr.shape[0]):
            if arr[i] == arr[i]:
                s += arr[i]
                sq += arr[i] * arr[i]
                count += 1.0

    if count > 0:
        s /= count
        sq /= count

    return s, sq

# pythran export had_rebin_pythran_3672_to_1836(float32[][][], float32[][][][], int, int, int, float32, float32)


def had_rebin_pythran_3672_to_1836(mem, hadens, ens, startyear,
                                   endyear, nc_miss_val, miss_val):

    index = hadens.shape[0] - hadens.shape[0] - 1
    for itime in range((startyear - 1850) * 12, (endyear - 1850 + 1) * 12):
        index += 1
        for ilat in range(hadens.shape[2]):
            for ilon in range(hadens.shape[3]):
                if ilon < 18:
                    ishift = 18
                else:
                    ishift = -18
                sum = nc_miss_val - nc_miss_val
                n = nc_miss_val - nc_miss_val
                jj = (ilon + ishift) * 2
                for i in range(2):
                    ii = ilat * 2 + i
                    for j in range(2):
                        x = mem[itime, ii, jj + j]
                        if(x != nc_miss_val):
                            sum += x
                            n += 1.0
                if n > 0:
                    hadens[ens, index, ilat, ilon] = sum / n
                else:
                    hadens[ens, index, ilat, ilon] = miss_val

    return


def getindex(alldates, somedates, index):

    n = somedates.shape[0]
    m = alldates.shape[0]
    index[0] = 0
    jsave = 0
    for j in range(n):
        for k in range(index[jsave], m):
            if alldates[k] == somedates[j]:
                index[j] = k
                jsave = j
                break

    return


def find_gstatindex(glons, glats, lons, lats, gstatindex):

    for l in range(lats.shape[0]):
        ilon = numpy.floor(lons[l] / (360. / glons.shape[0]))
        ilat = numpy.floor((lats[l] + 90.) / (180. / glats.shape[0]))
        gstatindex[ilat, ilon, 0] += 1
        gstatindex[ilat, ilon, gstatindex[ilat, ilon, 0]] = l

    return


def getindex2(alldates, somedates, index):

    n = somedates.shape[0]
    m = alldates.shape[0]
    index[0] = 0
    jsave = 0
    for j in range(n):
        for k in range(index[jsave], m):
            if alldates[k] == somedates[j]:
                index[j] = k
                jsave = j
                break
            if alldates[k] > somedates[j]:
                index[j] = k
                jsave = j + 1
                break
    if index[jsave] == 0:
        index[jsave] = m - 1

    return


# pythran export copystride(float[][][][][], float[][][], int[], int, int, int[], float)


def copystride(a, b, index, n, m, pindex, fmiss_val):
    # expand time dimension,  reduce p dimension

    for l in range(b.shape[0]):
        for k in range(pindex.shape[0]):  # reduce p dimension
            for j in range(b.shape[2]):
                if b[l, pindex[k], j] !=  fmiss_val:
                    a[n, m, l, k, index[j]] = b[l, pindex[k], j]

    return

# pythran export expandandadd(float[][][][][], float[][][][], int[][][], int[], int, float[][][][][], float)


def expandandadd(b, ref, index, pindex, iens, a, sign):
    # copy strides from up to 8 dimensional arrays
    pm = pindex[0]
    for l in range(b.shape[0]):
        for k in range(b.shape[2]):
            for m in range(pindex.shape[0]):
                if sign == 0.0:
                    for j in range(b.shape[4]):
                        a[l, k, m, j] = b[l, iens, k, pindex[m], j]
                else:

                    if ref.shape[2] != a.shape[2]:
                        pm = pindex[m]+m-m
                    else:
                        pm = pindex[m]+m-pindex[m]
                    for j in range(b.shape[4]-1):
                        if index[l, iens, j+1] == 0:
                            if j == 0:
                                index[l, iens, j+1] = a.shape[3]
                                for i in range(index[l, iens, j+1]):
                                    a[l, k, m, i] = ref[l, k, pm, i]
                                index[l, iens, j+1] = 0
                            else:
                                index[l, iens, j+1] = a.shape[3]
                                for i in range(index[l, iens, j],
                                               index[l, iens, j+1]):
                                    a[l, k, m, i] = ref[l, k, pm, i]
                                index[l, iens, j+1] = 0
                            break
                        else:
                            for i in range(index[l, iens, j],
                                           index[l, iens, j+1]):
                                if ref[l, k, pm, i] == ref[l, k, pm, i]:
                                    a[l, k, m, i] = ref[l, k, pm, i] + \
                                        sign*b[l, iens, k, pindex[m], j]
                                else:
                                    a[l, k, m, i] = numpy.nan

    return


def tcost(dists, slopes, cost):
    # calculate trend heterogeneity cost function

    id = 0
    goodstats = 0
    tcost = 0.
    for l in range(slopes.shape[0]):
        cost[l] = 0.
    for l in range(slopes.shape[0]):
        if slopes[l] == slopes[l]:
            goodstats += 1
            for k in range(l, slopes.shape[0]):
                if slopes[l] == slopes[l] and slopes[k] == slopes[k]:
                    s = (slopes[l]-slopes[k])*dists[id]
                    cost[l] += s*s
                    cost[k] += s*s
                    if l == 24 and k == 29:
                        print dists[id], slopes[l], slopes[k]
                id += 1
            if goodstats > 0:
                tcost += cost[l]
        else:
            id += slopes.shape[0]-l

    if goodstats > 0:
        tcost /= goodstats
        for l in range(slopes.shape[0]):
            cost[l] /= goodstats

    return tcost


def thin(t, index, n):

    ni = t.shape[0]
    if n < 2:
        for i in range(t.shape[0]):
            index[i] = i
    else:
        ni = t.shape[0]/n
        for i in range(ni):
            for j in range(n):
                index[i] = i*n
                if t[i*n+j] == t[i*n+j]:
                    index[i] = i*n+j
                    break

    return ni


def rmean(t, tret, tmean, index, runmean):

    ni = t.shape[0]
    good = runmean-runmean
    if runmean < 2:
        for i in range(ni):
            tret[i] = t[i]
    else:

        for j in range(ni):
            if t[j] == t[j]:
                index[good] = j
                good += 1
            else:
                tret[j] = t[j]

        if good > runmean+2:
            i = runmean/2
            if runmean % 2 == 1:
                tmean[i] = 0.
                for k in range(-runmean/2+1, runmean/2):
                    tmean[i] += t[index[i+k]]
                tmean[i] /= runmean

                for i in range(runmean/2+1, good-runmean/2-1):
                    tmean[i] = (tmean[i-1]*runmean +
                                t[index[i+runmean/2]] -
                                t[index[i-runmean/2]])/runmean

            else:

                i = runmean/2
                tmean[i] = 0.
                for k in range(-runmean/2, runmean/2-1):
                    tmean[i] += t[index[i+k]]
                tmean[i] /= runmean

                for i in range(runmean/2+1, good-runmean/2-1):
                    tmean[i] = (tmean[i-1]*runmean +
                                t[index[i+runmean/2-1]] -
                                t[index[i-runmean/2]])/runmean

            for i in range(good):
                tret[index[i]] = tmean[i]
        else:
            for i in range(good):
                tret[index[i]] = t[index[i]]

    return tret


def snhtmov(t, tsa, snhtparas, index, count, tmean, tsquare):

    n = snhtparas[0]
    max_miss = snhtparas[1]
    ninc = snhtparas[2]

    ni = t.shape[0]
    good = 0
    for j in range(ni):
        if t[j] == t[j]:
            index[good] = j
            good += 1
        count[j] = good

    if good > n-2*max_miss:
        rm = n/2
        i = rm/2
        for k in range(n/2, ni, ninc):
            if count[k]-count[k-n/2] > n/2-max_miss:
                tmean[k] = 0.
                tsquare[k] = 0.
                for i in range(count[k-n/2], count[k]):
                    x = t[index[i]]
                    tmean[k] += x
                    tsquare[k] += x*x
                tmean[k] /= count[k]-count[k-n/2]
                tsquare[k] /= count[k]-count[k-n/2]
        fak = n/2/ninc  # approximate length of SNHT window with fak*ninc
        for k in range(n/2, ni, ninc):
            if count[k]-count[k-n/2] > n/2-max_miss and \
               count[k+fak*ninc]-count[k] > n/2-max_miss:
                m = (tmean[k]+tmean[k+fak*ninc])/2.
                tdiv = tsquare[k]-tmean[k]*tmean[k]
                if tdiv > 0:
                    tsa[k] = n/2*((tmean[k]-m)*(tmean[k]-m) +
                                  (tmean[k+fak*ninc]-m) *
                                  (tmean[k+fak*ninc]-m))/math.sqrt(tdiv)
                else:
                    tsa[k] = 0.
                for i in range(ninc):
                    tsa[k+i] = tsa[k]

    return


def snhtmov2(t, tsa, snhtparas, index, count, tmean, tsquare):

    n = snhtparas[0]
    max_miss = snhtparas[1]
    ninc = snhtparas[2]

    ni = t.shape[0]
    good = 0
    for j in range(ni):
        tmean[j] = numpy.nan
        tsquare[j] = numpy.nan
        tsa[j] = numpy.nan
    tmean[0] = 0.
    tsquare[0] = 0.
    for j in range(ni):
        count[j] = 0
        if t[j] == t[j]:
            index[good] = j
            if good > 0:
                tmean[good] = tmean[good-1]+t[j]
                tsquare[good] = tsquare[good-1]+t[j]*t[j]
            else:
                tmean[good] = t[j]
                tsquare[good] = t[j]*t[j]
            good += 1
        if good > 0:
            count[j] = good-1

    if good > n-2*max_miss:
        rm = n/2
        for k in range(rm-max_miss, ni-(rm-max_miss)):
            xm = k-rm
            if xm < 0:
                xm = 0
            xp = k+rm
            if xp > ni-1:
                xp = ni-1
            if count[k]-count[xm] > rm-max_miss and \
               count[xp]-count[k] > rm-max_miss:
                x = (tmean[count[k]]-tmean[count[xm]])/(count[k]-count[xm])
                y = (tmean[count[xp]]-tmean[count[k]])/(count[xp]-count[k])
                xy = (tmean[count[xp]]-tmean[count[xm]])/(count[xp]-count[xm])

                sig = (tsquare[count[xp]]-tsquare[count[xm]])/(count[xp] -
                                                               count[xm])
                if(sig > xy*xy):
                    sig = math.sqrt(sig-xy*xy)
                    tsa[k] = ((count[k]-count[xm])*(x-xy)*(x-xy) +
                              (count[xp]-count[k])*(y-xy)*(y-xy))/sig
                else:
                    tsa[k] = 0.

                print xm, k, xp, tsa[k], sig, count[k], \
                    count[xp]-count[k], count[k]-count[xm]
    return


def zonaltrends(gslopes):

    s = gslopes.shape
    zslopes = numpy.zeros([s[0], s[3]])
    for k in range(s[0]):
        for ip in range(s[3]):
            zslopes[k, ip] = 0.
            b = 0
            for ipar in range(s[2]):
                mask = ~numpy.isnan(gslopes[k, :, ipar, ip])
                b += sum(mask)
                zslopes[k, ip] += sum(gslopes[k, mask, ipar, ip])
            if(b > 0):
                zslopes[k, ip] /= b
            else:
                zslopes[k, ip] = numpy.nan
    return zslopes


def belttrends(zslopes, belts):
    beltslopes = numpy.zeros([belts.shape[0], zslopes.shape[1]])
    for k in range(belts.shape[0]):
        for ip in range(zslopes.shape[1]):
            hilf = zslopes[belts[k, 0]:belts[k, 1], ip]
            mask = ~numpy.isnan(hilf)
            b = sum(mask)
            if(b > 0):
                beltslopes[k, ip] = sum(hilf[mask])/b
            else:
                beltslopes[k, ip] = numpy.nan
    return beltslopes




def rgb(r, g, b):
    return numpy.asarray([r, g, b], dtype=numpy.float)


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

