#!/usr/bin/env python
import traceback
import sys,glob
import os.path
import copy
import threading
from multiprocessing import Pool

from scipy.io import netcdf

#from ecmwfapi import ECMWFDataServer
import numpy
from gribapi import *
from datetime import date
import netCDF4
import time
from numba import *
from rasotools.utils import *
import matplotlib.pylab as plt
#from readCHUAN import read_ERACLIM_UAInv
import shutil
from functools import partial
#from pythransplit import esplit
#import fallocate

def netcdf4speedtest(fn):

    tt=time.time()
    ftn='/home/srvx7/leo/fastscratch/ei6/001001/001001_t'+'.nc'
    shutil.copy(ftn,fn)
    print(time.time()-tt)
    f=netCDF4.Dataset(ftn,'r')
    fo=netCDF4.Dataset(fn,'r+')
    var=f.variables['temperatures']
    for iv in range(20):
        fo.createVariable('var{}'.format(iv),var.dtype,var.dimensions,fill_value=numpy.nan)
        fo.variables['var{}'.format(iv)][:]=var[:]
        print(time.time()-tt)

        var=f.variables['temperatures']
        for i in list(fo.variables.keys()):
            if 'var10' not in i:
                continue
            for j in var.ncattrs():
                if j!='_FillValue' and j!='scale_factor' and j!='add_offset':
                    if i!='flags':
                        io=i
#			if str(i)=='temperatures' or i==refinput['parameter']:
#			    io=refinput['parameter']
#			    setattr(fo.variables[io],j,getattr(var,j))
                        try:
                            setattr(fo.variables[i],j,getattr(var,j))
#			    if io==refinput['parameter'] and j=='long_name':
#				setattr(fo.variables[io],j,refinput['parameter'])				
                            if j=='missing_value':
                                if fo.variables[io].dtype in (numpy.dtype('float32'),numpy.dtype('float64')):
                                    setattr(fo.variables[io],j,numpy.nan)
                                if fo.variables[io].dtype in (numpy.dtype('int32'),numpy.dtype('int64')):
                                    setattr(fo.variables[io],j,-999)
#			    if refinput['pk']!='t':
#				if io==refinput['parameter'] and j=='valid_range':
#				    setattr(fo.variables[io],j,[-100.,100.])
#				if io==refinput['parameter'] and j=='units':
#				    setattr(fo.variables[io],j,'m/s')
                        except KeyError:
                            pass

                        if i=='datum' and j=='units':
                            setattr(fo.variables[i],j,'days since 1900-01-01 00:00:00')
    f.close()
    print(time.time()-tt)
    print(os.path.getsize(fn),os.path.getsize(ftn))


    return

@njit
def sphdist_pre(x0,lat1,y0,lon1,z0,length):
    d=math.pi/180.
#    x0=numpy.cos(lat0*d)*numpy.cos(lon0*d)
    x1=numpy.cos(lat1*d)*numpy.cos(lon1*d)
#    y0=numpy.cos(lat0*d)*numpy.sin(lon0*d)
    y1=numpy.cos(lat1*d)*numpy.sin(lon1*d)
#    z0=numpy.sin(lat0*d)
    z1=numpy.sin(lat1*d)
    sp=numpy.empty(length)
    for i in range(length):
        sp[i]=x0[i]*x1+y0[i]*y1+z0[i]*z1
    spm=numpy.argmax(sp)
    if sp[spm]>1.0:
        sp[spm]=1.0
    ang=numpy.arccos(sp[spm])/d #numpy.sqrt(x0*x0+y0*y0+z0*z0)/numpy.sqrt(x1*x1+y1*y1+z1*z1))/d

    return (ang,spm)

@njit
def sphdist(lat0,lat1,lon0,lon1):
    d=math.pi/180.
    x0=numpy.cos(lat0*d)*numpy.cos(lon0*d)
    x1=numpy.cos(lat1*d)*numpy.cos(lon1*d)
    y0=numpy.cos(lat0*d)*numpy.sin(lon0*d)
    y1=numpy.cos(lat1*d)*numpy.sin(lon1*d)
    z0=numpy.sin(lat0*d)
    z1=numpy.sin(lat1*d)
    sp=x0*x1+y0*y1+z0*z1
    spm=numpy.argmax(sp)
    if sp[spm]>1.0:
        sp[spm]=1.0
    ang=numpy.arccos(sp[spm])/d #numpy.sqrt(x0*x0+y0*y0+z0*z0)/numpy.sqrt(x1*x1+y1*y1+z1*z1))/d

    return (ang,spm)

@njit
def sphdistall(lat0,lat1,lon0,lon1):
    d=math.pi/180.
    x0=numpy.cos(lat0*d)*numpy.cos(lon0*d)
    x1=numpy.cos(lat1*d)*numpy.cos(lon1*d)
    y0=numpy.cos(lat0*d)*numpy.sin(lon0*d)
    y1=numpy.cos(lat1*d)*numpy.sin(lon1*d)
    z0=numpy.sin(lat0*d)
    z1=numpy.sin(lat1*d)
    sp=x0*x1+y0*y1+z0*z1
    spm=numpy.argmax(sp)

    sp[sp>1.0]=1.0	
    ang=numpy.arccos(sp)/d

    return (ang,spm)

def read_wban(fn):
    f=open(fn,'r')
    rdata=f.read().split('\n')
    f.close()
    wbanids=[]
    wbanlons=[]
    wbanlats=[]
    wmoids=[]
    wmonames=[]
    for r in rdata[2:]:
        if len(r)>0:
            wbanids.append(r[10:15])
            lat=float(r[179:182])+float(r[183:185])/60.+float(r[186:188])/3600
            wbanlats.append(lat)
            lon=float(r[189:193])+float(r[194:196])/60.+float(r[197:200])/3600
            wbanlons.append(lon)
            wmoids.append(r[16:21])
            wmonames.append(r[130:160])

    return wbanids,wbanlats,wbanlons,wmoids,wmonames

@njit(cache=True)
def cscore1(tgrid,fv,offset,missval,corr,tidx,im,i,corr2):
    if (i-offset)%2!=0:
        return
    vs=fv.shape
    for j in range(tidx.shape[0]):
        if i//2>=tidx[j]:
            ik=j
    for j in range(vs[1]):
        for k in range(vs[2]):
            for l in range(vs[3]):
                if fv[i//2,j,k,l]!=missval:
                    if corr.size>1:
                        tgrid[i,vs[1]-j-1,k,l]=fv[i//2,j,k,l]-corr[im+ik,offset,vs[1]-j-1,k,l]-corr2[im+ik,offset,vs[1]-j-1,k,l]
                    else:
                        tgrid[i,vs[1]-j-1,k,l]=fv[i//2,j,k,l]-corr2[im+ik,offset,vs[1]-j-1,k,l]			
                else:
                    tgrid[i,vs[1]-j-1,k,l]=numpy.nan


@njit(cache=True)
def cscore2(tgrid,fv,offset,missval,im,i,corr2):
    if (i-offset)%2==0:
        vs=fv.shape
        for j in range(tidx.shape[0]):
            if i//2>=tidx[j]:
                ik=j
        for j in range(vs[1]):
            for k in range(vs[2]):
                for l in range(vs[3]):
                    if fv[i//2,j,k,l]!=missval:
                        tgrid[i,vs[1]-j-1,k,l]=fv[i//2,j,k,l]-corr2[im+ik,offset,vs[1]-j-1,k,l]
                    else:
                        tgrid[i,vs[1]-j-1,k,l]=numpy.nan

#@njit(cache=True)
def cscore3(tgrid,fv,offset,missval,i):
    if (i-offset)%2==0:
        vs=fv.shape
        for j in range(vs[1]):
            for k in range(vs[2]):
                tgrid[i,vs[1]-j-1,k,:]=fv[i//2,j,k,:]
                mask=fv[i//2,j,k,:]==missval
                tgrid[i,vs[1]-j-1,k,mask]=numpy.nan
#		for l in range(vs[3]):
                    #if fv[i/2,j,k,l]!=missval:
                        #tgrid[i,vs[1]-j-1,k,l]=fv[i/2,j,k,l]
                    #else:
                        #tgrid[i,vs[1]-j-1,k,l]=numpy.nan

#@njit(parallel=True)
def csplit(tgrid,fv,offset,missval,corr,tidx,im,corr2):

    vs=fv.shape
    if im==-1:
        for i in prange(offset,vs[0]*2):
            cscore2(tgrid,fv,offset,missval,i,corr2)
            #for j in range(vs[1]):
                #for k in range(vs[2]):
                    #for l in range(vs[3]):
                        #if fv[i/2,j,k,l]!=missval:
                            #tgrid[i,vs[1]-j-1,k,l]=fv[i/2,j,k,l]
                        #else:
                            #tgrid[i,vs[1]-j-1,k,l]=numpy.nan
    else:
        for i in prange(offset,vs[0]*2):
            cscore1(tgrid,fv,offset,missval,corr,tidx,im,i,corr2)
            #for j in range(tidx.shape[0]):
                #if i/2>=tidx[j]:
                    #ik=j
            #for j in range(vs[1]):
                #for k in range(vs[2]):
                    #for l in range(vs[3]):
                        #if fv[i/2,j,k,l]!=missval:
                            #tgrid[i,vs[1]-j-1,k,l]=fv[i/2,j,k,l]-corr[im+ik,offset,vs[1]-j-1,k,l]
                        #else:
                            #tgrid[i,vs[1]-j-1,k,l]=numpy.nan

    return

def dsplit(tgrid,fv,offset,missval,corr,tidx,im):

    vs=fv.shape
    if im==-1:
        for i in range(offset,vs[0]*2,2):
#	    for j in range(vs[1]):
#		for k in range(vs[2]):
#		    for l in range(vs[3]):
            h=fv[i/2,:,:,:].astype(numpy.float32)
#	    if i/2==0:
#		idx=numpy.where(h==missval)
#	    print numpy.sum(h==miss_val)
            h[h==missval]=numpy.nan
            tgrid[i,::-1,:,:]=h #fv[i/2,:,:,:]
    else:

        for i in range(offset,vs[0]*2,2):
            for j in range(tidx.shape[0]):
                if i/2>=tidx[j]:
                    ik=j
#	    for j in range(vs[1]):
#		for k in range(vs[2]):
#		    for l in range(vs[3]):
            h=fv[i/2,:,:,:].astype(numpy.float32)
#		if i/2==0:
#		    idx=numpy.where(h==missval)
#		print numpy.sum(h==miss_val)
            h[h==missval]=numpy.nan

#			if fv[i/2,j,k,l]!=missval:
            tgrid[i,::-1,:,:]=h-corr[im+ik,offset,::-1,:,:]
#			else:
#			    tgrid[i,::-1,:,:]=numpy.nan

    return

@njit(cache=True)
def latlon_to_station(ifield,ofield,ofieldref,hi,hw,idatum,odatum,hours,fcstep):


    fcshift=0
    if fcstep==12:
        fcshift=1

    iold=0
    i=0
    ot=0
    while ot<odatum.shape[0]:
        while odatum[ot]<idatum[0] and ot<odatum.shape[0]-1:
            ot+=1
        if odatum[ot]<idatum[0]:
            break   
        while idatum[i]<odatum[ot] and i<idatum.shape[0]-1:
            i+=1
        if i==idatum.shape[0]-1 and idatum[i]<odatum[ot]:
            break
        for ih in range(ofield.shape[0]):
            h=hours[ih,ot]
            it=i
            itn=i
            if h==-999:
                h=ih*12
            else:
                if h%3 !=0:
                    pass
                if h<0 or h>23:
                    continue

#	    ohv=hi[h,0]-fcshift
            #print otxv,otxn
            for il in range(ofield.shape[1]):
                if ofieldref[ih,il,ot]==ofieldref[ih,il,ot]:
                    ofield[ih,il,ot]=hw[h,0]*ifield[hi[h,0],il,it+hi[h,2]]+hw[h,1]*ifield[hi[h,1],il,it]	
        ot+=1

    return 

#@njit(cache=False,parallel=True)
def sgrid_to_station(ifield,ofield,ofieldref,xi,yi,lons,lats,weights,hi,hw,datum,hours,tx,md,stride,preverse,fcstep):
    for s in range(len(ofield)):
        pgrid_to_station(ifield,ofield[s],ofieldref[s],xi,yi,lons[s],lats[s],weights[s],hi,hw,datum[s],hours[s],tx,md,stride,preverse,fcstep)

    return 
#@njit(cache=False,parallel=True)
#def grid_to_station(ifield,ofield,ofieldref,xi,yi,xoo,yoo,weights,hi,hw,datum,hours,tx,md,stride,preverse,fcstep):
    #its=[]
    #for it in range(datum.shape[0]):
        #if datum[it]-1>=tx and datum[it]-1<tx+md+1:         #Assumption: datum has FORTRAN index, tx has C index thus -1 !!!
            #its.append(it)
    #for it in prange(its[0],its[-1]+1):
        #pgrid_to_station(ifield,ofield,ofieldref,xi,yi,xoo,yoo,weights,hi,hw,datum,hours,tx,md,stride,preverse,fcstep,it)
@njit(cache=True)
def pgrid_to_station(ifield,ofield,ofieldref,xi,yi,xoo,yoo,weights,hi,hw,datum,hours,tx,md,stride,preverse,fcstep):

    dxi=360./xi.shape[0]
    if yi.shape[0]%2==1:
        dyj=180./(yi.shape[0]-1)
        idyj=numpy.int(yi.shape[0]/2-numpy.floor(yoo/dyj))
    else:
        idyj=0
        while yi[idyj]>yoo:
            idyj+=1
            if idyj==yi.shape[0]:
                break

    xo=xoo
    yo=yoo

    weights_initialized=weights[0]!=0 or weights[1]!=0 or weights[2]!=0 or weights[3]!=0

    if idyj==yi.shape[0]:
        sp=True
        idyj-=1
    else:
        sp=False
    idyjm1=idyj-1
    dyj=yi[idyjm1]-yi[idyj]
    if xo<0.:
        xo+=360.
    idxi=numpy.int(numpy.floor(xo/dxi)+1)
    idxim1=idxi-1
    if idxi==xi.shape[0]:
        idxi=0
    if idxim1<0:
        idxim1+=xi.shape[0]
    if not weights_initialized:
        dyp=-(yi[idyj]-yo)
        dypm1=-(yo-yi[idyjm1])
        dxp=xi[idxi]-xo
        if dxp <-dxi:
            dxp+=360.
        elif dxp<0.:
            dxp=0.
        if dxp >=360:
            dxp-=360.
        dxpm1=xo-xi[idxim1]
        if dxpm1 <-dxi:
            dxpm1+=360
        elif dxpm1<0.:
            dxpm1=0.
        if dxpm1 >=360:
            dxpm1-=360.

        if sp:
            dyp=dyj
            dypm1=0.
        weights[0]=dyp/dyj*dxp/dxi
        weights[1]=dypm1/dyj*dxp/dxi
        weights[2]=dyp/dyj*dxpm1/dxi
        weights[3]=dypm1/dyj*dxpm1/dxi
        xo=weights[0]+weights[1]+weights[2]+weights[3]
        weights/=xo #numpy.sum(weights)
#        if sum(weights[:])<0.99999 or sum(weights[:])>1.00001:            
#            print xo,yo,dxi,dyj,weights[:],sum(weights[:])

    ia=int(ifield.shape[0]/md)

    fcshift=0
    if fcstep==12:
        fcshift=1
    for it in range(datum.shape[0]):
        if datum[it]-1>=tx and datum[it]-1<tx+md+fcshift:         #Assumption: datum has FORTRAN index, tx has C index thus -1 !!!
            for ih in range(ofield.shape[0]):
                h=hours[ih,it]
                if h==-999:
                    h=ih*12
                else:
                    if h==23:
                        h=0
                    if h<0 or h>23:
                        continue
                otxv=ia*(datum[it]-1-tx+hi[h,2])+hi[h,0]-fcshift
                if otxv<0 or otxv>=ifield.shape[0]:
                    continue
                if otxv<ifield.shape[0]-1:
                    otxn=otxv+1
                if preverse==1:
                    for il in range(ofield.shape[1]):
                        ilr=ofield.shape[1]-il-1
                        if ofieldref[ih,il,it]==ofieldref[ih,il,it]:
                            ofield[ih,il,it]=hw[h,0]*(weights[0]*ifield[otxv,ilr,idyjm1,idxim1]+weights[1]*ifield[otxv,ilr,idyj,idxim1]+ \
                                                      weights[2]*ifield[otxv,ilr,idyjm1,idxi]+weights[3]*ifield[otxv,ilr,idyj,idxi])+\
                                hw[h,1]*(weights[0]*ifield[otxn,ilr,idyjm1,idxim1]+weights[1]*ifield[otxn,ilr,idyj,idxim1]+ \
                                         weights[2]*ifield[otxn,ilr,idyjm1,idxi]+weights[3]*ifield[otxn,ilr,idyj,idxi])			        
                else:
                    for il in range(ofield.shape[1]):
                        if ofieldref[ih,il,it]==ofieldref[ih,il,it]:
                            ofield[ih,il,it]=hw[h,0]*(weights[0]*ifield[otxv,il,idyjm1,idxim1]+weights[1]*ifield[otxv,il,idyj,idxim1]+ \
                                                      weights[2]*ifield[otxv,il,idyjm1,idxi]+weights[3]*ifield[otxv,il,idyj,idxi])+\
                                hw[h,1]*(weights[0]*ifield[otxn,il,idyjm1,idxim1]+weights[1]*ifield[otxn,il,idyj,idxim1]+ \
                                      weights[2]*ifield[otxn,il,idyjm1,idxi]+weights[3]*ifield[otxn,il,idyj,idxi])	


    return 

@njit(cache=False)
def grid_to_gpsstation(ifield,ofield,ofieldref,xi,yi,xoo,yoo,weights,hi,hw,datum,refdatum,hours,tx,md,stride,preverse,fcstep):

    dxi=360./xi.shape[0]

#        if sum(weights[:])<0.99999 or sum(weights[:])>1.00001:            
#            print xo,yo,dxi,dyj,weights[:],sum(weights[:])

    ia=int(ifield.shape[0]/md)
    itref=0
    for it in range(datum.shape[0]):
#	print datum[it],tx,tx+md
        fcshift=0
        if fcstep==12:
            fcshift=1

        if datum[it]-1>=tx and datum[it]-1<tx+md+fcshift:  #Assumption: datum has FORTRAN index, tx has C index thus -1 !!!
            while(refdatum[itref]<datum[it] and itref<refdatum.shape[0]):
                itref+=1
            if refdatum[itref]!=datum[it]:
                continue

            for ih in range(ofield.shape[0]):
#                 otx=numpy.int(stride*2*(datum[it]-1-tx)+stride*ih)
                h=hours[ih,it]
                if h==-999:
                    h=ih*12
                else:
                    if h==23:
                        h=0
                    if h<0 or h>23:
                        continue
#		print it,h
                otxv=ia*(datum[it]-1-tx+hi[h,2])+hi[h,0]-fcshift
                if otxv<0 or otxv>=ifield.shape[0]:
                    continue
                if otxv<ifield.shape[0]-1:
                    otxn=otxv+1
                #print otxv,otxn
                if preverse==1:
                    for il in range(ofield.shape[1]):
                        ilr=ofield.shape[1]-il-1
                        xo=xoo[ih,il,it]
                        yo=yoo[ih,il,it]
                        if xo!=xo or yo!=yo:
                            continue
                        if yi.shape[0]%2==1:
                            dyj=180./(yi.shape[0]-1)
                            idyj=numpy.int(yi.shape[0]/2-numpy.floor(yo/dyj))
                        else:
                            idyj=0
                            while yi[idyj]>yo:
                                idyj+=1
                                if idyj==yi.shape[0]:
                                    break
                        if idyj==yi.shape[0]:
                            sp=True
                            idyj-=1
                        else:
                            sp=False
                        idyjm1=idyj-1
                        dyj=yi[idyjm1]-yi[idyj]
                        if xo<0.:
                            xo+=360.
                        idxi=numpy.int(numpy.floor(xo/dxi)+1)
                        idxim1=idxi-1
                        if idxi==xi.shape[0]:
                            idxi=0
                        if idxim1<0:
                            idxim1+=xi.shape[0]
#			if not weights_initialized:
                        dyp=-(yi[idyj]-yo)
                        dypm1=-(yo-yi[idyjm1])
                        dxp=xi[idxi]-xo
                        if dxp <-dxi:
                            dxp+=360.
                        elif dxp<0.:
                            dxp=0.
                        if dxp >=360:
                            dxp-=360.
                        dxpm1=xo-xi[idxim1]
                        if dxpm1 <-dxi:
                            dxpm1+=360
                        elif dxpm1<0.:
                            dxpm1=0.
                        if dxpm1 >=360:
                            dxpm1-=360.

                        if sp:
                            dyp=dyj
                            dypm1=0.
                        weights[0]=dyp/dyj*dxp/dxi
                        weights[1]=dypm1/dyj*dxp/dxi
                        weights[2]=dyp/dyj*dxpm1/dxi
                        weights[3]=dypm1/dyj*dxpm1/dxi
                        xo=weights[0]+weights[1]+weights[2]+weights[3]
                        weights/=xo #numpy.sum(weights)
                        if ofieldref[ih,il,itref]==ofieldref[ih,il,itref]:
                            ofield[ih,il,it]=hw[h,0]*(weights[0]*ifield[otxv,ilr,idyjm1,idxim1]+weights[1]*ifield[otxv,ilr,idyj,idxim1]+ \
                                                      weights[2]*ifield[otxv,ilr,idyjm1,idxi]+weights[3]*ifield[otxv,ilr,idyj,idxi])+\
                                hw[h,1]*(weights[0]*ifield[otxn,ilr,idyjm1,idxim1]+weights[1]*ifield[otxn,ilr,idyj,idxim1]+ \
                                         weights[2]*ifield[otxn,ilr,idyjm1,idxi]+weights[3]*ifield[otxn,ilr,idyj,idxi])			        
                else:
                    for il in range(ofield.shape[1]):
                        xo=xoo[ih,il,it]
                        yo=yoo[ih,il,it]
                        if xo!=xo or yo!=yo or yo>90. or yo<-90.:
                            continue
                        if yi.shape[0]%2==1:
                            dyj=180./(yi.shape[0]-1)
                            idyj=numpy.int(yi.shape[0]/2-numpy.floor(yo/dyj))
                        else:
                            idyj=0
                            while yi[idyj]>yo:
                                idyj+=1
                                if idyj==yi.shape[0]:
                                    break
                        if idyj==yi.shape[0]:
                            sp=True
                            idyj-=1
                        else:
                            sp=False
                        idyjm1=idyj-1
                        dyj=yi[idyjm1]-yi[idyj]
                        if xo<0.:
                            xo+=360.
                        idxi=numpy.int(numpy.floor(xo/dxi)+1)
                        idxim1=idxi-1
                        if idxi==xi.shape[0]:
                            idxi=0
                        if idxim1<0:
                            idxim1+=xi.shape[0]
#			if not weights_initialized:
                        dyp=-(yi[idyj]-yo)
                        dypm1=-(yo-yi[idyjm1])
                        dxp=xi[idxi]-xo
                        if dxp <-dxi:
                            dxp+=360.
                        elif dxp<0.:
                            dxp=0.
                        if dxp >=360:
                            dxp-=360.
                        dxpm1=xo-xi[idxim1]
                        if dxpm1 <-dxi:
                            dxpm1+=360
                        elif dxpm1<0.:
                            dxpm1=0.
                        if dxpm1 >=360:
                            dxpm1-=360.

                        if sp:
                            dyp=dyj
                            dypm1=0.
                        weights[0]=dyp/dyj*dxp/dxi
                        weights[1]=dypm1/dyj*dxp/dxi
                        weights[2]=dyp/dyj*dxpm1/dxi
                        weights[3]=dypm1/dyj*dxpm1/dxi
                        xo=weights[0]+weights[1]+weights[2]+weights[3]
                        weights/=xo #numpy.sum(weights)
                        if ofieldref[ih,il,itref]==ofieldref[ih,il,itref]:
                            ofield[ih,il,it]=hw[h,0]*(weights[0]*ifield[otxv,il,idyjm1,idxim1]+weights[1]*ifield[otxv,il,idyj,idxim1]+ \
                                                      weights[2]*ifield[otxv,il,idyjm1,idxi]+weights[3]*ifield[otxv,il,idyj,idxi])+\
                                hw[h,1]*(weights[0]*ifield[otxn,il,idyjm1,idxim1]+weights[1]*ifield[otxn,il,idyj,idxim1]+ \
                                      weights[2]*ifield[otxn,il,idyjm1,idxi]+weights[3]*ifield[otxn,il,idyj,idxi])	
                #if numpy.abs(ofield[ih,il,it])>1000.:


    return 

@njit(cache=False,parallel=True)
def quantilecheck(arr,sep=-1):
    if sep>0 and sep<arr.shape[2]:
        h1=arr[:,:,:sep]
        for ip in prange(h1.shape[1]):
            pquantilecheck(h1,ip)
        h2=arr[:,:,sep:]
        for ip in prange(h2.shape[1]):
            pquantilecheck(h2,ip)
        arr[:,:,:sep]=h1[:]
        arr[:,:,sep:]=h2[:]
    else:
        for ip in prange(arr.shape[1]):
            pquantilecheck(arr,ip)

@njit(cache=True)
def pquantilecheck(arr,ip):

    for ipar in range(arr.shape[0]):
#	for ip in range(arr.shape[1]):
        h=arr[ipar,ip,:]
        h=h[~numpy.isnan(h)]
        h=numpy.sort(h)
        l=len(h)
        if l>3:
            perc=numpy.empty(3)
            perc[1]=h[l//2]
            perc[0]=h[l//4]
            perc[2]=h[3*l//4]


            #perc=numpy.percentile(h[mask],[25.,50.,75.])
        ##perc=numpy.nanpercentile(arr[ipar,ip,:],[25.,50.,75.])
            pdmax=perc[1]+6*(perc[2]-perc[0])
            pdmin=perc[1]-6*(perc[2]-perc[0])
            for l in range(arr.shape[2]):
                if arr[ipar,ip,l]<pdmin or arr[ipar,ip,l]>pdmax:
                    arr[ipar,ip,l]=numpy.nan

                #pmax=h>pdmax
                #pmin=h<pdmin
                #arr[ipar,ip,pmin+pmax]=numpy.nan
            #except:
                #pass
#@njit
def quantilechecks(arr):

    for ipar in range(arr.shape[0]):
        for ip in range(arr.shape[1]):
            h=arr[ipar,ip,:]
            mask=~numpy.isnan(h)
            try:
                perc=numpy.percentile(h[mask],[25.,50.,75.])
            #perc=numpy.nanpercentile(arr[ipar,ip,:],[25.,50.,75.])
                pdmax=perc[1]+6*(perc[2]-perc[0])
                pdmin=perc[1]-6*(perc[2]-perc[0])
                pmax=h>pdmax
                pmin=h<pdmin
                arr[ipar,ip,pmin+pmax]=numpy.nan
            except:
                pass
#	    except TypeError:
#		pass
#	    n=numpy.sum(~numpy.isnan(arr[ipar,ip,:]))
#	    print ipar,ip,n,pdmax-pdmin,nmax*1.0/n,nmin*1.0/n
#            print ipar,ip,numpy.sum(pmax+pmin), 'values removed'

    return

def read_temperature(stats,st,fn,tidx,temperaturesdict,hours,tempname='temperatures'):

    indexmax=temperaturesdict[tempname].shape[2]

    if 'mdatum' not in list(stats.keys()):
        try:
            f=netCDF4.Dataset(fn,'r')
            f.set_auto_mask(False)
            s=getattr(f.variables['datum'],'units').split()[2].split('-')[0]
            shift=tidx[(int(s)-1900)*12]+1
            if len(fn.split('/')[-1])==11:
                shift -=1
            if f.variables['datum']==2:  
                stats['mdatum']=f.variables['datum'][0]+shift
                print(fn,f.variables['datum'][0].shape[0])
            else:
                stats['mdatum']=f.variables['datum']+shift
                print(fn,f.variables['datum'].shape[0])
            if numpy.abs(f.variables['lat'][:].flatten()[0]-stats['lat'])>1.0  or numpy.abs(f.variables['lon'][:].flatten()[0]-stats['lon'])>1.0:
                print('LAT LON INCONSISTENT:',st,fn,f.variables['lat'][:].flatten()[0],stats['lat'],f.variables['lon'][:].flatten()[0],stats['lon'])
                if st=='096645':
                    print(stats['lat'],stats['lon'],'is correct and used')
                else:
                    f.close()
                    return
            stats['mhours']=f.variables['hours'][:]
            for var in [tempname,'fg_dep','an_dep']:
                try:
                    miss_val=getattr(f.variables[var],'missing_value')
                    vals=f.variables[var][:]
                    vals[vals==miss_val]=numpy.nan
                    vals[vals==9.96920996839e+36]=numpy.nan
                except:
                    vals=f.variables[var][:]
                    pass
                if 'presat' in fn or '1761' in fn or '1759' in fn:
                    it=numpy.where(f.variables['datum'][0,:]>20819+365)[0]
                    if len(it)>0:
                        vals[:,:,it[0]:]=numpy.nan
                if numpy.nanmax(vals)>900 or numpy.nanmin(vals)<-900:
                    print(fn,numpy.sum(vals>900)+numpy.sum(vals<-900),' spurious values, set to NaN')
                    vals[vals>900]=numpy.nan
                    vals[vals<-900]=numpy.nan
                stats['m'+var]=vals
            if 'merged' not in fn and ':' not in fn:
## departures in merged files are defined as bg-obs, not obs-bg, thus the sign change
                if numpy.sum(~numpy.isnan(f.variables['bias'][:]))>0:
                    b=f.variables['bias'][:]
                else:
                    b=numpy.zeros(stats['m'+tempname].shape)
            else:
                b=numpy.zeros(stats['m'+tempname].shape)

            stats['m'+tempname]-=b
            stats['mfg_dep']=b-stats['mfg_dep']    
            stats['man_dep']=b-stats['man_dep']
            #if len(fn.split('/')[-1])!=11:
                #print fn.split('/')[-1],': Changing sign of departures'
                #stats['mfg_dep']=-stats['mfg_dep']
                #stats['man_dep']=-stats['man_dep']


            stats['msource']=numpy.zeros(indexmax,numpy.dtype('|S8'))
            try:
                stats['msource'][stats['mdatum']]=numpy.asarray(f.variables['source'][0],dtype='|S8')
            except:
                stats['msource'][stats['mdatum']]='BUFRDATA'

            stats['mstations']=numpy.zeros((indexmax),numpy.dtype('|S44')) # can hold 5strings,4 commas
            try:
                fstring=[str(f.variables['odbstatid'][0])]
                if fstring[0][0]=='0':
                    fstring[0]=fstring[0][1:]
                while fstring[-1][0]=='0':
                    fstring.append(fstring[-1][1:])

                if len(fstring)>1:   
                    stats['mstations'][stats['mdatum']]=','.join(fstring) # caution: id is not always same as file name		
                else:
                    stats['mstations'][stats['mdatum']]=fstring[0]

            except:
                fstring=fn.split('/')[-2]
                if ':' not in fstring and len(fstring)==6 and fstring[0]=='0':
                    fstring=fstring[1:]

                print('odbstatid not found, using '+fstring+' as a proxy')
                stats['mstations'][stats['mdatum']]=fstring # caution: id is not always same as file name

            f.close()
            hours[:,stats['mdatum']]=stats['mhours']
            temperaturesdict[tempname][:,:,stats['mdatum']]=stats['m'+tempname]
            temperaturesdict['fg_dep'][:,:,stats['mdatum']]=stats['mfg_dep']
            temperaturesdict['an_dep'][:,:,stats['mdatum']]=stats['man_dep']
        except:
            p=False
            pass
    else:

#   merge
        ntimes=list()
        try:
            f=netCDF4.Dataset(fn,'r')
            if 'datum' in list(f.variables.keys()):
                if f.variables['datum'].ndim==2:
                    print(fn,f.variables['datum'][0].shape[0],f.variables['lat'][:].flatten()[0],f.variables['lon'][:].flatten()[0])
                else:
                    print(fn,f.variables['datum'].shape[0],f.variables['lat'][:].flatten()[0],f.variables['lon'][:].flatten()[0])
                    
                f.set_auto_mask(False)
                s=getattr(f.variables['datum'],'units').split()[2].split('-')[0]
                shift=tidx[(int(s)-1900)*12]+1
                if len(fn.split('/')[-1])==11:
                    shift -=1
                if f.variables['datum'].ndim==2:
                    times=f.variables['datum'][0]+shift
                else:
                    times=f.variables['datum']+shift
                if st in ('027612','021965','043279','085469') or st>'104000':
                    if numpy.min(times)<0 or numpy.max(times)>=indexmax:
                        bd=numpy.where(numpy.logical_or(times<0,times>=indexmax))[0] 
                        print(fn,': Bad time index found: ',numpy.min(times),numpy.max(times))
                        try:
                            if times[bd+1]-times[bd-1]==2:
                                times[bd]=(times[bd+1]+times[bd-1])/2
                                print(fn,': Bad time index fixed')
                            elif bd==0:
                                times[bd]=times[bd+1]-1
                                print(fn,': Bad time index fixed')
                            else:
                                return
                        except:
                            bd=numpy.where(numpy.logical_and(times>=0,times<indexmax))[0] 
                            times=times[bd]
                            return


                #hours[:,stats['mdatum']]=stats['mhours'][:]
                #hh=copy.copy(f.variables['hours'][:])
                #for ipar in range(2):
                    #mask=numpy.logical_and(hh[ipar,:]>-1,hh[ipar,:]<24)
                    #if sum(mask)>0:
                        #nhours[ipar,times[mask]]=hh[ipar,mask]

                #idx=numpy.where(nhours[1,:]==4)[0]
                #if idx.shape[0]>0:
                    #print' wrong hour'
                t=time.time()
                #mdathilf=numpy.asarray(numpy.where(ntimes>0)[0],dtype=numpy.int32)

                overwrite=numpy.zeros(times.shape[0],dtype=numpy.bool)
                #if '1759' not in fn:
                    #try:
                        #overwrite[:]=numpy.logical_and(stats['msource'][times]=='NCARUA20',f.variables['source'][0]!='NCARUA20')
                    #except KeyError:
                        #pass
                bad=False
                for var in [tempname,'fg_dep','an_dep']: #temperaturesdict.keys():
                    try:
                        mergevar(stats,fn,f.variables,times,hours,var,temperaturesdict[var],overwrite,tempname=tempname,thresh=0.1)
                    except (KeyError,ValueError) as merr:
                        if 'm'+var not in list(stats.keys()):
                            print('Variable '+var+'not found in file '+fn+', adding empty variable')
                            stats['m'+var]=numpy.empty((2,16,stats['mdatum'].shape[0]),dtype=numpy.float32)
                            stats['m'+var].fill(numpy.nan)
                        if 'm'+var=='m'+tempname:
                            bad=True
                        print(merr)

                tmerge=time.time()-t
                if tmerge>0.1:
                    print('merge',tmerge)
                if(bad):
                    f.close()
                    return
#                stats['mdatum']=



#mfg_dep not always available -> no check		
                #tt=time.time()
                #quantilecheck(stats['mfg_dep'])
                #delmask=numpy.isnan(stats['mfg_dep'])
                #stats['m'+tempname][delmask]=numpy.nan
                #try:
                    #stats['man_dep'][delmask]=numpy.nan
                #except:
                    #stats['man_dep']=numpy.empty(stats['m'+tempname].shape,dtype=numpy.float32)
                    #stats['man_dep'].fill(numpy.nan)

                #print 'qc:',time.time()-tt

            f.close()
        except (RuntimeError, IOError,AttributeError) as nerr:
            print(fn,nerr)
            p=False
            pass

def merge_gps(fn,stats):
    if 'wet' in fn:
        syll='gpswet'
        vn='wet'
    else:
        syll='gps'
        vn=''
    try:
        with netCDF4.Dataset(fn,'r') as f:

            for d in ['gpstemperatures','gpslat','gpslon']:
                stats['gps'+vn+d[3:]]=numpy.empty(f.variables[d].shape,dtype=numpy.float32)
                stats['gps'+vn+d[3:]][:]=f.variables[d][:]
            stats[syll+'hours']=numpy.empty(f.variables['gpshours'].shape,dtype=numpy.int)
            stats[syll+'hours'][:]=f.variables['gpshours'][:]
            stats[syll+'datum']=f.variables['datum'][0,:]
            for ds in ['erai']: #,'era5','jra55']:
                stats[ds+'_'+syll+'temperatures']=numpy.empty(f.variables['gpstemperatures'].shape,dtype=numpy.float32)
                stats[ds+'_'+syll+'temperatures'][:]=numpy.nan
                stats[ds+'_fg'+syll+'dep']=numpy.empty(stats['mtemperatures'].shape,dtype=numpy.float32)
                stats[ds+'_fg'+syll+'dep'][:]=numpy.nan
    except:
        for ds in ['erai']: #,'era5','jra55']:
            stats[ds+'_fg'+syll+'dep']=numpy.empty(stats['mtemperatures'].shape,dtype=numpy.float32)
            stats[ds+'_fg'+syll+'dep'][:]=numpy.nan

    return

def mergevar(stats,fn,variables,times,hours,var,temperatures,overwrite,tempname='temperatures',thresh=0.):
# times contains dates (as days since 19000101) of data to be merged in
# stats['mdatum'] contains dates (as days since 19000101) of existing data
# mdathilf contains dates (as days since 19000101) of both existing and merged data
#	    t1=time.time()
    try:
        miss_val=getattr(variables[var],'missing_value')
    except:
        print(fn,' NO missing value specified')
        pass
#            temperatures[:,:,stats['mdatum']]=stats['m'+var]
#	    t=time.time()
    vals=variables[var][:]
    vhours=variables['hours'][:]
#            vals[vals==miss_val]=numpy.nan
    if 'presat' in fn or '1761' in fn or '1759' in fn:
        it=numpy.where(variables['datum'][0,:]>20819+365)[0]
        if len(it)>0:
            vals[:,:,it[0]:]=numpy.nan
    vals[numpy.logical_or(vals>400.,vals<-100.)]=numpy.nan
#	    print'mergeread',time.time()-t
#            if 'merged' not in fn and ':' not in fn and 'ai_bfr' not in fn  and 'igra' not in fn:
    if numpy.sum(~numpy.isnan(variables['bias'][:,11,:]))>0:
        b=variables['bias'][:]
        if var==tempname:
            vals-=b
        elif var=='fg_dep':
            vals=b-vals
        elif var=='an_dep':
            vals=b-vals
    #if len(fn.split('/')[-1])!=11 and var!=tempname:
        #print fn.split('/')[-1],': Changing sign of departures'
##		vals=-vals
        #vals=-vals
# departures in merged files are defined as bg-obs, not obs-bg, thus the sign change

#	    print'mergeread1',time.time()-t
#            hilf=temperatures[:,:,times]
    if '3188' not in fn:
        pass
    if numpy.sum(~numpy.isnan(temperatures[:,:,times]))==0:  # new data - easy merge
        temperatures[:,:,times]=vals
        if var==tempname:
            try:
                stats['msource'][times]=numpy.asarray(variables['source'][0],dtype='|S8')
            except:
                stats['msource'][times]='BUFRDATA'
#		overwrite=vals.shape[2]
        inconsistent=0
    else:
#		overwrite=numpy.zeros(times.shape[0],dtype=numpy.bool)
#		if '1759' not in fn:
#		    overwrite[:]=numpy.logical_and(stats['msource'][times]=='NCARUA20',variables['source'][0]!='NCARUA20')
        ov=numpy.copy(overwrite)
#		print 'before', numpy.sum(~numpy.isnan(temperatures[0,7,:]))
        inconsistent=reccheck(temperatures,hours,times,vals,vhours,ov,thresh)
#		print 'after', numpy.sum(~numpy.isnan(temperatures[0,7,:]))
        #reccheck(temperatures,times,vals,thresh)
        if numpy.sum(hours>23):
            print(fn, numpy.sum(hours>23),'invalid values in hours')
        if var==tempname:
            try:
                stats['msource'][times[ov]]=variables['source'][0]
            except:
                stats['msource'][times[ov]]='BUFRDATA'
        idx=numpy.where(inconsistent>thresh)
        if numpy.sum(inconsistent[idx])>1:
            print(var,len(idx[0]),' records inconsistent:',fn)#.split('/')[-2], times[idx],inconsistent[idx]
#	    print 'mergevar: ',fn,numpy.sum(overwrite)

##	    mask=numpy.isnan(hilf)
##	    hilf[mask]=vals[mask]
    #print numpy.sum(mask),numpy.sum(numpy.isnan(hilf))
##            temperatures[:,:,times]=hilf
    if var==tempname:
#		stats['m'+var]=temperatures[:,:,mdathilf]

        tx=time.time()
        try:
            fstring=[str(variables['odbstatid'][0])]
            if fstring[0][0]=='0':
                fstring[0]=fstring[0][1:]
            while fstring[-1][0]=='0':
                fstring.append(fstring[-1][1:])


        except:
            fstring=fn.split('/')[-2]
            if ':' not in fstring and len(fstring)==6 and fstring[0]=='0':
                fstring=fstring[1:]
            fstring=[fstring]

        save="x"
        for i in range(len(times)):
            if stats['mstations'][times[i]]!=save:
                save=stats['mstations'][times[i]]
            else:
                stats['mstations'][times[i]]=sold
                continue
            if stats['mstations'][times[i]]=='':
                if len(fstring)<2:
                    stats['mstations'][times[i]]=fstring[0]
                else:
                    stats['mstations'][times[i]]=','.join(fstring)
            else: 
                slist=stats['mstations'][times[i]].decode('utf-8').split(',')
                for s in fstring:
                    if s not in slist:
                        slist.append(s)

                if len(slist)>1:   
                    stats['mstations'][times[i]]=','.join(slist) # caution: id is not always same as file name		
                else:
                    stats['mstations'][times[i]]=slist[0]
            sold=stats['mstations'][times[i]]

        print('mergeread2',time.time()-tx)
#	    print time.time()-t1
    return
@njit()
def reccheck(temperatures,hours,datum,vals,vhours,overwrite,thresh):
    #overwrite=numpy.zeros(datum.shape[0],dtype=boolean)
    inconsistent=numpy.zeros(datum.shape[0])
    #overwrite=numpy.zeros(datum.shape[0],dtype=numpy.bool)
    #inconsistent=numpy.zeros(datum.shape[0],dtype=numpy.bool)
    for it in range(datum.shape[0]):
        sq=0.
        dc=0.001
        for ipar in range(temperatures.shape[0]):  
            vc=0
            tc=0
            vg=0
            if vhours[ipar,it]>=0 and vhours[ipar,it]<24:
                for ip in range(temperatures.shape[1]):

                    diff=temperatures[ipar,ip,datum[it]]-vals[ipar,ip,it]
                    if diff==diff:
                        sq+=diff*diff
                        dc+=1
                        vg+=1
                        tc+=1
                    else:
                        if vals[ipar,ip,it]==vals[ipar,ip,it]:
                            temperatures[ipar,ip,datum[it]]=vals[ipar,ip,it]
                            vc+=1
                            vg+=1
                        elif temperatures[ipar,ip,datum[it]]==temperatures[ipar,ip,datum[it]]:
                            tc+=1
            else:
                if vhours[ipar,it]>23:
                    print('invalid')

            if tc<3 and vg>2:
                temperatures[ipar,:,datum[it]]=vals[ipar,:,it]
                hours[ipar,datum[it]]=vhours[ipar,it]
        #	    days[:,:,datum[it]]=True
                overwrite[it]=True
            else:
        #	    print it
                pass
        inconsistent[it]=numpy.sqrt(sq/dc)
        #if inconsistent[it]>0.1:
            #print it,inconsistent[it]

    return inconsistent





def add_latlonfeedback(ncpath,varname,stats,st,fcstep=0):


    hi=[]
    hw=[]
    dmd=4
    for i in range(0,dmd):
        for j in range(24/dmd):
            hi.append((i,(i+1)%dmd,-((i+1)/dmd)))
            hw.append((1.0-j*dmd/24.,j*dmd/24.))
    hw=numpy.asarray(hw)
    hi=numpy.asarray(hi,dtype='int')

    tt=time.time()
    if type(stats[st]["lon"]) is list:
        lon=stats[st]["lon"][0]
        lat=stats[st]["lat"][0]
    else:
        lon=stats[st]["lon"]
        lat=stats[st]["lat"]
    if lat!=lat or lon!=lon:
        return
    try:
        fn=ncpath+'/ref_{}_{}_t.nc'.format(int(lat*100),int(lon*100))
        f=netCDF4.Dataset(fn)
    except:
        print(fn,' trying..')
        lf=False
        for l in range(int(lat*100-2),int(lat*100+3)):
            for k in range(int(lon*100-2),int(lon*100+3)):
                try:
                    fnx=ncpath+'/ref_{}_{}_t.nc'.format(l,k)
                    f=netCDF4.Dataset(fnx)
                    print('found!',fnx)
                    lf=True
                    break
                except:
                    pass
            if lf:
                break
        if not lf:
            print(fn,' not found!')
            return
    if varname in list(f.variables.keys()):
        latlon_to_station(f.variables[varname][:],stats[st][varname],stats[st]['mtemperatures'],
                          hi,hw,f.variables['datum'][0,:],stats[st]['mdatum'],stats[st]['mhours'],fcstep)

            #plt.figure(figsize=(12,5))
            #plt.subplot(1,2,1)
            #plt.xlabel('days since 19000101')
            #plt.ylabel('Temperature [K]')
            #plt.legend()
            #plt.subplot(1,2,2)
            #plt.plot(smd,stats[st]['erai_fg'+syll+'dep'][0,5,:],'x',label=''+syll+'-erai_fg'+syll+'')
            #plt.plot(smd,stats[st]['erai_fgtemperatures'][0,5,:]-stats[st]['mtemperatures'][0,5,:],label='obs-erai_fg')
            #plt.xlabel('days since 19000101')
            #plt.ylabel('Temperature [K]')
            #plt.legend()


    print(time.time()-tt)
    return

def refdiff(ce20cpath,filepattern,start,stop,tidx,fieldsperday=2,fcstep=0,mrange=list(range(1,13)),tgrid=None,tempname='temperatures',iens=0):

    t=time.time()
    try:
        if 'av' in filepattern:
            mdiff=numpy.load(ce20cpath+tempname+'_{}-{}_av.npy'.format(start,stop))
        else:
            mdiff=numpy.load(ce20cpath+tempname+'_{}-{}_{}.npy'.format(start,stop,iens))
    except:
        il=-1
        ini=False
        for iy in range(start,stop):
            for im in mrange:
                t2=time.time()
                il+=1
                l=(iy-1900)*12+im-1
                tx=tidx[l]
                try:
                    md=tidx[l+1]-tidx[l]
                except IndexError:
                    md=31 # December assumed

                fn=ce20cpath+filepattern.format(iy,im)
                print(fn)
                if not os.path.isfile(fn):
                    continue
                f=open(fn,'r')

                gid=0
                d=0
                h=0
                ip=0
                while True:
                    gid = grib_new_from_file(f)
                    if gid is None:
                        break
                    hh=grib_get(gid,'time')
                    if hh%1200!=0:
                        grib_release(gid)
                        continue
                    if ip==0 and h==0 and d==0:
                        inc=-1+2*grib_get(gid,'jScansPositively')
                        ss=grib_get(gid,'step')
                        Ni=grib_get(gid,'Ni')
                        Nj=grib_get(gid,'Nj')
                        yi=grib_get(gid,'latitudeOfFirstGridPointInDegrees')+inc*numpy.arange(Nj)*grib_get(gid,'jDirectionIncrementInDegrees')
                        xi=grib_get(gid,'longitudeOfFirstGridPointInDegrees')+numpy.arange(Ni)*grib_get(gid,'iDirectionIncrementInDegrees')

                    if tgrid is None:
                        tgrid=numpy.empty((31,2,16,Nj,Ni),dtype=numpy.float32)
                        mdiff=numpy.empty((stop-start,12,2,16,Nj,Ni),dtype=numpy.float32)
                        ini=True

                    if ss==12:
                        if hh==1200:
                            tgrid[d,h-1,ip,:,:]=numpy.reshape(grib_get_values(gid),(Nj,Ni))
                        else:
                            tgrid[d,h+1,ip,:,:]=numpy.reshape(grib_get_values(gid),(Nj,Ni))
                    else:
                        tgrid[d,h,ip,:,:]=numpy.reshape(grib_get_values(gid),(Nj,Ni))

                    grib_release(gid)
                    ip+=1
                    if ip==16:
                        h+=1
                        ip=0
                        if h==2:
                            d+=1
                            h=0

                f.close()

                if d/md not in [1]:
                    print('grib file incomplete')
                    exit()
                mdiff[iy-start,im-1,:,:,:,:]=numpy.mean(tgrid,axis=0)



                print('after grib:',time.time()-t2)
        if 'av' in filepattern:
            numpy.save(ce20cpath+tempname+'_{}-{}_av'.format(start,stop),mdiff)
        else:
            numpy.save(ce20cpath+tempname+'_{}-{}_{}'.format(start,stop,iens),mdiff)
    print('refdiff',time.time()-t)
    return mdiff

def make_ens(ce20cpath,filepattern,start,stop,tempname='temperatures',iens=0):

    t=time.time()
    il=-1
    ini=False
    for iy in range(start,stop):
        for im in range(1,13):
            t2=time.time()

            fno=ce20cpath+(filepattern+'.av.grb').format(iy,im)
            fo=open(fno,'w')
            gids=[]
            av=[]
            for ie in range(10):

                fn=ce20cpath+filepattern.format(iy,im)+'.'+'{0}'.format(ie)+'.grb'
                print(fn)
                if not os.path.isfile(fn):
                    break
                f=open(fn,'r')

                k=0
                while True:
                    gid = grib_new_from_file(f)
                    if gid is None:
                        break
                    if ie==0:
                        gids.append(grib_clone(gid))
                        av.append(grib_get_values(gid))
                    elif ie==9:
                        av[k]+=grib_get_values(gid)
                        grib_set_values(gids[k],av[k]/10)
                        grib_write(gids[k],fo)
                        grib_release(gids[k])
                    else:
                        av[k]+=grib_get_values(gid)
                    k+=1

                    grib_release(gid)
                f.close()
            fo.close()
            print(time.time()-t2)
    print(time.time()-t)

    return

@njit(nogil=True)
def add_togrid(tgrid,x,d,ip,Nj,Ni):
    l=0
    for j in range(Nj):
        for i in range(Ni):
            tgrid[d,ip,j,i]=x[l]
            l+=1

def add_feedback_mp(ce20cpath,filepattern,varname,start,stop,tidx,stats,fieldsperday=4,fcstep=0,
                    mrange=list(range(1,13)),tgrid=None,sylls=[],tempname='temperatures',
                 corr=None,cstart=None,cstop=None,varname2=None):

#    ce20cpath='/nas/srvx8/leo/CERA20C/'
    t=time.time()

    global mystats
    mystats=stats

    func = partial(grread,ce20cpath,filepattern,varname,start,stop,tidx,varname2,fieldsperday,fcstep,sylls,tempname,
                   corr,cstart,cstop)
    p=Pool(12)
    out=p.map(func,mrange)
#    p=Pool(12)
#    out=list(map(func,mrange))

    tt=time.time()
    i=0
    for s in out[0][0]:
        for m in mrange:
            stats[s][varname][out[m-1][1][i]]=out[m-1][2][i]
            if varname2 is not None:
                stats[s][varname2][out[m-1][3][i]]=out[m-1][4][i]
        i+=1

    p.close()
    p.terminate()
    p.join()

    print('remap',time.time()-tt)    
    gpsname=[]
    for syll in sylls:
        gpsname+=['erai_fg'+syll+'dep'] # return departure string since this departure is already finished

    if len(gpsname)>0:
        return gpsname
    else:
        if varname2 is not None:
            return [varname,varname2]
        else:
            return [varname]


def grread(ce20cpath,filepattern,varname,start,stop,tidx,varname2,fieldsperday,fcstep,sylls,tempname,corr,cstart,cstop,m):   

# expects dictionary stats

    slist=[]
    lats=[]
    lons=[]
    svar=[]
    svar2=[]
    stemp=[]
    sweights=[]
    sdatum=[]
    shours=[]
    for st in list(mystats.keys()):
        if st in ['odbfile', 'header']:
            continue
        if 'm'+tempname in list(mystats[st].keys()):
            if type(mystats[st]["lon"]) is list:
                lon=mystats[st]["lon"][0]
                lat=mystats[st]["lat"][0]
            else:
                lon=mystats[st]["lon"]
                lat=mystats[st]["lat"]
            if lat!=lat or lon!=lon:
                continue
            svar.append(numpy.empty(mystats[st]['m'+tempname].shape,dtype=numpy.float32))
            svar[-1].fill(numpy.nan)
            if varname2 is not None :
                svar2.append(numpy.empty(mystats[st]['m'+tempname].shape,dtype=numpy.float32))
                svar2[-1].fill(numpy.nan)

        else:
            continue

        stemp.append(mystats[st]['m'+tempname])
        mystats[st]['weights'][:]=0.
        sweights.append(mystats[st]['weights'])
        sdatum.append(mystats[st]['mdatum'])
        shours.append(mystats[st]['mhours'])
        slist.append(st)
        lats.append(lat)
        lons.append(lon)	


    print('grread m',m)
    t=time.time()
    tgrid=None
    il=-1
    ini=False
    for iy in range(start,stop):
        for im in [m]:
            t2=time.time()
            il+=1
            l=(iy-1900)*12+im-1
            tx=tidx[l]
            try:
                md=tidx[l+1]-tidx[l]
            except IndexError:
                md=31 # December assumed
            fn=ce20cpath+filepattern.format(iy,im)
            print(fn)
            if not os.path.isfile(fn):
                continue
            f=open(fn,'r')

            gid=0
            d=0
            h=0
            ip=0
            while True:
                gid = grib_new_from_file(f)
                if gid is None:
                    break
                #hh=grib_get(gid,'time')
                #if hh%600!=0:
                    #grib_release(gid)
                    #continue
                if ip==0 and h==0 and d==0:
                    inc=-1+2*grib_get(gid,'jScansPositively')
                    Ni=grib_get(gid,'Ni')
                    Nj=grib_get(gid,'Nj')
                    yi=grib_get(gid,'latitudeOfFirstGridPointInDegrees')+inc*numpy.arange(Nj)*grib_get(gid,'jDirectionIncrementInDegrees')
                    xi=grib_get(gid,'longitudeOfFirstGridPointInDegrees')+numpy.arange(Ni)*grib_get(gid,'iDirectionIncrementInDegrees')

                    if tgrid is None:
                        tgrid=numpy.empty((31*fieldsperday,16,Nj,Ni),dtype=numpy.float32)
                    if varname2 is not None:
                        tgrid2=numpy.empty((31*fieldsperday,16,Nj,Ni),dtype=numpy.float32)
                    ini=True
                x=grib_get_values(gid)
                add_togrid(tgrid,x,d,ip,Nj,Ni)

                grib_release(gid)
                ip+=1
                if ip==16:
                    d+=1
                    ip=0
                    #if h==2:
                        #d+=1
                        #h=0
            f.close()

            if d/md not in [2,4]:
                print('grib file incomplete, skipped')
                continue


            try:
                hi=[]
                hw=[]
                dmd=d//md
                for i in range(0,dmd):
                    for j in range(24//dmd):
                        hi.append((i,(i+1)%dmd,-((i+1)//dmd)))
                        hw.append((1.0-j*dmd/24.,j*dmd/24.))
                hw=numpy.asarray(hw)
                hi=numpy.asarray(hi,dtype='int')

                print('after grib:',time.time()-t2)
                if len(svar2)>0:
                    sgrid_to_station(tgrid,tuple(svar2),tuple(stemp),xi,yi,
                                     tuple(lons),
                                                tuple(lats),
                                                tuple(sweights),hi,hw,
                                                tuple(sdatum),
                                                tuple(shours),
                                                tx,
                                                md,numpy.int(1),numpy.int(0),numpy.int(fcstep))

                if corr is not None and iy>=cstart and iy<cstop:
                    fpd=fieldsperday
                    tgrid[0:fpd*md:fpd,:,:,:]-=corr[im-1,0,:,:,:] # assumed averaged over years
                    if fpd==2:
                        tgrid[1:fpd*md:fpd,:,:,:]-=corr[im-1,1,:,:,:] # assumed averaged over years
                    else:
                        tgrid[2:fpd*md:fpd,:,:,:]-=corr[im-1,1,:,:,:] # assumed averaged over years
                        c=0.5*(corr[im-1,0,:,:,:]+corr[im-1,1,:,:,:])
                        tgrid[1:fpd*md:fpd,:,:,:]-=c # assumed averaged over years
                        tgrid[3:fpd*md:fpd,:,:,:]-=c # assumed averaged over years
                    print('after correction:',time.time()-t2)

                sgrid_to_station(tgrid,tuple(svar),tuple(stemp),xi,yi,
                                 tuple(lons),
                                    tuple(lats),
                                    tuple(sweights),hi,hw,
                                    tuple(sdatum),
                                    tuple(shours),
                                    tx,
                                    md,numpy.int(1),numpy.int(0),numpy.int(fcstep))
                for syll in sylls:
                    for st in list(mystats.keys()):
                        if st in ['odbfile', 'header'] or 'erai_'+syll+'temperatures' not in list(mystats[st].keys()):
                            continue

                        grid_to_gpsstation(tgrid,mystats[st]['erai_'+syll+'temperatures'],mystats[st]['m'+tempname],xi,yi,
                                           mystats[st][syll+'lon'],
                                    mystats[st][syll+"lat"],
                                    mystats[st]["weights"],hi,hw,
                                    mystats[st][syll+'datum'],
                                    mystats[st]['mdatum'],
                                    mystats[st][syll+'hours'],
                                    tx,
                                    md,numpy.int(1),numpy.int(0),numpy.int(fcstep))

                        #print 'gs2',time.time()-tt
                        mi=numpy.where(numpy.logical_and(mystats[st]['mdatum']>=tx,mystats[st]['mdatum']<tx+md ))[0]
                        gi=numpy.where(numpy.logical_and(mystats[st][syll+'datum']>=tx,mystats[st][syll+'datum']<tx+md ))[0]

                        #print len(mi),len(gi),time.time()-tt
                        l=0
                        lg=[]
                        sgd=mystats[st][syll+'datum']
                        smd=mystats[st]['mdatum']
                        if len(gi)>0:
                            for i in range(mi.shape[0]):
                                while l<gi.shape[0]-1 and sgd[gi[l]]<smd[mi[i]]:
                                    l+=1
                                if sgd[gi[l]]==smd[mi[i]] and numpy.sum(~numpy.isnan(mystats[st][syll+'temperatures'][:,:,gi[l]]))>0:
                                    mystats[st]['erai_fg'+syll+'dep'][:,:,mi[i]]=mystats[st]['erai_'+syll+tempname][:,:,gi[l]]-mystats[st][syll+'temperatures'][:,:,gi[l]]
                                    lg.append(i)


                print(iy*100+im)
                print(time.time()-t,time.time()-t2)
            except IOError:
                print(fn+' not found')
                pass
    idx=[]
    vals=[]
    for s in svar:	   
        idx.append(numpy.where(~numpy.isnan(s)))
        vals.append(s[idx[-1]])
    if len(svar2)>0:
        idx2=[]
        vals2=[]
        for s in svar2:	   
            idx2.append(numpy.where(~numpy.isnan(s)))
            vals2.append(s[idx2[-1]])
        return slist,idx,vals,idx2,vals2
    else:
        return slist,idx,vals

def add_feedback(ce20cpath,filepattern,varname,start,stop,tidx,stats,fieldsperday=4,fcstep=0,
                 mrange=list(range(1,13)),tgrid=None,sylls=[],tempname='temperatures',
                 corr=None,cstart=None,cstop=None,varname2=None):

#    ce20cpath='/nas/srvx8/leo/CERA20C/'
    t=time.time()
    slist=[]
    lats=[]
    lons=[]
    svar=[]
    svar2=[]
    stemp=[]
    sweights=[]
    sdatum=[]
    shours=[]
    for st in list(stats.keys()):
        if st in ['odbfile', 'header']:
            continue
        if 'm'+tempname in list(stats[st].keys()):
            if varname not in list(stats[st].keys()):
                stats[st][varname]=numpy.empty(stats[st]['m'+tempname].shape,dtype=numpy.float32)
                stats[st][varname].fill(numpy.nan)
            if varname2 is not None :
                if varname2 not in list(stats[st].keys()):
                    stats[st][varname2]=numpy.empty(stats[st]['m'+tempname].shape,dtype=numpy.float32)
                    stats[st][varname2].fill(numpy.nan)

            if type(stats[st]["lon"]) is list:
                lon=stats[st]["lon"][0]
                lat=stats[st]["lat"][0]
            else:
                lon=stats[st]["lon"]
                lat=stats[st]["lat"]
            if lat!=lat or lon!=lon:
                continue
        else:
            continue

        stemp.append(stats[st]['m'+tempname])
        svar.append(stats[st][varname])
        if varname2 is not None:
            svar2.append(stats[st][varname2])
        stats[st]['weights'][:]=0.
        sweights.append(stats[st]['weights'])
        sdatum.append(stats[st]['mdatum'])
        shours.append(stats[st]['mhours'])
        slist.append(st)
        lats.append(lat)
        lons.append(lon)
    il=-1
    ini=False
    for iy in range(start,stop):
        if corr is not None:
            if iy<cstart or iy>=cstop:
                continue
        for im in mrange:
            t2=time.time()
            il+=1
            l=(iy-1900)*12+im-1
            tx=tidx[l]
            try:
                md=tidx[l+1]-tidx[l]
            except IndexError:
                md=31 # December assumed
#            fn=ce20cpath+'CERA20C{0}{1:02}.grb'.format(iy,im)
            fn=ce20cpath+filepattern.format(iy,im)
            print(fn)
            if not os.path.isfile(fn):
                continue
            f=open(fn,'r')

            gid=0
            d=0
            h=0
            ip=0
            while True:
#		gid = grib_new_from_index(iid)
                gid = grib_new_from_file(f)
                if gid is None:
                    break
                #hh=grib_get(gid,'time')
                #if hh%600!=0:
                    #grib_release(gid)
                    #continue
                if ip==0 and h==0 and d==0:
                    inc=-1+2*grib_get(gid,'jScansPositively')
                    Ni=grib_get(gid,'Ni')
                    Nj=grib_get(gid,'Nj')
                    yi=grib_get(gid,'latitudeOfFirstGridPointInDegrees')+inc*numpy.arange(Nj)*grib_get(gid,'jDirectionIncrementInDegrees')
                    xi=grib_get(gid,'longitudeOfFirstGridPointInDegrees')+numpy.arange(Ni)*grib_get(gid,'iDirectionIncrementInDegrees')
#		pid=grib_get(gid,'paramId')
#		if pid==130:
                if tgrid is None:
                    tgrid=numpy.empty((31*fieldsperday,16,Nj,Ni),dtype=numpy.float32)
                    if varname2 is not None:
                        tgrid2=numpy.empty((31*fieldsperday,16,Nj,Ni),dtype=numpy.float32)
                    ini=True
#		tt=time.time()
                x=grib_get_values(gid)
                add_togrid(tgrid,x,d,ip,Nj,Ni)

#                tgrid[d,ip,:,:]=numpy.reshape(grib_get_values(gid),(Nj,Ni))
#		print time.time()-t
                #print d,grib_get(gid,'date'),hh
                grib_release(gid)
                ip+=1
                if ip==16:
                    d+=1
                    ip=0
#		    if h==2:
#			d+=1
#			h=0
            #siz=os.path.getsize(fn)
            #m=fallocate.posix_fadvise(iid, 0, 0,fallocate.POSIX_FADV_DONTNEED) 
#            grib_index_release(iid)
            f.close()

            if d/md not in [2,4]:
                print('grib file incomplete, skipped')
                continue


            try:
                hi=[]
                hw=[]
                dmd=d//md
                for i in range(0,dmd):
                    for j in range(24//dmd):
                        hi.append((i,(i+1)%dmd,-((i+1)//dmd)))
                        hw.append((1.0-j*dmd/24.,j*dmd/24.))
                hw=numpy.asarray(hw)
                hi=numpy.asarray(hi,dtype='int')

                print('after grib:',time.time()-t2)
                if varname2 is not None:
                    sgrid_to_station(tgrid,tuple(svar2),tuple(stemp),xi,yi,
                                     tuple(lons),
                                                tuple(lats),
                                                tuple(sweights),hi,hw,
                                                tuple(sdatum),
                                                tuple(shours),
                                                tx,
                                                md,numpy.int(1),numpy.int(0),numpy.int(fcstep))

                if corr is not None and iy>=cstart and iy<cstop:
                    fpd=fieldsperday
                    tgrid[0:fpd*md:fpd,:,:,:]-=corr[im-1,0,:,:,:] # assumed averaged over years
                    if fpd==2:
                        tgrid[1:fpd*md:fpd,:,:,:]-=corr[im-1,1,:,:,:] # assumed averaged over years
                    else:
                        tgrid[2:fpd*md:fpd,:,:,:]-=corr[im-1,1,:,:,:] # assumed averaged over years
                        c=0.5*(corr[im-1,0,:,:,:]+corr[im-1,1,:,:,:])
                        tgrid[1:fpd*md:fpd,:,:,:]-=c # assumed averaged over years
                        tgrid[3:fpd*md:fpd,:,:,:]-=c # assumed averaged over years
                    print('after correction:',time.time()-t2)

                sgrid_to_station(tgrid,tuple(svar),tuple(stemp),xi,yi,
                                 tuple(lons),
                                    tuple(lats),
                                    tuple(sweights),hi,hw,
                                    tuple(sdatum),
                                    tuple(shours),
                                    tx,
                                    md,numpy.int(1),numpy.int(0),numpy.int(fcstep))
                for syll in sylls:
                    for st in list(stats.keys()):
                        if st in ['odbfile', 'header'] or 'erai_'+syll+'temperatures' not in list(stats[st].keys()):
                            continue

                    #if 'erai_'+syll+'temperatures' in stats[st].keys():
                        grid_to_gpsstation(tgrid,stats[st]['erai_'+syll+'temperatures'],stats[st]['m'+tempname],xi,yi,
                                           stats[st][syll+'lon'],
                                    stats[st][syll+"lat"],
                                    stats[st]["weights"],hi,hw,
                                    stats[st][syll+'datum'],
                                    stats[st]['mdatum'],
                                    stats[st][syll+'hours'],
                                    tx,
                                    md,numpy.int(1),numpy.int(0),numpy.int(fcstep))

                        #print 'gs2',time.time()-tt
                        mi=numpy.where(numpy.logical_and(stats[st]['mdatum']>=tx,stats[st]['mdatum']<tx+md ))[0]
                        gi=numpy.where(numpy.logical_and(stats[st][syll+'datum']>=tx,stats[st][syll+'datum']<tx+md ))[0]

                        #print len(mi),len(gi),time.time()-tt
                        l=0
                        lg=[]
                        sgd=stats[st][syll+'datum']
                        smd=stats[st]['mdatum']
                        if len(gi)>0:
                            for i in range(mi.shape[0]):
                                while l<gi.shape[0]-1 and sgd[gi[l]]<smd[mi[i]]:
                                    l+=1
                                if sgd[gi[l]]==smd[mi[i]] and numpy.sum(~numpy.isnan(stats[st][syll+'temperatures'][:,:,gi[l]]))>0:
                                    stats[st]['erai_fg'+syll+'dep'][:,:,mi[i]]=stats[st]['erai_'+syll+tempname][:,:,gi[l]]-stats[st][syll+'temperatures'][:,:,gi[l]]
                                    lg.append(i)


                        #plt.figure(figsize=(12,5))
                        #plt.subplot(1,2,1)
                        #plt.plot(smd[mi],stats[st]['erai_fgtemperatures'][0,5,mi],label='erai_fg')
                        #plt.plot(smd[mi],stats[st]['mtemperatures'][0,5,mi],label='rsobs')
                        #plt.plot(sgd[gi],stats[st][syll+'temperatures'][0,5,gi],'x',label=syll+'ret')
                        #plt.plot(sgd[gi],stats[st]['erai_'+syll+'temperatures'][0,5,gi],'+',label='erai_fg'+syll+'')
                        #plt.xlabel('days since 19000101')
                        #plt.ylabel('Temperature [K]')
                        #plt.legend()
                        #plt.subplot(1,2,2)
                        #plt.plot(smd,stats[st]['erai_fg'+syll+'dep'][0,5,:],'x',label=syll+'-erai_fg'+syll+'')
                        #plt.plot(smd,stats[st]['erai_fgtemperatures'][0,5,:]-stats[st]['mtemperatures'][0,5,:],label='obs-erai_fg')
                        #plt.xlabel('days since 19000101')
                        #plt.ylabel('Temperature [K]')
                        #plt.legend()

                        #maskg=stats[st]['erai_'+syll+'temperatures'][0,5,gi[0]:gi[-1]+1]==stats[st]['erai_'+syll+'temperatures'][0,5,gi[0]:gi[-1]+1]
                        #print 'erai'+syll+' vs '+syll+'',sum(maskg),numpy.corrcoef(stats[st]['erai_'+syll+'temperatures'][0,5,gi[0]:gi[-1]+1][maskg],
                                                #stats[st][syll+'temperatures'][0,5,gi[0]:gi[-1]+1][maskg])[0,1]
                        #maskm=stats[st]['erai_fgtemperatures'][0,5,mi[0]:mi[-1]+1]==stats[st]['erai_fgtemperatures'][0,5,mi[0]:mi[-1]+1]
                        #print 'eraifg vs obs',sum(maskm),numpy.corrcoef(stats[st]['erai_fgtemperatures'][0,5,mi[0]:mi[-1]+1][maskm],
                                                #stats[st]['mtemperatures'][0,5,mi[0]:mi[-1]+1][maskm])[0,1]
                        #maskgm=stats[st]['erai_fg'+syll+'dep'][0,5,mi[0]:mi[-1]+1]==stats[st]['erai_fg'+syll+'dep'][0,5,mi[0]:mi[-1]+1]
                        #print syll+'dep vs fgdep',sum(maskgm),numpy.corrcoef(stats[st]['erai_fg'+syll+'dep'][0,5,mi[0]:mi[-1]+1][maskgm],
                                                #stats[st]['erai_fgtemperatures'][0,5,mi[0]:mi[-1]+1][maskgm]-stats[st]['mtemperatures'][0,5,mi[0]:mi[-1]+1][maskgm])[0,1]

                        #aux=stats[st]['mtemperatures'][0,5,mi[0]:mi[-1]+1]
                        #aux=numpy.take(aux,lg)		    
                        #masktm=stats[st][syll+'temperatures'][0,5,gi[0]:gi[-1]+1]==stats[st][syll+'temperatures'][0,5,gi[0]:gi[-1]+1]
                        #print syll+' vs obs',sum(masktm),numpy.corrcoef(stats[st][syll+'temperatures'][0,5,gi[0]:gi[-1]+1][masktm],
                                                #aux[masktm])[0,1]
                        #aux=stats[st]['erai_fgtemperatures'][0,5,mi[0]:mi[-1]+1]
                        #aux=numpy.take(aux,lg)
                        #masktm=stats[st][syll+'temperatures'][0,5,gi[0]:gi[-1]+1]==stats[st][syll+'temperatures'][0,5,gi[0]:gi[-1]+1]
                        #print 'jrafg vs jra'+syll+'',sum(masktm),numpy.corrcoef(stats[st]['erai_'+syll+'temperatures'][0,5,gi[0]:gi[-1]+1][masktm],
                                                #aux[masktm])[0,1]
                        #plt.show()
                        #print 'done'

                #if numpy.sum(stats[st][varname]==0.)>0 or numpy.nanmax(numpy.abs(stats[st][varname]))>1.e10:
                    #print 'spurious value',st,varname,tx,iy,im
                    #exit()
                print(iy*100+im)
                print(time.time()-t,time.time()-t2)
            except IOError:
                print(fn+' not found')
                pass
    gpsname=[]
    for syll in sylls:
        gpsname+=['erai_fg'+syll+'dep'] # return departure string since this departure is already finished

    if len(gpsname)>0:
        return gpsname
    else:
        if varname2 is not None:
            return [varname,varname2]
        else:
            return [varname]

@njit(cache=True)
def reduce_20cr(tin,tout,mask,stride):
    for it in range(0,tin.shape[0],stride):
        il=0
        ith=it/stride
        for ip in mask:
            if ip==1:
                for ilat in range(tin.shape[2]):
                    for ilon in range(tin.shape[3]):
                        tout[ith,il,ilat,ilon]=(tin[it,ip,ilat,ilon]+tin[it,ip+1,ilat,ilon])/2.  # 925 hPa
            else:
                for ilat in range(tin.shape[2]):
                    for ilon in range(tin.shape[3]):
                        tout[ith,il,ilat,ilon]=tin[it,ip,ilat,ilon]

            il+=1	
    return

class myThread (threading.Thread):
    def __init__(self, threadID, name, sdict,tdict,iy,im):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.sdict=sdict
        self.tdict=tdict
        self.iy=iy
        self.im=im
    def run(self):
        print("Starting " + self.name)
        add_feedback_mt_core(self.threadID,self.sdict,self.tdict,self.iy,self.im)
        print("Exiting " + self.name)

def add_feedback_mt(ce20cpath,filepattern,varname,start,stop,tidx,stats,fieldsperday=4,fcstep=0,
                    mrange=list(range(1,13)),tgrid=None,sylls=[],tempname='temperatures',
                 corr=None,cstart=None,cstop=None,varname2=None):

    tdict=(ce20cpath,filepattern,varname,start,stop,tidx,stats,fieldsperday,fcstep,
           mrange,sylls,tempname,corr,cstart,cstop,varname2)

    t=time.time()
    slist=[]
    lats=[]
    lons=[]
    svar=[]
    svar2=[]
    stemp=[]
    sweights=[]
    sdatum=[]
    shours=[]
    for st in list(stats.keys()):
        if st in ['odbfile', 'header']:
            continue
        if 'm'+tempname in list(stats[st].keys()):
            if varname not in list(stats[st].keys()):
                stats[st][varname]=numpy.empty(stats[st]['m'+tempname].shape,dtype=numpy.float32)
                stats[st][varname].fill(numpy.nan)
            if varname2 is not None :
                if varname2 not in list(stats[st].keys()):
                    stats[st][varname2]=numpy.empty(stats[st]['m'+tempname].shape,dtype=numpy.float32)
                    stats[st][varname2].fill(numpy.nan)

            if type(stats[st]["lon"]) is list:
                lon=stats[st]["lon"][0]
                lat=stats[st]["lat"][0]
            else:
                lon=stats[st]["lon"]
                lat=stats[st]["lat"]
            if lat!=lat or lon!=lon:
                continue
        else:
            continue

        stemp.append(stats[st]['m'+tempname])
        svar.append(stats[st][varname])
        if varname2 is not None:
            svar2.append(stats[st][varname2])
        stats[st]['weights'][:]=0.
        sweights.append(stats[st]['weights'])
        sdatum.append(stats[st]['mdatum'])
        shours.append(stats[st]['mhours'])
        slist.append(st)
        lats.append(lat)
        lons.append(lon)


    il=-1
    ini=False
    tgrid=[None]*12
    tgrid2=[None]*12


    for iy in range(start,stop):
        mthread=[]
        if corr is not None:
            if iy<cstart or iy>=cstop:
                continue

        for im in mrange:
            if not ini:
                fn=ce20cpath+filepattern.format(iy,im)
                if not os.path.isfile(fn):
                    continue
                f=open(fn,'r')

                gid = grib_new_from_file(f)
                Ni=grib_get(gid,'Ni')
                Nj=grib_get(gid,'Nj')
                grib_release(gid)
                f.close()
                ini=True

            if tgrid[im-1] is None:
                tgrid[im-1]=numpy.empty((31*fieldsperday,16,Nj,Ni),dtype=numpy.float32)
                if varname2 is not None:
                    tgrid2[im-1]=numpy.empty((31*fieldsperday,16,Nj,Ni),dtype=numpy.float32)
                    ini=True

        for im in mrange:
            sdict=(slist,lats,lons,svar,svar2,stemp,sweights,sdatum,shours,tidx,tgrid,tgrid2)
            mthread.append(myThread(1, "Thread{}".format(im), sdict,tdict,iy,im))
            mthread[im-1].start()

        #for im in mrange:
        mthread[11].join()
        print(time.time()-t)
    gpsname=[]
    for syll in sylls:
        gpsname+=['erai_fg'+syll+'dep'] # return departure string since this departure is already finished

    if len(gpsname)>0:
        return gpsname
    else:
        if varname2 is not None:
            return [varname,varname2]
        else:
            return [varname]


def add_feedback_mt_core(threadid,sdict,tdict,iy,im):

    slist,lats,lons,svar,svar2,stemp,sweights,sdatum,shours,tidx,tgrid,tgrid2=sdict
    ce20cpath,filepattern,varname,start,stop,tidx,stats,fieldsperday,fcstep,mrange,sylls,tempname,corr,cstart,cstop,varname2=tdict    

    t2=time.time()
    #il+=1
    l=(iy-1900)*12+im-1
    tx=tidx[l]
    try:
        md=tidx[l+1]-tidx[l]
    except IndexError:
        md=31 # December assumed
#            fn=ce20cpath+'CERA20C{0}{1:02}.grb'.format(iy,im)
    fn=ce20cpath+filepattern.format(iy,im)
    print(fn)
    if not os.path.isfile(fn):
        return
    f=open(fn,'r')

    gid=0
    d=0
    h=0
    ip=0
    while True:
#		gid = grib_new_from_index(iid)
        gid = grib_new_from_file(f)
        if gid is None:
            break
        if ip==0 and h==0 and d==0:
            inc=-1+2*grib_get(gid,'jScansPositively')
            Ni=grib_get(gid,'Ni')
            Nj=grib_get(gid,'Nj')
            yi=grib_get(gid,'latitudeOfFirstGridPointInDegrees')+inc*numpy.arange(Nj)*grib_get(gid,'jDirectionIncrementInDegrees')
            xi=grib_get(gid,'longitudeOfFirstGridPointInDegrees')+numpy.arange(Ni)*grib_get(gid,'iDirectionIncrementInDegrees')

            x=grib_get_values(gid)
        add_togrid(tgrid[im-1],x,d,ip,Nj,Ni)

        grib_release(gid)
        ip+=1
        if ip==16:
            d+=1
            ip=0
    f.close()

    if d/md not in [2,4]:
        print('grib file incomplete, skipped')
        return


    try:
        hi=[]
        hw=[]
        dmd=d/md
        for i in range(0,dmd):
            for j in range(24/dmd):
                hi.append((i,(i+1)%dmd,-((i+1)/dmd)))
                hw.append((1.0-j*dmd/24.,j*dmd/24.))
        hw=numpy.asarray(hw)
        hi=numpy.asarray(hi,dtype='int')

        print('after grib:',time.time()-t2)
        if varname2 is not None:
            sgrid_to_station(tgrid[im-1],tuple(svar2),tuple(stemp),xi,yi,
                                     tuple(lons),
                                                tuple(lats),
                                                tuple(sweights),hi,hw,
                                                tuple(sdatum),
                                                tuple(shours),
                                                tx,
                                                md,numpy.int(1),numpy.int(0),numpy.int(fcstep))

        if corr is not None and iy>=cstart and iy<cstop:
            fpd=fieldsperday
            tgrid[im-1][0:fpd*md:fpd,:,:,:]-=corr[im-1,0,:,:,:] # assumed averaged over years
            if fpd==2:
                tgrid[im-1][1:fpd*md:fpd,:,:,:]-=corr[im-1,1,:,:,:] # assumed averaged over years
            else:
                tgrid[im-1][2:fpd*md:fpd,:,:,:]-=corr[im-1,1,:,:,:] # assumed averaged over years
                c=0.5*(corr[im-1,0,:,:,:]+corr[im-1,1,:,:,:])
                tgrid[im-1][1:fpd*md:fpd,:,:,:]-=c # assumed averaged over years
                tgrid[im-1][3:fpd*md:fpd,:,:,:]-=c # assumed averaged over years
            print('after correction:',time.time()-t2)

        sgrid_to_station(tgrid[im-1],tuple(svar),tuple(stemp),xi,yi,
                                 tuple(lons),
                                    tuple(lats),
                                    tuple(sweights),hi,hw,
                                    tuple(sdatum),
                                    tuple(shours),
                                    tx,
                                    md,numpy.int(1),numpy.int(0),numpy.int(fcstep))
        for syll in sylls:
            for st in list(stats.keys()):
                if st in ['odbfile', 'header'] or 'erai_'+syll+'temperatures' not in list(stats[st].keys()):
                    continue

            #if 'erai_'+syll+'temperatures' in stats[st].keys():
                grid_to_gpsstation(tgrid[im-1],stats[st]['erai_'+syll+'temperatures'],stats[st]['m'+tempname],xi,yi,
                                           stats[st][syll+'lon'],
                                    stats[st][syll+"lat"],
                                    stats[st]["weights"],hi,hw,
                                    stats[st][syll+'datum'],
                                    stats[st]['mdatum'],
                                    stats[st][syll+'hours'],
                                    tx,
                                    md,numpy.int(1),numpy.int(0),numpy.int(fcstep))

                #print 'gs2',time.time()-tt
                mi=numpy.where(numpy.logical_and(stats[st]['mdatum']>=tx,stats[st]['mdatum']<tx+md ))[0]
                gi=numpy.where(numpy.logical_and(stats[st][syll+'datum']>=tx,stats[st][syll+'datum']<tx+md ))[0]

                #print len(mi),len(gi),time.time()-tt
                l=0
                lg=[]
                sgd=stats[st][syll+'datum']
                smd=stats[st]['mdatum']
                if len(gi)>0:
                    for i in range(mi.shape[0]):
                        while l<gi.shape[0]-1 and sgd[gi[l]]<smd[mi[i]]:
                            l+=1
                        if sgd[gi[l]]==smd[mi[i]] and numpy.sum(~numpy.isnan(stats[st][syll+'temperatures'][:,:,gi[l]]))>0:
                            stats[st]['erai_fg'+syll+'dep'][:,:,mi[i]]=stats[st]['erai_'+syll+tempname][:,:,gi[l]]-stats[st][syll+'temperatures'][:,:,gi[l]]
                            lg.append(i)

        print(iy*100+im)
        print(time.time()-t2)
    except IOError:
        print(fn+' not found')
        pass

@njit(parallel=True)
def preduce_20cr(tin,tout,mask,stride):
    for it in prange(tin.shape[0]):
        reduce_20cr_map(tin,tout,mask,stride,it)   
    return

@njit(cache=True)
def reduce_20cr(tin,tout,mask,stride):
    for it in range(0,tin.shape[0],stride):
        il=0
        ith=it/stride
        for ip in mask:
            if ip==1:
                for ilat in range(tin.shape[2]):
                    for ilon in range(tin.shape[3]):
                        tout[ith,il,ilat,ilon]=(tin[it,ip,ilat,ilon]+tin[it,ip+1,ilat,ilon])/2.  # 925 hPa
            else:
                for ilat in range(tin.shape[2]):
                    for ilon in range(tin.shape[3]):
                        tout[ith,il,ilat,ilon]=tin[it,ip,ilat,ilon]

            il+=1	
    return

@njit(cache=True)
def reduce_20cr_map(tin,tout,mask,stride,it):
#    for it in range(0,tin.shape[0],stride):
    il=0
    ith=it/stride
    for ip in mask:
        if ip==1:
            for ilat in range(tin.shape[2]):
                for ilon in range(tin.shape[3]):
                    tout[ith,il,ilat,ilon]=(tin[it,ip,ilat,ilon]+tin[it,ip+1,ilat,ilon])/2.  # 925 hPa
        else:
            for ilat in range(tin.shape[2]):
                for ilon in range(tin.shape[3]):
                    tout[ith,il,ilat,ilon]=tin[it,ip,ilat,ilon]

        il+=1	
    return


#
def add_feedback_netcdf(path,filepattern,filepar,varname,start,stop,tidx,stats,st,sylls=[],tempname='temperatures',
                        corr=None,cstart=None,cstop=None,varname2=None,corr2=None):
    il=0
    t=time.time()
    values=None
    for iy in range(start,stop):
        mcheck=len(filepattern.split('{'))
        if mcheck==2:  # NOAA 20CR
            l=(iy-1900)*12
            tx=tidx[l]
            md=tidx[l+12]-tidx[l]
            fn=path+filepattern.format(iy)

            try:
                t2=time.time()
                f=netCDF4.Dataset(fn)
                #f=netcdf.netcdf_file(fn,'r')

                level=f.variables['level'][:]
                mask=numpy.asarray((0,1,3,6,10,12,14,15,16,17,18,19,20,21,22,23))
                vs=f.variables[filepar].shape
                if il==0:
                    tgrid=numpy.empty((366*4,mask.shape[0],vs[2],vs[3]),dtype=numpy.float32)
                    #values=numpy.empty((365*4,vs[1],vs[2],vs[3]),dtype=numpy.float32)
                    #values2=numpy.empty((366*4,vs[1],vs[2],vs[3]),dtype=numpy.float32)
                    #v32=numpy.empty((366*4,vs[1],vs[2],vs[3]),dtype=numpy.float32)
                    il=1

                #for it in range(vs[0]):
                #print time.time()-t2
                #if vs[0]==1460:
                    #values[:]=f.variables[filepar][:]#.astype(numpy.float32)
                    #preduce_20cr(values,tgrid,mask,1)
                #else:
                    #values2[:]=f.variables[filepar][:]#.astype(numpy.float32)
                    #preduce_20cr(values2,tgrid,mask,1)
                #print time.time()-t2
                preduce_20cr(f.variables[filepar][:],tgrid,mask,1)


                #reduce_20cr(f.variables[filepar][:],tgrid,mask,1)
                print(time.time()-t2)
#		tgrid[:,1,:,:]=(tgrid[:,1,:,:]+tgrid[:,2,:,:])/2.
#		tgrid=f.variables[filepar][:]
                xi=f.variables['lon'][:]
                yi=f.variables['lat'][:]
                print(time.time()-t2)
                hi=[]
                hw=[]
                dmd=f.variables[filepar].shape[0]/md
                for i in range(0,dmd):
                    for j in range(24/dmd):
                        hi.append((i,(i+1)%dmd,-((i+1)/dmd)))
                        hw.append((1.0-j*dmd/24.,j*dmd/24.))
                hw=numpy.asarray(hw)
                hi=numpy.asarray(hi,dtype='int')
                for st in list(stats.keys()):

                    stats[st]['weights'][:]=0.
                    if 'm'+tempname in list(stats[st].keys()):
                        pgrid_to_station(tgrid,stats[st][varname],stats[st]['m'+tempname],xi,yi,
                                         numpy.float(stats[st]["lon"]),
                                    numpy.float(stats[st]["lat"]),
                                    stats[st]["weights"],hi,hw,
                                    stats[st]['mdatum'],
                                    stats[st]['mhours'],
                                    tx,
                                    md,numpy.int(1),numpy.int(1),numpy.int(0))
                        #if numpy.sum(stats[st][varname]==0.)>0 or numpy.nanmax(numpy.abs(stats[st][varname]))>1.e10:
                            #print 'spurious value',st,varname,tx,iy,im
                            #exit()

                print(time.time()-t2)
                print(iy)
                f.close()
            except (RuntimeError,IOError) as e:
                print(fn+' not found')
                pass
        else:   # JRA55
            inc=1
            if mcheck==4:
                inc=3
            for im in range(1,13,inc):
                l=(iy-1900)*12+im-1
                tx=tidx[l]
                md=tidx[l+inc]-tidx[l]
                fnotfound=False
                for ih in range((inc+1)//2):
                    if mcheck==6:
                        fn=path+filepattern.format(iy,im,md)
                    else:
                        fn=path+filepattern.format(iy,im,ih*12)

                    try:                    
#			f=netCDF4.Dataset(fn,'r')
                        f=netcdf.netcdf_file(fn,'r')
#			f.set_auto_mask(False)
                    except (RuntimeError,IOError) as e:
                        print(fn+' not found')
                        fnotfound=True
                        continue
                    if il==0:
                        vs=f.variables[filepar].shape
                        tgrid=numpy.empty([(vs[0]+2)*2,vs[1],vs[2],vs[3]],dtype=numpy.float32)
                        if varname2 is not None:
                            tgrid2=numpy.empty([(vs[0]+2)*2,vs[1],vs[2],vs[3]],dtype=numpy.float32)
                        il=1
                    if mcheck==6:
                        miss_val=f.variables[filepar].getncattr('_FillValue')
                        tgrid[:]=f.variables[filepar][:]
                        tgrid[tgrid==miss_val]=numpy.nan
                        xi=f.variables['g0_lon_3'][:]
                        yi=f.variables['g0_lat_2'][:]
                        preverse=0
                    else:
                        tt=time.time()
                        miss_val=f.variables[filepar]._get_missing_value() #f.variables[filepar].getncattr('_FillValue')
                        if values is None:
                            s=f.variables[filepar].shape
                            s=(92,s[1],s[2],s[3])
                            values=numpy.empty(s,dtype=numpy.float32)#.astype(numpy.float32

                        s0=f.variables[filepar].shape[0]    
                        values[:s0,:,:,:]=f.variables[filepar][:]#.astype(numpy.float32)
                        print('readnetcdf',time.time()-tt)
                        preverse=0
                        tt=time.time()
                        if corr is not None:
                            if iy*100+im>=cstart and iy*100+im<cstop:
                                csplit(tgrid,values,ih,miss_val,corr,tidx[l:l+3]-tidx[l],im-1,corr2)
                            else:
                                csplit(tgrid,values,ih,miss_val,numpy.empty((1,1,1,1,1),dtype=numpy.float32),tidx[l:l+3]-tidx[l],im-1,corr2)
                        else:
                            csplit(tgrid,values,ih,miss_val,numpy.empty((1,1,1,1,1),dtype=numpy.float32),tidx[l:l+3]-tidx[l],im-1,corr2)

                        if varname2 is not None:
                            csplit(tgrid2,values,ih,miss_val,numpy.empty((1,1,1,1,1),dtype=numpy.float32),tidx[l:l+3]-tidx[l],im-1,corr2)			    
                        print('csplit',time.time()-tt)
                        #if corr is not None:
                            #if iy>=cstart and iy<cstop:
                                #esplit(tgrid,values,ih,miss_val,corr,tidx[l:l+3]-tidx[l],im-1)
                                #if varname2 is not None:
                                    #esplit(tgrid2,values,ih,miss_val,numpy.empty((1,1,1,1,1),dtype=numpy.float32),tidx,-1)			    
                        #else:
                            #esplit(tgrid,values,ih,miss_val,numpy.empty((1,1,1,1,1),dtype=numpy.float32),tidx,-1)	
                        #print 'esplit',time.time()-tt
                        #print numpy.nanstd(tgrid),numpy.nanstd(tgrid2)
                        #if corr is not None:
                            #if iy>=cstart and iy<cstop:
                                #dsplit(tgrid,values,ih,miss_val,corr,tidx[l:l+3]-tidx[l],im-1)
                                #if varname2 is not None:
                                    #dsplit(tgrid2,values,ih,miss_val,numpy.empty((1,1,1,1,1),dtype=numpy.float32),tidx,-1)			    
                        #else:
                            #dsplit(tgrid,values,ih,miss_val,numpy.empty((1,1,1,1,1),dtype=numpy.float32),tidx,-1)
                        #del values
                        #print 'dsplit',time.time()-tt
#			print numpy.nanstd(tgrid),numpy.nanstd(tgrid2)
                        xi=f.variables['lon'][:].astype(numpy.float32)
                        yi=f.variables['lat'][:].astype(numpy.float32)
                    #siz=os.path.getsize(fn)
                    #m=fallocate.posix_fadvise(f, 0, siz,fallocate.POSIX_FADV_DONTNEED) 
                    f.close()

                if fnotfound:
                    continue

                hi=[]
                hw=[]
                dmd=2
                if mcheck==6:
                    dmd=4
                for i in range(0,dmd):
                    for j in range(24//dmd):
                        hi.append((i,(i+1)%dmd,-((i+1)//dmd)))
                        hw.append((1.0-j*dmd/24.,j*dmd/24.))
                hw=numpy.asarray(hw)
                hi=numpy.asarray(hi,dtype='int')
                if hi[12,2]==-1:
                    hi[12,2]=0
                for st in list(stats.keys()):

                    if 'm'+tempname in list(stats[st].keys()):
                        tt=time.time()
                        pgrid_to_station(tgrid,stats[st][varname],stats[st]['m'+tempname],xi,yi,
                                         numpy.float(stats[st]["lon"]),
                                    numpy.float(stats[st]["lat"]),
                                    stats[st]["weights"],hi,hw,
                                    stats[st]['mdatum'],
                                    stats[st]['mhours'],
                                    tx,
                                    md,numpy.int(2-inc/3),preverse,numpy.int(0))
                        if varname2 is not None:
                            pgrid_to_station(tgrid2,stats[st][varname2],stats[st]['m'+tempname],xi,yi,
                                             numpy.float(stats[st]["lon"]),
                                    numpy.float(stats[st]["lat"]),
                                    stats[st]["weights"],hi,hw,
                                    stats[st]['mdatum'],
                                    stats[st]['mhours'],
                                    tx,
                                    md,numpy.int(2-inc/3),preverse,numpy.int(0))
                        #print 'gs1',time.time()-tt
                        for syll in sylls:
                            #if 'jra55_'+syll+tempname in stats[st].keys():
                            for st in list(stats.keys()):
                                if st in ['odbfile', 'header'] or 'erai_'+syll+'temperatures' not in list(stats[st].keys()):
                                    continue
                                grid_to_gpsstation(tgrid,stats[st]['jra55_'+syll+tempname],stats[st][syll+tempname],xi,yi,
                                                   stats[st][syll+"lon"],
                                            stats[st][syll+"lat"],
                                            stats[st]["weights"],hi,hw,
                                            stats[st][syll+'datum'],
                                            stats[st]['mdatum'],
                                            stats[st][syll+'hours'],
                                            tx,
                                            md,numpy.int(2-inc/3),preverse,numpy.int(0))
                                print('gs2',time.time()-tt)
                                mi=numpy.where(numpy.logical_and(stats[st]['mdatum']>=tx,stats[st]['mdatum']<tx+md ))[0]
                                gi=numpy.where(numpy.logical_and(stats[st][syll+'datum']>=tx,stats[st][syll+'datum']<tx+md ))[0]

                                print(len(mi),len(gi),time.time()-tt)
                                l=0
                                lg=[]
                                sgd=stats[st][syll+'datum']
                                smd=stats[st]['mdatum']
                                if len(gi)>0:
                                    for i in range(mi.shape[0]):
                                        while l<gi.shape[0]-1 and sgd[gi[l]]<smd[mi[i]]:
                                            l+=1
                                        if sgd[gi[l]]==smd[mi[i]]:
                                            stats[st]['jra55_fg'+syll+'dep'][:,:,mi[i]]=stats[st]['jra55_'+syll+tempname][:,:,gi[l]]-stats[st][syll+tempname][:,:,gi[l]]
                                            lg.append(i)

                                    #plt.figure(figsize=(12,5))
                                    #plt.subplot(1,2,1)
                                    #plt.plot(smd[mi],stats[st]['jra55_fgtemperatures'][0,5,mi],label='jra55_fg')
                                    #plt.plot(smd[mi],stats[st]['mtemperatures'][0,5,mi],label='rsobs')
                                    #plt.plot(sgd[gi],stats[st][syll+'temperatures'][0,5,gi],'x',label=syll+'ret')
                                    #plt.plot(sgd[gi],stats[st]['jra55_'+syll+'temperatures'][0,5,gi],'+',label='jra55_fg'+syll+'')
                                    #plt.xlabel('days since 19000101')
                                    #plt.ylabel('Temperature [K]')
                                    #plt.legend()
                                    #plt.subplot(1,2,2)
                                    #plt.plot(smd,stats[st]['jra55_fg'+syll+'dep'][0,5,:],'x',label=syll+'-jra55_fg'+syll+'')
                                    #plt.plot(smd,stats[st]['jra55_fgtemperatures'][0,5,:]-stats[st]['mtemperatures'][0,5,:],label='obs-jra55_fg')
                                    #plt.xlabel('days since 19000101')
                                    #plt.ylabel('Temperature [K]')
                                    #plt.legend()

                                    maskg=stats[st]['jra55_'+syll+tempname][0,5,gi[0]:gi[-1]+1]==stats[st]['jra55_'+syll+tempname][0,5,gi[0]:gi[-1]+1]
                                    print('jra55'+syll+' vs '+syll+'',sum(maskg),numpy.corrcoef(stats[st]['jra55_'+syll+tempname][0,5,gi[0]:gi[-1]+1][maskg],
                                                                                                stats[st][syll+tempname][0,5,gi[0]:gi[-1]+1][maskg])[0,1])
                                    maskm=stats[st]['jra55_fgtemperatures'][0,5,mi[0]:mi[-1]+1]==stats[st]['jra55_fgtemperatures'][0,5,mi[0]:mi[-1]+1]
                                    print('jra55fg vs obs',sum(maskm),numpy.corrcoef(stats[st]['jra55_fgtemperatures'][0,5,mi[0]:mi[-1]+1][maskm],
                                                                                     stats[st]['mtemperatures'][0,5,mi[0]:mi[-1]+1][maskm])[0,1])
                                    maskgm=stats[st]['jra55_fg'+syll+'dep'][0,5,mi[0]:mi[-1]+1]==stats[st]['jra55_fg'+syll+'dep'][0,5,mi[0]:mi[-1]+1]
                                    print(syll+'dep vs fgdep',sum(maskgm),numpy.corrcoef(stats[st]['jra55_fg'+syll+'dep'][0,5,mi[0]:mi[-1]+1][maskgm],
                                                                                         stats[st]['jra55_fgtemperatures'][0,5,mi[0]:mi[-1]+1][maskgm]-stats[st]['mtemperatures'][0,5,mi[0]:mi[-1]+1][maskgm])[0,1])

                                    aux=stats[st]['mtemperatures'][0,5,mi[0]:mi[-1]+1]
                                    aux=numpy.take(aux,lg)		    
                                    masktm=stats[st][syll+tempname][0,5,gi[0]:gi[-1]+1]==stats[st][syll+tempname][0,5,gi[0]:gi[-1]+1]
                                    print(syll+' vs obs',sum(masktm),numpy.corrcoef(stats[st][syll+tempname][0,5,gi[0]:gi[-1]+1][masktm],
                                                                                    aux[masktm])[0,1])
                                    aux=stats[st]['jra55_fgtemperatures'][0,5,mi[0]:mi[-1]+1]
                                    aux=numpy.take(aux,lg)
                                    masktm=stats[st][syll+tempname][0,5,gi[0]:gi[-1]+1]==stats[st][syll+tempname][0,5,gi[0]:gi[-1]+1]
                                    print('jrafg vs jra'+syll+'',sum(masktm),numpy.corrcoef(stats[st]['jra55_'+syll+tempname][0,5,gi[0]:gi[-1]+1][masktm],
                                                                                            aux[masktm])[0,1])
                                    #plt.show()
                                    #print 'done'



                                    #if numpy.sum(stats[st][varname]==0.)>0 or numpy.nanmax(numpy.abs(stats[st][varname]))>1000:
                                        #print 'spurious value',st,varname,tx,iy,im
                                        #exit()
                print(iy*100+im)
                print(time.time()-t)
    if varname2 is not None:
        return [varname,varname2]
    else:
        return [varname]

def refdiff_netcdf(path,filepattern,filepar,start,stop,tidx,tempname='temperatures'):
    il=0
    t=time.time()
    try:
        mdiff=numpy.load(path+filepar+'_{}-{}.npy'.format(start,stop))
    except:
        for iy in range(start,stop):
            mcheck=len(filepattern.split('{'))
            if mcheck==2:  # NOAA 20CR
                l=(iy-1900)*12
                tx=tidx[l]
                md=tidx[l+12]-tidx[l]
                fn=path+filepattern.format(iy)

                try:
                    t2=time.time()
                    f=netCDF4.Dataset(fn,'r')

                    level=f.variables['level'][:]
                    mask=numpy.asarray((0,1,3,6,10,12,14,15,16,17,18,19,20,21,22,23))
                    vs=f.variables[filepar].shape
                    if il==0:
                        tgrid=numpy.empty((366*4,mask.shape[0],vs[2],vs[3]),dtype=numpy.float32)
                        il=1

                    reduce_20cr(f.variables[filepar][:],tgrid,mask,1)
                    print(time.time()-t2)

                    xi=f.variables['lon'][:]
                    yi=f.variables['lat'][:]
                    print(time.time()-t2)
                    hi=[]
                    hw=[]
                    dmd=f.variables[filepar].shape[0]/md

                    print(time.time()-t2)
                    print(iy)
                    f.close()
                except (RuntimeError,IOError) as e:
                    print(fn+' not found')
                    pass
            else:   # JRA55
                inc=1
                if mcheck==4:
                    inc=3
                for im in range(1,13,inc):
                    l=(iy-1900)*12+im-1
                    tx=tidx[l]
                    md=tidx[l+inc]-tidx[l]
                    fnotfound=False
                    for ih in range((inc+1)//2):
                        if mcheck==6:
                            fn=path+filepattern.format(iy,im,md)
                        else:
                            fn=path+filepattern.format(iy,im,ih*12)

                        try:                    
                            f=netCDF4.Dataset(fn,'r')
                            f.set_auto_mask(False)
                        except (RuntimeError,IOError) as e:
                            print(fn+' not found')
                            fnotfound=True
                            continue
                        if il==0:
                            vs=f.variables[filepar].shape
                            tgrid=numpy.empty([(vs[0]+2)*2,vs[1],vs[2],vs[3]],dtype=numpy.float32)
                            mdiff=numpy.empty((stop-start,12,2,vs[1],vs[2],vs[3]),dtype=numpy.float32)
                            il=1
                        if mcheck==6:
                            miss_val=f.variables[filepar].getncattr('_FillValue')
                            tgrid[:]=f.variables[filepar][:]
                            tgrid[tgrid==miss_val]=numpy.nan
                            xi=f.variables['g0_lon_3'][:]
                            yi=f.variables['g0_lat_2'][:]
                            preverse=0
                        else:
                            miss_val=f.variables[filepar].getncattr('_FillValue')
                            corr2=numpy.zeros((12, 2, 16, 320, 640),dtype=numpy.float32)
                            csplit(tgrid,f.variables[filepar][:],ih,miss_val,numpy.empty((1,1,1,1,1),dtype=numpy.float32),tidx[l:l+3]-tidx[l],im-1,corr2)
                            xi=f.variables['lon'][:]
                            yi=f.variables['lat'][:]
                            preverse=1
                        f.close()

                    if fnotfound:
                        continue
                    ix=0
                    for imm in range(im-1,im+2):
                        tgrsh=numpy.reshape(tgrid,(tgrid.shape[0]//2,2,tgrid.shape[1],tgrid.shape[2],tgrid.shape[3]))
                        mdiff[iy-start,imm,:,:,:,:]=numpy.mean(tgrsh[tidx[l+ix]-tidx[l]:tidx[l+ix+1]-tidx[l],:,:,:,:],axis=0)
                        ix+=1

                    print(iy*100+im)
                    print(time.time()-t)
        numpy.save(path+filepar+'_{}-{}'.format(start,stop), mdiff)

    print('refdiff_netcdf',time.time()-t)	
    return mdiff

@njit(cache=True)
def prtest(k):
    #s=str(k)
    print(('In prtest '+'ok'))

#################################################################################
#  Main Script
#################################################################################

#prtest(2)

def retrieve_fb_jra55():

    quantilecheck(numpy.random.randn(2,16,23000))
    monthdays=[31,28,31,30,31,30,31,31,30,31,30,31]  
    start=1958
    stop=2014
    plevs=numpy.asarray([10,20,30,50,70,100,150,200,250,300,400,500,700,850,925,1000])
    nlev=plevs.shape[0]

    with open('/home/srvx7/leo/tables/Pub9volA160415x.flatfile') as f:
        d=f.read().split('\n')
    wmonrs=[]
    wmonames=[]
    for dd in d:
        ddd=dd.split('\t')
        try: 
            wmonrs.append(ddd[5])
            wmonames.append(ddd[7])
        except:
            pass



    meta=read_ERACLIM_UAInv(os.path.expanduser('~/fastscratch/giub')+'/Inventory_ERACLIM_upperair_2.0.txt')
    for s in range(len(meta['WMO'])):
        if len(meta['WMO'][s])>5:
            meta['WMO'][s]=meta['WMO'][s][:5]  # some WMO Nrs have several instances
        if int(meta['WMO'][s])!=-999:
            meta['WMO'][s]='0'+meta['WMO'][s]

    meta['belongsto']=[]
    for m in meta['unique_record_ID']:
        meta['belongsto'].append('')
    figra=open(os.path.expanduser('~/tables/igra.station.list'))
    igrad=figra.read().split('\n')
    igrameta=dict(wmonr=[],wmoname=[])
    for ist in igrad:
        try:
            igrameta['wmonr'].append('{0:0>6}'.format(numpy.int(ist[4:9])))
            igrameta['wmoname'].append(ist[12:32])
        except:
            pass

    fmerged = open('/home/srvx7/leo/tables/mergedstations.t', 'r')
    mergedwban=('003132','011003','010003')
    rd=fmerged.read()
    rd=rd.split('\n')
    for i in range(len(rd)):
        rd[i]=rd[i][:5]+' '+rd[i][5:]
    sts=[]
    lats=[]
    lons=[]
    stsWMO=[]
    stnames=[]
    for l in rd:
        st=l[10:32].split()
        try:
            if len(st[0])<5:
                st[0]='0'+st[0]
            if len(st[0])<6:
                st[0]='0'+st[0]

            if st[0] in mergedwban:
                continue
            sts.append(st[0])
            try:
                idx=wmonrs.index(st[0][1:])
                stnames.append(wmonames[idx])
            except:
                stnames.append('missing')
            print(st[0],stnames[-1])
            stsWMO.append(st[0])
            lats.append(float(st[1]))
            lons.append(float(st[2]))
        except:
            pass

    lorig=len(sts)
    l=0
    mc=copy.copy(meta['WMO'])
    for m in mc:
        if m!='-999' and len(m)<6:
            meta['WMO'][l]='0'+m
        l+=1
    l=0

    wbanids,wbanlats,wbanlons,wmoids,wbannames=read_wban(os.path.expanduser('~/tables')+'/WBAN.TXT-2006jan')
    presatfiles=glob.glob(os.path.expanduser('~/fastscratch/presat')+'/*/*_t.nc')
    wbanverified=[]
    presatlats=[]
    presatlons=[]
    presatnrs=[]
    for s in presatfiles:
        presatnrs.append(s.split('/')[-1][:-5])
        sx=s.split('/')[-1][1:-5]
        f=netCDF4.Dataset(s)
        llon=f.variables['lon'][0]
        llat=f.variables['lat'][0]
        presatlons.append(llon)
        presatlats.append(llat)
        iwban=-1
        iwmo=-1
        iwmonew=-1
        ists=-1
        try:
            iwban=wbanids.index(sx)
        except:
            pass
        try:
            iwmo=wmoids.index(sx)
        except:
            pass
        try:
            iwmonew=wmonrs.index(sx)
        except:
            pass
        try:
            ists=sts.index('0'+sx)
        except:
            pass
        wbansph=-1
        wmosph=-1
        if iwban>-1:
            wbansph,ksave=sphdist(numpy.asarray([wbanlats[iwban]]),llat ,numpy.asarray([wbanlons[iwban]]),llon)
        if iwmo>-1:
            wmosph,ksave=sphdist(numpy.asarray([wbanlats[iwmo]]),llat ,numpy.asarray([wbanlons[iwmo]]),llon)
        if iwban>-1 and wbansph<1.:
            print(sx,' : ',iwban,iwmo,iwmonew,ists,wbansph,wmosph,llat,llon,wbannames[iwban])
            wbanverified.append(s.split('/')[-1][:-5])

    with open('/raid60/scratch/leo/scratch/Date_and_Station_files/Stationfiles_presat/2prot') as f:
        wbansuspect=f.read().split()
    for l in range(len(wbansuspect)):
        wbansuspect[l]=wbansuspect[l][:-6]

    #for m in meta['WMO']:
        #if m!='-999' and m not in sts[:lorig]:
            #sts.append(m)
            #stsWMO.append(m)
            #stnames.append(meta['Stationname'][l])
            #lons.append('{0:6.2f}'.format(meta['Lon_DegE'][l]))
            #lats.append('{0:5.2f}'.format(meta['Lat_DegN'][l]))
            #print m,meta['Lat_DegN'][l],meta['Lon_DegE'][l],meta['Stationname'][l]
        #else:
            #if m in sts:
                #meta['belongsto'][l]=m

        #l+=1



    l=0
    lfound=0
    lorphan=0
    lnew=0
    for m in meta['WMO']:
    #    if m=='-999':
        if meta['Lon_DegE'][l]!=-999.:
            minsph,ksave=sphdist(numpy.asarray(lats),meta['Lat_DegN'][l] ,numpy.asarray(lons),meta['Lon_DegE'][l])
            if minsph>0.3 and m not in sts:
                if m!='-999':
                    sts.append(m)
                else:
                    sts.append('{}'.format(100000+int(meta['unique_record_ID'][l])))
                stnames.append(meta['Stationname'][l])
                lons.append(meta['Lon_DegE'][l])
                lats.append(meta['Lat_DegN'][l])
                print(100000+int(meta['unique_record_ID'][l]),lats[-1],lons[-1],'new ',meta['Stationname'][l],\
                      'near',sts[ksave],stnames[ksave],lats[ksave],lons[ksave],minsph)
                lnew+=1
            else:
                print(100000+int(meta['unique_record_ID'][l]),meta['Stationname'][l],'belongs to',sts[ksave],stnames[ksave])
                meta['belongsto'][l]=sts[ksave]
                if stnames[ksave]=='':
                    stnames[ksave]=meta['Stationname'][l]
                lfound+=1
        else:
            print('Orphan ',meta['unique_record_ID'][l],meta['Lon_DegE'][l],meta['Stationname'][l])
            lorphan+=1

        l+=1
    print(lfound,' found ',lnew,' new ', lorphan,' orphans')

    l=0
    lfound=0
    lorphan=0
    lnew=0
    for m in meta['WMO']:
    #    if m=='-999':
        if meta['Lon_DegE'][l]!=-999.:
            minsph,ksave=sphdist(numpy.asarray(lats),meta['Lat_DegN'][l] ,numpy.asarray(lons),meta['Lon_DegE'][l])
            if minsph>0.3 and m not in sts:
                if m!='-999':
                    sts.append(m)
                else:
                    sts.append('{}'.format(100000+int(meta['unique_record_ID'][l])))
                stnames.append(meta['Stationname'][l])
                lons.append(meta['Lon_DegE'][l])
                lats.append(meta['Lat_DegN'][l])
                print(100000+int(meta['unique_record_ID'][l]),lats[-1],lons[-1],'new ',meta['Stationname'][l],\
                      'near',sts[ksave],stnames[ksave],lats[ksave],lons[ksave],minsph)
                lnew+=1
            else:
                print(100000+int(meta['unique_record_ID'][l]),meta['Stationname'][l],'belongs to',sts[ksave],stnames[ksave])
                meta['belongsto'][l]=sts[ksave]
                if stnames[ksave]=='':
                    stnames[ksave]=meta['Stationname'][l]
                lfound+=1
        else:
            print('Orphan ',meta['unique_record_ID'][l],meta['Lon_DegE'][l],meta['Stationname'][l])
            lorphan+=1

        l+=1
    print(lfound,' found ',lnew,' new ', lorphan,' orphans')
    #exit()
    stats=dict()
    indexmax=45000
    tidx=calcdays(19000101,(2020-1900)*12)-1
    l=-1
    t=time.time()

    temperatures=numpy.empty([2,16,indexmax],dtype=numpy.float32)

    sei=glob.glob('/home/srvx7/leo/fastscratch/ei6/0*/0*_t.nc')
    sbelong=[]
    for s in sts:
        sbelong.append(None)
    ni=0
    mi=0
    for s in sei:
        sx=s.split('/')[-1][:-5]
        #if sx=='001001':
            #print sx
        if sx not in sts and os.path.getsize(s)>30000:
    #	print sx,'with size',os.path.getsize(s),'not in sts'
            f=netCDF4.Dataset(s)
            llat=f.variables['lat'][0]
            llon=f.variables['lon'][0]
            ll=0
            minsph,ksave=sphdist(numpy.asarray(lats),llat,numpy.asarray(lons),llon)
            if minsph<0.3:
                sbelong[ksave]=sx
                mi+=1
            else:
                print(sx,' new station!',os.path.getsize(s))
                sts.append(sx)
                lats.append(llat)
                lons.append(llon)
                if sx[1:] in wmonrs:
                    idx=wmonrs.index(sx[1:])
                    stnames.append(wmonames[idx])
                else:
                    stnames.append('missing')
                sbelong.append(None)
                ni+=1
    print(ni,mi)
    nm=0
    ffns=open('mergeinputfiles','w')
    for st in sts:
        l+=1
        #if st!='001001' :
            #continue
        stsave=copy.copy(st)
        lat=float(lats[l])
        lon=float(lons[l])
        stats[st]={'lat':lat,'lon':lon,'weights':numpy.zeros(4)}
        fns=['/home/srvx7/leo/fastscratch/ei6/'+st+'/'+st+'_t.nc',
#	    '/home/srvx7/leo/scratch/ei/'+st[1:6]+'/feedbackmerged'+st[1:6]+'.nc',
]
        try:
            minsphs,ksave=sphdistall(numpy.asarray(presatlats),lat,numpy.asarray(presatlons),lon)
            mm=0
            nm=0
            for m in minsphs:
                if m<0.3:
                    fns.append('/home/srvx7/leo/fastscratch/presat/'+presatnrs[mm]+'/'+presatnrs[mm]+'_t.nc')
                    nm+=1
                mm+=1
            if nm>0:
                print(st,'appended',nm,'presat stations')
        except:
            pass

        if type(sbelong[l]) is str:
            fns.append('/home/srvx7/leo/fastscratch/ei6/'+sbelong[l]+'/'+sbelong[l]+'_t.nc')
        ll=0
        #if st=='007190':
            #print 'x'
        for a,b,c,d,e in zip(meta['WMO'],meta['unique_record_ID'],meta['Lat_DegN'],meta['Lon_DegE'],meta['belongsto']):
            if a==st or e==st:
                fns.append('/home/srvx7/leo/fastscratch/CHUAN/{}/feedbackmerged{}pl.nc'.format(100000+int(b),100000+int(b)))
                fns.append('/home/srvx7/leo/fastscratch/CHUAN/{}A/feedbackmerged{}Apl.nc'.format(100000+int(b),100000+int(b)))
                fns.append('/home/srvx7/leo/fastscratch/CHUAN/{}B/feedbackmerged{}Bpl.nc'.format(100000+int(b),100000+int(b)))
                print('appended',st,a,b,lat,lon,meta['Lat_DegN'][ll],meta['Lon_DegE'][ll],e)
                try:
                    f=netCDF4.Dataset(fns[-3],'r')
                    print(st,f.variables['lat'][0],f.variables['lon'][0],lat,lon,a,b,e)
                    print('')
                except:
                    pass
            else:
                if numpy.abs(lats[l]-c)<0.3 and numpy.abs(lons[l]-d)<1.0:
                    minsph,ksave=sphdist(numpy.asarray([lats[l]]),c,numpy.asarray([lons[l]]),d)
                    if minsph<0.3:
                        #fns.append('/home/srvx7/leo/fastscratch/CHUAN/{}/feedbackmerged{}pl.nc'.format(100000+int(b),100000+int(b)))
                        #fns.append('/home/srvx7/leo/fastscratch/CHUAN/{}A/feedbackmerged{}Apl.nc'.format(100000+int(b),100000+int(b)))
                        #fns.append('/home/srvx7/leo/fastscratch/CHUAN/{}B/feedbackmerged{}Bpl.nc'.format(100000+int(b),100000+int(b)))
                        print('small difference',a,st,b,c,float(lats[l]),d,float(lons[l]),e,minsph)
                        nm+=1
            ll+=1
    #    fns=['/home/srvx7/leo/fastscratch/ei6/'+st+'/feedbackmerged'+st+'.nc',]

        ffns.write('\n')
    #    if st=='007190':
        for fn in fns:
    #	if 'CHUAN' in fn:
    #	    print 'CHUAN'
            read_temperature(stats[st],st,fn,tidx,temperatures,indexmax)
            ffns.write(st+' : '+fn+'\n')    

        if 'mtemperatures' not in list(stats[st].keys()):
            print(st,': no temperature data found')
            continue
    #    else:
    #	print st,': some temperature data found'


        fnsurfs=[
            ('e20c','/home/srvx7/leo/fastscratch/e20c/'+st[1:6]+'/feedbackmerged'+st[1:6]+'.nc'),
            ('e20c','/home/srvx7/leo/fastscratch/presat/'+st+'/'+st+'_20c.nc'),
            ('n20cr','/home/srvx7/leo/fastscratch/n20cr/'+st+'/feedbackmerged'+st+'.nc'),
            ('n20cr','/home/srvx7/leo/fastscratch/presat/'+st+'/'+st+'_n20cr.nc')
        ]
        for fn in fnsurfs:
            try:
                f=netCDF4.Dataset(fn[1],'r')
                f.set_auto_mask(False)
                s=getattr(f.variables['datum'],'units').split()[2].split('-')[0]
                shift=tidx[(int(s)-1900)*12]
                miss_val=getattr(f.variables['e20c_0'],'missing_value')
    #            temperatures=numpy.empty([2,stats[st]['mtemperatures'].shape[1],indexmax],dtype=numpy.float32)
                temperatures[:]=numpy.nan
                vals=f.variables['e20c_0'][:]  
                vals[vals==miss_val]=numpy.nan
                merge=False
                try:
                    temperatures[:,:,stats[st]['mdatum']]=stats[st][fn[0]+'_andep']
                    merge=True
                except:
                    fnold=fn
                    pass
                temperatures[:,:,f.variables['datum'][0,:]+shift]=-vals #obs-bg
                stats[st][fn[0]+'_andep']=temperatures[:,:,stats[st]['mdatum']]
                if merge:
                    fnold=''
            except:
    #            try:
    #                temperatures=numpy.empty([2,stats[st]['mtemperatures'].shape[1],indexmax],dtype=numpy.float32)
    #            except:
                temperatures[:]=numpy.nan


        if 'mtemperatures' in list(stats[st].keys()):
            for stkey in ['jra55_antemperatures','jra55_fgtemperatures','erai_fgtemperatures',
                          'e20c_antemperatures','n20c_antemperatures','erapresat_antemperatures']:
                stats[st][stkey]=numpy.empty(stats[st]['mtemperatures'].shape,dtype=numpy.float32)
                stats[st][stkey][:]=numpy.nan
            for ie in range(10):
                stats[st]['ce20c{0}_antemperatures'.format(ie)]=numpy.empty(stats[st]['mtemperatures'].shape)
                stats[st]['ce20c{0}_antemperatures'.format(ie)][:]=numpy.nan

            # '+syll+' RO temperatures
            for syll in ['gpswet','gps']: #'gps'
                gpsi='/home/srvx7/leo/fastscratch/ei6/'+st+'/'+st+'_'+syll+'.nc'
                if os.path.isfile(gpsi):
                    merge_gps(gpsi,stats[st])
                else:
                    stats[st]['jra55_fg'+syll+'dep']=numpy.empty(stats[st]['mtemperatures'].shape,dtype=numpy.float32)
                    stats[st]['jra55_fg'+syll+'dep'][:]=numpy.nan
                    stats[st]['erai_fg'+syll+'dep']=numpy.empty(stats[st]['mtemperatures'].shape,dtype=numpy.float32)
                    stats[st]['erai_fg'+syll+'dep'][:]=numpy.nan

            print(st)

    #st='011035'
    ffns.close()       
    print(time.time()-t)
    nready=True ; fready=True ; aready=True; ce20cready=True; e20cready=True; ce20ceready=True; presatready=True; gpsready=True

    #print time.time()-t
    presatready=False
    if not presatready:
        add_feedback('/raid60/scratch/leo/scratch/ERApreSAT/','erapreSAT{0}{1:02}.130.grb','erapresat_antemperatures',1939,1967,tidx,stats)#1939,1967
        #stats[st]['ref']=numpy.empty(stats[st]['erapresat_antemperatures'].shape)
        #stats[st]['ref'][:]=stats[st]['erapresat_antemperatures'][:]

    #for st in sts:
        #if 'mtemperatures' in stats[st].keys() and 'erapresat_antemperatures' in stats[st].keys():
            #add_latlonfeedback(os.path.expandvars('$FSCRATCH/ei6/gridded'),'erapresat_antemperatures',stats,st,fcstep=0)
    #plt.plot(stats[st]['mdatum']/365.+1900,(stats[st]['erapresat_antemperatures']-stats[st]['mtemperatures'])[0,5,:])
    #plt.plot(stats[st]['mdatum']/365.+1900,(stats[st]['erapresat_antemperatures']-stats[st]['ref'])[0,5,:],'r-')
    #print 'diff:'
    print(time.time()-t)
    ce20cready=False
    if not ce20cready:
        add_feedback('/raid60/scratch/leo/scratch/CERA20C_ens/','CERA20C{0}{1:02}.130.0.grb','ce20c0_antemperatures',1930,2011,tidx,stats)#1930,2011

    print(time.time()-t)
    ce20ceready=False
    if not ce20ceready:
        for ie in range(1,10):
            add_feedback('/raid60/scratch/leo/scratch/CERA20C_ens/','CERA20C{0}{1:02}.130.'+'{0}'.format(ie)+'.grb',
                         'ce20c'+'{0}'.format(ie)+'_antemperatures',1930,2011,tidx,stats)#1930,2011

#tidx,stats,fieldsperday=4,fcstep=0,mrange=list(range(1,13)),tgrid=None,syll='gps'
    e20cready=False
    if not e20cready:
        add_feedback('/raid60/scratch/leo/scratch/ERA20C/','era20C.{0}{1:02}.130','e20c_antemperatures',1930,2011,tidx,stats)#1930,2011

    nready=False
    if not nready:
        add_feedback_netcdf('/raid60/scratch/leo/scratch/N20CV2c/','air.{0}.nc','air','n20c_antemperatures',1930,2015,tidx,stats,st)#1930,2015


    aready=False
    if not aready:
        add_feedback_netcdf('/raid60/scratch/leo/scratch/JRA55/','anl_p125.011_tmp.{0}{1:02}0100_{0}{1:02}{2}18.haimberger105010.nc',
                            'TMP_GDS0_ISBL','jra55_antemperatures',1958,2014,tidx,stats,st)#1958,2014

    fready=False
    if not fready:
        add_feedback_netcdf('/raid60/scratch/leo/scratch/JRA55/split/','test.fcst_mdl.011_tmp.reg_tl319.{0}{1:02}.{2:02}.grb.nc',
                            'TMP','jra55_fgtemperatures',1958,2017,tidx,stats,st)#1958,2017

    gpsready=False
    if not gpsready:
        for syll in ['gps','gpswet']:
            add_feedback('/raid60/scratch/leo/scratch/ERAI/','eraifc.{0}{1:02}.130',
                         'erai_fgtemperatures',2001,2017,tidx,stats,fieldsperday=2,syll=syll)#2001,2017
            add_feedback_netcdf('/raid60/scratch/leo/scratch/JRA55/split/','test.fcst_mdl.011_tmp.reg_tl319.{0}{1:02}.{2:02}.grb.nc',
                                'TMP','jra55_fgtemperatures',2001,2017,tidx,stats,st,syll=syll)#2001,2017

    print(time.time()-t)
    id0=daysbetween(19000101,'days since 1939-01-01 00:00:00')
    id1=daysbetween(19000101,'days since 1967-01-01 00:00:00')
    id2=daysbetween(19000101,'days since 1979-01-01 00:00:00')
    id2a=daysbetween(19000101,'days since 1977-01-01 00:00:00')
    id2b=daysbetween(19000101,'days since 1981-01-01 00:00:00')

    jd1=0
    jd2=id1
    jd2a=daysbetween(19000101,'days since 1965-01-01 00:00:00')
    jd2b=daysbetween(19000101,'days since 1969-01-01 00:00:00')

    for st in list(stats.keys()):
        if 'mtemperatures' in list(stats[st].keys()):
            flag=False
            if aready:
                fno='/home/srvx7/leo/fastscratch/ei6/'+st+'/feedbackmerged'+st+'.nc'
                try:
                    fo = netCDF4.Dataset(fno,"r")
                    stats[st]['omdatum']=fo.variables['datum'][:]
                    stats[st]['jra55_andep']=fo.variables['jra55_andep'][:]
                except:    
                    stats[st]['jra55_andep']=stats[st]['mtemperatures']-numpy.nan
    #                flag=True
                    pass
                try:
                    fo.close()
                except:
                    pass
            else:
                stats[st]['jra55_andep']=stats[st]['jra55_antemperatures']-stats[st]['mtemperatures']

            if fready:
                fno='/home/srvx7/leo/fastscratch/ei6/'+st+'/feedbackmerged'+st+'.nc'
                try:
                    fo = netCDF4.Dataset(fno,"r")
                    stats[st]['omdatum']=fo.variables['datum'][:]
                    stats[st]['jra55_fgdep']=fo.variables['jra55_fgdep'][:]
                except:    
                    stats[st]['jra55_fgdep']=stats[st]['mtemperatures']-numpy.nan
    #                flag=True
                    pass
                try:
                    fo.close()
                except:
                    pass
            else:
                stats[st]['jra55_fgdep']=stats[st]['jra55_fgtemperatures']-stats[st]['mtemperatures']

            if ce20cready:
                fno='/home/srvx7/leo/fastscratch/ei6/'+st+'/feedbackmerged'+st+'.nc'
                try:
                    fo = netCDF4.Dataset(fno,"r")
                    stats[st]['omdatum']=fo.variables['datum'][:]
                    stats[st]['ce20c0_andep']=fo.variables['ce20c0_andep'][:]
                except:    
                    stats[st]['ce20c0_andep']=stats[st]['mtemperatures']-numpy.nan
    #		flag=True
                    pass
                try:
                    fo.close()
                except:
                    pass
            else:
                stats[st]['ce20c0_andep']=stats[st]['ce20c0_antemperatures']-stats[st]['mtemperatures']

            if ce20ceready:
                fno='/home/srvx7/leo/fastscratch/ei6/'+st+'/feedbackmerged'+st+'.nc'
                try:
                    fo = netCDF4.Dataset(fno,"r")
                    stats[st]['omdatum']=fo.variables['datum'][:]
                    for ie in range(1,10):
                        stats[st]['ce20c{0}_andep'.format(ie)]=fo.variables['ce20c{0}_andep'.format(ie)][:]
                except:    
                    for ie in range(1,10):
                        stats[st]['ce20c{0}_andep'.format(ie)]=stats[st]['mtemperatures']-numpy.nan
    #		flag=True
                    pass
                try:
                    fo.close()
                except:
                    pass
            else:
                for ie in range(1,10):
                    stats[st]['ce20c{0}_andep'.format(ie)]=stats[st]['ce20c{0}_antemperatures'.format(ie)]-stats[st]['mtemperatures']

            if presatready:
                fno='/home/srvx7/leo/fastscratch/ei6/'+st+'/feedbackmerged'+st+'.nc'
                try:
                    fo = netCDF4.Dataset(fno,"r")
                    stats[st]['omdatum']=fo.variables['datum'][:]
                    stats[st]['erapresat_andep']=fo.variables['erapresat_andep'][:]
                    mask=numpy.logical_and(numpy.isnan(stats[st]['mfg_dep']),~numpy.isnan(stats[st]['erapresat_andep']))
                    ms=numpy.sum(mask)
                    if ms>0:
                        stats[st]['mfg_dep'][mask]=stats[st]['erapresat_andep'][mask]
                        print(time.time()-t)
                        print(st, 'late merge using preSAT an',ms,'values')
                except:    
                    stats[st]['erapresat_andep']=stats[st]['mtemperatures']-numpy.nan
    #		flag=True
                    pass
                try:
                    fo.close()
                except:
                    pass
            else:
                stats[st]['erapresat_andep']=stats[st]['erapresat_antemperatures']-stats[st]['mtemperatures']
                t=time.time()
                mask=numpy.logical_and(numpy.isnan(stats[st]['mfg_dep']),~numpy.isnan(stats[st]['erapresat_andep']))
                ms=numpy.sum(mask)
                if ms>0:
                    stats[st]['mfg_dep'][mask]=stats[st]['erapresat_andep'][mask]
                    print(time.time()-t)
                    print(st, 'late merge using preSAT an',ms,'values')


            if nready:
                fno='/home/srvx7/leo/fastscratch/ei6/'+st+'/feedbackmerged'+st+'.nc'
                try:
                    fo = netCDF4.Dataset(fno,"r")
                    stats[st]['omdatum']=fo.variables['datum'][:]
                    stats[st]['n20c_andep']=fo.variables['n20c_andep'][:]
                except:    
                    stats[st]['n20c_andep']=stats[st]['mtemperatures']-numpy.nan
    #		flag=True
                    pass
                try:
                    fo.close()
                except:
                    pass
            else:
                stats[st]['n20c_andep']=stats[st]['n20c_antemperatures']-stats[st]['mtemperatures']

            if e20cready:
                fno='/home/srvx7/leo/fastscratch/ei6/'+st+'/feedbackmerged'+st+'.nc'
                try:
                    fo = netCDF4.Dataset(fno,"r")
                    stats[st]['omdatum']=fo.variables['datum'][:]
                    stats[st]['e20c_andep']=fo.variables['e20c_andep'][:]
                except:    
                    stats[st]['e20c_andep']=stats[st]['mtemperatures']-numpy.nan
    #		flag=True
                    pass
                try:
                    fo.close()
                except:
                    pass
            else:
                stats[st]['e20c_andep']=stats[st]['e20c_antemperatures']-stats[st]['mtemperatures']

            if gpsready:
                fno='/home/srvx7/leo/fastscratch/ei6/'+st+'/feedbackmerged'+st+'.nc'
                try:
                    fo = netCDF4.Dataset(fno,"r")
                    stats[st]['omdatum']=fo.variables['datum'][:]
                    stats[st]['erai_fggpswetdep']=fo.variables['erai_fggpswetdep'][:]
                    stats[st]['erai_fggpsdep']=fo.variables['erai_fggpsdep'][:]
                except:    
                    stats[st]['erai_fggpsdep']=stats[st]['erai_fggpsdep']-numpy.nan
                    stats[st]['erai_fggpswetdep']=stats[st]['erai_fggpswetdep']-numpy.nan
    #		flag=True
                    pass
                try:
                    fo.close()
                except:
                    pass
            else:
                pass
                #stats[st]['e20c_andep']=stats[st]['e20c_antemperatures']-stats[st]['mtemperatures']


            fn='/home/srvx7/leo/fastscratch/ei6/'+st+'/'+st+'_t.nc'
            try:
                f = netCDF4.Dataset(fn,"r")
            except (RuntimeError,IOError) as e:
                try:
                    fn='/home/srvx7/leo/fastscratch/ei/'+st[1:6]+'/feedbackmerged'+st[1:6]+'.nc'
                    f = netCDF4.Dataset(fn,"r")
                except (RuntimeError,IOError) as e:
                    try:
                        fn='/home/srvx7/leo/fastscratch/presat/'+st+'/'+st+'_t.nc'
                        f = netCDF4.Dataset(fn,"r")
                    except (RuntimeError,IOError) as e:
                        try:
        #		    fn='/home/srvx7/leo/fastscratch/ei/01001/feedbackmerged01001.nc'
                            fn='/vgc/srvx7/leo/scratch/CHUAN/'+st+'/feedbackmerged'+st+'pl.nc'
                            f = netCDF4.Dataset(fn,"r")
                        except (RuntimeError,IOError) as e:
                            flag=True
                            continue
    #        if flag:
    #            continue

            fno='/home/srvx7/leo/fastscratch/ei6/'+st+'/feedbackmerged'+st+'.nc'
            try:
                fo = netCDF4.Dataset(fno,"w")
            except:
                os.mkdir('/home/srvx7/leo/fastscratch/ei6/'+st)
                fo = netCDF4.Dataset(fno,"w")

            for i in f.ncattrs():
                setattr(fo,i,getattr(f,i))

            try:
                setattr(fo,'Stationnname',stnames[sts.index(st)])
            except:
                setattr(fo,'Stationname','unknown')

            for i in list(f.dimensions.keys()):
                if i=='time':
                    fo.createDimension(i,stats[st]["mdatum"].shape[0])
                else:
                    try:
                        fo.createDimension(i,len(f.dimensions[i]))
                    except:
                        flag=True
                        continue
            if flag:
                continue

            vlist=[]
            for ie in range(10):
                vlist.append('ce20c{0}_andep'.format(ie))
            for i in list(f.variables.keys()):
                var=f.variables[i]
                if i=='datum':
                    fo.createVariable(i,var.dtype,var.dimensions)
                    fo.variables[i][:]=stats[st]["mdatum"][:]
                elif i=='hours':
                    fo.createVariable(i,var.dtype,var.dimensions)
                    fo.variables[i][:]=stats[st]["mhours"][:]
                elif i=='temperatures':
                    fo.createVariable(i,var.dtype,var.dimensions)
                    fo.variables[i][:]=stats[st]["mtemperatures"][:]
                elif i=='fg_dep':
                    fo.createVariable(i,var.dtype,var.dimensions)
                    if "mfg_dep" in list(stats[st].keys()):
                        fo.variables[i][:]=-stats[st]["mfg_dep"][:]
                    else:
                        fo.variables[i][:]=stats[st]["mtemperatures"][:]+numpy.nan
                elif i=='an_dep':
                    fo.createVariable(i,var.dtype,var.dimensions)
                    if "man_dep" in list(stats[st].keys()):
                        fo.variables[i][:]=-stats[st]["man_dep"][:]
                    else:
                        fo.variables[i][:]=stats[st]["mtemperatures"][:]+numpy.nan
                elif i in ('lon','lat'):
                    fo.createVariable(i,var.dtype,var.dimensions)
                    fo.variables[i][:]=stats[st][i]
                elif i in ('alt','press'):
                    fo.createVariable(i,var.dtype,var.dimensions)
                    fo.variables[i][:]=var[:]
                else:
                    if i!='fflags' and i!='bias' and i not in ['eijra_fgdep','jra55_fgdep','jra55_andep','e20c_andep','n20c_andep','erapresat_andep',
                                                               'jra55_fggpsdep','erai_fggpsdep','era5_fggpsdep','jra55_fggpswetdep','erai_fggpswetdep',
                                                               'era5_fggpswetdep']+vlist:
                        fo.createVariable(i,var.dtype,var.dimensions)

                for j in var.ncattrs():
                    if j!='_FillValue' and j!='scale_factor' and j!='add_offset':
                        if i!='fflags' and i!='bias':
                            setattr(fo.variables[i],j,getattr(var,j))
                            if i=='datum' and j=='units':
                                setattr(fo.variables[i],j,'days since 1900-01-01 00:00:00')

            temperatures[:]=numpy.nan
            temperatures[:,:,stats[st]["mdatum"][:]]=-stats[st]["mfg_dep"]

            # JRA55 -> ERA-Interim shift removal
            #        plt.figure(figsize=(12,5))
            #        plt.plot(stats[st]["mdatum"][:]/365.25+1900,rmeanw(-stats[st]["mfg_dep"][1,5,:],30))
            #	plt.savefig(st+'.eps')

            stats[st]['ce20cens_andep']=stats[st]['ce20c0_andep'][:]
            for ie in range(1,10):
                stats[st]['ce20cens_andep']+=stats[st]['ce20c{}_andep'.format(ie)][:]
            stats[st]['ce20cens_andep']/=10

            jrakey="ce20cens_andep"
            a,b,c,d=0,0,0,0
            try:
                a=numpy.where(stats[st]["mdatum"]>=id0)[0][0]
                b=numpy.where(stats[st]["mdatum"]>id2)[0][0]
                temperatures[:,:,stats[st]["mdatum"][a:b]]=-stats[st][jrakey][:,:,a:b]
                d=numpy.where(stats[st]["mdatum"]>id2b)[0][0]
                meanjra=numpy.nanmean(stats[st][jrakey][:,:,b:d],axis=2)
                meanei=numpy.nanmean(stats[st]["mfg_dep"][:,:,b:d],axis=2)
                count=numpy.sum(~numpy.isnan(stats[st]["mfg_dep"][:,:,b:d]+stats[st][jrakey][:,:,b:d]),axis=2)
                thresh=180
                for ipar in range(2):
                    for ip in range(stats[st][jrakey].shape[1]):
                        if count[ipar,ip]>thresh:            
                            temperatures[ipar,ip,stats[st]["mdatum"][:b]]+=(meanjra[ipar,ip]-meanei[ipar,ip])
            except:
                print(st,a,b, 'no data in this interval')

            stats[st]['eice20_fgdep']=-temperatures[:,:,stats[st]["mdatum"]]

    # JRA55 -> ERA-Interim shift removal
    #        plt.figure(figsize=(12,5))
    #        plt.plot(stats[st]["mdatum"][:]/365.25+1900,rmeanw(-stats[st]["mfg_dep"][1,5,:],30))
    #	plt.savefig(st+'.eps')
            temperatures[:]=numpy.nan
            temperatures[:,:,stats[st]["mdatum"][:]]=-stats[st]["mfg_dep"]

            jrakey="jra55_fgdep"
            a,b,c,d=0,0,0,0
            try:
                a=numpy.where(stats[st]["mdatum"]>=id1)[0][0]
                b=numpy.where(stats[st]["mdatum"]>id2)[0][0]
                temperatures[:,:,stats[st]["mdatum"][a:b]]=-stats[st][jrakey][:,:,a:b]
                d=numpy.where(stats[st]["mdatum"]>id2b)[0][0]
                meanjra=numpy.nanmean(stats[st][jrakey][:,:,b:d],axis=2)
                meanei=numpy.nanmean(stats[st]["mfg_dep"][:,:,b:d],axis=2)
                count=numpy.sum(~numpy.isnan(stats[st]["mfg_dep"][:,:,b:d]+stats[st][jrakey][:,:,b:d]),axis=2)
                thresh=180
                for ipar in range(2):
                    for ip in range(stats[st][jrakey].shape[1]):
                        if count[ipar,ip]>thresh:            
                            temperatures[ipar,ip,stats[st]["mdatum"][:b]]+=(meanjra[ipar,ip]-meanei[ipar,ip])
            except:
                print(st,a,b, 'no data in this interval')

    # ERA-Presat -> JRA55 shift removal
            a,b,c,d=0,0,0,0
            try:
                a=0
                b=numpy.where(stats[st]["mdatum"]>jd2)[0][0]
                c=numpy.where(stats[st]["mdatum"]>=jd2a)[0][0]
                meanjra=numpy.nanmean(stats[st]["mfg_dep"][:,:,c:b],axis=2)
                meanei=numpy.nanmean(stats[st][jrakey][:,:,c:b],axis=2)
                count=numpy.sum(~numpy.isnan(stats[st]["mfg_dep"][:,:,c:b]+stats[st][jrakey][:,:,c:b]),axis=2)
                thresh=180
                for ipar in range(2):
                    for ip in range(stats[st][jrakey].shape[1]):
                        if count[ipar,ip]>thresh:            
                            temperatures[ipar,ip,stats[st]["mdatum"][:b]]+=(meanjra[ipar,ip]-meanei[ipar,ip]) 
            except:
                print(st,c,b, 'no data in this interval')

            stats[st]['eijra_fgdep']=-temperatures[:,:,stats[st]["mdatum"]]

            for dep in ['erai_fggpswetdep','era5_fggpswetdep','jra55_fggpswetdep','eijra_fgdep','eice20_fgdep','ce20cens_andep','jra55_fgdep','jra55_andep','e20c_andep','n20c_andep','erapresat_andep',
                        'erai_fggpsdep','era5_fggpsdep','jra55_fggpsdep']+vlist:
                fo.createVariable(dep,numpy.float32,f.variables['temperatures'].dimensions,fill_value=numpy.nan)
                try:
    #                mask=numpy.isnan(stats[st][dep])
    #                stats[st][dep][mask]=-1.e30
                    fo.variables[dep][:]=stats[st][dep]
                except:
                    try:
                        if stats[st]['mdatum'].shape[0]>stats[st]['omdatum'].shape[1]:
                            fo.variables[dep][:,:,:stats[st]['omdatum'].shape[1]]=stats[st][dep]
                    except:
                        print(st,': could not write variable '+dep)

                for j in var.ncattrs():
                    if j!='_FillValue' and j!='scale_factor' and j!='add_offset':
                        setattr(fo.variables[i],j,getattr(var,j))





            fo.close()
            f.close()

    print(time.time()-t)




def add_feedback_i(ce20cpath,filepattern,varname,start,stop,tidx,stats,st):

#    ce20cpath='/nas/srvx8/leo/CERA20C/'
    t=time.time()
    il=-1
    ini=False
    for iy in range(start,stop):
        for im in range(1,13):
            t2=time.time()
            il+=1
            l=(iy-1900)*12+im-1
            tx=tidx[l]
            md=tidx[l+1]-tidx[l]
#            fn=ce20cpath+'CERA20C{0}{1:02}.grb'.format(iy,im)
            fn=ce20cpath+filepattern.format(iy,im)
#	    f=open(fn)

            if not os.path.isfile(fn):
                continue
            if os.path.isfile(fn+'.id'):
                iid = grib_index_read(fn+'.id')
            else:
                try:
                    iid = grib_index_new_from_file(fn,["paramId"])
                    grib_index_write(iid,fn+'.id')
                except:
                    print('index of ',fn,' could not be created, skipped')
                    continue

            try:
                l = grib_index_get(iid,'paramId')
            except GribInternalError:
                print('index of ',fn,' could not be read, skipped')
                continue

            grib_index_select(iid, 'paramId', '130')
            print('after gribindex:',time.time()-t2)
            gid=0
            d=0
            h=0
            ip=0
            while True:
                gid = grib_new_from_index(iid)
#		gid = grib_new_from_file(f)
                if gid is None:
                    break
                hh=grib_get(gid,'time')
                if hh%600!=0:
                    grib_release(gid)
                    continue
                if ip==0 and h==0 and d==0:
                    inc=-1+2*grib_get(gid,'jScansPositively')
                    Ni=grib_get(gid,'Ni')
                    Nj=grib_get(gid,'Nj')
                    yi=grib_get(gid,'latitudeOfFirstGridPointInDegrees')+inc*numpy.arange(Nj)*grib_get(gid,'jDirectionIncrementInDegrees')
                    xi=grib_get(gid,'longitudeOfFirstGridPointInDegrees')+numpy.arange(Ni)*grib_get(gid,'iDirectionIncrementInDegrees')
#		pid=grib_get(gid,'paramId')
#		if pid==130:
                if not ini:
                    tgrid=numpy.empty((31*4,16,Nj,Ni))
                    ini=True
#		tt=time.time()
                tgrid[d,ip,:,:]=numpy.reshape(grib_get_values(gid),(Nj,Ni))
#		print time.time()-t
                grib_release(gid)
                ip+=1
                if ip==16:
                    d+=1
                    ip=0
#		    if h==2:
#			d+=1
#			h=0
            #siz=os.path.getsize(fn)
            #m=fallocate.posix_fadvise(iid, 0, 0,fallocate.POSIX_FADV_DONTNEED) 
            grib_index_release(iid)

            if d/md not in [2,4]:
                print('grib file incomplete, skipped')
                continue


            try:
                hi=[]
                hw=[]
                dmd=d/md
                for i in range(0,dmd):
                    for j in range(24/dmd):
                        hi.append((i,(i+1)%dmd,-((i+1)/dmd)))
                        hw.append((1.0-j*dmd/24.,j*dmd/24.))
                hw=numpy.asarray(hw)
                hi=numpy.asarray(hi,dtype='int')

                print('after grib:',time.time()-t2)
                for st in list(stats.keys()):

                    if 'mtemperatures' in list(stats[st].keys()):
                        grid_to_station(tgrid,stats[st][varname],stats[st]['mtemperatures'],xi,yi,
                                        numpy.float(stats[st]["lon"]),
                                    numpy.float(stats[st]["lat"]),
                                    stats[st]["weights"],hi,hw,
                                    stats[st]['mdatum'],
                                    stats[st]['mhours'],
                                    tx,
                                    md,numpy.int(1),numpy.int(0))
                        #if numpy.sum(stats[st][varname]==0.)>0 or numpy.nanmax(numpy.abs(stats[st][varname]))>1.e10:
                            #print 'spurious value',st,varname,tx,iy,im
                            #exit()
                print(iy*100+im)
                print(time.time()-t,time.time()-t2)
            except IOError:
                print(fn+' not found')
                pass

if __name__ == '__main__':

    fn='/home/srvx7/leo/fastscratch/ei6/001001/001001_tx'+'.nc'
    #netcdf4speedtest(fn)
    retrieve_fb_jra55()
