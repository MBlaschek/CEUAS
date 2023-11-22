from numba import njit,jit,prange
import numpy as np
import math
#import matplotlib.pyplot as plt
from datetime import date
#from rasotools.anomaly import *
import datetime
import time
import os

miss_val=np.nan
# henderson moving average weights.
# first row=symmetric, other rows are increasingly asymmetric.
hmweights=np.asarray([
    [-0.01935,-0.02787,0,0.06549,0.14736,0.21434,0.24006,0.21434,0.14736,0.06549,0,-0.02787,-0.01935],
    [0,-0.01643,-0.02577,0.00127,0.06594,0.14698,0.21314,0.23803,0.21149,0.14368,0.06099,-0.00533,-0.03401],
    [0,0,-0.011,-0.02204,0.0033,0.06626,0.14559,0.21004,0.23324,0.20498,0.13547,0.05108,-0.01695],
    [0,0,0,-0.00814,-0.0202,0.00413,0.06608,0.14441,0.20784,0.23002,0.20076,0.13024,0.04483],
    [0,0,0,0,-0.01604,-0.02487,0.00267,0.06784,0.14939,0.21605,0.24144,0.2154,0.1481],
    [0,0,0,0,0,-0.04271,-0.03864,0.00182,0.0799,0.17436,0.25392,0.29223,0.2791],
    [0,0,0,0,0,0,-0.09187,-0.05812,0.01202,0.11977,0.2439,0.35315,0.42113]
])    


def Interp(datain,xin,yin,xout,yout,order=0):

    """
              Interpolates a 2D array onto a new grid (only works for linear grids), 
              with the Lat/Lon inputs of the old and new grid. Can perfom nearest
              neighbour interpolation or bilinear interpolation (of order 1)'

       This is an extract from the basemap module (truncated)
    """

    if order==0:
           interpolation='NearestNeighbour' 
    else:
           interpolation = 'Bilinear'
           
    # Mesh Coordinates so that they are both 2D arrays
 #   xout,yout = np.meshgrid(xout,yout)

   # compute grid coordinates of output grid.
    delx = xin[1:]-xin[0:-1]
    dely = yin[1:]-yin[0:-1]

    xcoords = (len(xin)-1)*(xout-xin[0])/(xin[-1]-xin[0])
    ycoords = (len(yin)-1)*(yout-yin[0])/(yin[-1]-yin[0])


    xcoords = np.clip(xcoords,0,len(xin)-1)
    ycoords = np.clip(ycoords,0,len(yin)-1)

    # Interpolate to output grid using nearest neighbour
    if interpolation == 'NearestNeighbour':
        xcoordsi = np.around(xcoords).astype(np.int32)
        ycoordsi = np.around(ycoords).astype(np.int32)
        dataout = datain[ycoordsi,xcoordsi]

    # Interpolate to output grid using bilinear interpolation.
    elif interpolation == 'Bilinear':
        xi = xcoords.astype(np.int32)
        yi = ycoords.astype(np.int32)
        xip1 = xi+1
        yip1 = yi+1
        xip1 = np.clip(xip1,0,len(xin)-1)
        yip1 = np.clip(yip1,0,len(yin)-1)
        delx = xcoords-xi.astype(np.float32)
        dely = ycoords-yi.astype(np.float32)
        dataout = (1.-delx)*(1.-dely)*datain[yi,xi] + \
                  delx*dely*datain[yip1,xip1] + \
                  (1.-delx)*dely*datain[yip1,xi] + \
                  delx*(1.-dely)*datain[yi,xip1]

    return dataout

    #@njit(cache=True,parallel=True)
def nncat(a,b,c):
#    for i in range(a.shape[0]):
#        for j in range(a.shape[1]):
#            for k in range(a.shape[2]):
#                for l in range(a.shape[3]):
    c[:,:,:a.shape[2],:]=a#[:,:,:,:]
#            for k in range(b.shape[2]):
#                for l in range(a.shape[3]):
    c[:,:,a.shape[2]:,:]=b#[:,:,:,:]
                    
    return 

#@njit(cache=True)
#pythran export numba_add(float[], float[])
def numba_add(x,y):
    z=x*y
    return z

@njit(cache=True)
def nmatmul(A,x,b):

    s=A.shape
    for i in range(s[1]):
        b[i]=0.
        for j in range(s[0]):
            b[i]+=A[i,j]*x[j]

@njit(cache=True)
def mmatmul(A,x,b):

    s=A.shape
    for i in range(s[1]):
        b[i]=0.
    for j in range(s[0]):
        for i in range(s[1]):
            b[i]+=A[i,j]*x[j]

@njit(cache=True)
def stationaverage(currentdata,currentdataav,avcount,minindex,maxindex,thresh):

    s=currentdata.shape
    for istat in range(s[0]):
        for ipar in range(s[1]):
            for ip in range(s[2]):
                for it in range(minindex[istat],maxindex[istat]+1):
                    if currentdata[istat,ipar,ip,it]==currentdata[istat,ipar,ip,it]:
                        currentdataav[ipar,ip,it]+=currentdata[istat,ipar,ip,it]
                        avcount[ipar,ip,it]+=1

    for ipar in range(s[1]):
        for ip in range(s[2]):
            for it in range(s[3]):
                if avcount[ipar,ip,it]>=thresh:
                    currentdataav[ipar,ip,it]/=avcount[ipar,ip,it]
                else:
                    currentdataav[ipar,ip,it]=miss_val


@njit
def nodatafound(feld,testindex=-1):
    notfound=True
    if abs(testindex)>=feld.shape[3]:
        testindex=-feld.shape[3]+1
    for i in range(feld.shape[0]):
        for j in range(feld.shape[2]):
            for k in range(feld.shape[4]):
                if feld[i,0,j,testindex,k]==feld[i,0,j,testindex,k]:
                    notfound=False
                    return notfound
    return notfound

def mktmppath(temppath):
    if not os.path.exists(temppath):
        tmppath=temppath.split('/')
        tp='/'
        for tmp in tmppath:
            tp+=tmp+'/'
            if not os.path.exists(tp):
                os.mkdir(tp)

def daysbetween(startdate,filedate):
    s=filedate.split()
    fileday=list(map(int,s[2].split('-')))
    refday=(startdate//10000,(startdate%10000)//100,startdate%100)
    diffd=datetime.date(*fileday)-datetime.date(*refday)

    return diffd.days

def toyearmonth(itime,units):

    if 'days since' in units:
        d=units.split()
        start=list(map(int,d[2].split('-')))
        otime=itime
        for i in range(itime.shape[0]):
            hilf=datetime.datetime(*start)+datetime.timedelta(days=int(itime[i]))
            otime[i]=hilf.year*100+hilf.month
    elif 'hours since' in units:
        d=units.split()
        start=list(map(int,d[2].split('-')))
        otime=itime
        for i in range(itime.shape[0]):
            hilf=datetime.datetime(*start)+datetime.timedelta(hours=int(itime[i]))
            otime[i]=hilf.year*100+hilf.month
    elif units=='yyyymm' or units=='YYYYMM':
        otime=itime
    elif units=='months-since-19790101':
        otime=itime
        for i in range(itime.shape[0]):
            otime[i]=197900+((itime[i]-1)/12)*100+np.mod(itime[i]-1,12)+1
    else:
        print('no valid time unit found, good luck')
        otime=itime
    return otime

def calcdays(startdate,mmax):

    days=np.zeros(mmax,dtype=np.int64)
    syear=startdate//10000
    sdate=date(syear,1,1)
    for k in range(mmax):
        delta=date(syear+k//12,k%12+1,1)-sdate
        days[k]=delta.days+1
    return days

def stats(feld,weights=None,dim=None,lang='ger',mima=None,short=None):

    if weights is not None:
        if dim is not None:
            pass
        else:
            dim=0
    else:
        weights=np.asarray([0])
        if dim is not None:
            pass
        else:
            dim=0

    arr=feld.ravel()
    
#    t=time.time()
    s,sq=statcore(feld,arr,weights,dim)
#    print 'statcore: ',time.time()-t
    
    sig=math.sqrt(np.abs(sq-s*s))
    rms=math.sqrt(np.abs(sq))
    if short is not None:
        if mima is not None and feld.size>0:
            smima='/%4.1f/%4.1f'%(np.nanmin(feld.flatten()),np.nanmax(feld.flatten()))
        else:
            smima=''
    
        if rms>10:
            if abs(s)<2.0:
                return "%3.1f/%3.0f/%3.0f" %( s,rms,sig)+smima
            else:
                return "%3.0f/%3.0f/%3.0f" %( s,rms,sig)+smima
                
        else:
            return "%4.1f/%4.1f/%4.1f" %( s,rms,sig)+smima
    else:
    
        if mima is not None:
            if sig>10:
                smima=', Min: %3.0f, Max: %3.0f'%(np.nanmin(feld.flatten()),np.nanmax(feld.flatten()))
            else:
                smima=', Min: %4.1f, Max: %4.1f'%(np.nanmin(feld.flatten()),np.nanmax(feld.flatten()))
                
        else:
            smima=''
    
        if rms>10:
            fs='Mittel: %3.0f, RMS: %3.0f, Sig: %3.0f'
        else:
            fs='Mittel: %4.1f, RMS: %4.1f, Sig: %4.1f'
            
        if lang=='ger':
            return fs %( s,rms,sig)+smima
        else:
            return "Mean: %4.1f, RMS: %4.1f, Sig: %4.1f" %( s,rms,sig)+smima

@njit(cache=True)
def statcore(feld,arr,weights,dim):
    s=0.
    sq=0.
    count=0.

    if weights.shape[0]!=1:

        stride=1
        for i in range(feld.ndim-1,dim,-1):
            stride*=feld.shape[i]

        k=0
        ki=0
        for i in range(arr.shape[0]):
            if arr[i]==arr[i]:
                s+=arr[i]*weights[k]
                sq+=arr[i]*arr[i]*weights[k]
                count+=weights[k]
            ki+=1
            if ki==stride:
                k+=1
                ki=0
                if k==weights.shape[0]:
                    k=0

    else:
        for i in range(arr.shape[0]):
            if arr[i]==arr[i]:
                s+=arr[i]
                sq+=arr[i]*arr[i]
                count+=1.0

    if count>0:             
        s/=count
        sq/=count

    return s,sq



@jit(cache=False,nopython=True)
def had_rebin_3672_to_1836(mem,hadens,ens,startyear,endyear,nc_miss_val,miss_val):
    for itime in prange((startyear-1850)*12,(endyear-1850+1)*12):
        phad_rebin_3672_to_1836(mem,hadens,ens,startyear,endyear,nc_miss_val,miss_val,itime)
        
@jit(cache=True,nopython=True)
def phad_rebin_3672_to_1836(mem,hadens,ens,startyear,endyear,nc_miss_val,miss_val,itime):

    index=hadens.shape[0]-hadens.shape[0]+itime-(startyear-1850)*12
#    for itime in range((startyear-1850)*12,(endyear-1850+1)*12):
    if itime>=mem.shape[0]:
        return
    for ilat in range(hadens.shape[2]):
        for ilon in range(hadens.shape[3]):
            if ilon<18:
                ishift=18
            else:
                ishift=-18
            sum=nc_miss_val-nc_miss_val
            n=nc_miss_val-nc_miss_val
            jj=(ilon+ishift)*2
            for i in range(2):
                ii=ilat*2+i
                for j in range(2):
                    x=mem[itime,ii,jj+j]
                    if(x>0.9*nc_miss_val):
                        sum+=x
                        n+=1.0
            if n>0:
                hadens[ens,index,ilat,ilon]=sum/n
            else:
                hadens[ens,index,ilat,ilon]=miss_val

    return
@njit(cache=True)
def rebin_72144_to_1836(mem,hadens,ens,global_startyear,startyear,endyear):

    for ip in range(mem.shape[0]):
        index=hadens.shape[0]-hadens.shape[0]-1
        for itime in range((startyear-global_startyear)*12,(endyear-global_startyear+1)*12):
            index+=1
            for ilat in range(hadens.shape[3]):
                for ilon in range(hadens.shape[4]):
                    sum=np.float32(0.)
                    n=np.float32(0.)
                    jj=ilon*4
                    for i in range(4):
                        ii=ilat*4+i
                        for j in range(4):
                            x=mem[ip,itime,ii,jj+j]
                            if x==x :
                                sum+=x
                                n+=1.0
                    if n>0:
                        hadens[ens,ip,index,ilat,ilon]=sum/n
                    else:
                        hadens[ens,ip,index,ilat,ilon]=np.nan

    return

@njit
def rebin_361576_to_1836(mem,hadens,ens,l,midx,global_startyear,startyear,endyear):

    for ip in range(hadens.shape[1]):
        index=hadens.shape[0]-hadens.shape[0]-1
        for itime in range(12): #(startyear-global_startyear)*12,(endyear-global_startyear+1)*12):
#            print ip,itime
            index+=1
            for ilat in range(hadens.shape[3]):
                for ilon in range(hadens.shape[4]):
                    sum=np.float32(0.)
                    n=np.float32(0.)
                    jj=ilon*16
                    for i in range(0,20):
                        ii=ilat*20+i
                        for j in range(0,16):
                            x=mem[itime,midx[ip],ii,jj+j]
                            if x==x :
                                sum+=x
                                n+=1.0
                    if n>0:
                        hadens[ens,ip,l*12+itime,ilat,ilon]=sum/n
                    else:
                        hadens[ens,ip,l*12+itime,ilat,ilon]=np.nan

    return
@njit(cache=True)
def rebin_256512_to_1836(mem,hadens72,hadens,ens,l,midx):

    #for ip in range(hadens.shape[1]):
#        index=hadens.shape[0]-hadens.shape[0]-1
#        for itime in range(12): #(startyear-global_startyear)*12,(endyear-global_startyear+1)*12):
#            print ip,itime
#            index+=1
    for h in hadens72,hadens:
        for ilat in range(h.shape[3]):
            for ilon in range(h.shape[4]):
    
                jstep=mem.shape[0]/h.shape[3]
                jstart=int(ilat*jstep)
                jstop=int((ilat+1)*jstep)
                istep=mem.shape[1]/h.shape[4]
                istart=int(ilon*istep)
                istop=int((ilon+1)*istep)
                for ip in range(h.shape[1]):
                    x=0.
                    n=0
                    for jj in range(jstart,jstop):     
                        for ii in range(istart,istop):
                            for it in range(mem.shape[2]):
                                y=mem[jj,ii,it,midx[ip]]
                                if y==y:
                                    x+=y
                                    n+=1
                    if n>0:
                        h[ens,ip,l,ilat,ilon]=x/n
                    else:
                        h[ens,ip,l,ilat,ilon]=np.nan
    print(hadens72[ens,5,l,0,0],hadens[ens,5,l,0,0])

    return

@njit
def rebin_361576_to_72144(mem,hadens,ens,l,midx,global_startyear,startyear,endyear):

    for ip in range(hadens.shape[1]):
        index=hadens.shape[0]-hadens.shape[0]-1
        for itime in range(12): #(startyear-global_startyear)*12,(endyear-global_startyear+1)*12):
#            print ip,itime
            index+=1
            for ilat in range(hadens.shape[3]):
                for ilon in range(hadens.shape[4]):
                    sum=np.float32(0.)
                    n=np.float32(0.)
                    jj=ilon*4
                    for i in range(0,5):
                        ii=ilat*5+i
                        for j in range(0,4):
                            x=mem[itime,midx[ip],ii,jj+j]
                            if x==x :
                                sum+=x
                                n+=1.0
                    if n>0:
                        hadens[ens,ip,l*12+itime,ilat,ilon]=sum/n
                    else:
                        hadens[ens,ip,l*12+itime,ilat,ilon]=np.nan

    return

@njit(cache=True)
def getindex (alldates,somedates,index):
#   import np
#   from scipy.stats.stats import nanmean
    n=somedates.shape[0]
    m=alldates.shape[0]
    index[0]=0
    jsave=0
    for j in range(n):
        for k in range(index[jsave],m):
            if alldates[k]==somedates[j] or alldates[k]==somedates[j]-13:
                index[j]=k
                jsave=j
                break

    return

@njit(cache=True)
def find_gstatindex(glons,glats,lons,lats,gstatindex):

    for l in range(lats.shape[0]):
        if np.isnan(lons[l]+lats[l]):
            continue
            ilon = glons.shape[0] // 2
            ilat = glats.shape[0] // 2
        else:
            ilon=int(np.floor(lons[l]/(360./glons.shape[0])))
            ilat=int(np.floor((lats[l]+90.)/(180./glats.shape[0])))
        gstatindex[ilat,ilon,0]+=1
        if gstatindex[ilat,ilon,0] < gstatindex.shape[2]:           
            gstatindex[ilat,ilon,gstatindex[ilat,ilon,0]]=l
        else:
            raise IndexError
    return

#@njit(cache=True)
def getindex2old (alldates,somedates,index):
#   import np
#   from scipy.stats.stats import nanmean
    n=somedates.shape[0]
    m=alldates.shape[0]
    index[0]=0
    jsave=0
    for j in range(n):
        for k in range(index[jsave],m):
            if alldates[k]==somedates[j]:
                index[j]=k
                jsave=j
                break
            if alldates[k]>somedates[j]:
                index[j]=k
                if j+1<n:
                    jsave=j+1
                else:
                    jsave=j
                break
 #   if jsave<n:
    if index[jsave] == 0:
        index[jsave]=m-1

    return
@njit(cache=True)
def getindex2 (alldates,somedates,index):
#   import np
#   from scipy.stats.stats import nanmean
    n=somedates.shape[0]
    m=alldates.shape[0]
    index[0]=0
    jsave=0
    for j in range(n):
        index[j]=int(somedates[j]/30.5)
        #for k in range(index[jsave],m):
            #if alldates[k]==somedates[j]:
                #index[j]=k
                #jsave=j
                #break
            #if alldates[k]>somedates[j]:
                #index[j]=k
                #if j+1<n:
                    #jsave=j+1
                #else:
                    #jsave=j
                #break
 #   if jsave<n:
 #   if index[jsave] == 0:
 #       index[jsave]=m-1

    return


@njit(cache=True)
def copystride4(a,b,index,n,m,pindex,fmiss_val):
#   expand time dimension, reduce p dimension

        for l in range(b.shape[0]):              
            for k in range(pindex.shape[0]): # reduce p dimension 
                pik=int(pindex[k])
                for j in range(index.shape[0]):
                    if b[l,pindex[k],j] != fmiss_val:
                        a[m,l,k,index[j]]=b[l,pindex[k],j]
        return
        
@njit(cache=True)
def copystride(a,b,index,n,m,pindex,fmiss_val):
#   expand time dimension, reduce p dimension

     for l in range(b.shape[0]):              
        for k in range(pindex.shape[0]): # reduce p dimension 
            pik=int(pindex[k])
            for j in range(index.shape[0]):
                if b[l,pik,j] != fmiss_val:
                    a[n,m,l,k,index[j]]=b[l,pik,j]
 
     return

@njit(parallel=True)
def picopy(a,b,c,iens,opindex):
#   expand time dimension, reduce p dimension
#   currentdata[:]=currentdatatm+tasks[key]["data"][:,iens,:,:,:].reshape(tmshape[0],tmshape[2],ps.shape[0],tmshape[4])[:,:,pindex,:]

    pindex=np.zeros_like(opindex)
    if opindex[0]>c.shape[3]-1:  
        pindex[0]=0
    else:
        pindex[:]=opindex[:]
    for m in range(c.shape[0]):              
        for l in range(c.shape[2]):              
            for k in range(pindex.shape[0]): # reduce p dimension 
#                pik=pindex[k]
                for i in range(c.shape[4]):
                    a[m,l,k,i]=b[m,l,k,i]+c[m,iens,l,pindex[k],i]

    return

@njit(cache=True)
def expand(b,index,pindex,a):
#   import np
#   from scipy.stats.stats import nanmean
#   copy strides from up to 8 dimensional arrays

    for l in range(b.shape[0]):              
        for k in range(b.shape[2]):
            for m in range(pindex.shape[0]):
                for j in range(b.shape[4]-1):
                    if index[l,j+1]==0:
                        break
                    for i in range(index[l,j],index[l,j+1]):
                        a[l,0,k,m,i]=b[l,0,k,pindex[m],j]

    return

@njit(cache=True)
def expand2(b,index,pindex,a):
#   import np
#   from scipy.stats.stats import nanmean
#   copy strides from up to 8 dimensional arrays

    for k in range(b.shape[0]):
        for m in range(pindex.shape[0]):
            for j in range(b.shape[2]-1):
                if index[j+1]==0:
                    break
                for i in range(index[j],index[j+1]):
                    a[k,m,i]=b[k,pindex[m],j]

    return


@njit(parallel=True)
def expandandadd(b,ref,index,opindex,iens,a,isign):
    
    sign=np.float32(isign)
    if b.shape[0]>a.shape[0]:
        print('expandandadd',b.shape[0],'>',a.shape[0],'exiting')
    for l in prange(b.shape[0]):  
        pexpandandadd(b,ref,index,opindex,iens,a,sign,l)

@njit(cache=True)
def pexpandandadd(b,ref,index,opindex,iens,a,sign,l):
#   import np
#   from scipy.stats.stats import nanmean
#   copy strides from up to 8 dimensional arrays
    pindex=np.zeros(opindex.shape[0],dtype=np.int64)
    pindex[:]=opindex[:]
    pm=int(pindex[0])
    if b.shape[3]==1:
        pindex[0]=0
    if abs(sign)<0.1:
#        for l in range(b.shape[0]):  
        for k in range(b.shape[2]):
            for m in range(pindex.shape[0]):
                ipm=pindex[m]
                for j in range(b.shape[4]):
                    a[l,k,m,j]=b[l,iens,k,ipm,j]
    else:
#        for l in range(b.shape[0]):  
        for k in range(b.shape[2]):
            for m in range(pindex.shape[0]):

                if ref.shape[2]!=a.shape[2]:
                    pm=pindex[m]#+m-m
                else:
                    pm=m#pindex[m]+m-pindex[m]
                for j in range(b.shape[4]-1):
                    if index[l,iens,j+1]==0:
                        if j==0:
                            index[l,iens,j+1]=a.shape[3]
                            for i in range(index[l,iens,j+1]):
                                a[l,k,m,i]=ref[l,k,pm,i]
                            index[l,iens,j+1]=0   
                        else:  
                            index[l,iens,j+1]=a.shape[3]
                            for i in range(index[l,iens,j],index[l,iens,j+1]):
                                a[l,k,m,i]=ref[l,k,pm,i]
                            index[l,iens,j+1]=0   
                        break
                    else:
                        ipm=pindex[m]
                        if sign==1.0:
                            for i in range(index[l,iens,j],min(index[l,iens,j+1],ref.shape[3])):
                                    a[l,k,m,i]=ref[l,k,pm,i]+b[l,iens,k,ipm,j]
                        else:
                            for i in range(index[l,iens,j],min(index[l,iens,j+1],ref.shape[3])):
                                    a[l,k,m,i]=ref[l,k,pm,i]-b[l,iens,k,ipm,j]

    return

@njit(cache=True)
def tdist(dists,lats,lons,weight):
#   import np
#   from scipy.stats.stats import nanmean
#   copy strides from up to 8 dimensional arrays

    if dists[0]==0: #nan?
        x=np.cos(lats*math.pi/180.)*np.cos(lons*math.pi/180.)
        y=np.cos(lats*math.pi/180.)*np.sin(lons*math.pi/180.)
        z=np.sin(lats*math.pi/180.)

        sdistp(dists, x, y, z)
        dists[:]=np.arccos(dists*(1.0-10*np.finfo(dists.dtype).eps))
#        print lats,lons,dists
        if weight!=0:
            dists[:]=np.exp(-dists*40./2/math.pi)*np.sqrt(800.)

    return x,y,z

@njit(cache=True)# ,fastmath=True)
def sdistp(dists,x,y,z):

    idx=np.zeros(x.shape[0],dtype=np.int32)
    for l in range(x.shape[0]-1):
        idx[l+1]=idx[l]+x.shape[0]-l
    for l in range(x.shape[0]):
        for k in range(l,x.shape[0]):
            dists[idx[l]+k-l]=x[l]*x[k]+y[l]*y[k]+z[l]*z[k]

    return

#sid is index of station to which distances are extracted
#dists are precalculated distances
@njit(cache=True)
def extract(dists,sid,lats=None,lons=None, dlen=None):

    if dlen is not None:
        n = np.int32(np.sqrt(dlen) *2)
    else:
        
        n=np.int32(np.sqrt(dists.shape[0]*2))
    idx=np.zeros(n,dtype=np.int32)
    stdists=np.zeros(n)
    for l in range(n-1):
        idx[l+1]=idx[l]+n-l
    m=0
    for l in range(sid):
        stdists[m]=dists[idx[l]+sid-l]
#        print('l {:6.2f},{:6.2f},{:6.2f},{:6.2f},{:d}'.format(lats[l],lons[l],lats[sid],lons[sid],np.int(stdists[m]*6370.)))
        m+=1
    for k in range(sid,n):
        stdists[m]=dists[idx[sid]+k-sid]
#        print('k {:6.2f},{:6.2f},{:6.2f},{:6.2f},{:d}'.format(lats[sid],lons[sid],lats[k],lons[k],np.int(stdists[m]*6370.)))
        m+=1

    return stdists

@njit(cache=True)
def sdist(dists,x,y,z):

    id=0
    for l in range(x.shape[0]):
        for k in range(l,x.shape[0]):
            dists[id]=x[l]*x[k]+y[l]*y[k]+z[l]*z[k]
            id+=1

    return
#@njit(cache=True)
def fdist(slat,slon,lat,lon):

#    if typeof(slon) is numba.types.misc.UnicodeType:
    #try:
        if type(slon) is str:
            if '.' in slon:
                flon=float(slon)
                flat=float(slat)
            elif ',' in slon:
                flon=float('.'.join(slon.split(',')))
                flat=float('.'.join(slat.split(',')))
            elif slon[-1] in 'WE':
                flon=float(slon[:3])+float(slon[4:7])/60+float(slon[-2:-1])/3600
                if slon[-1]=='W':
                    flon=-flon
                flat=float(slat[:2])+float(slat[3:6])/60+float(slat[-2:-1])/3600
                if slat[-1]=='S':
                    flat=-flat
            else:
                try:
                    
                    flon=float(slon[-9:-6])+float(slon[-5:-3])/60+float(slon[-2:])/3600
                    if slon[0]=='-':
                        flon=-flon
                    flat=float(slat[-8:-6])+float(slat[-5:-3])/60+float(slat[-2:])/3600
                    if slat[0]=='-':
                        flat=-flat
                except:
                    flon=np.nan
                    flat=np.nan
                
        else:
            flon=slon
            flat=slat
        
        
        lats=np.append(flat,lat)
        lons=np.append(flon,lon)
        
        x=np.cos(lats*math.pi/180.)*np.cos(lons*math.pi/180.)
        y=np.cos(lats*math.pi/180.)*np.sin(lons*math.pi/180.)
        z=np.sin(lats*math.pi/180.)
    
        dist=np.arccos((x[0]*x[1:]+y[0]*y[1:]+z[0]*z[1:])*0.9999999)
            
    #except Exception as e:
        #exc_type, exc_obj, exc_tb = sys.exc_info()
        #fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        #print(exc_type, fname, exc_tb.tb_lineno)
    
        return dist

@njit
def pcost(id,cost,slopes,dists,l):
    
        for k in range(slopes.shape[0]):
            if slopes[l]==slopes[l] and slopes[k]==slopes[k]:
                s=(slopes[l]-slopes[k])*dists[id[l]+k]
                cost[l]+=s*s
#                cost[k]+=s*s
#                    if l==24 and k==29:
#                        print dists[id],slopes[l],slopes[k]
#            id+=1      
        if goodstats>0:
            tcost+=cost[l]
    #else:
        #id+=slopes.shape[0]-l 
@njit        
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
                    if dists[id]!=dists[id]:
                        print('spurious',l,k,id)
                    else:
                        s = (slopes[l]-slopes[k])*dists[id]
                        cost[l] += s*s
                        cost[k] += s*s
                    #print((dists[id], slopes[l], slopes[k]))
                id += 1
            if goodstats > 0:
                tcost += cost[l]
        else:
            id += slopes.shape[0]-l

    if goodstats > 0:
        tcost /= goodstats
        for l in range(slopes.shape[0]):
            cost[l] /= goodstats

    print('tcost:',tcost)
    return tcost

@njit(cache=True,parallel=True)
def tcost2(dists,slopes,cost):
#   calculate trend heterogeneity cost function

    id=0
    goodstats=0
    tcost=0.
    cost[:]=0
    costk=np.zeros_like(cost)
    id=np.zeros(slopes.shape[0],np.int32)
    s=np.zeros_like(dists)
    for l in range(1,slopes.shape[0]):
        id[l]=id[l-1]+slopes.shape[0]-l
        if slopes[l]==slopes[l]:
            goodstats+=1
            
    for l in prange(slopes.shape[0]):
        if slopes[l]==slopes[l]:
            for k in range(l,slopes.shape[0]):
                if slopes[k]==slopes[k]:     
                    m=id[l]+k-l
                    s[m]=(slopes[l]-slopes[k])*dists[m]
                    s[m]*=s[m]
                    costk[k]+=s[m]
            cost[l]=np.nansum(s[id[l]:m+1])
            
    for l in range(slopes.shape[0]):
        cost[l]+=costk[l]
    tcost=np.nansum(cost)
    if goodstats>0:
        tcost/=goodstats
        for l in range(slopes.shape[0]):
            cost[l]/=goodstats


    return tcost

@njit(cache=True)
def thin(t,index,n):

    ni=t.shape[0]
    if n<2:
        for i in range(t.shape[0]):
            index[i]=i
    else:
        ni=t.shape[0]//n
        for i in range(ni):
            for j in range(n):
                index[i]=i*n
                if t[i*n+j]==t[i*n+j]:
                    index[i]=i*n+j
                    break

    return ni

@njit(cache=True)
def thin2(t,n):

    ni=t.shape[0]
    index=np.zeros(ni//n,dtype=int32)-1
    if n<2:
        for i in range(t.shape[0]):
            index[i]=i
    else:
        ni=t.shape[0]//n
        for i in range(ni):
            index[i]=i*n
            for j in range(n):
                idx=i*n+j
                if t[idx]==t[idx]:
                    index[i]=idx
                    break

    return index
def rmeanw(t,runmean):
    tmean=t.copy()
    index=np.zeros(tmean.shape[0],dtype='int')
    
    tret=rmean(t,tmean,index,runmean)
    tret[:runmean]=np.nan
    tret[-runmean:]=np.nan
    return tret

@njit(cache=True)
def rmean(t,tmean,index,runmean):

    tret=np.zeros(t.shape[0])
    ni=t.shape[0]
    good=runmean-runmean
    if runmean<2:
        for i in range(ni):
            tret[i]=t[i]
    else:

        for j in range(ni):
            tret[j]=np.nan
            if t[j]==t[j]:
                index[good]=j
                good+=1

        if good>runmean+2:
            i=runmean//2
            tmean[:]=np.nan
            if runmean%2==1:
                tmean[i]=0.
                for k in range(-runmean//2+1,runmean//2+1):
                    tmean[i]+=t[index[i+k]]
                tmean[i]/=runmean

                for i in range(runmean//2+1,good-runmean//2):
                    tmean[i]=(tmean[i-1]*runmean+t[index[i+runmean//2]])/(runmean+1)

            else:

                i=runmean//2
                tmean[i]=0.
                for k in range(-runmean//2,runmean//2):
                    tmean[i]+=t[index[i+k]]
                tmean[i]/=runmean

                for i in range(runmean//2+1,good-runmean//2-1):
                    tmean[i]=(tmean[i-1]*runmean+t[index[i+runmean//2-1]])/(runmean+1)

            for i in range(good):
                tret[index[i]]=tmean[i]
        else:
            for i in range(good):
                tret[index[i]]=t[index[i]]

    return tret
@njit
def monmean(var,tidx,thresh=15):
    out=np.empty((2,tidx.shape[0]))
    for i in range(tidx.shape[0]-1):
        mm=0.0
        mc=0
        for ipar in range(2):
            for j in range(tidx[i+1]-tidx[i]):
                v=var[ipar,tidx[i]+j]
                if v==v:
                    mm+=v
                    mc+=1
            if mc>thresh:
                out[ipar,i]=mm/mc
            else:
                out[ipar,i]=np.nan
    out[:,-1]=np.nan
    return out
            

@njit(cache=True)
def snhtmov(t,tsa,snhtparas,index,count,tmean,tsquare):

    n=snhtparas[0]
    max_miss=snhtparas[1]
    ninc=snhtparas[2]

    ni=t.shape[0]
    good=0
    for j in range(ni):
        if t[j]==t[j]:
            index[good]=j
            good+=1
        count[j]=good

    if good>n-2*max_miss:
        rm=n/2
        i=rm/2
        for k in range(n/2,ni,ninc):
            if count[k]-count[k-n/2]>n/2-max_miss:
                tmean[k]=0.
                tsquare[k]=0.
                for i in range(count[k-n/2],count[k]):
                    x=t[index[i]]
                    tmean[k]+=x
                    tsquare[k]+=x*x
                tmean[k]/=count[k]-count[k-n/2]
                tsquare[k]/=count[k]-count[k-n/2]
        fak=n/2/ninc # approximate length of SNHT window with fak*ninc
        for k in range(n/2,ni,ninc):
            if count[k]-count[k-n/2]>n/2-max_miss and count[k+fak*ninc]-count[k]>n/2-max_miss:
                m=(tmean[k]+tmean[k+fak*ninc])/2.
                tdiv=tsquare[k]-tmean[k]*tmean[k]
                if tdiv>0:
                    tsa[k]=n/2*((tmean[k]-m)*(tmean[k]-m)+(tmean[k+fak*ninc]-m)*(tmean[k+fak*ninc]-m))/math.sqrt(tdiv)
                else:
                    tsa[k]=0.
                for i in range(ninc):
                    tsa[k+i]=tsa[k]

    return

@njit(cache=True)
def snhtmov2(t,tsa,snhtparas,index,count,tmean,tsquare):

    n=snhtparas[0]
    max_miss=snhtparas[1]
    ninc=snhtparas[2]

    ni=t.shape[0]
    good=0
    #for j in range(ni):
        #tmean[j]=np.nan
        #tsquare[j]=np.nan
        #tsa[j]=np.nan
    tmean[0]=0.
    tsquare[0]=0.
    for j in range(ni):
        count[j]=0
        if t[j]==t[j]:
            index[good]=j
            if good>0:
                tmean[good]=tmean[good-1]+t[j]
                tsquare[good]=tsquare[good-1]+t[j]*t[j]
            else:
                tmean[good]=t[j]
                tsquare[good]=t[j]*t[j]      
            good+=1
        if good>0:
            count[j]=good-1

    if good>n-2*max_miss:
        rm=n//2
        for k in range(rm-max_miss,ni-(rm-max_miss)):
            xm=k-rm
            if xm<0:
                xm=0
            xp=k+rm
            if xp>ni-1:
                xp=ni-1
            if count[k]-count[xm]>rm-max_miss and count[xp]-count[k]>rm-max_miss:
                x=(tmean[count[k]]-tmean[count[xm]])/(count[k]-count[xm])
                y=(tmean[count[xp]]-tmean[count[k]])/(count[xp]-count[k])
                xy=(tmean[count[xp]]-tmean[count[xm]])/(count[xp]-count[xm])

                sig=(tsquare[count[xp]]-tsquare[count[xm]])/(count[xp]-count[xm])
                if(sig>xy*xy):
                    sig=math.sqrt(sig-xy*xy)
                    tsa[k]=((count[k]-count[xm])*(x-xy)*(x-xy)+(count[xp]-count[k])*(y-xy)*(y-xy))/sig
                else:
                    tsa[k]=0.

#                print xm,k,xp,tsa[k],sig,count[k],count[xp]-count[k],count[k]-count[xm]

    return

# compare three intervals to detect bad data episodes
@njit(cache=True)
def snhtmov3(t,tsa,snhtparas,index,count,tmean,tsquare,msum,mtmean,mtsquare,month):

    n=snhtparas[0]
    max_miss=snhtparas[1]
    ninc=snhtparas[2]

    ni=t.shape[0]
    good=0
    tmean[0]=0.
    tsquare[0]=0.
    msum[:,0]=0
    mtsquare[:,0]=0.
    mtmean[:,0]=0.
    #tsa[:]=np.nan
    for j in range(ni):
        tsa[j]=np.nan
        count[j]=0
        if t[j]==t[j]:
            index[good]=j
            if good>0:
                tmean[good]=tmean[good-1]+t[j]
                tsquare[good]=tsquare[good-1]+t[j]*t[j]
                mtmean[month[j],good]=mtmean[month[j],good-1]+t[j]
                mtsquare[month[j],good]=mtsquare[month[j],good-1]+t[j]*t[j]
                clast=j
            else:
                tmean[good]=t[j]
                tsquare[good]=t[j]*t[j]      
                mtmean[month[j],good]=t[j]
                mtsquare[month[j],good]=t[j]*t[j]
                cfirst=j
            msum[month[j],good]+=1
            good+=1
            for k in range(12):
                msum[k,good]=msum[k,good-1]
                mtmean[k,good]=mtmean[k,good-1]
                mtsquare[k,good]=mtsquare[k,good-1]
        if good>0:
            count[j]=good-1

    if good>n-3*max_miss:
        rm=n/3
        for k in range(cfirst+2*max_miss,clast-(rm-max_miss)):
            xm2=k-2*rm
            if xm2<0:
                xm2=0
            xm=k-rm
            if xm<0:
                xm=0
            xp=k+rm
            if xp>ni-1:
                xp=ni-1
            if count[xm]-count[xm2]>rm-max_miss and count[k]-count[xm]>rm-max_miss and count[xp]-count[k]>rm-max_miss:
                x=0.
                y=0.
                xy=0.
                l=0
                sig=0.
                for m in range(12):
                    ix=msum[m,count[k]]-msum[m,count[xm]]
                    iy=msum[m,count[xp]]-msum[m,count[k]]+msum[m,count[xm]]-msum[m,count[xm2]]
                    if ix>5 and iy>5:
                        x+=(mtmean[m,count[k]]-mtmean[m,count[xm]])/ix
                        y+=(mtmean[m,count[xp]]-mtmean[m,count[k]]+mtmean[m,count[xm]]-mtmean[m,count[xm2]])/iy
                        xy+=(mtmean[m,count[xp]]-mtmean[m,count[xm2]])/(ix+iy)
                        sig+=(mtsquare[m,count[xp]]-mtsquare[m,count[xm2]])/(ix+iy)
                        l+=1
                if l>3:
                    x/=l
                    y/=l
                    xy/=l
                    sig/=l
                #x=(tmean[count[k]]-tmean[count[xm]])/(count[k]-count[xm])
                #y=(tmean[count[xp]]-tmean[count[k]])+(tmean[count[xm]]-tmean[count[xm2]])
                #y/=(count[xp]-count[k]+count[xm]-count[xm2])
                #xy=(tmean[count[xp]]-tmean[count[xm2]])/(count[xp]-count[xm2])

                #sig=(tsquare[count[xp]]-tsquare[count[xm2]])/(count[xp]-count[xm2])
                if(sig>xy*xy and l>3):
                    sig=math.sqrt(sig-xy*xy)
                    tsa[k]=((count[k]-count[xm])*(x-xy)*(x-xy)+(count[xp]-count[k]+count[xm]-count[xm2])*(y-xy)*(y-xy))/sig
                    #if tsa[k]>1000:
                        #print k
                else:
                    tsa[k]=np.nan

#                print xm,k,xp,tsa[k],sig,count[k],count[xp]-count[k],count[k]-count[xm]

    return

@njit(cache=True)
def snhtmov4(t,tsa,snhtparas,index,count,tmean,tsquare,l,m,o):

    n=snhtparas[0]
    max_miss=snhtparas[1]
    ninc=snhtparas[2]

    ni=t.shape[3]
    good=0
    #for j in range(ni):
        #tmean[j]=np.nan
        #tsquare[j]=np.nan
        #tsa[j]=np.nan
    tmean[0]=0.
    tsquare[0]=0.
    for j in range(ni):
        count[j]=0
        if t[l,m,o,j]==t[l,m,o,j]:
            index[good]=j
            if good>0:
                tmean[good]=tmean[good-1]+t[l,m,o,j]
                tsquare[good]=tsquare[good-1]+t[l,m,o,j]*t[l,m,o,j]
            else:
                tmean[good]=t[l,m,o,j]
                tsquare[good]=t[l,m,o,j]*t[l,m,o,j]      
            good+=1
        if good>0:
            count[j]=good-1

    if good>n-2*max_miss:
        rm=n//2
        for k in range(rm-max_miss,ni-(rm-max_miss)):
            xm=k-rm
            if xm<0:
                xm=0
            xp=k+rm
            if xp>ni-1:
                xp=ni-1
            if count[k]-count[xm]>rm-max_miss and count[xp]-count[k]>rm-max_miss:
                x=(tmean[count[k]]-tmean[count[xm]])/(count[k]-count[xm])
                y=(tmean[count[xp]]-tmean[count[k]])/(count[xp]-count[k])
                xy=(tmean[count[xp]]-tmean[count[xm]])/(count[xp]-count[xm])

                sig=(tsquare[count[xp]]-tsquare[count[xm]])/(count[xp]-count[xm])
                if(sig>xy*xy):
                    sig=math.sqrt(sig-xy*xy)
                    tsa[l,m,o,k]=((count[k]-count[xm])*(x-xy)*(x-xy)+(count[xp]-count[k])*(y-xy)*(y-xy))/sig
                #else:
                    #tsa[l,m,o,k]=0.

#                print xm,k,xp,tsa[k],sig,count[k],count[xp]-count[k],count[k]-count[xm]

    return


@njit(cache=True)
def zonaltrends(gslopes,zslopes):

    s=gslopes.shape
#    zslopes=np.zeros([s[0],s[3]])
    for k in range(s[0]):
        for ip in range(s[3]):
            zslopes[k,ip]=0.
            b=0
            for ipar in range(s[2]):
                for j in range(s[1]):
                    if gslopes[k,j,ipar,ip]==gslopes[k,j,ipar,ip]:
                        zslopes[k,ip]+=gslopes[k,j,ipar,ip]
                        b+=1
#                mask=~np.isnan(gslopes[k,:,ipar,ip])
#                b+=sum(mask)
#                zslopes[k,ip]+=sum(gslopes[k,mask,ipar,ip])
            if(b>0):
                zslopes[k,ip]/=b
            else:
                zslopes[k,ip]=np.nan
    return zslopes

def belttrends(zslopes,belts):
    beltslopes=np.zeros([belts.shape[0],zslopes.shape[1]])
    weights=np.empty(zslopes.shape[0])
    for k in range(zslopes.shape[0]):
        weights[k]=np.cos((-85.+k*10)*np.pi/180.)
    for k in range(belts.shape[0]):
        wk=weights[belts[k,0]:belts[k,1]]
        for ip in range(zslopes.shape[1]):
            hilf=zslopes[belts[k,0]:belts[k,1],ip]
            mask=~np.isnan(hilf)
            b=sum(mask)
            if(b>0):
                beltslopes[k,ip]=sum(hilf[mask]*wk[mask])/np.sum(wk[mask])
            else:
                beltslopes[k,ip]=np.nan
    return beltslopes

@njit(cache=True)
def hmav13(x,y):
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
                y[i]=np.nan
        else:
            y[i]=np.nan 

    last=x.shape[0]
    for i in range(last-1,6,-1):
        if x[i]!=x[i]:
            y[i]=np.nan
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
                y[i]=np.nan
        else:
            y[i]=np.nan 

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
                y[i]=np.nan
        else:
            y[i]=np.nan 

    return y

def rgb(r,g,b):
    return np.asarray([r,g,b],dtype=np.float32)

