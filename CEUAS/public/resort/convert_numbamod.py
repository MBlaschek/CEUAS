#!/usr/bin/env
# coding: utf-8

from numba import njit, prange
import numpy as np
import sys,glob
import zipfile, os, time
import urllib3
from datetime import datetime, timedelta
import glob
sys.path.append(os.getcwd()+'/../cds-backend/code/')
sys.path.append(os.getcwd()+'/../harvest/code/')
sys.path.insert(0,os.getcwd()+'/../resort/rasotools-master/')
#import rasotools
#dir(rasotools)
import warnings
from functools import partial
import matplotlib.pylab as plt
warnings.filterwarnings('ignore')
import pickle
import ray
import rs_drift as rsd


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

@njit(cache=True)
def thin2(t,n):

    ni=t.shape[0]
    index=np.zeros(ni//n,dtype=np.int32)-1
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

# check if surface obs in unused IGRA

@njit
def is_surf3(oh, z, ori, mask, maskm1, imatch, iri, zi, ohi):
    for j in range(len(ori)):
        if mask[j]:
            jj = j
        elif maskm1[j]:
            jj = j - 1
        if mask[j] | maskm1[j] :
            if imatch[jj] == iri.shape[0] - 1:
                istop = iri.shape[0]
            else:
                istop = iri[imatch[jj]+1]
            for k in range(iri[imatch[jj]], istop):
                if ohi[k] == 1:
                    if j == len(ori) - 1:
                        jstop = z.shape[0]
                    else:
                        jstop = ori[j+1]
                    zmatch = np.searchsorted(z[ori[j]:jstop], zi[k])
                    if z[ori[j]+zmatch] == zi[k]:
                        l = 0
                        while ori[j]+zmatch+ l < z.shape[0] and z[ori[j]+zmatch+ l]  == z[ori[j]+zmatch]:
                            oh[ori[j]+zmatch+ l] = 1.0
                            l += 1
                    break
#zdiff = is_surf2(out['z_coordinate'], out['observed_variable'], out['observation_value'], evarno, oh, ori, station_elevation)
@njit
def is_surf2(z, ov, ovar, evarno, oh, ori, ps, station_height):
    
    zdiff = np.full_like(z, np.nan, dtype=np.float32)
    g = np.float32(9.80665)
    sv = np.array((40, 39, 58, 281, 41, 42))
    for i in range(len(ori)-1):
        sl = slice(ori[i], ori[i+1])
        ix = -1
        for j in range(sl.start, sl.stop):
            if np.any(evarno[j]==sv):
                ix = j
                break
            if (z[j] == z[j]) and (z[j]>z[ix]) and (z[j] != 100000.) and (z[j] != 85000.) and np.abs(z[j]-93215.53) > 0.1 and (np.abs(p_to_z_ifs(z[j]/100.) / g - station_height) < 500):
                #print(z[j], (np.abs(p_to_z_ifs(z[j]/100.) / g - station_height) < 200), np.abs(p_to_z_ifs(z[j]/100.) / g - station_height))
                ix = j
                
        if ix==-1:
                continue
                
        j = 0
        k = -1
        while ix + j < sl.stop and z[ix+j] == z[ix]:
            zdiff[ix+j] = 0.
            if ov[ix+j] == 117:
                zdiff[ix+j] = ovar[ix+j] / g - station_height
                k = j
            j += 1
        if k != -1:
            zdiff[ix:ix+j] = zdiff[ix+k]
    return zdiff
                  
        
    
@njit
def is_surf(z, ov, oh, min_press, idx):
    
    tg = tb = wg = wb = s2 = 0
    wind = np.array((106, 107, 139, 140))
    for ii in idx:
        j = 0
        if z[ii-1] < z[ii]:
            ii -= 1
        k = 0
        while (z[ii-j] == z[ii]):
            if z[ii-j] < min_press:
                j += 1
                continue
            if np.any(wind == ov[ii-j]):
                if oh[ii-j] != oh[ii-j]:
                    oh[ii-j] = 10.
                    wb += 1
                else:
                    wg += 1
            else:
                if oh[ii-j] != oh[ii-j]:
                    oh[ii-j] = 2.
                    tb += 1
                else:
                    tg += 1
            j += 1
        while (ii - j >= 0) and (z[ii-j] <= z[ii-j+1]):
            if oh[ii-j] == oh[ii-j]:
                #print(ii, 'two surface values, keeping that from the odb')
                oh[ii-j+1:ii] = np.nan
                break
            j += 1
    return wb, wg, tb, tg         

@njit
def rms(x):
    #x = np.array(x)
    mask = ~np.isnan(x)
    if np.any(mask):
        return np.sqrt(np.mean(x[mask]*x[mask])), np.sum(mask)
    else:
        return np.nan, 0
    
@njit
def make_vrindex(vridx,ridx,idx):
    l=0
    for i in range(1,len(idx)): # to set the recordindices
        if ridx[i]>ridx[i-1]:
            vridx[ridx[i-1]:ridx[i]]=l # next record after l
            l=i
    vridx[ridx[i]:]=len(idx) # next record for the last element is the len of the data


from numba.typed import List



        #if observed_variable[i]==ipar[38]:
            #relhum[j]=observation_value[i]
            #d_relhum[j]=observation_value[i]-departure[i]
            #fgd_relhum[j]=observation_value[i]-fg_departure[i]
        #elif observed_variable[i]==ipar[39]:
            #spechum[j]=observation_value[i]
            #d_spechum[j]=observation_value[i]-departure[i]
            #fgd_spechum[j]=observation_value[i]-fg_departure[i]
        #elif observed_variable[i]==ipar[34]:
            #dpd[j]=observation_value[i]
            #d_dpd[j]=observation_value[i]-departure[i]
            #fgd_dpd[j]=observation_value[i]-fg_departure[i]
        #elif observed_variable[i]==ipar[36]:
            #dewpoint[j]=observation_value[i]
            #d_dewpoint[j]=observation_value[i]-departure[i]
            #fgd_dewpoint[j]=observation_value[i]-fg_departure[i]
        #elif observed_variable[i]==ipar[85]:
            #temp[j]=observation_value[i]
            #d_temp[j]=observation_value[i]-departure[i]
            #fgd_temp[j]=observation_value[i]-fg_departure[i]
        #elif observed_variable[i]==ipar[104]:
            #uwind[j]=observation_value[i]
            #d_uwind[j]=observation_value[i]-departure[i]
            #fgd_uwind[j]=observation_value[i]-fg_departure[i]
        #elif observed_variable[i]==ipar[105]:
            #vwind[j]=observation_value[i]
            #d_vwind[j]=observation_value[i]-departure[i]
            #fgd_vwind[j]=observation_value[i]-fg_departure[i]
        #elif observed_variable[i]==ipar[106]:
            #wd[j]=observation_value[i]
            #d_wd[j]=observation_value[i]-departure[i]
            #fgd_wd[j]=observation_value[i]-fg_departure[i]
        #elif observed_variable[i]==ipar[107]:
            #ws[j]=observation_value[i]
            #d_ws[j]=observation_value[i]-departure[i]
            #fgd_ws[j]=observation_value[i]-fg_departure[i]
        #else:
            #pass


#@njit(boundscheck=False)
def qconvert(j,k,h,a_observation_value,a_conversion_flag,a_conversion_method,
             a_an_depar,a_fg_depar,
             temp,cdpddp,cdpdrh,crhdpd,cshrh,cshdpd,crhsh,cdpdsh,
             d_cdpddp,d_cdpdrh,d_cdpdsh,d_cshrh,d_cshdpd,d_crhsh,d_crhdpd,
             fgd_cdpddp,fgd_cdpdrh,fgd_cdpdsh,fgd_cshrh,fgd_cshdpd,fgd_crhsh,fgd_crhdpd):
    if h==ipar[34]:
        if cdpddp[k]==cdpddp[k]:
            a_observation_value[j]=cdpddp[k]
            a_an_depar[j]=cdpddp[k]-d_cdpddp[k]
            a_fg_depar[j]=cdpddp[k]-fgd_cdpddp[k]
            if (np.abs(cdpddp[k])>80) or (cdpddp[k]<0.01):
                a_observation_value[j]=np.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=2
        elif cdpdrh[k]==cdpdrh[k]:
            a_observation_value[j]=cdpdrh[k]
            a_an_depar[j]=cdpdrh[k]-d_cdpdrh[k]
            a_fg_depar[j]=cdpdrh[k]-fgd_cdpdrh[k]
            if (np.abs(cdpdrh[k])>80) or (cdpdrh[k]<0.01):
                a_observation_value[j]=np.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=3
        else:
            a_observation_value[j]=cdpdsh[k]
            a_an_depar[j]=cdpdsh[k]-d_cdpdsh[k]
            a_fg_depar[j]=cdpdsh[k]-fgd_cdpdsh[k]
            if (np.abs(cdpdsh[k])>80) or (cdpdsh[k]<0.01):
                a_observation_value[j]=np.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=4

    elif h==ipar[36]:
        if cdpdrh[k]==cdpdrh[k]:
            a_observation_value[j]=temp[k]-cdpdrh[k]
            a_an_depar[j]=(temp[k]-cdpdrh[k])-(temp[k]-d_cdpdrh[k])
            a_fg_depar[j]=(temp[k]-cdpdrh[k])-(temp[k]-fgd_cdpdrh[k])
            if (np.abs(cdpdrh[k])>80) or (cdpdrh[k]<0.01):
                a_observation_value[j]=np.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=3
        elif cdpdrh[k]==cdpdrh[k]:
            a_observation_value[j]=temp[k]-cdpddp[k]
            a_an_depar[j]=(temp[k]-cdpddp[k])-(temp[k]-d_cdpddp[k])
            a_fg_depar[j]=(temp[k]-cdpddp[k])-(temp[k]-fgd_cdpddp[k])
            if (np.abs(cdpddp[k])>80) or (cdpddp[k]<0.01):
                a_observation_value[j]=np.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=2
        else:
            a_observation_value[j]=temp[k]-cdpdsh[k]
            a_an_depar[j]=(temp[k]-cdpdsh[k])-(temp[k]-d_cdpdsh[k])
            a_fg_depar[j]=(temp[k]-cdpdsh[k])-(temp[k]-fgd_cdpdsh[k])
            if (np.abs(cdpdsh[k])>80) or (cdpdsh[k]<0.01):
                a_observation_value[j]=np.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=4

    elif h==ipar[38]:
        if crhdpd[k]==crhdpd[k]:
            a_observation_value[j]=crhdpd[k]
            a_an_depar[j]=crhdpd[k]-d_crhdpd[k]
            a_fg_depar[j]=crhdpd[k]-fgd_crhdpd[k]
            if (crhdpd[k]<0.) or (crhdpd[k]>1.03):
                a_observation_value[j]=np.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=2
        else: 
            a_observation_value[j]=crhsh[k]
            a_an_depar[j]=crhsh[k]-d_crhsh[k]
            a_fg_depar[j]=crhsh[k]-fgd_crhsh[k]
            if (crhsh[k]<0.) or (crhsh[k]>1.03):
                a_observation_value[j]=np.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=4

    elif h==ipar[39]:
        if cshrh[k]==cshrh[k]:
            a_observation_value[j]=cshrh[k]
            a_an_depar[j]=cshrh[k]-d_cshrh[k]
            a_fg_depar[j]=cshrh[k]-fgd_cshrh[k]
            if (cshrh[k]<0.) or (cshrh[k]>50.):
                a_observation_value[j]=np.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=3
        else:
            a_observation_value[j]=cshdpd[k]
            a_an_depar[j]=cshdpd[k]-d_cshdpd[k]
            a_fg_depar[j]=cshdpd[k]-fgd_cshdpd[k]
            if (cshdpd[k]<0.) or (cshdpd[k]>50.):
                a_observation_value[j]=np.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=2
    #else:
        #print('qconvert called with wrong variable')

    return

#@njit(boundscheck=True)
def wconvert(j,k,h,a_observation_value,a_conversion_flag,a_conversion_method,
             a_an_depar,a_fg_depar,
             cuwind,cvwind,cwd,cws,
             d_cwd,d_cws,
             fgd_cwd,fgd_cws):
    if h==ipar[104]:
        if cuwind[k]==cuwind[k]:
            a_observation_value[j]=cuwind[k]
            a_an_depar[j]=np.nan
            a_fg_depar[j]=np.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=1
    elif h==ipar[105]:
        if cvwind[k]==cvwind[k]:
            a_observation_value[j]=cvwind[k]
            a_an_depar[j]=np.nan
            a_fg_depar[j]=np.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=1
    elif h==ipar[106]:
        if cwd[k]==cwd[k]:
            a_observation_value[j]=cwd[k]
            a_an_depar[j]=d_cwd[k]
            a_fg_depar[j]=fgd_cwd[k]
            if (cwd[k]<0.) or (cwd[k]>360.):
                a_observation_value[j]=np.nan
            a_conversion_flag[j]=0
            a_conversion_method[j]=2
    elif h==ipar[107]:
        if cws[k]==cws[k]:
            a_observation_value[j]=cws[k]
            a_an_depar[j]=d_cws[k]
            a_fg_depar[j]=fgd_cws[k]
            a_conversion_flag[j]=0
            a_conversion_method[j]=2
    #else:
        #print('wconvert called with wrong variable')

@njit
def do_copy(a_obstab,obstab,j,i):
    # all written here - * will be overwritten if it's a converted variable
    a_obstab['date_time'][j]=obstab['date_time'][i]
    a_obstab['observation_id'][j]=obstab['observation_id'][i] # *
    #a_obstab['report_id'][j]=obstab['report_id'][i] 
    a_obstab['observation_value'][j]=obstab['observation_value'][i] # *
    a_obstab['observed_variable'][j]=obstab['observed_variable'][i] # *
    a_obstab['z_coordinate'][j]=obstab['z_coordinate'][i]
    a_obstab['z_coordinate_type'][j]=obstab['z_coordinate_type'][i]
    a_obstab['conversion_flag'][j]=obstab['conversion_flag'][i] # *
    a_obstab['conversion_method'][j]=obstab['conversion_method'][i] # *
    return

@njit
def do_fb_copy(a_loaded_feedback,loaded_feedback,j,i):
    # all written here - * will be overwritten if it's a converted variable
    a_loaded_feedback['fg_depar@body'][j]=loaded_feedback['fg_depar@body'][i] # *
    a_loaded_feedback['an_depar@body'][j]=loaded_feedback['an_depar@body'][i] # *
    a_loaded_feedback['biascorr@body'][j]=loaded_feedback['biascorr@body'][i]
    a_loaded_feedback['biascorr_fg@body'][j]=loaded_feedback['biascorr_fg@body'][i]
    return

@njit(cache=True)
def z_to_p_ifs(h): # geopotential height (m^/s^2) to pressure (Pa)
    a = 5.252368255329
    b = 44330.769230769
    c = 0.000157583169442
    ptro = 226.547172
    po = 1013.25
    g = 9.80665
 
    h /= g
    if h != h:

        p = h

    elif h > 11000.0:

        y = -c * (h - 11000.0)

        p = ptro * np.exp(y)

    else:

        y = 1.0 - h / b

        p = po * (y**a)

    return p * 100. # we want Pa

@njit(cache=True)
def p_to_z_ifs(p): # pressure (hPa) to geopotential height (m^/s^2) 
    a = 5.252368255329
    b = 44330.769230769
    c = 0.000157583169442
    ptro = 226.547172
    po = 1013.25
    g = 9.80665
 
    if p != p:

        h = p

    elif p < 226.5:
        
        y = np.log(p/ptro)

        h = y /(-c) + 11000.

    else:
        
        y = (p / po) ** (1. / a)
        
        h = (y - 1) * (-b)


    h *= g
    return h

@njit(cache=True)
def extract_tuv(out, lats, lons, obss, ovs, zs, tss, sl, station_elevation):
    
    idx = np.zeros(tss.shape[0], dtype=np.int32)
    k = 0
    pmax = 0
    for i in range(1, tss.shape[0]):
        if tss[i] != tss[i-1]:
            k += 1
            idx[k] = i
            if idx[k] - idx[k-1] > pmax:
                pmax = idx[k] - idx[k-1] 
    idx = idx[:k+1]
    #utss, idx = np.unique(tss, return_index=True)
    #pmax = 0
    #if idx.shape[0] > 1:
        #pmax = np.max(idx[1:]-idx[:-1])
    if tss.shape[0]-idx[-1] > pmax:
        pmax = tss.shape[0]-idx[-1] 
    
    temp = np.zeros(pmax)
    temp[:] = np.nan
    u = np.zeros(pmax)
    v = np.zeros(pmax)
    zu = np.zeros(pmax)
    gph = np.zeros(pmax)
    gph[:] = np.nan
    
    for i in range(len(idx)):
        temp[:] = np.nan
        u[:] = 0.
        v[:] = 0.
        zu[:] = 0.
        gph[:] = np.nan
        
        ii = idx[i]
        if i < len(idx) - 1:
            iip = idx[i+1]
        else:
            iip = tss.shape[0]
        
        k = -1
        zold = -10000
        m = 0
        while zu[0] == 0 and ii + m < iip:
            
            zu[0] = zs[ii+m]
            m += 1
        
        tz = False
        #zu = zs[ii:iip][obss[ii:iip]==139]
        #zu = np.unique(zs[ii:iip])
        #if ii == 27686:
            #print('x')
        for j in range(ii, iip):
            if zs[j] > zold :
                if not tz and k > -1: # neither z nor t are available
                    gph[k] = p_to_z_ifs(zs[j])
                tz = False
                
                k += 1
                zold = zs[j]
                zu[k] = zold
                m = 0
                while zu[k] == 0 and j + m < iip:
                    
                    zu[k] = zs[j+m]
                    m += 1
                    
            if obss[j] == 126:
                if ovs[j] > 150. and ovs[j] < 350:          
                    temp[k] = ovs[j]
                else:
                    temp[k] = np.nan
                if k < temp.shape[0] - 1:
                    
                    temp[k+1] = temp[k]
                tz = True
            elif obss[j] == 139:
                if ovs[j] == ovs[j] and np.abs(ovs[j]) < 150.:
                    u[k] = ovs[j]
                else:
                    u[k] = 0.
            elif obss[j] == 140 and np.abs(ovs[j]) < 150.:
                if ovs[j] == ovs[j]:
                    v[k] = ovs[j]
                else:
                    v[k] = 0.
            elif obss[j] == 117:
                if ovs[j] == ovs[j]:
                    gph[k] = ovs[j]
                    tz = True
                else:
                    pass #gph[k] = 0.
        
        if not tz: # neither z nor t are available
            gph[k] = p_to_z_ifs(zs[j])
        if zu[k] >= 1100.:
            zu[k] = zu[k-1] # probably a wrongly encoded missing surface pressure - replaced by next pressure above, likely ok for trajectory calculation
            
        timed = np.empty(k+1)
        timed.fill(np.nan)
        kstart = k // 2
        if kstart < 1:
            kstart = 1
        for i in range(kstart, k+1):
            if temp[i] != temp[i]:
                temp[i] = temp[i-1]
        if zu[k]<1100:   #np.any(~np.isnan(temp[k::-1])) and 
            if np.sum(np.isnan(temp[k::-1])>0.7*k):
                latd,lond,timed = rsd.numba_drift.trajectory(lats[ii],lons[ii], u[k::-1], v[k::-1],
                                                         zu[k::-1]*100, temp[k::-1])
            #print(zu[0], timed[k])
            #idy = np.where((zu[k::-1]==500) & (timed<0))[0]
            #idy = np.where(timed<0)[0]
            #if len(idy) >0:
                #x = 0
            isgph = False
            if np.all(np.isnan(timed[1:])):
                for i in range(k, -1, -1):
                    if gph[i] != gph[i] or np.abs(gph[i]) > 1.1 * p_to_z_ifs(zu[0]):
                        gph[i] = p_to_z_ifs(zu[i])
                latd,lond,timed = rsd.numba_drift.trajectory(lats[ii],lons[ii], u[k::-1], v[k::-1],
                                                         zu[k::-1]*100, temp[k::-1], zorig=gph[k::-1]/9.80665-0*station_elevation)
                isgph = True
        else:
            latd = np.zeros(zu[k::-1].shape[0])
            lond = np.zeros(latd.shape[0])
            timed = np.zeros(latd.shape[0])
            
        #idy = np.where((timed==0) & (zu[k::-1]==500.))[0]
        #idy = np.where(np.abs(latd)>10. )[0]
        #if np.any(np.abs(latd)>10.):
            
            #x = 0
        suspect = np.any((timed<1000) & (zu[k::-1]==500.) | (timed<0) | (timed>50000))
        #print(ii, iip, k, isgph,timed[zu[k::-1]==500.] )
        if suspect:
            x = 0
            if np.any(timed<0):
                timed[timed<0] = 0 # likely the result of unsorted geopotential values, data should be fixed.
            if np.any(timed>10000):
                if np.any(timed<10000):
                    
                    timed[timed>10000] = np.nanmax(timed[timed<10000])
                else:
                    timed[timed>10000] = 10000
            
        if len(latd) > 1:
            
            latd[-1] = latd[-2]
            lond[-1] = lond[-2]
            timed[-1] = timed[-2]
        else:
            if timed[0] !=timed[0]:
                timed[0] = 0
            
            
        #idy = np.where((timed<100) & (zu[k::-1]==500.) )[0]
        ##print(ii, iip, k, isgph,timed[zu[k::-1]==500.] )
        #if len(idy) >0 and len(timed) > 1:
            #x = 0

        zold = 1100.
        m = -1
        ostart = sl.start
        for j in range(iip-1, ii-1, -1):
            if zs[j] < zold:
                m += 1
                zold = zs[j]
            if m >= latd.shape[0]:
                out['latd'][ostart + j] = out['latd'][ostart + j + 1]  
                out['lond'][ostart + j] = out['lond'][ostart + j + 1] 
                out['timed'][ostart + j] = out['timed'][ostart + j + 1] 
                #print('m too large', m, latd.shape[0])
            else:
                

                if latd[m] == latd[m] and timed[m] == timed[m]:
                    out['latd'][ostart + j] = latd[m] 
                    out['lond'][ostart + j] = lond[m] 
                    out['timed'][ostart + j] = timed[m]
                else:
                    if ostart + j + 1 < out['latd'].shape[0]:
                        
                        out['latd'][ostart + j] = out['latd'][ostart + j + 1]  
                        out['lond'][ostart + j] = out['lond'][ostart + j + 1] 
                        out['timed'][ostart + j] = out['timed'][ostart + j + 1] 
            #if zs[j] ==500:
                #if out['timed'][ostart + j] == 0:
                    #x = 0
    
        #print(ii, iip, k, isgph,out['timed'][ostart + ii:ostart + iip][out['z_coordinate'][ostart + ii:ostart + iip]==85000] )
        #idy = np.where((out['timed'][ostart + ii:ostart + iip]<=0) & (out['z_coordinate'][ostart + ii:ostart + iip]==85000))[0]
        #if len(idy) >0:
            #x = 0
    #if False:    
        #plt.semilogy(u[:k+1], zu)
        #plt.semilogy(v[:k+1], zu)
        #plt.ylim(1000, 10)
        #plt.show()
        #plt.semilogy(temp[:k+1], zu)
        #plt.ylim(1000, 10)
        #plt.show()
        #plt.semilogy(latd[::-1], zu)
        #plt.semilogy(lond[::-1], zu)
        #plt.ylim(1000, 10)
        #plt.show()
        #plt.semilogy(timed[::-1], zu)
        #plt.ylim(1000, 10)
        #plt.show()
    
    if np.any(out['lond']>360.) or np.any(out['lond']<-360.):
        #print('wrong displacement',np.max(out['lond']),np.min(out['lond']) )
        out['lond'][:] = 0.
        out['latd'][:] = 0.
    return

@njit
def bisection(func, a,b, xref):

    fa = func(a, xref)
    fb = func(b, xref)
    mul = fa * fb
    if (mul >= 0 or np.isnan(mul)):
#        print("You have not assumed right a and b\n")
        return np.nan

    c = a
    while ((b-a) >= 0.01):

        # Find middle point
        c = (a+b)/2

        # Check if middle point is root
        fc = func(c, xref)
        if (fc == 0.0):
            break

        # Decide the side to repeat the steps
        if (fc*fa < 0):
            b = c
        else:
            a = c

    return c


@njit(cache=True)
def liquid(t):  # log of saturation pressure in Pa
    return -6096.9385 / t + 21.2409642 \
           - 2.711193e-2 * t \
           + 1.673952e-5 * t * t \
           + 2.433502 * np.log(t)

@njit(cache=True)
def liquid_hPa(t): # log of saturation pressure in hPa
    return -6096.9385 / t + 16.635794 \
           - 2.711193e-2 * t \
           + 1.673952e-5 * t * t \
           + 2.433502 * np.log(t)

@njit(cache=True)
def ice(t):
    return -6024.5282 / t + 29.32707 \
           + 1.0613868e-2 * t \
           - 1.3198825e-5 * t * t \
           - 0.49382577 * np.log(t)

@njit
def func(x, xref):
    return np.exp(liquid(x)) - xref

@njit
def func2(x, xref):
    return np.exp(ice(x)) - xref

@njit(cache=True)
def Sonntag(temp, over_water=True, over_ice=False):
    """Sonntag (1994): Sonntag, D., Advancements in the field of hygrometry,
    Meteorologische Zeitschrift, 3, 51-66, 1994. in Pa.
    Liquid: 273.16 < T < 373.15, Ice: 173.15 < T < 273.15

    ln ( ew )  = -6096.9385/t + 21.2409642
                 - 2.711193*10**(-2)*t
                 + 1.673952*10**(-5)*t**2
                 + 2.433502*math.log(t)

    ln ( ei )  = -6024.5282/t + 29.32707
                 + 1.0613868*10**(-2)*t
                 - 1.3198825*10**(-5)*t**2
                 - 0.49382577*math.log(t)

    Args:
        temp: air temperature in K
        liquid_only: use only water vapor over liquid water
        ice_only: use only water vapor over ice
        kwargs: dummy

    Returns:
         es : saturation water vapor pressure in Pa
    """


    if over_water:
        return np.exp(liquid(temp))
    elif over_ice:
        return np.exp(ice(temp))
    else:
        return np.where(temp < 273.16, np.exp(ice(temp)), np.exp(liquid(temp)))

@njit(cache=True)
def frostpoint_Sonntag(e):
    y = np.log(e / 611.213)
    return 273.15 + 13.7204 * y \
           + 7.36631e-1 * y * y \
           + 3.32136e-2 * y * y * y \
           + 7.78591e-4 * y * y * y * y

@njit(cache=True)
def dewpoint_Sonntag(e):
    """Calculation of Dewpoint after Sonntag 1994

    Parameters
    ----------
    e       Water vapor pressure [Pa]

    Returns
    -------
    dewpoint in [K]
    """
    #y = np.where(e > 0, np.log(e / 611.213), np.nan)
    y = np.log(e / 611.213)
    return 273.15 + 13.715 * y \
           + 8.4262e-1 * y * y \
           + 1.9048e-2 * y * y * y \
           + 7.8158e-3 * y * y * y * y

@njit(cache=True)
def vap2sh(e, p):
    """ Convert water vapor pressure to specific humidity
    Parameters
    ----------
    e      Water vapor [Pa]
    p      total air pressure [Pa]

    Returns
    -------
    specific humidity (1 = kg/kg)
    """
    #c = eps()  # Rd/Rv
    c = 2.8705e+2 / 4.6150e+2 #eps()  # Rd / Rv = 0.622
    pa = p - e  # dry air pressure
    return (e * c) / (e * c + pa)

@njit(cache=True)
def sh2vap(q, p):
    """Specific Humidity to Water vapor pressure
    formula derived from ratio of humid air vs. total air

    Parameters
    ----------
    q       specific humidity [kg/kg]
    p       total air pressure [Pa]

    Returns
    -------
    water vapor pressure [Pa]
    """
    #rd = constant(2.8705e+2, 'J/kg/K', 'gas constant air   ')
    #rv = constant(4.6150e+2, 'J/kg/K', 'gas constant H2O   ')
    
    
    c = 2.8705e+2 / 4.6150e+2 #eps()  # Rd / Rv = 0.622
    return (q * p) / (q + (1 - q) * c)

@njit(cache=True)
def ef(p, t=None, over_water=True, over_ice=False):
    """
    from
    Sugidachi, T. and Fujiwara, M.: Correction of the Stepwise Change Observed at 0C in
    Meisei RS2-91, RS-01G, and RS-06G Radiosonde Relative Humidity Profiles,
    Journal of the Meteorological Society of Japan. Ser. II, 91(3), 323-336,
    doi:10.2151/jmsj.2013-306, 2013.

    Args:
        t: air temperature K
        p: air pressure Pa

    Returns:
        f : enhancement factor for saturation water vapor pressure
            depending on air temperature and air pressure
    """
    # Vaisala / Buck 1981  / WMO 2008
    # 1.0016 + 3.15e-6*p - 0.074 / p  # kPa
    if over_water:
        return 1.0007 + (3.46e-6 * p / 100)
    if over_ice:
        return 1.0003 + (4.18e-6 * p / 100)
    if t is not None:
        return np.where(t < 273.16,
                        1.0003 + (4.18e-6 * p / 100),
                        1.0007 + (3.46e-6 * p / 100))

@njit(cache=True)
def svp(t, p):
    """
    Saturation water vapor pressure from Temperature
    The equations by Hyland and Wexler [4], the nearly identical equation by Wexler (1976, see reference below) and
    the equation by Sonntag [7] are the most commonly used equations among Radiosonde manufacturers
    and should be used in upper air applications to avoid inconsistencies.

    Known Saturation water Vapor Formulations:

    Bolton 1980
    Goff 1957, 1965        (180    - 273.15 / 273.15 - 373.15) 1957> WMO
    Goff and Gratch 1946   (184    - 273.15 / 273.15 - 373.15)
    Hyland and Wexler 1983 (173.16 - 273.15 / 273.15 - 473.15) Vaisala
    IAPWS 1995             (173.15 - 273.15 / 273.15 - 647   ) Wagner and Pruss + Wagner 1994
    Murphy and Koop 2005   (110    - 273.15 / 132      332   ) Ice formulation as well
    Sonntag 1990           (173.15 - 273.15 / 273.15 - 373.15)
    Sonntag 1994           (173.15 - 273.15 / 273.15 - 373.15)
    Wagner 1994            (190    - 273.15                  )
    Wagner and Pruss 1993  (190    - 273.15 / 273.15 - 647   )
    Wexler 1976            (                  273.15 - 373.15)
    Wright 1997

    Args:
        t: air temperature
        method: string or function
        p: pressure
        **kwargs: additional keywords passed to function

    Returns:
        es : saturation water vapor pressure in Pa
    """
    #if p is not None:
    f = ef(p)
    #else:
    #f = 1.
    return Sonntag(t) * f

@njit(cache=True)
def lvint(p, v, idx, ps):
    lp = np.log(p)
#    for j in range(v.shape[1]):
        
    for k in range(idx.shape[0]):
        j = idx[k]
        if np.sum(~np.isnan(v[:, j])) > 1:
            
            lastgoodpi = -1
            for i in range(0, v.shape[0]-1):
                if v[i, j] != v[i, j]:
                    if lastgoodpi == -1:
                        continue
                    kp =1
                    while i + kp < v.shape[0] and v[i+kp, j] != v[i+kp, j]:
                        kp += 1
                    if kp > 1:
                        if i + kp >= p.shape[0]:
                            continue
                        ti = np.searchsorted(ps, (p[lastgoodpi], p[i+kp]))
                        if ti[1] - ti[0] > 1:                           
                            continue
                    if lp[i+kp] - lp[lastgoodpi] > 0.:
                        
                        w0 = (lp[i] - lp[lastgoodpi]) / (lp[i+kp] - lp[lastgoodpi])
                    elif  lp[i+kp] - lp[lastgoodpi] < 0:
                        raise ValueError ('lvint: pressure not monotonic')
                    else:
                        w0 = 0.5
                    v[i, j] = (1 - w0) * v[lastgoodpi, j] + w0 * v[i+kp, j]
                else:
                    lastgoodpi = i
            ##if j==1 and v[i, j] == np.float32(29556.111):
                #print(i, j, v[i, j])
    return

@njit(cache=True)
def touv(wd, ws):
    wd[np.isnan(ws)] = np.nan
    wd[ws < 0.001] = 0.
    ws[np.isnan(wd)] = np.nan
    u = ws * np.cos(np.radians(270.-wd))
    v = ws * np.sin(np.radians(270.-wd))
    return u, v, wd, ws

@njit(cache=True)
def wswd(u, v):

    wd = - np.arctan2(v, u) * 180 / np.pi - 90.
    wd[wd<0] += 360.
    return np.sqrt(u ** 2 + v ** 2), wd 
    

@njit(cache=True, boundscheck=False)
def augment1(obstab, a_obstab, loaded_feedback, a_loaded_feedback,
            ri, ts, ps, humvar, wvar, zvar, fn):


    recordindex=np.empty(ri.shape[0],obstab['date_time'].dtype)
    recordtimestamp=np.empty(ts.shape[0],obstab['date_time'].dtype)
    rim = np.max(ri[1:]-ri[:-1])
    zlist = np.empty(rim, dtype=np.float32)
    ztype = np.empty(rim, dtype=np.int32)
    #allvar = np.concatenate(((0, 117, 126), humvar, wvar))
    allvar = np.zeros(11, dtype=np.int32)
    allvar[1:3] = np.array((117, 126))
    allvar[3:7] = humvar
    allvar[7:] = wvar
    alli = np.zeros(np.max(allvar) + 1, dtype=np.int32)
    
    
    for i in range(len(allvar)):
        alli[allvar[i]] = i

    j=-1 # augmented index
    jsave=0
    p= obstab['z_coordinate'][0]-1 # z_coordinate[0]-1
    rts= obstab['date_time'][0] # date_time[0]

    j=-1
    porig = -1
    #good=True
    spurious_z = 0
    z_coordinate = obstab['z_coordinate']
    g = 9.80665
    oovs = obstab['observation_value'].shape[0]
    k = 0
    ak = 0
    while k <oovs:
        i = k
        j = 0
        plist = np.empty(zlist.shape[0]+15)
        isave = np.empty(zlist.shape[0]+15, dtype=np.int32)
        isave[:] = -1
        while i < oovs and obstab['date_time'][i] ==rts:

            
            if z_coordinate[i]!=porig:
                porig = z_coordinate[i]
                zct = obstab['z_coordinate_type'][i]
                oov = obstab['observed_variable'][i]
                if zct == 1:
                    
                    plist[j]=z_coordinate[i]
                    ztype[j] = zct
                elif oov == 0:
                    plist[j]=z_coordinate[i]
                    ztype[j] = zct
                elif zct == 0: 
                    zlist[j] = z_coordinate[i]*g
                    plist[j] = z_to_p_ifs(zlist[j])
                    z_coordinate[i] = plist[j]
                    ztype[j] = 1
                elif zct == 2: 
                    zlist[j] = z_coordinate[i]
                    plist[j] = z_to_p_ifs(zlist[j])
                    z_coordinate[i] = plist[j]
                    ztype[j] = 1
                else:
                    spurious_z += 1
                    plist[j] = z_coordinate[i]
                    zlist[j] = np.nan
                    ztype[j] = zct
                isave[j] = i
                j += 1
            else:
                z_coordinate[i] = z_coordinate[i-1]
            i += 1
        
        
        
        idy = np.where((ps>np.min(plist[:j])) & (ps<np.max(plist[:j])))
        if len(idy[0]) > 0:
            plist[j:j+len(idy[0])] = ps[idy]
        idx = np.argsort(plist[:j+len(idy[0])])
        plist = plist[idx]
        isave = isave[idx]
                         
        
        fobs = np.empty((isave.shape[0], allvar.shape[0]), dtype=np.float32) # +15 to accomodate for additional standard pressure levels
        ffg_depar = np.empty_like(fobs) # +15 to accomodate for additional standard pressure levels
        fan_depar = np.empty_like(fobs)
        fbc = np.empty_like(fobs)
        fbc_fg = np.empty_like(fobs)
        fcf = np.empty((isave.shape[0], allvar.shape[0]), dtype=np.int32)
        fcm = np.empty_like(fcf)
        
        for f in fobs, ffg_depar, fan_depar, fbc: 
            f.fill(np.nan)
        
        i = k
        j = -1
        porig = -1
        while i < oovs and obstab['date_time'][i] ==rts:
          
            if z_coordinate[i]!=porig:
                j += 1
                porig = z_coordinate[i]
                if z_coordinate[i] != z_coordinate[i]:
                    ii = np.searchsorted(plist, porig)
                    if ii == plist.shape[0]:
                        ii -= 1
                else:
                    ii = j
                    #print('x', plist[-1], porig)
            vi = alli[obstab['observed_variable'][i]]
            fobs[ii, vi] = obstab['observation_value'][i]
            ffg_depar[ii, vi] = loaded_feedback['fg_depar@body'][i]
            fan_depar[ii, vi] = loaded_feedback['an_depar@body'][i]
            fbc[ii, vi] = loaded_feedback['biascorr@body'][i]
            fbc_fg[ii, vi] = loaded_feedback['biascorr_fg@body'][i]
            i += 1

        
        #for ii in range(isave.shape[0]):
            #if fobs[ii, 9] == fobs[ii, 9] and  fobs[ii, 10] == fobs[ii, 10]:
                
        mask =   ~np.isnan(fobs[:, 9]+   fobs[:, 10]) & (np.isnan(fobs[:, 7]+fobs[:, 8]))
        if np.any(mask):
            fobs[mask, 7], fobs[mask, 8] = touv(fobs[mask, 10] , fobs[mask, 9])
            fcm[mask] = 1
            fcf[mask] = 0
        
        # dpd present, rh missing
        mask =   ~np.isnan(fobs[:, 3]+   fobs[:, 2]) & (np.isnan(fobs[:, 5]))
        if np.any(mask):
            fobs[mask, 4] = fobs[:, 2] - fobs[mask, 3]
            fobs[mask, 5] = Sonntag(fobs[mask, 4]) / Sonntag(fobs[mask, 2]) *100.
            fcm[mask, 5] = 2
            fcf[mask, 5] = 0

        # td present, rh missing
        mask =   ~np.isnan(fobs[:, 4]+   fobs[:, 2]) & (np.isnan(fobs[:, 5]))
        if np.any(mask):
            fobs[mask, 5] = Sonntag(fobs[mask, 4]) / Sonntag(fobs[mask, 2]) *100.
            fcm[mask, 5] = 2
            fcf[mask, 5] = 0

        # sh present, rh missing
        mask =   ~np.isnan(fobs[:, 6]+   fobs[:, 2]) & (np.isnan(fobs[:, 5]))
        if np.any(mask):
            vpdata = sh2vap(fobs[mask, 6], plist[mask])
            fobs[mask, 5] = vpdata / svp(fobs[:, 2], p=plist[mask]) *100.
            fcm[mask, 5] = 4
            fcf[mask, 5] = 0
            
        #rh present, td, dpd missing
        mask =   ~np.isnan(fobs[:, 5]+   fobs[:, 2]) & (np.isnan(fobs[:, 3]) | np.isnan(fobs[:, 4]))
        if np.any(mask):
            vpdata = fobs[mask, 5] * np.exp(liquid(fobs[mask, 2])) /100. #Sonntag(fobs[mask, 2])
            fobs[mask, 3] = dewpoint_Sonntag(vpdata)
            fobs[mask, 4] = fobs[mask, 2] - fobs[mask, 3]
            fcm[mask, 3] = 3
            fcf[mask, 3] = 0
            fcm[mask, 4] = 3
            fcf[mask, 4] = 0

        # rh present, q missing
        mask =   ~np.isnan(fobs[:, 5]+   fobs[:, 2]) & (np.isnan(fobs[:, 6]))
        if np.any(mask):
            vpdata = fobs[mask, 5] * Sonntag(fobs[mask, 2]) / 100.
            fobs[mask, 6] = vap2sh(vpdata, plist[mask])
            fcm[mask, 6] = 3
            fcf[mask, 6] = 0


        for f in fobs, ffg_depar, fan_depar, fbc:
            lvint(plist, f[:isave.shape[0]], np.array((0, 1, 2, 4, 5, 7, 8)))
        
        #fobs[:, 6] =qs(fobs[:, 2], plist, fobs[:, 5])        
        #fobs[:, 5] =td(fobs[:, 2], plist, fobs[:, 4])        
        mask =   ~np.isnan(fobs[:, 7]+   fobs[:, 8]) & (np.isnan(fobs[:, 9]+fobs[:, 10]))
        if np.any(mask):
            fobs[:, 10], fobs[:, 9] = wswd(fobs[:, 7] , fobs[:, 8])
            fcm[mask, 9:] = 2
            fcf[mask, 9:] = 0


        fobs[:, 3] = fobs[:, 2] - fobs[:, 4]
        
        for j in range(len(plist)):
            for jj in range(11):
                if fobs[j, jj] == fobs[j, jj]:
                    a_obstab['observation_value'][ak] = fobs[j, jj]
                    
                    a_obstab['z_coordinate'][ak] = plist[j]
                    a_obstab['z_coordinate_type'][ak] = 1
                    a_obstab['observed_variable'][ak] = allvar[jj]        
                    a_obstab['date_time'][ak] = rts
                    a_obstab['conversion_flag'][ak] = 0
                    a_obstab['conversion_method'][ak] = 0
                    a_loaded_feedback['fg_depar@body'][ak] = ffg_depar[j, jj]
                    a_loaded_feedback['an_depar@body'][ak] = fan_depar[j, jj]
                    a_loaded_feedback['biascorr@body'][ak] = fbc[j, jj]
                    a_loaded_feedback['biascorr_fg@body'][ak] = fbc_fg[j, jj]
                    ak += 1
                
            
        k = i
        if i < oovs:
            
            rts=obstab['date_time'][i]


    return a_obstab, a_loaded_feedback, ak


@njit(cache=True, boundscheck=True)
def fill_obsid(avar,conversion_flag):
    for o in range(avar.shape[0]):
        if conversion_flag[o] == 0:
            for i in range(2):
                avar[o,i]=b'9'
    return avar

@njit(cache=True, boundscheck=True)
def fill_restdata(final, rest_data, ari, ri, fn, i): #,test):

    
    for l in range(ari.shape[0]-1):
        
        final[ari[l]:ari[l+1]] =  rest_data[ri[l]] 

    final[ari[-1]:] =  rest_data[ri[-1]] 

    return final

@njit(cache=True, boundscheck=True)
def fill_restdata_abs(final, rest_data, ari, ri, addedvar, absidx): #,test):

    
    k = -1
    ash =  ari.shape[0]
    for l in range(final.shape[0]):
        if k + 1 < ash and l >= ari[k+1]:
            k += 1
            #if k == ri.shape[0]:
                #k -= 1
        
            m =ari[k]
            mmax = ari[k+1] if k + 1 < ash else final.shape[0]
            while m < mmax and addedvar[m] == -1:
                m += 1
        final[l] = rest_data[addedvar[m]]
    
    final = final[absidx]
        
    return final

@njit(cache=True, boundscheck=True)
def fill_restdata_addedvar(final, rest_data, addedvar): #,test):

    
    for l in range(addedvar.shape[0]):
        if addedvar[l] != -1:
            final[l] = rest_data[addedvar[l]]
            

    return final


@njit(cache=True, boundscheck=True)
def fill_restdata_addedvar_abs(final, rest_data, addedvar, absidx, miss_val): #,test):

    
    for l in range(addedvar.shape[0]):
        if addedvar[absidx[l]] != -1:
            final[l] = rest_data[addedvar[absidx[l]]]
        else:
            final[l] = miss_val
            
    return final

@njit(cache=True, boundscheck=True)
def fill_restdata_obsid_abs(final, rest_data, addedvar, absidx, dum): #,test):

    m = 0
    #c = np.concatenate((rest_data[0][:4] , np.array(['9', '9'],dtype='S1') , np.string_([f'{m:014d}']).view('S1')))
    for l in range(addedvar.shape[0]):
        if addedvar[absidx[l]] != -1:
            final[l] = rest_data[addedvar[absidx[l]]]
        else:
            final[l] =dum[m]
            m += 1
            
    return final

@njit(cache=True, boundscheck=False)
def fill_restdata_old(final, rest_data, addedvar, j): #,test):

    for l in range(addedvar.shape[0]-1):
        diff=(addedvar[l+1,1]-addedvar[l,1])-(addedvar[l+1,0]-addedvar[l,0])
        final[addedvar[l,1]:addedvar[l+1,1]-diff]=rest_data[addedvar[l,0]:addedvar[l+1,0]]
        if diff>0:
            for i in range(diff):
                final[addedvar[l+1,1]-diff+i]=rest_data[addedvar[l+1,0]-1]
        #x=np.nansum(final[addedvar[l,1]:addedvar[l+1,1]]-test[addedvar[l,1]:addedvar[l+1,1]])
        #if ~np.isnan(x) and x!=0:
                        #print(addedvar[l:l+2])                             
                        #print('x')

    return final

def split(x): 
    return [i.encode() for i in x.decode()]

@njit(cache=True, boundscheck=False)
def polinl(XA,YA,X, Y):


#    REAL XA(4),YA(NFELD,4),X,Y(NFELD)
#    REAL XX1,XX2,XX3,XX4,X1X2

    XX1=X-XA[0]
    XX2=X-XA[1]
    XX3=X-XA[2]
    XX4=X-XA[3]
    X1X2=XA[0]-XA[1]
    X1X2=X1X2*X1X2*X1X2

    for I in range(Y.shape[0]):
        Y[I]=XX2*XX3*XX4/3.0*YA[I,0]-XX1*XX3*XX4*YA[I,1]+XX1*XX2*XX4*YA[I,2]-XX1*XX2*XX3/3.0*YA[I,3]
        Y[I]=Y[I]/X1X2/2.0

    #for I in range(Y.shape[0]):
        #if Y[I]>YA[I,1] and Y[I]> YA[I,2]:
            #Y[I]=np.max((YA[I,1],YA[I,2]))
        #if Y[I] < YA[I,1] and Y[I] < YA[I,2]:
            #Y[I]=np.min((YA[I,1],YA[I,2]))

    return Y

@njit
def is_sorted(a):
    
    for i in range(a.size-1):
        if a[i+1] < a[i] :
            return False
    return True

@njit(cache=True, boundscheck=True)
def trisimple(tera5, lons, lats, ts, xref, yref, secref, glon, glat, secs, ret):
    
    Y = np.zeros(4, dtype=np.float32)
    yval = np.zeros((1, 4), dtype=np.float32)
    pval = np.zeros(2, dtype=np.float32)
    dx = glon[1] - glon[0]
    dy = glat[1] - glat[0]
    w = np.zeros(2)
    
    #sold = 0
    for i in range(0, tera5.shape[2]//2):
        
        for it in range(2):

                YA = tera5[:,:,2*i+ it, 0]
                XA = glon[xref[i, 1]] + np.arange(-1, 3) *dx
                y = glat[yref[i, 1]] + np.arange(-1, 3)*dy
                yval[0, :] = polinl(XA,YA,lons[i], Y)
                
                xx = polinl(y, yval, lats[i], Y)
                #if xx[0] > tmax or xx[0] < tmin:
                    #print('fail',xx[0], tmax,tmin )
                pval[it] = xx[0] #polinl(tlat, yval, lats[il])          

        s =  secs[secref[i]]
        if secref[i]+1 < secs.shape[0]:
            
            s1 = secs[secref[i]+1]
            if (s1 - ts[i]) > (s1 - s) or (s1 - ts[i]) < 0:
                #raise ValueError(i,s1, s, ts[i])
                ret[i] = np.nan
                #sold = s1 - s
            else:
            
                w[0]=(s1-ts[i])/(s1 - s)
                w[1]=1.-w[0]
                if w[0] < 0 or w[0]> 1:
                    raise ValueError(str(w[0]))
                
                ret[i]=np.float32(pval[0]*w[0]+pval[1]*w[1])
        else:
            ret[i] = pval[0]
        
    return ret
    
#@njit(cache=False, boundscheck=False)
def triint(tera5, reatab, tlon, tlat, lons, lats, secs, ts, ori, obss,press, zs, p, dy, fn, ret):
    sh = tera5.shape
    tmax = np.max(tera5)
    tmin = np.min(tera5)
    pval = np.zeros((sh[2],sh[3]))
    pval.fill(np.nan)
    yval = np.zeros((1, 4))
    
    idx = np.zeros(len(reatab), dtype=np.int32)
    oldlats = -100.
    oldlons = -400.
    oldts = -1
    im = 0
    if ret == 0:
        return idx[:im]
    its = np.searchsorted(secs, ts) - 1
    its[its==-1] = 0
    if ret == 1:
        return idx[:im]
    iyref = np.searchsorted(-tlat, -lats) - 2
    hlon = tlon.copy()
    lons[lons>=359.9999] -= 360.
    if ret == 2:
        return idx[:im]
    if is_sorted(hlon):
        lons[lons > hlon[-2]] -= 360.
        lons[lons < hlon[0]] += 360.
        ixref = np.searchsorted(hlon, lons) - 2
             
    else:
        dx = np.max(tlon[1:] - tlon[:-1])
        hlon[:] = hlon[0] - 360 + np.arange(hlon.shape[0]) *dx
        lons[lons > hlon[-2]] -= 360.
        lons[lons < hlon[0]] += 360.
        ixref = np.searchsorted(hlon, lons) - 2

    if ret == 3:
        return idx[:im]
    if np.any(ixref < 0) or np.any(ixref == hlon.shape[0]-2):
        #print(fn)
        #if not np.any(lons==hlon[ixref+2]):
           #raise ValueError(fn) #print('fail')
        #else:
        ixref[ixref==-1] = 0
        if np.any(ixref < 0) or np.any(ixref == hlon.shape[0]-2):
            return idx[:im] #raise ValueError(fn) #print('fail')
    
    if ret ==4:    
        return idx[:im]    
    Y = np.zeros(4)
    yval = np.zeros((1, 4))
    
    for il in range(len(reatab)):
        if obss[il] == p:
            found = False
            for ip in range(sh[3]):
                if press[ip] == zs[il]:
                    found = True
                    break
            if not found:
                continue
            #if ts[il] != oldts:
                
                #its =np.searchsorted(secs, ts[il]) - 1
                #if its == -1:
                    #its = 0
                #oldts = ts[il]
            
            #if its + 1 >= tera5.shape[2]:
                #continue
            its1 = its[il] + 1
            if its1 == secs.shape[0]:
                continue
            w0=(secs[its1]-ts[il])/(secs[its1] - secs[its[il]])
            w1=1.-w0
            
            hh = hlon[ixref[il]:ixref[il]+4]
            
            for it in its[il], its1:
                if lats[il] > 90 - 2 * dy:
                    pval[it, ip] = tera5[0, 0, it, ip]          
                elif lats[il] < -90 + 2 * dy:
                    pval[it, ip] = tera5[-1, 0, it, ip]          
                else:
                    i = 0
                    
                    YA = tera5[iyref[il]:iyref[il]+4, ixref[il]:ixref[il]+4, it, ip]
                    yval[0, :] = polinl(hh,YA,lons[il], Y)
                    #for iy in range(4):
                        ##print(yval.dtype, tlon.dtype, tera5[iy, :, it, ip:ip+1].T.dtype, lons[il])
                        ##print(yval.shape, tlon.shape, tera5[iy, :, it, ip:ip+1].T.shape)
                        ##if ixref == 2:
                            ##print(ixref)
                        #if ixref[il] + 4 > tera5.shape[1]:
                            ##xx = polinl(hlon[0:0+4],tera5[iyref[il]+iy, 0:0+4, it, ip:ip+1].T,lons[il])
                            #raise ValueError #print('ixref was too large!')
                        #else:
                            ##YA = tera5[iyref[il]+iy, ixref[il]:ixref[il]+4, it, ip:ip+1].T
                            ##Y = np.zeros(YA.shape[0])
                            
                            #xx = polinl(hlon[ixref[il]:ixref[il]+4],tera5[iyref[il]+iy, ixref[il]:ixref[il]+4, it, ip:ip+1].T,lons[il], Y)
                        ##print(xx.shape, xx.dtype)
                        ##if xx[0] > tmax or xx[0] < tmin:
                            ##print('fail',xx[0], tmax,tmin )
                        #yval[i, iy] = xx[0] #polinl(tlon,tera5[iy, :, it, ip:ip+1].T,lons[il])
                    xx = polinl(tlat[iyref[il]:iyref[il]+4], yval, lats[il], Y)
                    #if xx[0] > tmax or xx[0] < tmin:
                        #print('fail',xx[0], tmax,tmin )
                    pval[it, ip] = xx[0] #polinl(tlat, yval, lats[il])          
    
            reatab[im]=pval[its[il],ip]*w0+pval[its1,ip]*w1
            idx[im] = il
            im += 1
    
    return idx[:im]


@njit(cache=True, boundscheck=True)
def add_fb(loaded_obstab,loaded_feedback,ref20CR,refera5an,refera5fc):
    i20=0
    iera=0
    for i in range(loaded_obstab['date_time'].shape[0]):
        if loaded_feedback['fg_depar@body'][i]!=loaded_feedback['fg_depar@body'][i]:
            if loaded_obstab['observation_value'][i]==loaded_obstab['observation_value'][i]:
                if refera5fc[i]==refera5fc[i]:
                    loaded_feedback['fg_depar@body'][i]=loaded_obstab['observation_value'][i]-refera5fc[i]
                    loaded_feedback['an_depar@body'][i]=loaded_obstab['observation_value'][i]-refera5an[i]
                    loaded_feedback['biascorr@body'][i]=0.
                    iera+=1
                elif ref20CR[i]==ref20CR[i]:
                    loaded_feedback['fg_depar@body'][i]=loaded_obstab['observation_value'][i]-ref20CR[i]
                    loaded_feedback['an_depar@body'][i]=loaded_obstab['observation_value'][i]-ref20CR[i]
                    loaded_feedback['biascorr@body'][i]=0.
                    i20+=1
        #if i%1000000==0:
            #print(i,i20,iera)

@njit(cache=True, boundscheck=False)
def make_vrindex(vridx,ridx): # this function is similar to np.unique with return_index=True, but it expands the index to the  
    # original array dimensions
    l=0
    for i in range(len(ridx)): # to set the recordindices
        if i == 0:
            l +=1
        else:
            if ridx[i]>ridx[i-1]:
                vridx[ridx[i-1] + 1 : ridx[i] + 1]=l # next record after l
                l += 1
            else:
                l += 1
    vridx[ridx[i] + 1 :]=l #len(idx) # next record for the last element is the len of the data

@njit(cache=True, boundscheck=True)
def augment2(obstab, a_obstab, loaded_feedback, a_loaded_feedback,
            ri, ts, ps, humvar, wvar, zvar, fn):


    rimp = np.argmax(ri[1:]-ri[:-1])
    u = np.unique(obstab['observed_variable'][ri[rimp]:ri[rimp+1]]).shape[0]
    rim = ri[rimp+1]-ri[rimp]
    #zlist = np.empty(ri[rim+1]-ri[rim], dtype=np.float32)

    allvar = np.zeros(11, dtype=np.int32)
    allvar[1:3] = np.array((117, 126))
    allvar[3:7] = humvar
    allvar[7:] = wvar
    allvarsurf = allvar.copy()
    alli = np.zeros(np.max(allvar) + 1, dtype=np.int32)
    a_obstab['observation_value'][:] = np.nan
    
    for i in range(len(allvar)):
        alli[allvar[i]] = i

    j=-1 # augmented index
    jsave=0
    p= obstab['z_coordinate'][0]-1 # z_coordinate[0]-1
    rts= obstab['date_time'][0] # date_time[0]

    j=-1
    porig = -1
    spurious_z = 0 
    z_coordinate = obstab['z_coordinate']
    ztype = obstab['z_coordinate_type']
    g = 9.80665
    oovs = obstab['observation_value'].shape[0]
    oidx = np.arange(oovs, dtype=np.int64)
    orig = np.empty(a_obstab.shape[0], dtype=np.int64)
    orig.fill(-1)
    ### adjust for relative humidity in %
    #idx = np.where((obstab['observed_variable']==138)&(obstab['observation_value']>2.0))
    #if len(idx) > 0:
        #obstab['observation_value'][idx] /=100.
    #idx = np.where((obstab['observed_variable']==138)&(obstab['observation_value']<0.))
    #if len(idx) > 0:
        #obstab['observation_value'][idx] = np.nan
       
    k = 0
    ak = 0
    #m = 0
    plc = 0
    tddmax = 0.
    while k < oovs:
        i = k
        j = 0
        plist = np.empty(int(rim*12/u)+ 16, dtype=np.float32)
        #ztype = np.empty(zlist.shape[0]+16000, dtype=np.int32) 
        isave = np.empty(int(rim*12/u)+16, dtype=np.int32)
        isave[:] = -1
        if rts == 2967624000:
            x = 0
        while i < oovs and obstab['date_time'][i] ==rts:

            #if i == 84:
                #x = 0
            if z_coordinate[i]!=porig or i == k:
                porig = z_coordinate[i]
                zct = obstab['z_coordinate_type'][i]
                oov = obstab['observed_variable'][i]
                if zct == 1:
                    
                    plist[j]=z_coordinate[i]
                    ztype[i] = zct
                elif oov == 0:
                    plist[j]=z_coordinate[i]
                    if z_coordinate[i] == z_coordinate[i]:
                        if ztype[i] == 2:
                            phi = z_coordinate[i]
                            plist[j] = z_to_p_ifs(phi)
                            z_coordinate[i] = plist[j]
                            
                            ztype[i] = 1
                    else:
                        ztype[i] = 1
                elif zct == 0: 
                    phi = z_coordinate[i]*g
                    plist[j] = z_to_p_ifs(phi)
                    z_coordinate[i] = plist[j]
                    #if obstab['observed_variable'][i] == 126:
                        #print('126', obstab['observation_value'][i] )
                    ztype[i] = 1
                elif zct == 2: 
                    phi = z_coordinate[i]
                    plist[j] = z_to_p_ifs(phi)
                    z_coordinate[i] = plist[j]
                    ztype[i] = 1
                    #if obstab['observed_variable'][i] == 126:
                        #print('126', obstab['observation_value'][i] )
                else:
                    spurious_z += 1
                    plist[j] = z_coordinate[i]
                    #zlist[j] = np.nan
                    ztype[i] = zct
                isave[j] = i
                j += 1
            else:
                z_coordinate[i] = z_coordinate[i-1]
                ztype[i] = ztype[i-1]
            if ztype[i] != 1:
                x = 0
            i += 1
        
        
        ix = np.where((obstab['z_coordinate_type'][k:i]!=1))[0]
        if len(ix) > 0:
            iy = np.where(obstab['observed_variable'][k + ix]!=0)[0] # probably a surface observation                
            obstab['observation_value'][k + ix] = np.nan
            if np.any(~np.isnan(obstab['z_coordinate'][k + ix])) and len(iy) > 0:
                raise ValueError(fn+' '+'z_coordinate_type '+str(len(ix))+', ' +str(obstab.shape[0]))
            else:
                pass
                #k = i
                #if i < obstab['date_time'].shape[0]:
                    
                    #rts = obstab['date_time'][i]
                #continue
                
            #return obstab, loaded_feedback, ak, orig[:ak],tddmax #addedvar[:l+1]
            #break
      
        
        if j == 0:
            k = i
            if i < obstab['date_time'].shape[0]:
                raise ValueError('augment2 - should never be there')
                    
                rts = obstab['date_time'][i]
            continue
        plj = plist[:j]
        pmax = np.max(plj)
        pmin = np.min(plj)
        
        
        for ip in range(ps.shape[0]):
            if ps[ip] > pmin and ps[ip] < pmax:
                if ps[ip] not in plj:
                    plist[j] = ps[ip]
                    j += 1
        
        idx = np.argsort(plist[:j])
        plist = plist[idx]
        isave = isave[idx]

        fobs = np.empty((isave.shape[0], allvar.shape[0]), dtype=np.float32) # +15 to accomodate for additional standard pressure levels
        fidx = np.empty((isave.shape[0], allvar.shape[0]), dtype=np.int64) # +15 to accomodate for additional standard pressure levels
        fobsvar = np.empty((isave.shape[0], allvar.shape[0]), dtype=np.int64) # +15 to accomodate for additional standard pressure levels
        fidx.fill(-1)
        fobsvar.fill(-1)
        ffg_depar = np.empty_like(fobs) # +15 to accomodate for additional standard pressure levels
        fan_depar = np.empty_like(fobs)
        fbc = np.empty_like(fobs)
        #fbc_fg = np.empty_like(fobs)
        fcf = np.zeros((isave.shape[0], allvar.shape[0]), dtype=np.int32)
        fcm = np.empty_like(fcf)
        int_min = -2147483647 - 1
        fcf[:] = int_min
        fcm[:] = int_min
        
        for f in fobs, ffg_depar, fan_depar, fbc: 
            f.fill(np.nan)
        
        i = k
        j = -1
        porig = -1
        while i < oovs and obstab['date_time'][i] ==rts:
          
            if z_coordinate[i]!=porig:
                j += 1
                porig = z_coordinate[i]
                if z_coordinate[i] == z_coordinate[i]:
                    ii = np.searchsorted(plist, porig)
                    if ii == plist.shape[0]:
                        ii -= 1
                else:
                    ii = j
                    
            if plist[ii] == porig or z_coordinate[i] != z_coordinate[i]:
                
                    ##print('x', plist[-1], porig)
                vi = alli[obstab['observed_variable'][i]]
                fobs[ii, vi] = obstab['observation_value'][i]
                fobsvar[ii, vi] = obstab['observed_variable'][i]
                fidx[ii, vi] = oidx[i]
                ffg_depar[ii, vi] = loaded_feedback['fg_depar@body'][i]
                fan_depar[ii, vi] = loaded_feedback['an_depar@body'][i]
                fbc[ii, vi] = loaded_feedback['biascorr@body'][i]
                #fbc_fg[ii, vi] = loaded_feedback['biascorr_fg@body'][i]
            
            i += 1
            
        foo = fobs.copy()
        
        ## convert ws,wd to u,v before vertical interpolation        
        mask =   np.isnan(fobs[:, 7]+fobs[:, 8])
        if np.any(mask):
            fobs[mask, 7], fobs[mask, 8], fobs[mask, 9], fobs[mask, 10] = touv(fobs[mask, 9] , fobs[mask, 10])
            fcm[mask] = 1
            fcf[mask] = 0

        for f in fobs, ffg_depar, fan_depar, fbc:
            lvint(plist, f[:isave.shape[0]], np.array((0, 1,2,3, 4, 5, 7, 8)), ps) # interpolate vertically all variables except specific humidity, wind speed and wind direction
        
                    
        tmask = ~np.isnan(fobs[:, 2]) & (fobs[:, 2] > 150) & (fobs[:, 2] < 350)


        mask = np.searchsorted(plist, 10000.)
        if mask < plist.shape[0] and plist[mask] == 10000.:# and ~np.isnan(fobs[mask, 4]):
            
            vals = fobs[mask, :].flatten()
            vpdata = fobs[mask, 5]* np.exp(liquid(fobs[mask, 2])) 
            td = bisection(func, fobs[mask, 2]-70., fobs[mask, 2]+1., vpdata)
            if np.abs(vals[4] - td) > tddmax:
                tddmax = np.abs(vals[4] - td)
        
        if np.sum(tmask) > 0:
            
            ## dpd present, rh missing
            #mask =   tmask & (~np.isnan(fobs[:, 3])) & (np.isnan(fobs[:, 5]))
            #if np.any(mask):
                #dp = fobs[mask, 2] - fobs[mask, 3]
                #fobs[mask, 5] = Sonntag(dp) / Sonntag(fobs[mask, 2]) 
                #fcm[mask, 5] = 2
                #fcf[mask, 5] = 0
    
            # dpd present, dp missing
            mask =  tmask & ( ~np.isnan(fobs[:, 3])) & (np.isnan(fobs[:, 4]))
            if np.any(mask):
                fobs[mask, 4] = fobs[mask, 2] - fobs[mask, 3]
                fcm[mask, 4] = 2
                fcf[mask, 4] = 0

            # dp present, dpd missing
            mask =  tmask & ( ~np.isnan(fobs[:, 4])) & (np.isnan(fobs[:, 3]))
            if np.any(mask):
                fobs[mask, 3] = fobs[mask, 2] - fobs[mask, 4]
                if np.any(fobs[mask, 3]<-0.2):
                    x = 0
                fcm[mask, 3] = 2
                fcf[mask, 3] = 0

            # td present, replace rh EVEN IF AVAILABLE
            mask =  tmask & ( ~np.isnan(fobs[:, 4])) & ( fobs[:, 4] < fobs[:, 2] + 0.1 ) #& (np.isnan(fobs[:, 5]))
            if np.any(mask):
                fobs[mask, 5] = Sonntag(fobs[mask, 4]) / Sonntag(fobs[mask, 2]) 
                fcm[mask, 5] = 2
                fcf[mask, 5] = 0
    
            # q present, replace rh EVEN IF AVAILABLE, only if DP,DPD are not available
            mask = tmask & (  ~np.isnan(fobs[:, 6])) & (plist > 0.) & (fcm[:, 5] != 2) #& (np.isnan(fobs[:, 5]))
            if np.any(mask):
                vpdata = sh2vap(fobs[mask, 6], plist[mask])
                fobs[mask, 5] = vpdata / svp(fobs[mask, 2], p=plist[mask])
                fcm[mask, 5] = 4
                fcf[mask, 5] = 0
                
    
            #rh present, td, dpd missing
            mask =  tmask & ( ~np.isnan(fobs[:, 5])) & np.isnan(fobs[:, 4]) & np.isnan(fobs[:, 3])
            #mask =  tmask & ( ~np.isnan(fobs[:, 5])) #&  np.isnan(fobs[:, 3])
            if np.any(mask):
                
                vpdata = fobs[mask, 5]* np.exp(liquid(fobs[mask, 2])) #* Sonntag(fobs[mask, 2])# * np.exp(liquid(fobs[mask, 2]))  #Sonntag(fobs[mask, 2])
                vpdata2 = fobs[mask, 5]* np.exp(ice(fobs[mask, 2])) #* Sonntag(fobs[mask, 2])# * np.exp(liquid(fobs[mask, 2]))  #Sonntag(fobs[mask, 2])
                tm = fobs[mask, 2] < 273.15
                vpdata[tm] = vpdata2[tm]
                l = 0
                for il in range(len(mask)):
                    if mask[il]:
                        #vpdata = fobs[mask, 5] * np.exp(liquid(fobs[mask, 2]))  #Sonntag(fobs[mask, 2])
                        if fobs[il, 2] < 273.15:
                            
                            vpdat = bisection(func, fobs[il, 2]-70., fobs[il, 2]+1., vpdata[l])
                        else:
                            vpdat = bisection(func, fobs[il, 2]-70., fobs[il, 2]+1., vpdata[l])
                        fobs[il, 4] = vpdat
                        fobs[il, 3] = fobs[il, 2] - vpdat
                        #fobs[il, 3] = vpdat
                        l += 1

                if np.any(fobs[mask, 3]<-0.2):
                    x = 0
                fcm[mask, 3] = 3
                fcf[mask, 3] = 0
                fcm[mask, 4] = 3
                fcf[mask, 4] = 0
    
            # rh present, q missing, #REPLACE even if available
            mask =  tmask & ( ~np.isnan(fobs[:, 5])) #& (np.isnan(fobs[:, 6])) & (plist > 0.)
            if np.any(mask):
                vpdata = fobs[mask, 5] * Sonntag(fobs[mask, 2])
                    
                if np.any(vpdata>0.2*plist[mask]):
                    #print('vapor pressure too large')
                    fobs[mask, 6] = np.nan
                else:
                    fobs[mask, 6] = vap2sh(vpdata, plist[mask])
                if np.any((fobs[mask, 6]<0)): #| (plist[mask]==50000.) &(fobs[mask, 6]>0.004)):# &(plist[mask]==50000.)):
                    raise ValueError('negative specific humidity')

                fcm[mask, 6] = 3
                fcf[mask, 6] = 0
    
            
        # convert vertically interpolated u,v into ws, wd
        mask =   ~np.isnan(fobs[:, 7]+   fobs[:, 8]) & np.isnan(foo[:, 9]+foo[:, 10])
        if np.any(mask):
            fobs[mask, 10], fobs[mask, 9] = wswd(fobs[mask, 7] , fobs[mask, 8])
            fcm[mask, 9:] = 2
            fcf[mask, 9:] = 0
            
            
        mask = ~np.isnan(foo)
        if np.all(np.isnan(fobs)):
            for j in range(len(plist)):
                if np.any(mask[j, 7:]):
                    mask[j, 7:] = True
                    
        good =np.any(mask) or np.any(~np.isnan(fobs))    
        for j in range(len(plist)):
            for jj in range(11):

                if not good or mask[j, jj] or (fobs[j, jj] == fobs[j, jj]): # not good copy all-nan records to keep consistency with record number
                    a_obstab['observation_value'][ak] = fobs[j, jj]
                    
                    a_obstab['z_coordinate'][ak] = plist[j]
                    a_obstab['z_coordinate_type'][ak] = 1
                    a_obstab['observed_variable'][ak] = allvar[jj]        
                    a_obstab['date_time'][ak] = rts
                    a_obstab['conversion_flag'][ak] = fcf[j, jj]
                    a_obstab['conversion_method'][ak] = fcm[j, jj]
                    a_loaded_feedback['fg_depar@body'][ak] = ffg_depar[j, jj]
                    a_loaded_feedback['an_depar@body'][ak] = fan_depar[j, jj]
                    a_loaded_feedback['biascorr@body'][ak] = fbc[j, jj]
                    #a_loaded_feedback['biascorr_fg@body'][ak] = fbc_fg[j, jj]
                    orig[ak] = fidx[j, jj]
                    ak += 1
            
        k = i
        if i < oovs:
            
            rts=obstab['date_time'][i]
            
    # add units
    
    ucodes = np.zeros(141, dtype=np.int32)
    ucodes[34] = 5
    ucodes[39] = 622
    ucodes[126] = 5
    ucodes[137] = 5
    ucodes[138] = 0
    ucodes[117] = 806
    ucodes[106] = 110
    ucodes[107] = 731
    ucodes[139] = 731
    ucodes[140] = 731
    
#    for i in range(ak):
    a_obstab['units'][:ak] = ucodes[a_obstab['observed_variable'][:ak]]

    #test code, do not delete
    #k = 0
    #for jj in range(a_obstab['observation_value'].shape[0]):#[:20]:
        
        #if orig[j] != -1:
            #if a_obstab['observation_value'][j] != obstab['observation_value'][orig[j]]:
                #print(j, a_obstab['date_time'][j], obstab['date_time'][orig[j]],
                    #a_obstab['z_coordinate'][j], obstab['z_coordinate'][orig[j]], 
                    #a_obstab['observation_value'][j], obstab['observation_value'][orig[j]])
            #else:
                #k += 1
    #if k != obstab['observation_value'].shape[0]:
        #print('too few values')

    
        
    return a_obstab[:ak], a_loaded_feedback[:ak], ak, orig[:ak],tddmax #addedvar[:l+1]

@njit(parallel=True)
def pfill(arr, v, n):
    if arr.shape[0] < 1000000:
        for i in range(arr.shape[0]):
            arr[i] = v
    else:
      
        chunk = arr.shape[0] // n
        idx = np.zeros(n+1, dtype=np.int64)
        for i in range(1, n):
            idx[i] = idx[i-1] + chunk
        idx[n] = arr.shape[0]
        for i in prange(n):
            for j in range(idx[i], idx[i+1]):
                arr[j] = v

@njit(parallel=True)
def pfilla(arr, v, n):
    if arr.shape[0] < 1000000:
        for i in range(arr.shape[0]):
            arr[i] = v[i]
    else:       
        chunk = arr.shape[0] // n
        idx = np.zeros(n+1, dtype=np.int64)
        for i in range(1, n):
            idx[i] = idx[i-1] + chunk
        idx[n] = arr.shape[0]
        for i in prange(n):
            for j in range(idx[i], idx[i+1]):
                arr[j] = v[j]



@njit(boundscheck=True, cache=True)
def nancheck(darr, farr, n, inan=-2147483648):
    if darr.shape[0] < 1000000:
        for i in range(farr.shape[0]):
            if farr[i] == farr[i]:
                darr[i] = farr[i]
            else:
                darr[i] = inan
    else:
        
        chunk = farr.shape[0] // n
        idx = np.zeros(n+1, dtype=np.int64)
        for i in range(1, n):
            idx[i] = idx[i-1] + chunk
        idx[n] = farr.shape[0]
        for i in prange(n):
            for j in range(idx[i], idx[i+1]):
                if farr[j] == farr[j]:
                    darr[j] = farr[j]
                else:
                    darr[j] = inan

        #for i in prange(farr.shape[0]):
            #if farr[i] == farr[i]:
                #darr[i] = farr[i]
            #else:
                #darr[i] = inan
    
@njit(boundscheck=False, cache=True)
def check_spurious(ov, ot, wind=False):
    
    elsum = 0
    olim = np.zeros((141, 2))
    olim[138, :] = np.array((0, 1.01))
    olim[117, :] = np.array((-10000., 500000.))
    olim[126, :] = np.array((150., 360.))
    olim[137, :] = np.array((120., 360.))
    olim[34, :] = np.array((-0.2, 70.))
    olim[39, :] = np.array((0., 0.1))
    if wind:
        olim[106, :] = np.array((0, 360.))
        olim[107, :] = np.array((0, 200.))
        olim[139, :] = np.array((-200, 200.))
        olim[140, :] = np.array((-200, 200.))
    
    rhsum = 0
    for i in range(ov.shape[0]):
        o = olim[ot[i], 1]
        if o > 0:
            
            if ov[i] < olim[ot[i], 0] or ov[i] > o:
                if ot[i] == 34:
                    x = 0
                if ot[i] == 138 :
                    if ov[i] > 1.1:
                        rhsum += 1
                        if ov[i] <= 101:
                            
                            ov[i] /= 100.
                        else:
                            ov[i] = np.nan
                    elif ov[i] < 0:
                        ov[i]=np.nan
                        #if rhsum > 0.01*ov.shape[0]:
                            #raise ValueError('RH in percent?')
                else:                   
                    ov[i] = np.nan
                elsum += 1
        #if rhsum > 0.01*ov.shape[0]:
            #raise ValueError(f'{rhsum} values above 1.1 - RH in percent?')
        else:
            pass  # no value check implemented
            
    return elsum, rhsum
    
    #idx = np.where(((ov<=0) | (ov>1.01)) &(ot==138) ) # rh
    #elsum += len(idx[0])
    #ov[idx] = np.nan
    #idx = np.where(((ov<0) | (ov>200.)) &(ot==107) ) # ws
    #elsum += len(idx[0])
    #ov[idx] = np.nan
    #idx = np.where(((ov<150) | (ov>360.)) &(ot==126) ) # temp
    #elsum += len(idx[0])
    #ov[idx] = np.nan
    #idx = np.where(((ov<120) | (ov>360.)) &(ot==137) ) # td
    #elsum += len(idx[0])
    #ov[idx] = np.nan
    ##print(len(idx[0]), loaded_obstab['observation_value'][idx])
    #idx = np.where(((ov<-0.2) | (ov>70.)) &(ot==34) ) # DPD
    #elsum += len(idx[0])
    #ov[idx] = np.nan
    #idx = np.where(((ov<0) | (ov>0.1)) &(ot==39) ) # q
    #elsum += len(idx[0])
    #ov[idx] = np.nan

################################################################
        
if __name__ == '__main__':



    ipar=np.zeros(140,dtype=np.int32)-2100000000 # 
    ipar[0]=0
    ipar[34]=34
    ipar[39]=39
    ipar[85]=126
    ipar[106]=106
    ipar[107]=107
    ipar[117]=117
    #ipar[]=136
    ipar[36]=137 #dp
    ipar[38]=138 #rh
    ipar[104]=139
    ipar[105]=140


