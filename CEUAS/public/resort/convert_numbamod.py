#!/usr/bin/env
# coding: utf-8

from numba import njit
import numpy as np
import sys,glob
import zipfile, os, time
import urllib3
from datetime import datetime, timedelta
import glob
sys.path.append(os.getcwd()+'/../cds-backend/code/')
sys.path.append(os.getcwd()+'/../harvest/code/')
sys.path.insert(0,os.getcwd()+'/../resort/rasotools-master/')
import rasotools
dir(rasotools)
import warnings
from functools import partial
import matplotlib.pylab as plt
warnings.filterwarnings('ignore')
import pickle
import ray


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
def liquid(t):
    return -6096.9385 / t + 21.2409642 \
           - 2.711193e-2 * t \
           + 1.673952e-5 * t * t \
           + 2.433502 * np.log(t)

@njit(cache=True)
def ice(t):
    return -6024.5282 / t + 29.32707 \
           + 1.0613868e-2 * t \
           - 1.3198825e-5 * t * t \
           - 0.49382577 * np.log(t)

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
def lvint(p, v, idx):
    lp = np.log(p)
#    for j in range(v.shape[1]):
        
    for k in range(idx.shape[0]):
        j = idx[k]
        if np.sum(~np.isnan(v[:, j])) > 1:
            
            for i in range(1, v.shape[0]-1):
                if v[i, j] != v[i, j]:
                    km =1
                    while i - km > -1 and v[i-km, j] != v[i-km, j]:
                        km += 1
                    kp =1
                    while i + kp < v.shape[0] and v[i+kp, j] != v[i+kp, j]:
                        kp += 1
                    if i + kp < v.shape[0] and  i - km > -1 :
                        if lp[i+kp] - lp[i-km] > 0.:
                            
                            w0 = (lp[i] - lp[i-km]) / (lp[i+kp] - lp[i-km])
                        elif  lp[i+kp] - lp[i-km] < 0:
                            raise ValueError ('lvint: pressure not monotonic')
                        else:
                            w0 = 0.5
                        v[i, j] = (1 - w0) * v[i-km, j] + w0 * v[i+kp, j]
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
    

#@njit(cache=True, boundscheck=False)
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
                ii = np.searchsorted(plist, porig)
                if ii == plist.shape[0]:
                    ii -= 1
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


@njit(cache=True, boundscheck=False)
def fill_obsid(avar,conversion_flag):
    for o in range(avar.shape[0]):
        if conversion_flag[o] == 0:
            for i in range(2):
                avar[o,i]=b'9'
    return avar

@njit(cache=True, boundscheck=False)
def fill_restdata(final, rest_data, addedvar): #,test):

    
    for l in range(addedvar.shape[0]):
        if addedvar[l] != -1:
            final[l] = rest_data[addedvar[l]]

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
def polinl(XA,YA,X):


#    REAL XA(4),YA(NFELD,4),X,Y(NFELD)
#    REAL XX1,XX2,XX3,XX4,X1X2

    Y = np.zeros(YA.shape[0])
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

@njit(cache=True, boundscheck=False)
def triint(tera5, reatab, tlon, tlat, lons, lats, secs, ts,obss,press, zs, p, dy, fn):
    sh = tera5.shape
    tmax = np.max(tera5)
    tmin = np.min(tera5)
    pval = np.zeros((sh[2],sh[3]))
    yval = np.zeros((1, 4))
    
    idx = np.zeros(len(reatab), dtype=np.int32)
    oldlats = -100.
    oldlons = -400.
    oldts = -1
    im = 0
    for il in range(len(reatab)):
        if obss[il] == p:
            found = False
            for ip in range(sh[3]):
                if press[ip] == zs[il]:
                    found = True
                    break
            if not found:
                continue
            if ts[il] != oldts:
                
                its =np.searchsorted(secs, ts[il]) - 1
                if its == -1:
                    its = 0
                oldts = ts[il]
            
            #if its + 1 >= tera5.shape[2]:
                #continue
            w0=(secs[its+1]-ts[il])/(secs[its+1] - secs[its])
            w1=1.-w0
            if lats[il] != oldlats: 
                iyref = np.searchsorted(-tlat, -lats[il]) - 2
                #if iyref == 2:
                    #print(iyref)
                oldlats = lats[il]
            if lons[il] != oldlons:
                if tlon[-1] - tlon[0] > 0:
                    ixref = np.searchsorted(tlon, lons[il]) - 2
                    hlon = tlon
                else:
                    hlon = tlon.copy()
                    delta = tlon[2] - tlon[1]
                    if delta < 0:
                        delta = tlon[3] - tlon[2]
                    hlon[1] = hlon[2] - delta
                    hlon[0] = hlon[1] - delta
                    hlon[-2] = hlon[-3] + delta
                    hlon[-1] = hlon[-2] + delta
                    if np.any(hlon[1:]-hlon[:-1]<0):
                        #print(' not monotonic - fail', fn)
                        raise ValueError
                        
                    if lons[il] < hlon[0]:
                        hlon -= 360.
                        if lons[il] < hlon[0]:
                            #print('fail lons>', lons[il] , hlon[0], fn)
                            raise ValueError
                    if lons[il] > hlon[-1]:
                        hlon += 360.
                        if lons[il] > hlon[-1]:
                            #print('fail lons<', lons[il] , hlon[-1], fn)
                            raise ValueError
                    
                    #if lons[il] > hlon[0]:
                        
                        #hlon[hlon<hlon[0]] += 360.
                    #else:
                        #hlon[hlon>hlon[-1]] -= 360.
                    ixref = np.searchsorted(hlon, lons[il]) - 2
                    #ixref = np.where(np.abs(hlon-lons[il]) <=0.5)[0][0] - 2
                if ixref <0:
                    ixref = 0
                oldlons = lons[il]
            
            for it in its, its + 1:
                if lats[il] > 90 - 2 * dy:
                    pval[it, ip] = tera5[0, 0, it, ip]          
                elif lats[il] < -90 + 2 * dy:
                    pval[it, ip] = tera5[-1, 0, it, ip]          
                else:
                    i = 0
                    for iy in range(4):
                        #print(yval.dtype, tlon.dtype, tera5[iy, :, it, ip:ip+1].T.dtype, lons[il])
                        #print(yval.shape, tlon.shape, tera5[iy, :, it, ip:ip+1].T.shape)
                        #if ixref == 2:
                            #print(ixref)
                        if ixref + 4 > tera5.shape[1]:
                            xx = polinl(hlon[0:0+4],tera5[iyref+iy, 0:0+4, it, ip:ip+1].T,lons[il])
                            #print('ixref was too large!')
                        else:
                            xx = polinl(hlon[ixref:ixref+4],tera5[iyref+iy, ixref:ixref+4, it, ip:ip+1].T,lons[il])
                        #print(xx.shape, xx.dtype)
                        #if xx[0] > tmax or xx[0] < tmin:
                            #print('fail',xx[0], tmax,tmin )
                        yval[i, iy] = xx[0] #polinl(tlon,tera5[iy, :, it, ip:ip+1].T,lons[il])
                    xx = polinl(tlat[iyref:iyref+4], yval, lats[il])
                    #if xx[0] > tmax or xx[0] < tmin:
                        #print('fail',xx[0], tmax,tmin )
                    pval[it, ip] = xx[0] #polinl(tlat, yval, lats[il])          
    
            reatab[im]=pval[its,ip]*w0+pval[its+1,ip]*w1
            idx[im] = il
            im += 1
    
    return idx[:im]


    @njit(cache=True, boundscheck=False)
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


