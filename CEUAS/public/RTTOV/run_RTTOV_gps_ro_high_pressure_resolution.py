# System information
import pyrttov
import numpy as np
import os, sys, glob
import xarray
import xarray as xr
print(sys.executable)
print(sys.version_info)
import pandas as pd
import pandas
sys.path.append(os.getcwd()+'/../cds-backend/code/')
import cds_eua4 as eua
import h5py
import pickle
import multiprocessing
from functools import partial
import time
import netCDF4
from datetime import date, timedelta

def ef(p, t=None, over_water=True, over_ice=False, **kwargs):
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

def HylandWexler(temp, over_water=True, over_ice=False, **kwargs):
    """Hyland and Wexler (1983), also in Wexler and Hyland (1983): Stated ranges 173.16 ?
    T < 273.16 for ice and 273.16 > 473.15 for liquid.

    Used by Vaisala

    ln( ew ) = -5800.2206/t + 1.3914993
               - 0.48640239 * 10**(-1)*t
               + 0.41764768 * 10**(-4)*t**2
               - 0.14452093 * 10**(-7)*t**3
               + 6.5459673*math.log(t)

    ln( ei ) = -5674.5359/t + 6.3925247
               - 0.96778430 * 10**(-2)*t
               + 0.62215701 * 10**(-6)*t**2
               + 0.20747825 * 10**(-8)*t**3
               - 0.94840240 * 10**(-12)*t**4
               + 4.1635019*math.log(t)

    Args:
        temp: air temperature in K
        liquid_only: use only water vapor over liquid water
        ice_only: use only water vapor over ice
        kwargs: dummy

    Returns:
         es : saturation water vapor pressure in Pa
    """

    def liquid(t):
        return np.exp(
            -5800.2206 / t + 1.3914993
            - 0.48640239e-1 * t
            + 0.41764768e-4 * t * t
            - 0.14452093e-7 * t * t * t
            + 6.5459673 * np.log(t)
        )

    def ice(t):
        return np.exp(
            -5674.5359 / t + 6.3925247
            - 0.96778430e-2 * t
            + 0.62215701e-6 * t * t
            + 0.20747825e-8 * t * t * t
            - 0.94840240e-12 * t * t * t * t
            + 4.1635019 * np.log(t)
        )

    if over_water:
        return liquid(temp)
    elif over_ice:
        return ice(temp)
    else:
        return np.where(temp < 273.16, ice(temp), liquid(temp))
    
def svp(t, method='HylandWexler', p=None, **kwargs):
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
    try:
        if callable(method):
            vpfunc = method
        else:
            vpfunc = eval(method)

        if p is not None:
            f = ef(p, **kwargs)
        else:
            f = 1.
        return vpfunc(t, **kwargs) * f
    except:
        import sys
        print("Functions: ", ", ".join([i for i in dir(sys.modules[__name__]) if i[0].upper() == i[0]]))
        
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
    rd = 287.05
    rv = 461.50
    c = rd/rv  # Rd/Rv = 0.6219934994582882
    pa = p - e  # dry air pressure
    return (e * c) / (e * c + pa)

def dp_sh(dp, press):
        vpdata = svp(dp, p=press)
        q = vap2sh(vpdata, press)
        return q

# @ray.remote       
def rttov_calc(tadata, humdata, pressdata, eradata, datedata, chan):
    # rttov_installdir = '/users/staff/leo/rttov13/'
    rttov_installdir = '/jetfs/scratch/uvoggenberger/rttov13/'

    # ------------------------------------------------------------------------
    # Set up the profile data
    # ------------------------------------------------------------------------

    # Declare an instance of Profiles
    nlevels = len(pressdata[0])
    nprofiles = len(tadata)
    myProfiles = pyrttov.Profiles(nprofiles, nlevels)

    # Associate the profiles and other data from example_data.h with myProfiles
    # Note that the simplecloud, clwscheme, icecloud and zeeman data are not mandatory and
    # are omitted here

    def expand2nprofiles(n, nprof):
        # Transform 1D array to a [nprof, nlevels] array
        outp = np.empty((nprof, len(n)), dtype=n.dtype)
        for i in range(nprof):
            outp[i, :] = n[:]
        return outp

#     dfsh = ascent.copy()
#     consthum =  pickle.load( open( "dfsh.p", "rb" ) )
#     pl = [1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000,92500,100000]
#     for i in range(len(pl)):
#         for j in range(len(ascent)):
#             if dfsh.index[j] == pl[i]:
#                 dfsh.ta.iloc[j] = consthum[i]
    
    myProfiles.GasUnits = 0
    myProfiles.P = pressdata # expand2nprofiles(pressdata/100., nprofiles) 
#     print(myProfiles.P)
    myProfiles.T = tadata # expand2nprofiles(tadata, nprofiles) 
#     print(myProfiles.T)
    myProfiles.Q = humdata # expand2nprofiles(humdata, nprofiles) 
#     myProfiles.Q = expand2nprofiles(np.array(dfsh.ta), nprofiles) 

#     print(myProfiles.Q)
    
    
    
    myProfiles.Angles = [[0, 0, 45, 180]] * nprofiles
    # satzen, satazi, sunzen, sunazi
    myProfiles.SurfType = [[0, 0]] * nprofiles
    # skin%surftype
    
    S2m = []
    Skin = []
    SurfGeom = []
    DateTimes = []
        
    for i in range(nprofiles):
        if np.array(eradata[i]['sp']).size > 1:
            S2m.append([float(eradata[i].sp[0])/100., float(eradata[i].t2m[0]), dp_sh(float(eradata[i].d2m[0]), float(eradata[i].sp[0])), float(eradata[i].u10[0]), float(eradata[i].v10[0]), 100000])
            Skin.append([float(eradata[i].skt[0]), 0, 0, 0, 3.0, 5., 15, 0.1, 0.3])
        else:
            S2m.append([float(eradata[i]['sp'])/100., float(eradata[i]['t2m']), dp_sh(float(eradata[i]['d2m']), float(eradata[i]['sp'])), float(eradata[i]['u10']), float(eradata[i]['v10']), 100000])
            Skin.append([float(eradata[i]['skt']), 0, 0, 0, 3.0, 5., 15, 0.1, 0.3])
        SurfGeom.append([float(eradata[i]['latitude']), float(eradata[i]['longitude']), 0.])
        dt=pd.to_datetime(datedata[i])
        DateTimes.append([dt.year, dt.month, dt.day, 0, 0, 0])
        
    myProfiles.S2m = S2m
    # s2m%p, s2m%t, s2m%q, s2m%u, s2m%v, s2m%wfetch
    myProfiles.Skin = Skin
    # (skin%t, skin%salinity, skin%foam_fraction, skin%snow_fraction skin%fastem(1:5)) --> fastem default =  3.0, 5., 15, 0.1, 0.3, 0
    myProfiles.SurfGeom = SurfGeom
    # (latitude, longitude, elevation)
    myProfiles.DateTimes = DateTimes
    # (year, month, day, hour, minute, second)

    
    # ------------------------------------------------------------------------
    # Set up Rttov instances for each instrument
    # ------------------------------------------------------------------------

    # Create three Rttov objects for three instruments
    msuRttov = pyrttov.Rttov()

    nchan_msu = len(chan)
    chan_list_msu = chan

    # Set the options for each Rttov instance:
    # - the path to the coefficient file must always be specified
    # - turn RTTOV interpolation on (because input pressure levels differ from
    #   coefficient file levels)
    # - set the verbose_wrapper flag to true so the wrapper provides more
    #   information
    # - enable solar simulations for SEVIRI
    # - enable CO2 simulations for HIRS (the CO2 profiles are ignored for
    #   the SEVIRI and MHS simulations)
    # - enable the store_trans wrapper option for MHS to provide access to
    #   RTTOV transmission structure
#     print("/rttov/rtcoef_rttov12/rttov7pred54L/rtcoef_noaa_14_msu.dat")
#     msuRttov.FileCoef = '{}/{}'.format(rttov_installdir,
#                                        "rtcoef_rttov12/rttov7pred54L/rtcoef_noaa_14_msu.dat")
    msuRttov.FileCoef = rttov_installdir+"rtcoef_rttov13/rttov13pred54L/rtcoef_noaa_14_msu.dat"
#     msuRttov.FileCoef = "/rttov/rtcoef_rttov12/rttov7pred54L/rtcoef_noaa_15_amsub.dat"
    msuRttov.Options.AddInterp = True
#     msuRttov.Options.AddSolar = True
#     msuRttov.Options.CO2Data = False
    msuRttov.Options.VerboseWrapper = True


    # Load the instruments: for HIRS and MHS do not supply a channel list and
    # so read all channels
#     try:
    msuRttov.loadInst(channels=chan_list_msu)
#     msuRttov.loadInst()
#     except pyrttov.RttovError as e:
#         sys.stderr.write("Error loading instrument(s): {!s}".format(e))
#         sys.exit(1)

    # Associate the profiles with each Rttov instance
    msuRttov.Profiles = myProfiles

    # ------------------------------------------------------------------------
    # Load the emissivity and BRDF atlases
    # ------------------------------------------------------------------------

    # Load the emissivity and BRDF atlases:
    # - load data for the month in the profile data
    # - load the IR emissivity atlas data for multiple instruments so it can be used for SEVIRI and HIRS
    # - SEVIRI is the only VIS/NIR instrument we can use the single-instrument initialisation for the BRDF atlas

#     irAtlas = pyrttov.Atlas()
#     irAtlas.AtlasPath = '{}/{}'.format(rttov_installdir, "emis_data")
#     irAtlas.loadIrEmisAtlas(ex.datetimes[1][0], ang_corr=True) # Include angular correction, but do not initialise for single-instrument

#     brdfAtlas = pyrttov.Atlas()
#     brdfAtlas.AtlasPath = '{}/{}'.format(rttov_installdir, "brdf_data")
#     brdfAtlas.loadBrdfAtlas(ex.datetimes[1][0], msuRttov) # Supply Rttov object to enable single-instrument initialisation
#     brdfAtlas.IncSea = False                                 # Do not use BRDF atlas for sea surface types

    # TELSEM2 atlas does not require an Rttov object to initialise
    mwAtlas = pyrttov.Atlas()
    mwAtlas.AtlasPath = '{}/{}'.format(rttov_installdir, "emis_data")
    mwAtlas.loadMwEmisAtlas(-1)

    # Set up the surface emissivity/reflectance arrays and associate with the Rttov objects
    surfemisrefl_msu = np.zeros((2,nprofiles,nchan_msu), dtype=np.float64)

#     seviriRttov.SurfEmisRefl = surfemisrefl_msu

    # ------------------------------------------------------------------------
    # Call RTTOV
    # ------------------------------------------------------------------------

    # Surface emissivity/reflectance arrays must be initialised *before every call to RTTOV*
    # Negative values will cause RTTOV to supply emissivity/BRDF values (i.e. equivalent to
    # calcemis/calcrefl TRUE - see RTTOV user guide)

    surfemisrefl_msu[:,:,:] = -1.

    # Call emissivity and BRDF atlases
    #try:
        ## Do not supply a channel list for SEVIRI: this returns emissivity/BRDF values for all
        ## *loaded* channels which is what is required
##         surfemisrefl_seviri[0,:,:] = irAtlas.getEmisBrdf(seviriRttov)
##         surfemisrefl_seviri[1,:,:] = brdfAtlas.getEmisBrdf(seviriRttov)
##         surfemisrefl_hirs[0,:,:] = irAtlas.getEmisBrdf(hirsRttov)
##         surfemisrefl_mhs[0,:,:] = mwAtlas.getEmisBrdf(mhsRttov)

        #surfemisrefl_msu[0,:,:] = mwAtlas.getEmisBrdf(msuRttov)


    #except pyrttov.RttovError as e:
        ## If there was an error the emissivities/BRDFs will not have been modified so it
        ## is OK to continue and call RTTOV with calcemis/calcrefl set to TRUE everywhere
        #sys.stderr.write("Error calling atlas: {!s}".format(e))

    # Call the RTTOV direct model for each instrument:
    # no arguments are supplied to runDirect so all loaded channels are
    # simulated
    try:
        msuRttov.runDirect()
    except pyrttov.RttovError as e:
        sys.stderr.write("Error running RTTOV direct model: {!s}".format(e))
        sys.exit(1)

    # ------------------------------------------------------------------------
    # Print out some of the output
    # ------------------------------------------------------------------------

    print
    print("SELECTED OUTPUT")
    print

#     print("MSU visible channel reflectances, channels 2-4")
#     for p in range(nprofiles):
#         print("Profile {:d}:".format(p))
#         for c in range(len(chan)):
#             print("  Ch #{:02d} refl={:f}".format(chan_list_msu[c],
#                                                   msuRttov.BtRefl[p, c]))
#         print
    return chan_list_msu, msuRttov.BtRefl
def monmean(data,days):
    sdays=pd.date_range(start='1900-01-01',end='2023-02-01',freq='MS')
    idays=np.array([(sdays[i]-sdays[0]).days for i in range(len(sdays))])
    montemp=[]
    good=[]
    gdays=[]
    for i in range(len(idays)-1):
        start,stop=np.searchsorted(days,idays[i:i+2]+1)
        if stop>start+1:
            d=data[:,:,start:stop]
            x=np.sum(~np.isnan(d),axis=2).reshape((data.shape[0],data.shape[1],1))
            if np.sum(x)>0:
                
                good.append(x.reshape((data.shape[0],data.shape[1],1)))
                montemp.append(np.nanmean(d,axis=2).reshape((data.shape[0],data.shape[1],1)))
                gdays.append(idays[i])
    
    return np.concatenate(montemp,axis=2),np.concatenate(good,axis=2),np.array(gdays)+1

def do_rt_gpsro(era_input,df_temp,df_press,df_time,df_good,chum,mask):
    
    tt = time.time()

    tadata = []
    daydata = []
    datedata = []
    humdata = []
    eradata = []
    chandata = []
    pressdata = []
    mondata = []

    tadata34 = []
    daydata34 = []
    datedata34 = []
    humdata34 = []
    eradata34 = []
    chandata34 = []
    pressdata34 = []
    mondata34 = []

    wholemon = []

    ngoodmon = []
    dgoodmon = []

    shpsave=-1
    ei34=dict()
    ei=dict()
    for k in era_input.keys():
        if k not in ['mon']:
            ei34[k]=None
            ei[k]=None
    
    for ih in [0,1]:

        tt=time.time()
        l=-1
        for yr in range(1950,2022,1):
            #print(yr,time.time()-tt)
            for imon in range(1,13):
                l+=1
                wholemon.append(yr*100+imon)
                mon=str(yr*100+imon)

                x=np.datetime64('{}-{:0>2}-{:0>2}'.format(yr,imon,1))

                idx=np.searchsorted(df_time,x)
                if idx==df_time.shape[0]:
                    continue
                #print(df_time[idx],x)
                if df_time[idx]!=x:
                    continue
                
                #pressure threshold check
                if np.nanmin(df_press) > 30 or np.nanmax(df_press) < 700:
                    print('pressure out of range!')
                    print(np.nanmin(df_press), np.nanmax(df_press))
                    continue

                if era_input['sp'][l]>0:  

                    # 3000 Pa - 70000 Pa
                    prof=df_temp[ih,:,idx]
                    prof = prof[np.logical_and(df_press >= 30 , df_press <=700)]
                    use_chum = chum[np.logical_and(df_press >= 30 , df_press <=700)]
                    use_df_press = df_press[np.logical_and(df_press >= 30 , df_press <=700)]
                            
                    if not any(np.isnan(prof)):
                        tadata34.append(prof)
                        humdata34.append(use_chum)
                        pressdata34.append(use_df_press)
                        datedata34.append(df_time[idx])
                        for k in ei34.keys():
                            if k in ['latitude','longitude']:
                                ei34[k]=era_input[k]
                            else:    
                                ei34[k]=era_input[k][l] 
                        eradata34.append(ei34)  
                        chandata34.append(34)
                        daydata34.append(ih)
                        mondata34.append(mon)
                        mask[ih,1:,idx]=True

                    # 5000 Pa - 85000 Pa
                    lprof=df_temp[ih,:,idx]
                    lprof = lprof[np.logical_and(df_press >= 50 , df_press <=850)]
                    use_chum = chum[np.logical_and(df_press >= 50 , df_press <=850)]
                    use_df_press = df_press[np.logical_and(df_press >= 50 , df_press <=850)]
                    if not any(np.isnan(lprof)):
                        #low_reduced_sh = reduced_sh[reduced_sh.index.isin(low_mon_mean.plev)]
                        tadata.append(lprof)
                        humdata.append(use_chum)
                        pressdata.append(use_df_press)
                        datedata.append(df_time[idx])
                        for k in ei.keys():
                            if k in ['latitude','longitude']:
                                ei[k]=era_input[k]
                            else:    
                                ei[k]=era_input[k][l] 
                        eradata.append(ei)  
                        chandata.append(2)
                        daydata.append(ih)
                        mondata.append(mon)
                        mask[ih,0,idx]=True

                #else:
                    #print(l,ih,idx,era_input['sp'][l])
                    #print('not found')
        
                else:
                    print('quality check: failed')
                    pass
                    #print(l,ih,idx,df_good[ih,3,idx],era_input['sp'][l])
        # print(ih,len(tadata))
    if len(tadata34)==0 or len(tadata)==0:
        print('nans found, returning')
        return
    print('before 34 c')
    #if statid!='91413':
        #return

    #for ih in 0,1:
        
        #indices[ih]=np.sort(np.array(list(set(indices34[0]+indices2[0]))))
    print(time.time()-tt)

    a34,b34 = rttov_calc(
        np.array(tadata34),
        np.array(humdata34),
        np.array(pressdata34),
        eradata34, # [x for x, y in zip(eradata34, chandata34) if y == 34],
        np.array(datedata34), # [np.array(chandata34) == 34],
        [3,4]
    )

    print('before 12 c')
    a2,b2 = rttov_calc(
        np.array(tadata),
        np.array(humdata),
        np.array(pressdata),
        eradata, # ([x for x, y in zip(eradata, chandata) if y == 2],
        np.array(datedata), # [np.array(chandata) == 2],
#         [2]
        [1,2]
    )

    #middle_index = len(wholemon)//2
    print(time.time()-tt)
    return a2,b2,a34,b34,mask


def do_rt(era_input,df_temp,df_press,df_time,df_good,chum,mask):
    
        tt = time.time()
    
        tadata = []
        daydata = []
        datedata = []
        humdata = []
        eradata = []
        chandata = []
        pressdata = []
        mondata = []

        tadata34 = []
        daydata34 = []
        datedata34 = []
        humdata34 = []
        eradata34 = []
        chandata34 = []
        pressdata34 = []
        mondata34 = []

        wholemon = []

        ngoodmon = []
        dgoodmon = []

        shpsave=-1
        ei34=dict()
        ei=dict()
        for k in era_input.keys():
            if k not in ['mon']:
                ei34[k]=None
                ei[k]=None
        
        for ih in [0,1]:

            tt=time.time()
            l=-1
            for yr in range(1950,2022,1):
                #print(yr,time.time()-tt)
                for imon in range(1,13):
                    l+=1
                    wholemon.append(yr*100+imon)
                    mon=str(yr*100+imon)

                    x=np.datetime64('{}-{:0>2}-{:0>2}'.format(yr,imon,1))

                    idx=np.searchsorted(df_time,x)
                    if idx==df_time.shape[0]:
                        continue
                    #print(df_time[idx],x)
                    if df_time[idx]!=x:
                        continue

                    if df_good[ih,3,idx] >= 3 and era_input['sp'][l]>0:  

                        prof=df_temp[ih,2:13,idx]
                                
                        if not any(np.isnan(prof)):
                            tadata34.append(prof)
                            humdata34.append(chum[0:11])
                            pressdata34.append(df_press[2:13])
                            datedata34.append(df_time[idx])
                            for k in ei34.keys():
                                if k in ['latitude','longitude']:
                                    ei34[k]=era_input[k]
                                else:    
                                    ei34[k]=era_input[k][l] 
                            eradata34.append(ei34)  
                            chandata34.append(34)
                            daydata34.append(ih)
                            mondata34.append(mon)
                            mask[ih,1:,idx]=True


                        lprof=df_temp[ih,3:14,idx]
                        if not any(np.isnan(lprof)):
                            #low_reduced_sh = reduced_sh[reduced_sh.index.isin(low_mon_mean.plev)]
                            tadata.append(lprof)
                            humdata.append(chum[1:])
                            pressdata.append(df_press[3:14])
                            datedata.append(df_time[idx])
                            for k in ei.keys():
                                if k in ['latitude','longitude']:
                                    ei[k]=era_input[k]
                                else:    
                                    ei[k]=era_input[k][l] 
                            eradata.append(ei)  
                            chandata.append(2)
                            daydata.append(ih)
                            mondata.append(mon)
                            mask[ih,0,idx]=True

                    #else:
                        #print(l,ih,idx,era_input['sp'][l])
                        #print('not found')
            
                    else:
                        print('quality check: failed')
                        pass
                        #print(l,ih,idx,df_good[ih,3,idx],era_input['sp'][l])
            # print(ih,len(tadata))
        if len(tadata34)==0 or len(tadata)==0:
            print('nans found, returning')
            return
        print('before 34 c')
        #if statid!='91413':
            #return

        #for ih in 0,1:
            
            #indices[ih]=np.sort(np.array(list(set(indices34[0]+indices2[0]))))
        print(time.time()-tt)

        a34,b34 = rttov_calc(
            np.array(tadata34),
            np.array(humdata34),
            np.array(pressdata34),
            eradata34, # [x for x, y in zip(eradata34, chandata34) if y == 34],
            np.array(datedata34), # [np.array(chandata34) == 34],
            [3,4]
        )

        print('before 12 c')
        a2,b2 = rttov_calc(
            np.array(tadata),
            np.array(humdata),
            np.array(pressdata),
            eradata, # ([x for x, y in zip(eradata, chandata) if y == 2],
            np.array(datedata), # [np.array(chandata) == 2],
    #         [2]
            [1,2]
        )

        #middle_index = len(wholemon)//2
        print(time.time()-tt)
        return a2,b2,a34,b34,mask

def calc_gridded(era_input,npzdict, chum, odir, pool,adj = None):

    try:
        tt=time.time()
        tups=[]
        df_press=npzdict['ps'][:]
        df_time=pd.to_datetime(npzdict['days'][:]-1,unit='D',origin='19000101').values
        df_time=df_time[600:600+npzdict['CR20v372'][0,0,600:,0,0].shape[0]]
        for j in range(era_input['t2m'].shape[1]):
            for i in range(era_input['t2m'].shape[2]):
                hilf=npzdict['CR20v372'][0,:,600:,j,i]
                ei=dict()
                for s in era_input.keys():
                    ei[s]=era_input[s][:852,j,i]
                ei['latitude']=-88.75+j*2.5
                ei['longitude']=1.25+i*2.5
                
                tups.append((i,j,np.concatenate((hilf.reshape((1,hilf.shape[0],hilf.shape[1])),hilf.reshape((1,hilf.shape[0],hilf.shape[1]))),axis=0),ei))
        df_good=~np.isnan(tups[0][2])*20
        func=partial(grid_rt,df_press,df_time,df_good,chum)
        print(time.time()-tt)
        result=list(pool.map(func,tups))
        print(time.time()-tt)
        bigarr=np.empty((1,result[0].shape[0],npzdict['CR20v372'].shape[2],era_input['t2m'].shape[1],era_input['t2m'].shape[2]))
        l=-1
        for j in range(era_input['t2m'].shape[1]):
            for i in range(era_input['t2m'].shape[2]):
                l+=1
                bigarr[0,:,600:600+result[0].shape[1]//2,j,i]=result[l][:,:result[0].shape[1]//2]
        print(time.time()-tt)
        
    except MemoryError:
        pass
    
    return bigarr

def grid_rt(df_press,df_time,df_good,chum,tup):
    
    i,j,df_temp,ei=tup        
    print(i, j)
    mask=np.zeros((2,3,df_time.shape[0]),dtype=bool)
    try:
        
        a2,b2,a34,b34 = do_rt(ei,df_temp,df_press,df_time,df_good,chum,mask)
    #def do_rt(era_input,df_temp,df_press,df_time,df_good,chum,mask):

    except MemoryError:
        print(i,j,'failed')
        return

    return np.concatenate((b2.T,b34.T),axis=0)

def calc_station(tup, chum, odir, adj = None, npzdict=None):
#     print(')
#     print(glob.glob('/rttov/rtcoef_rttov12/rttov7pred54L/*ssmt2*'))
    print('entering calc_station')
    tt = time.time()
    statid=tup[0]
    era_input=tup[1]
    statlist = statid
    print('starting calc station')
    if type(statlist)!=np.ndarray:
        
        statid = statlist.split('.nc')[0][-5:]

        try:
            os.makedirs(odir+"/"+statid+"/")
        except Exception as e:
            print('could not create outdir')
            print(str(e))
        print(statid)
    print('before first try')
    try:
        print("first try")
        adjstatlist = None
        mt='montemp'
        if adj == 'bg':
            statlist = 'bg'.join(statlist.split('bincorr'))
        elif adj in ('rharm', 'rharm_h'):
            statlist = (adj+'_').join(statlist.split('feedbackglobbincorrmon'))
            mt='ta' + adj.split('rharm')[1]
        elif adj in ('suny','sunyhom'):
            try:   
                statlist = '/users/staff/leo/fastscratch/SUNY/UA-HRD_stations_homo+raw-monthly-anomalyT_00-12Z_195801-202008.nc'
                with netCDF4.Dataset(statlist,'r') as h:
                    siteids=np.asarray(h['SiteID'][:].T.flatten().view('S11').astype('U11'))
                    slats=np.array(h['lat'])
                    slons=np.array(h['lon'])
                ssiteids=np.array([s.rstrip()[-5:] for s in siteids])
                try:
                    
                    idx=np.where(statid==ssiteids)[0][0]
                except:
                    try:
                        with netCDF4.Dataset(tup[0],'r') as h:
                            lon=np.array(h['lon'])[-1]
                            lat=np.array(h['lat'])[-1]
                        
                        idx=np.where(np.logical_and(np.abs(slats-lat)<2.0,np.abs(slons-lon)<2.0))[0][0]
                    except:
                        print(statid,' no IGRA match found')
                        return
                    
                
                statlist = '/users/staff/leo/fastscratch/SUNY/homo-raw-subdaily-station/'+siteids[idx].strip()+'.nc'
                if not os.path.isfile(statlist):
                    print(siteids[idx],' not found in daily files directory')
                    return
            except:
                print('suny file '+statid+' not found')
                return
            mt='rawT' if adj=='suny' else 'homoT'
        elif adj == 'rio':
            adjstatlist = glob.glob(statlist.split('feedback')[0]+'feedbackglobbincorrsave_rio24_0*.nc')
        elif adj == 'rit':
            adjstatlist = glob.glob(statlist.split('feedback')[0]+'feedbackglobbincorrsave_rit24_0*.nc')
        if adjstatlist != None:
            if len(adjstatlist)>0:
                
                adjstatlist=adjstatlist[0]
                with h5py.File(adjstatlist) as h:
                    adj_time=pd.to_datetime(h['datum'][0,:]-1,unit='D',origin='19000101').values
                    adj_corr=h['rasocorr'][:]
            else:
                print(statlist.split('feedback')[0]+'feedbackglobbincorrsave_rio24_0*.nc'+' : no match found')
                return
                
                
        print(statlist)
        try:
            
            if adj in ('suny','sunyhom'):
                spress=np.array((10.,20,30,50,70,100,150,200,250,300,400,500,700,850,925,1000))
                if 'monthly' in statlist:
                    with netCDF4.Dataset(statlist,'r') as h:
                        station_name=''
                        df_press=h['pres'][:]
                        idx=np.searchsorted(df_press,spress*100)
                        idx[idx==df_press.shape[0]]=df_press.shape[0]-1
                        
                        refdate=datetime.datetime(1900,1,1)
                        siteids=np.asarray(h['SiteID'][:].T.flatten().view('S11').astype('U11'))
                        slats=np.array(h['lat'])
                        slons=np.array(h['lon'])
                        ipx=np.arange(16,dtype='int')
                        ipx[15]=14
                        itxy=(h['time'][:]-190000)//100
                        itym=(h['time'][:]-190000)%100-1
                        itx=np.array(itxy*12+itym)
#                        itxfinal=np.searchsorted(itx,tasks[ishs[0]]['data'].shape[4])
#                        idx=np.zeros(siteids.shape[0],dtype=int)-999
                        T=np.array(h[mt][:])
                        T[T==-9999.]=np.nan
                        T=T+273.15
                        j=0
                        jj=0
                            
                        for ih in range(2):
                            if tasks[ish]['data'].shape[3]>pindex[-1]:
                                tasks[ish]['data'][idx[i],0,ih,pindex,itx[0]:itx[0]+itxfinal]=T[i,:itxfinal,ipx[pindex],ih]
                            else:
                                tasks[ish]['data'][idx[i],0,ih,:len(pindex),itx[0]:itx[0]+itxfinal]=T[i,:itxfinal,ipx[pindex],ih]
                else:
                    with h5py.File(statlist) as h:
                        station_name=''
                        df_press=h['pressure'][:]
                        idx=np.searchsorted(df_press,spress*100)
                        idx[idx==df_press.shape[0]]=df_press.shape[0]-1
                        
                        if df_press[0]>5000:
                            print(statid,' not enough pressure levels')
                            print(df_press)
                            return
                        refdate=datetime.datetime(1900,1,1)
                        hty=h['time'][:]//10000
                        htm=(h['time'][:]%10000)//100
                        htd=h['time'][:]%100
                        df_time=[datetime.datetime(hty[i],htm[i],htd[i]) for i in range(hty.shape[0])]
                        df_days=np.array([(df_time[i]-refdate).days for i in range(hty.shape[0])])+1
                        df_time=pd.to_datetime(df_days-1,unit='D',origin='19000101').values
                        
                        x=np.einsum('kli->ilk', h[mt][:])
                        x[x==-9999.]=np.nan
                        df_temp,df_good,df_days=monmean(x,df_days)
                        df_time=pd.to_datetime(df_days-1,unit='D',origin='19000101').values
                        df_temp=df_temp[:,idx,:]
                        df_temp=df_temp*0.1+273.15
                        df_good=df_good[:,idx,:]
                        for i in range(len(idx)):
                            if not any(df_press==spress[i]*100):
                                df_temp[:,i,:]=np.nan
                                df_good[:,i,:]=0
                        
                        df_press=spress
                        try:
                            
                            df_corr=h['rasocorrmon'][:]
                        except:
                            df_corr=np.zeros_like(df_temp)
                            
                        df_temp[df_temp>400]=np.nan
                    with h5py.File(tup[0]) as o:
                        
                        df_lat=o['lat'][:]
                        df_lon=o['lon'][:]
                    mask=np.zeros((2,3,df_time.shape[0]),dtype=bool)
            else:
                try:
                    
                    with h5py.File(statlist) as h:
                        station_name=h.attrs['Stationname']
                        df_press=h['press'][:]
                        x=h[mt][:]
                        if np.sum(~np.isnan(x))==0:
                            print(statlist,'only nans found')
                            return
                        if b'seconds' in h['datum'].attrs['units']:  
                            df_temp,df_good,df_days=monmean(x,np.asarray(h['datum'][0,:]/86400,dtype='int'))
                            df_time=pd.to_datetime(df_days-1,unit='D',origin='19000101').values
                            df_press/=100
                        else:
                            df_time=pd.to_datetime(h['datum'][0,:]-1,unit='D',origin='19000101').values
                            df_days=h['datum'][0,:]
                            df_temp=h[mt][:]
                            df_good=h['goodmon'][:]
                        try:
                            
                            df_corr=h['rasocorrmon'][:]
                        except:
                            df_corr=np.zeros_like(df_temp)
                        
                        if len(tup)==3:
                            idx=np.searchsorted(npzdict['days'],df_days)
                            idx[idx>=tup[2].shape[0]]=-1
                            df_hilf=df_temp.copy()
                            df_hilf=tup[2][:,:,idx]
                            for i in range(len(idx)):
                                if idx[i]==-1:
                                    df_hilf[:,:,i]=np.nan
                            df_hilf[np.isnan(df_temp)]=np.nan
                            df_temp=df_hilf[:]
                        df_temp[df_temp>400]=np.nan
                        df_lat=h['lat'][:]
                        df_lon=h['lon'][:]
                        mask=np.zeros((2,3,df_time.shape[0]),dtype=bool)
                except FileNotFoundError:
                    print(statlist,'not found')
                    return
                
        except MemoryError:
            print('file not found: ',statlist)
            return


        # #
        # ##
        # ###
        # ####
        # #####
        # if adj in [None, 'rharm','rio', 'rit','suny','sunyhom']:
        #     df_temp += df_corr
        # elif adj in ['raobcore', 'bg', 'rharm_h']:
        #     pass
        # else:
        #     print('not a valid adjustment: ', adj)
        #     return
        # #####
        # ####
        # ###
        # ##
        # #

    #     print('adj_df: ', adj_df)
    #     print('df: ',df)
        debug = False
        if adjstatlist != None:
            if debug:
                df_temp_orig = df_temp.copy()
            for it in range(0, adj_time.shape[0]-1):
                idb=np.searchsorted(df_time,adj_time[it:it+2])
                for ih in [0,1]:
                    for ip in range(df_press.shape[0]):
                        #print('adjusting: ', ip)
                        df_temp[ih,ip,idb[0]:idb[1]]-=adj_corr[ih,ip,it]

            if(debug):
                plt.subplot(1, 2, 1)    
                #plt.plot(df_time, df_temp_orig[0, 3, :])
                #plt.plot(df_time, df_temp[0, 3, :])
                #plt.plot(df_time, df_temp_orig[0, 3, :]-df_temp[0, 3, :])
                for i in range(2, 7):
                    
                    plt.plot(df_time, df_temp_orig[0, i, :]-df_temp[0, i, :])
                    pmask = (df_time >np.datetime64('1978-12-31')) & (df_time < np.datetime64('2007-01-01'))& ( ~np.isnan(df_temp[0, i, :]))
                    z = np.polyfit(np.float64(df_time[pmask])/1.e9/86400/365.25, df_temp[0, i, pmask], 1)
                    print(f'Trend 1979-2006:{i},{z[0] * 10:5.3f} K/10a')
                    if i == 3:
                        
                        plt.title(f'Trend 1979-2006:{z[0] * 10:5.3f} K/10a')
        # print("era_input", era_input)
        print("df_temp", df_temp)
        print("df_press", df_press)
        print(",df_time", df_time)
        print("df_good", df_good)
        print("chum", chum)
        print("mask", mask)
        try:
            a2,b2,a34,b34 = do_rt(era_input,df_temp,df_press,df_time,df_good,chum,mask)
        except Exception as e:
            print(str(e))
            print('nans found',statid,'returning')
            return

        try:
            os.makedirs(odir+"/"+statid+"/")
        except:
            pass

        with netCDF4.Dataset('/jetfs/scratch/uvoggenberger/bt_template.nc') as f:
            if adj in ['rio','rit']:
                fno=adjstatlist[:-9]+'_bt2_'+adjstatlist[-9:]
                fno='_'.join(fno.split('__'))
                fno='mon'.join(fno.split('save'))
            elif adj in ('raobcore','bg','rharm','rharm_h'):
                fno=statlist[:-9]+'_bt2_'+statlist[-9:]
                fno='_'.join(fno.split('__'))
            elif adj in ('suny','sunyhom'):
                fno='/'.join(tup[0].split('/')[:-1])+'/'+adj+'_bt2_'+tup[0][-9:]
            else:
                fno=''.join(statlist.split('corr'))[:-9]+'_bt2_'+statlist[-9:]
            print("fno", fno)
            with netCDF4.Dataset(fno,'w') as g:
                for attname in f.ncattrs():
                    if attname == 'Stationname':
                        try:
                            
                            setattr(g,attname,station_name)
                        except:
                            setattr(g,attname,b'unknown')
                            
                    else:
                        setattr(g,attname,getattr(f,attname))
            
                # To copy the dimension of the netCDF file
            
                for dimname,dim in f.dimensions.items():
                    if dimname == 'time':
                        gmask=np.sum(mask,axis=(0,1))>0
                        xmask=mask[:,:,gmask]
                        g.createDimension(dimname,xmask.shape[2])
                    else:
                        g.createDimension(dimname,len(dim))
            
            
                # To copy the variables of the netCDF file
                for varname,ncvar in f.variables.items():
                    
                    if varname in ['montemp', 'lat', 'lon', 'press', 'datum']:
                        var = g.createVariable(varname,ncvar.dtype,ncvar.dimensions)
                        #Proceed to copy the variable attributes
                        for attname in ncvar.ncattrs():  
                            setattr(var,attname,getattr(ncvar,attname))
                            
                var = g.createVariable('goodmon',int,ncvar.dimensions)
            
                #Finally copy the variable data to the new created variable
                g['lat'][:] = f['lat'][:][-1]
                g['lon'][:] = f['lon'][:][-1]
                g['press'][:] = [2,3,4]
                g['datum'][0,:] = df_days[gmask]
                print('df_days[gmask]', df_days[gmask][-1], df_days[gmask])
                
            
                fillvars = {}
                vars_to_write = ['montemp', 'goodmon']#, 'rasocorrmon', 'eracorrmon', 'ancorrmon']

                hilf=np.empty(g['montemp'].shape,g['montemp'].dtype)
                hilf.fill(np.nan)
                vals=np.sum(mask,axis=2) # needed to split b2 and b34
                #for ih in 0,1:
                    
                hilf[0,0,xmask[0,0,:]]=b2[:vals[0,0],1]
                hilf[0,1,xmask[0,1,:]]=b34[:vals[0,1],0]
                hilf[0,2,xmask[0,2,:]]=b34[:vals[0,2],1] # was [0,1] before -> wrong?
                hilf[1,0,xmask[1,0,:]]=b2[vals[0,0]:,1]
                hilf[1,1,xmask[1,1,:]]=b34[vals[0,1]:,0]
                hilf[1,2,xmask[1,2,:]]=b34[vals[0,2]:,1]
                
                if np.nanmax(hilf)>400. or np.nanmin(hilf)<150.:
                    print('spurious:',np.nanmax(hilf),np.nanmin(hilf)<150.)
                g['montemp'][:]=hilf[:]
                g['goodmon'][:]=0
                for ih in 0,1:
                    g['goodmon'][ih,0,:]=df_good[ih,3,gmask]
                    g['goodmon'][ih,1,:]=df_good[ih,2,gmask]
                    g['goodmon'][ih,2,:]=df_good[ih,2,gmask]
                
                if debug:    
                    plt.subplot(1, 2, 2)
                    mask = (g['datum'][0,:]/365.25 > 78) & (g['datum'][0,:]/365.25 < 107)& ( ~np.isnan(g['montemp'][0, 2, :]))
                    plt.plot(g['datum'][0,:]/365.25, g['montemp'][0, 2, :])
                    z = np.polyfit(g['datum'][0,mask]/365.25, g['montemp'][0, 2, mask], 1)
                    plt.title(f'Trend 1979-2006:{z[0] * 10:5.3f} K/10a')
                    plt.show()
            


        print('done: '+statid,time.time()-tt)
    except MemoryError as e:
        print(e,'nothing to calculate: '+statid)
        return
    return

def read_npz(statlist,shortname):
    
    ipath=os.path.dirname(statlist[0])[:-6]
    rdict={}
    with np.load('/jetfs/scratch/uvoggenberger/allsave.npz') as d:
        for k,v in d.items():
            rdict[k]=v[:]
        rdict['istat']=rdict['lats'].shape[0]
    with np.load('/jetfs/scratch/uvoggenberger/rss_1900_2020.npz') as d:
        for k,v in d.items():
            try:
                
                rdict[k]=v[:]
            except:
                pass

    return rdict
            
    rdict['data']=[np.empty((2,16,rdict['CR20v372'].shape[2]),dtype=np.float32) for i in range(rdict['lats'].shape[0])]
    
    l=-1
    for s in statlist:
        l+=1
        with h5py.File(s,'r') as f:
            
            rdict['data'][l].fill(np.nan)
            lati=int((f['lat'][-1]-90.)/(360/rdict['CR20v372'].shape[3]))
            loni=int((f['lon'][-1]+1.25)/(360/rdict['CR20v372'].shape[4]))
            rdict['lats'][l]=f['lat'][-1]
            rdict['lons'][l]=f['lon'][-1]
            if loni<0:
                loni+=rdict['CR20v372'].shape[4]
            rdict['data'][l][0,:,:]=rdict['CR20v372'][0,:,:,lati,loni].reshape((1,16,rdict['CR20v372'].shape[2]))
            rdict['data'][l][1,:,:]=rdict['data'][l][0,:,:]
            print(l)
        
                    #tasks[k]["msudata"]=d["msudata"].astype(np.float32)
                    #d.close()
    return rdict #lats,lons,days,ps,stnames,istat,data,index

def eragridded(yr):

    fill_values = {
        'u10': -0.4824346,
        'v10': 0.3942312,
        'd2m': 261.6807,
        't2m': 268.16467,
        'skt': 267.5933,
        'sp': 88436.48,
    }

    era_input={}
    tt=time.time()
    print(yr,time.time()-tt)
    with netCDF4.Dataset('/jetfs/scratch/uvoggenberger/era_land_monthly/era_'+str(yr)+'.nc') as h:
        hlat=h['latitude'][:]
        hlon=h['longitude'][:]

        era_all={}    
        #for k in 'u10', 'v10', 'd2m', 't2m', 'skt', 'sp':
            #era_all[k]=np.array(h[k][:])
        for targetlon in range(-180,185,5):
            for targetlat in range(-90,95,5):
                s = str(targetlat).zfill(4) + '_' + str(targetlon).zfill(4)
                era_input[s]={}
                idxlat=np.searchsorted(-hlat, -targetlat)
                if targetlon < 0: 
                    idxlon=np.searchsorted(hlon, targetlon + 360.)
                else:
                    idxlon=np.searchsorted(hlon, targetlon)
                if idxlon==h['longitude'].shape[0]:
                    idxlon=-1
                era_input[s]['latitude']=hlat[idxlat]
                era_input[s]['longitude']=hlon[idxlon]
                for k in 'u10', 'v10', 'd2m', 't2m', 'skt', 'sp':
                    era_input[s][k]=h[k][:,idxlat,idxlon]
                    era_input[s][k] = era_input[s][k].filled(fill_values[k])
                era_input[s]['mon']=yr*100+np.arange(1,13)
        print(yr,time.time()-tt)
    return era_input

def eraseries(statlist,yr):

    era_input={}
    tt=time.time()
    with netCDF4.Dataset('/jetfs/scratch/uvoggenberger/era_land_monthly/era_'+str(yr)+'.nc') as h:
        hlat=h['latitude'][:]
        hlon=h['longitude'][:]

        era_all={}    
        #for k in 'u10', 'v10', 'd2m', 't2m', 'skt', 'sp':
            #era_all[k]=np.array(h[k][:])
        for slong in statlist:
            s=slong.split('/')[-1].split('.')[0][-6:]
            with h5py.File(slong) as f:
                era_input[s]={}
                idxlat=np.searchsorted(-hlat,-f['lat'][0])
                if f['lon'][0]<0: 
                    idxlon=np.searchsorted(hlon,f['lon'][0]+360.)
                else:
                    idxlon=np.searchsorted(hlon,f['lon'][0])
                if idxlon==h['longitude'].shape[0]:
                    idxlon=-1
                era_input[s]['latitude']=hlat[idxlat]
                era_input[s]['longitude']=hlon[idxlon]
                for k in 'u10', 'v10', 'd2m', 't2m', 'skt', 'sp':
                    era_input[s][k]=h[k][:,idxlat,idxlon]
            era_input[s]['mon']=yr*100+np.arange(1,13)
        print(yr,time.time()-tt)
    return era_input

def era1836(yr):

    era_input={}
    hh={}
    tt=time.time()
    with netCDF4.Dataset('./era/era_'+str(yr)+'.nc') as h:
        hlat=np.array(h['latitude'][:])
        hlon=np.array(h['longitude'][:])

        cosj=np.cos(h['latitude'][:]*np.pi/180.)
        for k in 'u10', 'v10', 'd2m', 't2m', 'skt', 'sp':
            era_input[k]=np.zeros((h[k].shape[0],72,144))
            hh=np.flip(np.array(h[k][:]),axis=1)
            hh[hh==h[k].getncattr('missing_value')]=np.nan
            for j in range(len(cosj)):
                hh[:,j,:]*=cosj[j]
            #print(k,np.mean(hh[k])/np.mean(cosj))        

            rfak=h['u10'].shape[2]//era_input['u10'].shape[2]
            cosm=np.zeros(era_input[k].shape[1])
            for j in range(era_input[k].shape[1]):
                cosm[j]=np.mean(cosj[j*rfak+1:(j+1)*rfak+1])
                for i in range(era_input[k].shape[2]):                          
                    era_input[k][:,j,i]=np.mean(hh[:,j*rfak+1:(j+1)*rfak+1,i*rfak:(i+1)*rfak],axis=(1,2))/cosm[j]

            print(k,np.mean(era_input[k]))        

        print(yr,time.time()-tt)
    return era_input

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def calc_coord_pair_all_plev(coord_pair, df_ro, era_input, chum, odir, adj = None, npzdict=None):
    # print('entering calc_coord_pair')
    tt = time.time()
    
    # -> selecting form the dataframe 
    targetlat = float(coord_pair.split('_')[0])
    targetlon = float(coord_pair.split('_')[1])
    # print(targetlat, targetlon)
    try:
        df_target = df_ro.xs(targetlat, level='latitude_bins').xs(targetlon, level='longitude_bins')
    except:
        print('no data found')
        return

    # -> dataframe cleaning and extracting data
    df_target.dropna(subset=['temperature'], inplace=True)
    if len(df_target) < 10:
        return
    # -> selecting only pressure levels < 850 hPa
    df_target = df_target[df_target.pressure <= 85000]
    df_target.sort_values(by='pressure', inplace=True)
    # print(df_target)

    # -> puting data into the correct shape for further processing
    df_temp = [np.reshape((df_target.temperature.values), [-1,1]), np.reshape((df_target.temperature.values), [-1,1])]
    df_press = list(df_target.pressure.values/100)
    df_time = list(df_target.index[0])[:1]
    df_good = [np.reshape([30]*len(list(df_target.temperature.values)), [-1,1]), np.reshape([30]*len(list(df_target.temperature.values)), [-1,1])]
    mask =np.zeros((2,3,1),dtype=bool)
    era_input_target = era_input[coord_pair]
    chum = np.array([3.90, 5.90, 9.17, 20.30,
                        85.00, 1064.00 , 2475.60, 6631.20,
                        15468.00, 21684.00, 35328.00 , 44220.00]
                )/2.
    chum_plevs = [3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000]
    chum_all_plevs = np.interp(df_press, chum_plevs, chum)
    
    # -> RTTOV calculation
    try:
        # print("mask", mask)
        a2,b2,a34,b34,mask = do_rt_gpsro(era_input_target,np.array(df_temp),np.array(df_press),np.array(df_time),np.array(df_good),chum_all_plevs,np.array(mask))
        # print("mask", mask)
    except Exception as e:
        print(str(e))
        print('nans found',coord_pair,'returning')
        return

    # -> sorting and writing data to file
    with netCDF4.Dataset('/jetfs/scratch/uvoggenberger/bt_template.nc') as f:
        fno='/jetfs/scratch/uvoggenberger/rttov_out/gpsro_all_levels/singles/'+'bt2_'+coord_pair+'.nc'
        # print("fno", fno)
        with netCDF4.Dataset(fno,'w') as g:
            for attname in f.ncattrs():
                if attname == 'Stationname':
                    try:
                        
                        setattr(g,attname,station_name)
                    except:
                        setattr(g,attname,b'unknown')
                        
                else:
                    setattr(g,attname,getattr(f,attname))
        
            # To copy the dimension of the netCDF file
        
            for dimname,dim in f.dimensions.items():
                if dimname == 'time':
                    gmask=np.sum(mask,axis=(0,1))>0
                    xmask=mask[:,:,gmask]
                    g.createDimension(dimname,xmask.shape[2])
                else:
                    g.createDimension(dimname,len(dim))
        
        
            # To copy the variables of the netCDF file
            for varname,ncvar in f.variables.items():
                
                if varname in ['montemp', 'lat', 'lon', 'press', 'datum']:
                    var = g.createVariable(varname,ncvar.dtype,ncvar.dimensions)
                    #Proceed to copy the variable attributes
                    for attname in ncvar.ncattrs():  
                        setattr(var,attname,getattr(ncvar,attname))
                        
            var = g.createVariable('goodmon',int,ncvar.dimensions)
        
            # Finally copy the variable data to the new created variable
            g['lat'][:] = targetlat
            g['lon'][:] = targetlon
            g['press'][:] = [2,3,4]
            g['datum'][0,:] = (pd.to_datetime(df_target.index[0][0]) - pd.to_datetime("1900-1-1")).days
            
        
            fillvars = {}
            vars_to_write = ['montemp', 'goodmon']#, 'rasocorrmon', 'eracorrmon', 'ancorrmon']

            hilf=np.empty(g['montemp'].shape,g['montemp'].dtype)
            hilf.fill(np.nan)
            vals=np.sum(mask,axis=2) # needed to split b2 and b34
            
            #for ih in 0,1:
            hilf[0,0,xmask[0,0,:]]=b2[:vals[0,0],1]
            hilf[0,1,xmask[0,1,:]]=b34[:vals[0,1],0]
            hilf[0,2,xmask[0,2,:]]=b34[:vals[0,2],1]
            hilf[1,0,xmask[1,0,:]]=b2[vals[0,0]:,1]
            hilf[1,1,xmask[1,1,:]]=b34[vals[0,1]:,0]
            hilf[1,2,xmask[1,2,:]]=b34[vals[0,2]:,1]
            
            if np.nanmax(hilf)>400. or np.nanmin(hilf)<150.:
                print('spurious:',np.nanmax(hilf),np.nanmin(hilf)<150.)
            g['montemp'][:]=hilf[:]
            g['goodmon'][:]=0
            for ih in 0,1:
                g['goodmon'][ih,0,:]=30 # df_good[ih,3,gmask]
                g['goodmon'][ih,1,:]=30 # df_good[ih,2,gmask]
                g['goodmon'][ih,2,:]=30 # df_good[ih,2,gmask]            


        print('done: '+coord_pair,time.time()-tt)
    # except MemoryError as e:
    #     print(e,'nothing to calculate: '+statid)
        # return
    return

def calc_coord_pair(coord_pair, df_ro, era_input, chum, odir, adj = None, npzdict=None):
    # print('entering calc_coord_pair')
    tt = time.time()
    
    # -> selecting form the dataframe 
    targetlat = float(coord_pair.split('_')[0])
    targetlon = float(coord_pair.split('_')[1])
    # print(targetlat, targetlon)
    try:
        df_target = df_ro.xs(targetlat, level='latitude_bins').xs(targetlon, level='longitude_bins')
    except:
        print('no data found')
        return

    # -> dataframe cleaning and extracting data
    df_target.dropna(subset=['temperature'], inplace=True)
    if len(df_target) < 10:
        return
    target_pressure_list = []
    for tp in 1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000,92500,100000:
        target_pressure_list.append(find_nearest(df_target.pressure, tp))
    df_target = df_target[np.isin(df_target.pressure, target_pressure_list)]
    df_target.sort_values(by='pressure', inplace=True)
    # print(df_target)

    # -> puting data into the correct shape for further processing
    df_temp = [np.reshape((df_target.temperature.values), [-1,1]), np.reshape((df_target.temperature.values), [-1,1])]
    df_press = list(df_target.pressure.values/100)
    df_time = list(df_target.index[0])[:1]
    df_good = [np.reshape([30]*len(list(df_target.temperature.values)), [-1,1]), np.reshape([30]*len(list(df_target.temperature.values)), [-1,1])]
    mask =np.zeros((2,3,1),dtype=bool)
    era_input_target = era_input[coord_pair]
    chum = np.array([3.90, 5.90, 9.17, 20.30,
                        85.00, 1064.00 , 2475.60, 6631.20,
                        15468.00, 21684.00, 35328.00 , 44220.00]
                )/2.
    
    # -> RTTOV calculation
    try:
        # print("mask", mask)
        a2,b2,a34,b34,mask = do_rt(era_input_target,np.array(df_temp),np.array(df_press),np.array(df_time),np.array(df_good),chum,np.array(mask))
        # print("mask", mask)
    except Exception as e:
        print(str(e))
        print('nans found',coord_pair,'returning')
        return

    # -> sorting and writing data to file
    with netCDF4.Dataset('/jetfs/scratch/uvoggenberger/bt_template.nc') as f:
        fno='/jetfs/scratch/uvoggenberger/rttov_out/gpsro_all_levels/singles/'+'bt2_'+coord_pair+'.nc'
        # print("fno", fno)
        with netCDF4.Dataset(fno,'w') as g:
            for attname in f.ncattrs():
                if attname == 'Stationname':
                    try:
                        
                        setattr(g,attname,station_name)
                    except:
                        setattr(g,attname,b'unknown')
                        
                else:
                    setattr(g,attname,getattr(f,attname))
        
            # To copy the dimension of the netCDF file
        
            for dimname,dim in f.dimensions.items():
                if dimname == 'time':
                    gmask=np.sum(mask,axis=(0,1))>0
                    xmask=mask[:,:,gmask]
                    g.createDimension(dimname,xmask.shape[2])
                else:
                    g.createDimension(dimname,len(dim))
        
        
            # To copy the variables of the netCDF file
            for varname,ncvar in f.variables.items():
                
                if varname in ['montemp', 'lat', 'lon', 'press', 'datum']:
                    var = g.createVariable(varname,ncvar.dtype,ncvar.dimensions)
                    #Proceed to copy the variable attributes
                    for attname in ncvar.ncattrs():  
                        setattr(var,attname,getattr(ncvar,attname))
                        
            var = g.createVariable('goodmon',int,ncvar.dimensions)
        
            # Finally copy the variable data to the new created variable
            g['lat'][:] = targetlat
            g['lon'][:] = targetlon
            g['press'][:] = [2,3,4]
            g['datum'][0,:] = (pd.to_datetime(df_target.index[0][0]) - pd.to_datetime("1900-1-1")).days
            
        
            fillvars = {}
            vars_to_write = ['montemp', 'goodmon']#, 'rasocorrmon', 'eracorrmon', 'ancorrmon']

            hilf=np.empty(g['montemp'].shape,g['montemp'].dtype)
            hilf.fill(np.nan)
            vals=np.sum(mask,axis=2) # needed to split b2 and b34
            
            #for ih in 0,1:
            hilf[0,0,xmask[0,0,:]]=b2[:vals[0,0],1]
            hilf[0,1,xmask[0,1,:]]=b34[:vals[0,1],0]
            hilf[0,2,xmask[0,2,:]]=b34[:vals[0,2],1]
            hilf[1,0,xmask[1,0,:]]=b2[vals[0,0]:,1]
            hilf[1,1,xmask[1,1,:]]=b34[vals[0,1]:,0]
            hilf[1,2,xmask[1,2,:]]=b34[vals[0,2]:,1]
            
            if np.nanmax(hilf)>400. or np.nanmin(hilf)<150.:
                print('spurious:',np.nanmax(hilf),np.nanmin(hilf)<150.)
            g['montemp'][:]=hilf[:]
            g['goodmon'][:]=0
            for ih in 0,1:
                g['goodmon'][ih,0,:]=30 # df_good[ih,3,gmask]
                g['goodmon'][ih,1,:]=30 # df_good[ih,2,gmask]
                g['goodmon'][ih,2,:]=30 # df_good[ih,2,gmask]            


        print('done: '+coord_pair,time.time()-tt)
    # except MemoryError as e:
    #     print(e,'nothing to calculate: '+statid)
        # return
    return


###
##
# Run RTTOV for GPSRO data.
##
###

# -> Humidity parameters
# -> [3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000]
consthum = np.array([3.90, 5.90, 9.17, 20.30,
                     85.00, 1064.00 , 2475.60, 6631.20,
                     15468.00, 21684.00, 35328.00 , 44220.00]
               )/2.

# -> select input data
statlist = []
statlist = glob.glob(os.path.expandvars('/jetfs/scratch/uvoggenberger/RO_OPSv5.6.2_L2b_10x10/*.nc'))


i = None

pool = multiprocessing.Pool(processes=60)

tt=time.time()

##
# Create a prepared File with all the meta data from ERA5 for all months and coordinate pairs.
##
try:
    with open('/jetfs/scratch/uvoggenberger/era_gridded.pkl','rb') as f:
        era_input=pickle.load(f)

except:
    func=partial(eragridded)
    result_list = list(pool.map(eragridded,range(1950,2023)))

    era_input={}
    for s in result_list[0].keys():
        era_input[s]={}
        for k in result_list[0][s].keys():
            if k in ['latitude','longitude']:
                era_input[s][k]=result_list[0][s][k]
            else:
                era_input[s][k]=np.concatenate([result_list[i][s][k] for i in range(len(result_list))])

    with open('/jetfs/scratch/uvoggenberger/era_gridded.pkl','wb') as f:
        pickle.dump(era_input,f)

# print(time.time()-tt)


##
# Iterate through all adjustment types.
##

# -> general output directory
odir='/jetfs/scratch/uvoggenberger/rttov_out'
for i in ['gpsro'] :#,'rharm','sunyhom','20CRv3','suny']: #,'20CRv3',[None, 'raobcore', 'rio', 'rit']:

    # -> Create output directory
    if i == None:
        odir += "/rttov_out_unadj_testing"
    elif i == 'raobcore':
        odir += "/rttov_out_raobcore_testing"
    elif i == 'rio':
        odir += "/rttov_out_rio_testing"
    elif i == 'rit':
        odir += "/rttov_out_rit_testing"
    try:
        os.makedirs("./"+odir+"/")
    except:
        pass

    skeys=list(era_input.keys())

    tt=time.time()

    # -> listing all the coordinate pair names, for multiprocessing
    coord_pairs = []
    for targetlon in range(-175,180,5):
        for targetlat in range(-85,90,5):
            coord_pairs.append(str(targetlat).zfill(4) + '_' + str(targetlon).zfill(4))
# ### -> reduced coord pairs for testing!
#             break
#         break

    ro_files = glob.glob('/jetfs/scratch/uvoggenberger/RO_OPSv5.6.2_L2b_10x10/*')
### -> reduced list for testing!
    for ro_file in ro_files[:]:
        fno='/jetfs/scratch/uvoggenberger/rttov_out/gpsro_all_levels/'+'bt2_'+ro_file.split('/')[-1].split('.nc')[0]+'.nc'
        if len(glob.glob(fno)) > 0:
            continue
        df_ro = xr.open_dataset(ro_file).to_dataframe()

        # calc_coord_pair(coord_pairs[0], df_ro = df_ro, era_input = era_input,chum = consthum, adj = i, odir = odir)

        func=partial(calc_coord_pair_all_plev, df_ro = df_ro, era_input = era_input,chum = consthum, adj = i, odir = odir)
        result_list = list(pool.map(func, coord_pairs))

        files = glob.glob('/jetfs/scratch/uvoggenberger/rttov_out/gpsro_all_levels/singles/*.nc')
        vars = ['montemp','goodmon', 'lat', 'lon', 'press', 'datum']
        vars_to_cp = {}
        for var in vars:
            vars_to_cp[var] = []
        for file in files:
            # print(file)
            with netCDF4.Dataset(file) as f:
                for var in vars:
                    vars_to_cp[var].append(f[var][:].data)
        with netCDF4.Dataset(files[0]) as f:
            print("fno", fno)
            with netCDF4.Dataset(fno,'w') as g:
                for attname in f.ncattrs():
                    if attname == 'Stationname':
                        try:
                            setattr(g,attname,station_name)
                        except:
                            setattr(g,attname,b'unknown')
                            
                    else:
                        setattr(g,attname,getattr(f,attname))

                # To copy the dimension of the netCDF file

                for dimname,dim in f.dimensions.items():
                    diml = len(dim)
                    if dimname == 'station':
                        diml = len(vars_to_cp['montemp'])
                    g.createDimension(dimname,diml)

                vars_to_cp
                # To copy the variables of the netCDF file
                for varname,ncvar in f.variables.items():
                    if varname in ['lat', 'lon', 'press', 'datum',]:
                        var = g.createVariable(varname,ncvar.dtype,ncvar.dimensions)
                        #Proceed to copy the variable attributes
                        for attname in ncvar.ncattrs():  
                            setattr(var,attname,getattr(ncvar,attname))

                    if varname in ['montemp', 'goodmon']:
                        var = g.createVariable(varname,ncvar.dtype,('station',) + ncvar.dimensions)
                        #Proceed to copy the variable attributes
                        for attname in ncvar.ncattrs():  
                            setattr(var,attname,getattr(ncvar,attname))

                # Finally copy the variable data to the new created variable
                g['lat'][:] = vars_to_cp['lat']
                g['lon'][:] = vars_to_cp['lon']
                g['press'][:] = vars_to_cp['press'][0]
                g['datum'][:] = vars_to_cp['datum'][0]
                g['montemp'][:] = vars_to_cp['montemp']
                g['goodmon'][:] = vars_to_cp['goodmon']
        for del_f in glob.glob('/jetfs/scratch/uvoggenberger/rttov_out/gpsro_all_levels/singles/*.nc'):
            os.remove(del_f)
