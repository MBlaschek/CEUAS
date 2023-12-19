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
sys.path.append('/jetfs/home/uvoggenberger/CEUAS/CEUAS/public/cds-backend/code/')
import cds_eua4 as eua
# sys.path.append("/jetfs/home/uvoggenberger/CEUAS/CEUAS/public/resort/rasotools-master")
# from rasotools import met as met
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

def to_sh(vp=None, temp=None, rel_humi=None, dpd=None, press=None, method='HylandWexler', precision=6, **kwargs):
    """
    vp -(p)  sh
    rh (t,p) -> vp (p)-> sh
    dpd (t,p) -> vp (p)-> sh

    Args:
        vp: water vapor pressure
        temp: temperature
        rel_humi: rel. humidity
        dpd: dewpoint departure
        press: air pressure
        method:
        precision:

    Returns:
        spec_humi: specific humidity
    """
    from numpy import around
    from xarray import DataArray
    
    if not isinstance(temp, DataArray):
        raise ValueError("Requires a DataArray", type(temp))

    if rel_humi is None and dpd is None and vp is None:
        raise RuntimeError("Requires either r, dpd or vp for conversion")

    if rel_humi is not None and not isinstance(rel_humi, DataArray):
        raise ValueError("R Requires a DataArray", type(rel_humi))

    if dpd is not None and not isinstance(dpd, DataArray):
        raise ValueError("DPD Requires a DataArray", type(dpd))

    if vp is not None and not isinstance(vp, DataArray):
        raise ValueError("VP Requires a DataArray", type(vp))

    if press is None:
        raise RuntimeError("Conversion ?>Q requires a pressure variable as well")

    if vp is not None:
        qvar = vp.copy()
    else:
        qvar = temp.copy()

    if isinstance(press, str):
        if press in qvar.dims:
            press = qvar[press].values
            press = _conform(press, qvar.values.shape)
    elif isinstance(press, DataArray):
        press = press.values
    else:
        pass

    kwargs['tol'] = kwargs.get('tol', 0.1)  # Dewpoint minimization accuracy
    if vp is not None:
        # VP to Q
        qvar.values = vap2sh(vp.values, press)
        origin = 'vp,p'
    elif rel_humi is not None:
        # RH to Q
        vpdata = rel_humi.values * svp(temp.values, method=method, p=press, **kwargs)
        qvar.values = vap2sh(vpdata, press)
        origin = 't,rh,p'
    else:
        # DPD to Q
        vpdata = svp(temp.values - dpd.values, method=method, p=press, **kwargs)
        qvar.values = vap2sh(vpdata, press)
        origin = 't,dpd,p'

    r_att = {'origin': origin, 'esat': method, 'standard_name': 'specific_humidity',
             'long_name': 'specific humidity', 'units': 'kg/kg'}
    if press is not None:
        r_att['enhancement_factor'] = "yes"

    r_att['precision'] = precision
    qvar.attrs.update(r_att)
    qvar.values = around(qvar.values, decimals=precision)
    qvar.name = 'sh'
    return qvar

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
    """ set to 1 for sh !"""
    myProfiles.GasUnits = 1
    myProfiles.P = pressdata/100. # expand2nprofiles(pressdata/100., nprofiles) 
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
        dt=pd.to_datetime(datedata[i], format='%Y%m')
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

    # msuRttov.FileCoef = rttov_installdir+"rtcoef_rttov13/rttov13pred54L/rtcoef_noaa_14_msu.dat"
    msuRttov.FileCoef = rttov_installdir+"rtcoef_rttov13/rttov13pred54L/rtcoef_noaa_19_mhs.dat"
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
    # try:
    msuRttov.runDirect()
    # except pyrttov.RttovError as e:
    #     sys.stderr.write("Error running RTTOV direct model: {!s}".format(e))
    #     sys.exit(1)

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
            # set correct timerange -> 1970 up to date?
            for yr in range(1980,2022,1):
                #print(yr,time.time()-tt)
                for imon in range(1,13):
                    l+=1
                    wholemon.append(yr*100+imon)
                    mon=str(yr*100+imon)

                    x=str(np.datetime64('{}{:0>2}'.format(yr,imon)))

                    idx=np.searchsorted(df_time,x)
                    if idx==df_time.shape[0]:
                        continue
                    #print(df_time[idx],x)
                    if df_time[idx]!=x:
                        continue

                    if np.logical_and(era_input['sp'][l]>0 , np.all(df_good[ih,-9:-1,idx] >= 5)): # checks if all of the important levels are built means of at least 5 datetimes
                        prof_temp=df_temp[ih,-9:-1,idx]
                        prof_sh = np.array(to_sh(temp=xr.DataArray(prof_temp), press=xr.DataArray(np.array([20000, 25000, 30000, 40000, 50000, 70000, 85000, 92500])), rel_humi=xr.DataArray(chum[ih,-9:-1,idx])).values)
                        prof_sh[np.logical_or((prof_sh < 0.), (prof_sh > 1.03))] = np.nan

                        # # -> [20000,25000,30000,40000,50000,70000,85000]
                        # consthum = np.array([1064.00 , 2475.60, 6631.20,
                        #                     15468.00, 21684.00, 35328.00 , 44220.00, 48220.00]
                        # )/2.
                        # prof_sh = consthum
                                
                        if ((not any(np.isnan(prof_temp))) or (not any(np.isnan(prof_sh)))):
                            tadata34.append(prof_temp)
                            humdata34.append(prof_sh)
                            pressdata34.append(df_press[-9:-1])
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
                            mask[ih,0:,idx]=True


                        # lprof=df_temp[ih,3:14,idx]
                        # if not any(np.isnan(lprof)):
                        #     #low_reduced_sh = reduced_sh[reduced_sh.index.isin(low_mon_mean.plev)]
                        #     tadata.append(lprof)
                        #     humdata.append(chum[1:])
                        #     pressdata.append(df_press[3:14])
                        #     datedata.append(df_time[idx])
                        #     for k in ei.keys():
                        #         if k in ['latitude','longitude']:
                        #             ei[k]=era_input[k]
                        #         else:    
                        #             ei[k]=era_input[k][l] 
                        #     eradata.append(ei)  
                        #     chandata.append(2)
                        #     daydata.append(ih)
                        #     mondata.append(mon)
                        #     mask[ih,0,idx]=True

                    #else:
                        #print(l,ih,idx,era_input['sp'][l])
                        #print('not found')
            
                    else:
                        print('quality check: failed')
                        pass
                        #print(l,ih,idx,df_good[ih,3,idx],era_input['sp'][l])
            # print(ih,len(tadata))
        if len(tadata34)==0: # or len(tadata)==0:
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
            [3,4,5] # [3,4,5]
        )

    #     print('before 12 c')
    #     a2,b2 = rttov_calc(
    #         np.array(tadata),
    #         np.array(humdata),
    #         np.array(pressdata),
    #         eradata, # ([x for x, y in zip(eradata, chandata) if y == 2],
    #         np.array(datedata), # [np.array(chandata) == 2],
    # #         [2]
    #         [1,2]
    #     )

        #middle_index = len(wholemon)//2
        print(time.time()-tt)
        return a34,b34,mask


def eragridded(yr):

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
                era_input[s]['mon']=yr*100+np.arange(1,13)
        print(yr,time.time()-tt)
    return era_input

def eraseries(statlist,yr):
    print('eraseries running for: ', statlist, yr)

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
                try:
                    idxlat=np.searchsorted(-hlat,-f['observations_table']['latitude'][-1])
                    if f['observations_table']['longitude'][-1]<0: 
                        idxlon=np.searchsorted(hlon,f['observations_table']['longitude'][-1]+360.)
                    else:
                        idxlon=np.searchsorted(hlon,f['observations_table']['longitude'][-1])
                    if idxlon==h['longitude'].shape[0]:
                        idxlon=-1
                    era_input[s]['latitude']=hlat[idxlat]
                    era_input[s]['longitude']=hlon[idxlon]
                    for k in 'u10', 'v10', 'd2m', 't2m', 'skt', 'sp':
                        era_input[s][k]=h[k][:,idxlat,idxlon]
                except:
                    for k in 'u10', 'v10', 'd2m', 't2m', 'skt', 'sp':
                        era_input[s][k]=np.nan
                    era_input[s]['latitude']=hlat
                    era_input[s]['longitude']=hlon
            era_input[s]['mon']=yr*100+np.arange(1,13)
        print(yr,time.time()-tt)
    return era_input

def eraseries_station(slong,yr):
    print(slong, yr)
    era_input={}
    tt=time.time()
    with netCDF4.Dataset('/jetfs/scratch/uvoggenberger/era_land_monthly/era_'+str(yr)+'.nc') as h:
        hlat=h['latitude'][:]
        hlon=h['longitude'][:]

        era_all={}    
        #for k in 'u10', 'v10', 'd2m', 't2m', 'skt', 'sp':
            #era_all[k]=np.array(h[k][:])
        s=slong.split('/')[-1].split('.')[0][-6:]
        with h5py.File(slong) as f:
            era_input[s]={}
            try:
                idxlat=np.searchsorted(-hlat,-f['observations_table']['latitude'][-1])
                if f['observations_table']['longitude'][-1]<0: 
                    idxlon=np.searchsorted(hlon,f['observations_table']['longitude'][-1]+360.)
                else:
                    idxlon=np.searchsorted(hlon,f['observations_table']['longitude'][-1])
                if idxlon==h['longitude'].shape[0]:
                    idxlon=-1
                era_input[s]['latitude']=hlat[idxlat]
                era_input[s]['longitude']=hlon[idxlon]
                for k in 'u10', 'v10', 'd2m', 't2m', 'skt', 'sp':
                    era_input[s][k]=h[k][:,idxlat,idxlon]
            except:
                for k in 'u10', 'v10', 'd2m', 't2m', 'skt', 'sp':
                    era_input[s][k]=np.nan
                era_input[s]['latitude']=hlat
                era_input[s]['longitude']=hlon
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

def calc_coord_pair(station, era_input, odir, adj = None, npzdict=None):
    print(station)
    station_id = '0'+station.split('/')[-1].split('_')[0].split('-')[-1]
    if adj != None:
        fno=odir+'bt2_'+station_id+'_'+str(adj)+'.nc'
    else:
        fno=odir+'bt2_'+station_id+'.nc'
    if len(glob.glob(fno)) > 0:
        return
    # print('entering calc_coord_pair')
    tt = time.time()
    try:
        df = eua.CDMDataset(filename = station).to_dataframe(groups=['observations_table', 'advanced_homogenisation'], variables=['observed_variable', 'observation_value', 'date_time', 'z_coordinate', 'latitude', 'longitude','humidity_bias_estimate'])
    except:
        print('missing adv. hom.: ', station_id)
        return 
        ###
    # -> selecting form the dataframe 
    targetlat = float(df.latitude.iloc[-1])
    targetlon = float(df.longitude.iloc[-1])
    # print(targetlat, targetlon)
    df = df[df.date_time >= '1972']
    df = df[df.date_time <= '2022']
    if len(df) < 10:
        print("too little data")
        return
    df = df.rename({'latitude':'lat','longitude':'lon', 'date_time':'time', 'z_coordinate':'plev'}, axis='columns')

    ### maybe needs to be rh -> so humidity bias adjustment can be applied -> convert after with rasotools? YES! Still called sh, but is rh!                                                                                                                                                                                                                                                                            
    all_dfsh = df[df.observed_variable == 138]
    all_dfsh = all_dfsh.rename({'observation_value':'sh'}, axis='columns')

    all_dfta = df[df.observed_variable == 126]
    all_dfta = all_dfta.rename({'observation_value':'temperature'}, axis='columns')
    # all_dfta = all_dfta.drop(columns=['humidity_bias_estimate'])

    df_target = pd.merge(all_dfsh, all_dfta, on=["time", "plev", "lat", "lon"])
    df_target = df_target.rename({'plev':'pressure'}, axis='columns')

    if adj != None:
        try:
            df_target['sh'] = df_target['sh'] - np.nan_to_num(df_target['humidity_bias_estimate_y']) # add _y
        except:
            print('missing adv. hom.: ', station_id)
            return 

    # -> dataframe cleaning and extracting data
    df_target.dropna(subset=['temperature'], inplace=True)
    if len(df_target) < 10:
        print("too little data")
        return
    target_pressure_list = []
    std_plevs = [1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000,92500,100000]
    # for tp in std_plevs:
    #     target_pressure_list.append(find_nearest(df_target.pressure, tp))
    df_target = df_target[np.isin(df_target.pressure, std_plevs)]
    # multiple ascents - do not sorty by pressue -> would mix everything   # df_target.sort_values(by='pressure', inplace=True)

    # df_date = df_target.set_index('time', inplace=False)
    # df_month = df_date.resample('M').mean()
    df_target.index = pd.to_datetime(df_target['time'],format='%y-%m/%d %I:%M%p')
    
    ##
    # right now it does mean(day+night) twice -> can be split up
    ##

    # -> puting data into the correct shape for further processing
    df_target_temp = df_target.groupby(by=[df_target.index.month, df_target.index.year, 'pressure']).aggregate({"temperature":"mean"})
    time_for_df = np.char.add(np.array(df_target_temp.index.get_level_values(1)).astype('str') , np.char.zfill(np.array(df_target_temp.index.get_level_values(0)).astype('str'), 2))
    df_target_temp['time_for_sorting'] = time_for_df
    df_target_temp.sort_values(by=['time_for_sorting','pressure'], inplace=True)

    df_target_sh = df_target.groupby(by=[df_target.index.month, df_target.index.year, 'pressure']).aggregate({"sh":["count","mean"]})
    time_for_df = np.char.add(np.array(df_target_sh.index.get_level_values(1)).astype('str') , np.char.zfill(np.array(df_target_sh.index.get_level_values(0)).astype('str'), 2))
    df_target_sh['time_for_sorting'] = time_for_df
    df_target_sh.sort_values(by=['time_for_sorting','pressure'], inplace=True)

    df_temp = []
    df_sh = []
    df_good = []
    for pi in np.array(std_plevs)[:]:
        df_temp_pi = []
        df_sh_pi = []
        df_good_pi = []
        times = df_target_temp.time_for_sorting.drop_duplicates()
        try:
            xs_temp = df_target_temp.xs(pi, level='pressure')
            xs_sh = df_target_sh.xs(pi, level='pressure')
        
            for time_i in times:
                if time_i in list(xs_temp.time_for_sorting):
                    df_temp_pi.append(xs_temp[xs_temp.time_for_sorting == time_i].temperature.values[0])
                    df_sh_pi.append(xs_sh[xs_sh.time_for_sorting == time_i].sh['mean'].values[0])
                    df_good_pi.append(xs_sh[xs_sh.time_for_sorting == time_i].sh['count'].values[0])
                else:
                    df_temp_pi.append(np.nan)
                    df_sh_pi.append(np.nan)
                    df_good_pi.append(0)
            df_temp.append(df_temp_pi)
            df_sh.append(df_sh_pi)
            df_good.append(df_good_pi)
        except:
            df_temp.append([np.nan]*len(times))
            df_sh.append([np.nan]*len(times))
            df_good.append([0]*len(times))
    df_temp = [df_temp, df_temp]
    df_sh = [df_sh, df_sh]
    df_good = [df_good, df_good]
    # df_press = np.array(df_target_temp.index.get_level_values(2)) 
    df_press = std_plevs
    df_time = df_target_temp.time_for_sorting.drop_duplicates().values
    # df_good = np.full_like(df_temp, 30) # [np.reshape([30]*len(list(df_target.temperature.values)), [-1,1]), np.reshape([30]*len(list(df_target.temperature.values)), [-1,1])]
    mask =np.zeros((2,3,len(df_time)),dtype=bool)
    try:
        era_input_target = era_input[station_id]
    except:
        print("no era data")
        return
    
    
    # chum = df_target_sh.sh.values

    # chum = np.array([3.90, 5.90, 9.17, 20.30,
    #                     85.00, 1064.00 , 2475.60, 6631.20,
    #                     15468.00, 21684.00, 35328.00 , 44220.00]
    #             )/2.
    
    # -> RTTOV calculation
    # try:
    # print("mask", mask)
    try:
        a34,b34,mask = do_rt(era_input_target,np.array(df_temp),np.array(df_press),np.array(df_time),np.array(df_good),np.array(df_sh),np.array(mask))
    except:
        print('failed: ', station_id)
        return 

    # print("mask", mask)
    # except Exception as e:
    #     print(str(e))
    #     print('nans found',station_id,'returning')
    #     return

    # -> sorting and writing data to file
    with netCDF4.Dataset('/jetfs/scratch/uvoggenberger/bt_template.nc') as f:
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
            g['press'][:] = [3,4,5]
            g['datum'][0,:] = (pd.to_datetime(df_time[gmask],format='%Y%m') - pd.to_datetime("1900-1-1")).days
            g['datum'][1,:] = (pd.to_datetime(df_time[gmask],format='%Y%m') - pd.to_datetime("1900-1-1")).days
            g['datum'][2,:] = (pd.to_datetime(df_time[gmask],format='%Y%m') - pd.to_datetime("1900-1-1")).days
            g['datum'][3,:] = (pd.to_datetime(df_time[gmask],format='%Y%m') - pd.to_datetime("1900-1-1")).days
            
        
            fillvars = {}
            vars_to_write = ['montemp', 'goodmon']#, 'rasocorrmon', 'eracorrmon', 'ancorrmon']

            hilf=np.empty(g['montemp'].shape,g['montemp'].dtype)
            hilf.fill(np.nan)
            vals=np.sum(mask,axis=2) # needed to split b2 and b34
            
            #for ih in 0,1:
            # hilf[0,0,xmask[0,0,:]]=b2[:vals[0,0],1]
            hilf[0,0,xmask[0,0,:]]=b34[:vals[0,0],0]
            hilf[0,1,xmask[0,1,:]]=b34[:vals[0,1],1]
            hilf[0,2,xmask[0,2,:]]=b34[:vals[0,2],2]
            hilf[1,0,xmask[1,0,:]]=b34[vals[1,0]:,0]
            hilf[1,1,xmask[1,1,:]]=b34[vals[1,1]:,1]
            hilf[1,2,xmask[1,2,:]]=b34[vals[1,2]:,2]
            
            if np.nanmax(hilf)>400. or np.nanmin(hilf)<150.:
                print('spurious:',np.nanmax(hilf),np.nanmin(hilf)<150.)
            g['montemp'][:]=hilf[:]
            g['goodmon'][:]=0
            for ih in 0,1:
                try:
                    g['goodmon'][ih,0,:]= (np.nanmean(np.array(df_good)[ih,-9:-1,:], axis = 0)[gmask]).astype(int)
                    g['goodmon'][ih,1,:]= (np.nanmean(np.array(df_good)[ih,-9:-1,:], axis = 0)[gmask]).astype(int)
                    g['goodmon'][ih,2,:]= (np.nanmean(np.array(df_good)[ih,-9:-1,:], axis = 0)[gmask]).astype(int)      
                except:
                    g['goodmon'][ih,0,:]= 99999 
                    g['goodmon'][ih,1,:]= 99999 
                    g['goodmon'][ih,2,:]= 99999         


        print('done: '+station_id,time.time()-tt)
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

pool = multiprocessing.Pool(processes=30)

tt=time.time()

##
# Create a prepared File with all the meta data from ERA5 for all months and coordinate pairs.
##
statlist = glob.glob('/mnt/users/scratch/leo/scratch/converted_v13/long/*.nc')
try:
    with open('/jetfs/scratch/uvoggenberger/era_input.pkl','rb') as f:
        era_input=pickle.load(f)
        
except:
    for stat in statlist:
        func=partial(eraseries_station, stat)
        result_list = list(pool.map(func,range(1950,2023)))
        era_input={}
        for s in result_list[0].keys():
            era_input[s]={}
            for k in result_list[0][s].keys():
                if k in ['latitude','longitude']:
                    era_input[s][k]=result_list[0][s][k]
                else:
                    era_input[s][k]=np.concatenate([result_list[i][s][k] for i in range(len(result_list))])
        
        with open('/jetfs/scratch/uvoggenberger/era_stations/era_stations_'+stat.split('/')[-1].split('.')[0][-6:]+'.pkl','wb') as f:
            pickle.dump(era_input,f)
print(time.time()-tt)


##
# Iterate through all adjustment types.
##

# -> general output directory -> set with odir

station = glob.glob('/mnt/users/scratch/leo/scratch/converted_v13/long/*.nc')

# for unadj:

odir = '/jetfs/scratch/uvoggenberger/rttov_out/humidity/' # '/jetfs/scratch/uvoggenberger/rttov_out/humidity_adj/'
func=partial(calc_coord_pair, era_input = era_input, adj = None, odir = odir) # "hbe"
result_list = list(pool.map(func, station[:]))

# check: 080259, 043285
# for stat in station:
#     result_list = calc_coord_pair(stat, era_input = era_input, adj = "hbe", odir = odir)
#     break

# are "good" ok? Day night? Compare to NOAA. 
