# System information
import pyrttov
import numpy as np
import os, sys, glob
import xarray
print(sys.executable)
print(sys.version_info)
import pandas as pd
import pandas
pd.set_option('display.max_columns', None)
# pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)
sys.path.append(os.getcwd()+'/../cds-backend/code/')
import cds_eua3 as eua
import pickle
import multiprocessing
from functools import partial

sys.path.insert(0,os.getcwd()+'/../resort/rasotools-master/')
import rasotools
import xarray as xr

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
       
def rttov_calc(tadata, humdata, pressdata, eradata, datedata, chan):
    rttov_installdir = '/rttov/'

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
    
    myProfiles.GasUnits = 1 # kg/kg
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
        if np.array(eradata[i].sp).size > 1:
            S2m.append([float(eradata[i].sp[0])/100., float(eradata[i].t2m[0]), dp_sh(float(eradata[i].d2m[0]), float(eradata[i].sp[0])), float(eradata[i].u10[0]), float(eradata[i].v10[0]), 100000])
            Skin.append([float(eradata[i].skt[0]), 0, 0, 0, 3.0, 5., 15, 0.1, 0.3, 0])
        else:
            S2m.append([float(eradata[i].sp)/100., float(eradata[i].t2m), dp_sh(float(eradata[i].d2m), float(eradata[i].sp)), float(eradata[i].u10), float(eradata[i].v10), 100000])
            Skin.append([float(eradata[i].skt), 0, 0, 0, 3.0, 5., 15, 0.1, 0.3, 0])
        SurfGeom.append([float(eradata[i].latitude), float(eradata[i].longitude), 0.])
        DateTimes.append([datedata[i].year, datedata[i].month, datedata[i].day, 0, 0, 0])
        
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
    print(os.getcwd())
#     msuRttov.FileCoef = "/rttov/rtcoef_rttov12/rttov7pred54L/rtcoef_noaa_17_amsub.dat"
#     msuRttov.FileCoef = "/rttov/rtcoef_rttov12/rttov7pred54L/rtcoef_noaa_16_amsub.dat"
#     msuRttov.FileCoef = "/rttov/rtcoef_rttov12/rttov7pred54L/rtcoef_noaa_15_amsub.dat"
    msuRttov.FileCoef = "/rttov/rtcoef_rttov12/rttov7pred54L/rtcoef_noaa_19_mhs.dat"
    
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
    try:
        # Do not supply a channel list for SEVIRI: this returns emissivity/BRDF values for all
        # *loaded* channels which is what is required
#         surfemisrefl_seviri[0,:,:] = irAtlas.getEmisBrdf(seviriRttov)
#         surfemisrefl_seviri[1,:,:] = brdfAtlas.getEmisBrdf(seviriRttov)
#         surfemisrefl_hirs[0,:,:] = irAtlas.getEmisBrdf(hirsRttov)
#         surfemisrefl_mhs[0,:,:] = mwAtlas.getEmisBrdf(mhsRttov)

        surfemisrefl_msu[0,:,:] = mwAtlas.getEmisBrdf(msuRttov)


    except pyrttov.RttovError as e:
        # If there was an error the emissivities/BRDFs will not have been modified so it
        # is OK to continue and call RTTOV with calcemis/calcrefl set to TRUE everywhere
        sys.stderr.write("Error calling atlas: {!s}".format(e))

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
    print(os.getcwd())
    print

#     print("MSU visible channel reflectances, channels 2-4")
#     for p in range(nprofiles):
#         print("Profile {:d}:".format(p))
#         for c in range(len(chan)):
#             print("  Ch #{:02d} refl={:f}".format(chan_list_msu[c],
#                                                   msuRttov.BtRefl[p, c]))
#         print
    return chan_list_msu, msuRttov.BtRefl


def calc_station(statid, chum, adj = None):
    print(statid)
    statlist = glob.glob('/mnt/users/scratch/leo/scratch/converted_v7/*' + statid + '*_CEUAS_merged_v1.nc')

#     try:
    statlist = glob.glob('./igra_h_single_stats_new/*' + statid + '*.nc')
    ihdf = xarray.open_dataset(statlist[0]).to_dataframe()
    ihdf = ihdf[ihdf.datum > "1979"] 
    ihdf = ihdf.rename({'press':'plev'}, axis='columns')
    ihdf = ihdf[ihdf.plev.isin([3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000])]

    tadata = []
    daydata = []
    datedata = []
    humdata = []
    eradata = []
    chandata = []
    pressdata = []
    mondata = []

    wholemon = []

    for day in [True,False]:
        if day:
            dn_dfta = ihdf.query('hour == 0')
        else:
            dn_dfta = ihdf.query('hour == 1')

        for yr in range(1979,2022,1):
            for mon in range(int(str(yr)+'01'), int(str(yr)+'13'), 1):
                wholemon.append(mon)
                dfta = dn_dfta.loc[(dn_dfta['datum'].dt.year==int(str(mon)[:4])) & (dn_dfta['datum'].dt.month==int(str(mon)[-2:]))]
                df = dfta
                mon_mean = df.groupby(['plev']).aggregate({"ta":np.mean})
                reduced_sh = df.groupby(['plev']).aggregate({"rh":np.mean})
                for j in [3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000]:
                    if len(df[df.plev == j]) < 1:
                        mon_mean = mon_mean[mon_mean.index != j]
                        reduced_sh = reduced_sh[reduced_sh.index != j]

                if len(mon_mean) >= 9:
                    with xarray.open_dataset('./era/era_'+str(yr)+'.nc') as era:
                        era_input = era.sel(time = str(yr)+'-'+str(mon)[-2:]+'-01T00:00:00.000000000', latitude = dfta.lat.iloc[0], longitude = dfta.lon.iloc[0], method='nearest')
                    date = pd.to_datetime(float(era_input.time))
#                     high_mon_mean = mon_mean[mon_mean.index.isin([3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000])]
#                     high_hum_mean = reduced_sh[reduced_sh.index.isin([3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000])]
                    high_mon_mean = mon_mean[mon_mean.index.isin([10000,15000,20000,25000,30000,40000,50000,70000,85000])]
                    high_hum_mean = reduced_sh[reduced_sh.index.isin([10000,15000,20000,25000,30000,40000,50000,70000,85000])]
                    high_hum_mean = high_hum_mean[high_hum_mean.index.isin(high_mon_mean.index)]
                    high_mon_mean = high_mon_mean[high_mon_mean.index.isin(high_hum_mean.index)]

                    if len(high_mon_mean) >= 9:
                        plevs_to_check = high_mon_mean.index
#                         if((3000 in plevs_to_check and 5000 in plevs_to_check and 7000 in plevs_to_check and 10000 in plevs_to_check and 
                        if((10000 in plevs_to_check and 
                                15000 in plevs_to_check and 20000 in plevs_to_check and 25000 in plevs_to_check and 
                                30000 in plevs_to_check and 40000 in plevs_to_check and 50000 in plevs_to_check)):
                            if(not (70000 in plevs_to_check)):
                                s = high_mon_mean.iloc[0]
                                s.ta = np.nan
                                s.name = 70000
                                high_hum_mean = high_hum_mean.append(s)
                                high_mon_mean = high_mon_mean.append(s)
                            if(not (85000 in plevs_to_check)):
                                s = high_mon_mean.iloc[0]
                                s.ta = np.nan
                                s.name = 85000
                                high_hum_mean = high_hum_mean.append(s)
                                high_mon_mean = high_mon_mean.append(s)

                            tadata.append(high_mon_mean.ta)
#                                 humdata.append(high_hum_mean.hus)
                            inter_hum = np.array(rasotools.met.convert.to_sh(temp=xr.DataArray(high_mon_mean.ta), press=xr.DataArray(high_mon_mean.index), 
                                                                    rel_humi=xr.DataArray(high_hum_mean.rh/100.)).values)
                            print(inter_hum)
                            inter_hum[inter_hum <= 0] = np.nan
                            humdata.append(inter_hum)
                            pressdata.append(high_mon_mean.index)
                            datedata.append(date)
                            eradata.append(era_input)  
                            chandata.append(2)
                            daydata.append(day)
                            mondata.append(mon)
#                     else:
#                         print(mon)
#                         print(high_hum_mean)

    a2,b2 = rttov_calc(np.array(tadata)[np.array(chandata) == 2], np.array(humdata)[np.array(chandata) == 2], np.array(pressdata)[np.array(chandata) == 2], [x for x, y in zip(eradata, chandata) if y == 2], np.array(datedata)[np.array(chandata) == 2], [3,4,5])

    middle_index = len(wholemon)//2    

    testb2 =b2[np.array(daydata)[np.array(chandata) == 2] == True]
    date_out2 = np.array(mondata)[np.array(chandata) == 2][np.array(daydata)[np.array(chandata) == 2] == True]

    b_final = []
    for i in wholemon[:middle_index]:
        b = [[np.nan, np.nan, np.nan]]
        if i in date_out2:
            for j in range(len(b)):
                b[j] = testb2[date_out2 == i][0][j]
            b = testb2[date_out2 == i]
        b_final.append(b)

    try:
        os.makedirs("./rttov_out_hum/"+statid+"/")
    except:
        pass

    if adj == None:
        pickle.dump( b_final, open( "rttov_out_hum/"+statid+"/"+statid+"_day_refl.p", "wb" ) )
        pickle.dump( wholemon[:middle_index], open( "rttov_out_hum/"+statid+"/"+statid+"_day_dates.p", "wb" ) )
    else:
        pickle.dump( b_final, open( "rttov_out_hum/"+statid+"/"+adj+"_"+statid+"_day_refl.p", "wb" ) )
        pickle.dump( wholemon[:middle_index], open( "rttov_out_hum/"+statid+"/"+adj+"_"+statid+"_day_dates.p", "wb" ) )
    print('day done: '+statid)


    testb2 =b2[np.array(daydata)[np.array(chandata) == 2] == False]
    date_out2 = np.array(mondata)[np.array(chandata) == 2][np.array(daydata)[np.array(chandata) == 2] == False]

    b_final = []
    for i in wholemon[middle_index:]:
        b = [[np.nan, np.nan, np.nan]]
        if i in date_out2:
            b = testb2[date_out2 == i]
        b_final.append(b)

    if adj == None:
        pickle.dump( b_final, open( "rttov_out_hum/"+statid+"/"+statid+"_night_refl.p", "wb" ) )
        pickle.dump( wholemon[middle_index:], open( "rttov_out_hum/"+statid+"/"+statid+"_night_dates.p", "wb" ) )
    else:
        pickle.dump( b_final, open( "rttov_out_hum/"+statid+"/"+adj+"_"+statid+"_night_refl.p", "wb" ) )
        pickle.dump( wholemon[middle_index:], open( "rttov_out_hum/"+statid+"/"+adj+"_"+statid+"_night_dates.p", "wb" ) )
    print('night done: '+statid)
#     except:
#         print('nothing to calculate: '+statid)
#         if adj == None:
#             pickle.dump( np.array([[[np.nan, np.nan, np.nan]]]), open( "rttov_out_hum/"+statid+"/"+statid+"_day_refl.p", "wb" ) )
#             pickle.dump( [197901], open( "rttov_out_hum/"+statid+"/"+statid+"_day_dates.p", "wb" ) )
#         else:
#             pickle.dump( np.array([[[np.nan, np.nan, np.nan]]]), open( "rttov_out_hum/"+statid+"/"+adj+"_"+statid+"_day_refl.p", "wb" ) )
#             pickle.dump( [197901], open( "rttov_out_hum/"+statid+"/"+adj+"_"+statid+"_day_dates.p", "wb" ) )
#         print('nothing to calculate: '+statid)
        
#         if adj == None:
#             pickle.dump( np.array([[[np.nan, np.nan, np.nan]]]), open( "rttov_out_hum/"+statid+"/"+statid+"_night_refl.p", "wb" ) )
#             pickle.dump( [197901], open( "rttov_out_hum/"+statid+"/"+statid+"_night_dates.p", "wb" ) )
#         else:
#             pickle.dump( np.array([[[np.nan, np.nan, np.nan]]]), open( "rttov_out_hum/"+statid+"/"+adj+"_"+statid+"_night_refl.p", "wb" ) )
#             pickle.dump( [197901], open( "rttov_out_hum/"+statid+"/"+adj+"_"+statid+"_night_dates.p", "wb" ) )
#         print('nothing to calculate: '+statid)


if __name__ == '__main__': 
    try:
        os.makedirs("./rttov_out_hum/")
    except:
        pass
#     [3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000]
    consthum = np.array([6.46, 5.90, 9.17, 20.30,
                         285.00, 1464.00 , 2475.60, 6631.20,
                         15468.00, 21684.00, 35328.00 , 44220.00]
                   )/2.
    statlist = []
    stats = []
#     for i in ['11035', '10393', '72357', '50527']:
#         stats.append(glob.glob('/mnt/users/scratch/leo/scratch/converted_v7/*'+i+'*_CEUAS_merged_v1.nc')[0])
    stats = glob.glob('/mnt/users/scratch/leo/scratch/converted_v7/*_CEUAS_merged_v1.nc')
    for i in stats:
        statlist.append(i.split('-')[-1][:5])
        
    statlist = []
#     stats = glob.glob('./igra_h_single_stats_new/*')
    stats = glob.glob('./igra_h_single_stats_new/*11035*')
    for i in stats:
        statlist.append(i.split('_')[-1][:-3])
    for i in [None]:
        pool = multiprocessing.Pool(processes=40)
        func=partial(calc_station, chum = consthum, adj = i)
        result_list = list(pool.map(func, statlist))
        print(result_list)
