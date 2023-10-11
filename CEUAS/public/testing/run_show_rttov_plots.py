# /home/swd/jupyterhub/alpha/bin/rttov-launcher

# System information
# import pyrttov
import numpy as np
import os, sys, glob
import xarray
print(sys.executable)
print(sys.version_info)
import pandas as pd
# pd.set_option('display.max_columns', None)
# # pd.set_option('display.max_rows', None)
# pd.set_option('display.width', None)
# pd.set_option('display.max_colwidth', -1)

import matplotlib
import matplotlib.pylab as plt
import matplotlib.pyplot as maplt
matplotlib.rcParams.update({'font.size': 20})

import pickle

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
       
        
def rttov_calc(ascent, era, date):
    rttov_installdir = '/rttov/'

    # ------------------------------------------------------------------------
    # Set up the profile data
    # ------------------------------------------------------------------------

    # Declare an instance of Profiles
    nlevels = len(ascent.index)
    nprofiles = 1
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

    myProfiles.GasUnits = 1
    myProfiles.P = expand2nprofiles(np.array(ascent.index)/100., nprofiles) 
#     print(myProfiles.P)
    myProfiles.T = expand2nprofiles(np.array(ascent.ta), nprofiles) 
#     print(myProfiles.T)
    myProfiles.Q = expand2nprofiles(np.array(ascent.hus), nprofiles) 
#     print(myProfiles.Q)
    
    
    
    myProfiles.Angles = [[0, 0, 45, 180]]
    # satzen, satazi, sunzen, sunazi
    myProfiles.S2m = [[float(era.sp)/100., float(era.t2m), dp_sh(float(era.d2m), float(era.sp)), float(era.u10), float(era.v10), 100000]]
#     print(myProfiles.S2m)
    # s2m%p, s2m%t, s2m%q, s2m%u, s2m%v, s2m%wfetch
    myProfiles.Skin = [[float(era.skt), 0, 0, 0, 3.0, 5., 15, 0.1, 0.3, 0]] #float(era.snowc)/100., 3.0, 5., 15, 0.1, 0.3, 0]]
    # (skin%t, skin%salinity, skin%foam_fraction, skin%snow_fraction skin%fastem(1:5)) --> fastem default =  3.0, 5., 15, 0.1, 0.3, 0
    myProfiles.SurfType = [[0, 0]]
    # skin%surftype
    myProfiles.SurfGeom = [[df.lat.iloc[0], df.lon.iloc[0], 0.]]
    # (latitude, longitude, elevation)
    myProfiles.DateTimes = [[date.year, date.month, date.day, 0, 0, 0]]
    # (year, month, day, hour, minute, second)
    
    
    # ------------------------------------------------------------------------
    # Set up Rttov instances for each instrument
    # ------------------------------------------------------------------------

    # Create three Rttov objects for three instruments
    msuRttov = pyrttov.Rttov()

    nchan_msu = 3
    chan_list_msu = [2, 3, 4]

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
    msuRttov.FileCoef = "/rttov/rtcoef_rttov12/rttov7pred54L/rtcoef_noaa_14_msu.dat"
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
    mwAtlas.loadMwEmisAtlas(date.month)

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
    print

    print("MSU visible channel reflectances, channels 2-4")
    for p in range(nprofiles):
        print("Profile {:d}:".format(p))
        for c in range(3):
            print("  Ch #{:02d} refl={:f}".format(chan_list_msu[c],
                                                  msuRttov.BtRefl[p, c]))
        print
    return chan_list_msu, msuRttov.BtRefl




if __name__ == '__main__':

    statid = '11035'
    chan_list = []
    refl = []
    dates = []
    for yr in range(1979,2022,1):
        for mon in range(int(str(yr)+'01'), int(str(yr)+'13'), 1):
            try:
                # add file download once on srvx2
                dfsh = xarray.open_dataset(glob.glob('./rttov_cuon_input/'+statid+'/'+str(mon)+'/dest*spec*')[0]).hus.to_dataframe()
                dfta = xarray.open_dataset(glob.glob('./rttov_cuon_input/'+statid+'/'+str(mon)+'/dest*temp*')[0]).ta.to_dataframe()
                df = dfta.merge(dfsh, how='inner', on=['plev', 'time', 'lat', 'lon'])
                mon_mean = df.groupby(['plev']).aggregate({"ta":np.mean, "hus":np.mean})

                era = xarray.open_dataset('./era/era_'+str(yr)+'.nc')
                era_input = era.sel(time = str(yr)+'-'+str(mon)[-2:]+'-01T00:00:00.000000000', latitude = dfta.lat.iloc[0], longitude = dfta.lon.iloc[0], method='nearest')
                date = pd.to_datetime(float(era_input.time))

                a, b = rttov_calc(mon_mean, era_input, date)
                chan_list.append(a)
                refl.append(b)
                dates.append(mon)
            except:
                chan_list.append(np.nan)
                refl.append(np.nan)
                dates.append(mon)

    pickle.dump( refl, open( "11035_refl.p", "wb" ) )
    pickle.dump( dates, open( "11035_dates.p", "wb" ) )
    pickle.dump( chan_list, open( "11035_chan_list.p", "wb" ) )