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
       
def rttov_calc(ascent, hus, era, date):
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
    myProfiles.Q = expand2nprofiles(np.array(hus.hus), nprofiles) 
#     print(myProfiles.Q)
    
    
    
    myProfiles.Angles = [[0, 0, 45, 180]]
    # satzen, satazi, sunzen, sunazi
    myProfiles.S2m = [[float(era.sp)/100., float(era.t2m), dp_sh(float(era.d2m), float(era.sp)), float(era.u10), float(era.v10), 100000]]
#     print(myProfiles.S2m)
    # s2m%p, s2m%t, s2m%q, s2m%u, s2m%v, s2m%wfetch
    print(float(era.skt))
    myProfiles.Skin = [[float(era.skt), 0, 0, 0, 3.0, 5., 15, 0.1, 0.3, 0]] #float(era.snowc)/100., 3.0, 5., 15, 0.1, 0.3, 0]]
    # (skin%t, skin%salinity, skin%foam_fraction, skin%snow_fraction skin%fastem(1:5)) --> fastem default =  3.0, 5., 15, 0.1, 0.3, 0
    myProfiles.SurfType = [[0, 0]]
    # skin%surftype
    myProfiles.SurfGeom = [[float(era.latitude), float(era.longitude), 0.]]
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

def calc_station(statid, chum):
    print(statid)
    try:
        statlist = glob.glob('/mnt/users/scratch/leo/scratch/converted_v7/*' + statid + '*_CEUAS_merged_v1.nc')

        ###
    #     df = eua.CDMDataset(filename = statlist[0]).to_dataframe(groups=['observations_table'], variables=['observed_variable', 'observation_value', 'date_time', 'z_coordinate', 'latitude', 'longitude'])
        df = eua.CDMDataset(filename = statlist[0]).to_dataframe(groups=['observations_table', 'advanced_homogenisation'], variables=['observed_variable', 'observation_value', 'date_time', 'z_coordinate', 'latitude', 'longitude', 'RASE_bias_estimate'])
        ###

        df = df[df.z_coordinate.isin([1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000,92500,100000])]
        df = df[df.date_time > '1979']
        df = df.rename({'latitude':'lat','longitude':'lon', 'date_time':'time', 'z_coordinate':'plev'}, axis='columns')

        all_dfsh = df[df.observed_variable == 39]
        all_dfsh = all_dfsh.rename({'observation_value':'hus'}, axis='columns')
        dfsh = all_dfsh.groupby(['plev']).aggregate({"hus":np.mean})
        pl = [1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000,92500,100000]
        for i in range(len(pl)):
            for j in range(len(dfsh)):
                if dfsh.index[j] == pl[i]:
                    dfsh.hus.iloc[j] = chum[i]

        all_dfta = df[df.observed_variable == 85]
        all_dfta = all_dfta.rename({'observation_value':'ta'}, axis='columns')

        ##
        all_dfta.ta = all_dfta.ta - all_dfta.RASE_bias_estimate
        ##

        for day in [True,False]:
            if day:
                dn_dfta = all_dfta[all_dfta.time.dt.hour > 6][all_dfta.time.dt.hour <= 18]
            else:
                dn_dfta = all_dfta[all_dfta.time.dt.hour <= 6].append(all_dfta[all_dfta.time.dt.hour > 18]).sort_values(by='time')
            chan_list = []
            refl = []
            dates = []
            for yr in range(1979,2022,1):
                for mon in range(int(str(yr)+'01'), int(str(yr)+'13'), 1):
                    try:
                        dfta = dn_dfta.loc[(dn_dfta['time'].dt.year==int(str(mon)[:4])) & (dn_dfta['time'].dt.month==int(str(mon)[-2:]))]
                        df = dfta
                        mon_mean = df.groupby(['plev']).aggregate({"ta":np.mean})
                        print(mon)
                        reduced_sh = dfsh
                        for j in [1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000,92500,100000]:
                            if len(df[df.plev == j]) < 10:
                                print(j)
                                mon_mean = mon_mean[mon_mean.index != j]
                                reduced_sh = reduced_sh[reduced_sh.index != j]

                        if len(mon_mean) < 9:
                            print('too short')
                            chan_list.append(np.nan)
                            refl.append(np.array([[np.nan, np.nan, np.nan]]))
                            dates.append(mon)
                        else:
                            with xarray.open_dataset('./era/era_'+str(yr)+'.nc') as era:
                                era_input = era.sel(time = str(yr)+'-'+str(mon)[-2:]+'-01T00:00:00.000000000', latitude = dfta.lat.iloc[0], longitude = dfta.lon.iloc[0], method='nearest')
                            date = pd.to_datetime(float(era_input.time))
        #                             print(mon_mean)
        #                             print(reduced_sh)
        #                             print(era_input)
        #                             print(date)
                            a, b = rttov_calc(mon_mean,reduced_sh, era_input, date)
                            chan_list.append(a)

                            plevs_to_check = mon_mean.index
                            if(not (5000 in plevs_to_check and 7000 in plevs_to_check and 10000 in plevs_to_check and 
                                    15000 in plevs_to_check and 20000 in plevs_to_check and 25000 in plevs_to_check and 
                                    30000 in plevs_to_check and 40000 in plevs_to_check and 50000 in plevs_to_check and 
                                    70000 in plevs_to_check and 85000 in plevs_to_check)):
                                b = np.array([[np.nan, b[0][1], b[0][2]]])

                            if(not (3000 in plevs_to_check and 5000 in plevs_to_check and 7000 in plevs_to_check and 
                                    10000 in plevs_to_check and 15000 in plevs_to_check and 20000 in plevs_to_check and 
                                    25000 in plevs_to_check and 30000 in plevs_to_check and 40000 in plevs_to_check and 
                                    50000 in plevs_to_check and 70000 in plevs_to_check)):
                                b = np.array([[b[0][0], np.nan, np.nan]])

                            refl.append(b)
                            dates.append(mon)
                    except:
                        chan_list.append(np.nan)
                        refl.append(np.array([[np.nan, np.nan, np.nan]]))
                        dates.append(mon)

            try:
                os.makedirs("./rttov_out/"+statid+"/")
            except:
                pass
            if day:
                pickle.dump( refl, open( "rttov_out/"+statid+"/RASE_"+statid+"_day_refl.p", "wb" ) )
                pickle.dump( dates, open( "rttov_out/"+statid+"/RASE_"+statid+"_day_dates.p", "wb" ) )
                pickle.dump( chan_list, open( "rttov_out/"+statid+"/RASE_"+statid+"_day_chan_list.p", "wb" ) )
                print('done: '+statid)
            else:
                pickle.dump( refl, open( "rttov_out/"+statid+"/RASE_"+statid+"_night_refl.p", "wb" ) )
                pickle.dump( dates, open( "rttov_out/"+statid+"/RASE_"+statid+"_night_dates.p", "wb" ) )
                pickle.dump( chan_list, open( "rttov_out/"+statid+"/RASE_"+statid+"_night_chan_list.p", "wb" ) )
                print('done: '+statid)
    except:
        print('error: '+statid)


if __name__ == '__main__': 
#     multiprocessing.set_start_method('spawn', force=True)
    consthum =  pickle.load( open( "dfsh.p", "rb" ) )
    statlist = []
    stats = glob.glob('/mnt/users/scratch/leo/scratch/converted_v7/*_CEUAS_merged_v1.nc')
#     stats = glob.glob('/mnt/users/scratch/leo/scratch/converted_v7/*10393*_CEUAS_merged_v1.nc')
    for i in stats:
        statlist.append(i.split('-')[-1][:5])
#     calc_station(statlist[0])
    pool = multiprocessing.Pool(processes=30)
    func=partial(calc_station, chum = consthum)
    result_list = list(pool.map(func, statlist))
    print(result_list)