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
pd.set_option('display.max_columns', None)
# pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)
sys.path.append(os.getcwd()+'/../cds-backend/code/')
import cds_eua3 as eua
import pickle
import multiprocessing
from functools import partial
import time
import h5py
import netCDF4
# import ray
# import time
# ray.init(num_cpus=40)

# @ray.remote

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
    rttov_installdir = '/users/staff/leo/rttov13/'

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


def calc_station(tup, chum, odir, adj = None):
#     print('')
#     print(glob.glob('/rttov/rtcoef_rttov12/rttov7pred54L/*ssmt2*'))
#     print('')
    statid=tup[0]
    era_input=tup[1]
    statlist = statid
#     adjstatlist = glob.glob(statlist.split('feedbackmerged')[0]+'feedbackglobbincorrsave_rio24_0*.nc')
#     adjstatlist = glob.glob(statlist.split('feedbackmerged')[0]+'feedbackglobbincorrsave_rit24_0*.nc')
#     adjstatlist = glob.glob(statlist.split('feedbackmerged')[0]+'feedbackglobbincorrsave0*.nc')
#     feedbackglobbincorrmon

#     print(statlist.split('feedbackmerged')[0]+'feedbackglobbincorrsave0*.nc', adjstatlist)
#     if len(adjstatlist) == 0:
# #         print('no adj')
#         return
#     else:
#         adjstatlist = adjstatlist[0]
        
    statid = statlist.split('.nc')[0][-5:]
#     if len(glob.glob("./"+odir+"/"+statid)) > 0:
#         print('skipped', glob.glob("./"+odir+"/"+statid))
#         return
    try:
        os.makedirs("./"+odir+"/"+statid+"/")
    except:
        pass
    print(statid)

    try:

        adjstatlist = None
        if adj == 'bg':
            statlist = 'bg'.join(statlist.split('bincorr'))
        elif adj == 'rio':
            adjstatlist = glob.glob(statlist.split('feedback')[0]+'feedbackglobbincorrsave_rio24_0*.nc')
        elif adj == 'rit':
            adjstatlist = glob.glob(statlist.split('feedback')[0]+'feedbackglobbincorrsave_rit24_0*.nc')
        if adjstatlist != None:
            if len(adjstatlist)>0:
                
                adjstatlist=adjstatlist[0]
                with h5py.File(adjstatlist) as h:
                    adj_time=pd.to_datetime(h['datum'][0,:]-1,unit='d',origin='19000101').values
                    adj_corr=h['rasocorr'][:]
            else:
                print(statlist.split('feedback')[0]+'feedbackglobbincorrsave_rio24_0*.nc'+' : no match found')
                return
                
            #adj_df = xr.open_dataset(adjstatlist).to_dataframe()
            #adj_df.press = adj_df.press * 100.
            #adj_df = adj_df[adj_df.press.isin([3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000])]

        print(statlist)
        try:
            
            with h5py.File(statlist) as h:
                station_name=h.attrs['Stationname']
                df_time=pd.to_datetime(h['datum'][0,:]-1,unit='d',origin='19000101').values
                df_days=h['datum'][0,:]
                df_temp=h['montemp'][:]
                try:
                    
                    df_corr=h['rasocorrmon'][:]
                except:
                    df_corr=np.zeros_like(df_temp)
                    
                df_press=h['press'][:]
                df_good=h['goodmon'][:]
                df_temp[df_temp>400]=np.nan
                df_lat=h['lat'][:]
                df_lon=h['lon'][:]
                mask=np.zeros((2,3,df_time.shape[0]),dtype=bool)
        except:
            print('file not found: ',statlist)
            return


        #
        ##
        ###
        ####
        #####
        if adj in [None, 'bg','rio', 'rit']:
            df_temp += df_corr
        elif adj == 'raobcore':
            pass
        else:
            print('not a valid adjustment: ', adj)
            return
        #####
        ####
        ###
        ##
        #

    #     print('adj_df: ', adj_df)
    #     print('df: ',df)
        if adjstatlist != None:
            for it in range(1,adj_time.shape[0]-1):
                idb=np.searchsorted(df_time,adj_time[it-1:it+1])
                for ih in [0,1]:
                    for ip in range(df_press.shape[0]):
                        #print('adjusting: ', ip)
                        df_temp[ih,ip,idb[0]:idb[1]]-=adj_corr[ih,ip,it]


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

    #     print(all_dfta)

        shpsave=-1
        ei34=dict()
        ei=dict()
        for k in era_input.keys():
            if k not in ['mon']:
                ei34[k]=None
                ei[k]=None
        
        indices=[[],[]]        
        indices2=[[],[]]        
        indices34=[[],[]]    
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
                    #if imon<12:
                        #y=np.datetime64('{}-{:0>2}-{:0>2}'.format(yr,imon+1,1))-np.timedelta64(1,'s')
                    #else:
                        #y=np.datetime64('{}-{:0>2}-{:0>2}'.format(yr+1,1,1))-np.timedelta64(1,'s')
                    idx=np.searchsorted(df_time,x)
                    if idx==df_time.shape[0]:
                        continue
                    #print(df_time[idx],x)
                    if df_time[idx]!=x:
                        continue

                    indices[ih].append(idx)
                    if df_good[ih,3,idx] >= 3 and era_input['sp'][l]>0:                                

                        prof=df_temp[ih,2:13,idx]
                                
                        if not any(np.isnan(prof)):
                            #high_reduced_sh = reduced_sh[reduced_sh.index.isin(high_mon_mean.plev)]
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
            

            print(ih,len(tadata))
        if len(tadata34)==0 or len(tadata)==0:
            print('nans found, returning',statid)
            return
        print('before 34 c',statid)
        #if statid!='91413':
            #return

        for ih in 0,1:
            
            indices[ih]=np.sort(np.array(list(set(indices34[0]+indices2[0]))))

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

        #testb2 =b2[np.array(daydata)[np.array(chandata) == 2] == True]
        #testb34 = b34[np.array(daydata34)[np.array(chandata34) == 34] == True]
        #date_out2 = np.array(mondata)[np.array(chandata) == 2][np.array(daydata)[np.array(chandata) == 2] == True]
        #date_out34 = np.array(mondata34)[np.array(chandata34) == 34][np.array(daydata34)[np.array(chandata34) == 34] == True]

        try:
            os.makedirs("./"+odir+"/"+statid+"/")
        except:
            pass

        with netCDF4.Dataset('bt_template.nc') as f:
            if adj in ['rio','rit']:
                fno=adjstatlist[:-9]+'_bt2_'+adjstatlist[-9:]
                fno='_'.join(fno.split('__'))
                fno='mon'.join(fno.split('save'))
            elif adj=='raobcore':
                fno=statlist[:-9]+'_bt2_'+statlist[-9:]
            elif adj=='bg':
                fno=statlist[:-9]+'_bt2_'+statlist[-9:]
            else:
                fno=''.join(statlist.split('corr'))[:-9]+'_bt2_'+statlist[-9:]
            with netCDF4.Dataset(fno,'w') as g:
                for attname in f.ncattrs():
                    if attname == 'Stationname':
                        setattr(g,attname,station_name)
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
                
            
                fillvars = {}
                vars_to_write = ['montemp', 'goodmon']#, 'rasocorrmon', 'eracorrmon', 'ancorrmon']

                hilf=np.empty(g['montemp'].shape,g['montemp'].dtype)
                hilf.fill(np.nan)
                vals=np.sum(mask,axis=2) # needed to split b2 and b34
                #for ih in 0,1:
                    
                hilf[0,0,xmask[0,0,:]]=b2[:vals[0,0],1]
                hilf[0,1,xmask[0,1,:]]=b34[:vals[0,1],0]
                hilf[0,2,xmask[0,2,:]]=b34[:vals[0,1],1]
                hilf[1,0,xmask[1,0,:]]=b2[vals[0,0]:,1]
                hilf[1,1,xmask[1,1,:]]=b34[vals[0,1]:,0]
                hilf[1,2,xmask[1,2,:]]=b34[vals[0,2]:,1]
                    
                g['montemp'][:]=hilf[:]
                g['goodmon'][:]=0
                for ih in 0,1:
                    g['goodmon'][ih,0,:]=df_good[ih,3,gmask]
                    g['goodmon'][ih,1,:]=df_good[ih,2,gmask]
                    g['goodmon'][ih,2,:]=df_good[ih,2,gmask]
            

        #b_final = []
        #for i in wholemon[:middle_index]:
            #b = [[np.nan, np.nan, np.nan, np.nan]]
            #if i in date_out2:
                #b[0][0] = testb2[date_out2 == i][0][0]
                #b[0][1] = testb2[date_out2 == i][0][1]
            #if i in date_out34:
                #b[0][2] = testb34[date_out34 == i][0][0]
                #b[0][3] = testb34[date_out34 == i][0][1]
            #b_final.append(b)


        #pickle.dump( b_final, open( odir+"/"+statid+"/"+statid+"_day_refl.p", "wb" ) )
        #pickle.dump( dgoodmon, open( odir+"/"+statid+"/"+statid+"_day_goodmon.p", "wb" ) )
        #pickle.dump( wholemon[:middle_index], open( odir+"/"+statid+"/"+statid+"_day_dates.p", "wb" ) )
        print('day done: '+statid)


        #testb2 =b2[np.array(daydata)[np.array(chandata) == 2] == False]
        #testb34 = b34[np.array(daydata34)[np.array(chandata34) == 34] == False]
        #date_out2 = np.array(mondata)[np.array(chandata) == 2][np.array(daydata)[np.array(chandata) == 2] == False]
        #date_out34 = np.array(mondata34)[np.array(chandata34) == 34][np.array(daydata34)[np.array(chandata34) == 34] == False]

        #b_final = []
        #for i in wholemon[middle_index:]:
            #b = [[np.nan, np.nan, np.nan, np.nan]]
            #if i in date_out2:
                #b[0][0] = testb2[date_out2 == i][0][0]
                #b[0][1] = testb2[date_out2 == i][0][1]
            #if i in date_out34:
                #b[0][2] = testb34[date_out34 == i][0][0]
                #b[0][3] = testb34[date_out34 == i][0][1]
            #b_final.append(b)

        #pickle.dump( b_final, open( odir+"/"+statid+"/"+statid+"_night_refl.p", "wb" ) )
        #pickle.dump( ngoodmon, open( odir+"/"+statid+"/"+statid+"_night_goodmon.p", "wb" ) )
        #pickle.dump( wholemon[middle_index:], open( odir+"/"+statid+"/"+statid+"_night_dates.p", "wb" ) )
        print('night done: '+statid,time.time()-tt)
    except MemoryError as e:
        print(e,'nothing to calculate: '+statid)
        return
    return

if __name__ == '__main__': 
#     [3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000]
    consthum = np.array([3.90, 5.90, 9.17, 20.30,
                         85.00, 1064.00 , 2475.60, 6631.20,
                         15468.00, 21684.00, 35328.00 , 44220.00]
                   )/2.
    statlist = []
#     statlist = glob.glob('/mnt/ssdraid/scratch/leo/rise/1.0/exp02/*85469*/feedbackglobbincorrmon0*.nc')
#     statlist = glob.glob('/mnt/ssdraid/scratch/leo/rise/1.0/exp02/*17095*/feedbackglobbincorrmon0*.nc')

    statlist = glob.glob('/mnt/ssdraid/scratch/leo/rise/1.0/exp01/*/feedbackglobbincorrmon[01]*.nc')

    i = None
    

    def eraseries(statlist,yr):
        
        era_input={}
        tt=time.time()
        with netCDF4.Dataset('./era/era_'+str(yr)+'.nc') as h:
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
    
    pool = multiprocessing.Pool(processes=40)
    
    func=partial(eraseries, statlist[:])
    tt=time.time()
    result_list = list(pool.map(func,range(1950,2022)))
    era_input={}
    for s in result_list[0].keys():
        era_input[s]={}
        for k in result_list[0][s].keys():
            if k in ['latitude','longitude']:
                era_input[s][k]=result_list[0][s][k]
            else:
                era_input[s][k]=np.concatenate([result_list[i][s][k] for i in range(len(result_list))])
    print(time.time()-tt)

    odir='./'
    for i in ['bg',None, 'raobcore', 'rio', 'rit']:# [None, 'raobcore', 'rio', 'rit']:
        if i == None:
            odir = "rttov_out_unadj_testing"
        elif i == 'raobcore':
            odir = "rttov_out_raobcore_testing"
        elif i == 'rio':
            odir = "rttov_out_rio_testing"
        elif i == 'rit':
            odir = "rttov_out_rit_testing"
        try:
            os.makedirs("./"+odir+"/")
        except:
            pass
        skeys=list(era_input.keys())
        func=partial(calc_station, chum = consthum, adj = i, odir = odir)
        tt=time.time()
        result_list = list(pool.map(func, zip(statlist[:],[era_input[k] for k in skeys[:]])))
        print(time.time()-tt)
        print('finished')
        
    '''
    ['/rttov/rtcoef_rttov12/rttov7pred54L/rtcoef_noaa_8_msu.dat', 
    '/rttov/rtcoef_rttov12/rttov7pred54L/rtcoef_noaa_17_amsub.dat',
    '/rttov/rtcoef_rttov12/rttov7pred54L/rtcoef_noaa_16_amsub.dat',
    '/rttov/rtcoef_rttov12/rttov7pred54L/rtcoef_noaa_18_amsua.dat', 
    '/rttov/rtcoef_rttov12/rttov7pred54L/rtcoef_noaa_15_amsub.dat', 
    '/rttov/rtcoef_rttov12/rttov7pred54L/rtcoef_noaa_16_amsua.dat', 
    '/rttov/rtcoef_rttov12/rttov7pred54L/rtcoef_noaa_9_msu.dat',
    '/rttov/rtcoef_rttov12/rttov7pred54L/rtcoef_noaa_11_msu.dat',
    '/rttov/rtcoef_rttov12/rttov7pred54L/rtcoef_noaa_6_msu.dat', 
    '/rttov/rtcoef_rttov12/rttov7pred54L/rtcoef_noaa_14_msu.dat', 
    '/rttov/rtcoef_rttov12/rttov7pred54L/rtcoef_noaa_7_msu.dat', 
    '/rttov/rtcoef_rttov12/rttov7pred54L/rtcoef_noaa_19_amsua.dat',
    '/rttov/rtcoef_rttov12/rttov7pred54L/rtcoef_noaa_15_amsua.dat',
    '/rttov/rtcoef_rttov12/rttov7pred54L/rtcoef_noaa_12_msu.dat',
    '/rttov/rtcoef_rttov12/rttov7pred54L/rtcoef_noaa_17_amsua.dat',
    '/rttov/rtcoef_rttov12/rttov7pred54L/rtcoef_noaa_10_msu.dat', 
    '/rttov/rtcoef_rttov12/rttov7pred54L/rtcoef_noaa_5_msu.dat']
    '''