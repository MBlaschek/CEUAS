from __future__ import print_function

### pre-WMOBUFR: export ECCODES_DEFINITION_PATH=/home/erc/Work/mydefs:`codes_info -d` 
# WMOBUFR: export ECCODES_DEFINITION_PATH=/home/erc/Work/mydefsWMOBUFR:`codes_info -d` 



#from encode_TM309099 import bufr_encode # pre-WMOBUFR
from encode_TM309099WMO import bufr_encode # WMOBUFR
from bufr_utils import open_clean_bufr_files, close_clean_bufr_files, get_synoptic_date_and_time
import traceback
import sys, subprocess
import numpy as np
import datetime as dt
import zipfile
import xarray as xr
import pandas as pd
import os, copy
import cdsapi
from multiprocessing import Pool as multiPPool, RawArray as multiPArray, current_process as multiPCurrentProcess; nprocesses=8
from entropy import crawl_slice, function_sameprofileas

df_WMO_SCH = pd.read_csv('CUON_wmo_sch_mapping.tsv',sep='\t').set_index(['sch_id']).sort_index()
SCH_idx = df_WMO_SCH.index.get_level_values('sch_id').unique()

#data_acquisition_dir='CUON_DISPLACEMENT_NEW'
#data_acquisition_dir='CUON_DISPLACEMENT_RENEW'
data_acquisition_dir='CUON_DISPLACEMENT_FINAL'
#bufr_output_dir='CUON_DISPLACEMENT_FINAL' # pre-WMOBUFR
bufr_output_dir='CUON_WMO' # WMOBUFR

#timeoffset_name='_minus30minutes'; timeoffset_value=np.timedelta64(-1800,'s') # DEBUG
timeoffset_name=''; timeoffset_value=None

from io import StringIO 

class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio    # free up some memory
        sys.stdout = self._stdout

#def main():
if 1==1:
    if len(sys.argv)==2:
      yyyymmdd=int(sys.argv[1])
    else:
      #yyyymmdd=19800601
      yyyymmdd=20061001
    CUONPATH = '/ec/res4/scratch/erc/{0}/'.format(data_acquisition_dir)
    yyyymmdd_dt = dt.datetime.strptime(str(yyyymmdd),'%Y%m%d')
    print("==========================================================")
    print("CUON ENCODING FOR",yyyymmdd_dt)

    df_BIGcuon = None

    # Time window for BUFR encoding    
    min_datetime_exact  = np.datetime64(yyyymmdd_dt.strftime('%Y-%m-%d %H:%M:%S'))
    max_datetime_exact  = np.datetime64((yyyymmdd_dt + dt.timedelta(days= 1)).strftime('%Y-%m-%d %H:%M:%S'))

    # Time window for searching for duplicates
    min_datetime_margin = np.datetime64((yyyymmdd_dt - dt.timedelta(hours= 1)).strftime('%Y-%m-%d %H:%M:%S'))
    max_datetime_margin = np.datetime64((yyyymmdd_dt + dt.timedelta(hours=25)).strftime('%Y-%m-%d %H:%M:%S'))

    pkl_dup = "/ec/res4/scratch/erc/Leo/BUFRENCODE/perm/erc/ERA6BUFR/{0}/pickledup/{1}/pickledup_{2}.pkl.gz".format(data_acquisition_dir,yyyymmdd_dt.strftime('%Y'),yyyymmdd_dt.strftime('%Y%m%d'))
    pkl_dup_path = os.path.dirname(pkl_dup)
    if not os.path.exists(pkl_dup_path):
      os.makedirs(pkl_dup_path)

    for yyyymmdd_dt1 in [yyyymmdd_dt-dt.timedelta(days=1), yyyymmdd_dt, yyyymmdd_dt+dt.timedelta(days=1)]:
      pkl_file1 = "/ec/res4/scratch/erc/Leo/BUFRENCODE/perm/erc/ERA6BUFR/{0}/pickle/{1}/pickle_{2}.pkl.gz".format(data_acquisition_dir,yyyymmdd_dt1.strftime('%Y'),yyyymmdd_dt1.strftime('%Y%m%d'))
      pkl_path1 = os.path.dirname(pkl_file1)
      if not os.path.exists(pkl_path1):
        os.makedirs(pkl_path1)

      lread=1 # DEBUG
      if lread and not os.path.exists(pkl_file1):

        print(pkl_file1)
        stop

        # GET CUON
        cuon_vars = ['air_temperature','dew_point_temperature','wind_speed','wind_from_direction','geopotential']
        optn_vars = {'air_temperature':['RISE_bias_estimate','sonde_type','latitude_displacement','longitude_displacement','time_since_launch'], \
                     'dew_point_temperature':['humidity_bias_estimate'], \
                     'wind_speed':['wind_bias_estimate'], \
                     'wind_from_direction':['wind_bias_estimate'], \
                     }
        ren_vars = {'wind_speed':{'wind_bias_estimate':'wind_speed_bias_estimate'}, \
                    'wind_from_direction':{'wind_bias_estimate':'wind_from_direction_bias_estimate'}, \
                    'air_temperature':{'RISE_bias_estimate':'air_temperature_bias_estimate','latitude_displacement':'latitudeDisplacement','longitude_displacement':'longitudeDisplacement','time_since_launch':'timePeriod'}, \
                    'dew_point_temperature':{'humidity_bias_estimate':'dew_point_temperature_bias_estimate'}}
        df_cuon = None
        orph_list=[]; orph_counter=10000; orph_map=[] # ill-formed WIGOS IDs require special handling
        cuonfile_msk = CUONPATH + 'CUON_{{0}}_{0}.zip'.format(yyyymmdd_dt1.strftime('%Y%m%d'))
        CUONEXTR     = CUONPATH + '{0}/'.format(              yyyymmdd_dt1.strftime('%Y%m%d'))
        if not os.path.exists(CUONEXTR):
          os.makedirs(CUONEXTR)
        cuonfiles    = []
        for var1 in cuon_vars:
          df_cuon_var1 = None
          vars_to_save = {'air_temperature':['ta','RISE_bias_estimate','sonde_type','latitude_displacement','longitude_displacement','time_since_launch'],'dew_point_temperature':['dew_point_temperature','humidity_bias_estimate'],'wind_speed':['wind_speed','wind_bias_estimate'],'wind_from_direction':['wind_from_direction','wind_bias_estimate'],'geopotential':['geopotential']}[var1]
          cuonfile = cuonfile_msk.format(var1)
          if not os.path.exists(cuonfile):
            obsdataset = 'insitu-comprehensive-upper-air-observation-network'
            req_fmt = 'nc'
            obsret_args = {
            'variable': [var1],
            'format': req_fmt,
            'dummy': 'avoid_cach_'+str(dt.datetime.now()),
            "cdm":["observations_table/source_id","station_configuration/platform_type"],
            }
            if var1 in optn_vars.keys():
              obsret_args['optional'] = optn_vars[var1]
            else:
              obsret_args['optional'] = []
            obsret_args['optional'] += ['station_elevation']
            cdsapi_file = cuonfile+'-request.py'
            fcdsapi = open(cdsapi_file,'w')
            fcdsapi.write('import os\nimport cdsapi\nos.environ["CDSAPI_RC"]="/home/erc/.cdsapirc_test"\n')
            fcdsapi.write('c = cdsapi.Client(timeout=600,quiet=False,debug=True)\n')
            cargs = (obsdataset, {**obsret_args, **{'date': [yyyymmdd_dt1.strftime('%Y%m%d')]}}, cuonfile,)
            fcdsapi.write('c.retrieve("{0}",{4} {1},"date":[{2}] {5}, "{3}")\n'.format(obsdataset, ', '.join(['"{0}":{1}'.format(k, '"{0}"'.format(obsret_args[k]) if (type(obsret_args[k]) is str) else '{0}'.format(obsret_args[k])) for k in obsret_args.keys()]), yyyymmdd_dt1.strftime('%Y%m%d'), cuonfile, '{', '}'))
            fcdsapi.close()
            #c.retrieve(*cargs)
            cdsapi_stderr = subprocess.run(['python3',cdsapi_file], capture_output=True).stderr.strip().decode('utf-8')
            os.remove(cdsapi_file)
            lnodata = ("Try different pressure levels or variable" in cdsapi_stderr) or \
                      ("Something went wrong in the data provider service, check your query and try later" in cdsapi_stderr)
            if lnodata:
              print("No file for",var1,yyyymmdd_dt1.strftime('%Y%m%d'))
              os.system("touch {0}".format(cuonfile))
          else:
            print("file already acquired:",cuonfile)
 
          if not os.path.exists(cuonfile):
            print("ERROR... file not retrieved...")
            print(cdsapi_stderr)
            stop
   
          if os.path.getsize(cuonfile)>0:
            with zipfile.ZipFile( cuonfile, "r") as zip_ref:
              listOfFileNames = zip_ref.namelist()
              listOfFileNames_var = [x for x in listOfFileNames if var1 in x]
              zip_ref.extractall(members=listOfFileNames_var, path=CUONEXTR)
          else:
            listOfFileNames_var = []
          cuonfiles.append(cuonfile)
 
          for fic1 in np.sort(listOfFileNames_var):
            #print("reading",CUONEXTR+fic1)
            ds = xr.open_dataset(CUONEXTR+fic1); nds = ds.sizes['obs']
            vars_save_from_found = [x for x in vars_to_save if x in ds.variables] + ['source_id','station_elevation']
            # Check the WIGOS ID and handle missing ones
            wigos_statid1=os.path.basename(fic1).replace('dest_','').replace('_'+var1,'').replace('.nc','')
            wigos_statid1_elms = wigos_statid1.split('-')
            if ('coordinate-orphan' in wigos_statid1) or (len(wigos_statid1_elms)!=4):
              if wigos_statid1 not in orph_list:
                orph_counter-=1
                new_wigos_statid1 = '0-23001-2-orph{0:04d}'.format(orph_counter)
                orph_list.append(wigos_statid1)
                orph_map.append(new_wigos_statid1)
              else:
                new_wigos_statid1 = orph_map[orph_list.index(wigos_statid1)]
            else:
              new_wigos_statid1 = copy.copy(wigos_statid1)
            # Check station_id
            vars_save_extra_opt_args = {'wigos_statid':np.full((nds,), new_wigos_statid1, '<U32')}
            vars_save_for_identification = ['time','lat','lon','plev','report_id'] # it seems report_id now changes!!
            try:
              if ('station_id' not in ds.variables) or (ds['station_id'].values[:1].astype(str)[0]!=wigos_statid1):
                vars_save_extra_opt_args['station_id'] = np.full((nds,), new_wigos_statid1, 'S16')
              else:
                vars_save_for_identification.append('station_id')
            except:
              # if the station ID contains non-ASCII characters (rare, but happens, e.g. 19800604/dest_0-20999-0-UH....S_air_temperature.nc)
              # then we need to skip the file altogether
              print("ERROR... Non-ASCII character in station_id",str(ds['station_id'].values[0].decode('utf-8')),"skipping file",CUONEXTR+fic1)
              continue
            # Check platform_type
            try:
              platform_types = np.unique(xr.open_dataset(CUONEXTR+fic1, group='station_configuration')['platform_type'].values[:])
            except:
              platform_types = [0] # Fall-back: assume LAND
            if len(platform_types) != 1 :
              print("ERROR... non-unique platform type...")
              stop
            #DEBUG replace missing platform_type for SHIP
            #if platform_types[0]==-2147483648:
            #  platform_types[0]=2
            vars_save_extra_opt_args['platform_type'] = np.full((nds,), platform_types[0], int)
            # Load the data into a dataframe
            df_cuon_tmp = pd.DataFrame( { **{x:ds[x].values[:] for x in vars_save_for_identification+vars_save_from_found}, **vars_save_extra_opt_args } )

            # Roll back longitudes>=180
            wlon180 = np.where(df_cuon_tmp['lon']>=180.)[0]
            if len(wlon180)>0:
              ilon = list(df_cuon_tmp.columns).index('lon')
              df_cuon_tmp.iloc[wlon180,ilon] -= 360.

            # Ensure longitude_displacement is in the range [-180,180)
            if 'longitude_displacement' in vars_save_from_found:
              wlondispp180 = np.where((~np.isnan(df_cuon_tmp['longitude_displacement'].values[:])) & \
                                      (df_cuon_tmp['longitude_displacement'].values[:]>=180.))[0]
              if len(wlondispp180)>0:
                ilondisp = list(df_cuon_tmp.columns).index('longitude_displacement')
                df_cuon_tmp.iloc[wlondispp180,ilondisp] -= 360.
              wlondispm180 = np.where((~np.isnan(df_cuon_tmp['longitude_displacement'].values[:])) & \
                                      (df_cuon_tmp['longitude_displacement'].values[:]<-180.))[0]
              if len(wlondispm180)>0:
                ilondisp = list(df_cuon_tmp.columns).index('longitude_displacement')
                df_cuon_tmp.iloc[wlondispm180,ilondisp] += 360.

            if var1 in ren_vars.keys():
              df_cuon_tmp = df_cuon_tmp.rename(columns=ren_vars[var1])
            if df_cuon_var1 is None:
              df_cuon_var1 = df_cuon_tmp.copy()
            else:
              df_cuon_var1 = pd.concat([df_cuon_var1, df_cuon_tmp], axis=0) # add rows

            # Close file and remove
            ds.close()
            os.remove(CUONEXTR+fic1)

          # Form an indexed dataframe
          if df_cuon_var1 is None:
            empty_vars = ['wigos_statid','report_id','plev','time','lat','lon','station_id','source_id','station_elevation','platform_type']
            for var11 in vars_to_save:
              if var1 in ren_vars.keys() and var11 in ren_vars[var1].keys():
                empty_vars.append(ren_vars[var1][var11])
              else:
                empty_vars.append(var11)
            df_cuon_var1 = pd.DataFrame({k:[] for k in empty_vars})

          df_cuon_var1.set_index(['wigos_statid','report_id','plev','time','lat','lon','station_id'],inplace=True)
          #df_cuon_var1.set_index(['wigos_statid','time','plev','lat','lon','station_id'],inplace=True)
          if df_cuon is None:
            df_cuon = df_cuon_var1.copy()
          else:
            recon_vars = ['source_id','station_elevation','platform_type']
            df_cuon = df_cuon.join(df_cuon_var1.rename(columns={x:x+'NEW' for x in recon_vars}), how='outer')
            wnew_source_ID = np.where((df_cuon['source_idNEW'].values[:].astype(str)!='nan') & (df_cuon['source_id'].values[:].astype(str)=='nan'))[0]
            if len(wnew_source_ID)>0:
              isource_id = list(df_cuon.columns).index('source_id')
              df_cuon.iloc[wnew_source_ID, isource_id] = df_cuon['source_idNEW'].values[wnew_source_ID]
            wbad_source_ID = np.where((df_cuon['source_idNEW'].values[:].astype(str)!='nan') & (df_cuon['source_id'].values[:].astype(str)!='nan') & (df_cuon['source_idNEW'].values[:].astype(str)!=df_cuon['source_id'].values[:].astype(str)))[0]
            if len(wbad_source_ID)>0:
              print("Here we go... source_ids do not match...",df_cuon.iloc[wbad_source_ID][['source_id','source_idNEW']].values[:])
              stop
            df_cuon.drop(columns=['source_idNEW'], inplace=True)
            for recon_var1 in recon_vars:
              if recon_var1!='source_id':
                wnew_recon_var = np.where((~np.isnan(df_cuon[recon_var1+'NEW'].values[:])) & (np.isnan(df_cuon[recon_var1].values[:])))[0]
                if len(wnew_recon_var)>0:
                  irecon_var1 = list(df_cuon.columns).index(recon_var1)
                  df_cuon.iloc[wnew_recon_var, irecon_var1] = df_cuon[recon_var1+'NEW'].values[wnew_recon_var]
                wbad_recon_var = np.where((~np.isnan(df_cuon[recon_var1+'NEW'].values[:])) & (~np.isnan(df_cuon[recon_var1].values[:])) & (df_cuon[recon_var1+'NEW'].values[:]!=df_cuon[recon_var1].values[:]))[0]
                if len(wbad_recon_var)>0:
                  print("Here we go... ",recon_var," values do not match...",df_cuon.iloc[wbad_recon_var][[recon_var1,recon_var1+'NEW']].values[:])
                  stop
                df_cuon.drop(columns=[recon_var1+'NEW'], inplace=True)
 
        # RENAME COLUMNS AND QC
        wgoodlatlon = np.where((~np.isnan(df_cuon.index.get_level_values('lat').values[:])) & (~np.isnan(df_cuon.index.get_level_values('lon').values[:])) & (df_cuon.index.get_level_values('lat').values[:]>=-90.) & (df_cuon.index.get_level_values('lat').values[:]<=90.) & (df_cuon.index.get_level_values('lon').values[:]>=-180.) & (df_cuon.index.get_level_values('lon').values[:]<=180.))[0]
        df_cuon = df_cuon.iloc[wgoodlatlon]
        wgoodplev = np.where((~np.isnan(df_cuon.index.get_level_values('plev').values[:])) & (df_cuon.index.get_level_values('plev').values[:]>=50) & (df_cuon.index.get_level_values('plev').values[:]<125000))[0]
        df_cuon = df_cuon.iloc[wgoodplev]
        #wgoodwigosid = np.where(~np.char.endswith(df_cuon.index.get_level_values('wigos_statid').values[:].astype(str),'coordinate-orphan'))[0]
        #df_cuon = df_cuon.iloc[wgoodwigosid]
        df_cuon.rename(columns={'air_temperature_bias_estimate':'airTempBiasCorr','ta':'airTemperature','dew_point_temperature':'dewpointTemperature','dew_point_temperature_bias_estimate':'dewpointTempBiasCorr'},inplace=True)
        wbadbiasT = np.where((df_cuon['airTempBiasCorr'].values[:]>50.) | (df_cuon['airTempBiasCorr'].values[:]<-50.))[0]
        if len(wbadbiasT)>0: # T bias correction unreasonable
          iairTempBiasCorr = list(df_cuon.columns).index('airTempBiasCorr')
          df_cuon.iloc[wbadbiasT, iairTempBiasCorr] = np.nan
        wbadbiasdewT = np.where((df_cuon['dewpointTempBiasCorr'].values[:]>50.) | (df_cuon['dewpointTempBiasCorr'].values[:]<-50.))[0]
        if len(wbadbiasdewT)>0: # dewT bias correction unreasonable
          idewpointTempBiasCorr = list(df_cuon.columns).index('dewpointTempBiasCorr')
          df_cuon.iloc[wbadbiasdewT, idewpointTempBiasCorr] = np.nan
        df_cuon.rename(columns={'wind_from_direction_bias_estimate':'windDirectionBiasCorr','wind_from_direction':'windDirection'},inplace=True)
        wbaddd = np.where((df_cuon['windDirection'].values[:]<0.) | (df_cuon['windDirection'].values[:]>360.))[0]
        if len(wbaddd)>0: # negative or > 360 wind direction, not foreseen
          iwindDirection = list(df_cuon.columns).index('windDirection')
          df_cuon.iloc[wbaddd, iwindDirection] = np.nan
        df_cuon.rename(columns={'wind_speed':'windSpeed','wind_speed_bias_estimate':'windSpeedBiasCorr'},inplace=True)
        wbadff = np.where((df_cuon['windSpeed'].values[:]<0.) | (df_cuon['windSpeed'].values[:]>113.))[0]
        if len(wbadff)>0: # wind speed negative or unreasonably large
          iwindSpeed = list(df_cuon.columns).index('windSpeed')
          df_cuon.iloc[wbadff, iwindSpeed] = np.nan
        wbadT = np.where((df_cuon['airTemperature'].values[:]<100.) | (df_cuon['airTemperature'].values[:]>333.))[0]
        if len(wbadT)>0: # unreasonable temperature
          iairTemperature = list(df_cuon.columns).index('airTemperature')
          df_cuon.iloc[wbadT, iairTemperature] = np.nan
        wbaddewT = np.where((df_cuon['dewpointTemperature'].values[:]<100.) | (df_cuon['dewpointTemperature'].values[:]>350.))[0]
        if len(wbaddewT)>0: # unreasonable dew point temperature
          idewpointTemperature = list(df_cuon.columns).index('dewpointTemperature')
          df_cuon.iloc[wbaddewT, idewpointTemperature] = np.nan
        wokgeop = np.where(~np.isnan(df_cuon['geopotential'].values[:]))[0]
        geop_height = np.full( (len(df_cuon),), np.nan)
        zg0 = 9.80665 # WMO value for g
        geop_height[wokgeop] = df_cuon['geopotential'].values[wokgeop] / zg0
        df_cuon['nonCoordinateGeopotentialHeight'] = geop_height
        wbadz = np.where((df_cuon['nonCoordinateGeopotentialHeight'].values[:]>40.e3) | (df_cuon['nonCoordinateGeopotentialHeight'].values[:]<-500))[0] # -500 m to 40 km
        if len(wbadz)>0: # unreasonable geopotential height
          inonCoordinateGeopotentialHeight = list(df_cuon.columns).index('nonCoordinateGeopotentialHeight')
          df_cuon.iloc[wbadz, inonCoordinateGeopotentialHeight] = np.nan
        df_cuon.rename(columns={'station_elevation':'height'},inplace=True)
        # TODO specify if these are standard or significant levels ## THIS IS REVISED IN ENCODING BELOW (TO HANDLE THE CASE OF HIGH-RESOLUTION BUFR AND TO MARK SURFACE LEVELS)
        signifs = np.full((len(df_cuon),), 0., dtype=float)
        std_plevels = np.array([1000, 925, 850, 700, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10, 7, 5, 3, 2, 1])*100. # in Pa
        wstd = np.where(np.isin(df_cuon.index.get_level_values('plev').values[:],std_plevels))[0]
        signifs[wstd] += 2**(18-2)
        w_nonstd_t_present = np.where((~np.isin(df_cuon.index.get_level_values('plev').values[:],std_plevels)) & (~np.isnan(df_cuon['airTemperature'].values[:])))[0]
        signifs[w_nonstd_t_present] += 2**(18-5)
        w_nonstd_wind_present = np.where((~np.isin(df_cuon.index.get_level_values('plev').values[:],std_plevels)) & (~np.isnan(df_cuon['windDirection'].values[:])) & (~np.isnan(df_cuon['windSpeed'].values[:])))[0]
        signifs[w_nonstd_wind_present] += 2**(18-7)
        w_all_missing = np.where(signifs==0.)[0]
        signifs[w_all_missing] += 1
        df_cuon['extendedVerticalSoundingSignificance'] = signifs
        df_cuon = df_cuon.reset_index()
    
        ## CALCULATE DEWPOINT TEMPERATURE
        #coeff_a=17.625; coeff_b=243.04 # Alduchov and Eskridge (1996)
        ## Ts = (b × α(T,RH)) / (a - α(T,RH)) where a = 17.625 and b = 243.04 Celsius and alpha = ln(RH/100) + a*T/(b+T) with T in degrees C
        #woktq=np.where((~np.isnan(df_cuon['airTemperature'].values[:])) & (~np.isnan(df_cuon['hur'].values[:])))[0]
        #dewpoint = np.full( (len(df_cuon),), np.nan )
        #temp_C = df_cuon['airTemperature'].values[woktq] - 273.15
        #rh_100 = df_cuon['hur'].values[woktq]*100.
        #alpha = np.log(rh_100 / 100.) + coeff_a*temp_C/(coeff_b+temp_C)
        #dewpoint[woktq] = (coeff_b*alpha)/(coeff_a-alpha) + 273.15
        #df_cuon['dewpointTemperature'] = dewpoint
    
        # INDEX AND VERIFY UNICITY FOR HEADER ELEMENTS
        df_cuon = df_cuon.set_index(['wigos_statid','report_id','plev']).sort_index(ascending=[True,True,False])
        #df_cuon = df_cuon.set_index(['wigos_statid','time','plev']).sort_index(ascending=[True,True,False])
        n1 = len(df_cuon.index.unique())
        df_cuon['lat05'] = np.around(df_cuon.reset_index()['lat'].values[:]/5,decimals=1)*5.
        df_cuon['lon05'] = np.around(df_cuon.reset_index()['lon'].values[:]/5,decimals=1)*5.
        n2 = len(df_cuon.reset_index().set_index(['wigos_statid','report_id','plev','time','lat05','lon05','station_id']).index.unique())
        #n2 = len(df_cuon.reset_index().set_index(['wigos_statid','time','plev','lat05','lon05','station_id']).index.unique())
        if n1!=n2:
          print("WARNING... non-unique matching in indices...",n1,n2)
          print("example")
          for dimex in ['lat','lon']:
            nuq=df_cuon.reset_index().groupby(['wigos_statid','report_id','plev','time','station_id'])[dimex + '05'].nunique()
            #nuq=df_cuon.reset_index().groupby(['wigos_statid','time','plev','station_id'])[dimex + '05'].nunique()
            nuq_2=np.where(nuq>1)[0]
            if len(nuq_2)>0:
              print(dimex.upper()+" changes in the profile as follows:")
              #df_cuon.reset_index().set_index(['wigos_statid','report_id','plev','time','station_id']).loc[nuq.iloc[nuq_2].index][['lat', 'lon', 'airTemperature','dewpointTemperature','windDirection','windSpeed']].reset_index()
              df_cuon.reset_index().set_index(['wigos_statid','time','plev','station_id']).loc[nuq.iloc[nuq_2].index][['lat', 'lon', 'airTemperature','dewpointTemperature','windDirection','windSpeed']].reset_index()
    
        df_cuon = df_cuon.reset_index().rename(columns={'plev':'pressure'})
        #df_cuon = df_cuon.reset_index().set_index(['wigos_statid','time']).rename(columns={'plev':'pressure'})
        # END OF READING
  
        # SAVE DATA
        df_cuon.to_pickle(pkl_file1,protocol=3,compression="gzip")
        print("saved pickle in",pkl_file1)

        # Done reading all the data: Clean-up zip files and directory
        for cuonfile1 in cuonfiles:
          if os.path.exists(cuonfile1):
            os.remove(cuonfile1)
        if os.path.exists(CUONEXTR):
          os.rmdir(CUONEXTR)

      elif lread and os.path.exists(pkl_file1):
  
        df_cuon = pd.read_pickle(pkl_file1,compression="gzip")
        print("read from",pkl_file1,len(df_cuon))

      if len(df_cuon)>0:

        # Retain only data that will be encoded for that day plus the margin needed for duplicate checks
        wtime = np.where((df_cuon['time'].values[:]>=min_datetime_margin) & \
                         (df_cuon['time'].values[:]<=max_datetime_margin))[0]
        print("keeping",len(wtime))
        if df_BIGcuon is None:
          df_BIGcuon = df_cuon.iloc[wtime].copy()
        else:
          df_BIGcuon = pd.concat([df_BIGcuon, df_cuon.iloc[wtime].copy()], axis=0)

    if not os.path.exists(pkl_dup):
      # LOOK FOR DUPLICATES
      myindex_make = ['hour','lat01i','lon01i']
      time_axis='time'
      df_BIGcuon['hour'] = ((df_BIGcuon[time_axis].values[:]-min_datetime_margin).astype(int)//1e9//3600).astype(int) # TODO USE ONLY LAUNCH TIME
      df_BIGcuon['lat01i'] = np.around(df_BIGcuon['lat'],0).astype(int)
      df_BIGcuon['lon01i'] = np.around(df_BIGcuon['lon'],0).astype(int)
      # Roll back longitudes>=180
      wlon180 = np.where(df_BIGcuon['lon01i']>=180.)[0]
      if len(wlon180)>0:
        ilon = list(df_BIGcuon.columns).index('lon01i')
        df_BIGcuon.iloc[wlon180,ilon] -= 360
      df_BIGcuon = df_BIGcuon.set_index(myindex_make).sort_index()
      pairs_done = []; discarded_list = []; df_dup = {'wigos_statid':[], 'report_id':[]}
      function_call = {'nbcolumns':len(df_BIGcuon.columns)+3, 'args':(pairs_done,discarded_list,df_dup,), 'function':function_sameprofileas}
      crawl_slice(df_BIGcuon.reset_index(),df_BIGcuon.index.names,function_call)
      print("Number of profiles to remove",len(discarded_list))
      df_BIGcuon.reset_index(inplace=True)
      df_BIGcuon.drop(columns=['hour','lat01i','lon01i'], inplace=True)
      if len(df_dup['wigos_statid'])>0:
        df_dup = pd.DataFrame(df_dup)
        df_dup.to_pickle(pkl_dup, protocol=3, compression="gzip")
        print("saved duplicated checks in",pkl_dup)
      else:
        if os.path.exists(pkl_dup):
          os.system("rm -f {0}".format(pkl_dup))
        os.system("touch {0}".format(pkl_dup))
        print("empty duplicated checks",pkl_dup)
        df_dup = None
    else:
      if os.path.getsize(pkl_dup)>0:
        df_dup = pd.read_pickle(pkl_dup, compression="gzip")
      else:
        df_dup = None

    # SELECT ONLY DATA FOR THE DAY
    w_day = np.where((df_BIGcuon['time'].values[:]>=min_datetime_exact) & \
                     (df_BIGcuon['time'].values[:]< max_datetime_exact))[0]
    if len(w_day)>0:
      df_BIGcuon = df_BIGcuon.iloc[w_day]
    else:
      print("NO DATA FOUND FOR THE DAY... finishing here...")
      sys.exit(0)

    # SELECT ONLY DATA WHICH ARE NOT DUPLICATES
    if df_dup is not None:
      wwrite = np.where(np.count_nonzero( \
                  np.isin(df_BIGcuon[['wigos_statid','report_id']].values[:].astype(str), \
                          df_dup    [['wigos_statid','report_id']].values[:].astype(str)), axis=1) != 2 )[0]
      if len(wwrite)>0:
        df_BIGcuon = df_BIGcuon.iloc[wwrite]
      else:
        print("NO DATA LEFT AFTER APPLICATION OF DUPLICATE CHECK... THIS IS HIGHLY SUSPICIOUS... CHECK...!!!!")
        sys.exit(-1)
    else:
      print("No duplicates to worry about...")

    # NOW PROCEED TO WRITE
    df_write = df_BIGcuon.set_index(['wigos_statid','report_id'])

    # OFFSET DATA BACK IN TIME
    if timeoffset_value is not None:
      oldtime = df_BIGcuon['time'].values[:]
      newtime = oldtime + timeoffset_value
      df_BIGcuon['time'] = newtime
      print('BEFORE time offset',len(df_BIGcuon))
      # NEED TO REDO THE TIME SELECTION
      w_day = np.where((df_BIGcuon['time'].values[:]>=min_datetime_exact) & \
                       (df_BIGcuon['time'].values[:]< max_datetime_exact))[0]
      df_BIGcuon = df_BIGcuon.iloc[w_day]
      print('AFTER time offset',len(df_BIGcuon))

    # MAKE A SUMMARY
    pkl_sum = "/ec/res4/scratch/erc/Leo/BUFRENCODE/perm/erc/ERA6BUFR/{0}/picklesum/{1}/picklesum_{2}.pkl.gz".format(data_acquisition_dir,yyyymmdd_dt.strftime('%Y'),yyyymmdd_dt.strftime('%Y%m%d'))
    pkl_sum_path = os.path.dirname(pkl_sum)
    if not os.path.exists(pkl_sum_path):
      os.makedirs(pkl_sum_path)
    if not os.path.exists(pkl_sum):
      df_summary = None
      for utc1 in [0,6,12,18]:
        if utc1==0:
          wutc1 = np.where(df_write['time'].values<=np.datetime64(yyyymmdd_dt.strftime('%Y-%m-%d 03:00:00')))[0]
        elif utc1==6:
          wutc1 = np.where((df_write['time'].values> np.datetime64(yyyymmdd_dt.strftime('%Y-%m-%d 03:00:00'))) & \
                           (df_write['time'].values<=np.datetime64(yyyymmdd_dt.strftime('%Y-%m-%d 09:00:00'))))[0]
        elif utc1==12:
          wutc1 = np.where((df_write['time'].values> np.datetime64(yyyymmdd_dt.strftime('%Y-%m-%d 09:00:00'))) & \
                           (df_write['time'].values<=np.datetime64(yyyymmdd_dt.strftime('%Y-%m-%d 15:00:00'))))[0]
        elif utc1==18:
          wutc1 = np.where(df_write['time'].values> np.datetime64(yyyymmdd_dt.strftime('%Y-%m-%d 15:00:00')))[0]
        else:
          print("ERROR... utc1 unexpected",utc1)
          stop
        if len(wutc1)>0:
          ndata = df_write.iloc[wutc1].groupby(df_write.index.names).count()
          reordered_columns = ['pressure', 'time', 'lat', 'lat05', 'lon', 'lon05', 'station_id', 'source_id', 'height', 'platform_type', 'sonde_type', 'extendedVerticalSoundingSignificance', 'latitudeDisplacement', 'longitudeDisplacement', 'timePeriod', 'airTemperature', 'airTempBiasCorr', 'dewpointTemperature', 'dewpointTempBiasCorr', 'windSpeed', 'windSpeedBiasCorr', 'windDirection', 'windDirectionBiasCorr', 'geopotential', 'nonCoordinateGeopotentialHeight']
          if len(reordered_columns)!=len(ndata.columns):
            print("PROBLEM... reordered columns not the same length as columns...")
            print("reordered_columns",reordered_columns)
            print("ndata columns    ",ndata.columns)
            stop
          ndata = ndata[reordered_columns]
      
          df_ndata = pd.DataFrame({('Nobs','Total'):ndata.sum().astype(int)})
          df_ndata.index.name = 'variable'
          df_nprof = None
          for nquant1 in [0,1,5,50,95,99,100]:
            df1 = pd.DataFrame({('Nlev','P{0:03d}'.format(nquant1)):ndata.quantile(float(nquant1)/100.).astype(int)})
            df1.index.name = 'variable'
            df_ndata = pd.concat([ df_ndata, df1 ], axis=1)
            df_nprof = pd.DataFrame({('Nprof','Total'):np.count_nonzero(ndata.transpose(),axis=1)})
            df_nprof.index = ndata.columns; df_nprof.index.name = 'variable'
            for nlev_range1 in [[0,0],[1,3],[3,10],[10,50],[50,100],[100,200],[200,999]]:
              mycount = np.count_nonzero((ndata.transpose()>=nlev_range1[0]) & (ndata.transpose()<=nlev_range1[1]),axis=1)
              df1_nprof = pd.DataFrame({('by nlev','{0:03d}-{1:03d}'.format(*nlev_range1)):mycount})
              df1_nprof.index = ndata.columns; df1_nprof.index.name = 'variable'
              df_nprof = pd.concat([ df_nprof, df1_nprof], axis=1)
          df_ndata = pd.concat([ df_ndata, df_nprof], axis=1)
          df_ndata['UTC']=utc1
          df_ndata['date']=int(yyyymmdd_dt.strftime('%Y%m%d'))
          df_ndata = df_ndata.reset_index().set_index(['date','UTC','variable'])
          if df_summary is None:
            df_summary = df_ndata.copy()
          else:
            df_summary = pd.concat( [df_summary, df_ndata], axis=0)
      df_summary.to_pickle(pkl_sum, protocol=3, compression="gzip")
      print("saved summary in",pkl_sum)

    # POST-PROC QC
    n_levels = df_write.groupby(df_write.index)[df_write.columns[:1]].count()
    wsingle_level = np.where(n_levels==1)[0]
    if len(wsingle_level):
      print("WARNING... single level profiles...")
      print(df_write.loc[n_levels.iloc[wsingle_level].index]['timePeriod'])
    wbad_timePeriod = np.where(df_write['timePeriod'].values[:]>24575)[0]
    if len(wbad_timePeriod)>0:
      print("WARNING... bad timePeriod...")
      print(df_write.iloc[wbad_timePeriod][['station_id','time']])
      wgood_timePeriod = np.where(df_write['timePeriod'].values[:]<=24575)[0]
      df_write = df_write.iloc[wgood_timePeriod]

    # PREPARE BUFR ENCODING
    pnams = list(df_write.index.names)
    posis  = {x:pnams.index(x) for x in pnams}
    profs = df_write.index.unique()
    bufr_file_output_fmt = '/ec/res4/scratch/erc/Leo/BUFRENCODE/perm/erc/ERA6BUFR/{0}{1}'.format(bufr_output_dir,timeoffset_name)+'/{0}/{1}/{2}/{3}/{4}_prof{5}.bufr'

    # create directories, clean-up possibly pre-existing files, open new files for writing
    filearray = open_clean_bufr_files(bufr_file_output_fmt, yyyymmdd_dt, nprocesses)

    ### DEBUG ### ONLY FOR TESTING !!!!
    #w4036dd=np.where((df_cuon.index==profs[4036])& (~np.isnan(df_cuon['windDirection'].values[:])))[0]
    #df_cuon.iloc[w4036dd,list(df_cuon.columns).index('windDirectionBiasCorr')] = 90 # perform a 90-degree rotation clockwise (what came from the North now comes from the West) to see clear impact on u and v
    #w4036T=np.where((df_cuon.index==profs[4036])& (~np.isnan(df_cuon['airTemperature'].values[:])))[0]
    #df_cuon.iloc[w4036T ,list(df_cuon.columns).index('airTempBiasCorr')] = -1.55 # warm T by 1.55 Kelvin
    #w4036Td=np.where((df_cuon.index==profs[4036])& (~np.isnan(df_cuon['dewpointTemperature'].values[:])))[0]
        #df_cuon.iloc[w4036Td,list(df_cuon.columns).index('dewpointTempBiasCorr')] = 1.45 # cool Td by 1.45 Kelvin
        #w4036=np.where((df_cuon.index==profs[4036]))[0]
    #df_cuon.iloc[w4036,list(df_cuon.columns).index('height')]=1193 # assign station elevation for WMO ID 72364
    ### DEBUG ### ONLY FOR TESTING !!!!

    #for iprof1 in range(len(profs)):
    def single_call(iprof1):
      #print("run number",iprof1)
      try:
        iproc=multiPCurrentProcess()._identity[0]
      except:
        iproc=1
      prof1 = profs[iprof1]
      df_loc = df_write.loc[prof1].copy()
      obs_dt1 = df_loc['time'].iloc[0] # We pick the first time (lowest level)
      #obs_dt1 = prof1[posis['time']]
      synobs_dt = get_synoptic_date_and_time(obs_dt1)
      #report_id1 = prof1[posis['report_id']].decode('utf-8').strip() 
      source_id1 = df_loc['source_id'].values[:1].astype(bytes)[0].decode('utf-8').strip() # We pick the first source_id
      onebufr_filewb = filearray[iproc][synobs_dt.strftime('%Y%m%d%H')]
      oneprofile = {'header':{}, 'data':{}}
      oneprofile['header']['latitude'] = df_loc['lat'].values[0] # We pick the first lat/lon/elevation (lowest level)
      oneprofile['header']['longitude'] = df_loc['lon'].values[0]
      oneprofile['header']['height'] = df_loc['height'].values[0]
      oneprofile['header']['platform_type'] = df_loc['platform_type'].values[0]
      oneprofile['header']['year'] = obs_dt1.year
      oneprofile['header']['month'] = obs_dt1.month
      oneprofile['header']['day'] = obs_dt1.day
      oneprofile['header']['hour'] = obs_dt1.hour
      oneprofile['header']['minute'] = obs_dt1.minute
      oneprofile['header']['second'] = 0
      wigos_statid1_elms = prof1[posis['wigos_statid']].split('-')
      oneprofile['header']['wigosIdentifierSeries'], oneprofile['header']['wigosIssuerOfIdentifier'], oneprofile['header']['wigosIssueNumber'] = [int(x) for x in wigos_statid1_elms[:3]]
      oneprofile['header']['wigosLocalIdentifierCharacter'] = wigos_statid1_elms[3]
      # SPECIAL HANDLING OF RS IDENTIER TO ENSURE COMPATIBILITY WITH WMO IDENTIFIER
      try:
        ident_full = copy.copy(oneprofile['header']['wigosLocalIdentifierCharacter'])
        if ':' in ident_full:
          ident_full = ident_full.split(':')[-1]
        if ident_full.isdigit() and int(ident_full)<=99999:
          oneprofile['header']['stationNumber'] = int(ident_full) % 1000
          oneprofile['header']['blockNumber'] = int(ident_full) // 1000
          if int(ident_full)<10000: # format to 5-digit string
            ident_full = '{0:05d}'.format(int(ident_full))
      except:
        # we cannot encode this ill-formed WIGOS identifier, this shouldn't happen!
        print("ERROR... ILL-FORMED WIGOS IDENTIFIER...",'-'.join(wigos_statid1_elms))
        stop
      if len(ident_full)>8:
        print("WARNING ",oneprofile['header'],prof1,df_loc['station_id'].unique())
        ident_full = 'BAD'+ident_full[-5:]
      oneprofile['header']['shipOrMobileLandStationIdentifier'] = ident_full
          
      #print("encoding",'-'.join(wigos_statid1_elms),' '.join(['{0}:{1}'.format(k,oneprofile['header'][k]) for k in oneprofile['header'].keys() if k in ['wigosLocalIdentifierCharacter','stationNumber','blockNumber']]))
      sonde_type_uq = df_loc['sonde_type'].unique()
      sonde_type_uqok = []
      for sonde_type_uq1 in sonde_type_uq:
        if (type(sonde_type_uq1) is float) and (np.isnan(sonde_type_uq1)):
          continue
        sonde_type_uqok.append(sonde_type_uq1)
      if len(sonde_type_uqok)>1:
        print("ERROR ... non-unique sonde_type...")
        print("wigos_statid",wigos_statid1,"sonde_type",sonde_type_uqok)
        stop
      if len(sonde_type_uqok)==1 and type(sonde_type_uqok[0]) is bytes:
        sonde_type_utf = sonde_type_uqok[0].decode('utf-8').strip()
        try:
          if sonde_type_utf not in ['n','nan','NA']:
            sonde_type_flt = float(sonde_type_utf)
            sonde_type_int = int(sonde_type_flt)
            if (float(sonde_type_int)==sonde_type_flt) and \
               (sonde_type_int>=0) and (sonde_type_int<=255):
              oneprofile['header']['radiosondeType'] = sonde_type_int
            else:
              oneprofile['header']['sondeTypeDetail'] = int.from_bytes(sonde_type_uqok[0], 'little')
        except:
          oneprofile['header']['sondeTypeDetail'] = int.from_bytes(sonde_type_uqok[0], 'little')
        # REMINDER: TO DECODE sondeTypeDetail from a value x : int.to_bytes(int(x), length=4, byteorder='little').decode('utf-8')
        # try and map SCHROEDER description to WMO sonde type
        if ('radiosondeType' not in oneprofile['header'].keys()) and (sonde_type_utf in SCH_idx):
          oneprofile['header']['radiosondeType'] = df_WMO_SCH.loc[sonde_type_utf]['wmo_id']
          if type(oneprofile['header']['radiosondeType']) is pd.Series:
            oneprofile['header']['radiosondeType'] = oneprofile['header']['radiosondeType'].values[0]
          if 'sondeTypeDetail' in oneprofile['header'].keys():
            del oneprofile['header']['sondeTypeDetail'] # remove sondeTypeDetail if we have a match for the WMO type
      elif len(sonde_type_uqok)==1 and type(sonde_type_uqok[0]) is not bytes:
        print("UNEXPECTED SONDE_TYPE",sonde_type_uqok[0])
        stop
      for k in ['sondeTypeDetail','radiosondeType']:
        if k in oneprofile['header'].keys():
          if type(oneprofile['header'][k] is not int):
            oneprofile['header'][k] = int(oneprofile['header'][k])
      #for k in ['sondeTypeDetail','radiosondeType']:
      #  if k in oneprofile['header'].keys():
      #    print(sonde_type_uqok[0], k, oneprofile['header'][k])
      #oneprofile['header']['datasetSource'] = 'CUON   {0}'.format(int(report_id1[0]))
      oneprofile['header']['datasetSource'] = 'CUON{0:>4}'.format(source_id1[-4:])

      # HANDLING HIGH-RES SONDES - REMOVING SIGNIFICANT LEVEL MARKING
      if synobs_dt>=dt.datetime(2000,1,1): # we worry about high-resolution BUFR after this date only
        if ('timePeriod' in df_loc.columns) and ('extendedVerticalSoundingSignificance' in df_loc.columns):
          wts = np.where((~np.isnan(df_loc['timePeriod'].values[:])) & (~np.isnan(df_loc['extendedVerticalSoundingSignificance'].values[:])))[0]
          if len(wts)>50: # for fewer than 50 levels we do not bother
            times = np.unique(np.sort(df_loc['timePeriod'].values[wts]))
            deltat = (times[1:] -times[:-1])
            try:
              dt95 = np.percentile(deltat,95)
            except:
              dt95 = np.nan
              print("error when trying to get 95th percentile from",deltat)
            if (~np.isnan(dt95)) and (dt95 <= 30): # we set our threshold at 30 seconds -- this is approx 150 meters, at 5 m/s ascent speed
              print("high-res RS",oneprofile['header'])
              icol_extendedVerticalSoundingSignificance = list(df_loc.columns).index('extendedVerticalSoundingSignificance')
              for signifbit in [18-5,18-7]:
                wsignif = np.where((~np.isnan(df_loc['extendedVerticalSoundingSignificance'].values[:])) & (np.bitwise_and(df_loc['extendedVerticalSoundingSignificance'].values[:].astype(int),2**signifbit)>0))[0]
                if len(wsignif)>0:
                  df_loc.iloc[wsignif, icol_extendedVerticalSoundingSignificance] -= 2**signifbit
                  print("removed",len(wsignif),"bits for",signifbit)

      # MARKING SURFACE LEVELS
      if ('nonCoordinateGeopotentialHeight' in df_loc.columns) and (~np.isnan(oneprofile['header']['height'])) and ('extendedVerticalSoundingSignificance' in df_loc.columns):
        wZsfc = np.where((~np.isnan(df_loc['nonCoordinateGeopotentialHeight'].values[:])) & (~np.isnan(df_loc['pressure'].values[:])) & (np.abs(df_loc['nonCoordinateGeopotentialHeight'].values[:]/9.80665-oneprofile['header']['height'])<=10.))[0]
        if len(wZsfc)>0:
          icol_extendedVerticalSoundingSignificance = list(df_loc.columns).index('extendedVerticalSoundingSignificance')
          # pick the lowest level (highest pressure)
          isfc = np.argmax(df_loc['pressure'].values[wZsfc])
          df_loc.iloc[wZsfc[isfc], icol_extendedVerticalSoundingSignificance] += 2**(18-1)

      # STATIONS WITHIN 1 deg of NORTH OR SOUTH POLE AND REPORTING WIND
      if (abs(oneprofile['header']['latitude'])>=89.) and ('windDirection' in df_loc.columns) and ('windSpeed' in df_loc.columns) and ('latitudeDisplacement' in df_loc.columns) and ('longitudeDisplacement' in df_loc.columns):
        if np.count_nonzero((~np.isnan(df_loc['windDirection'].values[:])) | (~np.isnan(df_loc['windSpeed'].values[:])))>0:
          df_loc['latitudeDisplacement'] = np.nan # we cannot keep displacements for these locations because the wind should have been remapped
          df_loc['longitudeDisplacement'] = np.nan
          print("removed lat/lon displacement for station",oneprofile['header'])

      #oneprofile={'header':{'longitude':337.4-360., 'latitude':63.97, 'height':425.0, \
      #              'year':1997, 'month':7, 'day':1, 'hour':0, 'minute':0, 'second':0, \
      #              'wigosIdentifierSeries':0, 'wigosIssuerOfIdentifier':20000, 'wigosIssueNumber':0, 'wigosLocalIdentifierCharacter':'4018', 'blockNumber':4, 'stationNumber':18, 'shipOrMobileLandStationIdentifier':''}, \
      # 'data':{'pressure':[x[1] for x in T_P],'airTemperature':[x[0] for x in T_P],'airTemperatureCorrection':[0.01 for x in T_P]}}
      for k in [('pressure',float), \
                ('airTemperature',float), \
                ('airTempBiasCorr',float), \
                ('dewpointTemperature',float), \
                ('dewpointTempBiasCorr',float), \
                ('windDirection',int), \
                ('windDirectionBiasCorr',int), \
                ('windSpeed',float), \
                ('windSpeedBiasCorr',float), \
#('windSpeedCorrApplied',float), \
                ('nonCoordinateGeopotentialHeight',float), \
                ('extendedVerticalSoundingSignificance',int), \
                ('timePeriod',int), \
                ('latitudeDisplacement',float), \
                ('longitudeDisplacement',float), \
               ]:
        #newk0 = copy.copy(k[0]) # pre-WMOBUFR
        # WMOBUFR:
        if k[0] in ['airTempBiasCorr','dewpointTempBiasCorr','windDirectionBiasCorr','windSpeedBiasCorr']:
          newk0 = {'airTempBiasCorr':'airTemperature->differenceStatisticalValue', \
                   'dewpointTempBiasCorr':'dewpointTemperature->differenceStatisticalValue', \
                   'windDirectionBiasCorr':'windDirection->differenceStatisticalValue', \
                   'windSpeedBiasCorr':'windSpeed->differenceStatisticalValue'}[k[0]]
        else:
          newk0 = copy.copy(k[0])
        if k[0] in df_loc.columns:
          oneprofile['data'][newk0] = df_loc[k[0]].values[:]
        #if k=='windSpeedBiasCorr':
        #  print("FOUND WIND SPEED BIAS ...")
        #  stop

      #try:
      #print(oneprofile)
      bufr_encode(oneprofile, onebufr_filewb)
      #except CodesInternalError as err:
      #  traceback.print_exc(file=sys.stderr)
      #  #return 1


    # parallel processing
    ncalls=len(profs)
    print(ncalls,"to do")
    ##for iprof1 in range(len(profs)): # DEBUG
    ##  single_call(iprof1) # DEBUG
    with multiPPool(nprocesses) as p:
      p.map(single_call, list(range(ncalls)))
      p.close(); p.join()
      print("All tasks completed", flush=True)
    # close all files
    close_clean_bufr_files(filearray)

#if __name__ == "__main__":
#    sys.exit(main())
