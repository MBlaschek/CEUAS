# System information
import numpy as np
import os, sys, glob
print(sys.executable)
print(sys.version_info)
import pandas as pd
pd.set_option('display.max_columns', None)
# pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)
sys.path.append(os.getcwd()+'/../cds-backend/code/')
import cds_eua3 as eua
import pickle
import multiprocessing
from functools import partial
import netCDF4 as nc

def make_nc_file(name):
    
    statid = name[-5:]
    
    day_ref = glob.glob("./rttov_out/"+statid+"/*day_refl*")
    if len(day_ref) == 0:
        return
    day_re = pickle.load( open( day_ref[0], "rb" ) )
    
    day_dat = glob.glob("./rttov_out/"+statid+"/*day_dates*")
    day_da = pickle.load( open( day_dat[0], "rb" ) )
    
    night_ref = glob.glob("./rttov_out/"+statid+"/*night_refl*")
    if len(night_ref) == 0:
        return
    night_re = pickle.load( open( night_ref[0], "rb" ) )
    
    night_dat = glob.glob("./rttov_out/"+statid+"/*night_dates*")
    night_da = pickle.load( open( night_dat[0], "rb" ) )
            
            
    fn = './rttov_nc/rttov_'+name+'.nc'
    ds = nc.Dataset(fn, 'w', format='NETCDF4')
    
    # Dimensions:
    station = ds.createDimension('station', 1)
    numdat = ds.createDimension('numdat', 1)
    time = ds.createDimension('time', len(day_da))
    pressure = ds.createDimension('pressure', 3)
    hour = ds.createDimension('hour', 2)

    # Variables:
    lat = ds.createVariable('lat', 'f4', ('station',))
    lon = ds.createVariable('lon', 'f4', ('station',))
    press = ds.createVariable('press', 'f4', ('pressure',))
    montemp = ds.createVariable('montemp', 'f4', ('hour', 'pressure', 'time'))
    datum = ds.createVariable('datum', 'f4', ('hour', 'time',))
    
    # Attributes:

    lat.long_name = "station latitude"
    lat.units = "degrees_north"
    lat.axis = "Y"
    lat.valid_range = [-90.,90.]
    lat.missing_value = -999.
    
    lon.long_name = "station longitude"
    lon.units = "degrees_east"
    lon.axis = "X"
    lon.valid_range = [-180., 180.]
    lon.missing_value = -999.
    
    press.long_name = "MSU layers"
    press.units = ""
    press.axis = "Z"
    press.valid_range = [1, 4]
    press.missing_value = -999.

    datum.long_name = "datum"
    datum.units = "days since 1900-01-01 0:0:0"
    datum.axis = "T"
    datum.calendar = "gregorian"
    datum.missing_value = -999.
    
    montemp.long_name = "monthly_uncorrected_bt"
    montemp.units = "K"
    montemp.missing_value = -999.
    montemp.valid_range = [0., 400.]
    montemp.cell_methods = "time: monthly means"
    
    # Global Attributes:
    ds.Conventions = "CF-1.1" ;
    ds.title = "" ;
    ds.institution = "University of Vienna" ;
    ds.history = "17/09/20" ;
    ds.source = "ERA5, CUON" ;
    ds.references = "www.univie.ac.at/theoret-met/research/raobcore" ;
    ds.Stationname = str(statid) ;
                
    
    # populating Variables:
    
    datum[:] = [day_da,night_da]
    pressure = [2,3,4]
    
    day_re_new = []
    for i in day_re:
        if not isinstance(i, float):
            day_re_new.append(i[0])
        else:
            day_re_new.append(np.array([np.nan, np.nan, np.nan]))
    night_re_new = []
    for i in night_re:
        if not isinstance(i, float):
            night_re_new.append(i[0])
        else:
            night_re_new.append(np.array([np.nan, np.nan, np.nan]))

    montemp[:,:,:] = [np.transpose(day_re_new),np.transpose(night_re_new)]
    press[:] = [2,3,4]
    
    o_file = glob.glob('/mnt/users/scratch/leo/scratch/converted_v7/*'+statid+'*_CEUAS_merged_v1.nc')
    df = eua.CDMDataset(filename = o_file[0]).to_dataframe(groups=['observations_table'], variables=['latitude', 'longitude'])
    
    lat[:] = df.latitude.iloc[0]
    lon[:] = df.longitude.iloc[0]
    ds.close()
    

if __name__ == '__main__':
    
    statids = []
    stats = glob.glob('/mnt/users/scratch/leo/scratch/converted_v7/*_CEUAS_merged_v1.nc')
    for i in stats:
        statid = (i.split('/')[-1].split('_')[0])
        if not os.path.isfile('./rttov_nc/rttov_'+statid+'.nc'):
            statids.append(statid)
        else:
            print(statid)
    print(statids)
    for i in statids:
        try:
            make_nc_file(i)
        except:
            pass