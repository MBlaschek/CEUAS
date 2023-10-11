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
import datetime
import xarray

def to_days_since(time):
    reference = datetime.datetime(1900, 1, 1)
    fstart = datetime.datetime(int(str(time)[:4]),int(str(time)[-2:]),15)
    offset = fstart - reference
    return offset.days

def make_nc_file(name):
    
    statid = name
    
    day_ref = glob.glob("./rttov_out_ih/"+statid+"/IGRAH*"+statid+"*day_refl*")
    if len(day_ref) == 0:
        return
    day_re = pickle.load( open( day_ref[0], "rb" ) )
    
    day_dat = glob.glob("./rttov_out_ih/"+statid+"/IGRAH*"+statid+"*day_dates*")
    day_da = pickle.load( open( day_dat[0], "rb" ) )
    
    night_ref = glob.glob("./rttov_out_ih/"+statid+"/IGRAH*"+statid+"*night_refl*")
    if len(night_ref) == 0:
        return
    night_re = pickle.load( open( night_ref[0], "rb" ) )
    
    night_dat = glob.glob("./rttov_out_ih/"+statid+"/IGRAH*"+statid+"*night_dates*")
    night_da = pickle.load( open( night_dat[0], "rb" ) )
    
    
    fn = './rttov_nc/IH_rttov_'+name+'.nc'
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
    ds.references = "www.univie.ac.at/theoret-met/research/IH" ;
    ds.Stationname = str(statid) ;
                
    
    # populating Variables:
    dd = []
    nd = []
    for i in day_da:
        dd.append(to_days_since(i))
        
    for i in night_da:
        nd.append(to_days_since(i))
    
    datum[:] = dd
    pressure = [2,3,4]
    
    checkshape = np.shape([0,0,0])
    day_re_new = []
    for i in day_re:
        if not isinstance(i, float):
            if np.shape(i[0]) != checkshape:
                day_re_new.append(i)
            else:
                day_re_new.append(i[0])
        else:
            day_re_new.append(np.array([np.nan, np.nan, np.nan]))
    night_re_new = []
    for i in night_re:
        if not isinstance(i, float):
            if np.shape(i[0]) != checkshape:
                night_re_new.append(i)
            else:
                night_re_new.append(i[0])
        else:
            night_re_new.append(np.array([np.nan, np.nan, np.nan]))

    montemp[:,:,:] = [np.transpose(night_re_new),np.transpose(day_re_new)]
    press[:] = [2,3,4]
    
    o_file = glob.glob('igra_h_single_stats_new/*'+ statid +'*.nc')
    df = xarray.open_dataset(o_file[0]).to_dataframe()
    
    lat[:] = df.lat.iloc[0]
    lon[:] = df.lon.iloc[0]
    ds.close()
    

if __name__ == '__main__':
    
    statids = []
    stats = glob.glob('./rttov_out_ih/*')
    for i in stats:
        statid = (i.split('/')[-1])
        statids.append(statid)
    pool = multiprocessing.Pool(processes=15)
    func=partial(make_nc_file)
    result_list = list(pool.map(func, statids))
    print(result_list)