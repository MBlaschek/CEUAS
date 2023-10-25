# System information
import numpy as np
import os, sys, glob
import xarray
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


def to_seconds_since(time):
    """ return seconds since a Reference date applying time unit

    Args:
        time_elapsed: seconds or other time unit
        time_unit: 'seconds since 1900-01-01 00:00:00'
        reference: '1900-01-01 00:00:00'

    Returns:
        int : seconds since
    """
    reference = datetime.datetime(1900, 1, 1)
    date = datetime.datetime.strptime(str(time), '%Y-%m-%d')# %H:%M:%S')
    offset = (date-reference)
    return (offset.days * 24 * 60 * 60 + offset.seconds)

# def single_stat(i):
#     print(i)
#     stats = glob.glob('./igra_stats/'+i+'*.nc')
#     ds = xarray.open_dataset(stats[0], engine="netcdf4")
#     dafr = ds.to_dataframe()
#     for j in stats[1:]:
#         ds = xarray.open_dataset(j, engine="netcdf4")
#         df = ds.to_dataframe()
#         dafr = dafr.append(df)
#     dafr = dafr.sort_values(by=['report_timestamp', 'air_pressure'])
#     dafr = dafr.drop(columns=['index'])
#     rt = []
#     at = []
#     for i in range(len(dafr)):
#         at.append(to_seconds_since(dafr.actual_time[i]))
#         rt.append(to_seconds_since(dafr.index[i]))
#     dafr.index = rt
#     dafr.actual_time = at
#     ds = dafr.to_xarray()
#     ds.to_netcdf('./igra_single_stats/'+i+'.nc')
def single_stat(i):
    name = i
    stats = glob.glob('./igra_stats/'+i+'*.nc')
    ds = xarray.open_dataset(stats[0], engine="netcdf4")
    dafr = ds.to_dataframe()
    for j in stats[1:]:
        ds = xarray.open_dataset(j, engine="netcdf4")
        df = ds.to_dataframe()
        dafr = dafr.append(df)
    dafr = dafr.sort_values(by=['report_timestamp', 'air_pressure'])
    dafr = dafr.drop(columns=['index'])
    dafr = dafr[dafr.air_pressure.isin([1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000,92500,100000])]
    
#     yrs = glob.glob('/users/staff/a1400070/CEUAS/CEUAS/public/intercomparisons/IGRA_downloaded/*_H_*/*')
#     yrs.sort()
#     for i in range(len(yrs[:])):
#         startt = time.time()
#         print(yrs[i])
#         try:
#             ds = pd.read_csv(yrs[i], header=16)
#             if i == 0:
#                 ds70308 = ds[ds.station_name == 'USM00070308']
#                 ds70316 = ds[ds.station_name == 'USM00070316']
#                 ds70350 = ds[ds.station_name == 'USM00070350']
#                 ds91408 = ds[ds.station_name == 'PSM00091408']
#             else:
#                 ds70308 = ds70308.append(ds[ds.station_name == 'USM00070308'])
#                 ds70316 = ds70316.append(ds[ds.station_name == 'USM00070316'])
#                 ds70350 = ds70350.append(ds[ds.station_name == 'USM00070350'])
#                 ds91408 = ds91408.append(ds[ds.station_name == 'PSM00091408'])
#         except:
#             ds = pd.read_csv(yrs[i], header=15)
#             if i == 0:
#                 ds70308 = ds[ds.station_name == 'USM00070308']
#                 ds70316 = ds[ds.station_name == 'USM00070316']
#                 ds70350 = ds[ds.station_name == 'USM00070350']
#                 ds91408 = ds[ds.station_name == 'PSM00091408']
#             else:
#                 ds70308 = ds70308.append(ds[ds.station_name == 'USM00070308'])
#                 ds70316 = ds70316.append(ds[ds.station_name == 'USM00070316'])
#                 ds70350 = ds70350.append(ds[ds.station_name == 'USM00070350'])
#                 ds91408 = ds91408.append(ds[ds.station_name == 'PSM00091408'])
#         print(time.time()-startt)


    daterange = pd.date_range(start=dafr.index[0],end=dafr.index[-1]).to_pydatetime().tolist()
    pressurerange = [1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000,92500,100000]
    dafr_day = dafr[dafr.index.str[11:13] > '06']
    dafr_day = dafr_day[dafr_day.index.str[11:13] <= '18']
    dafr_night = dafr[dafr.index.str[11:13] > '18']
    dafr_night = dafr_night.append(dafr[dafr.index.str[11:13] <= '06'])
    as_list = dafr_day.index.tolist()
    for i in range(len(as_list)):
        if as_list[i][11:13] != "12":
            as_list[i] = as_list[i][:11]+"12"+as_list[i][13:]
    dafr_day.index = as_list
    as_list = dafr_night.index.tolist()
    for i in range(len(as_list)):
        if as_list[i][11:13] != "00":
            as_list[i] = as_list[i][:11]+"00"+as_list[i][13:]
    dafr_night.index = as_list
    dafr_day.index = pd.to_datetime(dafr_day.index ,format= '%Y-%m-%d %H:%M:%S+%f')
    dafr_night.index = pd.to_datetime(dafr_night.index ,format= '%Y-%m-%d %H:%M:%S+%f')

    night_datearray = []
    day_datearray = []
    dayadd = pd.Timedelta(hours =  12)
    for k in daterange:
        night_datearray.append(pd.Timestamp(k).tz_localize(None))
        day_datearray.append(pd.Timestamp(k).tz_localize(None) + dayadd)

    ta0 = []
    rh0 = []
    ew0 = []
    nw0 = []
    for j in pressurerange:
        print(j)
        fill_df = dafr_night[dafr_night.air_pressure == j]
        fill_df = fill_df.loc[~fill_df.index.duplicated(), :]
        final_df = fill_df.reindex(night_datearray, fill_value=np.nan)
        ta0.append(np.array(final_df.air_temperature))
        rh0.append(np.array(final_df.relative_humidity))
        ew0.append(np.array(final_df.eastward_wind_component))
        nw0.append(np.array(final_df.northward_wind_component))

    ta1 = []
    rh1 = []
    ew1 = []
    nw1 = []
    for j in pressurerange:
        print(j)
        fill_df = dafr_day[dafr_day.air_pressure == j]
        fill_df = fill_df.loc[~fill_df.index.duplicated(), :]
        final_df = fill_df.reindex(day_datearray, fill_value=np.nan)
        ta1.append(np.array(final_df.air_temperature))
        rh1.append(np.array(final_df.relative_humidity))
        ew1.append(np.array(final_df.eastward_wind_component))
        nw1.append(np.array(final_df.northward_wind_component))
    #    
    fn = './igra_h_single_stats_new/igra_h_'+name+'.nc'
    ds = nc.Dataset(fn, 'w', format='NETCDF4')

    # Dimensions:
    station = ds.createDimension('station', 1)
    time = ds.createDimension('time', len(night_datearray))
    pressure = ds.createDimension('pressure', 16)
    hour = ds.createDimension('hour', 2)

    # Variables:
    lat = ds.createVariable('lat', 'f4', ('station',))
    lon = ds.createVariable('lon', 'f4', ('station',))
    press = ds.createVariable('press', 'f4', ('pressure',))
    ta = ds.createVariable('ta', 'f4', ('hour', 'pressure', 'time'))
    rh = ds.createVariable('rh', 'f4', ('hour', 'pressure', 'time'))
    v = ds.createVariable('v', 'f4', ('hour', 'pressure', 'time'))
    u = ds.createVariable('u', 'f4', ('hour', 'pressure', 'time'))
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

    press.long_name = "pressure layer"
    press.units = "Pa"
    press.axis = "Z"
    press.valid_range = [1000, 100000]
    press.missing_value = -999.

    datum.long_name = "datum"
    datum.units = "seconds since 1900-01-01 0:0:0"
    datum.axis = "T"
    datum.calendar = "gregorian"
    datum.missing_value = -999.

    ta.long_name = "temperature"
    ta.units = "K"
    ta.missing_value = -999.
    ta.valid_range = [0., 400.]
    ta.cell_methods = "time"

    rh.long_name = "relative_humidity"
    rh.units = "1"
    rh.missing_value = -999.
    rh.valid_range = [0., 400.]
    rh.cell_methods = "time"

    u.long_name = "eastward_wind"
    u.units = "m/s"
    u.missing_value = -999.
    u.valid_range = [0., 400.]
    u.cell_methods = "time"

    v.long_name = "northward_wind"
    v.units = "m/s"
    v.missing_value = -999.
    v.valid_range = [0., 400.]
    v.cell_methods = "time"


    # Global Attributes:
    ds.Conventions = "CF-1.1"
    ds.title = ""
    ds.institution = "University of Vienna"
    ds.history = "17/09/20"
    ds.source = "IGRA_H"
    ds.Stationname = str(name)


    # populating Variables:

    dd = []
    for o in day_datearray:
        dd.append(to_seconds_since(str(o.date())))

    datum[:] = dd
    press[:] = [1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000,92500,100000]

    ta[:,:,:] = [ta0,ta1]
    rh[:,:,:] = [rh0,rh1]
    u[:,:,:] = [ew0,ew1]
    v[:,:,:] = [nw0,nw1]

    lat[:] = dafr.latitude.iloc[0]
    lon[:] = dafr.longitude.iloc[0]

    ds.close()
    
if __name__ == '__main__': 
    
    done = glob.glob('./igra_stats/*')
    stdone = []
    for i in done:
        st = i.split('/')[-1].split('_')[0]
        if not st in stdone:
            if not os.path.isfile('./igra_h_single_stats_new/igra_h_'+st+'.nc'):
                stdone.append(st)
                
#     stdone = ["GMM00010393"]
    pool = multiprocessing.Pool(processes=10)
    func=partial(single_stat)
    result_list = list(pool.map(func, stdone))
    print(result_list)