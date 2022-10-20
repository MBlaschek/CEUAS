#!/usr/bin/env
# coding: utf-8

import numpy as np
import os, glob, sys
import netCDF4
import xarray as xr
import time

import ray
ray.init(num_cpus=12)


@ray.remote
def process_file(fdir):
#     fdir = './era_hourly/era_land_hourly2000_1.nc'
    file = xr.open_dataset(fdir, decode_cf=False)
    dt = fdir.split('hourly')[-1].split('.nc')[0]
    yr = dt.split('_')[0]
    mn = dt.split('_')[1].zfill(2)
    
    intervall = 30
    try:
        for hour in [0,12]:

            if len(glob.glob('/users/staff/a1400070/scratch/era_land/era5_land_monthly_mean_surface_'+yr+'_'+mn+'_'+str(hour).zfill(2)+'.nc')) == 1:
                continue

            first = True
            for lat in range(-90,90,intervall):
                print('lat: ', lat, lat+intervall)
                for lon in range(0,360,intervall):
                        ds = file.sel(latitude=slice(lat+intervall,lat), longitude=slice(lon,lon+intervall))
                        out = xr.decode_cf(ds.where(ds['time']%24 == hour, drop=True))
                        out = out.resample(time="1MS").mean(dim="time")
        #                 out = out.assign_coords({'time': out_ds.time + hour*60*60*1000000000})
                        if first:
                            out_ds = out
                            first = False
                        else:
            #                 out_ds = xr.combine_by_coords([out_ds, out])
                            out_ds = xr.merge([out_ds, out])
            out_ds.to_netcdf('/users/staff/a1400070/scratch/era_land/era5_land_monthly_mean_surface_'+yr+'_'+mn+'_'+str(hour).zfill(2)+'.nc')
        return yr
    except:
        return 'ERROR: '+str(yr)
        
        
if __name__ == '__main__':
    file_list = []
    for i in glob.glob('./era_hourly/era_land_hourly*2022*.nc')[:]:
        file_list.append(process_file.remote(i))
    results = ray.get(file_list)
    print(results)
