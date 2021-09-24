import os,sys
import h5py as h5
import pandas as pd
import netCDF4 as nc

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker




"""
Reading the homogenized dataset Zhou at Al. 2020

The variable time has lenght:
   time(752)
since the paper describes the data between 1958-2018, but the files contain also data for 2019 and up until August 2020
so that 752 = 61*12 + (12+8)

Raw temperature shape
f_HOMO["rawT"].shape
(1184, 752, 15, 2)

"""

def get_SUNY_data(var, hour, plevel, file):
    '''
    Retrieves the temperature time series from a subdaily file e.g. homo-raw-subdaily-station/AUM00011035.nc
    (shape rawT(time, pressure, UTC) )
    NB only temperatures are available
            Parameters:
                    plevel (int):  index of the desired plevel, see below for mapping pascal-index

                                        Pressure: [  1000,   2000,   3000,   5000,   7000,  10000,  15000,
                                        20000,  25000,  30000,  40000,  50000,  70000,  85000,
                                        92500, 100000, 101300 ]
                                        
                                        50000: index 11
                                        70000: index 12
            
                    hour (int): index of the desired hour, see below for mapping hour-index
                    
                                    UTC: [0,12]
                                    0  : index 0
                                    12: index 1
                                    
                    var (int or str):variable number according to CDM convention (85: temp, 38: relative hum, 107: wind speed)
    '''
    
    file_data = nc.Dataset(file)        
    
    if var == 85 or var == '85':
        var_suny = 'rawT'
        
    if hour ==12:
        hour = 1
    pressure = file_data['pressure'][:]
    if plevel in pressure:
        
        index = np.where ( pressure == plevel )[0][0]
        
        obs = file_data[var_suny][:,index, hour].data
        times = file_data['time'][:]
        dt = pd.to_datetime( times, format='%Y%m%d' )
        df = pd.DataFrame( {"obs":obs, "time":dt, 'month': dt.month, 'year': dt.year} )
        # removing missing temperature values flagged with -9999.0
        cleaned_df = df.loc[df['obs'] > -999 ]
        return cleaned_df
    
    else:
        return 0
    



def get_monthly_counts(df):
    ''' Extract the monthly counts 
                Parameters:
                    df (pandas DataFrame): data 
    '''
    years, months, times, counts = [],[],[],[]
    
    for y in range(1900,2021):
        df_y = df.loc[ df['year']==y ]
        
        for m in range(1,13):
                            
            c = len (df_y.loc[ df_y['month'] == m ] )
    
            years.append(y)
            months.append(m)
            times.append(str(y)+str(m))
            counts.append(c)
            
            time_c = pd.to_datetime( times, format='%Y%m' )
            
    df_counts = pd.DataFrame(  { "year": years, "month": months, "time": time_c, "counts": counts } )    
    return df_counts

    
                     
                     
                     
    
    
    
def get_CUON_data(var, plev, hour , merged, file_data):
    f = h5.File(file_data, 'r')
    
    # if merged, will use sorted recordindex 
        
    if merged:
        if not str(var) in f["recordindices"].keys():
            print("No variable ", var, " for this file! ", file_data )            
            return 0
        
        else:
            ind = f["recordindices"][str(var)]
            imin, imax = ind[0], ind[-1]
    
            dt = f["observations_table"]["date_time"][imin:imax]
            red_df = pd.DataFrame(  { 'z':f["observations_table"]["z_coordinate"][imin:imax] }  )
            dt = pd.to_datetime(dt,  unit='s', origin=pd.Timestamp('1900-01-01'))
            red_df["month"] = dt.month
            red_df["year"] = dt.year        
            red_df["hour"] = dt.hour        
                
            red_df["time"] = dt               
            red_df["day"] = dt.day        
    
    else:
        dt = f["observations_table"]["date_time"]
        
        red_df = pd.DataFrame(  { 'z': f["observations_table"]["z_coordinate"], 
                                                    'observed_variable':  f["observations_table"]["observed_variable"],
                                                    'time': f["observations_table"]["date_time"] }   )

        red_df = red_df.loc[ (red_df["observed_variable"] == var) &  (red_df["z"] == plev) ]
        
        dt = red_df['time']
        
        dt = pd.to_datetime( dt,  unit='s', origin=pd.Timestamp('1900-01-01') )
        
        red_df["month"] = dt.dt.month
        red_df["year"] = dt.dt.year        
        red_df["hour"] = dt.dt.hour        
        red_df["day"] = dt.dt.day        
            
    red_df["time"] = dt        
            
    # Rounding to midnight/midday
    red_df = red_df.replace(   {'hour': {21:0, 22:0, 23: 0, 1:0, 2:0 , 3:0, 9:12, 10:12, 11:12, 13:12, 14:12, 15:12 }}    )
            
    red_df = red_df.loc[  red_df["hour"] == hour ]
    
    red_df = red_df.dropna()
    red_df = red_df.drop_duplicates(subset=['year', 'month', 'day','hour'], keep='last')
    
    return red_df



