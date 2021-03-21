""" Example to analyze Desroziers uncertainty from file """

import os,sys
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker


pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)

""" Test file """
f = '/raid60/scratch/federico/DESROZIERS_MARCH2021/0-20000-0-82930_CEUAS_merged_v0.nc'
station = f.split('/')[-1].split('_CEUAS')[0]

# NB do not use with very large station! Requires a lot of time to create dataframe 
df_a    = xr.open_dataset(f, engine = 'h5netcdf' , group = 'advanced_uncertainty', decode_times = True).to_dataframe()
df_obs = xr.open_dataset(f, engine = 'h5netcdf' , group = 'observations_table', decode_times = True).to_dataframe()

df_a['date_time'] = df_obs['date_time']

""" Selecting the variable temperature """
ind = np.where (df_obs['observed_variable'][:]==85)[0]
DF = df_a.iloc[ind]


""" Selecting a date """
date = '1986-08-22T12:00:00'
dt = np.datetime64(date)

DFF = DF.loc [DF['date_time'] == dt ]


def plot_profile(df, date, station ):
    
    std_plevs    = [1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 85000, 92500, 100000]
        
    temp, bias, des_30, press, adj = [], [], [], [], []
    des_30_p, des_30_m = [], []
    
    for p in std_plevs:
        ind = np.where ( df['z_coordinate'] == p )[0]
        if len(ind) > 0:
            
            t = df['observation_value'][ind].values[0]
            b = df['biascorr@body'][ind].values[0]
            d30 = df['desroziers_30'][ind].values[0]
            
            if not np.isnan(t) and not np.isnan(b) and not np.isnan(d30):
                temp.append(t )
                bias.append(b)
                des_30.append(d30)
                press.append(p/100)
                a = t-b
                adj.append(a)
                des_30_m.append(a - d30/2)
                des_30_p.append(a + d30/2)
                
    
    fs = 20
    fig, ax1= plt.subplots(figsize=(12,10) )        
          
    fig.suptitle('Station 0-20000-0-82390 - Profile '  + date , y = 0.94, fontsize = fs)
          
    ax1.tick_params(axis='both', which='major', labelsize=15)
    ax1.tick_params(axis='both', which='minor', labelsize=8)
    w = ax1.invert_yaxis() 
                  
    """ Adding error region """
    ax1.fill_betweenx(press, des_30_m, des_30_p, color = 'lime', alpha = 0.2)    
    ax1.plot(adj, press , label = 'Desroziers Unc. Band' , color = 'lime' , alpha = 0.2  )
                  
    ax1.set_ylabel( 'Pressure [hPa]'   , fontsize = fs )     
    ax1.set_xlabel( 'Temperature [K]' , fontsize = fs )          
    
    #ax1.plot(temp, press , label = 'Air Temperature' , color = 'red'  )
    ax1.scatter(temp, press, color = 'red'  , label = 'Air Temperature' )
    
    #ax1.plot(temp, press , label = 'Air Temperature Bias Adj.' , color = 'orange'  )
    ax1.scatter(adj, press, color = 'orange' , label = 'Air Temperature Bias Adj.' )
    
    
    #ax1.plot(des_30_m, press , label = 'Desroziers Uncertainty' , color = 'orange' , ls = '--')
    #ax1.plot(des_30_p, press , label = 'Desroziers Uncertainty' , color = 'orange' , ls = '--')
    
    
    ax1.grid(ls =":" , color = "lightgray")
    ax1.legend(fontsize = fs)
    os.system('mkdir Plots')
    plt.savefig('Plots/Desrozier_profile_' + date + '.png', dpi = 150 )
    
    
    
dummy = plot_profile(DFF, date, station)
0