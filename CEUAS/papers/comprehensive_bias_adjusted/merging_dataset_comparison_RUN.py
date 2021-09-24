import os,sys
import h5py as h5
import pandas as pd
import netCDF4 as nc

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker

from merging_dataset_comparison_UTILS import *


from multiprocessing import Pool
from functools  import partial


# to be imported:
#  get_SUNY_data(file_sub, hour, plevel)
#  get_monthly_counts(df)
#  get_CUON_data(file, var=85, merged=True, plev='', hour='' )



"""
# Example CUON station
cuon_vienna = '/mnt/users/scratch/leo/scratch/converted_v7/0-20001-0-11035_CEUAS_merged_v1.nc'
# Example SUNY station
subdaily_vienna = "/users/staff/leo/fastscratch/SUNY/homo-raw-subdaily-station/AUM00011035.nc"
subdaily_random = "/users/staff/leo/fastscratch/SUNY/homo-raw-subdaily-station/RSM00032215.nc"
#Example NCAR station
ncar_vienna = '/scratch/das/federico/MAY2021_HARVEST_secondary/ncar/0-20001-0-11035_ncar_harvested_uadb_trhc_11035.txt.nc' 
#Example IGRA station
igra_vienna = '/scratch/das/federico/MAY2021_HARVEST_secondary/igra2/0-20001-0-11035_igra2_harvested_AUM00011035-data.txt.nc' 
"""



# Some wrapping functions for multiprocesses and combination of output 

def wrapper( ds, var, plevel, hour, file):
    
    try:
        
        if 'SUNY' in ds:
            data = get_SUNY_data(var, hour, plevel, file)
        elif 'CUON' in ds:
            data = get_CUON_data(var, plevel, hour , True, file)
        else:
            data = get_CUON_data(var, plevel, hour , False, file)        
            
        if isinstance(data, pd.DataFrame):
            counts = get_monthly_counts(data)    
            return counts
        
    except:
        print("Fail !!! " , file )
        



def sum_counts(dic_df, dataset = "", var="", plevel="", hour=""):
    
    if dic_df:
        all_counts = {}
        for df in dic_df:
            if isinstance(df, pd.DataFrame):
                
                all_dt = df['time']
                for t,num in zip(all_dt, range(len(all_dt) ) ):
                    try:  
                        all_counts[t] = df['counts'][num]  + all_counts[t]
                    except:
                        all_counts[t] =  df['counts'][num] 
                        
        
            dates,counts = [],[]
            for dt in all_counts.keys():
                dates.append(dt)
                counts.append(all_counts[dt] )
        
            df = pd.DataFrame( { 'date':dates, 'counts':counts} ).to_csv( 'data/' + dataset +'_' + '_' + str(var) + '_' + str(plevel) + '_' + str(hour) + '_monthly_counts.csv' , sep = '\t')

    else:
        pass 


dbs = { 
                  'IGRA2': "/scratch/das/federico/MAY2021_HARVEST_secondary/igra2/" ,
            }

# Location of datasets 
dbs = { 'SUNY': "/users/staff/leo/fastscratch/SUNY/homo-raw-subdaily-station/" ,
              'IGRA2': "/scratch/das/federico/MAY2021_HARVEST_secondary/igra2/" ,
              'NCAR': "/scratch/das/federico/MAY2021_HARVEST_secondary/ncar/" ,
              'CUON': "/mnt/users/scratch/leo/scratch/converted_v7/" ,
              'ERA5_1': "//scratch/das/federico/MAY2021_HARVEST_secondary/era5_1//" ,
              'ERA5_2': "//scratch/das/federico/MAY2021_HARVEST_secondary/era5_2//" ,
              
        }





POOL = False
POOL = True
pool = Pool(30)

"""
# Loop over datasets 
for db in dbs.keys():  
   #files =  [dbs[db] + '/' + f  for f in os.listdir( dbs[db] ) if '.nc' in f  and '11035' in f ]
    files =  [dbs[db] + '/' + f  for f in os.listdir( dbs[db] ) if '.nc' in f  ]
    
    if db == 'NCAR':
        files  = [f for f in files if 'wind' not in f ]
    h = 0
    for v in [85]:
        if db == "SUNY" and v != 85:
            print("Skipping SUNY since no is data available except temperature")
            continue
        
        for p in [10000,50000,70000,85000]:
            print("+++++ Doing " , db , ' ' , h , ' ' , v , ' ' , p )
            
            if POOL:
                func = partial(wrapper, db, v, p, h )
                out = pool.map(func, files)
                                
                print("+++++ Done " , db , ' ' , h , ' ' , v , ' ' , p )
            else:
                out = []
                for file in files:
                    print('File: ', file )
                    counts_df = wrapper( db, v, p, h, file)
                    out.append(counts_df)
                    
            print("Aggregating results")
            c = sum_counts(out, dataset = db, var=v, plevel=p, hour=h)
 """
 
 
 

def plot(p, h):
    
    fs = 14
    
    dic = { 
            "CUON" : 'data/CUON__85_' + str(p)+ '_'  + h + '_monthly_counts.csv',
    
            "SUNY" : 'data/SUNY__85_' + str(p)+ '_'  + h + '_monthly_counts.csv',
        
            "IGRA" : 'data/IGRA2__85_' + str(p)+ '_'  + h + '_monthly_counts.csv',
            "ERA5_1" : 'data/ERA5_1__85_' + str(p)+ '_'  +  h + '_monthly_counts.csv',
            "ERA5_2" : 'data/ERA5_2__85_' + str(p)+ '_'  + h + '_monthly_counts.csv',
        
            "NCAR" : 'data/NCAR__85_' + str(p)+ '_'  + h + '_monthly_counts.csv',
    }
    
    colors = ['black','slateblue', 'orange', 'green', 'grey', 'gold']
    
    plt.subplots(figsize=(10,5) )
    
    for db,c in zip( dic.keys(), colors) :    
        df = pd.read_csv( dic[db] , sep = '\t')
        time, counts = df['date'] , df['counts']
        time = pd.to_datetime(time, infer_datetime_format=True)
        plt.plot(time, counts, label = db, color = c)
    
        
    plt.grid( ls=':' , color = 'lightgray')
    plt.legend(fontsize = fs-1)    
    plt.title("Global Monthly Temperature Records - " + str(p) + " Pa, h:" + h , fontsize=fs , y=1.02 )
    plt.xlim(pd.Timestamp('1940-01-01'),pd.Timestamp('2020-01-01') )
    
    #plt.yscale('log')
    plt.tight_layout()
    plt.savefig('Plots/merging_total_monthly_temperature_' + str(p) +'_' + h + '.png', dpi = 200)
    plt.show()


for p in  [10000,50000,70000,85000]:
    for h in [0,12]:
        
        h, p = str(h), str(p)
        plot(p, h)