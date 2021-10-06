import os,sys
import h5py as h5
import pandas as pd
import netCDF4 as nc

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker

# containing utilities for data reading and counting
from merging_dataset_comparison_UTILS import *

from multiprocessing import Pool
from functools  import partial


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


# /users/staff/a1400070/CEUAS/CEUAS/public/testing/igra_h_single_stats_new

# Some wrapping functions for multiprocesses and combination of output 

def wrapper( ds, var, plevel, hour, file):
    
    try:
        
        if 'SUNY' in ds:
            data_dic = get_SUNY_data(var, hour, plevel, file)
        elif "HARM" in ds:
            data_dic = get_HARM_data(plevel, hour, file)
        else:
            data_dic = get_CUON_data(var, plevel, hour , False, file)        
            
        if data_dic:
            counts = get_monthly_counts_dic(data_dic)    
            return counts
        else:
            print('No valid counts in the file ', file )
        
    except:
        print("Fail !!! " , file )
        



def sum_counts(dic_df, dataset = "", var="", plevel="", hour=""):
    
    all_counts = {} # holder for the results 
    
    for df in dic_df:
        try:  
            if not df.empty:
                all_dt   = df['time'][:]
                counts = df['counts'][:]
                # loop over the dt of the df of each station
                # if dt not a key of the dictionary, creates it
                for t,count in zip(all_dt, counts ):
                    tt =  t.strftime("%Y%m")
                    if tt in all_counts.keys():                        
                        all_counts[tt].append(count)
                    else:
                        all_counts[tt] =  []
                        all_counts[tt].append(count)
            else:
                print('Empty --- ')
                continue
        except:
            print('Exception ---')
            
    dates,counts = [],[]
        
    for dt in all_counts.keys():
        dates.append(dt)
        counts.append( sum(all_counts[dt] ) )
        
    df = pd.DataFrame( { 'date':dates, 'counts':counts} )
        
    df = df.sort_values(by = ['date'])
    df.to_csv( 'data/' + dataset +'_' + '_' + str(var) + '_' + str(plevel) + '_' + str(hour) + '_monthly_counts.csv' , sep = '\t' , index=False)
    return df 
    




#harm_dir = '/users/staff/a1400070/CEUAS/CEUAS/public/testing/igra_h_single_stats_new'
#files = os.listdir( harm_dir )
#f = harm_dir + '/' + files[3]
#dummy = get_HARM_data(f, plev = 50000, h=0)
    
# Location of datasets 
dbs = { 'SUNY': "/users/staff/leo/fastscratch/SUNY/homo-raw-subdaily-station/" ,
        
              'IGRA2': "/scratch/das/federico/MAY2021_HARVEST_secondary/igra2/" ,
              
              'NCAR': "/scratch/das/federico/MAY2021_HARVEST_secondary/ncar/" ,
              'CUON': "/mnt/users/scratch/leo/scratch/converted_v7/" ,
              'ERA5_1': "/scratch/das/federico/MAY2021_HARVEST_secondary/era5_1/" ,
              'ERA5_2': "/scratch/das/federico/MAY2021_HARVEST_secondary/era5_2/" ,
              
              'HARM': "/users/staff/a1400070/CEUAS/CEUAS/public/testing/igra_h_single_stats_new",
              
              'IGRA2_merg':          "/scratch/das/federico/MERGED_ONLY_IGRA_SEPT2021/" ,
              'IGRA2_merg_suny': "/scratch/das/federico/MERGED_ONLY_IGRA_SEPT2021/",
              'missing_suny':        "/scratch/das/federico/MERGED_ONLY_IGRA_SEPT2021/",

              'NCAR_merg':          "/scratch/das/federico/MERGED_ONLY_NCAR_SEPT2021/" ,

              
        }





# WHAT = plot or run 
WHAT = 'plot'
POOL = True
#POOL = False


if WHAT == 'run':
    pool = Pool(30)

    # Loop over datasets 
    for db in ['IGRA2_merg', 'CUON' ] :  #  'SUNY', 'IGRA2', 'IGRA2_merg' , 'IGRA2_merg_suny' , 'CUON', 'ERA5_1', 'ERA5_2' , 'HARM', "missing_suny" 
        
        if db == 'NCAR':
            files  = [dbs[db] + '/' + f  for f in os.listdir(dbs['NCAR'])  if 'wind' not in f ]

        elif db == 'IGRA2_merg_suny':
            
            # I pick only the files used by the SUNY analysis
            
            matching_CUON_SUNY_stations = [] # store the psuedo WIGOs if for the SUNY station in the CUON db
            
            SUNY_stations = [ s.split('.nc')[0] for s in os.listdir(dbs['SUNY']) if 'WMO' not in s  ]
            igra_harvested = [ f for f in os.listdir(dbs['IGRA2']  ) if '.nc' in f ]

            print('+++++ Extracting the common WIGOs for SUNY and CUON stations from the merged IGRA2 only merged dataset')
            for s in SUNY_stations:
                for i in igra_harvested:
                    if s in i:
                        stat = i.split('_igra2_')[0]
                        matching_CUON_SUNY_stations.append(stat)
                        
                        #print(0)
            print('+++++ Retrieving exact file names for SUNY and CUON stations from the merged IGRA2 only merged dataset')
            
            files = []
            for m in matching_CUON_SUNY_stations:  # there are 1161 matching SUNY-CUON stations (tot SUNY = 1185 , 14 WMO stations excluded )
                for s in [ f for f in os.listdir(dbs['IGRA2_merg']  ) if '.nc' in f ]:
                    if m in s:
                        #print(s)
                        files.append(dbs['IGRA2_merg'] + '/'  + s)
            files = list (np.unique(files) )  # removes double counting 
            
        elif db == 'missing_suny':
            FF = pd.read_csv('data/files_suny.csv')['files']
            files = [dbs[db] + '/' + f for f in os.listdir( dbs[db] )  if dbs[db] + '/' + f not in FF.values ]
            # df = pd.DataFrame (  {'files': files}).to_csv('files_not_in_suny.csv', index=False) 
            #print(0)
            
            
        else:
            files =  [dbs[db] + '/' + f  for f in os.listdir( dbs[db] ) if '.nc' in f  ]
            
        h = 0
        
        #files = files[:100] # TODO to make it smaller for tests 
        
        #files =  [ f for f in files if '10393' in f ]


        #files = list(pd.read_csv('data/files_suny.csv')['files'])
        #suny_miss = list(pd.read_csv('data/files_not_in_suny.csv')['files'])
        #files.extend(suny_miss)
        
        files = list (np.unique(files) ) 

        for v in [85]:
            if db == "SUNY" and v != 85:
                print("Skipping SUNY since no data available except temperature")
                continue
    
            for p in [10000,50000,70000,85000]:  # [10000,50000,70000,85000]
                print("+++++ Doing " , db , ' ' , h , ' ' , v , ' ' , p )
                
                if POOL:
                    func = partial(wrapper, db, v, p, h )
                    out = pool.map(func, files)
                    out = [o for o in out if isinstance(o, pd.DataFrame) ]                
                    print("+++++ Done " , db , ' ' , h , ' ' , v , ' ' , p )
                else:
                    out = []
                    for file in files:
                        print('File: ', file )
                        counts_df = wrapper( db, v, p, h, file)
                        out.append(counts_df)
                        
                print("Aggregating results")
                c = sum_counts(out, dataset = db, var=v, plevel=p, hour=h)
     
else:

    def plot(p, h, label=''):
        
        fs = 14  
        
        dic = { 
                        "CUON" : 'data/CUON__85_' + str(p)+ '_'  + h + '_monthly_counts.csv',
                        "SUNY"  : 'data/SUNY__85_' + str(p)+ '_'  + h + '_monthly_counts.csv',
                        "IGRA"                    : 'data/IGRA2__85_' + str(p)+ '_'  + h + '_monthly_counts.csv',
                        "IGRA_merg"          : 'data/IGRA2_merg__85_' + str(p)+ '_'  + h + '_monthly_counts.csv',
                        "NCAR_merg"          : 'data/NCAR_merg__85_' + str(p)+ '_'  + h + '_monthly_counts.csv',
                        
                        "IGRA_merg_suny" : 'data/IGRA2_merg_suny__85_' + str(p)+ '_'  + h + '_monthly_counts.csv',
                        "missing_suny"       : 'data/missing_suny__85_' + str(p)+ '_'  + h + '_monthly_counts.csv',
                        "HARM"  : 'data/HARM__85_' + str(p)+ '_'  + h + '_monthly_counts.csv',
                        
                        #"ERA5_1" : 'data/ERA5_1__85_' + str(p)+ '_'  +  h + '_monthly_counts.csv',
                        #"ERA5_2" : 'data/ERA5_2__85_' + str(p)+ '_'  + h + '_monthly_counts.csv',
                    
                        #"NCAR" : 'data/NCAR__85_' + str(p)+ '_'  + h + '_monthly_counts.csv',
                }

        colors = ['black','slateblue', 'orange', 'green', 'grey', 'gold', 'magenta']
        
        plt.subplots(figsize=(10,5) )
        
        ### runnign mean to smooth the curve
        
        #for k in ['obs_adj' , 'obs', 'fg_dep_adj' , "fg_dep" ]:
        #    mean = pd.Series( df[k]  )
        #    df[k+"_mean"] = mean.rolling(30, center=True).mean()   
            
            
        for db,c in zip( dic.keys(), colors) :    
            df = pd.read_csv( dic[db] , sep = '\t')
            df["mean"] = df['counts'].rolling(100, center=True).mean() 
            
            time, counts = df['date'] , df['mean']
            time = pd.to_datetime(time, format= "%Y%m")
            plt.plot(time, counts, label = db, color = c)
        
            
        plt.grid( ls=':' , color = 'lightgray')
        plt.legend(fontsize = fs-1)    
        plt.title("Global Monthly Temperature Records - " + str(p) + " Pa, h:" + h , fontsize=fs , y=1.02 )
        plt.xlim(pd.Timestamp('1940-01-01'),pd.Timestamp('2020-01-01') )
        plt.ylabel("Record counts / month "  , fontsize=fs )
        
        #plt.yscale('log')
        plt.tight_layout()
                
        plt.savefig('Plots/merging_total_monthly_temperature_' + str(p) +'_' + h + '_' + label + '.png', dpi = 200)
        plt.show()
        
    def plot_pretty(p, h, test_igra=True, label = ''):
        print("*** I am pretty-plotting ::: " , p , h  )
        fs = 14   
        
        dic = { 
                        "CUON (this work)" : 'data/CUON__85_' + str(p)+ '_'  + h + '_monthly_counts.csv',
                        "IGRA2"        : 'data/IGRA2_merg__85_' + str(p)+ '_'  + h + '_monthly_counts.csv',    
                        "NCAR"          : 'data/NCAR_merg__85_' + str(p)+ '_'  + h + '_monthly_counts.csv',
                        
                        "Zhou et al. 2021"  : 'data/SUNY__85_' + str(p)+ '_'  + h + '_monthly_counts.csv',
                        "Madonna et al. 2021"  : 'data/HARM__85_' + str(p)+ '_'  + h + '_monthly_counts.csv',
                        
                        #"ERA5_1" : 'data/ERA5_1__85_' + str(p)+ '_'  +  h + '_monthly_counts.csv',
                        #"ERA5_2" : 'data/ERA5_2__85_' + str(p)+ '_'  + h + '_monthly_counts.csv',
                    
                        #"NCAR" : 'data/NCAR__85_' + str(p)+ '_'  + h + '_monthly_counts.csv',
                }

        colors = ['black','slateblue', 'lightgreen', 'red', 'grey', 'gold', 'magenta']
        colors = ['royalblue','gold', 'lightgreen', 'salmon', 'lightblue']
        colors = ['royalblue','salmon', 'lightgreen', 'gold', 'lightblue']
        
        plt.subplots(figsize=(10,5) )
        
        ### runnign mean to smooth the curve
        
        #for k in ['obs_adj' , 'obs', 'fg_dep_adj' , "fg_dep" ]:
        #    mean = pd.Series( df[k]  )
        #    df[k+"_mean"] = mean.rolling(30, center=True).mean()   
            
            
        #
        # Smoothing procedure: using a rolling mean and then appling a cubic interpolation
        # 
        
        for db,c in zip( dic.keys(), colors) :    
            df = pd.read_csv( dic[db] , sep = '\t')
            time = pd.to_datetime(df['date'], format= "%Y%m")
            
            df["counts"] = df['counts'].rolling(50, center=False).mean() 
            
            df_r = df.set_index(time)
            
            # smoothing procedure, interpolating 
            # first setting the index to a datetime type
            # then interpolating with a certain frequency e.g. freq='15d' 
            
            index_hourly = pd.date_range(pd.Timestamp('1940-01-01'), pd.Timestamp('2020-12-31'), freq='20d')
            index_hourly = pd.to_datetime(index_hourly)
            df_smooth =  df_r.reindex(index=index_hourly).interpolate('cubic')
            
            # other option: using ONLY a rolling mean 
            #df["mean"] = df['counts'].rolling(30, center=False).mean() 
            #time, counts = df['date'] , df['mean']
            #time = pd.to_datetime(time, format= "%Y%m")
            
            time, counts = df_smooth.index , df_smooth['counts']
            plt.plot(time, counts, label = db, color = c , lw = 2)
        
            
        plt.grid( ls=':' , color = 'lightgray')
        plt.legend(fontsize = fs-2, loc = 'lower right')    
        
        if h=='0':
            h = '00'
            
        plt.title("Global Monthly Temperature Records - " + str(p) + " Pa, h:" + h , fontsize=fs , y=1.02 )
        plt.ylabel("Record counts / month "  , fontsize=fs )
        
        plt.xlim(pd.Timestamp('1940-01-01'),pd.Timestamp('2020-12-31') )
        
        #plt.yscale('log')
        plt.tight_layout()
        

        plt.savefig('Plots/merging_total_monthly_temperature_' + str(p) +'_' + h + '_' + label + '.png', dpi = 200)
        plt.show()        
        
        
    for p in  [70000,50000]: #10000,50000,70000,85000
        for h in [0,12]:    # [0,12]
            h, p = str(h), str(p)
            #plot(p, h)
            plot_pretty(p, h,  label = 'pretty')            
    
