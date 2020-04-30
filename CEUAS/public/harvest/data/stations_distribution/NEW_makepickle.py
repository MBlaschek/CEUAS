import os, sys
import xarray as xr
import netCDF4 as nc
import datetime
from datetime import datetime, timedelta
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import pandas as pd
# https://stackoverflow.com/questions/49204453/how-can-i-get-year-month-day-from-a-numpy-datetime64?noredirect=1&lq=1


import argparse
""" Parsing the string containing the list of stations to be merged """ 
parser = argparse.ArgumentParser(description="Make CDM compliant netCDFs")
parser.add_argument('--database' , '-d', 
                        help="Database"  ,
                  type = str)

args = parser.parse_args()
db = args.database

data_directories = {  'ncar'      : '/raid60/scratch/federico/Harvested_Databases/ncar'                         ,
                      'igra2'     : '/raid60/scratch/federico/Harvested_Databases/igra2'                        ,
                      #'era5_1'    : '/raid60/scratch/federico/Harvested_Database/era5_1'                             ,
                      'era5_2'    :'/raid60/scratch/leo/scratch/era5/odbs/2/era5_2' ,
                      'era5_1'    :'/raid60/scratch/leo/scratch/era5/odbs/1/era5_1' ,

                      'era5_1759' : '/raid60/scratch/federico/netCDF_converted_Jan2020/ready_for_merging/era5_1759'  ,
                      'era5_1761' : '/raid60/scratch/federico/netCDF_converted_Jan2020/ready_for_merging/era5_1761'  ,
                      'era5_3188' : '/raid60/scratch/federico/netCDF_converted_Jan2020/ready_for_merging/era5_3188'  ,     
                      
                      'bufr'      : '/raid60/scratch/federico/Harvested_Databases/bufr'                     }

  
for d,p in zip(  [db], [data_directories[db] ]  ):
    print(' I am  here')
    dic_all = {}    
    dic = ''
    
    print('Processing the dataset: ' , d , ' ' , p   ) 
    files  = [ k for k in os.listdir(p) if '.nc' in k ]  ### edit here
    
    for f in files:
        print(files.index(f) )
        
        try:             
          File = xr.open_dataset( p + '/' + f  , engine = 'h5netcdf'   )   
          rts = File['recordtimestamp'].data
        except:
          continue
        
        x = pd.to_datetime(rts)

        try:         
          dic =  dic.append(x)
        except:
          dic = x  
    dic_all[d]  = dic       
    np.save('NEW_NUMPYS/records_' + d  , dic_all,  allow_pickle=True )      
    print('Finished ***', d )
 
    
print('Finished ALL **** ')                
        
        
        
        
    
