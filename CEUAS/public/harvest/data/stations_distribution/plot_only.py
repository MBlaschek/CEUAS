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



### only the keys are actually used
data_directories = {  'ncar'      : '/raid60/scratch/federico/netCDF_converted_Jan2020/ncar'    ,
                                  'igra2'     : '/raid60/scratch/federico/netCDF_converted_Jan2020/igra2'               ,
                                  'era5_1'    :  '/raid60/scratch/federico/netCDF_converted_Jan2020/ready_for_merging/era5_1'          ,
                                  'era5_1759' : '/raid60/scratch/federico/netCDF_converted_Jan2020/ready_for_merging/era5_1759'  ,
                                  'era5_1761' : '/raid60/scratch/federico/netCDF_converted_Jan2020/ready_for_merging/era5_1761'  ,
                                  'era5_3188' : '/raid60/scratch/federico/netCDF_converted_Jan2020/ready_for_merging/era5_3188'  ,     
                      
                                  'bufr'      : '/raid60/scratch/federico/netCDF_converted_Jan2020//bufr'                     }






def plot(data_directories, direc = 'numpys/'):
  bins = 20
  fig,ax = plt.subplots()
  
  #plt.ylim(-0.5, 15.5)
  #plt.xlim(-0.5, 15.5)
  #plt.xlabel('Date')
  plt.ylabel('Records counts')  
  Min, Max = [], [] 
  
  # bin = 4 years
  BINS = { 'era5_1'    : 10 , 
           'era5_1759' : 11 ,
                  
           'era5_1761' : 10 ,
           'era5_3188' : 6  ,
           'igra2'     : 25 ,
           'ncar'      : 22 ,
           'bufr'      : 5  ,       }
  
  for d,p in data_directories.items() :  
         print('loading' , d , ' ' , p )
         #filenp = np.load( 'dic_records_TEST_' + d + '.npy' , allow_pickle = True).item() 
         filenp = np.load( direc + '/dic_records_' + d + '.npy' , allow_pickle = True).item()          
         data = filenp[d]
         #print( d , ' ' , min(data) , ' ' , max(data) , ' ' , len(data) )
         #Min.append(min(data))
         #Max.append(max(data))
         if d == 'era5_3188':
              plt.hist(data,  histtype='stepfilled', bins = BINS[d], stacked = False, label = '3188' , color = 'black', alpha = 0.7 , density = False , edgecolor='black')
         elif d == 'era5_1761':
              plt.hist(data,  histtype='stepfilled', bins = BINS[d], stacked = False, label = '1761' , color = 'yellow', alpha = 0.6 , density = False ,  edgecolor='black')     
         elif d == 'era5_1759':
              plt.hist(data,  histtype='stepfilled', bins = BINS[d] , stacked = False, label = '1759' , alpha = 0.7 , density = False , edgecolor='black')
         else:
              plt.hist(data,  histtype='stepfilled', bins = BINS[d] , stacked = False, label = d , alpha = 0.7 , density = False , edgecolor='black')
 
              
  #ax.yaxis.set_major_locator(plt.MaxNLocator(5))
  Min =  datetime.strptime('1915-01', '%Y-%m')  
  Max =  datetime.strptime('2020-01', '%Y-%m')  

  plt.xlim(Min, Max )
  
  #ax.xaxis.set_major_locator(plt.MaxNLocator(20))
  #plt.locator_params(nbins=10)
  plt.xticks( rotation=45)
  plt.legend(loc = 'upper left', fontsize =10)
  os.system('mkdir plots_records_counts')
  plt.savefig('plots_records_counts/MAY.png',  bbox_inches='tight')
  plt.close()  
  print('plotted')



a = plot(data_directories)
                
        
        
        
        
    
