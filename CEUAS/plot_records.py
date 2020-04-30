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

data_directories = {  'ncar'      : '/raid60/scratch/federico/Harvested_Database/ncar'                         ,
                      'igra2'     : '/raid60/scratch/federico/Harvested_Database/igra2'                        ,
                      #'era5_1'    : '/raid60/scratch/federico/Harvested_Database/era5_1'                             ,
                      'era5_2'    :'/raid60/scratch/leo/scratch/era5/odbs/2/era5_2' ,
                      'era5_1'    :'/raid60/scratch/leo/scratch/era5/odbs/1/era5_1' ,

                      'era5_1759' : '/raid60/scratch/federico/netCDF_converted_Jan2020/ready_for_merging/era5_1759'  ,
                      'era5_1761' : '/raid60/scratch/federico/netCDF_converted_Jan2020/ready_for_merging/era5_1761'  ,
                      'era5_3188' : '/raid60/scratch/federico/netCDF_converted_Jan2020/ready_for_merging/era5_3188'  ,     
                      
                      'bufr'      : '/raid60/scratch/federico/Harvested_Databases/bufr'                     }


def plot(data_directories, dic_all):
  bins = 20
  fig,ax = plt.subplots()
  
  #plt.ylim(-0.5, 15.5)
  #plt.xlim(-0.5, 15.5)
  #plt.xlabel('Date')
  plt.ylabel('Records counts')  
  Min, Max = [], [] 
  for d in dic_all:     
         
         filenp = np.load( 'dic_records_TEST_' + d + '.npy' , allow_pickle = True).item() 
         data = filenp[d]
         
         Min.append(min(data))
         Max.append(max(data))
         plt.hist(data,  histtype='stepfilled',  stacked = False, label = d , alpha = 0.4 , density = False)
         
  ax.xaxis.set_major_locator(plt.MaxNLocator(8))
  ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    
  plt.xlim(min(Min), max(Max) )
    
  plt.locator_params(nbins=10)
  plt.xticks( rotation=45)
  plt.legend(loc = 'upper left', fontsize =10)
  os.system('mkdir plots_records_counts')
  plt.savefig('plots_records_counts/record_counts_TRY.png',  bbox_inches='tight')
  plt.close()  

  
  

for d,p in zip(  [db], [data_directories[db] ]  ):
    print(' I am  here')
    dic_all = {}    
    dic = ''
    
    print('Processing the dataset: ' , d , ' ' , p   ) 
    files  = [ k for k in os.listdir(p) if '.nc' in k and 'wind' in k ]  ### edit here
    
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
    np.save('dic_records_newncar_' + d  , dic_all,  allow_pickle=True )      
          
    """
    File                        = nc.Dataset(  p + '/' + f  )   
    date_time              = File.variables['recordtimestamp'][:].data
    time_offset            = File.groups['observations_table']['date_time'].units
    time_offset_value  = time_offset.split('since ') [1]      
    time_offset_value  = datetime.strptime(time_offset_value, '%Y-%m-%d %H:%M:%S')

    if 'minutes' in  time_offset:
          date_time_delta = [ timedelta(minutes = float(i) )  + time_offset_value for i in date_time ]
    elif 'hours' in time_offset:
          date_time_delta = [ timedelta(hours = float(i) )      + time_offset_value  for i in date_time ]    
    elif 'seconds' in time_offset:                        
          date_time_delta = [ timedelta(seconds = float(i) )  + time_offset_value  for i in date_time ]   
          
         
    #rts  = [ i.replace(minute=0, second=0) for i in date_time_delta ]           
    #rts = [  i.strftime("%Y%m")   for i in date_time_delta ]
    #x = [np.datetime64(i) for i in ['2010-06-01T00:00:00.000000000', '2010-12-02T00:00:00.000000000']]
    #y = x.strftime('%Y%m') 
    #[str(i.date()) for i in x]
    
    np.save('dic_all_plots' , dic_all,  allow_pickle=True )      

    """
        
    #a = plot(d, dic[d])
        
 
#np.save('dic_all_plots' , dic_all,  allow_pickle=True )      

#a = plot(dic_all)
print('hello')    
#a = plot(d, dic[d])
                
        
        
        
        
    
