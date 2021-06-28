""" Make a pandas dataframe then a csv file to create the table of constraints for the CDS interface """


import xarray as xr
import pandas as pd
import numpy as np
import os,sys
import h5py as h5

from multiprocessing import Pool
from functools import partial
import datetime

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)

from tqdm import tqdm

import time 




#merged_files = [] 
"""
Eastward wind component	m s-1	Atmospheric eastward wind component at given pressure level
Geopotential height	m	Geopotential altitude from corrected pressure product
Northward wind component	m s-1	Atmospheric northward wind component at given pressure level
Relative humidity	%	Ratio of actual water vapor pressure to saturation water vapor pressure
Specific humidity	Kg/kg	Mass fraction of water vapor in atmospheric volume
Temperature	K	Atmospheric temperature at given pressure level
Wind from direction	Deg (0=N, 90=E)	Atmospheric wind direction at given pressure or height level, following the metorological convention, i.e. wind direction is wind from direction, it increases clockwise and starts with 0 deg if wind comes from North
Wind speed
"""

 

def make_csv_perstation(file):
        
     t = time.time()
                
     station = file.split('_CEUAS_')[0].split('/')[-1]
     
     h5f= h5.File(file,'r' ) 
     ri = h5f['recordindices']
     
     keys = ri.keys()
     
     for v in ['34','36','38','39','85','104','105','106','107','117']:  
             if v not in keys:
                     continue
             
             if os.path.isfile( out_dir + '/' + str(v) + '/' + str(v) + '_' + station + '.csv' ):
                     print('skip ' , v )
                     continue 
             
             print('doing:' , v)
             #index_dic[v] = {}
             mi,ma = h5f['recordindices'][v][0] , h5f['recordindices'][v][-1]
             #index_dic[v]['min'] = mi
             #index_dic[v]['max'] = ma
         
             z_coord = h5f['observations_table']['z_coordinate'][mi:ma]
             obs_val = h5f['observations_table']['observation_value'][mi:ma]
             dt = h5f['observations_table']['date_time'][mi:ma]      
         
         
             """ Filtering """                      
             dfs = []
             for p in   [1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000,92500,100000]:
                     #print(p)
                     res = {}
                     plevel_filter = np.where( (z_coord == p) )[0]
         
                     if not len(plevel_filter) >0:
                             continue
         
                     dt_f = dt[plevel_filter]  
                     z_coord_f = z_coord[plevel_filter] 
                     obs_val_f = obs_val[plevel_filter]
                           
                     nan_filter = np.where (obs_val_f != np.nan )        
                     
                     dt_f = dt_f[nan_filter]
                     z_coord_f = z_coord_f[nan_filter]
                     #indices = np.intersect1d (nan_filter, plevel_filter)
                 
                     res['date_time']=pd.to_datetime(dt_f, unit='s',origin=pd.Timestamp('1900-01-01'))      
                     res['z_coordinate'] = z_coord_f
                     
                     red = pd.DataFrame.from_dict(res)  
                     #red = red.drop_duplicates()
                     ydm_df = pd.DataFrame()
                     
                     if not red.empty:
                             ydm_df['year']  = red['date_time'].dt.year 
                             ydm_df['month'] = red['date_time'].dt.month 
                             ydm_df['day']   = red['date_time'].dt.day   
                             ydm_df['z_coordinate'] = z_coord_f   
                             
                     dfs.append(ydm_df)

             if len(dfs) >0:
                  DF = pd.concat(dfs)
             else:
                     continue 
             
             if not os.path.isdir(  out_dir  ):
                     os.mkdir(out_dir)
             
             if not os.path.isdir(  out_dir + '/' + str(v) + '/' ):
                     os.mkdir(out_dir + '/' + str(v) )
                 
             DF.to_csv( out_dir + '/' + str(v) + '/' + str(v) + '_' + station + '.csv', sep = '\t' , index = False)
                 
             print('TIME delta', time.time() - t )
             
                        
  
def make_csv_perstation(file):
                
     ttt=time.time()
     try:
          station = file.split('_CEUAS_')[0].split('/')[-1]

          print(station)
          processed = open('processed_stations.txt' , 'a')
          empty       = open('empty_stations.txt' , 'a')


          h5f=h5.File(file,'r')
          ri = h5f['recordindices']
          keys = ri.keys()

          for v in ['34','36','38','39','85','104','105','106','107','117']:  
               if v not in keys:
                    continue
               if os.path.isfile( out_dir + '/' + str(v) + '/' + str(v) + '_' + station + '.csv' ):
                    continue 
                          
               tt=time.time()
               mi,ma = h5f['recordindices'][v][0] , h5f['recordindices'][v][-1]
               #index_dic[v]['min'] = mi
               #index_dic[v]['max'] = ma

               #print(time.time()-tt)
               z_coord = h5f['observations_table']['z_coordinate'][mi:ma]
               obs_val = h5f['observations_table']['observation_value'][mi:ma]
               dt = h5f['observations_table']['date_time'][mi:ma]      

               """ Filtering """

               #print(time.time()-tt)
               plevel_filter = np.where( (z_coord == 1000) | (z_coord == 2000) | ( z_coord == 3000) | (z_coord == 5000 ) | (z_coord == 7000) | (z_coord == 10000) | 
                               (z_coord == 15000) |  (z_coord == 20000) |  (z_coord == 25000) |  (z_coord == 30000) |  (z_coord == 40000) |  (z_coord == 50000) | 
                               (z_coord == 70000) |  (z_coord == 85000) |  (z_coord == 92500) |  (z_coord == 100000) ) 
               #print(time.time()-tt)
               nan_filter = np.where (obs_val != np.nan )
               indices = np.intersect1d (nan_filter, plevel_filter)
 
               res = {}


               res['date_time']=pd.to_datetime(dt[indices], unit='s',origin=pd.Timestamp('1900-01-01'))      
               ydm_df=pd.DataFrame()
               if len(indices)>0:
                    ydm_df['year']  = res['date_time'].year 
                    ydm_df['month'] = res['date_time'].month 
                    ydm_df['day']   = res['date_time'].day 
                    ydm_df['z_coordinate']=z_coord[indices]

                    #ydm_df = red[['year', 'month' , 'day', 'z_coordinate']]
                    
                    if not os.path.isdir( out_dir ):
                            os.mkdir(out_dir)
                            
                    if not os.path.isdir( out_dir + '/' + str(v) + '/' ):
                                    os.mkdir(out_dir + '/' + str(v) )

                    ydm_df.to_csv( out_dir + '/' + str(v) + '/' + str(v) + '_' + station + '.csv', sep = '\t' , index = False)

                    processed.write(file + '\n' )
               else:
                    ydm_df.to_csv( out_dir + '/' + str(v) + '/' + str(v) + '_' + station + '_empty.csv', sep = '\t' , index = False)
          
     except MemoryError:
          pass

     return 


              
def combine_station(direc):
     
     files = [ direc+'/' + f for f in os.listdir(direc) ]
     lines = []
     
     for f in tqdm(files):
          lin = open(f,'r').readlines()
          lin=lin[1:]
          lines.extend(lin)

     
     unique = list(set(lines))
     return lines 
          
       
                

""" Main run """

#os.system('rm  -r /raid60/scratch/federico/PROVAAA/' )

out_dir = '/raid60/scratch/federico/constraints_form_JUNE2021/' 

merged_directory = '/raid60/scratch/leo/scratch/converted_v7/' 
merged_files = [ merged_directory + f for f in os.listdir(merged_directory) if '.nc' in f ]

#for m in merged_files:
#      make_csv_perstation(m)
        
#p      = Pool(25)
#func = partial(make_csv_perstation)    
#out   = p.map(func, merged_files)

        
all_lines = []

for v in ['34','36','38','39','85','104','105','106','107','117']:
     print(v)
     all_lines.extend( combine_station(out_dir + '/' + v ) )
     
a = open('constraints_databasev2_JUNE2021.txt', 'w')
a.write('year\tmonth\tday\tpressure\tvariable\n')
for l in all_lines:
     l = l.replace('\n', '\t'+v+'\n')
     a.write(l)
     
a.close()





