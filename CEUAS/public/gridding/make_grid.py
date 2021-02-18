""" Module to create a grid in lat-long boxes.
      1. Read the station_configuration from all the merged files
      2. Create a CSV (\t separated) file with all the station configurations
      3. Create a grid in lat-lon coordinates 
      4. Extract a list of stations included in each box
      
      The output file, called  stations_per_box_size_xxx
      is a numpy fe containing a dictionary that tells which primary ids of the stations are to be used inside
      each box
      
      """

import pandas as pd
import xarray as xr
import os,sys
from tqdm import tqdm
import numpy as np

from multiprocessing import Pool
from functools  import partial

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)

""" Database of the merged files """
merged_database = '/raid60/scratch/federico/DATABASE_JANUARY2021_FIXED_sensor'
file_list = os.listdir(merged_database)

""" Container for the list of station_configurations """
station_configuration = []


def get_sc(file):
    ''' Extract the station_configuration form the merged file '''
    f = merged_database + '/' + file
    sc = xr.open_dataset (f, engine='h5netcdf', group = 'station_configuration', decode_times = True ).to_dataframe()    
    first_sc = sc[0:1]
    first_sc['file']= f 
    print('Done ***' , f )
    return first_sc
    

if not os.path.isfile('all_merged_station_configurations.csv'):
    print('*** Will create the station_configuration file out of the merged files ***')
    p = Pool(30)
    func = partial(get_sc)
    out = p.map(func, file_list)
    sc = pd.concat(out)
    os.system('rm all_merged_station_configurations.csv')
    sc.to_csv('all_merged_station_configurations.csv' , sep = '\t')
    
else:
    print('*** Will read the station_configuration file ***')
    sc = pd.read_csv('all_merged_station_configurations.csv', delimiter = '\t')
    
    """ Fixing wrong convention for values > 180 degrees """
    lat, lon = sc['latitude'].values , sc['longitude'].values
    lat_fix, lon_fix = [ -360 + l if l > 180 else l for l in lat ] ,  [ -360 + lg if lg > 180 else lg for lg in lon ] 
    
    sc['latitude']   = lat_fix
    sc['longitude'] = lon_fix
    
    print(0)
print('*** DONE extracting all the station_configuration tables ')



def make_grid(size=10):
    """ Divide the globe into squared boxes of the given size (default 10 degrees) """
    
    boxes = {}
    for i, n in zip ( range(-180,180, size) , range(1,100000) ):     # longitude
        for j,m in zip( range(-180,180, size) , range (1,100000) ): # latitude 

            center_lat = i + size/2
            center_lon = j + size/2
            
            boxes[str(n)+'_'+str(m)] = {}
            boxes[str(n)+'_'+str(m)]['lat'] = []
            boxes[str(n)+'_'+str(m)]['lon'] = []            
            
            lat_min , lon_min  = center_lat - size/2 ,  center_lon - size/2
            lat_max , lon_max = center_lat + size/2  ,  center_lon + size/2 
            
            boxes[str(n)+'_'+str(m)]['lat'].append( float( "{:.2f}".format(lat_min) ) ) 
            boxes[str(n)+'_'+str(m)]['lat'].append( float( "{:.2f}".format(lat_max) ) ) 
            boxes[str(n)+'_'+str(m)]['lon'].append(float( "{:.2f}".format (lon_min) ) )
            boxes[str(n)+'_'+str(m)]['lon'].append(float( "{:.2f}".format (lon_max) ) )
    
    return boxes


""" Creating the gridded boxes """
boxes = make_grid(size=10)

def get_box(sc= '', boxes = ''):
    """ For each file in the station_configuration table, it associates the correct gridded box
          by comparing latitude and longitude """
    
    size = str(int(abs(boxes['1_1']['lat'][1] - boxes['1_1']['lat'][0]) ) ) # extracting the size of the box 
    found_boxes = {}
    
    for c,r in sc.iterrows():
        print('Done row::: ' , c)
        lat = r.latitude
        lon = r.longitude
        file = r.file
        primary = r.primary_id
        
        for k in boxes.keys():
            #print(' Checking box ::::' , k )
            if k not in found_boxes.keys():     
                found_boxes[k] = {}
                found_boxes[k]['lat'] = []
                found_boxes[k]['lon'] = []
                found_boxes[k]['files'] = []
                
                found_boxes[k]['lat'].append(boxes[k]['lat'][0] ) , found_boxes[k]['lat'].append(boxes[k]['lat'][1] )
                found_boxes[k]['lon'].append(boxes[k]['lon'][0] ), found_boxes[k]['lon'].append(boxes[k]['lon'][1] )
                
            try:
                if  lat < boxes[k]['lat'][0] or lat > boxes[k]['lat'][1]:
                    continue
                elif  lon < boxes[k]['lon'][0] or lon > boxes[k]['lon'][1]:
                    continue
                else:
                    found_boxes[k]['files'].append( primary )  # storing the primary id of the stations inisde each box 
            except:
                continue
            
    np.save('stations_per_box_size_' + size  , found_boxes )
    return found_boxes





if os.path.isfile('stations_per_box.npy'):
    print('*** Loading the gridded stations :::')
    assigned_boxes = np.load('stations_per_box.npy' , allow_pickle= True ).item()
else:
    assigned_boxes = get_box(sc= sc, boxes = boxes)
    
print('*** DONE extracting stations per box *** ')
   
   