import os,sys
import pandas as pd
import numpy as np
import h5py as h5 

from multiprocessing import Pool
from functools  import partial


merged = '/scratch/das/federico/MERGED_25FEB2022/'
files = [ f for f in os.listdir(merged) if '_bef' not in f and '.nc' in f ]


if not os.path.isdir('sensor_tables'):
    os.system('mkdir ' + 'sensor_tables')
    
def make_sensor_table(file):
    print(" *** " , file )
    f = h5.File(file, 'r')
        
        
    a = 0 
    # extracting indices
    indices = f['recordindex'][:]
    
    indices_p = [i for i in indices]
    
    
    # extracting sensor ids
    try:
        
        sensor_id = [ f['observations_table']['sensor_id'][i] for i in indices_p ] 
    except:
        print('Failed Sensor ::: ' , file  )
        return
        
    try:
        
        sensor_id = [ b''.join(f) for f in sensor_id ]
    
        # extracting times 
        date_time = f['recordtimestamp'][:]
        #date_time = [ b''.join(f) for f in date_time ]
        
        #dec = sensor_id[0].view('|S4') 
        #sensor_id = [ f[0] for f in dec ]
        
        df = pd.DataFrame( { 'date_time' : date_time , 'sensor_id' : sensor_id } )
        
        df = df.dropna()
        df = df.drop_duplicates(subset = ['sensor_id'])
        
        date_time_conv = pd.to_datetime( df['date_time'], unit='s',  origin=pd.Timestamp('1900-01-01') )
        df['date_time_conv'] = date_time_conv
        
        df.to_csv( 'sensor_tables/' + file.split('/')[-1].split('_')[0] + '_sensor.csv', sep = '\t' , index= True )
        
    except:
        print('Failed ::: ' , file  )
        return     
    
    return 0
    

POOL = True

if POOL:
    flist  = [ merged + '/' + f for f in files ]
    p = Pool(40)
    func = partial(make_sensor_table)
    out = p.map(func, flist)   
    
else:
    for f in files:
        path = merged + '/' + f
        d = make_sensor_table(path)
    





    
