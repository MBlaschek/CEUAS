import os,sys
import h5py as h5
import pandas as pd
import numpy as np

from multiprocessing import Pool
from functools  import partial

merged = '/scratch/das/federico/MERGED_25FEB2022/'
harvested = '/scratch/das/federico/COP2_HARVEST_FEB2022/'

if not os.path.isdir('coordinates'):
    os.system('mkdir coordinates')
    
def check_coord(limit, db, file, write = False):
        
    ot = h5.File(file, 'r')['observations_table']
    lat = ot['latitude'][:]
    lon = ot['longitude'][:]
    
    cons = True
    
    min_lat, max_lat = min(lat), max(lat)
    min_lon, max_lon = min(lon), max(lon)
    
    
    if ( abs( min_lat - max_lat ) > limit ):
        cons = False
        
    if ( abs( min_lon - max_lon ) > limit ):
        cons = False    
        
    if cons:
        print('Cons ' + file )
        print(file ,  ' LAT min-MAX :  ' , min_lat , '  ', max_lat )
        print(file,   'LON min-MAX :  ' , min_lon , '  ', max_lon  )
        
        if write:
            a = open('coordinates/' + db + '_consistent_coord.txt', 'a+')
            a.write(file + '\n')
    
    else:
        print('Inconsistent ****  ' + file )
        print( '\t\t LAT min-MAX :  ' ,  min_lat , '  ', max_lat  )
        print( '\t\t LON min-MAX :  ' ,  min_lon , '  ', max_lon    )
        
        if write:
            b = open('coordinates/' + db + '_inconsistent_coord.txt', 'a+')
            b.write(file + '\n')        
        
    
    
databases = [ 'era5_2', 'era5_1759', 'era5_1761', 'era5_3188', 'bufr', 'igra2', 'era5_1', 'ncar']

databases = [ 'merged' ]

databases = [ 'era5_2', 'era5_1759', 'era5_1761', 'era5_3188', 'bufr', 'igra2', 'era5_1', 'ncar' , 'merged' ]


databases = [ 'merged' ]

#databases = [ 'era5_2', 'era5_1759', 'era5_1761', 'era5_3188', 'bufr', 'igra2', 'era5_1', 'ncar']



databases = [ 'era5_2', 'era5_1759', 'era5_1761', 'era5_3188', 'bufr', 'igra2', 'era5_1', 'ncar' , 'merged' ]


# options 
POOL = True
CHECK_INCONSISTENT = False
DISCARD_CONSISTENT = True
WRITE = True


databases = [ 'era5_2', 'era5_1759', 'era5_1761', 'era5_3188', 'era5_1', 'ncar'  ]

for d in databases:
    
    print( 'Analysing database::: ' , d )
    try:
        inconsisten_stations = [ f.split('/')[-1].split('_')[0] for f in open(  'coordinates/merged_inconsistent_coord.txt', 'r').readlines() ]
    except:
        inconsisten_stations = []
    
    
    if d == 'merged':
        p = '/scratch/das/federico/MERGED_25FEB2022/' 
        if DISCARD_CONSISTENT:
            cons = [ f.replace('\n','') for f in open(  'coordinates/merged_consistent_coord.txt', 'r').readlines() ]            
            files = [ p + f for f in os.listdir(p ) if '.nc' in f and 'Sen' not in f and p+f not in cons ]
        else:
            files = [ p + f for f in os.listdir(p ) if '.nc' in f and 'Sen' not in f ]
            
    else:
        p = '/scratch/das/federico/COP2_HARVEST_FEB2022/' + d + '/'                
        files = [ p + '/' + f for f in os.listdir(p) if '.nc' in f  and '-1' not in f ]

    if CHECK_INCONSISTENT:
        
        ff = []
        for f in files:
            s = f.split('/')[-1].split('_')[0]
            if s in inconsisten_stations:
                ff.append(f)

        files = ff 

    if POOL:
        poo = Pool(40)
        
        func = partial(check_coord, 1, d, write = WRITE)
        out = poo.map(func, files)   
    else:
        for file  in files:
            dummy = check_coord(1, d, file, write = WRITE)
            c = 0
            
    a = 0
    

a = 0


#file = '/scratch/das/federico/MERGED_25FEB2022/0-20000-0-47104_CEUAS_merged_v1.nc'
#dummy = check_coord(1, 'merged', file, write = False)

"""
databases = [ 'era5_2', 'era5_1759', 'era5_1761', 'era5_3188', 'bufr', 'igra2', 'era5_1', 'ncar'  ]

# global check of consistentcy starting from one single wrongly merged file 
try:
    inconsisten_stations = [ f.split('/')[-1].split('_')[0] for f in open(  'coordinates/merged_inconsistent_coord.txt', 'r').readlines() ]
except:
    inconsisten_stations = []
        

for s in inconsisten_stations:
    print('**************** STATIONFAULTY' , s )
    for db in databases:
        
        files_check = []
        p = '/scratch/das/federico/COP2_HARVEST_FEB2022/' + db + '/'                
        files = [ p + '/' + f for f in os.listdir(p) if '.nc' in f  and '-1_' not in f ]
        
        # todo

        for f in files:
            stat = f.split('/')[-1].split('_')[0]
            if s in stat:
                check_coord(1, db, f, write = False)
        
    a = 0

"""
