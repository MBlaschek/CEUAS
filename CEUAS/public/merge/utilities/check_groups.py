""" Check if files contain source configuration, header table and station configuration groups """

import os,sys
import h5py as h5
import pandas as pd
import numpy as np

from multiprocessing import Pool
from functools  import partial



if not os.path.isdir('groups/'):
    os.system('mkdir groups ')
    
    
def check(file):
    
    try:
        f = h5.File(file, 'r')
        stat_c = f["station_configuration"]
        source_c = f["source_configuration"]
        header_t = f["header_table"]
        era5 = f["era5fb"]
        a = open('groups/correct_file_header_station_source_era5fb.txt', 'a+')
        a.write(file + '\n')
        
    except:
        failed = open('groups/failed_file_header_station_source_era5fb.txt', 'a+')
        failed.write(file + '\n')
        
        





merged = '/scratch/das/federico/MERGED_25FEB2022/'


files = [merged  + '/' + f for f in os.listdir(merged) if 'Sensor' not in f and '.nc' in f ]

# run 

REMOVE = False
POOL = True

if POOL:
    poo = Pool(40)
    func = partial(check)
    out = poo.map(func, files)       
else:
    for f in files:
        check(f)
    
    


if REMOVE:
    failed_files = [f.replace('\n','') for f in open( 'groups/failed_file_header_station_source_era5fb.txt', 'r').readlines() ] 
    for g in failed_files:
        os.system('rm ' + g )
