import pandas as pd
import numpy as np
import glob 
import os,sys

old_stat_conf = '/users/staff/federico/GitHub/CEUAS_master_SEPTEMBER2021/CEUAS/CEUAS/meta/inventory_comparison/code/station_configuration_OLD/'
new_stat_conf = '/users/staff/federico/GitHub/CEUAS_master_SEPTEMBER2021/CEUAS/CEUAS/meta/inventory_comparison/code/station_configuration/'

databases = [ 'era5_2', 'era5_1759', 'era5_1761', 'era5_3188', 'bufr', 'igra2', 'era5_1', 'ncar']

harvested = '/scratch/das/federico/COP2_HARVEST_FEB2022/'
merged = '/scratch/das/federico/MERGED_25FEB2022/'
for d in databases:
    
    new = pd.read_csv( new_stat_conf + '/' + d + '_station_configuration_extended.csv',  sep='\t' )
    old = pd.read_csv( old_stat_conf + '/' + d + '_station_configuration_extended.csv',    sep='\t' )
    
    p_old   = list(np.unique(old['primary_id'][:] ) ) 
    p_new = list(np.unique(new['primary_id'][:] ) )
    
    old_to_remove = [g for g in p_old if g not in p_new ]
    
    #print(old_to_remove)
    
    print(d , "  " , len(p_old) , " " , len(p_new) , " remove: " , len(old_to_remove) )
    
    for f in old_to_remove:
        files = glob.glob(harvested + '/' + d + '/' + f + "*" )
        print(files)
        for g in files:
            if os.path.isfile(g):
                a = 0
                os.system('rm ' + g )       
        
        files = glob.glob(merged +  '/' + f + "*"   )
        for g in files:
            if os.path.isfile(g):
                os.system('rm ' + g )         
                print(g)
                
    a = 0