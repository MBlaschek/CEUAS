# Remove files from harvested directory if the station id is not in the station_configuration 


import os,sys
import pandas as pd
import numpy as np
import glob

sc_d = '/users/staff/federico/GitHub/CEUAS_master_SEPTEMBER2021/CEUAS/CEUAS/public/harvest/data/station_configurations/'

f = '_station_configuration_extended.csv'

databases = [ 'era5_2', 'era5_1759', 'era5_1761', 'era5_3188', 'bufr', 'igra2', 'era5_1', 'ncar']

harvested = '/scratch/das/federico/COP2_HARVEST_FEB2022/'


wrong_stat_ids = {}


# cleaning harvested directory 
for d in databases:
    
    wrong_stat_ids[d] = []
    
    sc = sc_d + '/' + d + f 
    df = pd.read_csv(sc, sep = '\t')
    
    primary = df['primary_id'].values
    primary = list(np.unique(primary) )
    
    h_dir = harvested + '/' + d 
    
    files = [ h_dir + '/' + f for f in os.listdir(h_dir)  if '-1_' not in f  and 'unid' not in f ]
    
    for fi in files:
        stat = fi.split('/')[-1].split('_')[0] 
        if stat not in primary:
            wrong_stat_ids[d].append(fi)
            
    """
    stat_ids = [ f.split('_')[0] for f in os.listdir(h_dir) if '-1_' not in f  and 'unid' not in f ]
    stat_ids = list ( np.unique( stat_ids ) ) 
    wrong = [s for s  in stat_ids if s not in primary ]
    wrong_stat_ids[d] = wrong 
    """
    
    print(d , ' ' , len(wrong_stat_ids[d]) )
    
    a = 0
    
a = 0

for k in wrong_stat_ids.keys():
    files = wrong_stat_ids[k]
    for f in files:
        print(f)
        os.system('rm ' + f )
        

# cleaning merged directory 
sc = '/users/staff/federico/GitHub/CEUAS_master_SEPTEMBER2021/CEUAS/CEUAS/public/merge/CUON_station_configuration.csv'
df = pd.read_csv(sc, sep = '\t')

primary = df['primary_id'].values
primary = list(np.unique(primary) )

m_dir = '/scratch/das/federico/MERGED_25FEB2022/'
files = [ m_dir + '/' + f for f in os.listdir(m_dir)  if '.nc' in f   ]
wrong_merged = []
for fi in files:
    stat = fi.split('/')[-1].split('_')[0] 
    if stat not in primary:
        wrong_merged.append(fi)
        
a = 0