"""
Contain utilities for running the new version of the harvester with year_split 
"""

import os,sys
import pandas as pd
import numpy as np



def filelist_cleaner(lista, d=''):
    """ Removes unwanted files that might be present in the database directories """

    print('Cleaning the list of files to be converted')
    if d == 'ncar':
        cleaned = [ l for l in lista if '.nc' not in l and '.py' not in l ]
        
    if d == 'bufr':
        cleaned = [ l for l in lista if '.bfr' in l and 'era5.' in l and '+100-' not in l and 'undef' not in l and '99999' not in l]
        
    if d in ['era5_1759', 'era5_1761']:
        cleaned = [ l for l in lista if '.nc' not in l and '.conv.' in l and 'py' not in l and '.gz' in l and '.g.' not in l ] 
        
    if d =='igra2':
        cleaned = [ l for l in lista if '-data' in l and '.zip' not in l ]
        
    if d == 'era5_3188': #  era5.3188.conv._C:4687.gz
        cleaned = [ l for l in lista if '3188.conv' in l and '.nc' not in l and '.gz' in l]
        
    if d == 'era5_1':
        cleaned = [ f for f in lista if '.conv._' in f and '.nc' not in f and '.gz' not in f and '00000' not in f]
        
    if d == 'era5_2':
        cleaned = [ f for f in lista if '.conv._' in f and '.nc' not in f and '00000' not in f and '.gz' in f ]
        
    if d == 'amma' or d == 'giub' or d=='hara' or d=='npsound' or d=='shipsound':
        cleaned = [ f for f in lista if '.csv' in f ]

    return cleaned



def chunk_it(seq, num):
    """ Creates chunks of the input file list """
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out


def get_all_stations_files(db, station_kind, station_config_dir = '../data/station_configuration/'):
    """ Read the extended station configuration, to obtain the primary_id and the file name lists """
    
    if 'mobile' in db:
        station_kind = 'mobile'
        db = db.replace('_mobile' , '')
    if db in ['npsound' , 'shipsound']:
        station_kind = ''
        
    if station_kind in ['orphan']:
        sc = pd.read_csv( f'{station_config_dir}' + db + '_orphans_station_configuration_extended.csv', sep = '\t' )
    elif station_kind in ['mobile']:
        sc = pd.read_csv( f'{station_config_dir}' + db + '_mobile_station_configuration_extended.csv', sep = '\t' )
    else:
        sc = pd.read_csv( f'{station_config_dir}' + db + '_station_configuration_extended.csv', sep = '\t' )
        
    if not sc.empty:
        sc = sc[['primary_id' , 'file', 'latitude', 'longitude']]
        if db == 'giub':
            files = [f.replace('_reduced','') for f in sc.file ]
            sc['file'] = files

    ### filtering
    if station_kind == 'mobile':
        ids = [p for p in sc.primary_id.values if '20999' in p ]
        sc =sc.loc[sc.primary_id.isin(ids)]
    elif station_kind == 'orphan':
        ids = [p for p in sc.primary_id.values if '20999' in p ]
        sc =sc.loc[~sc.primary_id.isin(ids)]
        a =0        
        
    return sc


# a = get_all_stations_files('era5_1', 'regular')

def check_processed_stations(db, out_dir):
    """ Extract the list of stations that have already been processed with the list of years """
    
    processed  = {}
    
    if os.path.isdir(out_dir + '/' + db ):
        processed_stations = os.listdir(out_dir + '/' + db)
        
        for s in processed_stations:
            processed[s] = ''
            path = out_dir + '/' + db + '/' + s
            correctly_proc = [path+'/'+f for f in os.listdir(path ) if 'correctly' in f ]
            if len(correctly_proc) >0:
                
                lines = open(correctly_proc[0] , 'r').readlines()
                if len(lines) >0:
                    years = [l.split('\n')[0] for l in lines ]
            
                else:
                    years = []
            else:
                years = []
            processed[s] = years 
        
    return processed



'''
p = check_processed_stations(db, out_dir)

# at least one valid year
processed_stat = [s for s in p.keys() if len(p[s]) > 0 ]

all_files = get_all_stations_files(db)


print(0)
'''
