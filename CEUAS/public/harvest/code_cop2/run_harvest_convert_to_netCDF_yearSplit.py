import glob

from harvester_yearsplit_parameters import *  # datasets
from run_harvester_yearsplit_utils import *


#from check_correct import *

'''
# old run for large files without python3.8
tcsh
conda env list
conda activate py38
"""


"""
REMEMBER
- era5_1 and era5_1_mobile must use the monthly data files
- era5_2, despite having those files, can be run using directly the global text file 
'''



### make list of stations to run 


out_dir = '/scratch/das/federico/HARVEST_YEARLY_17AUG2023/'

all_stations  ={}


datasets = ['ncar' , 'igra2', 'era5_1759' , 'era5_1761']

#datasets = ['era5_1759' , 'era5_1761']

processes = 2

datasets = ['era5_2']

for db in datasets:
       
       ### all processed stations with years       
       #p = check_processed_stations(db, out_dir)
       #all_stations[db] = p
       
       #processed_stat = [s for s in p.keys() if len(p[s]) > 0 ]

       # total number of files in the original datasets 
       all_files = get_all_stations_files(db)
       
       files = list(all_files['file'].values)
       import random
       random.shuffle(files)       
       files=files[:10]
       
       chunks = chunk_it(files, processes)
       
       for c in chunks:
              print ('*** I am running CHUNK: ', chunks.index(c) , ' *** with: ' ,  len(c), ' files'  )
              c = str(','.join(c)).replace('#','')
              os.system('python3.8  harvest_convert_to_netCDF_yearSplit.py  -d ' + db + ' -o ' + out_dir + ' -f ' + c + '  -r ' + str(run_only_missing_stations) +   '  & ')
              
       print(0)

