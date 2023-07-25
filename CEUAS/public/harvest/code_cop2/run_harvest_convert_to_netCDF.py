import os,sys
import argparse
import numpy as np
import subprocess
import glob
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
- era5_2, despite having those files, can be run using directly the global ODB file 
'''

# on srvx1, srvx8 

# must be euqal to the db in the harvester script but cannot import it from the other file
# due to circular import i.e. out_dir is defined here and imported in the harvester
# but the db are better to be defined in the harvester
# so must check that these dictionaries are equal

db = { 'era5_1': '/mnt/users/scratch/leo/scratch/era5/odbs/1' ,
       'era5_1_mobile': '/mnt/users/scratch/leo/scratch/era5/odbs/1_mobile' ,
       
       
       'era5_2': '/mnt/users/scratch/leo/scratch/era5/odbs/2',
       'era5_2_mobile': '/mnt/users/scratch/leo/scratch/era5/odbs/2',
       
       
       'era5_3188': '/mnt/users/scratch/leo/scratch/era5/odbs/3188',
       'era5_1759': '/mnt/users/scratch/leo/scratch/era5/odbs/1759',
       'era5_1761': '/mnt/users/scratch/leo/scratch/era5/odbs/1761',
       'bufr': '/mnt/users/scratch/leo/scratch/era5/odbs/ai_bfr/',                                   
       'ncar': '/scratch/das/federico/databases_service2/UADB_20230109/', # check if extra '/' needed 

       'igra2': '/scratch/das/federico/databases_service2/IGRA2_20230106',
       
       'amma': '/scratch/das/federico/databases_service2/AMMA_BUFR/AMMA_split_csv/',

       'giub':'/scratch/das/federico/redo_GIUB_07072023_reduced/',

       'hara':'/scratch/das/federico/databases_service2/HARA-NSIDC-0008_csv/',
       'npsound' : '/scratch/das/federico/databases_service2/NPSOUND-NSIDC0060/NPSOUND_csv_converted',
       'shipsound' : '/scratch/das/federico/databases_service2/SHIPSOUND-NSIDC-0054/' ,

} 


#from run_harvest_convert_to_netCDF import db

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
       """ Creates sub sets of the input file list """
       avg = len(seq) / float(num)
       out = []
       last = 0.0
   
       while last < len(seq):
              out.append(seq[int(last):int(last + avg)])
              last += avg
   
       return out

# deine output directory "out_dir"


out_dir = '/scratch/das/federico/HARVEST_YEARSPLIT_25072023/'


""" Select the dataset to be processed """ 
datasets = ['era5_1', 'era5_2', 'era5_3188', 'era5_1759', 'era5_1761', 'ncar', 'igra2', 'bufr' ]
all_ds = ['era5_1',
          'era5_2',
          'era5_1_mobile',
          'era5_2_mobile',
          'era5_3188',
          'era5_1759',
          'era5_1761',
          'ncar',
          'igra2',
          'bufr',
          'amma',
          'npsound',
          'shipsound',
          'hara',
          'giub'
]


datasets = all_ds

datasets = ['bufr']

""" Numbe rof processed PER DATASET !!! """
processes = 20                                                                                                                                                                       

""" Check processed files """
check_missing = True
REDO = False


def rerun_list(f_list, processed_dir = '', input_dir = ''):
       try:
              processed = [ f.split('harvested_')[1].replace('.nc','') for f in os.listdir(processed_dir)  if 'harvested' in f ]

       except:
              processed = []
              
       f_list          = [ f.split(input_dir)[1] for f in f_list ]
       
       to_be_processed = []
       for file in f_list:
              file = file.replace('/','')
              if file in processed:
                     continue
              else:
                     to_be_processed.append(input_dir + '/' + file )
       #print('total to be processed: ', len(to_be_processed) )
       return to_be_processed

if not os.path.isdir('logs/'):
       os.mkdir('logs')
       
for d in datasets:
       processed_dir = out_dir + '/' + d 
       print ('DATASET IS', d )
       if d not in ['era5_1', 'era5_1_mobile', 'era5_2', 'era5_2_mobile']:
              files_list = [ db[d]+ '/' + f for f in os.listdir(db[d]) if os.path.isfile( db[d]+'/'+f ) ] # extracting the \
              f_list = [ f for f in files_list if os.path.getsize(f) > 1 ] # cleaning the list of files in the original database directories                                                               
              f_list = filelist_cleaner(f_list, d = d)
              f_list = [ f.replace('\n','')  for f in f_list ]
              #print("#### Number of files in the original directory::: " , len(f_list) , '      ', d)
              if check_missing :
                     f_list = rerun_list(f_list, processed_dir = processed_dir ,  input_dir =  db[d]  )              
                     #print(' Removed already processed #### ')
       else:
              processes = 3 # !!! there is already the pool splitting in the harvester, cannot be higher due to memory issues 
              # era5.conv.??????.82930.txt.gz
              if not os.path.isdir('files_list'):
                     os.mkdir('files_list')
                     
              if not os.path.isfile('files_list/' + d + '_files_list.txt'):  ### TODO, this is to speed ud reading the file names which take a lot of time
                     if d == 'era5_1':
                            flist=glob.glob('files_list/' + db[d] + "/era5.conv._*") # takes too long
                     elif d == 'era5_1_mobile':
                            #I cannot use the same approach cause I do not have single global station files...
                            flist = glob.glob(db[d] + "/era5.conv.??????._*")  
                            flist = np.unique( ['._' + f.split('._')[1] for f in flist if '_' in f  ])
                     #flist =[f for f in flist if '_40179' not in f  and '42147' not in f] # this file causes problems ???
                     a = open( 'files_list/' + d + '_files_list.txt','w')
                     for l in flist:
                            a.write(l + '\n')
                     a.close()                     
                
              else:
                     a = open('files_list/' + d + '_files_list.txt', 'r').readlines()
                     flist = [ f.replace('\n','').replace(' ','') for f in a]

              #odbs = flist
              #stat = [f.split('._')[1] for f in odbs ] # numerical codes for station e.g. '07630'
              
              # create names of merged files to check in the directory 
              if d == 'era5_1_mobile':
                     harvestable_stations = ['era5.conv.??????._' + s.split('._')[1] + '.nc' for s in flist ]
                     print(flist[:2], '\n\n', harvestable_stations[:2] )

              elif d=='era5_1':
                     harvestable_stations  = ['era5.conv.??????.'+ s.split('._')[1] + '.txt.gz.nc'  for s in flist ]

                     print(flist[:2], '\n\n', harvestable_stations[:2] )
              elif d in ['era5_2' ]:
                     harvestable_stations = ['era5.conv._' + s.split('._')[1] + '.gz.nc'   for s in flist ]
                     
              elif d in ['era5_2_mobile']:
                     
                     harvestable_stations = [ f.split('/')[-1] + '.gz'   for f in flist ]

                     
              if check_missing:
                     if os.path.isdir(processed_dir):
                            processed = [s.split('_harvested_')[1] for s in os.listdir(processed_dir) if 'harvested' in s]
                            print('PROCESSED ' , processed[:2] )
                            f_list = [f for f in harvestable_stations if f not in processed ]
                     else:
                            f_list = harvestable_stations
              else:
                     f_list = harvestable_stations
                     
              #if d in ['era5_1_mobile']:
              #       f_list = ['"' + db[d] + '/' + s.replace('.nc', '').replace('?.' , "?._")+ '"' for s in f_list ] # adjust names with "" symbols and remove .nc extension
              if d in ['era5_1']:
                     f_list = ['"' + db[d] + '/' + s.replace('.nc', '') + '"' for s in f_list ] 
              else:
                     f_list = [ db[d] + '/' + s.replace('.nc', '') for s in f_list ] 
                     


       # harvest random files per process 
       import random
       random.shuffle(f_list)
       f_list = [f for f in f_list if '#' not in f ]
       #f_list = f_list[:5]
       print('PROCESSING ::: ' , f_list )
       # divide into chunks
       chunks = chunk_it(f_list, processes)
       print('+++++++++ TOTAL NUMBER OF FILES to be reprocessed ::: ', len(f_list) )
       for c in chunks:
              print ('*** I am running CHUNK: ', chunks.index(c) , ' *** with: ' ,  len(c), ' files'  )
              c = str(','.join(c)).replace('#','')
              os.system('python3.8  harvest_convert_to_netCDF_yearSplit.py  -d ' + d + ' -o ' + out_dir + ' -f ' + c + ' & ')                                                                              

print('*** Finished with the parallel running ***')
