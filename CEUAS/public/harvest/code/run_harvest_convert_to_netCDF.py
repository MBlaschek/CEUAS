import os,sys
import argparse
from harvest_convert_to_netCDF import db
import numpy as np
import subprocess

#from check_correct import *
'''
tcsh
conda env list
conda activate py38
'''

""" Select the databases, split the files into the number of wanted processes and run in parallel """
def filelist_cleaner(lista, d=''):
       """ Removes unwanted files that might be present in the database directories """
       print('Cleaning the list of files to be converted')
       if d == 'ncar':
              cleaned = [ l for l in lista if '.nc' not in l ]
       if d == 'bufr':
              cleaned = [ l for l in lista if '.bfr' in l and 'era5.' in l and '+100-' not in l and 'undef' not in l]
       if d in ['era5_1759', 'era5_1761']:
              cleaned = [ l for l in lista if '.nc' not in l and '.conv.' in l and 'py' not in l and '.gz' in l and '.g.' not in l ] 
       if d =='igra2':
              cleaned = [ l for l in lista if '-data' in l ]
       if d == 'era5_3188': #  era5.3188.conv._C:4687.gz
              cleaned = [ l for l in lista if '3188.conv' in l and '.nc' not in l and '.gz' in l]
       if d == 'era5_1':
              cleaned = [ f for f in lista if '.conv._' in f and '.nc' not in f and '.gz' not in f and '00000' not in f]
       if d == 'era5_2':
              cleaned = [ f for f in lista if '.conv._' in f and '.nc' not in f and '00000' not in f and '.gz' in f ]
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


out_dir = '/raid60/scratch/federico/MARCH2021_WBAN_latswap_era5_1759/'

processes = 20 # number of process PER DATASET 



""" Select the dataset to be processed """ 
datasets = ['era5_1', 'era5_2', 'era5_3188', 'era5_1759', 'era5_1761', 'ncar', 'igra2', 'bufr' ]

datasets = ['era5_1759']



""" Check processed files """
check_missing = True
REDO = False

def rerun_list(f_list, processed_dir = '', split = '' , input_dir = ''):
       
       processed = [ f.split('harvested_')[1].replace('.nc','') for f in os.listdir(processed_dir)  ]
       f_list          = [ f.split(input_dir)[1] for f in f_list ]
       
       to_be_processed = []
       for file in f_list:
              file = file.replace('/','')
              if file in processed:
                     #print('skipping ' , file)
                     continue
              else:
                     to_be_processed.append(input_dir + '/' + file )
       #print('total to be processed: ', len(to_be_processed) )
       return to_be_processed



for d in datasets:
       processed_dir = out_dir + '/' + d 
       print ('DATASET IS', d )
       #files_list = [ db[d]['dbpath'] + '/' + f for f in os.listdir(db[d]['dbpath']) if os.path.isfile( db[d]['dbpath']+'/'+f ) ] # extracting the files list stores in the database path                   
       if d != 'era5_1':
              files_list = [ db[d]['dbpath'] + '/' + f for f in os.listdir(db[d]['dbpath']) if os.path.isfile( db[d]['dbpath']+'/'+f ) ] # extracting the \
              f_list = [ f for f in files_list if os.path.getsize(f) > 1 ] # cleaning the list of files in the original database directories                                                               
              f_list = filelist_cleaner(f_list, d = d)
              f_list = [ f.replace('\n','')  for f in f_list ]
              print("#### Number of files in the original directory::: " , len(f_list) , '      ', d)
              if check_missing :
                     f_list = rerun_list(f_list, processed_dir = processed_dir , split = 'ai_bfr' , input_dir =  db[d]['dbpath']  )              
                     #print(' Removed already processed #### ')
       else:
              processes = 5 # !!! there is already the pool splitting in the harvester !!!
              Dir = '/raid60/scratch/leo/scratch/era5/odbs/1/era5_1'
              # era5.conv.??????.82930.txt.gz
              stat =  [ f.split('era5.conv.')[1].split('.txt')[0] for f in os.listdir(Dir) if 'harvested' in f ]
              already_proc = [f.split('??????.')[1].split('.txt')[0] for f in os.listdir('/raid60/scratch/federico/TEST_MAY_ERA5/era5_1') ]
              stat = [s for s in stat if s not in already_proc ]
              f_list = ['"/raid60/scratch/leo/scratch/era5/odbs/1/era5.conv.??????.' + s + '.txt.gz' + '"' for s in stat]
       #f_list = f_list[:100]

     
            
       chunks = chunk_it(f_list, processes)
       
       print('+++++++++ TOTAL NUMBER OF FILES to be reprocessed ::: ', len(f_list) )
       for c in chunks:
              #print ('*** I am running CHUNK: ', chunks.index(c) , ' *** with: ' ,  len(c), ' files'  )
              c = str(','.join(c)).replace('#','')
        
              if d != 'era5_1':
                     print(1)
                     os.system('/opt/anaconda3/bin/python3  harvest_convert_to_netCDF.py  -d ' + d + ' -o ' + out_dir + ' -f ' + c + ' & ')  
              else:
                     print(0)
                     #os.system('python  harvest_convert_to_netCDF.py  -d ' + d + ' -o ' + out_dir + ' -f ' + c + ' & ')
              

print('*** Finished with the parallel running ***')
