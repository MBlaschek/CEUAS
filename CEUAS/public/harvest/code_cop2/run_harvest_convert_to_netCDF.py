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
'''


# on srvx1, srvx8 
db = { 'era5_1': '/mnt/users/scratch/leo/scratch/era5/odbs/1' ,
       'era5_2': '/mnt/users/scratch/leo/scratch/era5/odbs/2',
                               'era5_3188': '/mnt/users/scratch/leo/scratch/era5/odbs/3188',
                               'era5_1759': '/mnt/users/scratch/leo/scratch/era5/odbs/1759',
                               'era5_1761': '/mnt/users/scratch/leo/scratch/era5/odbs/1761',

                               'bufr': '/mnt/users/scratch/leo/scratch/era5/odbs/ai_bfr/',                                   
                               'ncar': '/scratch/das/federico/databases_service2/UADB_25012022/',
                               'igra2': '/scratch/das/federico/databases_service2/IGRA2_20211231/', 

                               } 


""" Select the databases, split the files into the number of wanted processes and run in parallel """
def filelist_cleaner(lista, d=''):
       """ Removes unwanted files that might be present in the database directories """
       print('Cleaning the list of files to be converted')
       if d == 'ncar':
              cleaned = [ l for l in lista if '.nc' not in l and '.py' not in l ]
       if d == 'bufr':
              cleaned = [ l for l in lista if '.bfr' in l and 'era5.' in l and '+100-' not in l and 'undef' not in l]
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

# deifne output directory "out_dir"
out_dir = '/scratch/das/federico/COP2_HARVEST_FEB2022/'


processes = 6 # number of process PER DATASET 


""" Select the dataset to be processed """ 
datasets = ['era5_1', 'era5_2', 'era5_3188', 'era5_1759', 'era5_1761', 'ncar', 'igra2', 'bufr' ]
datasets = ['era5_1761', 'era5_2', 'era5_3188', 'era5_1759' ]

datasets = ['era5_1']



""" Check processed files """
check_missing = True
REDO = False

def rerun_list(f_list, processed_dir = '', split = '' , input_dir = ''):
       try:
              processed = [ f.split('harvested_')[1].replace('.nc','') for f in os.listdir(processed_dir)  if 'harvested' in f ]

       except:
              processed = []
              
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
       print(to_be_processed)
       return to_be_processed


for d in datasets:
       processed_dir = out_dir + '/' + d 
       print ('DATASET IS', d )
       #files_list = [ db[d]['dbpath'] + '/' + f for f in os.listdir(db[d]['dbpath']) if os.path.isfile( db[d]['dbpath']+'/'+f ) ] # extracting the files list stores in the database path                   
       if d != 'era5_1':
              files_list = [ db[d]+ '/' + f for f in os.listdir(db[d]) if os.path.isfile( db[d]+'/'+f ) ] # extracting the \
              f_list = [ f for f in files_list if os.path.getsize(f) > 1 ] # cleaning the list of files in the original database directories                                                               
              f_list = filelist_cleaner(f_list, d = d)
              f_list = [ f.replace('\n','')  for f in f_list ]
              print("#### Number of files in the original directory::: " , len(f_list) , '      ', d)
              if check_missing :
                     f_list = rerun_list(f_list, processed_dir = processed_dir , split = 'ai_bfr' , input_dir =  db[d]  )              
                     #print(' Removed already processed #### ')
       else:
              processes = 6 # !!! there is already the pool splitting in the harvester, cannot be higher due to memory issues 
              # era5.conv.??????.82930.txt.gz
              odbs = glob.glob(db[d] + '/' + 'era5.conv._*')
              
              stat =  [f.split('._')[1] for f in odbs ]
              
              f_list = ['"' + db[d] + '/era5.conv.??????.' + s + '.txt.gz' + '"' for s in stat]
              if check_missing:
                     processed = [s.split('_harvested_')[1] for s in os.listdir(processed_dir) if 'harvested' in s]
                     print(processed[:10])
                     f_list = [f for f in f_list if f.split('/1/')[1].replace('"','')+'.nc' not in processed ]
                     print(f_list[:10])
            
       chunks = chunk_it(f_list, processes)
       
       print('+++++++++ TOTAL NUMBER OF FILES to be reprocessed ::: ', len(f_list) )
       for c in chunks:
              #print ('*** I am running CHUNK: ', chunks.index(c) , ' *** with: ' ,  len(c), ' files'  )
              c = str(','.join(c)).replace('#','')
        
              if d != 'era5_1':
                     print(" ")
                     os.system('python3.8  harvest_convert_to_netCDF.py  -d ' + d + ' -o ' + out_dir + ' -f ' + c + ' & ')  
              else:
                     print(" ")
                     #print(f_list[:10])
                     os.system('python3.8  harvest_convert_to_netCDF.py  -d ' + d + ' -o ' + out_dir + ' -f ' + c + ' & ')
              

print('*** Finished with the parallel running ***')
