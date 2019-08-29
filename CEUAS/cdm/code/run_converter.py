import os,sys
import argparse
from build_311c_cdmfiles_ALL_split import db
import numpy as np


""" Here I select the databases, split the files into the number of wanted processes and run in parallel """



def filelist_cleaner(lista, dataset=''):
       """ Removes unwanted files that might be present in the database directories """
       print('Cleaning the list of files to be converted')
       if dataset == 'ncar':
          cleaned = [ l for l in lista if '.nc' not in l ]
       if dataset == 'bufr':
          cleaned = [ l for l in lista if '.bfr' in l ]
       if 'era5' in dataset:
           cleaned = [ l for l in lista if '.nc' not in l and '.conv.' in l ]
       else:
           cleaned = lista

       cleaned = [c for c in cleaned if '.py' not in c ]
       return cleaned

def chunk_it(seq, num):
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out



def find_missing(processed_dir, original_dir, db = 'ncar'):
    """ Given input and  output directories, finds the files that have not been processed and return a list """
    miss = []
    proc = [ f for f in os.listdir(processed_dir) if '.py' not in f ]
    ori  = [ f for f in os.listdir(original_dir) if '.py' not in f ]

    for f in ori:
        if db == 'ncar': name = 'ch' + f + '.nc'
        if name not in proc:
            #print('file not processed! ', name )                                                                                                                                                                                                                                         
           miss.append(original_dir + '/' + f )
    print(' Total number of files missing::: ', str(len(miss)) )
    return miss


out_dir = '/raid60/scratch/federico/All_odb'
processes = 15

""" Select the dataset to be processed """ 
#datasets = db.keys() # select the dataset
#datasets = 'test'
datasets = ['era5_1']

check_missing = False

if datasets == 'all':
    datasets = db.keys()

if datasets == 'test':
    os.system('/opt/anaconda3/bin/python3 build_311c_cdmfiles_ALL_split.py -d test -o TEST_CHECK ')
 
else:
 if not check_missing:
   for d in datasets:

    files_list = [ db[d]['dbpath'] + '/' + f for f in os.listdir(db[d]['dbpath']) if os.path.isfile(db[d]['dbpath'] + '/' + f)] # extracting the files list stores in the database path                                                      \
    files_list = [ f for f in files_list if os.path.getsize(f) > 10 ] # cleaning the list of files in the original database directories                                                               
    print (files_list)                        
    files_list = filelist_cleaner(files_list, d)

    #chunks = [files_list[x:x+100] for x in range(0, len(files_list), processes)] # creates "processes                            
    
    chunks = chunk_it(files_list, processes)

    for c in chunks:
            
            c = str(','.join(c))
            #print(c)
            #input('hello')
            os.system('/opt/anaconda3/bin/python3 build_311c_cdmfiles_ALL_split.py -d ' + d + ' -o ' + out_dir + ' -f ' + c + ' & ')                      

 elif check_missing:
  for d in datasets:
         files_list = find_missing( out_dir + '/' + d , db[d]['dbpath'] , db = d )
         #print ('Processing the files: ', files_list )
         files_list = filelist_cleaner(files_list, d)
         chunks = chunk_it(files_list, processes)
         for c in chunks:
              c = str(','.join(c))
              out = out_dir + '_missing'
              os.system('/opt/anaconda3/bin/python3 build_311c_cdmfiles_ALL_split.py -d ' + d + ' -o ' + out + ' -f ' + c + ' & ')


print('\n\n\n *** Finished with the parallel running ***')
