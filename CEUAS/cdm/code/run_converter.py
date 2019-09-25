import os,sys
import argparse
from build_311c_cdmfiles_ALL_split import db
import numpy as np
import subprocess

from check_correct import *


""" Here I select the databases, split the files into the number of wanted processes and run in parallel """

def filelist_cleaner(lista, d=''):
       """ Removes unwanted files that might be present in the database directories """
       print('Cleaning the list of files to be converted')
       if d == 'ncar':
          cleaned = [ l for l in lista if '.nc' not in l ]
       if d == 'bufr':
          cleaned = [ l for l in lista if '.bfr' in l and 'era5.' in l and '+100-' not in l and 'undef' not in l]
       if 'era5' in d:
          cleaned = [ l for l in lista if '.nc' not in l and '.conv.' in l ]
       if d =='igra2':
          cleaned = [ l for l in lista if '-data' in l ]
       if d == 'era5_3188':
          cleaned = [ l for l in lista if '.conv.C' in l and '.nc' not in l ]
       if d == 'era5_1761':
          cleaned = [ l for l in lista if '.conv.' in l and '.nc' not in l ]

       cleaned = [c for c in cleaned if '.py' not in c and '.nc' not in c ]
       return cleaned

def chunk_it(seq, num):
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out


out_dir = '/raid60/scratch/federico/1_23Sept'
processes = 25

""" Select the dataset to be processed """ 
#datasets = db.keys() # select the dataset
#datasets = 'test' (string) or ['ncar'] (list)
datasets = ['era5_1'] 

check_missing = False

if datasets == 'all':
    datasets = db.keys()

if datasets == 'test':
    os.system('/opt/anaconda3/bin/python3 build_311c_cdmfiles_ALL_split.py -d test -o TEST_CHECK ')
 
else:
  for d in datasets:
    print ('DATASET IS', d )
    files_list = [ db[d]['dbpath'] + '/' + f for f in os.listdir(db[d]['dbpath']) if os.path.isfile( db[d]['dbpath']+'/'+f ) ] # extracting the files list stores in the database path                   
    f_list = [ f for f in files_list if os.path.getsize(f) > 1 ] # cleaning the list of files in the original database directories                                                               
    #files_list = open('/raid60/scratch/federico/ncar_MISS.txt', 'r').readlines()
    f_list = filelist_cleaner(f_list, d = d)
    f_list = [ f.replace('\n','')  for f in f_list ]
    if check_missing == True:
           f_list = finder(converted,source)

    #print('TOTAL NUMBER MISSING::: ', len(f_list) )     
    chunks = chunk_it(f_list, processes)
    print('+++++++++ NUMBER OF files', len(f_list) )
    #input('check')
    #print(f_list)
    for c in chunks:
           print ('*** I am running CHUNK: ', chunks.index(c) , ' ***' )
           print('I have ', len(c), ' files ' , )
           c = str(','.join(c)).replace('#','')
           
           #os.system('/opt/anaconda3/bin/python3  build_311c_cdmfiles_ALL_split.py -d ' + d + ' -o ' + out_dir + ' -f ' + c + ' & ')     
           comm = '/opt/anaconda3/bin/python3  build_311c_cdmfiles_ALL_split.py -d ' + d + ' -o ' + out_dir + ' -f ' + c + ' & '
           subprocess.call( comm, shell = True )


print('\n\n\n *** Finished with the parallel running ***')
