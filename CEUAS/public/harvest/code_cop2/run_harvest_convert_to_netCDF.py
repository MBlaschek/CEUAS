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
       'ncar': '/scratch/das/federico/databases_service2/UADB_25012022', # check if extra '/' needed 
       'igra2': '/scratch/das/federico/databases_service2/IGRA2_20211231',
       'era5_1_mobile': '/mnt/users/scratch/leo/scratch/era5/odbs/1_mobile' ,
       'era5_2_mobile': '',
       'amma': '/scratch/das/federico/databases_service2/AMMA_BUFR/',
} 


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
       if d == 'amma':
              cleaned = [ f for f in lista if '.' not in f ]
              
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
out_dir = '/scratch/das/federico/COP2_HARVEST_NOVEMBER2022/'


processes = 40 # number of process PER DATASET 


""" Select the dataset to be processed """ 
datasets = ['era5_1', 'era5_2', 'era5_3188', 'era5_1759', 'era5_1761', 'ncar', 'igra2', 'bufr' ]

all_ds = ['era5_1', 'era5_2', 'era5_1_mobile',
          'era5_3188', 'era5_1759', 'era5_1761',
          'ncar', 'igra2', 'bufr', 'amma' ]

datasets = ['era5_1_mobile']

datasets = all_ds

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
       if 'era5_1_' not in d:
              files_list = [ db[d]+ '/' + f for f in os.listdir(db[d]) if os.path.isfile( db[d]+'/'+f ) ] # extracting the \
              f_list = [ f for f in files_list if os.path.getsize(f) > 1 ] # cleaning the list of files in the original database directories                                                               
              f_list = filelist_cleaner(f_list, d = d)
              f_list = [ f.replace('\n','')  for f in f_list ]
              #print("#### Number of files in the original directory::: " , len(f_list) , '      ', d)
              if check_missing :
                     f_list = rerun_list(f_list, processed_dir = processed_dir ,  input_dir =  db[d]  )              
                     #print(' Removed already processed #### ')
       else:
              processes = 10 # !!! there is already the pool splitting in the harvester, cannot be higher due to memory issues 
              # era5.conv.??????.82930.txt.gz

              if not os.path.isfile('logs/' + d + '_files_list.txt'):  ### TODO, this is to speed ud reading the file names which take a lot of time
                     if d == 'era5_1':
                            flist=glob.glob('logs/' + db[d] + "/era5.conv._*") # takes too long
                     elif d == 'era5_1_mobile':
                            #I cannot use the same approach cause I do not have single global station files...
                            flist = glob.glob(db[d] + "/era5.conv.??????._*")  
                            flist = np.unique( ['._' + f.split('._')[1] for f in flist if '_' in f  ])
                     #flist =[f for f in flist if '_40179' not in f  and '42147' not in f] # this file causes problems ???
                     a = open( 'logs/' + d + '_files_list.txt','w')
                     for l in flist:
                            a.write(l + '\n')
                     a.close()
                
              else:
                     flist = [ f.replace('\n','') for f in open(d + '_files_list.txt').readlines() ]

              odbs = flist

              stat = [f.split('._')[1] for f in odbs ]
              if d == 'era5_1_mobile':
                     f_list = ['"' + db[d] + '/era5.conv.??????._' + s.replace('.txt.gz', '') + '.txt.gz' + '"' for s in stat ]
              else:
                     f_list = ['"' + db[d] + '/era5.conv.??????.' + s.replace('.txt.gz', '') + '.txt.gz' + '"' for s in stat ]
              if check_missing:
                     try:
                            processed = [s.split('_harvested_')[1] for s in os.listdir(processed_dir) if 'harvested' in s]
                            if d == 'era5_1':
                                   
                                   f_list = [f for f in f_list if f.split('/1/')[1].replace('"','')+'.nc' not in processed ]
                            else:
                                   f_list = [f for f in f_list if f.split('/1_mobile/')[1].replace('"','')+'.nc' not in processed ]
                                   
                     except:
                            pass

       # filter failed 
       # filt = [ f.replace('\n','') for f in open('logs/' + d + '_failed.txt').readlines() ]
       # f_list = [f for f in f_list if f not in filt ]
       print(len(f_list))
       #f_list = f_list[:5]
       if 'mobile'in d:
              f_list = [f for f in f_list if '_' in f ]
              # done to filter out aggregated yearly files...
       print("FLIST::: \n",   f_list )
       import random
       random.shuffle(f_list)
       chunks = chunk_it(f_list, processes)
       print('+++++++++ TOTAL NUMBER OF FILES to be reprocessed ::: ', len(f_list) )
       for c in chunks:
              print ('*** I am running CHUNK: ', chunks.index(c) , ' *** with: ' ,  len(c), ' files'  )
              c = str(','.join(c)).replace('#','')
              #os.system('python3.8  harvest_convert_to_netCDF.py  -d ' + d + ' -o ' + out_dir + ' -f ' + c + ' & ')                                                                              

print('*** Finished with the parallel running ***')
