import glob

from harvester_yearsplit_parameters import *  # datasets
from run_harvester_yearsplit_utils import *
from tqdm import tqdm
import random

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

all_stations  ={}

for db in datasets:

       # total number of files in the original datasets 
       file_df = get_all_stations_files(db)
 
       # select stations from a list 
       if stations:
              print("+++ Only running selected stations::: " , stations )
              file_df =file_df.loc[file_df.primary_id.isin(stations) ]                 

              
       ## extracting fully harvested files i.e. all available years were harvested  
       ### these files will not be checked anymore against the year during the harvesting 
       if skip_fully_harvested:
              if  os.path.isdir(out_dir+'/' + db):
                     stats = os.listdir(out_dir+'/' + db )
                     completed_files = []  # fully harvested files (all completed years)
                     
                     for s in tqdm(stats):
                            path= out_dir+'/' + db  + '/' + s 
                            files = [ f for f in os.listdir(path) if 'correctly_processed' in f ]
                            for f in files:
                                   f_name = f.split('harvested_')[1].split('_correctly')[0]
                                   lines = open(out_dir+'/' + db  + '/' + s + '/' + f , 'r').readlines()
                                   if lines[-1] == 'completed\n':
                                          fname = f_name.replace('.gz','')
                                          if db  =='era5_1':
                                                 fname = fname.replace('??????.','_').replace('.txt','')
                                          if db == 'era5_1_mobile':
                                                 fname = fname.replace('??????.','').replace('.txt','')
                                                 
                                          if db == 'era5_1761':
                                                 fname = fname.replace('._', '_').replace('.txt','')
                                                 
                                          completed_files.append(fname)
                                                 
                
                     file_df =file_df.loc[~file_df.file.isin(completed_files) ]                 
       
       '''
       df_red =file_df.loc[ ~file_df['primary_id'].isin(primary_proc)]
       file_df = df_red
       '''
       
       files = list(file_df['file'].values)
       all_f = list(file_df['file'].values)


       station_file = [  file_df['primary_id'].values[i] + '/' +  file_df['file'].values[i]  for i in range(len(file_df)) ] 
       files = station_file
       
       #primary_file_dic =dict(zip(  list(all_files['file'].values) ,  list(all_files['primary_id'].values) ))
       
       #files = [f for f in files if '11035' in f or '10393' in f or '9393' in f or '10395' in f or '06610' in f or '82930' in f ]

       #files = [f for f in all_f if '06610' in f  ]
       #random.shuffle(station_file)       
       #station_file=station_file[:10]

       #station_file = [ s for s in station_file if '10393' in s ]
       
       print(' ===== Running dataset ' + db + ' with ' + str(len(station_file) ) ) 
       #station_file = station_file[:20]
       chunks = chunk_it(station_file, processes)
       
       for c in chunks:
              print ('*** I am running CHUNK: ', chunks.index(c) , ' *** with: ' ,  len(c), ' files'  )
              c = str(','.join(c)).replace('#','')
              command = 'python3.8  harvest_convert_to_netCDF_yearSplit.py' 
              os.system( command + ' -d ' + db + ' -o ' + out_dir + ' -f ' + c + '  -r ' + str(run_only_missing_stations) +   '   & ' )
              
              
       a=0
