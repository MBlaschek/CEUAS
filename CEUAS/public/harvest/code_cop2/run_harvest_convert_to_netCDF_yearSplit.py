import glob
import sys
from run_harvester_yearsplit_utils import *
from tqdm import tqdm
import random
from random import shuffle
import argparse
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

parser = argparse.ArgumentParser(description="Crete station configuration table")

parser.add_argument('--station_config_dir' , '-s', 
                    help="Chose station configuration directory"  ,
                    type = str,
                    default = '.' )

parser.add_argument('--parameter_file' , '-p', 
                    help="Chose parameter file"  ,
                    type = str,
                    default = '.' )

parser.add_argument('--ceuas_dir' , '-c', 
                    help="Chose scrip directory"  ,
                    type = str,
                    default = '/srvfs/home/uvoggenberger/CEUAS/CEUAS/' )



args = parser.parse_args()
station_config_dir = args.station_config_dir
parameter_file = args.parameter_file
ceuas_dir = args.ceuas_dir

if parameter_file != '.': 
       print(parameter_file)
       sys.path.append(parameter_file.split('modded')[0])
       from modded_harvester_yearsplit_parameters import *
else:
       from harvester_yearsplit_parameters import *  

global automatic_call
automatic_call = False
if station_config_dir != '.':
       automatic_call = True


### make list of stations to run 

all_stations  ={}

# datasets = ['era5_1_mobile','era5_2_mobile']


for db in datasets:

       # total number of files in the original datasets 
       if automatic_call:
              file_df = get_all_stations_files(db, station_kind, station_config_dir)
              dsp = datasets_path[datasets[0]]
              available = [i for i in file_df.file if os.path.isfile(dsp+'/'+i)]
              file_df = file_df[file_df.file.isin(available)]
       else:
              file_df = get_all_stations_files(db, station_kind)
       
       # select stations from a list 
       if stations:
              print("+++ Only running selected stations::: " , stations )
              file_df =file_df.loc[file_df.primary_id.isin(stations) ]       

       # # # DEBUGGING: to run only a single station:
       # file_df = file_df[file_df['primary_id'].str.contains('11035|484848')]

       ## extracting fully harvested files i.e. all available years were harvested  
       ### these files will not be checked anymore against the year during the harvesting 
       if skip_fully_harvested:
              if  os.path.isdir(out_dir+'/' + db):
                     stats = os.listdir(out_dir+'/' + db )
                     completed_files = []  # fully harvested files (all completed years)
                     
                     for s in tqdm(stats, miniters=100):
                            if 'empty' in s: continue                            
                            path= out_dir+'/' + db  + '/' + s 
                            files = [ f for f in os.listdir(path) if 'correctly_processed' in f ]
                            for f in files:
                                   f_name = f.split('harvested_')[1].split('_correctly')[0]
                                   lines = open(out_dir+'/' + db  + '/' + s + '/' + f , 'r').readlines()
                                   if 'completed' in lines[-1] :
                                          fname = f_name.replace('.gz','')
                                          if db  =='era5_1':
                                                 fname = fname.replace('??????.','_').replace('.txt','')
                                          if db == 'era5_1_mobile':
                                                 fname = fname.replace('??????.','').replace('.txt','')
                                                 
                                          if db in ['era5_1761','era5_1759']:
                                                 fname = fname.replace('._', '.').replace('.txt','')
                                          if db in ['igra2']:
                                                 fname = fname.replace('-data.txt', '')
                                                 
                                          completed_files.append(fname)
                                                 
                
                     file_df =file_df.loc[~file_df.file.isin(completed_files) ]  
       
       '''
       df_red =file_df.loc[ ~file_df['primary_id'].isin(primary_proc)]
       file_df = df_red
       '''
       
       files = list(file_df['file'].values)
       all_f = list(file_df['file'].values)


       if db == 'era5_1_mobile': # fix for new files
              station_file = [  file_df['primary_id'].values[i] + '/' +  file_df['file'].values[i].replace('_', '')  for i in range(len(file_df)) ] 
       else:
              station_file = [  file_df['primary_id'].values[i] + '/' +  file_df['file'].values[i]  for i in range(len(file_df)) ] 
       files = station_file
       print(files)
       
       #primary_file_dic =dict(zip(  list(all_files['file'].values) ,  list(all_files['primary_id'].values) ))
       
       # files = [f for f in files if '11035' in f or '10393' in f or '9393' in f or '10395' in f or '06610' in f or '82930' in f ]
       # files = [f for f in files if '11035' in f]

       #files = [f for f in all_f if '06610' in f  ]
       #random.shuffle(station_file)       
       #station_file=station_file[:10]

       #station_file = [ s for s in station_file if '10393' in s ]
       #station_file = station_file[:5]
       
       print(' ===== Running dataset ' + db + ' with ' + str(len(station_file) ) ) 
       # station_file = station_file[:20]  ## to test with a subset
       # chunks = chunk_it(station_file, processes)
       shuffle(files)
       chunks = chunk_it(files, processes)
       
       for c in chunks:
              #print ('*** I am running CHUNK: ', chunks.index(c) , ' *** with: ' ,  len(c), ' files'  )
              c = str(','.join(c)).replace('#','')
              command = f'python {ceuas_dir}/public/harvest/code_cop2/harvest_convert_to_netCDF_yearSplit.py' # changed from python3.8
              com =  command + ' -d ' + db + ' -o ' + out_dir + ' -f  ' + c + '  -r ' + str(run_only_missing_stations) +   '  -k ' + station_kind + ' -pf ' + parameter_file + ' & '
              print(com)
              if len(c) > 0:
                     os.system( com )
              
              
       a=0


'''
era5_2 ['0-20500-0-3869/era5.conv._03869', '0-20500-0-93126/era5.conv._2:93126', '0-20500-0-3869/era5.conv._2:3869']
'''
