import glob

from harvester_yearsplit_parameters import *  # datasets
from run_harvester_yearsplit_utils import *
from tqdm import tqdm

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


out_dir = '/scratch/das/federico/HARVEST_TEST_exampleStations_new_write_dict/'
datasets = datasets_big


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
                                          if db == 'era5_1':
                                                 fname = fname.replace('??????.','_').replace('.txt','')
                                          if db == 'era5_1761':
                                                 fname = fname.replace('._', '_').replace('.txt','')
                                                 
                                          completed_files.append(fname)
                                                 
                
                     file_df =file_df.loc[~file_df.file.isin(completed_files) ]                 
       

       
       files = list(file_df['file'].values)
       all_f = list(file_df['file'].values)


       test_stations = ['0-20000-0-82930' , '0-20000-0-26258' ]
       station_file = [  file_df['primary_id'].values[i] + '/' +  file_df['file'].values[i]  for i in range(len(file_df)) ] 
       station_file = [f for f in station_file for s in test_stations if s in f ]

       print(' ===== Running dataset ' + db + ' with ' + str(len(station_file) ) ) 
       #station_file = station_file[:20]
       chunks = chunk_it(station_file, processes)
       
       for c in chunks:
              print ('*** I am running CHUNK: ', chunks.index(c) , ' *** with: ' ,  len(c), ' files'  )
              c = str(','.join(c)).replace('#','')
              command = 'python3.8  harvest_convert_to_netCDF_yearSplit.py' 
              os.system( command + ' -d ' + db + ' -o ' + out_dir + ' -f ' + c + '  -r ' + str(run_only_missing_stations) +   '   & ' )
              
