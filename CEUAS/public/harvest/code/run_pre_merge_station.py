""" Utility to run the pre-merging script, running in parallel. """

import os,sys
import pandas as pd
import numpy as np
from pre_merge_stations import dic_dataset_stationsfiles, read_input_file 


def chunk_it(seq, num):
    """ Split the number of files into a number=num of chunks. """
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out



""" Will read all the stations to be combined from the input file;
      each set will be processed separately.
      Will split the combining process into the number of processes selected with the processes variable and run in parallel. """


datasets = ['era5_3188']  # note that it will create a total of (processes x #datasets) number of parallel runs! 

out_dir_base = '/raid60/scratch/federico/Pre_Merged_File'
processes = 25


for d in datasets:
    print('*** Started combining the ' + d + ' dataset ')
    out_dir = out_dir_base + '/' + d
    stations_file = dic_dataset_stationsfiles[d]['stations_file']
    source_dir = dic_dataset_stationsfiles[d]['source_dir']
    
    dic = read_input_file(stations_file , source_dir)
    
    stations = list( dic.keys() ) # to test 
    
    chunks = chunk_it(stations, processes) # divide into smaller chunks for parallelization 
    
    print('+++ TOTAL NUMBER OF FILES TO COMBINE: ', len(stations) )
    print('+++ SAVING THE OUTPUT FILES IN : ', out_dir  )

    
    for c in chunks:            
        print ('*** I am running CHUNK: ', chunks.index(c) , ' *** with: ' ,  len(c), ' files \n'  )
        
        stations_string = ''
        for s in c:
            stations_string = stations_string + s + ','
            
        #c = str(','.join(c)).replace('#','')
        print         ( '/opt/anaconda3/bin/python3 pre_merge_stations.py  ' + ' -o ' + out_dir + ' -s ' + stations_string+ '  -d ' + d + ' & ' )
        os.system( '/opt/anaconda3/bin/python3 pre_merge_stations.py  ' + ' -o ' + out_dir + ' -s ' + stations_string+ '  -d ' + d + ' & ')    
    
    
