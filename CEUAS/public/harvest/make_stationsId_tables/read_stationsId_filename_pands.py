""" This tool creates a list of stations_id mapping to the source files. 
      For each dataset, it reads the txt files created by the script  extract_stationid_filename.py , 
      https://github.com/MBlaschek/CEUAS/blob/develop/CEUAS/cdm/code/merging/extract_stationid_filename.py
      which is a modified version of the script used by Leopold to extract the stations_configuration
      (that additionally deals with descondary_ids ).
      
      The script above is run for each dataset, and a csv file is created.
      
      This csv file is read by this tool, and each source file is retrieved for all the different stations_id. 
      
      Again, for each dataset, a 'merged' csv file is created.
      
      """

import os,sys
import pandas as pd
import numpy as np

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)



ds = ['era5_3188' , 'era5_1759' , 'era5_1761' , 'era5_1' , 'bufr', 'ncar' , 'igra2']

#ds = ['era5_3188' , 'era5_1759' ]



all_dic = {}
for d in ds:
    print('processing: ' , d)
    dic = {}
    
    F = 'stations_filenames_list/' + d + '_stationsId_filename.txt'
    df = pd.read_csv( F, delimiter = '\t' ) 
    df = df.sort_values(by = ['primary_id'])    
    
    unique_id= []
    primary_ids = np.array(df['primary_id'].values)
    for s in primary_ids:
        if s not in unique_id:
            unique_id.append(s)
    for s in unique_id:
        dic[s] = {}
        matching = df.loc[ df['primary_id'] == s ]
        
        dic[s]['source_files'] = list(matching['file'])
        dic[s]['station'] = matching['station_name'].values[0]
        
        #print(matching)   
        
        
    print(dic)
    
    all_dic[d] = dic   # storiung int he global dictionary 

    
print('a')


''' Writing in output files ''' 
for dataset in all_dic.keys():
    out = open('stations_filenames_list/' + dataset + '_summary_duplicated_stations.txt' , 'w')
    for d in all_dic[dataset].keys():
        line = str(d)  + '\t' + str(all_dic[dataset][d]['station']) + '\t' +  str(all_dic[dataset][d]['source_files']) + '\n' 
        print( line )
        out.write(line)
    
    
    
    

"""    
F = 'era5_3188_stationsId_filename.txt'

df = pd.read_csv( F, delimiter = '\t' ) 
df = df.sort_values(by = ['station_name'])



unique_stations = []
stations = np.array(df['station_name'].values)
for s in stations:
    if s not in unique_stations:
        unique_stations.append(s)
        
    
#stations = np.unique(stations)

dic = {}
for s in unique_stations:
    line = df.loc[df['station_name'] == s ]
    print(s , ' ' , line )
    
    
    
    
#print(stations)


#print('There is a total of ' , len(stations) , '  different cities ' )
#print (df)
#print (df)
"""
