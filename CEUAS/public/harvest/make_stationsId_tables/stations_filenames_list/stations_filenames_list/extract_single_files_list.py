""" Reads the stations which do not need to be combined,
    i.e. the data is contained in one single file per dataset, 
    ready to be merged with the other dataset """


import os,sys
import pandas as pd 

ds = ['era5_3188' , 'era5_1759' , 'era5_1761' , 'era5_1' , 'bufr', 'ncar' , 'igra2']


for d in ds:
    os.system('mkdir singlefile_stations')
    out = open( 'singlefile_stations/' + d + '_singlefile_stations.txt' , 'w')
    summary_file = d + '_summary_duplicated_stations.txt'
    df = pd.read_csv( summary_file, delimiter = '\t' ) 
    for index, row in df.iterrows():
        
        primary_id = row['primary_id']
        file = row['files'][:]        
        file = file.replace("'",'').replace('[','').replace(']','')
        file = file.split(',')
        
        if len(file) > 1:
            print (row)
        else:
            file_name = file[0]
            file_name = file_name.split('/')
            name = file_name[ len(file_name) -1 ]
            out.write(primary_id + '\t' + name + '\n' )
        
print('Done ***')