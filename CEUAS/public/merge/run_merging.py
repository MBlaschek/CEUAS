""" Utility to run the merging_cdm_netCDF.py script with parallel processes """

import os,sys

#from merging_cdm_netCDF_newfixes import data_directories, create_station_dics
from merging_all_last import data_directories 

def create_station_list(datadir):
    """ Finds all the stations in the data directories """
    all_stat = []
    for d,i in datadir.items():
        files = [ f for f in os.listdir(i) if '.nc' in f ]
        for F in files:
            stat_id = F.split('_'+d)[0]                
            if len(stat_id) > 5:
                if stat_id not in all_stat:                   
                    all_stat.append(stat_id)
                    
    return all_stat 
    
def chunk_it(seq, num):
    """ Creates sub sets of the input file list """
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out




lista = create_station_list(data_directories)


""" Select the number of parallel processes. Each will receive a list of stations to be merged """
processes = 3

chunks = chunk_it(lista, processes) # splitting into differnet list of stations

max_size = 500*10**6 # 5*10**6 = 5 MB


print('I split the input stations list ' , len(lista) , '  into : ', processes, '  different parallel processes \n \n ' )

for c in chunks:
    cc = str(','.join(c)).replace('#','')
    print(' I am running the chunk', chunks.index(c) , ' with  ' , len(c), '  stations ' )
    os.system('/opt/anaconda3/bin/python3 merging_all_last.py -s ' + cc + ' -m ' + str(max_size) + ' &' )


