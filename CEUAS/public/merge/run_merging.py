""" Utility to run the merging_cdm_netCDF.py script with parallel processes """

import os,sys
from merging_cdm_netCDF import data_directories,create_station_dics


def chunk_it(seq, num):
    """ Creates sub sets of the input file list """
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out


""" Select the number of parallel processes. Each will receive a list of stations to be merged """
processes = 6


all_stat_dics = create_station_dics(data_directories)

stations = list (all_stat_dics.keys())


chunks = chunk_it(stations,processes) # splitting into differnet list of stations

print(' Isplit the input stations list ' , len(stations) , '  into : ', processes, '  different parallel processes \n \n ' )
for c in chunks:
    cc = str(','.join(c)).replace('#','')
    print(' I am running the chunk', chunks.index(c) , ' with  ' , len(c), '  stations ' )
    os.system('/opt/anaconda3/bin/python3 merging_cdm_netCDF.py -s ' + cc + ' & ' )
