""" Utility to run the merging_cdm_netCDF.py script with parallel processes """

import os,sys

#from merging_cdm_netCDF_newfixes import data_directories, create_station_dics
from merging_all_last import data_directories 


""" Select the minimum and maximum size of the files to be processed """
max_size = 11*10**8     # 5*10**6 = 5 MB
min_size = 1*10**3



def create_station_list(datadir , processed_dir = '' ):
    """ Finds all the stations if, extracting from the files in the data directories """
    all_stat = []
    
    if processed_dir:
        processed_stat = [ f.split('_')[0] for f in os.listdir(processed_dir) ]
        
    for d,i in datadir.items():
        files = [ f for f in os.listdir(i) if '.nc' in f ]
        
        for F in files:
            stat_id = F.split('_'+d)[0]                
            if len(stat_id) > 5:
                if processed_dir:
                    if stat_id in processed_stat:
                        #print("Already processed: " , stat_id )
                        continue 
                   
                if stat_id not in all_stat:                   
                    all_stat.append(stat_id)
                    
    return all_stat 


def clean_list(stat_list, data_directories, min_size = '', max_size=''):
    """ Cleans the file list of all the file station according the to total dimension of the files in the list """
    stations_cleaned = []
       
    for s in stat_list:
           
        total_dim = []        
        for d,i in data_directories.items():
            files = os.listdir(i)
            for f in files:
                Id = f.split('_'+d)[0]
                if Id == s:           
                    total_dim. append( os.path.getsize (i + '/' + f) )
                    
        size = sum(total_dim)                 
        if size < min_size or size > max_size:
            #print('Skipping: size does not match ' , size/1000000000 )
            continue 
        stations_cleaned.append(s)
        #print('Keeping station: ' , s  , ' since the size is: ',  size/1000000000  )
        
    return stations_cleaned      
    
    
    
    
def chunk_it(seq, num):
    """ Creates sub sets of the input file list """
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out





lista = create_station_list(data_directories , processed_dir = '/raid60/scratch/federico/JUNE_TEST_MERGING_ALL' )
cleaned_list = clean_list(lista, data_directories, min_size = min_size, max_size = max_size  )
print('To be processed: ', len(cleaned_list), ' :::: files')

""" Select the number of parallel processes. Each will receive a list of stations to be merged """
processes = 6

chunks = chunk_it(lista, processes) # splitting into differnet list of stations



#print('I split the original input stations list ' , len(lista) , '  into : ', processes, '  different parallel processes \n \n ' )


for c in chunks:
    cc = str(','.join(c)).replace('#','')
    #print(' I am running the chunk', chunks.index(c) , ' with  ' , len(c), '  stations ' )
    os.system('/opt/anaconda3/bin/python3 merging_all_last.py -s ' + cc + ' -m ' + str(max_size) + ' -l ' + str(min_size) + '  &' )


