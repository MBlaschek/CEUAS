""" Utility to run the merging_cdm_netCDF.py script with parallel processes """

import os,sys

#from merging_cdm_netCDF_newfixes import data_directories, create_station_dics
from merging_cdm_netCDF import data_directories, out_dir


""" Select the minimum and maximum size of the files to be processed """
max_size = 100*10**9     # 5*10**6 = 5 MB
min_size = 1*10**4



def create_station_list(datadir , kind='mobile', processed_dir = '' ):
    """ Finds all the stations if, extracting from the files in the data directories """
    all_stat = []
    
    if processed_dir:
        try:
            processed_stat = [ f.split('_')[0] for f in os.listdir(processed_dir) ]
        except:
            print('+++ Could not find existing processed directory ')
            processed_stat = []
    for d,i in datadir.items():
        if kind == 'mobile':
            if d not in ['era5_1_mobile', 'era5_2_mobile']:
                continue
            files = [ f.split('_era5_')[0] for f in os.listdir(i) if '.nc' in f and '20999' in f ]
            all_stat.extend(files)
        
        else:
            files = [ f for f in os.listdir(i) if '.nc' in f and '20999' not in f ]

        
            for F in files:
                stat_id = F.split('_'+d)[0]                
                if len(stat_id) > 5:
                    if processed_dir:
                        if stat_id in processed_stat:
                            #print("Already processed: " , stat_id )
                            continue 
                        else:
                            if stat_id not in all_stat:                   
                                all_stat.append(stat_id)
    
    return all_stat 


def clean_list(stat_list, data_directories, min_size = '', max_size='', kind=''):
    """ Cleans the file list of all the file station according the to total dimension of the files in the list """
    stations_cleaned = []
    if kind == 'mobile':
        return stat_list
    
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
    import random
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    random.shuffle(seq)
    
    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out




kind = ''
kind = 'mobile'


""" Create station list to process, exclude processed stations """
lista = create_station_list(data_directories , kind = kind, processed_dir = out_dir )

""" Clean list wrt size """

cleaned_list = clean_list(lista, data_directories, min_size = min_size, max_size = max_size, kind=kind  )

print('To be processed: ', len(cleaned_list), ' :::: files')


""" Select the number of parallel processes. Each will receive a list of stations to be merged """
processes = 15

chunks = chunk_it(lista, processes) # splitting into differnet list of stations
#print('I split the original input stations list ' , len(lista) , '  into : ', processes, '  different parallel processes \n \n ' )


for c in chunks:
    cc = str(','.join(c)).replace('#','')
    print(' I am running the chunk', chunks.index(c) , ' with  ' , len(c), '  stations ' )
    print(c[:10])
    os.system('python3  merging_cdm_netCDF.py  -s ' + cc + ' -m ' + str(max_size) + ' -l ' + str(min_size) + '  -k ' + kind + '  &' )


