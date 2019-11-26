import os,sys
import pandas as pd
import numpy as np

""" Extracts all the stations from all the station_configuration files,
mapping them to the available source files """


""" For the igra2 dataset,
the primary station id is contained in the name of the source file.
This means that by reading the primary_id from th estation configuration file,
I loop over the possible file name in the igra2 databse directory, and I identify
which file name contains the desired primary_id.
This is done since the name cannot be recontructed by the primary_id number only
(they have other spceific alpha-numerical code)

One file is dubious, i.e. the same primary id maps onto 
two different file names:

94750 => ['ASM00094750-data.txt', 'USW00094750-data.txt']
"""



def find_igra2_file(n, files_list):
    """ Given a primary id in input, loops over the ifra2 files and finds a matching filename """        
    try:
        n = n.split('-')[3]
    except:
        return 0
    name_igra = []            
    for f in files_list:
        if n ==  f.split('-')[0][-5:]:
            name_igra.append(f)
            continue 
    #print('found ?' , len(name_igra) , ' ',  name_igra )        
    if len(name_igra) == 2:
        print ('Dubious file' , name_igra  )
        #dubious[n] = []
        #dubious[n].append(name_igra)
    elif len(name_igra) == 1:
        return name_igra[0]
    elif len(name_igra) == 0:
        print('igra2 not found ! +++ ')
        return np.nan 


def find_era5_1_file(n , files_list):
    """ Given a primary id in input, loops over the ifra2 files and finds a matching filename """        
    try:
        n = n.split('-')[3]
    except:
        return 0
    name = 'era5.conv._' + str(n)
    name_era5 = []            
    for f in files_list:
        if n ==  name:
            name_era5.append(f)
            continue 
    #print('found ?' , len(name_igra) , ' ',  name_igra )        
    if len(name_era5) == 2:
        print ('Dubious file' , name_era5  )
        #dubious[n] = []
        #dubious[n].append(name_igra)
    elif len(name_era5) == 1:
        print('found !!! ' ,  name , ' ****************************' )
        return name_era5[0]
    elif len(name_era5) == 0:
        print('era5_1 file not found ! ', name ,'  +++ ' )
        return np.nan 
    
    
        
        
        

ds = ['bufr'      , 
      'igra2'     ,
      'ncar'      ,
      'era5_1759' ,
      'era5_1761' ,
      'era5_1'    ,
      'era5_3188'   ]



# station_configurations list
stat_conf = [ d + '_station_configuration.dat' for d in ds ]



df = {}
""" Extracting primary ids from all the datasets and making a unique list """
all_primary = []
for d in ds:
    df[d] = pd.read_csv( 'station_configuration_'  +  d  +  '.dat' , sep='\t' )
    p = list( df[d]['primary_id'].values )
    all_primary = all_primary + p
    print('done station_conf +++  ' , d )    
uniques = np.unique (np.asarray (all_primary) ) 
    
    
""" Reading the filename lists """ 
df_names ={}
for d in ['bufr' ,  'ncar' ,  'era5_1759' , 'era5_1761' ,  'era5_3188' ]:
    df_names[d] = pd.read_csv(  d  +  '_stationIds_filename.dat' , sep='\t' , names = ['primary_id', 'secondary_id', 'file_name'] ) 
    print('done filename *** ' , d )    


""" Looping over the unique primary ids, find the filename , storing a dictionary """ 

map_dic = { 'u': [] , 'era5_1': [] , 'era5_3188':[] , 'era5_1759':[], 'era5_1761':[], 'igra2':[] , 'ncar': [] , 'bufr': [] }

for u in uniques:
    print( 'Processing   : ', list(uniques).index(u) )  
    map_dic['u'].append(u)    
    F = np.nan 
    for d, df in df_names.items(): # all ds except era5_1 and igra2
        primary = df['primary_id']
        file         = df['file_name']
        
        for p, f in zip (primary, file):
            if u ==p:
                F = f 
                continue
        map_dic[d].append(F)
        
    files_list = os.listdir('/raid60/scratch/federico/databases/IGRAv2') # getting the list of igra2 files             
    ig2 = find_igra2_file(u, files_list)
    map_dic['igra2'].append(ig2)
   
   
   
    files_list = os.listdir('/raid60/scratch/leo/scratch/era5/odbs/1') # getting the list of igra2 files      
    files_list = [ f for f in files_list if  '.conv._' in f  and '.gz' not in f and '.nc' not in f ]   
    era5_1 = find_era5_1_file(u, files_list)
    map_dic['era5_1'].append(era5_1)
        
        

""" Saving the dictionary as a DF """
df = pd.DataFrame.from_dict( map_dic  )


print('done done')