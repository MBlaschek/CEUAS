import os,sys
import pandas as pd 



def extract_latitude_mismatch_1759(v='era5_1759'):
    """ Creates list of files where a missing minus sign in latitude is very likely, from WBAN inventory """
    
    path = 'inventories/' +v 
    files = [path +'/'+ f for f in os.listdir( path ) if 'reduced' in f ]
    
    for file in files:
        df = pd.read_csv(file, sep = '\t', header = 0)
        print(0)
    
    
dummy = extract_latitude_mismatch_1759(v='era5_1759')

