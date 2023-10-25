""" Extract list of mobile stations from ERA2 file list """

from tqdm import  tqdm
import glob
import pandas as pd 
lista = open('era5_2_files_list.txt','r').readlines()

#stats = [ l.split('conv._')[1].strip('\n') for l in lista ]

'''
out = open('era5_2_stationtype.txt', 'w')

for s in tqdm(lista):
    file = s.replace('\n', '.gz')
    
    try:
        
        df = pd.read_csv(file, sep = '\t', nrows=5)
        rt = df.reportype.values[0]
    
        out.write(s.strip('\n') + ' \t' + str(rt) + '\n' )
    
    except:
        a = open('missing_era5_gzip.txt', 'a+')
        a.write(s)
        a.close()
        
    a = 0
out.close()
'''

df = pd.read_csv('era5_2_stationtype.txt', sep = '\t', names = ['file', 'type'])

land = open('era5_2_files_list.txt' , 'w')
mobile = open('era5_2_mobile_files_list.txt' , 'w')


for index, row in df.iterrows():
    file = row['file']
    typ = row['type']
    
    if typ in [16013 , 16022]:
        #print(file, typ)
        land.write(file + '\n')
    else:
        mobile.write(file + '\n')
        
land.close()
mobile.close()

### 16013, 16022 LAND 

### all other:  MOBILE 

print(0)