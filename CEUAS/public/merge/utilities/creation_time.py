import os,sys
import os.path, time
import time
from datetime import datetime
import glob
import numpy as np

### For harvested files 


db = ['era5_1','era5_2','era5_1759','era5_1761','era5_3188','bufr']
harv = '/scratch/das/federico/COP2_HARVEST_FEB2022/'
stations = open('redo_merging.txt' , 'w')
stat = []
lista = []
for d in db:

    files = [harv + '/' + d+'/'+f for f in os.listdir(harv + '/' + d) if '-1_' not in f ]
    for f  in files:
        #print("last modified: %s" % time.ctime(os.path.getmtime(f)))
        cftime = time.ctime(os.path.getctime(f))
        #print("created: %s" % cftime )
        created = datetime.strptime(cftime, "%a %b %d %H:%M:%S %Y")
        #print(date_of_created)
        m = created.month
        d = created.day
        h = created.hour
        if m==3 and d==14 or d==15:
            #print( created , ' ' , m , ' ' , d , ' ' , h )
            s = f.split('/')[-1].split('_')[0]
            print( created , ' ' , m , ' ' , d , ' ' , h , ' ' , s )
            stations.write(s + '\n')
            lista.append(s)

       
a = 0
lista = list (np.unique(lista) )

merged = '/scratch/das/federico/MERGED_25FEB2022/'
for l in lista:
    files = glob.glob( merged + '/' + l + '*')
    for f in files:
        if os.path.isfile(f):
            print(f)
            os.system('rm  ' + f )
            