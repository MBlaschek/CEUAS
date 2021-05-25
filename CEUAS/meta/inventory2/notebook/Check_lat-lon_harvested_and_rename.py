""" Check if harvested file contain multiple lat and lon.
      Rename the file accordingly using the identifier -20999- """
import subprocess
import pandas as pd
import os
import glob
import pandas as pd
import matplotlib.pylab as plt
import os,sys
import json
import xarray as xr
import numpy as np
import matplotlib.animation as animation
import h5py
from multiprocessing import Pool
from functools import partial

import numpy as np 

def run_fn(fn):

    if '-20999-' in fn: # nothing to do, already checked 
        return 0
    
    try:
        data = h5py.File(fn,'r')
        lats = np.unique(data['observations_table']['latitude'][:])
        lons = np.unique(data['observations_table']['longitude'][:])
        
        if  len(lats) ==1 and len(lons) ==1: # one single pair of cordinates 
            pass
        else:

            lat1, lat2 = min(lats), max(lats)
            lon1, lon2 = min(lons), max(lons)

            if (abs(lat1-lat2) < 0.5 and abs(lon1-lon2) < 0.5):
                pass

            else:
                fn_n = fn.replace('-20000-','-20999-').replace('-20001-','-20999-').replace('-20300-','-20999-').replace('-20400-','-20999-')
                fn_n = fn_n.replace('-20500-','-20999-').replace('-20600-','-20999-')
                
                os.system('mv ' + fn + '   ' + fn_n )
                print('moved ' + fn + '   ' + fn_n  )
                
    except:
        print("*** Fail: ", fn )
        pass
    
    return 0 


multi = True
for d in ['era5_1759','era5_3188','igra2', 'ncar', 'bufr']:

     harvested_dir = '/raid60/scratch/federico/MAY2021_HARVEST_secondary/' + d
     flist = os.listdir( harvested_dir )
     flist = [harvested_dir +'/'+f for f in flist if  '.nc' in f and '00000' not in f and '9999' not in f ]
     print(len(flist) , ' *** Files ')
     tot = len(flist)
 
 
     if multi:
          
          p=Pool(10)
          func=partial(run_fn)
          all_res=list(p.map(func,flist))
 
     else:
          all_res = []
          for f in flist:
               r = run_fn(f)
               all_res.append(r)
 

