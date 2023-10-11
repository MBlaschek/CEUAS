import numpy
import numpy as np
import pandas as pd
import sys, glob
import urllib3
import h5py
import cdsapi, zipfile, os, time
sys.path.append(os.getcwd()+'/../cds-backend/code/')
import cds_eua3 as eua
import warnings
import pickle
warnings.filterwarnings('ignore')

decades = [['19900101','19991231'],['20000101','20091231'],['20100101','20191231']]
decnames = ['1990s','2000s','2010s'] 
stats = glob.glob('/raid60/scratch/uli/converted_v5/*nc')
err = []
for i in range(len(decades)):
    print(i)
    print('\n')
    out = pd.DataFrame(columns=['lat', 'lon', 'obs', 'stat'])
    for o in stats:
        sid = o.split('/')[-1][:-19]
        print(sid)
        rmfiles = glob.glob('dest*.nc')
        for j in rmfiles:
            os.remove(j) 
        try:
            data = eua.vm_request_wrapper(vm_url='http://localhost:8000',overwrite=True, request={'variable':['temperature'], 'statid': sid, 'date': decades[i], 
                                                                                                  'pressure_level':50000})
            data = data.to_dataframe()
            da = data.dropna(axis=0, how='any')
            out = out.append({'lat': da.iloc[0].lat, 'lon': da.iloc[0].lon, 'obs': len(da), 'stat': da.iloc[0].station_id}, ignore_index=True)
        except:
            print('no data for this station')
        
    with open('/raid60/scratch/uli/cuon_values_500hpa_'+decnames[i]+'.p', 'wb') as handle:
        pickle.dump(out, handle, protocol=pickle.HIGHEST_PROTOCOL)
    print(err)
