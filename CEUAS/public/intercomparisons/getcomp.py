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
import shutil
import pickle
warnings.filterwarnings('ignore')

def request(rqdict, source, remove_file=True):
    t0 = time.time()

    c = cdsapi.Client()
    r = c.retrieve(
        source,rqdict)
    if True:
        r.download(target='download.zip')
        assert os.stat('download.zip').st_size == r.content_length, "Downloaded file is incomplete"
    z = zipfile.ZipFile('download.zip')
    z.extractall(path='./download/')
    z.close()
    print('Request took: ' + str(time.time() - t0) + ' seconds')
    
    files = glob.glob('./download/*.nc')
    
    if files[0].split('/')[-1].startswith('IGRA'):
        ds = xarray.open_dataset(files[0])            
        data = ds.to_dataframe()
        for i in files[1:]:
            ds = xarray.open_dataset(i)            
            data = data.append(ds.to_dataframe())

    else:
        data=eua.CDMDataset(files[0]).to_dataframe()
        for i in files[1:]:
            da = eua.CDMDataset(i).to_dataframe()
            data = data.append(da)
            
    os.remove('download.zip')
    if remove_file:
        try:
           shutil.rmtree('./download/')
        except:
           print('Error while deleting directory')

    return data

stations = pickle.load( open( "stations.p", "rb" ))
for i in stations.station_name:
    if not '/raid60/scratch/uli/IGRA_Data/COMP/COMP_'+str(i)+'_100.p' in glob.glob('/raid60/scratch/uli/IGRA_Data/COMP/*'):
        try:
            data = request({'variable':['temperature'],
                            'statid': i[-5:],
                            'period': '1979-01-01/2019-12-31',
                            'optional':'bias_estimate',
                            'pressure_level': 100,
                           }, 'insitu-comprehensive-upper-air-observation-network',
            )
        except:
            pass
        pickle.dump( data, open( '/raid60/scratch/uli/IGRA_Data/COMP/COMP_' + str(i) + '_100.p', "wb" ))
    else:
        print('-')