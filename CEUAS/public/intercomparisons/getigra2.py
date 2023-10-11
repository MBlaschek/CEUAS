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
import xarray
import pickle
warnings.filterwarnings('ignore')

def request(rqdict, source, sdir):
    t0 = time.time()

    c = cdsapi.Client()
    r = c.retrieve(
        source,rqdict)
    print('Request took: ' + str(time.time() - t0) + ' seconds')
    if True:
        r.download(target='download.zip')
        assert os.stat('download.zip').st_size == r.content_length, "Downloaded file is incomplete"
    z = zipfile.ZipFile('download.zip')
    z.extractall(path=sdir)
    z.close()
    return

# ---

igra = 'insitu-observations-igra-baseline-network'

# ---

files = glob.glob('tempdownload/*/*')
for j in files:
    print(j)
    ds = xr.open_dataset(i)
    data = ds.to_dataframe()
    data = data[data.air_pressure == 50000]
    data = data.dropna(axis=0, how='any')
    pickle.dump( out, open( '/raid60/scratch/uli/IGRA_Data/IGRA/IGRA_' + str(j) + '_500.p', "wb" ))
    print('done')


