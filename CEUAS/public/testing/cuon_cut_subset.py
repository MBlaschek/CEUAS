import h5py
import numpy as np
import pandas as pd
import os, glob, sys
import math
import ray

sys.path.append(os.getcwd()+'/../cds-backend/code/')
os.environ['PYTHONPATH'] = os.getcwd()+'/../cds-backend/code/'
import cds_eua4 as eua
from harvest_convert_to_netCDF import write_dict_h5

ray.init(num_cpus=10)

@ray.remote
def repack_file(inputfile):
    try:
        stat = inputfile.split('/')[-1] 
        targetfile = "/users/staff/uvoggenberger/scratch/converted_v13_windbe/"+stat

        mode='w'
        group = 'advanced_homogenisation'
        i = 'wind_bias_estimate'
            
        with eua.CDMDataset(inputfile) as file:
            ov_vars = file[group][i][:]
            alldict = pd.DataFrame({i:ov_vars})
            
        write_dict_h5(targetfile, alldict, group, {i: { 'compression': 'gzip' } }, [i]) 
        
    except:
        print(inputfile)
        return inputfile
    return 0

log = []
result_ids = []
for i in glob.glob('/mnt/users/scratch/leo/scratch/converted_v13/long/*.nc'):
    # log.append(repack_file(i))
    result_ids.append(repack_file.remote(i))
    
results = ray.get(result_ids)
ray.shutdown()
