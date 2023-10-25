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

# ray.init(num_cpus=3)
# ray.shutdown()

# @ray.remote
def repack_file(inputfile):
    try:
        stat = inputfile.split('/')[-1] 
        targetfile = "/users/staff/uvoggenberger/scratch/converted_v13_test_repacked/"+stat
        repackfile = "/users/staff/uvoggenberger/scratch/converted_v13_test_repacked_new/"+stat

        mode='a'
        group = 'header_table'
        i = 'product_version'
            
        with eua.CDMDataset(inputfile) as file:
            ov_vars = file[group][i][:]
            alldict = pd.DataFrame({i:ov_vars})
            
        with h5py.File(targetfile,  "a") as f:
            del f[group][i]
            
        write_dict_h5(targetfile, alldict, group, {i: { 'compression': 'gzip' } }, [i]) 
        command = 'h5repack -v GZIP=5 '+ targetfile + ' ' + repackfile
        log = os.system(command)
        
        command = 'rm '+ targetfile
        log = os.system(command)
        
    except:
        print(inputfile)
        return inputfile
    return 0

log = []
for i in glob.glob('/users/staff/uvoggenberger/scratch/converted_v13_test/*.nc'):
    log.append(repack_file(i))
