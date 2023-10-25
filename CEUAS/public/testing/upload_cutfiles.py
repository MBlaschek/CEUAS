import h5py
import numpy as np
import pandas as pd
import os, glob, sys
import math
import ray

ray.init(num_cpus=10)

@ray.remote
def repack_file(inputfile):
    command = 'scp '+ inputfile +' sis@136.156.154.104:/data/public/windbe/'
    out = os.system(command)
    print(inputfile, ' done')
    return out

log = []
result_ids = []
for i in glob.glob('/users/staff/uvoggenberger/scratch/converted_v13_windbe/*.nc'):
    result_ids.append(repack_file.remote(i))
    
results = ray.get(result_ids)
print(results)
ray.shutdown()
