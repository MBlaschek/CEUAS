import numpy as np
import pandas as pd
import os, glob, sys
import ray

ray.init(num_cpus=10)

@ray.remote
def repack_file(inputfile):
    command = 'scp '+ inputfile +' sis@136.156.154.104:/data/private/converted_v18/'
    out = os.system(command)
    
    file_path = '/users/staff/uvoggenberger/scratch/upload/' + inputfile.split('/0-')[-1].split('_C')[0] + '.txt'
    with open(file_path, "w") as text_file:
        text_file.write(str(out))
        text_file.write('\n')
        text_file.write('done')
        text_file.write('\n')

    print(inputfile, ' done', out)
    return out

log = []
result_ids = []
for i in glob.glob('/mnt/users/scratch/leo/scratch/converted_v18/long/*.nc')[10:20]:
    print(i)
    file_path = '/users/staff/uvoggenberger/scratch/upload/' + i.split('/0-')[-1].split('_C')[0] + '.txt'
    print(file_path)
    if len(glob.glob(file_path)) < 1:
        result_ids.append(repack_file.remote(i))
    # result_ids.append(repack_file.remote(i))
    
results = ray.get(result_ids)
print(results)
ray.shutdown()
