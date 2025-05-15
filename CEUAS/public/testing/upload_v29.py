import numpy as np
import pandas as pd
import os, glob, sys
import ray

# ray.init(num_cpu=5)

# @ray.remote
def repack_file(inputfile):
    command = 'scp '+ inputfile +' obs@136.156.140.94:/data/public/converted_v29/'
    out = os.system(command)
    
    file_path = '/mnt/users/scratch/uvoggenberger/upload/' + inputfile.split('/0-')[-1].split('_C')[0] + '.txt'
    with open(file_path, "w") as text_file:
        text_file.write(str(out))
        text_file.write('\n')
        text_file.write('done')
        text_file.write('\n')

    print(inputfile, ' done', out)
    return out


log = []
result_ids = []
for i in glob.glob('/mnt/users/scratch/leo/scratch/converted_v29/long/*.nc')[:]:
    print(i)
    file_path = '/mnt/users/scratch/uvoggenberger/upload/' + i.split('/0-')[-1].split('_C')[0] + '.txt'
    print(file_path)
    already_done = list(pd.read_csv('/mnt/users/scratch/uvoggenberger/upload/uploaded.txt', names=['files']).files)
    if not i.split('/')[-1] in already_done:
        if len(glob.glob(file_path)) < 1:
            repack_file(i)
        # result_ids.append(repack_file.remote(i))
    # result_ids.append(repack_file.remote(i))
    
# results = ray.get(result_ids)
# print(results)
# ray.shutdown()
