import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

import numpy as np
import glob

import ray

@ray.remote
def concat_files(sfd, out_dir):

    files = glob.glob(sfd + '/*')
    concat_list = []

    # collect files to concat:
    for fi in files:
        concat_list.append(pd.read_csv(fi, delimiter='\t'))
    
    # concat the data
    df = pd.concat(concat_list)
    
    # write to file:
    df.to_csv(out_dir + sfd.split('/')[-1] + '.csv', sep='\t')


if __name__ == '__main__':

    out_dir = '/users/staff/uvoggenberger/scratch/bufr_cnr/concated/'
    split_files_dirs = glob.glob('/users/staff/uvoggenberger/scratch/bufr_cnr/split_new_2/*')

    ray.init(num_cpus=40)

    result_ids = []
    for i in split_files_dirs[:]:
        result_ids.append(concat_files.remote(i, out_dir))
    
    results = ray.get(result_ids)
    ray.shutdown()

    # concat_files(split_files_dirs[0], out_dir)

