import glob
import pandas as pd 
import numpy as np
import h5py
import hdf5plugin
import ray
import multiprocessing
from functools import partial
import pickle
from collections import Counter


# @ray.remote
def get_source_count(file):
    with h5py.File(file) as fl:
        oid = fl['observations_table']['observation_id'][:]
        sid = fl['observations_table']['source_id'][:]

    df = pd.DataFrame.from_dict({'oid':np.array([''.join(row.astype(str)) for row in oid]), 'sid':np.array([''.join(row.astype(str)) for row in sid])})
    df = df.drop_duplicates('oid')
    
    return dict(df.value_counts('sid'))

pool = multiprocessing.Pool(processes=40)
year_counts = []

for year in range(2015,2025,1):
    print(year)
    files = glob.glob('/mnt/users/scratch/leo/scratch/converted_v27/'+str(year)+'/*.nc')
    
    result_list = list(pool.map(get_source_count, files[:]))

    total_counts = Counter()
    for d in result_list:
        total_counts.update(d)

    year_counts.append(total_counts)

    print(total_counts)
    


with open('20250204_stored_data_station_sources_2015_to_2024.p', 'wb') as file:
		pickle.dump(year_counts, file)