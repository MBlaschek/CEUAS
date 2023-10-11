import glob
import os
import pandas
import numpy
import h5py

log = []
files = glob.glob('/mnt/users/scratch/leo/scratch/converted_v13/long/*.nc')
for fn in files:
    with h5py.File(fn,'r') as f:
        try:
            a = f['header_table']['product_version'][0]
        except Exception as e:
            log.append(fn)
            print(e)

with open('h5_log.txt', 'w') as f:
    for j in log:
        f.write(j)
