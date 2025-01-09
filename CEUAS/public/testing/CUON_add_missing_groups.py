import h5py
import numpy as np
import pandas as pd
import os, glob, sys

sys.path.append(os.getcwd()+'/../cds-backend/code/')
os.environ['PYTHONPATH'] = os.getcwd()+'/../cds-backend/code/'
from harvest_convert_to_netCDF import write_dict_h5


def check_and_fix_file(file):
    """
    Find missing advanced_homogenisation groups and missing variables in it.
    """
    missing_vars = []
    target_vars = ['RAOBCORE_bias_estimate', 'RASE_bias_estimate', 'RICH_bias_estimate', 'RISE_bias_estimate', 'humidity_bias_estimate', 'wind_bias_estimate']
    with h5py.File(file,  "r") as f:
        if not 'advanced_homogenisation' in list(f.keys()):
            ov_vars = [np.nan] * len(f['observations_table']['date_time'][:])
            missing_vars = target_vars
        elif len(set(target_vars) - set(f['advanced_homogenisation'].keys())) > 0:
            ov_vars = [np.nan] * len(f['observations_table']['date_time'][:])
            missing_vars = list(set(target_vars) - set(f['advanced_homogenisation'].keys()))
        else:
            return 0

    for missing_var in missing_vars:
        alldict = pd.DataFrame({missing_var:ov_vars})
        write_dict_h5(file, alldict, 'advanced_homogenisation', {missing_var: { 'compression': 'gzip' } }, [missing_var])
    return 1


path = '/mnt/users/scratch/leo/scratch/converted_v25/long/'
path = '/mnt/users/scratch/uvoggenberger/test_add_group/'
all_files = glob.glob(path + "*.nc")

for file in all_files:
    result = check_and_fix_file(file)
    if result == 1:
        print(file, ": advanced_homogenisation has been added!")