import h5py
import numpy as np
import pandas as pd
import os, glob, sys

sys.path.append(os.getcwd()+'/../cds-backend/code/')
os.environ['PYTHONPATH'] = os.getcwd()+'/../cds-backend/code/'
from harvest_convert_to_netCDF import write_dict_h5, load_cdm_tables


def check_and_fix_file(file, attr_dict):
    """
    Find missing advanced_homogenisation groups and missing variables in it.
    """

    missing_vars = []
    target_vars = ['RAOBCORE_bias_estimate', 'RASE_bias_estimate', 'RICH_bias_estimate', 'RISE_bias_estimate', 'humidity_bias_estimate', 'wind_bias_estimate']

    missing_uncerts = []
    target_uncert = ['desroziers_30'] 

    with h5py.File(file,  "r") as f:
        if not 'advanced_homogenisation' in list(f.keys()):
            ov_vars = [0] * len(f['observations_table']['date_time'][:])
            missing_vars = target_vars
        elif len(set(target_vars) - set(f['advanced_homogenisation'].keys())) > 0:
            ov_vars = [0] * len(f['observations_table']['date_time'][:])
            missing_vars = list(set(target_vars) - set(f['advanced_homogenisation'].keys()))

        if not 'advanced_uncertainty' in list(f.keys()):
            ov_uncerts = [np.nan] * len(f['observations_table']['date_time'][:])
            missing_uncerts = target_uncert
        elif len(set(target_uncert) - set(f['advanced_uncertainty'].keys())) > 0:
            ov_uncerts = [np.nan] * len(f['observations_table']['date_time'][:])
            missing_uncerts = list(set(target_uncert) - set(f['advanced_uncertainty'].keys()))

    if len(missing_uncerts) == 0 and len(missing_vars) == 0:
        return 0

    for missing_var in missing_vars:
        alldict = pd.DataFrame({missing_var:ov_vars})
        write_dict_h5(file, alldict, 'advanced_homogenisation', {missing_var: { 'compression': 'gzip' } }, [missing_var], mode='a', attrs=attr_dict)

    for missing_uncert in missing_uncerts:
        alldict = pd.DataFrame({missing_uncert:ov_uncerts})
        write_dict_h5(file, alldict, 'advanced_uncertainty', {missing_uncert: { 'compression': 'gzip' } }, [missing_uncert], mode='a', attrs=attr_dict)

    return 1



# def find_bad_file(file):
#     with h5py.File(file,  "r") as f:
#         if not 'advanced_uncertainty' in list(f.keys()):
#             print(file)

# path = '/mnt/users/scratch/leo/scratch/converted_v25/long/'
# all_files = glob.glob(path + "*.nc")

# for i in all_files:
#     find_bad_file(i)



path = '/mnt/users/scratch/leo/scratch/converted_v25/long/'
path = '/mnt/users/scratch/uvoggenberger/test_add_group/'
all_files = glob.glob(path + "*.nc")
attr_file = glob.glob(path + "*11035*.nc")



attr_dict = {}

target_vars = ['RAOBCORE_bias_estimate', 'RASE_bias_estimate', 'RICH_bias_estimate', 'RISE_bias_estimate', 'humidity_bias_estimate', 'wind_bias_estimate']
target_uncerts = ['desroziers_30']
with h5py.File(attr_file[0],  "r") as f:
    for tv in target_vars:
        attr_dict[tv] = {}
        for atr in f['advanced_homogenisation'][tv].attrs:
            if atr != 'DIMENSION_LIST':
                attr_dict[tv][atr] = f['advanced_homogenisation'][tv].attrs[atr]
    for uc in target_uncerts:
        attr_dict[uc] = {}
        for atr in f['advanced_uncertainty'][uc].attrs:
            if atr != 'DIMENSION_LIST':
                attr_dict[uc][atr] = f['advanced_uncertainty'][uc].attrs[atr]


for file in all_files:
    result = check_and_fix_file(file, attr_dict)
    if result == 1:
        print(file, ": groups have been added!")