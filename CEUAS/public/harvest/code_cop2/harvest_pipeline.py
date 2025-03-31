import os
import sys
import subprocess
import pandas as pd
import numpy as np
import datetime
import glob


'''
change: set year in param file to the one from date check 
add: era5_1 download and processing
'''

# Create the directory to store the data for the current date:
global ceuas_dir
ceuas_dir = '/srvfs/home/uvoggenberger/CEUAS/CEUAS/'

global base_dir
base_dir = '/mnt/users/scratch/uvoggenberger/CUON_HARVEST'

global python_interpreter
python_interpreter = '/srvfs/home/uvoggenberger/micromamba/envs/uv12/bin/python'

datetime_now = datetime.datetime.now()
date_now = datetime_now.strftime('%Y%m')

global working_dir
working_dir = base_dir + '_' + date_now + '/'
os.system('mkdir -p ' + working_dir)

global table_dir
table_dir = f'{ceuas_dir}/meta/inventory_comparison_2/data/tables/'


def download_data_igra2(rm_zip=False):
    # Download data from NOAA
    igra_dir = working_dir + '/data/igra_data'
    if len(glob.glob(igra_dir + '/*.txt')) > 0:
        os.system('mkdir -p ' + igra_dir)
        os.system(f"wget -r -np -nH --cut-dirs=6 -A '*.zip' -P {igra_dir} https://www.ncei.noaa.gov/data/integrated-global-radiosonde-archive/access/data-y2d/")
        # Unzip the files
        os.system(f'unzip "{igra_dir}/*.zip" -d "{igra_dir}/"')
        # Remove the zip files
        if rm_zip:
            os.system(f'rm {igra_dir}/*.zip')
        # Get the list of files
        files = glob.glob(igra_dir + '/*.txt')
    else:
        print('Data already downloaded')
        files = glob.glob(igra_dir + '/*.txt')
    return files

def download_data_era5(rm_zip=False):
    
    # !/bin/bash
    # # Load the Teleport module
    # module load teleport
    # # Initiate Teleport login
    # tsh login
    # # Connect to the jump server
    # ssh -J lh4@jump.ecmwf.int lh4@hpc-login
    # doesn't work at all -> running multiple times for the same month??? 

    era5_dir = working_dir + '/data/era5_data'
    os.system('mkdir -p ' + era5_dir)

    # copy and edit: "/mnt/users/scratch/uvoggenberger/CUON_HARVEST_202503/data/splite5_1.ksh" also chmod +x on it, to make it executeable
    # run splite5_1 on the downloaded data and then run:
    # copy and edit: "/mnt/users/scratch/uvoggenberger/CUON_HARVEST_202503/code/odbgz_new.py" set month/year/... -> change paths to variable input!
    # copy and edit: "/mnt/users/scratch/uvoggenberger/CUON_HARVEST_202503/code/odbgz_new_mobile.py" set month/year/... -> change paths to variable input! 

    ll = glob.glob('./*')
    ls = [l for l in ll if not ".gz" in l]
    ls = [l for l in ls if len(l.split('.')) > 4]
    for i in ls:
        if len(glob.glob(i + '.gz')) < 1:
            os.system(f'cp {i} ../era5_1_mobile_data/')


def copy_tables_to_harvest():
    # Copy the tables to the harvest directory
    os.system(f'mkdir -p {working_dir}/data/tables/')
    os.system(f'cp -r {table_dir}* {working_dir}/data/tables/')
    os.system(f'cp {ceuas_dir}/meta/inventory_comparison_2/data/* {working_dir}/data/') 
    os.system(f"wget -O {working_dir}/data/tables/igra2-station-list.txt https://www.ncei.noaa.gov/data/integrated-global-radiosonde-archive/doc/igra2-station-list.txt")
    os.system(f'mkdir -p {working_dir}/code/file_list/')
    os.system(f'mkdir -p {working_dir}/code/inventories/')
    os.system(f'cp {working_dir}/data/tables/igra2-station-list.txt {working_dir}/code/file_list/')
    
def create_inventory(data_set):
    # Create the inventory
    input_dir = f'{working_dir}/data/{data_set}_data'
    analyze_inventory_functions = f'{ceuas_dir}/meta/inventory_comparison_2/code/analyze_inventory_functions.py'
    print(f"python {analyze_inventory_functions} -d {data_set} -w {working_dir} -i {input_dir}")
    os.system(f"module load odc; python {analyze_inventory_functions} -d {data_set} -w {working_dir} -i {input_dir}") # somehow selects the old dir -> fix this

def make_station_configuration(data_set):
    # Create the station configuration
    os.system(f'mkdir -p {working_dir}/code/station_configuration/')
    make_station_configuration_functions = f'{ceuas_dir}/meta/inventory_comparison_2/code/make_station_configuration.py'
    os.system(f"python {make_station_configuration_functions} -d {data_set} -w {working_dir}")
    o_file = glob.glob(f"{working_dir}/code/station_configuration/{data_set}_orphan*.csv")[0]
    os.system(f"cp {o_file} {o_file.replace('orphan', 'mobile')}")

def run_harvester(data_set, stat_kind='regular'):
    # Set up the parameter file:
    original_file_path = f'{ceuas_dir}/public/harvest/code_cop2/harvester_yearsplit_parameters.py'
    new_file_path = f'{working_dir}/code/modded_harvester_yearsplit_parameters.py'
    os.system(f'mkdir -p {working_dir}/harvest/')

    # for stat_kind in ['regular']: # , 'mobile', 'orphan']:
    #     if data_set == 'igra2' and stat_kind == 'mobile':
    #         continue

    # Read the content of the original file
    with open(original_file_path, 'r') as file:
        content = file.readlines()

    # Modify the content as needed
    modified_content = []
    for line in content:
        if "'igra2': " in line or "'igra2_mobile': " in line:
            line = line.split(': ')[0] + f': "{working_dir}/data/igra_data",\n'
        elif "'era5_1': " in line:
            line = line.split(': ')[0] + f': "{working_dir}/data/era5_1_data",\n'
        elif "'era5_1_mobile': " in line:
            line = line.split(': ')[0] + f': "{working_dir}/data/era5_1_mobile_data",\n'
        elif "datasets = " in line:
            line = line.split(' = ')[0] + f' = ["{data_set}"]\n'
        elif "out_dir = " in line:
            line = line.split(' = ')[0] + f' = "{working_dir}/harvest/"\n'
        elif "station_kind = " in line:
            line = line.split(' = ')[0] + f' = "{stat_kind}"\n'
        elif "processes = " in line:
            line = line.split(' = ')[0] + f' = 30\n'
        elif "min_year_to_process = " in line:
            line = line.split(' = ')[0] + f' = 2025\n'
        elif "max_year_to_process = " in line:
            line = line.split(' = ')[0] + f' = 2026\n'
        modified_content.append(line)

    # Write the modified content to a new file
    with open(new_file_path, 'w') as file:
        file.writelines(modified_content)

    print(f"Modified file has been saved to {new_file_path}")

    # Run the harvester
    os.system(f'mkdir -p {working_dir}/harvest/')
    harvester_functions = f'{ceuas_dir}/public/harvest/code_cop2/run_harvest_convert_to_netCDF_yearSplit.py'
    print(f"python {harvester_functions} -p {working_dir}/code/modded_harvester_yearsplit_parameters.py -s {working_dir}/code/station_configuration/ -c {ceuas_dir}")
    print()
    os.system(f"module load odc; python {harvester_functions} -p {working_dir}/code/modded_harvester_yearsplit_parameters.py -s {working_dir}/code/station_configuration/ -c {ceuas_dir}")
    # module load anaconda3 \n source activate uvn10 \n

if __name__ == '__main__':

    # download_data_igra2()
    # download_data_era5()
    # copy_tables_to_harvest()
    # create_inventory('igra2')
    # create_inventory('era5_1')
    # make_station_configuration('igra2')
    # make_station_configuration('era5_1')
    # run_harvester('igra2')
    # run_harvester('era5_1')
    # run_harvester('era5_1_mobile', stat_kind='mobile')
    # run_harvester('igra2_mobile', stat_kind='mobile')

    


