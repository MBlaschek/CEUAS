import os
import sys
import subprocess
import pandas as pd
import numpy as np
import datetime
import glob

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
    analyze_inventory_functions = f'{ceuas_dir}/meta/inventory_comparison_2/code/analyze_inventory_functions.py'
    os.system(f"python {analyze_inventory_functions} -d {data_set} -w {working_dir}") # somehow selects the old dir -> fix this

def make_station_configuration(data_set):
    # Create the station configuration
    os.system(f'mkdir -p {working_dir}/code/station_configuration/')
    make_station_configuration_functions = f'{ceuas_dir}/meta/inventory_comparison_2/code/make_station_configuration.py'
    os.system(f"python {make_station_configuration_functions} -d {data_set} -w {working_dir}")

def run_harvester(data_set):
    # Set up the parameter file:
    original_file_path = f'{ceuas_dir}/public/harvest/code_cop2/harvester_yearsplit_parameters.py'
    new_file_path = f'{working_dir}/code/modded_harvester_yearsplit_parameters.py'
    os.system(f'mkdir -p {working_dir}/harvest/')

    for stat_kind in ['regular']: # , 'mobile', 'orphan']:
        if data_set == 'igra2' and stat_kind == 'mobile':
            continue

        # Read the content of the original file
        with open(original_file_path, 'r') as file:
            content = file.readlines()

        # Modify the content as needed
        modified_content = []
        for line in content:
            if "'igra2': " in line or "'igra2_mobile': " in line:
                line = line.split(': ')[0] + f': "{working_dir}/data/igra_data",\n'
            elif "datasets = " in line:
                line = line.split(' = ')[0] + f' = ["{data_set}"]\n'
            elif "out_dir = " in line:
                line = line.split(' = ')[0] + f' = "{working_dir}/harvest/"\n'
            elif "station_kind = " in line:
                line = line.split(' = ')[0] + f' = "{stat_kind}"\n'
            elif "processes = " in line:
                line = line.split(' = ')[0] + f' = 30\n'
            elif "min_year_to_process = " in line:
                line = line.split(' = ')[0] + f' = 2024\n'
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
        os.system(f"python {harvester_functions} -p {working_dir}/code/modded_harvester_yearsplit_parameters.py -s {working_dir}/code/station_configuration/ -c {ceuas_dir}")
        # module load anaconda3 \n source activate uvn10 \n

if __name__ == '__main__':

    # download_data_igra2()
    # copy_tables_to_harvest()
    # create_inventory('igra2')
    # make_station_configuration('igra2')
    run_harvester('igra2')

