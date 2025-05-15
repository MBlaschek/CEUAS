import os
import sys
import getpass
import time
import psutil
import pandas as pd
import numpy as np
import datetime
import glob
import calendar
import h5py
import h5netcdf

'''
change: set year in param file to the one from date check 
add: era5_1 download and processing
'''

# Create the directory to store the data for the current date:

global user 
user = getpass.getuser()

global reference_file
global ceuas_dir
global python_interpreter
global base_dir
global rscratch
global refs

reference_file = '/mnt/users/scratch/leo/scratch/converted_v29/long/0-20001-0-11035_CEUAS_merged_v3.nc'
ceuas_dir = '/srvfs/home/uvoggenberger/CEUAS/CEUAS/'
base_dir = '/mnt/users/scratch/uvoggenberger/CUON_HARVEST'
python_interpreter = '/srvfs/home/uvoggenberger/micromamba/envs/uv12/bin/python'
rscratch = '/mnt/users/scratch/leo/scratch/'
refs = '/mnt/users/scratch/leo/scratch/converted_v13/rea/refs1940x.pkl'

global ecmwf_user
global ecmwf_output_dir

ecmwf_user = 'lh4'
ecmwf_output_dir = '/ec/res4/scratch/lh4/'

datetime_now = datetime.datetime.now()
global date_year
date_year = datetime_now.strftime('%Y')
global date_month
date_month = datetime_now.strftime('%m')

global download_year
if int(date_month) == 1:
    download_year = int(date_year) - 1
else:
    download_year = int(date_year)
global download_month
download_month = int(date_month) - 1

date_now = datetime_now.strftime('%Y%m')

global working_dir
working_dir = base_dir + '_' + date_now + '/'
os.system('mkdir -p ' + working_dir)

global table_dir
table_dir = f'{ceuas_dir}/meta/inventory_comparison_2/data/tables/'


sys.path.append(f"{ceuas_dir}/public/cds-backend/code/")
from harvest_convert_to_netCDF import write_dict_h5


def wait_for_python_processes(com_line_content = ''):
    while True:
        time.sleep(10)
        current_processes = []
        for p in psutil.process_iter():
            if p.username() == user:
                if np.any([com_line_content in pi for pi in p.cmdline()]): 
                    current_processes.append(p)
        # Check for running Python processes
        if len(current_processes) == 0:  # No Python processes are running
            print("No Python processes running. Proceeding with the script.")
            break

def download_data_igra2(rm_zip=False):
    # Download data from NOAA
    igra_dir = working_dir + '/data/igra_data'
    if len(glob.glob(igra_dir + '/*.txt')) < 1:
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

    import calendar

def days_in_month(year: int, month: int) -> int:
    if 1 <= month <= 12:
        return calendar.monthrange(year, month)[1]
    else:
        raise ValueError("Month must be between 1 and 12")

def download_data_era5(rm_zip=False):

    ## For now done manually, as there is a 2FA for the login to the server
    ## Don't forget to also download the gridded data -> necessary for later steps
    
    # module load teleport
    # ssh-agentreconnect -> or if that doesn't work tsh login
    # ssh lh4@hpc-login

    era5_dir = working_dir + '/data/era5_1_data'
    os.system('mkdir -p ' + era5_dir)

    with open(f'{ceuas_dir}/public/nrt_pipeline/request_era5_default.ksh', 'r') as file:
        content = file.readlines()
    # Modify the content as needed
    modified_content = []
    for line in content:
        if "MMM=" in line:
            line = f"MMM={str(download_month).zfill(2)} \n"
        if "YYY=" in line:
            line = f"YYY={str(download_year)} \n"
        if "DDD=" in line:
            line = f"DDD={days_in_month(download_year, download_month)} \n"
        modified_content.append(line)
    # Write the modified content to a new file
    with open(f'{working_dir}/data/request_era5.ksh', 'w') as file:
        file.writelines(modified_content)


    with open(f'{ceuas_dir}/public/nrt_pipeline/request_era5_gridded_1_default.ksh', 'r') as file:
        content = file.readlines()
    # Modify the content as needed
    modified_content = []
    for line in content:
        if "MMM=" in line:
            line = f"MMM={str(download_month).zfill(2)} \n"
        if "YYY=" in line:
            line = f"YYY={str(download_year)} \n"
        if "DDD=" in line:
            line = f"DDD={days_in_month(download_year, download_month)} \n"
        modified_content.append(line)
    # Write the modified content to a new file
    with open(f'{working_dir}/data/request_era5_gridded_1.ksh', 'w') as file:
        file.writelines(modified_content)


    with open(f'{ceuas_dir}/public/nrt_pipeline/request_era5_gridded_2_default.ksh', 'r') as file:
        content = file.readlines()
    # Modify the content as needed
    modified_content = []
    for line in content:
        if "MMM=" in line:
            line = f"MMM={str(download_month).zfill(2)} \n"
        if "YYY=" in line:
            line = f"YYY={str(download_year)} \n"
        if "DDD=" in line:
            line = f"DDD={days_in_month(download_year, download_month)} \n"
        modified_content.append(line)
    # Write the modified content to a new file
    with open(f'{working_dir}/data/request_era5_gridded_2.ksh', 'w') as file:
        file.writelines(modified_content)

   
    with open(f'{ceuas_dir}/public/nrt_pipeline/remote_creator.sh', 'r') as file:
        content = file.readlines()
    # Modify the content as needed
    modified_content = []
    for line in content:
        if "file_to_modify_1 " in line:
            line = line.replace("file_to_modify_1", f'{working_dir}/data/request_era5.ksh')
        if "file_to_modify_2 " in line:
            line = line.replace("file_to_modify_2", f'{working_dir}/data/request_era5_gridded_1.ksh')
        if "file_to_modify_3 " in line:
            line = line.replace("file_to_modify_3", f'{working_dir}/data/request_era5_gridded_2.ksh')
        if "harvest_dir " in line:
            line = line.replace("harvest_dir", f'{era5_dir}/')
        if "ecmwf_user " in line:
            line = line.replace("ecmwf_user", f'{ecmwf_user}')
        if "ecmwf_output_dir " in line:
            line = line.replace("ecmwf_output_dir", f'{ecmwf_output_dir}')
        modified_content.append(line)
    # Write the modified content to a new file
    with open(f'{working_dir}/data/get_era5.sh', 'w') as file:
        file.writelines(modified_content)

    os.system(f'chmod +x {working_dir}/data/get_era5.sh')
    os.system(f'{working_dir}/data/get_era5.sh')

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
    wait_for_python_processes(com_line_content = 'analyze_inventory_functions')

def make_station_configuration(data_set):
    # Create the station configuration
    os.system(f'mkdir -p {working_dir}/code/station_configuration/')
    make_station_configuration_functions = f'{ceuas_dir}/meta/inventory_comparison_2/code/make_station_configuration.py'
    os.system(f"python {make_station_configuration_functions} -d {data_set} -w {working_dir}")
    try:
        o_file = glob.glob(f"{working_dir}/code/station_configuration/{data_set}_orphan*.csv")[0]
        os.system(f"cp {o_file} {o_file.replace('orphans', 'mobile')}")
    except:
        pass

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
    wait_for_python_processes(com_line_content = 'harvest_convert_to_netCDF_yearSplit')

def set_up_merge():
    os.system(f'mkdir -p {working_dir}/merge/')
    os.system(f'mkdir -p {working_dir}/merge/merged_out/')

    os.system(f'mkdir -p {working_dir}/harvest/harvest_mobile/')
    os.system(f'mkdir -p {working_dir}/harvest/harvest_regular/')
    os.system(f'mv {working_dir}/harvest/era5_1 {working_dir}/harvest/harvest_regular/')
    os.system(f'mv {working_dir}/harvest/igra2 {working_dir}/harvest/harvest_regular/')
    os.system(f'mv {working_dir}/harvest/igra2_mobile {working_dir}/harvest/harvest_mobile/')
    os.system(f'mv {working_dir}/harvest/era5_1_mobile {working_dir}/harvest/harvest_mobile/')

    merge_params = f'{working_dir}/merge/modded_merging_yearly_parameters.py'
    os.system(f'cp {ceuas_dir}/public/merge/merging_yearly_parameters.py {merge_params}')

    with open(merge_params, 'r') as file:
        content = file.readlines()

    # Modify the content as needed
    modified_content = []
    for line in content:
        if "harvested_base_dir = " in line and not '    h' in line:
            line = line.split(' = ')[0] + f' = "{working_dir}/harvest/harvest"\n'
        elif "merged_out_dir = " in line:
            line = line.split(' = ')[0] + f' = "{working_dir}/merge/merged_out/"\n'
        elif "                                   " in line and not 'era5_1_' in line and not 'igra2' in line:
            line = ''
        modified_content.append(line)

    # Write the modified content to a new file
    with open(merge_params, 'w') as file:
        file.writelines(modified_content)

def run_merge(station_kind = "regular"):

    merge_params = f'{working_dir}/merge/modded_merging_yearly_parameters.py'
    with open(merge_params, 'r') as file:
        content = file.readlines()
    modified_content = []
    for line in content:
        if "station_kind = " in line:
            line = line.split(' = ')[0] + f' = "{station_kind}"\n'
        modified_content.append(line)
    with open(merge_params, 'w') as file:
        file.writelines(modified_content)

    os.system(f"python {ceuas_dir}/public/merge/merging_cdm_netCDF_yearSplit_SEP2023_pipeline.py -p {merge_params}  -max_y {date_year} -min_y {date_year}")
    
def run_resort():
    os.system(f'mkdir -p {working_dir}/resort/')
    os.system(f'mkdir -p {working_dir}/resort/long')
    os.system(f'cp {refs} {working_dir}/resort/')
    # set correct path to CUON statconf
    os.system(f"python {ceuas_dir}/public/resort/convert_and_resort_pipeline.py -i {working_dir}/merge/merged_out/ -w {working_dir} -c {ceuas_dir} -r {rscratch}")

def add_tables():
    path = f'{working_dir}/resort/{date_year}/'
    all_files = glob.glob(f'{path}/*.nc')
    attr_file = glob.glob(f'{reference_file}')

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

    # multiprocess here! 
    for file in all_files:
        result = check_and_fix_file(file, attr_dict)
        if result == 1:
            print(file, ": groups have been added!")

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



if __name__ == '__main__':

    # os.system(f'module load teleport; python3 -m teleport.login ') # export TSH_EXEC=/home/swd/manual/teleport/17.4.2/bin/tsh ; export TSH_PROXY=jump-17.ecmwf.int:443 ; export TSH_USERNAME=lh4 ; export TSH_PASSWORD=Ax9eZ3XHzW8T8Em ; 

    # Define the path to the marker file
    marker_file = os.path.join(working_dir, "download_complete.txt")

    # Check if the marker file exists
    if not os.path.exists(marker_file):
        print("Marker file not found. Running download functions...")
        download_data_igra2()
        download_data_era5()

        # Create the marker file to indicate completion
        with open(marker_file, "w") as f:
            f.write("Files prepared.\n")
        print("Marker file created.")

    print("Marker file found. Skipping download functions.")

    ## Call the following functions:

    copy_tables_to_harvest()
    create_inventory('igra2')
    create_inventory('era5_1')
    make_station_configuration('igra2')
    make_station_configuration('era5_1')
    run_harvester('igra2')
    run_harvester('era5_1')

    run_harvester('era5_1_mobile', stat_kind='mobile')
    run_harvester('igra2_mobile', stat_kind='mobile')

    set_up_merge()
    run_merge('regular')
    run_merge('mobile')
    run_merge('orphan')

    make_station_configuration("CUON")

    run_resort()

    add_tables()







