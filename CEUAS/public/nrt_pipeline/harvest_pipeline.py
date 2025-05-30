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
import hdf5plugin
import h5netcdf
from datetime import datetime, timedelta
import glob
import ray


####
# ADD LOGGING
# CHECK WAIT FOR PROCESSES
####

'''
This script is used to download data from NOAA and ECMWF, create inventories, run the harvester, merge the data, and resort it.
It is a part of the CEUAS project and is used to process data from the CUON network.
The script is divided into several functions, each responsible for a specific task.

The main functions are:
- download_data_igra2: Downloads data from NOAA and unzips it.
- download_data_era5: Downloads data from ECMWF and unzips it.
- create_inventory: Creates an inventory of the data.
- make_station_configuration: Creates a station configuration file.
- run_harvester: Runs the harvester to process the data.
- set_up_merge: Sets up the merging process.
- run_merge: Merges the data.
- run_resort: Resorts the data.
- add_tables: Adds tables to the data.

The script uses the following libraries:
- os: For file and directory operations.
- sys: For system-specific parameters and functions.
- getpass: For getting the username of the current user.
- time: For time-related functions.
- psutil: For process and system utilities.
- pandas: For data manipulation and analysis.
- numpy: For numerical operations.
- datetime: For date and time manipulation.
- glob: For file name pattern matching.
- calendar: For calendar-related functions.
- h5py: For reading and writing HDF5 files.
- h5netcdf: For reading and writing NetCDF files.

The script is designed to be run on a server with the necessary libraries and permissions.

'''
####       ####
# USER SETUP: #
#             #

global user 
user = getpass.getuser()

global ceuas_dir
global python_interpreter
global base_dir
# global rscratch
# global refs

ceuas_dir = '/srvfs/home/uvoggenberger/CEUAS/CEUAS/' # path to the CEUAS directory
base_dir = '/mnt/users/scratch/uvoggenberger/CUON_HARVEST' # path to the base directory for the harvest -> make sure to have enough disk space (100 GB per month)
python_interpreter = '/srvfs/home/uvoggenberger/micromamba/envs/uv12/bin/python' # path to the python interpreter

global ecmwf_user
global ecmwf_output_dir

ecmwf_user = 'lh4'
ecmwf_output_dir = '/ec/res4/scratch/lh4/'

#            #
##############  

global reference_file
reference_file = f'{ceuas_dir}/public/nrt_pipeline/0-20000-0-01107_CEUAS_merged_v3.nc' # A file for data structure reference.


###
# DATE SELECTION
auto_date = False # Set to True to automatically set the date to the previous month
selected_year = 2025
selected_month = 1
#
###

# GLOBAL VARIABLES ##############

global date_year
global date_month
global download_year
global download_month


if auto_date:
    datetime_now = datetime.datetime.now()
    date_year = datetime_now.strftime('%Y')
    date_month = datetime_now.strftime('%m')

    if int(date_month) == 1:
        download_year = int(date_year) - 1
        download_month = 12
        
    else:
        download_year = int(date_year)
        download_month = int(date_month) - 1
else:
    download_year = selected_year
    download_month = selected_month

date_year = download_year
date_month = download_month


date_now = str(download_year) + str(download_month).zfill(2)

global working_dir
os.system('mkdir -p ' + base_dir)
working_dir = base_dir + '/' + date_now + '/'
os.system('mkdir -p ' + working_dir)
os.system(f'mkdir -p {working_dir}/logs/')


global table_dir
table_dir = f'{ceuas_dir}/meta/inventory_comparison_2/data/tables/'


sys.path.append(f"{ceuas_dir}public/harvest/code_cop2/")
from harvest_convert_to_netCDF import write_dict_h5


def wait_for_python_processes(com_line_content = ''):
    zombie = 0
    while zombie < 5:
        time.sleep(10)
        current_processes = []
        for p in psutil.process_iter():
            if p.username() == user:
                try:
                    if np.any([com_line_content in pi for pi in p.cmdline()]): 
                        current_processes.append(p)
                except (psutil.ZombieProcess):
                    zombie += 1 
        # Check for running Python processes
        if len(current_processes) == 0:  # No Python processes are running
            print("No Python processes running. Proceeding with the script.")
            break

def download_data_igra2(rm_zip=False):
    # Download data from NOAA
    igra_dir = working_dir + '/data/igra2_data'
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
    
def monitor_odc_split():
    start_time = datetime.now()
    print("Waiting for a process with 'odc split' in its command line...")
    
    process = None
    while process is None:
        for proc in psutil.process_iter(attrs=['pid', 'cmdline', 'create_time']):
            try:
                if "odc split" in ' '.join(proc.info['cmdline']):
                    process = proc
                    break
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                continue

    print(f"Found process PID={process.pid}. Monitoring runtime...")

    try:
        while process.is_running():
            runtime = datetime.now() - start_time
            print(f"Runtime: {str(runtime).split('.')[0]}", end='\r')
            time.sleep(2)
    except psutil.NoSuchProcess:
        pass

    print(f"\nProcess PID={process.pid} has exited. Total runtime: {str(datetime.now() - start_time).split('.')[0]}")


def download_data_era5(rm_zip=False):  

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
        if "harvest_dir" in line:
            line = line.replace("harvest_dir", f'{era5_dir}/')
        if "ecmwf_user" in line:
            line = line.replace("ecmwf_user", f'{ecmwf_user}')
        if "ecmwf_out_dir" in line:
            line = line.replace("ecmwf_out_dir", f'{ecmwf_output_dir}')
        if "output.txt" in line:
            line = line.replace("output.txt", f'{working_dir}/logs/ecmwf_log.txt')
        modified_content.append(line)
    # Write the modified content to a new file
    with open(f'{working_dir}/data/get_era5.sh', 'w') as file:
        file.writelines(modified_content)

    os.system(f'chmod +x {working_dir}/data/get_era5.sh')
    os.system(f'{working_dir}/data/get_era5.sh')

    with open(f'{ceuas_dir}/public/nrt_pipeline/split.ksh', 'r') as file:
        content = file.readlines()
    # Modify the content as needed
    modified_content = []
    for line in content:
        if "YYYYMM" in line:
            line = line.replace("YYYYMM", f'{date_now}')
        if "PDIR" in line:
            line = line.replace("PDIR", f'{era5_dir}')
        modified_content.append(line)
    # Write the modified content to a new file
    with open(f'{working_dir}/data/split.ksh', 'w') as file:
        file.writelines(modified_content)
    
    os.system(f'chmod +x {working_dir}/data/split.ksh')
    os.system(f'module load odc; {working_dir}/data/split.ksh | tee {working_dir}/logs/splitting_era5_log.txt')

    ray.init(num_cpus=40)

    @ray.remote
    def odbmobile(fn):
        print(fn)
        return os.system(f"module load odc; odc sql --full-precision -q 'select *' -i "+f'"{fn}" | tr -d " "  | gzip   > "{fn}.gz"')
    
    @ray.remote
    def odb(fn):
        print(fn)
        return os.system(f"module load odc; odc sql --full-precision -q 'select *' -i "+f'"{fn}" | tr -d " "  | gzip   > "{fn}.gz"')
   
    fns=[]
    year='????'
    month='??'
    for patt in f'{era5_dir}/era5.conv.{year}{month}.[0-9]????',:
        print(patt)
        fns+=glob.glob(patt)
    
    futures = [odb.remote(fn) for fn in fns]
    results = ray.get(futures)
    
    fns=[]
    year='????'
    month='??'
    for patt in f'{era5_dir}/era5.conv.{year}{month}.[A-Z]????',\
        f'{era5_dir}/era5.conv.{year}{month}.????', f'{era5_dir}/era5.conv.{year}{month}.??????',\
        f'{era5_dir}/era5.conv.{year}{month}.???????',  f'{era5_dir}/era5.conv.{year}{month}.????????':
        print(patt)
        fns+=glob.glob(patt)
        fns = [fn for fn in fns if not '.gz' in fn]
    
    futures = [odbmobile.remote(fn) for fn in fns]
    results = ray.get(futures)

    ray.shutdown()



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
    analyze_inventory_functions = f'{ceuas_dir}/meta/inventory_comparison_2/code/analyze_inventory_functions_pipeline.py'
    print(f"python {analyze_inventory_functions} -d {data_set} -w {working_dir} -i {input_dir}")
    os.system(f"module load odc; python {analyze_inventory_functions} -d {data_set} -w {working_dir} -i {input_dir}") 
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
            line = line.split(': ')[0] + f': "{working_dir}/data/igra2_data",\n'
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
            line = line.split(' = ')[0] + f' = {date_year}\n'
        elif "max_year_to_process = " in line:
            line = line.split(' = ')[0] + f' = {str(int(date_year)+1)}\n'
        modified_content.append(line)

    # Write the modified content to a new file
    with open(new_file_path, 'w') as file:
        file.writelines(modified_content)

    print(f"Modified file has been saved to {new_file_path}")

    # Run the harvester
    harvester_functions = f'{ceuas_dir}/public/harvest/code_cop2/run_harvest_convert_to_netCDF_yearSplit.py'
    print(f"python {harvester_functions} -p {working_dir}/code/modded_harvester_yearsplit_parameters.py -s {working_dir}/code/station_configuration/ -c {ceuas_dir}")
    print()
    os.system(f"module load odc; python {harvester_functions} -p {working_dir}/code/modded_harvester_yearsplit_parameters.py -s {working_dir}/code/station_configuration/ -c {ceuas_dir} | tee {working_dir}/logs/harvester_log.txt")
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

    os.system(f"python {ceuas_dir}/public/merge/merging_cdm_netCDF_yearSplit_SEP2023_pipeline.py -p {merge_params}  -max_y {date_year} -min_y {date_year} | tee {working_dir}/logs/merging_log.txt")
    
def run_resort():
    os.system(f'mkdir -p {working_dir}/resort/')
    os.system(f'mkdir -p {working_dir}/resort/long')
    # os.system(f'cp {refs} {working_dir}/resort/')
    # set correct path to CUON statconf
    os.system(f"python {ceuas_dir}/public/resort/convert_and_resort_pipeline.py -i {working_dir}/merge/merged_out/ -w {working_dir} -c {ceuas_dir} | tee {working_dir}/logs/resort_log.txt") # -r {rscratch}

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
            ov_vars = [float(0)] * len(f['observations_table']['date_time'][:])
            missing_vars = target_vars
        elif len(set(target_vars) - set(f['advanced_homogenisation'].keys())) > 0:
            ov_vars = [float(0)] * len(f['observations_table']['date_time'][:])
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

    # Define the path to the marker file
    marker_file = os.path.join(working_dir, "download_complete.txt")

    # Check if the marker file exists
    if not os.path.exists(marker_file):
        print("Marker file not found. Running download functions...")
        download_data_era5()
        download_data_igra2()
        # Create the marker file to indicate completion
        with open(marker_file, "w") as f:
            f.write("Files prepared.\n")
        print("Marker file created.")

    print("Marker file found. Skipping download functions.")

    # Call the following functions:

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


