#!/usr/bin/env python
# By Ulrich Voggenberger
# Purpose: Create and compare adjustments with those from the CDS 
import sys
import time
import os
import logging

def download_data():
    """Download data from CDS and preprepare  for adjustment procedure
    """
#     import
    print('downloading and preparing data for processing')
    try:
        os.chdir('Converters')
        os.system('python3 from_cds_to_legacy.py')
    except:
        return 'downloading and converting error'
    
    return 0
    

def prepare_environment():
    """Prepare environment for adjustment procedure
    """
#     import
    print('loading modules')
    try:
        os.system('module load openmpi/3.1.6-intel-20.0.2')
        os.system('module load netcdf-fortran/4.5.3-intel-20.0.2')
        os.system('module load intel-parallel-studio/composer.2020.2-intel-20.0.2')
    except:
        return 'module loading error'
    
    print('setting stacksize to unlimited')
    try:
        os.system('limit stacksize unlimited')
    except:
        return 'stacksize setting error'
    
    print('compile fortran code')
    try:
        os.chdir('RISE_FORTRAN')
        os.system('make raso_correct_nc')
        os.system('setenv LD_LIBRARY_PATH /usr/local/lib')
    except:
        return 'fortran compiling error'
    
    print('fetching data')
    try:
        os.system('python3 import_structure.py')
    except:
        return
    
    return 0

def create_adj():
    """Run adjustment procedure and save to file
    """
#     import
    print('creating adjustments')
    try:
        os.chdir('RISE_FORTRAN')
        os.system('raso_correct_nc radcorpar06')
    except:
        return 'adjustment creating error'
    return 0

def compare_adj(station):
    """Compare created adjustment to downloaded adjustment
    """
    file = ('/Temperature_adjustment/0'+station+'/feedbackglobbincorrsave0'+station+'.nc')
    # TODO: compare with downloaded file
    
    return 0
    
if __name__ == "__main__":
    compare_list = ['15480']
    HOME = os.getcwd()
    print(HOME)
    
    status = prepare_environment()
    os.chdir(HOME)
    if status = 0:
        status = download_data()
    else:
        sys.exit(status)
    os.chdir(HOME)
    if status = 0:
        status = create_adj()
    else:
        sys.exit(status)
    os.chdir(HOME)
    if status = 0:
        for i in compare_list:
            status = compare_adj(i)
            if status != 0:
                break
    else:
        sys.exit(status)
    os.chdir(HOME)
    
    if status = 0:
        print('---')
        print('comparison successful')
        print('---')
    else:
        sys.exit(status)
    
    