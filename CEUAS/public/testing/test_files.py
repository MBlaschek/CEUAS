#!/usr/bin/env
# coding: utf-8

import numpy
import numpy as np
import pandas as pd
import sys, glob
import urllib3
import h5py
import cdsapi, zipfile, os, time
import warnings
import shutil
import xarray
from datetime import date
warnings.filterwarnings('ignore')
# import pycountry
sys.path.append(os.getcwd()+'/../cds-backend/code/')
import cds_eua3 as eua
import numba
import copy
import glob
from numba import njit
import pandas
import glob
import multiprocessing
from functools import partial

def testing(file):
    try:
        test = eua.CDMDataset(file)
        log = test.report_quality()
    except:
        log = {'Error': file}
    return log

if __name__ == '__main__':

    test_files = glob.glob('/raid60/scratch/leo/scratch/converted_v7/*.nc')

    pool = multiprocessing.Pool(processes=10)
    func=partial(testing)
    result_list = list(pool.map(func, test_files))
    log = []
    for i in result_list:
        for j in i.keys():
            log.append(j + ' : ' + str(i[j]) + '\n')
        log.append('\n')

    writelog=open(str(date.today()) + '_file_test_log.txt','w')
    writelog.writelines(log)
    writelog.close()