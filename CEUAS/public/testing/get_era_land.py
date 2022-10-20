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

import cdsapi
c = cdsapi.Client()
for yr in range(1977,1978,1):
    for mn in range(6,13,1):
        print(yr)
        r = c.retrieve(
            'reanalysis-era5-land',
            {
            'format': 'netcdf',
            'time': [
                '00:00', '12:00',
            ],
            'variable': [
                '10m_u_component_of_wind', '10m_v_component_of_wind', '2m_dewpoint_temperature',
                '2m_temperature', 'skin_temperature', 'snow_cover',
                'surface_pressure',
            ],
            'year': str(yr),
            'month': [ str(mn).zfill(2)], 
#                 '01', '02', '03',
#                 '04', '05', '06',
#                 '07', '08', '09',
#                 '10', '11', '12',
#             ],
            'day': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
                '13', '14', '15',
                '16', '17', '18',
                '19', '20', '21',
                '22', '23', '24',
                '25', '26', '27',
                '28', '29', '30',
                '31',
            ],
            },
            'download.nc')
        r.download(target='./era_hourly/era_land_hourly'+ str(yr)+'_'+str(mn)+'.nc')
        print('done')