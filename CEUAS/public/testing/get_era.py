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

# import cdsapi
# c = cdsapi.Client()
# for yr in range(1979,2022,1):
#     print(yr)
#     r = c.retrieve(
#         'reanalysis-era5-single-levels-monthly-means',
#         {
#             'format': 'netcdf',
#             'product_type': 'monthly_averaged_reanalysis',
#             'variable': [
#                 '10m_u_component_of_wind', '10m_v_component_of_wind', '2m_dewpoint_temperature',
#                 '2m_temperature', 'skin_temperature', 'soil_type',
#                 'surface_pressure',
#             ],
#             'year': yr,
#             'month': [
#                 '01', '02', '03',
#                 '04', '05', '06',
#                 '07', '08', '09',
#                 '10', '11', '12',
#             ],
#             'time': '00:00',
#         },
#         'download.nc')
#     r.download(target='era_'+ str(yr)+'.nc')
#     print('done')
    

import cdsapi
c = cdsapi.Client()
for yr in range(1950,2023,1):
    try:
        print(str(yr))
        r = c.retrieve(
            'reanalysis-era5-land-monthly-means',
            {
                'format': 'netcdf',
                'product_type': 'reanalysis-monthly-means-of-daily-means',
                'variable': [
                    '10m_u_component_of_wind', '10m_v_component_of_wind', '2m_dewpoint_temperature',
                    '2m_temperature', 'skin_temperature', 'soil_type',
                    'surface_pressure',
                ],
                'year': str(yr),
                'month': [
                    '01', '02', '03',
                    '04', '05', '06',
                    '07', '08', '09',
                    '10', '11', '12',
                ],
                'time': [
                    '00:00', '12:00',
                ],
            },
            'download.nc')
        r.download(target='./era/era_'+ str(yr)+'.nc')
        print('done')
    except:
        pass
    