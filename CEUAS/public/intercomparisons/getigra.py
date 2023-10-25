import numpy
import numpy as np
import pandas as pd
import sys, glob, os
import urllib3
import h5py
import cdsapi, zipfile, os, time
sys.path.append(os.getcwd()+'/../cds-backend/code/')
import cds_eua3 as eua
import warnings
import shutil
import xarray
import pickle
warnings.filterwarnings('ignore')
'download.csv-lev.zip'
def request(rqdict, source, sdir):
    t0 = time.time()

    c = cdsapi.Client()
    r = c.retrieve(
        source,rqdict, 'download.csv-lev.zip')
    print('Request took: ' + str(time.time() - t0) + ' seconds')
    if True:
        r.download(target='csv-lev.zip')
        assert os.stat('csv-lev.zip').st_size == r.content_length, "Downloaded file is incomplete"
    z = zipfile.ZipFile('csv-lev.zip')
    z.extractall(path=sdir)
    z.close()
    return

# ---

igra = 'insitu-observations-igra-baseline-network'

# ---


for i in range(1979,2021):
    print('IGRA', i)
    s2dir = './download/IGRA_' + str(i) + '/'
#     s2dir = './download/IGRA_H_' + str(i) + '/'
    if not os.path.isdir(s2dir):
        try:
            request({'source': 'global_radiosonde_archive',
#             request({'source': 'harmonized_global_radiosonde_archive',
                     'format': 'csv-lev.zip',
#                      'variable': [
#                          'air_temperature', 'eastward_wind_component', 'northward_wind_component',
#                          'relative_humidity',
#                      ],
                     'variable': [
                         'air_dewpoint_depression', 'air_temperature', 'geopotential_height',
                         'relative_humidity', 'wind_from_direction', 'wind_speed',
                     ],
                     'year': str(i),
                     'month': ['01','02','03','04','05','06','07','08','09','10','11','12'],
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
                    },igra, sdir=s2dir,
                   )
        except:
            pass

# def filterigra(yrs):
#     stations = pickle.load( open( "../intercomparisons/stations.p", "rb" ))
#     for j in [yrs]:
#         print(j)
#         file = glob.glob('/raid60/scratch/uli/IGRA_Data/IGRA/'+str(j)+'/*')
#         data = pd.read_csv(file[0], header=12, index_col=False)
#         data = data[data.report_timestamp > str(j-1)+'-12-31 23:59:59+00'][data.report_timestamp < str(j+1)+'-01-01 00:00:00+00']
#         for i in stations.station_name:
#             out = data[data.station_name == i]
#             out = out[out.air_pressure == 10000]
#             pickle.dump( out, open( '/raid60/scratch/uli/IGRA_Data/IGRA/IGRA_' + str(i) + '_' + str(j) + '_100.p', "wb" ))
#         print('done')

# import multiprocessing

# pool = multiprocessing.Pool(processes=3)
# result_list = pool.map(filterigra, [1982,1983,1984])        
        
        
# stations = pickle.load( open( "../intercomparisons/stations.p", "rb" ))

# for j in range(1979,2020,1):
#     print(j)
#     file = glob.glob('/raid60/scratch/uli/IGRA_Data/IGRA/'+str(j)+'/*')
#     data = pd.read_csv(file[0], header=12, index_col=False)
#     data = data[data.report_timestamp > str(j-1)+'-12-31 23:59:59+00'][data.report_timestamp < str(j+1)+'-01-01 00:00:00+00']
#     for i in stations.station_name:
#         out = data[data.station_name == i]
#         out = out[out.air_pressure == 10000]
#         pickle.dump( out, open( '/raid60/scratch/uli/IGRA_Data/IGRA/IGRA_' + str(i) + '_' + str(j) + '_100.p', "wb" ))
#     print(done)


for i in stations.station_name:
    try:
        data = request({'source': 'IGRA',
                        'variable':['air_temperature'],
                        'station_name': i,
                        'period': '1979-01-01/2019-12-31',
                       }, 'insitu-observations-igra-baseline-network',
        )
        print(i,'done')
    except:
        pass
    pickle.dump( data, open( '/raid60/scratch/uli/IGRA_Data/IGRA/IGRA_' + str(i) + '_100.p', "wb" ))
    
for i in stations.station_name:
    try:
        data = request({'source': 'IGRA_H',
                        'variable':['air_temperature'],
                        'station_name': i,
                        'period': '1979-01-01/2019-12-31',
                       }, 'insitu-observations-igra-baseline-network',
        )
        print(i,'done')
    except:
        pass
    pickle.dump( data, open( '/raid60/scratch/uli/IGRA_Data/IGRA_H/IGRA_H_' + str(i) + '_100.p', "wb" ))
