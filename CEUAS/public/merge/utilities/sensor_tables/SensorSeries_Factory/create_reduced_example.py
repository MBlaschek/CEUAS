import os,sys

""" Create example directory """
ex = 'example'


if not os.path.isdir(ex):
    os.mkdir(ex)
if not os.path.isdir(ex + '/data_plots'):
    os.mkdir(ex + '/data_plots')

    
for f in ['sensor_metadata_DASHBOARD.ipynb',
          'sensor_metadata_DASHBOARD.py',
          'sensor_configuration_all.csv',
          'igra2-metadata.txt',
          'modules']:

    os.system('cp -r ' + f + '  ' + ex )


""" copying the station data """
stations = ['0-20000-0-06610' ,
            '0-20001-0-10393' ,
            '0-20001-0-11035' ,
            '0-20000-0-82930' ]

import glob
for s in stations:
    files = glob.glob ('data_plots/*' + s + "*")
    for f in files:
        os.system('cp ' + f + ' ' + ex + '/data_plots/' )
