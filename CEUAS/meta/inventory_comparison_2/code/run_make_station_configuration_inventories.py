import os,sys


#os.system('python3.8 make_station_configuration.py -d CUON ')


#os.system('python3.8 make_station_configuration.py -d MERGE')



ds = [ 'era5_1', 'era5_2', 'era5_1759', 'era5_1761',
       'era5_3188', 'bufr', 'giub',
       'ncar', 'igra2',
       'amma', 'hara' ]

ds_mobile = ['era5_1_mobile' , 'era5_2_mobile', 'npsound' ]


for d in ds:
    os.system('python3.8 make_station_configuration.py -d ' + d )

'''
for dd in ds_mobile:
    os.system('python3.8 make_station_configuration.py -d ' + dd )
'''

