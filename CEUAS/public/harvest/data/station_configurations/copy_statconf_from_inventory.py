import os,sys

direc = '../../../../meta/inventory_comparison/code/station_configuration/'

databases = [ 'era5_2', 'era5_1759', 'era5_1761', 'era5_3188', 'bufr', 'igra2', 'era5_1', 'ncar']

for d in databases:
    f = direc + '/' + d + '_station_configuration_extended.csv'
    os.system('cp ' + f + '  . ' )
