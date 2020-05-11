""" This script reads the station_configurations files in the given directory,
and for each database it extracts the lat,lon, start and end dates of the obserations,
writing a separate file for each. """


import pandas as pd 
import sys, os



os.system('mkdir summaries')


""" Loop over the databases """
ds = ['era5_1','era5_2','era5_1759','era5_1761','era5_3188','bufr','ncar', 'igra2' ]
dir = '/raid8/srvx1/federico/GitHub/CEUAS_master_MAY/CEUAS/CEUAS/public/harvest/data/station_configurations/'

for d in ds:
    print(' *** Processing the database: ',  d )
    out = open('summaries/' + d + '_summary_distribution.dat','w')
    out.write('start\tend\tlat\tlon\n')
    f = dir + 'station_configuration_' + d + '.dat'
    df = pd.read_csv(f, delimiter = '\t')
    for index, row in df.iterrows():
        start, end = str(row['start_date']), str(row['end_date'])
        lat, lon = str(row['latitude']) , str(row['longitude'])
        out.write(start + '\t' + end + '\t' + lat + '\t' +  lon + '\n')
    out.close()
    print('*** Done with: ', d )
