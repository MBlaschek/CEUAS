import os,sys

os.system('python3.8 make_station_configuration.py -d era5_1 &')
os.system('python3.8 make_station_configuration.py -d bufr &')
os.system('python3.8 make_station_configuration.py -d era5_1759 &')
os.system('python3.8 make_station_configuration.py -d era5_1761 &')
os.system('python3.8 make_station_configuration.py -d igra2 &')
os.system('python3.8 make_station_configuration.py -d ncar &')
os.system('python3.8 make_station_configuration.py -d era5_3188 &')
os.system('python3.8 make_station_configuration.py -d era5_2 ')


os.system('python3.8 make_station_configuration.py -d CUON ')
os.system('python3.8 make_station_configuration.py -d MERGE ')


