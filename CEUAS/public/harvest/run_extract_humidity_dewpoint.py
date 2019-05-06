""" Module to extract dew point temperature or relative humidity from the netCDF files """

from extract_dewpoint_humidity import *

# Paths to the input files with temperature, dew point and relative humidity
input_netCDF_t  = 'out_netCDFs/1/10393/ERA5_1_10393_t.nc' #input file with temperature                                                                                                                                                                                                   
input_netCDF_sh = 'out_netCDFs/1/10393/ERA5_1_10393_sh.nc'
input_netCDF_rh = 'out_netCDFs/1/10393/ERA5_1_10393_rh.nc'
input_netCDF_dp = 'out_netCDFs/1/10393/ERA5_1_10393_dp.nc'

""" Loading the files """
netCDFs = netCDF_Files(file_t=input_netCDF_t  , file_rh = input_netCDF_rh  , file_dp = input_netCDF_dp)
dates = netCDFs.load_datum()
data  = netCDFs.load_data()
plevels = netCDFs.define_plevels()

""" Loading the data, datum (as lists) """
data_t  , data_rh  , data_dp = netCDFs.data_t ,  netCDFs.data_rh  , netCDFs.data_dp
datum_t , datum_rh , datum_dp = netCDFs.datum_t , netCDFs.datum_rh , netCDFs.datum_dp
datas = netCDFs.find_datums()
check_datum = netCDFs.check_datum() # true if datums are identical, false otherwise                                                                                                                                                                                                      
comb_file = netCDFs.create_new() # file to be filled with the combined data                                                                                                                                                                                                               

""" Loop over the pressure levels """
for p in plevels.keys():
    for h in [0,1]:
        print('processing the hour, pressure' , h , p )
        rh, dp, flagrh, flagdp = netCDFs.calc_missing(h=h, p=p)
        netCDFs.write_comb(p = p, h = h)


print("Done!")
