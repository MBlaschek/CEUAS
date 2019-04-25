# Introduction 
This README explains the content of the test_dir directory,
that includes the main scripts and relative information to:
1. convert the odb file to netCDF files
   (u,v wind , temperature-t, dew point temperature-dp, relative humidity-rh) 
2. plot the data from the netCDF files for u,v wind, t, dp, rh
3. process the netCDF files and calculate, when possible,
   the dp and rh from t,rh and t,dp if available,
   and create a 'combined' version of the file
4. make plots to check the combined values

### 1. readodbstationfiles_10393.py
This scripts reads the odb files for the station 10393 
and extracts the netCDF files.
The files are stored in the directory 'out_netCDFs'
The variables processed are temperature, relative humidity, dew point temp. and wind speed (u,v)
The data for each variable is stored in the files:
ERA5_1_10393_t.nc
ERA5_1_10393_rh.nc
ERA5_1_10393_dp.nc
ERA5_1_10393_u.nc
ERA5_1_10393_v.nc

### 2. plots_variables_10393.py
This scripts produces the plots for the temperature, relative humidity, wind (u,v) and dew point
form the netCDFs files, for each observation hour [0,1] and standard pressure level [0-15].
The plots are stored inside the 'Plots_10393' directory

### 3. 
