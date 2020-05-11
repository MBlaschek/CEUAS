# Make Plots
 Directory containing scripts and utilities to make pretty plots


## plot_stations_distributions.py

Script to produce nice distributions of the stations for differnet periods.

Data source:

- the file `summary_forplot.dat` 
  contains the columns: [lat, lon, start, end] i.e. coordinates and starting and ending dates of observations for the [ncar,igra2,bufr,era5_1, era5_1759, era5_1761, era5_3188] sets.
  The file is produced by the merging_dataset_netcdf.py script.

