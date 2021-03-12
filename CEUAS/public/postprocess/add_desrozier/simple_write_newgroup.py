# Modul Laden
import os,sys

sys.path.append('../../cds-backend/code/')
import cds_eua3 as eua



data = eua.CDMDataset('/raid60/scratch/uli/converted_v5/0-20000-0-11035_CEUAS_merged_v0.nc')

# Daten einlesen Temperatur
xyz = data.read_observed_variable(85, return_xarray=True)

# Daten schreiben neue Variable monkey in neuer gruppe adjust
data.write_observed_data('monkey',
                                                  ragged=xyz,  # input data
                                                  varnum=85,  # observed_variable to be aligned with
                                                  group='adjust',   # name of the new group
                                                  data_time='date_time',  # named datetime coordinate
                                                  data_plevs='z_coordinate'  # named pressure coordinate
                                                 )

                        
