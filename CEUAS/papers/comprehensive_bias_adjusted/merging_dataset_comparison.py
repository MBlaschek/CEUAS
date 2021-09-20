import os,sys
import h5py as h5
import pandas as pd
import netCDF4 as nc

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker




"""
Reading the homogenized dataset Zhou at Al. 2020

The variable time has lenght:
   time(752)
since the paper describes the data between 1958-2018, but the files contain also data for 2019 and up until August 2020
so that 752 = 61*12 + (12+8)

Raw temperature shape
f_HOMO["rawT"].shape
(1184, 752, 15, 2)


"""



#site_ids = [''.bjoin(i) for i in site_ids ]


"""
# might not be needed
def get_station_index(file_homo, station):
    '''
    Returns the index of the desired station in input
    from the file UA-HRD_stations_homo+raw-monthly-anomalyT_00-12Z_195801-202008.nc

            Parameters:
                    file (netCDF4 Dataset): input netCDF4 file
                    station_id (str): station identifier e.g. 
                                WIEN AUM00011035
                                LINDENBERG GMM00010393
                                PAYERNE SZM00006610
            Returns:
                    index (int): index of the station data as in the rawT variable
    '''
    
    site_ids = file_homo["SiteID"][:].data
    
    for i in range(1184):
        stat_id = site_ids[:,i]
        stat_id_str = b''.join(stat_id).decode('utf-8')
        if station == stat_id_str:
            return i
    
    print('*** Could not find the desired station, please check station parameter! ')
"""

def get_station_data(file_sub, hour, plevel):
    '''
    Retrieves the temperature time series from a subdaily file e.g. homo-raw-subdaily-station/AUM00011035.nc
    (shape rawT(time, pressure, UTC) )
            Parameters:
                    plevel (int):  index of the desired plevel, see below for mapping pascal-index

                                        Pressure: [  1000,   2000,   3000,   5000,   7000,  10000,  15000,
                                        20000,  25000,  30000,  40000,  50000,  70000,  85000,
                                        92500, 100000, 101300 ]
                                        
                                        50000: index 11
                                        70000: index 12
            
                    hour (int): index of the desired hour, see below for mapping hour-index
                    
                                    UTC: [0,12]
                                    0  : index 0
                                    12: index 1
    '''
    
    temp = file_sub['rawT'][:,hour,plevel].data
    times = file_sub['time'][:]
    dt = pd.to_datetime( times, format='%Y%m%d' )
    df = pd.DataFrame( {"temp":temp, "time":dt, 'month': dt.month, 'year': dt.year} )
    # removing missing temperature values flagged with -9999.0
    cleaned_df = df.loc[df['temp'] > -999 ]
    return cleaned_df


def get_monthly_counts(df):
    ''' Extract the monthly counts 
                Parameters:
                    df (pandas DataFrame): data 
    '''
    years, months, times, counts = [],[],[],[]
    
    for y in range(1900,2019):
        df_y = df.loc[ df['year']==y ]
        
        for m in range(1,13):
            c = len (df_y.loc[ df_y['month'] == m ] )
    
            years.append(y)
            months.append(m)
            times.append(str(y)+str(m))
            counts.append(c)
            
            time_c = pd.to_datetime( times, format='%Y%m' )
            
    df_counts = pd.DataFrame(  { "year": years, "month": months, "time": time_c, "counts": counts } )    
    return df_counts

    
                     
                     
                     
                
                     
def plot(counts_z, counts_cuon, counts_ncar, counts_igra, plevel, hour = "00", station=''):
    fs = 15
    
    plt.step(counts_cuon["time"], counts_cuon["counts"], where='mid',
            label= 'CUON (this work) '   , color = 'blue')
    
    plt.step(counts_z["time"], counts_z["counts"], where='mid',
            label= 'Zhou at Al. 2020 '  , color = 'orange')

    plt.step(counts_ncar["time"], counts_ncar["counts"], where='mid',
            label= 'NCAR' , color = 'lime')
    
    plt.step(counts_igra["time"], counts_igra["counts"], where='mid',
            label= 'IGRAv2' , color = 'red')
    
    plt.grid(color = "lightgray", ls=':')
    plt.ylabel('Counts/month at ' + str(hour) + " GMT" , fontsize = fs)
    plt.legend(fontsize = fs)
    
    plt.title("Monthly temperature record - station " + station +  ' at ' + str(plevel) + ' [hPa]' , fontsize = fs )

    from datetime import date
    
    plt.xlim( [ date(1950, 1, 1)  , date(2020, 1, 1) ] )
    plt.show()
    
    
    
    
def get_CUON_data(file, var=85, merged=True, plev='', hour='' ):
    f = h5.File(file, 'r')
    
    # if merged, will use sorted recordindex 
    if merged:
        ind = f["recordindices"]['85']
        imin, imax = ind[0], ind[-1]

        dt = f["observations_table"]["date_time"][imin:imax]
        red_df = pd.DataFrame(  { 'z':f["observations_table"]["z_coordinate"][imin:imax] }  )
        dt = pd.to_datetime(dt,  unit='s', origin=pd.Timestamp('1900-01-01'))
        red_df["month"] = dt.month
        red_df["year"] = dt.year        
        red_df["hour"] = dt.hour        
            
        red_df["time"] = dt               

    else:
        dt = f["observations_table"]["date_time"]
        
        red_df = pd.DataFrame(  { 'z': f["observations_table"]["z_coordinate"], 
                                                    'observed_variable':  f["observations_table"]["observed_variable"],
                                                    'time': f["observations_table"]["date_time"] }   )
        
        
        red_df = red_df.loc[ red_df["observed_variable"] == 85 ]
        dt = red_df['time']
        
        dt = pd.to_datetime(dt,  unit='s', origin=pd.Timestamp('1900-01-01'))
        red_df["month"] = dt.dt.month
        red_df["year"] = dt.dt.year        
        red_df["hour"] = dt.dt.hour        
        
    red_df["time"] = dt        
        
        # Rounding to midnight/midday
    red_df = red_df.replace(   {'hour': {23: 0, 22: 0, 1:0, 2:0 , 10:12, 11:12, 13:12, 14:12 }}    )
        
    red_df = red_df.loc[ (red_df["z"] == plev) &   (red_df["hour"] == hour) ]
    red_df = red_df.dropna()
        
        
     
    return red_df



#station = "AUM00011035"
#a = get_station_index(f_HOMO, station) 


# CUON station
cuon_vienna = '/mnt/users/scratch/leo/scratch/converted_v7/0-20001-0-11035_CEUAS_merged_v1.nc'
#f_homo = "/users/staff/leo/fastscratch/SUNY/UA-HRD_stations_homo+raw-monthly-anomalyT_00-12Z_195801-202008.nc"

# HOMO station
subdaily_vienna = "/users/staff/leo/fastscratch/SUNY/homo-raw-subdaily-station/AUM00011035.nc"
subdaily_random = "/users/staff/leo/fastscratch/SUNY/homo-raw-subdaily-station/RSM00032215.nc"

#NCAR station
ncar_vienna = '/scratch/das/federico/MAY2021_HARVEST_secondary/ncar/0-20001-0-11035_ncar_harvested_uadb_trhc_11035.txt.nc' 

#IGRA station
igra_vienna = '/scratch/das/federico/MAY2021_HARVEST_secondary/igra2/0-20001-0-11035_igra2_harvested_AUM00011035-data.txt.nc' 

""" Get subdaily RAW data """
plevel = 11 # 500 hPa



"""
random = nc.Dataset(subdaily_random)

data_0_clean, data_0 = get_station_data(random, plevel, 0)
counts_0 = get_monthly_counts(data_0_clean)   

data_12_clean, data_12 = get_station_data(random, plevel, 1)
counts_12 = get_monthly_counts(data_12_clean)       
dummy = plot(counts_0, counts_12, 500)
"""






#Get HOMO data 
print('***** Extracting ZHOU data')
vienna_sub = nc.Dataset(subdaily_vienna)
data = get_station_data(vienna_sub, plevel, 0)
counts_sub = get_monthly_counts(data)              

#Get CUON data 
print('***** Extracting CUON data')
data = get_CUON_data(cuon_vienna, var=85, merged=True, plev=50000, hour = 0)
counts_cuon = get_monthly_counts(data )              


#Get NCAR data
print('***** Extracting NCAR data')
data = get_CUON_data(ncar_vienna, var=85, merged=False, plev=50000, hour = 0)
counts_ncar = get_monthly_counts(data )   


#Get IGRA data
print('***** Extracting IGRA data')
data = get_CUON_data(igra_vienna, var=85, merged=False, plev=50000, hour = 0)
counts_igra = get_monthly_counts(data )   

# Making plot
dummy = plot(counts_sub, counts_cuon, counts_ncar, counts_igra, 
                        500, 
                        hour = '00')



print(0)