""" Create a single netCDF file contaning gridded data.
      structure of file:
      - dimensions
          lat
          lon
          time
          pleve
"""

import os,sys

import xarray as xr
import numpy as np
import pandas as pd



def initialize(variable = '' ):
    
    global Lat, Lon, Plev, Hour, Time, files, gridded_files_dir, hours, var, out_dir
    
    var = variable
    out_dir = '/raid60/scratch/federico/GRIDDED_FILES_FEB2021/'
    
    """ Directory with gridded files, list of files """
    gridded_files_dir = '/raid60/scratch/federico/GRIDDED_FILES_FEB2021/' + var + '/'
    files =  os.listdir(gridded_files_dir)
    
    """ Creating the latitude and longitude lists """
    lat = list(set( [ float(f.split('_')[5]) for f in os.listdir(gridded_files_dir) ] ) ) 
    lat.sort()
    lon = list(set( [ float(f.split('_')[6]) for f in os.listdir(gridded_files_dir) ] ) ) 
    lon.sort()
    
    
    Lat = np.array(lat)
    Lon = np.array(lon)
    Plev = np.array([1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 85000, 92500, 100000])
    Hour = np.array([0,12])
    
    """ Creating hour array (easier to find time in the df ) """
    df = xr.open_dataset( gridded_files_dir + '/' + files[0] , engine = 'h5netcdf', decode_times = True ).to_dataframe()
    q = pd.arrays.DatetimeArray( df['time'][:] )
    hours = q.hour
    df['hour'] = hours
    
    """ Total length of data """
    Time = df.loc [ (df['hour'] == 12) & (df['plev']==100000) ]['time'].values  # just need one plev per  hour, i.e. this is the list of distinct  time stamps



#Lat = Lat[:3]
#Lon = Lon[:5]

def make_xarray():
    res_average  = np.empty([len(Lat) , len(Lon), len(Hour) , len(Plev),  len(Time) ] )  # 2*16 is 2 hours x pressure levels , empty placeholder 
    res_anomaly = np.empty([len(Lat) , len(Lon), len(Hour) , len(Plev),  len(Time) ] )  # 2*16 is 2 hours x pressure levels , empty placeholder 
    
    #for lat in range(len(Lat)):
    #    for lon in range(len(Lon)):
    for lat in Lat:
        for lon in Lon:
            box_file = [f for f in files if '_' + str(lat) +'_'+str(lon) in f ][0]  # search for the grid file for those lat and lon 
            print(lat , ' ' , lon , ' ' , box_file )
            
            #if '24_20_lat' in box_file:
            #    print('check' , box_file )
            #else:
            #    continue
            #    print(0)
            df = xr.open_dataset( gridded_files_dir + '/' + box_file , engine = 'h5netcdf', decode_times = True ).to_dataframe()
            #q = pd.arrays.DatetimeArray( df['time'][:] )
            #hh = q.hour
            df['hour'] = hours  # it is the same for each dataframe     

            for p in range(len(Plev)):
                a = 0 # read here input data 
                for h in range(len(Hour)):
                    print('processing lat lon press hour : ',  lat,lon,p,h)
                    
                    df_red = df.loc[ (df['plev'] == Plev[p] )   & ( df['hour'] == Hour[h] )]
                    if var == 'ta':
                        values = df_red[var + '_average_bias'].values
                        ano = df_red[var + '_anomaly_bias'].values
                    else:
                        values = df_red[var + '_average'].values
                        ano = df_red[var + '_anomaly'].values                        

                    
                    res_average[np.where(Lat == lat),np.where(Lon == lon),h,p,:] = values
                    res_anomaly[np.where(Lat == lat),np.where(Lon == lon),h,p,:] = ano
                    
    return res_average, res_anomaly



def write_gridded_netCDF(res_average, res_anomaly):
    """ Writing the output netCDF files using xarray """
    da = xr.Dataset ( { var + '_average' : ( ['lat','lon','hour','pressure','time'] ,  res_average ), # variables 
                                     var + '_anomaly' : ( ['lat','lon','hour','pressure','time'] ,  res_anomaly ) },
                      
                                    coords = dict(   lat = Lat ,
                                               lon =  Lon,
                                               hour = Hour,
                                               pressure = Plev,
                                               time = Time,                  
                                               ),
                         )
    
    da.attrs['title'] = 'CEUAS gridded data'
    da.attrs['institution'] = 'University of Vienna',
    da.attrs['source']  = 'Institut fuer Meteorologie und Geophysik, leopold.haimberger@univie.ac.at'
    
    
    attr_dic = { 'ta' :                {'variable': 'Air temperature' , 'units': 'Kelvin [K]' },
                        'wind_speed' : {'variable': 'Wind speed' , 'units': 'meter per second [m/s]' },
                        'dew_point_temperature' : {'variable': 'Dew point temperature' , 'units': 'Kelvin [K]' },
                        'relative_humidity' : {'variable': 'relative humidity' , 'units': 'unitless [0-1]' } }
                        

    da.attrs['variable'] = attr_dic[var]['variable']
    da.attrs['units'] = attr_dic[var]['units']
    
    da.lat.attrs['variable'] = 'Latitude'
    da.lat.attrs['units'] = 'Degrees'
    
    da.lon.attrs['variable'] = 'Longitude'
    da.lon.attrs['units'] = 'Degrees'
    
    da.time.attrs['variable'] = 'Timestamp'
    #da.time.attrs['units'] = 'Degrees'
    
    da.hour.attrs['variable'] = 'Hour'
    da.hour.attrs['units'] = 'Hour from 00'
    
    da.pressure.attrs['variable'] = 'Pressure'
    da.pressure.attrs['units'] = 'Pascal [Pa]'
    
    from datetime import datetime
    
    now = datetime.now()
    dt = now.strftime("%d/%m/%Y %H:%M:%S")# dd/mm/YY H:M:S
    da.attrs['history'] = dt                    
    
                                                                 
    da.to_netcdf( out_dir + '/CEUAS_' + var +'_gridded.nc' )
    print('**Written file: ' ,   out_dir + '/CEUAS_' + var +'_gridded.nc *** '  )

""""
da = xr.Dataset ( { var + '_average' : ( ['lat','lon','hour','pressure','time'] ,  res_average ), # variables 
                                 var + '_anomaly' : ( ['lat','lon','hour','pressure','time'] ,  res_anomaly ) },
                  
                  #dims = ["lat","lon","hour","pressure","time"],
                  coords = dict(   lat = Lat ,
                                           lon =  Lon,
                                           hour = Hour,
                                           pressure = Plev,
                                           time = Time,                  
                                           ),
                     )
                                                             



da = xr.DataArray (data = res, name = 'ta',
                                dims = ["lat","lon","hour","pressure","time"],
                                coords = dict(   lat = Lat ,
                                                         lon =  Lon,
                                                         hour = Hour,
                                                         pressure = Plev,
                                                         time = Time,                  
                                                         ),
                                   )
"""                                                                                                       
                                           


    
    
""" Select here the variable to process.
      Will locate the proper directory and file list """

var = 'ta'
#Lat = Lat[:3]
#Lon = Lon[:5]
if __name__ == '__main__':

    dummy = initialize(variable = var)
    #Lat = Lat[10:]  ### uncomment for testing 
    #Lon = Lon[15:]    
    res_average, res_anomaly = make_xarray()
    write_gridded_netCDF(res_average, res_anomaly )





