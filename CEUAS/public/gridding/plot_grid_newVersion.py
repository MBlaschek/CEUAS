# fast gridded plot

import pandas as pd
import matplotlib.pylab as plt
import cartopy.crs as ccrs
import os,sys
import os
import json
import xarray as xr
import numpy as np
import shapely
from shapely.geometry import Point, Polygon

import matplotlib.animation as animation

import geopandas as gpd

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)
from tqdm import tqdm




def loading_WMO_regions_gpd():
    """ Getting the WMO regions json file """
    WMO_json = 'WMO_regions.json'
    if not os.path.isfile(WMO_json):
        os.system( 'wget https://cpdb.wmo.int/js/json/WMO_regions.json --no-check-certificate ')

    WMO =  gpd.read_file('WMO_regions.json')
    return WMO





def make_plot_gpd(WMO, sc = '', database = '', date = '', what = 'a', press = 85000):
    # https://stackoverflow.com/questions/59417997/how-to-plot-a-list-of-shapely-points                                                                                                                                                                                                   

    date_s = date
    date = np.datetime64(date)
    


    """ Loading from geopandas built-in methods """
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    #world = gpd.read_file('Europe_coastline.shp')
    #world = world.query( 'continent == "Europe"' )
    w = world.plot()
    WMO.plot( ax=w,  facecolor="none", edgecolor="lightgray", lw = 0.8)

    """ Saving coordinates to plot, and values """
    points_lat, points_lon, anomaly, average = [] , [] , [] , []
    lat_stations, lon_stations = [],[]
    none_lat, none_lon = [],[]
    
    """ Loop over files in database """
    db = os.listdir(database)    
    for f in tqdm(db):
        if 'dat' in f:
            no_lat, no_lon = f.split('_')[3] ,  f.split('_')[4]
            none_lat.append(float(no_lat))
            none_lon.append(float(no_lon))
            continue
        
        else:
            #loaded = xr.open_dataset(database + '/' + f).to_dataframe()
            #df_red = loaded.loc [ (loaded['date_time'] == date) & (loaded['z_coordinate'] == press ) 
            loaded = xr.open_dataset(database + '/' + f)
            stations = loaded['longitude'].attrs
            stat = stations['stations']
            if stat == 'None':
                continue
            else:
                for id in stat.split('_'):
                    ID = str(np.bytes_(id))
                    station = sc.loc [ sc['primary_id']== ID ]
                    try:
    
                        la = station['latitude'].values[0]
                        lo = station['longitude'].values[0]
                        if lo > 180:
                            lo = -360 + lo
                        
                        lat_stations.append(la)
                        lon_stations.append(lo)
                    except:
                        pass
            loaded = loaded.to_dataframe()
            df_red = loaded.loc [ (loaded['date_time'] == date) & (loaded['z_coordinate'] == press ) ]
            
            points_lat.append(df_red['latitude'].values)
            points_lon.append(df_red['longitude'].values)
            anomaly  .append(df_red['anomaly_bias'].values)
            average  .append(df_red['average_bias'].values)

        #print(points_lon, points_lat, anomaly, average)
   
    what = 'anomaly'   
    if what == 'anomaly':
        w = world.plot()
        WMO.plot( ax=w,  facecolor="none", edgecolor="lightgray", lw = 0.8)        
        plt.xlim([-180.,180.])
        plt.ylim([-90.,90.])        
        plt.scatter( lon_stations, lat_stations , c='lime',  s = 0.8, marker = 'o' )
        plt.scatter( none_lon, none_lat , c= 'lightgray' , s = 60, marker = 'x' )
        print(none_lat)
        plt.scatter( points_lon, points_lat , c= anomaly,  s = 60, marker = 's', cmap='bwr', alpha = 0.8, edgecolor = None)

        cbar = plt.colorbar(fraction=0.03, pad=0.03) # pad moves bar to left-right, fractions is the length of the bar        
        cbar.set_label('Temperature Anomaly over 20 years [K]')
        plt.clim(-5, 5)


        plt.title ('Climate Studies using Radiosonde Data - ' + date_s + ', p = ' +  str(press)[:-2] + ' [hPa] ' , fontsize = 9 )
        plt.savefig(out_dir + '/ClimateChange_' + date_s + '_' + what + '_' + '_plevel_' + str(press) + '.png', dpi= 250,   bbox_inches = 'tight' )
        plt.close()
        print('Done +++' , d  , '   ' , what )
        
        
    what = 'average'        
    
    if what == 'average':
        w = world.plot()
        WMO.plot( ax=w,  facecolor="none", edgecolor="lightgray", lw = 0.8)        
        plt.xlim([-180.,180.])
        plt.ylim([-90.,90.])        
        plt.scatter( lon_stations, lat_stations , c='black',  s = 1, marker = 'o')
        plt.scatter( points_lon, points_lat , c= average,  s = 60, marker = 's', cmap='rainbow' , alpha = 0.8, edgecolor = None )
        
        cbar = plt.colorbar(fraction=0.03, pad=0.03) # pad moves bar to left-right, fractions is the length of the bar        
        cbar.set_label('Average Temperature [K]')
        plt.clim(200, 330)


        plt.title ('Climate Studies using Radiosonde Data - ' + date_s + ', p = ' +  str(press)[:-2] + ' [hPa] ' , fontsize = 9)
        plt.savefig(out_dir + '/ClimateChange_' + date_s + '_' + what + '_' + '_plevel_' + str(press) + '.png', dpi= 250,   bbox_inches = 'tight' )
        plt.close()
        print('Done +++' , d  , '   ' , what )


""" Running """

out_dir = 'Plots_NEW/'
os.system('mkdir  ' + out_dir )
all_stat_conf = pd.read_csv('all_merged_station_configurations.csv', delimiter = '\t')



days = ['2019-06-15T12:00' , '2015-06-15T12:00' , '2010-06-15T12:00' , 
             '2019-12-15T12:00' , '2015-12-15T12:00' , '2010-12-15T12:00' , 
             '2000-06-15T12:00' ,  '2000-12-15T12:00' ]


days = ['2005-06-15T12:00' , '2005-12-15T12:00' ,  
             '1990-06-15T12:00' , '1990-12-15T12:00' , 
             '1980-06-15T12:00' ,  '1980-12-15T12:00' ,
             '1970-06-15T12:00' ,  '1970-12-15T12:00' ,
             '1960-06-15T12:00' ,  '1960-12-15T12:00' ,
             '1955-06-15T12:00' ,  '1955-12-15T12:00' ,
             '2000-06-15T12:00' ,  '2000-12-15T12:00'
             '2019-06-15T12:00' , '2019-12-15T12:00'
             ]


days = ['2005-06-15T12:00' , '2005-12-15T12:00' ,  
             '1980-06-15T12:00' ,  '1980-12-15T12:00' ,
             '1970-06-15T12:00' ,  '1970-12-15T12:00' ,
             '1960-06-15T12:00' ,  '1960-12-15T12:00' ,
             '2000-06-15T12:00' ,  '2000-12-15T12:00'
             '2019-06-15T12:00' , '2019-12-15T12:00'
             ]


WMO =  loading_WMO_regions_gpd()


plevels = [85000, 70000 ]
days = ['2000-06-15T12:00' , '1990-06-15T12:00' , '2010-06-15T12:00' , '1970-06-15T12:00'   ]
db = '/raid60/scratch/federico/GRIDDED_FILES_LATEST'

for d in days:
    for p in plevels:
        
        plotto = make_plot_gpd(WMO, sc = all_stat_conf, database = db , date = d , what = 'average' , press = p)        
        plotto = make_plot_gpd(WMO, sc = all_stat_conf, database = db , date = d , what = 'anomaly' , press = p )        


