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



db = '/raid8/srvx1/federico/GitHub/CEUAS_master_NOVEMBER2020/CEUAS/CEUAS/public/gridding/output_monthly'


def loading_WMO_regions_gpd():
    """ Getting the WMO regions json file """
    WMO_json = 'WMO_regions.json'
    if not os.path.isfile(WMO_json):
        os.system( 'wget https://cpdb.wmo.int/js/json/WMO_regions.json --no-check-certificate ')

    WMO =  gpd.read_file('WMO_regions.json')
    return WMO


def get_WMO_region(WMO, lat, lon):
    """ Extract the WMO region code given lat and lon """
    geoPoint = Point(lon,lat)

    for wmo_code,row in WMO.iterrows():
            geom = row.geometry
            for g in geom:
                if geoPoint.within(g):
                    #print('Point ' , Point , ' is in region: ', wmo_code , ' ' , row.Region_en  )                                                                                                                                                                                        
                    return wmo_code , row.Region_en





def make_plot_gpd(WMO, database = '', date = '', what = 'a'):
    # https://stackoverflow.com/questions/59417997/how-to-plot-a-list-of-shapely-points                                                                                                                                                                                                   

    plt.xlim([-180.,180.])
    plt.ylim([-90.,90.])

    """ Loading from geopandas built-in methods """
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    world = world.plot()

    """ Loading from geopandas built-in methods """
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    world = world.plot()
    WMO.plot( ax=world,  facecolor="none", edgecolor="lightgray", lw = 0.8)

    points_lat, points_lon, anomaly, average = [] , [] , [] , []
        
    for f in os.listdir(database):
        print (f)
        loaded = xr.open_dataset(database + '/' + f).to_dataframe()
        df_red = loaded.loc [ (loaded['date_time'] == date) & (loaded['z_coordinate'] == 85000 ) ]

        points_lat.append(df_red['latitude'].values)
        points_lon.append(df_red['longitude'].values)
        anomaly.append(df_red['anomaly_bias'].values)
        average.append(df_red['average_bias'].values)

        #print(points_lon, points_lat, anomaly, average)
        
    
    if what == 'a':
        plotto = plt.scatter( points_lon, points_lat , c= anomaly,  s = 50, marker = 's', cmap='rainbow', alpha = 0.7 )
        cbar = plt.colorbar(fraction=0.03, pad=0.03) # pad moves bar to left-right, fractions is the length of the bar        
        cbar.set_label('Temperature Anomaly over 20 years [K]')
        plt.clim(-5, 5)

    else:        
        plotto = plt.scatter( points_lon, points_lat , c= average,  s = 50, marker = 's', cmap='rainbow' , alpha = 0.7 )
        cbar = plt.colorbar(fraction=0.03, pad=0.03) # pad moves bar to left-right, fractions is the length of the bar        
        cbar.set_label('Average Temperature [K]')
        plt.clim(250, 330)
    

    date_string = 'June 2019'
    plt.title ('Climate Studies using Radiosonde Data ' + date_string , fontsize = 9)

    plt.xlim([-180.,180.])
    plt.ylim([-90.,90.])
    plt.savefig(out_dir + '/ClimateChange_' + date_string + '_' + what + '_average.png', dpi= 250,   bbox_inches = 'tight' )
    #plt.show()
    plt.close()
    print('Done +++' , d )


""" Running """


d = np.datetime64('2019-06-15T12:00:00.000000000')

out_dir = 'Plots/'
os.system('mkdir  ' + out_dir )
    
WMO =  loading_WMO_regions_gpd()
plotto = make_plot_gpd(WMO, database = db , date = d , what = 'no' )        
plotto = make_plot_gpd(WMO, database = db , date = d , what = 'a')        
