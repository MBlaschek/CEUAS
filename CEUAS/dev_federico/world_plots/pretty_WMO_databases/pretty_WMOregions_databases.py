""" Nice plotting utility to plot points given their
    geographical coordinates,
    e.g. to localize observing stations """



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
import numpy as np
import cartopy.crs as ccrs


pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)



def loading_WMO_regions_gpd():   
    """ Getting the WMO regions json file """
    WMO_json = 'WMO_regions.json'
    if not os.path.isfile(WMO_json):
        os.system( 'wget https://cpdb.wmo.int/js/json/WMO_regions.json --no-check-certificate ') 
   
    WMO =  gpd.read_file('WMO_regions.json')
    return WMO

                    
    
os.system('mkdir Pretty_Plots ')

def make_plot_gpd(WMO, start_date = '', end_date = '' , databases = ''):
    # https://stackoverflow.com/questions/59417997/how-to-plot-a-list-of-shapely-points
    
    f=plt.figure(figsize=(10,5))
    ax=plt.subplot(1,1,1,projection=ccrs.PlateCarree())
        
    start_date_string , end_date_string = np.datetime_as_string(start_date) ,  np.datetime_as_string(end_date) 
    
    plt.xlim([-180.,180.])
    plt.ylim([-90.,90.])
    #clb=f.colorbar(c,ax=ax,orientation='horizontal')
    #clb.ax.set_title('{:6.4f}'.format(np.mean(ds[k]).values/43200)+' '+units)
    
    """ Loading from geopandas built-in methods """
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))

    world = world.plot()
    
    WMO.plot( ax=world,  facecolor="none", edgecolor="lightgray", lw = 0.8)
    
    plt.title ('Active Radiosondes between ' + start_date_string + ' and ' + end_date_string , fontsize = 9)
        
    style_dic = { 'era5_1'    : { 'c' : 'navy' ,    'l': 'ERA5 1'  , 's': 4 } , 
                  'era5_2'    : { 'c' : 'lime' ,    'l': 'ERA5 2'  , 's': 4 } ,
                  'era5_3188' : { 'c' : 'black' ,   'l': '3188'    , 's': 9 } ,
                  'era5_1759' : { 'c' : 'magenta' , 'l': '1759'    , 's': 9 } ,
                  'era5_1761' : { 'c' : 'red' ,     'l': '1761'    , 's': 6 } ,
                  'bufr'      : { 'c' : 'orange' ,  'l': 'BUFR'    , 's': 6 } ,
                  'igra2'     : { 'c' : 'yellow' ,  'l': 'IGRA2'   , 's': 4 } ,
                  'ncar'      : { 'c' : 'cyan' ,    'l': 'NCAR'    , 's': 4 } ,
}

    for d in databases:
        print(' *** processing the database::: ', d )
        df = pd.read_csv('summaries/' + d + '_summary_distribution.dat', delimiter = '\t')
        LAT, LON = [], []
        for index, row in df.iterrows():
            start, end = np.datetime64(row['start']), np.datetime64(row['end'])
            #print (start, end, start_date, end_date )
            try:
                lat, lon = int(row['lat']) , int(row['lon'])
            except:
                pass

            if (start <= start_date) and (end >= start_date): # checking if the station is alive in the period
                LAT.append(lat)
                LON.append(lon)

        #lat = list(a['latitude'].astype(np.float32) )
        #lon = list(a['longitude'].astype(np.float32) )
        col = style_dic[d]['c'] 
        lab = style_dic[d]['l']
    
        if len(LAT) == 0:
            LAT.append(-999)
        if len(LON) == 0:            
            LON.append(-999) 
            plotto = plt.scatter( LON, LAT , color = col , label = lab + '[0]', s = 0.7 )
        else:
            #print(LON, '\n', LAT )
            plotto = plt.scatter( LON, LAT , color = col , label = lab + ' [' + str(len(LAT)) + ']', s = 0.7 )
 

    plt.xlim([-180.,180.])
    plt.ylim([-90.,90.])

    plt.legend(loc = 'lower left', fontsize = 6 ,ncol = 2)
    plt.savefig('Pretty_Plots/PROVA_' + start_date_string + '.png', dpi= 250,   bbox_inches = 'tight' )
    
    print(1, ' PLOT MADE !!!' )
        






""" Utility to read directly the station_configuration files """
def get_data_station_configuration(WMO, file= 'station_configuration_igra2.dat'):
        
    df = pd.read_csv(file ,  delimiter='\t', quoting=3, na_filter=False, comment='#')
    
    data = { 'latitude': [] , 'longitude':[], 'names':[], 'wmo':[] , 'start_date':[], 'end_date': [] , 'wmo_name':[]  }
    
    for index,r in df.iterrows():        
        try:
            
            pid      = r.primary_id
            name  = r.station_name
            lat       = float(r.latitude)
            lon      = float(r.longitude)
            min_t  = np.datetime64( r.start_date )
            max_t = np.datetime64( r.end_date )
            
            wmo, wmo_name = get_WMO_region(WMO, lat, lon) 
                
            data['latitude'].append(lat)
            data['longitude'].append(lon)
            data['names'].append(name)
            data['wmo'].append(int(wmo))
            data['start_date'].append(min_t)
            data['end_date'].append(max_t)
            data['wmo_name'].append(wmo_name)
                        
        except:
            print('Failed getting the wmo information')
            pass
        
    DF = pd.DataFrame.from_dict(data)
    DF.to_csv('Data.csv')
    return DF
    
    




date = '19xx-01-01'
Start = ['1900-01-01' , '1910-01-01' ,'1920-01-01' , '1930-01-01' , '1940-01-01' ,  '1950-01-01'  , '1960-01-01' , '1970-01-01' , '1980-01-01' , '1990-01-01' ,'2000-01-01' , '2010-01-01' ]
End   = ['1910-01-01' ,'1920-01-01' , '1930-01-01' , '1940-01-01' , '1950-01-01' ,  '1960-01-01' , '1970-01-01'  , '1980-01-01' , '1990-01-01' ,'2000-01-01' , '2010-01-01' , '2019-01-01']
 
 
Start = ['1950-01-01'  , '1960-01-01' , '1970-01-01' , '1980-01-01' , '1990-01-01' ,'2000-01-01' , '2010-01-01' ]
End   = ['1960-01-01' , '1970-01-01'  , '1980-01-01' , '1990-01-01' ,'2000-01-01' , '2010-01-01' , '2019-01-01']
 
 
#Start = ['1950-01-01'  , '1960-01-01' , '1970-01-01' , '1980-01-01' , '1990-01-01']
#End   = ['1960-01-01' , '1970-01-01'  , '1980-01-01' , '1990-01-01' ,'2000-01-01']


databases = ['era5_1','era5_2','era5_1759','era5_1761','era5_3188','bufr','ncar', 'igra2' ]

WMO = loading_WMO_regions_gpd()
for s , e in zip (Start, End):
    
    start = np.datetime64(s)
    end = np.datetime64(e)
    
    plotto = make_plot_gpd(WMO, start_date = start, end_date = end , databases = databases)

