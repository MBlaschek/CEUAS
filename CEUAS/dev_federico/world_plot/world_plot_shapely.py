""" Nice plotting utility to plot points given their
    geographical coordinates,
    e.g. to localize observing stations """


import netCDF4 as nc
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
                    
                    
    
    
def make_plot_gpd(WMO, df, start_date = '', end_date = '' ):
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
        
    style_dic = { 1 : { 'c' : 'navy' ,           'l': 'Africa'} , 
                          2 : { 'c' : 'lime' ,           'l': 'Asia' } ,
                          3 : { 'c' : 'red' ,            'l': 'South America' } ,
                          4 : { 'c' : 'gold' ,          'l': 'North America, Central America, Caribbean' } ,
                          5 : { 'c' : 'magenta' ,   'l': 'South-West Pacific' } ,
                          6 : { 'c' : 'cyan' ,         'l': 'Europe' }  }
    
    """
    style_dic = { 1 : { 'c' : 'blue' ,           'l': 'Africa'} , 
                          2 : { 'c' : 'lime' ,           'l': 'Asia' } ,
                          3 : { 'c' : 'red' ,            'l': 'South America' } ,
                          4 : { 'c' : 'gold' ,          'l': 'North America, Central America, Caribbean' } ,
                          5 : { 'c' : 'magenta' ,   'l': 'South-West Pacific' } ,
                          6 : { 'c' : 'cyan' ,         'l': 'Europe' }  }
    """

    df['start_date'] = pd.to_datetime(df['start_date'])
    df['end_date'] = pd.to_datetime(df['end_date'])
    LAT, LON, COL, LAB = [], [] , [] , []
    ims = []   

    for i in range(6):
        #print('Plotting WMO region: ' , i )
        a = df.loc[ (df['wmo'] == i ) & (df['start_date'] > start_date) & (df['end_date'] < end_date)   ]

        lat = list(a['latitude'].astype(np.float32) )
        lon = list(a['longitude'].astype(np.float32) )
        col = style_dic[i+1]['c'] 
        lab = style_dic[i+1]['l']

        if len(lat) == 0:
            lat.append(-999)
        if len(lon) == 0:            
            lon.append(-999)
  
        plotto = plt.scatter( LON, LAT , color = COL , label = LAB, s = 0.7 )
  
    """        
    for k in geoPoints.keys():
        xs = [point.x for point in  geoPoints[k] ]
        ys = [point.y for point in  geoPoints[k] ]        
        #plt.scatter( xs, ys , color = style_dic[k+1]['c'] , label = 'WMO region: ' +  style_dic[k+1]['l'], s = 1.5 )
        plt.scatter( xs, ys , color = style_dic[k+1]['c'] , label = style_dic[k+1]['l'], s = 1.6 )
    """
 
    plt.legend(loc = 'lower left', fontsize = 6 ,ncol = 2)
    plt.savefig('Map_shapely_' + start_date_string + '.png', dpi= 250,   bbox_inches = 'tight' )
    
    print(1)
        


 
def make_plot_anim(WMO, df, start_date = '', end_date = '' ):
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
        
    style_dic = { 1 : { 'c' : 'navy' ,           'l': 'Africa'} , 
                          2 : { 'c' : 'lime' ,           'l': 'Asia' } ,
                          3 : { 'c' : 'red' ,            'l': 'South America' } ,
                          4 : { 'c' : 'gold' ,          'l': 'North America, Central America, Caribbean' } ,
                          5 : { 'c' : 'magenta' ,   'l': 'South-West Pacific' } ,
                          6 : { 'c' : 'cyan' ,         'l': 'Europe' }  }
    
    """
    style_dic = { 1 : { 'c' : 'blue' ,           'l': 'Africa'} , 
                          2 : { 'c' : 'lime' ,           'l': 'Asia' } ,
                          3 : { 'c' : 'red' ,            'l': 'South America' } ,
                          4 : { 'c' : 'gold' ,          'l': 'North America, Central America, Caribbean' } ,
                          5 : { 'c' : 'magenta' ,   'l': 'South-West Pacific' } ,
                          6 : { 'c' : 'cyan' ,         'l': 'Europe' }  }
    """

    df['start_date'] = pd.to_datetime(df['start_date'])
    df['end_date'] = pd.to_datetime(df['end_date'])
    LAT, LON, COL, LAB = [], [] , [] , []
    ims = []   
    

    
    for i in range(6):
        #print('Plotting WMO region: ' , i )
        a = df.loc[ (df['wmo'] == i ) & (df['start_date'] > start_date) & (df['end_date'] < end_date)   ]

        lat = list(a['latitude'].astype(np.float32) )
        lon = list(a['longitude'].astype(np.float32) )
        col = style_dic[i+1]['c'] 
        lab = style_dic[i+1]['l']

        if len(lat) == 0:
            lat.append(-999)
        if len(lon) == 0:            
            lon.append(-999)
            
        
        
        ims.append( (plt.scatter( LON, LAT , color = COL , label = LAB, s = 0.7 ) ) , )
            
    print('Crated animation')
    fig2 = plt.figure()
    plt.show()
    im_ani = animation.ArtistAnimation(fig2, ims, interval=50, repeat_delay=3000, blit = False)
    im_ani.save(start_date_string + '_im.mp4', metadata={'artist':'Guido'} )
    
    

    
    
def prova_ani():
    fig2 = plt.figure()
    
    x = np.arange(-9, 10)
    y = np.arange(-9, 10).reshape(-1, 1)
    base = np.hypot(x, y)
    ims = []
    for add in np.arange(15):
        ims.append((plt.pcolor(x, y, base + add, norm=plt.Normalize(0, 30)),))
    
    im_ani = animation.ArtistAnimation(fig2, ims, interval=50, repeat_delay=3000,
                                       blit=False)
    # To save this second animation with some metadata, use the following command:
    im_ani.save('im.mp4')
    
    

'''
directory = '/raid60/scratch/federico/MERGED_20200130/'
def get_data_netCDF(directory):
    """ Extracting data from the station_configurations of the merged files """
    files = [ directory + f for f in os.listdir(directory) if '.nc' in f ][:300]

    data = {} 
    for d in ['era5_1' , 'era5_1759' , 'era5_1761' , 'era5_3188' , 'ncar_t' , 'ncar_w', 'igra2' , 'bufr']:
        data[d] = {}
        data[d]['latitude']   = []
        data[d]['longitude'] = []
        
    for f in files:    
        print(f)
        #for d in ['era5_1' , 'era5_1759' , 'era5_1761' , 'era5_3188' , 'ncar' , 'igra2' , 'bufr']:
        for d in ['era5_1','ncar_t', 'ncar_w']:
        
            try:
                sc = xr.open_dataset(f , engine = 'h5netcdf' , group = d + '_station_configuration')
                data[d]['latitude'].append( sc['latitude'].values[0] )
                data[d]['longitude'].append( sc['longitude'].values[0] )
                #sc.close()
            except:
                #print('failed')
                pass
        
    return data


if os.path.isfile('data.npy'):
    points = np.load('data.npy', allow_pickle = True).item()
else:
    points = get_data(directory)
'''



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
    
    
    
WMO = loading_WMO_regions_gpd()
if not os.path.isfile('Data.csv'):
    print('Reading the data from station_configuration')
    data = get_data_station_configuration(WMO)
else:
    print('Loading the data')
    data = pd.read_csv('Data.csv')
    

date = '19xx-01-01'
Start = ['1900-01-01' , '1910-01-01' ,'1920-01-01' , '1930-01-01' , '1940-01-01' ,  '1950-01-01'  , '1960-01-01' , '1970-01-01' , '1980-01-01' , '1990-01-01' ,'2000-01-01' , '2010-01-01' ]
End   = ['1910-01-01' ,'1920-01-01' , '1930-01-01' , '1940-01-01' , '1950-01-01' ,  '1960-01-01' , '1970-01-01'  , '1980-01-01' , '1990-01-01' ,'2000-01-01' , '2010-01-01' , '2020-01-01']
 
 
Start = ['1940-01-01' ,  '1950-01-01'  , '1960-01-01' , '1970-01-01' , '1980-01-01' , '1990-01-01' ,'2000-01-01' , '2010-01-01' ]
End   = ['1950-01-01' ,  '1960-01-01' , '1970-01-01'  , '1980-01-01' , '1990-01-01' ,'2000-01-01' , '2010-01-01' , '2020-01-01']
 
 
  
#plotto = make_plot_gpd(WMO, data , start_date = np.datetime64('1900-01-01'), end_date = np.datetime64('1910-01-01') )


prova_ani()
"""
for s , e in zip (Start, End):
    
    start = np.datetime64(s)
    end = np.datetime64(e)
    
    #plotto = make_plot_gpd(WMO, data , start_date = start, end_date = end )
    make_plot_anim(WMO, data , start_date = start, end_date = end )
"""
