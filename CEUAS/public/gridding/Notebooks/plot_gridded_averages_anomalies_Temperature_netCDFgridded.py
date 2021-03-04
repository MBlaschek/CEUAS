B""" Plotting utiity for netCDF gridded files """

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


import geopandas as gpd

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)
from tqdm import tqdm

from multiprocessing import Pool
from functools  import partial



""" Setting output directry """

out_dir = 'Plots/gridded/'
os.system('mkdir Plots')
os.system('mkdir  ' + out_dir )
def loading_WMO_regions_gpd():
    """ Getting the WMO regions json file """
    WMO_json = 'WMO_regions.json'
    if not os.path.isfile(WMO_json):
        os.system( 'wget https://cpdb.wmo.int/js/json/WMO_regions.json --no-check-certificate ')

    WMO =  gpd.read_file('WMO_regions.json')
    return WMO


# colorbar
import matplotlib.pylab as plt
import matplotlib.colors as mcolors
import numpy



def rgb(r,g,b):
    return tuple(numpy.asarray([r,g,b],dtype=numpy.float))

def make_colormap(seq):
        """Return a LinearSegmentedColormap
        seq: a sequence of floats and RGB-tuples. The floats should be increasing
        and in the interval (0,1).
        """
        se = [(None,) * 3, 0.0]
        for s in seq:
            se.append(s[0])
            se.append(s[1])#+ list(seq) +
            seq=se+[ (None,) * 3]
            cdict = {'red': [], 'green': [], 'blue': []}
        for i, item in enumerate(seq):
            if isinstance(item, float):
                r1, g1, b1 = seq[i - 1]
                r2, g2, b2 = seq[i + 1]
                cdict['red'].append([item, r1, r2])
                cdict['green'].append([item, g1, g2])
                cdict['blue'].append([item, b1, b2])
        return mcolors.LinearSegmentedColormap('CustomMap', cdict)

rgblist=["rgb(0,0,0.3)", "rgb(0,0,0.5)",
                 "rgb(0,0,0.7)", "rgb(0,0,0.9)", "rgb(0,0.15,1)",
                 "rgb(0,0.3,1)", "rgb(0,0.45,1)", "rgb(0,0.6,1)",
                 "rgb(0,0.75,1)", "rgb(0,0.85,1)", "rgb(0.2,0.95,1)",
                 "rgb(0.45,1,1)", "rgb(0.75,1,1)", "rgb(1,1,0)",
                 "rgb(1,0.9,0)", "rgb(1,0.8,0)", "rgb(1,0.7,0)",
                 "rgb(1,0.6,0)", "rgb(1,0.5,0)", "rgb(1,0.4,0)",
                 "rgb(1,0.3,0)", "rgb(1,0.15,0)", "rgb(0.9,0,0)",
                 "rgb(0.7,0,0)", "rgb(0.5,0,0)", "rgb(0.3,0,0)"]
rgblist2=zip([eval(rgblist[l]) for l in range(len(rgblist))],numpy.linspace(0,1,len(rgblist)))
                 
CM = make_colormap(rgblist2)


def read(file = '', var='ta'):
    
    df = xr.load_dataset(file, engine = 'h5netcdf').to_dataframe()
    df.reset_index(inplace=True)  
    return df
    
    


def make_plot_gpd(df, var, press, size, date):

    #def make_plot_gpd(WMO, sc = '', database = '', date = '', what = 'a', press = 85000, size = 10 ):
    # https://stackoverflow.com/questions/59417997/how-to-plot-a-list-of-shapely-points                                                                                                                                                                                                   

    date_s = date
    date = np.datetime64(date.split('T')[0])
    
    WMO =  loading_WMO_regions_gpd()

    """ Loading from geopandas built-in methods """
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    #world = gpd.read_file('Europe_coastline.shp')
    #world = world.query( 'continent == "Europe"' )
    w = world.plot()

    h = date_s.split('T')[1]
    if h[0]=='0':
        hour = 0
    else:
        hour = 12
    """ Select data from DF """
    reduced = df.loc[ (df['pressure']==press) & (df['time']==date) & (df['hour']==hour) ]

    points_lat = reduced['lat'].values
    points_lon = reduced['lon'].values
    anomaly = reduced[ var + '_anomaly'].values
    average = reduced[ var + '_average'].values
    


        
    def plot(what = ''):
        w = world.plot()
        WMO.plot( ax=w,  facecolor="none", edgecolor="lightgray", lw = 0.8)        
        plt.xlim([-180.,180.])
        plt.ylim([-95.,95.])   
        
        if size == 10:
            marker_size = 75
        elif size == 5:
            marker_size = 17
        
        x_t = list(range(-180,181,30))
        y_t = [ -90,-70,-50,-30,-10,0,10,30,50,70,90]
        y_t_l = ['90S','70S','50S','30S','10S','0','10N','30N','50N','70N','90N']
        x_t_l = ['180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E']
    
        plt.xticks(x_t, x_t_l, rotation=45, fontsize = 7)
        plt.yticks(y_t, y_t_l, fontsize = 7)
        
        if what == 'anomaly':
            plt.scatter( points_lon, points_lat , c= anomaly,  s = marker_size, marker = 's', cmap='bwr', alpha = 0.8, 
                     edgecolor = 'black' , linewidths=0.3 )
            cbar = plt.colorbar(fraction=0.03, pad=0.03) # pad moves bar to left-right, fractions is the length of the bar        
            cbar.set_label('Temperature Anomaly over 20 years [K]')
            plt.clim(-5, 5)
            
        elif what == 'average':
            plt.scatter( points_lon, points_lat , c= average,  s = marker_size, marker = 's', cmap= CM , alpha = 0.8,
                     edgecolor = 'black' , linewidths=0.3)
        
            cbar = plt.colorbar(fraction=0.03, pad=0.03) # pad moves bar to left-right, fractions is the length of the bar        
            cbar.set_label('Average Temperature [K]')
            plt.clim(200, 330)
        
        
        leg = plt.legend(fontsize = 7, loc = 'lower left', framealpha = 1, facecolor = 'lime')
        for lh in leg.legendHandles: 
            lh.set_alpha(1)

        plt.title ('Climate Studies using Radiosonde Data - ' + date_s + ', p=' +  str(press)[:-2] + ' [hPa] ' , fontsize = 8)
        plt.savefig(out_dir + '/ClimateChange_' + date_s + '_' + what + '_' + '_plevel_' + str(press) + '_gridsize_' + 
                    str(size) + '_netCDFgridded.png', dpi= 150,   bbox_inches = 'tight' )
        plt.show()
        plt.close()
        print('Done +++' , date , ' ' , what )
        
    dummy = plot(what='anomaly')
    dummy = plot(what='average')


# In[19]:


""" Running """



days = [ '1950-06-15T12:00' , '1950-12-15T12:00' , '1950-06-15T00:00' , '1950-12-15T00:00' ,
               '1960-06-15T12:00' , '1960-12-15T12:00' , '1960-06-15T00:00' , '1960-12-15T00:00' ,
               '1970-06-15T12:00' , '1970-12-15T12:00' , '1970-06-15T00:00' , '1970-12-15T00:00' ,
               '1980-06-15T12:00' , '1980-12-15T12:00' , '1980-06-15T00:00' , '1980-12-15T00:00' ,
               '1990-06-15T12:00' , '1990-12-15T12:00' , '1990-06-15T00:00' , '1990-12-15T00:00' ,
               '2000-06-15T12:00' , '2000-12-15T12:00' , '2000-06-15T00:00' , '2000-12-15T00:00' ,
               '2010-06-15T12:00' , '2010-12-15T12:00' , '2010-06-15T00:00' , '2010-12-15T00:00' ,
               
               
               '2019-06-15T12:00' , '2019-12-15T12:00' , '2019-06-15T00:00' , '2019-12-15T00:00' ,
               
               '2019-01-15T12:00' , '2019-02-15T12:00' , '2019-03-15T00:00' , '2019-04-15T00:00' ,
               '2019-05-15T12:00' , '2019-07-15T12:00' , '2019-08-15T00:00' , '2019-09-15T00:00' ,
               '2019-10-15T12:00' , '2019-11-15T12:00' ]


days = [ '1990-06-15T12:00' ,  '2019-06-15T12:00' ,  '2019-06-15T00:00' ]


#std_plevs    = [1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 85000, 92500, 100000]
plevels = [ 1000, 20000, 92500 ]


#out_dir = 'Plots/gridded/'
#os.system('mkdir Plots')
#os.system('mkdir  ' + out_dir )
all_stat_conf = pd.read_csv('../all_merged_station_configurations.csv', delimiter = '\t')

var = 'ta'
F = '/raid60/scratch/federico/GRIDDED_FILES_FEB2021/CEUAS_ta_gridded.nc/'
df = read(file = F, var='ta')


POOL = False
# POOL running
if POOL:
    p = Pool(40)
    for press in plevels:
            func = partial(make_plot_gpd, db_10 , var, press, 10)      
            out = p.map(func, days) 

else:     
    for d in days:   
        a = make_plot_gpd (df , var , 85000, 10, d) 


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




