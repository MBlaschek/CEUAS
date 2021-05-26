""" Analyze lat/lon of harvested files """

import subprocess
import pandas as pd
import os
import glob
import pandas as pd
import matplotlib.pylab as plt
import os,sys
import json
import xarray as xr
import numpy as np
import matplotlib.animation as animation
import geopandas as gpd
import h5py
from multiprocessing import Pool
from functools import partial

import numpy as np 


# In[6]:

"""
rpath='/raid60/scratch/leo/scratch/converted_v5/'
flist=glob.glob(rpath+'*0-99*CEUAS*.nc')
l=0
m=0
maxlist=[]
for fn in flist:
     with h5py.File(fn,'r') as f:
         lat=f['observations_table']['latitude'][:]
         lon=f['observations_table']['longitude'][:]
         lon[lon>180.]-=360
         ulat=np.unique(lat)
         ulon=np.unique(lon)
         l+=1
         if np.max(ulat)-np.min(ulat)>1 or np.max(ulon)-np.min(ulon)>1:
             m+=1
             print(l,m,fn,ulat,ulon)
maxlist.append((fn,np.float(np.max(ulat)),np.float(np.min(ulat)),np.float(np.max(ulon)),np.float(np.min(ulon))))
with open('laterr','w') as f:
     json.dump(maxlist,f)
"""


def run_fn(fn):
    res = { 'single' : { 'lat':[], 'lon':[], 'file':[] } ,   
                'multi'   : { 'lat':[], 'lon':[], 'file':[] } }
    try:
        data = h5py.File(fn,'r')
        lats = np.unique(data['observations_table']['latitude'][:])
        lons = np.unique(data['observations_table']['longitude'][:])
        
        
        if  len(lats) ==1 and len(lons) ==1: # one single pair of cordinates 
            res['single']['lat'].append( lats[0] )
            res['single']['lon'].append( lons[0] )
            res['single']['file'].append(fn)
        else:

            lat1, lat2 = min(lats), max(lats)
            lon1, lon2 = min(lons), max(lons)

            if (abs(lat1-lat2) < 0.5 and abs(lon1-lon2) < 0.5):
                res['single']['lat'].append(lat1)
                res['single']['lon'].append(lon1)
                res['single']['file'].append(fn)

            else:
                res['multi']['lat'].append([lat1,lat2])
                res['multi']['lon'].append([lon1,lon2])
                res['multi']['file'].append(fn)
                out = open(d+'_changing_coords_lat1lat2_lon1lon2.dat', 'a')
                lat1,lat2,lon1,lon2 = str(lat1), str(lat2), str(lon1), str(lon2)
                out.write(fn + '\t' + lat1 + '\t' + lat2 + '\t' + lon1 + '\t' + lon2 + '\n' )
        print('*** DONE ' , fn )
    except:
        print("*** Fail: ", fn )
        pass
    return res 



def plot(dataset, res):

    """ Getting the WMO regions json file """
    WMO_json = 'WMO_regions.json'
    if not os.path.isfile(WMO_json):
        os.system( 'wget https://cpdb.wmo.int/js/json/WMO_regions.json --no-check-certificate ')

    WMO =  gpd.read_file('WMO_regions.json')

    fs = 13
    plt.xlim([-180.,180.])
    plt.ylim([-90.,90.])
        #clb=f.colorbar(c,ax=ax,orientation='horizontal')                                                                                                                                                                                                                                     
        #clb.ax.set_title('{:6.4f}'.format(np.mean(ds[k]).values/43200)+' '+units)                                                                                                                                                                                                            

    """ Loading from geopandas built-in methods """
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    world = world.plot(figsize=(10,12))
    WMO.plot( ax=world,  facecolor="none", edgecolor="lightgray", lw = 0.5)

    plotto = plt.scatter( res['single']['lon'], res['single']['lat'] ,  
                          s = 3, color = 'red', 
                         label = 'Correct lat/lon '+ str(len(res['single']['lon'])) )
    #plotto = plt.scatter( res['single']['lon'], res['single']['lat'] ,  s = 0.7, color = 'red', label = 'Orphans ' )

    if 'lon' in res['multi'].keys(): # I have at least one 
        for lon,lat,num in zip(res['multi']['lon'],res['multi']['lat'], range(len(res['multi']['lat']))):
            if num == 1:
                plt.plot( [lon[0],lon[1]], [lat[0],lat[1]], 
                          color = 'lime', 
                         label = 'Changing lat/lon ' + str(len(res['multi']['lon'])))
                #plt.plot( [lon[0],lon[1]], [lat[0],lat[1]], color = 'blue' , label = 'Changing coords ' )

                plt.scatter( [lon[0],lon[1]], [lat[0],lat[1]], 
                             color = 'blue' ,s = 40, )

            else:
                plt.plot( [lon[0],lon[1]], [lat[0],lat[1]], color = 'lime')
                plt.scatter( [lon[0],lon[1]], [lat[0],lat[1]], 
                             color = 'blue', s = 40,)

    #cbar = plt.colorbar(fraction=0.03, pad=0.03) # pad moves bar to left-right, fractions is the length of the bar
    #cbar.set_label('Number of Records ')

    if dataset == 'ncar':
        ds = 'NCAR'
    elif dataset == 'bufr':
        ds = 'BUFR'
    elif dataset == 'igra2':
        ds = 'IGRA2'
    else:
        ds = 'ERA5 ' + dataset 
        
    tot = len(res['multi']['lon']) + len(res['single']['lon'])
    plt.title (str(tot) + ' Stations of the dataset ' + ds , fontsize = fs)

    plt.legend(fontsize = 12, loc = 'lower left')

    #plt.legend(loc = 'lower left', fontsize = 6 ,ncol = 2)
    plt.savefig('Plots/all_stations_map_' + dataset + '_harvested.png', dpi= 250,   bbox_inches = 'tight' )
    #plt.show()
        #plt.close()
    print('Done +++' , d )
    plt.close()





# looping over datasets 

multi = True
for d in ['bufr']:
#for d in ['era5_1759','era5_1761','era5_3188','era5_1','era5_2', 'igra2', 'ncar', 'bufr']:

     harvested_dir = '/raid60/scratch/federico/MAY2021_HARVEST/' + d
     flist = os.listdir( harvested_dir )
     flist = [harvested_dir +'/'+f for f in flist if  '.nc' in f and '00000' not in f and '9999' not in f ]
     print(len(flist) , ' *** Files ')
     tot = len(flist)
 
 
     if multi:
          
          p=Pool(25)
          func=partial(run_fn)
          all_res=list(p.map(func,flist))
 
     else:
          all_res = []
          for f in flist:
               r = run_fn(f)
               all_res.append(r)
               
     RES = { 'single' : { 'lat':[], 'lon':[], 'file':[] } ,   
                   'multi'  : { 'lat':[], 'lon':[], 'file':[] } }        #print(all_res)
 
     for r in all_res:
         if r['single']['lat']:
             RES['single']['lat'].append(r['single']['lat'][0])
             RES['single']['lon'].append(r['single']['lon'][0])
         if r['multi']['lat']:
             RES['multi']['lat'].append(r['multi']['lat'][0])
             RES['multi']['lon'].append(r['multi']['lon'][0])                
 
 
     np.save(d + '_data_lat_lon_from_harvested', RES )
     print(RES)
     plot(d,RES)



