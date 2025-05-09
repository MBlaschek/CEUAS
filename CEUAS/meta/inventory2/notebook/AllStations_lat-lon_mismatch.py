
""" Check if lat/lon change significantly in the files, causing the non-true-orphans problem """
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
    
from multiprocessing import Pool
from functools import partial

import numpy as np 


# In[6]:


def run_fn(fn):
        res = { 'single' : { 'lat':[], 'lon':[], 'file':[] } ,   
                    'multi'   : { 'lat':[], 'lon':[], 'file':[] } }
        qs='select distinct lat, lon'  #distinct,count(date) 

        try:
            rdata=subprocess.check_output(["odb","sql","-q",qs,"-i",fn,'--no_alignment'],stderr=subprocess.DEVNULL)
            data = rdata.decode('utf-8').split('\n')
            #print(rdata)
            #print(data)
            #print(data[1:-1])

            chunk = data[1:-1]

            #print(chunk)

            if  len(chunk) ==1:
                res['single']['lat'].append(float(chunk[0].split('\t')[0]))
                res['single']['lon'].append(float(chunk[0].split('\t')[1]))
                res['single']['file'].append(fn)
            else:

                lat1, lat2 = float(chunk[0].split('\t')[0]) , float(chunk[-1].split('\t')[0])
                lon1, lon2 = float(chunk[0].split('\t')[1]) , float(chunk[-1].split('\t')[1])

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
    
    tot = '0'
        
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

    if dataset == 'rda':
        ds = 'NCAR'
    elif dataset == 'ai_bfr':
        ds = 'BUFR'
    else:
        ds = 'ERA5 ' + dataset 
        
    tot = len(res['multi']['lon']) + len(res['single']['lon'])
    plt.title (str(tot) + ' Stations of the dataset ' + ds , fontsize = fs)
    
    plt.legend(fontsize = 12, loc = 'lower left')
    
    #plt.legend(loc = 'lower left', fontsize = 6 ,ncol = 2)
    plt.savefig('Plots/all_stations_map_' + dataset + '.png', dpi= 250,   bbox_inches = 'tight' )
    #plt.show()
        #plt.close()
    print('Done +++' , d )
    plt.close()





# looping over datasets 
odb_dir = '/raid60/scratch/leo/scratch/era5/odbs/'

datasets = { '1759' : 'era5.1759.conv.*'    ,
                     '1761' : 'era5.1761.conv.*'    ,
                     '1'    : 'era5.conv._*'               ,
                     '2'    : 'era5.conv._*'               ,
                     '3188' : 'era5.3188.conv.*'    }

for d in ['1759','1761','3188','1','2']:

    flist = glob.glob(odb_dir+'/'+ d + '/' + datasets[d])
    flist = [f for f in flist if '.gz' not in f and '.nc' not in f and '0000' not in f and '9999' not in f ]
    #print(flist)
    print(len(flist) , ' *** Files ')
    tot = len(flist)

    p=Pool(25)
    func=partial(run_fn)
    all_res=list(p.map(func,flist))

    RES = { 'single' : { 'lat':[], 'lon':[], 'file':[] } ,   
                  'multi'  : { 'lat':[], 'lon':[], 'file':[] } }        #print(all_res)
    
    for r in all_res:
        if r['single']['lat']:
            RES['single']['lat'].append(r['single']['lat'][0])
            RES['single']['lon'].append(r['single']['lon'][0])
        if r['multi']['lat']:
            RES['multi']['lat'].append(r['multi']['lat'][0])
            RES['multi']['lon'].append(r['multi']['lon'][0])                
        print(0)
            

    np.save(d + '_data_lat_lon', RES )
    print(RES)
    plot(d,RES)
