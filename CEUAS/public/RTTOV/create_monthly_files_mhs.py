import numpy
import numpy as np
import pandas
import pandas as pd
from numba import njit
import sys,glob
import zipfile, os, time
import urllib3
from datetime import datetime, timedelta
import glob
import h5py
import netCDF4 as nc
sys.path.append(os.getcwd()+'/../cds-backend/code/')
sys.path.append(os.getcwd()+'/../harvest/code/')
import rasotools
# from harvest_convert_to_netCDF_newfixes import write_dict_h5
import cds_eua4 as eua
# eua.logging_set_level(30)
import xarray as xr

import cdsapi, zipfile, os, time
#import schedule
import copy
from shutil import copyfile
import multiprocessing
import pickle

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.pylab as pylab
import warnings
warnings.filterwarnings('ignore')

from inspect import getmembers, isfunction


from IPython.display import Image
from IPython.core.display import HTML 
import rasotools


import matplotlib.pylab as plt
import matplotlib.colors as mcolors
import numpy

import ray
ray.init(num_cpus=40)

def rgb(r,g,b):
    return tuple(numpy.asarray([r,g,b],dtype=float))

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



# x,y=numpy.meshgrid(numpy.linspace(-3,3,101),numpy.linspace(-3,3,101))
# gauss=x/6+y/6
# clist=numpy.linspace(-1,1,26)
# plt.contourf(x,y,gauss,clist,cmap=cmnew)
# plt.colorbar()
# plt.show()
# print('ready')

def plt_trends(lla,pdict,var='_',bias='_', marker_size = 510, marker_shape = 's', alpha=0.8, cold_positive=False):


    rgblist=["rgb(0,0,0.3)", "rgb(0,0,0.5)",
    "rgb(0,0,0.7)", "rgb(0,0,0.9)", "rgb(0,0.15,1)",
    "rgb(0,0.3,1)", "rgb(0,0.45,1)", "rgb(0,0.6,1)",
    "rgb(0,0.75,1)", "rgb(0,0.85,1)", "rgb(0.2,0.95,1)",
    "rgb(0.45,1,1)", "rgb(0.75,1,1)", "rgb(1,1,0)",
    "rgb(1,0.9,0)", "rgb(1,0.8,0)", "rgb(1,0.7,0)",
    "rgb(1,0.6,0)", "rgb(1,0.5,0)", "rgb(1,0.4,0)",
    "rgb(1,0.3,0)", "rgb(1,0.15,0)", "rgb(0.9,0,0)",
    "rgb(0.7,0,0)", "rgb(0.5,0,0)", "rgb(0.3,0,0)"]

    if cold_positive:
        rgblist = np.flip(rgblist)
    rgblist2=zip([eval(rgblist[l]) for l in range(len(rgblist))],numpy.linspace(0,1,len(rgblist)))

    cmnew=make_colormap(rgblist2)

    params = {'legend.fontsize': 'x-large',
              'figure.figsize': (12, 8),
             'axes.labelsize': 'x-large',
             'axes.titlesize': 20,
             'xtick.labelsize':'medium',
             'ytick.labelsize':'medium'}
    pylab.rcParams.update(params)

    if 'scale' not in pdict.keys():
        pdict['scale']=2.0
    a = rasotools.plot._helpers.cost(lla[2],lla[1],lla[0])
    cost = np.sum(a)/len(a)

    ax = plt.axes([0, 0, 1, 1], projection=ccrs.PlateCarree())
    ax._autoscaleXon = False
    ax._autoscaleYon = False

    ax.add_feature(cfeature.OCEAN, zorder=0)
    ax.coastlines()

    plt.scatter(lla[2], lla[1], s=marker_size, alpha=alpha,
                c= lla[0],
                cmap=cmnew,
                vmin=-pdict['scale'],
                vmax=pdict['scale'],
                marker = marker_shape,
                edgecolor='k',)
    plt.colorbar(orientation='horizontal', label='Trend '+pdict['units'], shrink=0.9, pad=0.05)
    plt.tight_layout()
    plt.title('Brightness Temperature ' + str(pdict['start'])+'-'+str(pdict['stop'])+', '+str(pdict['pl'])+'\n'+'trend heterogeneity cost function: '+'{:.2f}'.format(cost)+'\n'+'number of stations: '+str(len(lla[0])))
        
    try:
        os.mkdir('plots_new')
    except:
        pass
#     plt.savefig('plots_new/cuon_'+names[bias]+'_'+str(pdict['start'])+'-'+str(pdict['stop'])+'_'+str(pdict['pl']), bbox_inches='tight')
    plt.show()
    plt.close()



pd.set_option('display.max_columns', None)
# pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
plt.rcParams['figure.figsize'] = [12, 8]

def plot_world_map(file, mission_channel, marker_size = 510, marker_shape = 's', alpha = 0.8, cold_positive=False):
    trends = pickle.load( open( file, "rb" ) )
    lats = []
    lons = []
    vals = []
    minl = []
    maxl =[]

    for i in trends:
        if trends[i] >= 15:
            maxl.append([i, trends[i]])
        elif trends[i] <= -15:
            minl.append([i, trends[i]])
        else:
            lats.append(float( i.split('_')[0]))#trends[i].lat))
            lons.append(float( i.split('_')[1]))#trends[i].lon))
            vals.append(float(trends[i].values))
    # print(' --- dropped ---')
    # print(maxl)
    # print(minl)
    # /10. for K/10a
    plt_trends(np.array([np.array(vals)/10., np.array(lats), np.array(lons)]), dict(var='temperature',pl=mission_channel,start='1992',stop='2022',units=r'K/10a'), marker_size = marker_size, marker_shape = marker_shape, alpha=alpha, cold_positive=cold_positive)
#     rasotools.plot.map.points(lon=np.array(lons), lat=np.array(lats), values=np.array(vals), vmin=-5, vmax=5)
    plt.show()
    plt.close()



### Convert amsub. data to gridded mean without interpolation

@ray.remote
def calc(df_mon, lat_range, lon_range, targetlat, targetlon):
    target_df = df_mon[np.logical_and(np.logical_and(df_mon.latitude >= lat_range[0], df_mon.latitude < lat_range[1]),np.logical_and(df_mon.longitude >= lon_range[0], df_mon.longitude < lon_range[1]))]
    r00 = np.nanmean(target_df['Ch3_BT'])
    r01 = np.nanmean(target_df['Ch4_BT'])
    r02 = np.nanmean(target_df['Ch5_BT'])
    r03 = len(target_df['Ch3_BT'])
    r04 = str(df_mon.iloc[0].Time.year) + '-' + str(df_mon.iloc[0].Time.month)
    return r00, r01, r02, r03, r04, targetlat, targetlon

######
for yr in [2010]: #range(2009,2018):
######
    lats = np.array(range(-8875,+9125, 250))/100.
    print(len(lats), lats)
    lons = np.array(range(125, 36000, 250))/100.
    print(len(lons), lons)

    time_series = {}
    for targetlon in lats:
        for targetlat in lons:
            time_series[str(targetlat) + '_' + str(targetlon)] = [[],[],[],[],[]]

    ######
    for imon in [1]: # range(1,13):
    ######
        print("Month: ", imon)
        fidu_files = glob.glob('/users/staff/uvoggenberger/scratch/fiduceo/dap.ceda.ac.uk/neodc/fiduceo/data/fcdr/microwave/v4.1/mhs/noaa19/'+str(yr)+'/'+str(imon).zfill(2)+'/*/*.0.1.nc')

        to_concat = []
        for fifi in fidu_files[:]:
            # print(fifi)
            with h5py.File(fifi) as fi:
                # print(fi.keys())
                vars = ['Ch3_BT', 'Ch4_BT', 'Ch5_BT', 'latitude', 'longitude', 'Time']
                to_df = {}
                for var in vars:
                    if var == 'Time':
                        to_df[var] = np.array([fi[var][:]]*90).flatten()
                    else:
                        to_df[var] = np.array(fi[var][:]).flatten()/100.
                    # print(var, len(to_df[var]))
            df = pd.DataFrame.from_dict(to_df)
    #         break
    #     break
    # break
            for chan in ['Ch3_BT', 'Ch4_BT', 'Ch5_BT']:
                df[chan].replace(655.35, np.nan, inplace=True)
            df.dropna(inplace=True)
            df.Time = pd.to_datetime(df.Time, unit='s')
            to_concat.append(df)

        df_mon = pd.concat(to_concat)
        ray_df_mon = ray.put(df_mon)
        result_ids = []
        for targetlon in lats:
            for targetlat in lons:
                lat_range = [targetlat - 1.25, targetlat + 1.25]
                lon_range = [targetlon - 1.25, targetlon + 1.25]
                result_ids.append(calc.remote(ray_df_mon, lat_range, lon_range, targetlat, targetlon))
        results = ray.get(result_ids)
        print()
        for ri in results:
            targetlat = ri[-2]
            targetlon = ri[-1]
            time_series[str(targetlat) + '_' + str(targetlon)][0].append(ri[0])
            time_series[str(targetlat) + '_' + str(targetlon)][1].append(ri[1])
            time_series[str(targetlat) + '_' + str(targetlon)][2].append(ri[2])
            time_series[str(targetlat) + '_' + str(targetlon)][3].append(ri[3])
            time_series[str(targetlat) + '_' + str(targetlon)][4].append(ri[4])
    
    pickle.dump( time_series, open( "/users/staff/uvoggenberger/scratch/RTTOV_output/mhs_data_"+str(yr)+".p", "wb" ) )

ray.shutdown()