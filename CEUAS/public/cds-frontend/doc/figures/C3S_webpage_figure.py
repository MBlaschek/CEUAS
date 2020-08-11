#!/usr/bin/env python
# coding: utf-8

import os
import numpy
import matplotlib
from matplotlib import pylab, mlab, pyplot
np = numpy
plt = pyplot
from IPython.core.pylabtools import figsize, getfigs
from pylab import *
from numpy import *

import cartopy as cpy
from matplotlib.colors import BoundaryNorm
import pandas as pd


idir = '../../../harvest/data/station_configurations/'


# In[8]:


files = os.listdir(idir)
print("\n".join(files))


# In[10]:


files_ordered = {}
for ifile in files:
    idata = pd.read_csv(idir + ifile,  delimiter='\t', quoting=3, comment='#')
    print(ifile, idata.shape[0], idata['observed_variables'].map(len).unique())
    files_ordered[ifile] = idata.shape[0]


# In[11]:


idata['primary_id']


# In[12]:


files_ordered = list(pd.Series(files_ordered).sort_values().index.values[::-1])


# In[13]:


print(files_ordered)


# In[14]:


# In[15]:


plt.rcParams['lines.linewidth'] = 2
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.major.size'] = 10
plt.rcParams['xtick.minor.size'] = 10
plt.rcParams['font.size'] = 15
plt.rcParams['figure.figsize'] = (19,10.6)
plt.style.use('seaborn-colorblind')


# ## Sources Map


projection = cpy.crs.PlateCarree()
ax = plt.axes(projection=projection)

ocean_facecolor='w'  # cpy.feature.COLORS['water']
land_facecolor='lightgray' # cpy.feature.COLORS['land']
if True:
    ax.add_feature(cpy.feature.OCEAN, zorder=0, facecolor=ocean_facecolor)

if True:
    ax.add_feature(cpy.feature.LAND, zorder=0, facecolor=land_facecolor)

if True:
    ax.coastlines()

if False:
    ax.add_feature(cpy.feature.LAKES, zorder=0)
    ax.add_feature(cpy.feature.RIVERS, zorder=1)

unique_stations = []
data_src = {}
# Sources:
for ifile in files_ordered:
    idata = pd.read_csv(idir + ifile,  delimiter='\t', quoting=3, comment='#')
    lon = idata['longitude'].values
    lat = idata['latitude'].values
    if any(lon > 180):
        lon = np.where(lon > 180, lon - 360, lon)
    idx = np.isfinite(lon) & np.isfinite(lat) & ((lon >= -180) & (lon <= 180)) & ((lat >= -90) & (lat <= 90)) & (idata['primary_id'].str.contains('0'))
    idata = idata[idx].drop_duplicates(subset='primary_id', keep='last')
    lon = idata['longitude'].values
    if any(lon > 180):
        lon = np.where(lon > 180, lon - 360, lon)
    lat = idata['latitude'].values
    ilabel = ifile.replace('station_configuration_','').replace('.dat','').upper().replace('_',' ') #+ ' (# %d)' % idx.sum()
    if '3188' in ilabel:
        ilabel = 'CHUAN'
    if 'BUFR' in ilabel:
        ilabel = 'ERA40'
    if 'ERA5' in ilabel:
        print(ilabel, 'to ECMWF')
        ilabel = 'ECMWF'
    
    
    if ilabel in data_src.keys():
        data_src[ilabel] = np.unique(data_src[ilabel] + list(zip(lon,lat)),axis=1).tolist()
    else:
        data_src[ilabel] = list(zip(lon,lat))
    
for ikey in data_src.keys():
    tmp = np.array(data_src[ikey])
    ilon, ilat = tmp[:,0], tmp[:,1]
    print(ikey, ilon.size)
    ax.scatter(ilon,ilat, s=40 , transform=cpy.crs.PlateCarree(), edgecolor='k', 
               label=ikey + ' (# %d)' % ilon.size)  # ontop
    unique_stations = list(np.unique(unique_stations + ["{:.1f}{:.1f}".format(i, j) for i,j in zip(ilon,ilat)]))

if True:
    try:
        gl = ax.gridlines(draw_labels=True, xlocs=None, ylocs=None,
                          linewidth=0.5, linestyle='--', color='k')
        gl.xformatter = cpy.mpl.gridliner.LONGITUDE_FORMATTER
        gl.yformatter = cpy.mpl.gridliner.LATITUDE_FORMATTER
        gl.xlabels_top = False
        gl.xlabels_bottom = False
        gl.ylabels_left = False
    except:
        ax.gridlines(draw_labels=False)

title = 'Locations of Ballon Records since 1905'
ax.set_title(title)
ax.legend(bbox_to_anchor=(0.5,-0.15), loc='lower center', ncol=4)
# savefig('C3S_webpage_logo.pdf', dpi=300, bbox_inches = 'tight',  pad_inches = 0)
savefig('CEUAS_network_v2.png', dpi=150, bbox_inches = 'tight',  pad_inches = 0)




