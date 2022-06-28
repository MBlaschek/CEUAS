#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
sys.path.append('modules')

from sensor_functions import *
from plot_functions_sensor import *

import glob 

import logging
logging.basicConfig()
logging.getLogger().setLevel(logging.DEBUG)

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)



# Select Station
station = '0-20000-0-82930'
merged = '/scratch/das/federico/CHACH_sensor'
merged = '/scratch/das/federico/CHACH'

clean_all = True 
if clean_all :
    files = glob.glob('data/*'+station + '*') + glob.glob('data_o/*'+station + '*')
    for f in files:
        os.remove(f)
      
        
fc = True



data_clean_all, all_sensor_station, data_all_wmo = get_data(merged, station, force_create=fc)





data_clean_all


save_fig = False
plot = Plot(station.split('_')[-1], save=save_fig)




# making the plots
series = plot.time_series( data_clean_all, label='')
table = plot.sensor_table( all_sensor_station)
wmo_bar = plot.wmo_bar_plot(data_all_wmo)


# In[7]:


#data_all_wmo


# In[8]:


#w = data_clean_all[data_clean_all.source=="WMO"]


# In[9]:


#data_clean_all


# In[10]:


#table = plot.sensor_table( all_sensor_station)
#table.show()


# In[11]:


series.show()


# In[12]:


table = plot.sensor_table( all_sensor_station)
table.show()


# In[13]:


wmo_bar.show()

