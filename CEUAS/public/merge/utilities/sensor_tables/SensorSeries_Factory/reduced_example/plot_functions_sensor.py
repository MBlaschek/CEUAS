#!/usr/bin/env python
# coding: utf-8

# In[31]:


import os,sys
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import logging

import pandas as pd
import numpy as np

import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly as plotly
import plotly.io as pio


from multiprocessing import Pool
from functools import partial

import json


from sensor_functions import * 


class Plot():
    """ Main class to hold plotting utilities """
    
    def __init__(self,station, save=False):
        if not os.path.isdir('plots'):
            os.mkdir('plots')
            
            
        self.station=station
        self.save = save
        self.fontsize = 11
        
    def time_series(self, data_df, label = ''):
        """ Creates a time series using also SNHT """
        #filter date
        #data_df = data_df.loc[data_df.date_time <= pd.Timestamp('1995-01-01')]

        station = self.station.split('-')[-1] if '-' in self.station else self.station 
        
        if not os.path.isdir('plots'):
            os.mkdir('plots')
        if not os.path.isdir('plots/html'):
            os.mkdir('plots/html')
        if not os.path.isdir('plots/png'):
            os.mkdir('plots/png')    
            
        if not os.path.isdir('data_plots'):
            os.mkdir('data_plots')
            
        if not os.path.isfile('data_plots/' + self.station + '_snht.csv' ): 
            logging.debug(" +++ Retrieving the SNHT json file from the local server")
            try:
                with open('/mnt/users/staff/leo/python/CEUAS/CEUAS/public/adjust/feedbackmerged0' + station + '_breakanalysis.json') as f:
                    d=json.load(f)
                    time = pd.to_datetime(d['days_since_1900'] , unit='d', origin=pd.Timestamp('1900-01-01') )
                    snht = pd.DataFrame( {'time': time , 'snht':d['tsasum'] } )
                    snht.to_csv('data_plots/' + self.station + '_snht.csv' , sep ='\t')
                    logging.debug(" +++ Saving the SNHT csv file into the local data directory ")


            except:
                return 
            
        else:
            logging.debug(" +++ Reading the SNHT csv file from the local data directory")
            snht = pd.read_csv('data_plots/' + self.station + '_snht.csv' , sep = '\t')
            
        symbols = {"IGRA2":'star', "WMO":'circle', "SCH":'square'}

        # Create figure with secondary y-axis
        subfig = make_subplots(specs=[[{"secondary_y": True}]])
        fig1 = px.line(snht, x="time", y="snht")
        fig2 = px.scatter(data_df, 
                          x="date_time", 
                          y="value", 
                          color="sensor_id",
                          hover_name="sensor_id", 
                          hover_data=["comment"],
                          symbol="source",
                          symbol_map= symbols )

        fig2.update_traces(yaxis="y2")
        
        subfig.update_layout(legend = dict(font = dict(size = self.fontsize-2, 
                                                       color = "black"))
                            )
        
        
        subfig.add_traces(fig1.data + fig2.data)
        subfig.layout.xaxis.title=""
        subfig.layout.yaxis.title="SNHT"
        #subfig.layout.yaxis2.title="Metadata Source"
        subfig.layout.yaxis2.title=""
        

        subfig.for_each_trace(lambda t: t.update(line=dict(color=t.marker.color)))

        subfig.update_layout(title='Sensors Time Series - ' + self.station + ' ' + label)
        #subfig.update_layout(width= 2000, height = 800)


        subfig.update_traces(marker=dict(size=11,
                                      line=dict(width=1.5,
                                                color='DarkSlateGrey')),
                          selector=dict(mode='markers'))



        igra2 = data_df.loc[data_df.source == 'IGRA2']
        for d in igra2.date_time:
            subfig.add_vline(x=d, line_width=3, line_dash="dash", line_color="green")



        subfig.update_layout(hovermode="x unified")

        subfig.update_layout(
        yaxis = dict(
        tickfont = dict(size= self.fontsize)),
        font=dict(
            size= self.fontsize,
            color="black"
            )
        )

        subfig.update_yaxes( ticktext= ['Schr.', 'WMO', 'IGRA2'],
                          tickvals= [1,2,3], secondary_y=True )

        
        
        if self.save:
            plotly.offline.plot(subfig, filename=  'plots/html/' + self.station + "_time_series.html" )
            pio.write_image(subfig, "plots/png/" + self.station + "_timeSeries_ku.png")
            

        return subfig

    def sensor_table(self, data):
        
        if len(data)==0:
            return False
        
        fig = go.Figure(data=[go.Table(
        header=dict(values=list(['Sensor','Source', 'Comment']),
                    fill_color='gold',
                    align='left',
                    font_size=self.fontsize+1),
        columnwidth = [1,1,8],
            
        cells=dict(values=[data.sensor_id, data.source, data.comment],
                   fill_color='aliceblue',
                   align='left',
                   font_size=self.fontsize,
                   line_color='darkslategray',
                   fill=dict(color=['paleturquoise', 'white']),
                   height=50)
            
                  
        ) ] )

        #fig.update_layout(height= 65*len(data))

        if self.save:
            plotly.offline.plot(fig, filename=  'plots/html/' + self.station + "_sensor_table.html" )
            pio.write_image(fig, "plots/png/" + self.station + "_sensor_table.png")

        return fig
    
    
    def wmo_bar_plot(self, data):
        
        data = data.loc[ data['source'] == 'WMO' ]
        nans_ind = np.where( (data.sensor_id != 'nan') & (data.sensor_id != '-922') ) [0]
        data_df_clean = data.iloc[nans_ind]


        df = data_df_clean.groupby(["sensor_id"]).count().sort_values(by=['value'])
        fig = px.bar(df, y=df.index, x='value', 
                     title = 'WMO codes counts',
                     orientation='h')

        return fig 
   
