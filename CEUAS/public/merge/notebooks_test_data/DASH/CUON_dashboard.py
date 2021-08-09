#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd

import os

import plotly
import plotly.express as px
import plotly.graph_objects as go
#import dash  # (version 1.12.0) pip install dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

import numpy as np 
from jupyter_dash import JupyterDash
import dash_bootstrap_components as dbc

from plots_maker import plot_properties

#import dash_daq as daq

from plots_maker import make_map, make_series, make_table_harvested, make_bars 

from data_reader import *
# importing df, df_map 


#app = dash.Dash(__name__ , external_stylesheets= [dbc.themes.CYBORG])
#app = JupyterDash(__name__, external_stylesheets= [dbc.themes.CYBORG])
app = JupyterDash(__name__, )


# In[ ]:





# In[2]:


""" Setting the APP LAYOUT 
What goes inside the app layout is your dash components,
with the graphs, layouts, checkboxes
anything that is listed here: https://dash.plotly.com/dash-core-components 

"""

std_plevs    = [1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 85000, 92500, 100000]

title = html.H1(
    className="title",
    children='CUON Station Dashboard',
)



#files_df = get_files_df()


app.layout = html.Div([
    # title of our web page
    html.H1("CUON Station Dashboard", style={'text-align': 'center'}),
    html.Br(),  # Br is a break i.e. a space in between
    # label and the key is what actually the user is going to see
    # values are integer since they come from the pandas dataframe and they are integer there
            
        
    html.Div([
            
            html.Div([
                html.Label('Input Primary Id'),
                dcc.Input(id="station", type="text", 
                          placeholder="0-20000-0-72507", 
                          style={'marginRight':'10px'},
                          #list = station_list,
                          #autoComplete = True,
                          #value = "0-20000-0-72507"
                         ),
            
                #html.Label('Select Station on the Map'),
                dcc.Graph(id='map', figure={}, 
                          clickData={'primary_id': '0-20000-0-82930'} ),
            ], className="subcontainer1" ),
                      
                
            html.Div([ 
                    html.Label('Select Variable'),
                
                    dcc.Dropdown(id="variable",
                         options=[
                         {"label": "Air Temperature"       , "value": 85  },
                         {"label": "Dew Point Temperature" , "value": 34  },
                         {"label": "Dew Point Depression " , "value": 36  },
                         {"label": "Relative Humidity"     , "value": 38  },
                         {"label": "u-wind"                , "value": 104 },
                         {"label": "v-wind"                , "value": 105 },
                         {"label": "Wind Direction"        , "value": 106 },
                         {"label": "Wind Speed"            , "value": 107 },
                         ],

                    multi=False,
                    value=85, #initial value in the drop down menu
                    style={'width': "98%"},
                    ),
        
                    html.Br(),
                    html.Label('Select Data Type'),
                    dcc.RadioItems(id='kind',
                    options=[ {'label': "Observation"     , "value": "observation_value" },
                              {'label': "An. Departure"   , "value": "an_depar@body" },
                              {'label': "Fg. Departure"   , "value": "fg_depar@body" },
                              {'label': "Bias Correction" , "value": "biascorr@body" },
                              {'label': "Desroziers Unc." , "value": "desroziers_180" },


                            ] ,
                    #labelStyle = {"display":"inline-block"},
                    labelStyle = {"display":"block"},
                    value = "observation_value", 
                    style={'width': "98%"}, ),

                ], className="container-kind"  )   ,
                             
                    
                               
            html.Div([
                    html.Label('Select Pressure Level'),
                    dcc.RadioItems(id='plevel',
                    options=[ {'label': str(i) + ' Pa', 'value': i} for i in std_plevs ] ,
                    value= std_plevs[0] ,
                    labelStyle={'display': 'block'},   
                            )], 
                    style={'width': "95%"},
                    className="container-plevel") ,
                    
        
            html.Div([
                    #html.Label('Time Series'),
                    dcc.Graph(id='series', figure={} ),                  
                    ], className="container-plot-series")  , 
   
            html.Div([
                    #html.Label('Time Series'),
                    dcc.Graph(id='bars', figure={} ),                  
                    ], className="container-plot-bar")  , 
                
                
            html.Div([dcc.Graph(id='table-harvested', figure={} ), 
                     ],
                className="container-harvested"
            ),
        
        
        
        
        
    ], className="dash-container", ),       
    
], style={'width': '100%',
          'display': 'inline-block',
          #'padding-left': '30px'
         },)


# In[3]:


# "an_depar@body" , "fg_depar@body" , "biascorr@body" , "sonde_type@conv"


# In[4]:


@app.callback(
    [Output(component_id='map', component_property='figure'), 
     Output(component_id='series', component_property='figure'),
     Output(component_id='table-harvested', component_property='figure'),
     Output(component_id='bars', component_property='figure'),

    ],
    
    [Input(component_id='variable', component_property='value'),
     Input(component_id='plevel'  , component_property='value'),
     Input(component_id='kind'    , component_property='value'),
     Input(component_id='station' , component_property='value'),
     Input(component_id='map'     , component_property='clickData'),
    ]
    
    )

    
def update_plots(variable, plevel, kind, station, clickData):

    """ Variable: meteo variable,
        kind: obs, fg_dep, an_dep, uncert, (+ bias,...) """
    
    if not station:
        station = "0-20000-0-82930"
        
    df , sources, counts = read_file(station)
    
    df_f = filter_data(df, plevel, variable)
    
    df_harvested = read_harvested(station)
    

    
    """ Creating the charts """
    map = make_map(df_map)
    series = make_series(df_f["date_time"], df_f[kind], variable, kind, station ) 
    table_harvested = make_table_harvested(df_harvested)
    bars = make_bars(sources, counts)
    
    
    return [map, series, table_harvested, bars]

    #return [series] # NB must always return a list even if you have one output only, due to @app definition 


# In[5]:


if __name__ == '__main__':
    #app.run_server(mode = 'inline', debug=True)
    import socket
    host = socket.gethostbyname(socket.gethostname())
    app.run_server(mode='inline',  port=8055, debug=True, host = host)
    # app.run_server(host='0.0.0.0', debug=True)


