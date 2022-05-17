#!/usr/bin/env python
# coding: utf-8

import sys
sys.path.append('modules')

# importing the modules for data reading and analysis
from sensor_functions import * 
from plot_functions_sensor import * 

import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)


warnings.simplefilter(action='ignore', category=FutureWarning)


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


#app = dash.Dash(__name__ , external_stylesheets= [dbc.themes.CYBORG])
app = JupyterDash(__name__, )


""" Setting the APP LAYOUT """

title = html.H1(
    className="title",
    children='Sensor Id Metadata Dashboard',
)

#files_df = get_files_df()


default_station = "0-20000-0-82930" # "0-20000-0-82930" is fast, "0-20000-0-06610"

app.layout = html.Div([
    # title of our web page
    html.H1("Sensor Id Metadata Dashboard", 
            style={'text-align': 'center'}),
    
    html.Br(),  # Br is a break i.e. a space in between
        
    
    html.Div([
            
        html.Div( children=[
            html.Label('Station Primary Id'),
            dcc.Input(id="input_station", 
                      type="text", 
                      placeholder= default_station, 
                      style={'marginLeft':'20px'},
                      value = default_station,
                      debounce=True, # to wait for click on Enter to start data retrieval
                     ),
            dcc.Loading( id="loading-1",
                          type="circle",
                          children=html.Div(id="loading-output-1")
            ),
        ]),

        html.Div( children=[
                    dcc.Graph(id='series', figure={},   # how to make two dcc close by: use inline-block and inside same html.Div as children
                                  style={'display': 'inline-block' , 
                                         'width': "71%",
                                         'marginRight':'-10px',
                                         'marginLeft':'-10px'}),
                    
                    dcc.Graph(id='wmo_table', figure={}, 
                                  style={'display': 'inline-block',
                                         'width': "28%",
                                         'marginLeft':'-10px'}),             
                    ],
                style={'width': '100%',
                       'display': 'inline-block',
                         },
                className="dash-container", ) ,
        
        
        html.Div([                  
                    dcc.Graph(id='sensor_table', figure={}, 
                                  style={'display': 'inline-block',
                                         'width': "100%",
                                         'marginLeft':'-10px'}),             
                    ],
                style={'width': '100%',
                       'display': 'inline-block',
                         },
                className="dash-container", ) ,

        ])
    ])



# setting a deafult station to be loaded first
station =     default_station


""" Calling the APP """

@app.callback(
    [ Output(component_id='series', component_property='figure'),
      Output(component_id='wmo_table', component_property='figure'),        
      Output(component_id='sensor_table', component_property='figure'),
    ],
    
    [ Input(component_id='input_station' , component_property='value')
    ]
    
)

    
def update_plots(station):
    """ Function to call back """ 

    # Extracting all the data
    data_clean_all, all_sensor_station, data_all_wmo = get_data(station, force_create=False)

    
    # Loading the chart class """
    plot = Plot(station.split('_')[-1], save=False)
    
    # Creating the plots
    sensor_table = plot.sensor_table( all_sensor_station)
    series = plot.time_series( data_clean_all, label='')
    wmo_table = plot.wmo_bar_plot(data_all_wmo)
    
    return [series, wmo_table, sensor_table]


if __name__ == '__main__':
    #app.run_server(mode = 'inline', debug=True)
    import socket
    host = socket.gethostbyname(socket.gethostname())
    #app.run_server(port=8058, debug=True, host = host)
    app.run_server(debug=True, host = host)                                        
    # app.run_server(host='0.0.0.0', debug=True)

