#!/usr/bin/env python
# coding: utf-8

from multiprocessing import Pool
from functools import partial

# importing the modules for data reading and analysis
from sensor_functions import * 
from plot_functions_sensor import * 

import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)


warnings.simplefilter(action='ignore', category=FutureWarning)



### Testing
# igra2_metadata.head(20)


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
#app = JupyterDash(__name__, external_stylesheets= [dbc.themes.CYBORG])
app = JupyterDash(__name__, )


""" Setting the APP LAYOUT 
What goes inside the app layout is your dash components,
with the graphs, layouts, checkboxes
anything that is listed here: https://dash.plotly.com/dash-core-components 

"""

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



"""
@app.callback(
        Output("loading-output-1", "value"),
    
        [ Input("input_station", "value") ]
)

    
"""

# Loading IGRA2 metadata 
ig = IgraMetaData()
igra2_metadata = ig.igra2_meta_df
    
# Loading sensor configuration
sensor = Sensor()

# Merged file source (if data not already available)
merged = '/scratch/das/federico/MERGED_APRIL2022'

all_stat = os.listdir(merged)
all_stat = [s.split('_')[0] for s in all_stat ]


def get_data(station):
    
    if not station:
        station = default_station
    
    if not os.path.isdir('data_plot'):
        os.mkdir('data_plot')
        
    if not (os.path.isfile('data_plot/' + station + '_data_clean_all.csv') 
            and os.path.isfile('data_plot/' + station + '_all_sensor_station.csv') ):

        logging.debug(" --- RETRIEVING --- data file: ")

        station = [c for c in all_stat if station in c ][0] 

        stat_igra2, stat_igra2_sonde = ig.get_igra_metadata(station)

        # Analyze data
        logging.debug(" --- ANALYZING --- data file: ")
        analyze = Analyze(sensor,merged,station)
        data_sch, data_wmo, data_df, data_wmo_clean, data_df_clean = analyze.analyze()

        data_clean_all = pd.concat([data_df_clean, stat_igra2_sonde])

        
        # extract unique sensor id table for the station
        all_sensor_station = analyze.get_all_sensors(data_clean_all)

        all_sensor_station.to_csv('data_plot/' + station + '_all_sensor_station.csv' , sep='\t')
        data_clean_all.to_csv('data_plot/' + station + '_data_clean_all.csv' , sep='\t')
        
        
    else:
        logging.debug(" --- READING --- data file: ")
        data_clean_all = pd.read_csv('data_plot/' + station + '_data_clean_all.csv', sep='\t')
        all_sensor_station = pd.read_csv('data_plot/' + station + '_all_sensor_station.csv', sep='\t')
        
            
    
    #print(data_df_clean_all.head(10) , all_sensor_station_df.head(10) )
    return [data_clean_all , all_sensor_station ]


station =     default_station
#data_clean_all , all_sensor_station = get_data(station)

### Testing
#station = default_station
#data_clean_all , all_sensor_station = get_data(station)
#data_clean_all
#all_sensor_station

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
    all_data , all_sensor = get_data(station)

    
    # Loading the chart class """
    plot = Plot(station.split('_')[-1], save=False)
    
    # Creating the plots
    sensor_table = plot.sensor_table( all_sensor)
    series = plot.time_series( all_data, label='')
    wmo_table = plot.wmo_bar_plot(all_data)
    
    return [series, wmo_table, sensor_table]



if __name__ == '__main__':
    #app.run_server(mode = 'inline', debug=True)
    import socket
    host = socket.gethostbyname(socket.gethostname())
    app.run_server(port=8058, debug=True, host = host)
    # app.run_server(host='0.0.0.0', debug=True)

