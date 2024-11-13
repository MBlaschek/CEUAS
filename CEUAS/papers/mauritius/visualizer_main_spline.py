import plotly.express as px  # (version 4.7.0)

import plotly.graph_objects as go
import plotly.express as px

from utils import load_data_df

import dash
from dash import dash_table
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
import dash_ag_grid as dag

import pandas as pd
from datetime import date

app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])

import numpy as np
from utils import *

import time


data_df = load_data_df( dir='spline_data', ts='20050207-10:01', var='temp' , sensor='Meisei')
data_json = data_df.to_json(date_format='iso', orient='split')


#data_df, sensors = load_data(dir = 'data_processed')
#data_json = data_df.to_json(date_format='iso', orient='split')


### https://www.youtube.com/watch?v=vqVwpL4bGKY
### How to style the app

"""
ROWS are defined first, and then we define the columns


NOTE 
the styling of the app, particularly the layout, is influenced by the css file!
"""

### GLOBAL VARIABLES
data_dir = 'spline_data'
all_files = [data_dir+'/'+f for f in os.listdir(data_dir)]

first_date = date(2005, 2, 7)



options_width = "100%"




app.layout = html.Div(children=[

                ### -------------------------- Title
                html.Div( children=[
                        html.H1("Mauritius (VACOAS 2005) Visualizer",
                            style={'color': 'blue',
                                   'fontSize': 40,
                                   }),
                        html.H3 ("University of Vienna, Department of Meteorology and Geophysics "),
                    ], className = 'title',
                       style = {'textAlign':'center',
                                'marginTop': '10px',
                                'marginBottom': '50px',
                                'fontSize': 15,
                                }
                ),

                html.Br(),
                html.Br(),

                ### -------------------------- Selection Box
                html.Div(children=[
                        html.Div(children=[
                            html.H2('Select the Data to Visualize',
                                    style={'paddingTop': '2rem',
                                           'fontSize': 30,
                                           }
                                    ),
                            html.Div( children=[
                                html.Label('Date',
                                           style={'paddingTop': '2rem',
                                                  'fontSize': 25,
                                                  }
                                           ),
                                dcc.DatePickerSingle(
                                        id='date-picker',
                                        min_date_allowed=date(2005, 2, 7),
                                        max_date_allowed=date(2005, 2, 25),
                                        initial_visible_month=date(2005, 2, 7),
                                        disabled_days = [date(2005, 2, 13)],
                                        date= first_date,
                                style={
                                        'font-size': '25',
                                        'width':'100%',
                                        'display': 'inline-block',
                                        'border-radius': '2px',
                                        'border': '1px solid #ccc',
                                        'color': '#333',
                                        'border-spacing': '0',
                                        'border-collapse': 'separate',
                                        #'text-align':'center',
                                        'padding-left': '3px',
                                        'padding-right': '3px',
                                        'paddingLeft': '2rem'
                                },

                                )
                                ]),

                            html.Div(children=[
                                html.Br(),
                                html.Label("Launching Time ",
                                            style={'paddingTop': '2rem',
                                            'fontSize': 25,
                                                            }
                                            )
                            ,
                                # gets updated after selecting the date
                                dcc.Dropdown(id='time-dropdown',
                                             style={'paddingLeft': '2rem',
                                                    'fontSize': 25,
                                                    }
                                             ),
                            ]),


                            html.Div(children=[
                                html.Br(),
                                html.Label("Variable",
                                           style={'paddingTop': '2rem',
                                                  'fontSize': 25,
                                                  }                                            ),
                                dcc.RadioItems(['Temperature', 'Humidity'], 'Temperature',
                                           id='variable',
                                               style={'paddingLeft': '2rem',
                                                      'fontSize': 25,
                                                      }
                                               ),
                            ]),

                            html.Div(children=
                                [
                                html.Br(),
                                    html.Label(
                                        "Click to export CSV data",
                                    style = {'paddingTop': '2rem',
                                             'fontSize': 25,
                                             },
                                    ),

                                    dbc.Button("Download CSV", id="download-button", n_clicks=0,
                                                size='lg', color='primary'),

                                    html.Div(id='container-button-timestamp'),

                                    dash_table.DataTable(id='dtable',
                                        columns=[{"name": i, "id": i} for i in data_df.columns ],
                                        #page_size=10 ### limit 10 rows per table, use pagination
                                        page_action='none',
                                        style_table={'height': '300px', 'overflowY': 'auto'},
                                        fixed_rows = {'headers': True},),

                                    dcc.Download(id='download_component'),

                                ]
                            ),



                        ], className = "four columns",
                           style = {
                                'padding':'2rem',
                                'margin':'1rem',
                                'boxShadow': '#e3e3e3 4px 4px 2px',
                                'border-radius': '10px',
                                'marginTop': '2rem',
                                }

                            ),

                            ### -------------------------- Plot Box
                        html.Div(
                            dcc.Graph(id='linechart', figure={},
                                  config={'scrollZoom': False})
                        , className = 'seven columns',
                            style={
                                   #'width':'60%'
                        }),
                ], className="twelve columns"),
            ]
)



#all_data = load_data(dir = 'data_processed')
#dates = [pd.to_datetime(k).date() for k in all_data.keys() ]
#print(dates)

### Available timestamps for date 2005-02-14  are:::  ['05:06:57' '10:08:04' '15:00:54' '18:00:36']


'''
""" Reading the dataframe """
@app.callback(
     [Output('data-dataframe', 'data')],
    Input('date-picker', 'date')
)

def read_dataframe(date):

    print('READ beginning:::' , data_df.head(), np.unique(data_df.sensor ) )

    #cleaned_df.to_json(date_format='iso', orient='split')
    a = data_df.to_json(date_format='iso', orient='split')
    print('Done putting to json')
    return a
'''

date_times = pd.read_csv('all_timestamps.csv' , sep='\t' )


""" Updating Dropdown """
@app.callback(
    [Output('time-dropdown', 'options'),
     Output('time-dropdown', 'value')],
    Input('date-picker', 'date'))

def dropdown_options(date):

    #print('Searching for date::: ', date )
    #date_times = pd.read_csv('all_timestamps.csv' , sep='\t' )
    ts = date_times.loc[date_times.date == date].time.values
    print('TIMESTAMPS ::: ' , ts )

    print('Available timestamps for date ', date , ' are::: ' , ts )

    return ts, ts[0]


'''
@app.callback(Output("loading-output-date", "children"),
          Input("date-picker", "date"))
def input_triggers_spinner(value):
    time.sleep(1)
    return value
'''




""" Plotting maps"""
@app.callback(    Output(component_id='linechart', component_property='figure'),
                  Output(component_id='dtable', component_property='data'),
    Input('date-picker', 'date'),
    Input('time-dropdown', 'value'),
    Input('variable', 'value'),

                  )

def plot_line_charts(start_date, time, variable):
    """ PLOTLY line chart """
    #print('Doing something :::' ,  data )

    """
    Available timestamps for date  2005-02-07  are:::  ['10:01:53' '14:57:50' '19:31:41']
    """

    print('DF for time stamp ::: ' , start_date, '  ', time[:5] )

    ### Adapting the variables to the file names
    time=time[:5]
    start_date = start_date.replace('-','')
    if variable == 'Temperature':
        var = 'temp'
        xlabel = 'Temperature [K]'

    else:
        var = 'hum'
        xlabel='Relative Humidity'

    files = [f for f in all_files if var in f and start_date in f and time in f ]
    sensors = [f.split('/')[1].split('_')[4].replace('.csv','') for f in files ]


    print('FILES ::: ' , files )
    print('Available SENSORS: ' , sensors )

    fig = go.Figure()


    # make DF
    all_df = []
    for f in files:
        df = pd.read_csv(f, sep='\t')
        dfr=df[::10]
        dfr['sensor'] = f.split('/')[1].split('_')[4].replace('.csv','')
        all_df.append(dfr)

    data_df = pd.concat(all_df)

    print(data_df.head(), ' ::: ' , len(data_df))


    fig = px.line(data_df, x='value', y='height', color='sensor')

    ylabel = "Height [m]"

    fig.update_layout(
        title={
            'text': variable + " - Date: " + start_date + ' / Launch: ' + time ,
             'font':dict(size=35),
             "yref":'paper'
        },

        width = 1200,
        height = 800,
        yaxis=dict(title = ylabel,
                   #autorange="reversed",
                   #range=[100000,0],
                   #type="log"
                   ),
        xaxis=dict(title= xlabel ),
        font=dict(
            #family="Courier New, monospace",
            size=22,
            #color="RebeccaPurple"
        )


    )


    data=data_df.to_dict("records")

    return fig,data






'''
@app.callback(
    Output("container-button-timestamp", "children"),
    Input("csv-button", "n_clicks"),
    prevent_initial_call=True
)
def displayClick(n_clicks):
    msg="I am doing something"
    return html.Div(msg)
'''

@app.callback(
    Output('download_component', "data"),
    Input('download-button', "n_clicks"),
    State('dtable', "derived_virtual_data"),
    prevent_initial_call=True,
)
def download_data(n_clicks, data):
    dff = pd.DataFrame(data)
    return dcc.send_data_frame(dff.to_csv, "data_csv.csv")



""" Launch the app """
if __name__ == '__main__':
    app.run_server(debug=True,
                   dev_tools_hot_reload=True,
   		   port='8050',
    ) # avoid auto reloadind = False


