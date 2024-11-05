import plotly.express as px  # (version 4.7.0)
import plotly.graph_objects as go
import plotly.express as px

from utils import load_data

import dash
from dash import dash_table
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
import pandas as pd
from datetime import date

app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])

import numpy as np

import time


data_df, sensors = load_data(dir = 'data_processed')
data_json = data_df.to_json(date_format='iso', orient='split')


### https://www.youtube.com/watch?v=vqVwpL4bGKY
### How to style the app

"""
ROWS are defined first, and then we define the columns


NOTE 
the styling of the app, particularly the layout, is influenced by the css file!
"""
first_date = date(2005, 2, 7)

options_width = "100%"

app.layout = html.Div(children=[

                ### dcc.Store stores the input data (whole dataframe)
                dcc.Store(id='data-dataframe',
                          data=data_json
                          ),

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





""" Updating Dropdown """
@app.callback(
    [Output('time-dropdown', 'options'),
     Output('time-dropdown', 'value')],
    Input('date-picker', 'date'))

def dropdown_options(date):

    #print('Searching for date::: ', date )
    date_times = pd.read_csv('all_timestamps.csv' , sep='\t' )

    ts = date_times.loc[date_times.date == date].time.values

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
    Input('date-picker', 'date'),
    Input('time-dropdown', 'value'),
    Input('data-dataframe', 'data'),
    Input('variable', 'value'),

                  )

def plot_line_charts(start_date, value, data, variable):
    """ PLOTLY line chart """
    #print('Doing something :::' ,  data )
    print('Restoring DF from json ::: ' )
    data = pd.read_json(data, orient='split')

    #print('READ IN GRAPH DATA aaaaa :::' ,  data )
    #print("DATE - " , start_date, ' TIME ' , value )
    #print('READ IN GRAPH DATA:::' ,  data.head() )

    sensors = np.unique(data.sensor )

    #print('READ IN GRAPH SENSORS:::' ,  data )

    #print("DATE - " , start_date, ' TIME ' , value )
    #print('KEYS: ' , all_data.keys() )
    #dates = [ k.split(' ')[0]for k in list( all_data.keys() ) ]
    #print('DATES ' , dates )
    #ddo = pd.to_datetime(start_date).date()

    date_string = str( pd.to_datetime(start_date).date() )

    #print('DATE IS ' , date_string )
    #boo = bool(date_string in dates)
    #print('date selected in dates? ' , boo )

    #print(' unique date times ' , np.unique (data.date.values )    )

    #all_data_date = all_data.loc[  all_data.date == pd.to_datetime( start_date, format="%Y-%m-%d",  utc=True ) ]

    data_ts = data.loc[  (data.date == start_date) & (data.time == value) ]

    #print('ALL TIMES FOR THIS DF ::: , ' ,  np.unique(data_ts.datetime.values) )
    #print('ALL TIMES FOR THIS DF ::: , ' ,  np.unique(data_ts.time ) )
    #print('ALL TIMES FOR THIS DF ::: , ' ,  np.unique(data_ts.date ) )


    fig = go.Figure()


    # Getting the sensors


    print('Available SENSORS: ' , sensors )

    '''
    fig.add_trace(go.Scatter(
        x=data_ts.temp,
        y=data_ts.press,
        color = 'sensor',
        connectgaps=True,
        # line=dict(color='blue', width=2),
        name='Temp'  # override default to connect the gaps
    ) )
    '''

    if variable == 'Temperature':
        x='temp'
        xlabel = 'Temperature [K]'
    else:
        x='hum'
        xlabel='Relative Humidity'

    ylabel = "Pressure [Pa]"

    fig = px.line(data_ts, x=x, y='press', color='sensor')

    fig.update_layout(
        title={
            'text': variable + " - Date: " + start_date + ' / Launch: ' + value ,
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



    return fig





""" Launch the app """
if __name__ == '__main__':
    app.run_server(debug=True,
                   dev_tools_hot_reload=True,
   		   port='8050',
    ) # avoid auto reloadind = False


