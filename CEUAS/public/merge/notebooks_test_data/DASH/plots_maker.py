import os,sys
import plotly.graph_objects as go
import plotly.express as px
import numpy as np
import pandas as pd 

plot_properties = {  85: { 'name' : "Air Temperature" , "u" : "[K]" },
                           
                     34: { 'name' : "Dew Point Temperature" , "u" : "[K]" },    
                     36: { 'name' : "Dew Point Depression"  , "u" : "[K]" },    
                     38: { 'name' : "Relative Humidity"     , "u" : "" },    

                     104: { 'name' : "u-wind Component"     , "u" : "[m/s]" },       
                     105: { 'name' : "v-wind Component"     , "u" : "[m/s]" },       
                     106: { 'name' : "Wind Direction"       , "u" : "" },       
                     107: { 'name' : "Wind Speed"           , "u" : "[m/s]" },       
                           
                           }
    



def make_series(times, values, variable, kind, station):
    """ Make time series plot """
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=times, y=values,
                    mode='lines+markers',
                            )
                 )
    
    
    fig.update_layout(
    title= kind + " , STATION " + station,
    #xaxis_title = "X Axis Title",
    yaxis_title= plot_properties[variable]["name"] + " " + plot_properties[variable]["u"],
    #legend_title="Legend Title",
        
    #font=dict(
    #    family="Courier New, monospace",
    #    size=18,
    #    color="RebeccaPurple"
    #)
)
    fig.update_xaxes(rangeslider_visible=True)
    return fig



def make_map(df):
    
    #df["color"] = "blue"
    
    
    fig = go.Figure()

    fig=px.scatter_geo(df,lon='longitude', lat='latitude',
                       opacity=0.8,
                       #projection="natural earth",
                       hover_data={'longitude':True,'latitude':True,'name':True, 'primary_id':True} )
                       
                       
    
    """
    for s in ["0-20000-0-06610" , "0-20001-0-10393" ,"0-20001-0-11035" ]:
        a = np.where( df.primary_id == s )[0]

        df.at[a, 'color'] = 'red'
        fig.add_trace(
            go.Scattergeo(lon=df["longitude"], 
                          lat=df["latitude"],
                          color = df["color"]
                       #projection="natural earth",
                       #hover_data={'longitude':True,'latitude':True,'name':True, 'primary_id':True} )
                        ) )
        
        #fig.update_traces(marker=dict(color="red" , size = 10 ) )
    """
        
    fig.update_layout(
    margin=dict(l=0, r=0, t=0, b=0),
    paper_bgcolor="lime",)
    
    return fig 


def make_table_harvested(df):
    """ Making a table oof the harvested files,
         res = {"file": [] , "counts":[], "dataset": [], "lat":[] , "lon":[], "min_date":[], "max_date": [] } """
    
    # print(df.columns)
    #  'dataset', 'lat', 'lon', 'min_date', 'max_date', 'counts', 'file'
    #df = df[ ['dataset', 'lat', 'lon', 'min_date', 'max_date', 'counts', 'file']]
    
    fig = go.Figure(data=[go.Table(
    header=dict(values= ['dataset', 'lat', 'lon', 'min_date', 'max_date', 'counts', 'file'],
                fill_color='paleturquoise',
                align='center',
                font_size=14),
        
        
    columnwidth = [80,70,70,100,100,100,400],
    cells=dict(
               fill_color='lavender',
               values=[df.dataset, df.lat, df.lon, df.min_date, df.max_date, df.counts, df.file],
               align='left',
               font_size=12,) 
    
    )
    ])
    fig.update_layout(title_text= 'Summary of surce files')
    return fig 


def make_bars(sources, counts):
    """ Create a bar chart for the distributions of source_id 
    (counting number of records)"""
    
    fig = go.Figure()

    df = pd.DataFrame.from_dict( {"Source": sources , "Counts": counts } )
        
    fig = px.bar(df, y='Counts', x='Source', color='Source',  text='Counts' )
    fig.update_traces(texttemplate='%{text:.2s}', textposition='outside')
    fig.update_layout(uniformtext_minsize=8, uniformtext_mode='hide')
    
    return fig