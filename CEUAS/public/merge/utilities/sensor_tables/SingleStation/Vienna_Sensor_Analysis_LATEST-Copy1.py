#!/usr/bin/env python
# coding: utf-8

# # Sensor Metadata Analysis for the CUON Dataset
# 
# This notebook provides a complete analysis of the sensor identification methods for a specific station.
# In particular, three different metadata sources are considered:
# 
# - WMO codes: the WMO codes are extracted from table at page A-398 of https://library.wmo.int/doc_num.php?explnum_id=10235 
# 
# - Schroeder's codes as can be found in XXX
# 
# - IGRA2 metadata that can be downloaded from https://www.ncei.noaa.gov/pub/data/igra/history/
# 
# 
# WMO codes are available in the 








#get_sensor_id_comments('VDT', sensor_conf)


# ## Analyze A Specific Station
# 
# Here we analyze a specific station.
# Reads data from the "/data" directory if data was already extracted in the past. 
# 
# Otherwise, the user must have access to the directory of the "merged"  files as indicated above.


def analyze_station(sensor_conf, station):
    
        
    def get_indices(data):
        data = data.reset_index()
        # find the indices where the sensor was replaced 
        indices = []
        last = ''
        for index, row in data.iterrows():
            sid = row.sensor_id
            #if sid =='nan':
            #    continue
            #print(index)
            if index ==0:
                indices.append(index)
                last = sid
            else:
                if sid == last:
                    continue
                else:
                    #print(sid , ' ' , last )
                    last = sid
                    #print(sid , ' ' , last )
                    indices.append(index)
        return indices
    
    # extracting station primary id 
    
    if not ( os.path.isfile('data/' + station + '_sch.csv') and os.path.isfile('data/' + station + '_wmo.csv')):
        
        print("Retrieving data from merged file")

        #file = [f for f in os.listdir(merged) if station in f and 'Sensor' not in f ][0]
        file = [f for f in os.listdir(merged) if station in f][0]

        station = file.split('/')[-1].split('_')[0]
        file = merged + '/' + file 

        f = h5.File(file, 'r')
        ts = f['recordtimestamp'][:]
        tsd = pd.to_datetime( ts, unit='s',  origin=pd.Timestamp('1900-01-01') )

        #index_minus = np.where(tsd <=  pd.Timestamp('1994-01-01')  )[0][-1]
        index_minus = 0   # change to start from a certain date onwards 
        
        #index_plus = np.where(tsd >  pd.Timestamp('1997-01-01')  )[0][0]
        index_plus = np.where(tsd <  pd.Timestamp('2013-01-01')  )[0][-1]


        # Extracting Schroeder 
        ind_obs_sch = list(f['recordindex'][:]) [index_minus:index_plus]
        i = np.take( f['observations_table']['sensor_id'][:].view('|S4') , ind_obs_sch) 
        ids_s = [s.decode('utf-8').replace('.0','').replace('.','') for s in i ]
        dic = {'date_time': tsd[index_minus:index_plus] , 'sensor_id': ids_s }

        data_sch = pd.DataFrame(dic)
        data_sch['value'] = 1

        # Extracting WMO
        ind_obs_wmo     = list(f['recordindex'][:]) [index_plus:]
        ind_obs_wmo_all = list(f['recordindex'][:]) # taking all WMOs
        
        
        wmoids = np.take(  f['observations_table']['sensor_id'][:].view('|S4') , ind_obs_wmo)
        wmoids = [s.decode('utf-8') for s in wmoids ]

        for s in np.unique(wmoids):
            print(s, '  ', type(s))
        
        dic_wmo = {'date_time':tsd[index_plus:] , 'sensor_id':wmoids }
        data_wmo = pd.DataFrame (dic_wmo)
        data_wmo['value'] = 2

        data_wmo.to_csv('data/' + station + '_wmo.csv' , sep = '\t') 
        data_sch.to_csv('data/' + station + '_sch.csv' , sep = '\t') 
            
        f.close()
            
    else:
        print("Loading existing data")
        data_wmo = pd.read_csv( 'data/' + station + '_wmo.csv' , sep = '\t')
        data_wmo['date_time'] = pd.to_datetime(data_wmo['date_time'] )

        data_sch = pd.read_csv('data/' + station + '_sch.csv' , sep = '\t') 
        data_sch['date_time'] = pd.to_datetime(data_sch['date_time'] )

        
    data_wmo['source'] = 'WMO'
    data_sch['source'] = 'SCH'

    # cleaning WMO data from nans 
    data_wmo_clean = data_wmo.loc[ (data_wmo.sensor_id != 'nan') & (data_wmo.sensor_id != '-922')].dropna( subset=['sensor_id'])
    data_wmo_clean.reset_index()
    
    #print(data_wmo_clean[data_wmo_clean.date_time >=  pd.Timestamp('1994-11-02') ][:20])
    
    indices_wmo_clean = get_indices(data_wmo_clean)
    #print(indices_wmo_clean)
    
    data_wmo_clean = data_wmo_clean.iloc[indices_wmo_clean]
    
    # getting only variation int he sensor_id indices 
    indices_sch = get_indices(data_sch)
    indices_wmo = get_indices(data_wmo)

    data_df = pd.concat( [data_sch.iloc[ list(indices_sch)], data_wmo. iloc[ list(indices_wmo)] ] )
    data_df_clean = pd.concat(  [data_sch.iloc[ list(indices_sch)], data_wmo_clean ]  ) # removed nans

    #unique_ids = np.unique(data['sensor_id'])
    comments = [ str(get_sensor_id_comments(str(i).replace(' ','').replace('.0',''), sensor_conf)) for i in data_df.sensor_id]

    data_df['comment'] = comments
    sid_clean = [str(i).replace('.0','')  for i in data_df.sensor_id]
    data_df['sensor_id'] = sid_clean    


    # Adding comments 
    for d in [data_df, data_df_clean]:
        comments = [ str(get_sensor_id_comments(str(i).replace(' ','').replace('.0',''), sensor_conf)) for i in d.sensor_id]
        d['comment'] = comments
        sid_clean = [str(i).replace('.0','')  for i in d.sensor_id]
        d['sensor_id'] = sid_clean    

        
    return data_sch, data_wmo, data_df, data_wmo_clean, data_df_clean 



# Extract all sensor data
data_sch, data_wmo, data_df, data_wmo_clean, data_df_clean = analyze_station(sensor_conf, station)



def make_time_series(data_df, station_name, label = ''):
    #filter date
    #data_df = data_df.loc[data_df.date_time <= pd.Timestamp('1995-01-01')]

    # converts categorical char values for sensor ids to integer values (on y-axis in the plots)
    data_df['values'] = pd.factorize( data_df.sensor_id)[0]


    fig = px.scatter(data_df, x="date_time", y="values", color="sensor_id",
                    hover_name="sensor_id", hover_data=["comment"]
                    )

    fig.update_layout(title='Sensors Time Series - ' + station_name + ' ' + label)
    #fig.update_yaxes( ticktext= ['Schroeder', 'WMO', 'IGRA2'],
    #                  tickvals= [1,2,3])

    fig.update_layout(width= 1800, height = 400)

    fig.update_layout(
        xaxis_title="Date of Sensor Replacement",
        yaxis_title="Sensor Metadata",
        legend_title="Sensor ID",
        font=dict(
            size=16,
            color="black"
        )
    )

    fig.update_traces(marker=dict(size=14,
                                  line=dict(width=2,
                                            color='DarkSlateGrey')),
                      selector=dict(mode='markers'))




    for d in stat_igra2.date_time:
        fig.add_vline(x=d, line_width=3, line_dash="dash", line_color="green")



    fig.update_layout(hovermode="x unified")

    fig.update_layout(
    yaxis = dict(
    tickfont = dict(size=16)))

    return fig
# fig = make_time_series(data_df, station_name) 


import json
from plotly.subplots import make_subplots

def make_time_series_2(data_df, station_name, label = ''):
    #filter date
    #data_df = data_df.loc[data_df.date_time <= pd.Timestamp('1995-01-01')]
    from plotly.subplots import make_subplots

    with open('/mnt/users/staff/leo/python/CEUAS/CEUAS/public/adjust/feedbackmerged070350_breakanalysis.json') as f:
        d=json.load(f)
        time = pd.to_datetime(d['days_since_1900'] , unit='d', origin=pd.Timestamp('1900-01-01') )
    
    snht = pd.DataFrame( {'time': time , 'snht':d['tsasum'] } )
    
    symbols = {"IGRA2":'star', "WMO":'circle', "SCH":'square'}

    #markers = {"IGRA2":0, "WMO":17, "SCH":22}
    #data_df['symbols'] = [markers[i] for i in data_df.source ]
    
    #print(data_df['symbols'])
    # Create figure with secondary y-axis
    subfig = make_subplots(specs=[[{"secondary_y": True}]])
    fig1 = px.line(snht, x="time", y="snht")
    fig2 = px.scatter(data_df, x="date_time", y="value", color="sensor_id",
                    hover_name="sensor_id", hover_data=["comment"],
                    symbol="source",
                    symbol_map= symbols )
    
    fig2.update_traces(yaxis="y2")



    

    subfig.add_traces(fig1.data + fig2.data)
    subfig.layout.xaxis.title=""
    subfig.layout.yaxis.title="SNHT"
    subfig.layout.yaxis2.title="Metadata Source"
    # recoloring is necessary otherwise lines from fig und fig2 would share each color
    # e.g. Linear-, Log- = blue; Linear+, Log+ = red... we don't want this
    subfig.for_each_trace(lambda t: t.update(line=dict(color=t.marker.color)))
    #subfig.show()

    subfig.update_layout(title='Sensors Time Series - ' + station_name + ' ' + label)

    subfig.update_layout(width= 1800, height = 550)


    subfig.update_traces(marker=dict(size=14,
                                  line=dict(width=2,
                                            color='DarkSlateGrey')),
                      selector=dict(mode='markers'))



    igra2 = data_df.loc[data_df.source == 'IGRA2']
    for d in igra2.date_time:
        subfig.add_vline(x=d, line_width=3, line_dash="dash", line_color="green")



    subfig.update_layout(hovermode="x unified")

    subfig.update_layout(
    yaxis = dict(
    tickfont = dict(size=16)),
    font=dict(
        size=16,
        color="black"
        )
    )

    subfig.update_yaxes( ticktext= ['Schroeder', 'WMO', 'IGRA2'],
                      tickvals= [1,2,3], secondary_y=True )
        
    return subfig


fig = make_time_series_2(data_df_clean, station_name) 

if not os.path.isdir('plots'):
    os.mkdir("plots")



import plotly as plotly
time_series_file =  "plots/" + station_name + "_timeSeries" 
plotly.offline.plot(fig, filename= time_series_file + ".html" )

import plotly.io as pio
pio.write_image(fig, "plots/" + station_name + "_timeSeries_ku.png")





fig = make_time_series(data_df_clean, station_name, label=' [Cleaned,Only WMO Variations]') 
#fig.show()

nans_ind = np.where( (data_df.sensor_id != 'nan') & (data_df.sensor_id != '-922') &  (data_df.source != 'IGRA2')  )[0]
data_sensor_clean = data_df.iloc[nans_ind]
data_sensor_clean_tab = data_sensor_clean.drop_duplicates(subset=['sensor_id'])


# ## Summary table of the sensor


fig = go.Figure(data=[go.Table(
    header=dict(values=list(['Sensor','Source', 'Comment']),
                fill_color='gold',
                align='left',
                font_size=20),
    columnwidth = [40,40,300],
    cells=dict(values=[data_sensor_clean_tab.sensor_id, data_sensor_clean_tab.source, data_sensor_clean_tab.comment],
               fill_color='aliceblue',
               align='left',
               font_size=16,
               height=30
              )),
])

fig.update_layout(width=1900, height=50*len(df_sensor))



# ## Counts and time series of the WMO sensors


df = data_sensor_clean.groupby(["sensor_id"]).count()
fig = px.bar(df, x=df.index, y='value', title = 'WMO codes counts')





