import pandas as pd
import os,sys 
import numpy as np

# default station if nothing has been selected
stat_def = '0-20000-0-82930'

database_path = '/raid60/scratch/federico/DATABASE_JULY2021_STD/'

station_list = [f.replace(".csv","") for f in os.listdir(database_path) if ".csv" in f and "station" not in f ]

station_map = database_path + '/stations_list_dash.csv'

df_map = pd.read_csv(station_map, sep = '\t' )



def get_bars_sources(df):
    """ Counts the data sources for the bar histogram """
    
    sources, counts = np.unique(df["source_id"], return_counts = True)
    sources = [ s.decode("utf-8") for s in sources ]
    return sources, counts


def read_file(station):
    """ Input from dashboard:
        - station
        - plevel
        - variable """

    if station not in station_list:
        station = stat_def
    try:
        file_name = database_path + '/' + station + '.pkl'
    except:
        file_name = database_path + '/' + stat_def + '.pkl'

        
    df = pd.read_pickle(file_name)
    df["date_time"] = pd.to_datetime(df['date_time'], unit='s',origin=pd.Timestamp('1900-01-01'))
    
    sources, counts = get_bars_sources(df)
    
    return df , sources, counts 


#df = read_file(station)

def filter_data(df, plevel, variable):
    """ Input from dashboard:
        - station
        - plevel
        - variable """
    
    df = df[ (df.z_coordinate == plevel) & (df.observed_variable == variable) ]
    
    return df


def read_harvested(station):
    file = database_path + "/" + station + "_summary_harvested_files.csv"
    df = pd.read_csv( file, sep = "\t" )

    return df



