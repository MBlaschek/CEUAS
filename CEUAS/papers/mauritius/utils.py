import os,sys
import pandas as pd
import numpy as np
from tqdm import tqdm



def load_data_old(dir='data_processed'):
    """ Retrieve the data and store a pandas dataframe """

    files = [ dir+'/'+f for f in os.listdir(dir) if 'fulldata' in f ]

    # will create a dict for each date-sensor pair
    all_df = {}

    for f in tqdm(files[:1]):
        df = pd.read_csv(f, sep='\t')
        sensor = f.split('/')[-1].split('_')[0]
        unique_ts = np.unique(df.datetime)
        for dt in unique_ts:
            if dt not in all_df.keys():
                all_df[dt] = {}

            df_dt = df.loc[df.datetime == dt ]

            dt_date = str( pd.to_datetime(dt).date() ) # change this !
            all_df[dt_date]= {}
            df =df[['press','hum','temp','datetime']]
            all_df[dt_date][sensor] = df

    return all_df


def load_data(dir='data_processed', size='reduced'):
    """ Retrieve the data and store a pandas dataframe """

    files = [ dir+'/'+f for f in os.listdir(dir) if size in f and '~' not in f  ]

    # will create a dict for each date-sensor pair
    all_df = []

    sensors = [  f.split('/')[-1].split('_')[0] for f in files]

    for f in tqdm(files[:]):  #TODO
        print('Extracting data from file::: ' , f )
        df = pd.read_csv(f, sep='\t')
        df =df.iloc[::5,:]  # TO DO change to exploit full data


        sensor = f.split('/')[-1].split('_')[0]

        dd =  pd.to_datetime(df['datetime']).dt.date  ### THIS IS VERY SLOW

        df['date'] = [d.strftime("%Y-%m-%d") for d in dd]   # -> format 2005-02-21

        tt =  pd.to_datetime(df['datetime']).dt.time
        df['time'] = [t.strftime("%H:%M:%S") for t in tt]   # -> format 15:00:54
        df = df[['press','hum','temp','datetime', 'date', 'time']]
        df['sensor'] = sensor


        all_df.append(df)

    DF = pd.concat(all_df)

    return DF, sensors


def load_data_single_timestamp(dir = 'data_processed', timestamp=''):
    """ Load the entire data for the given date """

    files = ['data_processed/single_ascent/' + f for f in os.listdir('data_processed/single_ascent/') if 'date' in f ]
