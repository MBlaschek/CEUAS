#!/usr/bin/env python
# coding: utf-8

# In[31]:


import os,sys
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas as pd
import numpy as np

from multiprocessing import Pool
from functools import partial

import h5py as h5

import urllib.request
pd.options.mode.chained_assignment = None 


from multiprocessing import Pool
from functools  import partial

import logging


if not os.path.isdir('data'):
    os.mkdir('data')
    
class IgraMetaData():
    """ Class holding the basic functionality to produce extract IGRA2 metadata information """
    
    def __init__(self):

        """ Convert the IGRA2 metadata from plain txt into a dataframe,
        for the select station WMO id """

        """
        IGRAID         1- 11   Character
        WMOID         13- 17   Integer
        NAME          19- 48   Character
        NAMFLAG       50- 50   Character
        LATITUDE      52- 60   Real
        LATFLAG       62- 62   Character
        LONGITUDE     64- 72   Real
        LONFLAG       74- 74   Character
        ELEVATION     76- 81   Real
        ELVFLAG       83- 83   Character
        YEAR          85- 88   Integer
        MONTH         90- 91   Integer
        DAY           93- 94   Integer
        HOUR          96- 97   Integer
        DATEIND       99- 99   Integer
        EVENT        101-119   Character
        ALTIND       121-122   Character
        BEFINFO      124-163   Character
        BEFFLAG      164-164   Character
        LINK         166-167   Character
        AFTINFO      169-208   Character
        AFTFLAG      209-209   Character
        REFERENCE    211-235   Character
        COMMENT      236-315   Character
        UPDCOM       316-346   Character
        UPDDATE      348-354   Character
        """
        
        # obtain the metadata file 
        if not os.path.isfile('igra2-metadata.txt'):
            url = 'https://www.ncei.noaa.gov/pub/data/igra/history/igra2-metadata.txt'
            urllib.request.urlretrieve(url, 'igra2-metadata.txt')


        names = ['IGRAID', 'WMOID', "NAME", "NAMFLAG" , "LATITUDE", "LATFLAG" , "LONGITUDE", "LONFLAG", "ELEVATION", "ELVFLAG", "YEAR",
                   "MONTH", "DAY", "HOUR","DATEIND","EVENT","ALTIND","BEFINFO","BEFFLAG","LINK","AFTINFO","AFTFLAG","REFERENCE","COMMENT","UPDCOM","UPDATE"] 

        colspecs = [(0,11),(12,17),(18,48),(49,50),(50,60),(61,62),(63,72),(73,74),(75,81),(82,83),(84,88),
                   (89,91),(92,94),(95,97),(98,99),(100,119),(120,122),(123,162),(163,164),(165,167),(168,207),(208,209),(210,234),(235,314),(315,346),(347,354)]


        self.igra2_meta_df = pd.read_fwf('igra2-metadata.txt', 
                            colspecs=colspecs, names=names,
                            ).astype(str)

    def get_igra_metadata(self,station):
        """ Extract the metadata from the igra df for the given station """ 
        df = self.igra2_meta_df
        
        station = station.split('-')[-1]
        # extracting all WMO ids
        wmos = [i if len(i) == 5 else '0'+i for i in df.WMOID]

        df['WMOID'] = wmos
        wmoid = station.split('-')[-1]
        stat = df.loc[df.WMOID == wmoid]

        # Extracting and converting dates 
        
        #print(stat)
        #months = np.unique(stat.MONTH)
        #print('MONTHS :::' , months)
        
        month = [i if len(i) == 2 else '0'+i for i in stat.MONTH]
        month = [m if int(m) <=12 else '01' for m in month]
        
        stat["MONTH"] = month
        stat["DATE"] = stat["YEAR"].astype(str) + stat["MONTH"].astype(str)
        stat["DATE"] = pd.to_datetime(stat['DATE'] , format='%Y%m' )

        
        # create a combined string with the most relevant information 
        update = stat[[ "EVENT","ALTIND","BEFINFO","BEFFLAG","LINK","AFTINFO","AFTFLAG" ]].agg(' , '.join, axis=1)

        u = [ ','.join(  [ v for v in c.split(',') if 'nan' not in v ]) for c in update  ]
        stat['UPDATE'] = u
        stat = stat[ ["DATE", "WMOID", "UPDATE", "REFERENCE", "COMMENT"] ]  # complete data for the station
        
        stat_igra2 = stat 
        # cleaning the station dataframe
        stat_igra2['value'] = 3
        stat_igra2['date_time'] = stat['DATE']
        stat_igra2['comment'] = stat['UPDATE']

        stat_igra2['sensor_id'] = 'IGRA2'
        stat_igra2['source'] = 'IGRA2'

        # Select only IGRA2 metadata relative to "SONDE" events 

        updates = list(stat_igra2.comment.values)
        ind = [ updates.index(i) for i in updates if 'SONDE' in i ]
        stat_igra2_sonde = stat_igra2.iloc[ind]
        
        stat_igra2 = stat_igra2[['value', 'date_time', 'comment', 'sensor_id', 'source']]
        stat_igra2_sonde = stat_igra2_sonde [['value', 'date_time', 'comment', 'sensor_id', 'source']]
        
        return stat_igra2, stat_igra2_sonde
    




class Sensor():
    """ Hold mehtods to extract the sensor_configuration table """ 
    
    def __init__(self):
        # sensor configuration
        sensor_conf = pd.read_csv('sensor_configuration_all.csv', sep='\t',  index_col=0)

        # add a converted column (to strings)
        sensor_id_s = []
        for s in sensor_conf.sensor_id.values:
            try:
                s = eval(s).decode('utf-8').replace(' ','')
            except:
                pass
            s = str(s)
            sensor_id_s.append(s)

        sensor_conf['sensor_id'] = sensor_id_s
        
        self.sensor_conf = sensor_conf
        
        
    def get_sensor_id_comments(self,sensor_id):
        """ Extracts the metadata realtive to a given sensor id from the sensor_configuration table """

        sensor_conf = self.sensor_conf
        s = sensor_id

        if s == 'NA':
            return 'NA'
        # list placeholders

        d = sensor_conf[sensor_conf['sensor_id'] == s ]
        if d.empty:
            s = s.replace('.0', '').replace('.', '')
            if len(s) == 2 and int(s) != 80:
                s = '1' + s
            elif len(s) ==2 and int(s) == 80:
                s = '80'
            d = sensor_conf[sensor_conf['sensor_id'] == s ]

        try:
            com = d.comments.values[0]
        except:
            com = 'NA'        

        try:
            com = eval(com).decode('utf-8')
        except:
            pass

        return com
    

class Analyze():
    
    def __init__(self,Sensor, merged, station):
        
        self.station = station
        self.Sensor = Sensor  # Sensor is a class 
        self.merged = merged
    
    def load_data(self):
        """ Load the data if existing or tries to read it from merged files """
        station = self.station
        
        lista = [f for f in os.listdir('data/') if station in f]
        #print('LISTA ::: ' , lista)
        if not (len(lista)>0):
            
            print("Retrieving data from merged file")
            merged = self.merged
            file = [f for f in os.listdir(merged) if station in f  and 'before' not in f ][0]

            station = file.split('/')[-1].split('_')[0]
            file = merged + '/' + file 

            f = h5.File(file, 'r')
            ts = f['recordtimestamp'][:]
            tsd = pd.to_datetime( ts, unit='s',  origin=pd.Timestamp('1900-01-01') )

            #index_minus = np.where(tsd <=  pd.Timestamp('1994-01-01')  )[0][-1]
            index_minus = 0   # change to start from a certain date onwards 

            #index_plus = np.where(tsd >  pd.Timestamp('1997-01-01')  )[0][0]
            index_plus = np.where(tsd <  pd.Timestamp('2013-01-01')  )[0][-1]

            ### Extracting Schroeder 
            ind_obs_sch = list(f['recordindex'][:]) [index_minus:index_plus]
            i = np.take( f['observations_table']['sensor_id'][:].view('|S4') , ind_obs_sch) 
            ids_s = [s.decode('utf-8').replace('.0','').replace('.','') for s in i ]
            dic = {'date_time': tsd[index_minus:index_plus] , 'sensor_id': ids_s }

            data_sch = pd.DataFrame(dic)
            data_sch['value'] = 1

            ### Extracting WMO
            ind_obs_wmo     = list(f['recordindex'][:]) [index_plus+1:]
            ind_obs_wmo_all = list(f['recordindex'][:]) # taking all WMOs


            wmoids = np.take(  f['observations_table']['sensor_id'][:].view('|S4') , ind_obs_wmo)
            wmoids = [s.decode('utf-8') for s in wmoids ]

            dic_wmo = {'date_time':tsd[index_plus+1:] , 'sensor_id':wmoids }
            data_wmo = pd.DataFrame (dic_wmo)
            data_wmo['value'] = 2

            data_wmo.to_csv('data/' + station + '_wmo.csv' , sep = '\t') 
            data_sch.to_csv('data/' + station + '_sch.csv' , sep = '\t') 

            f.close()

        else:
            print("Loading existing data")
            
            data_wmo = pd.read_csv( 'data/' + [f for f in lista if 'wmo' in f][0] , sep = '\t')
            data_wmo['date_time'] = pd.to_datetime(data_wmo['date_time'] )

            data_sch = pd.read_csv('data/' + [f for f in lista if 'sch' in f][0]  , sep = '\t') 
            data_sch['date_time'] = pd.to_datetime(data_sch['date_time'] )
        
    
        data_wmo['source'] = 'WMO'
        data_sch['source'] = 'SCH'
    
        # tidying some values 
        data_wmo['comment'] = [ str( self.Sensor.get_sensor_id_comments(str(i).replace(' ','').replace('.0',''))) for i in data_wmo.sensor_id]
        data_sch['comment'] = [ str( self.Sensor.get_sensor_id_comments(str(i).replace(' ','').replace('.0',''))) for i in data_sch.sensor_id]
            
        data_wmo['sensor_id'] = data_wmo['sensor_id'].astype(str)
        self.data_sch = data_sch
        self.data_wmo = data_wmo 
        
        
    def get_indices(self, data):
        """ # find the indices where the sensor was replaced 
        i.e. spots the change in the sensor ina  time series """

        data = data.reset_index()
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
                    last = sid
                    indices.append(index)
        return indices
 

    def clean_df(self, df):
        """ Clean the WMO dataframe from all nans """
        
        # cleaning WMO data from nans 
        data_wmo_clean = df.loc[ (self.data_wmo.sensor_id != 'nan') 
                               & (self.data_wmo.sensor_id != '-922')
                               & (self.data_wmo.sensor_id != -922)].dropna( subset=['sensor_id'] )                                                                          
                                                                                                                 
        data_wmo_clean.reset_index()

        return data_wmo_clean 
    
    
    
    def analyze(self):
        """ Extract information from the station file dataframe """
        
        self.load_data()
        
        # data from Schroeder and WMO 
        wmo = self.data_wmo
        sch = self.data_sch

        
        # cleaning the df from nans and nans values
        # need to be separated for WMO bar plot
        data_wmo_clean = self.clean_df(wmo)
        self.clean_data_wmo = data_wmo_clean
        indices_wmo_clean = self.get_indices(data_wmo_clean)
        
        #  self.get_indices simplifies the df by taking into accounts only updates in the sensors
        self.data_wmo_clean_simplified = data_wmo_clean.iloc [ indices_wmo_clean ] 

        # simmplified indices
        indices_sch = self.get_indices(sch)
        indices_wmo = self.get_indices(wmo)
        
        
        # data excluding nans, simplified sensors updates
        data_clean_simplified = pd.concat( [sch.iloc[ indices_sch], 
                                            data_wmo_clean. iloc[ indices_wmo_clean] ] ) 
                    

        self.data_all_cleaned_simplified = data_clean_simplified
        

    def get_all_sensors(self, df):
        """ Extract a small table with unqiue sensors and description """
        
        igra2 = df.loc[df['source'] == 'IGRA2']
        
        rest = df.loc[df['source'] != 'IGRA2']
        sensors,ind = np.unique( rest.sensor_id, return_index= True)
        df_sensor = rest.iloc[list(ind)] [['sensor_id', 'source', 'comment']]
        
        df_sensor = pd.concat([df_sensor, igra2])
        return df_sensor


    def simplify_day_night(self, df):
        df.reset_index()
        hours = pd.to_datetime( df.date_time).dt.hour
        moments = []
        for h in hours:
            h = int(h)
            if h <= 3 or h > 21:
                m = 'night'
            elif h > 3 and h <= 9:
                m = 'morning'
            elif h > 9 and h <=15:
                m = 'day'
            elif h > 15 and h <=21:
                m = 'evening'
                #print("EVENINGGGGG" , h )
            moments.append(m)

        df['moment'] = moments

        window = 10
        #dfr=df[window+1:]

        # fill vector with indices to keep
        indices_to_keep = list(range(window))

        for index in range(window, len(df)) :

            # idea: take "windows" previous records, extract only same moment of the day,
            # get list of unique sensor_ids
            # if list == one single id equal to the considered one, then do not save index

            dfr =  df[index-window:index]
            m = df.iloc[index].moment 

            sensor =  df.iloc[index].sensor_id

            dfr = dfr[ dfr.moment == m ]

            dfr_sensor = list( dfr[dfr.moment== m].sensor_id.values) 
            if len(dfr_sensor) ==0:
                indices_to_keep.append(index)

            else:

                if sensor == dfr_sensor[-1]:
                    continue

                else:
                    indices_to_keep.append(index)

        if len(indices_to_keep) >0:
            return df.iloc[indices_to_keep] 
        else:
            return df 
    
# Loading IGRA2 metadata 
ig = IgraMetaData()
#igra2_metadata = ig.igra2_meta_df
    
# Loading sensor configuration
sensor = Sensor()

# Merged file source (if data not already available)
merged = '/scratch/das/federico/MERGED_APRIL2022'

def get_data(station, force_create=False):
    """ Extract the data for plotting either from the database or from the reduced csv files 
    stored in data_plots directory """
    
    out_dir_data_plots = 'data_plots'
    
    if not os.path.isdir(out_dir_data_plots):
        os.mkdir(out_dir_data_plots)
        
    analyze = Analyze(sensor,merged,station)
            
    if not (os.path.isfile(out_dir_data_plots + '/' + station + '_data_clean_all.csv') 
            and os.path.isfile(out_dir_data_plots + '/' + station + '_all_sensor_station.csv')
            and os.path.isfile(out_dir_data_plots + '/' + station + '_clean_data_wmo.csv' )
            ) or force_create:

        try:
            all_stat = os.listdir(merged)
            all_stat = [s.split('_')[0] for s in all_stat ]
        except:
            logging.error("Cannot retrieve data! Please check accessibility to merged directory! ")
            sys.exit()
            
        if force_create:
            logging.debug(" --- FORCING CREATION --- data file: ")

        logging.debug(" --- RETRIEVING --- data file: ")

        stat_igra2, stat_igra2_sonde = ig.get_igra_metadata(station)

        # Analyze data
        logging.debug(" --- ANALYZING --- data file: ")

        
        analyze.analyze()
        
        data_all_cleaned_simplified = analyze.data_all_cleaned_simplified
        clean_data_wmo = analyze.clean_data_wmo
        
        # data_sch, data_wmo, data_df, data_wmo_clean, data_df_clean = analyze.analyze()
        
        # data_sch, data_wmo, data_df, data_wmo_clean_red, data_df_clean, data_wmo_clean = analyze.analyze()
        
        # concatenating WMO, Sch and IGRA2
        
        data_all_cleaned_simplified = analyze.simplify_day_night(data_all_cleaned_simplified)
        
        data_all_cleaned_simplified = pd.concat([data_all_cleaned_simplified, stat_igra2_sonde])
        
        # extract unique sensor id table for the station
        all_sensor_station = analyze.get_all_sensors(data_all_cleaned_simplified)

        all_sensor_station.to_csv(out_dir_data_plots          + '/' + station + '_all_sensor_station.csv' , sep='\t')
        data_all_cleaned_simplified.to_csv(out_dir_data_plots + '/' + station + '_data_clean_all.csv' , sep='\t')
        clean_data_wmo.to_csv(out_dir_data_plots              + '/' + station + '_clean_data_wmo.csv' , sep='\t')

    else:
        logging.debug(" --- READING --- data file: ")
        data_all_cleaned_simplified = pd.read_csv(out_dir_data_plots     + '/' + station + '_data_clean_all.csv', sep='\t')
        all_sensor_station = pd.read_csv(out_dir_data_plots + '/' + station + '_all_sensor_station.csv', sep='\t')
        clean_data_wmo = pd.read_csv(out_dir_data_plots     + '/' + station + '_clean_data_wmo.csv', sep='\t')
            
    
    #print(data_df_clean_all.head(10) , all_sensor_station_df.head(10) )
    all_sensor_station = all_sensor_station[ ['sensor_id', 'source', 'comment'] ] 
    all_sensor_station['comment'] = [str(s).replace('\n','').replace('\t','').replace('^', ' ').replace('>',']').replace('<','[')
                                     for s in all_sensor_station.comment]
    

    data_all_cleaned_simplified = data_all_cleaned_simplified [ [c for c in data_all_cleaned_simplified.columns if "Unnamed" not in c ] ]
    clean_data_wmo = clean_data_wmo [ [c for c in clean_data_wmo.columns if "Unnamed" not in c ] ]
    
    
    
    return [data_all_cleaned_simplified , all_sensor_station, clean_data_wmo ]

        





