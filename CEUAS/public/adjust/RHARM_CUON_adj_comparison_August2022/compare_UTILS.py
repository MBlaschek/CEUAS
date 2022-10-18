import numpy as np
import pandas as pd
import sys, glob
import urllib3
import h5py
import cdsapi, zipfile, os, time
import warnings
import shutil
import xarray
from datetime import date
warnings.filterwarnings('ignore')
import pycountry
sys.path.append(os.getcwd()+'/../../cds-backend/code/')
import cds_eua3 as eua
# import numbaprocess
import copy
import glob
from numba import njit
import pandas
import glob
import pickle
import matplotlib
import matplotlib.pyplot as plt
import h5py as h5
import netCDF4 as nc
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from multiprocessing import Pool
from functools  import partial






# source of the igra dataset
igra = 'insitu-observations-igra-baseline-network'

# cuon base dir
basedir = '/mnt/users/scratch/leo/scratch/converted_v8'

# SUNY dir
suny_dir = '/users/staff/leo/fastscratch/SUNY/homo-raw-subdaily-station/'


def request(rqdict, source, cds='_test', remove_file=True):
    """
    Run a request on the CDS server, and downloads the data
    
    Remember to UPDATE the file .cdsapirc on the server home directory with the CORRECT cds address!
    See e.g. https://cds-test.copernicus-climate.eu/cdsapp#!/dataset/insitu-observations-igra-baseline-network?tab=form 
    """
    if cds != '_test':
        cds = ''
    t0 = time.time()

    c = cdsapi.Client()
    print(c)
    r = c.retrieve(
        source, rqdict)
    print('Request took: ' + str(time.time() - t0) + ' seconds')
    if True:
        r.download(target='download.zip')
        assert os.stat('download.zip').st_size == r.content_length, "Downloaded file is incomplete"
    z = zipfile.ZipFile('download.zip')
    z.extractall(path='./download' + cds + '/' )
    z.close()
    return


# Fill the request form 
"""
out = request({
    'source': 'IGRA_H',# 'IGRA'
    'variable':['air_temperature'],
    'period':'2009-01-01/2010-12-31',
#     'format': 'csv-lev.zip', # no format = .nc
    'station_name':'AUM00011035',
},igra, remove_file=False)
"""


### TO DO -> might not be necessary anymore, it was just an additional test
# Here: select the name of the download directory according  to TEST or PUB
# You must change the API file .cdsapi accordingly before the download !!!

cds  = '_test'
cds = ''


def get_RHARM_IGRA_from_CDS(stations):
    """ Download HARM/IGRA data from the CDS (from both the TEST and PUBLIC cds)
        save csv files
        converts to pandas
        store in dictionaries """
    
    # container dictionaries for all the downloaded data in the form of pandas dataframe
    all_HARM_df = {} 
    all_IGRA_df = {}

    all_HARM_df_test = {}
    all_IGRA_df_test = {}

    stdplevs = [ i*100 for i in [10.0, 20.0, 30.0, 50.0, 70.0, 100.0, 150.0, 200.0, 250.0, 300.0, 400.0, 500.0, 700.0, 850.0, 925.0, 1000.0] ]

    # Download both public CDS and TEST RHARM
    # Loop over the station list, and download data from IGRA and RHARM as csv file.
    for cds in ['_test' , '']: 

        out_dir =  'download' + cds  ### CHECK that this line is here and not defined outside loop !!!

        # copy correct cdsapirc file 
        if 'test' in cds:
            os.system ( 'cp /users/staff/federico/cdsapirc_TEST.txt  /users/staff/federico/.cdsapirc' )
        else:
            os.system ('cp /users/staff/federico/cdsapirc_PUB.txt  /users/staff/federico/.cdsapirc' )

        for s in stations:
            out = request({
            'source': 'IGRA_H',# 'IGRA'
            'variable':['air_temperature'],
            'period':'2009-01-01/2010-12-31',
            'format': 'csv-lev.zip',
            'station_name': s,
            'archive_type': 'harmonized_global_radiosonde_archive',

            },igra, cds=cds, remove_file=False)


            files = glob.glob(out_dir+'/*.csv')


            file_dest = out_dir + '/' + s + '_RHARM_data.csv'
            for f in files:
                if '.csv' in f and 'data' not in f:
                    print( f, '        ' ,  file_dest)
                    os.rename( f, file_dest )

            # reading the data as pandas dataframe
            df = pd.read_csv(file_dest, skiprows=14)
            # storing only standard pressure levels
            df = df[df.air_pressure.isin(stdplevs ) ]
            df =df.sort_values(by=['report_timestamp'])
            # da = eua.CDMDataset('download/USM00070308_IGRA_data.nc').to_dataframe() ### example
            all_HARM_df[s] = df 

            if 'test' in cds:
                all_HARM_df_test[s] = df 

            else:
                all_HARM_df[s] = df 

            out = request({
            'source': 'IGRA',# 'IGRA'
            'variable':['air_temperature'],
            'period':'2009-01-01/2010-12-31',
            'format': 'csv-lev.zip',
            'station_name': s,
            'archive_type': 'global_radiosonde_archive',

            },igra, cds=cds, remove_file=False)

            filess = glob.glob(out_dir+'/*.csv')
            file_des = out_dir + '/' + s + '_IGRA_data.csv'
            for ff in filess:
                if '.csv' in ff and 'data' not in ff and 'HARM' not in ff :
                    os.rename( ff, file_des )

            dd = pd.read_csv(file_des, skiprows=14)
            dd = dd[dd.air_pressure.isin(stdplevs ) ] 
            dd = dd.sort_values(by=['report_timestamp'])

            if 'test' in cds:
                all_IGRA_df_test[s] = dd 

            else:
                all_IGRA_df[s] = dd       
                
    return all_HARM_df, all_IGRA_df, all_HARM_df_test, all_IGRA_df_test 





def get_HARM_IGRA_from_CSV(file):
    """ Read data from CSV (provided by Emanuele for few stations """
    
    df = pd.read_csv(file, sep = ',')
    
    return df 



def extract_CUON(file: str='' , year_filter: str=""):
    
    """ Extract CUON data from local database.
    Only standard plevels are selected. 
    year_filter : string with the minim year to consider e.g. '1987' --> TO DO IMPLEMENT
    """
    
    # loop over groups and variables
    var = {'observations_table':['z_coordinate' , 'date_time', 'observation_value', 'latitude' , 'longitude'], 
       'advanced_homogenisation': ['RAOBCORE_bias_estimate', 'RASE_bias_estimate', 'RICH_bias_estimate', 'RISE_bias_estimate'],
       'era5fb': ['biascorr@body', 'biascorr_fg@body']
      }


    f = h5.File(file, 'r')
    
    plev_indices = np.isin( f['observations_table']['z_coordinate'][:], stdplevs ) 
    ind_min, ind_max = f['recordindices']['126'][0], f['recordindices']['126'][-1]

    # placeholder dic
    d = {}
    
    # remove some data before 2009
    #year_filter = 109*365*60*60*24 + 1900
    #a = np.where( f['observations_table']['date_time'][ind_min:ind_max] > year_filter )[0][0]

    # boolean mask for plevels
    plev_indices = np.isin( f['observations_table']['z_coordinate'][ind_min:ind_max], stdplevs )

    for g in var.keys():
        for v in var[g]:
            d[v] = f[g][v][ind_min:ind_max][plev_indices] 
        
    
    df = pd.DataFrame.from_dict(d)
    df['date_time']  = pd.to_datetime( df['date_time'][:], unit='s',  origin=pd.Timestamp('1900-01-01') )
    
    if not os.path.isdir('data_CUON'):
        os.makedir('data_CUON')
        
    #df = df.loc [ (df.date_time >= min_date ) & (df.date_time <= max_date)  ]
    
    print(file.split('/')[-1].split('_')[0].split('-')[-1])
    df.to_csv( 'data_CUON/' + file.split('/')[-1].split('_')[0].split('-')[-1] + '.csv',  sep='\t' )
    return df




def extract_SUNY(file):
    """ Extract data from SUNY files,
    converting to a sorted pandas dataframe """
    
    # f['pressure'][:] -> [  1000,   2000,   3000,   5000,   7000,  10000,  15000, 20000,  25000,  30000,  40000,  50000,  70000,  85000, 92500, 100000, 101300]
    # f['UTC'][:] -> [0,12]
    f = nc.Dataset(file)
    pressure = f['pressure'][:]
    dic = { 'date_time':[], 'air_pressure':[], 'air_temperature':[], 'air_temperature_adj':[] , 'hour':[] }
    
    suny_dir = '/users/staff/leo/fastscratch/SUNY/homo-raw-subdaily-station/'

    for p in range(len(pressure)):  ### p is the index for the pressure value 
        for t in range(2): ### t is the index for the time value

            #print(len(f['time'][:].data) )
            #date_time = format='%Y%m%d'
            #print(f['time'][:].data[d] , ' ' , f['UTC'][:].data[t]  )
            hour = '00' if t==0 else '12'
            #dt = pd.to_datetime( str(f['time'][:].data[d]) + " " + hour , format='%Y%m%d %H' ) 

            dic['date_time']          .extend( f['time'][:].data  )
            #print(p, ' ', t ,  f['rawT'][:].data[:,p,t] )
            dic['air_temperature']    .extend( f['rawT'][:].data[:,p,t] )
            dic['air_temperature_adj'].extend( f['homoT'][:].data[:,p,t])

            dic['air_pressure']       .extend( len(f['time'][:].data ) * [f['pressure'][:].data[p]]  )
            
            hour = '00' if t==0 else '12'
            dic['hour']               .extend( len(f['time'][:].data ) * [ hour ] )
    
    df = pd.DataFrame.from_dict( dic ) #.sort_values( by=['date_time' , 'air_pressure'])
    dt = pd.to_datetime(  df['date_time'][:].astype(str)+ " " +  df['hour'][:].astype(str) , format='%Y%m%d %H' )
    
    df['date_time'] = dt
    df = df.sort_values( by=['date_time' , 'air_pressure'])

    if not os.path.isdir('data_SUNY'):
        os.mkdir('data_SUNY')
        


    df = df.loc[df["date_time"] <= pd.Timestamp('20181231') ]

        
    # correct original temp, which is in tenth of Celsius
    df['air_temperature'] = df['air_temperature']/10 + 273.15 
    df['air_temperature_adj'] = df['air_temperature_adj']/10 + 273.15 
    

    
    df.to_csv('data_SUNY/' + file.split('/')[-1].replace('.nc','') + '.csv' , sep='\t' )
    
    return df




### bug with timestamps https://github.com/matplotlib/matplotlib/issues/18158/  with matplotlib -> use plotly
def get_plevel_data_from_station(station, all_RHARM, all_CUON, all_SUNY, p, source= 'csv', cds='test'):
    
    """ Extract the data from the IGRA, RHARM, CUON dataframes for plotting.
        Adjustments of RHARM must be calculated by subtracting obs data RHARM-IGRA
        """
    
    # TO DO test switch NOT IMPLEMENTED
    if source != 'csv':
        return 0
        """
        if 'test' in cds:
            igra = all_IGRA_df_test
            rharm = all_HARM_df_test
        else:
            igra = all_IGRA_df
            rharm = all_HARM_df

        i = igra[station]
        igra_p = i[i.air_pressure==p]
        #igra_p = igra_p.reset_index()

        h = rharm[station]
        harm_p = h[h.air_pressure==p]
        """
    else:
        i = all_RHARM[station]
        igra_p = i[i.air_pressure==p]
        harm_p = igra_p
        
    c = all_CUON[station]
    cuon_p = c[c.z_coordinate == p] 
    
    s = all_SUNY[station]
    #print(s.head)
    suny_p = s[s.air_pressure == p]
    
    return igra_p, harm_p, cuon_p, suny_p



def get_rharm_adj(i,h,):
    """ Extract the RHARM adjustment.
        First join the dataframe along the timestamp-pressure axes, then calculated the difference between 
        the IGRA-RHARM temp. values.
        With this convention we consider [Adjusted = Raw - Bias] 
        Return the RHARM dataframe with an additional column """
    
    i = i[['report_timestamp', 'air_pressure', 'air_temperature']]
    i['air_temperature_i'] = i['air_temperature']
    
    h = h [['report_timestamp','air_pressure', 'air_temperature']]
    h['air_temperature_h'] = h['air_temperature']


    d = h.merge(i, how='outer')
    d['adj'] = d.air_temperature_i - d.air_temperature_h
    

        
    return d



def get_night_day_diff(df, what='IGRA', source=''):
    
    """ Extract day and night differences from a dataframe.
        Processing slightly different for IGRA/RHARM, CUON or SUNY.
        
        Have to return two different temp_diff since for RHAR/IGRA from CDS I process separately the two dataframes
        and then calculate the adj by subtracting RHARM-CUON,
         
        while for CUON and  RHAR/IGRA crom CSV I have the adj in the same dataframe """
    
    #print(df.columns, ' ' , what)
    # name of timestamp in the different dataframes 
    dt_col = 'date_time' if what != 'IGRA' else 'report_timestamp'

    df['date'] =  pd.to_datetime(df[dt_col])
    df['hour'] = df['date'].dt.hour
    
    df = df.reset_index()

    # Must correct for days 
    
    ### print( df.head(10) )
    night_hours = [0,1,2,3 ]
    nights = list(np.where( df.hour.isin(night_hours)  )[0] )
    df.loc[nights, 'hour'] = 0
    # ---    
    following_day_night_hours = [21,22,23]
    nights = list(np.where( df.hour.isin(following_day_night_hours))[0])
    df['date'][nights] = df['date'][nights] +  pd.Timedelta(1,unit='d')
    df.loc[nights, 'hour'] = 0
    # --- 
    day_hours = [11,12,13,14] 
    day = np.where( df.hour.isin(day_hours))[0]
    df.loc[day, 'hour'] = 12
    
    df['date'] = df['date'].dt.date

        
    iday   = df.loc[df.hour == 12 ]
    inight = df.loc[df.hour == 0 ]

    #iday = iday.drop_duplicates(subset=['date'])
    #inight = inight.drop_duplicates(subset=['date'])

    
    # df['adj'] = df.air_temperature_rharm - df.air_temperature_igra 

    days_set = set(iday.date) # make a set of the first list
    days_set.intersection_update(inight.date)
    #print(len(iday.date), len(inight.date), len(days_set) )
    
    day_set = list(days_set)
    df = df.loc[df['date'].isin(days_set)]
    
    #return day_set, df 

    days_set = df.date
    #print(df)
    #print(df[df.hour==12].head(10).air_temperature_igra)
    #print(df[df.hour==0].head(10).air_temperature_igra)
    #print(  df[df.hour==12].air_temperature_igra.values - df[df.hour==0].air_temperature_igra.values)

    df = df.drop_duplicates( subset =['date','hour'])
    
    days = df[df.hour==12].date 
    
    if what=='IGRA':
        if source != 'csv': # two deparates df for IGRA and RHARM, downloaded from CDS ### TO CHANGE, overly complicated :-(
            temp_diff =  df[df.hour==12].air_temperature.values  -  df[df.hour==0].air_temperature.values
            return days, temp_diff, temp_diff # 3rd value is dummy 
    
        else:
            # data taken from emanuele CSV 
            temp_diff_i =  df[df.hour==12].air_temperature_igra.values   - df[df.hour==0].air_temperature_igra.values
            temp_diff_h =  df[df.hour==12].air_temperature_rharm.values - df[df.hour==0].air_temperature_rharm.values
            return  days, temp_diff_i , temp_diff_h
        
    elif what == 'SUNY':
        temp_diff_raw =  df[df.hour==12].air_temperature.values  -  df[df.hour==0].air_temperature.values
        temp_diff_adj =  df[df.hour==12].air_temperature_adj.values  -  df[df.hour==0].air_temperature_adj.values

        return days, temp_diff_raw, temp_diff_adj
        
    else:
        temp_diff =  df[df.hour==12].observation_value.values  - df[df.hour==0].observation_value.values
        temp_diff_adj =  df[df.hour==12].observation_value.values - df[df.hour==12].RAOBCORE_bias_estimate.values
        temp_diff_adj = temp_diff_adj -(df[df.hour==0].observation_value.values-df[df.hour==0].RAOBCORE_bias_estimate.values)  

        
        return days, temp_diff, temp_diff_adj
    
    
    
    

    
def extract_data_from_df(station,  all_RHARM, all_CUON, all_SUNY, p, source='csv'):
    """ Extract and manipulate all data from the original dataframes for RHARM,IGRA,CUON,SUNY """
    
    # extracting data for plotting (IGRA, RHARM, CUON)
    # store data in dictionary
    data = {'CUON':{} , 'IGRA':{} , "RHARM":{} , "SUNY":{} }

        
    ##############  Extracting single plevel Igra, rHarm, Cuon, Suny
    i,h,c,s           = get_plevel_data_from_station(station, all_RHARM, all_CUON, all_SUNY, p, cds='') 
    data["CUON"]['df']  = c
    data["IGRA"]['df']  = i
    data["RHARM"]['df'] = h
    data["SUNY"]['df']  = s
    
    ############## calculating day-night temp difference
    if source != 'csv':
        i_date, i_tempdiff, dummy = get_night_day_diff(i)  # IGRA from CDS
        h_date, h_tempdiff, dummy = get_night_day_diff(h)  # RHARM from CDS
        data["IGRA"]['daynight']  = [i_date, i_tempdiff,  []] 
        data["RHARM"]['daynight'] = [h_date, h_tempdiff,  []] 
        
    else:
        ih_dates, temp_diff_i , temp_diff_h = get_night_day_diff(i, what='IGRA', source='csv')
        data["IGRA"]['daynight']  = [ih_dates, temp_diff_i , temp_diff_h] 
        data["RHARM"]['daynight'] =  data["IGRA"]['daynight']  # same, dummy 
        
    c_date, c_temp_diff, c_temp_diff_RAOBCORE = get_night_day_diff(c, what='CUON')
    data["CUON"]['daynight'] = [c_date, c_temp_diff, c_temp_diff_RAOBCORE]     

    s_date, s_tempdiff_raw, s_tempdiff_adj = get_night_day_diff(s, what="SUNY")  # RHARM from CDS    
    data["SUNY"]['daynight'] = [s_date, s_tempdiff_raw, s_tempdiff_adj]
    
    return data


def filter_day_night(i,h,c,s, time = 'day'):
    # filtering hours for day and night
    
    if time == 'day':
        i = i.loc[ i['report_timestamp'].dt.hour == 12 ] 
        h = h.loc[ h['report_timestamp'].dt.hour == 12 ] 
        c = c[c['date_time'].dt.hour.isin([10,11,12,13,14]) ] 
        s = s[s['date_time'].dt.hour == 12 ]
        
    elif time == 'night':
        i = i.loc[ i['report_timestamp'].dt.hour== 0 ] 
        h = h.loc[ h['report_timestamp'].dt.hour== 0 ] 
        c = c[c['date_time'].dt.hour.isin([22,23,0,1,2]) ] 
        s = s[s['date_time'].dt.hour == 0 ]
        
    return i,h,c,s


def calculate_day_night_adjustment(data='', time='night', source=''):
    """ Calculate adjustments """
    
    c,i,h,s = data["CUON"]['df'], data["IGRA"]['df'],  data["RHARM"]['df'], data["SUNY"]['df']
    
    ############## calculating adjustments
    i,h,c,s = filter_day_night(i,h,c,s, time=time)

    if source != 'csv':
        adj      = get_rharm_adj(i,h)
        #adj_test = get_rharm_adj( i_test,h_test)  ### this was for CDS-TEST, remove 
    else:
        adj = i
        #adj_test = i ### this was for CDS-TEST, remove 
        
    # return df ready for adjustments
    return c, s, i, adj     
        
    
    
    