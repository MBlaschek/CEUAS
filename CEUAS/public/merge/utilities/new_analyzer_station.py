
import os,sys
import netCDF4 as nc
import xarray as xr
import pandas as pd
import time
import numpy as np
from tqdm import tqdm
from collections import Counter
import datetime
import h5py as h5

import matplotlib
from matplotlib import ticker

#matplotlib.use('Agg')


import matplotlib.pylab as plt


""" Setting pandas printing output """
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)


""" Output Directory """
out_dir = os.getcwd() + '/Plots'
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)    
    
""" Class holding general properties of the plots """
class Common():
    def __init__(self):
        
        """ For plotting style """
        style_dic = {126 : { 'l': 'Air Temperature [K]'         , 'c': 'blue'      } , 
                     117: { 'l': 'Geopotential [m]'            , 'c': 'green'     } ,
                     139: { 'l': 'Wind u-component [m/s]'      , 'c': 'gold'      } ,
                     140: { 'l': 'Wind v-component [m/s]'      , 'c': 'limegreen' } ,
                     107: { 'l': 'Wind speed [m/s]'            , 'c': 'black'     } ,
                     106: { 'l': 'Wind from direction'         , 'c': 'cyan'      } ,    
                    }
    
        datasets_dic = {'era5_1'     : { 'l': 'ERA5 1'    , 'c': 'blue'      } ,
                        'era5_2'     : { 'l': 'ERA5 2'    , 'c': 'limegreen' } ,
                        'ncar'       : { 'l': 'NCAR'      , 'c': 'purple'    } ,
                        'igra2'      : { 'l': 'IGRA2'     , 'c': 'cyan'      } ,
                        'era5_1759'  : { 'l': 'ERA5 1759' , 'c': 'gray'      } ,
                        'era5_1761'  : { 'l': 'ERA5 1761' , 'c': 'pink'      } ,
                        'era5_3188'  : { 'l': 'ERA5 3188' , 'c': 'yellow'    } ,
                        'bufr'       : { 'l': 'BUFR'      , 'c': 'green'     } ,
                    }
        
        self.style_dic    = style_dic
        self.datasets_dic = datasets_dic 
        self.std_plevs    = [10, 20, 30, 50, 70, 100, 150, 200, 250, 300, 400, 500, 700, 850, 925, 1000]
        self.fontsize     = 15
        self.datasets     = ["era5_1" , "era5_2" , "era5_1759" , "era5_1761" , "era5_3188" , "bufr" , "ncar" , "igra2" ]
        
Common = Common()

""" Standard pressure levels in hPa """
std_plevs = Common.std_plevs




def analyze_headers(data= '', station = '', date_min='', date_max='', text = '' , text_save='', sources='' ):
        
    """ Extracts some information from the header table """
    timestamps_unconverted = data["record_timestamp_unconverted"]
    timestamps = data["record_timestamp"]
    
    #duplicates = data["duplicate"]
        
    datasets_dic = Common.datasets_dic
    style_dic    = Common.style_dic
    fs = Common.fontsize
    
    
    """ Plots and compare the distributions of time differences of subsequent records 
    (to check if delta time requirement during merging is respected) """
    def plot_timeshifts_distribution(timestamps):
        
        diff = []
        for dt,dtp in zip(timestamps, timestamps[1:]):
            try:
                d = (dtp - dt) / 3600.
                if d > 6:
                    continue
                diff.append(d)
            except:
                pass
        
    
        fig, ax = plt.subplots(figsize=(12,10) ) #    fig, axs = plt.subplots(nrows=2, ncols=3, constrained_layout=True , figsize=(15,10) )

        nbins = 24
        plt.hist(diff, nbins, histtype='bar', color = 'magenta')
        plt.xlim(0, 6)
        ax.grid(ls =":" , color = "lightgray")

        plt.xticks( fontsize = fs  )
        plt.ylabel("Record Counts", fontsize = fs )
        plt.xlabel("Time interval [bin = 1/4 hour]", fontsize = fs )


        plt.savefig(out_dir + '/'  + station + "_timeshifts_distribution_" + text_save + ".png", bbox_inches = 'tight' , dpi = 250 )     

        plt.show()
        plt.close()
        
        
    def plot_dataset_range(data, date_min='', date_max='' , text= '', station='' , text_save='', sources=''):
        """ Make a plot for the range of date time availabe for each dataset """
        
        sources = list(np.unique(data['source_id']) )
        results = {} 
        for s in sources:
            df_s = data.loc [data['source_id'] ==s ]
    
            results[s] = df_s
    
    
        if date_min:
            fig, ax = plt.subplots(figsize=(20,6) ) #    fig, axs = plt.subplots(nrows=2, ncols=3, constrained_layout=True , figsize=(15,10) )
        else:
            fig, ax = plt.subplots(figsize=(12,10) ) 
            
        fig.suptitle("Time Interval per Data Source for Station " + station + ' ' + text, y = 0.94, fontsize = fs )

        ticks, labels = [], [] 
        
        num_sources = len(results.keys() )        
        for r,i in zip(results.keys() , range(num_sources)):
            index = i+1
            ticks.append(index)
            label = Common.datasets_dic[r]['l']
            labels.append(label)
            y = np.empty( len(results[r]) ) 
            y.fill(index)
            plt.scatter(results[r]['record_timestamp'] , y, color = datasets_dic[r]['c'] , label = label )

        ax.grid(ls =":" , color = "lightgray")

        ax.set_yticks(ticks)
        ax.set_yticklabels(labels, fontsize = fs ,)
        plt.xticks( rotation = 45, fontsize = fs  )
        
        text = ''
        if date_min and date_max:
            plt.xlim(date_min, date_max)
            text = '_zoom_'
            
        
        plt.savefig(out_dir + '/'  + station + "_dataset_series_" + text_save + ".png", bbox_inches = 'tight' , dpi = 250 )     
        
        plt.show()
        plt.close()
        
        
        
    def plot_bar( counts = '', labels = '', rotation = False , text = '' , station = '' , colors = '', sources=''):
        """ Creates a simple bar plots """

        formatter = ticker.ScalarFormatter(useMathText=True)
        formatter.set_scientific(True) 
        formatter.set_powerlimits((-1,1)) 
        fig, ax = plt.subplots(figsize=(12,10) )
        ax.yaxis.set_major_formatter(formatter) 
        ax.grid(ls =":" , color = "lightgray")
        X = np.arange(len(counts))
        ax.bar(X, counts, color = colors )
        ax.set_xticks(X)
        ax.set_xticklabels(labels, fontsize = fs )
        if rotation:
            ax.set_xticklabels(labels, rotation = 45, fontsize = fs )

        plt.ylabel("Record Counts", fontsize = fs )
        plt.title('Record counts per data source for Station ' + str(station) + ' ' + text , y=1.02 , fontsize = fs )

        plt.savefig(out_dir + '/'  + station + "_recordcounts_datasets_distribution_global.pdf",
                bbox_inches = 'tight' , dpi = 250 )    
        plt.show()
        plt.close()  
    

    def plot_dataset_distribution(data, date_min='', date_max='', station='', text='' , text_save = 'zoom' , sources=''):
        """ Make a histogram for the distribution of date time availabe for each dataset """
        num_sources = len(data.keys() )
    
        """ Creating the list for the stacked histogram """
        ticks, labels, X, colors, counts = [], [], [], [], [] 
        Min, Max = [], [] 
        for rs in sources :
            d = data.loc[ data['source_id'] == rs]
            dt = d['record_timestamp']
            label = datasets_dic[rs]['l'] + "[" + str(len(dt)) + "]"
            counts.append(len(d['record_timestamp']))
            labels.append(label)
            colors.append(datasets_dic[rs]['c'] )
            X.append(dt)
            Min.append(min(dt))
            Max.append(max(dt))
            
        a = plot_bar( counts = counts, labels = labels, rotation = False , text = '' , station = station , colors = colors, sources=sources)     

        if date_min:
            fig, ax = plt.subplots(figsize=(20,6) ) #    fig, axs = plt.subplots(nrows=2, ncols=3, constrained_layout=True , figsize=(15,10) )
        else:
            fig, ax = plt.subplots(figsize=(15,10) ) 
            
        """ Extracting the min and Max years in the range """
        mm = np.array ( [min(Min)] , dtype = np.datetime64 )
        y_min = str(pd.to_datetime(mm).year[0])
        MM = np.array ( [max(Max)] , dtype = np.datetime64 )
        y_max = str(pd.to_datetime(MM).year[0])
        
        #print("min, Max" , y_min[0], y_max[0] )
        
        nbins = int(y_max) - int(y_min)
        plt.hist(X, nbins, histtype='bar', stacked=True, label=labels, color = colors )
    
        ax.grid(ls =":" , color = "lightgray")

        plt.xticks(rotation = 45 , fontsize = fs  )
        plt.yticks(fontsize = fs  )

        
        if date_min and date_max:
            plt.xlim(date_min, date_max)
    
        plt.ylim(0,4600)
        plt.ylabel("Record counts " , fontsize = fs)
        plt.legend(fontsize = fs  , loc = 'best' )
        plt.title("Records Distribution for station " + station + ' ' + text , y = 1.03, fontsize = fs  )
        plt.savefig(out_dir + '/'  + station + "_record_distribution_" + text_save + ".png", bbox_inches = 'tight' , dpi = 250 )     
        
        plt.show()
        plt.close()
        
    #p = plot_timeshifts_distribution(timestamps_unconverted)
    #p = plot_dataset_range(data , station=station , text = text, sources=sources)
    p = plot_dataset_distribution(data, date_min='', date_max='',  station=station , text = text, sources=sources )
    
    




def make_obs_tab(file, var="", resorted=True):
    f = h5.File( merged_file, 'r' )
    ot = f['observations_table']
    dic = {}

    if var:
        if resorted:
            ind_max, ind_min = max(f['recordindices'][var]) , min(f['recordindices'][var])
    
            for o in ['date_time' , 'observation_value' , 'observed_variable', 'z_coordinate', 'z_coordinate_type', 'source_id']:
                if o != 'source_id':
                    dic[o] = ot[o][ind_min:ind_max]
                else:
                    sid = [ b''.join(s).decode('utf-8') for s in ot[o][ind_min:ind_max] ] 
                    dic[o] = sid
        else:
            ind = np.where (  ot['observed_variable'] == var ) 
            
            for o in ['date_time' , 'observation_value' , 'observed_variable', 'z_coordinate', 'z_coordinate_type', 'source_id']:
                if o != 'source_id':
                    dic[o] = ot[o][ind_min:ind_max]
                else:
                    sid = [ b''.join(s).decode('utf-8') for s in ot[o][ind_min:ind_max] ] 
                    dic[o] = sid        
    else:
        
        for o in ['date_time' , 'observation_value' , 'observed_variable', 'z_coordinate', 'z_coordinate_type', 'source_id']:
            if o != 'source_id':
                dic[o] = ot[o][:]
            else:
                sid = [ b''.join(s).decode('utf-8') for s in ot[o][:] ] 
                dic[o] = sid    
        
            
    obs_df = pd.DataFrame( dic )
    
    return obs_df 
    
    
    
def plot_dataset_hist(data= "" ,  variable='' , p_levels= '' , station= '' , date_range= ['',''] , text='' , text_save = '', sources=''):
    """ Function to plot the distributions of datasets used in the merged file.
        Works with header_table or observations_table df (with header is faster) """

    """ Retrieving the style dictionaries """
    datasets_dic = Common.datasets_dic
    style_dic    = Common.style_dic
    fs = Common.fontsize - 2
    
    def count_data(data,sources):
        """ Counts data for the plot, make lables etc. """
        occ = dict(Counter(data["source_id"]) ) # Counter returns a dic with keys= all the items in the list, values=number of occurences per each item 
        for s in sources: #  putting back the sources with zero entries for nicer plots 
            if s not in list(occ.keys()):
                occ[s] = 0
       
        counts , labels, color = [], [] , []
        
        for source in Common.datasets:  # the double loops allows to plot the datasets always in the order defined in the list Common.all_sources
            for k,v in occ.items():
                if source == k:
                    labels.append(datasets_dic[k]["l"])
                    counts.append(v)
                    color .append( datasets_dic[k]["c"] )
                else:
                    continue 

        x = np.arange(len(occ)) # np.arange(3) = array([0, 1, 2])

        return counts, labels, color, x

   
    def plot_bar(data= '', ax = '', rotation = False , text = '' , station = '', sources=''):
        """ Creates a simple bar plots """
        counts, labels, color, x = count_data(data, sources)
        
        if not ax:
            fig, ax = plt.subplots(figsize=(10,7))
            
        formatter = ticker.ScalarFormatter(useMathText=True)
        formatter.set_scientific(True) 
        formatter.set_powerlimits((-1,1)) 
        ax.yaxis.set_major_formatter(formatter) 
    
        ax.grid(ls =":" , color = "lightgray")
        ax.bar(x, counts, color = color )
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize = fs )
        if rotation:
            ax.set_xticklabels(labels, rotation = 45, fontsize = fs )

    
    """ Global counts (all pressure levels) """
    #print("Plotting ::: Global counts (all pressure levels)")
    a = plot_bar(data=data , text= text, station = station, sources='')
    
    plt.ylabel("Data Counts (all records)", fontsize = fs )
    plt.title('Record counts per data source for Station ' + str(station) + ' ' + text , y=1.02 , fontsize = fs )

    plt.savefig( out_dir + '/' + station + "_datasets_distribution_global.pdf",
                bbox_inches = 'tight' , dpi = 250 )    
    plt.show()
    plt.close()    
        
    """ Global standard level counts (only available for observations_table ) """
    #std_plevs = [10, 20, 30, 50, 70, 100, 150, 200, 250, 300, 400, 500, 700, 850, 925, 1000]    
    std_plevs = Common.std_plevs
    
    standard_lev_data  = data.loc[ (data['z_coordinate'] == 10   ) |  
                                   (data['z_coordinate'] == 20   ) |
                                   (data['z_coordinate'] == 30   ) |
                                   (data['z_coordinate'] == 50   ) |
                                   (data['z_coordinate'] == 70   ) |
                                   (data['z_coordinate'] == 100  ) |
                                   (data['z_coordinate'] == 150  ) |
                                   (data['z_coordinate'] == 200  ) |
                                   (data['z_coordinate'] == 250  ) |
                                   (data['z_coordinate'] == 300  ) |
                                   (data['z_coordinate'] == 400  ) |
                                   (data['z_coordinate'] == 500  ) |
                                   (data['z_coordinate'] == 700  ) |
                                   (data['z_coordinate'] == 850  ) |
                                   (data['z_coordinate'] == 925  ) |
                                   (data['z_coordinate'] == 1000 )   ] 
    
    print("Printing ::: Global standard level counts (only available for observations_table ) ")
    a = plot_bar(data= standard_lev_data , station = station, text = text )
    plt.title('Data counts per data source for Station ' + str(station) + ' ' + text , y=1.02 , fontsize = fs )

    plt.ylabel("Counts (standard pressure levels)", fontsize = fs )
    plt.savefig(out_dir + '/'  + station + "_data_datasets_distribution_standard" + text_save + ".pdf", bbox_inches = 'tight' , dpi = 250 )      
    #plt.show()

    plt.close()
    
    """ Per-variable dataset distributions (all pressure levels) """
    print("Plotting ::: Per-variable dataset distributions (all pressure levels)")
    fig, axs = plt.subplots(nrows=2, ncols=3, constrained_layout=True , figsize=(18,10) )
    fig.suptitle('Data counts per variable and data source for Station ' + str(station) + ' ' + text , y=1.04 , fontsize = fs+3)

    for ax,v in zip(axs.flat, variables):
        v = int(v)
        ax.set_ylabel("Counts (all records) - " + style_dic[v]['l'] , fontsize = fs )
        data_v = data.loc[ (data['observed_variable'] == v) ]
        plot_bar(data= data_v , ax = ax , rotation = True , station = station, text = text )

    plt.savefig(out_dir + '/'  + station + "_data_datasets_distribution_perVariable_" + text_save + ".pdf", bbox_inches = 'tight' , dpi = 250 )   
    plt.show()
    plt.close()

    
    return 0


def plot_timeseries(resorted= True,
                    file='',
                    variables= '',
                    p_levels= '', 
                    station= '', 
                    date_range= ['',''], 
                    print_df = False, 
                    text = '',
                    sources=''):
        
    """ Function to plot time series """

    datasets_dic = Common.datasets_dic
    style_dic    = Common.style_dic
    fs           = Common.fontsize 

    date_min, date_max = date_range[0], date_range[1]
    
    for v in variables:
        v = str(v)
        data_v = make_obs_tab(file, var=v, resorted=resorted)

        for p in p_levels:
            
            data_p = data_v.loc[ (data_v["z_coordinate"] == p ) & ( data_v["z_coordinate_type"] == 1)  ] # filtering on various pressure levels

            #plt.figure(figsize=(20,10)) 
            fig, ax = plt.subplots(figsize=(16,9))
            ax.tick_params(axis='both', which='major', labelsize=15)
            ax.tick_params(axis='both', which='minor', labelsize=8)
        
            plt.grid(ls =":" , color = "lightgray")
            
            if  date_min and date_max:
                data_p = data_p.loc[ (data_p["date_time"] >= date_min ) & (data_p["date_time"] < date_max )  ] # option filtering on date range

                plt.xlim(date_min, date_max)
            
            source_datasets = np.unique(data_p["source_id"][:])
            print("All sources::: " , source_datasets )
    
            for s in sources:
            #for s in source_datasets:
                d = data_p.loc[ (data_p['source_id'] == s) ]
                x = d["date_time"].values
                y = d["observation_value"].values
                #S = s
                #print("Length for " , str(S) , " is: ", len(x) ) 
            
                num    = '[' + str(len(x)) + ']'
                legend =  datasets_dic[s]['l'] + ' ' + num 
                color  = datasets_dic[s]['c']
            
                #plt.plot(x, y, label = legend , color = color )
                plt.scatter(x, y, label = legend , color = color )
                
                """ Print here the dataframe """
                if print_df:
                    dates = [datetime.datetime(1959,12,14) , datetime.datetime(1959,12,15)]
                    df = d.loc [ (d["date_time"] >= dates[0]) & (d["date_time"] < dates[1])     ] 
                    pd.set_option('expand_frame_repr', False)

                    if not df.empty:
                        print( '  ' + s + '--------------------------------------------------------------------------------' , '\n', df)
        
            #plt.xlabel("Date Time", fontsize = fs)
            pressure = str(int(p/100) )
            v = int(v)
            plt.ylabel( style_dic[v]['l'] , fontsize = fs + 2 )
            
            #levels    = [100000, 50000, 1000]    

            #plt.text(0.73, 0.93, 'Plevel=' + pressure + ' [hPa]' , transform=ax.transAxes, color = 'red', fontsize = fs+5)
            
            plt.legend(loc = 'best', fontsize = fs)    

            
            plt.title("Station " + station + ' - Plevel=' + pressure + ' [hPa]' , fontsize = fs + 4, y = 1.03 )
            
            plt.savefig(out_dir + '/' + station + "_time_series_" + "_pressure_" + pressure + "_var_" + str(v) + text + ".pdf",
                         bbox_inches = 'tight' , dpi = 250 )
            plt.show()
            plt.close()
            print("*** Plot created for variable ", v , "  at pressure level " , str(p) )
        
        




data = {'path': '/mnt/users/scratch/leo/scratch/converted_v8/0-20000-0-58362_CEUAS_merged_v1.nc',
              'name': 'Chinese Stat',
              'station': '0-20000-0-58362'}
        
data = {'path': '/mnt/users/scratch/leo/scratch/converted_v8/0-20000-0-82930_CEUAS_merged_v1.nc',
              'name': 'Chinese Stat',
              'station': '0-20000-0-82930'}

def get_basic_data(dic):
    merged_file = data['path']
    station_name = data['name']
    station = data['station']
    


    
    
    # timestamps from header
    header = h5.File( merged_file, 'r' )['header_table']
    #report_timestamp = pd.to_datetime( header['report_timestamp'][:],  unit='s',  origin=pd.Timestamp('1900-01-01') ) 
    record_timestamp =  pd.to_datetime( header['record_timestamp'][:],  unit='s',  origin=pd.Timestamp('1900-01-01') ) 
    source_id = [ b''.join(s).decode('utf-8') for s in header['source_id'][:] ]


    header_df = pd.DataFrame ( { 'record_timestamp': record_timestamp ,
                               'record_timestamp_unconverted': header['record_timestamp'][:],
                               'source_id': source_id } )

    sources = list(np.unique( source_id ) )

    return merged_file, station_name, station, header, header_df, sources 


############################### extracting some basic data
merged_file, station_name, station, header, header_df, sources = get_basic_data(data)

###############################  headers 
#analyze_headers(data= header_df, station = station, date_min='', date_max='', text = '' , text_save='' , sources=sources )
        



############################### time series 

variables = ['126',"139","140","106","107","117"]
for v in variables:
    
    #a = make_obs_tab(merged_file, var=v)

    """
    plot_timeseries(variables= [ v ],
                         p_levels= [50000, 10000], 
                         file= merged_file,
                         station= station, 
                        date_range= ['',''], 
                        print_df = False, 
                        text = '',
                        sources=sources)

   """
    
############################### global histos 
a = make_obs_tab(merged_file, var=False)    
a = plot_dataset_hist(data= a,  variable=v , p_levels= '' , station= '' , date_range= ['',''] , text='' , text_save = '', sources='')
