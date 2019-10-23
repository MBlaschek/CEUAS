""" Merging the station configuration files """

import os,sys
import netCDF4 as nc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import argparse
from datetime import datetime, timedelta
import numpy.ma as ma

def to_integer(dt_time):
      return 10000*dt_time.year + 100*dt_time.month + dt_time.day

class netCDF():
      """ Class to handle the converted netCDF files """
      def __init__(self, file=''):
            self.f = file 
            self.df = nc.Dataset(file)
            self.groups = self.df.groups.keys()
            
      def stationConfiguration(): 
            stat_conf = elf.df.groups['station_configuration']
            self.station_configuration = stat_conf     
            self.primary_id        = stat_conf['primary_id']
            self.secondary_id    = stat_conf['secondary_id']
            self.lat                   = stat_conf['latitude']
            self.lon                  = stat_conf['longitude']
            self.station_name  = stat_conf['station_name']
            self.start_date       = stat_conf['start_date']
            self.end_date        = stat_conf['end_date']        
            
      
class stationConfigurationTable():
      """ Reading the CDM station_configuration tables in cvs format """

      def __init__(self, file=''):
          self.f = file 
          self.df = ''
      def readDF(self):
            """ Read the station_configuration as a Panda DF, and stor some variables """   
            
            try:
                  df = pd.read_csv(self.f, delimiter='\t', quoting=3, na_filter=False, comment='#')
                  self.df =  df
                  #self.primary_id        = df['primary_id']
                  #self.secondary_id    = df['secondary_id']
                  #self.lat                   = df['latitude']
                  #self.lon                  = df['longitude']
                  #self.station_name  = df['station_name']
                  #self.start_date       = df['start_date']
                  #self.end_date        = df['end_date']                  
                  
                  self.columns          = df.columns
            except IOError:
                print(' The station configuration file ', self.f , ' cannot be found!' ) 
                self.df = 0
              
      def extract_column(self, column):
            """ Extract the value of the variable in input from the dataframe """
            try:
                  return  self.df[column]
            except KeyError:
                  print('This columns does not exist in the pandas dataframe of the station_configuration file!')
                  print('The available columns are: ', df.columns )
      
      
class StationConfigurationAnalizer():
            """ Build a dictionary from the station_configuration"""
            
            def __init__(self, bufr = '', igra2 = '', ncar = '', era5_1 = '' , era5_1759 = '' , era5_1761 = '' , era5_3188 = ''):
                  """ Initilaize the station_configuration files, e.g. sc = stationConfiguration(file= stat_conf_file)  
                        and the files containing the primary.
                  , secondary and original file names """
                  
                  self.bufr = bufr
                  self.igra2 = igra2
                  self.ncar  = ncar
                  self.era5_1 = era5_1
                  self.era5_3188 = era5_3188
                  self.era5_1759 = era5_1759
                  self.era5_1761 = era5_1761 
                  
                  self.dic = '' 
                  self.all_primary_id = ''
                  self.igra2_file_finder = ''
                  self.bufr_file_finder           = pd.read_csv( 'stations_configurations/bufr_stationIds_filename.dat', sep='\t',  names = ['primary_id', 'secondary_id', 'file'] ) 
                  #self.igra2_file_finder = pd.read_csv('igra2_stationIds_filename.dat', 'r')
                  self.ncar_file_finder           = pd.read_csv('stations_configurations/ncar_stationIds_filename.dat' , sep='\t', names = ['primary_id', 'secondary_id', 'file'] )   
                  self.era5_3188_file_finder = pd.read_csv('stations_configurations/3188_stationIds_filename.dat' , sep='\t', names = ['primary_id', 'secondary_id', 'file'] )   
                  self.era5_1759_file_finder = pd.read_csv('stations_configurations/1759_stationIds_filename.dat' , sep='\t', names = ['primary_id', 'secondary_id', 'file'] )   
                  self.era5_1761_file_finder = pd.read_csv('stations_configurations/1761_stationIds_filename.dat' , sep='\t', names = ['primary_id', 'secondary_id', 'file'] )   
                  #self.era5_1_file_finder       = pd.read_csv('1_stationIds_filename.dat'       , sep='\t', names = ['primary_id', 'secondary_id', 'file'])   

                           
            def doDic(self):
                  """ Create a compact dictionary with the main data """
                  data = {}
                  datasets = { 'datasets' : [] , 'df': [] }
                  
                  for k,v in {'bufr':self.bufr, 'igra2':self.igra2, 'ncar': self.ncar, 'era5_1':self.era5_1, 'era5_1759':self.era5_1759, 'era5_1761': self.era5_1761, 'era5_3188':self.era5_3188   }.items():
                        datasets['datasets'].append(k)
                        datasets['df'].append(v)
                  
                  for dataset,df in zip( datasets['datasets'], datasets['df'] ):
                        data[dataset] = {}
                        data[dataset]['primary_id']       = list( df.df['primary_id'] )
                        data[dataset]['secondary_id']   = list( df.df['secondary_id'] )
                        data[dataset]['latitude']            = list( df.df['latitude'] )
                        data[dataset]['longitude']         = list( df.df['longitude'] )
                        data[dataset]['station_name']  = list(df.df['station_name'])
                        data[dataset]['start_date']        = list(df.df['start_date'])
                        data[dataset]['end_date']         = list(df.df['end_date'])                  
                  self.dic = data 
                  
            def extractUniqueStations(self):
                  """ Extract a list of unique stations primary_id from the datasets provided """
                  all_primary = list(self.igra2.df['primary_id']) 
                  all_primary = all_primary + list(self.bufr.df['primary_id']) 
                  all_primary = all_primary + list(self.ncar.df['primary_id']) 
                  
                  all_primary = all_primary + list(self.era5_1.df['primary_id']) 
                  all_primary = all_primary + list(self.era5_1759.df['primary_id']) 
                  all_primary = all_primary + list(self.era5_1761.df['primary_id']) 
                  all_primary = all_primary + list(self.era5_3188.df['primary_id']) 
                  
                  
                             
                  all_primary =  list( set(all_primary) )
                  
                  all_primary = [ p for p in  list( set(all_primary) ) if  isinstance(p, str)] # remove nans entries
                  all_primary = list( set(all_primary) )
                  self.all_primary_id = [ p for p in all_primary if p != '-1']
                  
            
            def extractDataFromDataset(self, primary_station= '', dataset='' ):
                  """ Extract the data by index in the corresponding list in the doDic() dictionary """
                  flag, station_name, lat, lon, start, end, original_file = '-1','-','-','-','-','-','-'
                  dic = self.dic[dataset]['primary_id']
                  if primary_station in dic:
                        flag = '1'
                        index = self.dic[dataset]['primary_id'].index(primary_station) 
                         
                        station_name     = str(self.dic[dataset]['station_name'][index])
                        lat, lon               = str(self.dic[dataset]['latitude'][index]), str(self.dic[dataset]['longitude'][index])
                        start, end           = str(self.dic[dataset]['start_date'][index]), str(self.dic[dataset]['end_date'][index] )          
                      
                        original_file = self.findOriginalFile(primary_station = primary_station, secondary_station = '', dataset = dataset)                              
                         
                  return flag, station_name, lat, lon, start, end , original_file
                      
            def f(self,stringa ):
                  """ Utility for nice printout in the summary file """
                  return "{:25}".format(stringa)
            
            def findOriginalFile(self, primary_station = '', secondary_station = '', dataset = ''):
                  """ Extract the original file(s) that were mapped to the primary station id 
                        The originale file names are extracted from the attribute file_finder
                        except for the igra2 datasets, for which the names are built directly
                        from the stations ids.
                  """ 
                  if dataset in ['era5_1']:  ######## TODO !!!! 
                         return '-'      
                  
                  found_file = ''
                  if dataset == 'bufr':
                        df =  self.bufr_file_finder
                  elif dataset == 'ncar':
                        df =  self.ncar_file_finder
                                             
                  elif dataset == 'era5_1759':
                        df = self.era5_1759_file_finder
                  elif dataset == 'era5_1761':
                        df = self.era5_1761_file_finder
                  elif dataset == 'era5_3188':
                        df = self.era5_3188_file_finder
                                          
                  if dataset == 'igra2':
                        igra2_path = '/raid60/scratch/federico/databases/IGRAv2/'
                        files = os.listdir(igra2_path)
                        found = '-1'
                        primary_station = primary_station.split('-')[3]
                        for f in files:
                              if primary_station in f:
                                    #print('found in igra2 file!')
                                    found_file = igra2_path + '/' +  f
                  
                  try:
                        find =  df.loc[ df['primary_id'] == primary_station ]  
                        found_file = find['file'].values[0]
                        return found_file
                  except:
                        return found_file 
                        
                       
            def basicSummary(self):
                  """ Write a txt summary out of the doDic() dictionary """
                  summary = open('summary_station_configuration.dat', 'w')
                  header = '#primary_id\tncar_flag\tigra2_flag\tbufr_flag\n'
                  summary.write(header)
                  
                  summary_forplot = open('summary_forplot.dat','w')
                  for i in self.all_primary_id:  # initialize some empty variable 
                        bufr_flag, ncar_flag, igra2_flag = '-1', '-1', '-1' 
                        ncar_stat, igra2_stat, bufr_stat = '','',''                        
                        ncar_lat, igra2_lat, bufr_lat, era5_1_lat , era5_1759_lat , era5_1761_lat , era5_3188_lat       = '-','-','-' , '-', '-', '-', '-'
                        ncar_lon, igra2_lon, bufr_lon, era5_1_lon , era5_1759_lon , era5_1761_lon , era5_3188_lon  = '-', '-', '-', '-', '-', '-', '-'
                        
                        
                        ncar_start, igra2_start, bufr_start, era5_1_start, era5_1759_start, era5_1761_start, era5_3188_start =  '-','-','-','-','-','-' ,'-'
                        bufr_end, ncar_end, igra2_end, era5_1_end, era5_1759_end, era5_1761_end, era5_3188_end = '-','-','-','-','-','-' ,'-'
                        
                        file_ncar , file_bufr , file_igra2, file_era5_1, file_era5_1759, file_era5_1761, file_era5_3188 = '-1', '-1' , '-1' , '-1',  '-1', '-1' , '-1' ,
                        
                        
                        bufr_flag , bufr_stat,  bufr_lat, bufr_lon , bufr_start, bufr_end, bufr_file             = self.extractDataFromDataset(dataset = 'bufr', primary_station = i)
                        ncar_flag , ncar_stat,  ncar_lat, ncar_lon , ncar_start, ncar_end, ncar_file         = self.extractDataFromDataset(dataset = 'ncar', primary_station = i)
                        igra2_flag , igra2_stat,  igra2_lat, igra2_lon , igra2_start, igra2_end, igra2_file = self.extractDataFromDataset(dataset = 'igra2', primary_station = i)
                        

                        era5_1_flag , era5_1_stat,  era5_1_lat, era5_1_lon , era5_1_start, era5_1_end, era5_1_file = self.extractDataFromDataset(dataset = 'era5_1', primary_station = i)
                        era5_1759_flag , era5_1759_stat,  era5_1759_lat, era5_1759_lon , era5_1759_start, era5_1759_end, era5_1759_file = self.extractDataFromDataset(dataset = 'era5_1759', primary_station = i)
                        era5_1761_flag , era5_1761_stat,  era5_1761_lat, era5_1761_lon , era5_1761_start, era5_1761_end, era5_1761_file = self.extractDataFromDataset(dataset = 'era5_1761', primary_station = i)
                        era5_3188_flag , era5_3188_stat,  era5_3188_lat, era5_3188_lon , era5_3188_start, era5_3188_end, era5_3188_file = self.extractDataFromDataset(dataset = 'era5_3188', primary_station = i)
                        

                        if ncar_flag == '1':
                              file_ncar = self.findOriginalFile(primary_station=i, dataset='ncar')
                        if bufr_flag == '1':
                              file_bufr = self.findOriginalFile(primary_station=i, dataset='bufr')                                 
                        if igra2_flag == '1':                              
                              file_igra2 = self.findOriginalFile(primary_station=i, dataset='igra2') 
                              
                        if era5_1759_flag == '1':                              
                              file_era5_1759 = self.findOriginalFile(primary_station=i, dataset='era5_1759') 
                        if era5_1761_flag == '1':                              
                              file_era5_1761 = self.findOriginalFile(primary_station=i, dataset='era5_1761')
                        if era5_3188_flag == '1':                              
                              file_era5_3188 = self.findOriginalFile(primary_station=i, dataset='era5_3188') 
                              
                        interesting = open('interesting_files.dat','a+')
                        
                        if ncar_flag == '1' and bufr_flag == '1' and igra2_flag == '1'  and era5_1759_flag == '1' and era5_1761_flag == '1' and era5_3188_flag == '1':
                              print ('fund' , str (self.all_primary_id.index(i)) )
                              interesting.write(i + ',' + file_ncar + ',' +  file_bufr + ',' + file_igra2 + ',' + file_era5_1759 + ',' + file_era5_1761 + ',' + file_era5_3188 + '\n')
                        
                        
                        l_flags = self.f(i) + '\t' +  self.f(ncar_flag) + '\t' + self.f(igra2_flag) + '\t' + self.f(bufr_flag) + '\t' +  self.f(bufr_flag) 
                        l_files = file_ncar + '\t' + file_igra2 + '\t' + file_bufr + '\n'  
                              
                        summary.write(l_flags)    
                        summary.write(l_files)
                        
                        ncar_l  = '# ncar\t\t'  + ncar_stat  + '\t\t' + ncar_lat   + '\t\t' + ncar_lon + '\t\t' + ncar_start   + '\t\t' + ncar_end  + '\n' 
                        bufr_l  = '# bufr\t\t'   + bufr_stat   + '\t\t' + bufr_lat   + '\t\t' + bufr_lon  + '\t\t' + bufr_start   + '\t\t' + bufr_end   + '\n' 
                        igra2_l = '# igra2\t\t' + igra2_stat + '\t\t' + igra2_lat + '\t\t' + igra2_lon + '\t\t' + igra2_start + '\t\t' + igra2_end + '\n' 
                              
                        summary.write(ncar_l)                                     
                        summary.write(bufr_l)       
                        summary.write(igra2_l)                                                          
                              
                        l_forplot =                   ncar_lat + ',' + ncar_lon  + ',' + ncar_start   + ','  + ncar_end  + ','
                        l_forplot = l_forplot + igra2_lat + ',' + igra2_lon + ',' + igra2_start + ','  + igra2_end + ','
                        l_forplot = l_forplot + bufr_lat   + ',' + bufr_lon   + ',' + bufr_start   + ',' + bufr_end   + ','

                        l_forplot = l_forplot + era5_1_lat       + ',' + era5_1_lon       + ',' + era5_1_start       + ',' + era5_1_end + ','
                        l_forplot = l_forplot + era5_1759_lat + ',' + era5_1759_lon + ',' + era5_1759_start + ',' + era5_1759_end + ','
                        l_forplot = l_forplot + era5_1761_lat + ',' + era5_1761_lon + ',' + era5_1761_start + ',' + era5_1761_end + ','
                        l_forplot = l_forplot + era5_3188_lat + ',' + era5_3188_lon + ',' + era5_3188_start + ',' + era5_3188_end + '\n'
                             
                        summary_forplot.write(l_forplot)
                  summary.close()
                              
                        #print(i , bufr_flag, ncar_flag, igra2_flag )

                        
class Merger():
      """ Class for the merging of the data from different netCDF files """
      def __init__(self ):
            self.data = {}
            self.databases = []
            #self.df = {}            # storing the dafarames 
            #self.DataDic = {}   # storing the info from the dataframe 
            #self.databases = ['ncar','igra2','bufr','era5_1','era5_1759','era5_1761','era5_3188']
            
      def InitializeData(self, datasets = {} ):
            """ Initialize dataset. 
                       Args:: datasets (dictionary where key+ dataset name e.g. bufr, igra2 etc. , and value is the path to the netCDF file """
            data = {}
            for k,v in datasets.items() :
                  data[k] = {}
                  ds =  nc.Dataset(v) 
                  data[k]['df'] = ds                
                  data[k]['dateindex'] = ds.variables['dateindex'][0,:]  # storing the dateindex 
                  data[k]['observations_table'] = {}
                  data[k]['observations_table']['date_time'] = ds.groups['observations_table']['date_time']   # storing the date_time 
                  
                  data[k]['source_file']      = ds.groups['source_configuration']['source_file'][0]
                  data[k]['product_code']  = ds.groups['source_configuration']['product_code'][0]
                  
                  
            self.databases = datasets.keys()
            self.data = data 
            
      def plot_styler(self):
            """ Building a style dictionary to make all the plots uniform """
            
            style_dic = {  'bufr'           : {'color':'cyan'           , 'size':8    , 'label':'BUFR'           , } ,
                                  'igra2'          : {'color':'orange'       , 'size' :6    , 'label':'IGRA2'         , } ,
                                  'ncar'           : {'color': 'magenta'  , 'size' :10  , 'label':'NCAR'           , } ,                      
                                  'era5_1'       : {'color': 'yellow'      , 'size' : 4  , 'label':'ERA5 1'          , } ,
                                  'era5_1759' : {'color': 'slateblue'  , 'size' : 3  , 'label':'ERA5 1759'    , } ,
                                  'era5_1761' : {'color': 'lime'          , 'size' : 2   , 'label':'ERA5 1761'   , } ,
                                  'era5_3188' : {'color': 'blue'          , 'size' : 1   , 'label':'ERA5 3188'   , } ,    }
            
            
            station , name        = self.data['ncar']['df'].groups['station_configuration']['primary_id'][0]      ,   self.data['ncar']['df'].groups['station_configuration']['station_name'][0]
            latitude, longitude = str (self.data['ncar']['df'].groups['station_configuration']['latitude'][0] ) ,    str (self.data['ncar']['df'].groups['station_configuration']['longitude'][0] )            
            
            style_dic_variable = {85   : {'label':'Temperature [C]'              , 'ylim' : [-40 , 40    ]  } ,
                                               107 : {'label': 'Wind Speed [m/s]'           , 'ylim' : [0     ,80    ]  }  ,
                                               106 : {'label': 'Wind Direction [degree]' , 'ylim' : [0     , 360 ]  }   ,
                                               38   : {'label': 'Relative Humidty'           , 'ylim'  : [0     , 1     ]  }   ,
                                               117 : {'label': 'Geopotential'                  , 'ylim'  : [0    , 100  ]  }   ,
                                               
                                               
                                               
                                               }  
            
            style_dic['station']     = station
            style_dic['name']       = name
            style_dic['latitude']    = latitude
            style_dic['longitude'] = longitude
            
            self.style_dic = style_dic
            self.style_dic_variable = style_dic_variable
            
            
      def PlotTimeDistribution(self):
            """ Script to plot the time distribution of the available data, by reading the date"""           
            # station info (statio ind, name, lat, lon)


            fig, ax = plt.subplots()            
            fig.set_size_inches(10, 4)
            plt.tick_params(
                axis='y',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                bottom=False,      # ticks along the bottom edge are off
                top=False,         # ticks along the top edge are off
                labelbottom=False)
            
            for d, position in zip ( self.databases, [1, 1.5,  2, 2.5, 3, 3.5, 4] ):
                  dates = self.data[d]['dateindex']
                  plt.scatter (dates, np.full( (len(dates)) , position ) , label = self.style_dic[d]['label'] , color =    self.style_dic[d]['color'] )

            plt.ylim(0,5)
            size = 13
            plt.text(19200000,4.5, 'Station Identifier: ' + self.style_dic['station'], fontsize = size)
            plt.text(19200000,4.1, 'Station Name: ' + self.style_dic['name'], fontsize = size)
            plt.text(19200000,3.7, 'Lat, Lon:  [' + self.style_dic['latitude'] + ','  + self.style_dic['longitude'] + ']', fontsize = size)
            

            #plt.yticks( (1, 1.5, 2 , 2.5 , 3, 3.5 , 4 ) ,('Tom', 'Dick', 'Harry', 'Sally', 'Sue', 'a', 'b') )
            
            plt.grid(linestyle = '--', color  ='lightgray', axis = 'x')
            plt.legend(fontsize = 10, loc = 'lower left', ncol = 2)
            #plt.grid()
            plt.xlabel('Observation Dates', fontsize =12)
            ax.ticklabel_format(useOffset=False, style='plain' )    
            ax.set_yticklabels(('', '', '','','','')) 
            plt.xticks(rotation = 45, fontsize = 8)
            plt.savefig('plot_directory/plots/time_series/' + self.style_dic['station'] + '_timeintervals.png', bbox_inches='tight', dpi = 200)

      def MakeDateTime(self, dataset='', shortener = 10000):   # only consider a small number of entries           
            """ Extracting the actual date_time from the time offset and the time deltas stored in ['observations_table']['date_time'] """
            date_times = self.data[dataset]['observations_table']['date_time']
            time_offset_units = date_times.units  
            time_offset           = time_offset_units.split('since ')[1].split(' ')[0]                             
            time_offset           = datetime.strptime(time_offset, '%Y-%m-%d')
            if shortener: 
                  date_time = set(date_times[:shortener])
            else:
                  date_time = set(date_times[:])
            
                  
            if 'minutes' in time_offset_units:
                  delta = [ timedelta(minutes = float(i) ) for i in date_time ]
            elif 'hours' in time_offset_units:
                  delta = [ timedelta(hours = float(i) )    for i in date_time ]
                  
            #date_time = [to_integer(i) for i in  [  time_offset + i  for i in delta  ] ]
            dt = [i for i in  [  time_offset + i  for i in delta  ] ]                              

            return dt      
            
            
            
            
      def Extract_variable_pressure_time(self, variables='', pressures = '' ):
            """ Given a dataset, extract the values of a variable for a certain pressure from the observation tables. 
                  Stores the data in a dictionary, separated """
            data = {}
            #for d in ['era5_3188']:
            
            for d in self.databases:
              
                  data[d] = {}
                  
                  obs_variable = self.data[d]['df'].groups['observations_table']['observed_variable'][:]
                  obs_values    = self.data[d]['df'].groups['observations_table']['observation_value'][:]
                  obs_id           = self.data[d]['df'].groups['observations_table']['observation_id'][:]
                  
                  p_levels         = self.data[d]['df'].groups['observations_table']['z_coordinate'][:]    
                  date_times    = self.data[d]['df'].groups['observations_table']['date_time']                                   
                  

                  data[d]['date_time']      = date_times[:] 
                   
                  # date_times are stored as a time interval from a certain starting point, which is given as an attribute of the date_time variable in the observation_table.
                  # This can be in minutes, hours, days etc. depending on the input data.
                  # We need to extract the time offset and add it to the entries store in the date_time obs table, to calculate the correct real date_time value 
                  time_offset_units = date_times.units  
                  time_offset           = time_offset_units.split('since ')[1].split(' ')[0]                             
                  time_offset           = datetime.strptime(time_offset, '%Y-%m-%d')
                  
                  print ('For the database' , d , ' the units of time are: ',  self.data[d]['df']['observations_table']['date_time'].units ) 
              
                  for v in variables:
                        data[d][v]={}
                        for p in pressures:
                              print('Processing: ' , d , ' ' , v , ' ' , p )
                              data[d][v][p] = {'observed_values': [] , 'indices': [] , 'date_time': [] , 'observation_id':[]  }
                              
                              pressure_indices  = np.where( p_levels == p)[0]   # extracting pressure levels matching with input pressure  
                              variables_indices = np.where( obs_variable == v)[0]  # extracting observed variable matching with input var                               
                              indices = list(set(pressure_indices).intersection(set(variables_indices))) # intersection
                                                           
                              data[d][v][p]['observed_values'] = np.take(obs_values, indices)
                              data[d][v][p]['observation_id']    = np.take(obs_id, indices)
                          
                              data[d][v][p]['indices']                = indices
                                                                                     
                              #pressure_indices    = np.where( p_levels == p)[0]  # extracting indices of the pressure levels where p == pressure 
                              #variables = np.take(obs_variable, pressure_indices)
                  
                              dt = np.take( date_times[:] , indices ) # extracting matchign date_time 
                              
                              if 'minutes' in time_offset_units:
                                    delta = [ timedelta(minutes = float(i) ) for i in dt ]
                              elif 'hour' in time_offset_units:
                                    delta = [ timedelta(hours = float(i) )    for i in dt ]
                                    
                              #date_time = [to_integer(i) for i in  [  time_offset + i  for i in delta  ] ]
                              dt = [i for i in  [  time_offset + i  for i in delta  ] ]                              
                              data[d][v][p]['date_time']  = dt

                  print('done with d: ', d)
                  
            return data
      
      
      def Merge(self, data = ''):
            """ Module to create a merged netCDF file from the input data. 
                  First step: creation of a list of available dates. 
                  If dates and time differ for less than the specified time delta, they are considered to represent the same observation """
            
            # Extract the dates and time. I  loop over all possible datetime of the various datasets,
            # which become a key of a dictionary. The value of the dic is a list that will contain the name of the datasets which have data for that specific date_time observation 
                       
            all_date_time = {} # list of all the dates and times 
            new_times = {}
            
            #for d in self.databases :            
            for d in ['era5_1761','igra2'] :
                  new_times[d] = []
                  
                  time_date_dataset = self.MakeDateTime(dataset = d , shortener = 1000)  # extracting the real date_time values 
                  #if d == 'igra2':
                  #      time_date_dataset = [datetime(1979, 1, 1, 0, 0) , datetime(1979, 1, 1, 0, 0) , datetime(1979, 1, 1, 0, 0) , datetime(1979, 1, 1, 0, 0) , datetime(1979, 1, 1, 0, 0) , datetime(1979, 1, 1, 0, 0)  ]
                  for td in time_date_dataset:
                        if td not in all_date_time.keys():
                              print('not found' , td )
                              try: 
                                    if d not in all_date_time[td]:
                                          all_date_time[td].append(d)
                              except KeyError:
                                    all_date_time[td] = []
                                    all_date_time[td].append(d)
                        else:
                              print(' found ' , td )
                              if d not in all_date_time[td]:                           
                                    all_date_time[td].append(d)                              
            print('hello')
            
      
      def PlotTimeSeries(self , data_dic = '', variable= '', pressure = ''):
            """ Make a time series plot for the variable, pressure levels specfified. Only identical dates are considered. """
            station = self.data['ncar']['df'].groups['station_configuration']['primary_id'][0] 
            self.plot_styler()
            fig, ax = plt.subplots()            
            fig.set_size_inches(15, 4)            
            #for d in self.databases:
            
            date_min, date_max = datetime.strptime('2100-01-01', '%Y-%m-%d') , datetime.strptime('1900-01-01', '%Y-%m-%d')
            min_obs, max_obs = 999999, 0.00001
            for d in self.databases: 
                
                  obs = data_dic[d][variable][pressure]['observed_values'] 
                  if d in ['era5_1', 'era5_3188','era5_1759','era5_1761', 'bufr'] and variable == 85:
                       obs = [ i - 273.15 for i in obs ]  # convert to Celsius 
                  print('database: ***** ', d , '   ', obs)
                  #time = np.datetime64(a[d][85][100000]['date_time'])
                  time = data_dic[d][variable][pressure]['date_time']
                  try:               
                        if min(time) < date_min :
                              date_min = min(time) 
                        if max(time) > date_max:
                              date_max = max(time)
                        if min(obs) < min_obs:
                              min_obs = min(obs)
                        if max(obs) > max_obs:
                              max_obs = max(obs)                              
                  except:
                        pass
                        
                  plt.scatter(time, obs , s = self.style_dic[d]['size'] , color = self.style_dic[d]['color'] , label = self.style_dic[d]['label'] )
            
            #fig.canvas.draw()
            #labels = [item.get_text() for item in ax.get_xticklabels()]                  
            #labels = [ i[:5] for i in labels ]                  
            #ax.set_xticklabels(labels)                  

            plt.title('Station  ' + self.style_dic['name'] + ' [ ' + self.style_dic['station'] + ' ]  for P = ' + str(pressure) + ' Pa', fontsize = 12 , y = 1.03)
            plt.grid(linestyle = '--', color = 'lightgray')

            plt.ylim(self.style_dic_variable[variable]['ylim'][0],    self.style_dic_variable[variable]['ylim'][1],  )
            plt.ylim(-50 , 50 )
            
            plt.ylim( min_obs - min_obs/10.  ,    max_obs - max_obs/10.   )
            
            plt.ylabel(self.style_dic_variable[variable]['label'])
            
            plt.legend(fontsize = 12, loc = 'lower left')
         
            plt.xlim( date_min , datetime(1975, 1, 1)       )                  
            plt.savefig('plot_directory/plots/time_series/' + str(variable) + '_'  + station + '_' + str(pressure) + '_timeseries_low.png', dpi = 200)
            plt.xlim( datetime(1975, 1, 1)   , date_max   )            
            plt.savefig('plot_directory/plots/time_series/' + str(variable) + '_' + station + '_' + str(pressure) + '_timeseries_high.png', dpi = 200)            
            fig.set_size_inches(20, 4)         
            plt.xlim( date_min, date_max  )
            plt.savefig('plot_directory/plots/time_series/' + str(variable) + '_' + station + '_' + str(pressure) + '_timeseries_all.png', dpi = 200)
            plt.close()

""" 
['primary_id', 'primary_id_scheme', 'record_number', 'secondary_id', 'secondary_id_scheme', 'station_name', 'station_abbreviation', 'alternative_name', 'station_crs', 'longitude', 'latitude',
'local_gravity', 'start_date', 'end_date', 'station_type', 'platform_type', 'platform_sub_type', 'operating_institute', 'operating_territory', 'city', 'contact', 'role', 'observing_frequency',
'reporting_time', 'telecommunication_method', 'station_automation', 'measuring_system_model', 'measuring_system_id', 'observed_variables', 'comment', 'optional_data', 'bbox_min_longitude',
'bbox_max_longitude', 'bbox_min_latitude', 'bbox_max_latitude', 'metadata_contact', 'metadata_contact_role']


['primary_id',  'secondary_id', 'station_name', 'alternative_name', 'station_crs', 'longitude', 'latitude',
'start_date', 'end_date',  'city', 
'observed_variables', ]
"""


make_summary = False

if make_summary:
      
      """ Read the station_configuration files """
      stat_conf_file = 'station_configuration_xxx.dat'
      bufr_sc   = stationConfigurationTable( file= 'stations_configurations/station_configuration_xxx.dat'.replace('xxx', 'bufr' )  )
      bufr_sc.readDF()
      igra2_sc = stationConfigurationTable( file= 'stations_configurations/station_configuration_xxx.dat'.replace('xxx', 'igra2' ) ) 
      igra2_sc.readDF()
      ncar_sc  = stationConfigurationTable( file= 'stations_configurations/station_configuration_xxx.dat'.replace('xxx', 'ncar' )  )
      ncar_sc.readDF()
      era5_1_sc  = stationConfigurationTable( file= 'stations_configurations/station_configuration_xxx.dat'.replace('xxx', 'era5_1' )  )
      era5_1_sc.readDF()
      era5_1759_sc  = stationConfigurationTable( file= 'stations_configurations/station_configuration_xxx.dat'.replace('xxx', 'era5_1759' )  )
      era5_1759_sc.readDF()
      era5_1761_sc  = stationConfigurationTable( file= 'stations_configurations/station_configuration_xxx.dat'.replace('xxx', 'era5_1761' )  )
      era5_1761_sc.readDF()
      era5_3188_sc  = stationConfigurationTable( file= 'stations_configurations/station_configuration_xxx.dat'.replace('xxx', 'era5_3188' )  )
      era5_3188_sc.readDF()
      
      

      """ Initialize the Analizer """
      analizer = StationConfigurationAnalizer (bufr = bufr_sc, 
                                                                         ncar = ncar_sc, 
                                                                         igra2 = igra2_sc, 
                                                                         era5_1 = era5_1_sc, 
                                                                         era5_1759 = era5_1759_sc ,  
                                                                         era5_1761 = era5_1761_sc, 
                                                                         era5_3188 = era5_3188_sc )
       

      analizer.extractUniqueStations()
      analizer.doDic()                                      # creates a compact dictionary from the station_configuration files 
      all_primary_id = analizer.all_primary_id # list of all the unique stations primary_ids, appearing at least in one of the datasets considered

      analizer.basicSummary()





""" Starting with the merging. Information stored in the summay file required (i.e. primary station id and filenames). """

data = { 'ncar'    : 'example_stations/ncar/chuadb_windc_22802.txt.nc'   ,
               'igra2'   : 'example_stations/igra2/chRSM00022802-data.txt.nc'  ,
               'bufr'     : 'example_stations/bufr/chera5.22802.bfr.nc'  ,
               
               'era5_1' : 'example_stations/era5_1/chera5.conv._22802.nc' , 
               
               'era5_1759' : 'example_stations/era5_1759/chera5.1759.conv.1:22802.nc' , 
               'era5_1761' : 'example_stations/era5_1761/chera5.1761.conv.1:22802.nc' , 
               'era5_3188' : 'example_stations/era5_3188/chera5.3188.conv.C:4648.nc' , 

               }



available_data = { 'ncar'    : 'example_stations/ncar/chuadb_windc_47646.txt.nc'   ,
               'igra2'   : 'example_stations/igra2/chJAM00047646-data.txt.nc'  ,
               'bufr'     : 'example_stations/bufr/chera5.47646.bfr.nc'  ,
               
               'era5_1' : 'example_stations/era5_1/chera5.conv._47646.nc' , 
               'era5_1759' : 'example_stations/era5_1759/chera5.1759.conv.1:47646.nc' , 
               'era5_1761' : 'example_stations/era5_1761/chera5.1761.conv.1:47646.nc' , 
               'era5_3188' : 'example_stations/era5_3188/chera5.3188.conv.C:5357.nc' , 

               }



#ncar  = nc.Dataset(ncar_test)
#igra2 = nc.Dataset(igra2_test)
#bufr   = nc.Dataset(bufr_test)

#ncar_di   = ncar.variables['dateindex'][0,:] # dateindex variable
#gra2_di  = igra2.variables['dateindex'][0,:] 
#bufr_di    = bufr.variables['dateindex'][0,:] 

#groups = ncar.groups.keys() # the groups keys are the same for all the files   odict_keys(['crs', 'era5fb', 'header_table', 'id_scheme', 'observations_table', 'observed_variable', 'station_configuration', 'station_configuration_codes', 'station_type'])
 

Merging = Merger()
""" Will create automatically the dictionary of the available data to consider """
Merging.InitializeData( datasets = available_data ) 

Merging.plot_styler()
Merging.PlotTimeDistribution()

variables = [85 ]
pressures = [100000, 500, 50000]


extracted_data = Merging.Extract_variable_pressure_time(variables= variables, pressures = pressures )            


for v in variables:
      for p in pressures:
            Merging.PlotTimeSeries( data_dic= extracted_data , variable = v , pressure = p )
      

merged = Merging.Merge(data = Merging.data )












# ncar_obs.variables = 


#analizer.findOriginalFile(dataset = 'ncar', primary_station = '0-20000-0-94711' )
#ncar_file_finder = analizer.ncar_file_finder

print('stop')


"""
['adjustment_id', 'advanced_assimilation_feedback', 'advanced_homogenisation', 'advanced_qc', 'advanced_uncertainty', 
'bbox_max_latitude', 'bbox_max_longitude', 'bbox_min_latitude', 'bbox_min_longitude', 'code_table', 'conversion_flag', 
'conversion_method', 'crs', 'data_policy_licence', 'date_time', 'date_time_meaning', 'exposure_of_sensor', 'latitude', 
'location_method', 'location_precision', 'longitude', 'numerical_precision', 'observation_duration', 'observation_height_above_station_surface', 
'observation_id', 'observation_value', 'observed_variable', 'original_code_table', 'original_precision', 'original_units', 'original_value', 
'processing_code', 'processing_level', 'quality_flag', 'report_id', 'secondary_value', 'secondary_variable', 'sensor_automation_status', 'sensor_id', 'source_id', 
'spatial_representativeness', 'traceability', 'units', 'value_significance', 'z_coordinate', 'z_coordinate_method', 'z_coordinate_type']


"""


"""


"""


