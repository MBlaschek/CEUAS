""" Merging the station configuration files """

import os,sys
import netCDF4 as nc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import argparse



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
                  #self.bufr_file_finder = ''
                  self.igra2_file_finder = ''
                  #self.ncar_file_finder = ''
                  self.bufr_file_finder = pd.read_csv( 'bufr_stationIds_filename.dat', sep='\t',  names = ['primary_id', 'secondary_id', 'file']) 
                  #self.igra2_file_finder = pd.read_csv('igra2_stationIds_filename.dat', 'r')
                  self.ncar_file_finder = pd.read_csv('ncar_stationIds_filename.dat' , sep='\t', names = ['primary_id', 'secondary_id', 'file'])   
                           
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
                  if dataset in ['era5_1', 'era5_1759', 'era5_1761' , 'era5_3188']:  ######## TODO !!!! 
                         return '-'      
                  
                  found_file = ''
                  if dataset == 'bufr':
                        df =  self.bufr_file_finder
                  elif dataset == 'ncar':
                        df =  self.ncar_file_finder
                                             
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
                        
                        

                        #print (  bufr_flag , bufr_stat,  bufr_lat, bufr_lon , bufr_start, bufr_end, bufr_file )
                        #print ( ncar_flag , ncar_stat,  ncar_lat, ncar_lon , ncar_start, ncar_end, ncar_file)
                          

                        if ncar_flag == '1':
                              file_ncar = self.findOriginalFile(primary_station=i, dataset='ncar')
                        if bufr_flag == '1':
                              file_bufr = self.findOriginalFile(primary_station=i, dataset='bufr')                                 
                        if igra2_flag == '1':                              
                              file_igra2 = self.findOriginalFile(primary_station=i, dataset='igra2') 
 
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
            self.df = {}            # storing the dafarames 
            self.DataDic = {}   # storing the info from the dataframe 
            
      def InitializeData(self, datasets = {'ncar' : '' , 'bufr' : '', 'igra2' : '' } ):
            """ Initialize dataset. 
                       Args:: datasets (dictionary where key+ dataset name e.g. bufr, igra2 etc. , and value is the path to the netCDF file """
            
            for k,v in datasets.items() :
                  self.df[k] = nc.Dataset(v)
                  self.DataDic[k] = {}
                  self.DataDic[k]['dateindex'] = self.df[k].variables['dateindex'][0,:]  # storing the dateindex 
                  self.DataDic[k]['observations_table'] = {}
                  self.DataDic[k]['observations_table']['date_time'] = self.df[k].groups['observations_table']['date_time']   # storing the date_time 
           
   
      def PlotTimeSeries(self):
            ncar = self.DataDic['ncar']['dateindex']
            bufr = self.DataDic['bufr']['dateindex']
            igra2 = self.DataDic['igra2']['dateindex']
            fig, ax = plt.subplots()            
            plt.scatter(ncar, np.full( (len(ncar)) , 1 ), label = 'ncar', color = 'blue')
            plt.scatter(bufr, np.full( (len(bufr)) , 2 ), label = 'bufr', color = 'lime')
            plt.scatter(igra2, np.full( (len(igra2)) , 3 ), label = 'igra2', color = 'cyan')
            plt.legend()
            #plt.grid()
            plt.xlabel('Observation Dates', fontsize =12)
            ax.ticklabel_format(useOffset=False, style='plain' )            
            plt.xticks(rotation = 45, fontsize = 8)
            plt.savefig('plots/test_timeseries.png', bbox_inches='tight', dpi = 200)


""" 
['primary_id', 'primary_id_scheme', 'record_number', 'secondary_id', 'secondary_id_scheme', 'station_name', 'station_abbreviation', 'alternative_name', 'station_crs', 'longitude', 'latitude',
'local_gravity', 'start_date', 'end_date', 'station_type', 'platform_type', 'platform_sub_type', 'operating_institute', 'operating_territory', 'city', 'contact', 'role', 'observing_frequency',
'reporting_time', 'telecommunication_method', 'station_automation', 'measuring_system_model', 'measuring_system_id', 'observed_variables', 'comment', 'optional_data', 'bbox_min_longitude',
'bbox_max_longitude', 'bbox_min_latitude', 'bbox_max_latitude', 'metadata_contact', 'metadata_contact_role']


['primary_id',  'secondary_id', 'station_name', 'alternative_name', 'station_crs', 'longitude', 'latitude',
'start_date', 'end_date',  'city', 
'observed_variables', ]
"""


make_summary = True

if make_summary:
      
      """ Read the station_configuration files """
      stat_conf_file = 'station_configuration_xxx.dat'
      bufr_sc   = stationConfigurationTable( file= 'station_configuration_xxx.dat'.replace('xxx', 'bufr' )  )
      bufr_sc.readDF()
      igra2_sc = stationConfigurationTable( file= 'station_configuration_xxx.dat'.replace('xxx', 'igra2' ) ) 
      igra2_sc.readDF()
      ncar_sc  = stationConfigurationTable( file= 'station_configuration_xxx.dat'.replace('xxx', 'ncar' )  )
      ncar_sc.readDF()
      
      era5_1_sc  = stationConfigurationTable( file= 'station_configuration_xxx.dat'.replace('xxx', 'era5_1' )  )
      era5_1_sc.readDF()
      era5_1759_sc  = stationConfigurationTable( file= 'station_configuration_xxx.dat'.replace('xxx', 'era5_1759' )  )
      era5_1759_sc.readDF()
      era5_1761_sc  = stationConfigurationTable( file= 'station_configuration_xxx.dat'.replace('xxx', 'era5_1761' )  )
      era5_1761_sc.readDF()
      era5_3188_sc  = stationConfigurationTable( file= 'station_configuration_xxx.dat'.replace('xxx', 'era5_3188' )  )
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

data = { 'ncar'    : '/raid60/scratch/federico/ConvertedDB/ncar/chuadb_windc_21946.txt.nc'   ,
               'igra2'   : '/raid60/scratch/federico/ConvertedDB/igra2/chRSM00021946-data.txt.nc' ,
               'bufr'     : 'chera5.21946.bfr.nc' ,
               'era5_1' : 'chera5.conv._21946.nc' }


#ncar  = nc.Dataset(ncar_test)
#igra2 = nc.Dataset(igra2_test)
#bufr   = nc.Dataset(bufr_test)

#ncar_di   = ncar.variables['dateindex'][0,:] # dateindex variable
#gra2_di  = igra2.variables['dateindex'][0,:] 
#bufr_di    = bufr.variables['dateindex'][0,:] 

#groups = ncar.groups.keys() # the groups keys are the same for all the files   odict_keys(['crs', 'era5fb', 'header_table', 'id_scheme', 'observations_table', 'observed_variable', 'station_configuration', 'station_configuration_codes', 'station_type'])
 

Merging = Merger()
Merging.InitializeData( datasets = data ) 

Merging.PlotTimeSeries()





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





