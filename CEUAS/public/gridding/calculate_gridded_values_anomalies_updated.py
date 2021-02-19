""" Module to calculate averages of monthly variables and anoalies,
      using stations inside different grid boxes.  """

import pandas as pd
import xarray as xr
import os,sys
from tqdm import tqdm
import numpy as np

from multiprocessing import Pool
from functools  import partial, partialmethod

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)



class Gridding():
    """ Main class to hold utilities to create gridded average monthly files and anomalies calculation """
    
    def __init__(self, out_dir = '', monthly_dir ='' , box = {'':''}, variable = '' ):
        
        self.out_dir = out_dir    # desired output directory 
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
        
        self.monthly_dir = monthly_dir # input directory with monthly files 
        #self.boxes = np.load(boxes, allow_pickle = True).item()
        self.std_plevs = [1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 85000, 92500, 100000]
        self.box_code = list(box.keys())[0]
        self.box = list(box.values())[0]
        self.variable = variable 
        
    def make_all_datetime(self, start_year=1900, end_year=2021):
        """ Create all date_time, on the day 15 of each month, at hours 00 and 12, 
              from start_year to end_year.
    
              In the output files, all these date_time will appear, eventually filled with
              nans if data not available. """
        
        date_time = []
        for y in range (start_year, end_year):
            for m in ['01' , '02' , '03' , '04' , '05' , '06' , '07' , '08' , '09', '10', '11', '12' ]:
                day = '15'
                for h in ['00','12']:
                    #DT = np.datetime64('1978-01-02T06:00:00.000000000')
                    timestamp = str(y) + '-' + m + '-' + day + 'T' + str(h) + ':00'
                    ts = np.datetime64(timestamp)
                    date_time.append(ts)
        date_time.sort()
        
        self.date_time = date_time         
        
    def find_file(self, primary):
        """ Find file in the monthly_files directory for a given primary id and given variable (e.g. temperature = ta)"""
        
        try:
            f = [f for f in os.listdir(self.monthly_dir) if primary in f and self.variable in f ][0]
            if not os.path.isfile(self.monthly_dir + '/' + f ):
                self.file = None
            else:
                self.file = self.monthly_dir + '/' + f
        except:
            self.file = None
            
            
    def make_empty(self):
        ''' Create empty structure of the results dictionary '''
        
        var = self.variable
        
        results = {}
        results['res'] = {}
        results['res']['empty'] = {}
        
        results['res']['empty'][var+'_average']         = ([np.nan] * len(self.date_time) * 16 )  
        results['res']['empty'][var+'_average_bias'] = ([np.nan] * len(self.date_time) * 16)
        results['res']['empty'][var+'_anomaly']        = ([np.nan] * len(self.date_time) * 16)
        results['res']['empty'][var+'_anomaly_bias']= ([np.nan] * len(self.date_time) * 16)
        results['res']['empty']['plev']             = ( self.std_plevs * len(self.date_time) )   
    
        results['res']['empty']['time'] = []
        for dt in self.date_time:
            results['res']['empty']['time'].extend( [dt] * 16 )  
                
        return results          


    def load_df(self):
        """ Load the dataframe from a monthly file. 
         ['index','lat','lon','plev','sonde_type','ta','ta_bias','ta_dep','ta_glob','ta_num_obs',
          'ta_num_obs_bias','ta_num_obs_dep','ta_num_obs_glob','ta_std_dev','ta_std_dev_bias',
          'ta_std_dev_dep', 'ta_std_dev_glob', 'time']
        """
        try:
            df = xr.open_dataset(self.file , engine = 'h5netcdf', decode_times = True ).to_dataframe()
            return df 
        except:
            print('*** ERROR: no available data in file: ' , self.file )
            return None
        
    def reduce_dataframe(self, obs_tab, min_days=10, hour = ''):
            """ Reduces the initial dataframe(obs_tab) removing variables and values that are not needed """

            obs_tab['year']    = pd.arrays.DatetimeArray (obs_tab['time'].values[:] ).year
            obs_tab['month'] = pd.arrays.DatetimeArray (obs_tab['time'].values[:] ).month
            obs_tab['hour']    = pd.arrays.DatetimeArray (obs_tab['time'].values[:] ).hour
            
            df_red = obs_tab.loc [ (obs_tab[self.variable + '_num_obs_dep'] > min_days) 
                                            &  (obs_tab['hour'] == hour ) ] 
            return df_red
        
        
    def calculate_anomaly(self, obs_tab, years_window=20, min_months=10, month ='', year ='', hour ='' ):
        """ main utility to calculate the monthly anomalies """
    
        """ Load the dataframe from a monthly file. 
         ['index','lat','lon','plev','sonde_type','ta','ta_bias','ta_dep','ta_glob','ta_num_obs',
          'ta_num_obs_bias','ta_num_obs_dep','ta_num_obs_glob','ta_std_dev','ta_std_dev_bias',
          'ta_std_dev_dep', 'ta_std_dev_glob', 'time'] """
          
        averages, averages_bias, anomalies, anomalies_bias, plevels = [],[],[],[],[]
        var = self.variable
        """ For each plevel, need to extract the df of the previous_climatology (depending on year_window),
              and the current one. """
        
        for p in self.std_plevs:
                # df of the previous 20 years of which I calculate the climatology 
                previous_climatology = obs_tab.loc [ (obs_tab['plev']==p) 
                                   & (obs_tab['month'] == month)
                                   & (obs_tab['hour'] == hour)
                                   & (obs_tab['year'] > year - years_window )
                                   & (obs_tab['year']  <= year - 1) ]
                            
                """ This reduced df contains only the monthly averages for the particular variable, pressure, time, 
                      that are calculated using more than the minimum required number of days,
                      see self.reduce_dataframe(min_days=10) """
    
                if len(previous_climatology) > min_months:
                    # dataframe of the year-month of which I want to calculate the anomaly  
                    current_df = obs_tab.loc [ (obs_tab['plev']==p) 
                                       & (obs_tab['month'] == month)
                                       & (obs_tab['hour'] == hour)
                                       & (obs_tab['year'] == year) ]
    
                    """ Reading the values of the current entry of the dataframe """
                    try:
                        
                        average = current_df[var + '_dep'].values[0]
                        average_bias =  current_df[var + '_bias'].values[0]
                    except:
                        average = current_df[var + '_dep'].values
                        average_bias =  current_df[var + '_bias'].values    
                    
                    """ Calculating the values of the previous climatology entries """                    
                    climatology_average = np.mean(previous_climatology[var + '_dep'].values)
                    climatology_average_bias = np.mean(previous_climatology[var + '_bias'].values)
                    
                    anomaly = average - climatology_average
                    anomaly_bias = average_bias - climatology_average_bias
    
                else:                
                    average, average_bias, anomaly, anomaly_bias = np.nan, np.nan, np.nan, np.nan 
                    
                averages.append(average) # no bias correction average
                averages_bias.append(average_bias) # bias corrected average
                anomalies.append(anomaly) # no bias correction anomaly 
                anomalies_bias.append(anomaly_bias) # bias corrected anomaly
                plevels.append(p) # pressure level 
                    
        return averages, averages_bias, anomalies, anomalies_bias, plevels 
                              
                                                       
    def calculate_box(self):
        """ Main utility to run the calculation over each gridded box.
              The loop is external to use paraller running. """
    
        bbox = self.box 
        var = self.variable 
        results = self.make_empty( ) # make a default empty result 
        self.results = results 
        
        #bbox = self.boxes[box] # analyzing a specific box
        if abs(bbox['lat'][0]) > 90 or abs(bbox['lon'][0]) > 180 :
            print('*** WARNING: wrong lat/lon' , bbox['lat'][0] , ' ' , bbox['lon'][0] )
            return False
    
        """ Check if there are station in the given box """
        if not bbox['files']:
            print('*** WARNING: no stations found inside this box *** ' , bbox )
            pass
        
        else:        
            #stations = [s for s in bbox['files'] if '10393' in s] 
            """ Loop over stations in the box """
            for station in bbox['files']:          # loop over each station file (if any)
            #for station in stations:          # loop over each station file (if any)
            
                print('*** Processing station *** ' , station , str(bbox['files'].index(station)) + '/' + str(len(bbox['files'])) , '  of box: ' , bbox  )
                station = station.replace("b'" , ''  ).replace("'",'')
                dummy_load = self.find_file(station) 
                
                if self.file:   # load dataframe if file exists 
                    df = self.load_df()
                    
                    if not isinstance(df, pd.DataFrame): # check if valid DataFrame exists
                        print('*** ERROR: invalid dataframe ***' , station )
                        continue

                    if df.empty: # check if DataFrame not empty
                        continue 
                    
                    for h in [0,12]:    # loop over the two possible hours 
                        red_df = self.reduce_dataframe(df, min_days=10, hour = h)
                        
                        if red_df.empty : 
                            continue
    
                        results['res'][station]  = {}
                        
                        results['res'][station][var + '_average']    , results['res'][station] [var + '_average_bias']  = [], []
                        results['res'][station][var + '_anomaly'] , results['res'][station][ var + '_anomaly_bias'] = [], []
                        
                        results['res'][station]['plev'] = []
                        results['res'][station]['time'] = []
                        
                        
                        all_dt = np.unique(df['time'])
                        for dt in self.date_time: # all possible date_time starting from 1900 until 12/2020
                            dt = pd.to_datetime(dt)
                            
                            results['res'][station]['time'].extend ( [dt] * 16 )      
    
                            if dt in all_dt : #  available date_time for this station, will calculate values
                                average, average_bias, anomaly, anomaly_bias, plevels = self.calculate_anomaly(red_df, years_window=20, min_months=10,
                                                                                                               month =dt.month , year = dt.year , hour = h )
                           
                                for p, a, ab, an, anb in zip (plevels, average, average_bias, anomaly, anomaly_bias ):
                                    results['res'][station][var + '_average'].append(a)  
                                    results['res'][station][var + '_average_bias'].append(ab)
                                    results['res'][station][var + '_anomaly'].append(an)
                                    results['res'][station][var + '_anomaly_bias'].append(anb)
                                    results['res'][station]['plev'].append(p)
                                
                                                                    
                            else:
                                results['res'][station][var + '_average']          .extend ( [np.nan] * 16 )  
                                results['res'][station][var + '_average_bias']  .extend ( [np.nan] * 16 )
                                results['res'][station][var + '_anomaly']        .extend ( [np.nan] * 16 )
                                results['res'][station][var + '_anomaly_bias'].extend ( [np.nan] * 16 )
                                results['res'][station]['plev'].extend ( self.std_plevs )   
                                
    
                else:
                    print('*** WARNING: file not found in the mothly directory created with the CDS database, skipping station ***' , station )
                    pass
                
                self.results = results 

                
        return True
    
    def  make_write_box_df(self):
        """ Combine data from different stations and save netcdf file output """
        var = self.variable        
        box = self.box 
        results = self.results
        
        lat, lon = box['lat'][0] , box['lon'][0]
        size      = abs( box['lat'][0] - box['lat'][1] )
        Lat , Lon = lat + size/2 , lon + size/2


        output_columns =  [ var+'_average' , var+'_anomaly' , var+'_average_bias' , var+'_anomaly_bias' ] # Variables to be used in the netCDF output files
        
        """ Creating the output dataframe """
        output_df = pd.DataFrame(columns = ['plev' , 'time', 'lat' , 'lon'] + output_columns )      
        
        stations = [ s for s in results['res'].keys() if 'empty' not in s ]

        """ Here, loop over the entire length of calculated  results,
              which are already sorted by time from 1900 to 2021, and by plevels.
              If no stations are available, use dummy nan results. """
        for column in  output_columns : # loop over each desired column
            # need to add a sum!!! 
            means = []
            if len(stations) >=1 :             
                for v in range( len(results['res']['empty']['time'] )  ) : # total length of the results for each station 
                    
                    values = [ results['res'][s][column][v] for s in stations ]
                    average = np.nanmean(values)
                    
                    if not average:
                        average = np.nan 
                    try:
                        means.append(average[0])
                    except:
                        means.append(average)

            else:
                means = results['res']['empty'][column]               
            
            output_df[column] = means 
        
        output_df['plev'] = results['res']['empty']['plev']
        output_df['time'] = results['res']['empty']['time']
       
        output_df['lat']            = ( [Lat] * len(results['res']['empty']['time'] ) )
        output_df['lon']           = ( [Lon] * len(results['res']['empty']['time'] ) )
        #output_df['grid_size']  = ( [size] * len(results['res']['empty']['time'] ) )

        """ Saving array with stations """
        if  len(stations) >=1:
            stat = '_'.join(stations)
        else:
            stat = 'None'
            

        
        #output_df['primary_ids']      = ( [stat] * len(results['res']['empty']['time'] ) )

        xarray = output_df.to_xarray()
        xarray['lat'].attrs['stations'] = stat 
        xarray['lon'].attrs['stations'] = stat 
        
        if len(stat) > 50:
            stat = stat.split('_')[0] + stat.split('_')[1] + stat.split('_')[1] + '_others'
        box_name = self.box_code + '_' + 'lat_lon_size_' + str(Lat) + '_' + str(Lon) + '_' + str(size) + '_' + stat        
        dummy = xarray.to_netcdf ( self.out_dir + '/' +  box_name + '.nc' , mode = 'w' )
        print('*** Written output file: ' , box_name )
        # xarray =output_df[43500:43600]
        # check::: 43536  1000         2018-05-15 12:00:00  []  
    
    
    def process_box(self):
        """ Wrapper for the operation to be perfomed onto each grid box """
        dummy = self.make_all_datetime()
        check = self.calculate_box()
        if check:
            dummy = self.make_write_box_df()
            
        return 




""" Initialize relevant variable for input/output """
monthly_file_dir = '/raid60/scratch/federico/MONTHLY_FEB2021/wind_speed/'
out_dir = '/raid60/scratch/federico/GRIDDED_FILES_FEB2021/wind_speed'

out_dir = '/raid60/scratch/federico/PROVA_WINDSPEED_FEB2021/'

boxes_file = 'stations_per_box_size_10.npy'
variable = 'wind_speed'

""" Numpy file containing the grid boxes """
boxes = np.load(boxes_file, allow_pickle = True).item()

remove = [] 


""" # for analsing one specific box 
for k in boxes.keys():
    lat, lon, latm, lonm = boxes[k]['lat'][0] , boxes[k]['lon'][0], boxes[k]['lat'][1], boxes[k]['lon'][1]
    if latm > 90 or lat < -90 or lon < -180 or lonm > 180 or  k != '24_20':
        remove.append(k)
"""  

    
print('*** Cleaning wrong boxes ***')

for r in remove:
    del boxes[r]




if __name__ == '__main__':

    MakeGrid = Gridding(out_dir = out_dir , monthly_dir = monthly_file_dir , variable = variable )

    """ List of all existing boxes to process, and removed the already processed ones  """
    #all_boxes = list(boxes.keys())
    
    processed_boxes = [f.split('_lat_lon_')[0] for f in os.listdir(out_dir)]
    processed = [b for b in list(boxes.keys()) if b in processed_boxes ]
    for p in processed:
        del boxes[p]
    
    print('*** Must process ' , len(boxes) , '  boxes ' ) 

    """ Select POOL=True for parallel processing """
    POOL = True

    if POOL:
        
        def run(out_dir,monthly_dir, variable, box):
            try:             
                MakeGrid = Gridding(out_dir = out_dir , monthly_dir = monthly_file_dir,  variable = variable, box = box)
                MakeGrid.process_box() 
            except:
                a = open('Failed_boxes.dat' , 'a')
                a.write(list(box.keys())[0] + '\n')
        """ List of boxes to process (list of dictionaries) """
        bboxes = [{b:boxes[b]} for b in boxes.keys() ]        
        
        #for b in bboxes:
        #    run(out_dir = out_dir , monthly_dir = monthly_file_dir, variable = variable, box = b)
        func = partial(run, out_dir , monthly_file_dir, variable )        
        p = Pool(40)
        out = p.map(func, bboxes)

    else:  
        for b in boxes.keys():
            box_dic = {b: boxes[b]}
            #box_dic = {'11_16': {'lat': [-80.0, -70.0], 'lon': [-30.0, -20.0], 'files': ["b'0-20000-0-89022'"]}}
            MakeGrid = Gridding(out_dir = out_dir , monthly_dir = monthly_file_dir, box = box_dic , variable = variable )
            a =  MakeGrid.process_box()    




