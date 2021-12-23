"""
    Module: analyze_inventory_functions.py
    
    Author:: Ambrogi Federico

"""
import os, sys
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

import numpy as np
import re 
import subprocess
import geopy.distance
import glob

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)

from multiprocessing import Pool
from functools  import partial

from eccodes import *


class Utils():
    """ Class holding a series of calculation utilites """
    def __init__(self, treshold = 10):
        """ threshold: maximum allowed distance in Km for points to be considered compatible """
        self.treshold = treshold
        
        # station_configuration columns and mapping
        self.inventory_to_statConf = { 'WIGOS_calc' : 'primary_id' ,
                'file_statid' : 'secondary_id',
                'station_name' : 'station_name',
                'latitude':'latitude',
                'longitude':'longitude',
                'start_date': 'start_date',
                'end_date':'end_date',
                'city':'city',
                }
        self.statConf_to_inventory = dict( zip(  self.inventory_to_statConf.values(),  self.inventory_to_statConf.keys()) )

        return None
    
    def distance(self, lat1,lon1,lat2,lon2):
        """ Return the Harvesine distance between two points """
        d = geopy.distance.geodesic( (lat1,lon1), (lat2,lon2) ).km
        return d     

    def degMinSec_to_decimal(self, coord):
        """ Converts lat and lon [lists] from degrees-minute-seconds to decimal """
        
        def dms2dec(dms):
            """ Split strings like -75 14 00N into degree, minute and seconds,
            convert into decimal with properly assigned sign """
            
            sign = -1 if re.search('[swSW]', dms) else 1
            
            dms = re.sub ( '[wWsSnNeE]', '', dms  )
            
            dms = [ d for d in dms.split(' ') if d ]
            #d = dms.split(' ')[0]
            #m = dms.split(' ')[1]
            #s = dms.split(' ')[2]
            
            d = dms[0]
            m = dms[1]
            s = dms[2]

            dec = float(d) + float(m)/60 + float(s)/3600
            
            return sign*round(dec, 2)
        
        coord_dec = map(dms2dec, coord)        

                
        return coord_dec
    
    def coord_consistency(self, lats='',lons='', threshold= 10 ):
        """ Check if the absolute variation in time of the coordinates is below a certain threshold in Km """
            
        check = self.distance( min(lats) , min(lons), max(lats), max(lons)  ) < threshold 
        return check
    
    
    
class Analyze():
    
    def __init__(self, dist_limit = 30, inv = '', data = '', utils='', cities = ''):
        """ parameters ::
                     dist_limit : maximum distance allowed for station matching
                     inv :: dataframe of a single the inventory
                     data :: data extracted from a single file from one of the data sets """
        
        self.dist_limit = dist_limit
        self.inv = inv
        self.data = data
        self.found_inv= { self.inv.name:{} }
        self.utils = utils 
        
        self.wigos = { 'IGRA2':'0-20300-0-'   ,
                                'WBAN': '0-20500-0-'  ,
                                'CHUAN':'0-20400-0-'  }
        
        self.best_match = ''
        self.cities = cities 
        
    def Make_WIGOS(self, df):
        """ Creates a WIGOS id for the file given the best matching inventory data """
        #df = self.matchingInv_dist_id
        # case OSCAR ds 
        if self.inv.name == 'OSCAR':
            
            if len(df) == 1: # only one matching, keep this 
                wigos = df['WIGOS'][0]
            else:
                # TO DO handle special case?                 
                wigos = df['WIGOS'][0]

                
        else:
            prefix = self.wigos[self.inv.name]
            if len(df) == 1: # only one matching, keep this 
                wigos = prefix + str(df['station_id'][0])
            else:
                # TO DO handle special case? 
                wigos = prefix + str(df['station_id'][0])
                 

        #print(0)                
        return wigos
        
        
    #def Distance(self, lat1,lon1,lat2,lon2):
    #  """ Return the Harvesine distance between two points """
    #    d = geopy.distance.geodesic( (lat1,lon1), (lat2,lon2) ).km
    #    return d 
        
    def FindCity(self):
        """ Find the closest city within a 100 km radius """
        
        lat, lon = self.data.lats[0] , self.data.lons[0]
        cities = self.cities
        cities = cities.loc [  (abs(cities['Latitude'] - lat) < 5 ) & (abs(cities['Longitude'] - lon) < 5) ] 
        if not cities.empty:
            cities = cities.reset_index()
            distances = [ round( self.utils.distance( lat, lon , cities['Latitude'][i] , cities['Longitude'][i] ), 1 ) for i in range(len(cities['Latitude']) ) ]
            ind_min = distances.index(min(distances)) # index of the closest city 
            
            city, distance = cities.iloc[ind_min] , distances[ind_min]
            self.best_match['city'] = city.City
            self.best_match['city_dist_km'] = distance
            self.best_match['city_lat'] = city.Latitude
            self.best_match['city_lon'] = city.Longitude
        
        else:
            self.best_match['city'] = 'None'
            self.best_match['city_dist_km'] = 'None'
            self.best_match['city_lat'] = 'None'
            self.best_match['city_lon'] = 'None'            
               
         

    def AddDistance(self):
        """ Calculates the distances between the point and the stations in the inventory """
        print("*** Adding distance ***")
        
        # skimming the inventory with close by stations only, to avoid useless calculation of distances
        df = self.inv
        red_inv = df.loc [  (abs(df['latitude'] - self.data.lats[-1]) < 2 ) & (abs(df['longitude'] - self.data.lons[-1]) < 2) ] 
        
        # i.e. at least one station found whihc si reasonably close
        if len(red_inv) >0:
            
            red_inv = red_inv.reset_index()
            
            distances = [ round( self.utils.distance( self.data.lats[-1], self.data.lons[-1] , red_inv['latitude'][i] , red_inv['longitude'][i] ), 1 ) for i in range(len(red_inv['latitude']) ) ]
    
            red_inv['distance_km'] = distances
            red_inv = red_inv.sort_values(by =['distance_km'])
        else:
            distances = []
            red_inv['distance_km'] = 'NAN'
            
        self.red_inv = red_inv 
        red_inv.name = self.inv.name     
        
        return distances
        
    def Closest(self):
        """ Find all the stations that are within the distance limit,
              or have the same station identifier,
              """
        print("*** Find Closest ***")

        if len(self.red_inv["distance_km"]) >0:
            
            closest_indices = np.searchsorted( self.red_inv["distance_km"], self.dist_limit )
            inv_closest = self.red_inv[0:closest_indices]
            inv_closest.name = self.inv.name 
            inv_closest = inv_closest.reset_index(drop=True)
            inv_closest = inv_closest.drop(columns=['index'])
            
            self.found_inv[self.inv.name ]['matching_dist'] = inv_closest
        else:
            inv_closest = pd.DataFrame()
            inv_closest.name = 'NONE'
            self.found_inv['NONE'] = {}
            self.found_inv['NONE']['matching_dist'] = inv_closest
            return 0
            
    def MatchingId(self):
        """ Check if there is a station id matching the input station id data.
             Note that statId is a list, most frequently with one element only,
             but might contain more if more ids are found in a single file e.g. for NCAR dataset.
             """
        
        print("*** Find MatchingId ***")
        
        if self.inv.name in ['CHUAN', 'WBAN']:
            stats = 'WMO_id'
     
        else:
            stats = 'station_id'
            
        if isinstance(self.data.statid, pd.Series):
            if self.data.dataset == 'igra2':
                try:
                    s = [eval(self.data.statid[0])]
                except:
                    s = [self.data.statid[0]]
                    
                    
            else:
                s = eval(self.data.statid[0])
            
            s = [str(s) for s in s ]
            stat_ids = [ f.split(':')[1]  if ':' in f else f for f in s ]
                
        df = self.inv.loc[ self.inv[stats].isin( stat_ids ) ] 
        
        # test the missing latitude sign for ERA5 1759 inventory 
        if len(df) ==0 and self.inv.name=='WBAN' and self.data.dataset:
            df_alt = self.inv.loc[ self.inv['station_id'].isin( stat_ids ) ]
            if len(df_alt) >0:      
                df = df_alt
                print(0)
            
            
            
        inv_matchingId = df
        inv_matchingId.name = self.inv.name 
        inv_matchingId = inv_matchingId.reset_index(drop=True)        
        
        lats = self.data.lats
        lons = self.data.lons
        
        distances = [ round( self.utils.distance( lats[-1], lons[-1] , inv_matchingId['latitude'][i] , inv_matchingId['longitude'][i] ), 1 ) for i in range(len(inv_matchingId['latitude']) ) ]

        inv_matchingId['distance_km'] = distances
        inv_matchingId['distance_km_minusLat'] = ['' for i in range(len(inv_matchingId))]  # must add to keep dataframes compatible in all cases
        
        if len(df) >0 and self.inv.name=='WBAN' and self.data.dataset:
            
            distances_swap = [ round( self.utils.distance( -lats[-1], lons[-1] , inv_matchingId['latitude'][i] , inv_matchingId['longitude'][i] ), 1 ) for i in range(len(inv_matchingId['latitude']) ) ]   
            inv_matchingId['distance_km_minusLat'] = distances_swap
            
        self.found_inv[self.inv.name ]['matching_id'] = inv_matchingId
        
        #if not inv_matchingId.empty:
        #    print('check')
               
    def MakeIgraInventory(self):
        """ Special case: for IGRA2 inventory only need to convert the inventory entries to the format of the other inventories """
        igra2 = self.inv
        df = pd.DataFrame()
        
        dic = { 'max_date': max_date , 
                'min_date':min_date, 
                'statid':statids, 
                'lats': [lats] , 
                'lons': [lons], 
                'file': file , 
                'file_path': self.file , 
                'db': self.dataset,
                'variables': str(variables) }

    def Combine(self):
        """ Concatenating the df for each inventory and removing duplicated entries """
        try:
            id_match = self.found_inv[self.inv.name ]['matching_id']
        except:
            id_match = pd.DataFrame()
        try:
            dist_match = self.found_inv[self.inv.name ]['matching_dist']
        except:
            dist_match = pd.DataFrame()
        
        if not ( id_match.empty and dist_match.empty ):  # case where I have both an ID match and a dist match 
            res = pd.concat ( [id_match, dist_match], ignore_index=True) 
            res = res.drop_duplicates()
        elif ( id_match.empty and not dist_match.empty ):
            res = dist_match
        elif ( not id_match.empty and dist_match.empty ):
            res = id_match
            
        else:
            res = id_match # empty case
            for c in res.columns:
                res[c] = [np.nan]            
            res['WIGOS_calc'] = ['None'] 
            res['inventory'] = self.inv.name 
            self.best_match = res 
            #return
                    
        wigos = self.Make_WIGOS(res)
        res['inventory'] = self.inv.name         
        res['WIGOS_calc'] = wigos 
        
        self.best_match = res 
                
        print(0)
           
        
        
    
class Data():
    """ Main class to hold utilities for reading original input data from 
         - ERA5 1,2,1759,1761,3188 
         - IGRA2
         - NCAR
         - BUFR 
         """
    
    def __init__(self, dataset='', utils = ''):
        
        self.dataset = dataset
        self.utils = utils 
        """
        self.datasets = { 'era5_1': '/mnt/users/scratch/leo/scratch/era5/odbs/1' ,
                                   'era5_2': '/mnt/users/scratch/leo/scratch/era5/odbs/2',
                                   'era5_3188': '/mnt/users/scratch/leo/scratch/era5/odbs/3188',
                                   'era5_1759': '/mnt/users/scratch/leo/scratch/era5/odbs/1759',
                                   'era5_1761': '/mnt/users/scratch/leo/scratch/era5/odbs/1761',
                                   
                                   'bufr': 'mnt/users/scratch/leo/scratch/era5/odbs/',                                   
                                   
                                   'ncar': '',
                                   'igra2': '',

                                   }    
        """
        self.odb_to_cdm = { 1    : 117 , # geopotential
                                   2    : 85          , # temperature K
                                   3    : 104        , # uwind m/s , upper air u component 
                                   4    : 105        ,  # vwind m/s
                                   7    : 39          ,  # specific humidity
                                   
                                  111 : 106    , # degree (angle) , wind from direction 
                                  112 : 107    , # m/s , wind force 
                                  29   : 38      , # relative humidity in %
                                  59    : 36      , # dew point (available in ODB files )
                                 
                          }
        
        
    def read_file(self, file):
        """ Wrapper for the functions to read input file from the datasets. 
             If a csv file exists holding the needed information, it is read.
             Otherwise it reads the original file with the appropriate functions. """
        
        # for igra2, file is indeed a pandas series 
        if self.dataset == 'igra2':
            self.file = file.station_id_igra
            self.series = file 
            
        else:
            self.file = file  

  
        # setting output directory holding csv files 
        out_dir = 'temp_data/' + self.dataset
        self.out_dir = out_dir 
        
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
            
        f = self.out_dir + '/' + self.file.split('/')[-1] + '.csv'
        
        if os.path.isfile( f ):         # read file if csv already exists
            df = pd.read_csv(f,  sep='\t', converters={'statid': str}  )
            
            self.statid = df['statid'].astype('str')
            if type(eval(df['lats'].astype('str').values[0])) == list:
                
                self.lats = eval(df['lats'].astype('str').values[0]) 
                self.lons = eval(df['lons'].astype('str').values[0]) 
            else:
                self.lats = [eval(df['lats'].astype('str').values[0]) ]   
                self.lons = [eval(df['lons'].astype('str').values[0]) ]
                
            self.min_date = df['min_date'].values[0]
            self.max_date = df['max_date'].values[0]
            self.consistent_coord = df['consistent_coord'].values[0]
            self.variables = df['variables'].values[0]
        
        else: # read from original files
            if 'era5_' in self.dataset :
                self.read_odb()
            elif 'ncar' in self.dataset:
                self.read_ncar()
            elif 'bufr' in self.dataset:
                self.read_bufr()
            elif 'igra2' in self.dataset:
                self.read_igra2()
        
    def read_odb(self):
        """ Extract data from a single ODB file 
        # odb headre -i file """
        
        file = self.file 
        #Extract unique lat and lon
        comm = "odb sql 'select distinct lat,lon'  -i FILE --no_alignment ".replace('FILE', file)  # ODB command in sql syntax 
        proc = subprocess.Popen(comm, stdout=subprocess.PIPE , shell=True)
        b = proc.stdout.read().decode('utf-8').split('\n')
        
        lats, lons  = [],[]
        for i in range (1,len(b)):
            if not b[i]:
                continue
            lats.append(  round(float(b[i].split('\t')[0])  , 3 ) )
            lons.append( round(float(b[i].split('\t')[1]) , 3 ) )
            #print(b[i])
            
        self.consistent_coord = True
        # define a max distance
        
        distances = [ round(self.utils.distance( lats[i], lons[i] , 
                                        lats[0], lons[0] , ), 1 ) for i in range(len(lats) ) ]         
        
        # maximum 20 km distance 
        if max(distances) > 20:
        
            #Extract all lat and lon if I have more than one value   
            # This is necessary to quantify inconsistent lat/lon 
            comm = "odb sql 'select lat,lon'  -i FILE --no_alignment ".replace('FILE', file)  # ODB command in sql syntax 
            proc = subprocess.Popen(comm, stdout=subprocess.PIPE , shell=True)
            b = proc.stdout.read().decode('utf-8').split('\n')
    
            lats_all, lons_all  = [],[]
            for i in range (1,len(b)):
                if not b[i]:
                    continue
                lats_all.append(  round(float(b[i].split('\t')[0])  , 3 ) )
                lons_all.append( round(float(b[i].split('\t')[1]) , 3 ) )
                #print(b[i])
                    
            ### TODO -> implement the funtion that does this 
            # group by equal lat/lon values 
            df = pd.DataFrame( {'lat': lats_all , 'lon' : lons_all } )
            grouped = df.groupby(df.columns.tolist(), as_index=False).size()
            # add distance column
            distances = [ round(self.utils.distance( grouped['lat'].values[i], grouped['lon'].values[i] , 
                                             grouped['lat'].values[0], grouped['lon'].values[0] ), 1 ) for i in range(len(grouped['lat'].values) ) ] 
    
            grouped['distance_km'] = distances
            
            compatible_coord = grouped.loc[grouped['distance_km'] <= self.utils.treshold ]
            frac = len(compatible_coord) / len(lats_all)
            if frac >= 0.99:
                self.consistent_coord = True
                pass
            else:
                self.consistent_coord = False
                
            grouped.to_csv( self.out_dir  + '/' + self.file.split('/')[-1] + '_coordinates.csv', sep = '\t' )        
        
                
        self.lats = lats
        self.lons = lons

        # extract statid 
        comm = "odb sql 'select distinct statid'  -i FILE --no_alignment ".replace('FILE', file)   
        proc = subprocess.Popen(comm, stdout=subprocess.PIPE , shell=True)
        b = proc.stdout.read().decode('utf-8').split('\n')
        statids = [ eval(c) for c in np.unique( b[1:-1] ) ] 
        
        # extract distinct variables
        comm = "odb sql 'select distinct varno'  -i FILE --no_alignment ".replace('FILE', file)   
        proc = subprocess.Popen(comm, stdout=subprocess.PIPE , shell=True)
        b = proc.stdout.read().decode('utf-8').split('\n')
        variables = [ eval(c) for c in np.unique( b[1:-1] ) ] 
        variables = [v for v in variables if v in self.odb_to_cdm.keys() ]
        variables = [self.odb_to_cdm[v] for v in variables]
        variables.sort()
        
        # extract min date, max date 
        comm = "odb sql 'select distinct MIN(date) '  -i FILE --no_alignment ".replace('FILE', file)   
        proc = subprocess.Popen(comm, stdout=subprocess.PIPE , shell=True)
        b = proc.stdout.read().decode('utf-8').split('\n')
        min_date = b[1]
        
        comm = "odb sql 'select distinct MAX(date) '  -i FILE --no_alignment ".replace('FILE', file)   
        proc = subprocess.Popen(comm, stdout=subprocess.PIPE , shell=True)
        b = proc.stdout.read().decode('utf-8').split('\n')
        
        max_date = b[1]
        
        dic = { 'max_date': max_date , 
                'min_date':min_date, 
                'statid':str(statids), # needed in case multiple ids 
                'lats': [lats] , 
                'lons': [lons], 
                'file': file , 
                'file_path': self.file , 
                'db': self.dataset,
                'variables': str(variables),
                'consistent_coord': str(self.consistent_coord) }
        
        # storing the 
        pd.DataFrame(dic).to_csv( self.out_dir  + '/' + self.file.split('/')[-1] + '.csv', sep = '\t' )        
        self.statid = statids
        self.lats = lats
        self.lons = lons 
        self.min_date = min_date 
        self.max_date = max_date
        self.variables = dic['variables']

        
    def check_consistent_coords(self, lats, lons):
        """ Given the entire list of lats and lons check if the values are consistent.
        If not, picks the most frequently occurring pairs, and checks how much of the data is compatible.
        If more than 99% of the data has similar coordinates, the file is considered correct.
        Otherwise it cannot be identified since coordinates are inconsistent in time."""
        
        # group by equal lat/lon values 
        df = pd.DataFrame( {'lat': lats , 'lon' : lons } )
        grouped = df.groupby(df.columns.tolist(), as_index=False).size()
        # add distance column
        distances = [ round(self.utils.distance( grouped['lat'].values[i], grouped['lon'].values[i] , 
                                         grouped['lat'].values[0], grouped['lon'].values[0] ), 1 ) for i in range(len(grouped['lat'].values) ) ] 

        grouped['distance_km'] = distances
        
        compatible_coord = grouped.loc[grouped['distance_km'] <= self.utils.treshold ]
        frac = len(compatible_coord) / len(lats)
        if frac >= 0.99:
            self.consistent_coord = True
            pass
        else:
            self.consistent_coord = False
            
        print(0)
        #grouped.to_csv( self.out_dir  + '/' + self.file.split('/')[-1] + '_coordinates.csv', sep = '\t' )  
        
        
    def read_bufr(self):
        """ Extract data from a single BUFR file """
        
        bufrfile = self.file 
        
        vdict={'111':'windDirection','112':'windSpeed','1':'geopotentialHeight',
               '2':'airTemperature','59':'dewpointTemperature','29':'relativeHumidity'}
                
        # check if valid bufr file, preliminary loop
        with open(bufrfile) as f:
            cnt = 0
            # loop over the messages in the file
            while 1:
                # get handle for message
                bufr = codes_bufr_new_from_file(f)
                if bufr is None:
                    break
                codes_release(bufr)
                cnt+=1
    
        if cnt==0:
            return None
        
        else:
            # storing here the complete data from file 
            data = {}
            var = ['latitude' , 'longitude', 'blockNumber', 'stationNumber', 'date', 'variables']
            for v in var: 
                data[v] = []

            # looping through buffer 
            with open(bufrfile) as f:
                # loop over the messages in the file
                for i in range(cnt):
                    # get handle for message
                    bufr = codes_bufr_new_from_file(f)
                    if i not in [0,cnt-1]:
                        codes_release(bufr)
                        continue
                    else:
                        codes_set(bufr, 'unpack', 0)
                        # get all the timePeriods
                        iterid = codes_bufr_keys_iterator_new(bufr)
        
                        # loop over the keys
                        #if codes_get_array(bufr,'dataSubCategory')[0] not in [91,101]:
                        #    codes_release(bufr)
                        #    continue
                        
                        # extracting date
                        while codes_bufr_keys_iterator_next(iterid):     
                            keyname = codes_bufr_keys_iterator_get_name(iterid)
                            
                            if keyname == 'localYear':
                                year = codes_get_array(bufr,keyname)[0]
                            if keyname == 'localMonth':
                                month = codes_get_array(bufr,keyname)[0]
                            if keyname == 'localDay':
                                day = codes_get_array(bufr,keyname)[0]
                        #datum='{}-{:0>2}-{:0>2}'.format(year, month, day  ) this is the format of the stat conf 
                        datum='{}{:0>2}{:0>2}'.format(year, month, day  )
                        
                        if i == 0:
                            data['start_date'] = datum
                        if i== (cnt-1):
                            data['end_date'] = datum
                            
                        # delete the key iterator
                        codes_bufr_keys_iterator_delete(iterid)
        
                        for v in ['latitude' , 'longitude', 'stationNumber', 'blockNumber']:
                            data[v].append(codes_get(bufr, v))
                            
                        #lat = codes_get(bufr, "latitude")
                        # odetype@hdr == 35
                        # expver@desc == ERA40
                        # ECMWF ? 
        
                        for varno in ['111','112','1','2','59','29']:
                            try:
                                a = codes_get_array(bufr, vdict[varno] )
                                data['variables'].append ( self.odb_to_cdm[int(varno)] ) 
                            except:
                                pass
        
                        codes_release(bufr)
        
        print(0)

        lats, lons = [round(l,3) for l in data['latitude']] , [round(c,4) for c in data['longitude']]
        df = pd.DataFrame( {'lat': lats , 'lon' : lons } )
        df = df.drop_duplicates().reset_index(drop=True) 
         
        statIds = [np.unique( ['{:0>2}{:0>3}'.format(data['blockNumber'][i], data['stationNumber'][i] ) for i in range(len(data['blockNumber'] ) )  ]  )[0]  ]
         
        self.consistent_coord = True
        # define a max distance
        
        distances = [ round(self.utils.distance( df.lat[i], df.lon[i] , 
                                        df.lat[0], df.lon[0] , ), 1 ) for i in range(len(df.lat) ) ]         
        
        
        # maximum 20 km distance 
        if max(distances) > 20:
            a = self.check_consistent_coords(lats, lons)
            print('check!')
            
        lats, lons = list(df.lat.values), list(df.lon.values)
        dic = { 'max_date': data['end_date'] , 
                 'min_date': data['start_date'], 
                 'statid': str(statIds), # needed in case multiple ids 
                 'lats': lats, 
                 'lons': lons, 
                 'file': self.file.split('/')[-1] , 
                 'file_path': self.file , 
                 'db': self.dataset,
                 'variables': str(list(np.unique(data['variables'])) ),
                 
                 'consistent_coord': str(self.consistent_coord) }
         
        pd.DataFrame(dic).to_csv( self.out_dir  + '/' + self.file.split('/')[-1] + '.csv', sep = '\t' )        
        self.statid = statIds
        self.lats = lats
        self.lons = lons 
        self.min_date = data['start_date'] 
        self.max_date = data['end_date']

        self.variables = dic['variables']
         

        
        return 0    
    
    def read_ncar(self):
        """ Extract data from a single ODB file """
        file = self.file 
        
        with open(file, 'rt') as infile:
            tmp = infile.read()                                                                                                                                                                                                                    
            data = tmp.splitlines()  
            
        stations_id = []
        dates = []
        lats, lons = [],[]
        
        for i, line in enumerate(data):
            if line[0] == 'H':
                usi      = int(line[2:14])  # unique station identifier
                                     
                ident    = int(line[15:21].replace(' ',''))# WMO
                if ident not in stations_id:
                    stations_id.append(ident)

                year     = line[38:42]  # year
                month    = "%02d" % int(line[43:45])
                day      = "%02d"  % int(line[46:48])
                
                if '99' in day:
                    continue
                        
                dates.append(year + month + day)
                
                lats.append( float(line[57:67]) )
                lons.append( float(line[68:78]) )
                
                stype    = int(line[86:88])
                
                #idate = datetime.strptime(year + month + day + hour + minutes, '%Y%m%d%H%M')
                
        df = pd.DataFrame( {'lat': lats , 'lon' : lons } )
        df = df.drop_duplicates().reset_index(drop=True) 
         
        statIds = [ np.unique(stations_id)[0]  ]
         
        self.consistent_coord = True        
        distances = [ round(self.utils.distance( df.lat[i], df.lon[i] , 
                                        df.lat[0], df.lon[0] , ), 1 ) for i in range(len(df.lat) ) ]         
        
        # maximum 20 km distance allowed
        if max(distances) > 20:
            a = self.check_consistent_coords(lats, lons)
            print('check!')
            
        lats, lons = list(df.lat.values), list(df.lon.values)
        
        dic = { 'max_date': max(dates), 
                 'min_date': min(dates), 
                 'statid': str(statIds), # needed in case multiple ids 
                 'lats': lats, 
                 'lons': lons, 
                 'file': self.file.split('/')[-1] , 
                 'file_path': self.file , 
                 'db': self.dataset,
                 'variables': str( ['38','85','106','107'] ),
                 
                 'consistent_coord': str(self.consistent_coord) }
         
        pd.DataFrame(dic).to_csv( self.out_dir  + '/' + self.file.split('/')[-1] + '.csv', sep = '\t' )        
        self.statid = statIds
        self.lats = lats
        self.lons = lons 
        self.min_date = min(dates)
        self.max_date = max(dates)
        self.variables = dic['variables']
        return 0    
    
    
    def read_igra2(self, file=''):
        """ Extract data from a pandas series extracted from igra2 stations_list """
        series = self.series
        variables = ['85','38','106','107']
        
        dic = { 'max_date': series.end_date , 
                'min_date':series.start_date, 
                'statid':str(series.station_id), # needed in case multiple ids 
                'lats': [series.latitude] , 
                'lons': [series.longitude], 
                'file': self.file , 
                'file_path':  self.file , 
                'db': self.dataset,
                'variables': str(variables),
                'consistent_coord': str('N/A') }  # NB coordinates should be checked against complete files, not stations lists !!! 
        
        # storing the 
        pd.DataFrame(dic).to_csv( self.out_dir  + '/' + self.file.split('/')[-1] + '.csv', sep = '\t' )     
        
        return 0        
    
    
class Inventory():
    """  Class holding the functionalities to read and homogenize 
    the OSCAR, IGRA2, WBAN and CHUAN inventories.
    For OSCAR and WBAN will convert lat and lon to decimal format. """
    
    def __init__(self, datadir ='', tablesdir = '../data/tables/', oscar ="", igra2= "", wban = '', chuan = '' , utils='' ):
        self.datadir = datadir + '/' + datadir
        self.oscar  = datadir + '/' + oscar
        self.igra2  = datadir + '/' + igra2
        self.wban  = datadir + '/' + wban
        self.chuan = datadir + '/' + chuan
        self.utils = utils 
        self.tables_dir = tablesdir
        self.inv = {}
        
    def readCities(self):
        """ Read the cities from the file worldcitiespop.csv to 
        Used to map the station to the closest city """
        df = pd.read_csv(self.tables_dir + '/worldcitiespop.csv')
        self.cities = df
        
    def readInventory(self):
        """ Read the source inventories (OSCAR,IGRA2,WBAN,CHUAN)
        and convert them to pandas dataframe with uniform names for common columns"""
         
        #############################################
        ### Reading IGRA2 data
        #############################################        
        igra2 = pd.read_fwf(self.igra2, widths=(11, 9, 10, 7, 4, 30, 5, 5, 7),
                        names=( 'station_id_igra', 'latitude', 'longitude', 'elevation', 'dummy',
                                      'station_name', 'start_date', 'end_date', 'records') )
        
        # extracting standard 5 digits station id 
        stat_id = [ str(s[6:]) for s in igra2['station_id_igra'] ]
        igra2['station_id'] = stat_id
        igra2 = igra2[['station_id_igra', 'station_name','latitude','longitude', 'station_id', 'start_date', 'end_date']]
        
        # clean for missing values in lat and lon (i.e. values < 0)
        igra2 = igra2.loc[ (igra2['latitude'] >= -180) & (igra2['longitude'] >= -180) ]
        igra2 = igra2.reset_index(drop=True)
        igra2['isRadio'] = 'True'
        igra2.name = 'IGRA2'
        
        #############################################
        ### Reading WBAN data
        #############################################        
        wban =pd.read_fwf(self.wban,widths=(10,5,6,17,22,39,31,31,8,10,10,10,8,7),
                        names=('dum1','station_id','WMO_id','dum0','Country','dum2','station_name','dum3','start_date','end_date',
                                      'latitude','longitude','dum4','Elev'),
                        skiprows=1 )
        
        wban = wban[ [ 'station_id', 'WMO_id' , 'station_name', 'latitude', 'longitude', 'start_date', 'end_date'] ]
        wban = wban.dropna(subset = ['latitude', "longitude"])
        #wban = wban.replace( {'00 00 00': 0 , '0  0  0':0} )
        
        # converting lat, lon to decimal format 
        lat_dec = list(self.utils.degMinSec_to_decimal(wban['latitude'] ))
        lon_dec = list(self.utils.degMinSec_to_decimal(wban['longitude'] ))
        
        wban['latitude'] = lat_dec
        wban['longitude'] = lon_dec

        wban['station_id'] = wban['station_id'].astype('int')
        wban['isRadio'] = None
        wban = wban.reset_index(drop=True)
        wban.name = 'WBAN'
        
        #############################################
        ### Reading OSCAR data
        #############################################
        oscar = pd.read_csv( self.oscar , sep='\t')
        oscar = oscar.rename(columns = {'StationId':'WIGOS', 'Longitude':'longitude', 
                                        'Latitude':'latitude' , 'StationName':'station_name',} )
        
        oscar = oscar [['WIGOS', 'station_name', 'latitude', 'longitude', 'ObsRems']]
        
        # converting lat, lon to decimal format 
        lat_dec = list(self.utils.degMinSec_to_decimal(oscar['latitude'] ))
        lon_dec = list(self.utils.degMinSec_to_decimal(oscar['longitude'] ))
        
        oscar['original_lat'] = oscar['latitude']
        oscar['latitude'] = lat_dec
    
        oscar['original_lon'] = oscar['longitude']
        oscar['longitude'] = lon_dec        
        radioFlag = oscar['ObsRems'].str.contains('Radio').astype(str)
        statids = [f.split('-')[-1] for f in oscar['WIGOS'] ]  # assume that th estation id is the last piece of the WIGOS id 
        oscar['station_id'] = statids
        oscar['start_date'] = 999
        oscar['end_date'] = 999
        oscar['isRadio'] = radioFlag 
        oscar.name = 'OSCAR'
        
        #############################################
        ### Reading CHUAN data
        #############################################        
        with open(self.chuan,'rb') as f:
            rdata=f.read().decode('ISO-8859-1').split('\n')
            z=[[] for _ in rdata[0].split('\t')]
            for r in rdata[1:-1]:
                li=r.split('\t')
#                print(len(li))
                if len(li)>=66:
                    for l in range(66): #len(li)):
                        z[l].append(li[l])

            chuan=dict(zip(rdata[0].split('\t'),z))
            chuan['Lat_DegN']  = [float(l.replace(',','.')) for l in chuan['Lat_DegN'] ]
            chuan['Lon_DegE'] = [float(l.replace(',','.')) for l in chuan['Lon_DegE']]

            chuan=pd.DataFrame(chuan)
            chuan=chuan.rename(columns = {'unique_record_ID':'station_id','WMO#':'WMO_id','Stationname':'station_name',
                                          'Lon_DegE':'longitude', 'Lat_DegN':'latitude', 'Alt_masl':'Elev'})

            # filtering missing values for coordinates
            chuan = chuan.loc[ (chuan['latitude'] > -90) & (chuan['longitude'] > -90)]
            chuan = chuan.reset_index(drop=True)
            chuan['station_id']= chuan['station_id'].astype('int')

            d = dict( zip( [ str(m) for m in range(1,10)] , [ '0'+ str(d) for d in range (1,10)   ] )  )
            chuan = chuan.replace( {'StartStationMonth': d} )
            chuan = chuan.replace( {'EndStationMonth': d} )
            chuan['start_date'] = chuan['StartStationYear'] + chuan['StartStationMonth']
            chuan['end_date'] = chuan['EndStationYear'] + chuan['EndStationMonth']
            
            chuan = chuan[['station_id', 'WMO_id' , 'station_name' , 'latitude' , 'longitude', 'start_date', 'end_date' ]]
            chuan['isRadio'] = None
            chuan.name = 'CHUAN'
        
        # saving csv files in dedicated inventories directory
        if not os.path.isdir('original_inventories'):
            os.mkdir('original_inventories')
        chuan.to_csv('original_inventories/chuan.csv' , sep = '\t')
        igra2.to_csv('original_inventories/igra2.csv' , sep = '\t')
        oscar.to_csv('original_inventories/oscar.csv' , sep = '\t')
        wban.to_csv('original_inventories/wban.csv' , sep = '\t')
               
        # storing df as instances
        self.inv["oscar"] = oscar
        self.inv["igra2"] = igra2
        self.inv["wban"] = wban
        self.inv["chuan"] = chuan
        
        print(0)
        
        

# initialize Utils
utils = Utils()

# Initialize inventory
inventory = Inventory( datadir ='../data/tables',  
                 oscar = "vola_legacy_report.txt" , 
                 igra2 = "igra2-station-list.txt"  , 
                 wban = "WBAN.TXT-2006jan"  , 
                 chuan = "Inventory_ERACLIM_upperair_2.1.txt",
                 utils = utils )

inventory.readInventory()
inventory.readCities()

################
### Analysis part 
################

#f = 'era5.conv._82930'
def wrapper(data, file):
    """ Wrapper function to run all the analysis steps for a single file """
    
    if isinstance(file, pd.Series):
        dataset = 'igra2'
        name_s = file.station_id_igra
    else:
        dataset = file.split('/')[-2]
        name_s = file
        
    print("Doing file: " , name_s )
    try:
        
        d = data.read_file(file=file)
      
        matching_inv = {}
        
        # checking and creating output directory 
        out_dir = 'inventories/' + data.dataset
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
            
        for i in inventories :
            print(' *** Analizing the inventory: ', i )
            analyze = Analyze( data = data, inv = inventory.inv[i], 
                               cities= inventory.cities, utils = utils)
            
            analyze.AddDistance()
            analyze.MatchingId()
            analyze.Closest()
            #matching_inv.append(analyze.found_inv)
            analyze.Combine()
            analyze.FindCity()
            matching_inv[i] = analyze.best_match    
        
        l = [matching_inv[f] for f in  matching_inv.keys() ]
        
        all_inventories = pd.concat(l)
        all_inventories = all_inventories.reset_index()
        
        # Extracting the best WIGOS
        #valid = np.where(all_inventories['latitude'].isfinite() )
        #df_red = all_inventories.drop_na(by=['latitude'])
        
        # removing non-matched files and inventories 
        df_red = all_inventories.dropna(subset=['latitude'])
        
        if df_red.empty:
            print("+++ No matching inventory for the file")
            out = open('inventories/' + dataset + '_unidentified.txt', 'a+')
            out.write(file + '\n')
            return 
        
        else:
            #if 'OSCAR' in list(df_red['inventory']):
            #    df_red = df_red.loc[df_red['inventory'] == 'OSCAR']
            best_wigos = df_red.iloc[0]['WIGOS_calc']
            if not data.consistent_coord:  # coordinates are inconsistent, skipping 
                best_wigos  = '0-20999-0-' + str(eval(data.statid[0])[0])
                a = open('inventories/' + dataset + '_inconsistent_coordinates.txt', 'a+')
                a.write(file + '\n')
                
            
        
        # Add file information     
        all_inventories['lat_file'] = str( [max(data.lats), min(data.lats)] )
        all_inventories['lon_file'] = str( [max(data.lons), min(data.lons)]  )
        all_inventories['file'] = data.file.split('/')[-1]
        all_inventories['file_dataset'] = data.dataset
        all_inventories['file_min_date'] = str(data.min_date)
        all_inventories['file_max_date'] = str(data.max_date)
        all_inventories['file_statid'] = data.statid[0]
        
        all_inventories['WIGOS_best'] = best_wigos # TODO to fix 
        all_inventories['variables'] = str(data.variables)  
        
        name = out_dir + '/' + data.file.split('/')[-1] + '_inventories.csv'
        
        all_inventories.to_csv( name, sep = '\t', index = False)
        
        all_inventories_red = all_inventories[ ['file_statid', 'station_id', 'station_name', 'latitude', 'longitude', 
                                                'original_lat', 'original_lon', 'distance_km', 'lat_file',
                                                'lon_file',  'inventory', 'WIGOS', 'WMO_id', 'WIGOS_calc', 'file_min_date', 
                                                'file_max_date', 'start_date', 'end_date', 'isRadio', 'WIGOS_best', 'city', 'variables'] ]
        
        all_inventories_red.to_csv( name.replace('.csv','_reduced.csv'), sep = '\t' , index = False )
        
        print("Done :::" , data.file )
        a = open('inventories/' + dataset + '_processed_all.txt', 'a+')
        a.write(name_s + '\n')       
        
    
    except:
        print("*** Cannot read file! ***")
        a = open('inventories/' + dataset +"_failed_files.txt", "a+")
        a.write(name_s + '\n')       
    
    
    

datasets = { 'era5_1': '/mnt/users/scratch/leo/scratch/era5/odbs/1' ,
                               'era5_2': '/mnt/users/scratch/leo/scratch/era5/odbs/2',
                               'era5_3188': '/mnt/users/scratch/leo/scratch/era5/odbs/3188',
                               'era5_1759': '/mnt/users/scratch/leo/scratch/era5/odbs/1759',
                               'era5_1761': '/mnt/users/scratch/leo/scratch/era5/odbs/1761',
                               
                               'bufr': '/mnt/users/scratch/leo/scratch/era5/odbs/ai_bfr/',                                   
                               'ncar': '/scratch/das/federico/databases/UADB',
                               'igra2': '', # dummy, use the igra2 station list file 

                               } 
    


if __name__ == '__main__':
    """ Parameters:
          - POOL: runs multiprocesses (default=30)
          - CHECK_MISSING: only runs missing files, otherwise will rerun and replace existing files """

    databases = ['era5_2']
        
    databases = [ 'era5_2', 'era5_1759', 'era5_1761', 'era5_3188', 'bufr', 'ncar', 'igra2', 'era5_1']
    
    #databases = [ 'era5_3188', 'bufr', 'ncar', 'igra2', 'era5_1']

    # list of inventories to consult  
    inventories  = [ 'oscar', 'igra2', 'wban', 'chuan']
    alldb = [ 'era5_2', 'era5_1759', 'era5_1761', 'era5_3188', 'bufr', 'ncar', 'igra2', 'era5_1']
    
    era5_block = [ 'era5_2', 'era5_1759', 'era5_1761', 'era5_3188']
    
    databases = [ 'igra2', 'ncar', 'bufr']
    
    databases = era5_block

    databases = alldb
    
    # enable multiprocesing
    POOL = True
    n_pool = 40
    # only process missing files 
    CHECK_MISSING = False  
    CHECK_FAILED = False
    # loop through each of the databases
    for db in databases:
        
        # getting only missing files, option CHECK_MISSING 
        if db == 'era5_1':
            flist=glob.glob(datasets[db] + "/era5.conv._*")
            #flist =[f for f in flist if '11035' in f ]
        elif db == 'era5_2':
            flist=glob.glob("/mnt/users/scratch/leo/scratch/era5/odbs/2/era5.conv._*")
            flist=[f for f in flist if '.gz' not in f and '.nc' not in f ]
            
        elif db == 'era5_1759':
            flist=glob.glob("/mnt/users/scratch/leo/scratch/era5/odbs/1759/era5.1759.conv.*")
            
        elif db == 'era5_1761':
            flist=glob.glob("/mnt/users/scratch/leo/scratch/era5/odbs/1761/era5.1761.conv.*")
            
        elif db == 'era5_3188':
            flist=glob.glob("/mnt/users/scratch/leo/scratch/era5/odbs/3188/era5.3188.conv.*")
            
        elif db == 'ncar':
            flist=glob.glob( datasets[db] + '/uadb*')
            
        elif db== 'igra2':
            f = inventory.inv['igra2']
            flist = [ f.iloc[n] for n in range(len(f)) ]
            
            # TODO special filterfor IGRA stations 
            #flist = [f for f in flist if f.station_id_igra == 'BLXUAE04838']
            #print(0)

        elif 'bufr' in db:
            flist=glob.glob(datasets[db] + '/'+'era5.*.bfr')
            print(0)
        # general cleaning of file lists    
        
        flist = [f for f in flist if '.gz' not in f and '.nc' not in f ]
        
        if CHECK_FAILED:
            failed = open(db.replace('era5_','') + '_failed_files.txt','r').readlines()
            flist = a
        #flist = [f  for f in flist if  'era5.1759.conv.2:81605' in f ]
        if db != 'igra2':   
            if CHECK_MISSING:
                flist_c = [f.split('/')[-1] for f in flist ]
                try: 
                    processed = [f.replace("_inventories.csv","") for f in os.listdir('inventories/'+db) ]
                except:
                    processed = []
                missing_stat = [f for f in flist_c if f not in processed ]
                flist = [datasets[db] + '/' + f for f in missing_stat ]
                
            # sort files from smallest 
            pairs = []
            for f in flist:
                # Get size and add to list of tuples.
                size = os.path.getsize(f)
                pairs.append((size, f))
        
            # Sort list of tuples by the first element, size.
            pairs.sort(key=lambda s: s[0])
            flist = [ s[1] for s in pairs ]
        
            
        data = Data(dataset=db, utils = utils )
        
        
        #failed = open('inventories/' + db.replace('era5_','').replace('ncar','UADB') + '_failed_files.txt','r').readlines()
        #failed = [f.replace('\n','') for f in failed]
        #flist = failed

        #unidentified = open('inventories/' + db.replace('era5_','') + '_unidentified.txt','r').readlines()
        #unidentified = [f.replace('\n','') for f in unidentified]
        
        #flist = [f for f in flist if f not in unidentified ]
        # TODO edit here to possibly filter file list
    
        if POOL: # running multiprocessing or single 
            p = Pool(n_pool)
            func = partial(wrapper,data)
            out = p.map(func, flist)   
                
        else:
            for f in flist:
                w = wrapper(data, f)


                
                
# example: bad coordinates '/mnt/users/scratch/leo/scratch/era5/odbs/1759/era5.1759.conv.1:57606'
 # example bad coordinates extreme: /mnt/users/scratch/leo/scratch/era5/odbs/1759/era5.1759.conv.1920 
 
 # TODO fix min_date era5.conv._84013 
 
"""
/mnt/users/scratch/leo/scratch/era5/odbs/1759/era5.1759.conv.6:100954
/mnt/users/scratch/leo/scratch/era5/odbs/1759/era5.1759.conv.6:102317
/mnt/users/scratch/leo/scratch/era5/odbs/1759/era5.1759.conv.1:67194q
"""
