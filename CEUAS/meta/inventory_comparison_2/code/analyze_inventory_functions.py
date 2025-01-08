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
from datetime import date, datetime,timedelta

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# pd.set_option('display.max_columns', None)
# pd.set_option('display.max_rows', None)
# pd.set_option('display.width', None)
# pd.set_option('display.max_colwidth', -1)

from multiprocessing import Pool
from functools  import partial

from eccodes import *

from tqdm import tqdm
"""
Definitions
- consistent_coordinates = 30 km (do not change over time)
- distance wrt station-inventory = 30 km 
- df_red = df_red.loc[df_red["distance_km"] <= 30 ]

"""
int_void = -2147483648


def check_read_file(file='', read= False):
    with open(file, 'rt') as infile:
        tmp = infile.read()  # alternative readlines (slower)                                                                                                                                                                                                                   
        data = tmp.splitlines()  # Memory (faster)  

    return data

def igra2_ascii_to_dataframe(file=''):
    """ Read an igra2 stationfile in ASCII format and convert to a Pandas DataFrame. 
        Adapted from https://github.com/MBlaschek/CEUAS/tree/master/CEUAS/data/igra/read.py 
        Variables used inside the DataFrame are already CDM compliant
        df.types

        Args:
             file (str): path to the igra2 station file

        Returns:
             Pandas DataFrame with cdm compliant column names
    """

    data = check_read_file(file=file, read=True)
    #source_file = [l for l in file.split('/') if '.txt' in l][0]
    read_data = [] #  Lists containing the raw data from the ascii file, and the observation dates
    """ Data to be extracted and stored from the igra2 station files 
        Some info is contained in the header of each ascent, some in the following data """

    """ Initialize the variables that can be read from the igra2 files """
    ident,year,month,day,hour,reltime,p_src,np_src,lat, lon = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan 
    lvltyp1,lvltyp2,etime,press,pflag,gph,zflag,temp,tflag,rh,dpdep,wdir,wspd = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan # initialize to zeros
    stations_id = []
    idate = np.nan
    count = 0
    head_count = 0

    obs_id = 0

    def make_release_time(date_time,  hour, release):
        """ build a sonde release time 
              ex 2019 02 20 00 2349 
              ex 2019 01 10 00 0011 
              They round the release time to the closest hour. 
              It can be the same day or the following !!!
              date_time = date_time python object, 
              date, time, release = original strings 
              """
        release_h  = int(release[:2])
        release_m = int(release[2:4])

        if release_h == 99:
            release_date_time = date_time 
            timestamp_flag = int_void #largest integer number int 64   

        else:
            timestamp_flag = 1 

            if release_m == 99:
                release_m = 0
            release_date_time = date_time.replace(hour= release_h, minute= release_m) 

            """ Here, I have to subtract one day to the release time stamp if the hour of the time stamp is in the night,
                  but the nominal time is reported at midnight hence in the following day. For example  2019 02 20 00 2349 from file VMM00048820 """
            if hour == '00':
                if release_h > 20:
                    release_date_time = release_date_time - timedelta(days=1)
                else:
                    pass

        return release_date_time, timestamp_flag 

    #data = data[3706329:]  # TO DO HERE TO DO CHANGE!!!
    #for i, line in enumerate(data):
    #      print(i , '  ' ,  line[0] )

    hash_lines = [line for i, line in enumerate(data) if line.strip().startswith("#")]
    for line in [hash_lines[0], hash_lines[-1]] :
        if line[0] == '#':
            head_count = head_count +1 
            # Info from the Header line of each ascent                                                                                                                                                                                                                   
            ident     = line[1:12]               # station identifier
            #ident     = ident[6:12]
            if ident not in stations_id:
                stations_id.append(ident)

            year      = line[13:17]               # year, months, day, hour of the observation
            month   = line[18:20]
            day       = line[21:23]
            hour      = line[24:26]               
            reltime  = line[27:31]            # release time of the sounding.
            numlev  = int(line[32:36])        # number of levels in the sounding == number of data recorded in the ascent
            p_src     = line[37:45]              # data source code for the pressure levels 
            np_src   = line[46:54]             # data source code for non-pressure levels
            lat         = int(line[55:62]) / 10000.  # latitude and longitude
            lon        = int(line[63:71]) / 10000.
            #observation_id = i
            if int(hour) == 99:
                time = reltime + '00'
            else:
                time = hour + '0000'

            if '99' in time:
                time = time.replace('99', '00')

            idate = datetime.strptime(year + month + day + time, '%Y%m%d%H%M%S') # constructed according to CDM

            ### making release time and its flag
            release_time , report_timeflag = make_release_time(idate, hour, reltime) # making the release time 

            ### report_meaning_of_timestamp	int	meaning_of_time_stamp:meaning	Report time - beginning, middle or end of reporting period
            ### 1	beginning	Date / time specified indicates the start of the period over which the observation was made.
            ### end	Date / time specified indicates the end of the period over which the observation was made.
            ### middle	Date / time specified indicates the middle of the period over which the observation was made.

            iday =  int(year + month + day)
            count = count + 1

        read_data.append ( ( 'IGRA2'.rjust(10), head_count,  idate, iday, ident, lat, lon, release_time, report_timeflag) )

    column_names_igra2 = [ 'source_id', 'report_id',  'record_timestamp' , 'iday', 'statid@hdr', 'lat@hdr', 'lon@hdr', 'report_timestamp', 'report_meaning_of_timestamp',]

    df = pd.DataFrame(data= read_data, columns= column_names_igra2)
    return df



class Utils():
    """ Class holding a series of calculation utilites """
    def __init__(self, treshold = 30):
        """ threshold: maximum allowed distance in Km for points to be considered compatible (does not change more than threshold over time ) """
        self.treshold = treshold
        
        # station_configuration columns and mapping
        """
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
        """
        
    def distance(self, lat1,lon1,lat2,lon2):
        """ Return the Harvesine distance between two points """
        d = geopy.distance.geodesic( (lat1,lon1), (lat2,lon2) ).km
        return d     

    def oscar_coord(self, oscar):
        """ Converts OSCAR coordinates to standard format, since it is not given in the range [-90,90] and [-180,180]"""
        vlat=[]
        vlon=[]
        
        for v in oscar['latitude'].values:
            dms=v[:-1].split()
            if v[-1]=='S':
                vlat.append(-(float(dms[0])+float(dms[1])/60+float(dms[2])/3600))
            else:
                vlat.append(float(dms[0])+float(dms[1])/60+float(dms[2])/3600)
        for v in oscar['longitude'].values:
            dms=v[:-1].split()
            if v[-1]=='W':
                vlon.append(-(float(dms[0])+float(dms[1])/60+float(dms[2])/3600))
            else:
                vlon.append(float(dms[0])+float(dms[1])/60+float(dms[2])/3600)
                
        oscar['latitude']=np.array(vlat)
        oscar['longitude']=np.array(vlon)        
            
        return oscar
            
    def degMinSec_to_decimal(self, coord):
        """ Converts lat and lon [lists] from degrees-minute-seconds to decimal """
                
                
        def dms2dec(dms):
            """ Split strings like -75 14 00N into degree, minute and seconds,
            convert into decimal with properly assigned sign """
            
            sign = -1 if re.search('[-swSW]', dms) else 1
            
            dms = re.sub ( '[wWsSnNeE]', '', dms  )            
            dms = [ d for d in dms.split(' ') if d ]

            d = dms[0]
            m = dms[1]
            s = dms[2]

            dec = abs(float(d)) + float(m)/60 + float(s)/3600
            return sign*round(dec, 2)
        
        coord_dec = map(dms2dec, coord)        

                
        return coord_dec
    
    
    #def coord_consistency(self, lats='',lons='' ):
    #    """ Check if the absolute variation in time of the coordinates is below a certain threshold in Km """
            
    #    #check = self.distance( min(lats) , min(lons), max(lats), max(lons)  ) < threshold 
    #    check = self.distance( min(lats) , min(lons), max(lats), max(lons)  ) < self.threshold 
        
    #    return check
    
    
class Analyze():
    
    def __init__(self, dist_limit = 30, inv = '', data = '', utils='', cities = ''):
        """ parameters ::
                     dist_limit : maximum distance allowed for station matching, default = 30km
                     inv :: dataframe of a single the inventory
                     data :: data extracted from a single file from one of the data sets """
        
        self.dist_limit = dist_limit
        self.inv = inv
        self.data = data
        self.found_inv= { self.inv.name: {} }
        self.utils = utils 
        
        self.mobile_datasets = ['npsound' , 'era5_1_mobile' , 'era5_2_mobile' , 'shipsound']
        
        
        self.wigos = { 'IGRA2'    :'0-20300-0-',
                       'CHUAN'    :'0-20400-0-',                       
                       'WBAN'     :'0-20500-0-',
                       
	               'SCHROEDER':'0-20600-0-',
                       
	               'WMO'      :'0-20700-0-',
                       
                       'AMMA' : '0-20800-0-',
                       'GIUB' : '0-20900-0-' ,
                       
                       'HARA' : '0-20111-0-' ,
        }
        
        self.best_match = ''
        self.cities = cities 
        
    def Make_WIGOS(self, df):
        """ Creates a pseudo WIGOS id for the file given the best matching inventory data.
              For OSCAR stations, we use the real WIGOS from the inventory.
              For th eother inventories, we create a pseudo WIGOS using the convention: 0-20x00-0-,
              see self.wigos in the initialization of the class """
        
        #df = self.matchingInv_dist_id
        # case OSCAR ds
        
        """
        if self.inv.name == 'OSCAR':
            
            if len(df) == 1: # only one matching, keep this 
                wigos = df['WIGOS'][0]
            else:
                # TO DO handle special case?        
                for i,text in enumerate(df['ObsRems'][:]) :
                    if ('Radiosonde' in text or 'Upper' in text):
                        wigos = df['WIGOS'][i]
                        break
                wigos = df['WIGOS'][0]
                
        else:
            prefix = self.wigos[self.inv.name]
            if len(df) == 1: # only one matching, keep this 
                wigos = prefix + str(df['station_id'][0])
            else:
                wigos = [prefix + str(df['station_id'].values[i])  for i in range(0, len(df)) ]
                 
        return wigos
        """
        
        if self.inv.name == 'OSCAR':
            wigos = df.WIGOS.values
            wigos = list(wigos)
            
        else:
            prefix = self.wigos[self.inv.name]
            if len(df) == 1: # only one matching, keep this 
                wigos = prefix + str(df['station_id'][0])
            else:
                wigos = [prefix + str(df['station_id'].values[i])  for i in range(0, len(df)) ]
                 
        return wigos

    def FindCity(self):
        """ Find the closest city within a 100 km radius """
        
        print('*** Find CIty *** ')
        if len(self.data.lats) < 1: # edit form "if not self.data.lats:""
            self.best_match['city'] = 'None'
            self.best_match['city_dist_km'] = 'None'
            self.best_match['city_lat'] = 'None'
            self.best_match['city_lon'] = 'None'               
            
        else:
            lat, lon = self.data.lats[0] , self.data.lons[0]
            cities = self.cities
            cities = cities.loc [  (abs(cities['Latitude'] - lat) < 5 ) & (abs(cities['Longitude'] - lon) < 5) ]  # skimming df 
            
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
        
        
        df['latitude'] = df['latitude'].astype(float)
        df['longitude'] = df['longitude'].astype(float)
        
        if len(self.data.lats) ==0 or 'mobile' in self.data.dataset: # empty dummy inventory for mobile stations 
            red_inv = pd.DataFrame(columns = df.columns )
            red_inv['distance_km'] = 'NAN'
            self.red_inv = red_inv 
            red_inv.name = self.inv.name                
            distances = []
            return distances
        
        
        red_inv = df.loc [  (abs(df['latitude']- self.data.lats[-1]) < 2 ) & (abs(df['longitude'] - self.data.lons[-1]) < 2) ] 
        
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
            
            return 
            
    def MatchingId(self):
        """ Check if there is a station id matching the input station id data.
             Note that statId is a list, most frequently with one element only,
             but might contain more if more ids are found in a single file e.g. for NCAR dataset.
             """
        
        print("*** Find MatchingId ***")
        
        if self.inv.name in ['CHUAN', 'WBAN', 'SCHROEDER', 'WMO', 'AMMA', 'GIUB', "HARA"]:
            stats = 'WMO_id'
     
        else:
            stats = 'station_id'
            
        if isinstance(self.data.statid, pd.Series):
            try:
                s = eval(self.data.statid[0])
            except:
                s = [self.data.statid[0]]
            stat_ids = s 
                    
        elif isinstance(self.data.statid, list):            
            s = [str(s) for s in self.data.statid ]
            stat_ids = [ f.split(':')[1]  if ':' in f else f for f in s ]
                
        if not isinstance(stat_ids, list):
            stat_ids = [stat_ids ]
            
        df = self.inv.loc[ self.inv[stats].isin( stat_ids ) ] 
        
        # test the missing latitude sign for ERA5 1759 inventory 
        # The HP is that the station id must match the id in the WBAN inventory, AND by adding a minus sing to the latitude the distance becomes small
        # Here I select the WBAN entry with the same station id only if all the other stations do not have compatible distance calculated
        # using original lat and lon 
        if len(df) ==0 and self.inv.name=='WBAN' and self.data.dataset=="era5_1759":
            df_alt = self.inv.loc[ self.inv['station_id'].isin( stat_ids ) ]
            if len(df_alt) >0:      
                df = df_alt
            else:
                print(0)

        inv_matchingId = df
        inv_matchingId.name = self.inv.name 
        inv_matchingId = inv_matchingId.reset_index(drop=True)        
        
        lats = self.data.lats
        lons = self.data.lons
        
        if len(lats) > 0: # edit from "if lats:"
            distances = [ round( self.utils.distance( lats[-1], lons[-1] , inv_matchingId['latitude'][i] , inv_matchingId['longitude'][i] ), 1 ) for i in range(len(inv_matchingId['latitude']) ) ]
    
            inv_matchingId['distance_km'] = distances
            inv_matchingId['distance_km_minusLat'] = ['' for i in range(len(inv_matchingId))]  # must add to keep dataframes compatible in all cases
        
        else: # edit from "else:"
            
            inv_matchingId['distance_km'] = ['' for i in range(len(inv_matchingId))]
            
        
        # if the lat in the WBAN inventory is positive, nothing to be done
        if len(df) >0 and self.inv.name=='WBAN' and self.data.dataset=="era5_1759":
            
            if inv_matchingId.latitude.values[0] >0:
                pass
            else: # try to add minus sign in latitude 
                #lat_pos_ind = np.where( np.array(lats) >0 )[0][0]
                #lat_pos = lats[lat_pos_ind]
                #lon_pos = lons[lat_pos_ind]
                
                #distances_swap = [ round( self.utils.distance( -lat_pos , lon_pos, inv_matchingId['latitude'][i] , inv_matchingId['longitude'][i] ), 1 ) for i in range(len(inv_matchingId['latitude']) ) ]   
                #inv_matchingId['distance_km_minusLat'] = distances_swap
                pass
            
        self.found_inv[self.inv.name ]['matching_id'] = inv_matchingId
        
               
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

        a=0
        
        
        
        
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
            if self.inv.name != 'wban':
                res = res.drop_duplicates(subset= [c for c in res.keys() if c != 'distance_km_minusLat' ]  )
            else:
                res = res.drop_duplicates()
            
        elif ( id_match.empty and not dist_match.empty ):
            res = dist_match
        elif ( not id_match.empty and dist_match.empty ):
            res = id_match
            
        else:
            res = id_match # empty case
            for c in res.columns:
                res[c] = [np.nan]            
            #if 'mobile' not in self.data.dataset and 'npsound' not in self.data.dataset and 'shipsound' not in self.data.dataset: # TODO simplify ?
            if self.data.dataset not in self.mobile_datasets:
                res['WIGOS_calc'] = ['None'] 
            else:
                res['WIGOS_calc'] = ['None'] 
                
                
            res['inventory'] = self.inv.name 
            self.best_match = res 
        
            #return
        res = res.reset_index()
        wigos = self.Make_WIGOS(res)
        res['inventory'] = self.inv.name         
        res['WIGOS_calc'] = wigos 
        
        self.best_match = res 
                


class Data():
    """ Main class to hold utilities for reading original input data from 
         - ERA5 1,2,1759,1761,3188 
         - IGRA2
         - NCAR
         - BUFR 
         """
    
    def __init__(self, dataset='', utils = '', inv=''):
        
        self.dataset = dataset
        self.utils = utils 
        self.inventory=inv
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
        
        # TO DO need to check this if it is still the same
        self.test_era5_1759_lat_mismatch_all = ['era5.1759.conv.2:82606', 'era5.1759.conv.2:82607', 'era5.1759.conv.2:80430', 'era5.1759.conv.2:82608', 'era5.1759.conv.2:80431', 'era5.1759.conv.2:80433', 
                                           'era5.1759.conv.2:80505', 'era5.1759.conv.2:80506', 'era5.1759.conv.2:80509', 'era5.1759.conv.2:80510', 'era5.1759.conv.2:80704', 
                                           'era5.1759.conv.2:81302', 'era5.1759.conv.2:81404', 'era5.1759.conv.2:81408', 'era5.1759.conv.2:81410', 'era5.1759.conv.2:81502', 
                                           'era5.1759.conv.2:81605', 'era5.1759.conv.2:50302', 'era5.1759.conv.2:51301', 'era5.1759.conv.2:51302', 'era5.1759.conv.2:51304', 
                                           'era5.1759.conv.2:50801', 'era5.1759.conv.2:52401', 'era5.1759.conv.2:50802', 'era5.1759.conv.2:60501', 'era5.1759.conv.2:60701', 
                                           'era5.1759.conv.2:61503', 'era5.1759.conv.2:60702', 'era5.1759.conv.2:61703', 'era5.1759.conv.2:61502', 'era5.1759.conv.2:62702', 
                                           'era5.1759.conv.2:50304', 'era5.1759.conv.2:61701', 'era5.1759.conv.2:50305', 'era5.1759.conv.2:80402', 'era5.1759.conv.2:62701', 
                                           'era5.1759.conv.2:50306', 'era5.1759.conv.2:80604', 'era5.1759.conv.2:80405', 'era5.1759.conv.2:81306', 'era5.1759.conv.2:50307', 
                                           'era5.1759.conv.2:80406', 'era5.1759.conv.2:81403', 'era5.1759.conv.2:50309', 'era5.1759.conv.2:80408', 'era5.1759.conv.2:50310', 
                                           'era5.1759.conv.2:81405', 'era5.1759.conv.2:80409', 'era5.1759.conv.2:50402', 'era5.1759.conv.2:81406', 'era5.1759.conv.2:80412', 
                                           'era5.1759.conv.2:81407', 'era5.1759.conv.2:50403', 'era5.1759.conv.2:81501', 'era5.1759.conv.2:50405', 'era5.1759.conv.2:80419', 
                                           'era5.1759.conv.2:81608', 'era5.1759.conv.2:50406', 'era5.1759.conv.2:80502', 'era5.1759.conv.2:81702', 'era5.1759.conv.2:51303', 
                                           'era5.1759.conv.2:80504', 'era5.1759.conv.2:61501', 'era5.1759.conv.2:80511', 'era5.1759.conv.2:82401', 'era5.1759.conv.2:80301', 
                                           'era5.1759.conv.2:80512', 'era5.1759.conv.2:82404', 'era5.1759.conv.2:80302', 'era5.1759.conv.2:80513', 'era5.1759.conv.2:82408', 
                                           'era5.1759.conv.2:83405', 'era5.1759.conv.2:80303', 'era5.1759.conv.2:83504', 'era5.1759.conv.2:80304', 'era5.1759.conv.2:80602', 
                                           'era5.1759.conv.2:83508', 'era5.1759.conv.2:80306', 'era5.1759.conv.2:80310', 'era5.1759.conv.2:80702', 'era5.1759.conv.2:80311', 
                                           'era5.1759.conv.2:81301', 'era5.1759.conv.2:80413', 'era5.1759.conv.2:81304', 'era5.1759.conv.2:80415', 'era5.1759.conv.2:81409', 
                                           'era5.1759.conv.2:80421', 'era5.1759.conv.2:81601', 'era5.1759.conv.2:80422', 'era5.1759.conv.2:81603', 'era5.1759.conv.2:80424', 
                                           'era5.1759.conv.2:81606', 'era5.1759.conv.2:80425', 'era5.1759.conv.2:81609', 'era5.1759.conv.2:80426', 'era5.1759.conv.2:80427', 
                                           'era5.1759.conv.2:82403', 'era5.1759.conv.2:80429', 'era5.1759.conv.2:82405', 'era5.1759.conv.2:82501', 'era5.1759.conv.2:82503', 
                                           'era5.1759.conv.2:82507', 'era5.1759.conv.2:80432', 'era5.1759.conv.2:50404', 'era5.1759.conv.2:80305', 'era5.1759.conv.2:80309', 
                                           'era5.1759.conv.2:80705', 'era5.1759.conv.2:81604']
        
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
        
        # to be filled afterwards
        # self.file = ''
        # self.dataset = ''
        # self.min_date = ''
        # self.max_date = ''
        # self.consistent_coord = ''
        # self.variables = ''
        #self.lats = ''
        #self.lons = ''
        
    def read_file(self, file):
        """ Wrapper for the functions to read input file from the datasets. 
             If a csv file exists holding the needed information, it is read.
             Otherwise it reads the original file with the appropriate functions. """
        
        # for igra2, file is indeed a pandas series 
        # if self.dataset == 'igra2':
        #     self.file = file.station_id_igra
        #     self.series = file 
            
        # else:
        ## IGRA2CHANGE
        self.file = file  

        f = 'temp_data/' + self.dataset + '/' + self.file.split('/')[-1] + '.csv'
        
        # NB DOES NOT WORK WELL WITH HARA, REDO TEMP FILES !!! 
        
        if os.path.isfile( f ):         # read file if csv already exists in temp_data directory 
            df = pd.read_csv(f,  sep='\t', converters={'statid': str}  ).astype('str')
            
            stat_id =  df['statid'].astype('str')
            
            if self.dataset == 'igra2':
                stat_id = [df['statid'].astype('str').values[0]]
            else:
                try:
                    stat_id = eval(df.statid[0])[0]
                except:
                    stat_id = df.statid[0]
                    
            if not isinstance(stat_id, list):
                try:
                    stat_id = [stat_id]
                except:
                    a=0
                    
            self.statid = stat_id
            
            """
            # fixing weird lats, lons values with point without digit 
            Lats, Lons =[],[]
            for lat,lon in zip( df['lats'].values, df['lons'].values ) :
                try:
                    l = float(lat)
                except:
                    lat = float(str(lat).replace('.','.0') )
                try:
                    a = float(lo)
                except:
                    lon = float(str(lon).replace('.','.0') )                    
                Lats.append(lat)
                Lons.append(lon)
            df['lats'] = Lats
            df['lons'] = Lons 
            """
            
            if 'mobile' in self.dataset:
                self.lats = []
                self.lons = []
                self.lats_all=[]
                self.lons_all = []
            else:
                try:
                    
                    if type(eval(df['lats'].values[0])) == list:
                        self.lats = eval(df['lats'].values[0]) 
                        self.lons = eval(df['lons'].values[0]) 
                    else:
                        self.lats = [eval(df['lats'].values[0]) ]   
                        self.lons = [eval(df['lons'].values[0]) ]
                except:
                    self.lats = []
                    self.lons = []
                    
                if 'lats_all' in df.columns:
                    if type(eval(df['lats_all'].values[0])) == list:
                            self.lats_all = eval(df['lats_all'].values[0]) 
                            self.lons_all = eval(df['lons_all'].values[0]) 
                    else:
                            self.lats_all = [eval(df['lats_all'].values[0]) ]   
                            self.lons_all = [eval(df['lons_all'].values[0]) ]
                        
                    
            if self.file.split("/")[-1] in self.test_era5_1759_lat_mismatch_all:
                self.lats = [-l for l in self.lats if l >0 ]
                
            self.min_date = df['min_date'].values[0]
            self.max_date = df['max_date'].values[0]
            try:
                self.consistent_coord = eval( df['consistent_coord'].values[0] )
            except:
                if self.dataset == 'igra2':
                    a = 0
                    
            self.variables = df['variables'].values[0]
            
            if 'frequency' in df.columns:
                self.frequency = df['frequency'].values[0]
            else:
                self.frequency = '[unknown]'
        
        else: # read from original files
            if 'era5_' in self.dataset and 'mobile' not in self.dataset :
                self.read_odb()
            elif self.dataset == 'era5_1_mobile':
                self.read_era5_1_mobile()
            elif self.dataset == 'era5_2_mobile':
                self.read_era5_1_mobile()                
            elif 'ncar' in self.dataset:
                self.read_ncar()
            elif 'bufr_cnr' in self.dataset :
                self.read_bufr_cnr()
            elif 'bufr' in self.dataset :
                self.read_bufr()
            elif 'maestro' in self.dataset :
                self.read_bufr_maestro()
            elif 'amma' in self.dataset:
                self.read_amma()
            elif 'giub' in self.dataset:
                self.read_giub()
            elif 'igra2' in self.dataset:
                self.read_igra2()
            elif 'hara' in self.dataset:
                self.read_hara()
            elif self.dataset in ['npsound' , 'shipsound']:
                self.read_shipsound_npsound()
   
                
    def read_era5_1_mobile(self):
        """ Works also for ERA5 2, might change th ename of the function """
        # Must read header from non mobile directory...
               
        file = self.file 
        statIds = file.split('.txt.gz')[0].split('._')[1]
        
        # need full list to extract max date 
        if self.dataset == 'era5_1_mobile':
            flist_all= [ f for f in glob.glob(datasets[self.dataset] + "/era5.conv.*._*.txt.gz" ) if f.split('._')[1].split('.txt.')[0] == statIds]  
            flist_all.sort()
            max_date = flist_all[-1].split('conv.')[1].split('._')[0]
            min_date = self.file.split('conv.')[1].split('._')[0]
            
        elif self.dataset == 'era5_2_mobile':
            flist_all= [ f for f in glob.glob(datasets[self.dataset] + "/era5.conv.*.*" ) if f.split('.gz')[0].split('.')[-1] == statIds]  
            flist_all.sort()
            max_date = flist_all[-1].split('conv.')[1].split('.')[0]
            min_date = flist_all[0].split('conv.')[1].split('.')[0]

        file = 'era5.conv._' + statIds 
        dic = { 'max_date': max_date, 
                 'min_date': min_date, 
                 'statid': str(statIds), # needed in case multiple ids 
                 'lats': '', 
                 'lons': '', 

                 'lats_all': '', 
                 'lons_all': '', 
                 
                 'file': file , 
                 'file_path': self.file , 
                 'db': self.dataset,
                 'variables': '',
                 'frequency': '',
                 'consistent_coord': 'False' }
        
        self.statid = [statIds]
        
        self.lats = []
        self.lons = [] 
        self.min_date = min_date
        self.max_date = max_date
        self.variables = dic['variables']
        self.consistent_coord = False
        self.file = file # update reduced file name 
        
        pd.DataFrame(dic, index=[0]).to_csv( 'temp_data/'  + self.dataset +  '/' + self.file.split('/')[-1] + '.csv', sep = '\t')

        
        
    def read_odb(self):
        """ Extract data from a single ODB file 
        # odb headre -i file """
        
        file = self.file 

        #a = open('check_odb','a+')
        #a.write(file + '\n')
        #a.close()

        """ In the first iteration, extract a list of lat and lon values to see if they change in time due to relocation.
        Then we calculate the distance between the first location and all the others to see if they are within a distance of 30 km maximum.
        If there are far away points, we check if the fractions of data far apart is small. We do so by grouping by distance. 
        """
        #1. Extract unique lat and lon        
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
        del b
        
        self.consistent_coord = True
        # define a max distance
        
        #2. Check if locations are sufficiently close, by calculating the distance of the first location with the rest 
        # I do this preliminary check for  ODB files sinc e this is quite slow
        
        distances = [ round(self.utils.distance( lats[i], lons[i] , 
                                        lats[0], lons[0] , ), 1 ) for i in range(len(lats) ) ]                 
        if max(distances) > 30:
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
            del b
            ### TODO -> implement the funtion that does this 
            # group by equal lat/lon values 
            df = pd.DataFrame( {'lat': lats_all , 'lon' : lons_all } )
            grouped = df.groupby(df.columns.tolist(), as_index=False).size()
            grouped = grouped.sort_values(by=['size'], ascending = False)
            # add distance column
            distances = [ round(self.utils.distance( grouped['lat'].values[i], grouped['lon'].values[i] , 
                                             grouped['lat'].values[0], grouped['lon'].values[0] ), 1 ) for i in range(len(grouped['lat'].values) ) ] 
    
            grouped['distance_km'] = distances
            
            frequency = str(list(grouped['size'].values)) # check if works
            self.frequency = frequency
            #TODO I implemented a separate function that does this, see later 
            """
            grouped
                lat     lon     size  distance_km
            0  28.05  9.63    20       0.0        
            1  31.10  121.25  4266     10252.9    
            2  31.17  121.43  2156598  10262.5
            """
            compatible_coord = grouped.loc[grouped['distance_km'] <= self.utils.treshold ]
            
            frac = sum(compatible_coord['size'].values) / sum(grouped['size'].values )

            if frac >= 0.99:
                self.consistent_coord = True
                pass
            else:
                self.consistent_coord = False
                
            
            grouped.to_csv( 'temp_data/' + self.dataset  + '/' + self.file.split('/')[-1] + '_coordinates.csv', sep = '\t' )        
                
            lats_good, lons_good= list(compatible_coord.lat.values) , list(compatible_coord.lon.values)
            
        else:
            lats_good, lons_good= lats, lons 
            frequency = '[1]'
            
        # saving list of all lats and lons data
        self.lats_all = lats
        self.lons_all = lons

        self.lats = lats_good
        self.lons = lons_good
        
        self.lat_most_freq = lats_good[0]
        self.lon_most_freq = lons_good[0]
        
        # extract statid 
        comm = "odb sql 'select distinct statid'  -i FILE --no_alignment ".replace('FILE', file)   
        proc = subprocess.Popen(comm, stdout=subprocess.PIPE , shell=True)
        b = proc.stdout.read().decode('utf-8').split('\n')
        statids = [ eval(c) for c in np.unique( b[1:-1] ) ] 
        del b 
        
        # extract distinct variables
        comm = "odb sql 'select distinct varno'  -i FILE --no_alignment ".replace('FILE', file)   
        proc = subprocess.Popen(comm, stdout=subprocess.PIPE , shell=True)
        b = proc.stdout.read().decode('utf-8').split('\n')
        variables = [ eval(c) for c in np.unique( b[1:-1] ) ] 
        variables = [v for v in variables if v in self.odb_to_cdm.keys() ]
        variables = [self.odb_to_cdm[v] for v in variables]
        variables.sort()
        del b
        
        # extract min date, max date 
        comm = "odb sql 'select distinct MIN(date) '  -i FILE --no_alignment ".replace('FILE', file)   
        proc = subprocess.Popen(comm, stdout=subprocess.PIPE , shell=True)
        b = proc.stdout.read().decode('utf-8').split('\n')
        min_date = b[1]
        del b
        
        comm = "odb sql 'select distinct MAX(date) '  -i FILE --no_alignment ".replace('FILE', file)   
        proc = subprocess.Popen(comm, stdout=subprocess.PIPE , shell=True)
        b = proc.stdout.read().decode('utf-8').split('\n')
        
        max_date = b[1]
        
        # TODO check if frequency works!! 
        dic = { 'max_date': max_date , 
                'min_date':min_date, 
                'statid':str(statids), # needed in case multiple ids 
                'lats': [lats_good] , 
                'lons': [lons_good], 
               
                'lats_all': [lats] , 
                'lons_all': [lons], 
                
                'frequency': frequency,
                
                'lats_most_freq':   lats_good[0], 
                'lons_most_freq':  lons_good[0], 

                'file': file , 
                'file_path': self.file , 
                'db': self.dataset,
                'variables': str(variables),
                'consistent_coord': str(self.consistent_coord) }
        
        # storing the 

        pd.DataFrame(dic).to_csv( 'temp_data/' + self.dataset + '/' + self.file.split('/')[-1] + '.csv', sep = '\t' )        
        self.statid = statids
        
        #self.lats = lats
        #self.lons = lons 
        
        self.min_date = min_date 
        self.max_date = max_date
        self.variables = dic['variables']

        
    def check_consistent_coords(self, lats, lons):
        """ Given the entire list of lats and lons check if the values are consistent.
        If not, picks the most frequently occurring pairs, and checks how much of the data is compatible.
        If more than 99% of the data has similar coordinates, the file is considered correct.
        Otherwise it cannot be identified since coordinates are inconsistent in time."""
        
        # saving all values of lat and lon 
        self.lats_all = lats
        self.lons_all = lons
        
        # group by equal lat/lon values 
        df = pd.DataFrame( {'lat': lats , 'lon' : lons } )
        grouped = df.groupby(df.columns.tolist(), as_index=False).size()
        grouped = grouped.sort_values(by=['size'], ascending = False)        
        
        ### If only one row, means all the coordinates are equal
        if len(grouped) == 1:
            self.lat_most_freq = lats[0]
            self.lon_most_freq = lons[0]
            self.lats = lats
            self.lons = lons
            self.consistent_coord = True    
            self.frequency= '[1]'
            return
        
        else:
            ### Case where the coordinate pairs are different 
            ### Must check for the distances in km of the various points and see how many are > 30 km (arbitrary threshold)
            # add distance column
            distances = [ round(self.utils.distance( grouped['lat'].values[i], grouped['lon'].values[i] , 
                                             grouped['lat'].values[0], grouped['lon'].values[0] ), 1 ) for i in range(len(grouped['lat'].values) ) ] 
    
            grouped['distance_km'] = distances
            
            compatible_coord = grouped.loc[grouped['distance_km'] <= self.utils.treshold ]
            
            # saving only compatible values 
            lats_good, lons_good= list(compatible_coord.lat.values) , list(compatible_coord.lon.values)
            self.lats = lats_good
            self.lons = lons_good
            frequency = str(list(grouped['size'].values))  # frequency of the different lat-lon pairs
            self.frequency = frequency  
            
            # best values = most frequent values 
            self.lat_most_freq = lats_good[0]
            self.lon_most_freq = lons_good[0]
            
            frac = sum(compatible_coord['size'].values) / sum(grouped['size'].values )
            if frac >= 0.99:
                self.consistent_coord = True
            else:
                self.consistent_coord = False

    def read_giub(self):
                
        lines = open(self.file,'r').readlines()
        dates = [l.split('\t')[1] for l in lines if 'date' not in l]
        dates = np.unique(dates)        
        
        # must grab the lat and lon from the inventory file CHUAN 
        statid = self.file.split('/')[-1].split('.txt')[0]
        if 'A' in statid or 'B' in statid:
            statid = statid.replace('A','').replace('B','')
            
        statIds = [  statid ]
        
        chuan = self.inventory.inv['chuan'].astype(str)
        loc = chuan.loc[chuan.station_id == statid ]
        lats = [float(loc.latitude.values[0])]
        lons = [float(loc.longitude.values[0])]

    
        dic = { 'max_date': max(dates), 
                 'min_date': min(dates), 
                 'statid': str(statIds), # needed in case multiple ids 
                 'lats': [lats], 
                 'lons': [lons], 

                 'lats_all': [lats], 
                 'lons_all': [lons], 
                 
                 'file': self.file.split('/')[-1] , 
                 'file_path': self.file , 
                 'db': self.dataset,
                 'variables': str( ['126','137','106','107','137', '138', '139', '140', '117'] ),
                 'frequency': [1],
                 'consistent_coord': 'True' }
        
        # set consistent_coord 
        self.consistent_coord = True
        
        pd.DataFrame(dic).to_csv( 'temp_data/'  + self.dataset +  '/' + self.file.split('/')[-1] + '.csv', sep = '\t' )        
        self.statid = statIds
        self.lats = lats
        self.lons = lons 
        self.min_date = min(dates)
        self.max_date = max(dates)
        self.variables = dic['variables']
        return 0    
            
        
        
    def read_amma(self):
        """ Reading AMMA station files writeen into csv station files """
        
        amma_file = self.file 
        
        # reading data via csv 
        data =pd.read_csv(self.file, sep='\t')
        dates = np.unique(data.typicalDate)        
        lats = np.unique(data.latitude.values)
        lons = np.unique(data.longitude.values)
        
        
        df = pd.DataFrame( {'lat': lats , 'lon' : lons } )
        df = df.drop_duplicates().reset_index(drop=True) 
         
        statIds = [ np.unique(data.statid.values)[0]  ]
          
        lats, lons = list(df.lat.values), list(df.lon.values)
        a = self.check_consistent_coords(lats, lons)
    
        dic = { 'max_date': max(dates), 
                 'min_date': min(dates), 
                 'statid': str(statIds), # needed in case multiple ids 
                 'lats': [lats], 
                 'lons': [lons], 

                 'lats_all': [lats], 
                 'lons_all': [lons], 
                 
                 'file': self.file.split('/')[-1] , 
                 'file_path': self.file , 
                 'db': self.dataset,
                 'variables': str( ['126','137','106','107','117'] ),
                 'frequency': self.frequency,
                 'consistent_coord': str(self.consistent_coord) }

        
        pd.DataFrame(dic).to_csv( 'temp_data/'  + self.dataset +  '/' + self.file.split('/')[-1] + '.csv', sep = '\t' )        
        self.statid = statIds
        self.lats = lats
        self.lons = lons 
        self.min_date = min(dates)
        self.max_date = max(dates)
        self.variables = dic['variables']
        return 0    
    
    
    
    def read_hara(self):
        """ Reading HARA station files written into csv station files """
        
        hara_file = self.file 
        
        # reading data via csv 
        data =pd.read_csv(self.file, sep='\t')
        dates = np.unique(data.date_time)     
        dates = [ ''.join(f.split('_')) for f in dates ]
        lats = data.lat.values
        lons = data.lon.values
        
        del data
        
        lats = [float(str(l).replace('\t','')) for l in lats ]
        lons = [float(str(l).replace('\t','')) for l in lons ]
        
        lats_all = np.unique(lats)
        lons_all = np.unique(lons)
        
        df = pd.DataFrame( {'lat': lats , 'lon' : lons } )
        df = df.drop_duplicates().reset_index(drop=True) 
         
        statIds = [ hara_file.split('/')[-1].split('_')[0] ]
          
        #lats, lons = list(df.lat.values), list(df.lon.values)
        
        a = self.check_consistent_coords(lats, lons)
    
        dic = { 'max_date': max(dates), 
                 'min_date': min(dates), 
                 'statid': str(statIds), # needed in case multiple ids 
                 'lats': [list(df.lat.values)], 
                 'lons': [list(df.lon.values)], 

                 'lats_all': [list(lats_all)], 
                 'lons_all': [list(lons_all)],   # actually the same thing...
                 
                 'file': self.file.split('/')[-1] , 
                 'file_path': self.file , 
                 'db': self.dataset,
                 'variables': str( ['126','137','106','107','117'] ),
                 'frequency': self.frequency,
                 'consistent_coord': str(self.consistent_coord) }

        
        pd.DataFrame(dic).to_csv( 'temp_data/'  + self.dataset +  '/' + self.file.split('/')[-1] + '.csv', sep = '\t' )        
        self.statid = statIds
        self.lats = lats
        self.lons = lons 
        self.min_date = min(dates)
        self.max_date = max(dates)
        self.variables = dic['variables']
        return 0    


    def read_bufr_cnr(self):
        """ Reading BURF CNR station files written into csv station files """
        
        bufr_cnr_file = self.file 
        
        # reading data via csv 
        data =pd.read_csv(self.file, sep='\t')
        dates = pd.to_datetime(data['date'].astype('string') + data['time'].astype('string').str.zfill(6), format='%Y%m%d%H%M%S')
        dates = np.unique(dates)     
        lats = data.latitude.values
        lons = data.longitude.values
        
        del data
                
        lats_all = np.unique(lats)
        lons_all = np.unique(lons)
        
        df = pd.DataFrame( {'lat': lats , 'lon' : lons } )
        df = df.drop_duplicates().reset_index(drop=True) 
         
        statIds = [ bufr_cnr_file.split('/')[-1].split('_')[0] ]
          
        #lats, lons = list(df.lat.values), list(df.lon.values)
        
        a = self.check_consistent_coords(lats, lons)
    
        dic = { 'max_date': max(dates), 
                 'min_date': min(dates), 
                 'statid': str(statIds), # needed in case multiple ids 
                 'lats': [list(df.lat.values)], 
                 'lons': [list(df.lon.values)], 

                 'lats_all': [list(lats_all)], 
                 'lons_all': [list(lons_all)],   # actually the same thing...
                 
                 'file': self.file.split('/')[-1] , 
                 'file_path': self.file , 
                 'db': self.dataset,
                 'variables': str( ['126','137','106','107','117'] ),
                 'frequency': self.frequency,
                 'consistent_coord': str(self.consistent_coord) }

        
        pd.DataFrame(dic).to_csv( 'temp_data/'  + self.dataset +  '/' + self.file.split('/')[-1] + '.csv', sep = '\t' )        
        self.statid = statIds
        self.lats = lats
        self.lons = lons 
        self.min_date = min(dates)
        self.max_date = max(dates)
        self.variables = dic['variables']
        return 0    



    def read_shipsound_npsound(self):
        """ Reading SHIPSOUND files written into csv station files """
               
        # reading data via csv 
        data =pd.read_csv(self.file, sep='\t')
        dates = np.unique(data.date_time)     
        dates = [ ''.join(f.split('_')) for f in dates ]
        lats = data.latitude.values
        lons = data.longitude.values
        
        del data
        
        lats = [float(str(l).replace('\t','')) for l in lats ]
        lons = [float(str(l).replace('\t','')) for l in lons ]
        
        lats_all = np.unique(lats)
        lons_all = np.unique(lons)
        
        df = pd.DataFrame( {'lat': lats , 'lon' : lons } )
        df = df.drop_duplicates().reset_index(drop=True) 
         
        statIds = [ self.file.split('/')[-1].split('.dat')[0] ]
          
        #lats, lons = list(df.lat.values), list(df.lon.values)
        
        a = self.check_consistent_coords(lats, lons)
    
        dic = { 'max_date': max(dates), 
                 'min_date': min(dates), 
                 'statid': str(statIds), # needed in case multiple ids 
                 'lats': [list(df.lat.values)], 
                 'lons': [list(df.lon.values)], 

                 'lats_all': [list(lats_all)], 
                 'lons_all': [list(lons_all)],   # actually the same thing...
                 
                 'file': self.file.split('/')[-1] , 
                 'file_path': self.file , 
                 'db': self.dataset,
                 'variables': str( ['126','137','106','107','117'] ),
                 'frequency': self.frequency,
                 'consistent_coord': str(self.consistent_coord) }

        
        pd.DataFrame(dic).to_csv( 'temp_data/'  + self.dataset +  '/' + self.file.split('/')[-1] + '.csv', sep = '\t' )        
        self.statid = statIds
        self.lats = lats
        self.lons = lons 
        self.min_date = min(dates)
        self.max_date = max(dates)
        self.variables = dic['variables']
        return 0  
    
    



    '''
    def read_npsound(self):
        """ Reading NPSOUND files written into csv station files """
               
        # reading data via csv 
        data =pd.read_csv(self.file, sep='\t')
        dates = np.unique(data.date_time)     
        dates = [ ''.join(f.split('_')) for f in dates ]
        lats = data.latitude.values
        lons = data.longitude.values
        
        del data
        
        lats = [float(str(l).replace('\t','')) for l in lats ]
        lons = [float(str(l).replace('\t','')) for l in lons ]
        
        lats_all = np.unique(lats)
        lons_all = np.unique(lons)
        
        df = pd.DataFrame( {'lat': lats , 'lon' : lons } )
        df = df.drop_duplicates().reset_index(drop=True) 
         
        statIds = [ self.file.split('/')[-1].split('.dat')[0] ]
          
        #lats, lons = list(df.lat.values), list(df.lon.values)
        
        a = self.check_consistent_coords(lats, lons)
    
        dic = { 'max_date': max(dates), 
                 'min_date': min(dates), 
                 'statid': str(statIds), # needed in case multiple ids 
                 'lats': [list(df.lat.values)], 
                 'lons': [list(df.lon.values)], 

                 'lats_all': [list(lats_all)], 
                 'lons_all': [list(lons_all)],   # actually the same thing...
                 
                 'file': self.file.split('/')[-1] , 
                 'file_path': self.file , 
                 'db': self.dataset,
                 'variables': str( ['126','137','106','107','117'] ),
                 'frequency': self.frequency,
                 'consistent_coord': str(self.consistent_coord) }

        
        pd.DataFrame(dic).to_csv( 'temp_data/'  + self.dataset +  '/' + self.file.split('/')[-1] + '.csv', sep = '\t' )        
        self.statid = statIds
        self.lats = lats
        self.lons = lons 
        self.min_date = min(dates)
        self.max_date = max(dates)
        self.variables = dic['variables']
        return 0  
    '''
    
        
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
                            try:
                                val = codes_get(bufr, v)
                            except:
                                val = 999
                            data[v].append(val)
                            
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
        
        lats, lons = [round(l,3) for l in data['latitude']] , [round(c,4) for c in data['longitude']]
        df = pd.DataFrame( {'lat': lats , 'lon' : lons } )
        df = df.drop_duplicates().reset_index(drop=True) 
         
        statIds = [np.unique( ['{:0>2}{:0>3}'.format(data['blockNumber'][i], data['stationNumber'][i] ) for i in range(len(data['blockNumber'] ) )  ]  )[0]  ]
         
        
        # maximum 30 km distance 
        a = self.check_consistent_coords(lats, lons)
            
        lats, lons = list(df.lat.values), list(df.lon.values)

        dic = { 'min_date': data['start_date'], 
                 'max_date': data['end_date'] , 
                 'statid': str(statIds), # needed in case multiple ids 
                 #'lats': lats, 
                 #'lons': lons, 
                 'lats' : [self.lats],
                 'lons': [self.lons],
                 
                 'lats_all': [self.lats_all] , 
                 'lons_all': [self.lons_all], 
                 
                 'frequency': self.frequency,
                 
                 'lats_most_freq':   self.lats[0], 
                 'lons_most_freq':  self.lons[0],                  
                 
                 'file': self.file.split('/')[-1] , 
                 'file_path': self.file , 
                 'db': self.dataset,
                 'variables': str(list(np.unique(data['variables'])) ),
                 'consistent_coord': str(self.consistent_coord) }

        pd.DataFrame(dic).to_csv( 'temp_data/' +self.dataset  + '/' + self.file.split('/')[-1] + '.csv', sep = '\t' )        
        self.statid = statIds
        self.min_date = data['start_date'] 
        self.max_date = data['end_date']
        self.variables = dic['variables']
        
        return 0    

    def read_bufr_maestro(self):
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
                        keynames = []
                        while codes_bufr_keys_iterator_next(iterid):     
                            keyname = codes_bufr_keys_iterator_get_name(iterid)
                            keynames.append(keyname)
                            
                            if keyname == 'typicalYear':
                                year = codes_get_array(bufr,keyname)[0]
                            if keyname == 'typicalMonth':
                                month = codes_get_array(bufr,keyname)[0]
                            if keyname == 'typicalDay':
                                day = codes_get_array(bufr,keyname)[0]
                        #datum='{}-{:0>2}-{:0>2}'.format(year, month, day  ) this is the format of the stat conf 
                        datum='{}{:0>2}{:0>2}'.format(year, month, day  )
                        
                        print(keynames)

                        if i == 0:
                            data['start_date'] = datum
                        if i== (cnt-1):
                            data['end_date'] = datum
                            
                        # delete the key iterator
                        codes_bufr_keys_iterator_delete(iterid)
        
                        for v in ['#1#latitude' , '#1#longitude', '#1#stationNumber', '#1#blockNumber']:
                            try:
                                val = codes_get(bufr, v)
                            except:
                                val = 999
                            data[v[3:]].append(val)
                            
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
        
        lats, lons = [round(l,3) for l in data['latitude']] , [round(c,4) for c in data['longitude']]
        df = pd.DataFrame( {'lat': lats , 'lon' : lons } )
        df = df.drop_duplicates().reset_index(drop=True) 
         
        statIds = [np.unique( ['{:0>2}{:0>3}'.format(data['blockNumber'][i], data['stationNumber'][i] ) for i in range(len(data['blockNumber'] ) )  ]  )[0]  ]
         
        
        # maximum 30 km distance 
        a = self.check_consistent_coords(lats, lons)
            
        lats, lons = list(df.lat.values), list(df.lon.values)

        dic = { 'min_date': data['start_date'], 
                 'max_date': data['end_date'] , 
                 'statid': str(statIds), # needed in case multiple ids 
                 #'lats': lats, 
                 #'lons': lons, 
                 'lats' : [self.lats],
                 'lons': [self.lons],
                 
                 'lats_all': [self.lats_all] , 
                 'lons_all': [self.lons_all], 
                 
                 'frequency': self.frequency,
                 
                 'lats_most_freq':   self.lats[0], 
                 'lons_most_freq':  self.lons[0],                  
                 
                 'file': self.file.split('/')[-1] , 
                 'file_path': self.file , 
                 'db': self.dataset,
                 'variables': str(list(np.unique(data['variables'])) ),
                 'consistent_coord': str(self.consistent_coord) }

        pd.DataFrame(dic).to_csv( 'temp_data/' +self.dataset  + '/' + self.file.split('/')[-1] + '.csv', sep = '\t' )        
        self.statid = statIds
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
                
                #NCAR gives longitude as [0,360]
                lo = float(line[68:78])
                if lo > 180:
                    lo = - (360-lo)
                lons.append( lo )
                                
                #idate = datetime.strptime(year + month + day + hour + minutes, '%Y%m%d%H%M')
                
        df = pd.DataFrame( {'lat': lats , 'lon' : lons } )
        df = df.drop_duplicates().reset_index(drop=True) 
         
        statIds = [ np.unique(stations_id)[0]  ]
         
        self.consistent_coord = True        
        distances = [ round(self.utils.distance( df.lat[i], df.lon[i] , 
                                        df.lat[0], df.lon[0] , ), 1 ) for i in range(len(df.lat) ) ]         
        
        a = self.check_consistent_coords(lats, lons)
        
        lats, lons = list(df.lat.values), list(df.lon.values)
    
        dic = { 'max_date': max(dates), 
                 'min_date': min(dates), 
                 'statid': str(statIds), # needed in case multiple ids 
                 'lats': [lats], 
                 'lons': [lons], 

                 'lats_all': [lats], 
                 'lons_all': [lons], 
                 
                 'file': self.file.split('/')[-1] , 
                 'file_path': self.file , 
                 'db': self.dataset,
                 'variables': str( ['38','85','106','107'] ),
                 'frequency': self.frequency,
                 'consistent_coord': str(self.consistent_coord) }

        
        pd.DataFrame(dic).to_csv( 'temp_data/'  + self.dataset +  '/' + self.file.split('/')[-1] + '.csv', sep = '\t' )        
        self.statid = statIds
        self.lats = lats
        self.lons = lons 
        self.min_date = min(dates)
        self.max_date = max(dates)
        self.variables = dic['variables']
        return 0    
    
    
    def read_igra2(self, file=''):
        """ Extract data from a pandas series extracted from igra2 stations_list """
        file = self.file
        df = igra2_ascii_to_dataframe(file)
        variables = ['85','38','106','107']
        
        lat = df['lat@hdr'][0]
        lon = df['lon@hdr'][0]
        # if lat < -180 or lon < -180 :
        #     lat, lon = '', '' 
        
        if not  (lat and lon):
            coord = False
        else:
            coord = True 
            
        dic = { 'max_date': df.iday[1], 
                'min_date':df.iday[0], 
                'statid':str(df['statid@hdr'][0]),
                'lats': [ lat ] , 
                'lons': [ lon ], 
                'file': self.file , 
                'file_path':  self.file , 
                'db': self.dataset,
                'variables': str(variables),
                'consistent_coord': coord }  # NB coordinates should be checked against complete files, not stations lists !!! 
        
        # storing the 
        pd.DataFrame(dic).to_csv( 'temp_data/' + self.dataset + '/'  + self.file.split('/')[-1] + '.csv', sep = '\t' )     
        
        self.statid = [df['statid@hdr'][0]]
        self.lats = [lat]
        self.lons = [lon]
        self.min_date = df.iday[0]
        self.max_date = df.iday[1]
        self.variables = str(variables)
        self.consistent_coord = coord 
        return 0       

    # def read_igra2(self, file=''):
    #     """ Extract data from a pandas series extracted from igra2 stations_list """
    #     series = self.series
    #     variables = ['85','38','106','107']
        
    #     ### must remove dummy lat and lon for mobing ships in igra2
    #     lat, lon = series.latitude, series.longitude
    #     if lat < -180 or lon < -180 :
    #         lat, lon = '', '' 
        
    #     if not  (lat and lon):
    #         coord = False
    #     else:
    #         coord = True 
            
    #     dic = { 'max_date': series.end_date , 
    #             'min_date':series.start_date, 
    #             'statid':str(series.station_id), # needed in case multiple ids 
    #             'lats': [ lat ] , 
    #             'lons': [ lon ], 
    #             'file': self.file , 
    #             'file_path':  self.file , 
    #             'db': self.dataset,
    #             'variables': str(variables),
    #             'consistent_coord': coord }  # NB coordinates should be checked against complete files, not stations lists !!! 
        
    #     # storing the 
    #     pd.DataFrame(dic).to_csv( 'temp_data/' + self.dataset + '/'  + self.file.split('/')[-1] + '.csv', sep = '\t' )     
        
    #     self.statid = [series.station_id]
    #     self.lats = [lat]
    #     self.lons = [lon]
    #     self.min_date = series.start_date
    #     self.max_date = series.end_date
    #     self.variables = str(variables)
    #     self.consistent_coord = coord 
    #     return 0   
    
    
class Inventory():
    """  Class holding the functionalities to read and homogenize 
    the OSCAR, IGRA2, WBAN and CHUAN inventories.
    For OSCAR and WBAN will convert lat and lon to decimal format. """
    
    def __init__(self, datadir ='', tablesdir = '/srvfs/home/uvoggenberger/CEUAS/CEUAS/meta/inventory_comparison_2/data/tables/',
                 oscar ="",
                 igra2= "",
                 wban = '',
                 chuan = '',
                 schroeder='',
                 wmo='',
                 amma='',
                 hara='',
                 utils='' ):
        #self.datadir = datadir + '/' + datadir # TO DO why this ?
        self.datadir = datadir + '/'      
        self.oscar  = datadir + '/' + oscar
        self.igra2  = datadir + '/' + igra2
        self.wban  = datadir + '/' + wban
        self.chuan = datadir + '/' + chuan
        self.schroeder = datadir + '/' + schroeder
        self.wmo = datadir + '/' + wmo
        self.amma = datadir + '/' + amma 
        self.hara = datadir + '/' + hara
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
        
        
        igra2 = igra2[['station_id_igra', 'station_name','latitude','longitude', 'elevation', 'station_id', 'start_date', 'end_date']]
        
        # clean for missing values in lat and lon (i.e. values < 0)
        #igra2 = igra2.loc[ (igra2['latitude'] >= -180) & (igra2['longitude'] >= -180) ] ## TO DO HERE TO DO 
        igra2 = igra2.reset_index(drop=True)
        igra2['isRadio'] = 'True'
        igra2.name = 'IGRA2'

        out = open('/srvfs/home/uvoggenberger/CEUAS/CEUAS/meta/inventory_comparison_2/code/file_list/igra2_files_list.txt' , 'w')
        
        for i in range(len(igra2)):
            stat = igra2.loc[i, 'station_id_igra' ] 
            lat = igra2.loc[i, 'latitude' ] 
            lon = igra2.loc[i, 'longitude' ] 
            
            if lat < -90 or lon < -180:
                k = 'mobile'
            else:
                k = 'regular'
                
            #print(stat + '\t' + k + '\n')
            out.write(stat + '\t' + k + '\n')
            
        #############################################
        ### Reading SCHROEDER data
        #############################################
    
        sch = pd.read_csv(self.schroeder, sep='\t') # columns:  'wmo', 'lat', 'lon', 'station', 'date_min', 'date_max' 
        sch = sch.rename(columns = {'wmo':'station_id','station':'station_name',
                                    'lon':'longitude', 'lat':'latitude'})    
        sch['WMO_id'] = sch.station_id
        sch['isRadio'] = None
        sch['elevation'] = -999  # not available
        sch = sch.reset_index(drop=True)
        sch.name='SCHROEDER'

        #############################################
        ### Reading WMO data
        #############################################
        
        wmo = pd.read_csv(self.wmo, sep='\t') # columns:  'wmo', 'lat', 'lon', 'station', 'date_min', 'date_max'
        wmo = wmo.rename(columns = {'id':'station_id','name':'station_name',
                                    'lons':'longitude', 'lats':'latitude'})
        wmo['WMO_id'] = wmo.station_id
        wmo['isRadio'] = None
        wmo = wmo.reset_index(drop=True)
        wmo.name='WMO'
        #wmo.sch['elevation'] = '' ### TO DO implement !!! it is available

 
        #############################################
        ### Reading AMMA data
        #############################################
        
        amma = pd.read_csv(self.amma, sep='\t') # columns:  'wmo', 'lat', 'lon', 'station', 'date_min', 'date_max'
        amma = amma.rename(columns = {'WMO station No.':'station_id','Station name':'station_name',
                                    'Lon':'longitude', 'Lat':'latitude'})
        amma['WMO_id'] = amma.station_id
        amma['isRadio'] = True
        amma = amma.reset_index(drop=True)
        amma['latitude'] = amma['conv_lat']
        amma['longitude'] = amma['conv_lon']
        cols = [f for f in amma.columns if f not in ['Unnamed: 0', 'GCOS',
       'ASECNA', 'Frequency per day', 'SOP frequency', 'IOP frequency',
       'Upgraded ground station', '2006 Success Rate (%)',
       'Best 2006 Monthly success (%)' ] ] 
        amma = amma[ cols ]
        amma.name='AMMA'
        
        #wmo.sch['elevation'] = '' ### TO DO implement !!! it is available
  
  
  
        #############################################
        ### Reading HARA data
        #############################################
        
        hara = pd.read_csv(self.hara, sep='\t') # columns:  'wmo', 'lat', 'lon', 'station', 'date_min', 'date_max'
        hara = hara.rename(columns = {'STATION_ID':'WMO_id','NAME':'station_name',
                                    'LON_(DEG_E)':'longitude', 'LAT_(DEG_N)':'latitude'})
        hara['WMO_id'] = hara.WMO_id.astype(str)
        ids = ['0'+ str(i) if len(str(i)) == 4 else i for i in hara.WMO_id ]
        hara['WMO_id'] = ids
        hara['station_id'] = ids       
        hara['isRadio'] = True
        hara = hara.reset_index(drop=True)
        cols = [c for c in hara.columns if c not in ['#_OF _SOUNDINGS', 'YEARS',] ]
        hara = hara [cols]
        hara.name='HARA'
        #wmo.sch['elevation'] = '' ### TO DO implement !!! it is available

        #############################################
        ### Reading WBAN data
        #############################################        
        wban =pd.read_fwf(self.wban,widths=(10,5,6,17,22,39,31,31,8,10,10,10,8,7),
                        names=('dum1','station_id','WMO_id','dum0','Country','dum2','station_name','dum3','start_date','end_date',
                                      'latitude','longitude','dum4','elevation'),
                        skiprows=1 )
        
        wban = wban[ [ 'station_id', 'WMO_id' , 'station_name', 'latitude', 'longitude', 'start_date', 'end_date', 'elevation'] ]
        wban = wban.dropna(subset = ['latitude', "longitude"])
        #wban = wban.replace( {'00 00 00': 0 , '0  0  0':0} )
        
        # converting lat, lon to decimal format 
        
        #l = list(self.utils.degMinSec_to_decimal(["-10 56 00"]) )
        wban['elevation']=wban['elevation']*0.3048
        wban['elevation'][wban['elevation']==-99999.0*0.3048]=np.nan
        
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
                                        'Latitude':'latitude' , 'StationName':'station_name', "Hp": 'elevation'} )
        
        oscar = oscar [['WIGOS', 'station_name', 'latitude', 'longitude', 'elevation', 'ObsRems' ]]
        
        # storing original lat and lon (using N,S,W,E notation )
        oscar['original_lat'] = oscar['latitude']
        oscar['original_lon'] = oscar['longitude']
        
        # converting lat, lon to decimal format  
        oscar = self.utils.oscar_coord(oscar) 
    
        # checking if station is a radiosonde 
        radioFlag = oscar['ObsRems'].str.contains('Radio').astype(str)
        statids = [f.split('-')[-1] for f in oscar['WIGOS'] ]  # assume that the station id is the last piece of the WIGOS id 
        oscar['station_id'] = statids
        oscar['start_date'] = 999  # not available in OSCAR 
        oscar['end_date'] = 999 # not available in OSCAR
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

            ## try with Stationname  field or use lat and lon to get country 
            chuan=dict(zip(rdata[0].split('\t'),z))
            chuan['Lat_DegN']  = [float(l.replace(',','.')) for l in chuan['Lat_DegN'] ]
            chuan['Lon_DegE'] = [float(l.replace(',','.')) for l in chuan['Lon_DegE']]

            chuan=pd.DataFrame(chuan)
            chuan=chuan.rename(columns = {'unique_record_ID':'station_id','WMO#':'WMO_id','Stationname':'station_name',
                                          'Lon_DegE':'longitude', 'Lat_DegN':'latitude', 'Alt_masl':'elevation'})

            # filtering missing values for coordinates
            chuan = chuan.loc[ (chuan['latitude'] > -90) & (chuan['longitude'] > -180)]
            chuan = chuan.reset_index(drop=True)
            chuan['station_id']= chuan['station_id'].astype('int')
            #chuan['station_id']= chuan['station_id'].astype('int')

            d = dict( zip( [ str(m) for m in range(1,10)] , [ '0'+ str(d) for d in range (1,10)   ] )  )
            chuan = chuan.replace( {'StartStationMonth': d} )
            chuan = chuan.replace( {'EndStationMonth': d} )
            chuan['start_date'] = chuan['StartStationYear'] + chuan['StartStationMonth']
            chuan['end_date'] = chuan['EndStationYear'] + chuan['EndStationMonth']
            
            chuan = chuan[['station_id', 'WMO_id' , 'station_name' , 'latitude' , 'longitude', 'start_date', 'end_date', 'elevation' ]]
            chuan['isRadio'] = None
            chuan.name = 'CHUAN'
        
        # saving csv files in dedicated inventories directory
        if not os.path.isdir('original_inventories'):
            os.mkdir('original_inventories')
        chuan.to_csv('original_inventories/chuan.csv' , sep = '\t')
        igra2.to_csv('original_inventories/igra2.csv' , sep = '\t')
        oscar.to_csv('original_inventories/oscar.csv' , sep = '\t')
        wban.to_csv('original_inventories/wban.csv' , sep = '\t')
        sch.to_csv('original_inventories/sch.csv' , sep = '\t')
        wmo.to_csv('original_inventories/wmo.csv' , sep = '\t')
        
               
        # storing df as instances
        self.inv["oscar"] = oscar
        self.inv["igra2"] = igra2
        self.inv["wban"] = wban
        self.inv["chuan"] = chuan
        self.inv["schroeder"] = sch
        self.inv["wmo"] = wmo
        self.inv["amma"] = amma
        self.inv["hara"] = hara 
                

# initialize Utils
utils = Utils()

# Initialize inventory
inventory = Inventory( datadir ='/srvfs/home/uvoggenberger/CEUAS/CEUAS/meta/inventory_comparison_2/data/tables',  
                 oscar = "vola_legacy_report.txt" , 
                 igra2 = "igra2-station-list.txt"  , 
                 wban = "WBAN.TXT-2006jan"  , 
                 chuan = "Inventory_ERACLIM_upperair_2.1.txt",
                 schroeder='schroeder_metadata.csv',
                 wmo='wmo_metadata.csv',        
                 amma = 'AMMA_campaign_digitized_metadata_extended.csv',
                 hara = 'HARA_metadata.csv',
                 utils = utils )

inventory.readInventory()
inventory.readCities()

################
### Analysis part 
################
def wrapper(data, file):
    """ Wrapper function to run all the analysis steps for a single file """
    
    dataset = data.dataset 
    if isinstance(file, pd.Series):
        #dataset = 'igra2'
        name_s = file.station_id_igra
    else:
        if 'UADB' in file:
            #dataset = 'ncar'
            name_s = file
        elif 'bfr' in file:
            #dataset = 'bufr'
            name_s = file
        #elif 'TEMP' in file or "PILOT" in file or 'DRIFT' in file:
        #    dataset = 'amma'
        #    name_s = file
        elif 'amma' in file:
            dataset = 'amma'
            name_s = file            
        else:
            #dataset = 'era5_' + file.split('/')[-2]
            name_s = file            
        dataset = data.dataset 
        if isinstance(file, pd.Series):
            #dataset = 'igra2'
            name_s = file.station_id_igra
        else:
            name_s = file
    print("Doing file: " , name_s )
 
    """
    ### creating inventories directory
    out_dir_inv = 'inventories/' + dataset 
    
    if not os.path.isdir(out_dir_inv):
        os.mkdir(out_dir_inv )
    if not os.path.isdir(out_dir_inv + '/logs'):
        os.mkdir(out_dir_inv + '/logs' )
    """
    try:
        print("Reading file::: " , file )
        d = data.read_file(file=file)
        print("Done Reading file::: " , file )
        
        # if data.dataset == 'igra2':
        #     file_name = file.station_id_igra
        # else:
        file_name = file.split("/")[-1]
        
        matching_inv = {}        
        for i in ["oscar","igra2","wban","chuan", 'schroeder', 'wmo', 'amma', 'hara'] :        
            print(' *** Analyzing the inventory: ', i )
            analyze = Analyze( data = data, inv = inventory.inv[i], 
                                cities= inventory.cities, utils = utils)
            
            analyze.AddDistance()
            
            # Find stations with matching ids from the inventory
            analyze.MatchingId()
            # Find closest station in the inventory 
            analyze.Closest()
            analyze.Combine()
            analyze.FindCity()
            matching_inv[i] = analyze.best_match    
            
        l = [matching_inv[f] for f in  matching_inv.keys() ]
        
        all_inventories = pd.concat(l)
        all_inventories = all_inventories.reset_index()
        
        # removing non-matched files and inventories 
        df_red = all_inventories.dropna(subset=['latitude'])
        df_red = df_red.sort_values(by=["inventory", "distance_km"])
        
        # getting internal station id 
        try:
            st = str(eval(data.statid[0][0])[0])
        except:
            st = str(data.statid[0])
            
        if data.dataset in ['igra2'] and not data.consistent_coord: # special case where igra2 has mobile stations
            df_red = pd.DataFrame (columns = df_red.columns )           
            
        if df_red.empty or data.dataset in ['npsound', 'era5_1_mobile' , 'era5_2_mobile' , 'shipsound' ] :
            
            if data.dataset in ['npsound', 'era5_1_mobile' , 'era5_2_mobile' , 'shipsound', 'igra2']:
                flag = 'mobile'
                wigos_pre = '0-20999-0-'
                print("+++ No matching inventory for the file (MOBILE station)")
                
            else:
                flag = 'NoCoordMatchNoIDMatch'
                wigos_pre = '0-20888-0-'
                print("+++ No matching inventory for the file (ORPHAN station)")
            
            all_inventories = {}
            
            try:
                df = pd.DataFrame.from_dict( {'lats':data.lats , 'lons':data.lons} )
                df =df.drop_duplicates()
                df = df.reset_index()
                all_inventories['lat_file'] = ','.join( [str(l) for l in df.lats.values ] )
                all_inventories['lon_file'] = ','.join( [str(l) for l in df.lons.values ] )
            
            except:
                all_inventories['lat_file']  = [-999]
                all_inventories['lon_file'] = [-999]

            all_inventories['file'] = data.file.split('/')[-1]
            all_inventories['file_dataset'] = data.dataset
            all_inventories['file_min_date'] = str(data.min_date)
            all_inventories['file_max_date'] = str(data.max_date)
            all_inventories['file_statid'] = data.statid[0]
            
            all_inventories['variables'] = str(data.variables)  
            
            all_inventories['kind'] = flag 

            a = open('inventories/' + dataset +'/logs/' + dataset + '_processed_all.txt', 'a+')
            a.write(name_s + '\t' + flag + '\t' + '' +  '\n')       

            df = pd.DataFrame( all_inventories, index=[0] )  ### case: fully orphans, no whatsoever matching 
            df['WIGOS_best'] = wigos_pre + st
            name = 'inventories/' + dataset + '/' + data.file.split('/')[-1] + '_' + flag + '.csv'            
            df.to_csv( name,  sep = '\t' , index = False )
            return 
        
        else:
            df_red = df_red.sort_values(by=['inventory','distance_km'])
            df_red['file_dataset'] = data.dataset
            df_red['file_min_date'] = str(data.min_date)
            df_red['file_max_date'] = str(data.max_date)
            df_red['variables'] = str(data.variables)  
            df_red['file'] = data.file.split('/')[-1]
            
            if not data.consistent_coord:  # coordinates are inconsistent, skipping 

                flag = 'inconsistentCoordinates' # to update the global file later 
                wigos_pre = '0-20666-0-'
                
                df_red['WIGOS_best'] = wigos_pre + st
                
                a = open( 'inventories/' + dataset + '/logs/' + dataset + '_' + flag + '.txt', 'a+')
                a.write(file_name + '\t' + str(data.lats_all) + '\t' + str(data.lons_all) + '\t' +  str(data.frequency)  + '\t' + flag +  '\n')
                
                name = 'inventories/' + dataset + '/' + data.file.split('/')[-1] + '_' + flag + '.csv'
                df_red['kind'] = flag
                df_red['file_statid'] = data.statid[0]
                
                df_red.to_csv(name, sep = '\t' , index = False )
                print("+++ No matching inventory for the file (INCONSISTENT COORDS station)")
                
                #return
                
            # Here I finally check that the distance between coordinates with positive lat are too large, 
            # while by adding the minus sing they become comaptible,
            # only valid for stations in the WBAN inventory with same station_id 
            if df_red.distance_km.values[0] > 30:
                try: 
                    df = df_red.loc[df_red["inventory"] == "WBAN"]
                    if not df.empty: 
                        
                        if df.distance_km_minusLat.values[0]:
                            
                            if d <= 30:
                                flag = 'era5-1769minusSignLatMismatch'
                                if not os.path.isfile('inventories/' + dataset + '/logs/' + dataset + '_' + flag +'.txt'):
                                    w = open('inventories/' + dataset + '/logs/' + dataset + '_' + flag + '.txt', 'w' )
                                    
                                    header = "\t".join( ["#station_id","wban_id","file_lat" , "wban_lat" , "file_lon" , "wban_lon" , "file","wigos", '\n'] ) 
                                    w.write(header)                            
                                    w.close()
                                    
                                a = open( 'inventories/' + dataset + '/logs/' + dataset + '_' + flag + '.txt', 'a+' )
                                
                                st = "\t". join( [ str(eval(data.statid[0])[0]) , str( int(df.station_id.values[0])) , str(data.lats[0]) ,  str(df['latitude'].values[0])  ,  str(data.lons[0]) , str(df['longitude'].values[0]) ,  file_name ,  flag , '\n'] ) 
                                a.write(st)

                        
                except:
                    print("PROBLEM==========" , file_name )
                    pass
                
                
        # inventories within 30km ===> IDENTIFIED STATIONS
        df_red = df_red.loc[df_red["distance_km"] <= 30 ]
        
        # Extracting the best wigos id 
        best_wigos = ''
        for i in ['OSCAR', 'IGRA2', 'WBAN', 'CHUAN', 'SCHROEDER', 'WMO', 'AMMA', 'HARA', 'GIUB']:
            if not best_wigos:
                d=df_red.loc[df_red.inventory == i ] 
                if not d.empty:
                    if i == "OSCAR":
                        if len(d[d.isRadio=="True"]) >0:  # best WIGOS would be the radiosonde 
                            best_wigos = d[d.isRadio=="True"].WIGOS.values[0]
                        else:
                            best_wigos = d.WIGOS.values[0]  # closest WIGOS otherwise (they are sorted by distance)
                    else:
                        best_wigos = d.WIGOS_calc.values[0]
            
        if isinstance(data.statid[0], str):
            stat = data.statid[0]
        elif isinstance(data.statid, list):
            stat = str(data.statid[0])
                
        else:
            try:
                stat = str(eval(data.statid[0].values[0])[0] )
            except:
                stat = str(eval(data.statid[0].values[0]))
        # 
        if data.dataset == 'era5_1759':
            try:
                oscar_dist = df_red.loc[df_red["inventory"] == 'OSCAR']["distance_km"].values[0]
            except:
                oscar_dist = 999
            try:
                wban_dist =   df_red.loc[df_red["inventory"] == 'WBAN']["distance_km"].values[0]
            except:
                wban_dist = 999
            
            if not df_red.empty:
                a = open( 'inventories/' + dataset + '/logs/' + dataset + '_CHECK_OSCAR_SignLatMismatch.txt', 'a+' )        
                #st = "\t". join( [ stat, str( df_red.station_id.values[0]) , str(data.lats[0]) ,  str(df_red['latitude'].values[0])  ,  str(data.lons[0]) , str(df_red['longitude'].values[0]) , str(oscar_dist), str(wban_dist),   file_name ,  '\n'] ) 
                st = "\t". join( [ stat, str(data.lats[0]) ,  str(df_red['latitude'].values[0])  ,  str(data.lons[0]) , str(df_red['longitude'].values[0]) , str(oscar_dist), str(wban_dist),   file_name ,  '\n'] ) 
                
                a.write(st)
                
        # Add file information     
        all_inventories['lat_file'] = str( [max(data.lats), min(data.lats)] )
        all_inventories['lon_file'] = str( [max(data.lons), min(data.lons)]  )
        all_inventories['file'] = data.file.split('/')[-1]
        all_inventories['file_dataset'] = data.dataset
        all_inventories['file_min_date'] = str(data.min_date)
        all_inventories['file_max_date'] = str(data.max_date)
        all_inventories['file_statid'] = data.statid[0]
        
        all_inventories['WIGOS_best'] = best_wigos  
        all_inventories['variables'] = str(data.variables)  
        
        name = 'inventories/' + dataset + '/' + data.file.split('/')[-1] + '_inventories.csv'
        
        all_inventories_red = all_inventories[ ['file_statid', 'station_id', 'station_name', 'latitude', 'longitude', 'elevation',
                                                'original_lat', 'original_lon', 'distance_km', 'lat_file',
                                                'lon_file',  'inventory', 'WIGOS', 'WMO_id', 'WIGOS_calc', 'file_min_date', 
                                                'file_max_date', 'start_date', 'end_date', 'isRadio', 'WIGOS_best', 'city', 'variables'] ]
        
        # NOTE: will save all inventories, even the ones whihc do not satisfy distance requirement 
        
        if len(df_red) >0 and data.consistent_coord:
            flag = 'identified' + '\t' + best_wigos 
            all_inventories['kind'] = 'identified'
            all_inventories.to_csv( name.replace('.csv','_identified.csv'), sep = '\t', index = False)
            all_inventories_red.to_csv( name.replace('.csv','_reduced.csv'), sep = '\t' , index = False )

        elif data.consistent_coord :      
            all_inventories['kind'] = 'noDistMatch'
            
            print("+++ No matching inventory wrt distance of 30 km  for the file but matching existing ID from inventory")
            out = open( 'inventories/' + dataset + '/logs/' + dataset + '_unidentified.txt', 'a+')
            out.write( name_s +   '\n' )
            
            flag = 'unidentified_nomatch' + '\t' + ' '
            
            all_inventories['WIGOS_best'] = '0-20777-0-' + st
            all_inventories.to_csv( name.replace('.csv', '_noDistMatch.csv'), sep = '\t', index = False)            
            all_inventories_red.to_csv( name.replace('.csv','_noDistMatch_reduced.csv'), sep = '\t' , index = False )
            
        # writing both identified and non-identified files to log
        a = open(  'inventories/' + dataset + '/logs/' + dataset + '_processed_all.txt', 'a+')        
        a.write(name_s + '\t' + str(data.lats[0]) + '\t' + str(data.lons[0]) + '\t' + flag + '\n' )       
        
        print("Done :::" , file_name )


    except:
        print("*** Cannot read file! ***" , name_s )
        a = open('inventories/' + dataset + '/logs/' + dataset +"_failed_files.txt", "a+")
        a.write(name_s + '\n')       
        
        a = open('inventories/' + dataset + '/logs/' + dataset + '_processed_all.txt', 'a+')
        a.write(name_s + '\t' + 'failed' + '\n')               
    
    
    
# missing files in the mismatch 
test_era5_1759_lat_mismatch_all = ['era5.1759.conv.2:82606', 'era5.1759.conv.2:82607', 'era5.1759.conv.2:80430', 'era5.1759.conv.2:82608', 'era5.1759.conv.2:80431', 'era5.1759.conv.2:80433', 
                                   'era5.1759.conv.2:80505', 'era5.1759.conv.2:80506', 'era5.1759.conv.2:80509', 'era5.1759.conv.2:80510', 'era5.1759.conv.2:80704', 
                                   'era5.1759.conv.2:81302', 'era5.1759.conv.2:81404', 'era5.1759.conv.2:81408', 'era5.1759.conv.2:81410', 'era5.1759.conv.2:81502', 
                                   'era5.1759.conv.2:81605', 'era5.1759.conv.2:50302', 'era5.1759.conv.2:51301', 'era5.1759.conv.2:51302', 'era5.1759.conv.2:51304', 
                                   'era5.1759.conv.2:50801', 'era5.1759.conv.2:52401', 'era5.1759.conv.2:50802', 'era5.1759.conv.2:60501', 'era5.1759.conv.2:60701', 
                                   'era5.1759.conv.2:61503', 'era5.1759.conv.2:60702', 'era5.1759.conv.2:61703', 'era5.1759.conv.2:61502', 'era5.1759.conv.2:62702', 
                                   'era5.1759.conv.2:50304', 'era5.1759.conv.2:61701', 'era5.1759.conv.2:50305', 'era5.1759.conv.2:80402', 'era5.1759.conv.2:62701', 
                                   'era5.1759.conv.2:50306', 'era5.1759.conv.2:80604', 'era5.1759.conv.2:80405', 'era5.1759.conv.2:81306', 'era5.1759.conv.2:50307', 
                                   'era5.1759.conv.2:80406', 'era5.1759.conv.2:81403', 'era5.1759.conv.2:50309', 'era5.1759.conv.2:80408', 'era5.1759.conv.2:50310', 
                                   'era5.1759.conv.2:81405', 'era5.1759.conv.2:80409', 'era5.1759.conv.2:50402', 'era5.1759.conv.2:81406', 'era5.1759.conv.2:80412', 
                                   'era5.1759.conv.2:81407', 'era5.1759.conv.2:50403', 'era5.1759.conv.2:81501', 'era5.1759.conv.2:50405', 'era5.1759.conv.2:80419', 
                                   'era5.1759.conv.2:81608', 'era5.1759.conv.2:50406', 'era5.1759.conv.2:80502', 'era5.1759.conv.2:81702', 'era5.1759.conv.2:51303', 
                                   'era5.1759.conv.2:80504', 'era5.1759.conv.2:61501', 'era5.1759.conv.2:80511', 'era5.1759.conv.2:82401', 'era5.1759.conv.2:80301', 
                                   'era5.1759.conv.2:80512', 'era5.1759.conv.2:82404', 'era5.1759.conv.2:80302', 'era5.1759.conv.2:80513', 'era5.1759.conv.2:82408', 
                                   'era5.1759.conv.2:83405', 'era5.1759.conv.2:80303', 'era5.1759.conv.2:83504', 'era5.1759.conv.2:80304', 'era5.1759.conv.2:80602', 
                                   'era5.1759.conv.2:83508', 'era5.1759.conv.2:80306', 'era5.1759.conv.2:80310', 'era5.1759.conv.2:80702', 'era5.1759.conv.2:80311', 
                                   'era5.1759.conv.2:81301', 'era5.1759.conv.2:80413', 'era5.1759.conv.2:81304', 'era5.1759.conv.2:80415', 'era5.1759.conv.2:81409', 
                                   'era5.1759.conv.2:80421', 'era5.1759.conv.2:81601', 'era5.1759.conv.2:80422', 'era5.1759.conv.2:81603', 'era5.1759.conv.2:80424', 
                                   'era5.1759.conv.2:81606', 'era5.1759.conv.2:80425', 'era5.1759.conv.2:81609', 'era5.1759.conv.2:80426', 'era5.1759.conv.2:80427', 
                                   'era5.1759.conv.2:82403', 'era5.1759.conv.2:80429', 'era5.1759.conv.2:82405', 'era5.1759.conv.2:82501', 'era5.1759.conv.2:82503', 
                                   'era5.1759.conv.2:82507', 'era5.1759.conv.2:80432', 'era5.1759.conv.2:50404', 'era5.1759.conv.2:80305', 'era5.1759.conv.2:80309', 
                                   'era5.1759.conv.2:80705', 'era5.1759.conv.2:81604']



### Directory containing the databases 
basedir = '/mnt/users/scratch/leo/scratch/era5/odbs/'
uvdir = '/mnt/users/scratch/uvoggenberger/'
datasets = {  'era5_1': basedir + '/1' ,
             'era5_1_mobile' : basedir + '1_mobile' ,
            'era5_2': basedir + '/2',
            'era5_2_mobile': basedir + '/2',

            'era5_3188': basedir +'/3188',
            'era5_1759': basedir +'/1759',
            'era5_1761': basedir + '/1761',
            'bufr': basedir + '/ai_bfr/',  
            'bufr_cnr': uvdir + 'bufr_cnr/concated',                                  
            'ncar': '/scratch/das/federico/databases_service2/UADB_22012024/',
            'amma': '/scratch/das/federico/databases_service2/AMMA_BUFR/AMMA_split_csv/' ,
            'igra2': '', # dummy, use the igra2 station list file 
            'igra2_mobile': '', # dummy, use the igra2 station list file 

            'hara': '/scratch/das/federico/databases_service2/HARA-NSIDC-0008_csv/',
            'maestro': uvdir + 'HARVEST_2025/MAESTRO_2024/IUSN74/',
            
            'giub': '/scratch/das/federico/databases_service2/GIUB_07072023/',
            'npsound' : '/scratch/das/federico/databases_service2/NPSOUND-NSIDC0060/NPSOUND_csv_converted' ,
            'shipsound' : '/scratch/das/federico/databases_service2/SHIPSOUND-NSIDC-0054/'
                               } 


def get_flist(db):
    """ Extarcts the list of files to be processed according to the dataset """
    if not os.path.isdir('file_list'):
        os.mkdir('file_list')
    
    if db == 'era5_1':
        if not os.path.isfile('file_list/era5_1_files_list.txt'):  ### TODO, this is to speed ud reading the file names which take a lot of time
            flist=glob.glob(datasets[db] + "/era5.conv._*") # takes too long 
            flist =[f for f in flist if '_40179' not in f  and '42147' not in f] # thse files cause problems ???
            a = open( 'file_list/era5_1_files_list.txt','w')
            for l in flist:
                a.write(l + '\n')
            a.close()
            
        else:
            flist = [ f.replace('\n','') for f in open(db + '_files_list.txt').readlines() ]
            
            
    elif db == 'era5_1_mobile':
        # here we do not have the full odb file so we cannot use the already implemented function... ok...
        
        if not os.path.isfile('file_list/era5_1_mobile_files_list.txt'):  ### TODO, this is to speed ud reading the file names which take a lot of time
            flist_all=glob.glob(datasets[db] + "/era5.conv.*._*.txt.gz") # takes too long 
            
            stations = list( np.unique( [f.split('.txt')[0].split('_')[-1] for f in flist_all ] ) )
                
            stations.sort()
            flist = []
            a = open( 'file_list/era5_1_mobile_files_list.txt','w')
            for station in stations:
                
                files = [f for f in flist_all if station in f.split('.txt.gz')[0].split('.')[-1] ]
                files.sort()
                file = files[0] 
                
                a.write(file + '\n')
                
                flist.append(file)
            a.close()        
        else:
            flist = [ f.replace('\n','') for f in open('file_list/'+db + '_files_list.txt').readlines() ]
            

    elif db in [ 'era5_2', 'era5_2_mobile' ]:
        # for this dataset, we must read the csv file and determine if it is a mobile or fixed station 
        
        ### ANALSYIS of ERA5 2 stations type (fixed, mobile)
        if not os.path.isfile('file_list/era5_2_files_list.csv' ):  ### TODO, this is to speed ud reading the file names which take a lot of time
            flist=glob.glob(datasets[db] + "/era5.conv._*") 
            flist = [f for f in flist if '.gz' not in f ]
            stations = list( np.unique( [f.split('.txt')[0].split('_')[-1] for f in flist ] ) )
            stations.sort()
                
            ### saving extra file information
            dic = {"file":[], 
                   "reportype": [],
                   'kind' :[] }
            
            for s in tqdm(stations):
                if s in ['','00000'] :
                    continue
                station_file = glob.glob(datasets[db] + "/era5.conv._" + s + '.gz')[0]

                
                df = pd.read_csv(station_file, sep = '\t' , usecols=['reportype'], nrows=50000 )    
                report_type = np.unique(df.reportype).astype(str)
                rt = ','.join(report_type)
                dic['file'].append(station_file)
                dic['reportype'].append(rt)
                
                rt_fixed = [ '16022' , '16013']
                rt_mobile = [r for r in report_type if r not in rt_fixed ]
                
                if len(rt_mobile) == 0:
                    dic['kind'].append('FIXED')
                else:
                    dic['kind'].append('MOBILE')
                    
            df = pd.DataFrame(dic)
            df.to_csv( 'file_list/era5_2_files_list.csv' , sep = '\t')

 
        else:
            df = pd.read_csv( 'file_list/era5_2_files_list.csv', sep='\t' )
            flist = list(df.file) 
            
        if 'mobile' in db:
            flist = list( df.loc[df.kind == 'MOBILE'].file )
        else:
            flist = list( df.loc[df.kind == 'FIXED'].file )

        flist=[f.replace('.gz', '')  for f in flist]            
            
            
    elif db == 'era5_1759':
        flist=glob.glob("/mnt/users/scratch/leo/scratch/era5/odbs/1759/era5.1759.conv.*")
        
    elif db == 'era5_1761':
        flist=glob.glob("/mnt/users/scratch/leo/scratch/era5/odbs/1761/era5.1761.conv.*")
        
    elif db == 'era5_3188':
        flist=glob.glob("/mnt/users/scratch/leo/scratch/era5/odbs/3188/era5.3188.conv.*")
        
    elif db == 'ncar':
        flist=glob.glob( datasets[db] + '/uadb*')
        
    elif db== 'igra2':
        flist=glob.glob("/mnt/users/scratch/uvoggenberger/HARVEST_2025/IGRA2_05012025/*.txt")
        # f = inventory.inv['igra2']
        # flist = [ f.iloc[n] for n in range(len(f)) ]
        
    elif 'bufr_cnr' in db:
        flist=glob.glob(datasets[db] + '/*') 

    elif 'bufr' in db:
        flist=glob.glob(datasets[db] + '/'+'era5.*.bfr')
        flist=[f for f in flist if 'undef' not in f ]

    elif 'maestro' in db:
        flist=glob.glob(datasets[db] + '/'+'*.bufr')

    elif 'amma' in db:
        flist=glob.glob(datasets[db] + '/'+'*')
        
    elif 'hara' in db:
        flist=glob.glob(datasets[db] + '/'+'*')      

    elif 'giub' in db:
        flist=glob.glob(datasets[db] + '/'+'*')         
        
    elif 'npsound' in db:
        flist=glob.glob(datasets[db] + '/'+'*')         
        #flist = [f for f in flist if '28' not in f and '31' not in f ]  # dont remember why? 
        
    elif 'shipsound' in db:
        flist=glob.glob(datasets[db] + '/'+'*.csv')         
        #flist = [f for f in flist if '.out' not in f and '31' not in f ]
        
    if db not in ['era5_1_mobile']:
        flist = [f for f in flist if '.gz' not in f and '.nc' not in f ]
        
    flist = [f for f in flist if '00000' not in f and '99999' not in f and '-1e+100' not in f ]
        
    return flist 
        
        
    
if __name__ == '__main__':
    """ Parameters:
          - POOL: runs multiprocesses (default=30)
          - CHECK_MISSING: only runs missing files, otherwise will rerun and replace existing files """

    all_db = [ 'era5_1', 'era5_2', 
               'era5_1759', 'era5_1761', 'era5_3188', 
               'bufr', 'ncar', 
               'igra2', 'igra2_mobile', 
               'amma', 
               'hara', 'npsound', 'shipsound',
               'giub', 
               'maestro',
               'era5_1_mobile' , 'era5_2_mobile']  # era5_1_mobile, era5_2_mobile 
    


    databases = all_db
    
    databases = ['era5_1759', 'era5_1761', 'era5_3188', 
               'bufr', 'bufr_cnr', 'ncar', 
               'igra2', 'igra2_mobile', 
               'amma', 
               'hara', 'npsound', 'shipsound',
               'giub', 'maestro',]
    
    databases = ['maestro']  
    
    # enable multiprocesing
    POOL = False 
    n_pool = 40
    CHECK_MISSING = False  
    CHECK_FAILED = False
    # loop through each of the databases
    for db in databases:
        
        if db == 'igra2':
            CHECK_MISSING = False
            
        ### Preparing directories 
        if not os.path.isdir( 'inventories/' + db + '/logs'):
            os.makedirs( 'inventories/' + db + '/logs' )
        if not os.path.isdir( 'temp_data/' + db + '/logs'):
                os.makedirs( 'temp_data/' + db + '/logs' )            
                                         
        ### Extracting files list to process 
        flist = get_flist(db)
        

        '''
        flist_n = [f.split('_inventories')[0] for f in os.listdir('inventories/' + db ) if 'noMatchingCoordNoIds'  in f]
        flist = [f for f in flist if f.split('/')[-1] in flist_n  ]
        '''
        
        # ### filtering mobile igra
        if db == 'igra2_mobile':
            mobile_igra_list = pd.read_csv('/srvfs/home/uvoggenberger/CEUAS/CEUAS/meta/inventory_comparison_2/code/file_list/igra2_files_list.txt' , sep = '\t', names = ['station', 'kind'] )
            mobile_igra_list = mobile_igra_list.loc[mobile_igra_list.kind == 'mobile']
            flist = [f for f in flist if f.station_id_igra in list(mobile_igra_list.station) ]
        
        
        if CHECK_FAILED:
            failed = open(db.replace('era5_','') + '_failed_files.txt','r').readlines()
            flist = a

        if db != 'igra2':   
            if CHECK_MISSING:
                flist_c = [f.split('/')[-1] for f in flist ]
                try: 
                    processed = [f.replace(".csv","") for f in os.listdir('temp_data/'+db)] 
                except:
                    processed = []
                missing_stat = [f for f in flist_c if f.replace('.csv','') not in processed ]
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
        
        data = Data(dataset=db, utils = utils, inv=inventory )
        

        ######  TODO edit here to possibly filter file list
        # flist = [f for f in flist if '94146' in f ]  # inconsitentCoord era5_1
        
        if db != 'igra2':
            skipped = [ 'era5.1759.conv.xxx'.replace('xxx', str(i))  for i in range(1900,2020,1) ]
            skipped = skipped + [ 'era5.1761.conv.xxx'.replace('xxx', str(i))  for i in range(1900,2020,1) ]
            flist = [f for f in flist if f.split('/')[-1] not in skipped ]
        
        # remove 0 size files:
        flist = [f for f in flist if os.path.getsize(f) != 0]
    
        # flist = flist[:100] ## edit! 
        if len(flist) ==1:
            POOL = False
            
        #fii = [f.replace('_coordinates','').replace('.csv','') for f in os.listdir('/users/staff/federico/GitHub/CEUAS_master_SEPTEMBER2021/CEUAS/CEUAS/meta/inventory_comparison/code/temp_data/' + db) ]
        #flist = [f for f in flist if f.split('/')[-1]  not in fii ]
        
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


# /mnt/users/scratch/leo/scratch/era5/odbs/2/era5.conv._94791


Problematic files:

- ERA5 3188
/mnt/users/scratch/leo/scratch/era5/odbs/3188/era5.3188.conv.C:4724
/mnt/users/scratch/leo/scratch/era5/odbs/3188/era5.3188.conv.C:5975
/mnt/users/scratch/leo/scratch/era5/odbs/3188/era5.3188.conv.C:5583
/mnt/users/scratch/leo/scratch/era5/odbs/3188/era5.3188.conv.C:5752
/mnt/users/scratch/leo/scratch/era5/odbs/3188/era5.3188.conv.C:5688
/mnt/users/scratch/leo/scratch/era5/odbs/3188/era5.3188.conv.
/mnt/users/scratch/leo/scratch/era5/odbs/3188/era5.3188.conv.C:4964

(missing lat and lon in the ODB)

- ERA 1761
/mnt/users/scratch/leo/scratch/era5/odbs/1761/era5.1761.conv.2:34014

- ERA 1759
/mnt/users/scratch/leo/scratch/era5/odbs/1759/era5.1759.conv.3:0N1
/mnt/users/scratch/leo/scratch/era5/odbs/1759/era5.1759.conv.6:99999
/mnt/users/scratch/leo/scratch/era5/odbs/1759/era5.1759.conv.2:34014  -> latitude not in [-90,90] range
/mnt/users/scratch/leo/scratch/era5/odbs/1759/era5.1759.conv.2:32001  -> latitude not in [-90,90] range

- BUFR
/mnt/users/scratch/leo/scratch/era5/odbs/ai_bfr/era5.undefundef.bfr

- ERA5_1 MOBILE
        era5.conv._98753_incosistentCoord.csv
        era5.conv._98501_inventories_noMatchingCoordNoIds.csv
        era5.conv._98222_inventories_noDistMatch.csv



# lat mismatch era5_1759 example
# '/mnt/users/scratch/leo/scratch/era5/odbs/1759/era5.1759.conv.2:21502'
"""

"""
Distances between the station coordinates and the inventories are calculated only if the lat and lon difference are below 2 degrees.
If no inventory station satisfies this condition, distances are not calculated.
The file: logs/
"""

### to try
# inconsistent coordinates
# era5.55299.bfr 


""" IGRA2 MOBILE stations
- they change location in time, hence they do not have defined lat and lo in the igra-station-list file 
/scratch/das/federico/databases_service2/IGRA2_20211231/
ZZV0000DBBH -98.8888 -998.8888 -998.8    METEOR                         1979 2005   3788
ZZV0000DBFM -98.8888 -998.8888 -998.8    MEERKATZE                      1985 1988    236
ZZV0000DBLK -98.8888 -998.8888 -998.8    POLARSTERN                     1985 2023  10383
ZZV0000DCMS -98.8888 -998.8888 -998.8    DOLORES                        1991 1993    530
ZZV0000DESI -98.8888 -998.8888 -998.8    VALDIVIA                       1993 1993    124
ZZV0000DFCG -98.8888 -998.8888 -998.8    SONNE                          2009 2015    254
ZZV0000DKIW -98.8888 -998.8888 -998.8    MERKUR PORTUGAL                1989 1993    735
ZZV0000DZHA -98.8888 -998.8888 -998.8    MANILA BAY                     1990 1991    672
ZZV0000EBUQ -98.8888 -998.8888 -998.8    ESPERANZA DEL MAR              2003 2007    719
ZZV0000EREB -98.8888 -998.8888 -998.8    BRIZ                           1973 1990   2917
ZZV0000EREI -98.8888 -998.8888 -998.8    OKEAN                          1973 1990   3291
ZZV0000ERES -98.8888 -998.8888 -998.8    VIKTOR BUGAYEV                 1974 1991   2110
ZZV0000FNOR -98.8888 -998.8888 -998.8    FORT ROYAL                     1986 2002   2494
ZZV0000FNOU -98.8888 -998.8888 -998.8    FORT FLEUR D'EPEE              1986 2002   2530
ZZV0000FNPH -98.8888 -998.8888 -998.8    FORT DESAIX                    1987 2003   2579
ZZV0000FNRS -98.8888 -998.8888 -998.8    DOUCE FRANCE                   1986 2003   2616
ZZV0000FQFL -98.8888 -998.8888 -998.8    CMA-CGM FORT SAINT LOUIS       2003 2006    522
ZZV0000FQFM -98.8888 -998.8888 -998.8    CMA-CGM FORT SAINT PIERRE      2003 2006    501
ZZV0000JBOA -98.8888 -998.8888 -998.8    KEIFU MARU                     1973 2000   4131
ZZV0000JCCX -98.8888 -998.8888 -998.8    CHOFU MARU                     1988 2009   1634
ZZV0000JDWX -98.8888 -998.8888 -998.8    KOFU MARU                      1989 2009   1243
ZZV0000JGQH -98.8888 -998.8888 -998.8    RYOFU MARU                     2000 2023   3133
ZZV0000JIVB -98.8888 -998.8888 -998.8    SEIFU MARU                     2000 2008    844
ZZV0000JNSR -98.8888 -998.8888 -998.8    MIRAI                          2000 2023   5430
ZZV0000KAOU -98.8888 -998.8888 -998.8    ROGER REVELLE (AWS)            2011 2012    581
ZZV0000KHRH -98.8888 -998.8888 -998.8    SEALAND DEVELOPER              2003 2004    398
ZZV0000KNIJ -98.8888 -998.8888 -998.8    MANULANI                       1985 1987    725
ZZV0000LDWR -98.8888 -998.8888 -998.8    POLARFRONT                     1990 2009  15239
ZZV0000LGKI -98.8888 -998.8888 -998.8    LANCE                          2015 2015    255
ZZV0000SWJS -98.8888 -998.8888 -998.8    PELJASPER                      2000 2003    439
ZZV0000UBLF -98.8888 -998.8888 -998.8    AKADEMIK KURCHATOV             1974 1979    142
ZZV0000UBNZ -98.8888 -998.8888 -998.8    AKADEMIK SHULEYKIN             1983 1991   1389
ZZV0000UFTA -98.8888 -998.8888 -998.8    VOLGONEFT-131                  2008 2012   1068
ZZV0000UHQS -98.8888 -998.8888 -998.8    AKADEMIK KOROLYOV              1974 1990   2966
ZZV0000UVMJ -98.8888 -998.8888 -998.8    VSEVOLOD BERYOZKIN             1977 1983    665
ZZV0000UZGH -98.8888 -998.8888 -998.8    PASSAT                         1974 1991   2079
ZZV0000V2LV -98.8888 -998.8888 -998.8    UAL TEXAS                      1993 1993    144
ZZV0000V2XM -98.8888 -998.8888 -998.8    SKOGAFOSS                      2003 2007    431
ZZV0000V2XO -98.8888 -998.8888 -998.8    LAGARFOSS                      2001 2002    192
ZZV0000VPHA -98.8888 -998.8888 -998.8    OOCL CHALLENGE                 1988 1992   1158
ZZV0000WAAH -98.8888 -998.8888 -998.8    SEALAND MOTIVATOR              2004 2005    361
ZZV0000WPKD -98.8888 -998.8888 -998.8    GALVESTON BAY                  2001 2005   1233
ZZV0000WTEC -98.8888 -998.8888 -998.8    RONALD H. BROWN (AWS)          2000 2020    876
ZZV0000ZDLG -98.8888 -998.8888 -998.8    BRANSFIELD                     1986 1993    288
ZZV0000ZSAF -98.8888 -998.8888 -998.8    S. A. AGULHAS                  1989 2008    325
ZZV000ASDK3 -98.8888 -998.8888 -998.8    ASDK3                          2009 2015   1929
ZZV000ASFR1 -98.8888 -998.8888 -998.8    ASFR1                          2007 2018   3128
ZZV000ASFR2 -98.8888 -998.8888 -998.8    ASFR2                          2007 2018   2610
ZZV000ASFR3 -98.8888 -998.8888 -998.8    ASFR3                          2010 2018   2280
ZZV000ASFR4 -98.8888 -998.8888 -998.8    ASFR4                          2010 2018   2258
ZZV000ELML7 -98.8888 -998.8888 -998.8    HORNBAY                        2000 2005   2309
ZZV000LADB2 -98.8888 -998.8888 -998.8    SKAUGRAN                       1992 1993    394
ZZV000LAJV4 -98.8888 -998.8888 -998.8    SKAUBRYN                       1991 1992    163
ZZV000OVYA2 -98.8888 -998.8888 -998.8    ARINA ARCTICA                  2000 2003    863
ZZV000OXTS2 -98.8888 -998.8888 -998.8    IRENA ARCTICA                  2000 2007    616
ZZV000OXVH2 -98.8888 -998.8888 -998.8    NAJA ARCTICA                   2003 2006    631
ZZV000OXYH2 -98.8888 -998.8888 -998.8    NUKA ARCTICA                   2000 2005   1339
ZZV000P3OO4 -98.8888 -998.8888 -998.8    DRAGON BAY                     1992 1993    148
ZZV000VSBV3 -98.8888 -998.8888 -998.8    CANMAR AMBASSADOR              1992 1992    250
ZZV000ZCBP6 -98.8888 -998.8888 -998.8    MISSISSAUGA EXPRESS            2000 2005   1069
ZZV00ASDE01 -98.8888 -998.8888 -998.8    ASDE01                         2006 2017   3234
ZZV00ASDE02 -98.8888 -998.8888 -998.8    ASDE02                         2006 2017   2737
ZZV00ASDE03 -98.8888 -998.8888 -998.8    ASDE03                         2006 2017   2831
ZZV00ASDE04 -98.8888 -998.8888 -998.8    ASDE04                         2006 2016   3046
ZZV00ASDE09 -98.8888 -998.8888 -998.8    ASDE09                         2007 2023    598
ZZV00ASDK01 -98.8888 -998.8888 -998.8    ASDK01                         2007 2021   2963
ZZV00ASDK02 -98.8888 -998.8888 -998.8    ASDK02                         2007 2017   2993
ZZV00ASDK03 -98.8888 -998.8888 -998.8    ASDK03                         2015 2017    745
ZZV00ASES01 -98.8888 -998.8888 -998.8    ASES01                         2008 2017   1839
ZZV00ASEU01 -98.8888 -998.8888 -998.8    ASEU01                         2006 2017   2124
ZZV00ASEU02 -98.8888 -998.8888 -998.8    ASEU02                         2007 2017   2288
ZZV00ASEU03 -98.8888 -998.8888 -998.8    ASEU03                         2007 2017   2118
ZZV00ASEU04 -98.8888 -998.8888 -998.8    ASEU04                         2006 2017   2028
ZZV00ASEU05 -98.8888 -998.8888 -998.8    ASEU05                         2006 2017   2300
ZZV00ASEU06 -98.8888 -998.8888 -998.8    ASEU06                         2011 2017   1404
ZZV00ASGB01 -98.8888 -998.8888 -998.8    ASGB01                         2007 2011    956
ZZXUAICE002 -98.8888 -998.8888 -998.8    NP02                           1950 1950    178
ZZXUAICE003 -98.8888 -998.8888 -998.8    NP03                           1954 1955    354
ZZXUAICE004 -98.8888 -998.8888 -998.8    NP04                           1954 1957   1089
ZZXUAICE005 -98.8888 -998.8888 -998.8    NP05                           1955 1956    345
ZZXUAICE006 -98.8888 -998.8888 -998.8    NP06                           1956 1959    853
ZZXUAICE007 -98.8888 -998.8888 -998.8    NP07                           1957 1959    680
ZZXUAICE008 -98.8888 -998.8888 -998.8    NP08                           1959 1961    890
ZZXUAICE009 -98.8888 -998.8888 -998.8    NP09                           1960 1961    269
ZZXUAICE010 -98.8888 -998.8888 -998.8    NP10                           1961 1964    873
ZZXUAICE011 -98.8888 -998.8888 -998.8    NP11                           1962 1963    327
ZZXUAICE012 -98.8888 -998.8888 -998.8    NP12                           1963 1965    674
ZZXUAICE013 -98.8888 -998.8888 -998.8    NP13                           1964 1967   1006
ZZXUAICE014 -98.8888 -998.8888 -998.8    NP14                           1965 1965    170
ZZXUAICE015 -98.8888 -998.8888 -998.8    NP15                           1966 1968    676
ZZXUAICE016 -98.8888 -998.8888 -998.8    NP16                           1968 1972   1356
ZZXUAICE017 -98.8888 -998.8888 -998.8    NP17                           1968 1969    482
ZZXUAICE019 -98.8888 -998.8888 -998.8    NP19                           1969 1973   1146
ZZXUAICE021 -98.8888 -998.8888 -998.8    NP21                           1972 1974    454
ZZXUAICE022 -98.8888 -998.8888 -998.8    NP22                           1974 1982   2862
ZZXUAICE026 -98.8888 -998.8888 -998.8    NP26                           1983 1986    824
ZZXUAICE028 -98.8888 -998.8888 -998.8    NP28                           1986 1988    915
ZZXUAICE030 -98.8888 -998.8888 -998.8    NP30                           1988 1990    576
ZZXUAICE031 -98.8888 -998.8888 -998.8    NP31                           1989 1991    717
"""