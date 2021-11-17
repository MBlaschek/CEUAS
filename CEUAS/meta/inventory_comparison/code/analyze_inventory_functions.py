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


class Utils():
    """ Class holding a series of calculation utilites """
    def __init__(self):
        return 0

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
    
    def coord_consistency(self, lats='',lons='', threshold=''):
        """ Check if the absolute variation in time of the coordinates is below a certain threshold """
        try:
            lats, lons =  eval(self.data.lats[0]) ,  eval(self.data.lons[0])
        except:
            lats, lons =  [self.data.lats[0]] ,  [self.data.lons[0] ]
            
        check = abs( min(lons) - max(lons)  ) < 0.2 and abs( min(lats) - max(lats) ) < 0.2    
        return check
    
    
class Analyze():
    
    def __init__(self, dist_limit = 30, inv = '', data = '', utils=''):
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
        
        
    def Distance(self, lat1,lon1,lat2,lon2):
        """ Return the Harvesine distance between two points """
        d = geopy.distance.geodesic( (lat1,lon1), (lat2,lon2) ).km
        return d 
        
    
    def AddDistance(self):
        """ Calculates the distances between the point and the stations in the inventory """
        print("*** Adding distance ***")
        # checking if the values of lat and lons do not change significantly over time 
        
        try:
            lats, lons = eval(self.data.lats[0]),  eval(self.data.lons[0])
        except:
            lats, lons = [self.data.lats[0]] ,  [self.data.lons[0]]
        
        self.data.lats_list = lats
        self.data.lons_list = lons 
        
        check = self.utils.coord_consistency(self, lats=lats, lons=lons, threshold=0.2)        
        
        if check:
            # filter far off stations 
            df = self.inv
            red_inv = df.loc [  (abs(df['latitude'] - lats[-1]) < 2 ) & (abs(df['longitude'] - lons[-1]) < 2) ] 
            red_inv = red_inv.reset_index()
            
            distances = [ round( self.Distance( lats[-1], lons[-1] , red_inv['latitude'][i] , red_inv['longitude'][i] ), 1 ) for i in range(len(red_inv['latitude']) ) ]
            
            #remove unecessary list of lats and los, only keep largest and smallest
            
            
        else:
            print("Incosistent Latitude and Longitude values !!!")
            raise NotImplementedError()
            # TODO implement checking for coordinates check as in the harvester 
            
        red_inv['distance_km'] = distances
        red_inv = red_inv.sort_values(by =['distance_km'])
        red_inv.name = self.inv.name 
        
        self.red_inv = red_inv 
        
        return distances
        
    def Closest(self):
        """ Find all the stations that are within the distance limit,
              or have the same station identifier,
              """
        print("*** Find Closest ***")

        closest_indices = np.searchsorted( self.red_inv["distance_km"], self.dist_limit )
        inv_closest = self.red_inv[0:closest_indices]
        inv_closest.name = self.inv.name 
        inv_closest = inv_closest.reset_index(drop=True)
        inv_closest = inv_closest.drop(columns=['index'])
        
        self.found_inv[self.inv.name ]['matching_dist'] = inv_closest
    
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
            
        df = self.inv.loc[ self.inv[stats].isin( self.data.statid ) ] 
        inv_matchingId = df
        inv_matchingId.name = self.inv.name 
        inv_matchingId = inv_matchingId.reset_index(drop=True)        
        
        lats = self.data.lats_list
        lons = self.data.lons_list
        
        '''
        try:   
            lats, lons = eval(self.data.lats[0]),  eval(self.data.lons[0])
        except:
            lats, lons = [self.data.lats[0]],  [self.data.lons[0]]
            
        self.data.lats_list = lats
        self.data.lons_list = lons
        '''
        
        distances = [ round( self.Distance( lats[-1], lons[-1] , inv_matchingId['latitude'][i] , inv_matchingId['longitude'][i] ), 1 ) for i in range(len(inv_matchingId['latitude']) ) ]
        inv_matchingId['distance_km'] = distances
        
        self.found_inv[self.inv.name ]['matching_id'] = inv_matchingId
               
               
    def Combine(self):
        """ Concatenating the found df for each inventory and removing duplicated entries """
        id_match = self.found_inv[self.inv.name ]['matching_id']
        dist_match = self.found_inv[self.inv.name ]['matching_dist']
        
        if not ( id_match.empty and dist_match.empty ):
            res = pd.concat ( [id_match, dist_match], ignore_index=True) 
            res = res.drop_duplicates()
        elif ( id_match.empty and not dist_match.empty ):
            res = dist_match
        elif ( not id_match.empty and dist_match.empty ):
            res = id_match
            
        else:
            res = id_match # empty case
            for c in res.columns:
                res[c] = ['nan']            
            res['WIGOS_calc'] = ['nan'] 
            res['inventory'] = self.inv.name 
            self.best_match = res 
            return
                    
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
    
    def __init__(self, dataset=''):
        
        self.dataset = dataset
        
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
        
    def read_file(self, file):
        """ Wrapper for the functions to read input file from the datasets. 
             If a csv file exists holding the needed information, it is read.
             Otherwise it reads the original file with the appropriate functions. """
        
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
            self.lats = df['lats']
            self.lons = df['lons'] 
            self.min_date = df['min_date']
            self.max_date = df['max_date']
        
        else: # read from original files
            if 'era5_' in self.dataset :
                self.read_odb()
            elif 'igra' in self.dataset:
                self.read_igra2()
            elif 'ncar' in self.dataset:
                self.read_ncar()
            elif 'bufr' in self.dataset:
                self.read_bufr()
        
        
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
            print(b[i])
            
        # extract statid 
        comm = "odb sql 'select distinct statid'  -i FILE --no_alignment ".replace('FILE', file)   
        proc = subprocess.Popen(comm, stdout=subprocess.PIPE , shell=True)
        b = proc.stdout.read().decode('utf-8').split('\n')
        statids = [ eval(c) for c in np.unique( b[1:-1] ) ] 
        
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
                'statid':statids, 
                'lats': [lats] , 
                'lons': [lons], 
                'file': file , 
                'file_path': self.file , 
                'db': self.dataset  }
        
        # storing the 
        pd.DataFrame(dic).to_csv( self.out_dir  + '/' + self.file.split('/')[-1] + '.csv', sep = '\t' )        
        self.statid = statids
        self.lats = lats
        self.lons = lons 
        self.min_date = min_date
        self.max_date = max_date
        
    def read_bufr(self, file =''):
        """ Extract data from a single BUFR file """
        return 0    
    
    def read_igra2(self, file ='' ):
        """ Extract data from a single IGRA2 file (txt format)"""
        return 0
    
    def read_ncar(self, file=''):
        """ Extract data from a single ODB file """
        return 0    
    
    
    
class Inventory():
    """  Class holding the functionalities to read and homogenize 
    the OSCAR, IGRA2, WBAN and CHUAN inventories.
    For OSCAR and WBAN will convert lat and lon to decimal format. """
    
    def __init__(self, datadir ='', oscar ="", igra2= "", wban = '', chuan = '' , utils='' ):
        self.datadir = datadir + '/' + datadir
        self.oscar  = datadir + '/' + oscar
        self.igra2  = datadir + '/' + igra2
        self.wban  = datadir + '/' + wban
        self.chuan = datadir + '/' + chuan
        self.utils = utils 
        self.inv = {}
        
    def readInventory(self):
        """ Read the source inventories (OSCAR,IGRA2,WBAN,CHUAN)
        and convert them to pandas dataframe with uniform names for common columns"""
         
        #############################################
        ### Reading IGRA2 data
        #############################################        
        igra2 = pd.read_fwf(self.igra2, widths=(11, 9, 10, 7, 4, 30, 5, 5, 7),
                        names=( 'station_id_igra', 'latitude', 'longitude', 'elevation', 'dummy',
                                      'station_name', 'start', 'end', 'records') )
        
        # extracting standard 5 digits station id 
        stat_id = [ s[6:] for s in igra2['station_id_igra'] ]
        igra2['station_id'] = stat_id
        igra2 = igra2[['station_id_igra', 'station_name','latitude','longitude', 'station_id']]
        
        # clean for missing values in lat and lon (i.e. values < 0)
        igra2 = igra2.loc[ (igra2['latitude'] >= -180) & (igra2['longitude'] >= -180) ]
        igra2 = igra2.reset_index(drop=True)
        igra2.name = 'IGRA2'
        
        
        #############################################
        ### Reading WBAN data
        #############################################        
        wban =pd.read_fwf(self.wban,widths=(10,5,6,17,22,39,31,31,8,10,10,10,8,7),
                        names=('dum1','station_id','WMO_id','dum0','Country','dum2','station_name','dum3','From','To',
                                      'latitude','longitude','dum4','Elev'),
                        skiprows=1 )
        
        wban = wban[ [ 'station_id', 'WMO_id' , 'station_name', 'latitude', 'longitude' ] ]
        wban = wban.dropna(subset = ['latitude', "longitude"])
        #wban = wban.replace( {'00 00 00': 0 , '0  0  0':0} )
        
        # converting lat, lon to decimal format 
        lat_dec = list(self.utils.degMinSec_to_decimal(self, wban['latitude'] ))
        lon_dec = list(self.utils.degMinSec_to_decimal(self, wban['longitude'] ))
        
        wban['latitude'] = lat_dec
        wban['longitude'] = lon_dec

        wban['station_id'] = wban['station_id'].astype('int')
        
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
        lat_dec = list(self.utils.degMinSec_to_decimal(self, oscar['latitude'] ))
        lon_dec = list(self.utils.degMinSec_to_decimal(self, oscar['longitude'] ))
        
        oscar['original_lat'] = oscar['latitude']
        oscar['latitude'] = lat_dec
    
        oscar['original_lon'] = oscar['longitude']
        oscar['longitude'] = lon_dec        
        
        statids = [f.split('-')[-1] for f in oscar['WIGOS'] ]  # assume that th estation id is the last piece of the WIGOS id 
        oscar['station_id'] = statids
        
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
            chuan = chuan[['station_id', 'WMO_id' , 'station_name' , 'latitude' , 'longitude' ]]
            
            # filtering missing values for coordinates
            chuan = chuan.loc[ (chuan['latitude'] > -90) & (chuan['longitude'] > -90)]
            chuan = chuan.reset_index(drop=True)
            chuan['station_id']= chuan['station_id'].astype('int')

            chuan.name = 'CHUAN'
        
        # saving csv files
        chuan.to_csv('chuan.csv' , sep = '\t')
        igra2.to_csv('igra2.csv' , sep = '\t')
        oscar.to_csv('oscar.csv' , sep = '\t')
        wban.to_csv('wban.csv' , sep = '\t')
               
        # storing df as instances
        self.inv["oscar"] = oscar
        self.inv["igra2"] = igra2
        self.inv["wban"] = wban
        self.inv["chuan"] = chuan
        
        print(0)
        
        
    
    


# initialize Utils
utils = Utils

# Initialize inventory
inventory = Inventory( datadir ='../data/tables',  
                 oscar = "vola_legacy_report.txt" , 
                 igra2 = "igra2-station-list.txt"  , 
                 wban = "WBAN.TXT-2006jan"  , 
                 chuan = "Inventory_ERACLIM_upperair_2.1.txt",
                 utils = utils )

inventory.readInventory()

################
### Analysis part 
################

#f = 'era5.conv._82930'
def wrapper(data, file):
    """ Wrapper function to run all the analysis steps for a single file """
    
    print("Doing file: " , file )
    d = data.read_file(file=file)
    matching_inv = {}
    
    # checking and creating output directory 
    out_dir = 'inventories/' + data.dataset
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
        
    for i in inventories :
        print(' *** Analizing the inventory: ', i )
        analyze = Analyze( data = data, inv = inventory.inv[i], utils = utils)
        analyze.AddDistance()
        analyze.MatchingId()
        analyze.Closest()
        #matching_inv.append(analyze.found_inv)
        analyze.Combine()
        matching_inv[i] = analyze.best_match    
    
    l = [matching_inv[f] for f in  matching_inv.keys() ]
    
    all_inventories = pd.concat(l)
    
    # Add file information     
    all_inventories['lat_file'] = str( [max(data.lats_list), min(data.lats_list)] )
    all_inventories['lon_file'] = str( [max(data.lons_list), min(data.lons_list)]  )
    all_inventories['file'] = data.file.split('/')[-1]
    all_inventories['file_dataset'] = data.dataset
    all_inventories['file_min_date'] = data.min_date[0]
    all_inventories['file_max_date'] = data.max_date[0]
    all_inventories['file_max_date'] = data.max_date[0]
    all_inventories['file_statid'] = data.statid[0]
    
    name = out_dir + '/' + data.file.split('/')[-1] + '_inventories.csv'
    
    all_inventories.to_csv( name, sep = '\t', index = False)
    
    all_inventories_red = all_inventories[ [ 'station_id', 'station_name', 'latitude', 'longitude', 'original_lat', 'original_lon', 'distance_km', 'lat_file',
       'lon_file',  'inventory', 'WIGOS', 'WMO_id', 'WIGOS_calc' ] ]
    
    all_inventories_red.to_csv( name.replace('.csv','_reduced.csv'), sep = '\t' , index = False )
    
    print("Done :::" , data.file )
    
    
    

datasets = { 'era5_1': '/mnt/users/scratch/leo/scratch/era5/odbs/1' ,
                               'era5_2': '/mnt/users/scratch/leo/scratch/era5/odbs/2',
                               'era5_3188': '/mnt/users/scratch/leo/scratch/era5/odbs/3188',
                               'era5_1759': '/mnt/users/scratch/leo/scratch/era5/odbs/1759',
                               'era5_1761': '/mnt/users/scratch/leo/scratch/era5/odbs/1761',
                               
                               'bufr': 'mnt/users/scratch/leo/scratch/era5/odbs/',                                   
                               
                               'ncar': '',
                               'igra2': '',

                               } 
    


if __name__ == '__main__':

    databases = ['era5_2']
    inventories  = [ 'oscar', 'igra2', 'wban', 'chuan']
    
    POOL = True
    CHECK_MISSING = True
    
    # loop through each of the databases
    for db in databases:
        
        # getting only missing files, option CHECK_MISSING 
        if db == 'era5_1':
            flist=glob.glob(datasets[db] + "/era5.conv._*")
        elif db == 'era5_2':
            flist=glob.glob("/mnt/users/scratch/leo/scratch/era5/odbs/2/era5.conv._*")
            flist=[f for f in flist if '.gz' not in f]
            
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
        
        data = Data(dataset=db)
    
        if POOL: # running multiprocessing or single 
            p = Pool(30)
            func = partial(wrapper,data)
            out = p.map(func, flist)   
            
        else:
            for f in flist:
                w = wrapper(data, f)

