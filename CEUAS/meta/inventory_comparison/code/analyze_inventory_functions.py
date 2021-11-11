"""
    Module: analyze_inventory_functions.py
    
    Author:: Ambrogi Federico

"""
import os, sys
import pandas as pd
import numpy as np
import re 
import subprocess
import geopy.distance


pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)


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
            
            d = dms.split(' ')[0]
            m = dms.split(' ')[1]
            s = dms.split(' ')[2]
            
            dec = float(d) + float(m)/60 + float(s)/3600
            
            return sign*round(dec, 2)
        
        coord_dec = map(dms2dec, coord)        
        
        return coord_dec
    
    

class Analyze():
    
    def __init__(self, dist_limit = 30, inv = '', data = '' ):
        """ parameters ::
                     dist_limit : maximum distance allowed for station matching
                     inv :: dataframe of a single the inventory
                     data :: data extracted from a single file from one of the data sets """
        
        self.dist_limit = dist_limit
        self.inv = inv
        self.data = data
        
    def Distance(self, lat1,lon1,lat2,lon2):
        """ Return the Harvesine distance between two points """
        d = geopy.distance.geodesic( (lat1,lon1), (lat2,lon2) ).km
        return d 
        
    def AddDistance(self):
        """ Calculates the distances between the point and the stations in the inventory """

        distances = [ round( self.Distance( self.data.lats[0], self.data.lons[0] , self.inv['latitude'][i] , self.inv['longitude'][i] ), 1 ) for i in range(len(self.inv['latitude']) )     ] 
        self.inv['distances_km'] = distances
        return distances
        
    def Closest(self):
        """ Find all the stations that are within the distance limit,
              or have the same station identifier,
              """
        df = self.inv.sort_values(by =['distance_km'])
        return df
    
    def MatchingId(self):
       """ Check if there is a station id matching the input station id data """
       
       df = self.inv.loc[ self.inv['station_id'] == self.data.statid ] 
       return df
       
    
class Data():
    """ Main class to hold utilities for reading original input data from 
         - ERA5 1,2,1759,1761,3188 
         - IGRA2
         - NCAR
         - BUFR 
         """
    
    def __init__(self, era5_1= '', era5_2='', era5_1759 = '',  era5_1761 = '', era5_3188 = '', bufr='', ncar= '', igra2='',   ):
        
        self.era5_1 = era5_1
        self.era5_2 = era5_2
        self.era5_1759 = era5_1759
        self.era_1761= era5_1761        
        self.era5_3188 = era5_3188
        self.bufr = bufr
        self.igra2 = igra2
        self.ncar = ncar
        
    def read_odb(self, file=''):
        """ Extract data from a single ODB file 
        # odb headre -i file """
        
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
                        names=( 'station_id', 'latitude', 'longitude', 'elevation', 'dummy',
                                      'station_name', 'start', 'end', 'records') )
        igra2 = igra2[['station_id', 'station_name','latitude','longitude']]
        igra2.name = 'IGRA2'
        
        
        #############################################
        ### Reading WBAN data
        #############################################        
        wban =pd.read_fwf(self.wban,widths=(10,5,6,17,22,39,31,31,8,10,10,10,8,7),
                        names=('dum1','station_id','WMO_id','dum0','Country','dum2','station_name','dum3','From','To',
                                      'latitude','longitude','dum4','Elev'),
                        skiprows=1 )
        
        wban = wban[ [ 'station_id', 'WMO_id' , 'station_name', 'latitude', 'longitude' ] ]
        wban.name = 'WBAN'
        
        
        #############################################
        ### Reading OSCAR data
        #############################################
        oscar = pd.read_csv( self.oscar , sep='\t')
        oscar = oscar.rename(columns = {'StationId':'WIGOS', 'Longitude':'longitude', 
                                        'Latitude':'latitude' , 'StationName':'station_name',} )
        
        oscar = oscar [['WIGOS', 'station_name', 'latitude', 'longitude', 'station_name', 'ObsRems']]
        
        # converting lat, lon to decimal format 
        lat_dec = list(self.utils.degMinSec_to_decimal(self, oscar['latitude'] ))
        lon_dec = list(self.utils.degMinSec_to_decimal(self, oscar['longitude'] ))
        
        oscar['original_latitude'] = oscar['latitude']
        oscar['latitude'] = lat_dec
    
        oscar['original_longitude'] = oscar['longitude']
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

# initialize Data
data = Data(era5_1 = '/mnt/users/scratch/leo/scratch/era5/odbs/1' )

# Initialize inventory
inventory = Inventory( datadir ='../data/tables',  
                 oscar = "vola_legacy_report.txt" , 
                 igra2 = "igra2-station-list.txt"  , 
                 wban = "WBAN.TXT-2006jan"  , 
                 chuan = "Inventory_ERACLIM_upperair_2.1.txt",
                 utils = utils )

inventory.readInventory()

# This must be done automatically for the various files in the dataset 
f = data.era5_1 + '/era5.conv._82930'
data.read_odb(file = f)
d = data.read_odb(file = f)




################
### Analysis part 
################

analyze =Analyze(data = data, inv = inventory.inv['oscar'])
analyze.AddDistance()
analyze.MatchingId()

print(inv.chuan)


