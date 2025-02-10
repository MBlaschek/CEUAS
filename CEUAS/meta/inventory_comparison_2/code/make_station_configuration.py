import os,sys
import pandas as pd
import urllib.request
import glob
from tqdm import tqdm
import numpy as np
from ast import literal_eval

from multiprocessing import Pool
from functools  import partial
import reverse_geocoder as reverse_geocoder 

# pd.set_option('display.max_columns', None)
# pd.set_option('display.max_rows', None)
# pd.set_option('display.width', None)
# pd.set_option('display.max_colwidth', -1)

import argparse


# the following are used to find the country given lat and lon (to fill territory information)
import requests
from shapely.geometry import mapping, shape
from shapely.prepared import prep
from shapely.geometry import Point
import pycountry 
import json


os.system('mkdir station_configuration')
    
# urls of the CDM station_configuration
# We need to use the table definition (not the table, which does not exist)
path='https://raw.githubusercontent.com/glamod/common_data_model/master/table_definitions/station_configuration.csv'
f=urllib.request.urlopen(path)
# reading column names and setting type to string, store in dictionary
col_names=pd.read_csv(f, delimiter='\t', quoting=3, nrows=0, comment='#')
columns_type_dic={col: str for col in col_names}
# need to repeat twice 
f=urllib.request.urlopen(path)

# reading station_configuration_columns
stat_conf=pd.read_csv(f, delimiter='\t', quoting=3, dtype=columns_type_dic ,na_filter=False,comment='#')
stat_conf_columns = stat_conf['element_name']

# mapping the column names from each inventory csv file, to the column names in the CDM
inventory_to_statConf = { 'WIGOS_best' : 'primary_id' ,
                'file_statid' : 'secondary_id',
                'station_name' : 'station_name',
                'latitude':'latitude',
                'longitude':'longitude',
                'file_min_date': 'start_date',
                'file_max_date':'end_date',
                'city':'city',
                'variables':'observed_variables'
                }

statConf_to_inventory = dict( zip( inventory_to_statConf.values(), inventory_to_statConf.keys()) )

### supplementary columns to fill in CUON
metadata_contact = 'L.Haimberger'
metadata_contact_role = '0'

### to decide 
primary_id_scheme = 0  # see Fabios list 
secondary_id_scheme = 3 # see Fabios list 
station_crs = 0
station_type = 1

platform_type = 0
platform_sub_type = 63

### to be copied when available from GRUAN stations 
operating_institute = 'NMS'
telecommunication_method = '' 
station_automation = '' 
measuring_system_model = '' 
measuring_system_id = ''

### save above variables in dictionary
special_columns = {'primary_id_scheme':primary_id_scheme,
                   'secondary_id_scheme':secondary_id_scheme,
                   'station_crs':station_crs,
                   'station_type':station_type,
                   'platform_type':platform_type,
                   'platform_sub_type':platform_sub_type, 
                   'operating_institute':operating_institute,
                   #'operating_territory':operating_territory,
                   'telecommunication_method':telecommunication_method,
                   'station_automation':station_automation,
                   'measuring_system_model':measuring_system_model,
                   'measuring_system_id':measuring_system_id, 
                   'metadata_contact_role':metadata_contact_role,
                   'metadata_contact':metadata_contact
                   }

# initialize empty dic as the stat_conf table 

def extract_inventory(ids):
    """ Converts primary_id to inventory list by checking the pseudo wigos ids
            self.wigos = { 'IGRA2'    :'0-20300-0-',
                       'CHUAN'    :'0-20400-0-',                       
                       'WBAN'     :'0-20500-0-',           
                       'SCHROEDER':'0-20600-0-',                    
                       'WMO'      :'0-20700-0-',
                       
                       'AMMA' : '0-20800-0-',
                       'GIUB' : '0-20900-0-' ,
                       
                       'HARA' : '0-20111-0-' ,
        }
    """ 
    
    invs = []
    for i in ids:
        i = str(i)
        if '-20300-' in i:
            inv = 'IGRA2'
        elif '-20400-' in i:
            inv = 'CHUAN'
        elif '-20500-' in i:
            inv = 'WBAN'
        elif '-20600-' in i:
            inv = 'SCHROEDER'
        elif '-20700-' in i:
            inv = 'WMO'            
        elif '-20800-' in i:
            inv = 'AMMA'
        elif '-20111-' in i:
            inv = 'HARA'
        else:
            inv = 'OSCAR'
        invs.append(inv)
        
    return invs


def iso_from_country_from_coordinates( lat='', lon='', json_data='',  cdm_sub_regions = '', file = ''):
    """ Extract the iso code / country from given lat and lon.
    Taken from: https://stackoverflow.com/questions/20169467/how-to-convert-from-longitude-and-latitude-to-country-or-city.
    
    Converts then the code to the CDM table convention 
    https://github.com/glamod/common_data_model/blob/master/tables/sub_region.dat
    
    So first check the country with the script, 
    then convert the proper name according to the CDM table sub_region.dat
    """

    countries = {}
    for feature in json_data["features"]:
        geom = feature["geometry"]
        country = feature["properties"]["ADMIN"]
        countries[country] = prep(shape(geom))
    
    #print(len(countries))
    def get_country(lon, lat):
        point = Point(lon, lat)
        for country, geom in countries.items():
            if geom.contains(point):
                return country
    
        return "unknown"
    
    retrieve_country = get_country(lon, lat)
    
    if retrieve_country == 'unknown':
            g = reverse_geocoder.search( (lat, lon) )
            alpha_3 = pycountry.countries.get(alpha_2= g[0]['cc']).alpha_3
            if alpha_3 == "WLF":
                alpha_3 = "WSM" # Wallis and Futuna  
            territory = cdm_sub_regions[cdm_sub_regions.alpha_3_code == alpha_3 ].sub_region.values[0]
            
            return territory
        
    else:
        if retrieve_country == 'Russia':
            retrieve_country = "Russian Federation"
        if retrieve_country == 'Iran':
            retrieve_country = 'Iran, Islamic Republic of'
        if retrieve_country == 'United Republic of Tanzania':            
            retrieve_country = 'Tanzania, United Republic of'
        if retrieve_country in  ['Somaliland']: # countries with geopolitical uncertainties 
            return ''
        if retrieve_country in ['Democratic Republic of Congo', 'Republic of Congo', 'Democratic Republic of the Congo']: # countries with geopolitical uncertainties 
            retrieve_country = 'Congo, The Democratic Republic of the'
        if retrieve_country == 'Syria':            
            retrieve_country = 'Syrian Arab Republic'        
        if retrieve_country == "Moldova":
            retrieve_country = "Moldova, Republic of"
        if retrieve_country == 'Guinea Bissau':            
            retrieve_country = 'Guinea-Bissau'      
        if retrieve_country == 'Macedonia':            
            retrieve_country = 'North Macedonia'                  
        if retrieve_country == 'Czech Republic':            
            retrieve_country = 'Czechia'    
        if retrieve_country == 'Republic of Serbia':            
            retrieve_country = 'Serbia'    
        if retrieve_country == 'Venezuela':            
            retrieve_country = 'Venezuela, Bolivarian Republic of'
        if retrieve_country == 'United States of America':
            retrieve_country = 'United States'
        if retrieve_country == 'Vietnam':
            retrieve_country = 'Viet Nam'              
        if retrieve_country == 'Ivory Coast':
            retrieve_country = "Côte d'Ivoire"        
        if retrieve_country == 'Bolivia':
            retrieve_country = "Bolivia, Plurinational State of"            
        """
        if retrieve_country ==    "French Southern and Antarctic Lands":
            retrieve_country = "France"      
        if retrieve_country ==    "North Korea":
                retrieve_country = "Korea, Democratic People's Republic of"                   
        if retrieve_country ==  "Taiwan":
            retrieve_country = "Taiwan, Province of China"                 
        if "US " in retrieve_country:
            retrieve_country = "United States"         
        if retrieve_country == "Brunei":
            retrieve_country = "Brunei Darussalam"      
        if retrieve_country == "United States Virgin Islands":
            retrieve_country = "Virgin Islands, U.S."   
        if retrieve_country == "Laos":
            retrieve_country = "Lao People's Democratic Republic"     
        """
        
        if "US " in retrieve_country:
            retrieve_country = "United States"      
            
        # ------------   
        # TODO replace statements above with the dictionary below, to keep same shape 
        terr_names = { "Aland": "Åland Islands",
                 "South Korea":"Korea, Republic of" ,
                 "Falkland Islands":"Falkland Islands (Malvinas)",
                 "Macao S.A.R": "Macao",
                 "Federated States of Micronesia":"Micronesia, Federated States of",
                 "East Timor": "Timor-Leste",
                 "Indian Ocean Territories":"British Indian Ocean Territory",
                 "Akrotiri Sovereign Base Area":"United Kingdom",
                 "Hong Kong S.A.R.": "Hong Kong",
                 "Saint Barthelemy":"Saint Barthélemy",
                 "Saint Helena": "Saint Helena, Ascension and Tristan da Cunha",
                 "Cape Verde": "Cabo Verde",
                 "The Bahamas": "Bahamas",
                 "Laos": "Lao People's Democratic Republic",
                 "United States Virgin Islands":"Virgin Islands, U.S.",
                 "Brunei":"Brunei Darussalam",
                 "Taiwan":"Taiwan, Province of China"    ,
                 "North Korea":"Korea, Democratic People's Republic of" ,
                 "French Southern and Antarctic Lands":"France"  ,
                 
        }
        
        if retrieve_country in terr_names.keys():
            retrieve_country = terr_names[retrieve_country]

            
        # Aland, South Korea, Falkland Islands, Macao S.A.R, Federated States of Micronesia, East Timor, Indian Ocean Territories, Akrotiri Sovereign Base Area
        # Palestine, Hong Kong S.A.R. , Saint Barthelemy,  Saint Helena, Cape Verde, The Bahamas, 
        if retrieve_country == "Palestine":
            alpha_3 = "PSE" # taken directly from the CDM codes     
        elif retrieve_country == "Kosovo":
            return ''
        elif retrieve_country == 'Dhekelia Sovereign Base Area': 
            alpha_3 = 'CYP'
        else:
            alpha_3 = pycountry.countries.get(name=retrieve_country).alpha_3
        if alpha_3 == "UZB":
            alpha_3 = "USB"
            
        try:
            territory = cdm_sub_regions[cdm_sub_regions.alpha_3_code == alpha_3 ].sub_region.values[0]
        except:
            territory = ''
            
    #print(territory, file )
    a = open('inventories/territory_code.csv', 'a+')
    a.write(file.split('/')[-1].split('_inventories')[0] + '\t' + str(territory) + '\n')
    a.close()
    
    return territory
    
def get_best_inventory(df):
    """ Select best inventory among the ones available.
    Give preference to upper-air / radiosonde station from OSCAR, then generic OSCAR,
    then IGRA2, then WBAN and finaly CHUAN inventory """
        
    df_best = df.loc[df.WIGOS_best == df.WIGOS_calc] # here, OSCAR thes best matching inventory
    if not df_best.empty:
        return df_best[:1]
    
    else: 
        for i in ['OSCAR', 'IGRA2', 'WBAN' , 'CHUAN', 'WMO', 'SCHROEDER', 'AMMA', 'HARA']:
            dfr = df[ df.inventory == i ]        
            #dfr = dfr.dropna( subset = ['latitude','longitude'])  
            if dfr.empty:
                continue
            if i=="OSCAR":         
                dfr_radio = dfr[dfr.isRadio == True ] # check if upperAir/radiosonde station flag 
                
                # select OSCAR radio station,
                # otherwise generic OSCAR station
                
                if not dfr_radio.empty:
                    return dfr_radio
                else:
                    return dfr[0:1]
                
            else:
                dfr = dfr.sort_values(by=['distance_km'] )     
                return dfr[0:1]  
            
    return dfr 

# directory containing each single inventory file for each station 
    
def make_inventory(v):
    """ Creates the inventory for the given database """
    
    print(' *** Creating the inventory for ' , v )
    if not os.path.isdir('station_configuration/logs'):
        os.system('mkdir station_configuration/logs')
        
    inv_path = 'inventories/' 

    files = glob.glob(inv_path + '/' + v + '/' + '*inventories_iden*')
    files = [f for f in files if 'all' not in f and 'Ids' not in f  ]
    
    # initialize each inventory dic
    stat_conf_dic = {}
    for c in stat_conf_columns:
        stat_conf_dic[c] = []

    cdm_sub_regions = pd.read_csv("../data/sub_region.dat", sep = '\t')
    # Retrieve or read json files with country borders
    if not os.path.isfile("../data/countries.geojson"):
        json_data = requests.get("https://raw.githubusercontent.com/datasets/geo-countries/master/data/countries.geojson").json()
    else:
        a = open('../data/countries.geojson')
        json_data = json.load(a)
    
    # holding extra variables that do not appear in the stat_conf but can be useful
    # will be saved in additional station_conf file 
    
    extra_vars = { 'distance_km': [],
                            'distance_km_minusLat': [],
                            'city_dist_km': [] ,
                            'elevation':[],
                            'file': [] }
    
    if os.path.isfile('inventories/territory_code.npy'):
        terr_dic = np.load('inventories/territory_code.npy' , allow_pickle=True ).item()
        
    else:
        terr_dic = {'file':[], 'terr':[]}
    
    #files = [f for f in files if '10393' in f or '10395' in f or '09393' in f ]
    #files = files[4000:]
    for file in tqdm(files):    
        if '26063' in file:
            print(file)
        if '26072' in file:
            print(file)
            
        df = pd.read_csv(file, sep = '\t', header = 0)
        df = df.dropna( subset = ['latitude','longitude'])
        df = df.loc[df['distance_km'] < 30 ]
        
        if df.empty:
            a = open('station_configuration/logs/unidentified_' + v + '.txt' , 'a+' )
            a.write(file + '\n')
            continue
        
        # select best inventory available        
        best_inv = get_best_inventory(df)
        #if best_inv.empty:
        #    continue
        
        for e in extra_vars.keys():
                if e != "elevation":
                    extra_vars[e].append(best_inv[e].values[0])                    
                else:
                    try:
                        df['elevation'] = df['elevation'].astype(float)
                    except:
                        df['elevation'] =  [str(i).replace(',', '.') for i in df.elevation.values]
                        df['elevation'] = df['elevation'].astype(float)
                        
                    #a = df.sort_values(by=['distance_km']) # already done above ?
                    #elev = [i.replace(',', '.') for i in df.elevation.values]
                    elev = [e for e in df.elevation.values]
                    
                    value = [e for e in elev if e and e > -999 ]
                    if len(value) >0:
                        extra_vars[e].append(value[0])
                    else:
                        extra_vars[e].append(-999) 


        for c in stat_conf_columns:
            # print(c, len(stat_conf_dic[c]))
            
            if c in special_columns.keys():  # specific values for columns 
                stat_conf_dic[c].append(special_columns[c])            

            elif c == "operating_territory":
                ### TO DO !!! remove comment 
                fname = file.split('/')[-1].split('_inventories')[0]  
                
                try: # try to use lookup table
                    terr = terr_dic[ fname ]
                    #print(' Extracted terrytory' , fname )
                except:
                    # this is slow due to the implementation of the function below
                    print('Calculating territory --- ')
                    terr = iso_from_country_from_coordinates(lat=best_inv.latitude.values[0], lon=best_inv.longitude.values[0] , 
                                                          cdm_sub_regions=cdm_sub_regions,
                                                          json_data = json_data, file = file)      
                    
                    terr_dic[  fname ] = terr
                    
                stat_conf_dic[c].append(terr) 
                
            elif c in statConf_to_inventory.keys():
                if c == 'secondary_id':
                    a = best_inv[ statConf_to_inventory[c] ].values[0]
                    #print(0)
                    
                elif c == 'start_date' or c == 'end_date':
                    a = best_inv[ statConf_to_inventory[c] ].values[0]
                    if type(a) == str:
                        if '-' in a:
                            a = ''.join(a.split('-'))
                    try:
                        b = float(a)
                    except:
                        #print(0)
                        continue
                    
                try:
                    value = best_inv[ statConf_to_inventory[c] ].values[0]
                    if str(value)=='999' or value==999:
                        stat_conf_dic[c].append('')
                        if c == 'start_date' or c == 'end_date':
                            print('check' ,   file )
                    else:
                        if 'city' in c:
                            stat_conf_dic[c].append( value.title() )
                        else:
                            stat_conf_dic[c].append(value)
                except:
                    print('wrong')
                    stat_conf_dic[c].append('')
            else:
                stat_conf_dic[c].append('')
                
    statconf_df = pd.DataFrame(stat_conf_dic)
    ind = statconf_df.index
    # adding record number
    statconf_df['record_number'] = list(range(1,len(statconf_df)+1))
    statconf_df.to_csv('station_configuration/' + v + '_station_configuration.csv' , sep= '\t' )           
    statconf_df.to_excel('station_configuration/' + v + '_station_configuration.xlsx' )
        

    for c in extra_vars.keys():
        statconf_df[c] = extra_vars[c]
    
    # adding the inventory name as source of the pseudo WIGOS id
    inv = extract_inventory(statconf_df.primary_id.values)
    statconf_df['inventory'] = inv     



    statconf_df.to_csv('station_configuration/' + v + '_station_configuration_extended.csv' , sep= '\t' )           
    statconf_df.to_excel('station_configuration/' + v + '_station_configuration_extended.xlsx' )

    np.save( 'inventories/territory_code', terr_dic, allow_pickle=True )
    print('*** Done with the inventory ',  v )








def make_CUON():
    ''' Create one single inventory out of the CUON datasets '''
    
    inv = []
    for s in glob.glob('station_configuration/*_station_configuration_extended*'):
        if 'CUON' in s or 'xl' in s or 'all'  in s:           
            continue
        df = pd.read_csv(s, sep='\t')
        print(s , ' ' , len(df))
        inv.append(df)
        
    # first I clean and extract the unique primary ids, then I will loop through them to extrac the combined information from them all 
    inventory_cuon = pd.concat(inv)
    inventory_cuon = inventory_cuon[ [ c for c in inventory_cuon.columns if c != 'record_number'  ] ] # remove record_number column 
    inventory_cuon = inventory_cuon.dropna(subset = ['primary_id']) # drop missing primary ids
    inventory_cuon = inventory_cuon.reset_index()
    
    all_ids = np.unique(inventory_cuon['primary_id'].astype(str))
    
    # storing combined CUON dictionary 
    combined_cuon= {}
    # reading all station configuration columns
    station_conf_columns = pd.read_csv('station_configuration/era5_1_station_configuration.csv', sep='\t').columns
    
    for c in station_conf_columns:
        if c in  'index' or 'Unnamed' in c :
            continue     
        combined_cuon[c] = []
    
    combined_cuon['elevation'] = []
    # looping through all the unique ids and combining data
    
    #all_ids = ['0-20001-0-10393']
    for index, stat in tqdm(enumerate(all_ids)):
        indices = np.where( inventory_cuon['primary_id'].astype(str) == stat )[0]
        df_i = inventory_cuon.iloc[indices]
    
        for c in list(station_conf_columns):
            if c in  'index' or 'Unnamed' in c :
                continue 

            if c == 'record_number':
                combined_cuon[c].append(index)
                index = index+1 
                continue
            
            values = [f.replace('[','').replace(']','') for f in df_i[c].astype(str).values if isinstance(f,str) ]
            #if c == 'observed_variables':
            #   values = np.unique( [f.split(',') for f in values ] ) 
            cleaned = []
                
            # extracting min and max date 
            if c == 'start_date':
                values = [v.replace('-','') for v in values ]
                try:
                    values = [ int ( min ( [float(s) for s in values if float(s) > 2021 ] ) ) ]
                except:
                    values = [ int ( min ( [float(s) for s in values ] ) ) ]
                    
            if c == 'end_date':
                values = [v.replace('-','') for v in values ]                
                try:
                    values = [ int ( max ( [float(s) for s in values ]) ) ]
                except:                    
                    values = [ int ( max ( [float(s) for s in values ] ) ) ]            
            #if c == 'observed_variables':
            #    print(0)
            for v in values:
                if c in ['observed_variables']:
                    lista = v.split(',')
                    for val in lista:
                        val = val.replace(' ', '')
                        try:
                            vv = str(literal_eval(val))
                        except:
                            vv =  str(val)          
                            
                        # converting old numberinf scheme for observed variables to new ones form the CDM (i.e. using profile observations)
                        if vv == '38': # relative humidity
                            vv = '138'
                        if vv == '85': # temperature
                            vv = '126'          
                        if vv == '104':  # eastward wind
                            vv = '139'  
                        if vv == '105': # northward wind
                                vv = '140'  
                        if vv == '36': # dew point
                            vv = '137'  
                                        
                        if vv not in cleaned:
                            cleaned.append( vv )
                            
                else:
                    try:
                        v =  str(literal_eval(v))
                    except:
                        pass
                    
                    if v is not None and v not in cleaned:
                        cleaned.append(str(v))
                    
            cleaned.sort()
           
            values = [f for f in cleaned if f is not None ]
    
            if isinstance(values, list):
                values = ','.join(values)
            elif isinstance(values, str):
                pass
            if values == 'nan':
                values = ''
                
            combined_cuon[c].append(values)
            #print(0)
        elev = str(df_i['elevation'].values[0]).replace(',','.')
        
        combined_cuon['elevation'].append(elev)
        
    d = pd.DataFrame(combined_cuon)
    combined_cuon['record_number'] = list(range(len(    combined_cuon['latitude'] ) ))
    d.to_csv('station_configuration/CUON_station_configuration_extended.csv', sep='\t')
    d.to_excel('station_configuration/CUON_station_configuration.xlsx' )


def make_ORPHAN():
    ''' Create combined inventory for all orphans, mobile and problematic stations  '''
    
    inv = []
    for s in glob.glob('station_configuration/*orphans_station_configuration_extended*'):
        if 'CUON' in s or 'xl' in s or 'all'  in s:
            continue
        df = pd.read_csv(s, sep='\t')
        print(s , ' ' , len(df))
        inv.append(df)
        
    # first I clean and extract the unique primary ids, then I will loop through them to extrac the combined information from them all 
    inventory_cuon = pd.concat(inv)
    inventory_cuon = inventory_cuon[ [ c for c in inventory_cuon.columns if c != 'record_number'  ] ] # remove record_number column 
    inventory_cuon = inventory_cuon.dropna(subset = ['primary_id']) # drop missing primary ids
    inventory_cuon = inventory_cuon.reset_index()
    
    all_ids = np.unique(inventory_cuon['primary_id'].astype(str))
    
    # storing combined CUON dictionary 
    combined_cuon= {}
    # reading all station configuration columns
    station_conf_columns = pd.read_csv('station_configuration/era5_1_orphans_station_configuration_extended.csv', sep='\t').columns
    
    for c in station_conf_columns:
        if c in  'index' or 'Unnamed' in c :
            continue     
        combined_cuon[c] = []
    
    combined_cuon['elevation'] = []
    # looping through all the unique ids and combining data
    
    #all_ids = ['0-20001-0-10393']
    for index, stat in tqdm(enumerate(all_ids)):
        indices = np.where( inventory_cuon['primary_id'].astype(str) == stat )[0]
        df_i = inventory_cuon.iloc[indices]
    
        for c in list(station_conf_columns):
            if c in  'index' or 'Unnamed' in c :
                continue 

            if c == 'record_number':
                combined_cuon[c].append(index)
                index = index+1 
                continue
            
            values = [f.replace('[','').replace(']','') for f in df_i[c].astype(str).values if isinstance(f,str) ]
            #if c == 'observed_variables':
            #   values = np.unique( [f.split(',') for f in values ] ) 
            cleaned = []
                
            # extracting min and max date 
            if c == 'start_date':
                values = [v.replace('-','') for v in values ]
                try:
                    values = [ int ( min ( [float(s) for s in values if float(s) > 2021 ] ) ) ]
                except:
                    values = [ int ( min ( [float(s) for s in values ] ) ) ]
                    
            if c == 'end_date':
                values = [v.replace('-','') for v in values ]                
                try:
                    values = [ int ( max ( [float(s) for s in values ]) ) ]
                except:                    
                    values = [ int ( max ( [float(s) for s in values ] ) ) ]            
            #if c == 'observed_variables':
            #    print(0)
            for v in values:
                if c in ['observed_variables']:
                    lista = v.split(',')
                    for val in lista:
                        val = val.replace(' ', '')
                        try:
                            vv = str(literal_eval(val))
                        except:
                            vv =  str(val)          
                            
                        # converting old numberinf scheme for observed variables to new ones form the CDM (i.e. using profile observations)
                        if vv == '38': # relative humidity
                            vv = '138'
                        if vv == '85': # temperature
                            vv = '126'          
                        if vv == '104':  # eastward wind
                            vv = '139'  
                        if vv == '105': # northward wind
                                vv = '140'  
                        if vv == '36': # dew point
                            vv = '137'  
                                        
                        if vv not in cleaned:
                            cleaned.append( vv )
                            
                else:
                    try:
                        v =  str(literal_eval(v))
                    except:
                        pass
                    
                    if v is not None and v not in cleaned:
                        cleaned.append(str(v))
                    
            cleaned.sort()
           
            values = [f for f in cleaned if f is not None ]
    
            if isinstance(values, list):
                values = ','.join(values)
            elif isinstance(values, str):
                pass
            if values == 'nan':
                values = ''
                
            combined_cuon[c].append(values)
            #print(0)
        
                    
        combined_cuon['elevation'].append(np.nan)
        
    d = pd.DataFrame(combined_cuon)
    combined_cuon['record_number'] = list(range(len(    combined_cuon['latitude'] ) ))
    d.to_csv('station_configuration/CUON_orphans_station_configuration_extended.csv', sep='\t')
    d.to_excel('station_configuration/CUON_orphans_station_configuration.xlsx' )



def make_inventory_orphans_mobile(v):
    """ Special case of mobile or oprhans station, where we will collect only the id given in the original file and extract the set of changing lat and lon,
    since it is not possible to identify an actual station; all other fields in the station_inventory will be neglected """

    print(' *** Creating the inventory for Dataset ' , v )

            
    inv_path = 'inventories/' 
    
    files = glob.glob(inv_path + '/' + v + '/' + '*')
    files = [f for f in files if 'inventories_iden' not in f and 'reduced' not in f and 'logs' not in f ]
    
    if v in ['era5_1' , 'era5_2']:
        files_m = glob.glob(inv_path + '/' + v + '_mobile' + '/' + '*')
        files_m = [f for f in files_m if 'inventories_iden' not in f and 'reduced' not in f and 'logs' not in f ]
        files.extend(files_m)
    
    
    # dic holding data
    stat_conf_dic = {}
    stat_conf_dic['file'] = []
    stat_conf_dic['kind'] = []
    
    for file in tqdm(files):    
        
        if '88888' in file:  # wrong file 
            continue 
        df = pd.read_csv(file, sep = '\t', header = 0)


        ### 
        if '0-20999-0' in df.WIGOS_best.values[0]:
            kind = 'MOBILE'
        elif '0-20888-0'in df.WIGOS_best.values[0]:
            kind = 'ORPHAN'
        elif '0-20777-0'in df.WIGOS_best.values[0]:
            kind = 'MATCHING_ID_NONMATCHING_COORD'
        elif '0-20666-0'in df.WIGOS_best.values[0]:
            kind = 'COORDINATE_ISSUES'
            
    
        # restric the columns, all the others will be empty
        stat_conf_columns_red =  ['primary_id' , 'secondary_id' , 'longitude' , 'latitude' , 'start_date' , 'end_date', 'observed_variables'] 
        # ['primary_id' , 'secondary_id' , 'longitude' , 'latitude' , 'start_date' , 'end_date' ]
        for c in stat_conf_columns:
            if c not in stat_conf_dic.keys():
                stat_conf_dic[c] = []
            if c not in stat_conf_columns_red:
                stat_conf_dic[c].append(' ')    
            
        stat_conf_dic['primary_id'].append (str(df.WIGOS_best.values[0]) ) 
        stat_conf_dic['secondary_id'].append( str(df.file_statid.values[0]).replace('\n','')  )
        
        stat_conf_dic['start_date'].append( str(df['file_min_date'].values[0]) )
        stat_conf_dic['end_date'].append( str(df['file_max_date'].values[0]) )
        try:
            
            stat_conf_dic['latitude'].append( str(df['lat_file'].values[0]) )
            stat_conf_dic['longitude'].append( str(df['lon_file'].values[0]) )
        except:
            stat_conf_dic['latitude'].append( np.nan )
            stat_conf_dic['longitude'].append( np.nan  )            
            
        stat_conf_dic['observed_variables'].append( str(df['variables'].values[0]) )

        stat_conf_dic['kind'].append( kind )
        stat_conf_dic['file'].append( str(df['file'].values[0]))
    statconf_df = pd.DataFrame(stat_conf_dic)
    # adding record number
    statconf_df['record_number'] = list(range(1,len(statconf_df)+1))
    statconf_df.to_csv('station_configuration/' + v + '_orphans_station_configuration_extended.csv' , sep= '\t' )           
    statconf_df.to_excel('station_configuration/' + v + '_orphans_station_configuration_extended.xlsx' )
        
        
    


def merge_inventories():
    
    inventory_cuon = pd.read_csv('station_configuration/CUON_station_configuration.csv', sep='\t')
    ids_cuon = np.unique(inventory_cuon['primary_id'].astype(str))
        
    # read IGRA, GRUAN
    igra_f = pd.read_excel('Station_configuration_HARM.xlsx', sheet_name='IGRA')
    igra_f = igra_f[ [ c for c in igra_f.columns if c != 'record_number'  ] ]
     
    gruan_f = pd.read_excel('Station_configuration_HARM.xlsx', sheet_name='GRUAN')
    gruan_f = gruan_f[ [ c for c in gruan_f.columns if c != 'record_number'  ] ]

    # fix a bug in gruan inventory, one station has a blank space in its ids so it is not located
    g_clean = [ s.replace(' ','').replace('1-620-2001-0507','1-620-2001-08507') for s in gruan_f.primary_id ]
    # GRACIOSA has two entries in the OSCAR both valid for upper-air but different WIGOS
    
    gruan_f['primary_id'] = g_clean 
    
    #ids_igra_f = igra_f['primary_id']
    
    # combining all the inventories
    combining_all = {}
    for c in inventory_cuon.columns:
        if c in  ['index','Unnamed: 0','record_index','record_number']:
            continue     
        combining_all[c] = []
        
    #inventory_cuon = inventory_cuon[:100]
    for index, row in tqdm( inventory_cuon.iterrows() ):     
        pi = row.primary_id 
        
        igra = igra_f.loc[igra_f['primary_id'] == pi ]
        igra = igra.reset_index()
        gruan = gruan_f.loc[gruan_f['primary_id'] == pi ]
        gruan = gruan.reset_index()
        
        for c in inventory_cuon.columns:
            if c in ['index', 'Unnamed: 0', 'record_index','record_number']:
                continue
                
            # renumbering and adding observed variables common to all datasets 
                            
            #if c == 'primary_id' and row[c] == '0-20000-0-63740':
            #    print(0)
            if len(igra) ==0 and len(gruan) ==0:
                s = row[c]
                if s == 'None':
                    s = ''
             
            else:
                s = row[c] 
                if s == 'None':
                    s = ''
    
                if c in ['secondary_id', 'metadata_contact', 'station_name', 'city']: # summing alls
                    
                    if len(igra) >0:
                        try:
                            if s !=  igra[c].values[0] and s != 'None':
                                s = s + ',' +  igra[c].values[0]
                        except:
                            s = s 
                    if len(gruan) >0:
                        try:
                            if s != gruan[c].values[0]:
                                s = s + ',' +  gruan[c].values[0]  # city is missing in GRUAN
                        except:
                            pass
    
                elif (c in special_columns.keys() or c in ['operating_territory']):
                    if len(gruan) >0:
                        s =  gruan[c].values[0]
                    elif len(igra) >0:
                        s =  igra[c].values[0]
                    else:
                        s = row[c]
                else:
                    s = s
                    
            s = str(s)
            if s != '' and s.split(',')[0] == '':
                s = s.split(',')[1]
            if s == 'nan':
                s = ''
            if c in ['secondary_id', 'station_name', 'city']:
                lista =  s.split(',') 
                lista = [l for l in lista if l != 'None' ]
                s = ','.join( list( np.unique( lista ) ) )
            if c == 'start_date': # taken from CUON 
                #lista =  s.split(',') 
                #lista = [l for l in lista if l != 'nan' ]        
                #s = min(lista)
                pass
            if c == 'end_date':
                #lista =  s.split(',') 
                #lista = [l for l in lista if l != 'nan' ]        
                #s = max(lista)
                pass
            if c in ['latitude', 'longitude']:
                if len(s) >8:
                    f = s.split(',')
                    s = f[-1]
            if c in ['metadata_contact']:
                lista = s.split(',')
                lista = list( np.unique( lista ) )
                              
                if len(lista) >1:
                    s = ','.join( list( np.unique( lista ) ) )
                else:
                    s =lista[0]
                    
            if c in ['observed_variables']: 
                s = s + ',57' # add pressure as variable, see CDM table observed_variables
                if len(gruan) >0:                        
                    for o in ['137','138']:
                        if o not in s.split(','):
                            s = s + ',' + o            
                            
            if s == 'nan':
                s = ''      
                    
            combining_all[c].append(s)     
    
    
    """ Add station only found in GRUAN """
    ids_gruan = gruan_f['primary_id'].values
    ids_cuon = inventory_cuon['primary_id'].values
    only_gruan = [s for s in ids_gruan if s not in ids_cuon]
    
    for station in only_gruan:     
        df = gruan_f.loc[gruan_f['primary_id'] == station ]
        for c in combining_all.keys():
            if c in ['index', 'Unnamed: 0', 'record_index']:
                continue
            elif c in ['elevation']: # remember, elevation is not a proper field in the station configuration table, I added it artifically to my station configs
                s = np.nan 
            elif c in ['start_date', 'end_date']:
                s = df[c].dt.strftime('%Y%m%d').values[0]
            else:
                s = df[c].values[0]
            
            if c == 'primary_id':
                s = s.replace(' ','')
            combining_all[c].append(s)
        
    # Saving the complete dataframe
    for c in combining_all.keys():
        print(c, len(combining_all[c] ) )
        
    # check record_number, elevation 
    df = pd.DataFrame(combining_all)    
    df = df.sort_values( by = 'primary_id' )  # sorting by primary_id
    
    
    df.insert(2, 'record_number', list(range(len(df))) )
    df.to_csv('station_configuration/all_combined_station_configuration.csv' , sep='\t')
    df.to_excel('station_configuration/all_combined_station_configuration.xlsx' )
    
    
    print('--- Finished with the global inventory --- ')
    

# define a list of operation to perform between  [ 'INVENTORY', CUON', 'MERGE']
TODO = ['INVENTORY', 'CUON', "MERGE"]
POOL = True   # NB when running the first time, it might be necessary tot urn off the POOL due to the function that extracts the territory which has another POOL inside
n_pool = 40

parser = argparse.ArgumentParser(description="Crete station configuration table")

parser.add_argument('--dataset' , '-d', 
                      help="Select the dataset"  ,
                      type = str,
                      default = 'era5_1' )

args = parser.parse_args()
v                = args.dataset



if v in  [ 'era5_2', 'era5_1759', 'era5_1761', 'era5_3188', 'bufr', 'maestro', 'ncar', 'igra2', 'era5_1', 'amma', 'hara', 'giub']:
    WHAT = 'INVENTORY'
    
elif v in ['CUON', 'MERGE', 'ORPHAN']:
    WHAT = v
    
elif v in ['npsound' , 'era5_1_mobile' , 'era5_2_mobile', 'shipsound']:
    WHAT = 'MOBILE'
    

if WHAT == 'INVENTORY':
        # inventories to process
        #inventories_all = [ 'era5_2', 'era5_1759', 'era5_1761', 'era5_3188', 'bufr', 'ncar', 'igra2', 'era5_1']
        #inventories = ['era5_1']

        inventories =  [v]
        
        if not POOL:
            for v in inventories:
                # dummy = make_inventory(v) # regular identified stations  ### TO DO the oprhans will replace the extended stat conf file! need to be fixed! 
                dummy = make_inventory_orphans_mobile(v) # orphans 
        else:
            p = Pool(n_pool)                
            func = partial(make_inventory) # regular identified stations 
            out = p.map(func, inventories)  # orphans 
            
            func = partial(make_inventory_orphans_mobile)
            out = p.map(func, inventories)       
            
            
elif WHAT == 'CUON':
        dummy = make_CUON()

elif WHAT == 'ORPHAN':
        dummy = make_ORPHAN()
        
elif WHAT== 'MERGE':
        dummy = merge_inventories()
        
elif WHAT == 'MOBILE':
    dummy = make_inventory_orphans_mobile(v)

    


