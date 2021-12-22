import os,sys
import pandas as pd
import urllib.request
import glob
from tqdm import tqdm

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

# mapping the column names from each inveotry csv file, to the column names in the CDM
inventory_to_statConf = { 'WIGOS_calc' : 'primary_id' ,
                'file_statid' : 'secondary_id',
                'station_name' : 'station_name',
                'latitude':'latitude',
                'longitude':'longitude',
                'start_date': 'start_date',
                'end_date':'end_date',
                'city':'city',
                }

statConf_to_inventory = dict( zip( inventory_to_statConf.values(), inventory_to_statConf.keys()) )

# supplementary columsn to fill in CUON
metadata_contact = 'L.Haimberger'
metadata_contact_role = '0'


# to decide 
primary_id_scheme = 0  # see Fabios list 
secondary_id_scheme = 1 # see Fabios list 
station_crs = 0
station_type = 1

platform_type = 0
platform_sub_type = 63

# to be copied when available from GRUAN stations 
operating_institute = ''
operating_territory = '' 
telecommunication_method = '' 
station_automation = '' 
measuring_system_model = '' 
measuring_system_id = ''


special_columns = {'primary_id_scheme':primary_id_scheme,
                   'secondary_id_scheme':secondary_id_scheme,
                   'station_crs':station_crs,
                   'station_type':station_type,
                   'platform_type':platform_type,
                   'platform_sub_type':platform_sub_type, 
                   'operating_institute':operating_institute,
                   'operating_territory':operating_territory,
                   'telecommunication_method':telecommunication_method,
                   'station_automation':station_automation,
                   'measuring_system_model':measuring_system_model,
                   'measuring_system_id':measuring_system_id, 
                   'metadata_contact_role':metadata_contact_role,
                   'metadata_contact':metadata_contact
                   }



print(stat_conf)





# initialize empty dic as the stat_conf table 

    
def get_best_inventory(df):
    """ Select best inventory among the ones available.
    Give preference to upper-air / radiosonde station from OSCAR, then generic OSCAR,
    then IGRA2, then WBAN and finaly CHUAN inventory """
    
    for i in ['OSCAR', 'IGRA2', 'WBAN' , 'CHUAN']:
        dfr = df[ df.inventory == i ]
        if not dfr.empty:
            if i == 'OSCAR':
                dfr_radio = dfr[dfr.isRadio == True ] # check if upperAir/radiosonde station flag 
                if not dfr_radio.empty:
                    return dfr_radio
                else:
                    return dfr[0:1]
            else:
                return dfr[0:1]            
            
            
            
# directory containing each single inventory file for each station 
inv_path = 'inventories/' 
            
# inventories to process
inventories = [ 'era5_2','bufr', 'ncar']
            
"""
for v in inventories:
    files = glob.glob(inv_path + '/' + v + '/' + '*inventories.csv*')
    
    # initialize each inventory dic
    stat_conf_dic = {}
    for c in stat_conf_columns:
        stat_conf_dic[c] = []
    
    for file in tqdm(files):    
        df = pd.read_csv(file, sep = '\t', header = 0)
        
        # select best inventory available        
        best_inv = get_best_inventory(df)
        
        for c in stat_conf_columns:
            
            if c in special_columns.keys():
                stat_conf_dic[c].append(special_columns[c])            
            
            elif c in statConf_to_inventory.keys():
                    try:
                        value = best_inv[ statConf_to_inventory[c] ].values[0]
                        if ('999' in str(value) or value==999 ):
                            stat_conf_dic[c].append('')
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
                
    df = pd.DataFrame(stat_conf_dic)
    ind = df.index
    # adding record number
    df['record_number'] = list(range(1,len(df)+1))
    df.to_csv('station_configuration/' + v + '_station_configuration.csv' , sep= '\t' )           
            
print(0)
"""


# crate one single inventory out of the CUON datasets
inv = []
for s in glob.glob('station_configuration/*_station_configuration*'):
    if 'CUON' in s:
        continue
    df = pd.read_csv(s, sep='\t')
    print(s , ' ' , len(df))
    inv.append(df)
    
inventory_cuon = pd.concat(inv)
inventory_cuon = inventory_cuon[ [ c for c in inventory_cuon.columns if c == 'record_number'  ]]

inventory_cuon = inventory_cuon.drop_duplicates( subset=['primary_id'])
inventory_cuon.to_csv('station_configuration/CUON_station_configuration.csv' , sep= '\t' )           

print(0)



"""
# merge IGRA, GRUAN
igra_f = pd.read_excel('Station_configuration_.xlsx', sheet_name='IGRA')
gruan_f = pd.read_excel('Station_configuration_.xlsx', sheet_name='GRUAN')

ids_igra_f = igra_f['primary_id']
ids_gruan_f = gruan_f['primary_id']
"""