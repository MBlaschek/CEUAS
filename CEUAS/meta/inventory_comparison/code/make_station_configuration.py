import os,sys
import pandas as pd
import urllib.request
import glob
from tqdm import tqdm
import numpy as np
from ast import literal_eval


pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)


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
                'file_min_date': 'start_date',
                'file_max_date':'end_date',
                'city':'city',
                'variables':'observed_variables'
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
inventories = [ 'era5_2', 'era5_1759', 'era5_1761', 'era5_3188', 'bufr', 'ncar', 'igra2', 'era5_1']

            
inventories = [ 'era5_2', 'era5_1759',]
#inventories = ['era5_1761', 'era5_3188', 'bufr',]
#inventories = [ 'ncar', 'igra2', 'era5_1']


'''
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
'''




# crate one single inventory out of the CUON datasets
inv = []
for s in glob.glob('station_configuration/*_station_configuration*'):
    if 'CUON' in s:
        continue
    df = pd.read_csv(s, sep='\t')
    print(s , ' ' , len(df))
    inv.append(df)
    
inventory_cuon = pd.concat(inv)

inventory_cuon = inventory_cuon[ [ c for c in inventory_cuon.columns if c != 'record_number'  ] ]
inventory_cuon = inventory_cuon.dropna(subset = ['primary_id'])
inventory_cuon = inventory_cuon.reset_index()

all_ids = np.unique(inventory_cuon['primary_id'].astype(str))

# storing combined cun dictionary 
combined_cuon= {}
for c in inventory_cuon.columns:
    if c in  'index' or 'Unnamed' in c :
         continue     
    combined_cuon[c] = []


all_ids = all_ids

"""
for i in tqdm(all_ids):
    indices = np.where( inventory_cuon['primary_id'].astype(str) == i )[0]
    df_i = inventory_cuon.iloc[indices]

    for c in list(inventory_cuon.columns):
        if c in  'index' or 'Unnamed' in c :
             continue 

        values = [f.replace('[','').replace(']','') for f in df_i[c].astype(str).values if isinstance(f,str) ]
        #if c == 'observed_variables':
        #   values = np.unique( [f.split(',') for f in values ] ) 
        cleaned = []
        if c == 'observed_variables':
            print(0)
        for v in values:
            if c in ['observed_variables']:
                lista = v.split(',')
                for val in lista:
                    val = val.replace(' ', '')
                    try:
                        vv = str(literal_eval(val))
                    except:
                        vv =  str(val)          
                        
                    if vv == '38':
                        vv = '138'
                    if vv == '85':
                        vv = '126'                    
                        
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
    
    

d = pd.DataFrame(combined_cuon)
d.to_csv('station_configuration/CUON_station_configuration.csv', sep='\t')



"""
   
combined_cuon = pd.read_csv('station_configuration/CUON_station_configuration.csv', sep='\t')
ids_cuon = np.unique(combined_cuon['primary_id'].astype(str))

print(0)



# merge IGRA, GRUAN
igra_f = pd.read_excel('Station_configuration_.xlsx', sheet_name='IGRA')
igra_f = igra_f[ [ c for c in igra_f.columns if c != 'record_number'  ] ]
 
gruan_f = pd.read_excel('Station_configuration_.xlsx', sheet_name='GRUAN')
gruan_f = gruan_f[ [ c for c in gruan_f.columns if c != 'record_number'  ] ]


ids_igra_f = igra_f['primary_id']



# combining all 
combining_all = {}
for c in inventory_cuon.columns:
    if c in  'index' or 'Unnamed' in c or c in 'record_index':
        continue     
    combining_all[c] = []
    

for index, row in tqdm( combined_cuon.iterrows() ):     
    pi = row.primary_id 
    
    #pi = '0-20001-0-10393'
    
    igra = igra_f.loc[igra_f['primary_id'] == pi ]
    igra = igra.reset_index()
    gruan = gruan_f.loc[gruan_f['primary_id'] == pi ]
    gruan = gruan.reset_index()
    
    if index > 100000:
        break 
    for c in inventory_cuon.columns:
        if c in ['index', 'Unnamed: 0', 'record_index']:
            continue
            
        if len(igra) ==0 and len(gruan) ==0:
            s = row[c]
            if s == 'None':
                s = ''
                
        else:
            s = row[c] 
            if s == 'None':
                s = ''
                
            if c in ['secondary_id', 'metadata_contact', 'station_name', 'city']: # summing alls
                if c == 'station_name':
                    print(0)

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

            elif c in ['start_date','end_date', 'observed_variables']: # CUON is most comprehensive
                    s = row[c]
    
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
        combining_all[c].append(s)     


""" Add station only found in GRUAN """
ids_gruan = gruan_f['primary_id'].values
ids_cuon = combined_cuon['primary_id'].values
miss = [s for s in ids_gruan if s not in ids_cuon]




    
    
""" Add station only found in GRUAN """
ids_gruan = gruan_f['primary_id'].values
ids_cuon = combined_cuon['primary_id'].values
only_gruan = [s for s in ids_gruan if s not in ids_cuon]

for station in only_gruan:     
    df = gruan_f.loc[gruan_f['primary_id'] == station ]
    for c in combining_all.keys():
        if c in ['index', 'Unnamed: 0', 'record_index']:
            continue
        s = df[c].values[0]
        combining_all[c].append(s)
    

# Saving the complete dataframe 
df = pd.DataFrame(combining_all)
df.insert(2, 'record_number', list(range(len(df))) )
df.to_csv('station_configuration/NEW_all_combined_station_configuration.csv' , sep='\t')


print('--- Finished with the global inventory --- ')

