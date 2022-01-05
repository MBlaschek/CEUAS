import os,sys
import pandas as pd
import urllib.request
import glob
from tqdm import tqdm
import numpy as np
from ast import literal_eval

from multiprocessing import Pool
from functools  import partial

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
        dfr = dfr.dropna( subset = ['latitude','longitude'])
        
        if not dfr.empty:
            if i == 'OSCAR':
                dfr_radio = dfr[dfr.isRadio == True ] # check if upperAir/radiosonde station flag 
                if not dfr_radio.empty:
                    return dfr_radio
                else:
                    return dfr[0:1]
            else:
                return dfr[0:1]  
            
    return dfr 
            
            
            
# directory containing each single inventory file for each station 
        

def make_inventory(v):
    """ Creates the inventory for the given database """
    
    print(' *** Creating the inventory for ' , v )
    
    inv_path = 'inventories/' 

    files = glob.glob(inv_path + '/' + v + '/' + '*inventories.csv*')
    
    # initialize each inventory dic
    stat_conf_dic = {}
    for c in stat_conf_columns:
        stat_conf_dic[c] = []
    
    for file in tqdm(files):    
        
        df = pd.read_csv(file, sep = '\t', header = 0)
        
        # select best inventory available        
        best_inv = get_best_inventory(df)
        
        if best_inv.empty:
            continue
        
        for c in stat_conf_columns:
            
            if c in special_columns.keys():
                stat_conf_dic[c].append(special_columns[c])            

            elif c in statConf_to_inventory.keys():
                if c == 'secondary_id':
                    a = best_inv[ statConf_to_inventory[c] ].values[0]
                    #print(0)
                    
                elif c == 'start_date' or c == 'end_date':
                    a = best_inv[ statConf_to_inventory[c] ].values[0]
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
                
    df = pd.DataFrame(stat_conf_dic)
    ind = df.index
    # adding record number
    df['record_number'] = list(range(1,len(df)+1))
    df.to_csv('station_configuration/' + v + '_station_configuration.csv' , sep= '\t' )           
        
    print('*** Done with the inventory ',  v )








def make_CUON():
    ''' Create one single inventory out of the CUON datasets '''
    
    inv = []
    for s in glob.glob('station_configuration/*_station_configuration*'):
        if 'CUON' in s:
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
    for c in inventory_cuon.columns:
        if c in  'index' or 'Unnamed' in c :
             continue     
        combined_cuon[c] = []
    
    # looping through all the unique ids and combining data
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
                
            # extracting min and max date 
            if c == 'start_date':
                try:
                    values = [ int ( min ( [float(s) for s in values if float(s) > 2021 ] ) ) ]
                except:
                    values = [ int ( min ( [float(s) for s in values ] ) ) ]
                    
            if c == 'end_date':
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


def merge_inventories():
    
    inventory_cuon = pd.read_csv('station_configuration/CUON_station_configuration.csv', sep='\t')
    ids_cuon = np.unique(inventory_cuon['primary_id'].astype(str))
        
    # read IGRA, GRUAN
    igra_f = pd.read_excel('Station_configuration_HARM.xlsx', sheet_name='IGRA')
    igra_f = igra_f[ [ c for c in igra_f.columns if c != 'record_number'  ] ]
     
    gruan_f = pd.read_excel('Station_configuration_HARM.xlsx', sheet_name='GRUAN')
    gruan_f = gruan_f[ [ c for c in gruan_f.columns if c != 'record_number'  ] ]

    ids_igra_f = igra_f['primary_id']
    
    # combining all the inventories
    combining_all = {}
    for c in inventory_cuon.columns:
        if c in  'index' or 'Unnamed' in c or c in 'record_index':
            continue     
        combining_all[c] = []
        
    #inventory_cuon = inventory_cuon[:100]
    for index, row in tqdm( inventory_cuon.iterrows() ):     
        pi = row.primary_id 
        
        #if index > 8500:
        #    break
        
        igra = igra_f.loc[igra_f['primary_id'] == pi ]
        igra = igra.reset_index()
        gruan = gruan_f.loc[gruan_f['primary_id'] == pi ]
        gruan = gruan.reset_index()
        
        for c in inventory_cuon.columns:
            if c in ['index', 'Unnamed: 0', 'record_index']:
                continue
                
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
                    #if c == 'station_name':
                    #    print(0)
    
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
    
                #elif c in ['start_date','end_date']: # CUON is most comprehensive
                #    try:
                #        s = row[c][:4]
                #    except:
                #        s = ''
                        
                elif c in ['observed_variables']: 
                    if len(gruan) >0:
                        
                        s = s + ',57'
                    
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
            if c == 'start_date' or c == 'end_date':
                s = df[c].dt.strftime('%Y%m%d').values[0]
            else:
                s = df[c].values[0]
            
            combining_all[c].append(s)
        
    # Saving the complete dataframe 
    df = pd.DataFrame(combining_all)
    df.insert(2, 'record_number', list(range(len(df))) )
    df.to_csv('station_configuration/all_combined_station_configuration.csv' , sep='\t')
    df.to_csv('station_configuration/all_combined_station_configuration' )
    
    
    print('--- Finished with the global inventory --- ')
    
    
    
    

WHAT = 'MERGE'
POOL = True
n_pool = 40

if WHAT == 'inventory':
    # inventories to process
    inventories_all = [ 'era5_2', 'era5_1759', 'era5_1761', 'era5_3188', 'bufr', 'ncar', 'igra2', 'era5_1']
    #inventories = ['era5_1761', 'era5_3188', 'bufr',]
    
    inventories = inventories_all
    
    if not POOL:
        for v in inventories:
            dummy = make_inventory(v)
    else:
        p = Pool(n_pool)                
        func = partial(make_inventory)
        out = p.map(func, inventories)       
    
    print(0)
    
elif WHAT == 'CUON':
    dummy = make_CUON()

elif WHAT== 'MERGE':
    dummy = merge_inventories()



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



'''
inventory_cuon = pd.read_csv('station_configuration/CUON_station_configuration.csv', sep='\t')
ids_cuon = np.unique(inventory_cuon['primary_id'].astype(str))

print(0)

# merge IGRA, GRUAN
igra_f = pd.read_excel('Station_configuration_HARM.xlsx', sheet_name='IGRA')
igra_f = igra_f[ [ c for c in igra_f.columns if c != 'record_number'  ] ]
 
gruan_f = pd.read_excel('Station_configuration_HARM.xlsx', sheet_name='GRUAN')
gruan_f = gruan_f[ [ c for c in gruan_f.columns if c != 'record_number'  ] ]


ids_igra_f = igra_f['primary_id']


# combining all the inventories
combining_all = {}
for c in inventory_cuon.columns:
    if c in  'index' or 'Unnamed' in c or c in 'record_index':
        continue     
    combining_all[c] = []
    

for index, row in tqdm( inventory_cuon.iterrows() ):     
    pi = row.primary_id 
    
    #pi = '0-20001-0-10393'
    if index > 8500:
        
        break
    
    igra = igra_f.loc[igra_f['primary_id'] == pi ]
    igra = igra.reset_index()
    gruan = gruan_f.loc[gruan_f['primary_id'] == pi ]
    gruan = gruan.reset_index()
    
    for c in inventory_cuon.columns:
        if c in ['index', 'Unnamed: 0', 'record_index']:
            continue
            
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
                #if c == 'station_name':
                #    print(0)

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

            elif c in ['start_date','end_date']: # CUON is most comprehensive
                try:
                    s = row[c][:4]
                except:
                    s = ''
                    
            elif c in ['observed_variables']: 
                s = row[c][:4]
                
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
        if c == 'start_date':
            lista =  s.split(',') 
            lista = [l for l in lista if l != 'nan' ]        
            s = min(lista)
        if c == 'end_date':
            lista =  s.split(',') 
            lista = [l for l in lista if l != 'nan' ]        
            s = max(lista)
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
                
        if s == 'nan':
            s = ''      
                
            #print(0)
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
        s = df[c].values[0]
        combining_all[c].append(s)
    

# Saving the complete dataframe 
df = pd.DataFrame(combining_all)
df.insert(2, 'record_number', list(range(len(df))) )
df.to_csv('station_configuration/all_combined_station_configuration.csv' , sep='\t')


print('--- Finished with the global inventory --- ')

'''





"""
# old loop for inventories wthout function 
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
                if c == 'secondary_id':
                    a = best_inv[ statConf_to_inventory[c] ].values[0]
                    print(0)
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



