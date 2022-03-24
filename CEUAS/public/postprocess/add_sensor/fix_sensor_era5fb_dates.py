""" Extract the wmo tables and creates a table for converting the old-to-new codes wrt the dates when it is applicable.
See table page A-398
https://library.wmo.int/doc_num.php?explnum_id=10235
"""

import pandas as pd
import numpy as np

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)


# Fix sensor_id codes from feedback wrt dates 

#file = '/users/staff/federico/GitHub/CEUAS_master_SEPTEMBER2021/CEUAS/CEUAS/public/merge/sensor_configuration_all.csv'
#df = pd.read_csv(file, sep = '\t')

wmo = pd.read_csv('data/table_BUFR_radiosonde.csv' , sep=',' , header=1 , names = ['date', 'table_1', 'sensor_id', 'comments'] )




dates = wmo.date.values
dates_n = []
old_sensor = []
new_sensor = []


def make_df_convertion(wmo):
    """ Creates a table for converting the old-to-new codes wrt the dates when it is applicable.
    See table page A-398
    https://library.wmo.int/doc_num.php?explnum_id=10235
    """    
    for index, row in wmo.iterrows():
        d = row.date
        if type(d) == str and '/' in d :
            print(d)
            s = d.split('/')
            d,m,y = s[0] , s[1], s[2]
            
            if len(y) <3:
                y = '20' + y 
            
            dt = pd.Timestamp( '-'.join([y,m,d]) ).date()
            dates_n.append(dt)
            old_sensor.append(row.table_1)
            new_sensor.append(row.sensor_id)
            
        else:
            pass    
        
    df = pd.DataFrame( {'date':dates_n , 'old_id': old_sensor , 'new_id': new_sensor}  )
    
    return df 
        
        
a = 0