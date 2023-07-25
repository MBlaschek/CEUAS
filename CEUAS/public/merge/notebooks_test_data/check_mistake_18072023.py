""" Check missing values for latitude,longitude 
      merged VS resorted files
      See email Paul 14/07/2023 
      
      AIM: check that the nans and misisng values for report ids come from steps subsequent to merging 
      i.e. during resorting """


import h5py as h5
import numpy as np
import pandas as pd 

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)




def make(file, var=126, date='1980-06-01'):
    """ Extract info from merged (pre-resorting) file """
    
    ff = h5.File(file, 'r')
    ts = pd.to_datetime(   ff['recordtimestamp'][:] , unit='s',  origin=pd.Timestamp('1900-01-01') )
    
    ind =  ff['recordindex'][:]
    
    df = pd.DataFrame.from_dict( {'inde': ind, 'timestamp': ts })
    
    
    # selecting variables
    variab = ['date_time' , 'observed_variable' , 'observation_value', 'latitude' , 'longitude', 'z_coordinate' , 'source_id' , 'report_id']
    
    # getting index for date
    ind = np.where(df.timestamp >= pd.Timestamp(date)) [0][0]
    indr = df.inde[ind]
    dic = {}
    for v in variab:
        if v in ['source_id', 'report_id']:
            d = [ b''.join(c) for c in  ff['observations_table'][v][indr : indr+1000 ]]
            dic[v] = d
        elif v =='date_time':
            dic[v]  = pd.to_datetime(   ff['observations_table'][v][indr : indr+1000  ], unit='s',  origin=pd.Timestamp('1900-01-01') )
            
        else:
            dic[v] = ff['observations_table'][v][indr : indr+1000  ]
      
    df = pd.DataFrame.from_dict( dic )
        
    # selecting variable, date 
    df =df[df.observed_variable == var]
    df =df.sort_values(by=['z_coordinate'])
    df =df[df.date_time == date]
    
    
    return df
    
    
file = '/scratch/das/federico/REMERGE_GIUB_12072023/0-20001-0-11035_CEUAS_merged_v1.nc'    
merged = make(file, var=126, date='1980-06-01')
    
"""
merged
     date_time  observed_variable  observation_value  latitude  longitude  z_coordinate  source_id           report_id
0   1980-06-01  126                245.500000         48.25     16.370001  500.0         b'era5_1'  b'100000000046090'
9   1980-06-01  126                237.500000         48.25     16.370001  760.0         b'era5_1'  b'100000000046090'
15  1980-06-01  126                236.300003         48.25     16.370001  1000.0        b'era5_1'  b'100000000046090'
24  1980-06-01  126                229.300003         48.25     16.370001  1370.0        b'era5_1'  b'100000000046090'
26  1980-06-01  126                223.899994         48.25     16.370001  2000.0        b'era5_1'  b'100000000046090'
31  1980-06-01  126                223.500000         48.25     16.370001  2050.0        b'era5_1'  b'100000000046090'
33  1980-06-01  126                220.899994         48.25     16.370001  3000.0        b'era5_1'  b'100000000046090'
43  1980-06-01  126                219.300003         48.25     16.370001  5000.0        b'era5_1'  b'100000000046090'
52  1980-06-01  126                218.899994         48.25     16.370001  5500.0        b'era5_1'  b'100000000046090'
58  1980-06-01  126                221.100006         48.25     16.370001  7000.0        b'era5_1'  b'100000000046090'
67  1980-06-01  126                221.899994         48.25     16.370001  7600.0        b'era5_1'  b'100000000046090'
68  1980-06-01  126                221.899994         48.25     16.370001  9000.0        b'era5_1'  b'100000000046090'
78  1980-06-01  126                224.100006         48.25     16.370001  10000.0       b'era5_1'  b'100000000046090'
"""
    
    
    
    
#############################3   RESORTED, VIENNA
file_res = '/mnt/users/scratch/leo/scratch/converted_v13/long/0-20001-0-11035_CEUAS_merged_v1.nc'    
ff = h5.File(file_res, 'r')


ts = pd.to_datetime(  ff['recordindices']['recordtimestamp'], unit='s',  origin=pd.Timestamp('1900-01-01') )

df_ind = pd.DataFrame.from_dict( {'inde': ff['recordindices']['126'][0:-1]  , 'timestamp': ts } )


date='1980-06-01'

ind = np.where(df_ind.timestamp >= pd.Timestamp(date)) [0][0]
indr = df_ind.inde[ind]


dic = {}
variab = ['date_time' , 'observed_variable' , 'observation_value', 'latitude' , 'longitude', 'z_coordinate' , 'source_id' , 'report_id']
for v in variab:
    if v in ['source_id', 'report_id']:
        d = [ b''.join(c) for c in  ff['observations_table'][v][indr : indr+1500 ]]
        dic[v] = d
    elif v == 'date_time':
            dic[v]  = pd.to_datetime(   ff['observations_table'][v][indr : indr+1500 ], unit='s',  origin=pd.Timestamp('1900-01-01') )
            
    else:
        dic[v] = ff['observations_table'][v][indr : indr+1500 ]
    
resorted = pd.DataFrame.from_dict( dic )
resorted =resorted.sort_values(by=['z_coordinate'])
resorted =resorted[resorted.date_time == '1980-06-01']


"""
resorted
   date_time  observed_variable  observation_value  latitude  longitude  z_coordinate  source_id                 report_id
0  1980-06-01  126                245.500000         48.25     16.370001  500.0         b'era5_1'  b'100000000104386'      
1  1980-06-01  126                240.783417        NaN       NaN         640.0         b'nnnnnn'  b'nnnnnnnnnnnnnnnnnnnnn'
2  1980-06-01  126                237.500000         48.25     16.370001  760.0         b'era5_1'  b'100000000104386'      
3  1980-06-01  126                237.010635        NaN       NaN         850.0         b'nnnnnn'  b'nnnnnnnnnnnnnnnnnnnnn'
4  1980-06-01  126                236.300003         48.25     16.370001  1000.0        b'era5_1'  b'100000000104386'      
5  1980-06-01  126                235.427902        NaN       NaN         1040.0        b'nnnnnn'  b'nnnnnnnnnnnnnnnnnnnnn'
6  1980-06-01  126                229.300003         48.25     16.370001  1370.0        b'era5_1'  b'100000000104386'      
7  1980-06-01  126                223.899994         48.25     16.370001  2000.0        b'era5_1'  b'100000000104386'      
"""