#!/usr/bin/env python
# coding: utf-8

# # Convert original ERA5 3188 files
# 
# TO DO
# - harvest original files /mnt/users/scratch/leo/scratch/giub2.1
# - replace ERA5 3188 before 1940 
# - after 1940 check if ERA5 2 are correct i.e. cases where no pressure and geopot nans
# 
# 

# In[1]:


import pandas as pd 
pd.set_option('mode.chained_assignment', None)

import os,sys
import numpy as np
from tqdm import tqdm

from multiprocessing import Pool
from functools import partial

import random

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)

from datetime import datetime 
def get_time():
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")    
    return current_time


# In[2]:


### QUICK EXAMPLE

# di = '/mnt/users/scratch/leo/scratch/giub2.1'
# f = '5518B.txt'
# file = di + '/' + f

# df = pd.read_csv(file, sep = '\t').astype(str)
# df
#columns = df.columns
#columns
#len(columns)
#for c in columns[:30]:
#    print(c)
    
# p1
# p1flag
# p1GPH
# p1GPHflag
# p1T
# p1Tflag
# p1WD
# p1WDflag
# p1WS
# p1WSflag
# p1U
# p1Uflag
# p1V
# p1Vflag
# p1RH
# p1RHflag
# p1DIFFTAU
# p1DIFFTAUflag
# p1SH
# p1SHflag

# initialize dictionary containing all data
all_data = {} 
for v in ['date' , 'time', 'temp' , 'gph' , 'wspeed' , 'wdir', 'uwind', 'vwind' , 'rh' , 'sh' , 'z_coordinate' , 'z_coordinate_type' , 'days' , 'months']:
    all_data[v] = [] 

def convert(out_dir, file, return_df=False):
    """ Converting the original txt files where each row contains all the variables measured at different z_coordinates as columns.
    Convert to a standard format where each row is a different level measurements, and columns are variables.
    Complication:
     the z coordinate is given as either pressure e.g. p1, p2... pMax or heigh e.g. h1,h2...h Max.
    Sometimes you have p, sometimes you have h.
    And the maximum number of levels is arbitrary i.e. Max=50 or Max=100.
    See implementation.
    """
    
    def create_date_time(month, day, year, hour, minutes):
        # creating dates

        if len(month) <2:
            month = '0' + month
        if len(day) <2:
            day = '0' + day
        date = year + '-' + month + '-' + day 

        # creating datetimes objects
        if len(hour) <2:
            hour = '0' + hour 
        if len(minutes) <2:
            minutes = '0' + minutes 
        time = hour + ':' + minutes
        #print(date ,  '   ' , time )   
        return date, time

    # read df with pandas 
    print('Reading df from original file === ' , file )
    df = pd.read_csv(file, sep = '\t').astype(str)
    print('Done Reading df === ' , file )
    
    # storing lists, drop dataframe
    dic_df = {}
    for c in df.columns:
        dic_df[c] = df[c].values
        
    total_rows = len(df)    
    columns = df.columns 
    del df
        
        
    max_level = 1
    if 'p1' in columns:
        for i in range(1,200):
            if 'p'+str(i) in columns:
                max_level=i
            else:
                break
    else:
        if 'h1' in columns:
            for i in range(1,200):
                if 'h'+str(i) in columns:
                    max_level=i
                else:
                    break
                    
    levels = list(range(1,max_level))
    # determine if levels in height or pressure 
    if 'p1' in columns:
            letter_level = 'p'
            z_type = '1'
    elif 'h1' in columns:  # use CDM convention for geo pot height in meters 
            letter_level = 'h'
            z_type = '0'
    

    
    for i in tqdm(range(total_rows)):

        z_coordinates  = [ dic_df[letter_level+str(lev)][i] for lev in levels ]
        z_types = [z_type for lev in levels ] 
        
        month = dic_df['month'][i]  
        day = dic_df['day'][i]  
        year = dic_df['year'][i]  
        hour = dic_df['hourUTC'][i]
        minutes = dic_df['minutes'][i] 
        
        date, time = create_date_time(month, day, year, hour, minutes)

        dates    = [date for l in levels ]
        times    = [time for l in levels ]
        days     = [day for lev in levels ]
        months = [month for lev in levels ]
        
        # extract data for all levels
        temp   = [ dic_df[letter_level + str(lev)+'T'][i] for lev in levels ]
        if z_type=='1':
            gph = [ dic_df[letter_level + str(lev)+'GPH'][i] for lev in levels ]
        else:
            gph = [-999 for lev in levels ] 

        wspeed = [ dic_df[letter_level + str(lev)+'WS'][i] for lev in levels ]
        wdir      = [ dic_df[letter_level + str(lev)+'WD'][i] for lev in levels ]
        uwind   = [ dic_df[letter_level + str(lev)+'U'][i]    for lev in levels ]
        vwind   = [ dic_df[letter_level + str(lev)+'V'][i]    for lev in levels ]
        rh         = [ dic_df[letter_level + str(lev)+'RH'][i]  for lev in levels ]
        sh         = [ dic_df[letter_level + str(lev)+'SH'][i]  for lev in levels ] 
        
        all_data['date'].extend(dates)
        all_data['time'].extend(times)

        all_data['z_coordinate'].extend(z_coordinates)
        all_data['z_coordinate_type'].extend(z_types)
        
        all_data['temp'].extend(temp)
        all_data['gph'].extend(gph) # gph in meters  
        all_data['wspeed'].extend(wspeed)
        all_data['wdir'].extend(wdir)
        all_data['uwind'].extend(uwind)
        all_data['vwind'].extend(vwind)
        all_data['rh'].extend(rh)
        all_data['sh'].extend(sh)

        all_data['days'].extend(days)
        all_data['months'].extend(months)
        
        
    # clean the dictionary from all emtpy observation i.e. when neither temp, not wspeed, wdir, uwind, udir are available 
    
    keep_indices = []
    for i in range( len(all_data['sh']) ):
        temp = all_data['temp'][i]
        winddir = all_data['wdir'][i]
        wspeed = all_data['wspeed'][i]
        uwind, vwind = all_data['uwind'][i], all_data['vwind'][i] 
        if '999' in temp and '999' in winddir and '999' in wspeed and '999' in vwind and '999' in uwind:
            pass
        else:
            keep_indices.append(i)
        
    all_data_valid = {}
    for c in all_data.keys():
        all_data_valid[c] = [ all_data[c][i] for i in keep_indices ]
        
        
    print('Creating dataframe for file :::' , file , '   ' ,   get_time() )
    all_data_df = pd.DataFrame.from_dict(all_data_valid)
    print('Done Creating dataframe for file :::' , file , '   ' ,   get_time() )
    
    # remove impossible dates  e.g. 31st of June 

    try:
        all_data_df['date'] = pd.to_datetime(all_data_df['date'] )
    except:
        ind = all_data_df.loc[ (all_data_df.months == '6') & (all_data_df.days == '31')  ]
        if len(ind) >0 :
            all_data_df = all_data_df.drop(index=ind.index)
        
        
    all_data_df = all_data_df.reset_index()
    
    all_data_df['gph'] = all_data_df['gph'].astype(float)
    valid = np.where(all_data_df.gph > 0 )[0] 

    all_data_df['geopotential'] = all_data_df['gph'] # gph in meters 
    all_data_df['geopotential'][valid] = all_data_df['geopotential'][valid].astype(float) * 9.80665 # convert height to geopotential 
        
    if not os.path.isdir(out_dir):
        os.system('mkdir ' + out_dir)
    out_name = out_dir + '/' + file.split('/')[-1] + '_converted_csv.csv'
    
    
    print('Writing dataframe to CSV for file :::' , file , '   ' ,   get_time() )
    all_data_df.to_csv(out_name,  sep = '\t' , index=False) 
    print('Done Writing dataframe to CSV for file :::' , file , '   ' ,   get_time() )


    print('Finished file::: ' , file, '   '  , get_time() )
    if return_df:
        return all_data_df

# Example single file - TEST 
di = '/mnt/users/scratch/leo/scratch/giub2.1'
out_dir = '/scratch/das/federico/databases_service2/giub2.1_ERA5_3188_14032023' 


'''
f = '5518B.txt'
file = di + '/' + f
file_df = pd.read_csv(file, sep = '\t').astype(str)

file_df.columns

df = convert(file , return_df = True)

dd = df.loc[df['date'] == pd.Timestamp('1925-05-14')]
dd = dd.loc[dd['time'] == '13:00']

ff = '/mnt/users/scratch/leo/scratch/giub2.1/4862.txt'
'''

'''
# quick test
out_dir = '/scratch/das/federico/databases_service2/giub2.1_ERA5_3188_14032023' 

f = '5066B.txt'q
file = di + '/' + f 
file_df = pd.read_csv(file, sep = '\t').astype(str)
df = convert(out_dir, file , return_df = True)
'''



'''
f = '4882.txt'
file = di + '/' + f 
file_df = pd.read_csv(file, sep = '\t').astype(str)
df = convert(out_dir, file , return_df = True)

a = 0
'''

'''
#f = '4963.txt'
f = '6015A.txt'
file = di + '/' + f

from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Before Time =", current_time)

df = convert(out_dir, file , return_df = True)

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("After Time =", current_time)

'''

### Running Multiproc 

def run(out_dir, exception, file):
    exception = False
    if exception:
        try:        
            convert(out_dir, file, return_df=False)
        except:
            a = open('Failed_ERA_3188-giub2.1.txt' , 'a+' )
            a.write(file + '\n')
    else:
        convert(out_dir, file)
        
out_dir = '/scratch/das/federico/databases_service2/giub2.1_ERA5_3188_02042023'

exception = False
func = partial(run, out_dir, exception)
p = Pool(12)        
POOL = True
skip_processed = True

if skip_processed:
    processed = [f.split('_converted')[0] for f in os.listdir(out_dir)]
    files = [di + '/' + f for f in os.listdir(di) if '.txt' in f and 'py' not in f and 'upper' not in f and f.split('/')[-1] not in processed ]
else:
    files = [di + '/' + f for f in os.listdir(di) if '.txt' in f and 'py' not in f and 'upper' not in f ]
    
files = [f for f in files  if '.tx' not in f  ]
random.shuffle(files)


#files = [ f for f in files if  '4702.txt' in f ]
print('Total Missing files ::: ' , len(files) )
#print(files)
if POOL:
    out = p.map(func, files)      
else:
    for f in files:
        convert(out_dir, f, False)
        
    
