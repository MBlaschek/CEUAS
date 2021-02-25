""" Analyze Payerne merged file """
import xarray as xr
import pandas as pd

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)



p = '/raid60/scratch/federico/DATABASE_JANUARY2021_FIXED_sensor/0-20000-0-06610_CEUAS_merged_v0.nc'
df = xr.open_dataset(p, engine = 'h5netcdf' , group = 'header_table', decode_times = True ).to_dataframe()

df = df[['report_timestamp', 'record_timestamp', 'report_id', 'source_id']]

source_con = xr.open_dataset(p, engine = 'h5netcdf' , group = 'source_configuration', decode_times = True ).to_dataframe()
files = source_con.index

df['source_file'] = files


rep_ids = df['report_id']
s = list(set(rep_ids))

print("There is a total of " , len(rep_ids) , ' but the distinct report_ids are: ' , len(s) )


import collections
duplicated = [item for item, count in collections.Counter(rep_ids).items() if count > 1] 
first = duplicated[0]
print("The first duplicated reportid is: " ,first )
print("*** Find it in the header_table ***" )
find_df = df.loc[ df['report_id'] == first ]
print(find_df)

for d in duplicated:
    find_df = df.loc[ df['report_id'] == d ]
    print(find_df)
    
print(0)
            
