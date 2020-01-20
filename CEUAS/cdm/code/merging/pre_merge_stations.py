""" Pre-merging utility for station data split into several files.
Example: MOSCOW ,  0-20000-0-27612 
 ['/raid60/scratch/leo/scratch/era5/odbs/3188/era5.3188.conv.C:4629', '/raid60/scratch/leo/scratch/era5/odbs/3188/era5.3188.conv.C:4567', '/raid60/scratch/leo/scratch/era5/odbs/3188/era5\
.3188.conv.C:5246']

"""

# variables ==> recordtimestamp, recordindex  (to be calculated new)


# groups ==>  to be updated: observations_table, header_table, era5fb 
#            ==> to be simply copied: 'station_configuration', 'crs' , 'observed_variable', 'units' , 'z_coordinate_type' , 'station_type'



import os,sys
import netCDF4 as nc
import pandas as pd
import xarray as xr 
import numpy as np

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)


class CombineNetCDF():
    """ Functionalitites to read netCDF files in input """
    
    def __init__(self, files ):
        self.files = files
        self.crs = ''
        self.observed_variables = ''
        self.units = ''
        self.z_coordinate_type = ''
        self.station_type = ''
        self.station_configuration = ''
        
        self.observations_table = ''
        self.era5fb = ''
        self.header_table = ''
        
        self.data = {}  #containing all the df to be combined
        self.combined = {}
        
    def read_netCDF(self):
        """ Reading the groups from the netCDF file using xarray """
        
        for f in self.files:     
            data = {}      
            for group in ['observations_table' , 'era5fb' , 'header_table' ]:    
                data[group] = xr.open_dataset(f , engine = 'h5netcdf' , group = group).to_dataframe() 

            if files.index(f) == 0:
                for group in ['station_configuration', 'crs' , 'observed_variable', 'units' , 'z_coordinate_type' , 'station_type']:    
                    self.combined[group] = xr.open_dataset(f1 , engine = 'h5netcdf' , group = group) # no need to save as panda, they will be copied as they in xarray are to the output file  
                    
            self.data[f] = data         

    def combine_df(self):
        """ Analize the dataframes for the  'observations_table' , 'era5fb' , 'header_table'  and combine """
        
        observations_tables, header_tables, era5fb_tables = [], [], [] 
        
        for k in  self.data.keys():
            observations_tables.append(self.data[k]['observations_table'] )
            header_tables.append(self.data[k]['header_table'] )
            era5fb_tables.append (self.data[k]['era5fb'])        
 
        # observations table
        observations_tables_combined = pd.concat(observations_tables)
        observations_tables_combined = observations_tables_combined.sort_values(by = ['date_time', 'z_coordinate' ] )    

        # header_table        
        header_tables_combined = pd.concat(header_tables)
        header_tables_combined = header_tables_combined.sort_values(by = ['record_timestamp' ] )    
        
        # era5fb     
        era5fb_tables_combined= pd.concat(era5fb_tables)        
        try:  # different sorting if the original source is in ODB vs all the rest of the formats 
            era5fb_tables_combined = era5fb_tables_combined.sort_values(by = ['report_timestamp' , 'vertco_reference_1@body' ] )        
        except:
            era5fb_tables_combined = era5fb_tables_combined.sort_values(by = ['date@hdr', 'time@hdr' , 'vertco_reference_1@body' ] )        
        
        self.combined['era5fb'] = era5fb_tables_combined.to_xarray()
        self.combined['header_table'] = header_tables_combined.to_xarray()
        self.combined['observations_table'] = observations_tables_combined.to_xarray()
                
                
        print('Done combining dataframes')



    def find_date_indices(self):
        """ Extracts the list of observation dates, and store the indices of the first and last observations. Copy from the netCDF_CDM_converter script """            
        datetime = self.combined['observations_table']['date_time'].values    
        
        date_times, indices, counts = np.unique(datetime, return_counts = True, return_index= True)

        print('check')
        # convert to date_time object 
        try:
            date_times = [  datetime.strptime(  str(int(i)) , '%Y%m%d%H%M') for i in date_times ]
        except:
            print('already date_time')
            pass
    
        return np.array(indices) , date_times,  counts  
        
        
    def write_netCDF(self, dataset = '', out_dir = '' ):

        station_id = self.combined['station_configuration']['primary_id'].values[0]         
        out_file = out_dir + '/' + dataset + '_' +  station_id + '.nc'
        
        for k in  self.combined.keys():
                self.combined[k].to_netcdf(out_file, format='netCDF4', engine='h5netcdf', mode='a' , group = k)  # writing the merged observations_table

        """ Writing the new recordtimestamps and recordindex """
        date_times, indices, counts = self.find_date_indices()
        di=xr.Dataset()

        di['recordtimestamps'] = ( {'recordtimestamps' : date_times.shape } , date_times )
        di['recordindex']          = ( {'recordindex' : indices.shape } , indices )
        di.to_netcdf(out_file, format='netCDF4', engine='h5netcdf', mode='a')
                
                
        print('Done writing ' , out_file )
        

""" Try """

if __name__ == '__main__':

    out_dir = os.getcwd()  # output directory  
    out_dir = '/raid60/scratch/federico/MERGED_FILES/'
       
    f1 = 'MOSCOW/chera5.3188.conv.C:4567.nc'
    f2 = 'MOSCOW/chera5.3188.conv.C:4629.nc'
    f3 = 'MOSCOW/chera5.3188.conv.C:5246.nc'

    files = [ f1, f2, f3 ]

    """ Analize """
    cb = CombineNetCDF(files = files)
    cb.read_netCDF ()
    cb.combine_df()
    cb.write_netCDF( dataset='era5_3188' , out_dir= out_dir )

    print('hello')

'''
ot1 = xr.open_dataset(f1 , engine = 'h5netcdf')

ot1 = xr.open_dataset(f1 , engine = 'h5netcdf' , group = 'observations_table').to_dataframe() #[['observed_variable', 'z_coordinate' , 'date_time' , 'observed_variable', 'observation_value'] ]
ot2 = xr.open_dataset(f2 , engine = 'h5netcdf' , group = 'observations_table').to_dataframe() #[['observed_variable', 'z_coordinate' , 'date_time' , 'observed_variable', 'observation_value'] ]
ot3 = xr.open_dataset(f3 , engine = 'h5netcdf' , group = 'observations_table').to_dataframe() #[['observed_variable', 'z_coordinate' , 'date_time' , 'observed_variable', 'observation_value'] ]


#['observed_variables', 'z_coordinate' , 'date_time' , 'observed_variable', 'observatio
#
#n_value']


s1, s2 , s3 = set(ot1['date_time']) , set(ot2['date_time']) , set(ot3['date_time'])

tot = list(s1) + list(s2) + list(s3) 

print ( len(set(tot)) ) 

concatenated = pd.concat ( [ot1, ot2, ot3 ] )
concatenated = concatenated.drop_duplicates()
'''

"""
1. concatenate obervation tables
2. concatenate feedback
3. sort both by date_time and pressure
4. find new record_index 


"""

