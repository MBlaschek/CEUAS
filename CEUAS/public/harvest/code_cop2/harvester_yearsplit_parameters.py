"""
parameters for harvester run
Will be copied inside the output directory so that the parameters will be stored

"""


### Datasets sources
database_service2_path = '/scratch/das/federico/databases_service2/'

datasets_path = { 'era5_1': '/mnt/users/scratch/leo/scratch/era5/odbs/1' ,
      'era5_2': '/mnt/users/scratch/leo/scratch/era5/odbs/2',
      'era5_3188': '/mnt/users/scratch/leo/scratch/era5/odbs/3188',
      'era5_1759': '/mnt/users/scratch/leo/scratch/era5/odbs/1759',
      'era5_1761': '/mnt/users/scratch/leo/scratch/era5/odbs/1761',
      
      'bufr': '/mnt/users/scratch/leo/scratch/era5/odbs/ai_bfr/',                                   
      'ncar': database_service2_path + 'UADB_20230109/', # check if extra '/' needed 

      #'igra2': database_service2_path + 'IGRA2_20230106',
      'igra2': database_service2_path + 'IGRA2_03012024',
                  
      'era5_1_mobile': '/mnt/users/scratch/leo/scratch/era5/odbs/1_mobile' ,
      'era5_2_mobile': '/mnt/users/scratch/leo/scratch/era5/odbs/2',
      
      'amma': database_service2_path + 'AMMA_BUFR/AMMA_split_csv/',
      
      # pre ERA5 3188
      'giub': database_service2_path + "GIUB_07072023",
      
      # ARCTIC datasets
      'hara': database_service2_path + 'HARA-NSIDC-0008_csv/',
      'npsound' : database_service2_path + 'NPSOUND-NSIDC0060/NPSOUND_csv_converted',
      'shipsound' : database_service2_path + 'SHIPSOUND-NSIDC-0054/' ,
      
      # Intercomparison
      'mauritius': database_service2_path + 'mauritius_2005/',
      'mauritius_digitized': database_service2_path + 'MAURITIUS/',
      'yangjiang' : database_service2_path + 'YANGJIANG/intercomparison_yangjiang_2010/'
} 


 ### Datasets to be provessed (list of datasets, e.g. datasets = ['era5_1', 'era5_2'] )
 # all_ds = list ( db.keys() )
datasets_big = ['era5_1', 'era5_2',
                'era5_3188', 'era5_1759', 'era5_1761',
                'ncar',
                'igra2',
                'bufr' ,
                'giub',
                'amma',
                'npsound', 'shipsound', 'hara']

datasets = datasets_big

datasets = ['npsound','hara','shipsound']

#datasets = ['npsound','shipsound','hara']


"""
era5_2 srvx1
era5_1 srvx8
"""

### output directory 

#out_dir = '/scratch/das/federico/HARVEST_YEARLY_16JAN2024/'

out_dir = '/scratch/das/federico/HARVEST_YEARLY_16JAN2024_full_harvest/'

### single_station - set to a primary id to run all the files mapping to the given primary_id e.g.   station = '0-20001-0-10393'
#stations = ['0-20001-0-10393' , '0-20001-0-11035' , '0-20000-0-70219' ]
stations = False

### Number of processes per dataset 
### WARNING: it multiplies the number of datasets chosen! 
processes = 5

### Pre-selection of already existing fully harvested files ( only checked by the running script)
skip_fully_harvested = True

### Only process missing stations (all years)
run_only_missing_stations = False

### Only process missing years per station 
check_missing_year = False

### Time range for harvesting
min_year_to_process = 1880 
max_year_to_process = 2024

