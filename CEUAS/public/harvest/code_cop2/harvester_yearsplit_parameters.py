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
      'ncar': database_service2_path + 'UADB_22012024/',

      #'igra2': database_service2_path + 'IGRA2_20230106',
      'igra2': database_service2_path + 'IGRA2_03012024',
                  
      'era5_1_mobile': '/mnt/users/scratch/leo/scratch/era5/odbs/1_mobile' ,
      'era5_2_mobile': '/mnt/users/scratch/leo/scratch/era5/odbs/2',
      
      'amma': database_service2_path + 'AMMA_BUFR/AMMA_split_csv_22FEB2024/', #/scratch/das/federico/databases_service2/AMMA_BUFR/AMMA_split_csv_22FEB2024
      
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


 ### Datasets to be processed (list of datasets, e.g. datasets = ['era5_1', 'era5_2'] )
datasets_big = ['era5_1', 'era5_2',
                #'era5_3188', 
                'era5_1759', 'era5_1761',
                'ncar',
                'igra2',
                #'bufr' ,
                'giub',
                'amma',
                'npsound', 'shipsound', 'hara']




#datasets = datasets_big
#out_dir = '/scratch/das/federico/HARVEST_YEARLY_16JAN2024_full_harvest/'
#station_kind = 'regular'



#datasets = datasets_big
#out_dir = '/scratch/das/federico/HARVEST_YEARLY_16JAN2024_full_harvest_orphan/'
#station_kind = 'orphan'


#datasets = ['era5_2' , ]
#out_dir = '/scratch/das/federico/HARVEST_YEARLY_16JAN2024_full_harvest/'
#station_kind = 'regular'


datasets = ['era5_1_mobile' , 'era5_2_mobile']
datasets = ['era5_2_mobile' ]

datasets = ['hara' ]

out_dir = '/scratch/das/federico/HARVEST_YEARLY_22FEB2024_amma/'
station_kind = 'regular'

### Select station kind [regular, orphan, mobile]

#station_kind = 'regular'
#station_kind = 'orphan'
#station_kind = 'mobile'

### single_station - set to a primary id to run all the files mapping to the given primary_id e.g.   station = '0-20001-0-10393'
#stations = ['0-20001-0-10393' , '0-20001-0-11035' , '0-20000-0-70219' ]
stations = False

### Number of processes per dataset 
### WARNING: it multiplies the number of datasets chosen! 
processes = 40

### Pre-selection of already existing fully harvested files ( only checked by the running script)
skip_fully_harvested = True

### Only process missing stations (all years)
run_only_missing_stations = False

### Only process missing years per station 
check_missing_year = False

### Time range for harvesting
min_year_to_process = 2016
max_year_to_process = 2017

