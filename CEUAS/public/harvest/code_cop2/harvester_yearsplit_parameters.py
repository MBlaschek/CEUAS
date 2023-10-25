"""
parameters for harvester run
Will be copied inside the output directory so that the parameters will be stored

"""


### Datasets sources

datasets_path = { 'era5_1': '/mnt/users/scratch/leo/scratch/era5/odbs/1' ,
      'era5_2': '/mnt/users/scratch/leo/scratch/era5/odbs/2',
      'era5_3188': '/mnt/users/scratch/leo/scratch/era5/odbs/3188',
      'era5_1759': '/mnt/users/scratch/leo/scratch/era5/odbs/1759',
      'era5_1761': '/mnt/users/scratch/leo/scratch/era5/odbs/1761',
      
      'bufr': '/mnt/users/scratch/leo/scratch/era5/odbs/ai_bfr/',                                   
      'ncar': '/scratch/das/federico/databases_service2/UADB_20230109/', # check if extra '/' needed 

      'igra2': '/scratch/das/federico/databases_service2/IGRA2_20230106',

      'era5_1_mobile': '/mnt/users/scratch/leo/scratch/era5/odbs/1_mobile' ,
      'era5_2_mobile': '',
      
      'amma': '/scratch/das/federico/databases_service2/AMMA_BUFR/AMMA_split_csv/',
      
      # pre ERA5 3188
      'giub': "/scratch/das/federico/redo_GIUB_07072023_reduced",
      
      # ARCTIC datasets
      'hara':'/scratch/das/federico/databases_service2/HARA-NSIDC-0008_csv/',
      'npsound' : '/scratch/das/federico/databases_service2/NPSOUND-NSIDC0060/NPSOUND_csv_converted',
      'shipsound' : '/scratch/das/federico/databases_service2/SHIPSOUND-NSIDC-0054/' ,
      
      'mauritius': '/scratch/das/federico/databases_service2/mauritius_2005/'     
} 


 ### Datasets to be provessed (list of datasets, e.g. datasets = ['era5_1', 'era5_2'] )
 # all_ds = list ( db.keys() )
datasets = ['era5_1', 'era5_2', 'era5_3188', 'era5_1759', 'era5_1761', 'ncar', 'igra2', 'bufr' ]
 
datasets = [ 'era5_1759', 'igra2' ]
datasets = [ 'era5_1759', 'era5_1761' ]


datasets = [ 'ncar' ]
 
### Select output directory 
out_dir = '/scratch/das/federico/HARVEST_YEARLY_17AUG2023/'
 
### Number of processes per dataset 
processes = 4

### Only process missing stations (all years)
run_only_missing_stations = False

### Only process missing years per station 
check_missing_year = False
