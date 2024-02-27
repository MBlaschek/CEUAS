### 
### Merging parameters
###

### Path to harvested file directory
harvested_base_dir = '/scratch/das/federico/HARVEST_YEARLY_16JAN2024_full_harvest'

### Output merged directory (create if does not exist)
merged_out_dir = '/scratch/das/federico/MERGED_YEARLY_14FEB2024'
merged_out_dir = '/scratch/das/federico/MERGED_YEARLY_20FEB2024'

merged_out_dir = '/scratch/das/federico/MERGED_YEARLY_26FEB2024_TEST_AMMA'

### Kind of stations to merge ['regular' , 'orphans', 'mobile'] (if set to 'mobile', will run only mobile stations, otherwise will run all types )
station_kind = 'regular'
    
### harvested directories
harvested_base_dir = '/scratch/das/federico/HARVEST_YEARLY_16JAN2024_full_harvest'
if station_kind in ['orphan' , 'mobile']:
    harvested_base_dir = harvested_base_dir + '_' + station_kind 
    
# can be changed singularly if needed 
data_directories   = { 'era5_1'       : harvested_base_dir + '/era5_1' ,
                       
                                   'era5_1_mobile'       : harvested_base_dir + '/era5_1_mobile' ,
                                   'era5_2_mobile'       : harvested_base_dir + '/era5_2_mobile' ,
                                   
                                   'era5_2'       : harvested_base_dir + '/era5_2' ,
                                   'era5_3188' : harvested_base_dir + '/era5_3188' ,
                                   'era5_1759' : harvested_base_dir + '/era5_1759' ,
                                   'era5_1761' : harvested_base_dir + '/era5_1761' ,
                                   'ncar'           : harvested_base_dir + '/ncar' ,
                                   'igra2'          : harvested_base_dir + '/igra2' ,
                                   'bufr'            : harvested_base_dir + '/bufr' , 
                                   #'amma'        : harvested_base_dir + '/amma' ,
                                   
                                   'amma': '/scratch/das/federico/HARVEST_YEARLY_22FEB2024_amma/amma', 
                                   
                                   'giub'           : harvested_base_dir + '/giub' ,

                                   'hara'           : harvested_base_dir + '/hara' ,
                                   
                                   }



### Creating merging list from the available files in the harvested directory. If FALSE, will read available file 
create_merging_list = False

### Runing exception
run_exception = False

### Adding sensor id
add_sensor = True

### stations_list (optional)
stations = False

### Check if merged file exixts, and remove it from list
check_missing_stations = True 

### Multiprocesses run
POOL = True
pool_number = 20

