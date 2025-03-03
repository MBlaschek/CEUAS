import glob
import os
### 
### Merging parameters
###

### Path to harvested file directory
#harvested_base_dir = '/scratch/das/federico/HARVEST_YEARLY_16JAN2024_full_harvest'

### Output merged directory (create if does not exist)
#merged_out_dir = '/scratch/das/federico/MERGED_YEARLY_06MAR2024_TEST'
merged_out_dir = '/mnt/users/scratch/leo/scratch/UH/MERGED_YEARLY_14FEB25'
merged_out_dir = '/mnt/scratch/scratch/leo/scratch/UH/MERGED_YEARLY_19FEB25'

### Kind of stations to merge ['regular' , 'orphans', 'mobile'] (if set to 'mobile', will run only mobile stations, otherwise will run all types )
station_kind = 'regular'
    
### harvested directories
#harvested_base_dir = '/scratch/das/federico/HARVEST_YEARLY_16JAN2024_full_harvest'
#harvested_base_dir = '/mnt/users/scratch/leo/scratch/FH'
#harvested_base_dir = '/mnt/users//scratch/uvoggenberger/CUON_HARVEST/harvest_test_rerun'
#harvested_base_dir = '/mnt/users/scratch/leo/scratch/UH/CUON_HARVEST2/harvest'
harvested_base_dir = '/mnt/users/scratch/leo/scratch/UH/CUON_HARVEST6/harvest'
if station_kind in ['regular', 'orphan' , 'mobile']:
    harvested_base_dir = harvested_base_dir + '_' + station_kind 

scpath = './'
#scpath = '/srvfs/home/uvoggenberger/CEUAS/CEUAS/public/harvest/data/station_configurations'
scpath = os.path.expandvars('$HOME/CEUAS/CEUAS/meta/inventory_comparison_2/code/station_configuration/')
# can be changed singularly if needed 
data_directories   = { 'era5_1'       : harvested_base_dir + '/era5_1' ,
                       
                                   'era5_1_mobile'       : harvested_base_dir + '/era5_1_mobile' ,
                                   'era5_2_mobile'       : harvested_base_dir + '/era5_2_mobile' ,
                                   'npsound': harvested_base_dir + '/npsound', 
                                   'shipsound': harvested_base_dir + '/shipsound', 
                                   
                                   'era5_2'       : harvested_base_dir + '/era5_2' ,
                                   #'era5_3188' : harvested_base_dir + '/era5_3188' ,
                                   'era5_1759' : harvested_base_dir + '/era5_1759' ,
                                   'era5_1761' : harvested_base_dir + '/era5_1761' ,
                                   'ncar'           : harvested_base_dir + '/ncar' ,
                                   'igra2'          : harvested_base_dir + '/igra2' ,
                                   'bufr'            : harvested_base_dir + '/bufr' , 
                                   'bufr_cnr'            : harvested_base_dir + '/bufr_cnr' , 
                                   'amma'        : harvested_base_dir + '/amma' ,
                                   'maestro'        : harvested_base_dir + '/maestro' ,

                                   #'giub'           : harvested_base_dir + '/giub' ,
                                   'giub'           : '/mnt/users/scratch/leo/scratch/fixedgiub' ,

                                   'hara'           : harvested_base_dir + '/hara' ,
                                   'woudc'           : harvested_base_dir + '/woudc' ,
                                   #'era5_2_restored'       : harvested_base_dir + '_orphan/era5_2/restored' ,
                                   #'era5_1_restored'       : harvested_base_dir + '_orphan/era5_1/restored' ,
                                   #'era5_1761_restored'       : harvested_base_dir + '_orphan/era5_1761/restored' ,
                                   #'era5_1759_restored'       : harvested_base_dir + '_orphan/era5_1759/restored' ,
                                   
                                   }

restored_base = '/mnt/users/scratch/leo/scratch/UH/CUON_HARVEST6/harvest_orphan/'
if station_kind in ['orphan']:
    data_directories = {}
    dnames = [os.path.basename(d) for d in glob.glob(restored_base+'*')]
    for d in dnames:
        data_directories[d+'_reduced'] = restored_base + d + '/reduced/'
        data_directories[d+'_unchanged'] = restored_base + d + '/unchanged/'
        
else:
    
    kl = list(data_directories.keys())
    for k in kl:
        if 'mobile' not in k and 'giub' not in k:
            data_directories[k+'_restored'] = restored_base + k + '/restored'
            data_directories[k+'_reduced'] = restored_base + k + '/reduced'

### Creating merging list from the available files in the harvested directory. If FALSE, will read available file 
create_merging_list = True

### Runing exception
run_exception = False

### Adding sensor id
add_sensor = True

### stations_list (optional)
stations = False

### Check if merged file exixts, and remove it from list
check_missing_stations = False 

### Multiprocesses run
#POOL = True
pool_number = 128

