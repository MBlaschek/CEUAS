### Merging parameters


### Path to harvested file directory
harvested_base_dir ='/scratch/das/federico/HARVEST_YEARLY_20OCT2023/'

harvested_base_dir = '/scratch/das/federico/HARVEST_YEARLY_10NOV2023_Vienna_Lin_Bethel/'


# can be changed singularly if needed 
data_directories   = { 'era5_1'       : harvested_base_dir + '/era5_1' ,
                                   'era5_1_mobile'       : harvested_base_dir + '/era5_1_mobile' ,
                                   'era5_2'       : harvested_base_dir + '/era5_2' ,
                                   'era5_3188' : harvested_base_dir + '/era5_3188' ,
                                   'era5_1759' : harvested_base_dir + '/era5_1759' ,
                                   'era5_1761' : harvested_base_dir + '/era5_1761' ,
                                   'ncar'           : harvested_base_dir + '/ncar' ,
                                   'igra2'          : harvested_base_dir + '/igra2' ,
                                   'bufr'            : harvested_base_dir + '/bufr' , 
                                   'amma'        : harvested_base_dir + '/amma' ,
                                   'giub'        : harvested_base_dir + '/giub' ,
                                   
                                   }

### Output merged directory (create if does not exist)
merged_out_dir = '/scratch/das/federico/MERGED_YEARLY_13NOV2023_ViennaCOMPLETE'

merged_out_dir = '/scratch/das/federico/MERGED_YEARLY_20NOV_full'

merged_out_dir = '/scratch/das/federico/MERGED_YEARLY_20NOV_checkPlots'

merged_out_dir = '/scratch/das/federico/MERGED_YEARLY_20NOV_FULL_VIENNA'


merged_out_dir = '/scratch/das/federico/MERGED_YEARLY_20NOV_FULL_VIENNA_newplots'


merged_out_dir = '/scratch/das/federico/MERGED_YEARLY_22NOV_FULL_checkDimensions'

### Years to be merged; default: all available years 
year = []

### Runing exception
run_exception = False

### Multiprocesses run
POOL = True
pool_number = 10
