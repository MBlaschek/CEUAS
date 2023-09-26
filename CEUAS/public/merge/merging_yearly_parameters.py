### Merging parameters


### Path to harvested file directory
harvested_base_dir = '/scratch/das/federico/HARVEST_YEARLY_18SEP2023/'

# can be changed singularly if needed 
data_directories   = { 'era5_1'       : harvested_base_dir + '/era5_1' ,
                                   'era5_1_mobile'       : harvested_base_dir + '/era5_1_mobile' ,
                                   #'era5_2'       : harvested_base_dir + '/era5_2' ,
                                   'era5_3188' : harvested_base_dir + '/era5_3188' ,
                                   'era5_1759' : harvested_base_dir + '/era5_1759' ,
                                   'era5_1761' : harvested_base_dir + '/era5_1761' ,
                                   'ncar'           : harvested_base_dir + '/ncar' ,
                                   'igra2'          : harvested_base_dir + '/igra2' ,
                                   'bufr'            : harvested_base_dir + '/bufr' , 
                                   'amma'        : harvested_base_dir + '/amma' ,
                                   }





### Output merged directory (create if does not exist)
merged_out_dir = '/scratch/das/federico/MERGED_YEARLY_20SEPT2023'

### Years to be merged; default: all available years 
year = []

###
run_exception = False

