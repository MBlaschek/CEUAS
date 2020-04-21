""" Utility to test the conversion of source files to CDM compliant netCDF files.
    We provide one group of files, one for each dataset, as an example. 
    Will run a separate process for each file listed in the dictionaries below.
    Please select an output directory name to store the output files.
"""
import os,sys

""" Example files, containing a small amount of data, for the same observation station. 
     Each key of the dictionary is the proper name of the dataset, while each value is the complete
     path to the file. In this case, they are found in the data/example_stations directory.
     
     Note that the the ncar dataset requires two distinct file:
     ncar (wind) : containing the full wind  data
     ncar_t (temperature): containing the full temperature data
"""

""" To test stations divided into multiple files, i.e. these files contain data for the same station from the 3188 dataset. """ 
multiple = {         'era5_3188_a' : '../data/example_stations/era5_3188/era5.3188.conv.C:4567' ,
                     'era5_3188_b' : '../data/example_stations/era5_3188/era5.3188.conv.C:5246' ,
                     'era5_3188_c' : '../data/example_stations/era5_3188/era5.3188.conv.C:4629' , }



small = { 'era5_1'    : '../data/example_stations/era5_1/era5.conv._82930.gz'           ,
          'era5_1759' : '../data/example_stations/era5_1759/era5.1759.conv._1:82930.gz' ,
          'ncar_t'    : '../data/example_stations/ncar/uadb_trhc_82930.txt'             ,
          'ncar'      : '../data/example_stations/ncar/uadb_windc_82930.txt'            ,
          'igra2'     : '../data/example_stations/igra2/BRM00082930-data.txt'           ,
          'bufr'      : '../data/example_stations/bufr/era5.82930.bfr'                  , 
          'era5_3188' : '../data/example_stations/era5_3188/era5.3188.conv._C:5246.gz'  ,

 }



michi = { 
          'ncar_t'    : '/raid60/scratch/federico/databases/UADB/uadb_trhc_53845.txt',
          'ncar'      :'/raid60/scratch/federico/databases/UADB/uadb_windc_53845.txt',

          'igra2'  : '/raid60/scratch/federico/databases/IGRAv2/CHM00053845-data.txt',

          'bufr'   : '/raid60/scratch/leo/scratch/era5/odbs/ai_bfr/era5.53845.bfr',  }




#0-20000-0-53845_era5_1_harvested_era5.conv.53845.txt.gz.nc


""" Select here the group of files """ 
selected = small

""" Select the output directory """
out_dir = 'output_test'


if __name__ == '__main__':

    for k,v in selected.items():
        k = k.replace('ncar_t','ncar').replace('_a','').replace('_b','').replace('_c','')
        print(k , '   ', v )
        os.system('/opt/anaconda3/bin/python3 harvest_convert_to_netCDF_newfixes_removeduplicates.py -f ' + v + ' -d ' + k + ' -o ' + out_dir + ' & ' )
        





    os.system('rm -r    /raid8/srvx1/federico/GitHub/CEUAS_master_FEB202/CEUAS/CEUAS/public/merge/example_stations' )
    os.system('cp -r output_test /raid8/srvx1/federico/GitHub/CEUAS_master_FEB202/CEUAS/CEUAS/public/merge/example_stations' )
    print(' Done copying!!!')










