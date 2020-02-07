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
                     'era5_3188_b' : './data/example_stations/era5_3188/era5.3188.conv.C:5246' ,
                     'era5_3188_c' : './data/example_stations/era5_3188/era5.3188.conv.C:4629' , }


small = { 'era5_1'    : '../data/example_stations/era5_1/era5.conv._82930'           , 
          'era5_1759' : '../data/example_stations/era5_1759/era5.1759.conv.1:63260'  ,
          'era5_1761' : '../data/example_stations/era5_1761/era5.1761.conv.2:32904'  ,

          'igra2'     :'../data/example_stations/igra2/BRM00082930-data.txt'         ,
          
          'ncar'      :'../data/example_stations/ncar/uadb_windc_82930.txt'          ,
          'ncar_t'    : '../data/example_stations/ncar/uadb_trhc_82930.txt'          ,

          'bufr'      :'/raid60/scratch/leo/scratch/era5/odbs/ai_bfr/era5.82930.bfr' , }




""" Select here the group of files """ 
selected = small

""" Select the output directory """
out_dir = 'output_test'


if __name__ == '__main__':

    for k,v in selected.items():
        k = k.replace('ncar_t','ncar').replace('_a','').replace('_b','').replace('_c','')
        print(k , '   ', v )
        os.system('/opt/anaconda3/bin/python3 harvest_convert_to_netCDF.py -f ' + v + ' -d ' + k + ' -o ' + out_dir + ' & ' )
        
















