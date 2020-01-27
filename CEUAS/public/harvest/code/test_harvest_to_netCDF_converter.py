""" Utility to test the conversion of source files to CDM compliant netCDF files.
    We provide, as an example, three groups of file. 
    Will run a separate process for each file listed in the dictionaries below.
    Please select an output directory name to store the output files.
"""

import os,sys

out_dir = 'output_test'

""" Example of TATENO station, containing a large amount of data """
big  =  { 'era5_1'    :'/raid60/scratch/leo/scratch/era5/odbs/1/era5.conv._47646'          , 
            'era5_1759' :'/raid60/scratch/leo/scratch/era5/odbs/1759/era5.1759.conv.1:47646' ,
            'era5_1761' :'/raid60/scratch/leo/scratch/era5/odbs/1761/era5.1761.conv.1:47646' ,
            'era5_3188' :'/raid60/scratch/leo/scratch/era5/odbs/3188/era5.3188.conv.C:5357'  ,
            'igra2'     :'/raid60/scratch/federico/databases/IGRAv2//JAM00047646-data.txt'   ,
            'ncar'      :'/raid60/scratch/federico/databases/UADB//uadb_windc_47646.txt'     , 
            'bufr'      :'/raid60/scratch/leo/scratch/era5/odbs/ai_bfr/era5.47646.bfr'       ,
            'ncar'      :'/raid60/scratch/federico/databases/UADB//uadb_trhc_47646.txt'  }


""" Exmaple of xxx station, containing a small amount of data """
small =  {       'era5_1'    :'/raid60/scratch/leo/scratch/era5/odbs/1/era5.conv._82930'           ,
                 'era5_1759' :'/raid60/scratch/leo/scratch/era5/odbs/1759/era5.1759.conv.1:82930'  ,
                 'igra2'     :'/raid60/scratch/federico/databases/IGRAv2/BRM00082930-data.txt'     ,
                 'ncar'      :'/raid60/scratch/federico/databases/UADB//uadb_windc_82930.txt'      ,
                 'bufr'      :'/raid60/scratch/leo/scratch/era5/odbs/ai_bfr/era5.82930.bfr'        ,
                 'ncar_t'    : '/raid60/scratch/federico/databases/UADB/uadb_trhc_82930.txt'          }


""" To test stations divided into multiple files """ 
multiple = { 'era5_3188_a' : '/raid60/scratch/leo/scratch/era5/odbs/3188/era5.3188.conv.C:4567' ,
             'era5_3188_b' : '/raid60/scratch/leo/scratch/era5/odbs/3188/era5.3188.conv.C:5246' ,
             'era5_3188_c' : '/raid60/scratch/leo/scratch/era5/odbs/3188/era5.3188.conv.C:4629' , }


""" Select here the group of files """ 
selected = small

for k,v in selected.items():
    k = k.replace('ncar_t','ncar').replace('_a','').replace('_b','').replace('_c','')
    print(k , '   ', v )
    os.system('/opt/anaconda3/bin/python3 harvest_to_netCDF_converter_leo+federico.py -f ' + v + ' -d ' + k + ' -o ' + O + ' & ' )
        
















