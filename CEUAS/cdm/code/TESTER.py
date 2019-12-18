import os,sys
""" This script runs 7 processes in parallel to convert the listed files. 
    Ueeful to create the single netCDF file for a single station, from all the datasets  """


""" Name of the output folder """
O = 'Big_Fixed'

""" Dictionary of files """
files_big  =  { 'era5_1'    :'/raid60/scratch/leo/scratch/era5/odbs/1/era5.conv._47646'          , 
            'era5_1759' :'/raid60/scratch/leo/scratch/era5/odbs/1759/era5.1759.conv.1:47646' ,
            'era5_1761' :'/raid60/scratch/leo/scratch/era5/odbs/1761/era5.1761.conv.1:47646' ,
            'era5_3188' :'/raid60/scratch/leo/scratch/era5/odbs/3188/era5.3188.conv.C:5357'  ,
            'igra2'     :'/raid60/scratch/federico/databases/IGRAv2//JAM00047646-data.txt'   ,
            'ncar'      :'/raid60/scratch/federico/databases/UADB//uadb_windc_47646.txt'     , 
            'bufr'      :'/raid60/scratch/leo/scratch/era5/odbs/ai_bfr/era5.47646.bfr'       ,
            'ncar'      :'/raid60/scratch/federico/databases/UADB//uadb_trhc_47646.txt'  }



files_small =  { 'era5_1'    :'/raid60/scratch/leo/scratch/era5/odbs/1/era5.conv._82930'     ,
            'era5_1759' :'/raid60/scratch/leo/scratch/era5/odbs/1759/era5.1759.conv.1:82930' ,
            'igra2'     :'/raid60/scratch/federico/databases/IGRAv2/BRM00082930-data.txt'   ,
            'ncar'      :'/raid60/scratch/federico/databases/UADB//uadb_windc_82930.txt'     ,
            'bufr'      :'/raid60/scratch/leo/scratch/era5/odbs/ai_bfr/era5.82930.bfr'       ,
}


era5 = { 'era5_1'    :'/raid60/scratch/leo/scratch/era5/odbs/1/era5.conv._47646'   }



""" Run the parallel processes """
for k,v in files_small.items():
    os.system('/opt/anaconda3/bin/python3 build_311c_cdmfiles_ALL_split.py -f ' + v + ' -d ' + k + ' -o ' + O + '    &   >> ciao.txt' )
        
