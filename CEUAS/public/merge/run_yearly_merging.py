import os,sys


ranges= [ 
           [1880, 1965] , 
           [1966,1990] , 
           [1991,2000] , 
           [2001,2010] ,            
           [2011, 2017],
           [2018, 2019],
           [2020, 2021],
           [2022, 2022],
]


ranges= [ 
           [2001,2004] ,           
           [2005,2008] ,            
           [2009,2010] ,            
           
]

for r in ranges:
    print('Running range ::: ' , r )
    os.system('python  merging_cdm_netCDF_yearSplit_SEP2023.py  -min_y ' + str(r[0])  + '  -max_y  ' + str(r[1] ) + ' & ')
    
    
