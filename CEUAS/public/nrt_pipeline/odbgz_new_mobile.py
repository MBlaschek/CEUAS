import os,sys,glob
import subprocess,psutil
from multiprocessing import Pool
fns=[]
year='????'
month='??'
def odbmobile(fn):
    print(fn)
#    os.system(f"mkdir /mnt/users/scratch/uvoggenberger/CUON_HARVEST_202503/data/era5_1_data/")
    os.system(f"odc sql --full-precision -q 'select *' -i "+f'"{fn}" | tr -d " "  | gzip   > "{fn}.gz"') # "{fn.replace('era5_1_data', 'era5_1_mobile_data')}.gz"')
    
def odb(fn):
    print(fn)
    os.system(f"odc sql --full-precision -q 'select *' -i "+f'"{fn}" | tr -d " "  | gzip   > "{fn}.gz"')

era_5_1_dir = ''

#for patt in f'/mnt/users/scratch/uvoggenberger/CUON_HARVEST_202503/data/era5_data/era5.conv.{year}{month}.[0-9]????',:
for patt in f'{era_5_1_dir}/era5.conv.{year}{month}.[A-Z]????',\
    f'{era_5_1_dir}/era5.conv.{year}{month}.????', f'{era_5_1_dir}/era5.conv.{year}{month}.??????',\
    f'{era_5_1_dir}/era5.conv.{year}{month}.???????',  f'{era_5_1_dir}/era5.conv.{year}{month}.????????':
#    if '.gz' in patt:
#        continue
    
    print(patt)
    fns+=glob.glob(patt)
    fns = [fn for fn in fns if not '.gz' in fn]
#for fn in fns:
#    if 'DAV' in fn:
#        print(f'x{fn}x')
#        odb(fn)
#exit()
   
P=Pool(40)
x=list(P.map(odbmobile,fns[:]))
