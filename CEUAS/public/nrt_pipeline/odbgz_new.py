import os,sys,glob
import subprocess,psutil
from multiprocessing import Pool
fns=[]
year='????'
month='??'
def odbmobile(fn):
    print(fn)
    os.system(f"odc sql --full-precision -q 'select *' -i "+f'"{fn}" | tr -d " "  | gzip   > "{fn}.gz"')
    
def odb(fn):
    print(fn)
    os.system(f"odc sql --full-precision -q 'select *' -i "+f'"{fn}" | tr -d " "  | gzip   > "{fn}.gz"')

era_5_1_dir = ''

for patt in f'{era_5_1_dir}/era5.conv.{year}{month}.[0-9]????',:
#for patt in f'era5.conv.{year}{month}.[A-Z]????',\
#    f'era5.conv.{year}{month}.????', f'era5.conv.{year}{month}.??????',\
#    f'era5.conv.{year}{month}.???????',  f'era5.conv.{year}{month}.????????':
    
    print(patt)
    fns+=glob.glob(patt)

#for fn in fns:
#    if 'DAV' in fn:
#        print(f'x{fn}x')
#        odb(fn)
#exit()
   
P=Pool(40)
x=list(P.map(odb,fns[:]))
