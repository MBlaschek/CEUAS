import os,sys,glob
import numpy as np
import h5py

'''
This script is to be called after the feedbackmerged files have been created in the working directory.
It creates the file mergedstations.t which is needed to run the RAOBCORE/RICH Fortran executables. 
Please adjust paths as appropriate
'''
ipath='/raid60/raid/home/srvx7/lehre/users/a1400070/CEUAS/CEUAS/public/adjust/Temperature_adjustment'
os.chdir(os.path.expanduser('~/tmp'))
cpath=os.getcwd()
opath=cpath
os.chdir(ipath)
mlist=glob.glob('[0-9]*/feedbackmerged*.nc')

os.chdir(cpath)
i=1
with open('mergedstations.t','w') as f:
    for m in mlist:
        with h5py.File(ipath+'/'+m,'r') as hf:
            
            if  hf['time'].shape[0]>365:
                
                #print(hf.keys())
                try:
                    lat=hf['lat'][0]
                    lon=hf['lon'][0]
                except:
                    continue
                unique_source_identifier=hf.attrs['unique_source_identifier']
    
                f.write('{:5d}   0 '.format(i)+m[:6]+' {:6.2f} {:7.2f} {}\n'.format(lat,lon,unique_source_identifier.decode())) 
                i+=1


print('ready')