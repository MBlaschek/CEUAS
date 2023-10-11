import matplotlib.pylab as plt
from netCDF4 import Dataset
import h5py
import numpy as np
import pandas as pd

statid='091165'
ps=np.array((10,20,30,50,70,100,150,200,250,300,400,500,700,850,925,1000))*100.
weights=pd.read_fwf('/users/staff/leo/fastscratch/MSUDaten/std_atmosphere_wt_function_chan_4.txt',skiprows=1,widths=(5,10,10,10,10,10,10,10,10,10,10),error_bad_lines=False)
x=weights.iloc[:-3].to_numpy()
x=np.asarray(x,dtype='float')
pweight=ps.copy()
for ip in range(pweight.shape[0]):
    l=0
    while ps[ip]<x[l][3]:
        l+=1
    #l
    #if ps[ip]<x[3][l] and ps[ip]>x[3][l+1]:
    pweight[ip]=(x[l][5]*(ps[ip]-x[l+1][3])+x[l+1][5]*(x[l][3]-ps[ip]))/(x[l][3]-x[l+1][3])

pweight /= np.sum(pweight[1:])

with Dataset(statid+'/feedbackglobbincorrmon_bt2_'+statid+'.nc') as f:
    x=np.nanmean(f.variables['montemp'][:1,2,:],axis=0)
    #plt.plot(f.variables['datum'][0,:]/365.25+1900,f.variables['montemp'][0,2,:],label='RAOBCORE BT')

with Dataset(statid+'/feedbackglobbinmon_bt2_'+statid+'.nc') as f:
    plt.plot(f.variables['datum'][0,:]/365.25+1900,np.nanmean(f.variables['montemp'][:1,2,:],axis=0)-x,label='ADJ BT')

#with Dataset(statid+'/feedbackglobbincorrmon'+statid+'.nc') as f:
    #plt.plot(f.variables['datum'][0,:]/365.25+1900,f.variables['montemp'][0,5,:],label='temperature at 100hPa 00GMT')
#with Dataset(statid+'/feedbackglobbincorrmon'+statid+'.nc') as f:
    #plt.plot(f.variables['datum'][0,:]/365.25+1900,f.variables['montemp'][1,5,:],label='temperature at 100hPa 12GMT')
#with Dataset(statid+'/feedbackglobbincorrmon'+statid+'.nc') as f:
    #x=f.variables['rasocorrmon'][:,2,:]
    #plt.plot(f.variables['datum'][0,:]/365.25+1900,f.variables['rasocorrmon'][0,5,:],label='temperature at 100hPa 00GMT')
with Dataset(statid+'/feedbackglobbincorrmon'+statid+'.nc') as f:
    plt.plot(f.variables['datum'][0,:]/365.25+1900,np.nanmean(f.variables['rasocorrmon'][:1,4,:],axis=0),label='ADJ 100 hPa')
    x=f.variables['rasocorrmon'][:]
    bt=np.zeros_like(x,shape=(x.shape[0],x.shape[2]))
    for ip in range(1,10):
        bt+=pweight[ip]*x[:,ip,:]
    plt.plot(f.variables['datum'][0,:]/365.25+1900,np.nanmean(bt[:1,:],axis=0),label='ADJ WEIGHTED')
    plt.title(statid)
    
    
'''
with Dataset(statid+'/feedbackmerged'+statid+'.nc') as f:
    plt.plot(f.variables['datum'][:]/365.25+1900,f.variables['temperatures'][0,0,:])

with h5py.File('/mnt/users/scratch/leo/scratch/converted_v8/0-20001-0-01001_CEUAS_merged_v1.nc') as f:
    start=f['recordindices']['126'][0]
    stop=f['recordindices']['126'][-1]
    idx=np.where(f['observations_table']['z_coordinate'][start:stop]==1000)
    plt.plot(f['observations_table']['date_time'][start:stop][idx]/86400/365.25+1900,f['observations_table']['observation_value'][start:stop][idx])
'''    

plt.legend()
plt.savefig('/users/staff/leo/bttest.png')
plt.show()
