# Modul Laden
import sys,os
sys.path.append(os.getcwd()+'/../../cds-backend/code/')
import cds_eua3 as eua
import xarray as xr
import numpy as np
import datetime
from numba import njit
import matplotlib.pylab as plt
import glob
from multiprocessing import Pool
import time
import h5py

@njit
def add_biasestimate(xyzv,xyzt,press,atime0,adj,adjpress):
    adjv=np.empty_like(xyzv)
    adjv.fill(np.nan)
    isec18=64800
    isec06=21600
    for ip in range(adjpress.shape[0]):
        idtold=0
        for it in range(atime0.shape[0]):
            if it==atime0.shape[0]-1:
                idt=press.shape[0]
            else:
                idt=np.searchsorted(xyzt,atime0[it+1])
            for idx in range(idtold,idt):
                if press[idx]==adjpress[ip]:
                    ih=xyzt[idx]%86400
                    if ih>=isec18 or ih<isec06:
                        adjv[idx]=adj[0,ip,it]
                    else:
                        adjv[idx]=adj[1,ip,it]
            idtold=idt
            

    return adjv

def add_adj(adjfile):
    
    statid=adjfile[-8:-3]
    adjustments={'raobcore':os.path.expandvars(adjfile),
                 'rich':os.path.expandvars('corrsave_rio24_'.join(adjfile.split('corrsave'))),
                 'rase':os.path.expandvars(('ERA5bc_RAOBCORE_v'+version+'_').join(adjfile.split('feedbackglobbincorrsave'))),
                 'rise':os.path.expandvars(('ERA5bc_RAOBCORE_v'+version+'_').join(adjfile.split('feedbackglobbincorrsave')))}
    adjname={'raobcore':'rasocorr',
                 'rich':'rasocorr',
                 'rase':'bias',
                 'rise':'richbias'}
    try:
        
        merged='feedbackmerged'.join(adjfile.split('feedbackglobbincorrsave'))
        with h5py.File(merged) as f:
            uid=f.attrs['unique_source_identifier'].decode()
        #ifile=glob.glob('/raid60/scratch/leo/scratch/converted_v5/*'+statid+'_CEUAS_merged_v1.nc')[0]
        try:
            
            ifile=glob.glob('/raid60/scratch/leo/scratch/converted_v7/'+uid+'_CEUAS_merged_v1.nc')[0]
        except:
            uid='2000?-?-'.join(uid.split(uid[2:10]))
            print(uid, 'chose more liberal uid pattern, should not be necessary')
            ifile=glob.glob('/raid60/scratch/leo/scratch/converted_v7/'+uid+'_CEUAS_merged_v1.nc')[0]
        data = eua.CDMDataset(ifile)
        #idx=np.where(np.logical_and(data.observations_table.z_coordinate[:]==10000,data.observations_table.observed_variable[:]==85))
        #plt.plot(data.observations_table.date_time[:][idx],data.era5fb['biascorr@body'][:][idx])        
        #plt.plot(data.observations_table.date_time[:][idx],data.adjust['RASE_bias_estimate'][:][idx])        
        #plt.plot(data.observations_table.date_time[:][idx],data.adjust['RISE_bias_estimate'][:][idx])        
        if 'advanced_homogenisation' in data.groups: 
            print('already has bias estimate')
            #return
    except Exception as e:
        print(e)
        return
    
    # Daten einlesen Temperatur
    try:
        
        xyz = data.read_observed_variable(85, return_xarray=True,date_time_in_seconds=True)
    except:
        print (statid, 'could not read temperature')
        return
    
    ref=np.datetime64(datetime.datetime(1900,1,1),'ns')
    xyzt=(xyz.date_time.values-ref).astype('long')//1000000000
    

    for k,v in adjustments.items():
        try:
            
            adjustments=xr.open_dataset(v,decode_times=False)
            if adjustments.datum.ndim==2:
                atime0=(adjustments.datum[0].values.astype(int)-1)*86400.
            else:
                atime0=(adjustments.datum.values.astype(int)-1)*86400.
           
            
            mask=adjustments[adjname[k]].values==-999.
            adjustments[adjname[k]].values[mask]=np.nan
            tt=time.time()
            adj=add_biasestimate(xyz.values,xyzt,xyz.z_coordinate.values,atime0,
                                 adjustments[adjname[k]].values,adjustments.press.values*100)
            print('add:',time.time()-tt)
            xyz.values=adj
            
            #idx=np.where(xyz.z_coordinate.values==50000)
            #plt.plot(xyzt[idx]/86400/365.25,adj[idx])
            #plt.plot(atime0/86400/365.25,adjustments.rasocorr.values[0,11,:])
            #plt.plot(atime0/86400/365.25,adjustments.rasocorr.values[1,11,:])
            #plt.show()
            # Daten schreiben neue Variable monkey in neuer gruppe adjust
            data.write_observed_data(k.upper()+'_bias_estimate',
                                     ragged=xyz,  # input data
                                     varnum=85,  # observed_variable to be aligned with
                                     group='advanced_homogenisation',   # name of the new group
                                     data_time='date_time',  # named datetime coordinate
                                     data_plevs='z_coordinate',  # named pressure coordinate
                                     attributes={'version':version}
                                    )
            print('write:',time.time()-tt)
        except Exception as e:
            print(v,e)
            
                        
version='1.8'                        
#mfiles=glob.glob(os.path.expandvars('$FSCRATCH/rise/1.0/exp06/*/feedbackglobbincorrsave*.nc'))
mfiles=glob.glob('/raid60/raid/home/srvx7/lehre/users/a1400070/CEUAS/CEUAS/public/adjust/Temperature_adjustment/*/feedbackglobbincorrsave*.nc')
#with h5py.File(os.path.expandvars('$FSCRATCH/rise/1.0/exp06/108087/feedbackmerged108087.nc')) as f:
    #print(f.attrs.keys()
P=Pool(40)
tt=time.time()
list(P.map(add_adj,mfiles))
print(time.time()-tt)
print('ready')
