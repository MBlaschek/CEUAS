import hdf5plugin
import h5py
import numpy as np
import os,glob,sys
import datetime

miss=glob.glob(os.path.expandvars('$RSCRATCH/converted_v11/long/*v1.nc'))
nmiss=[]
for m in miss:
    key=m.split('_CEUAS')[0].split('-')[-1]
    #print(key)
    #print(m)
    fns=[m]#glob.glob(os.path.expandvars('$RSCRATCH/converted_v11/long/*'+m+'*v1.nc'))
    #fns=glob.glob(os.path.expandvars('/mnt/scratch/scratch/federico/MERGED_FEB*/'+os.path.basename(m)))
    #fns=glob.glob(os.path.expandvars('/mnt/scratch/scratch/federico/COP2_HARVEST_JAN2023/*/*'+key+'*.nc'))
    #print(fns)

    ref=datetime.datetime(1900,1,1)
    dstart=datetime.datetime(1945,7,1)
    istart=(dstart-ref).total_seconds()
    for fn in fns:
        with h5py.File(fn,'r') as g:
                #print(g.filename,g.keys())
                idx=np.where(np.abs(g['recordindices/recordtimestamp'][:]-istart)<8640000)[0]
                if len(idx)>0:
                    try:
                        ridxu=g['recordindices']['139'][idx[0]:idx[0]+2]
                        ridxv=g['recordindices']['140'][idx[0]:idx[0]+2]
                        print('/'.join(m.split('/')[-2:]),ridxu)
                        if ridxu[0]==ridxu[-1]:
                            continue
                        idy=np.arange(len(ridxu),dtype=int)
                        if 'wind_bias_estimate' not in g['advanced_homogenisation'].keys():
                            continue
                        if np.abs( g['advanced_homogenisation/wind_bias_estimate'][ridxu[0]:ridxu[-1]][idy[0]])>0.1 or np.abs( g['advanced_homogenisation/wind_bias_estimate'][ridxv[0]:ridxv[-1]][idy[0]])>0.1:
                            print(ref+datetime.timedelta(seconds=int(g['observations_table/date_time'][ridxu[0]])),
                                  list(zip(g['observations_table/z_coordinate'][ridxu[0]:ridxu[-1]][idy],
                                      g['observations_table/z_coordinate_type'][ridxu[0]:ridxu[-1]][idy],
                                      g['observations_table/observation_value'][ridxu[0]:ridxu[-1]][idy],
                                      g['observations_table/observation_value'][ridxv[0]:ridxv[-1]][idy],
                                      g['advanced_homogenisation/wind_bias_estimate'][ridxu[0]:ridxu[-1]][idy],
                                      g['advanced_homogenisation/wind_bias_estimate'][ridxv[0]:ridxv[-1]][idy])))
                                  
                        #if np.isnan(g['observations_table/z_coordinate'][ridx[0]:ridx[-1]][idy[0]]):
                            #nmiss.append(g.filename)
                            #break
                    except :
                        print('could not read',m)
        
print(nmiss)
