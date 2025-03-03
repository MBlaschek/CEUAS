import os, sys, glob, shutil
import h5py
import numpy as np
import matplotlib.pylab as plt

def dim_attach(g, k):
    for v in g[k].keys(): #var_selection:
        l=0            
        try:
            fvv=g[k][v]
            if 'string' not in v and v not in ('index', 'hdrlen'):
                if 'index' in g[k].keys():
                    
                    fvv.dims[l].attach_scale(g[k]['index'])
                else:
                    fvv.dims[l].attach_scale(g[k]['hdrlen'])
                    
                #print(v,fvv.ndim,type(fvv[0]))
                if fvv.ndim==2 : #or type(fvv[0]) in [str,bytes,np.bytes_]:
                    #print(type(fvv[0]), fvv[0])
                    slen=fvv.shape[1] #sdict[v]
                    #slen=10
                    g[k][v].dims[1].attach_scale(g[k]['string{}'.format(slen)])
        except MemoryError as e:
            print(g.filename.split('/')[-1],k, e)
            pass

def fix(indir, outdir = '/mnt/users/scratch/leo/scratch/fixedgiub2'): 
    #outdir = '/mnt/users/scratch/leo/scratch/fixedgiub'
    iprob = 0
    #indir='/mnt/scratch/scratch/federico/HARVEST_YEARLY_16JAN2024_full_harvest/giub'
    #indir = '/users/staff/uvoggenberger/scratch/CUON_HARVEST/harvest_regular/giub'
    for fn in glob.glob(indir+'/*/*.nc'):
        #if '1956' not in fn:
            #continue
        fo = '/'.join(fn.split('/')[-2:])
        os.makedirs(os.path.dirname(outdir+'/'+fo), exist_ok=True)
        try:
            os.remove(os.path.dirname(outdir+'/'+fo))
        except:
            pass
        #shutil.copyfile(fn, outdir+'/'+fo)
        ri = [] # recordindex
        rts = [] #recordtimestamp
        tsl = [] # time since launch
        iiprob = 0
        with h5py.File(fn, 'r') as f:
            ts = f['observations_table']['date_time'][:]
            z = f['observations_table']['z_coordinate'][:]
            zt = f['observations_table']['z_coordinate_type'][:]
            ot = f['observations_table']['observed_variable'][:]
            rio = f['recordindex'][:]
            ztu = np.unique(zt)
            if len(ztu) == 1:
                if ztu[0] == 0:
                    
                    ri = np.concatenate((np.array([0]), np.where(((ts[1:]-ts[:-1])>=1800)|((z[1:]-z[:-1])<0))[0]+1))
                else:
                    ri = np.concatenate((np.array([0]), np.where(((ts[1:]-ts[:-1])>=1800)|((z[1:]-z[:-1])>0))[0]+1))
                    
            rts = ts[ri]
            tsl = np.zeros_like(ts)
            repeat = False
            for i in range(len(ri)-1):
                tsl[ri[i]:ri[i+1]] = ts[ri[i]:ri[i+1]] - rts[i]
                zvar = z[ri[i]:ri[i+1]]+ot[ri[i]:ri[i+1]] / 1000.
                zu, zi = np.unique(zvar, return_index=True)
                #zvar = list(zip(z[ri[i]:ri[i+1]], ot[ri[i]:ri[i+1]]))
                #zu,zi = np.unique(zvar, return_index=True)
                #print(len(zu), len(zvar))
                if len(zu) < len(zvar):
                    repeat = True
                    #print(os.path.basename(fn), 'problematic', i, np.max(ts[ri[i]:ri[i+1]]-ts[ri[i]]), len(ts))
                    iiprob += 1
                    if iiprob > 1:
                        iprob += 1
                        break
            if iiprob > 1:
                
                print(os.path.basename(fn), 'problematic', i, np.max(ts[ri[i]:ri[i+1]]-ts[ri[i]]), len(ts))
        #continue
        if not repeat:
            shutil.copyfile(fn, outdir+'/'+fo)
            if len(rio) != len(ri) or np.any((rio-ri)!=0):
            
                with h5py.File(outdir+'/'+fo, 'r+') as f:
                    ts = f['observations_table']['date_time'][:]
                    z = f['observations_table']['z_coordinate'][:]
                    zt = f['observations_table']['z_coordinate_type'][:]
                    ot = f['observations_table']['observed_variable'][:]
                    rio = f['recordindex'][:]
                    ri = np.concatenate((np.array([0]), np.where((ts[1:]-ts[:-1])>3600)[0]+1))
                    rts = ts[ri]
                    tsl = np.zeros_like(ts)
                    del f['recordindex']
                    del f['recordtimestamp']
                    
                    for i in range(len(ri)-1):
                        tsl[ri[i]:ri[i+1]] = ts[ri[i]:ri[i+1]] - rts[i]
                        ts[ri[i]:ri[i+1]] = rts[i]
                    tsl[ri[-1]:] = ts[ri[-1]:] - rts[-1]
                    ts[ri[-1]:] = rts[-1]
                    f['observations_table']['date_time'][:] = ts
                    #f['observations_table'].create_dataset('timed', data=tsl)
                    
                    f.create_dataset('recordindex', data=ri)
                    f.create_dataset('recordtimestamp', data=rts)
                    hi = np.searchsorted(f['header_table']['record_timestamp'][:], rts)
                    if np.std(rts - f['header_table']['record_timestamp'][:][hi]) != 0.:
                        raise ValueError(f.filename)
                    for v in f['header_table'].keys():
                        if v == 'hdrlen':
                            del f['header_table'][v]
                            f['header_table'].create_dataset('hdrlen', data=np.zeros(ri.shape[0], dtype='S1'))
                        elif 'string' not in v:
                            tmp = f['header_table'][v][:][hi]
                            del f['header_table'][v]
                            f['header_table'].create_dataset(v, data=tmp)
                            
                            
            
                print(os.path.basename(fn), 'modified')
                with h5py.File(outdir+'/'+fo, 'r+') as f:
                    f.create_dataset('index', data=np.zeros(ri.shape[0], dtype='S1'))
                    f['recordindex'].dims[0].attach_scale(f['index'])
                    f['recordtimestamp'].dims[0].attach_scale(f['index'])
                    #dim_attach(f, '/')
                    dim_attach(f, 'observations_table')
                    #dim_attach(f, 'header_table') # fails for some reason
                    
                print(os.path.basename(fn), 'modified')
            else:
                print(os.path.basename(fn), 'copied')
               
    print('problematic', iprob)        
 
 
 