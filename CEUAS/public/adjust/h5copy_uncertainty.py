import h5py
import os,glob
import numpy
import xarray as xr
from multiprocessing import Pool

#fns=glob.glob(os.path.expandvars('$RSCRATCH/converted_v5/old/*.nc'))
fns=glob.glob(os.path.expandvars('/data/public/adj/*.nc'))

print('input list length:',len(fns))
ll=0

def do_copy(fn):
    
#    fno=fn.split('.nc')[0]+'copy'+'.nc'
    fno=''.join(fn.split('old/'))
    fno=os.path.expandvars('/data/private/new/'+''.join(fn.split('/')[-1].split('_uncertainty')))
    k='advanced_uncertainty'
    #with xr.open_dataset(fno,k) as ds:
        #print(ds)
    with h5py.File(fn,'r') as f:
        try:
            
            with h5py.File(fno,'r+') as g:
        
                #f.copy(f[k],g,expand_refs=True, expand_external=True,expand_soft=True)
                if k not in g.keys():
                    
                    g.create_group(k)
                #for v in f[k].keys():
                for v in f[k].keys():
                    if 'num' not in v and 'desroziers' not in v:
                        if v=='index':
                            try:
                                
                                del g[k][v]
                            except:
                                pass
                            try:
                                
                                g[k].create_dataset_like(v,f[k][v])
                                g[k][v][:]=f[k][v][:]
                            except:
                                pass
                        continue
                    
                    try:
                        try:
                            del g[k][v]
                        except Exception as e:
                            #print(e)
                            pass
                        if 'num' in v:
                            g[k].create_dataset_like(v,f[k][v],dtype=numpy.int32,compression='gzip')
                            g[k][v][:]=f[k][v][:]
                        else:
                            g[k].create_dataset_like(v,f[k][v],compression='gzip')
                            g[k][v][:]=f[k][v][:]
                    except Exception as e:
                        print(e)
                        continue
                    #ll+=1
                    for a,av in f[k][v].attrs.items():
                        if a not in ('CLASS','NAME','REFERENCE_LIST','DIMENSION_LIST'):
                            #print(a,av)
                            g[k][v].attrs[a]=av
                
                for v in g[k].keys(): #var_selection:
                    l=0            
                    try:
                        fvv=g[k][v]
                        if 'string' not in v and v!='index':                    
                            g[k][v].dims[l].attach_scale(g[k]['index'])
                            #print(v,fvv.ndim,type(fvv[0]))
                            if fvv.ndim==2 or type(fvv[0]) in [str,bytes,numpy.bytes_]:
                                slen=sdict[v]
                                #slen=10
                                g[k][v].dims[1].attach_scale(g[k]['string{}'.format(slen)])
                    except Exception as e:
                        print(fn.split('/')[-1],e)
                        pass
                
        except Exception as e:
            print(fn.split('/')[-1],e)
            pass
        
        print(fn.split('/')[-1]+' copied')
                
                
                
            #for v in slist:
                #s='string{}'.format(v)
                #for a in ['NAME']:
                    #g[k][s].attrs[a]=numpy.bytes_('This is a netCDF dimension but not a netCDF variable.')
        
        
    #f.close()
    #g.close()
    
    #ds=xr.open_dataset(fno,k)
    #print(ds)

p=Pool(20)
list(p.map(do_copy,fns))
    
print('ready')