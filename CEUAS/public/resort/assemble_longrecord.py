
import hdf5plugin
import h5py
import os,glob
import numpy as np
import xarray as xr
import time
import matplotlib.pylab as plt
from datetime import datetime, timedelta
from numba import njit

from functools import partial
#from multiprocessing import Pool
import ray
from harvest_convert_to_netCDF import write_dict_h5

def write_var(k,fnsgood,fno,vi):
    

    ll=0
    tt=time.time()
    try:
        
        with h5py.File(fno,'r+') as g:
        
            for fn in fnsgood[:]:
            #    fno=fn.split('.nc')[0]+'copy'+'.nc'
                fno=''.join(fn.split('old/'))
                #with xr.open_dataset(fno,k) as ds:
                    #print(ds)
                with h5py.File(fn,'r') as f:
                    try:
                        
                        #for v in f[k].keys():
                        ub = f[k][vi].shape[0]
                        if k == 'recordindices' and vi != 'recordtimestamp':
                            ub -= 1
                        for v in [vi]:
                            if ll==0: 
                                #g[k].create_dataset_like(v,f[k][v],shape=(nptlens[-1],))
                                g[k][v][:nptlens[0]]=f[k][v][:ub]
                            else:
                                try:
                                    
                                    g[k][v][nptlens[ll-1]:nptlens[ll]]=f[k][v][:ub]
                                except:
                                    continue
                            ll+=1
                            print(ll,fnl,k, v, fn.split('/')[-1]+' copied',time.time()-tt)
        
                    except Exception as e:
                        print(fn.split('/')[-1],e)
                        pass
            print('x')
    except Exception as e:
        print('open',fn.split('/')[-1],e)
        pass
    return
        
# resorting the data
#
@njit
def make_vrindex(vridx,ridx): # this function is similar to np.unique with return_index=True, but it expands the index to the  
    # original array dimensions
    l=0
    for i in range(len(ridx)): # to set the recordindices
        if i == 0:
            l +=1
        else:
            if ridx[i]>ridx[i-1]:
                vridx[ridx[i-1] + 1 : ridx[i] + 1]=l # next record after l
                l += 1
            else:
                l += 1
    vridx[ridx[i] + 1 :]=l #len(idx) # next record for the last element is the len of the data

def pwrite_variable(fno, fps, rdict, ridict, tup):       

    tt = time.time()
    gname, vname = tup
    if vname == 'index' or 'string' in vname:
        return
    if vname == 'fg_depar@body':
        print(vname)
    rstarts = {s: 0 for s in sorted(list(ridict.keys())[:-2])}
    csum = 0
    for s in rstarts.keys():
        rstarts[s] = csum
#        csum += rsave[s]
        csum += np.sum([rdict[s][i].stop-rdict[s][i].start for i in range(len(rdict[s])) if rdict[s][i] is not None ])

    var = None #np.zeros(csum, dtype=int)
    l = 0
    for k, rs in sorted(rstarts.items()):
        irel = 0
        #for ifn in range(len(fns[:])):
            #with h5py.File(fns[ifn],'r') as f:
        for ifn in range(len(fps)):
                f = fps[ifn]
            
                if var is None:
                    if f[gname][vname].ndim == 1:
                        
                        var = np.empty_like(f[gname][vname][:], shape=(csum, ))
                        if var.dtype == np.dtype('int32') or var.dtype == np.dtype('int64'):
                            var[:] = -2147483648
                        else:
                            var[:] = np.nan
                    else:
                        if vname == 'sensor_id':                    
                            var = np.empty_like(f[gname][vname][:], shape=(csum, 4))
                        else:
                            var = np.empty_like(f[gname][vname][:], shape=(csum, f[gname][vname].shape[1]))
                        var[:] = b''
                        
                if rdict[k][ifn]:   #not none
                    try:
                        
                        chunk = f[gname][vname][rdict[k][ifn]]
                    except:
                        print('ERROR ', fno)
                    #print(f['observations_table']['observed_variable'][rdict[k][ifn].stop])
                    if vname == 'sensor_id':
                        
                        var[rs+irel:rs+irel+chunk.shape[0], :chunk.shape[1]] = chunk
                    else:
                        try:
                            
                            var[rs+irel:rs+irel+chunk.shape[0]] = chunk
                        except:
                            if chunk.shape[1] < var.shape[1]:
                                print( 'warning, var wider than chunk ...', fps[ifn], vname, k, var[0], chunk[0])
                                var[rs+irel:rs+irel+chunk.shape[0], :chunk.shape[1]] = chunk
                            else:                                
                                print( 'warning, cut chunk in width ...', fps[ifn], vname, k, var[0], chunk[0])
                                var[rs+irel:rs+irel+chunk.shape[0]] = chunk[:, :var.shape[1]]
                                print('cut chunk succeeded')

                    irel += chunk.shape[0] 
                    #if np.any(f['observations_table']['observed_variable'][rdict[k][ifn]]!=k):
                        #print(fns[ifn], k, ifn, f['observations_table']['observed_variable'][rdict[k][ifn]])
                    #if np.any(var[rs+irel-chunk.shape[0]:rs+irel]!=\
                              #h[gname][vname][rs+irel-chunk.shape[0]:rs+irel]):
                        #print(fns[ifn], k, ifn, h['recordindices'][str(k)][[0, -1]])
                    #else:
                        #l += 1
                        
    #for i in fbkeys:
        #print(i)
        #print(time.time()-tt)

        #if i in reduced_fbkeys:
            #continue
        #else: 
            #with eua.CDMDataset(fn) as data:
                #rest_data = data.era5fb[i][:]
            #if i in ['expver', 'source@hdr', 'source_id', 'statid@hdr']:
                #final = np.empty((addedvar[-1][1],len(rest_data[0])), dtype=rest_data[0].dtype)
            #else:
                #final = np.empty(addedvar[-1][1], dtype=rest_data[0].dtype)
            #ov_vars = fill_restdata(final, rest_data, addedvar, jj)

        #ov_vars = ov_vars[absidx]

        #if i == 'index':
            #pass
        #elif i in ['expver', 'source@hdr', 'source_id', 'statid@hdr']:
            #alldict = {i:np.asarray(ov_vars, dtype='S1')}
            #write_dict_h5(targetfile, alldict, 'era5fb', {i: { 'compression': 'gzip' } }, [i])
        #else:
            #alldict = pandas.DataFrame({i:ov_vars})
            #write_dict_h5(targetfile, alldict, 'era5fb', {i: { 'compression': 'gzip' } }, [i]) 
    fnov = fno[:-3]+'.'+vname+'@'+gname+'.nc'
    print(fnov.split('/')[-1], time.time() - tt)
    if var.ndim == 2:
        alldict = {vname:np.asarray(var, dtype='S1')}
#        write_dict_h5(fnov, alldict, gname, {vname: { 'compression': 'gzip' , 'chunksizes': ( min( [10000,var.shape[0] ] ), var.shape[1] ) } } , [vname], mode='w')
        write_dict_h5(fnov, alldict, gname, {vname: { 'compression': 'gzip' } } , [vname], mode='w')
    else:
        #alldict = pandas.DataFrame({i:ov_vars})
        alldict = {vname:var}
        write_dict_h5(fnov, alldict, gname, {vname: { 'compression': 'gzip' } }, [vname], mode='w')
    print(gname, vname, time.time() - tt)
    
ray_pwrite_variable = ray.remote(pwrite_variable)

def dim_attach(g, k):
    for v in g[k].keys(): #var_selection:
        l=0            
        try:
            fvv=g[k][v]
            if 'string' not in v and v!='index':                    
                g[k][v].dims[l].attach_scale(g[k]['index'])
                #print(v,fvv.ndim,type(fvv[0]))
                if fvv.ndim==2 or type(fvv[0]) in [str,bytes,np.bytes_]:
                    slen=fvv.shape[1] #sdict[v]
                    #slen=10
                    g[k][v].dims[1].attach_scale(g[k]['string{}'.format(slen)])
        except MemoryError as e:
            print(g.filename.split('/')[-1],k, e)
            pass

def h5concatenate(fkey, h=None) :
    
    tt = time.time()
    fns=glob.glob(os.path.expandvars(fkey))
    fns.sort()
    fnl=len(fns)
    
    dirs = fns[0].split('/')
    dirs[-2] = 'long'
    fno=os.path.expandvars('/'.join(dirs))
    if os.path.exists(fno+'.txt'):
        print(fno, 'already exists')
        wtime = os.path.getmtime(fno+'.txt')
        rtimes = np.array([os.path.getmtime(fn) for fn in fns])
        if np.any(rtimes>wtime):
            print(fno, 'older than some input')
        else:
            return
    
    vdict={'observations_table': ['observed_variable','observation_value', 'z_coordinate','date_time', 'index'],
            'recordindices': ['126', 'index', 'recordtimestamp']}
    vset = set()
    fps = []
    try:
        
        for fn in fns[:]:
            fps.append(h5py.File(fn,'r'))
            f = fps[-1]
    #        with h5py.File(fn,'r') as f:
            vset = vset.union(set(list(f['recordindices'].keys())[:-2]))
        sh = f['recordindices']['recordtimestamp'].shape[0] # just to check if recordtimestamp exists
    except Exception as e:
        print(fn, e)
        return
    ivset = np.array(list(vset), dtype=int)
    ivset.sort()
    rdict = {ivset[i]: [] for i in range(len(ivset))}
    rdict2 = {ivset[i]: [] for i in range(len(ivset))}
    ridict = {ivset[i]: [] for i in range(len(ivset))}
    rsave = {}
    rsave2 = {}
    risave = {}
    rts = []
    ris = []
    isum = 0
    #for fn in fns[:]:
        #with h5py.File(fn,'r') as f:
    for f in fps:
        #with h5py.File(fn,'r') as f:
        try:
            
            for i in ivset:
                rdict[i].append(None) 
                rdict2[i].append(None) 
                ridict[i].append(None)

# just for checking
            #vs, vidx= np.unique(f['observations_table']['observed_variable'][:], return_index=True)
            #for i in range(len(vs)):
                #if i < len(vs) - 1:
                    #rdict[vs[i]][-1] = slice(vidx[i], vidx[i+1])
                #else:
                    #rdict[vs[i]][-1] = slice(vidx[i], f['observations_table']['observed_variable'].shape[0])
                
                #if vs[i] in rsave.keys() :
                    
                    #rsave[vs[i]] += rdict[vs[i]][-1].stop - rdict[vs[i]][-1].start
                #else:
                    #rsave[vs[i]] = rdict[vs[i]][-1].stop - rdict[vs[i]][-1].start
                    
                    
            rts.append(f['recordindices']['recordtimestamp'][:])
            ris.append(f['recordindices']['index'][:-1])
            ivs = sorted(np.array(list(f['recordindices'].keys())[:-2], dtype=int))
            for i in ivs:
                ridict[i][-1]= f['recordindices'][str(i)][:]
                rdict[i][-1] = slice(ridict[i][-1][0], ridict[i][-1][-1])
                #rdict2[i][-1] = slice(ridict[i][-1][0], ridict[i][-1][-1])
                #if np.any(rdict2[i][-1] != rdict[i][-1]):
                    #print('x')

        except Exception as e:
            if(os.path.exists('err.log')):
                mode = 'a'
            else:
                mode = 'w'
            with open('err.log', mode) as fe:
                fe.write(f.filename+', no recordtimestamp')
                print(f.filename+' no recordtimestamp', e)
            return
                
            #print(fn.split('/')[-2])
        
    csum =0
    rt = np.zeros(np.sum([i.shape[0] for i in rts]), dtype=int)
    ri ={}#s: np.zeros(np.sum([ridict[s][l].shape[0] for l in ridict[s]]), dtype=int) for s in ridict.keys()}
    l = 0
    ival = 0
    for k, rs in ridict.items():
        
        csum = 0
        for s in rs:
            if s is not None:
                csum += s.shape[0]
        ri[k] = np.zeros(rt.shape[0]+1, dtype=int)
        irel = 0
        l = 0
        for i in range(len(rs)):
#            with h5py.File(fns[ifn],'r') as f:
                if rs[i] is not None:   #not none
                    chunk = rs[i] #[rdict[k][ifn]]
                    ri[k][irel:irel+chunk.shape[0]] = ival + chunk - chunk[0]
                    irel += chunk.shape[0] - 1
                    ival += chunk[-1]- chunk[0]
                else:
                    if irel == 0:
                        ri[k][irel:irel+rts[i].shape[0]] = ival
                    else:
                        ri[k][irel:irel+rts[i].shape[0]] = ival
                    irel += rts[i].shape[0]
                if h is not None:
                    if np.any(h['recordindices'][str(k)][:irel]!=ri[k][:irel]):
                        print(k,l,irel,ival, h['recordindices'][str(k)][:irel],ri[k][:irel])
                        raise AssertionError
                        #plt.plot(h['recordindices'][str(k)][:irel]-ri[k][:irel], label='h')
                        ##plt.plot(ri[k][:irel], label='ri')
                        #plt.title(str((k, l, irel)))
                        #plt.legend()
                        #plt.show()
                
                l+=1
        if ri[k][-1] == 0:
            ri[k][-1] =ri[k][-2]    

        
    ridict['recordtimestamp'] = np.concatenate(rts)
    ridict['index'] = np.concatenate(ris+[np.array([b''], dtype='S1')])
    
    if h is not None:
        for k in h['recordindices'].keys():
            if k not in ('index', 'recordtimestamp'):          
                print(k, np.where(h['recordindices'][k][:]-ri[int(k)][:]!=0)[0])
            else:
                print(k, np.where(h['recordindices'][k][:]!=ridict[k][:])[0])
            
    #rstarts = {s: 0 for s in sorted(rsave.keys())} #rvars[i]: np.sum(rdict[rvars[i]]) for i in range(len(rvars))}
    #s = list(rsave.keys())
    with h5py.File(fn, 'r') as g:
        
        vars = {'observations_table':list(g['observations_table'].keys()) ,
                'era5fb': list(g['era5fb'].keys()),}
    #vlist = [(k, i) for i in vars[k] for k in vars.keys()]
    gvlist = []
    for k,v in vars.items():
        for i in v:
            gvlist.append((k, i))
    
    print(fno, 'before write', time.time()-tt)
    #func = partial(pwrite_variable,fno, fns, rdict, ridict )
    func = partial(pwrite_variable,fno, fps, rdict, ridict )
    list(map(func, gvlist))
    for f in fps:
        f.close()
    #futures = [ray_pwrite_variable.remote( fno, fns, rdict, ridict, gv) for gv in gvlist]
    #obj_ref = ray.get(futures)
    
    print(fno, 'before copy', time.time() - tt)
    try:       
        os.remove(fno+'.txt') # this ensures that the file fno is regenerated in case something went wrong during write.
    except:
        pass
    # copy variable files to final file
    with h5py.File(fno, 'w') as g:
        for gv in gvlist:
            k, v = gv
            if(k not in g.keys()):                    
                g.create_group(k)
            rfile = fno[:-2] + '@'.join((v, k)) + '.nc'
            if v != 'index' and 'string' not in v:
                with h5py.File(rfile, 'r') as f:
                    for i in f[k].keys():
                        if i not in g[k]:
                            f.copy(k+'/'+i, g[k], i, without_attrs=True)
    
                            for a,av in f[k][i].attrs.items():
                                if a not in ('CLASS','NAME','REFERENCE_LIST','DIMENSION_LIST'):
                                    print(a,av)       
                                    g[k][i].attrs[a]=av
        for k in g.keys():
            dim_attach(g, k)
        
        for gr in ('header_table', 'source_configuration'):
            
            g.create_group(gr)
            start = 0
            noindex = True
            for fn in fns:
                with h5py.File(fn, 'r') as f:
                    for k in f[gr].keys():
                        fhk = f[gr][k]
                        if k not in g[gr].keys():
                            
                            if fhk.ndim == 1:
                                if 'string' not in k:
                                    g[gr].create_dataset_like(k, fhk, shape=ridict['recordtimestamp'].shape)
                                else:
                                    g[gr].create_dataset_like(k, fhk)
                                    
                            else:
                                g[gr].create_dataset_like(k, fhk, shape=(ridict['recordtimestamp'].shape[0], fhk.shape[1]))
                                
                        
                        if 'string' not in k:
                            stop = start + fhk[:].shape[0]
                            try:
                                
                                g[gr][k][start:stop] = fhk[:]
                            except:
                                print(start, stop, g[gr][k].shape)
                        if k == 'index':
                            noindex = False
                        
                start = stop
                
            if gr == 'header_table':
                # to satisfy need for second primary key relation to station_configuration table
                g[gr]['station_record_number'][:] = np.arange(ridict['recordtimestamp'].shape[0]) 
                # to fix missing or zero report_timestamp
                repts = g[gr]['report_timestamp'][:]
                invalid = np.where(repts<=0)
                if len(invalid[0]) > 0:
                    repts[invalid] = g[gr]['record_timestamp'][:][invalid]
                    print(fno, 'fixed', len(invalid[0]), 'timestamps')
                
                
            if noindex:
                #grlist = list(g[gr].keys())
                g[gr].create_dataset('index', data=np.full(ridict['recordtimestamp'].shape, b'', dtype='S1'))
                    
            dim_attach(g, gr)
        
                        
        rstarts = {s: 0 for s in sorted(list(ridict.keys())[:-2])}
        csum = 0
        for s in rstarts.keys():
            rstarts[s] = csum
    #        csum += rsave[s]
            csum += np.sum([rdict[s][i].stop-rdict[s][i].start for i in range(len(rdict[s])) if rdict[s][i] is not None ])
        
        with h5py.File(fn, 'r') as f:
            
            g.create_group('recordindices')
            var = None #np.zeros(csum, dtype=int)
            l = 0
#            rs = 0
            for k, rs in sorted(rstarts.items()):
                irel = 0
                var = np.full(ridict['index'].shape,-2147483648, dtype=np.int64)
                for ifn in range(len(ridict[k])):
                                
                        if rdict[k][ifn]:   #not none
                            csize = len(ridict[k][ifn])
                            #print(f['observations_table']['observed_variable'][rdict[k][ifn].stop])
                            var[irel:irel+csize] = rs + ridict[k][ifn] - ridict[k][ifn][0]
                            irel += csize - 1
                            rs += ridict[k][ifn][-1]- ridict[k][ifn][0]
                        else:
                            csize = len(rts[ifn])+1
                            var[irel:irel+csize] = rs
                            irel += csize - 1
                            

                g['recordindices'].create_dataset(str(k), data=var[:])
            g['recordindices'].create_dataset('index', data=ridict['index'])
            g['recordindices'].create_dataset('recordtimestamp', data=ridict['recordtimestamp'])
                
            dim_attach(g, 'recordindices')
            
            for gr in f.keys():
                if gr in ('dateindex', 'recordindices', 'header_table', 'source_configuration'):
                    continue                    
                if gr not in g.keys():
                    f.copy(gr, g, gr, without_attrs=True)
                    if gr == 'station_configuration' and 'record_number' not in g[gr].keys():
                        g[gr].create_dataset('record_number', data=np.arange(ridict['recordtimestamp'].shape[0]), dtype=ridict['recordtimestamp'].dtype)
                    grlist = list(f[gr].keys())
                    if 'index' not in grlist:
                        g[gr].create_dataset('index', data=np.array([b'']*f[gr][grlist[0]].shape[0], dtype='S1'))
                    
                dim_attach(g, gr)
        print('')
                    
           
        
    if h is not None:
        with h5py.File(fno, 'r') as g:
            for gr in h.keys():
                if gr in 'dateindex':
                    continue
                for k in h[gr].keys():
                    if k not in ['0', '34', '39', '106', '107', '126', '139', 'observation_value', 'fg_depar@body', 'date_time', 'report_id']:
                        continue
                    if h[gr][k].dtype != np.dtype('S1'):
                        
                        idx = np.where(np.logical_and(~np.isnan(h[gr][k]), ~np.isnan(g[gr][k])))[0]
                        if idx.shape[0] == 0:
                            print(gr, k)
                            continue
                    else:
                        idx = np.arange(h[gr][k].shape[0])
                    if idx[-1] +1 > g[gr][k].shape[0] or  idx[-1] +1 > h[gr][k].shape[0]:
                        print(g, k, 'x')
                    ggri = g[gr][k][:][idx]
                    hgri = h[gr][k][:][idx]
                    if g[gr][k].dtype == np.dtype('S1'):
                        crit = np.any(ggri != hgri)
                    else:
                        crit = np.any(np.abs(ggri - hgri) > 1.e-30)
                    if crit:
                        idy = np.where(ggri!=hgri)[0]
                        if g[gr][k].dtype != np.dtype('S1'):
                            iargmax = np.argmax(np.abs(ggri[idy] - hgri[idy]))
                            eps = 1.e-4
                            if np.abs(ggri[idy[iargmax]] - hgri[idy[iargmax]]) > eps * np.abs(ggri[idy[iargmax]]) or \
                               np.abs(ggri[idy[iargmax]] - hgri[idy[iargmax]]) > eps * np.abs(hgri[idy[iargmax]]):
                                print(os.path.basename(g.filename), gr, k, 'not equal', idx[idy[iargmax]], ggri[idy[iargmax]], hgri[idy[iargmax]])
                                plt.plot(g['observations_table']['date_time'][:][idx[idy]]/365.25/86400, np.array((ggri[idy[:]], hgri[idy[:]])).T)
                                plt.title(','.join((os.path.basename(fn), gr, k, str(g['observations_table']['observed_variable'][idx[idy[iargmax]]]))))
                                plt.savefig('_'.join((fn, gr, k, '.png')))
                                plt.close()
                        else:
                            print(os.path.basename(g.filename), gr, k, 'not equal', idx[idy[:2]], ggri[idy[:2]], hgri[idy[:2]])
                        x = 0
                    print(gr, k)
    
    vars =glob.glob(fno[:-2]+'*.nc')
    for v in vars:
        os.remove(v)
    with open(fno+'.txt', 'w') as f:
        f.write('done')
    print(fno, 'after copy', time.time() - tt)
        
    return

ray_h5concatenate = ray.remote(h5concatenate)

if __name__ == "__main__":

    #fns=glob.glob(os.path.expandvars('$RSCRATCH/converted_v5/old/*.nc'))
    #fkeys = os.path.expandvars('$RSCRATCH/converted_v10/[12]???/*26460*.nc')
    pwd = os.getcwd()
    os.chdir(os.path.expandvars('$RSCRATCH/converted_v10/'))
    longlist = glob.glob('[12]???/*11035*v1.nc')
    os.chdir(pwd)
    shlist = []
    shlist = [os.path.basename(ll) for ll in longlist]
    fkeys = np.unique(shlist)
    ray.init(num_cpus=20)
        
    for fks in fkeys:
        fkey = os.path.expandvars('$RSCRATCH/converted_v10/[12]???/') + fks
        try:
            
            fnc = glob.glob('/'.join(fkey.split('/')[:-2]) + '/' + fkey.split('/')[-1])[0]
        except:
            print(fkey, 'no match found')
            continue
        tt = time.time()
        h = h5py.File(fnc,'r')
        rvars = list(h['recordindices'].keys())
        #for i in range(len(rvars)):
            #if rvars[i].isdecimal():
                #rvars[i] = int(rvars[i])
        #rvars.sort()
        rdict = {rvars[i]: [] for i in range(len(rvars))}
        ridict = {rvars[i]: [] for i in range(len(rvars))}
        risave = {}
        
        h5concatenate(fkey)
    
    print('total', time.time() - tt)
    exit()
    func = partial(pwrite_variable, fno, fns, ridict)
    list(map(func, gvlist)) #'observations_table', 'observed_variable')      
    
    
    #hrsave = 0
    #for ik in irvars[:1]:
        #try:
            #k = str(ik)
            #hr = np.concatenate(ridict[k])
            #if len(hr) > 0:
                
                #plt.plot(h['recordindices'][k][:hr.shape[0]]-hr, label=k)
                #hrsave += hr[0]
        #except MemoryError:
            #pass
    #plt.legend()
    #plt.show()
    with h5py.File(fno,'w') as g:
        for k, vlist in vdict.items():
            
            tlens=[]
            fnsgood=[]
            for fn in fns[:]:
                try:
                    
                    with h5py.File(fn,'r') as f:
                        if k == 'recordindices':
                            tlens.append(f[k][vlist[0]].shape[0]-1)
                        else:
                            tlens.append(f[k][vlist[0]].shape[0])
                            
                        #for v in vlist:
                            #dum=f[k][v].shape[0]
                        fnsgood.append(fn)
                except:
                    continue
                    
            nptlens=np.cumsum(tlens,dtype=np.int)
    
        #with h5py.File(fno,'w') as g:
            #f.copy(f[k],g,expand_refs=True, expand_external=True,expand_soft=True)
            with h5py.File(fnsgood[0],'r') as f:
                if k not in g.keys():
                    
                    g.create_group(k)
                    for v in vlist:
                        if v not in g[k].keys():
                            g[k].create_dataset_like(v,f[k][v],shape=(nptlens[-1],))
            g.attrs['jumbo']='Y'
        
     
            func=partial(write_var,k,fnsgood,fno)
            dum=list(map(func, vlist))
        
        #with h5py.File(fno,'r+') as g:
            with h5py.File(fnsgood[-1],'r') as f:
                for a,av in f[k][v].attrs.items():
                    if a not in ('CLASS','NAME','REFERENCE_LIST','DIMENSION_LIST'):
                        print(a,av)
                        g[k][v].attrs[a]=av
            
                for v in g[k].keys(): #var_selection:
                    l=0            
                    try:
                        fvv=g[k][v]
                        if 'string' not in v and v!='index':                    
                            g[k][v].dims[l].attach_scale(g[k]['index'])
                            #print(v,fvv.ndim,type(fvv[0]))
                            if fvv.ndim==2 or type(fvv[0]) in [str,bytes,np.bytes_]:
                                slen=sdict[v]
                                #slen=10
                                g[k][v].dims[1].attach_scale(g[k]['string{}'.format(slen)])
                    except MemoryError as e:
                        print(fn.split('/')[-1],e)
                        pass
            
                   
                    
    print('ready', time.time-tt)