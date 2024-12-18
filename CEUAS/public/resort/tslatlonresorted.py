import os, glob, sys
#import hdf5plugin
import h5py
import numpy as np
import datetime
import time
#import ray
from numba import njit
import matplotlib.pylab as plt
import ray
import pickle
#os.chdir(os.path.expanduser('~/python/CEUAS/CEUAS/public/adjust'))
#sys.path.append(os.path.expanduser('~/python/CEUAS/CEUAS/public/common/Rasotools/rasotools/'))
#sys.path.append(os.path.expanduser('~/python/CEUAS/CEUAS/public/common/Rasotools/'))
from utils import tdist,extract,rmeanw
import copy
from tqdm import tqdm
import shutil

def fextractrs(fn):

    ref = datetime.datetime(1900, 1, 1)
    parlist = '106','107','126'
    dtemplate =  dict(id=[], lat=[], lon=[], ulatlon=[], rt=[])
    dtemplate0 =  dict(id=[], lat=[], lon=[], ulatlon=[], rt=[])
    for par in parlist:
        for prop in 'ws', 'zz', 'andep', 'rlen', 'pmin':         
            dtemplate[prop+par]=[np.nan]
            dtemplate0[prop+par]=[]

    daydict = {'fn': fn}
    plist = np.array([5000., 10000, 20000, 30000, 50000, 70000, 85000])
    ref = datetime.datetime(1900, 1, 1)
    rt = None
    try:
        
        with h5py.File(fn, 'r') as f:
            key = os.path.basename(fn.split('.nc')[0])
            #for k in range(46000*12):
                #daydict[k] = dict(sec=-3599+k*7200, id=[], lat=[], lon=[], ulatlon=[], rlen=[], pmin=[])
    
            p = 0
            rt = f['recordindices']['recordtimestamp'][:] 
            for k in range(len(rt)):
                dk = (rt[k] + 3600)//7200
                daydict[dk] = copy.deepcopy(dtemplate)
                    
            for par in parlist:
                    
                try:
                    ri = f['recordindices'][par][:]
                    tslice = slice(ri[0], ri[-1])
                    l = ri.shape[0] - 1
                    while ri[l] == ri[ri.shape[0] - 1] and l > 0:
                        l -= 1
                    ri = ri[:l+1] - ri[0]
                    
                    l=len(ri)
                    if l > 0:
                        
                        try:
                            
                            lat = f['observations_table']['latitude'][tslice][ri]
                        except:
                            print(fn, 'no index')
                            raise KeyError
                        lon = f['observations_table']['longitude'][tslice][ri]
                        z= f['observations_table']['z_coordinate'][tslice]
                        obs= f['observations_table']['observation_value'][tslice]
                        an_dep= f['era5fb']['an_depar@body'][tslice]
                        #obsvar = f['observations_table']['observed_variable'][:]
                        #f['observations_table']['longitude'][tslice].shape[0]
                        for k in range(len(ri)):
                            dk = (rt[k] + 3600)//7200
                            dd = daydict[dk]
                            if not dd['rt']:   

                                dd['rt'].append(rt[k])
                                dd['id'].append(key)
                                dd['lat'].append(lat[k])
                                dd['lon'].append(lon[k])
                                dd['ulatlon'].append(f'{lat[k]:.2f}_{lon[k]:.2f}')
                                
                            if k<l - 1:
                                dd['rlen'+par][0] = ri[k+1]-ri[k]
                            else:
                                dd['rlen'+par][0] = tslice.stop-tslice.start-ri[k]
                            zz = z[ri[k]:ri[k]+dd['rlen'+par][-1]]
                            if len(zz) > 0:                       
                                dd['pmin'+par][0] = np.nanmin(zz)
                            #else:
                                #dd['pmin'+par].append(np.nan)
                            ws = obs[ri[k]:ri[k]+dd['rlen'+par][-1]]
                            siglev = np.isin(zz, plist)
                            if np.sum(siglev) == 0:
                                siglev[:] = True
                            idx = np.where(siglev)[0]
                            dd['zz'+par][0] = zz[idx][:plist.shape[0]]
                            dd['ws'+par][0] = ws[idx][:plist.shape[0]]
                            dd['andep'+par][0] = np.nanstd(ws[idx][:plist.shape[0]])
                                    
                                
                    else:
                        print(ri)
                    print(key, par, ref+datetime.timedelta(seconds=int(rt[-1])), lat[-1], lon[-1], daydict[dk]['pmin'+par][-1])
                except KeyError:
                    if 'recordindices' in f.keys():
                        
                        print(key, 'no',par,'read failed', f['recordindices'].keys())
                    else:
                        print(key, 'no recordindices')
                    
    except Exception as e:
        print(fn, 'read failed', e)
    
    return daydict #ray.put(daydict)
        
    
ray_fextractrs = ray.remote(fextractrs)

def dim_attach(g, k):
    for v in g[k].keys(): #var_selection:
        l=0            
        try:
            fvv=g[k][v]
            if 'string' not in v and v!='index':                    
                g[k][v].dims[l].attach_scale(g[k]['index'])
                #print(v,fvv.ndim,type(fvv[0]))
                if len(fvv) > 0 and (fvv.ndim==2 or type(fvv[0]) in [str,bytes,np.bytes_]):
                    slen=fvv.shape[1] #sdict[v]
                    #slen=10
                    g[k][v].dims[1].attach_scale(g[k]['string{}'.format(slen)])
        except MemoryError as e:
            print(g.filename.split('/')[-1],k, e)
            pass

def writereduced(fn, daydict, ref):
    
    os.makedirs(os.path.dirname(fn)+'/reduced', exist_ok=True)
    for iters in range(3):
        
        try:
        
            with h5py.File(fn, 'r') as f:
                key = os.path.basename(fn.split('.nc')[0])
                #for k in range(46000*12):
                    #daydict[k] = dict(sec=-3599+k*7200, id=[], lat=[], lon=[], ulatlon=[], rlen=[], pmin=[])
        
                try:
                    
                    ri = f['recordindices']['107'][:]
                    tslice = slice(ri[0], ri[-1])
                    l = ri.shape[0] - 1
                    while ri[l] == ri[ri.shape[0] - 1] and l > 0:
                        l -= 1
                    #ri = ri[:l] - ri[0]
                    rt = f['recordindices']['recordtimestamp'][:]
                except:
                    return 0, 0
                
                keeplist = []
                i = 0
                #try:
                    
                for r in rt:
                    k = (r +3600) // 7200
                    try:
                        idx = daydict[k]['id'].index(key)
                        if not daydict[k]['dup'][idx]:
                            keeplist.append(i)
                    except Exception as e:
                        print(os.path.basename(fn), 'key', key, 'not found for', ref+datetime.timedelta(seconds=int(r)), ', but put to keeplist')
                        keeplist.append(i)
                        
                    i += 1
                #except:
                    #print(key, 'key not found!!!')
                    #return
                
                #print(fn, rt.shape[0], len(keeplist))
                if len(keeplist) ==0 or len(keeplist) / rt.shape[0] > 0.99:
                    #print(fn, rt.shape[0], len(keeplist))
                    if len(keeplist) == 0:
                        print(fn, 'full duplicate')
                        with h5py.File(os.path.dirname(fn)+'/reduced/'+os.path.basename(fn), 'w') as h:
                            h.attrs['status'] = 'full duplicate'
                    return rt.shape[0], len(keeplist)
                #else:
                    #print(fn, rt.shape[0], len(keeplist))
                    
                fslices = []
                nindex = 0
                parlist = []
                keylist = []
                for v in f['recordindices'].keys():
                    if v[0] in ('0','1','3'):
                        parlist.append(f['recordindices'][v][0])
                        keylist.append(v)
                pidx = np.argsort(parlist)
                keylist = np.array(keylist)[pidx]
                for v in keylist: #f['recordindices'].keys():
                    #if v[0] in ('0','1','3'):
                        wix = f['recordindices'][v][:]
                        for k in keeplist:
                            fslices.append(slice(wix[k], wix[k+1]))
                            nindex+=wix[k+1]-wix[k]
                
                os.makedirs(os.path.dirname(fn)+'/reduced', exist_ok=True)    
                with h5py.File(os.path.dirname(fn)+'/reduced/'+os.path.basename(fn), 'w') as h:
                    
                    for gr in ['era5fb', 'observations_table']:
                        g = h.create_group(gr)
                        for v in f[gr].keys():
                            if v == 'index':
                                g.create_dataset('index', data=np.empty(nindex, dtype='S1'), compression= 32015, compression_opts=(3, ))                       
                            elif 'string' in v:
                                xdata = f[gr][v][:]
                                g.create_dataset(v, data=xdata)
                            else:
                                shp = list(f[gr][v].shape)
                                shp[0] = nindex
                                if f[gr][v].chunks is None:
                                    chk = [0]
                                else:
                                    chk = list(f[gr][v].chunks)
                                chk[0] = np.min((chk[0], nindex))
                                chk = tuple(chk)
                                if chk[0] == 0:
                                    g.create_dataset_like(v, f[gr][v], shape=shp, chunks=None, compression= 32015, compression_opts=(3, ))
                                else:
                                    g.create_dataset_like(v, f[gr][v], shape=shp, chunks=chk, compression= 32015, compression_opts=(3, ))
                                vm = f[gr][v][:]
                                #gm = vm[:nindex]
                                #tt = time.time()
                                #i = 0
                                if(fslices):
                                    gm = np.concatenate([vm[s] for s in fslices])
                                    #for s in fslices:
                                        #ip = s.stop-s.start
                                        #gm[i:i+ip] = vm[s]
                                        #i += ip
                                    #print(time.time()-tt)
                                    g[v][:] = gm
                        for a,av in f[gr][v].attrs.items():
                            if a not in ('CLASS','NAME','REFERENCE_LIST','DIMENSION_LIST'):
                                #print(a,av)
                                h[gr][v].attrs[a]=av
                            
                        dim_attach(h, gr)
                        print(h.filename.split('/')[0], gr, 'written')
                        
                    for gr in ['recordindices']:
                        g = h.create_group(gr)
                        #for v in keylist: #f[gr].keys():
                            ##if v[0] in ('0', '1','3'):
                                #wix = f['recordindices'][v][:]
                                #for k in keeplist:
                                    #fslices.append(slice(wix[k], wix[k+1]))
                                    #nindex+=wix[k+1]-wix[k]
                        nindex = len(keeplist)
                        ref = 0
                        i = 0
                        g.create_dataset('index', data=np.empty(nindex+1, dtype='S1'), compression= 32015, compression_opts=(3, ))
                        v = 'recordtimestamp'
                        if f[gr][v].chunks is None:
                            g.create_dataset_like(v, f[gr][v], shape=(nindex, ), data=f[gr][v][:][keeplist])
                        else:
                            g.create_dataset_like(v, f[gr][v], shape=(nindex, ), data=f[gr][v][:][keeplist],
                                                    chunks=(np.min((nindex,f[gr][v].chunks[0])), ))
                        for v in keylist: #f[gr].keys():
                            #if v == 'index':
                                #g.create_dataset('index', data=np.empty(nindex, dtype='S1'), compression= 'gzip')                       
                            #elif v[0] in ('0','1','3'):
                            shp = list(f[gr][v].shape)
                            shp[0] = nindex + 1
                            if f[gr][v].chunks is None:
                                g.create_dataset_like(v, f[gr][v], shape=shp)
                            else:
                                g.create_dataset_like(v, f[gr][v], shape=shp, chunks=(np.min((shp[0],f[gr][v].chunks[0])), ))
                            wix = f[gr][v][:]
                            vs = [ref]
                            for k in keeplist:
                                vs.append(vs[-1]+(fslices[i].stop-fslices[i].start))
                                i += 1
                            ref = vs[-1]
        
                            g[v][:] = np.array(vs)
                            #else:
                                #g.create_dataset_like(v, f[gr][v], shape=(len(keeplist), ), data=f[gr][v][:][keeplist],
                                                      #chunks=(np.min((len(keeplist),f[gr][v].chunks[0])), ))
                            
                            for a,av in f[gr][v].attrs.items():
                                if a not in ('CLASS','NAME','REFERENCE_LIST','DIMENSION_LIST'):
                                    #print(a,av)
                                    h[gr][v].attrs[a]=av
                        dim_attach(h, gr)
                        print(h.filename.split('/')[0],gr, 'written')
                    glist = ['header_table']
                    if 'station_configuration' in f.keys():
                        glist.append('station_configuration')
                        if f['station_configuration']['latitude'].shape[0] != f['header_table']['latitude'].shape[0]:
                            print(fn, 'WARNING: station_configuration and header_table are not aligned')
                    else:
                        print('ERROR', fn, 'no station configuration')
                    for gr in glist:
                        g = h.create_group(gr)
                        nindex = len(keeplist)
                        for v in f[gr].keys():
                            if v == 'index':
                                g.create_dataset('index', data=np.empty(nindex, dtype='S1'), compression= 32015, compression_opts=(3, ))                       
                            elif 'string' in v:
                                xdata = f[gr][v][:]
                            #else:
                                #xdata = f[gr][j][headerslice][:]
                                g.create_dataset(v, data=xdata)
                            else:
                                shp = list(f[gr][v].shape)
                                shp[0] = nindex
                                chk = f[gr][v].chunks
                                if chk is not None:
                                    chk = list(chk)
                                    chk[0] = np.min((chk[0], nindex))
                                else:
                                    chk = [0]
                                chk = tuple(chk)
                                if chk[0] > 0:
                                    
                                    g.create_dataset_like(v, f[gr][v], shape=shp, chunks=chk, compression= 32015, compression_opts=(3, ))
                                else:
                                    g.create_dataset_like(v, f[gr][v], shape=shp, chunks=None, compression= 32015, compression_opts=(3, ))
                                if len(keeplist) > 0:
                                    
                                    gm = np.empty_like(f[gr][v][:len(keeplist) + 1])
                                    g[v][:] = f[gr][v][:][keeplist]
                        for a,av in f[gr][v].attrs.items():
                            if a not in ('CLASS','NAME','REFERENCE_LIST','DIMENSION_LIST'):
                                #print(a,av)
                                h[gr][v].attrs[a]=av
                            
                
                        dim_attach(h, gr)
                        print(h.filename.split('/')[0], gr, 'written')
                        
                    for k in f.keys():
                        if k not in ['era5fb', 'observations_table', 'recordindices'] + glist:
                            
                            h.create_group(k)
                            for v in f[k].keys():
                                f.copy(f[k][v],h[k],name=v,without_attrs=True)
                            for a in f[k].attrs.keys():
                                h[k].attrs[a]=f[k].attrs[a]
                            dim_attach(h, k)
                    
                    # check
                    fts = f['recordindices']['recordtimestamp'][:]
                    hts = h['recordindices']['recordtimestamp'][:]
                    d = np.searchsorted(fts, hts)
                    
                    for rv in f['recordindices'].keys():
                        if rv[0] in ['0', '1', '3']:
                            fri = f['recordindices'][rv][:]
                            hri = h['recordindices'][rv][:]
                            foov = f['observations_table']['observation_value'][:]
                            hoov = h['observations_table']['observation_value'][:]
                            fhr = f['header_table']['record_timestamp'][:]
                            hhr = h['header_table']['record_timestamp'][:]
                            for i in range(len(d)):
                                if np.any(foov[fri[d[i]]:fri[d[i]] + (hri[i+1]-hri[i])] !=hoov[hri[i]:hri[i+1]]) or \
                                   np.any(fhr[d[i]:d[i] + 1] !=
                                          hhr[i:i+1]):
                                    if ~np.any(np.isnan(foov[fri[d[i]]:fri[d[i]] + (hri[i+1]-hri[i])])):
                                        
                                        print(foov[fri[d[i]]:fri[d[i]] + (hri[i+1]-hri[i])],
                                              hoov[hri[i]:hri[i+1]])
                                        print(fhr[d[i]:d[i] + 1],
                                              hhr[i:i+1], )
                                        print('')
                                
            print(fn, rt.shape[0], len(keeplist), 'some unique values present, written to reduced/')
            return rt.shape[0], len(keeplist)
        except FileNotFoundError:
            print(fn, 'could not be written in', iters, 'attempts')
            time.sleep(0.1)
    raise FileNotFoundError
                

ray_writereduced = ray.remote(writereduced)

def score(k,p, match):
    parlist = '106','107','126'
    sc = [0, 0]
    for par in parlist:
        if k['pmin'+par][p[0]] < k['pmin'+par][p[1]]:  # higher levels is better
            sc[0] += 3
        elif k['pmin'+par][p[0]] > k['pmin'+par][p[1]]:
            sc[1] += 3
            
        if k['rlen'+par][p[0]] > k['rlen'+par][p[1]]: #more levels is better
            sc[0] += 1
        elif k['rlen'+par][p[0]] < k['rlen'+par][p[1]]:
            sc[1] += 1
        
        if k['andep'+par][p[0]] < k['andep'+par][p[1]]: #smaller andep is better
            sc[0] += 0.5
        elif k['andep'+par][p[0]] > k['andep'+par][p[1]]:
            sc[1] += 0.5
        
        if sc[0] == sc[1]:
            for ip in 0, 1:              
                if 'orphan' in  k['id'][p[ip]]:
                    sc[ip] -= 0.1
                    break
                
    #print(sc, match)
    #if sc[0] + sc[1] == 0:
        #x = 0
    #for par in '106', '126':       
        #for x in zip(k['zz'+par][p[0]], k["ws"+par][p[0]], k['ws'+par][p[1]]):
            #print(f'\t{x}')
        
    #if k['dup'][p[0]]:
        #sc[0] = -1
    #elif k['dup'][p[1]]:
        #sc[1] = -1
    
    return 0 if sc[0] < sc[1] else 1 

#else:
    #print(fn, 'full duplicate to '+str(fnids)+', can be removed')
    #try:
        #os.mkdir('remove')
    #except:
        #pass
    #shutil.copyfile(fn, 'remove/'+fn)
    
@njit
def makepairs(idx, ldist):

    pairs =np.zeros((ldist*(ldist+1)//2, 2), dtype=np.int32)
    k = 0
    for l in range(idx.shape[0]):
            
        isum = 0
        m = ldist - 1
        while (isum + ldist - 1 < idx[l]):
            isum += m
            m -= 1
        if ldist-1-m != idx[l]-isum:          
            pairs[k, 0] = ldist-1-m
            pairs[k, 1] = idx[l]-isum
            k += 1
            
    return pairs[:k]        

def finddups(k):
    
    ref = datetime.datetime(1900, 1, 1)
    parlist = '106','107','126'
    if len(k['pmin106'])>1 or len(k['pmin126'])>1:
        ttt = time.time()
        ldist = len(k['lon'])
        dists = np.zeros((ldist * (ldist + 1))//2)
        #ldist = ldist * (ldist - 1)
        #dists[0] = 0.
        tdist(dists, np.array(k['lat']), np.array(k['lon']), 0)
        dists = np.int32(dists*6370.)
        idx = np.where(dists<200)
        if len(idx[0]) > 2:
            print(ref+datetime.timedelta(seconds=int(k['rt'][0])))
            #print(k['sec']/86400/365.25)#, k['sec']-np.array(k['rt']))
            #print(k, dists[idx])
            l = 0
            ttt = time.time()
            pairs = makepairs(idx[0], ldist)
                #if pairs[-1][0] != 0:
                    #x = 0
            #print(time.time()-ttt)   
            dups = []
            threshs = {'106': 0.5,'107': 0.05,'126': 0.1}
            l = 0
            for p in pairs:
                
                if k['dup'][p[0]] or k['dup'][p[1]]:
                    continue
                sim = False
                for par in '107', '126':
                    if(sim):
                        #print(par, 'check no more needed')
                        continue
                    p0 = np.array(k['ws'+par][p[0]])
                    p1 = np.array(k['ws'+par][p[1]])
                    if par == '107':
                        pd0 = np.array(k['ws106'][p[0]])
                        pd1 = np.array(k['ws106'][p[1]])
                        
                    z0 = np.array(k['zz'+par][p[0]])
                    z1 = np.array(k['zz'+par][p[1]])
                    match = 0
                    if z0.ndim > 0 and z1.ndim > 0:
                        
                        idx = np.searchsorted(z0, z1)
                        for i in range(idx.shape[0]):
                            if idx[i] < p0.shape[0]:
                                if np.abs(p0[idx[i]]-p1[i]) < threshs[par]:
                                    if par == '126':                                       
                                        match += 1
                                    elif pd0.ndim > 0 and pd1.ndim > 0 and i <pd1.shape[0] and idx[i] < pd0.shape[0] and np.abs(pd0[idx[i]]-pd1[i]) < threshs['106']:
                                        match += 1                                        
                                else:
                                    if par == '107' and p0[idx[i]] !=0. and p1[i] !=0. :
                                        good = pd0.ndim > 0 and pd1.ndim > 0 and i <pd1.shape[0] and idx[i] < pd0.shape[0]
                                        if good and np.abs(p0[idx[i]]*1.852/3.6/p1[i]-1.0) < 0.01 and np.abs(pd0[idx[i]]-pd1[i]) < threshs['106']:
                                            match += 1
                                            if match > 3:
                                                print('KNOTS', match, z1[i], p0[idx[i]]*1.852/3.6, p1[i])
                                        elif good and np.abs(p0[idx[i]]/(p1[i]*1.852/3.6) -1) < 0.01 and np.abs(pd0[idx[i]]-pd1[i]) < threshs['106']:
                                            match += 1
                                            if match > 3:
                                                
                                                print('KNOTS', match, z1[i], p0[idx[i]], p1[i]*1.852/3.6)
                                            
                                
                    if match >3:
                        
                    #minlen = np.min((len(k['ws'+par][p[0]]),len(k['ws'+par][p[1]])))
                    #if np.sum(np.isin(k['zz'+par][p[0]],k['zz'+par][p[1]])) > 3:
                        ##if np.sum(np.isin(k['ws'][0],k['ws'][1])) > np.sum(~np.isin(k['ws'][0],k['ws'][1])):
                        #if np.sum(np.abs(k['ws'+par][p[0]][:minlen]-k['ws'+par][p[1]][:minlen])<0.5) > 3:# > np.sum(~np.isin(k['ws'][0],k['ws'][1])):
                            ##for ip in p:
                                ##print(k['id'][ip], dists[idx[0][l]], ref+datetime.timedelta(seconds=int(k['rt'][ip])),
                                      ##k['pmin'][ip], k['rlen'][ip], np.nanstd(k['andep'][ip]), score(k, p)) #k['zz'][ip], k['ws'][ip])
                            ##print('likely duplicates', score(k, p))
                            sim = True
                        #if not sim and par == '107':
                            #if np.sum(np.abs(k['ws'+par][p[0]][:minlen]*1.852/3.6-k['ws'+par][p[1]][:minlen])<0.5) > 3:
                                #sim = True
                                #print('KNOTS')
                            #if np.sum(np.abs(k['ws'+par][p[0]][:minlen]-k['ws'+par][p[1]][:minlen]*1.852/3.6)<0.5) > 3:
                                #sim = True 
                                #print('KNOTS')

                if sim:
                    k['dup'][p[score(k, p, match)]] = True
                    if p[0] == 0 and  k['dup'][0]:
                        x = 0
                    dups.append(p)
                else:
                    pass
                    #for par in '106','107', '126':       
                        #p0 = np.array(k['ws'+par][p[0]])
                        #p1 = np.array(k['ws'+par][p[1]])
                        #z0 = np.array(k['zz'+par][p[0]])
                        #z1 = np.array(k['zz'+par][p[1]])
                        #if z0.ndim > 0 and z1.ndim > 0:
                            #idx = np.searchsorted(z0, z1)
                            #match = 0
                            #for i in range(idx.shape[0]):
                                #if idx[i] < p0.shape[0]:
                                    #print('no',k['id'][p[0]][:-16], k['id'][p[1]][:-16], z1[i], p0[idx[i]], p1[i])
                    #for fn in fns:
                        #for ip in p:
                            #if k['id'][ip] in fn:
                                
                                #with h5py.File(fn) as f:
                                    #frt = f['recordindices']['recordtimestamp'][:]
                                    #fri = f['recordindices']['106'][:]
                                    #idt = np.where(frt==k['rt'][0])[0]
                                    #print(k['id'][ip][:-16], f['observations_table']['z_coordinate'][:][fri[idt[0]]:fri[idt[0]+1]],
                                          #f['observations_table']['observation_value'][:][fri[idt[0]]:fri[idt[0]+1]])
                    x = 0
                    
                l += 1

            if dups:
                #for d in dups:
                    #print(ref+datetime.timedelta(seconds=int(k['rt'][0])),
                          #k['id'][d[0]][:-16], k['dup'][d[0]],k['id'][d[1]][:-16], k['dup'][d[1]])
                print(len(dups), 'duplicates', end=',')
            
            print(time.time()-ttt)
                    
    return k['dup']

ray_finddups = ray.remote(finddups)

def do_year(year):
    
    parlist = '106','107','126'
    dtemplate =  dict(id=[], lat=[], lon=[], ulatlon=[], rt=[])
    dtemplate0 =  dict(id=[], lat=[], lon=[], ulatlon=[], rt=[])
    for par in parlist:
        for prop in 'ws', 'zz', 'andep', 'rlen', 'pmin':         
            dtemplate[prop+par]=[np.nan]
            dtemplate0[prop+par]=[]

    #fns = glob.glob(f'/mnt/users/scratch/leo/scratch/converted_v13/{year}/*v1.nc')[:] 
    fns = glob.glob(f'/mnt/users/scratch/leo/scratch/converted_v24/{year}/*v3.nc')[:]
    if len(fns) == 0:
        print(f'No files found: /mnt/users/scratch/leo/scratch/converted_v24/{year}/*v3.nc')
        return
    
    daydicts = []
    try:
            with open(f'daydictrs{fns[0].split("/")[-2]}.pkl', 'rb') as f:
                daydict = pickle.load(f)
            print(time.time()-tt)
    except:
        
        futures = []
        for fn in fns[:]:
            if True or '0-20000-0-02464' in fn:
                daydicts.append(fextractrs(fn))
                #futures.append(ray_fextractrs.remote(fn))
            
        #daydicts = ray.get(futures)
            
        daydict = {}
        #for k in range(46000*12):
            #daydict[k] = dict(sec=k*7200, id=[], lat=[], lon=[], ulatlon=[], rlen=[], pmin=[], ri=[], rt=[], ws=[], zz=[], andep=[], dup=[])
        for do in daydicts:
            #d = ray.get(do)
            d = do
            for k in d.keys():
                if k not in daydict.keys():
                    daydict[k] = copy.deepcopy(dtemplate0)
                    daydict[k]['dup'] = []
                    
                if k != 'fn':
                    
                    for p in d[k].keys(): #['id', 'lat', 'lon', 'ulatlon', 'rlen', 'pmin', 'ri', 'rt', 'ws', 'zz', 'andep']:
                        if d[k][p]:
                            daydict[k][p].append(d[k][p][0])
                    daydict[k]['dup'].append(False)
                    #daydict[k]['ws'].append(d[k]['ws'][0][:])
                    #if daydict[k]['id']:
                        
                        #print(k, daydict[k]['id'])
                else:
                    print(d['fn'])
        
        
        pot_dups = [len(daydict[k]['pmin106']) for k in daydict.keys()]
        
        with open(f'daydictrs{fns[0].split("/")[-2]}.pkl', 'wb') as f:
            pickle.dump(daydict, f)
    
    print('after read', f'{time.time()-tt:.3f}')
    ref = datetime.datetime(1900, 1, 1)
    #daydictref = ray.put(daydict)
        
    futures = []
    for k, v in daydict.items():
        #print(k['sec']/86400/365.25)#, k['sec']-np.array(k['rt']))
        if k != 'fn':
            dtk = ref+datetime.timedelta(seconds=int(k)*7200)
            if dtk > datetime.datetime(1980, 6, 1, 22) and  dtk < datetime.datetime(1980, 6, 2, 1):
                x = 0
            daydict[k]['dup'] = finddups(v)
        #futures.append(ray_finddups.remote(v))
        
    #dups = ray.get(futures)
    #l = 0
    #for k in daydict.keys():
        #daydict[k]['dup'] = dups[l]
        #l += 1
    print('after duplicate check', f'{time.time()-tt:.3f}')
    
    ndaydict ={}
    for k, v in daydict.items():
        ndaydict[k] = {'id': v['id'], 'dup': v['dup']}
        
    #obj =ray.put(ndaydict)                    
    print('after put', f'{time.time()-tt:.3f}')
    
    x = 0
    futures = []
    keepstats = []
    if True:
        for fn in fns:
            #if '0-20000-0-02464' in fn:
                keepstats.append(writereduced(fn, ndaydict, ref))
    else:
        for fn in fns:
            futures.append(ray_writereduced.remote(fn, obj, ref))
        
        keepstats = ray.get(futures)
    
    keepsum =0
    rtsum = 0
    full = 0
    for k in keepstats:
        rtsum += k[0]
        keepsum += k[1]
        if k[1] == 0:
            full += 1
    years.append(year)
    keepsums.append(keepsum)
    rtsums.append(rtsum)
    fulls.append(full)
    print('total records:', rtsum, 'kept records:', keepsum, 'full duplicates:', full)
    print('after write', time.time()-tt)
        
ray_do_year = ray.remote(do_year)

################################################################
        
if __name__ == '__main__':

    #fns = glob.glob('/mnt/scratch/scratch/federico/MERGED_FEB2023_MOBILE_mobile/*v1.nc')[:]
    fns = glob.glob('/mnt/users/scratch/leo/scratch/converted_v13/1980/0-20999-0-UH*v1.nc')[:] +\
        glob.glob('/mnt/users/scratch/leo/scratch/converted_v13/1980/20999*UH*v1.nc')[:]
    
    years=[]
    keepsums=[]
    rtsums=[]
    fulls=[]
    
    ray.init()

    
    tt = time.time()
    futures = []
    for year in range(2023, 1904, -1):
        
        if False:
            futures.append(ray_do_year.remote(year))
        else:
            do_year(year)
        
    ray.get(futures)
    
    for year in range(2005, 2004, -1):
        
        path = f'/mnt/users/scratch/leo/scratch/converted_v19/{year}/'
        os.makedirs(path+'original', exist_ok=True)
        fns = glob.glob(path+'reduced/*v3.nc')
        for fn in fns:
            try:
                if not os.path.exists(path+'original/'+os.path.basename(fn)):
                    
                    shutil.move(path+os.path.basename(fn),path+'original/'+os.path.basename(fn))
                shutil.copyfile(fn,path+os.path.basename(fn) )
            except:
                print(fn, 'copy failed')
        
    print('total records:', np.sum(rtsums), 'kept records:', np.sum(keepsums), 'full duplicates:', np.sum(fulls))
    plt.plot(years, np.array([rtsums, keepsums]))
    plt.savefig('dups.png')
    plt.show()
    
    ref = datetime.date(1900, 1, 1)
    ddlen = 0
    dlist = []
    ulist = []
    klist = []
    for k in range(len(daydict)):
        dlen = len(daydict[k]['id'])
        if dlen > 0:
            #print(ref+datetime.timedelta(days=k),dlen)
            ddlen += dlen
            klist.append(k)
            dlist.append(dlen)
            ulist.append(len(np.unique(daydict[k]['ulatlon'])))
            
    print('print',ddlen, time.time()-tt)
    plt.plot(np.array(klist)/325.25+1900, np.array(dlist))
    plt.plot(np.array(klist)/325.25+1900, np.array(ulist))
    plt.show()

