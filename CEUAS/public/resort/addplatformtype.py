import os,glob,shutil
import h5py
import numpy as np
import datetime
import time
from global_land_mask import globe
import ray

#tt = time.time()
#for lat in range(-90, 90):
    #for lon in range(-180, 180):
        
        #print(globe.is_land(lat, lon))
#print(time.time()-tt)

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

scplatform_type = {0: [16013, 16022, 16045, 16068], # land station
                   2:[16014, 16019, 16046, 16069], #ship 
                   8: [16021, 16075], # land vehicle
                   }
                   #41: [16086], # dropsonde (aircraft platform) should not be used
                   #10: [16087] # BUFR descent - so far no platform_type for this should not be used
os.chdir(os.path.expandvars('$RSCRATCH/converted_v13/long'))
fns=glob.glob('20999*.nc') + glob.glob('0-20999-0*.nc')

i=0
def addpf(fn, scplatform_type):
    
    #if '0-20' == fn[:4] or '20999' in fn[:5]:
        #continue
    droplist = []
    keeplist = []
    mismatch = 0
    with h5py.File(fn,'r+') as f:
        try:
            par = '107'
            wix = f['recordindices/'+par][:]
            wsix = slice(wix[0],wix[-1] )
        except:
            par = list(f['recordindices'].keys())[0]
            wix = f['recordindices/'+par][:]
            wsix = slice(wix[0],wix[-1] )
        
        pt = np.unique(f['era5fb']['reportype'][wsix])
        st = []
        for p in pt:
            for k, v in scplatform_type.items():
                
                if p > 0 and p in v:
#                    if k not in st:
                    st.append(k)
        #d,idx = np.unique(f['header_table/latitude'][:], return_index=True)
        #globe.is_land(f['header_table/latitude'][:], f['header_table/longitude'][:][idx])
        #if np.sum(globe.is_land(d, f['header_table/longitude'][:][idx]))/d.shape[0] < 0.5:
            #print(fn, 'ship')
        #else:
            #print(fn, 'land')
        if len(st) > 1:
            #print('ERROR 2 subtypes in 1 record', fn, pt)
            st = [st[-1]]
        elif len(st) ==0:
            try:
                
                hs = f['header_table']['primary_station_id'][0].view(f'S{f["header_table"]["primary_station_id"].shape[1]}')[0].split(b'/')[-1]
            except Exception as e:
                print(fn, e)
                hs = b''
                if len(pt) == 0:
                    pt = np.array([-1])
            if b'20999' == hs[:5]:
                if b'ncar' in hs:
                    hs = hs.split(b'-')[-1].split(b'.')[0].split(b'_')[-1].decode()
                elif b'uadb_' in hs:
                    print('x')
                elif b'npsound_' in hs:
                    print(fn.split('/')[-1],'unknown, assuming ship')
                    hs = 'x'
                elif b'bufr_' in hs:
                    print('x')
                    hs = hs.split(b'.')[-3].decode()
                elif b'conv' in hs:
                    print('x')
                    if b'txt' in hs:
                        hs = hs.split(b'.')[-4].decode()
                    else:
                        hs = hs.split(b'.')[-3].decode()
                        
                    if hs[0] == '_':
                        hs = hs[1:]
                    if ':' in hs:
                        hs = hs.split(':')[-1]
                    if hs == '':
                        print(fn.split('/')[-1],'unknown, assuming land')
                        hs = '0'
                    elif hs.isnumeric() and not (np.nanstd(f['header_table']['longitude'][:]) > 2.0 or np.nanstd(f['header_table']['latitude'][:]) > 2.0):
                        print(fn.split('/')[-1],'unknown, assuming land')
                        hs = '0'
                    else:
                        print(fn.split('/')[-1],'unknown, assuming ship', f['header_table']['latitude'][:], f['header_table']['longitude'][:])
                        hs = 'S'
                        
                            
                elif b'igra' in hs:
                    hs = hs.split(b'_')[-1].split(b'-')[0].decode()
                    if 'ZZ' in hs:
                        hs = 'Z'
                    else:
                        hs = '0'
                elif b'amma' in hs:
                    hs = hs.split(b'_')[-2].decode()
                elif b'shipsound' in hs:
                    hs = 'S'
                else:
                    hs = hs.split(b'-')[-1].split(b'.')[-4].decode()
                    if hs[0] == '_':
                        hs = hs[1:]
                    if ':' in hs:
                        hs = hs.split(':')[-1]
                    
                print(hs)
                if hs[0].isnumeric():
                    print('no subtype, substituting land')
                    st.append(0)
                else:
                    print('no subtype, substituting ship')
                    st.append(2)
                    
            else:
                hs = hs.split(b'-')[-1].split(b'_')[0].decode()
                
                if hs !='':
                    
                    if hs[0] == '_':
                        hs = hs[1:]
                    if ':' in hs:
                        hs = hs.split(':')[-1]
                print(hs)
                if hs.isnumeric():
                    if ('0-20'!=fn[4:] or '20999' in fn) and (f['header_table']['latitude'].shape[0] < 20 or np.nanstd(f['header_table']['longitude'][:]) > 1.0 or np.nanstd(f['header_table']['latitude'][:]) > 1.0):
                        print('no subtype, substituting ship')
                        st.append(2)
                    else:    
                        print('no subtype, weather ship?, substituting land ', f['header_table']['latitude'][-1], f['header_table']['longitude'][-1])
                        st.append(0)
                else:
                    if '0-20'!=fn[4:] and '20999' not in fn:
                        if (f['header_table']['latitude'].shape[0] > 20 and np.nanstd(f['header_table']['longitude'][:]) < 1.0 and np.nanstd(f['header_table']['latitude'][:]) < 1.0):
                            print('no subtype, new WIGOS, substituting land')
                            st.append(0)
                        else:
                            print('no subtype, substituting ship')
                            st.append(2)
                    else:
                        print('no subtype, substituting ship')
                        st.append(2)
        else:
            pass
        dlon = f['header_table/longitude'][:]
        dlat = f['header_table/latitude'][:]
        dlat[dlat<-90.] = -90.
        dlat[dlat>90.] = 90.
        idx = ~np.isnan(dlon) & ~np.isnan(dlat)
        if np.sum(idx) > 0:
            
            isl = np.sum(globe.is_land(dlat[idx], dlon[idx]))/dlat.shape[0]
            if len(pt) > 0 and pt[-1] < 0:
                if st[0] == 0 and isl < 0.1:
                    
                    print(fn, 'st', st, pt, isl, )
                    st[0] = 2
                    print('x')
                elif st[0] == 2 and isl > 0.5:
                    print(fn, 'st', st, pt, isl, )
                    st[0] = 0
                    print('x')
        try:
            
            f['header_table']['platform_type'][:] =st[0] 
            f['station_configuration']['platform_type'][:] =st[0]
            if f['station_configuration']['platform_type'][-1] != st[0]:
                print(fn.split('/')[-1], 'fixed platform_type', f['station_configuration']['platform_type'][-1], st[0])
        except Exception as e:
            try:
                
                f.create_group('station_configuration')
            except:
                pass
            try:
                
                f['station_configuration'].create_dataset_like('platform_type', f['header_table']['platform_type'])
                f['station_configuration'].create_dataset('index',
                                                          data=np.empty(f['header_table']['platform_type'].shape[0], dtype='S1'),
                                                          compression= 'gzip')                       

            except Exception as e:
                print(fn, e)
            dim_attach(f, 'station_configuration')   
            print(fn, e, ',fixed')


ray_addpf = ray.remote(addpf)

ray.init(num_cpus=40)

tt = time.time()
futures = []
for fn in fns:
    futures.append(ray_addpf.remote(fn, scplatform_type))
    #reduce_orphan(fn, glats, glons, gids)
obj_ref  = ray.get(futures)
               
print(time.time()-tt)
            
exit(0) 
for fn in fns[:]:
    droplist = []
    keeplist = []
    mismatch = 0
    with h5py.File(fn,'r') as f:
        maxdate=ithresh +1 #f['recordindices/recordtimestamp'][-1]
        
        if maxdate>ithresh:
            
            ix=np.searchsorted(f['recordindices/recordtimestamp'][:],ithresh)
            rix = f['recordindices/recordtimestamp'][:]#[ix:]
            try:
                par = '107'
                wix = f['recordindices/'+par][:]
                wsix = slice(wix[0],wix[-1] )#[ix:]
            except:
                par = list(f['recordindices'].keys())[0]
                wix = f['recordindices/'+par][:]
                wsix = slice(wix[0],wix[-1] )#[ix:]
                #wsix = f['recordindices/'+par][:-1]#[ix:]
            try:
                lat = f['observations_table']['latitude'][wsix]
                lon = f['observations_table']['longitude'][wsix]
            except:
                print('x')
            #lon=np.unique(f['observations_table']['longitude'][ix:])
            #lat=np.unique(f['observations_table']['latitude'][ix:])
            #print(np.max(lon)-np.min(lon))
            #print(np.max(lat)-np.min(lat))
            #if np.max(lat)-np.min(lat)>0.3 or np.max(lon)-np.min(lon)>0.3:
                #i+=1
            print(fn)
            if '0-20300' in fn:
                key = fn.split('_')[0].split('-')[-1]
            else:
                key = fn.split('_')[1]
                if 'era5' in key and 'conv' in key:
                    if '?' not in key:
                        try:
                            
                            key = fn.split('_')[2].split('.gz')[0].split(':')[1]
                        except:
                            key = ''
                    else:
                        key = key.split('.')[-3]
                if 'era5' in key and 'bfr' in key:
                    key = key.split('.')[1]
                    while len(key) < 5:
                        key = '0' + key
                if 'ZZ' in key:
                    key = key.split('-')[0][-5:]
                    while key[0] == '0':
                        key = key[1:]
                if 'uadb' in key:
                    key = fn.split('_')[3].split('.')[0]
                    
            
            print(key)
            fnids = glob.glob('0-20???-0-'+key+'_CEUAS_merged_v1.nc')
            threshs = []
            if len(fnids) > 0:
                threshs = [1.0] * len(fnids)
                if len(fnids) > 1:
                    
                    print(fnids)
            if 'ZZ' not in fn:
                ulats,idx = np.unique(lat, return_index=True)
                ulons = lon[idx]
                for i in range(len(ulats)):
                    udx = np.searchsorted(glats, ulats[i])
                    j = 0
                    while udx < glats.shape[0] and udx - j >= 0 and np.abs(glats[udx-j] -ulats[i]) < 0.3:
                        if np.abs(glons[udx-j] - ulons[i]) < 0.3 and gids[udx-j] not in fnids:
                            fnids.append(gids[udx-j])
                            threshs.append(0.3)
                        j += 1
                    k = 1
                    while udx + k <glats.shape[0] and np.abs(glats[udx+k] -ulats[i]) < 0.3 :
                        k += 1
                        if np.abs(glons[udx+k] - ulons[i]) < 0.3 and gids[udx+k] not in fnids:
                            fnids.append(gids[udx+k])
                            threshs.append(0.3)
                if len(fnids) > 0:
                    print(fnids)
                if len(fnids) > 1:
                    
                    print(fnids)
                    
            # now check candidates
            for fnid,thresh in zip(fnids, threshs):
                with h5py.File(fnid,'r') as g:
                    #iy=np.searchsorted(g['recordindices/recordtimestamp'][:],ithresh)
                    gk = len(list(g['recordindices'].keys()))
                             
                    riy = g['recordindices/recordtimestamp'][:]#[iy:]
                    try:
                        
                        wiy = g['recordindices/'+par][:]
                        wsiy = slice(wiy[0],wiy[-1]) #[ix:]
                    except:
                        continue
                    
                    #print(g['recordindices/recordtimestamp'][iy]-f['recordindices/recordtimestamp'][ix])
                    lony = g['observations_table/longitude'][wsiy] #[iy:])
                    laty = g['observations_table/latitude'][wsiy] #[iy:])
                    match = np.searchsorted(riy, rix)
                    match = match[match<riy.shape[0]]
                    foo = f['observations_table']['observation_value'][wsix]
                    goo = g['observations_table']['observation_value'][wsiy]
                    for m in range(len(match)):
                        if (match[m] < riy.shape[0]  and (riy[match[m]] - rix[m] <3700 ) and m < lat.shape[0]) or \
                           (match[m] < riy.shape[0]  and (riy[match[m]] - rix[m-1]<3700) and m < lat.shape[0]):
                            #if match[m] < riy.shape[0]  and (np.abs(riy[match[m]] - rix[m-1])<3700) and m < lat.shape[0]:
                                #print('x')
                            try:
                                
                                if np.abs(lony[match[m]]-lon[m]) < thresh and np.abs(laty[match[m]]-lat[m]) < thresh:
                                    if m == match.shape[0]-1:
                                        if len(foo) - wix[m] <= len(goo) - wiy[match[m]]:
                                            droplist.append(m)
                                    else:
                                        if wix[m+1]-wix[m] <= wiy[match[m+1]] - wiy[match[m]]:
                                        
                                            droplist.append(m)
                                        else:
                                            #print(f['observations_table']['z_coordinate'][wsix][wix[m]-wix[0]:wix[m+1]-wix[0]],
                                                  #g['observations_table']['z_coordinate'][wsiy][wiy[match[m]]- wiy[0]:wiy[match[m+1]]- wiy[0]])
                                            print('orphan record larger than original record', wix[m+1]-wix[m], wiy[match[m+1]] - wiy[match[m]])

                                    #if m < match.shape[0]-1 and m > match.shape[0]-4:
                                        
                                        #print(foo[wix[m]-wix[0]:wix[m+1]-wix[0]], goo[wiy[match[m]]-wiy[0]:wiy[match[m+1]]-wiy[0]])
                                        ##print(f['observations_table']['date_time'][wsix][wix[m]-wix[0]:wix[m+1]-wix[0]], g['observations_table']['date_time'][wsiy][wiy[match[m]] - wiy[0]:wiy[match[m+1]]- wiy[0]])
                                        #print(f['observations_table']['z_coordinate'][wsix][wix[m]-wix[0]:wix[m+1]-wix[0]],g['observations_table']['z_coordinate'][wsiy][wiy[match[m]]- wiy[0]:wiy[match[m+1]]- wiy[0]])
                                        #print('')
                                    if m == match.shape[0]-1:
                                        print(foo[wix[m]-wix[0]:], goo[wiy[match[m]]-wiy[0]:])
                                        #print(f['observations_table']['date_time'][wsix][wix[m]-wix[0]:], g['observations_table']['date_time'][wsiy][wiy[match[m]]-wiy[0]:])
                                        print(f['observations_table']['z_coordinate'][wsix][wix[m]-wix[0]:], g['observations_table']['z_coordinate'][wsiy][wiy[match[m]]-wiy[0]:])
                                        print('')
                                    
                                    if foo[m] != goo[match[m]]:
                                        if m < match.shape[0] - 1 and foo[wix[m+1]-wix[0] -1] != goo[wiy[match[m+1]]-wiy[0]-1]:
                                            mismatch += 1
                                        elif  m == match.shape[0] - 1:
                                            mismatch += 1
                                        
                                else:
                                    keeplist.append(m)
                            except Exception as e:
                                print(fn, fnid, e)
                                droplist.append(m)
                        else:
                            keeplist.append(m)
                    #print(lon, lony)
                    #print(lat, laty)
                    print(len(droplist), len(keeplist), mismatch)
            if len(droplist) == 0:
                shutil.copyfile(fn, 'reduced/'+fn)
            elif len(droplist) <rix.shape[0]:
                keeplist = []
                nindex = 0
                j = 0
                for i in range(rix.shape[0]):
                    if i == droplist[j]:
                        j += 1
                        if j == len(droplist):
                            break
                    else:
                        keeplist.append(i)
                fslices = []
                for v in f['recordindices'].keys():
                    if v == '0' or '1' in v:
                        wix = f['recordindices'][v][:]
                        for k in keeplist:
                            fslices.append(slice(wix[k], wix[k+1]))
                            nindex+=wix[k+1]-wix[k]
                        
                with h5py.File('reduced/'+fn, 'w') as h:
                    
                    for gr in ['era5fb', 'observations_table']:
                        g = h.create_group(gr)
                        for v in f[gr].keys():
                            if v == 'index':
                                g.create_dataset('index', data=np.empty(nindex, dtype='S1'), compression= 'gzip')                       
                            elif 'string' in v:
                                xdata = f[gr][v][:]
                                g.create_dataset(v, data=xdata)
                            else:
                                shp = list(f[gr][v].shape)
                                shp[0] = nindex
                                chk = list(f[gr][v].chunks)
                                chk[0] = np.min((f[gr][v].chunks[0], nindex))
                                chk = tuple(chk)
                                if chk[0] == 0:
                                    g.create_dataset_like(v, f[gr][v], shape=shp, chunks=None)
                                else:
                                    g.create_dataset_like(v, f[gr][v], shape=shp, chunks=chk)
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
                        print(gr, 'written')
                        
                    for gr in ['recordindices']:
                        g = h.create_group(gr)
                        for v in f[gr].keys():
                            if v == '0' or '1' in v:
                                wix = f['recordindices'][v][:]
                                for k in keeplist:
                                    fslices.append(slice(wix[k], wix[k+1]))
                                    nindex+=wix[k+1]-wix[k]
                        nindex = len(keeplist)
                        ref = 0
                        i = 0
                        for v in f[gr].keys():
                            if v == 'index':
                                g.create_dataset('index', data=np.empty(nindex, dtype='S1'), compression= 'gzip')                       
                            elif v == '0' or '1' in v:
                                shp = list(f[gr][v].shape)
                                shp[0] = nindex + 1
                                g.create_dataset_like(v, f[gr][v], shape=shp)
                                wix = f[gr][v][:]
                                vs = [ref]
                                for k in keeplist:
                                    vs.append(vs[-1]+(fslices[i].stop-fslices[i].start))
                                    i += 1
                                ref = vs[-1]

                                g[v][:] = np.array(vs)
                            else:
                                g.create_dataset_like(v, f[gr][v], shape=(len(keeplist), ), data=f[gr][v][:][keeplist])
                            
                            for a,av in f[gr][v].attrs.items():
                                if a not in ('CLASS','NAME','REFERENCE_LIST','DIMENSION_LIST'):
                                    #print(a,av)
                                    h[gr][v].attrs[a]=av
                        dim_attach(h, gr)
                        print(gr, 'written')
                    for gr in ['header_table']:
                        g = h.create_group(gr)
                        nindex = len(keeplist)
                        for v in f[gr].keys():
                            if v == 'index':
                                g.create_dataset('index', data=np.empty(nindex, dtype='S1'), compression= 'gzip')                       
                            elif 'string' in v:
                                xdata = f[gr][v][:]
                            #else:
                                #xdata = f[gr][j][headerslice][:]
                                g.create_dataset(v, data=xdata)
                            else:
                                shp = list(f[gr][v].shape)
                                shp[0] = nindex
                                chk = list(f[gr][v].chunks)
                                chk[0] = np.min((f[gr][v].chunks[0], nindex))
                                chk = tuple(chk)
                                if chk[0] > 0:
                                    
                                    g.create_dataset_like(v, f[gr][v], shape=shp, chunks=chk)
                                else:
                                    g.create_dataset_like(v, f[gr][v], shape=shp, chunks=None)
                                if len(keeplist) > 0:
                                    
                                    gm = np.empty_like(f[gr][v][:len(keeplist) + 1])
                                    g[v][:] = f[gr][v][:][keeplist, ]
                        for a,av in f[gr][v].attrs.items():
                            if a not in ('CLASS','NAME','REFERENCE_LIST','DIMENSION_LIST'):
                                #print(a,av)
                                h[gr][v].attrs[a]=av
                            
                
                        dim_attach(h, gr)
                        print(gr, 'written')
                        
                    for k in f.keys():
                        if k not in ['era5fb', 'observations_table', 'header_table', 'recordindices']:
                            
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
                        if '0' in rv or '1' in rv:
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
                                
                        
                print('some unique values present')
                
                
            
                   
                
                
                    
                    

    
        
print(i)
exit()
i=0
try:
    os.mkdir('orph')
except:
    pass
for fn in fns:
    i+=1
    wigos=f'0-23001-2-orph{i:04d}'
    fno='orph/'+wigos+'_CEUAS_merged_v1.nc'
    wigos=np.string_(wigos)
    shutil.copyfile(fn,fno)
    f=h5py.File(fno,'r+')
    #print(f['station_configuration/primary_id'][:])
    print(f['station_configuration/primary_id'].shape)
    l=f['station_configuration/index'].shape[0]
    del f['station_configuration/primary_id']
    f['station_configuration'].create_dataset('primary_id',data=(np.vstack([wigos]*l)).view('S1'))
    
    l=f['header_table/index'].shape[0]
    del f['header_table/primary_station_id']
    f['header_table'].create_dataset('primary_station_id',data=(np.vstack([wigos]*l)).view('S1'))
    print(fn,fno)
    
