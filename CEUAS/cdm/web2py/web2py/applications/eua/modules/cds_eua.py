import os,sys
import h5py #ickle as h5py
import numpy
from numba import njit
import xarray
import pandas as pd
import copy
import time

#@njit(cache=True)
def find_dateindex(y,x):
    """ creates the indices list from the dates, for quick access 
        nb the benchmark script will not work with these files since the definition of the array size is swapped i.e. (x.shape[0], 3)"""        


    #x=y#numpy.unique(y)
    z=numpy.zeros((3,x.shape[0]),dtype=numpy.int32)
    z-=1
    j=0
    for i in range(len(y)):
        m=y[i]
        if x[j]==y[i]:
            if z[1,j]==-1:
                z[1,j]=i
                #print(j,i)
            else:
                if z[2,j]<i:
                    z[2,j]=i
        elif x[j]<y[i]:
            j+=1
            if x[j]==y[i]:
                if z[1,j]==-1:
                    z[1,j]=i
                    #print(j,i)
                else:
                    if z[2,j]<i:
                        z[2,j]=i
            else:
                print('Error')
        else:
            j-=1
            if x[j]==y[i]:
                if z[1,j]==-1:
                    z[1,j]=i
                    #print(j,i)
                else:
                    if z[2,j]<i:
                        z[2,j]=i
            else:
                print('Error')
    z[0,:]=x
    return z

@njit
def find_dateindex_cg(y):
    
    x=numpy.unique(y)
    z=numpy.zeros((x.shape[0],3),dtype=numpy.int32)
    z-=1
    j=0
    for i in range(len(y)):
        m=y[i]
        if x[j]==y[i]:
            if z[j,0]==-1:
                z[j,0]=i
                #print(j,i)
            else:
                if z[j,1]<i:
                    z[j,1]=i
        elif x[j]<y[i]:
            j+=1
            if x[j]==y[i]:
                if z[j,0]==-1:
                    z[j,0]=i
                    #print(j,i)
                else:
                    if z[j,1]<i:
                        z[j,1]=i
            else:
                print('Error')
        else:
            j-=1
            if x[j]==y[i]:
                if z[j,0]==-1:
                    z[j,0]=i
                    #print(j,i)
                else:
                    if z[j,1]<i:
                        z[j,1]=i
            else:
                print('Error')
    z[:,2]=x
    return z

def do_copy(fd,f,k,idx,cut_dimension): # cuts vars and copies attributes of observation, feedback and header tables
    for v in f[k].keys():
        if f[k][v].ndim==1:
            if f[k][v].dtype!='S1':
                
                fd[k].create_dataset_like(v,f[k][v],shape=idx.shape,chunks=True)
                hilf=f[k][v][idx[0]:idx[-1]+1]
                fd[k][v][:]=hilf[idx-idx[0]]
            else:
                if v in [cut_dimension]:
                    
                    fd[k].create_dataset_like(v,f[k][v],shape=idx.shape,chunks=True)
                    hilf=f[k][v][idx[0]:idx[-1]+1]                                                
                    fd[k][v][:]=hilf[idx-idx[0]]
                else:
                    fd[k].create_dataset_like(v,f[k][v])
                    fd[k][v][:]=f[k][v][:]
                    pass
        else:
            fd[k].create_dataset_like(v,f[k][v],shape=(idx.shape[0],f[k][v].shape[1]),chunks=True)
            hilf=f[k][v][idx[0]:idx[-1]+1,:]
            fd[k][v][:]=hilf[idx-idx[0],:]
        for a in f[k][v].attrs.keys():
            if a not in ['DIMENSION_LIST','CLASS']:
                fd[k][v].attrs[a]=f[k][v].attrs[a]
    for v in f[k].keys():
        l=0
        for d in f[k][v].dims:
            if len(d)>0:
                fd[k][v].dims[l].attach_scale(fd[k][f[k][v].dims[l][0].name])
            l+=1

def cdmexplainer(rvars):
    #rfile='/fio/srvx7/leo/python/CEUAS/CEUAS/public/harvest/data/tables/'+'chera5.conv._01001.nc'
    rfile=os.path.expandvars('$EUA_ROOT/subdaily/v0.1/source/ERA5_1/obs/0-20000-0-'+'01001'+
                             '/eua_subdaily_v0.1_source_ERA5_1_obs_0-20000-0-'+'01001'+'_t.nc')
    print(rfile,rvars)           
    if len(rvars.keys())>=1:
        rvdict=copy.copy(rvars)
        rvdk=rvdict.keys()
        for k in rvdk:
            try:
                rvdict[k]=eval(rvdict[k])
                print(k,type(rvdict[k]),rvdict[k])
            except:
                print('could not evaluate '+rvdict[k])
        #for k in rvdk:
            #if type(rvdict[k]) is not list:
                #if rvdict[k] not in vdict.keys():
                    #rvdict[k]=[rvdict[k],rvdict[k]]
                #else:
                    #rvdict[k]=[vdict[rvdict[k]],vdict[rvdict[k]]]
        print(rfile)           
        #with xarray.open_dataset(rfile,group=rvdict['group']) as x: 
        with h5py.File(rfile) as f: 
            x=f[rvdict['group']]
            print(rfile,x)
            xk=list(x.keys())
#            if rvdict['group']=='era5fb':
#                xk=xk+list(x.dims)
            del rvdict['group']
            y=dict()
            rvdk=list(rvdk)[0]
            ki=xk.index(rvdk)
            xk.insert(0,xk.pop(ki))
            if rvdict:
                #idx=numpy.where(x[rvdk].values.astype(int)==int(rvdict[rvdk]))[0]
                if rvdk in ['index']:
                    idx=[numpy.int64(rvdict[rvdk])]
                else:
                    z=x[rvdk][:]
                    if z.ndim==1:
                        idx=numpy.where(x[rvdk][:]==numpy.int64(rvdict[rvdk]))[0]
                    else:
                        z=z.view('S{0:d}'.format(x[k].shape[1]))
                        idx=numpy.where(z[:,0].astype('int')==numpy.int64(rvdict[rvdk]))[0]
                for k in x.keys():
                    #print(k,x[k].shape)
                    if x[k].ndim>1:
                        #y[k]=x[k].values.view('S80')[idx[0]]
                        if x[k].shape[0]==x[rvdk].shape[0]:
                            
                            y[k]=x[k][idx[0]].view('S{0:d}'.format(x[k].shape[1]))
                            #print(y[k])
                            y[k]=numpy.asarray(y[k],dtype='unicode')
                        #else:
                            #y[k]=numpy.asarray([x[k]],dtype='unicode')
                    else:
                        if x[k].shape[0]==x[rvdk].shape[0]:
                            y[k]=numpy.asarray([x[k][idx[0]]])
                        else:
                            if type(x[k][0]) is numpy.bytes_:
                                #y[k]=x[k][:].view('S{0:d}'.format(x[k].shape[0]))
                                pass
                            else:
                                y[k]=x[k]
                                print(y[k])
            else:
                for k in xk:
                    y[k]=x[k].values.view('S80')
                    y[k]=numpy.asarray(y[k].reshape((y[k].shape[0])),dtype='unicode')
            #print(k,y[k].shape)
            df = pd.DataFrame(y,columns=xk)
            if df.shape[0]==1:
                df=df.T
            
        return df
            
        
    
''' Main routine for parsing the CDS request'''
def process(rvars):

    error=''
    rfile=''
    vdict={'temperature':2,'uwind':3,'vwind':4}
    rvkeys=rvars.keys()
    #cost=calculate_cost(rvars) # estimate size of output file
    if 'statid' in rvkeys:
        #rfile='applications/testajax/static/1/chera5.conv._'+rvars['statid']+'.nc'
        #rfile='/fio/srvx7/leo/python/CEUAS/CEUAS/public/harvest/data/tables/'+'chera5.conv._'+rvars['statid']+'.nc'
        rfile=os.path.expandvars('$EUA_ROOT/subdaily/v0.1/source/ERA5_1/obs/0-20000-0-'+rvars['statid']+
                                 '/eua_subdaily_v0.1_source_ERA5_1_obs_0-20000-0-'+rvars['statid']+'_t.nc')
        print(rfile)
        if len(rvkeys)>1:
            rvdict=copy.copy(rvars)
            del rvdict['statid']
            rvdk=rvdict.keys()
            for k in rvdk:
                try:
                    rvdict[k]=eval(rvdict[k])
                    print(k,type(rvdict[k]),rvdict[k])
                except:
                    print('could not evaluate '+rvdict[k])
            for k in rvdk:
                if type(rvdict[k]) is not list:
                    if rvdict[k] not in vdict.keys():
                        rvdict[k]=[rvdict[k],rvdict[k]]
                    else:
                        rvdict[k]=[vdict[rvdict[k]],vdict[rvdict[k]]]
                                    
            t=time.time()
            with h5py.File(rfile) as f:
                di=f['dateindex']
                #f['dat']=f['dateindex'][0,:]
                try:
                    sdate=rvdict.pop('date')
                    didx=numpy.where(numpy.logical_and(di[0,:]>=sdate[0],di[0,:]<=sdate[1]))[0]
                    didx=[didx[0],didx[-1]]
                except:
                    didx=[0,di.shape[1]-1]
                
                trange=[di[1,didx[0]],di[2,didx[1]]+1]    
                mask=numpy.ones(trange[1]-trange[0],dtype=numpy.bool)
                criteria={'variable':'era5fb/varno@body','level':'era5fb/vertco_reference_1@body'}
                ck=criteria.keys()
                t=time.time()
                for ck,cv in criteria.items():
                    if ck in rvdk:
                        mask=numpy.logical_and(mask,f[cv][trange[0]:trange[1]]>=rvdict[ck][0])
                        mask=numpy.logical_and(mask,f[cv][trange[0]:trange[1]]<=rvdict[ck][1])
                        
                print('mask:',time.time()-t)
                t=time.time()
                #y=f['era5fb/date@hdr'][mask]
                idx=numpy.where(mask)[0]+trange[0]
                print(f['era5fb/date@hdr'][idx[0]])
                z=f['recordindex'][:]
                zidx=numpy.where(numpy.logical_and(z>=trange[0],z<trange[1]))[0]
                
                dfile='applications/eua/static/tmp/dest.nc'
                with h5py.File(dfile,'w') as fd:
                    for k in f.keys():
                        if isinstance(f[k],h5py.Group):
                            fd.create_group(k)
                            if k in ['observations_table','era5fb']: # only obs, feedback fitting criteria (idx) is copied
                                do_copy(fd,f,k,idx,'obslen')
                            elif k in ['header_table']:  # only records fitting criteria (zidx) are copied
                                do_copy(fd,f,k,zidx,'hdrlen')

                            else: # groups that are simply copied
                                #print(k)
                                for v in f[k].keys():
                                    f.copy(f[k][v],fd[k],name=v,without_attrs=True)
                                for a in f[k].attrs.keys():
                                    fd[k].attrs[a]=f[k].attrs[a]
                                for v in f[k].keys():
                                    #print(k,v)
                                    #fd[k].create_dataset_like(v,f[k][v])
                                    for a in f[k][v].attrs.keys():
                                        print(k,v,a)
                                        if a not in ['DIMENSION_LIST','CLASS']:
                                            fd[k][v].attrs[a]=f[k][v].attrs[a]
                                    l=0
                                    for d in f[k][v].dims:
                                        if len(d)>0:
                                            fd[k][v].dims[l].attach_scale(fd[k][f[k][v].dims[l][0].name])
                                        l+=1
                        else:  # add vars in root group
                            #print(k)
                            fd.create_dataset_like(k,f[k])
                            fd[k][:]=f[k][:]
                            for a in f[k].attrs.keys():
                                    if a not in ['DIMENSION_LIST','CLASS']:
                                        fd[k].attrs[a]=f[k].attrs[a]
                    for k in f.keys():  # attach dims in root group
                        if not isinstance(f[k],h5py.Group):
                            l=0
                            for d in f[k].dims:
                                if len(d)>0:
                                    fd[k].dims[l].attach_scale(fd[f[k].dims[l][0].name])
                                l+=1
                    
                    #print(fd)

            print(time.time()-t)

        return dfile,''
    else:
       return rfile,'No station ID (statid) specified'
    
''' Main routine for parsing the CDS request'''
def oprocess(rvars):

    print('process: ',rvars)
    error=''
    rfile=''
    vdict={'temperature':2,'uwind':3,'vwind':4}
    rvkeys=rvars.keys()
    #cost=calculate_cost(rvars) # estimate size of output file
    if 'statid' in rvkeys:
        rfile='applications/eua/static/1/chera5.conv._'+rvars['statid']+'.nc'
        rfile='/fio/srvx7/leo/python/CEUAS/CEUAS/public/harvest/data/tables/'+'chera5.conv._'+rvars['statid']+'.nc'
        rfile=os.path.expandvars('$EUA_ROOT/subdaily/v0.1/source/ERA5_1/obs/0-20000-0-'+rvars['statid']+
                                 '/eua_subdaily_v0.1_source_ERA5_1_obs_0-20000-0-'+rvars['statid']+'_t.nc')
        print(rfile)
        if len(rvkeys)>1:
            rvdict=copy.copy(rvars)
            del rvdict['statid']
            rvdk=rvdict.keys()
            for k in rvdk:
                try:
                    rvdict[k]=eval(rvdict[k])
                    print(k,type(rvdict[k]),rvdict[k])
                except:
                    print('could not evaluate '+rvdict[k])
            for k in rvdk:
                if type(rvdict[k]) is not list:
                    if rvdict[k] not in vdict.keys():
                        rvdict[k]=[rvdict[k],rvdict[k]]
                    else:
                        rvdict[k]=[vdict[rvdict[k]],vdict[rvdict[k]]]
                        
            with h5py.File(rfile, 'r') as fs:
            
                di=fs['dateindex']
                if 'date' in rvdk:
                    istart=numpy.where(di[0,:]>=rvdict['date'][0])[0][0]
                    istop=numpy.where(di[0,:]<=rvdict['date'][1])[0][-1]
                ichunk=[di[1,istart],di[2,istop]]
                mask=numpy.ones(ichunk[1]-ichunk[0],dtype=numpy.bool)
                criteria={'date':'era5fb/date@hdr','variable':'era5fb/varno@body','level':'era5fb/vertco_reference_1@body'}
                ck=criteria.keys()
                t=time.time()
                for ck,cv in criteria.items():
                    if ck in rvdk:
                        mask=numpy.logical_and(mask,fs[cv][ichunk[0]:ichunk[1]]>=rvdict[ck][0])
                        mask=numpy.logical_and(mask,fs[cv][ichunk[0]:ichunk[1]]<=rvdict[ck][1])
                        
                print('mask:',time.time()-t)
                dfile='applications/eua/static/tmp/dest.nc'
                with h5py.File(dfile, 'w') as fd:
                    fsk=fs.keys()
                    for f in fsk:
                        if isinstance(fs[f],h5py.Group):
                            fd.create_group(f)
                            fsgk=fs[f].keys()
                            for fsd in fsgk:
                                print('mask:',fs[f][fsd].name,fs[f][fsd].shape,time.time()-t)
        #                        print(fs[f][fsd].name,fs[f][fsd].shape)
                                if 'expver' in fs[f][fsd].name:
                                    print('x')
                                if fs[f][fsd].shape[0]==di[2,-1]+1:
                                    if fs[f][fsd].ndim==1:
                                        z=fs[f][fsd][ichunk[0]:ichunk[1]]
                                        fd.create_dataset(fs[f][fsd].name, data=z[mask],compression=fs[f][fsd].compression,
                                                  compression_opts=fs[f][fsd].compression_opts) 
                                        #for k in fs[f][fsd].attrs.keys():
                                            #print(f,fsd,k,fs[f][fsd].attrs[k])
                                            #if type(fs[f][fsd].attrs[k]) is str:
                                                #fd[f][fsd].attrs.create(k,str.encode(fs[f][fsd].attrs[k]))
                                            #else:
                                                #fd[f][fsd].attrs.create(k,fs[f][fsd].attrs[k])
                                            #fd[f][fsd].attrs[k]=fs[f][fsd].attrs[k]
                                    else:
                                        #idx=numpy.where(mask)
                                        tt=time.time()
        #                                zz=numpy.zeros_like(fs[f][fsd])
                                        z=fs[f][fsd][ichunk[0]:ichunk[1],:]
                                        fd.create_dataset(fs[f][fsd].name, data=z[mask,:],compression=fs[f][fsd].compression,
                                                  compression_opts=fs[f][fsd].compression_opts) 
                                        #for k in fs[f][fsd].attrs.keys():
                                            #fd[f][fsd].attrs[k]=fs[f][fsd].attrs[k]
                                        print(time.time()-tt)
                                        print('x')
                                        
                                else:
                                    fd.create_dataset(fs[f][fsd].name, data=fs[f][fsd],compression=fs[f][fsd].compression,
                                                  compression_opts=fs[f][fsd].compression_opts) 
                        else:
                            print(fs[f].name,fs[f].shape)
                            tt=time.time()
                            y=fs['era5fb/date@hdr'][ichunk[0]:ichunk[1]][mask]
                            print(time.time()-tt)
                            z=find_dateindex(y,numpy.unique(y))
                            print('fd:',time.time()-tt)
                            #z=find_dateindex(y)
                            fd.create_dataset(f, data=z,compression=fs[f].compression,
                                                  compression_opts=fs[f].compression_opts) 
                            print('fd:',time.time()-tt)
                       
                
            #fs.close()
            #fd.close()

            return dfile,''

        else:
            return rfile,''
    else:
       return rfile,'No station ID (statid) specified'
    

    
if __name__ == '__main__':


    os.chdir(os.path.expandvars('$HOME/python/web2py'))
    print(('called with ',sys.argv[1],sys.argv[2]))
    rvars=eval(sys.argv[2])
    try:
        
        if '[0' in rvars['statid']:
            rvars['statid']="['"+rvars['statid'][1:7]+"']"
    except:
        
        pass
    
    df=cdmexplainer(rvars)   
    print(df)
    
    #rfile,error=process(rvars)
    #print(rfile,error)
    
    
