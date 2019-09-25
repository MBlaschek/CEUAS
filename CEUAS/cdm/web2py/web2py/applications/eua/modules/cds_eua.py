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
        with xarray.open_dataset(rfile,group=rvdict['group']) as x: 
            print(rfile,x)
            xk=list(x.keys())
            if rvdict['group']=='era5fb':
                xk=xk+list(x.dims)
            del rvdict['group']
            y=dict()
            rvdk=list(rvdk)[0]
            ki=xk.index(rvdk)
            xk.insert(0,xk.pop(ki))
            if rvdict:
                idx=numpy.where(x[rvdk].values.astype(int)==int(rvdict[rvdk]))[0]
                for k in x.keys():
                    print(x[k].values.shape)
                    if type(x[k].values[0]) in [bytes,numpy.bytes_]:
                        #y[k]=x[k].values.view('S80')[idx[0]]
                        y[k]=x[k].values[idx[0]]
                        print(y[k])
                        y[k]=numpy.asarray([y[k]],dtype='unicode')
                    else:
                        y[k]=numpy.asarray([x[k].values[idx[0]]])
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
    
    
