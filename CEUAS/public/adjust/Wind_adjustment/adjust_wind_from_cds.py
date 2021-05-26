#!/usr/bin/env python
import traceback
import sys,glob
import os.path

import numpy as np
from datetime import date
import netCDF4
import time
from numba import njit
sys.path.append('../rasotools/')
from utils import *
from anomaly import *
import matplotlib.pylab as plt
import scipy.stats
import f90nml
import xarray as xr
sys.path.append('../../cds-backend/code/')
import cds_eua3 as eua
import urllib3
import json
import h5py
import copy
from functools import partial
from multiprocessing import Pool

#@njit
def add_winddirbias(xyzu,xyzv,xyzd,xyzt,press,atime0,adj,adjpress):
    adjd=np.full_like(xyzv,np.nan)
    adju=np.full_like(xyzv,np.nan)
    adjv=np.full_like(xyzv,np.nan)
    adjbu=np.full_like(xyzv,np.nan)
    adjbv=np.full_like(xyzv,np.nan)
    ff=np.sqrt(xyzu**2+xyzv**2)

    idtold=0
    for it in range(atime0.shape[0]):
        if it==atime0.shape[0]-1:
            idt=press.shape[0]
        else:
            idt=np.searchsorted(xyzt,atime0[it+1])
        adjd[idtold:idt]=np.nanmean(adj[:,:,it])
        
        adju[idtold:idt]=np.cos((270-xyzd[idtold:idt])*np.pi/180)*ff[idtold:idt]
        adjv[idtold:idt]=np.sin((270-xyzd[idtold:idt])*np.pi/180)*ff[idtold:idt]
        adju[idtold:idt]=np.cos((270-xyzd[idtold:idt]+adjd[idtold:idt])*np.pi/180)*ff[idtold:idt]-adju[idtold:idt]
        adjv[idtold:idt]=np.sin((270-xyzd[idtold:idt]+adjd[idtold:idt])*np.pi/180)*ff[idtold:idt]-adjv[idtold:idt]

        idtold=idt
        #print(it,idt)
            

    return adjd,adju,adjv

@njit
def nanargmax(tsamax):
    
    absmax=0.
    argabsmax=0
    for it in range(tsamax.shape[0]):
        if tsamax[it]==tsamax[it]:
            if tsamax[it]>absmax:
                absmax=tsamax[it]
                argabsmax=it
    return argabsmax

@njit
def select_breaks(tsa,tsaint,thresh):
    breakidx=np.full(20,-1)
    
    tsamax=np.zeros(tsa.shape[2])
    tsamean=np.zeros(tsa.shape[2])
    tsacount=np.zeros(tsa.shape[2])
    for ih in range(tsa.shape[0]):
        for ip in range(tsa.shape[1]):
            for it in range(tsa.shape[2]):
                if tsa[ih,ip,it]==tsa[ih,ip,it]:
                    if tsa[ih,ip,it]>tsamax[it]:
                        tsamax[it]=tsa[ih,ip,it]
                    tsamean[it]+=tsa[ih,ip,it]
                    tsacount[it]+=1
                    
    for it in range(tsa.shape[2]):
        if tsacount[it]>0:
            tsamean[it]/=tsacount[it]
    #tsamax[:]=tsamean
    argabsmax=nanargmax(tsamax)
    #absmax=np.nanargmax(tsamax)
    j=0
    i0=np.int64(0)
    while(tsamax[argabsmax]>thresh) and j<breakidx.shape[0]:
        breakidx[j]=argabsmax
        istart=argabsmax-tsaint//3
        if istart<0:
            istart=0
        istop=argabsmax+tsaint//3
        if istop>tsamax.shape[0]:
            istop=tsamax.shape[0]
        print(argabsmax,tsamax[argabsmax])
        tsamax[istart:istop]=0.
        argabsmax=nanargmax(tsamax)
        j+=1
        
    
    return breakidx[:j+1]

def calc_dirshifts(ddeps,cdict,breakidx,tsaint,delta):
    
    dirshifts=np.zeros(breakidx.shape[0])
    ff=np.sqrt(cdict['u']['xrdq']['uwind'].values**2+cdict['v']['xrdq']['vwind'].values**2)
    #strongddeps=copy.deepcopy(ddeps)
    ddeps[ff<0.5]=np.nan
    ds=cdict['d']['xrdq']['winddirection'].values[:]
    #ds[ds>360]=np.nan
    #ds[ds<0]=np.nan
    
    ufg=cdict['u']['xrdq']['uwind'].values-cdict['u']['xrdq']['era5_fgdep'].values
    vfg=cdict['v']['xrdq']['vwind'].values-cdict['v']['xrdq']['era5_fgdep'].values
    cdict['u']['xrdq']['uwindbias']=cdict['u']['xrdq']['uwind'].copy(deep=True)
    cdict['v']['xrdq']['vwindbias']=cdict['v']['xrdq']['vwind'].copy(deep=True)
    #cdict['u']['xrdq']['uwind'].values[:]=np.cos((270-ds[:])*np.pi/180)*ff[:]
    #cdict['v']['xrdq']['vwind'].values[:]=np.sin((270-ds[:])*np.pi/180)*ff[:]
    #x=cdict['u']['xrdq']['uwind'].values[0,10,:]-cdict['u']['xrdq']['uwindbias'].values[0,10,:]

    cdict['d']['xrdq']['directionbias']=cdict['d']['xrdq']['winddirection'].copy(deep=True)
    cdict['d']['xrdq']['directionbias'].values[:]=0.
    #cdict['u']['xrdq']['uwindbias'].values[:]=0.
    #cdict['v']['xrdq']['vwindbias'].values[:]=0.
    for bi in range(breakidx.shape[0]-1,0,-1):
        istart=np.max((breakidx[bi-1],breakidx[bi]-tsaint))
        istop=np.min((breakidx[bi]+tsaint,ddeps.shape[2]))
        dirshifts[bi]=np.nanmean(ddeps[:,:,istart:breakidx[bi]-delta])-np.nanmean(ddeps[:,:,breakidx[bi]+delta:istop]) #:istop
        unc=np.sqrt(0.5*(np.nanstd(ddeps[:,:,istart:breakidx[bi]-delta])**2+np.nanstd(ddeps[:,:,breakidx[bi]+delta:istop])**2))
        count=np.sum(~np.isnan(ddeps[:,:,istart:istop]))
        if dirshifts[bi]!=dirshifts[bi]:
            dirshifts[bi]=0.
        adir=abs(dirshifts[bi])
        if adir>3.0 and adir>1.96*unc/np.sqrt(count):
            ddeps[:,:,:breakidx[bi]]-=dirshifts[bi]
            ds[:,:,:breakidx[bi]]+=dirshifts[bi]
            ds[ds>360.]-=360.
            ds[ds<0.]+=360.
            print('adjusting direction by{:5.2f}'.format(dirshifts[bi]),
                  'degrees at {:5.2f}'.format(1900+cdict['u']['xrdq']['datum'].values[breakidx[bi]]/365.25))
        else:
            dirshifts[bi]=0.
            continue
        # now calculate the increments
        cdict['d']['xrdq']['directionbias'].values[:,:,:breakidx[bi]]+=dirshifts[bi]
        cdict['u']['xrdq']['uwind'].values[:,:,:breakidx[bi]]=np.cos((270-ds[:,:,:breakidx[bi]])*np.pi/180)*ff[:,:,:breakidx[bi]]
        cdict['v']['xrdq']['vwind'].values[:,:,:breakidx[bi]]=np.sin((270-ds[:,:,:breakidx[bi]])*np.pi/180)*ff[:,:,:breakidx[bi]]
    cdict['u']['xrdq']['uwindbias'].values[:]-=cdict['u']['xrdq']['uwind'].values[:]
    cdict['v']['xrdq']['vwindbias'].values[:]-=cdict['v']['xrdq']['vwind'].values[:]                
        
        
    cdict['u']['xrdq']['era5_fgdep'].values=cdict['u']['xrdq']['uwind'].values-ufg
    cdict['v']['xrdq']['era5_fgdep'].values=cdict['v']['xrdq']['vwind'].values-vfg
    
    cdict['d']['xrdq']['winddirection'].values[:]=ds
    
    return dirshifts

@njit
def getpindex(press,plist):
    oindex=np.empty(len(press),dtype='int')
    pindex=np.empty(len(press),dtype='int')
    reverseindex=np.full(110000,-1,dtype='int')
    for i in range(len(plist)):
        reverseindex[plist[i]]=i
    l=0
    for i in range(press.shape[0]):
        if press[i]>0 and press[i]<110000:
            
            ri=reverseindex[press[i]]
            if ri!=-1:
                oindex[l]=i
                pindex[l]=ri
                l+=1
        else:
            pass
            #print(i,press[i])
    
    return oindex[:l],pindex[:l]
    

def homogenize_winddir(via_backend=False, fns=[]):
    
    lplot=False
#     with open(os.path.expanduser('~leo/python/hug2/config/active.json')) as f:
#         active=json.load(f)
#     ids=list(active.keys())
#     lats=np.asarray([active[x][2] for x in active.keys()])   
#     lons=np.asarray([active[x][3] for x in active.keys()])
#     starts=np.asarray([active[x][0] for x in active.keys()])   
#     stops=np.asarray([active[x][1] for x in active.keys()])
#     l=0
#     for i in range(lats.shape[0],lats.shape[0]):
#         idx=np.where(np.logical_and(np.abs(lats[i]-lats)<0.1,np.abs(lons[i]-lons)<0.1))[0]
#         if len(idx)>1:
#             fak=86400*365.25
#             print('duplicate {:s},{:s},{:4.0f},{:4.0f},{:4.0f},{:4.0f}'.format(ids[idx[0]],ids[idx[1]],
#                                                                                starts[idx[0]]/fak,starts[idx[1]]/fak,stops[idx[0]]/fak,stops[idx[1]]/fak))
#             try:
                
#                 with h5py.File('/raid60/scratch/leo/scratch/converted_v5/'+ids[idx[0]]+'_CEUAS_merged_v1.nc','r') as f:
#                     with h5py.File('/raid60/scratch/leo/scratch/converted_v5/'+ids[idx[1]]+'_CEUAS_merged_v1.nc','r') as g:
#                         try:
#                             print(f['observations_table']['latitude'][0],f['observations_table']['longitude'][0],
#                                   g['observations_table']['latitude'][1],g['observations_table']['longitude'][1])
#                             l+=1
#                         except:
                            
#                             print('table read error')
#             except:
#                 print('file open error')
                
#     print(l,' duplicates')
            
    
#     http = urllib3.PoolManager()
#     r = http.request('GET', 'http://early-upper-air.copernicus-climate.eu/statlist/?mindate=1900-01-01&enddate=2020-12-31')
# #     r = http.request('GET', 'http://srvx8.img.univie.ac.at:8002/statlist/?mindate=1900-01-01&enddate=2020-12-31')
#     fns=r.data.split(b'\n')
#     for i in range(len(fns)):
#         fns[i]=fns[i].split(b',')[0].decode()
    opath=os.path.expandvars('./')
    os.chdir(opath)
#     #fns=glob.glob('0?????/')
#     #fns=[fns[fns.index('0-20000-0-35229')]]
    tt=time.time()
    print(os.getcwd())
    if via_backend:
        func=partial(homogenize_station,opath,via_backend)
        adjusted=list(map(func,fns))
    else:    
        p=Pool(40)
        func=partial(homogenize_station,opath,via_backend)
        adjusted=list(p.map(func,fns))
        
    print(time.time()-tt,len(adjusted))

def homogenize_station(opath,via_backend,fnf):
    
    tt=time.time()
    lplot=False
    cdict={'d':{'wind_from_direction':'winddirection','obs_minus_bg':'era5_fgdep','bias_estimate':'bias_estimate','lat':'lat','lon':'lon','hours':'hours'},
           'u':{'ua':'uwind','obs_minus_bg':'era5_fgdep','bias_estimate':'bias_estimate','lat':'lat','lon':'lon','hours':'hours'},
           'v':{'va':'vwind','obs_minus_bg':'era5_fgdep','bias_estimate':'bias_estimate','lat':'lat','lon':'lon','hours':'hours'}}

    varnos={'d':106,'u':104,'v':105}
    plist=[1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000,92500,100000]
    fnu=[]
    fnd=[]
    dimsize_errors=0
    failedfiles=[]
    #for fnf in fns[:]:
    fn=fnf[-5:]
    
    prefix='0'
    if fn in fnu:
        print('duplicate '+fnf+', incrementing leading zero to 1')
        prefix='1'
    fnu.append(fn)
    fo=opath+prefix+fn+'/feedbackmergedwinddir'+prefix+fn+'.nc'
    try:
        
        mt=os.path.getmtime(fo)
        ts=time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(mt))
        #if ts.split(' ')[0][-2:]=='11':
            #continue
    except:
        pass
    
    variables=dict(zip(['d','u','v'],['wind_direction','u_component_of_wind','v_component_of_wind']))
    shortvariables=dict(zip(['d','u','v'],['winddirection','uwind','vwind']))
    
    if via_backend:
        print('cds data')
        
        try:        
            #data=eua.vm_request_wrapper({'variable': 'temperature', 'statid': fn, 'date':['20100101','20100102']}, 
                                        #overwrite=True,vm_url='http://srvx8.img.univie.ac.at:8002')
#             ifile=glob.glob('/raid60/scratch/leo/scratch/converted_v5/'+fnf+'_CEUAS_merged_v1.nc')[0]
            #try:
                
                #with h5py.File(ifile,'r') as f:
                    #if 'wind_direction_bias_estimate' in f['advanced_homogenisation'].keys():
                        #continue
            #except:
                #pass
            dir_file = glob.glob('downloaded/wind_downloaded_*/*'+fn+'*direction*')[0]
            u_file = glob.glob('downloaded/wind_downloaded_*/*'+fn+'*eastw*')[0]
            v_file = glob.glob('downloaded/wind_downloaded_*/*'+fn+'*north*')[0]
            data=dict(zip(['d','u','v'],[eua.CDMDataset(dir_file), eua.CDMDataset(u_file), eua.CDMDataset(v_file)]))

        except Exception as e:
            print(e)
            return
        
        try:
            
            for k,v in data.items():
                if 'cube' in cdict[k].keys():   
                    del cdict[k]['cube']
                    del cdict[k]['xrdq']
                cdict[k]['cube']=v.read_data_to_3dcube(list(cdict[k].keys()))
                dql={}
                dq=cdict[k]['cube']
                for kk in dq.keys():
                ###dq[k].rename_dims({'plev':'pressure'})
                    dql[cdict[k][kk]]=dq[kk].rename(cdict[k][kk])
                cdict[k]['xrdq']=xr.Dataset(dql).rename_dims({'press':'pressure','datum':'time'})
                cdict[k]['xrdq'].attrs['unique_source_identifier']=fnf

        except Exception as e:
            print(e)
            return
    
#     else:
#         try:
            
#             ifile=glob.glob('downloaded/wind_downloaded_*/'+fnf+'_CEUAS_merged_v1.nc')[0]
        
#             data=eua.CDMDataset(ifile)
#         except:
#             return
    
#         ndata=copy.deepcopy(cdict)
#         for para in ['d','u','v']:
#             #xyz = data.read_observed_variable(varnos[para], 
#                                                #pressure_level=[1000,2000,3000,5000,7000,10000,15000,20000,25000,
#                                                                  #30000,40000,50000,70000,85000,92500,100000],
#                                                #return_xarray=True,date_time_in_seconds=True)
#             try:
                
#                 istart=data['recordindices'][str(varnos[para])][0]
#                 istop=data['recordindices'][str(varnos[para])][-1]
#                 #ndata[para]={}
#                 ndata[para]['z_coordinate']=data['observations_table']['z_coordinate'][istart:istop].astype(np.int)
#                 ndata[para]['istart']=istart
#                 ndata[para]['istop']=istop
                
#             except:
#                 print(fnf,'no wind direction data')
#                 return
            
            
#             oindex,pindex=getpindex(ndata[para]['z_coordinate'],np.array(plist,dtype=np.int))
#             if len(pindex)==0:
#                 print(fnf,'no pressure on standard levels')
#                 return
            
#             ndata[para]['observation_value']=data['observations_table']['observation_value'][istart:istop][oindex]
#             ndata[para]['obs_minus_bg']=data['era5fb']['fg_depar@body'][istart:istop][oindex]
#             ndata[para]['date_time']=data['observations_table']['date_time'][istart:istop][oindex]
            
#             out=eua.daysx2(ndata[para]['date_time'],pindex,len(plist),ndata[para]['observation_value'])
    
#             ndata[para]['xrdq']=xr.Dataset()
#             ndims={'hour':2,'press':len(plist),'datum':len(out[5])}
#             ncoords={'hour':np.asarray([0,12]),'press':np.array(plist)/100.,'datum':out[5]}
#             ndata[para]['xrdq'][shortvariables[para]]=xr.DataArray(out[0],dims=ndims,coords=ncoords)
            
#             out=eua.daysx2(ndata[para]['date_time'],pindex,len(plist),ndata[para]['obs_minus_bg'])
#             ndata[para]['xrdq']['era5_fgdep']=xr.DataArray(out[0],dims=ndims)
#             ndata[para]['xrdq']['bias_estimate']=xr.DataArray(np.zeros_like(out[0]),dims=ndims)
            
#             ndata[para]['xrdq']['datum'].attrs['units']='days since 1900-01-01 00:00:00'
#             ndata[para]['xrdq']['press'].attrs['units']='hPa'
            
#             ndata[para]['xrdq']['lat']=xr.DataArray(np.asarray([data['observations_table']['latitude'][istop]]),name='lat',dims=('station'))
#             ndata[para]['xrdq']['lon']=xr.DataArray(np.asarray([data['observations_table']['latitude'][istop]]),name='lon',dims=('station'))
#             ndata[para]['xrdq']['hours']=xr.DataArray(out[6],name='hours',dims=('hour','datum'))
#             ndata[para]['xrdq']=ndata[para]['xrdq'].rename_dims({'press':'pressure','datum':'time'})
#             ndata[para]['xrdq'].attrs['unique_source_identifier']=fnf
#             #print(ndata[para]['xrdq'][shortvariables[para]])
        
#         cdict=ndata    
    try:
        
        tsa=np.full(cdict['d']['xrdq']['era5_fgdep'].shape,np.nan,dtype=np.float32)
        ddeps=np.full(cdict['d']['xrdq']['era5_fgdep'].shape,np.nan,dtype=np.float32)
        #ds=np.full_like(ddeps,np.nan)
        snhtparas=np.asarray([1460,650,30])
        index=np.zeros(tsa.shape[-1],dtype=np.int)
        count=np.zeros(tsa.shape[-1],dtype=np.int)
        tmean=np.zeros(tsa.shape[-1],dtype=np.float32)
        tsquare=np.zeros(tsa.shape[-1],dtype=np.float32)
        for i in range(1):
            
    
            for ih in range(cdict['d']['xrdq']['era5_fgdep'].shape[0]):
                for ip in range(cdict['d']['xrdq']['era5_fgdep'].shape[1]-2):
                    #v=-cdict['d']['xrdq']['era5_fgdep'].values[ih,ip,:]
                    vu=-cdict['u']['xrdq']['era5_fgdep'].values[ih,ip,:]
                    vv=-cdict['v']['xrdq']['era5_fgdep'].values[ih,ip,:]
                    o=cdict['d']['xrdq']['winddirection'].values[ih,ip,:]
                    ou=cdict['u']['xrdq']['uwind'].values[ih,ip,:]
                    ov=cdict['v']['xrdq']['vwind'].values[ih,ip,:]
                    ff=np.sqrt(ou**2+ov**2)
                    o2=270-np.arctan2(ov,ou)*180/np.pi
                    fu=ou-vu
                    fv=ov-vv
                    fd=270-np.arctan2(fv,fu)*180/np.pi
                    v=o2-fd
                    v[v>180]-=360
                    v[v<-180]+=360
                    v[np.abs(v)>90.]=np.nan # larger deviations than 90 degrees are highly unlikely to come from wrong north alignment, are therefore discarded
                    
                    #qs=np.nanquantile(v,[0.005,0.995])
                    #idx=np.where(np.logical_or(v<qs[0],v>qs[1]))
                    #print(qs,idx)
                    #v[idx]=np.nan
                    #o[idx]=np.nan
                    #ds[ih,ip,:]=o2[:]
                    v[ff<2.0]=np.nan
                    ddeps[ih,ip,:]=v[:]
                    tsa[ih,ip,:]=np.nan
                    snhtmov2(v, tsa[ih,ip,:], snhtparas, index, count, tmean, tsquare)    
                    if lplot and ip==10:
                        
                        plt.subplot(3,1,1)
                        plt.plot(cdict['d']['xrdq'].datum.values[:]/365.25,tsa[ih,ip,:],
                                 label='{} {} {:5.0f}'.format(ih,ip,np.nanmax(tsa[ih,ip,:])))
                        plt.subplot(3,1,2)
                        plt.plot(cdict['d']['xrdq'].datum.values[:]/365.25,rmeanw(vu,30),
                                 label='{} {} {:5.2f}'.format(ih,ip,np.nanstd(vu)))
                    
                    #hilf=xrdq['bias_estimate'].values[ih,ip,:]
                    #hilf[np.isnan(hilf)]=0.
                    #xrdq['era5_fgdep'].values[ih,ip,:]=-v-hilf
                        
                        # ignore missing fgdep, bias_estimate
                        #idx=np.where(np.logical_and(np.isnan(v),~np.isnan(o)))
                        ##print(len(idx[0]))
                        #xrdq['era5_fgdep'].values[ih,ip,idx]=0.
        
            breakpoints=select_breaks(tsa,1460,100.)
            breakpoints.sort()
            dirshifts=calc_dirshifts(ddeps,cdict,breakpoints,2*1460,160)
            if lplot:
                
                plt.subplot(3,1,3)
                plt.plot(cdict['d']['xrdq'].datum.values[:]/365.25,cdict['d']['xrdq']['directionbias'].values[:,10,:].T,
                         label='direction bias')
                
                plt.legend()
                plt.show()
            try:
                os.mkdir(opath+prefix+fn)
            except:
                pass
    except Exception as e:
        print(e)
        failedfiles.append(fnf+' early')
        
    try:
        
        cdict['d']['xrdq']['uwindbias']=cdict['u']['xrdq']['uwindbias']
        cdict['d']['xrdq']['vwindbias']=cdict['v']['xrdq']['vwindbias']
        cdict['d']['xrdq'].drop_vars('bias_estimate')
        cdict['d']['xrdq'].to_netcdf(path=fo, format='NETCDF4_CLASSIC')
    except Exception as e:
        print(e)
        dimsize_errors+=1
        failedfiles.append(fnf)
        return
    
    # Daten schreiben neue Variable monkey in neuer gruppe adjust

    if not via_backend:
        
        try:
            ranges={'d':[],'u':[],'v':[]}
            biasnames={'d':'directionbias','u':'uwindbias','v':'vwindbias'}

            data = eua.CDMDataset(ifile)
            xyz = data.read_observed_variable(varnos['u'], return_xarray=True,date_time_in_seconds=True)
            
            for k,v in varnos.items():
                ranges[k]=[data['recordindices'][str(v)][0],data['recordindices'][str(v)][-1]]
            datau=data['observations_table']['observation_value'][ranges['u'][0]:ranges['u'][-1]]
            datav=data['observations_table']['observation_value'][ranges['v'][0]:ranges['v'][-1]]
            datad=data['observations_table']['observation_value'][ranges['d'][0]:ranges['d'][-1]]
            
            adjd,adju,adjv=add_winddirbias(datau,
                                datav,datad,
                                data['observations_table']['date_time'][ranges['u'][0]:ranges['u'][-1]],
                                data['observations_table']['z_coordinate'][ranges['u'][0]:ranges['u'][-1]],
                                cdict['d']['xrdq']['datum'].values[:]*86400,
                                cdict['d']['xrdq']['directionbias'].values[:],
                                cdict['d']['xrdq']['press'].values[:])
            
            raggedd={'d':adjd,'u':adju,'v':adjv}
            for k in raggedd.keys():
            
                xyz.values[:]=raggedd[k]
                if k=='d':
                    
                
                    data.write_observed_data('wind_bias_estimate',
                                         ragged=xyz,  # input data
                                         varnum=varnos[k],  # observed_variable to be aligned with
                                         group='advanced_homogenisation',   # name of the new group
                                         data_time='date_time',  # named datetime coordinate
                                         data_plevs='z_coordinate',  # named pressure coordinate
                                         attributes={'version':version}
                                        )
                    data.close()
                else:
                    with h5py.File(ifile,'r+') as data:
                        
                        data.file['advanced_homogenisation']['wind_bias_estimate'][cdict[k]['istart']:cdict[k]['istop']]=raggedd[k]
                    
            print('write:',time.time()-tt)
            
            print('wrote '+fo)
        except Exception as e:
            print('writing back to merged file failed')
            failedfiles.append(fnf+'-merged')
    
    return fnf
    #print(dimsize_errors,failedfiles)



if __name__ == "__main__":


    plt.rcParams['lines.linewidth'] = 3
    tt=time.time()
    version='1.0'
#     via_backend=True
#     homogenize_winddir(via_backend, fns=['0-20001-0-11035', '0-20001-0-10393', '0-20001-0-70219', '0-20000-0-91413', '0-20000-0-35229', '0-20000-0-68994'])
    

    print((time.time()-tt))
