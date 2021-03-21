#!/usr/bin/env python
import traceback
import sys,glob
import os.path

import numpy as np
from datetime import date
import netCDF4
import time
from numba import njit
from rasotools.utils import *
from rasotools.anomaly import *
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

#@njit
def add_winddirbias(xyzu,xyzv,xyzt,press,atime0,adj,adjpress):
    adjd=np.full_like(xyzv,np.nan)
    adju=np.full_like(xyzv,np.nan)
    adjv=np.full_like(xyzv,np.nan)
    ff=np.sqrt(xyzu**2+xyzv**2)

    idtold=0
    for it in range(atime0.shape[0]):
        if it==atime0.shape[0]-1:
            idt=press.shape[0]
        else:
            idt=np.searchsorted(xyzt,atime0[it+1])
        adjd[idtold:idt]=np.nanmean(adj[:,:,it])
        adju[idtold:idt]=np.cos((270-adjd[idtold:idt])*np.pi/180)*ff[idtold:idt]
        adjv[idtold:idt]=np.sin((270-adjd[idtold:idt])*np.pi/180)*ff[idtold:idt]

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
    ds=cdict['d']['xrdq']['winddirection'].values
    ds[ds>360]=np.nan
    ds[ds<0]=np.nan
    
    ufg=-(cdict['u']['xrdq']['uwind'].values-cdict['u']['xrdq']['era5_fgdep'].values)
    vfg=-(cdict['v']['xrdq']['vwind'].values-cdict['v']['xrdq']['era5_fgdep'].values)
    cdict['u']['xrdq']['uwindbias']=cdict['u']['xrdq']['uwind'].copy(deep=True)
    cdict['v']['xrdq']['vwindbias']=cdict['v']['xrdq']['vwind'].copy(deep=True)
    cdict['d']['xrdq']['directionbias']=cdict['d']['xrdq']['winddirection'].copy(deep=True)
    cdict['d']['xrdq']['directionbias'].values[:]=0.
    cdict['u']['xrdq']['uwindbias'].values[:]=0.
    cdict['v']['xrdq']['vwindbias'].values[:]=0.
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
        cdict['u']['xrdq']['uwindbias'].values[:,:,:breakidx[bi]]-=np.cos((270-ds[:,:,:breakidx[bi]])*np.pi/180)*ff[:,:,:breakidx[bi]]
        cdict['v']['xrdq']['vwindbias'].values[:,:,:breakidx[bi]]-=np.sin((270-ds[:,:,:breakidx[bi]])*np.pi/180)*ff[:,:,:breakidx[bi]]
        cdict['u']['xrdq']['uwind'].values[:,:,:breakidx[bi]]-=cdict['u']['xrdq']['uwindbias'].values[:,:,:breakidx[bi]]
        cdict['v']['xrdq']['vwind'].values[:,:,:breakidx[bi]]-=cdict['v']['xrdq']['vwindbias'].values[:,:,:breakidx[bi]]                
        
        
    cdict['u']['xrdq']['era5_fgdep'].values=cdict['u']['xrdq']['uwind'].values-ufg
    cdict['v']['xrdq']['era5_fgdep'].values=cdict['v']['xrdq']['vwind'].values-vfg
    
    cdict['d']['xrdq']['winddirection'].values[:]=ds
    
    return dirshifts

def homogenize_winddir():
    
    lplot=False
    with open(os.path.expanduser('~leo/python/hug2/config/active.json')) as f:
        active=json.load(f)
    ids=list(active.keys())
    lats=np.asarray([active[x][2] for x in active.keys()])   
    lons=np.asarray([active[x][3] for x in active.keys()])
    starts=np.asarray([active[x][0] for x in active.keys()])   
    stops=np.asarray([active[x][1] for x in active.keys()])
    l=0
    for i in range(lats.shape[0],lats.shape[0]):
        idx=np.where(np.logical_and(np.abs(lats[i]-lats)<0.1,np.abs(lons[i]-lons)<0.1))[0]
        if len(idx)>1:
            fak=86400*365.25
            print('duplicate {:s},{:s},{:4.0f},{:4.0f},{:4.0f},{:4.0f}'.format(ids[idx[0]],ids[idx[1]],
                                                                               starts[idx[0]]/fak,starts[idx[1]]/fak,stops[idx[0]]/fak,stops[idx[1]]/fak))
            try:
                
                with h5py.File('/raid60/scratch/leo/scratch/converted_v5/'+ids[idx[0]]+'_CEUAS_merged_v1.nc','r') as f:
                    with h5py.File('/raid60/scratch/leo/scratch/converted_v5/'+ids[idx[1]]+'_CEUAS_merged_v1.nc','r') as g:
                        try:
                            print(f['observations_table']['latitude'][0],f['observations_table']['longitude'][0],
                                  g['observations_table']['latitude'][1],g['observations_table']['longitude'][1])
                            l+=1
                        except:
                            
                            print('table read error')
            except:
                print('file open error')
                
    print(l,' duplicates')
            
    
    http = urllib3.PoolManager()
    r = http.request('GET', 'http://srvx8.img.univie.ac.at:8002/statlist/?mindate=1900-01-01&enddate=2020-12-31')
    fns=r.data.split(b'\n')
    for i in range(len(fns)):
        fns[i]=fns[i].split(b',')[0].decode()
    opath=os.path.expandvars('$FSCRATCH/rise/1.0/exp06/')
    os.chdir(opath)
    #fns=glob.glob('0?????/')
    #fns=[fns[fns.index('0-20000-0-35229')]]
    tt=time.time()
    cdict={'d':{'wind_from_direction':'winddirection','obs_minus_bg':'era5_fgdep','bias_estimate':'bias_estimate','lat':'lat','lon':'lon','hours':'hours'},
           'u':{'ua':'uwind','obs_minus_bg':'era5_fgdep','bias_estimate':'bias_estimate','lat':'lat','lon':'lon','hours':'hours'},
           'v':{'va':'vwind','obs_minus_bg':'era5_fgdep','bias_estimate':'bias_estimate','lat':'lat','lon':'lon','hours':'hours'}}
    fnu=[]
    fnd=[]
    dimsize_errors=0
    failedfiles=[]
    for fnf in fns[633:]:
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
        try:        
            #data=eua.vm_request_wrapper({'variable': 'temperature', 'statid': fn, 'date':['20100101','20100102']}, 
                                        #overwrite=True,vm_url='http://srvx8.img.univie.ac.at:8002')
            data=dict(zip(['d','u','v'],eua.vm_request_wrapper({'variable': ['wind_direction','u_component_of_wind','v_component_of_wind'],  'optional':['obs_minus_bg','bias_estimate'],'statid': fn, 'pressure_level':[1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000,92500,100000]}, 
                                        overwrite=True,vm_url='http://srvx8.img.univie.ac.at:8002').values()))
            #udata=eua.vm_request_wrapper({'variable': 'u_component_of_wind', 'optional':['obs_minus_bg','bias_estimate'],'statid': fn, 'pressure_level':[1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000,92500,100000]}, 
                                        #overwrite=True,vm_url='http://srvx8.img.univie.ac.at:8002')
            #vdata=eua.vm_request_wrapper({'variable': 'v_component_of_wind', 'optional':['obs_minus_bg','bias_estimate'],'statid': fn, 'pressure_level':[1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000,92500,100000]}, 
                                        #overwrite=True,vm_url='http://srvx8.img.univie.ac.at:8002')
        except Exception as e:
            print(e)
            continue
        
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
            continue
        
        try:
            
            tsa=np.full(cdict[k]['xrdq']['era5_fgdep'].shape,np.nan,dtype=np.float32)
            ddeps=np.full(cdict[k]['xrdq']['era5_fgdep'].shape,np.nan,dtype=np.float32)
            #ds=np.full_like(ddeps,np.nan)
            snhtparas=np.asarray([1460,650,30])
            index=np.zeros(tsa.shape[-1],dtype=np.int)
            count=np.zeros(tsa.shape[-1],dtype=np.int)
            tmean=np.zeros(tsa.shape[-1],dtype=np.float32)
            tsquare=np.zeros(tsa.shape[-1],dtype=np.float32)
            for i in range(1):
                
        
                for ih in range(cdict[k]['xrdq']['era5_fgdep'].shape[0]):
                    for ip in range(cdict[k]['xrdq']['era5_fgdep'].shape[1]-2):
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
            continue
        
        # Daten schreiben neue Variable monkey in neuer gruppe adjust

        try:
            varnos={'d':106,'u':104,'v':105}
            ranges={'d':[],'u':[],'v':[]}
            biasnames={'d':'directionbias','u':'uwindbias','v':'vwindbias'}
            ifile=glob.glob('/raid60/scratch/leo/scratch/converted_v5/'+fnf+'_CEUAS_merged_v1.nc')[0]
            data = eua.CDMDataset(ifile)
            xyz = data.read_observed_variable(varnos['u'], return_xarray=True,date_time_in_seconds=True)
            
            for k,v in varnos.items():
                ranges[k]=[data['recordindices'][str(v)][0],data['recordindices'][str(v)][-1]]
            datau=data['observations_table']['observation_value'][ranges['u'][0]:ranges['u'][-1]]
            datav=data['observations_table']['observation_value'][ranges['v'][0]:ranges['v'][-1]]
            
            adjd,adju,adjv=add_winddirbias(datau,
                                datav,
                                data['observations_table']['date_time'][ranges['u'][0]:ranges['u'][-1]],
                                data['observations_table']['z_coordinate'][ranges['u'][0]:ranges['u'][-1]],
                                cdict['d']['xrdq']['datum'].values[:]*86400,
                                cdict['d']['xrdq']['directionbias'].values[:],
                                cdict['d']['xrdq']['press'].values[:])
            
            raggedd={'d':adjd,'u':adju,'v':adjv}
            for k in raggedd.keys():
            
                xyz.values[:]=raggedd[k]
                data.write_observed_data(variables[k]+'_bias_estimate',
                                     ragged=xyz,  # input data
                                     varnum=varnos[k],  # observed_variable to be aligned with
                                     group='advanced_homogenisation',   # name of the new group
                                     data_time='date_time',  # named datetime coordinate
                                     data_plevs='z_coordinate'  # named pressure coordinate
                                    )
            print('write:',time.time()-tt)
            
            print('wrote '+fo)
        except Exception as e:
            print('writing back to merged file failed')
            failedfiles.append(fnf+'-merged')
    
    print(dimsize_errors,failedfiles)



if __name__ == "__main__":


    plt.rcParams['lines.linewidth'] = 3
    
    homogenize_winddir()
    
    if len(sys.argv)<2:
        exp='exp06'
    else:
        exp=sys.argv[1]

    try:
        nml = f90nml.read('/home/srvx7/leo/fastscratch/rise/1.0/'+exp+'/'+'radcorpar')
    except:
        pass
    fgdepname=nml['rfpar']['fgdepname']
    #f=netCDF4.Dataset('/home/srvx7/raobcore/v1.5.1/export/ERA5/ERA5bc_RAOBCORE_v1.5_070219.nc','r')
    f=netCDF4.Dataset('/home/srvx7/raobcore/v1.5.1/export/ERA5/ERA5bc_RAOBCORE_v1.5_098646.nc','r')
    f.set_auto_mask(False)
    vals=f.variables['temperatures'][1,3,:]
    offsets=f.variables['datum'][0,:]
    print((vals[:4],offsets[0]))
    print((datetime.date(1900,1,1)+datetime.timedelta(days=int(offsets[0]))))
    mask=vals!=-999.
    #plt.plot(f.variables['datum'][0,mask]/365.25+1900,f.variables['bias'][0,5,mask])
    #plt.show()

    plpath=os.path.expandvars('$RS/plots/')
    try:
        os.mkdir(plpath)
    except:
        pass
    monthdays=[31,28,31,30,31,30,31,31,30,31,30,31]  
    start=1900
    stop=2014
    plevs=np.asarray([10,20,30,50,70,100,150,200,250,300,400,500,700,850,925,1000])
    nlev=plevs.shape[0]

    fmerged = open('/home/srvx7/leo/fastscratch/rise/1.0/'+exp+'/mergedstations.t', 'r').readlines()
    mlist=[]
    for l in fmerged:
        mlist.append(l.split()[2])

    stats=dict()
    lats=list()
    lons=list()
    indexmax=45000
    tidx=calcdays(19000101,(2015-1900)*12)-1
    l=0
    t=time.time()
    #forig=os.popen('ls /home/srvx7/leo/fastscratch/ei6/*/*_t.nc').read().split()
    forig=glob.glob('/home/srvx7/leo/fastscratch/rise/1.0/'+exp+'/*/feedbackmerged*.nc')
    #forig=glob.glob('/home/srvx7/leo/fastscratch/rise/1.0/'+exp+'/[0-9]*/feedbackmerged010046.nc')
    ic=0
    for line in forig:
    #    str=line[9:31]
    #    st=str.split()
        st=line.split('/')[-1][-9:-3]
        if st not in mlist:
            print((st,' nicht im merged'))
        else:
            print((st,'im merged'))
            ic+=1

    base = datetime.date(1900,1,1)

    il=0
    for line in forig:
    #    str=line[9:31]
    #    st=str.split()
        #st=line.split('/')[6]
        st=line.split('/')[-1][-9:-3]
    #    if st<='050000' or st>'090000':
        #if st!='027612':
            #continue
    #    if len(st)==5:
    #        st='0'+st
    #    st='070361'
        fns=['/home/srvx7/leo/fastscratch/rise/1.0/'+exp+'/'+st+'/feedbackmerged'+st+'.nc',
#	    '/home/srvx7/leo/scratch/stream1/RAOBCORE_RICH_v1.3_nc/feedbackmerged'+st[1:6]+'.nc',
'/home/srvx7/leo/fastscratch/rise/1.0/'+exp+'/'+st+'/'+'feedbackglobbincorrsave'+st+'.nc',
            '/home/srvx7/leo/fastscratch/rise/1.0/'+exp+'/'+st+'/'+'feedbackglobbincorrsave_rio24_'+st+'.nc'
            ]
        varlists=[['temperatures',fgdepname,'an_dep'],['rasocorr'],['rasocorr']]
    #    fns=['/home/srvx7/leo/fastscratch/ei6/'+st+'/feedbackmerged'+st+'.nc',]
        i=0
        stats[st]=dict()
        for fn in fns:
            read_temperature(stats[st],fn,tidx,varlist=varlists[i])
            i+=1

        if 'lat' not in list(stats[st].keys()) or 'mdatum' not in list(stats[st].keys()):
            print((st,'no data found,cycling ..'))
            continue
        s=stats[st]

        t2=time.time()
        elev=calc_elevangles(0.,np.asarray([s['lat']]),np.asarray([s['lon']]))
        for y in range(120):
            idx=datetime.date(1900+y,1,1)-base
            idx=idx.days
            elev[0,:,idx:idx+366]=elev[0,:,:366]
        print(('t2:',time.time()-t2))

        bdatum=[]
        try:
            if il==0:
                full=np.empty([2,s['mtemperatures'].shape[1],indexmax],dtype=np.float32)
                fullc=np.empty([2,s['mtemperatures'].shape[1],indexmax],dtype=np.float32)
                date_list = [base + datetime.timedelta(days=x) for x in range(0, indexmax)]
                date_listm = np.asarray([date_list[i].month for i in range(0, indexmax)])
                date_listy = np.asarray([date_list[i].year for i in range(0, indexmax)])
                il=1
            smask=np.abs(s['mtemperatures'])<400.
        except:
            print((st,' could not be read'))
            continue
        full.fill(np.nan)

        if np.sum(smask)>0:
            print((st, np.sum(smask), ' values found'))
        smask=np.abs(s['m'+fgdepname])>20.
        s['m'+fgdepname][smask]=np.nan
        s['mtemperatures'][smask]=np.nan
        s['mfg']=s['mtemperatures']+s['m'+fgdepname]
        full[:,:,s['mdatum']]=s['m'+fgdepname]
        #fullfg=np.empty([2,s['mtemperatures'].shape[1],indexmax],dtype=np.float32)
        #fullfg[:,:,s['mdatum']]=s['mfg']
        #fullt=np.empty([2,s['mtemperatures'].shape[1],indexmax],dtype=np.float32)
        #fullt[:,:,s['mdatum']]=s['mtemperatures']
        fullc[:]=full[:]
        allindices=[s['mdatum'][0]]
        try:
            for idx in s['index']:
                if idx>s['mdatum'][0] and idx!=20819:
                    allindices.append(idx)
            if s['index'][-1]<s['mdatum'][-1]:
                allindices+=[s['mdatum'][-1]]
            else:
                allindices[-1]=s['mdatum'][-1]
        except KeyError:
            print((st,'no break found'))
            allindices.append(s['mdatum'][-1])


    #    binsfg,goods,chunk,chunko=calc_bin_and_chunk(s['m'+fgdepname],fullfg,allindices)
    #    binst,goods,chunk,chunko=calc_bin_and_chunk(s['m'+fgdepname],fullt,allindices)

        bins=np.empty((len(allindices),12,s['m'+fgdepname].shape[0],s['m'+fgdepname].shape[1]))
        svar=np.empty((len(allindices),bins.shape[2],bins.shape[3]),dtype='bool')
        sratio=np.empty((len(allindices),bins.shape[2],bins.shape[3]),dtype='float32')
        scorr=np.empty((len(allindices),bins.shape[2],bins.shape[3]),dtype='float32')
        svar[:]=False
        weights=np.asarray((0.9,0.9,0.9,0.9,0.9,0.9,0.8,0.8,0.7,0.6,0.5,0.3,0.01,0.01,0.0,0.0))
        for i in range(1,len(allindices)):

            bins,goods,chunk,chunko=calc_bin_and_chunk(full,allindices,i,bins,weights)

            for k in range(chunk.shape[0]):
                for l in range(chunk.shape[1]):
                    if sum(goods[:,k,l]>5)>6:
                        c=chunk[k,l,:]-chunko[k,l,:]
                        mask=~np.isnan(c)
                        if sum(mask)>0:
                            selev=elev[0,k,allindices[i-1]:allindices[i]]
                            if mask.shape[0]>366:
                                sratio[i-1,k,l]=np.std(c[mask][:366])/np.std(chunk[k,l,mask][:366])
            #		    scorr[i-1,k,l]=np.corrcoef(c[mask],selev[mask])[0,1]
                                scorr[i-1,k,l]=scipy.stats.spearmanr(c[mask][:366],selev[mask][:366]).correlation
                            else:
                                sratio[i-1,k,l]=np.std(c[mask])/np.std(chunk[k,l,mask])
                                scorr[i-1,k,l]=scipy.stats.spearmanr(c[mask],selev[mask]).correlation
                            svar[i-1,k,l]=sratio[i-1,k,l]>0.1
                        else:
                            scorr[i-1,k,l]=np.nan
                            sratio[i-1,k,l]=np.nan
                    else:
                        scorr[i-1,k,l]=np.nan
                        sratio[i-1,k,l]=np.nan
                if sum(svar[i-1,k,:])>5:
                    svar[i-1,k,:]=True
                else:
                    svar[i-1,k,:]=False


            for k in range(chunk.shape[0]):
                for l in range(chunk.shape[1]-4):  # no RISE adjustment below 500 hPa
                    if sum(goods[:,k,l]>5)>6:
                        c=chunk[k,l,:]-chunko[k,l,:]
                        if svar[i-1,k,l]:
                            fullc[k,l,allindices[i-1]:allindices[i]]=chunk[k,l,:]-np.nanmean(chunk[k,l,:]-chunko[k,l,:])
    #		    print 'c/c',i,k,l,svar[k,l],np.nanstd(c)/np.nanstd(chunk[k,l,:])
                        if np.nanmax(np.abs(chunk[k,l,:]-chunko[k,l,:]))>5:
                            print((k,l,allindices[i-1]/365.25+1900,np.nanmax(np.abs(chunko[k,l,:]-chunk[k,l,:])),np.nanmean(chunk[k,l,:]-chunko[k,l,:])))
    #		else:
    #		    print k,l,sum(goods[:,k,l]>5)

    #    s['m'+fgdepname]=fullc[:,:,s['mdatum']]-full[:,:,s['mdatum']]
        if False:
            plt.figure(figsize=(10,5))
            for ipar in (0,1):
                plt.subplot(1,2,ipar+1)
                n=scorr.shape[0]*1.0
                for l in range(scorr.shape[0]-1):

                    try:
                        print((min(scorr[l,ipar,:12]),max(scorr[l,ipar,:12]),np.nanmean(sratio[l,ipar,2:8]),np.nanmean(np.abs(bins[l,:,ipar,2:8]))))
                        if np.nanmean(sratio[l,ipar,2:8])>0.1 and np.nanmean(np.abs(bins[l,:,ipar,2:8]))>0.2:
                            plt.semilogy(scorr[l,ipar,:12],plevs[:12],label=str(date_listy[allindices[l]]+1),color=[1-l/n,l/n,0.])
                    except:
                        pass
                plt.ylim(1000,10)
                plt.xlim(-1.5,1.)
                dpl=(1000,700,500,200,100,50,20,10)
                sdpl=[str(dp) for dp in dpl]
                plt.yticks(dpl,sdpl)
                plt.legend(loc='best',fontsize=10)
                plt.title(st+', Correlation {:0>2}GMT'.format(ipar*12))

            plt.savefig(plpath+st+'_corrprof.eps')
            plt.close()
        s['elev']=elev[0,:,s['mdatum']]
        print(('t2:',time.time()-t2))
        #continue


        if st in mlist:
            pass
        else:
            print((st,s['mtemperatures'].shape))
            pass


    #    lat=float(st[1])
    #    lon=float(st[2])
    #    s={'lat':float(st[1]),'lon':float(st[2])}
        if 'mbias' not in list(s.keys()):
            print((st+' no RASE bias'))
            continue
        if 'mtemperatures' not in list(s.keys()):
            print((st+' not plotted'))
            print((s['mbias'].shape))
            continue
        if 'mRAOBCORE' not in list(s.keys()):
            s['mRAOBCORE']=s['mtemperatures']-s['mtemperatures']
            print((st+' mRAOBCORE not found, zero correction supplied'))
        if 'mrasocorr' not in list(s.keys()) or 'mrasocorr' not in list(s.keys()):
            s['mrasocorr']=s['mtemperatures']-s['mtemperatures']
            print((st+' mrasocorr not found, zero correction supplied, cycle'))
            continue
    #	continue

        #try:
            #x2008=daysbetween(19000101,'days since 2008-01-01 00:00:00')        
            #x=np.where(np.logical_or(np.abs(s['mtemperatures']-
                                                        #s['mrasocorr'])>np.abs(s['mbias'])+5,
                                            #np.logical_and(np.isnan(s['mrasocorr']),s['mdatum']<x2008)))
            #s['mtemperatures'][x]=np.nan
        #except KeyError:
            #pass
        #plt.plot(1900+stats[st]['mdatum']/365.25,stats[st]['mtemperatures'][0,5,:])
        #plt.plot(1900+stats[st]['mdatum']/365.25,stats[st]['m'+fgdepname][0,5,:])
        #plt.plot(1900+stats[st]['mdatum']/365.25,stats[st]['mRAOBCORE'][0,5,:])
        #plt.show()
        t=time.time()
        if False:
            try:
                checkmean=-s['mtemperatures']+s['mRAOBCORE']+s['mbias']
                check=np.nanmean(checkmean[:,1:,:],axis=2)
                checkstd=np.nanstd(check)
                if checkstd<0.1:
                    s['newbias']=checkmean+s['mrasocorr']
                    s['newrichbias']=checkmean+s['mrichrasocorr']
                    x=np.where(s['mdatum']>=x2008)
                    s['newbias'][:,:,x]=s['mrasocorr'][:,:,x]
                    s['newrichbias'][:,:,x]=s['mrichrasocorr'][:,:,x]
                else:
                    if(not np.isnan(checkstd)):
                        plt.plot(check[0,:],15-np.arange(15))
                        plt.plot(check[1,:],15-np.arange(15))
                        plt.title(st+' {:5.2f}'.format(checkstd))
                        plt.savefig(plpath+st+'_profile.eps')
                        plt.close()
                        print(('ERAI bias and RAOBCORE v1.3 not consistent',checkstd))
                    s['newbias']=s['mrasocorr']
                    s['newrichbias']=s['mrichrasocorr']
            except KeyError:
                print('Had to use old bias')
                s['newbias']=s['mbias']
                s['newrichbias']=s['mbias']
                pass
        else:
            try:
                s['newbias']=s['mrasocorr']
                s['newrichbias']=s['mrichrasocorr']
            except:
                print('no bias corrections in file, cycling')
                continue

        try:
            x1979=daysbetween(19000101,'days since 1979-01-01 00:00:00')
            idx=np.where(s['mdatum']>x1979)
            if len(idx[0])>0:
                s['newrichbias'][:,:,idx[0]]-=fullc[:,:,s['mdatum'][idx]]-full[:,:,s['mdatum'][idx]]
                s['newbias'][:,:,idx[0]]-=fullc[:,:,s['mdatum'][idx]]-full[:,:,s['mdatum'][idx]]
            pass
        except:
            s['newrichbias']=-(fullc[:,:,s['mdatum']]-full[:,:,s['mdatum']])
            s['newbias']=-(fullc[:,:,s['mdatum']]-full[:,:,s['mdatum']])
            print((st,'no RICH adjustments'))
            pass

        fig=False
        pindex=[5]
    #    pindex=[0,1,2,3,4,5,6,7,8,9,10,11,12,13]
        nm=np.nanmax(np.abs(s['newrichbias']))
        if nm>7:
            print((st,'nm:',nm))
            pindex=[0,1,2,3,4,5,6,7,8,9,10,11,12,13]
    #	fig=True

        if fig and s['mbias'].shape[2]>10:
            for i in range(2):
                if sum(~np.isnan(s['m'+fgdepname][i,5,:]))<1000:
                    continue
                for ip in pindex:
                    plt.figure(figsize=(10,8))
                    plt.subplot(2,1,1)
                    try:
                        htime=np.arange(s['mdatum'][0],s['mdatum'][-1]+1)
                        htemp=np.empty(htime.shape[0])
                        htemp.fill(np.nan)
                        htemp2=np.empty(htime.shape[0])
                        htemp2.fill(np.nan)
                        htempfg=np.empty(htime.shape[0])
                        htempfg.fill(np.nan)
                        htempfg[:]=htemp[:]
                        htempfg[s['mdatum']-s['mdatum'][0]]=-s['m'+fgdepname][i,ip,:]
                        index=thin2(htempfg,30)

        #		index=thin2(s['mtemperatures'][i,ip,:],10)
                        mdi=1900+(htime[index])/365.25
                        plt.plot(mdi,rmeanw(htempfg,60)[index],'k',label='bgdep {:4.2f}'.format(np.nanstd(htempfg)),lw=2)
        #		plt.plot(1900+(htime)/365.25,rmeanw(htempfg,60),'k',label='bgdep',lw=2)
        #		plt.plot(mdi,rmeanw(-s['m'+fgdepname][i,ip,:],60)[index],'k',label='bgdep',lw=2)
                        #plt.plot(mdi,s['mtemperatures'][i,ip,index]-
                                #s['mRAOBCORE'][i,ip,index],'b',label='v1.3')
                        htemp[s['mdatum']-s['mdatum'][0]]=s['mrasocorr'][i,ip,:]
                        htemp[np.isnan(htempfg)]=np.nan
    #		    plt.plot(mdi,htemp[index],'r',label='v1.5')
                        plt.ylabel('bias estimate')
                        ax2=plt.twinx()
                        plt.plot(mdi,elev[0,i,s['mdatum'][0]:s['mdatum'][-1]+1][index],'b',label='elev',lw=0.5)
                        plt.ylim(-30.,90.)
                    except KeyError:
                        pass
            #	s['mbias'][i,ip,np.logical_and(s['mdatum']<x2008,
            #	                                            np.isnan(s['mtemperatures'][i,ip,index]))]=np.nan
                    ns='SN'
                    ew='WE'
                    ptitle=(st+', {0:0>2}GMT, {1}hPa, {2:5.2f}'+ns[s['lat']>0]+',{3:6.2f}'+ew[s['lon']>0]).format(i*12,plevs[ip],s['lat'],s['lon'])

                    htemp[s['mdatum']-s['mdatum'][0]]=s['mbias'][i,ip,:]
    #		plt.plot(mdi,htemp[index],'g',label='ERAI {:4.2f}'.format(np.nanstd(htempfg-htemp)))
                    plt.legend(loc='best',fontsize=12)
                    plt.title(ptitle)
                    plt.ylabel('Solar Elevation')
                    plt.xlim(np.floor(1900.+s['mdatum'][0]/365.25),2015)
    #		plt.xlim(1979,2015)
                    plt.subplot(2,1,2)
                    plt.plot(mdi,rmeanw(htempfg,60)[index],'k',label='bgdep {:4.2f}'.format(np.nanstd(htempfg)),lw=2)
                    try:
    #		    htemp[s['mdatum']-s['mdatum'][0]]=s['newbias'][i,ip,:]
    #		    plt.plot(mdi,htemp[index],'y',label='ERA5')#,lw=2)
                        htemp[s['mdatum']-s['mdatum'][0]]=s['newrichbias'][i,ip,:]
                        htemp2[s['mdatum']-s['mdatum'][0]]=s['newbias'][i,ip,:]
                        plt.plot(mdi,htemp[index],'m',label='ERA5_RICH {:4.2f}'.format(np.nanstd(htempfg-htemp)))#,lw=1)
    #		    plt.plot(mdi,htemp2[index],'r',label='ERA5_RAOB {:4.2f}'.format(np.nanstd(htempfg-htemp2)))#,lw=1)
                        startdate=datetime.datetime(1900,1,1,0,0,0)
                        mydate=datetime.datetime(1998,10,0o1,0,0,0)
                        ic=0
                        for ix in s['mdatum']:
                            sdate=startdate+datetime.timedelta(int(ix))
                            if sdate==mydate:
        #			print s['newrichbias'][0,:,ic]
        #			print s['newrichbias'][1,:,ic]
                                pass
                            ic+=1
                    except KeyError:
                        pass
                    plt.legend(loc='best',fontsize=12)
                    plt.title(ptitle)
                    plt.ylabel('bias estimate')
                    plt.xlim(np.floor(1900.+s['mdatum'][0]/365.25),2015)
    #		plt.xlim(1979,2015)

                    plt.ylim(-3,3)
                    pname=plpath+'/ERA5.'+st+'.{}.{:0>4}.eps'.format(i*12,plevs[ip])
                    plt.savefig(pname)
                    plt.close()
                    print((pname,time.time()-t))
                    pass


    #    if 'mtemperatures' in stats[st].keys():
    #        stats[st]['jra55_antemperatures']=stats[st]['mtemperatures']-np.nan
    #        stats[st]['jra55_fgtemperatures']=stats[st]['mtemperatures']-np.nan
    #        print st

        if 'mtemperatures' in list(s.keys()):
            flag=False
            fn='/home/srvx7/leo/fastscratch/rise/1.0/'+exp+'/'+st+'/feedbackmerged'+st+'.nc'
            f = netCDF4.Dataset(fn,"r")
            fno='/home/srvx7/leo/fastscratch/rise/1.0/'+exp+'/'+st+'/ERA5bc_RAOBCORE_v1.74_'+st+'.nc'
            fo = netCDF4.Dataset(fno,"w", format='NETCDF4_CLASSIC')

            for i in f.ncattrs():
                if i=='history':
                    setattr(fo,i,datetime.date.today().strftime("%Y/%m/%d"))
                elif i=='source':
                    setattr(fo,i,'RAOBCORE/RICH v1.7.4 + solar elevation dependency (from 197901 onward)' )
                elif i=='title':
                    setattr(fo,i,'Station daily temperature series with JRA55/CERA20C/ERApreSAT background departure statistics and RISE bias estimates' )
                else:
                    setattr(fo,i,getattr(f,i))
            for i in list(f.dimensions.keys()):
                if i=='time':
                    fo.createDimension(i,s["mdatum"].shape[0])
                else:
                    try:
                        fo.createDimension(i,len(f.dimensions[i]))
                    except:
                        flag=True
                        continue
            nalias=8    
            fo.createDimension('nalias',8)
            fo.createDimension('nchar',8)
            if flag:
                continue
            #nogos=['flags',u's_type', u'eijra_fgdep', u'jra55_fgdep', u'jra55_andep', u'e20c_andep', u'n20c_andep', u'ce20c_andep', u'erapresat_andep']
            tobecopied=['datum','hours','lat','lon','alt','press','temperatures','an_dep','fg_dep','source','mergedstats']
            for i in list(f.variables.keys()):
                var=f.variables[i]
                if i=='datum':
                    fo.createVariable(i,var.dtype,var.dimensions)
                    fo.variables[i][:]=s["mdatum"][:]
                elif i=='hours':
                    fo.createVariable(i,var.dtype,var.dimensions)
                    fo.variables[i][:]=s["mhours"][:]
                elif i=='temperatures':
                    fo.createVariable(i,var.dtype,var.dimensions)
                    s["mtemperatures"][np.isnan(s["mtemperatures"])]=-999.
                    fo.variables[i][:]=s["mtemperatures"][:]
                elif i=='fg_dep':
                    fo.createVariable(i,var.dtype,var.dimensions)
                    s["m"+fgdepname][np.isnan(s["m"+fgdepname])]=-999.
                    fo.variables[i][:]=s["m"+fgdepname][:]
                elif i=='an_dep':
                    s["newbias"][np.isnan(s["newbias"])]=999.
                    s["newrichbias"][np.isnan(s["newrichbias"])]=999.
                    try:
                        fo.createVariable('bias',var.dtype,var.dimensions)
                    except:
                        pass
                    fo.variables['bias'][:]=-s["newbias"][:]
                    fo.createVariable('richbias',var.dtype,var.dimensions)
                    fo.variables['richbias'][:]=-s["newrichbias"][:]
                elif i in ('lon','lat','alt','press'):
                    fo.createVariable(i,var.dtype,var.dimensions)
                    fo.variables[i][:]=var[:]
                elif i in ('source'):
#		    str_out = netCDF4.stringtochar(np.array(['test'], 'S4'))
                    fo.createVariable(i,'S1',('time','nchar'))
                    x=np.empty((var.shape[0]),dtype='S8')
                    x.fill('        ')
                    try:
                        svar=var[:]
                    except:
                        print('could not read source, supplying BUFRDATA')
                        svar=x
                    str_out = netCDF4.stringtochar(np.asarray(svar,'S8'))
                    fo.variables[i][:]=str_out
                elif i in ('mergedstats'):
                    fo.createVariable(i,'S1',('time','nalias','nchar'))
                    tt=time.time()
                    x=np.empty((var.shape[0],nalias),dtype='S8')
                    x.fill('        ')
                    try:
                        svar=var[:]
                        for k in range(svar.shape[0]):
                            l=svar[k].split(',')
                            #if len(l)>1:
                                #print 'l>1'
                            for m in range(len(l)):
                                x[k,m]=l[m]+' '*(8-len(l[m]))
                    except:
                        print('could not read mergedstats, filling with WMO number')
                        x.fill(st[1:]+'   ')

                    str_out = netCDF4.stringtochar(x)
                    fo.variables[i][:]=str_out
                    print(('mergevar:',time.time()-tt))
                else:
                    if i in tobecopied:
                        print((i,'some unknown variable'))
                        fo.createVariable(i,var.dtype,var.dimensions)

                for j in var.ncattrs():
                    if j!='_FillValue' and j!='scale_factor' and j!='add_offset':
                        if i in tobecopied:
                            if i=='an_dep':
                                setattr(fo.variables['bias'],j,getattr(var,j))
                                setattr(fo.variables['richbias'],j,getattr(var,j))
                            elif i=='datum' and j=='units':
                                setattr(fo.variables[i],j,'days since 1900-01-01 00:00:00')
                            else:
                                setattr(fo.variables[i],j,getattr(var,j))

            setattr(fo.variables['fg_dep'],'infos','obs-CERA20C/ERApreSAT ensemble average up to 196112; obs-ERA_Interim 1962 onwards')
            setattr(fo.variables['bias'],'long_name','RAOBCORE v1.6 RISE bias estimate')
            setattr(fo.variables['richbias'],'long_name','RICH v1.6 RICH bias estimate')
            setattr(fo.variables['source'],'info','Preferred source used for merged temperature record')
            setattr(fo.variables['source'],'valid_entries','BUFRDATA: (ECMWF holdings), NCARUA(20,21,24): sources from NCAR, CHUAN2.1: ERA-CLIM(2) digitized data')
            setattr(fo.variables['mergedstats'],'info','ODB statIDs that matched during merge - bias adjustments should be applied to those')


            fo.close()
            f.close()

    print((time.time()-t))
