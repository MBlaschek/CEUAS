#!/usr/bin/env python
import traceback
import sys,glob
import os.path

import numpy
from datetime import date
import netCDF4
import time
from numba import njit
from rasotools.utils import *
from rasotools.anomaly import *
import matplotlib.pylab as plt
import scipy.stats
import f90nml
sys.path.append('../../cds-backend/code/')
import cds_eua3 as eua
data=eua.vm_request_wrapper({'variable': 'temperature', 'statid': '01001'}, overwrite=True)

plt.rcParams['lines.linewidth'] = 3

def read_temperature(stats,fn,tidx,varlist=['temperatures','fg_dep','an_dep']):

    if 'mdatum' not in list(stats.keys()):
        try:
            f=netCDF4.Dataset(fn,'r')
            f.set_auto_mask(False)
            s=getattr(f.variables['datum'],'units').split()[2].split('-')[0]
            shift=tidx[(int(s)-1900)*12]
            stats['mdatum']=f.variables['datum'][0]+shift
            stats['mhours']=f.variables['hours'][:]
            stats['lat']=f.variables['lat'][0]
            stats['lon']=f.variables['lon'][0]
            for var in varlist:
                try:
                    miss_val=getattr(f.variables['temperatures'],'missing_value')
                    vals=f.variables[var][:]
                    vals[vals==miss_val]=numpy.nan
                    if var=='jracepre_fgdep':
                        vals2=f.variables['fg_dep'][:]
                        vals2[vals2==miss_val]=numpy.nan
                        idx=numpy.where(f.variables['datum'][0]>42735)[0]
                        if idx.shape[0]>0:
                            vals[:,:,idx]=-vals2[:,:,idx]
                except KeyError:
                    vals=f.variables['temperatures'][:]+numpy.nan
                if var not in ['fg_dep','an_dep']:
                    stats['m'+var]=vals
                else:
                    stats['m'+var]=-vals
            if 'merged' not in fn:
## departures in merged files are defined as bg-obs, not obs-bg, thus the sign change
                b=f.variables['bias'][:]
            else:
                b=stats['mtemperatures']-stats['mtemperatures']

            b[numpy.isnan(stats['mtemperatures'])]=numpy.nan   
            stats['mbias']=-b
            stats['mtemperatures']-=b
            stats['m'+fgdepname]=b-stats['m'+fgdepname]
            stats['man_dep']=b-stats['man_dep']

            f.close()
        except:
            p=False
            pass
    else:

#   merge
        ntimes=list()
        try:
            f=netCDF4.Dataset(fn,'r')
            f.set_auto_mask(False)
            if 'datum' in list(f.variables.keys()):
                s=getattr(f.variables['datum'],'units').split()[2].split('-')[0]
                shift=tidx[(int(s)-1900)*12]
                times=f.variables['datum'][0]+shift
                ntimes=numpy.zeros(indexmax,dtype=numpy.int)
#                ntimes[times]=1
                ntimes[stats['mdatum']]=1
                nhours=numpy.zeros([2,indexmax],dtype=numpy.int)-999
                try: 
                    nhours[:,times]=f.variables['hours'][:]
                except KeyError:
                    pass
                nhours[:,stats['mdatum']]=stats['mhours'][:]
                t=time.time()
                mdathilf=numpy.asarray(numpy.where(ntimes>0)[0],dtype=numpy.int)

                if 'mtemperatures' not in list(stats.keys()):
                    return
                for var in varlist:
                    mergevar(stats,fn,f.variables,times,mdathilf,var)

                #print 'merge',time.time()-t
                stats['mdatum']=mdathilf
                stats['mhours']=nhours[:,mdathilf]
            f.close()
        except:
            p=False
            pass

def mergevar(stats,fn,variables,times,mdathilf,var):
    temperatures=numpy.empty([2,stats['mtemperatures'].shape[1],indexmax],dtype=numpy.float32)
    miss_val=getattr(variables[var],'missing_value')
    vals=variables[var][:]

    if 'corr' in fn:
        timeoffset=daysbetween(19000101,variables["datum"].units)        
        index=variables['datum'][0,:]-1+timeoffset
        stats['index']=index
#                index=numpy.arange(len(index),dtype=int)   
        try:
            nc_miss_val=variables[var].getncattr('_FillValue')
        except:
            pass
        try:
            nc_miss_val=variables[var].getncattr('missing_value')
        except:
            pass


        temperatures[:]=0.
        hilf=variables[var][:]
        hilf[hilf==nc_miss_val]=numpy.nan
        expand2(hilf,index,numpy.arange(variables[var].shape[1],dtype='int'),temperatures)
#                copystride(temperatures,variables[var][:],index,0,0,
#                           numpy.arange(variables[var].shape[1],dtype='int'),nc_miss_val)
        vals=temperatures.copy()

    else:
        temperatures[:]=numpy.nan
        vals[vals==miss_val]=numpy.nan
        if 'merged' not in fn:
            try:
                b=variables['bias'][:]
            except KeyError:
                b=vals-vals
        else:
            b=vals-vals
    # departures in merged files are defined as bg-obs, not obs-bg, thus the sign change
        if var=='temperatures':
            vals-=b
        elif var=='fg_dep':
            vals=b-vals
        elif var=='an_dep':
            vals=b-vals

        temperatures[:,:,times]=vals
    if 'rio' in fn:
        stats['mrich'+var]=temperatures[:,:,mdathilf]
    else:
        stats['m'+var]=temperatures[:,:,mdathilf]


def calc_bin_and_chunk(full,allindices,i,bins,weights):

    chunko=full[:,:,allindices[i-1]:allindices[i]]
    chunk=numpy.empty(chunko.shape)
    chunk[:]=chunko[:]
    goods=numpy.empty((12,chunko.shape[0],chunko.shape[1]),dtype='int')
    for j in range(12):
        mask=date_listm[allindices[i-1]:allindices[i]]==(j+1)
        mask=numpy.where(date_listm[allindices[i-1]:allindices[i]]==(j+1))[0]
        bins[i-1,j,:,:]=numpy.nanmean(chunk[:,:,mask],axis=2)
        goods[j,:]=numpy.sum(~numpy.isnan(chunk[:,:,mask]),axis=2)
        bmask=goods[j,:]<=5
        bins[i-1,j,bmask]=numpy.nan
        for k in range(chunk.shape[0]):
            for l in range(chunk.shape[1]):
                if goods[j,k,l]>5:
                    chunk[k,l,mask]-=weights[l]*bins[i-1,j,k,l]
                else:
#			chunk[k,l,mask]=numpy.nan
                    chunko[k,l,mask]=numpy.nan

    return bins,goods,chunk,chunko


if __name__ == "__main__":


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
    plevs=numpy.asarray([10,20,30,50,70,100,150,200,250,300,400,500,700,850,925,1000])
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
        elev=calc_elevangles(0.,numpy.asarray([s['lat']]),numpy.asarray([s['lon']]))
        for y in range(120):
            idx=datetime.date(1900+y,1,1)-base
            idx=idx.days
            elev[0,:,idx:idx+366]=elev[0,:,:366]
        print(('t2:',time.time()-t2))

        bdatum=[]
        try:
            if il==0:
                full=numpy.empty([2,s['mtemperatures'].shape[1],indexmax],dtype=numpy.float32)
                fullc=numpy.empty([2,s['mtemperatures'].shape[1],indexmax],dtype=numpy.float32)
                date_list = [base + datetime.timedelta(days=x) for x in range(0, indexmax)]
                date_listm = numpy.asarray([date_list[i].month for i in range(0, indexmax)])
                date_listy = numpy.asarray([date_list[i].year for i in range(0, indexmax)])
                il=1
            smask=numpy.abs(s['mtemperatures'])<400.
        except:
            print((st,' could not be read'))
            continue
        full.fill(numpy.nan)

        if numpy.sum(smask)>0:
            print((st, numpy.sum(smask), ' values found'))
        smask=numpy.abs(s['m'+fgdepname])>20.
        s['m'+fgdepname][smask]=numpy.nan
        s['mtemperatures'][smask]=numpy.nan
        s['mfg']=s['mtemperatures']+s['m'+fgdepname]
        full[:,:,s['mdatum']]=s['m'+fgdepname]
        #fullfg=numpy.empty([2,s['mtemperatures'].shape[1],indexmax],dtype=numpy.float32)
        #fullfg[:,:,s['mdatum']]=s['mfg']
        #fullt=numpy.empty([2,s['mtemperatures'].shape[1],indexmax],dtype=numpy.float32)
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

        bins=numpy.empty((len(allindices),12,s['m'+fgdepname].shape[0],s['m'+fgdepname].shape[1]))
        svar=numpy.empty((len(allindices),bins.shape[2],bins.shape[3]),dtype='bool')
        sratio=numpy.empty((len(allindices),bins.shape[2],bins.shape[3]),dtype='float32')
        scorr=numpy.empty((len(allindices),bins.shape[2],bins.shape[3]),dtype='float32')
        svar[:]=False
        weights=numpy.asarray((0.9,0.9,0.9,0.9,0.9,0.9,0.8,0.8,0.7,0.6,0.5,0.3,0.01,0.01,0.0,0.0))
        for i in range(1,len(allindices)):

            bins,goods,chunk,chunko=calc_bin_and_chunk(full,allindices,i,bins,weights)

            for k in range(chunk.shape[0]):
                for l in range(chunk.shape[1]):
                    if sum(goods[:,k,l]>5)>6:
                        c=chunk[k,l,:]-chunko[k,l,:]
                        mask=~numpy.isnan(c)
                        if sum(mask)>0:
                            selev=elev[0,k,allindices[i-1]:allindices[i]]
                            if mask.shape[0]>366:
                                sratio[i-1,k,l]=numpy.std(c[mask][:366])/numpy.std(chunk[k,l,mask][:366])
            #		    scorr[i-1,k,l]=numpy.corrcoef(c[mask],selev[mask])[0,1]
                                scorr[i-1,k,l]=scipy.stats.spearmanr(c[mask][:366],selev[mask][:366]).correlation
                            else:
                                sratio[i-1,k,l]=numpy.std(c[mask])/numpy.std(chunk[k,l,mask])
                                scorr[i-1,k,l]=scipy.stats.spearmanr(c[mask],selev[mask]).correlation
                            svar[i-1,k,l]=sratio[i-1,k,l]>0.1
                        else:
                            scorr[i-1,k,l]=numpy.nan
                            sratio[i-1,k,l]=numpy.nan
                    else:
                        scorr[i-1,k,l]=numpy.nan
                        sratio[i-1,k,l]=numpy.nan
                if sum(svar[i-1,k,:])>5:
                    svar[i-1,k,:]=True
                else:
                    svar[i-1,k,:]=False


            for k in range(chunk.shape[0]):
                for l in range(chunk.shape[1]-4):  # no RISE adjustment below 500 hPa
                    if sum(goods[:,k,l]>5)>6:
                        c=chunk[k,l,:]-chunko[k,l,:]
                        if svar[i-1,k,l]:
                            fullc[k,l,allindices[i-1]:allindices[i]]=chunk[k,l,:]-numpy.nanmean(chunk[k,l,:]-chunko[k,l,:])
    #		    print 'c/c',i,k,l,svar[k,l],numpy.nanstd(c)/numpy.nanstd(chunk[k,l,:])
                        if numpy.nanmax(numpy.abs(chunk[k,l,:]-chunko[k,l,:]))>5:
                            print((k,l,allindices[i-1]/365.25+1900,numpy.nanmax(numpy.abs(chunko[k,l,:]-chunk[k,l,:])),numpy.nanmean(chunk[k,l,:]-chunko[k,l,:])))
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
                        print((min(scorr[l,ipar,:12]),max(scorr[l,ipar,:12]),numpy.nanmean(sratio[l,ipar,2:8]),numpy.nanmean(numpy.abs(bins[l,:,ipar,2:8]))))
                        if numpy.nanmean(sratio[l,ipar,2:8])>0.1 and numpy.nanmean(numpy.abs(bins[l,:,ipar,2:8]))>0.2:
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
            #x=numpy.where(numpy.logical_or(numpy.abs(s['mtemperatures']-
                                                        #s['mrasocorr'])>numpy.abs(s['mbias'])+5,
                                            #numpy.logical_and(numpy.isnan(s['mrasocorr']),s['mdatum']<x2008)))
            #s['mtemperatures'][x]=numpy.nan
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
                check=numpy.nanmean(checkmean[:,1:,:],axis=2)
                checkstd=numpy.nanstd(check)
                if checkstd<0.1:
                    s['newbias']=checkmean+s['mrasocorr']
                    s['newrichbias']=checkmean+s['mrichrasocorr']
                    x=numpy.where(s['mdatum']>=x2008)
                    s['newbias'][:,:,x]=s['mrasocorr'][:,:,x]
                    s['newrichbias'][:,:,x]=s['mrichrasocorr'][:,:,x]
                else:
                    if(not numpy.isnan(checkstd)):
                        plt.plot(check[0,:],15-numpy.arange(15))
                        plt.plot(check[1,:],15-numpy.arange(15))
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
            idx=numpy.where(s['mdatum']>x1979)
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
        nm=numpy.nanmax(numpy.abs(s['newrichbias']))
        if nm>7:
            print((st,'nm:',nm))
            pindex=[0,1,2,3,4,5,6,7,8,9,10,11,12,13]
    #	fig=True

        if fig and s['mbias'].shape[2]>10:
            for i in range(2):
                if sum(~numpy.isnan(s['m'+fgdepname][i,5,:]))<1000:
                    continue
                for ip in pindex:
                    plt.figure(figsize=(10,8))
                    plt.subplot(2,1,1)
                    try:
                        htime=numpy.arange(s['mdatum'][0],s['mdatum'][-1]+1)
                        htemp=numpy.empty(htime.shape[0])
                        htemp.fill(numpy.nan)
                        htemp2=numpy.empty(htime.shape[0])
                        htemp2.fill(numpy.nan)
                        htempfg=numpy.empty(htime.shape[0])
                        htempfg.fill(numpy.nan)
                        htempfg[:]=htemp[:]
                        htempfg[s['mdatum']-s['mdatum'][0]]=-s['m'+fgdepname][i,ip,:]
                        index=thin2(htempfg,30)

        #		index=thin2(s['mtemperatures'][i,ip,:],10)
                        mdi=1900+(htime[index])/365.25
                        plt.plot(mdi,rmeanw(htempfg,60)[index],'k',label='bgdep {:4.2f}'.format(numpy.nanstd(htempfg)),lw=2)
        #		plt.plot(1900+(htime)/365.25,rmeanw(htempfg,60),'k',label='bgdep',lw=2)
        #		plt.plot(mdi,rmeanw(-s['m'+fgdepname][i,ip,:],60)[index],'k',label='bgdep',lw=2)
                        #plt.plot(mdi,s['mtemperatures'][i,ip,index]-
                                #s['mRAOBCORE'][i,ip,index],'b',label='v1.3')
                        htemp[s['mdatum']-s['mdatum'][0]]=s['mrasocorr'][i,ip,:]
                        htemp[numpy.isnan(htempfg)]=numpy.nan
    #		    plt.plot(mdi,htemp[index],'r',label='v1.5')
                        plt.ylabel('bias estimate')
                        ax2=plt.twinx()
                        plt.plot(mdi,elev[0,i,s['mdatum'][0]:s['mdatum'][-1]+1][index],'b',label='elev',lw=0.5)
                        plt.ylim(-30.,90.)
                    except KeyError:
                        pass
            #	s['mbias'][i,ip,numpy.logical_and(s['mdatum']<x2008,
            #	                                            numpy.isnan(s['mtemperatures'][i,ip,index]))]=numpy.nan
                    ns='SN'
                    ew='WE'
                    ptitle=(st+', {0:0>2}GMT, {1}hPa, {2:5.2f}'+ns[s['lat']>0]+',{3:6.2f}'+ew[s['lon']>0]).format(i*12,plevs[ip],s['lat'],s['lon'])

                    htemp[s['mdatum']-s['mdatum'][0]]=s['mbias'][i,ip,:]
    #		plt.plot(mdi,htemp[index],'g',label='ERAI {:4.2f}'.format(numpy.nanstd(htempfg-htemp)))
                    plt.legend(loc='best',fontsize=12)
                    plt.title(ptitle)
                    plt.ylabel('Solar Elevation')
                    plt.xlim(numpy.floor(1900.+s['mdatum'][0]/365.25),2015)
    #		plt.xlim(1979,2015)
                    plt.subplot(2,1,2)
                    plt.plot(mdi,rmeanw(htempfg,60)[index],'k',label='bgdep {:4.2f}'.format(numpy.nanstd(htempfg)),lw=2)
                    try:
    #		    htemp[s['mdatum']-s['mdatum'][0]]=s['newbias'][i,ip,:]
    #		    plt.plot(mdi,htemp[index],'y',label='ERA5')#,lw=2)
                        htemp[s['mdatum']-s['mdatum'][0]]=s['newrichbias'][i,ip,:]
                        htemp2[s['mdatum']-s['mdatum'][0]]=s['newbias'][i,ip,:]
                        plt.plot(mdi,htemp[index],'m',label='ERA5_RICH {:4.2f}'.format(numpy.nanstd(htempfg-htemp)))#,lw=1)
    #		    plt.plot(mdi,htemp2[index],'r',label='ERA5_RAOB {:4.2f}'.format(numpy.nanstd(htempfg-htemp2)))#,lw=1)
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
                    plt.xlim(numpy.floor(1900.+s['mdatum'][0]/365.25),2015)
    #		plt.xlim(1979,2015)

                    plt.ylim(-3,3)
                    pname=plpath+'/ERA5.'+st+'.{}.{:0>4}.eps'.format(i*12,plevs[ip])
                    plt.savefig(pname)
                    plt.close()
                    print((pname,time.time()-t))
                    pass


    #    if 'mtemperatures' in stats[st].keys():
    #        stats[st]['jra55_antemperatures']=stats[st]['mtemperatures']-numpy.nan
    #        stats[st]['jra55_fgtemperatures']=stats[st]['mtemperatures']-numpy.nan
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
                    s["mtemperatures"][numpy.isnan(s["mtemperatures"])]=-999.
                    fo.variables[i][:]=s["mtemperatures"][:]
                elif i=='fg_dep':
                    fo.createVariable(i,var.dtype,var.dimensions)
                    s["m"+fgdepname][numpy.isnan(s["m"+fgdepname])]=-999.
                    fo.variables[i][:]=s["m"+fgdepname][:]
                elif i=='an_dep':
                    s["newbias"][numpy.isnan(s["newbias"])]=999.
                    s["newrichbias"][numpy.isnan(s["newrichbias"])]=999.
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
#		    str_out = netCDF4.stringtochar(numpy.array(['test'], 'S4'))
                    fo.createVariable(i,'S1',('time','nchar'))
                    x=numpy.empty((var.shape[0]),dtype='S8')
                    x.fill('        ')
                    try:
                        svar=var[:]
                    except:
                        print('could not read source, supplying BUFRDATA')
                        svar=x
                    str_out = netCDF4.stringtochar(numpy.asarray(svar,'S8'))
                    fo.variables[i][:]=str_out
                elif i in ('mergedstats'):
                    fo.createVariable(i,'S1',('time','nalias','nchar'))
                    tt=time.time()
                    x=numpy.empty((var.shape[0],nalias),dtype='S8')
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
