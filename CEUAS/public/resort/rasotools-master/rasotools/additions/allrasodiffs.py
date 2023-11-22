import os,glob,sys
import numpy
#from scipy.io import netcdf
from scipy.stats import linregress
from Magics.macro import *
import netCDF4
import time
from rasotools.additions.utils import *
#from rasotools.pythranutils import *
from rasotools.additions.anomaly import *
#from rasotools.panomaly import *
from rasotools.additions.set_trendmaps import *
#from mpl_toolkits.basemap import Basemap
#import matplotlib.pylab as plt
#import matplotlib.patches as mpatches
#import matplotlib.cm as cm
#import matplotlib.colors as colors
from rasotools.additions.dumpandspagarr import *
from rasotools.additions.define_datasets import *

from netCDF4 import Dataset
from numba import njit,prange
from rasotools.additions.allrasotrends import *


def allrasodiffs(path,tasks,plotproperties,intervals,days,stnames,interv,gstations=None,sats=0,init=False,daynight=None,ref="tmcorr"):

    pindex=numpy.asarray(plotproperties['pindex'],dtype='int')
    msupindex=numpy.asarray(plotproperties['msupindex'],dtype='int')
#    msups=plotproperties['msups']
    msups=['TLT','TMT','TTS','TLS']
    satlist=['rss','uah','star','wegc','merra2','20CRv3']
    reanlist=['wegc','merra2','20CRv3']

    os.chdir(path)
    startyear=tasks[0]["startdate"]//10000
    slowload=len(plotproperties['plotlist'])!=1 or plotproperties['plotlist'][0]!='stations'
    first=(intervals[0]-startyear)*12
    last=(intervals[1]-startyear+1)*12
    tolerance=intervals[2]
    t=time.time()
    if slowload:
        if interv==0 and 'satseries' in plotproperties['plotlist']:
            sats={'uah':dict(),'rss':dict()}
            if path[-1]=='/':
                spath='/'.join(path.split('/')[:-2])
            else:
                spath='/'.join(path.split('/')[:-1])

            sats['uah']['full'],sats['uah']['18']=read_uah(spath+'/../MSUUAHDaten/','5.6',startyear=startyear,endyear=2015,ens=0)
            sats['rss']['full'],sats['rss']['18']=read_rss(spath+'/../MSUDaten/','V3_3',startyear=startyear,endyear=2015,ens=0)

        hadtime=numpy.arange(last-first,dtype=float)
        hilf=hadtime-hadtime

    istat,lats,lons,ps,stnames=loadsondes(path,tasks,days,plotproperties,stnames,init,slowload,satlist)

    shortnames=[t['shortname'] for t in tasks]
    if 'suny' in shortnames or 'sunyhom' in shortnames:
        fill_suny(tasks,shortnames,stnames,lats,lons,ps,plotproperties['pindex'])

    if isinstance(stnames,list):
        stnames=numpy.asarray(stnames)
    #istnames=stnames.astype(numpy.int)
    tmshape=tasks[0]["data"].shape
    climatology=numpy.zeros(12,numpy.float64)

    slopes=[]
    sli=-1
    for k in range(len(tasks)):
        #if tasks[k]['shortname'] in sats.keys():
            #tasks[k]['msudata']=satst(lats,lons,sats[tasks[k]['shortname']]['full'])

        for inter in range(intervals.size//3):
            sli+=1
            slopes.append({ "shortname":tasks[k]["shortname"],                            "ens":tasks[k]["ens"],
                            "interval":intervals[:],
                            "data":numpy.empty([tmshape[0],len(tasks[k]["ens"]),tmshape[2]+1,
                                                pindex.shape[0]+msupindex.shape[0]],numpy.float64),
                            }
                          )
            slopes[sli]["data"].fill(numpy.nan)


    dists=numpy.zeros((tmshape[0]+1)*tmshape[0]//2,numpy.float64)
    tdist(dists,lats,lons,1) 


    ni=36
    nj=18
    glons=5.+10.*numpy.arange(ni)
    glats=-85.+10.*numpy.arange(nj)
    belts=numpy.asarray([[0,18],[11,18],[0,8],[7,11],[9,18],[0,9]])
    beltnames=["Globe",'NHEx','SHEx','Tropics','NH','SH']

    gstatindex=numpy.zeros([nj,ni,100],dtype=numpy.int32) # 20 is max stations per 10x10 degree box
    find_gstatindex(glons,glats,lons,lats,gstatindex)

    currentdatabg=numpy.zeros([tmshape[0],tmshape[2],pindex.shape[0],tmshape[4]],numpy.float32)
    currentdatatm=numpy.zeros([tmshape[0],tmshape[2],pindex.shape[0],tmshape[4]],numpy.float32)
    currentdata  =numpy.zeros([tmshape[0],tmshape[2],pindex.shape[0],tmshape[4]],numpy.float32)
    zero=         numpy.zeros([tmshape[0],tmshape[2],pindex.shape[0],tmshape[4]],numpy.float32)
    currentdatauah=numpy.zeros([tmshape[0],tmshape[2],msupindex.shape[0],tmshape[4]],numpy.float32)
    currentdatarss=numpy.zeros([tmshape[0],tmshape[2],msupindex.shape[0],tmshape[4]],numpy.float32)
    currentdatamsu=numpy.zeros([tmshape[0],tmshape[2],msupindex.shape[0],tmshape[4]],numpy.float32)
    jcurrentdata=numpy.zeros((1,1,1,1))

    orig=numpy.zeros(tmshape[4],dtype=numpy.float32)
    anomaly=numpy.zeros(tmshape[4],dtype=numpy.float32)
    gst=numpy.zeros(tmshape[4],dtype=numpy.float32)

    andeplist,andepstdlist,andeprmslist=make_andeplist()

    keytmcorr=-1
    keytm=-1
    for key in range(len(tasks)):
        if tasks[key]["shortname"] in ['uah','rss','star']:
            tasks[key]['cost']=numpy.empty((tasks[key]['msudata'].shape[1],3,pindex.shape[0]+msupindex.shape[0]))
        else:
            tasks[key]['cost']=numpy.empty((tasks[key]['data'].shape[1],3,pindex.shape[0]+msupindex.shape[0]))
        tasks[key]['cost'].fill(numpy.nan)
        if tasks[key]["shortname"]=="tmcorr":
            keytmcorr=key
        if tasks[key]["shortname"]=="tm":
            keytm=key        
        if tasks[key]["shortname"]=="rss":
            keyrss=key        
        if tasks[key]["shortname"]=="uah":
            keyuah=key        

    try:
        expandandadd(tasks[keyuah]["msudata"],currentdatauah,tasks[keytm]["index"],msupindex,0,currentdatauah,0.0)
        expandandadd(tasks[keyrss]["msudata"],currentdatarss,tasks[keytm]["index"],msupindex,0,currentdatarss,0.0)
    except:
        pass

    sli=-1
    nodata=[]
    for key in range(len(tasks)):
        print((tasks[key]["shortname"]))
        if tasks[key]["shortname"] in ["rio","rit","ce20c_andep"]:
            enslist=tasks[key]["ens"]
        else:
            enslist=[0]

        nodata.append(False)    
        for iens in enslist: #range(len(tasks[key]["ens"])):

            # replace RAOBCORE adjustments with unadjusted series

            sat=False
            t=time.time()
            dkey='data'
            if tasks[key]['shortname'] in ('rss','uah','star'):
                dkey='msudata'
            if nodatafound(tasks[key][dkey]):
                nodata[-1]=True
                break
            currentdata[:]=numpy.nan
            if tasks[key]["shortname"] == "tmcorr":
                picopy(currentdata,numpy.zeros_like(currentdata),tasks[key]['data'],iens,pindex)
                #expandandadd(tasks[key]["data"],tasks[key]["data"].reshape(tmshape[0],tmshape[2],tmshape[3],tmshape[4]),
                                #tasks[key]["index"],pindex,iens-iens,currentdata,0.0)
            elif tasks[key]["shortname"] == "tm":
                expandandadd(tasks[key]["data"],tasks[keytmcorr]["data"].reshape(tmshape[0],tmshape[2],tmshape[3],tmshape[4]),
                             tasks[key]["index"],pindex,iens-iens,currentdatatm,+1.0)
                currentdata[:]=currentdatatm[:]
            elif ('rio' in tasks[key]["shortname"] or 'rit' in tasks[key]["shortname"]) and 'corr' not in tasks[key]["shortname"] :  
                expandandadd(tasks[key]["data"],currentdatatm,tasks[key]["index"],pindex,iens,currentdata,-1.0)
            elif 'riocorr' in tasks[key]["shortname"] or 'ritcorr' in tasks[key]["shortname"]  :  
                expandandadd(tasks[key]["data"],currentdatatm-currentdatatm,tasks[key]["index"],pindex,iens,currentdata,-1.0)
            elif tasks[key]["shortname"] in ("rcorr",'riocorr','ritcorr'):
                expandandadd(tasks[key]["data"],currentdatatm-currentdatatm,tasks[key]["index"],pindex,iens,currentdata,-1.0)
            elif tasks[key]["shortname"] in ["eraibc","era5v2","era5v2429","era5v2rich","era5v3","era5v4","era5v5","era5v7","era5v7rich"]:
#                expandandadd(tasks[key]["data"],currentdatatm,tasks[key]["index"],pindex,iens,currentdata,-1.0)
                fak=1.0
                if tasks[key]["shortname"] in ["era5v2429"]:
                    fak=-1.0
                if ref=='zero':
                    #try:
                    #picopy(currentdata,currentdatatm-currentdatatm-currentdata2429,fak*tasks[key]["data"],iens,pindex)
                    picopy(currentdata,currentdatatm-currentdatatm,fak*tasks[key]["data"],iens,pindex)
                    #except:
                        #picopy(currentdata,currentdatatm-currentdatatm,fak*tasks[key]["data"],iens,pindex)
                        #currentdata2429=currentdata.copy()
                
                else:
                    picopy(currentdata,currentdatatm,fak*tasks[key]["data"],iens,pindex)

            elif tasks[key]["shortname"] in ["bg"]: 
                expandandadd(tasks[key]["data"],currentdata,tasks[keytm]["index"],pindex,iens-iens,currentdatabg,0.0)
                currentdata[:]=currentdatabg.copy()
            elif tasks[key]["shortname"] in andepstdlist or tasks[key]["shortname"] in andeprmslist or tasks[key]["shortname"]=='hadat': 
                expandandadd(tasks[key]["data"],zero,tasks[keytm]["index"],pindex,iens-iens,currentdata,0.0)
#                currentdatabg=tasks[key]["data"].reshape(tmshape[0],tmshape[2],ps.shape[0],tmshape[4])[:,:,pindex,:]
#                currentdata[:]=currentdatabg.copy()
            elif tasks[key]["shortname"] in ["comp-bg","comp2-bg","comp4-bg",]: 
                picopy(currentdata,currentdatabg,tasks[key]["data"],iens,pindex)
#                expandandadd(tasks[key]["data"],currentdatabg,tasks[keytm]["index"],pindex,0,currentdata,-1.0)
            elif tasks[key]["shortname"] in ["erai_fggpsdep","erai_fggpswetdep"]: 
                picopy(currentdata,currentdatabg,-tasks[key]["data"],iens,pindex)
#                expandandadd(tasks[key]["data"],currentdatabg,tasks[keytm]["index"],pindex,0,currentdata,-1.0)
            elif tasks[key]["shortname"] in andeplist or tasks[key]["shortname"] in ['ce20c_andep']: 
                if keytm>key:
                    expandandadd(tasks[keytm]["data"],tasks[keytmcorr]["data"].reshape(tmshape[0],tmshape[2],tmshape[3],tmshape[4]),
                                 tasks[keytm]["index"],pindex,iens,currentdatatm,+1.0)
                if ref=='zero':    
                    picopy(currentdata,zero,tasks[key]["data"],iens,pindex)
                else:
                    picopy(currentdata,currentdatatm,tasks[key]["data"],iens,pindex)
#		    currentdata[:]=currentdatatm[:]#+tasks[key]["data"].reshape(tmshape[0],tmshape[2],ps.shape[0],tmshape[4])[:,:,pindex,:]
            elif tasks[key]["shortname"] in ["ce20c0_andep"]:
#                currentdata[:]=currentdatatm+tasks[key]["data"][:,iens,:,:,:].reshape(tmshape[0],tmshape[2],ps.shape[0],tmshape[4])[:,:,pindex,:]
                picopy(currentdata,currentdatatm,tasks[key]["data"],iens,pindex)
            elif tasks[key]["shortname"] in ["fg_dep"]: 
                currentdata[:]=currentdatatm-tasks[key]["data"].reshape(tmshape[0],tmshape[2],ps.shape[0],tmshape[4])[:,:,pindex,:]
            elif tasks[key]["shortname"] == "bgdiff":  
                currentdata[:]=currentdatatm-currentdatabg
            elif tasks[key]["shortname"] in ["uah",'rss','star']:  
                #currentdata=numpy.asarray([0.])
                sat=True
            elif tasks[key]["shortname"] in ['suny','sunyhom', 'rharm', 'rharm_h', 'rharmbc']: 
                picopy(currentdata,currentdatatm-currentdatatm,tasks[key]["data"],iens,pindex)
            else:
                print(('key '+tasks[key]["shortname"] +' not processed'))
                return

            if slowload:
                expandandadd(tasks[key]["msudata"],currentdatamsu,tasks[keytm]["index"],msupindex,iens,currentdatamsu,0.0)

            if len(gstations)>0:
                try:
                    for stn in range(stnames.shape[0]):
                        if stnames[stn] not in gstations:
                            currentdata[stn,:,:,:]=numpy.nan
                except:
                    pass
                

            try:
                if daynight=='Day':
                    mask=numpy.logical_and(lons>-90.,lons<90.)
                    currentdata[mask,0,:,:]=numpy.nan
                    currentdata[numpy.logical_not(mask),1,:,:]=numpy.nan		    
                elif daynight=='Night':
                    mask=numpy.logical_and(lons>-90.,lons<90.)
                    currentdata[mask,1,:,:]=numpy.nan
                    currentdata[numpy.logical_not(mask),0,:,:]=numpy.nan
                else:
                    pass
            except:
                pass

            print(('expand',iens,time.time()-t))

            if intervals.ndim>1:
                intrange=intervals.shape[0]
            else:
                intrange=1
            for inter in range(intrange):
                interval=numpy.asarray(intervals[:],dtype=numpy.int32)
                stop=interval[1]
                start=interval[0]
                stime=(start-startyear)*12
                if stime <0:
                    continue
                itime=stime+numpy.arange((stop-start+1)*12)
                if iens==0:
                    sli+=1

                from fnmatch import fnmatch
                if len(plotproperties['sonde_mask']) > 0 and plotproperties['sonde_mask'][0] != b"":
                    try:
                        i = 0
                        for stn in range(stnames.shape[0]):
                            if stnames[stn] == '035229':
                                print('x')
                            u, bc = np.unique(tasks[key]['sonde_type'][i, itime], return_counts=True)
                            st = u[np.argmax(bc)].strip()
                            match = False
                            for p in plotproperties['sonde_mask']:
                                if fnmatch(st,p):
                                    match = True
                                    break
                                
                            if not match:
                                currentdata[stn,:, :, itime]=numpy.nan
                            i += 1
                    except:
                        pass
                t=time.time()

                s=slopes[sli]["data"]

                if not sat:
#		    jcurrentdata=numpy.concatenate((currentdata,currentdatamsu),axis=2)
#		    print 'if',time.time()-t
                    cs=currentdata.shape
                    if 'satseries' in plotproperties['plotlist']:
                        if jcurrentdata.shape[2]!=cs[2]+currentdatamsu.shape[2]:
                            jcurrentdata=numpy.zeros((cs[0],cs[1],cs[2]+currentdatamsu.shape[2],cs[3]))
                        ncat(currentdata,currentdatamsu,jcurrentdata)
                        anomalies=jcurrentdata[:]

                        jpindex=numpy.concatenate((pindex,msupindex))
                    else:
                        jcurrentdata=currentdata
                        anomalies=currentdata[:]
                        jpindex=pindex
                else:
                    jcurrentdata=currentdatamsu[:]
                    anomalies=currentdatamsu[:]
                    jpindex=msupindex
                    #anomalies_and_slopes(currentdatamsu,startyear,interval,tolerance,iens,itime,orig,anomaly,anomalies,climatology,
                                        #good,s)

                good=numpy.zeros([tmshape[0],tmshape[2],jpindex.shape[0]],numpy.int32)
                print((time.time()-t))
                if tasks[key]['shortname'] in andepstdlist or ref=='zero':
                    #jcurrentdata[jcurrentdata==0.0]=numpy.nan
                    calc_diffs(jcurrentdata,numpy.zeros(jcurrentdata.shape),startyear,interval,tolerance,iens,good,s)
                else:

                    if key<2 and (tasks[key]["shortname"]=='tm' and ref!='tm' or tasks[key]["shortname"]=='tmcorr' and ref!='tmcorr') and ref!='night':
                        tasks[key]['zslopes'][:]=numpy.nan
                        tasks[key]['beltslopes'][:]=numpy.nan
                        jcurrentdatatm=jcurrentdata.copy()
                        break
                    if tasks[key]["shortname"] == ref and inter==0:
                        jcurrentdatabg=jcurrentdata.copy()
                        break
                    else:
                        if ref=='night':
                            calc_diffs(jcurrentdata[:,::-1,:,:],jcurrentdata,startyear,interval,tolerance,iens,good,s)
                        else:			    
                            calc_diffs(jcurrentdata,jcurrentdatabg,startyear,interval,tolerance,iens,good,s)

#		anomaliesd_diffs(jcurrentdata,startyear,interval,int(tolerance),int(iens),itime,orig,anomaly,anomalies,
#		                       climatology,good,s)

    # remove India               
                mnames=read_mesural(os.path.expanduser('~/tables/MESURAL'))   

                mnames=[]
                for m in mnames:
                    idx=numpy.where(m==stnames)[0]
                    if len(idx)>0:
                        anomalies[idx[0],:,:,:]=numpy.nan
                        s[idx[0],:,:,:]=numpy.nan
                #mask=(istnames>61900)&(istnames<62000)
                #s[mask,:,:,:]=numpy.nan
                #anomalies[mask,:,:,:]=numpy.nan
                #mask=(istnames>48000)&(istnames<49000)
                #s[mask,:,:,:]=numpy.nan
                #anomalies[mask,:,:,:]=numpy.nan
                ##mask=(istnames<50000)|(istnames>60000)
                ##s[mask,:,:,:]=numpy.nan
                ##anomalies[mask,:,:,:]=numpy.nan
                #mask=(istnames>42000)&(istnames<44000)
                #s[mask,:,:,:]=numpy.nan
                #anomalies[mask,:,:,:]=numpy.nan

                if False:
                    ct=numpy.sum(~numpy.isnan(anomalies),axis=(0,1))
                    ct[ct>50]=50
                    if anomalies.shape[2]>1:
                        clevs=[0,2,5,10,20,50]
                        plt.contourf(numpy.arange(anomalies.shape[3])/12.+1900,ps[pindex],ct,levels=clevs,cmap='gist_rainbow')
                        plt.ylim(1000,10)
                        ax=plt.gca()
                        ax.set_yscale('log')
                        plt.yticks((1000,850,700,500,300,200,100,50,20,10),('1000','850','700','500','300','200','100','50','20','10'))
                        plt.ylabel('hPa')
                        plt.xlim(1950,2017)
                        plt.title('# of monthly means over China')
                        plt.colorbar()
                    else:    
                        plt.plot(numpy.arange(anomalies.shape[3])/12.+1900,ct.flatten())

                #mask=(istnames>48600)&(istnames<48900)
                #s[mask,:,:,:]=numpy.nan
                #anomalies[mask,:,:,:]=numpy.nan


#                if tasks[key]["shortname"]=='rcorr' or tasks[key]["shortname"]=='bgdiff':
#		    northof30N(tasks,key,corig,path,ps,pindex)
                for ipar in (0,1):
                    idy=numpy.where(~numpy.isnan(s[:,0,ipar,0]))[0]
                    dists2=numpy.zeros((idy.shape[0]+1)*(idy.shape[0]+2)//2)
                    tdist(dists2,lats[idy],lons[idy],0.)
                    idx=numpy.zeros(idy.shape,dtype=numpy.int32)
                    tlist=[]
                    for l in range(idy.shape[0]-2):
                        idx[l+1]=idx[l]+idy.shape[0]-l
                        m=numpy.argmin(dists2[idx[l]+1:idx[l+1]])
                        if dists2[idx[l]+m+1]*6370<50.:
                            print(('twins:',stnames[idy[l]],stnames[idy[l+1+m]],dists2[idx[l]+1+m]*6370))
                            tlist.append((idy[l],idy[l+1+m]))
                            if numpy.sum(~numpy.isnan(currentdata[idy[l],ipar,0,:]))>numpy.sum(~numpy.isnan(currentdata[idy[l+1+m],ipar,0,:])):
                                s[idy[l+1+m],0,ipar,0]=numpy.nan
                                #currentdata[idy[l+1+m],ipar,0,:]=numpy.nan
                            else:
                                s[idy[l],0,ipar,0]=numpy.nan
                                #currentdata[idy[l],ipar,0,:]=numpy.nan



                s[:,:,2,:]=numpy.nanmean(s[:,:,:2,:],axis=2)
                print(('slope',time.time()-t))
                if slowload or tasks[key]['shortname']==plotproperties['monthlytasks'][0] or 'gridded' in plotproperties['plotlist']:

                    if tasks[key]['shortname']=='rio':
                        print('rio')
                    t=time.time()

                    itime=numpy.asarray(itime,dtype=int)
                    try:
                        if gslopes.shape[3]!=jpindex.shape[0]:
                            gslopes=numpy.zeros([nj,ni,s.shape[2],jpindex.shape[0]],numpy.float64)
                    except:
                        gslopes=numpy.zeros([nj,ni,s.shape[2],jpindex.shape[0]],numpy.float64)

                    grid_diffs(s,good,int(tolerance),gstatindex,gslopes,iens)


                    print(('grid',time.time()-t))
                    g=gslopes.shape
                    t=time.time()

                    zanomalies=numpy.nanmean(numpy.nanmean(gslopes,axis=1),axis=1)
                    zs=zanomalies.shape
                    beltanomalies=numpy.zeros([len(beltnames),zs[1]])
                    for b in range(len(beltnames)):
                        beltanomalies[b,:]=numpy.nanmean(zanomalies[belts[b,0]:belts[b,1],:],axis=0)

    #		tasks[key]['beltanomalies'][:]=beltanomalies[:]
                    #plt.plot(1900.+itime/12.,beltanomalies[0,15,itime])
                    print((time.time()-t))

                pmax=3
                rprefix=''
                if ref=='night':
                    pmax=1
                    rprefix='rad'
                for ipar in range(pmax):
                    parstr="{0:0>2}".format(ipar*12)
                    for ip in range(jpindex.shape[0]):

                        if ip<pindex.shape[0] and jpindex.shape[0]!=msupindex.shape[0]:
                            pstr="{0:0>3}".format(ps[jpindex[ip]].astype(numpy.int32))
                            psuff=' hPa'
                        else:
                            pstr=msups[jpindex[ip]]
                            psuff=''

                        estr=''
                        if len(tasks[key]["ens"])>1:
                            estr="{0:0>2}".format(iens)

                        if plotproperties['map']=='Globe':
                            ppage = page(
                                layout='positional',  
                                page_x_length=29.7, 
                                page_y_length=21., 
                                page_id_line='off',
                                page_x_position=0., 
                                page_y_position=0.)

                            pmap=''
                        else:
                            ppage = page(
                                layout='positional',  
                                page_x_length=14.7, 
                                page_y_length=20., 
                                page_id_line='off',
                                page_x_position=0., 
                                page_y_position=0.)
                            pmap='_'+plotproperties['map']

                        nprefix=''
                        #if "gridded" in plotproperties["plotlist"]:
                            #nprefix+='gridded_'
                        if "stations" in plotproperties["plotlist"] or'gridded' in plotproperties["plotlist"] :
                            nprefix+='stations_'

                        oname=nprefix+'diff_'+tasks[key]["shortname"]+estr+'_'+ref+'_'+\
                            "{:4}".format(interval[0])+"-{:4}".format(interval[1])+'_'+\
                            pstr+'_'+parstr  
                        out = output({"output_name_first_page_number":"off",
                                      "output_formats":plotproperties["outputformat"], 
                                      'output_name': plotproperties['tmppath']+'/'+oname})
                            #define the cartesian projection

                        print(('stations output file: ',oname))

                        if tasks[key]['shortname'] in andepstdlist or tasks[key]['shortname'] in andeprmslist:
                            clev=numpy.linspace(0.0,4.0,21)
                            clev=numpy.append(clev,10.)
                        else:
                            clev=numpy.linspace(-1.2,1.2,25)
                            clev=numpy.append(clev,5.)
                            clev=numpy.append(-5.,clev)

                        if plotproperties['double']:
                            clev*=2

                        if "stations" in plotproperties["plotlist"] or 'cost' in plotproperties["plotlist"] or 'gridded' in plotproperties["plotlist"] :

                            scosts=numpy.zeros(s.shape[0])
                            tt=time.time()
                            cost=tcost(dists,s[:,iens,ipar,ip],scosts)
                            #print time.time()-tt
                            #cost2=tcost2(dists,s[:,iens,ipar,ip],scosts)
                            #print time.time()-tt
                            if ip<pindex.shape[0]:
                                tasks[key]['cost'][iens,ipar,ip]=cost

                            if ip==jpindex.shape[0]-1 and ipar==1:
                                print(('cost',time.time()-t))
                            cstr="{0:8.2f}".format(cost)
                            statstr="{:4}".format(sum(~numpy.isnan(s[:,iens,ipar,ip])))
                            if '_std' in tasks[key]['shortname']:
                                lines =["Difference Standard Deviation [K], "+tasks[key]["shortname"]+estr+
                                        ", {:4}".format(interval[0])+"-{:4}".format(interval[1])+', '+pstr+psuff]
                            elif '_corr' in tasks[key]['shortname']:
                                lines =["Anomaly Correlation, "+tasks[key]["shortname"]+estr+
                                        ", {:4}".format(interval[0])+"-{:4}".format(interval[1])+', '+pstr+psuff]
                            else:
                                if ref=='night':
                                    tref='12-00'
                                    lines =["Temperature Difference [K], "+tasks[key]["shortname"]+estr+',12h-00h, '+plotproperties['exps'][interv]+', '+
                                            "{:4}".format(interval[0])+"-{:4}".format(interval[1])+', '+pstr+psuff]
                                else:
                                    tref=ref
                                    lines =["Temperature Difference [K], "+tasks[key]["shortname"]+estr+'-'+ref+', '+plotproperties['exps'][interv]+', '+
                                            "{:4}".format(interval[0])+"-{:4}".format(interval[1])+', '+parstr+"h, "+pstr+psuff]
                            lines.append(str(statstr+' Stations, Cost: '+cstr+', '+plotproperties["version"]))



                            vec=s[0:istat,iens,ipar,ip].flatten()
                            mask=~numpy.isnan(vec)
                            if sum(mask)<2:
                                continue
                            inp=minput({'Input_values':vec[mask].astype(numpy.float64),  #.ravel?
                                        'Input_latitude_values':lats[mask].astype(numpy.float64),
                                        'Input_longitude_values':lons[mask].astype(numpy.float64),
                                        'Input_Type':'geographical'})

                        else:
                            lines =["Temperature Difference [K], "+tasks[key]["shortname"]+estr+', '+
                                    "{:4}".format(interval[0])+"-{:4}".format(interval[1])+', '+parstr+"h, "+pstr+psuff,
                                    plotproperties["version"]]

                        if "stations" in plotproperties["plotlist"]:
                            cind=numpy.argsort(-scosts[mask])
                            maxs=cind.shape[0]
                            if str(plotproperties['plotstatids'])=='False':
                                maxs=3
                                ilav=numpy.concatenate([[],numpy.asarray([-180.,-180.,-180.])+1.0])
                                ilov=numpy.concatenate([[],numpy.asarray([-20.,0.,20.])+1.0])
                            else:
                                ilav=numpy.concatenate([lats[mask][cind[0:maxs]].astype(numpy.float64),
                                                        numpy.asarray([-80.,-80.,-80.])+1.0])
                                ilov=numpy.concatenate([lons[mask][cind[0:maxs]].astype(numpy.float64),
                                                        numpy.asarray([-20.,0.,20.])+1.0])

                            projection,coastlines,title,legend,symb,symb2,symb3,symb4,cont=\
                                set_trendmaps(lines,clev,plotproperties,
                                              stnames=numpy.concatenate([stnames[mask][cind[:maxs]],stnames[mask][cind[:maxs]]]),
                                              slopes=numpy.concatenate([vec[mask][cind[:maxs]].astype(numpy.float64),vec[mask][cind[:maxs]]]),
                                              costs=numpy.concatenate([scosts[mask][cind[:maxs]],scosts[mask][cind[:maxs]]]),
                                              map=plotproperties['map'])

                            inp3=minput({'Input_latitude_values':ilav,
                                         'Input_longitude_values':ilov,
                                         'Input_Type':'geographical'})
                            inp4=minput({'Input_latitude_values':ilav-2.0,
                                         'Input_longitude_values':ilov+2.0,
                                         'Input_Type':'geographical'})
                            inp5=minput({'Input_latitude_values':ilav-4.0,
                                         'Input_longitude_values':ilov+2.0,
                                         'Input_Type':'geographical'})
                        else:
                            projection,coastlines,title,legend,symb,cont=set_trendmaps(lines,clev,plotproperties,
                                                                                       map=plotproperties['map'])

                        if "gridded" in plotproperties["plotlist"]:
                            glats=numpy.linspace(-94.999,94.999,nj+2)
                            glons=numpy.linspace(5,375,ni+2)

                            hilf=numpy.empty([nj+2,ni+2])
                            hilf[:]=numpy.nan

                            if ipar==0 and ip==9:
                                print((out.args['output_name']))
                            hilf[1:nj+1,1:ni+1]=gslopes[:,:,ipar,ip].reshape([nj,ni]).copy()
    #                        hilf[1:nj+1,1:ni+1]=hadmedslopes[0,:,:].reshape([nj,ni]).copy()
                            mask2=numpy.isnan(hilf)
                            hilf[mask2]=-1.e21

                            inp2=minput({'Input_field':hilf[:,1:ni+1],  #.ravel?
                                         'Input_field_initial_latitude':glats[0]-5,
                                         'Input_field_initial_longitude':glons[0]-5,
                                         'Input_field_latitude_step':(glats[1]-glats[0]),
                                         'Input_field_longitude_step':glons[1]-glons[0],
                                         'Input_Type':'geographical'})

                        #plot(out,projection,coastlines,inp,symb,title,legend,inp2,cont)
                            out.args['output_title']=os.getcwd()+'/'+out.args['output_name']
                            if len(vec[mask])>0 and "stations" in plotproperties["plotlist"] or "cost" in plotproperties["plotlist"]:
                                plot(out,projection,inp2,cont,inp,symb,title,legend,coastlines)
                            else:
                                plot(out,projection,inp2,cont,title,legend,coastlines)
                                
                        else:
                            if "stations" in plotproperties["plotlist"]:
                                if numpy.sum(mask)!=0:
                                    out.args['output_title']=os.getcwd()+'/'+out.args['output_name']
#				    if plotproperties['plotstatids']=='True':
                                    plot(out,projection,inp,symb,inp3,symb2,inp4,symb3,inp5,symb4,title,legend,coastlines)
#				    else:
                                        #print 'plot inp:',inp
#					plot(out,projection,inp,symb,title,legend,coastlines)
                                else:
                                    print((out.args['output_name']+': no valid data'))


                if nodata[-1]:
                    print((tasks[key]["shortname"]+': does not contain valid data, will be removed'))
                    continue

                if not slowload:
                    if tasks[key]['shortname']==plotproperties['monthlytasks'][0] and '':
                        return
                    continue
                t=time.time()
                s=gslopes.shape
                zslopes=numpy.empty([s[0],s[3]])
                zslopes.fill(numpy.nan)
                zslopes=zonaltrends(gslopes,zslopes)
                tasks[key]["zslopes"][iens,:,:pindex.shape[0]]=zslopes[:,:pindex.shape[0]]

                print(('zslopes',time.time()-t))
                if 'zonal' in plotproperties['plotlist'] and zslopes.shape[1]>=pindex.shape[0] and pindex.shape[0]>1:
                    out = output({"output_name_first_page_number":"off",
                                  "output_formats":plotproperties["outputformat"], 
                                  'output_name': plotproperties['tmppath']+'/'+'diffzonal_'+tasks[key]["shortname"]+estr+'-'+ref+'_'+
                                  "{:4}".format(interval[0])+"-{:4}".format(interval[1])})


                    mask=numpy.isnan(zslopes)

                    hilf=zslopes.T#[pindex,:]
                    hilf=hilf[:len(pindex),:]

                    lines =["Temperature Difference [K], "+tasks[key]["shortname"]+estr+'-'+ref+', '+plotproperties['exps'][interv]+', '+
                            "{:4}".format(interval[0])+"-{:4}".format(interval[1]),
                            stats(hilf[hilf==hilf],mima=1,short=1)]

                    hilf[numpy.isnan(hilf)]=-1.e21

                    projection,horizontal,vertical,title,legend,cont=set_zonalcross(lines,
                                                                                    clev,plotproperties)

                    ppage = page(
                        layout='positional',  
                        page_x_length=29.7, 
                        page_y_length=21., 
                        page_id_line='off',
                        page_x_position=0., 
                        page_y_position=0.)
                    out.args['output_title']=os.getcwd()+'/'+out.args['output_name']

                    xa=glats#[1:19]
                    ya=numpy.asarray(ps[numpy.asarray(pindex,dtype=int)],numpy.float64)
                    inputs=minput({'Input_field':hilf,  #.ravel?
                                   'Input_field_x_list':xa,
                                   'Input_field_y_list':ya,
                                   'Input_type':'cartesian'
                                   })
                    numpy.savez(out.args['output_name'],hilf=hilf,xa=xa,ya=ya)
                    plot(out,ppage,projection,horizontal,vertical,inputs,cont,title,legend)
    if "belts" in plotproperties['plotlist']:

        for ibelt in range(len(beltnames)):

            ppage = page(
                layout='positional',  
                page_x_length=29.7, 
                page_y_length=21., 
                page_id_line='off',
                page_x_position=0., 
                page_y_position=0.)

            linesplustext=[]
            legend_user_lines=[]
            lines =["Temperature Difference [K], "+beltnames[ibelt]+', ref='+ref+
                    ", {:4}".format(interval[0])+"-{:4}".format(interval[1]),plotproperties['version']+', '+plotproperties['exps'][0]]

            profsat=''
            if 'satseries' in plotproperties['plotlist']:
                profsat='sat'

            out = output({"output_name_first_page_number":"off",
                          "output_formats":plotproperties["outputformat"], 
                          'output_name': plotproperties['tmppath']+'/'+profsat+'diffsbelt_'+beltnames[ibelt]+'_'+
                          "{:4}".format(interval[0])+"-{:4}".format(interval[1])})

            ti=list(range(len(tasks)))
            li=[ta['shortname'] for ta in tasks]
            if 'rio' in li:
                ti.insert(0,ti.pop(li.index('rio')))
            if 'ce20c_andep' in li:
                ti.insert(0,ti.pop(li.index('ce20c_andep')))
            #try:
                #iref=li.index(ref)
                #ti.pop(iref)
                #li.pop(iref)
            #except:
                #pass

            for key in ti:

                for iens in tasks[key]["ens"]: #range(len(tasks[key]["ens"])):
                    if ibelt==0:
                        tasks[key]["beltslopes"][iens,:,:]=belttrends(tasks[key]["zslopes"][iens,:,:],belts)

                shade=1
                if numpy.sum(~numpy.isnan(tasks[key]["beltslopes"][0,:,:pindex.shape[0]]))>0 and \
                   tasks[key]["shortname"]!=ref and tasks[key]["shortname"]!='eijra_fgdep':
                    legname=copy.copy(tasks[key]["shortname"])
                    if tasks[key]["shortname"] not in [ "rio", 'ce20c_andep' ] or shade==0:
                        if 'andep' in legname:
                            legname=legname.split('_andep')[0]
                        for iens in tasks[key]["ens"]: #range(len(tasks[key]["ens"])):

                            if 'satseries' in plotproperties['plotlist']:
                                if tasks[key]["shortname"] != "uah" and tasks[key]["shortname"] != "rss":
                                    tasks[key]["beltslopes"][iens,ibelt,-4]=numpy.nan
                                linesplustext=linesplustext+addsattrend(tasks[key]["beltslopes"][iens,ibelt,-4:],[800.,550.,250.,90.],
                                                                        scol=tasks[key]["color"],iens=iens) # SAT
                                #amp=tasks[key]["beltslopes"][iens,ibelt,-3]/hadmedbeltslopes[ibelt,0]
                            else:
                                if tasks[key]["beltslopes"].shape[2]>=pindex.shape[0]:
                                    linesplustext=linesplustext+addprofile(tasks[key],ps[pindex],iens,ibelt) # RASO
                                else:
                                    linesplustext=linesplustext+addsattrend(tasks[key]["beltslopes"][iens,ibelt,:],[800.,550.,250.,90.],
                                                                            scol=tasks[key]["color"],iens=iens) # SAT
                    else:
                        if 'satseries' in plotproperties['plotlist']:
                            for iens in tasks[key]["ens"]: #range(len(tasks[key]["ens"])):
                                tasks[key]["beltslopes"][iens,ibelt,-4]=numpy.nan
                                linesplustext=linesplustext+addsattrend(tasks[key]["beltslopes"][iens,ibelt,-4:],[800.,550.,250.,90.],
                                                                        scol=tasks[key]["color"],iens=iens) # SAT
                                print((tasks[key]["beltslopes"][iens,ibelt,-4:]))
                        else:
                            linesplustext=linesplustext+add_ens_profile(tasks[key], ps[pindex], ibelt)

                    legend_user_lines=legend_user_lines+[legname]

            projection,horizontal,vertical,title,legend=set_belttrends(lines,legend_user_lines,plotproperties,beltinterval=[-1.8,0.6])

            if len(linesplustext)>0:
                out.args['output_title']=os.getcwd()+'/'+out.args['output_name']	
                plot(out,ppage,projection,horizontal,vertical,title,linesplustext,legend)
                print((out.args['output_name']+'.'+plotproperties["outputformat"][0]))
            else:
                print((out.args['output_name']+'.'+plotproperties["outputformat"][0]+ 'not created - no data found'))

    if 'stations' in plotproperties['plotlist'] or 'cost' in plotproperties['plotlist']:
        for key in range(len(tasks)):
            if key<2 and (tasks[key]["shortname"]=='tm' and ref!='tm' or tasks[key]["shortname"]=='tmcorr' and ref!='tmcorr') and ref!='night':
                tasks[key]['zslopes'][:]=numpy.nan
                tasks[key]['beltslopes'][:]=numpy.nan
        if slowload:
            makecostprofiles(tasks,nodata,[],interval,plotproperties,'diff',ps,pindex,iens,msups,msupindex)

    return           	

