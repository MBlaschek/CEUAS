import os,glob,sys
import numpy as np
import matplotlib
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
#from mpl_toolkits.basemap import Basemapnp.min
import matplotlib.pylab as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
import matplotlib.colors as colors
from rasotools.additions.dumpandspagarr import *
from rasotools.additions.define_datasets import *
import cartopy.crs as ccrs
import traceback
import copy
from functools import partial
from multiprocessing import Pool
import shutil
import h5py
import pandas as pd

from netCDF4 import Dataset
from numba import njit,prange


def prepare_ti(tasks,rlist=['rio','rit','ce20c_andep']):    
    ti=list(range(len(tasks)))
    li=[ta['shortname'] for ta in tasks]
    for r in rlist:
        if r in li:
            ti.insert(0,ti.pop(li.index(r)))

    return ti

@njit
def fcp(outVar,ghilf,idx,ip):
    for it in range(ghilf.shape[3]):
        for ilon in range(ghilf.shape[1]):
            for ilat in range(ghilf.shape[0]):
                outVar[it,idx,ilat,ilon]=ghilf[ghilf.shape[0]-ilat-1,(ilon+ghilf.shape[1]//2)%ghilf.shape[1],ip,it]

# as append_gridded, but up to 10 hPa
def append_gridded_10(ifile,ofile,ganomalies,days,ps,start=2006,stop=2015,version='1.9',append=True):
    #input file
    dsin = Dataset(ifile)
    dsin.set_auto_maskandscale(False)    

    #output file
    try:
        os.remove(ofile)
    except:
        pass
    dsout = Dataset(ofile, "w")
    dsout.set_auto_maskandscale(False)    
    #Copy dimensions
    for dname, the_dim in list(dsin.dimensions.items()):
        print((dname, len(the_dim)))
        if dname=='time':
            dsout.createDimension(dname, ganomalies.shape[4])
        elif dname=='pressure':
            dsout.createDimension(dname, ps.shape[0])
        else:
            dsout.createDimension(dname, len(the_dim))


    #Copy variables
    press=ps[::-1]
    t=time.time()
    ghilf=np.nanmean(ganomalies,axis=2)
    print(('toout:',ofile,time.time()-t))
    for v_name, varin in list(dsin.variables.items()):
        if v_name=='anomalies' and ganomalies.shape[4]==12:
            outVar = dsout.createVariable('climatology', varin.datatype, varin.dimensions)
        else:
            outVar = dsout.createVariable(v_name, varin.datatype, varin.dimensions)
        print((v_name,varin.datatype))
        if v_name=='time':
            outVar[:] = days[:]
        elif v_name=='pressure':
            outVar[:] = press[:]
        elif v_name=='anomalies':
            ov=np.empty(outVar.shape)
            ov.fill(np.nan)
            for ip in range(ps.shape[0]):
                if ps[ip] in press:
                    idx=np.where(press==ps[ip])[0][0]

                    fcp(ov, ghilf, idx, ip)

            ov[np.isnan(ov)]=-1.e30
#	    ov[696:696+varin.shape[0],:,:,:]=varin[:]
            outVar[:]=ov

            print(('toout:',time.time()-t))

        else:
            outVar[:] = varin[:]


        for attname in varin.ncattrs():
            print(attname)
            if attname=='units' and v_name=='time':
                setattr(outVar,attname,"days since 1900-01-01 00:00:00")
            elif v_name=='anomalies' and ganomalies.shape[4]==12 and attname in ['valid_min','valid_max','long_name']:
                setattr(outVar,'valid_min',160.)
                setattr(outVar,'valid_max',310.)
                setattr(outVar,'long_name','climatologies')
            else:
                setattr(outVar,attname,getattr(varin,attname))

    for attname in dsin.ncattrs():
        #print(attname)
        setattr(dsout,attname,getattr(dsin,attname))

    setattr(dsout,'history',str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
    tstring='{}-{}'.format(start,stop).join(getattr(dsin,'title').split('1981-2010'))
    if 'rio' in ofile:
        tstring='RICH-obs v1.5.1'+' member '+ ofile[3:5].join(tstring.split('RAOBCORE v1.5.1'))
    if 'rit' in ofile:
        tstring='RICH-tau v1.5.1'+' member '+ ofile[3:5].join(tstring.split('RAOBCORE v1.5.1'))
    if 'Control' not in os.getcwd():
        tstring=version.join(tstring.split('1.5.1'))
    if 'raw' not in ofile:
        tstring=''.join(tstring.split('raw '))
        
    if ganomalies.shape[4]==12:
        tstring=''.join(tstring.split('anomalies with respect to '))
    setattr(dsout,'title',str(tstring))
    print(tstring)
    #close the output file
    dsout.close()

def append_gridded(ifile,ofile,ganomalies,days,ps):
    #input file
    dsin = Dataset(ifile)
    dsin.set_auto_maskandscale(False)    

    #output file
    try:
        os.remove(ofile)
    except:
        pass
    dsout = Dataset(ofile, "w")
    dsout.set_auto_maskandscale(False)    
    #Copy dimensions
    for dname, the_dim in list(dsin.dimensions.items()):
        print((dname, len(the_dim)))
        if dname!='time':
            dsout.createDimension(dname, len(the_dim))
        else:
            dsout.createDimension(dname, ganomalies.shape[4])


    #Copy variables
    press=dsin.variables['pressure'][:]
    t=time.time()
    ghilf=np.nanmean(ganomalies,axis=2)
    print(('toout:',time.time()-t))
    for v_name, varin in list(dsin.variables.items()):
        outVar = dsout.createVariable(v_name, varin.datatype, varin.dimensions)
        print((v_name,varin.datatype))
        if v_name=='time':
            outVar[:] = days[:]
        elif v_name=='anomalies':
            ov=np.empty(outVar.shape)
            ov.fill(np.nan)
            for ip in range(ps.shape[0]):
                if ps[ip] in press:
                    idx=np.where(press==ps[ip])[0][0]

                    fcp(ov, ghilf, idx, ip)

            ov[np.isnan(ov)]=-1.e30
            ov[696:696+varin.shape[0],:,:,:]=varin[:]
            outVar[:]=ov

            print(('toout:',time.time()-t))

        else:
            outVar[:] = varin[:]


        for attname in varin.ncattrs():
            print(attname)
            if attname=='units' and v_name=='time':
                setattr(outVar,attname,"days since 1900-01-01 00:00:00")
            else:
                setattr(outVar,attname,getattr(varin,attname))

    for attname in dsin.ncattrs():
        print(attname)
        if attname=='history':
            setattr(dsout,attname,str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
        else:
            setattr(dsout,attname,getattr(dsin,attname))

    #close the output file
    dsout.close()


def read_hadCRUT4(path,prefix,startyear=1957,endyear=2013,ens=0):

    t1=time.time()
    try:
        if 'T.4' in prefix:
            f=np.load('CRUT4_'+str(startyear)+'_'+str(endyear)+'.npz')
        else:
            f=np.load('CRUT5_'+str(startyear)+'_'+str(endyear)+'.npz')
            
        hadmed=f["hadmed"]
        hadtem=f["hadtem"]
        hadens=f["hadens"]
    except:

        if 'T.4' in prefix:
            fn=path+prefix+".median.nc"
            fvar='temperature_anomaly'
            fevar='temperature_anomaly'
        else:
            fn=path+prefix+".analysis.anomalies.ensemble_mean.nc"
            fvar='tas_mean'
            fevar='tas'
        try:
            f = netCDF4.Dataset(fn,"r")
            f.set_auto_mask(False)
        except:
            print((fn+' not found'))

        x=f.variables[fvar][:]
        hadmed=np.zeros([1,(endyear-startyear+1)*12,x.shape[1]//2,x.shape[2]//2],dtype=np.float32)
        had_rebin_3672_to_1836(x,hadmed,0,startyear,endyear,np.float32(-1.e30),np.float32('nan'))

        fn=path+'CRUTEM.'+'.'.join(prefix.split('.')[1:])+".anomalies.nc"
        try:
            f = netCDF4.Dataset(fn,"r")
            f.set_auto_mask(False)
            x=f.variables['temperature_anomaly'][:]
            hadtem=np.zeros([1,(endyear-startyear+1)*12,x.shape[1]//2,x.shape[2]//2],dtype=np.float32)
            had_rebin_3672_to_1836(x,hadtem,0,startyear,endyear,np.float32(-1.e30),np.float('nan'))
        except:
            hadtem=np.zeros_like(hadmed)
            print((fn+' not found'))



#	had_rebin_pythran_3672_to_1836(x,hadmed,0,startyear,endyear,np.float(-1.e30),np.float('nan'))
        enslist=glob.glob(path+prefix+'*'+'anomalies.*[0-9].nc')
        hadens=np.zeros([len(enslist),hadmed.shape[1],x.shape[1]//2,x.shape[2]//2],dtype=np.float32)
        if ens!=0:
            for ens in range(len(enslist)):
                #fn=os.path.join(path,prefix+'.anomalies.'+'{0:1}'.format(ens+1)+".nc")
                fn=enslist[ens]
                try:
                    f = netCDF4.Dataset(fn,"r")
                except:
                    print((fn+' not found'))
                print(ens)
                x=f.variables[fevar][:]
                had_rebin_3672_to_1836(x,hadens,ens,startyear,endyear,np.float32(-1.e30),np.float32('nan'))

        if 'T.4' in prefix:
            np.savez('CRUT4_'+str(startyear)+'_'+str(endyear)+'.npz', hadmed=hadmed,hadtem=hadtem,hadens=hadens)
        else:
            np.savez('CRUT5_'+str(startyear)+'_'+str(endyear)+'.npz', hadmed=hadmed,hadtem=hadtem,hadens=hadens)
            

    print((time.time()-t1))

#    plt.contour(np.reshape(hadens[0,0,:,:],[18,36]))
    #sys.exit()
    return hadmed,hadtem,hadens

def read_rss(path,version,startyear=1957,endyear=2013,ens=0):

    t1=time.time()
    chlist=['TLT','TMT','TTS','TLS']
    chlist=['TMT','TTS','TLS']
    try:
        f=np.load('rss_'+str(startyear)+'_'+str(endyear)+'.npz')
        x=f["rssfull"]
        rss18=f["rss18"]
    except:

        l=1
        for ch in chlist:
            fn=path+'RSS_Tb_Anom_Maps_ch_'+ch+'_'+version+'.nc'
            try:
                f = netCDF4.Dataset(fn,"r")
                f.set_auto_mask(False)
            except:
                print((fn+' not found'))

            try:
                dstart=getattr(f.variables['months'],'units')
                syear=int(dstart.split()[2].split('-')[0])
                s=f.variables['brightness_temperature'].shape
                fill=getattr(f.variables['brightness_temperature'],'_FillValue')
                si=(syear-startyear)*12
                if 'x' not in locals():
                    x=np.empty([4,(endyear-startyear+1)*12,s[1],s[2]])
                    x[:,:,:,:]=np.nan
                x[l,si:(si+s[0]),:,:]=f.variables['brightness_temperature'][:(endyear-syear+1)*12,:,:]
            except NameError:
                if 'x' not in locals():
                    x=np.empty([4,(endyear-startyear+1)*12,s[1],s[2]])
                x[l,si:(si+s[0]),:,:].fill(np.nan)
            l+=1

        x[x==fill]=np.nan
    #    x=np.roll(x,s[2]/2,axis=3)

        s=x.shape
        rss18=np.empty([1,s[0],s[1],s[2]//4,s[3]//4])
        rebin_72144_to_1836(x,rss18,0,global_startyear=startyear,startyear=startyear,endyear=endyear)
        np.savez('rss_'+str(startyear)+'_'+str(endyear)+'.npz', rssfull=x,rss18=rss18)

    print(('read_rss: ',time.time()-t1))
#    plt.contour(np.reshape(x[3,1300,:,:],[72,144]))
#    plt.show()
    #sys.exit()
    return x,rss18

def read_star(path,fnlist,startyear=1957,endyear=2013,ens=0):

    t1=time.time()
    chlist=['TMT','TTS','TLS']
    chslist=['TMT','TUT','TLS']
    try:
        f=np.load('star_'+str(startyear)+'_'+str(endyear)+'.npz')
        x=f["starfull"]
        star18=f["star18"]
    except:

        l=1
        for ch,chs,fx in zip(chlist,chslist,fnlist):
            fn=path+fx
            try:
                f = netCDF4.Dataset(fn,"r")
                f.set_auto_mask(False)
            except:
                print((fn+' not found'))

            try:
                dstart=getattr(f.variables['time'],'units')
                syear=int(dstart.split()[2].split('-')[0])
                btvn='tcdr_MSU_AMSUA_ATMS_'+chs
                s=f.variables[btvn].shape
                m=int(f.variables['time'][0]+15)//30-1
                fill=getattr(f.variables[btvn],'_FillValue')
                si=(syear-startyear)*12+m
                smax=(endyear-startyear+1)*12

                if 'x' not in locals():
                    x=np.empty([4,(endyear-startyear+1)*12,s[1],s[2]])
                    x[:,:,:,:]=np.nan
                if si+s[0]>=smax:
                    x[l,si:smax,:,:]=f.variables[btvn][:smax-si,:,:]
                else:
                    x[l,si:si+s[0],:,:]=f.variables[btvn][:,:,:]

            except NameError:
                if 'x' not in locals():
                    x=np.empty([4,(endyear-startyear+1)*12,s[1],s[2]])
                    x[l,si:,:,:].fill(np.nan)
            l+=1

        x[x==fill]=np.nan
    #    x=np.roll(x,s[2]/2,axis=3)

        s=x.shape
        star18=np.empty([1,s[0],s[1],s[2]//4,s[3]//4])
        rebin_72144_to_1836(x,star18,0,global_startyear=startyear,startyear=startyear,endyear=endyear)
        np.savez('star_'+str(startyear)+'_'+str(endyear)+'.npz', starfull=x,star18=star18)

    print(('read_star: ',time.time()-t1))
#    plt.contour(np.reshape(x[3,1300,:,:],[72,144]))
#    plt.show()
    #sys.exit()
    return x,star18

def read_uah(path,version,startyear=1957,endyear=2013,ens=0):

    t1=time.time()
    chlist=['tlt','tmt','tts','tls']
    try:
        f=np.load('uah_'+str(startyear)+'_'+str(endyear)+'.npz')
        x=f["uahfull"]
        uah18=f["uah18"]
    except:
        l=0
        x=np.empty([4,(endyear-startyear+1)*12,72,144])
        x[:]=np.nan
        for ch in chlist:

            for year in range(1979,endyear+1):
                fn=path+ch+'monamg.'+str(year)+'_'+version
                try:
                    with open(fn,"r") as myfile:
    #		    d=myfile.read().split('\n')
                        data=(' -').join(myfile.read().split('-'))
                        data=data.split('\n')
                        dl=(len(data)-1)//649
                        for i in range(dl,-1,-1):
                            dum= data.pop(i*649)
                        data=' '.join(data)
                        si=(year-startyear)*12
                        c=np.fromstring(data,sep=' ')
                        cms=c.shape[0]//72//144
                        x[l,si:si+cms,:,:]=c.reshape([cms,72,144])
                except IOError:
                    #print fn+' not found'
                    continue
            l+=1

        x=np.roll(x,72,axis=3)
    #    x[x==-9999.]=np.nan
        x=x/100.
        s=x.shape
        uah18=np.empty([1,s[0],s[1],s[2]//4,s[3]//4])
        rebin_72144_to_1836(x,uah18,0,global_startyear=startyear,startyear=startyear,endyear=endyear)
        np.savez('uah_'+str(startyear)+'_'+str(endyear)+'.npz', uahfull=x,uah18=uah18)
    print(('read_uah: ',time.time()-t1))
#    plt.contour(np.reshape(x[3,1300,:,:],[72,144]))
#    plt.show()
    #sys.exit()
    return x,uah18

def read_merra2(path,version,ps,startyear=1979,endyear=2017,ens=0):

    t1=time.time()
    try:
        f=np.load('merra2_'+str(startyear)+'_'+str(endyear)+'.npz')
        merra272=f["merra272"]
        merra218=f["merra218"]
    except:
        l=0
#	x=np.empty([16,(endyear-startyear+1)*12,361,576])
#	x[:]=np.nan
        ini=True
        for year in range(1980,endyear+1):
            #for month in range(12):
            fn=path+'MERRA2_t_{:4}'.format(year)+version+'.nc'
            #mv='{:4}-{:0>2}-01'.format(year,month+1)
            try:
                f=netCDF4.Dataset(fn)
#		    print f.variables.keys()
            except IOError:
                print((fn+' not found'))
                continue

            l=(year-startyear)
            y=np.array(f.variables['ta'][:])
            y[y>f.variables['ta']._FillValue/10.]=np.nan
            midx=np.zeros(16,dtype=int)
            for ip in range(16):
                midx[ip]=np.where(f.variables['plev'][:]==ps[ip]*100.)[0]
            #for it in range(12):
                #for ip in range(14):
                    #x[ip,l+it,:,:]=y[it,midx[ip],:,:]
#		ganomalies=np.zeros([nj,ni,2,jpindex.shape[0],tmshape[4]],np.float32)

        ##    x[x==-9999.]=np.nan
            #x=x/100.
            if ini:
                s=y.shape
                merra218=np.empty([1,16,(endyear-startyear+1)*12,s[2]//20,s[3]//16])
                merra272=np.empty([1,16,(endyear-startyear+1)*12,s[2]//5,s[3]//4])
                ini=False
            rebin_361576_to_1836(y,merra218,0,l,midx,global_startyear=startyear,startyear=startyear,endyear=endyear)
            rebin_361576_to_72144(y,merra272,0,l,midx,global_startyear=startyear,startyear=startyear,endyear=endyear)

        merra272=np.reshape(merra272,merra272.shape[1:])		
    np.savez('merra2_'+str(startyear)+'_'+str(endyear)+'.npz',merra218=merra218,merra272=merra272)
    print(('read_merra2: ',time.time()-t1))
#    plt.contour(np.reshape(x[3,1300,:,:],[72,144]))
#    plt.show()
    #sys.exit()
    return merra272,merra218

def rebin_20CRmonth(ps,fn):
    
            #mv='{:4}-{:0>2}-01'.format(year,month+1)
    t1=time.time()
    merra218=np.empty([1,16,1,18,36])
    merra272=np.empty([1,16,1,18*4,36*4])
    try:
        f=h5py.File(fn)
    except IOError:
        print((fn+' not found'))
        merra218.fill(np.nan)
        merra272.fill(np.nan)
        return merra272,merra218

    #l=(year-startyear)*12+month
    #tt=time.time()
    #y=np.empty(f['t'].shape,dtype=np.float32)
    #for l in range(y.shape[0]):
        #y[y.shape[0]-l-1,:]=f['t'][l,:]
        #print(l)
    #print(time.time()-tt)
    #tt=time.time()
    y=np.flip(f['t'][:],axis=0)
    #print(time.time()-tt)
    #y[y>f.variables['ta']._FillValue/10.]=np.nan
    midx=np.zeros(16,dtype=int)
    for ip in range(16):
        midx[ip]=np.where(f['level'][:]==ps[ip])[0]

    #merra272=np.empty([1,16,1,72,144])
        #ini=False
    print(time.time()-t1)
    rebin_256512_to_1836(y,merra272,merra218,0,0,midx)
    #rebin_256512_to_72144(y,merra272,0,l,midx,global_startyear=startyear,startyear=startyear,endyear=endyear)
    print(fn.split('/')[-1],time.time()-t1)

    f.close()
    return merra272,merra218

def read_20CRv3(path,version,ps,startyear=1979,endyear=2017,ens=0):

    t1=time.time()
    try:
        f=np.load('20CRv3_'+str(startyear)+'_'+str(endyear)+'.npz')
        CR20v372=f["CR20v372"]
        CR20v318=f["CR20v318"]
        try:
            msuCR20v372=f["msuCR20v372"]
            msuCR20v372=np.reshape(msuCR20v372,msuCR20v372.shape[1:])		
        except:
            msuCR20v372=None
            pass
    except:

        fns=[]
        for year in range(startyear,endyear+1):
            for month in range(12):
                fns.append(path+'anl_meant_{}{:0>2}_TMP_pres.nc'.format(year,month+1))
                func=partial(rebin_20CRmonth,ps)
        with Pool(40) as p:
            rlist=list(map(func,fns))
        rlist72=[]
        rlist18=[]
        for r in rlist:
            rlist72.append(r[0])
            rlist18.append(r[1])
        CR20v372=np.concatenate(rlist72,axis=2)
        CR20v318=np.concatenate(rlist18,axis=2)
        #merra272=np.reshape(merra272,merra272.shape[1:])		
        np.savez('20CRv3_'+str(startyear)+'_'+str(endyear)+'.npz',CR20v318=CR20v318,CR20v372=CR20v372)#,merra272=merra272)
    print(('read_20CRv3: ',time.time()-t1))
    CR20v372=np.reshape(CR20v372,CR20v372.shape[1:])		
    return CR20v372,CR20v318,msuCR20v372#,merra218

#    plt.contour(np.reshape(x[3,1300,:,:],[72,144]))
#    plt.show()
    #sys.exit()

def read_wegc(path,version,startyear=2001,endyear=2017,ens=0):

    t1=time.time()
    try:
        f=np.load('wegc_'+str(startyear)+'_'+str(endyear)+'.npz')
        x=f["wegcfull"]
        wegc18=f["wegc18"]
    except:
        l=0
        x=np.empty([16,(endyear-startyear+1)*12,72,144])
        x[:]=np.nan
        mv='{:4}-{:0>2}-01'.format(2001,1)
        for year in range(2001,endyear+1):
            for month in range(12):
                fn=path+'2.5x2.5-'+mv+'_{:4}-{:0>2}-01-RO_OPSv'.format(year,month+1)+version+'_L2b.nc'
                mv='{:4}-{:0>2}-01'.format(year,month+1)
                try:
                    f=netCDF4.Dataset(fn)
#		    print f.variables.keys()
                except IOError:
                    try:
                        f=netCDF4.Dataset(fn[:-3]+'_interp'+fn[-3:])
                    except:
    #		    #print fn+' not found'
                        continue

                l=(year-startyear)*12+month-1
                y=np.roll(f.variables['dryTemperature'][0,:,:,:],72,axis=0)
                y[y>f.variables['dryTemperature']._FillValue/10.]=np.nan
                for ip in range(14):
                    x[ip,l,:,:]=y[:,:,ip].T
#		ganomalies=np.zeros([nj,ni,2,jpindex.shape[0],tmshape[4]],np.float32)

            ##    x[x==-9999.]=np.nan
                #x=x/100.
        s=x.shape
        wegc18=np.empty([1,s[0],s[1],s[2]//4,s[3]//4])
        rebin_72144_to_1836(x,wegc18,0,global_startyear=startyear,startyear=startyear,endyear=endyear)
        np.savez('wegc_'+str(startyear)+'_'+str(endyear)+'.npz', wegcfull=x,wegc18=wegc18)
    print(('read_wegc: ',time.time()-t1))
#    plt.contour(np.reshape(x[3,1300,:,:],[72,144]))
#    plt.show()
    #sys.exit()
    return x,wegc18

def read_wegc10(path,version,startyear=2006,endyear=2020,ens=0):

    t1=time.time()
    try:
        f=np.load('wegc_'+str(startyear)+'_'+str(endyear)+'.npz')
        x=f["wegcfull"]
        wegc18=f["wegc18"]
    except:
        l=0
        ps = np.array([10, 20, 30, 50, 70, 100, 150, 200, 250, 300, 400, 500, 700, 850, 925, 1000])*100.
        x=np.empty([16,(endyear-startyear+1)*12,18, 36])
        x[:]=np.nan
        for year in range(startyear,endyear+1):
            for month in range(12):
                mv='{:4}-{:0>2}-01'.format(year,month+1)
                mp1 = (month+2) % 13
                if mp1 == 0:
                    mp1 = 1
                mv1 ='{:4}-{:0>2}-01-'.format(year+(month+2)//13,mp1)
                fn=path+'/10x10-'+mv+'_'+mv1+ path.split('/')[-1][:-10]+version + '_L2b-data.nc'
                print(fn)
                try:
                    f=netCDF4.Dataset(fn)
#		    print f.variables.keys()
                except IOError:
                    try:
                        f=netCDF4.Dataset(fn[:-3]+'_interp'+fn[-3:])
                    except:
    #		    #print fn+' not found'
                        continue

                l=(year-startyear)*12+month-1
                y=np.roll(np.array(f.variables['dry_temperature'][0,:,:,:]),x.shape[-2],axis=0)
                y[y>f.variables['dry_temperature']._FillValue/10.]=np.nan
                gpsp = np.roll(np.array(f.variables['pressure'][0,:,:,:]),x.shape[-2],axis=0)
                gpsp[gpsp>f.variables['pressure']._FillValue/10.]=np.nan
                for ilat in range(x.shape[2]):
                    for ilon in range(x.shape[3]):
                        p = gpsp[ilat, ilon, :]
                        mask = ~np.isnan(p)
                        p = p[~np.isnan(p)][::-1]
                        t = y[ilat, ilon,mask][::-1]
                        plev = np.searchsorted(p, ps)
                        for ip in range(ps.shape[0]):
                            if plev[ip] >= p.shape[0] or plev[ip] == 0:
                                x[ip, l, ilat, ilon] = np.nan
                            else:
                                
                                x[ip,l,ilat,ilon]=(t[plev[ip]-1] * (p[plev[ip]] - ps[ip]) + \
                                    t[plev[ip]] * (ps[ip] - p[plev[ip]-1])) / \
                                    (p[plev[ip]] - p[plev[ip]-1])
                        
                        #if ilon == 0 and ilat == x.shape[2] // 3.5:
                            
                            #plt.semilogy(t, p)
                            #plt.semilogy(x[:, l, ilat, ilon], ps, 'o')
                            #plt.ylim(100000, 1)
                            #plt.show()
#		ganomalies=np.zeros([nj,ni,2,jpindex.shape[0],tmshape[4]],np.float32)

            ##    x[x==-9999.]=np.nan
                #x=x/100.
        #s=x.shape
        #wegc18=np.empty([1,s[0],s[1],s[2]//4,s[3]//4])
        #rebin_72144_to_1836(x,wegc18,0,global_startyear=startyear,startyear=startyear,endyear=endyear)
        wegc18 = x
        wegcfull = np.zeros((x.shape[0], x.shape[1], x.shape[2]*4, x.shape[3]*4))
        for ilat in range(x.shape[2]):
            for ilon in range(x.shape[3]):
                for i in range(4):
                    for j in range(4):
                        wegcfull[:, :, ilat*4+i,ilon*4+j ] = x[:, :, ilat, ilon]
        x = wegcfull
                        
        
        np.savez('wegc_'+str(startyear)+'_'+str(endyear)+'.npz', wegcfull=wegcfull,wegc18=wegc18)
    print(('read_wegc: ',time.time()-t1))
#    plt.contour(np.reshape(x[3,1300,:,:],[72,144]))
#    plt.show()
    #sys.exit()
    return x,wegc18

def make_corr(ipath,opath,statid):
    
    #for statid in ostnames: #[0:maxs0]:
    toofew=False
    flist=[]
    fn=statid
    #for d in tasks:
    if True:
        if not os.path.isfile(fn): # or '033345' not in fn:
            return


        f = netCDF4.Dataset(fn,"r")
        f.set_auto_mask(False)
        #if istat==0:
            #lats=np.empty(3000,np.float32)
            #lats[:]=np.nan
            #lons=lats.copy()
            #ps=f.variables['press'][:]

        try:
            dat=f.variables['datum'][0,:]
        except:
            return
        #print(istat,len(tasks[0]['data']))
        cstartdate=f.variables['datum'].getncattr('units').split()[2].split('-')
        startdate=datetime.date(int(cstartdate[0]),int(cstartdate[1]),int(cstartdate[2]))
        sd=19000101 #int(d['startdate'])
        offset=startdate-datetime.date(sd//10000,(sd-(sd//10000)*10000)//100,sd%100)
        dat=dat+offset.days
        dat[0]=1
        #istat+=1
        try:
            
            g=netCDF4.Dataset(opath+'/'+fn,"r")
        except:
            #print(opath+'/'+fn,'not found')
            return
        #igood+=1
        g.set_auto_mask(False)
        odat=g.variables['datum'][0,:]
        
        #fvarmax=f.variables['datum'][0,-1]
        #if fvarmax<23008 and fvarmax>gvmax:
            #gvmax=fvarmax
            #print(fn,gvmax)
        #continue
        
        cutdate=datetime.date(2010,1,1)-datetime.date(sd//10000,(sd-(sd//10000)*10000)//100,sd%100)
        cutdate=cutdate.days
        iidx=np.where(dat<cutdate)[0]
        oidx=np.where(odat>=cutdate)[0]
        ndat=np.empty_like(dat,shape=(len(iidx)+len(oidx)))
        nrdict={
            'rasocorr':np.empty_like(f.variables['rasocorr'][:],shape=(2,16,len(iidx)+len(oidx)))}
        if 'rasobreak' in f.variables.keys():
            
            nrdict['rasobreak']=np.empty_like(f.variables['rasobreak'][:],shape=(2,16,len(iidx)+len(oidx)))
        ndat[:len(iidx)]=dat[iidx]
        ndat[len(iidx):]=odat[oidx]
        
        for k,v  in nrdict.items():
            if len(oidx)>0:
                for l in range(v.shape[0]):
                    for m in range(v.shape[1]):
                        
                        v[l,m,iidx]=f.variables[k][l,m,iidx]+g.variables[k][l,m,oidx[0]-1]
                v[:,:,len(iidx):]=g.variables[k][:,:,oidx]
                print(fn,odat[oidx])
            else:
                v[:,:,iidx]=f.variables[k][:,:,iidx]
        
        #nrc[:,:,len(iidx):]=g.variables['rasocorr'][:,:,oidx]
        
        #cdict={}
        #corrs={'x':f.variables['rasocorr'][:],'y':g.variables['rasocorr'][:],'z':nrdict['rasocorr']}
        #dats={'x':dat,'y':odat,'z':ndat}
        #lw={'x':4,'y':2,'z':1}
        #for var in 'x','y','z':
            #cdict[var]=np.zeros((2,16,45000),np.float32)
            #for k in range(corrs[var].shape[2]-1,0,-1):
                #for l in range(corrs[var].shape[0]):
                    #for m in range(corrs[var].shape[1]):
                        #cdict[var][l,m,dats[var][k-1]:dats[var][k]]=corrs[var][l,m,k-1]
                
            #plt.plot(np.arange(45000)/365.25+1900,cdict[var][0,5,:],label=var,linewidth=lw[var])
        
        #plt.legend()
        ##plt.show()
        #plt.close()
        
        npath='ControlNew'.join(ipath.split('Control'))
        try:
            os.makedirs(npath+fn[-9:-3]+'/')
        except:
            pass
        with netCDF4.Dataset(npath+fn,'w') as h:
            ga=copy.copy(f.__dict__)
            ga['title']='Radiosonde temperature adjustments updated to 2019'
            ga['institution']='University of Vienna, Austria'
            ga['history']=datetime.datetime.now().strftime('%Y-%m-%d')
            ga['version']='1.5.1'
            h.setncatts(ga)
            for name,dimension in f.dimensions.items():
                if name=='time':
                    h.createDimension(name,ndat.shape[0])
                else:
                    h.createDimension(name,len(dimension))
            for name,variable in f.variables.items():
                x=h.createVariable(name,variable.datatype,variable.dimensions)
                if name in nrdict.keys():
                    h[name][:]=nrdict[name]
                elif name=='datum':
                    if f.variables[name].shape[0]==2:
                        h[name][:]=np.array((ndat,ndat))    
                    else:
                        h[name][:]=ndat   
                else:
                    h[name][:]=f[name][:]
                    
                ga=copy.copy(f[name].__dict__)
                if name=='datum':
                    ga['units']='1900'.join(ga['units'].split('1957'))
                h[name].setncatts(ga)
        
            #shutil.copy(opath+statid[-9:-3]+'/feedbackglobbincorrmon'+statid[-9:-3]+'.nc',npath+statid[-9:-3]+'/')
        if 'rio' not in statid and 'rit' not in statid:
            try:
                
                with netCDF4.Dataset(opath+statid[-9:-3]+'/feedbackglobbincorrmon'+statid[-9:-3]+'.nc','r') as f:
                    f.set_auto_mask(False)
                    with netCDF4.Dataset(npath+statid[-9:-3]+'/feedbackglobbincorrmon'+statid[-9:-3]+'.nc','w') as h:
                        h.set_auto_mask(False)
                        ga=copy.copy(f.__dict__)
                        ga['title']='Radiosonde temperature adjustments updated to 2019'
                        ga['institution']='University of Vienna, Austria'
                        ga['history']=datetime.datetime.now().strftime('%Y-%m-%d')
                        ga['version']='1.5.1'
                        h.setncatts(ga)
                        for name,dimension in f.dimensions.items():
                            h.createDimension(name,len(dimension))
                        for name,variable in f.variables.items():
                            if name =='rasocorrmon':
                                x=h.createVariable(name,variable.datatype,variable.dimensions)
                                mdays=f.variables['datum'][0,:]
                                nrdict[name]=copy.copy(variable[:])
                                for k in range(ndat.shape[0]):
                                    try:
                                        
                                        idx=np.where(mdays>ndat[k])[0][0]
                                        for i in range(nrdict[name].shape[0]):
                                            for j in range(nrdict[name].shape[1]):
                                                nrdict[name][i,j,idx:]=-nrdict['rasocorr'][i,j,k]
                                    except:
                                        pass
                                        
                                h[name][:]=nrdict[name]
                                idx=np.where(f.variables['montemp'][:]==-999.)
                                nrdict['montemp']=f.variables['montemp'][:]-f.variables['rasocorrmon'][:]+nrdict[name]
                                nrdict['montemp'][idx]=-999.
                                x=h.createVariable('montemp',f.variables['montemp'].datatype,f.variables['montemp'].dimensions)
                                h['montemp'][:]=nrdict['montemp']
                                
                                ga=copy.copy(f[name].__dict__)
                                h[name].setncatts(ga)
                                ga=copy.copy(f['montemp'].__dict__)
                                h['montemp'].setncatts(ga)
                            else:
                                if name !='montemp':
                                    x=h.createVariable(name,variable.datatype,variable.dimensions)
                                    
                                    h[name][:]=f[name][:]
                                
                                    ga=copy.copy(f[name].__dict__)
                                    h[name].setncatts(ga)
                        print('written')
                
            except Exception as e:
                print(e,'could not copy',opath+statid[-9:-3]+'/feedbackglobbincorrmon'+statid[-9:-3]+'.nc'+' to '+npath+statid[-9:-3]+'/')
        print('wrote',npath+fn)
                
        return


def append_stations(path,tasks,days,lats=[],lons=[],ps=[],stnames=[],minlen=24):

    sel='??'
    ipath=os.path.expandvars('$FSCRATCH/rise/1.0/Control/')
    os.chdir(os.path.expandvars('$FSCRATCH/rise/1.0/Control/'))
    #ostnames=glob.glob('00[1-9]???/feedbackglobbincorrsave_ri???_?'+sel+'???.nc')+glob.glob('0[1-9]????/feedbackglobbincorrsave_ri???_?'+sel+'???.nc')
    ostnames=glob.glob('00[1-9]???/feedbackglobbincorrsave?'+sel+'???.nc')+glob.glob('0[1-9]????/feedbackglobbincorrsave?'+sel+'???.nc')
    opath='exp00'.join(ipath.split('Control'))
    if len(stnames)==0:
        #stnames= glob.glob('00[1-9]???/feedbackglobbincorrsave_ri???_?'+sel+'???.nc')+glob.glob('0[1-9]????/feedbackglobbincorrsave_ri???_?'+sel+'???.nc')
        stnames= glob.glob('00[1-9]???/feedbackglobbincorrsave?'+sel+'???.nc')+glob.glob('0[1-9]????/feedbackglobbincorrsave?'+sel+'???.nc')
    pindex=np.arange(16,dtype=np.int32)

    tt=time.time()
    istat=0
    goodsts=[]
    found_tmcorr=False
    igood=0
    gvmax=0
    func=partial(make_corr,ipath,opath)
    p=Pool(30)
    list(map(func,ostnames))
    


    for k in range(len(tasks)):
        try:
            tasks[k]["data"]=np.asarray(tasks[k]["data"]) #[0:istat,:,:,:,:])
            tasks[k]["msudata"]=np.asarray(tasks[k]["msudata"]) #[0:istat,:,:,:,:]
        except:
            pass

    print((time.time()-tt))
    lats=lats[0:istat]
    lons=lons[0:istat]

    return istat,lats,lons,ps,goodsts

@njit(cache=True)
def monthly_average(mav,hilf,mstarts,thresh):
    
    
    for iens in range(hilf.shape[1]):
        for ih in range(hilf.shape[2]):
            for ip in range(hilf.shape[3]):
                for im in range(mstarts.shape[0]-1):
                    mv=0.
                    ic=0
                    for iday in range(mstarts[im],mstarts[im+1]):
                        h=hilf[0,iens,ih,ip,iday]
                        if h==h:
                            mv+=h
                            ic+=1
                    if ic>thresh:
                        mav[iens,ih,ip,im]=mv/ic
                    else:
                        mav[iens,ih,ip,im]=np.nan
    return
                        
def mon_mean(data,days,idays):
    montemp=[]
    good=[]
    gdays=[]
    astart,astop=np.searchsorted(idays+1,(days[0],days[-1]))
    for i in range(astart,astop):
        start,stop=np.searchsorted(days,idays[i:i+2]+1)
        if stop>start:
            d=data[:,:,start:stop]
            x=np.sum(~np.isnan(d),axis=2).reshape((data.shape[0],data.shape[1],1))
            if np.sum(x)>0:
                
                good.append(x.reshape((data.shape[0],data.shape[1],1)))
                montemp.append(np.nanmean(d,axis=2).reshape((data.shape[0],data.shape[1],1)))
                gdays.append(i)
    
    return np.concatenate(montemp,axis=2),np.concatenate(good,axis=2),np.array(gdays)

def readsuny(d,tasks,fn,dc,istat,statid,ens,pindex,ymds,idays,minlen,lats=None,lons=None,ps=None,stlongnames=None,monthly=False):
    tt=time.time()
    try:
        spress=np.array((10.,20,30,50,70,100,150,200,250,300,400,500,700,850,925,1000))
        with h5py.File(fn) as h:
            station_name=''
            df_press=h['pressure'][:]
            idx=np.searchsorted(df_press,spress*100)
            idx[idx==df_press.shape[0]]=df_press.shape[0]-1
            
            #refdate=datetime.datetime(1900,1,1)
            #hty=h['time'][:]//10000
            #htm=(h['time'][:]%10000)//100
            #htd=h['time'][:]%100
            #df_time=[datetime.datetime(hty[i],htm[i],htd[i]) for i in range(hty.shape[0])]
            #df_days=np.array([(df_time[i]-refdate).days for i in range(hty.shape[0])])+1
            df_days=np.searchsorted(ymds,h['time'][:])+1
            #df_time=pd.to_datetime(df_days-1,unit='d',origin='19000101').values
            
            mask=h['rawT'][:]==-9999.
            if d['shortname']=='suny':
                
                y=h['rawT'][:]*0.1
            else:
                
                y=h['homoT'][:]*0.1
                
            y[mask]=np.nan
            y=np.einsum('kli->ilk', y)
            x=np.empty_like(y,shape=(y.shape[0],spress.shape[0],y.shape[2]))
            x=y[:,idx,:]
            for i in range(len(idx)):
                if not any(df_press==spress[i]*100):
                    x[:,i,:]=np.nan
            
            df_press=spress
            x[x>400]=np.nan

            print(time.time()-tt)
            montemp,good,gdays=mon_mean(x,df_days,idays)
            print(time.time()-tt)
            #with h5py.File( '/users/staff/leo/fastscratch/SUNY/UA-HRD_stations_homo+raw-monthly-anomalyT_00-12Z_195801-202008.nc') as g:
                
                #print('x')
                
            
            if istat==0:
                sh=d['data'][1:]
                msh=d['msudata'][1:]
                d['data']=[np.empty(sh)]
                d['msudata']=[np.empty(msh)]
            else:
                d['data'].append(np.empty_like(d['data'][-1]))
                d['msudata'].append(np.empty_like(d['msudata'][-1]))
            d['data'][-1].fill(np.nan)
            d['msudata'][-1].fill(np.nan)
            for ih in 0,1:
                for ip in pindex:
                    d['data'][-1][0,ih,ip,gdays]=montemp[ih,ip,:]
            
            
                
    except MemoryError:
        return False,0,spress
        #with h5py.File(tup[0]) as o:
            
            #df_lat=o['lat'][:]
            #df_lon=o['lon'][:]
        #mask=np.zeros((2,3,df_time.shape[0]),dtype=bool)
    print(time.time()-tt)
    return True,1,spress
                   

def read_alltrends(path,tasks,days,satlist=[],lats=[],lons=[],ps=[],stnames=[],minlen=24,refseries='tmcorr',reftasks=[],plotproperties=None):

    # Activate if fill_suny is not used
    #with netCDF4.Dataset( '/users/staff/leo/fastscratch/SUNY/UA-HRD_stations_homo+raw-monthly-anomalyT_00-12Z_195801-202008.nc') as g:
        #print('x')
        #igralats=g['lat'][:]
        #igralons=g['lon'][:]
        #igraIDs=np.ascontiguousarray(g['SiteID'][:].T).view('|S11')
        #igrasIDs=np.ascontiguousarray(g['SiteID'][5:,:].T).view('|S6')

    
    if len(stnames)==0:
        stnames0= glob.glob('[0-9]?????')
        stnames=[]
        for s in stnames0:
            if '-' not in s:
                stnames.append(s)
        stnames.sort()
    #stnames=np.array(['068842'])
    pindex=np.arange(16,dtype=np.int32)
#    stnames=['001001','011035']
#    str=[' ']

    tt=time.time()
    alldays=pd.date_range(start='1900-01-01',end='2023-02-01',freq='d')
    ymds=alldays.year.values*10000+alldays.month.values*100+alldays.day.values
    sdays=pd.date_range(start='1900-01-01',end='2023-02-01',freq='MS')
    idays=np.array([(sdays[i]-sdays[0]).days for i in range(len(sdays))])
    istat=0
    goodsts=[]
    found_tmcorr=False
    for d in tasks:
        if d['shortname']==refseries or len(tasks) == 1:
            found_tmcorr=True

    if len(tasks) == 1:
        refseries = tasks[0]['shortname']
    l=-1
    for statid in stnames: #[0:maxs0]:
        
        if statid == '001010':
            x =0
        l+=1
        found=not found_tmcorr
        toofew=False
        flist=[]
        for d in tasks:
            flist.append(False)
            #if len(d["file"])>0 and toofew==False:
            if  toofew==False:
                for ens in d["ens"]:
                    try:
                        if type(d['relpath']) is list:
                            fpath=d['relpath'][0]
                        else:
                            fpath=d['relpath']
                        sids=0
                    except:
                        fpath=path
                        sids=0
                    if len(d["ens"])==1:
                        fn=os.path.join(fpath,statid[sids:],d["file"]+statid[sids:]+d["suff"]+".nc")
                        if "andepmon" in str(d["file"]):
                            #gn=''
                            gn=os.path.join(fpath,statid[sids:],d["msufile"]+d["suff"]+'_bt2_'+statid[sids:]+".nc")
                        elif 'SUNY' in d['file']:
                            # you should use this only if fill_suny is not used. 
                            try:
                                idx=np.where(igrasIDs.flatten()==np.string_(statid))
                                fn=d['file']+'/homo-raw-subdaily-station/'+igraIDs[idx][0][0].decode()+d["suff"]+".nc"
                                gn=os.path.join(fpath,statid[sids:],d["msufile"]+d["suff"]+'_bt2_'+statid+".nc")
                            except:
                                print('not found, looking for matching lat/lon')
                                idx=np.where(np.logical_and(np.abs(igralats-lats[l])<2.0,np.abs(igralons-lons[l])<2.0))[0]
                                if len(idx)>0:
                                    fn=glob.glob(d['file']+'/homo-raw-subdaily-station/*'+igraIDs[idx[0]][0].decode().strip()+d["suff"]+".nc")[0]
                                    gn=os.path.join(fpath,statid[sids:],d["msufile"]+d["suff"]+'_bt2_'+statid+".nc")
                                    print('found matching latitude')
                        else:
                            gn=os.path.join(fpath,statid[sids:],d["msufile"]+d["suff"]+'_bt2_'+statid[sids:]+".nc")
                    else:
                        if "ce20c" in d["file"]:
                            fn=os.path.join(fpath,statid[sids:],d["file"]+"{0:0>1}".format(ens)+d["suff"]+statid[sids:]+".nc")
                            gn='' #os.path.join(fpath,statid[sids:],d["msufile"]+"{0:0>1}".format(ens)+d["suff"]+'bt2'+statid[sids:]+".nc")
                        else:
                            fn=os.path.join(fpath,statid[sids:],d["file"]+d["suff"]+"{0:0>2}".format(ens)+'_'+statid[sids:]+".nc")
                            gn=os.path.join(fpath,statid[sids:],d["msufile"]+d["suff"]+"{0:0>2}".format(ens)+'_bt2_'+statid[sids:]+".nc")

#                    print d["shortname"]+' '+fn
                    if statid=='085201':
                        print (fn)
                    print((d['shortname'],fn))
                    if not os.path.isfile(fn): # or '033345' not in fn:
                        break


                    
                    if 'suny' in d['shortname']:
                        # this code is needed only if fill_suny is not used
                        dc=0
                        ifound,dc,press=readsuny(d,tasks,fn,dc,istat,statid,ens,pindex,ymds,idays,minlen,lats=lats[l],lons=lons[l],ps=ps,stlongnames=statid,monthly=True)
                        if os.path.isfile(gn):
                            g = netCDF4.Dataset(gn,"r")
                            dat=g.variables['datum'][0,:]
                            index=np.zeros(dat.shape[0],dtype=np.int32)
                            getindex(days,dat,index) 
                            v=g.variables[d['msuvar']][:]
                            if v.shape[1]==3:
                                w=np.empty((v.shape[0],4,v.shape[2]))
                                w[:,1:,:]=v
                                w[:,0,:]=np.nan
                                v=w
                            #copystride4(d["msudata"],v,index,istat,ens,
                                        #np.arange(v.shape[1]),nc_miss_val)
                            copystride4(d["msudata"][istat],v,index,0,ens,
                                        np.arange(v.shape[1]),np.nan)
                        else:
                            d["msudata"][istat][:,:,:,:]=np.nan
                        
                        flist[-1]=True
                        continue
                        

#                    t3=time.time()
                    f = netCDF4.Dataset(fn,"r")
                    f.set_auto_mask(False)
                    sat=False
                    try:
                        if gn != '':
                            g = netCDF4.Dataset(gn,"r")
                            g.set_auto_mask(False)
                            sat=True
                        else:
                            sat=False
                    except:
                        print(('no' +gn))
                        sat=False
#		    else:
#                    except:
#                        print fn+' not found'
#                        continue

                    if istat==0:
                        if d["shortname"]==refseries:
                            lats=np.empty(3000,np.float32)
                            lats[:]=np.nan
                            lons=lats.copy()
                            ps=f.variables['press'][:]
                            if f.variables['press'].getncattr('units')=='Pa':
                                ps=np.array(ps)/100
                        if ens==0:
                            for da in ('data','msudata'):
                                sh=[]
                                for dl in range(1,len(d[da])):
                                    sh.append(d[da][dl])
                                d[da]=[]
                                d[da].append(np.empty(sh,dtype=np.float32))
                                d[da][0].fill(np.nan)
                        d['sonde_type'] = []
                        #if len(days)==0:
                            #days=f.variables['datum'][0,:]
#		    print time.time()-t3
                    try:
                        dat=f.variables['datum'][0,:]
                        if dat.dtype==np.float32:
                            dat=np.array(dat,dtype=np.int64)
                            if 'seconds' in f.variables['datum'].getncattr('units'):
                                dat=dat//86400
                        
                        if 'sonde_type' in f.variables.keys():
                            sonde_type =f.variables['sonde_type'][:]
                        else:
                            print(d['shortname'],f.filepath().split('/')[-1], 'no sonde_type, try supplying from task 0')
                            sonde_type = None
                    except:
                        continue
                    print(istat,len(tasks[0]['data']))
                    cstartdate=f.variables['datum'].getncattr('units').split()[2].split('-')
                    startdate=datetime.date(int(cstartdate[0]),int(cstartdate[1]),int(cstartdate[2]))
                    sd=int(d['startdate'])
                    offset=startdate-datetime.date(sd//10000,(sd-(sd//10000)*10000)//100,sd%100)
                    dat=dat+offset.days
                    if d["shortname"]==refseries:
                        lats[istat]=f.variables['lat'][:]
                        lons[istat]=f.variables['lon'][:]
                        if abs(lats[istat])>90.0:
                            print(('Spurious lat',istat,lats[istat]))
                            break

                    nc_miss_val=np.float32(-999.)
                    if istat==10:
                        print(fn)

                    if len(d["index"])==1:
                        if dat.shape[0]<minlen:
                            if d["shortname"]==refseries:
                                toofew=True
                            print((statid+' series too short'))
                            f.close()
                            if sat:
                                g.close()
                            break
                        index=np.zeros(dat.shape[0],dtype=np.int32)
                        getindex(days,dat,index)
                        dv=d['var']
                        if dv not in list(f.variables.keys()):
                            for dk in list(f.variables.keys()):
                                if dv in dk:
                                    dv=dk
                                    break
                        #copystride4(d["data"],f.variables[dv][:],index,istat,ens,pindex,nc_miss_val)
                        if istat>0:
                            for da in ('data','msudata'):
                                d[da].append(np.empty_like(d[da][-1]))
                                d[da][-1].fill(np.nan)
                            #d['sonde_type'].append(numpy.zeros(d['data'][0].shape[3], dtype='|S4'))
                        if d['shortname']=='hadat':
                            d['data'][istat][0,:,0,:]=np.nan
                            copystride4(d["data"][istat],f.variables[dv][:],index,0,ens,np.asarray([15]),nc_miss_val)
                            for p in pindex:
                                d['data'][istat][0,:,p,:]=d['data'][istat][0,:,0,:]

                        else:
                            #d['data'][istat,:,:,:,:]=np.nan
                            if '_rms' in d['shortname']:
                                if "hilf" not in locals():
                                    hilf=np.empty_like(d["data"][istat])

                                hilf[:]=np.nan
                                dvm=dv.split('_std')[0]
                                copystride4(hilf,f.variables[dvm][:],index,0,ens,pindex,nc_miss_val)
                                copystride4(d["data"][istat],f.variables[dv][:],index,0,ens,pindex,nc_miss_val)
                                d["data"][istat][:]=np.sqrt(d["data"][istat][:]*d["data"][istat][:]+hilf[istat][:]*hilf[istat][:])
                            elif 'era5v' in d['shortname']:
                                x=np.empty((2,45000),dtype=np.float32)
                                for i in range(pindex.shape[0]):
                                    x.fill(np.nan)
                                    x[:,dat]=f.variables[dv][:,pindex[i],:]
                                    x[x==nc_miss_val]=np.nan
                                    d["data"][istat][0,:,pindex[i],:]=monmean(x,days,thresh=3)
                            else:
                                if 'seconds' in f.variables['datum'].getncattr('units'):
                                    s=d["data"][0].shape
                                    hilf=np.empty((1,s[0],s[1],s[2],d['ddata'][-1]))
                                    mstarts=calcdays(19000101,d['data'][istat].shape[3])
                                    hilf.fill(np.nan)
                                    index=dat
                                    copystride(hilf,np.array(f.variables[dv][:]),index,0,ens,pindex,nc_miss_val)
                                    monthly_average(d["data"][istat],hilf,mstarts,15)
                                else:
                                    try:
                                        plotproperties['monthresh']=5
                                        good=f.variables['goodmon'][:]
                                        mask=good<=plotproperties['monthresh']
                                        dum=f.variables[dv][:]
                                        dum[mask]=np.nan
                                    except:                                       
                                        dum=f.variables[dv][:]
                                    copystride4(d["data"][istat],dum,index,0,ens,pindex,nc_miss_val)
                                    d['sonde_type'].append(numpy.zeros(d['data'][0].shape[3], dtype='|S4'))
                                    d['sonde_type'][istat][index] = sonde_type


                        if sat:
                            if statid=='068842':
                                print(statid)
                            dat=g.variables['datum'][0,:]
                            index=np.zeros(dat.shape[0],dtype=np.int32)
                            getindex(days,dat,index)
                            try:
                                v=g.variables[d['msuvar']][:]
                                if v.shape[1]==3:
                                    w=np.empty((v.shape[0],4,v.shape[2]))
                                    w[:,1:,:]=v
                                    #w=np.empty((1,4,v.shape[2]))
                                    #w[:,1:,:]=np.nanmean(v,axis=0)
                                    w[:,0,:]=nc_miss_val
                                    v=w
                                v[v<150.]=np.nan
                                v[v>320.]=np.nan
                                copystride4(d["msudata"][istat],v,index,0,ens,
                                            np.arange(v.shape[1]),nc_miss_val)
                            except:
                                d["msudata"][istat][:,:,:,:]=np.nan
                               
                        else:
                            d["msudata"][istat][:,:,:,:]=np.nan
#			    if np.abs(np.nanmean(d["msudata"]))>1000.:
#				print np.nanmean(d["msudata"])

                    else:
                        index=np.zeros(dat.shape[0],dtype=np.int32)
                        getindex2(days,dat,index)   
                        d["index"][istat,ens,0:dat.shape[0]]=index
                        index=np.arange(dat.shape[0],dtype=np.int32)
                        #if type(d['data']) is list:
                            #d['data']=np.empty(d['data'],dtype=np.float32)
                        #if type(d['msudata']) is list:
                            #d['msudata']=np.empty(d['msudata'],dtype=np.float32)
                        if istat>0 and ens==0:
                            for da in ('data','msudata'):
                                if type(d[da][0]) is int:
                                    sh=[]
                                    for dl in range(1,len(d[da])):
                                        sh.append(d[da][dl])
                                    d[da].append(np.empty(sh,dtype=np.float32))   
                                    d[da][-1].fill(np.nan)
                                else:
                                    d[da].append(np.empty_like(d[da][-1]))
                                    d[da][-1].fill(np.nan)


                        d["data"][istat][ens,:,:,0:dat.shape[0]]=f.variables[d['var']][:]
                        try:
                            
                            d['sonde_type'].append(tasks[0]['sonde_type'][istat])
                        except:
                            corrmonpath = f.filepath().split('/')
                            corrmonpath[-1] = 'feedbackglobbincorrmon' + corrmonpath[-1][-9:]
                            corrmonpath = '/'.join(corrmonpath)
                            try:
                                
                                with Dataset(corrmonpath) as fcorr:
                                    cdat = np.array(fcorr.variables['datum'][0, :])
                                    cindex=np.zeros(cdat.shape[0],dtype=np.int32)
                                    getindex(days,cdat,cindex)
                                    d['sonde_type'].append(numpy.zeros(days.shape[0], dtype='|S4'))
                                    d['sonde_type'][istat][cindex] = np.array(fcorr.variables['sonde_type'][:])
                            except:
                                print('no metadata')
                                
                        if sat:
                            if statid=='068842':
                                print(statid)
                            dat=g.variables['datum'][0,:]
                            index=np.zeros(dat.shape[0],dtype=np.int32)
                            getindex(days,dat,index) 
                            v=g.variables[d['msuvar']][:]
                            if v.shape[1]==3:
                                w=np.empty((v.shape[0],4,v.shape[2]))
                                w[:,1:,:]=v
                                w[:,0,:]=nc_miss_val
                                v=w
                            #copystride4(d["msudata"],v,index,istat,ens,
                                        #np.arange(v.shape[1]),nc_miss_val)
                            copystride4(d["msudata"][istat],v,index,0,ens,
                                        np.arange(v.shape[1]),nc_miss_val)
                        else:
                            d["msudata"][istat][:,:,:,:]=np.nan

                            #if np.abs(np.nanmean(d["msudata"]))>1000.:
#				print np.nanmean(d["msudata"])

                    f.close()
                    if sat:
                        g.close()

                    flist[-1]=True
                    if d["shortname"]==refseries:
                        found=True
            else:
                print (d["name"]+' no data')
            if not flist[-1] and d["shortname"]==refseries:
                break
        #if len(tasks[0]['sonde_type'])>len(tasks[1]['sonde_type']):
            #print('one is missing')
        if found or len(tasks) == 1:
            for k in range(len(tasks)):
                if not flist[k] and tasks[k]['shortname'] not in satlist:
                    if istat==0:
                        sh=tasks[k]['data']
                        tasks[k]['data']=[]
                        tasks[k]['data'].append(np.empty(sh[1:],dtype=np.float32))
                        sh=tasks[k]['msudata']
                        tasks[k]['msudata']=[]
                        tasks[k]['msudata'].append(np.empty(sh[1:],dtype=np.float32))
                        tasks[k]['sonde_type'] = []
                        tasks[k]['sonde_type'].append(tasks[0]['sonde_type'][istat])
                    else:	
                        tasks[k]['data'].append(np.empty_like(tasks[k]['data'][-1]))
                        tasks[k]['msudata'].append(np.empty_like(tasks[k]['msudata'][-1]))
                        if istat < len(tasks[0]['sonde_type']):
                            tasks[k]['sonde_type'].append(tasks[0]['sonde_type'][istat])
                        else:
                            tasks[k]['sonde_type'].append(np.zeros_like(tasks[0]['sonde_type'][0]))
                            
                    tasks[k]['data'][istat][:,:,:,:]=np.nan
                    tasks[k]['msudata'][istat][:,:,:,:]=np.nan
            goodsts.append(statid)
#            print statid,istat,lons[istat]
            istat+=1

    for k in range(len(tasks)):
        #if not flist[k] and tasks[k]['shortname'] not in satlist:
            try:
                tasks[k]["data"]=np.asarray(tasks[k]["data"]) #[0:istat,:,:,:,:])
                tasks[k]["msudata"]=np.asarray(tasks[k]["msudata"]) #[0:istat,:,:,:,:]
                tasks[k]["sonde_type"]=np.asarray(tasks[k]["sonde_type"]) #[0:istat,:,:,:,:]
            except:
                pass

    print((time.time()-tt))
    lats=lats[0:istat]
    lons=lons[0:istat]

    return istat,lats,lons,ps,goodsts

@njit(parallel=False,cache=True)
def satst(lats,lons,btg,mmax):
    lati=((lats+90.)//2.5).astype(np.int32)
    loni=(lons//2.5).astype(np.int32)
    bts=np.empty((lats.shape[0],1,2,btg.shape[0],mmax),dtype=np.float32)
    for istat in range(lats.shape[0]):
        psatst(lati,loni,btg,bts,istat)

    return bts

@njit()
def psatst(lati,loni,btg,bts,istat):

    timin=min((bts.shape[4],btg.shape[1]))
    for ichan in range(btg.shape[0]):
#	for istat in range(lats.shape[0]):
        for it in range(timin):
            bts[istat,0,0,ichan,it]=btg[ichan,it,lati[istat],loni[istat]]
            bts[istat,0,1,ichan,it]=bts[istat,0,0,ichan,it]
        if btg.shape[1]<bts.shape[4]:
            bts[istat,0,:,ichan,btg.shape[1]:]=np.nan
    return

def loadsondes(path,tasks,days,plotproperties,stnames,init,slowload,satlist):
    t=time.time()

    newall=False
    try:
        npztime=os.path.getmtime(plotproperties['tmppath']+'/'+'allsave.npz')
    except:
        npztime=0
    for d in tasks:
        try:
            if slowload:
                fn=plotproperties['tmppath']+'/'+d["shortname"]+'.npz'
                npztime2=os.path.getmtime(fn)
            else:
                if len(d['ens'])<2:
                    vn='_data_{0:0>2}'.format(plotproperties['pindex'][0])
                else:
                    vn='_data_'+plotproperties['ens'][0]+'_'+'{0:0>2}'.format(plotproperties['pindex'][0])

                fn=plotproperties['tmppath']+'/'+d["shortname"]+vn+'.npz'   
                npztime2=os.path.getmtime(fn)
            raobcoretime=os.path.getmtime(stnames[0]+'/'+d['file']+stnames[0]+'.nc')
            if raobcoretime>npztime2:
                os.remove(fn)
            if d['shortname']=='tmcorr' and raobcoretime>npztime:
                try:
                    os.remove(plotproperties['tmppath']+'/'+'allsave.npz')
                except:
                    pass
        except Exception as e:
            print(e)
            pass
    try:
        k=0
        llist=[]
        d=np.load(plotproperties['tmppath']+'/'+'allsave.npz')
        lats=d["lats"]
        lons=d["lons"]
        days=d["days"]
        ps=d["ps"]
        stnames=d["stnames"]
        d.close()
        istat=lats.shape[0]
        if init:
            for k in range(len(tasks)):
                try:
                    if slowload:
                        vn=''
                        d=np.load(plotproperties['tmppath']+'/'+tasks[k]["shortname"]+'.npz')
                        tasks[k]["data"]=d["data"].astype(np.float32)
                        tasks[k]["index"]=d["index"]
                        tasks[k]["msudata"]=d["msudata"].astype(np.float32)
                        try:
                            
                            tasks[k]["sonde_type"]=d["sonde_type"]
                        except:
                            tasks[k]["sonde_type"]=np.zeros(tasks[k]['data'].shape[-1], dtype='S4')
                            
                        d.close()
                    else:
                        if len(tasks[k]['ens'])<2:
                            vn='_data_{0:0>2}'.format(plotproperties['pindex'][0])
                        else:
                            vn='_data_'+plotproperties['ens'][0]+'_'+'{0:0>2}'.format(plotproperties['pindex'][0])

                        d=np.load(plotproperties['tmppath']+'/'+tasks[k]["shortname"]+vn+'.npz')
                        sh=d['arr_0'].shape
                        tasks[k]["data"]=np.reshape(d["arr_0"],(sh[0],1,sh[1],1,sh[2]))
                        d.close()
                        d=np.load(plotproperties['tmppath']+'/'+tasks[k]["shortname"]+'_index.npz')
                        sh=d['arr_0'].shape
                        tasks[k]["index"]=d['arr_0']
                        d.close()
                        d=np.load(plotproperties['tmppath']+'/'+tasks[k]["shortname"]+'_sonde_type.npz')
                        sh=d['arr_0'].shape
                        tasks[k]["sonde_type"]=d['arr_0']
                        d.close()
                        print(('read: ',tasks[k]["shortname"],time.time()-t))
                except Exception as e:
                    if tasks[k]['shortname'] not in ['suny','sunyhom']:
                        llist+=[tasks[k]]
                    print((plotproperties['tmppath']+'/'+tasks[k]["shortname"]+vn+'.npz not found or outdated'))
                    print(('need to read from scratch: ',tasks[k]["shortname"],time.time()-t))

    except: # Exception as e:

        traceback.print_exc()
        refseries='tmcorr'
        for t in tasks:
            if t['shortname']=='uwind':
                refseries='uwind'
        istat,lats,lons,ps,stnames=read_alltrends(path,tasks,days,refseries=refseries,plotproperties=plotproperties)
        t=time.time()
        mktmppath(plotproperties['tmppath'])
        np.savez(plotproperties['tmppath']+'/'+'allsave.npz',
                    lats=lats,
                    lons=lons,
                    days=days,
                    ps=ps,
                    stnames=stnames

                    )
        newall=True

    if slowload:
        if len(llist)>0:
            istat,lats,lons,ps,stnames=read_alltrends(path,llist,days,satlist,lats,lons,ps,stnames)
        if newall:
            llist=copy.deepcopy(tasks)
        for k in range(len(llist)):
            if 'sonde_type' not in llist[k].keys():
                llist[k]["sonde_type"] = tasks[0]['sonde_type']
            np.savez(plotproperties['tmppath']+'/'+llist[k]["shortname"]+'.npz',
                        data=llist[k]["data"],
                        index=llist[k]["index"],
                        msudata=llist[k]["msudata"], 
                        sonde_type=llist[k]["sonde_type"]
                        )
            print(('save: ',llist[k]["shortname"],time.time()-t))
    else:
        if len(llist)>0:
            istat,lats,lons,ps,stnames=read_alltrends(path,llist,days,satlist,lats,lons,ps,stnames,reftasks=tasks, plotproperties=plotproperties)
        if newall:
            llist=copy.deepcopy(tasks)
        for k in range(len(llist)):
            sd=dict()
            for iens in llist[k]['ens']:
                for ip in range(ps.shape[0]):
                    if len(llist[k]['ens'])<2:
                        vn='_data_{0:0>2}'.format(ip)
                    else:
                        vn='_data_{0:0>2}_{1:0>2}'.format(iens,ip)
                    #sh=llist[k]["data"].shape
                    sd[vn]=llist[k]["data"][:,iens,:,ip,:]
                    np.savez(plotproperties['tmppath']+'/'+llist[k]["shortname"]+vn+'.npz',sd[vn])
            sd['index']=llist[k]["index"]
            np.savez(plotproperties['tmppath']+'/'+llist[k]["shortname"]+'_index.npz',sd['index'])
            np.savez(plotproperties['tmppath']+'/'+llist[k]["shortname"]+'_sonde_type.npz',llist[k]['sonde_type'])
            print(('save: ',llist[k]["shortname"],time.time()-t))

    return istat,lats,lons,ps,stnames

def read_mesural(fn):
    with open(fn) as f:
        rdata=f.read().split('\n')
    mnames=[]
    for r in rdata:
        if len(r)>0:
            if 'SONDES' in r:
                sn=r.split('SONDES')[1].split('MESURAL')[0]
                try:
                    if int(sn)<100000:
                        if '00'+sn[:-1].strip() not in mnames:
                            mnames.append('00'+sn[:-1].strip())
                    else:
                        if '0'+sn[:-1].strip() not in mnames:
                            mnames.append('0'+sn[:-1].strip())
                except:
                    #print r,sn
                    sn=r.split('SONDES')[1].split('METOX-')[0]
                    try:
                        if int(sn)<100000:
                            if '00'+sn[:-1].strip() not in mnames:
                                mnames.append('00'+sn[:-1].strip())
                        else:
                            if '0'+sn[:-1].strip() not in mnames:
                                mnames.append('0'+sn[:-1].strip())
                    except:
                        print((r,sn))
    return mnames

def read_viz(fn):
    with open(fn) as f:
        rdata=f.read().split('\n')
    mnames=[]
    oldsn=''
    viz=False
    for r in rdata:
        if len(r)>0:
            if 'SONDES' in r:
                sn=r.split('SONDES')[1][:8]
#		if not sn.isdigit():
#		    sn=sn[:5]
                try:
                    sn='{:0>6}'.format(int(sn))
                    if oldsn!=sn:
                        if viz:
                            mnames.append(oldsn)
                        oldsn=sn
                        viz=False
                    if 'VIZ' in r:
                        year=int(r[91:95])
                        if year<1988:
                            viz=True
                    else:
                        year=int(r[91:95])
                        if year<1988:
                            viz=False

                except:
                    print((r,sn))
    return mnames
#   18     0SONDES    10010VAISALA-RS80                         VAISALA-RS80          938   198600009912                 1984010199          1

def fill_suny(tasks,shortnames,stnames,lats,lons,ps,pindex):
    ishs=[]
    for s in ['suny','sunyhom']:
        try:
            ishs.append(shortnames.index(s))
            if type(tasks[ishs[-1]]['data']) is list:
                tasks[ishs[-1]]['data']=np.empty([tasks[0]['data'].shape[0]]+tasks[ishs[-1]]['data'][1:],dtype=np.float32)
                tasks[ishs[-1]]['msudata']=np.empty([tasks[0]['data'].shape[0]]+tasks[ishs[-1]]['msudata'][1:],dtype=np.float32)
                tasks[ishs[-1]]['data'].fill(np.nan)
                tasks[ishs[-1]]['msudata'].fill(np.nan)
        except:
            pass
        
    ifname=tasks[ishs[0]]['file']+'UA-HRD_stations_homo+raw-monthly-anomalyT_00-12Z_195801-202008.nc'
    with netCDF4.Dataset(ifname,'r') as f:
        siteids=np.asarray(f.variables['SiteID'][:].T.flatten().view('S11').astype('U11'))
        slats=np.array(f.variables['lat'])
        slons=np.array(f.variables['lon'])
        ipx=np.arange(16,dtype='int')
        ipx[15]=14
        itxy=(f.variables['time'][:]-190000)//100
        itym=(f.variables['time'][:]-190000)%100-1
        itx=np.array(itxy*12+itym)
        itxfinal=np.searchsorted(itx,tasks[ishs[0]]['data'].shape[4])
        for ish in ishs:
            
            idx=np.zeros(siteids.shape[0],dtype=int)-9999
            if shortnames[ish]=='suny':
                T=np.array(f.variables['rawT'][:])
            elif shortnames[ish]=='sunyhom':
                T=np.array(f.variables['homoT'][:])
            else:
                continue
            T[T==-9999.]=np.nan
            j=0
            jj=0
            for i in range(siteids.shape[0]):  
                
                try:
                    idx[i]=np.where(stnames==siteids[i].rstrip()[-6:])[0][0]
                    #print(idx[i],stnames[idx[i]],siteids[i].rstrip()[-6:])
                    j+=1
                except:
                    continue
                    #idx[i]=tasks[ish]['data'].shape[0]

                if False and idx[i]==tasks[ish]['data'].shape[0]:
                    #print(i,siteids[i],idx[i])
                    #print('missed ID, trying lat/lon')
                    #print(slats[i],slons[i])
                    #print(np.where(np.logical_and(np.abs(slats[i]-lats)<1.0,np.abs(slons[i]-lons)<1.0)))
                    try:
                        
                        idx[i]=np.where(np.logical_and(np.abs(slats[i]-lats)<2.0,np.abs(slons[i]-lons)<2.0))[0][0]
                        j+=1
                        if '0912' in stnames[idx[i]] or '0911' in stnames[idx[i]]:
                            
                            print('now found',i,idx[i],siteids[i],slats[i],lats[idx[i]],slons[i],lons[idx[i]])
                    except:
                        jj+=1
                        continue
                    
                for ih in range(2):
                    for ip in range(len(pindex)):
                        
                        tasks[ish]['data'][idx[i],0,ih,pindex[ip],itx[0]:itx[0]+itxfinal]=T[i,:itxfinal,ipx[pindex[ip]],ih]

            for i in range(siteids.shape[0]):  
                
                if idx[i]==-9999:
                    #print(i,siteids[i],idx[i])
                    #print('missed ID, trying lat/lon')
                    #print(slats[i],slons[i])
                    #print(np.where(np.logical_and(np.abs(slats[i]-lats)<1.0,np.abs(slons[i]-lons)<1.0)))
                    try:
                        
                        id=np.where(np.logical_and(np.abs(slats[i]-lats)<2.0,np.abs(slons[i]-lons)<2.0))[0]
                        for ix in id:
                            if ix not in idx:
                                idx[i]=ix
                                break
                        if idx[i]==-9999:
                            print('no match',siteids[i])
                            continue
                        j+=1
                        #if '0912' in stnames[idx[i]] or '0911' in stnames[idx[i]]:
                            
                            #print('now found',i,idx[i],siteids[i],slats[i],lats[idx[i]],slons[i],lons[idx[i]])
                    except:
                        jj+=1
                        continue
                    
                for ih in range(2):
                    for ip in range(len(pindex)):
                        
                        tasks[ish]['data'][idx[i],0,ih,pindex[ip],itx[0]:itx[0]+itxfinal]=T[i,:itxfinal,ipx[pindex[ip]],ih]
           
    np.savetxt('suny.idx',idx,fmt='%5d')    
    sdays=pd.date_range(start='1900-01-01',end='2023-02-01',freq='MS')
    idays=np.array([(sdays[i]-sdays[0]).days for i in range(len(sdays))])+1
    for ish in ishs:
        btlist=glob.glob(os.getcwd()+'/*/'+tasks[ish]['msufile']+'_bt2_??????.nc')
        for bt in btlist:
            try:
                
                sidx=np.where(stnames==bt[-9:-3])[0][0]
                if sidx not in idx:
                    continue
            except:
                print('no RAOBCORE match for '+bt)
                continue
            with h5py.File(bt) as f:
                tidx=np.searchsorted(idays,f['datum'][0,:])
                for ih in 0,1:     
                    for ic in 0,1,2:
                        
                        tasks[ish]['msudata'][sidx,0,ih,-3+ic,tidx]=f['montemp'][ih,ic,:]
        
    print('failed,good:',jj,j)
        

def allrasotrends(path,tasks,plotproperties,intervals,days,stnames,interv,hadmed=0,hadtem=0,hadens=0,sats=0,gstations=[],daynight='',init=False):

    t=time.time()

    #dummy=append_stations(path,tasks,days,lats=[],lons=[],ps=[],stnames=[],minlen=24)
    
    ppage = page(
        layout='positional',  
        page_x_length=29.7, 
        page_y_length=21., 
        page_id_line='off',
        page_x_position=0., 
        page_y_position=0.)

    pindex=np.asarray(plotproperties['pindex'],dtype='int')
    msupindex=np.asarray(plotproperties['msupindex'],dtype='int')
#    msups=plotproperties['msups']
    msups=np.asarray([800.,550.,250.,90.])
    msunames=['TLT','TMT','TTS','TLS']
    satlist=['rss','uah','star','wegc','merra2','20CRv3']
    reanlist=['wegc','merra2','20CRv3']
    if sats==0:
        sats={}
        for s in satlist:
            sats[s]=dict()
    snlist=[]
    for k in range(len(tasks)):
        snlist.append(tasks[k]['shortname']) 

    os.chdir(path)
    startyear=tasks[0]["startdate"]//10000
    first=(intervals[0]-startyear)*12
    last=(intervals[1]-startyear+1)*12
    endyear=2022
    if first<0:
        print((os.getcwd()+': Interval ',intervals,' not covered by available data starting at ',startyear))
        return hadmed,hadtem,hadens,sats
    tolerance=intervals[2]
    reanactive=False
    for r in reanlist:
        if r in plotproperties['monthlytasks']:
            reanactive=True
            break
    slowload=len(plotproperties['plotlist'])!=1 or plotproperties['plotlist'][0]!='stations' or reanactive 


    if init:
        if 'belts' in plotproperties['plotlist'] or 'beltanomalies' in plotproperties['plotlist'] or 'hadcrut4' in snlist or 'hadcrutem4' in snlist:
            hadpath=os.path.expanduser('~/fastscratch/rise/1.0/common/')
            if 'hadcrut4' in snlist:
                hadmed,hadtem,hadens=read_hadCRUT4(hadpath+'/../common/','HadCRUT.4.6.0.0',startyear=startyear,endyear=endyear,ens=1)
            else:
                hadmed,hadtem,hadens=read_hadCRUT4(hadpath+'/../common/','HadCRUT.5.0.1.0',startyear=startyear,endyear=endyear,ens=1)
            #hadmean=np.mean(hadens,axis=0)
        else:
            hadmed=np.zeros([1,(endyear-startyear+1)*12,18,36])
            hadtem=np.zeros([1,(endyear-startyear+1)*12,18,36])
    if slowload:
        if init:
            if 'wegc' in plotproperties['monthlytasks']:
#                sats['wegc']['full'],sats['wegc']['18']=read_wegc(os.path.expandvars('$RSCRATCH/GPS/wegz/pruned/'),'5.6.2',startyear=startyear,endyear=endyear,ens=0)
                sats['wegc']['full'],sats['wegc']['18']=read_wegc10(os.path.expanduser('~/fastscratch/RODaten/RO_OPSv5.6.2_L2b_10x10'),'daily',startyear=startyear,endyear=endyear,ens=0)
                sats['wegc']['18slopes']=np.zeros([1,16,hadmed.shape[2],hadmed.shape[3]],dtype=np.float32)
            if '20CRv3' in plotproperties['monthlytasks']:
                sats['20CRv3']['full'],sats['20CRv3']['18'],sats['20CRv3']['msu']=read_20CRv3(os.path.expandvars('$RSCRATCH/20CRv3/'),'',
                                                plotproperties['ps'],startyear=startyear,endyear=endyear,ens=0)
                sats['20CRv3']['18slopes']=np.zeros([1,16,hadmed.shape[2],hadmed.shape[3]],dtype=np.float32)
            if 'merra2' in plotproperties['monthlytasks']:
                sats['merra2']['full'],sats['merra2']['18']=read_merra2(os.path.expandvars('$RSCRATCH/MERRA2/'),'',
                                                                        plotproperties['ps'],startyear=startyear,endyear=endyear,ens=0)
                sats['merra2']['18slopes']=np.zeros([1,16,hadmed.shape[2],hadmed.shape[3]],dtype=np.float32)
            if 'uah' in plotproperties['monthlytasks']:
                sats['uah']['full'],sats['uah']['18']=read_uah('/'.join(path.split('/')[:-3])+'/MSUUAHDaten/','6.0',startyear=startyear,endyear=endyear,ens=0)
                sats['uah']['18slopes']=np.zeros([1,4,hadmed.shape[2],hadmed.shape[3]],dtype=np.float32)
            if 'rss' in plotproperties['monthlytasks']:
                sats['rss']['full'],sats['rss']['18']=read_rss('/'.join(path.split('/')[:-3])+'/MSUDaten/','V4_0',startyear=startyear,endyear=endyear,ens=0)
                sats['rss']['18slopes']=np.zeros([1,4,hadmed.shape[2],hadmed.shape[3]],dtype=np.float32)
            if 'star' in plotproperties['monthlytasks']:
                #starlist=['NESDIS-STAR_TCDR_MSU-AMSUA_V04R01_TMT_S197811_E202010_C20201125.nc',
                          #'NESDIS-STAR_TCDR_MSU-AMSUA_V04R01_TUT_S198101_E202010_C20201125.nc',
                          #'NESDIS-STAR_TCDR_MSU-AMSUA_V04R01_TLS_S197811_E202010_C20201125.nc']	    
                starlist=['NESDIS-STAR_TCDR_TMT_V05R00_S197811_E202309_C20231005.nc',
                          'NESDIS-STAR_TCDR_TUT_V05R00_S198101_E202309_C20231005.nc',
                          'NESDIS-STAR_TCDR_TLS_V05R00_S197812_E202309_C20231005.nc']	    
                sats['star']['full'],sats['star']['18']=read_star('/'.join(path.split('/')[:-3])+'/MSU_STAR/',starlist,startyear=startyear,endyear=endyear,ens=0)
                sats['star']['18slopes']=np.zeros([1,4,hadmed.shape[2],hadmed.shape[3]],dtype=np.float32)

        if 'belts' in plotproperties['plotlist']:
            hadmedslopes=np.zeros([1,hadmed.shape[2],hadmed.shape[3]],dtype=np.float32)
            hadtemslopes=np.zeros([1,hadtem.shape[2],hadtem.shape[3]],dtype=np.float32)
            hadslopes=np.zeros([hadens.shape[0],hadmed.shape[2],hadmed.shape[3]],dtype=np.float32)
            first=(intervals[0]-startyear)*12
            last=(intervals[1]-startyear+1)*12
            if first<0:
                print((os.getcwd()+': Interval ',intervals,' not covered by available data starting at ',startyear))
                return hadmed,hadtem,hadens,sats
            tolerance=intervals[2]
            hadtime=np.arange(last-first,dtype=float)
            hilf=hadtime-hadtime

            tx=time.time()
            for ilat in range(hadmed.shape[2]):
                for ilon in range(hadmed.shape[3]):
        #	    print typeof(np.asarray(hadtime,dtype=float)),typeof(np.asarray(np.reshape(hadmed[0,first:last,ilat,ilon],[hadtime.shape[0]]),dtype=float))
                    hadstop=hadtime.shape[0]	    
                    if last>hadmed.shape[1]:
                        hadstop=hadtime.shape[0]-(last-hadmed.shape[1])
                        print(('HADAT ends at ',hadmed.shape[1],' !!'))
                    hadmedslopes[0,ilat,ilon]=fastlinregress(np.asarray(hadtime[0:hadstop],dtype=float),
                                                             np.asarray(np.reshape(hadmed[0,first:last,ilat,ilon],[hadstop]).flatten(),dtype=float))*120.
                    hadtemslopes[0,ilat,ilon]=fastlinregress(np.asarray(hadtime[0:hadstop],dtype=float),
                                                             np.asarray(np.reshape(hadtem[0,first:last,ilat,ilon],[hadstop]).flatten(),dtype=float))*120.
                    for ich in range(1,4):
                        for s in satlist:
                            if s in plotproperties['monthlytasks']:
                                hilf[0:hadstop]=np.reshape(sats[s]['18'][0,ich,first:last,ilat,ilon],[hadstop])
                                sats[s]['18slopes'][0,ich,ilat,ilon]=fastlinregress(hadtime[0:hadstop],hilf[0:hadstop])*120.

            print((time.time()-t))

            allhad(hadens,hadtem,hadmed,hadslopes,hadtime,first,last,hadstop)
            #for iens in range(hadens.shape[0]):
                #for ilat in range(hadmed.shape[2]):
                    #for ilon in range(hadmed.shape[3]):
                        #hadslopes[iens,ilat,ilon]=fastlinregress(hadtime[0:hadstop],hadens[iens,first:last,ilat,ilon]
                                                                    #)*120.

            print(('had:',time.time()-t))
    if tasks[-1]['shortname']=='hadat':
        istat,lats,lons,ps,stnames=loadsondes(path,tasks[:-1],days,plotproperties,stnames,init,slowload,satlist)
    else:
        istat,lats,lons,ps,stnames=loadsondes(path,tasks,days,plotproperties,stnames,init,slowload,satlist)
        
    shortnames=[t['shortname'] for t in tasks]
    if 'suny' in shortnames or 'sunyhom' in shortnames:
        fill_suny(tasks,shortnames,stnames,lats,lons,ps,plotproperties['pindex'])
    t=time.time()
    if isinstance(stnames,list):
        stnames=np.asarray(stnames)
    #istnames=stnames.astype(np.int)
    tmshape=tasks[0]["data"].shape

    slopes=[]
    sli=-1
    for k in range(len(tasks)):
        if tasks[k]['shortname'] in list(sats.keys()):
            if tasks[k]['shortname'] in reanlist:
                tasks[k]['data']=satst(lats,lons,sats[tasks[k]['shortname']]['full'],tmshape[4])
                try:
                    tasks[k]['msudata']=satst(lats,lons,sats[tasks[k]['shortname']]['msu'],tmshape[4])
                except:
                    pass
                    
            else:
                tasks[k]['msudata']=satst(lats,lons,sats[tasks[k]['shortname']]['full'],tmshape[4])

        for inter in range(intervals.size//3):
            sli+=1
            slopes.append({ "shortname":tasks[k]["shortname"],                            "ens":tasks[k]["ens"],
                            "interval":intervals[:],
                            "data":np.empty([tmshape[0],len(tasks[k]["ens"]),tmshape[2]+1,
                                                pindex.shape[0]+msupindex.shape[0]],np.float32),
                            }
                          )
            slopes[sli]["data"].fill(np.nan)


    ni=36
    nj=18
    glons=5.+10.*np.arange(ni)
    glats=-85.+10.*np.arange(nj)
    belts=np.asarray([[0,18],[11,18],[0,8],[7,11],[9,18],[0,9]])
    beltnames=["Globe",'NHEx','SHEx','Tropics','NH','SH']

    gstatindex=np.zeros([nj,ni,100],dtype=np.int32) # 20 is max stations per 10x10 degree box
    find_gstatindex(glons,glats,lons,lats,gstatindex)

    currentdatabg=np.zeros([tmshape[0],tmshape[2],pindex.shape[0],tmshape[4]],np.float32)
    currentdatatm=np.zeros([tmshape[0],tmshape[2],pindex.shape[0],tmshape[4]],np.float32)
    currentdata  =np.zeros([tmshape[0],tmshape[2],pindex.shape[0],tmshape[4]],np.float32)
    #currentdatauah=np.zeros([tmshape[0],tmshape[2],msupindex.shape[0],tmshape[4]],np.float32)
    #currentdatarss=np.zeros([tmshape[0],tmshape[2],msupindex.shape[0],tmshape[4]],np.float32)
    currentdatamsu=np.zeros([tmshape[0],tmshape[2],msupindex.shape[0],tmshape[4]],np.float32)
    jcurrentdata=np.zeros((1,1,1,1))

    keytmcorr=-1
    keytm=-1
    shortnames=[]
    for key in range(len(tasks)):
        if 'stations' in plotproperties['plotlist'] or 'cost' in plotproperties['plotlist']:
            t=time.time()
            dists=np.zeros((tmshape[0]+1)*tmshape[0]//2,np.float32)
            x,y,z=tdist(dists,lats,lons,1)
            print((dists[:5]))
            print(('tdist',time.time()-t))

        shortnames.append(tasks[key]['shortname'])
        if tasks[key]['shortname'] in satlist and tasks[key]['shortname']  not in reanlist:
            tasks[key]['cost']=np.empty((tasks[key]['msudata'].shape[1],3,msupindex.shape[0]))
        else:
            tasks[key]['cost']=np.empty((tasks[key]['data'].shape[1],3,pindex.shape[0]+msupindex.shape[0]))
        tasks[key]['cost'].fill(np.nan)
        if tasks[key]["shortname"]=="tmcorr":
            keytmcorr=key
        if tasks[key]["shortname"]=="tm":
            keytm=key        
        #if tasks[key]["shortname"]=="rss":
            #keyrss=key        
            #expandandadd(tasks[keyrss]["msudata"][:,:,:,:currentdata.shape[3]],currentdatamsu,tasks[keytm]["index"],msupindex,0,currentdatamsu,0.0)
        #if tasks[key]["shortname"]=="uah":
            #keyuah=key        
            #expandandadd(tasks[keyuah]["msudata"][:,:,:,:currentdata.shape[3]],currentdatamsu,tasks[keytm]["index"],msupindex,0,currentdatamsu,0.0)
        #if tasks[key]["shortname"]=="star":
            #keystar=key        
            #expandandadd(tasks[keystar]["msudata"][:,:,:,:currentdata.shape[3]],currentdatamsu,tasks[keytm]["index"],msupindex,0,currentdatamsu,0.0)

    andeplist,andepstdlist,andeprmslist=make_andeplist()

    sli=-1
    nodata=[]
    cini=False
    print(('prep:', time.time()-t))
    for key in range(len(tasks)):
        print((tasks[key]["shortname"]))
        if tasks[key]["shortname"] in [ "rio", "rit","riocorr", "ritcorr","ce20c_andep"]:
            enslist=tasks[key]["ens"]
        else:
            enslist=[0]

        nodata.append(False)    
        for iens in enslist: #range(len(tasks[key]["ens"])):

            # replace RAOBCORE adjustments with unadjusted series
            sat=False
            t=time.time()
            dkey='data'
            ti=pindex[-1]
            if tasks[key]['shortname'] in satlist and tasks[key]['shortname'] not in reanlist:
                dkey='msudata'
                ti=-1
            if nodatafound(tasks[key][dkey],testindex=ti) and snlist[key]!='hadat':
                nodata[-1]=True
                break
            if tasks[key]["shortname"] == "tm":
                expandandadd(tasks[key]["data"],tasks[keytmcorr]["data"].reshape(tmshape[0],tmshape[2],tmshape[3],tmshape[4]),
                             tasks[key]["index"],pindex,iens-iens,currentdatatm,+1.0)
#                print currentdata[1171,0,0,:]-currentdatatm[1171,0,0,:]
                currentdata[:]=currentdatatm[:]
            elif tasks[key]["shortname"] == "milancorr":
#                expandandadd(tasks[key]["data"],tasks[keytm]["data"].reshape(tmshape[0],tmshape[2],tmshape[3],tmshape[4]),
#                             tasks[key]["index"],pindex,0,currentdata,+1.0)
#                print currentdata[1171,0,0,:]-currentdatatm[1171,0,0,:]
                ms=tasks[key]["data"].shape
                currentdata[:]=currentdatatm-(tasks[key]["data"].reshape((ms[0],ms[2],ms[3],ms[4])))[:,:,pindex,:]
            elif tasks[key]["shortname"] in ["erai_fggpsdep","erai_fggpswetdep"]: 
                picopy(currentdata,currentdataeibg,-tasks[key]["data"],iens,pindex)
                currentdata[:,:,:,:102*12]=np.nan
            elif tasks[key]["shortname"] == "eijra_fgdep": 
                picopy(currentdata,currentdatatm,tasks[key]["data"],iens,pindex)
                currentdataeibg=currentdata.copy()
            elif tasks[key]["shortname"] in andeplist:
#                currentdata[:]=currentdatatm+tasks[key]["data"].reshape(tmshape[0],tmshape[2],ps.shape[0],tmshape[4])[:,:,pindex,:]
                picopy(currentdata,currentdatatm,tasks[key]["data"],iens,pindex)
            elif tasks[key]["shortname"] in ["ce20c_andep"]:
#                currentdata[:]=currentdatatm+tasks[key]["data"][:,iens,:,:,:].reshape(tmshape[0],tmshape[2],ps.shape[0],tmshape[4])[:,:,pindex,:]
                picopy(currentdata,currentdatatm,tasks[key]["data"],iens,pindex)
            elif ('rio' in tasks[key]["shortname"] or 'rit' in tasks[key]["shortname"]) and 'corr' not in tasks[key]["shortname"] :  
                expandandadd(tasks[key]["data"],currentdatatm,tasks[key]["index"],pindex,iens,currentdata,-1.0)
            elif 'riocorr' in tasks[key]["shortname"] or 'ritcorr' in tasks[key]["shortname"]:  
                expandandadd(tasks[key]["data"],currentdatatm-currentdatatm,tasks[key]["index"],pindex,iens,currentdata,-1.0)

            elif 'era5v' in tasks[key]["shortname"] or 'eraibc' in tasks[key]["shortname"]:  
#		expandandadd(tasks[key]["data"],currentdatabg,tasks[key]["index"],pindex,iens,currentdata,-1.0)
                picopy(currentdata,currentdatatm,tasks[key]["data"],iens,pindex)

            elif tasks[key]["shortname"] == "rcorr":
                zero=currentdatatm-currentdatatm
                expandandadd(tasks[key]["data"],zero,tasks[key]["index"],pindex,iens,currentdata,-1.0)
            elif tasks[key]["shortname"] == "bg": 
                expandandadd(tasks[key]["data"],currentdata,tasks[keytm]["index"],pindex,iens-iens,currentdatabg,0.0)
                currentdata[:]=currentdatabg.copy()
            elif tasks[key]["shortname"] in ['suny','sunyhom']: 
                picopy(currentdata,currentdatatm-currentdatatm,tasks[key]["data"],iens,pindex)
            elif tasks[key]["shortname"] == "bgdiff":  
                currentdata[:]=currentdatatm-currentdatabg
            elif tasks[key]["shortname"] in reanlist:  
                picopy(currentdata,currentdatatm-currentdatatm,tasks[key]["data"],iens,pindex)
            elif tasks[key]["shortname"] in satlist:  
                #currentdata=np.asarray([0.])
                sat=True
            else:
                if snlist[key]!='hadat':
                    expandandadd(tasks[key]["data"],currentdata,tasks[keytm]["index"],pindex,iens-iens,currentdata,0.0)


            print(('expand1',iens,time.time()-t))
            if slowload and snlist[key]!='hadat':
                if tasks[key]["msudata"].shape[0]!=tasks[keytm]['msudata'].shape[0]:
                    tasks[key]["msudata"]=np.empty(tasks[keytm]['msudata'].shape,dtype=np.float32)
                    tasks[key]["msudata"].fill(np.nan)
                else:
                    if not cini:
                        cmask=np.isnan(tasks[keytm]['msudata'])
                        cini=True
#		    tasks[key]["msudata"][cmask]=np.nan
                ss=tasks[key]["msudata"].shape
                sm=currentdatamsu.shape
                if ss[0]==sm[0] and ss[4]==sm[3]:
                    currentdatamsu=tasks[key]["msudata"][:,iens,:,:,:]
                    if tasks[key]['shortname']=='tm':
                        tmindex=np.where(np.isnan(currentdatamsu))
                    else:
                        if tasks[key]['type']=='SAT':
                            currentdatamsu[tmindex]=np.nan
                        
                else:
                    expandandadd(tasks[key]["msudata"],currentdatamsu,tasks[keytm]["index"],msupindex,iens,currentdatamsu,0.0)
#		currentdatamsu[:,:,:,:(1979-tasks[key]['startdate']/10000)*12]=np.nan
#		expandandadd(tasks[key]["msudata"][:,:,:,:,:currentdatamsu.shape[3]],currentdatamsu,tasks[keytm]["index"],msupindex,iens,currentdatamsu,0.0)

            if len(gstations)>0:
                try:
                    for stn in range(stnames.shape[0]):
                        if stnames[stn] not in gstations:
                            currentdata[stn,:,:,:]=np.nan
                except:
                    pass

            print(('expand2',iens,time.time()-t))
            try:
                if daynight=='Day':
                    mask=np.logical_and(lons>-90.,lons<90.)
                    currentdata[mask,0,:,:]=np.nan
                    currentdata[np.logical_not(mask),1,:,:]=np.nan		    
                elif daynight=='Night':
                    mask=np.logical_and(lons>-90.,lons<90.)
                    currentdata[mask,1,:,:]=np.nan
                    currentdata[np.logical_not(mask),0,:,:]=np.nan
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
                interval=np.asarray(intervals[:],dtype=np.int32)
                stop=interval[1]
                start=interval[0]
                stime=(start-startyear)*12
                itime=stime+np.arange((stop-start+1)*12)
                if iens==0:
                    sli+=1

                t=time.time()

                s=slopes[sli]["data"]

                if not sat:
#		    jcurrentdata=np.concatenate((currentdata,currentdatamsu),axis=2)
#		    print 'if',time.time()-t
                    cs=currentdata.shape
                    if slowload:
                        if jcurrentdata.shape[2]!=cs[2]+currentdatamsu.shape[2]:
                            jcurrentdata=np.zeros((cs[0],cs[1],cs[2]+currentdatamsu.shape[2],cs[3]),dtype=np.float32)
                        print(('slope',time.time()-t))				   
                        ncat(currentdata,currentdatamsu,jcurrentdata)
                        print(('nslope',time.time()-t))				   
                    else:
                        jcurrentdata=currentdata
#		    print 'if',time.time()-t
                    anomalies=jcurrentdata #np.copy(jcurrentdata)

                    jpindex=np.concatenate((pindex,msupindex))
#		    print 'if',time.time()-t
                else:
                    print("should not be there")
#                    exit()
                    jcurrentdata=currentdatamsu[:]
#                    anomalies=jcurrentdatamsu #np.copy(currentdatamsu[:])
                    anomalies=np.copy(currentdatamsu[:])
                    jpindex=msupindex
#		    print 'else',time.time()-t
                    #anomalies_and_slopes(currentdatamsu,startyear,interval,tolerance,iens,itime,orig,anomaly,anomalies,climatology,
                                        #good,s)
                if 'good' not in locals():
                    good=np.zeros([tmshape[0],tmshape[2],jpindex.shape[0]],np.int32)
                    climatologies=np.zeros([tmshape[0],tmshape[2],jpindex.shape[0],12],np.float32)
                else:
                    if jpindex.shape[0]==good.shape[2]:
                        good[:]=np.zeros([tmshape[0],tmshape[2],jpindex.shape[0]],np.int32)
                        climatologies[:]=np.zeros([tmshape[0],tmshape[2],jpindex.shape[0],12],np.float32)
                    else:
                        good=np.zeros([tmshape[0],tmshape[2],jpindex.shape[0]],np.int32)
                        climatologies=np.zeros([tmshape[0],tmshape[2],jpindex.shape[0],12],np.float32)
                if slowload or tasks[key]['shortname']==plotproperties['monthlytasks'][0] or tasks[key]['shortname'] in ['tm','tmcorr']:

                    print(('slope',time.time()-t))				   
                    #good=np.zeros([tmshape[0],tmshape[2],jpindex.shape[0]],np.int)
                    #climatologies=np.zeros([tmshape[0],tmshape[2],jpindex.shape[0],12],np.float32)
#		    for x in jcurrentdata,startyear,interval,int(tolerance),int(iens),itime,anomalies,climatologies,good,s:               
#			print typeof(x)

                    if tasks[key]['data'].size>1:
                        if jcurrentdata.shape[2]>2:
                            pidx=np.array(list(np.arange(pindex.shape[0]))+list(msupindex+pindex.shape[0]))
                        else:
                            pidx=np.array([0])
                    else:
                        pidx=msupindex
                    if slowload:
                        #if 'anomalies' not in locals():
                            #anomalies=jcurrentdata[:,:,pidx,:]
                        #else:
                            #anomalies[:]=jcurrentdata[:,:,pidx,:]

                        anomalies_and_slopes(jcurrentdata[:,:,pidx,:],startyear,interval,int(tolerance),int(iens),itime,anomalies,
                                             climatologies,good,s)
                    else:
                        #if 'anomalies' not in locals():
                            #anomalies=np.copy(jcurrentdata)	
                        #else:
                            #anomalies[:]=jcurrentdata[:]
                        anomalies_and_slopes(jcurrentdata,startyear,interval,int(tolerance),int(iens),itime,anomalies,
                                             climatologies,good,s)

    # remove India             		    
#                    mnames=read_mesural(os.path.expanduser('~/tables/MESURAL'))   
                    mnames=read_viz(os.path.expanduser('~/tables/cards_meta.prn'))   

#		    mnames=['']
                    #for m in mnames:
                        #idx=np.where(m==stnames)[0]
                        #if len(idx)>0:
                            #anomalies[idx[0],:,:,:]=np.nan
                            #s[idx[0],:,:,:]=np.nan
                    idx=0
                    #for m in stnames:
                        #if m not in mnames:
                            #anomalies[idx,:,:,:]=np.nan
                            #s[idx,:,:,:]=np.nan
                        #idx+=1

                    #mask=(istnames>68000)|(istnames<60000)
                    #s[mask,:,:,:]=np.nan
                    #anomalies[mask,:,:,:]=np.nan
                    #mask=(istnames>20000)&(istnames<40000)
                    #mask=(istnames<50000)|(istnames>60000)
                    #s[mask,:,:,:]=np.nan
                    #anomalies[mask,:,:,:]=np.nan
                    #mask=(istnames<70000)|(istnames>75000)
                    #s[mask,:,:,:]=np.nan
                    #anomalies[mask,:,:,:]=np.nan
                    mask=(stnames>'061900')&(stnames<'062000')
                    s[mask,:,:,:]=np.nan
                    anomalies[mask,:,:,:]=np.nan
                    #mask=(istnames>84000)&(istnames<85000)
                    #s[mask,:,:,:]=np.nan
                    #anomalies[mask,:,:,:]=np.nan
                    mask=(stnames>'042000')&(stnames<'044000')
                    s[mask,:,:,:]=np.nan
                    anomalies[mask,:,:,:]=np.nan
                    #mask=(istnames>48600)&(istnames<48900)
                    #s[mask,:,:,:]=np.nan
                    #anomalies[mask,:,:,:]=np.nan

                    s[:,:,2,:]=np.nanmean(s[:,:,:2,:],axis=2)
                print(('slope',time.time()-t))
                if slowload or 'hadcrut4' in snlist or 'hadcrutem4' in snlist:
                    t=time.time()

                    itime=np.asarray(itime,dtype=int)
                    #try:
                        #if ganomalies.shape[3]!=jpindex.shape[0]:
                            #ganomalies=np.zeros([nj,ni,tmshape[2],jpindex.shape[0],tmshape[4]],np.float32)
                            #gslopes=np.zeros([nj,ni,tmshape[2],jpindex.shape[0]],np.float32)
                    #except:
                    ganomalies=np.zeros([nj,ni,2,jpindex.shape[0],tmshape[4]],np.float32)
                    gslopes=np.zeros([nj,ni,2,jpindex.shape[0]],np.float32)

                    sy=tasks[0]['startdate']//10000
                    if 'had' in snlist[key]:
                        if snlist[key]=='hadcrut4':
                            h=hadmed
                        else:
                            h=hadmed
                        for j in range(h.shape[2]):
                            for i in range(h.shape[3]):
                                for k in range(12):
                                    cg=h[0,(start-sy)*12+k:(stop-sy)*12+k:12,j,i]
                                    if np.sum(np.isnan(cg))<2:
                                        c=np.nanmean(cg)
                                        if ganomalies.shape[4]<=h.shape[1]:
                                            
                                            ganomalies[j,i,0,0,k::12]=h[0,k:ganomalies.shape[4]:12,j,i]-c
                                            
                                        else:
                                            ganomalies[j,i,0,0,k:h.shape[1]:12]=h[0,k::12,j,i]-c
                                            ganomalies[j,i,0,0,h.shape[1]+k::12]=np.nan
                                    else:
                                        ganomalies[j,i,0,0,k::12]=np.nan

                        ganomalies[:,:,1,0,:]=np.nan
                        for j in range(1,ganomalies.shape[3]):
                            ganomalies[:,:,0,j,:]=ganomalies[:,:,0,0,:]
                            ganomalies[:,:,1,j,:]=np.nan
                        ganomalies[gamask]=np.nan
                        tim=np.arange((stop-start+1)*12)/12.
                        istart=(start-1900)*12
                        istop=(stop+1-1900)*12
                        for l in range(ganomalies.shape[0]):
                            for k in range(ganomalies.shape[1]):
                                #for j in range(ganomalies.shape[3]):
                                if sum(np.isnan(ganomalies[l,k,0,0,istart:istop]))<tolerance*12:
                                    gslopes[l,k,0,:]=fastlinregress(tim,ganomalies[l,k,0,0,istart:istop])*10
                                    gslopes[l,k,0,:]=np.nan
                                else:
                                    gslopes[l,k,0,:]=np.nan
                        gslopes[:,:,1,:]=gslopes[:,:,0,:]
#			gslopes[:,:,2,:]=gslopes[:,:,0,:]

                    else:
                        grid_anomalies(anomalies,good,int(tolerance),gstatindex,ganomalies,gslopes,start,stop,itime)
                        if snlist[key]=='tmcorr':
#			    ganomalies[:,:,2,:,:]=np.nanmean(ganomalies,axis=2)
                            gamask=np.isnan(ganomalies)
                    if snlist[key] in satlist and snlist[key] not in reanlist:
                        pass
                        #ganomalies[gamask[:,:,:,-1:,:]]=np.nan
                    else:
                        ganomalies[gamask]=np.nan


                    print(('grid',time.time()-t))
                    g=ganomalies.shape
    #		zanomalies=np.zeros([g[0],g[3],g[4]])
    #		ndnanmean(ganomalies,zanomalies)
                    
                    if False:
                        tasks[key]['ganomalies'] = ganomalies.copy()
                        tasks[key]['anomalies'] = anomalies.copy()
                    t=time.time()
                    
                    if False and key == 2:
                        sdiff = tasks[1]['anomalies'][:, :, 4, :] - tasks[0]['anomalies'][:, :, 4, :]
                        idx = np.unravel_index(np.nanargmax(sdiff), sdiff.shape)
                        print(idx, sdiff[idx])
                        
                        for iii in range(7, 11):
                            
                            sdiff = tasks[1]['anomalies'] - tasks[0]['anomalies']
                            gdiff = tasks[1]['ganomalies'][iii:iii + 1] - tasks[0]['ganomalies'][iii:iii + 1]
                            idx = np.unravel_index(np.nanargmax(gdiff), gdiff.shape)
                            print(iii, idx, gdiff[idx])
                            gdm = np.nanmean(gdiff, axis=(0, 1, 2, 3))
                            plt.plot(gdm, label=str(iii))
                        plt.legend()
                        plt.show()

                    zanomalies=ndnanmean(ganomalies,2)
                    weights=np.empty(zanomalies.shape)
                    zslopes=np.zeros([nj,ni,1])
                    itime=np.arange((stop-start+1)*12)
                    istart=(start-tasks[key]['startdate']//10000)*12
                    istop=(stop-tasks[key]['startdate']//10000+1)*12
                    for k in range(zanomalies.shape[0]):
                        weights[k,:,:]=np.cos((-85.+k*10)*np.pi/180.)
                    if 'had' in tasks[key]['shortname']:
                        for k in range(zanomalies.shape[0]):
                            for l in range(ganomalies.shape[1]):
                                gslopes[k,l,:]=fastlinregress(itime,ganomalies[k,l,0,0,istart:istop])*120.
#		    tasks[key]['zslopes'][:]=np.nan
#		    tasks[key]['zslopes'][0,:,0]=zslopes[:,0]
                    #gdanomalies=np.nanmean(ganomalies,axis=2)
                    #zanomalies=np.nanmean(gdanomalies,axis=1)
                    zs=zanomalies.shape
                    beltanomalies=np.zeros([len(beltnames),zs[1],zs[2]])
                    hadmedanomalies=np.zeros([len(beltnames),zs[1]])
                    for b in range(len(beltnames)):
                        #if snlist[key]!='hadat':
                        beltanomalies[b,:,:]=np.nanmean(zanomalies[belts[b,0]:belts[b,1],:,:]*weights[belts[b,0]:belts[b,1],:,:],axis=0)/ \
                            np.nanmean(zanomalies[belts[b,0]:belts[b,1],:,:]/zanomalies[belts[b,0]:belts[b,1],:,:]*weights[belts[b,0]:belts[b,1],:,:],axis=0)
                        #if tasks[key]['shortname'] in satlist and tasks[key]['shortname'] not in ['wegc','merra2']:
                            #idx=np.where(np.sum(~np.isnan(zanomalies[belts[b,0]:belts[b,1],len(pindex)+msupindex[0],:]),axis=0)<=1)
                        #else:
                            #idx=np.where(np.sum(~np.isnan(zanomalies[belts[b,0]:belts[b,1],0,:]),axis=0)<=1)
                        #beltanomalies[b,:,idx[0]]=np.nan

                        #else:
                        #for z in range(zs[1]):
                            #beltanomalies[b,z,:]=np.nanmean(hadmed[0,:,belts[b,0]:belts[b,1],:],axis=(1,2))

                    if tasks[key]['shortname'] in satlist and tasks[key]['shortname'] not in reanlist:
                        #tasks[key]['beltanomalies'][iens,:,:,:]=np.empty((1,beltanomalies.shape[0],beltanomalies.shape[1]+pindex.shape[0],beltanomalies.shape[2]))
                        tasks[key]['beltanomalies'][iens,:,:,:]=np.empty((1,beltanomalies.shape[0],beltanomalies.shape[1],beltanomalies.shape[2]))
                        tasks[key]['beltanomalies'][iens,:,:,:]=np.nan
                        tasks[key]['beltanomalies'][iens,:,-4:,:]=beltanomalies[:]
                    else:
                        tasks[key]['beltanomalies'][iens,:,:,:]=np.copy(beltanomalies[:])
                    #plt.plot(1900.+itime/12.,beltanomalies[0,15,itime])
                    print((time.time()-t))

                    if ganomalies.shape[4]>10000 and jpindex.shape[0]>13 and 4 in jpindex and 0 in jpindex:
                        gclimatologies=np.zeros([nj,ni,tmshape[2],jpindex.shape[0],12],np.float32)
                        gcslopes=np.zeros([nj,ni,tmshape[2],jpindex.shape[0]],np.float32)
                        grid_anomalies(climatologies,good,0,gstatindex,gclimatologies,gcslopes,0,0,itime)
                        save_gridded(ganomalies,gclimatologies,tasks,key,days,ps,pindex,iens=iens,
                                     start=1991,stop=2020,version='1.9',append=False)

                estr=''	    
                clev=np.linspace(-1.2,1.2,25)
                clev=np.append(clev,10.)
                clev=np.append(-10.,clev)
                try:
                    if plotproperties['double']:
                        clev*=2
                except:
                    pass
                for ipar in range(3):
                    parstr="{0:0>2}".format(ipar*12)
                    if 'satseries' in plotproperties['plotlist']:
                        ipmax=jpindex.shape[0]
                    else:
                        ipmax=pindex.shape[0]

                    for ip in range(ipmax):

                        if sum(s[:,0,ipar,ip]==s[:,0,ipar,ip])==0:
                            continue

                        if ip<pindex.shape[0] and jpindex.shape[0]!=msupindex.shape[0]:
                            pstr="{0:0>3}".format(ps[jpindex[ip]].astype(np.int32))
                            psuff=' hPa'
                        else:
                            pstr=msunames[jpindex[ip]]
                            psuff=''

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
                        if "stations" in plotproperties["plotlist"]:
                            nprefix+='stations_'

                        oname=nprefix+'trends_'+tasks[key]["shortname"]+estr+'_'+\
                            "{:4}".format(interval[0])+"-{:4}".format(interval[1])+'_'+pstr+'_'+parstr    
                        out = output({"output_name_first_page_number":"off",
                                      "output_formats":plotproperties["outputformat"], 
                                      'output_name':plotproperties['tmppath']+'/'+oname})

                        if "stations" in plotproperties["plotlist"] or "cost" in plotproperties["plotlist"]:

                            #print 'stations output file: ',oname
                            #if plotproperties['exps'][0]=='exp07':
                                #xg=np.load('/vgc/srvx7/leo/scratch/rise/1.0/exp08/stnames.npy')
                            #elif plotproperties['exps'][0]=='exp08':
                                #xg=np.load('/vgc/srvx7/leo/scratch/rise/1.0/exp07/stnames.npy')
                            #else:
                                #xg=[]

                            #ii=0
                            #for st in stnames:
                                #if st not in xg:
                                    #s[ii,iens,ipar,ip]=np.nan
                                #ii+=1

                            scosts=np.zeros(s.shape[0])
                            cost=tcost(dists,s[:,iens,ipar,ip],scosts)

                            if ip==jpindex.shape[0]-1 and ipar==1:
                                print(('cost',time.time()-t))
                            cstr="{0:8.2f}".format(cost)
                            if ip<pindex.shape[0]:
                                tasks[key]['cost'][iens,ipar,ip]=cost
                        if "stations" in plotproperties["plotlist"]:
                            statstr="{:4}".format(sum(~np.isnan(s[:,iens,ipar,ip])))
                            lines =["Temperature Trends [K/10a], "+tasks[key]["shortname"]+estr+', '+plotproperties['exps'][interv]+', '+
                                    "{:4}".format(interval[0])+"-{:4}".format(interval[1])+', '+parstr+"h, "+pstr+psuff,
                                    str(statstr+' Stations, Cost: '+cstr+', '+plotproperties["version"]+', '+plotproperties["fgdepname"])]


                            vec=s[0:istat,iens,ipar,ip].flatten()
                            mask=~np.isnan(vec)
                            inp=minput({'Input_values':vec[mask].astype(np.float32),  #.ravel?
                                        'Input_latitude_values':lats[mask].astype(np.float32),
                                        'Input_longitude_values':lons[mask].astype(np.float32),
                                        'Input_Type':'geographical'})

                        else:
                            lines =["Temperature Trends [K/10a], "+tasks[key]["shortname"]+estr+', '+
                                    "{:4}".format(interval[0])+"-{:4}".format(interval[1])+', '+parstr+"h, "+pstr+psuff,
                                    str(plotproperties["version"])]

                        if "stations" in plotproperties["plotlist"]:
                            cind=np.argsort(-scosts[mask])
                            maxs=cind.shape[0]
                            if str(plotproperties['plotstatids'])=='False':
                                maxs=3
                                ilav=np.concatenate([[],np.asarray([-80.,-80.,-80.])+1.0])
                                ilov=np.concatenate([[],np.asarray([-20.,0.,20.])+1.0])
                            else:
                                ilav=np.concatenate([lats[mask][cind[0:maxs]].astype(np.float32),
                                                        np.asarray([-80.,-80.,-80.])+1.0])
                                ilov=np.concatenate([lons[mask][cind[0:maxs]].astype(np.float32),
                                                        np.asarray([-20.,0.,20.])+1.0])
                            
                            with open(plotproperties['tmppath']+'/'+'stnames_'+tasks[key]["shortname"]+estr+'_'+
                                  "{:4}".format(interval[0])+"-{:4}".format(interval[1])+'_'+pstr+'_'+parstr ,'w') as f:
                                f.write('statid,trend,cost\n')
                                for a,b,c in zip(stnames[mask][cind[:maxs]],
                                                 vec[mask][cind[:maxs]].astype(np.float32),
                                                 scosts[mask][cind[:maxs]]):
                                    f.write('{},{:5.2f},{:5.4f}\n'.format(a,b,c))
                                
                            projection,coastlines,title,legend,symb,symb2,symb3,symb4,cont=\
                                set_trendmaps(lines,clev,plotproperties,
                                              stnames=np.concatenate([stnames[mask][cind[:maxs]],stnames[mask][cind[:maxs]]]),
                                              slopes=np.concatenate([vec[mask][cind[:maxs]].astype(np.float32),vec[mask][cind[:maxs]]]),
                                              costs=np.concatenate([scosts[mask][cind[:maxs]],scosts[mask][cind[:maxs]]]),
                                              map=plotproperties['map'])


                            inp3=minput({'Input_latitude_values':ilav-0.5,
                                         'Input_longitude_values':ilov-1.3,
                                         'Input_Type':'geographical'})
                            inp4=minput({'Input_latitude_values':ilav-2.0,
                                         'Input_longitude_values':ilov+2.0,
                                         'Input_Type':'geographical'})
                            inp5=minput({'Input_latitude_values':ilav-4.0,
                                         'Input_longitude_values':ilov+2.0,
                                         'Input_Type':'geographical'})
                        else:
                            projection,coastlines,title,legend,symb,cont=set_trendmaps(lines,clev,plotproperties,map=plotproperties['map'])

                        if "gridded" in plotproperties["plotlist"]:
                            glats=np.linspace(-94.999,94.999,nj+2)
                            glons=np.linspace(5,375,ni+2)

                            hilf=np.empty([nj+2,ni+2])
                            hilf[:]=np.nan

                            if ipar==0 and ip==9:
                                print((out.args['output_name']))
                            if ipar==2:
                                hilf[1:nj+1,1:ni+1]=np.nanmean(gslopes[:,:,0:2,ip],axis=2).reshape([nj,ni]).copy()
                            else:
                                hilf[1:nj+1,1:ni+1]=gslopes[:,:,ipar,ip].reshape([nj,ni]).copy()
    #                        hilf[1:nj+1,1:ni+1]=hadmedslopes[0,:,:].reshape([nj,ni]).copy()
                            mask2=np.isnan(hilf)
                            hilf[mask2]=-1.e21

                            inp2=minput({'Input_field':hilf[:,1:ni+1],  #.ravel?
                                         'Input_field_initial_latitude':glats[0]-5,
                                         'Input_field_initial_longitude':glons[0]-5,
                                         'Input_field_latitude_step':(glats[1]-glats[0]),
                                         'Input_field_longitude_step':glons[1]-glons[0],
                                         'Input_Type':'geographical'})

                        #plot(out,projection,coastlines,inp,symb,title,legend,inp2,cont)
                            out.args['output_title']=os.getcwd()+'/'+out.args['output_name']
                            if len(vec[mask])>0:
                                plot(out,projection,inp2,cont,inp,symb,title,legend,coastlines)
                                print((out.args['output_name']+'.'+plotproperties["outputformat"][0]))

                        else:
                            if "stations" in plotproperties["plotlist"]:
                                if np.sum(mask)!=0:
                                    out.args['output_title']=os.getcwd()+'/'+out.args['output_name']
#				    if plotproperties['plotstatids']=='True':
                                    for ix in range(1):
                                        ptime=time.time()
                                        plot(out,projection,inp,symb,inp3,symb3,inp4,symb4,inp5,symb2,title,legend,coastlines)
                                        print((out.args['output_name']+'.'+plotproperties["outputformat"][0]))
                                    #del out,projection,inp,symb,inp3,symb3,inp4,symb4,inp5,symb2,title,legend,coastlines
                                        print(('ptime:',tasks[key]["shortname"],iens,ip,ipar,time.time()-ptime))    
                                    print('xtime')
                                #else:
                                        #print 'plot inp:',inp
#					plot(out,projection,inp,symb,title,legend,coastlines)
                                else:
                                    print((out.args['output_name']+': no valid data'))


                if not slowload:
                    if tasks[key]['shortname']==plotproperties['monthlytasks'][0]:
                        return hadmed,hadtem,hadens,sats
                    continue

                t=time.time()
                s=gslopes.shape
                zslopes=np.zeros([s[0],s[3]])
                zslopes=zonaltrends(gslopes,zslopes)
                if 'belts' in plotproperties['plotlist']:
                    hadmedzslopes=np.zeros([nj,1])
                    hadmedzslopes=zonaltrends(np.reshape(hadmedslopes[0,:,:],[nj,ni,1,1]),hadmedzslopes)
                    hadzslopes=np.zeros([hadens.shape[0],nj,1])
                    hilf=np.zeros([nj,1])
                    for jens in range(hadens.shape[0]):
                        hadzslopes[jens,:,:]=zonaltrends(np.reshape(hadslopes[jens,:,:],[nj,ni,1,1]),hilf)

                tasks[key]["zslopes"][iens,:,:]=zslopes.copy()
                print(('zslopes',time.time()-t))
                if 'zonal' in plotproperties['plotlist'] and zslopes.shape[1]>=pindex.shape[0] and pindex.shape[0]>1:
                    out = output({"output_name_first_page_number":"off",
                                  "output_formats":plotproperties["outputformat"], 
                                  'output_name': plotproperties['tmppath']+'/'+'trendszonal_'+tasks[key]["shortname"]+estr+'_'+
                                  "{:4}".format(interval[0])+"-{:4}".format(interval[1])})


                    hilf=zslopes.T#[pindex,:]
                    hilf=hilf[:len(pindex),:]

                    lines =["Temperature Trends [K/10a], "+tasks[key]["shortname"]+estr+', '+plotproperties['exps'][interv]+', '+
                            "{:4}".format(interval[0])+"-{:4}".format(interval[1]),
                            stats(hilf[hilf==hilf],mima=1,short=1)]

                    hilf[np.isnan(hilf)]=-1.e21

                    projection,horizontal,vertical,title,legend,cont=set_zonalcross(lines,clev,plotproperties)

                    ppage = page(
                        layout='positional',  
                        page_x_length=29.7, 
                        page_y_length=21., 
                        page_id_line='off',
                        page_x_position=0., 
                        page_y_position=0.)
                    out.args['output_title']=os.getcwd()+'/'+out.args['output_name']

                    xa=glats#[1:19]
                    ya=np.asarray(ps[np.asarray(pindex,dtype=int)],np.float64)

                    inputs=minput({'Input_field':hilf,  #.ravel?
                                   'Input_field_x_list':xa,
                                   'Input_field_y_list':ya,
                                   'Input_type':'cartesian'
                                   })

#		    if 'belt' in plotproperties['plotlist']:
                    plot(out,ppage,projection,horizontal,vertical,inputs,cont,title,legend)
                    print((out.args['output_name']+'.'+plotproperties["outputformat"][0]))


        if nodata[-1]:
            print((tasks[key]["shortname"]+': does not contain valid data, will be removed'))
            continue
    #for n in range(len(nodata)-1,0,-1):
        #del tasks[n]

    # Zonal SAT profiles

    if 'satseries' in plotproperties['plotlist'] and 'zonal' in plotproperties['plotlist']:

        for ip in range(jpindex.shape[0]):

            profsat='sat'
            if ip<pindex.shape[0]:
                suff='{:0>3}'.format(int(ps[pindex[ip]]))	    
            else:
                suff=msunames[msupindex[0]]
            out = output({"output_name_first_page_number":"off",
                          "output_formats":plotproperties["outputformat"], 
                          'output_name': plotproperties['tmppath']+'/'+profsat+'trendszonal_'+suff+
                          "_{:4}".format(interval[0])+"-{:4}".format(interval[1])
                          })

            ti=prepare_ti(tasks)

            legend_user_lines=[]
            linesplustext=[]
            for key in ti:
                if nodata[key]:
                    continue

                legname=copy.copy(tasks[key]["shortname"])
                soff=0
                if legname not in satlist:
                    soff=ip #pindex.shape[0]
                else:
                    soff=0
                if 'andep' in legname:
                    legname=legname.split('_andep')[0]
                if legname not in [ "rio","rit", 'ce20c' ] or shade==0:
                    for iens in tasks[key]["ens"]: #range(len(tasks[key]["ens"])):

                        linesplustext=linesplustext+addsatline(tasks[key],glats,iens,soff,pkey='zslopes') 
                        legend_user_lines=legend_user_lines+[legname]
                else:
                    linesplustext=linesplustext+add_ens_satline(tasks[key], glats,soff,pkey='zslopes')
                    legend_user_lines=legend_user_lines+[legname+"{0:0>2}".format(iens)]

            if ip<pindex.shape[0]:
                lines =["Zonal Mean Temperature Trends [K/10a], "+suff+"hPa, {:4}".format(interval[0])+"-{:4}".format(interval[1])]
            else:
                lines =["Zonal Mean Brightness Temperature Trends [K/10a], "+suff+", {:4}".format(interval[0])+"-{:4}".format(interval[1])]

            projection,horizontal,vertical,title,legend=set_zonalsat(lines,legend_user_lines,plotproperties,beltinterval=[-2.0,0.5])
            out.args['output_title']=os.getcwd()+'/'+out.args['output_name']	
            plot(out,ppage,projection,horizontal,vertical,linesplustext,title,legend)
            print((out.args['output_name']+'.'+plotproperties["outputformat"][0]))


    ti=prepare_ti(tasks,rlist=['rio','rit'])
    sps=[]
    jpindex=np.concatenate((pindex,msupindex))
    for ip in pindex:
        sps.append('{:4} hPa'.format(int(ps[ip])))
    for p in msunames:
        sps.append(p)
    if 'beltanomalies' in plotproperties['plotlist']:
        if plotproperties['set_ref_interval']=='True':
            sxl=35.7
        else:
            sxl = 29.7
        ppage = page(
            layout='positional',  
            super_page_x_length=sxl, 
            super_page_y_length=plotproperties['monthlyseriesheight']+3.0, 
            page_x_length=24.7, 
            page_y_length=plotproperties['monthlyseriesheight']+3.0, 
            page_id_line='off',
            page_x_position=0., 
            page_y_position=0.)
        for ref in plotproperties['reference']: #['tmcorr','zero']:
            iref=-1
            for t in range(len(tasks)):
                if tasks[t]['shortname']==ref:
                    iref=t
            for ibelt in range(len(beltnames)):
                if beltnames[ibelt]!=str(plotproperties['belt'][0]):
                    continue
                for ip in range(len(jpindex)):


                    linesplustext=[]
                    lraw=[]
                    traw=[]
                    legend_user_lines=[]
                    if ref=='zero':
                        lines =["Temperature Anomalies [K], "+beltnames[ibelt]+', '+sps[ip]+', '+
                                "({:4}".format(interval[0])+"-{:4} Clim), ".format(interval[1])+str('Version '+plotproperties['exps'][0])]
                    else:
                        lines =["Anomaly Differences [K] to "+tasks[iref]['shortname']+', '+beltnames[ibelt]+', '+sps[ip]+', '+
                                "({:4}".format(interval[0])+"-{:4} Clim), ".format(interval[1])+str('Version '+plotproperties['exps'][0])]

                    out = output({"output_name_first_page_number":"off",
                                  "output_formats":plotproperties["outputformat"], 
                                  'output_name': plotproperties['tmppath']+'/'+'beltanomalies_'+ref+'_'+beltnames[ibelt]+'_'+''.join(sps[ip].split(' hPa')).strip()+'_'+
                                  "{:4}".format(interval[0])+"-{:4}".format(interval[1])})

                    dp=len(jpindex)-len(msupindex)    
                    keys=[]
                    for key in ti:
                        if nodata[key] or key==iref or (
                            tasks[key]['shortname'] in satlist and tasks[key]['shortname'] not in reanlist and ip<dp):
                            continue

                        shade=1
                        legname=copy.copy(tasks[key]["shortname"])
                        if 'andep' in legname:
                            legname=legname.split('_andep')[0]
                        if 'fgdep' in legname:
                            legname=legname.split('_fgdep')[0]
                        if 'erai_fg' in legname and 'gps' in legname:
                            legname=legname.split('erai_fg')[1].split('dep')[0]
                        if legname=='tm':
                            legname='unadj'
                        if legname=='tmcorr':
                            legname='adj'
                        if legname=='sunyhom':
                            legname='sunyh'


                        if legname not in [ "rio","rit"] or shade==0:
                            for iens in tasks[key]["ens"]: #range(len(tasks[key]["ens"])):

#				if tasks[key]['shortname'] in satlist:
#				    linesplustext=linesplustext+add_mts(tasks,key,ip-dp,iens,ibelt,iref,scol=tasks[key]["color"]) # RASO
#				else:

                                hh=add_mts(tasks,key,ip,iens,ibelt,iref,plotproperties,scol=tasks[key]["color"])
                                linesplustext=linesplustext+hh[0]# RASO
                                lraw.append(hh[1])
                                traw.append(hh[2])
                                keys.append(key)				
                                if len(tasks[key]['ens'])>1:
                                    legend_user_lines=legend_user_lines+[legname+"{0:0>2}".format(iens)]
                                else:
                                    legend_user_lines=legend_user_lines+[legname]
                        else:
                            legend_user_lines=legend_user_lines+[legname]
                            hhilf=np.asarray([0])
                    dynrange=np.asarray(plotproperties["range"][:])/6.*1.05
                    gtdynrange=[np.asarray([-0.2,0.2])]
                    ggsls=[[]]
                    gpint=[plotproperties['plotinterval']]
                    if plotproperties['set_ref_interval'] == 'True':    
                        ggsls.append([])
                        gpint.append(plotproperties['trend_ref_interval'])
                        gtdynrange.append(np.asarray([-0.2,0.2]))
                    for il in range(len(lraw)):
                        print('here')
                        mask=np.abs(linesplustext[2*il].args['Input_y_values'])<(dynrange[1]-dynrange[0])*20.
                        for gsls,pint, tdynrange in zip(ggsls, gpint, gtdynrange):
#                        pint=plotproperties['plotinterval']
                            tmask=(traw[il]>=pint[0])*(traw[il]<pint[1]+1)*(~np.isnan(np.abs(lraw[il])))
                            if np.sum(tmask)>np.sum((traw[il]>=pint[0])&(traw[il]<pint[1]+1)) - 36:
                                print((np.std(traw[il][tmask]-1900.),np.nanstd(lraw[il][tmask])))
                                #st=linregress(traw[il][tmask]-1900.,lraw[il][tmask])
                            #st=np.polyfit(traw[il][tmask]-1900.,lraw[il][tmask], 1, cov=True)
                                st=fastlinregressnonanerr(traw[il][tmask]-1900.,lraw[il][tmask])
                                print((tasks[keys[il]]['shortname'] , np.sum(tmask), np.sum((traw[il]>=pint[0])&(traw[il]<pint[1]+1)), st))
                                alpha=alphaest(lraw[il][tmask]-np.mean(lraw[il][tmask]))
                                if (tasks[keys[il]]['shortname'] in ['tm','tmcorr'] or 'rio'  in tasks[keys[il]]['shortname']  or 'rit' in  tasks[keys[il]]['shortname'])  and ref!='zero':
                                    alpha=0.
                                neffcorr=np.sqrt((1+alpha)/(1-alpha))
                                gsls.append([])
                                gsls[-1].append(st[0]*10)
                                err=st[1]*1.96*neffcorr*10
                                gsls[-1].append(gsls[-1][0]-err)
                                gsls[-1].append(gsls[-1][0]+err)
    #gsls.append(fastlinregressnonan(l.args['Input_x_values'][tmask]-1900.,
    #					                l.args['Input_y_values'][tmask])*10)
                                if gsls[-1][2]>tdynrange[1]:
                                    tdynrange[1]=gsls[-1][2]*1.05
                                if gsls[-1][1]<tdynrange[0]:
                                    tdynrange[0]=gsls[-1][1]*1.05
    
                            else:
                                gsls.append([-1.e21]*3)
                            if sum(mask)>0:
                                m=np.amax(linesplustext[2*il].args['Input_y_values'][mask])
                                if m>dynrange[1]:
                                    dynrange[1]=m*1.05
                                m=np.amin(linesplustext[2*il].args['Input_y_values'][mask])
                                if m<dynrange[0]:
                                    dynrange[0]=m*1.05
                                print(tasks[keys[il]]['shortname'], dynrange)


                    print('here')
                    projection,horizontal,vertical,title,legend,null_data,null_line=monthlyseries(lines,legend_user_lines,plotproperties,dynrange=dynrange)
                    print('there')
#		    hadtrend=addhadtrend(hadmedbeltslopes[ibelt,0],hadbeltslopes[:,ibelt,0])
#		    legend_user_lines=legend_user_lines+["HadCRUT4"]
#                    tline=['Trends {}-{}'.format(plotproperties['plotinterval'][0],plotproperties['plotinterval'][1])]
                    tline=['Trends [K/10a]', '{}-{}'.format(plotproperties['plotinterval'][0],plotproperties['plotinterval'][1])]
                    ppage2,projection2,horizontal2,vertical2,title2,null_data2,null_line2=sbtrends(tline,plotproperties,dynrange=gtdynrange[0])
                    markers=sbmarkers(tasks,keys,ggsls[0])

                    if(plotproperties['set_ref_interval'] == 'True'):
                        
                        tline=['Trends [K/10a]', '{}-{}'.format(*plotproperties['trend_ref_interval'])]
                        ppage3,projection3,horizontal3,vertical3,title3,null_data3,null_line3=sbtrends(tline,plotproperties,dynrange=gtdynrange[1], second=True)
                        #markers=[]
                        #if len(gsls)>0:
                        markers3=sbmarkers(tasks,keys,ggsls[1])
                    out.args['output_title']=os.getcwd()+'/'+out.args['output_name']
                    try:
                        igood=0
                        for it in range(len(linesplustext)-2,-1,-2):
                            linesplustext[it].args['Input_x_values'][linesplustext[it].args['Input_x_values']>10000.]=-1.e21
                            linesplustext[it].args['Input_y_values'][linesplustext[it].args['Input_y_values']>10000.]=-1.e21
                            linesplustext[it].args['Input_x_values'][linesplustext[it].args['Input_x_values']<-10000.]=-1.e21
                            linesplustext[it].args['Input_y_values'][linesplustext[it].args['Input_y_values']<-10000.]=-1.e21
                            if any(linesplustext[it].args['Input_y_values']!=-1.e21):
                                #print(tasks[it//2]['shortname'],'good values')
                                igood +=1
                            else:
                                linesplustext.pop(it)
                                linesplustext.pop(it)
                                
                        if igood>0: #<len(linesplustext)//2:
                            #print(tasks[it//2]['shortname'],linesplustext[it].args['Input_y_values'],linesplustext[it].args['Input_y_values'])
                            cmd = [out,ppage,projection,horizontal,vertical,linesplustext,title,null_data,null_line,legend,
                                 ppage2,projection2,horizontal2,vertical2,title2,null_data2,null_line2,markers]
                            if plotproperties['set_ref_interval'] == 'True':
                                cmd += [ppage3,projection3,horizontal3,vertical3,title3,null_data3,null_line3,markers3]
                            plot(*cmd)
                        else:
                            print('NO PLOT PRODUCED', out.args['output_name']+'.'+plotproperties["outputformat"][0])
                            
                    except:
                        
                        print('NO PLOT PRODUCED', out.args['output_name']+'.'+plotproperties["outputformat"][0])
                    print((out.args['output_name']+'.'+plotproperties["outputformat"][0]))


    if 'belts' in plotproperties['plotlist']:
        hadmedbeltslopes=np.zeros([1,belts.shape[0],1])
        hadbeltslopes=np.zeros([hadens.shape[0],belts.shape[0],1])
        try:
            hadmedbeltslopes=belttrends(hadmedzslopes,belts)
        except IOError:
            print('no hadmedzslopes, returning ...')
            return hadmed,hadtem,hadens,sats
        for jens in range(hadens.shape[0]):
            hadbeltslopes[jens,:,:]=belttrends(np.reshape(hadzslopes[jens,:,:],[nj,1]),belts)

        for ibelt in range(len(beltnames)):

            linesplustext=[]
            legend_user_lines=[]
            lines =["Temperature Trends [K/10a], "+beltnames[ibelt]+', '+
                    "{:4}".format(interval[0])+"-{:4}".format(interval[1])+', '+
                    str(plotproperties['exps'][0])]

            profsat=''
            if 'satseries' in plotproperties['plotlist']:
                profsat='sat'

            out = output({"output_name_first_page_number":"off",
                          "output_formats":plotproperties["outputformat"], 
                          'output_name': plotproperties['tmppath']+'/'+profsat+'trendsbelt_'+beltnames[ibelt]+'_'+
                          "{:4}".format(interval[0])+"-{:4}".format(interval[1])})

            amplevels=np.where(np.logical_and(ps[pindex]>150,ps[pindex]<=400))[0]
            print(('amplevels:',amplevels))
            hadbeltq=np.quantile(hadbeltslopes[:,ibelt,0],[0.05,0.5,0.95])
            for key in ti:
                if nodata[key]:
                    continue
                if tasks[key]['shortname'] in satlist and 'had' not in tasks[key]['shortname'] and tasks[key]['shortname'] not in reanlist and start<1979:
                    continue
                #if tasks[key]['shortname'] in ['tmcorr']:
                    #continue
                if np.nansum(tasks[key]["zslopes"])==0:
                    continue

                for iens in tasks[key]["ens"]: #range(len(tasks[key]["ens"])):
                    if ibelt==0:
                        tasks[key]["beltslopes"][iens,:,:]=belttrends(tasks[key]["zslopes"][iens,:,:],belts)
                shade=1
                legname=copy.copy(tasks[key]["shortname"])
                if 'andep' in legname:
                    legname=legname.split('_andep')[0]
                if legname not in [ "rio","rit", 'ce20c' ] or shade==0:
                    for iens in tasks[key]["ens"]: #range(len(tasks[key]["ens"])):

                        if 'satseries' in plotproperties['plotlist']:
                            if legname not in satlist:
                                soff=plotproperties['msupindex'].shape[0]
                                tasks[key]["beltslopes"][iens,ibelt,-soff]=np.nan
                            linesplustext=linesplustext+addsattrend(tasks[key]["beltslopes"][iens,ibelt,-soff:],msups,
                                                                    scol=tasks[key]["color"],iens=iens) # SAT
                            amp=tasks[key]["beltslopes"][iens,ibelt,-soff+1]/hadbeltq
                        else:
                            if tasks[key]["beltslopes"].shape[2]>=pindex.shape[0]:
                                linesplustext=linesplustext+addprofile(tasks[key],ps[pindex],iens,ibelt) # RASO
                            else:
                                linesplustext=linesplustext+addsattrend(tasks[key]["beltslopes"][iens,ibelt,:],msups[msupindex],
                                                                        scol=tasks[key]["color"],iens=iens) # SAT
                            try:
                                if tasks[key]['shortname'] in ['uah','rss','star']:
                                    ampl=np.where(msupindex==1)[0][0]
                                else:
                                    ampl=amplevels.copy()
                                trmax=np.max(tasks[key]["beltslopes"][iens,ibelt,ampl].flatten())
                                #amp=trmax/hadmedbeltslopes[ibelt,0]
                                amp=trmax/hadbeltq
                            except:
                                amp=tasks[key]["beltslopes"][iens,ibelt,1]/hadbeltq
                        amps=" {0:3.1f}".format(amp[1])
                        if len(tasks[key]['ens'])>1:
                            legend_user_lines=legend_user_lines+[legname+"{0:0>2}".format(iens)+amps]
                        else:
                            legend_user_lines=legend_user_lines+[legname+amps]
                    tasks[key]['amps'][:,0,interv,ibelt]=amp
                else:
                    amp=np.zeros((tasks[key]['amps'].shape[0],tasks[key]["ens"].shape[0]))
                    if 'satseries' in plotproperties['plotlist']:
                        for iens in tasks[key]["ens"]: #range(len(tasks[key]["ens"])):
                            tasks[key]["beltslopes"][iens,ibelt,-4]=np.nan
                            amp[:,iens]=tasks[key]["beltslopes"][iens,ibelt,-3]/hadbeltq
                            linesplustext=linesplustext+addsattrend(tasks[key]["beltslopes"][iens,ibelt,-4:],[800.,550.,250.,90.],
                                                                    scol=tasks[key]["color"],iens=iens) # SAT
                            print((tasks[key]["beltslopes"][iens,ibelt,-4:]))
                    else:
                        linesplustext=linesplustext+add_ens_profile(tasks[key], ps[pindex], ibelt)
                        for iens in tasks[key]["ens"]: #range(len(tasks[key]["ens"])):
                            try:
                                if tasks[key]['shortname'] in ['uah','rss','star']:
                                    ampl=[1]
                                else:
                                    ampl=amplevels.copy()
                                amp[:,iens]=np.max(tasks[key]["beltslopes"][iens,ibelt,ampl])/hadbeltq
                            except:
                                amp[:,iens]=np.nan
                    amps=" {0:3.1f}".format(np.mean(amp,axis=1)[1])
                    tasks[key]['amps'][:,:,interv,ibelt]=amp[:]
                    legend_user_lines=legend_user_lines+[legname+amps]

                tasks[key]['had'][:,interv,ibelt]=np.quantile(hadbeltslopes[:,ibelt,0],[0.05,0.5,0.95])

            projection,horizontal,vertical,title,legend=set_belttrends(lines,legend_user_lines,plotproperties)
            hadtrend=addhadtrend(hadmedbeltslopes[ibelt,0],hadbeltslopes[:,ibelt,0])
            legend_user_lines=legend_user_lines+["HadCRUT5"]

            if len(linesplustext)>0:
                out.args['output_title']=os.getcwd()+'/'+out.args['output_name']
                ldat=False
                for il in range(0,len(linesplustext),2):
                    if len(linesplustext[il].args['Input_x_values'])>0:
                        ldat=True
                if ldat:
                    plot(out,ppage,projection,horizontal,vertical,linesplustext,hadtrend,title,legend)
                    print((out.args['output_name']+'.'+plotproperties["outputformat"][0]))
                else:
                    print((out.args['output_name']+'.'+plotproperties["outputformat"][0]+ 'empty linesplustext'))		    
            else:
                print((out.args['output_name']+'.'+plotproperties["outputformat"][0]+ 'not created - no data found'))


    if 'stations' in plotproperties['plotlist'] or 'cost' in plotproperties['plotlist']:
        makecostprofiles(tasks,nodata,satlist,interval,plotproperties,'trend',ps,pindex,iens,msups,msupindex)




    return hadmed,hadtem,hadens,sats

def makecostprofiles(tasks,nodata,satlist,interval,plotproperties,what,ps,pindex,iens,msups,msupindex):

    if what=='trend':
        whatlong='Trend'
    else:
        whatlong='Difference'
    out = output({"output_name_first_page_number":"off",
                  "output_formats":plotproperties["outputformat"], 
                  'output_name': plotproperties['tmppath']+'/'+what+'cost'+'_'+
                  "{:4}".format(interval[0])+"-{:4}".format(interval[1])})

    linesplustext=[]	
    legend_user_lines=[]
    lines =["Temperature "+whatlong+" Cost, "+
            "{:4}".format(interval[0])+"-{:4}".format(interval[1]),
            str(plotproperties['version']+', '+os.getcwd().split('/')[-1]+', '+plotproperties["fgdepname"])]
    ti=prepare_ti(tasks)

    shade=1
    for key in ti:
        if nodata[key]:
            continue
        legname=copy.copy(tasks[key]["shortname"])
        if 'andep' in legname:
            legname=legname.split('_andep')[0]
        if 'fgdep' in legname:
            legname=legname.split('_fgdep')[0]
        if legname not in [ "rio", 'rit', 'ce20c' ] or shade==0:
            if legname in satlist:
                linesplustext=linesplustext+addprofile(tasks[key],msups[msupindex],iens,0,pkey='cost') # RASO
            else:
                for iens in tasks[key]["ens"]: #range(len(tasks[key]["ens"])):

                    linesplustext=linesplustext+addprofile(tasks[key],ps[pindex],iens,2,pkey='cost') # RASO
        else:
            linesplustext=linesplustext+add_ens_profile(tasks[key], ps[pindex],2,pkey='cost')

        if linesplustext[-2].args['Input_x_values'].shape[0]<2:
            del linesplustext[-1]
            del linesplustext[-1]
            continue
        n=np.sum(~np.isnan(tasks[key]['cost'][:,2,:].flatten()))
        legname=legname+' {:d}'.format(np.int(np.sqrt(np.nansum(tasks[key]['cost'][:,2,:].flatten()**2/n))))
        if len(tasks[key]["ens"])>1:
            legend_user_lines=legend_user_lines+[legname+"{0:0>2}".format(iens)]
        else:
            legend_user_lines=legend_user_lines+[legname]

    maxcost=np.NaN
    for ta in tasks:
        maxcost=np.nanmax([maxcost,np.nanmax(ta['cost'][:,2,:].flatten())])
    if np.isnan(maxcost):
        print('All cost profiles are NaN, continuing')
    else:


        ppage = page(
            layout='positional',  
            page_x_length=29.7, 
            page_y_length=21., 
            page_id_line='off',
            page_x_position=0., 
            page_y_position=0.)
        projection,horizontal,vertical,title,legend=set_belttrends(lines,legend_user_lines,plotproperties,pkey='cost',beltinterval=[0.,float(maxcost)])

        out.args['output_title']=os.getcwd()+'/'+out.args['output_name']	
        if len(linesplustext)>0:
            plot(out,ppage,projection,horizontal,vertical,linesplustext,title,legend)
        else:
            plot(out,ppage,projection,horizontal,vertical,linesplustext,title) #,title,legend
        print((out.args['output_name']+'.'+plotproperties["outputformat"][0]))

def save_gridded(ganomalies,gclimatologies,tasks,key,days,ps,pindex,iens=0,start=2006,stop=2015,version='1.9',append=True):

    t=time.time()
    for im in range(ganomalies.shape[4]-25,ganomalies.shape[4]):
        f=plt.figure(figsize=(20,13))
        t1=time.time()
        #m =Basemap(llcrnrlon=-180.,llcrnrlat=-90.,urcrnrlon=180,urcrnrlat=90.)
        #m.drawmapboundary(fill_color=rgb(0.6,0.8,1))
        #parallels = (np.arange(-90.,91,90.))
        ## labels = [left,right,top,bottom]
        #m.drawparallels(parallels,labels=[True,True,True,True])
        #meridians = np.arange(-180.,181.,60.)
        #m.drawmeridians(meridians,labels=[True,True,True,True])
        #m.drawcoastlines(linewidth=0.25)
#        ax2 = f.add_axes([0.1,0.1,0.,0.8],projection=ccrs.PlateCarree())
        ax2=plt.axes(projection=ccrs.PlateCarree())
        ax2.coastlines()
        ax2.gridlines()
        ax2.plot([(-180.,180.),(-90.,90.)],'or',transform=ccrs.PlateCarree())

        ni=ganomalies.shape[1]
        nj=ganomalies.shape[0]
        glons=5.+10.*np.arange(ni)
        glats=-85.+10.*np.arange(nj)
        gglons=glons.copy()
        gglons[gglons>180]-=360
#        x, y = m(*np.meshgrid(gglons,glats))
        x, y = np.meshgrid(gglons,glats)

        clevs=(np.arange(21)-10)
        cmap=cm.get_cmap(name='hsv')
        cNorm  = colors.Normalize(vmin=clevs[0], vmax=clevs[-1])
        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmap)
        for ilon in range(glons.shape[0]):
            for ilat in range(glats.shape[0]):
                if ganomalies[ilat,ilon,0,3,im]==ganomalies[ilat,ilon,0,3,im] or ganomalies[ilat,ilon,1,3,im]==ganomalies[ilat,ilon,1,3,im]:
                    xy=list(zip([x[ilat,ilon]-5,x[ilat,ilon]+5,x[ilat,ilon]+5,x[ilat,ilon]-5,x[ilat,ilon]-5],
                                [y[ilat,ilon]-5,y[ilat,ilon]-5,y[ilat,ilon]+5,y[ilat,ilon]+5,y[ilat,ilon]-5]))
                    cl=clevs[clevs>np.nanmean(ganomalies[ilat,ilon,:,3,im])]
                    if cl.shape[0]==0:
                        cl=clevs.shape[0]-1
                    else:
                        cl=cl[0]
                    colorVal = scalarMap.to_rgba(cl)
                    poly = mpatches.Polygon( xy, facecolor=colorVal, edgecolor=[0.8,0.8,0.8],lw=0.3,transform=ccrs.PlateCarree())
                    plt.gca().add_patch(poly)
    #		m.contourf(x,y, ganomalies[:,:,0,3,im])
        ax=f.add_subplot(2,1,2)
        ax.set_position([0.13,0.87,0.77,0.04])
        ax.set_xlim([clevs[0],clevs[-1]])
    #		    ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        legname=copy.copy(tasks[key]["shortname"])
        ax.set_title(legname+', {0}{1:2>0}'.format(1900+im//12,np.mod(im,12)+1))
        for ipos in range(clevs.shape[0]-1):
            xy=list(zip([clevs[ipos],clevs[ipos+1],clevs[ipos+1],clevs[ipos],clevs[ipos]],[0,0,1,1,0]))
            colorVal = scalarMap.to_rgba(clevs[ipos])
            poly=mpatches.Polygon(xy, facecolor=colorVal, edgecolor=[0.8,0.8,0.8],lw=0.3)
            plt.gca().add_patch(poly)

        print((time.time()-t))
        f.savefig('test_'+tasks[key]["shortname"]+'{0:4}'.format(im)+'.eps')
        plt.close(f)
        print((time.time()-t))

    ipath=os.path.expandvars('$FSCRATCH/classic/sabiner/exp07/')
    #if tasks[key]["shortname"]=='tmcorr':
        #fn='raobcore15_gridded_20{}.nc'.format(ganomalies.shape[4]/12-101)
        #append_gridded(ipath+'raobcore15_gridded_2013.nc',fn,ganomalies,days,ps[pindex.astype(int)])
    #if tasks[key]["shortname"]=='tm':
        #fn='raobcore15_raw_gridded_20{}.nc'.format(ganomalies.shape[4]/12-101)
        #append_gridded(ipath+'raobcore_raw_gridded_2013.nc',fn,ganomalies,days,ps[pindex.astype(int)])
    #if tasks[key]["shortname"]=='rio24':
        #fn='rich15obs_mean_gridded_20{}.nc'.format(ganomalies.shape[4]/12-101)
        #append_gridded(ipath+'rich15obs_mean_gridded_2013.nc',fn,ganomalies,days,ps[pindex.astype(int)])

    vs=''.join(version.split('.'))
    if tasks[key]["shortname"]=='era5v7':

        fn='rase'+vs+'_gridded_20{}.nc'.format(ganomalies.shape[4]//12-101)
        append_gridded_10(ipath+'raobcore15_gridded_2013.nc',fn,ganomalies,days,ps[pindex.astype(int)],start=start,stop=stop,append=False)
    if tasks[key]["shortname"]=='era5v7rich':
        fn='rise'+vs+'_gridded_20{}.nc'.format(ganomalies.shape[4]//12-101)
        append_gridded_10(ipath+'raobcore15_gridded_2013.nc',fn,ganomalies,days,ps[pindex.astype(int)],start=start,stop=stop,append=False)
    if tasks[key]["shortname"]=='tmcorr':
        fn='raobcore'+vs+'_gridded_20{}.nc'.format(ganomalies.shape[4]//12-101)
        append_gridded_10(ipath+'raobcore15_gridded_2013.nc',fn,ganomalies,days,ps[pindex.astype(int)],start=start,stop=stop)
    if tasks[key]["shortname"]=='tm':
        fn='raobcore'+vs+'_raw_gridded_20{}.nc'.format(ganomalies.shape[4]//12-101)
        append_gridded_10(ipath+'raobcore_raw_gridded_2013.nc',fn,ganomalies,days,ps[pindex.astype(int)],start=start,stop=stop)
    if tasks[key]["shortname"]=='rio24':
        fn='rich'+vs+'obs_mean_gridded_20{}.nc'.format(ganomalies.shape[4]//12-101)
        append_gridded_10(ipath+'rich15obs_mean_gridded_2013.nc',fn,ganomalies,days,ps[pindex.astype(int)],start=start,stop=stop)

    if append==False:
        if tasks[key]["shortname"] in ['rio','rit']:
            fn=tasks[key]["shortname"]+'{:0>2}_gridded_20{}.nc'.format(iens,ganomalies.shape[4]//12-101)
        else:
            fn=tasks[key]["shortname"]+'_gridded_20{}.nc'.format(ganomalies.shape[4]//12-101)
            
        append_gridded_10(ipath+'raobcore15_gridded_2013.nc',fn,ganomalies,days,ps[pindex.astype(int)],start=start,stop=stop,append=append)


#    try:
    plot_gridded(fn)
    plot_gridded_trend(fn,[start,stop])
#    except:
#	pass

    if tasks[key]["shortname"]=='era5v7':
        fn='rase'+vs+'_gridded_clim_{}-{}.nc'.format(start,stop)
        append_gridded_10(ipath+'raobcore15_gridded_2013.nc',fn,gclimatologies,days[:12],ps[pindex.astype(int)],start=start,stop=stop,append=False)
    if tasks[key]["shortname"]=='era5v7rich':
        fn='rise'+vs+'_gridded_clim_{}-{}.nc'.format(start,stop)
        append_gridded_10(ipath+'raobcore15_gridded_2013.nc',fn,gclimatologies,days[:12],ps[pindex.astype(int)],start=start,stop=stop,append=False)
    if tasks[key]["shortname"]=='tmcorr':
        fn='raobcore'+vs+'_gridded_clim_{}-{}.nc'.format(start,stop)
        append_gridded_10(ipath+'raobcore15_gridded_2013.nc',fn,gclimatologies,days[:12],ps[pindex.astype(int)],start=start,stop=stop)
    if tasks[key]["shortname"]=='tm':
        fn='raobcore'+vs+'_raw_gridded_clim_{}-{}.nc'.format(start,stop)
        append_gridded_10(ipath+'raobcore_raw_gridded_2013.nc',fn,gclimatologies,days[:12],ps[pindex.astype(int)],start=start,stop=stop)
    if tasks[key]["shortname"]=='rio24':
        fn='rich'+vs+'obs_mean_gridded__clim_{}-{}.nc'.format(start,stop)
        append_gridded_10(ipath+'rich15obs_mean_gridded_2013.nc',fn,gclimatologies,days[:12],ps[pindex.astype(int)],start=start,stop=stop)

    try:
        plot_gridded(fn)
    except:
        pass

    return
    #		if tasks[key]["shortname"]=='tmcorr':
    #		    append_gridded('../exp07/rich15tau_mean_gridded_2013.nc','rich15tau_mean_gridded_2014.nc',
    #		                   ganomalies)


def plot_gridded(fn):
    t=time.time()
    f = netCDF4.Dataset(fn,"r")
    f.set_auto_mask(False)

    ganomalies=f.variables['anomalies'][:]
    print((f.variables['pressure']))
    gs=ganomalies.shape[0]
    for im in range(gs-13,gs,12):
        fig=plt.figure(figsize=(20,13))
        t1=time.time()
        #m =Basemap(llcrnrlon=-180.,llcrnrlat=-90.,urcrnrlon=180,urcrnrlat=90.)
        #m.drawmapboundary(fill_color=rgb(0.6,0.8,1))
        #parallels = (np.arange(-90.,91,90.))
        ## labels = [left,right,top,bottom]
        #m.drawparallels(parallels,labels=[True,True,True,True])
        #meridians = np.arange(-180.,181.,60.)
        #m.drawmeridians(meridians,labels=[True,True,True,True])
        #m.drawcoastlines(linewidth=0.25)
        ax2=plt.axes(projection=ccrs.PlateCarree())
        ax2.coastlines()
        ax2.gridlines()
        ax2.plot([(-180.,180.),(-90.,90.)],'or',transform=ccrs.PlateCarree())

        ni=ganomalies.shape[3]
        nj=ganomalies.shape[2]
        glons=f.variables['lon'][:]#5.+10.*np.arange(ni)
        glats=f.variables['lat'][:] #-85.+10.*np.arange(nj)
        gglons=glons.copy()
        gglons[gglons>180]-=360
        x, y = np.meshgrid(gglons,glats)

        clevs=(np.arange(21)-10)
        cmap=cm.get_cmap(name='jet')
        cNorm  = colors.Normalize(vmin=clevs[0], vmax=clevs[-1])
        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmap)
        for ilon in range(glons.shape[0]):
            for ilat in range(glats.shape[0]):
                if ganomalies[im,8,ilat,ilon]==ganomalies[im,8,ilat,ilon]:
                    xy=list(zip([x[ilat,ilon]-5,x[ilat,ilon]+5,x[ilat,ilon]+5,x[ilat,ilon]-5,x[ilat,ilon]-5],
                                [y[ilat,ilon]-5,y[ilat,ilon]-5,y[ilat,ilon]+5,y[ilat,ilon]+5,y[ilat,ilon]-5]))
                    cl=clevs[clevs>np.nanmean(ganomalies[im,8,ilat,ilon])]
                    if cl.shape[0]==0:
                        cl=clevs.shape[0]-1
                    else:
                        cl=cl[0]
                    colorVal = scalarMap.to_rgba(cl)
                    poly = mpatches.Polygon( xy, facecolor=colorVal, edgecolor=[0.8,0.8,0.8],lw=0.3,transform=ccrs.PlateCarree())
                    plt.gca().add_patch(poly)
    #		m.contourf(x,y, ganomalies[:,:,0,3,im])
        ax=fig.add_subplot(2,1,2)
        ax.set_position([0.13,0.87,0.77,0.04])
        ax.set_xlim([clevs[0],clevs[-1]])
    #		    ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.set_title(fn[:-3]+', {0}{1:2>0}'.format(1900+im//12,np.mod(im,12)+1))
        for ipos in range(clevs.shape[0]-1):
            xy=list(zip([clevs[ipos],clevs[ipos+1],clevs[ipos+1],clevs[ipos],clevs[ipos]],[0,0,1,1,0]))
            colorVal = scalarMap.to_rgba(clevs[ipos])
            poly=mpatches.Polygon(xy, facecolor=colorVal, edgecolor=[0.8,0.8,0.8],lw=0.3)
            plt.gca().add_patch(poly)

        print((time.time()-t))
        fig.savefig('test2_'+fn[:-3]+'{0:4}'.format(im)+'.eps')
        plt.close(fig)
        print((time.time()-t))
    return

def plot_gridded_trend(fn,interval):
    t=time.time()
    f = netCDF4.Dataset(fn,"r")
    f.set_auto_mask(False)

    ganomalies=f.variables['anomalies'][:]
    ganomalies[ganomalies==-1.e30]=np.nan
    ps=f.variables['pressure'][:]
    lats=f.variables['lat'][:]
    print(ps)
    startmonth=(interval[0]-1900)*12
    endmonth=(interval[1]-1900+1)*12

    iseries=np.arange(startmonth,endmonth)/120.
    gzmean=np.empty((ganomalies.shape[2],ganomalies.shape[1]))
    zmseries=iseries-iseries
    for ip in range(ganomalies.shape[1]):
        for j in range(ganomalies.shape[2]):
            for it in range(iseries.shape[0]):
                if it+startmonth>=ganomalies.shape[0]:
                    zmseries[it]=np.nan
                else:
                    zmseries[it]=np.nanmean(ganomalies[it+startmonth,ip,j,:])
            gzmean[j,ip]=fastlinregress(iseries, zmseries)

    fig=plt.figure(figsize=(20,13))
    clevs=(np.arange(21)-10)/10.
    cmap=cm.get_cmap(name='jet')
    cNorm  = colors.Normalize(vmin=clevs[0], vmax=clevs[-1])
    scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmap)
    plt.contourf(lats,ps,gzmean.T,clevs=clevs)	
    ax=plt.gca()
    ax.set_yscale('log')
    ax.set_ylim(850,20)
    ax.set_yticks((850,700,500,300,200,100,50,30,20))
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

    ax=fig.add_subplot(2,1,2)
    ax.set_position([0.13,0.91,0.77,0.04])
    ax.set_xlim([clevs[0],clevs[-1]])
    ax.get_yaxis().set_visible(False)
    ax.set_title(fn[:-3]+', {}-{}'.format(interval[0],interval[1]))
    for ipos in range(clevs.shape[0]-1):
        xy=list(zip([clevs[ipos],clevs[ipos+1],clevs[ipos+1],clevs[ipos],clevs[ipos]],[0,0,1,1,0]))
        colorVal = scalarMap.to_rgba(clevs[ipos])
        poly=mpatches.Polygon(xy, facecolor=colorVal, edgecolor=[0.8,0.8,0.8],lw=0.3)
        plt.gca().add_patch(poly)
    

    fig.savefig('ztrend_'+fn[:-3]+'{}-{}'.format(interval[0],interval[1])+'.eps')
    plt.close(fig)
    print((time.time()-t))
