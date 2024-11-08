#!/usr/bin/env python
import sys
import os.path
import glob

import subprocess
import urllib.request
import xarray as xr
import numpy
from datetime import date, datetime
import time
from multiprocessing import Pool
from netCDF4 import Dataset
import gzip
import pandas as pd
from functools import partial
from eccodes import *
import matplotlib.pylab as plt
#import cartopy.crs as ccrs
import argparse

pd.options.mode.chained_assignment = None

import code
#import path
#sys.path.append('/fio/srvx7/leo/python')
sys.path.append("/fio/srvx7/leo/python/Rasotools/")
from utils import *                                                                                                                                                

def now(time):
      """ For easily print a readable current time stamp """
      a = datetime.datetime.fromtimestamp(time).strftime('%Y-%m-%d %H:%M:%S')
      return  a

home= os.getcwd()
def plot_active(archives,period,path=''):

    if path=='':
        path=os.path.expandvars('$RSCRATCH/era5/odbs/')

    meta={}
    idx= {}

    f=plt.figure(figsize=(10,5))
    ax=plt.subplot(1,1,1,projection=ccrs.PlateCarree())
    for a in archives:
        meta[a]=pd.read_csv(path+a+'/meta.csv',sep='\t')

        #idx[a]=numpy.where(meta[a].start_date <'1950')

        try:
            idx[a]=numpy.where(numpy.logical_and.reduce((meta[a].start_date <period[1],meta[a].end_date >=period[0],
                                                     meta[a].longitude>-180,meta[a].latitude>=-90)))
        except:
            idx[a]=[]
            continue

        c=plt.scatter(meta[a]['longitude'].iloc[idx[a]],meta[a]['latitude'].iloc[idx[a]],transform=ccrs.PlateCarree(),s=9,label=a+' '+str(len(idx[a][0])))
        
    plt.xlim([-180.,180.])
    plt.ylim([-90.,90.])
    #clb=f.colorbar(c,ax=ax,orientation='horizontal')
    #clb.ax.set_title('{:6.4f}'.format(np.mean(ds[k]).values/43200)+' '+units)

    ax.coastlines()
    gl=ax.gridlines(draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_right = False

    if len(archives)>1:
        a='all'
        plt.legend()
    plt.title(a+', {}-{}'.format(period[0],period[1]))

    plt.savefig(path+a+'/{}-{}.png'.format(period[0],period[1]),bbox_inches='tight')
    plt.close()
    return

def plot_orphans(archives,path=''):

    if path=='':
        path=os.path.expandvars('$RSCRATCH/era5/odbs/')

    meta={}
    idx= {}

    f=plt.figure(figsize=(10,5))
    ax=plt.subplot(1,1,1,projection=ccrs.PlateCarree())
    for a in archives:
        try:
            
            meta[a]=pd.read_csv(path+a+'/orphans.csv',
                            names=['filename','identifier','data latitude','data longitude','latitude','longitude','from','to'])

        #idx[a]=numpy.where(meta[a].start_date <'1950')

            if a!='1759':
                c=plt.scatter(meta[a]['data longitude'][:],meta[a]['data latitude'][:],transform=ccrs.PlateCarree(),s=9,label=a+' '+str(len(meta[a]['data latitude'][:])))
            else:
                idx=numpy.where(abs(meta[a]['latitude'][:]+meta[a]['data latitude'][:])<0.5)
                idy=numpy.where(abs(meta[a]['latitude'][:]+meta[a]['data latitude'][:])>0.5)
                c=plt.scatter(meta[a]['data longitude'][idy[0]],meta[a]['data latitude'][idy[0]],transform=ccrs.PlateCarree(),s=9,label=a+' '+str(len(idy[0])))
                c=plt.scatter(meta[a]['longitude'][idx[0]],meta[a]['latitude'][idx[0]],transform=ccrs.PlateCarree(),s=25,marker='*',label=a+' '+str(len(idx[0])))
            
        except:
            pass
        
    plt.xlim([-180.,180.])
    plt.ylim([-90.,90.])

    ax.coastlines()
    gl=ax.gridlines(draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_right = False

    if len(archives)>1:
        a='all'
        plt.legend()
    plt.title(a+' orphans')

    plt.savefig(path+a+'/orphans.png',bbox_inches='tight')
    plt.close()
    return

def plot_nowigos(archives,path=''):

    if path=='':
        path=os.path.expandvars('$RSCRATCH/era5/odbs/')

    meta={}
    idx= {}

    f=plt.figure(figsize=(10,5))
    ax=plt.subplot(1,1,1,projection=ccrs.PlateCarree())
    for a in archives:
        try:
            
            meta[a]=pd.read_csv(path+a+'/meta.csv',sep='\t')

        #idx[a]=numpy.where(meta[a].start_date <'1950')

            idx=[]
            c=0
            for v in meta[a]['primary_id'].values:
                if '20000-0' not in v:
                    idx.append(c)
                c+=1
                    
            if len(idx)>0:
                
                pc=ax.scatter(meta[a]['longitude'][idx],meta[a]['latitude'][idx],transform=ccrs.PlateCarree(),
                          s=25,marker='*',label=a+' '+str(len(idx)))
            
        except:
            print('failed')
            pass
        
    plt.xlim([-180.,180.])
    plt.ylim([-90.,90.])

    ax.coastlines()
    gl=ax.gridlines(draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_right = False

    if len(archives)>1:
        a='all'
        plt.legend()
    plt.title(a+' stations with no WIGOS ID')

    plt.savefig(path+a+'/nowigos.png',bbox_inches='tight')
    plt.close()
    return

def plot_numbers(archives,path=''):

    if path=='':
        path=os.path.expandvars('$RSCRATCH/era5/odbs/')

    meta={}
    idx= {}

    f=plt.figure(figsize=(10,5))
    ax=plt.subplot(1,1,1)
    years=numpy.arange(1940,1980)
    for a in archives:
        try:
            
            meta[a]=pd.read_csv(path+a+'/meta.csv',sep='\t')

        #idx[a]=numpy.where(meta[a].start_date <'1950')

            nums=numpy.zeros(years.shape[0])
            for y in years:
                
                try:
                    idx[a]=numpy.where(numpy.logical_and.reduce((meta[a].start_date <str(y+1),meta[a].end_date >=str(y),
                                                             meta[a].longitude>-180,meta[a].latitude>=-90)))
                except:
                    idx[a]=[]
                    continue
                nums[y-years[0]]=len(idx[a][0])
            
            ax.step(years,nums,label=a)
        except:
            
            pass
            
    #plt.xlim([-180.,180.])
    #plt.ylim([-90.,90.])

    if len(archives)>1:
        a='all'
        plt.legend()
    plt.title('Active stations in year')

    plt.savefig(path+a+'/numbers.png',bbox_inches='tight')
    plt.close()
    return

def uai_trans(uaikeys,unified):
    trans=dict(zip(uaikeys,['']*len(uaikeys)))
    trans['ID']='station_id',
    trans['Source Link (Images)']=('source_uid','data_repository_ftp')
    trans['Source']=('source_name','originator_Organisation')
    trans['Original Country or Region']='territoryOfOriginOfData_Name',
    trans['Station Name']='station name',
    trans['Platform']='stationPlatformType',
    trans['Altitude']='elev',
    trans['Original_Latitude_Units']='stationDescription',
    trans['Original_Longitude_Units']='stationDescription',
    trans['Original_Altitude_Units']='stationDescription',
    trans['Last modified']='creationDate',
    trans['Start Station Year Month Day']='startOfOperation',
    #trans['End Station Year Month Day']='',
    trans['Start Record Year Month Day']='lot2start(vname)',
    trans['End Record Year Month Day']='lot2end(vname)',
    trans['WMO Region']='continent',
    trans['Latitude']='lat',
    trans['Longitude']='lon',
    trans['Current Station Name']='station name',
    trans['WMO_ID']='wmo_id','stationPlatformUniqueIdentifier', # WIGOS ID
    trans['Network1_name']='network_Names',
    trans['Network1_ID']='network_IDs',
    trans['Network2_name']='network_Names',
    trans['Network2_ID']='network_IDs',
    trans['Country or Region']='country',
    trans['Time Resolution']='obs_freq',
    trans['Observation Times']='obs_time_GMT',
    trans['Units']='measurementUnit_name','measurementUnit_abbreviation'  # 'lot2units(vname,units)',
    trans['Variable Instrument']='stationPlatformModel',
    trans['Type of Access (Public or Restricted)']='station_data_policy','useConstraints'
    trans['Data Owner']='pointOfContact_Organisation',
    trans['Data Owner Link or contact']='pointOfContact_MailAddress',
    trans['Data Provider']='principalInvestigator_Name',
    trans['Data Provider (Link or contact)']='principalInvestigator_OnlineResource',
    trans['State of Data Rescue']='statusOfObservation',
    trans['Project Status']='statusOfObservation',

    return trans

def odb_trans(odbkeys,unified):
    trans=dict(zip(odbkeys,['']*len(odbkeys)))
    trans['statid@hdr']='station_id',
    #trans['expver@desc']=('source_uid','data_repository_ftp')
    trans['expver']=('source_uid','data_repository_ftp')
    trans['source@hdr']='source_name',
    trans['collection_identifier@conv']='originator_Organisation',
    trans['obstype@hdr']='stationPlatformType',
    trans['codetype@hdr']='stationPlatformModel',
    trans['stalt@hdr']='elev',
    #trans['Original_Latitude_Units']='stationDescription',
    #trans['Original_Longitude_Units']='stationDescription',
    #trans['Original_Altitude_Units']='stationDescription',
    trans['lat@hdr']='lat',
    trans['lon@hdr']='lon',
    #trans['Last modified']='creationDate',
    #trans['Start Station Year Month Day']='startOfOperation',
    #trans['Start Record Year Month Day']='lot2start(vname)',
    #trans['End Record Year Month Day']='lot2end(vname)',
    #trans['WMO Region']='continent',
    #trans['Current Station Name']='station name',
    #trans['WMO_ID']='wmo_id','stationPlatformUniqueIdentifier', # WIGOS ID
    #trans['Network1_name']='network_Names',
    #trans['Network1_ID']='network_IDs',
    #trans['Network2_name']='network_Names',
    #trans['Network2_ID']='network_IDs',
    #trans['Country or Region']='country',
    #trans['Time Resolution']='obs_freq',
    #trans['Observation Times']='obs_time_GMT',
    #trans['Units']='measurementUnit_name','measurementUnit_abbreviation'  # 'lot2units(vname,units)',
    #trans['Variable Instrument']='stationPlatformModel',
    #trans['Type of Access (Public or Restricted)']='station_data_policy','useConstraints'
    #trans['Data Owner']='pointOfContact_Organisation',
    #trans['Data Owner Link or contact']='pointOfContact_MailAddress',
    #trans['Data Provider']='principalInvestigator_Name',
    #trans['Data Provider (Link or contact)']='principalInvestigator_OnlineResource',
    #trans['State of Data Rescue']='statusOfObservation',
    #trans['Project Status']='statusOfObservation',

    kdict=dict(zip(list(unified.columns),list(range(len(unified.columns)))))
    kdict['1']=kdict['Geopotential_height_start_year']
    kdict['2']=kdict['Temperature_start_year']
    kdict['29']=kdict['Relative_humidity_start_year']
    kdict['59']=kdict['Dew_Point_Temperature_start_year']
    kdict['111']=kdict['Wind_speed_start_year']
    kdict['112']=kdict['Wind_direction_start_year']

    return trans,kdict

vndict={'wind_speed':'Wind_speed','wind_from_direction':'Wind_direction',
        'dew_point_temperature':'Dew_Point_Temperature','air_temperature':'Temperature',
        'air_pressure':'Air_pressure','pressure_temp':'Air_pressure',
        'relative_humidity':'Relative_humidity','geopotential_height':'Geopotential_height'}

def uai_unified(trans,unified,uai,fips=None):

    kdict={}
    for k,v in trans.items():
        for vv in v:
            if '(' in vv:
                pass
            else:
                kdict[vv]=unified.columns.get_loc(vv)
    #kdict['Units']=unified.columns.get_loc('measurementUnit_name')
    #kdict['Units']=unified.columns.get_loc('measurementUnit_abbreviation')
    kdict['Original Units']=unified.columns.get_loc('measurementOriginalUnit_name')
    kdict['fips_code']=unified.columns.get_loc('fips_code')
    for k,v in vndict.items():
        vv=v+'_start_year'
        kdict[vv]=unified.columns.get_loc(vv)
        vv=v+'_end_year'
        kdict[vv]=unified.columns.get_loc(vv)

    vlist=list(trans.values())
    l=-1
    oldstationname='----'
    oldsource='----'
    ulist=[]
    row=['']*len(unified.columns)
    for i in range(len(uai)):
        if uai['Station Name'][i]!=oldstationname or uai['Source'][i]!=oldsource:
            oldstationname=uai['Station Name'][i]
            oldsource=uai['Source'][i]
            l+=1
            #unified.loc[l]=['']*100
            ulist.append(row[:])

            for k,v in trans.items():

                for vv in v:
                    if '(' in vv:
                        vname=uai['Variable Name'][i]
                        #units=uai['Units'][i]
                        try:
                            s=eval(vv)
                            ulist[-1][kdict[s]]=uai[k][i]
                        except:
                            print('No variable name for '+vname)
                    else:
                        if ulist[-1][kdict[vv]]=='':
                            ulist[-1][kdict[vv]]=uai[k][i]
                        else:
                            if vlist.count(v)>1:
                                if uai[k][i]!='' and uai[k][i]!=ulist[-1][kdict[vv]]:
                                    ulist[-1][kdict[vv]]+=' , '+uai[k][i]
                    if vv=='measurementUnit_name':
                        print(uai[k][i])
                    if vv=='country':
                        try:
                            cc=fips.at[numpy.where(fips['Name']==uai[k][i])[0][0],'FIPS 10-4']
                            ulist[-1][kdict['fips_code']]=cc
                        except:
                            print('country code for',uai[k][i],'not found')
                    if vv=='stationPlatformUniqueIdentifier':
                        if int(uai[k][i])>1000 and int(uai[k][i])<100000:
                            ulist[-1][kdict[vv]]='0-20000-0-'+uai[k][i]



            #print(unified.loc[l])
        else:
            for k,v in trans.items():
                for vv in v:
                    if vv=='measurementUnit_name':
                        print(uai[k][i])

                    if '(' in vv:
                        vname=uai['Variable Name'][i]
                        #units=uai['Units'][i]
                        try:
                            s=eval(vv)
                            if uai[k][i]!=ulist[-1][kdict[s]]:                                                   
                                if ulist[-1][kdict[s]]=='': 
                                    ulist[-1][kdict[s]]=uai[k][i]
                                else:
                                    ulist[-1][kdict[s]]+=' , '+uai[k][i]
                            #print(l,uai['Station Name'][i],s,ulist[-1][kdict[s]])
                        except:
                            print('No variable name for '+vname)
                    #else:
                        #ulist[-1][kdict[vv]]+=uai[k][i]
        print(i)
    unified=unified.append(pd.DataFrame(ulist,columns=unified.columns))           
    return unified

def lot2start(vname):
    return vndict[vname]+'_start_year'

def lot2end(vname):
    return vndict[vname]+'_end_year'

def lot2units(vname,units):
    return vndict[vname]

def do_row(row,kdict,meta,wmoid,nw,lat,lon,name) :
    try:
        if nw=='C' and 'Group' not in meta.columns:
            return []

        meta.name=name   
        if meta.name=='OSCAR':
            z=numpy.where(meta['StationId']=='0-20000-0-{:0>5}'.format(int(wmoid)))[0]
            if len(z)==0:
                dists=fdist(lat,lon,meta['Latitude'].values,meta['Longitude'].values)*180./math.pi
                z=numpy.nanargmin(dists)
                if dists[z]<0.5:
                    z=[z]
                    print(meta['StationId'][z[0]],meta['StationName'][z[0]])
                else:
                    z=[]
        elif meta.name=='IGRA2':
            z=numpy.where(meta['StationId']==wmoid)[0]
            if len(z)==0:
                dists=fdist(lat,lon,meta['Latitude'].values,meta['Longitude'].values)*180./math.pi
                z=numpy.nanargmin(dists)
                if dists[z]<0.5:
                    z=[z]
                    #print(meta['StationId'][z[0]],meta['StationName'][z[0]])
                else:
                    z=[]
        elif meta.name=='WBAN':
            z=numpy.where(meta['StationId']==wmoid)[0]
            if len(z)==0:
                dists=fdist(lat,lon,meta['Latitude'].values,meta['Longitude'].values)*180./math.pi
                z=numpy.nanargmin(dists)
                if dists[z]<0.5:
                    z=[z]
                    #print(meta['StationId'][z[0]],meta['StationName'][z[0]])
                else:
                    z=[]
                                        
        else:
            z=numpy.where(meta['StationId']==wmoid)[0]


        if len(z)==0:
            return []
        else:
            z=z[0]
            metaz=meta.loc[z]
            if fdist(metaz['Latitude'],metaz['Longitude'],lat,lon)*180./math.pi<1.0: 
                row[kdict['station name']]=metaz['StationName']
                row[kdict['stationDescription']]=metaz['StationName']
                #print(wmoid,': ',meta.name+' match')
                row[kdict['stationPlatformUniqueIdentifier']]=metaz['StationId']
                row[kdict['station name']]=metaz['StationName']
                row[kdict['stationDescription']]=metaz['StationName']
                if meta.name=='OSCAR': 
                    row[kdict['regionOfOriginOfData']]=metaz['RegionName']
                    row[kdict['territoryOfOriginOfData_Name']]=metaz['CountryArea']
                    row[kdict['territoryOfOriginOfData_ISO3CountryCode']]=metaz['CountryCode']

                if meta.name=='CHUAN': 
                    row[kdict['descriptionDataset']]=metaz['Source']
            else:
                #print('Position Mismatch:',name,wmoid,metaz['Latitude'],metaz['Longitude'],lat,lon)
                return []


    except KeyError as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        #print(exc_type, fname, exc_tb.tb_lineno)

    return [z]

def insert(row,vola,chuan,igrainv,wbaninv,transodb,trans,kdict,unified,rdata,varno,dmin,dmax,dcount,fi,fn):

    #try:

    ids=-1
    wmoid=''
    iwmoid=0
    if not fi:
        uai=dict(zip(rdata[0].split(),zip(rdata[1].split(),rdata[-2].split())))
        uai['statid@hdr']=list(uai['statid@hdr'])
        if uai['statid@hdr'][0]=="''":
            for i in range(len(uai['statid@hdr'])):
                uai['statid@hdr'][i]="'C:"+fn.split(':')[-1]
        vlist=list(transodb.values())

        #print('HERE', fn )
        for k,v in transodb.items():
            for vv in v:
                if row[kdict[vv]]=='':
                    try:
                        row[kdict[vv]]=uai[k][0]
                    except:
                        #print('skipping' , vv , ' ' , k )
                        row[kdict[vv]] = ''
                    #row[kdict[vv]]=uai[k][0]  
                else:
                    if vlist.count(v)>1:
                        if uai[k][0]!=row[kdict[vv]]:
                            row[kdict[vv]]+=' , '+uai[k][0]

            if k=='statid@hdr':   # extract information from vola data base
                if ':' in uai['statid@hdr'][0]:
                    try:
                        nw=uai['statid@hdr'][0][1:-1].split(':')[0]
                        if nw=='C':
                            wmoid=fn.split(':')[-1]
                            if wmoid[-1] in 'AB':
                                wmoid=wmoid[:-1]
                        else:
                            wmoid=uai['statid@hdr'][0][1:-1].split(':')[1]

                        lat=float(uai['lat@hdr'][0])
                        lon=float(uai['lon@hdr'][0])
                        iwmoid=int(wmoid)
                    except:
                        #print('error parsing wmo ID',wmoid,uai['lat@hdr'][0],uai['lon@hdr'][0])
                        break                    
                else:
                    nw=''
                    try:
                        lat=float(uai['lat@hdr'][0])
                        lon=float(uai['lon@hdr'][0])
                    except:
                        lat=numpy.nan
                        lon=numpy.nan
                    wmoid=uai['statid@hdr'][0]
                    if "'" in wmoid:
                        wmoid=wmoid[1:-1]
                    if len(wmoid)>3:
                        try:                                   # without the try, the station id e.g. W5331 fails. Is this fix correct?
                              iwmoid=int(wmoid)
                        except:
                              iwmoid = 0
                              
                #print('HERE 2 ', fn )
                if iwmoid>1000 and iwmoid<100000:
                    mi=do_row(row,kdict,vola,wmoid,nw,lat,lon,'OSCAR')
                    if not mi:
                        mi=do_row(row,kdict,igrainv,wmoid,nw,lat,lon,'IGRA2')
                        if not mi:
                            mi=do_row(row,kdict,wbaninv,wmoid,nw,lat,lon,'WBAN')
                            if not mi:
                                mi=do_row(row,kdict,chuan,wmoid,nw,lat,lon,'CHUAN')
                                if not mi:
                                        
                                    #print(wmoid,'not found, not identifiable!')
                                    wait_time=0.2
                                    #print('HERE 3 ', fn )
                                    #while not os.access(fn.split('/')[0]+'/'+'orphans.csv', os.W_OK ):
                                    #    if os.path.isfile(home + fn.split('/')[0]+'/'+'orphans.csv'):
                                    #        time.sleep(wait_time)
                                    orphans = home + '/' + fn.split('/')[0] + '_orphans.csv'
                                    orph = open(orphans, 'a')
                                    while not os.access(orphans, os.W_OK ):                                                                                                                                                                                                                  
                                        if os.path.isfile(home + fn.split('/')[0]+'/'+'orphans.csv'):                                                                                                                                                                                                                   
                                            time.sleep(wait_time)  
                                        else:
                                           orphans = home + '/' + fn.split('/')[0] + '_orphans.csv'   
                                           
                                    #print('HERE 4 ', fn )
                                    with open(orphans, 'a') as f:
                                        x=numpy.where(wbaninv.StationId==wmoid)
                                        if len(x[0])==0:
                                            f.write(fn+','+wmoid+',{:7.2f},{:7.2f},,,,\n'.format(lat,lon))
                                        else:
                                            f.write(fn+','+wmoid+',{:7.2f},{:7.2f},{:7.2f},{:7.2f},{},{}\n'.format(lat,lon,wbaninv.Latitude.values[x[0][0]],
                                                                                   wbaninv.Longitude.values[x[0][0]],
                                                                                   numpy.int(wbaninv.From.values[x[0][0]]),numpy.int(wbaninv.To.values[x[0][0]])))
                                            
                                    ids=-1
                                else:
                                    ids=4
                            else:
                                ids=2
                        else:
                            if igrainv['Network'].values[mi]=='M':
                                #print('IGRA2 WMO station')
                                ids=0
                            else:
                                ids=3
                    else:
                        ids=0
                else:
                    ids=-1
                row[kdict[varno]]=dmin
                row[kdict[varno]+1]=dmax


    #except Exception as e:
        #exc_type, exc_obj, exc_tb = sys.exc_info()
        #fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        #print(exc_type, fname, exc_tb.tb_lineno)

    cc=row[kdict['territoryOfOriginOfData_ISO3CountryCode']]
    return ids,wmoid,row[kdict['territoryOfOriginOfData_ISO3CountryCode']]

def odb_process(ci,unified,vola,chuan,igrainv,wbaninv,trans,fn):           
    #print('*** Processing ', fn, ' ' , now(time.time()) )
    transunified=[]
    fi=False
    try:
        for varno in ['111','112','1','2','59','29']:
            qs='select min(date),max(date) where varno='+varno #distinct,count(date)
            try:
                rdata=subprocess.check_output(["odb","sql","-q",qs,"-i",fn,'--no_alignment'],stderr=subprocess.DEVNULL)
                dmin,dmax=rdata.decode().split('\n')[1].split()
                dcount=0
                if dmin!='NULL':
                    qs='select * where varno='+varno+' and (date='+dmin+' or date='+dmax+')'
                    rdata=subprocess.check_output(["odb","sql","-q",qs,"-i",fn,'--no_alignment'],stderr=subprocess.DEVNULL).decode().split('\n')
                    if not ci:
                        ci=True
                        cols=rdata[0].split()
                        transodb,kdict=odb_trans(cols,unified)
                    if not fi:
                        transunified.append(row[:])
                    insert(transunified[-1],vola,chuan,igrainv,wbaninv,transodb,trans,kdict,unified,rdata,varno,dmin,dmax,dcount,fi,fn)
                    fi=True


            except subprocess.CalledProcessError as e:
                pass
                #print(e)
    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        #print(exc_type, fname, exc_tb.tb_lineno)
        exit()

    if len(transunified)==0:
        transunified.append(row[:])

    return transunified[0]

def scini(station_configuration,ids,wmoid,cwmoid,vola,transunified,kdict,dmax,dmin,cdm,cities,cc, fn ):

    try:
        station_configuration['longitude']=float(transunified[-1][kdict['lon']])
        station_configuration['latitude']=float(transunified[-1][kdict['lat']])
    except:
        station_configuration['longitude']=numpy.nan
        station_configuration['latitude']=numpy.nan
    station_configuration['station_name']=transunified[-1][kdict['station name']]
    if ids==0:
        station_configuration['primary_id']='0-20000-0-'+wmoid
        station_configuration['primary_id_scheme']=5
        #if wmoid!=cwmoid:
            
        station_configuration['secondary_id']=cwmoid
        station_configuration['secondary_id_scheme']=ids
    else:
        station_configuration['secondary_id']=cwmoid
        station_configuration['secondary_id_scheme']=ids
        dists=fdist(station_configuration['latitude'],station_configuration['longitude'],
                    vola['Latitude'].values,vola['Longitude'].values)*180./math.pi
        
        idy=numpy.argsort(dists)
        #idx=numpy.argmin(dists)
        i=0
        while dists[idy[i]]<0.5:
            idx=idy[i]
            if '0-724-0-0' in  vola['StationId'][idx]:
                i+=1
                continue
            station_configuration['primary_id']=vola['StationId'][idx]
            station_configuration['primary_id_scheme']=5
            station_configuration['station_name']=vola.StationName[idx]
            cc= vola['CountryCode'][idx]
            #print(wmoid,'colocated with WMO station',station_configuration['primary_id'])
            if 'Upper' in vola['ObsRems'].values[idx] or '0-2000' in vola['StationId'][idx]:
                #print(wmoid,'colocated with WMO station',station_configuration['primary_id']+', used this')
                break
            i+=1

    station_configuration['station_crs']=0
    station_configuration['local_gravity']=numpy.NaN
    try:
        if '-' in dmin:
            station_configuration['start_date']=datetime.datetime.strptime(dmin,'%Y-%m-%d')
            station_configuration['end_date']=datetime.datetime.strptime(dmax,'%Y-%m-%d')
        else:
            station_configuration['start_date']=datetime.datetime.strptime(dmin,'%Y%m%d')
            station_configuration['end_date']=datetime.datetime.strptime(dmax,'%Y%m%d')
    except:
        #print('invalid date encountered:',dmin,dmax)
        station_configuration['start_date']='NA'
        station_configuration['end_date']='NA' 
        
    station_configuration['station_type']=1 # Land Station
    station_configuration['platform_type']=0 # Land Station
    station_configuration['platform_sub_type']=63 
    station_configuration['operating_institute']='NMS' 
    try:
        station_configuration['operating_territory']=numpy.where(cdm['sub_region']['alpha_3_code']==cc)[0][0]
    except:
        pass
    try:
        if ~numpy.isnan(station_configuration['latitude']):
            dist=fdist(station_configuration['latitude'],station_configuration['longitude'],cities['Latitude'].values,cities['Longitude'].values)
            station_configuration['city']=cities.AccentCity.iloc[numpy.argmin(dist)]
        else:
            station_configuration['city']='NA' 
    except:
        station_configuration['city']='NA' 
    station_configuration['contact']='NA' 
    station_configuration['role']='NA' 
    try:
        lx=set(r.split('\t')[5] for r in rdata[:-1])                       
        station_configuration['observing_frequency']=len(lx)-1
        station_configuration['reporting_time']=lx[1:]
    except:
        station_configuration['observing_frequency']=0
        station_configuration['reporting_time']=-1
    station_configuration['reporting_time']=63 
    station_configuration['telecommunication_method']=29
    station_configuration['station_automation']=1 
    station_configuration['measuring_system_model']='NA' 
    station_configuration['measuring_system_id']='NA'
    station_configuration['observed_variables']=[] 
    station_configuration['comment']='NA' 
    station_configuration['bbox_min_longitude']=station_configuration['longitude']-0.5/numpy.cos(station_configuration['latitude']*numpy.pi/180.)
    station_configuration['bbox_min_latitude']=station_configuration['latitude']-0.5
    station_configuration['bbox_max_longitude']=station_configuration['longitude']+0.5/numpy.cos(station_configuration['latitude']*numpy.pi/180.) 
    station_configuration['bbox_max_latitude']=station_configuration['latitude']+0.5
    if station_configuration['bbox_min_longitude']<-180.:
        station_configuration['bbox_min_longitude']+=360.
    if station_configuration['bbox_max_longitude']>180.:
        station_configuration['bbox_max_longitude']-=360.
    if station_configuration['bbox_min_latitude']<-90.:
        station_configuration['bbox_min_latitude']=-90.
    if station_configuration['bbox_max_latitude']>90.:
        station_configuration['bbox_max_latitude']=90.

    station_configuration['metadata_contact']='S. Schroeder'
    station_configuration['metadata_contact_role']=0

    #code.interact(local = locals())
    #print(station_configuration['primary_id'])
    #print(station_configuration['secondary_id'])
    
    if fn != 0:
          
          db = fn.split("/")[0]
          out = open(home + '/' + db + '_summary.txt', 'a')
          try:
                Line = station_configuration['primary_id'] + '\t' + str(station_configuration['secondary_id']) + '\t' + station_configuration['station_name'] + '\t'
          except:
                Line = ' ' + '\t' + str(station_configuration['secondary_id']) + '\t' + str(station_configuration['station_name']) + '\t'
          Line = Line + str(station_configuration['latitude']) + '\t' + str(station_configuration['longitude']) + '\t' + fn + '\n'
          out.write(Line)
            
    return    

def odb_cdm(ci,cdm,cdmd,unified,vola,chuan,igrainv,wbaninv,trans,fn):           
    print('*** Processing the file odb_cdm ::: ', fn, ' ' , now(time.time()) )
    transunified=[]
    station_configuration=dict()
    for k in cdmd['station_configuration'].element_name.values:
        station_configuration[k]=list()
    fi=False
    #try:
    for varno in ['111','112','1','2','59','29']:
        qs='select min(date),max(date) where varno='+varno #distinct,count(date)
        try:
            #print('Checking the data ***', fn )  
            rdata=subprocess.check_output(["odb","sql","-q",qs,"-i",fn,'--no_alignment'],stderr=subprocess.DEVNULL)
            dmin,dmax=rdata.decode().split('\n')[1].split()
            #print ('Checked the data ***', fn , ' ' , dmin , ' ' , dmax )
            dcount=0
            if dmin!='NULL':
                if int(dmin)<19000000 or int(dmax)>21000000:
                    #print(fn,': MALFORMED DATE ERROR!!')
                    qs='select min(date),max(date) where varno='+varno +' and date>19000000' #distinct,count(date)
                    rdata=subprocess.check_output(["odb","sql","-q",qs,"-i",fn,'--no_alignment'],stderr=subprocess.DEVNULL)
                    dmin,dmax=rdata.decode().split('\n')[1].split()
                    #print('*** checked the rdata subprocess outside boundaries')

                qs='select * where varno='+varno+' and (date='+dmin+' or date='+dmax+')'
                rdata=subprocess.check_output(["odb","sql","-q",qs,"-i",fn,'--no_alignment'],stderr=subprocess.DEVNULL).decode().split('\n')
                #print('*** checked againg the rdata ', fn )
                if not ci:
                    ci=True
                    cols=rdata[0].split()
                    if 'C:' in fn:
                        idx=cols.index('lat@hdr')
                        if rdata[1].split('\t')[idx]=='NULL':
                            try:
                                idc=numpy.where(chuan['StationId']==fn.split(':')[-1])[0][0]
                                slat='.'.join(chuan['Latitude'][idc].split(','))
                                slon='.'.join(chuan['Longitude'][idc].split(','))
                                for i in range(len(rdata[1:])):
                                    if rdata[i+1]=='':
                                        break
                                    rl=rdata[i+1].split('\t')
                                    rl[idx]=slat[:]
                                    rl[idx+1]=slon[:]
                                    rdata[i+1]='\t'.join(rl)
                                #print(fn,' lat lon supplied from metadata on ',i,'lines')
                            except:
                                #print(fn,' lat lon could not be supplied from metadata')
                                pass
                    transodb,kdict=odb_trans(cols,unified)
                    #print('*** Done with the odb_trans ', fn )
                if not fi:
                    transunified.append(row[:])
                    #print('*** Doing the insert ', fn )
                    ids,wmoid,cc=insert(transunified[-1],vola,chuan,igrainv,wbaninv,transodb,trans,kdict,unified,rdata,varno,dmin,dmax,dcount,fi,fn)
                    #print('*** Done with the insert', fn )
                    if ':' in fn:  
                        fnl=fn.split(':')
                        cwmoid=fnl[-2][-1]+':'+fnl[-1]
                    else:
                        cwmoid=wmoid
                    #print ('*** Doing the scini ', fn )
                    scini(station_configuration,ids,wmoid,cwmoid,vola,transunified,kdict,dmax,dmin,cdm,cities,cc,fn) 
                    #if '.3188.' in fn:
                        #idx=numpy.where(chuan['WMO ID']==wmoid)[0]
                        #if len(idx)>0:
                            #station_configuration['secondary_id']=chuan['StationId'].values[idx]
                            #station_configuration['secondary_id_scheme']=4
                        #else:
                            #idy=numpy.where(chuan['StationId']==wmoid)[0]
                            #if len(idy)==0:
                                #print('CHUANstation not identified')
                            #else:
                                #station_configuration['secondary_id']=wmoid
                                #station_configuration['secondary_id_scheme']=4
                                
                    
                    fi=True
                station_configuration['observed_variables'].append(int(varno)) 


        except subprocess.CalledProcessError as e:      
             print(e)

    #except Exception as e:
        #exc_type, exc_obj, exc_tb = sys.exc_info()
        #fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        #print(exc_type, fname, exc_tb.tb_lineno)
        #exit()

    if len(transunified)==0:
        transunified.append(row[:])

    if station_configuration['observed_variables']==[] or numpy.isnan(station_configuration['latitude']):
        #print(fn,'insufficient data for assimilation')
        return None
    
    if 2 in station_configuration['observed_variables']:
        station_configuration['station_configuration_code']=13 # Radiosonde, taken from station_configuration_code cdm table
    else:
        station_configuration['station_configuration_code']=8 # Pilot, taken from station_configuration_code cdm table

    #print('*** DONE odb_cdm ', fn )
    return station_configuration,transunified[0] 

def read_bufr_stn_meta(varno,bufrfile):

    tx=time.time()
    vdict={'111':'windDirection','112':'windSpeed','1':'geopotentialHeight',
           '2':'airTemperature','59':'dewpointTemperature','29':'relativeHumidity'}
    vlog={'111':False,'112':False,'1':False,'2':False,'59':False,'29':False}
    transunified=[]
    fi=False

    with open(bufrfile) as f:
        cnt = 0
        # loop over the messages in the file
        while 1:
            # get handle for message
            bufr = codes_bufr_new_from_file(f)
            if bufr is None:
                break
            codes_release(bufr)
            cnt+=1

    if cnt==0:
        return None
    rdata=['expver@desc\tobstype@hdr\tcodetype@hdr\tstalt@hdr\tstatid@hdr\tdate@hdr\ttime@hdr\tlat@hdr\tlon@hdr\tsource@hdr\tcollection_identifier@conv\tobsvalue@body\tvarno@body\tvertco_type@body\tvertco_reference_1@body']
    with open(bufrfile) as f:
        # loop over the messages in the file
        for i in range(cnt):
            # get handle for message
            bufr = codes_bufr_new_from_file(f)
            if i not in [0,cnt-1]:
                codes_release(bufr)
                continue
            else:
                codes_set(bufr, 'unpack', 0)
                # get all the timePeriods
                iterid = codes_bufr_keys_iterator_new(bufr)

                # loop over the keys
                if codes_get_array(bufr,'dataSubCategory')[0] not in [91,101]:
                    #print codes_get_array(bufr,'dataSubCategory')[0]
                    codes_release(bufr)
                    continue
                d=dict()
                while codes_bufr_keys_iterator_next(iterid):

                    # print key name
                    keyname = codes_bufr_keys_iterator_get_name(iterid)
                    if keyname in ['localYear','localMonth','localDay','localHour']:
                        d[keyname]=codes_get_array(bufr,keyname)
                    #print (keyname,codes_get_array(bufr,keyname))
                    if keyname=='localHour':
                        datum='{}-{:0>2}-{:0>2}'.format(d['localYear'][0],d['localMonth'][0],d['localDay'][0])
                        hour='{:0>2}'.format(d['localHour'][0])
                        #print(datum)
                # delete the key iterator
                codes_bufr_keys_iterator_delete(iterid)

                if i==0:
                    lat = codes_get(bufr, "latitude")
                    lon = codes_get(bufr, "longitude")
                    alt = float(codes_get(bufr, "heightOfStation"))
                    blockNumber = codes_get(bufr, "blockNumber")
                    stationNumber = codes_get(bufr, "stationNumber")
                    statid = str(blockNumber*1000+stationNumber)
                    print(statid, "                      FFFFFFFFFFFFFF")
                    dmin=datum
                    rdata=['expver@desc\tobstype@hdr\tcodetype@hdr\tstalt@hdr\tstatid@hdr\tdate@hdr\ttime@hdr\tlat@hdr\tlon@hdr\tsource@hdr\tcollection_identifier@conv\tobsvalue@body\tvarno@body\tvertco_type@body\tvertco_reference_1@body']
                    rdata.append('\t'.join(['ERA40','{}'.format(codes_get_array(bufr,'dataSubCategory')[0]),'35',str(alt),
                                            '{:0>2}{:0>3}'.format(blockNumber,stationNumber),'x','x',
                                            '{:7.3f}'.format(lat),'{:8.3f}'.format(lon),'ECMWF','x','x','x','x','x']))
                    rdata.append('\t'.join(['ERA40','{}'.format(codes_get_array(bufr,'dataSubCategory')[0]),'35',str(alt),
                                            '{:0>2}{:0>3}'.format(blockNumber,stationNumber),'x','x',
                                            '{:7.3f}'.format(lat),'{:8.3f}'.format(lon),'ECMWF','x','x','x','x','x']))

                dmax=datum
                for varno in ['111','112','1','2','59','29']:
                    try:
                        airTemperature = codes_get_array(bufr, vdict[varno]) 
                        vlog[varno]=True
                    except:
                        pass

                codes_release(bufr)

        ci=False       
        station_configuration={}
        for k in cdmd['station_configuration'].element_name.values:
            station_configuration[k]=-1
        for varno in vlog.keys():

            if not ci:
                ci=True
                cols=['statid@hdr','expver@desc','source@hdr','collection_identifier@conv','obstype@hdr',
                      'codetype@hdr','stalt@hdr','lat@hdr','lon@hdr']
                transodb,kdict=odb_trans(cols,unified)
                dcount=0
            if not fi:
                transunified.append(row[:])
            ids,wmoid,cc=insert(transunified[-1],vola,chuan,igrainv,wbaninv,transodb,trans,kdict,unified,rdata,varno,dmin,dmax,dcount,fi,bufrfile)
            if not fi:
                # def scini(station_configuration,ids,wmoid,cwmoid,vola,transunified,kdict,dmax,dmin,cdm,cities,cc, fn):
                cwmoid = statid
                scini(station_configuration,ids,wmoid,cwmoid,vola,transunified,kdict,dmax,dmin,cdm,cities,cc,bufrfile)
                fi=True
            station_configuration['observed_variables'].append(int(varno)) 



    #print ('/'.join(bufrfile.split('/')[-1:]),cnt,"messages",time.time()-tx)

    if len(transunified)==0:
        transunified.append(row[:])

    return station_configuration,transunified[0]

def get_secondary_ncar(ncfile):
      """ This module extracts the list of secondary station ids inside the original ncar files.
          First it gtes the original file name, searches for it in the file directory given,
          then opens the original files and loops over the headers to get all the stations ids. """
      # UADB_trhc_074693.nc
      # UADB_windc_078730.nc
      # uadb_trhc_62378.txt
      # uadb_windc_70308.txt
      name = ncfile.split("/")[1]
      name = name.replace("UADB_","uadb_").replace("c_00","c_")
      name = name.replace("c_0","c_").replace('nc', 'txt')
      data_source = "/raid60/scratch/federico/databases/UADB"
      
      file_path = data_source + '/' + name 
      
      if os.path.isfile(file_path):
            print('*** Original ncar file found ')
            print('Hello')
            data = open(file_path, 'r').read().splitlines() 
            all_ids = []
            for i, line in enumerate(data):
                  if line[0] == 'H':
                      ident    = line[15:21].replace(' ','')# WMO
                      if ident not in all_ids:
                            all_ids.append(str(ident))
                            
            print("All_id: " , all_ids )
            if len(all_ids) >1:
                  All = ','.join(all_ids)
            else:
                  All = all_ids[0]
      else:
            print("*** Could not find the original file!")
            All = ''
      return All      

def read_rda_meta(ncfile):

    secondary = get_secondary_ncar(ncfile)
    
    tx=time.time()
    vdict={'windd':'111','winds':'112','gph':'1',
           'temp':'2','rhumi':'29'}
    vlog={'111':False,'112':False,'1':False,'2':False,'59':False,'29':False}
    transunified=[]
    fi=False
    #try:
    rdata=['expver@desc\tobstype@hdr\tcodetype@hdr\tstalt@hdr\tstatid@hdr\tdate@hdr\ttime@hdr\tlat@hdr\tlon@hdr\tsource@hdr\tcollection_identifier@conv\tobsvalue@body\tvarno@body\tvertco_type@body\tvertco_reference_1@body']
    nchdr='UADB_Station'.join(ncfile.split('UADB'))
    with Dataset(nchdr) as f:
        ds=xr.open_dataset(ncfile)
        d=dict()
        lat = f.variables['lat'][0]
        lon = f.variables['lon'][0]
        alt = f.variables['alt'][0]
        dmin=str(ds.date[0].values.astype('datetime64[D]'))
        dmax=str(ds.date[-1].values.astype('datetime64[D]'))
        rdata=['expver@desc\tobstype@hdr\tcodetype@hdr\tstalt@hdr\tstatid@hdr\tdate@hdr\ttime@hdr\tlat@hdr\tlon@hdr\tsource@hdr\tcollection_identifier@conv\tobsvalue@body\tvarno@body\tvertco_type@body\tvertco_reference_1@body']
        rdata.append('\t'.join(['RDA','101','35',str(alt),
                                nchdr[-8:-3],'x','x',
                                '{:7.3f}'.format(lat),'{:8.3f}'.format(lon),'NCAR','x','x','x','x','x']))

        for varno in ['windd','winds','gph','temp','rhumi']:
            try:
                if not numpy.isnan(numpy.nanmean(ds[varno])):
                    vlog[vdict[varno]]=True
            except:
                pass

        ci=False       
        station_conf={}
        for k in cdmd['station_configuration'].element_name.values:
            station_conf[k]=-1
        for varno in vlog.keys():

            if not ci:
                ci=True
                cols=['statid@hdr','expver@desc','source@hdr','collection_identifier@conv','obstype@hdr',
                      'codetype@hdr','stalt@hdr','lat@hdr','lon@hdr']
                transodb,kdict=odb_trans(cols,unified)
                dcount=0
            if not fi:
                transunified.append(row[:])
            ids,wmoid,cc=insert(transunified[-1],vola,chuan,igrainv,wbaninv,transodb,trans,
                                kdict,unified,rdata,varno,dmin,dmax,dcount,fi,ncfile)
            if not fi:
                cwmoid = secondary # list of extracted secondary ids 
                scini(station_conf,ids,wmoid,cwmoid, vola,transunified,kdict,dmax,dmin,cdm,cities,cc,ncfile)
                fi=True
            station_conf['observed_variables'].append(int(varno)) 



    #print ('/'.join(ncfile.split('/')[-1:]),time.time()-tx)


    if len(transunified)==0:
        transunified.append(row[:])

    return station_conf, transunified[0]

def read_igra_meta(line):

    tx=time.time()
    vdict={'windd':'111','winds':'112','gph':'1',
           'temp':'2','rhumi':'29'}
    vlog={'111':False,'112':False,'1':False,'2':False,'59':False,'29':False}
    transunified=[]
    fi=False
    #try:
    rdata=['expver@desc\tobstype@hdr\tcodetype@hdr\tstalt@hdr\tstatid@hdr\tdate@hdr\ttime@hdr\tlat@hdr\tlon@hdr\tsource@hdr\tcollection_identifier@conv\tobsvalue@body\tvarno@body\tvertco_type@body\tvertco_reference_1@body']
    d=dict()
    lat = line['Latitude']
    lon = line['Longitude']
    alt = line['Elev']
    dmin=str(line['From'])+'-01-01'
    dmax=str(line['To'])+'-12-31'
    rdata=['expver@desc\tobstype@hdr\tcodetype@hdr\tstalt@hdr\tstatid@hdr\tdate@hdr\ttime@hdr\tlat@hdr\tlon@hdr\tsource@hdr\tcollection_identifier@conv\tobsvalue@body\tvarno@body\tvertco_type@body\tvertco_reference_1@body']
    rdata.append('\t'.join(['IGRA','101','35',str(alt),
                            line['StationId'],'x','x',
                            '{:7.3f}'.format(lat),'{:8.3f}'.format(lon),'IGRA','x','x','x','x','x']))

    for varno in ['windd','winds','gph','temp','rhumi']:
        try:
            #if not numpy.isnan(numpy.nanmean(ds[varno])):
            vlog[vdict[varno]]=True
        except:
            pass

    ci=False       
    station_conf={}
    for k in cdmd['station_configuration'].element_name.values:
        station_conf[k]=-1
    for varno in vlog.keys():

        if not ci:
            ci=True
            cols=['statid@hdr','expver@desc','source@hdr','collection_identifier@conv','obstype@hdr',
                  'codetype@hdr','stalt@hdr','lat@hdr','lon@hdr']
            transodb,kdict=odb_trans(cols,unified)
            dcount=0
        if not fi:
            transunified.append(row[:])
        ids,wmoid,cc=insert(transunified[-1],vola,chuan,igrainv,wbaninv,transodb,trans,
                            kdict,unified,rdata,varno,dmin,dmax,dcount,fi,line['StationId'])
        if not fi:
            scini(station_conf,ids,wmoid,wmoid, vola,transunified,kdict,dmax,dmin,cdm,cities, cc, 0)
            fi=True
        station_conf['observed_variables'].append(int(varno)) 



    #print ('/'.join(line['StationId'].split('/')[-1:]),time.time()-tx)

    if len(transunified)==0:
        transunified.append(row[:])

    return station_conf,transunified[0]


if __name__ == '__main__':


    parser = argparse.ArgumentParser(description="Analyze the odb inventory")
    parser.add_argument('--database_dir' , '-d', 
                    help="Optional: path to the database directory. If not given, will use the files in the data directory" ,
                    default = '../data',
                    type = str)
    parser.add_argument('--auxtables_dir' , '-a', 
                    help="Optional: path to the auxiliary tables directory. If not given, will use the files in the data/tables directory" ,
                    default = '../data/tables/',
                    type = str)
    args = parser.parse_args()
    dpath = args.database_dir
    tpath = args.auxtables_dir

    #print ('THE DPATH IS', dpath)
    if not dpath:
        dpath = '/raid60/scratch/leo/scratch/era5/odbs/'
   
    dpath = '/raid60/scratch/leo/scratch/era5/odbs/' 
    if not tpath:
        tpath = '../data/tables/'
   
    print ('Analysing the databases stored in ', dpath)
    cdmpath='https://raw.githubusercontent.com/glamod/common_data_model/master/tables/'                                                                                                                                                                                               


    '''
    if len(sys.argv)==1:
        os.chdir('/'.join(sys.argv[0].split('/')[:-1]))
     
        dpath='../data/'
        tpath='../data/tables/'
        cdmpath='https://raw.githubusercontent.com/glamod/common_data_model/master/tables/'
    else:
        dpath=os.path.expandvars('$RSCRATCH/era5/odbs/')
        tpath=os.path.expanduser('~/tables/')
        #cdmpath=os.path.expandvars('$RSCRATCH/common_data_model/tables/')
    '''

    cities=pd.read_csv(tpath+'worldcitiespop.csv')
    #cities10000=cities.iloc[numpy.where(cities.Population>10000)[0]]
    #cities10000.to_csv('~/tables/worldcitiespop.csv')

    cdmtablelist=['id_scheme','crs','station_type','observed_variable','station_configuration_codes']        
    #os.chdir(os.path.expandvars('$RSCRATCH/common_data_model/tables'))
    #flist=glob.glob('*.dat')
    cdm=dict()
    for key in cdmtablelist:
        f=urllib.request.urlopen(cdmpath+key+'.dat')
        col_names=pd.read_csv(f,delimiter='\t',quoting=3,nrows=0)
        f=urllib.request.urlopen(cdmpath+key+'.dat')
        tdict={col: str for col in col_names}
        cdm[key]=pd.read_csv(f,delimiter='\t',quoting=3,dtype=tdict,na_filter=False)
#        print(cdm[key])
    #os.chdir(os.path.expandvars('$RSCRATCH/common_data_model/table_definitions'))
    #flist=glob.glob('*.csv')
    cdmd=dict()
    cdmtabledeflist=['id_scheme','crs','station_type','observed_variable','station_configuration']        
    for key in cdmtabledeflist:
        url='table_definitions'.join(cdmpath.split('tables'))+key+'.csv'
        f=urllib.request.urlopen(url)
        col_names=pd.read_csv(f,delimiter='\t',quoting=3,nrows=0,comment='#')
        f=urllib.request.urlopen(url)
        tdict={col: str for col in col_names}
        cdmd[key]=pd.read_csv(f,delimiter='\t',quoting=3,dtype=tdict,na_filter=False,comment='#')
#        print(cdmd[key])

    id_scheme={cdmd['id_scheme'].element_name.values[0]:[0,1,2,3,4,5,6],
               cdmd['id_scheme'].element_name.values[1]:['WMO Identifier','Volunteer Observing Ships network code',
                                                             'WBAN Identifier','ICAO call sign','CHUAN Identifier',
                                                             'WIGOS Identifier','Specially constructed Identifier']}
    cdm['id_scheme']=pd.DataFrame(id_scheme)
    cdm['id_scheme'].to_csv(tpath+'/id_scheme_ua.dat')
    cdm['crs']=pd.DataFrame({'crs':[0],'description':['wgs84']})
    cdm['crs'].to_csv(tpath+'/crs_ua.dat')
    cdm['station_type']=pd.DataFrame({'type':[0,1],'description':['Radiosonde','Pilot']})
    cdm['station_type'].to_csv(tpath+'/station_type_ua.dat')
    cdm['observed_variable']=pd.DataFrame({'variable':[1,2,3,4,7,29,59,111,112],
                                           'parameter_group':['geopotential','temperature','wind','wind','humidity',
                                                              'humidity','humidity','wind','wind'],
                                           'subdomain':['upper_air','upper_air','upper_air','upper_air','upper_air','upper_air','upper_air','upper_air','upper_air'],
                                           'name':['geopotential','air_temperature','eastward_wind','northward_wind','specific_humidity',
                                                   'relative_humidity','dew_point_temperature','wind speed','wind_from_direction'],
                                           'units':['m2 s-2','K','m s-1','m s-1', '1','1','K','m s-1','degree'],
                                           'description':['Geopotential is the sum of the specific gravitational potential energy relative to the geoid and the specific centripetal potential energy.',
                                                          'Air temperature is the bulk temperature of the air, not the surface (skin) temperature.',
                                                          '"Eastward" indicates a vector component which is positive when directed eastward (negative westward). Wind is defined as a two-dimensional (horizontal) \
                                                          air velocity vector, with no vertical component. (Vertical motion in the atmosphere has the standard name upward_air_velocity.)',
                                                                                                                                                                                          '"Northward" indicates a vector component which is positive when directed northward (negative southward). Wind is defined as a two-dimensional (horizontal) \
                                                          air velocity vector, with no vertical component. (Vertical motion in the atmosphere has the standard name upward_air_velocity.)' ,
                                                                                                                                                                                           '"specific" means per unit mass. Specific humidity is the mass fraction of water vapor in (moist) air.',
                                                          'actual water vapor pressure devided by saturation vapor pressure over water times hundred',
                                                          'Dew point temperature is the temperature at which a parcel of air reaches saturation upon being cooled at constant pressure and specific humidity.',
                                                          'Speed is the magnitude of velocity. Wind is defined as a two-dimensional (horizontal) air velocity vector, with no vertical component. \
                                                          (Vertical motion in the atmosphere has the standard name upward_air_velocity.) The wind speed is the magnitude of the wind velocity.',
                                                                                                                                                                                               'Wind is defined as a two-dimensional (horizontal) air velocity vector, with no vertical component. \
                                                          (Vertical motion in the atmosphere has the standard name upward_air_velocity.) In meteorological reports, \
                                                          the direction of the wind vector is usually (but not always) given as the direction from which it is blowing (wind_from_direction) \
                                                          (westerly, northerly, etc.). In other contexts, such as atmospheric modelling, it is often natural to give the direction\
                                                          in the usual manner of vectors as the heading or the direction to which it is blowing (wind_to_direction) (eastward, southward, etc.) \
                                                          "from_direction" is used in the construction X_from_direction and indicates the direction from which the velocity vector of X is coming.' 
                                                          ]})
    cdm['observed_variable'].to_csv(tpath+'/observed_variable_ua.dat')

    #os.chdir(os.path.expanduser('~/tables'))
    url='ftp://ftp.ncdc.noaa.gov/pub/data/igra/igra2-station-list.txt'
    urllib.request.urlretrieve(url,tpath+url.split('/')[-1])
    igrainv=pd.read_fwf(tpath+url.split('/')[-1],widths=(2,1,3,5,9,10,7,4,30,5,5,7),
                        names=('CC','Network','Code','StationId','Latitude','Longitude','Elev','dummy','StationName','From','To','Nrec'))

    #f=urllib.request.urlopen(url)
    #igrainv=pd.read_fwf(f,widths=(2,1,3,5,9,10,7,4,30,5,5,7),
                        #names=('CC','Network','Code','StationId','Latitude','Longitude','Elev','dummy','StationName','From','To','Nrec'))
    idx=numpy.where(igrainv.Latitude<-90.)
    igrainv.Latitude.iloc[idx]=numpy.nan
    igrainv.Longitude.iloc[idx]=numpy.nan
    igrainv.name='IGRA2'
    wbaninv=pd.read_fwf(tpath+'WBAN.TXT-2006jan',widths=(10,5,6,17,22,39,31,31,8,10,10,10,8,7),
                        names=('dum1','StationId','WMO ID','dum0','Country','dum2','StationName','dum3','From','To',
                               'Latitude','Longitude','dum4','Elev'))
    wblat=[]
    wblon=[]

    for slon,slat in zip(wbaninv.Longitude.values,wbaninv.Latitude.values):
        if type(slon)==float:
            wblon.append(slon)
            wblat.append(slat)
        else:
            llon=slon.split()
            llat=slat.split()
            wblon.append(abs(float(llon[0]))+float(llon[1])/60.+float(llon[2])/3600.)
            wblat.append(abs(float(llat[0]))+float(llat[1])/60.+float(llat[2])/3600.)
            if slon[0][0]=='-':
                wblon[-1]=-wblon[-1]
            if slat[0][0]=='-':
                wblat[-1]=-wblat[-1]
            
    wbaninv['Longitude'][:]=wblon[:]
    wbaninv['Latitude'][:]=wblat[:]    
    wbaninv['Elev'][:]=wbaninv['Elev'][:]*0.3048
    wbaninv['Elev'][wbaninv['Elev']==-99999.0*0.3048]=numpy.nan
    wbaninv.name='WBAN'

    lot2= pd.read_excel(tpath+'UA_Inventory_field_explanation.xlsx', sheet_name=None)
    unified=pd.DataFrame(columns=lot2['all_column_inventory']['Column_header'])
    col_names=pd.read_csv(os.path.expanduser(tpath+'uao_data.csv'),sep=';',nrows=0)
    trans=uai_trans(list(col_names),list(unified))

    if True:
        # TODO: replace with tpath tpath+'Inventory_ERACLIM_upperair_2.1.txt' 
        with open(os.path.expanduser(tpath+'Inventory_ERACLIM_upperair_2.1.txt'),'rb') as f:
            rdata=f.read().decode('ISO-8859-1').split('\n')
            z=[[] for _ in rdata[0].split('\t')]
            for r in rdata[1:-1]:
                li=r.split('\t')
#                print(len(li))
                if len(li)>=66:
                    for l in range(66): #len(li)):
                        z[l].append(li[l])

            chuan=dict(zip(rdata[0].split('\t'),z))
            chuan=pd.DataFrame(chuan)
            chuan=chuan.rename(columns = {'unique_record_ID':'StationId','WMO#':'WMO ID','Stationname':'StationName',
                                          'Lon_DegE':'Longitude', 'Lat_DegN':'Latitude', 'Alt_masl':'Elev'})
            chuan.name='CHUAN'
        url='https://oscar.wmo.int/oscar/vola/vola_legacy_report.txt'  
        urllib.request.urlretrieve(url,tpath+url.split('/')[-1])
        vola=pd.read_csv(tpath+'vola_legacy_report.txt',sep='\t')
        vlat=[]
        vlon=[]
        for v in vola['Latitude'].values:
            dms=v[:-1].split()
            if v[-1]=='S':
                vlat.append(-(float(dms[0])+float(dms[1])/60+float(dms[2])/3600))
            else:
                vlat.append(float(dms[0])+float(dms[1])/60+float(dms[2])/3600)
        for v in vola['Longitude'].values:
            dms=v[:-1].split()
            if v[-1]=='W':
                vlon.append(-(float(dms[0])+float(dms[1])/60+float(dms[2])/3600))
            else:
                vlon.append(float(dms[0])+float(dms[1])/60+float(dms[2])/3600)
        vola['Latitude']=numpy.array(vlat)
        vola['Longitude']=numpy.array(vlon)
        vola.name='OSCAR'
        os.chdir(dpath)
        row=['']*len(unified.columns)
        ci=False
        #func=partial(odb_process,ci,unified,vola,chuan,igrainv,wbaninv,trans)
        func=partial(odb_cdm,ci,cdm,cdmd,unified,vola,chuan,igrainv,wbaninv,trans)
        bfunc=partial(read_bufr_stn_meta,2)
        rfunc=partial(read_rda_meta) 
        tu=dict()
        p=Pool(20)

        dbs=['igra2','ai_bfr','rda','3188','1759','1761']
        dbs=['1759','1761']

        dbs=['3188']

        dbs=['ai_bfr'] 
        dbs = ["igra2"]        
        for odir in dbs:
            
            try: # what is this with the orphans ?
                with open(dpath+odir+'/orphans.csv','w'): 
                    pass
            except:
                pass
            if 'ai' in odir:
                flist=glob.glob(odir+'/'+'era5.*.bfr')
                #transunified=list(map(bfunc,flist))
                transunified=list(p.map(bfunc,flist))                                                                                                             
                print('******* Running ' , odir , '    ', len(flist) , ' files ')


            elif 'rda' in odir:
                flist=glob.glob(odir+'/'+'UADB_[tw]*.nc')
                flist.sort()
                glist=[]
                for i in range(len(flist)-1,-1,-1):
                    s=flist[i][-9:-3]
                    if s in glist:
                        del flist[i]
                    else:
                        glist.append(s)
                        
                #read_rda_meta(flist[0])
                #a = flist[:20]
                #for f in a:
                #      read_rda_meta(f)
                flist = flist[:20]      
                transunified=list(p.map(read_rda_meta, flist))
                    
            elif 'igra2' in odir:
                digrainv=igrainv.to_dict('records')
                #for d in digrainv:
                #      read_igra_meta(d)
                transunified=list(p.map(read_igra_meta,digrainv))

            else:
                if odir == '3188':
                    #code.interact(local = locals())
                    flist = glob.glob(odir+'/'+'era5.3188.*') #era5.3188.conv.C:4953
                    #flist = [ f for f in os.listdir(dpath + '/' + odir) if 'era5.3188.conv' in f ]
                elif odir == '1759' :
                    flist=glob.glob(odir+'/'+'era5.1759.conv.*') 
                    print(flist)
                elif odir == '1761':
                    flist=glob.glob(odir+'/'+'era5.1761.conv.*')
                elif odir == '2':
                    flist = glob.glob(odir+'/'+'era5.conv_*')
                elif odir == '1':
                    flist = glob.glob(odir+'/'+'era5.conv._*')

                hlist=[]
                
                for f in flist:
                    if '.nc' not in f and '.gz' not in f and '00000' not in f and '99999' not in f:
                        hlist.append(f)

                print('******* Running ' , len(hlist) , ' files ')
                #print(hlist[:100])
                #hlist = hlist[:10]
                #hlist =['3188/era5.3188.conv.1:6607']
                #print('HLIST', hlist )

                transunified=list(p.map(func,hlist))

            for i in range(len(transunified)-1,-1,-1):
                if type(transunified[i]) is not tuple:
                    del transunified[i]
            
            tucdmlist=[]
            tu311alist=[]
            for tucdm,tu311a in transunified:
                tucdmlist.append(tucdm)
                tu311alist.append(tu311a)
           
           
            #output in CDM format                
            transunified= list(filter(None, tucdmlist))           

            print('Processing the dataframe ****')
            v = {k: [dic[k] for dic in transunified] for k in transunified [0]}
            v['record_number']=list(range(len(v['record_number'])))
            transunified=v
            
            tu[odir]=pd.DataFrame(transunified)
            tu[odir].name= odir
            tu[odir].to_csv( home + '/' + odir + '_meta.csv',na_rep='NA',sep='\t',index=False) # changed to write in this directory 
            #os.system('mv' + home + '/summary.txt ' + home + '/' + odir + '_file_summary.txt' )
    
            #output in 311a legacy format                
            transunified= list(filter(None, tu311alist))     
    
            tu[odir]=pd.DataFrame(transunified,columns=unified.columns.values)
            tu[odir].name=odir
            tu[odir].to_csv( home + '/' + odir + '_meta_311a.csv',na_rep='NA',sep='\t',index=False)

        '''
        dbs=['3188','1759','1761']
        plot_orphans(dbs,path=dpath)
        dbs=['igra2','ai_bfr','rda','3188','1759','1761']
        plot_nowigos(dbs,path=dpath)
        plot_numbers(dbs,path=dpath)
        for dec in range(1940,2020,10):
            
            plot_active(dbs,[str(dec),str(dec+10)],path=dpath)

        print('ready')
        '''
