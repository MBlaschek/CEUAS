#!/usr/bin/env python
import sys
import os.path
import glob

import subprocess
import urllib.request
import xarray as xr
import numpy
import datetime
from datetime import date, datetime
import time
from multiprocessing import Pool
from netCDF4 import Dataset
import gzip
import pandas as pd
from functools import partial
#from utils import fdist
from eccodes import *
import matplotlib.pylab as plt
import cartopy.crs as ccrs
import argparse
sys.path.append("/fio/srvx7/leo/python/Rasotools/")
import math
#from numba.experimental import jitclass
#from numba import jitclass 
#from utils import *   


# need this since during the running, it moves to other directories 
home = os.getcwd()
out_dir = home +'/output_data/'

if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

def fdist(slat,slon,lat,lon):
    if type(slon) is str:
        if '.' in slon:
            flon=float(slon)
            flat=float(slat)
        elif ',' in slon:
            flon=float('.'.join(slon.split(',')))
            flat=float('.'.join(slat.split(',')))
        elif slon[-1] in 'WE':
            flon=float(slon[:3])+float(slon[4:7])/60+float(slon[-2:-1])/3600
            if slon[-1]=='W':
                flon=-flon
            flat=float(slat[:2])+float(slat[3:6])/60+float(slat[-2:-1])/3600
            if slat[-1]=='S':
                flat=-flat
        else:
            try:

                flon=float(slon[-9:-6])+float(slon[-5:-3])/60+float(slon[-2:])/3600
                if slon[0]=='-':
                    flon=-flon
                flat=float(slat[-8:-6])+float(slat[-5:-3])/60+float(slat[-2:])/3600
                if slat[0]=='-':
                    flat=-flat
            except:
                flon=numpy.nan
                flat=numpy.nan

    else:
        flon=slon
        flat=slat


    lats=numpy.append(flat,lat) # TODO changed this 
    lons=numpy.append(flon,lon)

    x=numpy.cos(lats*math.pi/180.)*numpy.cos(lons*math.pi/180.)
    y=numpy.cos(lats*math.pi/180.)*numpy.sin(lons*math.pi/180.)
    z=numpy.sin(lats*math.pi/180.)

    dist=numpy.arccos((x[0]*x[1:]+y[0]*y[1:]+z[0]*z[1:])*0.9999999)
    
    return dist


def get_secondary_nrda(file):
    """ This function extracts the complete list of secondary station ids from the original ncar files (not available from the netCDF) """
    ncardir = '/raid60/scratch/federico/databases/UADB/'
    tr = 'uadb_trhc_'
    wind = 'uadb_windc_'
    stat_in_name = file.split('_')[-1].replace('.nc','')[1:]
    # the station numbe rin the nc files alwayas have 6 digits while in the real txt files might have 4 or 5. Some exra 0s have been added.

    if 'trhc' in file:
        file_path = ncardir + tr 
    else:
        file_path = ncardir + wind 

    if stat_in_name == '10185':
        pass 
    if os.path.isfile(file_path + stat_in_name + '.txt' ):   
        file_path = file_path + stat_in_name + '.txt'
    elif  os.path.isfile(file_path + stat_in_name[1:] + '.txt'  ):
        file_path = file_path + stat_in_name[1:] + '.txt'
    else:
        print('File not found! ' , file )

    data = open(file_path , 'r').readlines()

    stations_id = []
    for i, line in enumerate(data):
        if line[0] == 'H':
            #usi      = int(line[2:14])  # unique station identifier
            ident    = line[15:21].replace(' ','')# WMO
            if ident not in stations_id:
                stations_id.append(ident)
    stat_ids = ','.join(stations_id)

    return stat_ids



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

def find_wban(lat,lon,stat_id):
    """ checks if station_id with  swapped lat and lon is found in the wban """
    z=numpy.where(wbaninv['StationId']==stat_id)[0]
    if len(z) >0:
        z = z[0]      
        dist = fdist(wbaninv['Latitude'][z],wbaninv['Longitude'][z],-lat,lon)*180./math.pi
        if dist < 0.5:
            return True
        else:
            return False
    else:
        return False

def do_row(row,kdict,meta,wmoid,iwmoid, nw,lat,lon,name, fn='', dataset='') :
    # meta is the data source (as a pd dataframe) of the metadata e.g. vola, wban etc. 

    inventory_codes = { 'OSCAR' : '20000',
                                      'IGRA2':'20300',
                                      'WBAN': '20500',
                                      'CHUAN':'20400',
                                      # unidentified: '20600',
                                      }

    try_wban = False 
    if dataset == '1759': # if True, will try for a latitude mismatch inside the OSCAR dataset 
        try_wban = find_wban(lat,lon, wmoid)

    try:

        z = []
        meta.name=name
        print('Testing the ' , meta.name , '  inventory ')
        # preliminary chekc if the station id is the in OSCAR dataset 
        if meta.name=='OSCAR':
            z=numpy.where(meta['StationId']=='0-20000-0-{:0>5}'.format(int(iwmoid)))[0]
        else:
            z=numpy.where(meta['StationId']==wmoid)[0]


        if len(z) !=0: # stationId found in the inventory 
            ind=z[0]
            metaz=meta.loc[ind]
            # check position mismatch: just by looking at the stationid, I might find wrong coordinates 

            dist = fdist(metaz['Latitude'],metaz['Longitude'],lat,lon)*180./math.pi 

            if dist > 1.0:
                print('StationId ' , wmoid , 'found in the ' , meta.name , ' report as ' , str(iwmoid) , ' but coordinates do not match')
                print('Lat, lon' , lat, lon , '    ', meta.name , ' lat,lon' , metaz['Latitude'],metaz['Longitude'] )
                found = False

            else: # if also coordinates match, keep the record found

                print(meta.name , ' coordinates match do_row for station: ' , wmoid, ' ' , meta['StationId'][z[0]], ' ' , meta['StationName'][z[0]])

                #print('Found: ', wmoid ,' inside the ', meta.name+' database with id' , meta['StationId'][z].values[0] )
                row[kdict['stationPlatformUniqueIdentifier']]=metaz['StationId']
                row[kdict['station name']]=metaz['StationName']
                row[kdict['stationDescription']]=metaz['StationName']
                # will use the file lat/lon and not the coords from the inventory item


                if meta.name == 'OSCAR':

                    row[kdict['regionOfOriginOfData']]=metaz['RegionName']
                    row[kdict['territoryOfOriginOfData_Name']]=metaz['CountryArea']
                    row[kdict['territoryOfOriginOfData_ISO3CountryCode']]=metaz['CountryCode']  
                    found = True

                    if 'Upper' not in vola['ObsRems'].values[ind] : # 
                        found = False           
                    else:
                        return ind  # if found = False, will not be used anyway 

                else:

                    if meta.name=='CHUAN': 
                        row[kdict['descriptionDataset']]=metaz['Source']

                    found = True
                    return ind


        if len(z)==0 or not found: # if no station id is found or coords do not match, check for matching coordinates, take the closest station 
            dists = fdist(lat, lon, meta['Latitude'],meta['Longitude'])*180./math.pi # !! DONOT invert the order of the lats and longs station / inventory !!! 
            print('*** Checking for matching coordinate for unknown station id in inventory: ' , meta.name )
            if meta.name =='OSCAR': # for OSCAR I have to loop over the possible stations since not all are upper air, and find the correct ones 
                idy=numpy.argsort(dists) # Returns the indices that would sort an array
                i=0
                oscar_found = False
                standard_wigos_idx = ''
                
                while dists[idy[i]]<0.5:
                    

                    
                    oscar_found = True
                    idx=idy[i]
                    
                    # here I select, among the matching entries, the one with the normal WIGOS id 
                    # e.g. I select Africa	Zambia	ZMB	0-20000-0-67575	67575	
                    # instead of Africa	Zambia	ZMB	0-894-1-mks0025
                    # that have the same lat and lon 
                    
                    if str(iwmoid ) in meta['StationId'][idx] and not standard_wigos_idx :
                        standard_wigos_idx = idx

                    row[kdict['station name']]=meta.StationName[idx]
                    row[kdict['stationDescription']]=meta.StationName[idx]
                    row[kdict['regionOfOriginOfData']]=meta.RegionName[idx]
                    row[kdict['territoryOfOriginOfData_Name']]=meta.CountryArea[idx]
                    row[kdict['territoryOfOriginOfData_ISO3CountryCode']]=meta.CountryCode[idx]  

                    #cc= vola['CountryCode'][idx]
                    #if 'Upper' in vola['ObsRems'].values[idx] or '0-2000' in vola['StationId'][idx]:
                    #print(vola.iloc[idx], '   ' ,  dists[idy[i]] )
                    
                    if 'Upper' in vola['ObsRems'].values[idx] :
                        print('Found an OSCAR upper air station with compatible coordinates ') 
                        row[kdict['stationPlatformUniqueIdentifier']]=meta.StationId[idx]                        
                        return idx

                    i+=1

                if oscar_found: # only reached if return idx is not executed 
                    print('Found an OSCAR station with compatible coordinates (but no Upper air/Radiosonde station)')                                
                    
                    # return the clostes entry with a regular WIGOS id, otherwise simply the closest entry 
                    if standard_wigos_idx:
                        row[kdict['stationPlatformUniqueIdentifier']]=meta.StationId[standard_wigos_idx]                        
                        return standard_wigos_idx
                    else:
                        row[kdict['stationPlatformUniqueIdentifier']]=meta.StationId[idy[0]]                                                
                        return idy[0]

                if try_wban: # it means that the id was found inside the wban with a latitude sign swap 
                    dists = fdist(-lat, lon, meta['Latitude'],meta['Longitude'])*180./math.pi # !! DONOT invert the order of the lats and longs station / inventory !!! 
                    idy=numpy.argsort(dists) # Returns the indices that would sort an array
                    i=0
                    oscar_found = False
                    while dists[idy[i]]<0.5:

                        print("Found a station in the WBAN with a latitude mismatch that also appears in the OSCAR dataset! ")

                        idx=idy[i]

                        row[kdict['station name']]=meta.StationName[idx]
                        row[kdict['stationPlatformUniqueIdentifier']]=meta.StationId[idx]
                        row[kdict['stationDescription']]=meta.StationName[idx]
                        row[kdict['regionOfOriginOfData']]=meta.RegionName[idx]
                        row[kdict['territoryOfOriginOfData_Name']]=meta.CountryArea[idx]
                        row[kdict['territoryOfOriginOfData_ISO3CountryCode']]=meta.CountryCode[idx]  

                        row[kdict['lat']] = -lat

                        if not os.path.isfile(out_dir + '/1759_WBAN_mismatch_OSCAR_found.dat'):
                            a = open(out_dir + '/1759_WBAN_mismatch_OSCAR_found.dat', 'w')
                            a.write('primary_id' + '\t' + 'file_lat' + '\t' + 'wban_lat' + '\t' + 'file_lon' + '\t' + 'wban_lon' + '\t' + 'file' + '\n')
                            a.close()
                            
                        a = open(out_dir + '/1759_WBAN_mismatch_OSCAR_found.dat', 'a+')
                        #a.write(fn + '\n')
                        a.write(meta.StationId[idx] + '\t' + str(lat) + '\t' + str(meta['Latitude'][idx]) + '\t' + 
                                                    str(lon) + '\t' + str(meta['Longitude'][idx]) + '\t' + fn.split('/')[1] + '\n')                        
                        a.close()
                        return idx

            else:

                z=numpy.nanargmin(dists)
                if dists[z]<0.5:
                    
                    z=[z]
                    metaz=meta.loc[z]
                    print(meta.name , ' coordinates match do_row  ' , wmoid, ' ' , metaz['StationId'].values[0], ' ' , metaz['StationName'].values[0])
                    if len(z) > 0:
                                          
                        row[kdict['station name']]=metaz['StationName'].values[0]
                        row[kdict['stationDescription']]=metaz['StationName'].values[0]
                        #print('Found: ', wmoid ,' inside the ', meta.name+' database with id' , meta['StationId'][z] )
                        row[kdict['stationPlatformUniqueIdentifier']]=metaz['StationId'].values[0]
                        row[kdict['station name']]=metaz['StationName'].values[0]
                        row[kdict['stationDescription']]=metaz['StationName'].values[0]

                        if meta.name=='CHUAN': 
                            row[kdict['descriptionDataset']]=metaz['Source'].values[0] 

                        return z[0]
                    else:
                            
                        return False

            # testing missing latitude; only works for era5_1759 
            #if meta.name == 'WBAN' and dataset == '1759':
            if  dataset == '1759' and meta.name == 'WBAN':
                z=numpy.where(meta['StationId']==wmoid)[0]
                if len(z) > 0:

                    z=numpy.where(meta['StationId']==wmoid)[0]

                    print('Testing possible missing latitude sign === ')
                    if lat > 0:
                        if fdist(metaz['Latitude'],metaz['Longitude'],-lat,lon)*180./math.pi < 0.5:
                            '0-' + inventory_codes[meta.name] + '0-' + str(iwmoid)
                            #print('StationId ' , wmoid , 'found in the ' , meta.name , ' report as ' , '0-20400-0-{:0>5}'.format(int(iwmoid)) , ' when swapping latitude')

                            print('StationId ' , wmoid , 'found in the ' , meta.name , ' report as ' , '0-' + inventory_codes[meta.name] + '0-' + str(iwmoid) , ' when swapping latitude')

                            print('! FOUND a latitude swap:  file -Lat, lon' , -lat, lon , '    ', meta.name , ' lat,lon' , metaz['Latitude'],metaz['Longitude'] )    
                            print(meta.name , ' coordinates match do_row for station: ' , wmoid, ' ' , meta['StationId'][z[0]], ' ' , meta['StationName'][z[0]])

                            if not os.path.isfile(out_dir + '/era5_1759_WBAN_latitude_mismatch.dat'):
                                a = open(out_dir + '/era5_1759_WBAN_latitude_mismatch.dat', 'w')
                                a.write('primary_id' + '\t' + 'file_lat' + '\t' + 'wban_lat' + '\t' + 'file_lon' + '\t' + 'wban_lon' + '\t' + 'file' + '\n')
                                a.close()

                            a = open(out_dir + '/era5_1759_WBAN_latitude_mismatch.dat', 'a+')
                            a.write(str(iwmoid) + '\t' + str(lat) + '\t' + str(metaz['Latitude']) + '\t' + 
                                                        str(lon) + '\t' + str(metaz['Longitude']) + '\t' + fn.split('/')[1] + '\n')

                            row[kdict['station name']]=metaz['StationName']
                            row[kdict['stationDescription']]=metaz['StationName']
                            row[kdict['stationPlatformUniqueIdentifier']]=metaz['StationId']
                            row[kdict['station name']]=metaz['StationName']
                            row[kdict['stationDescription']]=metaz['StationName']
                            row[kdict['lat']] = -lat
                            row[kdict['lon']] = lon                            
                            #row[kdict['latitude']]= -lat  # writing correct sign latitude

                            return ind
                        else:
                            found = False 
                    else:
                        found = False
                else:
                    found = False

    except KeyError as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno, 'FAILED do row: ' , wmoid)


def insert(row,vola,chuan,igrainv,wbaninv,transodb,trans,kdict,unified,rdata,varno,dmin,dmax,dcount,fi,fn):

    dataset= fn.split('/')[0]

    ids=-1
    wmoid=''
    iwmoid=0
    mi = False         

    if not fi:
        uai=dict(zip(rdata[0].split(),zip(rdata[1].split(),rdata[-2].split())))
        uai['statid@hdr']=list(uai['statid@hdr'])
        if uai['statid@hdr'][0]=="''":
            for i in range(len(uai['statid@hdr'])):
                uai['statid@hdr'][i]="'C:"+fn.split(':')[-1]
        vlist=list(transodb.values())

        for k,v in transodb.items():
            for vv in v:
                try:
                    if row[kdict[vv]]=='':
                        row[kdict[vv]]=uai[k][0]
                    else:
                        if vlist.count(v)>1:
                            if uai[k][0]!=row[kdict[vv]]:
                                row[kdict[vv]]+=' , '+uai[k][0]
                except:
                    pass
            if k=='statid@hdr':   # extract information from vola data base
                if ':' in uai['statid@hdr'][0]:
                    try:
                        nw=uai['statid@hdr'][0][1:-1].split(':')[0]
                        if nw=='C':
                            wmoid=fn.split(':')[-1]
                            if wmoid[-1] in 'AB':
                                wmoid=wmoid[:-1]
                            if wmoid[-1] in 'AB':
                                wmoid=wmoid[:-1]                                
                        else:
                            wmoid=uai['statid@hdr'][0][1:-1].split(':')[1]

                        lat=float(uai['lat@hdr'][0])
                        lon=float(uai['lon@hdr'][0])
                        iwmoid=int(wmoid)
                    except:
                        print('error parsing wmo ID',wmoid,uai['lat@hdr'][0],uai['lon@hdr'][0])
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
                            # error case: W5331. This will break the code. I then remove the non digits characters
                        try:                              
                            iwmoid=int(wmoid)
                        except:
                            import re
                            iwmoid = int(re.sub('\D', '', wmoid ))

                #if iwmoid>=1 and iwmoid<1000000:
                if True:
                
                    mi=do_row(row,kdict,vola,wmoid, iwmoid, nw,lat,lon,'OSCAR', dataset = dataset, fn = fn) # mi is the number of the row inside the OSCAR database for the station
                    if not mi:
                        mi=do_row(row,kdict,igrainv,wmoid, iwmoid,nw,lat,lon,'IGRA2')
                        if not mi:

                            mi=do_row(row,kdict,wbaninv,wmoid, iwmoid, nw,lat,lon,'WBAN', dataset = dataset, fn = fn )
                            if not mi:
                                mi=do_row(row,kdict,chuan,wmoid, iwmoid, nw,lat,lon,'CHUAN', fn = fn)
                                if not mi:
                                    print(wmoid,'not found, not identifiable!')
                                    #a = open(out_dir + '/' + fn.split('/')[0] + '_' + 'orphans_not_identified.dat','a')
                                    #a.write(fn + '\n' )
                                    filepath = out_dir + '/' + fn.split('/')[0] + '_orphans.csv'

                                    with open(filepath,'a') as f:
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
                                    print('FOUND CHUAN' , fn )                                                    
                            else:
                                ids=2
                        else:
                            if igrainv['Network'].values[mi]=='M':
                                print('IGRA2 WMO station , set ids=0 ')
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

    # cc=row[kdict['territoryOfOriginOfData_ISO3CountryCode']]
    #a = open(home + '/' + fn.split('/')[0] + '_all_files.dat' , 'a')
    #a.write(fn + ',' + wmoid + ',' + row[kdict['stationPlatformUniqueIdentifier']] + '\n')
    return ids, wmoid, row[kdict['territoryOfOriginOfData_ISO3CountryCode']] , row[kdict['stationPlatformUniqueIdentifier'] ] , mi 

def odb_process(ci,unified,vola,chuan,igrainv,wbaninv,trans,fn):           

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
                print(e, 'odb call failed')

    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno , 'odb_process fail' )
        a = open(out_dir + '/failed_odb.dat', 'a' )

        exit()

    if len(transunified)==0:
        transunified.append(row[:])

    return transunified[0]

def scini(station_configuration,ids,station_id,cwmoid,vola,transunified,kdict,dmax,dmin,cdm,cities,cc, mi, fn=''):

    try:
        station_configuration['longitude']=float(transunified[-1][kdict['lon']])
        station_configuration['latitude']=float(transunified[-1][kdict['lat']])
    except:
        station_configuration['longitude']=numpy.nan
        station_configuration['latitude']=numpy.nan
    station_configuration['station_name']=transunified[-1][kdict['station name']]

    station_configuration['secondary_id']=cwmoid
    station_configuration['primary_id_scheme']= ids

    #station_configuration['secondary_id_scheme']=ids

    """ If ids==0, the station id was found either in the OSCAR or in the IGRA inventory but from a  WMO network station """
    if ids==0:
        station_configuration['primary_id']='0-20000-0-'+ station_id
        if '0-2000' not in station_id:
            station_id =  '0-20000-0-' + station_id 
        station_configuration['primary_id']= station_id
        station_configuration['primary_id_scheme']=5

    elif ids == 3: # case of IGRA station not from a WMO network, use full IGRA id 
        igrainv_row = igrainv.iloc[mi]
        try:
            station_configuration['primary_id']= '0-20300-0-' + igrainv_row['igraId'].values[0]
        except:
            station_configuration['primary_id']= '0-20300-0-' + igrainv_row['igraId']

        station_configuration['primary_id_scheme']=3

    elif ids ==2: # case of WBAN station
        wbaninv_row = wbaninv.iloc[mi]
        try:
            station_configuration['primary_id']= '0-20500-0-' + str(wbaninv_row['StationId'].values[0])
        except:
            station_configuration['primary_id']= '0-20500-0-' + str(wbaninv_row['StationId'])                    

        station_configuration['primary_id_scheme']=2

        #vola,chuan,igrainv,wbaninv
        print(0)

    elif ids==4: # case of CHUAN station
        print ('CHUAN ************************** ' , ids, station_id, cwmoid )
        chuan_row = chuan.iloc[mi]
        station_configuration['primary_id']= '0-20400-0-' + str(chuan_row['StationId'])
        station_configuration['primary_id_scheme']=4
    else:

        print ('CANNOT FIND ************************** ' , ids, station_id, cwmoid )
        station_configuration['primary_id']= '0-20600-0-' + cwmoid
        b = open( out_dir + '/' + fn.split('/')[0] + '_20600.csv', 'a')
        b.write(fn + '\t' +  '0-20600-0-' + cwmoid + '\n')
        b.close()
        


    station_configuration['station_crs']=0
    station_configuration['local_gravity']=numpy.NaN
    try:
        if '-' in dmin:
            station_configuration['start_date']=datetime.strptime(dmin,'%Y-%m-%d')
            station_configuration['end_date']=datetime.strptime(dmax,'%Y-%m-%d')
        else:
            station_configuration['start_date']=datetime.strptime(dmin,'%Y%m%d')
            station_configuration['end_date']=datetime.strptime(dmax,'%Y%m%d')
    except:
        print('invalid date encountered:',dmin,dmax)
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

    return    

def odb_cdm(ci,cdm,cdmd,unified,vola,chuan,igrainv,wbaninv,trans,fn):           
    #print(fn)
    transunified=[]
    station_configuration=dict()
    for k in cdmd['station_configuration'].element_name.values:
        station_configuration[k]=list()
    fi=False
    #try:
    for varno in ['111','112','1','2','59','29']:
        qs='select min(date),max(date) where varno='+varno #distinct,count(date)
        try:
            rdata=subprocess.check_output(["odb","sql","-q",qs,"-i",fn,'--no_alignment'],stderr=subprocess.DEVNULL)
            dmin,dmax=rdata.decode().split('\n')[1].split()

            dcount=0
            if dmin!='NULL':
                if int(dmin)<19000000 or int(dmax)>21000000:
                    print(fn,': MALFORMED DATE ERROR!!')
                    qs='select min(date),max(date) where varno='+varno +' and date>19000000' #distinct,count(date)
                    rdata=subprocess.check_output(["odb","sql","-q",qs,"-i",fn,'--no_alignment'],stderr=subprocess.DEVNULL)
                    dmin,dmax=rdata.decode().split('\n')[1].split()
                qs='select * where varno='+varno+' and (date='+dmin+' or date='+dmax+')'
                rdata=subprocess.check_output(["odb","sql","-q",qs,"-i",fn,'--no_alignment'],stderr=subprocess.DEVNULL).decode().split('\n')
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
                                print(fn,' lat lon supplied from metadata on ',i,'lines')
                            except:
                                print(fn,' lat lon could not be supplied from metadata')
                                pass
                    transodb,kdict=odb_trans(cols,unified)
                if not fi:
                    transunified.append(row[:])
                    ids,wmoid,cc, primary_id, z =insert(transunified[-1],vola,chuan,igrainv,wbaninv,transodb,trans,kdict,unified,rdata,varno,dmin,dmax,dcount,fi,fn)
                    # ids,wmoid,cc, station_id =insert(transunified[-1],vola,chuan,igrainv,wbaninv,transodb,trans,kdict,unified,rdata,varno,dmin,dmax,dcount,fi,fn)

                    if ':' in fn:  
                        fnl=fn.split(':')
                        cwmoid=fnl[-2][-1]+':'+fnl[-1]
                    else:
                        cwmoid=wmoid

                    scini(station_configuration,ids,primary_id,cwmoid,vola,transunified,kdict,dmax,dmin,cdm,cities,cc,z, fn=fn)
                    


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
        print(fn,'insufficient data for assimilation')
        return None

    if 2 in station_configuration['observed_variables']:
        station_configuration['station_configuration_code']=13 # Radiosonde, taken from station_configuration_code cdm table
    else:
        station_configuration['station_configuration_code']=8 # Pilot, taken from station_configuration_code cdm table

    a = open(out_dir + '/' + fn.split('/')[0] + '_correctly_processed.dat' , 'a')
    #print('Done === ' , fn )
    a.write(fn + '\n')

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
            ids,wmoid,cc, primary_id, z =insert(transunified[-1],vola,chuan,igrainv,wbaninv,transodb,trans,kdict,unified,rdata,varno,dmin,dmax,dcount,fi,bufrfile)
            if not fi:
                scini(station_configuration,ids,primary_id,wmoid,vola,transunified,kdict,dmax,dmin,cdm,cities,cc, z, fn = bufrfile)
                fi=True
            station_configuration['observed_variables'].append(int(varno)) 

    #print ('/'.join(bufrfile.split('/')[-1:]),cnt,"messages",time.time()-tx)

    if len(transunified)==0:
        transunified.append(row[:])

    a = open(out_dir + '/' + bufrfile.split('/')[0] + '_correctly_processed.dat','a' )
    a.write(bufrfile + '\n' )

    return station_configuration,transunified[0]

def read_rda_meta(ncfile):

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
                transodb,kdict=odb_trans(cols,unified)  # this just maps the odb column names to the cdm column names 
                dcount=0
            if not fi:
                transunified.append(row[:])
            ids,wmoid,cc,primary_id, z=insert(transunified[-1],vola,chuan,igrainv,wbaninv,transodb,trans,
                                              kdict,unified,rdata,varno,dmin,dmax,dcount,fi,ncfile)
            if not fi:
                scini(station_conf,ids,primary_id,wmoid,vola,transunified,kdict,dmax,dmin,cdm,cities,cc,z, fn = ncfile)
                fi=True
            station_conf['observed_variables'].append(int(varno)) 

    #print ('/'.join(ncfile.split('/')[-1:]),time.time()-tx)

    if len(transunified)==0:
        transunified.append(row[:])

    secondary = get_secondary_nrda(ncfile)
    station_conf['secondary_id'] = secondary # updating the secondary with the complete list found in the file 
    return station_conf,transunified[0]

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

        ids,wmoid,cc,primary_id, z=insert(transunified[-1],vola,chuan,igrainv,wbaninv,transodb,trans,
                                          kdict,unified,rdata,varno,dmin,dmax,dcount,fi,line['StationId'])
        if not fi:
            scini(station_conf,ids,primary_id,wmoid,vola,transunified,kdict,dmax,dmin,cdm,cities,cc, z)
            fi=True
        station_conf['observed_variables'].append(int(varno)) 



    #print ('/'.join(line['StationId'].split('/')[-1:]),time.time()-tx)

    if len(transunified)==0:
        transunified.append(row[:])

    station_conf['secondary_id'] = line['igraId'] 
    station_conf['secondary_id_scheme'] = 3    
    return station_conf,transunified[0]


if __name__ == '__main__':

    print ('THE DPATH IS', dpath)
    if not dpath:
        dpath = '../data/'

    if not tpath:
        tpath = home + '/../data/tables/'

    print ('Analysing the databases stored in ', dpath)
    cdmpath='https://raw.githubusercontent.com/glamod/common_data_model/master/tables/'                                                                                                                                                                                               

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

    igrainv['igraId'] =  igrainv['CC'] + igrainv['Network'] + igrainv['Code'] + igrainv['StationId']
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

    lot2= pd.read_excel(tpath+'UA_Inventory_field_explanation.xlsx', sheet_name=None, engine ='openpyxl')
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
            chuan['Lat_DegN']  = [float(l.replace(',','.')) for l in chuan['Lat_DegN'] ]
            chuan['Lon_DegE'] = [float(l.replace(',','.')) for l in chuan['Lon_DegE']]

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
        p=Pool(30)
        
        dbs=['igra2','ai_bfr','rda','3188','1759','1761','1','2']
        dbs=['RI/nasa']
        dbs=['RI/Pangaea/COMP']


        PARALLEL = False  

        dbs=['ai_bfr','rda','3188','1761','2','1','igra2','1759']

    
        dbs=['igra2']
        
        transunified = []

        for odir in dbs: 
            #with open(dpath+odir+'/orphans.csv','w'):
            #    pass
            if 'ai' in odir:
                flist=glob.glob(odir+'/'+'era5.*.bfr')

                #flist = ['ai_bfr/era5.10106.bfr']

                flist = [f for f in flist if 'undef' not in f ] # clean some wrong files
                #proc = open(home + '/ai_bfr_correctly_processed.dat').readlines()
                #proc = [f.replace('\n','') for f in proc]

                if os.path.isfile(out_dir + '/' + odir +'_correctly_processed.dat'):
                    proc = open(out_dir + '/' + odir +'_correctly_processed.dat').readlines()
                    proc = [f.replace('\n','') for f in proc]
                else:
                    proc = []

                flist = [ f for f in flist if f not in proc ]
                #transunified=list(map(bfunc,flist))
                #flist = ['ai_bfr/era5.37260.bfr']
                #flist = ['ai_bfr/era5.70292.bfr']

                if PARALLEL:

                    transunified=list(p.map(bfunc,flist))
                else:
                    for f in flist:

                        transunified.append(read_bufr_stn_meta(2,f))

            elif 'rda' in odir:
                flist=glob.glob(odir+'/'+'UADB_[tw]*.nc')
                flist.sort()
                #f4list = [f for f in flist if '82599' in f ]
                glist=[]
                for i in range(len(flist)-1,-1,-1):
                    s=flist[i][-9:-3]
                    if s in glist:
                        del flist[i]
                    else:
                        glist.append(s)
                #transunified=list(p.map(read_rda_meta,flist))
                #transunified=[]
                print('running ', len(flist) , '  files')
                if not PARALLEL:
                    print('Running in single mode === ')                
                    for f in flist:
                        transunified.append(read_rda_meta(f))
                else:
                    transunified=list(p.map(read_rda_meta,flist))

            elif 'igra2' in odir:
                digrainv=igrainv.to_dict('records')
                transunified=list(map(read_igra_meta,digrainv))


            else:
                if odir == '1' or odir == '2':
                    flist=glob.glob(odir+'/'+'era5.conv._*')
                    #flist = ['1/era5.conv._73033']

                elif odir == '1759':
                    flist=glob.glob(odir+'/'+'era5.1759.conv.*')
                elif odir == '1761':
                    flist=glob.glob(odir+'/'+'era5.1761.conv.*')
                elif odir =='3188':
                    flist=glob.glob(odir+'/'+'era5.3188.conv.*')
                    #flist = [f for f in flist if '4581' in f ]

                flist =  [f for f in flist if 'gz' not in f and f != '3188/era5.3188.conv.' ]
                #flist =[f for  f in flist if '_10401' in f ]  # to test latitude mismatch in 1759

                #flist =[f for  f in flist if '2:61701' in f ]  # to test latitude mismatch in 1759
                #f = pd.read_csv(home + '/era5_1759_WBAN_latitude_mismatch_git.dat', delimiter = '\t')
                #flist = [ '1759/' + f for f in f['file'].values ]
                print(0)
                """
                if os.path.isfile(out_dir + '/' + odir +'_correctly_processed.dat'):
                    proc = open(out_dir + '/' + odir +'_correctly_processed.dat').readlines()
                    proc = [f.replace('\n','') for f in proc]
                else:
                    proc = []
                """
                proc = []
                
                #lista = pd.read_csv(home + '/output_data_SAVE/era5_1759_WBAN_latitude_mismatch.dat', delimiter = '\t')

                flist = [ f for f in flist if f not in proc ]

                hlist=[]
                for f in flist:
                    if '.nc' not in f and '.gz' not in f and '00000' not in f and '99999' not in f:
                        hlist.append(f)


                #hlist = hlist[:100]
                print('running ', len(hlist) , '  files')
                if not PARALLEL:
                    print('Running in single mode === ')
                    for f in hlist:
                        a = odb_cdm(ci,cdm,cdmd,unified,vola,chuan,igrainv,wbaninv,trans, f )
                        transunified.append(a)
                else:
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

            v = {k: [dic[k] for dic in transunified] for k in transunified [0]}
            v['record_number']=list(range(len(v['record_number'])))
            transunified=v

            tu[odir]=pd.DataFrame(transunified)
            tu[odir].name=odir
            tu[odir].to_csv(out_dir + '/' +  odir.replace('/','_') + '_meta.csv', na_rep='NA',sep='\t',index=False, mode = 'a')
            print('written stat_conf file')


