# -*- coding: utf-8 -*-
# this file is released under public domain and you can use without limitations

#########################################################################
## This is a simple controller
## - index is the default action of any application
##   here it forwards the request from a CDS server to further processing 
##   with python
##
##   For queries ask early-upper-air@copernicus-climate.eu
##   The C3S 311c Lot2 group
##   Vienna, 26 August 2019
## - user is required for authentication and authorization
## - download is for downloading files uploaded in the db (does streaming)
## - api is an example of Hypermedia API support and access control
#########################################################################

import os,sys
import urllib
import socket
#from gluon.debug import dbg
host=socket.gethostname()
print (host)
if 'srvx' in host:
    sys.path.append(os.path.expanduser('~leo/python/'))
else:
    sys.path.append('/data/private/soft/python/')
    os.environ["RSCRATCH"] = "/data/public/"
import cds_eua2 as eua
import pandas as pd
import xarray
import numpy
import hug
import h5py #ickle as h5py
import zipfile
import json
from multiprocessing import Pool
import glob
from functools import partial
from falcon import HTTPError, HTTP_400, HTTP_422
import copy
import time
from datetime import datetime, timedelta
import io
import matplotlib.pylab as plt
import subprocess

#@hug.exception(Exception)
#def handler(exception):
    #return str(exception)

def main():
    
    
    os.chdir(os.path.expandvars('$RSCRATCH/era5/odbs/merged'))
    z=numpy.zeros(1,dtype=numpy.int32)
    zidx=numpy.zeros(1,dtype=numpy.int)
    idx=numpy.zeros(1,dtype=numpy.int)
    trajectory_index=numpy.zeros(1,dtype=numpy.int)
    zz=eua.calc_trajindexfast(z,zidx,idx,trajectory_index)
    try:
        os.mkdir(os.path.expanduser('~/.tmp'))
    except:
        pass
    
    try:  
        with open(os.path.expanduser('~/.tmp/active.json')) as f:
            active=json.load(f)
    except:
        slist=glob.glob(os.path.expandvars('$RSCRATCH/era5/odbs/merged/'+'0-20000-0-?????_CEUAS_merged_v0.nc'))
        slnum=[i[-24:-19] for i in slist]
        
        volapath='https://oscar.wmo.int/oscar/vola/vola_legacy_report.txt'
        f=urllib.request.urlopen(volapath)
        col_names=pd.read_csv(f,delimiter='\t',quoting=3,nrows=0)
        f=urllib.request.urlopen(volapath)
        tdict={col: str for col in col_names}
        vola=pd.read_csv(f,delimiter='\t',quoting=3,dtype=tdict,na_filter=False)
 
        active={}
        
        for s,skey in zip(slist,slnum):
            try:
                
                with h5py.File(s,'r') as f:
                    print(skey)
                    try:
                        
                        #funits=f['recordtimestamp'].attrs['units']
                        funits='seconds since 1900-01-01 00:00:00'
                        active[skey]=[int(eua.secsince(f['recordtimestamp'][0],funits)),int(eua.secsince(f['recordtimestamp'][-1],funits)),
                                      float(f['observations_table']['latitude'][-1]),float(f['observations_table']['longitude'][-1])]
                        idx=numpy.where(vola.IndexNbr.values==skey)[0]
                        if len(idx)>0:
                            
                            active[skey].append(vola.CountryCode[idx[0]])
                        else:
                            active[skey].append('')
                            print('no key found for '+skey)
                    except KeyError:
                        print(skey+': a table is missing')
            except:
                print('file open error:'+s)
        with open(os.path.expanduser('~/.tmp/active.json'),'w') as f:
            json.dump(active,f)
            
        
    
    cf=eua.read_standardnames()
    
    #cdmtablelist=['id_scheme','crs','station_type','observed_variable','station_configuration_codes','units','sub_region']        
    cdmpath='https://raw.githubusercontent.com/glamod/common_data_model/master/tables/'                                                                                                                                                                                               
    cdmtablelist=['sub_region']        
    cdm=dict()
    for key in cdmtablelist:
        f=urllib.request.urlopen(cdmpath+key+'.dat')
        col_names=pd.read_csv(f,delimiter='\t',quoting=3,nrows=0)
        f=urllib.request.urlopen(cdmpath+key+'.dat')
        tdict={col: str for col in col_names}
        cdm[key]=pd.read_csv(f,delimiter='\t',quoting=3,dtype=tdict,na_filter=False)
        
    return active,cdm,cf

active,cdm,cf=main()

slnum=list(active.keys())
slist=[os.path.expandvars('$RSCRATCH/era5/odbs/merged/')+'0-20000-0-'+s+'_CEUAS_merged_v0.nc' for s in slnum]
#filelist=glob.glob('chera5.conv._?????.nc')

def check_body(body,cdm):
    required_keys=['variable']
    valid_keys=required_keys+['statid','fbstats','pressure_level','date','time','bbox','country']
    xor_keys=['statid','bbox','country']
    valid_ranges={}
    #valid_ranges['statid']=list(active.keys())# ['01001','99999']
    #if 'country' in body.keys():
        #valid_ranges['country']=['GLOBE','ALL']
        #for v in active.values():
            #if v[4] not in valid_ranges['country']:
                #valid_ranges['country'].append(v[4])
    #valid_ranges['bbox']=[[-180,180],[-90,90],[-180,180],[-90,90]]# luro
    valid_ranges['pressure_level']=['500','110000']
    valid_ranges['date']=['19000101','20301231']
    valid_ranges['time']=['0','24']
    valid_ranges['fbstats']=['obs_minus_an','obs_minus_bg','bias_estimate']
    valid_ranges['variable']=['temperature','u_component_of_wind','v_component_of_wind',
                              'wind_speed','wind_direction','relative_humidity',
                              'specific_humidity','dew_point_temperature'] #,'geopotential']
    
    try:
        
        bk=list(body.keys())
        for r in required_keys:
            if r not in bk:
                return ['Missing argument:', 'Argument '+r+' is required']
    
        for b in bk:
            if b not in valid_keys:
                return ['Invalid argument '+b+'.',' Valid arguments:'+str(valid_keys)]
            
        rxor=''
        for r in xor_keys:
            if r in bk:
                if len(rxor)==0:
                    
                    rxor=r
                else:
                    return ['Invalid station selection','Please do not specify both '+rxor+' and '+r]
        
        if 'country' in bk:
            if type(body['country']) is str:
                
                if body['country'].upper() in ('GLOBE','ALL'):
                    body['statid']=slnum
                    del body['country']
                else:
                    body['country']=[body['country']]
                    
            if 'country' in body.keys():
                vcountries=cdm['sub_region'].alpha_3_code.values
                body['statid']=[]
                for v in body['country']:
                    if v not in vcountries:
                        return ['Invalid station selection',v+' is not a valid country code']
                    for k,vv in active.items():
                        if vv[4]==v:
                            body['statid'].append(k)
                if len(body['statid'])==0:
                    return ['Invalid station selection','Countries '+str(body['country'])+' have no radiosondes']
                del body['country']
        elif 'bbox' in bk:
            if type(body['bbox']) is not list:
                return ['Invalid station selection','Bounding box: [lower, left, upper, right]']
            if len(body['bbox']) !=4:
                return ['Invalid station selection','Bounding box: [lower, left, upper, right]' ]
            try:
                for i in range(4):
                    body['bbox'][i]=float(body['bbox'][i])
            except:
                return ['Invalid station selection','Bounding box: [lower, left, upper, right] must be int or float']
            if body['bbox'][0]>body['bbox'][2] or body['bbox'][1]>body['bbox'][3]:
                return ['Invalid station selection','Bounding box requirements: lower<upper, left<right, -90<=lat<=90, -180<=lon<=360']
            if body['bbox'][0]<-90 or body['bbox'][0]>90 or body['bbox'][2]<-90 or body['bbox'][2]>90 or \
               body['bbox'][1]<-180 or body['bbox'][1]>360 or body['bbox'][3]<-180 or body['bbox'][3]>360 \
               or body['bbox'][3]-body['bbox'][1] >360:
                return ['Invalid station selection','Bounding box requirements: lower<upper, left<right, -90<=lat<=90, -180<=lon<=360']
            
            body['statid']=[]            
            for k,v in active.items():
                if v[2]>=body['bbox'][0] and v[2]<=body['bbox'][2]:
                    if body['bbox'][3]<=180:
                        if v[3]>=body['bbox'][1] and v[3]<=body['bbox'][3]:
                            body['statid'].append(k)
                    else: # rectangle crossing date line
                        if v[3]<0:
                            if v[3]>=body['bbox'][1]-360 and v[3]+360<=body['bbox'][3]:
                                body['statid'].append(k)
                        else:
                            if v[3]>=body['bbox'][1] and v[3]<=body['bbox'][3]:
                                body['statid'].append(k)
            if len(body['statid'])==0:
                return ['Invalid station selection','Bounding box '+str(body['bbox'])+' contains no radiosonde stations']
            del body['bbox']
        else:
            try:
                if type(body['statid'])==int:
                    body['statid']=['{:0>5}'.format(body['statid'])]
                elif body['statid']=='all':
                    body['statid']=slnum
                elif body['statid'][0]=='all':
                    body['statid']=slnum
                elif type(body['statid']) is not list:
                    body['statid']=[str(body['statid'])]
                else:
                    if type(body['statid'][0])==int:
                        for k in range(len(body['statid'])):
                            body['statid'][k]='{:0>5}'.format(body['statid'][k])
                    
                    
            except MemoryError:
                return ['Invalid station selection','Please specify either bbox, country or statid for station selection. Use "statid":"all" to select all stations']

        bk=list(body.keys())
        for v in bk:
            if v in valid_ranges:
                if type(body[v]) is not list:
                    body[v]=str(body[v])
                    if '-' in body[v]:
                        bvv=body[v].split('-')
                    else:
                        bvv=[body[v]]
                        body[v]=[body[v]]
                else:
                    for k in range(len(body[v])):
                        body[v][k]=str(body[v][k])
                    if len(body[v])==1 and '-' in body[v][0]:
                        bvv=body[v][0].split('-')
                    else:
                        bvv=body[v]
                    if v=='date':
                        d=int(bvv[0])
                        start=(datetime(year=d//10000,month=d%10000//100,day=d%100)-datetime(year=1900,month=1,day=1)).days
                        d=int(bvv[-1])
                        stop=(datetime(year=d//10000,month=d%10000//100,day=d%100)-datetime(year=1900,month=1,day=1)).days
                        if len(bvv)==stop-start+1:
                            body[v]=[bvv[0]+'-'+bvv[-1]]
                            bvv=[bvv[0],bvv[-1]]
                        
                for bv in bvv:
                    if v in ('pressure_level','date','time'):   
                        if int(bv)>int(valid_ranges[v][1]) or int(bv)<int(valid_ranges[v][0]):
                            return ['argument value(s) '+str(bvv)+' not valid.',
                                    'Valid values:'+str(valid_ranges[v])]
                    else:
                        if bv not in valid_ranges[v]:
                            return ['argument value(s) '+str(bv)+' not valid.',
                                    'Valid values:'+str(valid_ranges[v])]
                            
                print('body:',v,bvv[0],bvv[-1],len(body[v]))
                    #if len(bvv)>1:
                        #if body[v][0]==body[v][1]:
                            #body[v].pop()
                       
                #else:
                    #if type(body[v]) is not list:
                        #body[v]=[body[v],body[v]]
                    #for i in range(len(body[v])): #bv in body[v]:
                        #bv=body[v][i]
                        #print(bv)
                        #if int(bv) <valid_ranges[v][0] or int(bv)>valid_ranges[v][1]:
                            #return ['argument value(s) '+str(body[v])+' not valid.',
                                    #'Valid values:'+str(valid_ranges[v])]
                        #if v != 'statid':
                            #body[v][i]=int(body[v][i])
                    #print('body:',v,body[v][0],len(body[v]))
#                    if body[v][0]==body[v][1]:
#                        body[v].pop()
    except IOError:
        return ['general syntax error ',body]

    return ''

def makebodies(bodies,body,spv,bo,l):    

    for b in body[spv[l]]:

        if l<len(spv)-1:
            makebodies(bodies,body,spv,copy.copy(bo)+[b],l+1)
        else:
            bodies.append(dict(body))
            bn=copy.copy(bo)+[b]
            print('spv,bo:',spv,bn)
            for s,b in zip(spv,bn): 
                bodies[-1][s]=b
                print('makebodies',l,s,b)
    return        
    
def defproc(body,randdir,cdm):

    tt=time.time()
    error=check_body(body,cdm)
    print('body',body)
    if len(error)>0:
        raise HTTPError(HTTP_422,description=error)
        
    
    try:
        os.mkdir(randdir)
    except:
        pass

    bodies=[]
    spv=['statid','variable']
    bo=[]
    makebodies(bodies,body,spv,bo,0)         
    for k in range(len(bodies)-1,-1,-1):
        if 'date' in bodies[k].keys():
            #if bodies[k]['date'][0]>active[bodies[k]['statid']][1]//100 or bodies[k]['date'][-1]<active[bodies[k]['statid']][0]//100:
            if type(bodies[k]['date']) is not list:
                bodies[k]['date']=[bodies[k]['date']]
            dsec=[]
            dssold=''
            for ds in [bodies[k]['date'][0],bodies[k]['date'][-1]]:
                if '-' in ds:
                    if dssold=='':
                        dssold='-'
                        for dss in ds.split('-'):
                            d=int(dss)
                            dsec.append(((datetime(year=d//10000,month=d%10000//100,day=d%100)-datetime(year=1900,month=1,day=1))).days*86400)
                else:
                    d=int(ds)
                    dsec.append(((datetime(year=d//10000,month=d%10000//100,day=d%100)-datetime(year=1900,month=1,day=1))).days*86400)
            try:
                
                if dsec[0]>active[bodies[k]['statid']][1] or dsec[-1]+86399<active[bodies[k]['statid']][0]:
                    del bodies[k]
                    
            except KeyError:
                    del bodies[k]
        if 'time' in bodies[k].keys():
            if type(bodies[k]['time']) is not list:
                bodies[k]['time']=[bodies[k]['time']]
            tsec=[]
            tssold=''
            for ds in [bodies[k]['time'][0],bodies[k]['time'][-1]]:
                if '-' in ds:
                    if tssold=='':
                        tssold='-'
                        for dss in ds.split('-'):
                            d=int(dss)
                            tsec.append(d)
                else:
                    d=int(ds)
                    tsec.append(d)
            print('tsec:',tsec)
                
                                                                                                 
        
    
    if len(bodies)==0:
        raise HTTPError(HTTP_422,description=[body,'No selected station has data in specified date range'])
    print(bodies[:5])
    print('len:',len(bodies))
    
    #spvs={'statid':body.pop('statid'),'variable':body.pop('variable')}
    #for sv1 in body['statid']:,'variable']:
        #for sv1 in ['statid','variable']:
        
        #if type(body[splitvar]) is list:
            #for s in  body[splitvar]:
                #bodies.append(dict(body))
                #bodies[-1]['splitvar']=s
        #else:
            #if body[splitvar] == 'all':
                #statids=[]
                #for s in slist:
                    #bodies.append(dict(body))
                    #bodies[-1][splitvar]=s[-8:-3]
            #else:
                #bodies.append(dict(body))

    p=Pool(10)
    func=partial(eua.process_flat,randdir,cf)

    results=list(map(func,bodies))
    del p
    wpath=''
    for r in results:
        if r[0]!='':
            wpath=r[0]
            break
    if wpath=='':
        raise HTTPError(HTTP_422,description=[results[0][1],body])

    else:      
        rfile=os.path.dirname(wpath)+'/download.zip'

    print(results)
    print('wpath:'+wpath+';')
    with zipfile.ZipFile(rfile,'w') as f:
        for r in results:
            try:       
                if len(r[0])>0:
                    f.write(r[0],os.path.basename(r[0]))
                    
                os.remove(r[0])
            except:
                pass

    #for r in results:
        #os.remove(r[0])
        
    print('rfile:',rfile,'time: {:7.4f}s'.format(time.time()-tt))

    return rfile,''

def readandplot(rfile,body):

    with zipfile.ZipFile(rfile,'r') as a:
        for r in a.filelist:
            try:       
                with h5py.File(io.BytesIO(a.read(r)),'r') as hf:
                    print(r.filename,hf.keys())
                    qs='select date,time,vertco_reference_1,obsvalue,fg_depar@body,biascorr@body where date=20190101 and time>=40000 and time<50000 and varno=2'
                    #qs='select * where date=20190101 and time>=40000 and time<50000 and varno=2'
                    odbfile=os.path.expandvars('$RSCRATCH/era5/odbs/1/era5.conv.201901.10393')
                    rdata=subprocess.check_output(["odb","sql","-q",qs,"-i",odbfile,'--no_alignment'])
                    npos=rdata.index(b'\n')+1
                    xx=numpy.fromstring(rdata[npos:],dtype='float',sep='\t')
                    xx=xx.reshape((xx.shape[0]//6,6))
                    header=rdata[:npos].split()
                    check={}
                    for i in range(len(header)):
                        check[header[i]]=xx[:,i]
                    header=rdata[:npos].split()
                    datetime(1900,1,1)+timedelta(days=int(hf['time'][0]//86400),seconds=int(hf['time'][0]%86400))
                    if 'fbstats' in body.keys():
                        f,(ax1,ax2)=plt.subplots(1,2)
                        if type(body['fbstats']) is not list:
                            body['fbstats']=[body['fbstats']]
                        for fbs in body['fbstats']:
                            p=ax2.semilogy(hf[fbs],hf['plev'][:]/100,label=fbs)
                            
                        ax2.set_ylim(1030,5)
                        ax2.legend()
                        ax2.set_title(r.filename)
                    else:
                        f,ax1=plt.subplots(1,1)
                        
                    p=ax1.semilogy(hf['ta'],hf['plev'][:]/100,label='ta')
                    ax1.set_ylim(1030,5)
                    ax1.legend()
                    ax1.set_title(r.filename)
                    f.savefig(os.path.expanduser('~/plots/'+r.filename+'.png'))
                    plt.close()
                    
                    
            except MemoryError:
                pass
    return

def readandplot_ts(rfile,body):

    with zipfile.ZipFile(rfile,'r') as a:
        for r in a.filelist:
            try:       
                with h5py.File(io.BytesIO(a.read(r)),'r') as hf:
                    print(r.filename,hf.keys())
                    qs='select date,time,vertco_reference_1,obsvalue,fg_depar@body,biascorr@body where varno=2 and vertco_reference_1=10000'
                    #qs='select * where date=20190101 and time>=40000 and time<50000 and varno=2'
                    odbfile=os.path.expandvars('$RSCRATCH/era5/odbs/1/era5.conv.201901.10393')
                    rdata=subprocess.check_output(["odb","sql","-q",qs,"-i",odbfile,'--no_alignment'])
                    npos=rdata.index(b'\n')+1
                    rds=b'NaN'.join(rdata[npos:].split(b'NULL'))
                    xx=numpy.fromstring(rds,dtype='float',sep='\t')
                    xx=xx.reshape((xx.shape[0]//6,6))
                    header=rdata[:npos].split()
                    check={}
                    for i in range(len(header)):
                        check[header[i]]=xx[:,i]
                    header=rdata[:npos].split()
                    datetime(1900,1,1)+timedelta(days=int(hf['time'][0]//86400),seconds=int(hf['time'][0]%86400))
                    if 'fbstats' in body.keys():
                        f,(ax1,ax2)=plt.subplots(2,1)
                        if type(body['fbstats']) is not list:
                            body['fbstats']=[body['fbstats']]
                        for fbs in body['fbstats']:
                            mask=numpy.logical_and(hf['time'][:]%86400>=9*3600,hf['time'][:]%86400<15*3600)
                            p=ax2.plot(hf['time'][mask]/86400/365.25,hf[fbs][mask],label=fbs)
                            
                        #ax2.set_ylim(1030,5)
                        ax2.legend()
                        ax2.set_title(r.filename)
                    else:
                        f,ax1=plt.subplots(1,1)
                        
                    p=ax1.plot(hf['time'][mask]/86400/365.25,hf['dew_point_temperature'][mask],label='td')
                    #ax1.set_ylim(1030,5)
                    ax1.legend()
                    ax1.set_title(r.filename)
                    f.savefig(os.path.expanduser('~/plots/'+r.filename+'.png'))
                    plt.close()
                    
                    
            except MemoryError:
                pass
    return


@hug.get('/',output=hug.output_format.file)
def index(request=None,response=None):
    """
    index function requests get URI and converts into dictionary.
    Lists may be given via "[]"
    Allowed keys:
    statid (List) of strings
    bb= lat/lon rectangle (lower left upper right)
    variable (string) variable (one at a time)
    level (List) of levels in Pascal
    siglevs (bool) y/n
    glamod (bool)  y/n # add tables as hdf groups
    format (nc)
    """

    print(request.query_string)
    if '=' not in request.query_string:
        raise HTTPError(HTTP_422,description='A query string must be supplied')

    rs=request.query_string.split('&')
    body={}
    for r in rs:
        k,v=r.split('=')
        if '[' in v:
            vl=v[1:-1]
            if k in ['statid','variable','fbstats']:
                body[k]=vl.split(',')
            else:
                body[k]=list(numpy.fromstring(vl,sep=',',dtype='int'))
        else:
            body[k]=v

    randdir='{:012d}'.format(numpy.random.randint(100000000000))

    print(body)

    rfile,error=defproc(body,randdir,cdm)   

    if rfile=='':

        raise HTTPError(HTTP_422,description=error)

    response.set_header('Content-Disposition', 'attachment; filename='+os.path.basename(rfile))
    return rfile

@hug.post('/',output=hug.output_format.file)
def index(request=None,body=None,response=None):
    """
    index function requests get URI and converts into dictionary.
    Lists may be given via "[]"
    Allowed keys:
    statid (List) of strings
    bb= lat/lon rectangle (lower left upper right)
    variable (string) variable (one at a time)
    level (List) of levels in Pascal
    siglevs (bool) y/n
    glamod (bool)  y/n # add tables as hdf groups
    format (nc)
    """

    randdir='{:012d}'.format(numpy.random.randint(100000000000))

    print(json.dumps(body))

    rfile,error=defproc(body,randdir,cdm)   
    print(rfile,error)
    if rfile=='':
        raise HTTPError(HTTP_422,description=error)

    response.set_header('Content-Disposition', 'attachment; filename='+os.path.basename(rfile))

    return rfile


    
if __name__ == '__main__':

    #active,cdm,cf=main()
    body=eval(sys.argv[1])

    randdir=os.path.expandvars('$RSCRATCH/tmp/{:012d}'.format(100000000000))
    ret=defproc(body,randdir,cdm)
    print(ret)
    print('ready')
