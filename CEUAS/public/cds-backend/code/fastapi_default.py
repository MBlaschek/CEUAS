#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__version__ = '0.3'
__author__ = 'LH'
__status__ = 'dev'
__date__ = 'Mon Jul 13 10:29:21 CEST 2020'
__institute__ = 'UNIVIE'
__github__ = 'git@github.com:MBlaschek/CEUAS.git'
__doc__ = """
HUG Server v%s
Maintained by %s at %s
Github: %s [%s]
Updated: %s
""" % (__version__, __author__, __institute__, __github__, __status__, __date__)

"""
This is a simple controller
 - index is the default action of any application
   here it forwards the request from a CDS server to further processing
   with python

   For queries ask early-upper-air@copernicus-climate.eu
   The C3S 311c Lot2 group
   Vienna, 26 August 2019
 - user is required for authentication and authorization
 - download is for downloading files uploaded in the db (does streaming)
 - api is an example of Hypermedia API support and access control


License
This file is released under public domain and you can use without limitations

"""
import copy
import glob
import geopy
import csv
import json
import reverse_geocoder as rg
import pycountry
import logging
import os
import socket
import sys
import time
import zipfile
import urllib
import multiprocessing
from datetime import datetime, timedelta
from functools import partial
from multiprocessing import set_start_method, Pool
from typing import Union
from shutil import copyfile
import pickle
from itertools import product
import numpy as np
import h5py
import numpy
import pandas as pd
import xarray
import gzip

import asyncio
from fastapi import FastAPI
from fastapi.responses import FileResponse
from fastapi import Response
from fastapi import Request
app = FastAPI()
import subprocess

from pydantic import BaseModel
class Item(BaseModel):
    body: dict


#from typing import Any, Dict, AnyStr, List, Union
#JSONObject = Dict[AnyStr, Any]
#JSONArray = List[Any]
#JSONStructure = Union[JSONArray, JSONObject]



###############################################################################
#
# CONFIG
#
###############################################################################
config_file = 'hug.default.config.json'
global config
config = {'logger_name': 'upperair',
          'logger_level': 10,
          'logger_debug': 'hug.debug.log',
          'logger_dir': './logs',
          'src_path': '.',
          'data_dir': '.',
          'comp_dir': '.',
          'grid_dir': '.',
          'tmp_dir': '.',
          'config_dir': './config',
          'cds_eua_version': 4,
          'debug': False,
          'reload_pwd': 'reload'}

if os.path.isfile(config_file):
    new = json.load(open(config_file, 'r'))
    # expand ENVIRONMENTAL VARIABLES if any
    for ikey,ival in new.items():
        if 'dir' in ikey:
            new[ikey] = os.path.expandvars(ival)
            if '$' in new[ikey]:
                raise RuntimeError('configuration path expansion failed %s' %new[ikey])
    config.update(new)
else:
    print("Writing new config file:", config_file, "Adjust accordingly!")
    json.dump(config, open(config_file, 'w'))
    
if config['cds_eua_version'] == 3:
    import cds_eua3 as eua 
    # old variable codes
    CDS_EUA_VERSION = 3
elif config['cds_eua_version'] == 4:
    import cds_eua4 as eua
    CDS_EUA_VERSION = 4
elif config['cds_eua_version'] == 5:
    import cds_eua5 as eua
    CDS_EUA_VERSION = 5
    
print('CDS_EUA_VERSION: ', CDS_EUA_VERSION)
    
os.makedirs(config['logger_dir'], exist_ok=True)
os.makedirs(config['config_dir'], exist_ok=True)


###############################################################################
#
# LOGGING
#
###############################################################################
class MyFilter(object):
    def __init__(self, level):
        self.__level = level

    def filter(self, logrecord):
        return logrecord.levelno <= self.__level

request_logger = logging.getLogger('request_logger')
request_logger.setLevel(10)  # 10 Debug
r_ch = logging.FileHandler(config['logger_dir'] + '/request.log')
r_ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s | %(funcName)s - %(levelname)s - %(message)s')
r_ch.setFormatter(formatter)
request_logger.addHandler(r_ch)


logger = logging.getLogger(config['logger_name'])
logger.setLevel(config['logger_level'])  # 10 Debug
# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s | %(funcName)s - %(levelname)s - %(message)s')
# create console handler and set level to warning for stderr
if not config['debug']:
    ch = logging.StreamHandler()  # goes to std.err
    ch.setLevel(logging.ERROR)  # respond only to Debug and above
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    # create console handler and set level to info for stdout
    ch2 = logging.StreamHandler(sys.stdout)
    ch2.setLevel(logging.INFO)  # respond only to Debug and above
    ch2.setFormatter(formatter)
    logger.addHandler(ch2)
    ch2.addFilter(MyFilter(logging.INFO))

    try:
        ch3 = logging.FileHandler(config['logger_dir'] + '/' + config['logger_debug'])
    except PermissionError:
        ch3 = logging.FileHandler(config['logger_dir'] + '/' + config['logger_debug'].replace('.log', '.local.log'))
    ch3.setLevel(logging.DEBUG)  # respond only to Debug and above
    ch3.setFormatter(formatter)
    logger.addHandler(ch3)
    
else:
    # create console handler and set level to info for stdout
    ch2 = logging.StreamHandler(sys.stdout)
    ch2.setLevel(logging.DEBUG)  # respond only to Debug and above
    ch2.setFormatter(formatter)
    logger.addHandler(ch2)
    
for i, j in config.items():
    logger.debug('CONFIG %s : %s', i, j)

###############################################################################
#
# SPECIALS
#
###############################################################################
host = socket.gethostname()
logger.info("HUG started on %s", host)
# try:
#     # todo this part is deprecated / LEO delete?
#     if 'srvx' in host:
#         sys.path.append(os.path.expanduser('~leo/python/'))
#         config['data_dir'] = os.environ["RSCRATCH"]  # ?
# except:
#     pass

global constraints
try:
    constraints = '/tmp/constraints.csv'
    constraints = pd.read_csv(constraints)
    logger.info("constraints.csv read and ready")
except:
    logger.info("constraints.csv not found")

###############################################################################
#
# MAIN FUNCTIONS
#
###############################################################################

def makedaterange(vola: pd.DataFrame, itup: tuple, debug=False) -> dict:
    """ Read HDF5 radiosonde cdm backend file and extract datetime and geo information

    Args:
        vola: WMO Oscar Radiosonde list + Metadata
        itup: (filename, ID)

    Returns:
        dict : File Information:
            [ID] = [1.Date, last.Date, Lat, Lon, Country Code]
    """
    s, skey = itup  # filename, ID
    active = {}
    # creating a list for conversion between alpha_2 and alpha_3 countrycodes
    countrycodes = {}
    for country in pycountry.countries:
        countrycodes[country.alpha_2] = country.alpha_3
    countrycodes['XK']='XXK'
    try:

        with h5py.File(s, 'r') as f:
            logger.debug("SKEY %s", skey)
            try:

                # funits=f['recordtimestamp'].attrs['units']
                funits = 'seconds since 1900-01-01 00:00:00'
                # todo add list of available variables for future
                if CDS_EUA_VERSION == 2:
                    active[skey] = [int(eua.secsince(f['recordtimestamp'][0], funits)),
                                    int(eua.secsince(f['recordtimestamp'][-1], funits)),
                                    float(f['observations_table']['latitude'][-1]),
                                    float(f['observations_table']['longitude'][-1])]
                if CDS_EUA_VERSION >= 3:
                    #if isinstance(f['recordindices'], h5py._hl.group.Group):
                    if 'recordindices' in f.keys():
                        active[skey] = [int(eua.to_seconds_since(f['recordindices']['recordtimestamp'][0], funits)),
                                    int(eua.to_seconds_since(f['recordindices']['recordtimestamp'][-1], funits)),
                                    float(f['observations_table']['latitude'][-1]),
                                    float(f['observations_table']['longitude'][-1])]
                    else:
                        try:
                            
                            active[skey] = [int(eua.to_seconds_since(f['recordtimestamp'][0], funits)),
                                        int(eua.to_seconds_since(f['recordtimestamp'][-1], funits)),
                                        float(f['observations_table']['latitude'][-1]),
                                        float(f['observations_table']['longitude'][-1])]
                        except:
                            print('makedaterange: ',s.split('/')[-1],' recordtimestamp does not exist')
                            return active
                idx = numpy.where(vola.StationId.values == skey)[0]
                if len(idx) > 0:
                    active[skey].append(vola.CountryCode[idx[0]])
                else:
                    # if no country code available -> reverse geo search for them 
                    coordinates = (float(f['observations_table']['latitude'][-1]), float(f['observations_table']['longitude'][-1]))

                    #cc = rg.search(coordinates)[0]['cc']
                    cc='XK'
#                     if cc == 'XK':
#                         active[skey].append('XXK')
#                         logger.debug('reverse geo searche for: %s', skey)
#                     else:
                    # results are in alpha_2 country codes -> convert to alpha_3 like it is in the vola file
                    active[skey].append(countrycodes[cc])
                    logger.debug('reverse geo searche for: %s', skey)
            # add data directory for process_flat
                # active[skey].append(os.path.dirname(s))
                # add filepath
                active[skey].append(s)
            except KeyError as e:
#                 if debug:
#                     raise e
                logger.debug('%s : a table is missing', skey)
    except:
        if debug:
            raise 
        logger.error('file open error: %s', s)
    return active

def read_tstamps(fn):
    if '0-20100-0-01802' in fn:
        print(fn)
    with h5py.File(fn,'r') as f:
    #print(f[fk]['recordindices'].keys())
        fk=fn.split('/')[-1].split('_CEUAS_merged')[0]
        print(fk)
        rts=f['recordindices']['recordtimestamp'][:]
    return fk,rts

def pkl_initialize(config,slist=[]):

    #set_start_method('forkserver') 

##    flist=glob.glob(os.path.expandvars(config['data_dir'] + 'converted_v5/0-*-0-*_CEUAS_merged_v1.nc')) 
    #config['data_dir']='/raid60/scratch/leo/scratch/converted_v5'
    #config['comp_dir']='/raid60/scratch/leo/scratch/converted_v5'
    #slist = glob.glob(os.path.expandvars(config['data_dir'] + '/0-2000?-0-?????_CEUAS_merged_v1.nc'))
    #slist += glob.glob(os.path.expandvars(config['data_dir'] + '/0-20?00-0-*_CEUAS_merged_v1.nc'))
    #slist += glob.glob(os.path.expandvars(config['comp_dir'] + '/0-20?00-0-?????.nc'))
    #slist += glob.glob(os.path.expandvars(config['comp_dir'] + '/0-20?00-0-?????_CEUAS_merged_v0.nc'))
    #flist=glob.glob(os.path.expandvars(rpath+'*.nc'))
    fout=os.path.expandvars(config['config_dir']+'/h5link.pkl')
    
    tt=time.time()
    rtsdict={}
    i=0
    imem=0
    
    if not slist:
        
        try:
                
            with open(fout, 'rb') as f:
                rtskeys,rtsidx,rtsarr=pickle.load(f)
        except: 
            raise ValueError('cannot read '+fout)
    else:         
        #with h5py.File(fout,'r') as f:
        l=0
        #with Pool(10) as p:
        tup=map(read_tstamps,list(slist))
        rtsdict=dict(tup)
        
#         rtsdict = {}
#         for i in slist:
#             a,b = read_tstamps(i)
#             rtsdict[a]=b

        #for fn in flist:
            #with h5py.File(fn,'r') as f:
            ##print(f[fk]['recordindices'].keys())
                #fk=fn.split('/')[-1].split('_CEUAS_merged_v1.nc')[0]
                #rtsdict[fk]=f['recordindices']['recordtimestamp'][:]
                #print(fk,l)
                #l+=1
                #flinks[fk]=f[fk]
                #imem+=sys.getsizeof(rtsdict[fk])
                #print(i,imem,fk)
                #i+=1
            #with open('.pkl'.join(fout.split('.nc')), 'wb') as f:
                #pickle.dump(rtsdict, f, pickle.HIGHEST_PROTOCOL)
            
        rtsarr=numpy.concatenate(list(rtsdict.values()))      ## contains the record timestamps
        rtsidx=[0]+[len(rtsdict[v]) for v in rtsdict.keys()]  ## contains the indices where the timestamps of station i starts
        rtsidx=numpy.cumsum(rtsidx)
        rtskeys=list(rtsdict.keys())                          ## contains the station IDs
        with open(fout, 'wb') as f:
            pickle.dump((rtskeys,rtsidx,rtsarr), f, pickle.HIGHEST_PROTOCOL)
        
        x=0
    
    
    print('ready read',time.time()-tt) ; tt=time.time()
    return rtskeys,rtsidx,rtsarr,fout
    

def init_server(force_reload: bool = False, force_download: bool = False, debug:bool = False) -> tuple:
    """ Initialize Radiosonde Archive and read CDM Informations and CF Convention

    https://raw.githubusercontent.com/glamod/common_data_model/master/tables/

    Returns:
        dict : active stations
        dict : Common Data Model Name and Attributes
        dict : Climate Convention Naming
    """
    import os
    import json
    import urllib.request

    # os.chdir(os.path.expandvars('$RSCRATCH/era5/odbs/merged'))
    # wroot = os.path.expandvars('$RSCRATCH/era5/odbs/merged/tmp')
    # os.makedirs(wroot, exist_ok=True)
    z = numpy.zeros(1, dtype=numpy.int32)
    zidx = numpy.zeros(1, dtype=numpy.int)
    idx = numpy.zeros(1, dtype=numpy.int)
    trajectory_index = numpy.zeros(1, dtype=numpy.int)
    # What get's updated here?
    zz = eua.calc_trajindexfast(z, zidx, idx, trajectory_index)
    os.makedirs(config['config_dir'], exist_ok=True)
    #
    # Active Json:
    # [WIGOS ID] = [start time in seconds, end time in seconds, lat, lon, datadir]
    # seconds since 1900-01-01
    #
    active_file = config['config_dir'] + '/active.json'
    namelist_file  = config['config_dir'] + '/namelist.json'
    
#     config['data_dir']='/raid60/scratch/leo/scratch/converted_v5'
#     config['comp_dir']='/raid60/scratch/leo/scratch/comp'
    
    namelist = None
    active = None
    if os.path.isfile(active_file) and not force_reload:
        try:
            with open(active_file) as f:
                active = json.load(f)
            logger.info('Active Stations read. [%d]', len(active))
            
            rtskeys,rtsidx,rtsarr,fout=pkl_initialize(config)
            active['rtsarr']=rtsarr
            active['rtsidx']=rtsidx
            active['rtskeys']=rtskeys
        except Exception as e:
            logger.info('Active Stations will be created.')
            active = None
    if os.path.isfile(namelist_file) and not force_reload:
        try:
            with open(namelist_file) as f:
                namelist = json.load(f)
            logger.info('Namelist read. [%d]', len(namelist))
        except:
            namelist = None

    if active is None:
        #
        # find Merged Netcdf files and intercomparison files
        #
#        slist = glob.glob(os.path.expandvars(config['data_dir'] + '/0-2000?-0-?????_CEUAS_merged_v1.nc'))
#        slist += glob.glob(os.path.expandvars(config['data_dir'] + '/0-20?00-0-*_CEUAS_merged_v1.nc'))
        slist = glob.glob(os.path.expandvars(config['data_dir'] + '/0-*_CEUAS_merged_v1.nc'))
        slist += glob.glob(os.path.expandvars(config['comp_dir'] + '/0-20?00-0-?????.nc'))
        slist += glob.glob(os.path.expandvars(config['comp_dir'] + '/0-20?00-0-?????_CEUAS_merged_v0.nc'))
        # slnum = [i[-34:-19] for i in slist]
        slnum = [i.split('/')[-1].split('_')[0].replace('.nc','') for i in slist]
        volapath = 'https://oscar.wmo.int/oscar/vola/vola_legacy_report.txt'
        f = urllib.request.urlopen(volapath)
        col_names = pd.read_csv(f, delimiter='\t', quoting=3, nrows=0)
        # print(col_names)
        f = urllib.request.urlopen(volapath)
        tdict = {col: str for col in col_names}
        f = urllib.request.urlopen(volapath)
        vola = pd.read_csv(f, delimiter='\t', quoting=3, dtype=tdict, na_filter=False)
        # print (vola.iloc[0])
        # exit()
        active = {}
        func = partial(makedaterange, vola, debug=debug)
        #slist=[slist[3380]]
        #slnum=[slnum[3380]]
        if False:
#             with Pool(10) as p:
#                 sklist=list(p.map(func,zip(slist,slnum)))
#             with multiprocessing.get_context('spawn').Pool(10) as p:
            with Pool(9) as p:
                sklist=list(p.map(func,zip(slist,slnum)))
        else:
            sklist = list(map(func, zip(slist, slnum)))
        
        short_slist=[]
        for s,sl in zip(sklist,slist):
            if s:
                k = next(iter(s))
                if k in active.keys():
                    if active[k]!=s[k]:
                        print(active[k],'\n',s[k])
                else:
                    active[k] = s[k]
                    short_slist.append(active[k][-1])
        logger.info('Active Stations created. [%d]', len(active))
        
        try:
            with open(active_file, 'w') as f:
                json.dump(active, f)
        except Exception as e:
            logger.warning('Cannot write %s: %s', active_file, e)

        try:
            
            rtskeys,rtsidx,rtsarr,fout=pkl_initialize(config,slist=short_slist) # slist
            active['rtsarr']=rtsarr
            active['rtsidx']=rtsidx
            active['rtskeys']=rtskeys
            
        except Exception as e:
            logger.warning('Cannot write %s: %s', 'timestamp file h5link.pkl', e)
        
    
            
    if namelist is None:
        volapath = 'https://oscar.wmo.int/oscar/vola/vola_legacy_report.txt'
        f = urllib.request.urlopen(volapath)
        col_names = pd.read_csv(f, delimiter='\t', quoting=3, nrows=0)
        f = urllib.request.urlopen(volapath)
        tdict = {col: str for col in col_names}
        vola = pd.read_csv(f, delimiter='\t', quoting=3, dtype=tdict, na_filter=False)
        active_file = config['config_dir'] + '/active.json'
        act = json.load(open(active_file,"r"))

        namelist = {}
        for i in act:
            try:
                name = vola[vola['StationId']==i]['StationName'].iloc[0]
                name = name[:name.find(',')]
            except:
                if (i[:5] == '0-200'):
                    name = i
                elif (i[:5] == '0-201') or (i[:5] == '0-202'):
                    name = 'intercomparison campagne'
                elif (i[:5] == '0-203'):
                    name = 'xxxxxx'
                else:
                    name = 'missing station name'
            namelist[i]=name
        try:
            with open(namelist_file, 'w') as f:
                json.dump(namelist, f)
        except Exception as e:
            logger.warning('Cannot write %s: %s', namelist_file, e)
    #
    # Read CDM Definitions
    #
    cdm_file = config['config_dir'] + '/cf.json'
    if os.path.isfile(cdm_file) and not force_download:
        with open(cdm_file) as f:
            cf = json.load(f)
    else:
        cf = eua.read_standardnames()
        try:
            with open(cdm_file, 'w') as f:
                json.dump(cf, f)
        except Exception as e:
            logger.warning('Cannot write %s: %s', cdm_file, e)
    #
    # list of country codes -> used for country selection in check_body
    #
    cdmpath = 'https://raw.githubusercontent.com/glamod/common_data_model/master/tables/'
    cdmtablelist = ['sub_region']
    cdm = dict()
    for key in cdmtablelist:
        f = urllib.request.urlopen(cdmpath + key + '.dat')
        col_names = pd.read_csv(f, delimiter='\t', quoting=3, nrows=0)
        f = urllib.request.urlopen(cdmpath + key + '.dat')
        tdict = {col: str for col in col_names}
        cdm[key] = pd.read_csv(f, delimiter='\t', quoting=3, dtype=tdict, na_filter=False)

    
    return active, cdm, cf


###############################################################################
#
# Common Data Objects
#
# Active Station Dictionary
# WMO Regions Table for Country Codes
# CF Table Naming
###############################################################################


active, wmo_regions, cf = init_server()

# Active Station Numbers
slnum = list(active.keys())


# slist = [config['data_dir'] + '/0-20000-0-' + s + '_CEUAS_merged_v0.nc' for s in slnum]
# slist = [s[5] for _,s in active.items()]

# try:
# set_start_method("spawn")  # or fork ? not sure why, pickling?
# set_start_method("forkserver")  # fork is not threadsafe, unfortunately
# with multiprocessing.get_context('spawn').Pool(16) as P:
with Pool(9) as P:
    x=P.map(np.sin,np.arange(16))

# P=Pool(16) 
# x=P.map(np.sin,np.arange(16))
print(x)
# except RuntimeError:
#     pass

###############################################################################
#
# Functions
#
###############################################################################


def last_day_of_month(any_day: datetime) -> datetime:
    """ Get the last day of a month

    Args:
        any_day: Datetime

    Returns:
        datetime : last day of the month
    """
    next_month = any_day.replace(day=28) + timedelta(days=4)  # this will never fail
    return next_month - timedelta(days=next_month.day)


def to_valid_datetime(idate: str, as_string: bool = False) -> Union[str, datetime]:
    """ Convert date string or integer to a valid datetime (last day of the month)

    Args:
        idate: datetime like: 19990131
        as_string: return a string or the datetime object

    Returns:
        str : '19990131'
        datetime : datetime(1999,1,31)
    """
    d = 0
    try:
        d = int(idate)
        idate = datetime(year=d // 10000, month=d % 10000 // 100, day=d % 100)

    except ValueError:
        if (d % 100) > 28:
            idate = last_day_of_month(datetime(year=d // 10000, month=d % 10000 // 100, day=1))
            logger.debug('Date changed from %d to %s', d, idate.strftime('%Y%m%d'))
    if as_string:
        return idate.strftime('%Y%m%d')
    return idate


@app.get('/status')
def status_test(command=None) -> dict:
    """ Return the status of the hug server

    Returns:

    """
    import psutil
    hproc = None
    for proc in psutil.process_iter():
        if 'bin/hug' in " ".join(proc.cmdline()):
            hproc = proc
            break

    if hproc is not None:
        elapsed = datetime.now() - datetime.fromtimestamp(hproc.create_time())
        status_msg = {"version": __version__, "status": hproc.status(), "running": hproc.is_running(),
                      "available": str(elapsed), "memory": hproc.memory_percent(), "cpu": hproc.cpu_percent(),
                      "num_stations": len(slnum), "active": active}
        # psutil.disk_usage('/tmp/')  # '/data/private/',
        if command == config['reload_pwd']:
            if elapsed.total_seconds() > 120:
                # todo run this in background and wait until a request is finished before restarting
                init_server(force_reload=True, force_download=False)
                hproc.kill()  # this should end the hug server and cron should restart it
            else:
                status_msg['command'] = 'Sorry restarted only %d seconds ago [120 s]' % elapsed.total_seconds()
            return status_msg
        if command == 'restart':
            # raise RuntimeError("Restart requested ... killing myself softly.")
            hproc.kill()
            
        if command == 'cleanup':
            # search for request directories and results / remove them and report back
            pass

        if command == 'failed_requests':
            messages = {}
            with open(config['logger_dir'] + '/failed_requests.log', 'r') as f:
                tmp = f.read().splitlines()
                for iline in tmp:
                    idate = iline.split(' [')[0]
                    iline = iline.replace(idate, '').split('Message:')
                    messages[idate] = {'request': eval(iline[0].strip()[1:-1]), 'message': iline[1].strip()}
            status_msg["failed_requests"] = messages

        if command == 'finished_requests':
            messages = {}
            with open(config['logger_dir'] + '/finished_requests.log', 'r') as f:
                tmp = f.read().splitlines()
                for iline in tmp:
                    idate = iline.split(' [')[0]
                    iline = iline.replace(idate, '')
                    messages[idate] = {'request': eval(iline.strip()[1:-1])}
            status_msg["finsihed_requests"] = messages

        if command == 'running':
            return {'running': hproc.is_running()}

        return status_msg

    return {'error': "Can't find me....? :("}


def to_csv(flist: list, ofile: str = 'out.csv', name: str = 'variable'):
    """ Convert every file in flist to CSV

    Args:
        flist: list of files
        ofile: output filename
        name: variable filename inside zip file

    Returns:
        str: output filename (returned by the request)

    Profiling:
        on    Tue Jul 28 09:12:40 UTC 2020, 1143 temperature files, 1979-01, stdplevs

        Total time: 56.3918 s (load_dataset)
        Total time: 41.8788 s (zipping)
        Total time: 39.8347 s (open_dataset)
        File: /data/private/soft/python/hug/default.py
        Function: to_csv at line 366

        Line #      Hits         Time  Per Hit   % Time  Line Contents
        ==============================================================
           366                                           def to_csv(flist: list, ofile: str = 'out.csv'):
           377         1          2.0      2.0      0.0      statindex = 0
           378         1          1.0      1.0      0.0      dfs = []
           379      1144       1654.0      1.4      0.0      for fn in flist:
           380                                                   # changed from open_dataset to load_dataset
           381      1143   14480897.0  12669.2     36.4          ds = xarray.open_dataset(fn, drop_variables=['trajectory_label', 'trajectory_index', 'trajectory'])
           382      1143    7763230.0   6792.0     19.5          df = ds.to_dataframe()
           383      1143       9120.0      8.0      0.0          if 'primary_id' not in ds.attrs:
           384                                                       # /tmp//tmp//006691463272/dest_0-20000-0-53513_relative_humidity.nc
           385      1143    1462572.0   1279.6      3.7              df['statid'] = fn.split('/')[-1].split('_')[1]
           386                                                       # logger.warning('CSV no primary_id in %s', fn)
           387                                                       # continue
           388                                                   else:
           389                                                       df['statid'] = ds.attrs['primary_id']
           390      1143    1332445.0   1165.7      3.3          df['statindex'] = statindex
           391      1143       3120.0      2.7      0.0          dfs.append(df)
           392      1143       1151.0      1.0      0.0          statindex += 1
           393
           394         1     646098.0 646098.0      1.6      df = pd.concat(dfs, ignore_index=True)
           395         1         16.0     16.0      0.0      df.index.name = 'obs_id'
           396         1   14134345.0 14134345.0   35.5      df.to_csv(ofile)
           397         1          2.0      2.0      0.0      return ofile

    """
    statindex = 0
    dfs = []
    for fn in flist:
        # open_dataset (~20 s faster) than load_dataset
        # ds = xarray.open_dataset(fn, drop_variables=['trajectory_label', 'trajectory_index', 'trajectory'])
        logger.debug('Converting %s', fn)
        with eua.CDMDataset(fn) as ds:
            df = ds.to_dataframe()
        
#         to_be_removed = ['trajectory_index', 'trajectory']
#         for ivar in list(ds.variables):
#             if 'string' in ivar:
#                 to_be_removed.append(ivar)
            
#             if 'trajectory' in ds[ivar].dims and ivar not in list(ds.coords):
#                 report_id = ds[ivar].astype(object).sum(axis=1).astype(str)
#                 ds = ds.drop_vars(ivar)
#                 ds[ivar] = ('obs', report_id.values[ds.trajectory_index.values])  # todo obs ???
                
#             if ds[ivar].ndim > 1:
#                 tmp= ds[ivar].astype(object).sum(axis=1).astype(str)
#                 ds = ds.drop_vars(ivar)
#                 idim = tmp.dims[0]
#                 ds[ivar] = (idim, tmp)
        
    #         ds = ds.drop_vars(to_be_removed, errors='ignore')  # do not raise an error.
        
#         print(df)
#         timevar = df['time']
#         df.drop('time', axis=1)
#         df['time'] = timevar

        #
        # todo fix the primary_id in the NetCDF files
        #
#         if 'primary_id' not in ds.attrs:
#             # /tmp//tmp//006691463272/dest_0-20000-0-53513_relative_humidity.nc
#             df['statid'] = fn.split('/')[-1].split('_')[1]
#             # logger.warning('CSV no primary_id in %s', fn)
#             # continue
#         else:
#             df['statid'] = ds.attrs['primary_id']
        #
        #
        #
        
        
        att_ds = xarray.open_dataset(fn)   
        if 'ta' in df.keys():
            df['ta_units'] = att_ds['ta'].attrs['units']
        if 'va' in df.keys():
            df['va_units'] = att_ds['va'].attrs['units']
        if 'ua' in df.keys():
            df['ua_units'] = att_ds['ua'].attrs['units']
        if 'wind_speed' in df.keys():
            df['wind_speed_units'] = att_ds['wind_speed'].attrs['units']
        if 'wind_from_direction' in df.keys():
            df['wind_from_direction_units'] = att_ds['wind_from_direction'].attrs['units']
        if 'hur' in df.keys():
            df['hur_units'] = att_ds['hur'].attrs['units']
        if 'hus' in df.keys():
            df['hus_units'] = att_ds['hus'].attrs['units']
        if 'dew_point_temperature' in df.keys():
            df['dew_point_temperature_units'] = att_ds['dew_point_temperature'].attrs['units']
        if 'geopotential' in df.keys():
            df['geopotential_units'] = att_ds['geopotential'].attrs['units']
        df['plev_units'] = att_ds['plev'].attrs['units']
        df['statindex'] = statindex
        statindex += 1
        dfs.append(df)
        
    df = pd.concat(dfs, ignore_index=True)
    df.index.name = 'obs_id'
    #if '.zip' in ofile:
    #    df.to_csv(ofile, compression=dict(method='zip', archive_name=name + '.csv'), mode='a')  # might be 10x faster
    #else:
    df.to_csv(ofile)
    return ofile


###############################################################################


def check_body(observed_variable: list = None, variable: list = None, statid: list = None, product_type: str = None, 
               pressure_level: list = None,
               day: list = None, month: list = None, year: list = None, date: list = None, time: list = None, 
               bbox: list = None, country: str = None, area: list = None,
               format: str = None, period: list = None, optional: list = None, wmotable: dict = None,
               gridded: list = None, toolbox: str = None, cdm: list = None, da: bool = True, compression: str = None,
               pass_unknown_keys: bool = False, nodims: str = None, hdf: str = None, speed_test: str = None, 
               single_parallel: bool = False, single_csv: bool = True, single_csv_target: str = 'CDS_CUON_output_file',
               **kwargs) -> dict:
    """ Check Request for valid values and keys

    Args:
        variable: e.g. temperature, ...
        statid: e.g. 01001
        product_type: currently unused
        pressure_level: Pa, 500 - 110000 Pa values allowed
        date: '19990131' or ['19990101', '20000101'] as range
        time: 1,... or 0-23
        bbox: Bounding Box [lower left upper right]
        country: Country Code, e.g. DEU, AUT, USA, GBR'
        format: nc or csv
        period: ['19990101', '20000101'] see Notes
        wmotable: WMO Regions definition Table
        pass_unknown_keys: only for debugging and local use
        **kwargs:

    Returns:
        dict : Clean Request

    Notes:
        date and period are currently not really separated due to CDS filtering
    """
    d = {}
    '''
            FORMS:
            'air_temperature', 'dew_point_depression', 'dew_point_temperature',
            'eastward_wind_speed', 'geopotential', 'northward_wind_speed',
            'relative_humidity', 'specific_humidity', 'wind_from_direction',
            'wind_speed
            
            PUG:
            ["air_temperature", "zonal_wind", "meridional_wind",
            "wind_speed", "wind_direction", "air_relative_humidity", 
            "air_specific_humidity", "air_dewpoint"]
            
            PREVIEW:
            'air_temperature',
            'dew_point_temperature',
            'eastward_wind_speed',
            'northward_wind_speed',
            'relative_humidity',
            'specific_humidity',
            'wind_from_direction',
            'wind_speed'
            
    '''
    allowed_variables = ['temperature', 'air_temperature',
                         'u_component_of_wind', 'eastward_wind_speed',
                         'v_component_of_wind', 'northward_wind_speed',
                         'wind_speed',
                         'wind_direction', 'wind_from_direction',
                         'relative_humidity',
                         'specific_humidity',
                         'dew_point_temperature',
                         'geopotential', 'geopotential_height',
                         'dew_point_depression', 'dewpoint_departure','dewpoint_depression', 'dew_point_departure']
    #
    # single csv target name
    #
    d['single_csv_target']=single_csv_target
    
    #
    # Unknown keys ?
    #
    for ikey, ival in kwargs.items():
        logger.info('Requested key %s : %s unknown', ikey, ival)
        if pass_unknown_keys:
            d[ikey] = ival
    #
    # Product Type
    #
    # todo add product_type as a variable
    # possible values: [sounding, monthly, gridded]
    if product_type is not None:
        logger.warning('Not yet implemented: product_type : %s' % product_type)
    
    #
    # Direct Access
    #
    if da is not None:
        if da == 'False':
            da = False
        if da == 'True':
            da = True
        if not isinstance(da, bool):
            raise KeyError("Invalid type selected at da - only bool is valid: " + da)
        else:
            d['da'] = da

    #
    # Compression
    #
    if compression is not None:
        if compression in ['gzip','lzf']:
            d['compression'] = compression
    else:
        d['compression'] = 'gzip'
            
    
    #
    # Variable
    #
    if variable is None:
        if observed_variable is None:
            raise KeyError('Missing argument: variable')
        else:
            variable = observed_variable

    if not isinstance(variable, list):
        variable = [variable]

    for ivar in variable:
        if ivar not in allowed_variables:
            raise KeyError('Invalid variable selected: ' + ivar)
    for i in range(len(variable)):
        if variable[i] in ['dewpoint_departure', 'dew_point_departure','dewpoint_depression']:
            variable[i] = 'dew_point_depression'
        if variable[i] in ['geopotential_height']:
            variable[i] = 'geopotential'
        if variable[i] in ['air_temperature']:
            variable[i] = 'temperature'
        if variable[i] in ['eastward_wind_speed']:
            variable[i] = 'u_component_of_wind'
        if variable[i] in ['northward_wind_speed']:
            variable[i] = 'v_component_of_wind'
        if variable[i] in ['wind_from_direction']:
            variable[i] = 'wind_direction'
                
        d['variable'] = variable
    #
    # Optional
    #
    allowed_optionals = ['sonde_type', 'bias_estimate','obs_minus_an','obs_minus_bg', 'bias_estimate_method', 
                         'RISE_bias_estimate', 'RICH_bias_estimate', 'RASE_bias_estimate', 'RAOBCORE_bias_estimate',
                         'latitude_displacement', 'longitude_displacement', 'time_since_launch', 'true_time', 'station_elevation',
                         'RISE_1.8_bias_estimate', 'RICH_1.8_bias_estimate', 'RASE_1.8_bias_estimate', 'RAOBCORE_1.8_bias_estimate',
                         'desroziers_30', 'desroziers_60', 'desroziers_90', 'desroziers_180',
                         'wind_bias_estimate',
                         'humidity_bias_estimate', 'humidity_1.0_bias_estimate',
                        ]
    # bias_estimate_method : raobcore, rich, ...
    if optional is not None:
        if not isinstance(optional, list):
            if optional in allowed_optionals:
                d['optional'] = [optional]
            else:
                raise KeyError('Invalid optional selected: ' + optional)
        else:
            for iopt in optional:
                if iopt not in allowed_optionals:
                    raise KeyError('Invalid optional selected: ' + optional)
            d['optional'] = optional
#         new_opts = []
#         for opts in d['optional']:
#             if opts == 'RISE_bias_estimate':
#                 new_opts.append('RISE_1.8_bias_estimate')
#             elif opts == 'RICH_bias_estimate':
#                 new_opts.append('RICH_1.8_bias_estimate')
#             elif opts == 'RASE_bias_estimate':
#                 new_opts.append('RASE_1.8_bias_estimate')
#             elif opts == 'RAOBCORE_bias_estimate':
#                 new_opts.append('RAOBCORE_1.8_bias_estimate')
#             elif opts == 'humidity_bias_estimate':
#                 new_opts.append('humidity_1.0_bias_estimate')
#             else:
#                 new_opts.append(opts)
#         d['optional'] = new_opts

    #
    # toolbox
    #
    if toolbox is not None:
        if not toolbox == 'True':
            raise KeyError("Invalid type selected at toolbox - only string 'True' is valid: " + toolbox)
        try: 
            d['optional']
        except:
            pass
        else:
            if len(d['optional']) > 1:
                raise KeyError("toolbox 'True' is only valid if only one 'optional' is given.")  
        d['toolbox'] = True
        
    #
    # nodims
    #
    if nodims is not None:
        if not nodims == 'True':
            raise KeyError("Invalid type selected at nodims - only string 'True' is valid: " + nodims)
        try: 
            d['optional']
        except:
            pass
        else:
            if len(d['optional']) > 1:
                raise KeyError("nodims 'True' is only valid if only one 'optional' is given.")  
        d['nodims'] = True
    
    #
    # CDM
    #
    if cdm is not None:
        if not isinstance(cdm, list):
            if isinstance(cdm, str):
                d['cdm'] = [cdm]
            else:
                raise KeyError("Invalid type selected at CDM - only string is valid: " + cdm)
        else:
            d['cdm'] = cdm
            d['da'] = False
        print(d['cdm'])
#         if not ('observations_table/data_policy_licence' in d['cdm']):
#             d['cdm'].append('observations_table/data_policy_licence')
#             print(d['cdm'])
    else:
#         d['cdm'] = ['observations_table/data_policy_licence']
        d['cdm'] = None

    #
    # gridded [lower left upper right]
    #
    if gridded is not None:
        if not isinstance(gridded, (list, tuple)) or len(gridded) != 4:
            raise ValueError('Invalid selection, gridded: [lower left upper right]')
        try:
            for i in range(4):
                gridded[i] = float(gridded[i])
        except ValueError:
            raise ValueError('Invalid selection, gridded: [lower left upper right] must be int or float')

        if gridded[0] >= gridded[2] or gridded[1] >= gridded[3]:
            raise ValueError('Invalid selection, gridded: lower<upper [-90, 90], left<right [-180, 360]')

        if gridded[0] < -90 or gridded[0] > 90 or gridded[2] < -90 or gridded[2] > 90 or \
                gridded[1] < -180 or gridded[1] > 360 or gridded[3] < -180 or gridded[3] > 360 \
                or gridded[3] - gridded[1] > 360:
            raise ValueError('Invalid selection, gridded: lower<upper [-90, 90], left<right [-180, 360]')
        d['gridded'] = gridded
            
    #
    # Format
    #
    if format in ['fast_csv']:
        d['format'] = 'fast_csv'
    elif format is not None:
        if format not in ['nc', 'csv']:
            raise ValueError('Invalid format selected [nc, csv]: ' + format)
        d['format'] = format
    else:
        d['format'] = 'nc'
        
    #
    # hdf
    #
    if hdf is not None:
        if hdf == 'True':
            d['hdf'] = True
        else:
            d['hdf'] = False
    else:
        d['hdf'] = False
        
        
    #
    # speed_test
    #
    if speed_test is not None:
        d['speed_test'] = True
        
    #
    # only one of [statid, bbox, country]
    #
    if statid is None and bbox is None and country is None:
        statid = 'all'

    if sum((statid is not None, bbox is not None, country is not None)) != 1:
        raise RuntimeError('Invalid selection, Specify only one of statid: %s, bbox: %s and country: %s' % (
            statid, bbox, country))
        
    #
    # area to bbox
    #
    if area is not None:
        bbox = area
    else:
        d['bbox']=False
    
    #
    # Countries
    #
    if country is not None:
        statid = []
        if isinstance(country, str):
            if country.upper() in ('GLOBE', 'ALL'):
                statid = slnum  # get global variable (all sondes)
                country = []
            else:
                country = [country]
        # wmotable <-
        vcountries = wmotable['sub_region'].alpha_3_code.values
        for icountry in country:
            if icountry not in vcountries:
                raise ValueError('Invalid selection, %s is not a valid country code' % icountry)

            for k, vv in active.items():
                try:
                    if vv[4] == icountry:
                        statid.append(k)
                except:
                    raise ValueError(vv, k)
                    

        if len(statid) == 0:
            raise RuntimeError('Invalid selection, no Stations for Countries %s' % str(country))
    #
    # BBOX [lower left upper right]
    #
    elif bbox is not None:
        # converting from BBOX [upper left lower right] to [lower left upper right] :
        bbox = [float(bbox[2]), float(bbox[1]), float(bbox[0]), float(bbox[3])]
#         print(bbox)
        if not isinstance(bbox, (list, tuple)) or len(bbox) != 4:
            raise ValueError('Invalid selection, bounding box: [upper left lower right]')

        try:
            for i in range(4):
                bbox[i] = float(bbox[i])
        except ValueError:
            raise ValueError('Invalid selection, bounding box: [upper left lower right] must be int or float')

        if bbox[0] >= bbox[2] or bbox[1] >= bbox[3]:
            raise ValueError('Invalid selection, bounding box: lower<upper [-90, 90], left<right [-180, 360]')
        print(bbox)
        if bbox[0] < -90 or bbox[0] > 90 or bbox[2] < -90 or bbox[2] > 90 or bbox[1] < -180 or bbox[1] > 360 or bbox[3] < -180 or bbox[3] > 360 or bbox[3] - bbox[1] > 360:
            raise ValueError('Invalid selection, bounding box: lower<upper [-90, 90], left<right [-180, 360]')
        d['bbox']=bbox
        statid = []
        active_file = config['config_dir'] + '/active.json'
        bbact = json.load(open(active_file,"r"))
        for k, v in bbact.items():
            if bbox[0] <= float(v[2]) <= bbox[2]:
                if bbox[3] <= 180:
                    if bbox[1] <= float(v[3]) <= bbox[3]:
                        statid.append(k)
                else:
                    # rectangle crossing date line
                    if float(v[3]) < 0:
                        if float(v[3]) >= bbox[1] - 360 and float(v[3]) + 360 <= bbox[3]:
                            statid.append(k)
                    else:
                        if bbox[1] <= float(v[3]) <= bbox[3]:
                            statid.append(k)
        if len(statid) == 0:
            raise RuntimeError('Invalid selection, bounding box %s contains no radiosonde stations' % str(bbox))
    #
    # Stations
    #
    else:
        logger.info('start check statid with: %s', statid)
        # '0-20200-0-*': 
#         try:
        if statid == 'all' or statid == None:
            statid = slnum  # <- list of all station ids from init_server

        elif isinstance(statid, (str, int)):
            valid_id = None
            if('*' in statid):
                if statid[:3] == '0-2':
                    stats = []
                    pat=statid[:statid.index('*')]
                    for l in slnum: # -> searches all slnum for matching statids
                        if pat in l: 
                            stats.append(l)
                    valid_id = stats
                else:
                    stats = []
                    for s in ['0-20000-0-', '0-20001-0-', '0-20100-0-', '0-20200-0-', '0-20300-0-']:
                        m = s + statid
                        pat=m[:m.index('*')]
                        for l in slnum: # -> searches all slnum for matching statids
                            if pat in l: 
                                stats.append(l)
                        valid_id = stats

            else:
#                     if not ((len(statid) == 15) or (len(statid) == 5)):
#                         raise ValueError('statid %s of wrong size - please select statid without "0-20..."-prefix of 5 digits, or with "0-20..."-prefix of 15 digits' % str(statid))

                if statid[:2] == '0-' and statid in slnum:
                    valid_id = statid
                else:
                    for s in ['0-20000-0-', '0-20001-0-', '0-20100-0-', '0-20200-0-', '0-20300-0-']:
                        l = s + statid
                        if l in slnum:
                            valid_id = l
                            break

            if valid_id == None:
                raise ValueError('statid not available - please select an area, country or check your statid')

            # if wildcard was used, valid_id is already a list so it can be directly given to statid:
            if isinstance(valid_id, list):
                statid = valid_id
            else:
                statid = [valid_id]

        else:
            valid_id = None
            new_statid = []
            for k in statid:
                logger.info('preparing: %s', k)

                if('*' in k):
                    if k[:3] == '0-2':
                        stats = []
                        pat=k[:k.index('*')]
                        for l in slnum: # -> searches all slnum for matching statids
                            if pat in l: 
                                stats.append(l)
                        valid_id = stats
                    else:
                        stats = []
                        for s in ['0-20000-0-', '0-20001-0-', '0-20100-0-', '0-20200-0-', '0-20300-0-']:
                            m = s + k
                            pat=m[:m.index('*')]
                            for l in slnum: # -> searches all slnum for matching statids
                                if pat in l: 
                                    stats.append(l)
                        valid_id = stats

                else:
#                         if not ((len(k) == 15) or (len(k) == 5)):
#                             raise ValueError('statid %s of wrong size - please select statid without "0-20..."-prefix of 5 digits, or with "0-20..."-prefix of 15 digits' % str(statid))

                    if k[:3] == '0-2' and k in slnum:
                        valid_id = k
                    else:
                        for s in ['0-20000-0-', '0-20001-0-', '0-20100-0-', '0-20200-0-', '0-20300-0-']:
                            l = s + k
                            if l in slnum:
                                valid_id = l
                                break

                # if wildcard was used, valid_id is already a list so it can be directly given to new_statid:
                logger.info('final newstatid: %s   and valid_id: %s',new_statid, valid_id)
                if isinstance(valid_id, list):
                    for ns in valid_id:
                        new_statid.append(ns)
                else:
                    new_statid.append(valid_id)
                logger.info('after stat con: %s', new_statid)

            if valid_id == None:
                raise ValueError('statid not available - please select an area, country or check your statid')

            statid = [] 
            [statid.append(x) for x in new_statid if x not in statid] 
#         except Exception:
#             raise RuntimeError(
#                 'Invalid selection, specify either bbox, country or statid. Use "statid":"all" to select all ' \
#                 'stations ')
    d['statid'] = statid
    logger.info('final station selection: %s, %s', len(d['statid']), d['statid'])
    if len(d['statid']) == 1 and ((single_parallel == True) or (single_parallel == 'True')):
        # d['single_parallel'] = True
        d['single_parallel'] = False
        ##########
        ########## Turn back on!
        ##########
    #
    #
    # remove statids i
    
    #
    # single_parallel
    #
    if (single_parallel == False) or (single_parallel == 'False') or (len(d['statid']) > 1):
        d['single_parallel'] = False
#     if (single_parallel == True) or (single_parallel == 'True'):
#         d['single_parallel'] = True
#     else:
#         d['single_parallel'] = False
    print('single_parallel: ', d['single_parallel'])

    #
    # single_csv
    #
    if ((single_csv == True) or (single_csv == 'True')) and (d['format'] == 'csv'):
        d['single_csv'] = True
    else:
        d['single_csv'] = False
#     if (single_parallel == True) or (single_parallel == 'True'):
#         d['single_parallel'] = True
#     else:
#         d['single_parallel'] = False

    
    #
    #
    # Only pick one format for dates:
    date_not_yet_existing = True
    # prioritized order:, date, day/month/year
    # Period removed -> cds always converts period to date in this format: '19990101-20000101'
    #
    
    #
    # Date time selection
    # [DATE] or [START, END]
    #
    # todo only one date or two dates allowed at the moment by CDS
    if date is not None and date_not_yet_existing:
        # str, list (str, int)
        newdates = []
        # make a list
        if isinstance(date, (int, str)):
            date = [date]
        # BUG: date list is not sorted from CDS -> sort here
        date.sort()
        for idate in date:
            # convert to string
            if not isinstance(idate, str):
                idate = str(idate)
            # todo old style date range (not supported by CDS anymore ???)
            if '-' in idate:
                idate = idate.split('-')
                if idate[0] > idate[1]:
                    raise ValueError('starting date has to be before ending date: %s - %s' % (idate[0], idate[-1]))
                # check period dates (should not be out of range)
                if int(idate[0][-2:]) > 31 or int(idate[-1][-2:]) > 31:
                    raise ValueError('only valid dates allowed for date: %s' % idate)
                newdates.append(to_valid_datetime(idate[0], as_string=True))
                newdates.append(to_valid_datetime(idate[-1], as_string=True))
                # newdates.append('{}-{}'.format(to_valid_datetime(idate[0], as_string=True),
                #                                to_valid_datetime(idate[-1], as_string=True)))
            else:
                if int(idate[-2:]) > 31:
                    raise ValueError('only valid dates allowed for date: %s' % idate)
                try:
                    newdates.append(to_valid_datetime(idate, as_string=True))
                except:
                    raise ValueError('only valid dates allowed for date: %s' % idate)
        d['date'] = newdates
        date_not_yet_existing = False
    #
    # day/month/year selection
    #
    if year is not None and date_not_yet_existing:
        sday = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12',
                '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24',
                '25', '26', '27', '28', '29', '30', '31']
        smonth = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
        datelist = []
        newdates = []
        if isinstance(year, (int, str)):
            year = [year]
        if month is None:
            month = smonth
        elif isinstance(month, (int, str)):
            month = [month]
        if day is None:
            day=sday
        elif isinstance(day, (int, str)):
            day = [day]
        # for removal of e.g. 31.02.:
        for i in year:
            for j in month:
                if not (str(j).zfill(2) in smonth):
                    raise ValueError('not a valid month: %s' % j)
                for k in day:
                    if not (str(k).zfill(2) in sday):
                        raise ValueError('not a valid day: %s' % k)
                    try:
                        datetime.strptime(str(i)+str(j).zfill(2)+str(k).zfill(2), '%Y%m%d')
                        datelist.append(str(i)+str(j).zfill(2)+str(k).zfill(2))
                    except:
                        pass
        datelist.sort()
        newdates.append(to_valid_datetime(datelist[0], as_string=True))
        newdates.append(to_valid_datetime(datelist[-1], as_string=True))
        d['date'] = newdates
        date_not_yet_existing = False
        
    #
    # no given date:
    #
    if date is None and date_not_yet_existing:
        d['date'] = ['19000101', '21000101']
        date_not_yet_existing = False
    
    #
    # Pressure levels
    #
    if pressure_level is not None:
        if not isinstance(pressure_level, list):
            pressure_level = [pressure_level]

        for i in range(len(pressure_level)):
            # need to be string
            pressure_level[i] = str(pressure_level[i])
            try:
                pressure_level[i] = int(pressure_level[i])  # if not integer raises a ValueError
            except:
                raise TypeError('invalid selection, pressure_level allows only integer, ' + pressure_level[i])
            # pressure should be in Pa and between
            if pressure_level[i] < 500 or pressure_level[i] > 110000:
                raise ValueError(
                    'invalid selection, pressure_level out of range [5-1100 hPa]: %d Pa' % pressure_level[i])
            # need to be a string for processing
            pressure_level[i] = str(pressure_level[i])
        d['pressure_level'] = pressure_level
    #
    # times
    #
    if (time is not None) and (time != '') and (time != []):
        if not isinstance(time, list):
            time = [time]
        for i in range(len(time)):
            try:
                # need to be string
                time[i] = str(time[i])
                # in hours
                if int(time[i]) < 0 or int(time[i]) > 24:
                    raise ValueError('invalid selection, time out of range [0-24 h]: %d' % int(time[i]))
            except:
                raise ValueError('invalid selection, time allows only integer, ' + time[i])
        #d['da'] = False #leo
        d['time']=time

    return d


###############################################################################


def makebodiesold(bodies, body, spv, bo, l):
    """ Split Request by Variable and Station ID -> MP

    Args:
        bodies: list of requests
        body: request dictionary
        spv: variables to split by
        bo: tmp var ?
        l: 0

    Returns:
        bodies
    """
    # todo this fucntion is wired ... redesign needed ?
    for b in body[spv[l]]:
        if l < len(spv) - 1:
            makebodies(bodies, body, spv, copy.copy(bo) + [b], l + 1)
        else:
            #bodies.append(dict(body))
            bodies.append({})
            bn = copy.copy(bo) + [b]
            for k,v in body.items():
                if k not in spv:
                    bodies[-1][k]=v
            for s, b in zip(spv, bn):
                bodies[-1][s] = b
                logger.debug('makebodies %d %s %s', l, s, b)
    return

def makebodies(bodies, body, spv, bo, l):
    """ Split Request by Variable and Station ID -> MP

    Args:
        bodies: list of requests
        body: request dictionary
        spv: variables to split by
        bo: tmp var ?
        l: 0

    Returns:
        bodies
    """
    
    sbody={}
    for k,v in body.items():
        if k not in spv:
            sbody[k]=v
    
    tt=time.time()
    p=product(*[body[s] for s in spv])
    print(time.time()-tt)
    l=0
    for x in p:
        bodies.append(dict(sbody))
        for k,v in zip(spv,x):
            bodies[-1][k]=v
            #logger.debug('makebodies %d %s %s', l, k, v)
            l+=1
    print(time.time()-tt)
    return


def process_request(body: dict, output_dir: str, wmotable: dict, P, debug: bool = False) -> str:
    """ Main function of the hug server

    Args:
        body: request dictionary
        output_dir: Output directory
        wmotable: WMO regions definitions from init_server()
        debug: for debugging

    Returns:
        str : filename of zipped requested files
    """
    tt = time.time()
    #
    # Raises Errors will be handled by base_exception_handler of the hug server
    # in debug this will give a traceback
    #
    
    body = check_body(wmotable=wmotable, **body)
    #
    # check if the request isn't too long
    #
    time_len = ((int(body['date'][-1][:4])-int(body['date'][0][:4]))*12 + # variables and stations
                (int(body['date'][-1][4:6])-int(body['date'][0][4:6]))) # months
    if time_len == 0:
        time_len = 1
    lenprod = len(body['variable']) * len(body['statid']) * time_len
#     if len(body['variable']) == 1 and ((int(body['date'][-1][:4])-int(body['date'][0][:4]))*12 + (int(body['date'][-1][4:6])-int(body['date'][0][4:6]))) == 1:
#         logger.warning('Requesting too many elements - Exception: please split your request')
    request_logger.info('REQUEST ELEMENTS: %s', lenprod)
    if lenprod > 300000:
        # lenght restriction deactivated as long following line is out commented.
        raise RuntimeError('Request too large - please split your request')
#         logger.warning('Request very large - please split')
    #
    logger.debug('Cleaned Request %s', str(body))
    os.makedirs(output_dir, exist_ok=True)  # double check
    bodies = []
    spv = ['statid', 'variable']  # potential split up variables
    bo = []
    refdate = datetime(year=1900, month=1, day=1)
    print(' vor makebodies',time.time()-tt)
    if 'date' in body.keys():      
        start = (to_valid_datetime(body['date'][0]) - refdate).days * 86400
        ende = (to_valid_datetime(body['date'][-1]) - refdate).days * 86400+86399
    else:
        start=1
        ende=5000000000

    tt=time.time()
    if len(body['statid'])>1:     
        gdict2,lidx=eua.searchdate(active['rtsidx'], active['rtsarr'], start,ende)    
        gdict2=dict(zip([active['rtskeys'][l] for l in lidx],list(gdict2)))
        gd={}
        for b in body['statid']:
            try:
                
                gd[b]=gdict2[b]
            except:
                pass
        gdict2=gd
        body['statid']=gd

    else:
        idx=active['rtskeys'].index(body['statid'][0])
        gdict2,lidx=eua.searchdate(active['rtsidx'][idx:idx+2], active['rtsarr'], start,ende)    
        gdict2=dict(zip(body['statid'],list(gdict2)))
        
        

    print(time.time()-tt)
    makebodies(bodies, body, spv, bo, 0)  # List of split requests
    print('bodies: ', bodies)
    print('gdict2: ', gdict2)
    for k in range(len(bodies)):
        key=bodies[k]['statid']
        bodies[k]['filename']=active[key][5]
        bodies[k]['rtsidx']=gdict2[key]    #
    # Check all requests if dates are within active station limits
    #
    #activekeys=list(active.keys())
    #
    # seconds since Reference date
    #
    print('makebodies',time.time()-tt)
    
    #for k in range(len(bodies) - 1, -1, -1):
        ## date selection ? do all the stations have data in this period?
        #if 'date' in bodies[k].keys():
            ##
            ## Station Active ?
            ##
            
            #try: #if bodies[k]['statid'] in activekeys:
                #i=activekeys.index(bodies[k]['statid'])
                #idx=numpy.searchsorted(active['rtsarr'][active['rtsidx'][i]:active['rtsidx'][i+1]],(start,ende))
                #if idx[0]>=idx[1]: #rtsidx[i+1]-rtsidx[i]:
                ##if start > active[bodies[k]['statid']][1] or ende + 86399 < active[bodies[k]['statid']][0]:
                    #logger.debug('%s outside Index range', bodies[k]['statid'])
                    #del bodies[k]

                #else:
                    #logger.debug('%s Index[%d (%d / %d) %d]', bodies[k]['statid'],
                                 #active[bodies[k]['statid']][0],
                                 #start, ende,
                                 #active[bodies[k]['statid']][1])
                    ## add data path to request
                    ## input_dirs.append(active[bodies[k]['statid']][5])  # path from makedaterange (init_server)
                    ## file_paths.insert(0, active[bodies[k]['statid']][5])  # path from makedaterange (init_server)
                    #bodies[k]['filename'] = active[bodies[k]['statid']][5]
                    #bodies[k]['rtsidx']=idx[:]#(active['rtsidx'][i]+idx[0],active['rtsidx'][i]+idx[1])
            #except:
                #del bodies[k]
        #else:
            ## input_dirs.append(active[bodies[k]['statid']][5])  # path from makedaterange (init_server)
            #bodies[k]['filename'] = active[bodies[k]['statid']][5]
            #i=activekeys.index(bodies[k]['statid'])
            #idx=numpy.searchsorted(active['rtsarr'][active['rtsidx'][i]:active['rtsidx'][i+1]],(1,5000000000))
            #bodies[k]['rtsidx']=idx[:]#(active['rtsidx'][i]+idx[0],active['rtsidx'][i]+idx[1])
            ## file_paths.insert(0, active[bodies[k]['statid']][5])  # path from makedaterange (init_server)

    logger.debug('# requests %d', len(bodies))
    if len(bodies) == 0:
        raise RuntimeError('No selected station has data in specified date range: ' + str(body))

    # Make process_flat a function of only request_variables (dict)
    #
    # process_flat(outputdir: str, cftable: dict, debug:bool=False, request_variables: dict) -> tuple:
    # func = partial(eua.process_flat, output_dir, cf, input_dirs[0], debug)
#     if body['hdf']:
#         randdir = '{:012d}'.format(numpy.random.randint(100000000000))
#         rfile = os.path.dirname(output_dir) + '/request_'+randdir+'.nc'
#         with h5py.File(rfile, 'a') as merge:
#             for k in range(len(bodies)):
#                 name=bodies[k]['statid']
#                 merge.create_group(name.split('-')[-1])

# #         global OPEN_FILE
# #         OPEN_FILE =  h5py.File(rfile, 'a')
# #         for k in range(len(bodies)):
# #             name=bodies[k]['statid']
# #             OPEN_FILE.create_group(name.split('-')[-1])

#         func = partial(eua.process_flat, rfile, cf, debug)
# #         OPEN_FILE.close()
#     else:
#         func = partial(eua.process_flat, output_dir, cf, debug)

## DONT FORGET TO CHANGE IN PROCESS_FLAT AND READ_WRITE too!

    func = partial(eua.process_flat, output_dir, cf, debug)

    print('body', time.time()-tt)
    if 'gridded' in body:
        body['statid']=''
        print('body', body)
        results = [eua.process_flat(outputdir = output_dir, cftable = cf, debug = True, request_variables = body)]
    #
    # Smaller request?
    #
    
    elif debug or len(body['variable']) * len(body['statid'])<10:
        #
        # Single Threading
        #
        print('single threading')
        results = list(map(func, bodies))
    else:
        
        #
        # Multi Threading
        #
        print(time.time()-tt); tt=time.time()
        #with Pool(10) as p:
            # error with chunksize (from p.map to p.starmap)
        #bdist=active['rtsidx'][1:]-active['rtsidx'][:-1]
        #bidx=np.argsort(bdist)
        #for i in range(10):
            #print(np.sum(bdist[bidx][i::10]))
        
        #b2=[]
        #for b in reversed(bidx):
            #b2.append(bodies[b])
        with Pool(9) as Pl:
            results = list(Pl.map(func, bodies,chunksize=3))
        #results = list(map(func, bodies))
            # results = list(p.starmap(func, zip(input_dirs, [debug]*len(bodies), bodies), chunksize=1))
    #
    # Process the output 
    # todo catch Error Messages and store in a log file?
    #
    print(time.time()-tt,results)
    wpath = ''  # same as output_dir ?
    for r in results:
        print(r)
        if r[0] != '':
            wpath = r[0]
            break
            
    if body['hdf']:
        
        randdir = '{:012d}'.format(numpy.random.randint(100000000000))
        rfile = os.path.dirname(output_dir) + '/request_'+randdir+'.nc'
        print(rfile)
#         rfile = os.path.dirname(wpath) + '/download.nc'
        with h5py.File(rfile, 'w') as merge:
            for r in results:
                if(r[0]):                    
                    name = r[0].split('dest_')[-1].split('_')[0]
                    merge.create_group(name)
                    with h5py.File(r[0], 'r') as filetocopy:
                        for i in filetocopy.keys():
                            filetocopy.copy(i, merge[name], name=i)
        logger.debug('netcdfs merged [%d] to %s', len(results), rfile)
        
#         rfile = os.path.dirname(wpath) + '/download.nc'
#         with h5py.File(rfile, 'w') as merge:
#             for r in results:
#                 if(r[0]):                    
#                     name = r[0].split('dest_')[-1].split('_')[0]
#                     merge.create_group(name)
#                     with h5py.File(r[0], 'r') as filetocopy:
#                         for i in filetocopy.keys():
#                             filetocopy.copy(i, merge[name], name=i)
#         logger.debug('netcdfs merged [%d] to %s', len(results), rfile)
        
    else:
#     if not body['hdf']: 
        if wpath == '':
            raise RuntimeError('Error: %s (%s)' % (results[0][1], str(body)))
        else:
            rfile = os.path.dirname(wpath) + '/download.zip'
            print(rfile)

        logger.debug('wpath: %s; format %s Time %f', wpath, body['format'],time.time()-tt)

        if body['single_csv']:
            print('start writing single_csv')
            write_results = []
            for i_res in results:
                if i_res[0] != '':
                    write_results.append(i_res[0])
            # to place the stations in an order manner
            write_results.sort()
            ### <
#             combined_csv = pd.concat([pd.read_csv(f, header=14) for f in write_results], ignore_index=True)
            ### >
                
            geo_ll = []
            var_ll = []
            for geo in write_results:
#                 with gzip.open(geo, 'rt') as fd:
                with open(geo, 'rt') as fd:
                    reader = csv.reader(fd)
                    rows = [next(reader) for x in range(15)]
                    #rows = list(reader)
                    geo_row = []
                    exmn = rows[9][0]
                    exmn = exmn.split(' [South_West_North_East]')[0]
                    exmn = exmn.split('# Geographic area: ')[-1]
                    exmn =  exmn.split('_')
                    for row_ll in exmn:
                        geo_row.append(float(row_ll))
                    geo_ll.append(geo_row)
                    
                    var_ll.append(rows[10][0])
            var_ll = ", ".join(np.unique(var_ll))
            geo_ll = np.array(geo_ll)
#             results = [(''.join([i+'/' for i in rfile.split('/')[:-1]])+body['single_csv_target']+".csv.gz", '')]
            results = [(''.join([i+'/' for i in rfile.split('/')[:-1]])+body['single_csv_target']+".csv", '')]

#             with gzip.open(results[0][0], 'w') as file:
            with open(results[0][0], 'w') as file:
#                 with gzip.open(write_results[0], 'r') as f:
                with open(write_results[0], 'r') as f:
                    for i in range(15): ### < 15 >
                        if i == 10:
                            file.write(str('# Variables selected: ' + str(var_ll) + ' \n')) # str().encode()
                            f.readline()
                        elif i == 9:
                            file.write(str('# Geographic area: ' + str([np.min(geo_ll[:,0]),np.min(geo_ll[:,1]),np.max(geo_ll[:,2]),np.max(geo_ll[:,3])]) +
                                           ' [South_West_North_East] \n')) # str().encode()
                            f.readline()
                        else:
    #                         # for csv output str is needed
    #                         file.write(f.readline().decode())
                            # for gzip output byte is needed
                            file.write(f.readline())
#             combined_csv.to_csv(results[0][0], index=False, mode="a", compression='gzip')
            ### <
#             combined_csv.to_csv(results[0][0], index=False, mode="a")
            ### >
            
            print('pdconcat', time.time() - tt)
        
            with open(results[0][0]+'2', 'w') as fd:
                p1 = subprocess.Popen(["cat", results[0][0], results[0][0]], stdout=fd)
        
            print('input files: ', results)
            write_results.insert(0, results[0][0]+'2')
            p1 = subprocess.Popen(["tail", '-n', '+16', '-q']+write_results, stdout=subprocess.PIPE)
            with open(rfile, 'wb') as fout:
                print(['pigz', '-A', body['single_csv_target'] + '.csv', '--zip'])
                p2 = subprocess.run(['pigz', '-A', body['single_csv_target'] + '.csv', '--zip'], stdin=p1.stdout, stdout=fout)
            print('pdconcat',time.time() - tt, p2)
           
        else:
            if 'local_execution' in body.keys():
                return rfile

            #tt=time.time()

            with zipfile.ZipFile(rfile, 'w' ) as f:
                for r in results:
                    try:
                        if len(r[0]) > 0:
                            f.write(r[0], os.path.basename(r[0]))
                        if debug:
                            continue  # do not remove
                        os.remove(r[0])  # remove NetCDF file
                    except:
                        pass
            logger.debug('files compressed [%d] to %s', len(results), rfile)
            print('ZIPPING: ', tt - time.time())

        
#     else:
#         with zipfile.ZipFile(rfile, 'w', compression=zipfile.ZIP_DEFLATED) as f:
#             for v in body['variable']:
#                 ilist = glob.glob(output_dir + '/*' + v + '.nc')
#                 if len(ilist) > 0:
#                     ifile = to_csv(ilist, ofile=output_dir + '/' + v + '.csv')  # todo add correct name into zip
#                     f.write(ifile, os.path.basename(ifile))
#                     logger.debug('writing csv %s [%d] to %s', v, len(ilist), rfile)
#                     if debug:
#                         continue  # do not remove
#                     for i in ilist:
#                         os.remove(i)  # remove NetCDF file
#                     os.remove(ifile)  # remove csv

    logger.debug('Request-File: %s [Time: %7.4f s]', rfile, (time.time() - tt))
    request_logger.info('REQUEST TIME: %s', (time.time() - tt))
    return rfile


@app.get('/', response_class=FileResponse)
# def index(request:dict):
def index(request=None, response=None):
    """ Main Hug Index Function on get requests

    index function requests get URI and converts into dictionary.

    Args:
        request: dictionary

    Returns:

    """
#     logger.debug("GET %s", request.query_string)
#     if '=' not in request.query_string:
#         response.status = hug.HTTP_422
#         raise Exception('A query string must be supplied')

    try:
        rs = request.split('&') # rs = request.query_string.split('&')
        body = {}
        for r in rs:
            k, v = r.split('=')
            if '[' in v:
                vl = v[1:-1]
                if k in ['statid', 'variable', 'fbstats']:
                    body[k] = vl.split(',')
                else:
                    body[k] = list(numpy.fromstring(vl, sep=',', dtype='int'))
            else:
                body[k] = v
    except:
#         response.status = hug.HTTP_422
        raise Exception(request.query_string)
    logger.debug("BODY %s", body)

    randdir = '{:012d}'.format(numpy.random.randint(100000000000))
    logger.info("%s GET %s", randdir, str(body))
    tmpdir = config['tmp_dir'] + '/' + randdir
    try:
        # rfile = process_request(body, tmpdir, config['data_dir'], wmo_regions)
        rfile = process_request(body, tmpdir, wmo_regions, P, debug=config['debug'])
    except Exception as e:
        logger.error("%s GET FAILED, %s", randdir, e)
        with open(config['logger_dir'] + '/failed_requests.log', 'a+') as ff:
            ff.write('%s - %s [%s] Message: %s \n' % (str(datetime.now()), randdir, str(body), e))

        logger.info("%s GET FAILED %s", randdir, e)
        raise e
    logger.info("%s GET FINISHED", randdir)
    #
    # Write successful requests
    #
    with open(config['logger_dir'] + '/finished_requests.log', 'a+') as ff:
        ff.write('%s - %s [%s] \n' % (str(datetime.now()), randdir, str(body)))

#     response.set_header('Content-Disposition', 'attachment; filename=' + os.path.basename(rfile))
    return rfile


# @app.exception(Exception)
# def base_exception_handler(exception, response=None):
#     """ This captures any Exception from the Server

#     Args:
#         exception: hug.exceptions.Exception class
#         response: HTTP Response

#     Returns:
#         str :  Message
#     """
#     response.status = hug.HTTP_422  # Always the same Error ? or dependent on exception?
#     # one of the arguments could be the response status
#     return ",".join(exception.args)


@app.post('/', response_class=FileResponse)
async def index(body:Request, request=None, response=None):

    """ Main Hug index function for Post requests

    Args:
        body: dictionary of request

    Returns:

    Body:
     - variable         – e.g. temperature, ... as list or string
     - statid           – e.g. 01001 as list or string
     - product_type     – currently unused
     - pressure_level   – 500, ...
     - date             – '19990131' or ['19990101', '20000101'] as range
     - time             – 1,... or 0-23
     - fbstats          – only these are currently allowed: 'obs_minus_an', 'obs_minus_bg', 'bias_estimate'
     - bbox             – Bounding Box [lower left upper right]
     - country          – Country Code, e.g. DEU, AUT, USA, GBR
     - format           – 'nc' or 'csv'
     - period           – ['19990101', '20000101']

    """
    

    hostip = body.client.host
    clientip = body.client[0] + ':' + str(body.client[1])

    body_str = (await body.body()).decode()
    logger.info('REQUEST STRING: %s', body_str)
    
    request_logger.info('REQUEST STRING: %s', body_str)
    request_logger.info('body.client.host IP: %s', hostip)
    request_logger.info('body.client[0] IP: %s', clientip)

    if '&' in body_str:	
       	body = urllib.parse.parse_qs(body_str)
    else:
        body = json.loads(body_str)

    for i in body:
        if len(body[i]) == 1:
            body[i] = body[i][0]

    logger.info('ASYNC INDEX POST')
    randdir = '{:012d}'.format(numpy.random.randint(100000000000))
    logger.info("%s BODY %s %s", randdir,type(body),body)
    tmpdir = config['tmp_dir'] + '/' + randdir
    try:
        #rfile='/fio/srvx7/leo/x'
        logger.error("BODY: %s", body)
        rfile = process_request(body, tmpdir, wmo_regions, P, debug=config['debug'])
    except Exception as e:
        logger.error("%s POST FAILED, %s", randdir, e)
        with open(config['logger_dir'] + '/failed_requests.log', 'a+') as ff:
            ff.write('%s - %s [%s] Message: %s \n' % (str(datetime.now()), randdir, str(body), e))

        logger.info("%s POST FAILED %s", randdir, e)
        raise e
    logger.info("%s POST FINISHED", randdir)
    #
    # Write successful requests
    #
    with open(config['logger_dir'] + '/finished_requests.log', 'a+') as ff:
        ff.write('%s - %s [%s] \n' % (str(datetime.now()), randdir, str(body)))

#     response.set_header('Content-Disposition', 'attachment; filename=' + os.path.basename(rfile))
    return rfile


@app.get('/constraints/', response_class=FileResponse)
def index(request=None, response=None):
    """ Main Hug Index Function on get requests

    index function requests get URI and converts into dictionary.

    Args:
        request: dictionary
        response: str

    Returns:

    """
    logger.debug("GET %s", request.query_string)


    rfile='/data/public/tmp/constraints.csv'

    return rfile

@app.get('/variable_codes/', response_class=FileResponse)
def variable_codes(request=None, response=None):
    """ 
        Returns a json with all the variable codes
    """
    logger.debug("GET variable_codes")
    rfile='/data/public/tmp/variable_codes.json'
    with open(rfile, 'w') as fp:
        json.dump(eua.cdm_codes, fp)
#     response.set_header('Content-Disposition', 'attachment; filename=' + rfile)
    return rfile


@app.get('/variable_definition/', response_class=FileResponse)
def variable_definition(variable=None, response=None):
    """ 
        Returns a json with the variable definition
        Args:
            variable: str - variable of interest
    """
    logger.debug("GET variable_definition")
    if variable == None:
        rfile=config['config_dir']+'/cf.json'
#         response.set_header('Content-Disposition', 'attachment; filename=' + os.path.basename(rfile))
        return rfile
    else:
        cf = json.load(open(config['config_dir']+'/cf.json', 'r'))
        rfile='/data/public/tmp/variable_definition_'+variable+'.json'
        with open(rfile, 'w') as fp:
            json.dump(cf[variable], fp)
#         response.set_header('Content-Disposition', 'attachment; filename=' + os.path.basename(rfile))
        return rfile


@app.get('/station_configuration/', response_class=FileResponse)
def station_configuration(response=None):
    """ 
        Returns a the station configuration file for all merged stations
        Args:
            -
    """
    logger.debug("GET station_configuration")
    rfile = './CUON_station_configuration_extended.csv'
    return rfile


def datetime_to_seconds(dates, ref='1900-01-01T00:00:00'):
    """ from datetime64 to seconds since 1900-01-01 00:00:00"""
    return ((numpy.datetime64(dates) - numpy.datetime64(ref)) / numpy.timedelta64(1, 's')).astype(numpy.int64)

@app.get('/maplist/', response_class=FileResponse)
def mapdata(date=None, enddate=None, var=85, response=None,):
    """ Main Hug Index Function on get requests

    index function requests get URI and converts into dictionary.

    Args:
        request: dictionary
        response: str

    Returns:

    """
    active_file = config['config_dir'] + '/active.json'
    act = json.load(open(active_file,"r"))
    
    output_file = '/tmp/maplist_'+str(date)
    if (enddate is None) or (date == enddate):
        reqdate = date.split('-')
        interm_file = '/data/public/constraints_by_date/'+str(var)+'/'+str(var)+'_'+reqdate[0]+'_'+str(int(reqdate[1]))+'_'+str(int(reqdate[2]))+'.csv'
        copyfile(interm_file, output_file)
        with open(output_file) as f:
            lines = f.readlines()
        lines[0] = "station_name,longitude,latitude\n"
        with open(output_file, "w") as f:
            f.writelines(lines)
            
    else:    
        date = datetime_to_seconds(date)
        enddate = datetime_to_seconds(enddate)
        rows = []
        rows.append(['station_name', 'longitude', 'latitude'])
        for i in act:
            if (date >= act[i][0]) and (enddate <= act[i][1]):
                # renaming deactivated for now
                # name = namelist[i]
                name = i
                rows.append([name, act[i][3], act[i][2]])

        with open(output_file, 'w') as csvfile:  
            # creating a csv writer object  
            csvwriter = csv.writer(csvfile)  
            # writing the data rows  
            csvwriter.writerows(rows)

    return output_file


@app.post('/maplist_post/', response_class=FileResponse)
def mapdata(date=None, enddate=None, var=85, response=None,):
    """ Main Hug Index Function on get requests

    index function requests get URI and converts into dictionary.

    Args:
        request: dictionary
        response: str

    Returns:

    """
    active_file = config['config_dir'] + '/active.json'
    act = json.load(open(active_file,"r"))
    
    output_file = '/tmp/maplist_'+str(date)
    if (enddate is None) or (date == enddate):
        reqdate = date.split('-')
        interm_file = '/data/public/constraints_by_date/'+str(var)+'/'+str(var)+'_'+reqdate[0]+'_'+str(int(reqdate[1]))+'_'+str(int(reqdate[2]))+'.csv'
        copyfile(interm_file, output_file)
        with open(output_file) as f:
            lines = f.readlines()
        lines[0] = "station_name,longitude,latitude\n"
        with open(output_file, "w") as f:
            f.writelines(lines)
    else:    
        date = datetime_to_seconds(date)
        enddate = datetime_to_seconds(enddate)
        rows = []
        rows.append(['station_name', 'longitude', 'latitude'])
        for i in act:
            if (date >= act[i][0]) and (enddate <= act[i][1]):
                # renaming deactivated for now
                # name = namelist[i]
                name = i
                rows.append([name, act[i][3], act[i][2]])

        with open(output_file, 'w') as csvfile:  
            # creating a csv writer object  
            csvwriter = csv.writer(csvfile)  
            # writing the data rows  
            csvwriter.writerows(rows)

    return output_file

@app.get('/maplist2/', response_class=FileResponse)
def mapdata2(date=None, plev=None, var=None, response=None):
    """ Main Hug Index Function on get requests

    index function requests get URI and converts into dictionary.

    Args:
        request: dictionary
        response: str

    Returns:

    """
    
    yr, mn, dy = date.split('-')
    const = constraints[constraints.observed_variable == int(var)]
    const = const[const.z_coordinate == int(plev)]
    const = const[const.year == int(yr)]
    const = const[const.month == int(mn)]
    const = const[const.day == int(dy)]
    const = const.drop_duplicates(subset=['lat', 'lon'], keep='last')
    
    output_file = '/tmp/maplist_'+str(date)
    rows = []
    rows.append(['station_name', 'longitude', 'latitude'])
    for i in len(const):
        rw = const.iloc[i]
        rows.append([str(rw.station_name), float(rw.lon), float(rw.lat)])

    with open(output_file, 'w') as csvfile:  
        # creating a csv writer object  
        csvwriter = csv.writer(csvfile)  
        # writing the data rows  
        csvwriter.writerows(rows) 
            
    return output_file


@app.get('/statlist/', response_class=FileResponse)
def statdata(date:str=None, mindate:str=None, enddate:str=None, response:str=None):
    """ Main Hug Index Function on get requests

    index function requests get URI and converts into dictionary.

    Args:
        request: dictionary
        response: str

    Returns:

    """
    active_file = config['config_dir'] + '/active.json'
    act = json.load(open(active_file,"r"))
    
    output_file = '/tmp/maplist_'+str(date)
    
    if enddate is None:
        date = datetime_to_seconds(date)
        rows = []
        rows.append(['station', 'longitude', 'latitude'])
        for i in act:
            if (date >= act[i][0]) and (date <= act[i][1]):
#                 name = namelist[i]
                rows.append([i, act[i][3], act[i][2]])

        with open(output_file, 'w') as csvfile:  
            # creating a csv writer object  
            csvwriter = csv.writer(csvfile)  
            # writing the data rows  
            csvwriter.writerows(rows) 
            
    elif ((not enddate is None) and (not mindate is None)):
        mindate = datetime_to_seconds(mindate)
        enddate = datetime_to_seconds(enddate)
        rows = []
        rows.append(['station', 'longitude', 'latitude'])
        for i in act:
            if (((mindate >= act[i][0]) and (mindate <= act[i][1])) or
                ((enddate >= act[i][0]) and (enddate <= act[i][1])) or
                ((mindate <= act[i][0]) and (enddate >= act[i][1])) 
               ):
#                 name = namelist[i]
                rows.append([i, act[i][3], act[i][2]])

        with open(output_file, 'w') as csvfile:  
            # creating a csv writer object  
            csvwriter = csv.writer(csvfile)  
            # writing the data rows  
            csvwriter.writerows(rows)
            
    elif not enddate is None:
        date = datetime_to_seconds(date)
        enddate = datetime_to_seconds(enddate)
        rows = []
        rows.append(['station', 'longitude', 'latitude'])
        for i in act:
            if (date >= act[i][0]) and (enddate <= act[i][1]):
#                 name = namelist[i]
                rows.append([i, act[i][3], act[i][2]])

        with open(output_file, 'w') as csvfile:  
            # creating a csv writer object  
            csvwriter = csv.writer(csvfile)  
            # writing the data rows  
            csvwriter.writerows(rows)
            
    return output_file



if __name__ == '__main__':
    #active, wmo_regions, cf = init_server()
    #
    # Parse command line arguments for testing the server API
    #
    if len(sys.argv) == 1:
        print(
            "python default.py ""{'variable':['temperature'],'date':['20000101','20190131'], 'pressure_level': 500}""")
        sys.exit(0)
    body = eval(sys.argv[1])
    debug = body.pop('debug', False)
    if 'status' in body.keys():
        print(status_test(command=body['status']))
        sys.exit()
    #
    # Logging to DEBUG to std.out and hug.debug.local.log
    #
    logger.setLevel(10)
    for i in logger.handlers:
        i.setLevel(10)
    #
    # Specific directory for testing and clean it if necessary
    #
    randdir = os.path.expandvars('{:012d}'.format(100000000000))
    if os.path.isdir(randdir):
        for ifile in os.listdir(randdir):
            print(randdir + '/' + ifile)
            os.remove(randdir + '/' + ifile)

    logger.debug(str(body))
    #
    # Run the request
    #
    tmpdir = config['tmp_dir'] + '/' + randdir
#     with Pool(16) as P:
    ret = process_request(body, tmpdir, wmo_regions, P, debug=debug)
    logger.debug(str(ret))
