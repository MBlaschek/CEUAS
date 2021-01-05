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
import json
import logging
import os
import socket
import sys
import time
import zipfile
from datetime import datetime, timedelta
from functools import partial
from multiprocessing import set_start_method, Pool
from typing import Union

if False:
    import cds_eua2 as eua  # old version
    CDS_EUA_VERSION = 2
else:
    #sys.path.append(os.path.expanduser('~leo/python/'))
    import cds_eua3 as eua  # new version with CDMDataset class
    CDS_EUA_VERSION = 3

import h5py
import hug
import numpy
import pandas as pd
import xarray

try:
#    set_start_method("spawn")  # or fork ? not sure why, pickling?
    set_start_method("forkserver")  # or fork ? not sure why, pickling?
except RuntimeError:
    pass

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
          'tmp_dir': '.',
          'config_dir': './config',
          'debug': False,
          'reload_pwd': 'reload'}

if os.path.isfile(config_file):
    new = json.load(open(config_file, 'r'))
    # expand ENVIRONMENTAL VARIABLES if any
    for ikey,ival in new.items():
        if 'dir' in ikey:
            new[ikey] = os.path.expandvars(ival)
            if '$' in new[ikey]:
                raise RuntimeError('configuration path expansion failed %' % new[ikey])
    config.update(new)
else:
    print("Writing new config file:", config_file, "Adjust accordingly!")
    json.dump(config, open(config_file, 'w'))

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
try:
    # todo this part is deprecated / LEO delete?
    if 'srvx' in host:
        sys.path.append(os.path.expanduser('~leo/python/'))
        config['data_dir'] = os.environ["RSCRATCH"]  # ?
except:
    pass


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
                if CDS_EUA_VERSION == 3:
                    #if isinstance(f['recordindices'], h5py._hl.group.Group):
                    if 'recordindices' in f.keys():
                        active[skey] = [int(eua.to_seconds_since(f['recordindices']['recordtimestamp'][0], funits)),
                                    int(eua.to_seconds_since(f['recordindices']['recordtimestamp'][-1], funits)),
                                    float(f['observations_table']['latitude'][-1]),
                                    float(f['observations_table']['longitude'][-1])]
                    else:
                        active[skey] = [int(eua.to_seconds_since(f['recordtimestamp'][0], funits)),
                                        int(eua.to_seconds_since(f['recordtimestamp'][-1], funits)),
                                        float(f['observations_table']['latitude'][-1]),
                                        float(f['observations_table']['longitude'][-1])]
                idx = numpy.where(vola.StationId.values == skey)[0]
                if len(idx) > 0:
                    active[skey].append(vola.CountryCode[idx[0]])
                else:
                    active[skey].append('')
                    logger.debug('no key found for %s', skey)
                # add data directory for process_flat
                # active[skey].append(os.path.dirname(s))
                # add filepath
                active[skey].append(s)
            except KeyError as e:
                if debug:
                    raise e
                logger.error('%s : a table is missing', skey)
    except:
        if debug:
            raise 
        logger.error('file open error: %s', s)
    return active


def init_server(force_reload: bool = False, force_download: bool = False) -> tuple:
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
    active = None
    if os.path.isfile(active_file) and not force_reload:
        try:
            with open(active_file) as f:
                active = json.load(f)
            logger.info('Active Stations read. [%d]', len(active))
        except:
            active = None

    if active is None:
        #
        # find Merged Netcdf files and intercomparison files
        #
        slist = glob.glob(os.path.expandvars(config['data_dir'] + '/0-2000?-0-?????_CEUAS_merged_v0.nc'))
        slist += glob.glob(os.path.expandvars(config['data_dir'] + '/0-20?00-0-*_CEUAS_merged_v0.nc'))
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
        vola = pd.read_csv(f, delimiter='\t', quoting=3, dtype=tdict, na_filter=False)
        # print (vola.iloc[0])
        # exit()
        active = {}
        func = partial(makedaterange, vola)
        if False:
            with Pool(10) as p:
                sklist=list(p.map(func,zip(slist,slnum)))
        else:
            sklist = list(map(func, zip(slist, slnum)))
        
        for s in sklist:
            if s:
                k = next(iter(s))
                active[k] = s[k]
        logger.info('Active Stations created. [%d]', len(active))
        try:
            with open(active_file, 'w') as f:
                json.dump(active, f)
        except Exception as e:
            logger.warning('Cannot write %s: %s', active_file, e)
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


@hug.get_post('/status')
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
        # psutil.disk_usage('/data/public/')  # '/data/private/',
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
        flist: list fo files
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
           384                                                       # /data/public//tmp//006691463272/dest_0-20000-0-53513_relative_humidity.nc
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
        ds = xarray.open_dataset(fn)            
        
        to_be_removed = ['trajectory_index', 'trajectory']
        for ivar in list(ds.variables):
            if 'string' in ivar:
                to_be_removed.append(ivar)
            
            if 'trajectory' in ds[ivar].dims and ivar not in list(ds.coords):
                report_id = ds[ivar].astype(object).sum(axis=1).astype(str)
                ds = ds.drop_vars(ivar)
                ds[ivar] = ('obs', report_id.values[ds.trajectory_index.values])  # todo obs ???
                
            if ds[ivar].ndim > 1:
                tmp= ds[ivar].astype(object).sum(axis=1).astype(str)
                ds = ds.drop_vars(ivar)
                idim = tmp.dims[0]
                ds[ivar] = (idim, tmp)
        
        ds = ds.drop_vars(to_be_removed)
        df = ds.to_dataframe()
        #
        # todo fix the primary_id in the NetCDF files
        #
        if 'primary_id' not in ds.attrs:
            # /data/public//tmp//006691463272/dest_0-20000-0-53513_relative_humidity.nc
            df['statid'] = fn.split('/')[-1].split('_')[1]
            # logger.warning('CSV no primary_id in %s', fn)
            # continue
        else:
            df['statid'] = ds.attrs['primary_id']
        #
        #
        #
        df['statindex'] = statindex
        dfs.append(df)
        statindex += 1

    df = pd.concat(dfs, ignore_index=True)
    df.index.name = 'obs_id'
    #if '.zip' in ofile:
    #    df.to_csv(ofile, compression=dict(method='zip', archive_name=name + '.csv'), mode='a')  # might be 10x faster
    #else:
    df.to_csv(ofile)
    return ofile


###############################################################################


def check_body(variable: list = None, statid: list = None, product_type: str = None, pressure_level: list = None,
               day: list = None, month: list = None, year: list = None, date: list = None, time: list = None, 
               bbox: list = None, country: str = None, area: list = None,
               format: str = None, period: list = None, optional: list = None, wmotable: dict = None,
               pass_unknown_keys: bool = False,
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
    allowed_variables = ['temperature', 'u_component_of_wind', 'v_component_of_wind',
                         'wind_speed', 'wind_direction', 'relative_humidity',
                         'specific_humidity', 'dew_point_temperature']
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
    # Variable
    #
    if variable is None:
        raise KeyError('Missing argument: variable')
    else:
        if not isinstance(variable, list):
            variable = [variable]

        for ivar in variable:
            if ivar not in allowed_variables:
                raise KeyError('Invalid variable selected: ' + ivar)
        d['variable'] = variable
    #
    # Optional
    #
    allowed_optionals = ['sonde_type', 'bias_estimate','obs_minus_an','obs_minus_bg']
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
    # Check values:
    # ['obs_minus_an', 'obs_minus_bg', 'bias_estimate']
    
    #
    # Format
    #
    if format is not None:
        if format not in ['nc', 'csv']:
            raise ValueError('Invalid format selected [nc, csv]: ' + format)
        d['format'] = format
    else:
        d['format'] = 'nc'
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
                if vv[4] == icountry:
                    statid.append(k)

        if len(statid) == 0:
            raise RuntimeError('Invalid selection, no Stations for Countries %s' % str(country))
    #
    # BBOX [lower left upper right]
    #
    elif bbox is not None:
        # converting from BBOX [lower left upper right] to BBOX [upper left lower right]:
        bbox = [bbox[2], bbox[1], bbox[0], bbox[3]]
        if not isinstance(bbox, (list, tuple)) or len(bbox) != 4:
            raise ValueError('Invalid selection, bounding box: [upper left lower right]')

        try:
            for i in range(4):
                bbox[i] = float(bbox[i])
        except ValueError:
            raise ValueError('Invalid selection, bounding box: [upper left lower right] must be int or float')

        if bbox[0] >= bbox[2] or bbox[1] >= bbox[3]:
            raise ValueError('Invalid selection, bounding box: lower<upper [-90, 90], left<right [-180, 360]')

        if bbox[0] < -90 or bbox[0] > 90 or bbox[2] < -90 or bbox[2] > 90 or \
                bbox[1] < -180 or bbox[1] > 360 or bbox[3] < -180 or bbox[3] > 360 \
                or bbox[3] - bbox[1] > 360:
            raise ValueError('Invalid selection, bounding box: lower<upper [-90, 90], left<right [-180, 360]')
        statid = []
        for k, v in active.items():
            if bbox[0] <= v[2] <= bbox[2]:
                if bbox[3] <= 180:
                    if bbox[1] <= v[3] <= bbox[3]:
                        statid.append(k)
                else:
                    # rectangle crossing date line
                    if v[3] < 0:
                        if v[3] >= bbox[1] - 360 and v[3] + 360 <= bbox[3]:
                            statid.append(k)
                    else:
                        if bbox[1] <= v[3] <= bbox[3]:
                            statid.append(k)
        if len(statid) == 0:
            raise RuntimeError('Invalid selection, bounding box %s contains no radiosonde stations' % str(bbox))
    #
    # Stations
    #
    else:
        # '0-20200-0-*': 
        try:
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

                    if statid[:3] == '0-2' and statid in slnum:
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
                        if isinstance(valid_id, list):
                            new_statid = new_statid.extend(valid_id)
                        else:
                            new_statid.append(valid_id)
                            
                if valid_id == None:
                    raise ValueError('statid not available - please select an area, country or check your statid')
                    
                statid = [] 
                [statid.append(x) for x in new_statid if x not in statid] 
        except MemoryError:
            raise RuntimeError(
                'Invalid selection, specify either bbox, country or statid. Use "statid":"all" to select all ' \
                'stations ')
    d['statid'] = statid
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
        datelist = []
        newdates = []
        if isinstance(year, (int, str)):
            year = [year]
        if month is None:
            month = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
        elif isinstance(month, (int, str)):
            month = [month]
        if day is not None:
            day = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12',
                   '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24',
                   '25', '26', '27', '28', '29', '30', '31']
        elif isinstance(day, (int, str)):
            day = [day]
        # for removal of e.g. 31.02.:
        for i in year:
            for j in month:
                for k in day:
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
    if time is not None:
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
        d['time']=time

    return d


###############################################################################


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
    # todo this fucntion is wired ... redesign needed ?
    for b in body[spv[l]]:
        if l < len(spv) - 1:
            makebodies(bodies, body, spv, copy.copy(bo) + [b], l + 1)
        else:
            bodies.append(dict(body))
            bn = copy.copy(bo) + [b]
            for s, b in zip(spv, bn):
                bodies[-1][s] = b
                logger.debug('makebodies %d %s %s', l, s, b)
    return


def process_request(body: dict, output_dir: str, wmotable: dict, debug: bool = False) -> str:
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
    lenprod = (len(body['variable']) * len(body['statid']) * # variables and stations
               ((int(body['date'][-1][:4])-int(body['date'][0][:4]))*12 # years in months
                + (int(body['date'][-1][4:6])-int(body['date'][0][4:6])))) # months
#     if len(body['variable']) == 1 and ((int(body['date'][-1][:4])-int(body['date'][0][:4]))*12 + (int(body['date'][-1][4:6])-int(body['date'][0][4:6]))) == 1:
#         logger.warning('Requesting more than 500 elements - Exception: 1 variable and 1 month of every station')
    if lenprod > 30000:
        # lenght restriction deactivated as long following line is out commented.
        raise RuntimeError('Request too large - please split')
#         logger.warning('Request very large - please split')
    #
    logger.debug('Cleaned Request %s', str(body))
    os.makedirs(output_dir, exist_ok=True)  # double check
    bodies = []
    spv = ['statid', 'variable']  # potential split up variables
    bo = []
    refdate = datetime(year=1900, month=1, day=1)
    makebodies(bodies, body, spv, bo, 0)  # List of split requests
    #
    # Check all requests if dates are within active station limits
    #
    for k in range(len(bodies) - 1, -1, -1):
        # date selection ? do all the stations have data in this period?
        if 'date' in bodies[k].keys():
            #
            # seconds since Reference date
            #
            start = (to_valid_datetime(bodies[k]['date'][0]) - refdate).days * 86400
            ende = (to_valid_datetime(bodies[k]['date'][-1]) - refdate).days * 86400
            #
            # Station Active ?
            #
            if bodies[k]['statid'] in active.keys():
                if start > active[bodies[k]['statid']][1] or ende + 86399 < active[bodies[k]['statid']][0]:
                    logger.debug('%s outside Index range', bodies[k]['statid'])
                    del bodies[k]

                else:
                    logger.debug('%s Index[%d (%d / %d) %d]', bodies[k]['statid'],
                                 active[bodies[k]['statid']][0],
                                 start, ende,
                                 active[bodies[k]['statid']][1])
                    # add data path to request
                    # input_dirs.append(active[bodies[k]['statid']][5])  # path from makedaterange (init_server)
                    # file_paths.insert(0, active[bodies[k]['statid']][5])  # path from makedaterange (init_server)
                    bodies[k]['filename'] = active[bodies[k]['statid']][5]
            else:
                del bodies[k]
        else:
            # input_dirs.append(active[bodies[k]['statid']][5])  # path from makedaterange (init_server)
            bodies[k]['filename'] = active[bodies[k]['statid']][5]
            # file_paths.insert(0, active[bodies[k]['statid']][5])  # path from makedaterange (init_server)

    logger.debug('# requests %d', len(bodies))
    if len(bodies) == 0:
        raise RuntimeError('No selected station has data in specified date range: ' + str(body))

    # Make process_flat a function of only request_variables (dict)
    #
    # process_flat(outputdir: str, cftable: dict, debug:bool=False, request_variables: dict) -> tuple:
    # func = partial(eua.process_flat, output_dir, cf, input_dirs[0], debug)
    func = partial(eua.process_flat, output_dir, cf, debug)
    #
    # Smaller request?
    #
    if debug or len(body['variable']) * len(body['statid'])<10:
        #
        # Single Threading
        #
        results = list(map(func, bodies))
    else:
        #
        # Multi Threading
        #
        with Pool(10) as p:
            # error with chunksize (from p.map to p.starmap)
            results = p.map(func, bodies)
            # results = list(p.starmap(func, zip(input_dirs, [debug]*len(bodies), bodies), chunksize=1))
    #
    # Process the output 
    # todo catch Error Messages and store in a log file?
    #
    wpath = ''  # same as output_dir ?
    for r in results:
        if r[0] != '':
            wpath = r[0]
            break

    if wpath == '':
        raise RuntimeError('Error: %s (%s)' % (results[0][1], str(body)))
    else:
        rfile = os.path.dirname(wpath) + '/download.zip'

    logger.debug('wpath: %s; format %s Time %f', wpath, body['format'],time.time()-tt)

    if 'local_execution' in body.keys():
        return rfile

    if body['format'] == 'nc':
        with zipfile.ZipFile(rfile, 'w') as f:
            for r in results:
                try:
                    if len(r[0]) > 0:
                        f.write(r[0], os.path.basename(r[0]))
                    if debug:
                        continue  # do not remove
                    os.remove(r[0])  # remove NetCDF file
                except:
                    pass
        logger.debug('netcdfs compressed [%d] to %s', len(results), rfile)

    else:
        with zipfile.ZipFile(rfile, 'w', compression=zipfile.ZIP_DEFLATED) as f:
            for v in body['variable']:
                ilist = glob.glob(output_dir + '/*' + v + '.nc')
                if len(ilist) > 0:
                    ifile = to_csv(ilist, ofile=output_dir + '/' + v + '.csv')  # todo add correct name into zip
                    f.write(ifile, os.path.basename(ifile))
                    logger.debug('writing csv %s [%d] to %s', v, len(ilist), rfile)
                    if debug:
                        continue  # do not remove
                    for i in ilist:
                        os.remove(i)  # remove NetCDF file
                    os.remove(ifile)  # remove csv

    logger.debug('Request-File: %s [Time: %7.4f s]', rfile, (time.time() - tt))
    return rfile


@hug.get('/', output=hug.output_format.file)
def index(request=None, response=None):
    """ Main Hug Index Function on get requests

    index function requests get URI and converts into dictionary.

    Args:
        request: dictionary
        response: str

    Returns:

    """
    logger.debug("GET %s", request.query_string)
    if '=' not in request.query_string:
        response.status = hug.HTTP_422
        raise Exception('A query string must be supplied')

    try:
        rs = request.query_string.split('&')
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
        response.status = hug.HTTP_422
        raise Exception(request.query_string)

    randdir = '{:012d}'.format(numpy.random.randint(100000000000))
    logger.info("%s GET %s", randdir, str(body))
    tmpdir = config['tmp_dir'] + '/' + randdir
    try:
        # rfile = process_request(body, tmpdir, config['data_dir'], wmo_regions)
        rfile = process_request(body, tmpdir, wmo_regions, debug=config['debug'])
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

    response.set_header('Content-Disposition', 'attachment; filename=' + os.path.basename(rfile))
    return rfile


@hug.exception(Exception)
def base_exception_handler(exception, response=None):
    """ This captures any Exception from the Server

    Args:
        exception: hug.exceptions.Exception class
        response: HTTP Response

    Returns:
        str :  Message
    """
    response.status = hug.HTTP_422  # Always the same Error ? or dependent on exception?
    # one of the arguments could be the response status
    return ",".join(exception.args)


@hug.post('/', output=hug.output_format.file)
def index(request=None, body=None, response=None):
    """ Main Hug index function for Post requests

    Args:
        request: not used
        body: dictionary of request
        response:

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
    randdir = '{:012d}'.format(numpy.random.randint(100000000000))
    logger.info("%s POST %s", randdir, str(body))
    tmpdir = config['tmp_dir'] + '/' + randdir
    try:
        #rfile='/fio/srvx7/leo/x'
        rfile = process_request(body, tmpdir, wmo_regions, debug=config['debug'])
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

    response.set_header('Content-Disposition', 'attachment; filename=' + os.path.basename(rfile))
    return rfile


# @hug.get('/dataset')
# def dataset():
#     from pydap.wsgi.app import DapServer
#     return DapServer('/data/public/tmp/100000000000')  # maybe ??

# @hug.get('/dataset', output=hug.output_format.file)
# def opendap(request=None, response=None):
#     # todo does not work with pydap
#     #  it seems that
#     # from pydap.wsgi.app import DapServer
#     # application = DapServer('/data/public/tmp/100000000000')  # maybe ??
#     # should we return the DapServer ?
#     rfile = '/data/public/tmp/100000000000/dest_0-20000-0-70398_air_temperature.nc'
#     response.set_header('Content-Disposition', 'attachment; filename=' + os.path.basename(rfile))
#     return rfile

@hug.get('/constraints/', output=hug.output_format.file)
def index(request=None, response=None):
    """ Main Hug Index Function on get requests

    index function requests get URI and converts into dictionary.

    Args:
        request: dictionary
        response: str

    Returns:

    """
    logger.debug("GET %s", request.query_string)


    rfile='/data/public/constraints.csv'

    response.set_header('Content-Disposition', 'attachment; filename=' + os.path.basename(rfile))
    return rfile


if __name__ == '__main__':
    active, wmo_regions, cf = init_server()
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
    ret = process_request(body, tmpdir, wmo_regions, debug=debug)
    logger.debug(str(ret))
