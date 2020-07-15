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
import logging
import os
import socket
import sys
import time
import zipfile
from typing import Union
from datetime import datetime, timedelta
from functools import partial
from multiprocessing import set_start_method, Pool

import cds_eua2 as eua
import h5py
import hug
import numpy
import pandas as pd
import xarray
from falcon import HTTPError, HTTP_422

try:
    set_start_method("spawn")
except RuntimeError:
    pass

"""
Logging by external:
STDOUT > hug.log
STDERR > hug.err

hug.log (INFO)
hug.err (ERROR)
hug.debug.log (DEBUG)
finished_requests.log : Add successfully delivered requests
failed_requests.log : Add failed requests
"""


class MyFilter(object):
    def __init__(self, level):
        self.__level = level

    def filter(self, logrecord):
        return logrecord.levelno <= self.__level


logger = logging.getLogger('upperair')
logger.setLevel(logging.DEBUG)
# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s | %(funcName)s - %(levelname)s - %(message)s')
# create console handler and set level to warning for stderr
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
    ch3 = logging.FileHandler('logs/hug.debug.log')
except PermissionError:
    ch3 = logging.FileHandler('logs/hug.debug.local.log')
ch3.setLevel(logging.DEBUG)  # respond only to Debug and above
ch3.setFormatter(formatter)
logger.addHandler(ch3)

host = socket.gethostname()

logger.info("HUG started on %s", host)
if 'srvx' in host:
    sys.path.append(os.path.expanduser('~leo/python/'))
else:
    sys.path.append('/data/private/soft/python/')
    os.environ["RSCRATCH"] = "/data/public/"

global wroot
wroot = '.'


def makedaterange(vola: pd.DataFrame, itup: tuple) -> dict:
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
                active[skey] = [int(eua.secsince(f['recordtimestamp'][0], funits)),
                                int(eua.secsince(f['recordtimestamp'][-1], funits)),
                                float(f['observations_table']['latitude'][-1]),
                                float(f['observations_table']['longitude'][-1])]
                idx = numpy.where(vola.StationId.values == skey)[0]
                if len(idx) > 0:
                    active[skey].append(vola.CountryCode[idx[0]])
                else:
                    active[skey].append('')
                    logger.debug('no key found for %s', skey)
            except KeyError:
                logger.error('%s : a table is missing', skey)
    except:
        logger.error('file open error: %s', s)
    return active


def init_server() -> tuple:
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
    # todo make filepath more global and not defined in some function !!!
    os.chdir(os.path.expandvars('$RSCRATCH/era5/odbs/merged'))
    wroot = os.path.expandvars('$RSCRATCH/era5/odbs/merged/tmp')
    os.makedirs(wroot, exist_ok=True)
    z = numpy.zeros(1, dtype=numpy.int32)
    zidx = numpy.zeros(1, dtype=numpy.int)
    idx = numpy.zeros(1, dtype=numpy.int)
    trajectory_index = numpy.zeros(1, dtype=numpy.int)
    # What get's updated here?
    zz = eua.calc_trajindexfast(z, zidx, idx, trajectory_index)
    os.makedirs(os.path.expanduser('~/.tmp'), exist_ok=True)
    active_file = os.path.expanduser('~/.tmp/active.json')
    active = None
    if os.path.isfile(active_file):
        try:
            with open(active_file) as f:
                active = json.load(f)
            logger.info('Active Stations read. [%d]', len(active))
        except:
            active = None

    if active is None:
        slist = glob.glob(os.path.expandvars('0-20000-0-?????_CEUAS_merged_v0.nc'))
        slnum = [i[-34:-19] for i in slist]
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
        # with Pool(10) as p:
        # sklist=list(p.map(func,zip(slist,slnum)))
        sklist = list(map(func, zip(slist, slnum)))
        for s in sklist:
            if s:
                k = next(iter(s))
                active[k] = s[k]
        logger.info('Active Stations created. [%d]', len(active))
        with open(active_file, 'w') as f:
            json.dump(active, f)
    #
    # Read CDM Definitions
    #
    cdm_file = os.path.expanduser('~/.tmp/cf.json')
    if os.path.isfile(cdm_file):
        with open(cdm_file) as f:
            cf = json.load(f)
    else:
        cf = eua.read_standardnames()
        with open(cdm_file, 'w') as f:
            json.dump(cf, f)

    # cdmtablelist=['id_scheme','crs','station_type','observed_variable','station_configuration_codes','units','sub_region']
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
# CDM Table
# CF Table Naming
###############################################################################


active, cdm, cf = init_server()

# Active Station Numbers
slnum = list(active.keys())
slist = [os.path.expandvars('$RSCRATCH/era5/odbs/merged/') + '0-20000-0-' + s + '_CEUAS_merged_v0.nc' for s in slnum]


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


def to_csv(flist: list, ofile: str = 'out.csv'):
    """ Convert every file in flist to CSV

    Args:
        flist: list fo files
        ofile: output filename

    Returns:
        str: output filename (returned by the request)
    """
    statindex = 0
    dfs = []
    for fn in flist:
        ds = xarray.open_dataset(fn, drop_variables=['trajectory_label'])
        df = ds.to_dataframe()
        df['statid'] = ds.attrs['primary_id']
        df['statindex'] = statindex
        dfs.append(df)
        statindex += 1

    df = pd.concat(dfs, ignore_index=True)
    df.index.name = 'obs_id'
    df.to_csv(ofile)
    return ofile


###############################################################################


def check_body(variable: list = None, statid: list = None, product_type: str = None, pressure_level: list = None,
               date: list = None, time: list = None, fbstats=None, bbox: list = None, country: str = None,
               format: str = None, period: list = None, cdm: dict = None, pass_unknown_keys: bool = False,
               **kwargs) -> Union[dict, str]:
    """ Check Request for valid values and keys

    Args:
        variable: e.g. temperature, ...
        statid: e.g. 01001
        product_type: currently unused
        pressure_level: 500, ...
        date: '19990131' or ['19990101', '20000101'] as range
        time: 1,... or 0-23
        fbstats: only these are currently allowed: 'obs_minus_an', 'obs_minus_bg', 'bias_estimate'
        bbox: Bounding Box [lower left upper right]
        country: Country Code, e.g. DEU, AUT, USA, GBR
        format: nc or csv
        period: ['19990101', '20000101'] see Notes
        cdm: CDM definition Table
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
        return 'Missing argument: variable'
    else:
        if not isinstance(variable, list):
            variable = [variable]

        for ivar in variable:
            if ivar not in allowed_variables:
                return 'Invalid variable selected: ' + ivar
        d['variable'] = variable
    #
    # fb stats
    #
    if fbstats is not None:
        if not isinstance(fbstats, list):
            fbstats = [fbstats]
        for ifb in fbstats:
            if ifb not in ['obs_minus_an', 'obs_minus_bg', 'bias_estimate']:
                return 'Invalid fbstats variable selected: ' + ifb
    #
    # Format
    #
    if format is not None:
        if format not in ['nc', 'csv']:
            return 'Invalud format selected [nc, csv]: ' + format
        d['format'] = format
    else:
        d['format'] = 'nc'
    #
    # only one of [statid, bbox, country]
    #
    if statid is None and bbox is None and country is None:
        statid = 'all'

    if sum((statid is not None, bbox is not None, country is not None)) != 1:
        return 'Invalid selection, Specify only one of statid: %s, bbox: %s and country: %s' % (
            statid, bbox, country)
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
        vcountries = cdm['sub_region'].alpha_3_code.values
        for icountry in country:
            if icountry not in vcountries:
                return 'Invalid selection, %s is not a valid country code' % icountry

            for k, vv in active.items():
                if vv[4] == icountry:
                    statid.append(k)

        if len(statid) == 0:
            return 'Invalid selection, no Stations for Countries %s' % str(country)
    #
    # BBOX [lower left upper right]
    #
    elif bbox is not None:
        if not isinstance(bbox, (list, tuple)) or len(bbox) != 4:
            return 'Invalid selection, bounding box: [lower, left, upper, right]'

        try:
            for i in range(4):
                bbox[i] = float(bbox[i])
        except ValueError:
            return 'Invalid selection, bounding box: [lower, left, upper, right] must be int or float'

        if bbox[0] > bbox[2] or bbox[1] > bbox[3]:
            return 'Invalid selection, bounding box: lower<upper [-90, 90], left<right [-180, 360]'

        if bbox[0] < -90 or bbox[0] > 90 or bbox[2] < -90 or bbox[2] > 90 or \
                bbox[1] < -180 or bbox[1] > 360 or bbox[3] < -180 or bbox[3] > 360 \
                or bbox[3] - bbox[1] > 360:
            return 'Invalid selection, bounding box: lower<upper [-90, 90], left<right [-180, 360]'
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
            return 'Invalid selection, bounding box %s contains no radiosonde stations' % str(bbox)
    #
    # Stations
    #
    else:
        try:
            if statid == 'all':
                statid = slnum

            elif isinstance(statid, (str, int)):
                for s in ['0-20000-0-', '0-20001-0-']:
                    if isinstance(statid, int):
                        valid_id = s + '{:0>5}'.format(statid)
                    else:
                        if statid[:3] == '0-2':
                            valid_id = statid
                            break

                        valid_id = s + statid
                    if valid_id in slnum:
                        break
                statid = [valid_id]
            else:
                new_statid = []
                for k in range(len(statid)):
                    for s in ['0-20000-0-', '0-20001-0-']:
                        if isinstance(k, int):
                            valid_id = s + '{:0>5}'.format(k)
                        else:
                            if k[:3] == '0-2':
                                valid_id = k
                                break

                            valid_id = s + k
                        if valid_id in slnum:
                            break

                    new_statid.append(valid_id)
                statid = new_statid
        except MemoryError:
            return 'Invalid selection, specify either bbox, country or statid. Use "statid":"all" to select all ' \
                   'stations '
    d['statid'] = statid
    #
    # Date time selection
    # [DATE] or [START, END]
    #
    # todo only one date or two dates allowed at the moment by CDS
    if date is not None:
        # str, list (str, int)
        newdates = []
        # make a list
        if isinstance(date, (int, str)):
            date = [date]

        for idate in date:
            # convert to string
            if not isinstance(idate, str):
                idate = str(idate)
            # todo old style date range (not supported by CDS anymore ???)
            if '-' in idate:
                idate = idate.split('-')
                # check period dates (should not be out of range)
                newdates.append('{}-{}'.format(to_valid_datetime(idate[0], as_string=True),
                                               to_valid_datetime(idate[-1], as_string=True)))
            else:
                try:
                    newdates.append(to_valid_datetime(idate, as_string=True))
                except:
                    return 'only valid dates allowed for date: %s' % idate
        d['date'] = newdates
    #
    # Period [START, END] -> into date
    #
    # todo not forward by CDS -> to date [start end]
    if period is not None:
        if not isinstance(period, list):
            return 'invalid period selection, period [startdate, enddate], but %s' % str(period)

        for i in range(len(period)):
            period[i] = str(period[i])
        d['date'] = [to_valid_datetime(period[0], as_string=True), to_valid_datetime(period[-1], as_string=True)]
    #
    # Pressure levels
    #
    if pressure_level is not None:
        if not isinstance(pressure_level, list):
            pressure_level = [pressure_level]

        for i in range(len(pressure_level)):
            try:
                # need to be string
                pressure_level[i] = str(pressure_level[i])
                # in Pa
                if int(pressure_level[i]) < 500 or int(pressure_level[i]) > 110000:
                    return 'invalid selection, pressure_level out of range [50-1100 hPa]: %d' % int(
                        pressure_level[i]) / 100

            except:
                return 'invalid selection, pressure_level allows only integer, ' + pressure_level[i]
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
                    return 'invalid selection, time out of range [0-24 h]: %d' % int(time[i])
            except:
                return 'invalid selection, pressure_level allows only integer, ' + time[i]

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
    # todo this fucntion is wired ... redesign needed
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


def defproc(body: dict, wroot: str, randdir: str, cdm: dict, debug:bool=False) -> tuple:
    """ Main fucntion of the hug server

    Args:
        body: request dictionary
        wroot: station path root
        randdir: temporary directory for request
        cdm: CDM table definitions from init_server()
        debug: for debugging

    Returns:
        str : filename of zipped requested files
        str : message or error
    """
    tt = time.time()
    msg = check_body(cdm=cdm, **body)
    if isinstance(msg, str):
        return '', msg
    body = msg
    logger.debug('Cleaned Request %s', str(body))
    os.makedirs(wroot + '/' + randdir, exist_ok=True)
    bodies = []
    spv = ['statid', 'variable']  # potential split up variables
    bo = []
    refdate = datetime(year=1900, month=1, day=1)
    makebodies(bodies, body, spv, bo, 0)  # List of split requests
    #
    # Check all requests if dates are within active station limits
    #
    for k in range(len(bodies) - 1, -1, -1):
        deleted = False  # delete sonde from request ?
        if 'date' in bodies[k].keys():
            # seconds since Reference date
            start = (to_valid_datetime(bodies[k]['date'][0]) - refdate).days * 86400
            ende = (to_valid_datetime(bodies[k]['date'][-1]) - refdate).days * 86400
            # Station Active ?
            if bodies[k]['statid'] in active.keys():
                if start > active[bodies[k]['statid']][1] or ende + 86399 < active[bodies[k]['statid']][0]:
                    logger.debug('%s outside Index range', bodies[k]['statid'])
                    del bodies[k]
                    deleted = True

                else:
                    logger.debug('%s Index[%d (%d / %d) %d]', bodies[k]['statid'],
                                 active[bodies[k]['statid']][0],
                                 start, ende,
                                 active[bodies[k]['statid']][1])
            else:
                del bodies[k]
                deleted = True

    logger.debug('# requests %d', len(bodies))
    if len(bodies) == 0:
        return '', 'No selected station has data in specified date range: ' + str(body)

    func = partial(eua.process_flat, wroot, randdir, cf)

    if debug:
        #
        # Single Threading
        #
        results = list(map(func, bodies))
    else:
        #
        # Multi Threading
        #
        with Pool(10) as p:
            results = list(p.map(func, bodies, chunksize=1))

    wpath = ''
    for r in results:
        if r[0] != '':
            wpath = r[0]
            break

    if wpath == '':
        return '', 'Error: %s (%s)' % (results[0][1], str(body))
    else:
        rfile = os.path.dirname(wpath) + '/download.zip'

    logger.debug('wpath: %s; format %s', wpath, body['format'])

    if 'local_execution' in body.keys():
        return rfile, ''

    if body['format'] == 'nc':
        with zipfile.ZipFile(rfile, 'w') as f:
            for r in results:
                try:
                    # logger.debug('Zipping: %s in %s', r[0], os.path.basename(r[0]))
                    if len(r[0]) > 0:
                        f.write(r[0], os.path.basename(r[0]))
                    os.remove(r[0])
                except:
                    pass

    else:
        ofiles = []
        for v in body['variable']:
            ilist = glob.glob(os.path.dirname(wpath) + '/*' + v + '.nc')
            if len(ilist) > 0:
                ofiles.append(to_csv(ilist, v + '.csv'))
                logger.debug('writing csv %s [%d] to %s', v, len(ilist), rfile)
                for i in ilist:
                    os.remove(i)

        if len(ofiles) > 0:
            with zipfile.ZipFile(rfile, 'w', compression=zipfile.ZIP_DEFLATED) as f:
                for o in ofiles:
                    try:
                        f.write(o, os.path.basename(o))
                        logger.debug('writing %s to %s', o, rfile)
                        os.remove(o)
                    except:
                        pass

    logger.debug('Request-File: %s [Time: %7.4f s]', rfile, (time.time() - tt))
    return rfile, ''


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
        raise HTTPError(HTTP_422, title='malformed get request', description='A query string must be supplied')

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
        raise HTTPError(HTTP_422, title='malformed query string', description=request.query_string)

    randdir = '{:012d}'.format(numpy.random.randint(100000000000))
    logger.debug("GET BODY %s", str(body))
    wroot = os.path.expandvars('$RSCRATCH/tmp/')
    rfile, error = defproc(body, wroot, randdir, cdm)
    if rfile == '':
        logger.error("GET Request failed, %s", error)
        raise HTTPError(HTTP_422, title='malformed request', description=error)

    response.set_header('Content-Disposition', 'attachment; filename=' + os.path.basename(rfile))
    return rfile


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
    # todo add trigger for reloading active stations and restarting hug
    # todo add status request including uptime, ...
    randdir = '{:012d}'.format(numpy.random.randint(100000000000))
    request_id = id(body)  # or somthing better?
    logger.info("%s POST %s", request_id, str(body))
    wroot = os.path.expandvars('$RSCRATCH/tmp/')
    rfile, error = defproc(body, wroot, randdir, cdm)
    if rfile == '':
        logger.error("POST Request failed, %s", error)
        with open('./logs/failed_requests.log', 'a+') as ff:
            ff.write('%s - %s [%s] Message: %s \n' % (str(datetime.now()), request_id, str(body), error))

        logger.info("%s POST FAILED %s", request_id, error)
        raise HTTPError(HTTP_422, title='malformed request', description=error)
    logger.info("%s POST FINISHED", request_id)
    #
    # Write successful requests
    #
    with open('./logs/finished_requests.log', 'a+') as ff:
        ff.write('%s - %s [%s] \n' % (str(datetime.now()), request_id, str(body)))

    response.set_header('Content-Disposition', 'attachment; filename=' + os.path.basename(rfile))
    return rfile


if __name__ == '__main__':
    active, cdm, cf = init_server()
    #
    # Parse command line arguments for testing the server API
    #
    body = eval(sys.argv[1])
    debug = body.pop('debug', False)
    #
    # Logging to DEBUG to std.out and hug.debug.local.log
    #
    logger.setLevel(10)
    for i in logger.handlers:
        i.setLevel(10)
    #
    # Specific directory for testing
    #
    randdir = os.path.expandvars('{:012d}'.format(100000000000))
    if os.path.isdir(randdir):
        for ifile in os.listdir(randdir):
            print(randdir + '/' + ifile)
            os.remove(randdir + '/' + ifile)

    wroot = os.path.expandvars('$RSCRATCH/tmp/')
    logger.debug(str(body))
    #
    # Run the request
    #
    ret = defproc(body, wroot, randdir, cdm)
    idir = os.path.expandvars('$RSCRATCH/out')
    os.chdir(idir)
    logger.debug(str(ret))
