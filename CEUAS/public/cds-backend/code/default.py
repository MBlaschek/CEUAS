# -*- coding: utf-8 -*-
# this file is released under public domain and you can use without limitations

#########################################################################
# This is a simple controller
# - index is the default action of any application
#   here it forwards the request from a CDS server to further processing
#   with python
#
#   For queries ask early-upper-air@copernicus-climate.eu
#   The C3S 311c Lot2 group
#   Vienna, 26 August 2019
# - user is required for authentication and authorization
# - download is for downloading files uploaded in the db (does streaming)
# - api is an example of Hypermedia API support and access control
#########################################################################

import os
import socket
import sys
import urllib
import cds_eua2 as eua
import pandas as pd
import xarray
import numpy
import hug
import h5py  # pickle as h5py
import zipfile
import json
import glob
from functools import partial
from falcon import HTTPError, HTTP_422
import copy
import time
from datetime import datetime, timedelta
import logging
from multiprocessing import set_start_method, Pool

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

    def filter(self, logRecord):
        return logRecord.levelno <= self.__level


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

# @hug.exception(Exception)
# def handler(exception):
# return str(exception)
global wroot
wroot = '.'


def makedaterange(vola: pd.DataFrame, itup: tuple):
    s, skey = itup
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


def main():
    #    try:
    #        set_start_method("spawn")
    #    except RuntimeError:
    #        pass
    os.chdir(os.path.expandvars('$RSCRATCH/era5/odbs/merged'))
    #    os.chdir(os.path.expandvars('/raid60/scratch/federico/JUNE_TEST_MERGING_ALL'))
    wroot = os.path.expandvars('$RSCRATCH/era5/odbs/merged/tmp')
    try:
        os.mkdir(wroot)
    except:
        pass
    z = numpy.zeros(1, dtype=numpy.int32)
    zidx = numpy.zeros(1, dtype=numpy.int)
    idx = numpy.zeros(1, dtype=numpy.int)
    trajectory_index = numpy.zeros(1, dtype=numpy.int)
    zz = eua.calc_trajindexfast(z, zidx, idx, trajectory_index)
    try:
        os.mkdir(os.path.expanduser('~/.tmp'))
    except:
        pass

    try:
        with open(os.path.expanduser('~/.tmp/active.json')) as f:
            active = json.load(f)
    except:
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

        with open(os.path.expanduser('~/.tmp/active.json'), 'w') as f:
            json.dump(active, f)

    cf = eua.read_standardnames()

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


active, cdm, cf = main()

slnum = list(active.keys())
slist = [os.path.expandvars('$RSCRATCH/era5/odbs/merged/') + '0-20000-0-' + s + '_CEUAS_merged_v0.nc' for s in slnum]


# filelist=glob.glob('chera5.conv._?????.nc')
def last_day_of_month(any_day):
    next_month = any_day.replace(day=28) + timedelta(days=4)  # this will never fail
    return next_month - timedelta(days=next_month.day)


def to_valid_datetime(idate, as_string=False):
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


def check_body_new(variable: list = None, statid: list = None, product_type: str = None, pressure_level: list = None,
                   date: list = None, time: list = None, fbstats=None, bbox: list = None, country: str = None,
                   format: str = None, period: list = None, cdm: dict = None, pass_unknown_keys: bool = False,
                   **kwargs) -> str or dict:
    d = {}

    for ikey, ival in kwargs.items():
        logger.info('Requested key %s : %s unknown', ikey, ival)
        if pass_unknown_keys:
            d[ikey] = ival

    if variable is None:
        return 'Missing argument: variable'

    if sum((statid is None, bbox is None, country is None)) != 1:
        return 'Invalid selection, Specify only one of statid: %s, bbox: %s and country: %s' % (
            statid, bbox, country)

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
            return 'Invalid selection, specify either bbox, country or statid. Use "statid":"all" to select all stations'
        d['statid'] = statid

    if date is not None:
        # str, list (str, int)
        newdates = []
        if isinstance(date, str):
            date = [date]

        for idate in date:
            # if not isinstance(idate, list):
            #     idate = str(idate)
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

    if period is not None:
        if not isinstance(period, list):
            return 'invalid period selection, period [startdate, enddate], but %s' % str(period)
        for i in range(len(period)):
            period[i] = str(period[i])
        d['date'] = ['{}-{}'.format(to_valid_datetime(period[0], as_string=True),
                                   to_valid_datetime(period[-1], as_string=True))]

    if product_type is not None:
        logger.warning('Not yet implemented: product_type : %s' % product_type)

    return d


def check_body(body: dict, cdm: dict):
    """ Check dictionary from request if required keys are present and within ranges

    Args:
        body: request body
        cdm: CDM definitions table (here used only for Country codes)

    Returns:
        list: title, description
        str: empty
    """
    required_keys = ['variable']
    valid_keys = required_keys + ['statid', 'fbstats', 'pressure_level', 'date', 'time', 'bbox', 'country', 'format',
                                  'period']
    xor_keys = ['statid', 'bbox', 'country']
    valid_ranges = {'pressure_level': ['500', '110000'],
                    'date': ['19000101', '20301231'], 'time': ['0', '24'],
                    'fbstats': ['obs_minus_an', 'obs_minus_bg', 'bias_estimate'],
                    'variable': ['temperature', 'u_component_of_wind', 'v_component_of_wind',
                                 'wind_speed', 'wind_direction', 'relative_humidity',
                                 'specific_humidity', 'dew_point_temperature'],
                    'format': ['nc', 'csv']}

    try:
        bk = list(body.keys())
        for r in required_keys:
            if r not in bk:
                return ['Missing argument:', 'Argument ' + r + ' is required']

        for b in bk:
            if b not in valid_keys:
                logger.warning('Invalid argument key: %s', b)
                continue
                # return ['Invalid argument ' + b + '.', ' Valid arguments:' + str(valid_keys)]

        rxor = ''
        for r in xor_keys:
            if r in bk:
                if len(rxor) == 0:
                    rxor = r
                else:
                    return ['Invalid station selection', 'Please do not specify both ' + rxor + ' and ' + r]

        if 'date' in bk and 'period' in bk:
            return ['Invalid datetime selection', 'Please do not specify both date and period']

        if 'country' in bk:
            if type(body['country']) is str:

                if body['country'].upper() in ('GLOBE', 'ALL'):
                    body['statid'] = slnum
                    del body['country']
                else:
                    body['country'] = [body['country']]

            if 'country' in body.keys():
                vcountries = cdm['sub_region'].alpha_3_code.values
                body['statid'] = []
                for v in body['country']:
                    if v not in vcountries:
                        return ['Invalid station selection', v + ' is not a valid country code']
                    for k, vv in active.items():
                        if vv[4] == v:
                            body['statid'].append(k)
                if len(body['statid']) == 0:
                    return ['Invalid station selection', 'Countries ' + str(body['country']) + ' have no radiosondes']
                del body['country']
        elif 'bbox' in bk:
            if type(body['bbox']) is not list:
                return ['Invalid station selection', 'Bounding box: [lower, left, upper, right]']
            if len(body['bbox']) != 4:
                return ['Invalid station selection', 'Bounding box: [lower, left, upper, right]']
            try:
                for i in range(4):
                    body['bbox'][i] = float(body['bbox'][i])
            except:
                return ['Invalid station selection', 'Bounding box: [lower, left, upper, right] must be int or float']
            if body['bbox'][0] > body['bbox'][2] or body['bbox'][1] > body['bbox'][3]:
                return ['Invalid station selection',
                        'Bounding box requirements: lower<upper, left<right, -90<=lat<=90, -180<=lon<=360']
            if body['bbox'][0] < -90 or body['bbox'][0] > 90 or body['bbox'][2] < -90 or body['bbox'][2] > 90 or \
                    body['bbox'][1] < -180 or body['bbox'][1] > 360 or body['bbox'][3] < -180 or body['bbox'][3] > 360 \
                    or body['bbox'][3] - body['bbox'][1] > 360:
                return ['Invalid station selection',
                        'Bounding box requirements: lower<upper, left<right, -90<=lat<=90, -180<=lon<=360']

            body['statid'] = []
            for k, v in active.items():
                if v[2] >= body['bbox'][0] and v[2] <= body['bbox'][2]:
                    if body['bbox'][3] <= 180:
                        if v[3] >= body['bbox'][1] and v[3] <= body['bbox'][3]:
                            body['statid'].append(k)
                    else:  # rectangle crossing date line
                        if v[3] < 0:
                            if v[3] >= body['bbox'][1] - 360 and v[3] + 360 <= body['bbox'][3]:
                                body['statid'].append(k)
                        else:
                            if v[3] >= body['bbox'][1] and v[3] <= body['bbox'][3]:
                                body['statid'].append(k)
            if len(body['statid']) == 0:
                return ['Invalid station selection',
                        'Bounding box ' + str(body['bbox']) + ' contains no radiosonde stations']
            del body['bbox']
        else:
            try:
                suff = ['0-20000-0-', '0-20001-0-']
                if type(body['statid']) == int:
                    for s in suff:
                        bd = s + '{:0>5}'.format(body['statid'])
                        if bd in slnum:
                            break
                    body['statid'] = [bd]
                elif body['statid'] == 'all':
                    body['statid'] = slnum
                elif body['statid'][0] == 'all':
                    body['statid'] = slnum
                elif type(body['statid']) is not list:
                    bd = body['statid']
                    if body['statid'][:3] != '0-2':
                        for s in suff:
                            bd = s + body['statid']
                            if bd in slnum:
                                break
                    body['statid'] = [bd]
                else:
                    if type(body['statid'][0]) == int:
                        for k in range(len(body['statid'])):
                            for s in suff:
                                bd = s + '{:0>5}'.format(body['statid'][k])
                                if bd in slnum:
                                    break
                            body['statid'][k] = bd
                    else:
                        for k in range(len(body['statid'])):
                            if body['statid'][k][:3] != '0-2':
                                for s in suff:
                                    bd = s + body['statid'][k]
                                    if bd in slnum:
                                        break
                                body['statid'][k] = bd

            except MemoryError:
                return ['Invalid station selection',
                        'Please specify either bbox, country or statid for station selection. Use "statid":"all" to select all stations']

        refdate = datetime(year=1900, month=1, day=1)
        bk = list(body.keys())
        for v in bk:
            if v in valid_ranges:
                #
                # Convert to list
                #
                if isinstance(body[v], (str, int)):
                    body[v] = str(body[v])
                    body[v] = [body[v]]
                #
                # Convert to String
                #
                for k in range(len(body[v])):
                    body[v][k] = str(body[v][k])

                #
                # Check Valid Ranges
                #
                for bv in body[v]:
                    if v in ('pressure_level', 'time'):
                        try:
                            if int(bv) > int(valid_ranges[v][1]) or int(bv) < int(valid_ranges[v][0]):
                                return ['argument value(s) ' + str(bvv) + ' not valid.',
                                        'Valid values:' + str(valid_ranges[v])]
                        except:
                            return ['only integer arguments allowed for ' + v, str(bvv) + ' not valid.']
                    elif v == 'date':
                        continue
                    else:
                        if bv not in valid_ranges[v]:
                            return ['argument value(s) ' + str(bv) + ' not valid.',
                                    'Valid values:' + str(valid_ranges[v])]

                if v == 'date':
                    newdates = []
                    singledates = []
                    for k in range(len(body[v])):
                        # period
                        if '-' in body[v][k]:
                            bvv = body[v][k].split('-')
                            # check period dates (should not be out of range)
                            newdates.append('{}-{}'.format(to_valid_datetime(bvv[0], as_string=True),
                                                           to_valid_datetime(bvv[-1], as_string=True)))
                        else:
                            # try:
                            #     if int(body[v][k]) > int(valid_ranges[v][1]) or int(body[v][k]) < int(
                            #             valid_ranges[v][0]):
                            #         return ['argument value(s) ' + str(body[v][k]) + ' not valid.',
                            #                 'Valid values:' + str(valid_ranges[v])]
                            # except:
                            #     return ['only integer arguments allowed for ' + v, str(body[v][k]) + ' not valid.']
                            #
                            try:
                                singledates.append(to_valid_datetime(body[v][k], as_string=True))
                            except:
                                return ['only valid dates allowed for ' + v, str(body[v][k])]
                    #
                    # Continuous block of dates ?
                    #
                    if len(singledates) > 0:
                        start = (to_valid_datetime(singledates[0]) - refdate).days
                        stop = (to_valid_datetime(singledates[-1]) - refdate).days

                        # check if number of days between start and stop match length -> continous block
                        if len(singledates) == stop - start + 1:
                            newdates.append('{}-{}'.format(to_valid_datetime(singledates[0], as_string=True),
                                                           to_valid_datetime(singledates[-1], as_string=True)))
                            logger.debug('Dates to period: %s', newdates[-1])
                        else:
                            newdates.extend(singledates)

                    body[v] = newdates  # list(set(sorted(newdates)))   # unique list of dates
                if v == 'period':
                    body['date'] = '{}-{}'.format(to_valid_datetime(body[v][0], as_string=True),
                                                  to_valid_datetime(body[v][-1], as_string=True))
                    del body[v]
                logger.debug('%s %s %s [%d]', v, body[v][0], body[v][-1], len(body[v]))

    except IOError:
        logger.error("General syntax error %s", str(body))
        return ['general syntax error ', body]

    return ''


def makebodies(bodies, body, spv, bo, l):
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


def to_csv(flist: list, ofile: str = 'out.csv'):
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


def defproc(body: dict, wroot: str, randdir: str, cdm: dict):
    tt = time.time()
    error = check_body(body, cdm)
    if False:
        msg = check_body_new(cdm=cdm, **body)
        if isinstance(msg, str):
            return '', msg
        body = msg
    if len(error) > 0:
        return '', " ".join(error)

    logger.debug('Cleaned Request %s', str(body))
    try:
        os.mkdir(wroot + '/' + randdir)
    except:
        pass

    bodies = []
    spv = ['statid', 'variable']
    bo = []
    refdate = datetime(year=1900, month=1, day=1)
    makebodies(bodies, body, spv, bo, 0)
    for k in range(len(bodies) - 1, -1, -1):
        deleted = False
        if 'date' in bodies[k].keys():
            # if bodies[k]['date'][0]>active[bodies[k]['statid']][1]//100 or bodies[k]['date'][-1]<active[bodies[k]['statid']][0]//100:
            # if type(bodies[k]['date']) is not list:
            #     bodies[k]['date'] = [bodies[k]['date']]  # as list
            dsec = []
            dssold = ''
            try:
                for ds in [bodies[k]['date'][0], bodies[k]['date'][-1]]:
                    if '-' in ds:
                        if dssold == '':
                            dssold = '-'
                            for dss in ds.split('-'):
                                d = int(dss)
                                dsec.append((datetime(year=d // 10000, month=d % 10000 // 100,
                                                      day=d % 100) - refdate).days * 86400)
                    else:
                        d = int(ds)
                        dsec.append(
                            (datetime(year=d // 10000, month=d % 10000 // 100, day=d % 100) - refdate).days * 86400)
            except:
                return '', 'Invalid date specification' + str(bodies[k]['date'])

            if bodies[k]['statid'] in active.keys():
                if dsec[0] > active[bodies[k]['statid']][1] or dsec[-1] + 86399 < active[bodies[k]['statid']][0]:
                    logger.debug('%s outside Index range', bodies[k]['statid'])
                    del bodies[k]
                    deleted = True

                else:
                    logger.debug('%s Index[%d (%d / %d) %d]', bodies[k]['statid'],
                                 active[bodies[k]['statid']][0],
                                 dsec[0], dsec[-1],
                                 active[bodies[k]['statid']][1])
            else:
                del bodies[k]
                deleted = True

        if not deleted:
            if 'time' in bodies[k].keys():
                if type(bodies[k]['time']) is not list:
                    bodies[k]['time'] = [bodies[k]['time']]
                tsec = []
                tssold = ''
                try:
                    for ds in [bodies[k]['time'][0], bodies[k]['time'][-1]]:
                        if '-' in ds:
                            if tssold == '':
                                tssold = '-'
                                for dss in ds.split('-'):
                                    d = int(dss)
                                    tsec.append(d)
                        else:
                            d = int(ds)
                            tsec.append(d)
                except:
                    return '', 'Invalid time specification: ' + bodies[k]['time']
                logger.debug('tsec: %s', str(tsec))

    #    print(bodies[:5])
    logger.debug('# requests %d', len(bodies))

    if len(bodies) == 0:
        return '', 'No selected station has data in specified date range: ' + str(body)

    func = partial(eua.process_flat, wroot, randdir, cf)

    if False:
        results = list(map(func, bodies))
    else:
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

    try:
        oformat = body['format'][0]
    except:
        oformat = 'nc'
    logger.debug('wpath: %s; format %s', wpath, oformat)

    if False:
        #
        # just return a NetCDf file
        #
        if 'compression' in body.keys():
            rfile = []
            for r in results:
                if len(r[0]) > 0:
                    rfile.append(r[0])

            logger.debug('Request-File: %s [Time: %7.4f s]', rfile, (time.time() - tt))
            return rfile, ''

    if oformat == 'nc':
        with zipfile.ZipFile(rfile, 'w') as f:
            for r in results:
                try:
                    if len(r[0]) > 0:
                        f.write(r[0], os.path.basename(r[0]))
                    os.remove(r[0])
                except:
                    pass

    elif oformat == 'csv':
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
    randdir = '{:012d}'.format(numpy.random.randint(100000000000))
    request_id = id(body)  # or somthing better?
    logger.info("%s POST %s", request_id, str(body))
    wroot = os.path.expandvars('$RSCRATCH/tmp/')
    rfile, error = defproc(body, wroot, randdir, cdm)
    # os.makedirs('logs', exist_ok=True)
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
    active, cdm, cf = main()
    #
    # Parse command line arguments for testing the server API
    #
    body = eval(sys.argv[1])
    #
    # Logging to DEBUG
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
    # ofile=to_csv(idir,glob.glob('*temperature.nc'))
    logger.debug(str(ret))
