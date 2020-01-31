# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# This code has been developed in the service contract for C3S
# This is a simple controller
# - index is the default action of any application
#   here it forwards the request from a CDS server to further processing
#   with python
#
#   For queries ask early-upper-air@copernicus-climate.eu
#   The C3S 311c Lot2 group
# - user is required for authentication and authorization
# - download is for downloading files uploaded in the db (does streaming)
# - api is an example of Hypermedia API support and access control
#
# (c) University of Vienna, L. Haimberger, Vienna, Austria
# Copernicus Climate Change Service, 2020
# https://apps.ecmwf.int/datasets/licences/copernicus/
# email leopold.haimberger (at) univie.ac.at
# Created: Vienna, 26 August, 2019
# Last Modifed: 31 January, 2020
# -----------------------------------------------------------------------------

import os
import socket
import sys
import urllib

# from gluon.debug import dbg
host = socket.gethostname()
print(host)
if 'srvx' in host:
    sys.path.append(os.path.expanduser('~/python/'))
else:
    sys.path.append('/data/private/soft/python/')
    os.environ["RSCRATCH"] = "/data/private/"
import cds_eua as eua
import pandas as pd
import numpy
import hug
import h5py  # ickle as h5py
import zipfile
import json
from multiprocessing import Pool
import glob
from functools import partial
from falcon import HTTPError, HTTP_422
import copy
import time
from datetime import datetime


# @hug.exception(Exception)
# def handler(exception):
# return str(exception)

def main():
    os.chdir(os.path.expandvars('$RSCRATCH/era5/odbs/1'))
    z = numpy.zeros(1, dtype=numpy.int32)
    zidx = numpy.zeros(1, dtype=numpy.int)
    idx = numpy.zeros(1, dtype=numpy.int)
    trajectory_index = numpy.zeros(1, dtype=numpy.int)
    zz = eua.calc_trajindexfast(z, zidx, idx, trajectory_index)

    try:
        with open('activex.json') as f:
            active = json.load(f)
            slnum = active.keys()
            slist = [os.path.expandvars('$RSCRATCH/era5/odbs/1/') + 'chera5.conv._' + s + '.nc' for s in slnum]
    except:
        slist = glob.glob('/raid60/scratch/leo/scratch/era5/odbs/1/chera5.conv._?????.nc')
        slnum = [i[-8:-3] for i in slist]

        volapath = 'https://oscar.wmo.int/oscar/vola/vola_legacy_report.txt'
        f = urllib.request.urlopen(volapath)
        col_names = pd.read_csv(f, delimiter='\t', quoting=3, nrows=0)
        f = urllib.request.urlopen(volapath)
        tdict = {col: str for col in col_names}
        vola = pd.read_csv(f, delimiter='\t', quoting=3, dtype=tdict, na_filter=False)

        active = {}

        for s, skey in zip(slist, slnum):
            with h5py.File(s, 'r') as f:
                print(skey)
                try:

                    active[skey] = [int(f['recordtimestamp'][0]), int(f['recordtimestamp'][-1]),
                                    float(f['observations_table']['latitude'][-1]),
                                    float(f['observations_table']['longitude'][-1])]
                    idx = numpy.where(vola.IndexNbr.values == skey)[0]
                    if len(idx) > 0:

                        active[skey].append(vola.CountryCode[idx[0]])
                    else:
                        active[skey].append('')
                        print('no key found for ' + skey)
                except KeyError:
                    print(skey + ': a table is missing')
        with open('active.json', 'w') as f:
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


# filelist=glob.glob('chera5.conv._?????.nc')

def check_body(body, cdm):
    required_keys = ['variable']
    valid_keys = required_keys + ['statid', 'fbstats', 'pressure_level', 'date', 'time', 'bbox', 'country']
    xor_keys = ['statid', 'bbox', 'country']
    valid_ranges = {}
    valid_ranges['statid'] = [1001, 99999]
    valid_ranges['country'] = ['globe', 'all']
    valid_ranges['bbox'] = [[-180, 180], [-90, 90], [-180, 180], [-90, 90]]  # luro
    valid_ranges['pressure_level'] = [500, 100000]
    valid_ranges['date'] = [19000101, 20301231]
    valid_ranges['time'] = [0, 230000]
    valid_ranges['fbstats'] = ['obs_minus_an', 'obs_minus_bg', 'bias_estimate']
    valid_ranges['variable'] = ['temperature', 'u_component_of_wind', 'v_component_of_wind',
                                'wind_speed', 'wind_direction', 'relative_humidity',
                                'specific_humidity', 'dew_point_temperature', 'geopotential']

    try:

        bk = list(body.keys())
        for r in required_keys:
            if r not in bk:
                return 'argument ' + r + ' is missing'

        for b in bk:
            if b not in valid_keys:
                return 'argument name ' + b + ' is invalid. Valid arguments:' + str(valid_keys)

        rxor = ''
        for r in xor_keys:
            if r in bk:
                if len(rxor) == 0:

                    rxor = r
                else:
                    return 'Please do not specify both ' + rxor + ' and ' + r

        if 'country' in bk:
            vcountries = cdm['sub_region'].alpha_3_code.values
            body['statid'] = []
            for v in body['country']:
                if v not in vcountries:
                    return v + ' is not a valid country code'
                for k, vv in active.items():
                    if vv[4] == v:
                        body['statid'].append(k)
            if len(body['statid']) == 0:
                return 'Countries ' + str(body['country']) + ' have no radiosondes'
            del body['country']
        elif 'bbox' in bk:
            if type(body['bbox']) is not list:
                return 'Bounding box: [lower, left, upper, right]'
            if len(body['bbox']) != 4:
                return 'Bounding box: [lower, left, upper, right]'
            try:
                for i in range(4):
                    body['bbox'][i] = float(body['bbox'][i])
            except:
                return 'Bounding box: [lower, left, upper, right] must be int or float'
            if body['bbox'][0] > body['bbox'][2] or body['bbox'][1] > body['bbox'][3]:
                return 'Bounding box requirements: lower<upper, left<right, -90<=lat<=90, -180<=lon<=360'
            if body['bbox'][0] < -90 or body['bbox'][0] > 90 or body['bbox'][2] < -90 or body['bbox'][2] > 90 or \
                    body['bbox'][1] < -180 or body['bbox'][1] > 360 or body['bbox'][3] < -180 or body['bbox'][3] > 360 \
                    or body['bbox'][3] - body['bbox'][1] > 360:
                return 'Bounding box requirements: lower<upper, left<right, -90<=lat<=90, -180<=lon<=360'

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
                return 'Bounding box ' + str(body['bbox']) + ' contains no radiosondes'
            del body['bbox']
        else:
            try:
                if body['statid'] == 'all':
                    body['statid'] = slnum
            except:
                return 'Please specify either bbox, country or statid for station selection. Use "statid":"all" to select all stations'
        if body['variable'] == 'all':
            body['variable'] = ['temperature', 'u_component_of_wind', 'v_component_of_wind']

        bk = list(body.keys())
        for v in bk:
            if type(valid_ranges[v][0]) is str:
                if type(body[v]) is not list:
                    body[v] = [body[v], body[v]]
                for bv in body[v]:
                    if bv not in valid_ranges[v]:
                        return ['argument value(s) ' + str(body[v]) + ' not valid.',
                                'Valid values:' + str(valid_ranges[v])]
                print('bodyx:', v, body[v][0], len(body[v]))
                if body[v][0] == body[v][1]:
                    body[v].pop()

            else:
                if type(body[v]) is not list:
                    body[v] = [body[v], body[v]]
                for i in range(len(body[v])):  # bv in body[v]:
                    bv = body[v][i]
                    print(bv)
                    if int(bv) < valid_ranges[v][0] or int(bv) > valid_ranges[v][1]:
                        return ['argument value(s) ' + str(body[v]) + ' not valid.',
                                'Valid values:' + str(valid_ranges[v])]
                    if v != 'statid':
                        body[v][i] = int(body[v][i])
                print('body:', v, body[v][0], len(body[v]))
                if body[v][0] == body[v][1]:
                    body[v].pop()
    except IOError:
        return ['general syntax error ', body]

    return ''


def makebodies(bodies, body, spv, bo, l):
    for b in body[spv[l]]:

        if l < len(spv) - 1:
            makebodies(bodies, body, spv, copy.copy(bo) + [b], l + 1)
        else:
            bodies.append(dict(body))
            bn = copy.copy(bo) + [b]
            print('spv,bo:', spv, bn)
            for s, b in zip(spv, bn):
                bodies[-1][s] = b
                print('makebodies', l, s, b)


def defproc(body, randdir, cdm):
    tt = time.time()
    error = check_body(body, cdm)
    print('body', body)
    if len(error) > 0:
        return '', error

    try:
        os.mkdir(randdir)
    except:
        pass

    bodies = []
    spv = ['statid', 'variable']
    bo = []
    makebodies(bodies, body, spv, bo, 0)
    for k in range(len(bodies) - 1, -1, -1):
        if 'date' in bodies[k].keys():
            # if bodies[k]['date'][0]>active[bodies[k]['statid']][1]//100 or bodies[k]['date'][-1]<active[bodies[k]['statid']][0]//100:
            dsec = []
            for d in bodies[k]['date']:
                dsec.append(((datetime(year=d // 10000, month=d % 10000 // 100, day=d % 100) - datetime(year=1900,
                                                                                                        month=1,
                                                                                                        day=1))).days * 86400)
            if dsec[0] > active[bodies[k]['statid']][1] or dsec[-1] + 86399 < active[bodies[k]['statid']][0]:
                del bodies[k]

    print(bodies[:5])
    print('len:', len(bodies))

    # spvs={'statid':body.pop('statid'),'variable':body.pop('variable')}
    # for sv1 in body['statid']:,'variable']:
    # for sv1 in ['statid','variable']:

    # if type(body[splitvar]) is list:
    # for s in  body[splitvar]:
    # bodies.append(dict(body))
    # bodies[-1]['splitvar']=s
    # else:
    # if body[splitvar] == 'all':
    # statids=[]
    # for s in slist:
    # bodies.append(dict(body))
    # bodies[-1][splitvar]=s[-8:-3]
    # else:
    # bodies.append(dict(body))

    p = Pool(10)
    func = partial(eua.process_flat, randdir, cf)

    results = list(map(func, bodies))
    del p
    wpath = ''
    for r in results:
        if r[0] != '':
            wpath = r[0]
            break
    if wpath == '':
        raise HTTPError(HTTP_422, description=[results[0][1], body])

    else:
        rfile = os.path.dirname(wpath) + '/download.zip'

    print(results, 'wpath:' + wpath + 'x')
    with zipfile.ZipFile(rfile, 'w') as f:
        for r in results:
            try:
                if len(r[0]) > 0:
                    f.write(r[0], os.path.basename(r[0]))

                os.remove(r[0])
            except:
                pass

    # for r in results:
    # os.remove(r[0])

    print('rfile:', rfile, 'time: {:7.4f}s'.format(time.time() - tt))

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

    print(request.query_string)
    if '=' not in request.query_string:
        raise HTTPError(HTTP_422, description='A query string must be supplied')

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

    randdir = '{:012d}'.format(numpy.random.randint(100000000000))

    print(body)

    rfile, error = defproc(body, randdir)

    if rfile == '':
        raise HTTPError(HTTP_422, description=error)

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

    print(json.dumps(body))

    rfile, error = defproc(body, randdir, cdm)
    print(rfile, error)
    if rfile == '':
        raise HTTPError(HTTP_422, description=error)

    response.set_header('Content-Disposition', 'attachment; filename=' + os.path.basename(rfile))

    return rfile


if __name__ == '__main__':
    # active,cdm,cf=main()
    body = eval(sys.argv[1])

    randdir = '{:012d}'.format(100000000000)
    ret = defproc(body, randdir, cdm)
    print(ret)
    print('ready')
