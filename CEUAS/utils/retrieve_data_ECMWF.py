# -*- coding: utf-8 -*-

"""
Retrieve Data from ECMWF open access data



MARS Request:

retrieve,
class=ei,
date=19790101/19790201/19790301/19790401/19790501/19790601/19790701/19790801/19790901/19791001/19791101/19791201/19800101/19800201/19800301/19800401/19800501/19800601/19800701/19800801/19800901/19801001/19801101/19801201/19810101/19810201/19810301/19810401/19810501/19810601/19810701/19810801/19810901/19811001/19811101/19811201/19820101/19820201/19820301/19820401/19820501/19820601/19820701/19820801/19820901/19821001/19821101/19821201/19830101/19830201/19830301/19830401/19830501/19830601/19830701/19830801/19830901/19831001/19831101/19831201/19840101/19840201/19840301/19840401/19840501/19840601/19840701/19840801/19840901/19841001/19841101/19841201/19850101/19850201/19850301/19850401/19850501/19850601/19850701/19850801/19850901/19851001/19851101/19851201/19860101/19860201/19860301/19860401/19860501/19860601/19860701/19860801/19860901/19861001/19861101/19861201/19870101/19870201/19870301/19870401/19870501/19870601/19870701/19870801/19870901/19871001/19871101/19871201/19880101/19880201/19880301/19880401/19880501/19880601/19880701/19880801/19880901/19881001/19881101/19881201/19890101/19890201/19890301/19890401/19890501/19890601/19890701/19890801/19890901/19891001/19891101/19891201/19900101/19900201/19900301/19900401/19900501/19900601/19900701/19900801/19900901/19901001/19901101/19901201/19910101/19910201/19910301/19910401/19910501/19910601/19910701/19910801/19910901/19911001/19911101/19911201/19920101/19920201/19920301/19920401/19920501/19920601/19920701/19920801/19920901/19921001/19921101/19921201/19930101/19930201/19930301/19930401/19930501/19930601/19930701/19930801/19930901/19931001/19931101/19931201/19940101/19940201/19940301/19940401/19940501/19940601/19940701/19940801/19940901/19941001/19941101/19941201/19950101/19950201/19950301/19950401/19950501/19950601/19950701/19950801/19950901/19951001/19951101/19951201/19960101/19960201/19960301/19960401/19960501/19960601/19960701/19960801/19960901/19961001/19961101/19961201/19970101/19970201/19970301/19970401/19970501/19970601/19970701/19970801/19970901/19971001/19971101/19971201/19980101/19980201/19980301/19980401/19980501/19980601/19980701/19980801/19980901/19981001/19981101/19981201/19990101/19990201/19990301/19990401/19990501/19990601/19990701/19990801/19990901/19991001/19991101/19991201/20000101/20000201/20000301/20000401/20000501/20000601/20000701/20000801/20000901/20001001/20001101/20001201/20010101/20010201/20010301/20010401/20010501/20010601/20010701/20010801/20010901/20011001/20011101/20011201/20020101/20020201/20020301/20020401/20020501/20020601/20020701/20020801/20020901/20021001/20021101/20021201/20030101/20030201/20030301/20030401/20030501/20030601/20030701/20030801/20030901/20031001/20031101/20031201/20040101/20040201/20040301/20040401/20040501/20040601/20040701/20040801/20040901/20041001/20041101/20041201/20050101/20050201/20050301/20050401/20050501/20050601/20050701/20050801/20050901/20051001/20051101/20051201/20060101/20060201/20060301/20060401/20060501/20060601/20060701/20060801/20060901/20061001/20061101/20061201/20070101/20070201/20070301/20070401/20070501/20070601/20070701/20070801/20070901/20071001/20071101/20071201/20080101/20080201/20080301/20080401/20080501/20080601/20080701/20080801/20080901/20081001/20081101/20081201/20090101/20090201/20090301/20090401/20090501/20090601/20090701/20090801/20090901/20091001/20091101/20091201/20100101/20100201/20100301/20100401/20100501/20100601/20100701/20100801/20100901/20101001/20101101/20101201/20110101/20110201/20110301/20110401/20110501/20110601/20110701/20110801/20110901/20111001/20111101/20111201/20120101/20120201/20120301/20120401/20120501/20120601/20120701/20120801/20120901/20121001/20121101/20121201/20130101/20130201/20130301/20130401/20130501/20130601/20130701/20130801/20130901/20131001/20131101/20131201/20140101/20140201/20140301/20140401/20140501/20140601/20140701/20140801/20140901/20141001/20141101/20141201/20150101/20150201/20150301/20150401/20150501/20150601/20150701/20150801/20150901/20151001/20151101/20151201/20160101/20160201/20160301/20160401/20160501/20160601/20160701/20160801/20160901/20161001/20161101/20161201/20170101/20170201/20170301/20170401/20170501/20170601/20170701/20170801/20170901/20171001/20171101/20171201,
expver=1,
levtype=sfc,
param=151.128/167.128,
stream=moda,
type=an,
dataset=interim,
target=special/era_interim_mslp_moda_1979-2017.grb

retrieve,
class=ei,
date=19790101/19790201/19790301/19790401/19790501/19790601/19790701/19790801/19790901/19791001/19791101/19791201/19800101/19800201/19800301/19800401/19800501/19800601/19800701/19800801/19800901/19801001/19801101/19801201/19810101/19810201/19810301/19810401/19810501/19810601/19810701/19810801/19810901/19811001/19811101/19811201/19820101/19820201/19820301/19820401/19820501/19820601/19820701/19820801/19820901/19821001/19821101/19821201/19830101/19830201/19830301/19830401/19830501/19830601/19830701/19830801/19830901/19831001/19831101/19831201/19840101/19840201/19840301/19840401/19840501/19840601/19840701/19840801/19840901/19841001/19841101/19841201/19850101/19850201/19850301/19850401/19850501/19850601/19850701/19850801/19850901/19851001/19851101/19851201/19860101/19860201/19860301/19860401/19860501/19860601/19860701/19860801/19860901/19861001/19861101/19861201/19870101/19870201/19870301/19870401/19870501/19870601/19870701/19870801/19870901/19871001/19871101/19871201/19880101/19880201/19880301/19880401/19880501/19880601/19880701/19880801/19880901/19881001/19881101/19881201/19890101/19890201/19890301/19890401/19890501/19890601/19890701/19890801/19890901/19891001/19891101/19891201/19900101/19900201/19900301/19900401/19900501/19900601/19900701/19900801/19900901/19901001/19901101/19901201/19910101/19910201/19910301/19910401/19910501/19910601/19910701/19910801/19910901/19911001/19911101/19911201/19920101/19920201/19920301/19920401/19920501/19920601/19920701/19920801/19920901/19921001/19921101/19921201/19930101/19930201/19930301/19930401/19930501/19930601/19930701/19930801/19930901/19931001/19931101/19931201/19940101/19940201/19940301/19940401/19940501/19940601/19940701/19940801/19940901/19941001/19941101/19941201/19950101/19950201/19950301/19950401/19950501/19950601/19950701/19950801/19950901/19951001/19951101/19951201/19960101/19960201/19960301/19960401/19960501/19960601/19960701/19960801/19960901/19961001/19961101/19961201/19970101/19970201/19970301/19970401/19970501/19970601/19970701/19970801/19970901/19971001/19971101/19971201/19980101/19980201/19980301/19980401/19980501/19980601/19980701/19980801/19980901/19981001/19981101/19981201/19990101/19990201/19990301/19990401/19990501/19990601/19990701/19990801/19990901/19991001/19991101/19991201/20000101/20000201/20000301/20000401/20000501/20000601/20000701/20000801/20000901/20001001/20001101/20001201/20010101/20010201/20010301/20010401/20010501/20010601/20010701/20010801/20010901/20011001/20011101/20011201/20020101/20020201/20020301/20020401/20020501/20020601/20020701/20020801/20020901/20021001/20021101/20021201/20030101/20030201/20030301/20030401/20030501/20030601/20030701/20030801/20030901/20031001/20031101/20031201/20040101/20040201/20040301/20040401/20040501/20040601/20040701/20040801/20040901/20041001/20041101/20041201/20050101/20050201/20050301/20050401/20050501/20050601/20050701/20050801/20050901/20051001/20051101/20051201/20060101/20060201/20060301/20060401/20060501/20060601/20060701/20060801/20060901/20061001/20061101/20061201/20070101/20070201/20070301/20070401/20070501/20070601/20070701/20070801/20070901/20071001/20071101/20071201/20080101/20080201/20080301/20080401/20080501/20080601/20080701/20080801/20080901/20081001/20081101/20081201/20090101/20090201/20090301/20090401/20090501/20090601/20090701/20090801/20090901/20091001/20091101/20091201/20100101/20100201/20100301/20100401/20100501/20100601/20100701/20100801/20100901/20101001/20101101/20101201/20110101/20110201/20110301/20110401/20110501/20110601/20110701/20110801/20110901/20111001/20111101/20111201/20120101/20120201/20120301/20120401/20120501/20120601/20120701/20120801/20120901/20121001/20121101/20121201/20130101/20130201/20130301/20130401/20130501/20130601/20130701/20130801/20130901/20131001/20131101/20131201/20140101/20140201/20140301/20140401/20140501/20140601/20140701/20140801/20140901/20141001/20141101/20141201/20150101/20150201/20150301/20150401/20150501/20150601/20150701/20150801/20150901/20151001/20151101/20151201/20160101/20160201/20160301/20160401/20160501/20160601/20160701/20160801/20160901/20161001/20161101/20161201/20170101/20170201/20170301/20170401/20170501/20170601/20170701/20170801/20170901/20171001/20171101/20171201,
expver=1,
levtype=sfc,
param=144.128/228.128,
step=0-12,
stream=mdfa,
type=fc,
dataset=interim,
target=special/era_interim_total_precip_mdfa_1979-2017.grb

"""
# exp='{0}'.format(2367+(year-1908)/8)
# 'expver': '2367',
CERA20C = {
    'class': "ep",
    'dataset': "cera20c",
    'number': "0/1/2/3/4/5/6/7/8/9",
    'stream': "enda",
    'type': "an",
    'levtype': "pl",
    'levelist': '10/20/30/50/70/100/150/200/250/300/400/500/700/850/925/1000',
    'time': "00/06/12/18",
    'date': '19010101/to/19010131',
    'param': "u/v/T/q/r",
    'area': "G",
    'grid': "1/1",
    'target': 'CERA20C190101.grb'
}


def now():
    """
    Functions returns the current date string
    Returns
    -------
    str
        current date
    """
    import datetime
    return datetime.datetime.now().isoformat()


def request_4d_field(startdate, enddate, variables="t/q/o3", levellist="all", levtype="pl"):
    """
    MARS request for 4D Field of Temeprature (time x levels x lat x lon)

    Parameters
    ----------
    startdate
    enddate
    variables
    levellist
    levtype

    Returns
    -------
    mars_request
    """
    import pandas as pd
    date = pd.date_range(startdate, enddate, freq='D').to_period('M').unique().strftime("%Y%m")
    if len(date) > 1:

        print("Warning, requested Data might be too large for one request (2.? GB limit)")
        print("Number of Month: ", len(date))
        date = date[0] + '-' + date[-1]
    else:
        date = date[0]
    print(date)

    return {'class': 'ei',
            'stream': 'oper',
            'dataset': 'interim',
            'expver': '1',
            'resol': 'av',
            'grid': '128',
            'type': "an",
            'levtype': levtype,
            'levellist': levellist,
            'param': variables,
            'date': "%s/to/%s" % (startdate, enddate),
            'time': "00/12",
            'format': 'netcdf',
            'target': "lev/ERAINTERIM_LEV_%s.nc" % date}


def request_surface(startdate, enddate):
    """
    # MARS Table 128
    31 Sea-ice cover
    165 10m u-wind
    166 10m v-wind
    167 2m T
    134 Surface pressure
    235 Skin Temperature

    Parameters
    ----------
    startdate
    enddate

    Returns
    -------
    mars_request
    """
    import pandas as pd
    date = pd.date_range(startdate, enddate, freq='D').to_period('M').unique().strftime("%Y%m")
    if len(date) > 1:
        print("Warning, requested Data might be too large for one request (2.? GB limit)")
        date = date[0] + '-' + date[-1]

    else:
        date = date[0]

    return {'class': 'ei',
            'stream': 'oper',
            'dataset': 'interim',
            'expver': '1',
            'type': "an",
            'levtype': 'sfc',
            'resol': 'av',
            'grid': '128',
            'format': 'netcdf',
            'param': '134.128/165.128/166.128/167.128/235.128/31.128',
            'step': '0',
            'date': "%s/to/%s" % (startdate, enddate),
            'time': "00/12",
            'target': "sfc/ERAINTERIM_SFC_%s.nc" % date}


def retrieve(request, verbose=1, dryrun=False):
    """ Execute Mars request

    Parameters
    ----------
    request

    Returns
    -------
    filename
    """
    import os
    try:
        from ecmwfapi import ECMWFDataServer

    except ImportError as e:
        print("Install ECMWF webapi")
        raise e

    server = ECMWFDataServer()
    # server.verbose = True   # if an error occurs

    if not os.path.isdir(os.path.dirname(request['target'])):
        os.makedirs(os.path.dirname(request['target']))

    try:
        if dryrun:
            print(request)

        else:
            server.retrieve(request)

        if verbose > 0:
            print("Request was successful.")
            print(request['target'])
        return True

    except Exception as e:
        print(repr(e))
        return False


def retrieve_batch(startdate, enddate, request, pieces='M', dryrun=False):
    """ Batch process retrieval

    Modifiles 'date' and 'target' in request

    Parameters
    ----------
    startdate       str     e.g. '1979-01-01'
    enddate         str     e.g. '2014-12-31'
    request         dict    request_4d_field
    pieces          str     'M' : Monthly

    Returns
    -------
    status
    """
    import os
    import pandas as pd

    if 'date' not in request:
        print(request)
        raise RuntimeError('MARS request need date')

    if 'target' not in request:
        print(request)
        raise RuntimeError('MARS request need target')

    date = pd.date_range(startdate, pd.to_datetime(enddate) + pd.tseries.offsets.MonthEnd(1), freq='D')
    mdate = date.to_period(pieces).unique().to_timestamp()
    filedir = os.path.dirname(request['target'])
    ifile, iext = os.path.splitext(request['target'])
    if ifile[-1] != '_':
        ifile += '_'

    if not os.path.isdir(filedir):
        os.makedirs(filedir)

    n = len(mdate)
    success = {}
    for i, idate in enumerate(mdate):
        request['target'] = ifile + idate.strftime("%Y%m") + iext
        if i + 1 == n:
            idates = date[(date >= mdate[i])]  # & (date <=enddate)]

        else:
            idates = date[(date >= mdate[i]) & (date < mdate[i + 1])]

        request['date'] = "%s/to/%s" % (idates[0].strftime('%Y%m%d'), idates[-1].strftime('%Y%m%d'))
        print("[%04d / %04d] Retrieving: %s" % (i + 1, n, request['date']))
        status = retrieve(request, verbose=0, dryrun=dryrun)
        # print request['target']
        status = True
        success[request['date']] = status
        if not status:
            print("[%04d / %04d] ERROR: %s" % (i + 1, n, request['date']))

    return pd.Series(success)


def default_request(startdate, enddate, surface=True, level=True, batch=True, timestep='M', verbose=0, dryrun=False):
    req_sfc = request_surface(startdate, enddate)
    req_lev = request_4d_field(startdate, enddate)

    # START
    if batch:
        # Update Dates
        req_sfc['target'] = drop_dates(req_sfc['target'])
        req_lev['target'] = drop_dates(req_lev['target'])
        if surface:
            retrieve_batch(startdate, enddate, req_sfc, pieces=timestep, dryrun=dryrun)

        if level:
            retrieve_batch(startdate, enddate, req_lev, pieces=timestep, dryrun=dryrun)
    else:
        if surface:
            retrieve(req_sfc, verbose=verbose, dryrun=dryrun)

        if level:
            retrieve(req_lev, verbose=verbose, dryrun=dryrun)

    # DONE


def drop_dates(text):
    i = text.split('_')
    return "_".join(i[:-1]) + "." + i[-1].split('.')[-1]


def read_mars_request(req_file):
    import os

    if not os.path.isfile(req_file):
        raise IOError('File not found')

    with open(req_file, 'r') as f:
        req = f.read()

    req = req.splitlines()

    out = {}

    for iline in req:
        if 'retrieve' in iline:
            continue
        else:
            p = iline.split('=')
            out[p[0]] = p[1].replace(',', '')

    return out


def custom_request(startdate, enddate, req_file, batch=True, timestep='M', verbose=0, dryrun=False):
    req = read_mars_request(req_file)
    if verbose > 0:
        print(req)

    if batch:
        retrieve_batch(startdate, enddate, req, pieces=timestep, dryrun=dryrun)
    else:
        retrieve(req, verbose=verbose, dryrun=dryrun)


def usage():
    print(""" Retrieve ERA-Interim Data
retrieve_data_ECMWF.py -h -s [date] -e [date] -f [] -b -t [] -v -n

Options:
     -h             HELP
     -s [date]      Start date
     -e [date]      End date
     -f [file]      use mars request file
     -b             Batch, split into several tasks
     -t [M,A,D]     split into Monthly, Annual or daily
     -v             verbose
     -n             dryrun
""")


if __name__ == '__main__':
    import getopt, os, sys

    try:
        opts, args = getopt.getopt(sys.argv[1:], "s:e:f:bt:vn", ["sfc", "lev"])

    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit(2)

    batch = False
    surface = False
    level = False
    startdate = None
    enddate = None
    req_file = ""
    timestep = 'M'
    verbose = 0
    dryrun = False

    for opt, arg in opts:
        if opt in ("-s", "--start"):
            startdate = arg

        elif opt in ("-e", "--end"):
            enddate = arg

        elif opt in ("-f", "--file"):
            req_file = arg

        elif opt in ("-b", "--batch"):
            batch = True

        elif opt in ("--sfc"):
            surface = True

        elif opt in ("--lev"):
            level = True

        elif opt in ("-t", "--time"):
            timestep = arg

        elif opt in ("-v"):
            verbose = 1

        elif opt in ("-n"):
            dryrun = True
        else:
            usage()
            assert False, "unhandled option"

    if timestep not in ['M', 'A', 'D']:
        raise ValueError("Timestep not in M, A, D")

    # Check if DATES are there for some applications
    if batch or req_file == "":
        if startdate is None or enddate is None:
            usage()
            raise ValueError("Requires STARTDATE and ENDDATE")

    if req_file == "":
        default_request(startdate, enddate, surface=surface, level=level, timestep=timestep, batch=batch, dryrun=dryrun)
    else:
        custom_request(startdate, enddate, req_file, batch=batch, timestep=timestep, verbose=verbose, dryrun=dryrun)

    with open('retrieve_history', 'a') as f:
        f.write(now() + " ".join(sys.argv) + "\n")
