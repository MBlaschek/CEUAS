#!/usr/bin/env python3
import io
import os
import unittest
#import logging
import zipfile

import h5py

from default import defproc, active, cdm, cf

#
# global Informations for Tests
#
# active, cdm, cf = init_server()
# randdir = '{:012d}'.format(100000000000)
randdir = os.path.expandvars('{:012d}'.format(100000000000))
wroot = os.path.expandvars('$RSCRATCH/tmp/')
dzip = wroot + randdir + '/download.zip'

#logger = logging.getLogger('TestLog')
#logger.debug("DEFAULT ZIP: %s", dzip)

default_request = {
    # "date": ["19790101", "19790131"],
    "date": ["20000101", "20000131"],
    "variable": ["temperature"]
}


def update_request(*args, **kwargs):
    """ Allow update and copy of a dict with additonal keys

    Args:
        *args: dict or list(name, value)
        **kwargs: dict

    Returns:
        dict : updated
    """
    if len(args) == 1:
        kwargs.update(*args)
    else:
        kwargs.update({args[0]: args[1]})
    return kwargs


"""
Automated TestCase Builder ? not working at the moment

# test Case Builder Function
def make_request_case(request_body, expr, expected=True, msg=None):
    if msg is None:
        msg = str(request_body)

    class MyTestCase(unittest.TestCase):
        def test_request(self):
            ret = defproc(request_body, wroot, randdir, cdm)
            if expected:
                self.assertTrue(eval(expr), msg=msg)
            else:
                self.assertFalse(eval(expr), msg=msg)


# Dummy Test Case Wrapper
class RequestTest(make_request_case(default_request, 'ret==dzip', expected=True)):
    pass
"""


class BoundingBox(unittest.TestCase):
    bbox = {'bbox': ['40', '0', '50', '20']}

    def test_str_false(self):
        current = update_request({'bbox': "40, 0, 50, 20"}, **default_request)
        ret = defproc(current, wroot, randdir, cdm)
        self.assertFalse(ret[0] == dzip and ret[1] == '', msg=str(ret) + str(current))

    # todo add all definition failures

    def test_bbox_true(self):
        current = update_request(self.bbox, **default_request)
        ret = defproc(current, wroot, randdir, cdm)
        self.assertTrue(ret[0] == dzip and ret[1] == '', msg=str(ret) + str(current))

    def test_bbox_inside_true(self):
        current = update_request(self.bbox, **default_request)
        ret = defproc(current, wroot, randdir, cdm)
        with zipfile.ZipFile(ret[0], 'r') as f:
            x = f.namelist()
            hf = io.BytesIO(f.read(x[0]))
            h = h5py.File(hf, 'r')
            self.assertTrue('ta' in h.keys(), msg=str(ret) + str(current))


class Country(unittest.TestCase):
    country = {'country': 'DEU'}

    def test_country_true(self):
        current = update_request(self.country, **default_request)
        ret = defproc(current, wroot, randdir, cdm)
        self.assertTrue(ret[0] == dzip and ret[1] == '', msg=str(ret) + str(current))


class Statid(unittest.TestCase):
    station = {'statid': '01001'}

    def test_statid_true(self):
        current = update_request(self.station, **default_request)
        ret = defproc(current, wroot, randdir, cdm)
        self.assertTrue(ret[0] == dzip and ret[1] == '', msg=str(ret) + str(current))


class Datetime(unittest.TestCase):

    default_request = default_request.copy()
    default_request['statid'] = '1001'  # Jan Mayen

    def test_date_true(self):
        ret = defproc(default_request, wroot, randdir, cdm)
        self.assertTrue(ret[0] == dzip and ret[1] == '', msg=str(ret) + str(default_request))

    def test_time_true(self):
        current = update_request('time', 6, **default_request)
        ret = defproc(current, wroot, randdir, cdm)
        self.assertTrue(ret[0] == dzip and ret[1] == '', msg=str(ret) + str(current))


class PressureLevel(unittest.TestCase):
    def test_plevel_true(self):
        current = update_request('pressure_level', 500, **default_request)
        ret = defproc(current, wroot, randdir, cdm)
        self.assertTrue(ret[0] == dzip and ret[1] == '', msg=str(ret) + str(current))


class Variables(unittest.TestCase):
    def test_variable_true(self):
        ret = defproc(default_request, wroot, randdir, cdm)
        self.assertTrue(ret[0] == dzip and ret[1] == '', msg=str(ret) + str(default_request))

    def test_feedback_true(self):
        current = update_request('fbstats', ['biascorr'], **default_request)
        ret = defproc(current, wroot, randdir, cdm)
        self.assertTrue(ret[0] == dzip and ret[1] == '', msg=str(ret) + str(current))


#
# def test_dateplot(statstr):
#     evs = '{"statid":"10393","pressure_level":' + statstr + ',"variable":"temperature","fbstats":["obs_minus_bg","bias_estimate"]}'
#     body = eval(evs)
#     ret = defproc(body, wroot, randdir, cdm)
#     readandplot_ts('100000000000/download.zip', eval(evs))
#
#
# def test_dateplot2(statstr):
#     evs = '{"statid":"72764","pressure_level":' + statstr + ',"variable":"dew_point_temperature","fbstats":["obs_minus_bg","bias_estimate"]}'
#     body = eval(evs)
#     ret = defproc(body, wroot, randdir, cdm)
#     readandplot_ts('100000000000/download.zip', eval(evs))


if __name__ == '__main__':
    # for s in ["1002",'"01002"','["01002"]','[1002,1025]','["01002","01025"]','"all"']:
    # test_statid(s)
    # test_bbox2()
    # for s in ['["NOR","DEU"]','"NOR"','"ALL"']: #,'"XXX"']:
    # test_country(s)

    # for s in ['19571001','"19571001"','["19571001"]','[19571001,19571031]','["19571001-19571001"]','"19190101-20190101"','["19571001","19571031"]']: #,'"XXX"']:
    # test_date(s)

    # for d in ['["20140105"]','["20140102-20150105"]','["20130101","20140104","20140105","20140107"]']   :
    # for s in ['','0','"00"','["00"]','[0,12]','["0","12"]','"00-12"','"0-12"','"0-0"','"11-11"','["00-12"]','["21-03"]','["15-03"]','["09-15"]','["18-06"]','["6-18"]']:
    # test_time(s,d)

    # for d in ['','10000','"10000"','[10000]','["10000"]','[10000,20000]','["10000","20000"]']:
    # test_plevel(d)

    # test_dateplot2('50000')
    # test_dateplot('10000')
    #
    # for d in ['"temperature"', '["temperature"]',
    #           '["temperature","u_component_of_wind","v_component_of_wind","wind_speed","wind_direction","relative_humidity","specific_humidity","dew_point_temperature"]']:
    #     for fb in ['', '"obs_minus_an"', '["obs_minus_an"]', '["obs_minus_bg","bias_estimate"]']:
    #         test_variable(d, fb)
    #
    # print('ready')
    unittest.main(verbosity=2)
