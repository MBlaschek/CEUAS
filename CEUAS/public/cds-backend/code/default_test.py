#!/usr/bin/env python3
import io
import os
import unittest
# import logging
import zipfile

import h5py
import pandas as pd
# on import the init_server function will be executed
from default import process_request, wmo_regions, config, cf

#
# global Informations for Tests
#
randdir = os.path.expandvars('{:012d}'.format(100000000000))
output_dir = config['tmp_dir'] + '/' + randdir
dzip = output_dir + '/download.zip'

# logger = logging.getLogger('TestLog')
# logger.debug("DEFAULT ZIP: %s", dzip)

default_request = {
    "date": ["20000101", "20000131"],
    "variable": ["temperature"]
}


def clean_output_dir():
    list(map(os.unlink, (os.path.join(output_dir, f) for f in os.listdir(output_dir))))


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
        self.assertRaises(ValueError,
                          process_request,
                          current, output_dir, config['data_dir'], wmo_regions
                          )

    def test_float_false(self):
        current = update_request({'bbox': [40.0, 30, 30, 50.3]}, **default_request)
        self.assertRaises(ValueError,
                          process_request,
                          current, output_dir, config['data_dir'], wmo_regions
                          )

    def test_true(self):
        current = update_request(self.bbox, **default_request)
        ret = process_request(current, output_dir, config['data_dir'], wmo_regions)
        self.assertTrue(ret == dzip, msg=str(ret) + str(current))
        clean_output_dir()

    def test_data_true(self):
        current = update_request(self.bbox, **default_request)
        ret = process_request(current, output_dir, config['data_dir'], wmo_regions)
        with zipfile.ZipFile(ret, 'r') as f:
            x = f.namelist()
            hf = io.BytesIO(f.read(x[0]))
            h = h5py.File(hf, 'r')
            self.assertTrue('ta' in h.keys(), msg=str(ret) + str(current))
        clean_output_dir()


class Country(unittest.TestCase):
    country = {'country': 'DEU'}

    def test_deu_true(self):
        current = update_request(self.country, **default_request)
        ret = process_request(current, output_dir, config['data_dir'], wmo_regions)
        self.assertTrue(ret == dzip, msg=str(ret) + str(current))
        clean_output_dir()

    def test_atu_false(self):
        current = update_request({'country': 'ATU'}, **default_request)
        self.assertRaises(ValueError,
                          process_request,
                          current, output_dir, config['data_dir'], wmo_regions
                          )


# check more values or the data inside ?

class Statid(unittest.TestCase):
    station = {'statid': '10393'}

    def test_statid_true(self):
        current = update_request(self.station, **default_request)
        ret = process_request(current, output_dir, config['data_dir'], wmo_regions)
        self.assertTrue(ret == dzip, msg=str(ret) + str(current))
        clean_output_dir()

    def test_nostatid_false(self):
        current = dict(variable='temperature', date='20000101')
        ret = process_request(current, output_dir, config['data_dir'], wmo_regions)
        with zipfile.ZipFile(ret, 'r') as f:
            x = f.namelist()
            self.assertTrue(len(x) > 100, msg=str(ret) + str(current))
        clean_output_dir()

    def test_multi_statid_csv_true(self):
        # check for statindex
        pass


class Datetime(unittest.TestCase):
    default_request = default_request.copy()
    default_request['statid'] = '10393'  # Lindenberg

    def test_date_true(self):
        ret = process_request(self.default_request, output_dir, config['data_dir'], wmo_regions)
        self.assertTrue(ret == dzip, msg=str(ret) + str(self.default_request))
        clean_output_dir()

    def test_time_true(self):
        current = update_request('time', 6, **self.default_request)
        ret = process_request(current, output_dir, config['data_dir'], wmo_regions)
        self.assertTrue(ret == dzip, msg=str(ret) + str(current))
        clean_output_dir()

    def test_nodata_true(self):
        current = update_request('date', ['19000101', '19000131'], **self.default_request)
        self.assertRaises(RuntimeError,
                          process_request,
                          current, output_dir, config['data_dir'], wmo_regions
                          )


class PressureLevel(unittest.TestCase):
    default_request = default_request.copy()
    default_request['statid'] = '10393'  # Lindenberg

    def test_50000_true(self):
        current = update_request('pressure_level', 50000, **self.default_request)
        ret = process_request(current, output_dir, config['data_dir'], wmo_regions)
        self.assertTrue(ret == dzip, msg=str(ret) + str(current))
        clean_output_dir()

    def test_120000_false(self):
        current = update_request('pressure_level', 120000, **self.default_request)
        self.assertRaises(ValueError,
                          process_request,
                          current, output_dir, config['data_dir'], wmo_regions
                          )

    def test_type_false(self):
        current = update_request('pressure_level', 50003.2, **self.default_request)
        self.assertRaises(TypeError,
                          process_request,
                          current, output_dir, config['data_dir'], wmo_regions
                          )


class Variables(unittest.TestCase):
    # variables for check_body
    default_variables = ['temperature', 'u_component_of_wind', 'v_component_of_wind',
                         'wind_speed', 'wind_direction', 'relative_humidity',
                         'specific_humidity', 'dew_point_temperature']

    default_cdmname = {}  # this part is really confusing
    for ikey, ival in cf.items():
        if "odbcode" in ival.keys():
            default_cdmname[ival['cdsname']] = ikey
    default_request = default_request.copy()
    default_request['statid'] = '10393'  # Lindenberg

    def test_temperature_true(self):
        ret = process_request(self.default_request, output_dir, config['data_dir'], wmo_regions)
        self.assertTrue(ret == dzip, msg=str(ret) + str(self.default_request))
        clean_output_dir()

    def test_all_variables_true(self):
        current = update_request('variable', self.default_variables, **self.default_request)
        ret = process_request(current, output_dir, config['data_dir'], wmo_regions)
        with zipfile.ZipFile(ret, 'r') as f:
            x = f.namelist()
            for ivar in self.default_variables:
                with self.subTest(ivar=ivar):
                    self.assertTrue(any([self.default_cdmname[ivar] in ifile for ifile in x]),
                                    msg='Missing: %s (%s) %s' % (ivar, self.default_cdmname[ivar], str(x)))
        clean_output_dir()

    def test_variable_false(self):
        self.assertRaises(KeyError,
                          process_request,
                          dict(variable='air_temperature'), output_dir, config['data_dir'], wmo_regions
                          )

    def test_novariable_False(self):
        self.assertRaises(KeyError,
                          process_request,
                          dict(statid='1001'), output_dir, config['data_dir'], wmo_regions
                          )


class Feedback(unittest.TestCase):
    # todo check all variables and data
    default_variables = ['obs_minus_an', 'obs_minus_bg', 'bias_estimate']
    default_request = default_request.copy()
    default_request['statid'] = '10393'  # Lindenberg

    def test_an_true(self):
        current = update_request('fbstats', ['obs_minus_an'], **self.default_request)
        ret = process_request(current, output_dir, config['data_dir'], wmo_regions)
        self.assertTrue(ret == dzip, msg=str(ret) + str(current))
        clean_output_dir()

    def test_data_true(self):
        current = update_request('fbstats', self.default_variables, **self.default_request)
        ret = process_request(current, output_dir, config['data_dir'], wmo_regions)
        with zipfile.ZipFile(ret, 'r') as f:
            x = f.namelist()
            hf = io.BytesIO(f.read(x[0]))
            h = h5py.File(hf, 'r')
            for ivar in self.default_variables:
                with self.subTest(ivar=ivar):
                    self.assertTrue(ivar in h.keys(), msg=ivar + ' not found ' + str(h.keys()))
        clean_output_dir()


class Format(unittest.TestCase):
    default_request = default_request.copy()
    default_request['statid'] = '10393'  # Lindenberg

    def test_format_false(self):
        current = update_request('format', 'h5', **self.default_request)
        self.assertRaises(ValueError,
                          process_request,
                          current, output_dir, config['data_dir'], wmo_regions
                          )

    def test_format_nc_true(self):
        current = update_request('format', 'nc', **self.default_request)
        ret = process_request(current, output_dir, config['data_dir'], wmo_regions)
        with zipfile.ZipFile(ret, 'r') as f:
            x = f.namelist()
            hf = io.BytesIO(f.read(x[0]))
            h = h5py.File(hf, 'r')
            self.assertTrue('ta' in h.keys(), msg=str(ret) + str(current))
        clean_output_dir()

    def test_format_csv_true(self):
        current = update_request('format', 'csv', **self.default_request)
        ret = process_request(current, output_dir, config['data_dir'], wmo_regions)
        with zipfile.ZipFile(ret, 'r') as f:
            x = f.namelist()
            h = pd.read_csv(f.open(x[0]), index_col=0)
            self.assertTrue('ta' in h.columns, msg=str(ret) + str(current))
        clean_output_dir()


class Plotting(unittest.TestCase):

    def test_timeseries(self):
        # temperature at one level, make a long plot and compare to ?
        pass

    def test_feedback(self):
        #
        pass


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
