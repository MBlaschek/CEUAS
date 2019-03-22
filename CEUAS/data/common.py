# -*- coding: utf-8 -*-

"""
Common declarations:

variable names and meta data

can be used by any data process and employ common names

CDM ? CF
"""

std_plevels = [1000., 2000., 3000., 5000., 7000., 10000., 15000., 20000., 25000., 30000., 40000., 50000., 70000.,
               85000., 92500., 100000.]

era_plevels = [1000., 2000., 3000., 5000., 7000., 10000., 12500., 15000., 17500., 20000., 22500., 25000., 30000.,
               35000., 40000., 45000., 50000., 55000., 60000., 65000., 70000., 75000., 77500., 80000., 82500.,
               85000., 87500., 90000., 92500., 95000., 97500., 100000.]

# todo add more information
# what columns are required
# information on coordinates (pres, hour, time, date) ...

# t_fg_dep_era5  -> Obs - FG

metadata = {
    'temp': {
        'units': 'K',
        'standard_name': 'air_temperature'
    },
    'rhumi': {
        'units': '1',
        'standard_name': 'relative_humidity'
    },
    'dpd': {
        'units': 'K',
        'standard_name': 'dew_point_depression'
    },
    'windd': {
        'units': 'degree',
        'standard_name': 'wind_to_direction'
    },
    'winds': {
        'units': 'm/s',
        'standard_name': 'wind_speed'
    },
    'lon': {
        'units': 'degrees_east',
        'standard_name': 'longitude'
    },
    'lat': {
        'units': 'degrees_north',
        'standard_name': 'latitude'
    },
    'alt': {
        'units': 'm',
        'standard_name': 'altitude_above_sea_level'
    },
    'gph': {
        'units': 'm',
        'standard_name': 'geopotential_height'
    }
}


def get_cf_convention(varname, directory='./data'):
    import xmltodict
    import urllib
    import os

    if not os.path.isfile(directory + '/cf-standard-name-table.xml'):
        os.makedirs(directory, exist_ok=True)
        url = "http://cfconventions.org/Data/cf-standard-names/64/src/cf-standard-name-table.xml"
        try:
            urllib.request.urlretrieve(url, directory + '/cf-standard-name-table.xml')
            print("Downloaded: ", directory + '/cf-standard-name-table.xml')
        except:
            print("Error downloading")

    if os.path.isfile(directory + '/cf-standard-name-table.xml'):
        with open('cf-standard-name-table.xml') as fd:
            doc = xmltodict.parse(fd.read())
        # Use
        doc = doc['standard_name_table']['entry']

        # Search for variables ?? names and units

        return {'name': doc[0]['@id'], 'units': doc[0]['canonical_units'], 'grib': doc[0]['grib']}
