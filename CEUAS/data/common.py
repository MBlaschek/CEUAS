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


def to_cdm_coords(data):
    from xarray import DataArray, Dataset
    try:
        import cf2cdm
    except ImportError as e:
        print("Please install cfgrib")
        raise e

    if not isinstance(data, (DataArray, Dataset)):
        raise ValueError("Requires an xarray DataArray or Dataset")

    #
    # make sure there is no standard name on datetime coordinates -> error
    #
    for i, j in data.coords.items():
        if str(j.dtype) == 'datetime64[ns]':
            if 'standard_name' in j.attrs:
                del j.attrs['standard_name']
    return cf2cdm.translate_coords(data, cf2cdm.CDS)


def rename_to_cdm(data, variable, departure=None, cfunits=True, description=None, **kwargs):
    from xarray import DataArray

    if not isinstance(data, DataArray):
        raise ValueError("Requires an xarray DataArray")

    #
    # Convert units to CDM ?
    #
    if cfunits:
        data = to_cdm_coords(data)
    #
    # Search CF Convention for standard_name and units
    #
    attrs = get_cf_convention(variable)

    #
    # departure ?
    #
    if departure:
        #
        # name : ta_fg_dep_era5
        # standard_name : air_temperature_first_guess_departure_era5
        # description : first guess departure OBS - ERA5
        #
        attrs['name'] = attrs['name'] + departure

    #
    # Add description
    #
    if description is not None:
        attrs['description'] = description   # somehow related to departure

    #
    # Add all other keywords as attributes
    #
    attrs.update(kwargs)
    #
    # Check units
    #
    if 'units' in data.attrs:
        if attrs['units'] != data.attrs['units']:
            raise ValueError("Different units found!!!", data.attrs['units'], " to ", attrs['units'])
    #
    # Rename
    #
    data.name = attrs.pop('name', default=data.name)
    data.attrs.update(attrs)
    return data

#
# Standard naming dictionary based on CF Convention (might be faster)
#
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

#
# Download CF Convention and try to find he variable
#
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
        for ivar in doc:
            if varname in ivar['@id']:
                return {'name': ivar['@id'], 'units': ivar['canonical_units'], 'grib': ivar['grib']}
        return {'name': varname, 'units': None, 'grib': None}
