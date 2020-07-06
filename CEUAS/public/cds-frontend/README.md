# Access to the CDS via the cds python API

Related Milestones and Deliverables:


Type | Nr | Responsible | Nature | Title| Due | Status | File
---|---|---|---|---|---|---|---
Deliverable | DC3S311c_Lot2.3.1.1 | UNIVIE | Software, Report | First access to early upper air data base via CDS | March 2020 | June 2020 | code/* 

# Short Description

This directory contains a collection of python scripts and IPython notebook that access the Copernicus Climate Change Service (C3S) - Upper Air Service using the python module cdsapi.

At present, the cdsapi is directed to the sis-dev developer frontend. Once the service has been approved for public use the scripts can be used also via the standard server cds.copernicus.eu.


The front returns files, which are either

1. Zip files containing CF 1.7 compliant netCDF4 files (one per station). This is the default.
2. CSV files, which can be read in directly with spreadsheet applications such as Excel. 
2. JSON files containing error messages, if a HTTP error occurs

 Both file formats can be dealt with in the CDS toolbox. 

| Identifier       | All possible values                                          | Explanation                                                  |
| ---------------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| `date`           | `[YYYYMMDD,YYYYMMDD]`, `YYYYMMDD`, Integer or String         | List of dates of launches                        |
| `period`           | `[YYYYMMDD, YYYYMMDD]`, Integer or String         | Start and End of a period of dates of launches                        |
| `country`        | `[“ALL”,…,”USA”]`, String, Exclusive with `statid`, `bbox`, String | Country codes of stations to be selected according to WMO, see examples below.                    |
| `bbox`           | `[upper,left,lower,right]`, Float or String, Exclusive with `statid`, `country` | Boundaries of lat/lon rectangle to select stations           |
| `fbstats`        | `["obs_minus_bg","obs_minus_an","bias_estimate"]`            | ERA5 feedback information                                    |
| `pressure_level` | `[MMMM,…,NNNNN]`, `MMMM`, Integer or String                  | Pressure levels in Pascal. 16 standard pressure levels (10-1000 hPa) or significant levels (if omitted) |
| `statid`         | `[“SSSSS”]`, String, either WMO or WIGOS IDs. Special value “all”, Exclusive with `country`, `bbox` | WMO or WIGOS station ID                                      |
| `time`           | `[HHMMSS,HHMMSS]`                                            | List of times permitted.                                     |
| `variable`       | `[„air_temperature“, “zonal_wind“, “meridional_wind“, “wind_speed”, ”wind_direction”, ”air_relative_humidity”, ”air_specific_humidity”, "air_dewpoint"]`, String | Meteorological variables                                     |
| `format`         | `nc` or `csv`    | Output format |



# Installation

On the user side, python version 3 needs to be installed. The cdsapi module can be installed using:
```python
pip install cdsapi
```
Information on [pypi](https://pypi.org/project/cdsapi/) or at [cds]/https://cds.climate.copernicus.eu/api-how-to) 
To use the example notebook scripts, please install the following packages:
```python
pip install numpy pandas xarray matplotlib 
```

# Usage

There are examples for requests in the `Example.ipynb` Notebook and two simple requests are given here to show some functionality:

1. Request one sounding at a specific date and at one station

```python
import cdsapi
c = cdsapi.Client(url='https://sis-dev.climate.copernicus.eu/api/v2')   # at the moment this is not in the default catalogue
r = c.retrieve('insitu-comprehensive-upper-air-observation-network',
               {
                   'variable': ["air_temperature", "air_relative_humidity"],
                   'year': '2000',
                   'month':'02',
                   'day':'01',
                   'statid': '10393',
               }, 
               target='download.zip')
```

2. Request a timeseries of on station at a certain pressure level

```python
import cdsapi
c = cdsapi.Client(url='https://sis-dev.climate.copernicus.eu/api/v2')   # at the moment this is not in the default catalogue
r = c.retrieve('insitu-comprehensive-upper-air-observation-network',
               {
                   'variable': ["air_temperature", "air_relative_humidity"],
                   'period': ['19000101', '20201231'],
                   'statid': '10393',
                   'pressure_level': 500
               }, 
               target='download.zip')
```


# License

Generated using Copernicus Climate Change Service Information, 2020
[Copernicus Climate Change Service (C3S), 2020](https://apps.ecmwf.int/datasets/licences/copernicus/)

