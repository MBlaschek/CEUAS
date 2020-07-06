# Access to the CDS via the cds python API

Related Milestones and Deliverables:


Type | Nr | Responsible | Nature | Title| Due | Status | File
---|---|---|---|---|---|---|---
Deliverable | DC3S311c_Lot2.3.1.1 | UNIVIE | Software, Report | First access to early upper air data base via CDS | March 2020 | July 2020 | code/* 

# Short Description

This directory contains a collection of python scripts and an IPython notebook that access the Copernicus Climate Change Service (C3S) - Upper Air Service using the CDSAPI in Python.

At present, the cdsapi is directed to the sis-dev developer frontend. Once the service has been approved for public use the scripts can be used also via the standard server cds.copernicus.eu.

The front returns files, which are either

1. Zip file containing CF 1.7 compliant netCDF4 files (one per station). This is the default.
2. Zip file containing CSV files, which can be read in directly with spreadsheet applications such as Excel. 
3. JSON file containing error messages, if an HTTP error occurs

Both file formats can be dealt with in the CDS toolbox. 

| Identifier       | All possible values                                          | Explanation                                                  |
| ---------------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| `date`           | `[YYYYMMDD,YYYYMMDD]`, `YYYYMMDD`, Integer or String         | List of dates of launches                        |
| `period`           | `[YYYYMMDD, YYYYMMDD]`, Integer or String         | Start and End of a period of dates of launches                        |
| `country`        | `[“ALL”,…,”USA”]`, String, Exclusive with `statid`, `bbox` | Country codes of stations to be selected according to WMO, see examples below.                    |
| `bbox`           | `[upper,left,lower,right]`, Float or String, Exclusive with `statid`, `country` | Boundaries of lat/lon rectangle to select stations           |
| `fbstats`        | `["obs_minus_bg","obs_minus_an","bias_estimate"]`            | ERA5 feedback information                                    |
| `pressure_level` | `[MMMM,…,NNNNN]`, `MMMM`, Integer or String                  | Pressure levels in Pascal. 16 standard pressure levels (10-1000 hPa) or significant levels (if omitted) |
| `statid`         | `[“SSSSS”]`, String, either WMO or WIGOS IDs. Special value “all”, Exclusive with `country`, `bbox` | WMO or WIGOS station ID                                      |
| `time`           | `[HHMMSS,HHMMSS]`                                            | List of times permitted.                                     |
| `variable`       | `[„air_temperature“, “zonal_wind“, “meridional_wind“, “wind_speed”, ”wind_direction”, ”air_relative_humidity”, ”air_specific_humidity”, "air_dewpoint"]`, String | Meteorological variables                                     |
| `format`         | `nc` or `csv`    | Output format |



# Installation

On the user side, Python (advised version 3) needs to be installed. The cdsapi module can be installed using:
```python
pip install cdsapi
```
Information on [pypi](https://pypi.org/project/cdsapi/) or at [cds]/https://cds.climate.copernicus.eu/api-how-to) 
To use the example notebook scripts, please install the following packages:
```python
pip install numpy pandas xarray matplotlib 
```

# Usage

There are examples for different requests in the `Example.ipynb` Notebook and two simple requests are given here to show some functionality:

1. Request one sounding at a specific date and at one station:
We will ask for Lindenberg station on 1.February 2000. This will return a zip-File with two NetCDF files inside. One for air temperaure and one for air relative humidity. We specified not `pressure_level`, therefore all available pressure levels will be requested. We do not ask for a specific time, but all soundings on that day.
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
Executing this yields the first time (afterwards it is buffered for some time):
```python
2020-07-06 16:10:57,374 INFO Sending request to https://sis-dev.climate.copernicus.eu/api/v2/resources/insitu-comprehensive-upper-air-observation-network
2020-07-06 16:10:57,821 INFO Request is completed
2020-07-06 16:10:57,946 INFO Downloading http://***.***.***.***/cache-compute-0002/cache/data2/adaptor.comprehensive_upper_air.retrieve-1593518168.404673-17157-1-e369f43b-1db4-45bd-aa35-d909d7e46443.zip to download.zip (52.4K)
2020-07-06 16:10:58,072 INFO Download rate 418K/s
```
And unzipping:
```python
import zipfile
z = zipfile.ZipFile('download.zip')
z.extractall(path='./example_data/2')
z.close()
```
Finally we have two files in out directory: `example_data/2/dest_0-20000-0-10393_air_temperature.nc` and `example_data/2/dest_0-20000-0-10393_relative_humidity.nc`,
which we can read in with `Xarray` and have a look at the data:
```python
import xarray as xr
tdata = xr.load_dataset('example_data/2/dest_0-20000-0-10393_air_temperature.nc')
print(tdata)
<xarray.Dataset>
Dimensions:           (obs: 227, string10: 10, trajectory: 4)
Coordinates:
    lat               (obs) float32 52.22 52.22 52.22 ... 52.22 52.22 52.22
    lon               (obs) float32 14.12 14.12 14.12 ... 14.12 14.12 14.12
  * obs               (obs) int32 0 0 0 0 0 0 0 0 0 0 0 ... 0 0 0 0 0 0 0 0 0 0
    plev              (obs) float32 4260.0 5000.0 5180.0 ... 96100.0 100000.0
  * string10          (string10) |S1 b'' b'' b'' b'' b'' b'' b'' b'' b'' b''
    time              (obs) datetime64[ns] 2000-02-01T05:00:00 ... 2000-02-01T23:00:00
  * trajectory        (trajectory) int32 0 0 0 0
Data variables:
    ta                (obs) float32 200.3 203.7 nan 205.3 ... 278.2 278.2 nan
    trajectory_index  (obs) int32 0 0 0 0 0 0 0 0 0 0 0 ... 3 3 3 3 3 3 3 3 3 3
    trajectory_label  (trajectory, string10) |S1 b'0' b'0' b'0' ... b'9' b'2'
Attributes:
    primary_id:    0-20000-0-10393
    station_name:  LINDENBERG (10393-0)
    Conventions:   CF-1.7
    source:        radiosonde
    featureType:   trajectory
    history:       Created by Copernicus Early Upper Air Service Version 0, 3...
    license:       https://apps.ecmwf.int/datasets/licences/copernicus/
```
At the first look on `trajectory: 4` we know that there are 4 soundings on the requested day. A total of 227 observations are available. 

2. Request a timeseries of one station at a certain pressure level
We will ask for Lindenberg station on 500 hPa for all available dates and for air temperature and air relative humidity. This is a much larger request than the previous one.

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

