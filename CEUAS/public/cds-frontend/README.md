# Access to the CDS via the cds python API

Official Documentation:

Document| Abbreviation | Published
---|---|---
Product User Guide | PUG | May 2021
Algorithm Theoretical Basis Document | ATBD | June 2021

# Short Description

This directory contains a collection of Python scripts and an IPython notebook that access the Copernicus Climate Change Service (C3S) - Upper Air Service using the CDSAPI in Python.

At present, the cdsapi is directed to [sis-dev.climate.copernicus.eu]() - the developer frontend. Once the service has been approved for public use the scripts can be used also via the standard server [cds.copernicus.eu]().

The front returns files, which are either

1. Zip file containing CF 1.7 compliant netCDF4 files (one per station). This is the default.
2. Zip file containing CSV files, which can be read in directly with spreadsheet applications such as Excel. 
3. JSON file containing error messages, if an HTTP error occurs

Both file formats (csv, netCDF) can be dealt with in the CDS toolbox. 
A typical request should contain at least a `variable` and some of the other Identifiers as shown below.

| Identifier       | Request    | All possible values                                    | Type     | Explanation                                                  |
| ---------------- | ---------- | -------------------------------------------------------|--------- | ------------------------------------------------------------ |
| `variable`       | Mandatory | `[„air_temperature“, “zonal_wind“, “meridional_wind“, “wind_speed”, ”wind_direction”, "geopotential", ”air_relative_humidity”, ”air_specific_humidity”, "air_dewpoint", "dewpoint_depression"]` | String | Meteorological variables                                     |
| `country`        | Exclusive | `[“ALL”,…,”USA”]`| String|  Country codes of stations to be selected according to WMO, see examples below. Exclusive with `bbox` and `statid`|
| `bbox`           | Exclusive | `[upper,left,lower,right]`| Float or String | Boundaries of lat/lon rectangle to select stations|
| `statid`         | Exclusive | `[“SSSSS”]` | String | Allowed are: WMO or WIGOS station ID, “all”|
| `date`           | Optional | `[YYYYMMDD,YYYYMMDD]`, `YYYYMMDD`| Integer or String | List of dates of launches|
| `time`           | Optional  | `[HHMMSS,HHMMSS]` | Integer or String | List of times permitted.                                     |
| `year`           | Optional | `[YYYY,...,YYYY]`, `YYYY`| String | Years|
| `month`          | Optional | `[MM,..., MM]`, `MM`| String | Months|
| `day`            | Optional | `[DD,..., DD]`, `DD`| String | Days|
| `period`         | Optional | `[YYYYMMDD/YYYYMMDD]`| String | Start and End of a period of dates of launches|
| `optional`       | Optional |`[obs_minus_bg, ...]`            | String         | For allowed values see Optional Variables table below. |
| `pressure_level` | Optional |  `[MMMM,…,NNNNN]`, `MMMM` | Integer or String| Pressure levels in Pascal. 16 standard pressure levels (10-1000 hPa) or significant levels (if omitted) |
| `format`         | Optional  | `nc` or `csv`  | String  | Output format |
| `version`         | Optional  | `1.0` | String  | Version number (most recent version is chosen if this parameter is missing) |

### Optional Variables Table
Note: One variable must be selected above with the variables keyword. Depending on this, the optional variable then is in the same units as the main variable. 

| Optional | Description |
| --- | --- |
| `obs_minus_bg` | ERA5 observation minus first guess as stored in ERA5 ODB feedback data base|
| `obs_minus_an` | ERA5 observation minus analysis  as stored in ERA5 ODB feedback data base|
| `bias_estimate` | ERA5 variational bias estimate  as stored in ERA5 ODB feedback data base|
| `sonde_type` | Radiosonde types from S. Schroeder's VAPOR project or WMO BUFR radiosonde type codes (after 2013) as found in ERA5 ODB feedback data base|
| --- | --- |
| `desroziers_30` | Desroziers' (2005) observation uncertainty estimate, calculated from 30 day averages of ERA5 obs_minus_bg, obs_minus_an|
| `desroziers_60` | Desroziers' (2005) observation uncertainty estimate, calculated from 60 day averages of ERA5 obs_minus_bg, obs_minus_an|
| `desroziers_90` | Desroziers' (2005) observation uncertainty estimate, calculated from 90 day averages of ERA5 obs_minus_bg, obs_minus_an|
| `desroziers_180` | Desroziers' (2005) observation uncertainty estimate, calculated from 180 day averages of ERA5 obs_minus_bg, obs_minus_an|
| --- | --- |
| `RICH_bias_estimate` | Temperature Bias estimated using RICH (Haimberger et al. 2012) bias estimation method. Version indicated in variable attribute|
| `RAOBCORE_bias_estimate` | Temperature Bias estimated using RAOBCORE (Haimberger, 2007; Haimberger et al. 2012) bias estimation method. Version indicated in variable attribute|
| `RISE_bias_estimate` | Temperature Bias estimated using RICH, taking into account annual variation of solar elevation. Version indicated in variable attribute|
| `RASE_bias_estimate` | Temperature Bias estimated using RAOBCORE, taking into account annual variation of solar elevation. Version indicated in variable attribute|
| `humidity_bias_estimate` | Humidity bias estimates, using quantile matching applied to relative humidity background departures. Currently available for relative humidity, but estimates in other humidity variables will follow. Version indicated in variable attribute|
| `wind_bias_estimate` | Wind direction bias estimates using RAOBCORE (Gruber and Haimberger, 2008), available for wind_direction, u_component_of_wind and v_component_of_wind. Version indicated in variable attribute|
| --- | --- |



# Installation

On the user side, Python (advised version 3) needs to be installed. The cdsapi module can be installed using:
```python
pip install cdsapi
```
Information on how to install [pypi](https://pypi.org/project/cdsapi/) and create a `.cdsapirc` file with the valid access keys from your [cds account](https://cds.climate.copernicus.eu/api-how-to).
To use the example notebook scripts, please install the following packages:
```python
pip install numpy pandas xarray matplotlib 
```

# Usage

There are examples for different requests in the `Example.ipynb` Notebook and two simple requests are given here to show some functionality:

## Request one sounding at a specific date and at one station:
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

Finally we have two files in the directory: `example_data/2/dest_0-20000-0-10393_air_temperature.nc` and `example_data/2/dest_0-20000-0-10393_relative_humidity.nc`,
which we can read in with `xarray` and have a look at the data:
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

### Same Request but csv
A new feature is the CSV Output, so let's redo the request and set the format
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
                   'format':'csv'
               }, 
               target='download.zip')
``` 
Redo the unzipping and let's read the `csv` file with `pandas`:

```python
import pandas
tdata = pandas.read_csv('example_data/1/temperature.csv', index_col=0)
print(tdata)

          lat    lon     plev     ta                 time  trajectory_index           statid  statindex
obs_id                                                                                                 
0       52.22  14.12   4260.0  200.3  2000-02-01 05:00:00                 0  0-20000-0-10393          0
1       52.22  14.12   4260.0  200.3  2000-02-01 05:00:00                 0  0-20000-0-10393          0
2       52.22  14.12   4260.0  200.3  2000-02-01 05:00:00                 0  0-20000-0-10393          0
3       52.22  14.12   4260.0  200.3  2000-02-01 05:00:00                 0  0-20000-0-10393          0
4       52.22  14.12   4260.0  200.3  2000-02-01 05:00:00                 0  0-20000-0-10393          0
...       ...    ...      ...    ...                  ...               ...              ...        ...
9035    52.22  14.12  96100.0  278.2  2000-02-01 23:00:00                 3  0-20000-0-10393          0
9036    52.22  14.12  96100.0  278.2  2000-02-01 23:00:00                 3  0-20000-0-10393          0
9037    52.22  14.12  96100.0  278.2  2000-02-01 23:00:00                 3  0-20000-0-10393          0
9038    52.22  14.12  96100.0  278.2  2000-02-01 23:00:00                 3  0-20000-0-10393          0
9039    52.22  14.12  96100.0  278.2  2000-02-01 23:00:00                 3  0-20000-0-10393          0

[7880 rows x 8 columns]
```
Please note that there are two files, one for temperature and one for relative humdity. In order to separate stations the `statid` is present as well as the `statindex` refering to a number for each station. When multiple stations are selected, data will be appended for these stations to the corresponding variable file.

## Request a timeseries of one station at a certain pressure level
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

