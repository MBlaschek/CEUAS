[![GitHub version](https://badge.fury.io/gh/MBlaschek%2Frasotools.svg)](https://badge.fury.io/gh/MBlaschek%2Frasotools)

# Radiosonde Tools (RASOTOOLS)

Author: M. Blaschek
Last Update: Dec 2019
Status: Development

Version: 19.12

Read and bias adjust radiosonde data using background departure statistics from a reanalysis. The development and research has been conducted during [ERACLIM](https://www.ecmwf.int/en/research/projects/era-clim) (EU 2011-2013) , [ERACLIM2](https://www.ecmwf.int/en/research/projects/era-clim2) (EU 2014-2017) and [Copernicus Climate Change Service (C3S)](https://climate.copernicus.eu/) (EU 2018-*). In close collaboration with the [ECMWF](https://www.ecmwf.int/). Employing multiple [reanalysis](https://www.ecmwf.int/en/forecasts/datasets/browse-reanalysis-datasets) for bias adjustment such as ERA5, ERA-Interim, CERA-20C, ERA-20C as well as [JRA-55](https://jra.kishou.go.jp/JRA-55/index_en.html)



<img src="https://www.ecmwf.int/sites/default/files/styles/news_item_main_image/public/project_logos/ERA-Clim-logo-transparent.png?itok=OjZfu303" alt="eraclim" width="150px" /> <img src="https://www.ecmwf.int/sites/default/files/styles/news_item_main_image/public/project_logos/ERAClim2-logo%28203px%29.png?itok=QjDC6rd-"  alt="eraclim2" width="150px" /> <img src="https://climate.copernicus.eu/themes/custom/ce/logo.svg" alt="copernicus" width="150px" /><img src="https://www.ecmwf.int/sites/default/files/ECMWF_Master_Logo_RGB_nostrap.png" alt="ecmwf" width="150px" />



*This project has received funding from the European Union’s Sixth/Seventh Framework Programme for research, technological development and demonstration under grant agreement no 265229, no 607029.*


[TOC]

# Install

The code uses Python3 as a standard. 

Dependencies:

- Bottleneck
- scipy
- Cartopy
- tqdm
- numba
- pandas
- matplotlib
- numpy
- netCDF4

Optional:

- cfgrib
- eccodes

Download the Source code from GitHub

```shell
git clone https://github.com/MBlaschek/rasotools
```

or use the package on [PyPI rasotools](https://pypi.org/rasotools)

```shell
pip install rasotools
```

# How to Use?

The module can be imported like this

```python
>>> import rasotools as rt
>>> isonde = rt.open_radiosonde('example')   # Station in Vienna
Warning different idents:  AUM00011035 > 011035
```

The `Radiosonde` class provides routines to load and save data and attach data to an object. Most Routines in `rasotools` use [xarray](http://xarray.pydata.org/en/stable/index.html) Datasets or DataArrays to work.

```python
>>> isonde
Radiosonde (AUM00011035)
Data: 
IGRAv2     : <Dataset (9 vars [date(738), pres(32)])>
Global Attributes: 
ident      : <str (011035)>
source     : <str (NOAA NCDC)>
dataset    : <str (IGRAv2)>
levels     : <str (ERA-I 32 lower)>
processed  : <str (UNIVIE, IMG)>
libs       : <str (RT(0.2) NP(1.15.4) PD(0.23.4) XR(0.11.0))>
```

This shows a record from [IGRAv2](https://www.ncdc.noaa.gov/data-access/weather-balloon/integrated-global-radiosonde-archive) (NOAA) of the Vienna radiosonde. The original data (ASCII, table format, checkout my [igra](https://github.com/MBlaschek/igra) python3 module to read these tables) has been read and interpolated to standard pressure levels where necessary. 

```python
>>> isonde.data.IGRAv2                                                      
<xarray.Dataset>
Dimensions:  (date: 738, pres: 32)
Coordinates:
  * date     (date) datetime64[ns] 2016-01-01 ... 2016-12-31T12:00:00
  * pres     (pres) float64 1e+03 2e+03 3e+03 5e+03 ... 9.5e+04 9.75e+04 1e+05
Data variables:
    gph      (date, pres) float64 ...
    temp     (date, pres) float64 ...
    rhumi    (date, pres) float64 ...
    dpd      (date, pres) float64 ...
    windd    (date, pres) float64 ...
    winds    (date, pres) float64 ...
    numlev   (date) int64 ...
    lat      (date) float64 ...
    lon      (date) float64 ...
Attributes:
    ident:      011035
    source:     NOAA NCDC
    dataset:    IGRAv2
    levels:     ERA-I 32 lower
    processed:  UNIVIE, IMG
    libs:       RT(0.2) NP(1.15.4) PD(0.23.4) XR(0.11.0)

```



## Documentation

In `doc`  a Jupyter Notebook describes some functionality and shows how to use some functions. 



# History

**rasotools** is a continuous development at the University of Vienna for research purposes.  It is based on the need to handle radiosonde data on timeseries basis and bias adjust these timeseries. Similar projects are [*metpy*](https://unidata.github.io/MetPy/latest/tutorials/upperair_soundings.html), [*skewt*](https://pypi.org/project/SkewT/), [*PyIGRA*](https://github.com/retostauffer/PyIGRA) or [*Siphon*](https://unidata.github.io/siphon/latest/examples/upperair/IGRA2_Request.html)

There is also an IGRAv2 and UADB python module [*igra*](https://github.com/MBlaschek/igra) by the same author. At a later stage it might be useful to include some of these features for plotting an calculation of profiling indices.

# License

MIT License

Copyright (c) 2019 Michael Blaschek