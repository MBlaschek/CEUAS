from . import igra
from .chuan import chuan
"""
Module to read Data from different Archives

igra:
    IGRAv2 and UADB
    
odb_netcdf_format:
    read a netcdf converted odb file to standard radiosonde format (hour, date, pres) ???
    
"""

__doc__ = """
CEAUS Data Module
Read data from different sources:

Example IGRAv2 

Download Vienna Radiosonde to ./archive
>>> CEAUS.data.igra.download.station('AUM00011035', './archive/')

Read IGRAv2 Data to Xarray 
>>> data = CEAUS.data.igra.read.igra('AUM00011035', './archive/AUM00011035-data.txt.zip')
READ: 3724321
>>> data                                                                                                                
Out[:]: 
<xarray.Dataset>
Dimensions:  (date: 79123, pres: 16)
Coordinates:
  * date     (date) datetime64[ns] 1910-01-06T14:22:00 ... 2019-02-20T12:00:00
  * pres     (pres) float64 1e+03 2e+03 3e+03 5e+03 ... 8.5e+04 9.25e+04 1e+05
Data variables:
    gph      (date, pres) float64 nan nan nan nan ... 1.524e+03 846.0 205.0
    temp     (date, pres) float64 nan nan nan nan ... 264.2 270.0 276.3 283.8
    rhumi    (date, pres) float64 nan nan nan nan nan ... nan nan nan nan nan
    dpd      (date, pres) float64 nan nan nan nan nan ... 10.0 27.0 0.8 5.0 10.0
    windd    (date, pres) float64 nan nan nan nan ... 310.0 300.0 325.0 320.0
    winds    (date, pres) float64 nan nan nan nan nan ... 18.0 20.0 9.0 6.0 1.0
    numlev   (date) int64 5 11 9 4 4 8 7 6 12 ... 135 128 131 121 123 87 99 123
    lat      (date) float64 48.25 48.25 48.25 48.25 ... 48.25 48.25 48.25 48.25
    lon      (date) float64 16.36 16.36 16.36 16.36 ... 16.36 16.36 16.36 16.36
Attributes:
    ident:         AUM00011035
    source:        NOAA NCDC
    dataset:       IGRAv2
    processed:     UNIVIE, IMG
    interpolated:  to pres levs (#16)
"""
