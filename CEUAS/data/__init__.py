from . import igra

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
>>> 
"""
