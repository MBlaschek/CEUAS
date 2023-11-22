#import netCDF4
from numba import njit,prange
import numpy
import datetime
import copy
import string
#import ESMF
import math

class stationts:                           # format, scale, unit
    def __init__(self,
                 pindex,
                 statid,
                 shortname='',
                 name='',
                 startdate=19000101,
                 color='blue',
                 units='',
):
        self.statid=statid
        self.pindex=pindex
        self.startdate=startdate
        self.color=color
        self.units=units

class dailyts(stationts):
    
    def __init__(self,
                 rfpar,
                 pindex,
                 statid,
                 ens=[0],
                 dfile=[''],
                 dsuff=[''],
                 dvar='',
                 dindex=[0],
                 ddata=[0],
                 shortname='', 
                 name='', 
                 color='blue', 
                 units=''):
        
        stationts.__init__(self, pindex, statid, shortname=shortname, name=name, 
                          startdate=rfpar["startdate"], 
                          color=color, units=units)
        self.dfile=dfile

class monthlyts(stationts):
    
    def __init__(self,
                 rfpar,
                 pindex,
                 file='',
                 var='',
                 ens=[0],
                 index=[0],
                 data=[0],
                 ):
        self.file=file
    