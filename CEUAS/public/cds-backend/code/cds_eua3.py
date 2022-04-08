#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__version__ = '0.3'
__author__ = 'MB'
__status__ = 'dev'
__date__ = 'Thu Jul 16 09:25:53 UTC 2020'
__institute__ = 'UNIVIE'
__github__ = 'git@github.com:MBlaschek/CEUAS.git'
__doc__ = """
CDS_EUA Functions v%s
Maintained by %s at %s
Github: %s [%s]
Updated: %s
""" % (__version__, __author__, __institute__, __github__, __status__, __date__)

"""
Context
- This class CDMDataset with process_flat will replace the hug cds_eua2 function sets
- This class CDMDataset will be used in adjust and quality control
  
Performance
- HDF5 Netcdf files should be chunked or sorted by variable and a variable lie=ke recordindex , e.g. varindex could
  give a slice per Variable to speed up reading performance.
- HDF5 files could be split by variable, would improve reading performance dramatically, especially in MP
"""
import os
import glob
import time
from datetime import datetime, timedelta
import logging
import gzip
# import lz4.frame

import h5py  # most likely non standard on CDS
import numpy as np
import pandas as pd
import xarray as xr

# most likely non standard libraries on CDS
from numba import njit
import geopy
import json
import netCDF4
import multiprocessing

from functools import partial
from itertools import islice
from multiprocessing import set_start_method, Pool


# check codes from there
# https://github.com/glamod/common_data_model/blob/master/tables/observed_variable.dat
cdm_codes = {'temperature': 85, 'relative_humidity': 38, 'dew_point_temperature': 36, 'dew_point_depression': 34,
             'geopotential': 117, 'wind_direction': 106, 'wind_speed': 107, 'u_component_of_wind': 104,
             'v_component_of_wind': 105,
             'specific_humidity': 39}
#codes, to be used with old converted and harvested files (<=converted_v7)
ipar=[0]*140
ipar[0]=0
ipar[34]=34
ipar[39]=39
ipar[85]=85
ipar[106]=106
ipar[107]=107
ipar[117]=117
#ipar[]=136
ipar[36]=36 #dp
ipar[38]=38 #rh
ipar[104]=104
ipar[105]=105

cdm_codes = {'temperature': ipar[85], 'relative_humidity': ipar[38], 'dew_point_temperature': ipar[36], 
             'dew_point_depression': ipar[34],             
             'geopotential': ipar[117], 'wind_direction': ipar[106], 'wind_speed': ipar[107], 
             'u_component_of_wind': ipar[104],'v_component_of_wind': ipar[105],
             'specific_humidity': ipar[39]}
# get codes from there
# https://apps.ecmwf.int/odbgov/varno/
odb_codes = {'geopotential': 1, 'temperature': 2, 'u_component_of_wind': 3, 'v_component_of_wind': 4,
             'wind_direction': 111, 'wind_speed': 112, 'dew_point_temperature': 59, 'dew_point_depression': 299,
             'relative_humidity': 29, 'p': 999,
             'specific_humidity': 7}

std_plevs = np.asarray([10, 20, 30, 50, 70, 100, 150, 200, 250, 300, 400, 500, 700, 850, 925, 1000])

#
# This is a libary, logging must be controlled by the caller
#
logger = logging.getLogger('upperair.cdm') 

###############################################################################
#
# FUNCTIONS - NUMBA
# - calc_trajindexfast  - copy from cds_eua2 (not updated yet)
# - andisin_t           - unused, copy from cds_eua2 for times AND dates logic
# - is_sorted           - check if array is sorted
# - reverse_index       - match dates and pressures
# - andisin             - numpy.in1d
# - tohourday           - datetime to hours, days and day_shift
#
###############################################################################


#@njit(cache=True)
def searchdate(rtsidx,rtsarr,istart,istop=0):
    if istop==0:
        istop=istart
    gdict=np.zeros((rtsidx.shape[0],np.int(2)),dtype='int')
    lidx=np.zeros_like(rtsidx)
    l=0
    for i in range(len(rtsidx)-1):
        
        idx=np.searchsorted(rtsarr[rtsidx[i]:rtsidx[i+1]],(istart,istop))
        if idx[0]<idx[1]: #rtsidx[i+1]-rtsidx[i]:
#            if rtsarr[rtsidx[i]+idx]-istart<86400: 
            gdict[l]=idx
            lidx[l]=i
            l+=1
                
    return gdict[:l],lidx[:l]

#@njit
def calc_trajindexfast_rt(zidx,trajectory_index,tstart):

    for i in range(zidx.shape[0]-1):
        trajectory_index[zidx[i]-tstart:zidx[i+1]-tstart]=i

    return

def calc_trajindexfast(z, zidx, idx, trajectory_index):
    """ Calculate Trajectory Index """
    # zidx=numpy.zeros(z.shape[0],dtype=numpy.int32)
    z0 = zidx[0]
    j = 0
    l = 0
    i = 0
    for i in range(z.shape[0] - 1):
        jold = j
        #if i == 0:
            #trajectory_index[j] = l
            ##j += 1
        #else:
        nexti = i + 1
        if nexti < (z.shape[0] -1):
            while (z[nexti] == z[i]):
                nexti += 1
                if nexti==z.shape[0]:
                    break
        if nexti < (z.shape[0] -1):
            while (idx[j] >= z[i] and idx[j] < z[nexti]):
                print('trajectory_index[',j,']=', l)
                trajectory_index[j] = l
                j += 1
                if j == idx.shape[0]:
                    break
        if j > jold:
            zidx[l] = z0 + i
            l += 1
        if j == idx.shape[0]:
            break
    if j < idx.shape[0]:
        jold = j
        while (idx[j] >= z[i]): #and idx[j] < z[-1]):
            trajectory_index[j] = l
            zidx[l] = z0 + i
            j += 1
            #l += 1
            #i += 1
            if j == idx.shape[0]:
                break

           
    zidx = zidx[:l]

    return zidx#, trajectory_index

@njit(cache=True)
def calc_trajindexfastl(z, zidx, idx, trajectory_index):
    """ Calculate Trajectory Index 
Args:
    z : absolute index of all (variable)
    zidx : index of valid times in z
    idx : index of selected levels absolute
    trajectory_index : output
    """
    # zidx=numpy.zeros(z.shape[0],dtype=numpy.int32)
    z0 = 0 #zidx[0]
    j = 0
    l = 0
    i = 0
    for i in range(z.shape[0] - 1):
        jold = j
        while idx[j] >= z[i] and idx[j] < z[i + 1]:
            trajectory_index[j] = l
            j += 1
            if j == idx.shape[0]:
                break
        if j > jold:
            zidx[l] = z0 + i-1
            if zidx[l]<0:
                zidx[l]=0
            l += 1
        if j == idx.shape[0]:
            break

    if j < idx.shape[0]:
        if z.shape[0] > 1:
            i += 1
        jold = j
        while idx[j] >= z[i]:
            trajectory_index[j] = l
            j += 1
            if j == idx.shape[0]:
                break
        if j > jold:
            zidx[l] = z0 + i
            l += 1
    zidx = zidx[:l]
    return zidx

# @njit(cache=True)
# def calc_trajindexfast(z, zidx, idx, trajectory_index):
#     """ Calculate Trajectory Index 
# Args:
#     z : absolute index of all (variable)
#     zidx : index of valid times in z
#     idx : index of selected levels absolute
#     trajectory_index : output
#     """
#     # zidx=numpy.zeros(z.shape[0],dtype=numpy.int32)
#     z0 = zidx[0]
#     j = 0  # idx shape index
#     l = 0  # new unique change index -> zidx
#     i = 0
#     # Search for every recordindex
#     for i in range(z.shape[0] - 1):
#         jold = j
#         if i == 0:
#             trajectory_index[j] = l
#             j += 1
#         else:
#             nexti = i + 1
#             if nexti < (z.shape[0] -1):
#                 while (z[nexti] == z[i]):
#                     nexti += 1
#             if nexti < (z.shape[0] -1):
#                 while (idx[j] >= z[i] and idx[j] < z[nexti]):
#                     trajectory_index[j] = l
#                     j += 1
#                     if j == idx.shape[0]:
#                         break
#         # if a break as been found, write that position into zidx
#         if j > jold:
#             zidx[l] = z0 + i
#             l += 1
#         # at the end?
#         if j == idx.shape[0]:
#             break
    
#     if j < idx.shape[0]:
#         jold = j
#         while (idx[j] >= z[i]): #and idx[j] < z[-1]):
#             trajectory_index[j] = l
#             zidx[l] = z0 + i
#             j += 1
#             l += 1
#             i += 1
#             if j == idx.shape[0]:
#                 break
            
#     zidx = zidx[:l]

#     return zidx


@njit(cache=True)
def andisin_t(mask, x, v):
    """ numba version of numpy.in1d, for times"""
    jsave = 0
    for i in range(mask.shape[0]):
        if mask[i]:
            if x[i] == v[jsave]:
                mask[i] = True
            else:
                found = False
                for j in range(jsave, v.shape[0]):
                    if x[i] == v[j]:
                        found = True
                        jsave = j
                        break
                mask[i] = found
    return


@njit
def is_sorted(a) -> bool:
    """ check if an array is sorted

    Args:
        a:

    Returns:
        bool
    """
    for i in range(a.size - 1):
        if a[i + 1] < a[i]:
            return False
    return True


@njit(cache=True)
def reverse_index(mask, dates, plevs, idates, iplevs):
    """ Find all idates and iplevs in dates and plevs

    Args:
        mask: (input-output) integer index
        dates: dates in seconds
        plevs: pressure in Pa
        idates: dates in seconds
        iplevs: pressure in Pa

    Returns:
        mask

    Notes:
        Performance of numba is 10x faster than pure python.
    """
    # I assume that dates, idates and iplevs (per idate) are sorted !!!
    # plevs do not need to be sorted, but is faster
    k = 0
    for i in range(mask.shape[0]):
        # 1. Align dates (in seconds)
        if dates[i] < idates[k]:
            continue
        if dates[i] > idates[k]:
            for j in range(k, idates.shape[0]):
                if dates[i] > idates[k]:
                    k = k + 1
                else:
                    break
        # 2. Loop pressure levels
        for j in range(k, idates.shape[0]):
            if dates[i] == idates[j]:
                if plevs[i] == iplevs[j]:
                    mask[i] = j  # True
                    # print('MASK {:6} : {:6} {:4} {} == {}  {:6} == {:6}'.format(i, k, j, dates[i], idates[j], plevs[i], iplevs[j]))
                    break
                # loop j
            else:
                break
        # 3. End of dates to match
        if k >= (idates.shape[0] - 1):
            break


@njit(cache=True)
def andisin(mask, x, v):
    """ numba version of numpy.in1d(x, v)

    Args:
        mask: bool array
        x: array
        v: array

    Returns:
        mask
    """
    for i in range(mask.shape[0]):
        if mask[i]:
            found = False
            for j in range(v.shape[0]):
                if x[i] == v[j]:
                    found = True
                    break
            mask[i] = found
    return

# daysx2 puts a time array into an array of shape days x 2, where days is the number of days since 19000101
# and the second index is 0 for midnight (GMT+/- 3hrs) ascents and 1 for midday (GMT +/- 3hrs ascents).
# ascents between 3 and 9 and 15 and 21 are discarded
@njit(cache=True,boundscheck=True)
def daysx2(time,pindex,nplev,obs):
    dsecs=86400
    tofday=time%dsecs
    dindex=time//dsecs
    hindex=dindex-dindex #np.zeros(time.shape[0],dtype=np.int32)-1
    tevening=np.where(tofday>=6*dsecs//8)
    tmorning=np.where(tofday<2*dsecs//8)
    tmidday=np.where(np.logical_and(tofday>=2*dsecs/8,tofday<6*dsecs/8))
    dindex[tevening]+=1
    hindex[tmorning]=0
    hindex[tevening]=0
    hindex[tmidday]=1
    hgood=np.where(hindex>-1)
    gdays=np.unique(dindex)
    ghours=np.zeros((2,gdays.shape[0]),dtype=np.int32)
    
    rdays=np.zeros(np.max(gdays)+1,dtype=np.int32)-1
    for i in range(gdays.shape[0]): # reverse index
        rdays[gdays[i]]=i
        
    cobs=np.full((2,nplev,gdays.shape[0]),np.nan,dtype=obs.dtype)
    for i in range(time.shape[0]):
        ghours[hindex[i],rdays[dindex[i]]]=tofday[i]//3600
        cobs[hindex[i],pindex[i],rdays[dindex[i]]]=obs[i]
    
    return cobs,hindex,pindex,rdays[dindex],hgood,gdays,ghours
    

@njit(cache=True)
def tohourday(hours, days, datetimes, day_shift):
    """ convert datetimes to hours, days and day_shift

    Args:
        hours: hours from datetimes
        days: days from datetimes
        datetimes: input datetimes (in seconds)
        day_shift: 0,1, if hour like 23 should be valid for 00 next day

    Returns:
        hours, days, day_shift
    """
    ohilfold = -1
    for i in range(datetimes.shape[0]):
        if ohilfold == datetimes[i]:
            hours[i] = hours[i - 1]
            days[i] = days[i - 1]
            day_shift[i] = day_shift[i - 1]
        else:
            ohilfold = datetimes[i]
            hours[i] = ((datetimes[i] + 3599) // 3600) % 24
            if hours[i] == 0:
                day_shift[i] = datetimes[i] % 3600 != 0
            days[i] = datetimes[i] // 86400
            # if dshift[i]:
            # x=0
    return

@njit
def csvwrite(formatall,a,b,c,d,e,f,g):
    
    h=[item for sublist in zip(a,b,c,d,e,f,g) for item in sublist]
    s=formatall%tuple(h)
    
    return s

###############################################################################
#
# METADATA FUNCTIONS
# - read_standardnames      - Read CF Conventions and add some infos
# - get_attributes          - Return Metadata from read_standardnames by name, code
# - cds_to_cdm              - convert from CDS to CDM name
# - cdm_to_cds              - convert from CDM to CDS name
# - get_global_attributes   - return defined global NetCDF Attributes for this service
#
###############################################################################

def read_standardnames(url: str = None) -> dict:
    """ Read the Climate and Forcast Convention data

    Args:
        url : CF convention XML file
              default: http://cfconventions.org/Data/cf-standard-names/69/src/cf-standard-name-table.xml
    Returns:
        dict : Standard Names and Units for variables

    Notes:
            https://cfconventions.org/

    """
    import urllib.request
    import xml.etree.ElementTree as ET
    if url is None:
        url = 'http://cfconventions.org/Data/cf-standard-names/69/src/cf-standard-name-table.xml'

    response = urllib.request.urlopen(url).read()
    xmldocument = ET.fromstring(response)
    #
    # add all additional intercomparions variables
    #
    snames = ['platform_id', 'platform_name', 'latitude', 'longitude', 'time', 'air_pressure',
              'air_temperature', 'dew_point_temperature','dew_point_depression', 'relative_humidity', 'specific_humidity',
              'eastward_wind', 'northward_wind', 'wind_speed', 'wind_from_direction', 'geopotential',
              'trajectory_label', 'obs_minus_bg', 'obs_minus_an', 'bias_estimate', 'sonde_type',
              'sample_size', 'sample_error', 'report_id', 'reference_sonde_type', 
              'station_name', 
              'RISE_bias_estimate', 'RICH_bias_estimate', 'RASE_bias_estimate', 'RAOBCORE_bias_estimate', 
              'desroziers_30', 'desroziers_60', 'desroziers_90', 'desroziers_180',
              'wind_bias_estimate',
              'humidity_bias_estimate',
             ]

    cdmnames = ['header_table/primary_station_id', 'header_table/station_name', 'observations_table/latitude',
                'observations_table/longitude', 'observations_table/date_time', 'observations_table/z_coordinate']

    cdmnames += 10 * ['observations_table/observation_value']
    # todo at the moment this is hard coded here, what if JRA55 is requested?
    cdmnames += ['header_table/report_id', 'era5fb/fg_depar@body', 'era5fb/an_depar@body', 'era5fb/biascorr@body',
                 'observations_table/sensor_id',
                 'observations_table/secondary_value', 'observations_table/original_precision',
                 'observations_table/report_id', 'observations_table/reference_sensor_id',
                 'station_configuration/station_name',
                 'advanced_homogenisation/RISE_bias_estimate', 'advanced_homogenisation/RICH_bias_estimate',
                 'advanced_homogenisation/RASE_bias_estimate', 'advanced_homogenisation/RAOBCORE_bias_estimate',
                 'advanced_uncertainty/desroziers_30', 'advanced_uncertainty/desroziers_60', 'advanced_uncertainty/desroziers_90', 'advanced_uncertainty/desroziers_180',
                 'advanced_homogenisation/wind_bias_estimate', 
                 'advanced_homogenisation/humidity_bias_estimate',
                ]
    cf = {}
    for c, cdm in zip(snames, cdmnames):
        cf[c] = {'cdmname': cdm, 'units': 'NA', 'shortname': c}
        if c not in 'latitude longitude time air_pressure':
            cf[c]['coordinates'] = 'lat lon time plev'  # short names of dimensions
    l = 0
    for child in xmldocument:
        if 'id' not in child.keys():
            continue
        c = child.attrib['id']  # standard name
        if c in snames:
            i = snames.index(c)

            cf[c]['cdmname'] = cdmnames[i]
            if child[0].text is not None:
                cf[c]['units'] = child[0].text  # unit

            if child[2].text is not None:
                cf[c]['shortname'] = child[2].text  # shortname

            cf[c]['standard_name'] = c  # standard name

        l += 1
    cf['latitude']['shortname'] = 'lat'
    cf['longitude']['shortname'] = 'lon'
    cf['air_pressure']['shortname'] = 'plev'
    cf['time']['shortname'] = 'time'
    cf['bias_estimate']['cdsname'] = 'bias_estimate'
    cf['bias_estimate']['cdmcode'] = 0
    cf['bias_estimate']['odbcode'] = 0
    cf['obs_minus_bg']['cdsname'] = 'obs_minus_bg'
    cf['obs_minus_bg']['cdmcode'] = 0
    cf['obs_minus_bg']['odbcode'] = 0
    cf['obs_minus_an']['cdsname'] = 'obs_minus_an'
    cf['obs_minus_an']['cdmcode'] = 0
    cf['obs_minus_an']['odbcode'] = 0
    cf['air_temperature']['cdsname'] = 'temperature'
    cf['air_temperature']['cdmcode'] = cdm_codes['temperature']
    cf['air_temperature']['odbcode'] = 2
    cf['eastward_wind']['cdsname'] = 'u_component_of_wind'
    cf['eastward_wind']['cdmcode'] = cdm_codes['u_component_of_wind']
    cf['eastward_wind']['odbcode'] = 3
    cf['northward_wind']['cdsname'] = 'v_component_of_wind'
    cf['northward_wind']['cdmcode'] = cdm_codes['v_component_of_wind']
    cf['northward_wind']['odbcode'] = 4
    cf['wind_speed']['cdsname'] = 'wind_speed'
    cf['wind_speed']['cdmcode'] = cdm_codes['wind_speed']
    cf['wind_speed']['odbcode'] = 112
    cf['wind_from_direction']['cdsname'] = 'wind_direction'
    cf['wind_from_direction']['cdmcode'] = cdm_codes['wind_direction']
    cf['wind_from_direction']['odbcode'] = 111
    cf['relative_humidity']['cdsname'] = 'relative_humidity'
    cf['relative_humidity']['cdmcode'] = cdm_codes['relative_humidity']
    cf['relative_humidity']['odbcode'] = 29
    cf['specific_humidity']['cdsname'] = 'specific_humidity'
    cf['specific_humidity']['cdmcode'] = cdm_codes['specific_humidity']
    cf['specific_humidity']['odbcode'] = 7
    cf['dew_point_depression']['cdsname'] = 'dew_point_depression'
    cf['dew_point_depression']['cdmcode'] = cdm_codes['dew_point_depression']
    cf['dew_point_depression']['odbcode'] = 299
    cf['dew_point_temperature']['cdsname'] = 'dew_point_temperature'
    cf['dew_point_temperature']['cdmcode'] = cdm_codes['dew_point_temperature']
    cf['dew_point_temperature']['odbcode'] = 59
    cf['geopotential']['cdsname'] = 'geopotential'
    cf['geopotential']['cdmcode'] = cdm_codes['geopotential']
    cf['geopotential']['odbcode'] = 1
    cf['trajectory_label']['shortname'] = 'trajectory_label'
    cf['sonde_type']['shortname'] = 'sonde_type'
    cf['sample_size']['shortname'] = 'sample_size'
    cf['sample_error']['shortname'] = 'sample_error'
    cf['report_id']['shortname'] = 'report_id'
    cf['reference_sonde_type']['shortname'] = 'reference_sonde_type'
    cf['station_name']['shortname'] = 'station_name'
    cf['RISE_bias_estimate']['shortname'] = 'RISE_bias_estimate'
    cf['RICH_bias_estimate']['shortname'] = 'RICH_bias_estimate'
    cf['RASE_bias_estimate']['shortname'] = 'RASE_bias_estimate'
    cf['RAOBCORE_bias_estimate']['shortname'] = 'RAOBCORE_bias_estimate'
    cf['wind_bias_estimate']['shortname'] = 'wind_bias_estimate'
    cf['desroziers_30']['shortname'] = 'desroziers_30'
    cf['desroziers_60']['shortname'] = 'desroziers_60'
    cf['desroziers_90']['shortname'] = 'desroziers_90'
    cf['desroziers_180']['shortname'] = 'desroziers_180'
    cf['humidity_bias_estimate']['shortname'] = 'humidity_bias_estimate'
    return cf


def get_attributes(cdmname: str = None, cdsname: str = None, cdmcode: int = None, odbcode: int = None,
                   cf: dict = None, url: str = None, feedback: str = None) -> dict:
    """ Get attributes from CF Table for a code

    Args:
        cdmname: e.g. temperature
        cdsname: e.g. air_temperature
        cdmcode: 85, 38, 36, 34, 117, 106, 107, 104, 105, 39
        odbcode: 1, 2, 3, 4, 111, 112, 59, 299, 29, 999, 7
        cf: cf table read_standardnames()
        url: url for cf table
        feedback: feeedback variable: 'obs_minus_bg', 'obs_minus_an', 'bias_estimate'

    Returns:
        dict : Attributes

    Notes:
        ODB Codes
            geopotential : 1
            temperature : 2
            u_component_of_wind : 3
            v_component_of_wind : 4
            wind_direction : 111
            wind_speed : 112
            dew_point_temperature : 59
            dew_point_depression : 299
            relative_humidity : 29
            p : 999
            specific_humidity : 7

        CDM Codes
            temperature : 85
            relative_humidity : 38
            dew_point_temperature : 36
            dew_point_depression : 34
            geopotential : 117
            wind_direction : 106
            wind_speed : 107
            u_component_of_wind : 104
            v_component_of_wind : 105
            specific_humidity : 39
    Examples:
        >>> get_attributes(cdmname='temperature')
        {'air_temperature': {'cdmname': 'observations_table/observation_value',
          'units': 'K',
          'shortname': 'ta',
          'coordinates': 'lat lon time plev',
          'standard_name': 'air_temperature',
          'cdsname': 'temperature',
          'cdmcode': 85,
          'odbcode': 2}}
        >>> get_attributes(cdsname='air_temperature')
        {'temperature': {'cdmname': 'observations_table/observation_value',
          'units': 'K',
          'shortname': 'ta',
          'coordinates': 'lat lon time plev',
          'standard_name': 'air_temperature',
          'cdsname': 'temperature',
          'cdmcode': 85,
          'odbcode': 2}}
        >>> get_attributes(cdmcode=38)
        {'relative_humidity': {'cdmname': 'observations_table/observation_value',
          'units': '1',
          'shortname': 'hur',
          'coordinates': 'lat lon time plev',
          'standard_name': 'relative_humidity',
          'cdsname': 'relative_humidity',
          'cdmcode': 38,
          'odbcode': 29}}
    """
    if cf is None:
        cf = read_standardnames(url=url)

    if feedback is not None:
        if '@' in feedback:
            # currently the only accessable
            if feedback in ['fg_depar@body', 'an_depar@body', 'biascorr@body']:
                if 'fg_depar' in feedback:
                    feedback = 'obs_minus_bg'
                elif 'an_depar' in feedback:
                    feedback = 'obs_minus_an'
                else:
                    feedback = 'bias_estimate'
        if feedback not in cf.keys():
            return {}

    if cdmcode is not None:
        # 85 -> {air_temperature : {'units':'K', ...}}
        for ivar, iattrs in cf.items():
            if 'cdmcode' in iattrs.keys():
                if iattrs['cdmcode'] == cdmcode:
                    if feedback is None:
                        return {ivar: cf[ivar]}
                    else:
                        tmp = cf[feedback].copy()
                        tmp.update({'units': cf[ivar]['units'],
                                    'standard_name': cf[ivar]['standard_name'],
                                    'long_name': cf[feedback]['cdmname'].split('fb')[
                                                     0].upper() + ' reanalysis ' + feedback
                                    })
                        return {feedback: tmp}

    if odbcode is not None:
        # 2 -> {air_temperature : {'units':'K', ...}}
        for ivar, iattrs in cf.items():
            if 'odbcode' in iattrs.keys():
                if iattrs['odbcode'] == odbcode:
                    if feedback is None:
                        return {ivar: cf[ivar]}
                    else:
                        tmp = cf[feedback].copy()
                        tmp.update({'units': cf[ivar]['units'],
                                    'standard_name': cf[ivar]['standard_name'],
                                    'long_name': cf[feedback]['cdmname'].split('fb')[
                                                     0].upper() + ' reanalysis ' + feedback
                                    })
                        return {feedback: tmp}

    if cdmname is not None:
        cdmname = cdm_to_cds(cdmname)
        # feedback missing
        return {cdmname: cf[cdmname]}

    if cdsname is not None:
        return {cds_to_cdm(cdsname): cf[cdsname]}

    return {}


def cds_to_cdm(var, cf=None, url=None, ):
    """ CDS Name to CDM Name"""
    # air_temperature -> temperature
    if cf is None:
        cf = read_standardnames(url=url)
    if var in cf.keys():
        if 'cdsname' in cf[var].keys():
            return cf[var]['cdsname']
    return var


def cdm_to_cds(var, cf=None, url=None, ):
    """ CDM Name to CDS Name"""
    # temperature -> air_temperature
    if cf is None:
        cf = read_standardnames(url=url)
    for ivar, iattrs in cf.items():
        if 'cdsname' in iattrs.keys():
            if iattrs['cdsname'] == var:
                return ivar
    return var


def get_global_attributes(cf=None, url=None):
    # todo add more information, based on CF Table ?
    return {'Conventions': "CF-1.7", 'source': "radiosonde", 'featureType': "trajectory"}


###############################################################################
#
# BASIC HDF5 COPY Variable FUNCTION
#
###############################################################################

def do_cfcopy(fout, fin, group, idx, cf, dim0, compression, var_selection=None):
    """ Copy H5PY variables and apply subsetting (idx)

    Args:
        fout: output file
        fin: input file
        group: group
        idx: selection (mask)
        cf: cdm mapping of names
        dim0: record dimension name
        var_selection: variables

    """
    
    
    
    # cuts vars and copies attributes of observation, feedback and header tables
    tt = time.time()
    if not var_selection:
        var_selection = fin[group].keys()

    if not isinstance(var_selection, list):
        var_selection = [var_selection]

    for i in cf.keys():
        if i not in ['platform_id', 'platform_name']:
            if i in ['air_temperature', 'dew_point_temperature','dew_point_depression', 'relative_humidity', 'specific_humidity',
                     'eastward_wind', 'northward_wind', 'wind_speed', 'wind_from_direction', 'geopotential']:
                for fb in ['obs_minus_bg', 'obs_minus_an', 'bias_estimate']:
                    try:
                        cf[fb]['units'] = cf[i]['units']
                        cf[fb]['standard_name'] = i
                        cf[fb]['long_name'] = group.split('fb')[0].upper() + ' reanalysis ' + fb
                    except:
                        pass
    vlist = []
    mask=None
    for _, cfv in cf.items():
        logger.debug('CFCOPY Looking for: %s in %s', cfv['cdmname'], group)
        for v in var_selection:
            if mask is None and 'data_policy_licence' in var_selection:
                mask=fin[group]['data_policy_licence'][idx[0]:idx[-1] + 1] == 4
            if group + '/' + v == cfv['cdmname']:
                vlist.append(cfv['shortname'])
                try:
                    logger.debug('CFCOPY %s %s', v, vlist[-1])
                    if fin[group][v].ndim == 1:
#                         try:
                        if compression == 'lzf':
                            fout.create_dataset_like(vlist[-1], fin[group][v],
                                                     compression="lzf",
                                                     compression_opts = None,
                                                     shape=idx.shape,
                                                     chunks=True)
                        elif compression == 'gzip':
                            fout.create_dataset_like(vlist[-1], fin[group][v],
                                                     shape=idx.shape,
                                                     chunks=True)
                        hilf = fin[group][v][idx[0]:idx[-1] + 1]
                        if v in ['observation_value','obsvalue@body']:
                            hilf[mask]=np.nan
                        #hilf = hilf[idx[0]:idx[-1] + 1]  # use a min:max range
#                             hilf = fin[group][v][idx[0]:idx[-1] + 1]  # use a min:max range
                        if 'time' in v:
                            # convert time units
                            us = fin[group][v].attrs['units']
                            if b'hours' in us:
                                hilf = hilf * 3600  # hilf+=int(dh[0])
                            elif b'minutes' in us:
                                hilf = hilf * 60  # +int(dh[0])
                            elif b'seconds' in us:
                                hilf = hilf  # //60//60+int(dh[0])
                            elif b'days' in us:
                                hilf *= 24 * 3600
                        fout[vlist[-1]][:] = hilf[idx - idx[0]]  # but write just the subset corresponding to the variable
#                         except:
#                             logger.warning('not found: %s %s', group, v)
#                             pass
                    else:                        
                        if v == 'station_name':
                            s1 = fin[group][v].shape[1]
                            if compression == 'lzf':
                                fout.create_dataset_like(vlist[-1], fin[group][v],
                                                         compression="lzf",
                                                         compression_opts = None,
                                                         shape=(idx.shape[0], s1),
                                                         chunks=True)
                            elif compression == 'gzip':
                                fout.create_dataset_like(vlist[-1], fin[group][v],
                                                         compression=compression,
                                                         shape=(idx.shape[0], s1),
                                                         chunks=True)
                            sname = 'string{}'.format(s1)
                            if sname not in fout.keys():
                                if compression == 'lzf':
                                    fout.create_dataset(sname,
                                                        data=np.zeros(s1, dtype='S1'),
                                                        compression="lzf",
                                                        compression_opts = None,
                                                        chunks=True)
                                elif compression == 'gzip':
                                    fout.create_dataset(sname,
                                                        data=np.zeros(s1, dtype='S1'),
                                                        chunks=True)
                                fout[sname].attrs['NAME'] = np.string_(
                                    'This is a netCDF dimension but not a netCDF variable.')
                                fout[sname].make_scale(sname)
                            hilf = fin[group][v][0]
                            if hilf.shape[0] == 0:
                                print('x')
                            hilf = np.array([hilf]*idx.shape[0])
                            fout[vlist[-1]][:] = hilf
                        
                        else: 
                            s1 = fin[group][v].shape[1]
                            if compression == 'lzf':
                                fout.create_dataset_like(vlist[-1], fin[group][v],
                                                         compression="lzf",
                                                         compression_opts = None,
                                                         shape=(idx.shape[0], s1),
                                                         chunks=True)
                            elif compression == 'gzip':
                                fout.create_dataset_like(vlist[-1], fin[group][v],
                                                         shape=(idx.shape[0], s1),
                                                         chunks=True)
                            sname = 'string{}'.format(s1)
                            if sname not in fout.keys():
                                if compression == 'lzf':
                                    fout.create_dataset(sname,
                                                        data=np.zeros(s1, dtype='S1'),
                                                        compression="lzf",
                                                        compression_opts = None,
                                                        chunks=True)
                                elif compression == 'gzip':
                                    fout.create_dataset(sname,
                                                        data=np.zeros(s1, dtype='S1'),
                                                        chunks=True)
                                fout[sname].attrs['NAME'] = np.string_(
                                    'This is a netCDF dimension but not a netCDF variable.')
                                fout[sname].make_scale(sname)
                            hilf = fin[group][v][idx[0]:idx[-1] + 1, :]
                            hilf[mask, :] = np.nan
                            #hilf = hilf[idx[0]:idx[-1] + 1, :]
#                             hilf = fin[group][v][idx[0]:idx[-1] + 1, :]
                            if hilf.shape[0] == 0:
                                print('x')
                            fout[vlist[-1]][:] = hilf[idx - idx[0], :]
                            
                except Exception as e:
                    # todo fix for missing report_id SHOULD BE REMOVED
                    print(e)
                    hilf = np.zeros(shape=(idx.shape[0]), dtype='S10')
                    for i in range(hilf.shape[0]):
                        hilf[i] = '{:0>10}'.format(i)
                    fout.create_dataset(vlist[-1],
                                        data=hilf.view('S1'),
                                        shape=(idx.shape[0], 10),
                                        chunks=True)

                try:
                    for a in fin[group][v].attrs.keys():
                        if a not in ['DIMENSION_LIST', 'CLASS', 'external_table']:
                            if type(fin[group][v].attrs[a]) is str:
                                fout[vlist[-1]].attrs[a] = np.string_(fin[group][v].attrs[a])
                            else:
                                fout[vlist[-1]].attrs[a] = fin[group][v].attrs[a]

                    for a in cfv.keys():
                        if a not in ['shortname', 'odbcode', 'cdmcode']:
                            fout[vlist[-1]].attrs[a] = np.string_(cfv[a])
                        if a == 'units' and cfv[a] == 'NA':
                            fout[vlist[-1]].attrs[a] = np.string_('')
                        if a == 'units' and vlist[-1] == 'time':
                            ahilf = np.bytes_(fin[group][v].attrs[a])
                            fout[vlist[-1]].attrs[a] = ahilf
                except:
                    # quick fix should be removed
                    logger.warning('%s/%s has no attributes', group, v)
                l = 0
                for d in fout[cfv['shortname']].dims:
                    if len(d) > 0:
                        if l == 0:
                            fout[vlist[-1]].dims[l].attach_scale(fout[dim0])
                        else:
                            fout[vlist[-1]].dims[l].attach_scale(fout[sname])
                    l += 1
    tt = time.time() - tt
    if tt > 0.4:
        logger.warning('slow copy: %s %f s', group, tt)
        
def do_csvcopy(fout, fin, group, idx, cf, dim0, compression, var_selection=None):
    """ Copy H5PY variables and apply subsetting (idx)

    Args:
        fout: output file
        fin: input file
        group: group
        idx: selection (mask)
        cf: cdm mapping of names
        dim0: record dimension name
        var_selection: variables

    """
        
    # cuts vars and copies attributes of observation, feedback and header tables
    tt = time.time()
    if not var_selection:
        var_selection = fin[group].keys()

    if not isinstance(var_selection, list):
        var_selection = [var_selection]

    for i in cf.keys():
        if i not in ['platform_id', 'platform_name']:
            if i in ['air_temperature', 'dew_point_temperature','dew_point_depression', 'relative_humidity', 'specific_humidity',
                     'eastward_wind', 'northward_wind', 'wind_speed', 'wind_from_direction', 'geopotential']:
                for fb in ['obs_minus_bg', 'obs_minus_an', 'bias_estimate']:
                    try:
                        cf[fb]['units'] = cf[i]['units']
                        cf[fb]['standard_name'] = i
                        cf[fb]['long_name'] = group.split('fb')[0].upper() + ' reanalysis ' + fb
                    except:
                        pass
    vlist = []
    mask=None
    for _, cfv in cf.items():
        logger.debug('CFCOPY Looking for: %s in %s', cfv['cdmname'], group)
        for v in var_selection:
            if mask is None and 'data_policy_licence' in var_selection:
                mask=fin[group]['data_policy_licence'][idx[0]:idx[-1] + 1] == 4
            if group + '/' + v == cfv['cdmname']:
                vlist.append(cfv['shortname'])
                try:
                    logger.debug('CFCOPY %s %s', v, vlist[-1])
                    if fin[group][v].ndim == 1:
                        hilf = fin[group][v][idx[0]:idx[-1] + 1]
                        if v in ['observation_value','obsvalue@body']:
                            hilf[mask]=np.nan
                        if 'time' in v:
                            # convert time units
                            us = fin[group][v].attrs['units']
                            if b'hours' in us:
                                hilf = hilf * 3600  # hilf+=int(dh[0])
                            elif b'minutes' in us:
                                hilf = hilf * 60  # +int(dh[0])
                            elif b'seconds' in us:
                                hilf = hilf  # //60//60+int(dh[0])
                            elif b'days' in us:
                                hilf *= 24 * 3600
                        fout[v] = hilf[idx - idx[0]]  # but write just the subset corresponding to the variable
                    else:                        
                        if v == 'station_name':
                            s1 = fin[group][v].shape[1]
                            sname = 'string{}'.format(s1)
                            hilf = fin[group][v][0]
                            if hilf.shape[0] == 0:
                                print('x')
                            hilf = np.array([hilf]*idx.shape[0])
                            fout[v] = hilf
                        
                        else: 
                            s1 = fin[group][v].shape[1]
                            sname = 'string{}'.format(s1)
                            hilf = fin[group][v][idx[0]:idx[-1] + 1, :]
                            hilf[mask, :] = np.nan
                            #hilf = hilf[idx[0]:idx[-1] + 1, :]
#                             hilf = fin[group][v][idx[0]:idx[-1] + 1, :]
                            if hilf.shape[0] == 0:
                                print('x')
                            fout[v] = hilf[idx - idx[0], :]
                            
                except Exception as e:
                    # todo fix for missing report_id SHOULD BE REMOVED
                    print(e)
                    
    tt = time.time() - tt
    if tt > 0.4:
        logger.warning('slow copy: %s %f s', group, tt)


def copy_h5py_variable(fin:h5py.File, fout:h5py.File, variable:str):
    """Copy variable or group from one file to the other

    Args:
        fin ([type]): h5py File
        fout ([type]): h5py File
        variable (str): variable, e.g. /observations_table/observation_value

    Raises:
        IOError: [description]
    """
    # fin is open?
    if not (fin and fout):
        raise IOError("File status not opened. IN {} OUT {}".format(fin, fout))
    
    path = variable.split('/')
    variable = path[-1]
    current = ''
    gid = fout
    for ipath in path:
        if ipath != variable:
            # check group
            current += '/' + ipath
            try:
                gid = fout[current]
            except:
                gid = fout.create_group(current)
    # get variable
    fvar = fin[current + '/' + variable]
    if isinstance(fvar, h5py.Group):
        # recursive loop for Groups
        for ivar in fvar.keys():
            copy_h5py_variable(fin, fout, '{}/{}/{}'.format(current, variable, ivar))
        return
    
    # A Variable or Dimension ?
    if 'CLASS' in fvar.attrs.keys():
        # This is a dimension scale -> DIMENSION
        return
    
    # Create in output
    gid.create_dataset_like(variable, fvar, dtype=fvar.dtype, shape=fvar.shape, chunks=True)
    # trasfer data
    gid[variable][:] = fvar[:]
    for iatt, ival in fvar.attrs.items():
        if iatt in ['DIMENSION_LIST']:
            continue
        gid[variable].attrs[iatt] = ival
    #
    # Attach dimensions from file
    #
    try:
        idim = [ jdim[0].name for jdim in fvar.dims.keys() ]
    except:
        idim = None # this should happen for dimensions like index, string??
    
    if idim:
        for j,jdim in enumerate(idim):
            # /observations_table/index
            kdim = jdim.split('/')[-1]
            if kdim not in gid.keys():
                gid.create_dataset_like(kdim, fin[jdim])
                gid[kdim].make_scale(kdim)
                gid[variable].dims[j].attach_scale(gid[kdim])
    fout.flush() # Write

###############################################################################
#
# DATETIME FUNCTIONS
# - totimes                 - convert request time (hour) input to hour-array
# - to_seconds_since        - convert seconds with different time_unit to seconds since reference
# - seconds_to_datetime     - convert seconds to np.datetime64
# - last_day_of_month       - get last day of a month (request check)
# - convert_datetime        - convert to seconds or datetime64, from string
#
###############################################################################

def totimes(t: any) -> any:
    """ Input time to Array

    Args:
        t: e.g. '0-23', [0, 2, 6, 12]

    Returns:
        np.ndarray : list of hours
    """
    if type(t[0]) is str:
        if '-' in t[0]:
            ilist = np.array(t[0].split('-'), dtype=np.int32)
            if ilist[0] <= ilist[1]:
                ilist = np.arange(ilist[0], ilist[1] + 1)
            else:
                ilist = np.array(list(range(ilist[0], 24)) + list(range(ilist[1] + 1)), dtype=np.int32)
            out = ilist
        else:
            out = np.array(t, dtype=np.int32)
    else:
        out = np.array(t, dtype=np.int32)

    if np.min(out) < 0 or np.max(out) > 23:
        raise ValueError
    return out


def to_seconds_since(time_elapsed: int, time_unit: str, reference=None):
    """ return seconds since a Reference date applying time unit

    Args:
        time_elapsed: seconds or other time unit
        time_unit: 'hours since 1900-01-01 00:00:00'
        reference: '1900-01-01 00:00:00'

    Returns:
        int : seconds since
    """
    if reference is None:
        reference = datetime(1900, 1, 1)

    fstart = datetime.strptime(" ".join(time_unit.split(' ')[-2:]), '%Y-%m-%d %H:%M:%S')
    offset = fstart - reference
    offsets = offset.days * 24 * 3600 + offset.seconds
    fak = 1
    if 'hours' in time_unit:
        fak = 3600
    elif 'minutes' in time_unit:
        fak = 60

    secs = time_elapsed * fak + offsets
    return secs


def seconds_to_datetime(seconds, ref='1900-01-01'):
    """ from seconds to datetime64 """
    seconds = np.asarray(seconds)
    return pd.to_datetime(seconds, unit='s', origin=ref).values


def datetime_to_seconds(dates, ref='1900-01-01T00:00:00'):
    """ from datetime64 to seconds since 1900-01-01 00:00:00"""
    return ((dates - np.datetime64(ref)) / np.timedelta64(1, 's')).astype(np.int64)


# too slow
# @np.vectorize
# def seconds_to_datetime(seconds, ref=datetime(1900, 1, 1)):
#     return ref + timedelta(seconds=np.int(seconds))
# @np.vectorize
# def seconds_to_datetime(seconds, ref='1900-01-01T00:00:00'):
#     return np.datetime64(ref) + np.timedelta64(seconds, 's')


def last_day_of_month(any_day):
    """ last day of a month"""
    next_month = any_day.replace(day=28) + timedelta(days=4)  # this will never fail
    return next_month - timedelta(days=next_month.day)


# @np.vectorize
def convert_datetime(idate, return_string=False, return_seconds=False, reference=None, format=None):
    """ input date (str,int,...) to datetime64 or formatted string"""
    if isinstance(idate, str):
        # 19000101 or 01-01-1900 00:00:00
        if '-' in idate:
            if ':' in idate:
                for iformat in ['%Y-%m-%d %H:%M', '%Y-%m-%d %H:%M:%S']:
                    try:
                        idate = datetime.strptime(idate, iformat)
                        break
                    except ValueError:
                        pass
            else:
                idate = datetime.strptime(idate, '%Y-%m-%d')

    if not isinstance(idate, (datetime, np.datetime64)):
        try:
            d = int(idate)
            idate = datetime(year=d // 10000, month=d % 10000 // 100, day=d % 100)
        except ValueError:
            if (d % 100) > 28:
                idate = last_day_of_month(datetime(year=d // 10000, month=d % 10000 // 100, day=1))

    if isinstance(idate, np.datetime64):
        idate = idate.astype('M8[ms]').astype('O')

    if return_seconds:
        reference = reference if reference is not None else datetime(1900, 1, 1)
        return int((idate - reference).total_seconds())

    if return_string:
        if format is None:
            return idate.strftime('%Y%m%d')  # used by the request library
        else:
            return idate.strftime(format)
    return idate


###############################################################################
#
# Convert time, plev, obs to [time x plev] Cube
#
###############################################################################

def table_to_cube(time: np.ndarray, plev: np.ndarray, obs: np.ndarray, nplev: int = 16,
                  return_index_array: bool = False) -> tuple:
    """ Convert a ragged variable (table) to a data cube

    Args:
        time: datetime values
        plev: pressure index values
        obs: observation values
        nplev: number of pressure levels

    Returns:
        list : time indices, obs cube [time x plev]
    """
    xtime, jtime, itime = np.unique(time, return_index=True, return_inverse=True)
    data = np.full((xtime.size, nplev), np.nan, dtype=obs.dtype)
    data[itime, plev] = obs  # sortby date and plev after
    if return_index_array:
        indexes = np.full(data.shape, -1, dtype=np.int32)
        indexes[itime, plev] = np.arange(0, obs.size, dtype=np.int32)
        return jtime, indexes, data
    return jtime, data


###############################################################################
#
# MAIN FUNCTION for HUG
# - process_flat    - Take request, open file, and run read_write_request to
#                     produce CDM Frontend files. 1. Variable per request ?
#
###############################################################################

def process_flat(outputdir: str, cftable: dict, debug:bool, request_variables: dict) -> tuple:
    """ Process a station file with the requested variables

    Args:
        outputdir: output directory
        cftable: CF convention definitions table
        request_variables: request dictionary

    Returns:
        str : filename of results
        str : message or error
    """
    import os
    import glob
    # mimicks process_flat from cds_eua2
    msg = ''  # Message or error
    filename = ''  # Filename
    #print('request_variables: ', request_variables)
    try:
        if 'gridded' in request_variables:
            filename_out = outputdir + '/dest_gridded_' + str(request_variables['variable'][0]) + '.nc'
            request = request_variables
            # select via variable
            #
            # ToDo: EDIT FOR MULTIPLE VARIABLE REQUEST
            #
            
            config_file = 'hug.default.config.json'
            if os.path.isfile(config_file):
                config = json.load(open(config_file, 'r'))
            
            if 'temperature' in request['variable']:
                reqfile = config['grid_dir'] + '/CEUAS_ta_gridded.nc'
            elif 'relative_humidity' in request['variable']:
                reqfile = config['grid_dir'] + '/CEUAS_hur_gridded.nc'
            elif 'specific_humidity' in request['variable']:
                reqfile = config['grid_dir'] + '/CEUAS_hus_gridded.nc'
            elif 'wind_speed' in request['variable']:
                reqfile = config['grid_dir'] + '/CEUAS_wind_speed_gridded.nc'
            elif 'dew_point_temperature' in request['variable']:
                reqfile = config['grid_dir'] + '/CEUAS_dew_point_temperature_gridded.nc'
            
            with xr.load_dataset(reqfile) as f:
                # select via date
                print('selecting date: ', )
                if ('date' in request.keys()) and (len(request['date']) >= 1):
                    if len(request['date']) == 1:
                        odate = (request['date'][0])
                        odate = odate[:4]+'-'+odate[4:6] # +'-'+odate[6:]
#                         odate = pd.to_datetime(odate, format='%Y-%m-%d')
#                         data = f.where(f.time == odate, drop=True)
                        data = f.sel(time = odate)
                    else:
                        odate = (request['date'][0])
                        odate = odate[:4]+'-'+odate[4:6]+'-'+odate[6:]
                        edate = (request['date'][-1])
                        edate = edate[:4]+'-'+edate[4:6]+'-'+edate[6:]
                        print(odate, edate)
                        data = f.sel(time=slice(odate, edate))
                print('done: ', odate)
                        
                # select via pressure
                print('selecting date: ', )
                if ('pressure_level' in request.keys()) and (len(request['pressure_level']) > 0):
                    data =  data.where(data.pressure.isin([int(a) for a in request['pressure_level']]), drop=True)
                print('done: ', [int(a) for a in request['pressure_level']])
                    
                # select via time 
                print('selectin time: ')
                if ('time' in request.keys()) and (len(request['time']) == 1):
                    data =  data.where(data.hour == int(request['time'][0]), drop=True)
                else:
                    data =  data.where(data.hour == 12, drop=True)
                # select via coords
                print('selecting lat/lon: ')
                if len(request['gridded']) == 4 :
                    bounds = request['gridded']
                    data = data.where(data.lat >= bounds[0], drop=True).where(data.lat <= bounds[2], drop=True)
                    data = data.where(data.lon >= bounds[1], drop=True).where(data.lon <= bounds[3], drop=True)
#             except:
#                 logger.error('No gridded data available')
                print('done')
            print('cleanup output file')
#             print('rename dims')
#             data = data.rename_dims({'lat':'latitude','lon':'longitude'})
#             print('rename vars')
#             data = data.rename_vars({'ta_anomaly':'ta',}) #  'lat':'latitude', 'lon':'longitude'})
            print('drop vars')
            data = data.drop(['ta_average'])
            print('squeeze')
            data = data.squeeze(dim='hour', drop=True)
            data = data.squeeze(dim='pressure', drop=True)
            print('write to file: ', filename_out)
            data.to_netcdf(path=filename_out)
            print('done')

            
        else:
            # todo change this use the path
            statid = request_variables['statid']
#             statid = request_variables.pop('statid', None)

            if statid is None:
                logger.error('No station ID (statid) specified. %s', filename)
                raise ValueError('No station ID (statid) specified')

            filename = request_variables.pop('filename', None)
            if False:
                # old version -> not necessary anymore? will be added in check_body anyway
                if statid[:3] == '0-2':
                    suffix = ['']
                else:
                    suffix = ['0-20000-0-', '0-20300-0-', '0-20001-0-']

                for ss in suffix:
                    filename = os.path.expandvars(datadir + '/' + ss + statid + '_CEUAS_merged_v1.nc')  
                    # version as a variable
                    #filename = glob.glob(os.path.expandvars(datadir + '/' + ss + statid + '*.nc'))
                    #if len(filename) > 0:
                        #filename = filename[0]
                    #else:
                        #filename = ''

                    if os.path.isfile(filename):
                        break

            # Automaticaly adding variables for 20200- requests:
            if '0-20200-0' in statid:
                if 'optional' not in request_variables.keys():
                    request_variables['optional'] = []
                request_variables['optional'].extend(['reference_sonde_type', 'sample_size', 'sample_error']) 

            cdmnamedict = {}
            for igroup, v in cftable.items():
                if "odbcode" in v.keys():
                    cdmnamedict[v['cdsname']] = igroup

            # todo this could be changed to the cf.keys() -> cdm names of the variables
            # request_variables['variable'] is a list
            if request_variables['format'] in ['csv','fast_csv']: # in request_variables.keys():
                filename_out = outputdir + '/dest_' + statid + '_' + cdmnamedict[
                    request_variables['variable']] + '.csv.gz'
            ### for parallel writing into single file: 
#             elif request_variables['hdf']:
#                 filename_out = outputdir
            else:
                filename_out = outputdir + '/dest_' + statid + '_' + cdmnamedict[
                    request_variables['variable']] + '.nc'
            if debug: tt=time.time()
            # Make a subset of groups/variables to read (speed up)
            # Need to add station_configuration (required later) in read_write_request
            #
            if 'da' in request_variables:
                da=request_variables['da']
            else:
                da=True
            gdict = {
                'recordindices':[str(cdm_codes[request_variables['variable']]),'recordtimestamp'],
                'observations_table':['date_time','z_coordinate','observation_value','observed_variable','report_id'],
                'header_table':[],
                'station_configuration': ['station_name', 'primary_id']
            }
            if '0-20100-0' not in statid and '0-20200-0' not in statid:
                gdict["era5fb"]=[]
                gdict['advanced_uncertainty']=[]
                gdict['advanced_homogenisation']=[]
#                 gdict['advanced_homogenization']=[]
            if 'cdm' in request_variables:
                gdict['observations_table']=[]
                gdict['era5fb']=[]
                gdict['advanced_uncertainty']=[]
                gdict['advanced_homogenisation']=[]
                gdict['crs']=[]
                gdict['observed_variable']=[]
                gdict['sensor_configuration']=[]
                gdict['source_configuration']=[]
                gdict['station_configuration']=[]
                gdict['station_configuration_codes']=[]
                gdict['station_type']=[]
                gdict['units']=[]
                gdict['z_coordinate_type']=[]
                da=False
            with CDMDataset(filename=filename, groups=gdict,da=da) as data:
                if debug: print('x',time.time()-tt)
#                 try:
                #
                # parallel processing a single station file
                #
                if request_variables['single_parallel']:                
                    data.read_write_request_parallel(filename_out=filename_out,
                                            request=request_variables,
                                            cf_dict=cftable)
                else:
                    data.read_write_request(filename_out=filename_out,
                                            request=request_variables,
                                            cf_dict=cftable)
    
                
#                 except:
#                     pass
            if debug: 
                print(time.time()-tt)
                print('')

    except Exception as e:
    #except Exception as e:
        if debug:
            raise e
        logger.error('Exception %s occurred while reading %s', repr(e), filename)
        return '', 'Exception "{}" occurred while reading {}'.format(e, filename)

    return filename_out, msg


###############################################################################
#
# POST PROCESS FUNCTIONS
# - align_datetime          - Make datetime to day and standard hour [hour x day x plev]
# - stack_cube_by_hour      -
# - unstack_cube_by_hour    -
# - to_standard_launch_time -
# - level_interpolation     -
#
# Function is needed for post processing bias adjustments
#
###############################################################################


def stack_cube_by_hour(data: xr.DataArray, dim: str = 'time', hour: str = 'hour',
                       times: list = [0, 12]) -> xr.DataArray:
    """ Stack hour as a new dimension

    Args:
        data: data cube
        dim: datetime dimension
        hour: hour dimension
        times: nominal times

    Returns:
        xr.DataArray : (hour x time x plev)

    """
    # add a secret coordinate to undo stacking
    data.coords['stacking'] = (dim, np.zeros(data[dim].shape))
    data = dict(data.groupby(dim + '.hour'))
    for ikey in data.keys():
        data[ikey] = data[ikey].assign_coords(
            {dim: data[ikey][dim].to_index().to_period('D').to_timestamp().values})

    data = xr.concat(data.values(), dim=pd.Index(data.keys(), name=hour))
    data = data.reindex({hour: list(times)})
    return data


def unstack_cube_by_hour(data: xr.DataArray, dim: str = 'time', hour: str = 'hour') -> xr.DataArray:
    """ Unstack hour to time dimension

    Args:
        data: input data cube
        dim: datetime dimension
        hour: hour dimension

    Returns:
        xr.DataArray : (time x plev) cube
    """
    # use secret variable to unstack
    data = dict(data.groupby(hour))
    for ikey in data.keys():
        data[ikey] = data[ikey].assign_coords({dim: data[ikey][dim].values + np.timedelta64(ikey, 'h')})

    data = xr.concat(data.values(), dim=dim).sortby(dim)
    del data[hour]
    if 'stacking' in data.coords:
        index = (data['stacking']==0)
        data = data.isel({dim: index}).drop('stacking')

    return data


def align_datetime(data, dim: str = 'time', plev: str = 'plev', times: tuple = (0, 12), span: int = 3,
                   freq: str = '12h', apply:bool = False, **kwargs) -> xr.DataArray:
    """ Align dates and times according to times (0,12) and within a time span (e.g. 3 hours).
    Add a coordinate (e.g. standard_time) with nominal launch times according to times and span and freq.

    Args:
        data (DataArray, Dataset): Input data
        dim (str): datetime dimension
        plev (str): pressure level dimension
        times (tuple): sounding times
        span (int): plus minus times (smaller than freq/2)
        freq (str): frequency of output times
        apply (bool): apply new standard time?

    Returns:
        xr.DataArray : datetime standardized DataArray
    """
    import numpy as np
    from pandas import DatetimeIndex

    if not isinstance(data, (xr.DataArray, xr.Dataset)):
        raise ValueError('Requires a DataArray/Dataset', type(data))

    if dim not in data.dims:
        raise ValueError('Requires a datetime dimension', dim)

    if int(24 / (len(times) * 2)) < span:
        raise ValueError("Times and span do not agree!?", times, span)

    if int(24 / int(freq[:-1])) != len(times):
        raise ValueError("Times and freq do not match:", times, freq)

    if span > int(freq[:-1]) // 2:
        raise ValueError("Frequency and Span need to be consistent (span < freq/2): ", freq, span)

    dates = data[dim].values.copy()
    #
    # Convert all dates to standard_dates -> 0, 6, 12, 18 (times +/- span (3))
    #
    #_fix_datetime = np.vectorize(to_standard_launch_time)
    # newdates = _fix_datetime(dates, span=span)  # (time: 33%)
    newdates = to_standard_launch_time(dates, times, span=span)
    idx_std = DatetimeIndex(newdates).hour.isin(times)
    resolution = np.zeros(newdates.size)  # same size as dates
    #
    # check for duplicates in standard launch times
    #
    u, c = np.unique(newdates[idx_std], return_counts=True)
    conflicts = u[c > 1]
    if conflicts.size > 0:
        counts = _count_data(data, dim=dim, plev=plev)
        logger.warning("Conflicts %s [+/- %d]: %d (# %d) of %d", str(times), span, conflicts.size, c[c > 1].sum(),
                       idx_std.sum())
        for i in conflicts:
            indices = np.where(newdates == i)[0]  # numbers  (time: 45%)
            #
            # Count available data (DataArray or Dataset)
            #
            # slow
            # counts = data.isel(**{dim: indices}).count(plev).values
            # counts = _count_data(data.isel(**{dim: indices}), dim=dim, plev=plev)
            # slow end
            icounts = counts[indices]
            #
            # offsets to standard launch time
            #
            offset = np.abs((dates[indices] - i) / np.timedelta64(1, 'h'))
            j = np.argsort(offset)  # sort time offsets (first we want)
            jmax = np.argmax(icounts[j])  # check if counts from other time is larger
            if jmax != 0:
                #
                # there is a sounding with more level data (+/- 1 hour)
                #
                if (offset[j][0] + 1) <= offset[j][jmax]:
                    # ok close enough
                    jj = j.copy()
                    jj[j == 0] = jmax  # first pos is now at the position of the maximum
                    jj[j == jmax] = 0  # maximum is now first
                    j = jj
                #
                # there is a sounding with + 2 more levels
                #
                elif (icounts[j][0] + 2) <= icounts[j][jmax]:
                    # a lot more
                    jj = j.copy()
                    jj[j == 0] = jmax  # first pos is now at the position of the maximum
                    jj[j == jmax] = 0  # maximum is now first
                    j = jj
                else:
                    pass  # keep time sorting

            for m, k in enumerate(offset[j]):
                if m == 0:
                    continue  # this is the minimum

                if k == 0:
                    # same date/hour, but less data
                    newdates[indices[j][m]] += np.timedelta64(1, 'h')
                    resolution[indices[j][m]] = -1  # add hour
                else:
                    resolution[indices[j][m]] = -2  # revert back

    #
    # Select standard times (input) and resolved conflicts
    #
    idx_std = DatetimeIndex(newdates).hour.isin(times) & (resolution >= 0)  # this is the selection index
    resolution[idx_std] = 1
    u, c = np.unique(newdates[idx_std], return_counts=True)  # check only standard times
    conflicts = u[c > 1]
    if conflicts.size > 0:
        logger.warning("Conflicts remain: %d  Std: %d New: %d", conflicts.size, idx_std.sum(), newdates.size)
        if kwargs.get('debug', False):
            raise RuntimeError("Duplicates")

    nn = (resolution < 0).sum()  # number of resolved conflicts
    nx = (~idx_std).sum()  # number of non-standard times
    #
    # Standard datetime
    #
    jname = 'standard_{}'.format(dim)
    data[jname] = (dim, newdates)  # new coordinate
    data.attrs['info'] = 'Selected hours: %s' % str(times)
    data[jname].attrs.update({'long_name': 'nominal hour', 'units': 'h', 'info': str(times)})
    data[jname + '_flag'] = (dim, resolution.astype(int))
    data[jname + '_flag'].attrs.update({'units': '1', 'standard_name': 'flag_conflict_resolution',
                                        'info': '0: no conflict, 1: standard time, -1: datetime duplicate, less data, -2: lesser candidate'})

    logger.info("Modified dates: %d, Standard %s: %d / %d / %d", nn, str(times), idx_std.sum(), nx, newdates.size)
    if not all(data[jname].values == np.sort(data[jname].values)):
        logger.warning("Values are not sorted by %s", jname)
        # data = data.sortby(jname)
        pass
        
    if apply:
        data = data.swap_dims({dim: jname}).rename({dim: '{}_reported'.format(dim)}).rename({jname : dim})
        data[dim].attrs['info'] = jname
        data['{}_reported'.format(dim)].attrs['info'] = dim
    return data


def _count_data(data: xr.DataArray, dim: str = 'time', plev: str = 'plev') -> np.ndarray:
    """ Helper for align_datetime"""
    #
    # Count data per pressure level (if it is a dimension)
    #
    if plev in data.dims:
        #
        # Array format
        #
        if isinstance(data, xr.DataArray):
            return data.count(plev).values
        else:
            return data.count(plev).to_dataframe().sum(axis=1).values  # sum across variables

    elif data[dim].to_index().is_unique:
        #
        # has not pressure levels
        #
        if isinstance(data, xr.DataArray):
            return data.count(dim).values
        else:
            return data.count(dim).to_dataframe().sum(axis=1).values
    else:
        #
        # Table format
        #
        return data.groupby(dim).count().to_dataframe().max(axis=1).values


def to_standard_launch_time(dates: np.datetime64, std_times: list, span: int =3, debug: bool = False):
    """ Convert datetime64 to standard launch times

    Args:
        dates: input datetime64
        std_times: list of standard launch times, e.g. 0, 12
        span: allowed time departure, e.g. 3 for 3 hours (+/-)
        debug: debug flag
    """
    dateshour = dates.astype(np.datetime64(1, 'h'))  # truncate to hour
    hour = pd.DatetimeIndex(dateshour).hour.astype(int) # only the hour
    tcorr = np.zeros(hour.size, dtype=np.int)
    for itime in std_times:
        # +- (span) - > itime
        diff = (hour - itime)
        # within add the inverse, everywhere else 0
        tcorr = np.where(np.abs(diff) <= span, -1*diff, tcorr)
        if itime == 0:
            # check for values before midnight?
            tcorr = np.where(np.abs(diff -24) <= span, -1*(diff-24), tcorr)
    newdates = dateshour + (tcorr + np.timedelta64(0,'h'))
    newdates = np.where(tcorr == 0, dateshour, newdates) # fillback
    return newdates


def level_interpolation(idata: xr.DataArray, dim: str = 'time', method: str = 'linear', fill_value=None,
                        extrapolate: bool = False, **kwargs) -> xr.DataArray:
    """ Interpolate CDM variable per dimension (time, plev)

    Args:
        idata: Inputdata
        dim: coordinate of input, e.g. time, plev
        method: interpolation method, e.g. linear, log-linear for plev
        fill_value: Interpolation fill_value: 'extrapolate', None, np.nan
        extrapolate: fill missing values [False]

    Returns:
        xr.DataArray : Interpolate data, same as input

    Examples:
        Interpolate adjustments per level backwards and forward in time
        >>> data['hur_q'].groupby('plev').apply(level_interpolation)

        Interpolate adjustments per profile up and down in pressure using a log-linear interpolation
        >>> data['hur_q'].groupby('time').apply(level_interpolation, dim='plev', method='log-linear')

        Interpolate relative humidity per profile up and down in pressure using a log-linear interpolation, but no
        extrapolation
        >>> data['hur'].groupby('time').apply(level_interpolation, dim='plev', method='log-linear', extrapolate=False)

    """
    if not isinstance(idata, xr.DataArray):
        raise ValueError("Requires a DataArray, not ", type(idata))
    #
    # maybe check if dimensions are present
    #
    if idata.isnull().all().item():
        return idata

    idim = idata.dims[0]  # first dimension
    obs = idata[idim].values.copy()  # copy Dimension values
    #
    # swap first to interpolation dimension (dim)
    #
    idata = idata.swap_dims({idim: dim})
    #
    # Find duplicated values /not allowed in interpolation
    #
    itx = idata[dim].to_index()
    ittx = itx.duplicated(keep=False) & ~np.isfinite(idata.values)  # BUG duplicates in dim & not a value -> skip
    ittx = np.where(~np.isfinite(itx), True, ittx)  # BUG missing values can be in the coordinate (e.g. plev)
    #
    # some duplicates might remain
    #
    if itx[~ittx].duplicated().any():
        ittx[~ittx] = itx[~ittx].duplicated()  # duplicates with same values / should not happen

    # message(itx.size, ittx.sum(), mname='DUPLICATES', **update_kw('level', 1, **kwargs))
    # logger.debug('Int %s %d > %d', dim, itx.size, ittx.sum())
    if method == 'log-linear':
        #
        # Pressure levels to log
        #
        idata.values[~ittx] = idata[~ittx] \
            .assign_coords({'plev': np.log(idata['plev'].values[~ittx])}) \
            .interpolate_na(dim, method='linear', fill_value=fill_value) \
            .values
    else:
        idata[~ittx] = idata[~ittx].interpolate_na(dim, method=method, fill_value=fill_value)
    #
    # Extrapolate in dimension (backward and forward) or (up and down)
    #
    if extrapolate:
        idata = idata.bfill(dim).ffill(dim)
    #
    # revert to initial dimension
    #
    idata[idim] = (dim, obs)  # add coordinate / dimension
    # todo figure out if drop(obs) is requiredx
    return idata.swap_dims({dim: idim}).drop_vars('obs', errors = 'ignore')  # swap dimension back, remove obs coordinate


###############################################################################
#
# Request Functions for initiating CDMDataset
#
###############################################################################


def cds_request_wrapper(request: dict, request_filename: str = None, cds_dataset: str = None, cds_url: str = None,
                        overwrite: bool = False):
    """ Run a CDS Request and return a CDMDataset or a CDMDatasetList if more variables are requested

    Args:
        request: dictionary with CDS request
        request_filename: Output filename (zip)
        cds_dataset: CDS Dataset (insitu-comprehensive-upper-air-observation-network)
        cds_url: CDS Server URL (https://sis-dev.climate.copernicus.eu/api/v2)
        overwrite: flag if caching is disabled

    Returns:
        CDMDataset
        CDMDatasetList

    Examples:
            Request for Station 10393 Lindenberg the variable air_temperature for every available time and pressure
            level. The request will be forwarded to the CDS and to the VM, downloaded and loaded into the CDMDataset
            class.
            >>> data = cds_request_wrapper({'statid':'10393', 'variable':'air_temperature'})

    """
    import zlib
    import zipfile
    import cdsapi
    import urllib3
    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
    try:
        if request_filename is None:
            request_filename = '{}.zip'.format(zlib.adler32(bytes(repr(request), 'utf-8')))

        if not os.path.isfile(request_filename) or overwrite:
            client = cdsapi.Client(
                'https://sis-dev.climate.copernicus.eu/api/v2' if cds_url is None else cds_url)
            client.retrieve(
                cds_dataset if cds_dataset is not None else 'insitu-comprehensive-upper-air-observation-network',
                request,
                request_filename)
            # file is downloaded to cds_outputfile
        else:
            logger.info('Requested file exists: %s', request_filename)

        idir = os.path.dirname(request_filename) if '/' in request_filename else '.'
        os.makedirs(idir, exist_ok=True)
        with zipfile.ZipFile(request_filename, 'r') as f:
            files = f.namelist()
            f.extractall(idir + '/')
            for ifile in files:
                logger.debug('Extracting %s/%s', idir, ifile)
        if 'format' in request:
            if request['format'] == 'csv':
                # todo read multiple csv files and return a DataFrame
                raise NotImplementedError()
                
                
        files = ["{}/{}".format(idir, ifile) for ifile in files]
        if len(files) > 1:
            return CDMDatasetList(*files)
        return CDMDataset(filename=files[0])

    except Exception as e:
        logger.error('CDSAPI Request failed %s', str(request))
        raise e


def vm_request_wrapper(request: dict, request_filename: str = None, vm_url: str = None, overwrite: bool = False):
    """ Run a VM Request and return a CDMDataset or a CDMDatasetList if more variables are requested

    Args:
        request: dictionary with VM request
        request_filename: Output filename (zip)
        vm_url: VM Server URL (http://early-upper-air.copernicus-climate.eu)
        overwrite: flag if caching is disabled

    Returns:
        CDMDataset
        CDMDatasetList

    Examples:
            Request for Station 10393 Lindenberg the variable air_temperature for every available time and pressure
            level. The request will be forwarded to the CDS and to the VM, downloaded and loaded into the CDMDataset
            class.
            The same request as with the CDSAPI, but with the VM.
            >>> data = vm_request_wrapper({'statid':'10393', 'variable': ['temperature']})

    """
    import zlib
    import zipfile
    import requests
    import urllib3
    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
    try:
        if request_filename is None:
            if not 'hdf' in request:
                request_filename = '{}.zip'.format(zlib.adler32(bytes(repr(request), 'utf-8')))
            else:
                request_filename = 'out.nc'
        
        if not os.path.isfile(request_filename) or overwrite:
            r = requests.post('http://early-upper-air.copernicus-climate.eu' if vm_url is None else vm_url,
                              headers={'content-type': 'application/json'},
                              json=request,
                              stream=True)
            if r.status_code != requests.codes.ok:
                raise RuntimeError(r.text)

            with open(request_filename, 'wb') as f:
                f.write(r.content)

        idir = os.path.dirname(request_filename) if '/' in request_filename else '.'
        os.makedirs(idir, exist_ok=True)
        if not 'hdf' in request:
            with zipfile.ZipFile(request_filename, 'r') as f:
                files = f.namelist()
                f.extractall(idir + '/')
                for ifile in files:
                    logger.debug('Extracting %s/%s', idir, ifile)

#             if 'format' in request:
#                 if request['format'] == 'csv':
#                     # todo read multiple csv files and return a DataFrame
#                     raise NotImplementedError()

            files = ["{}/{}".format(idir, ifile) for ifile in files]
            if 'speed_test' in request:
                    return files[0]
            else:
                if len(files) > 1:
                    return CDMDatasetList(*files)
                if files[0][-4:] == '.csv':
                    return pd.read_csv(files[0])
                elif files[0][-7:] == '.csv.gz':
                    return pd.read_csv(files[0], compression='gzip')
                else:
                    return CDMDataset(filename=files[0])
        else:
            return request_filename

    except Exception as e:
        logger.error('VM Request failed %s', str(request))
        raise e


###############################################################################
#
# CDM Backend / Frontend Classesx
# 1. CDMVariable
#       * origin -> Link to H5PY file location
#       * Attributes, e.g.: shape : shape of H5PY Variable
# 2. CDMGroup
#       just a wrapper class for different printing
# 3. CDMDataset
#       Main class with H5PY File handlers
#       * file  : H5PY file handler (open)
#       * filename : string name of the file that is open
#       * groups : list of CDMGroups and CDMVariables
#       * hasgroups : True CDM Backend file / False CDM Frontend file
#                     depending on the file different variables are present
###############################################################################


class CDMVariable:
    """ CDMVariable(origin, name, ...)"""

    def __init__(self, filelink, name: str, **kwargs):
        self._name = name
        self.link = filelink
        self._data = filelink
        self._names = []
        for ikey, ival in kwargs.items():
            setattr(self, ikey, ival)
            self._names.append(ikey)

    def __setitem__(self, key, value):
        self._names.append(key)
        self.__setattr__(key, value)

    def __getitem__(self, item):
        return self._data[item]

    def __delitem__(self, item):
        # used for deleting attributes?
        if hasattr(self, item):
            self.__delattr__(item)  # private call to remove from Class
            self._names.remove(item)  # remove from list

    def __repr__(self):
        return self._name + " : " + " ".join(["{}".format(str(getattr(self, i))) for i in self._names])

    def keys(self):
        """ list of Attributes
        """
        return self._names

    def isArray(self):
        """ Loaded to Memory?
        """
        return False if 'HDF5' in str(self._data) else True

    def update(self, link=None, data=None):
        """ Update HDF5 File Reference, Link"""
        if link is not None:
            self.link = link

        if not self.isArray():
            self._data = link

        if data is not None:
            self._data = data

    def get_size(self):
        if hasattr(self.link.id, 'get_storage_size'):
            return self.link.id.get_storage_size()

class CDMGroup(CDMVariable):
    """ CDMGroup(origin, name, ...)"""
    shape = 0  # for visual representation, groups do not have a shape in h5py

    def __repr__(self):
        text = ''
        for i in self.keys():
            if i in ['shape']:
                continue
            istatus = 'L' if self[i].isArray() else ' '
            text += "\n{:_<50} :{:1}: {}".format(i, istatus, self[i].shape)
        return self._name + ":\n" + text

    def __getitem__(self, item):
        return self.__getattribute__(item)

    def __delitem__(self, item):
        # delete variable from Group
        if hasattr(self, item):
            delattr(self, item)  # remove from class
            self._names.remove(item) # remove from internal list
            del self.link[item]  # call h5py removal
            self.shape = len(self.keys())

class CDMDatasetList(dict):
    """ This is the dict subclass for CDMDataset with multiple variables
    """
    def __init__(self, *args):
        data = {}
        for ifile in args:
            name = (ifile.split('/')[-1]).replace('.nc','')
            data[name] = CDMDataset(filename=ifile)
        super().__init__(data)
        self.__dict__ = self
    def __repr__(self):
        text = "CDMDatasetList <{}>".format(len(self.keys()))
        for ikey in self.keys():
            text += "\n- {}\n  {}".format(ikey, "\n  ".join(str(self[ikey]).splitlines()))
            
        return text


def parallel_writing(filename, dims, trajectory_index, idx, zidx, compression, filename_out, request, globatts, cfcopy):
    time0 = time.time()
#     cfcopy, add_aux_vars = parallel_input
    key, value = list(cfcopy.items())[0]
    filename_out = filename_out[:-3] + '_' + key + '.nc'
    with CDMDataset(filename) as infile:
        with h5py.File(filename_out, 'w') as fout:
            #
            # Dimensions (obs, trajectory)
            #
#             if add_aux_vars == 0:
#                 for d, v in dims.items():
#                     fout.create_dataset(d, data=v)
#                     fout[d].attrs['NAME'] = np.string_('This is a netCDF dimension but not a netCDF variable.')
#                     fout[d].make_scale(d)  # resolves phony_dim problem
            
            if value == 0:
#                 for d, v in dims.items()[:1]:
                for d, v in islice(dims.items(), 0, 1):
                    fout.create_dataset(d, data=v)
                    fout[d].attrs['NAME'] = np.string_('This is a netCDF dimension but not a netCDF variable.')
                    fout[d].make_scale(d)  # resolves phony_dim problem
            elif value == 1:
#                 for d, v in dims.items()[1:]:
                for d, v in islice(dims.items(), 1, None):
                    fout.create_dataset(d, data=v)
                    fout[d].attrs['NAME'] = np.string_('This is a netCDF dimension but not a netCDF variable.')
                    fout[d].make_scale(d)  # resolves phony_dim problem

            elif value == 2:
                fout.create_dataset('trajectory_index', data=trajectory_index)
                fout['trajectory_index'].attrs['long_name'] = np.string_(
                    "index of trajectory this obs belongs to")
                fout['trajectory_index'].attrs['instance_dimension'] = np.string_("trajectory")
                fout['trajectory_index'].attrs['coordinates'] = np.string_("lat lon time plev")

                print('attrs written')
                
            elif value == 3:
                #
                # Station Configuration
                #
                # station_configuration
                if 'station_configuration' in infile.groups:
                    igroup = 'station_configuration'
                    cfcstationcon = {'station_name': 
                                     {
                                         'cdmname': 'station_configuration/station_name',
                                         'units': 'NA',
                                         'shortname': 'station_id',
                                         'coordinates': 'lat lon time plev',
                                         'standard_name': 'station_name'
                                     }
                                    } 
                    do_cfcopy(fout, infile.file, igroup, idx, cfcstationcon, 'obs', compression,
                              var_selection=['station_name'])
                    logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)

            else:           

                #
                # Variables based on cfcopy
                #
                if 'observations_table' in infile.groups:
                    igroup = 'observations_table'
                    do_cfcopy(fout, infile.file, igroup, idx, cfcopy, 'obs', compression,
                              var_selection=['observation_id', 'latitude', 'longitude', 'z_coordinate',
                                             'observation_value', 'date_time', 'sensor_id', 'secondary_value',
                                             'original_precision', 'reference_sensor_id', 'report_id','data_policy_licence']) #+cdm_obstab)
                    # 'observed_variable','units'
                    logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)
                    if 'nodims' in request.keys():
                        fout['time'].make_scale('time')
                #
                # Feedback Information
                #
                if 'era5fb' in infile.groups:
                    igroup = 'era5fb'
                    try:
                        do_cfcopy(fout, infile.file, igroup, idx, cfcopy, 'obs', compression,
                                  var_selection=['fg_depar@body', 'an_depar@body',
                                                 'biascorr@body']) #+cdm_eratab)
                        # ['vertco_reference_1@body','obsvalue@body','fg_depar@body'])
                        logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)
                    except KeyError as e:
                        print(e)
    #                     raise KeyError('{} not found in {} {}'.format(str(e), str(request['optional']), infile.name))
                #
                # Adjusted Values with ERA5
                # TODO: Names not set, group missing
                if 'adjera5' in infile.groups:
                    igroup = 'adjera5'
                    try:
                        do_cfcopy(fout, infile.file, igroup, idx, cfcopy, 'obs', compression,
                                  var_selection=['bias_estimate', 'bias_estimation_method'])
                        logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)
                    except KeyError as e:
                        print(e)
    #                     raise KeyError('{} not found in {} {}'.format(str(e), str(request['optional']), infile.name))

                if 'advanced_homogenisation' in infile.groups :
                    igroup = 'advanced_homogenisation'
                    varsel=[]
                    try:
                        for o in request['optional']:
                            if o in varseldict[request['variable']]:
                                varsel.append(o)
                    except:
                        pass

                    if varsel:       
                        try:
                            do_cfcopy(fout, infile.file, igroup, idx, cfcopy, 'obs', compression,
                                      var_selection=varsel)
                            logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)
                        except KeyError as e:
                            print(e)
                            #raise KeyError('{} not found in {} {}'.format(str(e), str(request['optional']), infile.name))


                #
                # advanced_uncertainty
                # 
                if 'advanced_uncertainty' in infile.groups:
                    print('advanced_uncertainty in infile.groups')
                    igroup = 'advanced_uncertainty'
                    try:
                        do_cfcopy(fout, infile.file, igroup, idx, cfcopy, 'obs', compression,
                                  var_selection=['desroziers_30', 'desroziers_60', 'desroziers_90', 'desroziers_180'])
                        logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)
                    except KeyError as e:
                        print(e)
    #                     raise KeyError('{} not found in {} {}'.format(str(e), str(request['optional']), infile.name))

                #
                # Header Information
                #
                if 'header_table' in infile.groups:
                    igroup = 'header_table'
                    # only records fitting criteria (zidx) are copied
                    # todo why is lon, lat not here?
                    do_cfcopy(fout, infile.file, igroup, zidx, cfcopy, 'trajectory', compression,
                              var_selection=['report_id'])
                    logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)
                    # ,'station_name','primary_station_id'])
                    # todo could be read from the observations_table


                #
                # Fix Attributes and Globals
                #
                if key == 'trajectory_label':
                    fout['trajectory_label'].attrs['cf_role'] = np.string_('trajectory_id')
                    fout['trajectory_label'].attrs['long_name'] = np.string_('Label of trajectory')
                for a, v in globatts.items():
                    fout.attrs[a] = np.string_(v)

                fout.attrs['history'] = np.string_(
                    'Created by Copernicus Early Upper Air Service Version 0, ' + datetime.now().strftime(
                        "%d-%b-%Y %H:%M:%S"))
                fout.attrs['license'] = np.string_('https://apps.ecmwf.int/datasets/licences/copernicus/')

                for i in fout.keys():
                    #print(i)
                    if 'nodims' in request.keys():
                        if (i == 'obs' or i == 'trajectory' or 'string' in i):
                            fout.__delitem__(i)
                    version = ''
                    if '.' in i:
                        vers = i.find('.')
                        version = i[(vers-2):(vers+2)] # '_x.x'
                        newname = i.replace(version, '')
                        fout[newname] = fout[i]
                        fout.__delitem__(i)
                        fout[newname].attrs['version'] = np.string_(version[1:]) # 'x.x'

                # removing variables with restricted access
                try: fout['observations_table'].__delitem__('observation_value')
                except: pass
                try: fout['era5fb'].__delitem__('obsvalue@body')
                except: pass

class CDMDataset:
    """ This is the main CDM Class for handling CDM files, both frontend and backend
    """

    # memory efficient, no duplicates
    # __slots__ = ['filename', 'file', 'groups', 'data']

    def __init__(self, filename: str = None, groups: dict = None, mode:str = 'r',da=False):
        """ Init Class CDMDataset with a filename, cds_request or vm_request

        Args:
            filename: path of NetCDF HDF5 backend of frontend file
            groups: read only these groups
            mode: h5py reading mode, e.g. r, w, r+, a

        Examples:
            Open a CDM Backend file
            >>> data = CDMDataset(filename='0-20000-0-10393_CEUAS_merged_v0.nc')

            Open a CDM frontend file
            >>> data = CDMDataset(filename='dest_0-20000-0-01001_air_temperature.nc')

        """
        if filename is None:
            raise ValueError('Specifiy either filename.')
 
        self.da=da
        if filename == 'empty':
            self.filename = ''
            self.name = ''
            self.file = None
            self.hasgroups = False
            self.groups = []
        else:
            self.filename = filename
            self.name = filename.split('/')[-1]  # just the name as in rname
            self.file = h5py.File(filename, mode)
            logger.debug("[OPEN] %s", self.filename)
            self.hasgroups = False
            self.groups = []
            if not da:
                self.inquire(groupdict=groups)  # Get variables and Groups
            else:
                self.groups=list(self.file.keys())

    def __getitem__(self, item):
        return self.__getattribute__(item)

    def __delitem__(self, item):
        # delete group or CDMFrontend variable
        if hasattr(self, item):
            delattr(self, item)  # remove it from class
            self.groups.remove(item) # remove from groups
            del self.file[item]  # remove from file
            self.file.flush()  # update file

    def __repr__(self):
        text = "File: " + repr(self.file)
#         text += "\nPermissions: %s%s (%s)" % ('R' if os.access(self.filename, os.R_OK) else '', 'W' if os.access(self.filename, os.W_OK) else '', os.getlogin())
        text += "\nFilesize: %.2f MB" % (self.file.id.get_filesize()/1024/1024)
        text += "\nFilename: " + self.filename
        text += "\n(G)roups/(V)ariables: \n"
        for igroup in self.groups:
            istatus = 'L' if self[igroup].isArray() else ' '
            text += "\n - {} | {:_<45} :{:1}: {}".format('G' if isinstance(self[igroup], CDMGroup) else 'V',
                                                         igroup, istatus, self[igroup].shape)
        return text

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def inquire(self, groupdict: dict = None):
        """ Read HDF5 file structures
        """
        groups = list(self.file.keys())
        if groupdict is not None:
            groups = [igroup for igroup in groups if igroup in groupdict.keys()]
        try:
            for igroup in groups:
                # Group or Variable
                new = False
                if igroup not in self.groups:
                    self.groups.append(igroup)
                    new = True

                if isinstance(self.file[igroup], h5py.Group):
                    if new:
                        jgroup = CDMGroup(self.file[igroup], igroup)
                    else:
                        jgroup = getattr(self, igroup)  # Get CDMGroup

                    jgroup.update(link=self.file[igroup])  # reconnect to Group, e.g. if reopened
                    varkeys=list(self.file[igroup].keys())

                    if groupdict is not None:
                        varkeys = groupdict[igroup] if len(groupdict[igroup]) > 0 else list(self.file[igroup].keys()) #varkeys
                    
                    for ivar in varkeys:
                        try:
                            if ivar not in jgroup.keys():
                                shape = self.file[igroup][ivar].shape
                                jgroup[ivar] = CDMVariable(self.file[igroup][ivar], ivar, shape=shape)
                            jgroup[ivar].update(link=self.file[igroup][ivar])
                        except:
                            # can fail if variable not present
                            pass
                    jgroup.shape = len(jgroup.keys())  # update never hurts
                    setattr(self, igroup, jgroup)  # attach to class
                    self.hasgroups = True
                else:
                    if new:
                        setattr(self, igroup, CDMVariable(self.file[igroup], igroup, shape=self.file[igroup].shape))
                    self[igroup].update(link=self.file[igroup])
                

        except Exception as e:
            logger.debug(repr(e))
            self.close()

    def close(self):
        """ Close H5py file
        """
        self.file.close()
        logger.debug("[CLOSED] %s", self.filename)

    def copy(self, filename:str, force:bool = False, variables:list = None, groups:list = None, fuzzy:bool=False):
        """Copy whole file (via shell) or variables or groups (via h5py)

        Args:
            filename (str): output filename
            force (bool, optional): Overwrite output file? Defaults to False.
            variables (list, optional): list of variables to copy. e.g. /observations_table/report_id. Defaults to None.
            groups (list, optional): list of groups to copy. Defaults to None.
            fuzzy (bool, optional): Search fuzzy for variable names. Defaults to False.

        Raises:
            OSError: [description]
            RuntimeError: [description]
        """
        import os
        import shutil

        if os.path.isfile(filename) and not force:
            raise OSError("File exists", filename)
        
        if filename == self.filename:
            raise RuntimeError("Filename is identical with current?")
        
        if variables is None and groups is None:
            shutil.copy(self.filename, filename)
            logger.info("copy complete [%s] to %s", self.name, filename)
        else:
            logger.warning("This is experimental. I hope you know waht you are doing.")
            if isinstance(variables, str):
                variables = [variables]
            if isinstance(groups, str):
                groups = [groups]
            
            with h5py.File(filename, 'w') as f:
                for ivar in variables:
                    copy_h5py_variable(self.file, f, ivar)
                for igroup in groups:
                    copy_h5py_variable(self.file, f, igroup)
            logger.info("Copied variables [%d] to %s", len(variables), filename)


    def reopen(self, mode: str = 'r', tmpdir: str='~/tmp', write_to_filename: str = None, strict:bool = False, **kwargs):
        """ Reopen the HDF5 file with different mode

        Args:
            mode: r, r+, w, w+, a
            tmpdir: temporary directory for writing
            write_to_filename: new filename (a copy of the object)
            strict: raise an Error if reopening fails (might not have been closed before)
            **kwargs: other options to h5py.File( )
        """
        import os
        import shutil

        self.file.close()

        if write_to_filename is not None:
            # cp file to that name and open with mode
            shutil.copy(self.filename, write_to_filename)
            logger.warning("reopen [%s] Copy file to %s", mode, write_to_filename)
            self.filename = write_to_filename
            self.name = self.filename.split('/')[-1]
            self.file = h5py.File(self.filename, mode=mode, **kwargs)

        else:
            try:
                self.file = h5py.File(self.filename, mode=mode, **kwargs)
            except OSError as e:
                if strict:
                    logger.error('reopen with', mode, self.filename)
                    raise e

                if '~' in tmpdir:
                    tmpdir = os.path.expanduser(tmpdir)

                if '$' in tmpdir:
                    tmpdir = os.path.expandvars(tmpdir)

                os.makedirs(tmpdir, exist_ok=True)
                if os.path.isdir(tmpdir):
                    shutil.copy(self.filename, '{}/{}'.format(tmpdir, self.name))
                    logger.warning("reopen [%s] Copy file to %s/%s", mode, tmpdir, self.name)
                    self.filename = '{}/{}'.format(tmpdir, self.name)
                    self.name = self.filename.split('/')[-1]
                    self.file = h5py.File(self.filename, mode=mode, **kwargs)
                else:
                    raise OSError('reopen with', mode, self.filename)
        logger.debug("reopen %s [%s]", self.filename, mode)
        self.inquire()  # Need to update links
    
    def storage_info(self, group:str=None, variable:str=None, limit:float = 0.1):
        # get file size and check 
        ffree = self.file.id.get_freespace() # space due to deleting
        fsize = self.file.id.get_filesize()
        text = "File: " + repr(self.file)
        text += "\nFilesize: %.2f MB" % (fsize/1024/1024)
        text += "\nFile freespace: %.2f MB (might need to repack)" % (ffree/1024/1024)
        text += "\nFilename: " + self.filename
        text += "\n(G)roups/(V)ariables: \n"
        if group is None:
            group = self.groups

        estimates = []
        textmsg = []
        for igroup in group:
            for ivar in self[igroup].keys():
                try:
                    if variable:
                        if ivar != variable:
                            continue
                    current = self[igroup][ivar].get_size()/1024/1024
                    if current < limit:
                        continue
                    shape = self[igroup][ivar].shape
                    estimates.append(current) 
                    textmsg.append("{:7.2f} MB| {:_<60} {}".format(current, igroup +'/'+ str(ivar), shape))
                except:
                    pass
        # sort
        estimates = np.array(estimates)
        textmsg = np.array(textmsg)
        idx = np.argsort(estimates)[::-1]  # reverse
        text += "\n".join(textmsg[idx].tolist())
        print(text)

    def select_variables(self, names: list, groups:list=None, fuzzy:bool = False):
        """Generate a list of existing variables based on names and groups
        Args:
            names: list of variables
            groups: list of groups
            fuzzy: partial match
        
        Returns:
            list: list of variables / groups
              (variable, group) or (variable)
        """
        if isinstance(names, str):
            names = [names]

        found = []
        if self.hasgroups:
            if groups is None:
                groups = self.groups

            for igroup in groups:
                for ivar in names:
                    if ivar in self[igroup].keys():
                        found.append((ivar, igroup))
                    else:
                        if fuzzy:
                            if any([ivar in i for i in self[igroup].keys()]):
                                found.append((ivar, igroup))
        else:
            for ivar in names:
                if ivar in self.keys():
                    found.append(ivar)
                else:
                    if fuzzy:
                        if any([ivar in i for i in self.keys()]):
                            found.append(ivar)
        return found

    def read_attributes(self, name: str, group: str = None, subset: list = None):
        """ Return all attributes or a subset for a variable

        Args:
            name: variable
            group: h5py group
            subset: list of attributes

        Returns:
            dict : Attributes of Variable
        """
        if group is not None:
            if not self.hasgroups:
                raise ValueError('CDM Frontend files do not have groups. Remove group keyword!', group)
            if group not in self.groups:
                raise ValueError('Group not found', group)
            if name not in self[group].keys():
                raise ValueError('Variable not found', name)
            name = "{}/{}".format(group, name)  # HDF5 access

        if subset is not None:
            if not isinstance(subset, list):
                subset = [subset]

        attributes = {}
        # can raise a KeyError if name not in file
        for ikey, ival in self.file[name].attrs.items():
            if subset is not None and ikey in subset:
                try:
                    attributes[ikey] = ival.decode()
                except AttributeError:
                    attributes[ikey] = ival
            else:
                if ikey not in ['DIMENSION_LIST']:
                    try:
                        attributes[ikey] = ival.decode()
                    except AttributeError:
                        attributes[ikey] = ival
        return attributes

    def load_variable_from_file(self, name, group: str = None, return_data: bool = False,
                                decode_byte_array: bool = True) -> list:
        """ Allow to load a variable from a group

        Args:
            name (str, list): name of the variable
            group (str): group
            return_data (bool): read data and return?
            decode_byte_array (bool): decode string arrays

        Returns:
            data: list of data variables requested
        """
        if group is not None:
            if group not in self.groups:
                raise ValueError('Group not found', group)

            if isinstance(name, str):
                name = [name]

            for iname in name:
                if iname not in self[group].keys():
                    raise ValueError('Variable not found', group, iname)

                if self[group][iname].isArray():
                    continue

                # TODO: use read_direct(array, source_sel=None, dest_sel=None) from h5py
                # according to documentation, this avoids internal copies ?
                if len(self[group][iname].shape) > 1 and decode_byte_array:
                    # concat byte-char-array
                    data = self[group][iname][()].astype(object).sum(axis=1).astype(str)
                else:
                    data = self[group][iname][()]  # get the full numpy array

                self[group][iname].update(data=data)
                # setattr(self[group], iname, CDMVariable(data, iname, shape=data.shape))

            if return_data:
                return [self[group][iname][()] for iname in name]
        else:
            if isinstance(name, str):
                name = [name]
            for iname in name:
                if self[iname].isArray():
                    continue
                if len(self[iname].shape) > 1 and decode_byte_array:
                    data = self[iname][()].astype(object).sum(axis=1).astype(str)
                else:
                    data = self[iname][()]  # Read data
                # setattr(self, iname, CDMVariable(data, iname, shape=data.shape))
                self[iname].update(data=data)
            if return_data:
                return [self[iname][()] for iname in name]

    def availability(self, varnum: int = None, date_time_name: str = 'date_time'):
        """ Read datetime and check how much data is available
        Returns:
            pd.DataFrame or pd.Series: counts per profile
        """
        if self.hasgroups:
            if 'recordtimestamp' in self.groups:
                timestamp = self.load_variable_from_file('recordtimestamp', return_data=True)[0]
            elif isinstance(self['recordindices'], CDMGroup):
                timestamp = self.load_variable_from_file('recordtimestamp', group='recordindices', return_data=True)[0]
            else:
                timestamp = self.load_variable_from_file(date_time_name, return_data=True)[0]

            if 'recordindex' in self.groups:
                recindex = self.load_variable_from_file('recordindex', return_data=True)[0]
                num_rec = np.diff(recindex).tolist() + [self['observations_table']['date_time'].shape[0] - recindex[-1]]
                return pd.Series(data=num_rec, index=seconds_to_datetime(timestamp), name='num_obs')
            elif 'recordindices' in self.groups:
                num_rec = {}
                for ikey in self['recordindices'].keys():
                    try:
                        ikey = int(ikey)
                        if varnum:
                            if ikey not in varnum:
                                continue
                        recindex = self.load_variable_from_file(str(ikey), group='recordindices', return_data=True)[0]
                        num_rec[ikey] = np.diff(recindex).tolist() + [0]
                    except:
                        pass
                return pd.DataFrame(data=num_rec, index=seconds_to_datetime(timestamp))
            else:
                num_rec = np.ones(timestamp.shape[0])
                return pd.Series(data=num_rec, index=seconds_to_datetime(timestamp), name='num_obs')
        else:
            timestamp = self.load_variable_from_file('time', return_data=True)[0]
            timestamp, num_rec = np.unique(timestamp, return_counts=True)
            return pd.Series(data=num_rec, index=seconds_to_datetime(timestamp), name='num_obs')
 

    def read_write_request(self, filename_out: str, request: dict, cf_dict: dict):
        """ This is the basic request used in the cds_eua2 script

        Args:
            filename_out: request output filename, e.g. : /data/public//tmp//049532808458/dest_0-20000-0-10393_air_temperature.nc
            request: Request dictionary, minimum: {'variable' : 'temperature'}
            cf_dict: CF Convention for Names (results from read_standardnames() Function)

        Examples:
            >>> data = CDMDataset(filename)
            >>> data.read_write_request('test_cdm_frontend.nc', {'statid': '01001', 'variable': ['temperature']})
            Write all temperature data to a CDM frontend file
        """
        # version of process_flat with Leo's optimizations
        if 'variable' not in request.keys():
            logger.error('No variable specified %s %s', str(request.keys()), self.name)
            raise ValueError('No variable specified')

        variable = request['variable']
        compression = request['compression']
        if isinstance(variable, list):
            if len(variable) > 1:
                # todo potentially launch here recursive ?
                raise RuntimeError('Requests need to be split up by variable')
            variable = variable[0]  # assumes a list ?

        if cf_dict is None:
            cf_dict = read_standardnames()

        cdsname = cdm_to_cds(variable, cf=cf_dict)
        cdmattrs = cf_dict[cdsname]  #
        cdmnum = cdmattrs['cdmcode']  # 85 for air_temperature
        time0 = time.time()
        #
        # todo future version will not used observed variable anymore / change here
        # HDF5 file work fast with slices
        # - trange      - datetime slice
        # - mask        - logical array applied after trange for time, plev, variable
        #

        trange, mask = self.read_observed_variable(cdmnum,
                                                   variable='observation_value',
                                                   dates=request.get('date', None),
                                                   plevs=request.get('pressure_level', None),
                                                   times=request.get('time', None),
                                                   observed_variable_name='observed_variable',
                                                   date_time_name='date_time',
                                                   date_time_in_seconds=False,
                                                   z_coordinate_name='z_coordinate',
                                                   group='observations_table',
                                                   dimgroup='observations_table',
                                                   return_index=True,
                                                   rtsindex=request.get('rtsidx',None)
                                                   )
        logger.debug('Datetime selection: %d - %d [%5.2f s] %s', trange.start,
                     trange.stop, time.time() - time0, self.name)
        tt=time.time() - time0
        print(tt)
        idx = np.where(mask)[0] + trange.start  # absolute integer index
        if len(idx) == 0:
            logger.warning('No matching data found %s', self.name)
            raise ValueError('No matching data found')  # add CDMname for logging
            #return

        logger.debug('Data found: %d %s', len(idx), self.name)
        #
        # Make Trajectory Information (lon, lat, profile id, ...)
        #
        trajectory_index = np.zeros_like(idx, dtype=np.int32)
        trajectory_index2 = np.zeros_like(idx, dtype=np.int32)
        if 'recordinindex' in self.groups:
            # unsorted indices in root
            recordindex = self.file['recordindex'][()]
        else:
            # sorted indices are in recordindices group / by variable
            recordindex = self.file['recordindices'][str(cdmnum)][()]  # values
            
        zidx = np.where(np.logical_and(recordindex >= trange.start, recordindex <= trange.stop))[0]
        #zidx2=np.zeros_like(zidx)
            
        #zidx2 = recordindex[zidx]
        #recordindex=np.unique(recordindex)
        #zidx=np.unique(recordindex)
        
        #rt=self.observations_table.date_time[trange][idx-trange.start]
        #for i in range(rt.shape[0]-1):
            #if rt[i+1]>rt[i]:
                #if i+1+trange.start not in recordindex:
                    #print(i)
        
        
        #rr=self.observations_table.report_id[trange][idx-trange.start]
        #rp=self.observations_table.z_coordinate[trange][idx-trange.start]
        #l=0
        #for a,b,c in zip(rt,rr,rp):
            #print(a,b.view('S11'),c)
            #l+=1
            #if l>100:
                #break
        
        
        #calc_trajindexfast_rt(recordindex,trajectory_index,trange.start)
        zidx = calc_trajindexfastl(recordindex, zidx, idx, trajectory_index)
        #
        # Dimensions and Global Attributes
        #
        tt=time.time() - time0
        print(tt)
        dims = {'obs': np.zeros(idx.shape[0], dtype=np.int32),
                'trajectory': np.zeros(zidx.shape[0], dtype=np.int32)}
        globatts = get_global_attributes()  # could put more infors there ?
        #
        # Definition of Variables to write
        #
        # Common Variables needed for a requested file
        #
        snames = ['platform_id', 'platform_name', 'observation_value', 'latitude',
                  'longitude', 'time', 'air_pressure', 'trajectory_label', 
                  'report_id', 'station_id']
        varseldict={}
        varseldict['temperature']=['RAOBCORE_bias_estimate', 'RASE_bias_estimate', 'RICH_bias_estimate', 'RISE_bias_estimate']
        varseldict['wind_direction']=['wind_bias_estimate']
        varseldict['u_component_of_wind']=['wind_bias_estimate']
        varseldict['v_component_of_wind']=['wind_bias_estimate']
        varseldict['relative_humidity']=['humidity_bias_estimate']
        # added report_id -> in observations_table, not to be confused with report_id from header_table -> trajectory_label
        logger.debug('Request-keys: %s', str(request.keys()))
        snames.append(cdsname)  # Add requested variable
        #
        # Add Optional Variables
        # Check read_standard_names for correct names otherwise it will be rejected
        #
        if 'optional' in request.keys():
            snames.extend(request['optional'])
            
        if 'bias_estimate_method' in request.keys():
            # use era5adj group
            pass
        
        # Copy Metadata -> used by do_cfcopy
        cfcopy = {}  # Copy of CDM Info Dictionary (read_standard_names())
        for ss in snames:
            try:
                cfcopy[ss] = cf_dict[ss]
            except:
                pass
        #
        # End Definition of Variables to write
        #
        logger.debug('Writing: %s', filename_out)
        tt=time.time() - time0
        print(tt)
        
#         rstcd = self.file['observations_table']['data_policy_licence'][idx[0]:idx[-1] + 1] == 4
#         zrstcd=recordindex[zidx+1]-idx[0]
#         rstcd=None
        
        if request['format'] in ['csv', 'fast_csv']:
            fout={}
            if 'observations_table' in self.groups:
                igroup = 'observations_table'
                do_csvcopy(fout, self.file, igroup, idx, cfcopy, 'obs', compression,
                          var_selection=['observation_id', 'latitude', 'longitude', 'z_coordinate',
                                         'observation_value', 'date_time', 'sensor_id', 'secondary_value',
                                         'original_precision', 'reference_sensor_id', 'report_id','data_policy_licence'])
                # 'observed_variable','units'
                logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)
            #
            # Feedback Information
            #
            if 'era5fb' in self.groups:
                igroup = 'era5fb'
                try:
                    do_csvcopy(fout, self.file, igroup, idx, cfcopy, 'obs', compression,
                              var_selection=['fg_depar@body', 'an_depar@body',
                                             'biascorr@body'])
                    # ['vertco_reference_1@body','obsvalue@body','fg_depar@body'])
                    logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)
                except KeyError as e:
                    raise KeyError('{} not found in {} {}'.format(str(e), str(request['optional']), self.name))
            #
            # Adjusted Values with ERA5
            # TODO: Names not set, group missing
            if 'adjera5' in self.groups:
                igroup = 'adjera5'
                try:
                    do_csvcopy(fout, self.file, igroup, idx, cfcopy, 'obs', compression,
                              var_selection=['bias_estimate', 'bias_estimation_method'])
                    logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)
                except KeyError as e:
                    raise KeyError('{} not found in {} {}'.format(str(e), str(request['optional']), self.name))

            if 'advanced_homogenisation' in self.groups :
                igroup = 'advanced_homogenisation'
                varsel=[]
                try:
                    for o in request['optional']:
                        if o in varseldict[request['variable']]:
                            varsel.append(o)
                except:
                    pass

                if varsel:       
                    try:
                        do_csvcopy(fout, self.file, igroup, idx, cfcopy, 'obs', compression,
                                  var_selection=varsel)
                        logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)
                    except KeyError as e:
                        print(e)
                        #raise KeyError('{} not found in {} {}'.format(str(e), str(request['optional']), self.name))
            #
            # advanced_uncertainty
            # 
            if 'advanced_uncertainty' in self.groups:
                print('advanced_uncertainty in self.groups')
                igroup = 'advanced_uncertainty'
                try:
                    do_csvcopy(fout, self.file, igroup, idx, cfcopy, 'obs', compression,
                              var_selection=['desroziers_30', 'desroziers_60', 'desroziers_90', 'desroziers_180'])
                    logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)
                except KeyError as e:
                    raise KeyError('{} not found in {} {}'.format(str(e), str(request['optional']), self.name))

            if request['format'] in ['fast_csv']:
                dellist = []
                for i in fout:
                    if(len(np.shape(fout[i])) > 1):
                        dellist.append(i)

                for i in dellist:
                    del fout[i]
    
                dtype = dict(names = list(fout.keys()), formats=[])
                lfout=[]
                for k in dtype['names']:
                    if len(fout[k].shape)==1:
                        pass #lfout.append(fout[k])
                    else:
                        fout[k]=fout[k].view('|S{}'.format(fout[k].shape[1])).flatten().astype(str) 
                    dtype['formats'].append(fout[k].dtype)
                headstr = ''
                formatstr = ''
                for n,d in zip(dtype['names'],dtype['formats']):
                    headstr = headstr+n+','
                    sd=str(d)
                    if 'int' in sd:
                        formatstr = formatstr+'%.0i,'
                    elif 'float' in sd:
                        formatstr = formatstr+'%.6f,'
                    else:
                        formatstr = formatstr+'"%.'+sd[2:]+'s",'
                        #formatstr = formatstr+'%.'+''+'s,'

                print(time.time()-time0)
                formatstrn=formatstr[:-1]+'\n'
                formatall=formatstrn*fout[dtype['names'][0]].shape[0]
#                formatall=formatstrn*lfout[0].shape[0]
#                 with open(filename_out,'w') as f:
#                 with lz4.frame.open(filename_out,'wt') as f:
                if filename_out is not None:
                    with gzip.open(filename_out,'wt',compresslevel=1) as f:
                        f.write(headstr[:-1]+'\n')
                        b=[item for sublist in zip(*fout.values()) for item in sublist]
                        f.write(formatall%tuple(b))
                else:
                    #b=[item for sublist in zip(*fout.values()) for item in sublist]
                    return headstr[:-1]+'\n'#+formatall%tuple(item for sublist in zip(*fout.values()) for item in sublist)
    
            else:
                #
                # Station Configuration
                #
                if 'station_configuration' in self.groups:
                    igroup = 'station_configuration'
                    cfcstationcon = {'station_name': 
                                     {
                                         'cdmname': 'station_configuration/station_name',
                                         'units': 'NA',
                                         'shortname': 'station_id',
                                         'coordinates': 'lat lon time plev',
                                         'standard_name': 'station_name'
                                     }
                                    } 
                    do_csvcopy(fout, self.file, igroup, idx, cfcstationcon, 'obs', compression,
                              var_selection=['station_name'])
                    logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)
                
                dtype = dict(names = list(fout.keys()), formats=[])
                lfout=[]
                print(fout)
                for k in dtype['names']:
                    if len(fout[k].shape)==1:
                        pass #lfout.append(fout[k])
                    else:
                        fout[k]=fout[k].view('|S{}'.format(fout[k].shape[1])).flatten().astype(str) 
                    dtype['formats'].append(fout[k].dtype)
                headstr = ''
                formatstr = ''
                for n,d in zip(dtype['names'],dtype['formats']):
                    headstr = headstr+n+','
                    sd=str(d)
                    if 'int' in sd:
                        formatstr = formatstr+'%.0i,'
                    elif 'float' in sd:
                        formatstr = formatstr+'%.6f,'
                    else:
                        formatstr = formatstr+'"%.'+sd[2:]+'s",'
                        #formatstr = formatstr+'%.'+''+'s,'

                print(time.time()-time0)
                formatstrn=formatstr[:-1]+'\n'
                formatall=formatstrn*fout[dtype['names'][0]].shape[0]
#                formatall=formatstrn*lfout[0].shape[0]
#                 with open(filename_out,'w') as f:
#                 with lz4.frame.open(filename_out,'wt') as f:

                tt=time.time()
                ss=[]
                for v in fout.values():
                    ss.append(f"{tuple(v)}".split(','))
                x=[','.join(t) for t in zip(*ss)]
                sss='\n'.join(x)
                print(len(sss),time.time()-tt)
                tt=time.time()
                b=[item for sublist in zip(*fout.values()) for item in sublist]
                s=formatall%tuple(b)
                print(len(sss),time.time()-tt)
                
#                s=csvwrite(formatall,*fout.values())

                if filename_out is not None:
                    with gzip.open(filename_out,'wt',compresslevel=1) as f:
                        f.write(headstr[:-1]+'\n')
                        b=[item for sublist in zip(*fout.values()) for item in sublist]
                        f.write(formatall%tuple(b))
                else:
                    #b=[item for sublist in zip(*fout.values()) for item in sublist]
                    return headstr[:-1]+'\n'#+formatall%tuple(item for sublist in zip(*fout.values()) for item in sublist)
        else:
            with h5py.File(filename_out, 'w') as fout:
                # todo future -> this could be replaced by a self.write_to_frontend_file(filename_out, )
                #
                # chose wich variable will be selected and copied into the outputfile 
                # from all the sources: observations_table, era5fb, header_table,...
                #
                # Dimensions (obs, trajectory)
                #
                for d, v in dims.items():
                    fout.create_dataset(d, data=v)
                    fout[d].attrs['NAME'] = np.string_('This is a netCDF dimension but not a netCDF variable.')
                    fout[d].make_scale(d)  # resolves phony_dim problem

                fout.create_dataset('trajectory_index', data=trajectory_index)
                fout['trajectory_index'].attrs['long_name'] = np.string_(
                    "index of trajectory this obs belongs to")
                fout['trajectory_index'].attrs['instance_dimension'] = np.string_("trajectory")
                fout['trajectory_index'].attrs['coordinates'] = np.string_("lat lon time plev")

                print('attrs written')

                #
                # Adding CDM
                #
                cdm_obstab = []
                cdm_eratab = []
                cdmlist = request.get('cdm', None)
                if cdmlist != None:
                    # removing variables with restricted access
                    if ('era5fb/obsvalue@body' in cdmlist) or ('observations_table/observation_value' in cdmlist):
                        try: cdmlist.remove('era5fb/obsvalue@body')
                        except: pass
                        try: cdmlist.remove('observations_table/observation_value')
                        except: pass
                    for i in cdmlist:
                        grp = i.split('/')[0]
                        if i == 'era5fb':
                            logger.debug('Full er5fb CDM-request disabled.')
                        elif grp == 'era5fb':
                            try: 
                                var = i.split('/')[1]
                                cdm_eratab.append(var)
                                cfcopy[var]={'cdmname': i, 'units': '', 'shortname': var, 'coordinates': 'lat lon time plev', 'standard_name': var, 'cdsname': var,}
                            except: pass
                            cdmlist.remove(i)

                        if i == 'observations_table':
                            logger.debug('Full observations_table CDM-request disabled.')
                        elif grp == 'observations_table':
                            try: 
                                var = i.split('/')[1]
                                cdm_obstab.append(var)
                                cfcopy[var]={'cdmname': i, 'units': '', 'shortname': var, 'coordinates': 'lat lon time plev', 'standard_name': var, 'cdsname': var,}
                            except: pass
                            cdmlist.remove(i)
                            # NOW USE THOSE for the do_cfcopy below
                    logger.debug('CDM - adding groups and variables: %s', str(cdmlist))
                    for cdmstring in cdmlist:
                        print(cdmstring)
                        cdmsplit = cdmstring.split('/')
    #                     try:
                        print(self.groups)
                        if cdmsplit[0] in self.groups:
                            print('group available')
                            # TODO: Add attributes and descriptions?
                            # whole group
                            if len(cdmsplit) == 1:
                                print('group copy')
                                fout.create_group(cdmsplit[0])
                                for cdmvar in self[cdmsplit[0]].keys():
                                    try:
                                        fout[cdmsplit[0]].create_dataset(cdmvar, data=self[cdmsplit[0]][cdmvar][:])
                                    except:
                                        pass
                            # single var of group 
                            if len(cdmsplit) == 2:
                                print('single var copy')
                                try:
                                    fout.create_group(cdmsplit[0])
                                except:
                                    pass # group alread exists?
                                try:
                                    fout[cdmsplit[0]].create_dataset(cdmsplit[1], data=self[cdmsplit[0]][cdmsplit[1]][:])
                                except:
                                    pass

                #
                # Variables based on cfcopy
                #
                if 'observations_table' in self.groups:
                    igroup = 'observations_table'
                    do_cfcopy(fout, self.file, igroup, idx, cfcopy, 'obs', compression,
                              var_selection=['observation_id', 'latitude', 'longitude', 'z_coordinate',
                                             'observation_value', 'date_time', 'sensor_id', 'secondary_value',
                                             'original_precision', 'reference_sensor_id', 'report_id','data_policy_licence']+cdm_obstab)
                    # 'observed_variable','units'
                    logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)
                    if 'nodims' in request.keys():
                        fout['time'].make_scale('time')
                #
                # Feedback Information
                #
                if 'era5fb' in self.groups:
                    igroup = 'era5fb'
                    try:
                        do_cfcopy(fout, self.file, igroup, idx, cfcopy, 'obs', compression,
                                  var_selection=['fg_depar@body', 'an_depar@body',
                                                 'biascorr@body']+cdm_eratab)
                        # ['vertco_reference_1@body','obsvalue@body','fg_depar@body'])
                        logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)
                    except KeyError as e:
                        raise KeyError('{} not found in {} {}'.format(str(e), str(request['optional']), self.name))
                #
                # Adjusted Values with ERA5
                # TODO: Names not set, group missing
                if 'adjera5' in self.groups:
                    igroup = 'adjera5'
                    try:
                        do_cfcopy(fout, self.file, igroup, idx, cfcopy, 'obs', compression,
                                  var_selection=['bias_estimate', 'bias_estimation_method'])
                        logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)
                    except KeyError as e:
                        raise KeyError('{} not found in {} {}'.format(str(e), str(request['optional']), self.name))

                if 'advanced_homogenisation' in self.groups :
                    igroup = 'advanced_homogenisation'
                    varsel=[]
                    try:
                        for o in request['optional']:
                            if o in varseldict[request['variable']]:
                                varsel.append(o)
                    except:
                        pass

                    if varsel:       
                        try:
                            do_cfcopy(fout, self.file, igroup, idx, cfcopy, 'obs', compression,
                                      var_selection=varsel)
                            logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)
                        except KeyError as e:
                            print(e)
                            #raise KeyError('{} not found in {} {}'.format(str(e), str(request['optional']), self.name))


                #
                # advanced_uncertainty
                # 
                if 'advanced_uncertainty' in self.groups:
                    print('advanced_uncertainty in self.groups')
                    igroup = 'advanced_uncertainty'
                    try:
                        do_cfcopy(fout, self.file, igroup, idx, cfcopy, 'obs', compression,
                                  var_selection=['desroziers_30', 'desroziers_60', 'desroziers_90', 'desroziers_180'])
                        logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)
                    except KeyError as e:
                        raise KeyError('{} not found in {} {}'.format(str(e), str(request['optional']), self.name))

                #
                # Header Information
                #
                if 'header_table' in self.groups:
                    igroup = 'header_table'
                    # only records fitting criteria (zidx) are copied
                    # todo why is lon, lat not here?
                    do_cfcopy(fout, self.file, igroup, zidx, cfcopy, 'trajectory', compression,
                              var_selection=['report_id'])
                    logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)
                    # ,'station_name','primary_station_id'])
                    # todo could be read from the observations_table
                #
                # Station Configuration
                #
                # station_configuration
                if 'station_configuration' in self.groups:
                    igroup = 'station_configuration'
                    cfcstationcon = {'station_name': 
                                     {
                                         'cdmname': 'station_configuration/station_name',
                                         'units': 'NA',
                                         'shortname': 'station_id',
                                         'coordinates': 'lat lon time plev',
                                         'standard_name': 'station_name'
                                     }
                                    } 
                    do_cfcopy(fout, self.file, igroup, idx, cfcstationcon, 'obs', compression,
                              var_selection=['station_name'])
                    logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)

                #
                # Fix Attributes and Globals
                #
                fout['trajectory_label'].attrs['cf_role'] = np.string_('trajectory_id')
                fout['trajectory_label'].attrs['long_name'] = np.string_('Label of trajectory')
                for a, v in globatts.items():
                    fout.attrs[a] = np.string_(v)

                fout.attrs['history'] = np.string_(
                    'Created by Copernicus Early Upper Air Service Version 0, ' + datetime.now().strftime(
                        "%d-%b-%Y %H:%M:%S"))
                fout.attrs['license'] = np.string_('https://apps.ecmwf.int/datasets/licences/copernicus/')


                for i in fout.keys():    
                    if 'toolbox' in request.keys() and not 'optional' in request.keys():
    #                     if i in ['wind_from_direction']:
    #                         fout['ta'] = fout[i]
    #                         fout.__delitem__(i)
                        print('no optionals')
                    elif 'toolbox' in request.keys():
                        if i in ['ta', 'hur', 'ua', 'va', 'wind_from_direction']:
                            oldunits = fout[i].attrs['units']
                            fout.__delitem__(i)
                            oldkey=(request['optional'][0])
                            fout[i]=fout[oldkey]
                            fout[i].attrs['units'] = oldunits
                            fout.__delitem__(oldkey)


                for i in fout.keys():
                    #print(i)
                    if 'nodims' in request.keys():
                        if (i == 'obs' or i == 'trajectory' or 'string' in i):
                            fout.__delitem__(i)
                    version = ''
                    if '.' in i:
                        vers = i.find('.')
                        version = i[(vers-2):(vers+2)] # '_x.x'
                        newname = i.replace(version, '')
                        fout[newname] = fout[i]
                        fout.__delitem__(i)
                        fout[newname].attrs['version'] = np.string_(version[1:]) # 'x.x'

                # removing variables with restricted access
                try: fout['observations_table'].__delitem__('observation_value')
                except: pass
                try: fout['era5fb'].__delitem__('obsvalue@body')
                except: pass

#             while True:
#                 try:
#                     with h5py.File(filename_out, 'a') as fout_base:
#                         if request['hdf']:
#                             fout = fout_base[request['statid'].split('-')[-1]]
#                         else:
#                             fout = fout_base
#                         # todo future -> this could be replaced by a self.write_to_frontend_file(filename_out, )
#                         #
#                         # chose wich variable will be selected and copied into the outputfile 
#                         # from all the sources: observations_table, era5fb, header_table,...
#                         #
#                         # Dimensions (obs, trajectory)
#                         #
#                         for d, v in dims.items():
#                             fout.create_dataset(d, data=v)
#                             fout[d].attrs['NAME'] = np.string_('This is a netCDF dimension but not a netCDF variable.')
#                             fout[d].make_scale(d)  # resolves phony_dim problem

#                         fout.create_dataset('trajectory_index', data=trajectory_index)
#                         fout['trajectory_index'].attrs['long_name'] = np.string_(
#                             "index of trajectory this obs belongs to")
#                         fout['trajectory_index'].attrs['instance_dimension'] = np.string_("trajectory")
#                         fout['trajectory_index'].attrs['coordinates'] = np.string_("lat lon time plev")

#                         print('attrs written')

#                         #
#                         # Adding CDM
#                         #
#                         cdm_obstab = []
#                         cdm_eratab = []
#                         cdmlist = request.get('cdm', None)
#                         if cdmlist != None:
#                             # removing variables with restricted access
#                             if ('era5fb/obsvalue@body' in cdmlist) or ('observations_table/observation_value' in cdmlist):
#                                 try: cdmlist.remove('era5fb/obsvalue@body')
#                                 except: pass
#                                 try: cdmlist.remove('observations_table/observation_value')
#                                 except: pass
#                             for i in cdmlist:
#                                 grp = i.split('/')[0]
#                                 if i == 'era5fb':
#                                     logger.debug('Full er5fb CDM-request disabled.')
#                                 elif grp == 'era5fb':
#                                     try: 
#                                         var = i.split('/')[1]
#                                         cdm_eratab.append(var)
#                                         cfcopy[var]={'cdmname': i, 'units': '', 'shortname': var, 'coordinates': 'lat lon time plev', 'standard_name': var, 'cdsname': var,}
#                                     except: pass
#                                     cdmlist.remove(i)

#                                 if i == 'observations_table':
#                                     logger.debug('Full observations_table CDM-request disabled.')
#                                 elif grp == 'observations_table':
#                                     try: 
#                                         var = i.split('/')[1]
#                                         cdm_obstab.append(var)
#                                         cfcopy[var]={'cdmname': i, 'units': '', 'shortname': var, 'coordinates': 'lat lon time plev', 'standard_name': var, 'cdsname': var,}
#                                     except: pass
#                                     cdmlist.remove(i)
#                                     # NOW USE THOSE for the do_cfcopy below
#                             logger.debug('CDM - adding groups and variables: %s', str(cdmlist))
#                             for cdmstring in cdmlist:
#                                 print(cdmstring)
#                                 cdmsplit = cdmstring.split('/')
#             #                     try:
#                                 print(self.groups)
#                                 if cdmsplit[0] in self.groups:
#                                     print('group available')
#                                     # TODO: Add attributes and descriptions?
#                                     # whole group
#                                     if len(cdmsplit) == 1:
#                                         print('group copy')
#                                         fout.create_group(cdmsplit[0])
#                                         for cdmvar in self[cdmsplit[0]].keys():
#                                             try:
#                                                 fout[cdmsplit[0]].create_dataset(cdmvar, data=self[cdmsplit[0]][cdmvar][:])
#                                             except:
#                                                 pass
#                                     # single var of group 
#                                     if len(cdmsplit) == 2:
#                                         print('single var copy')
#                                         try:
#                                             fout.create_group(cdmsplit[0])
#                                         except:
#                                             pass # group alread exists?
#                                         try:
#                                             fout[cdmsplit[0]].create_dataset(cdmsplit[1], data=self[cdmsplit[0]][cdmsplit[1]][:])
#                                         except:
#                                             pass

#                         #
#                         # Variables based on cfcopy
#                         #
#                         if 'observations_table' in self.groups:
#                             igroup = 'observations_table'
#                             do_cfcopy(fout, self.file, igroup, idx, cfcopy, 'obs', compression,
#                                       var_selection=['observation_id', 'latitude', 'longitude', 'z_coordinate',
#                                                      'observation_value', 'date_time', 'sensor_id', 'secondary_value',
#                                                      'original_precision', 'reference_sensor_id', 'report_id','data_policy_licence']+cdm_obstab)
#                             # 'observed_variable','units'
#                             logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)
#                             if 'nodims' in request.keys():
#                                 fout['time'].make_scale('time')
#                         #
#                         # Feedback Information
#                         #
#                         if 'era5fb' in self.groups:
#                             igroup = 'era5fb'
#                             try:
#                                 do_cfcopy(fout, self.file, igroup, idx, cfcopy, 'obs', compression,
#                                           var_selection=['fg_depar@body', 'an_depar@body',
#                                                          'biascorr@body']+cdm_eratab)
#                                 # ['vertco_reference_1@body','obsvalue@body','fg_depar@body'])
#                                 logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)
#                             except KeyError as e:
#                                 raise KeyError('{} not found in {} {}'.format(str(e), str(request['optional']), self.name))
#                         #
#                         # Adjusted Values with ERA5
#                         # TODO: Names not set, group missing
#                         if 'adjera5' in self.groups:
#                             igroup = 'adjera5'
#                             try:
#                                 do_cfcopy(fout, self.file, igroup, idx, cfcopy, 'obs', compression,
#                                           var_selection=['bias_estimate', 'bias_estimation_method'])
#                                 logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)
#                             except KeyError as e:
#                                 raise KeyError('{} not found in {} {}'.format(str(e), str(request['optional']), self.name))

#                         if 'advanced_homogenisation' in self.groups :
#                             igroup = 'advanced_homogenisation'
#                             varsel=[]
#                             try:
#                                 for o in request['optional']:
#                                     if o in varseldict[request['variable']]:
#                                         varsel.append(o)
#                             except:
#                                 pass

#                             if varsel:       
#                                 try:
#                                     do_cfcopy(fout, self.file, igroup, idx, cfcopy, 'obs', compression,
#                                               var_selection=varsel)
#                                     logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)
#                                 except KeyError as e:
#                                     print(e)
#                                     #raise KeyError('{} not found in {} {}'.format(str(e), str(request['optional']), self.name))


#                         #
#                         # advanced_uncertainty
#                         # 
#                         if 'advanced_uncertainty' in self.groups:
#                             print('advanced_uncertainty in self.groups')
#                             igroup = 'advanced_uncertainty'
#                             try:
#                                 do_cfcopy(fout, self.file, igroup, idx, cfcopy, 'obs', compression,
#                                           var_selection=['desroziers_30', 'desroziers_60', 'desroziers_90', 'desroziers_180'])
#                                 logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)
#                             except KeyError as e:
#                                 raise KeyError('{} not found in {} {}'.format(str(e), str(request['optional']), self.name))

#                         #
#                         # Header Information
#                         #
#                         if 'header_table' in self.groups:
#                             igroup = 'header_table'
#                             # only records fitting criteria (zidx) are copied
#                             # todo why is lon, lat not here?
#                             do_cfcopy(fout, self.file, igroup, zidx, cfcopy, 'trajectory', compression,
#                                       var_selection=['report_id'])
#                             logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)
#                             # ,'station_name','primary_station_id'])
#                             # todo could be read from the observations_table
#                         #
#                         # Station Configuration
#                         #
#                         # station_configuration
#                         if 'station_configuration' in self.groups:
#                             igroup = 'station_configuration'
#                             cfcstationcon = {'station_name': 
#                                              {
#                                                  'cdmname': 'station_configuration/station_name',
#                                                  'units': 'NA',
#                                                  'shortname': 'station_id',
#                                                  'coordinates': 'lat lon time plev',
#                                                  'standard_name': 'station_name'
#                                              }
#                                             } 
#                             do_cfcopy(fout, self.file, igroup, idx, cfcstationcon, 'obs', compression,
#                                       var_selection=['station_name'])
#                             logger.debug('Group %s copied [%5.2f s]', igroup, time.time() - time0)

#                         #
#                         # Fix Attributes and Globals
#                         #
#                         fout['trajectory_label'].attrs['cf_role'] = np.string_('trajectory_id')
#                         fout['trajectory_label'].attrs['long_name'] = np.string_('Label of trajectory')
#                         for a, v in globatts.items():
#                             fout.attrs[a] = np.string_(v)

#                         fout.attrs['history'] = np.string_(
#                             'Created by Copernicus Early Upper Air Service Version 0, ' + datetime.now().strftime(
#                                 "%d-%b-%Y %H:%M:%S"))
#                         fout.attrs['license'] = np.string_('https://apps.ecmwf.int/datasets/licences/copernicus/')


#                         for i in fout.keys():    
#                             if 'toolbox' in request.keys() and not 'optional' in request.keys():
#             #                     if i in ['wind_from_direction']:
#             #                         fout['ta'] = fout[i]
#             #                         fout.__delitem__(i)
#                                 print('no optionals')
#                             elif 'toolbox' in request.keys():
#                                 if i in ['ta', 'hur', 'ua', 'va', 'wind_from_direction']:
#                                     oldunits = fout[i].attrs['units']
#                                     fout.__delitem__(i)
#                                     oldkey=(request['optional'][0])
#                                     fout[i]=fout[oldkey]
#                                     fout[i].attrs['units'] = oldunits
#                                     fout.__delitem__(oldkey)


#                         for i in fout.keys():
#                             #print(i)
#                             if 'nodims' in request.keys():
#                                 if (i == 'obs' or i == 'trajectory' or 'string' in i):
#                                     fout.__delitem__(i)
#                             version = ''
#                             if '.' in i:
#                                 vers = i.find('.')
#                                 version = i[(vers-2):(vers+2)] # '_x.x'
#                                 newname = i.replace(version, '')
#                                 fout[newname] = fout[i]
#                                 fout.__delitem__(i)
#                                 fout[newname].attrs['version'] = np.string_(version[1:]) # 'x.x'

#                         # removing variables with restricted access
#                         try: fout['observations_table'].__delitem__('observation_value')
#                         except: pass
#                         try: fout['era5fb'].__delitem__('obsvalue@body')
#                         except: pass
#                 except:
#                     print('---')
#                     time.sleep(np.absolute(np.random.randn()/5.))
#                     continue
#                 break
                    
        logger.debug('Finished %s [%5.2f s]', self.name, time.time() - time0)
        tt=time.time() - time0
        print(tt)
    
        # FIN
                   
    def read_write_request_parallel(self, filename_out: str, request: dict, cf_dict: dict):
        """ This is the basic request used in the cds_eua2 script

        Args:
            filename_out: request output filename, e.g. : /data/public//tmp//049532808458/dest_0-20000-0-10393_air_temperature.nc
            request: Request dictionary, minimum: {'variable' : 'temperature'}
            cf_dict: CF Convention for Names (results from read_standardnames() Function)

        Examples:
            >>> data = CDMDataset(filename)
            >>> data.read_write_request('test_cdm_frontend.nc', {'statid': '01001', 'variable': ['temperature']})
            Write all temperature data to a CDM frontend file
        """
        # version of process_flat with Leo's optimizations
        if 'variable' not in request.keys():
            logger.error('No variable specified %s %s', str(request.keys()), self.name)
            raise ValueError('No variable specified')

        variable = request['variable']
        compression = request['compression']
        if isinstance(variable, list):
            if len(variable) > 1:
                # todo potentially launch here recursive ?
                raise RuntimeError('Requests need to be split up by variable')
            variable = variable[0]  # assumes a list ?

        if cf_dict is None:
            cf_dict = read_standardnames()

        cdsname = cdm_to_cds(variable, cf=cf_dict)
        cdmattrs = cf_dict[cdsname]  #
        cdmnum = cdmattrs['cdmcode']  # 85 for air_temperature
        time0 = time.time()

        trange, mask = self.read_observed_variable(cdmnum,
                                                   variable='observation_value',
                                                   dates=request.get('date', None),
                                                   plevs=request.get('pressure_level', None),
                                                   times=request.get('time', None),
                                                   observed_variable_name='observed_variable',
                                                   date_time_name='date_time',
                                                   date_time_in_seconds=False,
                                                   z_coordinate_name='z_coordinate',
                                                   group='observations_table',
                                                   dimgroup='observations_table',
                                                   return_index=True,
                                                   rtsindex=request.get('rtsidx',None)
                                                   )
        logger.debug('Datetime selection: %d - %d [%5.2f s] %s', trange.start,
                     trange.stop, time.time() - time0, self.name)
        tt=time.time() - time0
        print(tt)
        idx = np.where(mask)[0] + trange.start  # absolute integer index
        if len(idx) == 0:
            logger.warning('No matching data found %s', self.name)
            raise ValueError('No matching data found')  # add CDMname for logging
            #return

        logger.debug('Data found: %d %s', len(idx), self.name)
        #
        # Make Trajectory Information (lon, lat, profile id, ...)
        #
        trajectory_index = np.zeros_like(idx, dtype=np.int32)
        trajectory_index2 = np.zeros_like(idx, dtype=np.int32)
        if 'recordinindex' in self.groups:
            # unsorted indices in root
            recordindex = self.file['recordindex'][()]
        else:
            # sorted indices are in recordindices group / by variable
            recordindex = self.file['recordindices'][str(cdmnum)][()]  # values
            
        zidx = np.where(np.logical_and(recordindex >= trange.start, recordindex <= trange.stop))[0]
        zidx = calc_trajindexfastl(recordindex, zidx, idx, trajectory_index)
        #
        # Dimensions and Global Attributes
        #
        tt=time.time() - time0
        print(tt)
        dims = {'obs': np.zeros(idx.shape[0], dtype=np.int32),
                'trajectory': np.zeros(zidx.shape[0], dtype=np.int32)}
        globatts = get_global_attributes()  # could put more infors there ?
        #
        # Definition of Variables to write
        #
        # Common Variables needed for a requested file
        #
        snames = ['platform_id', 'platform_name', 'observation_value', 'latitude',
                  'longitude', 'time', 'air_pressure', 'trajectory_label', 
                  'report_id', 'station_id']
        varseldict={}
        varseldict['temperature']=['RAOBCORE_bias_estimate', 'RASE_bias_estimate', 'RICH_bias_estimate', 'RISE_bias_estimate']
        varseldict['wind_direction']=['wind_bias_estimate']
        varseldict['u_component_of_wind']=['wind_bias_estimate']
        varseldict['v_component_of_wind']=['wind_bias_estimate']
        varseldict['relative_humidity']=['humidity_bias_estimate']
        # added report_id -> in observations_table, not to be confused with report_id from header_table -> trajectory_label
        logger.debug('Request-keys: %s', str(request.keys()))
        snames.append(cdsname)  # Add requested variable
        #
        # Add Optional Variables
        # Check read_standard_names for correct names otherwise it will be rejected
        #
        if 'optional' in request.keys():
            snames.extend(request['optional'])
            
        if 'bias_estimate_method' in request.keys():
            # use era5adj group
            pass
        
        # Copy Metadata -> used by do_cfcopy
        cfcopy = {}  # Copy of CDM Info Dictionary (read_standard_names())
        for ss in snames:
            try:
                cfcopy[ss] = cf_dict[ss]
            except:
                pass
        #
        # End Definition of Variables to write
        #
        logger.debug('Writing: %s', filename_out)
        tt=time.time() - time0
        print(tt)
        
        ###
        # start splitting and parallel work here:
        ###
        
        # outputfiles:
        del cfcopy['platform_id']
        del cfcopy['platform_name']

        list_cfcopy = []
        for chunks in cfcopy:
            list_cfcopy.append({chunks:cfcopy[chunks]})
        list_cfcopy.append({'0':0})
        list_cfcopy.append({'1':1})
        list_cfcopy.append({'2':2})
        list_cfcopy.append({'3':3})
#         for o in range(len(list_cfcopy)):
#             list_cfcopy[o]=[list_cfcopy[o], o]
        
        ###pool###
        fname = str(self.filename)
        func = partial(parallel_writing, fname, dims, trajectory_index, idx, zidx, compression, filename_out, request, globatts)
        with multiprocessing.get_context('spawn').Pool(10) as P:
#             P = Pool(10)
            pool_out = list(P.map(func, list_cfcopy, chunksize=3))
            print(pool_out)
        filelist = glob.glob(filename_out[:-3]+'_*')
        
        ###
        # Concatination of single files
        ###
        
        with h5py.File(filename_out, 'w') as merge:
            print(filename_out)
            for r in filelist:
                print(r)
                with h5py.File(r, 'r') as filetocopy:
                    for i in filetocopy.keys():
                        if not i in merge.keys():
                            print(i)
                            filetocopy.copy(i, merge, name=i)
        logger.debug('netcdfs merged [%d] to %s', len(filelist), filename_out)
        logger.debug('Finished %s [%5.2f s]', self.name, time.time() - time0)
        tt=time.time() - time0
        print(tt)
    
        # FIN

    def make_datetime_slice(self, varnum:int = None, dates: list = None, date_time_name: str = 'date_time',
                            date_time_in_seconds: bool = False,
                            add_day_before: bool = False,
                            **kwargs) -> slice:
        """ Create a datetime slice for HDF5 optimal reading

        Args:
            varum: variable id number, e.g. 85 for temperature
            dates: date selection
            date_time_name: name of the datetime variable
            date_time_in_seconds: are dates in seconds ?
            add_day_before: do we have time in the request and need to make sure we have early launches?

        Returns:
            slice : datetime slice
        """
        if self.hasgroups:
            # CDM Backend file
            group = 'observations_table'
            date_time_name = date_time_name if date_time_name != 'date_time' else 'date_time'
        else:
            # CDM Frontend file
            group = None
            date_time_name = date_time_name if date_time_name != 'date_time' else 'time'

        if dates is not None:
            if 'recordtimestamp' in self.groups:
                timestamp = self.load_variable_from_file('recordtimestamp', return_data=True)[0]
                timestamp_units = self.read_attributes('recordtimestamp').get('units', None)
            elif 'recordindices' in self.groups:
                if isinstance(self['recordindices'], CDMGroup):
                    timestamp = self.load_variable_from_file('recordtimestamp', group='recordindices', return_data=True)[0]
                    timestamp_units = self.read_attributes('recordtimestamp', group='recordindices').get('units', None)
            else:
                # backup if recordtimestamp not present
                timestamp = self.load_variable_from_file(date_time_name, return_data=True)[0]
                timestamp_units = self.read_attributes(date_time_name).get('units', None)

            # Make sure units are accordingly
            if timestamp_units is None:
                timestamp_units = 'seconds since 1900-01-01 00:00:00'
            if 'seconds' not in timestamp_units:
                timestamp = to_seconds_since(timestamp, timestamp_units)

            if not date_time_in_seconds:
                # to seconds since 1900-01-01
                if len(dates) == 1:
                    dates = [convert_datetime(dates[0], return_seconds=True)]
                else:
                    dates = [convert_datetime(dates[0], return_seconds=True),
                             convert_datetime(dates[-1], return_seconds=True)]

            if add_day_before:
                logic = (timestamp >= (dates[0] - 86400)) & (timestamp <= (dates[-1] + 86399))
            else:
                logic = (timestamp >= (dates[0])) & (timestamp <= (dates[-1] + 86399))

            timeindex = np.where(logic)[0]

            if timeindex.shape[0] == 0:
                logger.warning('No data in time interval %s', self.name)
                raise ValueError('No data in specified time interval')

                
            if 'recordindex' in self.groups:
                # read correct var from group recordindices
                recordindex = self.load_variable_from_file('recordindex', return_data=True)[0]
                if timeindex[-1] < (recordindex.shape[0] - 1):
                    # within datetime range
                    trange = slice(recordindex[timeindex[0]], recordindex[timeindex[-1] + 1])
                else:
                    #
                    trange = slice(recordindex[timeindex[0]], self[group][date_time_name].shape[0])
                    
            elif 'recordindices' in self.groups:
                # update for sorted backend files
                try:
                    
                    recordindex = self.load_variable_from_file(str(varnum), group='recordindices', return_data=True)[0]
                except:
                    logger.warning(str(varnum)+' not found in recordindices of file ')
                    trange = slice(0,0)
                    return trange
                # last time index 
                if timeindex[-1] < (recordindex.shape[0] - 1):
                    # within datetime range
                    trange = slice(recordindex[timeindex[0]], recordindex[timeindex[-1] + 1])
                else:
                    #
                    trange = slice(recordindex[timeindex[0]], recordindex[-1])
#                 itx = np.isfinite(recordindex)
                
#                 if not np.all(itx):
#                     # no values:
#                     # return slice(0, None)
#                     raise ValueError("No data for variable")
                
#                 if timeindex[-1] < (recordindex.shape[0] - 1):
#                     # within datetime range
#                     trange = slice(int(recordindex[itx][timeindex[0]]), int(recordindex[itx][timeindex[-1] + 1]))
#                 else:
#                     #
#                     trange = slice(int(recordindex[itx][timeindex[0]]), int(recordindex[itx][-1]))
            else:
                trange = slice(timeindex[0], timeindex[-1] + 1)

            time_units = self.read_attributes(date_time_name, group=group).get('units', '')
            if timestamp_units != time_units:
                logger.warning('Timeunits missmatch? %s <> %s', timestamp_units, time_units)
                raise ValueError('Timeunits missmatch?', timestamp_units, time_units)
        else:
            # ALL
            if 'recordindices' in self.groups:
                recordindex = self.load_variable_from_file(str(varnum), group='recordindices', return_data=True)[0]
#                 itx = np.isfinite(recordindex)
#                 trange = slice(int(recordindex[itx][0]), int(recordindex[itx][-1]))  # all
                trange = slice(recordindex[0], recordindex[-1])
            elif group is not None:
                trange = slice(0, self[group][date_time_name].shape[0])
            
            else:
                trange = slice(0, self[date_time_name].shape[0])

        logger.debug('Datetime selection: %d - %d', trange.start, trange.stop)
        return trange

    def read_observed_variable(self, varnum: int,
                               variable: str = 'observation_value',
                               dates: list = None,
                               plevs: list = None,
                               times: list = None,
                               observed_variable_name: str = 'observed_variable',
                               date_time_name: str = 'date_time',
                               date_time_in_seconds: bool = False,
                               z_coordinate_name: str = 'z_coordinate',
                               group: str = 'observations_table',
                               dimgroup: str = 'observations_table',
                               return_coordinates: bool = False,
                               return_index: bool = False,
                               use_odb_codes: bool = False,
                               return_xarray: bool = False,
                               **kwargs
                               ):
        """ Read a variable from a CDM backend file
        Uses recordtimestamp and observed_variable for subsetting and z_coordinate as well

        Args:
            varnum: 85 or 2 for temperature
            variable: CDM variable name: observation_value
            dates: [start end]
            plevs: [plevs] in Pa
            times: [sounding times] in hours
            observed_variable_name:
            date_time_name:
            date_time_in_seconds:
            z_coordinate_name:
            group: group of variable
            dimgroup: group of date_time and z_coordinate and observed_variable
            return_coordinates: dates and pressure levels
            return_index: subset and logic array
            use_odb_codes: ODB Codes or CDM Codes
            return_xarray: convert to xarray object

        Returns:
            values
            values, dates, pressure
            trange, index
            trange, index, dates, pressure
            DataArray

        Examples:
            Read all Temperatures at 500 hPa from 2000 to 2019
            >>> tmp.read_observed_variable(85, dates=['2000-01-01','2019-12-31'], plevs=[50000])
            2020-07-15 12:58:02,518 - cdm | read_observed_variable - INFO - [READ] recordtimestamp: slice(6350527, 137146806, None)
            2020-07-15 12:59:13,630 - cdm | read_observed_variable - INFO - [READ] Observed variable 85
            2020-07-15 12:59:14,748 - cdm | read_observed_variable - INFO - [READ] pressure levels [50000]
            array([246.1 , 238.1 , 239.7 , ..., 239.92, 238.77, 237.7 ], dtype=float32)

        Notes:
            The CDM Backend files are based on ECMWF ODB format and are ragged arrays or tables. Variables are
            mixed, as they are available, like records. Hence to retrieve one variable the observed_variable needs
            to be used to subset the array in memory (faster)
        """
        
        if not self.da:        
            if not self.hasgroups:
                raise RuntimeError('This function only works with CDM Backend files')
    
            if dimgroup not in self.groups:
                raise ValueError('Missing group?', dimgroup)
    
            if group not in self.groups:
                raise ValueError('Missing group?', group)
    
            if not isinstance(variable, str):
                raise ValueError('(variable) Requires a string name, not ', str(variable))
    
            if not isinstance(varnum, int):
                raise ValueError('(varnum) Requires a integer number, not', str(varnum))
    
            if use_odb_codes:
                if varnum not in odb_codes.values():
                    raise ValueError('(varnum) Code not in ODB Codes', variable, str(odb_codes))
            else:
                if varnum not in cdm_codes.values():
                    raise ValueError('(varnum) Code not in CDM Codes', variable, str(cdm_codes))
    
            # TODO for the new files this variable is not mandatory anymore !!!
            if observed_variable_name not in self[dimgroup].keys():
                raise ValueError('Observed variable not found:', observed_variable_name, self[dimgroup].keys())
    
        if dates is not None:
            if not isinstance(dates, list):
                dates = [dates]

        if plevs is not None:
            if not isinstance(plevs, (list, np.ndarray)):
                plevs = [plevs]
            plevs=np.asarray(plevs,dtype=np.float32)
            
        
        if self.da and 'rtsindex' in kwargs.keys():
            svarnum=str(varnum)
            ssvar=self.file['recordindices'][svarnum]
            minday=0
            if times is not None:
                if int(times[0])>int(times[-1]):
                    if kwargs['rtsindex'][0]>0:
                        minday=-1
            if kwargs['rtsindex'][1]<ssvar.shape[0]:
                trange=slice(ssvar[kwargs['rtsindex'][0]+minday],ssvar[kwargs['rtsindex'][1]],None)
            else:
                trange=slice(ssvar[kwargs['rtsindex'][0]+minday],ssvar[-1],None)
                
        else:
            
            trange = self.make_datetime_slice(dates=dates, varnum=varnum, date_time_name=date_time_name,
                                          date_time_in_seconds=date_time_in_seconds,
                                          add_day_before=True if times is not None else False)
        
        if dates is not None:
            print('dates not None')
            if not self.da:
                print('da False')
                xdates = self[dimgroup][date_time_name][trange]
            else:
                if times is not None:
                    xdates = self.file[dimgroup][date_time_name][trange]
                
        else:
            print('dates None')
            xdates = self.load_variable_from_file(date_time_name, group=dimgroup, return_data=True)[0][trange]
        
        #
        # Observed Code
        #
        if dates is None and 'recordindex' in self.groups:
            self.load_variable_from_file(observed_variable_name, group=dimgroup)

        if False:
            # not really faster ?
            logic = np.ones(trange.stop - trange.start, dtype=np.bool)
            andisin(logic, self[dimgroup][observed_variable_name][trange], np.asarray([varnum], dtype=np.int32))
        elif 'recordindex' in self.groups:
            logic = (self[dimgroup][observed_variable_name][trange] == varnum)
        else:
            logic = np.ones(trange.stop-trange.start, dtype=np.bool)
            # todo apply trange before, to make a subset 
        logger.info('[READ] Observed variable %s', varnum)
        #
        # Pressure levels
        #
        xplevs = None
        if plevs is not None:
            if self.da:
                xplevs=self.file[dimgroup][z_coordinate_name][trange]
            else:
                if z_coordinate_name not in self[dimgroup].keys():
                    raise ValueError('Pressure variable not found:', z_coordinate_name, self[dimgroup].keys())
                p_attrs = self.read_attributes(z_coordinate_name, group=dimgroup)
                p_units = p_attrs.get('units', 'Pa')
                if p_units != 'Pa':
                    RuntimeWarning('Pressure variable wrong unit [Pa], but', p_units)
    
                if dates is None:
                    xplevs = self.load_variable_from_file(z_coordinate_name, group=dimgroup, return_data=True)[0][trange]
                else:
                    xplevs = self[dimgroup][z_coordinate_name][trange]
            if len(plevs) == 1:
                logic = logic & (xplevs == plevs[0])
            else:
                andisin(logic, xplevs, np.asarray(plevs))  # for cube that make 2s
                # logic = logic & (np.in1d(xplevs, plevs))
            logger.info('[READ] pressure levels %s', str(plevs))
        #
        # Times
        #
        if times is not None:
            # standard times [0-23], day before has been added
            times = totimes(times)  # not sure why this does notsort ?
            #times = np.sort(totimes(times))  # not sure why this does notsort ?
            hours = np.empty_like(xdates)
            date_shift = np.zeros_like(xdates, dtype=np.int32)  # earlier date, but would be
            days = np.empty_like(xdates)
            tohourday(hours, days, xdates, date_shift)
            # all previous day profiles should be included due to day before flag
            # match all times in hours (logical and)
            andisin(logic, hours, times)
            if times[0]>times[-1]:
                # use only late hours for the first day
                first_day = days[0]  # first day days[logic][0]
                first_day_logic = logic[days == first_day]
                last_day = days[-1]  # first day days[logic][0]
                last_day_logic = logic[days == last_day]
                # only true if lat times
                andisin(first_day_logic, hours[days == first_day], times[times >= times[0]])
                andisin(last_day_logic, hours[days == last_day], times[times < times[0]])
                logic[days == first_day] = first_day_logic  # modify logic
                logic[days == last_day] = last_day_logic  # modify logic
                logger.info('[READ] first day (%s) times %d/%d', seconds_to_datetime([first_day * 86400]),
                            first_day_logic.sum(), first_day_logic.size)

        if return_xarray:
            return_coordinates = True

        if return_coordinates:
            if xdates is None:
                if dates is None:
                    xdates = self.load_variable_from_file(date_time_name, group=dimgroup, return_data=True)[0][trange]
                else:
                    xdates = self[dimgroup][date_time_name][trange]

            if xplevs is None:
                if dates is None:
                    xplevs = self.load_variable_from_file(z_coordinate_name, group=dimgroup, return_data=True)[0][trange]
                else:
                    xplevs = self[dimgroup][z_coordinate_name][trange]

        if return_index:
            if return_coordinates:
                return trange, logic, xdates[logic], xplevs[logic]  # no trange here????
            return trange, logic

        if dates is None:
            self.load_variable_from_file(variable, group=group)

        if return_xarray:
            # todo add attributes from CF Table
            data = xr.DataArray(self[group][variable][trange][logic], dims=('obs'))
            logger.info('[READ] xarray ... %s', data.shape)
            if xdates is not None:
                data[date_time_name] = ('obs', seconds_to_datetime(xdates[logic]))
                logger.debug('[READ] datetime conversion %d', data[date_time_name].size)
            if xplevs is not None:
                data[z_coordinate_name] = ('obs', xplevs[logic])
                data[z_coordinate_name].attrs.update({'units': 'Pa', 'standard_name': 'air_pressure'})
            if len(data.coords) > 0:
                data = data.set_index(obs=list(data.coords))
            return data

        if return_coordinates:
            return self[group][variable][trange][logic], xdates[logic], xplevs[logic]

        return self[group][variable][trange][logic]

    def read_variable(self, name,
                      dates: list = None,
                      plevs: list = None,
                      times: list = None,
                      date_time_name: str = 'time',
                      date_time_in_seconds: bool = False,
                      z_coordinate_name: str = 'plev',
                      return_coordinates: bool = False,
                      return_index: bool = False,
                      return_xarray: bool = False):
        """ Read a variable from a CDM frontend file

        Args:
            name: name of the variable
            dates: [start end] datetime selection
            plevs: [pressure levels]
            date_time_name: Name of the datetime variable
            date_time_in_seconds: dates are in seconds since 1900-01-01 00:00:00
            z_coordinate_name: Name of the pressure level variable
            return_coordinates: add coordinate information
            return_index: no data only indices
            return_xarray: convert to xarray object including coordinates and attributes

        Returns:
            values
            values, dates, pressure
            trange, index
            trange, index, dates, pressure
            DataArray
        """

        if self.hasgroups:
            raise RuntimeError('Only for CDS frontend files')

        if not isinstance(name, str):
            raise ValueError('(variable) Requires a string name, not ', str(name))

        if name not in self.groups:
            raise ValueError('Variable not found', name)

        if dates is not None:
            if date_time_name not in self.groups:
                raise ValueError('date_time_name Variable not found', date_time_name)

            if not isinstance(dates, list):
                dates = [dates]

        if plevs is not None:
            if z_coordinate_name not in self.groups:
                raise ValueError('z_coordinate_name Variable not found', z_coordinate_name)

            if not isinstance(plevs, (list, np.ndarray)):
                plevs = [plevs]

        trange = self.make_datetime_slice(dates=dates, date_time_name=date_time_name,
                                          date_time_in_seconds=date_time_in_seconds)
        if dates is not None:
            xdates = self[date_time_name][trange]
        else:
            xdates = self.load_variable_from_file(date_time_name, return_data=True)[0]
        # trange = slice(None)
        # xdates = None
        # logic = np.ones(self[name].shape, dtype=np.bool)  # all True
        # if dates is not None:
        #     timestamp = self[date_time_name][()]
        #     d_attrs = self.read_attributes(date_time_name)
        #     time_units = d_attrs.get('units', 'seconds since 1900-01-01 00:00:00')
        #     if 'seconds' not in time_units:
        #         dates = to_seconds_since(timestamp, time_units)
        #     if not date_time_in_seconds:
        #         # to seconds since
        #         if len(dates) == 1:
        #             dates = [convert_datetime(dates[0], return_seconds=True)]
        #         else:
        #             dates = [convert_datetime(dates[0], return_seconds=True),
        #                      convert_datetime(dates[-1], return_seconds=True)]
        #     logic = (timestamp >= dates[0]) & (timestamp <= dates[-1])
        #     timeindex = np.where(logic)[0]
        #     trange = slice(timeindex[0], timeindex[-1])
        #     xdates = timestamp[trange]
        #     logger.info('[READ] %s : %s [%s]', date_time_name, str(trange), time_units)

        xplevs = None
        if plevs is not None:
            xplevs = self[z_coordinate_name][trange]
            p_attrs = self.read_attributes(z_coordinate_name)
            p_units = p_attrs.get('units', 'Pa')

            if len(plevs) == 1:
                logic = xplevs == plevs[0]
            else:
                # could be fast with andisin
                logic = np.in1d(xplevs, plevs)
            logger.info('[READ] %s : %s [%s]', z_coordinate_name, str(plevs), p_units)
        else:
            logic = np.ones(xdates.size, dtype=np.bool)

        if times is not None:
            # todo add time selection
            pass

        if return_xarray:
            return_coordinates = True  # need to be available for coordinates

        if return_coordinates:
            if xdates is None:
                try:
                    xdates = self[date_time_name][trange]  # read the data
                except:
                    pass
            if xplevs is None:
                try:
                    xplevs = self[z_coordinate_name][trange]  # read the data
                except:
                    pass

        # return only indices for repeated use, e.g. temperature + fbstats
        if return_index:
            if return_coordinates:
                return trange, logic, xdates[logic], xplevs[logic]
            return trange, logic

        if dates is None:
            self.load_variable_from_file(name)  # performance might be better if the whole is read and subset

        if return_xarray:
            # create a DataArray
            data = xr.DataArray(self[name][trange][logic], dims=('obs',), attrs=self.read_attributes(name))
            if xdates is not None:
                data[date_time_name] = ('obs', seconds_to_datetime(xdates[logic]))
                data[date_time_name].attrs.update(self.read_attributes(date_time_name))
            if xplevs is not None:
                data[z_coordinate_name] = ('obs', xplevs[logic])
                data[z_coordinate_name].attrs.update(self.read_attributes(z_coordinate_name))
            if len(data.coords) > 0:
                data = data.set_index(obs=list(data.coords))
            return data

        if return_coordinates:
            return self[name][trange][logic], xdates[logic], xplevs[logic]
        return self[name][trange][logic]

    def to_dataframe(self, groups=None, variables: list = None, date=None,
                     date_time_name: str = 'date_time',
                     date_is_index: bool = False,
                     decode_datetime: bool = True,
                     **kwargs):
        """ Convert variables to a DataFrame
        CDM Backend files need groups and variables
        CDM Frontend files do not

        Args:
            groups: list of groups to search in
            variables: names of variables
            date: datetime string selection
            date_time_name: name of datetime variable
            date_is_index: are date an index?
            **kwargs:

        Returns:
            DataFrame
        """
        # TODO change to search for variables and select whole groups
        # it is inconsistent to use variables + groups to select !!!
        # TODO add feature to select only one variable
        if self.hasgroups:
            if isinstance(groups, str):
                if groups == '/':
                    groups = None
                else:
                    groups = [groups]

            if isinstance(variables, str):
                variables = [variables]
        else:
            if variables is None:
                variables = self.groups  # all variables
                
            if date_time_name == 'date_time':
                date_time_name = 'time'
                
        if date is not None and 'recordindex' in self.groups:
            # todo missing feature, ask for varnum
            if isinstance(date, (str, int)):
                date = [date]

            if date_is_index:
                date = slice(date[0], date[-1] + 1)
            else:
                date = self.make_datetime_slice(dates=[date[0], date[-1]], date_time_name=date_time_name)
        else:
            date = slice(None)

        logger.info("Reading Profile on %s", str(date))
        data = {}

        # update / need to check
        # 
        # for ivar,igroup in self.select_variables(variables, groups=groups):
        #   data[ivar] = self[igroup][ivar]

        if groups is not None:
            # Backend file
            for igroup in groups:
                for ivar in variables:
                    if ivar in self[igroup].keys():
                        data[ivar] = self[igroup][ivar][date]
        else:
            # Frontend file
            for ivar in variables:
                if ivar in self.groups:
                    if ivar in ['trajectory_label', 'trajectory']:
                        data[ivar] = self[ivar][:]  # get data
                        continue
                    if 'string' in ivar:
                        continue
                    
                    data[ivar] = self[ivar][date]  # get data
        # raise RuntimeError()
        logger.debug('Read variables: %s', str(data.keys()))
        for ivar in list(data.keys()):
            # check if var still exists
            if ivar not in data.keys():
                continue

            # convert char arrays to string
            if len(data[ivar].shape) > 1:
                data[ivar] = data[ivar].astype(object).sum(axis=1).astype(str)
        
            # Trajectory info
            if ivar == 'trajectory_label':
                data['trajectory_label'] = data['trajectory_label'][data['trajectory_index']]
                try:
                    del data['trajectory_index']
                except:
                    pass
                try:
                    del data['trajectory']
                except: 
                    pass
            if ivar == date_time_name and decode_datetime:
                data[ivar] = seconds_to_datetime(data[ivar])
                
        return pd.DataFrame(data)

    def read_data_to_cube(self, variables: list, dates: list = None, plevs: list = None, feedback = None,
                          feedback_group = 'era5fb', **kwargs) -> dict:
        """ Read standard pressure levels and return a DataCube

        Args:
            variables: list of variables, e.g. temperature
            dates: [start, stop], str, int or datetime
            plevs: [list] in Pa or hPa
            feedback: list of feedback variables
            feedback_group: list of feedback groups (era5fb, adjust, ...)
            **kwargs:

        Optional Keywords:
            date_time_name: Name of the datetime variable
            z_coordinate_name: Name of the pressure level variable

        Returns:
            dict : {variable : xr.DataArray}
        """
        if len(variables) == 0:
            raise ValueError('Need a variables', str(cdm_codes.keys()))

        if isinstance(variables, str):
            variables = [variables]

        data = {}
        if plevs is not None:
            plevs = np.asarray(plevs)
            if any((plevs > 110000) | (plevs < 500)):
                raise ValueError('Pressure levels outside range [5, 1100] hPa')
        else:
            plevs = std_plevs * 100  # in Pa

        std_plevs_indices = np.zeros(1001, dtype=np.int32)  # hPa
        # in hPa
        for i, j in enumerate(plevs // 100):
            std_plevs_indices[j] = i

        # todo check if variable can be replaced inside HDF5 ?
        if self.hasgroups:
            #
            # check variables
            #
            varnum = []
            for ivar in variables:
                if ivar not in cdm_codes.keys():
                    raise ValueError('Variable not found', ivar)

                varnum.append({'varnum': cdm_codes[ivar], 'out_name': ivar,'variable': 'observation_value', 'group': 'observations_table', 'bkp_var': ivar})
                if feedback is not None:
                    if isinstance(feedback, str):
                        feedback = [feedback]
                    if isinstance(feedback_group, str):
                        feedback_group = [feedback_group]

                    for jvar, jgroup in self.select_variables(feedback, groups=feedback_group):
                        # jvar -> @body rename
                        # MB Update: new output name based on ivar and jvar
                        varnum.append({'varnum': cdm_codes[ivar],
                                       'out_name': '{}_{}'.format(ivar, jvar.split('@')[0] if '@' in jvar else jvar),
                                       'variable': jvar,
                                       'group': jgroup,
                                       'bkp_var': ivar})
            # multiprocessing of the requests?
            for ivarnum in varnum:
                logger.info('Reading ... %d  %s', ivarnum['varnum'], ivarnum['variable'])
                ivarnum.update(kwargs)
                # TIME SLICE, INDEX, SECONDS ARRAY, PRESSURE LEVELS
                trange, indices, secarray, pressure = self.read_observed_variable(dates=dates,
                                                                                  plevs=plevs,
                                                                                  return_coordinates=True,
                                                                                  return_index=True,
                                                                                  **ivarnum)
                logger.info('[CUBE] Variable Group %d %s %s', ivarnum['varnum'], str(trange), str(secarray.shape))
                obs = self[ivarnum['group']][ivarnum['variable']][trange][indices]
                #
                # to Cube
                #
                # requires hPa for indices
                itime, iobs = table_to_cube(secarray,
                                            std_plevs_indices[pressure.astype(np.int32) // 100],
                                            obs,
                                            nplev=plevs.size)
                logger.info('[CUBE] %s %s', ivarnum['variable'], iobs.shape)
                v_attrs = get_attributes(cdmcode=ivarnum['varnum'],
                                         feedback=ivarnum['variable'] if ivarnum['variable'] != 'observation_value' else None)
                if len(v_attrs) > 0:
                    v_attrs = v_attrs[list(v_attrs.keys())[0]]
                # Convert to Xarray [time x plev]
                data[ivarnum['out_name']] = xr.DataArray(iobs,
                                                        coords=(seconds_to_datetime(secarray[itime]), plevs),
                                                        dims=('time', 'plev'),
                                                        name=ivarnum['bkp_var'],
                                                        attrs=v_attrs,
                                                        )
                # todo add attributes for coordinates
        else:
            for ivar in variables:
                # Read Attributes Variable
                v_attrs = self.read_attributes(ivar)
                # why no trange?
                iobs, secarray, pressure = self.read_variable(ivar,
                                                              dates=dates,
                                                              plevs=plevs,
                                                              return_coordinates=True)

                itime, iobs = table_to_cube(secarray,
                                            std_plevs_indices[pressure.astype(np.int32) // 100],
                                            iobs)
                logger.info('[CUBE] %s %s', ivar, iobs.shape)
                # Convert to Xarray [time x plev]
                data[ivar] = xr.DataArray(iobs,
                                          coords=(seconds_to_datetime(secarray[itime]), plevs),
                                          dims=('time', 'plev'),
                                          name=ivar,
                                          attrs=v_attrs)
                # data[ivar].data['time'].attrs.update(self.read_attributes(kwargs.get('date_time_name', 'time')))
                # data[ivar].data['plev'].attrs.update(self.read_attributes(kwargs.get('z_coordinate_name', 'plev')))
        return data

    def convert_to_raobcore(self, variable: str, filename: str = None, dates: list = None, plevs: list = None,
                            times=[0, 12], span=3, freq='12h', feedback: list = None, feedback_group: str = 'era5fb',
                            source: str = 'RAOBCORE/RICH v1.7.2 + solar elevation dependency (from 197901 onward)',
                            title: str = 'Station daily temperature series with ERA5/NOAA_20CR background departure statistics and RISE bias estimates',
                            attrs: dict = None,
                            global_attrs: dict = None,
                            force:bool = False,
                            **kwargs):
        """Convert CDM backend file to RAOBCORE/RICH data format

        Args:
            variable (str): request variable, e.g. temperature
            filename (str, optional): output filename. Defaults to None.
            dates (list, optional): datetime selection. Defaults to None.
            plevs (list, optional): pressure level selection. Defaults to None.
            times (list, optional): time selection. Defaults to [0, 12].
            span (int, optional): allowed deviation from time. Defaults to 3.
            freq (str, optional): should correspond to times. Defaults to '12h'.
            feedback (list, optional): feedback variables. Defaults to None.
            feedback_group (str, optional): feedback groups. Defaults to 'era5fb'.
            source (str, optional): name of product. Defaults to 'RAOBCORE/RICH v1.7.2 + solar elevation dependency (from 197901 onward)'.
            title (str, optional): title of product. Defaults to 'Station daily temperature series with JRA55/CERA20C/ERApreSAT background departure statistics and RISE bias estimates'.
            attrs (dict, optional): attributes to variables. Defaults to None. e.g. {'temperatures' : {'units' : 'K'}}
            global_attrs (dict, optional): global attributes. Defaults to None.

        Raises:
            RuntimeError: not a backend file

        Returns:
            xr.Dataset: data cube [hour x pressure x time]

        Examples:
            Convert temperature to RAOBCORE FORMAT with ERA5 first guess feedback (era5fb group) and 
            RAOBCORE bias estimates (advanced_homogenisation group):

            >>> data.convert_to_raobcore('temperature', 'ERA5fc_bias_estimate_temperature.nc', 
                                        feedback=['fg_depar@body', 'RAOBCORE_bias_estimate'], 
                                        feedback_group=['era5fb', 'advanced_homogenisation'])

        """
        if not self.hasgroups:
            raise RuntimeError("Requires a backend file")
        
        if feedback:
            if not isinstance(feedback, list):
                feedback = [feedback]
        else:
            feedback = ['biascorr@body', 'fg_depar@body', 'an_depar@body']

        default_attrs = {
            "lat":  {'axis': 'Y', 'units': 'degrees_north', 'long_name': 'station latitude'},
            "lon": {'axis': 'X', 'units': 'degrees_east', 'long_name': 'station longitude'},
            "alt": {'axis': 'Z', 'units': 'm', 'long_name': 'station altitude'},
            "press": {"long_name": "pressure levels", "units": "hPa", "axis": "Z"},
            "datum": {"long_name": "datum", "axis": "T"},
            "hours": {"long_name": "launch time", "units": "hr", "valid_range": (0, 23)}
        }
        if attrs:
            default_attrs.update(attrs)

        default_global = {'Conventions': 'CF-1.1',
                            'title': title,
                            'institution': 'University of Vienna',
                            'history': '2020/02/15',
                            'source': source,
                            'references': 'www.univie.ac.at/theoret-met/research/raobcore',
                            'Stationnname': self.header_table.station_name[-1, :].tobytes().decode()}
        if global_attrs:
            default_global.update(global_attrs)

        # They will all have the same time and plev shapes
        # per variable of course (temp, hum, wind)
        tdata = self.read_data_to_cube(
                        variable, 
                        dates=dates, 
                        plevs=plevs, 
                        feedback=feedback, 
                        feedback_group=feedback_group, 
                        **kwargs)
        dim = 'time'
        plev = 'plev'
        # need standard times -> 0,12
        tdata[variable] = align_datetime(tdata[variable], times=times, span=span, freq=freq, dim=dim, plev=plev)
        icoord = 'standard_%s' % dim
        # Index for only standard times
        standard_index = np.where(tdata[variable][icoord + '_flag'].values == 1)[0]
        tdata = xr.Dataset(tdata)  # Convert to Dataset
        # 1. Select only standard times
        # 2. Swap dimension to new standard time
        # 3. Convert to a cube hour x time x pressure
        tdata = stack_cube_by_hour(tdata.isel(**{dim: standard_index}).swap_dims({'time': icoord}), dim=icoord, times=times)
        tdata = tdata.drop('stacking')  # drop stacking dimension, used to revert stacking
        # add hours variable
        # variable with the correct hours inside
        tdata['hours'] = tdata.time.dt.hour
        # reset time to index
        tdata['time'] = (icoord, range(tdata.standard_time.shape[0]))
        # rename pressure
        tdata = tdata.rename_dims({'plev': 'pressure'})
        # swap standard time back to time
        tdata = tdata.swap_dims({icoord: 'time'})
        # Cube hour x pressure x time
        tdata = tdata.transpose("hour", "pressure", "time")
        # pressure values in hPa
        tdata['press'] = ('pressure', tdata['plev'].values / 100)
        # datum variable with days
        tdata['datum'] = ('time', tdata[icoord].values)
        tdata['datum'] = tdata['datum'].expand_dims({'numdat': 1})
        # remove other variables and coordinates
        tdata = tdata.drop(icoord + '_flag')
        tdata = tdata.drop(icoord)
        tdata = tdata.drop('plev').drop('hour', dim=None).drop('time', dim=None)
        # ends up wrong in file
        # tdata = tdata.set_coords(['press', 'datum'])
        # Rename variables
        if 'temperature' in tdata.variables:
            tdata = tdata.rename_vars({'temperature': 'temperatures'})
        # ERA5 variables
        if '{}_fg_depar'.format(variable) in tdata.variables:
            tdata = tdata.rename_vars(
                {'{}_fg_depar'.format(variable): 'era5_fgdep'})
        if '{}_biascorr'.format(variable) in tdata.variables:
            tdata = tdata.rename_vars({'{}_biascorr'.format(variable): 'bias'})
        if '{}_an_depar'.format(variable) in tdata.variables:
            tdata = tdata.rename_vars(
                {'{}_an_depar'.format(variable): 'era5_andep'})
        # bias adjustment variables
        for ivar in tdata.variables:
            if '_bias_estimate' in ivar:
                tdata = tdata.rename_vars({ivar:'bias'})
                tdata['bias'].attrs['long_name'] = source
        
        # Add station information (latest)
        tdata['lat'] = ('station', [self.header_table.latitude[-1]])
        tdata['lon'] = ('station', [self.header_table.longitude[-1]])
        # TODO: alt does not appear in files ?
        #try:
        #    tdata['alt'] = ('station', [self.header_table.height_of_station_above_sea_level[-1] ])
        #except:
            # tdata['alt'] = ('station', [0])
        #    pass

        # Add attributes from the default ones
        for ikey, ival in default_attrs.items():
            try:
                tdata[ikey].attrs.update(ival)
            except:
                pass

        # Add Global Attributes
        tdata.attrs.update(default_global)
        # Remove some Attributes:
        for ivar in tdata:
            # remove cdm stuff
            for iatt in list(tdata[ivar].attrs.keys()):
                if iatt in ['cdmcode', 'odbcode']:
                    del tdata[ivar].attrs[iatt]

        # Output to file or return as Dataset 
        if filename is not None:
            if os.path.isfile(filename):
                # write only variables, no dims (assuming dims are all the same anyway)
                # Leo: netcdf_classic needs to be as format
                with xr.open_dataset(filename, engine='netcdf4', format='NETCDF4_CLASSIC',mode = 'a' if not force else 'w') as fopen:
                    for ivar in tdata.variables:
                        if ivar in fopen.variables:
                            logger.debug("Skipping {}".format(ivar))
                        else:
                            fopen[ivar] = tdata[ivar]
                            logger.info("Writing {}".format(ivar))
                    fopen.attrs.update(default_global)
                logger.info("Appended to {}".format(filename))
            else:
                # write a new file
                tdata.to_netcdf(filename)
                logger.info("Written to {}".format(filename))
        else:
            return tdata
     

    def read_data_to_3dcube(self, variables: list, dates: list = None, plevs: list = None, feedback: list = None,
                          feedback_group: str = 'era5fb', **kwargs) -> dict:
        """ Read standard pressure levels and return a DataCube

        Args:
            variables: list of variables, e.g. temperature
            dates: [start, stop], str, int or datetime
            plevs: [list] in Pa or hPa
            feedback: list of feedback variables
            feedback_group: group name of the feedback
            **kwargs:

        Optional Keywords:
            date_time_name: Name of the datetime variable
            z_coordinate_name: Name of the pressure level variable

        Returns:
            dict : {variable : xr.DataArray}
        """
        if len(variables) == 0:
            raise ValueError('Need a variables', str(cdm_codes.keys()))

        if isinstance(variables, str):
            variables = [variables]

        data = {}
        if plevs is not None:
            plevs = np.asarray(plevs)
            if any((plevs > 110000) | (plevs < 500)):
                raise ValueError('Pressure levels outside range [5, 1100] hPa')
        else:
            plevs = std_plevs * 100  # in Pa

        std_plevs_indices = np.zeros(1001, dtype=np.int32)  # hPa
        # in hPa
        for i, j in enumerate(plevs // 100):
            std_plevs_indices[j] = i

        # todo check if variable can be replaced inside HDF5 ?
        if self.hasgroups:
            #
            # check variables
            #
            varnum = []
            for ivar in variables:
                if ivar not in cdm_codes.keys():
                    raise ValueError('Variable not found', ivar)
                varnum.append(
                    {'varnum': cdm_codes[ivar], 'variable': 'observation_value', 'group': 'observations_table',
                     'bkp_var': ivar})
                if feedback is not None:
                    if isinstance(feedback, str):
                        feedback = [feedback]
                    for jvar in feedback:
                        # jvar -> @body rename
                        varnum.append({'varnum': cdm_codes[ivar], 'variable': jvar, 'group': feedback_group,
                                       'bkp_var': ivar})
            # multiprocessing of the requests?
            for ivarnum in varnum:
                logger.info('Reading ... %d  %s', ivarnum['varnum'], ivarnum['bkp_var'])
                ivarnum.update(kwargs)
                # TIME SLICE, INDEX, SECONDS ARRAY, PRESSURE LEVELS
                trange, indices, secarray, pressure = self.read_observed_variable(dates=dates,
                                                                                  plevs=plevs,
                                                                                  return_coordinates=True,
                                                                                  return_index=True,
                                                                                  **ivarnum)
                logger.info('[CUBE] Variable Group %d %s %s', ivarnum['varnum'], str(trange), str(secarray.shape))
                obs = self[ivarnum['group']][ivarnum['variable']][trange][indices]
                #
                # to Cube
                #
                # requires hPa for indices
                itime, iobs = table_to_cube(secarray,
                                            std_plevs_indices[pressure.astype(np.int32) // 100],
                                            obs,
                                            nplev=plevs.size)
                logger.info('[CUBE] %s %s', ivarnum['bkp_var'], iobs.shape)
                v_attrs = get_attributes(cdmcode=ivarnum['varnum'],
                                         feedback=ivarnum['variable'] if feedback is not None else None)
                if len(v_attrs) > 0:
                    v_attrs = v_attrs[list(v_attrs.keys())[0]]
                # Convert to Xarray [time x plev]
                data[ivarnum['bkp_var']] = xr.DataArray(iobs,
                                                        coords=(seconds_to_datetime(secarray[itime]), plevs),
                                                        dims=('time', 'plev'),
                                                        name=ivarnum['bkp_var'],
                                                        attrs=v_attrs,
                                                        )
                # todo add attributes for coordinates
        else:
            splevs=[1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000,92500,
                    100000]
            ldaysx2=False
            for ivar in variables:
                # Read Attributes Variable
                if ivar=='hours':
                    continue
                v_attrs = self.read_attributes(ivar)
                # why no trange?
                iobs, secarray, pressure = self.read_variable(ivar,
                                                              dates=dates,
                                                              plevs=plevs,
                                                              return_coordinates=True)
                if ldaysx2:
                    cobs=np.full(cobs.shape,np.nan,dtype=iobs.dtype)
                    cobs[hindex,pindex,dindex]=iobs[hgood]
                else:
                    #daysx2(secarray,
                                                             #std_plevs_indices[pressure.astype(np.int32) // 100],
                                                             #plevs.shape[0],iobs)
                    cobs,hindex,pindex,dindex,hgood,gdays,ghours=daysx2(secarray,
                                                             std_plevs_indices[pressure.astype(np.int32) // 100],
                                                             plevs.shape[0],iobs)

                #itime, iobs = table_to_3dcube(secarray,
                                            #std_plevs_indices[pressure.astype(np.int32) // 100],
                                            #iobs)
                logger.info('[CUBE] %s %s %s', ivar, iobs.shape,cobs.shape)
                # Convert to Xarray [time x plev]
                data[ivar] = xr.DataArray(cobs,
                                          coords=(np.array((0,12)),plevs/100,gdays),
                                          dims=('hour', 'press','datum'),
                                          name=ivar,
                                          attrs=v_attrs)
                data[ivar]['datum'].attrs['units']='days since 1900-01-01 00:00:00'
                data[ivar]['press'].attrs['units']='hPa'
            
            data['hours']=xr.DataArray(ghours,name='hours',dims=('hour','datum'),coords=(np.array((0,12)),gdays))
            data['lat']=xr.DataArray(np.asarray([self.lat[0]]),name='lat',dims=('station'))
            data['lon']=xr.DataArray(np.asarray([self.lon[0]]),name='lon',dims=('station'))
                # data[ivar].data['time'].attrs.update(self.read_attributes(kwargs.get('date_time_name', 'time')))
                # data[ivar].data['plev'].attrs.update(self.read_attributes(kwargs.get('z_coordinate_name', 'plev')))
        return data

    def trajectory_data(self):
        """ Return trajectory data"""
        # should return data that can be plotted on a map ?
        # todo finish this function
        # return trajectory data with lon, lat, label
        if self.hasgroups:
            # what to read here? station_configuration ?
            pass
        else:
            # read trajectory information, each profile has exectly one position, how much data is there per profile?
            pass

    def write_observed_data(self, name: str, varnum: int,
                            cube: xr.DataArray = None,
                            ragged=None,
                            group: str = None,
                            variable: str = 'observation_value',
                            dimgroup: str = 'observations_table',
                            data_time: str = 'time',
                            data_plevs: str = 'plev',
                            force_replace: bool = False,
                            attributes: dict = None,
                            interpolate: bool = False,
                            interpolate_datetime: bool = False,
                            extrapolate_time: bool = True,
                            extrapolate_plevs: bool = False,
                            **kwargs):
        """ Write a DataCube or a Table (DataArray, DataFrame) with Multiindex to CDM Backend file

        Args:
            name: name of variable to write to
            cube: Cube input data
            ragged: Table input data
            varnum: observed variable CDM Code
            group: new Group
            variable: Variable Name
            dimgroup: Variable group
            data_time: time dimension input name
            data_plevs: pressure dimension input name
            force_replace: overwrite string array due to shape missmatch?
            attributes: add Attributes to HDF variable
            interpolate: interpolate to all pressure levels available in the file
            interpolate_datetime: interpolate values backwards and forward in time?
            extrapolate_plevs: when interpolating, extrapolate outside of given values? constant
            extrapolate_time: when interpolating, extrapolate outside of given dates? constant

        Examples:
            Write an adjusted temperature to a CDM Backend file within a new group mytests that is align
            with 85 temperature. data is a ragged table
            >>> data = self.read_observed_variable(85)
            >>> write_observed_data('temperature_adjust', ragged=data, varnum=85, variable='observation_value',
            >>>               group='mytests')

        Notes:
            Interpolation
                When interpolating one should be careful. Interpolation for pressure is done using log-p.
                using interpolate_na from xarray

            Extrapolation
                using backward-fill and forward-fill for the dimension (time, plev), thus constant outside of range.
        """
        if not self.hasgroups:
            raise RuntimeError('This routine is intended for CDM Backend files')

        if dimgroup not in self.groups:
            raise ValueError('Dimgroup not found:', dimgroup)
        else:
            if variable not in self[dimgroup].keys():
                raise ValueError('Variable', variable, 'not in dimgroup', dimgroup)

        if cube is None and ragged is None:
            raise ValueError('Requires a Cube or a ragged Array')
        #
        #  CUBE (time x plev)
        #
        if cube is not None:
            if data_time not in cube.dims:
                raise ValueError('Datetime dimension not found', data_time)

            if data_plevs not in cube.dims:
                raise ValueError('Z coordinate not found', data_plevs)

            data = cube.stack(obs=(data_time, data_plevs)).dropna('obs')  # need to remove NaN
        #
        # Ragged Array (index (time x plev))
        #
        if ragged is not None:
            if isinstance(ragged, pd.DataFrame):
                ragged = ragged.to_xarray()

            # todo check this, assumes this is a Multi-Index Table, Xarray
            data = ragged  # is allready a table
            idim = data.dims[0]
            data = data.rename(**{idim: 'obs'})  # make sure the name is fixed
        #
        # Convert to ragged array, index should be sorted?
        #
        if not data.indexes['obs'].is_monotonic:
            logger.warning('Data not sorted !!!')
        #
        # Get data coordinate information
        #
        try:
            in_dates = data[data_time].values  # Probably np.datetime64
        except:
            raise ValueError('Datetime dimension not found', data_time)

        try:
            in_plevs = data[data_plevs].values.astype(np.int32)
        except:
            raise ValueError('Z coordinate not found', data_plevs)
        #
        # Make datetime slice
        #
        slice_dates = [in_dates[0], in_dates[-1]]  # for trange (slice)
        if isinstance(in_dates[0], np.datetime64):
            # to seconds since
            in_dates = ((in_dates - np.datetime64('1900-01-01T00:00:00')) / np.timedelta64(1, 's')).astype(np.int64)
        #
        # Reopen file for writing
        #
        if 'r' == self.file.mode:
            self.reopen(mode='r+', **kwargs)
        else:
            logger.debug('File mode: [%s] %s', self.file.mode, self.name)
        #
        # Read coordinate informations
        # - trange      - Make a time-slice
        # - mask        - a observed variable subset
        # - f_dates     - datetime[mask] in s
        # - f_plevs     - p-levels[mask]
        trange, mask, f_dates, f_plevs = self.read_observed_variable(varnum,
                                                                     dates=slice_dates,
                                                                     return_index=True,
                                                                     return_coordinates=True,
                                                                     **kwargs)
        #
        # Do we need to sort datetime from file ? (Potential src of error)
        #
        reverse_sort = False
        if not is_sorted(f_dates):
            if False:
                idx = np.argsort(f_dates)
                reverse_sort = True
                # todo see if reverse sort is necessary?
                f_dates = f_dates[idx]
                f_plevs = f_dates[idx]
                idx = np.argsort(idx)  # -> reverse sort index
            logger.warning('%s not sorted in file', data_time)
        #
        # Find reverse index for writing
        #
        match_index = np.full(f_dates.shape, -1, dtype=np.int32)  # Matching indices
        # Loop input and file dates/pressures -> Matches
        reverse_index(match_index, f_dates, f_plevs, in_dates, in_plevs)
        logic = (match_index > -1)
        if logic.sum() != in_dates.shape[0]:
            missing = in_dates[~np.in1d(in_dates, in_dates[match_index[logic]])]
            logger.warning('Not all input dates are found %s', missing.shape)
            logger.debug('Missing dates: %s', str(missing))
        else:
            logger.info('All dates found')

        values = data.values
        if interpolate and not np.issubdtype(values.dtype, str):
            mask = np.where(mask)[0]  # every observed value
            # create array [trange] -> all levels
            fvalues = np.full(f_dates.shape, np.nan, dtype=values.dtype)  # not compatible with interpolation
            # fill in data that is now sorted like in the file
            fvalues[logic] = values[match_index[logic]]
            logger.info('Interpolation: %d < %d', f_dates.shape[0], np.isfinite(values).sum())
            # interpolate data to missing parts
            farray = xr.DataArray(fvalues, dims=('obs'),
                                  coords={data_time: ('obs', f_dates), data_plevs: ('obs', f_plevs)})
            if interpolate_datetime:
                #
                # Interpolate back/forward in time
                # extrapolate=True
                #
                subset = farray[data_plevs].isin(std_plevs * 100)  # needs to be in Pa?
                farray[subset] = farray[subset].groupby(data_plevs).apply(level_interpolation,
                                                                          dim=data_time,
                                                                          extrapolate=extrapolate_time)
            #
            # Interpolate up/down levels
            #
            farray = farray.groupby(data_time).apply(level_interpolation, dim=data_plevs, method='log-linear',
                                                     extrapolate=extrapolate_plevs)
            values = farray.values
            logger.info('Post-Interpolation: %d < %d', f_dates.shape[0], np.isfinite(values).sum())
            #
            # logic / match_index are useless now (replace with dummys)
            #
            logic = np.ones(values.shape[0], dtype=np.bool)
            match_index = np.arange(values.shape[0], dtype=np.int)

        else:
            mask = np.where(mask)[0][logic]
        #
        # Get Group / New Group
        #
        if group in self.groups:
            gid = self.file[group]  # Append to existing group
        else:
            gid = self.file.create_group(group)  # New Group
        #
        # Variable new or old
        #
        if name in gid.keys():
            #
            # String Arrays are different
            #
            if np.issubdtype(values.dtype, str):
                n = len(values[0])
                m = gid[name].shape[0]
                # Shapes of file and input do not match?
                if gid[name].shape[1] != n:
                    if force_replace:
                        #
                        # Replace existing variable
                        #
                        logger.info('Removing %s/%s (%s) with (%d, %d)', group, name, gid[name].shape,
                                    values.shape, n)
                        gid.create_dataset_like('temp123',
                                                self.file[group][name],
                                                dtype='S{}'.format(n),
                                                shape=(m, n),
                                                chunks=True)  # Create a string array dataset
                        #
                        # Remove old variable and rename
                        #
                        del gid[name]  # remove in file
                        del self[group][name]  # remove in class
                        gid[name] = gid['temp123']  # overwrite
                        del gid['temp123']  # remove temporary dataset
                        #
                        # String dimensionwith with size (n)
                        #
                        sname = 'string{}'.format(n)
                        if sname not in gid.keys():
                            gid.create_dataset(sname, data=np.zeros(n, dtype='S1'), chunks=True)
                            gid[sname].attrs['NAME'] = np.string_(
                                'This is a netCDF dimension but not a netCDF variable.')
                            gid[sname].make_scale(sname)
                            gid[name].dims[1].attach_scale(gid[sname])
                        # TODO check if old string dimension is still used or not ?
                        #
                        # Attach dimensions from file
                        #
                        try:
                            # TODO there might be conflict if dimensions are used from a different group
                            # this is not accounted here (beacuse of our structure) 
                            # e.g. [/observations_table/index, /observations_table/string11]
                            idim = [jdim.name for jdim in self.file[dimgroup][variable].dims.keys()][0]
                        except:
                            idim = ''
                        idim = idim.split('/')[-1] if '/' in idim else idim
                        idim = idim if idim != '' else 'index'
                        if idim not in gid.keys():
                            gid[idim] = self.file[dimgroup][idim]
                            gid[idim].make_scale()
                        gid[name].dims[0].attach_scale(gid[idim])
                    else:
                        assert gid[name].shape[1] == n, 'Shapes do not match, force_replace=True to overwrite '
                    self.file.flush()  # write changes to file
                writeme = gid[name][()]  # get all data
                chardata = values.astype('S{}'.format(n)).view('S1').reshape((m, n))
                writeme[trange, :][mask, :] = chardata[match_index[logic], :]  # fill Array
            else:
                #
                # load exisiting data
                #
                writeme = self.load_variable_from_file(name, group=group, return_data=True)[0]
                writeme[trange][mask] = values[match_index[logic]]  # Write the data into the file
        else:
            #
            # create a new dataset
            #
            logger.info('Creating Dataset %s (%s/%s, %s)', name, dimgroup, variable,
                        self.file[dimgroup][variable].shape)
            if np.issubdtype(values.dtype, str):
                n = len(values[0])
                m = values.shape[0]
                chararray = values.astype('S{}'.format(n)).view('S1').reshape((m, n))
                gid.create_dataset(name,
                                   dtype='S{}'.format(n),
                                   shape=(self.file[dimgroup][variable].shape[0], n),
                                   chunks=True,
                                   compression='gzip')
                writeme = np.zeros((self.file[dimgroup][variable].shape[0], n), dtype='S1')
                writeme[trange, :][mask, :] = chararray[match_index[logic], :]  # fill Array
                #
                # String dimensionwith with size (n)
                #
                sname = 'string{}'.format(n)
                if sname not in gid.keys():
                    gid.create_dataset(sname, data=np.zeros(n, dtype='S1'), chunks=True)
                    gid[sname].attrs['NAME'] = np.string_('This is a netCDF dimension but not a netCDF variable.')
                    gid[sname].make_scale(sname)
                    gid[name].dims[1].attach_scale(gid[sname])
            else:
                #
                # Different fillvalues for float, integer
                #
                if np.issubdtype(values.dtype, int):
                    fillvalue = np.int32(-2147483648)
                else:
                    fillvalue = np.float32(np.nan)
                #
                # create a new dataset
                #
                gid.create_dataset_like(name, self.file[dimgroup][variable])  # Create a new dataset
                writeme = np.full(self.file[dimgroup][variable].shape,
                                  fillvalue,
                                  dtype=self.file[dimgroup][variable].dtype)
                writeme[trange][mask] = values[match_index[logic]]  # Write the data into the file
            #
            # Attach dimensions from file
            #
            try:
                idim = self.file[dimgroup][variable].dims[0].keys()[0]
            except:
                idim = ''
            idim = idim if idim != '' else 'index'
            if idim not in gid.keys():
                gid[idim] = self.file[dimgroup][idim]
                gid[idim].make_scale()
            gid[name].dims[0].attach_scale(gid[idim])

        #
        # Write and finish
        #
        gid[name][()] = writeme  # Write the new data
        #
        # Write Attributes
        #
        if attributes is not None:
            for ikey, ival in attributes.items():
                if isinstance(ival, str):
                    gid[name].attrs[ikey] = np.string_(ival)
                else:
                    gid[name].attrs[ikey] = ival
        self.file.flush()  # Put changes into the file
        self.inquire()  # update Groups and variable lists
        self[group][name].update(data=writeme)  # update class in memory
        logger.info('Finsihed writing %s to %s', name, self.name)

    def write_variable(self, name: str,
                       cube: xr.DataArray = None,
                       ragged=None,
                       data_time: str = 'time',
                       data_plevs: str = 'plev',
                       force_replace: bool = False,
                       attributes: dict = None,
                       interpolate: bool = False,
                       interpolate_datetime: bool = False,
                       extrapolate_time: bool = True,
                       extrapolate_plevs: bool = False,
                       **kwargs):
        """ Write a DataCube or a Table (DataArray, DataFrame) with Multiindex to CDM Frontend

        Args:
            name: name of variable to write to
            cube: Cube input data
            ragged: Table input data
            data_time: time dimension input name
            data_plevs: pressure dimension input name
            force_replace: overwrite string array
            attributes: data attributes to be written

        Optional Keywords:
            date_time_name: usually time
            z_coordinate_name: usually plev

        Examples:
            Write an adjusted temperature to a CDM Frontend file. data is a ragged table
            >>> data = self.read_variable('ta')
            >>> write_variable('temperature_adjust', ragged=data, variable='ta')
        """
        if self.hasgroups:
            raise RuntimeError('Only CDS frontend files')

        if cube is None and ragged is None:
            raise ValueError('Requires a Cube or a ragged Array')
        #
        #  CUBE (time x plev)
        #
        if cube is not None:
            if data_time not in cube.dims:
                raise ValueError('Datetime dimension not found', data_time)

            if data_plevs not in cube.dims:
                raise ValueError('Z coordinate not found', data_plevs)

            data = cube.stack(obs=(data_time, data_plevs)).dropna('obs')  # need to remove NaN
        #
        # Ragged Array (index (time x plev))
        #
        if ragged is not None:
            if isinstance(ragged, pd.DataFrame):
                ragged = ragged.to_xarray()

            # todo check this, assumes this is a Multi-Index Table, Xarray
            data = ragged  # is allready a table
            idim = data.dims[0]
            data = data.rename(**{idim: 'obs'})  # make sure the name is fixed

        if attributes is None:
            attributes = data.attrs.copy()
        #
        # Convert to ragged array, index should be sorted?
        #
        if not data.indexes['obs'].is_monotonic:
            logger.warning('Data not sorted !!!')
        #
        # Get data coordinate information
        #
        try:
            in_dates = data[data_time].values  # Probably np.datetime64
        except:
            raise ValueError('Datetime dimension not found', data_time)

        try:
            in_plevs = data[data_plevs].values.astype(np.int32)
        except:
            raise ValueError('Z coordinate not found', data_plevs)
        #
        # Make datetime slice
        #
        slice_dates = [in_dates[0], in_dates[-1]]  # for trange (slice)
        if isinstance(in_dates[0], np.datetime64):
            # to seconds since
            in_dates = ((in_dates - np.datetime64('1900-01-01T00:00:00')) / np.timedelta64(1, 's')).astype(np.int64)
        #
        # Reopen file for writing
        #
        if 'r' == self.file.mode:
            self.reopen(mode='r+', tmpdir=kwargs.get('tmpdir', '~/tmp'))
        else:
            logger.debug('File mode: [%s] %s', self.file.mode, self.name)
        #
        # Read coordinate information from file
        #
        idate_name = kwargs.get('date_time_name', 'time')
        trange = self.make_datetime_slice(dates=slice_dates, **kwargs)
        f_dates = self.load_variable_from_file(idate_name, return_data=True)[0]
        f_plevs = self.load_variable_from_file(kwargs.get('z_coordinate_name', 'plev'), return_data=True)[0]
        #
        # Do we need to sort datetime from file ? (Potential src of error)
        #
        reverse_sort = False
        if not is_sorted(f_dates):
            if False:
                idx = np.argsort(f_dates)
                reverse_sort = True
                # todo see if reverse sort is necessary?
                f_dates = f_dates[idx]
                f_plevs = f_dates[idx]
                idx = np.argsort(idx)  # -> reverse sort index
            logger.warning('%s not sorted in file', idate_name)
        #
        # Find reverse index for writing
        #
        # Match Input datetime and p-levs with the ones in the file
        match_index = np.full(f_dates.shape, -1, dtype=np.int32)  # Matching indices
        # Loop input and file dates/pressures -> Matches
        reverse_index(match_index, f_dates, f_plevs, in_dates, in_plevs)
        logic = (match_index > -1)
        # Might give 0 missing, because all dates have been found, but not all pressure levels (e.g. from a cube)
        if logic.sum() != in_dates.shape[0]:
            missing = in_dates[~np.in1d(in_dates, in_dates[match_index[logic]])]
            logger.debug('Not all input dates are found %s', missing.shape)
            logger.debug('Missing dates: %s', str(missing))
        else:
            logger.debug('All dates found %d', logic.sum())

        values = data.values
        if interpolate and not np.issubdtype(values.dtype, str):
            mask = np.arange(trange.start, trange.stop, dtype=np.int)
            # create array of file dimensions
            fvalues = np.full(f_dates.shape, np.nan, dtype=values.dtype)  # not compatible with interpolation
            # fill in data
            fvalues[logic] = values[match_index[logic]]
            logger.info('Interpolation: %d < %d', f_dates.shape[0], np.isfinite(values).sum())
            # interpolate data to missing parts
            farray = xr.DataArray(fvalues, dims=('obs'),
                                  coords={data_time: ('obs', f_dates), data_plevs: ('obs', f_plevs)})
            if interpolate_datetime:
                #
                # Interpolate back/forward in time
                # extrapolate=True
                #
                subset = farray[data_plevs].isin(std_plevs * 100)  # needs to be in Pa?
                farray[subset] = farray[subset].groupby(data_plevs).apply(level_interpolation,
                                                                          dim=data_time,
                                                                          extrapolate=extrapolate_time)
            #
            # Interpolate up/down levels
            #
            farray = farray.groupby(data_time).apply(level_interpolation,
                                                     dim=data_plevs,
                                                     method='log-linear',
                                                     extrapolate=extrapolate_plevs)
            values = farray.values
            logger.info('Post-Interpolation: %d < %d', f_dates.shape[0], np.isfinite(values).sum())
            #
            # logic / match_index are useless now (replace with dummys)
            #
            logic = np.ones(values.shape[0], dtype=np.bool)
            match_index = np.arange(values.shape[0], dtype=np.int)
        else:
            mask = np.where(logic)[0]
        #
        # String Arrays are different
        #
        if np.issubdtype(values.dtype, str):
            n = len(values[0])
            m = values.shape[0]
            if name in self.groups:
                # shapes of file and input do not match?
                if self[name].shape[1] != n:
                    if force_replace:
                        # TODO update with copy_h5py_variable code (which is much more general)
                        #
                        # Replace existing variable
                        #
                        logger.info('Removing %s (%s) with (%d,%d)', name, self[name].shape, m, n)
                        self.file.create_dataset_like('temp123',
                                                      self.file[name],
                                                      dtype='S{}'.format(n),
                                                      shape=(m, n))  # Create a string array dataset
                        #
                        # Remove old variable and rename
                        #
                        del self.file[name]  # remove in file
                        del self.file[name]  # remove in class
                        self.file[name] = self.file['temp123']  # overwrite
                        del self.file['temp123']  # remove temporary dataset
                        #
                        # Dimension names have been removed (phony_dim_x)
                        #
                        if False:
                            #
                            # String dimension with size (n)
                            #
                            sname = 'string{}'.format(n)
                            if sname not in self.file.keys():
                                self.file.create_dataset(sname, data=np.zeros(n, dtype='S1'), chunks=True)
                                self.file[sname].attrs['NAME'] = np.string_(
                                    'This is a netCDF dimension but not a netCDF variable.')
                                self.file[sname].make_scale(sname)
                                self.file[name].dims[1].attach_scale(self.file[sname])
                            # todo check if old string dimension is still used or not ?
                            #
                            # Attach dimensions from file
                            #
                            idim = self.file[name].dims[0].keys()[0]
                            idim = idim if idim != '' else 'obs'
                            if idim not in self.file.keys():
                                self.file[idim] = self.file[idim]
                                self.file[idim].make_scale()
                            self.file[name].dims[0].attach_scale(self.file[idim])
                    else:
                        assert self[name].shape[1] == n, 'Shapes do not match, force_replace=True to ' \
                                                         'overwrite '
                    self.file.flush()  # write changes to file
                writeme = self.file[name][()]  # get all data
                chardata = values.astype('S{}'.format(n)).view('S1').reshape((m, n))
                # writeme = np.full((self[name].shape[0], n), '', dtype='S1')
                writeme[mask, :] = chardata[match_index[logic], :]  # fill Array
            else:
                #
                # create a new dataset
                #
                chararray = values.astype('S{}'.format(n)).view('S1').reshape((m, n))
                self.file.create_dataset(name,
                                         dtype='S{}'.format(n),
                                         shape=(f_dates.shape[0], n),
                                         chunks=True,
                                         compression='gzip')
                writeme = np.zeros((f_dates.shape[0], n), dtype='S1')
                writeme[mask, :] = chararray[match_index[logic], :]  # fill Array
                #
                # Dimension names have been removed (phony_dim_x)
                # 
                if False:
                    #
                    # String dimensionwith with size (n)
                    #
                    sname = 'string{}'.format(n)
                    if sname not in self.file.keys():
                        self.file.create_dataset(sname, data=np.zeros(n, dtype='S1'), chunks=True)
                        self.file[sname].attrs['NAME'] = np.string_('This is a netCDF dimension but not a netCDF variable.')
                        self.file[sname].make_scale(sname)
                        self.file[name].dims[1].attach_scale(self.file[sname])
                    #
                    # Attach dimensions from file
                    #
                    idim = self.file[name].dims[0].keys()[0]
                    idim = idim if idim != '' else 'obs'
                    if idim not in self.file.keys():
                        self.file[idim] = self.file[idim]
                        self.file[idim].make_scale()
                    self.file[name].dims[0].attach_scale(self.file[idim])
        #
        # Just for numbers
        #
        else:
            #
            # Different fillvalues for float, integer
            #
            if np.issubdtype(values.dtype, int):
                fillvalue = np.int32(-2147483648)
            else:
                fillvalue = np.float32(np.nan)

            if name in self.groups:
                #
                # load exisiting data
                #
                logger.warning('Overwriting %s | %s', name, self.name)
                writeme = self.load_variable_from_file(name, return_data=True)[0]
            else:
                #
                # create a new dataset
                #
                logger.info('Creating Dataset %s (%s)', name, values.shape)
                self.file.create_dataset(name,
                                         dtype=values.dtype,
                                         shape=f_dates.shape,
                                         chunks=True,
                                         fillvalue=fillvalue,
                                         compression='gzip')  # Create a new dataset
                #
                # Dimension names have been removed (phony_dim_x)
                #
                if False:
                    #
                    # Attach dimensions from file
                    #
                    try:
                        idim = self.file[name].dims[0].keys()[0]
                    except:
                        idim = ''
                    idim = idim if idim != '' else 'obs'
                    if idim not in self.file.keys():
                        self.file[idim] = self.file[idim]
                        self.file[idim].make_scale()
                    self.file[name].dims[0].attach_scale(self.file[idim])
                writeme = np.full(f_dates.shape,
                                  fillvalue,
                                  dtype=values.dtype)

            writeme[mask] = values[match_index[logic]]  # Write the data
        #
        # Write and finish
        #
        self.file[name][()] = writeme  # Write the new data
        #
        # Write Attributes (can cause: can't open HDF5 attribute)
        #
        if attributes is not None:
            for ikey, ival in attributes.items():
                self.file[name].attrs[ikey] = np.string_(str(ival))
                #if isinstance(ival, str):
                #    self.file[name].attrs[ikey] = np.string_(ival)
                #else:
                #    self.file[name].attrs[ikey] = ival
        self.file.flush()  # Put changes into the file
        self.inquire()  # update Groups and variable lists
        self[name].update(data=writeme)  # update class in memory
        logger.info('Finsihed writing %s to %s', name, self.name)

#     def write_to_frontend_file(self, filename: str, variables: list, index: np.ndarray = None, cf_dict: dict = None,
#                                **kwargs):
#         """ Write a

#         Args:
#             filename:
#             variables:
#             index:
#             cf_dict:
#             **kwargs:

#         Examples:
#             >>> data = CDMDataset('testfile.nc')
#             >>> data.write_to_frontend_file('only_wind.nc', ['wind_speed', 'wind_from_direction'])

#         """
#         # search for variables if cdmname not in cf_dict
#         if not self.hasgroups:
#             raise NotImplementedError('Not yet implemented')

#         # make a decision which dimension requires what
#         # observations_table -> obs
#         # header_table -> trajectory
#         # stations_configuration -> trajectory

#         # trajectory index?
#         # do_cfcopy

#         # global attributes
#         pass

    #
    # finding country for coordinates:
    #
    def get_address_by_location(self, latitude, longitude, language="en"):
        """This function returns an address as raw from a location
        will repeat until success"""
        app = geopy.geocoders.Nominatim(user_agent="eua_cds3")
        # build coordinates string to pass to reverse() function
        coordinates = f"{latitude}, {longitude}"
        # sleep for a second to respect Usage Policy
        time.sleep(1)
        try:
            return app.reverse(coordinates, language=language).raw
        except:
            return self.get_address_by_location(latitude, longitude)
        

    def report_quality(self, filename: str = None, **kwargs):
        """ Compile a quality report

        Args:
            filename:

        Returns:

        Plans:
            Station : ID
            Region: EUR
            Country: AUT
            Location : lon, lat (most of the time)
                      other lon, lat combinations

            Merged Stations : IDs
            Distance between Stations : 1.2 km ID
            Temporal Coverage: [Start, End]
            before 1940: xxx
            before 1950: xxx
            before 1979: xxx
            after 1979: xxx
            Current Status: active / inactive ?
            Launch times: 00 (xxx), 12 (xxx), 6 (xxx)
            Observations: variables (xxx)
            Temperature Breakpoinst: from (adjust)
            Humidity Breakpoints:
            Wind Breakpoints:
            Temperature Outliers: (based on climatology)
            Humidity Outliers:
            Wind Outliers:
        """
        # run a quality control procedure to check if the CDM file is within the standard
        # run a duplicated controller
        # run a values consistency check
        # Calculate climatology and outlier statistics
        
        writerep = [] # list for writing report

        
        if not self.hasgroups:
            raise RuntimeError('Only available for CDM Backend files')

        report = {'Station': '', 'Country': '', 'State': '', 'Latitude': [], 'Longitude': [], 
                  'Temporal_Coverage':'', 'Ascents':'', 'Observations':'', 'b1940':'', 'b1950':'', 'b1979':'' , 'a1979':'', 'l00':'', 'l06':'', 'l12':'', 'l18':'', 
                  # 'Merged_Stations':[], 
                  'Merged_Observations': [],  'Merged_Sources':[],
                  'Variable_Code':[], 'Variables':[], 'Missing_Values': [], 'Min_Value':[], 'Max_Value':[], 
                  'PLvl':[85000,], 'PLvl_Observations':[], 'Outliers':[],}

        igroup = 'header_table'
        if igroup in self.groups:
            #
            # station id:
            #
            rdvar = None
            try:
                rdvar = self['header_table']['primary_station_id'][()][0]
            except:
                # not available right now - will be in 'primary_station_id' in the header_table or 'primary_id' in the station_configuration
                pass
            if rdvar is not None:
                wrtvar = ''
                for i in rdvar:
                    wrtvar = wrtvar + i.decode()
                report['Station'] = wrtvar
                writerep.append('station_id: '+ report['Station'])
            else:
                writerep.append('station_id: not readable')
            #
            # location:
            #
            lat = self['header_table']['latitude'][()][0]
            lon = self['header_table']['longitude'][()][0]
            writerep.append('latitude: '+ str(lat))
            report['Latitude'] = lat
            writerep.append('longitude: '+ str(lon))
            report['Longitude'] = lon
#             location = self.get_address_by_location(lat, lon)
#             report['Country'] = location['address']['country']
#             writerep.append('Country: '+ report['Country'])
#             report['State'] = location['address']['state']
#             writerep.append('State: '+ report['State'])
        else:
            writerep.append('header_table: not readable')
        #
        # time:
        #
        try:
            startt = datetime.utcfromtimestamp(self.recordindices.recordtimestamp[()][0]-2208988800)
            endt = datetime.utcfromtimestamp(self.recordindices.recordtimestamp[()][-1]-2208988800)
            if endt >= startt:
                report['Temporal_Coverage'] = (str(startt.year) + '-' + str(startt.month) + '-' + str(startt.day) + ' - ' + str(endt.year) + '-' + str(endt.month) + '-' + 
                                               str(endt.day))
                writerep.append('Temporal Coverage:' + report['Temporal_Coverage'])
        except:
            writerep.append('Temporal Coverage: not readable')

            
        igroup = 'observations_table'
        if igroup in self.groups:
            # converting self into dataframe
            data = self.to_dataframe('observations_table', ['observation_value', 'date_time', 'z_coordinate', 'observed_variable', 'report_id', 'source_id'], None, **kwargs)
            #
            # Ascents & Observations
            #
            try:
                asc = data.report_id.drop_duplicates() 
                report['Ascents'] = len(asc)
                writerep.append('Over all Ascents: ' + str(report['Ascents']))
                report['Observations'] = len(data.report_id)
                writerep.append('Over all Observations: ' + str(report['Observations']))
            except:
                writerep.append('Ascents and Observations: not readable')
            #
            # period:
            #
            try: 
                report['b1940'] = len(data[data.date_time.dt.year < 1940])
                report['b1950'] = len(data[data.date_time.dt.year < 1950])
                report['b1979'] = len(data[data.date_time.dt.year < 1979])
                report['a1979'] = len(data[data.date_time.dt.year >= 1979])
                writerep.append('Observations before 1940: ' + str(report['b1940']))
                writerep.append('Observations before 1950: ' + str(report['b1950']))
                writerep.append('Observations before 1979: ' + str(report['b1979']))
                writerep.append('Observations after 1979: ' + str(report['a1979']))
            except:
                writerep.append('Observations in Period: not readable')
            #
            # Launch times:
            #
#             try:
            report['l00'] = len(data[data.date_time.dt.hour >= 0][data.date_time.dt.hour < 6])
            report['l06'] = len(data[data.date_time.dt.hour >= 6][data.date_time.dt.hour < 12])
            report['l12'] = len(data[data.date_time.dt.hour >= 12][data.date_time.dt.hour < 18])
            report['l18'] = len(data[data.date_time.dt.hour >= 18][data.date_time.dt.hour < 24])
            writerep.append('Sonde starts at 00 UTC: ' + str(report['l00']))
            writerep.append('Sonde starts at 06 UTC: ' + str(report['l06']))
            writerep.append('Sonde starts at 12 UTC: ' + str(report['l12']))
            writerep.append('Sonde starts at 18 UTC: ' + str(report['l18']))
#             except:
#                 writerep.append('Sonde starts: not readable')
            #
            # Merge:
            #
            try:
#                 merge_stations = len(data.report_id.drop_duplicates())
#                 report['Merged_Stations'].append(merge_stations)
#                 writerep.append('Merged Stations: ' + str(merge_stations))
                for i in data.source_id.drop_duplicates():
                    merge_data = data[data.source_id == i]
                    merge_name = i
                    report['Merged_Observations'].append(len(merge_data))
                    report['Merged_Sources'].append(merge_name)
                    writerep.append('Observations merged from ' + merge_name + ': ' + str(len(merge_data)))
            except:
                writerep.append('Observations merged from: not readable')
            #
            # Value Check:
            #
            try:
                for i in data.observed_variable.drop_duplicates(): 
                    report['Variable_Code'].append(i)
                    vcheck = data[data.observed_variable == i]
                    report['Variables'].append(len(vcheck))
                    report['Missing_Values'].append(len(vcheck[vcheck.observed_variable == np.nan].observation_value)) # -> is it rly np.nan, or another value?
                    report['Min_Value'].append(np.nanmin(vcheck.observation_value))
                    report['Max_Value'].append(np.nanmax(vcheck.observation_value))                    
                writerep.append('Variable Codes: ' + str(report['Variable_Code']))
                writerep.append('Variables: ' + str(report['Variables']))
                writerep.append('Missing Values: ' + str(report['Missing_Values']))
                writerep.append('Minimum Values: ' + str(report['Min_Value']))
                writerep.append('Maximum Values: ' + str(report['Max_Value']))                
            except:
                writerep.append('Value Check: not readable')
        else:
            writerep.append('observations_table: not readable')
        #
        # check for Outliers
        #
        min_p = 25
        max_p = 75
        cut = 2
        skewed= False
        for j in report['PLvl']:
            for i in data.observed_variable.drop_duplicates():
                outl = 0
                da = data[data.observed_variable == i]
                da = da[da.z_coordinate == j]
                da = np.array(da.observation_value)

                q_min, q_max = np.nanpercentile(da, min_p), np.nanpercentile(da, max_p)
                cut_off = (q_max - q_min) * cut
                lower, upper = q_min-cut_off, q_max+cut_off

                if skewed==True:
                    q50 = np.nanpercentile(da, 50)
                    lower , upper = q_min-(q50-q_min)*cut ,  q_max+(q_max-q50)*cut  # the higher the cut, the more relaxed the contition for exclusion 

                median = np.nanmedian(da)
                cleaned, outliers = [],[]

                for d in np.asarray(da):
                    if not(d >= lower and d <= upper):
                        outl += 1
                report['Outliers'].append(outl)
                report['PLvl_Observations'].append(len(da))
            writerep.append('Pressurelevel: ' + str(j))
            writerep.append('Observations on Pressurelevel: ' + str(report['PLvl_Observations']))
            writerep.append('Outliers on Pressurelevel: ' + str(report['Outliers']))
                

#         igroup = 'observations_table'
#         if igroup in self.groups:
#             # Load groups (full arrays)
#             self.load_variable_from_file(['observed_variable', 'units', 'observation_value',
#                                           'z_coordinate', 'date_time'], group=igroup)
#             #
#             # observed_variables
#             #
#             varcodes = self[igroup]['observed_variable'][()]
            
#             unique_varcodes = np.unique(varcodes)
#             for ivar in unique_varcodes:
#                 _, idx = self.read_observed_variable(ivar, return_index=True)
#                 # Check units
#                 units = self[igroup]['units'][idx]
#                 # Check observation_value
#                 obs = self[igroup]['observation_value'][idx]
#                 # only on standard pressure levels ?

#         if 'era5fb' in self.groups:
#             pass

#         # if 'station_'
#         if filename is not None:
#             pass
        
        with open("quality_report_"+report['Station']+".txt", "w") as text_file:
            for i in writerep:
                print(i, file=text_file)
        
        return report

    def check_cdm_duplicates(self, groups, variables,
                             observed_variable_name: str = 'observed_variable',
                             date_time_name: str = 'date_time',
                             z_coordinate_name: str = 'z_coordinate',
                             **kwargs):
        """ Run a duplication check on a group and some variables

        Args:
            groups: CDMDataset group
            variables: CDMDataset group variables
            observed_variable_name:
            date_time_name:
            z_coordinate_name:
            **kwargs:

        Returns:

        Examples:
            >>> data = CDMDataset('0-20000-0-01001_CEUAS_merged_v0.nc')
            >>> data.check_cdm_duplicates('observations_table', ['observation_value', 'date_time', 'z_coordinate', 'observed_variable', 'source_id', 'z_coordinate_type'])

        """

        if not self.hasgroups:
            raise RuntimeError('Only for CDM Backend files')

        if observed_variable_name not in variables:
            variables.append(observed_variable_name)
        if date_time_name not in variables:
            variables.append(date_time_name)
        if z_coordinate_name not in variables:
            variables.append(z_coordinate_name)

        data = self.to_dataframe(groups, variables, None, **kwargs)
        logger.info("Evaluating for duplicates ... (pandas)")
        data = data[
            data.duplicated([date_time_name, z_coordinate_name, observed_variable_name], keep=False)].sort_values(
            [date_time_name, z_coordinate_name])
        return data

    def check_cdm_dimensions(self):
        """Check the CDM Backend file if it has proper dimensions attached to variables.
    
        """
        if not self.hasgroups:
            raise RuntimeError('Only for CDM Backend files')

        def check_dims(fin):
            collect = {}
            for ivar in fin.keys():
                if isinstance(fin[ivar], h5py.Group):
                    collect.update(check_dims(fin[ivar]))
                else:
                    try:
                        # get dimension names
                        collect[fin[ivar].name] = [ jdim[0].name for jdim in fin[ivar].dims.keys() ]
                    except:
                        # if error, than it is a dimension
                        collect[fin[ivar].name] = 'dim'
            return collect

        collect = check_dims(self.file)
        dims = []
        for ivar,ival in collect.items():
            if isinstance(ival, str):
                dims.append(ivar)
                
        dims = sorted(list(set(dims)))
        found = [False]*len(dims)
        for ivar,ival in collect.items():
            if isinstance(ival, list):
                for idim in ival:
                    if idim in dims:
                        found[dims.index(idim)] = True
            # TODO check if variables have proper dimensions set?
        
        for i,j in zip(dims, found):
            print("{:<70} : {}".format(i,j))
    # FIN class CDMDataset

"""
Speed test of h5py fancy indexing:

tmp = CDMDataset('01001')

%%time
    ...: idx = np.where(tmp.observations_table.observed_variable[()]==85)
    ...: tmp.observations_table.observation_value[()][idx], tmp.observations_table.date_time[()][idx], tmp.observations_table.z_coordinate[()][idx]
    ...:
    ...:
CPU times: user 18.8 s, sys: 5.17 s, total: 24 s
Wall time: 24 s


# insane indexing time

%%time
    ...: idx = np.where(tmp.observations_table.observed_variable[()]==85)
    ...: tmp.observations_table.observation_value[idx], tmp.observations_table.date_time[idx], tmp.observations_table.z_coordinate[idx]
    ...:
    ...:

"""
