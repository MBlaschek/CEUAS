# -*- coding: utf-8 -*-

__all__ = ['stationlist', 'datafiles']

"""
Examples:
>>> stations = rt.fun.station.read_igrav2_stationlist(None)
>>> stations = stations[stations.wmo != '']   # only WMO
>>> net = rasotools.network_from_stationlist('GLOBAL', stations)

"""


def stationlist(name, stations, ident='wmo', lon='lon', lat='lat', region=None, query=None, **kwargs):
    """ Setup a Radiosonde network based on a station list with coordinate informations

    Args:
        name (str):  Name
        stations (pd.DataFrame, dict): Station information DataFrame
        ident (str): column, key in stationlist
        lon (str): column, key in stationlist
        lat (str): column, key in stationlist
        region (tuple): (lon1, lon2, lat1, lat2)
        query (str): pandas eval on stations

    Returns:
        Network : Radiosonde network class
    """
    import numpy as np
    import pandas as pd
    from ..cls import Network
    from ..fun import message

    if not isinstance(stations, (dict, pd.DataFrame)):
        raise ValueError("requires a dict or DataFrame:", type(stations))

    new = Network(name)

    if not isinstance(stations, pd.DataFrame):
        stations = pd.DataFrame(stations)

    if query is not None:
        message("Stations: ", len(stations), **kwargs)
        stations = stations.query(query)
        message("Query: ", len(stations), **kwargs)

    nn = len(stations)
    stations = stations.dropna(subset=(ident, lon, lat), axis=0)
    if stations[ident].dtype != 'object':
        stations[ident] = stations[ident].map('{:06.0f}'.format)
    message("Dropping missing: ", nn - len(stations), **kwargs)
    # new.idents = list(map('{:06.0f}'.format, np.asarray(stations[ident])))
    new.idents = list(np.asarray(stations[ident]))
    new.lon = np.asarray(stations[lon])
    if any(new.lon > 180):
        message("Converting Longitudes to -180 180", **kwargs)
        new.lon = np.where(new.lon > 180, new.lon - 360., new.lon)

    new.lat = np.asarray(stations[lat])

    if region is not None:
        idx = (new.lon >= region[0]) & (new.lon <= region[1]) & (new.lat >= region[2]) & (new.lat <= region[3])
        if np.sum(idx) == 0:
            raise RuntimeError("No Radiosonde inside region:", region)
        message("Region: ", region, idx.size, np.sum(idx), **kwargs)
        new.idents = list(np.asarray(new.idents)[idx])
        new.lon = new.lon[idx]
        new.lat = new.lat[idx]
        stations = stations.iloc[idx]

    new.stations = stations.set_index(ident)
    return new


def datafiles(name, pattern, directory=None, **kwargs):
    from ..fun.station import build_stationlist_from_netcdf_archive
    #
    # find files
    # read meta information from files
    #
    stationlist = build_stationlist_from_netcdf_archive(pattern, directory=directory, **kwargs)
    new = stationlist(name, stationlist)

    return new
