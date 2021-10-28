#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__all__ = ['download_igrav2_stationlist', 'read_igrav2_stationlist', 'download_wmo_stationlist', 'read_wmo_statiolist',
           'build_stationlist_from_netcdf_archive']


def build_stationlist_from_netcdf_archive(pattern, directory=None, **kwargs):
    """

    Args:
        pattern (str): file pattern (ERA5_*.nc)
        directory: archive directory (rasotools.config.rasodir)
        **kwargs:

    Returns:
        DataFrame : station information (id, lon, lat, alt, start, end, toal, name)
    """
    import pandas as pd
    from .. import config
    from . import find_files, mp, levelup, message
    from .netcdf import view

    if directory is None:
        directory = config.rasodir

    message("Searching in", directory, " with ", pattern, **kwargs)
    files = find_files(directory, pattern)
    nn = len(files)
    if nn == 0:
        raise RuntimeError("No files found: ", directory, pattern)

    message("#" * 80, **kwargs)
    message("Example:", **kwargs)
    view(files[0])
    message("Station Information:", **kwargs)
    message(pd.Series(_get_netcdf_infos((0, files[0]), **kwargs)), **kwargs)
    message("#" * 80, **kwargs)
    # stations = []
    message("Starting Processing Pool(10) for ", nn, " files", **kwargs)
    stuff = mp.make_process_list(zip(range(nn), files), _get_netcdf_infos, **levelup(**kwargs))
    errors, stations = mp.execute_process(stuff,
                                          npp=kwargs.get('npp', 10),
                                          return_value=True,
                                          loop='loop' in kwargs.keys(),
                                          extend=kwargs.get('extend', False))

    stations = pd.DataFrame(stations)
    stations['src'] = list(map(lambda x: x.replace('.nc', '').split('_')[-2], stations['file']))
    stations = stations[~stations['src'].str.contains('Kopie')]  # remove something wired
    stations['id'] = stations.id.str.replace('-dataset.txt.gz', '')
    for i in range(stations.id.size):
        try:
            if len(stations.loc[i, 'id']) > 6:
                stations.loc[i, 'id'] = "%06d" % int(stations.loc[i, 'id'][-6:])
            else:
                stations.loc[i, 'id'] = "%06d" % int(stations.loc[i, 'id'])
        except:
            pass
    stations = stations.set_index(['id', 'src']).sort_index()
    return stations


def _get_netcdf_infos(arg, last_pos=True, debug=False,  **kwargs):
    import numpy as np
    import pandas as pd
    import xarray as xr
    #
    # recover arguments
    #
    j, ifile = arg
    #
    # Identifier
    #
    ident = ifile.split('/')

    if len(ident) >= 2:
        ident = ident[-2]
    else:
        ident = "SON%09d" % j
        print("unclear ident:", ifile, ident)

    infos = {
        'id': None,
        'lon': None,
        'lat': None,
        'alt': None,
        'start': None,
        'end': None,
        'total': None,
        'name': None
    }
    if 'gz' in ifile:
        import gzip
        g = gzip.open(ifile)
    else:
        g = ifile

    with xr.open_dataset(g) as ds:
        #
        # Attributes
        #
        for iatt, ival in ds.attrs.items():
            if iatt in ['id', 'ident', 'station_id', 'station_ident']:
                infos['id'] = ival
            if iatt in ['lon', 'longitude', 'stlon', 'station_lon']:
                infos['lon'] = ival
            if iatt in ['lat', 'latitude', 'stlat', 'station_lat']:
                infos['lat'] = ival
            if iatt in ['alt', 'stalt', 'station_alt']:
                infos['alt'] = ival
            if iatt in ['name', 'station_name', 'stname']:
                infos['name'] = ival
        #
        # Some information might be Variables (IGRA, UADB)
        #
        # if infos['lon'] is None:
        for i in ['lon', 'longitude', 'stlon']:
            if i in list(ds.data_vars):
                if last_pos:
                    infos['lon'] = float(ds[i].values[-1])
                else:
                    infos['lon'] = ds[i].copy()
        # if infos['lat'] is None:
        for i in ['lat', 'latitude', 'stlat']:
            if i in list(ds.data_vars):
                if last_pos:
                    infos['lat'] = float(ds[i].values[-1])
                else:
                    infos['lat'] = ds[i].copy()
        # if infos['alt'] is None:
        for i in ['alt', 'stalt']:
            if i in list(ds.data_vars):
                infos['alt'] = float(ds[i].values[-1])
        #
        # Min Year, Max Year
        #
        idate = None
        for i in ['time', 'date']:
            if i in list(ds.coords):
                infos['start'] = int(np.min(ds[i].dt.year))
                infos['end'] = int(np.max(ds[i].dt.year))
                idate = i
        #
        # Count events
        #
        if idate is None:
            infos['total'] = int(ds.count().to_array().max())
        else:
            infos['total'] = int(ds.count(idate).to_array().max())
    #
    # fix filenames
    #
    infos['file'] = ifile.split('/')[-1]
    #
    # fix ident
    #
    if infos['id'] is None:
        infos['id'] = ident
    #
    # Convert string to float
    #
    if isinstance(infos['lon'], str):
        try:
            infos['lon'] = float(infos['lon'].split()[0])
        except:
            pass

    if isinstance(infos['lat'], str):
        try:
            infos['lat'] = float(infos['lat'].split()[0])
        except:
            pass
    if isinstance(infos['alt'], str):
        try:
            infos['alt'] = float(infos['alt'].split()[0])
        except:
            pass

    try:
        tmp = pd.concat([infos['lon'].to_dataframe(), infos['lat'].to_dataframe()], axis=1)
        # infos['nloc'] = (tmp.groupby(['lon', 'lat']).size()).size
        new = []
        for i, m in list(tmp.groupby(['lon', 'lat'])):
            j = infos.copy()
            j['lon'] = m['lon'].iloc[-1]
            j['lat'] = m['lat'].iloc[-1]
            j['ptime'] = m.index.year[-1]
            j['pcount'] = m.size
            new.append(j)
        infos = new
        # infos['lon'], infos['lat'] = tmp.groupby(['lon', 'lat']).size().idxmax()
    except Exception as e:
        if debug:
            raise e
        infos['lon'], infos['lat'] = infos['lon'].item(-1),infos['lat'].item(-1)

    return infos


def read_igrav2_stationlist(stationfile, **kwargs):
    """Read IGRA Radiosondelist

    or download

    Parameters
    ----------
    new         bool
    ifile    str
    verbose     int

    Returns
    -------
    DataFrame
    """
    import os
    import numpy as np
    import pandas as pd
    from . import message

    if stationfile is None and os.path.isfile('results/igra2-station-list.csv'):
        try:
            return pd.read_csv('results/igra2-station-list.csv')
        except:
            pass
    if stationfile is None:
        try:
            from . import get_data
            stationfile = get_data("igra2-station-list.txt")
        except ImportError:
            message("Module Data not found: ax", **kwargs)
            pass

    message("Reading :", stationfile, **kwargs)
    try:
        infile = open(stationfile)
        tmp = infile.read()
        data = tmp.splitlines()[0:1]

    except IOError as e:
        message("File not found: " + stationfile, **kwargs)
        raise e
    else:
        infile.close()

    out = pd.DataFrame(columns=['id', 'wmo', 'lat', 'lon', 'alt', 'state', 'name', 'start', 'end', 'total'])

    for i, line in enumerate(data):
        id = line[0:11]

        try:
            id2 = "%06d" % int(line[5:11])  # substring

        except:
            id2 = ""

        lat = float(line[12:20])
        lon = float(line[21:30])
        alt = float(line[31:37])
        state = line[38:40]
        name = line[41:71]
        start = int(line[72:76])
        end = int(line[77:81])
        count = int(line[82:88])
        out.loc[i] = (id, id2, lat, lon, alt, state, name, start, end, count)

    out.loc[out.lon <= -998.8, 'lon'] = np.nan  # repalce missing values
    out.loc[out.alt <= -998.8, 'alt'] = np.nan  # repalce missing values
    out.loc[out.lat <= -98.8, 'lat'] = np.nan  # replace missing values
    out['start'] = out['start'].astype(int)
    out['end'] = out['end'].astype(int)
    out['total'] = out['total'].astype(int)
    out['name'] = out.name.str.strip()
    out = out.set_index('id')
    return out


def download_igrav2_stationlist():
    import urllib.request
    url = 'ftp://ftp.ncdc.noaa.gov/pub/data/igra/igra2-station-list.txt'

    try:
        urllib.request.urlretrieve(url, './results/igra2-station-list.txt')
    except Exception as e:
        print("Error: ", repr(e))
    print("Download completed")


def download_wmo_stationlist():
    import urllib.request
    url = 'https://oscar.wmo.int/oscar/vola/vola_legacy_report.txt'

    try:
        urllib.request.urlretrieve(url, './results/wmo-stations.txt')
    except Exception as e:
        print("Error: ", repr(e))
    print("Download completed: ./results/wmo-stations.txt")


def read_wmo_statiolist(stationfile, minimal=True, only_raso=True):
    """ Read WMO Radiosonde Station List

    ANTON(T)    : Antarctic Observing Network upper-air station (TEMP)
    GUAN        : GCOS Upper-Air Network station
    RBSN(T)     : Regional Basic Synoptic Network upper-air station (TEMP)
    RBSN(P)     : Regional Basic Synoptic Network upper-air station (PILOT)
    RBSN(ST)    : Regional Basic Synoptic Network surface and upper-air station (SYNOP/TEMP)
    RBSN(SP)    : Regional Basic Synoptic Network surface and upper-air station (SYNOP/PILOT)
    WN          : Upper-wind observations made by using navigation aids (NAVAID)
    WR          : Upper-wind observations made by radar
    WT          : Upper-wind observations made by radiotheodolite
    WTR         : Upper-wind observations made by radiotheodolite/radar composite method

    Args:
        stationfile (str): filename to read (wmo-stations.txt)
        minimal (bool): subset of columns
        only_raso (bool): only radiosonde stations

    Returns:
        pd.DataFrame : station list
    """
    import numpy as np
    import pandas as pd

    if stationfile is None:
        try:
            from .help import get_data
            stationfile = get_data('wmo-stations.txt')
        except ImportError:
            pass
    print("Reading :", stationfile)
    try:
        wd = pd.read_csv(stationfile, sep='\t')
    except IOError as e:
        print("Error missing file: ", stationfile)
        raise e

    sign = np.where(wd.Longitude.apply(lambda x: 'E' in x).values, 1, -1)
    wd.loc[:, 'Longitude'] = wd.Longitude.apply(
        lambda x: np.sum(np.float_(x[:-1].split()) / [1., 60., 3600.])).values * sign

    sign = np.where(wd.Latitude.apply(lambda x: 'S' in x).values, -1, 1)
    wd.loc[:, 'Latitude'] = wd.Latitude.apply(
        lambda x: np.sum(np.float_(x[:-1].split()) / [1., 60., 3600.])).values * sign

    # wd.ix[:, 'CountryArea'] = wd.CountryArea.apply(lambda x: x.split('/')[0])
    # wd.ix[:, 'RegionName'] = wd.RegionName.apply(lambda x: x.split('/')[0])
    wd = wd.rename(columns={'Longitude': 'lon', 'Latitude': 'lat', 'IndexNbr': 'id', 'StationName': 'name',
                            'Hp': 'alt', 'CountryArea': 'area', 'RegionName': 'region', 'StationId': 'wigos'})

    # require variables named: id, lon,lat,alt,name,count
    wd['id'] = np.where(np.isfinite(wd.id), wd.id.map('{:06.0f}'.format), '')  # wd.id.map('{:06.0f}'.format)
    rasotypes = ['RBSN(T)', 'RBSN(P)', 'RBSN(ST)', 'RBSN(SP)', 'GUAN', 'ANTON(T)', 'R']
    status = pd.concat([_wmo_obstype(wd.ObsRems, ityp) for ityp in rasotypes], axis=1, keys=rasotypes)

    print(status.sum())

    wd['raso'] = np.any([_wmo_obstype(wd.ObsRems, ityp) for ityp in rasotypes], axis=0)
    if only_raso:
        wd = wd[wd.raso]

    if minimal:
        wd = wd.set_index('id')
        wd = wd[['lon', 'lat', 'name', 'alt', 'area', 'region', 'wigos']].sort_index().drop_duplicates()

    return wd


def _wmo_obstype(data, typ):
    import pandas as pd
    out = []
    for i in data.values:
        if ',' in i:
            j = i.split(',')
            status = False
            for k in j:
                if k.strip() == typ:
                    status = True
            out += [status]
        else:
            if i.strip() == typ:
                out += [True]
            else:
                out += [False]
    return pd.Series(out, index=data.index)


if __name__ == "__main__":
    import os
    import sys
    import datetime

    doigra = False
    dowmo = False


    def usage():
        print("""
        Read radiosonde station list from IGRA v2 or WMO to pandas Dataframe
        Call: station.py -i -w [FILE]
        Options:
            -i      read IGRAv2
            -w      read WMO
            [FILE]  file of wmo or igrav2 list, alternative: download
        """)


    __name__ = os.getcwd().split('/')[-1]  # fixes realtive imports

    if len(sys.argv) < 2:
        usage()
        sys.exit(0)

    filename = "download"
    if len(sys.argv) == 3:
        if os.path.isfile(sys.argv[2]):
            filename = sys.argv[2]

    if sys.argv[1] not in ("-i", "-w"):
        usage()
        sys.exit(1)
    elif sys.argv[1] == "-i":
        doigra = True
    else:
        dowmo = True

    if filename == "download":
        if not os.path.isdir('./results/'):
            os.makedirs('./results/')

        if doigra:
            download_igrav2_stationlist()
            filename = './results/igra2-station-list.txt'
        if dowmo:
            download_wmo_stationlist()
            filename = './results/wmo-stations.txt'

    if not os.path.isfile(filename):
        raise IOError('File not found: %s' % filename)

    try:
        if doigra:
            out = read_igrav2_stationlist(filename)
        else:
            out = read_wmo_statiolist(filename)
        print(out)

    except Exception as e:
        print("Error:", repr(e))

    with open('job.execution', 'a') as f:
        f.write(datetime.datetime.now().isoformat(timespec='auto') + " ".join(sys.argv) + "\n")
