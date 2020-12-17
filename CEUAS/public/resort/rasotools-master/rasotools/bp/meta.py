# -*- coding: utf-8 -*-

__all__ = ['location_change', 'sondetype', 'metadata']


def location_change(lon, lat, dim='time', ilon=None, ilat=None, as_event=True, distance_threshold=10, window=180,
                    dates=None, **kwargs):
    """ Convert location change to breakpoint series

    Args:
        lon (DataArray): longitudes
        lat (DataArray): latitudes
        dim (str): datetime dimension
        ilon (float): location longitude
        ilat (float): location latitude
        as_event (bool): calculate location change as event (distance_threshold)
        distance_threshold (float): distance threshold for events
        window (int): event influence in days
        **kwargs:

    Returns:
        DataArray : distances between locations

    Notes:
        distance_threshold : degree 0.1 == 11.1 km, 0.01 == 1.1 km changes

    """
    import numpy as np
    from xarray import DataArray, full_like
    from ..fun.cal import distance

    # count occurence of each coordinate pair
    # use the most common (also the most recent?)
    # to estimate distance from,
    # era-interim has a distance of about 80km so only if larger it would make sense to split?
    if not isinstance(lon, DataArray):
        raise ValueError('requires a DataArray', type(lon))
    if not isinstance(lat, DataArray):
        raise ValueError('requires a DataArray', type(lat))

    lon = lon.copy()
    lat = lat.copy()
    lon = lon.bfill(dim)
    lat = lat.bfill(dim)

    dist = full_like(lon, 0, dtype=np.float)
    dist.name = 'distance'
    dist.attrs['units'] = 'km'
    fdistance = np.vectorize(distance)
    ishape = lon.values.shape
    lon = lon.values.flatten()
    lat = lat.values.flatten()
    if ilon is None and ilat is None:
        if lon.size > 1:
            # distance between more recent and less recent
            tmp = fdistance(lon[1:], lat[1:], lon[:-1], lat[:-1])
            tmp = np.append(tmp, tmp[-1])
        else:
            tmp = np.array([0])
        dist.values = tmp.reshape(ishape)
        dist.attrs['method'] = 'Backwards'
    else:
        tmp = fdistance(lon, lat, ilon, ilat)
        dist.values = tmp.reshape(ishape)
        dist.attrs['method'] = 'Point(%f E, %f N)' % (ilon, ilat)

    # Check for duplicates
    if dist[dim].to_index().duplicated().any():
        dist = dist.isel({dim: ~dist[dim].to_index().duplicated()})

    if as_event:
        dist.values = (dist > distance_threshold).astype(float) \
            .rolling(center=True, min_periods=1, **{dim: window}).mean() \
            .rolling(center=True, min_periods=1, **{dim: window}).sum()
        dist.attrs['threshold'] = distance_threshold
        dist.attrs['standard_name'] = 'location_change_point'
        if dates is not None:
            dist = dist.reindex(time=dates, method='nearest', fill_value=0)  # reindex
    else:
        dist.attrs['standard_name'] = 'distance'
        if dates is not None:
            dist = dist.reindex(time=dates, method='nearest').bfill(dim).ffill(dim)

    return dist


def sondetype(data, dim='time', window=30, as_event=True, dates=None, **kwargs):
    """ Make breakpoint series from sondetype changes

    Args:
        data (DataArray): Sondetype values
        dim (str): datetime dimension
        window (int): smoothing window
        **kwargs:

    Returns:
        DataArray : sondetype changes
    """
    import numpy as np
    from xarray import DataArray

    if not isinstance(data, DataArray):
        raise ValueError('requires a DataArray', type(data))

    data = data.copy()
    # Check for duplicates
    if data[dim].to_index().duplicated().any():
        data = data.isel({dim: ~data[dim].to_index().duplicated()})

    data = data.bfill(dim).rolling(center=True, min_periods=1, **{dim: window}).median().bfill(dim).ffill(dim)
    if as_event:
        idx = [slice(None)] * len(data.dims)
        axis = data.dims.index(dim)
        n = data.coords[dim].size
        idx[axis] = slice(0, n - 1)
        data.values[tuple(idx)] = (np.diff(data.values) != 0)
        idx[axis] = n - 1
        data.values[tuple(idx)] = 0
        data = data.rolling(center=True, **{dim: window}).mean().rolling(center=True, **{dim: window}).sum()
        data /= 2
        data = data.fillna(0)
        if dates is not None:
            data = data.reindex(time=dates, method='nearest', fill_value=0)  # reindex
    else:
        if dates is not None:
            data = data.reindex(time=dates, method='nearest').bfill(dim).ffill(dim)  # fill missing values
    return data


def metadata(ident, dates, dim='time', window=30, **kwargs):
    import numpy as np
    import pandas as pd
    from ..fun import message, get_data

    if not isinstance(ident, (str, int)):
        raise ValueError("Requires a string or int Radiosonde ID")

    if not isinstance(dates, pd.DatetimeIndex):
        raise ValueError("Requires a pandas.DatetimeIndex")

    data = pd.DataFrame(index=dates, columns=['day', 'igra_message', 'igra_event', 'wmo_sondetype', 'wmo_event'])
    data['day'] = data.index.to_period('D')
    data['igra_event'] = 0
    data['wmo_event'] = 0

    igra_id = False
    if isinstance(ident, str):
        if len(ident) > 6:
            igra_id = True
    else:
        ident = "%06d" % ident

    igra = pd.read_json(get_data('igrav2_metadata.json'))

    if igra_id:
        igra = igra[igra.igraid == ident]
    else:
        igra = igra[igra.wmoid == int(ident)]

    message('[META] IGRA Events:', igra.shape, **kwargs)
    if kwargs.get('verbose', 0) > 1:
        print(igra)

    igra = igra.drop_duplicates('date').set_index('date')
    igra.index = igra.index.tz_localize(tz=None)
    igra['message'] = igra['befinfo'] + ' ' + igra['link'] + ' ' + igra['aftinfo'] + ' ' + igra['comment']
    igra['message'] = igra['message'].str.strip()
    # Find closest match to given dates
    for idate in (igra.index.tz_localize(tz=None)):
        logic = abs(data.index - idate) <= pd.Timedelta(30, unit='d')
        jdate = np.argmin(abs(data.index - idate))
        message("[META]", idate, data.index[jdate], logic.sum(), data['day'].iloc[jdate], **kwargs)
        idx = np.where(data['day'] == data['day'].iloc[jdate])[0]
        data.iloc[idx, 2] = 1
        data.iloc[idx, 1] = igra.loc[idate, 'message']

    wmo = pd.read_json(get_data('wmo_metadata.json'))
    wmo = wmo[wmo.id.str.contains(ident)]
    message('[META] WMO Type changes', wmo.shape, **kwargs)
    data.index = data.index.tz_localize(tz='UTC')
    for row in wmo.iterrows():
        data.loc[slice(row[1]['start'], row[1]['stop']), 'wmo_sondetype'] = row[1]['c']

    data = data.drop('day', axis=1)
    data.index = data.index.tz_localize(tz=None)

    data['wmo_sondetype'] = data['wmo_sondetype'].replace(0, np.nan)
    data['wmo_sondetype'] = data['wmo_sondetype'].replace(-1, np.nan).astype(float)
    data.index.name = dim
    data = data.to_xarray()
    data['wmo_event'] = sondetype(data['wmo_sondetype'], dim=dim, window=window, dates=dates)
    data['igra_event'] = data['igra_event']\
        .rolling(center=True, min_periods=1, **{dim: window}).mean()\
        .rolling(center=True, min_periods=1, **{dim: window}).sum().fillna(0)
    return data


def wmo_code_table():
    """ Metadata from WMO Common Code Tables to binary and alphanumeric codes
    Table 3685
    Version 7.11.2018

    Returns:
        DataFrame : Radiosonde Code Table C2
    """
    import pandas as pd
    from ..fun import get_data

    return pd.read_csv(get_data('Common_C02_20181107_en.txt'))


def open_igra_metadata(filename):
    """ Read IGRAv2 _metadata file according to readme

    igra2-_metadata-readme.txt

    Documentation for IGRA Station History Information
    Accompanying IGRA Version 2.0.0b1
    August 2014

    Args:
        filename (str):  igra2-_metadata.txt

    Returns:
        DataFrame
    """
    import pandas as pd
    infos = """
    IGRAID         1- 11   Character
    WMOID         13- 17   Integer
    NAME          19- 48   Character
    NAMFLAG       50- 50   Character
    LATITUDE      52- 60   Real
    LATFLAG       62- 62   Character
    LONGITUDE     64- 72   Real
    LONFLAG       74- 74   Character
    ELEVATION     76- 81   Real
    ELVFLAG       83- 83   Character
    YEAR          85- 88   Integer
    MONTH         90- 91   Integer
    DAY           93- 94   Integer
    HOUR          96- 97   Integer
    DATEIND       99- 99   Integer
    EVENT        101-119   Character
    ALTIND       121-122   Character
    BEFINFO      124-163   Character
    BEFFLAG      164-164   Character
    LINK         166-167   Character
    AFTINFO      169-208   Character
    AFTFLAG      209-209   Character
    REFERENCE    211-235   Character
    COMMENT      236-315   Character
    UPDCOM       316-346   Character
    UPDDATE      348-354   Character
    """
    import numpy as np
    colspecs = []
    header = []
    types = {}
    for iline in infos.splitlines():
        if iline == '':
            continue
        ih = iline[0:11].strip().lower()
        header.append(ih)
        ii = int(iline[13:16]) - 1
        ij = int(iline[17:20])
        colspecs.append((ii, ij))
        it = iline[22:].strip()
        if it == 'Character':
            it = 'str'
        elif it == 'Real':
            it = 'float'
        else:
            it = 'int'
        types[ih] = it

    data = pd.read_fwf(filename, colspecs=colspecs, header=None, dtype=types, names=header)
    data = data.replace('nan', '')
    data['date'] = pd.to_datetime((data.year * 1000000 +
                                   np.where(data.month.values == 99, 6, data.month.values) * 10000 +
                                   np.where(data.day.values == 99, 15, data.day.values) * 100 +
                                   np.where(data.hour.values == 99, 0, data.hour.values)).apply(str), format='%Y%m%d%H')
    return data
