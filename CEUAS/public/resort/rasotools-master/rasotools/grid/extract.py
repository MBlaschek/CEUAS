# -*- coding: utf-8 -*-

__all__ = ['extract_locations']


def extract_locations(data, lon, lat, method='bilinear', raw=False, dim=None, debug=False, **kwargs):
    """ Extract location(s) from a DataArray

    Args:
        data : xr.DataArray
        lon : float, int, list
        lat : float, int, list
        method : str
            Interpolation Method: point, bilinear, distance
        raw : bool
        dim : str, list, Index
            concat along this new dimension
        debug : bool

    Returns:
        list or DataArray
    """
    import numpy as np
    from xarray import DataArray, concat
    from ..fun import message

    if not isinstance(data, DataArray):
        raise ValueError()

    if isinstance(lon, (int, float)):
        lon = [float(lon)]

    if isinstance(lat, (int, float)):
        lat = [float(lat)]

    if method not in ['point', 'bilinear', 'distance']:
        raise ValueError("Method unknown: point, bilinear, distance")

    data = data.copy()
    locations = []
    lon = np.array(lon)
    lat = np.array(lat)
    name_lon, name_lat = get_lonlat(list(data.dims))  # list
    if name_lat is None or name_lon is None:
        print(data.dims)
        raise ValueError("DataArray does not have named lon, lat coordinates")

    # need 0 to 360 deg
    xlons = data[name_lon].values
    if (xlons < 0).any():
        data[name_lon].values = np.where(xlons < 0, xlons + 360., xlons)
        message("Adjusting GRID Longitudes ...", **kwargs)

    # Input needs to be on the same lon
    if (lon < 0).any():
        lon = np.where(lon < 0, lon + 360., lon)
        message("Adjusting INPUT Longitudes ...", **kwargs)

    order = data.dims  # tuple
    ilon = order.index(name_lon)
    ilat = order.index(name_lat)
    lons = data[name_lon].values[:]  # copy
    lats = data[name_lat].values[:]  # copy
    # get both axis
    iaxis = sorted([ilon, ilat])

    iorder = [i for i in order if i not in [name_lat, name_lon]]
    dattrs = dict(data.attrs)
    newcoords = {i: data[i].copy() for i in iorder}

    if dim is None:
        dim = 'station'
    if len(iorder) == 0:
        iorder = [dim]

    mm = 0
    for jlon, jlat in zip(lon, lat):
        try:
            # indices and weights for location
            indices = [slice(None)] * len(order)
            #
            #  Point or Bilinear interpolation from rectangle
            #
            if method == 'point':
                ix, iy = get_pos(jlon, jlat, lons, lats)
                weights = np.array([1.])

            elif method == 'distance':
                ix, iy, weights = distance_weights(lons, lats, jlon, jlat)

            else:
                idx, points = rectangle(jlon, jlat, lons, lats)
                ix, iy, weights = bilinear_weights(jlon, jlat, idx, lons, lats)
                # Example for ERA-Interim GRID
                # plot_weights(ix, iy, weights, dataset.values[0, 0, :, :], lons, lats, jlon, jlat,
                #               lonlat=ilon < ilat, dx=-np.diff(lons).mean() / 2., dy=np.diff(lats).mean() / 2.)
            indices[ilon] = ix
            indices[ilat] = iy
            # convert to tuple ?
            indices = tuple(indices)
            #
            # get values
            #
            tmp = data.values[indices]
            if len(order) > 2:
                w = np.empty_like(tmp)
                w[::] = weights  # fill in
                weights = w

            if method != 'point':
                # Adjust Weights to shape and missing values
                weights = np.where(np.isnan(tmp), np.nan, weights)
                ws = np.nansum(weights, axis=iaxis[0], keepdims=True)
                w = np.where(ws == 0, np.nan, weights / ws)  # nansum changed !!
                weights = w

                tmp = np.nansum(tmp * weights, axis=iaxis[0])  # only the lower axis
            else:
                tmp = tmp * weights

            # results
            if not raw:
                dattrs['cell_method'] = "%s,%s: intp(%s)" % (name_lon, name_lat, method)
                newcoords[name_lon] = jlon - 360. if jlon > 180 else jlon
                newcoords[name_lat] = jlat
                newcoords[dim] = mm
                tmp = DataArray(tmp, coords=newcoords, dims=iorder, name=data.name, attrs=dattrs)
            locations.append(tmp.copy())
            mm += 1
        except Exception as e:
            if debug:
                raise e

            print(e)
            locations.append([])

    if len(locations) == 1:
        return locations[0]

    if dim is not None:
        locations = concat(locations, dim=dim)

    elif len(iorder) == 0:
        locations = concat(locations, dim='station')

    else:
        pass
    return locations


def get_lonlat(dims):
    "Look for lon lat -> return keys"
    jlon = None
    jlat = None
    for i in dims:
        if i.lower() in ['longitude', 'long', 'lon']:
            jlon = i

        if i.lower() in ['latitude', 'lati', 'lat']:
            jlat = i

    return jlon, jlat


def get_pos(ilon, ilat, lon, lat):
    import numpy as np
    # should be also valid for points in the range of -180 to 180
    if np.argmin([np.min(np.abs(lon - ilon)), np.min(np.abs(lon + 360. - ilon))]) == 0:
        return np.argmin(np.abs(lon - ilon)), np.argmin(np.abs(lat - ilat))
    else:
        return np.argmin(np.abs(lon + 360. - ilon)), np.argmin(np.abs(lat - ilat))


def rectangle(ilon, ilat, lons, lats, n=2):
    import numpy as np
    import itertools
    ix, iy = get_pos(ilon, ilat, lons, lats)
    if ix < 2 or ix > lons.size - 2:
        # wrap around for lon
        ilon = ilon if ilon <= 180. else ilon - 360.
        jx = np.abs(np.where(lons > 180., lons - 360., lons) - ilon).argsort()[slice(None, n)]
    else:
        jx = np.abs(lons - ilon).argsort()[slice(None, n)]
    jy = np.abs(lats - ilat).argsort()[slice(None, n)]
    jy = np.unique(np.where(jy > lats.size, lats.size, jy))  # upper limit
    return list(itertools.product(jx, jy)), list(itertools.product(lons[jx], lats[jy]))


def bilinear_weights(x, y, p, lons, lats):
    """Interpolate (x,y) from values associated with four points.

    The four points are a list of four triplets:  (x, y, value).
    The four points can be in any order.  They should form a rectangle.
    """
    # See formula at:  http://en.wikipedia.org/wiki/Bilinear_interpolation

    import numpy as np
    p = sorted(p)  # order points by x, then by y
    #
    (x1, y1), (_x1, y2), (x2, _y1), (_x2, _y2) = p
    # (lon1,lat1), (lon2,lat1), (lon1,lat2), (lon2,lat2)
    w = np.array([(lons[x2] - x) * (lats[y2] - y),
                  (x - lons[x1]) * (lats[y2] - y),
                  (lons[x2] - x) * (y - lats[y1]),
                  (x - lons[x1]) * (y - lats[y1])]) / ((lons[x2] - lons[x1]) * (lats[y2] - lats[y1]) + 0.0)
    if (w < 0).any():
        lons = np.where(lons > 180., lons - 360., lons)
        x = x if x <= 180. else x - 360.
        w = np.array([(lons[x2] - x) * (lats[y2] - y),
                      (x - lons[x1]) * (lats[y2] - y),
                      (lons[x2] - x) * (y - lats[y1]),
                      (x - lons[x1]) * (y - lats[y1])]) / ((lons[x2] - lons[x1]) * (lats[y2] - lats[y1]) + 0.0)

    # Make sure it is 1
    return np.array([x1, x2, x1, x2]), np.array([y1, y1, y2, y2]), w / w.sum()


def distance_weights(lons, lats, ilon, ilat):
    import numpy as np
    ix, iy = get_pos(ilon, ilat, lons, lats)
    ny = len(lats)
    nx = len(lons)
    dx = np.abs(lons - ilon)
    dy = np.abs(lats - ilat)
    if ix + 1 < nx:
        if dx[ix - 1] < dx[ix + 1]:
            idx = [ix - 1, ix]
        else:
            idx = [ix, ix + 1]
    else:
        if dx[ix - 1] < dx[0]:
            idx = [ix - 1, ix]
        else:
            idx = [ix, 0]

    if iy + 1 < ny:
        if dy[iy - 1] < dy[iy + 1]:
            idy = [iy - 1, iy]
            if iy - 1 < 0:
                idy = [iy, iy + 1]
        else:
            idy = [iy, iy + 1]
            if iy + 1 >= ny:
                idy = [iy - 1, iy]
    else:
        idy = [iy - 1, iy]

    dist = []
    for i in idx:
        for j in idy:
            dist.append(1 / distance(ilon, ilat, lons[i], lats[j]))
    dist = np.array(dist)
    return [idx[0], idx[0], idx[1], idx[1]], [idy[0], idy[1], idy[0], idy[1]], (dist / np.sum(dist))


def distance(lon, lat, lon0, lat0):
    """
    Calculates the distance between a point and an array of points

    Parameters
    ----------
    lon     Longitude Vector
    lat     Latitude Vector
    lon0    Longitude of Position
    lat0    Latitude of Position

    Returns
    -------
    numpy.array
    """

    import numpy as np
    lonvar, latvar = np.meshgrid(lon, lat)
    rad_factor = np.pi / 180.0  # for trignometry, need angles in radians
    # Read latitude and longitude from file into numpy arrays
    latvals = latvar * rad_factor
    lonvals = lonvar * rad_factor
    # ny,nx = latvals.shape
    lat0_rad = lat0 * rad_factor
    lon0_rad = lon0 * rad_factor
    # Compute numpy arrays for all values, no loops
    clat, clon = np.cos(latvals), np.cos(lonvals)
    slat, slon = np.sin(latvals), np.sin(lonvals)
    delX = np.cos(lat0_rad) * np.cos(lon0_rad) - clat * clon
    delY = np.cos(lat0_rad) * np.sin(lon0_rad) - clat * slon
    delZ = np.sin(lat0_rad) - slat;
    dist_sq = delX ** 2 + delY ** 2 + delZ ** 2  # Distance
    return np.squeeze(dist_sq)


#
#  DEBUGGING
#

def plot_weights(idx, idy, weights, data, lons, lats, jlon, jlat, lonlat=True, dx=0, dy=0):
    import matplotlib
    import matplotlib.pyplot as plt

    if not lonlat:
        data = data.T

    ix, iy = idx[0], idy[0]
    tmp = data[ix - 2:ix + 4, iy - 2:iy + 4]

    norm = matplotlib.colors.Normalize(vmin=tmp.min(), vmax=tmp.max())
    plt.figure()
    plt.pcolormesh(lons[ix - 2:ix + 4] + dx, lats[iy - 2:iy + 4] + dy, tmp.T, norm=norm)

    for i, j, k in zip(idx, idy, weights):
        plt.scatter(lons[i], lats[j], c=data[i, j], norm=norm, edgecolors='k')
        plt.text(lons[i], lats[j], "%6.2f\n(%6.2f %6.2f)\n%6.5f" % (data[i, j], lons[i], lats[j], k),
                 horizontalalignment='right')

    plt.plot(jlon, jlat, 'o', c='k')
    plt.colorbar()
