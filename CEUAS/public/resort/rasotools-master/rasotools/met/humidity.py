# -*- coding: utf-8 -*-
import numpy as np
from .esat import svp
from ..fun.constants import eps
#
# common problem with missing values
#
np.seterr(invalid='ignore')


__all__ = ['sh2ppmv', 'sh2vap', 'vap2sh', 'ppmv2pa', 'rh2vap', 'rh2ppmv', 'dewpoint', 'frostpoint']


def vap2ppmv(e, p):
    """ Convert water vapor pressure into PPMV

    Args:
        e (float): water vapor pressure
        p (float): air pressure

    Returns:
        e (float): PPMV
    """
    return 1e6 * (e / (p - e))


def sh2vap(q, p):
    """Specific Humidity to Water vapor pressure
    formula derived from ratio of humid air vs. total air

    Parameters
    ----------
    q       specific humidity [kg/kg]
    p       total air pressure [Pa]

    Returns
    -------
    water vapor pressure [Pa]
    """
    c = eps()  # Rd / Rv = 0.622
    return (q * p) / (q + (1 - q) * c)


def sh2ppmv(q, p):
    """Convert specific humidity [kg/kg] to parts per million by volume ppmv
    required for RTTOV

    Parameters
    ----------
    q       specific humidity [kg/kg]
    p       air pressure [Pa]

    Returns
    -------
    pressure [ppmv]
    """
    # c = con.eps() # Rd/Rv
    # Mixing Ratio
    # r = q / (1 - q)
    # aus e = r p /(r+c)
    # e = r * p / (r + c)
    # e = (q * p)/(q + (1 - q)*c)
    # Volumen-spez.
    e = sh2vap(q, p)
    # 10e6 * Pwater / (Ptotal - Pwater)
    return vap2ppmv(e, p)


def ppmv2pa(x, p):
    """Convert ppmv to Pa

    Parameters
    ----------
    x       Gas pressure [ppmv]
    p       total air pressure [Pa]

    Returns
    -------
    pressure [Pa]
    """
    return x * p / (1e6 + x)


def vap2sh(e, p):
    """ Convert water vapor pressure to specific humidity
    Parameters
    ----------
    e      Water vapor [Pa]
    p      total air pressure [Pa]

    Returns
    -------
    specific humidity (1 = kg/kg)
    """
    c = eps()  # Rd/Rv
    pa = p - e  # dry air pressure
    return (e * c) / (e * c + pa)


def vap2mix(e, p):
    """ convert water vapor pressure to mixing ratio

    Parameters
    ----------
    e       Water vapor (Pa)
    p       air pressure (Pa)

    Returns
    -------
    mixing ratio (1)
    """
    c = eps()  # Rd/Rv
    return c * e / (p - e)


def mix2vap(x, p):
    """ convert mixing ratio to water vapor pressure

    Parameters
    ----------
    x       Mixing ratio (1)
    p       air pressure (Pa)

    Returns
    -------
    water vapor pressure (Pa)
    """
    c = eps()
    return x * p / (c + x)


def rh2vap(r, t, method='HylandWexler', **kwargs):
    """Relative Humidity to Water vapor pressure
    Using Saturation Water vapor pressure function dependent on T

    Args:
        r: relative humidity in [0-1]
        t: Temperatur in [K]
        method: Saturation water vapor pressure function name

    Returns:
        e: water vapor pressure in [Pa]
    """
    return r * svp(t, method=method, **kwargs)


def rh2ppmv(r, t, p, method='HylandWexler', **kwargs):
    """Relative humidity to ppmv

    Parameters
    ----------
    r       relative humidity in [0-1]
    t       Temperatur in [K]
    p       air pressure in [Pa]
    esatfunc    saturation water vapor pressure function name

    Returns
    -------
    pressure in [ppmv]
    """
    e = r * svp(t, method=method, **kwargs)
    return vap2ppmv(e, p)


def frostpoint(e, method=None):
    try:
        if method is not None:
            if callable(method):
                fpfunc = method
            else:
                fpfunc = eval(method)

        else:
            fpfunc = eval('frostpoint_Magnus')
        return fpfunc(e)
    except:
        import sys
        print("Functions: ", ", ".join([i for i in dir(sys.modules[__name__]) if 'frostpoint_' in i]))



def frostpoint_Magnus(e):
    return np.where(e > 0,
                    273.15 + 272.62 * np.log(e / 611.2) / (22.46 - np.log(e / 611.2)),
                    np.nan)


def frostpoint_Sonntag(e):
    y = np.where(e > 0, np.log(e / 611.213), np.nan)
    return 273.15 + 13.7204 * y \
           + 7.36631e-1 * y * y \
           + 3.32136e-2 * y * y * y \
           + 7.78591e-4 * y * y * y * y


def frostpoint_MurphyKoop(e):
    """Frost Point Temperature from water vapor pressure
    after Murphy and Koop, 2005
    T > 115 K

    Input: e [Pa] Water Vapor Pressure
    Output: T_frost [K]
    """
    return np.where(e > 0, (1.814625 * np.log(e) + 6190.134) / (29.120 - np.log(e)), np.nan)


def compare_dewpoints():
    from .esat import reference_table

    data = reference_table()

    methods = ['dewpoint_Bolton', 'dewpoint_Boegel', 'dewpoint_AlduchovEskridge',
               'dewpoint_Sonntag', 'dewpoint_Magnus', 'dewpoint_ecmwf', 'dewpoint_JMA',
               'dewpoint_newton', 'dewpoint_minimize']

    for im in methods:
        try:
            data[im] = eval(im)(data['ew'], verbose=1)
            data[im] -= data['t']

        except Exception as e:
            print("Error: ", im, repr(e))
            data[im] = np.nan
            continue
    return data


def performance_dewpoints(esat=None):
    import timeit
    methods = ['dewpoint_Bolton', 'dewpoint_Boegel', 'dewpoint_AlduchovEskridge',
               'dewpoint_Sonntag', 'dewpoint_Magnus', 'dewpoint_ecmwf', 'dewpoint_JMA',
               'dewpoint_newton', 'dewpoint_minimize']

    for im in methods:
        try:
            def doit():
                func = eval(im)
                data = np.random.randint(1, 10000, 1000) / 10.
                _ = func(data, method=esat)

            t1 = timeit.Timer(doit)
            print(im, t1.timeit(100)/100.,'s')

        except Exception as e:
            print("Error: ", im, repr(e))
            continue


def dewpoint(e, method='dewpoint_Bolton', **kwargs):
    try:
        if callable(method):
            fpfunc = method
        else:
            fpfunc = eval(method)

        return fpfunc(e, **kwargs)
    except:
        import sys
        print("Functions: ", ", ".join([i for i in dir(sys.modules[__name__]) if 'dewpoint_' in i]))


def dewpoint_Bolton(e, **kwargs):
    """ Inverse Bolton Function to yield temperature
    Bolton 1981

    Input: e [Pa]
    Output: T [K]
    """
    a = 17.67
    b = 243.5
    return np.where(e > 0,
                    273.15 - b * np.log(e / 611.12) / (np.log(e / 611.12) - a),
                    np.nan)


def dewpoint_AlduchovEskridge(e, **kwargs):
    """ Alduchov and Eskridge(1996)

    Input: e [Pa]
    Output: T [K]
    """
    a = 17.625
    b = 243.04
    return np.where(e > 0,
                    273.15 - b * np.log(e / 611.12) / (np.log(e / 611.12) - a),
                    np.nan)


def dewpoint_Boegel(e, **kwargs):
    a = 611.21
    b = 18.729
    c = 257.87
    d = 227.3

    return 273.15 + np.where(e > 0,
                             d / 2. * (b - np.log(e / a) - (
                                     (b - np.log(e / a)) ** 2 - 4 * c * np.log(e / a) / d) ** 0.5),
                             np.nan)


def dewpoint_Sonntag(e, **kwargs):
    """Calculation of Dewpoint after Sonntag 1994

    Parameters
    ----------
    e       Water vapor pressure [Pa]

    Returns
    -------
    dewpoint in [K]
    """
    y = np.where(e > 0, np.log(e / 611.213), np.nan)
    return 273.15 + 13.715 * y \
           + 8.4262e-1 * y * y \
           + 1.9048e-2 * y * y * y \
           + 7.8158e-3 * y * y * y * y


def dewpoint_Magnus(e, r=None, t=None, **kwargs):
    """Dewpoint after Magnus Formulation

    Parameters
    ----------
    r       relative Humidity [0-1]
    t       Temperatur [K]

    Returns
    -------
    dewpoint in [K]
    """
    # alpha = np.log(r) + 17.27 * t / (237.7 + t)
    # return 273.15 + 237.7 * alpha / (17.27 - alpha)
    a = 17.27
    b = 237.7
    if r is not None and t is not None:
        t = t - 273.15  # requires Degrees
        return np.where(r > 0, 273.15 + b * (np.log(r) + a * t / (b + t))
                        / (a - (np.log(r) + a * t / (b + t))), np.nan)
    return np.where(e > 0,
                    273.15 - b * np.log(e / 611.12) / (np.log(e / 611.12) - a),
                    np.nan)


def dewpoint_ECMWF(e, t=None, rh=None, **kwargs):
    """ ECMWF IFS CYC31R1 Data Assimilation Documentation (Page 86-88)
    Derived Variables from TEMP
    Temperature and dew point are transformed into realtive humidity (RH) for
    TEMP observations, with a further transformation of RH into specific humidity (Q)
    for TEMP observations.

    Uses foeewmo for water only saturation water vapor pressure

    Parameters
    ----------
    t       temperature     [K]
    rh      rel. humidity   [1]

    Returns
    -------
    td      dew point       [K]
    """
    if t is not None:
        if rh is None:
            raise ValueError('requires t and rh or e')

        e = svp(t, method='FOEEWMO') * rh
    lnpart = np.log(e / 611.21)
    return (17.502 * 273.16 - 32.19 * lnpart) / (17.502 - lnpart)


def dewpoint_JMA(e, method='HylandWexler', tol=0.001, **kwargs):
    """ Algorithms used by electronic logbooks for the computation of dew point temperature

    OBSJMA Japan Meteorological Agency

    https://www.wmo.int/pages/prog/amp/mmop/JCOMM/OPA/SOT/documents/Dew-point-algorithm-OBSJMA.pdf

    Parameters
    ----------
    e           water vapor pressure    [Pa]
    method      saturation water vapor method
    tol         allowed tolerance in K

    Returns
    -------
    td          dew point       [K]
    """
    from ..fun.cal import fuzzy_equal
    verbose = kwargs.get('verbose', 0)

    @np.vectorize
    def dvp(ex):
        if not np.isfinite(ex):
            return np.nan
        n = 0
        t1 = 150.15
        t2 = 400.15
        d1 = svp(t1, method=method, **kwargs) - ex
        while not fuzzy_equal(abs(d1), 0, tol) and n < 100:
            d1 = svp(t1, method=method, **kwargs) - ex
            d2 = svp(t2, method=method, **kwargs) - ex
            if (d1 * d2 < 0):
                t2 = t2 - (t2 - t1) * 0.5
            else:
                t3 = t1
                t1 = t2
                t2 = t2 + (t2 - t3) * 0.5
            n += 1

        if verbose > 0:
            print("Warning Limit reached: ", n, t1, t2, svp(t1, method=method, **kwargs) - ex)
        return t1

    return dvp(e)


def dewpoint_minimize(e, method='HylandWexler', tol=0.001, **kwargs):
    from scipy import optimize

    verbose = kwargs.get('verbose', 0)
    if method is None:
        method = 'HylandWexler'
    @np.vectorize
    def dvp(ex):
        res = optimize.minimize_scalar(lambda x: abs(svp(x, method=method, **kwargs) - ex),
                                       bounds=(183, 383), options={'xatol': tol},
                                       method='bounded')
        if verbose > 0:
            # print(res)
            print("Min: %s %f %s %d" % (method, tol, res.success, res.nfev))
        return res.x

    return np.where(e > 0, dvp(e), np.nan)


def dewpoint_newton(e, method='HylandWexler', tol=0.001, maxiter=100, **kwargs):
    """Inverse Calculation of Dewpoint
    Calculate Dew point from Saturation water vapor pressure function, inverted by Newton Minimization
    """
    from scipy import optimize

    verbose = kwargs.get('verbose', 0)
    if method is None:
        method = 'HylandWexler'

    @np.vectorize
    def dvp(ex):
        res = optimize.newton(lambda x: abs(svp(x, method=method, **kwargs) - ex),
                              x0=dewpoint_Bolton(ex),
                              tol=tol, maxiter=maxiter)
        return res
    if verbose > 0:
        print("Newton: %s %f %d" % (method, tol, maxiter))

    return np.where(e > 0, dvp(e), np.nan)

#
# def humidity(x, t=288.15, p=101325, in_type='specific humidity', out_type='water vapor', nounits=False,
#              method=None, **kwargs):
#     """
#     Convert Humidity Data into other humidity xData
#     humidity:
#     - partial pressure of water vapor in Pa (not hPa nor mbar)
#     - specific humidity in kg/kg (not g/kg)
#     - mixing ratio in kg/kg (not g/kg)
#     - relative humidity in percent
#     - dew point temperature in K (not degree Celsius)
#     - virtual temperature in K (not degree Celsius)
#     """
#     in_unit = None
#     out_unit = None
#     ispandas = False
#     inputs = ['water vapor', 'specific humidity', 'relative humidity', 'virtual temperature', 'ppmv']
#     outputs = ['water vapor', 'specific humidity', 'relative humidity', 'mixing ratio', 'ppmv', 'virtual temperature',
#                'dewpoint']
#
#     if in_type not in inputs:
#         print "Unkown: %s" % in_type
#         print inputs
#
#     if out_type not in outputs:
#         print "Unkown: %s" % out_type
#         print outputs
#
#     # retain time information
#     if isinstance(x, pd.DataFrame) | isinstance(x, pd.Series):
#         # convert to numpy
#         time = x.index
#         if isinstance(x, pd.DataFrame):
#             columns = x.columns
#         x = np.ma.asarray(x)
#         ispandas = True
#     if isinstance(p, pd.DataFrame) | isinstance(p, pd.Series):
#         p = np.ma.asarray(p)
#     # Vectorize funktion or not?
#     x = np.ma.asarray(x)
#
#     c = 0.622  # Constant (Mv / Md)
#     if in_type in ('partial pressure', 'partial pressure of water vapor', 'water vapor'):
#         e = x
#         in_unit = 'Pa'
#     elif in_type == 'specific humidity':
#         q = x
#         e = (q * p) / (c + (1 - c) * q)
#         in_unit = 'kg/kg + Pa (p)'
#     elif in_type == 'mixing ratio':
#         r = x
#         e = (r * p) / (r + c)
#         in_unit = '1 + Pa (p)'
#     elif in_type == 'relative humidity':
#         RH = x
#         es = svp(t, method=method, **kwargs)
#         e = (RH / 100.) * es
#         in_unit = '% + K (t)'
#     elif in_type == 'virtual temperature':
#         Tv = x
#         e = p * (1 - t / Tv) / (1 - c)
#         in_unit = 'K + K (t)'
#     elif in_type == 'ppmv':
#         x = (x / 1e6)
#         e = x * p / (1 + x)
#         in_unit = 'ppmv + Pa (p)'
#     else:
#         print "[INPUT] Unknown Conversion"
#         return
#
#     if out_type in ('partial pressure', 'partial pressure of water vapor', 'water vapor'):
#         out = e
#         out_unit = 'Pa'
#     elif out_type == 'specific humidity':
#         Pd = p - e
#         out = (e * c) / (e * c + Pd)
#         out_unit = '+ Pa (p) = kg/kg '
#     elif out_type == 'mixing ratio':
#         Pd = p - e
#         out = (e / Pd) * c
#         out_unit = '+ Pa (p) = 1 '
#     elif out_type == 'ppmv':
#         Pd = p - e
#         out = (e / Pd) * 1e6
#         out_unit = '+ Pa (p) = ppmv'
#     elif out_type == 'relative humidity':
#         es = svp(t, method=method, **kwargs)
#         out = 100. * e / es
#         out_unit = '+ K (t) = %'
#     elif out_type == 'virtual temperature':
#         out = t / (1 - (e / p) * (1 - c))
#         out_unit = '+K (t) + Pa (p) = K'
#     elif out_type == 'dew point':
#         out = dewpoint(e, method=method)
#         out_unit = 'K'
#     else:
#         print "[OUTPUT] Unknown Conversion"
#         return
#     if ispandas:
#         if len(out.shape) > 1:
#             out = pd.DataFrame(out, index=time, columns=columns)
#         else:
#             out = pd.Series(out, index=time)
#
#     if nounits:
#         return out
#     else:
#         print "Conversion: %s => %s" % (in_unit, out_unit)
#         print in_type, ' => ', out_type
#         return out
