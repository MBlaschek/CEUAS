# -*- coding: utf-8 -*-


__all__ = ['temperature', 'specific_humidity', 'relative_humidity', 'water_vapor_pressure',
           'dew_point_depression']


def temperature(data, dim='plev', report=True, **kwargs):
    """ Quality control of Temperatures with RTTOV Profile limits

    Args:
        data (DataArray): temperatures [K]
        dim (str): pressure
        report (bool): return report
        **kwargs:

    Info:
        0 no_qc
        1 good_data
        2 outside_range
        3 dpd30

    Returns:
        DataArray : input data (unchanged)
        DataArray : qc [0, 1, 2, 3]

    """
    import numpy as np
    from xarray import DataArray
    from .. import config
    from ..fun import message, load_rttov

    if not isinstance(data, DataArray):
        raise ValueError("Requires a DataArray, ", type(data))

    if config.rttov_profile_limits is None:
        load_rttov()

    if 'units' in data.attrs:
        if data.attrs['units'] != 'K':
            raise RuntimeWarning("Temperature are not in Kelvin")

    # Valid Range: [183.15, 333.15] K
    rt = config.rttov_profile_limits
    pressure = data[dim].values
    tmin = rt.Temperature_min.values
    tmax = rt.Temperature_max.values
    pin = None
    if 'units' in data[dim].attrs:
        if data[dim].attrs['units'] == 'hPa':
            pin = rt.Pressure.values  # hPa to hPa
    if pin is None:
        pin = rt.Pressure.values * 100.  # hPa to Pa

    tmins = np.interp(np.log(pressure), np.log(pin), tmin, left=tmin.min(), right=tmin.max())
    tmins = np.broadcast_to(tmins, data.values.shape)
    tmaxs = np.interp(np.log(pressure), np.log(pin), tmax, left=tmax.min(), right=tmax.max())
    tmaxs = np.broadcast_to(tmaxs, data.values.shape)
    with np.errstate(invalid='ignore'):
        logic = (data.values < tmins) | (data.values > tmaxs)

    tcor = np.sum(logic)  # number of corrected values
    message("Temperatures %d [%d]" % (int(tcor), int(np.sum(np.isfinite(data.values)))), **kwargs)
    if 'QC' in data.attrs.keys():
        tcor += data.attrs['QC']

    if data.name is None:
        data.name = 't'
    # 0 no_qc 1 good_data 2 outside_range 3 dpd30
    flags = DataArray(data=np.isfinite(data.values), coords=data.coords, dims=data.dims,
                      name=data.name + '_qc',
                      attrs={'QC': tcor, 'flag_values': (0, 1, 2, 3),
                             'flag_meanings': "no_qc good_data outside_range dpd30",
                             'valid_range': (0, 3)})
    flags = flags.where(logic == 0, 2)
    # flags.values[logic] = 2  # outside range / overwritten

    if report:
        rep = qcreport(flags, data=data, min=tmins, max=tmaxs)

    data.attrs['QC'] = tcor
    data.attrs['ancilliary_variables'] = data.name + '_qc'
    if report:
        return data, flags, rep
    return data, flags


def specific_humidity(data, dim='plev', report=True, **kwargs):
    """ Quality control of spec. humidity with RTTOv Profile limits

    Args:
        data (DataArray): specifc humidity [kg/kg]
        dim (str): pressure
        report (bool): return report
        **kwargs:

    Info:
        0 no_qc
        1 good_data
        2 outside_range
        3 dpd30

    Returns:
        DataArray : input data (unchanged)
        DataArray : qc [0, 1, 2, 3]

    """
    import numpy as np
    from xarray import DataArray
    from .humidity import vap2sh
    from .. import config
    from ..fun import message, load_rttov

    if not isinstance(data, DataArray):
        raise ValueError("Requires a DataArray, ", type(data))

    if config.rttov_profile_limits is None:
        load_rttov()

    if 'units' in data.attrs:
        if data.attrs['units'] != 'kg/kg':
            raise RuntimeWarning("specific humidity is not in kg/kg")

    # Valid Range: [0, 1]  kg/kg
    rt = config.rttov_profile_limits
    pressure = data[dim].values
    qmin = vap2sh(rt.Water_vapour_min.values, rt.Pressure.values)
    qmax = vap2sh(rt.Water_vapour_max.values, rt.Pressure.values)
    pin = None
    if 'units' in data[dim].attrs:
        if data[dim].attrs['units'] == 'hPa':
            pin = rt.Pressure.values  # hPa to hPa
    if pin is None:
        pin = rt.Pressure.values * 100.  # hPa to Pa

    qmins = np.interp(np.log(pressure), np.log(pin), qmin, left=qmin.min(), right=qmin.max())
    qmins = np.broadcast_to(qmins, data.values.shape)
    qmaxs = np.interp(np.log(pressure), np.log(pin), qmax, left=qmax.min(), right=qmax.max())
    qmaxs = np.broadcast_to(qmaxs, data.values.shape)
    with np.errstate(invalid='ignore'):
        logic = (data.values < qmins) | (data.values > qmaxs)

    qcor = np.sum(logic)  # number of corrected values
    message("Spec. Humidity %d [%d]" % (int(qcor), int(np.sum(np.isfinite(data.values)))), **kwargs)
    if 'QC' in data.attrs.keys():
        qcor += data.attrs['QC']

    if data.name is None:
        data.name = 'q'

    # 0 no_qc 1 good_data 2 outside_range 3 dpd30
    flags = DataArray(data=np.isfinite(data.values), coords=data.coords, dims=data.dims,
                      name=data.name + '_qc',
                      attrs={'QC': qcor, 'flag_values': (0, 1, 2, 3),
                             'flag_meanings': "no_qc good_data outside_range dpd30",
                             'valid_range': (0, 3)})
    flags = flags.where(logic == 0, 2)
    # flags.values[logic] = 2  # outside range

    if report:
        rep = qcreport(flags, data=data, min=qmins, max=qmaxs)

    data.attrs['QC'] = qcor
    data.attrs['ancilliary_variables'] = data.name + '_qc'
    if report:
        return data, flags, rep
    return data, flags


def relative_humidity(data, report=True, **kwargs):
    """Quality control of rel. humidity with RTTOV Profile limits

    Args:
        data (DataArray): rel. humidity [1]
        report (bool): return report
        **kwargs:

    Info:
        0 no_qc
        1 good_data
        2 outside_range
        3 dpd30

    Returns:
        DataArray : input data (unchanged)
        DataArray : qc [0, 1, 2, 3]

    """
    import numpy as np
    from xarray import DataArray
    from ..fun import message

    if not isinstance(data, DataArray):
        raise ValueError("Requires a DataArray, ", type(data))

    if 'units' in data.attrs:
        if data.attrs['units'] != '1':
            raise RuntimeWarning("rel. humidity has units 1 (ratio)")

    r_absmin = 0
    r_absmax = 1.03  # 3 % plus
    # Valid Range: [ 0; 1] ratio
    with np.errstate(invalid='ignore'):
        logic = (data.values < r_absmin) | (data.values > r_absmax)

    rcor = np.sum(logic)  # number of corrected values
    message("Rel. Humidity %d [%d]" % (int(rcor), int(np.sum(np.isfinite(data.values)))), **kwargs)
    if 'QC' in data.attrs.keys():
        rcor += data.attrs['QC']

    if data.name is None:
        data.name = 'rh'

    # 0 no_qc 1 good_data 2 outside_range 3 dpd30
    flags = DataArray(data=np.isfinite(data.values), coords=data.coords, dims=data.dims,
                      name=data.name + '_qc',
                      attrs={'QC': rcor, 'flag_values': (0, 1, 2, 3),
                             'flag_meanings': "no_qc good_data outside_range dpd30",
                             'valid_range': (0, 3)})
    flags = flags.where(logic == 0, 2)
    # flags.values[logic] = 2  # outside range
    if report:
        rep = qcreport(flags, data=data, absmin=r_absmin, absmax=r_absmax)

    data.attrs['QC'] = rcor
    data.attrs['ancilliary_variables'] = data.name + '_qc'
    if report:
        return data, flags, rep
    return data, flags


def dew_point_depression(data, dim='time', iyear=1995, report=True, **kwargs):
    """Quality control of dewpoint depression with RTTOV Profile limits

    Args:
        data (DataArray): dewpoint depression [K]
        dim (str): time dimension
        iyear (int): year of DPD30
        report (bool): return report
        **kwargs:

    Info:
        0 no_qc
        1 good_data
        2 outside_range
        3 dpd30

    Returns:
        DataArray : input data (unchanged)
        DataArray : qc [0, 1, 2, 3]

    """
    import numpy as np
    from xarray import DataArray
    from ..fun import message

    if not isinstance(data, DataArray):
        raise ValueError("Requires a DataArray, ", type(data))

    if 'units' in data.attrs:
        if data.attrs['units'] != 'K':
            raise RuntimeWarning("Dewpoint depression is not in Kelvin")

    dpd_absmin = 0
    dpd_absmax = 80
    # Valid Range: [0, 80] K
    with np.errstate(invalid='ignore'):
        logic = (data.values < dpd_absmin) | (data.values > dpd_absmax)

    dpdcor = np.sum(logic)  # number of corrected values
    message("Dew point depression %d [%d]" % (int(dpdcor), int(np.sum(np.isfinite(data.values)))), **kwargs)
    if 'QC' in data.attrs.keys():
        dpdcor += data.attrs['QC']

    if data.name is None:
        data.name = 'dpd'

    # 0 no_qc 1 good_data 2 outside_range 3 dpd30
    # 0 also for NaN
    flags = DataArray(data=np.isfinite(data.values).astype(int), coords=data.coords, dims=data.dims,
                      name=data.name + '_qc',
                      attrs={'QC': dpdcor, 'flag_values': (0, 1, 2, 3),
                             'flag_meanings': "no_qc good_data outside_range dpd30",
                             'valid_range': (0, 3)})
    flags = flags.where(logic == 0, 2)
    # flags.values[logic] = 2  # outside range / overwritten
    #
    # DPD30
    #
    with np.errstate(invalid='ignore'):
        logic = data.where(data['{}.year'.format(dim)] <= iyear).groupby('{}.year'.format(dim)).apply(xhist, **kwargs)

    flags = flags.where(logic == 0, 3)
    # flags.values[logic.values.astype(bool)] = 3  # DPD 30

    if report:
        rep = qcreport(flags, data=data, absmin=dpd_absmin, absmax=dpd_absmax)
        rep.extend(qcreport(flags, iflag=3, data=data, flag_name='DPD30'))

    data.attrs['QC'] = dpdcor + logic.values.sum()
    data.attrs['ancilliary_variables'] = data.name + '_qc'
    if report:
        return data, flags, rep
    return data, flags


def water_vapor_pressure(data, dim='plev', report=True, **kwargs):
    """ Quality control of water vapor pressure with RTTOV Profile limits

    Args:
        data (DataArray): water vapor [Pa]
        dim (str): pressure
        report (bool): return report
        **kwargs:

    Info:
        0 no_qc
        1 good_data
        2 outside_range
        3 dpd30

    Returns:
        DataArray : input data (unchanged)
        DataArray : qc [0, 1, 2, 3]

    """
    import numpy as np
    from xarray import DataArray
    from .humidity import ppmv2pa
    from .. import config
    from ..fun import message, load_rttov

    if not isinstance(data, DataArray):
        raise ValueError("Requires a DataArray, ", type(data))

    if config.rttov_profile_limits is None:
        load_rttov()

    if 'units' in data.attrs:
        if data.attrs['units'] != 'Pa':
            raise RuntimeWarning("Water vapor are not in Pa")

    # Valid Range: [183.15, 333.15] K
    rt = config.rttov_profile_limits
    pressure = data[dim].values
    pfactor = 100  # Pa
    u1 = 'Pa'
    if 'units' in data[dim].attrs:
        u1 = data.attrs['units']
        if u1 == 'hPa':
            pfactor = 1

    pin = rt.Pressure.values * pfactor
    #
    # Check Water vapor Unit [hPa or Pa] ?
    #
    dfactor = 1  # Pa
    u2 = 'Pa'
    if 'units' in data.attrs:
        u2 = data.attrs['units']
        if u2 == 'hPa':
            dfactor = 0.01
    message("VP [{}]: {} , P [{}]: {}".format(u2, dfactor, u1, pfactor), **kwargs)
    #
    # Pressure is in hPa -> pfactor
    # Water vapor can be hPa or Pa -> dfactor
    #
    vpmin = ppmv2pa(rt.Water_vapour_min.values, rt.Pressure.values * pfactor) * dfactor
    vpmax = ppmv2pa(rt.Water_vapour_max.values, rt.Pressure.values * pfactor) * dfactor
    #
    #
    #
    vpmins = np.interp(np.log(pressure), np.log(pin), vpmin, left=vpmin.min(), right=vpmin.max())
    vpmins = np.broadcast_to(vpmins, data.values.shape)
    vpmaxs = np.interp(np.log(pressure), np.log(pin), vpmax, left=vpmax.min(), right=vpmax.max())
    vpmaxs = np.broadcast_to(vpmaxs, data.values.shape)
    with np.errstate(invalid='ignore'):
        logic = (data.values < vpmins) | (data.values > vpmaxs)

    vpcor = np.sum(logic)  # number of corrected values
    message("Water vapor pressure %d [%d]" % (int(vpcor), int(np.sum(np.isfinite(data.values)))), **kwargs)
    if 'QC' in data.attrs.keys():
        vpcor += data.attrs['QC']

    if data.name is None:
        data.name = 'vp'
    # 0 no_qc 1 good_data 2 outside_range 3 dpd30
    flags = DataArray(data=np.isfinite(data.values), coords=data.coords, dims=data.dims,
                      name=data.name + '_qc',
                      attrs={'QC': vpcor, 'flag_values': (0, 1, 2, 3),
                             'flag_meanings': "no_qc good_data outside_range dpd30",
                             'valid_range': (0, 3)})
    flags = flags.where(logic == 0, 2)
    # flags.values[logic] = 2  # outside range
    if report:
        rep = qcreport(flags, data=data, min=vpmins, max=vpmaxs)

    data.attrs['QC'] = vpcor
    data.attrs['ancilliary_variables'] = data.name + '_qc'
    if report:
        return data, flags, rep
    return data, flags


def wind_speed(data, report=True, **kwargs):
    import numpy as np
    from xarray import DataArray
    from ..fun import message

    if not isinstance(data, DataArray):
        raise ValueError("Requires a DataArray, ", type(data))

    if 'units' in data.attrs:
        if data.attrs['units'] not in ['m s-1', 'm/s']:
            raise RuntimeWarning("wind speed does not have units m/s (ms^-1)", data.attrs['units'])

    ws_absmin = 0
    ws_absmax = 150
    # Valid Range: [ 0; 150] m/s
    with np.errstate(invalid='ignore'):
        logic = (data.values < ws_absmin) | (data.values > ws_absmax)

    rcor = np.sum(logic)  # number of corrected values
    message("Wind Speed %d [%d]" % (int(rcor), int(np.sum(np.isfinite(data.values)))), **kwargs)
    if 'QC' in data.attrs.keys():
        rcor += data.attrs['QC']

    if data.name is None:
        data.name = 'ws'

    # 0 no_qc 1 good_data 2 outside_range 3 dpd30
    flags = DataArray(data=np.isfinite(data.values), coords=data.coords, dims=data.dims,
                      name=data.name + '_qc',
                      attrs={'QC': rcor, 'flag_values': (0, 1, 2, 3),
                             'flag_meanings': "no_qc good_data outside_range dpd30",
                             'valid_range': (0, 3)})
    flags = flags.where(logic == 0, 2)
    # flags.values[logic] = 2  # outside range
    if report:
        rep = qcreport(flags, data=data, absmin=ws_absmin, absmax=ws_absmax)

    data.attrs['QC'] = rcor
    data.attrs['ancilliary_variables'] = data.name + '_qc'
    if report:
        return data, flags, rep
    return data, flags


def wind_direction(data, report=True, **kwargs):
    import numpy as np
    from xarray import DataArray
    from ..fun import message

    if not isinstance(data, DataArray):
        raise ValueError("Requires a DataArray, ", type(data))

    if 'units' in data.attrs:
        if data.attrs['units'] != 'degree':
            raise RuntimeWarning("wind speed does not have units degree (degree)", data.attrs['units'])

    ws_absmin = 0
    ws_absmax = 360

    with np.errstate(invalid='ignore'):
        logic = (data.values < ws_absmin) | (data.values > ws_absmax)

    rcor = np.sum(logic)  # number of corrected values
    message("Wind Direction %d [%d]" % (int(rcor), int(np.sum(np.isfinite(data.values)))), **kwargs)
    if 'QC' in data.attrs.keys():
        rcor += data.attrs['QC']

    if data.name is None:
        data.name = 'wd'

    # 0 no_qc 1 good_data 2 outside_range 3 dpd30
    flags = DataArray(data=np.isfinite(data.values), coords=data.coords, dims=data.dims,
                      name=data.name + '_qc',
                      attrs={'QC': rcor, 'flag_values': (0, 1, 2, 3),
                             'flag_meanings': "no_qc good_data outside_range dpd30",
                             'valid_range': (0, 3)})
    flags = flags.where(logic == 0, 2)
    # flags.values[logic] = 2  # outside range
    if report:
        rep = qcreport(flags, data=data, absmin=ws_absmin, absmax=ws_absmax)

    data.attrs['QC'] = rcor
    data.attrs['ancilliary_variables'] = data.name + '_qc'
    if report:
        return data, flags, rep
    return data, flags


def qcreport(flags, iflag=2, data=None, min=None, max=None, absmin=None, absmax=None, flag_name=None):
    from numpy import where, char, array, datetime64
    # search data (logic), produce a list of flagged values
    events = []
    order = list(flags.dims)
    indices = where(flags == iflag)
    if indices[0].size > 0:
        events = array("[QC] " + flags.name + " | (")
        for i, idim in enumerate(order):
            if isinstance(flags[idim].values[0], datetime64):
                txt = idim + "='{}'"
            else:
                txt = idim + "={}"
            txt += ',' if i < len(order) - 1 else ')'
            events = char.add(events, array(list(map(txt.format, flags[idim].values[indices[i]]))))
        if flag_name is not None:
            events = char.add(events, array(" | " + flag_name))

        if data is not None:
            if min is not None:
                events = char.add(events, array(list(map(" | {:.4f} < ".format, min[indices]))))
            elif absmin is not None:
                events = char.add(events, array(" | {} <".format(absmin)))
            else:
                events = char.add(events, array(" | "))
            if data is not None:
                events = char.add(events, array(list(map("{:.4f}".format, data.values[indices]))))
            if max is not None:
                events = char.add(events, array(list(map(" < {:.4f}".format, max[indices]))))
            if absmax is not None:
                events = char.add(events, array("> {}".format(absmax)))
        events = list(events)
    return events


def combine_flags(a, b):
    a_attrs = {}
    b_attrs = {}
    if 'flag_values' in a.attrs.keys() and 'flag_meanings' in a.attrs.keys():
        a_attrs = dict(zip(a.attrs['flag_values'], a.attrs['flag_meanings'].split(' ')))
    if 'flag_values' in b.attrs.keys() and 'flag_meanings' in b.attrs.keys():
        b_attrs = dict(zip(b.attrs['flag_values'], b.attrs['flag_meanings'].split(' ')))
    #
    # Same flags + meanings ?
    #
    m_attrs = {}
    n = max(a_attrs.keys())
    for i, j in a_attrs.items():
        if i in b_attrs.keys():
            if j == b_attrs[i]:
                m_attrs[i] = j
            else:
                m_attrs[i] = j
                m_attrs[n + 1] = b_attrs[i]
                b = b.where(b != i, n + 1)
                n += 1
    for i, j in b_attrs.items():
        if i not in m_attrs.keys():
            m_attrs[i] = j

    # merged = a.combine_first(b)
    merged = a + b
    merged.attrs['flag_values'] = tuple(m_attrs.keys())
    merged.attrs['flag_meanings'] = " ".join(m_attrs.values())
    merged.attrs['valid_range'] = (0, n)
    return merged


def xhist(data, value=30, bins=None, excess=5, eps=0.5, **kwargs):
    """ Customized histogram

    Args:
        data:
        value:
        bins:
        normed:

    Returns:
        ndarray : boolean mask
    """
    from numpy import arange, linspace, isfinite, histogram, argmax, errstate
    from ..fun import message

    if bins is None:
        bins = arange(eps, value + value / 2. + eps, 2 * eps)
    elif isinstance(bins, int):
        bins = linspace(eps, value + value / 2. + eps, num=bins)
    else:
        pass
    out = data.copy()
    out.values[:] = 0
    out.values = out.values.astype(int)
    itx = isfinite(data.values)  # remove missing
    counts, divs = histogram(data.values[itx], bins=bins, density=True)
    #
    # is there enough data to make a histogram ?
    #
    n = len(counts)
    if itx.sum() < n * 2:
        return out

    #
    # mark values with anomalous high distribution! add 10% or 5%
    # bins: [ ) borders, right=False
    #
    jtx = argmax(histogram(value, bins=divs)[0])  # box of value in the current histogram
    crit = (excess / 100.) * (n / 100.)  # 5%  depends on sample size
    # Value is above critical
    #
    if counts[jtx] > crit:
        # Before and After are below critical
        a = counts[jtx - 1] if jtx - 1 >= 0 else 0
        b = counts[jtx + 1] if jtx + 1 < n else 0
        message("[%.2f pm %.2f] %.2f [%.2f] %.2f > %.2f" % (value, eps, a, counts[jtx], b, crit), **kwargs)
        if a < crit and b < crit:
            with errstate(invalid='ignore'):
                out.values = (itx & (data.values > (value - eps)) & (data.values < (value + eps))).astype(int)

    return out


def flagcount(data):
    from numpy import bincount, ravel
    counts = bincount(ravel(data.values).astype(int))
    names = [''] * len(counts)
    values = list(range(len(counts)))
    if 'flag_meanings' in data.attrs.keys():
        names = data.attrs['flag_meanings'].split(' ')
    if 'flag_values' in data.attrs.keys():
        values = data.attrs['flag_values']
    txt = ""
    j = 0
    for i in range(len(counts)):
        if i in values:
            txt += "{:2d} {} : {:6d}".format(j, names[j], counts[j])
            txt += "\n"
            j += 1

    return txt

