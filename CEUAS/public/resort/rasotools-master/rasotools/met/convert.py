# -*- coding: utf-8 -*-

__all__ = ['to_rh', 'to_vp', 'to_dpd', 'to_sh', 'to_dewpoint']


def to_rh(temp, dpd=None, spec_humi=None, press=None, method='HylandWexler', precision=6, **kwargs):
    """
    Convert dewpoint departure or specific humidity to relative humidity

    1. DPD to RH (T):
        + DPD to VP, T to VPsat
        + VP / VPsat

    2. Q to RH (Q, T, P):
        + Q to VP, T to VPsat
        + VP / VPsat

    Args:
        temp (DataArray): temperature
        dpd (DataArray): dewpoint depression
        spec_humi (DataArray): specific humidity
        press (DataArray,str): air pressure or dim name
        method (str): Saturation water vapor pressure formulation
        precision (int): decimal precision

    Returns:
        DataArray : relative humidity [1]
    """
    from numpy import around, errstate
    from xarray import DataArray
    from .esat import svp
    from .humidity import sh2vap
    from ..fun import message, leveldown

    if not isinstance(temp, DataArray):
        raise ValueError("Requires a DataArray", type(temp))

    if dpd is None and spec_humi is None:
        raise RuntimeError("Requires either dpd or q for conversion")

    if dpd is not None and not isinstance(dpd, DataArray):
        raise ValueError("Requires a DataArray", type(dpd))

    if spec_humi is not None:
        if not isinstance(spec_humi, DataArray):
            raise ValueError("Requires a DataArray", type(spec_humi))

        if press is None:
            raise RuntimeError("Conversion Q>RH requires a pressure variable as well")

    rvar = temp.copy()
    if press is not None:
        if isinstance(press, str):
            if press in rvar.dims:
                press = rvar[press].values
                press = _conform(press, rvar.values.shape)
        elif isinstance(press, DataArray):
            press = press.values
        else:
            pass

    if dpd is not None:
        # DPD to RH
        vpdata = svp(temp.values - dpd.values, method=method, p=press, **kwargs)
        rvar.values = vpdata / svp(temp.values, method=method, p=press, **kwargs)
        origin = 't,dpd' if press is None else 't,dpd,p'
    else:
        vpdata = sh2vap(spec_humi.values, press)
        rvar.values = vpdata / svp(temp.values, method=method, p=press, **kwargs)
        origin = 't,q,p'

    r_att = {'units': '1', 'standard_name': 'relativ_humidity', 'long_name': 'relative humidity',
             'esat': method, 'origin': origin}

    if press is not None:
        r_att['enhancement_factor'] = "yes"
    with errstate(invalid='ignore'):
        if ((rvar.values < 0) | (rvar.values > 1)).any():
            message("Warning relative humidiy outside of normal range", **leveldown(**kwargs))

    r_att['precision'] = precision
    rvar.attrs.update(r_att)
    rvar.values = around(rvar.values, decimals=precision)
    rvar.name = 'rh'
    return rvar


def to_vp(temp, dpd=None, td=False, rel_humi=None, spec_humi=None, press=None, method='HylandWexler', precision=9,
          **kwargs):
    """ Convert dewpoint departure, rel. humidity or specific humidity to water vapor pressure

    1. DPD to VP (T)
    2. Td to VP
    3. RH to VP (T)
    4. Q to VP (T, P)

    Args:
        temp (DataArray): temperature (input)
        rel_humi (DataArray):   rel. humidity (input)
        dpd:   Name of DPD (Input)
        td (bool): if temp is dewpoint or temperature
        spec_humi:  Name of Q (Input)
        press:   Name of pressure (Input)
        method:  Saturation water vapor pressure formulation

    Returns
    -------
        DataArray : water vapor pressure [Pa]
    """
    from numpy import around
    from xarray import DataArray
    from .esat import svp
    from .humidity import sh2vap

    if not isinstance(temp, DataArray):
        raise ValueError("Requires a DataArray", type(temp))

    if dpd is None and spec_humi is None and rel_humi is None and not td:
        raise RuntimeError("Requires either dpd, td, q or r for conversion")

    if dpd is not None and not isinstance(dpd, DataArray):
        raise ValueError("DPD Requires a DataArray", type(dpd))

    if spec_humi is not None and not isinstance(spec_humi, DataArray):
        raise ValueError("Q Requires a DataArray", type(spec_humi))

    if rel_humi is not None and not isinstance(rel_humi, DataArray):
        raise ValueError("R Requires a DataArray", type(rel_humi))

    if spec_humi is not None:
        if press is None:
            raise RuntimeError("Conversion Q>VP requires a pressure variable as well")

    vpvar = temp.copy()
    for iatt in list(vpvar.attrs.keys()):
        del vpvar.attrs[iatt]

    if press is not None:
        if isinstance(press, str):
            if press in vpvar.dims:
                press = vpvar[press].values
                press = _conform(press, vpvar.values.shape)
        elif isinstance(press, DataArray):
            press = press.values
        else:
            pass

    if dpd is not None:
        # DPD to VP
        vpvar.values = svp(temp.values - dpd.values, method=method, p=press)
        origin = 't,dpd' if press is None else 't,dpd,p'
    elif td:
        # Td to VP
        vpvar.values = svp(temp.values, method=method, p=press)
        origin = 'td' if press is None else 'td,p'

    elif rel_humi is not None:
        # RH to VP
        vpvar.values = svp(temp.values, method=method, p=press) * rel_humi.values
        origin = 't,rh' if press is None else 't,rh,p'

    else:
        # Q to VP
        vpvar.values = sh2vap(spec_humi.values, press)
        origin = 'q,p'

    r_att = {'origin': origin, 'esat': method, 'standard_name': 'water_vapor_pressure',
             'long_name': 'water vapor pressure', 'units': 'Pa'}

    if press is not None:
        r_att['enhancement_factor'] = "yes"

    r_att['precision'] = precision
    vpvar.attrs.update(r_att)
    vpvar.values = around(vpvar.values, decimals=precision)
    vpvar.name = 'vp'
    return vpvar


def to_dpd(temp, rel_humi=None, vp=None, spec_humi=None, press=None, svp_method='HylandWexler',
           dewp_method='dewpoint_Boegel', precision=2, **kwargs):
    """ Convert relative humidity, specific humidity or water vapor pressure to dewpoint departure

    1. RH to DPD (T):
        + RH to VP (T)
        + VP to DPD (T)

    2. VP to DPD (T):
        + VP to DPD (T)

    3. Q to DPD (T, P):
        + Q to VP (P)
        + VP to DPD (T)

    Args:
        temp: Name of temperature (input)
        rel_humi: Name of rel. humidity (Input)
        vp: Name of DPD (Input)
        spec_humi: Name of Q (Input)
        press: Name of pressure (Input)
        svp_method: Saturation water vapor pressure formulation
        dewp_method:
        precision (int)

    Returns:
        DataArray: dewpoint depression [K]
    """
    from numpy import around, errstate, NaN
    from xarray import DataArray
    from .esat import svp
    from .humidity import sh2vap, dewpoint
    from ..fun import message, leveldown

    if not isinstance(temp, DataArray):
        raise ValueError("Requires a DataArray", type(temp))

    if rel_humi is None and spec_humi is None and vp is None:
        raise RuntimeError("Requires either r, q or vp for conversion")

    if rel_humi is not None and not isinstance(rel_humi, DataArray):
        raise ValueError("R Requires a DataArray", type(rel_humi))

    if spec_humi is not None and not isinstance(spec_humi, DataArray):
        raise ValueError("Q Requires a DataArray", type(spec_humi))

    if vp is not None and not isinstance(vp, DataArray):
        raise ValueError("VP Requires a DataArray", type(vp))

    if spec_humi is not None:
        if press is None:
            raise RuntimeError("Conversion Q>DPD requires a pressure variable as well")

    dpdvar = temp.copy()
    for iatt in list(dpdvar.attrs.keys()):
        del dpdvar.attrs[iatt]

    if press is not None:
        if isinstance(press, str):
            if press in dpdvar.dims:
                press = dpdvar[press].values
                press = _conform(press, dpdvar.values.shape)
        elif isinstance(press, DataArray):
            press = press.values
        else:
            pass

    kwargs['tol'] = kwargs.get('tol', 0.1)  # Dewpoint minimization accuracy
    if rel_humi is not None:
        # RH to VP to DPD
        # if 'ECMWF' in method:
        #     dpdvar.values = temp.values - dewpoint_ecmwf(temp.values, rel_humi.values)
        # else:
        vpdata = rel_humi.values * svp(temp.values, method=svp_method, p=press, **kwargs)
        dpdvar.values = temp.values - dewpoint(vpdata, method=dewp_method, **kwargs)
        origin = 't,rh' if press is None else 't,rh,p'
    elif vp is not None:
        # VP to DPD
        # if method.lower() == 'ecmwf':
        #     dpdvar.values = temp.values - vp2td(vp.values)
        # else:
        dpdvar.values = temp.values - dewpoint(vp.values, method=dewp_method, **kwargs)
        origin = 't,vp'
    else:
        # Q to DPD
        # if method.lower() == 'ecmwf':
        #     dpdvar.values = temp.values - rh2td(temp.values, q2rh(spec_humi.values, temp.values, press.values))
        # else:
        try:
            vpdata = sh2vap(spec_humi.values, press.values)
        except:
            vpdata = sh2vap(spec_humi.values, press)
        dpdvar.values = temp.values - dewpoint(vpdata, method=dewp_method, **kwargs)
        origin = 't,q,p'

    r_att = {'origin': origin, 'svp': svp_method, 'dewp': dewp_method, 'standard_name': 'dew_point_depression',
             'long_name': 'dew point depression', 'units': 'K'}
    if press is not None:
        r_att['enhancement_factor'] = "yes"

    with errstate(invalid='ignore'):
        if (dpdvar.values < 0).any():
                # message("dew point depression outside range",dpdvar.values, **leveldown(**kwargs))
            dpdvar.values[dpdvar.values < -0.1] = NaN
            dpdvar.values[dpdvar.values < 0.0] = 0.0

    r_att['precision'] = precision
    dpdvar.values = around(dpdvar.values, decimals=precision)
    dpdvar.name = 'dpd'
    dpdvar.attrs.update(r_att)
    return dpdvar


def to_sh(vp=None, temp=None, rel_humi=None, dpd=None, press=None, method='HylandWexler', precision=6, **kwargs):
    """
    vp -(p)  sh
    rh (t,p) -> vp (p)-> sh
    dpd (t,p) -> vp (p)-> sh

    Args:
        vp: water vapor pressure
        temp: temperature
        rel_humi: rel. humidity
        dpd: dewpoint departure
        press: air pressure
        method:
        precision:

    Returns:
        spec_humi: specific humidity
    """
    from numpy import around
    from xarray import DataArray
    from .esat import svp
    from .humidity import vap2sh

    if not isinstance(temp, DataArray):
        raise ValueError("Requires a DataArray", type(temp))

    if rel_humi is None and dpd is None and vp is None:
        raise RuntimeError("Requires either r, dpd or vp for conversion")

    if rel_humi is not None and not isinstance(rel_humi, DataArray):
        raise ValueError("R Requires a DataArray", type(rel_humi))

    if dpd is not None and not isinstance(dpd, DataArray):
        raise ValueError("DPD Requires a DataArray", type(dpd))

    if vp is not None and not isinstance(vp, DataArray):
        raise ValueError("VP Requires a DataArray", type(vp))

    if press is None:
        raise RuntimeError("Conversion ?>Q requires a pressure variable as well")

    if vp is not None:
        qvar = vp.copy()
    else:
        qvar = temp.copy()

    if isinstance(press, str):
        if press in qvar.dims:
            press = qvar[press].values
            press = _conform(press, qvar.values.shape)
    elif isinstance(press, DataArray):
        press = press.values
    else:
        pass

    kwargs['tol'] = kwargs.get('tol', 0.1)  # Dewpoint minimization accuracy
    if vp is not None:
        # VP to Q
        qvar.values = vap2sh(vp.values, press)
        origin = 'vp,p'
    elif rel_humi is not None:
        # RH to Q
        vpdata = rel_humi.values * svp(temp.values, method=method, p=press, **kwargs)
        qvar.values = vap2sh(vpdata, press)
        origin = 't,rh,p'
    else:
        # DPD to Q
        vpdata = svp(temp.values - dpd.values, method=method, p=press, **kwargs)
        qvar.values = vap2sh(vpdata, press)
        origin = 't,dpd,p'

    r_att = {'origin': origin, 'esat': method, 'standard_name': 'specific_humidity',
             'long_name': 'specific humidity', 'units': 'kg/kg'}
    if press is not None:
        r_att['enhancement_factor'] = "yes"

    r_att['precision'] = precision
    qvar.attrs.update(r_att)
    qvar.values = around(qvar.values, decimals=precision)
    qvar.name = 'sh'
    return qvar


def to_dewpoint(vp=None, temp=None, rel_humi=None, spec_humi=None, press=None, svp_method='HylandWexler',
                dewp_method='dewpoint_Boegel', precision=2, **kwargs):
    """ calculate dewpoint
    VP -> Td
    T, RH (p= -> vp -> Td
    Q,P -> vp -> Td

    Returns:
        dewp: dewpoint
    """
    from numpy import around
    from xarray import DataArray
    from .esat import svp
    from .humidity import sh2vap, dewpoint
    from ..fun import message, leveldown

    if rel_humi is None and temp is None and vp is None and spec_humi is None:
        raise RuntimeError("Requires either rel. humidity, temp or vp for conversion")

    if rel_humi is not None and not isinstance(rel_humi, DataArray):
        raise ValueError("RH Requires a DataArray", type(rel_humi))

    if temp is not None and not isinstance(temp, DataArray):
        raise ValueError("TEMP Requires a DataArray", type(temp))

    if rel_humi is not None and temp is None:
        raise ValueError("requires TEMP and RH")

    if vp is not None and not isinstance(vp, DataArray):
        raise ValueError("VP Requires a DataArray", type(vp))

    if spec_humi is not None:
        if press is None:
            raise RuntimeError("Conversion Q>Td requires a pressure variable as well")

        dewp = spec_humi.copy()

    elif temp is not None:
        dewp = temp.copy()

    else:
        dewp = vp.copy()

    if press is not None:
        if isinstance(press, str):
            if press in dewp.dims:
                press = dewp[press].values
                press = _conform(press, dewp.values.shape)
        elif isinstance(press, DataArray):
            press = press.values
        else:
            pass

    if vp is not None:
        dewp.values = dewpoint(vp.values, method=dewp_method, **kwargs)
        origin = 'vp'
    elif temp is not None:
        vp = rel_humi.values * svp(temp.values, method=svp_method, p=press, **kwargs)
        dewp.values = dewpoint(vp, method=dewp_method, **kwargs)

        if ((temp.values - dewp.values) < 0).any():
            message("Dew point supersaturated", **leveldown(**kwargs))

        origin = 't,rh' if press is None else 't,rh,p'
    else:
        dewp.values = sh2vap(spec_humi.values, press)
        origin = 'q,p'

    r_att = {'origin': origin, 'svp': svp_method, 'dewp': dewp_method, 'standard_name': 'dew_point',
             'long_name': 'dew point', 'units': 'K'}
    if press is None:
        r_att['enhancement_factor'] = "yes"

    r_att['precision'] = precision
    dewp.values = around(dewp.values, decimals=precision)
    dewp.name = 'td'
    dewp.attrs.update(r_att)
    return dewp


def _conform(data, shape):
    """ Make numpy array conform to a certain shape

    Args:
        data:
        shape:

    Returns:

    """
    import numpy as np
    if not isinstance(data, np.ndarray):
        raise ValueError('Requires a numpy array')

    if not isinstance(shape, (tuple, list)):
        raise ValueError('Requires a tuple or list')

    data = data.copy()
    n = data.shape

    assert np.any([i in shape for i in n]), "Shapes do not allign?!"

    for i, j in enumerate(shape):
        if j not in n:
            data = np.expand_dims(data, axis=i)
    return data
