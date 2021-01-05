

def calculate_observation_errors(obs, analysis, firstguess, dim='time', n=730, min_periods=20, center=True, **kwargs):
    """ Calculate Observation Error statistics from first guess and analysis departures

    Args:
        obs (xarray.DataArray): Observations
        analysis (xarray.DataArray): Analysis
        firstguess (DataArray): First Guess
        dim (str): datetime dimension
        n (int): number of average (730 = 2 Years)
        min_periods (int): minimum required values
        center (bool): average alignment
        **kwargs:

    Returns:

    """
    import xarray
    #
    # Estimate obs errors from fg, an dep.
    #
    an_dep = obs - analysis
    an_dep.name = 'an_dep'
    fg_dep = obs - firstguess
    fg_dep.name = 'fg_dep'
    tmp = xarray.merge([an_dep, fg_dep])
    tmp = tmp.rolling(min_periods=min_periods, center=center, **{dim: n}).construct('pieces')
    return obserror(tmp, a='an_dep', b='fg_dep', dim='pieces', attrs=obs.attrs)


def obserror(data, a='t_an_dep', b='t_fg_dep', dim='time', attrs=None, **kwargs):
    import xarray as xr
    import numpy as np
    #
    #
    #
    if not isinstance(data, xr.Dataset):
        raise ValueError()
    #
    # figure out dimensions
    #
    jdims = list(data[a].dims)
    jdims.remove(dim)
    axis = list(data[a].dims).index(dim)
    if attrs is None:
        attrs = data[a].attrs.copy()
    #
    # Calculate product of DataArray's
    #
    tmp = np.sqrt(xr.apply_ufunc(np.nanmean, data[a]*data[b],
                                 input_core_dims=[data[a].dims],
                                 output_core_dims=[jdims],
                                 kwargs={'axis': axis}
                                 )
                  )
    #
    # Metadata
    #
    tmp.attrs.update(attrs)
    tmp.attrs.update({'standard_name': 'obs_error'})
    return tmp


def errorous(x, dim='time', n=730, min_periods=20, center=True, sigma=6, **kwargs):
    """ removes errorous values outside 4 sigma

    Args:
        x:
        dim:
        n:
        min_periods:
        center:

    Returns:

    """
    import xarray as xr
    import numpy as np
    if not isinstance(x, xr.DataArray):
        raise ValueError()

    m = x.rolling(min_periods=min_periods, center=center, **{dim: n}).mean().ffill(dim).bfill(dim)
    s = x.rolling(min_periods=min_periods, center=center, **{dim: n}).std().ffill(dim).bfill(dim)
    # print(m,s)
    x = x.where((~s.isnull()) & (~m.isnull()) & (x < (m + sigma*s)) & (x > (m - sigma*s)), np.nan)
    return x


def obserror_dpd(data, a='t_an_dep', b='t_fg_dep', c='rh_an_dep', d='rh_fg_dep', dim='time', attrs=None, **kwargs):
    import xarray as xr
    import numpy as np
    from .convert import to_dpd
    #
    #
    #
    if not isinstance(data, xr.Dataset):
        raise ValueError()
    #
    # figure out dimensions
    #
    jdims = list(data[a].dims)
    jdims.remove(dim)
    axis = list(data[a].dims).index(dim)
    if attrs is None:
        attrs = data[a].attrs.copy()
    #
    # Calculate DPD dep
    # sqrt( t_err^2 + td_err^2)

    # Error DPD ?
    # dpd = t - td(esat(T)*rh)
    # error prop for esat
    # e_err = 2.57775e6*sqrt((exp((35.004*t-9561.69)/(t-32.19))*t_1**2 )/ (t-32.19)**4)
    # sqrt( t_err^2 +   e_err^2 + (rh_err / abs(rh))^2 + e_err^2 * 216.218^2)
    #
    # Calculate product of DataArray's
    #
    tmp = np.sqrt(xr.apply_ufunc(np.nanmean, data[a]*data[b],
                                 input_core_dims=[data[a].dims],
                                 output_core_dims=[jdims],
                                 kwargs={'axis': axis}
                                 )
                  )
    #
    # Metadata
    #
    tmp.attrs.update(attrs)
    tmp.attrs.update({'standard_name': 'obs_error'})
    return tmp
