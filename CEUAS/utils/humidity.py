import numpy as np

"""
Three formulas for saturation water vapor pressure (e_sat)
Two dewpoint formulas

for ECMWF (ODBs) use dewpoint_ECMWF and FOEEWMO for e_sat

for IGRA and other sources use:
dewpoint_Bolton and either Boegel or Bolton for e_sat
"""


def dewpoint_ECMWF(e, **kwargs):
    """ ECMWF IFS CYC31R1 Data Assimilation Documentation (Page 86-88)
    Derived Variables from TEMP
    Temperature and dew point are transformed into realtive humidity (RH) for
    TEMP observations, with a further transformation of RH into specific humidity (Q)
    for TEMP observations.

    Uses foeewmo for water only saturation water vapor pressure

    Returns
    -------
    td      dew point       [K]
    """
    lnpart = np.log(e / 611.21)
    return (17.502 * 273.16 - 32.19 * lnpart) / (17.502 - lnpart)


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



def Boegel(t, over_water=True, over_ice=False, **kwargs):
    """ Boegel 1979
    Magnus formula used by Tetens and converted by Murray, and modified by Boegel, found in Buck, 1981

    Used by (mistake in formula):
    McCarthy, M. P., Titchner, H. A., Thorne, P. W., Tett, S. F. B., Haimberger, L. and
    Parker, D. E.: Assessing Bias and Uncertainty in the HadAT-Adjusted Radiosonde Climate Record,
    Journal of Climate, 21(4), 817-832, doi:10.1175/2007JCLI1733.1, 2008.

    Uses enhancement factor, according to McCarthy, 2008

    Args:
        t: temperatur in K

    Returns:
        es: saturation water vapor pressure in Pa
    """
    t = t - 273.15  # requires Degrees

    # wrong from McCarthy
    # return 611.21 * f * np.exp((18.729 - t / 227.3) * t / 257.87)
    if over_water:
        return 611.21 * np.exp((18.729 - t / 227.3) * t / (t + 257.87))
    elif over_ice:
        return 611.15 * np.exp((23.036 - t / 333.7) * t / (t + 279.82))
    else:
        return np.where(t < 0.01,
                        611.15 * np.exp((23.036 - t / 333.7) * t / (t + 279.82)),
                        611.21 * np.exp((18.729 - t / 227.3) * t / (t + 257.87))
                        )


def Bolton(t, **kwargs):
    """Bolton 1980

    tt = 273.16
    ew = 611.12 * math.exp(17.67*(t-273.15)/(t-tt+243.5) )

    Args:
        t: temperatur in K

    Returns:
        es: saturation water vapor pressure in Pa
    """
    tt = 273.16
    return 611.12 * np.exp(17.67 * (t - 273.15) / (t - tt + 243.5))



def FOEEWMO(t, **kwargs):
    """from IFS Documentation Cycle 31,
    Teten's formula for water only
    after Tiedtke 1993 (cloud and large-scale precipitation scheme)
    Uses the Tetens's formula
    Based on Buck 1981 & Alduchov and Eskridge 1996

    Args:
        t: air temperature K

    Returns:
        es : saturation water vapor in Pa
    """
    return 611.21 * np.exp(17.502 * (t - 273.16) / (t - 32.19))