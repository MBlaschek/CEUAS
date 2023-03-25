import numpy as np
import pandas as pd

__all__ = ['svp', 'comparison_chart', 'reference_table', 'Boegel', 'Bolton', 'Buck', 'Goff1957',
           'GoffGratch', 'FOEEWMO', 'HylandWexler', 'IAPWS', 'MurphyKoop', 'Sonntag', 'Wright']


def svp(t, method='HylandWexler', p=None, **kwargs):
    """
    Saturation water vapor pressure from Temperature
    The equations by Hyland and Wexler [4], the nearly identical equation by Wexler (1976, see reference below) and
    the equation by Sonntag [7] are the most commonly used equations among Radiosonde manufacturers
    and should be used in upper air applications to avoid inconsistencies.

    Known Saturation water Vapor Formulations:

    Bolton 1980
    Goff 1957, 1965        (180    - 273.15 / 273.15 - 373.15) 1957> WMO
    Goff and Gratch 1946   (184    - 273.15 / 273.15 - 373.15)
    Hyland and Wexler 1983 (173.16 - 273.15 / 273.15 - 473.15) Vaisala
    IAPWS 1995             (173.15 - 273.15 / 273.15 - 647   ) Wagner and Pruss + Wagner 1994
    Murphy and Koop 2005   (110    - 273.15 / 132      332   ) Ice formulation as well
    Sonntag 1990           (173.15 - 273.15 / 273.15 - 373.15)
    Sonntag 1994           (173.15 - 273.15 / 273.15 - 373.15)
    Wagner 1994            (190    - 273.15                  )
    Wagner and Pruss 1993  (190    - 273.15 / 273.15 - 647   )
    Wexler 1976            (                  273.15 - 373.15)
    Wright 1997

    Args:
        t: air temperature
        method: string or function
        p: pressure
        **kwargs: additional keywords passed to function

    Returns:
        es : saturation water vapor pressure in Pa
    """
    try:
        if callable(method):
            vpfunc = method
        else:
            vpfunc = eval(method)

        if p is not None:
            f = ef(p, **kwargs)
        else:
            f = 1.
        return vpfunc(t, **kwargs) * f
    except:
        import sys
        print("Functions: ", ", ".join([i for i in dir(sys.modules[__name__]) if i[0].upper() == i[0]]))


def reference_table():
    """Table 1 from Murphy and Koop (2005)
    Values of Saturation Water Vapor Pressure over liquid and ice.
    """
    print(""" This dataset has been calculated by Murphy and Koop (2005)
    and can be used as reference dataset for saturation water vapor pressure.
    IAPWS 97 Values give:
     T (K)   Es (Pa)
    213.15  1.0813475449
    243.15  38.005139487
    273.15  611.15347506
    303.15  4246.6883405
    333.15  19945.801925
    """)
    tx = np.asarray([-60, -30, 0, 30, 60]) + 273.15
    ex = np.asarray([1.0813475449, 38.005139487, 611.15347506, 4246.6883405, 19945.801925])
    t = np.asarray([150, 180, 210, 240, 273.15, 273.16, 300])
    ei = np.asarray([6.106e-6, 0.0053975, 0.70202, 27.272, 611.154, 611.657, np.nan])
    ew = np.asarray([1.562e-5, 0.011239, 1.2335, 37.667, 611.213, 611.657, 3536.8])
    return pd.DataFrame({'t': t, 'ei': ei, 'ew': ew, 'mix': np.where(t < 273.16, ei, ew)})


def comparison_chart(methods=None, over_water=True, over_ice=False, p=None, **kwargs):
    """ Make a comparison between Esat functions

    Returns
    -------

    """
    data = reference_table()
    if over_water:
        anomaly = 'ew'
    elif over_ice:
        anomaly = 'ei'
    else:
        anomaly = 'mix'

    print("Using %s as reference" % anomaly)
    if methods is None:
        methods = ['MurphyKoop']

    for im in methods:
        try:
            data[im] = svp(data['t'].values, method=im, p=p, **kwargs)
            data[im] -= data[anomaly]

        except Exception as e:
            print("Error: ", im, repr(e))
            data[im] = np.nan
            continue
    return data


def ef(p, t=None, over_water=True, over_ice=False, **kwargs):
    """
    from
    Sugidachi, T. and Fujiwara, M.: Correction of the Stepwise Change Observed at 0C in
    Meisei RS2-91, RS-01G, and RS-06G Radiosonde Relative Humidity Profiles,
    Journal of the Meteorological Society of Japan. Ser. II, 91(3), 323-336,
    doi:10.2151/jmsj.2013-306, 2013.

    Args:
        t: air temperature K
        p: air pressure Pa

    Returns:
        f : enhancement factor for saturation water vapor pressure
            depending on air temperature and air pressure
    """
    # Vaisala / Buck 1981  / WMO 2008
    # 1.0016 + 3.15e-6*p - 0.074 / p  # kPa
    if over_water:
        return 1.0007 + (3.46e-6 * p / 100)
    if over_ice:
        return 1.0003 + (4.18e-6 * p / 100)
    if t is not None:
        return np.where(t < 273.16,
                        1.0003 + (4.18e-6 * p / 100),
                        1.0007 + (3.46e-6 * p / 100))


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


def Buck(t, over_ice=False, over_water=True, **kwargs):
    """Buck 1981 Formulation
    pressure enhancement term

    For general use, equation ew and ei are recommended.
    Other equations are available as well.
    These equations correspond to ew1, ei2, fw3 and fi3 of the origional publication.

    ew = (1.000007 + (3.46*10**(-6)*p)) *
         611.21*math.exp( 17.502*(t-273.15) / (240.97+t-273.15) )

    ei = (1.000003 + (4.18*10**(-6)*p)) *
         611.15*math.exp( 22.452*(t-273.15) / (272.55+t-273.15) )

    Args:
        t: air temperaure in K
        p: air pressure in Pa

    Returns:
        es: saturation water vapor pressure in Pa
    """
    t = t - 273.15  # needs Celcius

    if over_water:
        return 611.21 * np.exp(17.502 * t / (240.97 + t))
    elif over_ice:
        return 611.15 * np.exp(22.452 * t / (272.55 + t))
    else:
        return np.where(t < 0.01,
                        611.15 * np.exp(22.452 * t / (272.55 + t)),
                        611.21 * np.exp(17.502 * t / (240.97 + t))
                        )


def Buck1996(t, over_ice=False, over_water=True, **kwargs):
    """Buck 1996 Formulation
    modified version of Buck 1981

    Args:
        t: air temperature K

    Returns:
        es : saturation water vapor in Pa
    """
    t = t - 273.15  # needs Celsius
    if over_water:
        return 611.21 * np.exp((18.678 - t / 234.5) * (t / (257.14 + t)))
    elif over_ice:
        return 611.15 * np.exp((23.036 - t / 333.7) * (t / (279.82 + t)))
    else:
        return np.where(t < 0.01,
                        611.15 * np.exp((23.036 - t / 333.7) * (t / (279.82 + t))),
                        611.21 * np.exp((18.678 - t / 234.5) * (t / (257.14 + t)))
                        )


def Goff1957(temp, over_ice=False, over_water=True, **kwargs):
    """Goff (1957): Units are converted here from atmospheres to Pa. Stated ranges 180 < T
    < 273.16 K for ice and 273.15 < T < 373.15 K for liquid with extension to 223 K.

    log (ew)  = math.log10(611.14) + 10.79574*(1 -tt/t )
                 - 5.0280*math.log10(t /tt)
                 + 1.50475 *10**(-4)*(1 -10*(-8.2969*(t/tt -1)))
                 + 0.42873 *10**(-3)*(10**(4.76955*(1-tt/t)) - 1)

    log (ei) = math.log10(611.14) - 9.096853*(tt/t - 1)
                 - 3.566506* math.log10(tt/t)
                 + 0.876812*(1 - t/tt)

    Args:
        temp: air temperature in K
        over_ice: use only water vapor over ice
        kwargs: dummy

    Returns:
         es : saturation water vapor pressure in Pa
    """

    def liquid(t):
        tt = 273.16
        return np.log10(611.14) + 10.79574 * (1 - tt / t) \
               - 5.0280 * np.log10(t / tt) \
               + 1.50475e-4 * (1 - 10 * (-8.2969 * (t / tt - 1))) \
               + 0.42873e-3 * (10 ** (4.76955 * (1 - tt / t)) - 1)

    def ice(t):
        tt = 273.16
        return np.log10(611.14) - 9.096853 * (tt / t - 1) \
               - 3.566506 * np.log10(tt / t) \
               + 0.876812 * (1 - t / tt)

    if over_water:
        return 10 ** liquid(temp)
    elif over_ice:
        return 10 ** ice(temp)
    else:
        return np.where(temp < 273.16, 10 ** ice(temp), 10 ** liquid(temp))


def Goff1965(temp, over_water=True, over_ice=False, **kwargs):
    """Goff (1965): Units are converted here from atmospheres to Pa.
    The range is the same as Goff (1957). This 1965 publication, taken from the proceedings of a conference in 1963, has also been frequently referenced as Goff 1963. It is a minor correction to Goff (1957).


    log( ew ) = math.log10(611.11) + 10.79586 *(1 - tt/t )
                - 5.02808*math.log10(t/tt)
                + 1.50474*10**(-4)*(1-10**(-8.29692*(t/tt -1)))
                + 0.42873*10**(-3)*(10**(4.76955*(1- tt/t)) - 1)

    log( ei ) = math.log10(611.11) - 9.096936*(tt/t - 1)
                - 3.56654*math.log10(tt/t )
                + 0.876817*(1- t/tt)

    Args:
        temp: air temperature in K
        liquid_only: use only water vapor over liquid water
        ice_only: use only water vapor over ice
        kwargs: dummy

    Returns:
         es : saturation water vapor pressure in Pa
    """

    def liquid(t):
        tt = 273.16
        return np.log10(611.11) + 10.79586 * (1 - tt / t) \
               - 5.02808 * np.log10(t / tt) \
               + 1.50474 * 10 ** (-4) * (1 - 10 ** (-8.29692 * (t / tt - 1))) \
               + 0.42873 * 10 ** (-3) * (10 ** (4.76955 * (1 - tt / t)) - 1)

    def ice(t):
        tt = 273.16
        return np.log10(611.11) - 9.096936 * (tt / t - 1) \
               - 3.56654 * np.log10(tt / t) \
               + 0.876817 * (1 - t / tt)

    if over_water:
        return 10 ** liquid(temp)
    elif over_ice:
        return 10 ** ice(temp)
    else:
        return np.where(temp < 273.16, 10 ** ice(temp), 10 ** liquid(temp))


def GoffGratch(temp, over_water=True, over_ice=False, **kwargs):
    """Goff and Gratch (1946): Units are converted here from atmospheres to Pa.
    Stated ranges: 184<T <273.16K for ice and 273.15 < T <373.15K for liquid.

    log( ew ) = -7.90298*( (373.16/t ) - 1) + 5.02808* math.log10(373.16/t)
                - 1.3816 *10**(-7)*(10**(11.344*(1-t/373.16)) - 1)
                + 8.1328*10**(-3)*(10**(-3.49149*(373.16/t -1)) - 1)
                + math.log10(101325)

    log( ei ) = -9.09718*((tt/t) - 1) - 3.56654*math.log10(tt/t)
                + 0.876793*(1 - (t /tt))
                + math.log10(610.71)

    Args:
        temp: air temperature in K
        liquid_only: use only water vapor over liquid water
        ice_only: use only water vapor over ice
        kwargs: dummy

    Returns:
         es : saturation water vapor pressure in Pa
    """

    def liquid(t):
        return -7.90298 * ((373.16 / t) - 1) \
               + 5.02808 * np.log10(373.16 / t) \
               - 1.3816 * 10 ** (-7) * (10 ** (11.344 * (1 - t / 373.16)) - 1) \
               + 8.1328 * 10 ** (-3) * (10 ** (-3.49149 * (373.16 / t - 1)) - 1) \
               + np.log10(101325)

    def ice(t):
        tt = 273.16
        return -9.09718 * ((tt / t) - 1) - 3.56654 * np.log10(tt / t) \
               + 0.876793 * (1 - (t / tt)) + np.log10(610.71)

    if over_water:
        return 10 ** liquid(temp)
    elif over_ice:
        return 10 ** ice(temp)
    else:
        return np.where(temp < 273.16, 10 ** ice(temp), 10 ** liquid(temp))


def FOEEWM(t, **kwargs):
    """from IFS Documentation Cycle 31,
    Teten's formula for mixed phases water, ice
    after Tiedtke 1993 (cloud and large-scale precipitation scheme)
    Based on Buck 1981 & Alduchov and Eskridge 1996

    Args:
        t: air temperature K

    Returns:
        es : saturation water vapor in Pa
    """
    # T Larger than 273.15 K > only water
    ew = 611.21 * np.exp(17.502 * (t - 273.16) / (t - 32.19))  # Liquid
    ei = 611.21 * np.exp(22.587 * (t - 273.16) / (t + 0.7))  # Ice
    e = np.where(t > 273.15, ew, ei)
    # combine ice and water for mixed clouds
    e = e + np.where((t >= 250.15) & (t <= 273.15), (ew - ei) * ((t - 250.16) / 23.) ** 2, 0)
    return e


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


def HylandWexler(temp, over_water=True, over_ice=False, **kwargs):
    """Hyland and Wexler (1983), also in Wexler and Hyland (1983): Stated ranges 173.16 ?
    T < 273.16 for ice and 273.16 > 473.15 for liquid.

    Used by Vaisala

    ln( ew ) = -5800.2206/t + 1.3914993
               - 0.48640239 * 10**(-1)*t
               + 0.41764768 * 10**(-4)*t**2
               - 0.14452093 * 10**(-7)*t**3
               + 6.5459673*math.log(t)

    ln( ei ) = -5674.5359/t + 6.3925247
               - 0.96778430 * 10**(-2)*t
               + 0.62215701 * 10**(-6)*t**2
               + 0.20747825 * 10**(-8)*t**3
               - 0.94840240 * 10**(-12)*t**4
               + 4.1635019*math.log(t)

    Args:
        temp: air temperature in K
        liquid_only: use only water vapor over liquid water
        ice_only: use only water vapor over ice
        kwargs: dummy

    Returns:
         es : saturation water vapor pressure in Pa
    """

    def liquid(t):
        mask = ~np.isnan(t)
        if np.sum(mask) > 0:
            
            tmask = t[mask]
            res = t.copy()
            res[mask] = np.exp(
                -5800.2206 / tmask + 1.3914993
                - 0.48640239e-1 *tmask
                + 0.41764768e-4 * tmask * tmask
                - 0.14452093e-7 * tmask * tmask * tmask
                + 6.5459673 * np.log(tmask)
            )
            return res
        else:
            return t

    def ice(t):
        return np.exp(
            -5674.5359 / t + 6.3925247
            - 0.96778430e-2 * t
            + 0.62215701e-6 * t * t
            + 0.20747825e-8 * t * t * t
            - 0.94840240e-12 * t * t * t * t
            + 4.1635019 * np.log(t)
        )

    if over_water:
        return liquid(temp)
    elif over_ice:
        return ice(temp)
    else:
        return np.where(temp < 273.16, ice(temp), liquid(temp))


def IAPWS(t, over_water=True, over_ice=False, **kwargs):
    """Formulation 1995 of IAPWS
    (The International Association for the Properties of Water and Steam)

    Liquid:
        Wagner and Pruss, 1993, 2002
    Ice:
        Wagner 1994

    Args:
        t: air temperature K

    Returns:
        es : saturation water vapor in Pa
    """
    if over_water:
        return WagnerPruss(t)
    elif over_ice:
        return Wagner1994(t)
    else:
        return np.where(t < 273.16, Wagner1994(t), WagnerPruss(t))


def MurphyKoop(temp, over_water=True, over_ice=False, **kwargs):
    """Murphy and Koop, 2005
    ln (ew ) = 54.842763 - 6763.22/t
              - 4.210 * math.log(t) + 0.000367*t
              + math.tanh( 0.0415 * (t - 218.8) )
              * (53.878 - 1331.22/t - 9.44523 * math.log(t) + 0.014025 * t)
    for 123 < T < 332 K.

    Args:
        temp: air temperature in K
        liquid_only: use only water vapor over liquid water
        ice_only: use only water vapor over ice
        kwargs: dummy

    Returns:
         es : saturation water vapor pressure in Pa
    """

    def liquid(t):
        return np.exp(
            54.842763 - 6763.22 / t
            - 4.210 * np.log(t)
            + 0.000367 * t
            + np.tanh(0.0415 * (t - 218.8))
            * (53.878 - 1331.22 / t - 9.44523 * np.log(t) + 0.014025 * t))

    def ice(t):
        return np.exp(9.550426 - 5723.265 / t + 3.53068 * np.log(t) - 0.00728332 * t)

    if over_water:
        return liquid(temp)
    elif over_ice:
        return ice(temp)
    else:
        return np.where(temp < 273.16, ice(temp), liquid(temp))


def Tetens(t, over_water=True, over_ice=False, **kwargs):
    """Tetens 1967 after Magnus in Murray, 1967

    Args:
        t: temperatur in K

    Returns:
        es: saturation water vapor pressure in Pa
    """
    if over_water:
        return 610.78 * np.exp(17.269 * (t - 273.16) / (t - 35.86))
    elif over_ice:
        return 610.78 * np.exp(21.875 * (t - 273.16) / (t - 7.66))
    else:
        return np.where(t > 273.15,
                        610.78 * np.exp(17.269 * (t - 273.16) / (t - 35.86)),
                        610.78 * np.exp(21.875 * (t - 273.16) / (t - 7.66)))


def Sonntag(temp, over_water=True, over_ice=False, **kwargs):
    """Sonntag (1994): Sonntag, D., Advancements in the field of hygrometry,
    Meteorologische Zeitschrift, 3, 51-66, 1994. in Pa.
    Liquid: 273.16 < T < 373.15, Ice: 173.15 < T < 273.15

    ln ( ew )  = -6096.9385/t + 21.2409642
                 - 2.711193*10**(-2)*t
                 + 1.673952*10**(-5)*t**2
                 + 2.433502*math.log(t)

    ln ( ei )  = -6024.5282/t + 29.32707
                 + 1.0613868*10**(-2)*t
                 - 1.3198825*10**(-5)*t**2
                 - 0.49382577*math.log(t)

    Args:
        temp: air temperature in K
        liquid_only: use only water vapor over liquid water
        ice_only: use only water vapor over ice
        kwargs: dummy

    Returns:
         es : saturation water vapor pressure in Pa
    """

    def liquid(t):
        return -6096.9385 / t + 21.2409642 \
               - 2.711193e-2 * t \
               + 1.673952e-5 * t * t \
               + 2.433502 * np.log(t)

    def ice(t):
        return -6024.5282 / t + 29.32707 \
               + 1.0613868e-2 * t \
               - 1.3198825e-5 * t * t \
               - 0.49382577 * np.log(t)

    if over_water:
        return np.exp(liquid(temp))
    elif over_ice:
        return np.exp(ice(temp))
    else:
        return np.where(temp < 273.16, np.exp(ice(temp)), np.exp(liquid(temp)))


def Sonntag1990(temp, over_water=True, over_ice=False, **kwargs):
    """Sonntag (1990)Sonntag (1990): Stated ranges are 173.15 < T < 273.16 K for ice and 173.15 < T <
    373.15 K for liquid, despite being based on Wexler (1976).

    ln ( ew )/100 = -6096.9385/t + 16.635794
                    - 2.711193*10**(-2)*t
                    + 1.673952*10**(-5)*t**2
                    + 2.433502*math.log(t)

    ln ( ei )/100 = 24.7219 - 6024.5282/t
                    + 1.0613868 * 10**(-2)*t
                    - 1.3198825 * 10**(-5)*t**2
                    - 0.49382577*math.log(t)

    Args:
        temp: air temperature in K
        liquid_only: use only water vapor over liquid water
        ice_only: use only water vapor over ice
        kwargs: dummy

    Returns:
         es : saturation water vapor pressure in Pa
    """

    def liquid(t):
        return -6096.9385 / t + 16.635794 \
               - 2.711193e-2 * t \
               + 1.673952e-5 * t * t \
               + 2.433502 * np.log(t)

    def ice(t):
        return 24.7219 - 6024.5282 / t \
               + 1.0613868e-2 * t \
               - 1.3198825e-5 * t * 2 \
               - 0.49382577 * np.log(t)

    if over_water:
        return 100 * np.exp(liquid(temp))
    elif over_ice:
        return 100 * np.exp(ice(temp))
    else:
        return np.where(temp < 273.16, 100 * np.exp(ice(temp)), 100 * np.exp(liquid(temp)))


def Wagner1994(t, **kwargs):
    """Wagner et al. (1994):

    ln( ei ) = math.log(611.657) - 13.9281690*(1-(tt/t)**1.5)
               + 34.7078238*(1-(tt/t)**1.25)

    190 < T < 273.16 K

    Args:
        t: temperatur in K

    Returns:
        es: saturation water vapor pressure in Pa
    """
    tt = 273.16
    return np.log(611.657) - 13.9281690 * (1 - (tt / t) ** 1.5) + 34.7078238 * (1 - (tt / t) ** 1.25)


def WagnerPruss(t, **kwargs):
    """Wagner and Pruss (1993):
    Updated in Wagner and Pruss (2002)
    IAPWS Formulation

    tc = 647.096 # K
    tau = 1 - t/tc
    ln( ew / 2.2064*10**7 ) = (tc/t)*( -7.85951783*tau
                                       + 1.84408259*tau**(1.5)
                                       - 11.7866497*tau**3
                                       + 22.6807411*tau**(3.5)
                                       - 15.9618719*tau**4
                                       + 1.80122502*tau**(7.5) )

    273.16 < T < 647 K

    Args:
        t: temperatur in K

    Returns:
        es: saturation water vapor pressure in Pa
    """
    tc = 647.096  # K
    tau = 1 - t / tc
    ew = (tc / t) * (-7.85951783 * tau + 1.84408259 * tau ** (1.5) -
                     11.7866497 * tau * tau * tau + 22.6807411 * tau ** (3.5) -
                     15.9618719 * tau * tau * tau * tau + 1.80122502 * tau ** (7.5))

    return np.exp(ew) * 2.2064e7


def Wexler(temp, over_water=False, over_ice=False, **kwargs):
    """Wexler (1976)

    ew = -2991.27290*t**(-2) - 6017.0128*t**(-1)
         + 18.87643854 - 0.028354721*t
         + 0.17838301*10**(-4)*t**2
         - 0.84150417*10**(-9)*t**3
         + 0.44412543*10**(-12)*t**4
         + 2.858487*math.log(t)

    ei = -5865.3696*t**(-1) + 22.241033
         + 0.013749042*t - 0.34031775*10**(-4)*t**2
         + 0.26967687*10**(-7)*t**3
         + 0.6918651 * math.log(t)

    Args:
        t: temperatur in K
        liquid_only: use only water vapor over liquid water
        ice_only: use only water vapor over ice
        kwargs: dummy

    Returns:
        es: saturation water vapor pressure in Pa
    """

    def liquid(t):
        return -2991.27290 / t / t - 6017.0128 / t \
               + 18.87643854 - 0.028354721 * t \
               + 0.17838301e-4 * t * t \
               - 0.84150417e-9 * t * t * t \
               + 0.44412543e-12 * t * t * t * t \
               + 2.858487 * np.log(t)

    def ice(t):
        return -5865.3696 / t + 22.241033 \
               + 0.013749042 * t \
               - 0.34031775e-4 * t * t \
               + 0.26967687e-7 * t * t * t \
               + 0.6918651 * np.log(t)

    if over_water:
        return np.exp(liquid(temp))
    elif over_ice:
        return np.exp(ice(temp))
    else:
        return np.where(temp < 273.16, np.exp(ice(temp)), np.exp(liquid(temp)))


def Wright(t, **kwargs):
    """Wright (1997), US Meteorological Handbook:
    Should be used for liquid only

    ew = 611.21 *math.exp(17.502 * (t -273.15) / (240.97 + t - 273.15) )

    Args:
        t: temperatur in K

    Returns:
        es: saturation water vapor pressure in Pa
    """
    return 611.21 * np.exp(17.502 * (t - 273.15) / (240.97 + t - 273.15))
