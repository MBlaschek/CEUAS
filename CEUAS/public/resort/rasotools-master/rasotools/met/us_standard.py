# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import datetime
import netCDF4 as nc

__all__ = ['atmospheres']


def atmospheres(filename=None, profile=None, levels=None, interpolate=False):
    """Standard Atmospheres

    Profile Types:
          TROPICAL                        1
          MIDLATITUDE SUMMER              2
          MIDLATITUDE WINTER              3
          SUBARCTIC SUMMER                4
          SUBARCTIC WINTER                5
          U. S. STANDARD,  1976           6

    from RTTOV v11
    
    Options:    filename
                profile          Profile Number
                interpolation    Interpolated to ERA-Interim Levels?

    Output:   profile , sonde_info
              as Pandas DataFrames
    """
    from .. import config

    if filename is None:
        filename = config.get_data('us_standard_76.H5')

    print("""
TROPICAL                        0
MIDLATITUDE SUMMER              1
MIDLATITUDE WINTER              2
SUBARCTIC SUMMER                3
SUBARCTIC WINTER                4
U. S. STANDARD,  1976           5
""")

    f = nc.Dataset(filename, 'r')

    info = pd.DataFrame()

    if profile is not None and isinstance(profile, str):
        profile = profile.upper()
        profs_infile = [f['PROFILES'][i]['ID'][:] for i in f['PROFILES'].groups.keys()]
        if profile in profs_infile:
            profile = profs_infile.index(profile) + 1
        else:
            raise ValueError("Could not find Profile: %s" % profile)

    if levels is None:
        levels = config.era_plevels

    if 0 < profile < 7:
        iobj = f['PROFILES']['%04d' % profile]
        print("Selecting: ", iobj['ID'][:])
        info = pd.concat([info, pd.DataFrame({'lon': float(iobj['LONGITUDE'][0]),
                                              'lat': float(iobj['LATITUDE'][0]),
                                              'time': datetime.datetime(*(list(iobj['DATE'][:]) + list(iobj['TIME'][:]))),
                                              'id': iobj['ID'][:]}, index=[int(profile)])])
        if interpolate:
            tnew = _interp_mod(iobj['T'][:], np.log(levels), np.log(iobj['P'][:] * 100.), False)
            qnew = _interp_mod(iobj['Q'][:], np.log(levels), np.log(iobj['P'][:] * 100.), False)
            out = pd.DataFrame({'p': levels, 't': tnew, 'q': qnew})

        else:
            out = pd.DataFrame({'p': iobj['P'][:] * 100., 't': iobj['T'][:], 'q': iobj['Q'][:]})
    else:
        out = {}
        for ikey, iobj in f['PROFILES'].groups.items():
            info = pd.concat([info, pd.DataFrame({'lon': float(iobj['LONGITUDE'][0]), 'lat': float(iobj['LATITUDE'][0]),
                                                  'time': datetime.datetime(*(list(iobj['DATE'][:]), list(iobj['TIME'][:]))),
                                                  'id': iobj['ID'][:]}, index=[int(ikey)])])
            if interpolate:
                tnew = _interp_mod(iobj['T'][:], np.log(levels), np.log(iobj['P'][:] * 100.), False)
                qnew = _interp_mod(iobj['Q'][:], np.log(levels), np.log(iobj['P'][:] * 100.), False)
                out[iobj['ID'].value] = pd.DataFrame({'p': levels, 't': tnew, 'q': qnew})
            else:
                out[iobj['ID'].value] = pd.DataFrame({'p': iobj['P'][:] * 100., 't': iobj['T'][:], 'q': iobj['Q'][:]})

    print("Units: p: %s t: %s q: %s " % ('Pa', iobj['T'].UNITS, iobj['Q'].UNITS))
    f.close()
    return out, info


# INTERPOLATION FUNCTION
def _interp_mod(values, pout, pin, dropna=True, min_values=4):
    """
    Modified numpy interp function

    Args:
      values     values to be interpolated
      pout       interpolation levels
      pin        input levels

    Keyword Args:
      dropna      [True]   Remove missing values before interpolation ?
      min_values  [4]      Minimum required levels for interpolation

    Results:
      values interpolated to levels
    """
    # remove NAN
    if dropna:
        ix = np.isnan(values)
        values = values[~ix]
        pin = pin[~ix]

    if np.sum(np.isfinite(values)) >= min_values:
        return np.interp(pout, pin, values, left=np.nan, right=np.nan)  # no extrapolation
    else:
        return np.full(len(pout), np.nan)

# download regression limits from RTTOV
# https://nwpsaf.eu/downloads/rtcoef_rttov12/PROFILES_ECMWF_83_2016_54L_regression_limits.txt?8447b5&8447b5
