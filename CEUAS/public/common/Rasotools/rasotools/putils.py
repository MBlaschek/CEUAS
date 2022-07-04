import numpy
import math
from datetime import date
import datetime
import time

# pythran export stationaverage(float[][][][], int)


def stationaverage(currentdata, thresh):

    s = currentdata.shape
    currentdataav = numpy.empty([s[1],s[2],s[3]])
    for it in range(s[3]):
        for ip in range(s[2]):
            for ipar in range(s[1]):
                mask = ~numpy.isnan(currentdata[:, ipar, ip, it])
                mask=mask.flatten()
                if sum(mask) >= thresh:
                    currentdataav[ipar, ip, it] = numpy.mean(
                        currentdata[mask, ipar, ip, it])
                else:
                    currentdataav[ipar, ip, it] = -1.e30
                    # currentdataav[ipar, ip, it] = numpy.nan

    return #currentdataav
