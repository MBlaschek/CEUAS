#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__version__ = '0.1'
__author__ = 'MB'
__status__ = 'dev'
__date__ = 'Mit Jan 15 15:30:59 CET 2020'
__institute__ = 'Univie, IMGW'
__github__ = 'git@github.com:MBlaschek/CEUAS.git'
__doc__ = """
Radiosonde UADB Reading Software v%s
Maintained by %s at %s
Github: %s [%s]
License: C3S
Updated: %s
""" % (__version__, __author__, __institute__, __github__, __status__, __date__)


def read_UADB(filename, odir, **kwargs):
    """ NCAR Upper Air Database
    This data is output from the NCAR Upper Air Database Project (UADB). The Composited
    UADB products (UADB-TRH,UADB-Wind) and Combined (UADB-TRHC,UADB-WindC)
    products have been output in this format. The output contains 1 descriptive header record (Table
    C-1) for each sounding, followed by the data records (Table C-2) for that sounding, with 1 line
    for each level. Note that each field in both the header and data records is separated by a space.

    Args:
        filename (str): filename
        **kwargs:

    Returns:


    Documentation:
            http://rda.ucar.edu/datasets/ds370.1/docs/uadb-format-ascii.pdf

        There is one header record for each sounding, followed by all the data records,
        one level per record.  The nlevs variable indicates how many data levels follow
        the header

        The values for all the flags below can be found at
        http://dss.ucar.edu/datasets/ds370.0/docs/format.html

        HEADER VARIABLES
            VAR            Type    Format     Description
            ---            ----    ------     -----------
            ht            char    a1         should alway be an H to identify header record
            uid           int*4   i12        sounding id, integer used to uniquely identify each sounding within a source
            idc           char    a6         station identifier, usually a WMO, WBAN or call sign
            idtype        int     i2         Flag identifying station identifier type
            isrc          int     i3         Flag indentifying source of data record
            dssrc         real    5.2        Version number of the source record (this will be 2.xx)
            datetype      int     i2         Flag identifying accuracy of the date
            iyr           int     i4         4 digit year
            imn           int     i2         2 digit month
            idy           int     i2         2 digit day
            ihr           int     i4         4 digit hour, 1 2 digits are the hour and the last 2 digits are the minutes
            lflg          int     i2         Flag identifying location type
            xlat          real    f8.2       latitude, in degrees and hundredths -90 to 90
            xlon          real    f8.2       longitude, in degrees and hundredths, 0-360E
            xelv          real    f5.0       elevation in whole meters
            iobt          int     i3         Flag indentifying observation type
            nlevs         int     i4         # of levels in the soudings
            cmpvrs        char    a6         Prodcut Version # (aa.bb.cc)

        After each header, there will be nlevs data records that follow.
        The variables for each level are put into a 1 dimension array


        DATA VARIABLES
            VAR            Type    Format     Description
            ---            ----    ------     -----------
            ltyp            int     i4        Flag indicating level type
                                              1:surf, 2:mand, 3:sig, 4:wind by Z,
                                              5:wind by P, 6=maxWind, 7:trop
                                              8:stdNCDC, 9:sigWind, 11:wind by Z & P
              p            real    f8.1       Presssure, in mbs and tenths
                                                 missing=-99999.0
              z            real    f8.1       geopotential , in mbs and tenths
                                                 missing=-99999.0
              t            real    f6.1       temperature, in degrees C and tenths
                                                 missing=-999.0
             rh            real    f6.1       relative humidity (0-100) , in tenths
                                                 missing=-999.0
             wd            real    f6.0       wind direction, in whole degrees
                                                 missing=-999
             ws            real    f6.1       wind speed, m/s and tenths
                                                 missing=-999.0
    """
    import datetime
    import zipfile
    import gzip
    import os
    import io
    import numpy as np
    import pandas as pd

    try:
        if not os.path.isfile(filename):
            raise IOError("File not Found! %s" % filename)

        if '.zip' in filename:
            archive = zipfile.ZipFile(filename, 'r')
            inside = archive.namelist()
            tmp = archive.open(inside[0])
            tmp = io.TextIOWrapper(tmp, encoding='utf-8')
            tmp = tmp.read()
            archive.close()
            data = tmp.splitlines()  # Memory (faster)
        elif '.gz' in filename:

            with gzip.open(filename, 'rt', encoding='utf-8') as infile:
                tmp = infile.read()  # alternative readlines (slower)
                data = tmp.splitlines()  # Memory (faster)
        else:
            with open(filename, 'rt') as infile:
                tmp = infile.read()  # alternative readlines (slower)
                data = tmp.splitlines()  # Memory (faster)

        raw = []
        headers = []
        dates = []
        nmiss = 0
        iprev = 0
        search_h = False
        i = 0
        for i, line in enumerate(data):
            if line[0] == 'H':
                try:
                    # Header
                    usi = int(line[2:14])  # unique station identifier
                    ident = line[15:21]  # WMO
                    idflag = int(line[22:24])  # id flag
                    d_src = int(line[25:28])  # source dataset
                    version = float(line[29:34])  # version
                    dateflag = int(line[35:37])  # date flag
                    year = line[38:42]  # year
                    month = "%02d" % int(line[43:45])
                    day = "%2d" % int(line[46:48])
                    hour = line[49:53]
                    locflag = int(line[54:56])  # Location Flag
                    lat = float(line[57:67])
                    lon = float(line[68:78])
                    ele = float(line[79:85])
                    stype = int(line[86:88])
                    numlev = int(line[89:93])
                    pvers = line[94:102]
                    #
                    # minutes and hours can have 99
                    #
                    if '99' in hour:
                        hour = hour.replace('99', '00')

                    if '99' in day:
                        search_h = True
                        continue

                    minutes = int(hour) % 100
                    hour = "%02d" % (int(hour) // 100)

                    if minutes > 60 or minutes < 0:
                        minutes = 0

                    elif minutes == 60:
                        minutes = 59

                    else:
                        pass
                    minutes = "%02d" % minutes
                    idate = datetime.datetime.strptime(year + month + day + hour + minutes, '%Y%m%d%H%M')
                    headers.append((idate, usi, numlev, lat, lon, ele, stype))
                    pday = int(day)
                    search_h = False

                except Exception as e:
                    print("Error: ", i, line, repr(e), "Skipping Block:")
                    search_h = True
                    iprev = i

            elif search_h:
                nmiss += 1
                continue  # Skipping block

            else:
                # Data
                ltyp = int(line[0:4])
                press = float(line[5:13])  # hPa
                gph = float(line[14:22])
                temp = float(line[23:29])  # degree
                rh = float(line[30:36])  # %
                wdir = float(line[37:43])
                wspd = float(line[44:50])  # m/s
                raw.append((press, gph, temp, rh, wdir, wspd))
                dates.append(idate)

        print("UADB Lines read:", i, "skipped:", nmiss, "Header:", len(headers), **kwargs)

        out = pd.DataFrame(data=raw, index=dates, columns=['pres', 'gph', 'temp', 'rhumi', 'windd', 'winds'])
        out = out.replace([-999.9, -9999, -999, -999.0, -99999.0, -99999.9], np.nan)
        # fix units
        out['pres'] *= 100.  # need Pa
        out.index.name = 'date'
        headers = pd.DataFrame(data=headers, columns=['date', 'uid', 'numlev', 'lat', 'lon', 'alt', 'stype']).set_index(
            'date')
        #
        # Combine header and Data ?
        #

        #
        # to Netcdf
        # Add Attributes / data ?
        #


    except Exception as e:
        return False, repr(e)
    #
    # All Well
    #
    return True, "Finished"


def usage():
    print("""
read_UADBto_netcdf.py -d [dir] -p [pat] -o [out]

Options:
    -d []       Input directory, e.g.: /tmp/data/ncar
    -p []       Filename pattern, e.g.: "*.zip"
    -o []       Output directory, default: /tmp/sondes
""")
    print(__doc__)


if __name__ == '__main__':
    import sys
    import os
    import getopt
    import glob
    import time

    idir = None
    odir = '/tmp/sondes'
    pat = '*.zip'

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hd:p:o:", ["help", "output", "pattern"])

    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt in "-d":
            idir = arg
        elif opt in ("-p", "--pattern"):
            pattern = arg
        elif opt in ("-o", "--output"):
            odir = arg
        elif opt in ('-h', '--help'):
            usage()
            assert False
        else:
            usage()
            assert False, "unhandled option"

    if not os.path.isdir(idir):
        raise IOError("no such directory", idir)

    if not os.path.isdir(odir):
        print("Creating directory:", odir)
        os.makedirs(odir)

    files = glob.glob(idir + pat, recursive=True)

    if len(files) == 0:
        raise RuntimeError("No files found:", idir, pat)

    t1 = time.time()
    #
    # work
    #
    n = len(files)
    for i, ifile in enumerate(files):
        ret, code = read_UADB(ifile, odir)
        print("{}/{} {} [{} : {}]".format(i, n, ifile, ret, code))
    #
    # FIN
    #
    t0 = time.time()
    print(t1 - t0, "s | Job Executed: ", " ".join(sys.argv))
