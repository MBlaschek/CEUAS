#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__version__ = '0.1'
__author__ = 'MB'
__status__ = 'dev'
__date__ = 'Mit Jan 15 15:30:59 CET 2020'
__institute__ = 'Univie, IMGW'
__github__ = 'git@github.com:MBlaschek/CEUAS.git'
__doc__ = """
Radiosonde IGRAv2 Reading Software v%s
Maintained by %s at %s
Github: %s [%s]
License: C3S
Updated: %s
""" % (__version__, __author__, __institute__, __github__, __status__, __date__)


def read_IGRA(filename, odir, **kwargs):
    """Read IGRA version 2 Data from NOAA

    Args:
        filename (str): Filename
        odir (str): Output directory
    Returns:
        bool : Status
        str : Message (Error)

    Info:
        Format Description of IGRA 2 Sounding Data Files

    ---------------------
    Notes:
    ---------------------

    2. Both types of files are updated once a day in the early morning Eastern
       Time. The latest observations usually become available within two
       calendar days of when they were taken.

    2. Data files are available for two different time spans:

       In subdirectory data-por, data files contain the full period of record.
       In subdirectory data-y2d, files only contain soundings from the current
         (or current and previous) year. For example, as of August 2016,
         the files in the data-y2d subdirectory begin with January 1, 2016.

    3. Each file in the data-por and data-y2d subdirectories contains the
       sounding data for one station.
       The name of the file corresponds to a station's IGRA 2 identifier (e.g.,
       "USM00072201-data.txt.zip"  contains the data for the station with the
       identifier USM00072201).

    3. Each sounding consists of one header record and n data
       records, where n (given in the header record) is the number of levels
       in the sounding.

    ---------------------
    Header Record Format:
    ---------------------

    ---------------------------------
    Variable   Columns  Type
    ---------------------------------
    HEADREC       1-  1  Character
    ID            2- 12  Character
    YEAR         14- 17  Integer
    MONTH        19- 20  Integer
    DAY          22- 23  Integer
    HOUR         25- 26  Integer
    RELTIME      28- 31  Integer
    NUMLEV       33- 36  Integer
    P_SRC        38- 45  Character
    NP_SRC       47- 54  Character
    LAT          56- 62  Integer
    LON          64- 71  Integer
    ---------------------------------

    These variables have the following definitions:

    HEADREC		is the header record indicator (always set to "#").

    ID		is the station identification code. See "igra2-stations.txt"
            for a complete list of stations and their names and locations.

    YEAR 		is the year of the sounding.

    MONTH 		is the month of the sounding.

    DAY 		is the day of the sounding.

    HOUR 		is the nominal or observation hour of the sounding (in UTC on
            the date indicated in the YEAR/MONTH/DAY fields). Possible
            valid hours are 00 through 23, and 99 = missing. Hours are
            given as provided by the data provider, and the relationship
            between this hour and the release time varies by data
            provider, over time, and among stations.

    RELTIME 	is the release time of the sounding in UTC. The format is
            HHMM, where HH is the hour and MM is the minute. Possible
            are 0000 through 2359, 0099 through 2399 when only the release
            hour is available, and 9999 when both hour and minute are
            missing.

    NUMLEV 		is the number of levels in the sounding (i.e., the number of
            data records that follow).

    P_SRC 		is the data source code for pressure levels in the sounding.
            It has 25 possible values:

            bas-data = British Antarctic Survey READER Upper-Air Data
            cdmp-amr = African Monthly Radiosonde Forms
                       digitized by the U.S. Climate Data Modernization
                       Program
            cdmp-awc = "African Wind Component Data" digitized from
                       Monthly Forms by the U.S. Climate Data
                       Modernization Program
            cdmp-mgr = "WMO-Coded Messages" for Malawi, digitized from
                       "Computer-Generated Forms" by the U.S. Climate
                       Data Modernization Program
            cdmp-zdm = Zambian "Daily UA MB Ascent Sheets" digitized by
                       the U.S. Climate Data Modernization Program
            chuan101 = Comprehensive Historical Upper Air Network (v1.01)
            erac-hud = ERA-CLIM Historical Upper Air Data
            iorgc-id = IORGC/JAMSTEC-Digitized data for Indonesia
            mfwa-ptu = West African Temperature-Humidity Soundings
                       digitized by Meteo-France
            ncar-ccd = C-Cards Radiosonde Data Set from NCAR
            ncar-mit = MIT Global Upper Air Data from NCAR
            ncdc6210 = NCDC Marine Upper Air Data (NCDC DSI-6210)
            ncdc6301 = NCDC U.S. Rawindsonde Data (NCDC DSI-6301)
            ncdc6309 = NCDC "NCAR-NMC Upper Air" (NCDC DSI-6309)
            ncdc6310 = NCDC "Global U/A Cards" (NCDC DSI-6310)
            ncdc6314 = Global Telecommunications System messages received
                       and processed at Roshydromet and archived at NCDC
                       (NCDC DSI-6314)
            ncdc6315 = NCDC "People's Republic of China Data" (NCDC DSI-6315)
            ncdc6316 = NCDC "Argentina National Upper Air Data" (NCDC
                       DSI-6316)
            ncdc6319 = NCDC "Korea National Upper Air Data" (NCDC DSI-6319)
            ncdc6322 = Global Telecommunications System messages received
                       at the Australian Bureau of Meteorology and
                       archived at NCDC (NCDC DSI-6322)
            ncdc6323 = NCDC "Australian U/A Thermo/Winds Merged" (NCDC
                       DSI-6323)
            ncdc6324 = NCDC "Brazil National Upper Air Data" (NCDC DSI-6324)
            ncdc6326 = NCDC "Global Upper Air Cards" (NCDC DSI-6326)
            ncdc6355 = Russian Ice Island upper air data  processed by
                       NCAR and archived at NCDC
            ncdc-gts = Global Telecommunications System (GTS) messages
                       received at NCDC from the National Centers for
                       Environmental Prediction
            ncdc-nws =  U.S. National Weather Service upper air data
                        received at NCDC in real-time
            ngdc-har = Historical Arctic radiosonde archive from the
                       National Geophysical Data Center
            usaf-ds3 = U.S. Air Force 14th Weather Squadron Upper Air
                       Data Set ( received in DS3 format)

    NP_SRC 		is the data source code for non-pressure levels in the
            sounding. These include levels whose vertical coordinate
            is only identified by height as well as surface levels without
            either pressure or height.
            NP_SRC has 15 possible values:

            cdmp-adp = "African Daily Pilot Balloon Ascent Sheets" digitized
                       by the U.S. Climate Data Modernization Program
            cdmp-awc = "African Wind Component Data" digitized from
                       "Monthly Forms" by the U.S. Climate Data
                       Modernization Program
            cdmp-us2 = "U.S. Winds Aloft digitized from "Daily Computation
                       Sheets" by the U.S. Climate Data Modernization
                       Program
            cdmp-us3 = "U.S. Winds Aloft" digitized from "Military Daily
                       Computation Sheets" by the U.S. Climate Data
                       Modernization Program
            cdmp-usm = U.S. pilot balloon observations digitized from
                       "Monthly Forms" by the U.S. Climate Data
                       Modernization Program
            chuan101 = Comprehensive Historical Upper Air Network (v1.01)
            erac-hud = ERA-CLIM Historical Upper Air Data
            mfwa-wnd = West African Winds Aloft digitized by Meteo-France
            ncdc6301 = NCDC U.S. Rawindsonde Data (NCDC DSI-6301)
            ncdc6309 = NCDC "NCAR-NMC Upper Air" (NCDC DSI-6309)
            ncdc6314 = Global Telecommunications System messages received
                       and processed at Roshydromet and archived at NCDC
                       (NCDC DSI-6314)
            ncdc-gts = Global Telecommunications System (GTS) messages
                       received at NCDC from the National Centers for
                       Environmental Prediction
            ncdc-nws =  U.S. National Weather Service upper air data
                        received at NCDC in real-time
            ngdc-har = Historical Arctic radiosonde archive from the
                       National Geophysical Data Center
            usaf-ds3 = U.S. Air Force 14th Weather Squadron Upper Air
                       Data Set (received in DS3 format)

    LAT 		is the Latitude at which the sounding was taken. For mobile
            stations, it is the latitude at the time of observation.
            For fixed stations, it is the same as the latitude shown
            in the IGRA station list regardless of the date of the
            sounding since no attempt was made to reconstruct the
            sounding-by-sounding location history of these stations.

    LON 		is the longitude at which the sounding was taken. For mobile
            stations, it is the longitude at the time of observation.
            For fixed stations, it is the same as the longitude shown
            in the IGRA station list regardless of the date of the
            sounding since no attempt was made to reconstruct the
            sounding-by-sounding location history of these stations.

    ---------------------
    Data Record Format:
    ---------------------

    -------------------------------
    Variable        Columns Type
    -------------------------------
    LVLTYP1         1-  1   Integer
    LVLTYP2         2-  2   Integer
    ETIME           4-  8   Integer
    PRESS          10- 15   Integer
    PFLAG          16- 16   Character
    GPH            17- 21   Integer
    ZFLAG          22- 22   Character
    TEMP           23- 27   Integer
    TFLAG          28- 28   Character
    RH             29- 33   Integer
    DPDP           35- 39   Integer
    WDIR           41- 45   Integer
    WSPD           47- 51   Integer
    -------------------------------

    These variables have the following definitions:

    LVLÞTYP1 	is the major level type indicator. It has the following
            three possible values:

            1 = Standard pressure level (for levels at 1000, 925, 850,
                700, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30,
                20, 10, 7, 5, 3, 2, and 1 hPa)
            2 = Other pressure level
            3 = Non-pressure level

    LVLÞTYP2 	is the minor level type indicator. It has the following
            three possible values:

            1 = Surface
            2 = Tropopause
            0 = Other

    ETIME		is the elapsed time since launch. The format is MMMSS, where
            MMM represents minutes and SS represents seconds, though
            values are not left-padded with zeros. The following special
            values are used:

            -8888 = Value removed by IGRA quality assurance, but valid
                    data remain at the same level.
            -9999 = Value missing prior to quality assurance.

    PRESS 		is the reported pressure (Pa or mb * 100, e.g.,
            100000 = 1000 hPa or 1000 mb). -9999 = missing.

    PFLAG 		is the pressure processing flag indicating what level of
            climatology-based quality assurance checks were applied. It
            has three possible values:

            blank = Not checked by any climatology checks. If data value
                    not equal to -9999, it passed all other applicable
                    checks.
            A     = Value falls within "tier-1" climatological limits
                    based on all days of the year and all times of day
                    at the station, but not checked by
                    "tier-2" climatology checks due to
                    insufficient data.
            B     = Value passes checks based on both the tier-1
                    climatology and a "tier-2" climatology specific to
                    the time of year and time of day of the data value.

    GPH 		is the reported geopotential height (meters above sea level).
            This value is often not available at variable-pressure levels.
            The following special values are used:

            -8888 = Value removed by IGRA quality assurance, but valid
                    data remain at the same level.
            -9999 = Value missing prior to quality assurance.

    ZFLAG 		is the  geopotential height processing flag indicating what
            level of climatology-based quality assurance checks were
            applied. It has three possible values:

            blank = Not checked by any climatology checks or flag not
                    applicable. If data value not equal to -8888 or -9999,
                    it passed all other applicable checks.
            A     = Value falls within "tier-1" climatological limits
                    based on all days of the year and all times of day
                    at the station, but not checked by
                    "tier-2" climatology checks due to insufficient data.
            B     = Value passes checks based on both the tier-1
                    climatology and a "tier-2" climatology specific to
                    the time of year and time of day of the data value.

    TEMP 		is the reported temperature (degrees C to tenths, e.g.,
            11 = 1.1°C). The following special values are used:

            -8888 = Value removed by IGRA quality assurance, but valid
                    data remain at the same level.
            -9999 = Value missing prior to quality assurance.

    TFLAG 		is the temperature processing flag indicating what
            level of climatology-based quality assurance checks were
            applied. It has three possible values:

            blank = Not checked by any climatology checks or flag not
                    applicable. If data value not equal to -8888 or -9999,
                    it passed all other applicable checks.
            A     = Value falls within "tier-1" climatological limits
                    based on all days of the year and all times of day
                    at the station, but not checked by "tier-2"
                    climatology checks due to insufficient data.
            B     = Value passes checks based on both the tier-1
                    climatology and a "tier-2" climatology specific to
                    the time of year and time of day of the data value.

    RH 		is the reported relative humidity (Percent to tenths, e.g.,
            11 = 1.1%). The following special values are used:

            -8888 = Value removed by IGRA quality assurance, but valid
                    data remain at the same level.
            -9999 = Value missing prior to quality assurance.
    DPDP 		is the reported dewpoint depression (degrees C to tenths, e.g.,
            11 = 1.1°C). The following special values are used:

            -8888 = Value removed by IGRA quality assurance, but valid
                    data remain at the same level.
            -9999 = Value missing prior to quality assurance.

    WDIR 		is the reported wind direction (degrees from north,
            90 = east). The following special values are used:

            -8888 = Value removed by IGRA quality assurance, but valid
                    data remain at the same level.
            -9999 = Value missing prior to quality assurance.

    WSPD 		is the reported wind speed (meters per second to tenths, e.g.,
            11 = 1.1 m/s). The following special values are used:

            -8888 = Value removed by IGRA quality assurance, but valid
                    data remain at the same level.
            -9999 = Value missing prior to quality assurance.
    """
    import datetime
    import gzip
    import zipfile
    import os
    import io
    import numpy as np
    import pandas as pd

    try:
        if not os.path.isfile(filename):
            raise IOError("File not found", filename)

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
        for i, line in enumerate(data):
            if line[0] == '#':
                # Header
                ident = line[1:12]
                year = line[13:17]
                month = line[18:20]
                day = line[21:23]
                hour = line[24:26]
                reltime = line[27:31]
                numlev = int(line[32:36])
                p_src = line[37:45]
                np_src = line[46:54]
                lat = int(line[55:62]) / 10000.
                lon = int(line[63:71]) / 10000.
                #
                #
                #
                if int(hour) == 99:
                    time = reltime + '00'
                else:
                    time = hour + '0000'
                #
                # minutes can be missing
                #
                if '99' in time:
                    time = time.replace('99', '00')

                idate = datetime.datetime.strptime(year + month + day + time, '%Y%m%d%H%M%S')
                headers.append((idate, ident, numlev, p_src.strip(), np_src.strip(), lat, lon))

            else:
                # Data
                lvltyp1 = int(line[0])  # 1-  1   integer
                lvltyp2 = int(line[1])  # 2-  2   integer
                etime = int(line[3:8])  # 4-  8   integer
                press = int(line[9:15])  # 10- 15   integer
                pflag = line[15]  # 16- 16   character
                gph = int(line[16:21])  # 17- 21   integer
                zflag = line[21]  # 22- 22   character
                temp = int(line[22:27]) / 10.  # 23- 27   integer
                tflag = line[27]  # 28- 28   character
                rh = int(line[28:33]) / 10.  # 30- 34   integer
                dpdp = int(line[34:39]) / 10.  # 36- 40   integer
                wdir = int(line[40:45])  # 41- 45   integer
                wspd = int(line[46:51]) / 10.  # 47- 51   integer

                raw.append((lvltyp1, lvltyp2, etime, press, pflag, gph, zflag, temp, tflag, rh, dpdp, wdir, wspd))
                dates.append(idate)

        print("IGRAv2 Lines read:", i, "Header count:", len(headers), **kwargs)
        c = ['ltyp1', 'ltyp2', 'etime', 'pres', 'pflag', 'gph', 'zflag', 'temp', 'tflag', 'rhumi', 'dpd', 'windd',
             'winds']
        out = pd.DataFrame(data=raw, index=dates, columns=c)
        out = out.replace([-999.9, -9999, -8888, -888.8], np.nan)  # known missing values by IGRAv2
        out.index.name = 'date'
        headers = pd.DataFrame(data=headers,
                               columns=['date', 'ident', 'numlev', 'p_src', 'np_src', 'lat', 'lon']).set_index(
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
read_IGRAv2_to_netcdf.py -d [dir] -p [pat] -o [out]

Options:
    -d []       Input directory, e.g.: /tmp/data/igra
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
        ret, code = read_IGRA(ifile, odir)
        print("{}/{} {} [{} : {}]".format(i, n, ifile, ret, code))
    #
    # FIN
    #
    t0 = time.time()
    print(t1 - t0, "s | Job Executed: ", " ".join(sys.argv))
