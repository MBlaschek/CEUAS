#  Writing data in TM309099 - high-resolution sondes

#  Starting point generated using ecCodes version: 2.23.0

# run with  strace -o trace.txt .... to see what path is being accessed

# export ECCODES_DEFINITION_PATH=/home/erc/Work/mydefsWMOBUFR:`codes_info -d`

from __future__ import print_function
import traceback
import sys, io
from eccodes import *
import numpy as np
from encode_rdb import encode_rdb, encode_profile

def bufr_encode(profile1, bufr_file1):
    ledition4=False
    if ledition4:
      ibufr = codes_bufr_new_from_samples('BUFR4')
    else:
      ibufr = codes_bufr_new_from_samples('BUFR3_local')
    ivalues = (1 ,)
    codes_set_array(ibufr, 'inputDelayedDescriptorReplicationFactor', ivalues)
    nlevs = len(profile1['data'][list(profile1['data'].keys())[0]])
    ivalues = (nlevs , 10*nlevs, 4*nlevs)
    codes_set_array(ibufr, 'inputExtendedDelayedDescriptorReplicationFactor', ivalues)
    if ledition4:
      codes_set(ibufr, 'edition', 4)
    else:
      codes_set(ibufr, 'edition', 3)
    codes_set(ibufr, 'masterTableNumber', 0)
    codes_set(ibufr, 'bufrHeaderCentre', 98)
    codes_set(ibufr, 'bufrHeaderSubCentre', 0)
    codes_set(ibufr, 'updateSequenceNumber', 0)
    codes_set(ibufr, 'dataCategory', 2) # see page https://www.nco.ncep.noaa.gov/sib/jeff/TableA_0_STDv22_LOC7.html
    if profile1['header']['platform_type'] > 8 or profile1['header']['platform_type'] < 0:
      subtype = {0:109,2:111,8:109}[0]
    else:
      subtype = {0:109,2:111,8:109}[profile1['header']['platform_type']] ## EDIT -> why can it now be nan, and didn't create an error before?
    #  subtype=111 # TEMP SHIP
    #  subtype=109 # TEMP LAND
    codes_set(ibufr, 'dataSubCategory', subtype)
    codes_set(ibufr, 'masterTablesVersionNumber', 31)
    codes_set(ibufr, 'localTablesVersionNumber', 6)

    rdbdata={}
    rdbdata['rdbType'] = 5
    rdbdata['rdbSubType'] = subtype
    rdbdata['ident'] = profile1['header']['shipOrMobileLandStationIdentifier']
    rdbdata['restricted'] = 0 # {1:1, 2:1,3:0, 4:0, 5:1, 6:1, 7:1, 8:1}[int(profile1['header']['datasetSource'][-1])] # TODO
    # CUON Data sources: 1:ERA5_1 2:ERA5_2 3:IGRA2 4:NCAR 5:BUFR 6:ERA5_1759 7:ERA5_1761 8:ERA5 3188
    for k in ['year','month','day','hour','minute','latitude','longitude']:
      rdbdata[k] = profile1['header'][k]
    if ledition4:
      codes_set(ibufr, 'typicalYear', rdbdata['year'])
      codes_set(ibufr, 'typicalMonth', rdbdata['month'])
      codes_set(ibufr, 'typicalDay', rdbdata['day'])
      codes_set(ibufr, 'typicalHour', rdbdata['hour'])
      codes_set(ibufr, 'typicalMinute', rdbdata['minute'])
      codes_set(ibufr, 'typicalSecond', 0)
    else:
      encode_rdb( codes_set, ibufr, rdbdata )

    codes_set(ibufr, 'numberOfSubsets', 1)
    codes_set(ibufr, 'observedData', 1)
    codes_set(ibufr, 'compressedData', 0)

    bitmapMask = []
    for ilev in range(nlevs):
      bitmapMask += [1, 1, 1, 1, 1, 1, 0, 0, 0, 0]
    codes_set_array (ibufr, 'inputDataPresentIndicator',bitmapMask)

    # Create the structure of the data section
    codes_set_array(ibufr, 'unexpandedDescriptors', [206064, 1099, 301150, 301111, 301128, 301113, 301114, 302049, 22043, 101000, 31002, 303056, 225000, 236000, 101000, 31002, 31031, 8024, 101000, 31002, 225255, 101000, 31001, 303051])                    
    # codes_set_array(ibufr, 'unexpandedDescriptors', [206064, 1231, 301150, 301111, 301128]) # , 301113, 301114, 302049, 22043, 101000, 31002, 303056, 225000, 236000, 101000, 31002, 31031, 8024, 101000, 31002, 225255, 101000, 31001, 303051])                    
    # [206064, 1231, 301150, 2231, 301111, 301128, 301113, 301114, 302049, 22043, 101000, 31002, 303056, 225000, 236000, 101000, 31002, 31031, 8024, 101000, 31002, 225255, 101000, 31001, 303051]
    # , 2231
          #           [ 206064, 1231] +
          #  [301150,2231,301111,301128,301113,301114,302049,22043] +
          #  [101000,31002,303056] + 
          #  #[4086,8042,207001,7004,10009,207000,5015,6015,12101,12103,11001,11002] + 
          #  #THIS WORKS: [222000,236000,101000,31002,31031,101000,31002,33007] +
          #  [225000,236000,101000,31002,31031,8024,101000,31002,225255] + # see /home/erc/ifs-source/erc_CY49R1_ISPD.IFS-3390.ED-85/odb/old_tools/bufr_add_bias.F
          #  [101000,31001,303051] )

    for k in ['wigosIdentifierSeries','wigosIssuerOfIdentifier','wigosIssueNumber','wigosLocalIdentifierCharacter','blockNumber','stationNumber','shipOrMobileLandStationIdentifier','height','year','month','day','hour','minute','latitude','longitude']:
      if k in profile1['header'].keys():
        codes_set(ibufr, k, profile1['header'][k])

    if 'radiosondeType' in profile1['header'].keys():
      codes_set(ibufr, 'radiosondeType', profile1['header']['radiosondeType'])
    else:
      codes_set(ibufr, 'radiosondeType', CODES_MISSING_LONG)

    codes_set(ibufr, 'solarAndInfraredRadiationCorrection', CODES_MISSING_LONG)
    codes_set(ibufr, 'trackingTechniqueOrStatusOfSystem', CODES_MISSING_LONG)
    codes_set(ibufr, 'measuringEquipmentType', CODES_MISSING_LONG)

    codes_set(ibufr, 'uniqueProductDefinition', profile1['header']['datasetSource'])
    # codes_set(ibufr, 'datasetSource', profile1['header']['datasetSource'])

#WMOBUFR    if 'sondeTypeDetail' in profile1['header'].keys():
#WMOBUFR      codes_set(ibufr, 'sondeTypeDetail', profile1['header']['sondeTypeDetail'])
#WMOBUFR    else:
#WMOBUFR      codes_set(ibufr, 'sondeTypeDetail', CODES_MISSING_LONG)

    codes_set(ibufr, 'radiosondeSerialNumber','')
    codes_set(ibufr, 'radiosondeAscensionNumber', CODES_MISSING_LONG)
    codes_set(ibufr, 'radiosondeReleaseNumber', CODES_MISSING_LONG)
    codes_set(ibufr, 'observerIdentification','')
    codes_set(ibufr, 'radiosondeCompleteness', CODES_MISSING_LONG)
    codes_set(ibufr, 'radiosondeConfiguration', CODES_MISSING_LONG)
    codes_set(ibufr, 'correctionAlgorithmsForHumidityMeasurements', CODES_MISSING_LONG)
    codes_set(ibufr, 'radiosondeGroundReceivingSystem', CODES_MISSING_LONG)
    codes_set(ibufr, 'radiosondeOperatingFrequency', CODES_MISSING_DOUBLE)
    codes_set(ibufr, 'balloonManufacturer', CODES_MISSING_LONG)
    codes_set(ibufr, 'balloonType', CODES_MISSING_LONG)
    codes_set(ibufr, 'weightOfBalloon', CODES_MISSING_DOUBLE)
    codes_set(ibufr, 'balloonShelterType', CODES_MISSING_LONG)
    codes_set(ibufr, 'typeOfGasUsedInBalloon', CODES_MISSING_LONG)
    codes_set(ibufr, 'amountOfGasUsedInBalloon', CODES_MISSING_DOUBLE)
    codes_set(ibufr, 'balloonFlightTrainLength', CODES_MISSING_DOUBLE)
    codes_set(ibufr, 'pressureSensorType', CODES_MISSING_LONG)
    codes_set(ibufr, 'temperatureSensorType', CODES_MISSING_LONG)
    codes_set(ibufr, 'humiditySensorType', CODES_MISSING_LONG)
    codes_set(ibufr, 'radome', CODES_MISSING_LONG)
    codes_set(ibufr, 'geopotentialHeightCalculation', CODES_MISSING_LONG)
    codes_set(ibufr, 'softwareVersionNumber','')
    codes_set(ibufr, 'reasonForTermination', CODES_MISSING_LONG)
    codes_set(ibufr, 'timeSignificance', CODES_MISSING_LONG)
    codes_set(ibufr, 'heightOfBarometerAboveMeanSeaLevel', CODES_MISSING_DOUBLE)
    codes_set(ibufr, 'heightOfStationGroundAboveMeanSeaLevel', CODES_MISSING_LONG)
    codes_set(ibufr, 'stationElevationQualityMarkForMobileStations', CODES_MISSING_LONG)
    codes_set(ibufr, '#1#verticalSignificanceSurfaceObservations', CODES_MISSING_LONG)
    codes_set(ibufr, 'cloudAmount', CODES_MISSING_LONG)
    codes_set(ibufr, 'heightOfBaseOfCloud', CODES_MISSING_DOUBLE)
    codes_set(ibufr, '#1#cloudType', CODES_MISSING_LONG)
    codes_set(ibufr, '#2#cloudType', CODES_MISSING_LONG)
    codes_set(ibufr, '#3#cloudType', CODES_MISSING_LONG)
    codes_set(ibufr, '#2#verticalSignificanceSurfaceObservations', CODES_MISSING_LONG)
    codes_set(ibufr, 'oceanographicWaterTemperature', CODES_MISSING_DOUBLE)
    codes_set(ibufr, 'differenceStatistics', 11) 

    level_kays = [('timePeriod',CODES_MISSING_LONG,int), \
                  ('extendedVerticalSoundingSignificance',CODES_MISSING_LONG,int), \
                  ('pressure',CODES_MISSING_LONG,int), \
                  ('nonCoordinateGeopotentialHeight',CODES_MISSING_DOUBLE,float), \
                  ('latitudeDisplacement',CODES_MISSING_DOUBLE,float),('longitudeDisplacement',CODES_MISSING_DOUBLE,float), \
                  ('airTemperature',CODES_MISSING_DOUBLE,float), \
                  ('airTemperature->differenceStatisticalValue',CODES_MISSING_DOUBLE,float), \
                  ('dewpointTemperature',CODES_MISSING_DOUBLE,float), \
                  ('dewpointTemperature->differenceStatisticalValue',CODES_MISSING_DOUBLE,float), \
                  ('windDirection',CODES_MISSING_LONG,int), \
                  ('windDirection->differenceStatisticalValue',CODES_MISSING_LONG,int), \
                  ('windSpeed',CODES_MISSING_DOUBLE,float), \
                  ('windSpeed->differenceStatisticalValue',CODES_MISSING_DOUBLE,float), \
                  ]

    encode_profile(codes_set, ibufr, level_kays, profile1['data'], nlevs1=nlevs)
    #data_kays = [] ; miss_kays = []
    #for k in level_kays:
    #  if k[0] in profile1['data'].keys():
    #    data_kays.append(k)
    #  else:
    #    miss_kays.append(k)
#
#    #print(data_kays)
#    #print(miss_kays)
#    for ilev in range(nlevs):
#      ient = '#{0}#'.format(ilev+1)
#      for k in data_kays:
#        codes_set(ibufr, ient+ k[0], check_value_before_bufr_encoding( profile1['data'][k[0]][ilev], k[1], k[2] ) )
#      for k in miss_kays:
#        codes_set(ibufr, ient+ k[0], k[1] )
    #codes_set(ibufr, '#2#timePeriod', CODES_MISSING_LONG)
    #codes_set(ibufr, '#2#extendedVerticalSoundingSignificance', CODES_MISSING_LONG)
    #codes_set(ibufr, '#2#pressure', CODES_MISSING_DOUBLE)
    #codes_set(ibufr, '#2#latitudeDisplacement', CODES_MISSING_DOUBLE)
    #codes_set(ibufr, '#2#longitudeDisplacement', CODES_MISSING_DOUBLE)
    #codes_set(ibufr, 'absoluteWindShearIn1KmLayerBelow', CODES_MISSING_DOUBLE)
    #codes_set(ibufr, 'absoluteWindShearIn1KmLayerAbove', CODES_MISSING_DOUBLE)

    # Encode the keys back in the data section
    codes_set(ibufr, 'pack', 1)

    if type(bufr_file1) is str:
      outfile = open(bufr_file1, 'wb')
    elif type(bufr_file1) is io.BufferedWriter:
      outfile = bufr_file1
    else:
      print("ERROR... unexpected type for the file: neither string nor io.BufferedWriter")
      stop

    codes_write(ibufr, outfile)
    if type(bufr_file1) is str:
      outfile.close()
      print ("Created output BUFR file '{0}'".format(bufr_file1))
    codes_release(ibufr)

"""
def main():
    T_P = np.array([[   281.85, 101800.  ],
       [   280.55, 100000.  ],
       [   278.75,  96700.  ],
       [   280.75,  96000.  ],
       [   280.95,  95000.  ],
       [   279.45,  92500.  ],
       [   277.15,  89200.  ],
       [   277.75,  88500.  ],
       [   276.75,  85000.  ],
       [   276.15,  81900.  ],
       [   273.75,  80500.  ],
       [   276.15,  77800.  ],
       [   274.95,  75300.  ],
       [   274.55,  72100.  ],
       [   272.65,  70000.  ],
       [   265.05,  62100.  ],
       [   262.85,  58000.  ],
       [   257.35,  53000.  ],
       [   258.35,  52500.  ],
       [   256.25,  50000.  ],
       [   243.85,  40000.  ],
       [   239.85,  37400.  ],
       [   234.25,  32800.  ],
       [   228.85,  30000.  ],
       [   222.05,  26500.  ],
       [   219.65,  25000.  ],
       [   218.25,  24100.  ],
       [   220.05,  20000.  ],
       [   218.45,  19000.  ],
       [   224.65,  16400.  ],
       [   224.05,  15000.  ],
       [   223.45,  10000.  ],
       [   224.85,   7000.  ],
       [   226.85,   5000.  ],
       [   226.45,   3340.  ],
       [   228.65,   3000.  ],
       [   229.25,   2800.  ]])
    oneprofile={'header':{'longitude':337.4-360., 'latitude':63.97, 'height':425.0, \
                  'year':1997, 'month':7, 'day':1, 'hour':0, 'minute':0, 'second':0, \
                  'wigosIdentifierSeries':0, 'wigosIssuerOfIdentifier':20000, 'wigosIssueNumber':0, 'wigosLocalIdentifierCharacter':'4018', 'blockNumber':4, 'stationNumber':18, 'shipOrMobileLandStationIdentifier':''}, \
     'data':{'pressure':[x[1] for x in T_P],'airTemperature':[x[0] for x in T_P],'airTemperatureBiasCorrection':[0.01 for x in T_P]}}
    try:
        bufr_encode(oneprofile)
    except CodesInternalError as err:
        traceback.print_exc(file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())

"""
