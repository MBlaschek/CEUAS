# Some routines useful for BUFR encoding
import numpy as np

def encode_rdb(codes_set1, ibfr1, rdbdata1):

    codes_set1(ibfr1, 'typicalYearOfCentury', rdbdata1['year'] % 100)
    codes_set1(ibfr1, 'typicalMonth', rdbdata1['month'])
    codes_set1(ibfr1, 'typicalDay', rdbdata1['day'])
    codes_set1(ibfr1, 'typicalHour', rdbdata1['hour'])
    codes_set1(ibfr1, 'typicalMinute', rdbdata1['minute'])

    codes_set1(ibfr1, 'rdbType', rdbdata1['rdbType'])
    codes_set1(ibfr1, 'oldSubtype', rdbdata1['rdbSubType'])

    codes_set1(ibfr1, 'localYear', rdbdata1['year'])
    codes_set1(ibfr1, 'localMonth', rdbdata1['month'])
    codes_set1(ibfr1, 'localDay', rdbdata1['day'])
    codes_set1(ibfr1, 'localHour', rdbdata1['hour'])
    codes_set1(ibfr1, 'localMinute', rdbdata1['minute'])
    codes_set1(ibfr1, 'localSecond', 0)

    codes_set1(ibfr1, 'section2Present', 1)

    codes_set1(ibfr1, 'rdbtimeDay', rdbdata1['day'])
    codes_set1(ibfr1, 'rdbtimeHour', rdbdata1['hour'])
    codes_set1(ibfr1, 'rdbtimeMinute', rdbdata1['minute'])
    codes_set1(ibfr1, 'rdbtimeSecond', 55)

    codes_set1(ibfr1, 'rectimeDay', rdbdata1['day'])
    codes_set1(ibfr1, 'rectimeHour', rdbdata1['hour'])
    codes_set1(ibfr1, 'rectimeMinute', rdbdata1['minute'])
    codes_set1(ibfr1, 'rectimeSecond', 55)

    codes_set1(ibfr1, 'restricted', rdbdata1['restricted'])

    codes_set1(ibfr1, 'correction1', 0)
    codes_set1(ibfr1, 'correction1Part', 0)
    codes_set1(ibfr1, 'correction2', 0)
    codes_set1(ibfr1, 'correction2Part', 0)
    codes_set1(ibfr1, 'correction3', 0)
    codes_set1(ibfr1, 'correction3Part', 0)
    codes_set1(ibfr1, 'correction4', 0)
    codes_set1(ibfr1, 'correction4Part', 0)
    codes_set1(ibfr1, 'qualityControl', 70)

    codes_set1(ibfr1, 'newSubtype', 0)

    codes_set1(ibfr1, 'localLatitude', rdbdata1['latitude'])
    codes_set1(ibfr1, 'localLongitude', rdbdata1['longitude'])
    codes_set1(ibfr1, 'ident', rdbdata1['ident'])


def check_value_before_bufr_encoding(data1,val_missing1,val_dtype1):
    if np.isnan(data1):
      return val_missing1
    elif val_dtype1==int:
      return int(data1)
    else:
      return data1

def encode_profile(codes_set1, ibfr1, level_kays1, profile_data1, nlevs1=None):

    data_kays = [] ; miss_kays = []
    for k in level_kays1:
      if k[0] in profile_data1.keys():
        data_kays.append(k)
      else:
        miss_kays.append(k)

    if nlevs1 is not None:
      for ilev1 in range(nlevs1):
        ient = '#{0}#'.format(ilev1+1)
        for k in data_kays:
          codes_set1(ibfr1, ient+ k[0], check_value_before_bufr_encoding( profile_data1[k[0]][ilev1], k[1], k[2] ) )
        for k in miss_kays:
          codes_set1(ibfr1, ient+ k[0], k[1] )
    else:
      for k in data_kays:
        codes_set1(ibfr1, k[0], check_value_before_bufr_encoding( profile_data1[k[0]], k[1], k[2] ) )
      for k in miss_kays:
        codes_set1(ibfr1, k[0], k[1] )

