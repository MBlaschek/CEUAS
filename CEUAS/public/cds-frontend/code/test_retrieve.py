# Copernicus in situ early upper air data
# frontend test script
#
# 5 June 2020
#

import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
import cdsapi
import time

tt=time.time()

c = cdsapi.Client(url='https://sis-dev.climate.copernicus.eu/api/v2',
                  key='please add your credentials here',
                  progress=True,retry_max=5)

l=0

#Comprehensive upper-air observation network from 1905 to present
#    'insitu-comprehensive-upper-air-observation-network',
'''
'''
# variables and feedback information
#
# dewpoint temperature not yet working
#
for variable in ["air_temperature",["air_temperature"],["zonal_wind","meridional_wind","wind_speed","wind_direction","air_relative_humidity","air_specific_humidity"]]:
        for fb in ['',"obs_minus_an",["obs_minus_an"],["obs_minus_bg","bias_estimate"]]:

            print(variable,fb)
            c.retrieve(
                'insitu-comprehensive-upper-air-observation-network',
                {
                    'variable': variable,
                    'date':['199602{:0>2}'.format(i) for i in range(1,30)],
                    'statid':'10393',
                    'pressure_level':'50000',
                    'fbstats':fb,
                },
                'download{:0>3}.zip'.format(l))
            l+=1

'''
'''
for d in ['',["20140105"],["20140102-20150105"],["20130101","20140104","20140105","20140107"]]   :
        if d=='':
           slist=['',"00"]
        else:
           slist= [[22,23,0,1,2,3],["22","23","00","01","02","03"],'','0',"00",["00"],[0,12],["0","12"],"00-12","0-12","0-0","11-11",["00-12"],["21-03"],["15-03"],["09-15"],["18-06"],["6-18"]]
       
        for s in slist:
        
            print(d,s)
            c.retrieve(
                'insitu-comprehensive-upper-air-observation-network',
                {
                    'variable': 'air_temperature',
                    'date':d,
                    'time':s,
                    'statid':'10393',
                    'pressure_level':'70000',
                },
                'download{:0>3}.zip'.format(l))
            l+=1
'''
'''
for d in ['','10000',"10000",[10000],["10000"],[10000,20000],["10000","20000"]]:
    for f in ['','csv','nc',['csv','nc']]:
            print(d,f)
            c.retrieve(
                'insitu-comprehensive-upper-air-observation-network',
                {
                    'variable': 'air_temperature',
                    'date':'19950101-19951231',
                    'time':['23','00','01'],
                    'statid':'10393',
                    'pressure_level':d,
                    'format':f,
                },
                'download{:0>3}.zip'.format(l))
            l+=1
'''
'''
# selecting radiosonde via bbox, country, (WMO) statid
#
# selection via WIGOS radiosonde type not yet working
#
for d in [
          {'statid':['0-20000-0-10393']},
          {'statid':'10393'},
          {'statid':'all'},
          {'country':["NOR","DEU"]},
          {'country':"NOR"},
          {'country':"ALL"},
          {"bbox":[60,0,40,20]},
          {"bbox":[60,0,-40.5,20]},
          {"bbox":[60,140,-40.5,280]},
          {"bbox":[90,-180,-90,180]},

          ]:
            print(d)
            c.retrieve(
                'insitu-comprehensive-upper-air-observation-network',
                dict(
                    variable= 'air_temperature',
                    date='19190101-20150131',
                    pressure_level=50000,
                    **d
                ),
                'download{:0>3}.zip'.format(l))
            l+=1

print (time.time()-tt,' seconds')

# I tried with year month day but failed because the frontend adds a zero to the year. The request below generates 199500101
#        'year': '1995',
#        'month':['{:0>2}'.format(i) for i in range(1,3)],
#        'day':['{:0>2}'.format(i) for i in range(1,10)],

# I tried dewpoint as variable. This does not seem a valid variable. 
