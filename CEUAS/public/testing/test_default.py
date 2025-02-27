import sys,os
sys.path.append(os.path.expanduser('~/python/hug2'))
from default import *
from cds_eua3 import *
import io
import h5py

active,cdm,cf=init_server()
randdir='{:012d}'.format(100000000000)
dzip=randdir+'/download.zip'

def test_bbox1():
    body=eval('{"bbox":[40,0,50,20],"date":["19790101","19790131"],"variable":"temperature"}')
    ret=process_request(body,randdir,cdm)
    assert ret[0]=='100000000000/download.zip' and ret[1]==''
    
def test_bbox2():
    body=eval('{"bbox":[-40,140,80,359],"date":["19190101-20130131"],"variable":"temperature"}')
    ret=process_request(body,randdir,cdm)
    with zipfile.ZipFile(dzip,'r') as f:
        x=f.namelist()
        hf = io.BytesIO(f.read(x[0]))
        h = h5py.File(hf,'r')
        assert 'ta' in h.keys()
        
    #assert 'dest_91938_air_temperature.nc' in x and 'dest_47991_air_temperature.nc' in x
    
def test_country(statstr):
    body=eval('{"country":'+statstr+',"date":["19190101-20090131"],"variable":"temperature"}')
    ret=process_request(body,randdir,cdm)
    assert ret=='100000000000/download.zip' #and ret[1]==''

def test_statid(statstr):
    if statstr=='all':
        evs='{"date":["19790101-19790131"],"variable":"temperature"}'
    else:
        evs='{"statid":'+statstr+',"date":["19190101-20190131"],"variable":"temperature"}'
        
    body=eval(evs)
    ret=process_request(body,randdir,cdm)
    assert ret=='100000000000/download.zip' # and ret[1]==''
    

def test_date(statstr):
    if statstr=='':
        evs='{"statid":"01001","variable":"temperature"}'
    else:
        evs='{"statid":"01001","date":'+statstr+',"variable":"temperature"}'
        
    body=eval(evs)
    ret=process_request(body,randdir,cdm)
    assert ret=='100000000000/download.zip'# and ret[1]==''

def test_dateplot(statstr):
    evs='{"statid":"10393","pressure_level":'+statstr+',"variable":"temperature","fbstats":["obs_minus_bg","bias_estimate"]}'
    body=eval(evs)
    ret=process_request(body,randdir,cdm)
    readandplot_ts('100000000000/download.zip',eval(evs))
    
def test_dateplot2(statstr):
    evs='{"statid":"72764","pressure_level":'+statstr+',"variable":"dew_point_temperature","fbstats":["obs_minus_bg","bias_estimate"]}'
    body=eval(evs)
    ret=process_request(body,randdir,cdm)
    readandplot_ts('100000000000/download.zip',eval(evs))
    
def test_plevel(statstr):
    if statstr=='':
        evs='{"statid":"10393","date":[20190101],"time":[5],"variable":"temperature"}'
    else: 
        evs='{"statid":"10393","date":[20190101],"time":[5],"pressure_level":'+statstr+',"variable":"temperature"}'
    body=eval(evs)
    ret=process_request(body,randdir,cdm)
    assert ret=='100000000000/download.zip' #and ret[1]==''
    
def test_variable(statstr,fbstr):
    if statstr=='':
        evs='{"statid":"10393","date":[20190101],"time":[5]}'
    else: 
        evs='{"statid":"10393","date":[20190101],"time":[5],"variable":'+statstr+'}'
    if fbstr!='':
        evs=evs[:-1]+',"fbstats":'+fbstr+'}'

    body=eval(evs)
    ret=process_request(body,randdir,cdm)
    assert ret=='100000000000/download.zip'# and ret[1]==''
    #readandplot('100000000000/download.zip',eval(evs))
    
def test_time(statstr,datestr):
    if statstr=='':
        evs='{"statid":"10393","date":'+datestr+',"variable":"temperature"}' 
    else:
        evs='{"statid":"10393","date":'+datestr+',"time":'+statstr+',"variable":"temperature"}' 
        
    body=eval(evs)
    ret=process_request(body,randdir,cdm)
    assert ret=='100000000000/download.zip' #and ret[1]==''
    

if __name__ == '__main__':

    #test for statid
    #for s in ["1001",'"01001"','["01001"]','[1001,1025]','["01001","01025"]','"all"']:
        #test_statid(s)

    test_bbox2()
    for s in ['["NOR","DEU"]','"NOR"','"ALL"']: #,'"XXX"']:
        test_country(s)
        
    for s in ['19571001','"19571001"','["19571001"]','[19571001,19571031]','["19571001-19571001"]','"19190101-20190101"','["19571001","19571031"]']: #,'"XXX"']:
        test_date(s)
    
    # test for "date" and "time"
    for d in ['["20140105"]','["20140102-20150105"]','["20130101","20140104","20140105","20140107"]']   :
        for s in ['[22,23,0,1,2,3]','["22","23","00","01","02","03"]','','0','"00"','["00"]','[0,12]','["0","12"]']: #,'"00-12"','"0-12"','"0-0"','"11-11"','["00-12"]','["21-03"]','["15-03"]','["09-15"]','["18-06"]','["6-18"]']:
            test_time(s,d)
    
    # test for "pressure_level"        
    for d in ['','10000','"10000"','[10000]','["10000"]','[10000,20000]','["10000","20000"]']:
        test_plevel(d)
    
    #test_dateplot2('50000')
    #test_dateplot('10000')
    
    # test for "variable" and "fbstats"
    for d in ['"temperature"','["temperature"]','["temperature","u_component_of_wind","v_component_of_wind","wind_speed","wind_direction","relative_humidity","specific_humidity","dew_point_temperature"]']:
        for fb in ['','"obs_minus_an"','["obs_minus_an"]','["obs_minus_bg","bias_estimate"]']:
            test_variable(d,fb)
        
    print('ready')