import numpy
import time
import datetime
import netCDF4
import matplotlib.pylab as plt
import os,sys,glob
sys.path.append("/fio/srvx7/leo/python/Rasotools/") # containing the functions/utilities

from rasotools.utils import *
from multiprocessing import Pool
#import odb

from eccodes import *
from functools import partial
from collections import OrderedDict
import subprocess
import json
import gzip
sys.path.append(".")
from retrieve_fb_jra55 import add_feedback
import copy


from collections import OrderedDict                                                     
import subprocess                                                                       
import json                                                                            
import gzip                                                                             
from retrieve_fb_jra55 import add_feedback                                              
import copy                                                                               

# /opt/anaconda3/bin/python project custom python executable
# PYTHONPATH=/usr/local:/fio/srvx7/leo/python/namelist_python-master:.:/fio/srvx7/leo/python/Rasotools


sys.path.append("/opt/anaconda3/lib/python3.7")
from rasotools.utils import *


plt.rcParams['lines.linewidth'] = 3

def stretrieve(countlist,mindays):

    stdict=OrderedDict()
    for c in countlist:
        with open(c) as cf:
            cl=cf.read().split('\n')
            for zeile in cl[1:]:
                sp=zeile.split("\t")
                try:
                    stdict[sp[0][1:6]]+=numpy.int(float(sp[1]))
                except:
                    try:
                        stdict[sp[0][1:6]]=numpy.int(float(sp[1]))
                    except:
                        pass
    for k in list(stdict.keys()):
        if stdict[k]<mindays:
            stdict.pop(k,None)

    return stdict

def read_odb(fn):

    conn = odb.connect('')
    c = conn.cursor()
    qs="select date, time, obstype, codetype, sonde_type, statid, lat, lon, stalt, vertco_reference_1, "
    qs+="varno, obsvalue@body, biascorr@body"#, fg_depar@body, an_depar@body "  
    qs+="WHERE obstype=5 and codetype=35 and (varno=2  or varno=3 or varno=4 or varno=29 or varno=7 ) "
    #qs+="and (vertco_reference_1=100000 or vertco_reference_1=92500 or vertco_reference_1=85000 or vertco_reference_1=70000 or vertco_reference_1=50000 or vertco_reference_1=40000 or vertco_reference_1=30000 or vertco_reference_1=25000 or vertco_reference_1=20000 or vertco_reference_1=15000 or vertco_reference_1=10000 or vertco_reference_1=7000 or vertco_reference_1=5000 or vertco_reference_1=3000 or vertco_reference_1=2000 or vertco_reference_1=1000) "
    #qs+="and (time >={} or time <={}) "
    #qs+="order by statid,date,time,vertco_reference_1"
    qs+=" from '"+fn+"';"
#    qs+="select date, time from '"+fn+"';"

    c.execute(qs)

    rdata=list(c.fetchall())
    print(rdata)

    return rdata

def read_odbsql(fn):

    qs="odb sql -q 'select date, time, obstype, codetype, sonde_type, statid, lat, lon, stalt, "
    qs+="vertco_reference_1, varno, obsvalue, biascorr, fg_depar@body, an_depar@body "
    qs+="WHERE obstype=5 and codetype=35 and (varno=2  or varno=3 or varno=4 or varno=29 or varno=7 ) "
    qs+="and (vertco_reference_1=100000 or vertco_reference_1=92500 or vertco_reference_1=85000 or vertco_reference_1=70000 or vertco_reference_1=50000 or vertco_reference_1=40000 or vertco_reference_1=30000 or vertco_reference_1=25000 or vertco_reference_1=20000 or vertco_reference_1=15000 or vertco_reference_1=10000 or vertco_reference_1=7000 or vertco_reference_1=5000 or vertco_reference_1=3000 or vertco_reference_1=2000 or vertco_reference_1=1000) "
    qs+="order by statid,date,time,vertco_reference_1' " 
    qs+="-i $SCRATCH/ectrans/era5/era5.conv."+fn[-6:]+"> /raid8/scratch/leo/scratch/Stationfiles_era5/"+fn[-6:]+".txt "

    if os.path.getsize(os.path.expandvars("$SCRATCH/ectrans/era5/era5.conv."+fn[-6:]))>100000:
        rdata=os.popen(qs).read()

    return 

def read_odbsql_stn(ipath,opath,statids):

    countlist=glob.glob(opath+'/*.count')
    alldata=''
    alldict=dict()
    for c in countlist:
        month=c[-12:-6]
        #toberead=False
        #with open(c) as cf:
            #cdata=cf.read().split('\n')
            #for zeile in cdata:
                #if statid in zeile:
                    #toberead=True
                    #break
#	if toberead:
        qs="select date, time, obstype, codetype, sonde_type, statid, int(lat*100), int(lon*100), int(stalt), "
        qs+="int(vertco_reference_1), int(varno), int(obsvalue*100), int(biascorr*100), int(fg_depar@body*100), int(an_depar@body*100) "
        qs+="WHERE varno=2"
        #and ("
        #st=""
        #for s in statids[:-1]:
            #st+="statid='"+s[0]+"' or "
        #st+="statid='"+statids[-1][0]+"'"
        #qs+=st+")
        qs+=" and (vertco_reference_1=100000 or vertco_reference_1=92500 or vertco_reference_1=85000 or vertco_reference_1=70000 or vertco_reference_1=50000 or vertco_reference_1=40000 or vertco_reference_1=30000 or vertco_reference_1=25000 or vertco_reference_1=20000 or vertco_reference_1=15000 or vertco_reference_1=10000 or vertco_reference_1=7000 or vertco_reference_1=5000 or vertco_reference_1=3000 or vertco_reference_1=2000 or vertco_reference_1=1000) "
        qs+="order by statid,date,time,vertco_reference_1" 

        t=time.time()
        if os.path.getsize(ipath+"/era5.conv."+month)>100000:
    #        rdata=os.popen(qs).read()
            rdata=subprocess.check_output(["odb","sql","-q",qs,"-i",ipath+"/era5.conv."+month])
            pass

        rdata=''.join(rdata.split(" "))
        rdata = rdata.split("'")
        strow=rdata[1::2]
        rdata = " ".join(rdata)
        rdata=rdata.split("NULL")
        rdata="NaN".join(rdata)
        rdata=rdata.split('\n')
        print(time.time()-t)
        for k,v in zip(strow,rdata[1:-1]):
            try:
                alldict[k].append(v)
            except:
                alldict[k]=[v]


        #with open(opath+'/'+month+'.json','w') as f:
            #json.dump(alldict,f)
        print(c.split('/')[-1],time.time()-t)
#	alldata+=rdata[rdata.index('\n'):]

    #f=open(opath+'001001.txt','r')
    #alldata=f.read()
    #f.close()

    for statid,v in alldict.items():
        readstationfiles_t(opath,os.path.expandvars('$FSCRATCH/ei6/0')+statid,'ERA5_0'+statid,'\n'.join(v))

    return

def par_read_odbsql_stn_withfeedback(varno,odbfile):

    #countlist=glob.glob(opath+'/*.count')
    alldata=''
    alldict=dict()


    if 'presat' in odbfile or '.2402.' in odbfile:
        qs="select statid, source,date, time, obstype, codetype, sonde_type, lat, lon, stalt, "
    else:
        qs='select statid, date, time, obstype, codetype, sonde_type, lat, lon, lsm@modsurf, '
    qs+="vertco_type, vertco_reference_1, varno, obsvalue@body, biascorr@body, fg_depar@body, an_depar@body "
    qs+="WHERE "
#    qs+="varno=2  and "
    qs+=" (vertco_reference_1=100000 or vertco_reference_1=92500 or vertco_reference_1=85000 or vertco_reference_1=70000 or vertco_reference_1=50000 or vertco_reference_1=40000 or vertco_reference_1=30000 or vertco_reference_1=25000 or vertco_reference_1=20000 or vertco_reference_1=15000 or vertco_reference_1=10000 or vertco_reference_1=7000 or vertco_reference_1=5000 or vertco_reference_1=3000 or vertco_reference_1=2000 or vertco_reference_1=1000) "

    t=time.time()
    sonde_type=True
    obstype=True
    if os.path.getsize(odbfile)>0:
#        rdata=os.popen(qs).read()
        try:
            rdata=subprocess.check_output(["odb","sql","-q",qs,"-i",odbfile,'--no_alignment'])
        except subprocess.CalledProcessError as e:
            qs=''.join(qs.split('sonde_type,'))
            try:
                rdata=subprocess.check_output(["odb","sql","-q",qs,"-i",odbfile,'--no_alignment'])
                sonde_type=False
            except subprocess.CalledProcessError as e:
                qs=''.join(qs.split('obstype,'))
                qs=''.join(qs.split('codetype,'))
                qs=''.join(qs.split('an_depar@body'))
                qs=''.join(qs.split('lsm@modsurf,'))
                try:
                    rdata=subprocess.check_output(["odb","sql","-q",qs,"-i",odbfile,'--no_alignment'])
                    obstype=False
                except subprocess.CalledProcessError as e:
                    print('odb failed!:'+' '+odbfile+' '+qs)
                    return alldict
    else:
        return alldict
    #print 'sql',time.time()-t
    tx=time.time()
    #if len(rdata.split('\n'))<3:
        #return alldict

    #rdata=''.join(rdata.split(" "))
    rdata = rdata.decode().split("'")
    if 'presat' in odbfile or '.2402.' in odbfile:
        stationrow=rdata[1::4]
        if len(stationrow)==0:
            return alldict
        sourcerow=rdata[3::4]
        header=rdata[0].split('\n')[0]
        rdata = "".join(rdata[4::4])
        cols=len(header.split('\t'))-2
        alldict['header']=header.split()[2:]
    else:
        stationrow=rdata[1::2]
        if len(stationrow)==0:
            return alldict
        sourcerow=['-']*len(stationrow)
        header=rdata[0].split('\n')[0]
        rdata = "".join(rdata[2::2])
        cols=len(header.split('\t'))-1
        alldict['header']=header.split()[1:]

    for h in range(len(alldict['header'])):
        alldict['header'][h]=alldict['header'][h].split('@')[0]
    rdata=rdata.split("NULL")
    rdata="NaN".join(rdata)

    # print time.time()-tx
    xdata=numpy.fromstring(rdata,sep='\t',dtype=numpy.float64)
    rdata=[]

    xdata=numpy.reshape(xdata,(xdata.shape[0]//cols,cols))
    #print time.time()-tx

    statid=stationrow[0]
    alldict[statid]=dict()
    ad=xdata

    h=alldict['header']
    idy=numpy.lexsort((ad[:,h.index('varno')],
                       ad[:,h.index('vertco_reference_1')],
                       ad[:,h.index('time')],
                       ad[:,h.index('date')]))
    alldict[statid]['data']=xdata[idy,:]
    alldict[statid]['source']=[sourcerow[0]] 
    alldict[statid]['odbfile']=odbfile
    alldict[statid]['odbstatid']=[stationrow[0]] 

    print(odbfile,time.time()-tx)

    return alldict

def par_read_odbascii_stn_withfeedback(varno,odbfile):

    #countlist=glob.glob(opath+'/*.count')
    alldata=''
    alldict=dict()


    t=time.time()
    try:
        with gzip.open(odbfile, 'rb') as f:
            rdata = f.read()   
        sonde_type=True
        obstype=True
    except:
        print('odb failed!:'+' '+odbfile)
        return alldict

    alldict['header']=['date', 'time', 'obstype', 'codetype', 'sonde_type', 'statid','lat', 'lon', 'stalt', 'vertco_reference_1', 
                       'varno', 'obsvalue', 'biascorr','fg_depar', 'an_depar','vertco_type','y','z']
    print('gzip',time.time()-t)
    tx=time.time()

    rdata=rdata.split("NULL")
    rdata="NaN".join(rdata)


    # print time.time()-tx
    xdata=numpy.fromstring(rdata,sep=' ',dtype=numpy.float64)
#    rdata=[]

    #cols=len(alldict['header'])
    cols=len(rdata[:rdata.index('\n')].split())

    xdata=numpy.reshape(xdata,(xdata.shape[0]/cols,cols))
    if cols==17:
        xd=numpy.zeros((xdata.shape[0],18))
        xd[:,:4]=xdata[:,:4]
        xd[:,5:]=xdata[:,4:]
        xdata=xd

    #if cols==17:
        #print cols,xdata[0,:]

    #else:
        #print cols,xdata[0,:]

    #print time.time()-tx

    if varno==2:
        idy=numpy.where(xdata[:,10]==varno)
    elif varno==3:
        idy=numpy.where(numpy.logical_or(xdata[:,10]==3,xdata[:,10]==4))

    ad=xdata[idy]
    xdata=0
    ad[:,15]=1.0 # set vertco_type (1=pressure levels)
    statid=odbfile[-15:-9]
    alldict[statid]=dict()

    alldict[statid]['data']=ad[:,:16]
    alldict[statid]['source']=['BUFRDATA'] 
    alldict[statid]['odbfile']=odbfile
    alldict[statid]['odbstatid']=[statid]

    print(odbfile,time.time()-tx)

    return alldict


def par_read_odbsql_stn_nofeedback(varno,odbfile):

    alldata=''
    alldict=dict()


    qs="select statid, source, date, time, obstype, codetype, lat, lon, stalt, "
    qs+="vertco_type,vertco_reference_1, varno, obsvalue"
    if varno==2:
        qs+=" where varno={} ".format(varno)
    elif varno==110:
        qs+=" where varno=111 or varno=112 ".format(varno)
    elif varno==7:
        #qs+=" where varno=7 ".format(varno)
        qs+=" where varno={} ".format(varno)
        
        
    t=time.time()
    sonde_type=True
    if os.path.getsize(odbfile)>0:
        try:
            rdata=subprocess.check_output(["odb","sql","-q",qs,"-i",odbfile,'--no_alignment'])
            #rdata=subprocess.check_output(["odb sql -q '" +qs+"' -i "+odbfile+' | '+'sed'+" 's/ *//g'"],shell=True)
        except subprocess.CalledProcessError as e:
            qs=''.join(qs.split('sonde_type,'))
            try:
                rdata=subprocess.check_output(["odb","sql","-q",qs,"-i",odbfile,'--no_alignment'])
                #rdata=subprocess.check_output(["odb sql -q '" +qs+"' -i "+odbfile+' | '+'sed'+" 's/ *//g'"],shell=True)
                sonde_type=False
            except subprocess.CalledProcessError as e:
                print('odb failed!:'+' '+odbfile+' '+qs)
                return alldict
    else:
        return alldict
    print('alldict FF:', alldict)
    #print 'sql',time.time()-t
    tx=time.time()

    #rdata=''.join(rdata.split(" "))
    rdata = rdata.decode().split("'")
    stationrow=rdata[1::4]
    if len(stationrow)==0:
        return alldict
    sourcerow=rdata[3::4]
    header=rdata[0]
    rdata = "".join(rdata[::4][1:])
    rdata=rdata.split("NULL")
    rdata="NaN".join(rdata)
    # print time.time()-tx
    xdata=numpy.fromstring(rdata,sep='\t',dtype=numpy.float64)
    rdata=[]

    cols=len(header.split('\t'))-2
    xdata=numpy.reshape(xdata,(xdata.shape[0]//cols,cols))
    #print time.time()-tx
    alldict['header']=header.split()[2:]
    for h in range(len(alldict['header'])):
        alldict['header'][h]=alldict['header'][h].split('@')[0]
    # 'header': ['date', 'time', 'obstype', 'codetype', 'lat', 'lon', 'stalt', 'vertco_type', 'vertco_reference_1', 'varno', 'obsvalue']
#    statid=stationrow[0]
    statid=odbfile.split('.')[-1]
    alldict[statid]=dict()
    ad=xdata
    h=alldict['header']
    idy=numpy.lexsort((ad[:,h.index('varno')],
                       ad[:,h.index('vertco_reference_1')],
                       ad[:,h.index('time')],
                       ad[:,h.index('date')]))
    alldict[statid]['data']=xdata[idy,:]
    alldict[statid]['source']=[sourcerow[0]] 
    alldict[statid]['odbstatid']=[stationrow[0]] 
    alldict[statid]['odbfile']=odbfile

    print(odbfile,time.time()-tx)
    #print('FF alldicts keys are: **** ', alldict.keys() , alldict['_67775'].keys() ) # e.g. for station 67775
    
    """ alldict is a dictionary with two keys. 
        The first one is called 'header', and contains the name of the 11 variables 
        ['date', 'time', 'obstype', 'codetype', 'lat', 'lon', 'stalt', 'vertco_type', 'vertco_reference_1', 'varno', 'obsvalue']
        contained in the data.
        The second key is the number identifying the station id, e.g. '_67693'.
        The second key maps to a list of data, each with 11 entries.
        alldict['_67693']['data'][5] will show the entries for the 6th data array in the list,
        for the 11 variables. 
        
        """
    return alldict

def par_read_igra_stn_nofeedback(varno,odbfile):

    alldata=''
    alldict=dict()


    t=time.time()
    sonde_type=True
    if os.path.getsize(odbfile)>0:
        try:
            with open(odbfile) as f:
                rdata=f.read().split('\n')
                hdr=rdata[0].split()
                l=0
                idx=[]
                h=[]
                rdata.pop()
                for r in rdata:
                    if r[0]=='#':
                        idx.append(l)
                        h.append(numpy.fromstring(r[6:36]+r[54:71],sep=' ',dtype=numpy.int))
                    l+=1
            for l in idx:
                rdata[l]=' '*50
            rdata='\n'.join(rdata)
            rdata=' '.join(rdata.split('A'))
            rdata=' '.join(rdata.split('B'))
            xdata=numpy.fromstring(rdata,sep=' ',dtype=numpy.float64)
            print(xdata.shape)


        except KeyError:
            print('igra read failed!:'+' '+odbfile)
            return alldict
    else:
        return alldict
    tx=time.time()


    cols=9
    xdata=numpy.reshape(xdata,(xdata.shape[0]/cols,cols))
    ydata=numpy.empty((xdata.shape[0]*2,11),numpy.float64)
    k=fill_ydata(xdata,ydata,numpy.asarray(idx),numpy.asarray(h),varno)

    #print time.time()-tx
    alldict['header']=['date', 'time', 'obstype', 'codetype', 'lat', 'lon', 'stalt', 'vertco_type', 'vertco_reference_1', 'varno', 'obsvalue']

#    statid=stationrow[0]
    statid=hdr[0][6:12]
    alldict[statid]=dict()

    alldict[statid]['data']=ydata[:k,:]
    alldict[statid]['source']=[hdr[7]] 
    alldict[statid]['odbstatid']=[statid] 
    alldict[statid]['odbfile']=odbfile

    print(odbfile,time.time()-tx)

    
    return alldict

@njit
def fill_ydata(xdata,ydata,idx,h,varno):

    c=(4,7,8)
    ll=0    
    k=0
    for i in range(h.shape[0]):
        l=ll
        if len(h[i])>6:
            for z in range(h[i,6]-1,-1,-1):
                if xdata[l+z,0]==10:
#                    for v in range(len(varno)):
                    ydata[k,0]=h[i,1]*10000+h[i,2]*100+h[i,3]
                    if h[i,4]!=99:
                        ydata[k,1]=h[i,4]*10000 #+h[i][5]
                    else:
                        if h[i,5]<2400:
                            ydata[k,1]=(h[i,5]/100)*10000 #+h[i][5]
                        else:
                            print('spurious hour')


                    ydata[k,2]=5.0
                    ydata[k,3]=35.0
                    ydata[k,4]=h[i,7]/10000.
                    ydata[k,5]=h[i,8]/10000.
                    ydata[k,6]=numpy.NaN
                    ydata[k,7]=1.0
                    ydata[k,8]=xdata[l+z,2]
                    ydata[k,9]=varno
                    if varno==2:
                        if xdata[l+z,c[0]]!=-9999 and xdata[l+z,c[0]]!=-8888:
                            ydata[k,10]=xdata[l+z,c[0]]/10.+273.15
                            k+=1
                    elif varno==111:
                        if xdata[l+z,c[1]]!=-9999 and xdata[l+z,c[1]]!=-8888:
                            ydata[k,10]=xdata[l+z,c[1]]
                            k+=1
                        if xdata[l+z,c[2]]!=-9999 and xdata[l+z,c[2]]!=-8888:
                            ydata[k,10]=xdata[l+z,c[2]]/10.
                            k+=1
                ll+=1

    if varno==111:
        k-=1
    return k



def par_read_bufr_stn_nofeedback(varno,bufrfile):

    '''       [  1.95301010e+07,   1.30000000e+04,   5.00000000e+00,
          3.50000000e+01,   5.59300000e+01,   3.75200000e+01,
          1.87000000e+02,   1.00000000e+00,   3.34000000e+04,
          2.00000000e+00,   2.22550003e+02]])
alldict['header']
['date', 'time', 'obstype', 'codetype', 'lat', 'lon', 'stalt', 'vertco_type', 'vertco_reference_1', 'varno', 'obsvalue']
len(alldict['header'])
11
    '''

    alldata=''
    alldict=dict()

    bufrlist=[]
    tx=time.time()
    try:
        f = open(bufrfile)
        cnt = 0
        # loop over the messages in the file
        while 1:
            # get handle for message
            bufr = codes_bufr_new_from_file(f)
            if bufr is None:
                break
            # we need to instruct ecCodes to expand all the descriptors
            # i.e. unpack the data section
            codes_set(bufr, 'unpack', 1)
            # get all the timePeriods
            #iterid = codes_bufr_keys_iterator_new(bufr)

            # loop over the keys
            #if codes_get_array(bufr,'dataSubCategory')[0]!=101:
    ##            print codes_get_array(bufr,'dataSubCategory')[0]
                #codes_release(bufr)
                #continue
            #while codes_bufr_keys_iterator_next(iterid):

                ## print key name
                #keyname = codes_bufr_keys_iterator_get_name(iterid)
                #print keyname,codes_get_array(bufr,keyname)

            ## delete the key iterator
            #codes_bufr_keys_iterator_delete(iterid)

            datum = float('19'+codes_get_array(bufr, "typicalDate")[0][2:])
            timePeriod = float(codes_get_array(bufr, "typicalTime")[0])
            pressure = codes_get_array(bufr, "pressure")
    #        extendedVerticalSoundingSignificance = codes_get_array(bufr, "extendedVerticalSoundingSignificance")
    #        geopotentialHeight = codes_get_array(bufr, "nonCoordinateGeopotentialHeight")
    #        latitudeDisplacement = codes_get_array(bufr, "latitudeDisplacement")
    #        longitudeDisplacement = codes_get_array(bufr, "longitudeDisplacement")
            if varno==2:
                airTemperature = codes_get_array(bufr, "airTemperature")
            elif varno==111:
            #dewpointTemperature = codes_get_array(bufr, "dewpointTemperature")
                windDirection = codes_get_array(bufr, "windDirection")
                windSpeed = codes_get_array(bufr, "windSpeed")
            else:
                print('unimplemented varno',varno)
                return alldict
            if cnt==0:
                lat = codes_get(bufr, "latitude")
                lon = codes_get(bufr, "longitude")
                alt = float(codes_get(bufr, "heightOfStation"))
                blockNumber = codes_get(bufr, "blockNumber")
                stationNumber = codes_get(bufr, "stationNumber")

            codes_release(bufr)
            #print 'station %d%d' % (blockNumber,stationNumber)

            #print 'timePeriod pressure geopotentialHeight latitudeDisplacement longitudeDisplacement airTemperature windDirection windSpeed significance'
            miss_val=-1.e100
            if varno==2:
                for i in range(0,len(airTemperature)):
                    if airTemperature[i]!=miss_val:
                        line=numpy.asarray((datum,timePeriod,5.,35.,lat,lon,alt,1.,pressure[i],2.0,airTemperature[i]))
                        bufrlist.append(line)#print pressure[i],airTemperature[i],windDirection[i],windSpeed[i]
                        cnt += 1
            else:
                miss_val=-1.e100
                for i in range(0,len(windDirection)):
                    if windSpeed[i]!=miss_val and windDirection[i]!=2147483647:
                        line=numpy.asarray((datum,timePeriod,5.,35.,lat,lon,alt,1.,pressure[i],111.0,windDirection[i]))
                        bufrlist.append(line)#print pressure[i],airTemperature[i],windDirection[i],windSpeed[i]
                        line=numpy.asarray((datum,timePeriod,5.,35.,lat,lon,alt,1.,pressure[i],112.0,windSpeed[i]))
                        bufrlist.append(line)#print pressure[i],airTemperature[i],windDirection[i],windSpeed[i]
                        cnt += 1

        f.close()
#        print '/'.join(bufrfile.split('/')[-1:]),cnt,"messages",time.time()-tx
    except:

        try:
            codes_release(bufr)
        except:
            pass
        try:
            f.close()
        except:
            pass
        return alldict

    if len(bufrlist)==0:
        return alldict
    alldict['header']=['date', 'time', 'obstype', 'codetype', 'lat', 'lon', 'stalt', 'vertco_type', 'vertco_reference_1', 
                       'varno', 'obsvalue']

    ad=numpy.asarray(bufrlist)
    h=alldict['header']
    statid=str(blockNumber*1000+stationNumber)
    alldict[statid]=dict()
    idy=numpy.lexsort((ad[:,h.index('varno')],
                       ad[:,h.index('vertco_reference_1')],
                       ad[:,h.index('time')],
                       ad[:,h.index('date')]))
    alldict[statid]['data']=ad[idy,:]
    alldict[statid]['source']=['BUFRDATA'] 
    alldict[statid]['odbstatid']=[statid]
    alldict[statid]['odbfile']=bufrfile

    print('/'.join(bufrfile.split('/')[-1:]),statid,cnt,"messages",time.time()-tx)

    return alldict

def read_odbsql_count(opath,ifile):

    ipath='/'.join(ifile.split('/')[:-1])
    month=ifile[-6:]
    qs="odb sql -q 'select statid,count(*) as num where varno=2 and vertco_reference_1=10000 order by statid' " 
    qs+="-i "+ipath+'/era5.conv.'+month+" > "+opath+"/"+month+".count "

    if os.path.getsize(ipath+'/era5.conv.'+month)>100000:
        rdata=os.popen(qs).read()

    return 

def profplot(statid,plevs,ts,tfgs,tbcs):
#    return
    plt.figure(figsize=(9,4))
    plt.subplot(1,2,2)
    h1,=plt.semilogy(numpy.nanmean(tfgs[0,:,:]+tbcs[0,:,:],axis=1),plevs/100.,label='00')
    h2,=plt.semilogy(-numpy.nanmean(tfgs[1,:,:]+tbcs[1,:,:],axis=1),plevs/100.,label='12')
    plt.xlabel('fg-bc [K]')
#    h3,=plt.semilogy(numpy.nanmean(tfgs[1,:,:],axis=1),plevs/100.,label='fg12')
#    h4,=plt.semilogy(-numpy.nanmean(tbcs[1,:,:],axis=1),plevs/100.,label='bc12')
    h3,=plt.semilogy(numpy.nanmean(ts[1,:,:]-ts[0,:,:],axis=1),plevs/100.,label='obs12-obs00')
    plt.ylim([1000,10])
    plt.legend(handles=[h1,h2,h3])
    plt.subplot(1,2,1)
    h1,=plt.semilogy(numpy.nanmean(tfgs[0,:,:],axis=1),plevs/100.,label='fg00')
    h2,=plt.semilogy(-numpy.nanmean(tbcs[0,:,:],axis=1),plevs/100.,label='bc00')
    h3,=plt.semilogy(numpy.nanmean(tfgs[1,:,:],axis=1),plevs/100.,label='fg12')
    h4,=plt.semilogy(-numpy.nanmean(tbcs[1,:,:],axis=1),plevs/100.,label='bc12')
#    h3,=plt.plot(lmeans[:,2],plevs/100.,label='fg')
    plt.ylim([1000,10])
    plt.xlim([-4,4])
    plt.ylabel('p [hPa]')
    plt.xlabel('[K]')
    plt.legend(handles=[h1,h2,h3,h4])
    plt.title('{}, {}h'.format(statid,0))
    plt.savefig('profs_{}_{}.ps'.format(statid,0))
    plt.close()

@jit(nopython=True)
def doqc(pmin,pmax,ip,uwind=0.,vwind=0.,ufg_dep=0.,vfg_dep=0.,uan_dep=0.,van_dep=0.):

    removed=0
    kept=0
    for ipar in range(uwind.shape[0]):
        for it in range(uwind.shape[2]):
            h=ufg_dep[ipar,ip,it]
            if h==h:
                if h<pmin or h>pmax or vfg_dep[ipar,ip,it]<pmin or vfg_dep[ipar,ip,it]>pmax:
                    ufg_dep[ipar,ip,it]=numpy.nan  
                    vfg_dep[ipar,ip,it]=numpy.nan  
                    uan_dep[ipar,ip,it]=numpy.nan  
                    van_dep[ipar,ip,it]=numpy.nan  
                    uwind[ipar,ip,it]=numpy.nan  
                    vwind[ipar,ip,it]=numpy.nan
                    removed+=1
                else:
                    kept+=1

    return removed,kept


@jit(nopython=True)
def mkwind(u,v,firsts,ps,mask,hours,source,uwind=0.,vwind=0.,ufg_dep=0.,vfg_dep=0.,uan_dep=0.,van_dep=0.):

    for iline in range(u.shape[0]):

        day=numpy.int(u[iline,0])
        idx=firsts[(day/10000-1900)*12+(day%10000)/100-1]+day%100-2 # days since... not Fortran convention
        if u[iline,1]>=60000 and u[iline,1]<=180000:
            ihour=1
        else:
            ihour=0

        ip=0
        found=False
        while ps[ip]<u[iline,9]:
            ip+=1
            if ip>=ps.shape[0]:
                break
            if ps[ip]==u[iline,9]:
                wp=0.
                found=True
                break
            else: # perhaps standard pl between this and the preceding significant level
                if ip>1:
                    if u[iline-1,9]<ps[ip-1] and u[iline,9]>ps[ip-1]:
                        if u[iline-1,9]-1>0:
                            ip-=1
                            wp=numpy.log(ps[ip]/u[iline-1,9])
                            wm=numpy.log(u[iline,9]/ps[ip])
                            found=True
                            break

        if not found:
            continue
        if wp==0.:
            uwind[ihour,ip,idx]=u[iline,11]
            vwind[ihour,ip,idx]=v[iline,11]
            ufg_dep[ihour,ip,idx]=u[iline,13]
            vfg_dep[ihour,ip,idx]=v[iline,13]
            uan_dep[ihour,ip,idx]=u[iline,14]
            van_dep[ihour,ip,idx]=v[iline,14]
        else:
            uwind[ihour,ip,idx]=(wp*u[iline,11]+wm*u[iline-1,11])/(wp+wm)
            vwind[ihour,ip,idx]=(wp*v[iline,11]+wm*v[iline-1,11])/(wp+wm)
            ufg_dep[ihour,ip,idx]=(wp*u[iline,13]+wm*u[iline-1,13])/(wp+wm)
            vfg_dep[ihour,ip,idx]=(wp*v[iline,13]+wm*v[iline-1,13])/(wp+wm)
            uan_dep[ihour,ip,idx]=(wp*u[iline,14]+wm*u[iline-1,14])/(wp+wm)
            van_dep[ihour,ip,idx]=(wp*v[iline,14]+wm*v[iline-1,14])/(wp+wm)


        hours[ihour,idx]=numpy.int(u[iline,1])/10000
        source[ihour,idx]=numpy.int(u[iline,15])
        mask[idx]=True


#@jit(nopython=True)
def mktemp(u,firsts,ps,mask,hours,source,uwind=0.,uwindbc=0.,ufg_dep=0.,uan_dep=0.):

    for iline in range(u.shape[0]):

        day=numpy.int(u[iline,0])
        idx=firsts[(day/10000-1900)*12+(day%10000)/100-1]+day%100-2 # days since... not Fortran convention
        if u[iline,1]>=60000 and u[iline,1]<=180000:
            ihour=1
        else:
            ihour=0

        ip=0
        found=False
        while ps[ip]<u[iline,9]:
            ip+=1
            if ip>=ps.shape[0]:
                break
            if ps[ip]==u[iline,9]:
                wp=0.
                found=True
                break
            else: # perhaps standard pl between this and the preceding significant level
                if ip>1:
                    if u[iline-1,9]<ps[ip-1] and u[iline,9]>ps[ip-1]:
                        if u[iline-1,9]-1>0:
                            ip-=1
                            wp=numpy.log(ps[ip]/u[iline-1,9])
                            wm=numpy.log(u[iline,9]/ps[ip])
                            found=True
                            break

        if not found:
            continue
        if wp==0.:
            uwind[ihour,ip,idx]=u[iline,11]/100. # /100 reverses multiplication with 100 during odb retrieve
            uwindbc[ihour,ip,idx]=u[iline,12]/100.
            ufg_dep[ihour,ip,idx]=u[iline,13]/100.
            uan_dep[ihour,ip,idx]=u[iline,14]/100.
        else:
            uwind[ihour,ip,idx]=(wp*u[iline,11]+wm*u[iline-1,11])/(wp+wm)/100.
            uwindbc[ihour,ip,idx]=(wp*u[iline,14]+wm*u[iline-1,12])/(wp+wm)/100.	    
            ufg_dep[ihour,ip,idx]=(wp*u[iline,13]+wm*u[iline-1,13])/(wp+wm)/100.
            uan_dep[ihour,ip,idx]=(wp*u[iline,14]+wm*u[iline-1,14])/(wp+wm)/100.	    

        hours[ihour,idx]=numpy.int(u[iline,1])/10000
        if u.shape[1]>15:
            source[ihour,idx]=numpy.int(u[iline,15])
        else:
            source[ihour,idx]=0

        mask[idx]=True



def readstationfiles_t(header,ncpath,prefix,statdata,varlist=['t','u','v','dp','h']):
    """ Reading the netCDF template file and filling with the data from the odb
    
        Args:   
             statdata: data from the odb file
             varlist : list of variables
                                        t   : temperature,
                                        u,v : u and v wind components
                                        dp  : dew point
                                        h   : specific humidity """
    #print('The nc path is FF', ncpath )
    #print('The prefix is FF', prefix )
    #print('The statdata[0]is FF', statdata[0] )
    
    t1=time.time()
    statid=statdata[0].replace('_','')

    # print('Statdata is FF:', statdata)
    # ('_98426', {'data': None, 'source': ['BUFRDATA'], 'odbstatid': ['98426'], 'odbfile': '/raid60/s...',
    
    
    #hlist=header.split()[2:]
    
    # FF Not used??? 
    #fields=['date','time','lat','lon','alt','vertco_type','vertco_reference_1','varno','obsvalue','biascor','fg_depar','an_depar']

    if 'ai_bfr' in prefix:
        while len(statid)<6:
            statid='0'+statid

    metadata=dict(t=dict(varname='temperatures',bias='tbiascorr',fg_dep='tfg_depar',an_dep='tan_depar',
                         datum='tdatum',hours='thours',units='K',long_name='Temperature',sonde_type='tsonde_type' , valid_range=[0,300]),
                  u=dict(varname='uwind',bias='ubiascorr',fg_dep='ufg_depar',an_dep='uan_depar',
                         datum='udatum',hours='uhours',units='m/s',long_name='Westerly wind component',sonde_type='usonde_type', valid_range=[-100.,100]),
                  v=dict(varname='vwind',bias='vbiascorr',fg_dep='vfg_depar',an_dep='van_depar',
                         datum='vdatum',hours='vhours',units='m/s',long_name='Southerly wind component',sonde_type='vsonde_type', valid_range=[-100.,100]) ,
                  dp=dict(varname='dp',bias='dpbiascorr',fg_dep='dpfg_depar',an_dep='dpan_depar',
                         datum='dpdatum',hours='dphours',units='???',long_name='Dew Point',sonde_type='dpsonde_type', valid_range=[-100.,100]) ,
                  h=dict(varname='h',bias='hbiascorr',fg_dep='hfg_depar',an_dep='han_depar',
                         datum='hdatum',hours='hhours',units='???',long_name='Specific Humidity',sonde_type='hsonde_type', valid_range=[-100.,100])                   
                  
                  
                  )
  
    ftn='/home/srvx7/leo/fastscratch/ei6/001001/001001_t'+'.nc'# template file
    f = netCDF4.Dataset(ftn,"r")
    #print ('Testing the template file FF', f , f.ncattrs())
    #print ('Variables of the template FF', f.variables.keys()) #['lat', 'lon', 'alt', 'press', 'datum', 'hours', 'temperatures', 'fg_dep', 'bias', 'flags', 'an_dep', 's_type']
    f.set_auto_maskandscale(False)
    
    print('ncpath: FF', ncpath+'/'+statid)
    for t in varlist: 
        print ('trying with the variable FF', t)
        #print ('the statdata is FF:', statdata)
        #print ('the metadata is FF:', metadata)
        
        """ 
        statdata[1].keys() = 'data', 'source', 'odbstatid', 'odbfile', 'nmax', 'mtemperatures', 'lat', 'lon', 'alt', 'obstype', 'weights', 'mdatum', 'mhours', 
        'tidx', 'sonde_type', 'wdidx', 'zref', 'tobsvalue', 'uobsvalue', 'vobsvalue', 'tbiascorr', 'ubiascorr', 'vbiascorr', 
        'tfg_depar', 'ufg_depar', 'vfg_depar', 'tan_depar', 'uan_depar', 'van_depar', 'temperatures', 'tdatum', 'tsonde_type', 'thours', 'tref']
        """
        try:
            if statdata[1][metadata[t]['datum']].shape[0]==0: # e.g. t = 't' , look for statdata[1]['tdatum'] where 'tdatum'=metadata['t']['datum'] 
                print(statid,metadata[t]['datum'], 'no data')
                continue
        except: # FF uncomment
            print('running in the exception FF ***')
            continue
        #except: 
        #    print('except')
        #    print('statdata[1] FF:', statdata[1])
        print('FF going through the try-except FF', )
        if not os.path.exists(ncpath+'/'+statid):
            os.makedirs(ncpath+'/'+statid)
        fno=ncpath+'/'+statid+'/'+prefix+statid+'_'+t+'.nc'
        print('the path to the output file is FF:', fno )
        fo = netCDF4.Dataset(fno,"w")
        print('written the new file FF', fno , fo.variables.keys())
        '''read and write global attributes'''
        for i in f.ncattrs(): # FF f.ncattrs=['Conventions', 'title', 'institution', 'history', 'source', 'references']
            if i=='history':
                setattr(fo,i,datetime.date.today().strftime("%Y/%m/%d"))
            elif i=='source':
                print('*** statdata[1][odbfile] FF', statdata[1]['odbfile'])
                setattr(fo,i,'Observation feedback, expID: '+statdata[1]['odbfile'].split('/')[-1] )
            elif i=='title':
                setattr(fo,i,'Station daily '+metadata[t]['varname']+' series' )
            else:
                setattr(fo,i,getattr(f,i))
                
        '''read and write dimensions of variables'''       
        tdim=statdata[1][metadata[t]['datum']].shape[0]
        #print('dimension keys of the template file FF:', f.dimensions.keys())
        for i in list(f.dimensions.keys()): # ['station', 'numdat', 'time', 'pressure', 'hour']
            if i=='time':
                fo.createDimension(i,tdim)
                if tdim !=statdata[1][metadata[t]['varname']].shape[2]:
                    print('DATUM DIMENSION NOT EQUAL NUMBER OF DAYS IN OBS')
                    continue

            else:
                try:
                    fo.createDimension(i,len(f.dimensions[i]))
                except:
                    flag=True
                    continue

        '''read and write data of variables'''     
        for i in list(f.variables.keys()): # f.variables.keys()= odict_keys(['lat', 'lon', 'alt', 'press', 'datum', 'hours', 'temperatures', 'fg_dep', 'bias', 'flags', 'an_dep', 's_type'])
            var=f.variables[i]          
            if i=='datum':
                fo.createVariable(i,var.dtype,var.dimensions)
                fo.variables[i][:]=statdata[1][metadata[t]['datum']][:]
            elif i=='lat':
                fo.createVariable(i,var.dtype,f.variables['datum'].dimensions)
                vhilf=numpy.asarray([statdata[1][i][0]]*tdim)
                fo.variables[i][:]=vhilf
            elif i=='lon':
                fo.createVariable(i,var.dtype,f.variables['datum'].dimensions)
                vhilf=numpy.asarray([statdata[1][i][0]]*tdim)
                fo.variables[i][:]=vhilf
            elif i=='alt':
                fo.createVariable(i,var.dtype,f.variables['datum'].dimensions)
                vhilf=numpy.asarray([statdata[1][i][0]]*tdim)
                fo.variables[i][:]=vhilf
            elif i=='hours':
                fo.createVariable(i,var.dtype,var.dimensions)
                fo.variables[i][:]=statdata[1][metadata[t]['hours']][:]
                fo.createVariable('source',str,var.dimensions)
                vhilf=numpy.asarray([[statdata[1]['source'][0]]*tdim,[statdata[1]['source'][0]]*tdim])
                fo.variables['source'][:]=vhilf
                fo.createVariable('odbstatid',str,('station'))
                vhilf=numpy.asarray(statdata[1]['odbstatid'])
                fo.variables['odbstatid'][:]=vhilf
            elif i=='temperatures':
                #print('writing the variable into temp, FF')
                #print('metadata[t][varname]:', metadata[t]['varname'])
                fo.createVariable(metadata[t]['varname'],var.dtype,var.dimensions,fill_value=numpy.nan)
                fo.variables[metadata[t]['varname']][:]=statdata[1][metadata[t]['varname']][:]
            elif i=='fg_dep':
                fo.createVariable(i,var.dtype,var.dimensions,fill_value=numpy.nan)
                if statdata[1][metadata[t]['fg_dep']].size>1:
                    fo.variables[i][:]=statdata[1][metadata[t]['fg_dep']][:]
            elif i=='an_dep':
                fo.createVariable(i,var.dtype,var.dimensions,fill_value=numpy.nan)
                if statdata[1][metadata[t]['an_dep']].size>1:
                    fo.variables[i][:]=statdata[1][metadata[t]['an_dep']][:]
            elif i=='bias':
                fo.createVariable(i,var.dtype,var.dimensions,fill_value=numpy.nan)
                if statdata[1][metadata[t]['bias']].size>1:
                    fo.variables[i][:]=statdata[1][metadata[t]['bias']][:]
            elif i=='press':
                fo.createVariable(i,var.dtype,var.dimensions)
                fo.variables[i][:]=var[:]
            elif i=='s_type':
                fo.createVariable(i,var.dtype,var.dimensions)
                fo.variables[i][:]=statdata[1][metadata[t]['sonde_type']][:]
            else:
                if i not in ('flags','lat','lon','alt'):
                    fo.createVariable(i,var.dtype,var.dimensions)

            '''read and write data of variables''' 
            iph= i in ['temperatures','fg_dep','an_dep','bias']
            for j in var.ncattrs(): # loop over the attributes of the data (read from template) and replace them 
                if j!='_FillValue' and j!='scale_factor' and j!='add_offset':
                    if i!='flags':
                        io=i
                        if i=='temperatures':
                            io=metadata[t]['varname'] # remember, t can be any variable e.g. uwind, but it will fill the field 'temperature' from the template
                        setattr(fo.variables[io],j,getattr(var,j))
                        if i=='temperatures' and j=='valid_range':
                            io=metadata[t]['varname']
                        #print('Set valid range to', metadata[t]['valid_range'])    
                        setattr(fo.variables[io],j,metadata[t]['valid_range'])
                        
                        if i=='temperatures' and j=='long_name':
                            setattr(fo.variables[io],j,metadata[t]['long_name'])				
                        if iph and j=='units':
                            setattr(fo.variables[io],j,metadata[t]['units'])				
                        if j=='missing_value':
                            if fo.variables[io].dtype in (numpy.dtype('float32'),numpy.dtype('float64')):
                                setattr(fo.variables[io],j,numpy.nan)
                            if fo.variables[io].dtype in (numpy.dtype('int32'),numpy.dtype('int64')):
                                setattr(fo.variables[io],j,-999)

                        #if i=='bias':
                            #setattr(fo.variables[i],j,getattr(var,j))
                        if i=='datum' and j=='units':
                            setattr(fo.variables[i],j,'days since 1900-01-01 00:00:00')
#	    setattr(fo.variables['source'],'missing_value','unknown')
        setattr(fo.variables['source'],'long_name','ODB2 source data set code')


        fo.close()
    f.close()
    print('wrote '+statid+':',time.time()-t1)


def topressure(gribpath,ps,tidx,stats):
    
    if not stats:
        return stats
    t=time.time()
    month=[]
    h=stats['header'] # FF ['date', 'time', 'obstype', 'codetype', 'lat', 'lon', 'stalt', 'vertco_type', 'vertco_reference_1', 'varno', 'obsvalue']
    hl=list()
    fgdepar=False
    for i in ['vertco_type','vertco_reference_1','varno','obsvalue','biascorr','fg_depar','an_depar','sonde_type']:
        try:
            hl.append(h.index(i))
            if i=='fg_depar':
                fgdepar=True
        except:
            hl.append(-1)

    for statid in list(stats.keys()):
        if statid not in ['header','odbfile']:
            s=stats[statid]
            if fgdepar:
                windno=3
            else:
                windno=111
            nmax=do_shape(tidx, s['data'],hl,varnos=(2,7,windno))
            s['nmax']=nmax+1 #len(numpy.unique(s['data'][:,0]))
            nmax=s['nmax']
            if nmax<2:
                continue
            for vn in ['mtemperatures']: #,'tref','uref','vref']:
                s[vn]=numpy.empty((2,ps.shape[0],nmax))
                s[vn][:]=numpy.nan

            s['lat']=[s['data'][0,h.index('lat')]]
            s['lon']=[s['data'][0,h.index('lon')]]
            s['alt']=[s['data'][0,h.index('lon')+1]]
            if s['lat']!=s['lat']:
                continue
            s['source']=[str(s['source'][0])]
            s['obstype']=s['data'][:,2]
            s['weights']=numpy.zeros(4)
            s['mdatum']=numpy.arange(nmax,dtype=numpy.int)
            s['mhours']=numpy.zeros((2,nmax),dtype=numpy.int)
            s['mhours'][:]=-999
            s['tidx']=numpy.where(s['data'][:,h.index('varno')]==2)[0]
            
            
            s['sonde_type']=numpy.zeros(nmax,dtype=numpy.int32)
            s['sonde_type'][:]=-999
            s['wdidx']=numpy.where(s['data'][:,h.index('varno')]==windno)[0]
            
            s['hidx']=numpy.where(s['data'][:,h.index('varno')]==7)[0] # FF humidity
            
            
            do_hours(tidx, s['mdatum'], s['mhours'],s['mtemperatures'], s['data'],s['sonde_type'], hl,varnos=(2,7,windno)) # ????
            
            if not fgdepar:
                for vn in ['zref']: #,'tref','uref','vref']:
                    nmax=s['mdatum'].shape[0]
                    s[vn]=numpy.empty((2,ps.shape[0],nmax))
                    s[vn][:]=numpy.nan
            else:
                s['zref']=numpy.empty((1,1,1))
            pass
            #print statid,s['mdatum'].shape
    if not fgdepar and gribpath is not None and 'ai_bfr' not in s['odbfile'] and 'igra' not in s['odbfile']:
        print('vor add',time.time()-t)
        tgrid=numpy.empty((124,16,181,360),dtype=numpy.float32)
        for vn in zip(['zref'],['129']): # ,'tref','uref','vref'],'130','131','132']):
            #add_feedback(gribpath,'erapreSAT{0}{1:02}.'+vn[1]+'.grb',
                            #vn[0],1949,1949+1,tidx,stats,fieldsperday=4,fcstep=0,tgrid=tgrid)#1930,2011
            if False:
                add_feedback(gribpath+'/../CERA20C_ens/','CERA20C{0}{1:02}.'+vn[1]+'.0.grb',
                             vn[0],1901,2010+1,tidx,stats,fieldsperday=4,fcstep=0,tgrid=tgrid)#1930,2011
            else:
                #add_feedback(gribpath+'/../CERA20C_ens/','CERA20C{0}{1:02}.'+vn[1]+'.0.grb',
                #             vn[0],1935,1935+1,tidx,stats,fieldsperday=4,fcstep=0,tgrid=tgrid)#1930,2011 # FF
                #add_feedback(gribpath+'/../CERA20C_ens/','CERA20C{0}{1:02}.'+vn[1]+'.0.grb',
                #             vn[0],1935,1979+1,tidx,stats,fieldsperday=4,fcstep=0,tgrid=tgrid)#1930,2011                
                add_feedback(gribpath+'/../CERA20C_ens/','CERA20C{0}{1:02}.'+vn[1]+'.0.grb',
                             vn[0],1980,1985,tidx,stats,fieldsperday=4,fcstep=0,tgrid=tgrid)#1930,2011 # FF

    for statid in list(stats.keys()):
        if statid not in ['header','odbfile']:
            s=stats[statid]
            if s['nmax']<2:
                continue
            for i in ['obsvalue','biascorr','fg_depar','an_depar']:
                for p in ['t','u','v','h','dp']: # FF added new variables h, dp
                    try:
                        s[p+i]=numpy.empty(s['mtemperatures'].shape)
                        s[p+i][:]=numpy.nan
                    except:
                        s[p+i]=numpy.empty((1,1,1))
                    pass
            #print statid
            
            if len(s['tidx'])>0:
                tgood=numpy.zeros(s['mdatum'].shape,dtype=numpy.bool)
                tgood=select_standardpressuredata(s['zref'],ps*100.,
                                                  s['mdatum'],s['data'],tgood,tidx,s['tidx'],2,hl,
                                                  s['tobsvalue'],
                                                  biascorr=s['tbiascorr'],
                                                  fg_depar=s['tfg_depar'],
                                                  an_depar=s['tan_depar'])
                if sum(tgood)<s['mdatum'].shape[0]:
                    for p in ['tobsvalue','tbiascorr','tfg_depar','tan_depar']:
                        try:
                            s[p]=s[p][:,:,tgood]
                        except:
                            pass
                s['temperatures']=s['tobsvalue']
                s['tdatum']=s['mdatum'][tgood]
                s['tsonde_type']=s['sonde_type'][tgood]
                s['thours']=s['mhours'][:,tgood]
                
                
                
            if len(s['hidx'])>0: # FF implemented new block for humidity h 
                tgood=numpy.zeros(s['mdatum'].shape,dtype=numpy.bool)
                tgood=select_standardpressuredata(s['zref'],ps*100.,
                                                  s['mdatum'],s['data'],tgood,tidx,s['hidx'],7,hl,
                                                  s['hobsvalue'],
                                                  biascorr=s['hbiascorr'],
                                                  fg_depar=s['hfg_depar'],
                                                  an_depar=s['han_depar'])
                if sum(tgood)<s['mdatum'].shape[0]:
                    for p in ['hobsvalue','hbiascorr','hfg_depar','han_depar']:
                        try:
                            s[p]=s[p][:,:,tgood]
                        except:
                            pass                    
                s['h']=s['hobsvalue']
                s['hdatum']=s['mdatum'][tgood]
                s['hsonde_type']=s['sonde_type'][tgood]
                s['hhours']=s['mhours'][:,tgood]
                    
                    
                
                
                
            if len(s['wdidx'])>0:
                ugood=numpy.zeros(s['mdatum'].shape,dtype=numpy.bool)
                ugood=select_standardpressuredata(s['zref'],ps*100.,
                                                  s['mdatum'],s['data'],ugood,tidx,s['wdidx'],windno,hl,                                                       
                                                  obs=s['uobsvalue'],
                                                  biascorr=s['ubiascorr'],
                                                  fg_depar=s['ufg_depar'],
                                                  an_depar=s['uan_depar'],
                                                  vobs=s['vobsvalue'],
                                                  vbiascorr=s['vbiascorr'],
                                                  vfg_depar=s['vfg_depar'],
                                                  van_depar=s['van_depar'])
                #if statid=='67197':
                    #print statid
                if sum(ugood)<s['mdatum'].shape[0]:
                    for p in ['uobsvalue','ubiascorr','ufg_depar','uan_depar','vobsvalue','vbiascorr','vfg_depar','van_depar']:
                        try:
                            s[p]=s[p][:,:,ugood]
                        except:
                            pass
                for vv in ['u','v']:
                    s[vv+'wind']=s[vv+'obsvalue']
                    s[vv+'datum']=s['mdatum'][ugood]
                    s[vv+'sonde_type']=s['sonde_type'][ugood]
                    s[vv+'hours']=s['mhours'][:,ugood]
            s['data']=None
            s['mtemperatures']=None
            s['tref']=None
            s['zref']=None

    print('nach select',time.time()-t)
    return stats
#@njit
def select_standardpressuredata(zref,ps,mdatum,data,good,tidx,didx,varno,hl,
                                obs,biascorr=numpy.empty((1,1,1)),
                                fg_depar=numpy.empty((1,1,1)),an_depar=numpy.empty((1,1,1)),
                                vobs=numpy.empty((1,1,1)),vbiascorr=numpy.empty((1,1,1)),
                                vfg_depar=numpy.empty((1,1,1)),van_depar=numpy.empty((1,1,1))
                                ):

    vday=0
    vhour=1
    vertco_type=hl[0]
    vcr=hl[1]
    vn=hl[2]

    #good=numpy.empty((mdatum.shape[0]),dtype=boolean)
    #good=numpy.empty((mdatum.shape[0]),dtype=numpy.bool)
    #good[:]=False

    irec=0
    m=-1
    idxold=-1
    for l in range(len(didx)):
        i=didx[l]
        if l>0:
            im1=didx[l-1]
        else:
            im1=didx[l]
        day=numpy.int(data[i,vday])
        year=day//10000
        if year<1900:
            print('Found an incorrect value for the year! Skipping', year)            
            continue
        if year> 2019:
            print('Found an incorrect value for the year! Skipping', year)
            continue 
        if data[i,vhour]>=60000 and data[i,vhour]<180000:
            ihour=1
            iday=0
        else:
            ihour=0
            iday=0
            if data[i,vhour]>=180000:
                iday=1

        idx=tidx[(year-1900)*12+(day%10000)//100-1]+day%100-1+iday#-tidx[(year-1900)*12] -1 # days since... not Fortran convention
        
        if idxold<idx:
            m+=1
            idxold=idx
            #if m>=mdatum.shape[0]: # can happen if launches are late and twice daily
                #ms=obs.shape
                #obs=numpy.concatenate((obs,numpy.zeros((ms[0],ms[1],1))),axis=2)
                #mdatum=numpy.concatenate((mdatum,numpy.zeros((1),dtype=numpy.int32)),axis=0)
                #good=numpy.concatenate((good,numpy.zeros((1),dtype=boolean)),axis=0)
                #mhours=numpy.concatenate((mhours,numpy.zeros((ms[0],1),dtype=numpy.int32)),axis=1)
                #sonde_type=numpy.concatenate((sonde_type,numpy.zeros((1),dtype=numpy.int32)),axis=0)
        if m<mdatum.shape[0]:
            mdatum[m]=idx
        #else:
            #print ('skip ',numpy.int(data[i,vday]))
            #continue
        if data[i,vertco_type]!=1: # height levels
            if data[i,vcr]!=data[i,vcr]:
                continue
            ip=0
            found=False
            while zref[ihour,ip,m]>data[i,vcr]:
                ip+=1
                if ip>=ps.shape[0]:
                    break
                if zref[ihour,ip,m]==data[i,vcr]:
                    wp=1.
                    wm=0.
                    found=True
                    break
                else: # perhaps standard pl between this and the preceding significant level
                    if ip>1:
                        if data[im1,vcr]<zref[ihour,ip,m] and data[i,vcr]>zref[ihour,ip,m]:
                            if data[im1,vcr]-1>0:
                                wp=numpy.log(zref[ihour,ip,m]/data[im1,vcr])
                                wm=numpy.log(data[i,vcr]/zref[ihour,ip,m])
                                found=True
                                break

            if not found:
                continue
            good[m]=True
            if data[i,vn]==111.0:
                if i >= len(data)-1 : continue # FF fixed for C:4843  
                u=numpy.cos((-data[i,hl[3]]+270.)*numpy.pi/180.)*data[i+1,hl[3]]
                v=numpy.sin((-data[i,hl[3]]+270.)*numpy.pi/180.)*data[i+1,hl[3]]
                um1=numpy.cos((-data[im1,hl[3]]+270.)*numpy.pi/180.)*data[im1+1,hl[3]]
                vm1=numpy.sin((-data[im1,hl[3]]+270.)*numpy.pi/180.)*data[im1+1,hl[3]]
                obs[ihour,ip,m]=(wp*u+wm*um1)/(wp+wm)
                vobs[ihour,ip,m]=(wp*v+wm*vm1)/(wp+wm)
            else:
                obs[ihour,ip,m]=(wp*data[i,hl[3]]+wm*data[im1,hl[3]])/(wp+wm)


        else:
            ip=0
            found=False
            while ps[ip]<data[i,vcr]:
                ip+=1
                if ip>=ps.shape[0]:
                    break

            if ip<ps.shape[0]:
                if ps[ip]==data[i,vcr]:
                    wp=1.
                    wm=0.
                    found=True
                else: # perhaps standard pl between this and the preceding significant level
                    if ip>1:

                        if data[im1,vcr]<ps[ip-1] and data[i,vcr]>ps[ip-1]:
                            if data[im1,vcr]-1>0:
                                ip-=1
                                wp=numpy.log(ps[ip]/data[im1,vcr])
                                wm=numpy.log(data[i,vcr]/ps[ip])
                                found=True

            if not found:
                continue
            good[m]=True
            if data[i,vn]==111.0:
                u=numpy.cos((-data[i,hl[3]]+270.)*numpy.pi/180.)*data[i+1,hl[3]]
                v=numpy.sin((-data[i,hl[3]]+270.)*numpy.pi/180.)*data[i+1,hl[3]]
                um1=numpy.cos((-data[im1,hl[3]]+270.)*numpy.pi/180.)*data[im1+1,hl[3]]
                vm1=numpy.sin((-data[im1,hl[3]]+270.)*numpy.pi/180.)*data[im1+1,hl[3]]
                obs[ihour,ip,m]=(wp*u+wm*um1)/(wp+wm)
                vobs[ihour,ip,m]=(wp*v+wm*vm1)/(wp+wm)
            elif data[i,vn]==3.0:
                obs[ihour,ip,m]=(wp*data[i,hl[3]]+wm*data[im1,hl[3]])/(wp+wm)
                vobs[ihour,ip,m]=(wp*data[i+1,hl[3]]+wm*data[im1+1,hl[3]])/(wp+wm)
                if hl[4]>-1:
                    biascorr[ihour,ip,m]=(wp*data[i,hl[4]]+wm*data[im1,hl[4]])/(wp+wm)
                    vbiascorr[ihour,ip,m]=(wp*data[i+1,hl[4]]+wm*data[im1+1,hl[4]])/(wp+wm)
                    fg_depar[ihour,ip,m]=(wp*data[i,hl[5]]+wm*data[im1,hl[5]])/(wp+wm)
                    vfg_depar[ihour,ip,m]=(wp*data[i+1,hl[5]]+wm*data[im1+1,hl[5]])/(wp+wm)
                    an_depar[ihour,ip,m]=(wp*data[i,hl[6]]+wm*data[im1,hl[6]])/(wp+wm)
                    van_depar[ihour,ip,m]=(wp*data[i+1,hl[6]]+wm*data[im1+1,hl[6]])/(wp+wm)

            else:
                obs[ihour,ip,m]=(wp*data[i,hl[3]]+wm*data[im1,hl[3]])/(wp+wm)
                if hl[4]>-1:
                    biascorr[ihour,ip,m]=(wp*data[i,hl[4]]+wm*data[im1,hl[4]])/(wp+wm)
                    fg_depar[ihour,ip,m]=(wp*data[i,hl[5]]+wm*data[im1,hl[5]])/(wp+wm)
                    an_depar[ihour,ip,m]=(wp*data[i,hl[6]]+wm*data[im1,hl[6]])/(wp+wm)


    return good

@njit
def do_shape(tidx,data,hl,varnos):

    m=-1
    ixold=-1
    for i in range(data.shape[0]):
        if data[i,hl[2]] not in varnos:
            continue
        md=int(data[i,0])
        year=md//10000
        if year<1900:
            continue
        month=md%10000//100
        day=md%100
        hour=int(data[i,1])//10000
        if hour>=6 and hour<18:
            ih=1
            iday=0
        else:
            ih=0
            iday=0
            if hour>=18:
                iday=1 # date shift

        ix=tidx[(year-1900)*12+month-1]+day-1+iday#-tidx[(year-1900)*12]	
        if ixold<ix:
            m+=1
            ixold=ix
    return m

@njit
def do_hours(tidx,mdatum,mhours,mtemperatures,data,sonde_type,hl,varnos):

    m=-1
    ixold=-1
    for i in range(data.shape[0]):
        if data[i,hl[2]] not in varnos:
            continue
        md=int(data[i,0])
        year=md//10000
        if year<1900:
            continue
        
        month=md%10000//100
        day=md%100
        hour=int(data[i,1])//10000
        if year==1980:
            a=0
        if hour>=6 and hour<18:
            ih=1
            iday=0
        else:
            ih=0
            iday=0
            if hour>=18:
                iday=1 # date shift

        ix=tidx[(year-1900)*12+month-1]+day-1+iday#-tidx[(year-1900)*12]	
        if ixold<ix:
            m+=1
            ixold=ix
        #if m>=mdatum.shape[0]: # can happen if launches are late and twice daily
            #ms=mtemperatures.shape
            #mtemperatures=numpy.concatenate((mtemperatures,numpy.zeros((ms[0],ms[1],1))),axis=2)
            #mdatum=numpy.concatenate((mdatum,numpy.zeros((1),dtype=numpy.int32)),axis=0)
            #mhours=numpy.concatenate((mhours,numpy.zeros((ms[0],1),dtype=numpy.int32)),axis=1)
            #sonde_type=numpy.concatenate((sonde_type,numpy.zeros((1),dtype=numpy.int32)),axis=0)

        mtemperatures[ih,:,m]=0.
        mdatum[m]=ix
        mhours[ih,m]=hour
        if hl[7]!=-1:
            if data[i,hl[7]]==data[i,hl[7]]:
                sonde_type[m]=int(data[i,hl[7]])

    return 

def odb2netcdf(gribpath,sodblist,varno,odbreader,idx,k):
    tidx=calcdays(19000101,(2023-1900)*12)-1    
    plevs=numpy.asarray([10,20,30,50,70,100,150,200,250,300,400,500,700,850,925,1000])

#    varno=2
#    func = partial(par_read_odbsql_stn_nofeedback,varno)
    #found=False
    #for s in sodblist[:idx[k+1]]:
        #if '002963' in s:
            #found=True
    #if not found:
        #return
        
    func = partial(odbreader,varno)
    alllist=list(map(func,sodblist[idx[k]:idx[k+1]]))
    
    alldicts=dict()
    while len(alllist)>0:
        a=alllist.pop()
        if len(list(a.keys()))==0:
            continue
        for b in list(a.keys()):
            if b=='header':
                alldicts[b]=a['header']
            else:
                statid=b
                alldicts[statid]=a[statid]
    referencevalues=topressure(gribpath,plevs,tidx,alldicts)

    if   varno==110:
        varlist=['u','v']   
    elif varno==2:
        varlist=['t']
    elif varno ==7:
        varlist = ['h']
        print('Processing the Specific Humidity FF ***')
        
        
    #global exp
    if exp == '1':
        print('FF the exp is', exp )
    #func = partial(readstationfiles_t,'header',os.path.expandvars('$FSCRATCH/ei6/'),'ERA5_'+exp+'_',varlist=varlist)
        #func = partial(readstationfiles_t,'header',os.path.expandvars('netCDF_'+exp+'/'+exp),'ERA5_'+exp+'_',varlist=varlist)
        func = partial(readstationfiles_t,'header',os.path.expandvars('1_humidity/'+exp),'ERA5_'+exp+'_',varlist=varlist)
        print('Done with FF', exp , varlist)
    else:
        func = partial(readstationfiles_t,'header',os.path.expandvars('DEBUG_netCDF_AllRes/'+exp), exp+'_',varlist=varlist)

    print('FF referencevalues.items()', referencevalues.items())
    list(map(func,list(referencevalues.items())))
    
    return



"""Dictionary mapping the name of the dataset to the storage directory
e.g. the dataset can be found at /raid60/scratch/leo/scratch/era5/odbs/1/ """

res_database = { '1'   : '/raid60/scratch/leo/scratch/era5/odbs/',
                 '1761': '/raid60/scratch/leo/scratch/era5/odbs/',
                 '1759': '/raid60/scratch/leo/scratch/era5/odbs/',
                 '3188': '/raid60/scratch/leo/scratch/era5/odbs/' 
                }

def run_converter(dataset='', single_stat= '', pool=1, varno=0, debug=False):
    """ Function converting the odb files to netCDF
    
        Args:
             dataset: number of the dataset, as defined in the res_database dictionary
             
             single_station: single station to be processed.
                             The absolute path to the 
                             If the varibale is defined, it will replace the dataset variable
                             and only the single station will be converted.
                             
             pull   : number of cores used in multithreading. 
                      If set to False, only one single process is executed
                      
             varno  : variable id number 
                      2:Temp, 111:wind dir,)
                      
             debug : if True, will print some additional information for debugging
             
        Example: 
             run_converter(dataset=['1','1759'], pull=False, varno=2)
             Will convert the datasets '1' and '1759', using a single core, for the Temperature variable
 
             run_converte(single_stat=/raid60/scratch/leo/scratch/era5/odbs/1/era5.conv.__01455, pool=False, varno=2 )
             Will convert the file specified
    """
    
    
    def process_list(sodblist, pool):
        """Looping over the files in the list
           Args:
                sodblist: the list of odb files to be processed
        """
        
        sodblist.sort(key=os.path.getsize) 
        s=numpy.cumsum(list(os.path.getsize(f) for f in sodblist))
        if debug: print('The sodblist is:', sodblist)
        stride=30
        idx=[0]
        k=1e9
        x=numpy.where(s>k)[0]
        while len(x)>0:
            k+=1e9
            idx.append(x[0])
            x=numpy.where(s>k)[0]
        idx2=list(range(len(idx)))  
        idx.append(len(s))
        if pool:
            p = Pool(pool)            
            if debug: print('Running using ', pool, ' cores.')
            if debug: print('The exp is:', exp)
            # TO DO: replace nofeedback with withfeedback for exp == 1
            func = partial(odb2netcdf,gribpath,sodblist,varno,par_read_odbsql_stn_nofeedback,idx)
            p.map(func,  idx2)
            if debug: print('Done with multithreading.')            
        else:   
            if debug: print('Running on single core.')            
            func = partial(odb2netcdf,gribpath,sodblist,varno,par_read_odbsql_stn_nofeedback,idx)
            list(map(func,  idx2))
            
    gribpath=os.path.expandvars('/raid60/scratch/leo/scratch/ERApreSAT/')

    sodblist = ''
    
    if dataset and not single_stat: # select all the stations in a database directory
        global exp
        exp = dataset
        path = res_database[dataset]
        ipath=os.path.expandvars(path + exp + '/')
        if debug: print('The selected database is:', ipath, ' the exp is', exp , ' d is ', dataset)
        
        
        # selecting different file names according to their naming conventions
        if dataset == '1761' or dataset == '1759': 
            sodblist=glob.glob(ipath+'/era5.????.conv.?:*[!.nc]') #ipath+'/era5.'+exp+'.conv.?:*[!.nc]'
        elif dataset == '3188':
            sodblist=glob.glob(ipath+'/era5.3188.conv.C:*[!.nc]') #ipath+'/era5.conv.C:4490*) 
        elif dataset == '1': 
            sodblist=glob.glob(ipath+'/*.conv._*')    
        
        '''
        # to process a subset of files
        #sodblist = [ el for el in sodblist if sodblist.index(el) in range(0,500)  ] # FF SLIMMER, processing only 5 files
        new = [ l.replace('/raid60/scratch/leo/scratch/era5/odbs/3188/era5.3188.conv.','') for l in sodblist ]
        processed = ['C:6146', 'C:6143', 'C:6155', 'C:6149', 'C:6158', 'C:6145', 'C:6157', 'C:6144', 'C:6154', 'C:6161', 'C:6164', 'C:6148', 'C:6168', 'C:6015A', 'C:6152', 'C:6160', 'C:6010A', 'C:6141', 'C:6014A', 'C:6137', 'C:6011A', 'C:6140', 'C:6159', 'C:6142', 'C:6136', 'C:6150', 'C:6162', 'C:6147', 'C:6132', 'C:6133', 'C:6167', 'C:6138', 'C:6165', 'C:6151', 'C:6130', 'C:6134', 'C:6156', 'C:6153', 'C:6131', 'C:6135', 'C:6139', 'C:6170', 'C:6012A', 'C:6169', 'C:6019A', 'C:6163', 'C:6067', 'C:6013A', 'C:6009A', 'C:5012A', 'C:6166', 'C:4529', 'C:4790', 'C:4528', 'C:4526', 'C:4527', 'C:4791', 'C:5853A', 'C:6017A', 'C:4525', 'C:5006A', 'C:4523', 'C:5856A', 'C:4524', 'C:4788', 'C:4793', 'C:5999', 'C:4530', 'C:4538', 'C:5012B', 'C:6003', 'C:6001', 'C:6018A', 'C:5003B', 'C:4963', 'C:4808', 'C:4806', 'C:6016A', 'C:4537', 'C:5991', 'C:6000', 'C:4536', 'C:4207A', 'C:5998', 'C:5782', 'C:5993', 'C:4535', 'C:5016B', 'C:6002', 'C:4704', 'C:5819', 'C:5995', 'C:4706', 'C:4534', 'C:4683', 'C:4663', 'C:4668', 'C:4672', 'C:4697', 'C:4731', 'C:5858A', 'C:5824', 'C:4532', 'C:6008', 'C:4703', 'C:4780', 'C:4740', 'C:4533', 'C:4976', 'C:5994', 'C:4531', 'C:4743', 'C:4755', 'C:5997', 'C:4750', 'C:4687', 'C:4775', 'C:4782', 'C:4779', 'C:4692', 'C:4778', 'C:5479', 'C:5992', 'C:4751', 'C:4721', 'C:4758', 'C:4702', 'C:4752', 'C:4774', 'C:4762', 'C:4616', 'C:5527', 'C:4207B', 'C:4701', 'C:4700', 'C:4681', 'C:4688', 'C:4667', 'C:4748', 'C:4655', 'C:4760', 'C:4629', 'C:4557', 'C:4696', 'C:5808', 'C:4781', 'C:5821', 'C:4665', 'C:4749', 'C:4759', 'C:4756', 'C:4565', 'C:4691', 'C:4605', 'C:5250A', 'C:4773', 'C:4567', 'C:4614', 'C:4678', 'C:4636', 'C:4662', 'C:4651', 'C:4563', 'C:4671', 'C:4627', 'C:4574', 'C:4611', 'C:4543', 'C:4604', 'C:4568', 'C:4552', 'C:4549', 'C:4595', 'C:4764', 'C:4637', 'C:4626', 'C:4577', 'C:4619', 'C:4625', 'C:4583', 'C:4682', 'C:4623', 'C:4640', 'C:4588', 'C:5004B', 'C:4674', 'C:4585', 'C:4630', 'C:4972', 'C:4561', 'C:4584', 'C:4542', 'C:4650', 'C:4569', 'C:4757', 'C:4649', 'C:4556', 'C:4564', 'C:4679', 'C:4657', 'C:4646', 'C:4771', 'C:4580', 'C:4612', 'C:4645', 'C:4560', 'C:4589', 'C:4591', 'C:4618', 'C:4550', 'C:4658', 'C:5531', 'C:4597', 'C:4754', 'C:4598', 'C:4777', 'C:4643', 'C:4648', 'C:4684', 'C:4669', 'C:4553', 'C:4644', 'C:4705', 'C:5812', 'C:4586', 'C:4628', 'C:4545', 'C:4622', 'C:4620', 'C:4766', 'C:4615', 'C:5828', 'C:4675', 'C:4582', 'C:4540', 'C:4555', 'C:4769', 'C:4590', 'C:4548', 'C:4559', 'C:4680', 'C:4763', 'C:4607', 'C:4566', 'C:5809', 'C:4558', 'C:4541', 'C:4783', 'C:5367A', 'C:4699', 'C:5006B', 'C:4547', 'C:4786', 'C:4728', 'C:4613', 'C:4994B', 'C:4621', 'C:4631', 'C:4652', 'C:4632', 'C:4581', 'C:4606', 'C:4647', 'C:4686', 'C:4659', 'C:4608', 'C:4578', 'C:4638', 'C:4570', 'C:4551', 'C:4546', 'C:4554', 'C:4713', 'C:4573', 'C:4593', 'C:4602', 'C:4571', 'C:4610', 'C:5007B', 'C:4641', 'C:4653', 'C:4587', 'C:4747', 'C:4617', 'C:4770', 'C:4735', 'C:4656', 'C:4767', 'C:4609', 'C:5853B', 'C:5510', 'C:4744', 'C:4670', 'C:4730', 'C:5492', 'C:4739', 'C:4596', 'C:4642', 'C:4603', 'C:4599', 'C:4654', 'C:4695', 'C:5996', 'C:4635', 'C:4693', 'C:4677', 'C:4601', 'C:4776', 'C:5391', 'C:4539', 'C:4738', 'C:4634', 'C:4736', 'C:4698', 'C:4772', 'C:5079', 'C:4690', 'C:5015', 'C:4855', 'C:5149', 'C:5856B', 'C:4742', 'C:6122', 'C:4666', 'C:5018B', 'C:4737', 'C:4746', 'C:4720', 'C:4592', 'C:4729', 'C:4572', 'C:4689', 'C:4544', 'C:6007', 'C:5519', 'C:4633', 'C:5516', 'C:4955', 'C:4732', 'C:4624', 'C:4792', 'C:4753', 'C:4594', 'C:4562', 'C:4639', 'C:4708', 'C:4660', 'C:4685', 'C:5508', 'C:4727', 'C:4745', 'C:5324', 'C:4765', 'C:4579', 'C:5118', 'C:5521', 'C:4741', 'C:5205', 'C:5307A', 'C:4676', 'C:4673', 'C:5272', 'C:4600', 'C:4768', 'C:4710', 'C:4796', 'C:5025', 'C:5125', 'C:4709', 'C:4761', 'C:5222A', 'C:5572', 'C:5984', 'C:5122A', 'C:4719', 'C:5000', 'C:5076', 'C:5349', 'C:4968', 'C:5386A', 'C:5120', 'C:4733', 'C:4717', 'C:4726', 'C:4714', 'C:6005', 'C:5225A', 'C:5341A', 'C:5031', 'C:5179', 'C:4723', 'C:5131A', 'C:4712', 'C:4716', 'C:5175A', 'C:5858B', 'C:4711', 'C:5063A', 'C:5528', 'C:6129A', 'C:5274A', 'C:5055A', 'C:5969B', 'C:6074', 'C:4794', 'C:5081', 'C:5509', 'C:5300A', 'C:5212A', 'C:5339A', 'C:4734', 'C:5351', 'C:5512', 'C:4724', 'C:5039', 'C:5126', 'C:5312A', 'C:5833', 'C:6004', 'C:4725', 'C:5522', 'C:4715', 'C:5353', 'C:5308', 'C:5314', 'C:5077', 'C:5504B', 'C:4722', 'C:5968B', 'C:5155', 'C:5185A', 'C:5243', 'C:5530', 'C:6129B', 'C:5113', 'C:5261', 'C:5505', 'C:5524', 'C:6126', 'C:5805', 'C:5218', 'C:5030A', 'C:5038', 'C:4707', 'C:5245', 'C:6125', 'C:5278A', 'C:4718', 'C:5533', 'C:5123', 'C:5294A', 'C:5132', 'C:5356A', 'C:5228A']       
        missing = [ n for n in new if n not in processed]
        
        #missing = ['C:4843'] # FF testing faulty one
        #print ('FF sodblist:', sodblist , 'num: ', len(sodblist))
        #print ('FF missing:', missing , 'num: ', len(missing))
        
        sodblist = ['/raid60/scratch/leo/scratch/era5/odbs/3188/era5.3188.conv.' + n for n in missing ]
        
        '''
        
        sodblist = sodblist[1000:1050]
        print('FF sodblist', sodblist)
        if debug:
            print('The selected odb list is: ',sodblist, 'made of ', len((sodblist)), ' data files. ***')
            print('The ipath variable is: ', ipath+'era5.conv.?:*[!.nc]')
        print('sodblist is: FF', sodblist)
        process_list(sodblist,pool)

    elif single_stat: # select a single station
        exp = single_stat.split('/')[len(single_stat.split('/'))-2] 
        print('The exp is', exp)
        sodblist=[single_stat] 
        process_list(sodblist,pool)

        
        
        
        
global exp    

if __name__ == "__main__":
    """ Call to the run_converter function 
    
        Example:
                 stat = '/raid60/scratch/leo/scratch/era5/odbs/1/era5.conv._98426' # selection of a single station
                 run_converter(dataset= False, single_stat= stat, pool=False, varno=2, debug = True )
                 
                 exp = '1'
                 run_converter(dataset= '1', single_stat= False, pool=20, varno=110, debug = True )  
                 
                 varno: 2(temp.), 110(u,v wind)
                 
                 
        Note that the single station option needs somes fixes (I did not fix where the output is saved)
        Remember to deactiovate the SLIMMER reduction
        Remember to select the correct range of years
        
    """
    
    stat = '/raid60/scratch/leo/scratch/era5/odbs/1/era5.conv._98426' # random station small size for single test
    #run_converter(dataset= ['1'], single_stat= False, pool=False, varno=2, debug = True )  
    
    #datasets = ['1','1761','1759']
    #stat = '/raid60/scratch/leo/scratch/era5/odbs/3188/era5.3188.conv.C:4751'
    stat = False
    #datasets = ['1','1761','1759','3188'] # 3188 gives: builtins.IndexError: i ndex 30533 is out of bounds for axis 0 with size 1476 when using pool
    
    datasets = ['1']
    variables = [7]
    
    for e in datasets:
        exp = e
        for v in variables:
            run_converter(dataset= e, single_stat= False, pool=False, varno= v, debug = True )   
            print('Finished with the database', e , ' **** for the variable: ', str(v))
            
    exit()



print('Done! *** ')    

'''
if __name__ == "__main__":

    gribpath=os.path.expandvars('/raid60/scratch/leo/scratch/ERApreSAT/')

     # choose '1' for the first set of observations
     
#    for exp in ['2824','2928','2504','2502','2506']:
    for exp in []: #'1']:#['3116']:#['e5']:#,'2502','2506']:
#	ipath=os.path.expandvars('$SCRATCH/ectrans/era5/'+exp+'/')
        ipath=os.path.expandvars('/raid60/scratch/leo/scratch/era5/odbs/'+exp+'/')
#        sodblist=glob.glob(ipath+'/*.conv._*')
        sodblist=glob.glob(ipath+'/*.conv._*') # select Vienna to try , or to make it faster use 91364 (below 1 MB)
        sodblist.sort(key=os.path.getsize)
        s=numpy.cumsum(list(os.path.getsize(f) for f in sodblist))

        stride=30
        p = Pool(20)
        idx=[0]
        k=1e9
        x=numpy.where(s>k)[0]
        while len(x)>0:
            k+=1e9
            idx.append(x[0])
            x=numpy.where(s>k)[0]
        idx2=list(range(len(idx)))   
        idx.append(len(s))
        #func = partial(odb2netcdf,None,sodblist,2,par_read_odbsql_stn_withfeedback,idx)
        #p.map(func,  idx2)
        func = partial(odb2netcdf,gribpath,sodblist,3,par_read_odbsql_stn_withfeedback,idx)
        p.map(func,  idx2)
        #list(map(func,  idx2))

    for exp in []:#['erapresat']:#['erapresat']:#['erapresat','erai']:
        ipath=os.path.expandvars('/raid60/scratch/leo/scratch/era5/odbs/'+exp+'/')
        sodblist=glob.glob(ipath+'/*.conv._*')
        sodblist.sort(key=os.path.getsize)
        s=numpy.cumsum(list(os.path.getsize(f) for f in sodblist))

        stride=30
        p=Pool(20)
        idx=[0]
        k=1e9
        x=numpy.where(s>k)[0]
        while len(x)>0:
            k+=1e9
            idx.append(x[0])
            x=numpy.where(s>k)[0]
        idx2=list(range(len(idx)))   
        idx.append(len(s))
        func = partial(odb2netcdf,gribpath,sodblist,2,par_read_odbsql_stn_withfeedback,idx)
        p.map(func,  idx2)
        #func = partial(odb2netcdf,gribpath,sodblist,3,par_read_odbsql_stn_withfeedback,idx)
        #p.map(func,  idx2)

    for exp in []:#['erai']:#['erapresat','erai']:
        ipath=os.path.expandvars('/raid60/scratch/leo/scratch/Date_and_Station_files/Stationfiles_complete/')
        sodblist=glob.glob(ipath+'/??????_t.txt.gz')
        sodblist.sort(key=os.path.getsize)
        s=numpy.cumsum(list(os.path.getsize(f) for f in sodblist))

        stride=30
        p=Pool(20)
        idx=[0]
        k=1e9
        x=numpy.where(s>k)[0]
        while len(x)>0:
            k+=1e9
            idx.append(x[0])
            x=numpy.where(s>k)[0]
        idx2=list(range(len(idx)))   
        idx.append(len(s))
        func = partial(odb2netcdf,gribpath,sodblist,2,par_read_odbascii_stn_withfeedback,idx)
        p.map(func,  idx2)
        #func = partial(odb2netcdf,gribpath,sodblist,3,par_read_odbascii_stn_withfeedback,idx)
        #p.map(func,  idx2)

    for exp in []:#['igra']:#['3188']: #'3188']:#['1759','2491','1761','1770']:
        ipath=os.path.expandvars('/raid60/scratch/leo/scratch/'+exp+'/')
        sodblist=glob.glob(ipath+'*.txt')
#        sodblist=glob.glob(ipath+'*57749*.txt')
        sodblist.sort(key=os.path.getsize)
        s=numpy.cumsum(list(os.path.getsize(f) for f in sodblist))

        stride=30
        p = Pool(25)
        idx=[0]
        k=1e9
        x=numpy.where(s>k)[0]
        while len(x)>0:
            k+=1e9
            idx.append(x[0])
            x=numpy.where(s>k)[0]
        idx2=list(range(len(idx)))   
        idx.append(len(s))
        func = partial(odb2netcdf,gribpath,sodblist,2,par_read_igra_stn_nofeedback,idx)
        p.map(func,  idx2)
        #func = partial(odb2netcdf,gribpath,sodblist,111,par_read_igra_stn_nofeedback,idx)
        #map(func,  idx2)

    #exit()
    for exp in ['3188','1759','1761']:#['3188']:#'1759','1761']: #,'3188']:#,'2491','1770']:
        ipath=os.path.expandvars('/raid60/scratch/leo/scratch/era5/odbs/'+exp+'/')
        sodblist=glob.glob(ipath+'/era5.'+exp+'.conv.?:*[!.nc]')
        #sodblist=glob.glob(ipath+'/era5.'+exp+'.conv.C:4490*')
        sodblist.sort(key=os.path.getsize)
        s=numpy.cumsum(list(os.path.getsize(f) for f in sodblist))

        stride=30
        p = Pool(25)
        idx=[0]
        k=1e9
        x=numpy.where(s>k)[0]
        while len(x)>0:
            k+=1e9
            idx.append(x[0])
            x=numpy.where(s>k)[0]
        idx2=list(range(len(idx)))   
        idx.append(len(s))
        #func = partial(odb2netcdf,gribpath,sodblist,2,par_read_odbsql_stn_nofeedback,idx)
        #p.map(func,  idx2)
        func = partial(odb2netcdf,gribpath,sodblist,111,par_read_odbsql_stn_nofeedback,idx)
        p.map(func,  idx2)
        #list(map(func,  idx2))

    #exit()
    for exp in ['ai_bfr']:#['ai_bfr']: 
        ipath=os.path.expandvars('/raid60/scratch/leo/scratch/ERA40/'+exp+'/')
#        sodblist=glob.glob(ipath+'/era5.'+exp+'.conv.?:*')
        sodblist=glob.glob(ipath+'/era5.*1135*.bfr')
#        sodblist=glob.glob(ipath+'/era5.62414.bfr')
        sodblist.sort(key=os.path.getsize)
        s=numpy.cumsum(list(os.path.getsize(f) for f in sodblist))

        stride=25
        p = Pool(25)
        idx=[0]
        dk=0.1e9
        k=dk
        x=numpy.where(s>k)[0]
        while len(x)>0:
            k+=dk
            idx.append(x[0])
            x=numpy.where(s>k)[0]
        idx2=list(range(len(idx)))   
        idx.append(len(s))
        func = partial(odb2netcdf,gribpath,sodblist,2,par_read_bufr_stn_nofeedback,idx)
        list(map(func,  idx2))
        #func = partial(odb2netcdf,gribpath,sodblist,111,par_read_bufr_stn_nofeedback,idx)
        #p.map(func,  idx2)

    exit()

'''
