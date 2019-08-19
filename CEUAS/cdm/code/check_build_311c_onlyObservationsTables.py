#!/usr/bin/env python
import sys
import os.path
import glob

import subprocess
import urllib.request
import xarray as xr
import numpy
#import h5pickle as h5py
import h5py
from datetime import date, datetime
import time
from multiprocessing import Pool
from netCDF4 import Dataset
import gzip
import pandas as pd
from functools import partial
#from rasotools.utils import *
#from eccodes import *
from numba import *
import matplotlib.pylab as plt
import cartopy.crs as ccrs
import argparse
import copy
from io import StringIO
import h5netcdf

"""  using fixed lenght strings (80 chars) 
Compression does not work well with normal strings 

"""
okinds={'varchar (pk)':numpy.dtype('|S80'),'varchar':numpy.dtype('|S80'),'numeric':numpy.float32,'int':numpy.int32,
       'timestamp with timezone':numpy.datetime64,
       'int[]*':list,'int[]':list,'varchar[]*':list,'varchar[]':list}


kinds={'varchar (pk)':str,'varchar':str,'numeric':numpy.float32,'int':numpy.int32,
       'timestamp with timezone':numpy.datetime64,
       'int[]*':list,'int[]':list,'varchar[]*':list,'varchar[]':list}



def make_datetime(dvar,tvar):
    """ Converts into date-time """
    dvari=dvar.values.astype(numpy.int)
    tvari=tvar.values.astype(numpy.int)
    df=pd.DataFrame({'year':dvar//10000,'month':(dvar%10000)//100,'day':dvar%100,
                        'hour':tvar//10000,'minute':(tvar%10000)//100,'second':tvar%100})
    dt=pd.to_datetime(df).values
    #print('dt id', dt)
    #input('FF check datetime')
    return dt


""" Translates the odb variables name  into cdm var """
cdmfb={'observation_value':'obsvalue@body',
       'observed_variable':'varno@body',
       'z_coordinate_type':'vertco_type@body',
       'z_coordinate':'vertco_reference_1@body',
       'date_time':[make_datetime,'date@hdr','time@hdr'],
       'longitude':'lon@hdr',
       'latitude':'lat@hdr'}
    
def read_all_odbsql_stn_withfeedback(odbfile):

    #countlist=glob.glob(opath+'/*.count')
    alldata=''

    alldict=xr.Dataset()
    t=time.time()
    sonde_type=True
    obstype=True
    if os.path.getsize(odbfile)>0:
        """ Read first the odb header to extract the column names and type """
        try:
            rdata=subprocess.check_output(["odb","header",odbfile])
            rdata=rdata.decode('latin-1').split('\n')
            columns=[]
            kinds=[]
            tdict={}
            for r in rdata[2:-2]:
                try:
                    print(r[:6])
                    if r[:6]=='Header':
                        break
                    else:    
                        columns.append(r.split('name: ')[1].split(',')[0])
                        kinds.append(r.split('type: ')[1].split(',')[0])
                        if kinds[-1]=='REAL':
                            tdict[columns[-1]]=numpy.float32
                        elif kinds[-1] in ('INTEGER','BITFIELD'):
                            if columns[-1]=='date@hdr':
                                tdict[columns[-1]]=numpy.int32
                            else: 
                                tdict[columns[-1]]=numpy.float32
                        else:
                            tdict[columns[-1]]=numpy.dtype('S8') # dict containng column name and type
                                     
                            
                except IndexError:
                    pass
        except:
            print('could not read odbfile '+odbfile)
            return alldict
        try:
            rdata=subprocess.check_output(["odb","sql","-q","select *","-i",odbfile,'--no_alignment']) # after reading the header it does the query
            # returns a byte string
            print('after odb:',time.time()-t)
            rdata=''.join(rdata.decode('latin-1').split("'")) # decoding the stirng into a unicode
            f=StringIO(rdata) # access the string like a file, return a file pointer to read the string with pandas
            # nb  if you have null values, reading of integer fails and are read as floats
            # to improve, you can convert the columns with nabs in the alldicts (pandas data frame) into int(np.nans)
            # date values are large so float 32 precision is not sufficient 
            alldict=pd.read_csv(f, delimiter='\t', quoting=3, comment='#',dtype=tdict) # example test: (21150x63 columns)
            del f,rdata

 
            """ alternative method to read the odb           
            if False:
                
                rdata='nan'.join(rdata.decode('latin-1').split('NULL'))
                rdata=''.join(rdata.split("'"))
                rdata='\t'.join(rdata.split('\n')[1:-1])
                rdata=tuple(rdata.split('\t'))
                
                print('after odb:',time.time()-t)
                #xdata=numpy.fromstring(rdata,sep='\t')
                cl=len(columns)
                rl=len(rdata)//cl
                #for k in range(cl):
                    #if kinds[k]=='REAL':
                        #alldict[columns[k]]=({'obslen':rl},numpy.empty(rl,dtype=numpy.float32))
                    #elif kinds[k] in ('INTEGER','BITFIELD'):
                        #alldict[columns[k]]=({'obslen':rl},numpy.empty(rl,dtype=numpy.int32))
                    #else:
                        #alldict[columns[k]]=({'obslen':rl},numpy.empty(rl,dtype='|S8'))
                #print(odbfile,time.time()-t)
                for k in range(cl):
                    if kinds[k]=='REAL':
                        alldict[columns[k]]=({'obslen':rl},numpy.float32(rdata[k::cl]))
                    elif kinds[k] in ('INTEGER','BITFIELD'):
                        rds=rdata[k::cl]
                        if 'nan' in rds: 
                            alldict[columns[k]]=({'obslen':rl},numpy.asarray(rds,dtype=float).astype(numpy.int32))
                        else:
                            alldict[columns[k]]=({'obslen':rl},numpy.asarray(rds,dtype=numpy.int32))
                    else:
                        alldict[columns[k]]=({'obslen':rl},numpy.asarray(rdata[k::cl],dtype='|S8'))
             """
       
            #input('FF check alldict')
            # alldict.head() # pandas method to print the columns name . alldict is now a panda dataframe with e.g. 22000*63 columns where the 63 columns are read from the odb file 
            
        except subprocess.CalledProcessError as e:
            print('odb failed!:'+' '+odbfile)
            return alldict

    print(odbfile,time.time()-t)
    #idy=numpy.lexsort((alldict['varno@body'],
                       #-alldict['vertco_reference_1@body'],
                       #alldict['time@hdr'],
                       #alldict['date@hdr']))
    #for k in alldict.columns:
        #alldict[k]=alldict[k][idy]

    print(odbfile,time.time()-t)

    """ may not be necessary to convert into x_array sicne you can write a pandas df into an HDF file """

    #print('FF alldict_to_xarray is: ', alldict.to_xarray )
   # input('verifica alldicts.to_xarray')

    return alldict.to_xarray()


def fromfb(fbv,cdmfb,cdmkind):
    """ 
         Args: 
               fbv    : xarray converted form the odb file, still with odb variable names
               cdmfb  : dictionary converting the odb variable names to CDM names 
               cdmkind: data type of the cdmfb
         Returns:
               values of the variables 

    """
    x=0
    # checks if the type of the variable is a list, so that it uses the function to extract the date time 

    #print('fbv, cdmfb, cdmkind are', fbv,cdmfb,cdmkind)

    ## 'date_time':[make_datetime,'date@hdr','time@hdr'],
    if type(cdmfb) is list:
        x=cdmfb[0](fbv[cdmfb[1]],fbv[cdmfb[2]]) # equals to make_datetime('date@hdr','time@hdr') where 'date@hdr','time@hdr' are read from the fbv xarray
        #print('x is', x )
        #input('FF check x time')
    else:
        if cdmfb=='varno@body':
            tr=numpy.zeros(113,dtype=int) 
            """ translates odb variables number to Lot3 numbering convention """
            tr[1]=117  # should change
            tr[2]=85
            tr[3]=104
            tr[4]=105
            tr[7]=39 #spec hum
            tr[29]=38 #relative hum
            tr[59]=36 # dew point
            tr[111]=106 #dd
            tr[112]=107  #ff
            
            tr[39]= 85 # 2m T
            tr[40]= 36 # 2m Td
            tr[41]= 104 #10m U
            tr[42]= 105  #10m V
            tr[58]= 38 # 2m rel hum
            
            x=tr[fbv[cdmfb].values.astype(int)] # reads the varno from the odb feedback and writes it into the variable id of the cdm 
        else:    
            x=fbv[cdmfb].values
    #print('check x', x )
    #input('check new X')    
    return x



def ttrans(cdmtype,kinds=kinds):
    """ convert the cdm types to numpy types """    
    nptype=numpy.float32
    try:
        nptype=kinds[cdmtype.strip()]
    except:
        print(cdmtype,'not found, using numpy.float32')
            
    return nptype

@njit
def find_dateindex(y):
    """ creates the indices list from the dates, for quick access 
        nb the benchmark script will not work with these files since the definition of the array size is swapped i.e. (x.shape[0], 3)"""        


    x=numpy.unique(y)
    z=numpy.zeros((3,x.shape[0]),dtype=numpy.int32)
    z-=1
    j=0
    for i in range(len(y)):
        m=y[i]
        if x[j]==y[i]:
            if z[1,j]==-1:
                z[1,j]=i
                #print(j,i)
            else:
                if z[2,j]<i:
                    z[2,j]=i
        elif x[j]<y[i]:
            j+=1
            if x[j]==y[i]:
                if z[1,j]==-1:
                    z[1,j]=i
                    #print(j,i)
                else:
                    if z[2,j]<i:
                        z[2,j]=i
            else:
                print('Error')
        else:
            j-=1
            if x[j]==y[i]:
                if z[1,j]==-1:
                    z[1,j]=i
                    #print(j,i)
                else:
                    if z[2,j]<i:
                        z[2,j]=i
            else:
                print('Error')
    z[0,:]=x
    return z

def odb_to_cdm(cdm, cdmd, fn):
    """ Convert the odb file into cdm compliant netCDF files
        input:
              fn   :: odb file name (e.g. era5.conv._10393)
              cdm  :: cdm tables (read with pandas)
              cdmd :: cdm tables definitions ("") """

    recl=0    
    t=time.time()

    fnl=fn.split('/')
    fnl[-1]='ch'+fnl[-1]
    fno=output_dir + '/' + fnl[-1] + '.nc' # creating an output file name e.g. chera5.conv._10393.nc  , try 01009 faster
    if not False:
        
        fbds=read_all_odbsql_stn_withfeedback(fn) # this is the xarray converted from the pandas dataframe         
        print(time.time()-t) # to check the reading of the odb

        fbencodings={}
        for d in fbds._variables.keys(): # variables stored in the xarray 
            if fbds.variables[d].dtype==numpy.dtype('float64'):
                if d!='date@hdr':             
                    fbencodings[d]={'dtype':numpy.dtype('float32'),'compression': 'gzip'} # probably dtype not neccessary, but compression must be there
                else:
                    fbencodings[d]={'dtype':numpy.dtype('int32'),'compression': 'gzip'}               
            else:
                fbencodings[d]={'compression': 'gzip'}
                
        y=fbds['date@hdr'].values
        z=find_dateindex(y)
        di=xr.Dataset() 
        di['dateindex']=({'days':z.shape[1],'drange':z.shape[0]},z) # date, index of the first occurrance, index of the last 
    
        """ Empty dictionaries for hdf groups (one for each table) and their hdf encoding description """
        groups={}
        groupencodings={}
        
        """ Loop over every table_definition in the cdmd dictionary of pandas DF """
        for k in cdmd.keys(): # loop over all the table definitions e.g. ['id_scheme', 'crs', 'station_type', 'observed_variable', 'station_configuration', 'station_configuration_codes', 'observations_table', 'header_table'] 
            groups[k]=xr.Dataset() # create an  empy xarray
            groupencodings[k]={} # create a dict of group econding

            """ Loop over each row in the table """
            for i in range(len(cdmd[k])):  
                d=cdmd[k].iloc[i] # s
                if k in ('observations_table'): 
                    try:
                        groups[k][d.element_name]=({'hdrlen':fbds.variables['date@hdr'].shape[0]},
                                    fromfb(fbds._variables, cdmfb[d.element_name], ttrans(d.kind,kinds=okinds)) )
                        
                    except KeyError:
                        x=numpy.zeros(fbds.variables['date@hdr'].values.shape[0],dtype=numpy.dtype(ttrans(d.kind,kinds=okinds)))
                        x.fill(numpy.nan)
                        groups[k][d.element_name]=({'hdrlen':fbds.variables['date@hdr'].shape[0]},x)
                             
                try:
                    groups[k][d.element_name].attrs['external_table']=d.external_table # defining variable attributes that point to other tables (3rd and 4th columns)
                    groups[k][d.element_name].attrs['description']=d.description
                    groupencodings[k][d.element_name]={'compression': 'gzip'}
                except KeyError:
                    print('bad:',k,d.element_name)
                    pass

        #this writes the dateindex to the netcdf file. For faster access it is written into the root group
        di.to_netcdf(fno, format='netCDF4',engine='h5netcdf',mode='w')
        # add era5 feedback to the netcdf file. 
        # fbds is the 60 col. xarray
        fbds.to_netcdf(fno, format='netCDF4',engine='h5netcdf',encoding=fbencodings, group='era5fb', mode='a') # f
        
        for k in groups.keys():
            groups[k].to_netcdf(fno,format='netCDF4',engine='h5netcdf',encoding=groupencodings[k],group=k,mode='a') #
            
        print('sizes: in: {:6.2f} out: {:6.2f}'.format(os.path.getsize(fn)/1024/1024,
                                              os.path.getsize(fno)/1024/1024))
        del fbds
    

    # speed test for accessing dateindex
    #for k in range(10):
        
        #t=time.time()
        #with h5py.File(fno,'r') as f:
            #di=f['dateindex'][:]
            #idx=numpy.where(di[0,:]==19950101)[0]
            #print(numpy.nanmean(f['observations_table']['observation_value'][di[1,idx[0]]:di[2,idx[0]]+1]))
        #print('ch',time.time()-t)
        #fng='cgera5'.join(fno.split('chera5'))
        #t=time.time()
        #with h5py.File(fng,'r') as g:
            #do=g['dateindex'][:]
            #idx=numpy.where(do[:,2]==19950101)[0]
            #print(numpy.nanmean(g['obsvalue@body'][do[idx[0],0]:do[idx[0],1]+1]))
        #print('cg',time.time()-t)
        
        
    # first implementation of inner join  - resolving numeric variable code 
    #with h5py.File(fno,'r') as f:
        #ext=f['observations_table']['observed_variable'].attrs['external_table']
        #lext=ext.split(':')
        
        #l=[]
        #lidx=[]
        #llen=len(f['observations_table']['observed_variable'][:])
        #for i in range(llen):
            #obv=f['observations_table']['observed_variable'][i]
            #if obv not in l:
                #idx=numpy.where(f[lext[0]][lext[1]][:]==obv)[0][0]
                #l.append(obv)
                #lidx.append(idx)
            #else:
                #idx=lidx[l.index(obv)]
            #for k in f[lext[0]].keys():
                #print(k,':',f[lext[0]][k][idx],)
            #print('value:',f['observations_table']['observation_value'][i])
        #print(f)

    
    print(fno,time.time()-t)

    return recl


def csvListFromUrls(url=''):
    """ Return a list of csv files, as fond in the url on the cdm GitHub """   
    urlpath = urlopen(url)
    string = urlpath.read().decode('utf-8')
    split = string.split(' ')
    csv_files_list = [m.replace('"','') for m in [n.split('title="')[1] for n in split if '.csv' in n and "title" in n] ] 
    return csv_files_list

if __name__ == '__main__':


    parser = argparse.ArgumentParser(description="Make CDM compliant netCDFs")
    parser.add_argument('--database_dir' , '-d', 
                    help="Optional: path to the database directory. If not given, will use the files in the data directory" ,
                    default = '../../cdm/data/',
                    type = str)
    parser.add_argument('--auxtables_dir' , '-a', 
                    help="Optional: path to the auxiliary tables directory. If not given, will use the files in the data/tables directory" ,
                    default = '../../cdm/data/',
                    type = str)
    parser.add_argument('--odbtables_dir' , '-odb', 
                    help="Optional: path to the odb tables directory. If not given, will use the files in the data/tables directory" ,
                    #default = '../../cdm/data/tables/',
                    default = '/raid60/scratch/leo/scratch/era5/odbs/',
                    type = str)

    parser.add_argument('--output_dir' , '-out',
                    help="Optional: path to the netcdf output directory" ,
                    default = '../../cdm/data/tables',
                    type = str)


    args = parser.parse_args()
    dpath = args.database_dir
    tpath = args.auxtables_dir
    output_dir = args.output_dir
    odbpath = args.odbtables_dir

    """
    print ('THE DPATH IS', dpath)
    if not dpath:
        dpath = '../cdm/code/data/'
   
    if not tpath:
        tpath = dpath+'/tables/'
   
    """
    print ('Analysing the databases: ')
    print ('dpath is ' , dpath)
    print ('tpath is ' , tpath )


    cdmpath='https://raw.githubusercontent.com/glamod/common_data_model/master/tables/' # cdm tables            

    
    """ Selecting the list of table definitions. Some of the entires do not have the corresponding implemented tables """
    cdmtabledeflist=['observations_table']
    cdm_tabdef = dict()
    for key in cdmtabledeflist:
        url='table_definitions'.join(cdmpath.split('tables'))+key+'.csv' # https://github.com/glamod/common_data_model/tree/master/table_definitions/ + ..._.dat 
        f=urllib.request.urlopen(url)
        col_names=pd.read_csv(f,delimiter='\t',quoting=3,nrows=0,comment='#')
        f=urllib.request.urlopen(url)
        tdict={col: str for col in col_names}
        cdm_tabdef[key]=pd.read_csv(f,delimiter='\t',quoting=3,dtype=tdict,na_filter=False,comment='#')
        
    
    """ Selecting the list of tables. 'station_configuration_codes','observations_table','header_table' are not implemented in the CDM GitHub"""        
    cdm_tab=dict() # dictionary where each key is the name of the cdm table, and the value is read from the .dat file    
    
    cdm_tab['observations_table']=pd.read_csv(tpath+'/table_definitions/observations_table.csv',delimiter='\t',quoting=3,comment='#')





    """ Up to here, we have two different dictionaries.
    cdm_tab: contains the dictionary for the cdm tables
    cdm_tabdef: contains the dictionary for the cdm tables definitions
    """


    #func=partial(odb_to_cdm, cdm=cdm_tab, cdmd=cdm_tabdef)
    func=partial(odb_to_cdm, cdm_tab, cdm_tabdef)

    #bfunc=partial(read_bufr_stn_meta,2)
    #rfunc=partial(read_rda_meta) 
    tu=dict()
    p=Pool(25)
    dbs=['1'] #,'igra2','ai_bfr','rda','3188','1759','1761']
    for odir in dbs: 
        cdm_tab['station_configuration']=pd.read_csv(odbpath+'/'+odir+'/station_configuration.dat',delimiter='\t',quoting=3,dtype=tdict,na_filter=False,comment='#')
        if 'ai' in odir:
            pass
        elif 'rda' in odir:
            pass
        elif 'igra2' in odir:
            pass

        else:
            #flist=glob.glob(odbpath+odir+'/'+'era5.conv.*01009.nc.gz')
            flist=glob.glob(odbpath+odir+'/'+'era5.conv._01009')
            transunified=list(map(func,flist))

