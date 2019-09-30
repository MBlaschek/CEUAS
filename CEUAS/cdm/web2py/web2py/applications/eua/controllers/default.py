# -*- coding: utf-8 -*-
# this file is released under public domain and you can use without limitations

#########################################################################
## This is a simple controller
## - index is the default action of any application
##   here it forwards the request from a CDS server to further processing 
##   with python
##
##   For queries ask early-upper-air@copernicus-climate.eu
##   The C3S 311c Lot2 group
##   Vienna, 26 August 2019
## - user is required for authentication and authorization
## - download is for downloading files uploaded in the db (does streaming)
## - api is an example of Hypermedia API support and access control
#########################################################################

import os,sys
sys.path.append(os.path.expanduser('~/python/'))
sys.path.append('/data/private/soft/python/')
from gluon.debug import dbg
import cds_eua as eua
import pandas as pd
import xarray
import numpy

def index():
    """
    example action using the internationalization operator T and flash
    rendered by views/default/index.html or views/generic.html

    if you need a simple wiki simply replace the two lines below with:
    return auth.wiki()
    """
    #print(dict(request).keys())
    #return 'Wort'
    rvars=dict(request.vars)
    print(rvars)
    if not rvars:
        return  HTML(BODY(H1(T('This is the backend for the Copernicus CDS Early Upper air service.')),_style="color: red;"), 
                     H1(T('It is in alpha stage.')), H1(T(' Once approved, it should be accessed via')),
                     H1(A('https://cds.climate.copernicus.eu',_href='https://cds.climate.copernicus.eu'))).xml() # .xml to serialize
        
        #return dict(message=XML(cdmd.to_html(table_id='First',render_links=True)))
        #return cdmd.to_html(table_id='First',render_links=True)
    print(rvars.keys())
    #if 'format' in rvars.keys():
        #if rvars['format']=='html':
            #cdmd=pd.read_csv('applications/eua/static/observations_table.csv',delimiter='\t',quoting=3,na_filter=False,comment='#')
            #return dict(message=XML(cdmd.to_html(table_id='First',render_links=True)))
            
    rfile=''
    if 'group' not in rvars:
        rfile,error=eua.process(rvars)
    else:
        cdmd=eua.cdmexplainer(rvars)
        return dict(content=XML(cdmd.to_html(table_id='First',escape=False)))#,index=False
        
    #rfile='applications/testajax/static/images/chera5.conv._'+rvars['name']+'.nc'
    print(rfile)
    if(os.path.isfile(rfile)):
        if 'format' in rvars.keys():
            if rvars['format']=='html':
                vars=['date_time','observation_id','latitude','longitude','z_coordinate','observed_variable','observation_value','units']
                if 'columns' in rvars.keys():  
                    print (rvars['columns'])
                    rvars['columns']=eval(rvars['columns'])
                    print (rvars['columns'])
                    
                    vars=vars+rvars['columns']
                with xarray.open_dataset(rfile,group='observations_table') as f:
                    print(f.keys())
                    x=xarray.Dataset()
                    for v in vars:
                        print(v)
                        if v in ['date_time']:
                            x[v]=f[v]
                            #print(x[v].values)
                            #x[v].values=pd.to_datetime(f[v].values,origin='1979-01-01',unit='h').values
                            #print(v,x[v])
                        else:
                            x[v]=f[v]
                            print(v,'type',type(f[v].values[0]))
                            if type(f[v].values[0]) in [ bytes,numpy.bytes_ ]:
                                x[v].values=x[v].values.astype(str)
                print('vor to_dataframe')
                cdmd=x.to_dataframe()
                print('nach to_dataframe')
                cdmd['observed_variable'] = cdmd['observed_variable'].apply(lambda y: '<a href="?group=observed_variable&variable={0}">{0}</a>'.format(y))
                if 'units' in cdmd.columns:
                    cdmd['units'] = cdmd['units'].apply(lambda y: '<a href="?group=units&units={0}">{0}</a>'.format(y))
                if 'report_id' in cdmd.columns:
                    cdmd['report_id'] = cdmd['report_id'].apply(lambda y: '<a href="?group=header_table&report_id={0}">{0}</a>'.format(y))
                if 'observation_id' in cdmd.columns:
                    cdmd['observation_id'] = cdmd['observation_id'].apply(lambda y: '<a href="?group=era5fb&index={0}">{0}</a>'.format(y))
                    
                print('vor return')
                return dict(message=XML(cdmd.to_html(table_id='First',escape=False,index=False)))#
        else:
            return response.stream(rfile,65536,attachment=False)
    else:
        
        return dict(aaa='Could not create netcdf file',aareason=error,aarequest='',**rvars)
#response.stream('applications/testajax/static/images/chera5.conv._01009.nc',attachment=False)

#def manifest():
    #return response.stream('applications/testajax/static/manifest.txt')

#def manifest2():
    #return response.stream('applications/testajax/static/manifest2.txt')
def user():
    """
    exposes:
    http://..../[app]/default/user/login
    http://..../[app]/default/user/logout
    http://..../[app]/default/user/register
    http://..../[app]/default/user/profile
    http://..../[app]/default/user/retrieve_password
    http://..../[app]/default/user/change_password
    http://..../[app]/default/user/manage_users (requires membership in
    use @auth.requires_login()
        @auth.requires_membership('group name')
        @auth.requires_permission('read','table name',record_id)
    to decorate functions that need access control
    """
    return dict(form=auth())


@cache.action()
def download():
    """
    allows downloading of uploaded files
    http://..../[app]/default/download/[filename]
    """
    return response.download(request, db)


def call():
    """
    exposes services. for example:
    http://..../[app]/default/call/jsonrpc
    decorate with @services.jsonrpc the functions to expose
    supports xml, json, xmlrpc, jsonrpc, amfrpc, rss, csv
    """
    return service()


@auth.requires_login() 
def api():
    """
    this is example of API with access control
    WEB2PY provides Hypermedia API (Collection+JSON) Experimental
    """
    from gluon.contrib.hypermedia import Collection
    rules = {
        '<tablename>': {'GET':{},'POST':{},'PUT':{},'DELETE':{}},
        }
    return Collection(db).process(request,response,rules)
