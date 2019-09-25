# -*- coding: utf-8 -*-
# this file is released under public domain and you can use without limitations

#########################################################################
## This is a sample controller
## - index is the default action of any application
## - user is required for authentication and authorization
## - download is for downloading files uploaded in the db (does streaming)
## - api is an example of Hypermedia API support and access control
#########################################################################

import matplotlib.pylab as plt
import numpy

def index():
    return dict()

def first():
    return 'Hi there'

def second():
    return request

def third():
    n=eval(request.vars.name)
    x=numpy.arange(n)
    y=x*x
    plt.plot(x,y)
#    dbg.set_trace()
    
    s=sum(x)
    plt.title('{0}'.format(s))
    plt.savefig('applications/myapp/static/images/test.png')
    plt.clf()
    return IMG(_src=URL('static','images/test.png'),id='testfig', _onclick="jQuery(this).slideToggle()")

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



def call():
    """
    exposes services. for example:
    http://..../[app]/default/call/jsonrpc
    decorate with @services.jsonrpc the functions to expose
    supports xml, json, xmlrpc, jsonrpc, amfrpc, rss, csv
    """
    return service()
