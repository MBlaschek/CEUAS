#!/opt/anaconda/bin/python


import sys
if sys.version_info[0]!=4:
    import matplotlib.pylab as plt
import os
import numpy
import matplotlib.patches as mpatches
import rasotools.additions.anomaly as rst
from scipy import linalg
from Magics.macro import *

#
from netCDF4 import Dataset
import rasotools.additions.allrasotrends as rsta
import rasotools.additions.allrasodiffs as rstd
import rasotools.additions.statanalysis as rsts
import rasotools.additions.stationts as rst
import rasotools.additions.utils as rstu
import rasotools.additions.dumpandspagarr as rsdump
import rasotools.additions.obstat as obstat
from rasotools.additions.set_trendmaps import *
#import rasotools.anomaly as rstan
#
#f=Dataset('rich_gridded.nc','r')
#f=Dataset('001001/feedbackglobbincorrmon001001.nc','r')
import json
import time
import string
#from   namelist import *
from datetime import date
import glob
import f90nml
#import cherrypy
try:
    from numba import njit,config,threading_layer
    #config.THREADING_LAYER = 'tbb'
except:
    from numba import njit,config
    


from rasotools.additions.anomaly import *
from rasotools.additions.define_datasets import *
import copy


def compare_profiles(plist,varlist,namelist,clist,plotproperties,cvar,beltname='',had=numpy.nan,hadens=numpy.zeros(100)):

    linesplustext=[]
    legend_user_lines=[]

    if plotproperties['mode'][0]!='monthlydiff':
        title='Temperature Trends [K/10a]'
        mode='Trend'
        interval=plotproperties['intervals'][0,:]
    else:
        title='Temperature Difference [K]'	
        mode='Diff'
        interval=plotproperties['diffintervals'][0,:]


    lines =[title+", "+beltname+', '+cvar+','
            "{:4}".format(interval[0])+"-{:4}".format(interval[1])]
    ps=plotproperties['ps']
    pindex=plotproperties['pindex']# numpy.concatenate((plotproperties['pindex'],plotproperties['msupindex']))

    profsat=''
    if 'satseries' in plotproperties['plotlist']:
        profsat='sat'

    ppage = page(
        layout='positional',  
        page_x_length=29.7, 
        page_y_length=21., 
        page_id_line='off',
        page_x_position=0., 
        page_y_position=0.)

    out = output({"output_name_first_page_number":"off",
                  "output_formats":plotproperties["outputformat"], 
                  'output_name': plotproperties['tmppath']+'/'+profsat+mode.lower()+'belt_'+beltname+'_'+cvar+'_'+
                  "{:4}".format(interval[0])+"-{:4}".format(interval[1])})

    amplevels=numpy.where(numpy.logical_and(ps[pindex]>100,ps[pindex]<=300))[0]

    shade=1
    for i in range(len(namelist)):
        if 'andep' in namelist[i]:
            namelist[i]=namelist[i].split('_andep')[0]
        if 'fgdep' in namelist[i]:
            namelist[i]=namelist[i].split('_fgdep')[0]

    amp=numpy.zeros(len(varlist))#tasks[key]["ens"].shape[0])
    for i in range(len(varlist)):
        if 'satseries' in plotproperties['plotlist']:
            for iens in [0]: #tasks[key]["ens"]: #range(len(tasks[key]["ens"])):
                amp[iens]=varlist[i][-3]/had
                linesplustext=linesplustext+addsattrend(varlist[i][-4:],[800.,550.,250.,90.],
                                                        scol=clist[i],iens=iens) # SAT
        else:
            for iens in [0]: #tasks[key]["ens"]: #range(len(tasks[key]["ens"])):
                try:
                    if namelist[i].lower() in ['uah','rss','star']:
                        linesplustext=linesplustext+addsattrend(varlist[i][-4:],[800.,550.,250.,90.],
                                                                scol=clist[i],iens=iens) 
                        ampl=[1]
                    else:
                        linesplustext=linesplustext+addprofile(varlist[i], ps[pindex], 0,0,col=clist[i],pkey=beltname.lower())
                        ampl=amplevels.copy()
                    amp[iens]=numpy.max(varlist[i][ampl])/had
                except:
                    amp[iens]=numpy.nan
        amps=numpy.nanmean(amp[iens])
        if amps==amps:
            amps=" {0:3.1f}".format(numpy.nanmean(amp[iens]))
        else:
            amps=''
        legend_user_lines=legend_user_lines+[namelist[i]+amps]    

    pkey=beltname.lower()
    if pkey=='cost':
        bmax=0.
        for i in range(len(varlist)):
            bmax=numpy.nanmax((bmax,numpy.nanmax(varlist[i])))
        beltinterval=[0.,bmax]
    else:
        beltinterval=[-1.8,0.7]

    projection,horizontal,vertical,title,legend=set_belttrends(lines,legend_user_lines,plotproperties,beltinterval=beltinterval,pkey=beltname.lower())
    if had==had:
        hadtrend=addhadtrend(had,hadens)
        legend_user_lines=legend_user_lines+["HadCRUT4"]

    if len(linesplustext)>0:
        out.args['output_title']=os.getcwd()+'/'+out.args['output_name']	
        if had==had:
            plot(out,ppage,projection,horizontal,vertical,linesplustext,hadtrend,title,legend)
        else:
            plot(out,ppage,projection,horizontal,vertical,linesplustext,title,legend)

        print((out.args['output_name']+'.'+plotproperties["outputformat"][0]))
    else:
       print(( out.args['output_name']+'.'+plotproperties["outputformat"][0]+ 'not created - no data found'))
    return

def compare_tseries(tlist,varlist,namelist,clist,plotproperties,cvar,beltname='',had=numpy.nan):

    linesplustext=[]
    legend_user_lines=[]

    if plotproperties['mode'][0]!='monthlydiff':
        title='Temperature Trends [K/10a]'    
        mode='Trend'
        interval=plotproperties['intervals'][0,:]
    else:
        title='Temperature Difference [K]'	
        mode='Diff'
        interval=plotproperties['diffintervals'][0,:]


    lines =[title+", "+beltname+', '+cvar+','
            "{:4}".format(interval[0])+"-{:4}".format(interval[1])]
    ps=plotproperties['ps']
    pindex=plotproperties['pindex']# numpy.concatenate((plotproperties['pindex'],plotproperties['msupindex']))

    profsat=''
    if 'satseries' in plotproperties['plotlist']:
        profsat='sat'

    ppage = page(
        layout='positional',  
        page_x_length=29.7, 
        page_y_length=21., 
        page_id_line='off',
        page_x_position=0., 
        page_y_position=0.)

    out = output({"output_name_first_page_number":"off",
                  "output_formats":plotproperties["outputformat"], 
                  'output_name': plotproperties['tmppath']+'/'+profsat+mode.lower()+'beltanomalies_'+beltname+'_'+cvar+'_'+
                  "{:4}".format(interval[0])+"-{:4}".format(interval[1])})


    if len(linesplustext)>0:
        out.args['output_title']=os.getcwd()+'/'+out.args['output_name']	
        if had==had:
            plot(out,ppage,projection,horizontal,vertical,linesplustext,hadtrend,title,legend)
        else:
            plot(out,ppage,projection,horizontal,vertical,linesplustext,title,legend)

        print((out.args['output_name']+'.'+plotproperties["outputformat"][0]))
    else:
        print((out.args['output_name']+'.'+plotproperties["outputformat"][0]+ 'not created - no data found'))
    return

def compare_zmeans(latlist,varlist,namelist):
    return

def defcolours():
    colours_=dict()
    colours_["automatic"] = rgb(0., 0., 0.) #// NEEDS FIXING 
    colours_["none"] = rgb(-1., -1., -1.) 
    colours_["background"] = rgb(1., 1., 1.) 
    colours_["foreground"] = rgb(0., 0., 0.) 
    colours_["red"] = rgb(1.0000, 0.0000, 0.0000) 
    colours_["green"] = rgb(0.0000, 1.0000, 0.0000) 
    colours_["blue"] = rgb(0.0000, 0.0000, 1.0000)
    colours_["yellow"] = rgb(1.0000, 1.0000, 0.0000) 
    colours_["cyan"] = rgb(0.0000, 1.0000, 1.0000) 
    colours_["magenta"] = rgb(1.0000, 0.0000, 1.0000) 
    colours_["black"] = rgb(0.0000, 0.0000, 0.0000) 
    colours_["avocado"] = rgb(0.4225, 0.6500, 0.1950) 
    colours_["beige"] = rgb(0.8500, 0.7178, 0.4675) 
    colours_["brick"] = rgb(0.6000, 0.0844, 0.0300) 
    colours_["brown"] = rgb(0.4078, 0.0643, 0.0000) 
    colours_["burgundy"] = rgb(0.5000, 0.0000, 0.1727) 
    colours_["charcoal"] = rgb(0.2000, 0.2000, 0.2000) 
    colours_["chestnut"] = rgb(0.3200, 0.0112, 0.0000) 
    colours_["coral"] = rgb(0.9000, 0.2895, 0.2250) 
    colours_["cream"] = rgb(1.0000, 0.8860, 0.6700) 
    colours_["evergreen"] = rgb(0.0000, 0.4500, 0.2945) 
    colours_["gold"] = rgb(0.7500, 0.5751, 0.0750) 
    colours_["grey"] = rgb(0.7000, 0.7000, 0.7000) 
    colours_["khaki"] = rgb(0.5800, 0.4798, 0.2900) 
    colours_["kelly_green"] = rgb(0.0000, 0.5500, 0.1900) 
    colours_["lavender"] = rgb(0.6170, 0.4070, 0.9400) 
    colours_["mustard"] = rgb(0.6000, 0.3927, 0.0000) 
    colours_["navy"] = rgb(0.0000, 0.0000, 0.4000) 
    colours_["ochre"] = rgb(0.6800, 0.4501, 0.0680) 
    colours_["olive"] = rgb(0.3012, 0.3765, 0.0000) 
    colours_["peach"] = rgb(0.9400, 0.4739, 0.3788) 
    colours_["pink"] = rgb(0.9000, 0.3600, 0.4116) 
    colours_["rose"] = rgb(0.8000, 0.2400, 0.4335) 
    colours_["rust"] = rgb(0.7000, 0.2010, 0.0000) 
    colours_["sky"] = rgb(0.4500, 0.6400, 1.0000) 
    colours_["tan"] = rgb(0.4000, 0.3309, 0.2000) 
    colours_["tangerine"] = rgb(0.8784, 0.4226, 0.0000) 
    colours_["turquoise"] = rgb(0.1111, 0.7216, 0.6503) 
    colours_["violet"] = rgb(0.4823, 0.0700, 0.7000) 
    colours_["reddish_purple"] = rgb(1.0000, 0.0000, 0.8536) 
    colours_["purple_red"] = rgb(1.0000, 0.0000, 0.5000) 
    colours_["purplish_red"] = rgb(1.0000, 0.0000, 0.2730) 
    colours_["orangish_red"] = rgb(1.0000, 0.0381, 0.0000) 
    colours_["red_orange"] = rgb(1.0000, 0.1464, 0.0000) 
    colours_["reddish_orange"] = rgb(1.0000, 0.3087, 0.0000) 
    colours_["orange"] = rgb(1.0000, 0.5000, 0.0000) 
    colours_["yellowish_orange"] = rgb(1.0000, 0.6913, 0.0000) 
    colours_["orange_yellow"] = rgb(1.0000, 0.8536, 0.0000) 
    colours_["orangish_yellow"] = rgb(1.0000, 0.9619, 0.0000) 
    colours_["greenish_yellow"] = rgb(0.8536, 1.0000, 0.0000) 
    colours_["yellow_green"] = rgb(0.5000, 1.0000, 0.0000) 
    colours_["yellowish_green"] = rgb(0.1464, 1.0000, 0.0000) 
    colours_["bluish_green"] = rgb(0.0000, 1.0000, 0.5000) 
    colours_["blue_green"] = rgb(0.0000, 1.0000, 1.0000) 
    colours_["greenish_blue"] = rgb(0.0000, 0.5000, 1.0000) 
    colours_["purplish_blue"] = rgb(0.1464, 0.0000, 1.0000) 
    colours_["blue_purple"] = rgb(0.5000, 0.0000, 1.0000) 
    colours_["bluish_purple"] = rgb(0.8536, 0.0000, 1.0000) 
    colours_["purple"] = rgb(1.0000, 0.0000, 1.0000) 
    colours_["white"] = rgb(1.0000, 1.0000, 1.0000) 
    colours_["undefined"] = rgb(-1., -1., -1.)
    return colours_

@njit
def afill(	a,pindex):
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            for k in range(a.shape[2]):
                for l in range(pindex.shape[0]):
                    for m in range(a.shape[4]):
                        a[i,j,k,pindex[l],m]=numpy.nan


def break_analysis(exp,rfpar,plotproperties,stnames):

    files=glob.glob(plotproperties['tmppath']+'/??????/found_breaks*')
    totalbreaks=0
    abreakts=numpy.zeros(rfpar["mmax"],dtype=int)
    rbreakts=numpy.zeros(rfpar["mmax"],dtype=int)
    refdate=datetime.date(1900,1,1)
    for fn in files:
        with open(fn) as f:
            data=f.read().split('\n')
        for d in data:
            di=d.split()
            if 'accepted' in d:
                bdate=di[1]
                print(d)
    #		    diff=datetime.date(int(bdate[:4]),int(bdate[4:6]),int(bdate[6:8]))-refdate
                if 'sabiner' in exp or 'Control' in exp:
                    abreakts[int(int(bdate)/30.5)]+=1
                else:
                    abreakts[(int(bdate[:4])-1900)*12+int(bdate[4:6])]+=1
            if 'rejected' in d:
                bdate=di[1]
    #		    diff=datetime.date(int(bdate[:4]),int(bdate[4:6]),int(bdate[6:8]))-refdate
                if 'sabiner' in exp or 'Control' in exp:
                    rbreakts[int(int(bdate)/30.5)]+=1
                else:
                    rbreakts[(int(bdate[:4])-1900)*12+int(bdate[4:6])]+=1
#		    rbreakint(bdate[:4]ts[diff.days]+=1

    med=numpy.median(abreakts[abreakts>0])
    x=rmeanw(abreakts*1.0,3)
    print(('med: ',med))
    peaklist=[]
    for i in range(2,abreakts.size-2):
        if abreakts[i]==numpy.max(abreakts[i-2:i+2]) and abreakts[i]>3*med:
            peaklist.append((i,x[i]))


    ab,=plt.plot(int(rfpar['startdate'])//10000+numpy.arange(rfpar["mmax"])/12.,rmeanw(abreakts*1.0,3),label='accepted')
    rb,=plt.plot(int(rfpar['startdate'])//10000+numpy.arange(rfpar["mmax"])/12.,rmeanw(rbreakts*1.0,3),label='rejected')
    for p in peaklist:
        plt.text(int(rfpar['startdate'])//10000+p[0]/12.,p[1],str(rfpar['startdate']//10000+p[0]//12)+"{:0>2}".format(numpy.mod(p[0],12)+1))
    plt.legend()
    plt.title(exp+' '+str(sum(abreakts))+' breaks')
    plt.xlim(1940,2016)
    print(('writing '+exp+'/breakanalysis.ps'))
    plt.savefig(plotproperties['tmppath']+'/breakanalysis.ps')
    plt.close()

def call_allrasotrends_presat(settingsf='settings',pardict={}):

    magicscolours=defcolours()
#    print('settings:',settingsf, 'pardict:',pardict)
    home=os.getenv('HOME')
    plottasks_file=open(home+'/python/allrasotrends/'+settingsf+'.json')
    settings=json.load(plottasks_file)#,encoding='latin-1')
    plottasks_file.close()
    if type(pardict) is str:
        try:
            
            with open(os.path.expandvars(os.path.expanduser(pardict))) as clf:
                pardict = json.load(clf)
        except:
            pass
        
    if pardict:
        for k,v in list(pardict.items()):
            print((k,v))
            settings['plotproperties'][k]=v
            
        with open(home+'/args.json', 'w') as clf:
            json.dump(pardict, clf)

    try:
        
        settings['plotproperties']['sonde_mask'] = np.string_(settings['plotproperties'][settings['plotproperties']['sonde_mask']])
    except:
        settings['plotproperties']['sonde_mask'] = []
        
    if type(settings['plotproperties']['plotstatids']) is list:
        settings['plotproperties']['plotstatids']=settings['plotproperties']['plotstatids'][0]

    settings['plotproperties']['contour_shade_colour_list']= ["rgb(0,0,0.3)", "rgb(0,0,0.5)",
                                                              "rgb(0,0,0.7)", "rgb(0,0,0.9)", "rgb(0,0.15,1)",
                                                              "rgb(0,0.3,1)", "rgb(0,0.45,1)", "rgb(0,0.6,1)",
                                                              "rgb(0,0.75,1)", "rgb(0,0.85,1)", "rgb(0.2,0.95,1)",
                                                              "rgb(0.45,1,1)", "rgb(0.75,1,1)", "rgb(1,1,0)",
                                                              "rgb(1,0.9,0)", "rgb(1,0.8,0)", "rgb(1,0.7,0)",
                                                              "rgb(1,0.6,0)", "rgb(1,0.5,0)", "rgb(1,0.4,0)",
                                                              "rgb(1,0.3,0)", "rgb(1,0.15,0)", "rgb(0.9,0,0)",
                                                              "rgb(0.7,0,0)", "rgb(0.5,0,0)", "rgb(0.3,0,0)"]

    f=open(os.path.expanduser('~/tables/guanstationslist.txt'))
    data=f.read().split('\n')
    gstations=[]
    for d in data:
        gstations.append('0'+d[:5])

    gstations=[]
    daynight='' # 'Night' 'Day'
    nens=16


    try:
        tmppath=copy.copy(settings['plotproperties']['tmppath'])
    except:
        tmppath='fastscratch'

    dplotproperties={}    
    dtasks={}
    iexp=-1
    for exp in settings['plotproperties']['exps']:
        iexp+=1
        path=home+'/fastscratch/rise/'+settings['plotproperties']['version']+'/'+exp+'/'    
        settings['plotproperties']['tmppath']=home+'/'+tmppath+'/rise/'+settings['plotproperties']['version']+'/'+exp+'/'

        try:
            nml = f90nml.read(path+'/'+'radcorpar')
        except:
            continue
        days=rstu.calcdays(nml['rfpar']['startdate'],nml['rfpar']['mmax'])

        dplotproperties[str(exp)]=copy.deepcopy(settings['plotproperties'])
        dplotproperties[str(exp)]['path']=str(path)
        plotproperties=dplotproperties[str(exp)]
        plotproperties['exps']=[str(exp)]
        try:
            plotproperties['logp']=eval(plotproperties['logp'])
        except:
            pass
        try:
            plotproperties['double']=eval(plotproperties['double'])
        except:
            pass
        try:
            plotproperties['version']=nml['rfpar']['version']
            plotproperties['fgdepname']=nml['rfpar']['fgdepname']
        except:
            plotproperties['fgdepname']='fg_dep'
            plotproperties['version']=path.split('/')[-2]


        ps=numpy.asarray((10,20,30,50,70,100,150,200,250,300,400,500,700,850,925,1000),dtype='float32')
        msups=numpy.asarray((550,220,90,880),dtype='float32')
        plotproperties['ps']=ps
        plotproperties['msups']=ps
        plotproperties['initial_size']=5000
        plotproperties['pindex']=numpy.asarray(plotproperties['pindex'],dtype='int')
        plotproperties['msupindex']=numpy.asarray(plotproperties['msupindex'],dtype='int')
        plotproperties['jpindex']=numpy.concatenate((plotproperties['pindex'],plotproperties['msupindex']))
        mode=plotproperties['mode']
        intervals=numpy.asarray(plotproperties['intervals'],dtype=numpy.int32)
        plotproperties['intervals']=intervals.copy()
        if sys.version_info[0]==2:
            strtype=[str,str]
        else:
            strtype=[str]
            
        if type(plotproperties['diffintervals']) in strtype:
            if plotproperties['diffintervals']=='all':
                diffintervals=[]
                for y in range(1939,2016,2):
                    tol=3
                    if y<1979:
                        tol=4
                    if y<1958:
                        tol=6
                    diffintervals.append([y,y+1,tol])
                diffintervals=numpy.asarray(diffintervals,dtype=numpy.int32)

            else:
                diffintervals=eval(plotproperties['diffintervals'])
        else:    
            diffintervals=numpy.asarray(plotproperties['diffintervals'],dtype=numpy.int32)
        plotproperties['diffintervals']=diffintervals.copy()
        plotinterval=numpy.asarray(plotproperties['plotinterval'],dtype=numpy.int32)
        plotproperties['plotinterval']=plotinterval


        stnames=eval(plotproperties['stnames'])
        if type(stnames[0]) is not str:
            print('please specify stnames as list of strings')
            return
        elif '-' in stnames[0]:
            plotproperties['super']=True
            plotproperties['superstring']=copy.copy(stnames[0])
            smin=int(stnames[0].split('-')[0])
            smax=int(stnames[0].split('-')[1])
            l=glob.glob(path+'??????/feedbackmerged??????.nc')
            nstnames=[]
            for ll in l:
                s=ll.split('/')[-1][-9:-3]
                try:
                    
                    ns=int(s)
                except:
                    continue
                if ns>=smin and ns<=smax:
                    nstnames.append(s)
            stnames=nstnames
        else:
            if plotproperties['super']:
                plotproperties['superstring']=stnames[0]+'-'+stnames[-1]



        print(('stnames after evaluation:',stnames))
        rfpar=nml['rfpar']
        rfpar['brmax']=100
        tmtest=rst.dailyts(rfpar,plotproperties['pindex'],stnames[0])

        # order is important tmcorr and tm must always be the first two

        for m in mode:
            if 'breaks'==m:

                break_analysis(exp,rfpar,plotproperties,stnames)

            elif 'dump'==m:

                tasks=[]
                tasknames=['rcorr','rio']
                tasks=tasks+define_ds(tasknames,rfpar,plotproperties,stnames,path)

                t=time.time()
                alldd=rsdump.readraobcoredump(path,stnames)
                rsdump.dumpanalysis(path,tasks,plotproperties,alldd)
                print((time.time()-t))

            elif 'spag' in m:

                tasks=[]
                tasknames=['rcorr']+plotproperties['spagtasks']

                tasks=tasks+define_ds(tasknames,rfpar,plotproperties,stnames,path)

                t=time.time()
                with open(path+'mergedstations.json','r') as f:
                    jdict=json.load(f) #,encoding='utf8')

                rcorrlist=list()
                for d in tasks:
                    d['dindex']=numpy.empty([d['ddata'][0],d['ddata'][4]],dtype='int')
                    d['ddata']=numpy.empty([d['ddata'][0],d['ddata'][1],d['ddata'][2],ps.shape[0],d['ddata'][4]])
                    rsdump.readrcorr(d,path,stnames,plotproperties=plotproperties)
                print(( 'vor readspagarr',time.time()-t))
                for spt in plotproperties['spagtasks']:
                    print((path,stnames,spt,[int(plotproperties['ens'][0])]))
                    spagdict=rsdump.readspagarr(path,stnames,spt,[int(plotproperties['ens'][0])])
                    print(('vor spaganalysis',spt,stnames,tasknames,time.time()-t))
    #                rsdump.spaganalysis(path,tasks,stnames,plotproperties,spagdict)
                    if m=='spag':
#			rsdump.spaganalysis_allbreaks(path,tasks,stnames,plotproperties,spagdict,spt,str(nml['rfpar']['startdate']/10000))
#		    if m=='spagmap':
                        rsdump.spagplusmaps(path,tasks,stnames,plotproperties,spagdict,jdict,spt,str(nml['rfpar']['startdate']//10000))
                        #rsdump.spagmap_allbreaks(path,tasks,stnames,plotproperties,spagdict,jdict,spt,str(nml['rfpar']['startdate']/10000))
                print((time.time()-t))

            elif 'monthlydiff'==m:
                t=time.time()
                try:
                    del tasks
                except:
                    pass
                dtasks[exp]=[]
                try:
                    tasknames=copy.copy(plotproperties['monthlytasks'])
                    for k in range(len(tasknames)):
                        tasknames[k]=str(tasknames[k])
                    if ('erai_fggpsdep' in tasknames or 'erai_fggpswetdep' in tasknames) and 'eijra_fgdep' not in tasknames:
                        tasknames=['eijra_fgdep']+tasknames
                    else:
                        try:
                            del tasknames[tasknames.index('eijra_fgdep')]
                            tasknames=['eijra_fgdep']+tasknames
                        except:
                            pass

                    if 'tm' not in tasknames:
                        if 'tmcorr' in tasknames:
                            tasknames.insert(1,'tm')
                        else:
                            tasknames=['tm']+tasknames
                    if 'tmcorr' not in tasknames:
                        tasknames=['tmcorr']+tasknames
                    else:
                        tasknames.remove('tmcorr')
                        tasknames=['tmcorr']+tasknames

                    if plotproperties['reference'][0] =='bg' and 'bg' not in tasknames:
                        tasknames.insert(2,'bg')
                    if 'tmcorr' in plotproperties['monthlytasks']:
                        tasknames.append('tmcorr')
                    if 'tm' in plotproperties['monthlytasks']:
                        tasknames.append('tm')

                except IOError:
                    print('monthly: no tasks defined')

                dtasks[exp]=dtasks[exp]+define_ds(tasknames,rfpar,plotproperties,stnames,path)
                tasks=dtasks[exp]

                print('initializing arrays...')
                tt=time.time()
                #for interv in  range(diffintervals.shape[0]): #range(diffintervals.shape[0]):
                interv=0
                print(('Interval:',diffintervals[interv], 'ini: ',time.time()-tt))
                rstd.allrasodiffs(path,tasks,plotproperties,diffintervals[interv],
                                  days,stnames,interv,gstations=gstations,daynight=daynight,
                                  ref=plotproperties['reference'][0],init=True)

            elif 'monthly'==m:
                t=time.time()
                try:
                    del tasks
                except:
                    pass
                dtasks[exp]=[]
                #tasknames=['tmcorr','tm','bgpresat','rio','ce20c0_andep','ce20cens_andep','e20c_andep','n20c_andep']#,bgpresat,ce20c0_andep]#,e20c_andep,n20c_andep,rio24,riomean,rss,uah] # SATSeries am Schluss!!#,bgdiff,rcorr
                try:
                    tasknames=plotproperties['monthlytasks']
                    for k in range(len(tasknames)):
                        tasknames[k]=str(tasknames[k])
                    if ('erai_fggpsdep' in tasknames or 'erai_fggpswetdep' in tasknames) and 'eijra_fgdep' not in tasknames:
                        tasknames=['eijra_fgdep']+tasknames
                    else:
                        try:
                            del tasknames[tasknames.index('eijra_fgdep')]
                            tasknames=['eijra_fgdep']+tasknames
                        except:
                            pass

                    if 'tm' not in tasknames:
                        if 'tmcorr' in tasknames:
                            tasknames.insert(1,'tm')
                        else:
                            tasknames=['tm']+tasknames
                    if 'tmcorr' not in tasknames:
                        tasknames=['tmcorr']+tasknames
                    else:
                        tasknames.remove('tmcorr')
                        tasknames=['tmcorr']+tasknames

                except IOError:
                    print('monthly: no tasks defined')

                dtasks[exp]=dtasks[exp]+define_ds(tasknames,rfpar,plotproperties,stnames,path)
                tasks=dtasks[exp]

                for ta in tasks:
                    ta['amps']=numpy.zeros((3,numpy.asarray(ta['ens']).shape[0],intervals.shape[0],6)) # trend amplification factors
                    ta['btm']=numpy.zeros((3,numpy.asarray(ta['ens']).shape[0],intervals.shape[0],6)) # belt trend maxima in upper troposphere
                    ta['had']=numpy.zeros((3,intervals.shape[0],6)) 


                goodinterval=0
                for interv in range(intervals.shape[0]):
                    print(('Interval:',intervals[interv]))
                    if goodinterval==0:
                        try:
                            hadmed,hadtem,hadens,sats=rsta.allrasotrends(path,tasks,plotproperties,intervals[interv],
                                                                         days,stnames,interv,gstations=gstations,daynight=daynight,init=True)
                            goodinterval=1
                        except MemoryError:
                            pass
                    else:
                        rsta.allrasotrends(path,tasks,plotproperties,intervals[interv],
                                           days,stnames,interv,hadmed,hadtem,hadens,sats,gstations=gstations,daynight=daynight)

                print((time.time()-t))
                if sys.version_info[0]<4:
                    if 'belts' in plotproperties['plotlist'] and plotproperties['pindex'].shape[0]>1:
                        centers=0.5*(intervals[:,0]+intervals[:,1])
                        beltnames=["Globe",'NHEx','SHEx','Tropics','NH','SH']
                        ti=list(range(len(tasks)))
                        li=[ta['shortname'] for ta in tasks]
                        if 'rio' in li:
                            ti.insert(0,ti.pop(li.index('rio')))
                        if 'ce20c_andep' in li:
                            ti.insert(0,ti.pop(li.index('ce20c_andep')))
    
                        for ibelt in range(len(beltnames)):
                            entries=[]
                            h=[]
                            for ii in ti:
                                ta=tasks[ii]
                                amps=ta['amps'][:,:,:,ibelt]
                                col=magicscolours[str(ta['color'])]
                                lab=ta['name']
                                if '_andep' in lab:
                                    lab=lab.split('_andep')[0]
    
                                if amps.shape[1]>1:
                                    for iens in range(amps.shape[0]):
                                        for iq in range(ta['had'].shape[0]):
                                            
                                            amps[iq,iens,ta['had'][iq,:,ibelt]<0.05]=numpy.nan
    
                                    amin=numpy.min(amps,axis=1)
                                    amax=numpy.max(amps,axis=1)
                                    amean=numpy.mean(amps,axis=1)
                                    l=0
                                    for iint in centers[0:-1]:
                                        di=centers[l+1]-centers[l]
                                        xy=list(zip([iint,iint+di,iint+di,iint,iint],[amin[l],amin[l+1],amax[l+1],amax[l],amin[l]]))
                                        l+=1
                                        mp=mpatches.Polygon( xy, facecolor=col, edgecolor=col,label=lab,alpha=0.7)
                                        plt.gca().add_patch(mp)
                                    if 'mp' in locals():
                                        h.append(mp)
                    #                plt.plot(intervals[:,0],amean,'o-')
                                    entries.append(lab)
                                else:
                                    for iq in range(ta['had'].shape[0]):
                                        amps[iq,0,ta['had'][iq,:,ibelt]<0.05]=numpy.nan
                                    pl=plt.plot(centers,amps[1,0,:],'o-',ms=5,lw=2,color=col,
                                                markerfacecolor=col,markeredgecolor=col,label=lab)
                                    h.append(pl[0])
                                    entries.append(lab)
    
    
                            plt.xlim([intervals[0,0]+9,intervals[-1,1]-9])
                            plt.ylabel('Amplification Factor')
                            plt.ylim([-0.5,4])
                            ax=plt.gca()
                            ax2= ax.twinx()
                            ax2.plot(centers,ta['had'][1,:,ibelt]-ta['had'][1,:,ibelt],':')
                            pl=plt.plot(centers,ta['had'][1,:,ibelt],'k-',label='HadCRUT5 trend')
                            h.append(pl[0])
                            ax2.set_ylim([-0.5*0.6,4*0.6])
                            ax2.set_ylabel('K/10a')
                            amin=ta['had'][0,:,ibelt]
                            amax=ta['had'][2,:,ibelt]
                            l=0
                            for iint in centers[0:-1]:
                                di=centers[l+1]-centers[l]
                                xy=list(zip([iint,iint+di,iint+di,iint,iint],
                                            [amin[l],amin[l+1],amax[l+1],amax[l],amin[l]]))
                                l+=1
                                mp=mpatches.Polygon( xy, facecolor=col, edgecolor=col,label=lab,alpha=0.7)
                                ax2.add_patch(mp)
                            #if 'mp' in locals(): # not needed
                                #h.append(mp)
                    #        plt.legend(entries)
                            plt.legend(handles=h,loc=2)
                            profsat=''
                            if 'satseries' in plotproperties['plotlist']:
                                profsat='sat'
                                
                            plt.tight_layout()
                            plt.savefig(plotproperties['tmppath']+'/'+profsat+'amplifications'+beltnames[ibelt]+'{}-{}'.format(intervals[0,0],intervals[-1,1])+'.pdf')
                            plt.close()
                        # trends only - no amplification   
                        # not yet implemented
                        pass
            elif 'daily' in m:
# 81: "[\"094302\",\"094312\",\"094430\",\"094430\",\"094510\",\"094637\",\"094647\",\"094659\",\"094711\"]"
# 80: "[\"094120\",\"094150\",\"094203\",\"094294\",\"094299\",\"094374\",\"094403\",\"094461\",\"094578\",\"094638\",\"094672\",\"094802\",\"094821\",\"094866\",\"094910\",\"094975\",\"094996\"]"
                tasks=[]
                tasknames=['obs','bgdep','eijra_fgdep','e20c_andep','n20c_andep','ce20c0_andep','jra55_fgdep']#,bgpresat,ce20c0_andep]#,e20c_andep,n20c_andep,rio24,riomean,rss,uah] # SATSeries am Schluss!!#,bgdiff,rcorr
                tasknames=plotproperties['dailytasks']#,bgpresat,ce20c0_andep]#,e20c_andep,n20c_andep,rio24,riomean,rss,uah] # SATSeries am Schluss!!#,bgdiff,rcorr

                corrlist=['rcorr']
                for i in range(32):
                    corrlist.append('era5v{:0>1}'.format(i))
                    corrlist.append('rio{:0>2}'.format(i))
                    corrlist.append('rit{:0>2}'.format(i))
                for corr in corrlist:
                    for ta in tasknames:
                        tx=ta.split(corr)[0].split('-')[0]
                        if corr in ta and tx!=ta:
                            if corr not in tasknames:
                                plotproperties['dailytasks']=[corr]+plotproperties['dailytasks']
                            if tx not in tasknames and tx!='':
                                plotproperties['dailytasks']=[tx]+plotproperties['dailytasks']

                tasknames=plotproperties['dailytasks']
                print(('dailytasks:',tasknames))
                tasks=tasks+define_ds(tasknames,rfpar,plotproperties,stnames,path)

                for ta in tasks:
                    ta['ddata']=numpy.empty(ta['ddata'],dtype=numpy.float32)
                    ta['ddata'].fill(numpy.nan)
                    ta['parameter']=plotproperties['parameters'][0]

            #    for interv in range(intervals.shape[0]):
                for interv in range(1):
                    if intervals.ndim>1:
                        plotproperties["intervals"]=numpy.asarray(intervals[interv],dtype=numpy.float32)
                    t=time.time()
            #        dumpdict=readraobcoredump(path,stnames)
                    if m=='dailyscatter':
                        rsts.statscatter(path,tasks,plotproperties,nml,stnames)
                    else:
                        rsts.statanalysis(path,tasks,plotproperties,nml,stnames,minlen=10)

            #        rsts.statanalysis(path,tasks,plotproperties,stnames,minlen=360,ref='../exp05/')
                    print(('interval: ',time.time()-t))
            elif 'obstat'==m:
                tasks=[]
                tasknames=['obs','bgdep','eijra_fgdep','e20c_andep','n20c_andep','ce20c0_andep','jra55_fgdep']
                tasknames=plotproperties['dailytasks']

                corrlist=['rcorr']
                for i in range(32):
                    corrlist.append('rio{:0>2}'.format(i))
                    corrlist.append('rit{:0>2}'.format(i))
                for corr in corrlist:
                    for ta in tasknames:
                        tx=ta.split(corr)[0]
                        if corr in ta and tx!=ta:
                            if corr not in tasknames:
                                plotproperties['dailytasks']=[corr]+plotproperties['dailytasks']
                            if tx not in tasknames and tx!='':
                                plotproperties['dailytasks']=[tx]+plotproperties['dailytasks']

                tasknames=plotproperties['dailytasks']
                print(('dailytasks:',tasknames))
                tasks=tasks+define_ds(tasknames,rfpar,plotproperties,stnames,path)

                for ta in tasks:
                    ta['ddata']=numpy.empty(ta['ddata'])

                for interv in range(1):
                    if intervals.ndim>1:
                        plotproperties["intervals"]=numpy.asarray(intervals[interv],dtype=numpy.float)
                    t=time.time()
                    obstat.obstat(path,tasks,plotproperties,nml,stnames,minlen=360)

                    print(('interval: ',time.time()-t))

            else:
                print(('no valid mode:',mode))

            print('finished')

    print('exp finished')
    return
    if mode[0] in ['monthly','monthlydiff']:
        settings['plotproperties']['refexp']='Control'
    #    settings['plotproperties']['cvar']='tmcorr'
        plist=[]
        varlist=[]
        namelist=[]
        clist=[]
        expnames=dict(Control='RC1.5',exp01='EIJRAPRE',exp02='JRACE20C',exp04='JRACEPRE',
                      rss='RSS',uah='UAH',star='STAR',wegc='WegC',tm='RAW',tmcorr='RAOBCORE',rio='RICH',rit='RICHT',bg='BG',
                      merra='MERRA2',era5='ERA5')
        for i in range(32):
            expnames['rio{:0>2}'.format(i)]='RIO{:0>2}'.format(i)
            expnames['rit{:0>2}'.format(i)]='RIT{:0>2}'.format(i)
        beltnames=['Globe','NH','SH','Tropics','NHEx','SHEx','Cost']
        ibelt=beltnames.index(plotproperties['belt'][0])

        j=1
        if len(settings['plotproperties']['exps'])>0:
            ie=0
            for exp in settings['plotproperties']['exps']:
                for it in range(len(dtasks[exp])):
                    if ie==0:
                        plist.append(ps[plotproperties['pindex']])
                        if beltnames[ibelt]=='Cost':
                            varlist.append(dtasks[exp][it]['cost'][0,2,:])
                        else:
                            varlist.append(dtasks[exp][it]['beltslopes'][0,ibelt,:])

                        if numpy.sum(~numpy.isnan(varlist[-1]))==0:
                            varlist.pop()
                            continue
                        if dtasks[exp][it]['shortname']==settings['plotproperties']['cvar']:
                            namelist.append(expnames[exp])
                        else:
                            namelist.append(expnames[dtasks[exp][it]['shortname']])
                        clist.append(dtasks[exp][it]['color'])
                    else:
                        if dtasks[exp][it]['shortname']==settings['plotproperties']['cvar']:
                            plist.append(ps[plotproperties['pindex']])
                            if beltnames[ibelt]=='Cost':
                                varlist.append(dtasks[exp][it]['cost'][0,2,:])		
                            else:
                                varlist.append(dtasks[exp][it]['beltslopes'][0,ibelt,:])		

                            if numpy.sum(~numpy.isnan(varlist[-1]))==0:
                                varlist.pop()
                                continue
                            namelist.append(expnames[exp])
                            rms=[]
                            mgk=list(magicscolours.keys())
                            for k in range(len(mgk)):
                                rms.append(numpy.mean((magicscolours[dtasks[exp][it]['color']]-magicscolours[mgk[k]])**2))
                            idx=numpy.argsort(rms)
                            while rms[idx[j]]==0:
                                j+=1
                            clist.append(mgk[idx[j]])
                            j+=1
                ie+=1
                magicscolours[dtasks[exp][it]['color']]

            if 'belts' in plotproperties['plotlist']:
                if 'cost' in plotproperties['plotlist']:
                    compare_profiles(plist,varlist,namelist,clist,plotproperties,expnames[settings['plotproperties']['cvar']],beltname=plotproperties['belt'][0])
                else:    
                    compare_profiles(plist,varlist,namelist,clist,plotproperties,expnames[settings['plotproperties']['cvar']],beltname=plotproperties['belt'][0],had=tasks[0]['had'][0,ibelt])
            #if 'cost' in plotproperties['plotlist']:
                #compare_profiles()
            if 'beltanomalies' in plotproperties['plotlist']:
                compare_tseries()
            if 'zmeans' in plotproperties['plotlist']:
                compare_zmeans()


@njit
def ffff(b):
    a = numpy.zeros(5,dtype=int64)
    return a


if __name__ == '__main__':
 
    #os.chdir('/fio/srvx7/leo/fastscratch/rise/1.0/exp04')
    #for f in glob.glob('??????/feedbackmerged??????.nc'):
        #rsts.add_meanandspread(f)
    t=time.time()
    #try:
    if True:
        if len(sys.argv) == 2:
            pardict = {}
            sys.argv[2] = ''
        elif len(sys.argv)==3:
            if '.json' in sys.argv[2]:
                pardict = sys.argv[2]
            else:
                pardict=eval(sys.argv[2])
                if '[0' in pardict['stnames']:
                    pardict['stnames']="['"+pardict['stnames'][1:7]+"']"
        print('called with ',sys.argv)
        call_allrasotrends_presat(settingsf=sys.argv[1],pardict=pardict)
    #except KeyError:
        #setf='settings'
        #print('called with ',setf)
        #call_allrasotrends_presat(settingsf=setf)

    print(('grand total:',time.time()-t))
    pass
