#from Magics.macro import *
from rasotools.additions.utils import *
import netCDF4
import time
#from mpl_toolkits.basemap import Basemap
#import  matplotlib.pylab as plt
#import matplotlib.patches as mpatches
#import matplotlib.cm as cm
#import matplotlib.colors as colors



def readspagarr(path,stnames,spt,ens):

    allspag=dict()
    for st in stnames:
        allspag[st]=dict()
        for it in ['1','2']:
            allspag[st][it]=list()
            for mem in ens:
                if spt=='rio':
                    spd={'stname':st,'type':'obs','iens':mem}
                else:
                    spd={'stname':st,'type':'','iens':mem}
                fn=path+'/'+st+'/'+spd['type']+'spagarr_'+it+'_{0:0>2}_'.format(mem)+st
                try:
                    f=open(fn,'rb')
                    recl = numpy.fromfile(f, dtype='uint32', count=1)
                    tbi,l,pmax,parmax = numpy.ndarray.byteswap(numpy.fromfile(f, dtype='int32', count=4),True)
                    spd['tbi'],spd['l'],spd['pmax'],spd['parmax']=tbi,l,pmax,parmax
                    recl = numpy.fromfile(f, dtype='uint32', count=2)
                    spd['tbindex']=numpy.ndarray.byteswap(numpy.fromfile(f, dtype='int32', count=tbi),True)
                    recl = numpy.ndarray.byteswap(numpy.fromfile(f, dtype='uint32', count=2))
                    spd['wmonrs']=numpy.ndarray.byteswap(numpy.fromfile(f, dtype='int32', count=tbi*l),True)
                    spd['wmonrs']=numpy.transpose(numpy.reshape(spd['wmonrs'],[l,tbi]))
                    recl = numpy.ndarray.byteswap(numpy.fromfile(f, dtype='uint32', count=2))
                    spd['xdarr']=numpy.ndarray.byteswap(numpy.fromfile(f, dtype='float32', count=parmax*pmax*tbi),True)
                    spd['xdarr']=numpy.transpose(numpy.reshape(spd['xdarr'],[parmax,pmax,tbi]))
                    spd['xdarr'][spd['xdarr']==-999.]=numpy.nan
                    recl = numpy.ndarray.byteswap(numpy.fromfile(f, dtype='uint32', count=2))
                    n = numpy.ndarray.byteswap(numpy.fromfile(f, dtype='int32', count=1),True)
                    spd['n'] = n
                    recl = numpy.fromfile(f, dtype='uint32', count=2)

                    spd['spagindex']=numpy.ndarray.byteswap(numpy.fromfile(f, dtype='int32', count=n),True)
                    recl = numpy.ndarray.byteswap(numpy.fromfile(f, dtype='uint32', count=2))
                    spd['spagarrc']=numpy.ndarray.byteswap(numpy.fromfile(f, dtype='float32', count=n),True)
                    spd['spagarrc'][spd['spagarrc']==-999.]=numpy.nan

                    allspag[st][it].append(spd.copy())

                    #spagarr=numpy.empty([l*pmax*parmax*tbi],dtype='float32')
                    #spagarr.fill(numpy.nan)
                    #spagarr[spd['spagindex']]=spd['spagarrc']
                    #spagarr=numpy.reshape(spagarr,[l,parmax,pmax,tbi])
                except:
                    pass
                #for i in range(tbi):
                    #print spd['tbindex'][i],1900+spd['tbindex'][i]/365.25
                    #for j in range(spagarr.shape[0]):
                        #plt.plot(spagarr[j,0,:,i],numpy.arange(pmax),'k')
                        #plt.ylim([15,0])
                        #plt.plot(spd['xdarr'][i,:,0],numpy.arange(pmax),'r')
                    #plt.show()

    return allspag

def readraobcoredump(path,stnames):

    alldd=dict()
    for st in stnames:
        dd=dict()
        fn=path+'/'+st+'/'+st+'.dump'
        f=open(fn,'rb')
        recl = numpy.fromfile(f, dtype='uint32', count=1)
        nmax,il,pmax,parmax,probmax = numpy.ndarray.byteswap(numpy.fromfile(f, dtype='int32', count=5),True)
        dd['nmax'],dd['il'],dd['pmax'],dd['parmax'],dd['probmax']=nmax,il,pmax,parmax,probmax
        recl = numpy.fromfile(f, dtype='uint32', count=2)
        dd['goodindex']=numpy.ndarray.byteswap(numpy.fromfile(f, dtype='int32', count=il),True)
        recl = numpy.ndarray.byteswap(numpy.fromfile(f, dtype='uint32', count=2))
        dd['tsaofmean']=numpy.ndarray.byteswap(numpy.fromfile(f, dtype='float32', count=il*probmax),True)
        dd['tsaofmean']=numpy.transpose(numpy.reshape(dd['tsaofmean'],[probmax,il]))
        dd['tsaofmean'][dd['tsaofmean']==-999.]=numpy.nan
        recl = numpy.ndarray.byteswap(numpy.fromfile(f, dtype='uint32', count=2))
        dd['aprioriprobs']=numpy.ndarray.byteswap(numpy.fromfile(f, dtype='float32', count=il),True)
        dd['aprioriprobs'][dd['aprioriprobs']==-999.]=numpy.nan
        recl = numpy.ndarray.byteswap(numpy.fromfile(f, dtype='uint32', count=2))
        dd['breakmeanprobs']=numpy.ndarray.byteswap(numpy.fromfile(f, dtype='float32', count=il*probmax),True)
        dd['breakmeanprobs']=numpy.transpose(numpy.reshape(dd['breakmeanprobs'],[probmax,il]))
        dd['breakmeanprobs'][dd['breakmeanprobs']==-999.]=numpy.nan
        recl = numpy.ndarray.byteswap(numpy.fromfile(f, dtype='uint32', count=2))
        dd['chosenbreaks']=numpy.ndarray.byteswap(numpy.fromfile(f, dtype='int32', count=30),True)
        alldd[st]=dd.copy()

    return alldd

def readrcorr(d,path,stnames,plotproperties=None):

    dens=d['ens']
    try:
        if len(dens)>1:
            dens=[int(plotproperties['ens'][0])]
    except:
        pass
    rcorrdict={}
    d['stnames']=list()
    l=0
    for statid in stnames:
        for ens in dens:
            if "dfile" in d:
                prefix=d["dfile"][:]
            else:
                prefix=d["file"][:]

            if "prefix" in d:
                relpath=d["prefix"][:]
            else:
                relpath="./"

            if "dsuff" in d:
                suffix=d["dsuff"][:]
            else:
                suffix=['']

            fn=list()
            for ifile in range(len(prefix)):
                if 'rio' not in d['shortname'] and 'rit' not in d['shortname']:
                    fn=fn+[os.path.join(path,relpath,statid,prefix[ifile]+statid+suffix[ifile]+".nc")]
                else:
                    fn=fn+[os.path.join(path,relpath,statid,prefix[ifile]+"{0:0>2}".format(ens)+'_'+statid+suffix[ifile]+".nc")]

                if os.path.isfile(fn[ifile]):


                    try:
                        f = netCDF4.Dataset(fn[ifile],"r")
                        f.set_auto_mask(False)
                        dvshape=f.variables[d['dvar']].shape
                        d['dindex'][l,:dvshape[2]]=f.variables['datum'][0,:]
                        d['ddata'][l,0,:,:,:dvshape[2]]=f.variables[d['dvar']][:]
                        d['stnames'].append(statid)
                        l+=1
                    except:
                        pass


    return 

def dumpanalysis(path,tasks,plotproperties,alldd):

    try:
        dummy=plotproperties['dummy']
    except:
        dummy=''
    for stname in list(alldd.keys()):
        t=time.time()
        mktmppath(plotproperties['tmppath']+'/'+stname)
        oname=plotproperties['tmppath']+'/'+stname+'/'+stname+'_breakanalysis'+dummy
        poutput = output(output_name = oname, 
                         output_formats = plotproperties["outputformat"],
                         output_title= oname,
                         output_name_first_page_number = "off", 
                         super_page_x_length=21.,super_page_y_length=29., )

        plotlist=[poutput]
        dd=alldd[stname]
        gi=dd['goodindex']
        nmax=dd['nmax']
        y=numpy.empty(nmax)
        y.fill(numpy.nan)
        x=numpy.floor(tasks[0]['startdate']//10000)+numpy.arange(nmax)/365.25

        ypos=0.
        dypos=2.5
        y[gi]=dd['aprioriprobs']
        plotlist=plotlist+dumpseries(x,y,dd['chosenbreaks'],plotproperties,title='aprioriprobs',xlab='on')

        titles=['strat00','strat12','trop00','trop12','rad','fg','stratcomp00','stratcomp12','tropcomp00','tropcomp12']
        l=0
        for ipar in range(dd["probmax"]-3,-1,-1):
            for par in ['tsaofmean','breakmeanprobs']:
                y.fill(numpy.nan)
                y[gi]=dd[par][:,ipar]
                axpos='right'
                col='black'
                if par=='tsaofmean':
                    axpos='left'
                    col='red'
                    ypos+=dypos
                print((titles[len(titles)-1-l]))
                plotlist=plotlist+dumpseries(x,y,dd['chosenbreaks'],plotproperties,title=titles[len(titles)-1-l],
                                             ypos=ypos,axpos=axpos,col=col,xlab='off')
            l+=1

        plot(plotlist)
        print((oname,time.time()-t))

def dumpseries(x,y,chosenbreaks,plotproperties,title='',ypos=0.,axpos='left',col='red',xlab='on'):

    index=numpy.zeros(y.shape[0],dtype=numpy.int32)
    n=thin(y,index,10)
    y[numpy.isnan(y)]=-1.e21

    data = minput(Input_x_values=x[index[0:n]],
                  Input_y_values=y[index[0:n]],
                  Input_y_missing_value=-1.e21
                  )

    ppage = page(
        layout='positional',  
        page_x_length=21., 
        page_y_length=2., 
        page_id_line='off',
        page_x_position=0., 
        page_y_position=ypos
    )

    dynrange=[0.,1.]
    y=y[y>-1.e21]
    if(len(y)>0):
        m=max(y)
        if(m>dynrange[1]):
            dynrange[1]=m
        m=min(y)
        if(m<dynrange[0]):
            dynrange[0]=m

    projection = mmap(
        subpage_map_projection='cartesian',
        subpage_x_axis_type='regular',
        subpage_y_axis_type='regular',   
        subpage_x_position=1.5,
        subpage_y_position=1.0,
        subpage_x_length=17.5,
        subpage_y_length=2.0,
        subpage_x_max=float(plotproperties["plotinterval"][1]),
        subpage_x_min=float(plotproperties["plotinterval"][0]),
        subpage_y_max=dynrange[1],
        subpage_y_min=dynrange[0],
    )

    # Vertical axis

    vertical = maxis(
        axis_orientation='vertical',
        axis_grid='off',
        axis_type='regular',
        axis_tick_label_height=0.4,
        axis_grid_reference_level=0.,
        axis_grid_reference_thickness=1,
        axis_grid_reference_colour='black',
        axis_title='on',
        axis_position=axpos,
        axis_title_text='',
        axis_title_height=0.4,
    )

    # Horizontal axis

    horizontal = maxis(
        axis_orientation='horizontal',
        axis_type='regular',
        axis_tick_label_height=0.4,
        axis_tick_label_colour='charcoal',
        axis_tick_label=xlab,
        axis_grid='off',
        axis_grid_colour='charcoal',
        axis_grid_thickness=1,
        axis_grid_line_style='dash',
    )

    text=mtext({
        "text_lines" : [title],
        "text_html" : "true",
        "text_colour" : "black",
        "text_font_size" : 0.4,
        "text_mode" : "positional",
        "text_box_x_position": 1.7,
        "text_box_y_position": 2.4,
        "text_box_x_length": 5.,
        "text_box_y_length": 0.4,
        "text_border": "off",
        "text_justification": "left",})

    graph = mgraph( legend='off' ,
                    graph_line_colour=col, 
                    graph_line_thickness= 3,
                    )

    chosen = minput(Input_x_values=1900.+chosenbreaks[chosenbreaks>0]/365.25,
                    Input_y_values=chosenbreaks[chosenbreaks>0]-chosenbreaks[chosenbreaks>0]+0.95,
                    )
#    print 'chosen:',chosenbreaks[chosenbreaks>0]
    symb=msymb(symbol_marker_index = 27,
               #symbol_advanced_table_min_value = -2.4,
               #symbol_advanced_table_max_value = 0.4,
               legend = "off",
               symbol_type = "marker",)

    return [ppage,projection,horizontal,vertical,data,graph,chosen,symb,text]

def spaganalysis(path,tasks,stnames,plotproperties,spagdict):

    global context
    global actions
    ps=plotproperties['ps']
    for stname in list(spagdict.keys()):
        mktmppath(plotproperties['tmppath']+'/'+stname)
        for it in list(spagdict[stname].keys()):
            for iens in range(len(spagdict[stname][it])):
                spd=spagdict[stname][it][iens]
                l,pmax,parmax,tbi=spd['l'],spd['pmax'],spd['parmax'],spd['tbi']
                spagarr=numpy.empty([l*pmax*parmax*tbi],dtype='float32')
                spagarr.fill(numpy.nan)
                spagarr[spd['spagindex']]=spd['spagarrc']
                spagarr=numpy.reshape(spagarr,[l,parmax,pmax,tbi])
                wmonrs=spd['wmonrs']
                for ib in range(1,spd['tbindex'].shape[0]-1):
                    ym=toyearmonth(numpy.asarray([spd['tbindex'][ib]],dtype='int'),'days since 1900-01-01')[0]
                    out = output({"output_name_first_page_number":"off",
                                  "output_formats":plotproperties["outputformat"], 
                                  'output_name': plotproperties['tmppath']+'/'+'spag_'+stname+'_'+it+'_{0}'.format(ym),
                                  'page_id_line':'off', 
                                  'super_page_x_length':29.,'super_page_y_length':14., 
})

                    glist=[out]
                    legend_user_lines=['SPAG','RAOB','RICH']
                    validprof=False
                    if ym==196912:
                        pass
                    for ipar in range(parmax):
                        lines =["Break profile estimates, "+stname,'{0}'.format(ym)+',{0:0>2}h'.format(ipar*12)]
                        glist+=set_spagarr(spagdict,lines,xlim=[-4,4],xpos=14.*ipar,legend_user_lines=legend_user_lines) 
                        leg=glist[-1]
                        del glist[-1]
#			for i in range(spagarr.shape[0]):
                        tlist=list()
                        for ispag in range(spagarr.shape[0]):
                            if sum(numpy.isnan(spagarr[ispag,ipar,:,ib]))<pmax-1:
                                tlist.append('{0:0>6}'.format(wmonrs[ib,ispag]))

                        sv=spagarr[:,ipar,:,ib].flatten()
                        pps=numpy.tile(ps,spagarr.shape[0]).flatten()
                        glist+=add_mprofile(sv,pps,'black',tlist=tlist,thickness=1)
                        #glist+=add_mprofile(spagarr[i,ipar,:,ib],ps,'black')
                        #t=time.time()
                        #plt.plot(sv)
                        #plt.savefig('spag_'+stname+'_'+it+'_{0}'.format(ym)+'.eps')
                        #plt.close()
                        #print 'plt:',time.time()-t,sv.size

                        for d in tasks:
#			    ld=d['stnames'].index(stname)
                            ld=stnames.index(stname)
                            irb=numpy.asarray(numpy.where(abs(spd['tbindex'][ib]-d['dindex'][ld,:])<=1))
                            if irb.size>0:
                                if irb[0][0]>0:
                                    glist+=add_mprofile(d['ddata'][ld,iens,ipar,:,irb[0][0]]-d['ddata'][ld,iens,ipar,:,irb[0][0]-1],
                                                        ps,d['color'],thickness=12)
                                    validprof=True
                        glist+=[leg]

                        if validprof:
                            t=time.time()
                            plot(glist)
                            try:
                                del context
                                del actions
                            except:
                                pass
                            print((time.time()-t,sv.size))

    return

def spaganalysis_allbreaks(path,tasks,stnames,plotproperties,spagdict,spt,startyear):

    global context
    global actions
    ps=plotproperties['ps']
    if 'dummy' not in list(plotproperties.keys()):
        plotproperties['dummy']=''
    #print spagdict.keys(),plotproperties['iteration']

    for stname in list(spagdict.keys()):
        for it in plotproperties['iteration']:
            if len(spagdict[stname][it])==0:
                continue
            mktmppath(plotproperties['tmppath']+'/'+stname)
            #print spagdict[stname][it][0].keys()
            for iens in range(len(spagdict[stname][it])):
                spd=spagdict[stname][it][iens]
                #print 'spd',spd
                l,pmax,parmax,tbi=spd['l'],spd['pmax'],spd['parmax'],spd['tbi']
                spagarr=numpy.empty([l*pmax*parmax*tbi],dtype='float32')
                spagarr.fill(numpy.nan)
                spagarr[spd['spagindex']]=spd['spagarrc']
                spagarr=numpy.reshape(spagarr,[l,parmax,pmax,tbi])
                wmonrs=spd['wmonrs']
                for ipar in range(len(plotproperties['time'])):
                #for ib in range(1,spd['tbindex'].shape[0]-1):
                    out = output({"output_name_first_page_number":"off",
                                  "output_formats":plotproperties["outputformat"], 
                                  'output_name': plotproperties['tmppath']+'/'+stname+'/spag_'+stname+'_'+spt+plotproperties['ens'][0]+'_'+it+'_'+plotproperties['time'][ipar]+plotproperties['dummy'],
                                  'page_id_line':'off', 
                                  'super_page_x_length':14.,'super_page_y_length':float((spd['tbindex'].shape[0]-1)/2)*7.+2., 
})

                    print(('out:',out.args['output_name']))  
                    glist=[out]
                    if spd['type']=='obs':
                        legend_user_lines=['SPAG','RAOB','RICH']
                    else:
                        legend_user_lines=['SPAG','RAOB','RICH-tau']
                    validprof=False
                    #for ipar in range(parmax):
                    for ib in range(1,spd['tbindex'].shape[0]-1):
                        ym=toyearmonth(numpy.asarray([spd['tbindex'][ib]],dtype='int'),'days since '+startyear+'-01-01')[0]
                        lines =["Break profile estimates, "+stname+', {0:0>2}h, '.format(ipar*12)+plotproperties['ens'][0]]
                        superpagetitle=ib==spd['tbindex'].shape[0]-2

                        glist+=set_spagarr(spagdict,lines,ym,xlim=[-4,4],xpos=7.*numpy.mod(ib-1,2),ypos=7.*((ib-1)/2),pxl=6.,pyl=7.,
                                           legend_user_lines=legend_user_lines,superpagetitle=superpagetitle) 
                        leg=glist[-1]
                        del glist[-1]
#			for i in range(spagarr.shape[0]):
                        tlist=list()
                        try:
                            if plotproperties['tlist']:
                                for ispag in range(spagarr.shape[0]):
                                    if sum(numpy.isnan(spagarr[ispag,ipar,:,ib]))<pmax-1:
                                        tlist.append('{0:0>6}'.format(wmonrs[ib,ispag]))
                        except:
                            pass

                        splist=[]
                        for k in range(spagarr.shape[0]):
                            if sum(spagarr[k,ipar,:,ib]==spagarr[k,ipar,:,ib]):
                                splist.append(spagarr[k,ipar,:,ib])
                        sv=numpy.asarray(splist).flatten()
                        if sv.size>0:
                        #sv=spagarr[:,ipar,:,ib].flatten()
                            pps=numpy.tile(ps,len(splist)).flatten()
                            glist+=add_mprofile(sv,pps,'black',tlist=tlist,pyl=7.0,thickness=1,superpage=superpagetitle)
                        #glist+=add_mprofile(spagarr[i,ipar,:,ib],ps,'black')
                        #t=time.time()
                        #plt.plot(sv)
                        #plt.savefig('spag_'+stname+'_'+it+'_{0}'.format(ym)+'.eps')
                        #plt.close()
                        #print 'plt:',time.time()-t,sv.size

                        idt=0
                        for d in tasks:
#			    ld=d['stnames'].index(stname)
                            ld=stnames.index(stname)
                            irb=numpy.asarray(numpy.where(abs(spd['tbindex'][ib]-d['dindex'][ld,:])<=1))
                            if irb.size>0:
                                if irb[0][0]>0:
                                    if sv.size==0 and idt==0:
                                        glist+=add_mprofile(d['ddata'][ld,iens,ipar,:,irb[0][0]]-d['ddata'][ld,iens,ipar,:,irb[0][0]-1],
                                                            ps,'black',tlist=tlist,pyl=7.0,thickness=1,superpage=superpagetitle)
                                        idt=1

                                    glist+=add_mprofile(d['ddata'][ld,iens,ipar,:,irb[0][0]]-d['ddata'][ld,iens,ipar,:,irb[0][0]-1],
                                                        ps,d['color'],pyl=7.0,thickness=12,superpage=superpagetitle)
                                    validprof=True
                                    glist+=[leg]

                        if validprof:
                            if superpagetitle:
                                legend = mlegend( {"legend":"on", 
                                                   "legend_text_font_size":0.5, 
                                                   "legend_text_colour":"black",
                                                   "legend_box_mode":"positional",
                                                   "legend_box_x_position":-1.0,
                                                   "legend_box_y_position":7.0-0.4,
                                                   "legend_box_x_length":6.0,
                                                   "legend_box_y_length":1.2,
                                                   "legend_text_composition":"user_text_only",
                                                   "legend_user_lines":legend_user_lines
                                                   })
                                glist+=[legend]

                    t=time.time()
                    plot(glist)
                            #try:
                                #del context
                                #del actions
                            #except:
                                #pass
                    print((time.time()-t,sv.size))

    return

def spagplusmaps(path,tasks,stnames,plotproperties,spagdict,jdict,spt,startyear):

    global context
    global actions
    ps=plotproperties['ps']
    if 'dummy' not in list(plotproperties.keys()):
        plotproperties['dummy']=''
    #print spagdict.keys(),plotproperties['iteration']

    for stname in list(spagdict.keys()):
        mktmppath(plotproperties['tmppath']+'/'+stname)
        for it in plotproperties['iteration']:
            if len(spagdict[stname][it])==0:
                continue
            #print spagdict[stname][it][0].keys()
            for iens in range(len(spagdict[stname][it])):
                spd=spagdict[stname][it][iens]
                #print 'spd',spd
                l,pmax,parmax,tbi=spd['l'],spd['pmax'],spd['parmax'],spd['tbi']
                spagarr=numpy.empty([l*pmax*parmax*tbi],dtype='float32')
                spagarr.fill(numpy.nan)
                spagarr[spd['spagindex']]=spd['spagarrc']
                spagarr=numpy.reshape(spagarr,[l,parmax,pmax,tbi])
                wmonrs=spd['wmonrs']
                for ipar in range(len(plotproperties['time'])):
                #for ib in range(1,spd['tbindex'].shape[0]-1):
                    out = output({"output_name_first_page_number":"off",
                                  "output_formats":plotproperties["outputformat"], 
                                  'output_name': plotproperties['tmppath']+'/'+stname+'/spag_'+stname+'_'+spt+plotproperties['ens'][0]+'_'+it+'_'+plotproperties['time'][ipar]+plotproperties['dummy'],
                                  'page_id_line':'off', 
                                  'super_page_x_length':14.,'super_page_y_length':(spd['tbindex'].shape[0]-1)*4.+2., 
})

                    print(('out:',out.args['output_name']))  
                    glist=[out]
                    if spd['type']=='obs':
                        legend_user_lines=['SPAG','RAOB','RICH']
                    else:
                        legend_user_lines=['SPAG','RAOB','RICH-tau']
                    validprof=False
                    #for ipar in range(parmax):
                    for ib in range(1,spd['tbindex'].shape[0]-1):
                        ym=toyearmonth(numpy.asarray([spd['tbindex'][ib]],dtype='int'),'days since '+startyear+'-01-01')[0]
                        lines =["Break profile estimates, "+stname+', {0:0>2}h, '.format(ipar*12)+it+', '+spt+plotproperties['ens'][0]]
                        superpagetitle=ib==spd['tbindex'].shape[0]-2

                        if True:
                            glist+=set_spagarr(spagdict,lines,ym,xlim=[-4,4],xpos=0.,ypos=4.*ib,pxl=5.,pyl=4.,
                                               legend_user_lines=legend_user_lines,superpagetitle=False) 
                            leg=glist[-1]
                            del glist[-1]

                            tlist=list()
                            try:
                                if plotproperties['tlist']:
                                    for ispag in range(spagarr.shape[0]):
                                        if sum(numpy.isnan(spagarr[ispag,ipar,:,ib]))<pmax-1:
                                            tlist.append('{0:0>6}'.format(wmonrs[ib,ispag]))
                                else:
                                    tlist=[]

                            except:
                                tlist=[]
                                pass

                            splist=[]
                            for k in range(spagarr.shape[0]):
                                if sum(spagarr[k,ipar,:,ib]==spagarr[k,ipar,:,ib]):
                                    splist.append(spagarr[k,ipar,:,ib])
                            sv=numpy.asarray(splist).flatten()
                            if sv.size>0:
                                pps=numpy.tile(ps,len(splist)).flatten()
                                glist+=add_mprofile(sv,pps,'black',tlist=tlist,pyl=7.0,thickness=1,superpage=superpagetitle)

                            idt=0
                            for d in tasks:
                                ld=stnames.index(stname)
                                irb=numpy.asarray(numpy.where(abs(spd['tbindex'][ib]-d['dindex'][ld,:])<=1))
                                if irb.size>0:
                                    if irb[0][0]>0:
                                        if sv.size==0 and idt==0:
                                            glist+=add_mprofile(d['ddata'][ld,iens,ipar,:,irb[0][0]]-d['ddata'][ld,iens,ipar,:,irb[0][0]-1],
                                                                ps,'black',tlist=tlist,pyl=7.0,thickness=1,superpage=superpagetitle)
                                            idt=1

                                        glist+=add_mprofile(d['ddata'][ld,iens,ipar,:,irb[0][0]]-d['ddata'][ld,iens,ipar,:,irb[0][0]-1],
                                                            ps,d['color'],pyl=7.0,thickness=12,superpage=superpagetitle)
                                        validprof=True
                                        glist+=[leg]

                            if validprof:
                                if superpagetitle:
                                    legend = mlegend( {"legend":"on", 
                                                       "legend_text_font_size":0.5, 
                                                       "legend_text_colour":"black",
                                                       "legend_box_mode":"positional",
                                                       "legend_box_x_position":-1.0,
                                                       "legend_box_y_position":4.0-0.4,
                                                       "legend_box_x_length":5.0,
                                                       "legend_box_y_length":1.2,
                                                       "legend_text_composition":"user_text_only",
                                                       "legend_user_lines":legend_user_lines
                                                       })
                                    glist+=[legend]

                        glist+=set_mapspagarr(spagdict,lines,ym,xlim=[-4,4],xpos=6.0,ypos=4.0*ib+0.4,pxl=7.0,pyl=4.0,
                                              legend_user_lines=legend_user_lines,superpagetitle=superpagetitle) 


                        #tlist=[stname]
                        idx=jdict['snrs'].index(stname)
                        lonlist=[jdict['lons'][idx]]
                        latlist=[jdict['lats'][idx]]
                        try:
#			    if plotproperties['tlist']:
                            for ispag in range(spagarr.shape[0]):
                                if sum(numpy.isnan(spagarr[ispag,ipar,:,ib]))<pmax-1:
                                    ts='{0:0>6}'.format(wmonrs[ib,ispag])
                                    idx=jdict['snrs'].index(ts)
                                    lonlist.append(jdict['lons'][idx])
                                    latlist.append(jdict['lats'][idx])
                        except:
                            pass

                        clev=numpy.linspace(-1.6,0.4,21)
                        clev=numpy.append(clev,10.)
                        clev=numpy.append(-10.,clev)
                        symb=msymb(symbol_advanced_table_selection_type = "list",
                                   symbol_marker_index = 27,
                                   symbol_table_mode = "advanced",
                                   legend = "off",
                                   symbol_type = "marker",
                                   symbol_outline='on',
                                   symbol_outline_colour='black',
                                   symbol_advanced_table_interval = 0.2,
                                   symbol_advanced_table_level_list = clev,
                                   symbol_advanced_table_reference_level= -1.6,
                                   symbol_advanced_table_level_count=21,
                                   symbol_advanced_table_level_tolerance=3,
                                   symbol_advanced_table_height_max_value = 0.6,
                                   symbol_advanced_table_height_min_value = 0.6,
                                   symbol_advanced_table_height_method = 'calculate',
                                   symbol_advanced_table_colour_direction = "clockwise",
                                   symbol_advanced_table_min_level_colour = "reddish_purple",
                                   symbol_advanced_table_max_level_colour = "red_orange")
                        vals=numpy.zeros(len(latlist))
                        vals[0]=1.0
                        inp=minput({'Input_values':vals,  #.ravel?
                                    'Input_y_values':numpy.asarray(latlist,dtype=numpy.float64),
                                    'Input_x_values':numpy.asarray(lonlist,dtype=numpy.float64),
                                    'Input_Type':'geographical'})

                        coastlines = mcoast(map_coastline_resolution='low')
                        if len(latlist)>0:
                            glist=glist+[inp,symb,coastlines]

                    t=time.time()
                    plot(glist)
                            #try:
                                #del context
                                #del actions
                            #except:
                                #pass
                    print((time.time()-t))#,sv.size

    return


# plot sonde locations for breaks
def spagmap_allbreaks(path,tasks,stnames,plotproperties,spagdict,jdict,spt,startyear):

    global context
    global actions
    ps=plotproperties['ps']
    if 'dummy' not in list(plotproperties.keys()):
        plotproperties['dummy']=''
    #print spagdict.keys(),plotproperties['iteration']

    for stname in list(spagdict.keys()):
        mktmppath(plotproperties['tmppath']+'/'+stname)
        for it in plotproperties['iteration']:
            if len(spagdict[stname][it])==0:
                continue
            #print spagdict[stname][it][0].keys()
            for iens in range(len(spagdict[stname][it])):
                spd=spagdict[stname][it][iens]
                #print 'spd',spd
                l,pmax,parmax,tbi=spd['l'],spd['pmax'],spd['parmax'],spd['tbi']
                spagarr=numpy.empty([l*pmax*parmax*tbi],dtype='float32')
                spagarr.fill(numpy.nan)
                spagarr[spd['spagindex']]=spd['spagarrc']
                spagarr=numpy.reshape(spagarr,[l,parmax,pmax,tbi])
                wmonrs=spd['wmonrs']
                for ipar in range(len(plotproperties['time'])):
                #for ib in range(1,spd['tbindex'].shape[0]-1):
                    out = output({"output_name_first_page_number":"off",
                                  "output_formats":plotproperties["outputformat"], 
                                  'output_name': plotproperties['tmppath']+'/'+stname+'/spagmap_'+stname+'_'+spt+plotproperties['ens'][0]+'_'+it+'_'+plotproperties['time'][ipar]+plotproperties['dummy'],
                                  'page_id_line':'off', 
                                  'super_page_x_length':14.,'super_page_y_length':float((spd['tbindex'].shape[0]-1)/2)*4.0+1.2, 
})

                    print(('out:',out.args['output_name']))  
                    glist=[out]
                    if spd['type']=='obs':
                        legend_user_lines=['SPAG','RAOB','RICH']
                    else:
                        legend_user_lines=['SPAG','RAOB','RICH-tau']
                    validprof=False
                    #for ipar in range(parmax):
                    for ib in range(1,spd['tbindex'].shape[0]-1):
                        ym=toyearmonth(numpy.asarray([spd['tbindex'][ib]],dtype='int'),'days since '+startyear+'-01-01')[0]
                        lines =["Neighbor locations, "+stname+', {0:0>2}h, {1}, '.format(ipar*12,it)+spt+plotproperties['ens'][0]]
                        superpagetitle=ib==spd['tbindex'].shape[0]-2

                        glist+=set_mapspagarr(spagdict,lines,ym,xlim=[-4,4],xpos=7.*numpy.mod(ib-1,2),ypos=4.0*((ib-1)/2),pxl=7.0,pyl=4.0,
                                              legend_user_lines=legend_user_lines,superpagetitle=superpagetitle) 
                        leg=glist[-1]
#			del glist[-1]
#			for i in range(spagarr.shape[0]):

                        tlist=[stname]
                        idx=jdict['snrs'].index(stname)
                        lonlist=[jdict['lons'][idx]]
                        latlist=[jdict['lats'][idx]]
                        try:
#			    if plotproperties['tlist']:
                            for ispag in range(spagarr.shape[0]):
                                if sum(numpy.isnan(spagarr[ispag,ipar,:,ib]))<pmax-1:
                                    tlist.append('{0:0>6}'.format(wmonrs[ib,ispag]))
                                    idx=jdict['snrs'].index(tlist[-1])
                                    lonlist.append(jdict['lons'][idx])
                                    latlist.append(jdict['lats'][idx])
                        except:
                            pass

                        clev=numpy.linspace(-1.6,0.4,21)
                        clev=numpy.append(clev,10.)
                        clev=numpy.append(-10.,clev)
                        symb=msymb(symbol_advanced_table_selection_type = "list",
                                   symbol_marker_index = 27,
                                   symbol_table_mode = "advanced",
                                   legend = "off",
                                   symbol_type = "marker",
                                   symbol_outline='on',
                                   symbol_outline_colour='black',
                                   symbol_advanced_table_interval = 0.2,
                                   symbol_advanced_table_level_list = clev,
                                   symbol_advanced_table_reference_level= -1.6,
                                   symbol_advanced_table_level_count=21,
                                   symbol_advanced_table_level_tolerance=3,
                                   symbol_advanced_table_height_max_value = 0.6,
                                   symbol_advanced_table_height_min_value = 0.6,
                                   symbol_advanced_table_height_method = 'calculate',
                                   symbol_advanced_table_colour_direction = "clockwise",
                                   symbol_advanced_table_min_level_colour = "reddish_purple",
                                   symbol_advanced_table_max_level_colour = "red_orange")
                        vals=numpy.zeros(len(latlist))
                        vals[0]=1.0
                        inp=minput({'Input_values':vals,  #.ravel?
                                    'Input_x_missing_value':-1.e21,
                                    'Input_y_missing_value':-1.e21,
                                    'Input_missing_value':-1.e21,
                                    'Input_latitude_values':numpy.asarray(latlist,dtype=numpy.float64),
                                    'Input_longitude_values':numpy.asarray(lonlist,dtype=numpy.float64),
                                    'Input_Type':'geographical'})

                        coastlines = mcoast(map_coastline_resolution='low')
                        if len(latlist)>0:
                            glist=glist+[inp,symb,coastlines]
                        splist=[]
                        for k in range(spagarr.shape[0]):
                            if sum(spagarr[k,ipar,:,ib]==spagarr[k,ipar,:,ib]):
                                splist.append(spagarr[k,ipar,:,ib])
                        sv=numpy.asarray(splist).flatten()
                        if sv.size>0:
                        #sv=spagarr[:,ipar,:,ib].flatten()
                            pps=numpy.tile(ps,len(splist)).flatten()
#			    glist+=add_mprofile(sv,pps,'black',tlist=tlist,pyl=7.0,thickness=1,superpage=superpagetitle)


                    t=time.time()
                    plot(glist)
                    print((time.time()-t,sv.size))

    return

def set_spagarr(allspag,lines,ym,xlim=[-4,4],xpos=1.5,ypos=1.5,pxl=13.,pyl=10.,legend_user_lines=[],superpagetitle=False):

    ppage = page(
        layout='positional',  
        page_x_length=pxl, 
        page_y_length=pyl, 
        page_id_line='off',
        page_x_position=1.5+xpos, 
        page_y_position=1.0+ypos
    )
    projection = mmap(
        subpage_map_projection='cartesian',
        subpage_x_axis_type='regular',
        subpage_y_axis_type='logarithmic',
        subpage_x_min=numpy.float(xlim[0]),
        subpage_x_max=numpy.float(xlim[1]),
        subpage_y_min=1000.,
        subpage_y_max=20.,
        subpage_y_length=pyl-1.2,
        subpage_x_length=pxl-1.,
        subpage_y_position=1.0 ,
    )
    # Vertical axis
    vertical = maxis(
        axis_orientation='vertical',
        axis_grid='off',
        axis_type='logarithmic',
        axis_tick_label_height=0.5,
        axis_tick_label_colour='black',
        axis_grid_colour='black',
        axis_grid_thickness=1,
        axis_grid_reference_line_style='solid',
        axis_grid_reference_thickness=1,
        axis_grid_line_style='dash',
        axis_title='on',
        axis_title_text='Pressure',
        axis_title_height=0.5,
    )
    # Horizontal axis
    horizontal = maxis(
        axis_orientation='horizontal',
        axis_type='regular',
        axis_tick_label_height=0.5,
        axis_tick_label_colour='black',
        axis_grid='on',
        axis_grid_colour='black',
        axis_grid_thickness=1,
        axis_grid_line_style='dash',
        axis_title='on',
        axis_title_text='K',
        axis_title_height=0.5,
    )

    plist=[ppage,projection,horizontal,vertical]
    if superpagetitle:

        title = mtext(
            text_lines= lines,
            text_html= 'true',
            text_colour= 'black',
            text_font_size= 0.6,
            text_mode = 'positional',
            text_box_x_position= 0.,
            text_box_y_position= pyl+0.1,
            text_box_x_length= pxl,
            text_box_y_length= 0.6,
            text_border= "off",
            text_justification = "left"    )

        plist+=[title]

    tym = mtext(
        text_lines= ['{0}'.format(ym)],
        text_html= 'true',
        text_colour= 'black',
        text_font_size= 0.6,
        text_mode = 'positional',
        text_box_x_position= 0.5,
        text_box_y_position= pyl-0.7,
        text_box_x_length= pxl,
        text_box_y_length= 0.6,
        text_border= "off",
        text_justification = "left"    )

    plist+=[tym]


    return plist

def set_mapspagarr(allspag,lines,ym,xlim=[-4,4],xpos=1.5,ypos=1.5,pxl=13.,pyl=10.,legend_user_lines=[],superpagetitle=False):

    ppage = page(
        layout='positional',  
        page_x_length=pxl, 
        page_y_length=pyl, 
        page_id_line='off',
        page_x_position=0.2+xpos, 
        page_y_position=0.2+ypos
    )
    projection =mmap({"subpage_map_projection":"cylindrical", 
                      "subpage_y_length":pyl-0.2,
                      "subpage_x_length":pxl-0.4,
                      "subpage_y_position":1.0 ,
                      "subpage_lower_left_longitude": -180.,
                      "subpage_lower_left_latitude":-90.,
                      "subpage_upper_right_longitude":180.,
                      "subpage_upper_right_latitude":90.,
                      "page_id_line":"off",
                      "map_grid":"off",
                      "map_label":"off",
                      "map_grid_latitude_increment":30.,
                      "map_grid_longitude_increment":60.,
                      "map_label_height":0.4})


    plist=[ppage,projection]
    if superpagetitle:

        title = mtext(
            text_lines= lines,
            text_html= 'true',
            text_colour= 'black',
            text_font_size= 0.6,
            text_mode = 'positional',
            text_box_x_position= 0.,
            text_box_y_position= pyl+0.5,
            text_box_x_length= pxl,
            text_box_y_length= 0.6,
            text_border= "off",
            text_justification = "left"    )

        plist+=[title]

    tym = mtext(
        text_lines= ['{0}'.format(ym)],
        text_html= 'true',
        text_colour= 'black',
        text_font_size= 0.6,
        text_mode = 'positional',
        text_box_x_position= 0.1,
        text_box_y_position= pyl-0.05,
        text_box_x_length= pxl,
        text_box_y_length= 0.6,
        text_border= "off",
        text_justification = "left"    )

    plist+=[tym]


    return plist

def add_mprofile(prof,ps,colour,graphdict={},pyl=13.,tlist=[],thickness=1,legend_user_lines='',superpage=False):

    mask=numpy.isnan(prof)
    prof[mask]=-1.e21
    data = minput(Input_x_values=prof,
                  Input_y_values=ps,
                  Input_x_missing_value=-1.e21
                  )
    lg="off"
    if superpage:
        lg="on"
    graph = mgraph( graph_line_colour= colour, #task["color"], 
                    graph_line_thickness= thickness,
                    graph_x_suppress_above=10.,
                    graph_x_suppress_below=-10.,
                    legend=lg,*graphdict
                    )

    plist=[data,graph]

    tlist=[]
    if tlist:

        if len(tlist)>20:
            dtlist=[]
            tstr=''
            for t in tlist:
                if tstr=='':
                    tstr+=t
                else:
                    tstr+=', '+t
                    dtlist.append(tstr)
                    tstr=''
            if tstr!='':
                dtlist.append(tstr)
            txp=3.5
        else:
            dtlist=tlist
            txp=4

        tl = mtext(
            text_lines= dtlist,
            text_line_height_ratios= [1.0]*len(dtlist),
            text_html= 'true',
            text_colour= 'black',
            text_font_size= 0.3,
            text_mode = 'positional',
            text_box_x_position= txp,
            text_box_y_position= 1.0,
            text_box_x_length= 2.,
            text_box_y_length= 0.02*len(dtlist),
            text_border= "off")

        plist+=[tl]

    return plist

def add_spagprofiles(task,ps,ibelt):


    esiz=task["beltslopes"].shape[0]
    estart=numpy.round(esiz*0.025)
    estop=numpy.round(esiz*0.975)
    if estop==esiz:
        estop=esiz-1
    ens05=numpy.zeros(ps.shape[0])
    ens95=numpy.zeros(ps.shape[0])

    for ip in range(ps.shape[0]):
        enstrends=task["beltslopes"][:,ibelt,ip].flatten()
        index=numpy.argsort(enstrends)
        ens05[ip]=enstrends[index[estart]]
        ens95[ip]=enstrends[index[estop]]
    if(numpy.nanstd(ens05-ens95)<1.e-7):
        ens95=ens05+0.01
    mask=numpy.isnan(ens05)
    if sum(mask)!=0:
        ens05[mask]=-1.e21
        ens95[mask]=-1.e21
    data = minput(Input_x_values=ens05,
                  Input_x2_values=ens95,
                  Input_y_values=ps,
                  Input_y2_values=ps,
                  Input_x_missing_value=-1.e21
                  )

    graph = mgraph( graph_shade='on',
                    graph_shade_colour= task["color"], 
                    graph_line_colour= task["color"], 
                    graph_shade_style= "area_fill",
                    graph_type='area',
                    graph_x_suppress_above=10.,
                    graph_x_suppress_below=-10.,
                    legend="off",
                    )
    return [data,graph]

#def addhadtrend(medtrend,enstrends):


    #p=numpy.asarray([990.])
    #data = minput(Input_x_values=medtrend.flatten(),
                    #Input_y_values=p,
                    #)

    #symb = msymb(   symbol_marker_index= 28,
                    #symbol_colour='black',
                    #symbol_height= 0.5,
                    #symbol_type="marker",
                    #legend="on",
                    #legend_user_text='HadCRUT4'
                    #)
    #estart=2
    #estop=97
    #index=numpy.argsort(enstrends.flatten())
    #ens05=numpy.asarray([enstrends[index[estart]],enstrends[index[estart]],enstrends[index[estop]],
                            #enstrends[index[estop]],enstrends[index[estop]],enstrends[index[estop]]])

    #mask=numpy.isnan(ens05)
    #if sum(mask)!=0:
        #ens05[mask]=-1.e21
    #dataens=minput(Input_x_values=ens05,
                    #Input_y_values=numpy.asarray([970.,1010.,990.,990.,970.,1010.]),
                    #Input_x_missing_value=-1.e21
                    #)
    #graph = mgraph( graph_line_colour= 'black', 
                    #graph_line_style = "solid",
                    #graph_line_thickness = 2,
                    #legend="off",
                    #)
    #return [data,symb,dataens,graph]



