import os,glob,sys
import numpy
from rasotools.additions.anomaly import *
#from Magics.macro import *
import netCDF4
from rasotools.additions.utils import *
from rasotools.additions.dailyseries import *
from rasotools.additions.statanalysis import *
from rasotools.additions.set_trendmaps import *
import time
import copy
from collections import OrderedDict
#import matplotlib.pyplot as plt
#from matplotlib.patches import Polygon

base = datetime.date(1900,1,1)
date_list = [base + datetime.timedelta(days=x) for x in range(0, 45000)]
date_listm = numpy.asarray([date_list[i].month for i in range(0, 45000)])
date_listy= [(datetime.date(x,1,1)-base).days for x in range(1900,2020)]

def set_bias(lines,legend_user_lines,plotproperties,obprop,sxpos=0.0,pkey='beltslopes'):

    beltinterval=obprop['bi']
    dsx=obprop["dxpos"]

    ptype='logarithmic'
    try:
        if plotproperties['logp']=='False':
            ptype='regular'
    except:
        pass

    apl=[-2.0,-1.0,0.0]
    haxistype='regular'
    haxpl=[0,1]
    if pkey!='cost':
        try:
            haxpl=plotproperties['haxpl']
            haxistype=plotproperties['haxistype']
        except:
            pass

    dsxx=dsx
    if sxpos==0.0:
        dsxx+=1.5
    ppage=page(
        page_y_length=12.5,
        page_y_position=0. ,
        page_x_position=sxpos ,
        page_x_length=dsxx,
        page_frame='off',
        layout='positional',
        page_id_line='off'
)
    projection = mmap(
        subpage_map_projection='cartesian',
        subpage_x_axis_type='regular',
        subpage_y_axis_type=ptype,
        subpage_x_min=beltinterval[0],
        subpage_x_max=beltinterval[1],
        subpage_y_min=1000.,
        subpage_y_max=20.,
        subpage_y_length=10.,
        subpage_y_position=1.5 ,
        subpage_x_position=dsxx-dsx,
        subpage_x_length=dsx,
    )
    # Horizontal axis
#    axis_tick_position_list=numpy.asarray(haxpl,dtype='float32'),
    hxd=dict(
        axis_orientation='horizontal',
        axis_type=haxistype,
        axis_tick_label_height=0.4,
        axis_tick_label_colour='black',
        axis_grid='on',
        axis_grid_colour='black',
        axis_grid_thickness=1,
        axis_grid_line_style='dot',
        axis_title="on",
        axis_grid_reference_level=0.,
        axis_title_text=obprop["axtitle"])

    if haxistype=='position_list':
        hxl=[];
        for h in haxpl:
            hxl.append('{:4.1f}'.format(h))
        hxd['axis_tick_label_list']=hxl

    horizontal = maxis(hxd)

    plist=[ppage,projection,horizontal]
    
    if sxpos==0.0:
        # Vertical axis
        vertical = maxis(
            axis_orientation='vertical',
            axis_grid='on',
            axis_type=ptype,
            axis_tick_label_height=0.4,
            axis_tick_label_colour='black',
            axis_grid_colour='black',
            axis_grid_thickness=1,
            #axis_grid_reference_line_style='dot',
            #axis_grid_reference_thickness=1,
            axis_grid_line_style='dot',
            axis_title='on',
            axis_title_text='Pressure',
            axis_title_height=0.6,
        )
    
        plist+=[vertical]
    else:
        vertical=maxis(
            axis_orientation='vertical',
            axis_grid='on',
            axis_grid_colour='black',
            axis_grid_thickness=1,
            #axis_grid_reference_line_style='dot',
            #axis_grid_reference_thickness=1,
            axis_grid_line_style='dot',
            axis_type=ptype,
            axis_tick_label='off',
            axis_title='off',
            axis_tick='off')
        legend = mlegend( {"legend": "off"})
        plist+=[vertical]
       
       
       
    #sxpos+=dsx
    return plist

def addxprofile(task,ps,iens,ibelt,obprop,pkey="zslopes"):


    hilf=task[pkey][iens,ibelt+obprop["ioff"],:ps.shape[0]].flatten()
#    if sum(hilf==hilf)==0:
#        return []
    mask=numpy.isnan(hilf)
    hilf[mask]=-1.e21
    data = minput(Input_x_values=hilf,
                  Input_y_values=ps,
                  Input_x_missing_value=-1.e21
                  )
    if pkey=='cost':
        supp=[0.,10000.]
    else:
        supp=[-10.,10.]

    graph = mgraph( graph_line_colour= task["color"], 
                    graph_line_thickness= 12,
#                    graph_x_suppress_above=supp[1],
#                    graph_x_suppress_below=supp[0],
                    legend=obprop["legend"],
                    )
    return [data,graph]

def add_xens_profile(task,ps,ibelt,pkey='beltslopes'):


    esiz=task[pkey].shape[0]
    estart=numpy.round(esiz*0.025)
    estop=numpy.round(esiz*0.975)
    if estop==esiz:
        estop=esiz-1
    ens05=numpy.zeros(ps.shape[0])
    ens95=numpy.zeros(ps.shape[0])

    for ip in range(ps.shape[0]):
        enstrends=task[pkey][:,ibelt,ip].flatten()
        index=numpy.argsort(enstrends)
        ens05[ip]=enstrends[index[estart]]
        ens95[ip]=enstrends[index[estop]]
#    if sum(ens05==ens05)==0:
#        return []
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

    if pkey=='cost':
        supp=[0.,10000.]
    else:
        supp=[-10.,10.]
    graph = mgraph( graph_shade='on',
                    graph_shade_colour= task["color"], 
                    graph_line_colour= task["color"], 
                    graph_shade_style= "area_fill",
                    graph_type='area',
#                    graph_x_suppress_above=supp[1],
#                    graph_x_suppress_below=supp[0],
                    legend="on",
                    )
    return [data,graph]

def do_obstat(ddata,obsdata,plotproperties,startdate=1900):
    out=numpy.zeros((18,ddata.shape[1]))
    first=date_listy[plotproperties['diffintervals'][0][0]-startdate]
    last=date_listy[plotproperties['diffintervals'][0][1]+1-startdate]


    xdata=ddata[:,0,first:last]
    for ip in range(out.shape[1]):
        xdata[:]=ddata[:,ip,first:last]
        for ipar in range(2):
            out[ipar,ip]=numpy.nanmean(xdata[ipar:])
        out[2,ip]=numpy.nanmean(xdata)

        for ipar in range(2):
            out[9+ipar,ip]=numpy.nanstd(xdata[ipar,:])
        out[11,ip]=numpy.nanstd(xdata)
        
        for ipar in range(2):
            out[3+ipar,ip]=numpy.sqrt(out[9+ipar,ip]**2+out[ipar,ip]**2)
        out[5,ip]=numpy.sqrt(out[11,ip]**2+out[2,ip]**2)
        
        odata=obsdata[0,:,ip,first:last]
        for ipar in range(2):
            mask=~numpy.isnan(xdata[ipar,:])*~numpy.isnan(odata[ipar,:])
            try:
                out[12+ipar,ip]=numpy.corrcoef(xdata[ipar,mask]+odata[ipar,mask],odata[ipar,mask])[0][1]
            except:
                out[12+ipar,ip]=numpy.nan
        out[14,ip]=numpy.nanmax(out[12:14,ip])
            
            
        for ipar in range(2):
            out[6+ipar,ip]=numpy.sum(~numpy.isnan(xdata[ipar,:]))/100.
        out[8,ip]=numpy.sum(~numpy.isnan(xdata))/100.
        
    return out
@njit(cache=True)
def do_obstrends(out,ddata,obsdata,iens,first,last):
    #out=numpy.empty((6,ddata.shape[1]))
    #out.fill(numpy.nan)

    xdata=ddata[iens,:,0,first:last]+obsdata[0,:,0,first:last]
    xtime=numpy.arange(last-first)/365.25
    d24=numpy.empty(last-first)
    for ip in range(ddata.shape[2]):
        xdata[:]=ddata[iens,:,ip,first:last]+obsdata[0,:,ip,first:last]
        d24[:]=0.
        for i in range(last-first):
            l=0
            if xdata[0,i]==xdata[0,i]:
                d24[i]+=xdata[0,i]
                l+=1
            if xdata[1,i]==xdata[1,i]:
                d24[i]+=xdata[1,i]
                l+=1
            if l>0:
                d24[i]/=l
            else:
                d24[i]=numpy.nan
        for ipar in range(2):
            l=0
            for i in range(last-first):
                if xdata[ipar,i]==xdata[ipar,i]:
                    l+=1
            if l>(last-first)*0.8:
                out[iens,ipar,ip]=fastlinregress(xtime, xdata[ipar,:])*10.
            else:
                out[iens,ipar,ip]=numpy.nan
        out[iens,2,ip]=numpy.nan
        l=0
        for i in range(last-first):
            if d24[i]==d24[i]:
                l+=1
        if l>(last-first)*0.8:
            out[iens,2,ip]=fastlinregress(xtime, d24)*10.
    
    return

def obstat(path,tasks,plotproperties,nml,stnames,minlen=360,ref=''):

    tt=time.time()
    print(('I am in '+path+ ' should be in ',plotproperties["exps"]))
    pindex=plotproperties["pindex"]
    try:
        dummy=plotproperties['dummy']
    except:
        dummy=''

    istat,lats,lons,ps,goodsts,stlongnames=read_dailyseries(path,tasks,nml,stnames,pindex,minlen=minlen)
    if istat==0:
        return
    diff=''
    if ref != '':
        tasks2=copy.deepcopy(tasks)
        istat2,lats2,lons2,ps2,goodsts2,stlongnames2=read_dailyseries(path,tasks2,nml,stnames,pindex,minlen=minlen,ref=ref)
        diff='_diff'

        for j in range(len(tasks)):
            tasks[j]['ddata']-=tasks2[j]['ddata']

    if isinstance(goodsts,list):
        goodsts=numpy.asarray(goodsts)
    istnames=goodsts.astype(numpy.int32)
    tmshape=tasks[0]["ddata"].shape

    data=OrderedDict()
    if 'trends' in plotproperties['plotlist']:
        
        dint=plotproperties['intervals'][0:2]
    else:
        dint=plotproperties['diffintervals'][0:2]
       
    for d in tasks:
        ens=0
        # replace RAOBCORE adjustments with unadjusted series
        t=time.time()
        if str(d["shortname"]) in ("tm",'obs',"an20cr",'rad'):
            data['current']=dailyanomalydriver(d,ens,dint)
            if d["shortname"] in ('rad'):
                data['current'][0,0,:,:]=data['current'][0,1,:,:]-data['current'][0,0,:,:]
                data['current'][0,1,:,:]=numpy.nan
            else:
                data[d["shortname"]]=data['current'].copy()                
        elif d["shortname"] == "tmcorr":
            if ens==0:
                expandandadd(d["ddata"],data['current'],d["dindex"],
                             numpy.arange(data['current'].shape[2],dtype=numpy.int64),ens,data['current'],-1.0)
        elif d["shortname"] in ("an20c","andep","bgdep","era5v2","eraibc","bgcorrection","jra55_andep","jra55_fgdep","eijra_fgdep"):  
            data["current"]=d["ddata"][:,ens,:,:,:]
            if d["shortname"] in ["bgdep",'andep']:
                data["current"]=-data["current"]
            data[d["shortname"]]=data['current'].copy()
        elif d["shortname"] in ("e20c_andep","ce20c0_andep","erapresat_andep","n20c_andep","ce20cens_mean","ce20cens_spread"):  
            data["current"]=d["ddata"][:,ens,:,:,:]
            data[d["shortname"]]=data['current'].copy()
        elif d["shortname"] == "jra55":  
            data["current"]=data['current']+data['jra55_andep']
        elif d["shortname"] in ("20cr-e20c","presat-e20c"):
            data["current"]=dailyanomalydriver(d,ens,dint)-data['e20c']              
        elif d["shortname"] == "aninc":
            data["current"]=d["ddata"][:,ens,:,:,:]-data['bgdep']
        elif d["shortname"] in ("e20c-ce20c"):
            data["current"]=d["ddata"][:,ens,:,:,:]-data['ce20c']
        elif d["shortname"] in ["rcorr"] or 'rio' in d["shortname"][0:3] or 'rit' in d["shortname"][0:3]:
            for ens in d["ens"]:
                if ens==0:
                    dds=d["ddata"].shape
                    nullref=numpy.zeros((dds[0],dds[2],dds[3],nml['rfpar']['nmax']))
                    data['current']=copy.copy(nullref)
                    data['accum']=copy.copy(nullref)
                expandandadd(d["ddata"],data['current'],d["dindex"],numpy.arange(nullref.shape[2],dtype=numpy.int64),
                             ens,data['accum'],1.0)
                data['current'][:]=data['accum'][:]
            data[d["shortname"]]=-data['accum'].copy()                
            try:
                data['current']/=-numpy.array(d["ens"]).size
                data[d["shortname"]]/=numpy.array(d["ens"]).size
            except:
                pass
            print((d['shortname'],d["ddata"][0,0,0,d['ddata'].shape[3]-1,:]))

        elif 'bgdepr' in d['shortname']:
            for bgc in ["bgdeprio","bgdeprit","bgdeprcorr"]: 
                if bgc in d['shortname']:#,"eijradeprcorr","eijradepriocorr","eijradepritcorr"]:
                    expandandadd(d["ddata"],data['bgdep'],d["dindex"],
                                 numpy.arange(data['bgdep'].shape[2],
                                              dtype=numpy.int64),ens,data['current'],-1.0)
        else:
            print((d["shortname"],' not implemented'))


#            print 'expand',ens,time.time()-t
#            plotlist=list()
        if istnames.shape[0]==1 or not plotproperties["super"]:
            pass
        else:
            data['currentav']=numpy.zeros(data['current'].shape[1:4])
            avcount=data['currentav'].copy()
            d['minindex'][:]=0
            d['maxindex'][:]=data['current'].shape[3]-1
            stationaverage(data['current'],data['currentav'],avcount,d['minindex'],d['maxindex'],5)
            print(('saverage,',time.time()-tt))
            data[d['shortname']][0,:,:,:]=data['currentav']
    
    interval=numpy.asarray(plotproperties['intervals'],dtype=numpy.int)
    null=numpy.zeros(data['current'].shape)
    for d in tasks:    
        first=date_listy[interval[0]-d['startdate']//10000]
        last=date_listy[interval[1]+1-d['startdate']//10000]
        for iens in [0]:
            if 'trends' in plotproperties['plotlist']:
                if d['shortname']=='obs':
                    do_obstrends(d["beltslopes"],null,data['obs'],iens,first,last)
                else:
                    do_obstrends(d["beltslopes"],data[d['shortname']],data['obs'],iens,first,last)
                   
            else:
                d["zslopes"][iens,:,0:pindex.shape[0]]=do_obstat(data[d['shortname']][iens,:,:],data['obs'],plotproperties,d['startdate']//10000)
                    


    rmax=6.0/1.02
    bmax=4.0/1.02
    tmax=2.0/1.02
    for d in tasks:
        rmax=max([rmax,numpy.nanmax(d["zslopes"][0,5,:pindex.shape[0]])])
        bmax=max([bmax,numpy.nanmax(numpy.abs(d["zslopes"][0,2,:pindex.shape[0]]))])
        tmax=max([tmax,numpy.nanmax(numpy.abs(d["beltslopes"][0,2,:pindex.shape[0]]))])
        d['zslopes'][numpy.isnan(d['zslopes'])]=-1.e21
    rmax*=1.02
    bmax*=1.02
    tmax*=1.02
    print(('max:',rmax,bmax,tmax,time.time()-tt))
    
    obprop=dict(trends={"bi":[-tmax,tmax],"dxpos":6.0,"legend":"off","axtitle":"bias [K]","ioff":0},
                bias={"bi":[-bmax,bmax],"dxpos":6.0,"legend":"off","axtitle":"bias [K]","ioff":0},
                corr={"bi":[-0.2,1.],"dxpos":6.0,"legend":"off","axtitle":"corr []","ioff":12},
                rms={"bi":[0.0,rmax],"dxpos":6.0,"legend":"off","axtitle":"rms [K]","ioff":3},
                std={"bi":[0.0,rmax],"dxpos":6.0,"legend":"off","axtitle":"std [K]","ioff":9},
                count={"bi":[0.,numpy.max(tasks[0]["zslopes"][0,8,:pindex.shape[0]])],"dxpos":3.0,"legend":"off","axtitle":"count [100]","ioff":6})

    for k in list(obprop.keys()):
        if k==plotproperties['plotlist'][-1]:
            obprop[k]['legend']='on'
    ti=list(range(len(tasks)))
    li=[ta['shortname'] for ta in tasks]
    if 'rio' in li:
        ti.insert(0,ti.pop(li.index('rio')))
    if 'rit' in li:
        ti.insert(0,ti.pop(li.index('rit')))
    if 'ce20c_andep' in li:
        ti.insert(0,ti.pop(li.index('ce20c_andep')))
    
    mktmppath(plotproperties['tmppath']+'/'+goodsts[0])
    
    if 'trends' in plotproperties['plotlist']:
        plot_obstrends(goodsts,stnames,lats,lons,dummy,tasks,plotproperties,obprop,ti,ps,pindex)
    else:
        plot_obstat(goodsts,stnames,lats,lons,dummy,tasks,plotproperties,obprop,ti,ps,pindex)

#    plot_obstrends()
    
    print((time.time()-tt))

    return
    
def plot_obstrends(goodsts,stnames,lats,lons,dummy,tasks,plotproperties,obprop,ti,ps,pindex):
    spage = page(
        layout='positional',  
        page_x_length=29.7, 
        page_y_length=21., 
        page_id_line='off',
        page_x_position=0., 
        page_y_position=0.)
    interval=numpy.asarray(plotproperties['intervals'],dtype='int')
    for ipar in range(3):
        if '{:0>2}'.format(ipar*12) not in plotproperties['time']:
            continue
        linesplustext=[]
        legend_user_lines=[]
        lines =["Temperature Trends [K/10a], {:0>2}h".format(ipar*12)+', '+
                "{:4}".format(interval[0])+"-{:4}".format(interval[1])+', '+goodsts[0],
                str('Version '+plotproperties['version']+', '+plotproperties['exps'][0])]

        profsat=''
        if 'satseries' in plotproperties['plotlist']:
            profsat='sat'

        out = output({"output_name_first_page_number":"off",
                      "output_formats":plotproperties["outputformat"], 
                      'output_name': plotproperties['tmppath']+'/'+goodsts[0]+'/'+profsat+'trends_'+goodsts[0]+
                      "_{:4}".format(interval[0])+"-{:4}".format(interval[1])+'_{:0>2}'.format(ipar*12)+dummy})

        plotlist=[spage,out]

        for key in ti:
            shade=1
            legname=copy.copy(tasks[key]["shortname"])
            if 'andep' in legname:
                legname=legname.split('_andep')[0]
            if legname not in [ "rio","rit", 'ce20c' ] or shade==0:
                for iens in [0]:#tasks[key]["ens"]: 
 
                    if 'satseries' in plotproperties['plotlist']:
                        if legname != "uah" and legname != "rss":
                            tasks[key]["zlopes"][iens,ibelt,-4]=numpy.nan
                        linesplustext=linesplustext+addsattrend(tasks[key]["zslopes"][iens,ibelt,-4:],msups,
                                                                scol=tasks[key]["color"],iens=iens) # SAT
                    else:
                        if tasks[key]["zslopes"].shape[2]>=pindex.shape[0]:
                            linesplustext=linesplustext+addprofile(tasks[key],ps[pindex],iens,ipar) # RASO
                        else:
                            linesplustext=linesplustext+addsattrend(tasks[key]["beltslopes"][iens,ibelt,:],[800.,550.,250.,90.],
                                                                    scol=tasks[key]["color"],iens=iens) # SAT
 
                        legend_user_lines=legend_user_lines+[legname]
            else:
                if 'satseries' in plotproperties['plotlist']:
                    for iens in tasks[key]["ens"]: #range(len(tasks[key]["ens"])):
                        tasks[key]["beltslopes"][iens,ibelt,-4]=numpy.nan
                        linesplustext=linesplustext+addsattrend(tasks[key]["beltslopes"][iens,ibelt,-4:],[800.,550.,250.,90.],
                                                                scol=tasks[key]["color"],iens=iens) # SAT
                        print((tasks[key]["beltslopes"][iens,ibelt,-4:]))
                else:
                    linesplustext=linesplustext+add_ens_profile(tasks[key], ps[pindex], ibelt)
                legend_user_lines=legend_user_lines+[legname]

        projection,horizontal,vertical,title,legend=set_belttrends(lines,legend_user_lines,plotproperties,pkey='cost',
                                                                   beltinterval=obprop['trends']['bi'])
        plotlist=plotlist+[projection,horizontal,vertical]+linesplustext+[title,legend]


        try:
            plot(plotlist)
            print((out.args['output_name']+'.'+str(plotproperties["outputformat"][0])))
        except IOError:
            print(('no time series available:', stnames))
    
def plot_obstat(goodsts,stnames,lats,lons,dummy,tasks,plotproperties,obprop,ti,ps,pindex):
    
    spl=1.5
    for p in plotproperties['plotlist']:
        spl+=obprop[p]['dxpos']+0.3
    if spl<17.5:
        spl=17.5
    spage=page(super_page_y_length=14.5,
               super_page_x_length=spl+0.5)
    for ipar in range(3):
        if '{:0>2}'.format(ipar*12) not in plotproperties['time']:
            continue
        sxpos=0.0
#	if '{:0>2}'.format(ipar*12) in plotproperties['time']:

    #Setting the output
        interval=plotproperties['diffintervals'][0]
        if len(goodsts)<2:
            oname=plotproperties['tmppath']+'/'+goodsts[0]+'/obstat_'+goodsts[0]+"_{:4}".format(interval[0])+"-{:4}".format(interval[1])+'_{0:0>2}'.format(ipar*12)+dummy
        else:
            oname=plotproperties['tmppath']+'/'+goodsts[0]+'/obstat_'+goodsts[0]+'-'+max(goodsts)+"_{:4}".format(interval[0])+"-{:4}".format(interval[1])+'_{0:0>2}'.format(ipar*12)+dummy

        print(('output file: ',oname))
        poutput = output(output_name = oname, 
                             output_formats = plotproperties["outputformat"],
                             output_title= oname,
                             output_name_first_page_number = "off")

        ns='N0S'[1-int(lats[0]/abs(lats[0]))]
        ew='E0W'[1-int(lons[0]/abs(lons[0]))]
        lines =["T Departure Statistics, "+'{:0>2}h, '.format(ipar*12)+stnames[0]+', {:3.1f}{} {:3.1f}{}'.format(abs(lats[0]),ns,abs(lons[0]),ew),
                "{:4}".format(interval[0])+"-{:4}".format(interval[1])+', '+plotproperties['exps'][0]]

        plotlist=[spage,poutput]


        legend_user_lines=[]
        for pl in plotproperties['plotlist']:
            linesplustext=[]
            for key in ti:
                if tasks[key]['shortname']=='obs':
                    continue
                shade=1
                legname=copy.copy(tasks[key]["shortname"])
                if 'andep' in legname:
                    legname=legname.split('_andep')[0]
                if legname not in [ "rio","rit", 'ce20c' ] or shade==0:
                    for iens in tasks[key]["ens"]: 
    
                        if 'satseries' in plotproperties['plotlist']:
                            if legname != "uah" and legname != "rss":
                                tasks[key]["zlopes"][iens,ibelt,-4]=numpy.nan
                            linesplustext=linesplustext+addsattrend(tasks[key]["zslopes"][iens,ibelt,-4:],msups,
                                                                    scol=tasks[key]["color"],iens=iens) # SAT
                        else:
                            if tasks[key]["zslopes"].shape[2]>=pindex.shape[0]:
                                linesplustext=linesplustext+addxprofile(tasks[key],ps[pindex],iens,ipar,obprop[pl]) # RASO
                            else:
                                linesplustext=linesplustext+addsattrend(tasks[key]["beltslopes"][iens,ibelt,:],[800.,550.,250.,90.],
                                                                        scol=tasks[key]["color"],iens=iens) # SAT
    
                            legend_user_lines=legend_user_lines+[legname]
                else:
                    if 'satseries' in plotproperties['plotlist']:
                        for iens in tasks[key]["ens"]: #range(len(tasks[key]["ens"])):
                            tasks[key]["beltslopes"][iens,ibelt,-4]=numpy.nan
                            linesplustext=linesplustext+addsattrend(tasks[key]["beltslopes"][iens,ibelt,-4:],[800.,550.,250.,90.],
                                                                    scol=tasks[key]["color"],iens=iens) # SAT
                            print((tasks[key]["beltslopes"][iens,ibelt,-4:]))
                    else:
                        linesplustext=linesplustext+add_ens_profile(tasks[key], ps[pindex], ibelt)
                    legend_user_lines=legend_user_lines+[legname]
    
 #               plotlist=plotlist+set_bias(lines,legend_user_lines,plotproperties,pl,sxpos,beltinterval=bi[pl])+linesplustext
            plotlist=plotlist+set_bias(lines,legend_user_lines,plotproperties,obprop[pl],sxpos,pkey='cost')+linesplustext
            dxi=0.0
            if sxpos==0.0:
                sxpos+=obprop[pl]["dxpos"]+1.8
                dxi=1.8
            else:
                sxpos+=obprop[pl]["dxpos"]+0.3
                dxi=0.3

        ppage=page(
            page_y_length=1.5,
            page_y_position=13.0 ,
            page_x_position=0. ,
            page_x_length=sxpos,
            page_frame='off',
            layout='positional',
            page_id_line='off',
        )

        xlegend = mlegend( {"legend": "on", 
                           "legend_border":"off",
                           "legend_text_font_size":0.5, 
                           "legend_text_colour":"black",
                           "legend_box_mode":"positional",
                           "legend_box_x_position":-sxpos+obprop[pl]["dxpos"]+1.0+dxi,
                           "legend_box_y_position":11.5,
                           "legend_box_x_length":spl-2.0,
                           "legend_box_y_length":1.5,
                           "legend_text_composition":"user_text_only",
                           "legend_column_count":4,
                           "legend_user_lines":legend_user_lines
        })
        title = mtext(
            text_lines= lines,
            text_html= 'true',
            text_colour= 'black',
            text_font_size= 0.5,
            text_mode = 'positional',
            text_box_x_position= 1.5,
            text_box_y_position= 0.,
            text_box_x_length= sxpos-2.0,
            text_box_y_length= 1.5,
            text_border= "off",
            text_justification = "left"    )
        try:
            plot(plotlist+[xlegend,ppage,title])
            print((poutput.args['output_name']+'.'+str(plotproperties["outputformat"][0])))
        except IOError:
            print(('no time series available:', stnames))



