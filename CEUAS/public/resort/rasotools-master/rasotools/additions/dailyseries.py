# importing Magics module

from Magics.macro import *
from rasotools.additions.utils import *
from rasotools.additions.anomaly import *
import time

ref = 'timeserie'
base = datetime.date(1900,1,1)
date_list = [base + datetime.timedelta(days=x) for x in range(0, 45000)]
date_listm = numpy.asarray([date_list[i].month for i in range(0, 45000)])


def dailyseries(task,currentdata,istname,lat,lon,istat,ens,ipar,ps,ip,plotproperties,longname=''):

    global base,date_list,date_listm
    # Setting the cartesian view

    if 'units' not in list(task.keys()):
        task['units']=plotproperties['units']
    pindex=plotproperties["pindex"]
    snhtparas=numpy.asarray(list(plotproperties["snht"].values()))
    hilf=currentdata[istat,ipar,ip,:].flatten()
    hilf2=hilf.copy()

#	hilf2=numpy.random.randn(hilf2.shape[0]) # for testing snht
    tmean=numpy.zeros(hilf.shape[0])
    tmean.fill(numpy.nan)
    index=numpy.zeros(currentdata.shape[3],dtype=numpy.int32)

    xtime=numpy.floor(task['startdate']//10000)+numpy.arange(currentdata.shape[3])/365.25
    tmask=numpy.logical_and(xtime>plotproperties['plotinterval'][0],xtime<plotproperties['plotinterval'][1]+1)
    slope=fastlinregress(xtime[tmask], hilf[tmask])*10.


    cs=currentdata.shape[3]
    if len(date_list)!=currentdata.shape[3]: #task['startdate']/10000!=base.year: # assume 19570101
        cstart=datetime.date(task['startdate']//10000,1,1)
        offset=(cstart-base).days
    else:
        offset=0
        #date_list = [base + datetime.timedelta(days=x) for x in range(0, currentdata.shape[3])]
        #date_listm = numpy.asarray([date_list[i].month for i in range(0, currentdata.shape[3])])

    try:
        season=False
        if 'seasonal' in plotproperties['dailytrends']:	
            season=True
            maskwinter=numpy.logical_and(tmask,numpy.logical_or(date_listm[offset:offset+cs]<3,date_listm[offset:offset+cs]>11))
            masksummer=numpy.logical_and(tmask,numpy.logical_and(date_listm[offset:offset+cs]>5,date_listm[offset:offset+cs]<9))

            slopesummer=fastlinregress(xtime[masksummer], hilf[masksummer])*10.
            slopewinter=fastlinregress(xtime[maskwinter], hilf[maskwinter])*10.
    except:
        pass
#	t=time.time()
#	print time.time()-t
    ohilf=numpy.empty(hilf.shape[0])
    ohilf[:]=hilf[:]
    omask=numpy.logical_and(xtime>plotproperties['plotinterval'][0],xtime<plotproperties['plotinterval'][1])
    if task['shortname'] in ('obs') and 'obsrmean' in list(plotproperties.keys()):
        hilf=rmean(hilf,tmean,index,numpy.int32(plotproperties["obsrmean"]))
    if task['shortname'] in ('rcorr'):
        pass
    else:
        hilf=rmean(hilf,tmean,index,numpy.int32(plotproperties["rmean"]))
#	print time.time()-t
    std=numpy.nanstd(ohilf[omask])
    n=thin(hilf,index,plotproperties["thin"])
#	print time.time()-t
    hilf[numpy.isnan(hilf)]=-1.e21
    hhilf=hilf[index[0:n]]
    xt2=xtime[index[0:n]]
    mask=numpy.logical_and(xt2>plotproperties['plotinterval'][0],xt2<plotproperties['plotinterval'][1])
    xt2=xt2[mask]	
    hhilf=hhilf[mask]
    N=numpy.sum(~numpy.isnan(ohilf[omask]))
    
    if sum(hhilf>-1.e21)>1:
        data = minput(Input_x_values=xt2,
                      Input_y_values=hhilf,
                      Input_y_missing_value=-1.e21
                      )
    else:
        return []

    ppage = page(
        layout='positional',  
        page_x_length=21., 
        page_y_length=plotproperties['dailyseriesheight'], 
        page_id_line='off',
        page_x_position=0., 
        page_y_position=(currentdata.shape[2]-1-ip)*plotproperties['dailyseriesheight']
    )

    dynrange=plotproperties["range"][:]
    print(('dynrange:',dynrange))
    hhilf=hhilf[hhilf>-1.e21]
    if(len(hhilf)>0):
        m=max(hhilf)
        print(('max',m))
        if(m*1.05>dynrange[1]):
            dynrange[1]=m*1.05
        m=min(hhilf)
        print(('min',m))
        if(m*1.05<dynrange[0]):
            dynrange[0]=m*1.05
        print(('dynrange:',dynrange))

    projection = mmap(
        subpage_map_projection='cartesian',
        subpage_x_axis_type='regular',
        subpage_y_axis_type='regular',   
        subpage_x_position=1.5,
        subpage_y_position=plotproperties['dailyseriesheight']*4./25.,
        subpage_x_length=17.5,
        subpage_y_length=plotproperties['dailyseriesheight']*17./25.,
        subpage_x_max=float(plotproperties["plotinterval"][1]),
        subpage_x_min=float(plotproperties["plotinterval"][0]),
        subpage_y_max=dynrange[1],
        subpage_y_min=dynrange[0],
    )

    # Vertical axis

    if type(task['units']) is list:
        task['units']=task['units'][0]
        
    vertical = maxis(
        axis_orientation='vertical',
        axis_grid='off',
        axis_type='regular',
        axis_tick_label_height=0.4,
        axis_grid_reference_level=0.,
        axis_grid_reference_thickness=1,
        axis_grid_reference_colour='black',
        axis_title='on',
        axis_title_text=task['units'],
        axis_title_height=0.4,
    )
    print(('units: ',task['units']))

    # Horizontal axis

    horizontal = maxis(
        axis_orientation='horizontal',
        axis_type='regular',
        axis_tick_label_height=0.4,
        axis_tick_label_colour='charcoal',
        axis_grid='off',
        axis_grid_colour='charcoal',
        axis_grid_thickness=1,
        axis_grid_line_style='dash',
    )

#	lines = ['Departure', ]

    switch='off'
    if ip==0:
        switch='on'
    graph = mgraph( legend=switch ,
                    graph_line_colour="black", 
                    graph_line_thickness= 5,
                    )

    sl=[]
    try:
        if 'annual' in plotproperties['dailytrends']:
            h=slope*((xtime[tmask])-numpy.mean(xtime[tmask]))/10.
            datat = minput(Input_x_values=[xtime[tmask][0],xtime[tmask][-1]],
                           Input_y_values=[h[0],h[-1]],
                           Input_y_missing_value=-1.e21
                           )
            grapht = mgraph( legend='off' ,
                             graph_line_colour="blue", 
                             graph_line_thickness= 5,
                             )
            sl=[datat,grapht]
    except:
        pass
    #if season:
        #datatw = minput(Input_x_values=xtime[maskwinter],
                        #Input_y_values=slopewinter*(xtime[maskwinter]-numpy.mean(xtime[maskwinter]))/10.,
                        #Input_y_missing_value=-1.e21
                        #)
        #graphtw = mgraph( legend='off' ,
                          #graph_line_colour="blue", 
                          #graph_line_thickness= 5,
                          #)
        #datats = minput(Input_x_values=xtime[masksummer],
                        #Input_y_values=slopesummer*(xtime[masksummer]-numpy.mean(xtime[masksummer]))/10.,
                        #Input_y_missing_value=-1.e21
                        #)
        #graphts = mgraph( legend='off' ,
                          #graph_line_colour="blue", 
                          #graph_line_thickness= 5,
                          #)
        #sl=sl+[datatw,graphtw,datats,graphts]

    legend_user_lines=['<font colour="black" size="0.4"> Difference</font>']
    legend = mlegend({"legend": "on", 
                      "legend_text_font_size":0.4, 
                      "legend_text_colour":"black",
                      "legend_box_mode":"positional",
                      "legend_box_x_position":2.0,
                      "legend_box_y_position":plotproperties['dailyseriesheight']*21.5/25.,
                      "legend_box_x_length":6.0,
                      "legend_box_y_length":plotproperties['dailyseriesheight']*2./25.,
                      "legend_column_count":1,
                      "legend_text_composition":"user_text_only",
                      "legend_user_lines":legend_user_lines
                      }) 
    line=['{0}hPa'.format(int(ps[pindex[ip]]))]
    text=mtext({
        "text_lines" : line,
        "text_html" : "true",
        "text_colour" : "black",
        "text_font_size" : 0.4,
        "text_mode" : "positional",
        "text_box_x_position": 2.0,
        "text_box_y_position":plotproperties['dailyseriesheight']*18.2/25.,
        "text_box_x_length": 10.,
        "text_box_y_length": 0.4,
        "text_border": "off",
        "text_justification" : "left"})

    line=['N={0:d} Std:{1:.2f}'.format(N,std)]
    stext=mtext({
        "text_lines" : line,
        "text_html" : "true",
        "text_colour" : "black",
        "text_font_size" : 0.4,
        "text_mode" : "positional",
        "text_box_x_position": 2.0,
        "text_box_y_position":plotproperties['dailyseriesheight']*4.5/25.,
        "text_box_x_length": 10.,
        "text_box_y_length": 0.4,
        "text_border": "off",
        "text_justification" : "left"})

    if season:
        line=['Trend {0}-{1}: {2:.2f} S{3:.2f} W{4:.2f} K/10a'.format(int(plotproperties['plotinterval'][0]),
                                                                      int(plotproperties['plotinterval'][1]),slope,slopesummer,slopewinter)]
    else:
        line=['Trend {0}-{1}: {2:.2f} K/10a'.format(int(plotproperties['plotinterval'][0]),
                                                    int(plotproperties['plotinterval'][1]),slope)]
    ttext=mtext({
        "text_lines" : line,
        "text_html" : "true",
        "text_colour" : "black",
        "text_font_size" : 0.4,
        "text_mode" : "positional",
        "text_box_x_position": 14.0,
        "text_box_y_position":plotproperties['dailyseriesheight']*18.2/25.,
        "text_box_x_length": 5.,
        "text_box_y_length": 0.4,
        "text_border": "off",
        "text_justification" : "right"})

    if ip==0:
        twelveminuszero=''
        if 'rad' in task['shortname']:
            twelveminuszero=', 12h-00h'
        if str(longname)=='missing':
            longname=''
        if type(istname)==numpy.bytes_:
            istname=istname.decode('utf-8')
        tline=['{0}, {1} {2}, {3:6.2f}N,{4:7.2f}E, {5:0>2}h'.format(str(task["name"])+twelveminuszero,
                    istname,str(longname),lat,lon,ipar*12)+', '+plotproperties['exps'][0]]
        print((lat,lon,istname,type(istname)))
    else:
        tline=['']

    title=mtext({
        "text_lines" : tline,
        "text_html" : "true",
        "text_colour" : "black",
        "text_font_size" : 0.5,
        "text_mode" : "positional",
        "text_box_x_position": 1.0,
        "text_box_y_position": plotproperties['dailyseriesheight']*25./25.,
        "text_box_x_length": 10.,
        "text_box_y_length": 0.5,
        "text_border": "off",
        "text_justification" : "left"})

    null_line=mgraph(legend='off',graph_line_colour='black')
    null_data = minput(Input_x_values=numpy.array(plotproperties["plotinterval"],dtype=numpy.float32),
                       Input_y_values=numpy.array([0,0],dtype=numpy.float32),
                       )

    rlist=[ppage,projection, vertical, horizontal, data, graph,legend,title,null_data,null_line]
    if not ("obs" in task["shortname"] or task["shortname"] in ["rcorr"] or task["shortname"][:3] in ['rio','rit']):
        tsa=numpy.zeros(tmean.shape[0])
        tsa.fill(numpy.nan)
        count=index.copy()
        #t=time.time()
#		snhtmov(hilf2,tsa,snhtparas,index,count,tmean,hilf)
#		print time.time()-t
        snhtmov2(hilf2,tsa,snhtparas,index,count,tmean,hilf)
#		print time.time()-t
        tsa[numpy.isnan(tsa)]=-1.e21
        index2=index.copy()
        n=thin(tsa,index2,plotproperties["thin"])

#		xt2=numpy.concatenate((xtime[index2[0:n]],xtime[index2[0:n]][::-1]))
        try:
            shadesnht=plotproperties['shadesnht']=='True'
        except:
            shadesnht=False
#		yt2=numpy.concatenate((tsa[index2[0:n]],tsacrop[::-1]-1))
        xt2=xtime[index2[0:n]]
        mask=numpy.logical_and(xt2>plotproperties['plotinterval'][0],xt2<plotproperties['plotinterval'][1])
        xt2=xt2[mask]	
        yt2=tsa[index2[0:n]][mask]
        yt2[numpy.isnan(yt2)]=-1.e21
        if sum(yt2>-1.e21)>1:
            xt2=xt2[yt2>-1.e20]
            if shadesnht:
                tsacrop=tsa[index2[0:n]][:]
                tsacrop[tsacrop>50.]=50.
                yt3=tsacrop[mask]
                yt3[numpy.isnan(yt3)]=-1.e21
                yt3=yt3[yt2>-1.e20]		
                yt2=yt2[yt2>-1.e20]
    
                data2 = minput(Input_x2_values=xt2,
                               Input_y2_values=yt3,
                               Input_x_values=xt2,
                               Input_y_values=yt2+2,		              
                               Input_y_missing_value=-1.e21,
                               Input_x_missing_value=-1.e21
                               )
            else:
                yt2=yt2[yt2>-1.e20]
                data2 = minput(Input_x_values=xt2,
                               Input_y_values=yt2+2,		              
                               Input_y_missing_value=-1.e21,
                               Input_x_missing_value=-1.e21
                               )
 
            projection2 = mmap(
                subpage_y_max=2.5*numpy.max([20.,numpy.max(tsa[index2[0:n]])]),
                subpage_y_min=-1.,
            )

            vertical2 = maxis(
                axis_orientation='vertical',
                axis_grid='off',
                axis_position='right',
                axis_tick_label_height=0.4,
                axis_title='on',
                axis_title_text='SNHT',
                #			axis_title_font='arial',
                axis_title_height=0.4,
            )
            if shadesnht:
                gsh='on'
                gt='area'
                glt=2
            else:
                gsh='off'
                gt='curve'
                glt=5
            graph2 = mgraph( legend=switch ,
                             graph_line_colour="red",
                             graph_line="on",
                             graph_type=gt,
                             graph_shade=gsh,
                             graph_shade_style='area_fill',
                             graph_shade_colour='red',
                             graph_line_thickness= glt,
                             graph_y_suppress_above=5000.,
                             graph_y_suppress_below=-10.,
                             )


            legend_user_lines=['<font colour="black" size="0.4"> SNHT</font>' ]
            legend2 = mlegend({"legend": "on", 
                               "legend_text_font_size":0.4, 
                               "legend_text_colour":"black",
                               "legend_box_mode":"positional",
                               "legend_box_x_position":8.0,
                               "legend_box_y_position":plotproperties['dailyseriesheight']*21.5/25,
                               "legend_box_x_length":6.0,
                               "legend_box_y_length":plotproperties['dailyseriesheight']*2./25,
                               "legend_column_count":1,
                               "legend_text_composition":"user_text_only",
                               "legend_user_lines":legend_user_lines
                               }) 

            rlist=rlist+[ppage,projection2,vertical2,data2,graph2,legend2]

    return [ppage,projection, vertical, horizontal, data, graph,legend,title,null_data,null_line]+sl+rlist+[text,stext,ttext]



if __name__ == "__main__":
    #settings of the png output 
    output = output(
        output_formats = ['png'],
        output_name = "precip_timeserie",
        output_name_first_page_number = "off"
    )

    plot(output, msl_timeserie(title_size = 0.8))

