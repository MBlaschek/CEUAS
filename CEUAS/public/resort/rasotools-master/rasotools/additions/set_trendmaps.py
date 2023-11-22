from Magics.macro import *
from rasotools.additions.utils import *

def set_trendmaps(lines,clev,plotproperties,stnames=None,slopes=None,costs=None,map='Globe'):
    
    contour_shade_colour_list=plotproperties['contour_shade_colour_list']
    if 'contour_shade_colour_method' not in plotproperties:
        cscm='list'
    else:
        cscm=plotproperties['contour_shade_colour_method']
    

    if map=='SEAsia':
        projection =mmap({"subpage_map_projection":"cylindrical", 
                          "subpage_y_length": 14.,
                          "subpage_y_position": 1.5 ,
                          "subpage_lower_left_longitude": 60.,
                          "subpage_lower_left_latitude":00.,
                          "subpage_upper_right_longitude":180.,
                          "subpage_upper_right_latitude":60.,
                          "page_id_line":"off",
                          "map_grid_latitude_increment":20.,
                          "map_grid_longitude_increment":30.,
                          "map_boundaries":'on',
                          "map_label_height":0.4})
    elif map=='Europe':
        projection =mmap({"subpage_map_projection":"cylindrical", 
                          "subpage_y_length": 14.,
                          "subpage_y_position": 1.5 ,
                          "subpage_lower_left_longitude":-50.,
                          "subpage_lower_left_latitude":30.,
                          "subpage_upper_right_longitude":60.,
                          "subpage_upper_right_latitude":85.,
                          "page_id_line":"off",
                          "map_grid_latitude_increment":20.,
                          "map_grid_longitude_increment":30.,
                          "map_boundaries":'on',
                          "map_label_height":0.4})
    elif map=='NAmerica':
        projection =mmap({"subpage_map_projection":"cylindrical", 
                          "subpage_y_length": 14.,
                          "subpage_y_position": 1.5 ,
                          "subpage_lower_left_longitude": -180.,
                          "subpage_lower_left_latitude":10.,
                          "subpage_upper_right_longitude":-40.,
                          "subpage_upper_right_latitude":80.,
                          "page_id_line":"off",
                          "map_grid_latitude_increment":30.,
                          "map_grid_longitude_increment":30.,
                          "map_boundaries":'on',
                          "map_label_height":0.4})
    elif map=='NP' or map=='SP':
        if map=='NP':
            hem='north'
            clat=90.
        else:
            hem='south'
            clat=-90.
        projection =mmap({"subpage_map_projection":"polar_stereographic", 
                          "subpage_y_length": 14.5,
                          "subpage_x_length": 14.5,
                          "subpage_y_position": 0.5 ,
                          "subpage_map_area_definition_polar": "centre",
                          "subpage_map_hemisphere":hem, 
                          "subpage_map_scale":50e6,
                          "subpage_map_centre_latitude":clat,
                          "subpage_map_vertical_longitude":0.,
                          "page_id_line":"off",
                          "map_grid_latitude_increment":20.,
                          "map_grid_longitude_increment":90.,
                          "map_label":'off',
                          "map_label_height":0.4})
    else: # Globe Cylindrical
        projection =mmap({"subpage_map_projection":"cylindrical", 
                          "subpage_y_length": 14.,
                          "subpage_y_position": 1.5 ,
                          "subpage_lower_left_longitude": -180.,
                          "subpage_lower_left_latitude":-90.,
                          "subpage_upper_right_longitude":180.,
                          "subpage_upper_right_latitude":90.,
                          "page_id_line":"off",
                          "map_grid_latitude_increment":30.,
                          "map_grid_longitude_increment":60.,
                          "map_label_height":0.4})
        

    #coastlines = mcoast(map_coastline_land_shade= 'on',
        #map_coastline_land_shade_colour= 'grey',
        #map_coastline_sea_shade= 'on',
        #map_coastline_sea_shade_colour= 'white',
    coastlines = mcoast(map_coastline_resolution='low')
    if map not in ('NP','SP'):
        txl=26.7
        lxl=26.7
        txp=1.5
        llf=1
    else:
        txl=15.5
        txp=0.5
        lxl=13.5
        llf=2
        
    title = mtext({
        "text_lines" : lines,
        "text_html" : "true",
        "text_colour" : "black",
        "text_font_size" : 0.6,
        "text_mode" : "positional",
        "text_box_x_position": txp,
        "text_box_y_position": 17.5,
        "text_box_x_length": txl,
        "text_box_y_length": 1.5,
        "text_border": "on",
        "text_justification" : "left"})

    legend = mlegend( {"legend": "on", 
                       "legend_display_type": "continuous",
                       "legend_label_frequency":llf,
                       "legend_text_font_size":0.5, 
                       "legend_text_colour":"black",
                       "legend_box_mode":"positional",
                       "legend_box_x_position":1.5,
                       "legend_box_y_position":15.5,
                       "legend_box_x_length":lxl,
                       "legend_box_y_length":1.2
                       })

    symb=msymb(symbol_advanced_table_selection_type = "list",
               symbol_marker_index = 27,
               symbol_table_mode = "advanced",
               legend = "on",
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
               symbol_advanced_table_colour_method = cscm,
               symbol_advanced_table_colour_list = contour_shade_colour_list,
               
               symbol_advanced_table_colour_direction = "clockwise",
               symbol_advanced_table_min_level_colour = "#ebf5fb",
               symbol_advanced_table_max_level_colour = "#1b4f72"
               )

    cont=mcont(contour_shade='on',
               contour_shade_method='area_fill',
               contour_shade_technique='grid_shading',
               contour_level_selection_type='list',
               contour_level_list=clev,
               contour_label='off',
               contour='off',
               contour_grid_shading_position='bottom_left',
               contour_shade_colour_list=contour_shade_colour_list,
               contour_shade_colour_method='list'               
               #contour_shade_colour_direction = "clockwise",
               #contour_shade_min_level_colour = "reddish_purple",
               #contour_shade_max_level_colour = "red_orange"
               )
    
    if stnames is not None:
        #tl=''
        #for tx in stnames[-1:0:-1]:
            #tl+=tx+'/'
        #tl+=stnames[0]
            
        symb2=msymb(symbol_table_mode = "off",
                   symbol_colour='yellow',#'rgb(1.0,1.0,1.0)',
                   legend = "off",
                   symbol_type = "text",
                   symbol_text_list=list(stnames),
                   symbol_text_position='right',
                   symbol_text_font_size=0.4,
                   symbol_text_font_colour='black',
                   symbol_outline='off',
                   symbol_outline_colour='yellow',
                   symbol_marker_index=27
                   )
        
        tcosts=list()
        for i in range(costs.shape[0]):
            tcosts.append('{0:4.1f}'.format(costs[i]))
        symb3=msymb(symbol_advanced_table_selection_type = "list",
                   symbol_table_mode = "off",
                   symbol_height = 0.01,
                   symbol_colour='yellow',#rgb(1.0,1.0,1.0)',
                   legend = "off",
                   symbol_type = "text",
                   symbol_text_list=tcosts,
                   symbol_text_position='right',
                   symbol_text_font_size=0.2,
                   symbol_text_font_colour='black',
                   symbol_outline='off',
                   symbol_outline_colour='yellow',
                   symbol_marker_index=27
                   )
    
        tslopes=list()
        for i in range(slopes.shape[0]):
            tslopes.append('{0:4.1f}'.format(slopes[i]))
        symb4=msymb(symbol_advanced_table_selection_type = "list",
                   symbol_table_mode = "off",
                   symbol_colour='yellow',#rgb(1.0,1.0,1.0)',
                   legend = "off",
                   symbol_type = "text",
                   symbol_text_list=tslopes,
                   symbol_text_position='right',
                   symbol_text_font_size=0.2,
                   symbol_text_font_colour='black',
                   symbol_outline='off',
                   symbol_outline_colour='yellow',
                   symbol_marker_index=27
                   )

    if stnames is not None:
        return projection,coastlines,title,legend,symb,symb2,symb3,symb4,cont
    else:
        return projection,coastlines,title,legend,symb,cont

def set_zonalsat(lines,legend_user_lines,plotproperties,beltinterval=[-2.0,0.5],pkey=''):

    
    apl=[-2.0,-1.0,0.0]
    vaxistype='regular'
    vaxpl=[0,1]
    if pkey!='cost':
        try:
            vaxpl=plotproperties['haxpl']
            vaxistype=plotproperties['haxistype']
        except:
            pass
    
    projection = mmap(
        subpage_map_projection='cartesian',
        subpage_x_min=-90.,
        subpage_x_max=90.,
        subpage_y_min=-1.,
        subpage_y_max=1.,
        subpage_y_length=14.,
        subpage_y_position=1.5 ,
    )
    # Vertical axis
    vxd=dict(
        axis_orientation='vertical',
        axis_type=vaxistype,
        axis_tick_label_height=0.4,
        axis_tick_label_colour='black',
        axis_tick_position_list=numpy.asarray(vaxpl,dtype='float32'),
        axis_grid='on',
        axis_grid_colour='black',
        axis_grid_thickness=1,
        axis_grid_line_style='dot')
    
    if vaxistype=='position_list':
        hxl=[];
        for h in vaxpl:
            hxl.append('{:4.1f}'.format(h))
        vxd['axis_tick_label_list']=hxl
        
    vertical = maxis(vxd)
    # Horizontal axis
    horizontal = maxis(
        axis_orientation='horizontal',
        axis_type='regular',
        axis_tick_label_height=0.4,
        axis_tick_label_colour='black',
        axis_grid='on',
        axis_grid_colour='black',
        axis_grid_thickness=1,
        axis_grid_line_style='dash',
    )
    
    if len(legend_user_lines)==1:
        legend_userlines=legend_user_lines*2
    legend = mlegend( {"legend": "on", 
                       "legend_text_font_size":0.5, 
                       "legend_text_colour":"black",
                       "legend_box_mode":"positional",
                       "legend_box_x_position":0.0,
                       "legend_box_y_position":15.5,
                       "legend_box_x_length":29.7,
                       "legend_box_y_length":1.8,
                       "legend_text_composition":"user_text_only",
                       "legend_column_count":4,
                       "legend_user_lines":legend_user_lines
                       })

    title = mtext(
        text_lines= lines,
        text_html= 'true',
        text_colour= 'black',
        text_font_size= 0.6,
        text_mode = 'positional',
        text_box_x_position= 1.5,
        text_box_y_position= 17.5,
        text_box_x_length= 20.,
        text_box_y_length= 1.5,
        text_border= "off",
        text_justification = "left"    )

    return projection,horizontal,vertical,title,legend    
def set_zonalcross(lines,clev,plotproperties):
    
    
    contour_shade_colour_list=plotproperties['contour_shade_colour_list']
    logp=plotproperties['logp'] 
    try:
        ymax=plotproperties['ymax']
    except:
        ymax=10.
        print('ymax not defined!!')
    if 'contour_shade_colour_method' not in plotproperties:
        cscm='list'
    else:
        cscm=plotproperties['contour_shade_colour_method']
    
    
    xtype='regular'
    if logp==True:
        xtype='logarithmic'
    
    print(('zlogp',type(logp),logp,xtype))
    
    projection = mmap(
        subpage_map_projection='cartesian',
        subpage_x_min=-90.,
        subpage_x_max=90.,
        subpage_y_min=1000.,
        subpage_y_max=ymax,
        subpage_y_axis_type = xtype,
        subpage_y_length=14.,
        subpage_y_position=1.5 ,
    )
    # Vertical axis
    vertical = maxis(
        axis_orientation='vertical',
        axis_grid='on',
        axis_type=xtype,
        axis_tick_label_height=0.4,
        axis_tick_label_colour='black',
        axis_grid_colour='black',
        axis_grid_thickness=1,
        axis_grid_reference_line_style='solid',
        axis_grid_reference_thickness=1,
        axis_grid_line_style='dash',
        axis_title='on',
        axis_title_text='Pressure',
        axis_title_height=0.6,
    )
    # Horizontal axis
    horizontal = maxis(
        axis_orientation='horizontal',
        axis_type='regular',
        axis_tick_label_height=0.4,
        axis_tick_label_colour='black',
        axis_grid='on',
        axis_grid_colour='black',
        axis_grid_thickness=1,
        axis_grid_line_style='dash',
    )
    
    legend = mlegend( {"legend": "on", 
                       "legend_display_type": "continuous",
                       "legend_text_font_size":0.5, 
                       "legend_text_colour":"black",
                       "legend_box_mode":"positional",
                       "legend_box_x_position":0.0,
                       "legend_box_y_position":15.7,
                       "legend_box_x_length":29.7,
                       "legend_box_y_length":0.9
                       })

    title = mtext(
        text_lines= lines,
        text_html= 'true',
        text_colour= 'black',
        text_font_size= 0.6,
        text_mode = 'positional',
        text_box_x_position= 1.5,
        text_box_y_position= 17.5,
        text_box_x_length= 20.,
        text_box_y_length= 1.5,
        text_border= "on",
        text_justification = "left"    )
    
    cont=mcont(contour_shade="on",
            #contour_shade_technique='cell_shading',
            contour_shade_method='area_fill',
            #contour_shade_technique='grid_shading',
            contour_level_selection_type='list',
            contour_level_list=clev,
            contour_label='off',
            contour='off',
            contour_grid_shading_position='bottom_left',
            contour_shade_colour_list=contour_shade_colour_list,
            contour_shade_colour_method=cscm,
            #contour_shade_colour_direction = "clockwise",
            #contour_shade_min_level_colour = "blue",
            #contour_shade_max_level_colour = "red",
            ##contour_shade_min_level_colour = "ebf5fb",
            ##contour_shade_max_level_colour = "1b4f72",
           legend='on'
            )
    #cont=mcont()

    return projection,horizontal,vertical,title,legend,cont

def set_belttrends(lines,legend_user_lines,plotproperties,pkey='beltslopes',beltinterval=[-1.8,0.7],haxistype=''):

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
            if haxistype is 'regular':
                pass
            else:
                haxpl=plotproperties['haxpl']
                haxistype=plotproperties['haxistype']
        except:
            pass
        
    projection = mmap(
        subpage_map_projection='cartesian',
        subpage_x_position=2.0,
        subpage_x_axis_type='regular',
        subpage_y_axis_type=ptype,
        subpage_x_min=beltinterval[0],
        subpage_x_max=beltinterval[1],
        subpage_y_min=1000.,
        subpage_y_max=20.,
        subpage_y_length=14.,
        subpage_y_position=1.5 ,
    )
    # Vertical axis
    vertical = maxis(
        axis_orientation='vertical',
        axis_grid='on',
        axis_type=ptype,
        axis_tick_label_height=0.6,
        axis_tick_label_colour='black',
        axis_grid_colour='black',
        axis_grid_thickness=1,
        axis_grid_reference_line_style='solid',
        axis_grid_reference_thickness=1,
        axis_grid_line_style='dash',
        axis_title='on',
        axis_title_text='Pressure',
        axis_title_height=0.6,
    )
    # Horizontal axis
    hxd=dict(
        axis_orientation='horizontal',
        axis_type=haxistype,
        axis_tick_label_height=0.6,
        axis_tick_label_colour='black',
        axis_tick_position_list=numpy.asarray(haxpl,dtype='float32'),
        axis_grid='on',
        axis_grid_colour='black',
        axis_grid_thickness=1,
        axis_grid_line_style='dot')
    
    if haxistype=='position_list':
        hxl=[];
        for h in haxpl:
            hxl.append('{:4.1f}'.format(h))
        hxd['axis_tick_label_list']=hxl
        
    horizontal = maxis(hxd)
    
    if len(legend_user_lines)==1:
        legend_user_lines=legend_user_lines*2
    legend = mlegend( {"legend": "on", 
                       "legend_text_font_size":0.5, 
                       "legend_text_colour":"black",
                       "legend_box_mode":"positional",
                       "legend_box_x_position":0.0,
                       "legend_box_y_position":15.5,
                       "legend_box_x_length":29.7,
                       "legend_box_y_length":1.8,
                       "legend_text_composition":"user_text_only",
                       "legend_column_count":4,
                       "legend_user_lines":legend_user_lines
                       })

    title = mtext(
        text_lines= lines,
        text_html= 'true',
        text_colour= 'black',
        text_font_size= 0.7,
        text_mode = 'positional',
        text_box_x_position= 1.5,
        text_box_y_position= 17.5,
        text_box_x_length= 20.,
        text_box_y_length= 1.5,
        text_border= "off",
        text_justification = "left"    )
    
    return projection,horizontal,vertical,title,legend

def monthlyseries(lines,legend_user_lines,plotproperties,dynrange=[-1.,1.]):

    projection = mmap(
        subpage_map_projection='cartesian',
        subpage_x_axis_type='regular',
        subpage_y_axis_type='regular',   
        subpage_x_position=1.4,
        subpage_y_position=1.0,
        subpage_x_length=22.8,
        subpage_y_length=plotproperties['monthlyseriesheight'],
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
        axis_tick_label_height=0.6,
        axis_grid_reference_level=0.,
        axis_grid_reference_thickness=1,
        axis_grid_reference_colour='black',
        axis_title='off',
        axis_title_text=plotproperties['units'][0],
        axis_title_height=0.6,
    )

    # Horizontal axis

    horizontal = maxis(
        axis_orientation='horizontal',
        axis_type='regular',
        axis_tick_label_height=0.6,
        axis_tick_label_colour='black',
        axis_grid='off',
        axis_grid_colour='black',
        axis_grid_thickness=1,
        axis_grid_line_style='dash',
    )

#	lines = ['Departure', ]

    #switch='off'
    #if ip==0:
        #switch='on'
    #graph = mgraph( legend='on' ,
                    #graph_line_colour="black", 
                    #graph_line_thickness= 5,
                    #)

    if len(legend_user_lines)==1:
        legend_user_lines=legend_user_lines*2
    legend = mlegend( {"legend": "on", 
                       "legend_text_font_size":0.5, 
                       "legend_text_colour":"black",
                       "legend_box_mode":"positional",
                       "legend_box_x_position":0.0,
                       "legend_box_y_position":plotproperties['monthlyseriesheight']+0.4,
                       "legend_box_x_length":25.0,
                       "legend_box_y_length":1.8,
                       "legend_text_composition":"user_text_only",
                       "legend_column_count":8,
                       "legend_user_lines":legend_user_lines
                       })

    title = mtext(
        text_lines= lines,
        text_html= 'true',
        text_colour= 'black',
        text_font_size= 0.6,
        text_mode = 'positional',
        text_box_x_position= 1.5,
        text_box_y_position= plotproperties['monthlyseriesheight']+2.,
        text_box_x_length= 20.,
        text_box_y_length= 1.5,
        text_border= "off",
        text_justification = "left"    )

    null_line=mgraph(legend='off',graph_line_colour='black')
    null_data = minput(Input_x_values=numpy.array(plotproperties["plotinterval"],dtype=numpy.float32),
                       Input_y_values=numpy.array([0,0],dtype=numpy.float32),
                       )

    return [projection, horizontal, vertical,title, legend,null_data,null_line]

def sbtrends(lines,plotproperties,dynrange=[-1.,1.], second=False):

#    global base,date_list,date_listm
    # Setting the cartesian view

    if second:
        ppage = page(
            layout='positional',  
            page_x_length=4.9, 
            page_y_length=plotproperties['monthlyseriesheight']+3.0, 
            page_id_line='off',
            page_x_position=30.4, 
            page_y_position=0.)
    else:
    
        ppage = page(
            layout='positional',  
            page_x_length=4.9, 
            page_y_length=plotproperties['monthlyseriesheight']+3.0, 
            page_id_line='off',
            page_x_position=24.7, 
            page_y_position=0.)
    

    projection = mmap(
        subpage_map_projection='cartesian',
        subpage_x_axis_type='regular',
        subpage_y_axis_type='regular',   
        subpage_x_position=1.0,
        subpage_y_position=1.0,
        subpage_x_length=4.0,
        subpage_y_length=plotproperties['monthlyseriesheight'],
        subpage_x_max=4.,
        subpage_x_min=0.,
        subpage_y_max=dynrange[1],
        subpage_y_min=dynrange[0],
    )

    # Vertical axis

    vertical = maxis(
        axis_orientation='vertical',
        axis_grid='off',
        axis_type='regular',
        axis_tick_label_height=0.6,
        axis_grid_reference_level=0.,
        axis_grid_reference_thickness=1,
        axis_grid_reference_colour='black',
        axis_title='off',
        axis_title_text='K/10a',
        axis_title_height=0.6,
    )

    # Horizontal axis

    horizontal = maxis(
        axis_orientation='horizontal',
        axis_type='position_list',
        axis_tick_position_list=[0.6,2.,3.4],
        axis_tick_label_height=0.6,
        axis_tick_label_colour='black',
        axis_tick_label_type='label_list',
        axis_tick_label_list=['ROB','REA','SAT'],
        axis_grid='off',
        axis_grid_colour='black',
        axis_grid_thickness=1,
        axis_grid_line_style='dash',
    )

    title = mtext(
        text_lines= lines,
        text_html= 'true',
        text_colour= 'black',
        text_font_size= 0.6,
        text_mode = 'positional',
        text_box_x_position= 1.0,
        text_box_y_position= plotproperties['monthlyseriesheight']+1.3,
        text_box_x_length= 4.2,
        text_box_y_length= 1.5,
        text_border= "off",
        text_justification = "left"    )

    null_line=mgraph(legend='off',graph_line_colour='black')
    null_data = minput(Input_x_values=numpy.array([0,4],dtype=numpy.float32),
                       Input_y_values=numpy.array([0,0],dtype=numpy.float32),
                       )

    return [ppage,projection, horizontal, vertical,title,null_data,null_line]#+sl+rlist+[text,stext,ttext]

def sbmarkers(tasks,keys,gsls):
    
    print('sbmarkers')
    pos=['ROB','REA','SAT']
    sh=dict()
    for p in pos:
        sh[p]=-0.1
    idx=[15,17,18]
    l=[]
    for k in range(len(keys)):
        
        tk=tasks[keys[k]]['type']
        sh[tk]+=0.1
        xpos=float(pos.index(tk))*1.4+0.6+sh[tk]
        print((k,tasks[keys[k]]['shortname'],xpos,gsls[k][0]))
        data = minput(Input_x_values=[xpos],
          Input_y_values=[gsls[k][0]],
          Input_y_missing_value=-1.e21
          )
    
        graph = msymb( symbol_type= 'marker',
                       symbol_marker_mode='index',
                       symbol_marker_index=idx[pos.index(tk)],
                       symbol_colour=tasks[keys[k]]['color'], 
                       symbol_height=0.8, 
                        legend="off"
                        )
        l+=  [data,graph]  
        
        if(True):
            data = minput(Input_x_values=[xpos-0.2,xpos+0.2],
              Input_y_values=([gsls[k][2],gsls[k][2]]),
              Input_y_missing_value=-1.e21
              )
        
            graph = mgraph(graph_line_colour= tasks[keys[k]]['color'], 
                        graph_line_thickness= 2
                            )
            
            l+=  [data,graph]  
            data = minput(Input_x_values=[xpos,xpos],
              Input_y_values=[gsls[k][2],gsls[k][1]],
              Input_y_missing_value=-1.e21
              )
        
            graph = mgraph(graph_line_colour= tasks[keys[k]]['color'], 
                        graph_line_thickness= 2
                            )
            
            l+=  [data,graph]  
            data = minput(Input_x_values=[xpos-0.2,xpos+0.2],
              Input_y_values=[gsls[k][1],gsls[k][1]],
              Input_y_missing_value=-1.e21
              )
        
            graph = mgraph(graph_line_colour= tasks[keys[k]]['color'], 
                        graph_line_thickness= 2
                            )
            
            l+=  [data,graph]  
        
    return l

def add_mts(tasks,i,ip,iens,ibelt,iref,plotproperties,scol=None,supp=[-5.,5.]):
    
    if scol is not None:
        col=scol
    else:
        col=task["color"]
    if iref==-1:
        ds=tasks[0]['beltanomalies'].shape[2]-tasks[i]['beltanomalies'].shape[2]
        raw=tasks[i]['beltanomalies'][iens,ibelt,ip-ds,:]
        hilf=rmeanw(raw,3)
    else:
        ds=tasks[iref]['beltanomalies'].shape[2]-tasks[i]['beltanomalies'].shape[2]
        raw=tasks[i]['beltanomalies'][iens,ibelt,ip-ds,:]-tasks[iref]['beltanomalies'][0,ibelt,ip,:]
        if tasks[i]['shortname']=='eraibc':
            hilf=rmeanw(raw,13)
        else:
            hilf=rmeanw(raw,3)
       
    hilf[hilf!=hilf]=-1.e21
    time=tasks[i]['startdate']//10000+1./24.+numpy.arange(hilf.shape[0])/12.
    mask=numpy.logical_and(time>=plotproperties['plotinterval'][0],time<plotproperties['plotinterval'][1])
    data = minput(Input_x_values=time[mask],
      Input_y_values=hilf[mask],
      Input_x_missing_value=-1.e21
      )

    graph = mgraph( graph_line_colour= col, 
                    graph_line_thickness= 6,
                    #graph_y_suppress_above=supp[1],
                    #graph_y_suppress_below=supp[0],
                    legend="on",
                    )
    return [data,graph],raw,time

def addprofile(task,ps,iens,ibelt,pkey="beltslopes",col='black'):
    
    if type(task) is dict:
        hilf=task[pkey][iens,ibelt,:].flatten()
        col=task["color"]
    else:
        hilf=task.flatten()
        col=col
#    if sum(hilf==hilf)==0:
#        return []
    mask=~numpy.isnan(hilf[:ps.shape[0]])
    #hilf[mask]=-1.e21
    data = minput(Input_x_values=hilf[:ps.shape[0]][mask],
      Input_y_values=ps[mask],
      Input_x_missing_value=-1.e21
      )
    if pkey=='cost':
        supp=[0.,10000.]
    else:
        supp=[-10.,10.]

    graph = mgraph( graph_line_colour= col, 
                    graph_line_thickness= 12,
                    graph_x_suppress_above=supp[1],
                    graph_x_suppress_below=supp[0],
                    legend="on",
                    )
    return [data,graph]

def addsatline(task,lat,iens,ilayer,pkey='zslopes'):
    
    hilf=task[pkey][iens,:,ilayer].flatten()
#    if sum(hilf==hilf)==0:
#        return []
    mask=numpy.isnan(hilf)
    hilf[mask]=-1.e21
    data = minput(Input_y_values=hilf,
      Input_x_values=lat,
      Input_y_missing_value=-1.e21
      )
    if pkey=='cost':
        supp=[0.,10000.]
    else:
        supp=[-10.,10.]

    graph = mgraph( graph_line_colour= task["color"], 
                    graph_line_thickness= 12,
                    graph_y_suppress_above=supp[1],
                    graph_y_suppress_below=supp[0],
                    legend="on",
                    )
    return [data,graph]

def add_ens_profile(task,ps,ibelt,pkey='beltslopes'):
    
    
    esiz=task[pkey].shape[0]
    estart=numpy.int(numpy.round(esiz*0.025))
    estop=numpy.int(numpy.round(esiz*0.975))
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
                    graph_x_suppress_above=supp[1],
                    graph_x_suppress_below=supp[0],
                    legend="on",
                    )
    return [data,graph]

def add_ens_satline(task,lat,ilayer,pkey='zslopes'):
    
    
    esiz=task[pkey].shape[0]
    estart=numpy.round(esiz*0.025)
    estop=numpy.round(esiz*0.975)
    if estop==esiz:
        estop=esiz-1
    ens05=numpy.zeros(lat.shape[0])
    ens95=numpy.zeros(lat.shape[0])
    
    for ip in range(lat.shape[0]):
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
    data = minput(Input_y_values=ens05,
                  Input_y2_values=ens95,
                  Input_x_values=lat,
                  Input_x2_values=lat,
      Input_y_missing_value=-1.e21
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
                    graph_y_suppress_above=supp[1],
                    graph_y_suppress_below=supp[0],
                    legend="on",
                    )
    return [data,graph]

def addhadtrend(medtrend,enstrends):
    
    
    p=numpy.asarray([990.,990.])
    data = minput(Input_x_values=[medtrend.flatten()[0],medtrend.flatten()[0]],
      Input_y_values=p,
      )

    symb = msymb(   symbol_marker_index= 28,
                    symbol_colour='black',
                    symbol_height= 0.5,
                    symbol_type="marker",
                    legend="on",
                    legend_user_text='HadCRUT5'
                    )
    estart=2
    estop=97
    index=numpy.argsort(enstrends.flatten())
    ens05=numpy.asarray([enstrends[index[estart]],enstrends[index[estart]],enstrends[index[estop]],
                         enstrends[index[estop]],enstrends[index[estop]],enstrends[index[estop]]])
    
    mask=numpy.isnan(ens05)
    if sum(mask)!=0:

        ens05[mask]=-1.e21
    dataens=minput(Input_x_values=ens05,
                 Input_y_values=numpy.asarray([970.,1010.,990.,990.,970.,1010.]),
                 Input_x_missing_value=-1.e21
                 )
    graph = mgraph( graph_line_colour= 'black', 
                    graph_line_style = "solid",
                    graph_line_thickness = 2,
                    legend="off",
                    )
    return [data,symb,dataens,graph]

def addsattrend(medtrend,p,scol,iens=0,enstrends=None):
    
    
    x=medtrend.flatten()[:]
#    if sum(medtrend==medtrend)==0:
#        return []
    x[numpy.isnan(x)]=-1.e21
    data = minput(Input_x_values=x,
      Input_y_values=p,
      )

    leg='on'
    if iens>0:
        leg='off'
        
    symb = msymb(   symbol_marker_index= 28,
                    symbol_colour=scol,
                    symbol_height= 0.5,
                    symbol_type="marker",
                    legend=leg,
                    )
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
    return [data,symb]

def northof30N(tasks,key,corig,path,ps,pindex):
    ylim=dict()
    ylim['rcorr']=[-1,3]
    ylim['bgdiff']=[-2,5]
    tit={'rcorr':'adjustment','bgdiff':'obs-bg'}
    xlim={'05':[1958,2013],'14':[1935,2015]}
    xstart={'05':1958,'14':1900}
    mask=(lats<30.)
    s[mask,:,:,:]=numpy.nan
    anomalies[mask,:,:,:]=numpy.nan
    rmsanomalies=numpy.empty([corig.shape[3],corig.shape[2]])
    meananomalies=numpy.empty([corig.shape[3],corig.shape[2]])
    countanomalies=numpy.empty([corig.shape[3],corig.shape[2]],dtype='int')
    for ip in range(corig.shape[2]):
        for it in range(corig.shape[3]):
            countanomalies[it,ip]=numpy.sum(~numpy.isnan(corig[:,:,ip,it]))
            meananomalies[it,ip]=numpy.nanmean(corig[:,:,ip,it])
            rmsanomalies[it,ip]=numpy.nanstd(corig[:,:,ip,it])
            rmsanomalies[it,ip]=numpy.sqrt(rmsanomalies[it,ip]*rmsanomalies[it,ip]+meananomalies[it,ip]*meananomalies[it,ip])
    
        font = {'family' : 'sans-serif',
                'weight' : 'normal',
                'size'   : 10}
    
        plt.rc('font', **font)
    
        plt.figure(figsize=[20/2.54,8/2.54])
        lw=3
        plt.subplot(2,1,1)    
        l1,=plt.plot(xstart[path[-2:]]+1./24+numpy.arange(anomalies.shape[3])/12.,
                     rmsanomalies[:,ip],label='rms '+tit[tasks[key]["shortname"]],linewidth=lw)
        l2,=plt.plot(xstart[path[-2:]]+1./24+numpy.arange(anomalies.shape[3])/12.,
                     meananomalies[:,ip],'r',label='mean '+tit[tasks[key]["shortname"]],linewidth=lw)
        plt.plot(xlim[path[-2:]],[0.,0.],'k')
        if ip >0:
            plt.ylim(ylim[tasks[key]["shortname"]])
        else:
            yl=ylim[tasks[key]["shortname"]]
            plt.ylim([yl[0],yl[1]+1])
    
        plt.xlim(xlim[path[-2:]])
        plt.ylabel('K')
        plt.legend(handles=[l1,l2],frameon=False)
        plt.title('Radiosondes north of 30N, {0} hPa'.format(int(ps[pindex[ip]])))
        plt.subplot(2,1,2)    
        plt.plot(xstart[path[-2:]]+1./24+numpy.arange(anomalies.shape[3])/12.,
                 countanomalies[:,ip],'g',label='Data count',linewidth=lw)
        plt.xlim(xlim[path[-2:]])
        plt.ylabel('Data Count')
        plt.savefig(tasks[key]["shortname"]+'_'+path[-2:]+'_{0}'.format(int(ps[pindex[ip]]))+'.eps')
