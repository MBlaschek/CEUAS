import os,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import pandas as pd
from datetime import date, datetime
import cartopy.crs as ccrs




# d = '1965-03-04'
# datetime_obj = datetime.strptime(d, '%Y-%m-%d')

os.system('mkdir out_plots')

def plot_properties(flag = ''):
    
    f=plt.figure(figsize=(10,5))

    
    if flag != 'pretty':
        
        ax=plt.subplot(1,1,1,projection=ccrs.Robinson())  
        '''
        ax=plt.subplot(1,1,1,projection=ccrs.PlateCarree())  
        gl.xlabels_top = False
        gl.ylabels_right = False        
        gl=ax.gridlines(draw_labels=True , ls = ':', color = 'lightgray')
        plt.grid(ls = '--', color = 'lightgray')  
        '''        
        ax.stock_img()        


    elif flag == 'pretty': 
        ax = plt.axes(projection=ccrs.Mollweide())
        ax.stock_img()

    ax.coastlines(color = 'black' , linewidth=0.3)

    #ax.coastline()
    #plt.xlim([-180.,180.])
    #plt.ylim([-90.,90.])

    """ dictionary defining the style of the points """   
    d = {'ncar'     : ['darkorange' , 4 , 'o'] , 
         'bufr'     : ['red'    , 3 , 'o'] ,
         'igra2'    : ['mediumpurple'   , 6 , 'o'] ,
         'era5_1'   : ['lime' , 3  , 'o'] , 
         'era5_1759': ['blue'   , 13 , 'o'] , 
         'era5_1761': ['cyan'   , 8 , 'o'] ,
         'era5_3188': ['yellow' , 2 , 'o'] ,
 } 




    
    return d

def clean_data(lat, lon, start, end, date_min = '', date_max=''):
    """ Extracts only the info to be plotted, between the chosen date_min and date_max"""
    
    date_min = datetime.strptime( date_min , '%Y-%m-%d')
    date_max= datetime.strptime( date_max , '%Y-%m-%d') # convert to datetime objects
    
    latc, lonc= [], []
    
    for i in range(len(lat)):
        st = start[i]  
        en = end[i]
        if start[i] == '-':  continue 
        
        st = datetime.strptime( st , '%Y-%m-%d')
        en =  datetime.strptime( en , '%Y-%m-%d')
        if st > date_max : continue
        if en < date_min: continue
        
        lon_check = float(lon[i])
        if lon_check > 180:
            lon_check = lon_check - 360 
        
        latc.append(float(lat[i]))
        lonc.append(lon_check)

    return latc, lonc 
    
       
def makePlot(start_date= '', end_date = '', file = '', flag = ''):
    """ Main function to plot maps. 
         [start, end]: time range for the data 
         file+ summary of all the station_configuration files merged together """

    print('Plotting :::' , start_date, ' ' , end_date , ' ' , flag )

    # check the column names order in the merging_netcdf.py script     
    col_names = ['ncar_lat', 'ncar_lon', 'ncar_start', 'ncar_end',    
                             'igra2_lat', 'igra2_lon', 'igra2_start', 'igra2_end'  ,    
                             'bufr_lat', 'bufr_lon', 'bufr_start', 'bufr_end', 
     
                             'era5_1_lat', 'era5_1_lon', 'era5_1_start', 'era5_1_end',      
                             'era5_1759_lat', 'era5_1759_lon', 'era5_1759_start', 'era5_1759_end', 
                             'era5_1761_lat', 'era5_1761_lon', 'era5_1761_start', 'era5_1761_end',
                             'era5_3188_lat', 'era5_3188_lon', 'era5_3188_start', 'era5_3188_end' ]
                    

    tab =pd.read_csv(file, delimiter=',',  na_filter=False, comment='#' , names = col_names )
       
       
    dic = plot_properties(flag = flag)
    
    for d in ['era5_1759', 'era5_1761', 'igra2', 'ncar', 'era5_3188', 'bufr', 'era5_1']:

        lat, lon =  clean_data( tab[ d + '_lat'], tab[d + '_lon'],  tab[d + '_start'] , tab[d + '_end'] , date_min = start_date, date_max= end_date )  # cleaning the station to be plotted                                                                                    

        if len(lon) <1: continue 

        color  = dic[d][0]
        marker = dic[d][2]
        size   = dic[d][1]
        lab = d.replace('era5_1','ERA5 1').replace('era5_3188','ERA5 3188').replace('era5_1759','ERA5 1759').replace('era5_1761','ERA5 1761').replace('igra2','IGRA2').replace('ncar','NCAR').replace('bufr','BUFR')         

        plt.scatter (lon, lat,   color = color, transform=ccrs.PlateCarree(central_longitude = 0.0), s = size, label = lab + ' [' + str(len(lon)) + ']' )
       
    plt.title('Data Availability from ' + start_date + ' to ' + end_date , y = 1.04 )   
    plt.legend(loc = 'lower left', ncol = 1, fontsize = 12)   
    plt.savefig('out_plots/map_' + start_date + '_' + end_date + '_' + flag + '.png', dpi = 300 , transparent=True )
    plt.close()
     
       
       

############################# Plotting part       
       
for f in ['', 'pretty']:
    
    a = makePlot(start_date = '1940-01-01', end_date = '1950-01-01',   file = 'summary_forplot.dat', flag = f )
   
    '''
    a = makePlot(start_date = '1900-01-01', end_date = '1920-01-01',   file = 'summary_forplot.dat', flag = f )
    a = makePlot(start_date = '1920-01-01', end_date = '1940-01-01',   file = 'summary_forplot.dat', flag = f ) 
    a = makePlot(start_date = '1940-01-01', end_date = '1950-01-01',   file = 'summary_forplot.dat', flag = f ) 
    a = makePlot(start_date = '1955-01-01', end_date = '1960-01-01',   file = 'summary_forplot.dat', flag = f )
    a = makePlot(start_date = '1965-01-01', end_date = '1970-01-01',   file = 'summary_forplot.dat', flag = f )
    a = makePlot(start_date = '1975-01-01', end_date = '1980-01-01',   file = 'summary_forplot.dat', flag = f )
    a = makePlot(start_date = '1980-01-01', end_date = '1990-01-01',   file = 'summary_forplot.dat', flag = f )
    a = makePlot(start_date = '1990-01-01', end_date = '2000-01-01',   file = 'summary_forplot.dat', flag = f )
    a = makePlot(start_date = '2000-01-01', end_date = '2010-01-01',   file = 'summary_forplot.dat', flag = f )
    a = makePlot(start_date = '2010-01-01', end_date = '2020-01-01',   file = 'summary_forplot.dat', flag = f )
    '''