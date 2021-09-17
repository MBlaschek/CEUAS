""" Produce a Hovmoeller plot """

import os,sys
import h5py as h5
import pandas as pd
import netCDF4 as nc

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker

import matplotlib as mpl
mpl.rcParams['agg.path.chunksize'] = 1000000



"""
dic = { 85 : {'x':'Temperature [K]'       , 'y': 'Pressure [hPa]' , 'ylim': (180,340), 'ylimdes': (0,10) } ,

        106 : {'x': 'Wind from Direction'  , 'y': 'Pressure [hPa]', 'ylim': (180,340), 'ylimdes': (0,14) },
        107 : {'x': 'Wind Speed [m/s] '    , 'y': 'Pressure [hPa]', 'ylim': (0,100), 'ylimdes': (0,30) },
        104 : {'x': 'Wind u-component [m/s] '    , 'y': 'Pressure [hPa]', 'ylim': (-50,100), 'ylimdes': (0,30) },
        105 : {'x': 'Wind v-component [m/s] '    , 'y': 'Pressure [hPa]', 'ylim': (-50,100), 'ylimdes': (0,30) },

        34 : {'x': 'Dew Point Depr. [K] '    , 'y': 'Pressure [hPa]', 'ylim': (-5,60), 'ylimdes': (0,30) },
        36 : {'x': 'Dew Point Temp. [K]  '    , 'y': 'Pressure [hPa]', 'ylim': (180,340), 'ylimdes': (0,30) },
        38 : {'x': 'Relative Humidity '       , 'y': 'Pressure [hPa]', 'ylim': (-1.5,1.5), 'ylimdes': (0,1.5) },
        39 : {'x': 'Specific Humidity '      , 'y': 'Pressure [hPa]', 'ylim': (0,0.01), 'ylimdes': (0,0.001) },
          }
"""





def get_CUON_wind(file, merged=True, plev=True, hour='' ):
    """ Dedicated module to analyze wind data from merged file """
    
    f = h5.File(file, 'r')
    
    ind = f["recordindices"]['104']
    imin_uwind, imax_uwind = ind[0], ind[-1]

    ind = f["recordindices"]['105']
    imin_vwind, imax_vwind = ind[0], ind[-1]
    
    ind = f["recordindices"]['106']
    imin_dir, imax_dir = ind[0], ind[-1]
    
    ind = f["recordindices"]['107']
    imin_speed, imax_speed = ind[0], ind[-1]
    
    dt = f["observations_table"]["date_time"][imin_uwind:imax_uwind]
    time_c = dt/(365.25*60*60*24) + 1900
    dt = pd.to_datetime(dt,  unit='s', origin=pd.Timestamp('1900-01-01'))
    
    #calc_speed =  np.sqrt (  np.square( f["observations_table"]["observation_value"][imin_uwind:imax_uwind]) + np.square(f["observations_table"]["observation_value"][imin_vwind:imax_vwind]) ) 
    #calc_dir = np.arctan2( f["observations_table"]["observation_value"][imin_uwind:imax_uwind] , f["observations_table"]["observation_value"][imin_vwind:imax_vwind]  )
    #calc_dep_dir = np.arctan2(f["era5fb"]["an_depar@body"][imin_uwind:imax_uwind] , f["observations_table"]["observation_value"][imin_vwind:imax_vwind]  )
    
    red_df = pd.DataFrame(  { 'z'     : f["observations_table"]["z_coordinate"][imin_uwind:imax_uwind] ,
                              
                                                'fg_uwind'    : f["era5fb"]["fg_depar@body"][imin_uwind:imax_uwind],
                                                'an_uwind'   : f["era5fb"]["an_depar@body"][imin_uwind:imax_uwind],
                                                'obs_uwind' : f["observations_table"]["observation_value"][imin_uwind:imax_uwind],
                                                
                                                'fg_vwind'    : f["era5fb"]["fg_depar@body"][imin_vwind:imax_vwind],
                                                'an_vwind'   : f["era5fb"]["an_depar@body"][imin_vwind:imax_vwind],
                                                'obs_vwind' : f["observations_table"]["observation_value"][imin_vwind:imax_vwind],

                                                'fg_dir'    : f["era5fb"]["fg_depar@body"][imin_dir:imax_dir],
                                                'an_dir'   : f["era5fb"]["an_depar@body"][imin_dir:imax_dir],
                                                'obs_dir' : f["observations_table"]["observation_value"][imin_dir:imax_dir],
                                                
                                                'fg_speed'    : f["era5fb"]["fg_depar@body"][imin_speed:imax_speed],
                                                'an_speed'   : f["era5fb"]["an_depar@body"][imin_speed:imax_speed],
                                                'obs_speed' : f["observations_table"]["observation_value"][imin_speed:imax_speed],                                                
                                                
                                            }
                                        )

    red_df["month"] = dt.month
    red_df["year"] = dt.year        
    red_df["hour"] = dt.hour        
        
    red_df["time"] = time_c         
    
    if hour<1:
        print("*** Extracting hour 00 ")        
        red_df = red_df.loc[ (red_df['hour'] > 21) | (red_df['hour'] < 3 ) ]
    else:
        print("*** Extracting hour 12 ")        
        red_df = red_df.loc[ (red_df['hour'] > 9) & (red_df['hour'] < 15 ) ]        
        
    if plev:
        print("*** Extracting standard pressure level")
        red_df = red_df.loc[ (red_df['z'] == 1000 )  |  ( red_df['z'] == 2000 ) | ( red_df['z'] == 3000 ) | ( red_df['z'] == 5000 ) | ( red_df['z'] == 7000 ) | ( red_df['z'] == 10000 ) |
                                         (red_df['z'] == 15000 )  |
                                         (red_df['z'] == 20000 )  |  ( red_df['z'] == 25000 ) | ( red_df['z'] == 30000 ) |
                                         (red_df['z'] == 40000 )  |
                                         (red_df['z']== 50000 ) | ( red_df['z'] == 70000 ) | ( red_df['z'] == 85000 ) |
                                         (red_df['z'] == 92500 )  |  ( red_df['z']== 100000 )  ]

    # [1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 85000, 92500, 100000]




    return red_df




def get_CUON_data(file, var="85", merged=True, std_plev=True, hour='', rolling_mean=True ):
    f = h5.File(file, 'r')
    
    ind = f["recordindices"][var]
    imin, imax = ind[0], ind[-1]

    z = f["observations_table"]["z_coordinate"][imin:imax] 
    
    dt = f["observations_table"]["date_time"][imin:imax]
    time_c = dt/(365.25*60*60*24) + 1900
    dt = pd.to_datetime(dt,  unit='s', origin=pd.Timestamp('1900-01-01'))

    obs = f["observations_table"]["observation_value"][imin:imax]
    bias = f["era5fb"]["biascorr@body"][imin:imax]
    fg_dep_adj = f["era5fb"]["fg_depar@body"][imin:imax]

    if var == "85": # bias is available only for temperature 
         
        obs_adj = obs - bias
        fg_dep = fg_dep_adj + bias # TODO


        """
        for p in [1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 85000, 92500, 100000]:
            
            ind_z = np.where ( z == p )[0]
            
            #plt.plot(time_c[ind_z], fg_dep_adj[ind_z]        , label = 'fg adj')
            #plt.plot(time_c[ind_z], fg_dep[ind_z]  , label = 'fg')
            
            plt.plot(time_c[ind_z], fg_dep[ind_z]  - fg_dep_adj[ind_z]        , label = 'bias calc')
            plt.plot(time_c[ind_z], bias[ind_z]       , label = 'bias')
            
            
            plt.ylim(-5,5)
            plt.legend()
            plt.title('Pressure ' + str(p ) )
            #plt.show()
            os.system('mkdir Plots/bias/check_bias/')
            plt.savefig('Plots/bias/check_bias/check_bias_' + str(p) + '.png')
            plt.close()
        """
        
        print(0)
    else: # no bias is available 
        obs_adj = obs
        fg_dep = fg_dep_adj
        
        
    red_df = pd.DataFrame(  { 'z'            : z,
                              
                                                'fg_dep_adj'    : fg_dep_adj,
                                                'fg_dep'           : fg_dep ,                                                 
                                                'bias'       : bias ,
                                                'obs'        : obs,
                                                'obs_adj'        : obs_adj,
                                                'bias_calc'     : fg_dep - fg_dep_adj,
                                            }
                                        )
    
    red_df["month"] = dt.month
    red_df["year"] = dt.year        
    red_df["hour"] = dt.hour        
        
    red_df["time"] = time_c         
    
    if hour<1:
        print("*** Extracting hour 00 ")        
        red_df = red_df.loc[ (red_df['hour'] > 21) | (red_df['hour'] < 3 ) ]
    else:
        print("*** Extracting hour 12 ")        
        red_df = red_df.loc[ (red_df['hour'] > 9) & (red_df['hour'] < 15 ) ]        
        
    if std_plev:
        print("*** Extracting standard pressure level")
        red_df = red_df.loc[ (red_df['z'] == 1000 )  |  
                                         (red_df['z'] == 2000 ) |
                                         (red_df['z'] == 3000 ) |
                                         (red_df['z'] == 5000 ) |
                                         (red_df['z'] == 7000 ) |
                                         (red_df['z'] == 10000 ) |
                                         (red_df['z'] == 15000 )  |
                                         (red_df['z'] == 20000 )  | 
                                         (red_df['z'] == 25000 ) | 
                                         (red_df['z'] == 30000 ) |
                                         (red_df['z'] == 40000 )  |
                                         (red_df['z']== 50000 ) | 
                                         (red_df['z'] == 70000 ) | 
                                         (red_df['z'] == 85000 ) |
                                         (red_df['z'] == 92500 )  |  
                                         (red_df['z']== 100000 )  ]

    # [1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 85000, 92500, 100000]

    # Applying rolling mean 
    # first I have to split ds by pressure level, apply the rolling mean, then re-add the ds 
    dfs = []
    if rolling_mean:
    
        for p in  [1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 85000, 92500, 100000]:
            d = red_df.loc[red_df['z'] == p]
            
            for k in ['obs_adj' , 'obs', 'fg_dep_adj' , "fg_dep" ]:
                dep_series = pd.Series( d[k]  )
                d[k] = dep_series.rolling(60, center=True).mean()          
            
            dfs.append(d)
        
        red_df = pd.concat( dfs )
        
    return red_df






def make_hovmoller(time, plevel, value, variable="Temperature", hour=0,
                                     what='' ,
                                     station='Vienna 0-20001-0-11035' , 
                                     z_label = "Obs-Bkg") :
    
    print("Making the plot for ", station , ' ' , variable )
    fs = 15
    
    # setting bar color
    if "Temp" in variable:
        col = 'coolwarm'
        levels = np.linspace(-3,8,12)
        nums = [ 1*i for i in range(-3,9) ]  # must be = number of levels +2
        
        if 'Bethel' in station:
            col = 'seismic'
            levels = np.linspace(-2,2,100)
            nums = [ 0.5*i for i in range(-4,4) ]  # must be = number of levels +2
            
    elif "Hum" in variable:
        col = 'cividis'
        levels = np.linspace(-0.5,0.5,11)
        nums = [ 0.1*i for i in range(-12,12) ]  
        
    elif "speed" in variable:
        col = 'GnBu_r'
        levels = np.linspace(-5,10,16)
        nums = [ 1*i for i in range(-5,12) ]  
        
    elif "direction" in variable:
        col = 'afmhot'
        levels = np.linspace(-20,70,10)
        nums = [ 10*i for i in range(-3,8) ]  
        
    fig,ax = plt.subplots(figsize = (10,5))

    contours = ax.tricontourf(time, plevel/100, value, levels = levels, cmap = col, extend = 'min')
    #contours = ax.tricontourf(time, plevel, value, cmap = col, extend = 'min')

    if hour ==0:
        h_n = '00'
    else:
        h_n = '12'
    plt.title(variable + ' - Station ' + station + " " + h_n + "Z" , fontsize = fs, y = 1.03)
    plt.ylabel('Pressure [hPa] ', fontsize = fs)
    
    plt.gca().invert_yaxis()
    
    # Color bar properties
    cbar = fig.colorbar(contours)
    cbar.set_ticks(nums)
    cbar.set_label(z_label, fontsize = fs )
    
    # Set axis labels
    plevels_toshow = [10, 50, 70, 100, 150, 300, 500, 700, 850, 925, 1000]
    ax.set_yticks( plevels_toshow ) 
    ax.set_yticklabels( plevels_toshow )
    
    if not os.path.isdir("Plots/bias/"):
        os.system("mkdir Plots/bias/")
        
    stat_name = station.replace(" ","-").replace("/","").replace("(","").replace(")","")
    
    if not os.path.isdir("Plots/bias/" + stat_name  ):
            os.system("mkdir Plots/bias/" + stat_name)
            
    plt.savefig("Plots/bias/"+ stat_name + '/' + what + '_' + variable.replace(" ","-").replace("/","") + "_" +  station.replace("/","")  + '_' + str(hour) + '_' + z_label + '.png' , dpi=200)
    plt.close()
    
#data = data[10000:50000]



stat_dic = {
                   #"Vienna (0-20001-0-11035)": '/mnt/users/scratch/leo/scratch/converted_v7/0-20001-0-11035_CEUAS_merged_v1.nc'         , 
                   #"Aktobe(0-20000-0-35229)" : '/mnt/users/scratch/leo/scratch/converted_v7/0-20000-0-35229_CEUAS_merged_v1.nc'         ,
                   #"Payerne (0-20000-0-06610)" : '/mnt/users/scratch/leo/scratch/converted_v7/0-20000-0-06610_CEUAS_merged_v1.nc'       ,
                   #"Lindenberg (0-20001-0-10393)" : '/mnt/users/scratch/leo/scratch/converted_v7/0-20001-0-10393_CEUAS_merged_v1.nc'  ,
                   "Bethel (0-20000-0-70219)" : '/mnt/users/scratch/leo/scratch/converted_v7/0-20000-0-70219_CEUAS_merged_v1.nc'          ,

 }


#0-20000-0-70219_CEUAS_merged_v1.nc

def make_grid_adjustment(df, station = '', hour ='0', all_levels= False):
    """ Plot the grid of temperature adjustments for a station on all pressure level """
    if all_levels:
        std =      [1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 85000, 92500, 100000]
        lab = 'all'
        fig,axs = plt.subplots(len(std), figsize = (10,20))
        
    else:
        std =      [1000, 5000, 10000, 20000, 50000, 85000, 100000] # 10
        lab = 'small'
        fig,axs = plt.subplots(len(std), figsize = (10,15))


    if hour ==0:
        hour = '00'
    else:
        hour = '12'
        
    for p,num in zip(std, range(len(std)) ):

        r = df.loc[df['z'] == p]
        time = r['time']
        bias = r['bias']
        
        axs[num].plot(time, bias, label = str(int(p/100)) + ' [hPa]')
        axs[num].legend()
        axs[num].set_ylim([-3,3])
        
        if num < len(std)-1:
            axs[num].set_xticks([])
            
        axs[num].set_ylabel("Bias [K]")
        
    axs[0].set_title("Temperature Bias Adjustment " + station + ' - '  + str(hour) + 'Z', y =1.02)
    
    stat_name = station.replace(" ","-").replace("/","").replace("(","").replace(")","")
    
    plt.savefig("Plots/bias/"+ stat_name + '/bias-grid_' +  station.replace("/","")  + '_' + str(hour) + '_' + lab +  '.png' , dpi=200)

    return 0
    
    
var_dic = {
                   "Temperature [K]"       : '85' ,            
                   'Relative Humidity'      : '38',
                   'Dew Point Temp. [K]  ': '34',
    
                  "Wind speed [m/s]"     : '107', 
                  "Wind direction"          : '106',              
                  "Wind speed u-comp [m/s]"     : '104', 
                  "Wind speed v-comp [m/s]"     : '105', 
                  
                  }


"""
var_dic = {                
                  "Wind direction"          : '106', 
                  "Wind speed [m/s]"     : '107', 

                  }
"""



var_dic = {                
    "Temperature [K]"       : '85' ,            


                  }




def get_direction_dep(red_df):
    """ Calculates the background u and v component, and the direction of the wind.
          Extracts then the obs-bkg using the observed direction from u and v. 
          
          fg = obs - bkg
          bkg = obs - fg 
          
          # The meteorological convention for winds is that U component is positive for a west to east flow (eastward wind) 
          # and the V component is positive for south to north flow (northward wind).
          
          #https://confluence.ecmwf.int/pages/viewpage.action?pageId=133262398
          
          """
    
    u, v = red_df["obs_uwind"] ,  red_df["obs_vwind"]
    u_dep, v_dep = red_df["fg_uwind"] ,  red_df["fg_vwind"]
    
    u_bkg, v_bkg = u-u_dep , v-v_dep
    
    obs_dir = np.arctan2(v, u) * 180 / np.pi
    bkg_dir = np.arctan2(v_bkg, u_bkg) * 180 / np.pi
    
    dep_dir = obs_dir - bkg_dir
    
    red_df["calc_dir_dep"] = dep_dir
    red_df["obs_dir"] = obs_dir 
    red_df["bkg_dir"] = bkg_dir 
    
    dep_series = pd.Series( dep_dir  )
    red_df["rolled_mean"] = dep_series.rolling(30, center=True).mean()  
    
    return red_df

def make_plot_direction(red_df, plevel = 20000, station=''):
    
    red_df = red_df.loc[ data_wind["z"] == 20000 ]
 
    fig, axs = plt.subplots(3, figsize=(15,10))
    
    axs[0].set_title("Observations, Background and Departures of wind direction - Station " + station + ' - ' + str(plevel) + ' [hPa]', y =1.03)
    
    obs_dir = red_df["obs_dir"] 
    bkg_dir = red_df["bkg_dir"] 
    dep_dir = red_df["calc_dir_dep"]
    
    axs[0].plot( red_df['time'], obs_dir, label = 'Obs. Direction', color = 'blue')
    axs[0].plot( red_df['time'], bkg_dir, label = 'Bkg. Direction ', color = 'orange')
    axs[0].legend()
    
    axs[1].plot( red_df['time'], dep_dir, label = 'Dep. Obs-Bkg', color = 'red')
    axs[1].legend()
    
    axs[2].plot( red_df['time'], red_df['rolled_mean'], label = 'Obs-Bkg Rolling Mean - 30 Days ', color = 'lime')
    axs[2].legend()    
    
    for i in [0,1,2]:
              axs[i].legend(loc = 'upper left')   
              axs[i].grid(color='lightgray', ls=":")
              axs[i].set_ylabel("Wind direction [degree]")
              
    axs[2].set_ylim([-50,100])
    plt.savefig("Plots/bias/" + station + "_direction_test_departure.png" , dpi = 200)
    
    plt.close()        
    
    # Calculate rolling mean

    
def check_some_levels(data, station=''):
    
    """ Make a grid plot with departures and bias adjustments """
    
    levels = [2000,20000,50000,70000,85000,100000]
    
    fig, ax = plt.subplots(len(levels)*2, figsize=(10,20))

    data["bias_calc"] = data['fg_dep'] - data['fg_dep_adj']
    
    f = data[ ['bias','bias_calc','fg_dep','fg_dep_adj','z','time'] ]
    for p,n in zip(levels, [2*i for i in range(len(levels)) ] ): 

        d = f.loc[f['z'] == p ] 
        
        #ax[n].plot(d['time']      , d['fg_dep']                            , label = 'Obs-Bkg ' + str(int(p/100)) + ' [hPa]'  , color = 'lime'     )
        #ax[n].plot(d['time']      , d['fg_dep_adj']                      , label = 'Obs-Bkg ' + str(int(p/100)) + ' [hPa]'  , color = 'orange')
        #ax[n+1].plot(d['time'] , d['bias']                                 , label = 'Bias'                                                     , color = 'black'  )   
        #ax[n+1].plot(d['time'] , d['fg_dep'] - d['fg_dep_adj']  , label = 'Bias diff '                                              , color = 'blue'     )
    
        ax[n].plot(d['time']      , d['fg_dep'] - d['fg_dep_adj']                      , label = 'Bias diff  ' + str(int(p/100)) + ' [hPa]'  , color = 'orange')
        ax[n+1].plot(d['time'] , d['bias']                                 , label = 'Bias'                                                     , color = 'black'  )       
    
        ax[n+1].set_ylim([-1,3])
        
        ax[n].set_xticks([])
        if n <= max([2*i for i in range(len(levels)) ]):
            ax[n].set_xticks([])

        ax[n].legend(loc= 'upper left')
        ax[n+1].legend(loc= 'upper left')

    """
    for i in range(11):
        ax[i].set_xticks([])
        ax[i].legend(loc= 'upper left')
        
    for i in [0,2,4,6,8,10]:
        
        ax[i].set_ylim([-1,3])
        
    ax[7].legend(loc= 'upper left')
    """
    
    plt.savefig("Plots/bias/check_bias/" + station + "_test_levels.png" , dpi = 200)
    print(0)
    
     
for  station in stat_dic.keys():
    file = stat_dic[station]
    
    for v in var_dic.keys():
    
        var = var_dic[v]
        
        for h in [0,12]:
            
            # special case: wind direction 
            if var == '106':
                
                data_wind = get_CUON_wind('/mnt/users/scratch/leo/scratch/converted_v7/0-20000-0-35229_CEUAS_merged_v1.nc', plev=True, hour = h)
                
                data_wind = get_direction_dep(data_wind)
                if data_wind.empty: 
                    continue
                data_wind = data_wind.dropna()
                
                # making separated direction plots 
                make_plot_direction(data_wind, plevel = 20000, station=station)

                # extracting data for Hovmoller plot
                time = data_wind["time"]
                plevel = data_wind["z"]
                value = data_wind["rolled_mean"]
            
                d = make_hovmoller(time, plevel, value, variable="Wind direction", station=station, hour = h)
            
            else:
                rm = False
                data = get_CUON_data(file, var=var, std_plev=True, hour = h , rolling_mean=rm)
                
                data = data.dropna()
                
                if data.empty:
                    print('----- Empty Dataframe :-( ')
                    continue 
                
                time = data["time"]
                plevel = data["z"]
                value = data["fg_dep"]
                bias = data["bias"]
                
                check_some_levels(data, station=station)
                
                make_hovmoller(time, plevel, bias, variable=v, hour=h,  station=station,  what = 'Bias'   , z_label='Bias')
                make_hovmoller(time, plevel, value, variable=v, hour=h,  station=station,  what = 'Temp' , z_label='Obs-Bkg (No Bias Adj.)')

                if var == '85': 
                    value = data["fg_dep_adj"]
                    make_hovmoller(time, plevel, value, variable=v, hour=h,  station=station, z_label='Obs-Bkg (Bias Adjusted)')
                    
                    make_grid_adjustment(data,  station=station, hour=h)
                    make_grid_adjustment(data,  station=station, hour=h, all_levels= True)

print(0)