import h5py
import matplotlib.pylab as plt
import numpy
import glob
import os,sys
import matplotlib.gridspec as gridspec
import xarray as xr 
import pandas as pd
import numpy as np
from tqdm import tqdm

import matplotlib
matplotlib.use('Agg')

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)


fn = '/raid8/srvx1/federico/GitHub/CEUAS_master_FEBRUARY2021/CEUAS/CEUAS/public/merge/PROVA/0-20000-0-94120_CEUAS_merged_v0_twohours.nc'
f = h5py.File(fn,'r')
nm = fn.split('/')[-1]
SIZE = 1.5
#print(list(f['observations_table'].keys()))
fo=f['observations_table']


os.system('mkdir Plots')
    
def plot_comparison(p, zoom= ''):
        
        fig = plt.figure(constrained_layout = True, figsize=(10,10))
        
        gs = gridspec.GridSpec(3, 1, figure = fig) # (number of rows, number of columns )
        ax = fig.add_subplot(gs[0,:])
        print(p)
        
        mask=numpy.logical_and(fo['z_coordinate'][:]==p, fo['observed_variable'][:]==85)
        #print(fn,numpy.sum(numpy.logical_and(fo['date_time'][mask]/365.25/86400>=60,fo['date_time'][mask]/365.25/86400<61)))
        
        
        time = [1900 + f for f in fo['date_time'][mask]/365.25/86400 ] 
        plt.scatter(time,fo['observation_value'][mask], label = nm, color = 'blue' , s = SIZE)
        ax.legend()
        if zoom:
                
                ax.set_xlim(1955,1965)
        else:
                
                ax.set_xlim(1945,1980)
                
        #ax.set_ylim(180,300)
        
        ax.set_title(' p=' + str(p) + ' [Pa]' , y = 1.02 , fontsize = 14)    
        
        fns=glob.glob('/raid60/scratch/leo/scratch/era5/odbs/2/era5_2/*94120*.nc')    
        for fn,i in zip(fns,[1,2]):
            ax = fig.add_subplot(gs[i,:])
            
            n = fn.split('/')[-1]
            fg=h5py.File(fn,'r')
    
            foo=fg['observations_table']
            mask=numpy.logical_and(foo['z_coordinate'][:]==p, foo['observed_variable'][:]==85)
            #print(fn,numpy.sum(numpy.logical_and(fo['date_time'][mask]/365.25/86400>=60,fo['date_time'][mask]/365.25/86400<61)))
            time = [1900 + f for f in foo['date_time'][mask]/365.25/86400 ]    
            plt.scatter(time,foo['observation_value'][mask] , label = n, s = SIZE)
            ax.legend()
            
            if zoom:
                ax.set_xlim(1955,1965)
            else:
                ax.set_xlim(1945,1980)
            #ax.set_ylim(180,300)
            
        if zoom:
                
                plt.savefig('Plots/era5_merged_zoom_' + str(p) + '.png', dpi = 200)    
        else:
                
                plt.savefig('Plots/era5_merged_' + str(p) + '.png', dpi = 200)    
                
        plt.close()

""""
for p in [1000, 5000, 90000, 10000, 50000, 90000, 92500, 100000]:
        plot_comparison(p)
        plot_comparison(p, zoom= True)
"""

def plot_records(zoom = False):
        
        fig = plt.figure(constrained_layout = True, figsize=(10,10))
        
        gs = gridspec.GridSpec(3, 1, figure = fig) # (number of rows, number of columns )
        ax = fig.add_subplot(gs[0,:])
        
        records_1, plev_1, time_1 = [],[],[]
        
        ind_min = f['recordindex']
        ind_max = [ ind_min[i+1] for i in range (0, len(ind_min)-2)  ]
        ts = f['recordtimestamp']
        zcoo = fo['z_coordinate']
        for d, imin, imax in zip( tqdm(ts) , ind_min, ind_max ) :
                if (d/365.25/86400 + 1900)< 1950 or ( d/365.25/86400 + 1900) > 1980:
                        continue
                
                records_1.append(  imax-imin )
                time_1.append( 1900 +d/365.25/86400 )
                plev_1.append( len (list (set(zcoo[imin:imax]) ) ) ) 
                
        print('Plotting first plot ')
        plt.scatter(time_1, records_1 , label = 'Total Records ', color = 'blue' , s = SIZE)
        plt.scatter(time_1, plev_1      , label = 'Pressure Levels', color = 'cyan' , s = SIZE)
        plt.scatter([-9999], [-9999]   , label = 'Merged File ' , color = 'cyan' , s = SIZE)
        ax.set_ylim(10,300)
                
        ax.legend()
        if zoom:
                ax.set_xlim(1955,1965)
        else:
                ax.set_xlim(1945,1980)
                
        ax.grid(ls =":" , color = "lightgray")
        ax.set_title('Records and Pressure Levels Counts ' , y = 1.02 , fontsize = 14)   
        fns=glob.glob('/raid60/scratch/leo/scratch/era5/odbs/2/era5_2/*94120*.nc')  
        #f2 = b'0-20000-0-94120_era5_2_harvested_era5.conv._1:94120.gz.nc'
        #f1 = b'0-20000-0-94120_era5_2_harvested_era5.conv._94120.gz.nc'    
        
        #ff1 = '/raid60/scratch/federico/HARVESTED_JAN2021/era5_2/' + f1.decode('utf-8')
        #ff2 = '/raid60/scratch/federico/HARVESTED_JAN2021/era5_2/' + f2.decode('utf-8')
        
        fns = [ff1,ff2] 
        
        for fn,i in zip(fns,[1,2]):
                
                Records, Plev, Time = [],[], []
                
                ax = fig.add_subplot(gs[i,:])
                n = fn.split('/')[-1]
                
                fg=h5py.File(fn,'r')
                foo = fg['observations_table']
                ts   = fg['recordtimestamp']                
                ind_min = fg['recordindex']
                
                ind_max = [ ind_min[i+1] for i in range (0, len(ind_min)-2)  ]
                zcoo = foo['z_coordinate']
                
                for d, imin, imax in zip( tqdm(ts) , ind_min, ind_max ) :
                        if (d/365.25/86400 + 1900)< 1950 or ( d/365.25/86400 + 1900) > 1980:
                                continue
                        
                        Records.append(  imax-imin )
                        Time.append( 1900 +d/365.25/86400 )
                        Plev.append( len (list (set(zcoo[imin:imax]) ) ) ) 
    
                        
                print('Plotting ', str(i) ,  ' plot ')
                        
                plt.scatter(Time, Records , label = 'Total Records ', color = 'red' , s = SIZE)
                plt.scatter(Time, Plev   , label = 'Pressure Levels', color = 'orange' , s = SIZE)
                plt.scatter([-9999], [-9999]   , label = n , color = 'cyan' , s = SIZE)    
            
                ax.legend()
                
                if zoom:
                        ax.set_xlim(1955,1965)
                else:
                        ax.set_xlim(1945,1980)
                ax.set_ylim(10,300)
                ax.grid(ls =":" , color = "lightgray")
            
        if zoom:
                plt.savefig('Plots/era5_merged_records_zoom' + '.png', dpi = 200)    
        else:
                plt.savefig('Plots/era5_merged_records' + '.png', dpi = 200)    
                
        plt.close()
        
        

        
        
        
def plot_time_distr():
    
    """ I analyze the time distirbution for the two different sources in the merged file """
    fn='/raid60/scratch/federico/MERGING_JAN2021_FIXED/0-20000-0-94120_CEUAS_merged_v0.nc'
    fn ='/raid8/srvx1/federico/GitHub/CEUAS_master_FEBRUARY2021/CEUAS/CEUAS/public/merge/PROVA/0-20000-0-94120_CEUAS_merged_v0_twohours.nc'
    
    time = xr.open_dataset(fn , engine = 'h5netcdf' , decode_times = True )['recordtimestamp'].values
    s = xr.open_dataset(fn,engine = 'h5netcdf' , group = 'source_configuration', decode_times = True )
    source = s['source_file'].values
    
    for s in list(set(source)):
        print(s)
        
    f2 = b'0-20000-0-94120_era5_2_harvested_era5.conv._1:94120.gz.nc'
    f1 = b'0-20000-0-94120_era5_2_harvested_era5.conv._94120.gz.nc'
    
    dic = {'time': time , 'source':source}
    df = pd.DataFrame.from_dict(dic)    
    
    t1 = df.loc [df['source'] == f1]['time']
    t2 = df.loc [df['source'] == f2]['time']
    
    y1 = np.empty( len(t1) )
    y1.fill( 1 )
    
    y2 = np.empty( len(t2) )
    y2.fill( 2 )
    
    plt.scatter(t1, y1, label = f1.decode('utf-8').split('harvested_')[1]  + '[' + str(len(t1))+ '] Merged' , color = 'orange')
    plt.scatter(t2, y2, label = f2.decode('utf-8').split('harvested_')[1]  + '[' + str(len(t2))+ '] Merged', color = 'lime')
    
    f2 = b'0-20000-0-94120_era5_2_harvested_era5.conv._1:94120.gz.nc'
    f1 = b'0-20000-0-94120_era5_2_harvested_era5.conv._94120.gz.nc'    
    
    ff1 = '/raid60/scratch/federico/HARVESTED_JAN2021/era5_2/' + f1.decode('utf-8')
    ff2 = '/raid60/scratch/federico/HARVESTED_JAN2021/era5_2/' + f2.decode('utf-8')
    
    time_ff1 = xr.open_dataset(ff1 , engine = 'h5netcdf' , decode_times = True )['recordtimestamp'].values
    time_ff2 = xr.open_dataset(ff2 , engine = 'h5netcdf' , decode_times = True )['recordtimestamp'].values
    yy1 = np.empty( len(time_ff1) )
    yy1.fill(3)
    yy2 = np.empty( len(time_ff2) )
    yy2.fill(4)
    plt.scatter(time_ff1, yy1, label = f1.decode('utf-8').split('harvested_')[1] + '[' + str(len(time_ff1))+ '] harvested', color = 'red', s = SIZE)
    plt.scatter(time_ff2, yy2, label = f2.decode('utf-8').split('harvested_')[1]  + '[' + str(len(time_ff2))+ '] harvested' , color = 'green', s = SIZE)
    
    
    plt.legend()
    plt.savefig('Plots/time_intervals.png')
    
    
    
    

def makedf(f):
    s = xr.open_dataset(f,engine = 'h5netcdf' , group = 'observations_table', decode_times = True )
    time = s['date_time'].values
    plev = s['z_coordinate'].values
    ztype = s['z_coordinate_type'].values
    var = s['observed_variable'].values
    val = s['observation_value'].values
    
    dic = {'time': time , 'plev':plev, 'ztype':ztype , 'var':var , 'val': val }
    df = pd.DataFrame.from_dict(dic)          
    return df

"""
def plot_fast():

    print('Extracted DF')
    

    
    
    def reduce(df, p):
        return df.loc [ (df['plev'] == p) & (df['var'] == 85) ]
    
    merged = fn

    ff1 = '/raid60/scratch/federico/HARVESTED_JAN2021/era5_2/0-20000-0-94120_era5_2_harvested_era5.conv._1:94120.gz.nc'
    ff2 = '/raid60/scratch/federico/HARVESTED_JAN2021/era5_2/0-20000-0-94120_era5_2_harvested_era5.conv._94120.gz.nc'
    
    m   = makedf(merged)
    df1 = makedf(ff1)
    df2 = makedf(ff2)
    

    
    def add_subplot(df, fig, ax, i = 0, l = 'ciao', color = 'black'):
        
        ax = fig.add_subplot(gs[i,:])
        ax.set_ylim(160,330)
        plt.scatter(df['time'], df['val'], label = l , color = color, s = SIZE)        
        ax.legend()
        ax.set_xlim(np.datetime64('1940-01-01'),np.datetime64('1980-01-01'))
        
    for p in [1000, 5000, 10000, 50000, 92500, 100000]:
        fig = plt.figure(constrained_layout = True, figsize=(10,8))
        
        gs = gridspec.GridSpec(3, 1, figure = fig) # (number of rows, number of columns )        
        print(p)
        m_red = reduce(m, p)
        df1_red = reduce(df1, p)
        df2_red = reduce(df2, p )
        
        ax = fig.add_subplot(gs[0,:])


        add_subplot(m_red, fig, ax, i = 0, l='Merged', color = 'blue')
        add_subplot(df1_red,fig,  ax, i = 1, l='era5.conv._1:94120.gz.nc', color = 'lime')
        add_subplot(df2_red, fig, ax, i = 2, l='era5.conv._94120.gz.nc', color = 'red')
    
        plt.savefig('Plots/fast_pandas_' + str(p) + '.png')
        plt.close()
    


        
        
        
        
#plot_fast()
"""






merged = '/raid8/srvx1/federico/GitHub/CEUAS_master_FEBRUARY2021/CEUAS/CEUAS/public/merge/PROVA/0-20000-0-94120_CEUAS_merged_v0_twohours.nc'

ff1 = '/raid60/scratch/federico/HARVESTED_JAN2021/era5_2/0-20000-0-94120_era5_2_harvested_era5.conv._1:94120.gz.nc'
ff2 = '/raid60/scratch/federico/HARVESTED_JAN2021/era5_2/0-20000-0-94120_era5_2_harvested_era5.conv._94120.gz.nc'

a = '/raid60/scratch/leo/scratch/era5/odbs/2/era5_2/'

# I checked and the two harvested files seem identical, so OK
g2 = a + '0-20000-0-94120_era5_2_harvested_era5.conv._94120.gz.nc'
g1 = a + '0-20000-0-94120_era5_2_harvested_era5.conv._1:94120.gz.nc'




if __name__ == '__main__':

        #plot_time_distr()

        #plot_records()  
        #plot_records(zoom = True )  
        df1 = makedf(ff1)
        df2 = makedf(g1)
        
        df_merged = makedf(merged)
        
        0















