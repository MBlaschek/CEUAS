""" Analyse and plot the u,v wind components time series from netCDF files.

author: Ambrogi Federico

This files gives retrives basic information and creates the plots
for the u,v wind compontents from the nectCDF files.

The netCDF files were converted from the odb format to net CDF using the 
version of script 'readodbstationfiles.py'.
Databases included: ['1', 3188','1759','1761'] in the root dir '/raid60/scratch/leo/scratch/era5/odbs/'

netCDF file root directory:  '/raid8/srvx1/federico/odb_netCDF/$FSCRATCH/ei6'

"""
import netCDF4
import matplotlib.pylab as plt
import os,sys
import os.path
import matplotlib.gridspec as gridspec
from pylab import rcParams


comb_file = 'station10393_test/ERA5_1_10393_calc_sh_rh_dp.nc'
t_file = 'station10393_test/ERA5_1_10393_t.nc'
dp_file = 'station10393_test/ERA5_1_10393_dp.nc'
rh_file = 'station10393_test/ERA5_1_10393_rh.nc'
c = netCDF4.Dataset(comb_file) 
#var = f.variables
#print(var)

plt.plot()


#plt.tight_layout()
#ax1.xaxis.set_major_formatter(plt.NullFormatter())

def plot_comb(p,h):
     
     dic_prop = {  'temp':{'xlab':'Year', 'ylab':'Temperature [K]'         ,'leg':'Temperature',       'ax':[1990,2019,  200 ,300] ,'c':'orange'       } ,
                   'rh':{'xlab':'Year', 'ylab':'Relative Humidity ' ,'leg':'Relative Humidity',   'ax':[1990,2019,  0   , 1 ] ,'c':'lightgray'    } ,                       
                   'dp':{'xlab':'Year', 'ylab':'Dew Point [K]'           ,'leg':'Upper Air Dew Point', 'ax':[1990,2019,  150 ,300] ,'c':'mediumpurple'} ,
                           }     
     fnt_size = 10
     rcParams['figure.figsize']= 10, 13
     
     gs = gridspec.GridSpec(5,1)
     
     datum = [ 1900 + d for d in c.variables['datum'][0,:]/365.25 ] 
     
     ct = c.variables['t_comb'][h,p,:]
     cdp   = c.variables['dew_point_temperature_comb'][h,p,:]
     crh   = c.variables['relative_humidity_comb'][h,p,:]
     cfdp  = c.variables['dp_flag'][h,p,:]
     cfrh  = c.variables['rh_flag'][h,p,:]     
     
     t =  netCDF4.Dataset(t_file).variables['temperatures'][h,p,:]
     dp =  netCDF4.Dataset(dp_file).variables['dp'][h,p,:]
     rh =  netCDF4.Dataset(rh_file).variables['rh'][h,p,:]

     # temp *******************************************************************************
     ax0 = plt.subplot(gs[0])
     plt.title('Station 10393 - Combination test for h=' + str(h) + ' and p=' + str(p))
     
     plt.plot(datum    , t  , label = 'Temperature'       , color = 'red'   , linestyle = '--' , lw = 1.2 )     
     plt.scatter(datum , ct , label = 'Temperature Comb.' , color = 'black' , s = 4 )     

     plt.axis(dic_prop['temp']['ax'])
     plt.xlabel(dic_prop['temp']['xlab'], fontsize = fnt_size)
     plt.ylabel(dic_prop['temp']['ylab'], fontsize = fnt_size)
     plt.legend(loc = 'lower left', fontsize = 10)
     
     # rh *******************************************************************************
   
     ax0 = plt.subplot(gs[1])
     plt.plot(datum    , rh  , label = 'Relative Hum.', color = 'gold'   , linestyle = '--' , lw = 1.2 )   
     plt.scatter(datum , crh , label = 'Relative Hum. Comb.' , color = 'black' , s = 1 )     

     plt.axis(dic_prop['rh']['ax'])
     plt.xlabel(dic_prop['rh']['xlab'], fontsize = fnt_size)
     plt.ylabel(dic_prop['rh']['ylab'], fontsize = fnt_size)
     plt.legend(loc = 'lower left', fontsize = 10)   

     ax0 = plt.subplot(gs[2])
     obs, meas , dobs, dmeas = [],[],[],[]
     for v,f,d in zip(crh,cfrh,datum):
          if f == 1 or f == '1':
               obs.append(f)
               dobs.append(d)
          elif f == 2 or f == '2':
               meas.append(f)
               dmeas.append(d)               
               
     plt.scatter(dobs  , obs  , color = 'black', s = 2 )     
     plt.scatter(dmeas , meas , color = 'blue', s = 2 , marker = '^')     
         
     plt.scatter(-1000,-1000, label = 'Measured'  , color = 'black' , s = 1  )    
     plt.scatter(-1000,-1000, label = 'Calculated',  color = 'blue' , s = 1 , marker = '^')    
     
     plt.axis([1990,2019,  0 ,3])
     plt.xlabel(dic_prop['rh']['xlab'], fontsize = fnt_size)
     plt.ylabel('1=Obs,2=Calc', fontsize = fnt_size)
     plt.legend(loc = 'lower left', fontsize = 10)
     
     # dp *******************************************************************************
     ax0 = plt.subplot(gs[3])
     plt.plot(datum    , dp  , label = 'Dew Point Temp.'       , color = 'cyan'   , linestyle = '--' , lw = 1.2 )     
     plt.scatter(datum , cdp , label = 'Dew Point Temp. Comb.' , color = 'black' , s = 4 )     

     plt.axis(dic_prop['dp']['ax'])
     plt.xlabel(dic_prop['dp']['xlab'], fontsize = fnt_size)
     plt.ylabel(dic_prop['dp']['ylab'], fontsize = fnt_size)
     plt.legend(loc = 'lower left', fontsize = 10)
     
     
     ax0 = plt.subplot(gs[4])
     
     obs, meas , dobs, dmeas = [],[],[],[]
     
     for v,f,d in zip(cdp,cfdp,datum):
          if f == 1 or f == '1':
               obs.append(f)
               dobs.append(d)
          elif f == 2 or f == '2':
               meas.append(f)
               dmeas.append(d)               
               
     plt.scatter(dobs  , obs  , color = 'black', s = 2 )     
     plt.scatter(dmeas , meas , color = 'blue', s = 2 , marker = '^')     
           
               
     plt.scatter(-1000,-1000, label = 'Measured'  , color = 'black' , s = 1  )    
     plt.scatter(-1000,-1000, label = 'Calculated',  color = 'blue' , s = 1 , marker = '^')    
     
     plt.axis([1990,2019,  0 ,3])
     plt.xlabel(dic_prop['rh']['xlab'], fontsize = fnt_size)
     plt.ylabel('1=Obs,2=Calc', fontsize = fnt_size)
     plt.legend(loc = 'lower left', fontsize = 10)
     
     
     
     plt.savefig('PLOTS/Combination/Combination_test_' + str(h) + '_' + str(p) + '.pdf' ,  bbox_inches='tight' )
     plt.close()
     
     
#os.mkdir('PLOTS/Combination/')
for p in range(0,16):

     for h in [0,1]:
          print('I am plotting the h , p' , h , p )
          plot_comb(p,h)
         
 
print('done')


