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

class netCDF:
     """" Class containing the main properties and functionalities to handle netCDF files
          Args:
               res_dir : root directory containing the results 
               file_name : base name
               var : variable type. i.e. 'temp','uwind','vwind' , 'hum'       
     """
     #global database_dir 
     def __init__(self, database_dir, res, file_name, hour , plevel):
          """ Information identifying the file 
              Args:
                   database_dir: root directory where the results for each dataset are stored e.g /raid8/srvx1/federico/odb_netCDF/$New_Results
                   res: specific datased considered e.g. 1759 
                   file_name: base file name stripped from extension and specific variable e.g. ERA5_1_16087_          
          """
          self.database = database_dir
          self.res = res
          self.file_name  = file_name
          self.uwind = ''
          self.vwind = ''
          self.temp  = ''
          self.sh    = ''
          self.hour = hour
          self.plevel = plevel
          
          
          self.datum_uwind = ''
          self.datum_vwind = ''
          self.datum_temp = ''         
          self.datum_sh = ''
          
          self.variables = ''
          
     def load(self, var=['temp','uwind','vwind','rh','dp' ,'sh']):
          """ Loading the u,v wind comp. files.
              The 'datum' and 'variables' should be identical """
          #print ('DATABASE is', self.database)
          file_dir = self.database + '/' + self.res + '/' + self.file_name
          print ('base dir:' , file_dir )
          print ('res dir:'  , self.res )
          print ('file name:', self.file_name )
          
          # example: /raid8/srvx1/federico/odb_netCDF/$New_Results/1759/1:87418
          paths = { 'dir': file_dir , 
                    'uwind': 'u.nc' , 'vwind':'v.nc' , 'temp':'t.nc', 'sh':'sh.nc' , 'dp':'dp.nc' , 'rh':'rh.nc' ,
                    
                    } 
          
          #uwind_path = file_dir + 'u.nc'

          data_loaded = {}
                
          for v in var:
               data_loaded[v] = {}
               data_loaded[v]['datum'] = []
               data_loaded[v]['data']  = [] # initialising empty dictionary
               
               path = file_dir +  paths[v] # extracting the path from the dictionary
               if os.path.isfile(path):
                    f = netCDF4.Dataset(path) 
                    data_loaded[v]['datum'] = [ 1900 + d for d in f.variables['datum'][0,:]/365.25 ] # converting the datum in years  
                    data_loaded[v]['data']  = f.variables[v.replace('temp','temperatures')][self.hour,self.plevel,:]
               else: raise ValueError('netCDF files:', path, ' not found!!!')
   
          self.datum_uwind = data_loaded['uwind']['datum']
          self.datum_vwind = data_loaded['vwind']['datum']
          self.datum_temp  = data_loaded['temp' ]['datum']
          self.datum_sh    = data_loaded['sh'   ]['datum']
          self.datum_rh    = data_loaded['rh'   ]['datum']
          self.datum_dp    = data_loaded['dp'   ]['datum']
          
          
          
          self.uwind = data_loaded['uwind']['data']
          self.vwind = data_loaded['vwind']['data']
          self.temp  = data_loaded['temp' ]['data']
          self.sh    = data_loaded['sh'   ]['data'] 
          self.rh    = data_loaded['rh'   ]['data']                    
          self.dp    = data_loaded['dp'   ]['data']          
          
               

     def printInfo(self):
          print('Basic info for *** ', self.file_name, '\n')
          print('The available variables are: ')
          for v in self.variables:
               print(v)

     def analyser(self):
          """Module that analyses the netCDF file"""
          print('Will do something')
          
          

class Plotter():
     """ Class containing the basic functionalities 
     for plotting the u,v wind components """
     def __init__(self, netCDF, var=''):
          self.file  = netCDF
          self.name = netCDF.file_name
          self.var     = var

          self.x_data = ''
          self.y_data = ''
          
          '''
          """ wind components """          
          self.datum_u = netCDF.datum_u
          self.datum_v = netCDF.datum_u
          self.uwind = netCDF.uwind
          self.vwind = netCDF.vwind
          """ temperature  """          
          self.datum_temp = netCDF.datum_u
          self.temp       = netCDF.uwind
          '''
          
     
     def load_xy(self):
          ''' Loading the x and y data to plot and analyze, according the chosen variable'''
          var = self.var
          print('The var is: ', var )
          if var == 'temp':
               self.x_data = self.file.datum_temp
               self.y_data = self.file.temp
          elif var == 'uwind':
               self.x_data = self.file.datum_uwind
               self.y_data = self.file.uwind
          elif var == 'vwind':
               self.x_data = self.file.datum_vwind
               self.y_data = self.file.vwind        
          elif var == 'sh':
               self.x_data = self.file.datum_sh
               self.y_data = self.file.sh       
               self.y_data = self.file.vwind        
          elif var == 'rh':
               self.x_data = self.file.datum_rh
               self.y_data = self.file.rh                      
          elif var == 'dp':
               self.x_data = self.file.datum_dp
               self.y_data = self.file.dp                     
               
               
          else:
               raise ValueError('Unknown variable:', var, ' !!!')
               
     
     def create_outDir(self,out_path):
          if not os.path.isdir(out_path):
               print('Creating the PLOT output directory: ', out_path)
               os.mkdir(out_path)
               
     
     def style_dic(self):          
          dic_prop = { 'uwind':{'xlab':'Year', 'ylab':'Speed [m/s]'             ,'leg':'uwind'      ,         'ax':[1920,2019, -60  ,60 ] ,'c':'blue'         } ,
                       'vwind':{'xlab':'Year', 'ylab':'Speed [m/s]'             ,'leg':'vwind'      ,         'ax':[1920,2019, -60  ,60 ] ,'c':'cyan'         } ,
                        'temp':{'xlab':'Year', 'ylab':'Temperature [K]'         ,'leg':'Temperature',         'ax':[1920,2019,  200 ,300] ,'c':'orange'       } ,
                          'sh':{'xlab':'Year', 'ylab':'Specific Humidity [???]' ,'leg':'Specific Humidity',   'ax':[1920,2019,  0   , 50] ,'c':'lime'         } ,
                          'rh':{'xlab':'Year', 'ylab':'Relative Humidity [???]' ,'leg':'Relative Humidity',   'ax':[1920,2019,  0   , 1 ] ,'c':'lightgray'    } ,                       
                           'dp':{'xlab':'Year', 'ylab':'Dew Point [K]'           ,'leg':'Upper Air Dew Point','ax':[1920,2019,  200 ,300] ,'c':'mediumpurple'} ,
                       
                           }
          return dic_prop
               
     def plotter_prop(self, xlabel = True):
          var = self.var
          dic_prop = self.style_dic()
          fnt_size = 12
          plt.title('Station: ' + res_name.split("_")[2] + ' - Hour: ' + str(self.file.hour) + ' Press. Lev: ' + str(self.file.plevel) )
          #ax_lim = min(self.)
          minimum = min(self.x_data)-5 
          axes_lim = dic_prop[var]['ax']
          limits = [minimum, axes_lim[1] , axes_lim[2] , axes_lim[3]  ]
          if var =='temp': 
               plt.text(minimum+2, 290, 'Dataset: ' + res_dir, fontsize = 12 , color = 'gray')
               #plt.text(minimum+2, 280, 'Station: ' + res_name.split("_")[2] , fontsize = 12 , color = 'red')            
          print(limits)
          plt.axis(limits)
          if xlabel: plt.xlabel(dic_prop[var]['xlab'], fontsize = fnt_size)
          plt.ylabel(dic_prop[var]['ylab'], fontsize = fnt_size)
          plt.legend(loc = 'lower left', fontsize = 10)
          
     def plotter(self, out_dir='', save = True, xlabel = True):
          self.load_xy()                     
          x = self.x_data 
          y = self.y_data
          
          print('x: ', x , 'y: ', y)
          dic_prop = self.style_dic()
          plt.plot( x , y , label = dic_prop[self.var]['leg'], linestyle = '-', color = dic_prop[self.var]['c'] )
          self.plotter_prop(xlabel = xlabel)
          plt.grid(linestyle = ':')                    
          if save: 
               save_path = out_dir + '/' + self.name + '_' + self.var +'.pdf'
               self.create_outDir(out_path=out_dir)
               plt.savefig(save_path,  bbox_inches='tight')
               plt.close()
 
 
                     
#netCDF_file = netCDF(test)     
#netCDF_file.load()
#print ('the file is', netCDF_file.path , netCDF_file.variables)


''' Reading the list of dataset contained in the database directory,
storing the absolute path '''



""" Producing single plots for each variable """
'''
for d in dataset:
     res_dir = database_dir+'/'+d # e.g. /raid8/srvx1/federico/odb_netCDF/$New_Results/1759/
     res_name = os.listdir(res_dir)[0].replace('u.nc','').replace('v.nc','').replace('t.nc','')

     for v in vars:
          print('processing: ', res_name) # e.g. ERA5_1759_2:11903_
          netCDF_file = netCDF(database_dir,d,res_name) 
          
          
          netCDF_file.load(v)     
          netCDF_file.printInfo()
      
          Plot = Plotter(netCDF_file, var=v )
          Plot.plotter(out_dir= out_dir)
'''  

from pylab import rcParams
rcParams['figure.figsize']= 10, 13# first is the height, second is the width


#database_dir = '/raid8/srvx1/federico/odb_netCDF/netCDF_1_1759_1761/1'


#dataset = [ res for res in os.listdir(database_dir)]
out_dir = os.getcwd() + '/CIAO/' 

os.system('mkdir Plots_Station_10393')

database_dir = '/raid8/srvx1/federico/odb_netCDF/Ready_ForHumidity/1/'
#database_dir = '/raid8/srvx1/federico/odb_netCDF/redo_10393/1/'
dataset = ['10393']
for d in dataset:
     for h in [0,1]:
          for p in range(0,16):
               gs = gridspec.GridSpec(6,1)
     
               res_dir = database_dir+'/'+d # e.g. /raid8/srvx1/federico/odb_netCDF/$New_Results/1759/
               res_name = os.listdir(res_dir)[0].replace('u.nc','').replace('v.nc','').replace('t.nc','')     
            
               ax0 = plt.subplot(gs[0])
            
               netCDF_file = netCDF(database_dir,d,res_name, hour = h , plevel = p)
               netCDF_file.load()     
               
               Plot = Plotter(netCDF_file, var='temp' )
               Plot.plotter(out_dir= out_dir, save = False , xlabel = False)
               plt.tight_layout()
               ax0.xaxis.set_major_formatter(plt.NullFormatter())
            
            
               ax1 = plt.subplot(gs[1])
               netCDF_file.load()     
               Plot = Plotter(netCDF_file, var='uwind' )
               Plot.plotter(out_dir= out_dir, save = False , xlabel = False)  
               plt.tight_layout()
               ax1.xaxis.set_major_formatter(plt.NullFormatter())
                    
               ax2 = plt.subplot(gs[2])
               netCDF_file.load()     
               Plot = Plotter(netCDF_file, var='vwind' )
               Plot.plotter(out_dir= out_dir, save = False , xlabel = False) 
               plt.tight_layout()
               ax2.xaxis.set_major_formatter(plt.NullFormatter())  
            
               ax3 = plt.subplot(gs[3])
               netCDF_file.load()     
               Plot = Plotter(netCDF_file, var='sh' )
               Plot.plotter(out_dir= out_dir, save = False , xlabel = False) 
               plt.tight_layout()
               ax3.xaxis.set_major_formatter(plt.NullFormatter())
            
               ax3 = plt.subplot(gs[4])
               netCDF_file.load()     
               Plot = Plotter(netCDF_file, var='rh' )
               Plot.plotter(out_dir= out_dir, save = False , xlabel = False) 
               plt.tight_layout()
               ax3.xaxis.set_major_formatter(plt.NullFormatter())  
               
               ax5 = plt.subplot(gs[5])
               netCDF_file.load()     
               Plot = Plotter(netCDF_file, var='dp' )
               Plot.plotter(out_dir= out_dir, save = False , xlabel = True) 
               plt.tight_layout()          
               
               plt.savefig('Plots_Station_10393/' + netCDF_file.file_name + '_' + str(h) + '_' + str(p) +'_' +'_variables.pdf' , bbox_inches='tight' )    
                

print('done')


'''
test = database_dir +'/16087/ERA5_1_16087_u.nc'
test_false = 'dsgh'

test_netcdf(file=test, print_info=True)

for var,val in vars.items():

    with netCDF4.Dataset( sample1 + var + '.nc') as f:
        plt.plot(f.variables['datum'][:]/365.25,f.variables[val][0,12,:] , 
                 label='ERA_1', color = 'blue')
    
    with netCDF4.Dataset(sample2 + var +  '.nc') as f:
        plt.plot(f.variables['datum'][0,:]/365.25,f.variables[val][0,12,:],
                label ='ERA5_erapresat', color = 'green')
    
    with netCDF4.Dataset(sample3 + var + '.nc') as f:
        plt.plot(f.variables['datum'][0,:]/365.25,f.variables[val][0,12,:],
                 label = 'feedbackmerged' , color = 'orange' )

    print(f.variables.keys())
    plt.ylabel(var)
    plt.xlim(20,120)
    plt.xlabel('Year')
    plt.legend(loc = 'lower left')
    plt.savefig(var + '.png',  bbox_inches='tight')
    plt.close()


print('*** Finished')
'''
