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

database_dir = '/raid8/srvx1/federico/odb_netCDF/$FSCRATCH/ei6'

class netCDF:
     """" Class containing the main properties and 
     functionalities to handle netCDF files
     """
     #global database_dir 
     def __init__(self, res_dir, file_name):
          self.database = database_dir
          self.res_dir = res_dir
          self.file  = file_name # pure name of the file with no path nor file extension , e.g. ERA5_1_16087_
          self.variables = ''
          self.uwind = ''
          self.vwind = ''
          self.datum_u = ''
          self.datum_v = ''
         
     def load(self):
          """ Loading the u,v wind comp. files.
              The 'datum' and 'variables' should be identical """
          
          #print ('DATABASE is', self.database)
          
          uwind_path = self.database + '/' + self.res_dir + '/' + self.file + 'u.nc'
          vwind_path = self.database + '/' + self.res_dir + '/' + self.file + 'v.nc'          
          #print('paths: ',uwind_path,vwind_path)
          
          if os.path.isfile(uwind_path) and os.path.isfile(vwind_path):
               f_u = netCDF4.Dataset(uwind_path) 
               f_v = netCDF4.Dataset(vwind_path)
               
               self.variables = f_u.variables.keys()
               
               self.datum_u = f_u.variables['datum'][0,:]/365.25
               datum_v = f_v.variables['datum'][0,:]/365.25
               self.datum_u = [ 1900 + d for d in datum_v ]
               
               self.uwind = f_u.variables['uwind'][0,12,:]
               self.vwind = f_v.variables['vwind'][0,12,:]
          
          else: 
               raise ValueError('netCDF files:', uwind_path, ' not found!!!')
      
     def printinfo(self):
          print('Basic info for *** ', self.file, '\n')
          print('The variables are: ')
          for v in self.variables:
               print(v)
          print('Datum: ', self.datum_u) 
          print('uwind: ', self.uwind) 
          print('vwind: ', self.vwind) 
          

class Plotter():
     """ Class containing the basic functionalities 
     for plotting the u,v wind components """
     def __init__(self, netCDF):
          self.file  = netCDF
          self.name = netCDF.file
          self.x = netCDF.datum_u
          self.uwind = netCDF.uwind
          self.vwind = netCDF.vwind
          self.out_dir = os.getcwd()+'/PLOTS'
     
     def create_out_fold(self):
          if not os.path.isdir(self.out_dir):
               print('Creating the PLOT output directory: ', self.out_dir)
               os.mkdir(self.out_dir)
               
     def plotter_prop(self):
          fontsize = 12
          plt.xlim(1920,2019)
          plt.ylim(-50,50)
          plt.xlabel('Year', fontsize = fontsize)
          plt.ylabel('Component Speed [m/s]', fontsize = fontsize)
          plt.legend(loc = 'lower left', fontsize = 12)
          
          
     def plotter(self):
          plt.plot(self.x , self.uwind , label = 'uwind', linestyle = '-', color = 'blue' )
          save_path = self.out_dir + '/' + self.name + 'uwind.pdf'
          plt.savefig(save_path,  bbox_inches='tight')
          self.plotter_prop()
          plt.close()
          
          plt.plot(self.x , self.vwind , label = 'vwind', linestyle = '-', color = 'green' )
          self.plotter_prop()
          save_path = self.out_dir + '/' + self.name + 'vwind.pdf'
          plt.savefig(save_path,  bbox_inches='tight')          
          plt.close()
          
#netCDF_file = netCDF(test)     
#netCDF_file.load()
#print ('the file is', netCDF_file.path , netCDF_file.variables)



''' Reading the list of dataset contained in the database directory,
storing the absoulte path '''

dataset = [ res for res in os.listdir(database_dir)]#[:100]

for d in dataset:
    
     res_name = os.listdir(database_dir+'/'+d)[0].replace('u.nc','').replace('v.nc','')
     print('processing: ', res_name)
     netCDF_file = netCDF(d,res_name)
     netCDF_file.load()     
     netCDF_file.printinfo()
 
     Plot = Plotter(netCDF_file)
     Plot.create_out_fold()
     Plot.plotter()
      

      

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