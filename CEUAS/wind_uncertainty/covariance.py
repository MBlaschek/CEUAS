""" Module for extracting the covariance matrix """


import os,sys
import netCDF4
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
#import matplotlib.gridspec as gridspec



# see definition in plotter script 
class netCDF:
     def __init__(self, file=""):
          self.file = file

     def read_data(self, file = '', var = ['fg_dep', 'an_dep']):
         
         data_loaded = {}  # Dictioary containing the extracted information                                                                                                                                                                                                              
         for v in var:
                  data_loaded[v] = {}
                  data_loaded[v]['datum'] = []
                  data_loaded[v]['data' ] = [] # contains either fg_Dep or an_Dep                                                                                                                                                                                                    

                  if os.path.isfile(file):
                      f = netCDF4.Dataset(file)
                      data_loaded[v]['datum'] = [ 1900 + d for d in f.variables['datum'][0,:]/365.25 ] # converting the datum in years                                                                                           
                      data_loaded[v]['data']  = f.variables[v][:,:,:]
                  else: raise ValueError('netCDF files:', path, ' not found!!!')
         return data_loaded 


uwind_file = 'data/ERA5_1_10393_u.nc'

t_file = 'data/ERA5_1_10393_t.nc'


class Covariance:
     def __init__(self,netCDF):
         self.data = netCDF
         return

     def running_mean(self,x, N):
          """ Retunrs the running mean of x, averaged over N """
          cumsum = np.cumsum(np.insert(x, 0, 0)) 
          return (cumsum[N:] - cumsum[:-N]) / float(N)

     def calc_cov(self, array_x, array_y):
          """ Calculate the (cross) covariance matrix for the two arrays x and y """
          cov = np.empty([len(array_x),len(array_x)], dtype = float) # initialize an empty 16*16 matrix (16 pressure levels)
          for x in range(len(array_x)):
               for y in range(len(array_y)):
                    entry = array_x[x] * array_y[y]
                    cov[x,y] = entry
          return cov

     def select_ijentry(self, matrices = '', i = '' , j = ''):
          """ Extract a list of entries of the covariance matrices in the list, 
          e.g. all the "ith,jth" entries of all the matrices.
          Return the values of the entries in a list """
          lista = [ m[i, j] for m in matrices ]
          return lista

     def extract_matrix_list(self, an_dep='', fg_dep='', datums = '', hours=[0,1]):
         """ Extracts a dictionary containing the covariance matrices for a specific hour, one for each date in datums """

         matrices = {'0':'', '12':''}

         for h in hours:
             lista = [] 
             for date in range(len(datums)):
                  an = an_dep[h,:,date]
                  fg = fg_dep[h,:,date]
                  corrMatrix = self.calc_cov(an,fg)
                  lista.append(corrMatrix)
             matrices[str(h)] = lista
         return matrices

     def average_matrix(self,matrices_list='', N='', matrix_index=100):
         """ Given the list of caov. matrices, returns an averaged matrix (for the index matrix_index in the list)
         where each entry is the mean calculated with the N following matrices. """

         averaged = np.empty(16, 16, dtype = float) # initialize an empty 16*16 matrix (16 pressure levels)                                                                                                                            
         new_mlist = matrices_list[matrix_index:N]  # selecting the subset of matrices
         for x in range(16):
              for y in range(16):
                   indices_list = self.select_ijentry(matrices=new_mlist , i=x , j=y)
                   averaged[x,y] = sum(indices_list)/N  
                    
         return averaged
 


class Plotter:

     def plot_prop(self, var='', fg_p = '' , an_p = ''):
         self.font = 15

         self.fg_p = fg_p #  pressure level of the first guess dep.
         self.an_p = an_p #  pressure level of the analysis dep.

         self.pretty_pressure = [10,20,30,50,70,100,150,200,250,300,400,500,700,850,925,1000]
         self.var = var
         self.var_dics = { 'temp': { 'units': 'K' , 'name':'Air Temperature' , 'x_range': [] , 'y_range': []  } }  


     def initialize_dirs(self):
         if not os.path.isdir('plots'):
              os.mkdir('plots')
         if not os.path.isdir('plots/covariances'):
              os.mkdir('plots/covariances')
         if not os.path.isdir('plots/series'):
              os.mkdir('plots/series')
         if not os.path.isdir('plots/histo'):
              os.mkdir('plots/histo')

     def cov_plot(self, matrix, station="", hour = "", date="" , averaged = "" ):
          """ Basic plot for the correlation matrix """

          var = self.var 

          fig,ax = plt.subplots()

          if date and not averaged:
             title = "Stat: " + station + ', H: ' + hour + ', Date: ' + date + ', ' + var
             filename = 'Cov_' + station + '_hour_' + hour.replace(':','') + '_date_' + date + ', ' +var
 
          elif averaged and not date:
               title = "Station: " + station + ', Hour: ' + hour + ', Averaged: ' + averaged + ', ' +var
               filename ='Cov_' + station + '_hour_' + hour.replace(':','') + '_averaged_' + averaged + ', ' +var

          plt.title(title, y=1.03, fontsize = self.font)

          num = len(matrix[0,:])
          Num = range(num)

          vmin, vmax = -10, 10
          color_map= plt.imshow(matrix, interpolation= 'nearest', cmap = 'RdYlBu', vmin = vmin, vmax = vmax ) # nearest serves for discreete grid  # cmaps blue, seismic            
          plt.ylim(-0.5, 15.5)
          plt.xlim(-0.5, 15.5)
          plt.xticks(Num, Num)
          plt.xlabel('Pressure level', fontsize = FONT)
          plt.yticks(Num, Num)
          plt.xlabel('Pressure level', fontsize = FONT)
          bar = plt.colorbar()
          bar.ax.set_ylabel("Covariance", fontsize = FONT)

          #  Creating text values                                                                                                                              
          for i in Num:
               for j in Num:
                    value = '{0:.2f}'.format(matrix[i,j])
                    text = ax.text( j,i, value , ha = 'center' , va = 'center', color = 'black', fontsize = 5)

          plt.savefig('plots/' + name + '.pdf', bbox_inches='tight')
          plt.close()

     def histo(self, X="", colors="", labels="" , bins = ''):
          plt.title('Estimated observation errors for the temperature')

          weights = []
          for x in X:
              w = []
              for e in range(len(x)):
                   w.append(1./len(x)*100)
              weights.append(w)

          print('weights', weights)
          plt.hist(X, bins, histtype='stepfilled',  stacked = False, color = C , label = L , alpha = 0.7 , density = True, weights = weights)

          plt.text(1, 0.95, "pressure(an,fg)="+ str(self.an_p) + ',' + str(self.fg_p) , fontsize= self.font)
          plt.grid(linestyle= ':', color = 'lightgray', lw = 1.2 )
          plt.legend(loc = 'upper right', fontsize = self.font - 3)
          plt.ylabel('Numbers / '+str(bins), fontsize = self.font)
          plt.ylim(0, 10)
          plt.xlim(0, 2,2)
          plt.xlabel(r'Errors [' + self.var_dics[self.var]['units'] + ']', fontsize= self.font)
          plt.savefig('plots/histo/' + self.var + '_anp_' + str(self.an_p) + '_fgp_' + str(self.fg_p) + '.pdf',  bbox_inches='tight')
          plt.close()
          
     

     def time_series(self, means='', datums = '', labels = '', colors = ''):
 
          plt.title('Time series for ' + self.var + ' for (an_dep,fg_dep)=(' + str(self.an_p)+ ',' + str(self.fg_p) + ')' )
          for (m,l,c) in zip(means,labels,colors):
               print(len(m))
               plt.plot( datums[:(len(m))], m, color =c, label =l)
               plt.legend(fontsize = self.font , loc = 'upper right')

          plt.ylim(0,1.2)
          plt.grid(linestyle= ':', color = 'lightgray', lw = 1.2 )
          plt.ylabel('Error [' + self.var_dics[self.var]['units'] + ']', fontsize= self.font) 
          plt.xlabel('Dates', fontsize = self.font)
          plt.savefig('plots/series/timeseries_' + self.var + '_anp_' + str(self.an_p) + '_fgp_' + str(self.fg_p) + '.pdf'  , bbox_inches = 'tight')
          plt.close()








""" Dirs, definitions, select datasets """
base_dir = 'data/'
file_dic = {'temp':'ERA5_1_10393_t.nc' , 'uwind': 'ERA5_1_10393_u.nc'  }
#variables = ['temp','uwind','vwind']
variables = ['temp']
stations = ['Lindenberg']

netCDF = netCDF()
Cov = Covariance(netCDF)
Plot = Plotter()
Plot.initialize_dirs()

matrices_dict = ''
""" Extracting the full (cross)covariance matrices """
for s in stations:

     for v in variables:

          input_file = base_dir + file_dic[v]

          data = netCDF.read_data(file = input_file) #  Loading the necCDF file
        
          andep = data['an_dep']['data'] #  Extracting first guess and analysis departures
          fgdep = data['fg_dep']['data'] 

          print('*** Check the data shape for fg_dep and an_dep: ', andep.shape, fgdep.shape )

          datums = data['fg_dep']['datum']

          # Extract the list of matrices, for all the days of observation """
          print("*** Extracting the matrices """)
          matrices_dict= Cov.extract_matrix_list(an_dep=andep, fg_dep=fgdep, datums=datums, hours=[0])
          
          # slimmed down
          #matrices_dict= Cov.extract_matrix_list(an_dep=andep[:1000], fg_dep=fgdep[:1000], datums=datums[:1000], hours=[0])

          """ Looping over the pressure levels """
          plevels_i = [11]
          plevels_j = [11]

          for i in plevels_i:
               for j in plevels_j:
                    print("*** Processing the i,j entries: ", i , j )
                    means = Cov.select_ijentry(matrices = matrices_dict['0'], i = i , j = j) 
                    means = [ m for m in means if not np.isnan(m) ] 

                    runningmean_30  = Cov.running_mean(means , 30 )
                    #print(runningmean_30)
                    runningmean_60  = Cov.running_mean(means , 60 )
                    runningmean_90  = Cov.running_mean(means , 90 )
                    runningmean_180 = Cov.running_mean(means , 180 )
                    runningmean_365 = Cov.running_mean(means , 365 )

                    X = [runningmean_30,runningmean_60,runningmean_90,runningmean_180,runningmean_365]

                    #C = ['slateblue', 'cyan', 'lime', 'orange', 'gold']                                                                                                                                                                                                            
                    C = ['yellow', 'orange', 'lime', 'cyan', 'slateblue']
                    L = ['Desroziers(1m)', 'Desroziers(2m)', 'Desroziers(3m)', 'Desroziers(6m)', 'Desroziers(1y)']

                    C = ['blue', 'orange', 'green', 'red', 'slateblue'] # color similar to Michi
                    Plot.plot_prop(var = v, fg_p = i, an_p = j)                     
                    Plot.time_series(means=X, datums = datums, labels = L, colors = C)
                    #cov_values = [ x for x in means  if not np.isnan(x) ]  #  Excluding 'nan' values since they give problems with hostograms
                   
                    """ Cleaning the running means from the nan values """
                    rm_1m = [x for x in runningmean_30 if not np.isnan(x) ] 
                    rm_2m = [x for x in runningmean_60 if not np.isnan(x) ] 
                    rm_3m = [x for x in runningmean_90 if not np.isnan(x) ]
                    rm_6m = [x for x in runningmean_180 if not np.isnan(x) ]
                    rm_1y = [x for x in runningmean_365 if not np.isnan(x) ]

                    """ Plotting the histograms """

                    X = [ rm_1m, rm_2m, rm_3m, rm_6m, rm_1y ]
                    bins = 150
                    Plot.histo(X= X, colors = C, labels = L, bins = bins)





'''
print(X)
plt.title('Estimated observation errors for the temperature')
Bins = 50
FONTSIZE = 15
n, bins, patches = plt.hist(X, Bins, histtype='stepfilled' ,  stacked = False, color = C , label = L , alpha = 0.7 , normed = True)
plt.text(1, 0.28,"pressure(an,fg)=(0,0)", fontsize= FONTSIZE)
plt.grid(linestyle= ':', color = 'lightgray', lw = 1.2 )
plt.legend(loc = 'upper right', fontsize = FONTSIZE - 3)
plt.ylabel('Numbers / '+str(Bins), fontsize = FONTSIZE)
plt.ylim(0, 0.30)
plt.xlabel(r'Errors [m/s]', fontsize=15)
plt.savefig('plots/Temperature.pdf',  bbox_inches='tight')
plt.close()
'''





'''
""" Plotting some covariances """
for hour in [0,1]:
    hour_name = str(hour).replace('0','00:00').replace('1','12:00')
    for date in [0,100,1000,5000,7000]:

        date_name = str(datums[date])

        print("*** I am extracting the covariances for the date ", date_name , " and hour ", hour_name )

        an = andep[hour,:,date]
        fg = fgdep[hour,:,date]

        corr_matrix = calc_cov(an,fg)  # Extracting the matrix

        cov_plot(corr_matrix , station = station, hour = hour_name , date = date_name)  # Plotting
'''

'''
# Example
date = str (data['fg_dep']['datum'][0]) 
hour = 0
Hour = '00'
an = andep[hour,:,0]
fg = fgdep[hour,:,0]
corr_matrix = calc_cov(an,fg)
print('Date is: ', date )
cov_plot(corr_matrix , station = 'Lindenberg', hour = Hour , date = date)
print('andep, fgdep', an, fg)
'''
        


