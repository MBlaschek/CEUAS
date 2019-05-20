""" Module for extracting the covariance matrix 

    Author: Ambrogi Federico, federico.ambrogi@univie.ac.at
   
"""


import os,sys
import netCDF4
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import argparse
import pandas as pd
#import matplotlib.gridspec as gridspec

""" Running mode: forcing the creation of the numpy files """
parser = argparse.ArgumentParser(description="Force the creation/update of the numpy files storing the covariances")
parser.add_argument('--force' , '-f', 
                    help="Optional: select True to force the creation of the numpy files storing the covariances" ,
                    default = 'False',
                    type = bool)
args = parser.parse_args()
force = args.force


# see definition in plotter script 
class netCDF:
     """ Class containig the basic information of the netCDF files """
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






class Covariance:
     """ Class with the main methods for computing and analising covariances """
     def __init__(self,netCDF):
         self.data = netCDF
         return

     def remove_outliers(self, data, min= '', max= '', cut= ''):
          """ Replacing the outlier values with nan (to keep vector of same length)"""
          q25, q75 = np.nanpercentile(data,25), np.nanpercentile(data,75)
          iqr = q75 - q25
          cut_off = iqr * cut
          lower, upper = q25-cut_off, q75+cut_off
          #print(q25, q75 , 'percentiles') 
          #input('check')
          cleaned = []
 
          for d in data:
              if d > lower and d < upper:
                   cleaned.append(d)
              else:
                   cleaned.append(np.nan)         
          #print('cleaned is', cleaned)
          return cleaned

     def running_mean_old(self, data = '', n='', datums=''):
          """ Returns the running mean of x, averaged over N  
              Data is reduced by removing outliers """
          #print(len(data), len(datums))
          #input('ciao')
          not_nans         = [ x for x in data if not np.isnan(x) ]
          not_nans_indexes = [ data.index(x) for x in data if not np.isnan(x) ]
     
          datums_nans      = [ datums[i] for i in not_nans_indexes ] # extracting not nans values and corresponding index in datums (for plotting) 

          #print (not_nans[:100], datums_nans[:100])
          #input('check values')
          cumsum = np.cumsum(np.insert(not_nans, 0, 0))
          means = (cumsum[n:] - cumsum[:-n]) / float(n)
          means = np.sqrt(means)
          return means , datums_nans
    
     def running_mean(self,data='',n=''): 
          df = pd.DataFrame( {'mean': data } ) 
          mean = df.rolling(n, center = True, min_periods=1).mean()
          #print('mean in rolling is', mean['mean'] )
          conv = np.array(mean['mean'])
          return conv 

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

         matrices = {}

         for h in hours:
             lista = [] 
             for date in range(len(datums)):
                  an = an_dep[h,:,date]
                  fg = fg_dep[h,:,date]
                  corrMatrix = self.calc_cov(an,fg)
                  lista.append(corrMatrix)
             matrices[str(h)] = lista
         return matrices

     def average_matrix(self,matrices_list='', N='', matrix_index=""):
         """ Given the list of caov. matrices, returns an averaged matrix (for the index matrix_index in the list)
         where each entry is the mean calculated with the N following matrices. 
         The average is computed considering only the actual 'n' entries with real values """

         averaged = np.empty([16, 16], dtype = float) # initialize an empty 16*16 matrix (16 pressure levels)                                                                                                                            
         new_mlist = matrices_list[matrix_index:matrix_index+N]  # selecting the subset of matrices from index=matrix_index to matrix_index+N
         for x in range(16):
              for y in range(16):
                   valueslist = self.select_ijentry(matrices=new_mlist , i=x , j=y)
                   valueslist_nooutliers = self.remove_outliers(valueslist, min= 0.25, max= 0.75, cut= 1.5) 
                   not_nan = [ x for x in valueslist_nooutliers if not np.isnan(x) ]
                   #print(not_nan)
                   N = len(not_nan)
                   try:
                        averaged[x,y] = sum(not_nan)/N  
                   except ZeroDivisionError:
                        averaged[x,y] = np.nan
         return averaged
 


class Plotter:

     def plot_prop(self, var='', fg_p = '' , an_p = '' , hour = ''):
         self.font = 15

         self.fg_p = fg_p #  pressure level of the first guess dep.
         self.an_p = an_p #  pressure level of the analysis dep.

         self.pretty_pressure = [10,20,30,50,70,100,150,200,250,300,400,500,700,850,925,1000]
         self.var = var
         self.hour = hour
         self.var_dics = { 'temp' : { 'units': 'K'     , 'name':'Air Temperature'  , 'x_range': [] , 'y_range': []  } ,
                           'uwind': { 'units': '[m/s]' , 'name':'u wind Component' , 'x_range': [] , 'y_range': [] } ,  
                           'vwind': { 'units': '[m/s]' , 'name':'v wind Component' , 'x_range': [] , 'y_range': [] } ,
                           'dp'   : { 'units': '[K]'   , 'name':'Dew Point Temp.'  , 'x_range': [] , 'y_range': [] } ,
                           'rh'   : { 'units': ''      , 'name':'Relative Hum.   ' , 'x_range': [] , 'y_range': [] } ,
                           }



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

          if not averaged:
             title = "Stat: " + station + ', H: ' + hour + ', Date: ' + date + ', ' + var
             filename = 'Cov_' + station + '_hour_' + hour.replace(':','') + '_date_' + date + '_' +var
 
          elif averaged :
               title = var.replace('temp','Temp.') + " , Stat: " + station + ', H: ' + str(hour) + ', Date: ' + str(date)
               filename ='Cov_' + station + '_hour_' + str(hour).replace(':','') + '_averaged_' + str(date) + '_' + var 

          plt.title(title, y=1.03, fontsize = self.font-2)

          num = len(matrix[0,:])
          Num = range(num)

          vmin, vmax = -5, 5
          color_map= plt.imshow(matrix, interpolation= 'nearest', cmap = 'RdYlBu', vmin = vmin, vmax = vmax ) # nearest serves for discreete grid  # cmaps blue, seismic            
          plt.ylim(-0.5, 15.5)
          plt.xlim(-0.5, 15.5)
          plt.xticks(Num, Num)
          plt.xlabel('Pressure level', fontsize = self.font)
          plt.yticks(Num, Num)
          plt.xlabel('Pressure level', fontsize = self.font)
          bar = plt.colorbar()
          bar.ax.set_ylabel("Covariance", fontsize = self.font)

          #  Creating text values                                                                                                                              
          for i in Num:
               for j in Num:
                    value = '{0:.2f}'.format(matrix[i,j])
                    text = ax.text( j,i, value , ha = 'center' , va = 'center', color = 'black', fontsize = 5)

          if not os.path.isdir('plots/covariances/'+station): os.mkdir('plots/covariances/'+station)
          plt.savefig('plots/covariances/' + station + '/' + filename + '.pdf', bbox_inches='tight')
          plt.close()

     def histo(self, X="", colors="", labels="" , bins = ''):
          plt.title('Estimated errors for ' + self.var + ' for (an_dep,fg_dep)=(' + str(self.an_p)+ ',' + str(self.fg_p) + ')' , y=1.03)

          weights = []
          for x in X:
              w = []
              for e in range(len(x)):
                   w.append(1./len(x)*100)
              weights.append(w)

          plt.hist(X, bins, histtype='stepfilled',  stacked = False, color = C , label = L , alpha = 0.7 , density = True)

          plt.grid(linestyle= ':', color = 'lightgray', lw = 1.2 )
          plt.legend(loc = 'upper right', fontsize = self.font - 3)
          plt.ylabel('Numbers / '+str(bins), fontsize = self.font)
          plt.ylim(0, 15)
          plt.xlim(0, 2.1 )
          plt.xlabel(r'Errors [' + self.var_dics[self.var]['units'] + ']', fontsize= self.font)
          plt.savefig('plots/histo/histo_' + self.var + '_hour_' + self.hour + '_anp_' + str(self.an_p) + '_fgp_' + str(self.fg_p) + '.pdf',  bbox_inches='tight')
          plt.close()
          
     def time_series(self, means='', datums = '', labels = '', colors = '', interval=50):
 
          plt.title('Time series for ' + self.var + ' for (an_dep,fg_dep)=(' + str(self.an_p)+ ',' + str(self.fg_p) + ')' , y=1.03)
          for (m,d,l,c) in zip(means,datums,labels,colors):
               #print('plots: datums , m ', len(d), len(m) )
               #input('prova')
               plt.plot( d[:(len(m))][::interval], m[::interval], color =c, label =l)
          plt.legend(fontsize = self.font-3 , loc = 'upper left')

          plt.ylim(0, 1.2)
          plt.grid(linestyle= ':', color = 'lightgray', lw = 1.2 )
          plt.ylabel('Error [' + self.var_dics[self.var]['units'] + ']', fontsize= self.font) 
          plt.xlabel('Dates', fontsize = self.font)
          plt.savefig('plots/series/timeseries_hour_' + self.hour + '_' + self.var + '_anp_' + str(self.an_p) + '_fgp_' + str(self.fg_p) + '.pdf'  , bbox_inches = 'tight')
          plt.close()


""" Dirs, definitions, select datasets """
base_dir = 'data/'
file_dic = {'temp':'ERA5_1_10393_t.nc' , 'uwind': 'ERA5_1_10393_u.nc'  , 'vwind' : 'ERA5_1_10393_v.nc' , 'dp':'ERA5_1_10393_dp.nc' , 'rh':'ERA5_1_10393_rh.nc'}
#variables = ['temp','uwind','vwind']
variables = ['temp', 'uwind','vwind','rh','dp']
stations = ['Lindenberg']

netCDF = netCDF()
Cov = Covariance(netCDF)
Plot = Plotter()
Plot.initialize_dirs()


# *********************************************
# Extracting the cross-covarinace matrices
# ********************************************* 

""" Initialize and empty dictionary of all the matrices """
matrices = {}

""" Extracting the full (cross)covariance matrices, and store in dictionary """

if not os.path.isfile('covariance_matrices.npy') and force != True :
     print(' *** Extracting the covariance matrices and storing in a numpy dictionary \n')
     for s in stations:
          print('  *** Processing the Station: ', s , '\n')
          matrices[s]= {}
          for v in variables:
               print('   *** Analysing the variable: ', v , '\n')
               matrices[s][v]={}
               input_file = base_dir + file_dic[v]
               data = netCDF.read_data(file = input_file) #  Loading the necCDF file
               andep = data['an_dep']['data'] #  Extracting first guess and analysis departures
               fgdep = data['fg_dep']['data'] 

               print('*** Check the data shape for fg_dep and an_dep: ', andep.shape, fgdep.shape )

               datums = data['fg_dep']['datum']

          # Extract the list of matrices, for all the days of observation """
               print("*** Extracting the values::: """)
          #matrices_dict[station][v] = Cov.extract_matrix_list(an_dep=andep, fg_dep=fgdep, datums=datums, hours=[0])
          # slimmed down
               #all_matrices = Cov.extract_matrix_list(an_dep=andep[:100], fg_dep=fgdep[:100], datums=datums[:100], hours=[0,1])
               all_matrices = Cov.extract_matrix_list(an_dep=andep, fg_dep=fgdep, datums=datums, hours=[0,1])
               matrices[s][v]['0']   = all_matrices['0']
               matrices[s][v]['1']   = all_matrices['1']
               matrices[s]['datums'] = datums
     np.save('covariance_matrices', matrices)

else:
     print('*** Loading the covariance_matrices.npy dictionary') 
     matrices = np.load('covariance_matrices.npy', allow_pickle = True).item()


print('The matrices are', matrices.keys())

print('The matrices are', matrices['Lindenberg'].keys())

'''
""" Plotting the covariance matrices """
for s in stations:
     datums = matrices[s]['datums'] 

     for v in variables:
          Plot.plot_prop(var=v)                                                                                                                                                                                                                                        
          for h in [0,1]: 
              lista = matrices[s][v][str(h)]
              #print(lista)
              for m_index in [100,500,1000,1500,2500,5000,7500,8500]:
                   print('Producing hte plots for', m_index)
                   averaged = Cov.average_matrix(matrices_list= lista, N=365, matrix_index= m_index) # considering the m_indexth matrix, and the 365 matrices that follows
                   plotting = Plot.cov_plot(averaged, station=s , hour = h, date= datums[m_index] , averaged = True )

'''

variables = ['temp','uwind','vwind']
variables = ['temp']
#variables = ['uwind','vwind']

""" Plotting the time series of the running mean """
for s in stations:
     datums = matrices[s]['datums']
     for v in variables:
          matrices_dic = matrices[s][v]
          for h in [0,1]:
               """ Looping over the pressure levels """
               plevels_i = range(16)
               plevels_j = range(16)

               plevels_i = [3,8,11]
               plevels_j = [3,8,11]  # 11 is 500 hPa

               plevels_i = [11]
               plevels_j = [11]  # 11 is 500 hPa
               #for i in plevels_i:
               #     for j in plevels_j:
               for i,j in zip(plevels_i,plevels_j): # no off-diagonal elements

                         print("*** Processing the i,j entries: ", i , j )
                         values = Cov.select_ijentry(matrices = matrices_dic[str(h)], i = i , j = j) 

                         values_cleaned = Cov.remove_outliers(values, min= 0.25, max= 0.75, cut= 1.5 )
            
                         means, dates = [],[]
                         for n in [30,60,90,180,365]:
                              runningmean, date  = Cov.running_mean_old(data= values_cleaned, n= n, datums= datums)
                              means.append(runningmean)
                              dates.append(date)

                         C = ['yellow', 'orange', 'lime', 'cyan', 'slateblue']
                         L = ['Desroziers (1m)', 'Desroziers (2m)', 'Desroziers (3m)', 'Desroziers (6m)', 'Desroziers (1y)']
                         C = ['steelblue', 'orange', 'forestgreen', 'red', 'slateblue'] # color similar to Michi

                         
                         """ Plotting the time series """
                         Plot.plot_prop(var = v, fg_p = i, an_p = j , hour = str(h))                     
                         Plot.time_series(means= means, datums= dates, labels= L, colors= C ,interval= 25)

                         """ Plotting the histograms """
                         bins = 30
                         Plot.histo(X= means, colors = C, labels = L, bins = bins)









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
        


