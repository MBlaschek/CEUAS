""" Module for extracting the covariance matrix 

    Author: Ambrogi Federico, federico.ambrogi@univie.ac.at

When calculating the error,
to check: first calculate the average then take the sqrt of the covariances or 
first calculate the sqrt then calculate the average ?
   
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
          """ Replacing the outlier values with np.nan (to keep vector of same length) 
              Returns: cleaned             = list of values without outliers (filled with nans replacing outliers)
                       outliers            = list of outlier values 
                       lower,upper, median = outliers delimiter and median values """
          q25, q75 = np.nanpercentile(data,25), np.nanpercentile(data,75)
          cut_off = (q75 - q25) * cut
          lower, upper = q25-cut_off, q75+cut_off
    
          median = np.nanmedian(data)
          cleaned, outliers = [],[]
          for d in data:
              if d >= lower and d <= upper or np.isnan(d):
                   cleaned.append(d)
                   outliers.append(np.nan)                   
              else:
                   cleaned.append(np.nan)
                   outliers.append(d) 
          return cleaned, outliers, lower, upper, median

     def running_mean_old(self, data = '', n='', datums=''):
          """ Returns the running mean of x, averaged over N  
              Data is reduced by removing outliers """
          not_nans         = [ x for x in data if not np.isnan(x) and x > -1000 ] #-1000 are outlier values
          not_nans_indexes = [ data.index(x) for x in data if not np.isnan(x) ]
          datums_nans      = [ datums[i] for i in not_nans_indexes ] # extracting not nans values and corresponding index in datums (for plotting) 
          cumsum = np.cumsum(np.insert(not_nans, 0, 0))
          means = (cumsum[n:] - cumsum[:-n]) / float(n)
          means = np.sqrt(np.absolute(means)) #  !!! check if correct !!!
          return means , datums_nans
    
     """
     def running_mean(self,data='',n=''): 
          df = pd.DataFrame( {'mean': data } ) 
          mean = df.rolling(n, center = True, min_periods=1).mean()
          #print('mean in rolling is', mean['mean'] )
          conv = np.array(mean['mean'])
          return conv 
     """

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
         """ Extracts a dictionary containing the covariance matrices, one for each date in datums """

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
         """ Given the list of cov. matrices and the index of the chosen matrix, returns an averaged matrix.
         The average is calculated using the N following matrices. 
         Only the the entries with actual values are used, nans values are skipped """

         averaged = np.empty([16, 16], dtype = float) # initialize an empty 16*16 matrix (for the 16 standard pressure levels)                                                                                                                            
         new_mlist = matrices_list[matrix_index:matrix_index+N]  # selecting the subset of matrices from index=matrix_index to matrix_index+N
         for x in range(16):
              for y in range(16):
                   valueslist = self.select_ijentry(matrices=new_mlist , i=x , j=y)
                   cleaned, outliers, lower, upper, median = self.remove_outliers(valueslist, min= 0.25, max= 0.75, cut= 1.5) 
                   #print('valueslist', valueslist, cleaned)
                   not_nan = [ x for x in cleaned if not np.isnan(x) ]
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
         self.var_dics = { 'temp' : { 'units': '[K]'     , 'name':'Air Temperature'  , 'x_range': [] , 'y_range': [] } ,
                           'uwind': { 'units': '[m/s]' , 'name':'u wind Component' , 'x_range': [] , 'y_range': [] } ,  
                           'vwind': { 'units': '[m/s]' , 'name':'v wind Component' , 'x_range': [] , 'y_range': [] } ,
                           'dp'   : { 'units': '[K]'   , 'name':'Dew Point Temp.'  , 'x_range': [] , 'y_range': [] } ,
                           'rh'   : { 'units': ''      , 'name':'Relative Hum.'    , 'x_range': [] , 'y_range': [] } ,
                           'speed': { 'units': ''      , 'name':'Wind Speed'       , 'x_range': [] , 'y_range': [] } ,
                           'direction': { 'units': ''  , 'name':'Wind Direction'       , 'x_range': [] , 'y_range': [] } ,

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
         if not os.path.isdir('plots/outliers'):
              os.mkdir('plots/outliers')
         if not os.path.isdir('plots/outliers/zoom'):
              os.mkdir('plots/outliers/zoom')

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
          plt.xlabel('Pressure level an_dep [hPa]', fontsize = self.font)
          plt.yticks(Num, Num)
          plt.xlabel('Pressure level fg_dep [hPa]', fontsize = self.font)
          ax.set_xticklabels(labels = self.pretty_pressure, fontsize = self.font-1)
          ax.set_yticklabels(labels = self.pretty_pressure, fontsize = self.font-1)

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

     def outliers (self, outliers='', cleaned = '', datums = '', bins = 25, median = '', lower = '', upper = ''):
          """ Plot the cumulative distributions of the number of outliers and true data """
          plt.title('Outlier for ' + self.var + ' for (an_dep,fg_dep)=(' + str(self.an_p)+ ',' + str(self.fg_p) + ')' , y=1.03)
 
          #median = '{:2.f}'.format(median,4)
          #lower = '{:2.f}'.format(lower,4)
          #upper = '{:2.f}'.format(upper,4)

          plt.text(0.1 , 1790 , 'Median [{:.2f}]'.format(median) , fontsize = self.font -3 )
          plt.text(0.1 , 1700 , 'Lower quart. [{:.2f}]'.format(lower) , fontsize = self.font -3 )
          plt.text(0.1 , 1610 , 'Upper quart. [{:.2f}]'.format(upper) , fontsize = self.font -3)

          cleaned = [ c for c in cleaned if not np.isnan(c)]
          colors = ['slateblue','limegreen']
          labels = ['Outliers ['+ str(len(outliers)) + ']', 'Data [' + str(len(cleaned)) + ']' ]
          plt.hist([outliers,cleaned], bins, histtype='stepfilled',  stacked = False, color = colors , label = labels , alpha = 0.7 , density = False)

          #plt.axvline(x= median, color = 'blue', ls = '-' , label = 'Median')
          #plt.axvline(x= lower , color = 'cyan', ls = '--' , label = 'Lower quart.')
          #plt.axvline(x= upper , color = 'cyan', ls = '--' , label = 'Upper quart.')
 
          plt.legend(fontsize = self.font-3 , loc = 'upper right')

          plt.xlim(0, 3)
          plt.ylim(0, 1900)
          plt.grid(linestyle= ':', color = 'lightgray', lw = 1.2 )
          plt.ylabel('Number ', fontsize= self.font)
          plt.xlabel('Covariance', fontsize = self.font)
          plt.savefig('plots/histo/outliers_hour_' + self.hour + '_' + self.var + '_anp_' + str(self.an_p) + '_fgp_' + str(self.fg_p) + '.pdf'  , bbox_inches = 'tight')
          plt.close()

     def outliers_series(self, corr='', out='', datums='', N= '', interval= ''):
          """ Plot the numbe rof outliers and true data vs date """

          fig = plt.figure(figsize = (11,8))

          ax1 = fig.add_subplot(211)
          plt.title('Outliers series over '+ str(N) +  ' days for ' + self.var + ', (an_dep,fg_dep)=(' + str(self.an_p)+ ',' + str(self.fg_p) + ')' , y=1.03)
          ax1.plot(datums, corr, label= 'True data' , color = 'cyan' )
          ax1.plot(datums, out , label= 'Outliers' , color = 'lime' )
          ax1.set_ylim(0, int(max(corr)+max(corr)/10) )
          ax1.grid(linestyle= ':', color = 'lightgray', lw = 1.2 )
          ax1.set_ylabel('Data size', fontsize= self.font)
          ax1.legend(loc = 'upper right', fontsize= self.font -3)
      

          ratios = 100* (np.array(out) / ( np.array(corr) + np.array(out) ))
          ax2 =fig.add_subplot(212, sharex= ax1)
          ax2.plot(datums, ratios ,color = 'blue' )
          ax2.set_ylim(0,50)
          ax2.yaxis.set_major_locator(plt.MaxNLocator(5))
          ax2.grid(linestyle= ':', color = 'lightgray', lw = 1.2 )
          ax2.set_ylabel('% Outliers', fontsize = self.font -3 )
          ax2.set_xlabel('Dates', fontsize = self.font)
          plt.savefig('plots/outliers/outliers_timeseries_days_' + str(N) + '_hour_' + self.hour + '_' + self.var + '_anp_' + str(self.an_p) + '_fgp_' + str(self.fg_p) + '.pdf'  , bbox_inches = 'tight')
          plt.close()

     def outliers_example(self, corr='', out='', date='', N= '', lower= '', upper= '', median = ''):
          plt.title( self.var + ' for (an_dep,fg_dep)=(' + str(self.an_p)+ ',' + str(self.fg_p) + ') averaged on N='+str(N) + ' days' , y=1.03)
          #plt.plot(date, corr, label = 'True Values', color = 'cyan' )
          #plt.plot(date, out, label = 'Outliers', color ='lime' )

          plt.scatter(date, corr, label = 'True Values', color = 'cyan' , s = 3)
          plt.scatter(date, out, label = 'Outliers', color ='lime'      , s = 3)
          X= [min(date), max(date)]
          plt.plot(X, [lower,lower]  , label = 'Lower' , color ='blue'  , ls = '--' )
          plt.plot(X, [upper,upper]  , label = 'Upper' , color ='red'   , ls = '--' )
          plt.plot(X, [median,median], label = 'Median', color ='black' , ls = '--' )


          plt.legend(fontsize = self.font-6 , loc = 'upper right', ncol = 2)
          plt.grid(linestyle= ':', color = 'lightgray', lw = 1.2 )

          plt.ylabel('Error ' + self.var_dics[self.var]['units'] , fontsize= self.font)

          plt.xlabel('Date', fontsize = self.font)
          plt.xticks(rotation=45)
          Y_max = max ( max(out) +1 , max(corr)+1  )
          Y_min = min ( min(out) -1 , min(corr)-1  )
         
          plt.xlim(min(date)-1/365 , max(date)+1/365 )
      
          if not (np.isnan(Y_max) and np.isnan(Y_min) ): plt.ylim(Y_min, Y_max )
          plt.savefig('plots/outliers/zoom/outliers_' + str(N) + '_date_' + str(min(date)) + '_hour_' + self.hour + '_' + self.var + '_anp_' + str(self.an_p) + '_fgp_' + str(self.fg_p) + '.pdf'  , bbox_inches = 'tight')
          plt.close()

     def time_series(self, means='', datums = '', labels = '', colors = '', interval=50):
 
          plt.title('Time series for ' + self.var + ' for (an_dep,fg_dep)=(' + str(self.an_p)+ ',' + str(self.fg_p) + ')' , y=1.03)
          for (m,d,l,c) in zip(means,datums,labels,colors):
               plt.plot( d[:(len(m))][::interval], m[::interval], color =c, label =l)
          plt.legend(fontsize = self.font-3 , loc = 'upper left')
          plt.ylim(0, 1.2)
          plt.grid(linestyle= ':', color = 'lightgray', lw = 1.2 )
          plt.ylabel('Error ' + self.var_dics[self.var]['units'],  fontsize= self.font) 
          plt.xlabel('Dates', fontsize = self.font)
          plt.savefig('plots/series/timeseries_hour_' + self.hour + '_' + self.var + '_anp_' + str(self.an_p) + '_fgp_' + str(self.fg_p) + '.pdf'  , bbox_inches = 'tight')
          plt.close()


""" Dirs, definitions, select datasets """
cov_file = 'covariance_matrices.npy'
variables = ['temp', 'uwind','vwind','speed','direction']
stations = ['Lindenberg']

netCDF = netCDF()
Cov = Covariance(netCDF)
Plot = Plotter()
Plot.initialize_dirs()

matrices = np.load(cov_file,  allow_pickle = True).item()
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
              for m_index in [500, 5000]: # select here the indexes of the matrices 
                   print('Producing the plots for', m_index)
                   averaged = Cov.average_matrix(matrices_list= lista, N=365, matrix_index= m_index) # considering the m_indexth matrix, and the 365 matrices that follows
                   plotting = Plot.cov_plot(averaged, station=s , hour = h, date= datums[m_index] , averaged = True )
'''



variables = ['temp','uwind','vwind','speed','direction']
variables = ['temp','uwind','vwind', 'speed']

#variables = ['temp']
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

               plevels_i = [0,8,15]
               plevels_j = [0,8,15]  # 11 is 500 hPa

               #plevels_i = [11]
               #plevels_j = [11]  # 11 is 500 hPa
               #for i in plevels_i:
               #     for j in plevels_j:
               for i,j in zip(plevels_i,plevels_j): # no off-diagonal elements

                    Plot.plot_prop(var = v, fg_p = i, an_p = j , hour = str(h))

                    print("*** Processing the i,j entries: ", i , j )
                    values = Cov.select_ijentry(matrices = matrices_dic[str(h)], i = i , j = j)

                    """ The errors are the sqrt of the covariances; note that you must take the abs since they might be negative """
                    values = [ abs(np.sqrt(v)) for v in values ]

                    '''
                    """ Outliers analysis """
                    for N in [30,60,180,365]:
                         print ("*** Outliers analysis  for ", v , ' levels ', i,j )
                         correct, excl = [], []

                         for d in datums:
                              data = values[datums.index(d):datums.index(d)+N]
                              values_cleaned, outliers, lower, upper, median = Cov.remove_outliers(data, min= 0.25, max= 0.75, cut= 1.5 )
                              correct.append( len([x for x in values_cleaned if not np.isnan(x)] ) )
                              excl.append( len([y for y in outliers if not np.isnan(y)] ) )
                              if datums.index(d) in [700,7000]:
                                   index = datums.index(d)
                                   date = datums[index:index+N]
                                   print ('d,corr,out', data, values_cleaned, outliers, date)
                                   Plot.outliers_example(corr=values_cleaned, out= outliers, date= date, N= N, lower= lower , upper= upper , median = median)

                         dat = datums[:len(correct)]
                         Plot.outliers_series(corr= correct, out= excl, datums= dat, N= N, interval= 50)
                    '''

                    values_cleaned, outliers, lower, upper, median = Cov.remove_outliers(values, min= 0.25, max= 0.75, cut= 1.5 )
 
                    """ Plotting outliers total distribution """                    
                    Plot.outliers(outliers = np.sqrt(outliers), cleaned = np.sqrt(values_cleaned), lower = lower, upper = upper, median = median, bins = 200 ) 
                                         
                    means, dates = [], []
                    for n in [30,60,90,180,365]:
                         runningmean, date  = Cov.running_mean_old(data= values_cleaned, n= n, datums= datums)
                         means.append(runningmean)
                         dates.append(date)
                         
                    C = ['yellow', 'orange', 'lime', 'cyan', 'slateblue']
                    L = ['Desroziers (1m)', 'Desroziers (2m)', 'Desroziers (3m)', 'Desroziers (6m)', 'Desroziers (1y)']
                    C = ['steelblue', 'orange', 'forestgreen', 'red', 'slateblue'] # color similar to Michi
     
                    """ Plotting the error time series """
                    #Plot.plot_prop(var = v, fg_p = i, an_p = j , hour = str(h))                     
                    Plot.time_series(means= means, datums= dates, labels= L, colors= C ,interval= 25)

                    """ Plotting the error distributions """
                    bins = 30
                    Plot.histo(X= means, colors = C, labels = L, bins = bins)
               





