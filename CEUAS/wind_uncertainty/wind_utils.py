""" Module containing utilities for the wind uncertainty analysis

    Author: Ambrogi Federico, federico.ambrogi@univie.ac.at

    class DataHandler
    class Covariance
    class Plotter
"""



import os,sys
import netCDF4
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import argparse
import cftime 
from scipy.stats import norm  



""" Some colors for pretty printout """ 
red    = '\033[91m' 
cend   = '\033[0m'
blue   = '\033[34m'
green  = '\033[92m'
yellow = '\033[33m'

  
class DataHandler:
     """ Class containig the basic information of the netCDF files.
         All the netCDF files must contain the observations, and analysis and background departures """
     def __init__(self, file=""):
          self.file = file

     def remove_outliers(self, data, min_p= 25, max_p= 75, cut= '', skewed= False):
          """ Finds outliers, and replace them with np.nan (to keep vector of same length)                                                                                                                                                                                              
              Returns: cleaned             = list of values without outliers (filled with nans replacing outliers)                                                                                                                                                                   
                       outliers            = list of outlier values                                                                                                                                                                                                               
                       lower,upper, median = outliers delimiter and median values """
          data_c = [ d for d in data if d ]
          q25, q75 = np.nanpercentile(data_c, min_p), np.nanpercentile(data_c, max_p)
          cut_off = (q75 - q25) * cut
          lower, upper = q25-cut_off, q75+cut_off

          if skewed==True:
               q50 = np.nanpercentile(data_c, 50)
               lower , upper = q25-(q50-q25)*cut ,  q75+(q75-q50)*cut

          median = np.nanmedian(data_c)
          cleaned, outliers = [],[]

          for d in np.asarray(data):
              if d >= lower and d <= upper:
                   cleaned.append(d)
                   outliers.append(np.nan)
              elif np.isnan(d):
                   cleaned.append(np.nan)
                   outliers.append(np.nan)
              else:
                   cleaned.append(np.nan)
                   outliers.append(d)
          return cleaned, outliers, lower, upper, median

     def read_data(self, file = '', var = ['fg_dep', 'an_dep']):
         
         data_loaded = {}  # Dictionary containing the extracted information                                                                                                                                                                                                              
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

     def running_mean_old(self, data = '', n='', datums=''):
          """ Returns the running mean of x, averaged over N  
              Data is reduced by removing outliers """
          not_nans         = [ x for x in data if not np.isnan(x) and x > -1000 ] # FF maybe not needed anymopre the > -1000 (from old version?)
          not_nans_indexes = [ data.index(x) for x in data if not np.isnan(x) ]
          datums_nans      = [ datums[i] for i in not_nans_indexes ] # extracting not nans values and corresponding index in datums (for plotting) 
          cumsum = np.cumsum(np.insert(not_nans, 0, 0))
          means = (cumsum[n:] - cumsum[:-n]) / float(n)
          #means = np.sqrt(np.absolute(means)) # CHECK !!!
          return means , datums_nans
    
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

     def gauss_fit(seld, data=''):
          """ Fits the data assuming a normal distributions, returns mean and standard deviation """
          mean, std = norm.fit(data)
          return mean, std

     def average_matrix(self,matrices_list='', N='', matrix_index=""):
         """ Given the list of cov. matrices and the index of the chosen matrix, returns an averaged matrix.
         The average is calculated using the N following matrices. 
         Only the the entries with actual values are used, nans values are skipped """

         averaged = np.empty([16, 16], dtype = float) # initialize an empty 16*16 matrix (for the 16 standard pressure levels)                                                                                                                            
         new_mlist = matrices_list[matrix_index:matrix_index+N]  # selecting the subset of matrices from index=matrix_index to matrix_index+N
         for x in range(16):
              for y in range(16):
                   valueslist = self.select_ijentry(matrices=new_mlist , i=x , j=y)
                   cleaned, outliers, lower, upper, median = self.data.remove_outliers(valueslist, cut= 1.5) 
                   not_nan = [ x for x in cleaned if not np.isnan(x) ]
                   N = len(not_nan)
                   try:
                        averaged[x,y] = sum(not_nan)/N  
                   except ZeroDivisionError:
                        averaged[x,y] = np.nan
         return averaged
 


class Plotter:

     def __init__(self):
          self.pretty_pressure_dic = { '0':'10' , '1':'20' , '2':'30' , '3':'50' , '4' :'70' , '5':'100',
                                       '6':'150', '7':'200', '8':'250', '9':'300', '10':'400',
                                       '11':'500', '12':'700' , '13':'850', '14':'925', '15':'1000' }

     def plot_prop(self, var='', fg_p = '' , an_p = '' , hour = ''):
         self.font = 15

         self.fg_p = fg_p #  pressure level of the first guess dep.
         self.an_p = an_p #  pressure level of the analysis dep.

         self.pretty_pressure_dic = { '0':'10' , '1':'20' , '2':'30' , '3':'50' , '4' :'70' , '5':'100',
                                      '6':'150', '7':'200', '8':'250', '9':'300', '10':'400',
                                      '11':'500', '12':'700' , '13':'850', '14':'925', '15':'1000' }
                                      
         self.pretty_pressure = [10,20,30,50,70,100,150,200,250,300,400,500,700,850,925,1000]
         self.var = var
         #self.pretty_var = self.var.replace('uwind','Wind u-component').replace('vwind','Wind v-component').replace('temp', 'Temperature').replace('speed','Wind Speed').replace('direction','Wind Direction')
         self.hour = hour
         self.var_dics = { 'temp' : { 'units': '[K]'         , 'name':'Air_Temperature'  , 'x_range': [] , 'y_range': [] } ,
                           'uwind': { 'units': '[m/s]'       , 'name':'Wind_u-component' , 'x_range': [] , 'y_range': [] } ,  
                           'vwind': { 'units': '[m/s]'       , 'name':'Wind_v-component' , 'x_range': [] , 'y_range': [] } ,
                           'dp'   : { 'units': '[K]'         , 'name':'Dew_Point_Temp.'  , 'x_range': [] , 'y_range': [] } ,
                           'rh'   : { 'units': ''            , 'name':'Relative_Hum.'    , 'x_range': [] , 'y_range': [] } ,
                           'speed': { 'units': '[m/s]'       , 'name':'Wind_Speed'       , 'x_range': [] , 'y_range': [] } ,
                           'direction': { 'units': 'Degree'  , 'name':'Wind_Direction'   , 'x_range': [] , 'y_range': [] } ,  }

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

     def date_prettyfier(self, date):
         """ Give the date in text format day/motnh/year from the datum format """
         units = 'days since 1900-01-01 00:00'
         date = date * 365.25
         date = cftime.num2date(date, units)
         pretty_date = str(date.day)+'/'+str(date.month)+'/'+str(date.year-1900)         
         return pretty_date

     def cov_plot(self, matrix, station="", hour = "", date="" , averaged = "" ):
          """ Basic plot for the correlation matrix """
          var = self.var_dics[self.var]['name'] 
          fig,ax = plt.subplots()
          date = self.date_prettyfier(date)
          hour = str(hour).replace('0','00:00').replace('1','12:00')
          if not averaged:
             title    = "Stat: " + station + ', H: '  + hour + ', Date: '    + date + ', ' + var
             filename = 'Cov_'   + station + '_hour_' + hour.replace(':','') + '_date_'    + str(date).replace('/','') + '_' +var
 
          elif averaged :
               title = var.replace('temp','Temp.') + " , Stat: " + station + ', H: ' + str(hour) + ', Date: ' + str(date)
               filename ='Cov_' + station + '_hour_' + str(hour).replace(':','') + '_averaged_' + str(date).replace('/','') + '_' + var 

          plt.title(title.replace('_', ' ' ), y=1.03, fontsize = self.font-2)

          num = len(matrix[0,:])
          Num = range(num)

          vmin, vmax = -3, 3
          if self.var == 'direction': 
               vmin, vmax = -10, 10
          color_map= plt.imshow(matrix, interpolation= 'nearest', cmap = 'RdYlBu', vmin = vmin, vmax = vmax ) # nearest serves for discreete grid  # cmaps blue, seismic            
          plt.ylim(-0.5, 15.5)
          plt.xlim(-0.5, 15.5)
          plt.xticks(Num, Num)
          plt.xlabel('Pressure level an_dep [hPa]', fontsize = self.font-2)
          plt.yticks(Num, Num)
          plt.ylabel('Pressure level fg_dep [hPa]', fontsize = self.font-2)
          ax.set_xticklabels(labels = self.pretty_pressure, fontsize = self.font-4, rotation=45)
          ax.set_yticklabels(labels = self.pretty_pressure, fontsize = self.font-4)

          bar = plt.colorbar()
          bar.ax.set_ylabel("Covariance", fontsize = self.font)
                                                                                                                 
          for i in Num: # creating text labels
               for j in Num:
                    value = '{0:.2f}'.format(matrix[i,j])
                    text = ax.text( j,i, value , ha = 'center' , va = 'center', color = 'black', fontsize = 5)

          if not os.path.isdir('plots/covariances/'+station): os.mkdir('plots/covariances/'+station)
          plt.savefig('plots/covariances/' + station + '/' + filename + '.png', bbox_inches='tight', dpi = 200)
          plt.close()

     def err_histo(self, X="", colors="", labels="" , bins = "", station = ""):
          hour = str(self.hour).replace('0','00:00').replace('1','12:00')

          pressure = self.pretty_pressure_dic[str(self.an_p)]
          var = self.var_dics[self.var]['name']
          plt.title(var + ' errors, Stat: ' + station + ', H: '  + hour + ', P: ' + pressure + ' [hPa]' , y=1.03)

          plt.hist(X, bins, histtype='stepfilled',  stacked = False, color = colors , label = labels , alpha = 0.4 , density = True)

          plt.grid(linestyle= ':', color = 'lightgray', lw = 1.2 )
          plt.legend(loc = 'upper right', fontsize = self.font - 4)
          plt.ylabel('Normalized Counts' , fontsize = self.font)

          min_x , max_x = 100, 0          
          for x in X:
              m = min(x)
              if m < min_x: min_x = m
              M = max(x)
              if M > max_x: max_x = M
          if self.var =='temp':
               y_max = 25
          elif self.var =='direction':
               y_max = 5
          else: y_max = 10

          plt.ylim(0, y_max)
          plt.xlim(min_x - min_x/3. , max_x + max_x/3.)
          plt.xlabel(r'Errors ' + self.var_dics[self.var]['units'] , fontsize= self.font)
          plt.savefig('plots/histo/histo_' + self.var + '_hour_' + self.hour + '_anp_' + str(self.an_p) + '_fgp_' + str(self.fg_p) + '.png',  bbox_inches='tight')
          plt.close()
     '''
     def outliers (self, outliers='', cleaned = '', datums = '', bins = 25, median = '', lower = '', upper = ''):
          """ Plot the cumulative distributions of the number of outliers and true data """
          plt.title('Outlier for ' + self.var + ' for (an_dep,fg_dep)=(' + str(self.an_p)+ ',' + str(self.fg_p) + ')' , y=1.03)

          plt.text(0.1 , 1790 , 'Median [{:.2f}]'.format(median) , fontsize = self.font -3 )
          plt.text(0.1 , 1700 , 'Lower quart. [{:.2f}]'.format(lower) , fontsize = self.font -3 )
          plt.text(0.1 , 1610 , 'Upper quart. [{:.2f}]'.format(upper) , fontsize = self.font -3)

          cleaned = [ c for c in cleaned if not np.isnan(c)]
          colors = ['slateblue','limegreen']
          labels = ['Outliers ['+ str(len(outliers)) + ']', 'Data [' + str(len(cleaned)) + ']' ]
          plt.hist([outliers,cleaned], bins, histtype='stepfilled',  stacked = False, color = colors , label = labels , alpha = 0.7 , density = False)
 
          plt.legend(fontsize = self.font-3 , loc = 'upper right')

          plt.xlim(0, 3)
          plt.ylim(0, 1900)
          plt.grid(linestyle= ':', color = 'lightgray', lw = 1.2 )
          plt.ylabel('Number ', fontsize= self.font)
          plt.xlabel('Cross-Covariance', fontsize = self.font)
          plt.savefig('plots/histo/outliers_hour_' + self.hour + '_' + self.var + '_anp_' + str(self.an_p) + '_fgp_' + str(self.fg_p) + '.png'  , bbox_inches = 'tight')
          plt.close()

     def outliers_series(self, corr='', out='', datums='', N= '', interval= '', flag=''):
          """ Plot the number of outliers and true data vs date """

          #pressure = self.pretty_pressure_dic[str(self.an_p)]
          #var = self.var_dics[self.var]['name']
          #hour = str(self.hour).replace('0','00:00').replace('1','12:00')
          #plt.title(var + ' time series, Stat: ' + station + ', H: '  + hour + ', P: ' + pressure + ' [hPa]' , y=1.03)

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
          plt.savefig('plots/outliers/' + flag + '_outliers_timeseries_days_' + str(N) + '_hour_' + self.hour + '_' + self.var + '_anp_' + str(self.an_p) + '_fgp_' + str(self.fg_p) + '.png'  , bbox_inches = 'tight')
          plt.close()
     '''
     def outliers_example(self, corr= '', out= '', date= '', N= '', lower= '', upper= '', median= '', flag='', upper_s = '', lower_s= '' , station = '', what = ''):

          pressure = self.pretty_pressure_dic[str(self.an_p)]
          var = self.var_dics[self.var]['name']
          hour = str(self.hour).replace('0','00:00').replace('1','12:00')
          plt.title(var + ' ' + what + ' Outliers - Stat: ' + station + ', H: '  + hour + ', P: ' + pressure + ' [hPa]' , y=1.03)

          corr_ = [ n for n in corr if not np.isnan(n) ]
          out_ = [ n for n in out if not np.isnan(n) ]

          num_a = '{:.1f}'.format(len(corr_)/len(out_ + corr_) * 100)
          num_o = '{:.1f}'.format(len(out_)/len(out_ + corr_) * 100)

          plt.scatter(date, corr, label = 'Accepted [' + num_a + '%]', color = 'cyan' , s = 3)
          plt.scatter(date, out, label = 'Outliers [' + num_o + '%]', color ='black'      , s = 3)
          X= [min(date), max(date)]

          plt.plot(X, [lower,lower]  , label = 'Lower' , color ='blue'  , ls = '--' )
          plt.plot(X, [upper,upper]  , label = 'Upper' , color ='red'   , ls = '--' )

          # adding the upper and lower values for skewed distributions
          plt.plot(X, [lower_s, lower_s]  , label = 'Lower Skewed' , color ='blue'  , ls = '-' )
          plt.plot(X, [upper_s, upper_s]  , label = 'Upper Skewed' , color ='red'   , ls = '-' )

          plt.plot(X, [median,median], label = 'Median [' + '{:.1f}'.format(median) + ']', color ='black' , ls = '--' )

          plt.legend(fontsize = self.font-6 , loc = 'upper right', ncol = 2)
          plt.grid(linestyle= ':', color = 'lightgray', lw = 1.2 )

          plt.ylabel('Departure ' + self.var_dics[self.var]['units'] , fontsize= self.font)

          plt.xlabel('Date', fontsize = self.font)
          plt.xticks(rotation=45)

          out_c =  [ n for n in out  if not np.isnan(n)] 
          corr_c = [ n for n in corr if not np.isnan(n)] 

          plt.xlim(min(date)-1/365 , max(date)+1/365 )
     
          plt.ylim(-10, 10)
     
          plt.savefig('plots/outliers/outliers_' + flag + '_' + str(N) + '_date_' + str(min(date)) + '_hour_' + self.hour + '_' + self.var + '_anp_' + str(self.an_p) + '_fgp_' + str(self.fg_p) + '.png'  , bbox_inches = 'tight')
          plt.close()





     def time_series(self, means='', datums = '', labels = '', colors = '', interval=30, station=''):
          pressure = self.pretty_pressure_dic[str(self.an_p)]
          var = self.var_dics[self.var]['name']
          hour = str(self.hour).replace('0','00:00').replace('1','12:00')

          plt.title(var + ' time series, Stat: ' + station + ', H: '  + hour + ', P: ' + pressure + ' [hPa]' , y=1.03)

          for (m,d,l,c) in zip(means,datums, labels, colors):
               plt.plot( d[:(len(m))][::interval], m[::interval], color =c, label =l )
               y_min, y_max = min(m)-min(m), max(m) 
               plt.ylim(y_min - y_min/2. , y_max + y_max/2.  )

          plt.legend(fontsize = self.font-5 , loc = 'upper left')
          plt.grid(linestyle= ':', color = 'lightgray', lw = 1.2 )
          plt.ylabel('Error ' + self.var_dics[self.var]['units'],  fontsize= self.font) 
          plt.xlabel('Dates', fontsize = self.font)
          plt.savefig('plots/series/timeseries_hour_' + self.hour + '_' + self.var + '_anp_' + str(self.an_p) + '_fgp_' + str(self.fg_p) + '.png'  , bbox_inches = 'tight')
          plt.close()

