B""" Module for extracting the covariance matrix 

    Author: Ambrogi Federico, federico.ambrogi@univie.ac.at

When calculating the error,
to check: first calculate the average then take the sqrt of the covariances or 
first calculate the sqrt then calculate the average ?
   
"""


from wind_utils import *






""" Running mode: forcing the creation of the numpy files """
parser = argparse.ArgumentParser(description="Tool to analyse and create the plots for the covariance matrices and errors")
parser.add_argument('--cov' , '-c', 
                    help="Optional: select True to plot covariance matrices" ,
                    default = False,
                    type = bool)

parser.add_argument('--errors' , '-e',
                    help="Optional: select True to plot time series and distributions of errors " ,
                    default = False,
                    type = bool)

parser.add_argument('--outliers' , '-o',
                    help="Optional: select True to analyse outlier values of the errors" ,
                    default = False,
                    type = bool)

args = parser.parse_args()
plot_cov      = bool(args.cov)
plot_errors   = bool(args.errors)
plot_outliers = bool(args.outliers)
  
if (not plot_cov and not plot_errors and not plot_outliers):
     print(red + '*** ERROR: no plot selected. Please select at least one plot to create\n' + cend )
     sys.exit()

""" Dirs, definitions, select datasets """
cov_file = 'data/covariance_matrices.npy'
variables = ['temp', 'uwind','vwind','speed','direction']
stations = ['Lindenberg']

netCDF = DataHandler()
Cov = Covariance(netCDF)
Plot = Plotter()
Plot.initialize_dirs()

matrices = np.load(cov_file,  allow_pickle = True).item()
print('The matrices are', matrices.keys())
print('The matrices are', matrices['Lindenberg'].keys())

#variables = ['temp','uwind','vwind','speed','direction']
#variables = ['temp','uwind','vwind', 'speed']

variables = ['temp','uwind','vwind', 'speed','direction']

gauss_fit = {}
""" Plotting the time series of the running mean """
for s in stations:
     print ('***** Processing the station: ', s , '\n')
     datums = matrices[s]['datums']
     for v in variables:
          print ('**** Variable: ', v )
          Plot.plot_prop(var=v)                                                                                                                                                                                                                                                          
          matrices_dic = matrices[s][v]
          gauss_fit[v] = {}

          for h in [0,1]:
               gauss_fit[v][h]= {}
               print ('*** Hour: ', h )



               if plot_cov == True:
                    lista = matrices[s][v][str(h)]
                    print(blue + '*** Plotting the covariance matrices\n' + cend)
                    for m_index in [500, 1000, 2000, 5000, 8000]: # select here the indexes of the matrices
                         print('   Processing matrix number: ', str(m_index) )                                                                                                                                                                      
                         averaged = Cov.average_matrix(matrices_list= lista, N=365, matrix_index= m_index) # considering the m_indexth matrix, and the 365 matrices that follows                                                                                  
                         plotting = Plot.cov_plot(averaged, station=s , hour = h, date= datums[m_index] , averaged = True )                                                    


               if plot_errors == True:
               
                    """ Looping over the pressure levels """
                    plevels_i = range(16)
                    plevels_j = range(16)
                    

                    

                    #plevels_i = [0,8,11,15]
                    #plevels_j = [0,8,11,15]  # 11 is 500 hPa

                    for i,j in zip(plevels_i,plevels_j): # no off-diagonal elements

                         gauss_fit[v][h][j]={}

                         Plot.plot_prop(var = v, fg_p = i, an_p = j , hour = str(h))

                         print("*** Processing the i,j entries: ", i , j )
                         values = Cov.select_ijentry(matrices = matrices_dic[str(h)], i = i , j = j)

                         """ The errors are the sqrt of the covariances; note that you must take the abs since they might be negative """
                         values = [ abs(np.sqrt(v)) for v in values ]
                         values_cleaned, outliers, lower, upper, median = netCDF.remove_outliers(values, cut= 1.5 )

                         if plot_outliers:
                              """ Outliers analysis """                    
                              for N in [30,60,180,365]:
                                   print (green + '*** Plotting the outliers analysis\n' + cend)
                                   correct, excl = [], []

                                   for d in datums:
                                        data = values[datums.index(d):datums.index(d)+N]
                                        values_cleaned, outliers, lower, upper, median = Cov.remove_outliers(data, cut= 1.5 )
                                        correct.append( len([x for x in values_cleaned if not np.isnan(x)] ) )
                                        excl.append( len([y for y in outliers if not np.isnan(y)] ) )
                                        if datums.index(d) in [700,7000]:
                                             index = datums.index(d)
                                             date = datums[index:index+N]
                                             print ('d,corr,out', data, values_cleaned, outliers, date)
                                             Plot.outliers_example(corr=values_cleaned, out= outliers, date= date, N= N, lower= lower , upper= upper , median = median)

                                        dat = datums[:len(correct)]
                                        Plot.outliers_series(corr= correct, out= excl, datums= dat, N= N, interval= 50)
                    
 
                              """ Plotting outliers total distribution """                    
                              Plot.outliers(outliers = np.sqrt(outliers), cleaned = np.sqrt(values_cleaned), lower = lower, upper = upper, median = median, bins = 200 ) 


                         if plot_errors:             
                              """ Errors analysis """
                              print(yellow+ '*** Plotting the error analysis \n' + cend )
                              means, dates = [], []
                              for n in [30,60,90,180,365]:
             
                                   gauss_fit[v][h][j][n]={}

                                   runningmean, date  = Cov.running_mean_old(data= values_cleaned, n= n, datums= datums)
                                   means.append(runningmean)
                                   dates.append(date)

                                   mean,std = Cov.gauss_fit(runningmean)
                                   gauss_fit[v][h][j][n]['mean'] = mean
                                   gauss_fit[v][h][j][n]['std'] = std

                         
                              C = ['yellow', 'orange', 'lime', 'cyan', 'slateblue']
                              L = ['Desroziers (1m)', 'Desroziers (2m)', 'Desroziers (3m)', 'Desroziers (6m)', 'Desroziers (1y)']
                              C = ['steelblue', 'orange', 'forestgreen', 'red', 'slateblue'] # color similar to Michi
     
                              """ Plotting the error time series """
                              Plot.time_series(means= means, datums= dates, labels= L, colors= C ,interval= 25, station=s)
                              
                              bins = 30
                              Plot.err_histo(X= means, colors = C, labels = L, bins = bins, station=s)
               



np.save('Gauss_analysis', gauss_fit)




