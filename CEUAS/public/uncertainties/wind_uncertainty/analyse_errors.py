""" Module for analysisng the covariance matrices
    and extract the Desroziers errors 

    Author: Ambrogi Federico, federico.ambrogi@univie.ac.at
"""


from wind_utils import *


""" Running mode: analyse the covariance matrices and/or the errors """
parser = argparse.ArgumentParser(description="Tool to analyse and create the plots for the covariance matrices and errors")
parser.add_argument('--cov' , '-c', 
                    help="Optional: select True to plot covariance matrices" ,
                    default = False,
                    type = bool)

parser.add_argument('--errors' , '-e',
                    help="Optional: select True to plot time series and distributions of errors " ,
                    default = False,
                    type = bool)


args = parser.parse_args()
plot_cov      = bool(args.cov)
plot_errors   = bool(args.errors)
  
if (not plot_cov and not plot_errors ):
     print(red + '*** ERROR: no plot selected. Please select at least one plot to create\n' + cend )
     sys.exit()

""" Dirs, definitions, select datasets """
cov_file = 'data/covariance_matrices.npy'
variables = ['temp', 'uwind','vwind','speed','direction','rh']
stations = ['Lindenberg']

netCDF = DataHandler()
Cov = Covariance(netCDF)
Plot = Plotter()
Plot.initialize_dirs()

matrices = np.load(cov_file,  allow_pickle = True).item()
print('The matrices are', matrices.keys())
print('The matrices are', matrices['Lindenberg'].keys())

gauss_fit = {}

""" Selecting the time averages """
time_averages = [30,60,90,180,365]


""" For pretty plots """
labels = ['Desroziers (1m)', 'Desroziers (2m)', 'Desroziers (3m)', 'Desroziers (6m)', 'Desroziers (1y)']
colors = ['steelblue', 'orange', 'forestgreen', 'red', 'slateblue'] # color similar to Michi







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
                    """ Creating the covariance matrices plots """
                    lista = matrices[s][v][str(h)]
                    print(blue + '*** Plotting the covariance matrices\n' + cend)
                    for m_index in [500, 1000, 2000, 5000, 8000]: # select here the indexes of the matrices
                         print('   Processing matrix with index: ', str(m_index) , ' for hour: ', str(h) )                                                                                                                                                                   
                         averaged = Cov.average_matrix(matrices_list= lista, N=365, matrix_index= m_index) # considering the m_indexth matrix, and the 365 matrices that follows                                         
                         plotting = Plot.cov_plot(averaged, station=s , hour = h, date= datums[m_index] , averaged = True )                                                    
               
               
               if plot_errors == True:           
                    """ Errors analysis """
                    plevels_i = range(16)
                    plevels_j = range(16)


                    print(yellow+ '*** Plotting the error analysis \n' + cend )

                    for i,j in zip (plevels_i,plevels_j):
                         
                         print('Extracting errors for pressure levels: ', i , j )
                         means, dates = [], []  

                         Plot.plot_prop(var = v, fg_p = i, an_p = j , hour = str(h)) # initialize global plot properties

                         values = Cov.select_ijentry(matrices = matrices_dic[str(h)], i = i , j = j)
                         values_cleaned, outliers, lower, upper, median = netCDF.remove_outliers(values, cut= 1.5 )
                         #values_cleaned = list(np.sqrt(values_cleaned) )
                         values_cleaned = values
                         gauss_fit[v][h][j]={}

                         for n in time_averages : # time intervals defined in the list above
                     
                              print('Averaging over time period: ', n )
                              gauss_fit[v][h][j][n]={}

                              
                              runningmean, date  = Cov.running_mean_old(data= values_cleaned, n= n, datums= datums)
                              means.append(runningmean)
                              dates.append(date)

                              mean,std = Cov.gauss_fit(runningmean)
                              gauss_fit[v][h][j][n]['mean'] = mean
                              gauss_fit[v][h][j][n]['std'] = std
                              
                              """ Plotting the error time series """

                         Plot.time_series(means= means, datums= dates, labels= labels, colors= colors ,interval= 25, station=s)
                              
                         Plot.err_histo(X= means, colors= colors, labels= labels, bins= 30, station=s)
               


""" Extract mean and standard deviation of the distributions of errors """
if not os.path.isfile('data/means_std.npy'):
     np.save('data/means_std', gauss_fit) 




