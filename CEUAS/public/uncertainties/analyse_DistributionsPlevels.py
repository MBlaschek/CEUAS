""" Plotting the errors distributions for different pressure levels 

    Author: Ambrogi Federico, federico.ambrogi@univie.ac.at
   
"""


from uncertainties_utils import *

""" Dirs, definitions, select datasets """
cov_file = 'data/covariance_matrices.npy'
variables = ['temp', 'uwind','vwind','speed','direction', 'rh']
variables = ['temp']
stations = ['Lindenberg']

netCDF = DataHandler()
Cov = Covariance(netCDF)
Plot = Plotter()
Plot.initialize_dirs()

matrices = np.load(cov_file,  allow_pickle = True).item()
print('The matrices are', matrices.keys())
print('The matrices are', matrices['Lindenberg'].keys())



""" pressure levels to be considered """

LEVELS = [0,3,5,8,11,13,15]



res = {}

'''
""" Plotting the time series of the running mean """
for s in stations:
     print ('***** Processing the station: ', s , '\n')
     datums = matrices[s]['datums']
     for v in variables:
          res[v] = {}

          print ('**** Variable: ', v )
          Plot.plot_prop(var=v)                                                                                                                                                                                                                                                          
          matrices_dic = matrices[s][v]
          for h in [1]:

                    plevels_i = LEVELS
                    plevels_j = LEVELS   # 11 is 500 hPa

                    for i,j in zip(plevels_i,plevels_j): # no off-diagonal elements

                         res[v][i] = {}
                         Plot.plot_prop(var = v, fg_p = i, an_p = j , hour = str(h))

                         print("*** Processing the i,j entries: ", i , j )
                         values = Cov.select_ijentry(matrices = matrices_dic[str(h)], i = i , j = j)
                         values_cleaned, outliers, lower, upper, median = netCDF.remove_outliers(values, cut= 1.5 )

                         """ The errors are the sqrt of the covariances; note that you must take the abs since they might be negative """
                         #values_cleaned = [ abs(np.sqrt(v)) for v in values ]

                                           
                         means, dates = [], []
                         for n in [60,365]:
                                   res[v][i][n] = {}
                                   res[v][i][n]['means'] = []
                               
                                   runningmean, date  = Cov.running_mean_old(data= values_cleaned, n= n, datums= datums)
                                   res[v][i][n]['means'].append(runningmean)
                            

print('the results is', res)
np.save('results_for_distributions',res)
'''

res = np.load('data/plevels_distributions.npy', allow_pickle = True).item()


### PLOTTING
colors = ['blue','slateblue','lime','gold','orange','magenta','cyan']
labels_f = [10,20,30,50,70,100,150,200,250,300,400,500, 700,850,925, 1000]
labels = [10,      50,   100,        250,        500, 700,         1000] 

bins = 25
FONT = 13
range = [0.,3.]
dic = {'uwind':'[m/s]' , 'vwind':'[m/s]' , 'speed':'[m/s]', 'temp':'[K]' , 'direction':'[Degree]', 'rh':''}

variables = ['temp','speed']
h = 1
for v in variables:
     for c,i in zip(colors, LEVELS ):
          data_60  = res[v][h][i][60]['means']
          data_365 = res[v][h][i][365]['means'] 
          plt.hist( data_365, histtype = 'stepfilled', range = range, bins = bins , color= c , alpha = 0.2, density = True)
          plt.hist( data_365, histtype = 'step', range = range ,bins = bins , color= c , ls = '-', density = True, alpha=0.6)


     v_n = v.replace('temp','Temperature').replace('uwind','Wind u-component').replace('vwind','Wind v-component')

     for c,l in zip (colors, labels):
           plt.plot([-100,-50], [-100,-50], color = c, label = str(l) + ' [hPa]' )

     plt.title(v_n + ' errors comparison for different pressure levels', fontsize = FONT, y = 1.03)
     plt.text(0.2 , 9.2 , 'Desroziers 1 year ', fontsize = FONT )

     plt.legend(loc = 'upper right', fontsize = FONT-4)
     plt.xlabel('Error ' + dic[v], fontsize = FONT)
     plt.ylabel('Normalized Counts', fontsize = FONT)
     plt.xlim(0.,3)
     plt.ylim(0,10)
     plt.grid(linestyle= ':', color = 'lightgray', lw = 1.2 )
     plt.savefig('plots/plevels_histo/'+ v + '_Errors_Comparison_365.png', bbox_inches = 'tight')
     plt.close()

   
     for c,i in zip(colors, LEVELS ):
          data_60  = res[v][h][i][60]['means']
          data_365 = res[v][h][i][365]['means']
          plt.hist( data_60, histtype = 'stepfilled', bins = bins , color= c , alpha = 0.2, density = True)
          plt.hist( data_60, histtype = 'step', bins = bins , color= c , ls = '-', density = True, alpha=0.6)

     for c,l in zip (colors, labels):
           plt.plot([-100,-50], [-100,-50], color = c, label = str(l) + ' [hPa]' )

     plt.text(0.2 , 9.2 , 'Desroziers 2 months', fontsize = FONT )
     v_n = v.replace('temp','Temperature').replace('uwind','Wind u-component').replace('vwind','Wind v-component')
     plt.title(v_n + ' errors comparison for different pressure levels', fontsize = FONT, y = 1.03)
     plt.legend(loc = 'upper right', fontsize = FONT-4)
     plt.xlabel('Error ' + dic[v], fontsize = FONT)
     plt.ylabel('Normalized Counts', fontsize = FONT)
     plt.xlim(0., 3)
     plt.ylim(0, 10)
     plt.grid(linestyle= ':', color = 'lightgray', lw = 1.2 )
     plt.savefig('plots/plevels_histo/'+ v + '_Errors_Comparison_60.png', bbox_inches = 'tight')
     plt.close()























