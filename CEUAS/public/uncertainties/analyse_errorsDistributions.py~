""" Module for plotting the erorrs for each pressure level 

    Author: Ambrogi Federico, federico.ambrogi@univie.ac.at

    For each pressure level it plots the mean and standard deviation 
    of the distributions of the Desroziers errors, averaged over 1 month and 1 year

"""


import os,sys
import netCDF4
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
from wind_utils import *


""" Dirs, definitions, select datasets """
#cov_file = 'data/covariance_matrices.npy'
variables = ['temp', 'uwind','vwind','speed','direction']

stations = ['Lindenberg']

#netCDF = netCDF()
#Cov = Covariance(netCDF)
#Plot = Plotter()
#Plot.initialize_dirs()

#matrices = np.load(cov_file,  allow_pickle = True).item()
#print('The matrices are', matrices.keys())
#print('The matrices are', matrices['Lindenberg'].keys())





variables = ['temp','uwind','vwind','speed','direction','rh']

LEVELS = [0,3,5,8,11,13,15]



#mv meres = np.load('data/mean_std_distributions.npy', allow_pickle = True).item()


### PLOTTING
colors = ['blue','slateblue','lime','gold','orange','magenta','cyan']
labels_f = [10,20,30,50,70,100,150,200,250,300,400,500, 700,850,925, 1000]
labels = [10,      50,   100,        250,        500, 700,         1000] 

bins = 50
FONT = 13
#range = [0.4,2]
dic = {'uwind':'[m/s]' , 'vwind':'[m/s]' , 'speed':'[m/s]', 'temp':'[K]' , 'direction':'[Degree]' , 'rh':''}


""" Preparing the plots directory """
path = 'plots/mean_std_plevels'
if not os.path.os.path.exists(path):
     os.makedirs(path)

""" Loading the numpy file """
try:
     gauss = np.load('data/means_std.npy', allow_pickle=True).item()
except:
     print('The file data/means_std.npy cannot be found!')
     sys.exit()


""" Plotting """
'''
for v in variables:
    h0 = gauss[v][0]
    h1 = gauss[v][1]
    for m,c in zip([30, 60, 365], ['blue','cyan','lime']):

         fig,ax = plt.subplots(figsize=(5,7))
         FONT = 13
       
         means0, means1, stds0, stds1 = [],[],[], []
      
         levels = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]


         for p in levels:
             means0.append( h0[p][m]['mean'] )
             stds0.append ( h0[p][m]['std' ] )
             means1.append( h1[p][m]['mean'] )
             stds1.append ( h1[p][m]['std' ] )
      
         means0 = list(reversed(means0) )
         means1 = list(reversed(means1) )
         stds0  = list(reversed(stds0)  )
         stds1  = list(reversed(stds1)  )

         lab_f = list(reversed(labels_f))

         plt.yticks(labels_f, labels_f, fontsize = 7)
         plt.ylim(1030,0)
         #ax.set_yticklabels(labels = labels_f, fontsize = FONT)

         M = str(m).replace('30','1 month').replace('60', '2 months').replace('365', '1 year')
         
         plt.title(v + ' error mean values for a Gaussian distribution, Desroziers ' + M, y=1.03)
         plt.plot(means0, lab_f, color  = c , label = 'Z00', ls = ":")
         plt.plot(means1, lab_f, color  = c , label = 'Z12', ls = '--')

         for m1,s1,p in zip(means1,stds1,lab_f):
             x , y  = [m1-s1/2 , m1+s1/2] , [p , p]
             plt.plot(x,y, color = c, ls = '--')

         for m0,s0,p in zip(means0,stds0,lab_f):
              x , y  = [m0-s0/2 , m0+s0/2] , [p , p]
              plt.plot(x,y, color = c, ls = ':')

         plt.xlabel('Error ' + dic[v] , fontsize = FONT )
         plt.ylabel('Pressure [hPa]'  , fontsize = FONT)

         plt.grid(linestyle= ':', color = 'lightgray', lw = 1.2 )
         plt.legend(loc = 'best')        
         plt.savefig(path + '/' + v + '_averages_gaussian_'+str(m)+'.png', bbox_inches = 'tight' , dpi = 180 ) 
         plt.close()
'''


#### new version, new plots

for v in variables:
    print('Plotting the errors vs peressure level for the variable: ',  v )
    h0 = gauss[v][0]
    h1 = gauss[v][1]
    for a in [0]:

         fig,ax = plt.subplots(figsize=(5,7))
         FONT = 13

         means0_30, means1_30, stds0_30, stds1_30 = [],[],[], []
         means0_365, means1_365, stds0_365, stds1_365 = [],[],[], []

         levels = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]

         for p in levels:
             means0_30.append( h0[p][30]['mean'] )
             stds0_30.append ( h0[p][30]['std' ] )
             means1_30.append( h1[p][30]['mean'] )
             stds1_30.append ( h1[p][30]['std' ] )


             means0_365.append( h0[p][365]['mean'] )
             stds0_365.append ( h0[p][365]['std' ] )
             means1_365.append( h1[p][365]['mean'] )
             stds1_365.append ( h1[p][365]['std' ] )


         means0_30 = list(reversed(means0_30) )
         means1_30 = list(reversed(means1_30) )
         stds0_30  = list(reversed(stds0_30)  )
         stds1_30  = list(reversed(stds1_30)  )

         means0_365 = list(reversed(means0_365) )
         means1_365 = list(reversed(means1_365) )
         stds0_365  = list(reversed(stds0_365)  )
         stds1_365  = list(reversed(stds1_365)  )


         lab_f = list(reversed(labels_f))

         plt.yticks(labels_f, labels_f, fontsize = 7)
         plt.ylim(1030,0)

         plt.title(v+' Error means', y=1.03)

         plt.plot(means0_30,  lab_f, color  = 'blue' , ls = '-' )
         plt.plot(means0_365, lab_f, color  = 'red'  , ls = '-' )

         plt.plot(means1_30,  lab_f, color  = 'blue' ,  ls = ':' )
         plt.plot(means1_365, lab_f, color  = 'red'  ,  ls = ':' )

         plt.plot([1,1.1],[-5,-5], color = 'black' , ls = ':'  , label = '12GTM' )
         plt.plot([1,1.1],[-5,-5], color = 'black' , ls = '-' , label = '00GMT' )
         plt.plot([1,1.1],[-5,-5], color = 'blue'  , ls = '--'  , label = 'Desroziers (1m)' )
         plt.plot([1,1.1],[-5,-5], color = 'red'   , ls = '--' , label = 'Desroziers (1y)' )


         
         for m1,s1,p in zip(means1_30,stds1_30,lab_f):
             x , y  = [m1-s1/2 , m1+s1/2] , [p , p]
             plt.plot(x,y, color = 'blue', ls = ':')

         for m0,s0,p in zip(means0_30,stds0_30,lab_f):
              x , y  = [m0-s0/2 , m0+s0/2] , [p , p]
              plt.plot(x,y, color = 'blue', ls = '-')


         for m0,s0,p in zip(means1_365,stds1_365,lab_f):
              x , y  = [m0-s0/2 , m0+s0/2] , [p , p]
              plt.plot(x,y, color = 'red', ls = ':')

         for m0,s0,p in zip(means0_365,stds0_365,lab_f):
              x , y  = [m0-s0/2 , m0+s0/2] , [p , p]
              plt.plot(x,y, color = 'red', ls = '-')



         
         plt.xlabel('Error ' + dic[v] , fontsize = FONT )
         plt.ylabel('Pressure [hPa]'  , fontsize = FONT)

         plt.grid(linestyle= ':', color = 'lightgray', lw = 1.2 )
         plt.legend(loc = 'best', fontsize = 9)
         plt.savefig(path + '/'+ v + '_all_errors.png', bbox_inches = 'tight' , dpi = 180 )
         plt.close()






























