""" Module for extracting the covariance matrix 

    Author: Ambrogi Federico, federico.ambrogi@univie.ac.at
    
    Return: numpy file called covariance_matrices.npy
            Contains a python dictionary storing the covariance matrices as numpy arrays
            If required, it creates plots to analyse the outliers, stored in plots/outliers

    Usage: python extract_covaraince_analyseOutliers.py [-f True][-o True]
           default = False

"""


#import os,sys
#import netCDF4
#import numpy as np
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pylab as plt
#import argparse
from wind_utils import *


""" Running mode: forcing the creation of the numpy files with [-f True] option , or outliers analysis with [-o True]"""
parser = argparse.ArgumentParser(description="Force the creation/update of the numpy files storing the covariances")
parser.add_argument('--force' , '-f', 
                    help="Optional: select True to force the creation of the numpy files storing the covariances" ,
                    default = False,
                    type = bool)
parser.add_argument('--outliers' , '-o',
                    help="Optional: select True to analyse outlier values of the errors" ,
                    default = False,
                    type = bool)

args = parser.parse_args()
force = args.force
plot_outliers = bool(args.outliers)


""" File that will contain the numpy dictionary.
    If exists, the code will not do anything unless the [force] option is chosen
    Otherwise, it will create the file at the specified location """
cov_file = 'data/covariance_matrices.npy'


""" Path to the the netCDF files"""
base_dir = 'data/'
file_dic = {'temp'     : 'ERA5_1_10393_t.nc'  , 
            'uwind'    : 'ERA5_1_10393_u.nc'  , 
            'vwind'    : 'ERA5_1_10393_v.nc'  , 
            'dp'       : 'ERA5_1_10393_dp.nc' , 
            'rh'       : 'ERA5_1_10393_rh.nc' ,
            'speed'    : 'ERA5_1_10393_speed.nc', 
            'direction': 'ERA5_1_10393_direction.nc'  }

""" Variables and observation station """
variables = ['temp', 'uwind','vwind', 'speed','direction','rh']
stations = ['Lindenberg']


""" Select the time average for the Desroziers method """
time_average = [30,180,365]

""" Select the index of the matrices to analyse (mapping to the datum list) """
matrix_indeces = [700,7000] 

""" Initialising the classes from thewind_utils.py module """
netCDF = DataHandler()
Plot = Plotter()
Plot.initialize_dirs()

Cov = Covariance(netCDF)


# *********************************************
# Extracting the covariance matrices
# ********************************************* 

""" Initialize and empty dictionary of all the matrices """
matrices = {}


""" Extracting the full (cross)covariance matrices, and store in dictionary """
if not os.path.isfile(cov_file) or force == True :
     print(red + '*** Extracting the covariance matrices and storing in a numpy dictionary' + cend)
     for s in stations:
          print(blue + '** Processing the Station: ', s  + cend)
          matrices[s]= {}
          for v in variables:
               print(green + '*Analysing the variable: ', v , '\n' + cend)
               matrices[s][v]={}
               input_file = base_dir + file_dic[v]
               data = netCDF.read_data(file = input_file) #  Loading the necCDF file

               andep = np.zeros([2,16,len(data['an_dep']['data'][1,1,:])]) # empty matrix
               fgdep = np.zeros([2,16,len(data['an_dep']['data'][1,1,:])])

               datums = data['fg_dep']['datum']

               for x in [0,1]: #  loop over the 00GMT(0) and 12GMT(1)   
                    for y in range(16): #  loop over the 16 standard pressure levels
                         print('Removing the outlier values for hour:', x , ' plevel:' , y)

                         cleaned_an, outliers, lower, upper, median = netCDF.remove_outliers(data['an_dep']['data'][x,y,:], cut= 1.5 )  # removing outliers
                         cleaned_fg, outliers, lower, upper, median = netCDF.remove_outliers(data['fg_dep']['data'][x,y,:], cut= 1.5 )      
                         andep[x,y,:]=cleaned_an
                         fgdep[x,y,:]=cleaned_fg 


                         if plot_outliers == True:
                              Plot.plot_prop(var= v, fg_p= str(y), an_p= str(y), hour= str(x))
                              """ Performs the outlier analysis ona subset of an_dep and fg_dep,  produce plots """
                              for N in time_average: # [30,60,180,365] # select the time period for the Desroziers average 
                                   for d in ['an_dep', 'fg_dep']:
                                        print (green + '*** Plotting the ' , d , ' outliers analysis taking ndays = ' + str(N) + ' \n' + cend)
                                        correct, excl , correct_s, excl_s= [], [], [], []
 
                                        for index in matrix_indeces: 

                                             Data = data[d]['data'][x,y,:][index:index+N] # selecting only an arbitrary subset of the data (i.e. for indexes in above list)                                             
                                             dates = datums[index:index+N]     

                                             values_cleaned, outliers, lower, upper, median = netCDF.remove_outliers(Data, cut= 1.5 )
                                             values_cleaned_s, outliers_s, lower_s, upper_s, median_s = netCDF.remove_outliers(Data, cut= 1.5 , skewed= True)

                                             correct.append( len([x for x in values_cleaned if not np.isnan(x)] ) )
                                             excl.append( len([y for y in outliers if not np.isnan(y)] ) )

                                             correct_s.append( len([x for x in values_cleaned_s if not np.isnan(x)] ) )
                                             excl_s.append( len([y for y in outliers_s if not np.isnan(y)] ) )

                                             Plot.outliers_example(corr=values_cleaned, out= outliers, what = d, date= dates, N= N, lower= lower , upper= upper , median= median, lower_s= lower_s , upper_s= upper_s, station = s)
                              print('Finished plotting outlier examples\n')
                 
               print(blue + " Extracting the covariances """ + cend)

               all_matrices = Cov.extract_matrix_list(an_dep=andep, fg_dep=fgdep, datums=datums, hours=[0,1])
               matrices[s][v]['0']   = all_matrices['0']
               matrices[s][v]['1']   = all_matrices['1']
               matrices[s]['datums'] = datums
     np.save('data/covariance_matrices', matrices)
     print('Saved the covariance matrices in data/covariance_matrices.npy\n')

else:
     print(red + '*** Loading the covariance_matrices.npy dictionary' + cend) 
     matrices = np.load(cov_file,  allow_pickle = True).item()


