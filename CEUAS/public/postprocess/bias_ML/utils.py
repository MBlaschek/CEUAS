import os,sys
import h5py as h5
import pandas as pd 
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import ticker

import seaborn as sb

from scipy import stats

from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split


from statsmodels.tsa.vector_ar.var_model import VAR





def get_CUON_data(file, var="85", merged=True, std_plev=True ):
    '''
    Retrieves the data from a merged CUON file, using h5py, and creates a pandas dataframe.
            Parameters:
                    std_plev(Bool):  return only data on standard pressure levels
            
                    merged(Bool): tells if the file is of merged type (otherwise standard harvested netCDF)
                                   NB simple harvested files still to be implemented
                                  
                    var (int or str): variable number according to CDM convention (85: temp, 38: relative hum, 107: wind speed)
            Return:
                    pandas dataframe
    '''
    
    stat = file.split('/')[-1].split("-")[3]
    name = stat + '_' + var + "_" + '.pkl'
                    
    print(name)
    
    if os.path.isfile(name):
        print("---Loading saved pickle file ", name )
        df = pd.read_pickle(name)
        return df
                                     
    else:     
                                     
        f = h5.File(file, 'r')
        var = str(var)

        ind = f["recordindices"][var]
        imin, imax = ind[0], ind[-1]

        z = f["observations_table"]["z_coordinate"][imin:imax] 

        dt = f["observations_table"]["date_time"][imin:imax]
        dt = pd.to_datetime(dt,  unit='s', origin=pd.Timestamp('1900-01-01'))

        obs = f["observations_table"]["observation_value"][imin:imax]
        bias = f["era5fb"]["biascorr@body"][imin:imax]
        fg_dep_adj = f["era5fb"]["fg_depar@body"][imin:imax]

        an_dep_adj = f["era5fb"]["an_depar@body"][imin:imax]


        if var == "85": # bias is available only for temperature 

            obs_adj = obs - bias
            fg_dep = fg_dep_adj + bias
            an_dep = an_dep_adj + bias 
            red_df = pd.DataFrame(  { 'z'            : z,
                                                    'fg_dep_adj'    : fg_dep_adj,
                                                    'fg_dep'           : fg_dep ,    

                                                    'an_dep_adj'    : an_dep_adj,
                                                    'an_dep'        : an_dep,                                

                                                    'bias'       : bias ,
                                                    'obs'        : obs,
                                                    'obs_adj'        : obs_adj,
                                                    'bias_calc'     : fg_dep - fg_dep_adj,


                                                }
                                            )

        red_df["time"] = dt            

        # selecting standard pressure level
        if std_plev:
            print("*** Extracting standard pressure level")
            red_df = red_df.loc[ (red_df['z'] == 1000 )  |  
                                             (red_df['z'] == 2000 ) |
                                             (red_df['z'] == 3000 ) |
                                             (red_df['z'] == 5000 ) |
                                             (red_df['z'] == 7000 ) |
                                             (red_df['z'] == 10000 ) |
                                             (red_df['z'] == 15000 )  |
                                             (red_df['z'] == 20000 )  | 
                                             (red_df['z'] == 25000 ) | 
                                             (red_df['z'] == 30000 ) |
                                             (red_df['z'] == 40000 )  |
                                             (red_df['z']== 50000 ) | 
                                             (red_df['z'] == 70000 ) | 
                                             (red_df['z'] == 85000 ) |
                                             (red_df['z'] == 92500 )  |  
                                             (red_df['z']== 100000 )  ]
            
        print("---Saving pickle file " , name )
        red_df.to_pickle(name)
        return red_df


    
def split_train_test_data(df, plevel = 85000, all_level=True):
    '''
    Split the input dataframe into train and test data.
    If all_level == False, it will use the series at the plevel indicated to split into train and test, 
    and use the last 30% of the data as test.
    If all_level == True, it will use all the other plevels as train set and the given plevel as test set.
    
            Parameters:
                    plev(int):  select the given pressure level. The last 30% of the series is used as a test.
                    all_level(bool): if True, it will consider the time series of the entire level plevel give as

            Return:
                    pandas dataframe/series for train and test data
    '''    
    
    train_labels = ['fg_dep', 'fg_dep_adj','obs_adj','obs', 'an_dep', 'an_dep_adj', 'time'] 
    target_labels =  ['bias']

    
    # selecting one single pressure level
    df_z = df.loc[ df['z'] == plevel ]

    # use only one plevel data and split into training-testing 70/30 %
    if not all_level:
        X_train, X_test, y_train, y_test = train_test_split( df_z[train_labels], df_z[target_labels],  test_size = 0.30, random_state = 42, shuffle=False )

    # use the series of the specified plevel as a test, all the rest as training 
    else:
        X_test = df_z[train_labels]
        y_test = df_z[['bias']]
        
        X_train_all = df.loc[ df['z'] != plevel ]
        X_train = X_train_all[train_labels]
        y_train = X_train_all[['bias']]
        
        
    return X_train, X_test, y_train, y_test 
        
        
def make_plot(X_train_time , y_train, results, title='', zoom=False ):
    """ Creates a plot of the comparison between test data and predictions """
    
    fs = 15
    plt.figure(figsize=(20,10))
    plt.plot( X_train_time['time'].values , y_train['bias'].values, label = 'Observed Bias ')
    plt.plot( results['time'].values , results['test'].values, label = 'Test Bias ')
    plt.plot( results['time'].values , results['prediction_rf'].values, label = 'Random Forest Bias ')
    plt.plot( results['time'].values , results['prediction_dt'].values, label = 'Decision Tree Bias ')

    plt.title(title)
    plt.ylabel("Bias [K]")
    plt.legend(fontsize = fs)
    
    plt.xlim(pd.Timestamp('19490101'), pd.Timestamp('20200101') )
    
    if zoom:
            plt.xlim(pd.Timestamp('20120101'), pd.Timestamp('20200101') )

    plt.show()
    
    
def make_seaborn_scattermatrix(df, title=''):
    plot = sb.pairplot(df, corner=True)
    plot.fig.suptitle(title)