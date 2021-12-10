"""
.. module:: data_preparation
   :synopsis: Module for reading and preparing CUON data for ML applications.
.. moduleauthor:: Federico Ambrogi <federico.ambrogi88@gmail.com>
"""

import h5py 
import logging
import pandas as pd
import numpy as np



# create logger for messages
logger = logging.getLogger('--LOG--')
logger.setLevel(logging.INFO)

class CUON_data():
     """ Main class to read CUON data """
     
     def __init__(self, station_id ="", file= ""):
          """ Iniziale file """
          
          self.file = file
          self.station_id = station_id
          
          
     def read_CUON_data(self, pandas = False, round_time = False, var = 85):
          """ Read a CUON data file
               Parameters:
                    pands (bool): flag to return pandas dataframes if True. Return a dictionary of numpy arrays otherwise
                    round_time(int): time interval witin wich hours are rouded to standard values 00 and 12 .
                                               If False, will extract exact times only.
                                               Range: [0,6] hours 
                    
          """
          # Extract standard plevels
          h5_file  = h5py.File(self.file , 'r')
          
          # Chunking for the correct variable indices
          ind = h5_file["recordindices"][str(var)]
          imin, imax = ind[0], ind[-1]
          

          # Extracting only standar plevels data
          plev = h5_file["observations_table"]["z_coordinate"][imin:imax][:]
          
          logging.info('*** Extracting standard plevels ')          
          ind_plev = np.where( (plev == 1000 )   | (plev == 2000 )   | (plev == 3000 )   | (plev == 5000 )   | (plev == 7000 )   |
                                                 (plev == 10000 ) | (plev == 15000 ) | (plev == 20000 ) | (plev == 25000 ) | (plev == 30000 ) | (plev == 40000 ) |
                                                 (plev == 50000 ) | (plev == 70000 )  | (plev == 85000 ) | (plev == 92500 ) | (plev == 100000 ) ) [0]
          
          # Extract times          
          dt = h5_file["observations_table"]["date_time"][:][imin:imax][ind_plev]
          dt = pd.to_datetime( dt,  unit='s', origin=pd.Timestamp('1900-01-01') )
          
 
          
          #hours = np.where( hours >= 18, 0, hours  )
          #hours = np.where( (hours <=6), 0, hours )
          #hours = np.where( (hours >=6) & (hours <=18) , 12, hours)
          
          days = dt.day   
          
          # Rounding hours and extracting indices 
          if not round_time: 
               logging.warning('*** No time rounding was selected, only exact 00 and 12 hours will be kept!')
               ind_h = np.where( ( hours[:] ==0) | (hours[:]  ==12 ) )[0]
               
          else:
               logging.info('*** Rounding times using a %f time window', round_time)
               
               lim_inf_00, lim_sup_00 = 24-round_time , round_time 
               lim_inf_12, lim_sup_12 = 12-round_time , 12+round_time 
               
               hours = dt.hour
               ind_h = np.where( ( hours[:] <= lim_sup_00) | 
                                             (  (hours[:]  >= lim_inf_12 ) &   (hours[:]  <= lim_sup_12 )   ) |
                                             (  (hours[:]  >= lim_inf_00 ) ) ) [0]
               
               hours =hours[ind_h]
               hours = np.where( hours >= 18, 0, hours   )
               hours = np.where( (hours  <=6), 0, hours  )
               hours = np.where( (hours  >=6) & (hours  <=18) , 12, hours )
               print(0)
               
          # Check finite values for data
          obsval_ind =  np.isfinite( h5_file["observations_table"]["observation_value"][:][imin:imax][ind_plev][ind_h] )
          
          data = {}
          
          for col in ['observation_value', 'z_coordinate','date_time']:
               data[col] = h5_file["observations_table"][col][imin:imax][ind_plev][ind_h][obsval_ind]
         
          for col in ['biascorr@body', 'fg_depar@body','an_depar@body']:
               data[col] = h5_file["era5fb"][col][imin:imax][ind_plev][ind_h][obsval_ind]          
          
          data['months'] = dt[ind_h][obsval_ind].month
          data['years'] = dt[ind_h][obsval_ind].year
          data['days'] = dt[ind_h][obsval_ind].day
          data['hours'] = hours[obsval_ind]
          
          # create a pandas DF if option selected      
          if pandas:
               data = pd.DataFrame( data )
          
          return data
     
     


cuon_vienna = '/mnt/users/scratch/leo/scratch/converted_v7/0-20001-0-11035_CEUAS_merged_v1.nc'
statid_vienna = 'Vienna(0-20001-0-11035)'

CUON = CUON_data(station_id= statid_vienna, file = cuon_vienna )

data = CUON.read_CUON_data(pandas=True, round_time=2)

print(0)



