""" Script testing the convertion of humidity/dew point 
    Derives the functionylitis from the extract_Humidity_Temperature script

    Author: Ambrogi Federico

    Usage : (python) test_humidity.py 

    Input : data files contained in the station10393_test directory
            for temperature, dew point, specific and relative humidity    

    Output: print on screen the values of variables
            creates plots in the directory PLOTS 

"""

from extract_Humidity_Temperature import netCDF_Files , humidity_dewPoint

import matplotlib
matplotlib.use('Agg')

import matplotlib.pylab as plt
import os,sys 

input_netCDF_t  = 'station10393_test/ERA5_1_10393_t.nc'  #input file with temperature        
input_netCDF_sh = 'station10393_test/ERA5_1_10393_sh.nc'
input_netCDF_rh = 'station10393_test/ERA5_1_10393_rh.nc'
input_netCDF_dp = 'station10393_test/ERA5_1_10393_dp.nc'



# Loading the netCDFs files
netCDFs = netCDF_Files(file_t  = input_netCDF_t , 
                       file_rh = input_netCDF_rh , 
                       file_sh = input_netCDF_sh , 
                       file_dp = input_netCDF_dp )

dates   = netCDFs.load_datum()
data    = netCDFs.load_data()
plevels = netCDFs.define_plevels()

""" Check if the datum arrays are equivalent for all the input variables,
    and creates a dictionary for each date """
check_datum = netCDFs.check_datum()
datas       = netCDFs.find_datums()



""" Initialize the class """
hum_dp = humidity_dewPoint()


""" Checking the results """

ratios = {}
for p in plevels.keys():

  ratios[p] = {}


  press = plevels[p]  # warning: formulas require p in pascal and not hPa, and here is in hPa

  print ('Analysing the pressure level p: ', p , ' corresponding to a pressure of ', press , ' hPa ' )

  for h in [0,1]:
    ratios[p][h] = []

    print('Analysing the hour: ', h )
    temp = netCDFs.data_t [h,p,:]
    dew  = netCDFs.data_dp[h,p,:]
    spec = netCDFs.data_sh[h,p,:]
    rel  = netCDFs.data_rh[h,p,:]
   
    for  t,dp,sh,rh in zip(temp,dew,spec,rel):
         check_t , check_rh , check_sh , check_dp = netCDFs.check_value( t=t , rh=rh , sh=sh , dp=dp )
         flag = bool(check_t and check_rh and check_sh and check_dp) # check == False if the 4 values are defined
         flag_tobechecked = bool( not check_t and not check_rh and not check_sh and not check_dp) # if nan, the value is False
         #print ('flag is', flag_tobechecked , ' for values ',  check_t , check_rh , check_sh , check_dp  )
         if  flag_tobechecked:
             #print('found all values!' ,  t,dp,sh,rh)
             p_pa = press * 100
             calc_rh = hum_dp.specific_to_relative( t=t , sh=sh , pressure= p_pa )
             #print('Relative Humidity. Calc: ', calc_rh , ' Odb: ', rh )         
             ratio = calc_rh / rh        
             ratios[p][h].append(ratio)
             #raw_input('check ')
         #if not check:
         #    print(' t,dp,sh,rh' ,  t,dp,sh,rh , ' which bools are: ', check_t , check_rh , check_sh , check_dp)
         #    raw_input('check missing values')
         #elif check:
         #    print('found all values!' ,  t,dp,sh,rh)
         

os.system('mkdir PLOTS')
def fast_plot(dic = ratios):
    for p in plevels.keys():
        for h in [0,1]:
            print('Making plot for ', p , ' ' , h )
            plt.grid()
            plt.title('Station 10393')
            plt.plot(dic[p][h] , label = 'Pressure: ' + str(p) + ' , hour: ' + str(h) )
            plt.ylabel('Ratio Calculated/Meas. Relative Humidity ')
            plt.savefig('PLOTS/ratios_' + str(p) + '_' + str(h) + '.pdf' , bbox ='tight')
            plt.close()

a = fast_plot(dic = ratios)

