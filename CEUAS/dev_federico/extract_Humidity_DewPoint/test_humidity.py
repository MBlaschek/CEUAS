""" Script testing the convertion of humidity/dew point 
    Derives the functionylities from  the "extract_Humidity_Temperature.py" script

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


""" Defining the input files for temperature, specific humidity, relative humidity, dew point """
input_netCDF_t  = 'station10393_test/ERA5_1_10393_t.nc'  #input file with temperature        
input_netCDF_sh = 'station10393_test/ERA5_1_10393_sh.nc'
input_netCDF_rh = 'station10393_test/ERA5_1_10393_rh.nc'
input_netCDF_dp = 'station10393_test/ERA5_1_10393_dp.nc'

""" Loading the netCDF files """
netCDFs = netCDF_Files(file_t  = input_netCDF_t , 
                       file_rh = input_netCDF_rh , 
                       file_sh = input_netCDF_sh , 
                       file_dp = input_netCDF_dp )

""" Loading the variables """

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
ratios_f = {}
ratios_h = {}

for p in plevels.keys():
#for p in [15]:
  ratios_f[p] = {}
  ratios_h[p] = {}

  press = plevels[p]  # warning: formulas require p in pascal and not hPa, and here is in hPa

  print ('Analysing the pressure level p: ', p , ' corresponding to a pressure of ', press , ' hPa ' )

  for h in [0,1]:
    ratios_f[p][h] = []
    ratios_h[p][h] = []

    print('Analysing the hour: ', h )
    temp = netCDFs.data_t [h,p,:]
    dew  = netCDFs.data_dp[h,p,:]
    spec = netCDFs.data_sh[h,p,:]
    rel  = netCDFs.data_rh[h,p,:]
   

    for  t,dp,sh,rh in zip(temp,dew,spec,rel):
          check_t , check_rh , check_sh , check_dp = netCDFs.check_value( t=t , rh=rh , sh=sh , dp=dp )
          print('check', check_t , check_rh , check_sh , check_dp )
          input('')
          flag = bool(check_t and check_rh and check_sh and check_dp) # check == False if the 4 values are defined
          flag_tobechecked = bool (not check_t and not check_rh and not check_sh and not check_dp) # if nan, the value is False
          
          if  flag_tobechecked:
            p_pa = press * 100
                        
            ''' using the sat. vapor formula at T=T(dew point) '''   
            
            # ##################################### dsljkmfnsdfjmn to do! should be 1 not zero!!!!
            calc_f = hum_dp.vapor_FOEEWMO(dp)/hum_dp.vapor_FOEEWMO(t)  # using the sat. vapor at dew point      
            ratiof = calc_f / rh        
            ratios_f[p][h].append(ratiof)
                        
            ''' using the specific to relative humidity conv. formula '''
            calc_h = hum_dp.specific_to_relative( t=t, sh=sh, pressure= p_pa )
            ratioh = calc_h / rh
            ratios_h[p][h].append(ratioh) # using the convertion formula      

            r2 = calc_h / calc_f
            
            #print ('checking the values: ', calc_h , calc_f , r2 , ' for t and dp and sh ', t , dp, sh )
            #print ('checking the ratios: ', ratiof , ratioh , rh , ' for t and dp and sh ', t , dp, sh )


""" Making polots to compare the values of the calculated relative humidity 
    (using directly the vapor_FOEEWMO formula with the dew point, or the  )
    with the value stored in the odb files """
os.system('mkdir PLOTS/')
def fast_plot(dic = '' , what = '' , c = ''):
    os.system('mkdir PLOTS/'+what)
    
    for p in plevels.keys():
    #for p in [15]:
  
        for h in [0,1]:
            print('Making plot for ', p , ' ' , h )
            plt.grid(linestyle = ':', color = 'lightgray')
            plt.title('Station 10393 - Pressure: ' + str(p) + ' , hour: ' + str(h) )
            plt.plot(dic[p][h] , label = what , color = c )
            plt.ylabel('Ratio Calculated/Meas. Relative Humidity')
            plt.xlabel('Observation Day ')
            #plt.ylim(0.95 , 1.05)
            plt.legend(fontsize = 12)
            plt.savefig('PLOTS/'+ what + '/ratios_'+ what + '_' + str(p) + '_' + str(h) + '.pdf' , bbox ='tight')
            plt.close()

for c,d,l in zip( ['blue','lime'], [ratios_f , ratios_h ], ['FOEEWMO','sh2rl']):
  a = fast_plot( dic = d , what = l , c = c )

