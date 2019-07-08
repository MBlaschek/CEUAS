""" Analyse and plot the u,v wind components time series from netCDF files.

author: Ambrogi Federico

This files gives retrives basic information and creates the plots
for the u,v wind compontents from the nectCDF files.

The netCDF files were converted from the odb format to net CDF using the 
version of script 'readodbstationfiles.py'.
Databases included: ['1', 3188','1759','1761'] in the root dir '/raid60/scratch/leo/scratch/era5/odbs/'

netCDF file root directory:  '/raid8/srvx1/federico/odb_netCDF/$FSCRATCH/ei6'

"""
import netCDF4
import matplotlib.pylab as plt
import os,sys
import os.path
import matplotlib.gridspec as gridspec








file_u, file_v = "ERA5_1_10393_u.nc" , "ERA5_1_10393_v.nc"

u = netCDF4.Dataset(file_u)
v = netCDF4.Dataset(file_v)


print (v.variables.keys())

datum = [ 1900 + d for d in u.variables['datum'][0,:]/365.25 ] 

'''
interval = 10
for p in range(15):
    u_an = u.variables['fg_dep'][0,p,:]
    u_fg = u.variables['an_dep'][0,p,:]

    plt.plot(datum[::interval] , u_an[::interval] , label = 'an_dep', color = 'blue' )
    plt.plot(datum[::interval] , u_fg[::interval] , label = 'fg_dep', color = 'lime' )

    plt.ylabel('u-wind departues')                                                                                                                                                                                                                                                  
    plt.ylim(-20,20)                                                                                                                                                                                                                                                                 
    plt.xlabel('Date')                                                                                                                                                                                                                                                                 
    plt.legend(loc = 'lower left')                                                                                                                                                                                                                                                    
    plt.savefig('plots/' +  str(p) + '_uwind_departures_series.png',  bbox_inches='tight')                                                                                                                                                                                                                
                
    plt.close()   


    v_an = v.variables['fg_dep'][0,p,:]
    v_fg = v.variables['an_dep'][0,p,:]

    plt.plot(datum[::interval] , v_an[::interval] , label = 'an_dep', color = 'blue' )
    plt.plot(datum[::interval] , v_fg[::interval] , label = 'fg_dep', color = 'lime' )

    plt.ylabel('v-wind departues for p=' + str(p) )
    plt.ylim(-20,20)

    plt.xlabel('Date')
    plt.legend(loc = 'lower left')
    plt.savefig('plots/' + str(p) + '_vwind_departures_series.png',  bbox_inches='tight')

    plt.close()
'''


index = 4809 # index of the year 2005 in datum
bins = 50
for p in [0,5,10]:
    u_an = u.variables['fg_dep'][0,p,:]
    u_fg = u.variables['an_dep'][0,p,:]
    
    an_before, an_after = list(u_an[:4808]) , list(u_an[4809:])
    fg_before, fg_after = list(u_fg[:4808]) , list(u_fg[4809:])

   
    print (an_before)
    print (an_after)

    plt.hist(an_before, histtype = 'step', bins = bins, color = 'blue', label = 'an_dep <2005 p='+str(p) , ls = '-')
    plt.hist(an_after,  histtype = 'step', bins = bins, color ='blue',    label ='an_dep >2005 p='+str(p) , ls = '--')

    plt.hist(fg_before, histtype = 'step', bins = bins, color = 'lime', label = 'fg_dep <2005 p='+str(p) , ls = '-')
    plt.hist(fg_after,  histtype = 'step', bins = bins, color ='lime',    label ='fg_dep >2005 p='+str(p) , ls = '--')



    plt.title('u-wind departures distributions for p=' + str(p))
    plt.xlabel('Departure')
    plt.xlim(-15,15)
    plt.legend(fontsize = 7, loc = 'best')
    plt.savefig('plots/_vwind_departures_distr_p'+str(p)+'.png',  bbox_inches='tight')
    plt.close()



''' 
datum = [ 1900 + d for d in u.variables['datum'][0,:]/365.25 ]





    print(f.variables.keys())
    plt.ylabel(var)
    plt.xlim(20,120)
    plt.xlabel('Year')
    plt.legend(loc = 'lower left')
    plt.savefig(var + '.png',  bbox_inches='tight')
    plt.close()
'''

