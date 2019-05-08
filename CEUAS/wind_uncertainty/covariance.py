""" Module for extracting the covariance matrix """


import os,sys
import netCDF4
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
#import matplotlib.gridspec as gridspec



# see definition in plotter script 
class netCDF:
     def __init__(self, file=""):
          self.file = file

     def read_data(self, file = '', var = ['fg_dep', 'an_dep']):
         
         data_loaded = {}  # Dictioary containing the extracted information                                                                                                                                                                                                              
         for v in var:
                  data_loaded[v] = {}
                  data_loaded[v]['datum'] = []
                  data_loaded[v]['data' ] = [] # contains either fg_Dep or an_Dep                                                                                                                                                                                                    

                  if os.path.isfile(file):
                      f = netCDF4.Dataset(file)
                      data_loaded[v]['datum'] = [ 1900 + d for d in f.variables['datum'][0,:]/365.25 ] # converting the datum in years                                                                                           
                      data_loaded[v]['data']  = f.variables[v][:,:,:]
                  else: raise ValueError('netCDF files:', path, ' not found!!!')
         return data_loaded 


uwind_file = 'data/ERA5_1_10393_u.nc'



def running_mean(x, N):
    cumsum = numpy.cumsum(numpy.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / float(N)



def calc_cov(array_x, array_y):
    cov = np.empty([len(array_x),len(array_x)], dtype = float) # initialize an empty 16*16 matrix (16 pressure levels)
    for x in range(len(array_x)):
        for y in range(len(array_y)):
            entry = array_x[x] * array_y[y]
            #print (array_x[x] , array_y[y] , entry , x, y)
            cov[x,y] = entry
    #print('*** Covariance matrix:', cov )
    return cov


print('*** Loading the file> ', uwind_file )
data = netCDF().read_data(file = uwind_file)

andep = data['an_dep']['data'] 
fgdep = data['fg_dep']['data'] 

print('*** Check the data shape for fg_dep and an_dep: ', andep.shape, fgdep.shape )



#test_x , test_y = [1,2,3] , [-1,9,0]
#matrix = calc_cov(test_x , test_y)


def cov_plot(matrix, station="", hour = "", date=""):
    """ Basic plot for the correlation matrix """

    fig,ax = plt.subplots()
    FONT = 15

    plt.title("Station: " + station + ', Hour: ' + hour + ', Date: ' + date , y=1.03)
    
    num = len(matrix[0,:])
    Num = range(num)
    vmin, vmax = -10, 10
    color_map= plt.imshow(matrix, interpolation= 'nearest', cmap = 'RdYlBu', vmin = vmin, vmax = vmax ) # nearest serves for discreete grid
    #color_map.set_cmap('Blues_r')
    #color_map.set_cmap('seismic')
    plt.ylim(-0.5, 15.5)
    plt.xlim(-0.5, 15.5)
    plt.xticks(Num, Num)
    plt.xlabel('Pressure level', fontsize = FONT)
    plt.yticks(Num, Num)
    plt.xlabel('Pressure level', fontsize = FONT)
    bar = plt.colorbar()
    bar.ax.set_ylabel("Covariance", fontsize = FONT)

    #  Creating text values
    for i in Num:
        for j in Num:
            value = '{0:.2f}'.format(matrix[i,j])
            #print(value)
            #input('next')
            text = ax.text( j,i, value , ha = 'center' , va = 'center', color = 'black', fontsize = 5)

    name = 'Cov_' + station + '_hour_' + hour.replace(':','') + '_date_' + date 
    plt.savefig('plots/' + name + '.pdf', bbox_inches='tight')
    plt.close()

os.system('mkdir plots')




datums = data['fg_dep']['datum']

station = 'Lindenberg'

""" Plotting some covariances """
for hour in [0,1]:
    hour_name = str(hour).replace('0','00:00').replace('1','12:00')
    for date in [0,100,1000,5000,7000]:

        date_name = str(datums[date])

        print("*** I am extracting the covariances for the date ", date_name , " and hour ", hour_name )

        an = andep[hour,:,date]
        fg = fgdep[hour,:,date]

        corr_matrix = calc_cov(an,fg)  # Extracting the matrix

        cov_plot(corr_matrix , station = station, hour = hour_name , date = date_name)  # Plotting
 
'''
# Example
date = str (data['fg_dep']['datum'][0]) 
hour = 0
Hour = '00'
an = andep[hour,:,0]
fg = fgdep[hour,:,0]
corr_matrix = calc_cov(an,fg)
print('Date is: ', date )
cov_plot(corr_matrix , station = 'Lindenberg', hour = Hour , date = date)
print('andep, fgdep', an, fg)
'''
        


