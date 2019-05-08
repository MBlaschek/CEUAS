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

t_file = 'data/ERA5_1_10393_t.nc'


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


def running_mean(x,N):
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)


'''
x = np.random.random(100000)

res_1 = running_mean(x, 1000)
print('res1', res_1, len(res_1))
res_2 = running_mean(x, 2000)
print('res2', res_2,  len(res_2))
res_3 = running_mean(x, 100)
print('res3', res_3,  len(res_3))
'''

print('*** Loading the file> ', t_file )
data = netCDF().read_data(file = t_file)

andep = data['an_dep']['data'] 
fgdep = data['fg_dep']['data'] 

print('*** Check the data shape for fg_dep and an_dep: ', andep.shape, fgdep.shape )



#test_x , test_y = [1,2,3] , [-1,9,0]
#matrix = calc_cov(test_x , test_y)

def extract_entry(matrices = '', pressure_i = '' , pressure_j = ''):
    """ Extract a list of entries of the covariance matrices , e.g. all the "ith,jth" entries of all the matrices """
    lista = [ m[pressure_i, pressure_j] for m in matrices ] 
    return lista
        




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



matrices = {'0': '' , '12':'' }

""" Extract the list of matrices, for all the days of observation """
for hour in [0]:
    lista = []
    h = str(hour)                                                                                                                                                                                                                                                                   
    for date in range(len(datums)):
        an = andep[hour,:,date]
        fg = fgdep[hour,:,date]
        corrMatrix = calc_cov(an,fg) 
        lista.append(corrMatrix)
    matrices[h] = lista
#print ('*** The list of matrices are:', matrices )

entry_1_1 = [ x for x in extract_entry(matrices = matrices['0'], pressure_i = 14 , pressure_j = 14) if not np.isnan(x) ]


#a = [ x for x in entry_1_1 if not np.isnan(x)  ]
#print('The entries are: ', entry_1_1 )
#a =  running_mean(a , 30 )
#print(a)
#input('')

""" Extracting the running means for various time intervals """
rm_1m = [x for x in running_mean(entry_1_1 , 30 ) if not np.isnan(x) ] 
rm_2m = [x for x in running_mean(entry_1_1 , 60 ) if not np.isnan(x) ] 
rm_3m = [x for x in running_mean(entry_1_1 , 90 ) if not np.isnan(x) ]
rm_6m = [x for x in running_mean(entry_1_1 , 180) if not np.isnan(x) ]
rm_1y = [x for x in running_mean(entry_1_1 , 365) if not np.isnan(x) ]





""" Plotting the histograms """

X = [ rm_1m, rm_2m, rm_3m, rm_6m, rm_1y ]
C = ['slateblue', 'cyan', 'lime', 'orange', 'gold'] 
L = ['Desroziers(1m)', 'Desroziers(2m)', 'Desroziers(3m)', 'Desroziers(6m)', 'Desroziers(1y)']


print(X)
plt.title('Estimated observation errors for the temperature')
Bins = 50
FONTSIZE = 15
n, bins, patches = plt.hist(X, Bins, histtype='stepfilled' ,  stacked = False, color = C , label = L , alpha = 0.7 , normed = True)
plt.text(1, 0.28,"pressure(an,fg)=(0,0)", fontsize= FONTSIZE)
plt.grid(linestyle= ':', color = 'lightgray', lw = 1.2 )
plt.legend(loc = 'upper right', fontsize = FONTSIZE - 3)
plt.ylabel('Numbers / '+str(Bins), fontsize = FONTSIZE)
plt.ylim(0, 0.30)
plt.xlabel(r'Errors [m/s]', fontsize=15)
plt.savefig('plots/Temperature.pdf',  bbox_inches='tight')
plt.close()






'''
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
        


