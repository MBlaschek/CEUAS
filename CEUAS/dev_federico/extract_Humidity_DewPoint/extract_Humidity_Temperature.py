""" Extract humidity and DewPoints from the temperature in the netCDF files 
    
    Author: Ambrogi Federico , federico.ambrogi@univie.ac.at   
    
    for CDM compliant naming, see the COnventions.pdf file on GitHub
    """


from humidity import *
from shutil import copyfile
import netCDF4
import numpy as np

'''
_t.nc variables: ['lat', 'lon', 'alt', 'press', 'datum', 'hours', 'source', 'odbstatid', 'temperatures', 'fg_dep', 'bias', 'an_dep', 's_type']
_h.nc variables: ['lat', 'lon', 'alt', 'press', 'datum', 'hours', 'source', 'odbstatid', 'h', 'fg_dep', 'bias', 'an_dep', 's_type']
'''
class netCDF_Files:
    """ Class containing the basic information of the netCDF file
        to be processed """
    
    def __init__(self, file_t='' , file_rh='' , file_sh='', file_dp=''):
        
        self.file_t  = netCDF4.Dataset(file_t ) 
        self.file_rh = netCDF4.Dataset(file_rh)
        self.file_sh = netCDF4.Dataset(file_sh )         
        self.file_dp = netCDF4.Dataset(file_dp) 
        
        self.old_file = file_t
        self.new_file = file_t.replace('_t.nc','_calc_sh_rh_dp.nc')
        
    def define_plevels(self):
        """ Create a dict with the pressure level and correpsondind index in the data arrays """
        plevs=np.asarray([10,20,30,50,70,100,150,200,250,300,400,500,700,850,925,1000])
        numbs=np.asarray(range(16))
        dict_plevels = {}
        for p,n in zip(plevs,numbs):
            dict_plevels[n]=p        
        self.plevels = dict_plevels
        return dict_plevels
        
    def load_datum(self):
        """ Loading the observed values as ([2,:] arrays)"""                
        self.datum_t  = [ 1900 + d for d in self.file_t. variables['datum'][0,:]] 
        self.datum_rh = [ 1900 + d for d in self.file_rh.variables['datum'][0,:]] 
        self.datum_sh = [ 1900 + d for d in self.file_sh.variables['datum'][0,:]]         
        self.datum_dp = [ 1900 + d for d in self.file_dp.variables['datum'][0,:]] 
        
    def load_data(self):
        """ Loading the observed values as ([2,:] arrays)"""        
        self.data_t  = self.file_t.variables  ['temperatures']#[:,:,:] (2,16,number)=(2hours,16pressures,observationDays)
        self.data_rh = self.file_rh.variables ['rh'][:,:,:]
        self.data_sh = self.file_sh.variables ['sh'][:,:,:]        
        self.data_dp = self.file_dp.variables ['dp'][:,:,:]
       
    def create_new(self):
        """ Copy the temperature file into the new one containing spec.hum, rel.hum and dew point 
            the new_var dic is contains the names of the new variables as keys and the variable range as values 
        """        
        copyfile(self.old_file , self.new_file)
        out_file = netCDF4.Dataset(self.new_file , "w")
        self.file_sh_rh_dp = out_file 
        
        # using the tmep field as a reference variable
        var = self.file_t.variables['temperatures']
        
        
        new_var = { 'specific_humidity'    : [-100 , 100] , 
                    'relative_humidity'    : [-100 , 100] , 
                    'dew_point_temperature': [0    , 300] ,
                    'specific_humidity'    : [0    , 300]    } # CDM compliant
        
        self.new_vars = new_var
        
        for k,v in new_var.items():
            out_file.createDimension(i,len(self.file_t.dimensions[var]))
            out_file.createVariable(k,  var.dtype , var.dimensions, fill_value=numpy.nan)
            #fo.variables[metadata[t]['varname']][:]=statdata[1][metadata[t]['varname']][:]  # to do later?          
            setattr(out_file.variables[k], j , v) # setting the variable range
        
        self.file_sh_rh_dp = out_file 
       
    def fill_missing(self, var = '' , values = [] ):
        """ Filling the new files with the values. 
            var is the name of the variable 
            values is a 2d numpy array for the two observation hours """
        if var not in self.new_vars.keys():
            raise ValueError ('Variable not recognized!')                     
        self.out_file.variables[var][:]= values[:]         
            
    def find_datums(self):
        ''' Loops over the datums of the temperature, humidities and dew points available 
            and creates a dictionary for each separate entry 
            (in principle different variables might have different datums)
            Returns:
                    dictionary with days as keys, empty values
            '''
        dic = {}
        for datum in [self.datum_dp , self.datum_rh , self.datum_sh , self.datum_t]:
            for day in datum:
                if day not in dic.keys():
                    dic[day] = '' # initialize an empty dictionary
        self.datum_dics = dic
        return dic
        
    def check_datum(self):
        ''' Check if the datum sets are the same for all the variables. 
            In this case, there is no need for further checks.
            Returns:
                    True or False'''
        if self.datum_dp.sort() == self.datum_rh.sort() and self.datum_rh.sort() == self.datum_t.sort() and self.datum_dp.sort() == self.datum_t.sort():
            print('*** The datum lists of all the variables are identical')
            return True
        else:
            print('*** The datum lists of all the variables are different: proceeding with one by one checks')            
            return False
    
    def check_value(self, t='', rh='', sh='', dp=''):
        ''' Check if the values are correct, and not numpys "nan" 
            Returns four boolean values '''
        check_t  = np.isnan(t) 
        check_rh = np.isnan(rh)        
        check_sh = np.isnan(sh)        
        check_dp = np.isnan(dp)        
        return check_t , check_rh , check_sh , check_dp
        
    def check_values(self, fast= True, p=''):
        ''' Checks if the values of some variables are missing 
            If fast=True  , the datum lists have passed the check, and they are identical
            if fast=False , the datum are different and require further analysis 
            '''
        if fast:
            for i in [0,1]: # 2 hours measurements (00:00 . 12:00)
                for t,dp,sh,rh in zip( self.data_t[i,p,:] , self.data_dp[i,p,:] , self.data_sh[i,p,:] , self.data_rh[i,p,:] ):
                    print(self.check_singlevalue (t=t, dp=dp, rh=rh, sh=sh ) )
        
        
    def clean_close(self):
        self.file_sh_rh_dp.close()
        print('Finished with processing the file', self.new_file )
        
       
       
       
       
class humidity_dewPoint:
    """ Module to extract humitiy and dew point from the temperature """
    
    def __init__(self):
                
        self.vap_FOEEWMO = ''
        self.vap_Bolton  = ''
        
        self.specific_h = ''        
        self.relative_h = ''
        
        
    def specific_to_relative(self, t='', sh='', pressure=''):
        try:
            v = sh2rh_ecmwf(sh, t, pressure)
            self.sh2rh_ecmwf = v
            return v
        except:
            print('Cannot convert specific humidity to relative! skipping')
            pass            
        
    def vapor_FOEEWMO(self, t):
        try:
            v = FOEEWMO(t)
            self.vap_FOEEWMO = v
            return v
        except:
            print('Cannot calculate the saturation vapor! skipping')
            pass
        
    def vapor_Bolton(self, t):
        try:
            v = Bolton(t)
            self.vap_Bolton = v
            return v
        except:
            print('Cannot calculate the saturation vapor! skipping')
            pass
        
        
    def vapor_to_specificHumidity(self, t = '', vap = '' , press= ''):
        
        if not vap:
            vap = self.vap_FOEEWMO
        if not t:
            t = self.t      
        if not press:
            p = self.pressure
        
        try: 
            sh = vp2sh(vap, press)
            self.specific_h = sh
            return sh
        except:
            print('Cannot calculate the specific humidity from vapor and pressure! skipping')
            pass            
    
    def relative_to_specific(self, t = '', pressure = '', sh = '' ):       
        try:     
            rel_h =sh2rh_ecmwf(sh, t, press)
            self.relative_h = rel_h
            return rel_h
        except:
            print('Cannot calculate the relative humidity! skipping')
            pass 
    
    
    
class fill_new:
    """ Fills the new files with the missing information """    
    
    def __init__(self):
        pass



'''
input_netCDF_t  = 'ERA5_1_10393_t.nc' #input file with temperature 
input_netCDF_sh = 'ERA5_1_10393_sh.nc'
input_netCDF_rh = 'ERA5_1_10393_rh.nc'
input_netCDF_dp = 'ERA5_1_10393_dp.nc'
'''

# ##################################################################################################################


    
    

'''
netCDFs = netCDF_Files(file_t=input_netCDF_t  , file_rh = input_netCDF_rh , file_sh = input_netCDF_sh , file_dp = input_netCDF_dp)

dates = netCDFs.load_datum()
data  = netCDFs.load_data()
plevels = netCDFs.define_plevels()



""" Loading the data, datum (as lists) """
data_t  , data_rh  , data_sh  , data_dp = netCDFs.data_t ,netCDFs.data_rh, netCDFs.data_sh, netCDFs.data_dp
datum_t , datum_rh , datum_sh, datum_dp = netCDFs.datum_t , netCDFs.datum_rh , netCDFs.datum_sh , netCDFs.datum_dp




#data_t shape (2, 16, 9916)
datas = netCDFs.find_datums()

check_datum = netCDFs.check_datum() # true if datums are identical, false otherwise

netCDFs.check_values(fast = check_datum)


#for t,s,r,d in zip(datum_t , datum_rh , datum_sh, datum_dp):
    #print('checking the datums ', t,s,r,d )

a = netCDFs.check_datum()

sh  = np.empty([2,1] , dtype = float) 
'''


'''
#for pressure in 
for temp in data_t[:,:]:
    print('For the temperature ', temp ) 
    
    h_dp = humidity_dewPoint(t = temp, h = '',  dp = '')
    foe = h_dp.vapor_FOEEWMO(temp)
    bolt = h_dp.vapor_Bolton(temp)
    
    sh = h_dp.vapor_to_specificHumidity (t = temp, vap = foe , press = 1000 )
    rh = h_dp.specific_to_relative      (t = temp, sh  = sh,   press = 1000 )
    
    
    print(' the saturation water vapor in Pa is: ', foe, ' the specific humidity is:' , sh , ' the relative humidity is:' , rh)
    
    #input('continue')
'''
    


'''
todo: 
- initialize file i.e. copy from temperature
- loop over temp
- extract sh, rh, dp ad 2d numpy arrays for the two observation hours
- fill the copied file with the new variables
- close
'''
    
    
#print(len(data_t), len(data_h), len(data_dp) )
print(' i am printing')

#print('temp is', data_t)
#print('datum is', datum_t)
    
'''  
class humidity_dewPoint:
    
    def __init__(self):
       '''
    

