""" Extract humidity and DewPoints from the temperature in the netCDF files 
    
    Author: Ambrogi Federico , federico.ambrogi@univie.ac.at   
    
    """


from humidity import *

import netCDF4


class netCDF_Files:
    """ Class containing the basic information of the netCDF file
        to be processed """
    
    def __init__(self, file_t='' , file_h='' , file_dp=''):
        self.file_t  = netCDF4.Dataset(file_t ) 
        self.file_h  = netCDF4.Dataset(file_h ) 
        self.file_dp = netCDF4.Dataset(file_dp) 
        self.data_t = ''
        
    def load_datum(self):
        """ Loading the observed values as ([2,:] arrays)"""                
        self.datum_t  = [ 1900 + d for d in self.file_t.variables['datum'][0,:]/365.25 ] 
        self.datum_h  = [ 1900 + d for d in self.file_h.variables['datum'][0,:]/365.25 ] 
        self.datum_dp = [ 1900 + d for d in self.file_dp.variables['datum'][0,:]/365.25 ] 
        
    def load_data(self):
        """ Loading the observed values as ([2,:] arrays)"""        
        self.data_t  =  self.file_t.variables['temperatures'][:,12,:]
        #print(self.file_t.variables['temperatures'][:,12,:])
        self.data_h  = self.file_h.variables['h'] [:,12,:]
        self.data_dp = self.file_dp.variables['dp'][:,12,:]
        
       
class humidity_dewPoint:
    """ Module to extract humitiy and dew point from the temperature """
    
    def __init__(self, t = '', h = '', dp = '' , pressure = 1000 ):
        self.t = t
        self.h = h
        self.dp = dp
                
        self.pressure = pressure
        
        self.vap_FOEEWMO = ''
        self.vap_Bolton  = ''
        
        self.specific_h = ''        
        self.relative_h = ''
        
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
    
    def specific_to_relative(self, t = '', press = '', sh = '' ):
        if not sh:
            sh = self.specific_h
        if not t:
            t = self.t      
        if not press:
            p = self.pressure
          
        try:     
            rel_h =sh2rh_ecmwf(sh, t, press)
            self.relative_h = rel_h
            return rel_h
        except:
            print('Cannot calculate the relative humidity! skipping')
            pass
    
    
        
input_netCDF_t  = 'ERA5_1_10393_t.nc' #input file with temperature 
input_netCDF_h  = 'ERA5_1_10393_h.nc'
input_netCDF_dp = 'ERA5_1_10393_dp.nc'



netCDFs = netCDF_Files(file_t=input_netCDF_t  , file_h = input_netCDF_h , file_dp = input_netCDF_dp)

dates = netCDFs.load_datum()
data  = netCDFs.load_data()


""" Loading the data, datum (as lists) """
data_t  , data_h, data_dp =  netCDFs.data_t , netCDFs.data_h , netCDFs.data_dp
datum_t, datum_h, datum_dp = netCDFs.datum_t , netCDFs.datum_h , netCDFs.datum_dp


for temp in data_t[0,:]:
    print('For the temperature ', temp ) 
    
    h_dp = humidity_dewPoint(t = temp, h = '',  dp = '')
    foe = h_dp.vapor_FOEEWMO(temp)
    bolt = h_dp.vapor_Bolton(temp)
    
    sh = h_dp.vapor_to_specificHumidity(t = temp, vap = foe , press = 1000 )
    rh = h_dp.specific_to_relative      (t = temp, sh  = sh,   press = 1000 )
    
    
    print(' the saturation water vapor in Pa is: ', foe, ' the specific humidity is:' , sh , ' the relative humidity is:' , rh)
    
    #input('continue')
    
    
    
    
    
#print(len(data_t), len(data_h), len(data_dp) )
print(' i am printing')

#print('temp is', data_t)
#print('datum is', datum_t)
    
'''  
class humidity_dewPoint:
    
    def __init__(self):
       '''
    

