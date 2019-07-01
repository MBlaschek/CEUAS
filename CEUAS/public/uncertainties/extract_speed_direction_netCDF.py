""" Module to extract the wind speed and direction from the netCDF files with u and v components 

    Author :: Federico Ambrogi, federico.ambrogi@univie.ac.at
    Usage  :: python extract_speed_direction_netCDF.py
    Output :: netCDF files with wind speed and wind direction and analysis/background departures """



import matplotlib.gridspec as gridspec                                                                                          
import netCDF4 
import numpy as np
import datetime
import os,sys
from shutil import copyfile


class windAnalysis:

    def __init__(self, file_u= "" , file_v= ""):
        """ Stores the netCDF files of u and v wind components """
        self.file_u = netCDF4.Dataset(file_u)
        self.file_v = netCDF4.Dataset(file_v)
        
    def read_data(self):
        """ Read the data from the netCDF files for u and v wind components,
            calculate departure values """
            
        self.datum   = self.file_u.variables['datum'][:]

        self.obs_u   = self.file_u.variables['uwind'][:,:,:]
        self.andep_u = self.file_u.variables['an_dep'][:,:,:]
        self.fgdep_u = self.file_u.variables['fg_dep'][:,:,:]
        self.an_u    = self.obs_u - self.andep_u   
        self.fg_u    = self.obs_u - self.fgdep_u 

        self.obs_v   = self.file_v.variables['vwind'][:,:,:]
        self.andep_v = self.file_v.variables['an_dep'][:,:,:]
        self.fgdep_v = self.file_v.variables['fg_dep'][:,:,:]
        self.an_v    = self.obs_v - self.andep_v                                                                                                                                                                                       
        self.fg_v    = self.obs_v - self.fgdep_v

        self.file_u.close()
        self.file_v.close()

    def calc_square(self, u_comp= '', v_comp=''):
        """ Extract the sqrt of the sum of the squares of the given vectors """ 
        return np.sqrt( np.power(u_comp,2) + np.power(v_comp, 2) )

    def extract_direction(self, u_comp="", v_comp="" ):
        """ Extract the direction from the u and v components """
        directions = np.zeros([2,16,len(u_comp[1,1,:])])

        for x in [0,1]:
            for y in range(16):
                for z in range(len(u_comp[1,1,:])): 
                    a = 270 - (180/np.pi)*np.arctan2(v_comp[x,y,z], u_comp[x,y,z])
                    if a > 360:
                        a = a -360
                    elif a == 360:
                        a = 0
                    #print('checking direction! ', u_comp[x,y,z], ' ' , v_comp[x,y,z], ' ' , a )
                    directions[x,y,z] = a

        #for u,v,a_r,r in zip (u_comp[1,10,:], v_comp[1,10,:], angle_rot[1,10,:], angle[1,10,:]):
        #    print("u_comp, v_comp, angle_rot, angle", u, v, a_r, r)

        return directions


    def dep_angle(self,x,y):
        """ Extract the smalles angular difference between two angles.
            Gives minus sign if dep > obs.
            The case of exactl 180 degree diff. is ambiguous 
            Example: obs=20  , an=10 -> dep=10                                                                         
                     obs=10  , an=20 -> dep=-10                                                                                                                                     
                     obs=10  , an=350 -> dep=-20                                                                                                                                                                
                     obs=350 , an=10 -> dep=+20 """
                 
        dep = min( 360 - abs(x-y) , abs(x-y))
        if x<y:
            dep = -dep 
        return dep


    def calc_direction_dep(self):
        """ Extract departures from the observed and analysis/backgrund values """

        self.obs_dir     = self.extract_direction(u_comp= self.obs_u, v_comp= self.obs_v )
        self.dir_an      = self.extract_direction(u_comp= self.an_u , v_comp= self.an_v)
        self.dir_fg      = self.extract_direction(u_comp= self.fg_u , v_comp= self.fg_v)

        dir_andep = np.zeros([2,16,len(self.obs_dir[1,1,:])])
        dir_fgdep = np.zeros([2,16,len(self.obs_dir[1,1,:])])

        for x in [0,1]:
            for y in range(16):
                for z in range(len(self.obs_dir[1,1,:])):
                    obs = self.obs_dir[x,y,z]
                    an  = self.dir_an[x,y,z]
                    fg  = self.dir_fg[x,y,z]
                    obs_an = self.dep_angle(obs,an)
                    obs_fg = self.dep_angle(obs,fg)

                    dir_andep[x,y,z]=obs_an
                    dir_fgdep[x,y,z]=obs_fg

        self.dir_andep   = dir_andep                                                                                                                                                                                                                                    
        self.dir_fgdep   = dir_fgdep
        
    def calc_speed_dep(self):
        """ Extract first guess and analysis departures for the wind speed and direction """
        self.obs_speed   = self.calc_square( u_comp= self.obs_u, v_comp= self.obs_v )
        self.speed_an    = self.calc_square( u_comp= self.an_u, v_comp= self.an_v )
        self.speed_andep = self.obs_speed - self.speed_an 
        self.speed_fg    = self.calc_square( u_comp= self.fg_u, v_comp= self.fg_v )
        self.speed_fgdep = self.obs_speed - self.speed_fg


    def create_netCDF(self, template = "", out_dir = "", file_name= ''):
        """ Create netCDF files for the wind speed and direction, as well as the first guess and analysis departures """

        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)

        temp = netCDF4.Dataset(template) #  Using a template

        file_out_s = out_dir+ '/' +file_name+ '_speed.nc'
        file_out_d = out_dir+ '/' +file_name+ '_direction.nc'

        copyfile(template, file_out_d)
        copyfile(template, file_out_s)

        out_s = netCDF4.Dataset(file_out_s, 'w') # output speed netCDF file
        out_d = netCDF4.Dataset(file_out_d, 'w') # output direction file

        for i in temp.ncattrs(): # FF f.ncattrs=['Conventions', 'title', 'institution', 'history', 'source', 'references']                                                                                                                                               
            if i=='history':
                setattr(out_s, i, datetime.date.today().strftime("%Y/%m/%d"))
                setattr(out_d, i, datetime.date.today().strftime("%Y/%m/%d"))

            elif i=='source':
                setattr(out_s,i, 'Calculated wind speed and departures from the u,v netCDF files for '     + file_name )
                setattr(out_d,i, 'Calculated wind direction and departures from the u,v netCDF files for ' + file_name )


        for dim in list(temp.dimensions.keys() ):
            out_s.createDimension(dim, len(temp.dimensions[dim]) )
            out_d.createDimension(dim, len(temp.dimensions[dim]) )


        for i in list(temp.variables.keys()):
            var = temp.variables[i]

            if i=='datum':
                print(var)
                out_s.createVariable(i, var.dtype, var.dimensions)
                out_s.variables[i][:]= self.datum[:]
                out_d.createVariable(i, var.dtype, var.dimensions)
                out_d.variables[i][:]= self.datum[:]

            if i=='fg_dep':
                out_s.createVariable(i, var.dtype, var.dimensions, fill_value=np.nan)
                #print('sizes', out_s.variables[i].shape , self.speed_fgdep.shape ) 
                out_s.variables[i][:,:,:]= self.speed_fgdep[:,:,:]
                out_d.createVariable(i, var.dtype, var.dimensions, fill_value=np.nan)
                out_d.variables[i][:,:,:]= self.dir_fgdep[:,:,:]

            if i=='an_dep':
                out_s.createVariable(i, var.dtype, var.dimensions, fill_value=np.nan)
                out_s.variables[i][:,:,:]= self.speed_andep[:,:,:]

                out_d.createVariable(i, var.dtype, var.dimensions, fill_value=np.nan)
                out_d.variables[i][:,:,:]= self.dir_andep[:,:,:]

        out_s.close()
        out_d.close()
        print('*** Wrote the files for the wind speed and direction in the %s directory \n', out_dir)





# Running 
base_dir = 'data/'
file_dic = {'temp':'ERA5_1_10393_t.nc' , 'uwind': 'ERA5_1_10393_u.nc'  , 'vwind' : 'ERA5_1_10393_v.nc' }                                                                                                                          

file_u = base_dir + file_dic['uwind']
file_v = base_dir + file_dic['vwind']

wind = windAnalysis(file_u = file_u, file_v = file_v)
wind.read_data()
wind.calc_speed_dep()
wind.calc_direction_dep()

x,y,z = 1 , 12, 500

a = wind.fg_v[x,y,z]
b = wind.obs_v[x,y,z]
c = wind.fgdep_v[x,y,z]

print('shape is', a.shape )
print(a,b,c)

a = wind.fg_u[x,y,z]
b = wind.obs_u[x,y,z]
c = wind.fgdep_u[x,y,z]

print('shape is', a.shape )
print(a,b,c)

print(file_u)
wind.create_netCDF(template = file_u, out_dir = 'data', file_name= 'ERA5_1_10393')



