""" Module to extract the wind speed and direction from the netCDF files with u and v components 

    Author:: Federico Ambrogi, federico.ambrogi@univie.ac.at

    Usage :: python extract_speed_direction_netCDF.py

    Output:: netCDF files with wind speed and wind direction and analysis,background departures """


import matplotlib.gridspec as gridspec                                                                                                                                                                                                                                          
import netCDF4 
import numpy as np
import datetime
import os,sys
from shutil import copyfile



""" Dirs, definitions, select datasets """
base_dir = 'data/'
file_dic = {'temp':'ERA5_1_10393_t.nc' , 'uwind': 'ERA5_1_10393_u.nc'  , 'vwind' : 'ERA5_1_10393_v.nc' , 'dp':'ERA5_1_10393_dp.nc' , 'rh':'ERA5_1_10393_rh.nc'}
#variables = ['temp','uwind','vwind']                                                                                                                                                                                                                                                  
variables = ['temp', 'uwind','vwind','rh','dp']
stations = ['Lindenberg']



class windAnalysis:

    def __init__(self, file_u= "" , file_v= ""):
        """ Stores the netCDF files of u and v wind components """
        self.file_u = netCDF4.Dataset(file_u)
        self.file_v = netCDF4.Dataset(file_v)
        
    def read_data(self):
        
        self.obs_u   = self.file_u.variables['uwind'][:,:,:]
        self.andep_u = self.file_u.variables['an_dep'][:,:,:]
        self.fgdep_u = self.file_u.variables['fg_dep'][:,:,:]
        #print(self.obs_u )
        self.an_u    = self.obs_u - self.andep_u  # extracting the analysis and first guess values 
        self.fg_u    = self.obs_u - self.fgdep_u 

        self.obs_v   = self.file_v.variables['vwind'][:,:,:]
        self.andep_v = self.file_v.variables['an_dep'][:,:,:]
        self.fgdep_v = self.file_v.variables['fg_dep'][:,:,:]
        self.an_v    = self.obs_v - self.andep_v  # extracting the analysis and first guess values                                                                                                                                                                                     
        self.fg_v    = self.obs_v - self.fgdep_v
        self.file_u.close()
        self.file_v.close()

    def calc_square(self, u_comp= '', v_comp=''):
        return np.sqrt( np.power(u_comp,2) + np.power(v_comp, 2) )

    def extract_direction(self, u_comp="", v_comp="" ):
        """ Extract the direction from the u and v components """
        angle_rot = np.zeros([2,16,len(u_comp[1,1,:])])
        angle = 270 - (180/np.pi)*np.arctan2(u_comp, v_comp)

        for x in [0,1]:
            for y in range(16):
                for z in range(len(u_comp[1,1,:])): 
                    a = angle[x,y,z]
                    if a > 360:
                        a = a -360
                    angle_rot[x,y,z] = a

        return angle_rot
        
        
    def calc_speed_dir_dep(self):
        """ Extracting first guess and analysis departures for wind speed and directions """
        self.obs_speed   = self.calc_square( u_comp= self.obs_u, v_comp= self.obs_v )

        self.speed_an    = self.calc_square( u_comp= self.an_u, v_comp= self.an_v )
        self.speed_andep = self.obs_speed - self.speed_an 
        self.speed_fg    = self.calc_square( u_comp= self.fg_u, v_comp= self.fg_v )
        self.speed_fgdep = self.obs_speed - self.speed_fg

        self.obs_dir     = self.extract_direction(u_comp= self.obs_u, v_comp= self.obs_v ) 
        self.dir_an      = self.extract_direction(u_comp= self.an_u , v_comp= self.an_v)
        self.dir_andep   = self.obs_dir - self.dir_an
        self.dir_fg      = self.extract_direction(u_comp= self.fg_u , v_comp= self.fg_v)
        self.dir_fgdep   = self.obs_dir - self.dir_fg

    def create_netCDF(self, template = "", out_dir = "", file_name= ''):

        os.system('rm -r extracted')
        os.system('mkdir extracted')

        temp = netCDF4.Dataset(template)
        file_out_s = out_dir +'/'+ file_name+'_speed.nc'
        file_out_d = out_dir +'/'+ file_name+'_direction.nc'

        copyfile(template, file_out_d)
        copyfile(template, file_out_s)

        out_s = netCDF4.Dataset(file_out_s, 'w') # output speed netCDF file
        out_d = netCDF4.Dataset(file_out_d, 'w') #  output direction file


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


            if i=='fg_dep':
                out_s.createVariable(i, var.dtype, var.dimensions, fill_value=np.nan)
                print('sizes', out_s.variables[i].shape , self.speed_fgdep.shape ) 
                out_s.variables[i][:,:,:]= self.speed_fgdep[:,:,:]

                out_d.createVariable(i, var.dtype, var.dimensions, fill_value=np.nan)
                out_d.variables[i][:]= self.dir_fgdep[:]

            if i=='an_dep':
                out_s.createVariable(i, var.dtype, var.dimensions, fill_value=np.nan)
                out_s.variables[i][:]= self.speed_andep[:]

                out_d.createVariable(i, var.dtype, var.dimensions, fill_value=np.nan)
                out_d.variables[i][:]= self.dir_andep[:]



        out_s.close()
        out_d.close()
        print('*** Wrote the files for the wind speed and direction \n')










file_u = base_dir + file_dic['uwind']
file_v = base_dir + file_dic['vwind']

wind = windAnalysis(file_u = file_u, file_v = file_v)
wind.read_data()
wind.calc_speed_dir_dep()

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
wind.create_netCDF(template = file_u, out_dir = 'extracted', file_name= 'prova')



