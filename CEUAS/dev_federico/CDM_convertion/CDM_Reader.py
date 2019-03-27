""" Module for parsing the CDM compliant data and information

    Author      :: Ambrogi Federico , federico.ambrogi@univie.ac.at
    Data Source :: https://github.com/glamod/common_data_model/
    Directory   :: /raid60/scratch/copernicus/common_data_model

    See definitions in  'Common Data Model for in situ observations' document ("cdm_latest.pdf")
    See observation variables in   'observed_variables.dat'
    
    vars in netCDF files: ['date','time','obstype','codetype','lat','lon','stalt','vertco_type','vertco_reference_1','varno','obsvalue']
"""    

import os,sys
import path


""" path to the directory where the github repository was cloned 
    path to the cvs file containing the information of the observations 
    path to the cvs file containing the information of the variables """

git_root_dir = '/raid60/scratch/copernicus/common_data_model'

Observation_file = os.path.join(git_root_dir,'table_definitions','observations_table.csv')

Variables_file = os.path.join(git_root_dir,'tables','observed_variable.dat')





""" Dictionary of variables to be included in the netCDF files
    key   = generic name (appearing e.g. in the readodbstations.py code),
    value = CDM name as found in the the table 
    Note that this MUST be checked and defined by hand """

observation_table = { 'lat'                 : 'latitude'      ,
                      'long'                : 'longitude'     , 
                      'vertco_type'         : 'z_coordinate'  , 
            #'vertco_reference_1'  : 'xx' , 
            #'stalt'               : 'xx' , 
            #'date'                : 'xx' , 
            #'time'                : 'xx' , 
            #'obstype'             : 'xx' , 
            #'codetype'            : 'xx' ,
               }


observed_variable = {'wind'     : 'wind speed', 
                     'pressure' : 'air pressure', 
                     'dewpoint' : 'dew point temperature' , 
                     'humidity' : 'specific humidity', 
                     #'varno'    : 'observed_variable',
                     #'obsvalue' : 'observation_value' }
                     }


def extractObsInfo(lines, var):
    """Loops over the lines of the observation cvs file, checks if contain the wanted variable,
       returns the information of the units (or type of observation) and definition of variable 
       
       Args: 
            lines : lines of the input csv file containing the observation definition 
            var   : variable to be searched for (values of the observation_table dictionary) 
       Returns     
            units (or format of the variable)
            definition of the variable
    """
    
    for l in lines:    
        a = l.split('\t')
        print(a )
        unit, definition = '', ''
        if var == a[0]: 
            unit, definition = a[2] , a[3]
            print('Check ', unit, definition)            
            return unit, definition
        else:
            continue 

def extractVarsInfo(lines, var):
    """ Same as  extractObsInfo but for the observed variables      
       Args: 
            lines : lines of the input csv file containing the variables definition 
            var   : variable to be searched for (values of the observed_variable dictionary) 
       Returns     
            units (or format of the variable)
            definition of the variable
    """
    
    for l in lines:    
        a = l.split('\t')
        print(var, a[4] )
        unit, definition, number = '', '' , ''
        
        if var == a[4]: 
            unit, definition, number = a[5] , a[6] ,  a[0]
            print('Check ', unit, definition)            
            return unit, definition , number
        else:
            continue
        
        



'''
""" Reading the observation file and creating the dictionary """
Dic_Obs = {}
lines = open(Observation_file).readlines()
for O in observation_table.keys():
    Dic_Obs[O] = {}
    unit, definition = extractObsInfo(lines, observation_table[O])
    Dic_Obs[O] ['units']      = unit
    Dic_Obs[O] ['definition'] = definition
print('Check the Obs Dictionary', Dic_Obs)

'''

""" Reading the observation file and creating the dictionary """
Dic_Var = {}
lines = open(Variables_file).readlines()
for O in observed_variable.keys():
    Dic_Var[O] = {}
    unit, definition, number = extractVarsInfo(lines, observed_variable[O])
    Dic_Var[O] ['units']      = unit
    Dic_Var[O] ['definition'] = definition
    Dic_Var[O] ['number'] = definition
    
print('Check the Obs Dictionary', Dic_Var)





'''
Data = { 'observation_table' : { 'lat'   : { 'name' :'latitude'   , 
                                             'def'  :'Latitude of station, -90 to 90 (or other as defined by station_crs'        ,
                                             'units': '' , }      , 
                                 'long' : { 'name' :'longitude'   , 
                                             'def'  :'Longitude of station, -180.0 to 180 (or others as defined by station_crs)' ,
                                             'units': '' , }      ,
                                 
                                 'vertco_type': { 'name' : 'z_coordinate_type'   , 
                                                  'def'  : 'Type of z coordinate' ,
                                                  'units': '' , }      ,
                                 
                                 'vertco_reference_1': { 'name' :'???'   ,
                                                         'def'  :'???' ,
                                                         'units': '' , }      ,
                                 
                                 'stalt': {'name'  : 'observation_height_above_station_surface' ,
                                           'def'   : 'Height of sensor above local ground or sea surface. Positive values for above surface (e.g. sondes), negative for below (e.g. xbt). For visualobservations, height of the visual observing platform.' ,
                                           'units' : '' } ,

                                 'date': { 'name'  : '' ,
                                           'def'   : '' ,
                                           'units' : '' , } ,

                                 'time': { 'name'  : '' ,
                                           'def'   : '' ,
                                           'units' : '' , } ,

                                 'obstype': { 'name'  : '' ,
                                              'def'   : '' ,
                                              'units' : '' , } ,

                                 'codetype': { 'name'  : '' ,
                                               'def'   : '' ,
                                               'units' : '' , } ,

                                } , 
                          
         'observed_variable' : { 'wind'     : { 'name': 'wind'     ,
                                                'def': 'Speed is the magnitude of velocity. Wind is defined as a two-dimensional (horizontal) air velocity vector,  with no vertical component. (Vertical motion in the atmosphere has the standard name upward air velocity.) The wind speed is the magnitude of the wind velocity. Lot 1 uses ff  - WMO abbrev.',
                                                'units': 'm s-1'}      ,
                                 'pressure' : { 'name': 'pressure' ,
                                                'def': 'pressure of air column at specified height',
                                                'units':'Pa'}      ,

                                 'dewpoint': {'name': 'dew point temperature' ,
                                              'def' : 'Dew point temperature is the temperature at which a parcel of air reaches saturation upon being cooled at constant pressure and specific humidity.',
                                              'units': 'K', }                 ,
                                 
                                 'humidity': {'name': 'specific humidity'     ,
                                              'def' : 'specific means per unit mass. Specific humidity is the mass fraction of water vapor in (moist) air.',
                                              'units': 'g kg-1' } ,

                                 'varno' :  { 'name': 'observed_variable'     ,
                                              'def' : 'The variable being observed / measured.',
                                              'units': 'int' } ,

                                 'obsvalue' :  { 'name': 'observation_value'     ,
                                                 'def' : 'The observed value.',
                                                 'units': 'numeric' } ,



                                 }
}
'''


odb_vars_numbering = { 'Temperature'            : 2 , 
                       'Wind speed'             : 110,
                       'Wind u-component speed' : -999 , 
                       'Wind v-component speed' : -999 ,
                       'Specific Humidity'      : 7 ,
                       'Dew Point Temp.'        : 59
                      }




'''
top =     [r'\begin{table}[!htbp] ',
           r'\footnotesize',
           r'\begin{center}',
           r'\renewcommand{\arraystretch}{1.3}',
           r'\begin{tabular}{  l p{1.5in} l p{3.0in} } ',
           r'\toprule \toprule', '\n']
 

bottom = [r'\bottomrule \bottomrule',
              r'\end{tabular}',
              r'\end{center}',
              r'\caption{Definition of naming convention, description and units for the variables contained in the netCDF files.}',
              r'\label{CDM}',
              r'\end{table}' ]





def printLatexTab(top = '' , bottom = '' , data = '' , outname = '' ):
    """ prints a latex style table for text editing 
        out is the name of the output tex file
        must be either ObsTab or VarTab, according to the data chosen"""
    out = open('Tex/' + outname + '.tex','w')
 
    for t in top:
        out.write(t + '\n')

    
    tabType = outname.replace('ObsTab','Observation Table').replace('VarTab','Variables Table')    
    
    out.write(r'\multicolumn{4}{c}{ ' + tabType + r' } \toprule \toprule \\' + '\n' )
    out.write(r'\textbf{Variable} & \textbf{CDM Name} & \textbf{Units} & \textbf{Description}  \\ \toprule' + '\n')
    for k in data.keys():
        #print k , D[k]['name'] , D[k]['def'] , D['units'] , '\n'
        n   = data[k]['name'].replace('_','\_')
        d   = data[k]['def'] .replace('_above',' above').replace('_','\_')
        un  = r'$[$' + data[k]['units'] + r'$]$'


        print(k , n , d, un )
        line = k.replace('_','\_') + ' & ' + n + ' & ' + un + ' & ' + d + r'\\ ' + '\n'
        out.write(line)
   
    for b in bottom:
        out.write(b + '\n')

    out.close()


obs = Data['observation_table']
var = Data['observed_variable']

a = printLatexTab(top = top, bottom = bottom , data = obs , outname = 'ObsTab')
b = printLatexTab(top = top, bottom = bottom , data = var , outname = 'VarTab')


'''








