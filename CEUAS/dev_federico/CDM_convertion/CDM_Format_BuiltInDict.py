""" Data structure and information as defined in the (preliminary) CDM format 

    Author:: Ambrogi Federico , federico.ambrogi@univie.ac.at
    Source:: https://github.com/glamod/common_data_model/

    See definition in  'Common Data Model for in situ observations' document ("cdm_latest.pdf")
    See variables in   'observed_variables.dat' """    



# current variables in netCDF files:
# ['date','time','obstype','codetype','lat','lon','stalt','vertco_type','vertco_reference_1','varno','obsvalue']

# observation_Tables: data taken from  Tab 4.2 of the "cdm_latest.pdf" 

'''
need to clarify what is the variable 
observation_heihgt_above.... (page 20) vs height_of_station_above_sea_level (page 18) ???
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
                                                'def': 'Speed is the magnitude of velocity. Wind is defined as a two-dimensional (horizontal) air velocity vector,  with no vertical component. '
                                                       '(Vertical motion in the atmosphere has the standard name upward air velocity.) The wind speed is the magnitude of the wind velocity. Lot 1 uses ff  - WMO abbrev.',
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



odb_vars_numbering = { 'Temperature'            : 2 , 
                       'Wind speed'             : 110,
                       'Wind u-component speed' : -999 , 
                       'Wind v-component speed' : -999 ,
                       'Specific Humidity'      : 7 ,
                       'Dew Point Temp.'        : 59

                      }


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
        print k 
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











