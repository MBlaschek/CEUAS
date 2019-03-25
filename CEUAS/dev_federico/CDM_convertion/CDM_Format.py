B""" Data structure and information as defined in the (preliminary) CDM format 

    Author:: Ambrogi Federico , federico.ambrogi@univie.ac.at
    Source:: https://github.com/glamod/common_data_model/

    See definition in  'Common Data Model for in situ observations' document
    See variables in   'observed_variables.dat' """    





Data = { 'observation_table' : { 'lat'   : { 'name' :'latitude'   , 
                                             'def'  :'Latitude of station, -90 to 90 (or other as defined by station_crs'        ,
                                             'units': '' , }      , 
                                 'long' : { 'name' :'longitude'   , 
                                             'def'  :'Longitude of station, -180.0 to 180 (or others as defined by station_crs)' ,
                                             'units': '' , }      ,
                                },
                          
                          
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
                                              'units': 'g kg-1' }
                                 }
}


def Print_Latex():
    out = open('CDM_Table.txt','w')
    top = [r'\begin{table}[!htbp] ',
           r'\footnotesize',
           r'\begin{center}',
           r'\renewcommand{\arraystretch}{1.3}',
           r'\begin{tabular}{ l l l p{3.5in}} ',
           r'\toprule ',
           r'\textbf{Variable} & \textbf{CDM Name} & \textbf{Units} & \textbf{Description}  \\ \toprule \toprule' ]

    bottom = [r'\bottomrule \\bottomrule',
              r'\end{tabular}',
              r'\end{center}',
              r'\caption{Definition of naming convention, description and units for the variables contained in the netCDF files.}',
              r'\label{CDM}', 
              r'\end{table}' ]
 
    for t in top:
        out.write(t + '\n')
       
 
   
    D = Data['observed_variable']
    for k in D.keys():
        #print k , D[k]['name'] , D[k]['def'] , D['units'] , '\n'
        print k 
        n   = D[k]['name'] 
        d   = D[k]['def'] 
        un  = '[' + D[k]['units'] + ']'
        line = k + ' & ' + n + ' & ' + un + ' & ' + d + '\\ ' + '\n'
        out.write(line)
   
    for b in bottom:
        out.write(b + '\n')



a = Print_Latex()












