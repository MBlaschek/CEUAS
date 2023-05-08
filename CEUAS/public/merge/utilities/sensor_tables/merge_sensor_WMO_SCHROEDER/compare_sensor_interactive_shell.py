import pandas as pd 
import numpy as np
import os,sys


""" Reading the file with all the pairs from last Schoreder sensor and irst WMO sensor 
    ['wmo','wmo_comment', 'sch' , 'sch_comment' , 'min_date' , 'max_date' , 'station'] """



df = pd.read_csv('sensor_from_sch_to_wmo.csv', sep = '\t')



if os.path.isfile('dictionary_sensor_mapping.npy'):
    all_sch = np.load('dictionary_sensor_mapping.npy', allow_pickle=True).item()

else:
    all_sch = {}
    all_sch['stations'] = []



"""
Loop over all the unique WMO sensor, and extract a reduced dataframe
Check for all possible Sch sensor
if they have been analyzed already or not
If yes, skip
if not, an interactive shell will ask to select if compatible or not  
"""


unique_wmo = list(np.unique( df.wmo.values ))

comp  = { 'wmo':[] , 'sch':[] , "wmo_c":[], "sch_c" : [] , "stations":[] }
wrong = { 'wmo':[] , 'sch':[] , "wmo_c":[], "sch_c" : [] , "stations":[]}
un    = { 'wmo':[] , 'sch':[] , "wmo_c":[], "sch_c" : [] , "stations":[] }

all_res = {'compatible':comp, 'uncertain': un, 'wrong': wrong }


proc = {}

for uw in unique_wmo:
    
    # unspecified WMO sensor 
    if uw in ['90']: 
        continue
                
    dfr = df.loc[df.wmo == uw ]
    print('\nLoop on WMO code ' , uw )
    
    sch = dfr.sch.values

    if not uw in proc.keys():
        proc[uw] = []
    
    # extracting and looping over all unique Sch ids for same WMO id
    unique_sch = np.unique(sch)
    
    cont_all=True

    break_all = False
    
    if not break_all:
        pass
    else:
        break
        
    while cont_all:
        
        for us in unique_sch:

            
            dfrs = dfr.loc[dfr.sch == us]
            dfrs = dfrs [['sch','sch_comment','wmo','wmo_comment','station']]
            
            print("Matching Schroeder sensors" , dfrs )
            
            sch_com = dfrs.sch_comment.values
            wmo_com = dfrs.wmo_comment.values
            stat = dfrs.station.values
                            
            print(uw, us, cont_all)
                

            cs =  dfrs.sch_comment.values[0]
            cw = dfrs.wmo_comment.values[0]
            st = list(dfrs.station.values)
            
            # set of unidentified sensor, too little information           
            if us in ['XXX','???','RPDE','ZP9','V??','M__','W_R','ZMW' , 'W__' , 'ZKG', 'RZD']: 
                continue

            print('********************************\n')
            print('Station: ', st , '\n')
            print('Sch: ' , us )
            print('WMO: ' , uw )
            print('Sch com:' , cs , '\n' , "WMO com: " , cw)
            print('\n********************************')

            ask = "Compatible? Type Y(compatible) N(wrong) U(uncertain).  \n *** 'E' for ending analysis ---> "

            # interactive shell question 
            answer = input(ask)


            def update_res(d, st, s, cs, w, cw):
                """ Update the result dictionaries according to user answer """
                proc[w].append(s)

                all_res[d]['wmo'].append(w)
                all_res[d]['wmo_c'].append(cw)
                all_res[d]['sch'].append(s)
                all_res[d]['sch_c'].append(cs)

                return True

            # terminating option
            if answer == 'E':
                #cont_all=False
                #break_all= True
                np.save('dictionary_sensor_mapping', all_res, allow_pickle=True)
                sys.exit()
                break

            elif answer in ['Y','y','yes','Yes','YES']:
                d = 'compatible'
                cont = update_res(d, st, us, cs, uw, cw)
                continue

            elif answer in ['n','N','no','No',"NO"]:
                d = 'wrong'
                cont = update_res(d, st, us, cs, uw, cw)
                continue

            elif answer in ['u', "U", 'un', 'Un', 'UN']:
                d = 'uncertain'
                cont = update_res(d, st, us, cs, uw, cw)
                continue

            else:
                print('Not valid option!')
                cont= True
                answer = input(ask)
                if answer == 'E': 
                    #cont_all=False
                    #break_all= True
                    np.save('dictionary_sensor_mapping', all_res, allow_pickle=True)
                    sys.exit()
                    break
                    
        cont_all = False                                 

    
print(all_res)

np.save('dictionary_sensor_mapping', all_res, allow_pickle=True)
