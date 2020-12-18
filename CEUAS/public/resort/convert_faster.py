#!/usr/bin/env
# coding: utf-8

import numpy
import numpy as np
import pandas
import pandas as pd
from numba import njit
import sys,glob
import zipfile, os, time
import urllib3
from datetime import datetime, timedelta
import glob
import h5py
sys.path.append(os.getcwd()+'/../cds-backend/code/')
sys.path.append(os.getcwd()+'/../harvest/code/')
from harvest_convert_to_netCDF_newfixes import write_dict_h5
import cds_eua3 as eua
eua.logging_set_level(30)
import xarray as xr
import cdsapi, zipfile, os, time
import copy
from shutil import copyfile
import multiprocessing
sys.path.append(os.getcwd()+'/../resort/rasotools-master/')
import rasotools
import warnings
warnings.filterwarnings('ignore')

opath='/raid60/scratch/uli/converted_v2/'
# if there are nan values in the pressure level - we will just sort without any converting!
def do_resort(fn):
    targetfile = opath+fn.split('/')[-1]  
    
    with h5py.File(fn, 'r') as file:
        with h5py.File(targetfile, 'w') as newfile:
            groups = []
            for i in file.keys():
                if type(file[i]) == h5py._hl.group.Group:
                    newfile.create_group(i)
                    groups.append(i)
                elif i == 'recordindex' or i == 'recordtimestamp':
                    pass
                else:
                    newfile.create_dataset(i, data=file[i][:])
            for i in groups:
                if(i == 'recordindices' or i == 'observations_table' or i == 'era5fb'):
                    pass
                else:
                    for j in file[i].keys():
                        newfile[i].create_dataset(j, data=file[i][j][:])
    
    data =  eua.CDMDataset(fn)
    allvars = data.observations_table.observed_variable[()]
    allvars.sort()
    allvars = numpy.unique(allvars)
    #
    ri = data.recordindex[()]
#     print('recordindex: ', len(ri))
    rt = data.recordtimestamp[()]
    keys = data.observations_table.keys()[:-1]
    fbkeys = data.era5fb.keys()[:-1]
    # dropping all keys, where dimensions won't work - just help variabels for dimensions
    pops = []
    for i in range(len(keys)):
        if 'string' in keys[i]:
            pops.append(keys[i])
    for i in pops: keys.remove(i)
    pops = []
    for i in range(len(fbkeys)):
        if 'string' in fbkeys[i]:
            pops.append(fbkeys[i])
    for i in pops: fbkeys.remove(i)

    recordindices = [[] for i in range(len(allvars))]
    recordtimestamps = [[] for i in range(len(allvars))]

    # output variables (from observations_table)
    ov = []
    for o in keys:
        ov.append([[] for i in range(len(allvars))])
    fb = []
    for o in fbkeys:
        fb.append([[] for i in range(len(allvars))])
    #
    # loading the observed_variables
    #
    obsv = data.observations_table.observed_variable[:]
    #
    # resorting the data
    #
#     print('resort:start')
    @njit
    def make_vrindex(vridx,ridx,idx):
        l=0
        for i in range(1,len(idx)): # to set the recordindices
            if ridx[i]>ridx[i-1]:
                vridx[ridx[i-1]:ridx[i]]=l # next record after l
                l=i
        vridx[ridx[i]:]=len(idx) # next record for the last element is the len of the data


    tt=time.time()

    ridxall=np.zeros(obsv.shape[0],dtype=np.int64) # reverse index - index of the record index
    j=-1
    for j in range(len(ri)-1):
        ridxall[ri[j]:ri[j+1]]=j
    j+=1
    ridxall[ri[j]:]=j # for the last elemenet
    ridx=[]
    vridx=[]
    absidx=[]
    abscount=0
    for j in range(len(allvars)):
        idx=np.where(obsv==allvars[j])[0] # index of all elements form certain variable j
#         print(j,len(idx),',',end='')
        vridx.append(np.zeros(ri.shape[0],dtype=np.int64)) # all zeros in lenght of record index
        ridx=ridxall[idx] # ridxall where variable is j
        make_vrindex(vridx[-1],ridx,idx)
        vridx[-1]+=abscount # abscount for stacking the recordindex

        absidx.append(copy.copy(idx)) # why copy? - to make sure it's not just the ref. - maybe ok without the cp
        abscount+=len(idx)

#     print('')
    #
    # finishing the sorting 
    #
    absidx=np.concatenate(absidx)
    #
    # recordtimestamps are only necessary once
    #
    recordtimestamps = recordtimestamps[0]
    #
    # targetfile has to be a copy of the original file
    #
#     targetfile = '/raid60/scratch/uli/converted_v2/'+fn.split('/')[-1]# 0-20000-0-63894_CEUAS_merged_v0.nc'
    if os.path.isfile(targetfile):
        mode='r+'
    else:
        mode='w'
#     print()
#     print('writing '+targetfile)

    for i in range(len(keys)):
        ov_vars = data.observations_table[keys[i]][:]
        ov_vars = ov_vars[absidx]
        if keys[i] == 'index':
            pass
        elif keys[i] == 'observation_id' or keys[i] == 'report_id' or keys[i] == 'sensor_id' or keys[i] == 'source_id':
            alldict = {keys[i]:np.asarray(ov_vars, dtype='S1')}
            write_dict_h5(targetfile, alldict, 'observations_table', {keys[i]: { 'compression': 'gzip' } }, [keys[i]])
        else:
            alldict = pandas.DataFrame({keys[i]:ov_vars})
            write_dict_h5(targetfile, alldict, 'observations_table', {keys[i]: { 'compression': 'gzip' } }, [keys[i]])  

    for i in range(len(fbkeys)):
        fb_vars = data.era5fb[fbkeys[i]][:]
        fb_vars = fb_vars[absidx]
        if fbkeys[i] == 'index' or fbkeys[i] == 'string6' or fbkeys[i] == 'string7' or fbkeys[i] == 'string10':
            pass
        elif fbkeys[i] == 'expver' or fbkeys[i] == 'source@hdr' or fbkeys[i] == 'source_id' or fbkeys[i] == 'statid@hdr':
            alldict = {fbkeys[i]:np.asarray(fb_vars, dtype='S1')}
            write_dict_h5(targetfile, alldict, 'era5fb', {fbkeys[i]: { 'compression': 'gzip' } }, [fbkeys[i]])
        else:
            alldict = pandas.DataFrame({fbkeys[i]:fb_vars})
            write_dict_h5(targetfile, alldict, 'era5fb', {fbkeys[i]: { 'compression': 'gzip' } }, [fbkeys[i]]) 
    #
    # writing the recordindices and recordtimestamp.
    #       
    recordindices=vridx
    for i in range(len(recordindices)):
        testvar = pandas.DataFrame({str(allvars[i]):recordindices[i]})
        write_dict_h5(targetfile, testvar, 'recordindices', {str(allvars[i]): { 'compression': None } }, [str(allvars[i])]) 

    write_dict_h5(targetfile, {'recordtimestamp':rt}, 'recordindices', {'recordtimestamp': { 'compression': None } }, ['recordtimestamp'])

    print('elapsed:',time.time()-tt)



def convert_missing(fn, destination: str = opath):
    tt=time.time()
    nanlist = [float('nan'), np.nan, 0, -2147483648]
    with eua.CDMDataset(fn) as data:
        arrayconverter = data.to_dataframe(groups='observations_table', variables=['observed_variable'])
        arrayconverter = arrayconverter.observed_variable.head(1).to_xarray()
        rt = data.recordtimestamp[:]
        keys = data.observations_table.keys()
        keys = [x for x in keys if not x.startswith('string')]
        keys.remove('index')
        keys.remove('shape')
        obskeys = keys
        obstab_writetofile = [[] for i in range(len(obskeys))]
        
        keys = data.era5fb.keys()
        keys = [x for x in keys if not x.startswith('string')]
        keys.remove('index')
        keys.remove('shape')
        fg_depar = keys.index('fg_depar@body')
        depar = keys.index('an_depar@body')
        biascorr = keys.index('biascorr@body')
        fg_biascorr = keys.index('biascorr_fg@body')
        fbkeys = keys
        fb_writetofile = [[] for i in range(len(fbkeys))]
        
        recidxlen = len(data.recordindex[:])
        
        addtorecordindex = [] # will be filled with the count of how many variables have been added before !! addtorecordindex[0] has to be added to recordindex [1] !!
        addedvarscount = 0 # will grow with every variable added
        onlyone = False
        if recidxlen == 1:
            recidxlen = 2
            onlyone = True
        
        # loading data:
        loaded_data = [[]]*len(obskeys)
        for o in range(len(obskeys)):
            loaded_data[o] = np.asarray(data.observations_table[obskeys[o]][:])
            
        loaded_fb = [[]]*len(fbkeys)
        for o in range(len(fbkeys)):
            loaded_fb[o] = np.asarray(data.era5fb[fbkeys[o]][:])
            
        recordindex = data.recordindex[:]
        # --->
    tt=time.time()
    tttp=0.
    for i in range(recidxlen-1):
        obstab = [[]]*len(obskeys)
        fb = [[]]*len(fbkeys)
#         if i%10==0 :
#             print(i,time.time()-tt)

        for o in range(len(obskeys)):
            if onlyone:
                obstab[o] = loaded_data[o][recordindex[i]:]
            else:
                obstab[o] = loaded_data[o][recordindex[i]:recordindex[i+1]]
            if obskeys[o] == 'observed_variable':
                obsvar = obstab[o]
            elif obskeys[o] == 'z_coordinate':
                plev = obstab[o]
            elif obskeys[o] == 'observation_value':
                obsval = obstab[o]
            elif obskeys[o] == 'z_coordinate_type':
                plevtype = obstab[o]
                
        if np.isnan(plev).any():
            do_resort(fn)
            return

        for o in range(len(fbkeys)):
            if onlyone:
                fb[o] = loaded_fb[o][recordindex[i]:]
            else:
                fb[o] = loaded_fb[o][recordindex[i]:recordindex[i+1]]

        raso_t = arrayconverter.copy()
        raso_p = arrayconverter.copy()
        raso_raso_rh = arrayconverter.copy()
        raso_sh = arrayconverter.copy()
        raso_rh = arrayconverter.copy()
        raso_dep_t = arrayconverter.copy()
        raso_dep_rh = arrayconverter.copy()
        raso_dep_sh = arrayconverter.copy()
        raso_dpd = arrayconverter.copy()
        for j in np.unique(plev):
            convertedfrom = []
            #select = plev == j
            select=np.where(plev==j)[0]
            obsvr = obsvar[select]
            obsvl = obsval[select]
            ptype = plevtype[select]
            
            #
            # Converting to dewpointdepression
            #
            valid_conversion_found = False

            if not (34 in obsvr) and (38 in obsvr) and not valid_conversion_found:
                if (85 in obsvr) and (ptype[obsvr == 85] == 1):
                    raso_t.values = obsvl[obsvr == 85]
                    #t = arrayconverter.copy()
                    raso_p.values = np.array([j])
                    #p = arrayconverter.copy()
                    raso_rh.values = obsvl[obsvr == 38]
                    #rh = arrayconverter.copy()
                    dpd = rasotools.met.convert.to_dpd(temp=raso_t,press=raso_p,rel_humi=raso_rh)
                    
                    if not np.isnan(dpd).any():

                        idx=np.where(obsvr==85)[0][0]
                        iselect=select[idx]
                        for o in range(len(obskeys)):
                            # write the new variable everywhere into the observationstable
                            if obskeys[o] == 'observed_variable':
                                obstab_writetofile[o].append(34)
                            elif obskeys[o] == 'observation_value':
                                obstab_writetofile[o].append(dpd[0])
                            elif obskeys[o] == 'conversion_method':
                                obstab_writetofile[o].append(3) # 'from_relative_humidity'
                            elif obskeys[o] == 'conversion_flag':
                                obstab_writetofile[o].append(0)
                            else:
#                                obstab_writetofile[o].append(obstab[o][select][obsvr == 85][0])
                                obstab_writetofile[o].append(obstab[o][iselect])

                        for o in range(len(fbkeys)):
                            # write the new variable everywhere into the fb
                            #if ((o == depar) or (o == fg_depar)) and (not (fb[o][select][obsvr == 85][0] in nanlist)) and (not (fb[o][select][obsvr == 38][0] in nanlist)):
                            if ((o == depar) or (o == fg_depar)) and (not (fb[o][iselect] in nanlist)) and (not (fb[o][select][obsvr == 38][0] in nanlist)):
                                raso_dep_t.values = obsvl[obsvr == 85] - fb[o][iselect]
                                #dep_t = arrayconverter.copy()
                                raso_dep_rh.values = obsvl[obsvr == 38] - fb[o][select][obsvr == 38]
                                #dep_rh = arrayconverter.copy()
                                fb_writetofile[o].append(dpd.values[0] - rasotools.met.convert.to_dpd(temp=raso_dep_t,press=raso_p,rel_humi=raso_dep_rh).values[0])
                            elif ((o == depar) or (o == fg_depar) or (o == biascorr) or (o == fg_biascorr)):
                                fb_writetofile[o].append(float('nan'))
                            else:
                                fb_writetofile[o].append(fb[o][iselect])
                        addedvarscount += 1
                        convertedfrom.append(38)
                        valid_conversion_found = True


            if not (34 in obsvr) and (39 in obsvr) and not valid_conversion_found:
                if (85 in obsvr) and (ptype[obsvr == 85] == 1): 
                    raso_t.values = obsvl[obsvr == 85]
#                     t = arrayconverter.copy()
                    raso_p.values = [j]
#                     p = arrayconverter.copy()
                    raso_sh.values = obsvl[obsvr == 39]
#                     sh = arrayconverter.copy()
                    dpd = rasotools.met.convert.to_dpd(temp=raso_t,press=raso_p,spec_humi=raso_sh)
                    
                    if not np.isnan(dpd).any():
                        
                        idx=np.where(obsvr==85)[0][0]
                        iselect=select[idx]
                        for o in range(len(obskeys)):
                            # write the new variable everywhere into the observationstable
                            if obskeys[o] == 'observed_variable':
                                obstab_writetofile[o].append(34)
                            elif obskeys[o] == 'observation_value':
                                obstab_writetofile[o].append(dpd[0])
                            elif obskeys[o] == 'conversion_method':
                                obstab_writetofile[o].append(4) # 'from_specific_humidity'
                            elif obskeys[o] == 'conversion_flag':
                                obstab_writetofile[o].append(0)
                            else:
                                obstab_writetofile[o].append(obstab[o][iselect])

#                         dep_t = arrayconverter.copy()
#                         dep_sh = arrayconverter.copy()
                        for o in range(len(fbkeys)):
                            # write the new variable everywhere into the fb
                            if ((o == depar) or (o == fg_depar)) and (not (fb[o][iselect] in nanlist)) and (not (fb[o][select][obsvr == 39][0] in nanlist)):
                                raso_dep_t.values = obsvl[obsvr == 85] - fb[o][iselect]
                                #dep_t = arrayconverter.copy()
                                raso_dep_sh.values = obsvl[obsvr == 39][0] - fb[o][select][obsvr == 39]
                                #dep_sh = arrayconverter.copy()
                                fb_writetofile[o].append(dpd.values[0] - rasotools.met.convert.to_dpd(temp=raso_dep_t,press=raso_p,spec_humi=raso_dep_sh).values[0])
                            elif ((o == depar) or (o == fg_depar) or (o == biascorr) or (o == fg_biascorr)):
                                fb_writetofile[o].append(float('nan'))
                            else:
                                fb_writetofile[o].append(fb[o][iselect])
                        addedvarscount += 1
                        convertedfrom.append(39)
                        valid_conversion_found = True

            if not (34 in obsvr) and (36 in obsvr) and not valid_conversion_found:
                if 85 in obsvr: 
                    t = obsvl[obsvr == 85]
                    dp = obsvl[obsvr == 36]
                    dpd = t-dp
                    
                    if not np.isnan(dpd).any():

                        for o in range(len(obskeys)):
                            # write the new variable everywhere into the observationstable
                            if obskeys[o] == 'observed_variable':
                                obstab_writetofile[o].append(34)
                            elif obskeys[o] == 'observation_value':
                                obstab_writetofile[o].append(dpd[0])
                            elif obskeys[o] == 'conversion_method':
                                obstab_writetofile[o].append(2) # 'from_dewpoint'
                            elif obskeys[o] == 'conversion_flag':
                                obstab_writetofile[o].append(0)
                            else:
                                obstab_writetofile[o].append(obstab[o][iselect])

                        for o in range(len(fbkeys)):
                            # write the new variable everywhere into the fb
                            if ((o == depar) or (o == fg_depar) or (o == biascorr) or (o == fg_biascorr)):
                                fb_writetofile[o].append(float('nan'))
                            else:
                                fb_writetofile[o].append(fb[o][iselect])
                        addedvarscount += 1
                        convertedfrom.append(36)
                        valid_conversion_found = True

            #
            # Converting to relative humidity
            #
            valid_conversion_found = False

            if not (38 in obsvr) and (39 in obsvr) and not valid_conversion_found:
                if (85 in obsvr) and (ptype[obsvr == 85] == 1): 
                    raso_t.values = obsvl[obsvr == 85]
#                     t = arrayconverter.copy()
                    raso_p.values = [j]
#                     p = arrayconverter.copy()
                    raso_sh.values = obsvl[obsvr == 39]
#                     sh = arrayconverter.copy()
                    rh = rasotools.met.convert.to_rh(temp=raso_t, spec_humi=raso_sh, press=raso_p)
                    
                    if not np.isnan(rh).any():

                        for o in range(len(obskeys)):
                            # write the new variable everywhere into the observationstable
                            if obskeys[o] == 'observed_variable':
                                obstab_writetofile[o].append(38)
                            elif obskeys[o] == 'observation_value':
                                obstab_writetofile[o].append(rh[0])
                            elif obskeys[o] == 'conversion_method':
                                obstab_writetofile[o].append(4) # 'from_specific_humidity'
                            elif obskeys[o] == 'conversion_flag':
                                obstab_writetofile[o].append(0)
                            else:
                                obstab_writetofile[o].append(obstab[o][iselect])

#                         dep_t = arrayconverter.copy()
#                         dep_sh = arrayconverter.copy()
                        for o in range(len(fbkeys)):
                            # write the new variable everywhere into the fb
                            if ((o == depar) or (o == fg_depar)) and (not (fb[o][iselect] in nanlist)) and (not (fb[o][select][obsvr == 39][0] in nanlist)):
                                raso_dep_t.values = obsvl[obsvr == 85] - fb[o][iselect]
                                #dep_t = arrayconverter.copy()
                                raso_dep_sh.values = obsvl[obsvr == 39] - fb[o][select][obsvr == 39]
                                #dep_sh = arrayconverter.copy()
                                fb_writetofile[o].append(rh.values[0] - rasotools.met.convert.to_rh(temp=raso_dep_t,press=raso_p,spec_humi=raso_dep_sh).values[0])
                            elif ((o == depar) or (o == fg_depar) or (o == biascorr) or (o == fg_biascorr)):
                                fb_writetofile[o].append(float('nan'))
                            else:
                                fb_writetofile[o].append(fb[o][iselect])
                        addedvarscount += 1
                        convertedfrom.append(39)
                        valid_conversion_found = True

            if not (38 in obsvr) and (36 in obsvr) and not valid_conversion_found:
                if (85 in obsvr): 
                    t = obsvl[obsvr == 85]
                    dp = obsvl[obsvr == 36]
                    raso_dpd.values = t-dp
#                     arrayconverter.values = dpd
#                     dpd = arrayconverter
                    raso_t.values = t
#                     t = arrayconverter
                    rh = rasotools.met.convert.to_rh(temp=raso_t,dpd=raso_dpd)
                    
                    if not np.isnan(rh).any():

                        for o in range(len(obskeys)):
                            # write the new variable everywhere into the observationstable
                            if obskeys[o] == 'observed_variable':
                                obstab_writetofile[o].append(38)
                            elif obskeys[o] == 'observation_value':
                                obstab_writetofile[o].append(rh[0])
                            elif obskeys[o] == 'conversion_method':
                                obstab_writetofile[o].append(2) # 'from_dewpoint'
                            elif obskeys[o] == 'conversion_flag':
                                obstab_writetofile[o].append(0)
                            else:
                                obstab_writetofile[o].append(obstab[o][iselect])

                        for o in range(len(fbkeys)):
                            # write the new variable everywhere into the fb
                            if ((o == depar) or (o == fg_depar) or (o == biascorr) or (o == fg_biascorr)):
                                fb_writetofile[o].append(float('nan'))
                            else:
                                fb_writetofile[o].append(fb[o][iselect])
                        addedvarscount += 1 
                        convertedfrom.append(36)
                        valid_conversion_found = True

            if not (38 in obsvr) and (34 in obsvr) and not valid_conversion_found:
                if (85 in obsvr): 
                    t = obsvl[obsvr == 85]
                    raso_t.values = t
#                     t = arrayconverter.copy()
                    raso_dpd.values = t - obsvl[obsvr == 34]
#                     dpd = arrayconverter.copy()
                    rh = rasotools.met.convert.to_rh(temp=raso_t, dpd=raso_dpd)
                    
                    if not np.isnan(rh).any():

                        for o in range(len(obskeys)):
                            # write the new variable everywhere into the observationstable
                            if obskeys[o] == 'observed_variable':
                                obstab_writetofile[o].append(38)
                            elif obskeys[o] == 'observation_value':
                                obstab_writetofile[o].append(rh[0])
                            elif obskeys[o] == 'conversion_method':
                                obstab_writetofile[o].append(3) # 'from_dewpointdepression'
                            elif obskeys[o] == 'conversion_flag':
                                obstab_writetofile[o].append(0)
                            else:
                                obstab_writetofile[o].append(obstab[o][iselect])

                        for o in range(len(fbkeys)):
                            # write the new variable everywhere into the fb
                            if ((o == depar) or (o == fg_depar) or (o == biascorr) or (o == fg_biascorr)):
                                fb_writetofile[o].append(float('nan'))
                            else:
                                fb_writetofile[o].append(fb[o][iselect])
                        addedvarscount += 1
                        convertedfrom.append(34)
                        valid_conversion_found = True

            #
            # Converting to specific humidity
            #
            valid_conversion_found = False

            if not (39 in obsvr) and (38 in obsvr) and not valid_conversion_found:
                if (85 in obsvr) and (ptype[obsvr == 85] == 1): 
                    raso_t.values = obsvl[obsvr == 85]
#                     t = arrayconverter.copy()
                    raso_p.values = [j]
#                     p = arrayconverter.copy()
                    raso_rh.values = obsvl[obsvr == 38]
#                     rh = arrayconverter.copy()
                    sh = rasotools.met.convert.to_sh(temp=raso_t, press=raso_p, rel_humi=raso_rh)
                    
                    if not np.isnan(sh).any():

                        for o in range(len(obskeys)):
                            # write the new variable everywhere into the observationstable
                            if obskeys[o] == 'observed_variable':
                                obstab_writetofile[o].append(39)
                            elif obskeys[o] == 'observation_value':
                                obstab_writetofile[o].append(sh[0])
                            elif obskeys[o] == 'conversion_method':
                                obstab_writetofile[o].append(2) # 'from_relative_humidity'
                            elif obskeys[o] == 'conversion_flag':
                                obstab_writetofile[o].append(0)
                            else:
                                obstab_writetofile[o].append(obstab[o][iselect])

                        dep_t = arrayconverter.copy()
                        dep_rh = arrayconverter.copy()
                        for o in range(len(fbkeys)):
                            # write the new variable everywhere into the fb
                            if ((o == depar) or (o == fg_depar)) and (not (fb[o][iselect] in nanlist)) and (not (fb[o][select][obsvr == 38][0] in nanlist)):
                                raso_dep_t.values = obsvl[obsvr == 85] - fb[o][iselect]
                                #dep_t = arrayconverter.copy()
                                raso_dep_rh.values = obsvl[obsvr == 38] - fb[o][select][obsvr == 38]
                                #dep_rh = arrayconverter.copy()
                                fb_writetofile[o].append(sh.values[0] - rasotools.met.convert.to_sh(temp=raso_dep_t, press=raso_p, rel_humi=raso_dep_rh).values[0])
                            elif ((o == depar) or (o == fg_depar) or (o == biascorr) or (o == fg_biascorr)):
                                fb_writetofile[o].append(float('nan'))
                            else:
                                fb_writetofile[o].append(fb[o][iselect])
                        addedvarscount += 1   
                        convertedfrom.append(38)
                        valid_conversion_found = True

            if not (39 in obsvr) and (34 in obsvr) and not valid_conversion_found:
                if (85 in obsvr) and (ptype[obsvr == 85] == 1): 
                    raso_t.values = obsvl[obsvr == 85]
#                     t = arrayconverter.copy()
                    raso_dpd.values = obsvl[obsvr == 34]
#                     dpd = arrayconverter.copy()
                    raso_p.values = np.array([j])
#                     p = arrayconverter.copy()
                    sh = rasotools.met.convert.to_sh(dpd=raso_dpd, press=raso_p, temp=raso_t)
                    
                    if not np.isnan(sh).any():

                        for o in range(len(obskeys)):
                            # write the new variable everywhere into the observationstable
                            if obskeys[o] == 'observed_variable':
                                obstab_writetofile[o].append(39)
                            elif obskeys[o] == 'observation_value':
                                obstab_writetofile[o].append(rh[0])
                            elif obskeys[o] == 'conversion_method':
                                obstab_writetofile[o].append(3) # 'from_dewpointdepression'
                            elif obskeys[o] == 'conversion_flag':
                                obstab_writetofile[o].append(0)
                            else:
                                obstab_writetofile[o].append(obstab[o][iselect])

                        for o in range(len(fbkeys)):
                            # write the new variable everywhere into the fb
                            if ((o == depar) or (o == fg_depar) or (o == biascorr) or (o == fg_biascorr)):
                                fb_writetofile[o].append(float('nan'))
                            else:
                                fb_writetofile[o].append(fb[o][iselect])
                        addedvarscount += 1
                        convertedfrom.append(34)

            if not (39 in obsvr) and (36 in obsvr) and not valid_conversion_found:
                if (85 in obsvr) and (ptype[obsvr == 85] == 1): 
                    t = obsvl[obsvr == 85]
                    dp = obsvl[obsvr == 36]
                    raso_dpd.values = t-dp
#                     dpd = arrayconverter.copy()
                    raso_t.values = t
#                     t = arrayconverter.copy()
                    raso_p.values = [j]
#                     p = arrayconverter.copy()
                    sh = rasotools.met.convert.to_sh(dpd=raso_dpd, press=raso_p, temp=raso_t)
                    
                    if not np.isnan(sh).any():

                        for o in range(len(obskeys)):
                            # write the new variable everywhere into the observationstable
                            if obskeys[o] == 'observed_variable':
                                obstab_writetofile[o].append(39)
                            elif obskeys[o] == 'observation_value':
                                obstab_writetofile[o].append(rh[0])
                            elif obskeys[o] == 'conversion_method':
                                obstab_writetofile[o].append(4) # 'from_dewpoint'
                            elif obskeys[o] == 'conversion_flag':
                                obstab_writetofile[o].append(0)
    #                             elif obskeys[o] == 'observation_id':
    #                                 obstab_writetofile[o].append(np.asarray([b'']))
                            else:
                                obstab_writetofile[o].append(obstab[o][iselect])

                        for o in range(len(fbkeys)):
                            # write the new variable everywhere into the fb
                            if ((o == depar) or (o == fg_depar) or (o == biascorr) or (o == fg_biascorr)):
                                fb_writetofile[o].append(float('nan'))
                            else:
                                fb_writetofile[o].append(fb[o][iselect])
                        addedvarscount += 1
                        convertedfrom.append(36)
                        valid_conversion_found = True

            #
            # Converting to wind components
            #

            if not (104 in obsvr) or not (105 in obsvr):
                if 106 in obsvr and 107 in obsvr:
                    wd = obsvl[obsvr == 106]
                    ws = obsvl[obsvr == 107]
                    u = ws * np.cos(np.radians(wd))
                    v = ws * np.sin(np.radians(wd))

                    if not (104 in obsvr) and not np.isnan(u).any():
                        for o in range(len(obskeys)):
                            # write the new variable everywhere into the observationstable
                            if obskeys[o] == 'observed_variable':
                                obstab_writetofile[o].append(104)
                            elif obskeys[o] == 'observation_value':
                                obstab_writetofile[o].append(u[0])
                            elif obskeys[o] == 'conversion_method':
                                obstab_writetofile[o].append(2) # 'from_speed_and_direction'
                            elif obskeys[o] == 'conversion_flag':
                                obstab_writetofile[o].append(0)
#                                 elif obskeys[o] == 'observation_id':
#                                     obstab_writetofile[o].append(np.asarray([b'']))
                            else:
                                obstab_writetofile[o].append(obstab[o][select][obsvr == 106][0])

                        for o in range(len(fbkeys)):
                            # write the new variable everywhere into the fb
                            if ((o == depar) or (o == fg_depar) or (o == biascorr) or (o == fg_biascorr)):
                                fb_writetofile[o].append(float('nan'))
                            else:
                                fb_writetofile[o].append(fb[o][select][obsvr == 106][0])
                        addedvarscount += 1

                    if not (105 in obsvr) and not np.isnan(v).any():
                        for o in range(len(obskeys)):
                            # write the new variable everywhere into the observationstable
                            if obskeys[o] == 'observed_variable':
                                obstab_writetofile[o].append(105)
                            elif obskeys[o] == 'observation_value':
                                obstab_writetofile[o].append(v[0])
                            elif obskeys[o] == 'conversion_method':
                                obstab_writetofile[o].append(2) # 'from_speed_and_direction'
                            elif obskeys[o] == 'conversion_flag':
                                obstab_writetofile[o].append(0)
#                                 elif obskeys[o] == 'observation_id':
#                                     obstab_writetofile[o].append(np.asarray([b'']))
                            else:
                                obstab_writetofile[o].append(obstab[o][select][obsvr == 106][0])

                        for o in range(len(fbkeys)):
                            # write the new variable everywhere into the fb
                            if ((o == depar) or (o == fg_depar) or (o == biascorr) or (o == fg_biascorr)):
                                fb_writetofile[o].append(float('nan'))
                            else:
                                fb_writetofile[o].append(fb[o][select][obsvr == 106][0])
                        addedvarscount += 1
                    convertedfrom.append(106)
                    convertedfrom.append(107)

            #
            # Converting to windspeed and winddirection
            #

            if not (106 in obsvr) or not (107 in obsvr):
                if 104 in obsvr and 105 in obsvr:
                    u = obsvl[obsvr == 104][0]
                    v = obsvl[obsvr == 105][0]
                    ws = np.sqrt(u ** 2 + v ** 2)
                    wd = 90 - np.arctan2(-v, -u) * 180 / np.pi - 180.
                    wd = np.where(wd > 0., wd, 360.+wd)

                    if not (106 in obsvr) and not np.isnan(wd).any():
                        for o in range(len(obskeys)):
                            # write the new variable everywhere into the observationstable
                            if obskeys[o] == 'observed_variable':
                                obstab_writetofile[o].append(106)
                            elif obskeys[o] == 'observation_value':
                                obstab_writetofile[o].append(wd)
                            elif obskeys[o] == 'conversion_method':
                                obstab_writetofile[o].append(2) # 'from_wind_components'
                            elif obskeys[o] == 'conversion_flag':
                                obstab_writetofile[o].append(0)
#                                 elif obskeys[o] == 'observation_id':
#                                     obstab_writetofile[o].append(np.asarray([b'']))
                            else:
                                obstab_writetofile[o].append(obstab[o][select][obsvr == 104][0])

                        for o in range(len(fbkeys)):
                            # write the new variable everywhere into the fb
                            if ((o == depar) or (o == fg_depar)) and (not (fb[o][select][obsvr == 104][0] in nanlist)) and (not (fb[o][select][obsvr == 105][0] in nanlist)):
                                dep_u = u - fb[o][select][obsvr == 104][0]
                                dep_v = v - fb[o][select][obsvr == 105][0]
                                dep_wd = 90 - np.arctan2(-dep_v, -dep_u) * 180 / np.pi - 180.
                                dep_wd = np.where(dep_wd > 0., dep_wd, 360.+dep_wd)
                                fb_writetofile[o].append(float(wd - dep_wd))
                            elif ((o == depar) or (o == fg_depar) or (o == biascorr) or (o == fg_biascorr)):
                                fb_writetofile[o].append(float('nan'))
                            else:
                                fb_writetofile[o].append(fb[o][select][obsvr == 104][0])
                        addedvarscount += 1

                    if not (107 in obsvr) and not np.isnan(ws).any():
                        for o in range(len(obskeys)):
                            # write the new variable everywhere into the observationstable
                            if obskeys[o] == 'observed_variable':
                                obstab_writetofile[o].append(107)
                            elif obskeys[o] == 'observation_value':
                                obstab_writetofile[o].append(ws)
                            elif obskeys[o] == 'conversion_method':
                                obstab_writetofile[o].append(2) # 'from_wind_components'
                            elif obskeys[o] == 'conversion_flag':
                                obstab_writetofile[o].append(0)
#                                 elif obskeys[o] == 'observation_id':
#                                     obstab_writetofile[o].append(np.asarray([b'']))
                            else:
                                obstab_writetofile[o].append(obstab[o][select][obsvr == 104][0])

                        for o in range(len(fbkeys)):
                            # write the new variable everywhere into the fb
                            if ((o == depar) or (o == fg_depar)) and (not (fb[o][select][obsvr == 104][0] in nanlist)) and (not (fb[o][select][obsvr == 105][0] in nanlist)):
                                dep_u = u - fb[o][select][obsvr == 104][0]
                                dep_v = v - fb[o][select][obsvr == 105][0]
                                dep_ws = np.sqrt(dep_u ** 2 + dep_v ** 2)
                                fb_writetofile[o].append(float(ws - dep_ws))
                            if ((o == depar) or (o == fg_depar) or (o == biascorr) or (o == fg_biascorr)):
                                fb_writetofile[o].append(float('nan'))
                            else:
                                fb_writetofile[o].append(fb[o][select][obsvr == 104][0])
                        addedvarscount += 1
                    convertedfrom.append(104)
                    convertedfrom.append(105)

            #
            # add all the non converted variables into the list
            #
            ttt=time.time()
            for o in range(len(obskeys)):
                # write everything what was already in the file
                if obskeys[o] == 'conversion_flag' and  len(convertedfrom) > 0:
                    wrt = obstab[o][select]
                    for pp in convertedfrom:
                        wrt[obsvr == pp] = 0
                    obstab_writetofile[o].extend(wrt)
                else:
                    obstab_writetofile[o].extend(obstab[o][select])
            for o in range(len(fbkeys)):
                # write everything what was already in the file
                fb_writetofile[o].extend(fb[o][select])

        # adjusting the recordindex
#             tttp+=time.time()-ttt
        addtorecordindex.append(addedvarscount)
#         if i%10==0:
#             print('writetofile',tttp)
#             tttp=0.
        
    ri = recordindex
    ri = np.asarray([0] + (addtorecordindex)) + np.asarray(ri) # [0] before the array!
          
#     for ow in range(len(obstab_writetofile)):
#         print(obskeys[ow], len(obstab_writetofile[ow]))
    
    # sorting:
    print('start sorting')
    targetfile = destination+fn.split('/')[-1]  
    
    with h5py.File(fn, 'r') as file:
        with h5py.File(targetfile, 'w') as newfile:
            groups = []
            for i in file.keys():
                if type(file[i]) == h5py._hl.group.Group:
                    newfile.create_group(i)
                    groups.append(i)
                elif i == 'recordindex' or i == 'recordtimestamp':
                    pass
                else:
                    newfile.create_dataset(i, data=file[i][:])
            for i in groups:
                if(i == 'recordindices' or i == 'observations_table' or i == 'era5fb'):
                    pass
                else:
                    for j in file[i].keys():
                        newfile[i].create_dataset(j, data=file[i][j][:])
    
#     data =  eua.CDMDataset(fn)
    allvars = np.asarray(obstab_writetofile[obskeys.index('observed_variable')])
    allvars.sort()
    allvars = numpy.unique(allvars)
    #
#     ri = data.recordindex[()]
#     print('recordindex: ', len(ri))
#     rt = data.recordtimestamp[()]
    keys = obskeys # data.observations_table.keys()[:-1]
    fbkeys = fbkeys # data.era5fb.keys()[:-1]
    # dropping all keys, where dimensions won't work - just help variabels for dimensions
    pops = []
    for i in range(len(keys)):
        if 'string' in keys[i]:
            pops.append(keys[i])
    for i in pops: keys.remove(i)
    pops = []
    for i in range(len(fbkeys)):
        if 'string' in fbkeys[i]:
            pops.append(fbkeys[i])
    for i in pops: fbkeys.remove(i)

    recordindices = [[] for i in range(len(allvars))]
    recordtimestamps = [[] for i in range(len(allvars))]

    # output variables (from observations_table)
    ov = []
    for o in keys:
        ov.append([[] for i in range(len(allvars))])
    fb = []
    for o in fbkeys:
        fb.append([[] for i in range(len(allvars))])
    #
    # loading the observed_variables
    #
    obsv = np.asarray(obstab_writetofile[obskeys.index('observed_variable')]) # data.observations_table.observed_variable[:]
    #
    # resorting the data
    #
#     print('resort:start')
    @njit
    def make_vrindex(vridx,ridx,idx):
        l=0
        for i in range(1,len(idx)): # to set the recordindices
            if ridx[i]>ridx[i-1]:
                vridx[ridx[i-1]:ridx[i]]=l # next record after l
                l=i
        vridx[ridx[i]:]=len(idx) # next record for the last element is the len of the data

    ridxall=np.zeros(obsv.shape[0],dtype=np.int64) # reverse index - index of the record index
    j=-1
    for j in range(len(ri)-1):
        ridxall[ri[j]:ri[j+1]]=j
    j+=1
    ridxall[ri[j]:]=j # for the last elemenet
    ridx=[]
    vridx=[]
    absidx=[]
    abscount=0
    for j in range(len(allvars)):
        idx=np.where(obsv==allvars[j])[0] # index of all elements form certain variable j
#         print(j,len(idx),',',end='')
        vridx.append(np.zeros(ri.shape[0],dtype=np.int64)) # all zeros in lenght of record index
        ridx=ridxall[idx] # ridxall where variable is j
        make_vrindex(vridx[-1],ridx,idx)
        vridx[-1]+=abscount # abscount for stacking the recordindex

        absidx.append(copy.copy(idx)) # why copy? - to make sure it's not just the ref. - maybe ok without the cp
        abscount+=len(idx)

    # finishing the sorting 
    #
    absidx=np.concatenate(absidx)
#     try:
#         absidx=np.concatenate(absidx)
#     except:
#         absidx = absidx[0]
    #
    # recordtimestamps are only necessary once
    #
    recordtimestamps = recordtimestamps[0]
    #
    # targetfile has to be a copy of the original file
    #
    print('elapsed converting: ',time.time()-tt)
    tt=time.time()
    if os.path.isfile(targetfile):
        mode='r+'
    else:
        mode='w'
#     print()
#     print('writing '+targetfile)
    
    for i in range(len(keys)):
        ov_vars = np.asarray(obstab_writetofile[i]) # data.observations_table[keys[i]][:]
        ov_vars = ov_vars[absidx]
        if keys[i] == 'index':
            pass
        elif keys[i] == 'observation_id' or keys[i] == 'report_id' or keys[i] == 'sensor_id' or keys[i] == 'source_id':
            alldict = {keys[i]:np.asarray(ov_vars, dtype='S1')}
            write_dict_h5(targetfile, alldict, 'observations_table', {keys[i]: { 'compression': 'gzip' } }, [keys[i]])
        else:
            alldict = pandas.DataFrame({keys[i]:ov_vars})
            write_dict_h5(targetfile, alldict, 'observations_table', {keys[i]: { 'compression': 'gzip' } }, [keys[i]])  

    for i in range(len(fbkeys)):
        fb_vars = np.asarray(fb_writetofile[i]) # data.era5fb[fbkeys[i]][:]
        fb_vars = fb_vars[absidx]
        if fbkeys[i] == 'index' or fbkeys[i] == 'string6' or fbkeys[i] == 'string7' or fbkeys[i] == 'string10':
            pass
        elif fbkeys[i] == 'expver' or fbkeys[i] == 'source@hdr' or fbkeys[i] == 'source_id' or fbkeys[i] == 'statid@hdr':
            alldict = {fbkeys[i]:np.asarray(fb_vars, dtype='S1')}
            write_dict_h5(targetfile, alldict, 'era5fb', {fbkeys[i]: { 'compression': 'gzip' } }, [fbkeys[i]])
        else:
            alldict = pandas.DataFrame({fbkeys[i]:fb_vars})
            write_dict_h5(targetfile, alldict, 'era5fb', {fbkeys[i]: { 'compression': 'gzip' } }, [fbkeys[i]]) 
    #
    # writing the recordindices and recordtimestamp.
    #       
    recordindices=vridx
    for i in range(len(recordindices)):
        testvar = pandas.DataFrame({str(allvars[i]):recordindices[i]})
        write_dict_h5(targetfile, testvar, 'recordindices', {str(allvars[i]): { 'compression': None } }, [str(allvars[i])]) 

    write_dict_h5(targetfile, {'recordtimestamp':rt}, 'recordindices', {'recordtimestamp': { 'compression': None } }, ['recordtimestamp'])

    print('elapsed writing:',time.time()-tt)
    

files = glob.glob('/raid60/scratch/federico/MERGED_DATABASE_OCTOBER2020_sensor/*.nc')
# print(files[:10])

# convert_missing(files[6020])
# convert_missing('/raid60/scratch/federico/MERGED_DATABASE_OCTOBER2020_sensor/0-20000-0-11035_CEUAS_merged_v0.nc')

if __name__ == '__main__':
#    pool = multiprocessing.Pool(processes=20)
#    result_list = pool.map(convert_missing, files[100:1000])
    result_list = list(map(convert_missing, files[0:1000]))
    print(result_list)