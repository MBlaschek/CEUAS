#!/usr/bin/env
# coding: utf-8

import h5py
import numpy as np
import pandas as pd
import os, glob, sys
import math

sys.path.append(os.getcwd()+'/../cds-backend/code/')
os.environ['PYTHONPATH'] = os.getcwd()+'/../cds-backend/code/'
import cds_eua4 as eua
import trajectory as trj

import ray
ray.init(num_cpus=20)



@ray.remote
def write_trj(stat):
#     stat = '11035'
#     try:
    # check if output already exists:
    checkfile = glob.glob('/mnt/users/staff/a1400070/scratch/trajectory_files_new/*'+str(stat)+'*.nc')
    if len(checkfile) > 0:
        return 1

    # checking if it was written to file already
#         with h5py.File(checkfile[0], 'r') as inputfile:
#             if ( 'displaced_time' in inputfile['advanced_homogenisation']):
#                 return 1

#             if not ('advanced_homogenisation' in inputfile):
#                 inputfile.create_group('advanced_homogenisation')

    # read converted file and open it
    statlist = glob.glob('/mnt/users/scratch/leo/scratch/converted_v8/*' + stat + '*_CEUAS_merged_v1.nc')
    file = eua.CDMDataset(filename = statlist[0])

    # check if all variables are available:
    if (not ('126' in file.recordindices.keys())) or (not ('139' in file.recordindices.keys())) or (not ('140' in file.recordindices.keys())):
        return 2

#     try:
#         igra_file = glob.glob('/scratch/das/federico/COP2_HARVEST_APRIL2022/igra2/*'+stat+'*')[0]
#         i_file = eua.CDMDataset(filename = igra_file)
#         igra_file_avail = True
#     except:
#         igra_file_avail = False
#     print('igra_file_avail:', igra_file_avail)
    
    try:
#         y = file.observations_table.z_coordinate_type[:]
#         x = file.observations_table.z_coordinate[:]
#         x = x[y==1]
#         x = x[~np.isnan(x)]
#         x = (np.sort(np.unique(x))[-100:])
#         plev_threshold = np.median(x) - 1000


        statlen = len(file.observations_table.observed_variable[:])
        latd = np.full(statlen, np.nan)
        lond = np.full(statlen, np.nan)
        timed = np.full(statlen, np.nan)
        ttime = np.full(statlen, np.nan)

        slat = file.observations_table.latitude[0]
        slon = file.observations_table.longitude[0]
    except:
        return 4
    print('ascents: ', len(file.header_table.report_id[:])-1)
    for i in range(len(file.header_table.report_id[:])-1):
        var_recidx = {}
        try:
            for j in file.recordindices.keys():
                if j not in ['index', 'recordtimestamp']:
                    var_recidx[j]=[file.recordindices[j][i], file.recordindices[j][i+1]]
        except:
            return 5
        
        t_idx_s = var_recidx['126'][0] #file.recordindices['126'][i]
        t_idx_e = var_recidx['126'][1] #file.recordindices['126'][i+1]
        u_idx_s = var_recidx['139'][0] #file.recordindices['139'][i]
        u_idx_e = var_recidx['139'][1] #file.recordindices['139'][i+1]
        v_idx_s = var_recidx['140'][0] #file.recordindices['140'][i]
        v_idx_e = var_recidx['140'][1] #file.recordindices['140'][i+1]

        # -----------

        if (u_idx_s == u_idx_e) or (u_idx_s == u_idx_e) or (t_idx_s == t_idx_e):
            # replace with nan filling
            continue
        
        # check for z_coordinate_type
        z_coordinate_type = file.observations_table.z_coordinate_type[t_idx_s:t_idx_e]
        if len(np.where(z_coordinate_type == 1)[0]) < 3:
            continue

        # -----------

        repid = file.header_table.report_id[i]
        u = file.observations_table.observation_value[u_idx_s:u_idx_e]
        v = file.observations_table.observation_value[v_idx_s:v_idx_e]
        t = file.observations_table.observation_value[t_idx_s:t_idx_e]
        z_coordinate_t = file.observations_table.z_coordinate[t_idx_s:t_idx_e]
        z_coordinate_u = file.observations_table.z_coordinate[u_idx_s:u_idx_e]
        z_coordinate_v = file.observations_table.z_coordinate[v_idx_s:v_idx_e]

        # find shortest array
        z_coords = [z_coordinate_t, z_coordinate_u, z_coordinate_v]
        shortest_zc = 0

        for k in range(len(z_coords)-1):
            if len(z_coords[k+1]) < len(z_coords[shortest_zc]):
                shortest_zc = k+1

        # if shortest array < 7 -> skip
        if len(z_coords[shortest_zc]) < 3:
            # replace with nan filling
            continue

        # -----------

        # check if it's neccessary to interpolate
        check = False
        if (len(np.array(z_coordinate_t)) == len(np.array(z_coordinate_u))):
            if not (np.array(z_coordinate_t) == np.array(z_coordinate_u)).all():
                check = True
        else:
            check = True

        if check:
            u_new = []
            v_new = []

            for k in range(len(z_coordinate_t)):
                if z_coordinate_u[0] < z_coordinate_t[k] < z_coordinate_u[-1]:
                    u_new.append(np.interp(z_coordinate_t[k], z_coordinate_u, u))
                    v_new.append(np.interp(z_coordinate_t[k], z_coordinate_v, v))
                else:
                    u_new.append(np.nan)
                    v_new.append(np.nan)

            u = np.array(u_new)
            v = np.array(v_new)

        # -----------
        input_df = pd.DataFrame({'t':t, 'u':u, 'v':v, 'p':z_coordinate_t, 'idx':np.array(range(t_idx_s, t_idx_e))})
        # flip for ascending order
        input_df = input_df.dropna().iloc[::-1].reset_index()
        if len(input_df) < 3:
            continue
#         print(input_df)



        # check if best possible time is selected
        date_time = file.observations_table.date_time[t_idx_s] #dt_date = pd.to_datetime(date_time, unit='s', origin='1900-01-01')
        '''
        if not (file.era5fb.reportype[t_idx_s] == 16045):
            if igra_file_avail:
                #if not already igra data -> select igra datetime
                if int(repid[0]) != 3:
                    dups = file.header_table.duplicates[i]
                    dups = dups[dups != b'']
                    dups = dups[dups != b',']
                    #iterate through all duplicates
                    for j in range(0,int((len(dups)/11))):
                        #if there is an igra duplicate:
                        if int((dups[(j*11):((j+1)*11)])[0]) == 3:
                            save_id = 0
                            deci = 1
                            a = ((dups[(j*11):((j+1)*11)]))[1:]
                            for o in np.flip(a):
                                save_id += int(o)* deci
                                deci = deci*10
                            date_time = i_file.recordtimestamp[save_id]  
        #                     dt_date = pd.to_datetime(date_time, unit='s', origin='1900-01-01')

        #                     print(i)
        #                     print((dups[(j*11):((j+1)*11)]))
        #                     print(save_id)
        #                     print('igra_time', pd.to_datetime(i_file.recordtimestamp[save_id], unit='s', origin='1900-01-01'))
        #                     print('cuon time', dt_date)
        #                     print()


            #     print('t: ', len(t), 'u: ', len(u), 'v: ', len(v), )

        # ----------------
        '''
        phys_model = trj.trajectory(lat=slat, lon=slon, temperature=np.array(input_df.t), u=np.array(input_df.u), v=np.array(input_df.v), pressure=np.array(input_df.p))
        ###
        #
        # check for phys_model[4] -> do not write to output if input_df.p[0] < 850
        # check for station height or mean hPa of first level and compare to that? or min(p)?
        #
        ###
    #     input_df['latitude_displacement'] = phys_model[0]
    #     input_df['longitude_displacement'] = phys_model[1]
    #     input_df['displaced_time'] = np.array(phys_model[4]) + date_time
    #     input_df['displaced_time_dt'] = pd.to_datetime(input_df['displaced_time'], unit='s', origin='1900-01-01')
    #     display(input_df)
    
        helper = list(input_df.idx)
#         print(helper, np.array(phys_model[0]))
        latd[helper] = np.array(phys_model[0])
        lond[helper] = np.array(phys_model[1])
        timed[helper] = np.array(phys_model[4])
        ttime[helper] = np.array(phys_model[4])+date_time
        
#         for j in var_recidx:
#             helper = np.array(range(var_recidx[j][0], var_recidx[j][1]))
#             print('helper', helper)
#             if len(helper) == len(phys_model[0]):
#                 latd[helper] = np.array(phys_model[0])
#                 lond[helper] = np.array(phys_model[1])
# #                 if input_df.p[0] > plev_threshold:
#                 timed[helper] = np.array(phys_model[4])
#                 ttime[helper] = np.array(phys_model[4])+date_time

    file.close()

#         output_file = glob.glob('/mnt/users/staff/a1400070/scratch/converted_v8/*'+stat+'*.nc')[0]
    output_file = '/mnt/users/staff/a1400070/scratch/trajectory_files_new/trajectory_'+str(stat)+'.nc'
#         with h5py.File(output_file, 'r+') as newfile:
    with h5py.File(output_file, 'w') as newfile:
        newfile.create_group('advanced_homogenisation')
        newfile['advanced_homogenisation'].create_dataset('latitude_displacement', data=latd)
        newfile['advanced_homogenisation']['latitude_displacement'].attrs['version'] = b'1.0'
        newfile['advanced_homogenisation'].create_dataset('longitude_displacement', data=lond)
        newfile['advanced_homogenisation']['longitude_displacement'].attrs['version'] = b'1.0'
        newfile['advanced_homogenisation'].create_dataset('time_since_launch', data=timed)
        newfile['advanced_homogenisation']['time_since_launch'].attrs['version'] = b'1.0'
        newfile['advanced_homogenisation'].create_dataset('true_time', data=ttime)
        newfile['advanced_homogenisation']['true_time'].attrs['version'] = b'1.0'

    return 0 
#     except Exception as e:
#         print(e)
#         return 3

        
if __name__ == '__main__':
    
    file_list = []
    result_ids = []
    for i in glob.glob('/mnt/users/scratch/leo/scratch/converted_v8/*.nc')[:]:
        file_list.append(i.split('_CEUAS_merged_v1.nc')[0][-5:])
#         print(write_trj(file_list[-1]))
        result_ids.append(write_trj.remote(file_list[-1]))
                    
    results = ray.get(result_ids)
    print(results)
