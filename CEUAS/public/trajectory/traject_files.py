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
from harvest_convert_to_netCDF import write_dict_h5

import ray
ray.init(num_cpus=30)



@ray.remote
def write_trj(stat):
    test_counter = 0
    try:
        # check if output already exists:
        targetfile = '/mnt/users/staff/uvoggenberger/scratch/converted_v11/trajectory_files_new/trajectory_'+str(stat.split('/')[-1])
        checkfile = glob.glob(targetfile)
        # if target file already exists
        if len(checkfile) > 0:
            # if input is older than target
            if os.path.getmtime(stat) < os.path.getmtime(targetfile):
                return 1
        
        

        # read converted file and open it
        file = eua.CDMDataset(filename = stat)

        # check if all variables are available:
        # if (not ('126' in file.recordindices.keys())) or (not ('139' in file.recordindices.keys())) or (not ('140' in file.recordindices.keys())):
            # return 2

    #     try:
    #         igra_file = glob.glob('/scratch/das/federico/COP2_HARVEST_APRIL2022/igra2/*'+stat+'*')[0]
    #         i_file = eua.CDMDataset(filename = igra_file)
    #         igra_file_avail = True
    #     except:
    #         igra_file_avail = False
    #     print('igra_file_avail:', igra_file_avail)

        # creating fillable output variables
        try:
            statlen = len(file.observations_table.observed_variable[:])
            latd = np.full(statlen, np.nan)
            lond = np.full(statlen, np.nan)
            timed = np.full(statlen, np.nan)
            ttime = np.full(statlen, np.nan)

            slat = file.observations_table.latitude[0]
            slon = file.observations_table.longitude[0]
        except:
            return 4
        if ('126' in file.recordindices.keys()) and ('139' in file.recordindices.keys()) and ('140' in file.recordindices.keys()):
            print('ascents: ', len(file.header_table.report_id[:])-1)
            for i in range(len(file.header_table.report_id[:])-1):
                var_recidx = {}
                try:
                    for j in file.recordindices.keys():
                        if j not in ['index', 'recordtimestamp']:
                            var_recidx[j]=[file.recordindices[j][i], file.recordindices[j][i+1]]
                except:
                    return 5

                t_idx_s = var_recidx['126'][0]
                t_idx_e = var_recidx['126'][1]
                u_idx_s = var_recidx['139'][0]
                u_idx_e = var_recidx['139'][1]
                v_idx_s = var_recidx['140'][0]
                v_idx_e = var_recidx['140'][1]

                # -----------

                if (u_idx_s == u_idx_e) or (u_idx_s == u_idx_e) or (t_idx_s == t_idx_e):
                    # replace with nan filling
                    continue

                # check for z_coordinate_type
                z_coordinate_type = file.observations_table.z_coordinate_type[t_idx_s:t_idx_e]
                if len(np.where(z_coordinate_type == 1)[0]) < 3:
                    continue

                # create variables needed for dataframe
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

                # if shortest array < 3 -> skip
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
                        if z_coordinate_u[0] <= z_coordinate_t[k] <= z_coordinate_u[-1]:
                            u_new.append(np.interp(z_coordinate_t[k], z_coordinate_u, u))
                            v_new.append(np.interp(z_coordinate_t[k], z_coordinate_v, v))
                        else:
                            u_new.append(np.nan)
                            v_new.append(np.nan)

                    u = np.array(u_new)
                    v = np.array(v_new)

                # done collecting data for calculation
                input_df = pd.DataFrame({'t':t, 'u':u, 'v':v, 'p':z_coordinate_t, 'idx':np.array(range(t_idx_s, t_idx_e))})

                # clean input array and check for outliers:
                input_df.drop(input_df[input_df.t < 172].index, inplace=True)
                input_df.drop(input_df[input_df.t > 372].index, inplace=True)
                input_df.drop(input_df[input_df.u > 150].index, inplace=True)
                input_df.drop(input_df[input_df.v > 150].index, inplace=True)

                # flip for ascending order
                input_df = input_df.dropna().iloc[::-1].reset_index()

                # last check if not to little levels # 3 levels is minimum -> maybe more should be needed
                # if changed - also change for shortest array!
                if len(input_df) < 3:
                    continue

                # check for lowest level pressure - is it too high above the ground?
                # whats the lowest levels pressure in Pa?
                p_lowest_level = input_df.p.iloc[0] # Pa
                t_lowest_level = input_df.t.iloc[0] # K
                # what is the acepted pressure range around that, for a given station height?
                # delta_p = g/(R*T) * p_0 * delta_h
                # delta_h = delta_p * R * T / (g * p_0)
                R = 287 # J/(kg K)
                g = 9.80 # m/s²

                # low_msl_pressure = 99000 # Pa
                # mean_msl_pressure = 101300 # Pa
                high_msl_pressure = 103000 # Pa

                # z_low = (low_msl_pressure - p_lowest_level) * R * t_lowest_level / g / low_msl_pressure
                # z_mean = (mean_msl_pressure - p_lowest_level) * R * t_lowest_level / g / mean_msl_pressure
                z_high = (high_msl_pressure - p_lowest_level) * R * t_lowest_level / g / high_msl_pressure

                # print('maximal height above ground', z_high, p_lowest_level)

                station_z = file.observations_table.station_elevation[input_df.idx.iloc[0]]
                # if the first observation is more than 1500 m above ground, the ascent is invalid (we check for high and low pressure)
                if np.logical_and(((z_high - station_z) > 1500), ((z_high - station_z) > 1500)):
                    continue


                # check if best possible time is selected # skipped for now
                date_time = file.observations_table.date_time[t_idx_s] 

                # if dt format is needed convert with:
                # dt_date = pd.to_datetime(date_time, unit='s', origin='1900-01-01')
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
                '''

                # calculate trajectory
                phys_model = trj.trajectory(lat=slat, lon=slon, temperature=np.array(input_df.t), u=np.array(input_df.u), v=np.array(input_df.v), pressure=np.array(input_df.p))

                # helper knows where to write the data
                helper = list(input_df.idx)

                # filling output variables with calculated data
                latd[helper] = np.array(phys_model[0])
                lond[helper] = np.array(phys_model[1])
                timed[helper] = np.array(phys_model[4])
                ttime[helper] = np.array(phys_model[4])+date_time

                # # stopper for some tests:
                # if test_counter > 100:
                #     break
                # else:
                #     test_counter += 1


        file.close()


        # writing to homogenisations only file

        
        mode='w'
        group = 'advanced_homogenisation'

        i = 'latitude_displacement'
        ov_vars = latd
        alldict = pd.DataFrame({i:ov_vars})
        write_dict_h5(targetfile, alldict, group, {i: { 'compression': 'gzip' } }, [i]) 

        i = 'longitude_displacement'
        ov_vars = lond
        alldict = pd.DataFrame({i:ov_vars})
        write_dict_h5(targetfile, alldict, group, {i: { 'compression': 'gzip' } }, [i]) 

        i = 'time_since_launch'
        ov_vars = timed
        alldict = pd.DataFrame({i:ov_vars})
        write_dict_h5(targetfile, alldict, group, {i: { 'compression': 'gzip' } }, [i]) 

        i = 'true_time'
        ov_vars = ttime
        alldict = pd.DataFrame({i:ov_vars})
        write_dict_h5(targetfile, alldict, group, {i: { 'compression': 'gzip' } }, [i]) 


        # writing to input file

    #     targetfile = '/mnt/users/staff/uvoggenberger/scratch/hum_adj_2022/'+str(stat.split('/')[-1])
    #     mode='r+'
    #     group = 'advanced_homogenisation'

    #     i = 'latitude_displacement'
    #     ov_vars = latd
    #     alldict = pd.DataFrame({i:ov_vars})
    #     write_dict_h5(targetfile, alldict, group, {i: { 'compression': 'gzip' } }, [i]) 

    #     i = 'longitude_displacement'
    #     ov_vars = lond
    #     alldict = pd.DataFrame({i:ov_vars})
    #     write_dict_h5(targetfile, alldict, group, {i: { 'compression': 'gzip' } }, [i]) 

    #     i = 'time_since_launch'
    #     ov_vars = timed
    #     alldict = pd.DataFrame({i:ov_vars})
    #     write_dict_h5(targetfile, alldict, group, {i: { 'compression': 'gzip' } }, [i]) 

    #     i = 'true_time'
    #     ov_vars = ttime
    #     alldict = pd.DataFrame({i:ov_vars})
    #     write_dict_h5(targetfile, alldict, group, {i: { 'compression': 'gzip' } }, [i]) 



        return 0 
    except Exception as e:
        print(e)
        return 3

        
if __name__ == '__main__':
    
    file_list = []
    result_ids = []
    for i in glob.glob('/mnt/users/scratch/leo/scratch/converted_v11/long/*.nc')[:]:
        # file_list.append(i.split('_CEUAS_merged_v1.nc')[0][-5:])
        file_list.append(i)
        result_ids.append(write_trj.remote(file_list[-1]))
                    
    results = ray.get(result_ids)
    print(results)
