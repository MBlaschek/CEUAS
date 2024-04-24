#!/usr/bin/env python

import os
import sys
import glob
import warnings

import numpy as np
import pickle
import pandas as pd
import h5py
import ray
import hdf5plugin
from tqdm import tqdm
from datetime import datetime
import argparse
import json

sys.path.append(os.path.expanduser("~")+'/uvpy/')
import uvfunctions as uvf
import uvplot as uvp
import uvtests as uvt


warnings.filterwarnings("ignore")

def show_trend_map(file,label, c_bar, multiplier = 1, c_bar_red_top = True, cbar_limit = 2, save_at=''):
    sys.path.append(os.getcwd()+'/../resort/rasotools-master/')
    import rasotools
    import hdf5plugin
    
    plev = file.split('_')[-7]
    sdate = file.split('_')[-4]
    edate = file.split('_')[-3]
    good_results = {}
    good_results['lat'] = []
    good_results['lon'] = []
    good_results['st'] = []
    good_results['at'] = []
    good_results['label'] = []
    
    results = pickle.load(open(file, 'rb'))
    for i in results:
        badflag = 0
        for j in i:
            try:
                if np.isnan(j):
                    badflag = 1
            except:
                pass
        if badflag == 0:
            good_results['lat'].append(i[1])
            good_results['lon'].append(i[2])
            good_results['st'].append(i[3][0])
            good_results['at'].append(i[4][0])
            good_results['label'].append(i[0].split('/')[-1].split('_')[0])
    da = pd.DataFrame.from_dict(good_results)
    statnum = len(da)

    a = rasotools.plot._helpers.cost(np.asarray(da.lon), np.asarray(da.lat), np.asarray(da.st))
    cost = np.sum(a)/len(a)
    ua_save_at = save_at[:-4]  + '_unadjusted.png'
    print(ua_save_at)
    fig1 = uvp.world_map_mpl(da.lat, da.lon, da['st'], label + ' trend unadjusted \n '+str(plev)+'_'+sdate+'_'+edate+' \n heterogeneity cost: ' + str(cost), inp_vmin=-cbar_limit, inp_vmax=cbar_limit, invert_cbar=c_bar_red_top, cbar_label = c_bar,save_at = ua_save_at)

    a = rasotools.plot._helpers.cost(np.asarray(da.lon), np.asarray(da.lat), np.asarray(da['at']))
    cost = np.sum(a)/len(a)
    a_save_at = save_at[:-4] + '_adjusted.png'
    print(a_save_at)
    fig1 = uvp.world_map_mpl(da.lat, da.lon, da['at'], label + ' trend adjusted \n '+str(plev)+'_'+sdate+'_'+edate+' \n heterogeneity cost: ' + str(cost), inp_vmin=-cbar_limit, inp_vmax=cbar_limit, invert_cbar=c_bar_red_top, cbar_label = c_bar, save_at = a_save_at)



def datetime_to_seconds(dates, ref="1900-01-01T00:00:00"):
    """from datetime64 to seconds since 1900-01-01 00:00:00"""
    return ((dates - np.datetime64(ref)) / np.timedelta64(1, "s")).astype(np.int64)


def seconds_to_datetime(seconds, ref="1900-01-01"):
    """from seconds to datetime64"""
    seconds = np.asarray(seconds)
    return pd.to_datetime(seconds, unit="s", origin=ref)


def trend_station(i, dt_from, dt_to, variable, adjustment, pressure):
    """i ... station file path"""
    # t0 = time.time()

    sys.path.append(os.getcwd() + "/../resort/rasotools-master/")
    import rasotools
    import hdf5plugin

    df_dict = {}
    sout = []
    aout = []
    stats = []
    lats = []
    lons = []
    try:
        with h5py.File(i, "r") as file:
            rts = file["recordindices"]["recordtimestamp"][:]
            idx = np.where(np.logical_and((rts >= dt_from), (rts <= dt_to)))[0]
            plevs = [pressure]

            idx_d = {}
            var_d = {'ta':'126', 'rh':'138', 'u':'139', 'v':'140', 'dp': '137', 'sh':'39', 'wd':'106'} # gp:117

            idx_d[variable] = file["recordindices"][var_d[variable]][idx]

            masks = {}
            for j in idx_d:
                masks[j] = file["observations_table"]["z_coordinate"][
                    idx_d[j][0]:idx_d[j][-1]
                ]
                masks[j] = np.isin(masks[j], plevs)
                # masks[i] = np.isfinite(masks[i])

            mask = masks[variable]
            t_idx = idx_d[variable]

            df_dict["observation_value"] = list(file["observations_table"]["observation_value"][t_idx[0]:t_idx[-1]][mask])
            df_dict["z_coordinate"] = list(file["observations_table"]["z_coordinate"][t_idx[0]:t_idx[-1]][mask])
            df_dict["date_time"] = seconds_to_datetime(
                list(
                    file["observations_table"]["date_time"][t_idx[0]:t_idx[-1]][mask]
                )
            )
            df_dict["latitude"] = list(
                file["observations_table"]["latitude"][t_idx[0]:t_idx[-1]][mask]
            )
            df_dict["longitude"] = list(
                file["observations_table"]["longitude"][t_idx[0]:t_idx[-1]][mask]
            )
            df_dict["temperature_bias_estimate"] = list(
                file["advanced_homogenisation"][adjustment][
                    t_idx[0]:t_idx[-1]
                ][mask]
            )
            temp = pd.DataFrame.from_dict(df_dict)

            # drop the 1st and 99th percentile
            q_01 = temp["observation_value"].quantile(
                0.01
            )  # get the 1st percentile value
            q_99 = temp["observation_value"].quantile(
                0.99
            )  # get the 99th percentile value
            temp = temp[
                (temp["observation_value"] > q_01) & (temp["observation_value"] < q_99)
            ]

            if len(temp) > 0:
                # temp.sort_values('time')
                temp["time"] = pd.to_datetime(temp["date_time"])
                temp["lat"] = np.array([temp.latitude.iloc[-1]] * len(temp))
                temp["lon"] = np.array([temp.longitude.iloc[-1]] * len(temp))
                temp["adjusted"] = (
                    temp["observation_value"] - temp["temperature_bias_estimate"]
                )
                temptime = temp.time
                if len(temp) >= 19 * 365 and len(np.unique(temptime.dt.year)) > 19:
                    # print('enough data')
                    xa = temp.set_index(["lat", "lon", "time"]).to_xarray()
                    # and do it twice for the adjusted values too!
                    # print('data prep: ', time.time()-t2, ' s')
                    # t3 = time.time()
                    out = rasotools.met.time.trend(
                        xa.observation_value, only_slopes=True, method="polyfit"
                    ).to_dataframe(
                        name="out"
                    )  # leave methode empty for default -> robust Theil–Sen estimator for trend
                    sout.append(float(out.out.iloc[-1]) * 3650)
                    # print('trend ua: ', time.time() - t3, ' s')
                    # t4 = time.time()
                    out_adj = rasotools.met.time.trend(
                        xa.adjusted, only_slopes=True, method="polyfit"
                    ).to_dataframe(
                        name="out_adj"
                    )  # leave methode empty for default -> robust Theil–Sen estimator for trend
                    aout.append(float(out_adj.out_adj.iloc[-1]) * 3650)
                    # print('trend adj: ', time.time() - t4, ' s')
                    # t5 = time.time()
                else:
                    sout = np.nan
                    aout = np.nan
            else:
                sout = np.nan
                aout = np.nan
            try:
                lats = temp.latitude.iloc[-1]
                lons = temp.longitude.iloc[-1]
                stats = i
            except:
                lats = np.nan
                lons = np.nan
                stats = i
    except:
        lats = np.nan
        lons = np.nan
        stats = i
        sout = np.nan
        aout = np.nan

    # print('over all: ', time.time()-t0, ' s')
    # print(sout, aout)
    return [stats, lats, lons, sout, aout]


def to_iterator(obj_ids):
    while obj_ids:
        done, obj_ids = ray.wait(obj_ids)
        yield ray.get(done[0])

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Create trend plots form CUON data")
    parser.add_argument('--variable' , '-v', 
                        help="Select the variable to calculate trend of. Available options: ta (temperature), dp (dewpoint), rh (relative humidity), u (u component of wind), v (v component of wind), wd (wind direction). If not selected, the script will run the ta."  ,
                        default = 'ta',
                        type = str)

    parser.add_argument('--startdate' , '-sd', 
                        help='Trend start date. Provide string in following format: "1994-01-01""' ,
                        default = '1958-01-01',
                        type = str)    

    parser.add_argument('--enddate' , '-ed',
                        help = 'Trend end date. Provide string in following format: "2023-12-31"', 
                        default = '1987-12-31',
                        type = str)    

    parser.add_argument('--input' , '-i',
                        help = "Path to the input directory, where *.nc files are."  ,
                        default = "/mnt/users/scratch/leo/scratch/converted_v19/long/",
                        type = str)    

    parser.add_argument('--output' , '-o',
                        help = "Path to the output directory, where notebooks are created."  ,
                        default = './cuon_trends',
                        type = str)   

    parser.add_argument('--adjustment' , '-a',
                        help = "Select adjustment for trend calculation. Available options: default (for default adjustment), RAOBCORE_bias_estimate, RASE_bias_estimate, RICH_bias_estimate, RISE_bias_estimate (for ta)"  ,
                        default = 'default',
                        type = str) 
    
    parser.add_argument('--threads' , '-nt',
                        help = "Thread number for multiprocessing. Default = 20"  ,
                        default = 40,
                        type = int)  
    
    parser.add_argument('--pressure' , '-p',
                        help = "Select the pressure in Pa on which the trends should be calculated. Default = 10000"  ,
                        default = 10000,
                        type = int)

    parser.add_argument('--cbar_limit' , '-cbl',
                        help = "Select the integer limit for the colorbar. Default = 2"  ,
                        default = 2,
                        type = float)
    
    parser.add_argument('--force_rerun' , '-fr',
                        help = "Force a rerun, even if data with same timecode exists. Set to 1 for forced rerun. Default = 0"  ,
                        default = 0,
                        type = int)
    
    parser.add_argument('--version' , '-vr',
                        help = "Version to keep track of trends. Default = 0"  ,
                        default = "V0",
                        type = str)
    
    args = parser.parse_args()
    input = args.input 
    out_dir = args.output
    variable = args.variable
    start_date = args.startdate
    end_date = args.enddate
    threads = args.threads
    adjustment = args.adjustment
    pressure = args.pressure
    cbl = args.cbar_limit
    fr = args.force_rerun
    version = args.version


    ## Input files:
    files = glob.glob(input + "*.nc")

    ## Creating output directory:
    if not os.path.isdir(out_dir):
        os.system('mkdir ' + out_dir )
       
    ## select adjustment:humidity_bias_estimate
    default_adjustment = {'ta':'RAOBCORE_bias_estimate', 'rh':'humidity_bias_estimate', 'u':'wind_bias_estimate', 'v':'wind_bias_estimate', 'dp': 'humidity_bias_estimate', 'sh':'humidity_bias_estimate', 'wd':'wind_bias_estimate'}
    if adjustment == 'default':
        adjustment = default_adjustment[variable]

    ## add identifier for creation:
    write_data_to = out_dir
    if write_data_to[-1] != '/':
        write_data_to += '/'
    for add_path in ['polyfit_trends', '_', variable, '_', str(pressure), '_', adjustment, '_', start_date, '_', end_date, '_', 'trendfile', '_', version, '.p']:
        write_data_to += add_path

    ## run trend calculation if file not already there
    print(start_date, end_date)
    print(variable, adjustment)
    print(pressure, out_dir)
    if not os.path.isfile(write_data_to) or fr == 1:
        print('creating files')

        dt_from = datetime_to_seconds(np.datetime64(start_date))
        dt_to = datetime_to_seconds(np.datetime64(end_date))


        test_r = ray.remote(trend_station)
        ray.init(num_cpus=threads)
        results = []
        obj_ids = [test_r.remote(i,dt_from, dt_to, variable, adjustment, pressure) for i in files]
        for x in tqdm(to_iterator(obj_ids), total=len(obj_ids)):
            results.append(x)
        
        pickle.dump(results, open(write_data_to, "wb"))
        ray.shutdown()

        
    else:
        print('trend file already exists')
    print('----')
    units = {'ta':'K', 'rh':'1', 'u':'m/s', 'v':'m/s', 'dp': 'K', 'sh':'g/kg', 'wd':'°'}
    cbar_red_equals_negative = {'ta':False, 'rh':True, 'u':False, 'v':False, 'dp': True, 'sh':True, 'wd':False}

    trend_notebook_info = {}
    trend_notebook_info['file'] = write_data_to
    trend_notebook_info['variable'] = variable
    trend_notebook_info['label'] = variable + ' ' + str(pressure)
    trend_notebook_info['c_bar'] = '['+units[variable]+']/10a'
    trend_notebook_info['c_bar_red_top'] = cbar_red_equals_negative[variable]
    

    # with open("trend_notebook_info.json", "w") as write_file:
    #     json.dump(trend_notebook_info, write_file)

    save_at = write_data_to[:-2] + '.png'
    show_trend_map(trend_notebook_info['file'], trend_notebook_info['label'], trend_notebook_info['c_bar'], c_bar_red_top=trend_notebook_info['c_bar_red_top'], cbar_limit = cbl, save_at=save_at)

