   
"""Genereate trend data for multiple stations and save to pickle."""

# import time
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

warnings.filterwarnings("ignore")


def datetime_to_seconds(dates, ref="1900-01-01T00:00:00"):
    """from datetime64 to seconds since 1900-01-01 00:00:00"""
    return ((dates - np.datetime64(ref)) / np.timedelta64(1, "s")).astype(np.int64)


def seconds_to_datetime(seconds, ref="1900-01-01"):
    """from seconds to datetime64"""
    seconds = np.asarray(seconds)
    return pd.to_datetime(seconds, unit="s", origin=ref)


def trend_station(i, dt_from, dt_to):
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
    # multiprocessing here
        with h5py.File(i, "r") as file:
            rts = file["recordindices"]["recordtimestamp"][:]
            idx = np.where(np.logical_and((rts >= dt_from), (rts <= dt_to)))[0]
            plevs = [70000]

            idx_d = {}
            # var_d = {'air_temperature':'126', 'relative_humidty':'138', 'geopotential':'117', 'eastward_wind_speed':'139', 'northward_wind_speed':'140', 'dew_point': '137', 'specific_humidity':'39'}
            var_d = {'dew_point_temperature':'137'}
            # var_d = {"relative_humidty": "138"}

            for j in var_d:
                idx_d[j] = file["recordindices"][var_d[j]][idx]

            masks = {}
            for j in idx_d:
                masks[j] = file["observations_table"]["z_coordinate"][
                    idx_d[j][0]:idx_d[j][-1]
                ]
                masks[j] = np.isin(masks[j], plevs)
                # masks[i] = np.isfinite(masks[i])

            # mask = masks['dew_point_temperature']
            # t_idx = idx_d['dew_point_temperature']
            mask = masks["dew_point_temperature"]
            t_idx = idx_d["dew_point_temperature"]

            df_dict["observation_value"] = list(
                file["observations_table"]["observation_value"][t_idx[0]:t_idx[-1]][
                    mask
                ]
            )
            df_dict["z_coordinate"] = list(
                file["observations_table"]["z_coordinate"][t_idx[0]:t_idx[-1]][mask]
            )
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
            df_dict["humidity_bias_estimate"] = list(
                file["advanced_homogenisation"]["humidity_bias_estimate"][
                    t_idx[0]:t_idx[-1]
                ][mask]
            )
            # print('reading data: ', time.time()-t0, ' s')
            # t1 = time.time()
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
            # print('to dataframe: ', time.time()-t1, ' s')
            # t2 = time.time()
            # display(temp)

            if len(temp) > 0:
                # temp.sort_values('time')
                temp["time"] = pd.to_datetime(temp["date_time"])
                temp["lat"] = np.array([temp.latitude.iloc[-1]] * len(temp))
                temp["lon"] = np.array([temp.longitude.iloc[-1]] * len(temp))
                temp["adjusted"] = (
                    temp["observation_value"] - temp["humidity_bias_estimate"]
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
                    sout.append(float(out.iloc[-1] * 3650))
                    # print('trend ua: ', time.time() - t3, ' s')
                    # t4 = time.time()
                    out_adj = rasotools.met.time.trend(
                        xa.adjusted, only_slopes=True, method="polyfit"
                    ).to_dataframe(
                        name="out_adj"
                    )  # leave methode empty for default -> robust Theil–Sen estimator for trend
                    aout.append(float(out_adj.iloc[-1] * 3650))
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


# files =  glob.glob("/mnt/users/scratch/leo/scratch/converted_v19/long/*11035*.nc")
# out = trend_station(files[0], datetime_to_seconds(np.datetime64("1994-01-01")), datetime_to_seconds(np.datetime64("2023-12-31")))
# print(out)

files = glob.glob("/mnt/users/scratch/leo/scratch/converted_v19/long/*.nc")
# files =  glob.glob('/users/staff/uvoggenberger/scratch/humtest/*.nc')


# dt_from = datetime_to_seconds(np.datetime64("1973-01-01"))
# dt_to = datetime_to_seconds(np.datetime64("2002-12-31"))
# dt_from = datetime_to_seconds(np.datetime64("1958-01-01"))
# dt_to = datetime_to_seconds(np.datetime64("1987-12-31"))
for dts, dte in [("1994-01-01", "2023-12-31"), ("1973-01-01", "2002-12-31")]: # , ("1958-01-01", "1987-12-31")
    
    print(dts, dte)
    print('----')
    
    dt_from = datetime_to_seconds(np.datetime64(dts))
    dt_to = datetime_to_seconds(np.datetime64(dte))

    test_r = ray.remote(trend_station)
    ray.init(num_cpus=20)

    # result_ids = []
    # for i in files:
    #     result_ids.append(test_r.remote(i, dt_from, dt_to))
    # results = ray.get(result_ids)

    results = []
    obj_ids = [test_r.remote(i,dt_from,dt_to) for i in files]
    for x in tqdm(to_iterator(obj_ids), total=len(obj_ids)):
        results.append(x)

    # pickle.dump(results, open("polyfit_trends_dewpoint_700hPa_1958_1988_Trend_20230717.p", "wb"))
    # pickle.dump(results, open("polyfit_trends_dewpoint_700hPa_1973_2003_Trend_20230717.p", "wb"))
    pickle.dump(results, open("/users/staff/uvoggenberger/scratch/CUON_trends/polyfit_trends_dewpoint_700hPa_"+dts+"_"+dte+"_Trend_20240422.p", "wb"))
    ray.shutdown()
