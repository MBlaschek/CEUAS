import h5py
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import glob


# def z_to_p_ifs(h): # geopotential height (m^/s^2) to pressure (Pa)
#     a = 5.252368255329
#     b = 44330.769230769
#     c = 0.000157583169442
#     ptro = 226.547172
#     po = 1013.25
#     g = 9.80665
 
#     h /= g
#     if h != h:

#         p = h

#     elif h > 11000.0:

#         y = -c * (h - 11000.0)

#         p = ptro * np.exp(y)

#     else:

#         y = 1.0 - h / b

#         p = po * (y**a)

#     return p * 100. # we want Pa


# # @ray.remote
# def process_station_file(file):

#     with h5py.File(file, 'r') as fl:
#         x = fl['observations_table']['observed_variable'][:]
#         val = fl['observations_table']['observation_value'][:]
#         p = fl['observations_table']['z_coordinate'][:]
#         rid_raw = fl['observations_table']['report_id'][:]
#         lat = fl['observations_table']['latitude'][:]
#         lon = fl['observations_table']['longitude'][:]

#     # Convert report_id to strings
#     report_ids = ["".join(byte.decode() for byte in rid).strip() for rid in rid_raw]

#     # Group by report_id
#     from collections import defaultdict
#     report_obs = defaultdict(list)

#     for var, v, pr, rid, la, lo in zip(x, val, p, report_ids, lat, lon):
#         report_obs[rid].append((var, v, pr, la, lo))

#     results_temp = []
#     results_pilot = []
#     results_pseudo_pilot = []

#     for rid, obs in report_obs.items():
#         observed_vars = [o[0] for o in obs]
#         has_126 = 126 in observed_vars
        
#         # Open a wider window for pressure conversion:
#         obs_100_500 = [o for o in obs if 10000 <= o[2] <= 50000]

#         u_vals = [o[1] for o in obs_100_500 if o[0] == 139]
#         v_vals = [o[1] for o in obs_100_500 if o[0] == 140]
#         g_vals = [o[1] for o in obs_100_500 if o[0] == 117]

#         # Ensure all data is of the same shape:
#         if len(g_vals) != len(u_vals) or len(g_vals) == 0:
#             continue
        
#         # Convert geopotential height to pressure using IFS formula
#         p_vals = [z_to_p_ifs(i_g) for i_g in g_vals]

#         # Recreate obs with only u_vals, v_vals, lat, lon, and use p_vals as pressure
#         obs_uv = []
#         i_g_v = 0
#         i_g_u = 0
#         for idx, (var, v, pr, la, lo) in enumerate(obs_100_500):
#             if var in [139]:
#                 obs_uv.append((var, v, p_vals[i_g_u], la, lo))
#                 i_g_u += 1
#             if var in [140]:
#                 obs_uv.append((var, v, p_vals[i_g_v], la, lo))
#                 i_g_v += 1
        
#         # Filter pressure as usual:
#         obs_200_400 = [o for o in obs_uv if 27000 <= o[2] <= 32000]
#         if not obs_200_400:
#             continue

#         u_vals = [o[1] for o in obs_200_400 if o[0] == 139]
#         v_vals = [o[1] for o in obs_200_400 if o[0] == 140]

#         if not u_vals or not v_vals:
#             continue

#         u_mean = np.mean(u_vals)
#         v_mean = np.mean(v_vals)
#         lat_m = obs_200_400[0][3] 
#         lon_m = obs_200_400[0][4] 
#         result = (lat_m, lon_m, u_mean, v_mean)

#         if has_126:
#             results_temp.append(result)
#         else:
#             results_pilot.append(result)

#         # Pseudo pilot: use TEMP data but recalculate pressure using IFS formula
#         # IFS formula: p = 1013.25 * exp(-z/7000) [hPa], z in meters
#         # Only for TEMP data (has_126)
#         if has_126:
#             obs_temp = [o for o in obs if o[0] in (139, 140)]
#             if obs_temp:
#                 pseudo_obs = []
#                 for o in obs_temp:
#                     # o[2] is pressure in hPa, need to estimate z from p
#                     # Invert IFS formula: z = -7000 * ln(p / 1013.25)
#                     p_hpa = o[2]
#                     z = -7000 * np.log(p_hpa / 1013.25)
#                     pseudo_p = 1013.25 * np.exp(-z / 7000)
#                     if 270 <= pseudo_p <= 320:  # pressure in hPa
#                         pseudo_obs.append(o)
#                 if pseudo_obs:
#                     u_vals_pseudo = [o[1] for o in pseudo_obs if o[0] == 139]
#                     v_vals_pseudo = [o[1] for o in pseudo_obs if o[0] == 140]
#                     if u_vals_pseudo and v_vals_pseudo:
#                         u_mean_pseudo = np.mean(u_vals_pseudo)
#                         v_mean_pseudo = np.mean(v_vals_pseudo)
#                         lat_m_pseudo = pseudo_obs[0][3]
#                         lon_m_pseudo = pseudo_obs[0][4]
#                         results_pseudo_pilot.append((lat_m_pseudo, lon_m_pseudo, u_mean_pseudo, v_mean_pseudo))

#     # Aggregate results
#     def aggregate(results):
#         if not results:
#             return [np.nan, np.nan, np.nan, np.nan]
#         lat = [entry[0] for entry in results]
#         lon = [entry[1] for entry in results]
#         u = [entry[2] for entry in results]
#         v = [entry[3] for entry in results]
#         return [np.mean(lat), np.mean(lon), np.mean(u), np.mean(v)]

#     results_pilot = aggregate(results_pilot)
#     results_temp = aggregate(results_temp)
#     results_pseudo_pilot = aggregate(results_pseudo_pilot)

#     return [results_temp, results_pilot, results_pseudo_pilot]


# files = glob.glob("/mnt/users/scratch/leo/scratch/converted_v29/1985/*11035*.nc", recursive=True)[:]

# process_station_file(files[0])  # Example call to test the function

def z_to_p_ifs(h): # geopotential height (m^/s^2) to pressure (Pa)
    a = 5.252368255329
    b = 44330.769230769
    c = 0.000157583169442
    ptro = 226.547172
    po = 1013.25
    g = 9.80665
 
    h /= g
    if h != h:

        p = h

    elif h > 11000.0:

        y = -c * (h - 11000.0)

        p = ptro * np.exp(y)

    else:

        y = 1.0 - h / b

        p = po * (y**a)

    return p * 100. # we want Pa


def process_station_file_pseudo_pilot(file):

    statid = file.split('/')[-1].split('.')[0]
    with h5py.File(file, 'r') as fl:
        x = fl['observations_table']['observed_variable'][:]
        val = fl['observations_table']['observation_value'][:]
        p = fl['observations_table']['z_coordinate'][:]
        rid_raw = fl['observations_table']['report_id'][:]
        lat_z = fl['observations_table']['latitude'][-1]
        lon_z = fl['observations_table']['longitude'][-1]

    # Convert report_id to strings
    report_ids = ["".join(byte.decode() for byte in rid).strip() for rid in rid_raw]

    # Group by report_id
    from collections import defaultdict
    report_obs = defaultdict(list)

    for var, v, pr, rid in zip(x, val, p, report_ids): #, lat, lon):
        report_obs[rid].append((var, v, pr)) #, la, lo))

    results_temp = []
    results_pilot = []
    results_pseudo_pilot = []

    for rid, obs in report_obs.items():
        observed_vars = [o[0] for o in obs]
        if not 126 in observed_vars:
            continue
        
        # Open a wider window for pressure conversion:
        obs_100_500 = [o for o in obs if 10000 <= o[2] <= 50000]

        u_vals = [o[1] for o in obs_100_500 if o[0] == 139]
        v_vals = [o[1] for o in obs_100_500 if o[0] == 140]
        g_vals = [o[1] for o in obs_100_500 if o[0] == 117]
        p_vals_orig = [o[2] for o in obs_100_500 if o[0] == 117]

        # Ensure all data is of the same shape:
        if len(g_vals) != len(u_vals) or len(g_vals) == 0:
            continue
        
        # Convert geopotential height to pressure using IFS formula
        p_vals = [z_to_p_ifs(i_g) for i_g in g_vals]

        # Recreate obs with only u_vals, v_vals, lat, lon, and use p_vals as pressure
        obs_uv = []
        i_g_v = 0
        i_g_u = 0
        for idx, (var, v, pr) in enumerate(obs_100_500): # , la, lo
            if var in [139]:
                obs_uv.append((var, v, p_vals[i_g_u])) # , la, lo
                i_g_u += 1
            if var in [140]:
                obs_uv.append((var, v, p_vals[i_g_v])) # , la, lo
                i_g_v += 1
                
        # Filter pressure as usual:
        obs_200_400 = [o for o in obs_uv if 27000 <= o[2] <= 32000]
        if not obs_200_400:
            continue

        u_vals = [o[1] for o in obs_200_400 if o[0] == 139]
        v_vals = [o[1] for o in obs_200_400 if o[0] == 140]

        if not u_vals or not v_vals:
            continue

        u_mean = np.mean(u_vals)
        v_mean = np.mean(v_vals)
        lat_m = lat_z # obs_200_400[0][3] 
        lon_m = lon_z #obs_200_400[0][4] 
        result = (lat_m, lon_m, u_mean, v_mean)

        results_temp.append(result)

    lat = []
    lon = []
    u = []
    v = []
    if len(results_temp) <= 10:
        lat.append(np.nan) 
        lon.append(np.nan) 
        u.append(np.nan) 
        v.append(np.nan) 
    else:
        for entry in results_temp:
            lat.append(entry[0]) 
            lon.append(entry[1]) 
            u.append(entry[2]) 
            v.append(entry[3]) 
    results_temp = [np.mean(lat), np.mean(lon), np.mean(u), np.mean(v), statid]    

    return [results_temp]

files = glob.glob("/mnt/users/scratch/leo/scratch/converted_v29/1985/0-20999-0-XXBB_CEUAS_merged_v3.nc", recursive=True)[:]
process_station_file_pseudo_pilot(files[0])