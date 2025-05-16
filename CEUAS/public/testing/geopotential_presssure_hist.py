import glob
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd 
pd.options.mode.chained_assignment = None
import numpy as np
import h5py
import hdf5plugin
import ray
import pickle
import multiprocessing
from functools import partial
import matplotlib.pyplot as plt
import os, sys

from collections import Counter
from scipy.interpolate import interp1d


sys.path.append(os.getcwd()+'/../../public/cds-backend/code/')
import cds_eua4 as eua

#!module load odc
import pyodc
# import codc

from scipy.stats import linregress

import csv

# Ray setup with 40 cores:
ray.init(num_cpus=40)


def find_closest_date(dates, target):
    return dates.loc[(dates - target).abs().idxmin()]

def calculate_pressure_layered(p0, L, T, z, R=287.053, g=9.80665):
    """
    Calculate the pressure (p) layer by layer for a given atmosphere model.
    
    Parameters:
    - p0: float, pressure at the start of the first layer (in Pascals)
    - L: list of float, temperature lapse rates for each layer (in K/m)
    - z: list of float, thickness of each layer (in meters)
    - T: list of float, temperature at the start of each layer (in Kelvin)
    - R: float, specific gas constant for dry air (default: 287.053 J/(kg·K))
    - g: float, acceleration due to gravity (default: 9.80665 m/s²)
    
    Returns:
    - list of float, pressures at the end of each layer
    """
    p_out = [] 

    exponent = g / (L * R) 

    for i in range(len(z)):
        if i == 0:
            p_out.append(p0)
        else:
            if np.abs(L[i]) > 0.001:  # If L is non-zero, use the polytropic formula
                p_next = p_out[i-1] * (1 - (L[i] * (z[i] - z[i-1])) / (T[i])) ** exponent[i]
            else:  # If L close to 0, use the isothermal layer formula
                p_next = p_out[i-1] * np.exp(-g * (z[i] - z[i-1]) / (R * T[i-1]))
            p_out.append(p_next)
    
    return p_out

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

def p_to_z_ifs(p): # pressure (hPa) to geopotential height (m^/s^2) 
    a = 5.252368255329
    b = 44330.769230769
    c = 0.000157583169442
    ptro = 226.547172
    po = 1013.25
    g = 9.80665
 
    if p != p:

        h = p

    elif p < 226.5:
        
        y = np.log(p/ptro)

        h = y /(-c) + 11000.

    else:
        
        y = (p / po) ** (1. / a)
        
        h = (y - 1) * (-b)


    h *= g
    return h



def geopotential_to_pressure(P1, Z1, Z2, Tm):
    """Convert geopotential height and temperature to pressure.
    
    Parameters:
        P1 (float): Pressure at height Z1 (Pa)
        Z1 (float): Geopotential height at P1 (m)
        Z2 (float): Geopotential height at P2 (m)
        Tm (float): Mean temperature between Z1 and Z2 (K)
    
    Returns:
        float: Estimated pressure at Z2 (Pa)
    """
    g = 9.80665  # Gravity (m/s²)
    R = 287.05   # Specific gas constant for dry air (J/kg·K)
    
    # Compute pressure at Z2 using the hypsometric equation
    P2 = P1 * np.exp(-g * (Z2 - Z1) / (R * Tm))
    
    return P2

def compute_pressure_profile(temperature, geopotential, P_surface):
    """Compute pressure profile from temperature and geopotential arrays.
    
    Parameters:
        temperature (array-like): Temperature at different heights (K)
        geopotential (array-like): Corresponding geopotential heights (m)
        P_surface (float): Surface pressure (Pa)
    
    Returns:
        np.ndarray: Computed pressure values at each geopotential height.
    """
    temperature = np.array(temperature)
    geopotential = np.array(geopotential)
    
    pressures = np.zeros_like(geopotential)  # Initialize pressure array
    pressures[0] = P_surface  # Set surface pressure
    
    for i in range(1, len(geopotential)):
        # Compute mean temperature in the layer
        Tm = (temperature[i] + temperature[i - 1]) / 2.0  
        
        # Compute pressure at height Z[i]
        pressures[i] = geopotential_to_pressure(pressures[i - 1], geopotential[i - 1], geopotential[i], Tm)
    
    return pressures

from collections import defaultdict

# Function to process individual file
def find_pilot(file):
    r_temp = 0
    r_pilot = 0
    with h5py.File(file, 'r') as fl:
        x = fl['era5fb']['reportype'][:]
        xx = np.unique(x, return_counts=True)
        pilot = 0
        temp = 0

        if 16068 in xx[0]:
            pilot += int(xx[1][xx[0] == 16068])
        if 16013 in xx[0]:
            pilot += int(xx[1][xx[0] == 16013])

        if 16022 in xx[0]:
            temp += int(xx[1][xx[0] == 16022])
        if 16045 in xx[0]:
            temp += int(xx[1][xx[0] == 16045])
        
    if pilot >= 50:
        r_pilot = 1
    if temp >= 50:
        r_temp = 1
    return [r_temp, r_pilot]

# Function to process individual file -> second option
def find_pilot(file):
    with h5py.File(file, 'r') as fl:
        x = fl['observations_table']['observed_variable'][:]
        y = fl['observations_table']['report_id'][:]
        
    # Group observed variables by report_id
    report_vars = defaultdict(set)

    for var, rid in zip(x, y):
        rid = b"".join(rid).decode("utf-8").strip()
        report_vars[rid].add(var)

    temp = 0
    pilot = 0

    for vars in report_vars.values():
        if 126 in vars:
            temp += 1
        else:
            pilot += 1

    # Convert counts into binary indicators (as in your original version)
    r_temp = int(temp >= 50)
    r_pilot = int(pilot >= 50)

    return [r_temp, r_pilot]

# Ray-remote version of yearly processing
@ray.remote
def process_year(year):
    print(f'Starting {year}!')
    temp_sum = 0
    pilot_sum = 0
    files = glob.glob(f'/mnt/users/scratch/leo/scratch/converted_v29/{year}/*.nc')
    
    for file in files:
        r_temp, r_pilot = find_pilot(file)
        temp_sum += r_temp
        pilot_sum += r_pilot
    
    return year, temp_sum, pilot_sum

# Process years in parallel
years = list(range(1900, 2026))
futures = [process_year.remote(year) for year in years]
results = ray.get(futures)

# Sort results by year
results.sort()

# Extract counts
temp_counts = [res[1] for res in results]
pilot_counts = [res[2] for res in results]

# Save results to CSV
with open("yearly_station_counts_obsvar.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Year", "Temp_Station_Count", "Pilot_Station_Count"])
    writer.writerows(results)
