import glob
import os
import numpy as np 
import pandas as pd
import xarray as xr
import re
import math

### Extract data into yearwise station files
# for year in range(1962,2024):
#     print(year)
#     out_dir = "/mnt/users/scratch/uvoggenberger/woudc_obsstore/stat_data/"+str(year)
#     if not os.path.exists(out_dir):
#         os.makedirs(out_dir)

#     files = glob.glob('/mnt/users/scratch/uvoggenberger/woudc_obsstore/*'+ str(year) +'*.nc')
#     year_df = []
#     for fi in files:
#         ds = xr.open_dataset(fi)
#         df = ds.to_dataframe()
#         year_df.append(df)

#     ydf = pd.concat(year_df)
#     sel_df = ydf.drop_duplicates(['longitude|header_table','latitude|header_table'])
#     for ll in range(len(sel_df)):
#         adf = sel_df.iloc[ll]
#         target_lon = adf['longitude|header_table']
#         target_lat = adf['latitude|header_table']
#         stat_df = ydf[np.logical_and(ydf['longitude|header_table'] == target_lon, ydf['latitude|header_table'] == target_lat)]
#         stat_df = stat_df.sort_values(by=['report_timestamp', 'observed_variable'])
#         stat_df.to_csv(out_dir+'/woudc_'+str(target_lat)+'_'+str(target_lon)+'.csv', index=False)
    

### Put files together into station files
# Directory with your data files
data_dir = '/mnt/users/scratch/uvoggenberger/woudc_obsstore/stat_data'
out_dir = '/mnt/users/scratch/uvoggenberger/woudc_obsstore/concat_stations/'

# Regex to extract coordinates from file names
coord_pattern = re.compile(r'woudc_([-+]?\d+\.\d+)_([-+]?\d+\.\d+)\.csv')

# Function to calculate the Haversine distance between two lat/lon pairs
def haversine(lat1, lon1, lat2, lon2):
    # Radius of Earth in kilometers
    R = 6371.0
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    delta_phi = math.radians(lat2 - lat1)
    delta_lambda = math.radians(lon2 - lon1)
    a = math.sin(delta_phi / 2)**2 + math.cos(phi1) * math.cos(phi2) * math.sin(delta_lambda / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    return R * c

# Parse filenames and store coordinates
files = []
for root, _, filenames in os.walk(data_dir):
    for filename in filenames:
        match = coord_pattern.search(filename)
        if match:
            lat, lon = float(match.group(1)), float(match.group(2))
            files.append((os.path.join(root, filename), lat, lon))

# Set distance threshold to 30 km
distance_threshold = 30.0

# Group files into clusters within the 30 km radius
clusters = []
visited = set()

for i, (file1, lat1, lon1) in enumerate(files):
    if file1 in visited:
        continue
    
    # Start a new cluster
    cluster = [file1]
    visited.add(file1)
    
    for j in range(i + 1, len(files)):
        file2, lat2, lon2 = files[j]
        
        # Check if file2 is within 30 km of file1
        if file2 not in visited and haversine(lat1, lon1, lat2, lon2) <= distance_threshold:
            cluster.append(file2)
            visited.add(file2)
    
    # Add the cluster to the list of clusters
    clusters.append(cluster)

# Print results
for idx, cluster in enumerate(clusters, start=1):
    print(f"Cluster {idx} (files within 30 km of each other):")
    concat_out = []
    for file in cluster:
        print(f"  {file}")
        concat_out.append(pd.read_csv(file))
    print()
    file_to_write = pd.concat(concat_out)
    file_to_write = file_to_write.sort_values(by=['report_timestamp', 'observed_variable'])
    file_to_write.to_csv(out_dir+file.split('/')[-1], index=False)
