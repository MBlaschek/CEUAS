# Import necessary libraries
import requests
import time
import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import transform

# Function to retrieve altitude from an online API based on latitude and longitude
def get_altitude(latitude, longitude):
    url = f"https://api.open-meteo.com/v1/elevation?latitude={latitude}&longitude={longitude}"
    response = requests.get(url)
    data = response.json()
    return data["elevation"][0]

# Load the data from CSV file into a DataFrame
# Replace 'test_orph.csv' with the actual file path if needed
df = pd.read_csv('/srvfs/home/uvoggenberger/CEUAS/CEUAS/meta/inventory_comparison_2/code/station_configuration/CUON_orphans_station_configuration_extended.csv', delimiter='\t') # test_orph.csv', delimiter=',')   # /srvfs/home/uvoggenberger/CEUAS/CEUAS/meta/inventory_comparison_2/code/station_configuration/CUON_orphans_station_configuration_extended.csv', delimiter='\t') # test_orph ','
if not 'calculated_elevation' in df.keys():
    df['calculated_elevation'] = np.full(len(df), np.nan)

# Initialize an empty list to store calculated elevations
calc_elevation = []

# Initialize counters for API request tracking
counter = 0  # Tracks requests between pauses
counter_total = 0  # Tracks total requests

# Iterate over each row in the DataFrame
for i in range(len(df)):
    i_df = df.iloc[i]  # Get the row as a Series
    
    # Check if 'calculated_elevation' is missing and needs to be filled
    if np.isnan(i_df.calculated_elevation):
        
        # Pause requests every 598 API calls to avoid rate limits
        if counter == 598:
            time.sleep(61)
            counter = 0  # Reset the counter after sleeping

        # Process latitude: handle comma separators and remove any parentheses
        if not pd.isna(i_df.latitude):
            if ',' in i_df.latitude:
                i_df.latitude = i_df.latitude.split(',')[0].replace('(', '')
        
        # Process longitude: handle comma separators and remove any parentheses
        if not pd.isna(i_df.longitude):
            if ',' in i_df.longitude:
                i_df.longitude = i_df.longitude.split(',')[0].replace('(', '')

        # Check if latitude and longitude are valid numbers
        if not (np.isnan(float(i_df.latitude)) or np.isnan(float(i_df.longitude))):
            try:
                # Get elevation from the API
                ele = get_altitude(float(i_df.latitude), float(i_df.longitude))
                counter += 1  # Increment request counter for each successful API call
                counter_total += 1  # Increment total request counter
            except:
                # Handle failed API calls by assigning NaN to elevation
                ele = np.nan
        else:
            ele = np.nan  # Assign NaN if coordinates are invalid

        # Append the calculated elevation (or NaN) to the list
        calc_elevation.append(ele)
        
        # Print original elevation and the calculated elevation for comparison
        print(i_df.elevation, calc_elevation[-1])
    else:
        # Use existing calculated elevation if available
        calc_elevation.append(i_df.calculated_elevation)

    # Stop processing if the total request limit is reached
    if counter_total > 4998:
        print('total_requests == 5000, continue next time at: ', i)
        break

# Add the calculated elevations to the DataFrame as a new column
df['calculated_elevation'] = calc_elevation

# Save the updated DataFrame back to the CSV file
df.to_csv('test_orph.csv', index=False)