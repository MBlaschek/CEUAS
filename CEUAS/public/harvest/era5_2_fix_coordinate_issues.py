import numpy as np
import pandas as pd
import glob


# Read data from file 'withwrongcoords'
with open('withwrongcoords.txt', 'r') as file:
    lines = file.readlines()  # Read all lines from the file

# Clean up the lines (remove newline characters) and print initial length
lines = [line.strip() for line in lines]
print(f"Initial number of lines: {len(lines)}")

# Remove duplicate lines
lines = np.unique(lines)
print(f"Number of unique lines: {len(lines)}")

# Prepare data to write by cleaning each line
to_write = []
for line in lines:
    cleaned_line = line.split(' check')[0]  # Remove everything after 'check'
    cleaned_line = cleaned_line.replace('_', '').replace('[', '').replace(']', '').replace(',', '')  # Remove undesired characters
    to_write.append(cleaned_line)

# Further process the cleaned data
final = []
for cleaned in to_write:
    parts = cleaned.split(' ')  # Split line into parts
    if len(parts) > 4:  # If there are more than 4 parts, concatenate the first three and the sixth part
        to_conc = f"{parts[0]} {parts[1]} {parts[2]} {parts[5]}"
    else:
        to_conc = cleaned  # If less than 4 parts, keep the line as it is
    final.append(to_conc)

# Remove duplicates again
final = np.unique(final)

# Write cleaned data to a new file 'cleaned_era5_2_latlon.txt'
with open('cleaned_era5_2_latlon.txt', 'w') as file:
    file.writelines([line + '\n' for line in final])  # Add newline after each line for proper file formatting

# Load the cleaned data into a pandas DataFrame
df = pd.read_csv('cleaned_era5_2_latlon.txt', names=['lat/lon', 'statid', 'year', 'value'], delimiter=' ')
df = df.drop_duplicates()  # Remove duplicate rows

# Initialize dictionary to store the final output
outd = {'lat/lon':[], 'file_lat':[], 'value':[]}

# Iterate over each row in the dataframe
for i in range(len(df)):
    row = df.iloc[i]
    ll = row['lat/lon']  # Latitude/Longitude
    year = row['year']
    sid = row['statid'].replace('0-', '')  # Clean up the statid by removing '0-'
    val = row['value'] * -1  # Invert the value (it seems negative correction is needed)
    
    print(f"Processing: Year={year}, StatID={sid}, Value={val}")

    # Search for relevant files that match the statid in the filename
    for file in glob.glob(f'/mnt/users/scratch/leo/scratch/era5/odbs/2/era5.conv*{sid}*.gz'):
        # Load the gzipped file
        df_fi = pd.read_csv(file, compression='gzip', delimiter='\t')

        # Extract the relevant value from the file using the lat/lon and @hdr suffix
        ll_to_check = df_fi[ll[:3] + "@hdr"].values  # Assuming 'll' is the correct column

        # Check if the value matches (within a tolerance of 0.001)
        if (np.abs(ll_to_check - (val * -1)) < 0.001).any():
            # Append matching results to the output dictionary
            outd['file_lat'].append(file.split('/')[-1])  # Extract file name from full path
            outd['value'].append(val)
            outd['lat/lon'].append(ll)

# Convert the output dictionary to a pandas DataFrame
odf = pd.DataFrame.from_dict(outd)

# Rename columns to match the desired output
odf.columns = ['lat/lon', 'file', 'value']  # Rename columns

# Add an index column, which will be used as 'index'
odf['index'] = odf.index  # The index of the DataFrame will be added as a column

# Reorder the columns to ensure the correct order
odf = odf[['index', 'lat/lon', 'file', 'value']]  # Reorder columns as requested

# Save the DataFrame to a CSV file
odf.to_csv('era5_2_lat_lon_mismatch.csv', index=False)  # Save CSV without the default DataFrame index
