"""
HDF5 File Validation Script

This script performs a series of tests on HDF5 files to validate their structure, data integrity, 
and content thresholds. It uses `pytest` for parameterized tests and logs the results of failed tests.

Usage:
- Run the tests using `pytest` with the following example commands:
  pytest -v merged_test.py --log-file=merge_test_log.log --log-file-level=INFO
  pytest -n 40 -v merged_test.py --log-file=merge_test_log.log --log-file-level=INFO

Tests performed:
1. File structure validation: Check for the presence of expected groups and variables.
2. Threshold validation: Ensure data values fall within defined acceptable ranges.
3. Source file validation: Verify the existence and consistency of linked source files.
4. (Optional) Data type validation: Check the consistency of data types against predefined expectations.

Constants:
- `DATA_DIRECTORY`: Directory containing HDF5 files for validation.
- `HARVEST_DIRECTORY`: Directory for additional resources (unused here).
- Threshold values for various meteorological and observational parameters.
- Expected group and variable names in HDF5 files.
"""


import h5py
import os
import numpy as np
import pandas as pd
import pytest
import glob
import hdf5plugin


DATA_DIRECTORY = "/mnt/users/scratch/leo/scratch/converted_v23/long/"
HARVEST_DIRECTORY = "/mnt/users/scratch/uvoggenberger/CUON_HARVEST/"


# Threshold constants for data validation
T_MIN, T_MAX = 150.0, 350.0  # Temperature [K]
RH_MIN, RH_MAX = 0.0, 1.0  # Relative Humidity [1]
SH_MIN, SH_MAX = 0.0, 0.040  # Specific Humidity [kg/kg]
DP_MIN, DP_MAX = 150.0, 350.0  # Dewpoint Temperature [K]
DPD_MIN, DPD_MAX = 0.0, 100.0  # Dewpoint Depression [K]
WD_MIN, WD_MAX = 0.0, 360.0  # Wind Direction [°]
WS_MIN, WS_MAX = 0.0, 150.0  # Wind Speed [m/s]
U_MIN, U_MAX = -150.0, 150.0  # Eastward Wind Speed [m/s]
V_MIN, V_MAX = -150.0, 150.0  # Northward Wind Speed [m/s]
GP_MIN, GP_MAX = 0.0, 500000.0  # Geopotential [m²/s²]
P_MIN, P_MAX = 0.0, 110000.0  # Pressure [Pa]



EXPECTED_GROUPS = ['advanced_homogenisation', 'crs', 'era5fb', 'header_table', 'observations_table', 
                   'recordindices', 'sensor_configuration', 'source_configuration', 'station_configuration', 
                   'station_configuration_codes', 'station_type', 'units', 'z_coordinate_type']
EXPECTED_VARIABLES = {
    "advanced_homogenisation": ['RAOBCORE_bias_estimate', 'RASE_bias_estimate', 'RICH_bias_estimate', 'RISE_bias_estimate', 'humidity_bias_estimate', 'wind_bias_estimate'],
    "crs": ['crs', 'crs_len', 'description'],
    "era5fb": ['albedo@modsurf', 'an_depar@body', 'an_depar@surfbody_feedback', 'an_sens_obs@body', 
               'andate', 'antime', 'biascorr@body', 'biascorr_fg@body', 'bufrtype@hdr', 'class', 
               'codetype@hdr', 'collection_identifier@conv', 'date@hdr', 'datum_anflag@body', 
               'datum_event1@body', 'datum_rdbflag@body', 'datum_sfc_event@surfbody_feedback', 
               'datum_status@body', 'datum_status@surfbody_feedback', 'eda_spread@errstat', 
               'entryno@body', 'expver', 'fg_depar@body', 'fg_depar@offline', 'fg_depar@surfbody_feedback', 
               'fg_error@errstat', 'final_obs_error@errstat', 'groupid@hdr', 'lat@hdr', 'lon@hdr', 
               'lsm@modsurf', 'lsm@surfbody_feedback', 'numtsl@desc', 'obs_error@errstat', 'obstype@hdr', 
               'obsvalue@body', 'orography@modsurf', 'ppcode@conv_body', 'qc_pge@body', 'report_event1@hdr', 
               'report_rdbflag@hdr', 'report_status@hdr', 'reportype', 'seaice@modsurf', 'sensor@hdr', 
               'seqno@hdr', 'snow_density@surfbody_feedback', 'snow_depth@modsurf', 
               'snow_depth@surfbody_feedback', 'sonde_type@conv', 'source@hdr', 'stalt@hdr', 'statid@hdr', 
               'station_type@conv', 'stream', 'subtype@hdr', 'time@hdr', 
               'timeseries_index@conv', 'timeslot@timeslot_index', 'tsfc@modsurf', 'type', 
               'unique_identifier@conv', 'varbc_ix@body', 'varno@body', 'vertco_reference_1@body', 
               'vertco_reference_2@body', 'vertco_type@body', 'windspeed10m@modsurf'],
    "header_table": ['application_area', 'comments', 'crs', 'duplicate_status', 'duplicates', 
                     'events_at_station', 'height_of_station_above_local_ground', 
                     'height_of_station_above_sea_level', 'height_of_station_above_sea_level_accuracy', 
                     'history', 'instrument', 'latitude', 'location_accuracy', 'location_method', 
                     'location_quality', 'longitude', 'number_of_pressure_levels', 'observing_programme', 
                     'owner', 'platform_sub_type', 'platform_type', 'primary_station_id', 
                     'primary_station_id_scheme', 'processing_codes', 'processing_level', 'product_name', 
                     'product_version', 'profile_id', 'record_timestamp', 'references', 'region', 
                     'report_duration', 'report_id', 'report_meaning_of_timestamp', 'report_quality', 
                     'report_synoptic_time', 'report_time_accuracy', 'report_time_quality', 
                     'report_time_reference', 'report_timestamp', 'report_type', 'sea_level_datum', 
                     'source_id', 'source_record_id', 'station_course', 'station_heading', 'station_name', 
                     'station_record_number', 'station_speed', 'station_type', 'sub_region'],
    "observations_table": ['adjustment_id', 'advanced_assimilation_feedback', 'advanced_homogenisation', 
                           'advanced_qc', 'advanced_uncertainty', 'bbox_max_latitude', 'bbox_max_longitude', 
                           'bbox_min_latitude', 'bbox_min_longitude', 'code_table', 'conversion_flag', 
                           'conversion_method', 'crs', 'data_policy_licence', 'date_time', 'date_time_meaning', 
                           'exposure_of_sensor', 'latd', 'latitude', 'location_method', 
                           'location_precision', 'lond', 'longitude', 'numerical_precision', 
                           'observation_duration', 'observation_height_above_station_surface', 'observation_id', 
                           'observation_value', 'observed_variable', 'original_code_table', 
                           'original_precision', 'original_units', 'original_value', 'processing_level', 
                           'quality_flag', 'report_id', 'secondary_value', 'secondary_variable', 
                           'sensor_automation_status', 'sensor_id', 'source_id', 'spatial_representativeness', 
                           'station_elevation', 'timed', 'traceability', 
                           'units', 'value_significance', 'z_coordinate', 'z_coordinate_method', 
                           'z_coordinate_type'],
    "recordindices": ['0', '106', '107', '117', '126', '137', '138', '139', '140', '34', '39', 'recordtimestamp'],
    "sensor_configuration": ['calibration_date', 'calibration_status', 'comments', 'date_end', 'date_start', 
                             'observing_method', 'optional_data', 'sampling_strategy', 'sensor_id'],
    "source_configuration": ['source_file'],
    "station_configuration": ['alternative_name', 'bbox_max_latitude', 'bbox_max_longitude', 
                              'bbox_min_latitude', 'bbox_min_longitude', 'city', 'comment', 'contact', 
                              'elevation', 'end_date', 'latitude', 'local_gravity', 'longitude', 
                              'measuring_system_id', 'measuring_system_model', 'metadata_contact', 
                              'metadata_contact_role', 'observed_variables', 'observing_frequency', 
                              'operating_institute', 'operating_territory', 'optional_data', 'platform_sub_type', 
                              'platform_type', 'primary_id', 'primary_id_scheme', 'record_number', 
                              'reporting_time', 'role', 'secondary_id', 'secondary_id_scheme', 'start_date', 
                              'station_abbreviation', 'station_automation', 'station_crs', 'station_name', 
                              'station_type', 'telecommunication_method'],
    "station_configuration_codes": ['abbreviation', 'code_value', 'description', 'field_id', 'field_name', 
                                    'station_configuration_codes_len'],
    "station_type": ['description', 'station_type_len', 'type'],
    "units": ['abbreviation', 'base_units', 'name', 'units', 'units_len'],
    "z_coordinate_type": ['description', 'type', 'z_coordinate_type_len'],
}

# Load the data type table
data_types = pd.read_csv('encodings.txt', delimiter='\t', names=['variable', 'group', 'dtype'])

# 1. Setup and file selection
def get_all_files(directory, file_extension):
    """Helper function to get all files with a specific extension in a directory."""
    file_paths = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(file_extension):
                file_paths.append(os.path.join(root, file))
    return file_paths

hdf5_files = get_all_files(DATA_DIRECTORY, ".nc")# [:400] #! Debugging limit!


# 2. Remove old log files
@pytest.mark.parametrize("hdf5_file_path", hdf5_files)
def remove_old_log(hdf5_file_path):
    if os.path.exists(hdf5_file_path):
        os.remove(hdf5_file_path)


# 3. File Structure Tests
@pytest.mark.parametrize("hdf5_file_path", hdf5_files)
def test_dataset_structure_hdf5(hdf5_file_path):
    """Test that all required datasets and groups are present in an HDF5 file."""
    
    # Initialize a logger to capture any errors during the test
    logger = []
    
    # Define the log file name dynamically based on the current file being tested
    log_name = "./logs/log_" + hdf5_file_path.split("/")[-1] + ".log"

    # Open the HDF5 file in read-only mode
    with h5py.File(hdf5_file_path, "r") as f:
        # Iterate through the list of expected groups to validate their existence
        for group in EXPECTED_GROUPS:
            try:
                # Check if the group exists in the file
                assert group in f.keys(), f"[{group}] is missing from file."
                
                # If the group exists, validate that all expected variables within the group are present
                for variable in EXPECTED_VARIABLES[group]:
                    assert variable in f[group].keys(), f"[{group}][{variable}] is missing from file."
            
            # Capture any assertion errors and log them for review
            except AssertionError as e:
                logger.append(f"File {hdf5_file_path} failed structure test: {e}")
    
    # Write the error log to a file for the current HDF5 file
    with open(log_name, "a") as file:
        file.writelines(line + "\n" for line in logger)


# 4. Threshold Tests
@pytest.mark.parametrize("hdf5_file_path", hdf5_files)
def test_temperature_thresholds_hdf5(hdf5_file_path):
    """
    Ensure temperature values in the HDF5 file are within defined realistic bounds.
    
    Parameters:
    - hdf5_file_path: Path to the HDF5 file being tested.
    """
    # Initialize a logger to capture any errors or issues
    logger = []
    # Generate a log file specific to the current HDF5 file
    log_name = "./logs/log_" + hdf5_file_path.split("/")[-1] + ".log"
    
    # Open the HDF5 file for reading
    with h5py.File(hdf5_file_path, "r") as f:
        # Check if the temperature variable ("126") exists in the "recordindices" group
        if "126" in f["recordindices"].keys():
            try:
                # Extract the indices for temperature data
                t_idx = f["recordindices"]["126"]
                
                # Use the extracted indices to retrieve the temperature values
                temperature_data = f["observations_table"]["observation_value"][t_idx[0]:t_idx[-1]]
                
                # Validate that all temperature values are above the minimum threshold (T_MIN)
                assert np.nanmin(temperature_data) >= T_MIN, "Temperature below minimum threshold."
                
                # Validate that all temperature values are below the maximum threshold (T_MAX)
                assert np.nanmax(temperature_data) <= T_MAX, "Temperature above maximum threshold."
            
            # Handle and log any assertion errors
            except AssertionError as e:
                logger.append(f"File {hdf5_file_path} failed structure test: {e}")
    

        # rh
        if "138" in f["recordindices"].keys():
            try:
                t_idx = f["recordindices"]["138"]
                humidity_data = f["observations_table"]["observation_value"][t_idx[0]:t_idx[-1]]
                assert np.nanmin(humidity_data) >= RH_MIN, "Relative humidity below minimum threshold."
                assert np.nanmax(humidity_data) <= RH_MAX, "Relative humidity above maximum threshold."
            except AssertionError as e:
                logger.append(f"File {hdf5_file_path} failed structure test: {e}")

        #sh 
        if "39" in f["recordindices"].keys():
            try:
                t_idx = f["recordindices"]["39"]
                humidity_data = f["observations_table"]["observation_value"][t_idx[0]:t_idx[-1]]
                assert np.nanmin(humidity_data) >= SH_MIN, "Specific humidity below minimum threshold."
                assert np.nanmax(humidity_data) <= SH_MAX, "Specific humidity above maximum threshold."
            except AssertionError as e:
                logger.append(f"File {hdf5_file_path} failed structure test: {e}")

        #dp 
        if "137" in f["recordindices"].keys():
            try:
                t_idx = f["recordindices"]["137"]
                humidity_data = f["observations_table"]["observation_value"][t_idx[0]:t_idx[-1]]
                assert np.nanmin(humidity_data) >= DP_MIN, "Dewpoint temperature below minimum threshold."
                assert np.nanmax(humidity_data) <= DP_MAX, "Dewpoint temperature above maximum threshold."
            except AssertionError as e:
                logger.append(f"File {hdf5_file_path} failed structure test: {e}")

        # #dpd  
        # if "?" in f["recordindices"].keys():
        #     try:
        #         t_idx = f["recordindices"]["?"]
        #         humidity_data = f["observations_table"]["observation_value"][t_idx[0]:t_idx[-1]]
        #         assert np.min(humidity_data) >= DPD_MIN, "Dewpoint depression below minimum threshold."
        #         assert np.max(humidity_data) <= DPD_MAX, "Dewpoint depression above maximum threshold."
        #     except AssertionError as e:
        #         logger.append(f"File {hdf5_file_path} failed structure test: {e}")

        # wdir
        if "106" in f["recordindices"].keys():
            try:
                t_idx = f["recordindices"]["106"]
                wind_data = f["observations_table"]["observation_value"][t_idx[0]:t_idx[-1]]
                assert np.nanmin(wind_data) >= WD_MIN, "Wind direction below minimum threshold."
                assert np.nanmax(wind_data) <= WD_MAX, "Wind direction above maximum threshold."
            except AssertionError as e:
                logger.append(f"File {hdf5_file_path} failed structure test: {e}")

        # wspeed
        if "107" in f["recordindices"].keys():
            try:
                t_idx = f["recordindices"]["107"]
                wind_data = f["observations_table"]["observation_value"][t_idx[0]:t_idx[-1]]
                assert np.nanmin(wind_data) >= WS_MIN, "Wind speed below minimum threshold."
                assert np.nanmax(wind_data) <= WS_MAX, "Wind speed above maximum threshold."
            except AssertionError as e:
                logger.append(f"File {hdf5_file_path} failed structure test: {e}")

        # u
        if "139" in f["recordindices"].keys():
            try:
                t_idx = f["recordindices"]["139"]
                wind_data = f["observations_table"]["observation_value"][t_idx[0]:t_idx[-1]]
                assert np.nanmin(wind_data) >= U_MIN, "Eastward wind speed below minimum threshold."
                assert np.nanmax(wind_data) <= U_MAX, "Eastward wind speed above maximum threshold."
            except AssertionError as e:
                logger.append(f"File {hdf5_file_path} failed structure test: {e}")

        # v
        if "138" in f["recordindices"].keys():
            try:
                t_idx = f["recordindices"]["138"]
                wind_data = f["observations_table"]["observation_value"][t_idx[0]:t_idx[-1]]
                assert np.nanmin(wind_data) >= V_MIN, "Northward wind speed below minimum threshold."
                assert np.nanmax(wind_data) <= V_MAX, "Northward wind speed above maximum threshold."
            except AssertionError as e:
                logger.append(f"File {hdf5_file_path} failed structure test: {e}")

        # gp 
        if "117" in f["recordindices"].keys():
            try:
                t_idx = f["recordindices"]["117"]
                geopotential_data = f["observations_table"]["observation_value"][t_idx[0]:t_idx[-1]]
                assert np.nanmin(geopotential_data) >= GP_MIN, "Geopotential below minimum threshold."
                assert np.nanmax(geopotential_data) <= GP_MAX, "Geopotential above maximum threshold."
            except AssertionError as e:
                logger.append(f"File {hdf5_file_path} failed structure test: {e}")

        # p
        try:
            pressure_data = f["observations_table"]["z_coordinate"][:]
            assert np.nanmin(pressure_data) >= P_MIN, "Pressure below minimum threshold."
            assert np.nanmax(pressure_data) <= P_MAX, "Pressure above maximum threshold."
        except AssertionError as e:
            logger.append(f"File {hdf5_file_path} failed structure test: {e}")


    with open(log_name, "a") as file:
        file.writelines(line + "\n" for line in logger)


# 5. Source Tests
@pytest.mark.parametrize("hdf5_file_path, dataset_name, expected_dtype", hdf5_files)
def test_source(hdf5_file_path):
    """Check that each dataset has a working link to its source file."""
    
    # Log file setup for capturing errors
    log_name = "./logs/log_" + hdf5_file_path.split("/")[-1] + ".log"
    logger = []

    # Open the main HDF5 file to validate source links
    with h5py.File(hdf5_file_path, "r") as f:
        # Extract timestamps from the main file
        f_timestamps = f['recordindices']['recordtimestamp'][:]
        already_checked_files = []

        # Verify all source files listed in the "source_configuration" group
        for i in np.unique(f["source_configuration"]['source_file'][:], axis=0):
            # Decode and clean up the source file path
            i = i.tobytes().decode('utf-8').replace(' ', '')
            print(i)
            already_checked_files.append(i)

            # Check if the source file exists
            try:
                assert os.path.exists(i), f"Expected {i} to exist, file not found."
            except AssertionError as e:
                logger.append(f"File {hdf5_file_path} failed source file test: {e}")
            
            # Validate that all record timestamps from the source file match with the main file
            with h5py.File(i, "r") as f_source:
                f_source_rts = f_source['recordtimestamp'][:]
                for fsrts in f_source_rts:
                    idx = np.searchsorted(f_timestamps, fsrts)
                    # Compute the minimum difference between timestamps
                    if (idx > 0) and (idx < (len(f_timestamps) - 1)):
                        result = np.min(f_timestamps[idx - 1:idx + 1] - fsrts)
                    elif idx == 0:
                        result = np.min(f_timestamps[idx:idx + 1] - fsrts)
                    else:
                        result = np.min(f_timestamps[idx - 1:idx] - fsrts)
                    
                    # Ensure all timestamps are within an acceptable range (e.g., 7200 seconds)
                    try:
                        assert np.abs(result) <= 7200, f"Expected all record timestamps from source file to exist. {idx} are missing timestamps from {i}."
                    except AssertionError as e:
                        logger.append(f"File {hdf5_file_path} failed source file test: {e}")

        # Validate data from harvested files not listed as sources
        # Extract the station ID from the main file name
        stat_id = hdf5_file_path.split('/')[-1].split('_')[0]
        # Find all directories in the harvest directory matching the station ID
        matching_dirs = glob.glob(HARVEST_DIRECTORY + 'harvest*/*/' + stat_id)
        for md in matching_dirs:
            # Identify matching harvested files
            matching_files = glob.glob(md + '/*.nc')
            for mf in matching_files:
                # Skip already-checked files
                if mf not in already_checked_files:
                    print(mf)
                    with h5py.File(mf, "r") as f_source:
                        f_source_rts = f_source['recordtimestamp'][:]
                        for fsrts in f_source_rts:
                            idx = np.searchsorted(f_timestamps, fsrts)
                            # Compute timestamp differences
                            if (idx > 0) and (idx < (len(f_timestamps) - 1)):
                                result = np.min(f_timestamps[idx - 1:idx + 1] - fsrts)
                            elif idx == 0:
                                result = np.min(f_timestamps[idx:idx + 1] - fsrts)
                            else:
                                result = np.min(f_timestamps[idx - 1:idx] - fsrts)
                            
                            # Ensure no additional data exists in harvested files
                            try:
                                assert np.abs(result) <= 7200, f"The file {mf} contains data, which is not in final file {stat_id}."
                            except AssertionError as e:
                                logger.append(f"File {mf} failed source file test: {e}")

    # Write all logged errors to a log file for further debugging
    with open(log_name, "a") as file:
        file.writelines(line + "\n" for line in logger)


# # 6. Data Type Tests
# @pytest.mark.parametrize("hdf5_file_path", hdf5_files)
# def test_data_type_hdf5(hdf5_file_path):
#     """Check that each dataset has the expected data type in HDF5."""
#     log_name = "./logs/log_"+hdf5_file_path.split("/")[-1]+".log"
#     logger = []
#     with h5py.File(hdf5_file_path, "r") as f:
#         for group in f.keys():
#             for variable in f[group].keys():
#                 compare_to = data_types[np.logical_and(data_types.group == group, data_types.variable == variable)]
#                 if len(compare_to) > 0:
#                     try:
#                         actual_dtype = f[group][variable].dtype

#                         type_name = compare_to.dtype.values[0].split("'")[1]
#                         expected_dtype = np.dtype(getattr(np, type_name.split('.')[-1]))
                        
#                         assert actual_dtype == expected_dtype, f"Expected {expected_dtype} for {group}/{variable}, got {actual_dtype}."
#                     except AssertionError as e:
#                         logger.append(f"File {hdf5_file_path} failed dtype test: {e}")
#     with open(log_name, "a") as file:
#         file.writelines(line + "\n" for line in logger)


