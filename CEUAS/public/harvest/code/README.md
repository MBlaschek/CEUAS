### Harvesting utilities
## Here we provide a brief description on how to run the data harvesting utilities




# harvest_to_netCDF_converter_leo+federico.py
Main script that converts an input file (form any of the dataset) to CDM compliant netCDF.
It is called by the scripts test_xxx and run_xxx , see below.
Example usage:
python3 [-f][FILE] [-o][OUTPUT DIRECTORY] [-d][DATASET NAME]
python3 -f /raid60/scratch/leo/scratch/era5/odbs/1/era5.conv._82930.gz -d era5_1 -o output_dir


# test_harvest_to_netCDF_converter_leo+federico.py
Runs the script harvest_to_netCDF_converter_leo+federico.py on a selection of files,
as a test. 
The output directory name and the grousp of file can be edited in the script.

Example usage:
python3 test_harvest_to_netCDF_converter_leo+federico.py

# run_xxx
Runs the harvest_to_netCDF_converter_leo+federico.py on the entire datasets available,
that can be specified inside the script.
The user can select the output directory and the number of files to be processed in parallel for each
dataset.

Example usage:
python3 run_harvest_to_netCDF_converter_leo+federico.py 