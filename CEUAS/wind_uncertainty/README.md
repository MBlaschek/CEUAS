# Estimation of the wind uncertainties

  Author: Federico Ambrogi, federico.ambrogi@univie.ac.at



## Introduction
We describe here the scripts and the necessary steps to extract the covariance matrices
and reproduce the results presented in the deliverable xxx
The required input files are stored in the *data* directory:
- ERA5_1_10393_t.nc
- ERA5_1_10393_u.nc
- ERA5_1_10393_v.nc
These are netCDF files, converted from the ECMWF odb files, containing the observation data for temperature, u- and v-compononents of the wind,
for the Lindenberg station (taken as example in the document).

All the other intermediate files processed during the various steps of the analysis will be also 
stored in the *data* directory.

## Workflow
The workflow is the following:
1. extract_speed_direction_netCDF.py 
2. extract_covariance_analyseOutliers.py
3. analyse_errors.py 

The module wind_utils.py contains classes and variables that are used vy the other scripts.

### Extracting the wind speed and direction from the u- and v-components
Calling the script
` python extract_speed_direction_netCDF.py` 
will generate netCDF files containg the wind speed and direction, to be found in the *data* directory.
These files are needed for the following analysis.

### Extracting the covariance matrices, and analyse outliers.
Since the Desroziers diagnostics is based on the analysis of time avergaes of (cross)covariances of the vectors analysis departure
 and background departures, it is is convenient to extract such matrices and store them in a numpy file.
This is done with
`python extract_covariance_analyseOutliers.py`
which creates the file *covariance_matrices.npy* in the *data* directory.
It will create the file, by default, only if not found in the *data* directory.
With the option
`python extract_covariance_analyseOutliers.py [-f True] `
it forces the creation of the file (note that it will replace the old file, if it exists).
With the option
`python extract_covariance_analyseOutliers.py [-o True] `
the code will creates plots to analyse the outliers, stored in *plots/outliers*
Note that the outliers analysis is not presented in the document, and can be performed
only if creating the matrix.

### Extract the Desroziers errors
BY calling
` python analyse_errors.py [-e True][-c True]`
the covariance matric plots and the errors distribution plots will be created inside the *plots* directory.
Note that at least of the two option must be flagged as True, otherwise the script will quit.
It will also create the numpy file *data/means_std.npy* which contains
the standard deviation and means obtained from the distribution of errors
calculated with the Desroziers statistics.

### Analyse the errror distributions for all the pressure levels
`analyse_errorsDistributions.py`
will create the plots of the mean values and standard deviation of the errors, 
evaluated by *analyse_errors.py* and stored in the *data/means_std.npy* file
The results are stored in *plots/mean_std_plevels*



 