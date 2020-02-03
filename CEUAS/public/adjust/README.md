# Homogenisation Software

Related Milestones and Deliverables:

Type | Nr | Responsible | Nature | Title| Due | Status | File
---|---|---|---|---|---|---|---
Milestone | MC3S311c_Lot2.2.1.2 | UNIVIE | Homogenization software on GitHub | Demonstrate run on Linux platform, using CDM compliant input | December 2019 | January 2020 | raso_adj_cdm_v0.py
Deliverable | DC3S311c_Lot2.2.4.2| UNIVIE | Software, Report | Homogeneity adjustment software + Adjustments plus user guide and documentation | February 2021 | ontime | 
Milestone | MC3S311c_Lot2.2.1.3 | UNIVIE | Final version of homogenization software on GitHub | Demonstration session, technical report | March 2021 | ontime | 
Deliverable | DC3S311c_Lot2.3.2.2 | UNIVIE | Report, Software | Reproduction of bias adjustment calculation in CDS | March 2021 | ontime |

[Use Interactive Notebook on COLAB (Google)](https://colab.research.google.com/github/MBlaschek/CEUAS/blob/master/CEUAS/public/cds-backend/Example_Homogenization.ipynb)

# Short description

Here the homogenization software employed in the Copernicus Climate Change Service (C3S) - Upper Air Data Service is provided. Please find more details in the script (`raso_adj_cdm_v0.py`) or the example Notebook (`Example_Homogenization.ipynb`). The IPython notebook is intended for demonstration purposes and example usage of the Homogenization software for Radiosonde humidity data and delivers information on the following:
* Downloading of data
* Investigating the data
* Running the Homogenization software
* Investigating results

Here we provide a brief description on the [installation](#Installation), [how to use the script](#How-to-use?) and the [license](#License).

# Installation
The script `raso_adj_cdm_v0.py` is in pure python code with imports of standard packages as follows:

* sys
* os
* numpy
* pandas
* xarray
* numba
* datetime
* getopt
* glob
* matplotlib

```bash
pip install sys os numpy pandas xarray numba datetime getopt glob matplotlib
```

## Homogenization software RISE
The code is provided in order to document the current state, but it is not advisable to use this code at the moment as it is very specific to a certain system and not yet ported for public use. However, a python version that consolidates with the humidity adjustment software is planed and will be delivered at the latest by the end of the contract.

# How to use?

The script can be executed as a script or imported like a module (there is no difference in functionality):

```bash
python raso_adj_cdm_v0.py -h
[MAIN] Executing ...

Run standardized radiosonde homogenisation software on CDM compliant file

raso_adj_cdm_v0.py -h -f [file] -o [name] 

Options:
    -h              Help
    --help      
    -f []           Input CDM compliant file
    --file []       
    -o []           Output name
    --output []
    
Optional Keyword Options:
    --thres []          Threshold value for SNHT, default: 50
    --window []         Moving Window for SNHT, default: 1470 (in days, 4 years)
    --missing []        Maximum allowed missing values in window, default: 600 (in days)
    --min_levels []     Minimum required levels for significant breakpoint, default: 3
    --dist []           Minimum distance between breakpoints, default: 730 (in days, 2 years)
    --sample_size []    Minimum sample size for statistics, default: 130 (in days)
    --borders []        Breakpoint zone, default: 90 (in days)
    --ratio []          Use ratio instead of differences, default: 0 (not)

    --logfile []        Write messages to a log file

Experimental Keyword Options:
    --donotwrite          Returns xarray Dataset
    --interpolate_missing Interpolate Adjustments to non-standard times and pressure levels
    
```

A common way to employ this script would be to download some radiosonde station data from the CDS (or at the moment only from the backend).
This example of a retrieval can be viewed in the Notebook (`Example_Homogenization.ipynb`) or in the example Notebook (`Example.ipynb`) from Deliverable [DC3S311c_Lot2.3.1.1](https://github.com/MBlaschek/CEUAS/tree/master/CEUAS/public/cds-backend). The following retrieves station data for 72357 (Norman USA) for all significant levels with reanalysis feedback information:

```python
import requests, zipfile, io, os, time
t0 = time.time()
#
# This is for Python 3+
#
options = {"statid": "72357",
           "date": [19900101, 20181231],
           "variable": ["temperature", "relative_humidity"],
           "fbstats": ["obs_minus_bg", "obs_minus_an", "bias_estimate"]
           }
#
# Download significant levels
#
if False:
    #
    # Standard pressure levels only (faster retrieval)
    #
    options["pressure_level"] = [1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 85000, 92500, 100000]
#                            
# 
r = requests.post('http://early-upper-air.copernicus-climate.eu',
                  headers={'content-type': 'application/json'},
                  json=options,
                  stream=True)
#
# Check for known Error Message
#
if r.status_code != requests.codes.ok:
    print(r.text)
    raise RuntimeError("Something is not correct with the request")
#
# Create directory
#
if not os.path.isdir('./example_data'):
    os.makedirs('./example_data')
#
# Unzip the Data Stream
#
try:
    z = zipfile.ZipFile(io.BytesIO(r.content))
    print("Retreived files: ")
    print(z.namelist())
    z.extractall(path='./example_data')
except:
    print(r.text)
    print("Error in request")    
print("Time elapsed: ", time.time()-t0, "s")
```

Following this data retrieval, a default command would be to employ the standard homogeneity adjustment procedure on the files in `example_data` and interpolate to non-standard dates and times as well as non-standard pressure levels. At the moment only the variable `hur` (relative humidity) will be adjusted.

```bash
python raso_adj_cdm_v0.py -f "example_data/*[!_out].nc" --adddate --interpolate_missing

[2020-01-30T09:22:37.272539] [INFO] Multiple input files:  example_data/*[!_out].nc  | # 2
[2020-01-30T09:22:37.530517] [INPUT] example_data/dest_72357_air_temperature.nc ta
[2020-01-30T09:22:37.760866] [INPUT] example_data/dest_72357_relative_humidity.nc hur
[2020-01-30T09:22:37.809099] [CONVERT] Converting to DataCube ...
[2020-01-30T09:22:37.811751] [BKP] Observations obs: 1346329
[2020-01-30T09:22:37.811818] [BKP] Trajectories obs: 1346329,trajectory: 21418
[2020-01-30T09:22:37.871434] [TIME] Dates: 1990-01-01 00:00:00  -  2018-12-31 23:04:00 Duplicates: 1324911
[2020-01-30T09:22:37.899379] [CONVERT] Selecting only standard pressure levels
[2020-01-30T09:22:40.110501] [CONVERT] Done
[2020-01-30T09:22:40.399273] [CONVERT] Standard time calculated, Duplicates resolved: 53
[2020-01-30T09:22:40.485993] [CONVERT] Converting to day-night Array [hour x time x pressure]
[2020-01-30T09:22:40.486258] [CHECK] Departures found ta obs_minus_an
[2020-01-30T09:22:40.486298] [CHECK] Departures found ta bias_estimate
[2020-01-30T09:22:40.486326] [CHECK] Departures found hur obs_minus_an
[2020-01-30T09:22:40.486356] [CHECK] Departures found hur bias_estimate
[2020-01-30T09:22:41.514256] [SNHT] Test statistics calculted ta_obs_minus_an_snht
[2020-01-30T09:22:42.564257] Breaks: [3000]
[2020-01-30T09:22:42.570787] [DETECT] Test statistics calculated ta_obs_minus_an_snht_breaks
[2020-01-30T09:22:42.582691] [SNHT] Test statistics calculted hur_obs_minus_an_snht
[2020-01-30T09:22:42.708884] Breaks: [3067 6126 8654]
[2020-01-30T09:22:42.715653] Breaks: [3067 6125 8653]
[2020-01-30T09:22:42.716585] [DETECT] Test statistics calculated hur_obs_minus_an_snht_breaks
[2020-01-30T09:22:42.723885] Breakpoints for  hur_obs_minus_an_snht_breaks
[2020-01-30T09:22:42.723941] [     idx] [     end] [    peak] [   start] [ #]
[2020-01-30T09:22:42.723978] [    3067] 1996-11-18 1998-05-31 1999-12-14 1120
[2020-01-30T09:22:42.724007] [    6126] 2005-06-02 2006-11-13 2008-01-06  920
[2020-01-30T09:22:42.724033] [    8654] 2012-09-27 2013-10-24 2015-07-12 1010
[2020-01-30T09:22:42.724072] [ADJUST] Breakpoints:  3
[2020-01-30T09:22:42.724573] [ADJUST] Sample size: 10 N-Q: 13
[2020-01-30T09:22:42.773872] Breakpoints for  hur_obs_minus_an_snht_breaks
[2020-01-30T09:22:42.773928] [     idx] [     end] [    peak] [   start] [ #]
[2020-01-30T09:22:42.773962] [    3067] 1996-10-10 1998-05-31 1999-10-20 1104
[2020-01-30T09:22:42.773999] [    6125] 2005-06-07 2006-10-31 2007-12-25  903
[2020-01-30T09:22:42.774025] [    8653] 2012-09-07 2013-10-23 2015-03-28  924
[2020-01-30T09:22:42.774054] [ADJUST] Breakpoints:  3
[2020-01-30T09:22:42.774598] [ADJUST] Sample size: 10 N-Q: 13
[2020-01-30T09:22:42.804510] [ALIGN] hour x time -> datetime 
[2020-01-30T09:22:42.862105] [REVERSE] Selecting variables
[2020-01-30T09:22:42.862210] [REVERSE] time to obs
[2020-01-30T09:22:43.210428] [TRAJ] Add trajectory information
[2020-01-30T09:22:43.239575] [OUT] Writing ...
[2020-01-30T09:22:43.240041] [OUTPUT] example_data/dest_72357_relative_humidity_out.nc bias_estimate, obs_minus_an, obs_minus_bg, hur_q, obs_minus_an_snht, obs_minus_an_snht_breaks, trajectory_index, trajectory_label, hur
[2020-01-30T09:23:04.717266] [MAIN] Finished
```

In a similar fashion the script can be imported as a module and used like this:
```python
import raso_adj_cdm_v0
raso_adj_cdm_v0.main(ifile="example_data/*[!_out].nc", adddate=True)
```

or return output to the python environment:
```python
import raso_adj_cdm_v0
data = raso_adj_cdm_v0.main(ifile="example_data/*[!_out].nc", adddate=True, donotwrite=True)
```
```
[2020-01-30T09:23:05.061040] [INFO] Multiple input files:  example_data/*[!_out].nc  | # 2
[2020-01-30T09:23:05.336528] [INPUT] example_data/dest_72357_air_temperature.nc ta
[2020-01-30T09:23:05.572736] [INPUT] example_data/dest_72357_relative_humidity.nc hur
[2020-01-30T09:23:05.625143] [CONVERT] Converting to DataCube ...
[2020-01-30T09:23:05.628032] [BKP] Observations obs: 1346329
[2020-01-30T09:23:05.628096] [BKP] Trajectories obs: 1346329,trajectory: 21418
[2020-01-30T09:23:05.687884] [TIME] Dates: 1990-01-01 00:00:00  -  2018-12-31 23:04:00 Duplicates: 1324911
[2020-01-30T09:23:05.715448] [CONVERT] Selecting only standard pressure levels
[2020-01-30T09:23:08.133193] [CONVERT] Done
[2020-01-30T09:23:08.405790] [CONVERT] Standard time calculated, Duplicates resolved: 53
[2020-01-30T09:23:08.493634] [CONVERT] Converting to day-night Array [hour x time x pressure]
[2020-01-30T09:23:08.493726] [CHECK] Departures found ta obs_minus_an
[2020-01-30T09:23:08.493765] [CHECK] Departures found ta bias_estimate
[2020-01-30T09:23:08.493796] [CHECK] Departures found hur obs_minus_an
[2020-01-30T09:23:08.493833] [CHECK] Departures found hur bias_estimate
[2020-01-30T09:23:08.506615] [SNHT] Test statistics calculted ta_obs_minus_an_snht
[2020-01-30T09:23:08.618187] Breaks: [3000]
[2020-01-30T09:23:08.623898] [DETECT] Test statistics calculated ta_obs_minus_an_snht_breaks
[2020-01-30T09:23:08.635725] [SNHT] Test statistics calculted hur_obs_minus_an_snht
[2020-01-30T09:23:08.748921] Breaks: [3067 6126 8654]
[2020-01-30T09:23:08.755645] Breaks: [3067 6125 8653]
[2020-01-30T09:23:08.756559] [DETECT] Test statistics calculated hur_obs_minus_an_snht_breaks
[2020-01-30T09:23:08.763884] Breakpoints for  hur_obs_minus_an_snht_breaks
[2020-01-30T09:23:08.763936] [     idx] [     end] [    peak] [   start] [ #]
[2020-01-30T09:23:08.763973] [    3067] 1996-11-18 1998-05-31 1999-12-14 1120
[2020-01-30T09:23:08.764002] [    6126] 2005-06-02 2006-11-13 2008-01-06  920
[2020-01-30T09:23:08.764029] [    8654] 2012-09-27 2013-10-24 2015-07-12 1010
[2020-01-30T09:23:08.764060] [ADJUST] Breakpoints:  3
[2020-01-30T09:23:08.764585] [ADJUST] Sample size: 10 N-Q: 13
[2020-01-30T09:23:08.800837] Breakpoints for  hur_obs_minus_an_snht_breaks
[2020-01-30T09:23:08.800889] [     idx] [     end] [    peak] [   start] [ #]
[2020-01-30T09:23:08.800924] [    3067] 1996-10-10 1998-05-31 1999-10-20 1104
[2020-01-30T09:23:08.800954] [    6125] 2005-06-07 2006-10-31 2007-12-25  903
[2020-01-30T09:23:08.800980] [    8653] 2012-09-07 2013-10-23 2015-03-28  924
[2020-01-30T09:23:08.801010] [ADJUST] Breakpoints:  3
[2020-01-30T09:23:08.801242] [ADJUST] Sample size: 10 N-Q: 13
[2020-01-30T09:23:08.831652] [ALIGN] hour x time -> datetime 
[2020-01-30T09:23:08.907494] [REVERSE] Selecting variables
[2020-01-30T09:23:08.907598] [REVERSE] time to obs
[2020-01-30T09:23:09.249673] [TRAJ] Add trajectory information
[2020-01-30T09:23:09.280101] [OUT] Returning data
```
```python
data
<xarray.Dataset>
Dimensions:                       (obs: 1346329, trajectory: 21418)
Coordinates:
    lat                           (obs) float32 35.23 35.23 ... 35.18083
    lon                           (obs) float32 -97.47 -97.47 ... -97.43778
    plev                          (obs) float64 980.0 1e+03 ... 9.25e+04 1e+05
    time                          (obs) datetime64[ns] 1990-01-01 ... 2018-12-31T23:04:00
Dimensions without coordinates: obs, trajectory
Data variables:
    ta_bias_estimate              (obs) float32 -1.340332 -1.340332 ... nan
    ta_obs_minus_an               (obs) float32 -0.173956 0.280334 ... nan
    ta_obs_minus_bg               (obs) float32 -0.039262 0.427206 ... nan
    ta                            (obs) float32 218.9 219.3 nan ... 279.75 nan
    hur_bias_estimate             (obs) float32 nan nan nan nan ... nan nan nan
    hur                           (obs) float32 nan nan nan ... 0.654979 nan
    hur_obs_minus_an              (obs) float32 nan nan nan ... -0.091665 nan
    hur_obs_minus_bg              (obs) float32 nan nan nan ... -0.050518 nan
    hur_q                         (obs) float32 nan nan nan nan ... 0.0 0.0 nan
    hur_obs_minus_an_snht         (obs) float64 0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0
    hur_obs_minus_an_snht_breaks  (obs) float64 0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0
    ta_trajectory_index           (obs) int32 0 0 0 0 ... 21417 21417 21417
    ta_trajectory_label           (trajectory) |S5 b'00551' ... b'23680'
    hur_trajectory_index          (obs) int32 0 0 0 0 ... 21417 21417 21417
    hur_trajectory_label          (trajectory) |S5 b'00551' ... b'23680'
```

For further examples we refer to the Notebook `Example_Homogenization.ipynb`

# License

Generated using Copernicus Climate Change Service Information, 2020
[Copernicus Climate Change Service (C3S), 2020](https://apps.ecmwf.int/datasets/licences/copernicus/)

