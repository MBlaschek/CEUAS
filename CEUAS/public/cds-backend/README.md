# Access to the CDS 

Related Milestones and Deliverables:


Type | Nr | Responsible | Nature | Title| Due | Status | File
---|---|---|---|---|---|---|---
Deliverable | DC3S311c_Lot2.3.1.1 | UNIVIE | Software, Report | First access to early upper air data base via CDS | September 2019 | January 2020 | code/* 

# Short Description

The code supplied here allows to run a python web-server (hug) that will be a backend for the Copernicus Climate Change Service (C3S) - Upper Air Service.

The backend expects HTTP POST requests, where the query string must be in JSON format. Table 1 describes the allowed keys and values of the requests. HTTP GET requests may work as well but are accepted only for debugging.

 

The backend returns files, which are either

1. Zip files containing CF 1.7 compliant netCDF4 files (one per station). The default name is `download.zip`.
2. JSON files containing error messages, if a HTTP error occurs

 Both file formats can be dealt with in the CDS toolbox. 

| Identifier       | All possible values                                          | Explanation                                                  |
| ---------------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| `date`           | `[YYYYMMDD,YYYYMMDD]`, `YYYYMMDD`, Integer or String         | Range of dates of radiosonde launches                        |
| `country`        | `[“CCC”,…,”DDD”]`, String, Exclusive with `statid`, `bbox`, String | Country codes of stations to be selected                     |
| `bbox`           | `[lower,left,upper,right]`, Float or String, Exclusive with `statid`, `country` | Boundaries of lat/lon rectangle to select stations           |
| `fbstats`        | `["obs_minus_bg","obs_minus_an","bias_estimate"]`            | ERA5 feedback information                                    |
| `pressure_level` | `[MMMM,…,NNNNN]`, `MMMM`, Integer or String                  | Pressure levels in Pascal. 16 standard pressure levels (10-1000 hPa) or significant levels (if omitted) |
| `statid`         | `[“SSSSS”]`, String, Special value “all”, Exclusive with `country`, `bbox` | WMO or WIGOS station ID                                      |
| `time`           | `[HHMMSS,HHMMSS]`                                            | List of times permitted.                                     |
| `variable`       | `[„temperature“, “u_component_of_wind“, “v_component_of_wind“, “wind_speed”, ”wind_direction”, ”relative_humidity”, ”specific_humidity”]`, String | Meteorological variables                                     |

Here we provide a brief description on the [installation](#Installation), [how to use the script](#How-to-use?) and the [license](#License).
[Use Interactive Notebook on COLAB (Google)](https://colab.research.google.com/github/MBlaschek/CEUAS/blob/master/CEUAS/public/cds-backend/Example.ipynb)


# Installation

There are a few steps to setup the whole server:

1. install hug
2. download the hug server scripts
3. upload data `code/upload_to_vm.ksh`



Ad 1) install the required packages on the server

```bash
pip install hug
pip install falcon pandas xarray numba numpy h5py
```

Ad 2) get the latest version of the code from GitHub:

* [CDS Module for hug](https://raw.githubusercontent.com/MBlaschek/CEUAS/master/CEUAS/public/cds-backend/code/cds_eua.py)

* [hug default launch script](https://raw.githubusercontent.com/MBlaschek/CEUAS/master/CEUAS/public/cds-backend/code/default.py)

Ad 3) Execute the upload [script](https://raw.githubusercontent.com/MBlaschek/CEUAS/master/CEUAS/public/cds-backend/code/upload_to_vm.ksh) for the local merged data base to the Virtual Machine (VM):

```bash
./upload_to_vm.ksh
```



## Start the backend

Launch the [hug](https://www.hug.rest/) server using default parameters:

```bash
hug -f default.py
/#######################################################################\
          `.----``..-------..``.----.
         :/:::::--:---------:--::::://.
        .+::::----##/-/oo+:-##----:::://
        `//::-------/oosoo-------::://.       ##    ##  ##    ##    #####
          .-:------./++o/o-.------::-`   ```  ##    ##  ##    ##  ##
             `----.-./+o+:..----.     `.:///. ########  ##    ## ##
   ```        `----.-::::::------  `.-:::://. ##    ##  ##    ## ##   ####
  ://::--.``` -:``...-----...` `:--::::::-.`  ##    ##  ##   ##   ##    ##
  :/:::::::::-:-     `````      .:::::-.`     ##    ##    ####     ######
   ``.--:::::::.                .:::.`
         ``..::.                .::         EMBRACE THE APIs OF THE FUTURE
             ::-                .:-
             -::`               ::-                   VERSION 2.6.0
             `::-              -::`
              -::-`           -::-
\########################################################################/

 Copyright (C) 2016 Timothy Edmund Crosley
 Under the MIT License

Serving on :8000...
```

# How to use?

One way to access the data from the backend can be a Linux tool called `curl` to download a zip file. 
The request (`--data ...`) is identical in the python request and also the retrieved file is identical. The retrieved file needs to be unzipped.

```bash
curl -H "Content-Type: application/json" -X POST --digest --data '{"statid":"11035","date":[20000101,20000101],"pressure_level":[1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 85000, 92500, 100000],"variable":["temperature","relative_humidity"],"fbstats":["obs_minus_bg","obs_minus_an","bias_estimate"]}' -o example_data/download.zip http://early-upper-air.copernicus-climate.eu
```

```bash
unzip -o example_data/download.zip
Archive:  example_data/download.zip
 extracting: example_data/dest_11035_air_temperature.nc  
 extracting: example_data/dest_11035_relative_humidity.nc  
```

or as a python request:

```python
import requests, zipfile, io, os, time
t0 = time.time()
# 
r = requests.post('http://early-upper-air.copernicus-climate.eu',
                  headers={'content-type': 'application/json'},
                  json={"statid": "11035",
                        "date": [20000101, 20000101],
                        "pressure_level": [1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000, 40000,
                                           50000, 70000, 85000, 92500, 100000],
                        "variable": ["temperature", "relative_humidity"],
                        "fbstats": ["obs_minus_bg", "obs_minus_an", "bias_estimate"]
                        },
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
    print("Unzipping retrieved files: to ./exampla_data")
    print(z.namelist())
    z.extractall(path='./example_data')
    z.close()
except:
    print(r.text)
    print("Error in request")
print("Time elapsed: ", time.time()-t0, "s")
```

# License

Generated using Copernicus Climate Change Service Information, 2020
[Copernicus Climate Change Service (C3S), 2020](https://apps.ecmwf.int/datasets/licences/copernicus/)

