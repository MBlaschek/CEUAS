# Data harvesting

It is explained here how to download Radiosonde data from different sources.

Type | Nr | Responsible | Nature | Title| Due | Status | File
---|---|---|---|---|---|---|---
Deliverable | DC3S311c_Lot2.1.2.1 | UNIVIE | Software, Report | C3S Upper air data harvest and upload toolbox v0|  Nov 2019 | delayed | code/ 

Table of Content:
- [Download radiosonde data from different sources](#download-radiosonde-data-from-different-sources)
  * [Requirements](#requirements)
  * [IGRAv2](#igrav2)
  * [UADB](#uadb)
  * [MARS ODB](#MARS-ODB)
- [Convert to Netcdf](#Convert-to-CDM-complaint-NetCDF)
- [License](#License)

# Tasks

In order to build the radiosonde archive follow these steps:
* Download radiosonde data from different sources
  * IGRAv2 (not credentials needed)
    * `code/download/download_IGRAv2.sh`
  * UADB (credentials needed [NCAR](https://rda.ncar.edu)
    * `code/download/download_UADB.csh`
  * MARS - ERA5, ECMWF (credentials needed)
* Read / convert data to NetCDF4
* Produce Archive Overview

# Download radiosonde data from different sources
In order to download the data from all sources (**Attention: credentials are needed for some**) follow the specific information here or in the Deliverable Document. In the Chapters below on can find the individual download procedures per data source. In order to work, it is highly recommended to use a Linux environment, a stable and fast internet connection and have enough disk space available. 

## Requirements
Please check the following before you execute:
* Each script has the needed credentials (modify script)
* Tools (wget, csh, bash, odb, mars)

## IGRAv2
The Integrated Global Radiosonde Archive ([IGRA](https://doi.org/10.7289/V5X63K0Q)) version 2 was released in 2016 and updated since then. There are about 2700 stations in that archive covering the time 1905 to present. Please find the description at the data site or in the `pre_merge_stations.py` routine.

In order to download the whole IGRA v2 Archive execute the following script. Note if you want to change the default download location, please change the `DATADIR` variable inside the script:
```bash
bash download_IGRAv2.sh
```
Example output, downloading the zipped ascii files:
```bash
--2020-01-17 13:23:15--  https://www1.ncdc.noaa.gov/pub/data/igra/data/data-por/ACM00078861-data.txt.zip
Verbindungsaufbau zu www1.ncdc.noaa.gov|205.167.25.178|:443... verbunden.
HTTP Anforderung gesendet, warte auf Antwort... 200 OK
Länge: 6731308 (6,4M) [application/zip]
In »»ACM00078861-data.txt.zip«« speichern.

100%[======================================>] 6.731.308   3,12M/s   in 2,1s    

2020-01-17 13:23:18 (3,12 MB/s) - »»ACM00078861-data.txt.zip«« gespeichert [6731308/6731308]

--2020-01-17 13:23:18--  https://www1.ncdc.noaa.gov/pub/data/igra/data/data-por/AEM00041217-data.txt.zip
Verbindungsaufbau zu www1.ncdc.noaa.gov|205.167.25.178|:443... verbunden.
HTTP Anforderung gesendet, warte auf Antwort... 200 OK
Länge: 22943202 (22M) [application/zip]
In »»AEM00041217-data.txt.zip«« speichern.

100%[======================================>] 22.943.202  7,62M/s   in 2,9s    

2020-01-17 13:23:21 (7,62 MB/s) - »»AEM00041217-data.txt.zip«« gespeichert [22943202/22943202]

--2020-01-17 13:23:21--  https://www1.ncdc.noaa.gov/pub/data/igra/data/data-por/AEXUAE05467-data.txt.zip
Verbindungsaufbau zu www1.ncdc.noaa.gov|205.167.25.178|:443... verbunden.
HTTP Anforderung gesendet, warte auf Antwort... 200 OK
Länge: 105036 (103K) [application/zip]
In »»AEXUAE05467-data.txt.zip«« speichern.

100%[======================================>] 105.036      300K/s   in 0,3s  
...
```

The script creates the following structure and source files:
```bash
ll /tmp/data/igra | head
total 22G
-rw-r--r--. 1 user users  6,5M  8. Dez 13:48 ACM00078861-data.txt.zip
-rw-r--r--. 1 user users   22M  8. Dez 13:48 AEM00041217-data.txt.zip
-rw-r--r--. 1 user users  103K  8. Dez 13:48 AEXUAE05467-data.txt.zip
```
## UADB
The NCAR Upper Air Database ([UADB](https://rda.ucar.edu/datasets/ds370.1/)) was released in 2014 and updated since then. There are about 1100 stations in that archive covering the time 1920 to present. Please find the description at the data site or in the `pre_merge_stations.py` routine.

In order to download the whole UADB Archive execute the following script. The script requires credentials (username (email) and password) from the user and the file `UADB.files.list`. The needed information can be found on [rda.ucar.edu](https://rda.ucar.edu/index.html?hash=data_user&action=register). Note if you want to change the default download location, please change the `DATADIR` variable inside the script:
```bash
csh download_UADB.csh [USERNAME] [PASSWORD]
```
Example output, downloading the ascii files:
```bash
wget --no-check-certificate -N --load-cookies auth.rda_ucar_edu http://rda.ucar.edu/data/ds370.1/uadb_trhc_10185.txt ...
--2020-01-21 15:32:27--  http://rda.ucar.edu/data/ds370.1/uadb_trhc_10185.txt
Auflösen des Hostnamen »rda.ucar.edu«.... 128.117.181.113
Verbindungsaufbau zu rda.ucar.edu|128.117.181.113|:80... verbunden.
HTTP Anforderung gesendet, warte auf Antwort... 301 Moved Permanently
Platz: https://rda.ucar.edu/data/ds370.1/uadb_trhc_10185.txt[folge]
--2020-01-21 15:32:27--  https://rda.ucar.edu/data/ds370.1/uadb_trhc_10185.txt
Verbindungsaufbau zu rda.ucar.edu|128.117.181.113|:443... verbunden.
HTTP Anforderung gesendet, warte auf Antwort... 200 OK
Länge: 1834826 (1,7M) [text/plain]
In »»uadb_trhc_10185.txt«« speichern.

100%[======================================>] 1.834.826   1,22M/s   in 1,4s    

2020-01-21 15:32:29 (1,22 MB/s) - »»uadb_trhc_10185.txt«« gespeichert [1834826/1834826]

wget --no-check-certificate -N --load-cookies auth.rda_ucar_edu http://rda.ucar.edu/data/ds370.1/uadb_trhc_1020.txt ...
--2020-01-21 15:32:29--  http://rda.ucar.edu/data/ds370.1/uadb_trhc_1020.txt
Auflösen des Hostnamen »rda.ucar.edu«.... 128.117.181.113
Verbindungsaufbau zu rda.ucar.edu|128.117.181.113|:80... verbunden.
HTTP Anforderung gesendet, warte auf Antwort... 301 Moved Permanently
Platz: https://rda.ucar.edu/data/ds370.1/uadb_trhc_1020.txt[folge]
--2020-01-21 15:32:30--  https://rda.ucar.edu/data/ds370.1/uadb_trhc_1020.txt
Verbindungsaufbau zu rda.ucar.edu|128.117.181.113|:443... verbunden.
HTTP Anforderung gesendet, warte auf Antwort... 200 OK
Länge: 260767 (255K) [text/plain]
In »»uadb_trhc_1020.txt«« speichern.

100%[======================================>] 260.767      442K/s   in 0,6s    

2020-01-21 15:32:31 (442 KB/s) - »»uadb_trhc_1020.txt«« gespeichert [260767/260767]
...
```
The script creates the following structure and source files:
```bash
ll /tmp/data/ncar | head
insgesamt 175G
-rw-r--r--. 1 user users   21M  1. Mär 2019  uadb_trhc_1001.txt
-rw-r--r--. 1 user users  249K  1. Mär 2019  uadb_trhc_10034.txt
-rw-r--r--. 1 user users  168M 15. Mai 2019  uadb_windc_10035.txt
-rw-r--r--. 1 user users   22M 15. Mai 2019  uadb_windc_10046.txt
```
## MARS ODB 

The MARS archive is a huge data collection including satellite and in-situ observations from around the world and resides at the ECMWF. The Observation Data Base (ODB) can be retrieved from the ECMWF assuming valid credentials. The Copernicus Early Upper Air archive consists of three sources available from ECMWF; these are

1. The ERA-CLIM(2) digitized data in ODB format (experiment identifier 3188)

2. The NCAR upper air data sets in ODB format (experiment identifiers 1759 and 1761)

3. The ERA5 analysis feedback archive

The ODB files can be retrieved with the following MARS (https://www.ecmwf.int/en/forecasts/datasets/archive-datasets) jobs:

 ```bash
retrieve, 
	type=an, 
	class=ea, 
	expver=${EXP}, 
	stream=oper, 
	date=${YYY}${MMM}01/to/${YYY}${MMM}${DDD}, 
	time=00/12, 
	decade=${DEC}, 
	reportype=16013/16022/16045/16068, 
	type=ofb,
	target='era5.conv.${YYY}${MMM}'
 ```

`EXP` can be `3188 (CHUAN)`, `1759 (NCAR)` or `1761 (NCAR)`. If one wants to retrieve the ERA5 feedback data, one can use experiment numbers `3645` (1958-), `3647` (1968-), `3649` (1972-) and `3651` (1950-) as well as `1` (1979-). The monthly ODB files must be rearranged into station series, using ODB and cat commands. A scripts for this can be found  in `code/download/splite5_1.ksh`. This is a relatively time consuming process.

The script creates the following structure and source files:
```bash
ll /tmp/data/era5 | head
insgesamt 175G
-rw-r--r--. 1 user users   21M  1. Mär 2019  uadb_trhc_1001.txt
-rw-r--r--. 1 user users  249K  1. Mär 2019  uadb_trhc_10034.txt
-rw-r--r--. 1 user users  168M 15. Mai 2019  uadb_windc_10035.txt
-rw-r--r--. 1 user users   22M 15. Mai 2019  uadb_windc_10046.txt
```

The ODB files need to be converted to ASCII:

```bash
odb sql -q 'select *' -i ${file} | tr -d " "  | gzip  > ${file}.gz
```

The gzipped ASCII files can be read by `pandas.read_csv`.

# Convert to CDM compliant NetCDF

The script `harvest_convert_to_netCDF.py` reads the data from the different data sources and converts them into a CDM complaint netCDF.

???



The final output looks like this:

???



# License

