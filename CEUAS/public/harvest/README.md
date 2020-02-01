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
ls -lh /tmp/data/igra | head
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
ls -lh /tmp/data/ncar | head
total 175G
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

The ODB files need to be converted to ASCII:

```bash
odb sql -q 'select *' -i ${file} | tr -d " "  | gzip  > ${file}.gz
```

The script creates the following structure and source files:
```bash
ls -lh /tmp/data/era5/odbs/1 | head
total 935G
-rw-r--r--. 1 user users  4,9G 11. Jän 17:27 era5.conv._01001.gz
-rw-r--r--. 1 user users  1,7K 11. Jän 16:00 era5.conv._01002.gz
-rw-r--r--. 1 user users  5,8K 11. Jän 16:00 era5.conv._01003.gz
-rw-r--r--. 1 uesr users  3,3G 11. Jän 16:56 era5.conv._01004.gz
-rw-r--r--. 1 uesr users  5,3K 11. Jän 16:00 era5.conv._01006.gz
```

The gzipped ASCII files can be read by `pandas.read_csv`.

# Convert to CDM compliant NetCDF

The script `harvest_convert_to_netCDF.py` reads the data from the different data sources and converts them into a CDM complaint netCDF.

In order to test the harvesting, one would of course have to retrieve all the original data from ECMWF. MARS requests for retrieving them are available in the software repository as well. To execute them, an ECMWF member state user ID is currently required. The conversion into netCDF4 is done with the script `harvest_to_netCDF_converter.py`, which is called by the run script `run_converter.py`.

To test it, it is best to clone the whole repository and to run the  `run_converter.py` script there. The necessary test data are fetched from the Internet or are provided in the repository.

For running the scripts we recommend a Linux environment with the [Anaconda Python 3](https://www.anaconda.com/distribution/) distribution installed. It contains the necessary packages (such as `xarray, pandas, netCDF4, h5py`). 

Additionally, [eccodes](([https://confluence.ecmwf.int//display/ECC/ecCodes+Home](https://confluence.ecmwf.int/display/ECC/ecCodes+Home)) and the ECMWF [ODB API](https://confluence.ecmwf.int/display/ODBAPI/ODB+API+Home)  must be installed on the machine to be able to decode the BUFR and ODB data bases. 

The final output looks like this:

```bash
ls -lh /tmp/data/
drwxr-xr-x. 124 user users  23M 31. Jän 10:58 era5_1
-rw-r--r--. 1 user users  275M  2. Jun 2019  era5_1/ccera5.conv._01001.nc
-rw-r--r--. 1 user users  172K  2. Jun 2019  era5_1/ccera5.conv._01002.nc
-rw-r--r--. 1 user users  176K  2. Jun 2019  era5_1/ccera5.conv._01003.nc
drwxr-xr-x.   3 user users 564K 27. Jän 17:17 era5_1759
-rw-r--r--. 1 user users  191M 19. Jun 2019  era5_1759/era5.1759.conv.1:1001.nc
-rw-r--r--. 1 user users   20M 19. Jun 2019  era5_1759/era5.1759.conv.1:10034.nc
-rw-r--r--. 1 user users  235M 19. Jun 2019  era5_1759/era5.1759.conv.1:10035.nc
drwxr-xr-x.   4 user users 1,1M 27. Jän 22:37 1761
drwxr-xr-x.   6 user users 152K 27. Jän 13:22 3188
drwxr-xr-x.   3 user users 8,1M 23. Jän 23:18 3645
drwxr-xr-x.   3 user users 3,7M 12. Jän 09:41 3647
drwxr-xr-x.   3 user users 8,0M 12. Jän 09:55 3649
drwxr-xr-x.   3 user users 8,5M 12. Jän 10:32 3651
drwxr-xr-x.   2 user users 4,0K  4. Okt 10:55 igra
-rw-r--r--. 1 user users   15M 29. Jän 20:02 chGMXUAC00001-data.txt.nc
-rw-r--r--. 1 user users  2,9M 29. Jän 20:02 chGMXUAC00002-data.txt.nc
-rw-r--r--. 1 user users  5,0M 29. Jän 20:02 chGMXUAC00006-data.txt.nc
drwxr-xr-x.   2 user users 4,0K  4. Okt 10:55 ncar
-rw-r--r--. 1 user users  3,0G 30. Jän 12:32 UADB_trhc_010035.nc
-rw-r--r--. 1 user users  388M 30. Jän 10:07 UADB_trhc_010046.nc
-rw-r--r--. 1 user users  719M 30. Jän 12:33 UADB_trhc_010113.nc
```

### NCDUMP of example CDM file

An example file looks like this:

```bash
netcdf \0-20000-0-01001_era5.conv._01001 {
dimensions:
	record = 33688 ;
	days = 3 ;
	drange = 14906 ;
variables:
	int dateindex(days, drange) ;
	int recordindex(record) ;
	int64 recordtimestamp(record) ;
		string recordtimestamp:units = "seconds since 1900-01-01 00:00:00" ;

group: crs {
  dimensions:
  	crs_len = 1 ;
  	string80 = 80 ;
  variables:
  	int64 crs(crs_len) ;
  	char description(crs_len, string80) ;
  } // group crs

group: era5fb {
  dimensions:
  	index = 182438469 ;
  	string10 = 10 ;
  	string7 = 7 ;
  variables:
  	float an_depar@body(index) ;
  	int andate(index) ;
  	int antime(index) ;
  	float biascorr@body(index) ;
  	int date@hdr(index) ;
  	float fg_depar@body(index) ;
  	char index(index) ;
  	float lat@hdr(index) ;
  	float lon@hdr(index) ;
  	float obsvalue@body(index) ;
  	int reportype(index) ;
  	int seqno@hdr(index) ;
  	float sonde_type@conv(index) ;
  	char source@hdr(index, string10) ;
  	char statid@hdr(index, string7) ;
  	int time@hdr(index) ;
  	int varno@body(index) ;
  	float vertco_reference_1@body(index) ;
  } // group era5fb

group: header_table {
  dimensions:
  	hdrlen = 33688 ;
  	string5 = 5 ;
  	string80 = 80 ;
  	string7 = 7 ;
  variables:
  	double application_area(hdrlen) ;
  		string application_area:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	char comments(hdrlen, string80) ;
  		string comments:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	int crs(hdrlen) ;
  		string crs:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	int duplicate_status(hdrlen) ;
  		string duplicate_status:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	double duplicates(hdrlen) ;
  		string duplicates:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	double events_at_station(hdrlen) ;
  		string events_at_station:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	float height_of_station_above_local_ground(hdrlen) ;
  		height_of_station_above_local_ground:_FillValue = NaNf ;
  		string height_of_station_above_local_ground:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	float height_of_station_above_sea_level(hdrlen) ;
  		height_of_station_above_sea_level:_FillValue = NaNf ;
  		string height_of_station_above_sea_level:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	float height_of_station_above_sea_level_accuracy(hdrlen) ;
  		height_of_station_above_sea_level_accuracy:_FillValue = NaNf ;
  		string height_of_station_above_sea_level_accuracy:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	char history(hdrlen, string80) ;
  		string history:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	char instrument(hdrlen, string80) ;
  		string instrument:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	float latitude(hdrlen) ;
  		latitude:_FillValue = NaNf ;
  		string latitude:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	float location_accuracy(hdrlen) ;
  		location_accuracy:_FillValue = NaNf ;
  		string location_accuracy:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	int location_method(hdrlen) ;
  		string location_method:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	int location_quality(hdrlen) ;
  		string location_quality:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	float longitude(hdrlen) ;
  		longitude:_FillValue = NaNf ;
  		string longitude:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	int number_of_pressure_levels(hdrlen) ;
  		string number_of_pressure_levels:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	double observing_programme(hdrlen) ;
  		string observing_programme:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	char owner(hdrlen, string80) ;
  		string owner:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	int platform_sub_type(hdrlen) ;
  	int platform_type(hdrlen) ;
  	char primary_station_id(hdrlen, string7) ;
  		string primary_station_id:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	int primary_station_id_scheme(hdrlen) ;
  		string primary_station_id_scheme:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	double processing_codes(hdrlen) ;
  		string processing_codes:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	int processing_level(hdrlen) ;
  		string processing_level:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	char product_name(hdrlen, string80) ;
  		string product_name:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	int product_version(hdrlen) ;
  		string product_version:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	char profile_id(hdrlen, string80) ;
  		string profile_id:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	int64 record_timestamp(hdrlen) ;
  		string record_timestamp:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	char references(hdrlen, string80) ;
  		string references:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	int region(hdrlen) ;
  		string region:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	int report_duration(hdrlen) ;
  		string report_duration:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	char report_id(hdrlen, string5) ;
  		string report_id:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	int report_meaning_of_timestamp(hdrlen) ;
  		string report_meaning_of_timestamp:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	int report_quality(hdrlen) ;
  		string report_quality:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	float report_synoptic_time(hdrlen) ;
  		report_synoptic_time:_FillValue = NaNf ;
  		string report_synoptic_time:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	float report_time_accuracy(hdrlen) ;
  		report_time_accuracy:_FillValue = NaNf ;
  		string report_time_accuracy:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	int report_time_quality(hdrlen) ;
  		string report_time_quality:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	int report_time_reference(hdrlen) ;
  		string report_time_reference:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	int64 report_timestamp(hdrlen) ;
  		string report_timestamp:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	int report_type(hdrlen) ;
  		string report_type:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	int sea_level_datum(hdrlen) ;
  		string sea_level_datum:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	char source_id(hdrlen, string80) ;
  		string source_id:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	char source_record_id(hdrlen, string80) ;
  		string source_record_id:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	float station_course(hdrlen) ;
  		station_course:_FillValue = NaNf ;
  		string station_course:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	float station_heading(hdrlen) ;
  		station_heading:_FillValue = NaNf ;
  		string station_heading:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	char station_name(hdrlen, string80) ;
  	int station_record_number(hdrlen) ;
  		string station_record_number:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	float station_speed(hdrlen) ;
  		station_speed:_FillValue = NaNf ;
  		string station_speed:coordinates = "platform_type platform_sub_type station_name station_type" ;
  	int station_type(hdrlen) ;
  	int sub_region(hdrlen) ;
  		string sub_region:coordinates = "platform_type platform_sub_type station_name station_type" ;
  } // group header_table

group: id_scheme {
  dimensions:
  	id_scheme_len = 7 ;
  	string80 = 80 ;
  variables:
  	char description(id_scheme_len, string80) ;
  	int64 scheme(id_scheme_len) ;
  } // group id_scheme

group: observations_table {
  dimensions:
  	index = 182438469 ;
  	string5 = 5 ;
  	string8 = 8 ;
  	string80 = 80 ;
  variables:
  	char adjustment_id(index, string80) ;
  	int advanced_assimilation_feedback(index) ;
  	int advanced_homogenisation(index) ;
  	int advanced_qc(index) ;
  	int advanced_uncertainty(index) ;
  	float bbox_max_latitude(index) ;
  	float bbox_max_longitude(index) ;
  	float bbox_min_latitude(index) ;
  	float bbox_min_longitude(index) ;
  	int code_table(index) ;
  	int conversion_flag(index) ;
  	int conversion_method(index) ;
  	int crs(index) ;
  	int data_policy_licence(index) ;
  	int64 date_time(index) ;
  		date_time:units = "seconds since 1900-01-01 00:00:00" ;
  	int date_time_meaning(index) ;
  	int exposure_of_sensor(index) ;
  	char index(index) ;
  	float latitude(index) ;
  	int location_method(index) ;
  	float location_precision(index) ;
  	float longitude(index) ;
  	float numerical_precision(index) ;
  	int observation_duration(index) ;
  	float observation_height_above_station_surface(index) ;
  	char observation_id(index, string8) ;
  	float observation_value(index) ;
  	int64 observed_variable(index) ;
  	int original_code_table(index) ;
  	float original_precision(index) ;
  	int original_units(index) ;
  	float original_value(index) ;
  	int processing_level(index) ;
  	int quality_flag(index) ;
  	char report_id(index, string5) ;
  	int secondary_value(index) ;
  	int secondary_variable(index) ;
  	int sensor_automation_status(index) ;
  	char sensor_id(index, string80) ;
  	char source_id(index, string80) ;
  	int spatial_representativeness(index) ;
  	int traceability(index) ;
  	int units(index) ;
  	int value_significance(index) ;
  	float z_coordinate(index) ;
  	int z_coordinate_method(index) ;
  	int z_coordinate_type(index) ;
  } // group observations_table

group: observed_variable {
  dimensions:
  	observed_variable_len = 120 ;
  	string80 = 80 ;
  variables:
  	char description(observed_variable_len, string80) ;
  	char domain(observed_variable_len, string80) ;
  	char name(observed_variable_len, string80) ;
  	char parameter_group(observed_variable_len, string80) ;
  	char sub_domain(observed_variable_len, string80) ;
  	char units(observed_variable_len, string80) ;
  	char variable(observed_variable_len, string80) ;
  } // group observed_variable

group: source_configuration {
  } // group source_configuration

group: station_configuration {
  dimensions:
  	station_configuration_len = 1 ;
  	string15 = 15 ;
  	string1 = 1 ;
  	string2 = 2 ;
  	string6 = 6 ;
  	string18 = 18 ;
  	string10 = 10 ;
  	string3 = 3 ;
  	string8 = 8 ;
  	string24 = 24 ;
  	string12 = 12 ;
  variables:
  	char alternative_name(station_configuration_len, string2) ;
  	double bbox_max_latitude(station_configuration_len) ;
  		bbox_max_latitude:_FillValue = NaN ;
  	double bbox_max_longitude(station_configuration_len) ;
  		bbox_max_longitude:_FillValue = NaN ;
  	double bbox_min_latitude(station_configuration_len) ;
  		bbox_min_latitude:_FillValue = NaN ;
  	double bbox_min_longitude(station_configuration_len) ;
  		bbox_min_longitude:_FillValue = NaN ;
  	char city(station_configuration_len, string8) ;
  	char comment(station_configuration_len, string2) ;
  	char contact(station_configuration_len, string2) ;
  	char end_date(station_configuration_len, string10) ;
  	double latitude(station_configuration_len) ;
  		latitude:_FillValue = NaN ;
  	char local_gravity(station_configuration_len, string2) ;
  	double longitude(station_configuration_len) ;
  		longitude:_FillValue = NaN ;
  	char measuring_system_id(station_configuration_len, string2) ;
  	char measuring_system_model(station_configuration_len, string2) ;
  	char metadata_contact(station_configuration_len, string12) ;
  	int64 metadata_contact_role(station_configuration_len) ;
  	char observed_variables(station_configuration_len, string24) ;
  	int64 observing_frequency(station_configuration_len) ;
  	char operating_institute(station_configuration_len, string3) ;
  	char operating_territory(station_configuration_len, string2) ;
  	char optional_data(station_configuration_len, string2) ;
  	int64 platform_sub_type(station_configuration_len) ;
  	int64 platform_type(station_configuration_len) ;
  	char primary_id(station_configuration_len, string15) ;
  	char primary_id_scheme(station_configuration_len, string1) ;
  	char record_number(station_configuration_len, string2) ;
  	int64 reporting_time(station_configuration_len) ;
  	char role(station_configuration_len, string2) ;
  	char secondary_id(station_configuration_len, string6) ;
  	char secondary_id_scheme(station_configuration_len, string2) ;
  	char start_date(station_configuration_len, string10) ;
  	char station_abbreviation(station_configuration_len, string2) ;
  	int64 station_automation(station_configuration_len) ;
  	int64 station_crs(station_configuration_len) ;
  	char station_name(station_configuration_len, string18) ;
  	int64 station_type(station_configuration_len) ;
  	int64 telecommunication_method(station_configuration_len) ;
  } // group station_configuration

group: station_configuration_codes {
  dimensions:
  	station_configuration_codes_len = 44 ;
  	string80 = 80 ;
  variables:
  	char abbreviation(station_configuration_codes_len, string80) ;
  	char code_value(station_configuration_codes_len, string80) ;
  	char description(station_configuration_codes_len, string80) ;
  	char field_id(station_configuration_codes_len, string80) ;
  	char field_name(station_configuration_codes_len, string80) ;
  } // group station_configuration_codes

group: station_type {
  dimensions:
  	station_type_len = 2 ;
  	string80 = 80 ;
  variables:
  	char description(station_type_len, string80) ;
  	int64 type(station_type_len) ;
  } // group station_type

group: units {
  dimensions:
  	units_len = 139 ;
  	string80 = 80 ;
  variables:
  	char abbreviation(units_len, string80) ;
  	char base_units(units_len, string80) ;
  	char name(units_len, string80) ;
  	char units(units_len, string80) ;
  } // group units

group: z_coordinate_type {
  dimensions:
  	z_coordinate_type_len = 2 ;
  	string80 = 80 ;
  variables:
  	char description(z_coordinate_type_len, string80) ;
  	int64 type(z_coordinate_type_len) ;
  } // group z_coordinate_type
}
```



# License

Durre, Imke; Xungang, Yin; Vose, Russell S.; Applequist, Scott; Arnfield, Jeff. (2016) Integrated Global Radiosonde Archive (IGRA), Version 2. NOAA National Centers for Environmental Information. DOI:10.7289/V5X63K0Q [15.12.2019]. This work is licensed under [U.S. Goverment Work](http://www.usa.gov/publicdomain/label/1.0/)

Research Data Archive/Computational and Information Systems Laboratory/National Center for Atmospheric Research/University Corporation for Atmospheric Research, 2014: NCAR Upper Air Database, 1920-ongoing. Research Data Archive at the National Center for Atmospheric Research, Computational and Information Systems Laboratory, Boulder, CO. [Available online at http://rda.ucar.edu/datasets/ds370.1/.] Accessed† 15 05 2019. This work is licensed under a [Creative Commons Attribution 4.0 International License](http://creativecommons.org/licenses/by/4.0/).

Generated using Copernicus Climate Change Service Information, 2020
[Copernicus Climate Change Service (C3S), 2020](https://apps.ecmwf.int/datasets/licences/copernicus/)