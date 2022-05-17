# Sensor Metadata Inventory

This repository contains the scripts and the metadata sources used to build the inventory of the radiosonde metadata for the Comprehensive Upper Air Network (CUON) dataset.

This inventory represents an extensive collection of metadata describing the sensors used in radiosonde experiments. 


## Metadata Sources
In this Section we briefly describe the origin of the metadata consulted for the identification of the sensors.

Three main metadata sources were considered to build the inventory:

- the Schroeder metadata from xxx.

- the IGRA2 historical metadata at https://www.ncei.noaa.gov/pub/data/igra/history/ 

- WMO sensor identifiers contained in the ODB files from ERA5 reanalysis. The codes are extracted from the table at page A-398 of https://library.wmo.int/doc_num.php?explnum_id=10235 

### WMO Metadata

Describe problem with duplicated codes for different sensors.

### Schoreder Metadata



### IGRA2 Metadata

Reading the file from https://www.ncei.noaa.gov/pub/data/igra/history/igra2-metadata.txt

See https://www.ncei.noaa.gov/pub/data/igra/history/igra2-metadata-readme.txt for explanations.

The IGRA metadata consist of three principal types of radiosonde station
history information.

1. The primary source is Gaffen (1996):  Records from this source are 
   identified by an update date of 99/1996.

2. Updates to records for GUAN and RATPAC stations (Free et al. 2005) for 
   the period since Gaffen (1996):  Most of these records are identified 
   by a source attribution of Hammer (personal communication) and are 
   based on information received since 2004 through communication between 
   the NCDC GCOS Lead Center and GCOS focal points at Hydrological and 
   Meteorological Service Centers (HMSCs) in WMO member countries. 
   Additional sources of updated information include Joyce (personal 
   communication) and NOAA/NWS.

3. Events supplied by Steve Schroeder of Texas A&M University: 
   Identified by a source attribution of Schroeder (2007), these events 
   cover stations in the Russian Federation, India, Japan, China 
   (including Hong Kong and Taiwan), and Antarctic stations operated by 
   those countries.

Note that often, only the year of the IGRA2 metadata is available, and not the month. By default, when the month is missing, it is set to January.


### Combining Metadata

Sensor data were extracted based on the primary identification number of the station data. Most of the stations possess a proper WIGOS identifier, which was built in the form 0-2000x-0-yyyyy, where x is a numerical value while y can be general (numerical or literal).

WIGOS ids of this type are based on WMO ids, and data of such identifiers can be extracted from both the IGRA2 and the Schroeders inventory.

As for the WMO metadata, they do not rely on WMO ids but data is contained inside the ODB files that can be downloaded from the MARS archive XXXX . Frequently the stations are identified by numerical code which do not correspond to WIGOS/WMO ids, hence at times only WMO data is available.


## Workflow
The idea of this project is to extract all the available information concenring the sensor identification for a given radiosonde station, from the three available sources. Once this is available, we also provide a graphical interface to visualize the data.

The extraction of the metadata involves several steps.

- WMO data are extracted from the 'era5fb' table which is contained in the merged files constituting the CUON dataset.

- Schroeders identifiers are also extracted during the merging phase, in the following way.
All the relevant data for a station with mathcing WMO identifier is extracted; data usually comprise a date and a sensor id, where the date corresponds to the time when the sensor was replaced. As an approximation, we extend backward the information of the first sensor recorded, meaning that if data from a given station is recorded before the date of the first sensor recorded in Schroeder's data, we assume that the same sensor had been used from the beginning. Furthermore, we assume that the same sensor is being used until a new sensor replacement was recorded.

- Schoreder's data is assumed to be the most accurate metadata available, hence to be preferred with respect to WMO data, until the year 2013. This means that, when assigning a sensor identifier for a specific date, if Schroeder's data was available (before 2013), then Schroeder's id was stored. After 2013, only WMO data is considered.


- Valuable information can also be retireved in the IGRA2 archives, although not as extensive as Schoreders'


## Data Visualization

Sensor Metadata can be conveniently visualized in different ways.
The charts are produced using the \verb{Plotly} library, which can be embedded in Jupyter notebooks.

Charts can be saved as "html"" files, that can be visualized in common web browsing software to take advantage of the dashboard interactive functionalities, as well as static png images.

Two types of charts are realized: plain tables and time series.

Tables summarises the content of all the different sensors that can be retrieved in the three sources (WMO, Schroeders' and IGRA2). No time information regarding the introduction of the sensor is stored here, so the order of appearance of the sensor in the table is not meaningful. The description of the sensor is also provided.

Time series gives the complete picture of the evolution of the usage of different sensor in the station. Three different markers are used to indicate when sensor events (i.e. introduction of a new sensor) happened as recorded in the three sources.

The following points are worth some explanation.
- WMO data is very extensive since it is provided for every record available in the original ERA5 data, eventually stored as NULL values. These have been removed. A high density of points can be explained since, frequently, different sonde models were used in different times of the days, so every change in the sensor is being tracked. This applies also to Schroeder's data but in a more limited way.

IGRA2 data is more limited, and it is represented both by markers and by vertical lines.

Finally, the solid green line represent the Standard Normal Homogeneity Test statistics which is used to determine break-points in time series and calculate bias adjustments. 


### Interactive Dashboard

An interactive dashboard can be launched via a python file called *sensor_metadata_DASHBOARD.py* or, alternatively, by a jupyter notebook called *sensor_metadata_DASHBOARD.ipynb*.

In particular we suggest to download the directory called "example" since it contains a self-contained working example. The scripts were tested using a dedicated anaconda environment based on python 3.8. The environment can be reproduced by installing the dependencies described in the file "requirements_Dash.txt".

The example can be tested only on four sample stations, namely 

- Lindenberg (0-20001-0-10393)

- Payerne (0-20000-06610)

- Vienna (0-20001-0-11035)

- Cachimbo (0-20000-0-82930)

By default, Cachimbo is loaded. Stations can be selected by typing the WIGOS identifier in the cell on top left of the dashboard. For now, if a wrong WIGOS is provided, the software will break; in this case, please reload the dashboard.


A third chart is also provided, which gives the total count of the WMO sensor identifiers. Many stations report unreasonable sensor models for few records, which most likely are due to erroneous data processing or data recording. 


