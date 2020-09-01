# Copernicus Climate Change - Early Upper Air Service (CEUAS)

> Service Contract C3S 311c Lot2

Collaborators:

* L. Haimberger
* M. Blaschek
* F. Ambrogi
* U. Voggenberger



Status: Development
Date: 01-2020

Host: UNIVIE



Contents:
* [cdm](#CDM)
* [meta](#Metadata)
* [public/cds-backend](#CDS-Backend)
* [public/harvest](#Harvest)
* [public/merge](#Merge)
* [public/adjust](#Adjust)
* [public/uncertainties](#Uncertainties)

# CDM

Information can be found in `CEUAS/cdm`.

Downloads and create the CDM tables, ready as input for the netCDF files.

# Metadata

Information can be found in `CEUAS/meta`.

This generates the station inventory from multiple sources and produces a consolidated list.

# CDS Backend

Information can be found in `CEUAS/public/cds-backend`.

This runs the backend server for the CDS. An example on how to use can be found here: [Example.ipynb](https://github.com/MBlaschek/CEUAS/blob/master/CEUAS/public/cds-backend/Example.ipynb)

# Harvest

Information can be found in `CEUAS/public/harvest`.

The scripts there allow to download the source datasets and produces CDM compliant netCDF files.

# Merge

Information can be found in `CEUAS/public/merge`.

The scripts there allow to merge the CDM compliant netCDF files from the Harvest step and produce a combined dataset.

# Adjust

Information can be found in `CEUAS/public/adjust`.

The scripts there allow to run the humidity homogenization on the CDM compliant files delivered via the CDS backend. An example can be found here: [Example Homogenization.ipynb](https://github.com/MBlaschek/CEUAS/blob/master/CEUAS/public/adjust/Example_Homogenization.ipynb)

# Uncertainties

Information can be found in `CEUAS/public/uncertainties`.

The scripts there allow to calculate estimates of upper air observation errors from Desrozier's method using background departures.

