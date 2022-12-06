#!/bin/bash
# -----------------------------------------------------------------------------
# Bash script to retrieve all IGRAv2 files (zip)
# Uses wget
# No credentials are needed
#
# Specify DATADIR according to your system
#
# (c) University of Vienna, M. Blaschek, Vienna, Austria
# Copernicus Climate Change Service, 2020
# https://apps.ecmwf.int/datasets/licences/copernicus/
# email michael.blaschek (at) univie.ac.at
# Created: Vienna, 15 August, 2017
# Last Accessed: 15 January, 2020
# -----------------------------------------------------------------------------
# Data directory
#
DATADIR=/tmp/data/igra
mkdir -p $DATADIR
cd $DATADIR
#
# DOWNLOAD from NCDC / NOAA
#
wget -r -np -nd -R "index.html*" -e robots=off https://www1.ncdc.noaa.gov/pub/data/igra/data/data-por
cd -
