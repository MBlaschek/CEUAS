#!/bin/bash
#################################################################
# Bash script to retrieve all IGRAv2 files (zip)
# Uses wget
# No credentials are needed
# 
# Specify DATADIR according to your system
#################################################################
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
