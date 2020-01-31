#!/bin/csh
# -----------------------------------------------------------------------------
# Csh Script to retrieve all online Data file of 'ds370.1',
# uadb-trhc: 82.42 GB, uadb-windc: 92.52 GB
# total 174.94G. The Script needs a filelist (UADB.files.list)
# This script uses 'wget' to download data and
# tar to extract the archive.
#
# Credential are needed, register at rda.ucar.edu (free)
#
# Call with Email Adress and Password to start download
# Specify DATADIR according to your system
#
# (c) University of Vienna, M. Blaschek, Vienna, Austria
# Copernicus Climate Change Service, 2020
# https://apps.ecmwf.int/datasets/licences/copernicus/
# email michael.blaschek (at) univie.ac.at
# Created: Vienna, 15 August, 2017
# Last Accessed: 15 January, 2020
# -----------------------------------------------------------------------------

set userid = $1
set pswd = $2

if(x$pswd == x && `env | grep RDAPSWD` != '') then
 set pswd = $RDAPSWD
endif

if (x$pswd == x) | (x$userid == x) then
 echo
 echo Usage: $0 YourUsername YourPassword
 echo
 exit 1
endif

set v = `wget -V |grep 'GNU Wget ' | cut -d ' ' -f 3`
set a = `echo $v | cut -d '.' -f 1`
set b = `echo $v | cut -d '.' -f 2`
if(100 * $a + $b > 109) then
 set opt = 'wget --no-check-certificate'
else
 set opt = 'wget'
endif

set opt1 = '-O Authentication.log --save-cookies auth.rda_ucar_edu --post-data'
set opt2 = "email=$userid&passwd=$pswd&action=login"
$opt $opt1="$opt2" https://rda.ucar.edu/cgi-bin/login
set opt1 = "-N --load-cookies auth.rda_ucar_edu"

set opt2 = "$opt $opt1 http://rda.ucar.edu/data/ds370.1/"
#
# read filelist from file
#
set filelist = `cat UADB.files.list`
#
# Datadir
#
set DATADIR=/tmp/data/ncar
mkdir -vp $DATADIR
cd $DATADIR
#
# DOWNLOAD
#
while($#filelist > 0)
 set syscmd = "$opt2$filelist[1]"
 echo "$syscmd ..."
 $syscmd
 shift filelist
end
#
# change back
#
cd ..
#
#
#
rm -f auth.rda_ucar_edu Authentication.log
exit 0
