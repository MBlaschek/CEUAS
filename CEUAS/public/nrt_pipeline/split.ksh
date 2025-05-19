#!/bin/ksh
#
# This script takes all monthly odb files in a directory, splits them 
# according to statid@hdr and then concatenates the split files into 
# time series of observations with the same statid@hdr.
#
# While the script is fairly generic, the expids and  
#
# maximum open file limit should be set generously by OS (>=5000)
#
# Leopold Haimberger, 23 January 2020
#

ulimit -a

echo 'NOTE NOTE NOTE NOTE'
echo 'NOTE: script filter_odbgz.ksh must be executed afterwards'
echo 'NOTE NOTE NOTE NOTE'

pdir=PDIR
n=0
#for dir in 1759 1761 3188 ; do   # original input from NCAR/CHUAN
#for dir in 3647 3649 3651 ; do # expids for ERA5.1 and the pre-1979 streams
#for dir in 13651 ; do # expids for ERA5.1 and the pre-1979 streams
#for dir in 2929 ; do # expids for ERA5.1 and the pre-1979 streams
for dir in 1 ; do   # main expid for ERA5
  cd $pdir
#  \rm *.conv*._*
  for file in `ls *.conv.YYYYMM` ; do
     echo $file

     n=$(($n+1))
     if [ $n -lt 12 ] ; then
        odc split -no_verification -maxopenfiles 5000 $file "$file.{statid@hdr}"
     else
	      odc split -no_verification -maxopenfiles 5000 $file "$file.{statid@hdr}" 
        n=0
     fi


  done

  exit

  nproc=`ps -ef | grep 'odb split' | wc -l`
  while [ $nproc -gt 1 ] ; do
    sleep 10
    nproc=`ps -ef | grep 'odb split' | wc -l`
  done
 
  yy=1986
  n=0
  while [ $yy -lt 1987 ] ; do
    for file in `ls era5.conv.${yy}??.* 2>/dev/null` ; do

      fil=`echo $file | cut -c 1-9`
      suff=`echo $file |  cut -c 18-25`
      echo $yy,$fil $suff, $file
    
      if [ $n -lt 20 ] ; then
#        cat $file >> $fil._$suff &
        n=$(($n+1))
      else
#        cat $file >> $fil._$suff 
        n=0
      fi
 

      #ls -l  $fil._$suff
    done
    echo $yy
    yy=$(($yy+1))
  done
done

touch split_done.txt
