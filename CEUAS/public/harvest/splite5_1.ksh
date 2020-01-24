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

pdir=`pwd`
n=0
for dir in 1759 1761 3188 ; do   # original input from NCAR/CHUAN
#for dir in 2929 3645 3647 3649 3651 # expids for ERA5.1 and the pre-1979 streams
#for dir in 1 ; do   # main expid for ERA5
  cd $pdir/$dir
  \rm *.conv.??????._*
  for file in `ls *.conv.??????` ; do
     echo $file

     n=$(($n+1))
     if [ $n -lt 12 ] ; then
        odb split -maxopenfiles 5000 $file "$file.{statid@hdr}" &
     else
        odb split -maxopenfiles 5000 $file "$file.{statid@hdr}" 
        n=0
     fi
      
  done

  nproc=`ps -ef | grep 'odb split' | wc -l`
  while [ $nproc -gt 1 ] ; do
    sleep 10
    nproc=`ps -ef | grep 'odb split' | wc -l`
  done
 
  yy=1900
  while [ $yy -lt 2020 ] ; do
    for file in `ls *.conv.${yy}??.* 2>/dev/null` ; do

      fil=`echo $file | cut -c 1-9`
      suff=`echo $file |  cut -c 18-23`
      echo $yy,$fil $suff, $file
    
      cat $file >> $fil._$suff

      #ls -l  $fil._$suff
    done
    echo $yy
    yy=$(($yy+1))
  done
done
