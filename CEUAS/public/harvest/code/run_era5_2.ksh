#!/bin/ksh

istart=0
n=0
idir=/raid60/scratch/leo/scratch/era5/odbs/2/
for file in `ls ${idir}/era5.conv._*.gz` ; do

    
    #statnr=`echo "${file#${file%?????}}"` # $file | cut -c 12-17`

    #mfile="${idir}/era5.conv.??????.${statnr}.txt.gz"
    #fil=

    echo $file


    n=`cat /proc/loadavg |cut -f1 -d"."`  #`ps u  | grep 'newfixes' | wc -l`
    m=`grep -i memfree /proc/meminfo | awk '{print $2}'`
    echo 'load and mem' $n $m
    while [[ $n -gt 30 || $m -lt 30178336 ]] ; do
      sleep 5
      n=`cat /proc/loadavg |cut -f1 -d"."`  #`ps u  | grep 'newfixes' | wc -l`
      m=`grep -i memfree /proc/meminfo | awk '{print $2}'`
      echo 'load and mem' $n $m
    done

    time python3 harvest_convert_to_netCDF_newfixes.py -f ${file} -d era5_2 -o ${idir} &

    if [ $istart -lt 5 ] ; then
       sleep 5
    else
       sleep 1
    fi
    istart=$(($istart+1))
    n=$(($n+1))
  
done
