#!/bin/ksh
#
# script for converting odb station files into gzipped ASCII format
#
# Limits for n should be adjusted depending on the server capacity
#
#
# Leopold Haimberger, 23 January 2020
#

for file in $(ls era5.*conv.* | grep -v '.nc'); do
  echo $file
  fil=$(echo $file | cut -c 1-14)
  suff=$(echo $file | cut -c 16-23)
  echo $fil,$suff

  n=$(ps u | grep 'odb sql' | wc -l)
  echo $n
  rm ${fil}._$suff.gz
  if [ ! -f ${fil}._$suff.gz ]; then
    while [[ $n -gt 20 ]]; do
      sleep 5
      n=$(ps u | grep 'odb sql' | wc -l)
      echo $n
    done

    time odb sql -q 'select *' -i ${file} | tr -d " " | gzip >${fil}._${suff}.gz &
  fi

done
