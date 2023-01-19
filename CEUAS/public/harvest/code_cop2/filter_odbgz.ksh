#!/bin/ksh

n=0
y=2013
while [ $y -lt 2020 ] ; do

for m in 01 02 03 04 05 06 07 08 09 10 11 12 ; do

echo $y$m

for file in `ls era5.conv.$y$m.?????` ; do

    n=`ps u  | grep 'filter_odbgz.py' | wc -l`
    echo $n
    while [[ $n -gt 30 ]] ; do
      sleep 2
      n=`ps u  | grep 'filter_odbgz.py' | wc -l`
      echo $n
    done

    if [ $y -ge 2013 ] ; then
       (time odb sql -q 'select *' -i ${file} | tr -d " "  > ${file}.txt ; python3 ~/python/filter_odbgz.py ${file} ; pigz -f ${file}.txt ) &
    else
       time odb sql -q 'select *' -i ${file} | tr -d " " | gzip > ${file}.txt.gz &
    fi
       n=$(($n+1))
#    else
#      (time odb sql -q 'select *' -i ${file} | tr -d " "  > ${file}.txt ; python3 ~/python/filter_odbgz.py ${file} ) 
#       n=0
#    fi
  
  
done
done
y=$(($y+1))
done
