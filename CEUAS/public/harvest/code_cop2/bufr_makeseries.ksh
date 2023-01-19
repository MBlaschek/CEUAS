#!/bin/ksh

n=0
for file in `ls ERA40*.bfr` ; do
  YYMM=`echo $file | cut -c 9-14`
  for file in `ls ${YYMM}.101.*.bfr` ; do
    STATID=`echo $file | awk -F . '{print $3}' -`
    #echo $YYMM$STATID
    #ls ${YYMM}.101.${STATID}.bfr
    if [ $n -lt 10 ]; then 
	cat $YYMM.101.${STATID}.bfr >> era5.${STATID}.bfr&
        n=$(($n+1))
    else
	cat $YYMM.101.${STATID}.bfr >> era5.${STATID}.bfr
        n=0
    fi

  done
done
  
