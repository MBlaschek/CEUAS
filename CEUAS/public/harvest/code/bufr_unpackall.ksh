#!/bin/ksh

n=0
for file in `ls ERA40*.bfr` ; do
  YYMM=`echo $file | cut -c 9-14`
  echo $YYMM

cat >rules.$YYMM<<EOF
if (dataSubCategory ==101 || dataSubCategory == 91) {
write "${YYMM}.[dataSubCategory].[blockNumber][stationNumber].bfr";
}
EOF
if [ $n -lt 30 ] ; then 
   bufr_filter rules.$YYMM ERA40_ai${YYMM}.bfr &
   n=$(($n+1))
else
   bufr_filter rules.$YYMM ERA40_ai${YYMM}.bfr
   n=0
fi

done
