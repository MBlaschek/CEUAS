#!/bin/ksh

ulimit -n 4096

#for file in `ls *.conv.??????` ; do

#    odb split -maxopenfiles 4000 $file "$file.{statid@hdr}"
#done

yy=1900
while [ $yy -lt 1980 ] ; do
mm=1
while [ $mm -lt 13 ] ; do
for file in `ls *.conv.${yy}${mm}.* 2>/dev/null` ; do

    fil=`echo $file | cut -c 1-15`
    suff=`echo $file |  cut -c 21-30`
    echo $yy,$fil$suff,$file
    cat $file >> $fil$suff
    exit
done
mm=$(($mm+1))
done
yy=$(($yy+1))
done
