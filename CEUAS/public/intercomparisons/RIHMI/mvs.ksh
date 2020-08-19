#!/bin/ksh

for file in `ls` ; do
  echo $file
#  mv $file/$file/* $file/
  rmdir $file/$file
done
