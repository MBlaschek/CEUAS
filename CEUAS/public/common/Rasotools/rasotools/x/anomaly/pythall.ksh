#!/bin/ksh

for file in `ls *.py` ; do

   pythran $file 2>${file}.log &

done
