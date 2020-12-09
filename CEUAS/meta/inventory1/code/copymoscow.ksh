#!/bin/ksh


cp $RSCRATCH/era5/odbs/3188/era5.3188.conv.1:27612 ../data/3188/
cp $RSCRATCH/era5/odbs/1759/era5.1759.conv.1:27612 ../data/1759/
cp $RSCRATCH/era5/odbs/1761/era5.1761.conv.1:27612 ../data/1761/
cp $RSCRATCH/era5/odbs/ai_bfr/era5.27612.bfr ../data/ai_bfr/
cp $RSCRATCH/era5/odbs/rda/*27612*.nc ../data/rda/

cp $RSCRATCH/era5/odbs/3188/meta.csv ../data/3188/meta_full.csv
cp $RSCRATCH/era5/odbs/1759/meta.csv ../data/1759/meta_full.csv
cp $RSCRATCH/era5/odbs/1761/meta.csv ../data/1761/meta_full.csv
cp $RSCRATCH/era5/odbs/ai_bfr/meta.csv ../data/ai_bfr/meta_full.csv
cp $RSCRATCH/era5/odbs/rda/meta.csv ../data/rda/meta_full.csv

cp $RSCRATCH/era5/odbs/all/*.png ../data/all/full/
