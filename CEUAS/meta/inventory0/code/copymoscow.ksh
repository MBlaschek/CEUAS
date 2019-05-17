#!/bin/ksh


for f in '4525' '4526' '4527' '4528' '4529' '4530' '4531' '4532' '4533' '4534' '4567' '4629' '5195A' '5246A' '5195B' '5246B' 
do

  cp $RSCRATCH/era5/odbs/3188/era5.3188.conv.C:${f} ../data/3188/
done

cp $RSCRATCH/era5/odbs/1759/era5.1759.conv.1:27612 ../data/1759/
cp $RSCRATCH/era5/odbs/1761/era5.1761.conv.1:27612 ../data/1761/
cp $RSCRATCH/era5/odbs/ai_bfr/era5.27612.bfr ../data/ai_bfr/
cp $RSCRATCH/era5/odbs/rda/*27612*.nc ../data/rda/

cp $RSCRATCH/era5/odbs/3188/meta.csv ../data/3188/meta_full.csv
cp $RSCRATCH/era5/odbs/1759/meta.csv ../data/1759/meta_full.csv
cp $RSCRATCH/era5/odbs/1761/meta.csv ../data/1761/meta_full.csv
cp $RSCRATCH/era5/odbs/ai_bfr/meta.csv ../data/ai_bfr/meta_full.csv
cp $RSCRATCH/era5/odbs/rda/meta.csv ../data/rda/meta_full.csv
