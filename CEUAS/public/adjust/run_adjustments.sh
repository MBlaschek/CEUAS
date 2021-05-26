#!/bin/csh -f
setenv LD_LIBRARY_PATH /usr/local/lib
# export LD_LIBRARAY_PATH="usr/local/lib"
module load openmpi/3.1.6-intel-20.0.2
module load netcdf-fortran/4.5.3-intel-20.0.2
module load intel-parallel-studio/composer.2020.2-intel-20.0.2
limit stacksize unlimited
# ulimit -s unlimited
setenv OMP_NUM_THREADS 1

# download data from cds
python gettemp.py
# compiling fortran code
cd RISE_FORTRAN
\rm *.mod *.o
make raso_correct_nc
cd ..
# creating adjustments
cd Temperature_adjustment/
# copy station data for further use
python import_structure.py
# create RAOBCORE adjusments
../RISE_FORTRAN/raso_correct_nc radcorpar06
ls -l */*corrsave??????.nc | wc -l
# create RISE adjustments
../RISE_FORTRAN/raso_correct_nc radcorpar06_24
ls -l */*corrsave*rio24*.nc | wc -l
cd ..
# write adjustments to outputfiles
cd Converters
python from_cds_to_legacy.py
cd ..
# create RICH and RASE adjustments
cd Converters
python ../Converters/add_solarangle_adjustments.py
ls -l */ERA5bc_*.nc | wc -l
cd ..
