#!/bin/bash
source /etc/profile.d/modules.sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
module load anaconda3/2020.07-gcc-8.3.1
module load netcdf-fortran/4.5.3-intel-20.0.2
module load intel-parallel-studio/composer.2020.2-intel-20.0.2
# unlimit stacksize
ulimit -s unlimited

export OMP_NUM_THREADS 1
cd Temperature_adjustment
python3 import_structure.py
../RISE_FORTRAN/raso_correct_nc radcorpar06
#../RISE_FORTRAN/raso_correct_nc radcorpar06_24