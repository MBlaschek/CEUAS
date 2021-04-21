#!/bin/bash
source /etc/profile.d/modules.sh
module load anaconda3/2020.07-gcc-8.3.1
module load netcdf-fortran/4.5.3-intel-20.0.2
module load intel-parallel-studio/composer.2020.2-intel-20.0.2
# unlimit stacksize

export LD_LIBRARY_PATH=/usr/local/lib
cd Temperature_adjustment
python3 import_structure.py
../RISE_FORTRAN/raso_correct_nc radcorpar06