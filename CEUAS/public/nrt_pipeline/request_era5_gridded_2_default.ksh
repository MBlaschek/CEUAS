#!/bin/ksh
#SBATCH --qos=el
#SBATCH --job-name=monmeanall
#SBATCH --output=monmeanall.$(hostname).%j.out
#SBATCH --error=monmeanall.$(hostname).%j.out
#SBATCH --mail-type=FAIL
#SBATCH --time=7-00:00:00

# Constants
export jobname='jobname'
export host=$(hostname)
export jobid=$SLURM_JOB_ID

YYY=2025
DEC=$((($YYY / 10) * 10))
EXP=5
PARAMS="T/Q/U/V/Z"
LEVELS="10/20/30/50/70/100/150/200/250/300/400/500/700/850/925/1000"
GRID="0.25/0.25"
TIMES="06/18"
STEPS="3/6/9/12"

cd $SCRATCH

MMM=01
DDD=30

cat marsecheader_SLURM > jobpl${EXP}${YYY}${MMM}
cat >> jobpl${EXP}${YYY}${MMM} <<EOF

cd $SCRATCH
cat > marsserapl.${EXP}.${YYY}${MMM} << *eof
retrieve, type=fc, class=ea, expver=${EXP}, stream=oper, param=${PARAMS},
date=${YYY}${MMM}01/to/${YYY}${MMM}${DDD}, time=${TIMES}, step=${STEPS},
decade=${DEC},
level=${LEVELS},
grid=${GRID},
target='era5fc.0.25.${YYY}${MMM}.[param]'
*eof

mars marsserapl.${EXP}.${YYY}${MMM}
(scp era5fc.0.25.${YYY}${MMM}.* leo@aurora.img.univie.ac.at:/mnt/users/scratch/leo/scratch/era5/gridded/ ; rm era5fc.0.25.${YYY}${MMM}.*) &
EOF

ksh jobpl${EXP}${YYY}${MMM}
rm jobpl${EXP}${YYY}${MMM}
done
