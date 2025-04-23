#!/bin/ksh
#SBATCH --qos=el
#SBATCH --job-name=monmeanall
#SBATCH --output=monmeanall.%j.out
#SBATCH --error=monmeanall.%j.out
#SBATCH --mail-type=FAIL
#SBATCH --time=7-00:00:00

# Set environment variables
export jobname='jobname'
export host=host
export jobid=jobid

# Define year, month, and experiment
YYY=2025
MMM=03
EXP=5

# Number of days in March
DDD=31

# Calculate decade
DEC=$(( ($YYY / 10) * 10 ))

# Create mars request script
cd $SCRATCH

cat marsecheader_SLURM > jobconv${EXP}${YYY}${MMM}
cat >> jobconv${EXP}${YYY}${MMM} <<EOF

cd $SCRATCH
cat > marsseraconv.${EXP}.${YYY}${MMM} << *eof
retrieve, type=an, class=ea, expver=${EXP}, stream=oper,
date=${YYY}${MMM}01/to/${YYY}${MMM}${DDD}, time=00/12,
decade=${DEC},
reportype=16013/16022/16045/16068/16014/16019/16046/16069/16021/16075,
type=ofb,
target="${SCRATCH}/era5.conv.${YYY}${MMM}".
*eof

mars marsseraconv.${EXP}.${YYY}${MMM}
EOF

# Run the job script
ksh jobconv${EXP}${YYY}${MMM}

