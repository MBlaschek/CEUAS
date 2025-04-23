#!/bin/bash

module load teleport
source /srvfs/home/uvoggenberger/CEUAS/CEUAS/public/nrt_pipeline/ssh-agent.sh
ssh_agentreconnect
python3 -m teleport.login
scp file_to_modify_1 lh4@hpc-login:/home/lh4/jobs/job_1.ksh
scp file_to_modify_2 lh4@hpc-login:/home/lh4/jobs/job_2.ksh
scp file_to_modify_3 lh4@hpc-login:/home/lh4/jobs/job_3.ksh
ssh -vv lh4@hpc-login 'chmod +x jobs/job_1.ksh ; jobs/job_1.ksh ; chmod +x jobs/job_2.ksh ; jobs/job_2.ksh ; chmod +x jobs/job_3.ksh ; jobs/job_3.ksh'
# scp era5.conv.202503.* uvoggenberger@aurora.img.univie.ac.at:/mnt/users/scratch/uvoggenberger/CUON_HARVEST_202504/data/era5_1_data/
