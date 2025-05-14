#!/bin/bash

module load teleport
source /srvfs/home/uvoggenberger/CEUAS/CEUAS/public/nrt_pipeline/ssh-agent.sh
ssh_agentreconnect
python3 -m teleport.login
ssh -vv ecmwf_user@hpc-login 'mkdir jobs'
scp file_to_modify_1 ecmwf_user@hpc-login:/home/ecmwf_user/jobs/job_1.ksh
scp file_to_modify_2 ecmwf_user@hpc-login:/home/ecmwf_user/jobs/job_2.ksh
scp file_to_modify_3 ecmwf_user@hpc-login:/home/ecmwf_user/jobs/job_3.ksh
ssh -vv ecmwf_user@hpc-login 'chmod +x jobs/job_1.ksh ; jobs/job_1.ksh ; chmod +x jobs/job_2.ksh ; jobs/job_2.ksh ; chmod +x jobs/job_3.ksh ; jobs/job_3.ksh' | tee output.txt
scp ecmwf_user@hpc-login:ecmwf_out_dir* harvest_dir
ssh -vv ecmwf_user@hpc-login 'rm ecmwf_out_dir*'

