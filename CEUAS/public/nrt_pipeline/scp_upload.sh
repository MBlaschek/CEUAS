#!/bin/bash
ssh -J proxy@136.156.142.19  obs@dbdataset-cci2-0000 'mkdir /mnt/public/202502/'
scp /mnt/users/scratch/uvoggenberger/CUON_HARVEST_202502/resort/2025/*.nc -J proxy@136.156.142.19  obs@dbdataset-cci2-0000:/mnt/public/202502/
scp /srvfs/home/uvoggenberger/CEUAS/CEUAS/public/nrt_pipeline/ingestion_script.sh -J proxy@136.156.142.19  obs@dbdataset-cci2-0000:/mnt/public/202502/
ssh -J proxy@136.156.142.19  obs@dbdataset-cci2-0000 'chmod +x /mnt/public/202502/remote_ingestion_script.sh; /mnt/public/202502/remote_ingestion_script.sh'
