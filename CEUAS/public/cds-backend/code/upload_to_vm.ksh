#!/bin/ksh

cd ~/python/web2py/applications
tar -cvf eua.tar eua
scp eua.tar  sis@136.156.133.52:/data/private/soft/web2py/applications/
ssh sis@136.156.133.52 "cd /data/private/soft/web2py/applications && tar -xvf eua.tar"
ssh sis@136.156.133.52 "cd /data/private/soft/web2py/applications/eua/static&& rm monthly&& rm subdaily && ln -s /data/public/eua/subdaily . && ln -s /data/public/eua/monthly ."
scp ~/python/cds_eua.py sis@136.156.133.52:python/
scp ~/python/web2py/routes.py sis@136.156.133.52:/data/private/soft/web2py/

exit
cd $RSCRATCH/eua
for file in 01001 01009 ; do 
  cp /fio/srvx7/leo/python/CEUAS/CEUAS/public/harvest/code/../data/tables/chera5.conv._${file}.nc /raid60/scratch/leo/scratch/eua/subdaily/v0.1/source/ERA5_1/obs/0-20000-0-${file}/eua_subdaily_v0.1_source_ERA5_1_obs_0-20000-0-${file}_t.nc 
done

for d in subdaily monthly ; do
  tar -cvf ${d}.tar ${d}
  ssh sis@136.156.133.52 "cd /data/public/eua && rm -r ${d}/*"
  scp ${d}.tar sis@136.156.133.52:/data/public/eua/
  ssh sis@136.156.133.52 "cd /data/public/eua && tar -xvf ${d}.tar"
done


