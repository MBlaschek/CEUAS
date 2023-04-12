   637	15:02	grep 63705 *.csv | grep 50000.0
   638	15:03	grep 47744 *.csv | grep 50000.0
   639	15:03	grep 74574 *.csv | grep 50000.0
   640	15:03	grep 52652 *.csv | grep 50000.0
   641	15:03	grep 52652 *.csv | grep 30000.0
   642	15:05	ls /mnt/scratch/scratch/federico/COP2_HARVEST_APRIL2022/era5_1/*63705*
   643	15:05	grep 11035 *.csv | grep 50000.0
   645	15:06	ls /mnt/scratch/scratch/federico/COP2_HARVEST_APRIL2022/era5_1/*47744*
   646	15:06	ls /mnt/scratch/scratch/federico/COP2_HARVEST_APRIL2022/era5_1/*74547*
   647	15:07	ls /mnt/scratch/scratch/federico/COP2_HARVEST_APRIL2022/era5_1/*74646*
   648	15:07	ls /mnt/scratch/scratch/federico/COP2_HARVEST_APRIL2022/era5_1/*74650*
   649	15:07	ls /mnt/scratch/scratch/federico/COP2_HARVEST_APRIL2022/era5_1/*52652*
   650	15:08	ls /mnt/scratch/scratch/federico/COP2_HARVEST_APRIL2022/era5_1/*55664*
   651	15:08	ls /mnt/scratch/scratch/federico/COP2_HARVEST_APRIL2022/era5_1/*56173*
   652	23:11	ls -lrt $SCRATCH/ectrans
   655	23:14	ls -lrt $RSCRATCH/ectrans/era5/odbs
   658	23:23	ls -lrt $RSCRATCH/ectrans
   659	23:23	ls -lrt $RSCRATCH/ectrans/*1979*
   660	23:25	cd $RSCRATCH/ectrans/
   661	23:25	h | grep odb
   665	9:23	ls -lrt *.conv.*12
   666	9:23	ls -lrt *.conv.198*
   669	9:43	ls -lrt *.conv.1980*
   670	9:43	ls -lrt *.conv.1979*
   673	9:48	ls -lrt *.conv.1959*
   674	9:50	ls -lrt *.conv.1958*
   676	10:49	ls -lrt *.conv.2010*
   680	8:34	du
   681	8:34	du -h
   683	9:32	df -h .
   689	15:43	history
   691	15:45	ls -lrt *.conv.*
   692	15:49	odb sql -q 'select statid,obstype,codetype,bufrtype,subtype where vertco_reference_1=50000 and varno=2 and date=19970701' -i era5.conv.199707
   694	15:55	ls -lrt
   695	15:55	python checkschroed.py
   699	16:03	cd
   700	16:03	python checkschroed.py | less
   701	16:18	back
   702	16:37	cd ~/tables
   703	16:38	mkdir world-weather-records
   704	16:38	cd world-weather-records/
   705	16:38	wget -1_era5_2_harvested_era5.conv._6:11012.gz.nc
   706	16:39	wget https://www.ncei.noaa.gov/data/world-weather-records/series-6/access/data/WWR_v01_global_array_1961-1970.txt
   707	16:39	wget https://www.ncei.noaa.gov/data/world-weather-records/series-7/access/data/WWR_v01_global_array_1911-1980.txt
   708	16:40	wget https://www.ncei.noaa.gov/data/world-weather-records/series-7/access/data/WWR_v01_global_array_1971-1980.txt
   709	16:41	wget https://www.ncei.noaa.gov/data/world-weather-records/series-8/access/data/WWR_v01_global_array_1981-1990.txt
   710	16:42	wget https://www.ncei.noaa.gov/data/world-weather-records/series-9/access/data/WWR_v01_global_array_1991-2000.txt
   711	16:43	wget https://www.ncei.noaa.gov/data/world-weather-records/series-9/doc/WWR_v02_global_StnLst.txt
   712	16:44	wget https://www.ncei.noaa.gov/data/world-weather-records/series-10/doc/StnLstWWR_Region00_2001-2010.txt
   713	16:44	wget https://www.ncei.noaa.gov/data/world-weather-records/series-10/doc/StnLstWWR_Region01_2001-2010.txt
   714	16:44	wget https://www.ncei.noaa.gov/data/world-weather-records/series-10/doc/StnLstWWR_Region02_2001-2010.txt
   715	16:44	wget https://www.ncei.noaa.gov/data/world-weather-records/series-10/doc/StnLstWWR_Region03_2001-2010.txt
   716	16:44	wget https://www.ncei.noaa.gov/data/world-weather-records/series-10/doc/StnLstWWR_Region04_2001-2010.txt
   717	16:44	wget https://www.ncei.noaa.gov/data/world-weather-records/series-10/doc/StnLstWWR_Region05_2001-2010.txt
   718	16:44	wget https://www.ncei.noaa.gov/data/world-weather-records/series-10/doc/StnLstWWR_Region06_2001-2010.txt
   719	16:45	wget https://www.ncei.noaa.gov/data/world-weather-records/series-11/doc/StnLstWWR_Region06_2011-2016.txt
   720	16:45	wget https://www.ncei.noaa.gov/data/world-weather-records/series-11/doc/StnLstWWR_Region05_2011-2016.txt
   721	16:45	wget https://www.ncei.noaa.gov/data/world-weather-records/series-11/doc/StnLstWWR_Region04_2011-2016.txt
   722	16:45	wget https://www.ncei.noaa.gov/data/world-weather-records/series-11/doc/StnLstWWR_Region03_2011-2016.txt
   723	16:45	wget https://www.ncei.noaa.gov/data/world-weather-records/series-11/doc/StnLstWWR_Region02_2011-2016.txt
   724	16:45	wget https://www.ncei.noaa.gov/data/world-weather-records/series-11/doc/StnLstWWR_Region01_2011-2016.txt
   725	16:45	wget https://www.ncei.noaa.gov/data/world-weather-records/series-11/doc/StnLstWWR_Region00_2011-2016.txt
   726	16:46	ls
   727	16:47	cat *.txt | cut -c 1-5 | grep 11012
   728	16:47	cat *.txt | cut -c 1-5
   729	16:47	cat *.txt | cut -c 0-5
   730	16:47	cat *.txt | cut -c 1-7
   731	16:48	cat *.txt | cut -c 1-7 | grep 11012
   732	18:39	h
   733	18:39	h > readme.txt
