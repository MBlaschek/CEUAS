import glob
import numpy as np
import pandas as pd

if False:
    tlist=glob.glob('/mnt/users/scratch/leo/scratch/era5/odbs/1_mobile/new/*.gz')
    nlist=glob.glob('/mnt/users/scratch/uvoggenberger/CUON_HARVEST/harvest_test_float32_source/era5_1_mobile/*/*.nc')
    olist=glob.glob('/mnt/users/scratch/leo/scratch/UH/CUON_HARVEST/harvest_mobile/era5_1_mobile/*/*.nc')

else:
    nlist=glob.glob('/mnt/users/scratch/uvoggenberger/CUON_HARVEST_2024/harvest_regular/era5_1/*/*.nc')
    nlist=nlist+glob.glob('/mnt/users/scratch/uvoggenberger/CUON_HARVEST_NEW/harvest_orphan/era5_1/*/*.nc')
    nlist=nlist+(glob.glob('/mnt/users/scratch/uvoggenberger/CUON_HARVEST_NEW/additions/era5_1/*/*.nc'))
    # nlist+=glob.glob('/mnt/users/scratch/uvoggenberger/CUON_HARVEST_NEW/harvest_test_float32_source_check/era5_1/*/*.nc')
    olist=glob.glob('/mnt/users/scratch/leo/scratch/UH/CUON_HARVEST/harvest_regular/era5_1/*/*.nc')
    olist=olist+glob.glob('/mnt/users/scratch/leo/scratch/UH/CUON_HARVEST/harvest_orphan/era5_1/*/*.nc')
    tlist=glob.glob('/mnt/users/scratch/leo/scratch/era5/odbs/1/new/*.gz')

#print(nlist[0])
#exit()
nlines=['era5.conv.'+n.split('_era5_1')[0][-4:]+'??.'+n.split('.')[-3]+'.gz.nc' for n in nlist]
# olines=['era5.conv.'+n.split('_era5_1')[0][-4:]+'??.'+n.split('.')[-3]+'.gz.nc' for n in olist]
tlines=[n.split('/')[-1][:14]+'??'+n.split('/')[-1][16:]+'.nc' for n in tlist if '.txt.gz' not in n]
tlines=np.unique(tlines)

#x=['11938', '17196', '17516', '22008', '24122', '24947', '25403', '26075', '27713', '28951', '30557', '36872', '70165', '71126', '72413', '76526', '84384', '96739']
#for xx in x:
#    for n in nlist:
#        if xx in n:
#            print(xx,'/'.join(n.split('/')[-3:]))
#            break
#exit()
#print(tlines)
#print(nlines[0])
#exit()

#with open('new23') as f:
#    nlines=f.read().split('\n')
#with open('old23') as f:
#    olines=f.read().split('\n')

#nlines=[n.split('/')[-1] for n in nlines]
#olines=[n.split('/')[-1] for n in olines]
#print(nlines)
#print(olines)
print('not in tlines:')
osum=0
# for n in tlines:
#     if n not in olines:
# #        print(n)
#         osum+=1

missing_stats = []
for year in range(1880,2024):
    tlines_y = [n for n in tlines if '.'+str(year)+'??' in n]
    nlines_y = [n for n in nlines if '.'+str(year)+'??' in n]

    print('not in nlines:')
    nsum=0
    for n in tlines_y:
        if n not in nlines_y:
            missing_stats.append(n[:-3].split('.')[-2])
            print(n[:-3])
            nsum+=1

print('fehlt in alt:',osum,'fehlt in neu:',nsum)
print(np.unique(missing_stats))

sc = pd.read_csv('/srvfs/home/uvoggenberger/CEUAS/CEUAS/public/harvest/data/station_configurations/era5_1_station_configuration_extended.csv', sep='\t')
ms = np.unique(missing_stats)
pid_list = []
for i in ms:
    pid_list = pid_list + [j for j in sc.primary_id if i in j] + [list(sc[sc.secondary_id == j].primary_id.values)[0] for j in sc.secondary_id if i in j]

print(np.unique(pid_list))
print(*np.unique(pid_list), sep=',')




['0-124-0-73110', '0-20000-0-01002', '0-20000-0-01003', '0-20000-0-01006',
 '0-20000-0-01008', '0-20000-0-01009', '0-20000-0-06242', '0-20000-0-07257',
 '0-20000-0-15615', '0-20000-0-16400', '0-20000-0-16723', '0-20000-0-17292',
 '0-20000-0-17300', '0-20000-0-25745', '0-20000-0-26069', '0-20000-0-28676',
 '0-20000-0-30925', '0-20000-0-32287', '0-20000-0-37958', '0-20000-0-38606',
 '0-20000-0-38705', '0-20000-0-40437', '0-20000-0-41929', '0-20000-0-42260',
 '0-20000-0-47157', '0-20000-0-47404', '0-20000-0-47590', '0-20000-0-47756',
 '0-20000-0-47777', '0-20000-0-51133', '0-20000-0-51288', '0-20000-0-57411',
 '0-20000-0-59046', '0-20000-0-59559', '0-20000-0-60018', '0-20000-0-60402',
 '0-20000-0-62008', '0-20000-0-62176', '0-20000-0-62417', '0-20000-0-64758',
 '0-20000-0-65599', '0-20000-0-67875', '0-20000-0-68816', '0-20000-0-71123',
 '0-20000-0-71742', '0-20000-0-71807', '0-20000-0-72281', '0-20000-0-72364',
 '0-20000-0-72489', '0-20000-0-72776', '0-20000-0-78862', '0-20000-0-83362',
 '0-20000-0-83779', '0-20000-0-91557', '0-20000-0-91929', '0-20000-0-94377',
 '0-20000-0-95721', '0-20000-0-98333', '0-20001-0-57516', '0-20300-0-96773']