import glob
import numpy as np

if False:
    tlist=glob.glob('/mnt/users/scratch/leo/scratch/era5/odbs/1_mobile/new/*.gz')
    nlist=glob.glob('/mnt/users/scratch/uvoggenberger/CUON_HARVEST/harvest_test_float32_source/era5_1_mobile/*/*.nc')
    olist=glob.glob('/mnt/users/scratch/leo/scratch/UH/CUON_HARVEST/harvest_mobile/era5_1_mobile/*/*.nc')

else:
    nlist=glob.glob('/mnt/users/scratch/uvoggenberger/CUON_HARVEST/harvest_test_float32_source/era5_1/*/*.nc')
    nlist=nlist+glob.glob('/mnt/users/scratch/uvoggenberger/CUON_HARVEST/harvest_test_float32_source_orphan/era5_1/*/*.nc')
    nlist+=glob.glob('/mnt/users/scratch/uvoggenberger/CUON_HARVEST/harvest_test_float32_source_check/era5_1/*/*.nc')
    olist=glob.glob('/mnt/users/scratch/leo/scratch/UH/CUON_HARVEST/harvest_regular/era5_1/*/*.nc')
    olist=olist+glob.glob('/mnt/users/scratch/leo/scratch/UH/CUON_HARVEST/harvest_orphan/era5_1/*/*.nc')
    tlist=glob.glob('/mnt/users/scratch/leo/scratch/era5/odbs/1/new/*.gz')

#print(nlist[0])
#exit()
nlines=['era5.conv.'+n.split('_era5_1')[0][-4:]+'??.'+n.split('.')[-3]+'.gz.nc' for n in nlist]
olines=['era5.conv.'+n.split('_era5_1')[0][-4:]+'??.'+n.split('.')[-3]+'.gz.nc' for n in olist]
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
for n in tlines:
    if n not in olines:
#        print(n)
        osum+=1

print('not in nlines:')
nsum=0
for n in tlines:
    if n not in nlines:
        print(n[:-3])
        nsum+=1

print('fehlt in alt:',osum,'fehlt in neu:',nsum)
print(len(tlines))