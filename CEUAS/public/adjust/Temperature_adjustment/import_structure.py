import glob,os,sys,shutil

ipath='/ssdraid/scratch/leo/rise/1.0/exp06/'
opath='.'

ifiles=['radcorpar06','radcorpar06_24','RS_IRI_CORRTABLES','rtcoef_noaa_14_msu.dat','rtcoef_noaa_16_amsua.dat','mergedstations.t','vapor.instruments.1']

for f in ifiles:
    shutil.copy(ipath+f,opath)

try:
    os.mkdir('../common')
except:
    pass
shutil.copy(ipath+'../common/TSK19582020_s.nc','../common/')
shutil.copy(ipath+'../common/FAL_1979000000','../common/')
shutil.copy(ipath+'../common/AAL_1979000000','../common/')

# assumption is that script from_cds_to_legacy has been run so that feedbackmerged... files have already been generated..
for f in glob.glob(ipath+'*/feedbackmerged??????.nc'):
    print(f,f[-9:-3])
    try:
        os.mkdir(f[-9:-3])
        shutil.copy(f,f[-9:-3])
    except:
        pass
    break

