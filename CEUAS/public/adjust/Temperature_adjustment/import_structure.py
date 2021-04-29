import glob,os,sys,shutil

ipath=''
opath='.'

ifiles=['radcorpar06','radcorpar06_24','RS_IRI_CORRTABLES','rtcoef_noaa_14_msu.dat','rtcoef_noaa_16_amsua.dat','mergedstations.t','vapor.instruments.1']

try:
    for f in ifiles:
        shutil.copy(ipath+f,opath)
except:
    pass

try:
    os.mkdir('../common')
except:
    pass
#shutil.copy(ipath+'../common/TSK19582020_s.nc','../common/')
#shutil.copy(ipath+'../common/FAL_1979000000','../common/')
#shutil.copy(ipath+'../common/AAL_1979000000','../common/')
try:
    shutil.copy(ipath+'../common/HadCRUT.4.6.0.0.median.nc','../common/')
    shutil.copy(ipath+'../common/LSM_2000010100_s.nc','../common/')
    shutil.copy(ipath+'../common/ELEV2000010100_s.nc','../common/')
    start='1958'
    end='2020'
    shutil.copy(ipath+'../common/TSK'+start+end+'_s.nc','../common/')
    #shutil.copy(ipath+'../common/LSP_'+start+end+'_s.nc','../common/')
    shutil.copy(ipath+'../common/CI__'+start+end+'_s.nc','../common/')
except:
    pass

# assumption is that script from_cds_to_legacy has been run so that feedbackmerged... files have already been generated..
i = 0
for f in glob.glob(ipath+'*/feedbackmerged??????.nc'):
    print(f,f[-9:-3])
    try:
        os.mkdir(f[-9:-3])
        shutil.copy(f,f[-9:-3])
    except:
        pass
    i += 1
#     if i > 10:
#         break

