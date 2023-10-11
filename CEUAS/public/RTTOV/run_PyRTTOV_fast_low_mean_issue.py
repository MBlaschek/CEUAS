# System information
import numpy as np
import os, sys, glob
import xarray
import xarray as xr
print(sys.executable)
print(sys.version_info)
import pandas as pd
import pandas
pd.set_option('display.max_columns', None)
# pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)
sys.path.append(os.getcwd()+'/../cds-backend/code/')
import cds_eua3 as eua
import pickle
import multiprocessing
from functools import partial

def calc_station(statid, chum, odir, adj = None):
    statlist = statid
    statid = statlist.split('.nc')[0][-5:]
    print(statid)
    
    fbmean = glob.glob('/mnt/ssdraid/scratch/leo/rise/1.0/exp02/*'+str(statid)+'*/feedbackglobbincorrmon0*.nc')[0]
    dfm = xr.open_dataset(fbmean).to_dataframe()
    dfm.press = dfm.press * 100.
    dfm = dfm[dfm.press.isin([3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000])]
    dfm = dfm.loc[dfm.index.get_level_values('numdat') == 0]
    dfm = dfm.reset_index()    
    dfm = dfm.rename({'time':'time_idx'}, axis='columns')
    dfm = dfm.rename({'datum':'time', 'press':'plev'}, axis='columns')

    df = xr.open_dataset(statlist).to_dataframe()
    df.press = df.press * 100.
    ###
    df = df[df.press.isin([3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000])]
    df = df[df.datum > '1950']
    df = df.reset_index()
    df = df.rename({'time':'time_idx'}, axis='columns')
    df = df.rename({'datum':'time', 'press':'plev'}, axis='columns')

    all_dfta = df
    all_dfta = all_dfta.rename({'temperatures':'ta'}, axis='columns')

    for day in [True,False]:
        if day:
            dn_dfta = all_dfta.loc[all_dfta.hour == 12]
            m_dt = 1
        else:
            dn_dfta = all_dfta.loc[all_dfta.hour == 0]
            m_dt = 0
        goodmon = []
        for yr in range(1995,1996,1):
            for mon in range(int(str(yr)+'01'), int(str(yr)+'7'), 1):
                
                dfta = dn_dfta.loc[(dn_dfta['time'].dt.year==int(str(mon)[:4])) 
                                   & (dn_dfta['time'].dt.month==int(str(mon)[-2:]))
                                  ]
                idfm =  dfm.loc[(dfm['time'].dt.year==int(str(mon)[:4])) 
                                & (dfm['time'].dt.month==int(str(mon)[-2:])) 
                                & (dfm['hour']==m_dt) 
                               ]
                df = dfta
                mon_mean = df.groupby(['plev']).aggregate({"ta":np.mean})
                print('calc_mean: ',mon_mean)
                print('file_mean: ', idfm)

    
    return



if __name__ == '__main__': 
    odir = "rttov_out_unadj_testing"
    try:
        os.makedirs("./"+odir+"/")
    except:
        pass
#     [3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000]
    consthum = np.array([3.90, 5.90, 9.17, 20.30,
                         285.00, 1464.00 , 2475.60, 6631.20,
                         15468.00, 21684.00, 35328.00 , 44220.00]
                       )/2.
    statlist = []
    statlist = glob.glob('/mnt/ssdraid/scratch/leo/rise/1.0/exp02/*/feedbackmerged*.nc')
    statlist = glob.glob('/mnt/ssdraid/scratch/leo/rise/1.0/exp02/*11035*/feedbackmerged*.nc')

#     statlist = glob.glob('/mnt/ssdraid/scratch/leo/rise/1.0/exp02/*11035*/feedbackmerged*.nc')
#     stats = glob.glob('/mnt/users/scratch/leo/scratch/converted_v7/*68842*_CEUAS_merged_v1.nc')
#     for i in stats:
#         statlist.append(i.split('-')[-1][:5])
#     calc_station(statlist[0])
#     for i in ["RISE", "RASE"]: #[None, "RAOBCORE", "RICH", "RISE", "RASE"]:
    i = None
    
#     for j in statlist:
#         print(j)
# #         try:
#         calc_station(j, chum = consthum, adj = i, odir = odir)
# #             print('done')
# #         except:
# #             print('failed')

    pool = multiprocessing.Pool(processes=1)
    func=partial(calc_station, chum = consthum, adj = i, odir = odir)
    result_list = list(pool.map(func, statlist[:]))
    print(result_list)