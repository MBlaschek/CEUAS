import numpy
import numpy as np
import pandas as pd
import sys, glob
import urllib3
import h5py
import cdsapi, zipfile, os, time
sys.path.append(os.getcwd()+'/../cds-backend/code/')
import cds_eua3 as eua
import warnings
import shutil
import pickle
warnings.filterwarnings('ignore')

ihfiles = glob.glob('/raid60/scratch/uli/IGRA_HUM/IGRA_H*/*.nc')
ifiles = glob.glob('/raid60/scratch/uli/IGRA_HUM/IGRA_*/*.nc')
for i in ihfiles:
    ifiles.remove(i)
    
for i in ifiles:
    if i == ifiles[0]:
        ds = xarray.open_dataset(i)            
        data = ds.to_dataframe()
    else:
        ds = xarray.open_dataset(i)            
        adddata = ds.to_dataframe()
        data = data.append(adddata, ignore_index=True)
display(data)
idata = data
    
# for i in ihfiles:
#     if i == ihfiles[0]:
#         ds = xarray.open_dataset(i)            
#         data = ds.to_dataframe()
#     else:
#         ds = xarray.open_dataset(i)            
#         adddata = ds.to_dataframe()
#         data = data.append(adddata, ignore_index=True)
# display(data)
# ihdata = data


stations = idata pickle.load( open( "stations.p", "rb" ))
sout = []
start = 1979
end = 2019
intervall = end - start
for i in range(len(stations.station_name)):
    files = glob.glob('/raid60/scratch/uli/IGRA_Data/COMP/COMP_' + str(stations.station_name.iloc[i]) + '*.p')
#     print(files)
    if len(files) > 0:
        temp = idata[idata.station_name == str(stations.station_name.iloc[i])]
        temp = temp[temp.plev == 30000]
#         temp = temp[temp.time > '1979-01-01 00:00:00'][temp.time < '2020-01-01 00:00:00']
        if len(temp) > 0:
            temp.sort_values('time')
            temp.time = pd.to_datetime(temp['time'])
            temp.lat = numpy.array([temp.lat.iloc[0]]*len(temp))
            temp.lon = numpy.array([temp.lon.iloc[0]]*len(temp))
            temptime = temp.time
            if len(temp) >= 17*365 and len(numpy.unique(temptime.dt.year)) > 17:
                temp.ta = temp.ta - numpy.nan_to_num(temp.bias_estimate)
                xa = temp.set_index(['lat', 'lon', 'time']).to_xarray()
                out = rasotools.met.time.trend(xa.ta,only_slopes=True).to_dataframe(name='out')
                sout.append(float(out.iloc[-1] *3650))
            else:
                sout.append(np.nan)
        else: sout.append(np.nan)
    else: sout.append(np.nan)
stations.station_name = sout
pickle.dump( stations, open( "new_CUON_BC_100hPa_2000_2019_Trend.p", "wb" ))



stations = pickle.load( open( "stations.p", "rb" ))
for i in stations.station_name:
    if not '/raid60/scratch/uli/IGRA_Data/COMP/COMP_'+str(i)+'_100.p' in glob.glob('/raid60/scratch/uli/IGRA_Data/COMP/*'):
        try:
            data = request({'variable':['temperature'],
                            'statid': i[-5:],
                            'period': '1979-01-01/2019-12-31',
                            'optional':'bias_estimate',
                            'pressure_level': 100,
                           }, 'insitu-comprehensive-upper-air-observation-network',
            )
        except:
            pass
        pickle.dump( data, open( '/raid60/scratch/uli/IGRA_Data/COMP/COMP_' + str(i) + '_100.p', "wb" ))
    else:
        print('-')