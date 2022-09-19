import glob, os, sys
import xarray as xr
import numpy as np
import pandas as pd
import multiprocessing
from functools import partial
import pickle

sys.path.append(os.getcwd()+'/../cds-backend/code/')
import cds_eua4 as eua

import netCDF4
#First import the netcdf4 library
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/

from datetime import datetime, timedelta
def to_days_since(time_elapsed: int, reference=None):
    if reference is None:
        reference = datetime(1900, 1, 1)
    time_elapsed = time_elapsed*100 + 2
    fstart = [datetime.strptime(str(date), '%Y%m%d') for date in time_elapsed]
    out = [int((dt - reference).days) for dt in fstart]
    return out

def write_bt_file(stat, adj):
    print(stat, adj)
#     folder_date = '20220906'
    if adj == 'unadj':
        adjid = '_unadj_testing' # '_unadj_20220906'
        outname = 'feedbackglobbinmon_bt2_0' #  unadj'
    elif adj == 'raobcore':
        adjid = '_raobcore_testing' # '_raobcore_20220906'
        outname = 'feedbackglobbincorrmon_bt2_0' # raobcore
    elif adj == 'rio':
        adjid = '_rio_testing' # '_raobcore_20220906'
        outname = 'feedbackglobbincorrmon_rio24_bt2_0' # rio
    elif adj == 'rit':
        adjid = '_rit_testing' # '_raobcore_20220906'
        outname = 'feedbackglobbincorrmon_rit24_bt2_0' # rit
        
    dirid = './rttov_out'+adjid+'/'

    
#     try:
    outdir = './bt_nc_testing/0'+str(stat)+'/'
    with eua.CDMDataset(glob.glob('/users/staff/a1400070/scratch/converted_v8/*'+str(stat)+'*_CEUAS_merged_v1.nc')[0]) as conv:
        lat = conv.observations_table.latitude[-1]
        lon = conv.observations_table.longitude[-1]
        try:
            station_name = conv.station_configuration.station_name[:]
            station_name = ''.join([i.decode() for i in station_name[0]])
        except:
            station_name = str(stat)
    print(station_name)

    dr = glob.glob(dirid + str(stat) + '/*day_refl.p')[0]
    nr = glob.glob(dirid + str(stat) + '/*night_refl.p')[0]
    dd = glob.glob(dirid + str(stat) + '/*day_dates.p')[0]
    nd = glob.glob(dirid + str(stat) + '/*night_dates.p')[0]
    dg = glob.glob(dirid + str(stat) + '/*day_goodmon.p')[0]
    ng = glob.glob(dirid + str(stat) + '/*night_goodmon.p')[0]

    out = {}
    for i in range(0,4):
        out['nrl_'+str(i+1)] = []
        out['drl_'+str(i+1)] = []
    with open(dr, "rb") as input_file:
        drl = pickle.load(input_file)
        for i in range(len(drl)):
            for j in range(0,4):
                out['drl_'+str(j+1)].append(drl[i][0][j])
    with open(nr, "rb") as input_file:
        nrl = pickle.load(input_file)
        for i in range(len(nrl)):
            for j in range(0,4):
                out['nrl_'+str(j+1)].append(nrl[i][0][j])
    with open(dd, "rb") as input_file:
        out['ddl'] = pickle.load(input_file)
    with open(nd, "rb") as input_file:
        out['ndl'] = pickle.load(input_file)
    with open(dg, "rb") as input_file:
        out['dgm'] = pickle.load(input_file)
    with open(ng, "rb") as input_file:
        out['ngm'] = pickle.load(input_file)
    odf = pd.DataFrame.from_dict(out)
#     print(odf)
#     odf = odf.dropna()
    for col in odf:
        if odf[col].isnull().all():
            return 1
    if len(odf) < 1:
        print('all nan input')
        return 1
    odf = odf.reset_index()
#     print(odf)

#         adj_check = False
#         if (len(glob.glob('./rttov_adj_out/' + str(stat) + '/*day_refl.p')) > 0):
#             adj_check = True
#         if adj_check:
#             a_dr = glob.glob('./rttov_adj_out/' + str(stat) + '/*day_refl.p')[0]
#             a_nr = glob.glob('./rttov_adj_out/' + str(stat) + '/*night_refl.p')[0]
#             a_dd = glob.glob('./rttov_adj_out/' + str(stat) + '/*day_dates.p')[0]
#             a_nd = glob.glob('./rttov_adj_out/' + str(stat) + '/*night_dates.p')[0]
#             a_dg = glob.glob('./rttov_unadj_out/' + str(stat) + '/*day_goodmon.p')[0]
#             a_ng = glob.glob('./rttov_unadj_out/' + str(stat) + '/*night_goodmon.p')[0]
#             a_out = {}
#             for i in range(0,4):
#                 a_out['nrl_'+str(i+1)] = []
#                 a_out['drl_'+str(i+1)] = []
#             with open(a_dr, "rb") as input_file:
#                 drl = pickle.load(input_file)
#                 for i in range(len(drl)):
#                     for j in range(0,4):
#                         a_out['drl_'+str(j+1)].append(drl[i][0][j])
#             with open(a_nr, "rb") as input_file:
#                 nrl = pickle.load(input_file)
#                 for i in range(len(nrl)):
#                     for j in range(0,4):
#                         a_out['nrl_'+str(j+1)].append(nrl[i][0][j])
#             with open(dd, "rb") as input_file:
#                 a_out['ddl'] = pickle.load(input_file)
#             with open(nd, "rb") as input_file:
#                 a_out['ndl'] = pickle.load(input_file)
#             with open(dg, "rb") as input_file:
#                 a_out['dgm'] = pickle.load(input_file)
#             with open(ng, "rb") as input_file:
#                 a_out['ngm'] = pickle.load(input_file)
#             a_odf = pd.DataFrame.from_dict(a_out)
#             a_odf = a_odf.dropna()
#             if len(a_odf) < 1:
#                 adj_check = False

    out = {}
    # only using channel 2,3,4:
    for i in range(2,5):
        if i == 2:
            out['montemp'] = []
            out['goodmon'] = []
#                 out['rasocorrmon'] = []
#                 out['eracorrmon'] = []
#                 out['ancorrmon'] = []
            out['press'] = []
            out['datum'] = []
            out['hour'] = []
            out['lat'] = []
            out['lon'] = []
        for j in range(len(odf.ddl)):
            for k in ['n','d']:
                out['montemp'].append(odf[odf.ddl == odf.ddl[j]][k+'rl_'+str(i)].iloc[0])
#                 print(odf[k+'gm'])
#                 print(odf[k+'gm'][j])
                out['goodmon'].append(odf[k+'gm'][j])
#                     if adj_check:
#                         if  (len(a_odf[a_odf.ddl == j])) < 1:
#                             out['rasocorrmon'].append(np.nan)
#                         else:
#                             out['rasocorrmon'].append(out['montemp'][-1] - (a_odf[a_odf.ddl == j][k+'rl_'+str(i)].iloc[0]))
#                     else:
#                         out['rasocorrmon'].append(np.nan)
#                     out['eracorrmon'].append(np.nan)
#                     out['ancorrmon'].append(np.nan)
                out['press'].append(i)
                out['datum'].append(odf.ddl[j])
                out['lat'].append(lat)
                out['lon'].append(lon)      
                if k == 'n':
                    out['hour'].append(0)
                else:
                    out['hour'].append(1)

    resorted_df = pd.DataFrame.from_dict(out)

    # Read en existing NetCDF file and create a new one
    # f is going to be the existing NetCDF file from where we want to import data
    # and g is going to be the new file.

    try:
        os.makedirs(outdir)
    except:
        pass

    f=Dataset('bt_template.nc','r') # r is for read only
    g=Dataset(outdir+outname+str(stat)+'.nc','w') # w if for creating a file
                                          # if the file already exists it  
                                          # file will be deleted 


    # To copy the global attributes of the netCDF file  

    for attname in f.ncattrs():
        if attname == 'Stationname':
            setattr(g,attname,station_name)
        else:
            setattr(g,attname,getattr(f,attname))

    # To copy the dimension of the netCDF file

    for dimname,dim in f.dimensions.items():
        if dimname == 'time':
            g.createDimension(dimname,len(resorted_df.datum.unique()))
        else:
            g.createDimension(dimname,len(dim))


    # To copy the variables of the netCDF file
    for varname,ncvar in f.variables.items():
        
        if varname in ['montemp', 'lat', 'lon', 'press', 'datum']:
            var = g.createVariable(varname,ncvar.dtype,ncvar.dimensions)
            #Proceed to copy the variable attributes
            for attname in ncvar.ncattrs():  
                setattr(var,attname,getattr(ncvar,attname))
                
    var = g.createVariable('goodmon',int,ncvar.dimensions)

    #Finally copy the variable data to the new created variable
    g['lat'][:] = resorted_df.lat.iloc[-1]
    g['lon'][:] = resorted_df.lon.iloc[-1]
    g['press'][:] = [2,3,4]
    g['datum'][0,:] = to_days_since(np.array(resorted_df.datum.unique()))

    fillvars = {}
    vars_to_write = ['montemp', 'goodmon']#, 'rasocorrmon', 'eracorrmon', 'ancorrmon']
    for i in vars_to_write:
        fillvars[i] = np.empty([len(resorted_df.hour.unique()), len(resorted_df.press.unique()), len(resorted_df.datum.unique())])

    u_hours = resorted_df.hour.unique()
    u_press = resorted_df.press.unique()
    u_datum = resorted_df.datum.unique()
    for i in range(len(u_hours)):
        for j in range(len(u_press)):
            for k in range(len(u_datum)):
                data = resorted_df[(resorted_df.hour == u_hours[i]) & (resorted_df.press == u_press[j]) & (resorted_df.datum == u_datum[k])]
                for l in vars_to_write:
                    fillvars[l][i,j,k] = data[l]

    for i in vars_to_write:
        g[i][:] = fillvars[i]


    f.close()
    g.close()
    return 0
    
#     except Exception as err:
#         print(stat, err)
#         return [stat, err]


if __name__ == '__main__': 
    for o in [ 'rio', 'rit', 'unadj', 'raobcore']: # 'rio', 'rit', 'unadj', 'raobcore',
        fs = glob.glob('./rttov_out_'+o+'_testing/*')
        print(fs)
#         fs = glob.glob('./rttov_out_'+i+'_20220906/*')

        sid = []
        for i in fs:
            if len(glob.glob(i+'/*')) > 0:
                sid.append(i.split('/')[-1])

        pool = multiprocessing.Pool(processes=40)
        func=partial(write_bt_file, adj = o)
        result_list = list(pool.map(func, sid[:]))
        print(result_list)

