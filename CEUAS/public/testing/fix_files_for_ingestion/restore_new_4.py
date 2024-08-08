import h5py
import numpy
import numpy as np
import pandas
import pandas as pd
import glob
import os
# import harvest_convert_to_netCDF_yearSplit as harvest	
import ray

long_string_len = 100  # maximum allowed length of strings in header and observations_table
fixed_string_len = 20  # maximum allowed length of strings in header and observations_table
id_string_length = 10 # maximum allowed length of strings for observation_id and report_id in header and observations_table


### nan for integer  number 
int_void = -2147483648


""" Possible variable types as listed in the CMD tables """
okinds={'varchar (pk)':np.dtype('|S' + str(fixed_string_len) ),
        'varchar':np.dtype('|S' + str(fixed_string_len) ), 
               'numeric':np.float32, 
               'int':np.int32,
               #'timestamp with timezone':np.datetime64,
               'timestamp with timezone':np.int64,
               'int[]*':list,
               'int[]':list,
               'varchar[]*':list,
               'varchar[]':list}

""" Variable types to be used in the compressed netCDF files """
kinds={'varchar (pk)':str,
       'varchar(pk)':str,
             'varchar ':str,
             'varchar':str,
             ' varchar':str,             
             'numeric':np.float32,
             'numeric ':np.float32,             
             'int':np.int32,
             'int ':np.int32,             
             'int(pk)' : np.int32,
             'timestamp with timezone':np.float32,             
             'int[]*':list,
             'int[]':list,
             'varchar[]*':list,
             'varchar[]':list}

gkinds={'varchar (pk)':numpy.dtype('|S'+ str(fixed_string_len)  ),
        'varchar':numpy.dtype('|S'+ str(fixed_string_len)  ),
               'numeric':numpy.float32,

               'int':numpy.int32,
               'timestamp with timezone':numpy.datetime64,
               'int[]*':numpy.int32,
               'int[]':numpy.int32,
               'varchar[]*':numpy.dtype('|S'+ str(long_string_len)  ),
               'varchar[]':numpy.dtype('|S'+ str(long_string_len)  )}



def write_dict_h5(dfile, f, k, fbencodings, var_selection=[], mode='a', attrs={}): 
    """ Writes each separate variable from the observation or feedback tables inot netcdf using h5py.
          f is a pandas dataframe with one column, one for each variable
          k is either 'era5fb' or 'observations_table'
          fbencodings is the encodings of variable types, e.g. {'observations_id': { 'compression': 'gzip' } }
    """

    #attrs=  {'date_time':('units','seconds since 1900-01-01 00:00:00')}
    #attrs = {'observation_id': ('description', 'unique ID for observation'), 'report_id': ('description', 'Link to header information') , 'date_time':('units','seconds since 1900-01-01 00:00:00') }

    with h5py.File(dfile,mode) as fd:
        try:
            fd.create_group(k)
            index=numpy.zeros (f[list(f.keys())[0]].shape[0], dtype='S1')
            fd[k].create_dataset('index', data=index)
        except:
            pass
        if not var_selection:
            var_selection=list(f.keys())

        string10=numpy.zeros(fixed_string_len,dtype='S1')
        sdict={}
        slist=[]

        for v in var_selection:
            if v in fd[k].keys():
                continue

            #if v == 'source_file':
            #    a = 0
            #if v in [ 'report_event1@hdr' , 'report_rdbflag@hdr' , 'datum_anflag@body', 'datum_event1@body', 'datum_rdbflag@body', 'index' , 'varbc_ix@body']:
            #    a=0
                #continue 
            if v in [ 'index']:
                    continue
                    
            if type(f[v]) == pd.core.series.Series:
                fvv=f[v].values
            else:
                fvv=f[v]

            try:
                if fvv.dtype ==pd.Int64Dtype(): ### TO DO 
                    continue
            except:
                pass

            # this is needed to remove the <NA> types from pandas in the array or it crashes, due to type Int64 NAN in pandas 
            if v in [ 'report_event1@hdr' , 'report_rdbflag@hdr' , 'datum_anflag@body', 'datum_event1@body', 'datum_rdbflag@body', 'index' , 'varbc_ix@body']:
                fvv = np.array( [int_void if pd.isna(i) else i for i in fvv   ] )
            
            if type(fvv[0]) not in [str,bytes,numpy.bytes_]:  ### HORRIBLE HANDLING of types, dtypes, strings, bytes... 
                #print(v, '  ', type(fvv[0]) , '  ' , fvv.dtype )

                if fvv.dtype !='S1':
                    #if fvv.dtype == "Int64":
                    #    0
                    #vtype = np.int32
                    #else:
                    #    vtype = fvv.dtype
                    
                    try:
                        fd[k].create_dataset(v,fvv.shape,fvv.dtype,compression=fbencodings[v]['compression'], chunks=True)
                    except:
                        #fd[k].create_dataset(v,fvv.shape,'int32',compression=fbencodings[v]['compression'], chunks=True)  
                        fd[k].create_dataset(v,fvv.shape,fvv.dtype,compression='gzip', chunks=True)

                    try:
                        fd[k][v][:]=fvv[:]
                    except:
                        fd[k][v][:] = np.empty( (len( fvv)) )

                    if attrs:    #  attrs={'date_time':('units','seconds since 1900-01-01 00:00:00')}
                        try:
                            if v in attrs.keys():
                                for kk,vv in attrs[v].items():
                                    if type(vv) is str:  
                                        fd[k][v].attrs[kk]=numpy.bytes_(vv)
                                    else:
                                        fd[k][v].attrs[kk]=vv
                            else:
                                for kk in attrs.keys():
                                    if kk == 'type':
                                        continue
                                    vv = attrs[kk]
                                    fd[k][v].attrs[kk]=vv

                        except:
                            for kk in attrs.keys():
                                if kk == 'type':
                                    continue
                                vv = attrs[kk]
                                fd[k][v].attrs[kk]=vv

                    if v in ['date_time','report_timestamp','record_timestamp']:
                        fd[k][v].attrs['units']=numpy.bytes_('seconds since 1900-01-01 00:00:00')                            #print (  fk, ' ' , v , ' ' ,   ) 

                else:
                    fd[k].create_dataset(v,fvv.shape,fvv.dtype,compression=fbencodings[v]['compression'], chunks=True)
                    fd[k][v][:]=fvv[:]
                    slen=fvv.shape[1]
                    sdict[v]=slen
                    if slen not in slist:
                        slist.append(slen)
                        try:
                            fd[k].create_dataset( 'string{}'.format(slen),  data=string10[:slen]  )
                        except:
                            pass               
                    if v in attrs.keys():
                        fd[k][v].attrs['description']=numpy.bytes_(attrs[v]['description'])
                        fd[k][v].attrs['external_table']=numpy.bytes_(attrs[v]['external_table'])

            else:
                sleno=len(fvv[0])
                slen=sleno
                try:
                    slen=int(fvv.dtype.descr[0][1].split('S')[1])
                except:  
                    slen=15

                sdict[v]=slen
                if slen not in slist:
                    slist.append(slen)
                    try:
                        fd[k].create_dataset( 'string{}'.format(slen),  data=string10[:slen]  )
                    except:
                        pass               
                try:

                    fd[k].create_dataset(v,data=fvv.view('S1').reshape(fvv.shape[0],slen),compression=fbencodings[v]['compression'],chunks=True)
                except KeyError:
                    fd[k].create_dataset(v,data=fvv.view('S1').reshape(fvv.shape[0],slen),compression='gzip',chunks=True)

                # if v in attrs.keys():
                #     fd[k][v].attrs['description']     =numpy.bytes_(attrs[v]['description'])
                #     fd[k][v].attrs['external_table']=numpy.bytes_(attrs[v]['external_table'])       
                if attrs:    #  attrs={'date_time':('units','seconds since 1900-01-01 00:00:00')}
                        try:
                            if v in attrs.keys():
                                for kk,vv in attrs[v].items():
                                    if type(vv) is str:  
                                        fd[k][v].attrs[kk]=numpy.bytes_(vv)
                                    else:
                                        fd[k][v].attrs[kk]=vv
                            else:
                                for kk in attrs.keys():
                                    if kk == 'type':
                                        continue
                                    vv = attrs[kk]
                                    fd[k][v].attrs[kk]=vv
                        except:
                            for kk in attrs.keys():
                                if kk == 'type':
                                    continue
                                vv = attrs[kk]
                                fd[k][v].attrs[kk]=vv


            #variables_dic[v] = f[v].values.dtype

        for v in fd[k].keys(): #var_selection:
            l=0      
            if 'string'  in v or v== 'index' :                    
                continue 
            try:
                if type(f[v]) == pd.core.series.Series:
                    fvv=f[v].values
                else:
                    fvv=f[v]
                fd[k][v].dims[l].attach_scale(fd[k]['index'])
                #print(v,fvv.ndim,type(fvv[0]))
                if fvv.ndim==2 or type(fvv[0]) in [str,bytes,numpy.bytes_]:
                    slen=sdict[v]
                    #slen=10
                    fd[k][v].dims[1].attach_scale(fd[k]['string{}'.format(slen)])
            except:
                pass

        #i=4        
        for v in slist:
            s='string{}'.format(v)
            for a in ['NAME']:
                fd[k][s].attrs[a]=numpy.bytes_('This is a netCDF dimension but not a netCDF variable.')
            #i+=1

    return

# # @ray.remote
# def redo_file(file, outdir = '/data/public/converted_v19/'):
#     import hdf5plugin
#     print(file)

#     targetfile = outdir + '/' + file.split('/')[-1]
#     # if os.path.isfile(targetfile):
#     #     return
#     # try: 
#     with h5py.File(file, 'r') as fi:
#         for table in fi.keys():
#             if table == 'station_configuration':
#                 sc_len = len(fi[table]['elevation'])
#                 ht_len = len(fi['header_table']['latitude'][:])
#                 if sc_len == 1:
#                     print('!')
#                 for var in fi[table].keys():
#                     attrs = {}
#                     if len(fi[table][var].attrs.keys()) > 0:
#                         attrs[var] = {}
#                         for attr in fi[table][var].attrs:
#                             if attr != 'DIMENSION_LIST':
#                                 attrs[var][attr] = fi[table][var].attrs[attr]

#                 # else:
#                     # for var in fi[table].keys():
#                     #     attrs = {}
#                     #     if len(fi[table][var].attrs.keys()) > 0:
#                     #         attrs[var] = {}
#                     #         for attr in fi[table][var].attrs:
#                     #             if attr != 'DIMENSION_LIST':
#                     #                 attrs[var][attr] = fi[table][var].attrs[attr]

#                     #     data = fi[table][var][:]
#                     #     if data.ndim > 1: # decode is quite slow -> faster possible?
#                     #         data = [''.join([b.decode() for b in a]) for a in data]
#                     #     df_to_write = pd.DataFrame({var: data})
#                     #     write_dict_h5(targetfile, df_to_write, table, {var: { 'compression': 'gzip'} }, [var], attrs=attrs) #  ,'compression_opts': 4
                    
#         file_path = './done/' + file.split('/')[-1] + '.txt'
#         with open(file_path, "w") as text_file:
#             text_file.write('done')
#     # except:
#     #     file_path = './failed/' + file.split('/')[-1] + '.txt'
#     #     with open(file_path, "w") as text_file:
#     #         text_file.write('failed')


# files = glob.glob('/mnt/users/scratch/leo/scratch/converted_v19/long/*89056*.nc')# data/public/converted_v19_irregular/*.nc')

# # ray.init(num_cpus=20)
# #print(files)
# # result_ids = []
# # for file in files:
# #     result_ids.append(redo_file.remote(file))

# # results = ray.get(result_ids)
# # ray.shutdown()

# for file in files:
#     redo_file(file, outdir='.')

def find_adjacent_duplicates_indices(arr):
    indices = []
    for i in range(len(arr) - 1):
        if arr[i] == arr[i + 1]:
            indices.append(i)
    return indices

@ray.remote
def redo_file(file, attrs_encodings, outdir = '/data/public/converted_v19_dtype'):
    import hdf5plugin
    print(file)
    #{'observations_id': { 'compression': 'gzip' } ,'compression_opts': 4 }

    targetfile = outdir + '/' + file.split('/')[-1]
    if os.path.isfile('./done_4/' + file.split('/')[-1] + '.txt'):
        return
    # try: 
    with h5py.File(file, 'r') as fi:
        for table in fi.keys():
            print(table)
            if table == 'station_configuration':
                ht_len = len(fi['header_table']['latitude'][:])
                for var in fi[table].keys():
                    try:
                        attrs = attrs_encodings[table][var]
                    except:
                        attrs = {}
                        if len(fi[table][var].attrs.keys()) > 0:
                            attrs[var] = {}
                            for attr in fi[table][var].attrs:
                                if attr != 'DIMENSION_LIST':
                                    attrs[var][attr] = fi[table][var].attrs[attr]
                    if var == 'elevation':
                        data = fi['header_table']['height_of_station_above_sea_level'][:]
                    elif var in ['latitude', 'longitude']:
                        if ('0-20999-0' in file):
                            if ('107' in list(fi['recordindices'].keys())) and ('126' in list(fi['recordindices'].keys())):
                                wind_idx = fi['recordindices']['107'][1:] -1 # get all the indices of the first elements -> take the last and do -1
                                temp_idx = fi['recordindices']['126'][1:] -1 # get all the indices of the first elements -> take the last and do -1
                                if np.all(np.array(fi['observations_table'][var][:])[temp_idx] == np.array(fi['observations_table'][var][:])[wind_idx]):
                                    data = np.array(fi['observations_table'][var][:])[temp_idx]
                                else:
                                    wd = np.array(fi['observations_table'][var][:])[wind_idx]
                                    td = np.array(fi['observations_table'][var][:])[temp_idx]
                                    fadi = find_adjacent_duplicates_indices(wd) 
                                    if len(fadi) >= 1:
                                        swap_idx = np.array(fadi) + 1
                                        wd[swap_idx] = td[swap_idx]
                                    data = wd
                            elif ('107' in list(fi['recordindices'].keys())):
                                wind_idx = fi['recordindices']['107'][1:] -1
                                data = np.array(fi['observations_table'][var][:])[wind_idx]
                            elif ('126' in list(fi['recordindices'].keys())):
                                temp_idx = fi['recordindices']['126'][1:] -1
                                data = np.array(fi['observations_table'][var][:])[temp_idx]
                        else:
                            data = fi['header_table'][var][:]
                    elif var == 'record_number':
                        data = fi['header_table']['station_record_number'][:]
                    elif var == 'primary_id':
                        data = fi['header_table']['primary_station_id'][:]
                    else:
                        read_data = fi[table][var][0]
                        if read_data.ndim > 0: # decode is quite slow -> faster possible?
                            data = np.array([b''.join(read_data)] * ht_len).astype('|S'+str(len(read_data)))
                        else:
                            data = np.array([read_data] * ht_len)

                    try:
                        df_to_write = {var: data} # .astype(attrs_encodings[table][var]['type'])}
                    except:
                        df_to_write = {var: data} # pd.DataFrame({var: data})
                    write_dict_h5(targetfile, df_to_write, table, {var: { 'compression': 'gzip'} }, [var], attrs=attrs) #  ,'compression_opts': 4

            else:
                for var in fi[table].keys():
                    try:
                        attrs = attrs_encodings[table][var]
                    except:
                        attrs={}
                        if len(fi[table][var].attrs.keys()) > 0:
                            attrs[var] = {}
                            for attr in fi[table][var].attrs:
                                if attr != 'DIMENSION_LIST':
                                    attrs[var][attr] = fi[table][var].attrs[attr]

                    if ('0-20999-0' in file) and (table == 'header_table') and (var in ['latitude', 'longitude']):
                        if ('107' in list(fi['recordindices'].keys())) and ('126' in list(fi['recordindices'].keys())):
                            wind_idx = fi['recordindices']['107'][1:] -1 # get all the indices of the first elements -> take the last and do -1
                            temp_idx = fi['recordindices']['126'][1:] -1 # get all the indices of the first elements -> take the last and do -1
                            if np.all(np.array(fi['observations_table'][var][:])[temp_idx] == np.array(fi['observations_table'][var][:])[wind_idx]):
                                data = np.array(fi['observations_table'][var][:])[temp_idx]
                            else:
                                wd = np.array(fi['observations_table'][var][:])[wind_idx]
                                td = np.array(fi['observations_table'][var][:])[temp_idx]
                                fadi = find_adjacent_duplicates_indices(wd) 
                                if len(fadi) >= 1:
                                    swap_idx = np.array(fadi) + 1
                                    wd[swap_idx] = td[swap_idx]
                                data = wd
                        elif ('107' in list(fi['recordindices'].keys())):
                            wind_idx = fi['recordindices']['107'][1:] -1
                            data = np.array(fi['observations_table'][var][:])[wind_idx]
                        elif ('126' in list(fi['recordindices'].keys())):
                            temp_idx = fi['recordindices']['126'][1:] -1
                            data = np.array(fi['observations_table'][var][:])[temp_idx]
                    else:
                        data = fi[table][var][:]
                    if data.ndim > 1: # decode is quite slow -> faster possible?
                        data = np.array([b''.join(bytestri) for bytestri in data]).astype('|S'+str(len(data[0])))
                    try:
                        df_to_write = {var: data} #.astype(attrs_encodings[table][var]['type'])}
                    except:
                        df_to_write = {var: data} # pd.DataFrame({var: data})
                    write_dict_h5(targetfile, df_to_write, table, {var: { 'compression': 'gzip'} }, [var], attrs=attrs) #  ,'compression_opts': 4

        table = 'advanced_uncertainty'
        ot_len = len(fi['observations_table']['z_coordinate'][:])
        attrs = {}
        for var in ['180', '90', '60','30']:
            description = "Desroziers uncertainty v 1.0 - " + str(var) + " days window"
        
            var = 'desroziers_' + str(var)
            attrs[var] = {}
            attrs[var]['description'] = description

            data = np.array([np.nan] * ot_len)
            try:
                df_to_write = {var: data} # .astype(attrs_encodings[table][var]['type'])}
            except:
                df_to_write = {var: data} # pd.DataFrame({var: data})
            write_dict_h5(targetfile, df_to_write, table, {var: { 'compression': 'gzip'} }, [var], attrs=attrs) #  ,'compression_opts': 4

        for var in ['uncertainty_method', 'uncertainty_type']:
            
            attrs[var] = {}
            attrs[var]['description'] = 'General information about the uncertainty.'

            data = np.array([5] * ot_len)
            try:
                df_to_write = {var: data} # .astype(attrs_encodings[table][var]['type'])}
            except:
                df_to_write = {var: data} # pd.DataFrame({var: data})
            write_dict_h5(targetfile, df_to_write, table, {var: { 'compression': 'gzip'} }, [var], attrs=attrs)

    file_path = './done_4/' + file.split('/')[-1] + '.txt'
    with open(file_path, "w") as text_file:
        text_file.write('done')


files = glob.glob('/data/public/converted_v19_irregular/*00-0-42314*.nc') # data/public/converted_v19_irregular/*.nc')
attrs_encodings = np.load('dic_type_attributes.npy', allow_pickle=True).item()


ray.init(num_cpus=15)# , object_store_memory=(10**9)*65)

#files = glob.glob('/data/public/converted_v19_irregular/*00-0-26258*.nc')
#print(files)
result_ids = []
for file in files:
    result_ids.append(redo_file.remote(file, attrs_encodings))

results = ray.get(result_ids)
ray.shutdown()
