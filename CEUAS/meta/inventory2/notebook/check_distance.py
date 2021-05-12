""" Check distance between all the stations identified in a inventory to avoid having duplicates """

import pandas as pd
import os
import numpy
import math 
import numpy as np

def fdist(slat,slon,lat,lon):


    if type(slon) is str:
        if '.' in slon:
            flon=float(slon)
            flat=float(slat)
        elif ',' in slon:
            flon=float('.'.join(slon.split(',')))
            flat=float('.'.join(slat.split(',')))
        elif slon[-1] in 'WE':
            flon=float(slon[:3])+float(slon[4:7])/60+float(slon[-2:-1])/3600
            if slon[-1]=='W':
                flon=-flon
            flat=float(slat[:2])+float(slat[3:6])/60+float(slat[-2:-1])/3600
            if slat[-1]=='S':
                flat=-flat
        else:
            try:

                flon=float(slon[-9:-6])+float(slon[-5:-3])/60+float(slon[-2:])/3600
                if slon[0]=='-':
                    flon=-flon
                flat=float(slat[-8:-6])+float(slat[-5:-3])/60+float(slat[-2:])/3600
                if slat[0]=='-':
                    flat=-flat
            except:
                flon=numpy.nan
                flat=numpy.nan

    else:
        flon=slon
        flat=slat


    lats=numpy.append(flat,lat) # TODO changed this 
    lons=numpy.append(flon,lon)

    x=numpy.cos(lats*math.pi/180.)*numpy.cos(lons*math.pi/180.)
    y=numpy.cos(lats*math.pi/180.)*numpy.sin(lons*math.pi/180.)
    z=numpy.sin(lats*math.pi/180.)

    dist=numpy.arccos((x[0]*x[1:]+y[0]*y[1:]+z[0]*z[1:])*0.9999999)

    return dist



def load_df(dataset):
    """ Reading the station_configuration as a dataframe """
    ds = '../code/output_data/' + dataset 
    df = pd.read_csv(ds + '_meta.csv', sep = '\t')
    df = df[['primary_id', 'secondary_id', 'latitude', 'longitude', 'station_name']]
    return df


#for i in ['1','2','1759', '1761', '3188', 'igra2', 'bufr', 'ncar']:
    
for i in ['1','2','1759', '1761', '3188', 'bufr', 'ncar']:
    print("Dataset: " , i)
    df = load_df(i)
    
    for s, num in zip(df['primary_id'], range(len(df)) ):
        if '0-20' not in s:
            continue
        
        tested_stat = df.loc[num]        
        lat, lon = tested_stat.latitude , tested_stat.longitude
        
        different_stations = df.loc [df['primary_id'] != s ]

        dists = fdist(lat, lon, different_stations['latitude'], different_stations['longitude'])*180./math.pi # !! DONOT invert the order of the lats and longs station / inventory !!! 
        idy=numpy.argsort(dists)

        for y in idy:
            if dists[y] > 0.:
                #print('OK!')
                continue
                
            else:
                stat = different_stations.loc[y+2]
                dist_check = fdist(lat, lon, stat['latitude'], stat['longitude'] )*180./math.pi
                
                print('Found a station too close!',  dists[y], ' ' , dist_check, '\n', stat , '\n' , tested_stat  )
                print(0)