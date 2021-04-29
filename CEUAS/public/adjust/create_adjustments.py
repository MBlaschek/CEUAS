#!/usr/bin/env python
# By Ulrich Voggenberger
# Purpose: Create and compare adjustments with those from the CDS 
import logging
import pandas as pd
import numpy as np
import sys, zipfile, os, time
import matplotlib.pyplot as plt
import glob
import datetime
import urllib3
import cdsapi
# import xarray
import shutil
import warnings
import pickle
import h5py
# import netCDF4
sys.path.append(os.getcwd()+'/../cds-backend/code/')
import cds_eua3 as eua

def compare_adj(station):
    """Compare created adjustment to downloaded adjustment
    """
    file = ('Temperature_adjustment/0'+station+'/feedbackglobbincorrsave0'+station+'.nc')
    data = eua.CDMDataset(file)
    
    f = open('Temperature_adjustment/0'+station+'/found_breaks0'+station, "r")
    breaks = f.read()
    breaksfound = 0
    breakdates = []
    while breaksfound >= 0:
        breaksfound = breaks.find('break   ', breaksfound) + 1
        if breaksfound == 0:
            break
        breakdates.append(breaks[breaksfound+7:breaksfound+15])
    press = data.press[:]
    nightadj = data.rasocorr[0,:]
    dayadj = data.rasocorr[1,:]
    
    c = cdsapi.Client()
    r = c.retrieve('insitu-comprehensive-upper-air-observation-network',
                   {'variable': 'temperature',
                    'optional':['obs_minus_bg','bias_estimate','RISE_bias_estimate', 'RICH_bias_estimate', 'RASE_bias_estimate', 'RAOBCORE_bias_estimate'],
                    'statid': station,
                    'pressure_level':[10,20,30,50,70,100,150,200,250,300,400,500,700,850,925,1000]
                   }
                  )
    r.download(target='download.zip')
    assert os.stat('download.zip').st_size == r.content_length, "Downloaded file is incomplete"
    z = zipfile.ZipFile('download.zip')
    z.extractall(path='./tocompare/cds_'+station)
    z.close()
    cdsfile = glob.glob('./tocompare/cds_'+ station +'/*.nc')
    cdsdata=eua.CDMDataset(cdsfile[0]).to_dataframe()
    daydata = cdsdata[cdsdata.time.dt.hour > 6][cdsdata.time.dt.hour <= 18]
    nightdata = cdsdata[cdsdata.time.dt.hour <= 6].append(cdsdata[cdsdata.time.dt.hour > 18])
    print(str(cdsdata.time.iloc[0]))#.year)+str(cdsdata.time.iloc[0].month)+str(cdsdata.time.iloc[0].day))
    print(breakdates)
    
    breakdates.append(str(cdsdata.time.iloc[0]))#.year)+str(cdsdata.time.iloc[0].month)+str(cdsdata.time.iloc[0].day))
    dates = breakdates
    dates.reverse()
    print(dates)
    
    for i in range(len(dates)-1):
        print(dates[i])
        for j in range(len(press)):
            try:
                print(press[j])
                d_cdsd = daydata[daydata.plev == press[j]*100][daydata.time >= dates[i]][daydata.time < dates[i+1]]
                d_adjd = dayadj[j][i]
                print('day')
                print('cds',d_cdsd.RAOBCORE_bias_estimate.iloc[0])
                print('calc', d_adjd)

                n_cdsd = nightdata[nightdata.plev == press[j]*100][nightdata.time >= dates[i]][nightdata.time < dates[i+1]]
                n_adjd = nightadj[j][i]
                print('night')
                print('cds', n_cdsd.RAOBCORE_bias_estimate.iloc[0])
                print('calc', n_adjd)
                print('')
#                 if d_adjd == 0.872276 or n_adjd == 0.872276:
#                     print('+++++++++++++++++++++')
#                     print(j, i)
            except:
                pass
        print('')
        print('')
        print('')
    return 0

def plot_adj(station):
    """Compare created adjustment to downloaded adjustment
    """
    file = ('Temperature_adjustment/0'+station+'/feedbackglobbincorrsave0'+station+'.nc')
    data = eua.CDMDataset(file)
    
    f = open('Temperature_adjustment/0'+station+'/found_breaks0'+station, "r")
    breaks = f.read()
    breaksfound = 0
    breakdates = []
    while breaksfound >= 0:
        breaksfound = breaks.find('break   ', breaksfound) + 1
        if breaksfound == 0:
            break
        breakdates.append(breaks[breaksfound+7:breaksfound+15])
    press = data.press[:]
    nightadj = data.rasocorr[0,:]
    dayadj = data.rasocorr[1,:]
    
    c = cdsapi.Client()
    r = c.retrieve('insitu-comprehensive-upper-air-observation-network',
                   {'variable': 'temperature',
                    'optional':['obs_minus_bg','bias_estimate','RISE_bias_estimate', 'RICH_bias_estimate', 'RASE_bias_estimate', 'RAOBCORE_bias_estimate'],
                    'statid': station,
                    'pressure_level':[10,20,30,50,70,100,150,200,250,300,400,500,700,850,925,1000]
                   }
                  )
    r.download(target='download.zip')
    assert os.stat('download.zip').st_size == r.content_length, "Downloaded file is incomplete"
    z = zipfile.ZipFile('download.zip')
    z.extractall(path='./tocompare/cds_'+station)
    z.close()
    cdsfile = glob.glob('./tocompare/cds_'+ station +'/*.nc')
    cdsdata=eua.CDMDataset(cdsfile[0]).to_dataframe()
    daydata = cdsdata[cdsdata.time.dt.hour > 6][cdsdata.time.dt.hour <= 18]
    nightdata = cdsdata[cdsdata.time.dt.hour <= 6].append(cdsdata[cdsdata.time.dt.hour > 18])
    print(str(cdsdata.time.iloc[0]))#.year)+str(cdsdata.time.iloc[0].month)+str(cdsdata.time.iloc[0].day))
    print(breakdates)
    
    breakdates.append(str(cdsdata.time.iloc[0]))#.year)+str(cdsdata.time.iloc[0].month)+str(cdsdata.time.iloc[0].day))
    dates = breakdates
    dates.reverse()
    print(dates)
    
    for i in range(len(dates)-1):
        print(dates[i])
        for j in range(len(press)):
            try:
                print(press[j])
                d_cdsd = daydata[daydata.plev == press[j]*100][daydata.time >= dates[i]][daydata.time < dates[i+1]]
                d_adjd = dayadj[j][i]
                print('day')
                print('cds',d_cdsd.RAOBCORE_bias_estimate.iloc[0])
                print('calc', d_adjd)

                n_cdsd = nightdata[nightdata.plev == press[j]*100][nightdata.time >= dates[i]][nightdata.time < dates[i+1]]
                n_adjd = nightadj[j][i]
                print('night')
                print('cds', n_cdsd.RAOBCORE_bias_estimate.iloc[0])
                print('calc', n_adjd)
                print('')
#                 if d_adjd == 0.872276 or n_adjd == 0.872276:
#                     print('+++++++++++++++++++++')
#                     print(j, i)
            except:
                pass
        print('')
        print('')
        print('')
    return 0

def plt_adj(dates, cdata, adata, plev, station, time):
    dates_list = [datetime.datetime.strptime(date, '%Y%m%d').date() for date in dates[1:]]
    fig, ax = plt.subplots()
    ax.scatter(dates_list, cdata[str(plev)], label='CDS')
    ax.scatter(dates_list, adata[str(plev)], label='Calc')

    ax.set(xlabel='breaks', ylabel='Adjustment',
           title= str(station)+' Temperature '+str(plev)+' hPa '+time)
    ax.legend()
    ax.grid()

    fig.savefig(str(station)+'_Temperature_'+str(plev)+'_hPa_'+time+'.png')
    plt.show()
    return 0

def get_adj(station):
    import datetime
    """Compare created adjustment to downloaded adjustment
    """
    file = ('/home/srvx11/lehre/users/a1400070/CEUAS/CEUAS/public/adjust/Temperature_adjustment/0'+station+'/feedbackglobbincorrsave0'+station+'.nc')
    data = eua.CDMDataset(file)
    
    f = open('/home/srvx11/lehre/users/a1400070/CEUAS/CEUAS/public/adjust/Temperature_adjustment/0'+station+'/found_breaks0'+station, "r")
    breaks = f.read()
    breaksfound = 0
    breakdates = []
    while breaksfound >= 0:
        breaksfound = breaks.find('break   ', breaksfound) + 1
        if breaksfound == 0:
            break
        breakdates.append(breaks[breaksfound+7:breaksfound+15])
    press = data.press[:]
    nightadj = data.rasocorr[0,:]
    dayadj = data.rasocorr[1,:]
    
    c = cdsapi.Client()
    r = c.retrieve('insitu-comprehensive-upper-air-observation-network',
                   {'variable': 'temperature',
                    'optional':['obs_minus_bg','bias_estimate','RISE_bias_estimate', 'RICH_bias_estimate', 'RASE_bias_estimate', 'RAOBCORE_bias_estimate'],
                    'statid': station,
                    'pressure_level':[10,20,30,50,70,100,150,200,250,300,400,500,700,850,925,1000]
                   }
                  )
    r.download(target='download.zip')
    assert os.stat('download.zip').st_size == r.content_length, "Downloaded file is incomplete"
    z = zipfile.ZipFile('download.zip')
    z.extractall(path='./tocompare/cds_'+station)
    z.close()
    cdsfile = glob.glob('./tocompare/cds_'+ station +'/*.nc')
    cdsdata=eua.CDMDataset(cdsfile[0]).to_dataframe()
    daydata = cdsdata[cdsdata.time.dt.hour > 6][cdsdata.time.dt.hour <= 18]
    nightdata = cdsdata[cdsdata.time.dt.hour <= 6].append(cdsdata[cdsdata.time.dt.hour > 18])
    print(str(cdsdata.time.iloc[0]))#.year)+str(cdsdata.time.iloc[0].month)+str(cdsdata.time.iloc[0].day))
    print(breakdates)
    
    breakdates.append(str(cdsdata.time.iloc[0]))#.year)+str(cdsdata.time.iloc[0].month)+str(cdsdata.time.iloc[0].day))
    dates = breakdates
    dates.reverse()
    print(dates)
    
    calc_day_out = {}
    calc_night_out = {}
    cds_day_out = {}
    cds_night_out = {}
    for i in [10,20,30,50,70,100,150,200,250,300,400,500,700,850,925,1000]:
        calc_day_out[str(i)]=[]
        calc_night_out[str(i)]=[]
        cds_day_out[str(i)]=[]
        cds_night_out[str(i)]=[]
    for i in range(len(dates)-1):
        for j in range(len(press)):
            try:
                d_adjd = dayadj[j][i]
                calc_day_out[str(int(press[j]))].append(d_adjd)
            except:
                calc_day_out[str(int(press[j]))].append(np.nan) 
            try:
                d_cdsd = daydata[daydata.plev == press[j]*100][daydata.time >= dates[i]][daydata.time < dates[i+1]]
                cds_day_out[str(int(press[j]))].append(d_cdsd.RAOBCORE_bias_estimate.iloc[0])
            except:
                cds_day_out[str(int(press[j]))].append(np.nan)
            try:
                n_cdsd = nightdata[nightdata.plev == press[j]*100][nightdata.time >= dates[i]][nightdata.time < dates[i+1]]
                cds_night_out[str(int(press[j]))].append(n_cdsd.RAOBCORE_bias_estimate.iloc[0])
            except:
                cds_night_out[str(int(press[j]))].append(np.nan)
            try:
                n_adjd = nightadj[j][i]
                calc_night_out[str(int(press[j]))].append(n_adjd)
            except:
                calc_night_out[str(int(press[j]))].append(np.nan)
                
    return dates, calc_day_out, calc_night_out, cds_day_out, cds_night_out

if __name__ == "__main__":
#     print('preparing the environment')
#     os.system('bash prep.sh')
    
#     print('download data from cds for further processing')
#     os.system('bash dldata.sh')

#     print('create adjustments')
#     os.system('bash adjust.sh')

#     glob.glob
    compare_list = ['11035']
    HOME = os.getcwd()
    print(HOME)
    
    for i in compare_list:
        status = compare_adj(i)
        if status != 0:
            break
    
    if status == 0:
        print('---')
        print('comparison successful')
        print('---')
    else:
        sys.exit(status)
    
    