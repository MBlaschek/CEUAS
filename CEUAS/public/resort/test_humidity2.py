#!/usr/bin/env python
# coding: utf-8

# # Check humidity conversions and adjustments

# In[1]:


import os,sys,glob
import numpy as np
import h5py
import gzip
from datetime import datetime,timedelta
print(os.getcwd())
from convert_numbamod import Sonntag,sh2vap,vap2sh,svp,liquid,bisection,func,rms
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from numba import njit, prange, set_num_threads
from tqdm import tqdm
from functools import partial
from multiprocessing import Pool
import time
import pickle
import ray
sys.path.append(os.getcwd()+'/../common/Rasotools/rasotools')
from anomaly import janomaly, fastlinregress
from utils import tdist, tcost
#import rasotools

#get_ipython().run_line_magic('matplotlib', 'inline')


# In[2]:


WIGOS='0-20001-0-71957'
WIGOS='0-20000-0-72357'
WIGOS='0-20001-0-11035'
WIGOS='0-20000-0-12843'
WIGOS='0-20001-0-10393'
WIGOS='0-20000-0-72357'
WIGOS='0-20001-0-71957'
YEAR='2011'
RVERSION='25'


# In[3]:


paths=dict(hte5='/mnt/users/scratch/leo/scratch/era5/odbs/1/new/',
htigra='/mnt/scratch/scratch/federico/databases_service2/IGRA2_03012024/',
he5=f'/mnt/users/scratch/uvoggenberger/CUON_HARVEST/harvest_regular/era5_1/{WIGOS}',
higra=f'/mnt/users/scratch/uvoggenberger/CUON_HARVEST/harvest_regular/igra2/{WIGOS}',
#higra=f'/mnt/users/scratch/uvoggenberger/CUON_HARVEST/harvest_regular/ncar/{WIGOS}',
merged=f'/mnt/users/scratch/leo/scratch/UH/MERGED_YEARLY_20NOV2024_REGULAR/{WIGOS}/',
resorted=f'/mnt/users/scratch/leo/scratch/converted_v{RVERSION}/long/')

files=dict(hte5='',
htigra='',
he5=f'{WIGOS}_{YEAR}_*.gz.nc',
higra=f'{WIGOS}_{YEAR}_igra2_harvested_*-data.txt.nc',
merged=f'{WIGOS}_{YEAR}_CEUAS_merged_v3.nc',
resorted=f'{WIGOS}_CEUAS_merged_v3.nc')


@njit#(parallel=True)
def nnanmean(a, am,t, mingood, quant):
    #am = np.zeros((1, a.shape[0], a.shape[1]))
    #ac = np.zeros((a.shape[0], a.shape[1]), dtype=np.int32)
    for k in prange(a.shape[0]):
        for j in range(a.shape[1]):
            ac = 0
            for i in range(a.shape[2]):
                if a[k, j, i] == a[k, j, i]:
                    am[t, k, j] += a[k, j, i]
                    ac += 1
            if ac >= mingood:
                am[t, k, j] /= ac
            else:
                am[t, k, j] = np.nan
    return

@njit#(parallel=True)
def nnanrms(a, am,t, mingood, quant):
    #am = np.zeros((1, a.shape[0], a.shape[1]))
    #ac = np.zeros((a.shape[0], a.shape[1]), dtype=np.int32)
    for k in prange(a.shape[0]):
        for j in range(a.shape[1]):
            ac = 0
            for i in range(a.shape[2]):
                if a[k, j, i] == a[k, j, i]:
                    am[t, k, j] += a[k, j, i] * a[k, j, i]
                    ac += 1
            if ac >= mingood:
                am[t, k, j] = np.sqrt(am[t, k, j]/ac)
            else:
                am[t, k, j] = np.nan
    return

@njit#(parallel=True)
def nnanquant(a, am,t, mingood, quant=0.5):
    #am = np.zeros((1, a.shape[0], a.shape[1]))
    #ac = np.zeros((a.shape[0], a.shape[1]), dtype=np.int32)
    for k in prange(a.shape[0]):
        for j in range(a.shape[1]):
            mask = ~np.isnan(a[k, j, :])
            if np.sum(mask) > mingood:
                am[t, k, j] = np.quantile(np.abs(a[k, j, mask]), quant)
                #am[t, k, j] = np.quantile(a[k, j, mask], quant)
            else:
                am[t, k, j] = np.nan
                
    return

##dtn = np.arange(20, dtype=np.int64)
#plevs = np.arange(9)
#pparams = np.array(['dt', 'fg_dep','fg_off','fg_off_corr', 'adj', 'obs', 'obs_adj'])
#params = pparams[1:]
#ld = list(zip(pparams, (np.int64, np.float32, np.float32, np.float32, np.float32, np.float32, np.float32)))
#deps = np.empty((len(params), len(plevs), dtn.shape[0]), dtype=np.float32)
#ts = np.empty((1, len(params), len(plevs)), dtype=np.float32)
##deps['dt'][0, :] = dtn
#deps[:] =1.

#set_num_threads(len(params))
#x = nnanmean(deps, ts, 0, 30)
#print(x)
                
def rgb(r,g,b):
    return tuple(np.asarray([r,g,b],dtype=np.float32))

def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    se = [(None,) * 3, 0.0]
    for s in seq:
        se.append(s[0])
        se.append(s[1])#+ list(seq) +
        seq=se+[ (None,) * 3]
        cdict = {'red': [], 'green': [], 'blue': []}
        for i, item in enumerate(seq):
            if isinstance(item, float):
                r1, g1, b1 = seq[i - 1]
                r2, g2, b2 = seq[i + 1]
                cdict['red'].append([item, r1, r2])
                cdict['green'].append([item, g1, g2])
                cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)

rgblist=["rgb(0,0,0.3)", "rgb(0,0,0.5)",
"rgb(0,0,0.7)", "rgb(0,0,0.9)", "rgb(0,0.15,1)",
"rgb(0,0.3,1)", "rgb(0,0.45,1)", "rgb(0,0.6,1)",
"rgb(0,0.75,1)", "rgb(0,0.85,1)", "rgb(0.2,0.95,1)",
"rgb(0.45,1,1)", "rgb(0.75,1,1)", "rgb(1,1,0)",
"rgb(1,0.9,0)", "rgb(1,0.8,0)", "rgb(1,0.7,0)",
"rgb(1,0.6,0)", "rgb(1,0.5,0)", "rgb(1,0.4,0)",
"rgb(1,0.3,0)", "rgb(1,0.15,0)", "rgb(0.9,0,0)",
"rgb(0.7,0,0)", "rgb(0.5,0,0)", "rgb(0.3,0,0)"]
rgblist2=zip([eval(rgblist[l]) for l in range(len(rgblist))],np.linspace(0,1,len(rgblist)))

cmnew=make_colormap(rgblist2)

sdict=dict(
    rs90 = ['71', '72', '73', '74','78','V9a','V9t' ],
    rs92 = ['70', '79', '80', '81','114','VNT','VnG'  ],
    rs80 = ['37', '52', '60', '61', '62', '63' , '66', '67','VC_'  ],
    rs41 = ['123', '124', '125', '141', '142'],
    graw = ['17'],
    m10= ['177'],
    m20= ['163'],
    viz  = ['49','ZKB'],
    vschroeder = ['VN', 'VP'],
)
    
sdict_ref = ray.put(sdict)

def plt_trends(da):
    pltparams = {'legend.fontsize': 'x-large',
              'figure.figsize': (12, 8),
             'axes.labelsize': 'x-large',
             'axes.titlesize': 20,
             'xtick.labelsize':'medium',
             'ytick.labelsize':'medium'}
    plt.rcParams.update(pltparams)

    statnum = len(da['lon'])
    
    dists=np.zeros(np.max(((statnum+1)*statnum//2, 1)),np.float32)
    x,y,z=tdist(dists,da['lat'], da['lon'],1)
    scosts=np.zeros(statnum)
    cost=tcost(dists,da['trend'],scosts)

    plt.figure(figsize=(10, 7))
    ax = plt.axes([0.1,0.05,0.8,0.8], projection=ccrs.PlateCarree())
    ax.set_extent([-180,180, -90, 90], ccrs.PlateCarree())
    ax.add_feature(cfeature.OCEAN, zorder=0)
    ax.coastlines()

    reduced = np.asarray(da['fn'])
    cm = ax.scatter(da['lon'], da['lat'], s=40, alpha=1,
                c= da['trend'],
                cmap=cmnew,
                vmin=da['ranges'][da['plev']][da['v']][0],
                vmax=da['ranges'][da['plev']][da['v']][1],
                transform=ccrs.PlateCarree())
    ax.gridlines(draw_labels=True)
    plt.colorbar(cm, orientation='horizontal', label=f'Trend [{da['units']}/10a]', shrink=0.9, pad=0.05)
    #plt.tight_layout()
    os.makedirs('plots_new', exist_ok=True)
    pfname = f'{da['wpath']}/plots_new/{da['version']}/{da['vname'][da['v']]}_{da['param']}_trend_{da['interval'][0]}-{da['interval'][1]}_{da['plev']}.png'

    mtrend = np.full((18, 36), np.nan)
    i = 0
    for ilat in range(-85, 90, 10):
        j = 0
        for ilon in range(-175, 180, 10):
            idx = np.where((da['lat']>=ilat-5)& (da['lat']<ilat+5)&(da['lon']>=ilon-5)& (da['lon']<ilon+5))[0]
            if len(idx) > 0:
                mtrend[i, j] = np.nanmean(da['trend'][idx])
            j += 1
        i += 1
    zmtrend = np.nanmean(mtrend, axis=1)
    weight = np.cos(np.arange(-85, 90, 10)*np.pi/180.)
    mask = ~np.isnan(zmtrend)
    gtrend = np.mean(zmtrend[mask]*weight[mask]) / np.mean(weight[mask])
    print(f'gtrend: {gtrend:5.3f} [{da['units']}/10a], # {np.sum(mask)}, {da['plev']} hPa')
    plt.title(f'{da['vname'][da['v']]}-{da['param']}, {da['interval'][0]}-{da['interval'][1]}, {da['plev']} hPa'+'\n'\
              +'trend heterogeneity cost function: '+'{:.2f}'.format(cost)+'\n'+'# stations: '+str(statnum)+f' global mean trend: {gtrend:5.3f}, # {np.sum(mask)}', y=1.1)
        
    try:
        
        plt.savefig(pfname, bbox_inches='tight')
    except Exception as e:
        print('could not write', pfname)
    plt.close()

ray_plt_trends = ray.remote(plt_trends)

def plt_rmsdevs(da):
    pltparams = {'legend.fontsize': 'x-large',
              'figure.figsize': (12, 8),
             'axes.labelsize': 'x-large',
             'axes.titlesize': 20,
             'text.usetex' : False,
             'backend': 'agg',
             'xtick.labelsize':'medium',
             'ytick.labelsize':'medium'}
    plt.rcParams.update(pltparams)

    statnum = len(da['lon'])
    
    dists=np.zeros(np.max(((statnum+1)*statnum//2, 1)),np.float32)
    x,y,z=tdist(dists,da['lat'], da['lon'],1)
    scosts=np.zeros(statnum)
    cost=tcost(dists,da['rmsdevs'],scosts)
    qsv = np.arange(0.85, 1.01, 0.05)
    qs = np.nanquantile(da['rmsdevs'], qsv)
    #idx = np.argsort(da['rmsdevs'])
    ##da['ranks'] = idx / len(idx) *100
    if qs.size<qsv.size:
        return 
    
    flog = np.floor(np.log10(qs[3]))
    if True: #flog>=0:
        if da['title'] == 'mean':
            
            da['ranges'][da['plev']][da['v']][0] = -np.floor(qs[3]/(10**flog)) *10**flog
        else:
            da['ranges'][da['plev']][da['v']][0] = 0
            
        da['ranges'][da['plev']][da['v']][1] = np.floor(qs[3]/(10**flog)) *10**flog
    else:
        da['ranges'][da['plev']][da['v']][0] = -np.floor(qs[3]/(10**flog)) *10**flog
        da['ranges'][da['plev']][da['v']][1] = np.floor(qs[3]/(10**flog)) *10**flog
    
    #da['ranges'][da['plev']][da['v']][0] = 85
    #da['ranges'][da['plev']][da['v']][1] = 105
        

    plt.figure(figsize=(10, 7))
    ax = plt.axes([0.1,0.05,0.8,0.8], projection=ccrs.PlateCarree())
    ax.set_extent([-180,180, -90, 90], ccrs.PlateCarree())
    ax.add_feature(cfeature.OCEAN, zorder=0)
    ax.coastlines()

    reduced = np.asarray(da['fn'])
    cm = ax.scatter(da['lon'], da['lat'], s=40, alpha=1,
                c= da['rmsdevs'],
                cmap=cmnew,
                vmin=da['ranges'][da['plev']][da['v']][0],
                vmax=da['ranges'][da['plev']][da['v']][1],
                transform=ccrs.PlateCarree())
    ax.gridlines(draw_labels=True)
    plt.title(f'{da['vname'][da['v']]}-{da['reportype']}{da['title']}-{da['param']}, {da['interval'][0]}-{da['interval'][1]}, {da['plev']} hPa, {statnum} stations', y=1.1)
    cbar = plt.colorbar(cm, orientation='horizontal', label=f'{da['units']}', shrink=0.9, pad=0.1)
    #pos = cbar.ax.get_position()
    #ax2 = cbar.ax.twiny()
    #ax2.set_xlim([85, 105])
    
    ## resize the colorbar (otherwise it overlays the plot)
    #pos.y0 +=0.05
    #cbar.ax.set_position(pos)
    #ax2.set_position(pos)
    #ax2.set_xticklabels([f'{q / 10 ** flog:.2f}' for q in qs]+[f'RMS / 10^{int(flog)} {da['units']}'])
    #print(da['plev'])
    if da['plev'] == 700:
        xx = 0
    
    #plt.tight_layout()
    os.makedirs(f'{da['wpath']}/plots_rmsnew/{da['version']}', exist_ok=True)
    pfname = f'{da['wpath']}/plots_rmsnew/{da['version']}/{da['vname'][da['v']]}_{da['reportype']}{da['title']}_{da['param']}_{da['interval'][0]}-{da['interval'][1]}_{da['plev']}.png'

    try:
        
        plt.savefig(pfname, bbox_inches='tight')
    except Exception as e:
        print('could not write', pfname, e)
    plt.close()
    return

ray_plt_rmsdevs = ray.remote(plt_rmsdevs)

def plt_rmsdevquantiles(da):
    pltparams = {'legend.fontsize': 'x-large',
              'figure.figsize': (12, 8),
             'axes.labelsize': 'x-large',
             'axes.titlesize': 20,
             'text.usetex' : False,
             'backend': 'agg',
             'xtick.labelsize':'medium',
             'ytick.labelsize':'medium'}
    plt.rcParams.update(pltparams)

    statnum = len(da['lon'])
    
    #dists=np.zeros(np.max(((statnum+1)*statnum//2, 1)),np.float32)
    #x,y,z=tdist(dists,da['lat'], da['lon'],1)
    #scosts=np.zeros(statnum)
    #cost=tcost(dists,da['rmsdevs'],scosts)
    qsv = np.arange(0.85, 1.01, 0.025)
    qs = np.nanquantile(da['rmsdevs'], qsv)
    idx = np.argsort(da['rmsdevs'])
    #da['ranks'] = idx / len(idx) *100
    if qs.size<qsv.size:
        return qsv * np.nan
    
    flog = np.floor(np.log10(qs[3]))
    if np.isnan(flog):
        return qsv * np.nan
    
    if flog>=0:
        da['ranges'][da['plev']][da['v']][1] = np.floor(qs[3]/(10**flog)) *10**flog
    else:
        da['ranges'][da['plev']][da['v']][1] = np.floor(qs[3]/(10**flog)) *10**flog
    
    da['ranges'][da['plev']][da['v']][0] = 85
    da['ranges'][da['plev']][da['v']][1] = 105
        

    plt.figure(figsize=(10, 7))
    ax = plt.axes([0.1,0.05,0.8,0.8], projection=ccrs.PlateCarree())
    ax.set_extent([-180,180, -90, 90], ccrs.PlateCarree())
    ax.add_feature(cfeature.OCEAN, zorder=0)
    ax.coastlines()

    reduced = np.asarray(da['fn'])
    cm = ax.scatter(da['lon'][idx], da['lat'][idx], s=40, alpha=1,
                c= np.arange(idx.shape[0])/idx.shape[0]*100,
                cmap=cmnew,
                vmin=da['ranges'][da['plev']][da['v']][0],
                vmax=da['ranges'][da['plev']][da['v']][1],
                transform=ccrs.PlateCarree())
    ax.gridlines(draw_labels=True)
    plt.title(f'{da['vname'][da['v']]}-{da['reportype']}{da['title']}-{da['param']}, {da['interval'][0]}-{da['interval'][1]}, {da['plev']} hPa, {statnum} stations', y=1.1)
    cbar = plt.colorbar(cm, orientation='horizontal', label=f'Percentile [%]', shrink=0.9, pad=0.1)
    pos = cbar.ax.get_position()
    ax2 = cbar.ax.twiny()
    ax2.set_xlim([85, 105])
    
    # resize the colorbar (otherwise it overlays the plot)
    pos.y0 +=0.05
    cbar.ax.set_position(pos)
    ax2.set_position(pos)
    ax2.set_xticklabels([f'{q / 10 ** flog:.2f}' for q in qs] + [f'RMS / 10^{int(flog)} {da['units']}'])
    #print(da['plev'])
    if da['plev'] == 700:
        xx = 0
    
    #plt.tight_layout()
    os.makedirs(f'{da['wpath']}/plots_rmsnew/{da['version']}', exist_ok=True)
    pfname = f'{da['wpath']}/plots_rmsnew/{da['version']}/{da['vname'][da['v']]}_{da['reportype']}q{da['title']}_{da['param']}_{da['interval'][0]}-{da['interval'][1]}_{da['plev']}.png'

    try:
        
        plt.savefig(pfname, bbox_inches='tight')
    except Exception as e:
        print('could not write', pfname, e)
    plt.close()
    return

ray_plt_rmsdevquantiles = ray.remote(plt_rmsdevquantiles)

@njit(cache=True)
def touv(wd, ws):
    wd[np.isnan(ws)] = np.nan
    wd[ws < 0.001] = 0.
    ws[np.isnan(wd)] = np.nan
    u = ws * np.cos(np.radians(270.-wd))
    v = ws * np.sin(np.radians(270.-wd))
    return u, v, wd, ws

@njit(cache=True)
def wswd(u, v):

    wd = - np.arctan2(v, u) * 180 / np.pi - 90.
    wd[wd<0] += 360.
    return np.sqrt(u ** 2 + v ** 2), wd 
    

def read_adj(starts,plevs,factors, bedict, vtup,params,mingood,sdict,reportype, fn):
    tt = time.time()
    wigos=fn.split('/')[-1].split('_CEUAS')[0]
    with h5py.File(fn) as f:
        fdicts={'fn':wigos, 'lat':f['header_table']['latitude'][-1], 'lon':f['header_table']['longitude'][-1], 'sdict': sdict,'reportype': reportype,}
        #if (fdicts['lon']) < - 30.:
            #return None
        fri = f['recordindices']
        ril = list(fri.keys())
        istarts=np.searchsorted(fri['recordtimestamp'][:],starts)
        g = 0
        for i in range(len(istarts)-1): 
            if istarts[i+1] - istarts[i] > 0:
                fdicts[starts[i]] = {}
                g += 1
        if g < mingood:
            print(f'{os.path.basename(fn)}, {len(fdicts) - 3}, {time.time()-tt:5.4f}')
            return None
        
        fo = f['observations_table']
        if sdict is not None:
            rk = fri[ril[1]][[0, -1]]
            sondetype = fo['sensor_id'][rk[0]:rk[1]].view('S4')
            usondetype = np.unique(sondetype)
            st = []
            for sdk in sdict.keys():
                if 'rs' in sdk:
                    st = st + sdict[sdk]
            st = np.sort(np.array(st, dtype='S4'))

            grst = []
            for rst in usondetype:
                if rst in st:
                    grst.append(rst)
            if len(grst) == 0:
                return None
        
        foo = fo['observation_value']#[:]
        fe5 = f['era5fb']['fg_depar@body']
        foz = fo['z_coordinate']
        try:
            
            fah = f['advanced_homogenisation']
        except:
            print(fn, 'no homogeneity adjustments')
            fah = None
        #print(wigos)
        metadict = {}
        repdict = {}
        oris = {}
        orstarts = {}
        orstops = {}
        rstarts = {}
        rstops = {}
        for var in vtup: #'137','39', '138':#,'126','139', '140':
            
            try:
                oris[var]=fri[var][istarts[0]:istarts[-1] + 1]
                ori = oris[var]
                orstarts[var],orstops[var]=ori[[0,-1]]
            except:
                continue
            
            rstarts[var] = {}
            rstops[var] = {}

            ii = -1
            for istart, istop in zip(istarts[:-1], istarts[1:]):
                ii += 1
                if istop - istart == 0:
                    continue
                try:
#                    ri=f['recordindices'][var][istart:istop + 1]
                    ri=ori[istart-istarts[0]:istop-istarts[0] + 1]
                    if ri[-1] - ri[0] == 0:
                        continue
                    rstarts[var][istart],rstops[var][istart]=ri[[0,-1]]
                except:
                    continue

        for var in vtup: #'137','39', '138':#,'126','139', '140':
            
            try:
                #ori=fri[var][istarts[0]:istarts[-1] + 1]
                #orstart,orstop=ori[[0,-1]]
                ori = oris[var]
                orstart,orstop=orstarts[var],orstops[var]            
            except:
                continue
            
            oz=foz[orstart:orstop]
            odt=fo['date_time'][orstart:orstop]
            ofg_dep=fe5[orstart:orstop]
            ofg_off=f['era5fb']['fg_depar@offline'][orstart:orstop]
            oov = foo[orstart:orstop]
            
            if sdict is not None:
                ometa = fo['sensor_id'][orstart:orstop].view('S4')
            if reportype is not None:
                oreptype = f['era5fb']['reportype'][orstart:orstop]
            #for g in  ofg_dep, ofg_off:
                #for p in plevs:
                    #zmask = oz == p * 100
                    #if g[zmask].shape[0] > 30:
                        #x = np.nanpercentile(g[zmask], (0.1, 99.9))
                        #mask = (g<=x[0])|(g>=x[1])
                        #if np.sum(zmask&mask) > 0:
                            
                            #g[mask&zmask] = np.nan
                            
            try:
                
                obe = fah[bedict[var]][orstart:orstop]
            except:
                obe = np.zeros(oz.shape[0])
            #whstd = np.nanstd(obe)
            #if var == '106' and whstd > 0:
                
                #print('windhom:', fn, np.nanstd(obe))

            ii = -1
            tri = None
            for istart, istop in zip(istarts[:-1], istarts[1:]):
                ii += 1
                if istop - istart == 0:
                    continue
                try:
#                    ri=f['recordindices'][var][istart:istop + 1]
                    ri=ori[istart-istarts[0]:istop-istarts[0] + 1]
                    #rstart,rstop=ri[[0,-1]]
                    #if rstop - rstart == 0:
                        #continue
                    rstart, rstop = rstarts[var][istart], rstops[var][istart]
                except:
                    continue
                # check if records have particular sonde type
                #ttt = time.time()
                if istart not in metadict.keys() and sdict is not None:
                    meta = ometa[rstart-orstart:rstop-orstart][:, 0]
                    imeta = np.searchsorted(st, meta)
                    imeta[imeta==st.shape[0]] = st.shape[0] - 1
                    metadict[istart] = np.sum(meta==st[imeta]) > len(imeta) *0.9
                else:
                    x = 0
                if var in ('106', '107', '139', '140') and reportype is not None:
                    if '126' in fri.keys():
                        if tri is None:
                            tri = fri['126'][istarts[0]:istarts[-1] + 1]
                        tgood = tri[istart-istarts[0]:istop-istarts[0]+1]
                        repdict[istart] = np.ones(rstop-rstart, dtype=bool)
                        for i in range(istart, istop):
                            if tgood[i-istart+1] - tgood[i-istart] > 0:
                                #print(rep[ri[i-istart] - orstart:ri[i-istart+1] - orstart])
                                #print(odt[ri[i-istart] - orstart:ri[i-istart+1] - orstart])
                                repdict[istart][ri[i-istart] - rstart:ri[i-istart+1] - rstart] = False
                        if 16022 in reportype:
                            repdict[istart] = ~repdict[istart]
                            
                    #rep = oreptype[rstart-orstart:rstop-orstart]
                    #irep = np.searchsorted(reportype, rep)
                    #irep[irep==reportype.shape[0]] -= 1 
                    #repdict[istart] = rep==reportype[irep]
                else:
                    repdict[istart] = np.ones(rstop-rstart, dtype=bool)
                    
                #print('meta', time.time()-ttt)
                if sdict is not None and not metadict[istart]:
                    continue
                
                
                #if reportype is not None and not repdict[istart]:
                    #continue

                z=oz[rstart-orstart:rstop-orstart]
                dt=odt[rstart-orstart:rstop-orstart]
                fg_dep=ofg_dep[rstart-orstart:rstop-orstart]
                fg_off=ofg_off[rstart-orstart:rstop-orstart]
                ov=oov[rstart-orstart:rstop-orstart]

                fg_off_corr=fg_off.copy()
    
                try:
                    if var=='126':
                        try:
                            hbe = obe[rstart-orstart:rstop-orstart]
                            mask = np.isnan(hbe) & ~np.isnan(fg_off)
                            hbe[mask] = fah['RICH_bias_estimate'][rstart:rstop][mask]
                            fg_dep += f['era5fb']['biascorr@body'][rstart:rstop]
                            fg_off_corr[:]=fg_off-hbe #obe[rstart-orstart:rstop-orstart]
                            #mask = np.isnan(fg_off_corr) & ~np.isnan(fg_off)
                            #fg_off_corr[mask]=fg_off-fah['RICH_bias_estimate'][rstart:rstop]
                                
                        except:
                            pass 
                    elif var=='137':
                        try:
                            rqi=fri['39'][istart:istop+1]
                            rqstart,rqstop=rqi[[0,-1]]
                            rti=fri['126'][istart:istop+1]
                            rtstart,rtstop=rti[[0,-1]]
                        except:
                            print(fn, 'there should be q and T if there is dewpoint')
                            continue
                    
                        zt=foz[rtstart:rtstop]
                        zq=foz[rqstart:rqstop]
                        #ov = foo[rstart:rstop]
                        qov = foo[rqstart:rqstop]
                        qfgd = fe5[rqstart:rqstop]
                        try:
                            hbe=obe[rstart-orstart:rstop-orstart]
                            hbe[np.isnan(hbe)]=0.
                        except:
                            hbe=np.zeros(rstop-rstart)
                        try:
                            rbe=fah['RISE_bias_estimate'][rtstart:rtstop]
                        except:
                            rbe=np.zeros(rtstop-rtstart)
                        #if np.any(np.isnan(rbe)):
                            #fdicts[starts[ii]][var]={}
                            #print(fn, 'spurious adjustments')
                            #continue
                        fg_off_corr=fg_off.copy()
                        fg_dep[:]=np.nan
                        #print(istop-istart)
                        for i in range(istop-istart):
                            zti=zt[rti[i]-rtstart:rti[i+1]-rtstart]
                            zqi=zq[rqi[i]-rqstart:rqi[i+1]-rqstart]
                            if len(zqi) == 0:
                                continue
                            idx=np.searchsorted(zti,z[ri[i]-rstart:ri[i+1]-rstart])
                            idx[idx==zti.shape[0]] = zti.shape[0] - 1
                            qidx=np.searchsorted(zqi,z[ri[i]-rstart:ri[i+1]-rstart])
                            mask = qidx==zqi.shape[0]
                            qidx[mask] = zqi.shape[0] - 1
                            #tera=f['observations_table']['observation_value'][rti[i]:rti[i+1]][idx]-f['era5fb']['fg_depar@body'][rti[i]:rti[i+1]][idx]
                            #qera=f['observations_table']['observation_value'][rqi[i]:rqi[i+1]][qidx]-f['era5fb']['fg_depar@body'][rqi[i]:rqi[i+1]][qidx]
                            qera=qov[rqi[i]-rqstart:rqi[i+1]-rqstart][qidx]-qfgd[rqi[i]-rqstart:rqi[i+1]-rqstart][qidx]
                            qera[zqi[qidx]!=z[ri[i]-rstart:ri[i+1]-rstart]] = np.nan
                            dpera = np.full(ri[i+1]-ri[i], np.nan)
                            zi = z[ri[i]-rstart:ri[i+1]-rstart]
                            if np.any(~np.isnan(qera)):
                                
                                try:
                                    pidxq = np.searchsorted(zqi, plevs*100)
                            
                                    pidxq[pidxq==zqi.shape[0]] = zqi.shape[0] - 1
                                    mask = (zqi[pidxq]-plevs*100)==0
                                    qera[pidxq[~mask]] = np.nan
                                    vpdata = sh2vap(qera[pidxq][mask], zqi[pidxq][mask])
                                    for j in range(vpdata.shape[0]):
                                        vpdata[j] = bisection(func, 150, 350, vpdata[j])
                                    pidx = np.searchsorted(zi, zqi[pidxq[mask]])
                                    pidx[pidx==zi.shape[0]] = zi.shape[0] - 1
                                    dpera[pidx]=vpdata
                                except:
                                    print(i, fn, 'vpdata failed')
                            
                            try:
                                
                                if zti.shape[0] == zi.shape[0]:
                                    
                                    fg_off_corr[ri[i]-rstart:ri[i+1]-rstart] -= rbe[rti[i]-rtstart:rti[i+1]-rtstart][idx] + hbe[ri[i]-rstart:ri[i+1]-rstart]#[qidx]
                                elif zti.shape[0] < zi.shape[0]:
                                        
                                    fg_off_corr[ri[i]-rstart:ri[i+1]-rstart][:zti.shape[0]] -= rbe[rti[i]-rtstart:rti[i+1]-rtstart][idx[:zti.shape[0]]] + hbe[ri[i]-rstart:ri[i+1]-rstart][:zti.shape[0]]
                                else:
                                    fg_off_corr[ri[i]-rstart:ri[i+1]-rstart] -= rbe[rti[i]-rtstart:rti[i+1]-rtstart][idx] + hbe[ri[i]-rstart:ri[i+1]-rstart]#[qidx]
                            except Exception as e:
                                raise ValueError(fn)
                        
                            #fg_dep[ri[i]-rstart:ri[i+1]-rstart]=f['observations_table']['observation_value'][ri[i]:ri[i+1]]-dpera
                            fg_dep[ri[i]-rstart:ri[i+1]-rstart]=ov[ri[i]-rstart:ri[i+1]-rstart]-dpera
                        #print(rms(fg_dep),rms(f['observations_table']['observation_value'][ri[i]:ri[i+1]]),np.nanmean(fg_off))
                    elif var=='138':
                        try:
                            hbe = obe[rstart-orstart:rstop-orstart]
                            if np.any(~np.isnan(fg_off)):
                                fg_off_corr[:]=fg_off-hbe
                            else:
                                fg_off_corr[:]=fg_dep-hbe
                        except:
                            pass
                    elif var=='39':
                        try:
                            hbe = obe[rstart-orstart:rstop-orstart]
                            if np.any(~np.isnan(fg_off)):
                                fg_off_corr[:]=fg_off-hbe
                            else:
                                fg_off_corr[:]=fg_dep-hbe
                        except:
                            pass
                        x = 0
                    elif var=='139' or var == '140' or var == '34':
                        hbe = obe[rstart-orstart:rstop-orstart]
                        try:
                            fg_off_corr[:]=fg_off-hbe
                        except:
                            pass
                    elif var== '106':
                      
                        #ufg = foo[rstart:rstop]
                        try:
                            rui=fri['139'][istart:istop+1]
                            rustart,rustop=rui[[0,-1]]
                            rvi=fri['140'][istart:istop+1]
                            rvstart,rvstop=rvi[[0,-1]]
                            rwsi=fri['107'][istart:istop+1]
                            rwsstart,rwsstop=rwsi[[0,-1]]
                        except:
                            print(fn, 'there should be u and v if there is winddirection')
                            continue
                    
                        #zu=f['observations_table']['z_coordinate'][rustart:rustop] 
                        #zv=f['observations_table']['z_coordinate'][rvstart:rvstop]
                        uvo = foo[rustart:rustop]
                        vvo = foo[rvstart:rvstop]
                        wso = foo[rwsstart:rwsstop]
                        _, wd = wswd(uvo- fe5[rustart:rustop], vvo- fe5[rvstart:rvstop])  # archived ERA5 wd background departure needs to be calculated from u,v
                        
                        speed_thresh = 2.0
                        try:
                            if ov.shape[0] == rustop - rustart and wso.shape[0] == rustop - rustart:
                                ov[wso<speed_thresh] = np.nan
                                fg_dep = ov - wd
                                fg_dep[fg_dep>180.] -= 360
                                fg_dep[fg_dep<-180.] += 360
                            else:
                                
                                zu=foz[rustart:rustop]
                                zws=foz[rwsstart:rwsstop]
                                j = 0
                                for i in range(istop-istart):
                                    zui=zu[rui[i]-rustart:rui[i+1]-rustart]
                                    if len(zui) == 0:
                                        continue
                                    zi = z[ri[i]-rstart:ri[i+1]-rstart]
                                    if rui[i+1] - rui[i] != ri[i+1] - ri[i] or rwsi[i+1] - rwsi[i] != ri[i+1] - ri[i]:
                                        x = 0
                                        zwsi=zws[rwsi[i]-rwsstart:rwsi[i+1]-rwsstart]
                                        idx=np.searchsorted(zi, zui)
                                        idx[idx==zi.shape[0]] = zi.shape[0] - 1
                                        mask = zi[idx] == zui
                                        idy=np.searchsorted(zwsi, zi)
                                        idy[idy==zwsi.shape[0]] = zwsi.shape[0] - 1
                                        ymask = zwsi[idy] == zi
                                        
                                        fgd = np.full(ri[i+1]-ri[i], np.nan)
                                        ufgm = ov[ri[i]-rstart:ri[i+1]-rstart][idx[mask]]
                                        #wsom = ~np.isnan(wso[rwsi[i]-rwsstart:rwsi[i+1]-rwsstart][idy[ymask]])
                                        try:
                                            
                                            ufgm[wso[rwsi[i]-rwsstart:rwsi[i+1]-rwsstart][idy[ymask]]<speed_thresh] = np.nan
                                            fgd[idx[mask]] = ufgm - wd[rui[i]-rustart:rui[i+1]-rustart][mask]
                                        except:
                                            print(fn,i, ri[i], ufgm.shape, mask.shape, 'ws,wd,u,v do not match')
                                    else:
                                        fgd = ov[ri[i]-rstart:ri[i+1]-rstart]
                                        fgd[wso[rwsi[i]-rwsstart:rwsi[i+1]-rwsstart] < speed_thresh] = np.nan
                                        fgd -= wd[rui[i]-rustart:rui[i+1]-rustart]
                                    fgd[fgd>180.] -= 360.
                                    fgd[fgd<-180.] += 360
                                    fg_dep[ri[i]-rstart:ri[i+1]-rstart] = fgd
                                    
                            fg_dep[~np.isin(z, plevs*100)] = np.nan
                                
                                
                        except MemoryError as e:
                            raise ValueError(fn) #,rustart, 'wd departure calculation failed')
                            pass
                        
                        
                        try:
                            hbe=obe[rstart-orstart:rstop-orstart]
                            hbe[np.isnan(hbe)]=0.
                        except:
                            hbe=np.zeros(rstop-rstart)

                        fg_off_corr = fg_off - hbe
                        fg_off_corr[fg_off_corr>180.] -= 360.
                        fg_off_corr[fg_off_corr<-180.] += 360
                        
                        #fg_off[np.isnan(fg_dep)] = np.nan
                        #fg_off_corr[np.isnan(fg_dep)] = np.nan

                        if False: #and rms(fg_off)[0] / rms(fg_dep)[0] > 1.1:
                            
                            plt.subplot(2, 1, 1)
                            plt.plot(fg_dep.flatten(), '.', label=f'fg_dep {rms(fg_dep)[0]:.2f} {rms(fg_dep)[1]}')
                            plt.plot(fg_off.flatten(), 'o', label=f'fg_off  {rms(fg_off)[0]:.2f} {rms(fg_off)[1]}')
                            plt.plot(fg_off_corr.flatten(), '.', label=f'fg_off_corr  {rms(fg_off_corr)[0]:.2f} {rms(fg_off_corr)[1]}')
                            plt.legend()
                            plt.title(f'{fn.split('/')[-1]} {dt[0] / 86400 / 365.25 + 1900}')
                            #plt.subplot(3, 1, 2)
                            #plt.plot(fg_off.flatten()-fg_off_corr.flatten(), 'o', label=f'fg_off  {rms(fg_off-fg_off_corr)[0]:.2f}')
                            #plt.legend()
                            #plt.subplot(2, 1, 2)
                            #plt.plot(wso.flatten(), '.', label=f'wso {rms(wso)[0]:.2f} {rms(wso)[1]}')
                            #plt.legend()
                            #plt.tight_layout()
                            plt.show()
                            x = 0
                        
                    elif var== '107':
                        #fg = f['observations_table']['observation_value'][rstart:rstop]
                        
                        #try:
                            #fg_off_corr[:]=fg_off-obe[rstart-orstart:rstop-orstart]
                        #except:
                        try:
                            hbe=obe[rstart-orstart:rstop-orstart]
                            hbe[np.isnan(hbe)]=0.
                        except:
                            hbe=np.zeros(rstop-rstart)

                        pass
                    
                    dtn = np.unique(dt)
                    if var == '106':  # ensure that there as many wind speeds as wind directions - needed for later check of windspeed to eliminate direction departures at low wind speeds
                        try:
                            
                            if len(dtn) != len(np.unique(fo['date_time'][rstarts['107'][istart]:rstops['107'][istart]])):
                                continue
                        except:
                            continue
                    if var == '107':  # ensure that there as many wind speeds as wind directions - needed for later check of windspeed to eliminate direction departures at low wind speeds
                        try:
                            
                            if len(dtn) != len(np.unique(fo['date_time'][rstarts['106'][istart]:rstops['106'][istart]])):
                                continue
                        except:
                            continue
                            
                    #if starts[ii] == 2637532800:
                        
                        #print(starts[ii], var, len(dtn), orstop-orstart)
#                    ld = list(zip(('fg_dep','fg_off','fg_off_corr', 'adj', 'obs', 'obs_adj'), (np.float32, np.float32, np.float32, np.float32, np.float32, np.float32)))
                    #ld = (np.float32, np.float32, np.float32, np.float32, np.float32, np.float32))
                    deps = np.full((len(params), len(plevs), dtn.shape[0]), np.nan, dtype=np.float32)
                    #deps['dt'][0, :] = dtn
                    #deps[['fg_dep','fg_off','fg_off_corr', 'adj', 'obs', 'obs_adj']][:] =np.nan

                    fdicts[starts[ii]][var]={'dt':dtn, 'params' : params}
                    for ip in range(len(plevs)):
                        idx=np.where((z==plevs[ip]*100)&(~np.isnan(fg_off))&repdict[istart]) [0]  #&(~np.isnan(fg_dep))
                        nidx = np.searchsorted(dtn, dt[idx])
                        #print(plev,len(idx))
                        mix = reportype is not None and 16999 in reportype
                        if mix:
                            idy=np.where((z==plevs[ip]*100)&(~np.isnan(fg_off))&(~repdict[istart])) [0]
                            nidy = np.searchsorted(dtn, dt[idy])
                            if len(idx)>0 and len(idy) > 0:
                                deps[0, ip,nidx] = np.nanmean(fg_dep[idx]) - np.nanmean(fg_dep[idy])
                                deps[1,ip,nidx] = np.nanmean(fg_off[idx]) - np.nanmean(fg_off[idy])
                                deps[2,ip,nidx] = np.nanmean(fg_off_corr[idx]) - np.nanmean(fg_off_corr[idy])
                                deps[3,ip, nidx] = np.nanmean(ov[idx]) - np.nanmean(ov[idy]) # obs
                                deps[4,ip, nidx] = hbe[idx] #fg_off[idx] - fg_off_corr[idx] #adj
                                deps[5,ip, nidx] = ov[idx] - hbe[idx] #deps[4][ip, nidx] - deps[3][ip, nidx] #obs-adj
                        else:
                            if len(idx)>0:
                                deps[0, ip,nidx] = fg_dep[idx]
                                deps[1,ip,nidx] = fg_off[idx]
                                deps[2,ip,nidx] = fg_off_corr[idx]
                                deps[3,ip, nidx] = ov[idx] # obs
                                deps[4,ip, nidx] = hbe[idx] #fg_off[idx] - fg_off_corr[idx] #adj
                                deps[5,ip, nidx] = ov[idx] - hbe[idx] #deps[4][ip, nidx] - deps[3][ip, nidx] #obs-adj
                    
                                 
                    fdicts[starts[ii]][var]['deps'] = deps * factors[var]
                except MemoryError as e:
                    raise ValueError(fn+' '+str(e))
    for istart in starts: #'137','39', '138':#,'126','139', '140':
        try:
            
            if fdicts[istart]['106']['deps'].shape[2] != fdicts[istart]['107']['deps'].shape[2]:
                
                print(f'{os.path.basename(fn)}, {istart}, {fdicts[istart]['106']['deps'].shape[2]},{fdicts[istart]['107']['deps'].shape[2]},{time.time()-tt:5.4f}')
        except Exception as e:
            #print(e)
            pass

    print(f'{os.path.basename(fn)}, {len(fdicts) - 3}, {mingood}, {time.time()-tt:5.4f}')
    if len(fdicts) -3 < mingood:
        return None
    fdicts_ref = ray.put(fdicts)
    return fdicts_ref

ray_read_adj = ray.remote(read_adj)

@njit(boundscheck=True)
def nf(f, ddv, val):
    
    #fd = f[d]
    for ip in range(f.shape[0]):
        for it in range(f.shape[1]):
            if f[ip, it] == f[ip, it]:
                ddv[ip, val[ip]] = f[ip, it]
                val[ip] += 1
    return #ddv,val

@njit(boundscheck=True)
def nf2(f,ddv, val):
    
    #fd = f[d]
    if f.shape[1] == 0:
        return
    for ip in range(f.shape[0]):
        for it in range(f.shape[1]):
            if f[ip, it] == f[ip, it]:
                ddv[ip, val[ip]] = f[ip, it]
                val[ip] += 1
    return #ddv,val

@njit(boundscheck=True)
def nf3(f,k, ddv, val):
    
    #fd = f[d]
    if f.shape[2] == 0:
        return
    for ip in range(f.shape[1]):
        for it in range(f.shape[2]):
            if f[0, ip, it] == f[0, ip, it] or f[1, ip, it] == f[1, ip, it] or f[2, ip, it] == f[2, ip, it]:
                ddv[ip, val[ip]] = f[k, ip, it]
                val[ip] += 1
    return #ddv,val

@njit#(parallel=True)
def anomtrend(dt, ts, trend, mingood):
    for ip in range(ts.shape[0]):
        for par in prange(ts.shape[1]):
            #janomaly(ts[ip, par,:])
            g = 0
            for i in range(ts.shape[2]):
                if ts[ip, par, i] == ts[ip, par, i]:
                    g += 1
            #print(ip, par, g)
            if g >= mingood:                
                trend[ip, par] = fastlinregress(dt, ts[ip, par, :])
            else:
                trend[ip, par] = np.nan
    return

def do_trend(ff, var, params, tint, plevs, glist, mingood=20, vgood=30):
    ts = np.zeros((len(glist), len(params), len(plevs)))
    dt = []
    t = 0
    
    trend = np.full((plevs.shape[0], params.shape[0]),  np.nan)
    #fk = list(ff.keys())[3:]
    #if var not in ff[fk[0]].keys():
        
        #return trend
    for istart in glist:
        #if type(istart) == np.int64:
            
            #if istart >= tint[0] and istart < tint[1]:
               # if var in ff[istart].keys(): #and par in ff[istart][var].dtype.fields.keys():
                    
                    #ts.append(np.nanmean(ff[istart][var][par], axis=1))
                    #if np.any(ff[istart][var]['params'] != params):
                        #raise ValueError(ff['fn']+' params mismatch')
                    if var in ff[istart].keys():
                        
                        a = ff[istart][var]['deps']
                        if a.shape[-1] > vgood:
                            nnanmean(a,ts,t, vgood, quant=None)
                            t += 1
                            dt.append(istart)
            #if istart >= tint[1]:
                #break
    
    #return trend
    if len(dt) < mingood:
        return trend
    
    dt = np.array(dt)
    
    
#def janomaly( series,startyear,interval,anomaly,climatology):
#def fastlinregress(x,y):
    ts = ts[:t, :, :].T
    anomtrend(dt/86400/365/10, ts, trend, mingood) #np.concatenate(ts).T
    
    
    #for par in range(len(params)):
        #for ip in range(len(plevs)):
            #mask = ~np.isnan(ts[ip, par, :])
            #if np.sum(mask) != 0:
        
                #trend[ip, par] = np.polyfit(dt[mask]/86400/365/10, ts[ip,par, mask], 1)[0]
    
    return trend

def do_func(ff, var, params, tint, plevs, glist, func, mingood=20, vgood=30, quant=None):
    ts = np.zeros((len(glist), len(params), len(plevs)))
    dt = []
    t = 0
    
    trend = np.full((plevs.shape[0], params.shape[0]),  np.nan)

    for istart in glist:
        if var in  ff[istart].keys():
            
            a = ff[istart][var]['deps']
            if a.shape[-1] > vgood:
                func(a,ts,t, vgood, quant)
                t += 1
                dt.append(istart)
                

    if len(dt) < mingood:
        return trend
    
    trend = np.nanmean(ts[:t, :3, :], axis=0).T
    #anomtrend(dt/86400/365/10, ts, trend, mingood) #np.concatenate(ts).T
        
    return trend

def xdo_rmsdev(ff, var, params, tint, plevs, glist, mingood=20, vgood=30):
    ts = np.zeros((len(glist), len(params), len(plevs)))
    dt = []
    t = 0
    
    trend = np.full((plevs.shape[0], params.shape[0]),  np.nan)

    for istart in glist:
        if var in  ff[istart].keys():
            
            a = ff[istart][var]['deps']
            if a.shape[-1] > vgood:
                nnanrms(a,ts,t, vgood)
                t += 1
                dt.append(istart)
                

    if len(dt) < mingood:
        return trend
    
    trend = np.nanmean(ts[:t, :3, :], axis=0).T
    #anomtrend(dt/86400/365/10, ts, trend, mingood) #np.concatenate(ts).T
        
    return trend

def xdo_meandev(ff, var, params, tint, plevs, glist, mingood=20, vgood=30):
    ts = np.zeros((len(glist), len(params), len(plevs)))
    dt = []
    t = 0
    
    trend = np.full((plevs.shape[0], params.shape[0]),  np.nan)

    for istart in glist:
        if var in  ff[istart].keys():
            
            a = ff[istart][var]['deps']
            if a.shape[-1] > vgood:
                nnanmean(a,ts,t, vgood)
                t += 1
                dt.append(istart)
                

    if len(dt) < mingood:
        return trend
    
    trend = np.nanmean(ts[:t, :3, :], axis=0).T
    #anomtrend(dt/86400/365/10, ts, trend, mingood) #np.concatenate(ts).T
        
    return trend


def xdo_quant(ff, var, params, tint, plevs, glist, quant, mingood=20, vgood=30):
    ts = np.zeros((len(glist), len(params), len(plevs)))
    dt = []
    t = 0
    
    trend = np.full((plevs.shape[0], params.shape[0]),  np.nan)

    for istart in glist:
        a = ff[istart][var]['deps']
        if a.shape[-1] > vgood:
            nnanquant(a,ts,t,quant, vgood)
            t += 1
            dt.append(istart)

    if len(dt) < mingood:
        return trend
    
    trend = np.nanmean(ts[:t, :3, :], axis=0).T
    #anomtrend(dt/86400/365/10, ts, trend, mingood) #np.concatenate(ts).T
        
    return trend

   
def trendsperstation(ff, starts, vtup, params, plevs, max_miss, vgood, intervals):
    
    trends = {}
    coords = {}
    for interval in intervals :
        if interval not in trends.keys():
            trends[interval] = {}
            coords[interval] = {}
            ystart = np.int64((datetime(interval[0], 1, 1) - datetime(1900, 1, 1)).total_seconds())
            ystop = np.int64((datetime(interval [1]+1, 1, 1) - datetime(1900, 1, 1)).total_seconds())
            slist = [start for start in starts if start >= ystart and start < ystop]
            mingood = np.int32((interval[1]+1-interval[0])*12 -max_miss)
        glist = [k for k in ff.keys() if k in slist]
        if len(glist) == 0 or len(glist) < mingood:
            #print(interval, len(glist), mingood)
            continue
        #tt = time.time()
        for v in vtup:
            g = 0
            for k in glist:
                if v in ff[k].keys():
                    g += 1
            if g < mingood:
                continue
            if v not in trends[interval].keys():
                trends[interval][v] = {}
                coords[interval][v] = []
            coords[interval][v].append({'lat': ff['lat'], 'lon': ff['lon']})
            for param in params:
                if param not in trends[interval][v].keys():
                    trends[interval][v][param] = []
            x = do_trend(ff, v, params, (ystart, ystop), plevs, glist, mingood=mingood, vgood=vgood)
            for i in range(len(params)):                       
                trends[interval][v][params[i]].append(x[:, i])
        #print(time.time()-tt)
        x = 0
    return trends, coords

ray_trendsperstation = ray.remote(trendsperstation)

def funcdevperstation(ff, starts, vtup, params, plevs, vgood, func, intervals, quant=None):
    
    trends = {}
    coords = {}
    for interval in intervals :
        if interval not in trends.keys():
            trends[interval] = {}
            coords[interval] = {}
            ystart = np.int64((datetime(interval[0], 1, 1) - datetime(1900, 1, 1)).total_seconds())
            ystop = np.int64((datetime(interval [1]+1, 1, 1) - datetime(1900, 1, 1)).total_seconds())
            slist = [start for start in starts if start >= ystart and start < ystop]
            mingood = np.int32(len(slist) *0.9)
        glist = [k for k in ff.keys() if k in slist]
        if len(glist) < mingood:
            continue
        #tt = time.time()
        for v in vtup:
            g = 0
            for k in glist:
                if v in ff[k].keys():
                    g += 1
            if g < mingood:
                continue
            if v not in trends[interval].keys():
                trends[interval][v] = {}
                coords[interval][v] = []
            coords[interval][v].append({'lat': ff['lat'], 'lon': ff['lon']})
            for param in params:
                if param not in trends[interval][v].keys():
                    trends[interval][v][param] = []
            x = do_func(ff, v, params, (ystart, ystop), plevs, glist, func, mingood=mingood, vgood=vgood, quant=quant)
            for i in range(3):                       
                trends[interval][v][params[i]].append(x[:, i])
        #print(time.time()-tt)
        x = 0
    return trends, coords

ray_funcdevperstation = ray.remote(funcdevperstation)

def plot_hovmoellers(fdicts_ref, ref, starts, plevs, vname, vtup, rversion):
    tt = time.time()
    fdicts = ray.get(fdicts_ref)
    ratios={'fg_off':dict(zip(vtup,[[] for v in vtup])),'fg_off_corr':dict(zip(vtup,[[] for v in vtup])),'off':dict(zip(vtup,[[] for v in vtup]))}
    for start in starts:
        deps={'fg_dep':{},'fg_off':{},'fg_off_corr':{}}
        vals={'fg_dep':{},'fg_off':{},'fg_off_corr':{}}
        anan = np.array([np.nan])
        m=0
        for d in deps:
            for var in vtup: #
                deps[d][var]=np.full((plevs.shape[0], 200000), np.nan)
                vals[d][var] = np.zeros(plevs.shape[0], dtype=np.int32)
        for ff in fdicts:
            
            if start in ff.keys():
                for var in vtup: #
                    if var in ff[start].keys():
                        f = ff[start][var]['deps']
                        if f.shape[-1] == 0:
                            continue
                        k = 0
                        for d in deps:                                    
                            #fff = f[k].copy()
                            #if var == '106':
                                #fff[ff[start]['107']['deps'][k]<5.0] = np.nan
                            #if k == 2 and var == '106':
                                #fff[fff > 180.] -= 360.
                                #fff[fff < -180.] += 360.
                            #nf2(f[k],deps[d][var],vals[d][var] )
                            nf3(f,k, deps[d][var],vals[d][var] )
                            k += 1
                            
        #print('start0:', start, time.time()-tt)
        for var in vtup: #
            l = 0
            r = {}
            #v = vals['fg_dep'][var]
            for d in deps:
                l += 1
                ddv = deps[d][var]
                v = vals[d][var]
                r[d] = np.full(len(plevs), np.nan)
                r[d+'2'] = np.full(len(plevs), np.nan)
                for ip in range(len(plevs)):
                    if v[ip] > 30:
                        samp = ddv[ip, :v[ip]]
                        x1 = rms(samp)
                        if x1[1] > 30:                                         
                            r[d][ip] = x1[0]
                        if l >= 2:
                            samp[np.isnan(deps['fg_dep'][var][ip, :v[ip]])] = np.nan
                            x2 = rms(samp)
                            v[ip] = x2[1]
                            if x2[1] > 30:                                         
                                r[d+'2'][ip] = x2[0]
                
                rref=r['fg_dep']
                #vref = vals['fg_dep'][var]
                #print(r)
                if l>1:
                    if var in ['137','138','39']:
                        r[d][plevs<100.]=np.nan
                        r[d+'2'][plevs<100.]=np.nan
                    
                    ratios[d][var].append(r[d+'2']/rref*100)
                    #for ip in range(len(plevs)):
                        #if vref[ip] > 30:
                            #if vref[ip] / v[ip] < 0.2 or vref[ip] / v[ip] > 5:
                                #ratios[d][var][-1][ip] = np.nan
                        #else:    
                            #ratios[d][var][-1][ip] = np.nan
                    
                    if np.any((ratios[d][var][-1]>130)|(ratios[d][var][-1]<70)):
                        idx = np.where((ratios[d][var][-1]>130)|(ratios[d][var][-1]<70))
                        print(ref+timedelta(seconds=int(start)),d,var,plevs[idx], ratios[d][var][-1][idx])
                        if np.any(ratios[d][var][-1]>500):
                            print(ref+timedelta(seconds=int(start)),d,var,plevs[idx], ratios[d][var][-1][idx])
                            x=0
                                  
                        #print(y,quart,d,var,ratios['fg_off'][var][-1])
                    rref2=r['fg_off'] #deps['fg_off'][var][:, 0]
                    if l > 2:
                        ratios['off'][var].append(r[d]/rref2*100)
                    

        print('start:', start, time.time()-tt)

    m=0
    fig = plt.figure(figsize=(10, len(vtup)*1.5))
    years = ((ref + timedelta(seconds=int(starts[0]))).year, (ref + timedelta(seconds=int(starts[-1]))).year)
    qm = (len(starts) - 1) // (years[1] - years[0])
    for var in vtup: #
        l=0
        for d in ratios.keys():
            m+=1
            l+=1
            ax=fig.add_subplot(len(vtup),3,m)
            ratios[d][var]=np.array(ratios[d][var])
            mgx,mgp=np.meshgrid(years[0]+np.arange((years[1]-years[0])*qm)/qm,plevs)
            #print(var,d,ratios[d][var].shape,mgx.shape)            
            cm = ax.contourf(mgx,mgp,ratios[d][var][:-1, :].T,np.linspace(75,115,9),cmap='gist_ncar')
            #plt.contour(mgx,mgp,ratios[d][var].T,[80,100],'k')
            from matplotlib.ticker import FormatStrFormatter

            ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
            ax.set_ylabel('Pressure / hPa')
            ax.set_ylim(1000,50)
            if d == 'off':
                ax.set_title(f'fg_{d}_corr/fg_off, {vname[var]}, {rms(ratios[d][var][:-1, :].flatten())[0]:5.2f}')
            else:
                ax.set_title(f'{d}/ERA5 fg, {vname[var]}, {rms(ratios[d][var][:-1, :].flatten())[0]:5.2f}')
            ax2=plt.colorbar(cm)
            ax2.set_label('%', rotation=270)

    #plt.subplot(len(vtup),2,2*len(vtup))
    fig.tight_layout()
    fig.savefig(f'hov_{rversion}_{'-'.join(vtup)}_{years[0]}-{years[1]}.png')
    plt.close()

def getsize(fn, starts, sdict, mingood, vgood):
    
    with h5py.File(fn) as f:
        fri = f['recordindices']
        ril = list(fri.keys())
        rt = fri['recordtimestamp']
        if len(rt) > mingood * vgood:            
            ts = rt[:]
            idx = np.searchsorted(ts, starts)
            idx = idx[(idx<ts.shape[0])&(idx>0)]
            if len(idx) > mingood:
                s = np.sum(idx[1:]-idx[:-1]>=vgood)
                if s >= mingood:
                    #rk = fri[ril[1]][[0, -1]]
                    #sondetype = np.unique(f['observations_table']['sensor_id'][rk[0]:rk[1]].view('S4'))
                    #st = []
                    #for sdk in sdict.keys():
                        #if 'rs' in sdk:
                            #st = st + sdict[sdk]
                    #for rst in sondetype:
                        #if rst.decode() in st:
                            #return fn
                        
                    
                    return fn
                else:
                    return None
        else:
            return None
        
ray_getsize = ray.remote(getsize)

def quantileplot(allqs, goods, da , i): #wpath, interval, vtup, params, plevs, vname, units
    
    pltparams = {'legend.fontsize': 'small',
              'figure.figsize': (12, 8),
             'axes.labelsize': 'medium',
             'axes.titlesize': 'large',
             'text.usetex' : False,
             'backend': 'agg',
             'xtick.labelsize':'medium',
             'ytick.labelsize':'medium'}
    plt.rcParams.update(pltparams)
    for v in da['vtup']:
        j = 0
        plt.figure(figsize=(6, 3*2))
        for param in da['params'][:3]:
            ax = plt.subplot(3,1,j+1 )
            profs = np.vstack(allqs[i:i+len(da['plevs'])])
            if da['title'] == 'mean':
                plt.semilogy(profs[da['pplevs'], -3:], da['plevs'][da['pplevs']], label=[f'{int(da['qsv'][-3 + i]*100)}' for i in range(3)])
            else:
                plt.semilogy(profs[da['pplevs'], -3:], da['plevs'][da['pplevs']], label=[f'{da['qsv'][-3 + i]*100:.1f}' for i in range(3)])
            handles, labels = ax.get_legend_handles_labels() 
            
            plt.ylim(1000, da['plevs'][da['pplevs']][0])
            plt.ylabel('p / hPa')
            if j == 2:                        
                plt.xlabel(f'{da['vnames'][v]} / {da['units'][v]}')
            ax2 = ax.twiny()
            ax2.plot(goods[i:i+len(da['plevs'])][da['pplevs']], da['plevs'][da['pplevs']], color=(0.7, 0.7, 0.7), label='count')
            ax2.set_xlim(0, np.max(goods[i:i+len(da['plevs'])][da['pplevs']])*1.3)
            handles2, labels2 = ax2.get_legend_handles_labels() 
            plt.legend(handles+handles2, labels+labels2, loc='upper right')
            #ax2.xlabel('# of stations')
            reptype = ''
            if da['reportype'] is not None :
                reptype = da['reportype']
            plt.title(f'{da['vnames'][v]}-{da['title']}-{param} {reptype} quantiles, {da['interval'][0]}-{da['interval'][1]}')
            i += len(da['plevs'])
            j += 1
        plt.tight_layout()
        plt.savefig(f'{da['wpath']}/{da['vnames'][v]}-{da['title']}_{reptype}quantiles_{da['interval'][0]}-{da['interval'][1]}')
        plt.close()
    x = 0
    
ray_quantileplot = ray.remote(quantileplot)

def sub(sdict):
    
    ref=datetime(1900,1,1)
    plevs=np.array((10, 20, 30, 50,70, 100,150,200,250, 300,400, 500,700,850,925, 1000))
    ranges={plev: {'137':(-5,5),'138':(-10, 10),'126':(-2, 2),'139':(-2, 2),'140':(-2, 2),'39':(-1, 1), '107': (-2, 2),'106': (-10, 10),'34': (-5, 5),} for plev in plevs}
    for plev in plevs[plevs<500]:
        ranges[plev]['138'] = [-4,4 ]
        ranges[plev]['39'] = [-0.1,0.1]
    for plev in plevs[plevs<150]:
        ranges[plev]['138'] = [-2,2 ]
        ranges[plev]['39'] = [-0.01,0.01]
    rmsranges = {}
    for plev in plevs:
        rmsranges[plev] = {}
        for k in ranges[plev].keys():
            rmsranges[plev][k] = [0, -ranges[plev][k][0] + ranges[plev][k][1]]
            
        
    factors={'137':1,'138':100,'126':1,'139':1,'140':1,'39':1000, '107': 1,'106': 1,'34': 1,'0': 1,}
    units={'137':'K','138':'%','126':'K','139':'m/s','140':'m/s','39':'g/kg', '107': 'm/s','106': 'deg','34': 'K',}
    vname={'137':'Td','138':'RH','126':'T','139':'u-wind','140':'v-wind','39':'q', '107': 'windspeed','106': 'winddir','34': 'T-Td'}
    bedict={'137':'humidity_bias_estimate','138':'humidity_bias_estimate','126':'RISE_bias_estimate','139':'wind_bias_estimate',
            '140':'wind_bias_estimate','39':'humidity_bias_estimate', '107': 'wind_bias_estimate','106': 'wind_bias_estimate', '34': 'humidity_bias_estimate',}
    vtup= '107','106' , '140', '139','126'#,'34','39','138', '107'#,'137'#,'137','139','140'
    #vtup= '126','34','39','138' #,'137'#,'137','139','140'
        
    params = np.array(('fg_dep','fg_off', 'fg_off_corr', 'obs', 'adj', 'obs_adj'))
    set_num_threads(len(params))
    vgood = 5
     
    versions = '29',# '25'
    intervals = [(1905 + t, 1925 + t) for t in range(65, 111, 5)]
    intervals = [(1979, 2019)]
    rmsintervals = intervals
    #rmsintervals = [(y, y + 1) for y in range(1940, 1981, 5)] #[(1957, 1959), (1986, 1987), (2006, 2007), (2016, 2017), (2023, 2024)]
    years=(np.min((np.min(np.array(intervals)), np.min(np.array(rmsintervals)))), np.max((np.max(np.array(intervals)), np.max(np.array(rmsintervals))))+ 1)
    #years=(np.min(np.array(rmsintervals)), np.max(np.array(rmsintervals)) + 1)
    #years = (1920, 1960)
    pilot = np.array((16013, 16068,))
    temp = np.array((16022, 16045))
    
    ncpu = 80
    mp = False
    if(mp):
        P=Pool(80)
    else:
        ray.init(address='131.130.157.5:6379', ignore_reinit_error=True) #(num_cpus=80, ignore_reinit_error=True)
    
    tt = time.time()
    for rversion in versions: #
        resorted=f'/mnt/users/scratch/leo/scratch/converted_v{rversion}/long/'
        #fns=glob.glob(f"{resorted}/{WIGOS}*.nc")
        #fns=glob.glob(f"{resorted}/0-2000[01]-[0-7]*v3.nc")
        tt = time.time()
        #fns=np.array(glob.glob(f"{resorted}/*-0-94[4-9]*v3.nc"))
        fns=np.array(glob.glob(f"{resorted}/*v3.nc"))
        sizes = [os.path.getsize(fn) for fn in fns]
        fnx = fns[np.argsort(sizes)[::-1]]
        
        #fns=[fn for fn in fnsa if '20999'  in fn]
    
        efiles={}
        starts = []
        qm = 12  # 1 for annual, 4 for quarterly, 12 for monthly
        for y in range(*years):
            for quart in range(1, 13, 12//qm):
    
                if qm == 4:
                    qq = quart - 1
                    if qq == 0:
                        
                        yy=y-1
                        qq=12
                else:
                    yy=y
                    qq=quart
                dstart=datetime(yy,qq,1)
                qq = qq + 12 // qm
                if qq > 12:
                    qq -= 12
                    yy += 1
                dstop=datetime(yy,qq,1)
                #print(dstart,dstop)
                start=(dstart-ref).total_seconds()
                stop=(dstop-ref).total_seconds()
                starts.append(np.int64(start))
        starts.append(np.int64(stop))
        starts = np.array(starts)
        starts_ref = ray.put(starts)
    
        #print(resorted,WIGOS,fns)
        #print(rversion,start,stop)
        mingood = 24 # 19 years (assuming min interval is 21 years)
        
        fns = []
        futures = []
        print('vor len(fns)',len(fns), time.time()-tt)
        if len(fnx) < 5:
            for fn in fnx:
                fns.append(getsize(fn, starts, sdict, mingood, vgood))
            
        else:
            futures = [ray_getsize.remote(fn, starts_ref, sdict_ref, mingood, vgood) for fn in fnx]
            fns = ray.get(futures)
        fns = [fn for fn in fns if fn is not None]#[:49]# and '0-20000-0-64387' in fn]

        print('len(fns)',len(fns), time.time()-tt)
        
        if(mp):
            pfunc=partial(read_adj,starts,plevs, factors, bedict, vtup, params, mingood)
            idx = 0
            fdicts=list(P.map(pfunc,fns[idx:]))
        else:
            futures = []
            fdicts = []
            if len(fns) < 50:
                for fn in fns:
                    fdicts.append(read_adj(starts, plevs, factors, bedict, vtup, params,mingood,None, None, fn))
            else:
                for fn in fns:
                    futures.append(ray_read_adj.remote(starts_ref, plevs,factors, bedict, vtup, params,mingood,None,None, fn))
                fdicts = ray.get(futures)

        fdicts_ref = [ff for ff in fdicts if ff is not None]
        reportype = ray.get(fdicts_ref[0])['reportype']
        reptype = ''
        if reportype is not None :
            if 16022 in reportype or 16045 in reportype:
                reptype = 'TEMP'
            else:
                if 16999 in reportype:
                    reptype = 'DIFF'
                else:
                    reptype = 'PILOT'

        sdict = ray.get(fdicts_ref[0])['sdict']
        print()
        print(len(fdicts), time.time()-tt)
        tt = time.time()
        npz = np.zeros((plevs.shape[0], 200000))
        hash = dict(zip(plevs, np.arange(plevs.shape[0])))
        if False:
            plot_hovmoellers(fdicts_ref, ref, starts, plevs, vname, vtup, rversion)
            
        futures = []
        trends = []
        for ff in fdicts_ref:
            futures.append(ray_trendsperstation.remote(ff, starts_ref,vtup,params, plevs, mingood, vgood, intervals))
            #trends.append(trendsperstation(ray.get(ff), starts,vtup,params, plevs, mingood, vgood, intervals))
        
        trends = ray.get(futures)
        print('Trend calc', time.time()-tt)
        
        futures = []
        rmsdevs = []
        alldevs = {'mean': {'func': nnanmean, 'quant': None,},'rms': {'func': nnanrms, 'quant': None,}, '80': {'func': nnanquant, 'quant': 0.8,}}
        for d in alldevs.keys():
            alldevs[d]['devs']=[]
            for ff in fdicts_ref:
                futures.append(ray_funcdevperstation.remote(ff, starts_ref,vtup,params, plevs, vgood,alldevs[d]['func'], rmsintervals, quant=alldevs[d]['quant']))
                #alldevs[d]['devs'].append(funcdevperstation(ray.get(ff), starts,vtup,params, plevs, vgood,alldevs[d]['func'], rmsintervals, quant=alldevs[d]['quant']))
        
            alldevs[d]['devs'] = ray.get(futures)
        
        print('rmsdev calc', time.time()-tt)
        
        futures = []
        for interval in intervals :
            for v in vtup:
                for param in params:
                    tr = np.array([t[0][interval][v][param] for t in trends if v in t[0][interval].keys()]).T
                    vlats = np.array([t[1][interval][v][0]['lat'] for t in trends if v in t[0][interval].keys()])
                    vlons = np.array([t[1][interval][v][0]['lon'] for t in trends if v in t[0][interval].keys()])
                    if tr.ndim < 3:
                        continue
                    for ip in range(len(plevs)):                        
                        mask = ~np.isnan(tr[ip,0, :])
                        if np.sum(mask) > 0:
                            
                            da = {'wpath': os.path.expanduser('~/CEUAS/CEUAS/public/resort'),'version' : rversion, 'param': param,'v' : v, 'vname': vname, 'units': units[v], 'trend': tr[ip,0, mask], 'interval' : interval, 'plev': plevs[ip],
                                  'lat': vlats[mask], 'lon': vlons[mask], 'fn': [], 'ranges' : ranges}
                            #plt_trends(da)
                            futures.append(ray_plt_trends.remote(da))
                        
        ray.get(futures)
                        
        print('Trend plot', time.time()-tt)
        futures = []
        #allqs = []
        #goods = []
        pplevs = np.searchsorted(plevs, (20, 50, 100, 150, 200, 300, 500, 700, 850, 1000))
        for d in alldevs.keys():
            
            qsv = np.arange(0.85, 1.01, 0.025)
            qsv[-1] -= 0.005
            alldevs[d]['allqs'] = []
            alldevs[d]['goods'] = []
            for interval in rmsintervals :
                for v in vtup:
                    for param in params[:3]:
                        tr = np.array([t[0][interval][v][param] for t in alldevs[d]['devs'] if v in t[0][interval].keys()]).T
                        if tr.ndim < 3:
                            continue
                        
                        alldevs[d]['goods'].append(np.sum(~np.isnan(tr), axis=2).flatten())
                        vlats = np.array([t[1][interval][v][0]['lat'] for t in alldevs[d]['devs'] if v in t[0][interval].keys()])
                        vlons = np.array([t[1][interval][v][0]['lon'] for t in alldevs[d]['devs'] if v in t[0][interval].keys()])
                        for ip in range(len(plevs)):                        
                            mask = ~np.isnan(tr[ip,0, :])
                            #if np.sum(mask) > 0:
                                
                            da = {'wpath': os.path.expanduser('~/CEUAS/CEUAS/public/resort'),'version' : rversion, 'param': param,'v' : v, 'vname': vname, 'units': units[v],
                                  'rmsdevs': tr[ip,0, mask], 'interval' : interval, 'plev': plevs[ip], 'vtup': vtup, 'plevs': plevs, 'vnames': vname, 
                                  'lat': vlats[mask], 'lon': vlons[mask], 'fn': [], 'ranges' : rmsranges, 'title': d, 'sdict': sdict,'reportype': reptype,}
    
                            if d =='mean':
                                qsv = np.array((0, 0, 0, 0, 0.05, 0.5, 0.95))
                                nq = np.nanquantile(da['rmsdevs'], qsv)
                                da['qsv'] = qsv
                            else:    
                                nq = np.nanquantile(da['rmsdevs'], qsv)
                                da['qsv'] = qsv
                            if nq.size == qsv.size:
                                
                                alldevs[d]['allqs'].append(nq)
                            else:
                                alldevs[d]['allqs'].append(qsv+np.nan)
                            
                            if ip in pplevs:
                                #plt_rmsdevs(da)
                                futures.append(ray_plt_rmsdevs.remote(da))
                            if ip in pplevs:
                                #plt_rmsdevquantiles(da)
                                futures.append(ray_plt_rmsdevquantiles.remote(da))
                                
                            
            ray.get(futures)
            print('RMS plot', time.time()-tt)
            alldevs[d]['goods'] = np.concatenate(alldevs[d]['goods'])
            allqsdict = {}
            i = 0
            futures = []
            goods_ref = ray.put(alldevs[d]['goods'])
            allqs_ref = ray.put(alldevs[d]['allqs'])
            da['wpath'] = os.path.expanduser(f'~/CEUAS/CEUAS/public/resort/plots_rmsnew/{rversion}/')
            da['units'] = units
            da['params'] = params
            da['pplevs'] = np.arange(5, len(plevs), dtype=np.int32)
            for interval in rmsintervals :
                #quantileplot(alldevs[d]['allqs'], goods, da, i)
                da['interval'] = interval
                futures.append(ray_quantileplot.remote(allqs_ref, goods_ref, da, i))
                i += len(vtup) *3 * len(plevs)
                
            ray.get(futures)
        
            
                    
        
                        
        print('quantile plot', time.time()-tt)
        ray.internal.free(fdicts_ref, starts_ref)
        print('x')
        
        #spmean = np.zeros(12)
        #for im in range(1, 13):
            #d = {}
            #p = {'u': '131','v': '132',}
            #for c in 'u', 'v':
                
                #with h5py.File(f'/mnt/users/scratch/leo/scratch/era5/gridded/era5t.1990{im:0>2d}.{p[c]}.nc') as f:
                    #weight = np.cos(f['latitude'][:]*np.pi/180.) / np.mean(np.cos(f['latitude'][:]*np.pi/180.))
                    #d[c] = f[c][:, :, 0, 7]*f[c].attrs['scale_factor']+f[c].attrs['add_offset']
                    #print(im, c, f[c].shape,f['level'][7], np.mean(d[c]))
                
            #sp =np.mean(np.sqrt(d['u'] *d['u'] +d['v'] *d['v']), axis=1)
            #spmean[im-1] = np.mean(sp*weight)
        #print(f'mean wind speed: {np.mean(spmean):.2f}')
                
        
    return



sub(sdict)
# In[10]:

print(os.getcwd())

print(len(fns))
for ef in efiles.items():
    print(ef)


# In[ ]:





# In[ ]:


deps.keys(),ratios.keys(),d,var,ratios[d][var].shape


# In[ ]:


fns[-3:],len(fns)


# In[ ]:


for i in range(len(fdicts)):                
    for var in '137','138','126','139': 
        if var not in  fdicts[i].keys():
            continue
        plt.figure(figsize=(6,6))
        l=1
        for plev in plevs:
    #        idx=np.where((z==plev*100)&(~np.isnan(fg_dep))&(~np.isnan(fg_off)))[0]
    #        if len(idx)>50:
    #            print(f'{var} {plev}hPa: online {rms(fg_dep[idx]):.3f}, offline {rms(fg_off[idx]):.3f}, adj {rms(fg_off_corr[idx]):.3f}')
    #        else:
    #            continue
            try:
                #print(fdicts[i][var].keys())
                dt,fg_dep,fg_off,fg_off_corr=fdicts[i][var][plev]
            except Exception as e:
            #    print(fdicts['fn'], e)
                continue
            idx=np.arange(dt.shape[0])
    
            plt.subplot(5,1,l)
            plt.plot(dt[idx]/86400/365.25,fg_off_corr[idx],'.',label=f'off_corr {rms(fg_off_corr[idx]):.3f}')
            plt.plot(dt[idx]/86400/365.25,fg_off[idx],'.',label=f'off {rms(fg_off[idx]):.3f}, on {rms(fg_dep[idx]):.3f}')
            plt.plot(dt[idx]/86400/365.25,fg_off[idx]-fg_off_corr[idx],'.',label=f'bias {rms(fg_off[idx]-fg_off_corr[idx]):.3f}')
            yl=plt.gca().get_ylim()
            plt.text(dt[idx[1]]/86400/365.25,yl[0]+0.8*(yl[1]-yl[0]),f'{plev}hPa')
            plt.legend(bbox_to_anchor=(1, 1), loc='upper left', ncol=1)
            if l==1:
                plt.title(f'{var} {fdicts[i]['fn']}')
                
            l+=1
        plt.tight_layout()


# In[ ]:


plt.plot(dt2[:-1]/86400/365.25,dt2[1:]-dt2[:-1])
plt.plot(np.array([dt2[0],dt2[-1]])/86400/365.25,[0,0])
plt.ylim(-22000,22000)


# In[ ]:


idx=np.searchsorted(dt2,dt1)
idy=idx-1
plt.plot(dt2[1:]-dt2[:-1])
#idx[idx==dt2.shape[0]]-=1

#plt.plot(dt1,dt2[idx]-dt1)
mask=dt2[idy]-dt1<-7200
#idy=np.where(mask)[0]
#print(len(idy))
#plt.plot(dt1[mask]/86400/365.25,dt2[idx[mask]]-dt1[mask])
idz=np.where(dt2[idx[mask]]-dt1[mask]>7200)[0]
print(dt1[mask][idz]/86400/365.25)
#plt.plot(dt1[mask]/86400/365.25,dt2[idx[mask]]-dt1[mask])


# In[ ]:


fn=glob.glob(f"{paths['higra']}/{files['higra']}")[0]
with h5py.File(fn) as f:
    #print(f['source_configuration']['source_file'][0])
    files['htigra']=f['source_configuration']['source_file'][0].decode().strip()
    print(f"x{files['htigra']}x")
    
fn=glob.glob(f"{paths['he5']}/{files['he5']}")[0]
with h5py.File(fn) as f:
    #print(f['source_configuration']['source_file'][0])
    files['hte5']=f['source_configuration']['source_file'][0].decode().strip()
    files['hte5']=f'{YEAR}12'.join(files['hte5'].split('??????'))
    print(files['hte5'])
print(len('#AUM00011035 2023 12 30 12 1131   60 ncdc-gts           482486   163564, 21     0  99480'))
print(len('21 -9999  98270'))


# In[ ]:


fn=f"{paths['htigra']}/{files['htigra']}"
print(fn)
with open(fn) as f:
    data=f.read().split('\n')
    hlist=[]
    slist=[]
    i=0
    while i< len(data):
        if data[i]:
            if data[i][0]=='#':
                if YEAR==data[i][13:17]:
                    hlist.append(data[i])
                    i+=1
                    found=False
                    while data[i] and data[i][0]!='#':
                        try:
                            if int(data[i][9:15])==10000:
                                slist.append(data[i])
                                found=True
                        except Exception  as e:
                            print(i,data[i],e)
                            
                        i+=1
                    if not found:
                        hlist.pop()
                    #print(f'{hlist[-1]},\n{slist[-1]}')
        i+=1
vals={'htigra':{'datetime':[datetime(int(hlist[0][13:17]),int(hlist[0][18:20]),int(hlist[0][21:23]),int(hlist[0][24:26]))],
                                    'T':[int(slist[0][22:27])/10+273.15],'DPD':[int(slist[0][34:39])/10],'RH':[int(slist[0][28:33])/1000]}}
vals['htigra']['Td']=[vals['htigra']['T'][0]-vals['htigra']['DPD'][0]]
vals['htigra']['RHcalc']=[(Sonntag(np.array([vals['htigra']['Td'][0]]))/Sonntag(np.array([vals['htigra']['T'][0]])))[0]]
print(vals)


# In[ ]:


print(len(hlist),len(slist))
vals={'htigra':{'datetime':[],'T':[],'DPD':[],'RH':[],'Td':[],'RHcalc':[]}}
for i in range(0,len(hlist)):
    vals['htigra']['datetime'].append(datetime(int(hlist[i][13:17]),int(hlist[i][18:20]),int(hlist[i][21:23]),int(hlist[i][24:26])))
    if slist[i][22:27] in ['-8888','-9999']: 
        vals['htigra']['T'].append(np.nan)
    else:
        vals['htigra']['T'].append(int(slist[i][22:27])/10+273.15)
    if slist[i][34:39] in ['-8888','-9999']: 
        vals['htigra']['DPD'].append(np.nan)
    else:
        vals['htigra']['DPD'].append(int(slist[i][34:39])/10),
    if slist[i][28:33] in ['-8888','-9999']: 
        vals['htigra']['RH'].append(np.nan)
    else:
        vals['htigra']['RH'].append(int(slist[i][28:33])/1000)
    vals['htigra']['Td'].append(vals['htigra']['T'][i]-vals['htigra']['DPD'][i])
    vals['htigra']['RHcalc'].append((Sonntag(np.array([vals['htigra']['Td'][i]]))/Sonntag(np.array([vals['htigra']['T'][i]])))[0])

for k in vals['htigra'].keys():
    if k!='datetime':
        vals['htigra'][k]=np.array(vals['htigra'][k])
plt.figure(figsize=(12,5))
plt.subplot(1,2,1)
#plt.plot(vals['htigra']['datetime'],vals['htigra']['RH'])
plt.plot(vals['htigra']['datetime'],vals['htigra']['RHcalc'])
plt.plot(vals['htigra']['datetime'],vals['htigra']['RHcalc']-vals['htigra']['RH'])
plt.title('RH calc - RH, '+WIGOS+', '+YEAR+', 100 hPa')
#plt.ylim(0,1)
plt.subplot(1,2,2)
plt.plot(vals['htigra']['datetime'],vals['htigra']['DPD'])
plt.title('DPD, '+WIGOS+', '+YEAR+', 100 hPa')
mask=np.where(~np.isnan(vals['htigra']['RHcalc']))
print(np.corrcoef(np.array(vals['htigra']['RHcalc'])[mask],np.array(vals['htigra']['DPD'])[mask])[0,1])
mask=np.where(~np.isnan(vals['htigra']['RH']))
print(np.corrcoef(np.array(vals['htigra']['RH'])[mask],np.array(vals['htigra']['DPD'])[mask])[0,1])
#print(vals['htigra']['datetime'])


# In[ ]:


fns=glob.glob(f"{paths['hte5']}/{files['hte5']}")
for fn in fns[:1]:
    print(fn)
    with open('odcheader.txt') as f:
              data=f.read().split('\n')
    vars=[d.split(':')[1].split(',')[0].strip() for d in data[:-1]]
    zc=vars.index('vertco_reference_1@body')
    varnoc=vars.index('varno@body')
    obsc=vars.index('obsvalue@body')
    print(varnoc)
    with gzip.open(fn,'rt') as f:
        data=f.read().split('\n')
        hlist=[]
        slist=[]
        nobs=0
        for i in range(len(data)):
            if data[i]:
                l=data[i].split('\t')
                print('x'+l[varnoc]+'x'+str(l[obsc]))
                if l[varnoc]=='29' and float(l[zc])==10000.0:
                    hlist.append(float(l[obsc]))
                    #slist.append(data[i+1])
                    print(f'{hlist[-1]}')
                    nobs+=1
                    if(nobs>50):
                        break


# In[ ]:


def running_mean(x, N):
    nsum=np.cumsum(~np.isnan(np.insert(x, 0, 0)))
    x[np.isnan(x)]=0.
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / (nsum[N:]-nsum[:-N])


# In[ ]:


fn=glob.glob(f"{paths['resorted']}/{files['resorted']}")[0]
print(fn)
plevs=np.array([50,100,150,200,250,300,500,700])
plevs=np.array([200,300])
ref=datetime(1900,1,1)
names={'126':'T','137':'Td','138':'RH'}
def plot_station(fn):
    wigos=fn.split('/')[-1].split('_')[0]
    ref=datetime(1900,1,1)
    with h5py.File(fn) as f:
        #print(f['observations_table'].keys())
        #print(np.unique(f['observations_table']['observation_height_above_station_surface'][:]))
        #return
        trange=1995,2024
        rt=f['recordindices']['recordtimestamp'][:]
        hrec=f['header_table']['record_timestamp'][:]
        hrep=f['header_table']['report_timestamp'][:]

    plt.plot(hrec/86400/365.25,hrep-hrec)

plot_station(fn)


# In[102]:




# In[440]:


fn=glob.glob(f"{paths['resorted']}/{files['resorted']}")[0]
print(fn)
plevs=np.array([50,100,150,200,250,300,500,700])
plevs=np.array([200,300])
ref=datetime(1900,1,1)
names={'126':'T','137':'Td','138':'RH'}
def plot_station(fn):
    wigos=fn.split('/')[-1].split('_')[0]
    ref=datetime(1900,1,1)
    with h5py.File(fn) as f:
        #print(f['observations_table'].keys())
        #print(np.unique(f['observations_table']['observation_height_above_station_surface'][:]))
        #return
        trange=1990,2024
        rt=f['recordindices']['recordtimestamp'][:]
        idt=np.searchsorted(rt,((datetime(trange[0],1,1)-ref).total_seconds(),(datetime(trange[1],1,1)-ref).total_seconds()))
        names={'126':'T','137':'Td','138':'RH'}
        for par in '138','137',:
            tdslice=slice(f['recordindices'][par][idt[0]],f['recordindices'][par][idt[1]])
            sens=f['observations_table']['sensor_id'][tdslice].view('S4')[:,0]
            print(sens.shape)
            
            vtype,vi=np.unique(sens,return_index=True)
            print(vtype)
            #dpdslice=slice(f['recordindices']['34'][idt[0]],f['recordindices']['34'][idt[1]])
            try:
                hbe=f['advanced_homogenisation']['humidity_bias_estimate'][tdslice]
            except:
                hbe=np.zeros_like(f['observations_table']['observation_value'][tdslice])
                
            fgdep_o=f['era5fb']['fg_depar@offline'][tdslice]
            print('fgdep_o',np.nanstd(fgdep_o))
            obs=f['observations_table']['observation_value'][tdslice]
            print('minobs',par,np.nanmin(obs))
            if par=='137':
                obs[(obs-hbe)<160.]=np.nan
            ts=f['observations_table']['date_time'][tdslice]
            ztd=f['observations_table']['z_coordinate'][tdslice]
            #zdpd=f['observations_table']['z_coordinate'][tdslice]
            print(len(obs),len(ts))
            plt.figure(figsize=(16,5))
            ip=0
            for pl in plevs:
                idx=np.where((ztd==pl*100) & (sens!=b'nan') )[0]
                #print('bc',np.sum(~np.isnan(hbe[idx])),np.sum(~np.isnan(obs[idx])))
                plt.subplot(1,4,1)
                N=60
                rm=running_mean(obs[idx],N)
                plt.plot(ts[idx[N//2:-N//2+1]]/86400/365.25+1900,rm,label=f'{pl},{np.nanstd(obs[idx]):.3f}')
                plt.title(f'{wigos} {names[par]}')
                if ip==0:
                    for v in vi:
                        key=f['observations_table']['sensor_id'][tdslice][v].view('S4')[0].decode()
                        kk='?'
                        for k in sdict.keys():
                            if key in sdict[k]:
                                kk=k
                        if par=='137':
                            plt.text(ts[v]/86400/365.25+1900,ts[v]-ts[v]+np.min(rm)+0.05*(np.max(rm)-np.min(rm)),kk,rotation=90)
                        else:
                            plt.text(ts[v]/86400/365.25+1900,-0.05,kk,rotation=90)
                            plt.plot(ts[v]/86400/365.25+1900,-0.06,'ok')
                #plt.ylim(np.min(rm),np.max(rm))
                ylimunadj=plt.gca().get_ylim()
                plt.legend(loc='upper right')
                plt.subplot(1,4,2)
                rm=running_mean(fgdep_o[idx],N)
                plt.plot(ts[idx[N//2:-N//2+1]]/86400/365.25+1900,rm,label=f'fgdep, {pl},{np.nanstd(fgdep_o[idx]):.3f}')
                plt.title(f'{wigos} {names[par]} obs-ERA5bg')
                plt.legend()
                ylim=plt.gca().get_ylim()
                plt.subplot(1,4,3)
                #idx=np.where(ztd==pl*100)[0]
                rm=running_mean(hbe[idx],N)
                #print('bc',np.sum(~np.isnan(hbe[idx])),np.sum(~np.isnan(rm)))
                plt.plot(ts[idx[N//2:-N//2+1]]/86400/365.25+1900,rm,label=f'bias {pl},{np.nanstd(hbe[idx]):.3f}')
                #ylim=plt.gca().get_ylim()
                #print(ylim)
                #if np.max(np.abs(ylim))<0.1:
                #plt.ylim(ylim)
                plt.legend()
                plt.title(f'{wigos} {names[par]} bias')
                plt.subplot(1,4,4)
                rm=running_mean(obs[idx]-hbe[idx],N)
                print('obs-hbe',np.nanmin(obs[idx]-hbe[idx]),np.nanmin(rm))
                plt.title(f'{wigos} {names[par]} adjusted')
                plt.plot(ts[idx[N//2:-N//2+1]]/86400/365.25+1900,rm,label=f'{pl},{np.nanstd(obs[idx]-hbe[idx]):.3f}')
                plt.ylim(ylimunadj)
                plt.legend(loc='upper right')
                ip+=1

plot_station(fn)
    #print(f['observations_table']['z_coordinate'][idx])


# In[431]:


rversion='25'
resorted=f'/mnt/users/scratch/leo/scratch/converted_v{rversion}/long/'
WIGOS='0-20000-0-72357'
fns=glob.glob(f"{resorted}/{WIGOS}*.nc")
fn=fns[0]
print(fn)
with h5py.File(fn) as f:
    trange=1995,2024
    rt=f['recordindices']['recordtimestamp'][:]
    idt=np.searchsorted(rt,((datetime(2006,7,3,11)-ref).total_seconds()))
    print(idt,ref+timedelta(seconds=int(rt[idt])))
    tslice=slice(f['recordindices']['126'][idt],f['recordindices']['126'][idt+1])
    rhslice=slice(f['recordindices']['138'][idt],f['recordindices']['138'][idt+1])
    tdslice=slice(f['recordindices']['137'][idt],f['recordindices']['137'][idt+1])
    zt=f['observations_table']['z_coordinate'][tslice]
    zrh=f['observations_table']['z_coordinate'][rhslice]
    iz=np.searchsorted(zt,zrh)
    fbt=f['era5fb']['obsvalue@body'][tslice][iz]
    fbtdep=f['era5fb']['fg_depar@body'][tslice][iz]
    fbrh=f['era5fb']['obsvalue@body'][rhslice]
    fbrhdep=f['era5fb']['fg_depar@body'][rhslice]
    print(len(fbt),len(fbrh))
    vpdata = (fbrh-fbrhdep)* np.exp(liquid(fbt-fbtdep)) 
            #vpdata2 = fobs[mask, 5]* np.exp(ice(fobs[mask, 2]))
    fbtd=np.empty_like(fbrh)
    obstd=np.empty_like(fbrh)
    for i in range(len(fbtd)):
        fbtd[i] = bisection(func, fbt[i]-fbtdep[i]-70., fbt[i]-fbtdep[i]+1., vpdata[i])
    if not np.any(~np.isnan(f['era5fb']['obsvalue@body'][tdslice])):
        vpdata = (fbrh)* np.exp(liquid(fbt)) 
        for i in range(len(fbtd)):
            obstd[i] = bisection(func, fbt[i]-70., fbt[i]+1., vpdata[i])
        
        fbtddep=obstd-fbtd
    else:
        obstd=f['era5fb']['obsvalue@body'][tdslice]
        fbtddep=obstd-fbtd
    #print(fbtddep)
    plt.figure(figsize=(12,5))
    spd={'126':1,'137':1,'138':2}
    for par in '126','137','138':
        
        tdslice=slice(f['recordindices'][par][idt],f['recordindices'][par][idt+1])
        print(tdslice)
        plt.subplot(1,2,spd[par])
        plt.semilogy(f['observations_table']['observation_value'][tdslice],f['observations_table']['z_coordinate'][tdslice],label=names[par]+' obs')
        if par!='137':
            plt.semilogy(f['era5fb']['obsvalue@body'][tdslice]-f['era5fb']['fg_depar@body'][tdslice],
                         f['observations_table']['z_coordinate'][tdslice],'-x',label=names[par]+' ERA5')
        else:
            plt.semilogy(obstd-fbtddep,
                         f['observations_table']['z_coordinate'][tdslice],'-x',label=names[par]+' ERA5')
        #plt.semilogy(f['observations_table']['observation_value'][tdslice]-f['advanced_homogenisation']['humidity_bias_estimate'][tdslice],
        #             f['observations_table']['z_coordinate'][tdslice],'-o',label=names[par]+' adjusted obs')
        #plt.semilogy(f['era5fb']['obsvalue@body'][tdslice]-f['advanced_homogenisation']['humidity_bias_estimate'][tdslice],
        #             f['observations_table']['z_coordinate'][tdslice])
        plt.ylim(101000,5000)
        plt.ylabel('pressure [Pa]')
        plt.title(f'{WIGOS}, {names[par]},{(ref+timedelta(seconds=int(rt[idt])))}')
        plt.legend()


# In[12]:


fns=glob.glob(f"{paths['resorted']}/*0-71*v3.nc")
print(len(fns))
vai=[]
for r in sdict.items():
    if 'rs' in r[0]:
        vai+=r[1]
print(vai)
for fn in fns:
    with h5py.File(fn,'r') as f:
        rt=f['recordindices']['recordtimestamp'][:]
        idt=np.searchsorted(rt,((datetime(1995,1,1,0)-ref).total_seconds()))
        par='126'
        try:
            tdslice=slice(f['recordindices'][par][idt],f['recordindices'][par][-1])
            idx=np.where(f['observations_table']['z_coordinate'][tdslice]==50000)[0][::-1] # we want last occurence
            vtype,vi,vc,=np.unique(f['observations_table']['sensor_id'][tdslice][idx].view('S4'),return_index=True,return_counts=True)
            ti=f['observations_table']['date_time'][tdslice][idx][vi]
            #print(ri.shape)
            i=0
            lastvai={'time':-1,'count':0,'type':''}
            lastother={'time':-1,'count':0,'type':''}
            for v in vtype:
                vv=v.decode()
                if vv in vai:
                    if ti[i]>lastvai['time']:
                        lastvai['time']=ti[i]
                        lastvai['count']=vc[i]
                        lastvai['type']=vv
                else:
                    if ti[i]>lastother['time']:
                        lastother['time']=ti[i]
                        lastother['count']=vc[i]
                        lastother['type']=vv
                i+=1
            for l in lastvai,lastother:
                if l['time']>-1 and l['type']!='nan':
                    print(fn.split('/')[-1],l['type'],l['count'],ref+timedelta(seconds=int(l['time'])))
            #break
        except KeyError as e:
            #print(fn.split('/')[-1],e)
            pass
            




# In[13]:


from tqdm import tqdm
fns=glob.glob(f"{paths['resorted']}/*v3.nc")
print(len(fns))
vai=[]
for r in sdict.items():
    if 'rs' in r[0]:
        vai+=r[1]
print(vai)

@njit
def calc_tds(fbt,vpdata):
    td=np.empty_like(fbt)
    for i in range(len(fbt)):
        td[i] = bisection(func, fbt[i]-70., fbt[i]+1., vpdata[i])
    return td

def ptt(tt,label=''):
    print(f'{label} {time.time()-tt:.4f}')
    
def vaisala(rtype,vlist=[]):
    return rtype in vlist or rtype[0]=='V'
    
def nonvaisala(rtype,vlist=[]):
    return rtype not in vlist and rtype[0]!='V'
    
def all (rtypes,vlist=[]):
    return True

    
def readstats(fn,plev=150,selector=vaisala,vlist=vai,startyear=1990,endyear=2015,min_ascents=50,debug=False): #sdict['rs90']+sdict['rs92']+sdict['rs41'] #sdict['rs92']
#for fn in fns: #tqdm(fns[:300]):
    stats=dict(vmeans=[],vfgdepmeans=[],vfgdepadjmeans=[],vfgdeprms=[],vfgdepadjrms=[],vfn=[],vlat=[],vlon=[])
    wigos=fn.split('/')[-1].split('_')[0]
    ref=datetime(1900,1,1)
    tt=time.time()
    with h5py.File(fn,'r') as f:
        par='137'
        if par not in f['recordindices'].keys():
            return {}
        rt=f['recordindices']['recordtimestamp'][:]
        idts=np.searchsorted(rt,((datetime(startyear,1,1,0)-ref).total_seconds()))
        idte=np.searchsorted(rt,((datetime(endyear+1,1,1,0)-ref).total_seconds()))
        try:
            tdslice=slice(f['recordindices']['137'][idts],f['recordindices']['137'][idte])
            tslice=slice(f['recordindices']['126'][idts],f['recordindices']['126'][idte])
            rhslice=slice(f['recordindices']['138'][idts],f['recordindices']['138'][idte])
            tidx=np.where(f['observations_table']['z_coordinate'][tslice]==plev*100)[0] # we want last occurence
            tdidx=np.where(f['observations_table']['z_coordinate'][tdslice]==plev*100)[0] # we want last occurence
            rhidx=np.where(f['observations_table']['z_coordinate'][rhslice]==plev*100)[0] # we want last occurence
            #ptt(tt,wigos)
            ttime=f['observations_table']['date_time'][tslice][tidx]
            rhtime=f['observations_table']['date_time'][rhslice][rhidx]
            #tdtime=f['observations_table']['date_time'][tdslice][tdidx]
            #if len(rhtime)<=len(tdtime):
            it=np.searchsorted(ttime,rhtime)
            #else:
            #    it=np.searchsorted(ttime,tdtime)
                
            notfound=np.any(it==ttime.shape[0])
            if notfound:
                ity=np.where(it<ttime.shape[0])[0]
                print(wigos,'rh time not found',len(ity),len(rhtime))
                if len(ity)>0:
                    rhidx=rhidx[ity]
                    tdidx=tdidx[ity]
                    it=it[ity]
                #return {}
            #print(len(ttime),len(rhtime),it,tslice,rhslice,idt)
            sens=f['observations_table']['sensor_id'][tdslice][tdidx].view('S4')
            #print(np.unique(sens,return_counts=True))
            adj=f['advanced_homogenisation']['humidity_bias_estimate'][tdslice][tdidx]
            rh=f['observations_table']['observation_value'][rhslice][rhidx]
            rhfg=rh-f['era5fb']['fg_depar@body'][rhslice][rhidx]
            t=f['observations_table']['observation_value'][tslice][tidx][it]
            tfg=t-f['era5fb']['fg_depar@body'][tslice][tidx][it]
            #ptt(tt,wigos)
            vpdata = (rhfg)* np.exp(liquid(tfg))
            tdfg=calc_tds(tfg,vpdata)
            vpdata = (rh)* np.exp(liquid(t))
            tdobs=calc_tds(t,vpdata)
            tdfgdep=tdfg-tdobs
            #ptt(tt,wigos+' bisect')
            
            vadj=[]
            vtdfgdep=[]
            vtdfgdepadj=[]
            vi=[]
            for i in range(len(tdidx)):
                if selector(sens[i][0].decode(),vlist): #sdict['rs90']+or sens[i][0].decode()[:2]=='VN' or sens[i][0].decode()[:2]=='V9':
                    try:
                        vtdfgdepadj.append(tdfgdep[i]-adj[i])
                        vtdfgdep.append(tdfgdep[i])
                        vadj.append(adj[i])
                        vi.append(i)
                    except:
                        #print('f',i,len(tdidx))
                        pass
            #ptt(tt,wigos+' select vaisala')
            if len(vadj)>min_ascents:
                nm=np.nanmean(vadj)
                mtdfgdep=np.nanmean(vtdfgdep)
                rmstdfgdep=rms(np.array(vtdfgdep))
                madjtdfgdep=np.nanmean(vtdfgdepadj)
                rmsadjtdfgdep=rms(np.array(vtdfgdepadj))
                if debug:
                    vi=np.array(vi)
                    print(np.unique(sens))
                    plt.plot(rhtime/86400/365.25+1900,tdobs)
                    plt.plot(rhtime/86400/365.25+1900,t)
                    plt.plot(rhtime[vi]/86400/365.25+1900,vtdfgdep,'o',label=f'fgdep,{len(vi)},{mtdfgdep:.3f}')
                    plt.legend()
                if not np.isnan(nm):
                    stats['vmeans'].append(nm)
                    stats['vfgdepmeans'].append(mtdfgdep)
                    stats['vfgdepadjmeans'].append(madjtdfgdep)
                    stats['vfgdeprms'].append(rmstdfgdep)
                    stats['vfgdepadjrms'].append(rmsadjtdfgdep)
                    
                    stats['vfn'].append(wigos)
                    stats['vlat'].append(f['header_table']['latitude'][0])
                    stats['vlon'].append(f['header_table']['longitude'][0])
                #print(f"{fn.split('/')[-1]}: Mean Td bias estimate, {np.nanmean(vadj):.2f}, {len(vadj)}")   
        except Exception as ex:
            if type(ex).__name__!='KeyError':
                print(type(ex).__name__, ex.args)
            return {}
            #vtype,vi,vc,=np.unique(f['observations_table']['sensor_id'][tdslice][idx].view('S4'),return_index=True,return_counts=True)
            #ti=f['observations_table']['date_time'][tdslice][idx][vi]
            #svtype=[v.decode() for v in vtype]
        return stats
            


# In[14]:


from multiprocessing import Pool
P=Pool(40)
statslist=list(P.map(readstats,fns[:]))
statslist=[s for s in statslist if s and s['vmeans']]
stats={}
stats['vaisala']={}
if statslist:
    stats['vaisala'][150]=statslist[0]
    for k in stats['vaisala'][150].keys():
        stats['vaisala'][150][k] =stats['vaisala'][150][k]+[s[k][0] for s in statslist[1:]]
print(len(stats['vaisala'][150]['vmeans']))


# In[15]:


stats={}
stats['vaisala']={}
if statslist:
    stats['vaisala'][150]=statslist[0]
    for k in stats['vaisala'][150].keys():
        stats['vaisala'][150][k] =stats['vaisala'][150][k]+[s[k][0] for s in statslist[1:]]
print(len(stats['vaisala'][150]['vmeans']))


# In[16]:


def plstats(stats,idx,title=''):
    for k in stats.keys():
        stats[k]=np.array(stats[k])
    plt.figure(figsize=(12,5))
    plt.subplot(1,2,1)
    plt.plot(stats['vmeans'][idx],label=f'adj {rms(stats['vmeans'][idx]):.3f}')
    plt.plot(stats['vfgdepmeans'][idx],label=f'fgdep {rms(stats['vfgdepmeans'][idx]):.3f}')
    plt.plot(stats['vfgdepadjmeans'][idx],label=f'fgdepadj {rms(stats['vfgdepadjmeans'][idx]):.3f}')
    plt.title(title)
    plt.legend()
    plt.subplot(1,2,2)
    plt.plot(stats['vfgdeprms'][idx],label=f'fgdep {rms(stats['vfgdeprms'][idx]):.3f}')
    plt.plot(stats['vfgdepadjrms'][idx],label=f'fgdepadj {rms(stats['vfgdepadjrms'][idx]):.3f}')
    plt.legend()
plstats(stats['vaisala'][150],slice(0,len(stats['vaisala'][150]['vmeans'])),title='Means and rms departures, All Vaisala sondes')


# In[17]:


plt.hist(stats['vaisala'][150]['vmeans'],bins=20)


# In[18]:


idy=np.where(np.array(stats['vaisala'][150]['vmeans'])<-0.5)[0]
plstats(stats['vaisala'][150],idy,title='Vaisala sondes with adjustments <=-0.5K at 150 hPa')
idy=np.where(np.array(stats['vaisala'][150]['vmeans'])>-0.5)[0]
plstats(stats['vaisala'][150],idy,title='Vaisala sondes with adjustments >-0.5K at 150 hPa')

idz=np.where(np.abs(stats['vaisala'][150]['vmeans'][idy])>10)
problematic=np.array(stats['vaisala'][150]['vfn'])[idy][idz]
print(len(idy),problematic)


# In[19]:


glob.glob(f"{paths['resorted']}/{files['resorted']}")[0]
problematic='0-20000-0-94120',
for w in problematic:
    fn=f"{paths['resorted']}/{w}_CEUAS_merged_v3.nc"
    plot_station(fn)


# In[24]:


from functools import partial
rs80 =['52','37',  ] # '60', '61','62','63' ,  '66', '67'

fns=[ fn for fn in fns if '47122' not in fn and '06458' not in fn]
#fns=[fn for fn in fns if '0-94120' in fn]
stats={}
for sele in vaisala,: #vaisala,nonvaisala,all:
    sn=sele.__name__
    stats[sn]={}
    for plev in 150,: #100,150,200,250,300:
        vfunc=partial(readstats,selector=sele,vlist=rs80+sdict['rs90']+sdict['rs92']+sdict['rs41'] ,plev=plev,startyear=1979,endyear=2015,min_ascents=30) #
        statslist=list(P.map(vfunc,fns[:]))
        statslist=[s for s in statslist if s and s['vmeans']]
        if statslist:
            stats[sn][plev]=statslist[0]
            for k in stats[sn][plev].keys():
                stats[sn][plev][k] =stats[sn][plev][k]+[s[k][0] for s in statslist[1:]]
        print(sn,plev,len(stats[sn][plev]['vfgdepadjrms']))
    


# In[21]:


for k in stats.keys():
    print(k,stats[k].keys())


# In[29]:


gstats=[]
for thresh in 0.5,:
    for sele in vaisala,: #vaisala,nonvaisala,all:
        sn=sele.__name__
        for plev in 150,:#stats[sn].keys():
            idy=np.where(np.array(stats[sn][plev]['vfgdepadjrms'])-np.array(stats[sn][plev]['vfgdeprms'])<thresh)[0]
            #idy=np.where(np.array(stats[sn][plev]['vmeans'])>thresh)[0]
            plstats(stats[sn][plev],idy,title=f'{sn} sondes with rms difference\n due to adjustments <={thresh}K at {plev} hPa')
            
            #idy=np.where(np.array(stats[sn][plev]['vmeans'])>-1000.)[0]
            #plstats(stats[sn][plev],idy,title=f'{sn} sondes with adjustments >{thresh}K at {plev} hPa')
            
            idz=np.where(stats[sn][plev]['vfgdepadjmeans'][idy]-stats[sn][plev]['vfgdepmeans'][idy]>3)[0]
            problematic=np.array(stats[sn][plev]['vfn'])[idy][idz]
            print(len(idy),len(idz))
            print(problematic)
lstr='['
l=0
for i in range(len(stats[sn][plev]['vfn'])):
    lstr+='"'+stats[sn][plev]['vfn'][i]+'",'
    l+=1
    if l%6==5:
        lstr+='\n'
lstr+=']'
print(lstr)
import cartopy.crs as ccrs
plt.figure()
ax=plt.subplot(1,1,1,projection=ccrs.PlateCarree())
ax.coastlines(lw=0.5)
c=plt.scatter(stats[sn][plev]['vlon'],stats[sn][plev]['vlat'],transform=ccrs.PlateCarree(),s=9)
plt.title('Vaisala sondes with working adjustments')

#ax=plt.subplot(1,1,1,projection=ccrs.PlateCarree())
#c=plt.scatter(lons,lats,transform=ccrs.PlateCarree(),s=9,c=ps)


# In[ ]:


stats[sn][plev].keys()


# In[23]:


all.__name__


# In[ ]:




