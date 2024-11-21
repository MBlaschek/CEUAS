# Compute entropy of a system

import numpy as np
import pandas as pd
from scipy.special import xlogy
import matplotlib.pyplot as plt; plt.ion()
import matplotlib.colors as mcolors
import pyproj; geoid = pyproj.Geod(ellps="WGS84")
import itertools
import cartopy.crs as ccrs

# System is defined by:
# - axes (index of a dataframe)
# - expected deviation for each axis

#entropies = None


def create_distance_matrix_from_dimension(xvals0): #,ijs0,jis0):
  #print("we are here 1")
  xvals1, xvals2 = np.meshgrid( xvals0, xvals0 )
  #print("we are here 2")
  Dist_matrix0 = np.abs(xvals1 - xvals2)
  #print("we are here 3")
  #print(xvals0)
  return Dist_matrix0

# Compute entropy from a dataframe
def entropy(df,dims,mysigmas,gains=None):
  n = len(df)
  Dmatrix = np.zeros((n,n,))
  ndims = len(dims)-1
  for idim1, dim1 in enumerate(dims):
    Dist_matrix = create_distance_matrix_from_dimension(df[dim1].values[:])
    Dmatrix = Dmatrix + (Dist_matrix / mysigmas[idim1])
  Dmatrix = Dmatrix / float(ndims)
  Smatrix = np.exp( - Dmatrix / mysigmas[-1] )
  d1 = xlogy(Smatrix, Smatrix)
  d2 = xlogy(1.-Smatrix, 1.-Smatrix)
  Hmatrix = (-1./np.log(2.)) * ( d1 + d2 ) # this is the binary entropy -- maximum value is 0.5
  E = np.sum(Hmatrix) / (float(n)**2.)
  if gains is not None:
    entropy_withouts = []
    # remove one observation each time
    if n>2:
      if type(gains) is int:
        gain_list = [gains]
      elif type(gains) is slice:
        gain_list = list(range(gains.stop)[gains])
      else:
        gain_list = list(gains)
        stop
      for j in gain_list:
        if j==0:
          entropy_without1 = np.sum(Hmatrix[1:,1:])
        elif j==n-1:
          entropy_without1 = np.sum(Hmatrix[:-1,:-1])
        else:
          entropy_without1 =  np.sum(Hmatrix[   :j,   :j])
          entropy_without1 += np.sum(Hmatrix[j+1: ,j+1: ])
          entropy_without1 += np.sum(Hmatrix[   :j,j+1: ])
          entropy_without1 += np.sum(Hmatrix[j+1: ,   :j])
        entropy_withouts.append(entropy_without1)
      entropy_gains = (np.array(entropy_withouts) /(float(n-1)**2.)) * (-1.)
      entropy_gains += E
      return E, entropy_gains
    else:
      return E, None
  else:
    return E

def function_entropy(df0, index_center_call0, index_names0, entropy_dims0, sigmas0, dfent0):

  # count deciles
  def count_deciles(vals0):
    decile_counts = []
    for k in range(11):
      kmin=float(k-1)*0.1
      kmax=float(k)*0.1
      decile_counts.append(np.count_nonzero((vals0>=kmin) & (vals0<kmax)))
    return decile_counts

  try:
    iindex = dfent0.index.get_loc(index_center_call0)
  except:
    iindex = -1
  if iindex!=-1:
    # total entropy
    E, entropy_gains = entropy(df0,entropy_dims0,sigmas0,gains=df0.set_index(index_names0).index.get_loc(index_center_call0))
    dfent0.iloc[iindex,0] = E # Assign entropy value
    dfent0.iloc[iindex,1] = len(df0) # Assign number of elements
    if entropy_gains is not None:
      dfent0.iloc[iindex,2] = entropy_gains
      #decs = count_deciles(entropy_gains)
      #for idec1,dec1 in enumerate(decs):
      #  dfent0.iloc[iindex,3+idec1] = dec1
    if np.count_nonzero(dfent0['nobs']>0) % 1000 == 0:
      print("+1000")

def get_RS_profile_value(df0,myvars0,orphvar0,sourcevar0):

  def rm_var(arr,rmvars):
    return [x for x in arr if x not in rmvars]

  # Return a "value" for the profile -- so the results do not depend on the order in which the comparisons are made
  n = len(df0)
  conds = {}
  for var1 in myvars0:
    #print(var1)
    conds[var1] = (~np.isnan(df0[var1]))
    #print(conds[var1])
  prev_cond = np.full( (n,), False)
  nvalue = 0
  for combvars in [myvars0, rm_var(myvars0, ['dewpointTempBiasCorr']), rm_var(myvars0, ['airTempBiasCorr','dewpointTempBiasCorr']), rm_var(myvars0, ['airTempBiasCorr','dewpointTempBiasCorr','dewpointTemperature']), rm_var(myvars0, ['airTempBiasCorr','dewpointTempBiasCorr','dewpointTemperature','airTemperature'])]:
    nvars = len(combvars)
    #print(combvars)
    cond_0 = np.full( (n,), True)
    for combvar1 in combvars:
      cond_0 = cond_0 & conds[combvar1]
    cond_app = (~prev_cond) & cond_0
    ngood = np.count_nonzero(cond_app)
    #print("ngood",ngood)
    if ngood>0:
      nvalue += ngood*(nvars**2)
    prev_cond = (prev_cond | cond_app)
    #print("new value",nvalue)
  lorph = np.count_nonzero(np.char.find(df0[orphvar0].values[:].astype(str),'orph')>=0)>0
  if lorph:
    nvalue -= 100 # penalize such profiles
  lncar = np.count_nonzero(np.char.find(df0[sourcevar0].values[:].astype(str),'ncar')>=0)>0
  if lncar:
   nvalue -= 10 # penalize such profiles
  ligra = np.count_nonzero(np.char.find(df0[sourcevar0].values[:].astype(str),'igra')>=0)>0
  if ligra:
   nvalue += 15 # give such profiles more points
  #print("orphan!",nvalue)
  #print("final value",nvalue)
  return nvalue

def RobustStandardDeviation(myarray):
  # compute robust standard deviation
  mask_ok = np.where(~np.isnan(myarray))[0]
  if mask_ok.size==0:
    return np.nan
  else:
    try:
      myarray_val = myarray[mask_ok]
    except:
      myarray_val = myarray.values[mask_ok]
    mymedian_val = np.median(myarray_val)
    RobStdDev = np.median( np.abs(myarray_val-mymedian_val) ) * 1.4826
    return RobStdDev

def PercentileAbsoluteDeviation(myarray1, myarray2, q1):
  # compute Percentile Absolute Deviation
  mask_ok = np.where((~np.isnan(myarray1)) & (~np.isnan(myarray2)))[0]
  if mask_ok.size==0:
    return np.nan
  else:
    try:
      myarray_val1 = myarray1[mask_ok]
    except:
      myarray_val1 = myarray1.values[mask_ok]
    try:
      myarray_val2 = myarray2[mask_ok]
    except:
      myarray_val2 = myarray2.values[mask_ok]
    mymedian_val1 = np.median(myarray_val1)
    mymedian_val2 = np.median(myarray_val2)
    PercAbsDev = np.percentile( np.abs( (myarray_val1-mymedian_val1) - (myarray_val2-mymedian_val2) ), q1 )
    return PercAbsDev

def get_DRIFTER_trajectory_value(df0,myvars0,bounds0):

  # Return a "value" for the trajectory -- so the results do not depend on the order in which the comparisons are made
  nvalue = 0
  ivar = len(myvars0)
  for var1 in myvars0:
    if type(var1) is list:
      wp = np.full( (len(df0),), True)
      for var0 in var1:
        wok = wp | ( (~np.isnan(df0[var0].values[:])) & (df0[var0].values[:]>bounds0[var0][0]) & (df0[var0].values[:]<bounds0[var0][1]))
      nvalue += np.count_nonzero(wok) * (ivar**2)
    else:
      nvalue += np.count_nonzero((~np.isnan(df0[var1].values[:])) & (df0[var1].values[:]>bounds0[var1][0]) & (df0[var1].values[:]<bounds0[var1][1])) * (ivar**2)
    ivar -= 1
  return nvalue
  
def function_sameprofileas(df0, index_center_call0, index_names0, pairs_done0, discarded_list0, df_dup0):

  #try:
  #  iindex = dfsameas0.index.get_loc(index_center_call0)
  #except:
  #  iindex = -1
  #if iindex!=-1:
  if 1==1:
    df0_profiles = df0.set_index(['wigos_statid','report_id']).index.unique()
    nprofiles = len(df0_profiles)
    if (nprofiles>1):
      #print("nprofiles",nprofiles)
      #vals1 = df0[index_names0[0]].values[:]
      #cond = np.isin(dfsameas0.index.get_level_values(index_names0[0]), vals1)
      #for index_names01 in index_names0[1:]:
      #  vals1 = df0[index_names01].values[:]
      #  cond = cond & np.isin(dfsameas0.index.get_level_values(index_names01), vals1)
      #iindex = np.where(cond) [0]
      #if (len(df0)!=len(dfsameas0.iloc[iindex])):
      #  print(df0)
      #  print(dfsameas0.iloc[iindex])
      #  stop
      #sameas_vals = dfsameas0.iloc[iindex]['sameas'].values[:]
      #discard_vals= dfsameas0.iloc[iindex]['discard'].values[:]
      #print(df0indexed.groupby(df0indexed.index.names).count()[['hour']].rename(columns={'hour':'npoints'}))
      # now pick one pair at a time and compute distance
      for ijpair in itertools.combinations(range(nprofiles),2):
        statid1,statid2 = [x for x in sorted( [df0_profiles[ijpair[0]], df0_profiles[ijpair[1]]] )] # note SORTING
        if (statid1 in discarded_list0) or (statid2 in discarded_list0):
          continue # Data were found to be inferior to another profile, so skip
        if (statid1,statid2) in pairs_done0:
          continue # Each pair of profiles only needs to be checked once
        pairs_done0.append((statid1,statid2))
        statid_pair = [statid1, statid2]
        dfs = []
        for ipos in [0,1]:
          iindex = np.where((df0['wigos_statid'].values[:]==statid_pair[ipos][0]) & (df0['report_id'].values[:]==statid_pair[ipos][1]))[0]
          dfs.append(df0.iloc[iindex].copy())
        if len(dfs[0])<3 or len(dfs[1])<3:
          continue # we do not have much to deal with... too few points!
        # are the data both from the ERA5 source?
        lbothERA5_1 = (dfs[0]['source_id'].values[0]==b'era5_1') and (dfs[1]['source_id'].values[0]==b'era5_1')
        #iindex1 = np.where((df0['wigos_statid'].values[:]==df0_profiles[ijpair[0]][0]) & (df0['hour'].values[:]==df0_profiles[ijpair[0]][1]))[0]
        #iindex2 = np.where((df0['wigos_statid'].values[:]==df0_profiles[ijpair[1]][0]) & (df0['hour'].values[:]==df0_profiles[ijpair[1]][1]))[0]
        ##if (type(iindex1) is not slice) or (type(iindex1) is not slice):
        ##  continue # at least one of these two is a single-level profile...
        #n1 = len(discard_vals[iindex1])
        #n2 = len(discard_vals[iindex2])
        #if n1<3 or n2<3:
        #  continue # we do not have much to deal with... too few points!
        #if list(np.unique(discard_vals[iindex1]))!=[0] or list(np.unique(discard_vals[iindex2]))!=[0]:
        #  continue # Data were found to be inferior to another profile, so skip
        cols_compare = ['airTemperature','dewpointTemperature','windSpeed','windDirection','airTempBiasCorr','dewpointTempBiasCorr','wigos_statid','source_id']
        #df1 = df0indexed.loc[df0_profiles[ijpair[0]]].reset_index().set_index(['pressure'])[cols_compare]
        #df2 = df0indexed.loc[df0_profiles[ijpair[1]]].reset_index().set_index(['pressure'])[cols_compare]
        df1 = dfs[0].set_index(['pressure'])[cols_compare]
        df2 = dfs[1].set_index(['pressure'])[cols_compare]
        #df1.to_pickle("df1.pkl")
        #df2.to_pickle("df2.pkl")
        if type(df1) is not pd.DataFrame or type(df2) is not pd.DataFrame:
          continue
        #print(df1)
        #print(df2)
        # COLUMNS TO CONSIDER: ['hour', 'lat01i', 'lon01i', 'wigos_statid', 'pressure', 'lat', 'lon', 'airTemperature', 'airTempBiasCorr', 'sonde_type', 'source_id', 'height', 'platform_type', 'dewpointTemperature', 'dewpointTempBiasCorr', 'windSpeed', 'windSpeedBiasCorr', 'latitude_displacement', 'longitude_displacement', 'windDirection', 'windDirectionBiasCorr', 'geopotential', 'nonCoordinateGeopotentialHeight', 'extendedVerticalSoundingSignificance', 'lat05', 'lon05']
        # start with temperature
        dfinn = pd.merge(df1,df2,how='inner',left_index=True,right_index=True,suffixes=('_1','_2'))
        #print(dfinn[['airTemperature_1','airTemperature_2']])
        #print(dfinn[['windSpeed_1','windSpeed_2']])
        # start with temperature
        var0 = 'airTemperature_{0}'; var1=var0.format(1); var2=var0.format(2)
        wtemp = np.where((~np.isnan(dfinn[var1].values[:])) & (~np.isnan(dfinn[var2].values[:])))[0]; ntemp=len(wtemp)
        if ntemp>0:
          temp1 = np.around(dfinn.iloc[wtemp][var1].values[:],2)
          temp2 = np.around(dfinn.iloc[wtemp][var2].values[:],2)
          ntempsame = np.count_nonzero(np.isin(np.around(np.abs(temp1-temp2),2),[0.,0.01,0.05,0.1]))
          ntempdiff = ntemp-ntempsame
        else:
          ntempsame = 0
        var0s = 'windSpeed_{0}';     var1s=var0s.format(1); var2s=var0s.format(2)
        var0d = 'windDirection_{0}'; var1d=var0d.format(1); var2d=var0d.format(2)
        wwind = np.where((~np.isnan(dfinn[var1s].values[:])) & (~np.isnan(dfinn[var2s].values[:])) & (~np.isnan(dfinn[var1d].values[:])) & (~np.isnan(dfinn[var2d].values[:])))[0]; nwind=len(wwind)
        if nwind>0:
          wnds1 = np.around(dfinn.iloc[wwind][var1s].values[:],1)
          wnds2 = np.around(dfinn.iloc[wwind][var2s].values[:],1)
          wndd1 = np.around(dfinn.iloc[wwind][var1d].values[:],0)
          wndd2 = np.around(dfinn.iloc[wwind][var2d].values[:],0)
          nwindsame_msms  = np.count_nonzero(((np.around(wnds1,0)>0) | (np.around(wnds2,0)>0)) & np.isin(np.around(np.abs(wnds1-wnds2),0),[0.]) & np.isin(np.around(np.abs(wndd1-wndd2),0),[0.,5.,10.,350.,355.,360.]))
          nwindsame_knot1 = np.count_nonzero(((np.around(wnds1,0)>0) | (np.around(wnds2,0)>0)) & np.isin(np.around(np.abs(wnds1*1852/3600.-wnds2),0),[0.]) & np.isin(np.around(np.abs(wndd1-wndd2),0),[0.,5.,10.,350.,355.,360.]))
          nwindsame_knot2 = np.count_nonzero(((np.around(wnds1,0)>0) | (np.around(wnds2,0)>0)) & np.isin(np.around(np.abs(wnds1-wnds2*1852/3600.),0),[0.]) & np.isin(np.around(np.abs(wndd1-wndd2),0),[0.,5.,10.,350.,355.,360.]))
          nwindsame = nwindsame_msms + nwindsame_knot1 + nwindsame_knot2
          nwinddiff = nwind - nwindsame
        else:
          nwindsame = 0
        #if ((ntemp<=5) and (ntempsame>=3)) or \
        #   ((ntemp>=6) and (ntempsame>=ntemp-3)) or \
        #   ((nwind<=5) and (nwindsame>=3)) or \
        #   ((nwind>=6) and (nwindsame>=nwind-3)):
        #  # then we need to decide which one we keep
        if (ntemp>=3 and ntempsame>=3 and ntempdiff<0.5*ntemp) or \
           (ntemp<3 and nwind>=3 and nwindsame>=3 and nwinddiff<0.5*nwind):
          pref=-1
          if (nwind>0) and ((nwindsame_knot1>=3) or (nwindsame_knot2>=3)):
            if (nwindsame_knot1==nwindsame_knot2): # profiles have the same number of suspicious wind matches after knots to m/s conversion
              pref=-1
            else:
              if (nwindsame_knot1>nwindsame_knot2): # profile 1 has more suspicious winds than profile 2
                pref=2
                reason="knot"
              else:
                pref=1
                reason="knot"
          if pref==-1:
            profile_value_1 = get_RS_profile_value(df1,cols_compare[:-2],cols_compare[-2],cols_compare[-1])
          #print(df2[cols_compare])
          #print(ntemp,ntempsame,nwind,nwindsame)
            profile_value_2 = get_RS_profile_value(df2,cols_compare[:-2],cols_compare[-2],cols_compare[-1])
            #print(dfinn.iloc[wtemp][[var1,var2]])
            #print(dfsameas0.iloc[iindex]['sameas'])
            #print(iindex)
            # then we need to sort out and decide which one we prefer to keep
            # first check: which one has more humidity data?
          #nhum1 = np.count_nonzero(~np.isnan(dfinn['dewpointTemperature_1'].values[:]))
          #nhum2 = np.count_nonzero(~np.isnan(dfinn['dewpointTemperature_2'].values[:]))
          #nwin1 = np.count_nonzero((~np.isnan(dfinn['windSpeed_1'].values[:])) & (~np.isnan(dfinn['windDirection_1'].values[:])))
          #nwin2 = np.count_nonzero((~np.isnan(dfinn['windSpeed_2'].values[:])) & (~np.isnan(dfinn['windDirection_2'].values[:])))
          #nbco1 = np.count_nonzero(~np.isnan(dfinn['airTempBiasCorr_1'].values[:]))
          #nbco2 = np.count_nonzero(~np.isnan(dfinn['airTempBiasCorr_2'].values[:]))
          #lorp1 = np.count_nonzero(np.char.find(dfinn['wigos_statid_1'].values[:].astype(str),'orph')>=0)>0
          #lorp2 = np.count_nonzero(np.char.find(dfinn['wigos_statid_2'].values[:].astype(str),'orph')>=0)>0
          #pref = None
          #if nhum1<nhum2 or (nhum1==nhum2 and nwin1<nwin2) or (nhum1==nhum2 and nwin1==nwin2 and nbco1<nbco2) or (nhum1==nhum2 and nwin1==nwin2 and nbco1==nbco2 and lorp1):
          #print(profile_value_1,profile_value_2)
            if profile_value_2>profile_value_1:
              pref=2
              reason="value"
          #elif nhum1>nhum2 or (nhum1==nhum2 and nwin1>nwin2) or (nhum1==nhum2 and nwin1==nwin2 and nbco1>nbco2) or (nhum1==nhum2 and nwin1==nwin2 and nbco1==nbco2 and lorp2):
            elif profile_value_1>profile_value_2:
              pref=1
              reason="value"
            else:
          #  #print(discard_vals)
              #print(df1[cols_compare])
              #print(df2[cols_compare])
              pref=1 # RANDOM!!
              reason="random"
          #discard_vals[{1:iindex1,2:iindex2}[pref]] += 1
          #print("discarded",{1:n1,2:n2}[pref])
          #print(sameas_vals[iindex1])
          #print(sameas_vals[iindex2])
          #sameas_vals[iindex1]=np.char.add(sameas_vals[iindex1],[str(ijpair[1]) for x in range(n1)])
          #sameas_vals[iindex2]=np.char.add(sameas_vals[iindex2],[str(ijpair[0]) for x in range(n2)])
        else:
          pref=-1
        if pref in [1,2]:
          discarded_list0.append({1:statid2,2:statid1}[pref])
          df_dup0['wigos_statid'].append(discarded_list0[-1][0])
          df_dup0['report_id'].append   (discarded_list0[-1][1])
          #print(df1[cols_compare])
          #print(df2[cols_compare])
          #print(statid_pair,'won',discarded_list0[-1],reason,ntemp,ntempsame,nwind,nwindsame,lbothERA5_1)
      #print(sameas_vals)
      #print(discard_vals)
      #dfsameas0.iloc[iindex,0] = discard_vals

def function_sametrajectoryas(df0, index_center_call0, index_names0, dfFULL0, pairs_done0, discarded_list0, cross_stats0):

  bounds = {'pressureReducedToMeanSeaLevel':[85000.+10.,105000.-10.],'nonCoordinatePressure':[85000.+10.,105000.-10.],'oceanographicWaterTemperature':[268.15+0.5,309.15-0.5]}
  if 1==1:
    #print(df0.set_index(['statid']).index.unique())
    #stop
    df0_trajectories = df0.set_index(['statid']).index.unique()
    ntrajectories = len(df0_trajectories)
    if (ntrajectories>1):
      #print("number of trajectories",ntrajectories)
      # now pick one pair at a time and compute distance
      for ijpair in itertools.combinations(range(ntrajectories),2):
        statid1,statid2 = [x.strip() for x in sorted( ['{0:>12s}'.format(df0_trajectories[ijpair[0]]), '{0:>12s}'.format(df0_trajectories[ijpair[1]])] )] # note SORTING
        if (statid1 in discarded_list0) or (statid2 in discarded_list0):
          continue # Data were found to be inferior to another trajectory, so skip
        if (statid1,statid2) in pairs_done0:
          continue # Each pair of buoy only needs to be checked once
        pairs_done0.append((statid1,statid2))
        statid_pair = [statid1, statid2]
        dfs = []
        for ipos in [0,1]:
          iindex = np.where(dfFULL0['statid'].values[:]==statid_pair[ipos])[0]
          dfs.append(dfFULL0.iloc[iindex].copy())
        if len(dfs[0])<3 or len(dfs[1])<3:
          continue # we do not have much to deal with... too few points!
        # are the IDs a near-match?
        lsameID = (len(statid1)==5 and len(statid2) in [7,8] and statid1[:2]==statid2[:2] and statid2[2:4]=='00' and statid1[2:5]==statid2[4:7]) or \
                  (len(statid1)==7 and len(statid2)==8       and statid1[:7]==statid2[:7])
        # exact comparison for MSLP and SP
        lsameP = False
        cols_compare = ['datetime']
        ntimesame = {x:0 for x in ['p','sst']}
        dest={'nonCoordinatePressure':'p','pressureReducedToMeanSeaLevel':'p'}
        for var_choice in ['nonCoordinatePressure','pressureReducedToMeanSeaLevel']:
          dest1 = dest[var_choice]
          nok12 = [np.count_nonzero(~np.isnan(dfs[ipos][var_choice].values[:])) for ipos in [0,1]]
          if nok12==0 or nok12==0:
            continue
          df1 = dfs[0].set_index([var_choice])[cols_compare]
          df2 = dfs[1].set_index([var_choice])[cols_compare]
          if type(df1) is not pd.DataFrame or type(df2) is not pd.DataFrame:
            continue
          # we only rely on a comparison of the times for the exact, same observed values
          dfinn = pd.merge(df1,df2,how='inner',left_index=True,right_index=True,suffixes=('_1','_2'))
          if len(dfinn)==0:
            continue
          var0 = cols_compare[0]+'_{0}'; var1=var0.format(1); var2=var0.format(2)
          wtime = np.where((~np.isnan(dfinn.index.values[:])) & (dfinn.index.values[:]>bounds[var_choice][0]) & (dfinn.index.values[:]<bounds[var_choice][1]) & (~np.isnat(dfinn[var1].values[:])) & (~np.isnat(dfinn[var2].values[:])))[0]; ntime=len(wtime)
          if ntime>0:
            delta_time = dfinn.iloc[wtime][var1].values[:] - dfinn.iloc[wtime][var2].values[:]
            ntimesame[dest1] += np.count_nonzero(np.abs(delta_time.astype(int)/1e9)<600) # 10 minutes
        lsameP = (ntimesame['p']>100)
        # look at agreements between SSTs, tracks
        cols_compare=[['oceanographicWaterTemperature'],['lat','lon']]
        lsameSST = False
        lsameTRACK=False
        for cols_compare1 in cols_compare:
          woks = []
          for ipos in [0,1]:
            woks.append(np.where(eval(' & '.join(["(~np.isnan(dfs[{0}]['{1}']))".format(ipos,col1_compare1) for col1_compare1 in cols_compare1])))[0])
          if len(woks[0])>=1 and len(woks[1])>=1:
            df1 = dfs[0].iloc[woks[0]].set_index(['hour']).sort_index()[cols_compare1]
            df2 = dfs[1].iloc[woks[1]].set_index(['hour']).sort_index()[cols_compare1]
            dfinn = pd.merge(df1, df2, \
                             how='inner', left_index=True, right_index=True, suffixes=('_1','_2'))
            if len(dfinn)>0:

              if len(cols_compare1)==1: # SST comparison
                var0 = cols_compare1[0]+'_{0}'; var1=var0.format(1); var2=var0.format(2)
                SSTinput_medd  = [np.abs(np.median(dfinn[var1])-np.median(dfinn[var2]))]
                SSTinput_rstd  = [RobustStandardDeviation(dfinn[var1]-dfinn[var2])]
                SSTinput_pad95 = [PercentileAbsoluteDeviation(dfinn[var1],dfinn[var2],95.)]
                if len(dfinn)>=2:
                  SSTinput_medd  += [np.abs(np.median(dfinn[var1].values[1:  ])-np.median(dfinn[var2].values[ :-1])), \
                                     np.abs(np.median(dfinn[var1].values[ :-1])-np.median(dfinn[var2].values[1:  ]))]
                  SSTinput_rstd  += [RobustStandardDeviation(dfinn[var1].values[1:  ]-dfinn[var2].values[ :-1]), \
                                     RobustStandardDeviation(dfinn[var1].values[ :-1]-dfinn[var2].values[1:  ])]
                  SSTinput_pad95 += [PercentileAbsoluteDeviation(dfinn[var1].values[1:  ],dfinn[var2].values[ :-1],95.), \
                                     PercentileAbsoluteDeviation(dfinn[var1].values[ :-1],dfinn[var2].values[1:  ],95.)]
                SSTMEDD  = np.min(np.array( SSTinput_medd  ))
                SSTRSTD  = np.min(np.array( SSTinput_rstd  ))
                SSTPAD95 = np.min(np.array( SSTinput_pad95 ))
                SSTNI = len(dfinn); SSTN1 = len(woks[0]); SSTN2 = len(woks[1])
                lsameSST = ((np.around(SSTMEDD,2)<=0.05) and ((np.around(SSTRSTD,2)<=0.07) or (np.around(SSTPAD95,2)<=0.10))) or \
                           ((np.around(SSTMEDD,2)<=0.07) and ((np.around(SSTRSTD,2)<=0.07) or (np.around(SSTPAD95,2)<=0.07))) or \
                           ((np.around(SSTMEDD,2)<=0.10) and ((np.around(SSTRSTD,2)<=0.04) or (np.around(SSTPAD95,2)<=0.07)))
              else: # TRACK comparison
                az120,az210,dist_m0 = geoid.inv(dfinn['lon_1'].values[:], dfinn['lat_1'].values[:], \
                                                dfinn['lon_2'].values[:], dfinn['lat_2'].values[:])
                dist_km = dist_m0*1.e-3 # convert to km
                dist01km = np.percentile(dist_km, 1.)
                dist05km = np.percentile(dist_km, 5.)
                dist50km = np.percentile(dist_km,50.)
                TRACKNI = len(dfinn); TRACKN1 = len(woks[0]); TRACKN2 = len(woks[1])
                lsameTRACK=((dist05km<1.0) and (dist50km<5.0)) or ((dist05km<0.2) and (dist50km<7.0))
        if lsameTRACK:
          if lsameID or (len(statid1)==5 and len(statid2)==8 and statid2[:3]=='EXD'):
            pref=2
          else:
            cols_value = [['pressureReducedToMeanSeaLevel','nonCoordinatePressure'],'oceanographicWaterTemperature']
            trajectory_value_1 = get_DRIFTER_trajectory_value(dfs[0],cols_value,bounds)
            trajectory_value_2 = get_DRIFTER_trajectory_value(dfs[1],cols_value,bounds)
            if trajectory_value_2>trajectory_value_1:
              pref=2
            elif trajectory_value_1>trajectory_value_2:
              pref=1
            else:
              pref=2 # UNCLEAR -- but we need to pick one!
            #print(trajectory_value_1,trajectory_value_2,pref)
        else: # buoys differ
          pref = -1
        if pref in [1,2]:
          discarded_list0.append({1:statid2,2:statid1}[pref])
        cross_stats0['lsameI'].append(lsameID)
        cross_stats0['lsameP'].append(lsameP)
        cross_stats0['lsameT'].append(lsameSST)
        cross_stats0['lsameX'].append(lsameTRACK)
        cross_stats0['statid1'].append(statid1)
        cross_stats0['statid2'].append(statid2)
        cross_stats0['pref_statid'].append(pref)

def crawl_slice(df0,dims0,function_call0,index_center0=(),index_names0=[]):
  if len(df0)<2:
    return
  if len(dims0)>0:
    dim1 = dims0[0]
    xval = df0[dim1].values[:]
    vmin, vmax = np.min(xval), np.max(xval)
    check_lon = (dim1=='lon01i') and (vmin==-180) and (vmax==179) # then we need to handle this circular coordinate
    if check_lon: # perform roll-over left and right for min and max
      j1=vmin; wleft=np.where(xval==vmax)[0]; wright=np.where((xval>=(j1)) & (xval<=(j1+1)))[0]
      crawl_slice(pd.concat([ df0.iloc[wleft], df0.iloc[wright] ], axis=0), \
           dims0[1:], function_call0, index_center0=index_center0+(j1,), index_names0=index_names0+[dim1])
      j1=vmax; wleft=np.where((xval>=(j1-1)) & (xval<=(j1)))[0]; wright=np.where(xval==vmin)[0]
      crawl_slice(pd.concat([ df0.iloc[wleft], df0.iloc[wright] ], axis=0), \
           dims0[1:], function_call0, index_center0=index_center0+(j1,), index_names0=index_names0+[dim1])
    else: # perform min and max without roll-over
      j1=vmin; wmid=np.where((xval>=(j1-1)) & (xval<=(j1+1)))[0]
      crawl_slice(df0.iloc[wmid], \
           dims0[1:], function_call0, index_center0=index_center0+(j1,), index_names0=index_names0+[dim1])
      if vmax!=vmin:
        j1=vmax; wmid=np.where((xval>=(j1-1)) & (xval<=(j1+1)))[0]
        crawl_slice(df0.iloc[wmid], \
           dims0[1:], function_call0, index_center0=index_center0+(j1,), index_names0=index_names0+[dim1])
    # Process indices in between min and max
    xl=vmin+1; xu=vmax
    for j1 in range(xl,xu):
      #print(xl,j1,xu)
      wmid=np.where((xval>=(j1-1)) & (xval<=(j1+1)))[0]
      crawl_slice(df0.iloc[wmid], \
           dims0[1:], function_call0, index_center0=index_center0+(j1,), index_names0=index_names0+[dim1])
  else:
    # we are now in the cake!
    if len(index_names0)>1:
      index_center0_call = index_center0
    else:
      index_center0_call = index_center0[0]
    if len(df0.columns)==function_call0['nbcolumns']:
      function_call0['function'](df0, index_center0_call, index_names0, *function_call0['args'])

"""
# ISPD
fic = "/perm/erc/ERA6BUFR/ISPD/pickle/1980/07/pickle_19800701.pkl"; time_axis="datetime"
df = pd.read_pickle(fic).reset_index()
df['hour'] = ((df[time_axis]-np.min(df[time_axis])).astype(int)//1e9//3600).astype(int)
df['lat01i'] = np.around(df['lat']).astype(int)
df['lon01i'] = np.around(df['lon']).astype(int)
df = df.set_index(['hour','lat01i','lon01i']).sort_index()

entropy_dims = ['lat','lon',time_axis,'pressureReducedToMeanSeaLevel']
sigmas = [df[dim1].std() for dim1 in entropy_dims]+[1.]
#sigmas = [25., 90., pd.Timedelta('0 days 06:00:00'), obs_var]

wok=np.where((~np.isnan(df['pressureReducedToMeanSeaLevel'])) & (df.index.get_level_values('hour')>=11) & (df.index.get_level_values('hour')<=13))[0]
dfloc = df.iloc[wok].sort_index()
dfentropy = dfloc[['nonCoordinatePressure','pts']].rename(columns={'nonCoordinatePressure':'entropy','pts':'nobs'})
dfentropy['entropy'] = np.nan
dfentropy['entropy_gain'] = np.nan
for k in range(11):
  dfentropy['entropy_g{0:02d}'.format(k)] = np.nan
dfentropy['nobs']    = 0
function_call = {'nbcolumns':len(df.columns)+3, 'args':(entropy_dims,sigmas,dfentropy,), 'function':function_entropy}
dfsmall = crawl_slice(dfloc.reset_index(),dfloc.index.names,function_call)

went=np.where((dfentropy['nobs']>0) & (~np.isnan(dfentropy['entropy_gain'])))[0]
dfentloc = dfentropy.iloc[went]

fig1 = plt.figure()
xlats=dfentloc.index.get_level_values('lat01i')
xlons=dfentloc.index.get_level_values('lon01i')

ax1 = fig1.add_subplot(2,2,1,projection=ccrs.PlateCarree())
scat1 = ax1.scatter(xlons, xlats, c=dfentloc['nobs'], cmap=plt.get_cmap('RdYlBu_r'), s=0.2)
ax1.coastlines()
plt.colorbar(scat1,ax=ax1)

ax2 = fig1.add_subplot(2,2,2,projection=ccrs.PlateCarree())
scat2 = ax2.scatter(xlons, xlats, c=dfentloc['entropy'], cmap=plt.get_cmap('RdYlBu_r'), s=0.2)
ax2.coastlines()
plt.colorbar(scat2,ax=ax2)

ax3 = fig1.add_subplot(2,2,3,projection=ccrs.PlateCarree())
scat2 = ax3.scatter(xlons, xlats, c=dfentloc['entropy_gain'], cmap=plt.get_cmap('RdYlBu_r'), s=0.2)
ax3.coastlines()
plt.colorbar(scat2,ax=ax3)

ax4 = fig1.add_subplot(2,2,4)
ax4.hist2d(dfentloc['nobs'], dfentloc['entropy'], cmin=1)

fig2 = plt.figure()
mycmap = plt.get_cmap('RdYlBu_r')
nbounds = [1,2,5,10]
mynorm = mcolors.BoundaryNorm(nbounds, mycmap.N, extend='max')
ncols=6
for k in range(ncols):
  col_k = 'entropy_gain'; val_k_low = (k-1)*0.1; val_k_high = k*0.1
  w_k = np.where((dfentloc[col_k]>=val_k_low) & (dfentloc[col_k]<val_k_high))[0]
  ax_k = fig2.add_subplot(2,ncols//2,k+1, projection=ccrs.PlateCarree())
  scat_k = ax_k.scatter(xlons[w_k], xlats[w_k], color='red', s=0.2) # =dfentloc[col_k].values[w_k], cmap=mycmap, norm=mynorm, s=0.2)
  ax_k.coastlines()
  #plt.colorbar(scat_k,ax=ax_k)
  ax_k.set_title(('{0:4.1f}'+u'\u2264'+' '+u'\u0394'+'E_n <{1:4.1f}').format(val_k_low,val_k_high))
  ax_k.set_extent([-180,180,-90,90])
fig2.subplots_adjust(left=0.01,right=0.99,bottom=0.01,top=0.97,hspace=0.1,wspace=0.1)

#fig1.savefig("/home/erc/Work/k1ss_entropy.png")

# TESTING:
t_df = pd.DataFrame({'xaxis':[0. for x in range(10)],'zvalue':[-1,-1,-1,-1,-1,1,1,1,1,1],'phony_dim':[-1 for x in range(2)]+[0 for x in range(6)]+[1 for x in range(2)]}).set_index(['phony_dim']).sort_index()
t_entropy_dims = ['xaxis','zvalue']
t_sigmas = [t_df[dim1].std() for dim1 in t_entropy_dims]+[1.]
t_sigmas [0]= 1.
t_sigmas [1]= 1.

t_dfentropy = t_df.copy().drop(columns=['xaxis']).rename(columns={'zvalue':'entropy'})
t_dfentropy['entropy'] = np.nan
for k in range(11):
  t_dfentropy['entropy_g{0:02d}'.format(k)] = np.nan
t_dfentropy['nobs']    = 0
t_dfsmall = crawl_slice(t_df.reset_index(),t_df.index.names,t_entropy_dims,t_sigmas,t_dfentropy)
"""

"""
# CUON
fic = "/perm/erc/ERA6BUFR/CUON/pickle_19800601.pkl"; time_axis="time"
df = pd.read_pickle(fic).reset_index()
df['hour'] = ((df[time_axis]-np.min(df[time_axis])).astype(int)//1e9//3600).astype(int)
df['lat01i'] = np.around(df['lat']).astype(int)
df['lon01i'] = np.around(df['lon']).astype(int)
# Roll back longitudes>=180
wlon180 = np.where(df['lon01i']>=180.)[0]
if len(wlon180)>0:
  ilon = list(df.columns).index('lon01i')
  df.iloc[wlon180,ilon] -= 360
myindex_make = ['hour','lat01i','lon01i']
df = df.set_index(myindex_make).sort_index()

df_index_make = df.index.names+['station_id',time_axis]
#dfprofiles = pd.concat([ df.groupby(df_index_make).first()[['lat','lon']], df.groupby(df_index_make).count()[['pressure']].rename(columns={'pressure':'npoints'})], axis=1 )

#dfloc = df.iloc[ np.where((df.index.get_level_values('hour')>=11) & (df.index.get_level_values('hour')<=13) & (np.isin(df.index.get_level_values('lat01i'),[45,46,47])) & (np.isin(df.index.get_level_values('lon01i'),[172,173,174])))[0] ].reset_index().set_index(myindex_make)
#dfloc = df.iloc[ np.where((df.index.get_level_values('hour')>=3) & (df.index.get_level_values('hour')<=9) & (np.isin(df.index.get_level_values('lat01i'),[43,44,45])) & (np.isin(df.index.get_level_values('lon01i'),[28,29,30])))[0] ].reset_index().set_index(myindex_make)
#dfloc = df.iloc[ np.where((df.index.get_level_values('hour')>=3) & (df.index.get_level_values('hour')<=9) & (np.isin(df.index.get_level_values('lat01i'),[67])) & (np.isin(df.index.get_level_values('lon01i'),[27])))[0] ].reset_index().set_index(myindex_make)
dfloc = df.iloc[ np.where((df.index.get_level_values('hour')>=3) & (df.index.get_level_values('hour')<=9)) [0] ].reset_index().set_index(myindex_make).sort_index()
#dfloc = df.iloc[ np.where((df.index.get_level_values('hour')>=3) & (df.index.get_level_values('hour')<=9) & (np.isin(df.index.get_level_values('lat01i'),[66])) & (np.isin(df.index.get_level_values('lon01i'),[2])))[0] ].reset_index().set_index(myindex_make)
#dfprofilesloc = dfprofiles.iloc[  np.where((dfprofiles.index.get_level_values('hour')>=11) & (dfprofiles.index.get_level_values('hour')<=13))[0] ].reset_index().set_index(myindex_make)

dfsameas = dfloc[['wigos_statid']].copy()
dfsameas['discard'] = 0
dfsameas['nobs']   = 0
dfsameas.drop(columns=['wigos_statid'],inplace=True)
function_call = {'nbcolumns':len(dfloc.columns)+3, 'args':(dfsameas,), 'function':function_sameprofileas}
dfsmall = crawl_slice(dfloc.reset_index(),dfloc.index.names,function_call)
"""
