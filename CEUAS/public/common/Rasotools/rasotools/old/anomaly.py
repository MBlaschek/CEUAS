
def manomaly( series,startyear,interval,anomaly,climatology):
#   import numpy
   from scipy.stats.stats import nanmean
   import numexpr as ne
   
   n=len(series)
   for j in range(0,12):
      hanomaly=series[j:n:12]
      
#      clnan=hanomaly[interval[0]-startyear:interval[1]-startyear]
#      clmask=numpy.ma.masked_array(clnan,numpy.isnan(clnan))
#      climatology[j]=numpy.mean(clmask)
      climatology[j]=nanmean(hanomaly[interval[0]-startyear:interval[1]-startyear])
      anomaly[j:n:12]=hanomaly-climatology[j]
   return

