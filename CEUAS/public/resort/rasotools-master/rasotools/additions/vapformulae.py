
from numba import *
import numpy
import matplotlib.pyplot as plt
from scipy import stats,optimize
#import utils
import math

@njit
def muko(T):
   ew = 54.842763 - 6763.22 / T - 4.21*numpy.log(T) + 0.000367*T \
      + numpy.tanh(0.0415 *(T - 218.8))* (53.878 - 1331.22 / T - 9.44523*numpy.log(T) + 0.014025*T)
    

   return numpy.exp(ew)

@njit
def mukodiff(T,e):
   ew = 54.842763 - 6763.22 / T - 4.21*numpy.log(T) + 0.000367*T \
      + numpy.tanh(0.0415 *(T - 218.8))* (53.878 - 1331.22 / T - 9.44523*numpy.log(T) + 0.014025*T)
   
   return abs(e-ew)

def dpd(T,ew,rh):
   e=rh+ew
   
   Td=optimize.newton(mukodiff,273.0,args=(e,))
   
   return T-Td
@njit
def fdpd(T,ew,rh):   
   e=rh+ew
   x0=273.0
   maxiter=50
   tol=1.e-7
   # Secant method
   p0 = x0
   if x0 >= 0:
       p1 = x0*(1 + 1e-4) + 1e-4
   else:
       p1 = x0*(1 + 1e-4) - 1e-4
   q0 = mukodiff(p0,e) #func(*((p0,) + args))
   q1 = mukodiff(p1,e) #func(*((p1,) + args))
   for iter in range(maxiter):
       if q1 == q0:
           if p1 != p0:
               print("Tolerance reached") 
           return T-(p1 + p0)/2.0
       else:
           p = p1 - q1*(p1 - p0)/(q1 - q0)
       if abs(p - p1) < tol:
           return T-p
       p0 = p1
       q0 = q1
       p1 = p
       q1 = mukodiff(p1,e)
       

###############################################################
        
if __name__ == '__main__':

   T=numpy.arange(1000)/10.+233.15
   ew=numpy.log(muko(T))
   
   #plt.plot(T,muko(T))
   print(muko(273.15))
   #print dpd(T,0.4)
   #plt.show()
   rh=numpy.arange(96)+5
   rh=numpy.log(rh/100.)
   dpds=numpy.empty((T.shape[0],rh.shape[0]))
   n=10
   for irh in range(0,rh.shape[0]):
      print(irh)
      for it in range(0,T.shape[0]):
         dpds[it,irh]=fdpd(T[it],ew[it],rh[irh])
   
   plt.subplot(1,2,1)
   rh=numpy.exp(rh)
   plt.contourf(rh,T,dpds,levels=numpy.arange(400)/10.)
   plt.colorbar()
   plt.subplot(1,2,2)
   plt.plot(rh,dpds[600])
   plt.plot(rh,dpds[300])
   plt.plot(rh,dpds[200])
   plt.show()
   print('')
         
