import numpy as np

def mod(a, p):
     return a-int(a/p)*p

def snhteqsamp2y(test,ref,ni,istartorig,istoporig,maxlen,increment,miss_val,max_miss,critical_dates,ncritical,tsa,plus,minus,prms,mrms,pcount,mcount):
"""
snhtmov calculates running mean SNHT test parameter tsa
by dividing interval in two pieces and then calculating
means, rms of subsamples
means, rms are calculated recursively (= fairly efficiently)

Args:
    test=tested series
    ref=reference series
    ni size of test,ref
    istart= start index of analysis window
    istop=  end index of analysis window
    maxlen = length of averaging window
    miss_val = missing value indicator
    max_miss = maximum number of missing values in one of the two subsamples
    plus,minus = means of second, first subsample
    prms,mrms = rms of second, first subsample
    pcount,mcount = number of values used for calculating averages in second, first subsample
    
the SNHT test statistic tsa is returned
L. Haimberger, 2.7.2014
"""

#     integer         :: ni,maxlen,max_miss,istart,istop,jstart,jend,jstartlen,jendlen,i,k,l,jbin(12),psub,msub
#     integer,intent(in) :: istartorig,istoporig
#     real(kind=JPRM) :: miss_val

#     real(kind=JPRM),intent(in) :: test(ni),ref(ni)
#     integer,intent(in) :: bin(ni),ncritical,critical_dates(ncritical)
#     real(kind=JPRM) :: diff(ni),plus(ni),minus(ni),prms(ni),mrms(ni),mean(ni,12),square(ni,12)
#     real(kind=JPRM) :: tsa(ni),xm,xp,x,y,xy,sig,xx,yy,cratio,xpi,xmi

#     integer         :: pcount(ni),mcount(ni),gcount(ni,12),gindex(ni),ibin,pmon(12),mmon(12),stat
#     integer         :: increment,m2,j,ipc,imc,ipindex(maxlen/2,12),imindex(maxlen/2,12)

    istart=istartorig
    istop=istoporig
    istart=1
    istop=ni
    if(ni == 45000):
        istart=mod(20819,increment)+1

    if(istart < 1) or (istop > ni) :
        print('istart,istop ',istart,istop,' must be in range ',1,ni)
        return

    m2=maxlen/2
    tsa = []

    j=0
    mean = []
    square = []
    gcount = []
    for i in range(0, ni):
        if (test[i] != miss_val and ref[i] != miss_val):
            j=j+1
#             gindex(j)=i
            diff=test[i]-ref[i]
            if (j>1) :
                mean.append(mean[-1]+diff)
                square.append(square[-1]+diff*diff)
#                 mean(j)=mean(j-1)+diff
#                 square(j)=square(j-1)+diff*diff
            else:
                mean.append(diff)
                square.append(diff*diff)
#                 mean(j)=diff
#                 square(j)=diff*diff
        gcount.append(j)

    if(j < 2*(m2-miss_val)): 
        plus=miss_val
        minus=miss_val
        prms=miss_val
        mrms=miss_val
        pcount=0
        mcount=0
        return

    for k in range(m2-max_miss, ni-max_miss):
        xm=k-m2
        if (xm  <  1):
            xm=1

        xp=k+m2
        if (xp  >  ni):
            xp=ni

        pcount[k]=gcount[xp]-gcount[k]
        mcount[k]=gcount[k]-gcount[xm]

        if (gcount[k]-gcount[xm]  >  m2-max_miss) and  (count[xp]-gcount[k]  >  m2-max_miss):
            x=(mean[gcount[k]]-mean[gcount[xm]])/(gcount[k]-gcount[xm])
            y=(mean[gcount[xp]]-mean[gcount[k]])/(gcount[xp]-gcount[k])
            xy=(mean[gcount[xp]]-mean[gcount[xm]])/(gcount[xp]-gcount[xm])

            sig=(square[gcount[xp]]-square[gcount[xm]])/(gcount[xp]-gcount[xm])

            if (sig > 0):
                sig=np.sqrt(sig-xy*xy)
                tsa.append(((gcount[k]-gcount[xm])*(x-xy)*(x-xy)+(gcount[xp]-gcount[k])*(y-xy)*(y-xy))/sig) 
#                 tsa[k]=((gcount[k]-gcount[xm])*(x-xy)*(x-xy)+(gcount[xp]-gcount[k])*(y-xy)*(y-xy))/sig
            else:
                sig=0.0
                tsa.append(0.)
#                 tsa[k] = 0.
            plus[k]=y
            minus[k]=x
            prms[k]=(square[gcount[xp]]-square[gcount[k]])/(gcount[xp]-gcount[k])
            mrms[k]=(square[gcount[k]]-square[gcount[xm]])/(gcount[k]-gcount[xm])

    return
