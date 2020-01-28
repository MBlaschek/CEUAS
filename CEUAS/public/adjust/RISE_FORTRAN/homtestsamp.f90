module homtestsamp

contains

!! snhtmov calculates running mean SNHT test parameter tsa
!! by dividing interval in two pieces and then calculating
!! means, rms of subsamples
!! means, rms are calculated recursively (= fairly efficiently)
!!
!! test=tested series
!! ref=reference series
!! ni size of test,ref
!! istart= start index of analysis window
!! istop=  end index of analysis window
!! maxlen = length of averaging window
!! miss_val = missing value indicator
!! max_miss = maximum number of missing values in one of the two subsamples
!! plus,minus = means of second, first subsample
!! prms,mrms = rms of second, first subsample
!! pcount,mcount = number of values used for calculating averages in second, first subsample
!!
!! the SNHT test statistic tsa is returned
!!
!! L. Haimberger, 14.10.2004
!!
subroutine snhtsamp(test,ref,rcpara,istart,istop,tsa,plus,minus,prms,mrms,pcount,mcount)

use rfmod

implicit none

type(rasocor_namelist) :: rcpara

integer         :: istart,istop,jstart,jend,jstartlen,jendlen,k

real(kind=JPRM) :: test(rcpara%nmax),ref(rcpara%nmax),diff(rcpara%nmax),plus(rcpara%nmax),minus(rcpara%nmax),prms(rcpara%nmax),mrms(rcpara%nmax)
real(kind=JPRM) :: tsa(rcpara%nmax)
real(kind=JPRM) :: qquer,rms,sigq,z1,z2,ti,tti,tm,tmm,tp,tpp,eps
real(kind=JPRM) :: pmonsum(12),mmonsum(12),pmonsqsum(12),mmonsqsum(12)

integer         :: pcount(rcpara%nmax),mcount(rcpara%nmax)
integer         :: m2,j,ipc,imc,ipnc,imnc,ipoc,imoc,jmi,jm2i,ipcold,imcold
logical         :: ini,index(rcpara%nmax),ip(rcpara%snht_maxlen/2),im(rcpara%snht_maxlen/2),pmask(rcpara%snht_maxlen/2),mmask(rcpara%snht_maxlen/2)
logical         :: ymask(12),montouched(12)
integer         :: pmoncount(12),mmoncount(12),imon


if(istart .lt. 1 .or. istop .gt. rcpara%nmax) then
  print*,'istart,istop ',istart,istop,' must be in range ',1,rcpara%nmax
  return
endif

m2=rcpara%snht_maxlen/2
if(istart .eq. 1 .and. istop .eq. rcpara%nmax) then
tsa=rcpara%miss_val
else
!!tsa(istart+m2:istop-m2+1)=rcpara%miss_val
endif

diff=rcpara%miss_val

where(test .ne. rcpara%miss_val .and. ref .ne. rcpara%miss_val) diff=test-ref

index=diff .ne. rcpara%miss_val


if(istart .eq. 1 .and. istop .eq. rcpara%nmax) then
  plus=rcpara%miss_val ; minus=rcpara%miss_val ; prms=rcpara%miss_val ; mrms=rcpara%miss_val
  pcount=0 ; mcount=0
else
!!  plus(istart+m2:istop-m2+1)=rcpara%miss_val ; minus(istart+m2:istop-m2+1)=rcpara%miss_val ; prms(istart+m2:istop-m2+1)=rcpara%miss_val ; mrms(istart+m2:istop-m2+1)=rcpara%miss_val
!!  pcount(istart+m2:istop-m2+1)=0 ; mcount(istart+m2:istop-m2+1)=0
endif

imcold=0
ipcold=0
pmoncount=0
mmoncount=0
do j=istart-rcpara%max_miss,istop-rcpara%snht_maxlen+rcpara%max_miss+1,rcpara%snht_increment

  jstart=j
  if(jstart .lt. 1) jstart=1 
  jstartlen=j+m2-jstart
  jend=j+rcpara%snht_maxlen-1
  if(jend .gt. rcpara%nmax) jend=rcpara%nmax
  jendlen=m2+jend-(j+rcpara%snht_maxlen-1)


  ip(1:jstartlen)=diff(jstart:j+m2-1) .ne. rcpara%miss_val
  ipc=count(ip(1:jstartlen))
  imc=0.
  if(ipc .gt. 0) then 
    im(1:jendlen)=diff(j+m2:jend) .ne. rcpara%miss_val
    imc=count(im(1:jendlen))
  endif

  montouched=.true.
    where(pmoncount .gt. 3 .and. mmoncount .gt. 3) montouched=.false.
    do k=1,rcpara%snht_increment
      montouched(rcpara%month(j-1+k-rcpara%snht_increment))=.true.
      montouched(rcpara%month(j+m2-1+k-rcpara%snht_increment))=.true.
      montouched(rcpara%month(jend-1+k-rcpara%snht_increment))=.true.
    enddo
!!  montouched=.true.

  if(ipc*imc .gt. 0) then 

!! calculate monthly averages
   do imon=1,12
     if(montouched(imon)) then
       pmoncount(imon)=0
       mmoncount(imon)=0
       pmonsum(imon)=0.
       mmonsum(imon)=0.
       pmonsqsum(imon)=0.
       mmonsqsum(imon)=0.
       do k=1,jstartlen
         if(ip(k) .and. rcpara%month(jstart-1+k) .eq. imon) then
           pmonsum(imon)=pmonsum(imon)+diff(jstart-1+k)
           pmonsqsum(imon)=pmonsqsum(imon)+diff(jstart-1+k)*diff(jstart-1+k)
           pmoncount(imon)=pmoncount(imon)+1
         endif
       enddo
    
       do k=1,jendlen
         if(im(k) .and. rcpara%month(j+m2-1+k) .eq. imon) then
           mmonsum(imon)=mmonsum(imon)+diff(j+m2-1+k)
           mmonsqsum(imon)=mmonsqsum(imon)+diff(j+m2-1+k)*diff(j+m2-1+k)
           mmoncount(imon)=mmoncount(imon)+1
         endif
       enddo
    
!!    pmask(1:jstartlen)=ip(1:jstartlen) .and. rcpara%month(jstart:j+m2-1) .eq. imon
!!    pmonsum(imon)=sum(diff(jstart:j+m2-1),pmask(1:jstartlen))
!!    pmoncount(imon)=count(pmask(1:jstartlen))

!!    mmask(1:jendlen)=im(1:jendlen) .and. rcpara%month(j+m2:jend) .eq. imon
!!    mmonsum(imon)=sum(diff(j+m2:jend),mmask(1:jendlen))
!!    mmoncount(imon)=count(mmask(1:jendlen))
     endif
   enddo
  endif

  ymask=pmoncount .gt. 3 .and. mmoncount .gt. 3
  ipc=count(ymask)
  if(ipc .gt. 3) then
    imc=ipc
    tp=sum(pmonsum/pmoncount,ymask)
    tm=sum(mmonsum/mmoncount,ymask)
    tpp=sum(pmonsum*pmonsum/pmoncount/pmoncount,ymask)
    tmm=sum(mmonsum*mmonsum/mmoncount/mmoncount,ymask)


    qquer=(tp+tm)/(ipc+imc)
    rms=(tpp+tmm)/(ipc+imc)
    sigq=sqrt(rms-qquer*qquer)  
       
    if(sigq .gt. 0.001) then
      z1=(tp/ipc-qquer)/sigq
      z2=(tm/imc-qquer)/sigq
      tsa(j+m2:j+m2+rcpara%snht_increment-1)=ipc*z1*z1+imc*z2*z2
      plus(j+m2:j+m2+rcpara%snht_increment-1)=tp/ipc
      minus(j+m2:j+m2+rcpara%snht_increment-1)=tm/imc
      prms(j+m2:j+m2+rcpara%snht_increment-1)=tpp/ipc
      mrms(j+m2:j+m2+rcpara%snht_increment-1)=tmm/imc
    endif
    pcount(j+m2:j+m2+rcpara%snht_increment-1)=ipc
    mcount(j+m2:j+m2+rcpara%snht_increment-1)=imc
  endif
enddo


return
end subroutine snhtsamp

end module homtestsamp
