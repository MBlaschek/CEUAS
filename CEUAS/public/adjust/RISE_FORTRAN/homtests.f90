module homtests

use rfmod

interface
  subroutine compute_segmentation_mean(diff,ndiff,kmax,lmin,lmax,vh,j_est,t_est)
    real*8,intent(in) :: diff(ndiff)
    integer*4 :: ndiff,kmax,lmin,lmax
    logical vh
    real*8 :: j_est(kmax)
    real*8:: t_est(kmax,kmax)
  end subroutine compute_segmentation_mean
end interface

contains



! snhtmov calculates running mean SNHT test parameter tsa
! by dividing interval in two pieces and then calculating
! means, rms of subsamples
! means, rms are calculated recursively (= fairly efficiently)
!
! test=tested series
! ref=reference series
! ni size of test,ref
! istart= start index of analysis window
! istop=  end index of analysis window
! maxlen = length of averaging window
! miss_val = missing value indicator
! max_miss = maximum number of missing values in one of the two subsamples
! plus,minus = means of second, first subsample
! prms,mrms = rms of second, first subsample
! pcount,mcount = number of values used for calculating averages in second, first subsample
!
! the SNHT test statistic tsa is returned
!
! L. Haimberger, 7.7.2004
!
subroutine snhtmov(test,ref,ni,istart,istop,maxlen,increment,miss_val,max_miss,tsa,plus,minus,prms,mrms,pcount,mcount)

implicit none

integer         :: ni,maxlen,max_miss,istart,istop,jstart,jend,jstartlen,jendlen,i
real(kind=JPRM) :: miss_val,inn

real(kind=JPRM) :: test(ni),ref(ni),diff(ni),plus(ni),minus(ni),prms(ni),mrms(ni)
real(kind=JPRM) :: tsa(ni),innov(increment),plusshort(maxlen/2),minusshort(maxlen/2)
real(kind=JPRM) :: qquer,rms,sigq,z1,z2,ti,tti,tm,tmm,tp,tpp,eps

integer         :: pcount(ni),mcount(ni)
integer         :: increment,m2,j,ipc,imc,ipnc,imnc,ipoc,imoc,jmi,jm2i
logical         :: ini,index(ni),ip(maxlen/2),im(maxlen/2),ipn(increment),imn(increment),ipo(increment),imo(increment)


if(istart .lt. 1 .or. istop .gt. ni) then
  print*,'istart,istop ',istart,istop,' must be in range ',1,ni
  return
endif

m2=maxlen/2
if(istart .eq. 1 .and. istop .eq. ni) then
tsa=miss_val
else
!tsa(istart+m2:istop-m2+1)=miss_val
endif

diff=miss_val

where(test .ne. miss_val .and. ref .ne. miss_val) diff=test-ref

index=diff .ne. miss_val

!if(count(index) .lt. maxlen) return


if(istart .eq. 1 .and. istop .eq. ni) then
  plus=miss_val ; minus=miss_val ; prms=miss_val ; mrms=miss_val
  pcount=0 ; mcount=0
else
!  plus(istart+m2:istop-m2+1)=miss_val ; minus(istart+m2:istop-m2+1)=miss_val ; prms(istart+m2:istop-m2+1)=miss_val ; mrms(istart+m2:istop-m2+1)=miss_val
!  pcount(istart+m2:istop-m2+1)=0 ; mcount(istart+m2:istop-m2+1)=0
endif

ini=.true.

do j=istart-max_miss,istop-maxlen+max_miss+1,increment

  jstart=j
  if(jstart .lt. 1) jstart=1
  jstartlen=j+m2-jstart
  jend=j+maxlen-1
  if(jend .gt. ni) jend=ni
  jendlen=m2+jend-(j+maxlen-1)

  if(jstart .eq. 1 .or. jend .eq. ni .or. ini) then
    ip(1:jstartlen)=diff(jstart:j+m2-1) .ne. miss_val
    ipc=count(ip(1:jstartlen))
    im(1:jendlen)=diff(j+m2:jend) .ne. miss_val
    imc=count(im(1:jendlen))
  else
    jm2i=j+m2-increment
    ipn=diff(jm2i:j+m2-1) .ne. miss_val
    ipnc=count(ipn)
    
    jmi=j+maxlen-increment
    imn=diff(jmi:jend) .ne. miss_val
    imnc=count(imn)

    ipo=diff(jstart-increment:jstart-1) .ne. miss_val
    ipoc=count(ipo)

    imo=ipn
    imoc=ipnc

    imc=imc+imnc-imoc
    ipc=ipc+ipnc-ipoc
  endif

  if(maxlen/2-ipc .lt. max_miss .and. maxlen/2-imc .lt. max_miss) then
       
       if(jstart .eq. 1 .or. jend .eq. ni .or. ini) then
         where(ip(1:jstartlen)) plusshort(1:jstartlen)=diff(jstart:j+m2-1)
         where(im(1:jendlen)) minusshort(1:jendlen)=diff(j+m2:jend)
         tp=sum(plusshort(1:jstartlen),ip(1:jstartlen))
         tm=sum(minusshort(1:jendlen),im(1:jendlen))
         tpp=sum(plusshort(1:jstartlen)*plusshort(1:jstartlen),ip(1:jstartlen))
         tmm=sum(minusshort(1:jendlen)*minusshort(1:jendlen),im(1:jendlen))
         ini=.false.
       else
         if(ipnc .gt. 0) then
           where(ipn) innov=diff(jm2i:j+m2-1)
           ti=sum(innov,ipn)
           tti=sum(innov*innov,ipn)
!where, sum construct is faster          
!           ti=0.
!           tti=0.
!           do i=1,increment
!             if(ipn(i)) then
!               inn=diff(jm2i+i-1)
!               ti=ti+inn
!               tti=tti+inn*inn
!             endif
!           enddo
           tp=tp+ti
           tpp=tpp+tti
           tm=tm-ti
           tmm=tmm-tti
         endif
         if(ipoc .gt. 0) then
           where(ipo) innov=diff(jstart-increment:jstart-1)
           tp=tp-sum(innov,ipo)
           tpp=tpp-sum(innov*innov,ipo)
!           do i=1,increment
!             if(ipo(i)) then
!               inn=diff(jstart-increment+i-1)
!               tp=tp-inn
!               tpp=tpp-inn*inn
!             endif
!           enddo
         endif
         if(imnc .gt. 0) then
           where(imn) innov=diff(jmi:jend)
           tm=tm+sum(innov,imn)
           tmm=tmm+sum(innov*innov,imn)
!           do i=1,increment
!             if(imn(i)) then
!               inn=diff(jmi+i-1)
!               tm=tm+inn
!               tmm=tmm+inn*inn
!             endif
!           enddo
         endif
         eps=0.001
       endif

       qquer=(tp+tm)/(ipc+imc)
       rms=(tpp+tmm)/(ipc+imc)
       sigq=sqrt(rms-qquer*qquer)  
       
     if(sigq .gt. 0.001) then
       z1=(tp/ipc-qquer)/sigq
       z2=(tm/imc-qquer)/sigq
       tsa(j+m2:j+m2+increment-1)=ipc*z1*z1+imc*z2*z2
       plus(j+m2:j+m2+increment-1)=tp/ipc
       minus(j+m2:j+m2+increment-1)=tm/imc
       prms(j+m2:j+m2+increment-1)=tpp/ipc
       mrms(j+m2:j+m2+increment-1)=tmm/imc
     endif
 else
   ini=.true.
 endif
   pcount(j+m2:j+m2+increment-1)=ipc
   mcount(j+m2:j+m2+increment-1)=imc

enddo


!if(any(plus .ne. miss_val .and. abs(plus) .gt. 10.)) then
!  write(*,*) 'snhtmov: plus invalid'
!  call exit(1)
!endif
!if(any(minus .ne.miss_val .and. abs(minus) .gt. 10.)) then
!  write(*,*) 'snhtmov: minus invalid'
!  call exit(1)
!endif

return
end subroutine snhtmov

! snhtmov calculates running mean SNHT test parameter tsa
! by dividing interval in two pieces and then calculating
! means, rms of subsamples
! means, rms are calculated recursively (= fairly efficiently)
!
! test=tested series
! ref=reference series
! ni size of test,ref
! istart= start index of analysis window
! istop=  end index of analysis window
! maxlen = length of averaging window
! miss_val = missing value indicator
! max_miss = maximum number of missing values in one of the two subsamples
! plus,minus = means of second, first subsample
! prms,mrms = rms of second, first subsample
! pcount,mcount = number of values used for calculating averages in second, first subsample
!
! the SNHT test statistic tsa is returned
!
! L. Haimberger, 7.7.2004
!
subroutine snhteqsamp(test,ref,ni,istart,istop,maxlen,increment,miss_val,max_miss,tsa,plus,minus,prms,mrms,pcount,mcount)

implicit none

integer         :: ni,maxlen,max_miss,istart,istop,jstart,jend,jstartlen,jendlen,i,k,stat
real(kind=JPRM) :: miss_val

real(kind=JPRM) :: test(ni),ref(ni),diff(ni),plus(ni),minus(ni),prms(ni),mrms(ni)
real(kind=JPRM) :: tsa(ni),innov(increment),plusshort(maxlen/2),minusshort(maxlen/2)
real(kind=JPRM) :: qquer,rms,sigq,z1,z2,ti,tti,tm,tmm,tp,tpp,eps

integer         :: pcount(ni),mcount(ni),ibin,pmon(12),mmon(12)
integer         :: increment,m2,j,ipc,imc,ipnc,imnc,ipoc,imoc,jmi,jm2i,ipcall,imcall
logical         :: ini,ip(maxlen/2),im(maxlen/2),iph(maxlen/2),imh(maxlen/2),ipn(increment),imn(increment),ipo(increment),imo(increment)
real(kind=JPRM),save,allocatable :: bin(:)


if(.not. allocated(bin)) then
  allocate(bin(ni),stat=stat)
  do i=1,ni
    bin(i)=modulo(i*1.,365.)
    bin(i)=floor(bin(i)/365.*12+1)
  enddo
endif

istart=1
istop=ni
if(istart .lt. 1 .or. istop .gt. ni) then
  print*,'istart,istop ',istart,istop,' must be in range ',1,ni
  return
endif

m2=maxlen/2
if(istart .eq. 1 .and. istop .eq. ni) then
tsa=miss_val
else
!tsa(istart+m2:istop-m2+1)=miss_val
endif

diff=miss_val

where(test .ne. miss_val .and. ref .ne. miss_val) diff=test-ref
do while(diff(istart) .eq. miss_val .and. istart .lt. istop)
  istart=istart+1
enddo
do while(diff(istop) .eq. miss_val .and. istop .gt. istart)
  istop=istop+1
enddo

if(istart .eq. ni) return

if(istart .eq. 1 .and. istop .eq. ni) then
  plus=miss_val ; minus=miss_val ; prms=miss_val ; mrms=miss_val
  pcount=0 ; mcount=0
else
  plus=miss_val ; minus=miss_val ; prms=miss_val ; mrms=miss_val
  pcount=0 ; mcount=0
!  plus(istart+m2:istop-m2+1)=miss_val ; minus(istart+m2:istop-m2+1)=miss_val ; prms(istart+m2:istop-m2+1)=miss_val ; mrms(istart+m2:istop-m2+1)=miss_val
!  pcount(istart+m2:istop-m2+1)=0 ; mcount(istart+m2:istop-m2+1)=0
endif


ini=.true.

do j=istart-max_miss,istop-maxlen+max_miss+1,increment

  jstart=j
  if(jstart .lt. 1) jstart=1
  jstartlen=j+m2-jstart
  jend=j+maxlen-1
  if(jend .gt. ni) jend=ni
  jendlen=m2+jend-(j+maxlen-1)

  if(jstart .eq. 1 .or. jend .eq. ni .or. ini) then
    ip=.false.
    ipc=0
!    do i=1,jstartlen
!      if(j+m2+i-1 .le. ni) then
!        ip(i)=diff(jstart+i-1) .ne. miss_val .and. diff(j+m2+i-1) .ne. miss_val
!        ipc=ipc+1
!      endif
!    enddo
    im=.false.
    imc=0
!    do i=1,jendlen
!      if(j+i-1 .gt. 0) then
!        im(i)=diff(j+m2+i-1) .ne. miss_val .and. diff(j+i-1) .ne. miss_val
!        imc=imc+1
!      endif
!    enddo

    ip(1:jstartlen)=diff(jstart:j+m2-1) .ne. miss_val
    ipc=count(ip(1:jstartlen))
    im(1:jendlen)=diff(j+m2:jend) .ne. miss_val
    imc=count(im(1:jendlen))
    ipcall=ipc
    imcall=imc

  else
  endif

  if(maxlen/2-ipc .lt. max_miss .and. maxlen/2-imc .lt. max_miss) then
       
    do ibin=1,12
      iph(1:jstartlen)=ip(1:jstartlen) .and. bin(jstart:j+m2-1) .eq. ibin
      pmon(ibin)=count(iph(1:jstartlen))
      imh(1:jendlen)=im(1:jendlen) .and. bin(j+m2:jend) .eq. ibin
      mmon(ibin)=count(imh(1:jendlen))

      if(pmon(ibin) .lt. mmon(ibin)) then
        i=0
        k=0
        do while(i .lt. mmon(ibin)-pmon(ibin))
          do while( .not. imh(jendlen-k) .and. k .lt. jendlen-1) 
            k=k+1
          enddo
          imh(jendlen-k)=.false.
          im(jendlen-k)=.false.
          i=i+1
        enddo
      endif
      if(mmon(ibin) .lt. pmon(ibin)) then
        i=0
        k=0
        do while(i .lt. pmon(ibin)-mmon(ibin))
          do while( .not. iph(1+k) .and. k .lt. jstartlen-1) 
             k=k+1
          enddo
          iph(1+k)=.false.
          ip(1+k)=.false.
          i=i+1
        enddo
      endif
    enddo 
    ipc=count(ip(1:jstartlen))
    imc=count(im(1:jendlen))
!    do ibin=1,12
!      pmon(ibin)=count(ip(1:jstartlen) .and. bin(jstart:j+m2-1) .eq. ibin)
!      mmon(ibin)=count(im(1:jendlen) .and. bin(j+m2:jend) .eq. ibin)
!    enddo
       if(jstart .eq. 1 .or. jend .eq. ni .or. ini) then
         where(ip(1:jstartlen)) plusshort(1:jstartlen)=diff(jstart:j+m2-1)
         where(im(1:jendlen)) minusshort(1:jendlen)=diff(j+m2:jend)
         tp=sum(plusshort(1:jstartlen),ip(1:jstartlen))
         tm=sum(minusshort(1:jendlen),im(1:jendlen))
!         write(12,'(I6,2F5.2)') j,tp/ipc,tm/imc
         tpp=sum(plusshort(1:jstartlen)*plusshort(1:jstartlen),ip(1:jstartlen))
         tmm=sum(minusshort(1:jendlen)*minusshort(1:jendlen),im(1:jendlen))
         ini=.false.
    if(j .eq. 31) then
!      write(*,*) j
    endif
       else
       endif

       qquer=(tp+tm)/(ipc+imc)
       rms=(tpp+tmm)/(ipc+imc)
       sigq=sqrt(rms-qquer*qquer)  
       
     if(sigq .gt. 0.001) then
       z1=(tp/ipc-qquer)/sigq
       z2=(tm/imc-qquer)/sigq
       tsa(j+m2:j+m2+increment-1)=ipc*z1*z1+imc*z2*z2
       plus(j+m2:j+m2+increment-1)=tp/ipc
       minus(j+m2:j+m2+increment-1)=tm/imc
       prms(j+m2:j+m2+increment-1)=tpp/ipc
       mrms(j+m2:j+m2+increment-1)=tmm/imc
     endif
 else
   ini=.true.
 endif
   ini=.true.
   pcount(j+m2:j+m2+increment-1)=ipcall
   mcount(j+m2:j+m2+increment-1)=imcall

enddo

!where(diff .eq. miss_val) tsa=miss_val


return
end subroutine snhteqsamp

! snhtmov calculates running mean SNHT test parameter tsa
! by dividing interval in two pieces and then calculating
! means, rms of subsamples
! means, rms are calculated recursively (= fairly efficiently)
!
! test=tested series
! ref=reference series
! ni size of test,ref
! istart= start index of analysis window
! istop=  end index of analysis window
! maxlen = length of averaging window
! miss_val = missing value indicator
! max_miss = maximum number of missing values in one of the two subsamples
! plus,minus = means of second, first subsample
! prms,mrms = rms of second, first subsample
! pcount,mcount = number of values used for calculating averages in second, first subsample
!
! the SNHT test statistic tsa is returned
!
! L. Haimberger, 7.7.2004
!
subroutine snhteqsamp2(test,ref,ni,istartorig,istoporig,maxlen,increment,miss_val,max_miss,critical_dates,ncritical,tsa,plus,minus,prms,mrms,pcount,mcount,bin) !

implicit none

integer         :: ni,maxlen,max_miss,istart,istop,jstart,jend,jstartlen,jendlen,i,k,l
integer,intent(in) :: istartorig,istoporig
real(kind=JPRM) :: miss_val

real(kind=JPRM),intent(in) :: test(ni),ref(ni)
integer,intent(in) :: bin(ni)
real(kind=JPRM) :: diff(ni),plus(ni),minus(ni),prms(ni),mrms(ni)
real(kind=JPRM) :: tsa(ni),innov(increment),plusshort(maxlen/2),minusshort(maxlen/2)
real(kind=JPRM) :: qquer,rms,sigq,z1,z2,ti,tti,tm,tmm,tp,tpp,eps

integer         :: pcount(ni),mcount(ni),ibin,pmon(12),mmon(12),ncritical,critical_dates(ncritical),stat
integer         :: increment,m2,j,ipc,imc,ipnc,imnc,ipoc,imoc,jmi,jm2i,ipcall,imcall,ipindex(maxlen/2,12),imindex(maxlen/2,12)
logical         :: ini,ip(maxlen/2),im(maxlen/2),ipn(increment),imn(increment),ipo(increment),imo(increment)

!!$real(kind=JPRM),save,allocatable :: bin(:)
!!$
!!$
!!$if( .not. allocated(bin)) then
!!$  allocate(bin(ni),stat=stat)
!!$  do i=1,ni
!!$    bin(i)=modulo(i*1.,365.)
!!$    bin(i)=floor(bin(i)/365.*12+1)
!!$  enddo
!!$endif

istart=istartorig
istop=istoporig
istart=1
istop=ni
if(ni .eq. 45000) istart=mod(20819,increment)+1

if(istart .lt. 1 .or. istop .gt. ni) then
  print*,'istart,istop ',istart,istop,' must be in range ',1,ni
  return
endif

m2=maxlen/2
!if(istart .eq. 1 .and. istop .eq. ni) then
tsa=miss_val
!else
!tsa(istart+m2:istop-m2+1)=miss_val
!endif

i=1
do while(i .lt. ni .and. (test(i) .eq. miss_val .or. ref(i) .eq. miss_val))
 diff(i)=miss_val
 i=i+1
enddo
istart=i
do i=istart,ni
  if(test(i) .eq. miss_val .or. ref(i) .eq. miss_val) then
    diff(i)=miss_val
  else
    diff(i)=test(i)-ref(i)
    istop=i
  endif
enddo
!!$where(test .ne. miss_val .and. ref .ne. miss_val) diff=test-ref
!!$do while(diff(istart) .eq. miss_val .and. istart .lt. istop)
!!$  istart=istart+1
!!$enddo
!!$do while(diff(istop) .eq. miss_val .and. istop .gt. istart)
!!$  istop=istop-1
!!$enddo

if(istart .eq. ni) return

!if(istart .eq. 1 .and. istop .eq. ni) then
!  plus=miss_val ; minus=miss_val ; prms=miss_val ; mrms=miss_val
!  pcount=0 ; mcount=0
!else
  plus=miss_val ; minus=miss_val ; prms=miss_val ; mrms=miss_val
  pcount=0 ; mcount=0
!endif


ini=.true.
do j=istart-max_miss,istop-maxlen+max_miss+1,increment

  jstart=j
  if(jstart .lt. 1) jstart=1
  jstartlen=j+m2-jstart
  jend=j+maxlen-1
  if(jend .gt. ni) jend=ni
  jendlen=m2+jend-(j+maxlen-1)

  if(jstart .eq. 1 .or. jend .eq. ni .or. ini) then
    ip=.false.
    ipc=0
    im=.false.
    imc=0

    ip(1:jstartlen)=diff(jstart:j+m2-1) .ne. miss_val
    do l=1,ncritical
      if(jstart .le. critical_dates(l) .and. j+m2-1 .gt. critical_dates(l)) then
        ip(1:critical_dates(l)-jstart+1)=.false.
      endif
    enddo 
    ipc=count(ip(1:jstartlen))
    im(1:jendlen)=diff(j+m2:jend) .ne. miss_val
    do l=1,ncritical
      if(j+m2 .le. critical_dates(l) .and. jend .gt. critical_dates(l)) then
        im(jendlen-(jend-critical_dates(l))+1:jendlen)=.false.
      endif
    enddo 
    imc=count(im(1:jendlen))
    ipcall=ipc
    imcall=imc

  else
  endif

  if(maxlen/2-ipc .lt. max_miss .and. maxlen/2-imc .lt. max_miss) then
       
!    do ibin=1,12
!      l=0
      pmon=0
      mmon=0
      do i=1,jstartlen
        if(ip(i)) then ! .and. jstart+i-1 .lt. ni.and. jstart+i-1 .gt. 0 ) then
          ibin=bin(jstart+i-1)
          pmon(ibin)=pmon(ibin)+1
          ipindex(pmon(ibin),ibin)=i
        endif         
      enddo
     
      do i=1,jendlen
        if(im(i)) then ! .and. j+m2+i-1 .le. ni .and. j+m2+i-1 .gt. 0) then
          ibin=bin(j+m2+i-1)
          mmon(ibin)=mmon(ibin)+1
          imindex(mmon(ibin),ibin)=i
        endif       
      enddo

      
    do ibin=1,12
      l=mmon(ibin)-pmon(ibin)
      if(l .gt. 0) then 
        do i=1,l
          im(imindex(mmon(ibin),ibin))=.false.
          mmon(ibin)=mmon(ibin)-1
          imc=imc-1
        enddo
      endif
      
      l=pmon(ibin)-mmon(ibin)
      if(l .gt. 0) then 
        do i=1,l
          ip(ipindex(i,ibin))=.false.
          pmon(ibin)=pmon(ibin)-1
          ipc=ipc-1
        enddo
      endif

    enddo 

!    write(*,'(I6,24I4)') j+m2,mmon,pmon
!    ipc=sum(pmon) !count(ip(1:jstartlen))
!    imc=sum(mmon) !count(im(1:jendlen))
         tp=0.
         tpp=0.
         tm=0.
         tmm=0.
       if(jstart .eq. 1 .or. jend .eq. ni .or. ini) then

         do i=1,jstartlen
            if(ip(i)) then
              tp=tp+diff(jstart+i-1)
              tpp=tpp+diff(jstart+i-1)*diff(jstart+i-1)
            endif
         enddo
         do i=1,jendlen
            if(im(i)) then
              tm=tm+diff(j+m2+i-1)
              tmm=tmm+diff(j+m2+i-1)*diff(j+m2+i-1)
            endif
         enddo
         ini=.false.
       else
       endif

       if(ipc .gt. 0 .and. imc .gt. 0) then
         qquer=(tp+tm)/(ipc+imc)
         rms=(tpp+tmm)/(ipc+imc)
         if(rms .gt. qquer*qquer) then
           sigq=sqrt(rms-qquer*qquer)  
         else
           sigq=0
         endif
       else
         sigq=0
       endif
       
     if(sigq .gt. 0.001) then
       z1=(tp/ipc-qquer)/sigq
       z2=(tm/imc-qquer)/sigq
       tsa(j+m2:j+m2+increment-1)=ipc*z1*z1+imc*z2*z2
       plus(j+m2:j+m2+increment-1)=tp/ipc
       minus(j+m2:j+m2+increment-1)=tm/imc
       prms(j+m2:j+m2+increment-1)=tpp/ipc
       mrms(j+m2:j+m2+increment-1)=tmm/imc
     endif
 else
   ini=.true.
 endif
   ini=.true.
   pcount(j+m2:j+m2+increment-1)=ipc !all
   mcount(j+m2:j+m2+increment-1)=imc !all

enddo

do l=1,ncritical
  tsa(critical_dates(l)-maxlen/2+max_miss:critical_dates(l)+maxlen/2-max_miss)=miss_val
  plus(critical_dates(l)-maxlen/2+max_miss:critical_dates(l)+maxlen/2-max_miss)=miss_val
  minus(critical_dates(l)-maxlen/2+max_miss:critical_dates(l)+maxlen/2-max_miss)=miss_val
  prms(critical_dates(l)-maxlen/2+max_miss:critical_dates(l)+maxlen/2-max_miss)=miss_val
  mrms(critical_dates(l)-maxlen/2+max_miss:critical_dates(l)+maxlen/2-max_miss)=miss_val
enddo
!where(diff .eq. miss_val) tsa=miss_val


return
end subroutine snhteqsamp2


! meaneqsamp calculates running mean SNHT test parameter tsa
! by dividing interval in two pieces and then calculating
! means, rms of subsamples
! means, rms are calculated recursively (= fairly efficiently)
!
! test=tested series
! ref=reference series
! ni size of test,ref
! istart= start index of analysis window
! istop=  end index of analysis window
! maxlen = length of averaging window
! miss_val = missing value indicator
! max_miss = maximum number of missing values in one of the two subsamples
! plus,minus = means of second, first subsample
! prms,mrms = rms of second, first subsample
! pcount,mcount = number of values used for calculating averages in second, first subsample
!
! the SNHT test statistic tsa is returned
!
! L. Haimberger, 7.7.2004
!
subroutine meaneqsamp(test,ref,ni,istart,left_maxlen,right_maxlen,increment,miss_val,max_miss,tsa,plus,minus,prms,mrms,pcount,mcount)

implicit none

integer         :: ni,left_maxlen,right_maxlen,max_miss,istart,jstart,jend,jstartlen,jendlen,i,k,minlen
real(kind=JPRM) :: miss_val

real(kind=JPRM),intent(in) :: test(ni),ref(ni)
real(kind=JPRM) :: diff(ni),plus,minus,prms,mrms
real(kind=JPRM) :: tsa,innov(increment),plusshort(left_maxlen),minusshort(right_maxlen)
real(kind=JPRM) :: qquer,rms,sigq,z1,z2,ti,tti,tm,tmm,tp,tpp,eps

integer         :: pcount,mcount,ibin,pmon(12),mmon(12)
integer         :: increment,j,ipc,imc,ipnc,imnc,ipoc,imoc,jmi,jm2i,mmonhilf,pmonhilf,stat
logical         :: ini,ip(left_maxlen),im(right_maxlen),iph(left_maxlen),imh(right_maxlen),ipn(increment),imn(increment),ipo(increment),imo(increment)
real(kind=JPRM),save,allocatable :: bin(:)

if( .not. allocated(bin)) then
  allocate(bin(ni),stat=stat)   
  do i=1,ni
    bin(i)=modulo(i*1.,365.)
    bin(i)=floor(bin(i)/365.*12+1)
  enddo
endif

if(istart .lt. 1-left_maxlen .or. istart .gt. ni) then
  print*,'istart ',istart,' must be in range ',1,ni
  return
endif

tsa=miss_val
diff=miss_val

where(test .ne. miss_val .and. ref .ne. miss_val) diff=test-ref


  plus=miss_val ; minus=miss_val ; prms=miss_val ; mrms=miss_val
  pcount=0 ; mcount=0


ini=.true.

do j=istart,istart,increment

  jstart=j
  if(jstart .lt. 1) jstart=1
  jstartlen=j+left_maxlen-jstart
  jend=j+left_maxlen+right_maxlen-1
  if(jend .gt. ni) then
    jend=ni
    jendlen=ni-j
  else
    jendlen=jend-j-left_maxlen+1
  endif

  if(jstart .eq. 1 .or. jend .eq. ni .or. ini) then
    ip=.false.
    ipc=0
    im=.false.
    imc=0

    ip(1:jstartlen)=diff(jstart:j+left_maxlen-1) .ne. miss_val
    ipc=count(ip(1:jstartlen))
    im(1:jendlen)=diff(j+left_maxlen:jend) .ne. miss_val
    imc=count(im(1:jendlen))

  else
  endif

  minlen=min(left_maxlen,right_maxlen)
  if(minlen-ipc .lt. max_miss .and. minlen-imc .lt. max_miss) then
       
    do ibin=1,12
      iph(1:jstartlen)=ip(1:jstartlen) .and. bin(jstart:j+left_maxlen-1) .eq. ibin
      pmon(ibin)=count(iph(1:jstartlen))
      imh(1:jendlen)=im(1:jendlen) .and. bin(j+left_maxlen:jend) .eq. ibin
      mmon(ibin)=count(imh(1:jendlen))
      mmonhilf=mmon(ibin)*left_maxlen/right_maxlen
      pmonhilf=pmon(ibin)*right_maxlen/left_maxlen

      
      if(pmonhilf .lt. mmon(ibin)) then
        i=0
        k=0
        
        do while(i .lt. mmon(ibin)-pmonhilf)
          do while( .not. imh(jendlen-k) .and. k .lt. jendlen-1) 
            k=k+1
          enddo
          imh(jendlen-k)=.false.
          im(jendlen-k)=.false.
          i=i+1
        enddo
      endif
      if(mmonhilf .lt. pmon(ibin)) then
        i=0
        k=0
        do while(i .lt. pmon(ibin)-mmonhilf)
          do while( .not. iph(1+k) .and. k .lt. jstartlen-1) 
             k=k+1
          enddo
          iph(1+k)=.false.
          ip(1+k)=.false.
          i=i+1
        enddo
      endif
    enddo 
    ipc=count(ip(1:jstartlen))
    imc=count(im(1:jendlen))
    if(minlen-ipc .lt. max_miss .and. minlen-imc .lt. max_miss) then
    do ibin=1,12
      pmon(ibin)=count(ip(1:jstartlen) .and. bin(jstart:j+left_maxlen-1) .eq. ibin)
      mmon(ibin)=count(im(1:jendlen) .and. bin(j+left_maxlen:jend) .eq. ibin)
    enddo
       if(jstart .eq. 1 .or. jend .eq. ni .or. ini) then
         where(ip(1:jstartlen)) plusshort(1:jstartlen)=diff(jstart:j+left_maxlen-1)
         where(im(1:jendlen)) minusshort(1:jendlen)=diff(j+left_maxlen:jend)
         tp=sum(plusshort(1:jstartlen),ip(1:jstartlen))
         tm=sum(minusshort(1:jendlen),im(1:jendlen))
!         write(12,'(I6,2F5.2)') j,tp/ipc,tm/imc
         tpp=sum(plusshort(1:jstartlen)*plusshort(1:jstartlen),ip(1:jstartlen))
         tmm=sum(minusshort(1:jendlen)*minusshort(1:jendlen),im(1:jendlen))
         ini=.false.
       else
       endif

       qquer=(tp+tm)/(ipc+imc)
       rms=(tpp+tmm)/(ipc+imc)
       sigq=sqrt(rms-qquer*qquer)  
       
     if(sigq .gt. 0.001) then
       z1=(tp/ipc-qquer)/sigq
       z2=(tm/imc-qquer)/sigq
       tsa=ipc*z1*z1+imc*z2*z2
       plus=tp/ipc
       minus=tm/imc
       prms=tpp/ipc
       mrms=tmm/imc
     endif
    else
      ini=.true.
    endif
 else
   ini=.true.
 endif
   ini=.true.
   pcount=ipc
   mcount=imc

enddo


return
end subroutine meaneqsamp

! meaneqsamp calculates running mean SNHT test parameter tsa
! by dividing interval in two pieces and then calculating
! means, rms of subsamples
! means, rms are calculated recursively (= fairly efficiently)
!
! test=tested series
! ref=reference series
! ni size of test,ref
! istart= start index of analysis window
! istop=  end index of analysis window
! maxlen = length of averaging window
! miss_val = missing value indicator
! max_miss = maximum number of missing values in one of the two subsamples
! plus,minus = means of second, first subsample
! prms,mrms = rms of second, first subsample
! pcount,mcount = number of values used for calculating averages in second, first subsample
!
! the SNHT test statistic tsa is returned
!
! L. Haimberger, 7.7.2004
!
subroutine meaneqsamp2(test,ref,ni,istart,left_maxlen,right_maxlen,increment,miss_val,max_miss,critical_dates,ncritical,tsa,plus,minus,prms,mrms,pcount,mcount,bin) !

implicit none

integer         :: ni,left_maxlen,right_maxlen,max_miss,istart,jstart,jend,jstartlen,jendlen,i,k,l,minlen,offset
real(kind=JPRM) :: miss_val,fac

real(kind=JPRM),intent(in) :: test(ni),ref(ni)
integer,intent(in) :: bin(ni)
real(kind=JPRM) :: diff(ni),plus,minus,prms,mrms
real(kind=JPRM) :: tsa,innov(increment),plusshort(left_maxlen),minusshort(right_maxlen)
real(kind=JPRM) :: qquer,rms,sigq,z1,z2,ti,tti,tm,tmm,tp,tpp,eps

integer         :: pcount,mcount,ibin,pmon(12),mmon(12),ncritical,critical_dates(ncritical),stat
integer         :: increment,j,ipc,imc,ipnc,imnc,ipoc,imoc,jmi,jm2i,mmonhilf,pmonhilf,ipindex(left_maxlen,12),imindex(right_maxlen,12)
logical         :: ini,ip(left_maxlen),im(right_maxlen),ipn(increment),imn(increment),ipo(increment),imo(increment),lcut

!!$real(kind=JPRM),save,allocatable :: bin(:)
!!$
!!$
!!$if( .not. allocated(bin)) then
!!$  allocate(bin(ni),stat=stat)
!!$  do i=1,ni
!!$    bin(i)=modulo(i*1.,365.)
!!$    bin(i)=floor(bin(i)/365.*12+1)
!!$  enddo
!!$endif

if(istart .lt. 1-left_maxlen .or. istart .gt. ni) then
  print*,'istart ',istart,' must be in range ',1,ni
  return
endif

tsa=miss_val
diff=miss_val



  plus=miss_val ; minus=miss_val ; prms=miss_val ; mrms=miss_val
  pcount=0 ; mcount=0

l=0
do i=1,ni
if(test(i) .ne. miss_val .and. ref(i) .ne. miss_val) then
  diff(i)=test(i)-ref(i)
  l=l+1
endif
enddo
if(l .eq. 0) return
  
where(test .ne. miss_val .and. ref .ne. miss_val) diff=test-ref

ini=.true.

do j=istart,istart,increment

  jstart=j
  if(jstart .lt. 1) jstart=1
  jstartlen=j+left_maxlen-jstart
  jend=j+left_maxlen+right_maxlen-1
  if(jend .gt. ni) jend=ni
  jendlen=jend-j-left_maxlen+1

! breakpoint may be blurred -> avoid data surrounding breakpoint.
  minlen=min(left_maxlen,right_maxlen)
  offset=180
  ipc=0
  imc=0
  do while(offset .ge. 0 .and.  minlen-0.5*ipc .gt. max_miss )

    lcut=.false.
    ip(1:jstartlen)=diff(jstart:j+left_maxlen-1) .ne. miss_val
    if(offset .ne. 0) then
      im(1:1+offset)=.false.
      im(jendlen-offset:jendlen)=.false.
    endif
    ipc=count(ip(1:jstartlen))

    if(minlen-0.5*ipc-max_miss .gt. 30) then
      offset=offset-(minlen-0.5*ipc-max_miss)
    else
      offset=offset-30
    endif

    offset=offset-30

  enddo
  offset=180
  do while(offset .ge. 0 .and. minlen-0.5*imc .gt. max_miss)
    im(1:jendlen)=diff(j+left_maxlen:jend) .ne. miss_val
    if(offset .ne. 0) then
      im(1:1+offset)=.false.
      im(jendlen-offset:jendlen)=.false.
    endif
    imc=count(im(1:jendlen))

     if( minlen-0.5*imc-max_miss .gt. 30) then
      offset=offset-(minlen-0.5*imc-max_miss)
    else
      offset=offset-30
    endif

  enddo

  if(minlen-ipc .lt. max_miss .and. minlen-imc .lt. max_miss) then
      
      ibin=-1000
      pmon=0
      do i=1,jstartlen
        if(ip(i) ) then !.and. jstart+i-1 .le. ni .and.  jstart+i-1 .gt. 0
          ibin=bin(jstart+i-1)
          pmon(ibin)=pmon(ibin)+1
          ipindex(pmon(ibin),ibin)=i
        endif         
      enddo
     
      mmon=0
      do i=1,jendlen
        if(im(i)  ) then !.and. j+left_maxlen+i-1 .le. ni .and. j+left_maxlen+i-1 .gt. 0 
         ibin=bin(j+left_maxlen+i-1)
          mmon(ibin)=mmon(ibin)+1
          imindex(mmon(ibin),ibin)=i
        endif       
      enddo
      if(ibin .eq. -1000) then
        stop 'ibin problem'
      endif
!      mmonhilf=sum(mmon)/12
!      pmonhilf=sum(pmon)/12
      fac=sum(mmon)*1.0/sum(pmon)
      
    do ibin=1,12
      l=mmon(ibin)-floor(fac*pmon(ibin))
      if(l .gt. 1) then 
        do i=1,l
          im(imindex(mmon(ibin),ibin))=.false.
          mmon(ibin)=mmon(ibin)-1
          imc=imc-1
        enddo
      endif
      
      l=pmon(ibin)-floor(1./fac*mmon(ibin))
      if(l .gt. 1) then 
        do i=1,l
          ip(ipindex(i,ibin))=.false.
          pmon(ibin)=pmon(ibin)-1
          ipc=ipc-1
        enddo
      endif


    enddo 

!    ipc=count(ip(1:jstartlen))
!    imc=count(im(1:jendlen))

!write(*,*) 'Max:',istart+left_maxlen,left_maxlen,right_maxlen,ipc,imc

    if(minlen-ipc .lt. max_miss .and. minlen-imc .lt. max_miss) then

!if(fac .gt. 10 .or. fac .lt. 0.1) then
!  write(*,*) 'fac problem',fac,sum(mmon),sum(pmon)
!endif

!    do ibin=1,12
!      pmon(ibin)=count(ip(1:jstartlen) .and. bin(jstart:j+left_maxlen-1) .eq. ibin)
!      mmon(ibin)=count(im(1:jendlen) .and. bin(j+left_maxlen:jend) .eq. ibin)
!    enddo
       if(jstart .eq. 1 .or. jend .eq. ni .or. ini) then

         tp=0.
         tpp=0.
         do i=1,jstartlen
            if(ip(i)) then
              tp=tp+diff(jstart+i-1)
              tpp=tpp+diff(jstart+i-1)*diff(jstart+i-1)
            endif
         enddo
         tm=0.
         tmm=0.
         do i=1,jendlen
            if(im(i)) then
              tm=tm+diff(j+left_maxlen+i-1)
              tmm=tmm+diff(j+left_maxlen+i-1)*diff(j+left_maxlen+i-1)
            endif
         enddo
         ini=.false.
       else
       endif

     if(ipc*imc .gt. 0) then
       qquer=(tp+tm)/(ipc+imc)
       rms=(tpp+tmm)/(ipc+imc)
       if(rms-qquer*qquer .gt. 0) then
         sigq=sqrt(rms-qquer*qquer) 
       else
         sigq=0.
       endif 
     else
       sigq=0.
     endif 
       
     if(sigq .gt. 0.001) then
       z1=(tp/ipc-qquer)/sigq
       z2=(tm/imc-qquer)/sigq
       tsa=ipc*z1*z1+imc*z2*z2
       plus=tp/ipc
       minus=tm/imc
       prms=tpp/ipc
       mrms=tmm/imc
     endif
    else
      ini=.true.
    endif
 else
   ini=.true.
 endif
   ini=.true.
   pcount=ipc
   mcount=imc

enddo


return
end subroutine meaneqsamp2

! medeqsamp calculates running mean SNHT test parameter tsa
! by dividing interval in two pieces and then calculating
! means, rms of subsamples
! means, rms are calculated recursively (= fairly efficiently)
!
! test=tested series
! ref=reference series
! ni size of test,ref
! istart= start index of analysis window
! istop=  end index of analysis window
! maxlen = length of averaging window
! miss_val = missing value indicator
! max_miss = maximum number of missing values in one of the two subsamples
! plus,minus = means of second, first subsample
! prms,mrms = rms of second, first subsample
! pcount,mcount = number of values used for calculating averages in second, first subsample
!
! the SNHT test statistic tsa is returned
!
! L. Haimberger, 19.11.2013
!
subroutine meaneqsamp2_med(test,ref,ni,istart,left_maxlen,right_maxlen,increment,miss_val,max_miss,critical_dates,ncritical,tsa,plus,minus,prms,mrms,pcount,mcount,bin) !

use sort

implicit none

integer         :: ni,left_maxlen,right_maxlen,max_miss,istart,jstart,jend,jstartlen,jendlen,i,k,l,minlen,offset,ll
real(kind=JPRM) :: miss_val,fac

real(kind=JPRM),intent(in) :: test(ni),ref(ni)
integer,intent(in) :: bin(ni)
real(kind=JPRM) :: diff(ni),plus,minus,prms,mrms
real(kind=JPRM) :: tsa,innov(increment),plusshort(left_maxlen),minusshort(right_maxlen)
real(kind=JPRM) :: qquer,rms,sigq,z1,z2,ti,tti,tm,tmm,tp,tpp,eps

integer         :: pcount,mcount,ibin,pmon(12),mmon(12),ncritical,critical_dates(ncritical),stat
integer         :: increment,j,ipc,imc,ipnc,imnc,ipoc,imoc,jmi,jm2i,mmonhilf,pmonhilf,ipindex(left_maxlen,12),imindex(right_maxlen,12)
logical         :: ini,ip(left_maxlen),im(right_maxlen),ipn(increment),imn(increment),ipo(increment),imo(increment),lcut

integer,target,allocatable:: lindex(:)
real:: ldists(ni)
real,target,allocatable:: lldists(:)

!!$real(kind=JPRM),save,allocatable :: bin(:)
!!$
!!$
!!$if( .not. allocated(bin)) then
!!$  allocate(bin(ni),stat=stat)
!!$  do i=1,ni
!!$    bin(i)=modulo(i*1.,365.)
!!$    bin(i)=floor(bin(i)/365.*12+1)
!!$  enddo
!!$endif

if(istart .lt. 1-left_maxlen .or. istart .gt. ni) then
  print*,'istart ',istart,' must be in range ',1,ni
  return
endif

tsa=miss_val
diff=miss_val



  plus=miss_val ; minus=miss_val ; prms=miss_val ; mrms=miss_val
  pcount=0 ; mcount=0

l=0
do i=1,ni
if(test(i) .ne. miss_val .and. ref(i) .ne. miss_val) then
  diff(i)=test(i)-ref(i)
  l=l+1
endif
enddo
if(l .eq. 0) return
  
where(test .ne. miss_val .and. ref .ne. miss_val) diff=test-ref

ini=.true.

do j=istart,istart,increment

  jstart=j
  if(jstart .lt. 1) jstart=1
  jstartlen=j+left_maxlen-jstart
  jend=j+left_maxlen+right_maxlen-1
  if(jend .gt. ni) jend=ni
  jendlen=jend-j-left_maxlen+1

! breakpoint may be blurred -> avoid data surrounding breakpoint.
  minlen=min(left_maxlen,right_maxlen)
  offset=180
  ipc=0
  imc=0
  do while(offset .ge. 0 .and.  minlen-0.5*ipc .gt. max_miss )

    lcut=.false.
    ip(1:jstartlen)=diff(jstart:j+left_maxlen-1) .ne. miss_val
    if(offset .ne. 0) then
      im(1:1+offset)=.false.
      im(jendlen-offset:jendlen)=.false.
    endif
    ipc=count(ip(1:jstartlen))

    if(minlen-0.5*ipc-max_miss .gt. 30) then
      offset=offset-(minlen-0.5*ipc-max_miss)
    else
      offset=offset-30
    endif

    offset=offset-30

  enddo
  offset=180
  do while(offset .ge. 0 .and. minlen-0.5*imc .gt. max_miss)
    im(1:jendlen)=diff(j+left_maxlen:jend) .ne. miss_val
    if(offset .ne. 0) then
      im(1:1+offset)=.false.
      im(jendlen-offset:jendlen)=.false.
    endif
    imc=count(im(1:jendlen))

     if( minlen-0.5*imc-max_miss .gt. 30) then
      offset=offset-(minlen-0.5*imc-max_miss)
    else
      offset=offset-30
    endif

  enddo

  if(minlen-ipc .lt. max_miss .and. minlen-imc .lt. max_miss) then
      
      ibin=-1000
      pmon=0
      do i=1,jstartlen
        if(ip(i) ) then !.and. jstart+i-1 .le. ni .and.  jstart+i-1 .gt. 0
          ibin=bin(jstart+i-1)
          pmon(ibin)=pmon(ibin)+1
          ipindex(pmon(ibin),ibin)=i
        endif         
      enddo
     
      mmon=0
      do i=1,jendlen
        if(im(i)  ) then !.and. j+left_maxlen+i-1 .le. ni .and. j+left_maxlen+i-1 .gt. 0 
         ibin=bin(j+left_maxlen+i-1)
          mmon(ibin)=mmon(ibin)+1
          imindex(mmon(ibin),ibin)=i
        endif       
      enddo
      if(ibin .eq. -1000) then
        stop 'ibin problem'
      endif
!      mmonhilf=sum(mmon)/12
!      pmonhilf=sum(pmon)/12
      fac=sum(mmon)*1.0/sum(pmon)
      
    do ibin=1,12
      l=mmon(ibin)-floor(fac*pmon(ibin))
      if(l .gt. 1) then 
        do i=1,l
          im(imindex(mmon(ibin),ibin))=.false.
          mmon(ibin)=mmon(ibin)-1
          imc=imc-1
        enddo
      endif
      
      l=pmon(ibin)-floor(1./fac*mmon(ibin))
      if(l .gt. 1) then 
        do i=1,l
          ip(ipindex(i,ibin))=.false.
          pmon(ibin)=pmon(ibin)-1
          ipc=ipc-1
        enddo
      endif


    enddo 

!    ipc=count(ip(1:jstartlen))
!    imc=count(im(1:jendlen))

!write(*,*) 'Max:',istart+left_maxlen,left_maxlen,right_maxlen,ipc,imc

    if(minlen-ipc .lt. max_miss .and. minlen-imc .lt. max_miss) then

       if(jstart .eq. 1 .or. jend .eq. ni .or. ini) then

         tp=0.
         tpp=0.
         ll=0
         do i=1,jstartlen
            if(ip(i)) then
              ll=ll+1
              ldists(ll)=diff(jstart+i-1)+ll*1.d-12
              tpp=tpp+diff(jstart+i-1)*diff(jstart+i-1)
            endif
         enddo
         allocate(lldists(ll),lindex(ll))
         lldists=ldists(1:ll)
!!$OMP CRITICAL
         call qsort(lldists,lindex)
!!$OMP END CRITICAL
         plus=lldists(ll/2+1)
         tp=sum(lldists)/ll
!!$         if(abs(tp-plus) .gt. 5.) then
!!$           write(*,*) 'mean,median ',tp,plus,ll,sqrt(tpp/ll)
!!$         endif
         deallocate(lldists,lindex)
         tm=0.
         tmm=0.
         ll=0
         do i=1,jendlen
            if(im(i)) then
              ll=ll+1
              ldists(ll)=diff(j+left_maxlen+i-1)+ll*1.d-12
              tmm=tmm+diff(j+left_maxlen+i-1)*diff(j+left_maxlen+i-1)
            endif
         enddo
         allocate(lldists(ll),lindex(ll))
         lldists=ldists(1:ll)
         if(ll .eq. 1262) then
           write(*,*) lldists
         endif
!!$OMP CRITICAL
         call qsort(lldists,lindex)
!!$OMP END CRITICAL
         minus=lldists(ll/2+1)
         tm=sum(lldists)/ll
!!$         if(abs(tm-minus) .gt. 5.) then
!!$           write(*,*) 'mean,median ',tm,minus,ll,sqrt(tmm/ll)
!!$         endif
         deallocate(lldists,lindex)

         ini=.false.
       else
       endif

     if(ipc*imc .gt. 0) then
       qquer=tp+tm
       rms=(tpp+tmm)/(ipc+imc)
       if(rms-qquer*qquer .gt. 0) then
         sigq=sqrt(rms-qquer*qquer) 
       else
         sigq=0.
       endif 
     else
       sigq=0.
     endif 
       
     if(sigq .gt. 0.001) then
       z1=(plus-qquer)/sigq
       z2=(minus-qquer)/sigq
       tsa=ipc*z1*z1+imc*z2*z2
       prms=tpp/ipc
       mrms=tmm/imc
     endif
    else
      ini=.true.
    endif
 else
   ini=.true.
 endif
   ini=.true.
   pcount=ipc
   mcount=imc

enddo


return
end subroutine meaneqsamp2_med

! meaneqsamp calculates running mean SNHT test parameter tsa
! by dividing interval in two pieces and then calculating
! means, rms of subsamples
! means, rms are calculated recursively (= fairly efficiently)
!
! test=tested series
! ref=reference series
! ni size of test,ref
! istart= start index of analysis window
! istop=  end index of analysis window
! maxlen = length of averaging window
! miss_val = missing value indicator
! max_miss = maximum number of missing values in one of the two subsamples
! plus,minus = means of second, first subsample
! prms,mrms = rms of second, first subsample
! pcount,mcount = number of values used for calculating averages in second, first subsample
!
! the SNHT test statistic tsa is returned
!
! L. Haimberger, 7.7.2004
!
subroutine meaneqsamp3old(test,ref,ni,istart,left_maxlen,right_maxlen,increment,miss_val,min_count,critical_dates,ncritical,plus,minus,pcount,mcount,pms,mms) !

implicit none

integer,intent(in)         :: ni,left_maxlen,right_maxlen
integer                    :: max_miss,istart,jstart,jend,jstartlen,jendlen,i,k,l,minlen,pcount,mcount,min_count,offset,maxoff,mcfak,stat
real(kind=JPRM),intent(in) :: miss_val

real(kind=JPRM),intent(in) :: test(ni),ref(ni)
real(kind=JPRM) :: diff(ni),plus,minus,pms,mms,max
!real(kind=JPRM),optional ::pcorr,mcorr
real(kind=JPRM) :: innov(increment),plusshort(left_maxlen),minusshort(right_maxlen)
real(kind=JPRM) :: qquer,rms,sigq,z1,z2,ti,tti,tm,tmm,tp,tpp,eps,fak

integer         :: ibin,pmon(12),mmon(12),ncritical,critical_dates(ncritical),fillmax
integer         :: increment,j,ipc,imc,ipcsave,imcsave,ipnc,imnc,ipoc,imoc,jmi,jm2i,mmonhilf,pmonhilf,ipindex(left_maxlen,12),imindex(right_maxlen,12)
logical         :: ini,ip(left_maxlen),im(right_maxlen),ipn(increment),imn(increment),ipo(increment),imo(increment),lcut
real(kind=JPRM),save,allocatable :: bin(:)


if( .not. allocated(bin)) then
  allocate(bin(ni),stat=stat)
  do i=1,ni
    bin(i)=modulo(i*1.,365.)
    bin(i)=floor(bin(i)/365.*12+1)
  enddo
endif

  plus=miss_val ; minus=miss_val
  pms=miss_val ; mms=miss_val
!  pcorr=miss_val ; mcorr=miss_val
  pcount=0 ; mcount=0

if(istart .lt. 1-left_maxlen .or. istart .gt. ni) then
  print*,'istart ',istart,' must be in range ',1,ni
  return
endif

diff(istart:istart+left_maxlen+right_maxlen)=miss_val

max=0.
fillmax=0
do i=istart+fillmax+1,istart+left_maxlen+right_maxlen-fillmax-1
  if(test(i) .ne. miss_val) then
    j=0    
    if(ref(i) .ne. miss_val) then
      diff(i)=test(i)-ref(i)
    else
      do while(j .lt. fillmax+1 .and. diff(i) .eq. miss_val)
        j=j+1
        if(ref(i+j) .ne. miss_val) then
          diff(i)=test(i)-ref(i+j)
        else   
          if(ref(i-j) .ne. miss_val) diff(i)=test(i)-ref(i-j) 
        endif     
      enddo   
    endif 
  endif
!  if(diff(i) .ne. miss_val) then
!    if(abs(diff(i)) .gt. max) max=abs(diff(i))
!  endif
enddo


!where(test .ne. miss_val .and. ref .ne. miss_val) diff=test-ref



maxoff=180
mcfak=3
ini=.true.

do j=istart,istart,increment

  jstart=j
  if(jstart .lt. 1) jstart=1
  jstartlen=j+left_maxlen-jstart
  jend=j+left_maxlen+right_maxlen-1
  if(jend .gt. ni) jend=ni
  jendlen=jend-j-left_maxlen+1

    ipc=0
    do i=1,jstartlen
      if(diff(jstart+i-1) .ne. miss_val) then
        ip(i)=.true.
        ipc=ipc+1
      else
        ip(i)=.false.
      endif
    enddo
  offset=0
  do while(offset .le. maxoff .and. ipc .ge. mcfak*min_count)

    if(ip(1+offset)) then
      ip(1+offset)=.false.
      ipc=ipc-1
    endif
    if(ip(jstartlen-offset)) then
      ip(jstartlen-offset)=.false.
      ipc=ipc-1
    endif
    offset=offset+1
  enddo

    im(1:jendlen)=diff(j+left_maxlen:jend) .ne. miss_val
    imc=count(im(1:jendlen))
  offset=0
  do while(offset .le. maxoff .and. imc .ge. mcfak*min_count)
    if(im(1+offset)) then
      im(1+offset)=.false.
      imc=imc-1
    endif
    if(im(jendlen-offset)) then
      im(jendlen-offset)=.false.
      imc=imc-1
    endif
    offset=offset+1
  enddo

   if(ipc .gt. min_count .and. imc .gt. min_count) then  
      pmon=0
      do i=1,jstartlen
        if(ip(i)) then
          ibin=bin(jstart+i-1)
          pmon(ibin)=pmon(ibin)+1
          ipindex(pmon(ibin),ibin)=i
        endif         
      enddo
     
      mmon=0
      do i=1,jendlen
        if(im(i) ) then
         ibin=bin(j+left_maxlen+i-1)
          mmon(ibin)=mmon(ibin)+1
          imindex(mmon(ibin),ibin)=i
        endif       
      enddo

    fak=sum(pmon)*1.0_JPRM/(sum(mmon)*1.0_JPRM)
    ipcsave=ipc
    imcsave=imc
    do ibin=1,12
      mmonhilf=mmon(ibin)*fak
      pmonhilf=pmon(ibin)/fak

!      if(.false.) then      
      l=mmon(ibin)-pmonhilf
      if(l .gt. 0) then 
        do i=1,l
          im(imindex(mmon(ibin),ibin))=.false.
          mmon(ibin)=mmon(ibin)-1
          imc=imc-1
        enddo
      endif

      l=pmon(ibin)-mmonhilf
      if(l .gt. 0) then 
        do i=1,l
          ip(ipindex(i,ibin))=.false.
          pmon(ibin)=pmon(ibin)-1
          ipc=ipc-1
        enddo
      endif
!      endif

    enddo 

!    ipc=count(ip(1:jstartlen))
!    imc=count(im(1:jendlen))
!    ipc=0
!    imc=0
!    do i=1,jstartlen
!      if(ip(i)) ipc=ipc+1
!    enddo
!    do i=1,jendlen
!      if(im(i)) imc=imc+1
!    enddo
!write(*,*) 'Max:',istart+left_maxlen,max,left_maxlen,right_maxlen,ipc,imc

!    if(minlen-ipc .lt. max_miss .and. minlen-imc .lt. max_miss) then
!   if(ipc .gt. min_count .and. imc .gt. min_count) then 
!   pmon=0
!   mmon=0 
!    do ibin=1,12
!      pmon(ibin)=0
!      mmon(ibin)=0
!      do i=1,jstartlen
!        if(ip(i)) pmon(bin(jstart+i-1))=pmon(bin(jstart+i-1))+1    
!      enddo
!      do i=1,jendlen
!        if(im(i)) mmon(bin(j+left_maxlen-1+i))=mmon(bin(j+left_maxlen-1+i))+1   
!      enddo
!      pmon(ibin)=count(ip(1:jstartlen) .and. bin(jstart:j+left_maxlen-1) .eq. ibin)
!      mmon(ibin)=count(im(1:jendlen) .and. bin(j+left_maxlen:jend) .eq. ibin)
!    enddo
       if(jstart .eq. 1 .or. jend .eq. ni .or. ini) then

         tp=0.
         tpp=0.
         do i=1,jstartlen
            if(ip(i)) then
              tp=tp+diff(jstart+i-1)
              tpp=tpp+diff(jstart+i-1)*diff(jstart+i-1)
!debug              if(abs(diff(jstart+i-1)) .gt. 100.) then
!                 write(*,*) diff(jstart+i-1),jstart+i-1,' error in meaneqsamp3'
!              endif
            endif
         enddo
         tm=0.
         tmm=0.
         do i=1,jendlen
            if(im(i)) then
              tm=tm+diff(j+left_maxlen+i-1)
              tmm=tmm+diff(j+left_maxlen+i-1)*diff(j+left_maxlen+i-1)
!              if(abs(diff(j+left_maxlen+i-1)) .gt. 100.) then
!                 write(*,*) diff(j+left_maxlen+i-1),j+left_maxlen+i-1,' error in meaneqsamp3'
!              endif
            endif
         enddo
         ini=.false.
       else
       endif

       qquer=(tp+tm)/(ipc+imc)
       
       plus=tp/ipc
       minus=tm/imc
       pms=tpp/ipc
       mms=tmm/imc
!    else
!      ini=.true.
!    endif
 else
   ini=.true.
 endif
   ini=.true.
   pcount=ipcsave
   mcount=imcsave

enddo


return
end subroutine meaneqsamp3old

! meaneqsamp calculates running mean SNHT test parameter tsa
! by dividing interval in two pieces and then calculating
! means, rms of subsamples
! means, rms are calculated recursively (= fairly efficiently)
!
! test=tested series
! ref=reference series
! ni size of test,ref
! istart= start index of analysis window
! istop=  end index of analysis window
! maxlen = length of averaging window
! miss_val = missing value indicator
! max_miss = maximum number of missing values in one of the two subsamples
! plus,minus = means of second, first subsample
! prms,mrms = rms of second, first subsample
! pcount,mcount = number of values used for calculating averages in second, first subsample
!
! the SNHT test statistic tsa is returned
!
! L. Haimberger, 7.7.2004
!
subroutine meaneqsamp3(test,ref,ni,istart,left_maxlen,right_maxlen,increment,miss_val,min_count,mfak,critical_dates,ncritical,plus,minus,pcount,mcount,pms,mms,bin) !

implicit none

integer,intent(in)         :: ni,left_maxlen,right_maxlen
integer                    :: max_miss,istart,jstart,jend,jstartlen,jendlen,i,k,l,minlen
integer,intent(out)        :: pcount,mcount
integer                    :: min_count,offset,maxoff
real(kind=JPRM),intent(in) :: miss_val,mfak

real(kind=JPRM),intent(in) :: test(ni),ref(ni)
integer,intent(in) :: bin(ni)
real(kind=JPRM) :: diff(ni)
real(kind=JPRM),intent(out) :: plus,minus,pms,mms
real(kind=JPRM) :: pvar,mvar,psi,msi,pssi,mssi,max,pfak
!real(kind=JPRM),optional ::pcorr,mcorr
real(kind=JPRM) :: innov(increment),plusshort(left_maxlen),minusshort(right_maxlen)
real(kind=JPRM) :: qquer,rms,sigq,z1,z2,ti,tti,tm,tmm,tp,tpp,eps,fak
real(kind=JPRM) :: psum(12),msum(12),pssum(12),mssum(12),val

integer         :: ibin,pmon(12),mmon(12),ncritical,critical_dates(ncritical),fillmax,mthresh,gm,stat
integer         :: increment,j,ipc,imc,ipcsave,imcsave,ipnc,imnc,ipoc,imoc,jmi,jm2i,mmonhilf,pmonhilf,ipindex(left_maxlen,12),imindex(right_maxlen,12)
logical         :: ini,ip(left_maxlen),im(right_maxlen),ipn(increment),imn(increment),ipo(increment),imo(increment),lcut



  plus=miss_val ; minus=miss_val
  pms=miss_val ; mms=miss_val
  pcount=0 ; mcount=0
if(istart .lt. 1 .or. istart+left_maxlen+right_maxlen .gt. ni) then
  print*,'istart ',istart,' must be in range ',1,ni
  return
endif

diff(istart:istart+left_maxlen+right_maxlen)=miss_val

max=0.
fillmax=1
do i=istart+fillmax,istart+left_maxlen+right_maxlen-fillmax
  if(test(i) .ne. miss_val) then
    j=0    
    if(ref(i) .ne. miss_val) then
      diff(i)=test(i)-ref(i)
    else
      do while(j .lt. fillmax .and. diff(i) .eq. miss_val)
        j=j+1
        if(ref(i+j) .ne. miss_val) then
          diff(i)=test(i)-ref(i+j)
        else   
          if(ref(i-j) .ne. miss_val) diff(i)=test(i)-ref(i-j) 
        endif     
      enddo   
    endif 
  endif
!  if(diff(i) .ne. miss_val) then
!    if(abs(diff(i)) .gt. max) max=abs(diff(i))
!  endif
enddo


!where(test .ne. miss_val .and. ref .ne. miss_val) diff=test-ref



maxoff=90
pfak=0.7
mthresh=5

  j=istart
  jstart=j
  if(jstart .lt. 1) jstart=1
  jstartlen=j+left_maxlen-jstart
  jend=j+left_maxlen+right_maxlen-1
  if(jend .gt. ni) jend=ni
  jendlen=jend-j-left_maxlen+1

    ipc=0
    do i=1,jstartlen
      if(diff(jstart+i-1) .ne. miss_val) then
        ip(i)=.true.
        ipc=ipc+1
      else
        ip(i)=.false.
      endif
    enddo
  offset=0
  do while(offset .le. maxoff .and. ipc .gt. pfak*min_count+2)

    if(ip(1+offset)) then
      ip(1+offset)=.false.
      ipc=ipc-1
    endif
    if(ip(jstartlen-offset)) then
      ip(jstartlen-offset)=.false.
      ipc=ipc-1
    endif
    offset=offset+1
  enddo

  if(ipc .gt. pfak*min_count) then

    im(1:jendlen)=diff(j+left_maxlen:jend) .ne. miss_val
    imc=count(im(1:jendlen))
    offset=0
    do while(offset .le. maxoff .and. imc .gt. mfak*min_count+2)
      if(im(1+offset)) then
        im(1+offset)=.false.
        imc=imc-1
      endif
      if(im(jendlen-offset)) then
        im(jendlen-offset)=.false.
        imc=imc-1
      endif
      offset=offset+1
    enddo
  endif
!write(*,*) 'eqsamp3:',ipc,imc,min_count
if(ipc .gt. pfak*min_count .and. imc .gt. mfak*min_count) then  
      pmon=0
      do i=1,jstartlen
        if(ip(i)) then ! .and. jstart+i-1 .le. ni .and. jstart+i-1 .gt. 0 
          ibin=bin(jstart+i-1)
          pmon(ibin)=pmon(ibin)+1
          ipindex(pmon(ibin),ibin)=i
        endif         
      enddo
     
      mmon=0
      do i=1,jendlen
        if(im(i)) then ! .and. j+left_maxlen+i-1 .le. ni .and. j+left_maxlen+i-1 .gt. 0 
         ibin=bin(j+left_maxlen+i-1)
          mmon(ibin)=mmon(ibin)+1
          imindex(mmon(ibin),ibin)=i
        endif       
      enddo

    plus=0
    minus=0
    pms=0
    mms=0
    psum=0
    pssum=0
    msum=0
    mssum=0
    pcount=0
    mcount=0
    pvar=0
    mvar=0

    gm=0
    do ibin=1,12
      if(pmon(ibin) .gt. mthresh .and. mmon(ibin) .gt. mthresh) then
      gm=gm+1
      do i=1,pmon(ibin)
        val=diff(jstart+ipindex(i,ibin)-1)
        psum(ibin)=psum(ibin)+val
        pssum(ibin)=pssum(ibin)+val*val
      enddo

      do i=1,mmon(ibin)
        val=diff(j+left_maxlen+imindex(i,ibin)-1)
        msum(ibin)=msum(ibin)+val
        mssum(ibin)=mssum(ibin)+val*val
      enddo
 
      pcount=pcount+pmon(ibin)
      mcount=mcount+mmon(ibin)

      psi=psum(ibin)/pmon(ibin)
      msi=msum(ibin)/mmon(ibin)
      plus=plus+psi
      minus=minus+msi

      pssi=pssum(ibin)/pmon(ibin)
      mssi=mssum(ibin)/mmon(ibin)
      pms=pms+pssi
      mms=mms+mssi

      pvar=pvar+pssi-psi*psi
      mvar=mvar+mssi-msi*msi

      endif
    enddo
    
       
    if(gm .gt. 0) then
       plus=plus/gm
       minus=minus/gm
       mvar=mvar/gm
       pvar=pvar/gm
!       pms=pms/gm
!       mms=mms/gm
       pms=plus*plus+pvar
       mms=minus*minus+mvar
    else
       plus=miss_val
       minus=miss_val
       pms=miss_val
       mms=miss_val
    endif
 endif


return
end subroutine meaneqsamp3


! meaneqsamp calculates running mean SNHT test parameter tsa
! by dividing interval in two pieces and then calculating
! means, rms of subsamples
! means, rms are calculated recursively (= fairly efficiently)
!
! test=tested series
! ref=reference series
! ni size of test,ref
! istart= start index of analysis window
! istop=  end index of analysis window
! maxlen = length of averaging window
! miss_val = missing value indicator
! max_miss = maximum number of missing values in one of the two subsamples
! plus,minus = means of second, first subsample
! prms,mrms = rms of second, first subsample
! pcount,mcount = number of values used for calculating averages in second, first subsample
!
! the SNHT test statistic tsa is returned
!
! L. Haimberger, 7.7.2004
!
subroutine meaneqsamp3orig(test,ref,ni,istart,left_maxlen,right_maxlen,increment,miss_val,min_count,critical_dates,ncritical,plus,minus,pcount,mcount,pms,mms) !

implicit none

integer,intent(in)         :: ni,left_maxlen,right_maxlen
integer                    :: max_miss,istart,jstart,jend,jstartlen,jendlen,i,k,l,minlen,pcount,mcount,min_count,offset,maxoff,pfak,mfak
real(kind=JPRM),intent(in) :: miss_val

real(kind=JPRM),intent(in) :: test(ni),ref(ni)
real(kind=JPRM) :: diff(ni),plus,minus,pms,mms,max
!real(kind=JPRM),optional ::pcorr,mcorr
real(kind=JPRM) :: innov(increment),plusshort(left_maxlen),minusshort(right_maxlen)
real(kind=JPRM) :: qquer,rms,sigq,z1,z2,ti,tti,tm,tmm,tp,tpp,eps,fak
real(kind=JPRM) :: psum(12),msum(12),pssum(12),mssum(12),val

integer         :: ibin,pmon(12),mmon(12),ncritical,critical_dates(ncritical),fillmax,mthresh,gm,stat
integer         :: increment,j,ipc,imc,ipcsave,imcsave,ipnc,imnc,ipoc,imoc,jmi,jm2i,mmonhilf,pmonhilf,ipindex(left_maxlen,12),imindex(right_maxlen,12)
logical         :: ini,ip(left_maxlen),im(right_maxlen),ipn(increment),imn(increment),ipo(increment),imo(increment),lcut
real(kind=JPRM),save,allocatable :: bin(:)


if( .not. allocated(bin)) then
  allocate(bin(ni),stat=stat)
  do i=1,ni
    bin(i)=modulo(i*1.,365.)
    bin(i)=floor(bin(i)/365.*12+1)
  enddo
endif

  plus=miss_val ; minus=miss_val
  pms=miss_val ; mms=miss_val
!  pcorr=miss_val ; mcorr=miss_val
  pcount=0 ; mcount=0

if(istart .lt. 1-left_maxlen .or. istart .gt. ni) then
  print*,'istart ',istart,' must be in range ',1,ni
  return
endif

diff(istart:istart+left_maxlen+right_maxlen)=miss_val

max=0.
fillmax=0
do i=istart+fillmax,istart+left_maxlen+right_maxlen-fillmax
  if(test(i) .ne. miss_val) then
    j=0    
    if(ref(i) .ne. miss_val) then
      diff(i)=test(i)-ref(i)
    else
      do while(j .lt. fillmax .and. diff(i) .eq. miss_val)
        j=j+1
        if(ref(i+j) .ne. miss_val) then
          diff(i)=test(i)-ref(i+j)
        else   
          if(ref(i-j) .ne. miss_val) diff(i)=test(i)-ref(i-j) 
        endif     
      enddo   
    endif 
  endif
!  if(diff(i) .ne. miss_val) then
!    if(abs(diff(i)) .gt. max) max=abs(diff(i))
!  endif
enddo


!where(test .ne. miss_val .and. ref .ne. miss_val) diff=test-ref



maxoff=90
pfak=1
mthresh=5

  j=istart
  jstart=j
  if(jstart .lt. 1) jstart=1
  jstartlen=j+left_maxlen-jstart
  jend=j+left_maxlen+right_maxlen-1
  if(jend .gt. ni) jend=ni
  jendlen=jend-j-left_maxlen+1

    ipc=0
    do i=1,jstartlen
      if(diff(jstart+i-1) .ne. miss_val) then
        ip(i)=.true.
        ipc=ipc+1
      else
        ip(i)=.false.
      endif
    enddo
  offset=0
  do while(offset .le. maxoff .and. ipc .gt. pfak*min_count+2)

    if(ip(1+offset)) then
      ip(1+offset)=.false.
      ipc=ipc-1
    endif
    if(ip(jstartlen-offset)) then
      ip(jstartlen-offset)=.false.
      ipc=ipc-1
    endif
    offset=offset+1
  enddo

    im(1:jendlen)=diff(j+left_maxlen:jend) .ne. miss_val
    imc=count(im(1:jendlen))
  offset=0
  do while(offset .le. maxoff .and. imc .gt. mfak*min_count+2)
    if(im(1+offset)) then
      im(1+offset)=.false.
      imc=imc-1
    endif
    if(im(jendlen-offset)) then
      im(jendlen-offset)=.false.
      imc=imc-1
    endif
    offset=offset+1
  enddo

if(ipc .gt. pfak*min_count .and. imc .gt. mfak*min_count) then  
      pmon=0
      do i=1,jstartlen
        if(ip(i)) then ! .and. jstart+i-1 .le. ni .and. jstart+i-1 .gt. 0 
          ibin=bin(jstart+i-1)
          pmon(ibin)=pmon(ibin)+1
          ipindex(pmon(ibin),ibin)=i
        endif         
      enddo
     
      mmon=0
      do i=1,jendlen
        if(im(i)) then ! .and. j+left_maxlen+i-1 .le. ni .and. j+left_maxlen+i-1 .gt. 0 
         ibin=bin(j+left_maxlen+i-1)
          mmon(ibin)=mmon(ibin)+1
          imindex(mmon(ibin),ibin)=i
        endif       
      enddo

    plus=0
    minus=0
    pms=0
    mms=0
    psum=0
    pssum=0
    msum=0
    mssum=0
    pcount=0
    mcount=0
    gm=0
    do ibin=1,12
      if(pmon(ibin) .gt. mthresh .and. mmon(ibin) .gt. mthresh) then
      gm=gm+1
      do i=1,pmon(ibin)
        val=diff(jstart+ipindex(i,ibin)-1)
        psum(ibin)=psum(ibin)+val
        pssum(ibin)=pssum(ibin)+val*val
      enddo

      do i=1,mmon(ibin)
        val=diff(j+left_maxlen+imindex(i,ibin)-1)
        msum(ibin)=msum(ibin)+val
        mssum(ibin)=mssum(ibin)+val*val
      enddo
 
      pcount=pcount+pmon(ibin)
      mcount=mcount+mmon(ibin)

      plus=plus+psum(ibin)/pmon(ibin)
      minus=minus+msum(ibin)/mmon(ibin)

      pms=pms+pssum(ibin)/pmon(ibin)
      mms=mms+mssum(ibin)/mmon(ibin)
      endif
    enddo
       
    if(gm .gt. 0) then
       plus=plus/gm
       minus=minus/gm
       pms=pms/gm
       mms=mms/gm
    else
       plus=miss_val
       minus=miss_val
       pms=miss_val
       mms=miss_val
    endif
 endif


return
end subroutine meaneqsamp3orig

! snhtmov calculates running mean SNHT test parameter tsa
! by dividing interval in two pieces and then calculating
! means, rms of subsamples
! means, rms are calculated recursively (= fairly efficiently)
!
! test=tested series
! ref=reference series
! ni size of test,ref
! istart= start index of analysis window
! istop=  end index of analysis window
! maxlen = length of averaging window
! miss_val = missing value indicator
! max_miss = maximum number of missing values in one of the two subsamples
! plus,minus = means of second, first subsample
! prms,mrms = rms of second, first subsample
! pcount,mcount = number of values used for calculating averages in second, first subsample
!
! the SNHT test statistic tsa is returned
!
! L. Haimberger, 2.7.2014
!
subroutine snhteqsamp2y(test,ref,ni,istartorig,istoporig,maxlen,increment,miss_val,max_miss,critical_dates,ncritical,tsa,plus,minus,prms,mrms,pcount,mcount,bin) !

implicit none

integer         :: ni,maxlen,max_miss,istart,istop,jstart,jend,jstartlen,jendlen,i,k,l
integer,intent(in) :: istartorig,istoporig
real(kind=JPRM) :: miss_val

real(kind=JPRM),intent(in) :: test(ni),ref(ni)
integer,intent(in) :: bin(ni),ncritical,critical_dates(ncritical)
real(kind=JPRM) :: diff(ni),plus(ni),minus(ni),prms(ni),mrms(ni),mean(ni),square(ni)
real(kind=JPRM) :: tsa(ni),xm,xp,x,y,xy,sig

integer         :: pcount(ni),mcount(ni),gcount(ni),gindex(ni),ibin,pmon(12),mmon(12),stat
integer         :: increment,m2,j,ipc,imc,ipindex(maxlen/2,12),imindex(maxlen/2,12)

istart=istartorig
istop=istoporig
istart=1
istop=ni
if(ni .eq. 45000) istart=mod(20819,increment)+1

if(istart .lt. 1 .or. istop .gt. ni) then
  print*,'istart,istop ',istart,istop,' must be in range ',1,ni
  return
endif

m2=maxlen/2
tsa=miss_val

j=0
do i=1,ni
 if (test(i) .ne. miss_val .and. ref(i) .ne. miss_val) then
  j=j+1
  gindex(j)=i
  diff(i)=test(i)-ref(i)
  if (j>1) then
    mean(j)=mean(j-1)+diff(i)
    square(j)=square(j-1)+diff(i)*diff(i)
  else
    mean(j)=diff(i)
    square(j)=diff(i)*diff(i)
  endif
 endif
 gcount(i)=j
enddo

if(j .lt. 2*(m2-miss_val)) return

  plus=miss_val ; minus=miss_val ; prms=miss_val ; mrms=miss_val
  pcount=0 ; mcount=0

do k=m2-max_miss,ni-max_miss

  xm=k-m2
  if (xm .lt. 1) xm=1

  xp=k+m2
  if (xp .gt. ni) xp=ni

  pcount(k)=gcount(xp)-gcount(k)
  mcount(k)=gcount(k)-gcount(xm)

  if (gcount(k)-gcount(xm) .gt. m2-max_miss .and. &
     gcount(xp)-gcount(k) .gt. m2-max_miss) then
     
     
     x=(mean(gcount(k))-mean(gcount(xm)))/(gcount(k)-gcount(xm))
     y=(mean(gcount(xp))-mean(gcount(k)))/(gcount(xp)-gcount(k))
     xy=(mean(gcount(xp))-mean(gcount(xm)))/(gcount(xp)-gcount(xm))

     sig=(square(gcount(xp))-square(gcount(xm)))/(gcount(xp)-gcount(xm))
     
     if (sig .gt. 0) then
       sig=sqrt(sig-xy*xy)
       tsa(k)=((gcount(k)-gcount(xm))*(x-xy)*(x-xy)+(gcount(xp)-gcount(k))*(y-xy)*(y-xy))/sig
     else
       sig=0.0
       tsa(k)=0.
     endif

     plus(k)=y
     minus(k)=x
     prms(k)=(square(gcount(xp))-square(gcount(k)))/(gcount(xp)-gcount(k))
     mrms(k)=(square(gcount(k))-square(gcount(xm)))/(gcount(k)-gcount(xm))
     
  endif


enddo

return
end subroutine snhteqsamp2y

! snhtmov calculates running mean SNHT test parameter tsa
! by dividing interval in two pieces and then calculating
! means, rms of subsamples
! means, rms are calculated recursively (= fairly efficiently)
!
! test=tested series
! ref=reference series
! ni size of test,ref
! istart= start index of analysis window
! istop=  end index of analysis window
! maxlen = length of averaging window
! miss_val = missing value indicator
! max_miss = maximum number of missing values in one of the two subsamples
! plus,minus = means of second, first subsample
! prms,mrms = rms of second, first subsample
! pcount,mcount = number of values used for calculating averages in second, first subsample
!
! the SNHT test statistic tsa is returned
!
! L. Haimberger, 2.7.2014
!
subroutine snhteqsamp2z(test,ref,ni,istartorig,istoporig,maxlen,increment,miss_val,max_miss,critical_dates,ncritical,tsa,plus,minus,prms,mrms,pcount,mcount,bin) !

implicit none

integer         :: ni,maxlen,max_miss,istart,istop,jstart,jend,jstartlen,jendlen,i,k,l,jbin(12),psub,msub
integer,intent(in) :: istartorig,istoporig
real(kind=JPRM) :: miss_val

real(kind=JPRM),intent(in) :: test(ni),ref(ni)
integer,intent(in) :: bin(ni),ncritical,critical_dates(ncritical)
real(kind=JPRM) :: diff(ni),plus(ni),minus(ni),prms(ni),mrms(ni),mean(ni,12),square(ni,12)
real(kind=JPRM) :: tsa(ni),xm,xp,x,y,xy,sig,xx,yy,cratio,xpi,xmi

integer         :: pcount(ni),mcount(ni),gcount(ni,12),gindex(ni),ibin,pmon(12),mmon(12),stat
integer         :: increment,m2,j,ipc,imc,ipindex(maxlen/2,12),imindex(maxlen/2,12)

istart=istartorig
istop=istoporig
istart=1
istop=ni
if(ni .eq. 45000) istart=mod(20819,increment)+1

if(istart .lt. 1 .or. istop .gt. ni) then
  print*,'istart,istop ',istart,istop,' must be in range ',1,ni
  return
endif

m2=maxlen/2
tsa=miss_val

jbin=0
do i=1,ni
 if (test(i) .ne. miss_val .and. ref(i) .ne. miss_val) then
  jbin(bin(i))=jbin(bin(i))+1
  j=jbin(bin(i))
  diff(i)=test(i)-ref(i)
  if (j>1) then
    mean(j,bin(i))=mean(j-1,bin(i))+diff(i)
    square(j,bin(i))=square(j-1,bin(i))+diff(i)*diff(i)
  else
    mean(j,bin(i))=diff(i)
    square(j,bin(i))=diff(i)*diff(i)
  endif
 endif
 gcount(i,:)=jbin(:)
enddo

if(sum(jbin) .lt. 2*(m2-miss_val)) return

  plus=miss_val ; minus=miss_val ; prms=miss_val ; mrms=miss_val
  pcount=0 ; mcount=0

do k=m2-max_miss,ni-max_miss,increment

  xm=k-m2
  if (xm .lt. 1) xm=1

  xp=k+m2
  if (xp .gt. ni) xp=ni

  j=sum(gcount(k,:))
  pcount(k)=sum(gcount(xp,:))-j
  mcount(k)=j-sum(gcount(xm,:))
  msub=0
  psub=0

  if (pcount(k) .gt. m2-max_miss .and. &
     mcount(k) .gt. m2-max_miss) then
     
     cratio=float(mcount(k))/pcount(k)
    
     x=0
     y=0
     xx=0
     yy=0
     do i=1,12
       xmi=xm
       xpi=xp
       msub=gcount(k,i)-gcount(xm,i)-cratio*(gcount(xp,i)-gcount(k,i))
       if (msub .ge. 1) then
         do while(gcount(xmi,i)-gcount(xm,i) .lt. msub .and. xmi .lt. k-(m2-max_miss))
           xmi=xmi+mod(i+12-bin(xmi),12)*30+msub
         enddo
       endif
       psub=gcount(xp,i)-gcount(k,i)-(gcount(k,i)-gcount(xm,i))/cratio
       if (psub .ge. 1) then
         do while(gcount(xp,i)-gcount(xpi,i) .lt. psub .and. xpi .gt. k+(m2-max_miss))
           xpi=xpi-mod(bin(xpi)+12-i,12)*30-psub
         enddo
       endif
       if(psub .lt. 1) psub=0
       if(msub .lt. 1) msub=0

      if (gcount(k,i) .gt. 0) then
          x=x+mean(gcount(k,i),i)
          y=y-mean(gcount(k,i),i)
          xx=xx+square(gcount(k,i),i)
          yy=yy-square(gcount(k,i),i)
       endif
       if (gcount(xmi,i) .gt. 0) then
         x=x-mean(gcount(xmi,i),i)
        xx=xx-square(gcount(xmi,i),i)
       endif 
       if (gcount(xpi,i) .gt. 0) then
       y=y+mean(gcount(xpi,i),i)
       yy=yy+square(gcount(xpi,i),i)
       endif
     enddo
     xy=(x+y)/(pcount(k)+mcount(k)-msub-psub)
     sig=(xx+yy)/(pcount(k)+mcount(k)-msub-psub)
     x=x/(mcount(k)-msub)
     y=y/(pcount(k)-psub)
     xx=xx/(mcount(k)-msub)
     yy=yy/(pcount(k)-psub)

     
     if (sig .gt. 0) then
       sig=sqrt(sig-xy*xy)
       tsa(k)=((mcount(k)-msub)*(x-xy)*(x-xy)+(pcount(k)-psub)*(y-xy)*(y-xy))/sig
     else
       sig=0.0
       tsa(k)=0.
     endif

     plus(k:k+increment-1)=y
     minus(k:k+increment-1)=x
     prms(k:k+increment-1)=yy
     mrms(k:k+increment-1)=xx
     tsa(k:k+increment-1)=tsa(k)
     
  endif
  pcount(k:k+increment-1)=pcount(k)-psub
  mcount(k:k+increment-1)=mcount(k)-msub


enddo

do l=1,ncritical
  tsa(critical_dates(l)-maxlen/2+max_miss:critical_dates(l)+maxlen/2-max_miss)=miss_val
  plus(critical_dates(l)-maxlen/2+max_miss:critical_dates(l)+maxlen/2-max_miss)=miss_val
  minus(critical_dates(l)-maxlen/2+max_miss:critical_dates(l)+maxlen/2-max_miss)=miss_val
  prms(critical_dates(l)-maxlen/2+max_miss:critical_dates(l)+maxlen/2-max_miss)=miss_val
  mrms(critical_dates(l)-maxlen/2+max_miss:critical_dates(l)+maxlen/2-max_miss)=miss_val
enddo

return
end subroutine snhteqsamp2z

end module homtests
