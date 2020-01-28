module rfcomp_1

use homtests
use rfcor

contains


subroutine winsor_mean(rcpara,spagarr1,dists,istat,index,spagarr,cspagarr,cspagarrl,distarr,distarrl,tbi)

implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer istat,i,ibl,ib,l

integer index(rcpara%statmax)
real :: dists(rcpara%statmax),distc(rcpara%cachemax),distarr(rcpara%cachemax),distarrl(rcpara%cachemax)

integer :: tbi,tbisave,ipmax,iparmax,ip,ipar
real(kind=JPRM) :: spaghilf(rcpara%cachemax),darrc(rcpara%brmax,rcpara%pmax,rcpara%parmax),xdrarrc(rcpara%brmax,rcpara%pmax,rcpara%parmax),tbindexn(rcpara%brmax),darr1(rcpara%brmax,rcpara%pmax,rcpara%parmax),darrcl(rcpara%brmax,rcpara%pmax,rcpara%parmax)


real(kind=JPRM)  ::   spagarr(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax),spagarr1(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax)
real(kind=JPRM)  ::   cspagarr(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax),cspagarrl(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax)

integer         :: icache(rcpara%statmax+1),icistat,istati

real meanhilf,maxhilf,meanhilf1,maxhilf1,dhilf,xdr,edr,swap,dhilfc,xdrc,edrc,swapd
integer enn,caunt(rcpara%pmax,rcpara%parmax),m,fpr,idx,ennc,spagindex(rcpara%brmax,rcpara%pmax,rcpara%parmax),sub,sub1,gson,n_critical,critical_dates(1),spagindexl(rcpara%brmax,rcpara%pmax,rcpara%parmax)
real :: spaghilf1(rcpara%cachemax),hilf,hilf13

do ipar=1,rcpara%parmax
  do ib=2,tbi
    spaghilf=rcpara%miss_val
    caunt=0
    idx=0
    do ip=1,rcpara%pmax
      spaghilf=spagarr1(ib,ip,ipar,:)
      caunt(ip,ipar)=count(spaghilf .ne. rcpara%miss_val)
      if (caunt(ip,ipar) .gt. 10) then 

       l=1
      do istat=1,rcpara%cachemax
       if (spaghilf(istat) .ne. rcpara%miss_val) then
        spaghilf1(l)=spaghilf(istat)
        distc(l)=dists(index(istat))
        l=l+1
       endif
      enddo
       l=1
       do m=1,caunt(ip,ipar)
       do l=m,caunt(ip,ipar)
         if (spaghilf1(l) .gt. spaghilf1(m)) then
            swap=spaghilf1(l)
            swapd=distc(l)
            spaghilf1(l)=spaghilf1(m)
            distc(l)=distc(m)
            spaghilf1(m)=swap
            distc(m)=swapd
         endif
       enddo
       enddo
       fpr=int((caunt(ip,ipar)/100.)*5.+1/2)
       idx=caunt(ip,ipar)-2*fpr-1
       cspagarr(ib,ip,ipar,1:idx)=spaghilf1((fpr+1):(caunt(ip,ipar)-fpr-1))
       distarr(1:idx)=distc((fpr+1):(caunt(ip,ipar)-fpr-1))
 
      endif
   enddo
  enddo
enddo


cspagarrl=rcpara%miss_val
do ipar=1,rcpara%parmax
  do ib=2,tbi
    spaghilf=rcpara%miss_val
    caunt=0
    idx=0
    do ip=1,rcpara%pmax
      spaghilf=spagarr(ib,ip,ipar,:)
      caunt(ip,ipar)=count(spaghilf .ne. rcpara%miss_val)
      if (caunt(ip,ipar) .gt. 20) then 

       l=1
      do istat=1,rcpara%cachemax
       if (spaghilf(istat) .ne. rcpara%miss_val) then
        spaghilf1(l)=spaghilf(istat)
        distc(l)=dists(index(istat))
        l=l+1
       endif
      enddo
       l=1
       do m=1,caunt(ip,ipar)
       do l=m,caunt(ip,ipar)
         if (spaghilf1(l) .gt. spaghilf1(m)) then
            swap=spaghilf1(l)
            swapd=distc(l)
            spaghilf1(l)=spaghilf1(m)
            distc(l)=distc(m)
            spaghilf1(m)=swap
            distc(m)=swapd
         endif
       enddo
       enddo
       fpr=int((caunt(ip,ipar)/100.)*5.+1/2)
       idx=caunt(ip,ipar)-2*fpr-1
       cspagarrl(ib,ip,ipar,1:idx)=spaghilf1((fpr+1):(caunt(ip,ipar)-fpr-1))
       distarrl(1:idx)=distc((fpr+1):(caunt(ip,ipar)-fpr-1))
 
      endif
   enddo
  enddo
enddo

return
end subroutine winsor_mean


subroutine common_mean(rcpara,spagarr,darr,spagindex,obssig,ib,minst) !


implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer istat,i,ibl,tbi,ip,ipar,ib,l,minst

real(kind=JPRM)  ::   spagarr(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax),obssig(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax)
real(kind=JPRM) :: darr(rcpara%brmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: dhilf,edr,min,max
integer enn,spagindex(rcpara%brmax,rcpara%pmax,rcpara%parmax)


do ipar=1,rcpara%parmax
!  do ib=1,tbi
    do ip=1,rcpara%pmax
      dhilf=0
      enn=0
      edr=0
      l=0
stat: do istat=1,rcpara%cachemax
        if (spagarr(ib,ip,ipar,istat) .ne. rcpara%miss_val  .and. obssig(ib,ip,ipar,istat) .ne. rcpara%miss_val) then
          dhilf=dhilf+spagarr(ib,ip,ipar,istat)/obssig(ib,ip,ipar,istat)  
          edr=edr+1./obssig(ib,ip,ipar,istat)
          enn=enn+1
          l=l+1
        endif
       if(l .eq. minst) then 
!         spagarr(ib,ip,ipar,istat+1:rcpara%cachemax)=rcpara%miss_val
         exit stat
       endif
      enddo stat
      if (dhilf .ne. 0 .and. enn .ne. 0) then
        darr(ib,ip,ipar)=dhilf/edr
        spagindex(ib,ip,ipar)=enn 
      else     
        darr(ib,ip,ipar)=-999.
      endif
   enddo
!  enddo
enddo

return
end subroutine common_mean


subroutine rweight_mean(rcpara,wmonrs,spagarr,lspagarr,xdrarr,xdsarr,index,obssig,dists,ib,typ,minst,istatmax) !

type(rasocor_namelist),intent(in) :: rcpara

integer istat,i,ibl,l,lmin,lmax,istatmax

integer,intent(in) :: index(rcpara%statmax)
real :: dists(rcpara%statmax)

integer :: tbi,tbisave,ipmax,iparmax,ip,ipar,gindex(rcpara%statmax),wmonrs(rcpara%statmax)
real(kind=JPRM) :: xdrarr(rcpara%brmax,rcpara%pmax,rcpara%parmax),xdsarr(rcpara%brmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM)  ::   spagarr(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax),obssig(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax)
logical lspagarr(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax)
integer         :: icache(rcpara%statmax+1),icistat,istati,typ
real(kind=JPRM) meanhilf,maxhilf,meanhilf1,maxhilf1,dhilf,xdr,edr,swap,dhilfc,xdrc,edrc,swapd
integer enn,caunt(rcpara%pmax,rcpara%parmax),m,fpr,idx,ennc
real(kind=JPRM) :: spaghilf1(rcpara%cachemax),hilf,hilf13,weight(rcpara%cachemax,rcpara%pmax),radius(rcpara%pmax),max,min,w,xds

radius=3*(/200.,200.,200.,200.,200.,200.,500.,500.,500.,500.,500.,500.,500.,500.,500.,500./)

!radius=radius*3.

do istat=1,istatmax
fak=exp(-(dists(index(istat))*3.1415729/180.*6370.)/radius(7))
 do ip=1,rcpara%pmax
   weight(istat,ip)=fak
!  weight(istat)=1./dists(index(istat))
 enddo
enddo

do ipar=1,rcpara%parmax
!  do ib=1,tbi
    do ip=1,rcpara%pmax
     xdr=0
     edr=0
     l=0
stat: do istat=1,istatmax
        if(spagarr(ib,ip,ipar,istat) .ne. rcpara%miss_val .and. dists(index(istat)) .ne. 0 .and. obssig(ib,ip,ipar,istat) .ne. rcpara%miss_val .and. lspagarr(ib,ip,ipar,istat)) then
          l=l+1
          gindex(l)=istat
        endif
       if(l .eq. minst) then 
         exit stat
       endif
      enddo stat

max=-1.e30
min=1e30
lmin=0
lmax=0
do istat=1,l
  if(spagarr(ib,ip,ipar,gindex(istat)) .gt. max) then
    lmax=istat
    max=spagarr(ib,ip,ipar,gindex(istat)) 
  endif
  if(spagarr(ib,ip,ipar,gindex(istat)) .lt. min) then
    lmin=istat
    min=spagarr(ib,ip,ipar,gindex(istat)) 
  endif
enddo
if(l .gt. 3 .and. lmin .ne. 0) spagarr(ib,ip,ipar,gindex(lmin))=rcpara%miss_val
if(l .gt. 3 .and. lmax .ne. 0) spagarr(ib,ip,ipar,gindex(lmax))=rcpara%miss_val

l=0
stat2: do istat=1,istatmax
       if ((spagarr(ib,ip,ipar,istat) .ne. rcpara%miss_val) .and. dists(index(istat)) .ne. 0 .and. lspagarr(ib,ip,ipar,istat) .and. obssig(ib,ip,ipar,istat) .ne. 0) then
         w=weight(istat,ip)/obssig(ib,ip,ipar,istat)**2
         xdr=xdr+spagarr(ib,ip,ipar,istat)*w
         edr=edr+w
         l=l+1
       endif
       if(l .eq. minst) then 
!         spagarr(ib,ip,ipar,istat+1:rcpara%cachemax)=rcpara%miss_val
         exit stat2
       endif
      enddo stat2
     if (xdr .ne. 0 .and. edr .ne. 0) then
      xdrarr(ib,ip,ipar)=xdr/edr

l=0
xds=0.
ws=0.
stat5: do istat=1,istatmax
       if ((spagarr(ib,ip,ipar,istat) .ne. rcpara%miss_val) .and. dists(index(istat)) .ne. 0 .and. lspagarr(ib,ip,ipar,istat) .and. obssig(ib,ip,ipar,istat) .ne. 0) then
         w=weight(istat,ip)/obssig(ib,ip,ipar,istat)**2
         ws=ws+w*w
         xds=xds+(xdrarr(ib,ip,ipar)-spagarr(ib,ip,ipar,istat))**2*w
         l=l+1
       endif
       if(l .eq. minst) then 
!         spagarr(ib,ip,ipar,istat+1:rcpara%cachemax)=rcpara%miss_val
         exit stat5
       endif
      enddo stat5
      if(edr*edr-ws .gt. 0) then
        xdsarr(ib,ip,ipar)=sqrt(xds/(edr*edr-ws))
      else
        xdsarr(ib,ip,ipar)=rcpara%miss_val
      endif
!      write(*,'(3I3,F6.3,A20)') ib,ip,ipar,xdrarr(ib,ip,ipar),'used same level'
     else
! use also break estimates of level below current level
       l=0
       xdr=0
       edr=0
stat3: do istat=1,istatmax
       if ((spagarr(ib,ip,ipar,istat) .ne. rcpara%miss_val) .and. dists(index(istat)) .ne. 0 .and. obssig(ib,ip,ipar,istat) .ne. rcpara%miss_val .and. obssig(ib,ip,ipar,istat) .ne. 0.) then
         w=weight(istat,ip)/obssig(ib,ip,ipar,istat)**2
         xdr=xdr+spagarr(ib,ip,ipar,istat)*w
         edr=edr+w
         l=l+1
       endif
       if(l .eq. minst) then 
!         spagarr(ib,ip,ipar,istat+1:rcpara%cachemax)=rcpara%miss_val
         exit stat3
       endif
      enddo stat3
     if (xdr .ne. 0 .and. edr .ne. 0) then
      xdrarr(ib,ip,ipar)=xdr/edr
l=0
xds=0.
ws=0.
stat6: do istat=1,istatmax
       if ((spagarr(ib,ip,ipar,istat) .ne. rcpara%miss_val) .and. dists(index(istat)) .ne. 0 .and. lspagarr(ib,ip,ipar,istat) .and. obssig(ib,ip,ipar,istat) .ne. 0) then
         w=weight(istat,ip)/obssig(ib,ip,ipar,istat)**2
         ws=ws+w*w
         xds=xds+(xdrarr(ib,ip,ipar)-spagarr(ib,ip,ipar,istat))**2*w
         l=l+1
       endif
       if(l .eq. minst) then 
!         spagarr(ib,ip,ipar,istat+1:rcpara%cachemax)=rcpara%miss_val
         exit stat6
       endif
      enddo stat6
      if(edr*edr-ws .gt. 0) then
        xdsarr(ib,ip,ipar)=sqrt(xds/(edr*edr-ws))
      else
        xdsarr(ib,ip,ipar)=rcpara%miss_val
      endif

      write(*,'(I5,3I3,A21)') wmonrs(index(1)),ib,ip,ipar,' used one level below'
     else
      xdrarr(ib,ip,ipar)=-999.
      xdsarr(ib,ip,ipar)=-999.
     endif
     endif
    enddo
!  enddo
enddo


return
end subroutine rweight_mean


logical function goodest(rcpara,obssig,obscorrp,obscorrm,pcount,mcount,left_maxlen,right_maxlen,tpcount,tmcount,min_count,ib,ip,ipar,istat,wmonrs,index,innov,iter,ldbg)

implicit none
   
type(rasocor_namelist),intent(in) :: rcpara

integer,intent(in) :: min_count,ib,ip,ipar,istat,wmonrs(rcpara%statmax),index(rcpara%statmax),innov,iter
integer,intent(in) :: left_maxlen(rcpara%brmax,rcpara%parmax,rcpara%cachemax),right_maxlen(rcpara%brmax,rcpara%parmax,rcpara%cachemax)
integer,intent(in) :: tpcount(rcpara%brmax,rcpara%pmax,rcpara%parmax),tmcount(rcpara%brmax,rcpara%pmax,rcpara%parmax)
integer,intent(in) :: pcount(rcpara%brmax,rcpara%pmax,rcpara%parmax),mcount(rcpara%brmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM),intent(in) :: obssig(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax),obscorrp(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax),obscorrm(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax)
logical lhilf,ldbg
real(kind=JPRM) :: fak
integer :: mfak

lhilf=obssig(ib,ip,ipar,istat) .ne. rcpara%miss_val 
goodest=lhilf
return
fak=0.7
!mfak=3./iter

if(lhilf) then 
 
  if(ip .lt. 14) then
  lhilf=  &
        pcount(ib,ip,ipar) .gt. fak*tpcount(ib,ip,ipar) .and. mcount(ib,ip,ipar) .gt. fak*tmcount(ib,ip,ipar) .or. pcount(ib,ip,ipar) .ge. min_count  .and. mcount(ib,ip,ipar) .ge. min_count

!  if(innov .eq. 2) then
!    lhilf= lhilf .and. obssig(ib,ip,ipar,istat) .lt. 0.5
!  endif
  endif

!        pcount(ib,ip,ipar) .gt. 0*left_maxlen(ib,ipar,istat) .and. mcount(ib,ip,ipar) .gt. fak*right_maxlen(ib,ipar,istat) .or. &
!obssig(ib,ip,ipar,istat) .lt. 0.3 .or.
  if(innov .eq. 2) then
!    lhilf= lhilf .and. obssig(ib,ip,ipar,istat) .lt. 0.3
  else
!    lhilf= lhilf .and. (obssig(ib,ip,ipar,istat) .lt. 4*obssig(ib,ip,ipar,1) .or. obssig(ib,ip,ipar,istat) .lt. 1.0)
  endif   
  if(ldbg) then
    if(lhilf) write(*,'(A8,2I6,A3,I3,A3,2I3,A5,I2,F8.3,6I5)') 'accepted',wmonrs(index(1)),wmonrs(index(istat)),' ib',ib,' ip',ip,ipar,'innov',innov,&
obssig(ib,ip,ipar,istat),pcount(ib,ip,ipar),mcount(ib,ip,ipar),tpcount(ib,ip,ipar),tmcount(ib,ip,ipar),left_maxlen(ib,ipar,istat),right_maxlen(ib,ipar,istat)
    if(.not. lhilf ) write(*,'(A8,2I6,A3,I3,A3,2I3,A5,I2,F8.3,6I5,2F8.3)') 'rejected',wmonrs(index(1)),wmonrs(index(istat)),' ib',ib,' ip',ip,ipar,'innov',innov,obssig(ib,ip,ipar,istat),pcount(ib,ip,ipar),mcount(ib,ip,ipar),&
tpcount(ib,ip,ipar),tmcount(ib,ip,ipar),left_maxlen(ib,ipar,istat),right_maxlen(ib,ipar,istat),obssig(ib,ip,ipar,istat),obssig(ib,ip,ipar,1)
  endif
endif

goodest=lhilf

return
end function goodest


subroutine estimate(rcpara,ttm,cachehilf,left_maxlen,right_maxlen,tbindex,min_count,mfak,dists,wmonrs,wmolons,wmolats,index,obsspagarr,lobsspagarr,&
obssig,obscorrp,obscorrm,obsplus,obsminus,pcount,mcount,obspms,obsmms,trasocorrs,tpcount,tmcount,ib,ip,ipar,istat,iter,innov,try,ldebug) !

implicit none
   
type(rasocor_namelist),intent(in) :: rcpara

integer istat,i,j,statnr,ib,minst,l,min_count,n_critical,critical_dates(1),innov,iter,ip,ipar,try
integer,intent(in) :: wmonrs(rcpara%statmax)
real(kind=JPRM),intent(in) :: wmolons(rcpara%statmax),wmolats(rcpara%statmax),mfak
real,intent(in) :: dists(rcpara%statmax)
integer,intent(in) :: index(rcpara%statmax),tbindex(rcpara%brmax)

real(kind=JPRM),intent(in) :: ttm(rcpara%nmax,rcpara%pmax,rcpara%parmax),trasocorrs(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: obsspagarr(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax)
integer,intent(in) :: left_maxlen(rcpara%brmax,rcpara%parmax,rcpara%cachemax),right_maxlen(rcpara%brmax,rcpara%parmax,rcpara%cachemax)
real(kind=JPRM) :: obsplus(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax),obsminus(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax),obssig(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax)
real(kind=JPRM) :: obspms(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax),obsmms(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax)
real(kind=JPRM) :: obscorrm(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax),obscorrp(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax)
logical :: lobsspagarr(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax)
integer :: tpcount(rcpara%brmax,rcpara%pmax,rcpara%parmax),tmcount(rcpara%brmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM),intent(in) :: cachehilf(rcpara%nmax) 
integer :: pcount(rcpara%brmax,rcpara%pmax,rcpara%parmax),mcount(rcpara%brmax,rcpara%pmax,rcpara%parmax)
logical ldebug,ge

!             write(*,*) 'Station ',wmonrs(index(1)),istat,ib,ipar,ip
             call meaneqsamp3(ttm(:,ip,ipar),cachehilf,rcpara%nmax,tbindex(ib)-left_maxlen(ib,ipar,istat),left_maxlen(ib,ipar,istat),right_maxlen(ib,ipar,istat),rcpara%snht_increment,rcpara%miss_val,min_count,mfak,&
critical_dates,n_critical,obsplus(ib,ip,ipar,istat),obsminus(ib,ip,ipar,istat),pcount(ib,ip,ipar),mcount(ib,ip,ipar),obspms(ib,ip,ipar,istat),obsmms(ib,ip,ipar,istat),rcpara%month) !in file homtests.f90 line 1278 
 

             call spagstat(obsplus(ib,ip,ipar,istat),obsminus(ib,ip,ipar,istat),pcount(ib,ip,ipar),mcount(ib,ip,ipar),obspms(ib,ip,ipar,istat),obsmms(ib,ip,ipar,istat),obssig(ib,ip,ipar,istat),rcpara%miss_val) !in file rfcor.f90 line 26

             if (goodest(rcpara,obssig,obscorrp,obscorrm,pcount,mcount,left_maxlen,right_maxlen,tpcount,tmcount,min_count,ib,ip,ipar,istat,wmonrs,index,innov,iter,ldebug))  then
               ge=.true.

               if(innov .eq. 4) then
                  obsspagarr(ib,ip,ipar,istat)=-(obsplus(ib,ip,ipar,istat)-obsminus(ib,ip,ipar,istat))
               else
                  obsspagarr(ib,ip,ipar,istat)=-(obsplus(ib,ip,ipar,istat)-obsminus(ib,ip,ipar,istat))
               endif
               

      if(ldebug) write(*,'(A6,2I6.5,A6,I2,A4,I2,A6,I2,A4,I2,A4,I2,A5,I2,L2,2I5,4F8.3)') 'count ',wmonrs(index(1)),wmonrs(index(istat)),' innov',innov,' try',try,' iter ',iter,' ib ',ib,' ip ',ip,' ipar',ipar,ge,pcount(ib,ip,ipar),&
mcount(ib,ip,ipar),obsplus(ib,ip,ipar,istat),obsminus(ib,ip,ipar,istat),obsspagarr(ib,ip,ipar,istat),trasocorrs(tbindex(ib)+1,ip,ipar)-trasocorrs(tbindex(ib)-1,ip,ipar)

            else
             ge=.false.
!   .and. &
!                 obsspagarr(ib,ip+1,ipar,istat) .ne. obsspagarr(ib,ip+2,ipar,istat)
      if(ldebug) write(*,'(A6,2I6.5,A6,I2,A4,I2,A6,I2,A4,I2,A4,I2,A5,I2,L2,2I5,4F8.3)') 'count ',wmonrs(index(1)),wmonrs(index(istat)),' innov',innov,' try',try,' iter ',iter,' ib ',ib,' ip ',ip,' ipar',ipar,ge,pcount(ib,ip,ipar),&
mcount(ib,ip,ipar),obsplus(ib,ip,ipar,istat),obsminus(ib,ip,ipar,istat),obsspagarr(ib,ip,ipar,istat),trasocorrs(tbindex(ib)+1,ip,ipar)-trasocorrs(tbindex(ib)-1,ip,ipar)

            if(ip .lt. 8 .and. obsspagarr(ib,ip+1,ipar,istat) .ne. rcpara%miss_val .and. &
                 obsspagarr(ib,ip+1,ipar,istat) .ne. obsspagarr(ib,ip+2,ipar,istat) .and.&
                 trasocorrs(tbindex(ib)+1,ip,ipar)-trasocorrs(tbindex(ib)-1,ip,ipar) .ne. 0.) then
                obsspagarr(ib,ip,ipar,istat)=obsspagarr(ib,ip+1,ipar,istat)
                lobsspagarr(ib,ip,ipar,istat)=.false.
                pcount(ib,ip,ipar)=pcount(ib,ip+1,ipar)
                mcount(ib,ip,ipar)=mcount(ib,ip+1,ipar)
!                tpcount(ib,ip,ipar)=tpcount(ib,ip+1,ipar)
!                tmcount(ib,ip,ipar)=tmcount(ib,ip+1,ipar)
                obssig(ib,ip,ipar,istat)=obssig(ib,ip+1,ipar,istat)
                if ( .not. goodest(rcpara,obssig,obscorrp,obscorrm,pcount,mcount,left_maxlen,right_maxlen,tpcount,tmcount,min_count,ib,ip,ipar,istat,wmonrs,index,innov,iter,ldebug))  then
                  write(*,*) 'using lower level failed'
                  obsspagarr(ib,ip,ipar,istat)=rcpara%miss_val
                endif
              else
                obsspagarr(ib,ip,ipar,istat)=rcpara%miss_val
              endif
            endif

if(ip .eq. 3 .and. ipar .eq. 1 .and. wmonrs(index(1)) .eq. 52818 .and. innov.eq. 4 .and. ib .eq. 13) then
  write(*,*) 'at the place'
endif

return
end subroutine estimate




subroutine modify_tbindex(rcpara,wmonrs,index,istat,tbindex,tbi,lasts,tb,iter,meta_s,alarm,manual,ldebug) !


implicit none

type(rasocor_namelist),intent(in) :: rcpara
type(metadata) :: meta_s

integer istat,l,ib,i,iter

integer index(rcpara%statmax),statnr
integer :: goodsondes(33),istart
integer :: bindex(rcpara%brmax),bindexsave(rcpara%brmax),tbindex(rcpara%brmax),tbindexsave(rcpara%brmax),bi,tbi,tbisave,ipmax,iparmax,ip,ipar,rad,rad2
logical adj
integer,intent(in) :: wmonrs(rcpara%statmax)
integer,intent(in)         :: lasts(rcpara%brmax,rcpara%parmax)
integer,intent(in) :: alarm(10000,3)
character*6 cstatnr
logical,intent(in) :: tb,ldebug,manual



      statnr=wmonrs(index(istat))
      if (statnr.eq. 47678) then
          write(*,*) 'test'
      endif
      write(cstatnr,'(I6.6)') statnr

      do i=1,rcpara%nmax
        if(meta_s.cardsmeta_s(i,index(istat)) .eq. 1000) then
          call addbreak(rcpara,cstatnr,.true.,i,&
                        tbindex,tbisave,tbi,ldebug)
write(*,*) 'jump',cstatnr,todate(i,rcpara)
        endif
      enddo
      if(.false.) then
        rad=toindex(19991001,rcpara)
       adj=.false.
       do ib=1,tbi
         if(abs(tbindex(ib)-rad) .lt. rcpara%snht_maxlen/8) adj=.true.
       enddo
       if(adj) then
        i=1
        do while(tbindex(i) .lt. rad .and. i .le. tbi+1)
          i=i+1
        enddo
        tbindex(i:tbi-1)=tbindex(i+1:tbi)
        tbi=tbi-1
         if(ldebug) write(*,*) cstatnr,' break', rad,' removed from tbindex'
       endif
      endif

      call delbreak(rcpara,cstatnr,statnr .eq. 94998, &
                    (/toindex(19800101,rcpara),toindex(19850101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)

      call delbreak(rcpara,cstatnr,statnr .eq. 81405, &
                    (/toindex(19790101,rcpara),toindex(19820101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)

!      call delbreak(rcpara,cstatnr,.true.,(/18991,18993/)&
!                    ,tbindex,tbisave,tbi,ldebug)

if(rcpara%smooth_method .eq. 1) then 


      call addbreak(rcpara,cstatnr,statnr .eq. 63741, &
                    toindex(19740101,rcpara),&
                     tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 91765, &
                    (/toindex(19800101,rcpara),toindex(19900101,rcpara)/),&
                     tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .gt. 91300 .and. statnr .lt. 91500, &
                    (/toindex(19610101,rcpara),toindex(19710101,rcpara)/),&
                     tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 48698 , &
                    (/toindex(19610101,rcpara),toindex(19710101,rcpara)/),&
                     tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 60680 , &
                    (/toindex(19670101,rcpara),toindex(19690101,rcpara)/),&
                     tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 87623, &
                    (/toindex(20060901,rcpara),toindex(20070101,rcpara)/),tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .gt. 83700 .and. statnr .lt. 85800, &
                    (/toindex(19830101,rcpara),toindex(19870101,rcpara)/),tbindex,tbisave,tbi,ldebug)
!      call delbreak(rcpara,cstatnr,statnr .gt. 91000 .and. statnr .lt. 91500, &
!                    (/toindex(19620101,rcpara),toindex(19850101,rcpara)/),tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 85543, &
                    (/toindex(19620101,rcpara),toindex(19870101,rcpara)/),tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .gt. 83700 .and. statnr .lt. 83800, &
                    (/toindex(19890101,rcpara),toindex(19950101,rcpara)/),tbindex,tbisave,tbi,ldebug)
      call addbreak(rcpara,cstatnr,statnr .gt. 83000 .and. statnr .lt. 84800, &
                    toindex(19890101,rcpara),tbindex,tbisave,tbi,ldebug)

      call delbreak(rcpara,cstatnr,statnr .gt. 24100 .and. statnr .lt. 24400, &
                    (/toindex(19950101,rcpara),toindex(19990101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 29282, &
                    (/toindex(20030101,rcpara),toindex(20040101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 4320, &
                    (/toindex(19830101,rcpara),toindex(19840101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 91285, &
                    (/toindex(19980101,rcpara),toindex(20000101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 91165, &
                    (/toindex(19880101,rcpara),toindex(19900101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 89022, &
                    (/toindex(19960101,rcpara),toindex(19980101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 89022, &
                    (/toindex(20010101,rcpara),toindex(20030101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call addbreak(rcpara,cstatnr,statnr .eq. 61998, &
                    toindex(19950101,rcpara),tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 61998, &
                    (/toindex(19810101,rcpara),toindex(19910101,rcpara)/),tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 63985, &
                    (/toindex(19870101,rcpara),toindex(19950101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 63985, &
                    (/toindex(19980101,rcpara),toindex(20050101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 89532, &
                    (/toindex(19940101,rcpara),toindex(20030101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 89532, &
                    (/toindex(19880101,rcpara),toindex(19900101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call addbreak(rcpara,cstatnr,statnr .gt. 94000 .and. statnr .lt. 95000 &
                    .or. any(statnr .eq.(/91517,91680,96996,96441,96413,89611,89642,89571/)) , &
                    toindex(19880101,rcpara),tbindex,tbisave,tbi,ldebug)
      call addbreak(rcpara,cstatnr,statnr .eq. 91517 , &
                    toindex(19740601,rcpara),tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 91643, &
                    (/toindex(19930101,rcpara),toindex(19970101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 91643, &
                    (/toindex(19800101,rcpara),toindex(19820101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 89664, &
                    (/toindex(20020101,rcpara),toindex(20060101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 89664, &
                    (/toindex(19860101,rcpara),toindex(19960101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 89592, &
                    (/toindex(19850101,rcpara),toindex(19890101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 94975, &
                    (/toindex(19880101,rcpara),toindex(19990101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 94294, &
                    (/toindex(19850101,rcpara),toindex(19860301,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 81405, &
                    (/toindex(19820101,rcpara),toindex(19830101,rcpara)/),tbindex,tbisave,tbi,ldebug)
      call addbreak(rcpara,cstatnr,statnr .gt. 87000 .and. statnr .lt. 89000, &
                    toindex(19870101,rcpara),tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 71082, &
                    (/toindex(19930101,rcpara),toindex(19950101,rcpara)/),tbindex,tbisave,tbi,ldebug)
       call delbreak(rcpara,cstatnr,statnr .gt. 68000 .and. statnr .lt. 69000, &
                    (/toindex(19850601,rcpara),toindex(19900101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 68842, &
                    (/toindex(19970101,rcpara),toindex(19990101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
!      call delbreak(rcpara,cstatnr,statnr .eq. 68842, &
!                    (/toindex(19730101,rcpara),toindex(19740101,rcpara)/)&
!                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .gt. 68000 .and. statnr .lt. 69000, &
                    (/toindex(19920101,rcpara),toindex(19960101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .gt. 68000 .and. statnr .lt. 69000, &
                    (/toindex(19670101,rcpara),toindex(19760101,rcpara)/)&
                   ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 68994, &
                    (/toindex(19960101,rcpara),toindex(20060101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 68906, &
                    (/toindex(19820101,rcpara),toindex(19840101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 26063, &
                    (/toindex(19690101,rcpara),toindex(19700101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call addbreak(rcpara,cstatnr,statnr .eq. 26063, toindex(19680601,rcpara)&
                    ,tbindex,tbisave,tbi,ldebug)
!!$      call delbreak(rcpara,cstatnr,statnr .eq. 61902, &
!!$                    (/toindex(19860101,rcpara),toindex(19870101,rcpara)/)&
!!$                    ,tbindex,tbisave,tbi,ldebug)
!!$      call delbreak(rcpara,cstatnr,statnr .eq. 61902, &
!!$                    (/toindex(19860101,rcpara),toindex(19870101,rcpara)/)&
!!$                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 61902, &
                    (/toindex(19690101,rcpara),toindex(19710101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 82599, &
                    (/toindex(19780101,rcpara),toindex(19790101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call addbreak(rcpara,cstatnr,statnr .eq. 82599, &
                    toindex(19770101,rcpara),tbindex,tbisave,tbi,ldebug)

      call delbreak(rcpara,cstatnr,statnr .eq. 91217, &
                    (/toindex(19900101,rcpara),toindex(19920101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 91366 .or. statnr .eq. 91376, &
                    (/toindex(19850101,rcpara),toindex(19880101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 96471, &
                    (/toindex(19860101,rcpara),toindex(19880101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 70308, &
                    (/toindex(19880101,rcpara),toindex(19900101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 61996, &
                    (/toindex(19800101,rcpara),toindex(19840101,rcpara)/),tbindex,tbisave,tbi,ldebug)
!      call delbreak(rcpara,cstatnr,statnr .eq. 17280, &
!                    (/toindex(19840101,rcpara),toindex(19850101,rcpara)/),tbindex,tbisave,tbi,ldebug)
      call addbreak(rcpara,cstatnr,statnr .eq. 17030, &
                    toindex(19840501,rcpara),tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 91925, &
                    (/toindex(19820101,rcpara),toindex(19840101,rcpara)/),tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 91925, &
                    (/toindex(19900101,rcpara),toindex(19920101,rcpara)/),tbindex,tbisave,tbi,ldebug)
      call addbreak(rcpara,cstatnr,statnr .eq. 91925, &
                    toindex(19920801,rcpara),tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 91938, &
                    (/toindex(20020101,rcpara),toindex(20110101,rcpara)/),tbindex,tbisave,tbi,ldebug)
!      call delbreak(rcpara,cstatnr,statnr .eq. 91938, &
!                    (/toindex(19810101,rcpara),toindex(19830101,rcpara)/),tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 91938, &
                    (/toindex(19700101,rcpara),toindex(19720101,rcpara)/),tbindex,tbisave,tbi,ldebug)
!      call delbreak(rcpara,cstatnr,statnr .eq. 91938, &
!                    (/toindex(19620101,rcpara),toindex(19640101,rcpara)/),tbindex,tbisave,tbi,ldebug)

      call addbreak(rcpara,cstatnr,statnr .eq. 91938, &
                    toindex(19701101,rcpara),tbindex,tbisave,tbi,ldebug)
      call addbreak(rcpara,cstatnr,statnr .eq. 91938, &
                    toindex(19720301,rcpara),tbindex,tbisave,tbi,ldebug)

!      call delbreak(rcpara,cstatnr,statnr .gt. 91000 .and. statnr .lt. 92000, &
!                    (/toindex(20090101,rcpara),toindex(20120101,rcpara)/),tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 97180, &
                    (/toindex(19900101,rcpara),toindex(19970101,rcpara)/),&
                     tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 80222, &
                    (/toindex(19800101,rcpara),toindex(19920101,rcpara)/),&
                     tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 84628, &
                    (/toindex(19580101,rcpara),toindex(19690101,rcpara)/),&
                     tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 85442, &
                    (/toindex(19580101,rcpara),toindex(19800101,rcpara)/),&
                     tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 8579, &
                    (/toindex(19740101,rcpara),toindex(19800101,rcpara)/),&
                     tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 8579, &
                    (/toindex(19670101,rcpara),toindex(19680101,rcpara)/),&
                     tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 83746, &
                    (/toindex(19740101,rcpara),toindex(19750101,rcpara)/),&
                     tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 56146, &
                    (/toindex(20020101,rcpara),toindex(20070101,rcpara)/),&
                     tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 56146, &
                    (/toindex(19830101,rcpara),toindex(19840101,rcpara)/),&
                     tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 20744, &
                    (/toindex(20030101,rcpara),toindex(20060101,rcpara)/),&
                     tbindex,tbisave,tbi,ldebug)
!      call addbreak(rcpara,cstatnr,statnr .eq. 47058, &
!                    toindex(19940101,rcpara)&
!                    ,tbindex,tbisave,tbi,ldebug)


      l=1
      do while (alarm(l,1) .ne. 0)

        if(alarm(l,1) .eq. statnr) then  
!          if(alarm(l,3) .eq. 0) then 
             call delbreak(rcpara,cstatnr,statnr .eq. statnr, &
                    (/alarm(l,2)-rcpara%snht_maxlen/3,alarm(l,2)+rcpara%snht_maxlen/3/) ,tbindex,tbisave,tbi,ldebug)
             call addbreak(rcpara,cstatnr,statnr .eq. statnr, &
                    alarm(l,2) ,tbindex,tbisave,tbi,ldebug)
          exit
        endif
        l=l+1
      enddo

endif

if(.false.) then  
      call delbreak(rcpara,cstatnr,statnr .gt. 78000 .and. statnr .lt. 79000, &
                    (/toindex(19850101,rcpara),toindex(19900101,rcpara)/),tbindex,tbisave,tbi,ldebug)
!      call addbreak(rcpara,cstatnr,statnr .eq. 61998, &
!                    toindex(19890101,rcpara),tbindex,tbisave,tbi,ldebug)

!      call addbreak(rcpara,cstatnr,statnr .eq. 89642, &
!                    toindex(19850701,rcpara),tbindex,tbisave,tbi,ldebug)

!      call addbreak(rcpara,cstatnr,statnr .eq. 63985, &
!                    toindex(19970901,rcpara),tbindex,tbisave,tbi,ldebug)
!      call addbreak(rcpara,cstatnr,statnr .eq. 63985, &
!                    toindex(19880901,rcpara),tbindex,tbisave,tbi,ldebug)
!      call delbreak(rcpara,cstatnr,statnr .eq. 96996, &
!                    (/toindex(19860101,rcpara),toindex(19870101,rcpara)/)&
!                    ,tbindex,tbisave,tbi,ldebug)
!      call delbreak(rcpara,cstatnr,statnr .eq. 91408, &
!                    (/toindex(19990101,rcpara),toindex(20000101,rcpara)/)&
!                    ,tbindex,tbisave,tbi,ldebug)
!      call delbreak(rcpara,cstatnr,statnr .eq. 61641, &
!                    (/toindex(19750101,rcpara),toindex(19820101,rcpara)/)&
!                    ,tbindex,tbisave,tbi,ldebug)

      call delbreak(rcpara,cstatnr,statnr .eq. 91217, &
                    (/toindex(19600101,rcpara),toindex(19710101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 91212, &
                    (/toindex(19600101,rcpara),toindex(19710101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 91212, &
                    (/toindex(19730101,rcpara),toindex(19900101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 91217, &
                    (/toindex(19970101,rcpara),toindex(20000101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 91212, &
                    (/toindex(19970101,rcpara),toindex(20000101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
!      call delbreak(rcpara,cstatnr,statnr .eq. 91958, &
!                    (/toindex(19880101,rcpara),toindex(19890101,rcpara)/)&
!                    ,tbindex,tbisave,tbi,ldebug)
!      call delbreak(rcpara,cstatnr,statnr .eq. 83476, &
!                    (/toindex(19890101,rcpara),toindex(19920101,rcpara)/)&
!                    ,tbindex,tbisave,tbi,ldebug)
!      call delbreak(rcpara,cstatnr,statnr .eq. 61291, &
!                    (/toindex(19810101,rcpara),toindex(19830101,rcpara)/)&
!                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 83971, &
                    (/toindex(19990101,rcpara),toindex(20000101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)
      call delbreak(rcpara,cstatnr,statnr .eq. 83971, &
                    (/toindex(19890101,rcpara),toindex(19900101,rcpara)/)&
                    ,tbindex,tbisave,tbi,ldebug)


!      endif !.false.

endif

     tbindexsave=tbindex
     tbisave=tbi

     tbindex=0
     l=0
     do ib=1,tbi
       if (.not. any(tbindexsave(ib)-rcpara%snht_maxlen/2 .eq. lasts)) then
         if(l .eq. 0) then
           l=l+1
           tbindex(l)=tbindexsave(ib)
         else
!          if(tbindexsave(ib) .gt. tbindex(l)+rcpara%snht_maxlen/3) then
            l=l+1
            tbindex(l)=tbindexsave(ib)
!          else
!            if(ldebug) write(*,*) cstatnr,' break', todate(tbindexsave(ib),rcpara),' removed from tbindex, too close'
!          endif
         endif
       else
         if(ldebug) write(*,*) cstatnr,' initial adjustment', todate(tbindexsave(ib),rcpara),' removed from tbindex'
       endif
     enddo
     tbi=l
!
     tbindex(tbi+1)=rcpara%nmax
! create composite from nearest cmax stations
! index(istat) is station itself, therefore omitted

   end subroutine modify_tbindex

subroutine addbreak(rcpara,cstatnr,logic,rad,tbindex,tbisave,tbi,ldebug)

implicit none
type(rasocor_namelist),intent(in) :: rcpara

integer tbi,tbisave,tbindex(:),ib,i,rad
logical ldebug,logic,adj
character*6 cstatnr


      if(logic) then
       adj=.false.
       do ib=1,tbi
         if(abs(tbindex(ib)-rad) .lt. rcpara%snht_maxlen/4) adj=.true.
       enddo
       if(.not. adj) then
        i=1
        do while(tbindex(i) .lt. rad .and. i .le. tbi+1)
          i=i+1
        enddo
        tbindex(i+1:tbi+2)=tbindex(i:tbi+1)
        tbindex(i)=rad
        tbi=tbi+1
         if(ldebug) write(*,*) cstatnr,' break',  todate(tbindex(i),rcpara),rad,' added to tbindex'
       endif
         
      endif
end subroutine addbreak

subroutine delbreak(rcpara,cstatnr,logic,range,tbindex,tbisave,tbi,ldebug)

implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer tbi,tbisave,tbindex(tbi+1),i,range(2)
logical ldebug,logic
character*6 cstatnr

      if(logic) then
         
         tbisave=tbi
         i=1
         do while(i .le. tbi)
           if(tbindex(i) .gt. range(1) .and. tbindex(i) .lt. range(2)) then
             if(ldebug) write(*,*) cstatnr,' break', todate(tbindex(i),rcpara),tbindex(i),' removed from tbindex'
             tbindex(i:tbi)=tbindex(i+1:tbi+1)
             tbi=tbi-1
           else
             i=i+1
           endif
         enddo
      endif
end subroutine delbreak


end module rfcomp_1
