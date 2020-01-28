module correct_mr

use rfcor

contains

!! correct most recent part of series
!! do this after breakpoint correction since then a longer interval for more accurate
!! estimation can be used.
subroutine correct_mostrecent(rcpara,iname,midx,lasts,gcount,tfgm,stfgm,tgps,needs_composite,ini_correct,rasocorrs,rasobreaks,rasobreakuncertainties) !

implicit none

type(rasocor_namelist),intent(in) :: rcpara
type(rasocor_namelist) :: rcpara1

integer protunit,istart,istop,its,iunit,i,ib,iname,ilayer,l,idiv,j,ll,ii
integer ipar,ip,left_maxlen,right_maxlen,mindex,indext799


real(kind=JPRM) :: tfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax),stfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tgps(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: plus(20,rcpara%pmax,rcpara%parmax),minus(20,rcpara%pmax,rcpara%parmax),prms(20,rcpara%pmax,rcpara%parmax),mrms(20,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: plus2(20,rcpara%pmax,rcpara%parmax),minus2(20,rcpara%pmax,rcpara%parmax),prms2(20,rcpara%pmax,rcpara%parmax),mrms2(20,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: minussave(20,rcpara%pmax,rcpara%parmax),mrmssave(20,rcpara%pmax,rcpara%parmax)
real(kind=JPRM),intent(out) :: rasocorrs(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) rasobreaks(rcpara%nmax,rcpara%pmax,rcpara%parmax),rasobreakuncertainties(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: splus(20,rcpara%pmax,rcpara%parmax),sminus(20,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: breakprofile(rcpara%pmax,rcpara%parmax),breakprofile2(rcpara%pmax,rcpara%parmax),breakprofile3(rcpara%pmax,rcpara%parmax),f_val(rcpara%parmax)
real(kind=JPRM) :: dfak,bdev(rcpara%pmax),msave
integer         :: lasts(rcpara%brmax,rcpara%parmax),midx(rcpara%brmax,rcpara%parmax),gcount(rcpara%parmax),counts(rcpara%brmax,rcpara%pmax,rcpara%parmax),last1(20,3),swap(3),isave
integer(kind=JPRM) :: mcountmean(rcpara%nmax,rcpara%parmax,rcpara%parmax),icsave(20,rcpara%pmax,rcpara%parmax),ic(20,rcpara%pmax,rcpara%parmax),icsave2(20,rcpara%pmax,rcpara%parmax),ic2(20,rcpara%pmax,rcpara%parmax),needs_composite(:,:)
real(kind=JPRM) :: msig(20,rcpara%pmax,rcpara%parmax),psig(20,rcpara%pmax,rcpara%parmax),msigsave(20,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: msig2(20,rcpara%pmax,rcpara%parmax),psig2(20,rcpara%pmax,rcpara%parmax),msigsave2(20,rcpara%pmax,rcpara%parmax)
logical :: mask(rcpara%nmax),in_gap,ini_correct(rcpara%parmax),lprof3(20,rcpara%parmax),ldebug

ldebug=.false.

rcpara1=rcpara
rcpara1%nmax=20
rcpara1%smooth_method=0
msig=0.
mrms=0.
minussave=rcpara%miss_val
minus=rcpara%miss_val
minus2=rcpara%miss_val
ic=0
ic2=0
 
!!indext799=toindex(20060131,rcpara)

  lprof3=.false.
  ll=1
  do ipar=1,rcpara%parmax
    do l=1,gcount(ipar)
      lprof3(l,ipar)=lasts(l,ipar) .lt. rcpara%old .and. lasts(1,3-ipar) .ne. 1 .and. (.not. any(abs(lasts(l,ipar)-lasts(1:gcount(3-ipar),3-ipar)) .lt. rcpara%mean_maxlen/2))
      if(lasts(l,ipar) .gt. 1) then 
        last1(ll,:)=(/l,ipar,lasts(l,ipar)/)
        ll=ll+1
      endif
    enddo
  enddo
  lprof3=.false.
  
ll=ll-1
!! bubble sort last1
do i=1,ll
  do j=i,ll
    if(last1(j,3) .gt. last1(i,3)) then
      swap=last1(j,:)
      last1(j,:)=last1(i,:)
      last1(i,:)=swap
    endif
  enddo
enddo

  !!$ call omp_set_lock(omp_lp)
!!  write(*,'(I5,I9,A7,5I5)') iname,todate(lasts(l,ipar),rcpara),' last1 ',last1(1:5,3)
  !!$ call omp_unset_lock(omp_lp)

  do ipar=1,rcpara%parmax
  do i=1,gcount(ipar)
    if(lasts(i,ipar)+rcpara%snht_maxlen/2 .lt. rcpara%old-1) then
      lasts(i,ipar)=lasts(i,ipar)+rcpara%snht_maxlen/2
    else
      lasts(i,ipar)=rcpara%nmax-2
    endif
  enddo
  enddo

  !lasts(1,:)=rcpara%nmax-2 
  lasts(gcount(1)+1,1) =1
  lasts(gcount(2)+1,2) =1
  do ii=1,ll

  ipar=last1(ii,2)
  l=last1(ii,1)
!!  do ipar=1,rcpara%parmax
!! next line is to account for most values "sacrificed" when calculating lasts indetect_gaps
!!      do l=1,gcount(ipar)
!! this is to 

!!    if(lasts(l,ipar) .lt. rcpara%old .or. ini_correct(ipar)) then  !! <2001
        left_maxlen=rcpara%mean_maxlen
        right_maxlen=min(rcpara%nmax-midx(l,ipar),rcpara%mean_maxlen)
!!        right_maxlen=min(indext799-midx(l,ipar),rcpara%mean_maxlen)
      
        if(midx(l,ipar) .lt. rcpara%nmax) then 
  !!$ call omp_set_lock(omp_lp)
          write(*,*) 'midx:',l,ipar,midx(l,ipar)
  !!$ call omp_unset_lock(omp_lp)
          do ip=1,rcpara%pmax
              
              if (midx(l,ipar) .gt. rcpara%old-rcpara%snht_maxlen) then
!              if (last1(l,ipar) .gt. rcpara%old) then
! adjust with Reanalysis state
                 if (ip>9 .or. rcpara%initial_adjust(4:6).ne.'gps') then
                    call calc_gapbreak(rcpara,ip,ipar,midx(l,ipar),minval((/rcpara%nmax-2,last1(ii,3)+rcpara%snht_maxlen/2/)),lasts(l+1,ipar),tfgm,minus(l,ip,ipar),mrms(l,ip,ipar),msig(l,ip,ipar),ic(l,ip,ipar),rcpara%month) !in this file line 850
                    call calc_gapbreak(rcpara,ip,3-ipar,midx(l,ipar),minval((/rcpara%nmax-2,last1(ii,3)+rcpara%snht_maxlen/2/)),last1(ii+1,3)+rcpara%snht_maxlen/2,tfgm,minus2(l,ip,ipar),mrms2(l,ip,ipar),msig2(l,ip,ipar),ic2(l,ip,ipar),rcpara%month)  !in this file line 850
! adjust with GPS-RO state (if tfgm is gps departure)
                 else
                    call calc_gapbreak(rcpara,ip,ipar,midx(l,ipar),minval((/rcpara%nmax-2,last1(ii,3)+rcpara%snht_maxlen/2/)),lasts(l+1,ipar),tgps,minus(l,ip,ipar),mrms(l,ip,ipar),msig(l,ip,ipar),ic(l,ip,ipar),rcpara%month) !in this file line 850
                    if (ic(l,ip,ipar)<rcpara%snht_maxlen/2-rcpara%max_miss/2) then
                       isave=ic(l,ip,ipar)
                       msave=minus(l,ip,ipar)
                       if (msave.ne. rcpara%miss_val) then
                          write(*,*) msave
                       endif
                       call calc_gapbreak(rcpara,ip,ipar,midx(l,ipar),minval((/rcpara%nmax-2,last1(ii,3)+rcpara%snht_maxlen/2/)),lasts(l+1,ipar),tfgm,minus(l,ip,ipar),mrms(l,ip,ipar),msig(l,ip,ipar),ic(l,ip,ipar),rcpara%month)
                
                       write(*,'(A8,4I6,2F8.2)') 'few gps ',midx(l,ipar),ip,isave,ic(l,ip,ipar),msave,minus(l,ip,ipar)
                    endif
                    call calc_gapbreak(rcpara,ip,3-ipar,midx(l,ipar),minval((/rcpara%nmax-2,last1(ii,3)+rcpara%snht_maxlen/2/)),last1(ii+1,3)+rcpara%snht_maxlen/2,tgps,minus2(l,ip,ipar),mrms2(l,ip,ipar),msig2(l,ip,ipar),ic2(l,ip,ipar),rcpara%month)  !in this file line 850
                    if (ic2(l,ip,ipar)<rcpara%snht_maxlen/2-rcpara%max_miss/2) then
                       isave=ic2(l,ip,ipar) 
                       msave=minus2(l,ip,ipar)
                       call calc_gapbreak(rcpara,ip,3-ipar,midx(l,ipar),minval((/rcpara%nmax-2,last1(ii,3)+rcpara%snht_maxlen/2/)),last1(ii+1,3)+rcpara%snht_maxlen/2,tfgm,minus2(l,ip,ipar),mrms2(l,ip,ipar),msig2(l,ip,ipar),ic2(l,ip,ipar),rcpara%month)  !in this file line 850
!                    !minus(l,ip,ipar)= -minus(l,ip,ipar)
                    !minus2(l,ip,ipar)= -minus2(l,ip,ipar)
                       write(*,'(A9,4I6,2F8.2)') 'few gps2 ',midx(l,ipar),ip,isave,ic2(l,ip,ipar),msave,minus2(l,ip,ipar)
                    endif
                 endif
              else
! adjust with background departures from neighbouring series (double differencing approach)
                 call calc_gapbreak(rcpara,ip,ipar,midx(l,ipar),minval((/rcpara%nmax-2,last1(ii,3)+rcpara%snht_maxlen/2/)),lasts(l+1,ipar),stfgm,minus(l,ip,ipar),mrms(l,ip,ipar),msig(l,ip,ipar),ic(l,ip,ipar),rcpara%month) !in this file line 850
                 call calc_gapbreak(rcpara,ip,3-ipar,midx(l,ipar),minval((/rcpara%nmax-2,last1(ii,3)+rcpara%snht_maxlen/2/)),last1(ii+1,3)+rcpara%snht_maxlen/2,stfgm,minus2(l,ip,ipar),mrms2(l,ip,ipar),msig2(l,ip,ipar),ic2(l,ip,ipar),rcpara%month)  !in this file line 850
              endif
  !!$ call omp_set_lock(omp_lp)
                if(ldebug) write(*,'(I6,2I9,2I3,7F8.2,2I5,A30)') iname,todate(lasts(l,ipar),rcpara),todate(midx(l,ipar),rcpara), &
ip,ipar,minussave(l,ip,ipar),mrmssave(l,ip,ipar),msigsave(l,ip,ipar),minus(l,ip,ipar),minus2(l,ip,ipar),mrms(l,ip,ipar),msig(l,ip,ipar),ic(l,ip,ipar),ic2(l,ip,ipar),' initial adjustment calculated'
  !!$ call omp_unset_lock(omp_lp)
         enddo
!!      write(*,*) cname,ipar,l,midx(l,ipar),ini_correct,' correcting series between' ,lasts(l+1,ipar),lasts(l,ipar)
       
      splus(l,:,:)=0.
      sminus(l,:,:)=minus(l,:,:)
      breakprofile=rcpara%miss_val
      call calc_profile(rcpara1,iname,l,splus,sminus,psig,msig,breakprofile,f_val,rcpara1%smooth_method,rcpara%parmax) !in file rfcor.f90 line 9
      
      splus(l,:,:)=0.
      sminus(l,:,:)=minus2(l,:,:)
      breakprofile2=rcpara%miss_val
      call calc_profile(rcpara1,iname,l,splus,sminus,psig,msig,breakprofile2,f_val,rcpara1%smooth_method,rcpara%parmax) !in file rfcor.f90 line 9
      
      breakprofile3=rcpara%miss_val
      if(ii .ne. 1 .and. lprof3(l,ipar)) then
        where(breakprofile .ne. rcpara%miss_val .and. breakprofile2 .ne. rcpara%miss_val)
          breakprofile3=(breakprofile-breakprofile2)
        endwhere
        if(.not. any(breakprofile3 .ne. rcpara%miss_val)) breakprofile3=breakprofile
      else
        breakprofile3=breakprofile
      endif
      
!!      do ll=l+1,gcount(ipar)
!!        where(breakprofile3 .ne. rcpara%miss_val) minus(l+1:gcount(ipar),ipar,:)=breakprofile3(:,ipar)
!!      enddo
      do ip=1,rcpara%pmax
        if(breakprofile3(ip,ipar) .ne. rcpara%miss_val .and. (ini_correct(ipar) .and. ip.lt.14 )) then !.or. lasts(l,ipar).lt. rcpara%old)) then
     
!          do i=lasts(l+1,ipar),lasts(l,ipar)-1
          do i=1,lasts(l,ipar)-1
!!             rasocorrs(i,ip,ipar)=rasocorrs(i,ip,ipar)-rasocorrs(lasts(l,ipar)-1,ip,ipar)+(1.0-rcpara%plevs(ip)/850.)*breakprofile3(ip,ipar)
!             rasocorrs(i,ip,ipar)=rasocorrs(i,ip,ipar)-rasocorrs(lasts(l,ipar)-1,ip,ipar)+breakprofile3(ip,ipar)
             rasocorrs(i,ip,ipar)=rasocorrs(i,ip,ipar)+breakprofile3(ip,ipar)
             if(tfgm(i,ip,ipar) .ne. rcpara%miss_val) tfgm(i,ip,ipar)=tfgm(i,ip,ipar)+breakprofile3(ip,ipar)
             if(stfgm(i,ip,ipar) .ne. rcpara%miss_val) stfgm(i,ip,ipar)=stfgm(i,ip,ipar)+breakprofile3(ip,ipar)
             if(tgps(i,ip,ipar) .ne. rcpara%miss_val) tgps(i,ip,ipar)=tgps(i,ip,ipar)+breakprofile3(ip,ipar)
          enddo
!!$ call omp_set_lock(omp_lp)
          if(ldebug) then
             write(*,*) iname,ip,ipar,'lastsraso',todate(lasts(l+1,ipar),rcpara),todate(lasts(l,ipar)-1,rcpara),rasocorrs(lasts(l,ipar)-1,ip,ipar),' has been adjusted'
          endif
!!$ call omp_unset_lock(omp_lp)
          if(minus(l,ip,ipar) .ne. rcpara%miss_val ) then
            do i=lasts(l+1,ipar),lasts(l,ipar)-1
               rasobreaks(i,ip,ipar)=rasobreaks(i,ip,ipar)-rasobreaks(lasts(l,ipar),ip,ipar)-minus(l,ip,ipar)
               rasobreakuncertainties(i,ip,ipar)=rasobreakuncertainties(i,ip,ipar)+msig(l,ip,ipar)
            enddo
          endif
        else
!!$ call omp_set_lock(omp_lp)
            if (ini_correct(ipar) .or. ip.lt.14 .or. lasts(l,ipar).lt. rcpara%old) then 
            write(*,'(I6,2I3,A9,I9,I9,A30,3L3)') iname,ip,ipar,'lastsraso',todate(lasts(l+1,ipar),rcpara),todate(lasts(l,ipar)-1,rcpara),'has not been adjusted',&
                     ini_correct(ipar), ip.lt.14, lasts(l,ipar).lt. rcpara%old
            endif
!!$ call omp_unset_lock(omp_lp)
        endif
      enddo
      needs_composite(l,ipar)=rcpara%nmax-2
     else
!!$ call omp_set_lock(omp_lp)
       write(*,*) iname,ini_correct(ipar),lasts(l,ipar),' not initial adjusted'
!!$ call omp_unset_lock(omp_lp)
     endif
!!    if(any(breakprofile(:,ipar) .ne. rcpara%miss_val))  chosenbreaks(1)=rcpara%nmax-2 
!!     enddo
!!     endif
     ini_correct(ipar)=ini_correct(ipar) .and. any(needs_composite(1:gcount(ipar),ipar)<rcpara%old .and. needs_composite(1:gcount(ipar),ipar)>1)
  enddo
  plus(1,:,:)=rcpara%miss_val
  minus(1,:,:)=rcpara%miss_val

return
end subroutine correct_mostrecent

subroutine calc_gapbreak(rcpara,ip,ipar,midx,lasts,lasts1,tfgm,minus,mrms,msig,ic,bin) !

implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer midx,lasts,lasts1,i,ip,ipar,ic
real(kind=JPRM) :: tfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: minus,mrms,msig
!integer,save,allocatable :: bin(:)
integer :: bin(:)
real(kind=JPRM) :: msum(12),rmssum(12),fak
integer         :: mcount(12),mc
logical ldebug

ldebug=.false.
!!$if( .not. allocated(bin)) then
!!$  allocate(bin(rcpara%nmax))
!!$  do i=1,rcpara%nmax
!!$    bin(i)=modulo(i*1.,365.)
!!$    bin(i)=floor(bin(i)/365.*12+1)
!!$  enddo
!!$endif

  minus=0.
  mrms=0.
  msum=0.
  rmssum=0.
  mcount=0
  ic=0
  fak=2.0
  do i=midx,lasts
    if(tfgm(i,ip,ipar) .ne. rcpara%miss_val) then
      ic=ic+1
      msum(bin(i))=msum(bin(i))+tfgm(i,ip,ipar)
      mcount(bin(i))=mcount(bin(i))+1
      minus=minus+tfgm(i,ip,ipar)
      rmssum(bin(i))=rmssum(bin(i))+tfgm(i,ip,ipar)*tfgm(i,ip,ipar)
    endif
  enddo
  if(ic .le. rcpara%snht_maxlen/fak) then
    i=midx
!!    do while(i .gt. midx-rcpara%mean_maxlen/2. .and. ic .le. rcpara%snht_maxlen/2-rcpara%max_miss .and. i .gt. lasts1 .and. i .gt. 1)
    do while(ic .le. rcpara%snht_maxlen/fak .and. i .gt. lasts1 .and. i .gt. 1)
      i=i-1
      if(tfgm(i,ip,ipar) .ne. rcpara%miss_val) then
        ic=ic+1
        msum(bin(i))=msum(bin(i))+tfgm(i,ip,ipar)
        mcount(bin(i))=mcount(bin(i))+1
        minus=minus+tfgm(i,ip,ipar)
        rmssum(bin(i))=rmssum(bin(i))+tfgm(i,ip,ipar)*tfgm(i,ip,ipar)
      endif
    enddo
  endif

!!$ call omp_set_lock(omp_lp)
!!$ call omp_unset_lock(omp_lp)
  where(mcount .gt. 2)
    msum=msum/mcount
    rmssum=rmssum/mcount
  elsewhere
    msum=0.
    rmssum=0.
  endwhere
  mc=count(mcount .gt. 2)
  if(mc .gt. 1 .and. ic .gt. rcpara%snht_maxlen/fak-rcpara%max_miss) then
    minus=sum(msum)/mc
    mrms=sqrt(sum(rmssum)/mc)
    msig=sqrt(mrms*mrms-minus*minus)/sqrt(ic*1.0)
  else
    minus=rcpara%miss_val
    mrms=rcpara%miss_val
    msig=rcpara%miss_val
  endif
  if(ldebug)  write(*,'(A3,I6,I3,I2,I6,I6,I5,F4.0,I5,F9.4)') 'IC:',i,ip,ipar,midx,lasts1,ic,rcpara%snht_maxlen/fak-rcpara%max_miss,sum(mcount),minus


end subroutine calc_gapbreak

end module correct_mr
