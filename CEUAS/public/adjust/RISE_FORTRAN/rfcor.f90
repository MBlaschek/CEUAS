module rfcor
  use txtnc
  use rfmod
  use homtests
  use homtestsamp
  use rfcorio
  use sort
  use data_setting

  interface
     subroutine calc_profile(rcpara,iname,breakloc,plus,minus,rplus,rminus,breakprofile,f_val,kplus1,dim3) !
       use rfmod

       implicit none

       type(rasocor_namelist),intent(in) :: rcpara

       integer(kind=JPRM)        :: kplus1,dim3
       integer(kind=JPRM) :: breakloc,iname
       real(kind=JPRM),intent(in) :: plus(rcpara%nmax,rcpara%pmax,dim3),minus(rcpara%nmax,rcpara%pmax,dim3)
       real(kind=JPRM),intent(in) :: rplus(rcpara%nmax,rcpara%pmax,dim3),rminus(rcpara%nmax,rcpara%pmax,dim3)
       real(kind=JPRM) :: breakprofile(rcpara%pmax,dim3),f_val(dim3)
     end subroutine calc_profile
  end interface

contains

  subroutine  check_ini(meta_s,istat,wmonrs,rcpara,ini_correct)
    use rfmod

    implicit none

    type(rasocor_namelist),intent(in) :: rcpara
    type(metadata) :: meta_s
    integer :: goodsondes(28),istat,statnr,i
    integer :: wmonrs(:)
    logical :: ini_correct(rcpara%parmax,rcpara%statmax),ldebug

    !goodsondes=(/37,52,53,57,60,61,62,63,66,67,71,72,73,74,78,79,80,81,82,83,47,55,56,26,76,85,86,87/)
    !ini_correct(:,istat)=.not. any(imax .eq. goodsondes)
    ldebug=.False.
    statnr=wmonrs(istat)
    if(statnr .gt. 62000 .and. statnr .lt. 63000) ini_correct(:,istat)=.true.
    if(statnr .eq. 15614) ini_correct(:,istat)=.true.
    if(statnr .gt. 68000 .and. statnr .lt. 69000) ini_correct(:,istat)=.false.
    if(statnr .eq. 22550) ini_correct(:,istat)=.true.
    do i=rcpara%nmax,toindex(19750101,rcpara),-1

       if(meta_s.cardsmeta_s(i,istat) .gt. 0) then
          if(meta_s.trusted(meta_s.cardsmeta_s(i,istat))) then
             ini_correct(:,istat)=.false.
             if(ldebug) write(*,*) 'Trusted: ',wmonrs(istat),i,todate(i,rcpara),meta_s.rscodes(meta_s.cardsmeta_s(i,istat))
             exit
          endif
       endif
    enddo
    if(statnr .eq. 48407 .or. statnr .eq. 48568 .or. statnr .eq. 48565 .or. statnr .eq. 96315 &
         .or. statnr .eq. 85586 .or. statnr .eq. 48900 .or. statnr .eq. 78384 .or. statnr .eq. 78762 &
         .or. statnr .eq. 7645.or. statnr .eq.  85442 .or. statnr .eq. 85799 .or. statnr .eq. 44373) ini_correct(:,istat)=.true.
  end subroutine check_ini

subroutine load_richcorr(tccr,icistat,ipar,wmonrs,idx,rcpara, mrasocorrs,flogcache,rbindex,rbi)

type(rasocor_namelist),intent(in) :: rcpara
type(cacherecord) :: tccr(rcpara%cachemax,3)

logical :: logcache
logical,optional :: flogcache
integer :: icistat,ipar,idx,err
integer,intent(in) :: wmonrs(rcpara%statmax)
real(kind=JPRM) :: mrasocorrs(rcpara%mmax,rcpara%pmax,rcpara%parmax)
integer :: bindex(rcpara%mmax),bi
integer,optional :: rbindex(rcpara%mmax),rbi
character filename*100

!$ call omp_set_lock(omp_lp(icistat))
    if(present(flogcache)) then
      logcache=flogcache
    else
      logcache=tccr(icistat,1)%vals .ne. 0
    endif
!$ call omp_unset_lock(omp_lp(icistat))
      if(logcache) then

!$ call omp_set_lock(omp_lp(icistat))
      bi=tccr(icistat,1)%vals
      bindex(1:bi)=tccr(icistat,1)%index(1:bi)
      if (ipar==0) then
      mrasocorrs(1:bi,:,:)=tccr(icistat,1)%feld(1:bi,:,:)
      else
      mrasocorrs(1:bi,:,ipar)=tccr(icistat,1)%feld(1:bi,:,ipar)
      endif

!$ call omp_unset_lock(omp_lp(icistat))

     else

!$ call omp_set_lock(omp_lp(icistat))
      write(filename,'(I6.6,a,a,a,I6.6,a)') wmonrs(idx),'/feedbackglobbincorrsave_rit',rcpara%ens,'_',wmonrs(idx),'.nc'
      CALL read_sonde_corr_daily_nc(filename, rcpara,idx, err,mrasocorrs, bindex,bi)
      if(err .ne. 0 .or. any(isnan(mrasocorrs))) then
      write(filename,'(I6.6,a,I6.6,a)') wmonrs(idx),'/feedbackglobbincorrsave',wmonrs(idx),'.nc'
        CALL read_sonde_corr_daily_nc(filename, rcpara,idx, err,mrasocorrs, bindex,bi) 
        if(err .eq. 0) then 
          write(*,*) 'read ',filename,' for initial adjustment since RICH-adjusted version was not available or spurious'
        else
          stop 'could not find data for initial adjustment'
        endif
      endif
!$ call omp_unset_lock(omp_lp(icistat))

      if(bi .ge. 2) bindex(2:bi)=bindex(2:bi)-1
      
      if(present(flogcache)) then
        if(.not. flogcache) then
           rbindex=bindex
           rbi=bi
           return
        endif
      endif
!$ call omp_set_lock(omp_lp(icistat))
      tccr(icistat,1)%vals=bi
      tccr(icistat,1)%index=rcpara%nmax
      tccr(icistat,1)%index(1:bi)=bindex(1:bi)
      if(.not. allocated(tccr(icistat,1)%feld)) allocate(tccr(icistat,1)%feld(rcpara%brmax,rcpara%pmax,rcpara%parmax))
      tccr(icistat,1)%feld(1:bi,:,:)=mrasocorrs(1:bi,:,:)
!      if(any(isnan(tccr(icistat,1)%feld))) then
!        stop
!      endif
!$ call omp_unset_lock(omp_lp(icistat))
     endif
   return
    
end subroutine load_richcorr

  subroutine spagstat(plus,minus,pcount,mcount,pms,mms,sig,miss_val) !

    implicit none

    real(kind=JPRM),intent(in) :: plus,minus,pms,mms
    integer,intent(in) :: pcount,mcount
    real(kind=JPRM) :: psig,msig,sig,miss_val,sigma1,sigma2

    sig=miss_val
    if ((plus .ne. miss_val) .and. minus .ne. miss_val ) then

       sigma1=pms-plus*plus
       sigma2=mms-minus*minus
       if(sigma1 .gt. 0 .and. sigma2 .gt. 0) then
          sig=sqrt(sigma1/pcount+sigma2/mcount)
          !           write(*,'(a,7F8.3,2I6)') 'stddev',sig,sigma1,sigma2,pms,plus,mms,minus,pcount,mcount
          !         else
          !           if(pcount .gt. 600 .and. mcount .gt. 600) then
          !           write(*,'(a,6F8.3,2I6)') 'spurious radicand',sigma1,sigma2,pms,plus,mms,minus,pcount,mcount
          !           endif   
       endif
    endif
  end subroutine spagstat

  subroutine remove_signal(rcpara,cstatnr,tm,rasocorrs,istart,istop,daily,lrandom)

    use ifport,only: rand
    implicit none


    type(rasocor_namelist),intent(in) :: rcpara

    integer i,j,l,m,istart,istop,ip,ipar,err,iter,ib,statnr,tbi,k,iseed
    real(kind=JPRM),intent(in) :: tm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM)  :: rasocorrs(rcpara%nmax,rcpara%pmax,rcpara%parmax),rc2(rcpara%mmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM),allocatable  :: anom(:,:,:),tmmon(:,:,:)
    real(kind=JPRM),allocatable  :: hilf(:),hilf2(:)
    integer,allocatable         :: index(:),index2(:)
    integer :: mistart,mistop,mi,mrcindex(rcpara%mmax)
    integer,target :: rcindex(rcpara%mmax),lindex(4)
    real(kind=JPRM),allocatable  :: climate(:)
    real(kind=8),target,allocatable          :: rrcindex(:)
    real :: a(2,rcpara%pmax,rcpara%parmax),b(2,rcpara%pmax,rcpara%parmax),opt,optold
    character*80    :: filename
    character*6     :: cstatnr
    logical :: daily
    logical,intent(in),optional :: lrandom

    filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbincorrsave'//cstatnr//'.nc'
    read(cstatnr,*) statnr

    if(present(lrandom)) then 
       tbi=5
       do j=1,tbi
          rcindex(j)=istart+floor(rand(0)*(istop-istart))
          if(j .gt. 1 .and. istop-istart .gt. rcpara%snht_maxlen*2) then 
             do while (any(abs(rcindex(1:j-1)-rcindex(j)) .lt. rcpara%snht_maxlen/3))
                rcindex(j)=istart+floor(rand(0)*(istop-istart))
             enddo
          endif

       enddo
       print*,'rcindex:',rcindex(1:tbi)
       l=0
       do j=1,tbi
          if(rcindex(j) .gt. istart+(istop-istart)/8 .and. rcindex(j) .lt. istop-(istop-istart)/8) then
             l=l+1
          endif
       enddo
       if(l .gt. 1) then
          allocate(rrcindex(l))
          rrcindex=rcindex(1:l)
          call qsort(rrcindex)
          rcindex(1:l)=rrcindex
          deallocate(rrcindex)
       endif
    else
       CALL read_sonde_corr_daily_nc(filename, rcpara,0, err,rc2, rcindex,tbi)
       l=0
       do j=1,tbi
          if(rcindex(j) .gt. istart+(istop-istart)/8 .and. rcindex(j) .lt. istop-(istop-istart)/8) then
             l=l+1
          endif
       enddo
    endif
    if(l .eq. 0) then
       ib=(istart+istop)/2
       do j=1,tbi+1
          if(rcindex(j) .gt. ib)  then
             rcindex(j+1:tbi+2)=rcindex(j:tbi+1)   
             rcindex(j)=ib
             exit          
          endif
       enddo
       tbi=tbi+1
    endif

    if(daily) then
       allocate(anom(rcpara%nmax,rcpara%pmax,rcpara%parmax), &
            hilf(rcpara%nmax),hilf2(rcpara%nmax),&
            index(rcpara%nmax),index2(rcpara%nmax),&
            climate(365))  
       anom=tm 
       do ipar=1,rcpara%parmax
          do ip=1,rcpara%pmax
             call dailyanomaly(rcpara,tm(:,ip,ipar),climate,anom(:,ip,ipar),istart,istop)


             l=0
             m=1
             do i=istart,istop
                do j=1,tbi
                   if(rcindex(j) .eq. i) then
                      if(l .gt. 0) m=l
                   endif
                enddo
                if(anom(i,ip,ipar) .ne. rcpara%miss_val) then ! .and. (anom(istop+istart-i,ip,ipar) .ne. rcpara%miss_val .or. anom(istop+istart-i-1,ip,ipar) .ne. rcpara%miss_val)) then
                   l=l+1
                   index(l)=i
                   index2(l)=m
                   hilf(l)=anom(i,ip,ipar)
                   hilf2(l)=0
                endif
             enddo
             if(l .gt. 0) then 
                b(:,ip,ipar)=0.
                call linreg(dble(index(1:l)),dble(hilf(1:l)),a(:,ip,ipar),l)
                ! iterate until difference smaller than 0.01K/decade
                opt=2
                iter=0
                do while(abs(a(2,ip,ipar)-b(2,ip,ipar)) .gt. 0.01/3650 .and. iter .lt. 10)
                   do i=1,l
                      hilf2(i)=a(1,ip,ipar)+opt*a(2,ip,ipar)*index2(i)
                   enddo
                   call linreg(dble(index(1:l)),dble(hilf2(1:l)),b(:,ip,ipar),l)
                   if(b(2,ip,ipar) .ne. 0) then
                      optold=a(2,ip,ipar)/b(2,ip,ipar)
                   else
                      exit
                   endif
                   opt=opt*optold
                   iter=iter+1
                enddo
                do i=2,l
                   rasocorrs(index(i-1):index(i),ip,ipar)=hilf2(i)-hilf2(l)
                enddo
                do j=1,tbi
                   if(rcindex(j) .gt. index(1)) then
                      exit
                   endif
                enddo
                if(j .le. tbi) then 
                   rasocorrs(1:index(1)-1,ip,ipar)=rasocorrs(1:index(1)-1,ip,ipar)-(rasocorrs(index(1)-1,ip,ipar)-rasocorrs(index(1),ip,ipar))
                else
                   rasocorrs(1:index(1)-1,ip,ipar)=rasocorrs(index(1),ip,ipar)
                endif
             endif
          enddo
       enddo

    else ! monthly
       allocate(tmmon(rcpara%mmax,rcpara%pmax,rcpara%parmax), &
            anom(rcpara%mmax,rcpara%pmax,rcpara%parmax), &
            hilf(rcpara%mmax),hilf2(rcpara%mmax),&
            index(rcpara%mmax),index2(rcpara%mmax),&
            climate(12))
       call makemonthly(rcpara,tm,tmmon,5) !
       !      anom=tmmon
       do ipar=1,rcpara%parmax
          do ip=1,rcpara%pmax
             call anomaly(tmmon(:,ip,ipar),anom(:,ip,ipar),climate,rcpara%mmax,12,rcpara%miss_val,15) !in file rfmod.f90 line 2311
             if(.not. any(anom(:,ip,ipar) .ne. rcpara%miss_val)) anom(:,ip,ipar)=tmmon(:,ip,ipar)
             l=0
             m=1
             mistart=(rcpara%year(istart)-1957)*12+rcpara%month(istart)
             mistop=(rcpara%year(istop)-1957)*12+rcpara%month(istop)
             do i=mistart,mistop
                do j=1,tbi
                   mrcindex(j)=(rcpara%year(rcindex(j))-1957)*12+rcpara%month(rcindex(j))
                   if(mrcindex(j) .eq. i) then
                      if(l .gt. 0) m=l
                   endif
                enddo
                if(anom(i,ip,ipar) .ne. rcpara%miss_val   .and. anom(mistop+mistart-i,ip,ipar) .ne. rcpara%miss_val ) then
                   l=l+1
                   index(l)=i
                   index2(l)=m
                   hilf(l)=anom(i,ip,ipar)
                   hilf2(l)=0
                endif
             enddo
             if(l .gt. 0) then 
                call linreg(dble(index(1:l)),dble(hilf(1:l)),a(:,ip,ipar),l)
                ! iterate until difference smaller than 0.01K/decade
                opt=2
                iter=0
                b(:,ip,ipar)=0.
                do while(abs(a(2,ip,ipar)-b(2,ip,ipar)) .gt. 0.01/120 .and. iter .lt. 10)
                   do i=1,l
                      hilf2(i)=a(1,ip,ipar)+opt*a(2,ip,ipar)*index2(i)
                   enddo
                   call linreg(dble(index(1:l)),dble(hilf2(1:l)),b(:,ip,ipar),l)
                   if(b(2,ip,ipar) .ne. 0) then
                      optold=a(2,ip,ipar)/b(2,ip,ipar)
                   else
                      exit
                   endif
                   opt=opt*optold
                   iter=iter+1
                enddo
                ! from here on we again deal with daily data!
                do i=istart,istop
                   do j=1,tbi+1
                      if(rcindex(j) .eq. i) then
                         exit
                      endif
                   enddo
                   if(j .le. tbi) then
                      do k=1,l
                         if(index(k) .ge. mrcindex(j)) then
                            rasocorrs(i,ip,ipar)=hilf2(k)-hilf2(l)
                            exit
                         endif
                      enddo
                   else
                      rasocorrs(i,ip,ipar)=rasocorrs(i-1,ip,ipar)
                   endif
                   if(i .eq. istart) then
                      rasocorrs(i,ip,ipar)=hilf2(1)-hilf2(l)
                   endif
                enddo
                rasocorrs(1:istart-1,ip,ipar)=rasocorrs(1:istart-1,ip,ipar)+rasocorrs(istart,ip,ipar)-rasocorrs(istart-1,ip,ipar)
             endif
          enddo
       enddo
       deallocate(tmmon)
    endif !daily
    deallocate(anom, hilf,hilf2,index,index2,climate)

    return
  end subroutine remove_signal

  subroutine tomonth(rcpara,dindex,daily,monthly)
    implicit none
    type(rasocor_namelist),intent(in) :: rcpara

    integer i,j,l,m,iy,imod,imon,tgm,index
    real(kind=JPRM),intent(in) :: daily(rcpara%nmax)
    integer(kind=JPRM),intent(in) :: dindex(rcpara%nmax)
    real(kind=JPRM) :: monthly(rcpara%mmax)

    do imon=1,rcpara%mmax
       iy=1957+(imon-1)/12
       imod=mod(imon-1,12)+1
       monthly(imon)=0.
       tgm=0
       do while(rcpara%month(index) .eq. imod .and. index .lt. rcpara%nmax)
          if(daily(index) .ne. rcpara%miss_val) then
             monthly(imon)=monthly(imon)+daily(index)
             tgm=tgm+1
          endif
          index=index+1

       enddo

       if(tgm .gt. 0) then
          monthly(imon)=monthly(imon)/tgm
       else
          monthly(imon)=rcpara%miss_val
       endif
    enddo

  end subroutine tomonth

  subroutine linreg(x,y,k,n,weights)
    implicit none

    real :: x(n),y(n),k(2),xq,yq
    real,optional :: weights(n)
    integer:: n

    xq=sum(x)/n
    yq=sum(y)/n
    k(2)=dot_product(x-xq,y-yq)/dot_product(x-xq,x-xq)
    k(1)=yq-k(2)*xq

    return
  end subroutine linreg


  subroutine dailyanomaly(rcpara,original,climate,anom,istart,istop,thresh)

    implicit none

    type(rasocor_namelist),intent(in) :: rcpara

    integer i,j,l,monsave,istart,istop
    integer,optional:: thresh
    real(kind=JPRM),intent(in) :: original(rcpara%nmax)
    real(kind=JPRM)  :: climate(365),anom(rcpara%nmax)

    climate=0.
    anom=rcpara%miss_val

    do i=1,365
       l=0
       do j=istart/365+1,istop/365
          if(original(floor(i+(j-1)*365.25)) .ne. rcpara%miss_val) then
             climate(i)=climate(i)+original(floor(i+(j-1)*365.25))
             l=l+1
          endif
       enddo
       if(present(thresh)) then
          if(l .ge. thresh) then 
             climate(i)=climate(i)/l
          else
             climate(i)=rcpara%miss_val
          endif
       else
          if(l .gt. (istop-istart)/730) then 
             climate(i)=climate(i)/l
          else
             climate(i)=rcpara%miss_val
          endif
       endif
       do j=istart/365+1,istop/365
          if(original(floor(i+(j-1)*365.25)) .ne. rcpara%miss_val .and. climate(i) .ne. rcpara%miss_val) then
             anom(floor(i+(j-1)*365.25))=original(floor(i+(j-1)*365.25))-climate(i)
          endif
       enddo
    enddo
    !    do i=2,rcpara%nmax-1
    !      if(anom(i) .eq. rcpara%miss_val .and. anom(i+1) .ne. rcpara%miss_val .and. anom(i-1) .ne. rcpara%miss_val) then
    !        anom(i)=(anom(i+1)+anom(i-1))/2.
    !      endif
    !    enddo

    return
  end subroutine dailyanomaly

  subroutine dailyanomaly2(rcpara,original,climate,anom,istart,istop,thresh)

    implicit none

    type(rasocor_namelist),intent(in) :: rcpara

    integer i,j,l,monsave,istart,istop,ibin,m,k
    integer,optional:: thresh
    real(kind=JPRM),intent(in) :: original(rcpara%nmax)
    real(kind=JPRM)  :: climate(366),anom(rcpara%nmax)
    integer bin((istop-istart)/365,366)
    climate=0.
    anom=rcpara%miss_val


    l=1
    ibin=0
    bin=0
    do i=istart,istop
       if(rcpara%month(i) .eq. rcpara%month(istart) .and. rcpara%day(i) .eq. rcpara%day(istart)) then
          l=1
          ibin=ibin+1
       endif
       bin(ibin,l)=i 
       l=l+1
    enddo


    do l=1,365
       m=0
       do i=1,ibin
          k=bin(i,l)
          if(k .eq. 0) cycle
          if(original(k) .ne. rcpara%miss_val) then
             climate(l)=climate(l)+original(k)
             m=m+1
          endif
       enddo
       if(m .ge. thresh) then 
          climate(l)=climate(l)/m
          do i=1,ibin
             k=bin(i,l)
             if(k .eq. 0) cycle
             if(original(k) .ne. rcpara%miss_val) anom(k)=original(k)-climate(l) 
          enddo
       else
          climate(l)=rcpara%miss_val
          do i=1,ibin
             anom(bin(i,l))=rcpara%miss_val
          enddo
       endif

    enddo
    !    do i=2,rcpara%nmax-1
    !      if(anom(i) .eq. rcpara%miss_val .and. anom(i+1) .ne. rcpara%miss_val .and. anom(i-1) .ne. rcpara%miss_val) then
    !        anom(i)=(anom(i+1)+anom(i-1))/2.
    !      endif
    !    enddo

    return
  end subroutine dailyanomaly2

  subroutine expand(rcpara,monthly,daily) !

    implicit none

    type(rasocor_namelist),intent(in) :: rcpara

    integer i,l,monsave
    real(kind=JPRM),intent(in) :: monthly(:)
    real(kind=JPRM)  :: daily(:)

    daily=rcpara%miss_val
    l=1
    monsave=rcpara%month(1)
    do i=1,rcpara%nmax
       if(rcpara%month(i) .ne. monsave) then
          l=l+1
          if(l .gt. rcpara%mmax) exit
          monsave=rcpara%month(i)
       endif
       daily(i)=monthly(l)
    enddo

    return
  end subroutine expand

  subroutine detect_gaps(rcpara,iname,tfgm,midx,lasts,gcount) !

    implicit none

    type(rasocor_namelist),intent(in) :: rcpara

    integer mindex,l,ipar,iname
    real(kind=JPRM) :: dfak
    real(kind=JPRM) :: tfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM) :: tsa(rcpara%nmax,rcpara%parmax),plus(rcpara%nmax),minus(rcpara%nmax),prms(rcpara%nmax),mrms(rcpara%nmax),null(rcpara%nmax)
    integer :: philf(rcpara%nmax),mhilf(rcpara%nmax,rcpara%parmax),critical_dates(20)
    !real(kind=JPRM) :: stsa(rcpara%nmax,rcpara%parmax),splus(rcpara%nmax),sminus(rcpara%nmax),sprms(rcpara%nmax),smrms(rcpara%nmax)
    !integer :: sphilf(rcpara%nmax),smhilf(rcpara%nmax,rcpara%parmax)
    integer         :: lasts(rcpara%brmax,rcpara%parmax),midx(rcpara%brmax,rcpara%parmax),gcount(rcpara%parmax),maxlen,max_miss,mean_maxlen
    logical in_gap

    maxlen=rcpara%mean_maxlen
    !!maxlen=rcpara%snht_maxlen
    max_miss=rcpara%mean_maxlen/2-2*(rcpara%snht_maxlen/2-rcpara%max_miss)+10 !! so dass SNHT noch was liefern kann
    !!max_miss=rcpara%max_miss ! so dass SNHT noch was liefern kann
    mhilf=0
    null=0.
    do ipar=1,rcpara%parmax

       call snhteqsamp2(tfgm(:,8,ipar),null,rcpara%nmax,1,rcpara%nmax,maxlen,rcpara%snht_increment,rcpara%miss_val,max_miss,critical_dates,-1, & 
            tsa(:,ipar),plus,minus,prms,mrms,philf,mhilf(:,ipar),rcpara%month) !in file homtests.f90 line 654

       !  call snhteqsamp2(hilf,null,rcpara%nmax,1,rcpara%nmax,rcpara%snht_maxlen,rcpara%snht_increment,rcpara%miss_val,rcpara%max_miss,critical_dates,0, & 
       !    stsa(:,ipar),splus,sminus,sprms,smrms,sphilf,smhilf(:,ipar))

    enddo

    dfak=2.0
    midx=rcpara%nmax
    gcount=0
    do ipar=1,rcpara%parmax

       !! detect data gaps in time series that cannot be bridged by the homogenization algo rithm
       !! In these cases the time series before the gap is treated as a time series that ends
       !! before the gap.
       !! Even if a time series has no gaps, the period after the end of the time series is treated as a gap.
       !! Time series between gaps are treated as independent time series.

       !! minus, plus of the first few days (19570101-19570220) which contain no radiosonde data anyway,
       !! are used to store gap information

       mindex=rcpara%nmax
       in_gap=.true.
       l=0
       lasts(:,ipar)=0
       do while(mindex .gt. 1 .and. l .lt. rcpara%brmax)
          if(in_gap) then
             l=l+1
             gcount(ipar)=l
             do while(mhilf(mindex,ipar) .lt. maxlen/2-max_miss+10 .and. mindex .gt. 1)
                mindex=mindex-1
             enddo
             lasts(l,ipar)=mindex
             !!          do while(mhilf(mindex,ipar) .lt. dfak*rcpara%snht_maxlen/3 .and. mhilf(mindex,ipar) .gt. rcpara%snht_maxlen/2-rcpara%max_miss+1 .and. mindex .gt. 1)
             do while(lasts(l,ipar)-mindex .lt. rcpara%mean_maxlen .and. mhilf(mindex,ipar) .gt. maxlen/2-max_miss+1 .and. mindex .gt. 1)
                mindex=mindex-1
             enddo
             if(mindex .lt. lasts(l,ipar)-rcpara%mean_maxlen) then 
!!$ call omp_set_lock(omp_lp)
                write(*,*) 'adjusting lasts',lasts(l,ipar),mindex
!!$ call omp_unset_lock(omp_lp)
                if(mindex+rcpara%mean_maxlen/2 .lt. rcpara%nmax-1) then 
                   lasts(l,ipar)=mindex+rcpara%mean_maxlen/2
                else
                   lasts(l,ipar)=rcpara%nmax-2
                endif

             endif
             in_gap=.not. in_gap
             midx(gcount(ipar),ipar)=mindex
             !!          write(*,*) gcount(ipar),mindex
             !!          if(lasts(l,ipar) .gt. 1) write(*,*) ipar,' before gap',l,' :',todate(midx(l,ipar),rcpara),todate(lasts(l,ipar),rcpara)
          else
             do while((mhilf(mindex,ipar) .ge. (maxlen/2-max_miss)) .and. mindex .gt. 1)
                mindex=mindex-1
             enddo
             in_gap=.not. in_gap
             !!          lasts(l,ipar)=mindex
             !!          write(*,*) gcount(ipar),mindex
          endif
          if(mindex .eq. 1 .and. l .eq. 1) midx(gcount(ipar),ipar)=rcpara%nmax
       enddo
       !!      gcount(ipar)=1 ! adjust only last part of series
       gcount(ipar)=gcount(ipar)-1
    enddo

  end subroutine detect_gaps

  subroutine select_composite(rcpara,statnr,cmax,wmonrs,wmolons,wmolats,wmostats,dists,index) !


    implicit none

    type(rasocor_namelist),intent(in) :: rcpara

    integer istat,i,j,statnr,err,iunit,cmax,wmostats
    integer index(rcpara%statmax),wmonrs(rcpara%statmax)
    integer,target:: lindex(rcpara%statmax)
    real :: xswap,dists(rcpara%statmax)
    real,target:: ldists(rcpara%statmax)
    real(kind=JPRM) :: wmolons(rcpara%statmax),wmolats(rcpara%statmax)
    character*3 fstring

    !! dists already calculated?
    if(index(1) .ne. 0) return

    istat=0
    do i=1,wmostats !rcpara%statmax
       if(statnr .eq. wmonrs(i)) istat=i
       !  index(i)=i
    enddo

    !call sphdist(wmolats(istat),wmolons(istat),wmolats,wmolons,dists,rcpara%statmax) !in file rfmod.f90 line 2257
    call sphdist(wmolats(istat),wmolons(istat),wmolats,wmolons,dists,wmostats) !in file rfmod.f90 line 2257


    ldists=dists

    call qsort(ldists(1:wmostats),lindex(1:wmostats))

    index=lindex

!!$  do j=1,cmax
!!$   do i=wmostats-1,1,-1
!!$    if(dists(index(i)) .gt. dists(index(i+1))) then
!!$      xswap=index(i) 
!!$      index(i)=index(i+1)
!!$      index(i+1)=xswap
!!$    endif
!!$   enddo
!!$  enddo

    if(wmonrs(index(1)) .ne. statnr) then
       i=1
       xswap=index(i) 
       index(i)=index(i+1)
       index(i+1)=xswap
       if(wmonrs(index(1)) .ne. statnr) then
          write(*,*) 'problem in select_composite',statnr
       else
          write(*,*) 'problem solved in select_composite',statnr
       endif
    endif

!!!$ call omp_set_lock(omp_lp)
    !!  write(fstring,'(I4)') cmax
    !!  write(*,'('//fstring//'I6)') wmonrs(index(1:cmax))
    !!  write(*,'('//fstring//'F6.2)') dists(index(1:cmax))
!!!$    call omp_unset_lock(omp_lp)

    return
  end subroutine select_composite


  subroutine adjust_bg(cstatnr,istat,rcpara,tm,bgcorrs,densities,adjust,eracorrs,wmolons,wmolats,wmostats,iunit) !

    implicit none

    type(rasocor_namelist),intent(in) :: rcpara

    integer protunit,iunit,i,ib,iname,l
    integer ipar,ip,wmostats,istat

    real(kind=JPRM) :: wmolons(wmostats),wmolats(wmostats),x
    real(kind=JPRM) :: bgcorrs(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM),intent(in) :: tm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM) :: eracorrs(rcpara%nmax,rcpara%pmax,rcpara%parmax)

    real(kind=JPRM) ::densities(rcpara%nmax,wmostats)
    real(kind=JPRM) :: adjust(rcpara%pmax,rcpara%parmax),cosfak

    real(kind=JPRM) :: tsa(rcpara%brmax,rcpara%pmax,rcpara%probmax),angle(rcpara%brmax),fak(rcpara%brmax)
    integer(kind=JPRM) :: era40meta(rcpara%nmax)

    real(kind=JPRM) :: amp(rcpara%brmax,rcpara%pmax),hilf1d(rcpara%nmax),hilf2d(rcpara%nmax,rcpara%pmax),null(rcpara%nmax)
    integer :: err,mindex,xlev,new_maxlen,new_max_miss,j,left_maxlen,right_maxlen,left,right,philf(rcpara%nmax),mhilf(rcpara%nmax),mgood(rcpara%pmax,rcpara%parmax,3)

    real(kind=JPRM) :: iadjust,eadjust,vprof(rcpara%pmax)

    integer :: cb(rcpara%brmax),ncritical,initial,final

    character*6  :: cstatnr
    character*2    ::ch
    character :: filename*80


    !!Werte von vprof, cosfak frs paper
    vprof=(/0.8,0.8,0.8,0.8,0.8,1.0,1.0,1.0,1.0,0.8,0.5,0.3,0.3,1.0,1.0,1.0/)
    cosfak=(2.+2*cos(wmolats(istat)*3.1415729/180.)**2)/3.

    vprof=(/1.0,1.0,1.0,1.0,1.0,1.2,0.9,0.7,0.7,0.6,1.0,1.2,1.2,1.0,1.0,1.0/)
    vprof=(/1.0,1.0,1.0,1.0,1.0,1.4,0.9,0.7,0.7,0.8,1.0,1.0,1.0,1.0,1.0,1.0/)
    !!Werte von vprof, cosfak fuer 004
    !!vprof=vprof/vprof
    !!cosfak=1.0
!!$ call omp_set_lock(omp_lp)
    ! write(*,*) 'COSFAK:',cstatnr,cosfak,vprof
!!$ call omp_unset_lock(omp_lp)


    eracorrs=0.
    do ipar=1,rcpara%parmax
       do ip=1,rcpara%pmax
          do ib=1,1
             if(rcpara%bg_correction_factor .gt. 0) then
                eracorrs(rcpara%bg_initial(ib):rcpara%bg_final(ib),ip,ipar)=bgcorrs(rcpara%bg_initial(ib):rcpara%bg_final(ib),ip,ipar)*densities(rcpara%bg_initial(ib):rcpara%bg_final(ib),istat)*rcpara%bg_correction_factor*vprof(ip)*cosfak !! avoid overcorrection
                if(rcpara%bg_initial(ib) .gt. 364 .and. rcpara%bg_final(ib) .lt. 18627) then
                   iadjust=sum(bgcorrs(rcpara%bg_initial(ib)-364:rcpara%bg_initial(ib),ip,ipar)*densities(rcpara%bg_initial(ib)-364:rcpara%bg_initial(ib),istat)*1.0)/365*vprof(ip)*cosfak
                   if(ib .eq. 1) then
                      eadjust=sum(bgcorrs(rcpara%bg_final(ib):rcpara%bg_final(ib)+364,ip,ipar)*densities(rcpara%bg_final(ib):rcpara%bg_final(ib)+364,istat)*1.0)/365*vprof(ip)*cosfak
                   else
                      eadjust=0.
                   endif
                   eracorrs(1:rcpara%bg_final(ib),ip,ipar)=eracorrs(1:rcpara%bg_final(ib),ip,ipar)-eadjust
                   eracorrs(1:rcpara%bg_initial(ib),ip,ipar)=eracorrs(1:rcpara%bg_initial(ib),ip,ipar)+iadjust
                endif
                !!      else 
                !!        eracorrs(initial(ib):final(ib),ip,ipar)=0.
             endif
          enddo
          if(rcpara%innov .ne. 'ME' ) eracorrs(1:16072,ip,ipar)=eracorrs(1:16072,ip,ipar)+adjust(ip,ipar)
       enddo
    enddo

    return

  end subroutine adjust_bg

  subroutine adjust_bg_ei(cstatnr,istat,rcpara,tm,bgcorrs,densities,adjust,eracorrs,wmolons,wmolats,wmostats,iunit) !

    implicit none

    type(rasocor_namelist),intent(in) :: rcpara

    integer protunit,iunit,i,ib,iname,l,im
    integer ipar,ip,wmostats,istat

    real(kind=JPRM) :: wmolons(wmostats),wmolats(wmostats),x
    real(kind=JPRM) :: bgcorrs(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM),intent(in) :: tm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM) :: eracorrs(rcpara%nmax,rcpara%pmax,rcpara%parmax)

    real(kind=JPRM) ::densities(rcpara%nmax,wmostats)
    real(kind=JPRM) :: adjust(rcpara%pmax,rcpara%parmax,12),cosfak

    real(kind=JPRM) :: tsa(rcpara%brmax,rcpara%pmax,rcpara%probmax),angle(rcpara%brmax),fak(rcpara%brmax)
    integer(kind=JPRM) :: era40meta(rcpara%nmax)

    real(kind=JPRM) :: amp(rcpara%brmax,rcpara%pmax),hilf1d(rcpara%nmax),hilf2d(rcpara%nmax,rcpara%pmax),null(rcpara%nmax)
    integer :: err,mindex,xlev,new_maxlen,new_max_miss,j,left_maxlen,right_maxlen,left,right,philf(rcpara%nmax),mhilf(rcpara%nmax),mgood(rcpara%pmax,rcpara%parmax,3)

    real(kind=JPRM) :: iadjust,eadjust,vprof(rcpara%pmax),aprof(rcpara%pmax),adj

    integer :: cb(rcpara%brmax),ncritical,initial,final,avgint

    character*6  :: cstatnr
    character*2    ::ch
    character :: filename*80


    !!Werte von vprof, cosfak frs paper
    vprof=(/0.8,0.8,0.8,0.8,0.8,1.0,1.0,1.0,1.0,0.8,0.5,0.3,0.3,1.0,1.0,1.0/)
    cosfak=1! (2.+2*cos(wmolats(istat)*3.1415729/180.)**2)/3.

    vprof=(/1.0,1.0,1.0,1.0,1.0,1.2,0.9,0.7,0.7,0.6,1.0,1.2,1.2,1.0,1.0,1.0/)
    vprof=(/1.0,1.0,1.0,1.0,1.0,1.4,0.9,0.7,0.7,0.8,1.0,1.0,1.0,1.0,1.0,1.0/)
    vprof=(/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/)
    !aprof=(/1.0,1.0,1.0,0.8,0.5,0.5,0.8,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/)
    !aprof=(/1.0,1.0,1.0,1.0,0.9,0.5,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/)
    !aprof=(/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.5,0.0,0.0,0.0/)
    aprof=(/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/)
    aprof=(/1.0,0.8,1.0,1.0,0.7,0.4,0.6,0.9,1.0,1.0,1.0,1.0,1.0,0.5,1.0,1.0/)
    !!Werte von vprof, cosfak fuer 004
    !!vprof=vprof/vprof
    !!cosfak=1.0
!!$ call omp_set_lock(omp_lp)
    ! write(*,*) 'COSFAK:',cstatnr,cosfak,vprof
!!$ call omp_unset_lock(omp_lp)


    avgint=1440
    eracorrs=0.
    do ipar=1,rcpara%parmax
       do ip=1,rcpara%pmax
          do ib=1,1
             if(rcpara%bg_correction_factor .gt. 0) then
                eracorrs(rcpara%bg_initial(ib):rcpara%bg_final(ib),ip,ipar)=bgcorrs(rcpara%bg_initial(ib):rcpara%bg_final(ib),ip,ipar)*densities(rcpara%bg_initial(ib):rcpara%bg_final(ib),istat)*rcpara%bg_correction_factor*vprof(ip)*cosfak !! avoid overcorrection
                if(rcpara%bg_initial(ib) .gt. avgint .and. rcpara%bg_final(ib) .lt. 18627) then
                   iadjust=sum(bgcorrs(rcpara%bg_initial(ib)-avgint:rcpara%bg_initial(ib),ip,ipar)*densities(rcpara%bg_initial(ib)-avgint:rcpara%bg_initial(ib),istat)*1.0)/(avgint+1)*vprof(ip)*cosfak
                   if(ib .eq. 1) then
                      eadjust=sum(bgcorrs(rcpara%bg_final(ib):rcpara%bg_final(ib)+avgint,ip,ipar)*densities(rcpara%bg_final(ib):rcpara%bg_final(ib)+avgint,istat)*1.0)/(avgint+1)*vprof(ip)*cosfak
                   else
                      eadjust=0.
                   endif
                   eracorrs(1:rcpara%bg_final(ib),ip,ipar)=eracorrs(1:rcpara%bg_final(ib),ip,ipar)-eadjust
                   eracorrs(1:rcpara%bg_initial(ib),ip,ipar)=eracorrs(1:rcpara%bg_initial(ib),ip,ipar)+iadjust
                endif
                !!      else 
                !!        eracorrs(initial(ib):final(ib),ip,ipar)=0.
             endif
          enddo

       enddo
    enddo

    ib=toindex(rcpara%switchdate,rcpara)-1
    do ipar=1,rcpara%parmax
       do ip=1,rcpara%pmax
          adj=0
          l=0
          do im=1,12
             if(adjust(ip,ipar,im) /= rcpara%miss_val) then
                adj=adj+adjust(ip,ipar,im)
                l=l+1
             endif
          enddo
          if(l >0) adj=adj/l
          do i=1,ib
             eracorrs(i,ip,ipar)=eracorrs(i,ip,ipar)+aprof(ip)*adj
             !        eracorrs(i,ip,ipar)=aprof(ip)*adj
          enddo
       enddo
    enddo

    return


  end subroutine adjust_bg_ei


  subroutine make_composite(rcpara,statnr,cmax,wmonrs,wmolons,wmolats,wmostats,dists,index,ominuse40,adjust,stm,stfgm,stnum,&
       tfgmcr,tmcr,icache,meta_s,lasts,gcount,needs_composite,ini_correct,composite_exists,bad_intervals,tccr) !

    implicit none

    type(rasocor_namelist),intent(in) :: rcpara
    type(cacherecord) :: tfgmcr(rcpara%cachemax),tmcr(rcpara%cachemax)
    type(cacherecord),optional ::tccr(rcpara%cachemax,3)

    type(metadata) :: meta_s

    integer istat,i,j,statnr,err,iunit,cmax,imin,ib,minst,found,l,istart,istop,mstart,ipresatend
    integer,intent(in) :: wmostats,wmonrs(rcpara%statmax)
    integer index(rcpara%statmax),switchindex
    real :: dists(rcpara%statmax)
    real(kind=JPRM),intent(in) :: wmolons(rcpara%statmax),wmolats(rcpara%statmax)

    real(kind=JPRM) ,intent(in):: ominuse40(rcpara%ni,rcpara%nj,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM),intent(in) :: adjust(rcpara%pmax,rcpara%parmax,12)
    real(kind=JPRM) :: adjustlocal(rcpara%pmax,rcpara%parmax,12)

    real(kind=JPRM) :: mrasocorrs(rcpara%mmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM) :: tfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tanm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tbcm(rcpara%nmax,rcpara%pmax,rcpara%parmax),itfg12m(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    !    real(kind=JPRM),intent(in) :: ttestm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM) :: stfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax),stm(rcpara%nmax,rcpara%pmax,rcpara%parmax),stweight(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM) :: eracorrs(rcpara%nmax,rcpara%pmax,rcpara%parmax),rasocorrs(rcpara%nmax,rcpara%pmax,rcpara%parmax),weight,wsum,wcrit
    integer :: stnum(rcpara%nmax,rcpara%pmax,rcpara%parmax)

    integer mr,cachehilf(rcpara%nmax)

    integer :: bindex(rcpara%mmax),bi,ipmax,iparmax,ip,ipar,g_days,offset
    integer :: rbindex(rcpara%mmax),rbi
    real(kind=JPRM) :: hilfcorr(rcpara%mmax,rcpara%pmax,rcpara%parmax),stnumcrit(rcpara%brmax,rcpara%parmax)
    logical ex3,ex4,ldebug,lexit
    integer :: local_cardsmeta(rcpara%nmax,1,rcpara%nmeta)
    integer err2,err3,ios,ic(rcpara%parmax),omp_get_thread_num, needs_composite(:,:,:)
    character cstatnr*6,filename*120,routine_name*200,cf*2
    logical composite_exists(rcpara%parmax),logcache
    integer,intent(in)         :: lasts(rcpara%brmax,rcpara%parmax),gcount(rcpara%parmax)
    logical,intent(in) :: ini_correct(rcpara%parmax,rcpara%statmax)
    integer         :: llasts(rcpara%brmax,rcpara%parmax),lgcount(rcpara%parmax),lmidx(rcpara%brmax,rcpara%parmax),imax,icache(rcpara%statmax+1),icistat
    !    real(kind=JPRM) :: bgcorr(rcpara%nmax,rcpara%pmax,rcpara%parmax),densities(rcpara%nmax,wmostats)
    real(kind=JPRM) :: pi,stfgmsum
    integer:: bad_intervals(:,:),stfgsum,ncb(rcpara%parmax)

    TYPE(complete_station)::net_complete_stat

    TYPE(global_settings_raobcore) :: g_s,g_s_in !type used to store al

    ldebug=.false.
    call select_composite(rcpara,statnr,cmax,wmonrs,wmolons,wmolats,wmostats,dists,index) !in file rfcor.f90 line 206

    err3=0

    stnum=0
    stfgm=0.
    stm=0.
    stweight=0.
    wsum=0.

    if (wmolats(index(1)) .lt. -55. ) then
       minst=10
    else
       minst=20
    endif
    !! create composite from nearest cmax stations
    !! index(1) is station itself, therefore omitted
    istat=1
    found=0
    pi=acos(-1.0)
    do while(istat .lt. cmax) ! .and. abs(cos(wmolats(index(1))*pi/180.)-cos(wmolats(index(istat))*pi/180.)) .lt. 0.5) 

       istat=istat+1


       write(cstatnr,'(I6.6)') wmonrs(index(istat))

       filename=trim(rcpara%prefix)//cstatnr//'/feedbackmerged'//cstatnr//'.nc'
       inquire(file=filename,exist=ex3)

       filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbincorrsave'//cstatnr//'.nc'
       inquire(file=filename,exist=ex4)


       if(( wmonrs(index(istat)) .gt. 42000 .and. wmonrs(index(istat)) .lt. 44000 .or. wmonrs(index(istat)) .gt. 17000 .and. wmonrs(index(istat)) .lt. 17500 .or. wmonrs(index(istat)) .gt. 15000 .and. wmonrs(index(istat)) .lt. 16000)) then

          if(ldebug) write(*,*) cstatnr, 'Arabian sonde - not used'
          ex3=.false.
       endif
       call check_ini(meta_s,index(istat),wmonrs,rcpara,ini_correct)
       !write(*,*) wmonrs(index(1)),wmonrs(index(istat)),ini_correct(:,index(istat))

       !Wenn die betrachtete Station südlich ist und die andere auch, oder nicht südlich (514)
       if( ((wmolats(index(1)) .lt. -55. .and.  wmolats(index(istat)) .lt. -50.) .or. wmolats(index(1)) .ge. -55.) .and. any(.not. ini_correct(:,index(istat))) .and. ex3 .and. ex4) then 
          !  if( ((wmolats(index(1)) .lt. -55. .and.  wmolats(index(istat)) .lt. -50.) .or. wmolats(index(1)) .ge. -55.)  .and. ex3 .and. ex4) then 

          if(ldebug)    then
             write(*,*) 'stnumcrit:',statnr,':',wmonrs(index(istat)),ini_correct(:,index(istat)),istat,found,stnumcrit(1:3,:)
          endif




          !$ call omp_set_lock(omp_lp(rcpara%statmax+2))
          logcache=icache(index(istat)) .eq. 0
          !$ call omp_unset_lock(omp_lp(rcpara%statmax+2))

          if(logcache) then !524
             !$ call omp_set_lock(omp_lp(rcpara%statmax+2))
             icache(rcpara%statmax+1)=icache(rcpara%statmax+1)+1
             icache(index(istat))=icache(rcpara%statmax+1)
             !$ call omp_unset_lock(omp_lp(rcpara%statmax+2))

             iunit=500000
             !$    iunit=iunit+OMP_GET_THREAD_NUM()

             filename=trim(rcpara%prefix)//cstatnr//'/feedbackmerged'//cstatnr//'.nc'
             CALL read_odb_nc(filename,rcpara,index(istat),err3,tm,tfgm,bad_intervals=bad_intervals) !subroutine in read_txt_write_nc.f90
             do ipar=1,rcpara%parmax
                do ip=1,rcpara%pmax
                   do i=1,rcpara%nmax
                      if (isnan(tm(i,ip,ipar)) .or. isnan(tfgm(i,ip,ipar))  .or. abs(tfgm(i,ip,ipar))>20.) then
                         tm(i,ip,ipar)=rcpara%miss_val
                         tfgm(i,ip,ipar)=rcpara%miss_val
                      else
                         !                         tfgm(i,ip,ipar)=-tfgm(i,ip,ipar)
                      endif
                   enddo
                enddo
             enddo

             call eiminuse40(rcpara,wmolats(index(istat)),wmolons(index(istat)),ominuse40,adjustlocal)
             switchindex=toindex(rcpara%switchdate,rcpara)
             !     do ipar=1,rcpara%parmax
             !       do ip=1,rcpara%pmax
             !         do i=1,switchindex-1
             !           if(tfgm(i,ip,ipar) .ne. rcpara%miss_val) tfgm(i,ip,ipar)=tfgm(i,ip,ipar)-adjustlocal(ip,ipar,rcpara%month(i))
             !         enddo
             !       enddo
             !     enddo   

             !! only use most recent parts of time series, but no data before a gap 
             !call detect_gaps(rcpara,wmonrs(index(istat)),tfgm,lmidx,llasts,lgcount) !in file rfcor.f90 line 113


             filename=trim(trim(rcpara%prefix)//cstatnr//'/feedbackglobbincorrsave'//cstatnr//'.nc')

             CALL read_sonde_corr_daily_IFS_nc(filename, rcpara, index(istat),ios, hilfcorr, bindex,bi)

             if(ios .eq. 0) then  !if 567

                rasocorrs=0.
                do ipar=1,rcpara%parmax
                   do ip=1,rcpara%pmax
                      do i=1,bi-1
                         rasocorrs(bindex(i):bindex(i+1),ip,ipar)=-hilfcorr(i,ip,ipar)
                      enddo
                   enddo
                enddo

                icistat=icache(index(istat))
                !$ call omp_set_lock(omp_lp(icistat))
                l=0
                do i=1,rcpara%nmax
                   if((tfgm(i,12,1) .ne. rcpara%miss_val .or. tfgm(i,12,2) .ne. rcpara%miss_val) .and.  (tm(i,12,1) .ne. rcpara%miss_val .or.  tm(i,12,2) .ne. rcpara%miss_val)) then
                      l=l+1            
                      cachehilf(l)=i
                   endif
                enddo
                tfgmcr(icistat)%vals=l
                allocate(tfgmcr(icistat)%index(l))
                tfgmcr(icistat)%index=cachehilf(1:l)
                allocate(tfgmcr(icistat)%feld(l,rcpara%pmax,rcpara%parmax))
                allocate(tmcr(icistat)%feld(l,rcpara%pmax,rcpara%parmax))

                do ipar=1,rcpara%parmax
                   do ip=1,rcpara%pmax
                      do l=1,tfgmcr(icistat)%vals
                         i=tfgmcr(icistat)%index(l)

                         if(tfgm(i,ip,ipar) .ne. rcpara%miss_val .and. tm(i,ip,ipar) .ne. rcpara%miss_val.and. rasocorrs(i,ip,ipar) .ne. rcpara%miss_val ) then

                            tfgmcr(icistat)%feld(l,ip,ipar)=tfgm(i,ip,ipar)+rasocorrs(i,ip,ipar) !+eracorrs(i,ip,ipar)
                            tmcr(icistat)%feld(l,ip,ipar)=tm(i,ip,ipar)-rasocorrs(i,ip,ipar)
                            if(abs(tfgmcr(icistat)%feld(l,ip,ipar)).gt. 30.) then
                               write(*,*) cstatnr,l,i,ip,ipar, tfgmcr(icistat)%feld(l,ip,ipar),' make_composite: wrong temperature'
                               stop
                            endif
                         else
                            tfgmcr(icistat)%feld(l,ip,ipar)=rcpara%miss_val
                            tmcr(icistat)%feld(l,ip,ipar)=rcpara%miss_val
                         endif
                      enddo
                   enddo
                enddo

                write(*,*) wmonrs(index(istat)),' make written to cache with index',icache(index(istat))
                !$ call omp_unset_lock(omp_lp(icistat))

             else !576
                write(*,*) cstatnr,'corrections could not be read - not used'
             endif !! ios .gt. 0

          else
             icistat=icache(index(istat))


          endif !! icache 524

          ncb=1
          do ipar=1,rcpara%parmax
             l=1
             if (needs_composite(l,ipar,index(istat))>1 .and. needs_composite(l,ipar,index(istat))<rcpara%old) then
                ncb(ipar)=needs_composite(l,ipar,index(istat))
             endif
          enddo
          ic=0
          do ipar=1,rcpara%parmax
             if(.not. ini_correct(ipar,index(istat))) then
                do ib=1,gcount(ipar)
                   i=1
                   if (tfgmcr(icistat)%vals>0) then
                      do while(i.lt. tfgmcr(icistat)%vals .and.tfgmcr(icistat)%index(i) .lt. lasts(ib,ipar)-rcpara%mean_maxlen)
                         i=i+1
                      enddo
                      do while(i.lt. tfgmcr(icistat)%vals .and. tfgmcr(icistat)%index(i) .le. lasts(ib,ipar))
                         if(tfgmcr(icistat)%feld(i,8,ipar) .ne. rcpara%miss_val .and. tfgmcr(icistat)%index(i)>ncb(ipar) .and. .not. ini_correct(ipar,index(istat))) then
                           ic(ipar)=ic(ipar)+1
                         endif
                         i=i+1
                      enddo
                   endif
                enddo
             endif
          enddo

          if(ldebug)      write(*,'(2I6,2L2,A5,2I6,A5,2I6,A5,2I6)') wmonrs(index(1)),wmonrs(index(istat)),ini_correct(:,index(istat)),'last: ',lasts(1,:),' ic: ',ic,' ncb: ',ncb

          if(present(tccr)) then
             call load_richcorr(tccr,icistat,0,wmonrs,index(istat),rcpara, mrasocorrs,.false.,rbindex,rbi)
             do ipar=1,rcpara%parmax
                do ip=1,rcpara%pmax
                   do i=1,rbi-1
                      rasocorrs(rbindex(i):rbindex(i+1),ip,ipar)=-mrasocorrs(i,ip,ipar)
                   enddo
                enddo
             enddo
          else
             rasocorrs=0.
          endif

          weight=exp(-0.01-dists(index(istat))*3.1415729/180.*6370./rcpara%weight_distance)

          if(any(ic .gt. rcpara%mean_maxlen/2)) then !651
             found=found+1
             if(icistat .eq. 0) then
                write(*,'(2I6,A1)') wmonrs(index(1)),wmonrs(index(istat)),': ',' cache error'
             else
                if(ldebug)      write(*,'(2I6,A2,A7,I6)') wmonrs(index(1)),wmonrs(index(istat)),': ',' cache ',icistat
             endif

             mstart=1
             if(wmonrs(index(istat)) .gt. 7000 .and. wmonrs(index(istat)) .lt. 8000 .or. &
                  wmonrs(index(istat)) .gt. 91900 .and. wmonrs(index(istat)) .lt. 92000 .or. &
                  wmonrs(index(istat)) .gt. 91500 .and. wmonrs(index(istat)) .lt. 91600 .or. &
                  wmonrs(index(istat)) .gt. 60000 .and. wmonrs(index(istat)) .lt. 68000 .or. & 
                  wmonrs(index(istat)) .gt. 80000 .and. wmonrs(index(istat)) .lt. 89000 .or. & 
                  wmonrs(index(istat)) .gt. 50000 .and. wmonrs(index(istat)) .lt. 60000 ) then
                mstart=toindex(19870101,rcpara)
             endif
             if(wmonrs(index(istat)) .gt. 89000 .and. wmonrs(index(istat)) .lt. 90000 .or. & 
                wmonrs(index(istat)) .gt. 93000 .and. wmonrs(index(istat)) .lt. 95000 ) then
                mstart=toindex(19540101,rcpara)
             endif

             do ipar=1,rcpara%parmax
                do ip=1,rcpara%pmax-2
                   do l=1,tfgmcr(icistat)%vals
                      if(tfgmcr(icistat)%feld(l,ip,ipar) .ne. rcpara%miss_val .and. tmcr(icistat)%feld(l,ip,ipar) .ne. rcpara%miss_val) then
                         i=tfgmcr(icistat)%index(l)
                         if(i .ge. mstart) then
                            stfgm(i,ip,ipar)=stfgm(i,ip,ipar)+weight*(tfgmcr(icistat)%feld(l,ip,ipar)+rasocorrs(i,ip,ipar))
!!$                            if(tmcr(icistat)%feld(l,ip,ipar) .lt. 150. .or. tmcr(icistat)%feld(l,ip,ipar) .gt. 350.) then
!!$                               write(*,*) 'STOP: spurious tmcr' , l,ip,ipar,tmcr(icistat)%feld(l,ip,ipar)
!!$                               stop
!!$                            else
                               stm(i,ip,ipar)=stm(i,ip,ipar)+weight*(tmcr(icistat)%feld(l,ip,ipar)+rasocorrs(i,ip,ipar))
!!$                            endif
                            stnum(i,ip,ipar)=stnum(i,ip,ipar)+1
                            stweight(i,ip,ipar)=stweight(i,ip,ipar)+weight
                         endif
                      endif
                   enddo
                enddo
             enddo
             wsum=wsum+weight
             if(ldebug .and. (statnr .eq. 83229 .or. statnr .eq. 1152 ))      then 
                write(*,'(2I6,A2,2I6,F8.1,4F8.4,I7,F9.4,2I7)') wmonrs(index(1)),wmonrs(index(istat)),':: ',ic,dists(index(istat))*6370.*3.1415729/180.,weight,maxval(stfgm),sum(stfgm)/rcpara%pmax/rcpara%parmax/rcpara%nmax,maxval(tfgmcr(icistat)%feld),count(tfgmcr(icistat)%feld .ne. rcpara%miss_val),maxval(tmcr(icistat)%feld),count(tmcr(icistat)%feld .ne. rcpara%miss_val),count(stfgm .ne. 0)
                continue

             endif
          else
             if(ldebug)      write(*,*) cstatnr,' not used in composite'
          endif !! ic .gt. mean_maxlen/2 (651)

       else
       endif  !! ex3 (514)
       err=err3 !!*err1*err2

       lexit=.true.
       do ipar=1,rcpara%parmax
          do ib=1,gcount(ipar)
             imin=maxval((/1,lasts(ib,ipar)-rcpara%mean_maxlen/))
             if(lasts(ib,ipar) .lt. 2) then         
                stnumcrit(ib,ipar)=0.
             else
                stnumcrit(ib,ipar)=sum(stnum(imin:lasts(ib,ipar),6,ipar))*1.0/(lasts(ib,ipar)-imin+1)
                if  (stnumcrit(ib,ipar)<0.9*minst) then
                   lexit=.false.
                endif
             endif

             if(ldebug) then
                stfgmsum=0.
                stfgsum=0
                do i=imin,lasts(ib,ipar)
                   if (stweight(i,6,ipar).ne. 0.) then
                      stfgmsum=stfgmsum+stfgm(i,6,ipar)/stweight(i,6,ipar)
                      stfgsum=stfgsum+1
                   endif
                enddo
                if (stfgsum>0) then
                   write(*,'(a,1x,I6,I3,I3,I6,I6,F6.2,L2,F6.2)') 'stnumcrit:',wmonrs(index(1)),ib,ipar,imin,lasts(ib,ipar),stnumcrit(ib,ipar),lexit,stfgmsum/stfgsum   
                endif
             endif
          enddo
       enddo
       if(lexit) then
         exit
       endif

    enddo  !485

    write(*,'(I6,A22,I4,A6,2I9,F8.2)') wmonrs(index(1)),': stations searched: ',istat,' good:',count(stweight(:,4,1) .gt. wsum/2.),count(stweight(:,4,1) .lt. wsum/2.),wsum/2.
    if (istat .eq. wmostats) then 
       write(*,*) 'stations searched',istat .lt. cmax, found .lt. minst, any(stnumcrit .lt. 0.7*minst), wmolats(index(1)) .gt. -50., dists(index(istat))*3.1415729/180.*6370. .lt. 5000.
    endif



    do ipar=1,rcpara%parmax !713
       do ib=1,gcount(ipar) !714
          if(lasts(ib,ipar) .gt. 1) then !715
             imin=maxval((/1,lasts(ib,ipar)-rcpara%mean_maxlen/))
             do ip=1,rcpara%pmax !717
                wcrit=sum(stweight(imin:lasts(ib,ipar),ip,ipar))/(lasts(ib,ipar)-imin+1)
                istart=lasts(ib+1,ipar)+rcpara%snht_maxlen/2+1
                if(lasts(ib+1,ipar) .eq. 1) istart=1
                istop=lasts(ib,ipar)+rcpara%snht_maxlen/2
                if(ib .eq. 1) istop=rcpara%nmax
                do i=istart,istop !723
                   if(stfgm(i,ip,ipar) .ne. rcpara%miss_val .and. stweight(i,ip,ipar) .gt. 1.) then !*(i+1000.)/(lasts(ib,ipar)+rcpara%snht_maxlen/2+1000.)/2.) then  ! .and. stweight(i,ip,ipar) .gt. wsum/5.
                      stfgm(i,ip,ipar)=stfgm(i,ip,ipar)/stweight(i,ip,ipar)
                      stm(i,ip,ipar)=stm(i,ip,ipar)/stweight(i,ip,ipar)
                   else
                      stfgm(i,ip,ipar)=rcpara%miss_val
                      stm(i,ip,ipar)=rcpara%miss_val
                   endif
                enddo !723

                stnumcrit(1,1)=sum(stnum(imin:lasts(ib,ipar),ip,ipar))*1.0/(lasts(ib,ipar)-imin+1)
                if(ldebug)       write(*,'(A7,I6,I3,I9,I3,2I4,2F8.2)') 'rcomp: ',statnr,ipar,todate(lasts(ib,ipar),rcpara),ip,istat-1,found,stnumcrit(1,1),wcrit

             enddo !717
          endif !715
       enddo !714
       if(gcount(ipar) .eq. 0) then
          stfgm(:,:,ipar)=rcpara%miss_val
          stm(:,:,ipar)=rcpara%miss_val
       endif
    enddo !713
    if(ldebug)      write(*,'(I6,A6,2I9)') wmonrs(index(1)),': good',count(stfgm(:,4,1) .gt. rcpara%miss_val),count(stfgm(:,4,1) .le. rcpara%miss_val)
!!$    do ipar=1,rcpara%parmax
!!$       do ip=1,rcpara%pmax
!!$          do i=1,rcpara%nmax
!!$             if(stm(i,ip,ipar) .ne. rcpara%miss_val) then
!!$                if(stm(i,ip,ipar) .lt. 150. .or. stm(i,ip,ipar) .gt. 350) then
!!$                   write(*,*) wmonrs(index(1)),i,ip,ipar,stm(i,ip,ipar),'STOP: stm invalid '
!!$                   call exit(1)
!!$                endif
!!$             endif
!!$          enddo
!!$       enddo
!!$    enddo


    composite_exists=.true.

    return

  end subroutine make_composite

  subroutine save_state(iunit,filename,rcpara,tm,mean_of_tsa,tsa_of_mean,apriori_probs,meanbreak_probs,breakmean_probs,chosenbreaks) !

    implicit none

    type(rasocor_namelist),intent(in) :: rcpara

    character*80 filename
    integer iunit,err,chosenbreaks(rcpara%nmax),il,i
    integer goodindex(rcpara%nmax)
    !!real(kind=JPRM) :: tm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    !!real(kind=JPRM) :: tsa(rcpara%nmax,rcpara%pmax,rcpara%probmax)
    real(kind=JPRM) :: tsa_of_mean(rcpara%nmax,rcpara%probmax),mean_of_tsa(rcpara%nmax,rcpara%probmax),meanbreak_probs(rcpara%nmax,rcpara%probmax),breakmean_probs(rcpara%nmax,rcpara%probmax),apriori_probs(rcpara%nmax)
    real(kind=JPRM) :: tm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM),allocatable :: hilf(:,:)


    open(iunit,file=filename,form='unformatted',err=120)
    err=0

    il=1
    do i=1,rcpara%nmax
       if(any(tm(i,:,:) .gt. -999.)) then
          goodindex(il)=i
          il=il+1 
       endif
    enddo
    il=il-1
    allocate(hilf(il,rcpara%probmax))

    write(iunit) rcpara%nmax,il,rcpara%pmax,rcpara%parmax,rcpara%probmax

    write(iunit) goodindex(1:il)
    !!  hilf=mean_of_tsa(goodindex(1:il),:)
    !!  write(iunit) hilf
    hilf=tsa_of_mean(goodindex(1:il),:)
    write(iunit) hilf
    write(iunit) apriori_probs(goodindex(1:il))
    !!  hilf=meanbreak_probs(goodindex(1:il),:)
    !!  write(iunit) hilf
    hilf=breakmean_probs(goodindex(1:il),:)
    write(iunit) hilf
    write(iunit) chosenbreaks(1:30)

    deallocate(hilf)

    close(iunit)

    return
120 continue
!!$ call omp_set_lock(omp_lp)
    print*,'cannot write ',filename
!!$ call omp_unset_lock(omp_lp)

    return
  end subroutine save_state

  subroutine corr_stats(rcpara,tm,tfgm,stm,stfgm,protunit) !

    implicit none

    type(rasocor_namelist),intent(in) :: rcpara

    integer protunit,ib,bcount,ipar,ip,ic
    real(kind=JPRM) :: tm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax),stm(rcpara%nmax,rcpara%pmax,rcpara%parmax),stfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM) :: mean,meansq,sq(rcpara%nmax)
    logical :: mask(rcpara%nmax)
    do ipar=1,rcpara%parmax
       do ip=1,rcpara%pmax
          ic=count(tm(:,ip,ipar) .ne. rcpara%miss_val)
          if(ic .gt. 1) then 
             mask=tm(:,ip,ipar) .ne. rcpara%miss_val
             where(mask)
                sq=tm(:,ip,ipar)*tm(:,ip,ipar)
             endwhere
             meansq=sqrt(sum(sq,mask)/ic)
             mean=sum(tm(:,ip,ipar),mask)/ic

!!$ call omp_set_lock(omp_lp)
             if(meansq*meansq-mean*mean .gt. 0) then
                write(protunit,'(2I3,I6,3F8.3)') ipar,ip,ic,mean,meansq,sqrt(meansq*meansq-mean*mean)
             else
!!$ call omp_set_lock(omp_lp)
                write(protunit,'(2I3,I6,3F8.3)') ipar,ip,ic,0.,0.,0. !
             endif
!!$ call omp_unset_lock(omp_lp)
          else
!!$ call omp_set_lock(omp_lp)
             write(protunit,'(2I3,I6,3F8.3)') ipar,ip,ic,0.,0.,0. !
!!$ call omp_unset_lock(omp_lp)
          endif
       enddo
    enddo

  end subroutine corr_stats


  !! attributionfactors=1 if radiosonde
  !! attributionfactors=0 if ERA-40
  subroutine attribute_breaks(rcpara,cardsmeta,era40meta,probmeanbreaks,meanprobbreaks,chosenbreaks,attributionfactors,protunit)

    implicit none

    type(rasocor_namelist),intent(in) :: rcpara

    integer protunit,ib,bcount
    integer(kind=JPRM) :: meanprobbreaks(rcpara%nmax,rcpara%probmax),probmeanbreaks(rcpara%nmax,rcpara%probmax),chosenbreaks(rcpara%nmax)
    real(kind=JPRM) :: break_thresh(rcpara%probmax),attributionfactors(rcpara%nmax)
    integer cprobmean(rcpara%probmax),cmeanprob(rcpara%probmax)
    integer(kind=JPRM) :: cardsmeta(rcpara%nmax,1,rcpara%nmeta),era40meta(rcpara%nmax)

    bcount=count(chosenbreaks .gt. 0)
    attributionfactors=-999.

    do ib=1,bcount
       attributionfactors(ib)=1.
       if(era40meta(chosenbreaks(ib)) .gt. 0) then 
          attributionfactors(ib)=0.
       endif
    enddo

    return
  end subroutine attribute_breaks

  !! Bayes_prob combines a priori information from metadata with
  !! statistical information.
  !! The probability for a break in a particular interval is in most cases
  !! close to 1. What remains is to find the break location within that interval.
  !! For this purpose knowledge  from simulation experiments is used.
  subroutine bayes_prob(rcpara,apriori_probs,tsa,tsalocs,break_probs,protunit)

    implicit none

    type(rasocor_namelist),intent(in) :: rcpara

    integer IFAIL,ib,ip,bc,istart,istop,locint,i,protunit
    real(kind=JPRM) :: apriori_probs(rcpara%nmax),tsa(rcpara%nmax,rcpara%probmax)
    integer(kind=JPRM) :: tsalocs(rcpara%nmax,rcpara%probmax)
    real(kind=JPRM) :: break_probs(rcpara%nmax,rcpara%probmax),sig
    real xp,xm,G01EAF,prob

    break_probs=0.
    do ip=1,rcpara%probmax
       bc=count(tsalocs(:,ip) .ne. 0)
       do ib=1,bc
          !!     sig=-40.+rcpara%locsig/sqrt(tsa(tsalocs(ib,ip),ip))
          sig=rcpara%locsig/sqrt(tsa(tsalocs(ib,ip),ip))
          if(sig .lt. 2.) sig=2. !! avoid negative sig for extreme values of tsa
          locint=3*sig
          istart=max(1,tsalocs(ib,ip)-locint)
          istop=min(rcpara%nmax,tsalocs(ib,ip)+locint)
          do i=istart,istop
             xm=(i-tsalocs(ib,ip))/sig
             xp=(i-tsalocs(ib,ip)+1)/sig
             prob=G01EAF('L', xp, IFAIL)-G01EAF('L', xm, IFAIL)
             break_probs(i,ip)=prob
          enddo
       enddo
    enddo

    !! P(A)=aprioriprob
    !! P(B/A)=break_probs
    !! P(\A)=1-aprioriprob
    !! P(B/\A)=1-break_probs
    !! P(A/B)=P(A)P(B/A)/(P(A)P(B/A)+P(\A)P(B/\A))
    !!
    !! P(B/\A)=1/maxlen (uniform distribution of break location if there is no break)

    do ip=1,rcpara%probmax
       do i=1,rcpara%nmax
          break_probs(i,ip)=apriori_probs(i)*break_probs(i,ip)/(apriori_probs(i)*break_probs(i,ip)+(1.-apriori_probs(i))*1./rcpara%snht_maxlen)
       enddo
    enddo

    return
  end subroutine bayes_prob

  !! Bayes_prob combines a priori information from metadata with
  !! statistical information.
  !! The probability for a break in a particular interval is in most cases
  !! close to 1. What remains is to find the break location within that interval.
  !! For this purpose knowledge from simulation experiments is used.
  subroutine bayes_break_simple(rcpara,apriori_probs,tsa,break_probs,break_thresh,ipstart,ipend,protunit) !

    implicit none

    type(rasocor_namelist),intent(in) :: rcpara

    integer IFAIL,ib,ip,bc,istart,istop,locint,i,protunit,ipstart,ipend,ianf,iend,tsamaxloc,dim2
    real(kind=JPRM),intent(in) :: apriori_probs(:),tsa(rcpara%nmax,rcpara%pmax)
    real(kind=JPRM) :: tsahilf(rcpara%nmax),break_probs(:),sig,break_prob,redfak,sigfak,tsamax,break_thresh
    real xp,xm,G01EAF,prob,nprob
    logical :: reduced(rcpara%nmax)

    redfak=0.95
    sigfak=1.6
    break_probs=apriori_probs
    reduced=.false.
    do ip=ipend,ipstart,-1
       where(tsa(:,ip) .ne. rcpara%miss_val)
          tsahilf=tsa(:,ip)
       elsewhere
          tsahilf=0.
       endwhere
       if(ip .ne. 9) then 
          tsamax=maxval(tsahilf)
          do while(tsamax .gt. break_thresh)
             tsamaxloc=maxloc(tsahilf,1)
             ianf=max(1,tsamaxloc-rcpara%snht_maxlen/2)
             iend=min(rcpara%nmax,tsamaxloc+rcpara%snht_maxlen/2)
             where(.not. reduced(ianf:iend))
                tsahilf(ianf:iend)=tsahilf(ianf:iend)/(tsamax+1.)*break_thresh
                reduced(ianf:iend)=.true.
             endwhere
             tsamax=maxval(tsahilf)    
          enddo


          do i=1,rcpara%nmax
             if(tsahilf(i) .ne. 0.) then
                if(ipend .eq. ipstart) then
                   xp=-3.+5.5*tsahilf(i)/break_thresh
                   !!         xp=tsahilf(i)!/break_thresh
                else
                   xp=-3.+4.6*tsahilf(i)/break_thresh
                endif
                prob=G01EAF('L', xp, IFAIL)
                nprob=1.d0-prob
                break_prob=break_probs(i)*prob/(break_probs(i)*prob+(1.-apriori_probs(i))*nprob)
                !!         if(ipstart .ne. ipend) then 
                !!           if(break_prob .lt. redfak*break_probs(i)) break_prob=redfak*break_probs(i)
                !!           if(break_prob .lt. 0.1*apriori_probs(i)) break_prob=0.1*apriori_probs(i)
                !!         endif
                break_probs(i)=break_prob
             else
                if(ipend .eq. ipstart) break_probs(i)=0.
             endif
          enddo

       endif

       !! P(A)=aprioriprob
       !! P(B/A)=break_probs
       !! P(\A)=1-aprioriprob
       !! P(B/\A)=1-break_probs
       !! P(A/B)=P(A)P(B/A)/(P(A)P(B/A)+P(\A)P(B/\A))
       !!
       !! P(B/\A)=1/maxlen (uniform distribution of break location if there is no break)

    enddo

    return
  end subroutine bayes_break_simple



  subroutine locate_combined_breaks(rcpara,iname,prob,probbreaks,tsa_of_mean,mean_of_tsa,break_thresh,break_fak,apriori_probs,iter,protunit) !

    implicit none

    type(rasocor_namelist),intent(in) :: rcpara

    integer i,ii,istart,istop,ip,l,protunit,mindist,ihilf(2),ic1,ic2,mmindist
    integer(kind=JPRM) :: probbreaks(rcpara%nmax,rcpara%probmax),midx(rcpara%pmax,rcpara%parmax)
    integer, intent(in) :: iter,iname
    real(kind=JPRM),intent(in) :: prob(rcpara%nmax,rcpara%probmax)
    real(kind=JPRM) :: meanmax(rcpara%probmax),apriori_probs(rcpara%nmax)
    real(kind=JPRM) :: break_thresh(rcpara%probmax),break_fak
    real(kind=JPRM),intent(in) :: tsa_of_mean(rcpara%nmax,rcpara%probmax),mean_of_tsa(rcpara%nmax,rcpara%probmax)
    logical :: radbed(rcpara%nmax)
    meanmax=maxval(prob,1)
    probbreaks=0
    mindist=rcpara%snht_maxlen/4
    ic1=toindex(19730101,rcpara)
    ic2=toindex(19780101,rcpara)


    !! 00h and 12h break probabilities are compared together. Only break probabilities which are
    !! local maxima of both are corrected
    do ip=1,2
       l=0
       do i=1,rcpara%nmax !! only consider breaks well away from end of time series
          radbed(i)= apriori_probs(i) .gt. rcpara%ap_prob_default .or. tsa_of_mean(i,5) .eq. rcpara%miss_val .or. &
               tsa_of_mean(i,5) .gt. tsa_of_mean(i,6) .or. tsa_of_mean(i,6) .lt. rcpara%break_thresh_rad .and. &
               (tsa_of_mean(i,5) .gt. tsa_of_mean(i,12) .or. tsa_of_mean(i,12) .lt. rcpara%break_thresh_rad .or. tsa_of_mean(i,12) .eq. rcpara%miss_val)

          if(prob(i,ip) .gt. 0.7 .and. (prob(i,ip+6) .gt. break_thresh(ip) .or. tsa_of_mean(i,ip+6) .eq. rcpara%miss_val) .and. radbed(i)) then
             !    if(prob(i,ip) .gt. rcpara%break_thresh .and. iname .gt. 50000 .and. iname .lt. 60000 .or. prob(i,ip) .gt. 0.9 .and. (prob(i,ip+6) .gt. break_thresh(ip) .or. tsa_of_mean(i,ip+6) .eq. rcpara%miss_val) .and. radbed(i)) then

             istart=max(1,i-mindist)
             istop=min(rcpara%nmax,i+mindist)
             if(prob(i,ip) .eq. maxval(prob(istart:istop,1:2))) then
                if(l .eq. 0) then
                   l=l+1
                   probbreaks(l,ip)=i
                else
                   if(i-probbreaks(l,ip) .gt. mindist) then 
                      l=l+1
                      probbreaks(l,ip)=i
                   endif
                endif
             endif
          else
          endif
       enddo
    enddo

    !! 00h and 12h break probabilities are compared together. Only break probabilities which are
    !! local maxima of both are corrected
    do ip=3,4
       l=0
       do i=1,rcpara%nmax
          !    radbed(i)= apriori_probs(i) .gt. rcpara%ap_prob_default .or. tsa_of_mean(i,5) .eq. rcpara%miss_val .or. tsa_of_mean(i,5) .gt. tsa_of_mean(i,6)  &
          !.or. tsa_of_mean(i,6) .lt. rcpara%break_thresh_rad .and. (tsa_of_mean(i,5) .gt. tsa_of_mean(i,12) .or. tsa_of_mean(i,12) .lt. rcpara%break_thresh_rad .or. tsa_of_mean(i,12) .eq. rcpara%miss_val)

          if(prob(i,ip) .gt. 0.7 .and. (prob(i,ip+6) .gt. break_thresh(ip)  .or. tsa_of_mean(i,ip+6) .eq. rcpara%miss_val) .and. radbed(i)) then
             !    if(prob(i,ip) .gt. rcpara%break_thresh .and. iname .gt. 50000 .and. iname .lt. 60000 .or. prob(i,ip) .gt. 0.7 .and. (prob(i,ip+6) .gt. break_thresh(ip)  .or. tsa_of_mean(i,ip+6) .eq. rcpara%miss_val) .and. radbed(i)) then
             istart=max(1,i-mindist)
             istop=min(rcpara%nmax,i+mindist)
             if(prob(i,ip) .eq. maxval(prob(istart:istop,3:4))) then
                !! There may be equal values of prob if increment gt 1 has been used in snhteqsamp
                !! Use only first value in this case
                if(l .eq. 0) then
                   l=l+1
                   probbreaks(l,ip)=i
                else
                   if(i-probbreaks(l,ip) .gt. mindist) then 
                      l=l+1
                      probbreaks(l,ip)=i
                   endif
                endif
             endif
          else
          endif
       enddo
    enddo


    !! 00h and 12h break probabilities are compared together. Only break probabilities which are
    !! local maxima of both are corrected
    do ip=7,8
       l=0
       do i=1,rcpara%nmax !! only consider breaks well away from end of time series
          radbed(i)=apriori_probs(i) .gt. rcpara%ap_prob_default .or. tsa_of_mean(i,5) .eq. rcpara%miss_val .or. &
               tsa_of_mean(i,5) .gt. tsa_of_mean(i,6) .or. tsa_of_mean(i,6) .lt. rcpara%break_thresh_rad

          radbed(i)= .true. !!tsa_of_mean(i,5) .eq. rcpara%miss_val .or.  (tsa_of_mean(i,5) .gt. tsa_of_mean(i,12) .or. tsa_of_mean(i,12) .lt. rcpara%break_thresh_rad .or. tsa_of_mean(i,12) .eq. rcpara%miss_val)
          if(prob(i,ip) .gt. 0.7 .and. prob(i,ip-6) .gt. break_thresh(ip) .and. radbed(i)) then
             istart=max(1,i-mindist)
             istop=min(rcpara%nmax,i+mindist)
             if(prob(i,ip) .eq. maxval(prob(istart:istop,7:8))) then
                !! There may be equal values of prob if increment gt 1 has been used in snhteqsamp
                !! Use only first value in this case
                if(l .eq. 0) then
                   l=l+1
                   probbreaks(l,ip)=i
                else
                   if(i-probbreaks(l,ip) .gt. mindist) then 
                      l=l+1
                      probbreaks(l,ip)=i
                   endif
                endif
             endif
          else
          endif
       enddo
    enddo

    !! 00h and 12h break probabilities are compared together. Only break probabilities which are
    !! local maxima of both are corrected
    do ip=9,10
       l=0
       do i=1,rcpara%nmax
          radbed(i)=.true. !!tsa_of_mean(i,5) .eq. rcpara%miss_val .or.  (tsa_of_mean(i,5) .gt. tsa_of_mean(i,12) .or. tsa_of_mean(i,12) .lt. rcpara%break_thresh_rad .or. tsa_of_mean(i,12) .eq. rcpara%miss_val)
          if(prob(i,ip) .gt. 0.7 .and.  prob(i,ip-6) .gt. break_thresh(ip)  .and. radbed(i)) then
             istart=max(1,i-mindist)
             istop=min(rcpara%nmax,i+mindist)
             if(prob(i,ip) .eq. maxval(prob(istart:istop,9:10))) then
                !! There may be equal values of prob if increment gt 1 has been used in snhteqsamp
                !! Use only first value in this case
                if(l .eq. 0) then
                   l=l+1
                   probbreaks(l,ip)=i
                else
                   if(i-probbreaks(l,ip) .gt. mindist) then 
                      l=l+1
                      probbreaks(l,ip)=i
                   endif
                endif
             endif
          else
          endif
       enddo
    enddo

    !! breaks in radiation errors are treated separately. They have priority over breaks
    !! in either 00h bg-obs and 12h bg-obs
    do ip=5,5
       l=0
       do i=1,rcpara%nmax

          radbed(i)= tsa_of_mean(i,5) .gt. tsa_of_mean(i,6) .or. (tsa_of_mean(i,5) .gt. tsa_of_mean(i,12) .or. tsa_of_mean(i,12) .eq. rcpara%miss_val)

          if(prob(i,ip) .gt. break_thresh(ip) .and. prob(i,ip) .gt. break_fak*maxval(meanmax(5:5)) .and. radbed(i)) then
             mmindist=mindist
             if(tsa_of_mean(i,5) .gt. 100) mmindist=mindist/2
             istart=max(1,i-mmindist)
             istop=min(rcpara%nmax,i+mmindist)
             if(prob(i,ip) .eq. maxval(prob(istart:istop,5:5))) then
                istart=max(1,i-rcpara%snht_maxlen/3)
                istop=min(rcpara%nmax,i+rcpara%snht_maxlen/3)
                if(5*tsa_of_mean(i,ip) .gt. maxval(tsa_of_mean(istart:istop,1:4))) then
                   ii=i
                else
                   ihilf=maxloc(prob(istart:istop,1:4))
!!$ call omp_set_lock(omp_lp)
                   !          write(*,*) 'loc: ',iname,i,ihilf
!!$ call omp_unset_lock(omp_lp)
                   ii=istart+ihilf(1)-1
                endif
                !! There may be equal values of prob if increment gt 1 has been used in snhteqsamp
                !! Use only first value in this case
                if(l .eq. 0) then
                   l=l+1
                   probbreaks(l,ip)=ii
                else
                   if(i-probbreaks(l,ip) .gt. mmindist) then 
                      l=l+1
                      probbreaks(l,ip)=ii
                   endif
                endif
             else
                !!        if(5*tsa_of_mean(i,ip) .le. 0*maxval(tsa_of_mean(istart:istop,1:4))) then
                !!          write(*,*) 'rad_break overruled'
                !!        endif
             endif
          else
          endif
       enddo
    enddo

    !! breaks in bg 12h-00h difference are likely spurious. They are calculated only for diagnostic purposes.
    do ip=6,6
       l=0
       do i=1,rcpara%nmax
          if(tsa_of_mean(i,6) .gt. rcpara%break_thresh_rad) then
             istart=max(1,i-mindist)
             istop=min(rcpara%nmax,i+mindist)
             if(tsa_of_mean(i,6) .eq. maxval(tsa_of_mean(istart:istop,6))) then
                !! There may be equal values of prob if increment gt 1 has been used in snhteqsamp
                !! Use only first value in this case
                if(l .eq. 0) then
                   l=l+1
                   probbreaks(l,ip)=i
                else
                   if(i-probbreaks(l,ip) .gt. mindist) then 
                      l=l+1
                      probbreaks(l,ip)=i
                   endif
                endif
             endif
          endif
       enddo
    enddo

    if(l .gt. 0) then
       !!  write(*,*) 'spurious bg breaks',probbreaks(1:l,6),rcpara%year(probbreaks(1:l,6))*10000+rcpara%month(probbreaks(1:l,6))*100+rcpara%day(probbreaks(1:l,6))
    endif
    if(iname .eq. 59981) then
       do i=1,12
          write(*,'(A8,10I6)') 'locate:',probbreaks(1:10,i)
       enddo
    endif
    return
  end subroutine locate_combined_breaks


  subroutine select_breaks2(rcpara,meanprobbreaks,probmeanbreaks,cardsmeta,era40meta,chosenbreaks,protunit) !

    implicit none

    type(rasocor_namelist),intent(in) :: rcpara

    integer select_mode,protunit,i,k,j,istart,istop,ic,ii
    integer(kind=JPRM) :: meanprobbreaks(rcpara%nmax,rcpara%probmax),probmeanbreaks(rcpara%nmax,rcpara%probmax),chosenbreaks(rcpara%nmax)
    real(kind=JPRM) :: break_thresh(rcpara%probmax)
    integer cprobmean(rcpara%probmax),cmeanprob(rcpara%probmax),kstart,index(9)
    integer(kind=JPRM) :: cardsmeta(rcpara%nmax,1,rcpara%nmeta),era40meta(rcpara%nmax),mindist
    logical lnew,probmeanmask(rcpara%nmax),meanprobmask(rcpara%nmax)


    !! if bg has radiation error do not select breaks
    mindist=rcpara%snht_maxlen/4
    meanprobmask=.false.
    probmeanmask=.false.

    probmeanmask=probmeanmask .or. meanprobmask
    meanprobmask=probmeanmask

    kstart=1
    if(chosenbreaks(1).ne. 0)  kstart=2 !! if end of series is corrected above
    ic=1
    index=(/5,2,1,4,3,8,7,10,9/)
    do ii=1,9 !!rcpara%probmax
       i=index(ii)
       cmeanprob(i)=count(meanprobbreaks(:,i) .gt. 0)
       do j=1,cmeanprob(i)
          lnew=.true.
          do k=kstart,ic
             istart=chosenbreaks(k)-mindist
             istop=chosenbreaks(k)+mindist
             if(meanprobbreaks(j,i) .lt. istop .and. meanprobbreaks(j,i) .gt. istart) lnew=.false.
          enddo
          if(lnew .and. .not. meanprobmask(meanprobbreaks(j,i)) .and. meanprobbreaks(j,i) .ne. 1) then
             chosenbreaks(ic+1)=meanprobbreaks(j,i)
             ic=ic+1
          endif
       enddo
       cprobmean(i)=count(probmeanbreaks(:,i) .gt. 0)
       do j=1,cprobmean(i)
          lnew=.true.
          do k=kstart,ic
             istart=chosenbreaks(k)-mindist
             istop=chosenbreaks(k)+mindist
             if(probmeanbreaks(j,i) .lt. istop .and. probmeanbreaks(j,i) .gt. istart) lnew=.false.
          enddo
          if(lnew .and. .not. probmeanmask(probmeanbreaks(j,i)) .and. probmeanbreaks(j,i) .ne. 1) then 
             chosenbreaks(ic+1)=probmeanbreaks(j,i) 
             ic=ic+1
          endif
       enddo

    enddo

    if(chosenbreaks(1) .eq. 0) then 
       chosenbreaks(1:ic-1)=chosenbreaks(2:ic)
       chosenbreaks(ic)=0
    endif
    !!write(*,*) ic-1,'selected:',chosenbreaks(1:ic-1)

    return
  end subroutine select_breaks2


  subroutine calc_breakstats(rcpara,chosenbreaks,tm,tfgm,stm,stfgm,rasobreaks,rasobreakuncertainties,plus,minus,splus,sminus,rplus,rminus,plusmean,minusmean,protunit)

    implicit none

    type(rasocor_namelist),intent(in) :: rcpara

    integer protunit
    integer ipar,ib,ip,bcount,pindex,mindex
    integer(kind=JPRM) :: chosenbreaks(rcpara%nmax)
    real(kind=JPRM) :: tm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax),stm(rcpara%nmax,rcpara%pmax,rcpara%parmax),stfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM),intent(inout) :: rasobreaks(rcpara%nmax,rcpara%pmax,rcpara%parmax),rasobreakuncertainties(rcpara%nmax,rcpara%pmax,rcpara%parmax)

    real(kind=JPRM) :: plus(rcpara%nmax,rcpara%pmax,rcpara%parmax),minus(rcpara%nmax,rcpara%pmax,rcpara%parmax),prms(rcpara%nmax,rcpara%pmax,rcpara%parmax),mrms(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM) :: splus(rcpara%nmax,rcpara%pmax,rcpara%parmax),sminus(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM) :: rplus(rcpara%nmax,rcpara%pmax,1),rminus(rcpara%nmax,rcpara%pmax,1)
    real(kind=JPRM) :: plusmean(rcpara%nmax,rcpara%parmax),minusmean(rcpara%nmax,rcpara%parmax)
    real(kind=JPRM) :: breakprofile(rcpara%pmax,rcpara%parmax),f_val(rcpara%parmax)
    integer(kind=JPRM) :: cardsmeta(rcpara%nmax,1,rcpara%nmeta),era40meta(rcpara%nmax)
    real(kind=JPRM) :: sig

    bcount=count(chosenbreaks .gt. 0)

    do ib=1,bcount

       breakprofile=0.
       do ipar=1,rcpara%parmax
          do ip=1,rcpara%pmax
             if(count(plus(1:chosenbreaks(ib),ip,ipar) .ne. rcpara%miss_val) .gt. 0 .and. count(minus(chosenbreaks(ib):rcpara%nmax,ip,ipar) .ne. rcpara%miss_val) .gt. 0) then
                if(plus(chosenbreaks(ib),ip,ipar) .ne. 0. .or. minus(chosenbreaks(ib),ip,ipar) .ne. 0.) then 
                   pindex=chosenbreaks(ib)
                   do while(plus(pindex,ip,ipar) .eq. rcpara%miss_val .and. pindex .gt. 1)
                      pindex=pindex-1
                   enddo
                   mindex=chosenbreaks(ib)
                   do while(minus(mindex,ip,ipar) .eq. rcpara%miss_val .and. mindex .lt. rcpara%nmax)
                      mindex=mindex+1
                   enddo
                   breakprofile(ip,ipar)=plus(pindex,ip,ipar)-minus(mindex,ip,ipar)
                   rasobreaks(1:chosenbreaks(ib),ip,ipar)=rasobreaks(1:chosenbreaks(ib),ip,ipar)-breakprofile(ip,ipar)

                   sig=((prms(chosenbreaks(ib),ip,ipar)-plus(chosenbreaks(ib),ip,ipar)*plus(chosenbreaks(ib),ip,ipar)+mrms(chosenbreaks(ib),ip,ipar)-minus(chosenbreaks(ib),ip,ipar)*minus(chosenbreaks(ib),ip,ipar))/2.)
                   if(sig .gt. 0.)  sig=sqrt(sig)
                   rasobreakuncertainties(1:chosenbreaks(ib),ip,ipar)=sig
                else
                   !!        rasobreaks(1:chosenbreaks(ib),ip,ipar)=rcpara%miss_val
                   !!        rasobreakuncertainties(1:chosenbreaks(ib),ip,ipar)=rcpara%miss_val
                endif
                if(mindex .ne. chosenbreaks(ib) .or. pindex .ne. chosenbreaks(ib)) then
                   plus(pindex:mindex,ip,ipar)=0.
                   minus(pindex:mindex,ip,ipar)=0.
                endif
             endif
          enddo
       enddo
    enddo

    return

  end subroutine calc_breakstats

  subroutine adjust_series(rcpara,iter,cname,lon,chosenbreaks,ib,attributionfactors,tm,tfgm,stm,stfgm,tgps,cardsmeta,era40meta,eracorrs,rasocorrs,rasobreaks, &
       rasobreakuncertainties,delemask,plus,minus,splus,sminus,rplus,rminus,compplus,compminus,radplus,radminus,bgrplus,bgrminus,plusmean,minusmean,dailycrut2,midx,lasts,gcount,protunit,alt) !

    implicit none

    type(rasocor_namelist),intent(in) :: rcpara

    integer,intent(in) :: protunit,ib,iter
    integer ipar,ip,bcount,ic,im,cib,k,igood,ibad
    integer(kind=JPRM),intent(in) :: chosenbreaks(rcpara%nmax)
    integer(kind=JPRM) :: correctedbreaks(rcpara%nmax)
    integer         :: lasts(rcpara%brmax,rcpara%parmax),midx(rcpara%brmax,rcpara%parmax),gcount(rcpara%parmax),ibegin,l
    real(kind=JPRM) :: stm(rcpara%nmax,rcpara%pmax,rcpara%parmax),stfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tgps(rcpara%nmax,rcpara%pmax,rcpara%parmax),lon
    real(kind=JPRM),intent(inout) :: tfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM),intent(in) :: eracorrs(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM),intent(inout) :: rasocorrs(rcpara%nmax,rcpara%pmax,rcpara%parmax),rasobreaks(rcpara%nmax,rcpara%pmax,rcpara%parmax),rasobreakuncertainties(rcpara%nmax,rcpara%pmax,rcpara%parmax),delemask(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM),intent(in) :: attributionfactors(rcpara%nmax)

    real(kind=JPRM),intent(in) :: plus(rcpara%nmax,rcpara%pmax,rcpara%parmax),minus(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM),intent(in) :: compplus(rcpara%nmax,rcpara%pmax,rcpara%parmax),compminus(rcpara%nmax,rcpara%pmax,rcpara%parmax)

    !! splus is overwritten!
    real(kind=JPRM) :: splus(rcpara%nmax,rcpara%pmax,rcpara%parmax),sminus(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM),intent(in) :: rplus(rcpara%nmax,rcpara%pmax,1),rminus(rcpara%nmax,rcpara%pmax,1)
    real(kind=JPRM),intent(in) :: radplus(rcpara%nmax,rcpara%pmax,1),radminus(rcpara%nmax,rcpara%pmax,1)
    real(kind=JPRM),intent(in) :: bgrplus(rcpara%nmax,rcpara%pmax,1),bgrminus(rcpara%nmax,rcpara%pmax,1)
    real(kind=JPRM),intent(in) :: plusmean(rcpara%nmax,rcpara%parmax),minusmean(rcpara%nmax,rcpara%parmax)
    real(kind=JPRM) :: breakprofile(rcpara%pmax,rcpara%parmax),compprofile(rcpara%pmax,rcpara%parmax),f_val(rcpara%parmax),scale(rcpara%pmax),alt
    real(kind=JPRM) :: bgrbreakprofile(rcpara%pmax,1),brms(rcpara%parmax),bgbrms,rad(rcpara%pmax)
    real(kind=JPRM) :: rbhilf(rcpara%pmax,rcpara%parmax),rbuhilf(rcpara%pmax,rcpara%parmax),rbradhilf(rcpara%pmax,1)
    real(kind=JPRM) :: redfak(rcpara%pmax),dailycrut2(rcpara%nmax,rcpara%parmax)
    integer(kind=JPRM),intent(in) :: cardsmeta(rcpara%nmax,1,rcpara%nmeta),era40meta(rcpara%nmax)
    integer(kind=JPRM) :: iname,sigcount(rcpara%parmax),insigcount(rcpara%parmax),cibm1
    logical :: lsig
    character*6,intent(in) :: cname


    bcount=count(chosenbreaks .gt. 0)
    cib=chosenbreaks(ib)
    if(ib .gt. 1) then
       cibm1=chosenbreaks(ib-1)
    else
       cibm1=chosenbreaks(ib)+rcpara%mean_maxlen
    endif
    correctedbreaks=0
    bgrbreakprofile=rcpara%miss_val

!!$ call omp_set_lock(omp_lp)
    read(cname,'(I6)') iname
!!$ call omp_unset_lock(omp_lp)
    ic=0
    !!do ib=1,bcount

    !!attributionfactors(ib)=1.0

    rbhilf=0.
    rbuhilf=0.
    do ipar=1,rcpara%parmax
       do ip=1,rcpara%pmax
          if(sminus(cib,ip,ipar) .ne. rcpara%miss_val .and. splus(cib,ip,ipar) .ne. rcpara%miss_val)  then
             !        if(.false. .and. cib .lt. 9000 .and. ip .lt. 6 .and. (iname .gt. 47400 .and. iname .lt. 48000 .or. iname .gt. 20000 .and. iname .lt. 40000)) then
             !          splus(cib,ip,ipar)=0.8*splus(cib,ip,ipar)
             !          sminus(cib,ip,ipar)=0.8*sminus(cib,ip,ipar)
             !          rbhilf(ip,ipar)=-(splus(cib,ip,ipar)-sminus(cib,ip,ipar))
             !        else
             rbhilf(ip,ipar)=-(splus(cib,ip,ipar)-sminus(cib,ip,ipar))
             !        endif
             if(rminus(cib,ip,1) .ne. rcpara%miss_val .and. rplus(cib,ip,1) .ne. rcpara%miss_val .and. rminus(cib,ip,1)*rplus(cib,ip,1) .ne. 0)  then       
                rbuhilf(ip,ipar)=sqrt((rplus(cib,ip,1)*rplus(cib,ip,1)+rminus(cib,ip,1)*rminus(cib,ip,1))/2.)
             endif

             if(ipar .eq. 1 .and. sminus(cib,ip,2) .ne. rcpara%miss_val .and. splus(cib,ip,2) .ne. rcpara%miss_val) then 
                rbradhilf(ip,1)=(splus(cib,ip,2)-sminus(cib,ip,2))-(splus(cib,ip,ipar)-sminus(cib,ip,ipar))
             endif
          endif
       enddo
    enddo

    call calc_profile(rcpara,iname,cib,splus,sminus,rplus,rminus,breakprofile,f_val,rcpara%smooth_method,rcpara%parmax) !in this file line 9
    where(breakprofile .eq. rcpara%miss_val) breakprofile=0.

    call calc_profile(rcpara,iname,cib,bgrplus,bgrminus,rplus,rminus,bgrbreakprofile,f_val,rcpara%smooth_method,1) !in this file line 9
    where(bgrbreakprofile .eq. rcpara%miss_val) bgrbreakprofile=0.

    call calc_profile(rcpara,iname,cib,compplus,compminus,rplus,rminus,compprofile,f_val,rcpara%smooth_method,rcpara%parmax) !in this file line 9
    where(compprofile .eq. rcpara%miss_val) compprofile=0.

    !  if(iname .gt. 94000 .and. iname .lt. 95000 .or. any(iname .eq. (/91643,91517,91680,96996,96441,96413,89611,89642,89564,89571/)) .and. cib .gt.  toindex(19870101,rcpara) .and. cib .lt.  toindex(19890101,rcpara) ) then
    !   breakprofile(12,1)=0.15
    !  endif
    if (alt>2800.or. (iname>72000 .and. iname<75000 .and. todate(cib,rcpara)<19620101)) then !.or.(iname>89000 .and. iname<90000 .and. todate(cib,rcpara)<19920101)) then
       redfak=(/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.5,0.,0.,0./)
       write(*,*) 'reduce:', alt,iname,todate(cib,rcpara)
    else if (alt>1400 ) then
       redfak=(/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.5,0.,0./)
    else
       redfak=(/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.,0./)
    endif
    breakprofile(15:16,:)=0
    do ipar=1,rcpara%parmax
       if(splus(cib,15,ipar).ne. rcpara%miss_val .and.sminus(cib,15,ipar).ne. rcpara%miss_val)  breakprofile(15,ipar)=-(splus(cib,15,ipar)-sminus(cib,15,ipar))
       if(splus(cib,16,ipar).ne. rcpara%miss_val .and.sminus(cib,16,ipar).ne. rcpara%miss_val)  breakprofile(16,ipar)=-(splus(cib,16,ipar)-sminus(cib,16,ipar))
       if( abs(breakprofile(15,ipar)) .lt. 0.5 .or. sqrt((rplus(cib,15,1)*rplus(cib,15,1)+rminus(cib,15,1)*rminus(cib,15,1))/2.) .gt. abs(breakprofile(15,ipar))  &
            .or. 2.*abs(breakprofile(16,ipar)) .gt. abs(breakprofile(15,ipar)) .or. iter .eq. 1 .and. rcpara%maxiter .eq. 2) then
          do ip=1,rcpara%pmax
             breakprofile(ip,ipar)=breakprofile(ip,ipar)*redfak(ip)
          enddo
       else
!!$ call omp_set_lock(omp_lp)
          write(*,'(A6,I9,A19,3F8.3)') cname,todate(cib,rcpara),' Full surface break',breakprofile(15,ipar),sqrt((rplus(cib,15,1)*rplus(cib,15,1)+rminus(cib,15,1)*rminus(cib,15,1))/2.),breakprofile(16,ipar)
!!$ call omp_unset_lock(omp_lp)
          do ip=1,rcpara%pmax
             ! to activate cru4 uncomment this line and comment the following
             !        if(breakprofile(ip,ipar) .gt. rcpara%miss_val) breakprofile(ip,ipar)=breakprofile(ip,ipar)*redfak(ip)+breakprofile(15,ipar)*(1.-redfak(ip))

             if(breakprofile(ip,ipar) .gt. rcpara%miss_val) breakprofile(ip,ipar)=breakprofile(ip,ipar)*redfak(ip)+breakprofile(ip,ipar)*(1.-redfak(ip))
          enddo

       endif
    enddo

    if(any(breakprofile .ne. 0.) .or. count(f_val .gt. 1.) .gt. 0) then

       !!  if(cib .lt. 0) then 
       scale=rcpara%miss_val
       scale=1.0
       ! do ip=1,-rcpara%pmax
       !   if(breakprofile(ip,1) .gt. rcpara%miss_val .and.  breakprofile(ip,2) .gt. rcpara%miss_val) then
       !     scale(ip)=maxval(abs(breakprofile(ip,:)))
       !   else
       !     if(breakprofile(ip,1) .gt. rcpara%miss_val) scale(ip)=abs(maxval(breakprofile(ip,:)))
       !   endif
       !   if(bgrbreakprofile(ip,1) .gt. rcpara%miss_val .and. scale(ip) .gt. rcpara%miss_val) scale(ip)=1.0-abs(bgrbreakprofile(ip,1))/(scale(ip)+0.50)
       !   if(scale(ip) .lt. 0) scale(ip)=0.
       !enddo
!!$ call omp_set_lock(omp_lp)
       !    write(*,'(A5,I9,A8,2(I3,F9.3),14F9.3)') cname,todate(cib,rcpara),' scale: ',maxloc(scale(1:14)),maxval(scale(1:14)),minloc(scale),minval(scale),scale(1:14)
!!$ call omp_unset_lock(omp_lp)

       !  do ipar=1,rcpara%parmax
       !    where(breakprofile(:,ipar) .ne. rcpara%miss_val)
       !      breakprofile(:,ipar)=breakprofile(:,ipar)*scale
       !    endwhere
       !!    breakprofile(14:16,ipar)=breakprofile(14:16,ipar)/2.
       !  enddo
       !!  endif


       call test_significance(rcpara,iname,cib,breakprofile,bgrbreakprofile,compprofile,rbhilf,rbuhilf,rbradhilf,radplus,radminus,lsig,sigcount,insigcount) !in this file line 1496

       if(lsig) then 


          do ipar=1,rcpara%parmax
09           l=1
             do while(lasts(l,ipar) .gt. cib)
                l=l+1
             enddo
             if(l .eq. gcount(ipar)+1) then
                ibegin=1
             else
                ibegin=lasts(l,ipar)+rcpara%snht_maxlen/2.
             endif
             do ip=1,rcpara%pmax

                if(breakprofile(ip,ipar) .ne. 0.) then

                   if(sminus(cib,ip,ipar) .ne. rcpara%miss_val .and. splus(cib,ip,ipar) .ne. rcpara%miss_val)  then
                      !!           if(sminus(cib,ip,3-ipar) .eq. rcpara%miss_val .or. splus(cib,ip,3-ipar) .eq. rcpara%miss_val .or. &
                      !!              2*abs(rad(ip)) .gt. abs(breakprofile(ip,2)-breakprofile(ip,1)) .and. ip .lt. 16  &
                      !!              .or. abs(rad(ip)) .gt. abs(breakprofile(ip,2)-breakprofile(ip,1)) .or. cibm1 .lt. 0*5800) then

                      rasocorrs(ibegin:cib,ip,ipar)=rasocorrs(ibegin:cib,ip,ipar)-breakprofile(ip,ipar)
                      rasobreaks(ibegin:cib,ip,ipar)=rasobreaks(ibegin:cib,ip,ipar)+rbhilf(ip,ipar)
                      rasobreakuncertainties(ibegin:cib,ip,ipar)=rasobreakuncertainties(ibegin:cib,ip,ipar)+rbuhilf(ip,ipar)

                      where(tfgm(ibegin:cib,ip,ipar) .ne. rcpara%miss_val)
                         tfgm(ibegin:cib,ip,ipar)=tfgm(ibegin:cib,ip,ipar)-breakprofile(ip,ipar)
                      endwhere
                      where(tgps(ibegin:cib,ip,ipar) .ne. rcpara%miss_val)
                         tgps(ibegin:cib,ip,ipar)=tgps(ibegin:cib,ip,ipar)-breakprofile(ip,ipar)
                      endwhere
                      where(stfgm(ibegin:cib,ip,ipar) .ne. rcpara%miss_val)
                         stfgm(ibegin:cib,ip,ipar)=stfgm(ibegin:cib,ip,ipar)-breakprofile(ip,ipar)
                      endwhere

                      !!           else
!!$ call omp_set_lock(omp_lp)
                      !!        write(*,*)  cname,cib,ip,ipar,rad(ip),radplus(cib,ip,1),radminus(cib,ip,1),(breakprofile(ip,2)-breakprofile(ip,1)),'12h-00h of adjustment larger than 12h-00h of original series .. not adjusted'
!!$ call omp_unset_lock(omp_lp)           
                      !!           endif
                      !!          else
                   endif

                else
                endif
             enddo
             where(dailycrut2(ibegin:cib,ipar) .ne. rcpara%miss_val)
                dailycrut2(ibegin:cib,ipar)=dailycrut2(ibegin:cib,ipar)-breakprofile(14,ipar)
             endwhere
          enddo

          write(protunit,'(A6,A7,I10,A10,4I3)') cname,':break ',todate(cib,rcpara),' accepted',sigcount,insigcount

          do ip=1,rcpara%pmax
             write(protunit,'(I3,2F8.2)') ip,breakprofile(ip,:)
          enddo
          ic=ic+1
          correctedbreaks(ic)=cib

       else
          write(protunit,'(A6,A7,I10,A10,4I3)') cname,':break ',todate(cib,rcpara),' rejected',sigcount,insigcount
          do ip=1,rcpara%pmax
             write(protunit,'(I3,2F8.2)') ip,breakprofile(ip,:)
          enddo
       endif
    else
       write(protunit,'(A6,A7,I10,A10,3F7.3)') cname,':break ',todate(cib,rcpara),' rejected',f_val
       do ip=1,rcpara%pmax
          write(protunit,'(I3,2F8.2)') ip,breakprofile(ip,:)
       enddo
    endif
    !!enddo

    !!chosenbreaks=correctedbreaks
    !!print*,chosenbreaks(1:ic)

    return

  end subroutine adjust_series

  subroutine test_significance(rcpara,iname,cb,breakprofile,bgrbreakprofile,compprofile,rasobreak,rasobreakuncertainty,rasobreakrad,radplus,radminus,lsig,sigcount,insigcount) !

    implicit none

    type(rasocor_namelist),intent(in) :: rcpara

    integer i,ip,ipar,sigcount(rcpara%parmax),insigcount(rcpara%parmax),sigcountrad,insigcountrad,iname,cb
    logical lsig,ldebug
    real(kind=JPRM) :: sigs(rcpara%pmax,rcpara%parmax),rbreak,runc
    real(kind=JPRM) :: breakprofile(rcpara%pmax,rcpara%parmax),compprofile(rcpara%pmax,rcpara%parmax),bgrbreakprofile(rcpara%pmax,1),rasobreak(rcpara%pmax,rcpara%parmax),rasobreakrad(rcpara%pmax,1),rasobreakuncertainty(rcpara%pmax,rcpara%parmax),radbreak
    real(kind=JPRM),intent(in) :: radplus(rcpara%nmax,rcpara%pmax,1),radminus(rcpara%nmax,rcpara%pmax,1)

    ldebug=.false.
    lsig=.false.
    do ipar=1,rcpara%parmax
       sigcount(ipar)=0
       insigcount(ipar)=0
       do ip=1,14
          if(rcpara%mweights(ip) .ne. 0 .and. rasobreak(ip,ipar) .ne. 0. .and. breakprofile(ip,ipar) .ne. rcpara%miss_val ) then
             !! 
             !!       if(breakprofile(ip,ipar)*rasobreak(ip,ipar)*rasobreakuncertainty(i,ipar) .gt. 0) then
!!$ call omp_set_lock(omp_lp)

             if(ldebug)  write(*,'(I6,I9,2I3,3F7.2)') iname,todate(cb,rcpara),ip,ipar,abs(breakprofile(ip,ipar)),abs(breakprofile(ip,ipar))-rasobreakuncertainty(ip,ipar)*1.96,abs(bgrbreakprofile(ip,1))

!!$ call omp_unset_lock(omp_lp)
!!!          if(abs(breakprofile(ip,ipar)) .gt. 0.3 .and. abs(breakprofile(ip,ipar))-rasobreakuncertainty(ip,ipar)*1.0 .gt. 0 .and. abs(bgrbreakprofile(ip,1)) .lt. 2*abs(breakprofile(ip,ipar))) then

             if(abs(breakprofile(ip,ipar)) .gt. rcpara%sig_thresh .and. (todate(cb,rcpara) .gt. 19800101 .or. abs(compprofile(ip,ipar)) .gt. rcpara%sig_thresh .or. compprofile(ip,ipar) .eq. 0.) .and. abs(breakprofile(ip,ipar))-rasobreakuncertainty(ip,ipar)*1.96 .gt. 0) then
                !          if(abs(breakprofile(ip,ipar)) .gt. rcpara%sig_thresh .and. abs(breakprofile(ip,ipar))-rasobreakuncertainty(ip,ipar)*1.96 .gt. 0) then
                sigcount(ipar)=sigcount(ipar)+1
                if (ip .gt. 1 .and. abs(breakprofile(ip,ipar)) .gt. 5*rcpara%sig_thresh .and. abs(breakprofile(ip,ipar))-rasobreakuncertainty(ip,ipar)*2.97.gt. 0)        sigcount(ipar)=sigcount(ipar)+1

             else
                insigcount(ipar)=insigcount(ipar)+1
             endif
             !!       endif
          endif
       enddo
       if(ipar .eq. 1) then
          sigcountrad=0
          insigcountrad=0
          do ip=1,12
             if(rasobreak(ip,1) .ne. 0. .and. rasobreak(ip,2) .ne. 0 .and. breakprofile(ip,ipar) .ne. rcpara%miss_val ) then
                !!       if(breakprofile(ip,ipar)*rasobreak(ip,ipar)*rasobreakuncertainty(i,ipar) .gt. 0) then
                rbreak=abs(breakprofile(ip,2)-breakprofile(ip,1))
                if(radplus(cb,ip,1) .gt. rcpara%miss_val .and. radminus(cb,ip,1)  .gt. rcpara%miss_val) then 
                   radbreak=abs(radplus(cb,ip,1)-radminus(cb,ip,1))
                else
                   radbreak=999.
                endif
                runc=sqrt((rasobreakuncertainty(ip,1)**2+rasobreakuncertainty(ip,2)**2)/2.)*1.96
!!$ call omp_set_lock(omp_lp)
                if(ldebug)          write(*,'(I6,I9,I6,5F7.2)') iname,todate(cb,rcpara),ip,rbreak,rbreak-abs(bgrbreakprofile(ip,1)),abs(bgrbreakprofile(ip,1)),radbreak,runc
!!$ call omp_unset_lock(omp_lp)
!!!          if(rbreak.gt. 0.3 .and. rbreak-runc .gt. 0 .and.  2*radbreak .gt. rbreak) then
                if(rbreak.gt. rcpara%sig_thresh_rad .and. rbreak-runc .gt. 0) then
                   !!          if(rbreak.gt. 0.3 .and. rbreak-runc .gt. 0 &
                   !!.and. rbreak-abs(bgrbreakprofile(ip,1)) .gt. 0 .and. 2*radbreak .gt. rbreak) then
                   !!          if(rbreak.gt. 0.3 .and. rbreak-runc .gt. 0) then
                   sigcountrad=sigcountrad+1
                else
                   insigcountrad=insigcountrad+1
                endif
                !!       endif
             endif
          enddo
       endif
       lsig= lsig .or. (sigcount(ipar) .gt. 1 .or. sigcountrad .gt. 1)
       if(iname .gt. 50000 .and. iname .lt. 60000) lsig= lsig .or. (sigcount(ipar) .gt. 1) .or. (sigcountrad .gt. 0 .or. insigcountrad .eq. 0)
    enddo

    if(sigcount(1) .eq. 0 .and. sigcount(2) .eq. 0  .and. sigcountrad .eq. 0) lsig=.false.

    !!lsig=.true.
!!$ call omp_set_lock(omp_lp)
    write(*,'(A5,I6,I9,L2,6I3)') 'iname',iname,todate(cb,rcpara),lsig, sigcount,insigcount,sigcountrad,insigcountrad
!!$ call omp_unset_lock(omp_lp)

    return

  end subroutine test_significance


  !! average_rad calculates the average of 2 00h or 2 12h observations, to be compared with 
  !! observation between. Resulting average is better suited for finding radiation errors
  !! if only one obs is available, this one is returned as average
  !!
  !! L. Haimberger, 08072004
  !!
  !! ni=number of days
  !! t=unaveraged temperature series
  !! tav=averaged temperature series
  !! itime == 0 or 12 ; If 0 then average is made with i-1,i; If 12 then average is made with i,i+1
  subroutine average_rad(ni,t,tav,itime,miss_val) !

    implicit none

    integer ni,itime,i
    real(kind=JPRM) :: t(ni),tav(ni),miss_val

    tav=miss_val
    if(itime .eq. 0) then
       do i=2,ni
          if(t(i) .ne. miss_val .and. t(i-1) .ne. miss_val) then
             tav(i)=0.5*(t(i)+t(i-1))
          else
             if(t(i) .ne. miss_val) tav(i)=t(i)
             if(t(i-1) .ne. miss_val) tav(i)=t(i-1)
          endif
       enddo
    else
       do i=1,ni-1
          if(t(i) .ne. miss_val .and. t(i+1) .ne. miss_val) then
             tav(i)=0.5*(t(i)+t(i+1))
          else
             if(t(i) .ne. miss_val) tav(i)=t(i)
             if(t(i+1) .ne. miss_val) tav(i)=t(i+1)
          endif
       enddo
    endif

    return
  end subroutine average_rad


  subroutine apriori_prob(rcpara,cname,cardsmeta,era40meta,apriori_probs,apriori_probs_rad,iter) !

    implicit none

    type(rasocor_namelist),intent(in) :: rcpara

    integer iname,iter,i,j,mindist,istart,istop,ilatest,imax,jimax
    integer(kind=JPRM) :: cardsmeta(rcpara%nmax,1,rcpara%nmeta),era40meta(rcpara%nmax)
    real(kind=JPRM) :: apriori_probs(rcpara%nmax),apriori_probs_rad(rcpara%nmax),aphilf(rcpara%nmax),aphilf_rad(rcpara%nmax),rmax,fak(rcpara%snht_maxlen),half
    character*6 cname

    apriori_probs=rcpara%ap_prob_default
    !!apriori_probs(16072:rcpara%nmax)=0.
!!$ call omp_set_lock(omp_lp)
    read(cname,'(I6.6)') iname
!!$ call omp_unset_lock(omp_lp) 

    do j=1,rcpara%nmax
       if(era40meta(j) .gt. 0) then
          istart=maxval((/1,j-365/))
          istop=minval((/j+365,rcpara%nmax/))
          do i=istart,istop
             mindist=366
             if(apriori_probs(i) .lt. 0.1) then
                mindist=minval((/abs(j-i),mindist/))
                apriori_probs(i)=apriori_probs(i)*(1.-exp(-(mindist/90.)*(mindist/90.)))
             endif
          enddo
       endif
    enddo

    where(cardsmeta(:,1,2) .gt. 0 .and. apriori_probs .lt. 0.05) apriori_probs=0.96 !! sonde change
    where(cardsmeta(:,1,3) .gt. 0 .and. apriori_probs .lt. 0.05) apriori_probs=0.96 !! 0.5 !! radcor


!!!!!!!!!!!!! correct false metadata
    if(cname .eq. '70350') apriori_probs(11000:12000)=rcpara%ap_prob_default !!Kodiak
    if(cname .eq. '70316') apriori_probs(11000:12000)=rcpara%ap_prob_default !!Kodiak
    !!if(cname .eq. '22217') apriori_probs(15700:15800)=ap_prob_default*10 !Kandalaksa
    if(cname .eq. '02527') apriori_probs(4000:4030)=0.99
    if(cname .eq. '61902') apriori_probs(11323:11433)=0.99
    if(iname .gt. 3000 .and. iname .lt. 4000) apriori_probs(toindex(19710101,rcpara))=0.96
    if(iname .gt. 3000 .and. iname .lt. 4000) apriori_probs(toindex(19600101,rcpara))=0.96
    if(iname .eq. 3322 ) apriori_probs_rad(toindex(19660101,rcpara))=0.99
    if(iname .eq. 6447 ) apriori_probs_rad(toindex(19860101,rcpara):toindex(19890101,rcpara))=0.1
    if(iname .eq. 6447 ) apriori_probs(toindex(19860101,rcpara):toindex(19890101,rcpara))=0.1
    if(iname .gt. 4200 .and. iname .lt. 4400) apriori_probs(toindex(19690101,rcpara))=0.96
    if(iname .eq. 3496) apriori_probs(toindex(19710101,rcpara))=rcpara%ap_prob_default
    if(iname .gt. 20000 .and. iname .lt. 40000) where(apriori_probs(4400:4750) .lt. 0.02) apriori_probs(4400:4750)=0.1
    if(iname .eq. 34560) apriori_probs(toindex(19650101,rcpara))=0.96
    if(iname .eq. 35361) apriori_probs(toindex(19740101,rcpara))=0.96
    if(iname .eq. 10238) apriori_probs(toindex(19901001,rcpara))=0.96
    if(iname .eq. 94975) apriori_probs(toindex(19870101,rcpara))=0.96
    if(iname .eq. 47138) apriori_probs(toindex(19790101,rcpara))=0.96
    if(iname .eq. 47058) then
       apriori_probs(toindex(20090101,rcpara))=0.96
       apriori_probs(toindex(19690101,rcpara))=0.96
    endif
    if(iname .eq. 91958) then
       apriori_probs(toindex(19700101,rcpara))=0.96
       apriori_probs(toindex(19860101,rcpara))=0.96
    endif
    if(iname .gt. 91000 .and. iname .lt. 92000) then
       !   apriori_probs(toindex(19800101,rcpara):toindex(19840101,rcpara))=0
    endif
    if(iname .eq. 94120) apriori_probs(toindex(19630101,rcpara))=0.96
    if(iname .gt. 7000 .and. iname .lt. 7600) then
       apriori_probs(5420:5450)=rcpara%ap_prob_default+0.9
       apriori_probs(5450:5510+366)=rcpara%ap_prob_default

       !!  apriori_probs(4384:4384+366)=ap_prob_default
       !!  apriori_probs(4384+182)=0.8
    endif
    if(iname .eq. 7145 ) apriori_probs(toindex(19801001,rcpara))=0.96
    if(iname .gt. 7000 .and. iname .lt. 8000) apriori_probs(14780)=0.96
    if(iname .gt. 70000 .and. iname .lt. 71000) then
       apriori_probs_rad(toindex(19811001,rcpara))=0.96
       apriori_probs_rad(toindex(19830601,rcpara))=0.96
    endif


    !! sonde change from GTS reports are considered certain
    istart=toindex(19900101,rcpara)
    if(iname .eq. 70219) istart=toindex(19960101,rcpara)
    istop=0
    if(any(cardsmeta(istart:rcpara%nmax,1,7) .gt. 0)) then
       loop:  do i=istart,rcpara%nmax
          if(cardsmeta(i,1,7) .gt. 0) then
             istop=i
             exit loop
          endif
       enddo loop
       if(istop .gt. 0) then
          !    where(cardsmeta(istop:rcpara%nmax,1,3) .eq. 0)
          !      apriori_probs(istop:rcpara%nmax)=0.001
          !    endwhere
       endif
    endif


    apriori_probs_rad=apriori_probs
    where(apriori_probs_rad .lt. 0.1)
       apriori_probs_rad=0.1
    endwhere
!!$if(iname .eq. 12120) then
!!$  apriori_probs_rad(toindex(19860101,rcpara):toindex(19940101,rcpara))=0.01
!!$  apriori_probs_rad(toindex(19880101,rcpara))=0.96
!!$  apriori_probs_rad(toindex(19920101,rcpara))=0.96
!!$  apriori_probs(toindex(19860101,rcpara):toindex(19940101,rcpara))=0.1*rcpara%ap_prob_default
!!$  apriori_probs(toindex(19880101,rcpara))=0.96
!!$  apriori_probs(toindex(19920101,rcpara))=0.96
!!$endif
    if(iname .gt. 94000 .and. iname .lt. 95000) then
       apriori_probs(istart:rcpara%nmax)=0.005
       apriori_probs_rad(istart:rcpara%nmax)=0.005
    endif


    if(iname .gt. 68000 .and. iname .lt. 69000) then
       apriori_probs(toindex(19690101,rcpara):toindex(19710101,rcpara))=rcpara%ap_prob_default
       apriori_probs_rad(toindex(19690101,rcpara):toindex(19710101,rcpara))=rcpara%ap_prob_default
       apriori_probs(toindex(19710101,rcpara):toindex(19720101,rcpara))=0.9
    endif

    !! sonde change from GTS reports are considered certain
    if(any(cardsmeta(:,1,7) .gt. 0)) then
       where(cardsmeta(:,1,7) .gt. 0) 
          apriori_probs=0.96
          apriori_probs_rad=0.96
       endwhere

    endif

    !!if(iname .gt. 3000) then
    apriori_probs(toindex(19990901,rcpara):toindex(20000101,rcpara))=rcpara%ap_prob_default !! don't trust metad
    apriori_probs_rad(toindex(19990901,rcpara):toindex(20000101,rcpara))=rcpara%ap_prob_default !! don't trust metadata from 199910, seem to be generic

    half=20_JPRM
    if(rcpara%use_meta .eq. 2) half=200_JPRM

    aphilf=apriori_probs
    aphilf_rad=apriori_probs_rad
    fak=0.
    write(*,*) rcpara%snht_maxlen
    do i=1,rcpara%snht_maxlen/2+1
       fak(i)=exp(-((i-1)*half/rcpara%snht_maxlen)**2)
    enddo
    !! broaden metadata peaks
    do i=1,rcpara%nmax

       if(apriori_probs(i) .gt. 0.1) then
          istart=maxval((/1,i-rcpara%snht_maxlen/4/))
          istop=minval((/rcpara%nmax,i+rcpara%snht_maxlen/4/))
          imax=istart
          rmax=-1.0
          do j=istart,istop
             if(aphilf(j) .gt. rmax) then
                rmax=aphilf(j)
                imax=j
             endif
             if(imax .lt. i .and. apriori_probs(j) .eq. rmax) imax=j
          enddo
          do j=istart,istop
             jimax=abs(j-imax)
             if(jimax .eq. 0) jimax=1
             if(jimax .gt. rcpara%snht_maxlen/2) jimax=rcpara%snht_maxlen/2
             apriori_probs(j)=maxval((/apriori_probs(j),aphilf(imax)*fak(jimax)/))
             apriori_probs_rad(j)=maxval((/apriori_probs_rad(j),aphilf_rad(imax)*fak(jimax)/))
          enddo
       endif
    enddo


    if(iname .eq. 34247) apriori_probs_rad(toindex(19681001,rcpara):toindex(19681201,rcpara))=0.999
    if(iname .eq. 34247) apriori_probs_rad(toindex(19690101,rcpara):toindex(19710101,rcpara))=0.0
    if(iname .eq. 34247) apriori_probs(toindex(19690101,rcpara):toindex(19681201,rcpara))=0.999
    if(iname .eq. 34247) apriori_probs(toindex(19690101,rcpara):toindex(19710101,rcpara))=0.0

    if(iname .gt. 72000 .and. iname .lt. 75000) then 
       where(apriori_probs(2557:11323) .le. rcpara%ap_prob_default .and.  apriori_probs_rad(2557:11323) .le. 0.1)
          apriori_probs(2557:11323)=0.
          apriori_probs_rad(2557:11323)=0.
       endwhere
    endif

    if(iname .eq. 76654) then
       apriori_probs(toindex(19910101,rcpara):toindex(19940101,rcpara))=0.0
       apriori_probs(toindex(19831201,rcpara):toindex(19840201,rcpara))=0.9
       apriori_probs_rad(toindex(19910101,rcpara):toindex(19940101,rcpara))=0.0
    endif

    if(iname .eq. 85543) then
       apriori_probs(toindex(19770101,rcpara):toindex(19800101,rcpara))=0.0
       apriori_probs_rad(toindex(19860101,rcpara):toindex(19880501,rcpara))=0.0
       apriori_probs(toindex(19860101,rcpara):toindex(19880501,rcpara))=0.0
    endif


    if(rcpara%use_meta .eq. 0) then
       apriori_probs=rcpara%ap_prob_default
       apriori_probs_rad=rcpara%ap_prob_default*5
       do j=1,rcpara%nmax
          if(era40meta(j) .gt. 0) then
             istart=maxval((/1,j-365/))
             istop=minval((/j+365,rcpara%nmax/))
             do i=istart,istop
                mindist=366
                if(apriori_probs(i) .lt. 0.1) then
                   mindist=minval((/abs(j-i),mindist/))
                   apriori_probs(i)=apriori_probs(i)*(1.-exp(-(mindist/180.)*(mindist/180.)))
                endif
             enddo
          endif
       enddo
    endif
where(apriori_probs .le. rcpara%ap_prob_default)
apriori_probs=rcpara%ap_prob_default
apriori_probs_rad=0.1
endwhere
if(iname .gt. 47800 .and. iname .lt. 48000 .or. iname .gt. 70000 .and. iname .lt. 75000) then
apriori_probs=rcpara%ap_prob_default
apriori_probs_rad=0.1
endif
    if(rcpara%use_meta .eq. 2) then
       where(apriori_probs .gt. 0.45)
          apriori_probs=0.999
          apriori_probs_rad=0.999
       elsewhere
          apriori_probs=0.
          apriori_probs_rad=0.
       endwhere
    endif

    if(iname .eq. 10238) apriori_probs(toindex(19901001,rcpara))=0.96
    if(iname .eq. 81405) apriori_probs(toindex(19870101,rcpara):toindex(19960101,rcpara))=rcpara%ap_prob_default
    if(iname .eq. 20744) apriori_probs(toindex(20020101,rcpara):toindex(20060101,rcpara))=0.0
    if(iname .eq. 54135) apriori_probs(toindex(19670101,rcpara):toindex(19710101,rcpara))=0.
    if(iname .eq. 54135) apriori_probs_rad(toindex(19670101,rcpara):toindex(19710101,rcpara))=0.
    if(iname .eq. 47058) then
        apriori_probs(toindex(19770201,rcpara):toindex(19840201,rcpara))=0.0
        apriori_probs_rad(toindex(19770201,rcpara):toindex(19840201,rcpara))=0.0
    endif

    return
  end subroutine apriori_prob

  subroutine apriori_prob_schroeder(rcpara,cname,cardsmeta,era40meta,apriori_probs,apriori_probs_rad,iter) !

    implicit none

    type(rasocor_namelist),intent(in) :: rcpara

    integer iname,iter,i,j,mindist,istart,istop,ilatest,imax,jimax
    integer(kind=JPRM) :: cardsmeta(rcpara%nmax),era40meta(rcpara%nmax),pdist(rcpara%nmax),mdist(rcpara%nmax)
    real(kind=JPRM) :: apriori_probs(rcpara%nmax),apriori_probs_rad(rcpara%nmax),aphilf(rcpara%nmax),aphilf_rad(rcpara%nmax),rmax,fak(rcpara%snht_maxlen),half
    character*6 cname

    read(cname,'(I6.6)') iname

    !if( .false.) then 

    !apriori_probs=rcpara%ap_prob_default
    !!apriori_probs(16072:rcpara%nmax)=0.
!!$ call omp_set_lock(omp_lp)
!!$ call omp_unset_lock(omp_lp) 

    do j=1,rcpara%nmax
       if(era40meta(j) .gt. 0) then
          istart=maxval((/1,j-365/))
          istop=minval((/j+365,rcpara%nmax/))
          do i=istart,istop
             mindist=366
             if(apriori_probs(i) .lt. 0.1) then
                mindist=minval((/abs(j-i),mindist/))
                apriori_probs(i)=apriori_probs(i)*(1.-exp(-(mindist/180.)*(mindist/180.)))
             endif
          enddo
       endif
    enddo
    !endif

    mdist(1)=1000
    pdist(rcpara%nmax)=0
    do i=2,rcpara%nmax
       mdist(i)=mdist(i-1)+1
       if(cardsmeta(i-1) .ne. 0) mdist(i)=0
    enddo
    do i=rcpara%nmax-1,1,-1
       pdist(i)=pdist(i+1)+1
       if(cardsmeta(i+1) .ne. 0) pdist(i)=0
    enddo
    do i=2,rcpara%nmax
       if(cardsmeta(i) .gt. 0 .and. mdist(i) .gt. 92 .and. apriori_probs(i) .lt. 0.05)  apriori_probs(i)=0.96
    enddo
    !where(cardsmeta .gt. 0 .and. apriori_probs .lt. 0.05) apriori_probs=0.96 !! sonde change


    !if( .false.) then 

    apriori_probs_rad=apriori_probs
    where(apriori_probs_rad .lt. 0.1)
       apriori_probs_rad=0.1
    endwhere

    !endif


    half=20_JPRM
    if(rcpara%use_meta .eq. 2) half=200_JPRM

    aphilf=apriori_probs
    aphilf_rad=apriori_probs_rad
    fak=0.
    write(*,*) rcpara%snht_maxlen
    do i=1,rcpara%snht_maxlen/2+1
       fak(i)=exp(-((i-1)*half/rcpara%snht_maxlen)**2)
    enddo
    !! broaden metadata peaks
    do i=1,rcpara%nmax

       if(apriori_probs(i) .gt. 0.1) then
          istart=maxval((/1,i-rcpara%snht_maxlen/4/))
          istop=minval((/rcpara%nmax,i+rcpara%snht_maxlen/4/))
          imax=istart
          rmax=-1.0
          do j=istart,istop
             if(aphilf(j) .gt. rmax) then
                rmax=aphilf(j)
                imax=j
             endif
             if(imax .lt. i .and. apriori_probs(j) .eq. rmax) imax=j
          enddo
          do j=istart,istop
             jimax=abs(j-imax)
             if(jimax .eq. 0) jimax=1
             if(jimax .gt. rcpara%snht_maxlen/2) jimax=rcpara%snht_maxlen/2
             apriori_probs(j)=maxval((/apriori_probs(j),aphilf(imax)*fak(jimax)/))
             apriori_probs_rad(j)=maxval((/apriori_probs_rad(j),aphilf_rad(imax)*fak(jimax)/))
          enddo
       endif
    enddo



    if(rcpara%use_meta .eq. 0) then
       apriori_probs=rcpara%ap_prob_default
       apriori_probs_rad=rcpara%ap_prob_default*5
       do j=1,rcpara%nmax
          if(era40meta(j) .gt. 0) then
             istart=maxval((/1,j-365/))
             istop=minval((/j+365,rcpara%nmax/))
             do i=istart,istop
                mindist=366
                if(apriori_probs(i) .lt. 0.1) then
                   mindist=minval((/abs(j-i),mindist/))
                   apriori_probs(i)=apriori_probs(i)*(1.-exp(-(mindist/180.)*(mindist/180.)))
                endif
             enddo
          endif
       enddo
    endif

    if(iname .eq. 10238) apriori_probs(toindex(19921001,rcpara))=0.96
    if(iname .gt. 55000 .and. iname .lt. 56300) apriori_probs(toindex(19780101,rcpara):toindex(19800101,rcpara)) =0
    if(iname .gt. 55000 .and. iname .lt. 56300) apriori_probs_rad(toindex(19780101,rcpara):toindex(19800101,rcpara)) =0

    if(rcpara%use_meta .gt. 1 .and. (iname .gt. 70000 .and. iname .lt. 75000 .or. iname .gt. 91000 .and. iname .lt. 91300)) then
       where(apriori_probs .gt. 2*rcpara%ap_prob_default)
          apriori_probs=0.999
          apriori_probs_rad=0.999
       elsewhere
          apriori_probs=0.005
          apriori_probs_rad=0.02
       endwhere
    endif
    if(iname .eq. 20744) apriori_probs(toindex(20020101,rcpara):toindex(20120101,rcpara))=0.0
    if(iname .eq. 47058) apriori_probs(toindex(19730101,rcpara):toindex(19950101,rcpara))=0.0


    return
  end subroutine apriori_prob_schroeder


  subroutine bg_corrx(bgcorrs,bgcorrdates,tm,tfgm,eracorrs,istat,rcpara,wmostats)

    type(rasocor_namelist),intent(in) :: rcpara
    integer istart,istop,i,ip,ipar,iunit,ips,ios,cstart,cstop,istat,wmostats
    real(kind=JPRM) :: bgcorrs(wmostats,rcpara%pmax,rcpara%parmax,rcpara%bgbrmax)
    integer(kind=JPRM) :: bgcorrdates(2,rcpara%bgbrmax)

    real(kind=JPRM) :: tfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tm(rcpara%nmax,rcpara%pmax,rcpara%parmax),corr(rcpara%parmax)
    real(kind=JPRM),intent(out) :: eracorrs(rcpara%nmax,rcpara%pmax,rcpara%parmax)

    eracorrs=0.
    j=1
    do while(j .le. rcpara%bgbrmax .and. bgcorrs(istat,1,1,j) .ne. rcpara%miss_val)

       istart=1
       do i=1,rcpara%nmax
          if(rcpara%year(i)*10000+rcpara%month(i)*100+rcpara%day(i) .eq. bgcorrdates(1,j)) istart=i
          if(rcpara%year(i)*10000+rcpara%month(i)*100+rcpara%day(i) .eq. bgcorrdates(2,j)) istop=i
       enddo

       do it=istart,istop

          eracorrs(it,:,:)=eracorrs(it,:,:)+bgcorrs(istat,:,:,j)

       enddo

       j=j+1

    enddo

    where(tfgm .ne. rcpara%miss_val) 
       tfgm=tfgm+eracorrs
    endwhere

    return

  end subroutine bg_corrx

  subroutine adjust_series_comp(rcpara,iter,cname,lon, chosenbreaks,ib,attributionfactors,tm,tfgm,stm,stfgm,tgps,cardsmeta,era40meta,eracorrs,rasocorrs,rasobreaks,&
       rasobreakuncertainties,delemask,plus,minus,splus,sminus,rplus,rminus,compplus,compminus,radplus,radminus,bgrplus,bgrminus,plusmean,minusmean,dailycrut2,midx,lasts,gcount,protunit,alt) !

    implicit none

    type(rasocor_namelist),intent(in) :: rcpara

    integer,intent(in) :: protunit,ib,iter
    integer ipar,ip,bcount,ic,im,cib,k,igood,ibad
    integer(kind=JPRM),intent(in) :: chosenbreaks(rcpara%nmax)
    integer(kind=JPRM) :: correctedbreaks(rcpara%nmax)
    integer         :: lasts(rcpara%brmax,rcpara%parmax),midx(rcpara%brmax,rcpara%parmax),gcount(rcpara%parmax),ibegin,l
    real(kind=JPRM) :: stm(rcpara%nmax,rcpara%pmax,rcpara%parmax),stfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tgps(rcpara%nmax,rcpara%pmax,rcpara%parmax),lon
    real(kind=JPRM),intent(inout) :: tfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM),intent(in) :: eracorrs(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM),intent(inout) :: rasocorrs(rcpara%nmax,rcpara%pmax,rcpara%parmax),rasobreaks(rcpara%nmax,rcpara%pmax,rcpara%parmax),rasobreakuncertainties(rcpara%nmax,rcpara%pmax,rcpara%parmax),delemask(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM),intent(in) :: attributionfactors(rcpara%nmax)

    real(kind=JPRM),intent(in) :: plus(rcpara%nmax,rcpara%pmax,rcpara%parmax),minus(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM),intent(in) :: compplus(rcpara%nmax,rcpara%pmax,rcpara%parmax),compminus(rcpara%nmax,rcpara%pmax,rcpara%parmax),alt

    !! splus is overwritten!
    real(kind=JPRM) :: splus(rcpara%nmax,rcpara%pmax,rcpara%parmax),sminus(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM),intent(in) :: rplus(rcpara%nmax,rcpara%pmax,1),rminus(rcpara%nmax,rcpara%pmax,1)
    real(kind=JPRM),intent(in) :: radplus(rcpara%nmax,rcpara%pmax,1),radminus(rcpara%nmax,rcpara%pmax,1)
    real(kind=JPRM),intent(in) :: bgrplus(rcpara%nmax,rcpara%pmax,1),bgrminus(rcpara%nmax,rcpara%pmax,1)
    real(kind=JPRM),intent(in) :: plusmean(rcpara%nmax,rcpara%parmax),minusmean(rcpara%nmax,rcpara%parmax)
    real(kind=JPRM) :: breakprofile(rcpara%pmax,rcpara%parmax),compprofile(rcpara%pmax,rcpara%parmax),f_val(rcpara%parmax),scale(rcpara%pmax), compbreakprofile(rcpara%pmax,rcpara%parmax)
    real(kind=JPRM) :: bgrbreakprofile(rcpara%pmax,1),brms(rcpara%parmax),bgbrms,rad(rcpara%pmax)
    real(kind=JPRM) :: rbhilf(rcpara%pmax,rcpara%parmax),rbuhilf(rcpara%pmax,rcpara%parmax),rbradhilf(rcpara%pmax,1)
    real(kind=JPRM) :: redfak(rcpara%pmax),dailycrut2(rcpara%nmax,rcpara%parmax)
    integer(kind=JPRM),intent(in) :: cardsmeta(rcpara%nmax,1,rcpara%nmeta),era40meta(rcpara%nmax)
    integer(kind=JPRM) :: iname,sigcount(rcpara%parmax),insigcount(rcpara%parmax),cibm1
    logical :: lsig
    character*6,intent(in) :: cname



    bcount=count(chosenbreaks .gt. 0)
    cib=chosenbreaks(ib)
    if(ib .gt. 1) then
       cibm1=chosenbreaks(ib-1)
    else
       cibm1=chosenbreaks(ib)+rcpara%mean_maxlen
    endif
    correctedbreaks=0
    bgrbreakprofile=rcpara%miss_val

!!$ call omp_set_lock(omp_lp)
    read(cname,'(I5)') iname
!!$ call omp_unset_lock(omp_lp)
    ic=0
    !!do ib=1,bcount

    !!attributionfactors(ib)=1.0

    rbhilf=0.
    rbuhilf=0.
    do ipar=1,rcpara%parmax
       do ip=1,rcpara%pmax
          if(sminus(cib,ip,ipar) .ne. rcpara%miss_val .and. splus(cib,ip,ipar) .ne. rcpara%miss_val .and. rminus(cib,ip,1) .ne. rcpara%miss_val .and. rplus(cib,ip,1) .ne. rcpara%miss_val .and. rminus(cib,ip,1)*rplus(cib,ip,1) .ne. 0)  then
             !        if(.false. .and. cib .lt. 9000 .and. ip .lt. 6 .and. (iname .gt. 47400 .and. iname .lt. 48000 .or. iname .gt. 20000 .and. iname .lt. 40000)) then
             !          splus(cib,ip,ipar)=0.8*splus(cib,ip,ipar)
             !          sminus(cib,ip,ipar)=0.8*sminus(cib,ip,ipar)
             !          rbhilf(ip,ipar)=-(splus(cib,ip,ipar)-sminus(cib,ip,ipar))
             !        else
             rbhilf(ip,ipar)=-(splus(cib,ip,ipar)-sminus(cib,ip,ipar))
             !        endif
             rbuhilf(ip,ipar)=sqrt((rplus(cib,ip,1)*rplus(cib,ip,1)+rminus(cib,ip,1)*rminus(cib,ip,1))/2.)

             if(ipar .eq. 1 .and. sminus(cib,ip,2) .ne. rcpara%miss_val .and. splus(cib,ip,2) .ne. rcpara%miss_val) then 
                rbradhilf(ip,1)=(splus(cib,ip,2)-sminus(cib,ip,2))-(splus(cib,ip,ipar)-sminus(cib,ip,ipar))
             endif
          endif
       enddo
    enddo

    call calc_profile(rcpara,iname,cib,splus,sminus,rplus,rminus,breakprofile,f_val,rcpara%smooth_method,rcpara%parmax) !in this file line 9
    where(breakprofile .eq. rcpara%miss_val) breakprofile=0.
    !bgrbreakprofile ist hier sbgrbreakprofile
    call calc_profile(rcpara,iname,cib,bgrplus,bgrminus,rplus,rminus,bgrbreakprofile,f_val,rcpara%smooth_method,1) !in this file line 9
    where(bgrbreakprofile .eq. rcpara%miss_val) bgrbreakprofile=0.

    call calc_profile(rcpara,iname,cib,compplus,compminus,rplus,rminus,compprofile,f_val,rcpara%smooth_method,rcpara%parmax) !in this file line 9
    where(compprofile .eq. rcpara%miss_val) compprofile=0.

    !  if (count(compprofile .ne. 0.) .ge. 8) then

    redfak=(/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.8,0.2,0.,0.,0./)
    compprofile(15:16,:)=0
    do ipar=1,rcpara%parmax
       if(compplus(cib,15,ipar).ne. rcpara%miss_val .and. compminus(cib,15,ipar).ne. rcpara%miss_val) then
          compprofile(15,ipar)=-(compplus(cib,15,ipar)-compminus(cib,15,ipar))
       endif
       if(compplus(cib,16,ipar).ne. rcpara%miss_val .and. compminus(cib,16,ipar).ne. rcpara%miss_val) then
          compprofile(16,ipar)=-(compplus(cib,16,ipar)-compminus(cib,16,ipar))
       endif
       if( abs(compprofile(15,ipar)) .lt. 0.5 .or. sqrt((rplus(cib,15,ipar)*rplus(cib,15,ipar)+rminus(cib,15,ipar)*rminus(cib,15,ipar))/2.) .gt. abs(compprofile(15,ipar))  &
            .or. 2.*abs(compprofile(16,ipar)) .gt. abs(compprofile(15,ipar)) .or. iter .eq. 1 .and. rcpara%maxiter .eq. 2) then
          do ip=1,rcpara%pmax
             compprofile(ip,ipar)=compprofile(ip,ipar)*redfak(ip)
          enddo
       else
!!$ call omp_set_lock(omp_lp)
          write(*,'(A5,I9,A18,3F8.3)') cname,todate(cib,rcpara),'Full surface break',compprofile(15,ipar),sqrt((rplus(cib,15,ipar)*rplus(cib,15,ipar)+rminus(cib,15,ipar)*rminus(cib,15,ipar))/2.),compprofile(16,ipar)
!!$ call omp_unset_lock(omp_lp)
          do ip=1,rcpara%pmax
             if(compprofile(ip,ipar) .gt. rcpara%miss_val) then
                compprofile(ip,ipar)=compprofile(ip,ipar)*redfak(ip)+compprofile(15,ipar)*(1.-redfak(ip))
             endif
          enddo

       endif
    enddo

    if(any(compprofile .ne. 0.) .or. count(f_val .gt. 1.) .gt. 0) then

       call test_significance(rcpara,iname,cib,breakprofile,bgrbreakprofile,compprofile,rbhilf,rbuhilf,rbradhilf,radplus,radminus,lsig,sigcount,insigcount) !in this file line 1623

       if(lsig) then 


          do ipar=1,rcpara%parmax
09           l=1
             do while(lasts(l,ipar) .gt. cib)
                l=l+1
             enddo
             if(l .eq. gcount(ipar)+1) then
                ibegin=1
             else
                ibegin=lasts(l,ipar)+rcpara%snht_maxlen/2.
             endif
             do ip=1,rcpara%pmax

                if(compprofile(ip,ipar) .ne. 0.) then

                   if(compminus(cib,ip,ipar) .ne. rcpara%miss_val .and. compplus(cib,ip,ipar) .ne. rcpara%miss_val)  then
                      rasocorrs(ibegin:cib,ip,ipar)=rasocorrs(ibegin:cib,ip,ipar)-compprofile(ip,ipar)
                      rasobreaks(ibegin:cib,ip,ipar)=rasobreaks(ibegin:cib,ip,ipar)+rbhilf(ip,ipar)
                      rasobreakuncertainties(ibegin:cib,ip,ipar)=rasobreakuncertainties(ibegin:cib,ip,ipar)+rbuhilf(ip,ipar)
                   endif

                   where(tfgm(ibegin:cib,ip,ipar) .ne. rcpara%miss_val)
                      tfgm(ibegin:cib,ip,ipar)=tfgm(ibegin:cib,ip,ipar)-compprofile(ip,ipar)
                   endwhere
                   where(tgps(ibegin:cib,ip,ipar) .ne. rcpara%miss_val)
                      tgps(ibegin:cib,ip,ipar)=tgps(ibegin:cib,ip,ipar)-compprofile(ip,ipar)
                   endwhere
                   where(stfgm(ibegin:cib,ip,ipar) .ne. rcpara%miss_val)
                      stfgm(ibegin:cib,ip,ipar)=stfgm(ibegin:cib,ip,ipar)-compprofile(ip,ipar)
                   endwhere
                endif
             enddo
             where(dailycrut2(ibegin:cib,ipar) .ne. rcpara%miss_val)
                dailycrut2(ibegin:cib,ipar)=dailycrut2(ibegin:cib,ipar)-compprofile(14,ipar)
             endwhere
          enddo

!!$ call omp_set_lock(omp_lp)
          write(protunit,'(A6,I6,A10,4I3)') 'break ',cib,' accepted',sigcount,insigcount

          do ip=1,rcpara%pmax
             write(protunit,'(I3,2F8.2)') ip,compprofile(ip,:)
          enddo
!!$ call omp_unset_lock(omp_lp)
          ic=ic+1
          correctedbreaks(ic)=cib

       else
!!$ call omp_set_lock(omp_lp)
          write(protunit,'(A6,I6,A10,4I3)') 'break ',cib,' rejected',sigcount,insigcount
          do ip=1,rcpara%pmax
             write(protunit,'(I3,2F8.2)') ip,compprofile(ip,:)
          enddo
!!$ call omp_unset_lock(omp_lp)
       endif
    else
!!$ call omp_set_lock(omp_lp)
       write(protunit,'(A6,I6,A10,3F7.3)') 'break ',cib,' rejected',f_val
       do ip=1,rcpara%pmax
          write(protunit,'(I3,2F8.2)') ip,compprofile(ip,:)
       enddo
!!$ call omp_unset_lock(omp_lp)
    endif
    !  endif


    return

  end subroutine adjust_series_comp

end module rfcor
