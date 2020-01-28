module rcorrect
contains

  subroutine raso_correct(rcpara,wmonrs,wmolats,wmolons,wmonames,wmostats,indexmeta,hilfmeta,ilmeta,meta_s,era40meta,bgcorrs,bgcorrs_e,bgcorrs_w,&
       densities,ominuse40,ominuse40_an,crut2,ini_correct,tfgmcr,tmcr,icache,needs_composite,statnr,gstat,istat,iter,skin,bad_intervals) !

    use rfmod
    use rfcor
    use rfcorio 
    use correct_breaks2
    use txtnc
    use tosat

    implicit none

    type(rasocor_namelist),intent(in) :: rcpara
    type(cacherecord) :: tfgmcr(rcpara%cachemax),tmcr(rcpara%cachemax)
    type(ecskin) :: skin
    type(metadata) :: meta_s

    integer         :: maxlen,i,ic,j,k,select_mode,iunit,ip,ipar,l,switchindex,lorig,ipresatend,ii,jj
    integer tdiff,tbi,tbindex(rcpara%mmax)
    integer,intent(in) :: wmostats,iter,istat,gstat,statnr,ilmeta
    real(kind=JPRM) :: locsig,break_fak,cachehilf(rcpara%nmax)

    real(kind=JPRM) :: diff(rcpara%nmax),plus(rcpara%nmax),minus(rcpara%nmax),prms(rcpara%nmax),mrms(rcpara%nmax),tsa(rcpara%nmax)
    integer         :: pcount(rcpara%nmax),mcount(rcpara%nmax)
    integer, intent(in) :: indexmeta(3000000),hilfmeta(3000000)
    integer         :: index(rcpara%statmax),omp_get_thread_num,imax,icache(rcpara%statmax+1)

    real(kind=JPRM) :: tfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tfg(rcpara%nmax,rcpara%pmax,rcpara%parmax),tanm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM) :: tm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tbcm(rcpara%nmax,rcpara%pmax,rcpara%parmax),fflags(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM) :: e20cm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM) :: stfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax),stanm(rcpara%nmax,rcpara%pmax,rcpara%parmax),stm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM) :: stbcm(rcpara%nmax,rcpara%pmax,rcpara%parmax),eracorrs(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM) :: rasocorrs(rcpara%nmax,rcpara%pmax,rcpara%parmax),remcorrs(rcpara%nmax,rcpara%pmax,rcpara%parmax),mrasocorrs(rcpara%mmax,rcpara%pmax,rcpara%parmax),trasocorrs(rcpara%nmax,rcpara%pmax,rcpara%parmax),delemask(rcpara%nmax,rcpara%pmax,rcpara%parmax),tgps(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    !!real(kind=JPRM) :: bgcorrs(rcpara%statmax,rcpara%pmax,rcpara%parmax,rcpara%brmax)
    real(kind=JPRM),intent(in)  :: bgcorrs(rcpara%nmax,rcpara%pmax,rcpara%parmax),bgcorrs_e(rcpara%nmax,rcpara%pmax,rcpara%parmax),bgcorrs_w(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    !real(kind=JPRM)  :: ancorrs(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    !!integer(kind=JPRM) :: bgcorrdates(2,rcpara%brmax)
    real(kind=JPRM) :: itfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax),itanm(rcpara%nmax,rcpara%pmax,rcpara%parmax),itm(rcpara%nmax,rcpara%pmax,rcpara%parmax),itfg12m(rcpara%nmax,rcpara%pmax,rcpara%parmax),newbc(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM)  :: solarangles(rcpara%nmax,rcpara%parmax)
    real(kind=JPRM),intent(in)  :: densities(rcpara%nmax,wmostats)
    integer(kind=JPRM)  :: tnumbers(rcpara%nmax,rcpara%pmax,rcpara%parmax),cardsmeta(rcpara%nmax,1,rcpara%nmeta),stgroups(rcpara%nmem,rcpara%ngroup)
    integer(kind=JPRM),intent(in)  :: era40meta(rcpara%nmax)
    integer,allocatable :: hours(:,:)
    real :: dists(rcpara%statmax)

    real(kind=JPRM) :: apriori_probs(rcpara%nmax),apriori_probs_rad(rcpara%nmax)
    real(kind=JPRM) :: break_thresh(rcpara%probmax),break_thresh_prob(rcpara%probmax)
    real(kind=JPRM) :: mweights(rcpara%pmax),tsaweights(rcpara%pmax),vtprcorr(rcpara%pmax)
    real(kind=JPRM) :: density(wmostats)

    real(kind=JPRM),intent(in) :: ominuse40(rcpara%ni,rcpara%nj,rcpara%pmax,rcpara%parmax,12),ominuse40_an(rcpara%ni,rcpara%nj,rcpara%pmax,rcpara%parmax,12)
    real(kind=JPRM) :: adjust(rcpara%pmax,rcpara%parmax,12),adjust_an(rcpara%pmax,rcpara%parmax,12)
    real(kind=JPRM) :: crut2(72,36,rcpara%mmax),dailycrut2(rcpara%nmax,rcpara%parmax)

    integer,intent(in) :: wmonrs(rcpara%statmax)
    logical :: ini_correct(rcpara%parmax,rcpara%statmax),logcache
    logical :: composite_exists(rcpara%parmax)
    real(kind=JPRM),intent(in) :: wmolats(rcpara%statmax),wmolons(rcpara%statmax)
    REAL (kind=JPRM):: lat, lon , stype(rcpara%nmax),alt
    character*50 :: wmonames(rcpara%statmax)
    !!character*64 :: stgroupnames(:)
    character binprefix*20,startdate*8,cdatum*8,cstatnr*6, cstatnr2*5, filename*100, cgstat*5
    integer err,err1,err2,err3,cerr
    integer prob_method,protunit,blunit,goodsondes(28),id
    logical ex,ex1,ex2,ex3
    integer rtype(rcpara%nmax),mr,cmax,icistat,iosmeta,tindex(rcpara%nmax)
    integer         :: lasts(rcpara%brmax,rcpara%parmax),midx(rcpara%brmax,rcpara%parmax),gcount(rcpara%parmax),lastsave(rcpara%brmax,rcpara%parmax),needs_composite(:,:,:)
    logical :: iczwisch,surfonly

    integer iindex,iindexsave,imon,imod,tgm,iy
    real(kind=JPRM) :: tmmon(rcpara%mmax,rcpara%pmax,rcpara%parmax),tsum

    real(kind=JPRM),allocatable,dimension(:,:,:,:) :: arr
    INTEGER,allocatable,dimension(:,:) :: hdatum,hhours

    integer:: maxatts,vals,remerr
    integer,allocatable:: datum(:,:)!,hours(:,:)
    integer:: bad_intervals(:,:)

!!$ call omp_set_lock(omp_lp)
    write(*,*) 'statnr ',statnr
!!$ call omp_unset_lock(omp_lp)

    !!allocate(cardsmeta(nmax,1,nmeta)) ! wmostats waeren alle stationen (zu viele)

    cardsmeta=0
    remerr=1
    err1=1

    !!if(gstat .eq. 0) return
    iunit=20
    !$ ic=omp_get_thread_num()
    !$ iunit=iunit+ic ! this should avoid clashes in parallel runs

!!$ call omp_set_lock(omp_lp)
    write(cstatnr,'(I6.6)') statnr

!!$ call omp_unset_lock(omp_lp)
    if(rcpara%innov .eq. 'MO') THEN
       filename=trim(rcpara%prefix)//cstatnr//'/feedbackmerged'//cstatnr
       inquire(file=filename,exist=ex3)
       !!  filename=trim(rcpara%prefix)//'feedbackglobbinsaveera'//cstatnr
       !!  inquire(file=filename,exist=ex1)
       !!  filename=trim(rcpara%prefix)//'feedbackglobbinsaveoper'//cstatnr
       !!  inquire(file=filename,exist=ex2)
       ex=ex3 !!.or. ex1 .or. ex2
    else if(rcpara%innov .eq. 'NN' .or. rcpara%innov .eq. 'NE' .or. rcpara%innov .eq. 'NR') then
       !  filename=trim(rcpara%prefix)//cstatnr//'/'//cstatnr//'_t.nc'
       filename=trim(rcpara%prefix)//cstatnr//'/feedbackmerged'//cstatnr//'.nc'
       inquire(file=filename,exist=ex3)
       filename=trim(rcpara%prefix)//cstatnr//'/'//cstatnr//'_t.nc'
       inquire(file=filename,exist=ex1)
       !!  filename=trim(rcpara%prefix)//'feedbackglobbinsaveera'//cstatnr
       !!  inquire(file=filename,exist=ex1)
       !!  filename=trim(rcpara%prefix)//'feedbackglobbinsaveoper'//cstatnr
       !!  inquire(file=filename,exist=ex2)
       ex=ex3 .or. ex1! .or. ex2
    else 
!!$ call omp_set_lock(omp_lp)
       write(*,*) 'innov must be either U or E or M or EO or MO or NN'
!!$ call omp_unset_lock(omp_lp)
       call exit(1)
    endif

    if(.not. ex) then
!!$ call omp_set_lock(omp_lp)
       print*,trim(filename),' not found'
!!$ call omp_unset_lock(omp_lp)
       return
    endif

    iosmeta=1
    if(iter .gt. 1) then
!!$ call omp_set_lock(omp_lp2)
    endif
    if(rcpara%innov .eq. 'MO') THEN
       if(ex3) then
          !!      read_sonde_oper calculates adjust .. do not remove call
          filename='dummy'
          solarangles(1,1)=1.0
          call read_sonde_oper(iunit,filename,rcpara,wmolats(istat),wmolons(istat),tm,tanm,tfgm,tbcm,solarangles,ominuse40,adjust,rtype,err2) !in file rfcorio.f90 line 611

          !! ERA/ECMWF data must be read anyway to get GTS sonde types
          filename=trim(rcpara%prefix)//cstatnr//'/metadata'//cstatnr
          iosmeta=1
!          open(iunit,file=filename,form='unformatted',status='old',action='read',iostat=iosmeta)
          if(iosmeta .eq. 0) then
             read(iunit) rtype
             read(iunit) apriori_probs
             read(iunit) apriori_probs_rad
             close(iunit)
!!$ call omp_set_lock(omp_lp)
             write(*,*) statnr,' read metadata file'
!!$ call omp_unset_lock(omp_lp)
          else   
             filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbinsaveera'//cstatnr
             solarangles(1,1)=1.
             call read_sonde(iunit,filename,rcpara,tm,tanm,tfgm,tbcm,solarangles,rtype,err1) !in file rfcorio.f90 line 373
             filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbinsaveoper'//cstatnr
             solarangles(1,1)=1.0
             call read_sonde_oper(iunit,filename,rcpara,wmolats(istat),wmolons(istat),tm,tanm,tfgm,tbcm,solarangles,ominuse40,adjust,rtype,err2) !in file rfcorio.f90 line 611
          endif
          filename=trim(rcpara%prefix)//cstatnr//'/feedbackmerged'//cstatnr
          call read_igrasonde_daily_new(iunit,filename,rcpara%nmax,rcpara%pmax,rcpara%parmax,tm,tanm,tfgm,itfg12m,rcpara%miss_val,rcpara%snht_maxlen/2,err3) !in file rfmod.f90 line 2044
          call rtypeinfo(rtype,rcpara,mr,istat,imax) !in file rfcorio.f90 line 847


       endif
       err=err3 !!*err1*err2
    else if(rcpara%innov .eq. 'NN' .or. rcpara%innov .eq. 'NE' .or. rcpara%innov .eq. 'NR') THEN
       if(ex) then
          !    filename=trim(rcpara%prefix)//cstatnr//'/feedbackmerged'//cstatnr
          !    call read_igrasonde_daily_new(iunit,filename,rcpara%nmax,rcpara%pmax,rcpara%parmax,tm,tanm,tfgm,itfg12m,rcpara%miss_val,rcpara%snht_maxlen/2,err3)
          !    if(err3 .eq. 2) return
          !! ERA/ECMWF data must be read anyway to get GTS sonde types
          iosmeta=1
          filename=trim(rcpara%prefix)//cstatnr//'/metadata'//cstatnr
          !    open(iunit,file=filename,form='unformatted',status='old',action='read',iostat=iosmeta)
          if(iosmeta .eq. 0) then
             read(iunit) rtype
             read(iunit) apriori_probs
             read(iunit) apriori_probs_rad
             close(iunit)
!!$ call omp_set_lock(omp_lp)
             write(*,*) statnr,' read metadata file'
!!$ call omp_unset_lock(omp_lp)
          else   
             filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbinsaveera'//cstatnr
             solarangles(1,1)=1.
!                   call read_sonde(iunit,filename,rcpara,stm,stanm,stfgm,tbcm,solarangles,rtype,err1) !in file rfcorio.f90 line 373
             filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbinsaveoper'//cstatnr
             solarangles(1,1)=1.0
!                   call read_sonde_oper(iunit,filename,rcpara,wmolats(istat),wmolons(istat),stm,stanm,stfgm,tbcm,solarangles,ominuse40,adjust,rtype,err2) !in file rfcorio.f90 line 611
          endif
          !    call rtypeinfo(rtype,rcpara,mr,istat,imax) !in file rfcorio.f90 line 847

          filename=trim(rcpara%prefix)//cstatnr//'/feedbackmerged'//cstatnr//'.nc'
          if(rcpara%innov .eq. 'NN') then
            CALL read_odb_nc(filename,rcpara,istat,err3,tm,tfgm,tanm=tanm,stname=wmonames(istat),bad_intervals=bad_intervals,tgps=tgps,alt=alt) !subroutine in read_txt_write_nc.f90
!!$            where (isnan(tm) .or. abs(tfgm)>20. .or. tm>330. .or. tm<170.)
!!$               tm=rcpara%miss_val
!!$               tfgm=rcpara%miss_val
!!$               tanm=rcpara%miss_val
!!$            endwhere
            do ipar=1,rcpara%parmax
               do ip=1,rcpara%pmax
                  do i=1,rcpara%nmax
                     if (isnan(tm(i,ip,ipar)) .or. abs(tfgm(i,ip,ipar))>20.) then
                         tm(i,ip,ipar)=rcpara%miss_val
                         tfgm(i,ip,ipar)=rcpara%miss_val
                         tanm(i,ip,ipar)=rcpara%miss_val
                     else
!                         tfgm(i,ip,ipar)=-tfgm(i,ip,ipar)
!                         tanm(i,ip,ipar)=-tanm(i,ip,ipar)
                     endif
                  enddo
                enddo
            enddo

          else
            CALL read_odb_nc(filename,rcpara,istat,err3,tm,tfgm,e20c0=tanm,stype=stype,stname=wmonames(istat),tgps=tgps) !subroutine in read_txt_write_nc.f90
             tfgm=tanm
          endif

!          do i=1,rcpara%nmax
             !     if(rtype(i) .ne. -999. .or. stype(i) .ne. -999.) print*,i,rtype(i),stype(i)
!          enddo

          rtype=floor(stype)
          call rtypeinfo(rtype,rcpara,mr,istat,imax) !in file rfcorio.f90 line 847

          maxatts=4
          vals=rcpara%nmax

          !allocate(datum(rcpara%nmax,1),hours(rcpara%nmax,rcpara%parmax))
          !call txttonc(rcpara,istat,filename,(/tm,tfgm,tbcm,fflags,tanm/),(/'temperatures','fg_dep      ','bias        ','flags       ','an_dep      '/),vals,5,datum, hours, 0._JPRM,0._JPRM,0)     
          !    CALL read_odb_nc(filename,rcpara,istat,err3,tm,tfgm,tbcm,fflags,tanm) !subroutine in read_txt_write_nc.f90

          call eiminuse40(rcpara,wmolats(istat),wmolons(istat),ominuse40,adjust)
          call eiminuse40(rcpara,wmolats(istat),wmolons(istat),ominuse40_an,adjust_an)

          switchindex=toindex(rcpara%switchdate,rcpara)
          do ipar=1,rcpara%parmax
             do ip=1,rcpara%pmax
                do i=1,rcpara%nmax
                   if(tfgm(i,ip,ipar) .ne. rcpara%miss_val .and. tm(i,ip,ipar) .ne. rcpara%miss_val) then
                      tfg(i,ip,ipar)=tm(i,ip,ipar)-tfgm(i,ip,ipar)
                      if(tanm(i,ip,ipar) .ne. rcpara%miss_val) then
                         tanm(i,ip,ipar)=tm(i,ip,ipar)-tanm(i,ip,ipar)
                      else
                         tanm(i,ip,ipar)=rcpara%miss_val
                      endif
                      if(tgps(i,ip,ipar) .ne. rcpara%miss_val) then
                         tgps(i,ip,ipar)=tfgm(i,ip,ipar)-tgps(i,ip,ipar)
                      else
                         tgps(i,ip,ipar)=rcpara%miss_val
                      endif
                   else
                      tanm(i,ip,ipar)=rcpara%miss_val
                      tfg(i,ip,ipar)=rcpara%miss_val
                      tgps(i,ip,ipar)=rcpara%miss_val
                   endif
                enddo
             enddo
          enddo
       endif
       err=err1*err3 !!*err1*err2
    endif ! rcpara%innov

    if(err .ne. 0) then
!!$ call omp_set_lock(omp_lp)
       write(*,*) 'found but could not read ',filename
!!$ call omp_unset_lock(omp_lp)
       !  call exit(1)
       return
    endif
    if(iter .gt. 1) then
!!$ call omp_unset_lock(omp_lp2)
    endif

    if(rcpara%nmax .eq. 23010) then
      call read_cards_meta(0,rcpara,wmostats,statnr,wmonrs,cardsmeta,indexmeta,hilfmeta,ilmeta,err) !in file rfcorio.f90 line 2081

      call create_meta(cardsmeta,rtype,statnr,rcpara,iter) !in file rfcorio.f90 line 891
    endif
    
    call check_ini(meta_s,istat,wmonrs,rcpara,ini_correct)

    protunit=100
    !$ protunit=protunit+omp_get_thread_num()

!!$ call omp_set_lock(omp_lp)
    open(protunit,file=trim(rcpara%prefix)//cstatnr//'/found_breaks'//cstatnr,form='formatted')
!!$ call omp_unset_lock(omp_lp)

!    needs_composite(:,:,istat)=.false.
    call detect_gaps(rcpara,statnr,tfgm,midx,lasts,gcount) !in file rfcor.f90 line 113
!    lastsave=lasts
    if (.not. any(gcount>20)) then
       needs_composite(:,:,istat)=lasts(1:20,:)
    else
       write(*,*) statnr,' too many gaps',gcount
    endif
!!$    do ipar=1,rcpara%parmax
!!$       if(gcount(ipar) .gt. size(needs_composite,1)) then 
!!$          write(*,*) statnr,' too many gaps',ipar,gcount(ipar)
!!$          stop
!!$       endif
!!$       do i=1,gcount(ipar)
!!$          needs_composite(i,ipar,istat)=(lasts(i,ipar) .lt. rcpara%old .or. ini_correct(ipar,istat))
!!$          !! benutze first guess nur wenn beide sehr rezent sind
!!$          !!    if(lasts(i,ipar) .ge. rcpara%old .and. lasts(1,ipar-3) .lt. rcpara%old) needs_composite(i,ipar,istat)=.true.
!!$       enddo
!!$       if(lasts(1,ipar) .gt. 0 .and.(lasts(1,ipar) .lt. rcpara%old)) ini_correct(ipar,istat)=.true.
!!$    enddo


!!$ call omp_set_lock(omp_lp)
    write(*,'(6I6,A17,2L2,A6,I6)') wmonrs(istat),imax,lasts(1,1),lasts(1,2),count(tfgm(:,10,1) .ne. rcpara%miss_val),count(tfgm(:,10,2) .ne. rcpara%miss_val),' initial_correct:',ini_correct(:,istat),' imax: ',imax, rcpara%old
!!$ call omp_unset_lock(omp_lp)

    lorig=0 ! original number of days with data
    do i=1,rcpara%nmax
       if((tfgm(i,12,1) .ne. rcpara%miss_val .or. tfgm(i,12,2) .ne. rcpara%miss_val) .and.  (tm(i,12,1) .ne. rcpara%miss_val .or.  tm(i,12,2) .ne. rcpara%miss_val)) then
          lorig=lorig+1 
          tindex(lorig)=i           
       endif
    enddo

    index=0
    composite_exists=.false.
    stm=rcpara%miss_val
    stfgm=rcpara%miss_val
    if(any(gcount .gt. 0)) then !257
       !!  if(iter .lt. 3 .and. .not. any(needs_composite(:,:,istat))) then 
       if(iter .lt. 3 .and. any(.not. ini_correct(:,istat)) ) then !259
          cerr=1
          filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbincomp'//cstatnr
          !!      call read_sonde_daily(iunit,filename,rcpara,cerr,stm,stfgm)
          if(cerr .ne. 0) then
             cmax=wmostats
             if(iter .eq. 2) then !265
                write(*,*) 'gcount if',statnr,gcount
                call make_composite(rcpara,statnr,cmax,wmonrs,wmolons,wmolats,wmostats,dists,index,ominuse40,adjust,stm,stfgm,tnumbers,&
                     tfgmcr,tmcr,icache,meta_s,lasts,gcount,needs_composite,ini_correct,composite_exists,bad_intervals) !in file rfcor.f90 line 426
                if(any(composite_exists)) then
                   filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbincomp'//cstatnr
                   where(stm(:,:,2) .ne. rcpara%miss_val .and. stm(:,:,1) .ne. rcpara%miss_val)
                      tfgm(:,:,1)=stm(:,:,2)-stm(:,:,1)
                   elsewhere
                      tfgm(:,:,1)=0.
                   endwhere
                   !call write_sonde_daily(iunit,filename,rcpara,err,stm,stfgm)
                endif
                close(protunit)
                return !?????
             else !265
                   do ipar=1,rcpara%parmax
                      l=0
                      do  while(lasts(l,ipar)>1)
                         if(lasts(l,ipar) .lt. rcpara%old) then 
                            tfgm(l:lasts(1,ipar),:,ipar)=rcpara%miss_val
                            tm(l:lasts(1,ipar),:,ipar)=rcpara%miss_val
                         endif
                         l=l+1
                      enddo
                   enddo
             endif !265
          endif !263
          !!    else
          !!      stm=rcpara%miss_val
          !!      stfgm=rcpara%miss_val
       else !iter .lt. 3 .and. any (.not. ini_correct) (259)
          if(iter .lt. rcpara%maxiter) then
             close(protunit)
             return
          endif
          cerr=1
          !!    if(.not. ini_correct(istat))  call read_sonde_daily(iunit,filename,rcpara,cerr,stm,stfgm)
          if(cerr .ne. 0 .and. rcpara%maxiter .gt. 1) then !302
             cmax=wmostats
             composite_exists=.false.
!!$ call omp_set_lock(omp_lp)
             write(*,*) 'gcount else',statnr,gcount
!!$ call omp_unset_lock(omp_lp)
             call make_composite(rcpara,statnr,cmax,wmonrs,wmolons,wmolats,wmostats,dists,index,ominuse40,adjust,stm,stfgm,tnumbers,&
                  tfgmcr,tmcr,icache,meta_s,lasts,gcount,needs_composite,ini_correct,composite_exists,bad_intervals)!in file rfcor.f90 line 426
          else
             composite_exists=.true.
          endif !302
       endif !259
    endif !257

    if(any(gcount .gt. 0) .and. (any(.not. ini_correct(:,istat)) .and. iter .lt. 3 .or. any(composite_exists)) .or. iter .eq. rcpara%maxiter) then !315



       call adjust_bg_ei(cstatnr,istat,rcpara,tm,bgcorrs,densities,adjust_an,eracorrs,wmolons,wmolats,wmostats,iunit) 

       filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbguncorrmon'//cstatnr//'.nc'
       if(ex .and. rcpara%extended_output .eq. 'Y' .and. count(tm(:,8,:) .gt. rcpara%miss_val) .gt. 100) then
          call write_sonde_monthly_nc(filename,rcpara,tfg,istat,err,wmolons(istat),wmolats(istat),tanm,stname=wmonames(istat)) !in file rfcorio.f90 line 1823
       endif
       where(tanm .ne. rcpara%miss_val)
          tanm=tanm+eracorrs
       endwhere

       call adjust_bg_ei(cstatnr,istat,rcpara,tm,bgcorrs,densities,adjust,eracorrs,wmolons,wmolats,wmostats,iunit) 

       if(any(tfg .lt. -700. .and. tfg .ne. rcpara%miss_val) ) stop 'STOP: spurious tfg values'
       where(tfg .ne. rcpara%miss_val)
          tfg=tfg+eracorrs
       endwhere

       filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbgmon'//cstatnr//'.nc'
       if(ex .and. rcpara%extended_output .eq. 'Y' .and. count(tm(:,8,:) .gt. rcpara%miss_val) .gt. 100) then
          call write_sonde_monthly_nc(filename,rcpara,tfg,istat,err,wmolons(istat),wmolats(istat),tanm,stname=wmonames(istat)) !in file rfcorio.f90 line 1823
       endif

!       filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobanmon'//cstatnr//'.nc'
       if(iunit .gt. 0 .and. rcpara%extended_output .eq. 'Y') then

!          call write_sonde_monthly_nc(filename,rcpara,tanm,istat,err,wmolons(istat),wmolats(istat),tanm) !in file   endif


!!$ call omp_set_lock(omp_lp)
!          write(*,*) wmonrs(istat),iter,any(gcount .gt. 0),.not. any(needs_composite(1,:,istat)) ,iter .eq. 1,any(composite_exists),' write correction and composite'
!!$ call omp_unset_lock(omp_lp)
iosmeta=1

          if(iosmeta .ne. 0) then
             call apriori_prob(rcpara,cstatnr,cardsmeta,era40meta,apriori_probs,apriori_probs_rad,iter) !in file rfcor.f90 line 1615


             call apriori_prob_schroeder(rcpara,cstatnr,meta_s.cardsmeta_s(:,istat),era40meta,apriori_probs,apriori_probs_rad,iter) !in file rfcor.f90 line 1615
             filename=trim(rcpara%prefix)//cstatnr//'/metadata'//cstatnr
             open(iunit,file=filename,form='unformatted',iostat=iosmeta)
             if(iosmeta .eq. 0) then
                write(iunit) rtype
                write(iunit) apriori_probs
                write(iunit) apriori_probs_rad
                close(iunit)
             endif
!!$ call omp_set_lock(omp_lp)
             write(*,*) statnr,' created and wrote metadata file'
!!$ call omp_unset_lock(omp_lp)
          endif
       endif


       rasocorrs=0.
       delemask=0.
       id=0
       do jj=1,2
        do ii=1,rcpara%statmax
         if (.not. ini_correct(jj,ii)) id=id+1
        enddo
       enddo
       write(*,*) 'Trusted 0012 stations:',id

       if(rcpara%innov .NE. 'NR' .or. remerr .eq. 1) THEN
          call correct_break(cstatnr,istat,wmolons(istat),wmolats(istat),rcpara,tm,tfgm,solarangles,stm,stfgm,tnumbers,tfg,tgps,rasocorrs,eracorrs,crut2,&
               dailycrut2,delemask,cardsmeta,era40meta,apriori_probs,apriori_probs_rad,iter,needs_composite(:,:,istat),ini_correct(:,istat),midx,lasts,gcount,protunit,iunit,alt) !in file correct_breaks2.f90 line 34
       ELSE
          call correct_break(cstatnr,istat,wmolons(istat),wmolats(istat),rcpara,tm,tfgm,solarangles,stm,stfgm,tnumbers,tfg,tgps,rasocorrs,eracorrs,crut2,&
               dailycrut2,delemask,cardsmeta,era40meta,apriori_probs,apriori_probs_rad,iter,needs_composite(:,:,istat),ini_correct(:,istat),midx,lasts,gcount,protunit,iunit,alt) !in file correct_breaks2.f90 line 34
       endif


          !$ call omp_set_lock(omp_lp(rcpara%statmax+2))
          logcache=icache(istat) .eq. 0
          if(logcache) then
             icache(rcpara%statmax+1)=icache(rcpara%statmax+1)+1
             icache(istat)=icache(rcpara%statmax+1)
          endif
             icistat=icache(istat)
             !$ call omp_unset_lock(omp_lp(rcpara%statmax+2))

             !$ call omp_set_lock(omp_lp(icistat))

           if(logcache) then
             l=0
             do i=1,rcpara%nmax
                if((tfgm(i,12,1) .ne. rcpara%miss_val .or. tfgm(i,12,2) .ne. rcpara%miss_val) .and.  (tm(i,12,1) .ne. rcpara%miss_val .or.  tm(i,12,2) .ne. rcpara%miss_val)) then
                   l=l+1            
!haim                   tfgmcr(icistat)%index(l)=i
                   cachehilf(l)=i
                endif
             enddo
             tfgmcr(icistat)%vals=l
             allocate(tfgmcr(icistat)%index(l))
             tfgmcr(icistat)%index=cachehilf(1:l)
             allocate(tfgmcr(icistat)%feld(l,rcpara%pmax,rcpara%parmax))
             allocate(tmcr(icistat)%feld(l,rcpara%pmax,rcpara%parmax))
           endif

             do ipar=1,rcpara%parmax
                do ip=1,rcpara%pmax
                   do l=1,tfgmcr(icistat)%vals
                      i=tfgmcr(icistat)%index(l)
                      if(tfgm(i,ip,ipar) .ne. rcpara%miss_val .and. tm(i,ip,ipar) .ne. rcpara%miss_val) then
                         tfgmcr(icistat)%feld(l,ip,ipar)=tfgm(i,ip,ipar)+rasocorrs(i,ip,ipar)+eracorrs(i,ip,ipar)
                         tmcr(icistat)%feld(l,ip,ipar)=tm(i,ip,ipar)-rasocorrs(i,ip,ipar)
!!$                         if(tmcr(icistat)%feld(l,ip,ipar).gt. 350. .or. tmcr(icistat)%feld(l,ip,ipar).lt. 150.) then
!!$write(*,*) cstatnr,' STOP: raso_correct: wrong temperature'
!!$stop
!!$endif
                      else
                         tfgmcr(icistat)%feld(l,ip,ipar)=rcpara%miss_val
                         tmcr(icistat)%feld(l,ip,ipar)=rcpara%miss_val
                      endif
                   enddo
                enddo
             enddo
             write(*,*) statnr,' written to cache with index',icache(istat)
             !$ call omp_unset_lock(omp_lp(icistat))


       if(iter .gt. 1 .and. any(composite_exists) .and. rcpara%extended_output .eq. 'Y') then
          filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbincompmon'//cstatnr//'.nc'
          where(tm .ne. rcpara%miss_val .and. stfgm .ne. rcpara%miss_val .and. tfgm .ne. rcpara%miss_val)
             stanm=tm+tfgm-stfgm
          elsewhere
             stanm=rcpara%miss_val
             stfgm=rcpara%miss_val
             stm=rcpara%miss_val
          endwhere
          call write_sonde_monthly_nc(filename,rcpara,stm,istat,err,lon, lat, stfgm,stanm,stname=wmonames(istat)) !in file read_txt_write_nc line 796

          where(tm .eq. rcpara%miss_val)
             stm=rcpara%miss_val
             stfgm=rcpara%miss_val
          endwhere
          filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbincomp'//cstatnr//'.nc'
          allocate(arr(tfgmcr(icistat)%vals,rcpara%pmax,rcpara%parmax,2),hours(tfgmcr(icistat)%vals,rcpara%parmax))
          do i=1,tfgmcr(icistat)%vals
             arr(i,:,:,1)=stm(tfgmcr(icistat)%index(i),:,:)
             arr(i,:,:,2)=stfgm(tfgmcr(icistat)%index(i),:,:)
          enddo
          CALL txttonc(rcpara,istat,filename,arr,(/'comp-obs    ','comp-bg     '/),tfgmcr(icistat)%vals, 2,tfgmcr(icistat)%index, hours, wmolats(istat),wmolons(istat),0,stname=wmonames(istat))
          deallocate(hours,arr)
       endif

       filename=trim(rcpara%prefix)//cstatnr//'/feedbackbgcorr'//cstatnr//'.nc'
       allocate(arr(lorig,rcpara%pmax,rcpara%parmax,1),hours(lorig,rcpara%parmax),stat=err)
       do i=1,lorig
          arr(i,:,:,1)=eracorrs(tindex(i),:,:)
       enddo
       CALL txttonc(rcpara,istat,filename,arr,(/'bg correction'/),lorig, 1,tindex(1:lorig), hours, wmolats(istat),wmolons(istat),0,stname=wmonames(istat))
       deallocate(hours,arr)

         write(*,'(a,i1,1x,I6,1x,A,1x,i6,5F13.4)') 'raso_correct ',iter,wmonrs(istat),'::',tfgmcr(icistat)%vals,sum(rasocorrs)

       filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbincorrsave'//cstatnr//'.nc'
       where(rasocorrs .eq. rcpara%miss_val) rasocorrs=0.

       trasocorrs=rasocorrs
       !    call write_sonde_corr_daily(iunit,filename,rcpara,err,rasocorrs,eracorrs-eracorrs)
       lat=wmolats(istat)
       lon=wmolons(istat) 
       call write_sonde_corr_daily_nc(filename,rcpara,istat,err,-rasocorrs,wmolons(istat),wmolats(istat),eracorrs,stname=wmonames(istat)) !in file read_txt_write_nc line 471

       ! Activate for initialization of "removal of signal" experiments
       if(remerr .ne. 0 .and. iter .ge. rcpara%maxiter .and. rcpara%innov .eq. 'NR') THEN
          call remove_signal(rcpara,cstatnr,tm,remcorrs,toindex(19790101,rcpara),toindex(20061231,rcpara),.false.,.true.)

          filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbinremsave'//cstatnr//'.nc'
          !      trasocorrs=rasocorrs
          where(remcorrs .eq. rcpara%miss_val) remcorrs=0.
          !      call write_sonde_corr_daily(iunit,filename,rcpara,err,rasocorrs,eracorrs-eracorrs)
          lat=wmolats(istat)
          lon=wmolons(istat) 
          call write_sonde_corr_daily_nc(filename,rcpara,istat,err,remcorrs,wmolons(istat),wmolats(istat),eracorrs,stname=wmonames(istat)) !in file read_txt_write_nc line 471

          !    filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbincorrsave'//cstatnr//'.nc'
          !    call write_sonde_corr_daily_nc(filename,rcpara,istat,err,trasocorrs,wmolons(istat),wmolats(istat),eracorrs) !in file read_txt_write_nc line 471
          where(trasocorrs .ne. rcpara%miss_val)
             rasocorrs=-trasocorrs
          endwhere
       endif

!       where(rasocorrs .ne. rcpara%miss_val)
!          rasocorrs=-rasocorrs
!       endwhere
       !!print*,'write ',filename
       filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbincorrmon'//cstatnr//'.nc'
       where(tfgm .ne. rcpara%miss_val .and. tm .ne. rcpara%miss_val .and.delemask .ne. rcpara%miss_val)
          tm=tm+rasocorrs
       elsewhere
          tm=rcpara%miss_val
          rasocorrs=rcpara%miss_val
          eracorrs=rcpara%miss_val
       endwhere
       where(dailycrut2 .ne. rcpara%miss_val .and. tm(:,11,:) .ne. rcpara%miss_val)
          tm(:,16,:)=dailycrut2
       elsewhere
          tm(:,16,:)=rcpara%miss_val
       endwhere
       ic=0
       do k=1,rcpara%parmax
          do j=1,rcpara%pmax
             do i=1,rcpara%nmax
                if(tm(i,j,k) .ne. rcpara%miss_val .and. tm(i,j,k) .lt. -700. ) then
                   ic=ic+1

!!$ call omp_set_lock(omp_lp)
                   write(*,*) statnr,' invalid values ',i,j,k,tm(i,j,k) 
!!$ call omp_unset_lock(omp_lp)
                endif
             enddo
          enddo
       enddo
       if(ic .gt. 0)  call exit(1)
       !call write_sonde_monthly(iunit,filename,rcpara,tm,err,rasocorrs,eracorrs) !in file rfcorio.f90 line 1823
       call write_sonde_monthly_nc(filename,rcpara,tm,istat,err,wmolons(istat),wmolats(istat),rasocorrs,eracorrs,stname=wmonames(istat))
       filename=filename(1:len(trim(filename))-9)//'bt'//cstatnr//'.nc'
       if(remerr .ne. 0 .and. iter .ge. rcpara%maxiter .and. rcpara%innov .eq. 'NR') THEN

          do k=1,rcpara%parmax
             do j=1,rcpara%pmax
                ic=0
                tsum=0.
                do i=1,rcpara%nmax
                   if(tm(i,j,k) .ne. rcpara%miss_val ) then
                      ic=ic+1
                      tsum=tsum+tm(i,j,k)
                   endif
                enddo
                do i=1,rcpara%nmax
                   if(tm(i,j,k) .ne. rcpara%miss_val)   tm(i,j,k)=tsum/ic+rasocorrs(i,j,k)
                enddo
             enddo
          enddo
       endif
       call write_sonde_monthly_bt_nc(filename,rcpara,tm,istat,err,wmolons(istat),wmolats(istat),4,skin,rasocorrs,tfg,tanm,stname=wmonames(istat))
       filename=filename(1:len(trim(filename))-9)//'2'//cstatnr//'.nc'
       call write_sonde_monthly_bt_nc(filename,rcpara,tm,istat,err,wmolons(istat),wmolats(istat),2,skin,rasocorrs,tfg,tanm,stname=wmonames(istat))
       !in file read_txt_write_nc.f90 line 799

    endif

    close(protunit)

    return
  end subroutine raso_correct

end module rcorrect
