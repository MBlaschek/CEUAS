module comp_c

contains

  subroutine comp_correct(rcpara,wmonrs,wmolats,wmolons,wmonames,wmostats,meta_s,era40meta,bgcorrs,bgcorrs_e,bgcorrs_w,densities,&
       ominuse40,crut2,ini_correct,tfgmcr,tmcr,tccr,lastscache,gcountcache,tbicache,icache,needs_composite,statnr,gstat,istat,ex,iter,alarm,skin,bad_intervals) !

    use rfmod
    use rfcor
    use rfcorio
    use correct_breaks2
    use rfcomp2
    use rfcomp_1
    use tosat

    implicit none

    type(rasocor_namelist),intent(in) :: rcpara
    type(cacherecord) :: tfgmcr(rcpara%cachemax),tmcr(rcpara%cachemax),tccr(rcpara%cachemax,3)
    type(ecskin)  :: skin
    type(metadata) :: meta_s

    integer         :: maxlen,i,select_mode,iunit,ip,ipar,l,switchindex,bi
    integer tdiff,tbi,tbindex(rcpara%mmax)
    integer,intent(in) :: wmostats,iter,istat,gstat,statnr
    real(kind=JPRM) :: locsig,break_fak
    REAL (kind=JPRM):: lat, lon 

    real(kind=JPRM) :: diff(rcpara%nmax),plus(rcpara%nmax),minus(rcpara%nmax),prms(rcpara%nmax),mrms(rcpara%nmax)
    integer(kind=JPRM) :: philf(rcpara%nmax),mhilf(rcpara%nmax)
    real(kind=JPRM) :: tsa(rcpara%nmax,rcpara%pmax,rcpara%parmax,2),tsa2(rcpara%nmax,rcpara%pmax,rcpara%parmax,2),tsa3(rcpara%nmax,rcpara%pmax,rcpara%parmax,2),null(rcpara%nmax)
    real(kind=JPRM) :: rasobreaks(rcpara%nmax,rcpara%pmax,rcpara%parmax), rasobreakuncertainties(rcpara%nmax,rcpara%pmax,rcpara%parmax)

    integer         :: pcount(rcpara%nmax),mcount(rcpara%nmax),critical_dates(10),istart
    integer         :: index(rcpara%statmax),omp_get_thread_num,imax,icache(rcpara%statmax+1)
    logical         :: ex(rcpara%statmax)

    integer         :: lastscache(rcpara%brmax,rcpara%parmax,rcpara%cachemax),gcountcache(rcpara%parmax,rcpara%cachemax)
    real(kind=JPRM) :: stfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax),stanm(rcpara%nmax,rcpara%pmax,rcpara%parmax),stm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM) :: stfgmsave(rcpara%nmax,rcpara%pmax,rcpara%parmax),tfgmsave(rcpara%nmax,rcpara%pmax,rcpara%parmax),tgpssave(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM) :: tfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tanm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tbcm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tgps(rcpara%nmax,rcpara%pmax,rcpara%parmax),rasocorrhomd(rcpara%nmax,rcpara%pmax,rcpara%parmax,4)
    real(kind=JPRM) :: tfgmorig(rcpara%nmax,rcpara%pmax,rcpara%parmax),stfgmorig(rcpara%nmax,rcpara%pmax,rcpara%parmax),stmorig(rcpara%nmax,rcpara%pmax,rcpara%parmax),tgpsorig(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM) :: mrasocorrs(rcpara%mmax,rcpara%pmax,rcpara%parmax),trasocorrs(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM) :: eracorrs(rcpara%nmax,rcpara%pmax,rcpara%parmax),rasocorrs(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM),intent(in)  :: bgcorrs(rcpara%nmax,rcpara%pmax,rcpara%parmax),bgcorrs_e(rcpara%nmax,rcpara%pmax,rcpara%parmax),bgcorrs_w(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM),allocatable,dimension(:,:,:,:) :: arr
    integer,allocatable :: hours(:,:)
    real(kind=JPRM)  :: solarangles(rcpara%nmax,rcpara%parmax)
    real(kind=JPRM),intent(in)  :: densities(rcpara%nmax,wmostats)
    integer(kind=JPRM),intent(in)  :: era40meta(rcpara%nmax)

    real :: dists(rcpara %statmax)

    real(kind=JPRM),intent(in) :: ominuse40(rcpara%ni,rcpara%nj,rcpara%pmax,rcpara%parmax,12)
    real(kind=JPRM) :: adjust(rcpara%pmax,rcpara%parmax)

    real(kind=JPRM) :: crut2(72,36,rcpara%mmax),dailycrut2(rcpara%nmax,rcpara%parmax),rimax,ramax,rawmax

    integer,intent(in) :: wmonrs(rcpara%statmax)
    logical :: ini_correct(rcpara%parmax,rcpara%statmax),logcache,ini_correct2(rcpara%parmax)
    logical :: composite_exists(rcpara%parmax)
    real(kind=JPRM),intent(in) :: wmolats(rcpara%statmax),wmolons(rcpara%statmax)
    character*50,intent(in) :: wmonames(rcpara%statmax)
    !character*64 :: stgroupnames(:)
    character binprefix*20,startdate*8,cdatum*8,cstatnr*6,filename*100,cgstat*5,cj*3
    integer err,err1,err2,err3,cerr
    integer prob_method,protunit,blunit,ipresatend
    logical ex1,ex2,ex3
    integer rtype(rcpara%nmax),mr,cmax,icistat,iosmeta,j,thread_num,riloc(1),raloc(1),rawloc(1),ml(3)
    integer         :: lasts(rcpara%brmax,rcpara%parmax),midx(rcpara%brmax,rcpara%parmax),gcount(rcpara%parmax),lastsave(rcpara%brmax,rcpara%parmax),midxsave(rcpara%brmax,rcpara%parmax),gcountsave(rcpara%parmax)
    integer         :: tbicache(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax)
    integer          :: alarm(10000,3),needs_composite(:,:,:),needs_compositesave(rcpara%brmax,rcpara%parmax)
    logical :: iczwisch,lrem,lgalert,lalert
    character*5 :: names(4)
    integer:: bad_intervals(:,:)
    integer(kind=JPRM)  :: tnumbers(rcpara%nmax,rcpara%pmax,rcpara%parmax)


    ldebug=.false.
    null=0.

    write(*,*) 'statnr ',statnr

    iunit=20
    !$ thread_num=omp_get_thread_num()
    !$ iunit=iunit+thread_num ! this should avoid clashes in parallel runs

    write(cstatnr,'(I6.6)') statnr

    select case(rcpara%innov)
    case('RI')
       filename=trim(rcpara%prefix)//cstatnr//'/feedbackmerged'//cstatnr//'.nc'
    case('RE')
       filename=trim(rcpara%prefix)//cstatnr//'/feedbackmerged'//cstatnr//'.nc'
    case('RO')
       filename=trim(rcpara%prefix)//cstatnr//'/feedbackmerged_lo'//cstatnr//'.nc'
    case default
       write(*,*) 'innov must be either RI or RO but is ',rcpara%innov
       call exit(1)
    end select

    if(.not. ex(istat) ) then
       print*,trim(filename),' not found'
       return
    endif

    iosmeta=1
    if(iter .gt. 1) then
    endif
    !select case(rcpara%innov)
    !case('RI')
    !    filename=trim(rcpara%prefix)//cstatnr//'/feedbackmerged'//cstatnr//'.nc'
!!$ call omp_set_lock(omp_lp(rcpara%statmax+2))
    !$OMP CRITICAL
    !    logcache=icache(istat) .eq. 0
    logcache=.true.
!!$ call omp_unset_lock(omp_lp(rcpara%statmax+2))
    !$OMP END CRITICAL
!!$ call omp_set_lock(omp_lp(istat))
    if(logcache) then
       !      icache(rcpara%statmax+1)=icache(rcpara%statmax+1)+1
       !      icache(istat)=icache(rcpara%statmax+1)

       tfgm=rcpara%miss_val
       tm=tfgm
       eracorrs=tfgm
       if(rcpara%innov .ne. 'RE') then
          !MERRA      CALL read_odb_nc(filename,rcpara,istat,err3,tm,tfgm,tbcm=tbcm) !subroutine in read_txt_write_nc.f90
          CALL read_odb_nc(filename,rcpara,istat,err3,tm,tfgm,tanm=tanm,stname=wmonames(istat),bad_intervals=bad_intervals,tgps=tgps) !subroutine in read_txt_write_nc.f90
          do ipar=1,rcpara%parmax
             do ip=1,rcpara%pmax
                do i=1,rcpara%nmax
                   if (isnan(tm(i,ip,ipar)) .or. abs(tfgm(i,ip,ipar))>20.) then
                      tm(i,ip,ipar)=rcpara%miss_val
                      tfgm(i,ip,ipar)=rcpara%miss_val
                      tanm(i,ip,ipar)=rcpara%miss_val
                      tgps(i,ip,ipar)=rcpara%miss_val
                   else
                      if (tgps(i,ip,ipar) .ne. rcpara%miss_val) then
                         tgps(i,ip,ipar)=tgps(i,ip,ipar)-tfgm(i,ip,ipar)
                      endif
                      !                         tfgm(i,ip,ipar)=-tfgm(i,ip,ipar)
                      !                         tanm(i,ip,ipar)=-tanm(i,ip,ipar)
                   endif

                enddo
             enddo
          enddo
       else
       endif
       filename=trim(rcpara%prefix)//cstatnr//'/feedbackbgcorr'//cstatnr//'.nc'
       CALL read_bgcorr_nc(filename,rcpara,istat,err3,eracorrs,(/'bg correctio',''/)) !subroutine in read_txt_write_nc.f90
       if(rcpara%innov .ne. 'RI') eracorrs=0.
       !  if(any(abs(tm) .gt. 1.e10).or. any(abs(tfgm) .gt. 1.e10)) then
       !    stop 'invalid value'
       !  endif
       filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbinremsave'//cstatnr//'.nc'
       CALL read_sonde_corr_daily_nc(filename, rcpara,0, err,mrasocorrs, tbindex,tbi)
       if(err .eq. 0) then
          do ipar=1,rcpara%parmax
             do ip=1,rcpara%pmax
                do i=1,tbi
                   trasocorrs(tbindex(i)+1:tbindex(i+1),ip,ipar)=mrasocorrs(i,ip,ipar)         
                enddo
                trasocorrs(1,ip,ipar)=mrasocorrs(1,ip,ipar) 
                where(tm(:,ip,ipar) .ne. rcpara%miss_val .and. tfgm(:,ip,ipar) .ne. rcpara%miss_val )
                   tm(:,ip,ipar)=tm(:,ip,ipar)-trasocorrs(:,ip,ipar)
                   tfgm(:,ip,ipar)=tfgm(:,ip,ipar)+trasocorrs(:,ip,ipar)!-eracorrs(:,ip,ipar)
                endwhere
             enddo
          enddo
          lrem=.true.
       else
          where(tfgm .ne. rcpara%miss_val .and. eracorrs .ne. rcpara%miss_val )
             tfgm=tfgm-eracorrs
          elsewhere
             tfgm=rcpara%miss_val
          endwhere
          lrem=.false.
       endif

       do i=1,rcpara%nmax
          if(.not. any(tm(i,12,:) .ne. rcpara%miss_val) .or. .not. any(tfgm(i,12,:) .ne. rcpara%miss_val)) then
             tm(i,:,:)=rcpara%miss_val
             tfgm(i,:,:)=rcpara%miss_val
          endif
       enddo

       icistat=icache(istat)
    else

       if(ldebug) then
          filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbincorrsave'//cstatnr//'.nc'
          CALL read_sonde_corr_daily_nc(filename, rcpara,0, err,mrasocorrs, tbindex,tbi)
          if(err .eq. 0) then
             do ipar=1,rcpara%parmax
                do ip=1,rcpara%pmax
                   do i=1,tbi
                      rasocorrs(tbindex(i)+1:tbindex(i+1),ip,ipar)=mrasocorrs(i,ip,ipar)         
                   enddo
                   rasocorrs(1,ip,ipar)=mrasocorrs(1,ip,ipar) 
                enddo
             enddo
          endif
       endif
       err3=0
       tfgm=rcpara%miss_val
       tm=tfgm
       !$OMP CRITICAL
       icistat=icache(istat)
       do ipar=1,rcpara%parmax
          do ip=1,rcpara%pmax
             do i=1,tfgmcr(icistat)%vals
                l=tfgmcr(icistat)%index(i)
                tm(l,ip,ipar)=tmcr(icistat)%feld(i,ip,ipar)
                tfgm(l,ip,ipar)=tfgmcr(icistat)%feld(i,ip,ipar)
             enddo
          enddo
       enddo
       !$OMP END CRITICAL
       if(any(abs(tm) .gt. 1.e10) .or. any(abs(tfgm) .gt. 1.e10)) then
          write(*,*) 'invalid value'
          call abort
       endif
    endif
!!$ call omp_unset_lock(omp_lp(istat))
    call eiminuse40(rcpara,wmolats(istat),wmolons(istat),ominuse40,adjust)
    switchindex=toindex(rcpara%switchdate,rcpara)
    !    call rtypeinfo(rtype,rcpara,mr,istat,imax) !in file rfcorio.f90 line 847
    err=err3 !*err1*err2
    !end select
    if(err .ne. 0) then
       write(*,*) 'found but could not read ',filename
       return
    endif

    blunit=300
    !$ blunit=blunit+thread_num
    filename=trim(rcpara%prefix)//'raso_blacklist' 
    !call read_raso_blacklist(blunit,filename,rcpara,statnr,tm,tfgm,err) !in file rfcorio.f90 line 245

    call check_ini(meta_s,istat,wmonrs,rcpara,ini_correct)
    protunit=100
    !$ protunit=protunit+thread_num


    call detect_gaps(rcpara,statnr,tfgm,midx,lasts,gcount) !in file rfcor.f90 line 113
    if (.not. any(gcount>20)) then
       needs_composite(:,:,istat)=lasts(1:20,:)
    else
       write(*,*) statnr,' too many gaps',gcount
    endif

    composite_exists=.false.
    if(any(gcount .gt. 0)) then 

       cmax=wmostats
       composite_exists=.false.
       if(ldebug)  then
          write(*,*) 'gcount',statnr,gcount
       endif

       call make_hom_composite2(rcpara,statnr,cmax,wmonrs,wmolons,wmolats,wmostats,dists,index,ominuse40,adjust,crut2,rtype,solarangles,tm,tfgm,tgps,tfgmcr,&
            tmcr,tccr,icache,meta_s,lastsave,midx,gcount,lastscache,gcountcache,tbicache,needs_composite,ini_correct,composite_exists,rasocorrhomd,eracorrs,ex,iter,thread_num,lrem,alarm,bad_intervals) !in file rfcomp2.f90 line 13

       if(iter .eq. rcpara%maxiter) then
          names=(/'nm___','rit'//rcpara%ens,'rgm__','rio'//rcpara%ens/)

          tfgmorig=tfgm
          tgpsorig=tgps

          do j=2,4,2
             tfgm=tfgmorig
             tgps=tgpsorig
             istart=floor((rcpara%mmax/12-2)*365.25) 
             write(*,*) 'late break after',istart
             do ipar=1,rcpara%parmax
                do i=istart,rcpara%nmax
                   if(rasocorrhomd(i,8,ipar,j) .ne. rasocorrhomd(i-1,8,ipar,j) .and. rasocorrhomd(i-1,8,ipar,j) .ne. 0) then
                      ini_correct(ipar,istat)=.true.
                      write(*,*) cstatnr,'late break'
                      exit
                   endif
                enddo
             enddo
             !! correct most recent part of series
             !! do this after breakpoint correction since then a longer interval for more accurate
             !! estimation can be used.
             call make_composite(rcpara,statnr,cmax,wmonrs,wmolons,wmolats,wmostats,dists,index,ominuse40,adjust,stm,stfgm,tnumbers,&
                  tfgmcr,tmcr,icache,meta_s,lasts,gcount,needs_composite,ini_correct,composite_exists,bad_intervals,tccr) !
            
             
             write(cj,'(A1,I1,A1)') '_',j,'_'
             if(iter .eq. rcpara%maxiter .and. any(composite_exists) .and. rcpara%extended_output .eq. 'Y') then
                filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbincompmon'//cj//cstatnr//'.nc'
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
                filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbincomp'//cj//cstatnr//'.nc'
                allocate(arr(tfgmcr(icistat)%vals,rcpara%pmax,rcpara%parmax,2),hours(tfgmcr(icistat)%vals,rcpara%parmax))
                do i=1,tfgmcr(icistat)%vals
                   arr(i,:,:,1)=stm(tfgmcr(icistat)%index(i),:,:)
                   arr(i,:,:,2)=stfgm(tfgmcr(icistat)%index(i),:,:)
                enddo
                CALL txttonc(rcpara,istat,filename,arr,(/'comp-obs    ','comp-bg     '/),tfgmcr(icistat)%vals, 2,tfgmcr(icistat)%index, hours, wmolats(istat),wmolons(istat),0,stname=wmonames(istat))
                deallocate(hours,arr)
             endif

             if (rcpara%initial_adjust(1:3).eq. 'all') then
                ini_correct2=.true.
             else
                ini_correct2=ini_correct(:,istat)
             endif

             do ipar=1,rcpara%parmax
               do i=1,gcount(ipar)
                 do ip=1,rcpara%pmax
                   rasocorrhomd(1:lasts(i,ipar)+rcpara%snht_maxlen/2,ip,ipar,j)=rasocorrhomd(1:lasts(i,ipar)+rcpara%snht_maxlen/2,ip,ipar,j)-rasocorrhomd(lasts(i,ipar)+rcpara%snht_maxlen/2,ip,ipar,j)
                enddo
               enddo
             enddo
             needs_compositesave= needs_composite(:,:,istat)
             lastsave=lasts
             midxsave=midx
             gcountsave=gcount
             where (stfgm .ne. rcpara%miss_val .and. tfgm .ne. rcpara%miss_val )
                stfgmsave=tfgm-stfgm-rasocorrhomd(:,:,:,j)
!                stfgmsave=tfgm-rasocorrhomd(:,:,:,j)
             elsewhere
                stfgmsave=rcpara%miss_val
             endwhere
             where (tfgm .ne. rcpara%miss_val)
                tfgmsave=tfgm-rasocorrhomd(:,:,:,j)
             elsewhere
                tfgmsave=rcpara%miss_val
             endwhere
             where (tgps .ne. rcpara%miss_val)
                tgpssave=tgps-rasocorrhomd(:,:,:,j)
             elsewhere
                tgpssave=rcpara%miss_val
             endwhere
             trasocorrs=-rasocorrhomd(:,:,:,j)
             call correct_mostrecent(rcpara,wmonrs(istat),midxsave,lastsave,gcountsave,tfgmsave,stfgmsave,tgpssave,needs_compositesave,ini_correct2,trasocorrs,rasobreaks,rasobreakuncertainties) !in this file line681
              rasocorrhomd(:,:,:,j)=-trasocorrs

!             rasocorrhomd(:,:,:,j-1)=rasocorrhomd(:,:,:,j)+trasocorrs

!             if (j==2) then
!                rasocorrhomd(:,:,:,j)=-trasocorrs
!             else
!                rasocorrhomd(:,:,:,j)=-trasocorrs+(rasocorrhomd(:,:,:,j-1)-rasocorrhomd(:,:,:,j-3))
!             endif


             if (j==4) then
                needs_composite(:,:,istat)=needs_compositesave
                tbi=1
                tbindex(1)=1
                l=1
                do i=1,rcpara%nmax
                   if (rasocorrhomd(i-1,9,1,2).ne. rasocorrhomd(i,9,1,2) .or. rasocorrhomd(i-1,9,2,2).ne. rasocorrhomd(i,9,2,2)) then
                      l=l+1
                      tbindex(l)=i
                      tbi=l
                   endif
                enddo
    if(icistat .gt. 0) then
!$ call omp_set_lock(omp_lp(icistat))
                lastscache(:,:,icistat)=needs_composite(:,:,istat)
                gcountcache(:,icistat)=gcountsave
       
       tccr(icistat,2)%vals=tbi
       tccr(icistat,3)%vals=tbi
       tccr(icistat,2)%index(1:tbi)=tbindex
       tccr(icistat,3)%index(1:tbi)=tbindex
       do bi=1,tbi
!          if(iter .eq. 2) then
             !      tccr(icistat,1)%feld(bi,:,:)=rasocorrhomd(tbindex(bi),:,:,2)
!          else
             tccr(icistat,1)%vals=0
             tccr(icistat,2)%feld(bi,:,:)=rasocorrhomd(tbindex(bi),:,:,4)
             tccr(icistat,3)%feld(bi,:,:)=rasocorrhomd(tbindex(bi),:,:,2)
!          endif
       enddo
       tccr(icistat,2)%feld(tbi+1,:,:)=0
       tccr(icistat,3)%feld(tbi+1,:,:)=0

!$ call omp_unset_lock(omp_lp(icistat))
    else
       write(*,*) 'error ',wmonrs(index(1))
    endif
            endif

             filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbincorrsave_'//names(j)//'_'//cstatnr//'.nc'
             call write_sonde_corr_daily_nc(filename,rcpara,istat,err,rasocorrhomd(:,:,:,j),wmolons(istat),wmolats(istat),stname=wmonames(istat)) !in file rfcorio.f90 line 1175
             filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbincorrsave_'//names(j)//'bt'//cstatnr//'.nc'
             where(tm .ne. rcpara%miss_val)
                tbcm=tm-rasocorrhomd(:,:,:,j)
             elsewhere
                tbcm=rcpara%miss_val
             endwhere
#ifdef RTTOV
             call write_sonde_monthly_bt_nc(filename,rcpara,tbcm,istat,err,wmolons(istat),wmolats(istat),4,skin,stname=wmonames(istat))
             filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbincorrsave_'//names(j)//'bt2'//cstatnr//'.nc'
             call write_sonde_monthly_bt_nc(filename,rcpara,tbcm,istat,err,wmolons(istat),wmolats(istat),2,skin,stname=wmonames(istat))
#endif
             if(ldebug) then
                ! check adjustment quality
                lgalert=.false.
                do ipar=1,rcpara%parmax
                   do ip=4,14
                      call snhteqsamp2(tfgm(:,ip,ipar),null,rcpara%nmax,1,rcpara%nmax,rcpara%snht_maxlen,rcpara%snht_increment,rcpara%miss_val,rcpara%max_miss,critical_dates,0, & 
                           tsa(:,ip,ipar,j/2),plus,minus,prms,mrms,philf,mhilf,rcpara%month)

                      where(tfgm(:,ip,ipar) .ne. rcpara%miss_val)
                         tbcm(:,ip,ipar)=tfgm(:,ip,ipar)-rasocorrhomd(:,ip,ipar,j)
                      elsewhere
                         tbcm(:,ip,ipar)=rcpara%miss_val
                      endwhere
                      call snhteqsamp2(tbcm(:,ip,ipar),null,rcpara%nmax,1,rcpara%nmax,rcpara%snht_maxlen,rcpara%snht_increment,rcpara%miss_val,rcpara%max_miss,critical_dates,0, & 
                           tsa2(:,ip,ipar,j/2),plus,minus,prms,mrms,philf,mhilf,rcpara%month)

                      where(tfgm(:,ip,ipar) .ne. rcpara%miss_val)
                         tbcm(:,ip,ipar)=tfgm(:,ip,ipar)-rasocorrs(:,ip,ipar)
                      elsewhere
                         tbcm(:,ip,ipar)=rcpara%miss_val
                      endwhere

                      call snhteqsamp2(tbcm(:,ip,ipar),null,rcpara%nmax,1,rcpara%nmax,rcpara%snht_maxlen,rcpara%snht_increment,rcpara%miss_val,rcpara%max_miss,critical_dates,0, & 
                           tsa3(:,ip,ipar,j/2),plus,minus,prms,mrms,philf,mhilf,rcpara%month)

                      lalert=.false.
                      do i=1,rcpara%nmax
                         if(tsa(i,ip,ipar,j/2) .ne. rcpara%miss_val .and. tsa2(i,ip,ipar,j/2) .ne. rcpara%miss_val) then
                            if(tsa2(i,ip,ipar,j/2) .gt. 300 .and. .not. lalert) then
                               rawmax=maxval(tsa(:,ip,ipar,j/2))
                               rawloc=maxloc(tsa(:,ip,ipar,j/2))
                               ramax=maxval(tsa3(:,ip,ipar,j/2))
                               raloc=maxloc(tsa3(:,ip,ipar,j/2))
                               rimax=maxval(tsa2(:,ip,ipar,j/2))
                               riloc=maxloc(tsa2(:,ip,ipar,j/2))
                               write(*,'(A5,I3,I2,I2,A12,3(I6,I11,F7.1))') cstatnr,ip,ipar,j,'SNHT alarm ',riloc(1),todate(riloc(1),rcpara),rimax,rawloc(1),todate(rawloc(1),rcpara),rawmax,raloc(1),todate(raloc(1),rcpara),ramax
                               lalert=.true.
                               lgalert=.true.
                            endif
                         endif
                      enddo
                   enddo
                enddo

                if(lgalert) then 
                   ml=maxloc(tsa2(:,:,:,j/2))
                   write(*,'(A5,I3,A16,3F7.1,3I6)') cstatnr,j,'SNHT max alarm ',tsa2(ml(1),ml(2),ml(3),j/2),tsa(ml(1),ml(2),ml(3),j/2),tsa3(ml(1),ml(2),ml(3),j/2),ml
                endif

             endif
          enddo
       else
          names=(/'n1___','n1ra_','rg1__','rg1ra'/)
          do j=2,4,2
             filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbincorrsave_'//names(j)//'_'//cstatnr//'.nc'
             !      call write_sonde_corr_daily_nc(filename,rcpara,istat,err,rasocorrhomd(:,:,:,j),wmolons(istat),wmolats(istat)) !in file rfcorio.f90 line 1175
          enddo
       endif


    endif


    return
  end subroutine comp_correct


end module comp_c
