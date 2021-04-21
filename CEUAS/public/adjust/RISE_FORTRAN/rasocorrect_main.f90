subroutine  rasocorrect_main(rcpara)

  use rfmod, only: rasocor_namelist,cacherecord,JPRB
  use rfcor 
  use lgribread
  use txtnc
  use rcorrect
  !  use rcorrect_lo
  use comp_c
  use tosat
  use netcdf

  implicit none

  integer,parameter :: nens=10

  type(rasocor_namelist) rcpara,rcpara1,rcparas(nens)
  type(ecskin) :: gskin,skin

  type(metadata) :: meta_s

  integer index,indexmax,statnr,wmostats,cwmostats,max_miss,i,j,k,l,istat,gstat,probmax,ilmeta,ip,ipar
  integer err,mmax,first,last,estatmax,estatindex,ewmostats,wmonr,istatmax,iwmostats
  integer iyear,imonth,iday,itime,cstatmax,imon,bi
  integer llla,lllo,urla,urlo
  integer ni,nj,nk,ic,iens,indstat,iswap

  integer(kind=JPIM),allocatable:: stgroups(:,:),hilfmeta(:),indexmeta(:),cardsmeta(:,:,:),ocardsmeta(:,:,:),era40meta(:)
  character startdate*8,cdatum*8,cstatnr*6,cstatnr2*5, cgstat*5, filename*120, oper*1,cens*2

  integer, allocatable:: wmonrs(:),ewmonrs(:),iwmonrs(:),cwmonrs(:),wmoheights(:)
  logical, allocatable:: ini_correct(:,:)
  real(kind=JPRM), allocatable:: wmolats(:),wmolons(:),ewmolats(:),ewmolons(:),iwmolats(:),iwmolons(:),cwmolats(:),cwmolons(:)
  character*50, allocatable:: wmonames(:),ewmonames(:),iwmonames(:),cwmonames(:)
  character*64, allocatable :: stgroupnames(:)
  real(kind=JPRM), allocatable:: bgcorrs(:,:,:,:),bgcorrs_w(:,:,:,:),bgcorrs_e(:,:,:,:)
  real(kind=JPRM), allocatable:: ifs_rasocorrs(:,:,:,:),ifs_rasobreaks(:,:,:,:),ifs_rasobreakuncertainties(:,:,:,:)
  real(kind=JPRM), allocatable:: ifs_rasocorrsNASA(:,:,:,:),ifs_rasobreaksNASA(:,:,:,:),ifs_rasobreakuncertaintiesNASA(:,:,:,:)
  integer(kind=JPIM), allocatable:: ifs_index(:,:),ifs_breakindex(:,:),needed(:),icache(:)
  integer(kind=JPIM), allocatable:: ifs_indexNASA(:,:),ifs_breakindexNASA(:,:)
  logical lcorr,omp_test_lock
  real(kind=JPRM), allocatable  ::   dists(:,:),density(:),densities(:,:),dsum(:)

  type(cacherecord),allocatable :: tfgmcr(:),tmcr(:),tccr(:,:)

  real(kind=JPRM) :: bias(16),h,han
  real(kind=JPRM),allocatable :: ominuse40(:,:,:,:,:),ominuse40p(:,:,:,:,:),ominuse40_an(:,:,:,:,:),ominuse40_anp(:,:,:,:,:),crut2(:,:,:)
  real           ,allocatable :: omhilf(:)
  real(kind=JPRM),allocatable :: crplus(:,:,:,:),crminus(:,:,:,:)

  real(kind=JPRM),allocatable :: bt(:,:),rswap(:)
  real(kind=JPRM),allocatable :: feldbox(:,:,:,:,:),feldclimbox(:,:,:,:,:)
  real(kind=JPRM),allocatable :: feldboxbg(:,:,:,:,:),feldclimboxbg(:,:,:,:,:)
  real(kind=JPRM),allocatable :: feldboxan(:,:,:,:,:),feldclimboxan(:,:,:,:,:)

  real(kind=JPRM),allocatable :: teclimgens(:,:,:,:,:),teclimgsens(:,:,:,:),teclimghilf(:,:,:,:,:),teclimgshilf(:,:,:,:)
  real(kind=JPRM),allocatable :: talat(:),talon(:)
  integer(kind=JPRM),allocatable :: bad_intervals(:,:)
  integer*4 :: ecsiz(7),ecensn,tastart
  integer*1 :: i2 
  integer :: lmmax,ipmax,iparmax,ii,jj,im,iplus,iminus,ish

  character*10 :: termin
  integer ::iter,iunit,iunit2,ios,cnum,iprof,ncid
  logical,allocatable :: ex(:)
  logical ex3,ex4
  character*2 CM,HH
  character*8 :: cpar(2)
  character*20 :: NameOfRoutine='rasocorrect_main'
#ifdef RTTOV
  Type(profile_Type), Allocatable   :: profiles(:)    
  Integer(Kind=jpim)                :: rttov_errorstatus  ! rttov error return code
  Integer(Kind=jpim)                :: nchanprof
  Integer(Kind=jpim) :: alloc_status(20)
  Integer(Kind=jpim) :: asw
  Integer(Kind=jpim) :: errorstatus
  Character (len=80) :: errMessage
#endif
  logical*1, allocatable :: used(:,:,:,:)
  integer,allocatable :: needs_composite(:,:,:)


  integer,allocatable         :: lastscache(:,:,:),gcountcache(:,:)
  integer,allocatable         :: tbicache(:,:,:,:)

  integer :: nprof,nlevels
  integer :: p, my_rank, ierr,status,alarm(10000,3)
  integer :: istart,istop,inc,ioffset,fgc
  real(kind=8),target,allocatable          :: firstgap(:)
  real(kind=8):: fswap
  integer,target,allocatable          :: fgindex(:)
  ! Initialize MPI, learn local rank and total number of processors.


  !!  wmostats=2881
  !!  filename=trim(rcpara%prefix)//'feedbackglobbinbgcorr'
  !!  allocate(bgcorrs(rcpara%nmax,rcpara%pmax,rcpara%parmax,1),densities(rcpara%nmax,wmostats),dsum(rcpara%nmax))
  !!  call read_sonde_daily(iunit,filename,rcpara,err,bgcorrs)

  istatmax=1536
  cstatmax=963
  estatmax=3070
  allocate(stgroups(rcpara%nmem,rcpara%ngroup),stgroupnames(rcpara%ngroup))
  allocate(wmonrs(rcpara%statmax),wmolons(rcpara%statmax),wmolats(rcpara%statmax),wmonames(rcpara%statmax),wmoheights(rcpara%statmax),ini_correct(rcpara%parmax,rcpara%statmax))
  allocate(ewmonrs(estatmax),ewmolons(estatmax),ewmolats(estatmax),ewmonames(estatmax))
  allocate(iwmonrs(istatmax),iwmolons(istatmax),iwmolats(istatmax),iwmonames(istatmax))
  allocate(cwmonrs(rcpara%statmax),cwmolons(rcpara%statmax),cwmolats(rcpara%statmax),cwmonames(rcpara%statmax))
  allocate(era40meta(rcpara%nmax))
  allocate(cardsmeta(rcpara%nmax,1,rcpara%nmeta),indexmeta(3000000),hilfmeta(3000000))
  allocate(meta_s.cardsmeta_s(rcpara%nmax,rcpara%statmax),meta_s.rsicodes(5000),meta_s.rscodes(5000),meta_s.trusted(5000))
  allocate(ifs_rasocorrs(rcpara%mmax,rcpara%pmax,rcpara%parmax,estatmax),ifs_index(rcpara%mmax,estatmax))
  allocate(ifs_rasobreaks(rcpara%mmax,rcpara%pmax,rcpara%parmax,estatmax),ifs_rasobreakuncertainties(rcpara%mmax,rcpara%pmax,rcpara%parmax,estatmax),ifs_breakindex(rcpara%mmax,estatmax))
  allocate(ifs_rasocorrsNASA(rcpara%brmax,rcpara%pmax,rcpara%parmax,estatmax),ifs_indexNASA(rcpara%brmax,estatmax))
  allocate(ifs_rasobreaksNASA(rcpara%brmax,rcpara%pmax,rcpara%parmax,estatmax),ifs_rasobreakuncertaintiesNASA(rcpara%brmax,rcpara%pmax,rcpara%parmax,estatmax),ifs_breakindexNASA(rcpara%brmax,estatmax))
  allocate(ex(rcpara%statmax))
  allocate(rswap(rcpara%ni))

  !do i=1,rcpara%statmax+3
  do i=1,8000
     !$ call omp_init_lock(omp_lp(i))
  enddo

  filename='./mergedstations.t'
  call erastations(wmonrs,wmolats,wmolons,rcpara%statmax,wmostats,filename) !in file rfcorio.f90 line 305
  rcpara%cachemax=wmostats+1


  filename=trim(rcpara%prefix)//'country.bin'
  !call read_country(20,filename,rcpara%nmem,rcpara%ngroup,stgroups,stgroupnames,err)

  iunit=20
  statnr=1001

  call read_cards_meta_schroeder(iunit,rcpara,statnr,wmonrs,wmolons,wmolats,wmostats,meta_s,err) !in file rfcorio.f90 line 2081
  if(err .ne. 0) call exit(1)

  ini_correct=.true.
  inquire(file='ini_correct',exist=ex(1))
  if(ex(1)) then
     open(iunit,file='ini_correct',action='read',form='unformatted')
     read(iunit) ini_correct
     close(iunit)
  endif
  if(err .ne. 0) then
    write(*,*) 'could not read ini_correct'
  endif !call exit(1)

  call read_bad_intervals('bad_intervals',rcpara,bad_intervals)

#ifdef RTTOV

  !msu setup
  nrttovid=1
  nprof=rcpara%ni*rcpara%nj
  errorstatus     = 0_jpim
  alloc_status(:) = 0_jpim
  allocate (instrument(3,nrttovid),stat= alloc_status(1))
  instrument(1,1)=1 ! platform NOAA
  instrument(2,1)=14 ! sat ID 14
  instrument(3,1)=1  ! instrument MSU
  !  instrument(2,1)=16 ! sat ID 16
  !  instrument(3,1)=3  ! instrument AMSU
  nchannels=3
  input_chan(1:nchannels)=(/2,3,4/) ! MSU 2,3,4
  !  input_chan(1:nchannels)=(/7,9/) ! AMSU 7,9
  input_ems(1:nchannels)=(/0.0,0.0,0.0/)
  
  call msu_fwd_init(instrument,0_jpim,0_jpim,0.0_jprm,input_chan,input_ems,rcpara%pmax-5)
#endif

  !  if(rcpara%innov .ne. 'RI') then

  allocate(crut2(72,36,rcpara%mmax))
  !  filename='../common/hadcrut3.dat'
  !  call read_hadCRUT3(iunit,filename,crut2,rcpara,'../common/abstem3.dat') !in file rfcorio.f90 line 172
  filename='../common/HadCRUT.4.6.0.0.median.nc'
  call read_hadCRUT4(filename,crut2,rcpara) !in file rfcorio.f90 line 172
  where(crut2 .lt. -1.e29) crut2=rcpara%miss_val

  !  endif

  oper="D"
  termin="0000000000"
  ni=rcpara%ni
  nj=rcpara%nj
  nk=rcpara%pmax
  allocate(ominuse40(rcpara%ni,rcpara%nj,rcpara%pmax,rcpara%parmax,12),ominuse40p(rcpara%ni,rcpara%nj,rcpara%pmax,rcpara%parmax,12),omhilf(ni*nj*nk))
  allocate(crplus(rcpara%ni,rcpara%nj,rcpara%pmax,rcpara%parmax),crminus(rcpara%ni,rcpara%nj,rcpara%pmax,rcpara%parmax))
  allocate(ominuse40_an(rcpara%ni,rcpara%nj,rcpara%pmax,rcpara%parmax,12))
  write(*,*) rcpara%ni,rcpara%nj,rcpara%pmax,rcpara%parmax

#ifdef RTTOV
  if(rcpara%startdate .eq. 19570101) then
     call ecskin_init(skin,1958,1957+(rcpara%mmax-1)/12,rcpara%startdate,rcpara%mmax,crut2,'../common/','_s',rcpara%miss_val)
  else
     call ecskin_init(skin,1958,rcpara%startdate/10000+(rcpara%mmax-1)/12,rcpara%startdate,rcpara%mmax,crut2,'../common/','_s',rcpara%miss_val)

  endif
#endif
  ominuse40_an=0.
  if (rcpara%switchdate .eq. 19790100 ) THEN  ! defacto disabled



     filename="../common/FAL_1979000000"
     open(iunit,file=filename,form='unformatted',status='old')
     read(iunit) ni,nj,nk,ipar,imon
     read(iunit) ominuse40
     close(iunit)
     filename="../common/AAL_1979000000"
     open(iunit,file=filename,form='unformatted',status='old')
     read(iunit) ni,nj,nk,ipar,imon
     read(iunit) ominuse40_an
     close(iunit)
  endif

  allocate(bgcorrs(rcpara%nmax,rcpara%pmax,rcpara%parmax,1),bgcorrs_e(rcpara%nmax,rcpara%pmax,rcpara%parmax,1),bgcorrs_w(rcpara%nmax,rcpara%pmax,rcpara%parmax,1),densities(rcpara%nmax,wmostats),dsum(rcpara%nmax))
  if(rcpara%bg_correction_factor .gt. 0) then
     selectcase(rcpara%bghom)
     case(1)  
        filename=trim(rcpara%prefix)//'feedbackglobbinbgcorr_hom'
     case(2)  
        filename=trim(rcpara%prefix)//'feedbackglobbinbgcorr_eqarea'
     case default
        filename=trim(rcpara%prefix)//'feedbackglobbinbgcorr'
     end select

     err=0
     if(rcpara%innov .ne. 'RI'.and. rcpara%innov .ne. 'RO'.and. rcpara%innov .ne. 'RE') call read_sonde_daily(iunit,filename,rcpara,err,bgcorrs) !in file rfcorio.f90 line 1009
     if(err .gt. 0) call exit(1)
     where(bgcorrs .eq. rcpara%miss_val) bgcorrs=0.
     filename=trim(rcpara%prefix)//'feedbackglobbinbgcorr_e'
     !  call read_sonde_daily(iunit,filename,rcpara,err,bgcorrs_e) !in file rfcorio.f90 line 1009
     if(err .gt. 0) call exit(1)
     filename=trim(rcpara%prefix)//'feedbackglobbinbgcorr_w'
     !  call read_sonde_daily(iunit,filename,rcpara,err,bgcorrs_w) !in file rfcorio.f90 line 1009
     if(err .gt. 0) call exit(1)
     rcpara1=rcpara
     call rcparapmaxparmax(rcpara1,wmostats,1) !in file rfmod.f90 line 206

     if(rcpara%statmax .ne. 3070) then
        filename=trim(rcpara%prefix)//'feedbackglobbindensities'
     else
        filename=trim(rcpara%prefix)//'feedbackglobbindensities_all'
     endif
     if(rcpara%innov .ne. 'RI' .and. rcpara%innov .ne. 'RE' .and. rcpara%innov .ne. 'LO' .and. rcpara%innov .ne. 'RO') then
        call read_sonde_daily(iunit,filename,rcpara1,err,densities) !in file rfcorio.f90 line 1009
        if(err .gt. 0) call exit(1)
        dsum=maxval(densities,dim=2)
        !$omp parallel do private(i)
        do istat=1,wmostats
           do i=1,rcpara%nmax
              if(dsum(i) .gt. 0) then 
                 !! weg am 27.10.2006      densities(i,istat)=rcpara%bg_correction_factor-densities(i,istat)/dsum(i)
                 densities(i,istat)=1.2-densities(i,istat)/dsum(i)
              else 
                 densities(i,istat)=0.
              endif
           enddo
        enddo
        !$omp end parallel do
        dsum=sum(densities,dim=2)/wmostats
        !$omp parallel do private(i,istat)
        do istat=1,wmostats
           do i=1,rcpara%nmax
              if(dsum(i) .gt. 0) then 
                 densities(i,istat)=densities(i,istat)/dsum(i)
              endif
           enddo
        enddo
        !$omp end parallel do
     endif

  endif


  allocate(lastscache(rcpara%brmax,rcpara%parmax,rcpara%cachemax),gcountcache(rcpara%parmax,rcpara%cachemax))
  allocate(tbicache(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax))

  allocate(icache(rcpara%statmax+1))
  allocate(tfgmcr(rcpara%cachemax),tmcr(rcpara%cachemax),tccr(rcpara%cachemax,3))
  do i=1,rcpara%cachemax
     !haim     allocate(tfgmcr(i)%index(rcpara%nmax))
     do j=1,3
        allocate(tccr(i,j)%index(rcpara%brmax))
        allocate(tccr(i,j)%feld(rcpara%brmax,rcpara%pmax,rcpara%parmax))
        tccr(i,j)%index=rcpara%nmax !! later exit criterion
        tccr(i,j)%vals=0
        tccr(i,j)%feld=0.
     enddo
     !haim     tfgmcr(i)%index=rcpara%nmax !! later exit criterion
     !! index for tmcr is identical, i.e. not needed
  enddo
  icache=0
  lastscache=0
  tbicache=0
  gcountcache=0

  if(rcpara%innov .eq. 'RI'.or. rcpara%innov .eq. 'RE'.or. rcpara%innov .eq. 'RO') then
     ex=.FALSE.
     do istat=wmostats,1,-1
        write(cstatnr,'(I6.6)') wmonrs(istat)
        IF (INDEX(cstatnr(1:1),'0') .NE. 0) THEN 
           cstatnr2=cstatnr(2:5)
        ELSE
           cstatnr2=cstatnr
        END IF
        if( rcpara%innov .eq. 'RO') then
           filename=trim(rcpara%prefix)//cstatnr//'/feedbackmerged_lo'//cstatnr//'.nc'
        else
           filename=trim(rcpara%prefix)//cstatnr//'/feedbackmerged'//cstatnr//'.nc'
        endif
        inquire(file=filename,exist=ex4)

        if(ex4) then
           filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbincorrsave'//cstatnr//'.nc'
           inquire(file=filename,exist=ex4)
        endif

        ex(istat)=ex4
     enddo
     !  deallocate(ominuse40_an) 

     allocate(firstgap(wmostats),fgindex(wmostats))

     alarm=0

     call read_alarm(iunit,'alarm',alarm)

     allocate(needs_composite(20,rcpara%parmax,rcpara%statmax))
     needs_composite=0

     open(iunit,file='tbicachedump',form='unformatted',iostat=ios,status='old')
     if(ios .eq. 0) then
        !    read(iunit) tbicache
        close(iunit)
     endif

     cnum=0
     open(iunit,file='comp_tmcachedump',form='unformatted',iostat=ios,status='old')
     if(ios .eq. 0) then
        read(iunit) i
        read(iunit) icache
        do istat=1,rcpara%statmax
           if(icache(istat) .gt. 0) then
              read(iunit) tfgmcr(icache(istat))%vals
              tmcr(icache(istat))%vals=tfgmcr(icache(istat))%vals
              allocate(tfgmcr(icache(istat))%index(tmcr(icache(istat))%vals))
              allocate(tfgmcr(icache(istat))%feld(tmcr(icache(istat))%vals,rcpara%pmax,rcpara%parmax))
              allocate(tmcr(icache(istat))%feld(tmcr(icache(istat))%vals,rcpara%pmax,rcpara%parmax))
              read(iunit) tfgmcr(icache(istat))%index
              read(iunit) tfgmcr(icache(istat))%feld  
              read(iunit) tmcr(icache(istat))%feld  

              do j=1,3 
                 read(iunit) tccr(icache(istat),j)%vals
                 if(tccr(icache(istat),j)%vals .gt. 0) then
                    tccr(icache(istat),j)%index=rcpara%nmax
                    read(iunit) tccr(icache(istat),j)%index(1:tccr(icache(istat),j)%vals)
                    !                 allocate(tccr(icache(istat),j)%feld(rcpara%brmax,rcpara%pmax,rcpara%parmax))
                    read(iunit) tccr(icache(istat),j)%feld  
                 endif
              enddo
           endif
        enddo
        close(iunit)
        cnum=icache(rcpara%statmax+1)
     endif

     iens=1
     rcparas(iens)=rcpara

     iter=1
     !     do while(iter .lt. rcpara%maxiter)
     !        iter=iter+1
     istart=1
     istop=wmostats
     inc=1


     write(*,*) rcpara%anfsonde,rcpara%endsonde
     !$OMP PARALLEL DO PRIVATE(statnr,gstat,istat,i) SCHEDULE(dynamic,1)
     do istat=istart,istop,inc
        !        do istat=1,wmostats

        if(.not. any(.not. ini_correct(:,istat))) cycle
        statnr=wmonrs(istat)

        gstat=0
        do i=1,rcpara%ngroup
           if(any(stgroups(:,i) .eq. statnr)) gstat=i
        enddo

        if(statnr .ge. rcpara%anfsonde .and. statnr .le. rcpara%endsonde  .and. abs(wmolats(istat)) .lt. 200.)  then !  .and. (wmolats(istat)) .gt. 60. ) then ! ) 

           call comp_correct(rcparas(iens),wmonrs,wmolats,wmolons,wmonames,wmostats,meta_s,era40meta,bgcorrs,bgcorrs_e,bgcorrs_w,densities,&
                ominuse40,crut2,ini_correct,tfgmcr,tmcr,tccr,lastscache,gcountcache,tbicache,icache,needs_composite,statnr,gstat,istat,ex,iter,alarm,skin,bad_intervals) !in file comp_correct.f90
        endif

     enddo
     !$OMP END PARALLEL DO

     ! then adjust the others
     !$OMP PARALLEL DO PRIVATE(statnr,gstat,istat,i) SCHEDULE(dynamic,1)
     do istat=istart,istop,inc
        !        do istat=1,wmostats

        if( any(.not. ini_correct(:,istat))) cycle
        statnr=wmonrs(istat)

        gstat=0
        do i=1,rcpara%ngroup
           if(any(stgroups(:,i) .eq. statnr)) gstat=i
        enddo

        if(statnr .ge. rcpara%anfsonde .and. statnr .le. rcpara%endsonde  .and. abs(wmolats(istat)) .lt. 200.)  then !  .and. (wmolats(istat)) .gt. 60. ) then ! ) 

           call comp_correct(rcparas(iens),wmonrs,wmolats,wmolons,wmonames,wmostats,meta_s,era40meta,bgcorrs,bgcorrs_e,bgcorrs_w,densities,&
                ominuse40,crut2,ini_correct,tfgmcr,tmcr,tccr,lastscache,gcountcache,tbicache,icache,needs_composite,statnr,gstat,istat,ex,iter,alarm,skin,bad_intervals) !in file comp_correct.f90
        endif

     enddo
     !$OMP END PARALLEL DO

     if (rcpara%maxiter > 1) then
     do istat=1,wmostats
        firstgap(istat)=rcpara%nmax
        do ipar=1,2
           l=1
           do while (needs_composite(l,ipar,istat)>1)
              firstgap(istat)=minval((/dble(needs_composite(l,ipar,istat)),firstgap(istat)/))
              l=l+1
           enddo
        enddo
     enddo
     fgc=count(firstgap<45000)
     if (fgc .gt. 5) then
        call qsort(firstgap,fgindex)
     else
        l=1
        do istat=1,wmostats
           if (firstgap(istat)<45000) then
              fgindex(l)=istat
              l=l+1
           endif
        enddo
     endif

     iter=2
     !$OMP PARALLEL DO PRIVATE(statnr,gstat,istat,i) SCHEDULE(dynamic,1)
     do indstat=wmostats,1,-1
        istat=fgindex(indstat)
        if (istat .eq. 0) cycle 
        statnr=wmonrs(istat)

        if(statnr .ge. rcpara%anfsonde .and. statnr .le. rcpara%endsonde) then
           call comp_correct(rcparas(iens),wmonrs,wmolats,wmolons,wmonames,wmostats,meta_s,era40meta,bgcorrs,bgcorrs_e,bgcorrs_w,densities,&
                ominuse40,crut2,ini_correct,tfgmcr,tmcr,tccr,lastscache,gcountcache,tbicache,icache,needs_composite,statnr,gstat,istat,ex,iter,alarm,skin,bad_intervals) !in file comp_correct.f90
        endif

     enddo !do 398
     !$OMP END PARALLEL DO

     endif

     if(.false. .and. cnum .lt. icache(rcpara%statmax+1)) then  ! overwrite cachefile only if cache has changed
        open(iunit,file='comp_tmcachedump',form='unformatted')
        write(iunit) rcpara%statmax
        write(iunit) icache
        do istat=1,rcpara%statmax
           if(icache(istat) .gt. 0) then 
              write(iunit) tfgmcr(icache(istat))%vals
              write(iunit) tfgmcr(icache(istat))%index
              write(iunit) tfgmcr(icache(istat))%feld  
              write(iunit) tmcr(icache(istat))%feld  
              do j=1,3 
                 write(iunit) tccr(icache(istat),j)%vals
                 if( tccr(icache(istat),j)%vals .gt. 0) then
                    write(iunit) tccr(icache(istat),j)%index(1:tccr(icache(istat),j)%vals)
                    write(iunit) tccr(icache(istat),j)%feld  
                 endif
              enddo
           endif
        enddo
        close(iunit)
     endif

     !  open(iunit,file='tbicachedump',form='unformatted')
     !  write(iunit) tbicache
     !  close(iunit)
!     open(iunit,file='ini_correct',form='unformatted')
!     write(iunit) ini_correct
!     close(iunit)
     deallocate(needs_composite)

     stop
  else if (rcpara%innov .eq. 'NN' .or. rcpara%innov .eq. 'NR' .or. rcpara%innov .eq. 'NE') THEN
     !  deallocate(ominuse40_an) 

     allocate(needs_composite(20,rcpara%parmax,rcpara%statmax))
     allocate(firstgap(wmostats),fgindex(wmostats))


     iter=1
     ini_correct=.true.
     !$OMP PARALLEL DO PRIVATE(statnr,gstat,istat,i) SCHEDULE(dynamic,1)
     do istat=wmostats,1,-1 !do 398
        !        do istat=1,wmostats !do 398

        statnr=wmonrs(istat)

        gstat=0
        do i=1,rcpara%ngroup
           if(any(stgroups(:,i) .eq. statnr)) gstat=i
        enddo

        if(statnr .ge. rcpara%anfsonde .and. statnr .le. rcpara%endsonde) then

           call raso_correct(rcpara,wmonrs,wmolats,wmolons,wmonames,wmostats,indexmeta,hilfmeta,ilmeta,meta_s,era40meta,bgcorrs,bgcorrs_e,bgcorrs_w,densities,&
                ominuse40,ominuse40_an,crut2,ini_correct,tfgmcr,tmcr,icache,needs_composite,statnr,gstat,istat,iter,skin,bad_intervals) !in file raso_correct.f90 
        endif

     enddo !do 398
     !$OMP END PARALLEL DO

     do istat=1,wmostats
        firstgap(istat)=rcpara%nmax
        do ipar=1,2
           l=1
           do while (needs_composite(l,ipar,istat)>1)
              firstgap(istat)=minval((/dble(needs_composite(l,ipar,istat)),firstgap(istat)/))
              l=l+1
           enddo
        enddo
     enddo
     fgc=count(firstgap<45000)
     if (count(firstgap<45000) .gt. 1) then
!        call qsort(firstgap,fgindex)
!! bubble sort firstgap
        do i=1,wmostats
          fgindex(i)=i
        enddo
        do i=1,wmostats
          do j=i,wmostats
            if(firstgap(j) .lt. firstgap(i)) then
               fswap=firstgap(j)
               firstgap(j)=firstgap(i)
               firstgap(i)=fswap
               iswap=fgindex(j)
               fgindex(j)=fgindex(i)
               fgindex(i)=iswap
            endif
          enddo
        enddo
     else
        do istat=1,wmostats
           if (firstgap(istat)<45000) then
              fgindex(1)=istat
           endif
        enddo
     endif

     iter=3
     !$OMP PARALLEL DO PRIVATE(statnr,gstat,istat,i) SCHEDULE(dynamic,1)
     do indstat=wmostats,1,-1
        istat=fgindex(indstat)
        if (istat .eq. 0) cycle 
        statnr=wmonrs(istat)

        if(statnr .ge. rcpara%anfsonde .and. statnr .le. rcpara%endsonde) then

           call raso_correct(rcpara,wmonrs,wmolats,wmolons,wmonames,wmostats,indexmeta,hilfmeta,ilmeta,meta_s,era40meta,bgcorrs,bgcorrs_e,bgcorrs_w,densities,&
                ominuse40,ominuse40_an,crut2,ini_correct,tfgmcr,tmcr,icache,needs_composite,statnr,gstat,istat,iter,skin,bad_intervals) !in file raso_correct.f90 
        endif

     enddo !do 398
     !$OMP END PARALLEL DO

     open(iunit,file='ini_correct',form='unformatted')
     write(iunit) ini_correct
     close(iunit)
     deallocate(needs_composite)

  else if (rcpara%innov .eq. 'LO') THEN
     !  deallocate(ominuse40_an) 

     allocate(needs_composite(20,rcpara%parmax,wmostats))
     needs_composite=.false.


     iter=0
     do while(iter .lt. rcpara%maxiter) !do 393
        iter=iter+1
        if(iter .eq. 2) iter=iter+1


        !$OMP PARALLEL DO PRIVATE(statnr,gstat,istat,i) SCHEDULE(dynamic,1)
        do istat=wmostats,1,-1 !do 398
           !        do istat=1,wmostats !do 398

           statnr=wmonrs(istat)

           gstat=0
           do i=1,rcpara%ngroup
              if(any(stgroups(:,i) .eq. statnr)) gstat=i
           enddo

           if(statnr .ge. rcpara%anfsonde .and. statnr .le. rcpara%endsonde) then

              !              call raso_correct_lo(rcpara,wmonrs,wmolats,wmolons,wmonames,wmostats,indexmeta,hilfmeta,ilmeta,meta_s,era40meta,bgcorrs,bgcorrs_e,bgcorrs_w,densities,&
              !                   ominuse40,ominuse40_an,crut2,ini_correct,tfgmcache,tmcache,tfgmcr,tmcr,icache,needs_composite,statnr,gstat,istat,iter,skin,bad_intervals) !in file raso_correct.f90 
           endif

        enddo !do 398
        !$OMP END PARALLEL DO
     enddo !do 393
     open(iunit,file='ini_correct',form='unformatted')
     write(iunit) ini_correct
     close(iunit)
     deallocate(needs_composite)

  else if(rcpara%innov .eq. 'MO') then
     !  deallocate(ominuse40_an) 

     allocate(needs_composite(10,rcpara%parmax,rcpara%statmax))
     needs_composite=.false.


     iter=0
     do while(iter .lt. 3)
        iter=iter+1
        if(iter .eq. 2) iter=iter+1
        !! this loop may be parallelized
        !$OMP PARALLEL DO PRIVATE(statnr,gstat,istat,i) SCHEDULE(dynamic,1)
        do istat=1,wmostats

           statnr=wmonrs(istat)

           gstat=0
           do i=1,rcpara%ngroup
              if(any(stgroups(:,i) .eq. statnr)) gstat=i
           enddo

           if(statnr .ge. rcpara%anfsonde .and. statnr .le. rcpara%endsonde) then

              call raso_correct(rcpara,wmonrs,wmolats,wmolons,wmonames,wmostats,indexmeta,hilfmeta,ilmeta,meta_s,era40meta,bgcorrs,bgcorrs_e,bgcorrs_w,densities,&
                   ominuse40,ominuse40_an,crut2,ini_correct,tfgmcr,tmcr,icache,needs_composite,statnr,gstat,istat,iter,skin,bad_intervals) !in file raso_correct.f90 
           endif

        enddo
        !$OMP END PARALLEL DO
     enddo
     open(iunit,file='ini_correct',form='unformatted')
     write(iunit) ini_correct
     close(iunit)
     deallocate(needs_composite)

  else if(rcpara%innov .eq. 'MI') then
     allocate(needed(rcpara%statmax))

     needed=0
     allocate(used(rcpara%nmax,rcpara%pmax,rcpara%parmax,1))
     used=.false.

     !! this loop may be parallelized
!!$!OMP PARALLEL DO PRIVATE(statnr,gstat,istat,i) SCHEDULE(dynamic,1)
     do istat=1,iwmostats

        statnr=iwmonrs(istat)

        if(statnr .ge. rcpara%anfsonde .and. statnr .le. rcpara%endsonde) then !! .or. statnr .eq. 01001 .or. statnr .eq. 48900 .or. statnr .eq. 94120) then

           call raso_correct_igra_ei(rcpara,wmonrs,wmolats,wmolons,wmonames,wmostats,iwmonrs,iwmolats,iwmolons,iwmonames,iwmostats,istatmax,ominuse40,used,needed,statnr,istat) !in file raso_correct_igra.f90 

        endif

     enddo
!!$!OMP END PARALLEL DO

     open(20,file='mergedstations.t')
     l=0
     do istat=1,iwmostats
        l=l+1
        write(iunit,'(I5,I4,1X,I5,F8.2,F8.2,I8)') l,0,iwmonrs(istat),iwmolats(istat),iwmolons(istat),0
     enddo

     do istat=1,wmostats
        if(wmonrs(istat) .ne. 0 .and. .not. any(iwmonrs .eq. wmonrs(istat)) .and. .not. any(wmonrs(1:istat-1) .eq. wmonrs(istat))) then 
           write(*,*) l,istat,needed(istat),wmonrs(istat)
           l=l+1

           if(needed(istat) .eq. 0) then
              write(iunit,'(I5,I4,1X,I5,F8.2,F8.2,I8)') l,istat,wmonrs(istat),wmolats(istat),wmolons(istat),0
              statnr=wmonrs(istat)
              call raso_correct_igra_ei(rcpara,wmonrs,wmolats,wmolons,wmonames,wmostats,wmonrs,wmolats,wmolons,wmonames,wmostats,rcpara%statmax,ominuse40,used,needed,statnr,istat) !in file raso_correct_igra.f90 
           else
              write(iunit,'(I5,I4,1X,I5,F8.2,F8.2,I8)') l,istat,wmonrs(istat),wmolats(istat),wmolons(istat),needed(istat)
           endif
        endif
     enddo

     close(20)

     stop
  else if(rcpara%innov .eq. 'MN') then
     allocate(needed(rcpara%statmax))

     needed=0
     !!allocate(used(rcpara%nmax,pmax,parmax,2*rcpara%statmax/3))
     allocate(used(rcpara%nmax,rcpara%pmax,rcpara%parmax,1))
     used=.false.

     !! this loop may be parallelized
!!$!OMP PARALLEL DO PRIVATE(statnr,gstat,istat,i) SCHEDULE(dynamic,1)
     do istat=1,iwmostats

        statnr=iwmonrs(istat)

        if(statnr .ge. rcpara%anfsonde .and. statnr .le. rcpara%endsonde) then !! .or. statnr .eq. 01001 .or. statnr .eq. 48900 .or. statnr .eq. 94120) then

           call raso_correct_igra_NASA(rcpara,wmonrs,wmolats,wmolons,wmonames,wmostats,iwmonrs,iwmolats,iwmolons,iwmonames,iwmostats,istatmax,ominuse40,used,needed,statnr,istat,'IGRA') !in file raso_correct_igra_NASA
        endif

     enddo
!!$!OMP END PARALLEL DO

     open(20,file='mergedstationsNASA.t')
     l=0
     do istat=1,iwmostats
        l=l+1
        write(iunit,'(I5,I4,1X,I5,F8.2,F8.2,I8)') l,0,iwmonrs(istat),iwmolats(istat),iwmolons(istat),0
     enddo

     do istat=1,wmostats
        if(wmonrs(istat) .ne. 0 .and. .not. any(iwmonrs .eq. wmonrs(istat)) .and. .not. any(wmonrs(1:istat-1) .eq. wmonrs(istat))) then 
           write(*,*) l,istat,needed(istat),wmonrs(istat)
           l=l+1

           if(needed(istat) .eq. 0) then
              write(iunit,'(I5,I4,1X,I5,F8.2,F8.2,I8)') l,istat,wmonrs(istat),wmolats(istat),wmolons(istat),0
              statnr=wmonrs(istat)
              call raso_correct_igra_NASA(rcpara,wmonrs,wmolats,wmolons,wmonames,wmostats,wmonrs,wmolats,wmolons,wmonames,wmostats,rcpara%statmax,ominuse40,used,needed,statnr,istat,'IGRA') !in file raso_correct_igra_NASA
           else
              write(iunit,'(I5,I4,1X,I5,F8.2,F8.2,I8)') l,istat,wmonrs(istat),wmolats(istat),wmolons(istat),needed(istat)
           endif
        endif
     enddo

     close(20)

     stop
  endif


  !!generate correction table to be read into IFS
  do istat=1,wmostats

     estatindex=0
     esm:  do j=1,estatmax
        if(wmonrs(istat) .eq. ewmonrs(j) .or.  wmolons(istat) .eq. ewmolons(j) .and. wmolats(istat) .eq. ewmolats(j)) then 
           estatindex=j
           exit esm
        endif
     enddo esm
     if(estatindex .eq. 0) then
        if(ldebug) write(*,*)  'no match found for station ',wmonrs(istat)
     else

!!$ call omp_set_lock(omp_lp)
        write(cstatnr,'(I6.6)') wmonrs(istat)
!!$ call omp_unset_lock(omp_lp)

        filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbincorrsave'//cstatnr//'.nc'
        err=1
        !  call read_sonde_corr_daily_IFS(20,filename,rcpara,estatindex,err,ifs_rasocorrs,ifs_index)
        call read_sonde_corr_daily_IFS_nc(filename,rcpara,istat,err,ifs_rasocorrs(:,:,:,estatindex),ifs_index(:,estatindex),bi,ifs_rasobreaks(:,:,:,estatindex),ifs_rasobreakuncertainties(:,:,:,estatindex)) !in file read_txt_write_nc.f90 line 694
        !  filename=trim(rcpara%prefix)//'feedbackglobbincorrsaveNASA'//cstatnr
        !  err=1
        !  call read_sonde_corr_daily_IFS(20,filename,rcpara,estatindex,err,ifs_rasocorrsNASA,ifs_indexNASA)

        !  filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbinbreaksave'//cstatnr//'.nc'
        err=1
        !  call read_sonde_corr_daily_IFS(20,filename,rcpara,estatindex,err,ifs_rasobreaks,ifs_breakindex,ifs_rasobreakuncertainties)
        !  call read_sonde_corr_daily_IFS_nc(filename,rcpara,err,ifs_rasobreaks(:,:,:,estatindex),ifs_breakindex(:,estatindex),ifs_rasobreakuncertainties(:,:,:,estatindex)) !in file read_txt_write_nc.f90 line 688

     endif

  enddo


  filename=trim(rcpara%prefix)//'leobiascor.t'
  call write_leobiascor_IFS_table(20,filename,rcpara,estatmax,ifs_rasocorrs,ifs_index,ewmonrs,ewmolats,ewmolons,err) !in file rfcorio.f90 line 1391

  !filename=trim(rcpara%prefix)//'biascorNASA.t'
  !call write_leobiascor_IFS_table(20,filename,rcpara,estatmax,ifs_rasocorrsNASA,ifs_indexNASA,ewmonrs,ewmolats,ewmolons,err)

  filename=trim(rcpara%prefix)//'leobiasbreaks.t'
  call write_leobiascor_IFS_table(20,filename,rcpara,estatmax,ifs_rasobreaks,ifs_breakindex,ewmonrs,ewmolats,ewmolons,err) !in file rfcorio.f90 line 1391

  filename=trim(rcpara%prefix)//'leobiasbreakuncertainties.t'
  call write_leobiascor_IFS_table(20,filename,rcpara,estatmax,ifs_rasobreakuncertainties,ifs_breakindex,ewmonrs,ewmolats,ewmolons,err) !in file rfcorio.f90 line 1391


  do i=1,rcpara%statmax+2
     !$ call omp_destroy_lock(omp_lp(i))
  enddo

  wmonr=94120
  iyear=186
  imonth=04
  iday=12
  itime=00
  !!call read_leobiascor_IFS_table(iunit,filename,wmonr,iyear,imonth,iday,itime,plevs,rcpara%pmax,lcorr,bias,err)
  !!write(*,*) lcorr,err
  !!write(*,'(16F6.2)') bias

  stop
end subroutine rasocorrect_main

