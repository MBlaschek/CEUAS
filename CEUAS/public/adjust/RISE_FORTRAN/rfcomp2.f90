module rfcomp2

  use rfmod
  use homtests
  use homtestsamp
  use rfcorio
  use rfcor
  use rfcomp_1
  use correct_mr

contains



  !obstau ist 2 fuer rich-tau, 4 fuer rich-obs
  !xpar ist 3-par wenn zu wenige nachbarn zum Zeitpunkt par da sind, sonst par
  subroutine ch2(cachehilf,tfgmcrvals,tfgmcrindex,tmcrfeld,tccrindex,tccrfeld,rasocorrhomd,rcpara,tbindex,ib,ip,ipar,istat,left_maxlen,right_maxlen,icistat,iter,sign)

    use rfmod
    implicit none
    integer :: i,l,m,first,last,isave
    integer,intent(in) :: ib,ip,ipar,istat,icistat,iter
    type(rasocor_namelist),intent(in) :: rcpara
    integer,intent(in) :: tbindex(rcpara%mmax)
    integer,intent(in) :: left_maxlen(rcpara%brmax,rcpara%parmax,rcpara%cachemax),right_maxlen(rcpara%brmax,rcpara%parmax,rcpara%cachemax)

    integer,intent(in) :: tfgmcrvals,tfgmcrindex(tfgmcrvals),tccrindex(rcpara%brmax)
    real(kind=JPRM),intent(in) ::tmcrfeld(tfgmcrvals),tccrfeld(rcpara%brmax)
    real(kind=JPRM) :: cachehilf(rcpara%nmax) ,sign
    real(kind=JPRM),intent(in) :: rasocorrhomd(rcpara%nmax)


    cachehilf=rcpara%miss_val
    first=tbindex(ib)-left_maxlen(ib,ipar,istat)
    last=tbindex(ib)+right_maxlen(ib,ipar,istat)
    if(last .gt. rcpara%nmax) then
       last=rcpara%nmax
    endif
    !write(*,*) 'first:',first,'last',last
    cachehilf(first:last)=rcpara%miss_val

    !$ call omp_set_lock(omp_lp(icistat))
    m=1
!!$OMP CRITICAL
    do i=1,tfgmcrvals
       l=tfgmcrindex(i)
       if(l .ge. first) exit
    enddo
    isave=i
    if (ib.eq.24 .and. ip.eq. 4 .and. ipar .eq. 1) then
      write(*,*) 'test here'
    endif
    if(iter .eq. 2 .and. tccrindex(1) .ne. rcpara%nmax) then
       do i=isave,tfgmcrvals
          l=tfgmcrindex(i)
          if(l .gt. last) exit
          if(tmcrfeld(i) .ne. rcpara%miss_val) then 
             cachehilf(l)=tmcrfeld(i)+rasocorrhomd(l)
             do while (tccrindex(m) .lt. l)
                m=m+1
             enddo
             if(tccrfeld(m) .ne. rcpara%miss_val) then
                cachehilf(l)=cachehilf(l)-tccrfeld(m) !-
             endif
          endif
       enddo
       i=0
    else
       do i=isave,tfgmcrvals
          l=tfgmcrindex(i)
          if(l .gt. last) exit
          if(tmcrfeld(i) .ne. rcpara%miss_val) then 
             cachehilf(l)=tmcrfeld(i)+rasocorrhomd(l)
          endif
       enddo
    endif
!!$OMP END CRITICAL
    !$ call omp_unset_lock(omp_lp(icistat))

    return
  end subroutine ch2

  subroutine make_hom_composite2(rcpara,statnr,cmax,wmonrs,wmolons,wmolats,wmostats,dists,index,ominuse40,adjust,crut2,rtype,solarangles,ttm,ttfgm,tgps,tfgmcr,tmcr,tccr,&
       icache,meta_s,lasts,midx,gcount,lastscache,gcountcache,tbicache,needs_composite,ini_correct,composite_exists,rasocorrhomd,eracorrs,ex,iter,thread_num,lrem,alarm,bad_intervals) !

    implicit none

    type(rasocor_namelist),intent(in) :: rcpara
    type(metadata) :: meta_s

    integer istat,i,j,statnr,err,iunit,cmax,imin,ib,minst,found(rcpara%brmax,rcpara%parmax),pfound(rcpara%brmax,rcpara%pmax,rcpara%parmax),l,istart,istop,left_lengtha,right_lengtha,rib,typ,k,min_count,bimax,isample,ifail(rcpara%parmax),imiss(rcpara%parmax),nbi,lcount,bisave,ipp,thread_num,bi2,istatmax,itccr,l1,l2,iorig,isubs,isave
    integer,intent(in) :: wmostats,wmonrs(rcpara%statmax),iter
    integer index(rcpara%statmax),sammem(rcpara%brmax,rcpara%cachemax),sammem1(rcpara%brmax,rcpara%cachemax),dunit

    real :: dists(rcpara%statmax),hilfdiff(rcpara%pmax),hilfdiff1(rcpara%pmax),distc(rcpara%cachemax),distarr(rcpara%cachemax),distarrl(rcpara%cachemax)
    real(kind=JPRM),intent(in) :: wmolons(rcpara%statmax),wmolats(rcpara%statmax)

    real(kind=JPRM) ,intent(in):: ominuse40(rcpara%ni,rcpara%nj,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM),intent(in) :: adjust(rcpara%pmax,rcpara%parmax)
    real(kind=JPRM) :: adjustlocal(rcpara%pmax,rcpara%parmax)
    real(kind=JPRM),intent(in)  :: solarangles(rcpara%nmax,rcpara%parmax)

    real(kind=JPRM) :: a,tfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM) :: ttm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tfg(rcpara%nmax,rcpara%pmax,rcpara%parmax),tanm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tbcm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM) :: ttfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax),e20cm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tfgmorig(rcpara%nmax,rcpara%pmax,rcpara%parmax),teststat(rcpara%brmax,rcpara%pmax,rcpara%parmax),tgps(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM) :: mrasocorrs(rcpara%mmax,rcpara%pmax,rcpara%parmax),xrasocorrs(rcpara%nmax,rcpara%pmax,rcpara%parmax)

    real(kind=JPRM) :: stfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax),stm(rcpara%nmax,rcpara%pmax,rcpara%parmax),stanm(rcpara%nmax,rcpara%pmax,rcpara%parmax),rasobreaks(rcpara%nmax,rcpara%pmax,rcpara%parmax),rasobreakuncertainties(rcpara%nmax,rcpara%pmax,rcpara%parmax)

    real(kind=JPRM) :: eracorrs(rcpara%nmax,rcpara%pmax,rcpara%parmax),rasocorrs(rcpara%mmax,rcpara%pmax,rcpara%parmax),weight,wsum,wcrit,trasocorrs(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM) :: rasocorrhomd(rcpara%nmax,rcpara%pmax,rcpara%parmax,4),r1,r2
    integer rtype(rcpara%nmax),mr,rad

    integer :: bindex(rcpara%mmax),bindex2(rcpara%mmax),bindexsave(rcpara%mmax),tbindex(rcpara%mmax),tbindex2(rcpara%mmax),mtbindex(rcpara%mmax),xbindex(rcpara%mmax),tbindexsave(rcpara%mmax),iana(rcpara%brmax),bi,tbi,tbi2,tbisave,ipmax,iparmax,ip,ipar,xbi,mtbi
    integer         :: tbicache(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax)
    real(kind=JPRM) :: thilfcorr(rcpara%mmax,rcpara%pmax,rcpara%parmax),tfghilfcorr(rcpara%mmax,rcpara%pmax,rcpara%parmax),stnumcrit(rcpara%brmax,rcpara%parmax),darr(rcpara%brmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM) :: xdrarr(rcpara%brmax,rcpara%pmax,rcpara%parmax),xdsarr(rcpara%brmax,rcpara%pmax,rcpara%parmax),xdrarr2(rcpara%brmax,rcpara%pmax,rcpara%parmax),obsxdrarr(rcpara%brmax,rcpara%pmax,rcpara%parmax),obsxdsarr(rcpara%brmax,rcpara%pmax,rcpara%parmax),obsdarr(rcpara%brmax,rcpara%pmax,rcpara%parmax),spaghilf(rcpara%cachemax)

    logical ex3,first,ex4,ex5,lbreak,adj,ldebug,jump,lrem
    logical :: ex(rcpara%statmax)
    integer err2,err3,ios,rios,ic(rcpara%parmax),omp_get_thread_num
    character cstatnr*6,filename*100,citer*1,ch*4
    logical composite_exists(rcpara%parmax),logcache
    integer,intent(in)         :: lasts(rcpara%brmax,rcpara%parmax),midx(rcpara%brmax,rcpara%parmax),gcount(rcpara%parmax)
    integer         :: lastscache(rcpara%brmax,rcpara%parmax,rcpara%cachemax),gcountcache(rcpara%parmax,rcpara%cachemax)
    integer         :: sidx(rcpara%brmax,rcpara%pmax,rcpara%parmax),needs_composite(:,:,:)
    logical,intent(in) :: ini_correct(rcpara%parmax,rcpara%statmax)
    real(kind=JPRM)  ::   tauspagarr(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax),obsspagarr(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax),null(rcpara%nmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM)  ::   tauspagarr2(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax),obsspagarr2(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax)

    integer left_maxlen(rcpara%brmax,rcpara%parmax,rcpara%cachemax),right_maxlen(rcpara%brmax,rcpara%parmax,rcpara%cachemax)

    type(cacherecord) :: tfgmcr(rcpara%cachemax),tmcr(rcpara%cachemax),tfgmc(rcpara%cachemax),tccr(rcpara%cachemax,3)

    integer         :: llasts(rcpara%brmax,rcpara%parmax),lgcount(rcpara%parmax),lmidx(rcpara%brmax,rcpara%parmax),imax,icache(rcpara%statmax+1),icistat,istati,chosenbreaks(30)

    real(kind=JPRM) :: tsa(rcpara%brmax,rcpara%pmax,rcpara%parmax),plus(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax),obsplus(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax)
    real(kind=JPRM) :: obssig(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax),tausig(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax)
    real(kind=JPRM) :: obssig2(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax),tausig2(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax)
    real(kind=JPRM) :: obspms(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax),obsmms(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax),taupms(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax)
    real(kind=JPRM) :: taumms(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax),tauplus(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax),minus(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax)
    real(kind=JPRM) :: obsminus(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax),tauminus(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax)
    real(kind=JPRM) :: prms(rcpara%brmax,rcpara%pmax,rcpara%parmax),mrms(rcpara%brmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM) :: obscorrm(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax),obscorrp(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax)
    real(kind=JPRM) :: taucorrm(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax),taucorrp(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax)

    real(kind=JPRM) :: cachetmp(rcpara%nmax),cachehilf(rcpara%nmax), cachehilf2(rcpara%nmax) ,mfak(rcpara%pmax)
    real(kind=JPRM) :: crut2(72,36,rcpara%mmax),dailycrut2(rcpara%nmax,rcpara%parmax),ratio

    integer :: pcount(rcpara%brmax,rcpara%pmax,rcpara%parmax),mcount(rcpara%brmax,rcpara%pmax,rcpara%parmax),cacheindex(rcpara%nmax),iplevs(rcpara%pmax)
    integer :: pcount2(rcpara%brmax,rcpara%pmax,rcpara%parmax),mcount2(rcpara%brmax,rcpara%pmax,rcpara%parmax)
    integer :: tpcount(rcpara%brmax,rcpara%pmax,rcpara%parmax),tmcount(rcpara%brmax,rcpara%pmax,rcpara%parmax),obsspagarrcount(rcpara%brmax,rcpara%pmax,rcpara%parmax)
    real(kind=JPRM) meanhilf,maxhilf,meanhilf1,maxhilf1,dhilf,xdr,edr,swap,dhilfc,xdrc,edrc,swapd
    integer enn,m,fpr,idx,ennc,spagindex(rcpara%brmax,rcpara%pmax,rcpara%parmax),sub,sub1,gson,n_critical,critical_dates(1),spagindexl(rcpara%brmax,rcpara%pmax,rcpara%parmax),halt,max_miss,max_miss1,breakdates(rcpara%brmax),full(rcpara%brmax,rcpara%parmax)
    real(kind=JPRM):: hilf,hilf13,sigma1,sigma2,sigma11,sigma21,sigma111,sigma211,pc,mc,diff,sqdiff
    real(kind=JPRM) :: apriori_probs(rcpara%nmax),profiles(rcpara%pmax,rcpara%parmax)
    logical :: lobsspagarr(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax),ltauspagarr(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax)
    logical :: lobsspagarr2(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax),ltauspagarr2(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax)
    integer ipstart,ipstop,inc,ig,igmax
    logical :: manual,lday,use_night,good,downwards,surfaceonly
    integer :: alarm(10000,3),ipresatend,eistart
    integer :: bad_intervals(:,:)

    character*2 cform


    manual=.false.
    ldebug=.false.
    lobsspagarr=.true.
    ltauspagarr=.true.
    lobsspagarr2=.true.
    ltauspagarr2=.true.


    write(cstatnr,'(I6.6)') statnr
    write(citer,'(I1)') iter
    dunit=0
    if(dunit .gt. 0) then 
       !$ dunit=dunit+omp_get_thread_num()
       open(dunit,file=cstatnr//'/'//'bstats_'//cstatnr//'_'//citer//'_'//rcpara%ens)
    endif

    index=0
    call select_composite(rcpara,statnr,cmax,wmonrs,wmolons,wmolats,wmostats,dists,index) !in file rfcor.f90 line 206
    prms=rcpara%miss_val
    mrms=rcpara%miss_val
    wsum=0.
    tauspagarr=rcpara%miss_val
    obsspagarr=rcpara%miss_val
    tauspagarr2=rcpara%miss_val
    obsspagarr2=rcpara%miss_val
    obsspagarrcount=0
    plus=rcpara%miss_val
    minus=rcpara%miss_val
    obsdarr=rcpara%miss_val
    xdrarr=rcpara%miss_val
    obsxdrarr=rcpara%miss_val
    minst=rcpara%minst
    tbindex=0

    min_count=rcpara%snht_maxlen/2-rcpara%rimax_miss
    if(iter .gt. 1)  then
       minst=3*minst
       min_count=min_count
    endif
    !  mfak=1
    !  write(*,'(a,14F6.2)') 'MFAK: ',mfak(1:14)


    iunit=200000
    iunit=iunit+thread_num
    hilf13=0.
    rasocorrhomd(:,:,:,2)=0.
    rasocorrhomd(:,:,:,4)=0.

    n_critical=1
    critical_dates=(/0/)
    teststat=rcpara%miss_val
    null=0.
    spagindex=rcpara%miss_val

    write(cstatnr,'(I6.6)') statnr

    filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbinremsave'//cstatnr//'.nc'
    CALL read_sonde_corr_daily_nc(filename, rcpara,index(1), err,rasocorrs, mtbindex,mtbi)
    if(err .ne. 0) mtbi=0
    filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbincorrsave'//cstatnr//'.nc'
    !    call read_firstiter(rcpara,iunit,filename,tbindex,tbi,trasocorrs,ios) !in file rfcorio.f90 line 18
    CALL read_sonde_corr_daily_nc(filename, rcpara,index(1), err,rasocorrs, tbindex,tbi)


    if(mtbi .gt. tbi) then
       do i=1,mtbi
          if(.not. any(tbindex .eq. mtbindex(i))) then
             call addbreak(rcpara,cstatnr,.true., &
                  mtbindex(i),tbindex,tbisave,tbi,ldebug)
          endif
       enddo
    endif
    do ipar=1,rcpara%parmax
       do ip=1,rcpara%pmax
          do i=1,tbi
             trasocorrs(tbindex(i)+1:tbindex(i+1),ip,ipar)=rasocorrs(i,ip,ipar)         
          enddo
          trasocorrs(1,ip,ipar)=rasocorrs(1,ip,ipar) 
       enddo
    enddo


    filename=filename(1:len(trim(filename))-12)//'monbt'//cstatnr//'.nc'
    inquire(file=filename,exist=ex3)
    if( .not. ex3) then 
       do ipar=1,rcpara%parmax
          do ip=1,rcpara%pmax
             do i=1,rcpara%nmax
                if(ttfgm(i,ip,ipar) .ne. rcpara%miss_val .and. ttm(i,ip,ipar) .ne. rcpara%miss_val .and. trasocorrs(i,ip,ipar) .ne. rcpara%miss_val ) then
                   tm(i,ip,ipar)=ttm(i,ip,ipar)-trasocorrs(i,ip,ipar)
                   tfg(i,ip,ipar)=tm(i,ip,ipar)+ttfgm(i,ip,ipar)
                else
                   tfg(i,ip,ipar)=rcpara%miss_val
                   tm(i,ip,ipar)=rcpara%miss_val
                endif
             enddo
          enddo
       enddo

       !     call write_sonde_monthly_bt_nc(filename,rcpara,tm,index(1),err,wmolons(index(1)),wmolats(index(1)),4,-trasocorrs,tfg,tfg)
       !     filename=filename(1:len(trim(filename))-8)//'2'//cstatnr//'.nc'
       !     call write_sonde_monthly_bt_nc(filename,rcpara,tm,index(1),err,wmolons(index(1)),wmolats(index(1)),2,-trasocorrs,tfg,tfg)

       !     filename=filename(1:len(trim(filename))-11)//cstatnr//'.nc'
       !     call write_sonde_monthly_nc(filename,rcpara,tm,index(1),err,wmolons(index(1)),wmolats(index(1)),-trasocorrs,tfg)

    endif

    do ipar=1,rcpara%parmax
       call expand(rcpara,crut2(floor((wmolons(index(1))+180.)/5)+1,floor((89.999-wmolats(index(1)))/5)+1,:) ,dailycrut2(:,ipar)) !in file rfcor.f90 line 88
    enddo
   

    tbindexsave=tbindex
    tbisave=tbi
    call modify_tbindex(rcpara,wmonrs,index,1,tbindex,tbi,needs_composite(:,:,index(1)),.true.,iter,meta_s,alarm,manual,.false.) !ldebug in file rfcomp_1.f90 line 321

    do ib=1,tbi
       write(*,*) cstatnr,' tbindex:',ib,tbindex(ib),todate(tbindex(ib),rcpara)
    enddo

    do ipar=1,rcpara%parmax
       do ib=2,tbi
          right_maxlen(ib,ipar,1)=tbindex(ib+1)-tbindex(ib)
          if(iter .eq. 2) right_maxlen(ib,ipar,1)=rcpara%nmax-2-tbindex(ib)
          right_maxlen(ib,ipar,1)=rcpara%nmax-2-tbindex(ib)
          if(any(tbindex(ib-1)-rcpara%snht_maxlen/2.eq. lasts)) then
             left_maxlen(ib,ipar,1)=tbindex(ib)-tbindex(ib-2)
          else
             left_maxlen(ib,ipar,1)=tbindex(ib)-tbindex(ib-1)
          endif
          if(iter .le. 1) then
             if(right_maxlen(ib,ipar,1) .gt. rcpara%mean_maxlen/2) right_maxlen(ib,ipar,1)=rcpara%mean_maxlen/2
             if(left_maxlen(ib,ipar,1) .gt. rcpara%mean_maxlen/2) left_maxlen(ib,ipar,1)=rcpara%mean_maxlen/2
          else
             if(right_maxlen(ib,ipar,1) .gt. rcpara%mean_maxlen) right_maxlen(ib,ipar,1)=rcpara%mean_maxlen
             if(left_maxlen(ib,ipar,1) .gt. rcpara%mean_maxlen) left_maxlen(ib,ipar,1)=rcpara%mean_maxlen
          endif

          mfak=2.0*rcpara%mean_maxlen/330/20
          do i=1,rcpara%pmax-2
             mfak(i)=mfak(i)*(1+(2*i-1.0)/rcpara%pmax)
             if(mfak(i) .gt. 2*mfak(1)) mfak(i)=2*mfak(1)
          enddo
          ! at end of series there may be only few data after break - in this case reduce mfak
          do i=tbindex(ib)+right_maxlen(ib,ipar,1),tbindex(ib),-1
             if(ttfgm(i,8,ipar) .ne. rcpara%miss_val) exit
          enddo
          if(tbindex(ib)+1.5*mfak(8)*min_count .gt. i .and. tbindex(ib)+min_count .lt. i ) then
             mfak(:)=1.
          endif
          mfak=1.0
          do ip=1,rcpara%pmax-2
             call meaneqsamp3(ttfgm(:,ip,ipar),null(:,ip,ipar),rcpara%nmax,tbindex(ib)-left_maxlen(ib,ipar,1),left_maxlen(ib,ipar,1),right_maxlen(ib,ipar,1),&
                  rcpara%snht_increment,rcpara%miss_val,min_count,mfak(ip),critical_dates,n_critical,plus(ib,ip,ipar,1),minus(ib,ip,ipar,1),tpcount(ib,ip,ipar),tmcount(ib,ip,ipar),prms(ib,ip,ipar),mrms(ib,ip,ipar),rcpara%month) !in file homtests.f90 line 1278 


             call spagstat(plus(ib,ip,ipar,1),minus(ib,ip,ipar,1),tpcount(ib,ip,ipar),tmcount(ib,ip,ipar),prms(ib,ip,ipar),mrms(ib,ip,ipar),obssig(ib,ip,ipar,1),rcpara%miss_val) !in file rfcor.f90 line 26
             tausig(ib,ip,ipar,1)=obssig(ib,ip,ipar,1)
             if(dunit .gt. 0) write(dunit,'(A7,I4,I2,I6,I3,I2,I5,I5,3F9.3,F8.1)') cstatnr//': ',1,0,tbindex(ib),ip,ipar,tpcount(ib,ip,ipar),tmcount(ib,ip,ipar),plus(ib,ip,ipar,1),minus(ib,ip,ipar,1),obssig(ib,ip,ipar,1),0.

             if(ldebug) write(*,'(A6,I6.5,I6,3I3,2I5,2F8.3,2I5,2F8.3)') 'tcount ',wmonrs(index(1)),tbindex(ib),ib,ip,ipar,tpcount(ib,ip,ipar),tmcount(ib,ip,ipar),plus(ib,ip,ipar,1),minus(ib,ip,ipar,1),left_maxlen(ib,ipar,1),right_maxlen(ib,ipar,1),trasocorrs(tbindex(ib)+1,ip,ipar)-trasocorrs(tbindex(ib)-1,ip,ipar),obssig(ib,ip,ipar,1)


          enddo
       enddo
    enddo


    call read_IRI(iunit,'RS_IRI_CORRTABLES',rcpara,'MKIII-RS80',iplevs,profiles,err) !in file rfcorio.f90 line 71

    found=0
    full=0.
    stnumcrit=0.

    bimax=0
    do ib=1,tbi
       if(tbindex(ib) .lt. rcpara%nmax-3) bimax=bimax+1
    enddo

    if(bimax .gt. rcpara%brmax) then
       write(*,*),statnr,bimax,' too large'
    endif
    if(iter .eq. 1) then
!!$OMP CRITICAL
       !$ call omp_set_lock(omp_lp(rcpara%statmax+2))
       if(icache(index(1)) .gt. 0) tbicache(:,:,:,icache(index(1)))=0
       !$ call omp_unset_lock(omp_lp(rcpara%statmax+2))
!!$OMP END CRITICAL
    endif

!h    tbindexsave=tbindex ! needed later on for adjustment of most recent time series

    iana=0
    istat=1
    istatmax=0
    do ib=bimax,2,-1

  !     if(any(tbindex(ib)-rcpara%snht_maxlen/2 .eq. lasts)) then
  !        if(ib .ne. bimax) tbindex(ib)=tbindex(ib+1)
  !        cycle
  !     endif

       istat=0 ! include test station to get data into cache
       istati=0
       halt=0
       eistart=toindex(19790101,rcpara)
       open(11,form='formatted')

       do while(halt .eq. 0 .and. (any(found(ib,:) .lt. minst ))  .and. (any(found(ib,:) .lt. maxval((/1,minst/2/)) .or. dists(index(istat+1))*3.1415729/180.*6370. .lt. rcpara%weight_distance) &
            .or. abs(wmolats(index(istat+1))-wmolats(index(1))) .lt. 30.) .and. wmonrs(index(istat+1)) .ne. 0 .and. (istat .lt. wmostats-1) )

          !abs(wmolats(index(istat+1))-wmolats(index(1))) .lt. 30.   
          istat=istat+1
          if(istat.gt.istatmax) istatmax=istat


          !    if( wmonrs(index(1)) .gt. 91700 .and. wmonrs(index(1)) .lt. 92000 .and. wmonrs(index(istat)) .gt. 91700 .and. wmonrs(index(istat)) .lt. 92000 .and. istat .ne. 1 )  then
          !      cycle
          !    endif

          if( istat .eq. 1 .or.  &
               (.not. manual .or. wmonrs(index(istat)).lt. 89000 .or. wmonrs(index(istat)).gt.90000 .or. wmolats(index(1)).lt. -50.)  .and. &
               (wmonrs(index(istat)).lt. 42000 .or. wmonrs(index(istat)).gt. 44000)  .and. &
               (wmonrs(index(istat)).lt. 48500 .or. wmonrs(index(istat)).gt. 49000)  .and. &
               (wmonrs(index(istat)).lt. 47000 .or. wmonrs(index(istat)).gt. 47100)  .and. &
               (tbindex(ib) .gt. eistart .or. (wmonrs(index(1)) .lt. 20000 .or. wmonrs(index(1)).gt. 40000) .or.  (wmonrs(index(istat)).lt. 20000 .or. wmonrs(index(istat)).gt. 40000))  .and. &
               (tbindex(ib) .gt. eistart .or. (wmonrs(index(1)) .lt. 50000 .or. wmonrs(index(1)).gt. 60000) .or.  (wmonrs(index(istat)).lt. 50000 .or. wmonrs(index(istat)).gt. 60000))  .and. &
               ( (wmonrs(index(1)) .lt. 47000 .or. wmonrs(index(1)).gt. 48000) .or.  (wmonrs(index(istat)).lt. 50000 .or. wmonrs(index(istat)).gt. 60000))  .and. &
               (.not. manual .or. .not. any(wmonrs(index(istat)) .eq. (/89009,96996,72486/))) .and. &
               dists(index(istat))*3.1415729/180.*6370. .lt. 5*rcpara%weight_distance&
               ) then

             write(cstatnr,'(I6.6)') wmonrs(index(istat))

             if( ex(index(istat))  ) then

                iana(ib)=iana(ib)+1

                !$ call omp_set_lock(omp_lp(rcpara%statmax+2))
                logcache=icache(index(istat)) .eq. 0
                if(logcache) then
                   icache(rcpara%statmax+1)=icache(rcpara%statmax+1)+1
                   icache(index(istat))=icache(rcpara%statmax+1)

                   icistat=icache(index(istat))
                   if (icistat.eq.9) then
                      write(*,*) icistat
                   endif
                   !$ call omp_unset_lock(omp_lp(rcpara%statmax+2))


                   !$ call omp_set_lock(omp_lp(icistat))

                   filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbincorrsave'//cstatnr//'.nc'
                   !      call read_firstiter(rcpara,iunit,filename,bindex,bi,rasocorrs,rios) !in file rfcorio.f90 line 18
                   CALL read_sonde_corr_daily_nc(filename, rcpara,index(istat), rios,rasocorrs, bindex,bi)


                   filename=trim(rcpara%prefix)//cstatnr//'/x'//cstatnr
                   !    call read_sonde_oper(iunit,filename,rcpara,wmolats(index(istat)),wmolons(index(istat)),tm,tanm,tfgm,tbcm,solarangles,ominuse40,adjustlocal,rtype,err2) !in file rfcorio.f90 line 611

                   if(rcpara%innov .eq. 'RI' .or. rcpara%innov .eq. 'RE') then 
                      filename=trim(rcpara%prefix)//cstatnr//'/feedbackmerged'//cstatnr//'.nc'
                   else
                      filename=trim(rcpara%prefix)//cstatnr//'/feedbackmerged_lo'//cstatnr//'.nc'
                   endif
                   if(rcpara%innov .eq. 'RE') then 
                      CALL read_odb_nc(filename,rcpara,index(istat),err3,tm,tfgm,e20c0=tanm)
                      tfgm=tanm
                   else
                      ! tbcm required for MERRA
                      !MERRA      CALL read_odb_nc(filename,rcpara,istat,err3,tm,tfgm,tbcm=tbcm) !subroutine in read_txt_write_nc.f90
                      CALL read_odb_nc(filename,rcpara,index(istat),err3,tm,tfgm,tanm=tanm,bad_intervals=bad_intervals)
                      do ipar=1,rcpara%parmax
                         do ip=1,rcpara%pmax
                            do i=1,rcpara%nmax
                               if (isnan(tm(i,ip,ipar)) .or. isnan(tfgm(i,ip,ipar)) .or. abs(tfgm(i,ip,ipar))>20.) then
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
                      !MERRA      where(tfgm .ne. rcpara%miss_val .and. tbcm .ne. rcpara%miss_val)
                      !       tfgm=-(tfgm-0*tbcm)
                      !      endwhere
                   endif
                   !subroutine in read_txt_write_nc.f90
                   !    call eiminuse40(rcpara,wmolats(istat),wmolons(istat),ominuse40,adjust)


                   filename=trim(rcpara%prefix)//cstatnr//'/feedbackbgcorr'//cstatnr//'.nc'
                   !    call read_sonde_daily(iunit,filename,rcpara,err,eracorrs) !in file rfcorio.f90 line 1009
                   !     eracorrs=0.
                   CALL read_bgcorr_nc(filename,rcpara,index(istat),err3,eracorrs,(/'bg correctio',''/)) !subroutine in read_txt_write_nc.f90
                   if(rcpara%innov .ne. 'RI') eracorrs=0. 

                   if(lrem) then
                      filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbinremsave'//cstatnr//'.nc'
                      CALL read_sonde_corr_daily_nc(filename, rcpara,0, err,mrasocorrs, xbindex,xbi)
                      if(err .eq. 0) then
                         do ipar=1,rcpara%parmax
                            do ip=1,rcpara%pmax
                               do i=1,xbi
                                  xrasocorrs(xbindex(i)+1:xbindex(i+1),ip,ipar)=-mrasocorrs(i,ip,ipar)         
                               enddo
                               xrasocorrs(1,ip,ipar)=-mrasocorrs(1,ip,ipar) 
                               do i=1,rcpara%nmax
                                  if(tm(i,ip,ipar) .ne. rcpara%miss_val .and. tfgm(i,ip,ipar) .ne. rcpara%miss_val ) then
                                     tm(i,ip,ipar)=tm(i,ip,ipar)+xrasocorrs(i,ip,ipar)
                                     tfgm(i,ip,ipar)=tfgm(i,ip,ipar)+xrasocorrs(i,ip,ipar)
                                  endif
                               enddo
                            enddo
                         enddo
                      endif
                   endif !lrem

                   istati=icistat
                   l=0
                   do i=1,rcpara%nmax
                      if((tfgm(i,12,1) .ne. rcpara%miss_val .or. tfgm(i,12,2) .ne. rcpara%miss_val) .and.  (tm(i,12,1) .ne. rcpara%miss_val .or.  tm(i,12,2) .ne. rcpara%miss_val)) then
                         l=l+1            
                         cachetmp(l)=i
                      endif
                   enddo
                   tfgmcr(icistat)%vals=l
                   if(l .eq. 0) then
                      print*, 'vals is zero, cycle'
                      !         icache(rcpara%statmax+1)=icache(rcpara%statmax+1)-1
                      allocate(tfgmcr(icistat)%index(1))
                      tfgmcr(icistat)%index=1
                      allocate(tfgmcr(icistat)%feld(1,rcpara%pmax,rcpara%parmax))
                      allocate(tmcr(icistat)%feld(1,rcpara%pmax,rcpara%parmax))
                      !$ call omp_unset_lock(omp_lp(icistat))
                      cycle
                      !call abort
                   endif
                   allocate(tfgmcr(icistat)%index(l))
                   tfgmcr(icistat)%index=cachetmp(1:l)
                   allocate(tfgmcr(icistat)%feld(l,rcpara%pmax,rcpara%parmax))
                   allocate(tmcr(icistat)%feld(l,rcpara%pmax,rcpara%parmax))

                   do j=2,3
                      tccr(icistat,j)%vals=0
                      tccr(icistat,j)%index(1:bi)=bindex(1:bi)
                      !       allocate(tccr(icistat,j)%feld(rcpara%brmax,rcpara%pmax,rcpara%parmax))
                      tccr(icistat,j)%feld=rcpara%miss_val
                   enddo
                   tccr(icistat,2)%vals=bi ! other two are initalized to bi in second iteration!
                   !     tfgmcr(icistat)%feld=rcpara%miss_val
                   do ipar=1,rcpara%parmax
                      do ip=1,rcpara%pmax-2
                         do l=1,tfgmcr(icistat)%vals
                            i=tfgmcr(icistat)%index(l)
                            if(tfgm(i,ip,ipar) .ne. rcpara%miss_val .and. tm(i,ip,ipar) .ne. rcpara%miss_val) then
                               if(eracorrs(i,ip,ipar) .ne. rcpara%miss_val) then
                                  tfgmcr(icistat)%feld(l,ip,ipar)=tfgm(i,ip,ipar)-eracorrs(i,ip,ipar)!+rasocorrs(i,ip,ipar)
                               else
                                  !           tfgmcr(icistat)%feld(l,ip,ipar)=tfgm(i,ip,ipar)
                                  tfgmcr(icistat)%feld(l,ip,ipar)=rcpara%miss_val 
                               endif
                               tmcr(icistat)%feld(l,ip,ipar)=tm(i,ip,ipar)!-rasocorrs(i,ip,ipar)

                               !         tmpcr(icistat)%feld(l,ip,ipar)=tm(i,ip,ipar)
                            else
                               tfgmcr(icistat)%feld(l,ip,ipar)=rcpara%miss_val
                               tmcr(icistat)%feld(l,ip,ipar)=rcpara%miss_val

                               !         tmpcr(icistat)%feld(l,ip,ipar)=rcpara%miss_val
                            endif
                         enddo
                      enddo
                   enddo

                   !    if (any(isnan(tfgmcr(icistat)%feld))) then
                   !      call abort()
                   !    endif
                   llasts=0
                   lgcount=0
                   call detect_gaps(rcpara,wmonrs(index(istat)),tfgm,lmidx,llasts,lgcount) !in file rfcor.f90 line 113
                   lastscache(:,:,icistat)=llasts
                   gcountcache(:,icistat)=lgcount


                   write(*,*) wmonrs(index(istat)),' make written to cache with index',icache(index(istat))
                   !$ call omp_unset_lock(omp_lp(icistat))


                else
                   icistat=icache(index(istat))
                   !$ call omp_unset_lock(omp_lp(rcpara%statmax+2))

                endif ! icache

                if(istat .eq. 1) cycle

                !$ call omp_set_lock(omp_lp(icistat))
                bi=tccr(icistat,2)%vals
                bindex=rcpara%nmax
                bindex(1:bi)=tccr(icistat,2)%index(1:bi)
                !$ call omp_unset_lock(omp_lp(icistat))


                write(cstatnr,'(I6.6)') wmonrs(index(istat))
                filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbincorrsave'//cstatnr//'.nc'
                !      CALL read_sonde_corr_daily_nc(filename, rcpara,index(istat), rios,rasocorrs, bindex2,bi2)
                call modify_tbindex(rcpara,wmonrs,index,istat,bindex,bi,lastscache(:,:,icistat),.false.,iter,meta_s,alarm,manual,ldebug) !in file rfcomp_1.f90 line 321

                par: do ipar=1,rcpara%parmax

                   if (ib .gt. 1 .and. ib .lt. rcpara%mmax .and. found(ib,ipar) .lt. minst) then
                      if(any(tbindex(ib-1)-rcpara%snht_maxlen/2.eq. lasts)) then
                         left_lengtha=tbindex(ib)-tbindex(ib-2)
                      else
                         left_lengtha=tbindex(ib)-tbindex(ib-1)
                      endif
                      !        right_lengtha=tbindex(ib+1)-tbindex(ib)
                      right_lengtha=rcpara%nmax-2-tbindex(ib)

                      if(iter .le. 1) then
                         if(right_lengtha .gt. rcpara%mean_maxlen/2) right_lengtha=rcpara%mean_maxlen/2
                         if(left_lengtha .gt. rcpara%mean_maxlen/2) left_lengtha=rcpara%mean_maxlen/2
                      else
                         if(right_lengtha .gt. rcpara%mean_maxlen) right_lengtha=rcpara%mean_maxlen
                         if(left_lengtha .gt. rcpara%mean_maxlen) left_lengtha=rcpara%mean_maxlen
                      endif

                      right_maxlen(ib,ipar,istat)=right_lengtha
                      left_maxlen(ib,ipar,istat)=left_lengtha
                      !        if(iter .eq. 0) then
                      if(ib .eq. 12 .and. wmonrs(index(istat)) .eq. 47991) then
                         write(*,*) 'Test'
                      endif
                      do rib=1,bi
                         if(bindex(rib) .gt. tbindex(ib)-left_lengtha .and. bindex(rib) .le. tbindex(ib)) then
                            left_maxlen(ib,ipar,istat)=tbindex(ib)-bindex(rib)
                            if(ldebug) write(*,'(A10,2I7,6I6)') 'leftmaxlen',wmonrs(index(1)),wmonrs(index(istat)),ib,ipar,istat,tbindex(ib),bindex(rib),left_maxlen(ib,ipar,istat)
                         endif
                      enddo
                      if(istat<30 .and. left_maxlen(ib,ipar,istat)<0.5*left_lengtha) then 
                         cycle
                      endif
                      do rib=bi,1,-1
                         if(bindex(rib) .ge. tbindex(ib) .and. bindex(rib) .lt. tbindex(ib)+right_lengtha) then

                            right_maxlen(ib,ipar,istat)=bindex(rib)-tbindex(ib)
                            if(ldebug) write(*,'(A11,2I7,6I6)') 'rightmaxlen',wmonrs(index(1)),wmonrs(index(istat)),ib,ipar,istat,tbindex(ib),bindex(rib),left_maxlen(ib,ipar,istat)
                         endif
                      enddo
                      !         endif 

                      if(ib .eq. 12 .and. wmonrs(index(istat)) .eq. 47991) then
                         write(*,*) 'Test'
                      endif
                      mfak=2.0*rcpara%mean_maxlen/min_count/20
                      do i=1,rcpara%pmax-2
                         mfak(i)=mfak(i)*(1+(2*i-1.0)/rcpara%pmax)
                         if(mfak(i) .gt. 2*mfak(1)) mfak(i)=2*mfak(1)
                      enddo
                      ! at end of series there may be only few data after break - in this case reduce mfak
                      do i=tbindex(ib)+right_maxlen(ib,ipar,1),tbindex(ib),-1
                         if(ttfgm(i,8,ipar) .ne. rcpara%miss_val) exit
                      enddo
                      if(tbindex(ib)+1.5*mfak(8)*min_count .gt. i .and. tbindex(ib)+min_count .lt. i ) then
                         mfak(:)=1.
                      endif

                      ifail=0
                      if(rcpara%downwards) then
                         ipstart=1
                         ipstop=rcpara%pmax-2
                         inc=1
                      else
                         ipstart=rcpara%pmax-2
                         ipstop=1
                         inc=-1
                      endif

                      use_night=.false.
                      if(ipar .eq. 2 .and. (statnr .gt. 80000 .and. statnr .lt. 89000 .or. statnr .eq. 61901 .or. statnr .eq. 81405 )) use_night=.true.
                      if(ipar .eq. 1 .and. (statnr .gt. 94000 .and. statnr .lt. 95000 )) use_night=.true.

                      pp:      do ip=ipstart,ipstop,inc

                         if(obsspagarrcount(ib,ip,ipar) .lt. minst .and. plus(ib,ip,ipar,1) .ne. rcpara%miss_val .and. minus(ib,ip,ipar,1) .ne. rcpara%miss_val) then

                            if(.not. use_night) then

                               if(iter .eq. 2) then
                                  right_maxlen(ib,ipar,istat)=right_lengtha
                                  left_maxlen(ib,ipar,istat)=left_lengtha
                                  do rib=1,bi
                                     if(tbicache(rib,ip,ipar,icistat) .gt. tbindex(ib)-left_maxlen(ib,ipar,istat) .and. tbicache(rib,ip,ipar,icistat) .le. tbindex(ib)) then
                                        left_maxlen(ib,ipar,istat)=tbindex(ib)-tbicache(rib,ip,ipar,icistat)
                                        if(ldebug) write(*,'(2I6.5,3I3,I6,A14,I5)') wmonrs(index(1)),wmonrs(index(istat)),ib,rib, ip, tbicache(rib,ip,ipar,icistat), ' 1 leftmaxlen ', left_maxlen(ib,ipar,istat)    
                                     endif
                                  enddo
                                  do rib=bi,1,-1
                                     if(tbicache(rib,ip,ipar,icistat) .ge. tbindex(ib) .and. tbicache(rib,ip,ipar,icistat) .lt. tbindex(ib)+right_maxlen(ib,ipar,istat)) then

                                        right_maxlen(ib,ipar,istat)= tbicache(rib,ip,ipar,icistat)-tbindex(ib)
                                        if(ldebug) write(*,'(2I6.5,3I3,I6,A14,I5)') wmonrs(index(1)),wmonrs(index(istat)),ib,rib, ip, tbicache(rib,ip,ipar,icistat), ' 1 rightmaxlen ', right_maxlen(ib,ipar,istat)    
                                     endif
                                  enddo
                               endif

                               if(ldebug) write(*,'(A8,2I6.5,2I3,4I5)') 'maxlens ',wmonrs(index(1)),wmonrs(index(istat)),ib,ipar,left_maxlen(ib,ipar,istat),right_maxlen(ib,ipar,istat),left_lengtha,right_lengtha

                               good=.false.
                               if (left_maxlen(ib,ipar,istat) .gt. min_count .and. right_maxlen(ib,ipar,istat) .gt. mfak(ip)*min_count) then 

                                  call  ch2(cachehilf,tfgmcr(icistat)%vals,tfgmcr(icistat)%index,tfgmcr(icistat)%feld(:,ip,ipar),tccr(icistat,3)%index,&
                                       tccr(icistat,3)%feld(:,ip,ipar),rasocorrhomd(:,ip,ipar,2),rcpara,tbindex,ib,ip,ipar,istat,left_maxlen,right_maxlen,icistat,iter,-1.0_JPRM)

                                  call estimate(rcpara,ttfgm,cachehilf,left_maxlen,right_maxlen,tbindex,min_count,mfak(ip),dists,wmonrs,wmolons,wmolats,index,&
                                       tauspagarr,ltauspagarr,tausig,taucorrp,taucorrm,tauplus,tauminus,pcount,mcount,taupms,taumms,trasocorrs,tpcount,tmcount,ib,ip,ipar,istat,iter,2,1,ldebug) !in file rfcomp_1.f90 line 272
                                  if(dunit .gt. 0) write(dunit,'(I6.6,A2,I4,I2,I6,I3,I2,I5,I5,3F9.3,F8.1)') wmonrs(index(1)),': ',istat,2,tbindex(ib),ip,ipar,pcount(ib,ip,ipar),mcount(ib,ip,ipar),plus(ib,ip,ipar,istat),minus(ib,ip,ipar,istat),tausig(ib,ip,ipar,istat),dists(index(istat))*3.1415729/180.*6370

                                  good=goodest(rcpara,tausig,taucorrp,taucorrm,pcount,mcount,left_maxlen,right_maxlen,tpcount,tmcount,min_count,ib,ip,ipar,istat,wmonrs,index,2,iter,.false.)


                                  if (good) then

                                     call  ch2(cachehilf2,tfgmcr(icistat)%vals,tfgmcr(icistat)%index,tmcr(icistat)%feld(:,ip,ipar),tccr(icistat,2)%index,&
                                          tccr(icistat,2)%feld(:,ip,ipar),rasocorrhomd(:,ip,ipar,4),rcpara,tbindex,ib,ip,ipar,istat,left_maxlen,right_maxlen,icistat,iter,1.0_JPRM)

                                     call estimate(rcpara,ttm,cachehilf2,left_maxlen,right_maxlen,tbindex,min_count,mfak(ip),dists,wmonrs,wmolons,wmolats,index,obsspagarr,lobsspagarr,&
                                          obssig,obscorrp,obscorrm,obsplus,obsminus,pcount,mcount,obspms,obsmms,trasocorrs,tpcount,tmcount,ib,ip,ipar,istat,iter,4,1,ldebug) !in file rfcomp_1.f90 line 272


                                     if(dunit .gt. 0) write(dunit,'(I6.6,A2,I4,I2,I6,I3,I2,I5,I5,3F9.3,F8.1)') wmonrs(index(1)),': ',istat,4,tbindex(ib),ip,ipar,pcount(ib,ip,ipar),mcount(ib,ip,ipar),plus(ib,ip,ipar,istat),minus(ib,ip,ipar,istat),obssig(ib,ip,ipar,istat),dists(index(istat))*3.1415729/180.*6370.
                                     if (ib==5 .and. ip==5 .and. ipar==2) then
                                        print*,wmonrs(index(istat))
                                        write(ch,'(I4)') left_maxlen(ib,ipar,istat)+right_maxlen(ib,ipar,istat)+1
                                        write(11,'(5I7)') wmonrs(index(1)), wmonrs(index(istat)), tbindex(ib),left_maxlen(ib,ipar,istat),right_maxlen(ib,ipar,istat)
                                        write(11,'('//ch//'F8.2)') ttfgm( tbindex(ib)-left_maxlen(ib,ipar,istat):tbindex(ib)+right_maxlen(ib,ipar,istat),ip,ipar)
                                        write(11,'('//ch//'F8.2)') cachehilf( tbindex(ib)-left_maxlen(ib,ipar,istat):tbindex(ib)+right_maxlen(ib,ipar,istat))
                                        write(11,'('//ch//'F8.2)') ttm( tbindex(ib)-left_maxlen(ib,ipar,istat):tbindex(ib)+right_maxlen(ib,ipar,istat),ip,ipar)
                                        write(11,'('//ch//'F8.2)') cachehilf2( tbindex(ib)-left_maxlen(ib,ipar,istat):tbindex(ib)+right_maxlen(ib,ipar,istat))
                              
                                     endif
                                  endif

                               endif

                            endif ! use_night

                            !                  if(abs(amod(real(wmolons(index(istat))-wmolons(index(1))),360.)) .lt. 20. .and. (wmolons(index(1)) .gt. -90. .and. wmolons(index(1)) .lt. 90 .and. ipar .eq. 2 .or. &
                            !(wmolons(index(1)) .lt. -90. .or. wmolons(index(1)) .gt. 90) .and. ipar .eq. 1)) then
                            if(iter .eq. 2) then 
                               right_maxlen(ib,ipar,istat)=right_lengtha
                               left_maxlen(ib,ipar,istat)=left_lengtha
                               do rib=1,bi
                                  if(tbicache(rib,ip,3-ipar,icistat) .gt. tbindex(ib)-left_maxlen(ib,ipar,istat) .and. tbicache(rib,ip,3-ipar,icistat) .le. tbindex(ib)) then
                                     left_maxlen(ib,ipar,istat)=tbindex(ib)-tbicache(rib,ip,3-ipar,icistat)
                                     if(ldebug) write(*,'(2I6.5,2I3,I6,A14,I5)') wmonrs(index(1)),wmonrs(index(istat)),ib,rib, tbicache(rib,ip,3-ipar,icistat), ' 2 leftmaxlen ', left_maxlen(ib,ipar,istat)    
                                  endif
                               enddo
                               do rib=bi,1,-1
                                  if(tbicache(rib,ip,3-ipar,icistat) .ge. tbindex(ib) .and. tbicache(rib,ip,3-ipar,icistat) .lt. tbindex(ib)+right_maxlen(ib,ipar,istat)) then

                                     right_maxlen(ib,ipar,istat)= tbicache(rib,ip,3-ipar,icistat)-tbindex(ib)
                                     if(ldebug) write(*,'(2I6.5,2I3,I6,A14,I5)') wmonrs(index(1)),wmonrs(index(istat)),ib,rib, tbicache(rib,ip,3-ipar,icistat), ' 2 rightmaxlen ', right_maxlen(ib,ipar,istat)    
                                  endif
                               enddo
                            endif

                            if(ldebug) write(*,'(A9,2I6.5,2I3,4I5)') 'maxlens2 ',wmonrs(index(1)),wmonrs(index(istat)),ib,ipar,left_maxlen(ib,ipar,istat),right_maxlen(ib,ipar,istat),left_lengtha,right_lengtha

                            if (.false. .and. left_maxlen(ib,ipar,istat) .gt. min_count  .and. right_maxlen(ib,ipar,istat) .gt. mfak(ip)*min_count) then 

                               call  ch2(cachehilf,tfgmcr(icistat)%vals,tfgmcr(icistat)%index,tfgmcr(icistat)%feld(:,ip,3-ipar),tccr(icistat,3)%index,&
                                    tccr(icistat,3)%feld(:,ip,3-ipar),rasocorrhomd(:,ip,ipar,2),rcpara,tbindex,ib,ip,ipar,istat,left_maxlen,right_maxlen,icistat,iter,-1.0_JPRM)

                               call estimate(rcpara,ttfgm,cachehilf,left_maxlen,right_maxlen,tbindex,min_count,mfak(ip),dists,wmonrs,wmolons,wmolats,index,&
                                    tauspagarr2,ltauspagarr2,tausig2,taucorrp,taucorrm,tauplus,tauminus,pcount2,mcount2,taupms,taumms,trasocorrs,tpcount,tmcount,ib,ip,ipar,istat,iter,2,2,ldebug) !in file rfcomp_1.f90 line 272

                               if(dunit .gt. 0) write(dunit,'(I6.6,A2,I4,I2,I6,I3,I2,I5,I5,3F9.3,F8.1)') wmonrs(index(1)),': ',istat,2,tbindex(ib),ip,ipar,pcount(ib,ip,ipar),mcount(ib,ip,ipar),plus(ib,ip,ipar,istat),minus(ib,ip,ipar,istat),tausig2(ib,ip,ipar,istat),dists(index(istat))*3.1415729/180.*6370.

                               if (goodest(rcpara,tausig2,taucorrp,taucorrm,pcount2,mcount2,left_maxlen,right_maxlen,tpcount,tmcount,min_count,ib,ip,ipar,istat,wmonrs,index,2,iter,.false.)) then

                                  call  ch2(cachehilf,tfgmcr(icistat)%vals,tfgmcr(icistat)%index,tmcr(icistat)%feld(:,ip,3-ipar),tccr(icistat,2)%index,&
                                       tccr(icistat,2)%feld(:,ip,3-ipar),rasocorrhomd(:,ip,ipar,4),rcpara,tbindex,ib,ip,ipar,istat,left_maxlen,right_maxlen,icistat,iter,1.0_JPRM)

                                  call estimate(rcpara,ttm,cachehilf,left_maxlen,right_maxlen,tbindex,min_count,mfak(ip),dists,wmonrs,wmolons,wmolats,index,obsspagarr2,&
                                       lobsspagarr2,obssig2,obscorrp,obscorrm,obsplus,obsminus,pcount2,mcount2,obspms,obsmms,trasocorrs,tpcount,tmcount,ib,ip,ipar,istat,iter,4,2,ldebug) !in file rfcomp_1.f90 line 272
                                  if(dunit .gt. 0) write(dunit,'(I6.6,A2,I4,I2,I6,I3,I2,I5,I5,3F9.3,F8.1)') wmonrs(index(1)),': ',istat,4,tbindex(ib),ip,ipar,pcount(ib,ip,ipar),mcount(ib,ip,ipar),plus(ib,ip,ipar,istat),minus(ib,ip,ipar,istat),obssig2(ib,ip,ipar,istat),dists(index(istat))*3.1415729/180.*6370.

                               endif
                            endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            ! used -105 since RS takes an hour (15 degrees) to reach Stratosphere

                            lday=(statnr .lt. 70000 .or. statnr .gt. 75000) .and. wmolons(index(1)) .gt. -105. .and. wmolons(index(1)) .lt. 75 .and. ipar .eq. 2 .or. &
                                 (wmolons(index(1)) .lt. -105. .or. wmolons(index(1)) .gt. 75) .and. ipar .eq. 1

                            !            if((lday .or. use_night .or. &
                            !               0.8*minval((/pcount2(ib,ip,ipar),mcount2(ib,ip,ipar)/)) .gt. minval((/pcount(ib,ip,ipar),mcount(ib,ip,ipar)/)) ) &
                            if(lday .and. obsspagarr2(ib,ip,ipar,istat) .ne. rcpara%miss_val) then
                               obsspagarr(ib,ip,ipar,istat)=obsspagarr2(ib,ip,ipar,istat)
                               tauspagarr(ib,ip,ipar,istat)=tauspagarr2(ib,ip,ipar,istat)
                               obssig(ib,ip,ipar,istat)=obssig2(ib,ip,ipar,istat)
                               tausig(ib,ip,ipar,istat)=tausig2(ib,ip,ipar,istat)
                               lobsspagarr(ib,ip,ipar,istat)=lobsspagarr2(ib,ip,ipar,istat)
                               ltauspagarr(ib,ip,ipar,istat)=ltauspagarr2(ib,ip,ipar,istat)
                               use_night=.true.
                               if(ldebug) write(*,'(6I6,F8.2,a10)') wmonrs(index(istat)),ib,ip,ipar,minval((/pcount2(ib,ip,ipar),mcount2(ib,ip,ipar)/)),minval((/pcount(ib,ip,ipar),mcount(ib,ip,ipar)/)),obsspagarr2(ib,ip,ipar,istat),'12h shift'
                            endif

                            if(obsspagarr(ib,ip,ipar,istat) .eq. rcpara%miss_val)  then
                               ifail(ipar)=ifail(ipar)+1 
                            else
                               if(obsspagarr(ib,ip,ipar,istat) .ne. rcpara%miss_val .and. lobsspagarr(ib,ip,ipar,istat)) obsspagarrcount(ib,ip,ipar)=obsspagarrcount(ib,ip,ipar)+1
                            endif
                            if(.not. rcpara%downwards .and. ifail(ipar) .gt. 1) exit pp
                            if(rcpara%downwards .and. ifail(ipar) .gt. 1) exit pp



                         endif   !plus/minus
                         if (ip.eq.6.and. ib.eq.7.and. ipar.eq.2) then
                            !               write(*,'(A6,I6,2F8.3)')'ospag:',wmonrs(index(istat)),obsspagarr(ib,ip,ipar,istat),obsspagarr2(ib,ip,ipar,istat)
                         endif
                      enddo pp !ip   
                   endif

                enddo par ! ipar


                ip=1
                ifail=0
                imiss=0
                do ipar=1,rcpara%parmax
                   do ip=1,rcpara%pmax-4
                      if(plus(ib,ip,ipar,1) .ne. rcpara%miss_val .and. minus(ib,ip,ipar,1) .ne. rcpara%miss_val) then
                         if((obsspagarr(ib,ip,ipar,istat) .eq. rcpara%miss_val .or. .not. lobsspagarr(ib,ip,ipar,istat)) .and. obsspagarrcount(ib,ip,ipar) .lt. minst) then
                            ifail(ipar)=ifail(ipar)+1
                         endif
                      else 
                         imiss(ipar)=imiss(ipar)+1
                      endif
                   enddo
                   if(ifail(ipar) .lt. 2 .and. imiss(ipar) .lt. rcpara%pmax-4 .and. found(ib,ipar) .lt. minst) then
                      found(ib,ipar)=found(ib,ipar)+1
                      full(ib,ipar)=istat
                   endif
                   if(imiss(ipar) .ge. rcpara%pmax-4) then
                      found(ib,ipar)=2*minst
                   endif
                enddo

                if(ldebug) write(*,'(A8,I6.6,I3,I10,A5,A6,4I3,A7,2I4,A6,F9.1,2F5.2)') 'Station ',wmonrs(index(1)),ib,todate(tbindex(ib),rcpara),' ref ',cstatnr,imiss,ifail,' found ',found(ib,:),' dist', dists(index(istat))*3.1415729/180.*6370.

             else
                !      write(*,*) wmonrs(index(istat)),'feedbackmerged or feedbackglobbincorsave not found'

                if(wmonrs(index(istat)) .eq. 0) then
                   write(*,*) istat,'error'
                endif
             endif
             istati=icache(rcpara%statmax+1) 
             if (istati .eq. rcpara%cachemax-1) halt=1


          endif


       enddo !while

       close(11)

       tbindex(tbi+1)=0

       if(found(ib,1) .eq. 0 .or. found(ib,2) .eq. 0) then
          write(*,*) wmonrs(index(1)),' No buddy found'
       endif

       call rweight_mean(rcpara,wmonrs,obsspagarr,lobsspagarr,obsxdrarr,obsxdsarr,index,obssig,dists,ib,typ,minst,istatmax) !in file rfcomp_1.f90 line 165


       call rweight_mean(rcpara,wmonrs,tauspagarr,ltauspagarr,xdrarr,xdsarr,index,tausig,dists,ib,typ,minst,istatmax) !in file rfcomp_1.f90 line 165

       do ipar=1,rcpara%parmax
          do ip=1,rcpara%pmax-2
             if(obsxdsarr(ib,ip,ipar) .ne. rcpara%miss_val) then
                r1=trasocorrs(tbindex(ib)+1,ip,ipar)
                r2=trasocorrs(tbindex(ib)-1,ip,ipar)
                if(ldebug) write(*,'(A5,I6.6,A2,I4,I2,I3,I6,I3,I2,3F9.3,F8.1)') 'obs: ',wmonrs(index(1)),'::',0,4,ib,tbindex(ib),ip,ipar,obsxdrarr(ib,ip,ipar)+(r2-r1),obsxdsarr(ib,ip,ipar),dists(index(istat))*3.1415729/180.*6370.
                if(dunit .gt. 0) write(dunit,'(I6.6,A2,I4,I2,I6,I3,I2,3F9.3,F8.1)') wmonrs(index(1)),'::',0,4,tbindex(ib),ip,ipar,obsxdrarr(ib,ip,ipar)+(r2-r1),obsxdsarr(ib,ip,ipar),dists(index(istat))*3.1415729/180.*6370.
                if(dunit .gt. 0) write(dunit,'(I6.6,A2,I4,I2,I6,I3,I2,3F9.3,F8.1)') wmonrs(index(1)),'::',0,2,tbindex(ib),ip,ipar,xdrarr(ib,ip,ipar)+(r2-r1),xdsarr(ib,ip,ipar),dists(index(istat))*3.1415729/180.*6370.
                if(dunit .gt. 0) write(dunit,'(I6.6,A2,I4,I2,I6,I3,I2,3F9.3,F8.1)') wmonrs(index(1)),'::',0,0,tbindex(ib),ip,ipar,obsxdrarr(ib,ip,ipar)-xdrarr(ib,ip,ipar),rcpara%miss_val,dists(index(istat))*3.1415729/180.*6370.
             endif
          enddo
       enddo


       ! Derzeit nur kurze Intervalle
       ! 1 ist normales mittel tau nm
       ! 2 ist normales mittel obs nmra
       ! 3 ist radiusgewichtetes mittel tau rgm
       ! 4 ist radiusgewichtetes mittel obs rgm
       isave=0
       do i=1,tbisave
         if(tbindexsave(i)==tbindex(ib)) then
            isave=i
            exit
         endif
       enddo
       do isample=2,4,2
          do ipar=1,rcpara%parmax
             do ip=1,rcpara%pmax


                ipp=ip
                !     if(ip .eq. 14) ipp=13
                select case(isample)
                case(1)
                   hilf=darr(ib,ipp,ipar)
                case(2)
                   if(xdrarr(ib,ip,ipar) .eq. rcpara%miss_val .and. rcpara%rimax_miss < rcpara%max_miss ) then
                     hilf=trasocorrs(tbindex(ib)+1,ip,ipar)-trasocorrs(tbindex(ib)-1,ip,ipar)
!                     hilf=-(plus(ib,ip,ipar,1)-minus(ib,ipp,ipar,1))
                   else
                     hilf=xdrarr(ib,ipp,ipar)
                   endif
                case(3)
                   hilf=obsdarr(ib,ipp,ipar)
                case(4)
                   if(obsxdrarr(ib,ip,ipar) .eq. rcpara%miss_val  .and. rcpara%rimax_miss<rcpara%max_miss ) then
                     hilf=trasocorrs(tbindex(ib)+1,ip,ipar)-trasocorrs(tbindex(ib)-1,ip,ipar)
                   else
                     if(ipp .ge. 13) then
                       hilf=xdrarr(ib,ipp,ipar)
                     else
                       hilf=obsxdrarr(ib,ipp,ipar)
                     endif
                   endif
                case(5)
                case(6)
                case(7)
                case(8)
                end select
                !
                ! Vorsicht bei 700/850 hPa - Skalieren falls RAOBCORE-Korrektur viel kleiner ist
                if(hilf .ne. rcpara%miss_val .and. abs(hilf) .gt. 20) then
                   write(*,*) 'spurious hilf'
                   write(*,'(A8,I6.6,I7,3I3,2I5,2F8.3)') 'Station ',wmonrs(index(1)),breakdates(ib)/100,iter,ip,ipar,tpcount(ib,ip,ipar),tmcount(ib,ip,ipar),plus(ib,ip,ipar,1)-minus(ib,ip,ipar,1),hilf
                   hilf=rcpara%miss_val
                endif
                if(ip .gt. 12 .and. hilf .ne. rcpara%miss_val) then 
                   if(abs(plus(ib,ip,ipar,1)-minus(ib,ip,ipar,1)) .gt. 0) then
                      ratio=abs(hilf)/abs(plus(ib,ip,ipar,1)-minus(ib,ip,ipar,1))
                      if(abs(hilf) .gt. 0.1 .and. ratio .gt. 3) then
                         !          hilf=hilf/ratio*2.
                         write(*,*) ib,ip,ipar,ratio,' hilf scaled'
                      endif
                   endif
                endif

                ! South American fix for surface temperatures
                if(rcpara%smooth_method .eq. 1 .and. any(wmonrs(index(1)) .eq. (/85543,68263,68242,68424,68442,91610/)) .and. ip .ge. 12 ) then
                   if(tbindex(ib) .gt. toindex(19790101,rcpara) .and. tbindex(ib) .lt. toindex(20000101,rcpara)) then
                      if(hilf .ne. rcpara%miss_val) then
                         hilf=hilf*(14-ip)/3
                      endif
                   endif
                endif
                ! 
                ! zu starke korrektur australischer Stationen
                if(rcpara%smooth_method .eq. 1 .and.( (wmonrs(index(1)) .gt. 94000 .and. wmonrs(index(1)) .lt. 95000 .or. wmonrs(index(1)) .gt. 48600 .and. &
                     wmonrs(index(1)) .lt. 48600 .or. any(wmonrs(index(1)) .eq. (/91643,91517,91680,89611,89642,89564,89571/))))) then
                   !ip .gt. 5 .and. 
                   if(tbindex(ib) .gt. toindex(19860101,rcpara) .and. tbindex(ib) .lt. toindex(19891201,rcpara) .and. tbindex(ib-1) .lt. toindex(19870101,rcpara)) then

                      if(profiles(ip,ipar)*10 .ne. rcpara%miss_val) then
                         if(hilf .gt. profiles(ip,1)) hilf=profiles(ip,1)
                         darr(ib,ip,ipar)=hilf
                         xdrarr(ib,ip,ipar)=hilf
                         obsdarr(ib,ip,ipar)=hilf
                         obsxdrarr(ib,ip,ipar)=hilf
                      endif
                   endif
                endif

                if (hilf .ne. rcpara%miss_val) then
                   rasocorrhomd(1:tbindex(ib),ip,ipar,isample)=-hilf+rasocorrhomd(tbindex(ib)+1,ip,ipar,isample)
                else
                   rasocorrhomd(1:tbindex(ib),ip,ipar,isample)=rasocorrhomd(tbindex(ib)+1,ip,ipar,isample)
                endif
             enddo
          enddo
       enddo


    enddo !ib 

    tbindex=tbindexsave

    write(cform,'(I2.2)') bimax-1
    isubs=0
    iorig=0
    do ib=2,bimax
       ifail=0
       breakdates(ib)=todate(tbindex(ib),rcpara)
       do ipar=1,rcpara%parmax
          if(full(ib,ipar) .eq. 0) full(ib,ipar)=istat
          !if(iter .eq. 1) then
          do ip=1,rcpara%pmax-2
             !    if((plus(ib,ip,ipar,1) .ne. rcpara%miss_val .or. minus(ib,ip,ipar,1) .ne. rcpara%miss_val) .and. obsxdrarr(ib,ip,ipar) .eq. rcpara%miss_val) then
             if(obsxdrarr(ib,ip,ipar) .eq. rcpara%miss_val) then
                if (plus(ib,ip,ipar,1).ne. rcpara%miss_val .and. minus(ib,ip,ipar,1) .ne. rcpara%miss_val) then
                   if (rcpara%rimax_miss==rcpara%max_miss) then
                      write(*,'(A8,I6.6,I7,3I3,2I5,F8.3,A25)') 'Station ',wmonrs(index(1)),breakdates(ib)/100,iter,ip,ipar,tpcount(ib,ip,ipar),tmcount(ib,ip,ipar),plus(ib,ip,ipar,1)-minus(ib,ip,ipar,1),' level adjustment failed'
                      if(icache(index(1)) .gt. 0) then
                         !$ call omp_set_lock(omp_lp(icache(index(1))))
!!$omp critical
                         tbicache(ib,ip,ipar,icache(index(1)))=tbindex(ib)
!!$omp end critical
                         !$ call omp_unset_lock(omp_lp(icache(index(1))))
                      else
                         write(*,*) 'Warning: tbicache not filled',cstatnr
                      endif
                   else
                      obsxdrarr(ib,ip,ipar)=plus(ib,ip,ipar,1)-minus(ib,ip,ipar,1)
                      xdrarr(ib,ip,ipar)=plus(ib,ip,ipar,1)-minus(ib,ip,ipar,1)
                      write(*,'(A8,I6.6,I7,3I3,2I5,F8.3,A29)') 'Station ',wmonrs(index(1)),breakdates(ib)/100,iter,ip,ipar,tpcount(ib,ip,ipar),tmcount(ib,ip,ipar),plus(ib,ip,ipar,1)-minus(ib,ip,ipar,1),' level adjustment substituted'
                      ifail(ipar)=ifail(ipar)-1
                      isubs=isubs+1
                   endif
                endif
                ifail(ipar)=ifail(ipar)+1
             else
               iorig=iorig+1
             endif
!!$             if(icache(index(1)) .gt. 0) then
!!$                !$ call omp_set_lock(omp_lp(icache(index(1))))
! !$omp critical
!!$                tbicache(ib,ip,ipar,icache(index(1)))=tbindex(ib)
! !$omp end critical
!!$                !$ call omp_unset_lock(omp_lp(icache(index(1))))
!!$             else
!!$                write(*,*) 'Warning: tbicache not filled',cstatnr
!!$             endif
          enddo
          !  endif
       enddo
    enddo

    if(bimax .gt. 1) then
       write(*,'(A8,I6.6,2I3,A16,'//cform//'I7)') 'Station ',wmonrs(index(1)),tbi,iter,' final breaks   ',breakdates(2:bimax)/100
       write(*,'(A8,I6.6,2I3,A16,'//cform//'I7)') 'Station ',wmonrs(index(1)),tbi,iter,' final analysed ',iana(2:bimax)
       write(*,'(A8,I6.6,2I3,A16,'//cform//'I7)') 'Station ',wmonrs(index(1)),tbi,iter,' final found  00',found(2:bimax,1)
       write(*,'(A8,I6.6,2I3,A16,'//cform//'F7.0)') 'Station ',wmonrs(index(1)),tbi,iter,' final dists  00', dists(index(full(2:bimax,1)))*3.1415729/180.*6370. 
       write(*,'(A8,I6.6,2I3,A16,'//cform//'I7)') 'Station ',wmonrs(index(1)),tbi,iter,' final found  12',found(2:bimax,2)
       write(*,'(A8,I6.6,2I3,A16,'//cform//'F7.0)') 'Station ',wmonrs(index(1)),tbi,iter,' final dists  12', dists(index(full(2:bimax,2)))*3.1415729/180.*6370. 
       write(*,'(A8,I6.6,2I3,A16,2I4)') 'Station ',wmonrs(index(1)),tbi,iter,' final orig subs', iorig,isubs 
    endif


    igmax=maxval(gcount)
    xdrarr2=rcpara%miss_val
    !if(any(gcount .ne. 0)) then
    if(.false. .and. iter .eq. rcpara%maxiter .and. any(gcount .ne. 0)) then
       !! correct most recent part of series
       !! do this after breakpoint correction since then a longer interval for more accurate
       !! estimation can be used.

       found=0
       tauspagarr2=rcpara%miss_val
       ltauspagarr=.true.
       istatmax=0
       do ipar=1,rcpara%parmax
          do ig=1,gcount(ipar)

             if(lasts(ig,ipar) .ne. 1 .and. lasts(ig,ipar)+rcpara%snht_maxlen/2-1 .le. rcpara%nmax) then
                do ip=1,rcpara%pmax-2
                   rasocorrhomd(1:lasts(ig,ipar)+rcpara%snht_maxlen/2-1,ip,ipar,2)=rasocorrhomd(1:lasts(ig,ipar)+rcpara%snht_maxlen/2-1,ip,ipar,2)-rasocorrhomd(lasts(ig,ipar)+rcpara%snht_maxlen/2-1,ip,ipar,2)
                enddo
             else
                if(lasts(ig,ipar)+rcpara%snht_maxlen/2-1 .gt. rcpara%nmax) then 
                   print*, statnr,'lasts invalid',lasts(1:5,:)
                endif
             endif

             istat=1 
             do while ((found(ig,ipar) .lt. minst )  .and. dists(index(istat+1))*3.1415729/180.*6370. .lt. 2*rcpara%weight_distance .and. wmonrs(index(istat+1)) .ne. 0 .and. istat .lt. rcpara%cachemax )

                istat=istat+1
                if(istatmax .lt. istat) istatmax=istat

                !$ call omp_set_lock(omp_lp(rcpara%statmax+2))
                icistat=icache(index(istat))
                !$ call omp_unset_lock(omp_lp(rcpara%statmax+2))
                if((.not. ini_correct(ipar,index(istat))) .and. icistat .gt. 0) then

                   call load_richcorr(tccr,icistat,0,wmonrs,index(istat),rcpara,mrasocorrs)
!!$!$ call omp_set_lock(omp_lp(icistat))
!!$    logcache=tccr(icistat,1)%vals .ne. 0
!!$!$ call omp_unset_lock(omp_lp(icistat))
!!$      if(logcache) then
!!$
!!$!$ call omp_set_lock(omp_lp(icistat))
!!$      bi=tccr(icistat,1)%vals
!!$      bindex(1:bi)=tccr(icistat,1)%index(1:bi)
!!$      mrasocorrs(1:bi,:,:)=tccr(icistat,1)%feld(1:bi,:,:)
!!$
!!$!$ call omp_unset_lock(omp_lp(icistat))
!!$
!!$     else
!!$
!!$      write(filename,'(I6.6,a,a,a,I6.6,a)') wmonrs(index(istat)),'/feedbackglobbincorrsave_rit',rcpara%ens,'_',wmonrs(index(istat)),'.nc'
!!$      CALL read_sonde_corr_daily_nc(filename, rcpara,index(istat), err,mrasocorrs, bindex,bi)
!!$      if(err .ne. 0 .or. any(isnan(mrasocorrs))) then
!!$      write(filename,'(I6.6,a,I6.6,a)') wmonrs(index(istat)),'/feedbackglobbincorrsave',wmonrs(index(istat)),'.nc'
!!$        CALL read_sonde_corr_daily_nc(filename, rcpara,index(istat), err,mrasocorrs, bindex,bi) 
!!$        if(err .eq. 0) then 
!!$          write(*,*) 'read ',filename,' for initial adjustment since RICH-adjusted version was not available or spurious'
!!$        else
!!$          stop 'could not find data for initial adjustment'
!!$        endif
!!$      endif
!!$
!!$      if(bi .ge. 2) bindex(2:bi)=bindex(2:bi)-1
!!$
!!$!$ call omp_set_lock(omp_lp(icistat))
!!$      tccr(icistat,1)%vals=bi
!!$      tccr(icistat,1)%index=rcpara%nmax
!!$      tccr(icistat,1)%index(1:bi)=bindex(1:bi)
!!$      if(.not. allocated(tccr(icistat,1)%feld)) allocate(tccr(icistat,1)%feld(rcpara%brmax,rcpara%pmax,rcpara%parmax))
!!$      tccr(icistat,1)%feld(1:bi,:,:)=mrasocorrs(1:bi,:,:)
!!$!      if(any(isnan(tccr(icistat,1)%feld))) then
!!$!        stop
!!$!      endif
!!$!$ call omp_unset_lock(omp_lp(icistat))
!!$     endif


                   right_maxlen(1,ipar,istat)=rcpara%snht_maxlen
                   left_maxlen(1,ipar,istat)=lasts(ig,ipar)-midx(ig,ipar)
                   do ip=1,rcpara%pmax-2

                      cachehilf=rcpara%miss_val
                      call  ch2(cachehilf,tfgmcr(icistat)%vals,tfgmcr(icistat)%index,tfgmcr(icistat)%feld(:,ip,ipar),tccr(icistat,1)%index,&
                           tccr(icistat,1)%feld(:,ip,ipar),rasocorrhomd(:,ip,ipar,2),rcpara,lasts(ig,ipar),1,ip,ipar,istat,left_maxlen,right_maxlen,icistat,iter,-1.0_JPRM)

                      diff=0.
                      sqdiff=0.
                      ic=0
                      do i=lasts(ig,ipar)-left_maxlen(1,ipar,istat),lasts(ig,ipar)
                         if(cachehilf(i) .ne. rcpara%miss_val .and. ttfgm(i,ip,ipar) .ne. rcpara%miss_val) then
                            hilf=ttfgm(i,ip,ipar)-cachehilf(i)
                            diff=diff+hilf
                            sqdiff=sqdiff+hilf*hilf
                            ic(ipar)=ic(ipar)+1
                            !              if(isnan(diff+sqdiff)) then
                            !                  call abort()
                            !              endif
                         endif
                      enddo
                      if(ic(ipar) .gt. rcpara%snht_maxlen/8) then  !1/2 year
                         !              if(isnan(diff/ic(ipar))) then
                         !                  call abort()
                         !              endif
                         tauspagarr(tbi+ig,ip,ipar,istat)=diff/ic(ipar)
                         sqdiff=sqdiff/ic(ipar)-diff*diff/ic(ipar)/ic(ipar)
                         if(sqdiff .lt. 0) then
                            tausig(tbi+ig,ip,ipar,istat)=0.
                         else
                            tausig(tbi+ig,ip,ipar,istat)=sqrt(sqdiff)
                         endif
                      endif
                   enddo

                   if(count(tauspagarr(tbi+ig,:,ipar,istat) .ne. rcpara%miss_val) .gt. 5) found(ig,ipar)=found(ig,ipar)+1
                endif

             enddo

          enddo
       enddo

       do ig=1,igmax
          !  call rweight_mean(rcpara,wmonrs,tauspagarr2,ltauspagarr,xdrarr2,xdsarr,index,tausig,dists,ig,typ,minst,istatmax) !in file rfcomp_1.f90 line 165
          call rweight_mean(rcpara,wmonrs,tauspagarr,ltauspagarr,xdrarr,xdsarr,index,tausig,dists,tbi+ig,typ,minst,istatmax) !in file rfcomp_1.f90 line 165
       enddo

    endif


    tbindex2=tbindex
    tbi2=tbi+igmax
    tbindex2(tbi2)=rcpara%nmax
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (.false.) then
!!!!!!!!!!!!!!!
       ! assemble adjustment time series
       do j=2,4,2
          do ipar=1,rcpara%parmax
             do ip=1,rcpara%pmax
                do i=1,gcount(ipar)
                   if(i .eq. 1 .and. lasts(i,ipar) .gt. rcpara%old) then    
                      tbindex2(tbi+i)=rcpara%nmax-2
                      !      if(xdrarr(tbi+i,ip,ipar) .ne. rcpara%miss_val) then
                      !        if(ini_correct(ipar,index(1))) rasocorrhomd(1:rcpara%nmax-2,ip,ipar,j)=rasocorrhomd(1:rcpara%nmax-2,ip,ipar,j)+xdrarr(tbi+i,ip,ipar)
                      !      else
                      if(ini_correct(ipar,index(1))) rasocorrhomd(1:rcpara%nmax-2,ip,ipar,j)=rasocorrhomd(1:rcpara%nmax-2,ip,ipar,j)+trasocorrs(rcpara%nmax-2,ip,ipar)
                      !      endif
                   else
                      ! Correction also of most recent time series
                      tbindex2(tbi+i)=lasts(i,ipar)+rcpara%snht_maxlen/2-1
                      if(xdrarr(tbi+i,ip,ipar) .ne. rcpara%miss_val) then
                         rasocorrhomd(1:lasts(i,ipar)+rcpara%snht_maxlen/2-1,ip,ipar,j)=rasocorrhomd(1:lasts(i,ipar)+rcpara%snht_maxlen/2-1,ip,ipar,j)-rasocorrhomd(lasts(i,ipar)+rcpara%snht_maxlen/2-1,ip,ipar,j)+xdrarr(tbi+i,ip,ipar)
                      else
                         rasocorrhomd(1:lasts(i,ipar)+rcpara%snht_maxlen/2-1,ip,ipar,j)=rasocorrhomd(1:lasts(i,ipar)+rcpara%snht_maxlen/2-1,ip,ipar,j)-rasocorrhomd(lasts(i,ipar)+rcpara%snht_maxlen/2-1,ip,ipar,j)+trasocorrs(lasts(i,ipar)+rcpara%snht_maxlen/2,ip,ipar)
                      endif
                      rasocorrhomd(lasts(i,ipar)+rcpara%snht_maxlen/2-1,ip,ipar,j)=rasocorrhomd(lasts(i,ipar)+rcpara%snht_maxlen/2,ip,ipar,j)  
                   endif
                enddo
             enddo
          enddo
       enddo

       !ccccccccccccccccccccccccc
       xdrarr2=rcpara%miss_val
       !if(any(gcount .ne. 0)) then
       if(iter .eq. rcpara%maxiter .and. tbi .ne. 0) then
          !! correct most recent part of series
          !! do this after breakpoint correction since then a longer interval for more accurate
          !! estimation can be used.

          found=0
          tauspagarr2=rcpara%miss_val
          ltauspagarr=.true.
          sidx=0
          do ipar=1,rcpara%parmax
             do ib=tbi,2,-1
                if(tbindex(ib) .lt. rcpara%nmax-2) then
                   do ip=1,7
                      if(.not. any(sidx(ib:tbi,ip,ipar) .gt. 0 .and. sidx(ib:tbi,ip,ipar) .lt. tbindex(ib))) then
                         if(rasocorrhomd(tbindex(ib),ip,ipar,4) .eq. rasocorrhomd(tbindex(ib)+1,ip,ipar,4)) then
                            if(any(rasocorrhomd(tbindex(ib),8:14,ipar,4)-rasocorrhomd(tbindex(ib)+1,8:14,ipar,4) .ne. 0.)) then
                               ig=1
                               do while(lasts(ig,ipar).gt. tbindex(ib))
                                  ig=ig+1
                               enddo
                               pcount(ib,ip,ipar)=count(ttfgm(lasts(ig,ipar):tbindex(ib),ip,ipar) .ne. rcpara%miss_val)
                               if(pcount(ib,ip,ipar) .gt. rcpara%snht_maxlen/2) then
                                  ic(ipar)=0
                                  il:        do i=tbindex(ib),lasts(ig,ipar),-1
                                     if(ttfgm(i,ip,ipar) .ne. rcpara%miss_val) ic(ipar)=ic(ipar)+1
                                     if(ic(ipar) .gt. rcpara%snht_maxlen/2) then
                                        sidx(ib,ip,ipar)=i
                                        exit il
                                     endif
                                  enddo il
                                  write(*,'(a,8I6)') 'HL-gaps ',wmonrs(index(1)),ib,ip,ipar,lasts(ig,ipar),tbindex(ib),sidx(ib,ip,ipar),pcount(ib,ip,ipar)
                               endif
                            endif
                         endif
                      endif
                   enddo
                endif
             enddo
          enddo


          do ipar=1,rcpara%parmax
             do ib=2,tbi
                pfound(ib,:,ipar)=0
                istat=1 
                if(any(sidx(ib,1:7,ipar) .ne. 0)) then
                   stat: do while (any(pfound(ib,1:7,ipar) .lt. minst)  .and. dists(index(istat+1))*3.1415729/180.*6370. .lt. 2*rcpara%weight_distance .and. wmonrs(index(istat+1)) .ne. 0 .and. istat .lt. rcpara%cachemax )


                      istat=istat+1
                      if(istatmax .lt. istat) istatmax=istat


                      !$ call omp_set_lock(omp_lp(rcpara%statmax+2))
                      icistat=icache(index(istat))
                      !$ call omp_unset_lock(omp_lp(rcpara%statmax+2))
                      if((.not. ini_correct(ipar,index(istat))) .and. icistat .gt. 0) then

                         call load_richcorr(tccr,icistat,ipar,wmonrs,index(istat),rcpara,mrasocorrs)

                         do ip=1,7
                            if(sidx(ib,ip,ipar) .ne. 0 .and. pfound(ib,ip,ipar) .lt. minst) then
                               right_maxlen(ib,ipar,istat)=rcpara%snht_maxlen
                               left_maxlen(ib,ipar,istat)=tbindex(ib)-sidx(ib,ip,ipar)

                               cachehilf=rcpara%miss_val
                               call  ch2(cachehilf,tfgmcr(icistat)%vals,tfgmcr(icistat)%index,tfgmcr(icistat)%feld(:,ip,ipar),tccr(icistat,1)%index,&
                                    tccr(icistat,1)%feld(:,ip,ipar),rasocorrhomd(:,ip,ipar,2)-rasocorrhomd(tbindex(ib),ip,ipar,2),rcpara,tbindex,ib,ip,ipar,istat,left_maxlen,right_maxlen,icistat,iter,-1.0_JPRM)

                               diff=0.
                               sqdiff=0.
                               ic=0
                               do i=sidx(ib,ip,ipar),tbindex(ib)
                                  if(i.lt. 1) then
                                     write(*,*) 'error'
                                  endif
                                  if(cachehilf(i) .ne. rcpara%miss_val .and. ttfgm(i,ip,ipar) .ne. rcpara%miss_val) then
                                     hilf=ttfgm(i,ip,ipar)-cachehilf(i)
                                     diff=diff+hilf
                                     sqdiff=sqdiff+hilf*hilf
                                     ic(ipar)=ic(ipar)+1
                                     !              if(isnan(diff+sqdiff)) then
                                     !                  stop
                                     !              endif
                                  endif
                               enddo
                               ! write(*,'(a,5I6)') 'IC: ',wmonrs(index(1)),ib,ip,ipar,ic(ipar) 
                               if(ic(ipar) .gt. rcpara%snht_maxlen/8) then  !1/2 year
                                  !              if(isnan(diff/ic(ipar))) then
                                  !                  stop
                                  !              endif
                                  tauspagarr2(ib,ip,ipar,istat)=diff/ic(ipar)
                                  sqdiff=sqdiff/ic(ipar)-diff*diff/ic(ipar)/ic(ipar)
                                  if(sqdiff .lt. 0) then
                                     tausig(ib,ip,ipar,istat)=0.
                                  else
                                     tausig(ib,ip,ipar,istat)=sqrt(sqdiff)
                                  endif
                               endif

                               if(tauspagarr2(ib,ip,ipar,istat) .ne. rcpara%miss_val) pfound(ib,ip,ipar)=pfound(ib,ip,ipar)+1
                            endif
                         enddo ! ip

                      endif ! ini_correct
                   enddo stat ! istat
                endif ! sidx
             enddo ! ib
          enddo ! ipar


          ! assemble adjustment time series
          do ib=2,tbi
             if(any(sidx(ib,:,:) .ne. 0)) then
                call rweight_mean(rcpara,wmonrs,tauspagarr2,ltauspagarr,xdrarr2,xdsarr,index,tausig,dists,ib,typ,minst,istatmax) !in file rfcomp_1.f90 line 165
             endif
          enddo

          !endif

          do j=2,4,2
             do ipar=1,rcpara%parmax
                do ip=1,7
                   do ib=2,tbi
                      ! Correction also of most recent time series
                      if(xdrarr2(ib,ip,ipar) .ne. rcpara%miss_val) then
                         rasocorrhomd(1:tbindex(ib),ip,ipar,j)=rasocorrhomd(1:tbindex(ib),ip,ipar,j)-rasocorrhomd(tbindex(ib),ip,ipar,j)+xdrarr2(ib,ip,ipar)
                      endif
                   enddo
                enddo
             enddo
          enddo

       endif !iter

    endif !.false.

    !ccccccccccccccccccccccccc

    rasocorrhomd(rcpara%nmax,:,:,:)=0.

    write(cstatnr,'(I6.6)') statnr
    if (tbindex(tbi) .eq. rcpara%nmax-2) then
       tauspagarr(tbi:tbi2,:,:,:)=tauspagarr(tbi+1:tbi2+1,:,:,:)
       xdrarr(tbi:tbi2,:,:)=xdrarr(tbi+1:tbi2+1,:,:)
       obsxdrarr(tbi:tbi2,:,:)=xdrarr(tbi+1:tbi2+1,:,:)
    endif


    filename=trim(rcpara%prefix)//cstatnr//'/spagarr_'//citer//'_'//rcpara%ens//'_'//cstatnr
    call savespagarr(rcpara,iunit,filename,tauspagarr,tbindex2,tbi2,xdrarr,wmonrs,dists,index,istatmax) !in file rfcorio.f90 line 121

    filename=trim(rcpara%prefix)//cstatnr//'/obsspagarr_'//citer//'_'//rcpara%ens//'_'//cstatnr
    call savespagarr(rcpara,iunit,filename,obsspagarr,tbindex2,tbi2,obsxdrarr,wmonrs,dists,index,istatmax) !in file rfcorio.f90 line 121


    if(icache(index(1)) .gt. 0) then
       !$ call omp_set_lock(omp_lp(icache(index(1))))
       tccr(icache(index(1)),2)%vals=tbi
       tccr(icache(index(1)),3)%vals=tbi
       tccr(icache(index(1)),2)%index(1:tbi)=tbindex
       tccr(icache(index(1)),3)%index(1:tbi)=tbindex
       do bi=1,tbi
          if(iter .eq. 2) then
             !      tccr(icache(index(1)),1)%feld(bi,:,:)=rasocorrhomd(tbindex(bi),:,:,2)
          else
             tccr(icache(index(1)),1)%vals=0
             tccr(icache(index(1)),2)%feld(bi,:,:)=rasocorrhomd(tbindex(bi),:,:,4)
             tccr(icache(index(1)),3)%feld(bi,:,:)=rasocorrhomd(tbindex(bi),:,:,2)
          endif
       enddo
       tccr(icache(index(1)),2)%feld(tbi+1,:,:)=0
       tccr(icache(index(1)),3)%feld(tbi+1,:,:)=0

       !$ call omp_unset_lock(omp_lp(icache(index(1))))
    else
       write(*,*) 'error ',wmonrs(index(1))
    endif


    if(dunit .gt. 0) close(dunit)

    return

  end subroutine make_hom_composite2

end module rfcomp2
