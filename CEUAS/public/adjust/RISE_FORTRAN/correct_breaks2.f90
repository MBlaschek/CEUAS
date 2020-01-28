!!filename=trim(rcpara%prefix)//'igrafgbin'//igcstatnr
!!call write_igrasonde_daily_new(iunit,'igrafgbin'//igcstatnr,rcpara%nmax,rcpara%pmax,rcpara%parmax,itm,itanm,itfgm,itfg12m,rcpara%miss_val,ierr)


module correct_breaks2

use homtests
use homtestsamp
use rfcorio
use rfcor
use txtnc
use correct_mr

contains


!! correct_break is the main routine governing the correction
!!  
!! Leopold Haimberger, 08072004
!!
!! ni=days
!! pmax=pressure levels
!! parmax=launches/day (usually 2)
!! tm=observed T
!! tfgm=bg-obs T
!! solarangles=solar elevation in degrees
!! stm=observed composite T
!! stfgm=bg-obs composite
!! tnumbers=number of obs used for composite
!! cardsmeta=CARDS metadata
!! era40meta=ERA-40 metadata
!!
!! externals: bayes_break,apriori_prob,calc_mean_temp,snhtmov
!! externals: G08AHF (Mann Whitney Test) G08CDF (Komolgoroff Smirnoff Test)

subroutine correct_break(cname,istat,lon,lat,rcpara,tmorig,tfgmorig,solarangles,stm, stfgmorig, tnumbers,tfg,tgps,rasocorrs,eracorrs,crut2,dailycrut2,&
delemask,cardsmeta, era40meta,apriori_probs,apriori_probs_rad,iteration,needs_composite,ini_correct,midx,lasts,gcount,protunit,iunit,alt) !

implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer protunit,ic,istart,istop,its,iunit,i,ib,iname,ilayer,l,idiv,k
integer ipar,ip,iter,isurf


integer,intent(in) :: iteration,istat
real(kind=JPRM) :: tm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tfgm(rcpara%nmax, rcpara%pmax,rcpara%parmax),stfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: tfg(rcpara%nmax,rcpara%pmax,rcpara%parmax),tgps(rcpara%nmax,rcpara%pmax,rcpara%parmax),alt
real(kind=JPRM),intent(in) :: tmorig(rcpara%nmax,rcpara%pmax,rcpara%parmax), tfgmorig(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: weights(rcpara%pmax)
real(kind=JPRM)  :: solarangles(rcpara%nmax,rcpara%parmax)
real(kind=JPRM),intent(in) :: lon,lat,eracorrs(rcpara%nmax,rcpara%pmax, rcpara%parmax), stm(rcpara%nmax,rcpara%pmax,rcpara%parmax), stfgmorig(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM),intent(out) :: rasocorrs(rcpara%nmax,rcpara%pmax,rcpara%parmax), delemask(rcpara%nmax,rcpara%pmax,rcpara%parmax)
integer(kind=JPRM) :: tnumbers(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: plus(rcpara%nmax,rcpara%pmax,rcpara%parmax), minus(rcpara%nmax,rcpara%pmax,rcpara%parmax),prms(rcpara%nmax,rcpara%pmax,rcpara%parmax),mrms(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: splus(rcpara%nmax,rcpara%pmax,rcpara%parmax), sminus(rcpara%nmax,rcpara%pmax,rcpara%parmax),sprms(rcpara%nmax,rcpara%pmax,rcpara%parmax),smrms(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: rplus(rcpara%nmax,rcpara%pmax,2), rminus(rcpara%nmax,rcpara%pmax,2) !!,rprms(rcpara%nmax,rcpara%pmax,1),rmrms(rcpara%nmax,rcpara%pmax,1)
real(kind=JPRM) :: radplus(rcpara%nmax,rcpara%pmax,2), radminus(rcpara%nmax,rcpara%pmax,2) !!,rprms(rcpara%nmax,rcpara%pmax,1),rmrms(rcpara%nmax,rcpara%pmax,1)
real(kind=JPRM) :: bgrplus(rcpara%nmax,rcpara%pmax,1), bgrminus(rcpara%nmax,rcpara%pmax,1), bgrprms(rcpara%nmax,rcpara%pmax,1), bgrmrms(rcpara%nmax,rcpara%pmax,1)
real(kind=JPRM) :: sbgrplus(rcpara%nmax,rcpara%pmax,1), sbgrminus(rcpara%nmax,rcpara%pmax,1),sbgrprms(rcpara%nmax,rcpara%pmax,1),sbgrmrms(rcpara%nmax,rcpara%pmax,1)
real(kind=JPRM) :: plusmean(rcpara%nmax,rcpara%parmax), minusmean(rcpara%nmax,rcpara%parmax),prmsmean(rcpara%nmax,rcpara%parmax),mrmsmean(rcpara%nmax,rcpara%parmax)
real(kind=JPRM) :: splusmean(rcpara%nmax,rcpara%parmax), sminusmean(rcpara%nmax,rcpara%parmax),sprmsmean(rcpara%nmax,rcpara%parmax),smrmsmean(rcpara%nmax,rcpara%parmax) 
real(kind=JPRM) :: rplusmean(rcpara%nmax,1),rminusmean(rcpara%nmax,1), rprmsmean(rcpara%nmax,1),rmrmsmean(rcpara%nmax,1),tsahilf(rcpara%nmax),tsahilfx(rcpara%nmax),tsahilf2(rcpara%nmax)
real(kind=JPRM) :: bgrplusmean(rcpara%nmax,1),bgrminusmean(rcpara%nmax,1), bgrprmsmean(rcpara%nmax,1),bgrmrmsmean(rcpara%nmax,1)
real(kind=JPRM) :: tmean(rcpara%nmax,rcpara%parmax,2), tfgmean(rcpara%nmax,rcpara%parmax,2),stfgmean(rcpara%nmax,rcpara%parmax,2), fgmean(rcpara%nmax,rcpara%parmax,2),stmean(rcpara%nmax,rcpara%parmax,2), tmhilf(rcpara%nmax,rcpara%parmax)
real(kind=JPRM) :: compplus(rcpara%nmax,rcpara%pmax,rcpara%parmax), compminus(rcpara%nmax,rcpara%pmax,rcpara%parmax)

integer(kind=JPRM) :: pcount(rcpara%nmax,rcpara%pmax,rcpara%parmax), mcount(rcpara%nmax,rcpara%pmax,rcpara%parmax),spcount(rcpara%nmax,rcpara%pmax, rcpara%parmax),smcount(rcpara%nmax,rcpara%pmax,rcpara%parmax)
integer(kind=JPRM) :: rpcount(rcpara%nmax,rcpara%pmax,1), rmcount(rcpara%nmax,rcpara%pmax,1),bgrpcount(rcpara%nmax,rcpara%pmax,1), bgrmcount(rcpara%nmax,rcpara%pmax,1),sbgrpcount(rcpara%nmax,rcpara%pmax,1), sbgrmcount(rcpara%nmax,rcpara%pmax,1)
integer(kind=JPRM) :: pcountmean(rcpara%nmax,rcpara%parmax,rcpara%parmax), mcountmean(rcpara%nmax,rcpara%parmax,rcpara%parmax)
!!,sfgpcountmean(rcpara%nmax,rcpara%parmax,rcpara%parmax),sfgmcountmean(rcpara%nmax,rcpara%parmax,rcpara%parmax)
!!integer(kind=JPRM) :: spcountmean(rcpara%nmax,rcpara%parmax),smcountmean(rcpara%nmax,rcpara%parmax)
!!integer(kind=JPRM) :: rpcountmean(rcpara%nmax,1),rmcountmean(rcpara%nmax,1),bgrpcountmean(rcpara%nmax,1),bgrmcountmean(rcpara%nmax,1)

real(kind=JPRM) :: tsa(rcpara%nmax,rcpara%pmax,rcpara%probmax), tsa_of_mean(rcpara%nmax,rcpara%probmax),tsa_of_meanx(rcpara%nmax,rcpara%probmax),mean_of_tsa(rcpara%nmax,rcpara%probmax)
!!real(kind=JPRM) :: ks(rcpara%nmax,rcpara%pmax,rcpara%probmax),ksmean(rcpara%nmax,rcpara%probmax),meanks(rcpara%nmax,rcpara%probmax)
!!real(kind=JPRM) :: mw(rcpara%nmax,rcpara%pmax,rcpara%probmax),mwmean(rcpara%nmax,rcpara%probmax),meanmw(rcpara%nmax,rcpara%probmax)
integer(kind=JPRM) :: mean_of_tsabreaks(rcpara%nmax,rcpara%probmax), tsa_of_meanbreaks(rcpara%nmax,rcpara%probmax)
integer(kind=JPRM) :: meanprobbreaks(rcpara%nmax,rcpara%probmax), probmeanbreaks(rcpara%nmax,rcpara%probmax),probbreaks(rcpara%nmax,rcpara%probmax),chosenbreaks(rcpara%nmax)
real(kind=JPRM) :: breakmean_probs(rcpara%nmax,rcpara%probmax), meanbreak_probs(rcpara%nmax,rcpara%probmax),break_probs(rcpara%nmax,rcpara%probmax)

real(kind=JPRM) :: apriori_probs(rcpara%nmax),apriori_probs_rad(rcpara%nmax), attributionfactors(rcpara%nmax)
integer(kind=JPRM) :: cardsmeta(rcpara%nmax,1,rcpara%nmeta),era40meta(rcpara%nmax), noaa4start,noaa4end

real(kind=JPRM) :: hilfdn(rcpara%nmax,rcpara%parmax),hilf1d(rcpara%nmax), hilf2d(rcpara%nmax,rcpara%pmax),null(rcpara%nmax)
real(kind=JPRM) :: tsa_break_threshs(rcpara%probmax), prob_break_threshs(rcpara%probmax), bprof(1,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: breakprofile(rcpara%pmax,rcpara%parmax),break_thresh,f_val(rcpara%parmax)

real(kind=JPRM) :: crut2(72,36,rcpara%mmax),dailycrut2(rcpara%nmax,rcpara%parmax), manom(rcpara%mmax,rcpara%parmax),tmmon(rcpara%mmax,rcpara%pmax,rcpara%parmax),tmonmean(rcpara%mmax),bgmonmean(rcpara%mmax)
real(kind=JPRM) :: mbganom(rcpara%mmax,rcpara%parmax),bgmon(rcpara%mmax,rcpara%pmax,rcpara%parmax),diff(rcpara%mmax,rcpara%parmax),dailybgdiff(rcpara%nmax,rcpara%parmax),bgdiff(rcpara%mmax,rcpara%parmax)

real(kind=JPRM) :: rasobreaks(rcpara%nmax,rcpara%pmax,rcpara%parmax), rasobreakuncertainties(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: pfill,mfill,swap,mcsave,climate(365,rcpara%parmax), anom(rcpara%nmax,rcpara%parmax),bgclimate(365,rcpara%parmax),bganom(rcpara%nmax,rcpara%parmax)

integer :: err,mindex,xlev,new_maxlen,new_max_miss,j,left_maxlen, right_maxlen, left, right,philf(rcpara%nmax),mhilf(rcpara%nmax),philfx(rcpara%nmax),mhilfx(rcpara%nmax),mgood(rcpara%pmax,rcpara%parmax,3)
integer :: critical_dates(20),ncritical,cbi
integer         :: lasts(rcpara%brmax,rcpara%parmax), midx(rcpara%brmax,rcpara%parmax)
integer,intent(in) :: gcount(rcpara%parmax)
integer :: lim(2,2),good(2),gmin,gmax,indext799

integer :: tsamaxloc(rcpara%probmax)
real(kind=JPRM) :: tsamax(rcpara%probmax)

integer :: needs_composite(:,:)
logical :: use_composite,risky(rcpara%pmax)

logical :: mask(rcpara%nmax),ini_correct(rcpara%parmax),in_gap,ini_correct2(rcpara%parmax)

character*6  :: cname
character*2    ::ch
character :: filename*80

!segmentation_mean variables
real*8::ddiff(rcpara%nmax)
logical:: vh=.true.,ldebug
real*8:: j_est(rcpara%brmax)
real*8:: t_est(rcpara%brmax*rcpara%brmax)
integer :: gindex(rcpara%nmax),ngood

ldebug=.false.

critical_dates=0
ncritical=0
do i=1,rcpara%nmax
  if(era40meta(i) .gt. 0) then
    ncritical=ncritical+1
    critical_dates(ncritical)=i
  endif
enddo
!!ncritical=0

indext799=rcpara%nmax !!toindex(20060131,rcpara)

null=0.

!!$ call omp_set_lock(omp_lp)
read(cname,'(I6)') iname
!!$ call omp_unset_lock(omp_lp)

!!if(iname .gt. 70000 .and. iname .lt. 71000 .or.  iname .gt. 72000 .and. iname .lt. 75000) then
!!  where(apriori_probs_rad(8000:10000) .eq. 0.1)
!!    apriori_probs_rad(8000:10000)=0.
!!  endwhere
     !!$ call omp_set_lock(omp_lp)
!!       write(*,*) cname,' apriori_probs_rad set to zero'
     !!$ call omp_unset_lock(omp_lp)
!!endif

!! no metadata
!!apriori_probs_rad=rcpara%ap_prob_default*5.
!!apriori_probs=rcpara%ap_prob_default

if(iname .eq. 62378) then
  iname=iname
endif

iter=0
 chosenbreaks(1)=1

tfgm=tfgmorig
tm=tmorig
where(tfgmorig .ne. rcpara%miss_val)
  tfgm=tfgmorig+eracorrs !backgroundkorrektur
endwhere

where(tfg .ne. rcpara%miss_val .and. eracorrs .ne. rcpara%miss_val)
  tfg=tfg !-eracorrs
elsewhere
  tfg=rcpara%miss_val
endwhere

if(any(stfgmorig .ne. rcpara%miss_val .and. abs(stfgmorig) .gt. 25.)) then
  !!$ call omp_set_lock(omp_lp)
  write(*,*) cname,': before adjust stfgmorig invalid',maxval(stfgmorig),minval(stfgmorig)
  !!$ call omp_unset_lock(omp_lp)
endif
!! correct most recent part of series
!! do this after breakpoint correction since then a longer interval for more accurate
!! estimation can be used.
if(iter .eq. 0 ) then
!!  call detect_gaps(rcpara,iname,tfgm,midx,lasts,gcount)
  if(any(gcount .gt. 0)) then
      where(tfgm .ne. rcpara%miss_val .and. stfgmorig .ne. rcpara%miss_val)
         stfgm=tfgm-stfgmorig
      elsewhere
         stfgm=rcpara%miss_val
      endwhere
  endif
endif

!!where(tfgmorig .ne. rcpara%miss_val)
!!  tfg=tm+tfgm
!!elsewhere
!!  tfg=rcpara%miss_val
!!endwhere

!! exclude NOAA-4 period
i=1
do while(rcpara%year(i)*10000+rcpara%month(i)*100+rcpara%day(i) .ne. 19750101) 
  i=i+1
enddo
noaa4start=i-1
do while(rcpara%year(i)*10000+rcpara%month(i)*100+rcpara%day(i) .ne. 19760901)   
  i=i+1
enddo
noaa4end=i-1
!!write(*,*) 'Noaa-4 period',todate(noaa4start,rcpara),todate(noaa4end,rcpara),' excluded'
!!tfgm(noaa4start:noaa4end,:,:)=rcpara%miss_val

rasocorrs=0.
rasobreakuncertainties=0.
rasobreaks=0.
!!tsa=rcpara%miss_val
!!plus=rcpara%miss_val
!!minus=rcpara%miss_val
!!prms=rcpara%miss_val
!!mrms=rcpara%miss_val
!!pcount=0
!!mcount=0

!call expand(rcpara,crut2(floor((lon+180.)/5)+1,floor((90-lat)/5)+1,:),dailycrut2)

do while(iter .lt. 1 .and. chosenbreaks(1) .ne. 0) !! correct until no further break is found


if(iter .eq. 0) then
  istart=1
  istop=rcpara%nmax
  ic=1
else
!   if(.false. .and. rcpara%prob_method .ne. 2) then
!     ic=count(chosenbreaks .gt. 0)
!     istart=minval(chosenbreaks(1:ic))-rcpara%snht_maxlen
!     if(istart .lt. 1) istart=1
!     istop=maxval(chosenbreaks(1:ic))+rcpara%snht_maxlen+rcpara%snht_increment
!     if(istop .gt. rcpara%nmax) istop=rcpara%nmax
!   endif
endif

call makemonthly(rcpara,tm,tmmon,9) !in file rfcorio.f90 line 1783
call makemonthly(rcpara,tfg,bgmon,9) !in file rfcorio.f90 line 1783

do ipar=1,rcpara%parmax
  call anomaly(tmmon(:,14,ipar),manom(:,ipar),tmonmean,rcpara%mmax,12,rcpara%miss_val,5) !in file rfmod.f90 line 2311
  call anomaly(bgmon(:,14,ipar),mbganom(:,ipar),bgmonmean,rcpara%mmax,12,rcpara%miss_val,5) !in file rfmod.f90 line 2311
  if (.not. any(tmmon(:,14,ipar) .ne. rcpara%miss_val)) then
    call anomaly(tmmon(:,14,ipar),manom(:,ipar),tmonmean,rcpara%mmax,12,rcpara%miss_val,5) !in file rfmod.f90 line 2311
    call anomaly(bgmon(:,14,ipar),mbganom(:,ipar),bgmonmean,rcpara%mmax,12,rcpara%miss_val,5) !in file rfmod.f90 line 23  
  endif
  i=floor((lon+180.)/5)+1
  j=floor((lat+90.)/5)+1
  where(manom(:,ipar) .ne. rcpara%miss_val .and. crut2(i,j,:) .ne. rcpara%miss_val)
    diff(:,ipar)=-(manom(:,ipar)-crut2(i,j,:))
  elsewhere
    diff(:,ipar)=rcpara%miss_val
  endwhere
  where(mbganom(:,ipar) .ne. rcpara%miss_val .and. crut2(i,j,:) .ne. rcpara%miss_val)
    bgdiff(:,ipar)=-(mbganom(:,ipar)-crut2(i,j,:))
  elsewhere
    bgdiff(:,ipar)=rcpara%miss_val
  endwhere
  call expand(rcpara,diff(:,ipar),dailycrut2(:,ipar)) !in file rfcor.f90 line 88
  call expand(rcpara,bgdiff(:,ipar),dailybgdiff(:,ipar)) !in file rfcor.f90 line 88
enddo



 chosenbreaks=0

break_probs=0.
breakmean_probs=0.
meanbreak_probs=0.

!! layer1 = 200-850 hPa
!! layer2 = 50-150 hPa
mcountmean=0
pcountmean=0
mcount=0
tsa_of_mean=rcpara%miss_val
tsa_of_meanx=rcpara%miss_val
lim(:,1)=(/2,7/)
lim(:,2)=(/10,13/)
call analyze_ts(iname,tfgm,tm,stm,tfg,stfgm,dailycrut2,lim,rcpara,protunit,apriori_probs,apriori_probs_rad,tsa_of_mean,breakmean_probs)

  if(iname .eq. 1001) then
    write(*,*) 'test',count(tm .ne. rcpara%miss_val),count(tfgm .ne. rcpara%miss_val),count(stm .ne. rcpara%miss_val),count(tfg .ne. rcpara%miss_val),count(stfgm .ne. rcpara%miss_val),count(dailycrut2 .ne. rcpara%miss_val),count(tsa_of_mean .ne. rcpara%miss_val)
  endif



mean_of_tsa=0.
breakmean_probs(:,6)=0.

!tsa_of_meanx=tsa_of_meanx-tsa_of_mean

prob_break_threshs=0.5
if(iname .eq. 40745) then
write(*,*) 'test'
endif
call locate_combined_breaks(rcpara,iname,breakmean_probs,probbreaks,tsa_of_mean,mean_of_tsa,prob_break_threshs,rcpara%break_fak-rcpara%break_fak,apriori_probs,iteration,protunit) !in file rfcor.f90 line 939
!! select breaks to be corrected
!if(iter .eq. 1) probbreaks(:,7:10)=probbreaks(:,1:4)

call select_breaks2(rcpara,probbreaks,probbreaks-probbreaks,cardsmeta,era40meta,chosenbreaks,protunit) !in file rfcor.f90 line 1154

ib=count(chosenbreaks .gt. 0)
!! bubble sort chosenbreaks
do i=1,ib
  do j=i,ib
    if(chosenbreaks(j) .gt. chosenbreaks(i)) then
      swap=chosenbreaks(j)
      chosenbreaks(j)=chosenbreaks(i)
      chosenbreaks(i)=swap
    endif
  enddo
enddo
if(iname .eq. 51463) then
write(*,*) 'test'
endif

if(ib .gt. 0) then
  write(ch,'(I2)') ib 
  write(*,'(A11,2F8.2,'//ch//'I9)') cname//' tsa:',maxval(tsa_of_mean(:,1:10)),maxval(breakmean_probs(:,1:10)),rcpara%year(chosenbreaks(1:ib))*10000+rcpara%month(chosenbreaks(1:ib))*100+rcpara%day(chosenbreaks(1:ib))
  write(*,'(A11,10I6,10I6)') cname//'bmp:',probbreaks(1,1:10),probbreaks(2,1:10)
else
  write(*,'(A6,2F8.2,a)') cname//' tsa:',maxval(tsa_of_mean(:,1:10)),maxval(breakmean_probs(:,1:10)),'no break'
endif


! if(.false.) then
!     call average_rad(rcpara%nmax,stfgm(:,4,1),hilf1d,0,rcpara%miss_val) !in file rfcor.f90 line 1583
!    !composit innovs
!     call snhteqsamp2(stfgm(:,4,2),hilf1d,rcpara%nmax,istart,istop, rcpara%snht_maxlen, rcpara%snht_increment, rcpara%miss_val, rcpara%max_miss, critical_dates, 0, & 
!     mean_of_tsa(:,1),plusmean(:,1),minusmean(:,1),prmsmean(:,1),mrmsmean(:,1),philf,mhilf) !in file homtests.f90 line 654
!     call average_rad(rcpara%nmax,tfgm(:,4,1),hilf1d,0,rcpara%miss_val) !in file rfcor.f90 line 1583
!    !Innovs
!     call snhteqsamp2(tfgm(:,4,2),hilf1d,rcpara%nmax,istart,istop,rcpara%snht_maxlen,rcpara%snht_increment,rcpara%miss_val,rcpara%max_miss,critical_dates,0, & 
!     mean_of_tsa(:,2),plusmean(:,1),minusmean(:,1),prmsmean(:,1),mrmsmean(:,1),philf,mhilf) !in file homtests.f90 line 654
!     call average_rad(rcpara%nmax,stfgmorig(:,4,1),hilf1d,0,rcpara%miss_val) !in file rfcor.f90 line 1583
!    !originale composit Innovs
!     call snhteqsamp2(stfgmorig(:,4,2),hilf1d,rcpara%nmax,istart,istop,rcpara%snht_maxlen,rcpara%snht_increment,rcpara%miss_val,rcpara%max_miss,critical_dates,0, & 
!     mean_of_tsa(:,3),plusmean(:,1),minusmean(:,1),prmsmean(:,1),mrmsmean(:,1),philf,mhilf) !in file homtests.f90 line 654
!     where(tmorig(:,4,:) .ne. rcpara%miss_val .and. stfgm(:,4,:) .ne. rcpara%miss_val)
!       stfgm(:,15,:)= tmorig(:,4,:)+stfgm(:,4,:)
!     elsewhere
!       stfgm(:,15,:)= rcpara%miss_val
!     endwhere
!     call average_rad(rcpara%nmax,stfgm(:,15,1),hilf1d,0,rcpara%miss_val) !in file rfcor.f90 line 1583
!     call snhteqsamp2(stfgm(:,15,2),hilf1d,rcpara%nmax,istart,istop,rcpara%snht_maxlen,rcpara%snht_increment,rcpara%miss_val,rcpara%max_miss,critical_dates,0, & 
!     mean_of_tsa(:,4),plusmean(:,1),minusmean(:,1),prmsmean(:,1),mrmsmean(:,1),philf,mhilf) !in file homtests.f90 line 654
!     where(mean_of_tsa .eq. rcpara%miss_val) mean_of_tsa=0.
! endif

do i=1,ib
  cbi=chosenbreaks(i)
  radplus(cbi,:,:)=rcpara%miss_val
  radminus(cbi,:,:)=rcpara%miss_val
  bgrplus(cbi,:,:)=rcpara%miss_val
  bgrminus(cbi,:,:)=rcpara%miss_val
  sbgrplus(cbi,:,:)=rcpara%miss_val
  sbgrminus(cbi,:,:)=rcpara%miss_val
  if(cbi .ne. rcpara%nmax-2) then 
!!    write(*,'(A6,4I9)') cname//':',cbi,left_maxlen,right_maxlen,new_max_miss
      if(i .eq. ib) then
        left=rcpara%mean_maxlen
      else
        left=min(cbi-chosenbreaks(i+1),rcpara%mean_maxlen)
      endif
      do while (any(cbi-critical_dates >0 .and. cbi-left-critical_dates <0) .and. left >rcpara%snht_maxlen/2) 
        left=left-10
!        write(*,*) 'left: ',istat,cbi,left
      enddo
      if(left <=rcpara%snht_maxlen/2) left=min(cbi-chosenbreaks(i+1),rcpara%mean_maxlen)
      right=rcpara%mean_maxlen
      do while (any(critical_dates-cbi >0 .and. critical_dates-cbi-right <0) .and. left >rcpara%snht_maxlen/2) 
        right=right-10
!        write(*,*) 'right: ',istat,cbi,right
      enddo
      if(right <=rcpara%snht_maxlen/2) right=min(rcpara%nmax-cbi,rcpara%mean_maxlen) 

      left_maxlen=min(left,rcpara%mean_maxlen)
      right_maxlen=min(rcpara%nmax-cbi,right)
!!      right_maxlen=min(indext799-cbi,rcpara%mean_maxlen)
      if(right_maxlen .lt. 1) right_maxlen=1
      good=0
    do ipar=1,rcpara%parmax
      do k=cbi-left_maxlen,cbi+right_maxlen
        if(k .gt. 0) then
          if(mcountmean(k,ipar,1) .gt. rcpara%snht_maxlen/2-rcpara%max_miss .and. good(ipar) .eq. 0) good(ipar) =k
        endif
      enddo
    enddo
    !!$ call omp_set_lock(omp_lp)
       if(ldebug) write(*,*) cname,cbi,todate(cbi,rcpara),'good: ',good,' min ',cbi-left_maxlen,rcpara%snht_maxlen/2-rcpara%max_miss,left_maxlen
    if(any(good .gt. cbi-left_maxlen)) then
      gmin=minval(good)
      gmax=maxval(good)
      if(gmin .ne. 0) then 
         if(cbi-gmax .gt. rcpara%snht_maxlen/2) then
           left_maxlen=min(cbi-gmax+(rcpara%snht_maxlen/2-rcpara%max_miss),left_maxlen)
           if(ldebug) write(*,*) cname,cbi,todate(cbi,rcpara),'good: new left_maxlen ',left_maxlen
         endif
      endif
    endif
      new_max_miss=min(left_maxlen-(rcpara%snht_maxlen/2-rcpara%max_miss),right_maxlen-(rcpara%snht_maxlen/2-rcpara%max_miss))
    !!$ call omp_unset_lock(omp_lp)
    if(rcpara%nmax .eq. 45000) then
!       if(cbi-left_maxlen-20820 .lt. 0) left_maxlen=cbi-20820
    endif
    do ipar=1,rcpara%parmax

      do ip=1,rcpara%pmax

        if(ipar .eq. 1) then
          call average_rad(rcpara%nmax,tfg(:,ip,1),hilf1d,0,rcpara%miss_val) !in file rfcor.f90 line 1583
	 !Temperatur bg levels
          call meaneqsamp2(tfg(:,ip,2),hilf1d,rcpara%nmax,cbi-left_maxlen,left_maxlen,right_maxlen,rcpara%snht_increment,rcpara%miss_val,new_max_miss,critical_dates,0, &
        tsa(cbi,ip,ipar),bgrplus(cbi,ip,ipar),bgrminus(cbi,ip,ipar),bgrprms(cbi,ip,ipar),bgrmrms(cbi,ip,ipar),bgrpcount(cbi,ip,ipar),bgrmcount(cbi,ip,ipar),rcpara%month) !in file homtests.f90 line 1050

          call average_rad(rcpara%nmax,stm(:,ip,1),hilf1d,0,rcpara%miss_val) !in file rfcor.f90 line 1583
	 !Temperatur composit levels
          call meaneqsamp2(stm(:,ip,2),hilf1d,rcpara%nmax,cbi-left_maxlen,left_maxlen,right_maxlen,rcpara%snht_increment,rcpara%miss_val,new_max_miss,critical_dates,0, &
        tsa(cbi,ip,ipar),sbgrplus(cbi,ip,ipar),sbgrminus(cbi,ip,ipar),sbgrprms(cbi,ip,ipar),sbgrmrms(cbi,ip,ipar),sbgrpcount(cbi,ip,ipar),sbgrmcount(cbi,ip,ipar),rcpara%month) !in file homtests.f90 line 1050
          if(abs(sbgrplus(cbi,ip,ipar)) .gt. 20. .and. sbgrplus(cbi,ip,ipar) .ne. rcpara%miss_val .or. abs(sbgrminus(cbi,ip,ipar)) .gt. 20. .and. sbgrminus(cbi,ip,ipar) .ne. rcpara%miss_val) then
        !!$ call omp_set_lock(omp_lp)
           write(*,*) iname,ip,ipar,'something is wrong with sbgrplus',sbgrplus(cbi,ip,ipar),sbgrminus(cbi,ip,ipar)
        !!$ call omp_unset_lock(omp_lp)
          endif  
          if(iname .lt. 50000 .or. iname .gt. 60000) then 
!!          if(cardsmeta(cbi,1,2) .eq. 0) then
            where(tm(:,ip,:) .gt. rcpara%miss_val)
              splus(:,ip,:)=tmorig(:,ip,:)-rasocorrs(:,ip,:)
            elsewhere
              splus(:,ip,:)=rcpara%miss_val
            endwhere
            call average_rad(rcpara%nmax,splus(:,ip,1),hilf1d,0,rcpara%miss_val) !in file rfcor.f90 line 1583
            radplus(cbi,ip,ipar)=-999.
            radminus(cbi,ip,ipar)=-999.
            call meaneqsamp2(splus(:,ip,2),hilf1d,rcpara%nmax,cbi-left_maxlen,left_maxlen,right_maxlen,rcpara%snht_increment,rcpara%miss_val,new_max_miss,critical_dates,0, &
        tsa(cbi,ip,ipar),radplus(cbi,ip,ipar),radminus(cbi,ip,ipar),sbgrprms(cbi,ip,ipar),sbgrmrms(cbi,ip,ipar),sbgrpcount(cbi,ip,ipar),sbgrmcount(cbi,ip,ipar),rcpara%month) !in file homtests.f90 line 1050
!!           endif
           endif
        endif

!        use_composite=lat .lt. -50. .or. any(abs(critical_dates-cbi) .lt. rcpara%snht_maxlen) .and. cbi .gt. rcpara%snht_maxlen+2 !!
        use_composite=any(abs(critical_dates-cbi) .lt. rcpara%snht_maxlen/2) .and. cbi .gt. rcpara%snht_maxlen+2 !!

!! .and. sbgrplus(cbi,ip,1) .ne. rcpara%miss_val .and. sbgrminus(cbi,ip,1).ne. rcpara%miss_val .and. abs(bgrplus(cbi,ip,1)-bgrminus(cbi,ip,1)) .ge. abs(sbgrplus(cbi,ip,1)-sbgrminus(cbi,ip,1))

          where(stfgm(:,ip,ipar) .ne. rcpara%miss_val)          
            hilf1d=stfgm(:,ip,ipar)
          elsewhere
            hilf1d=rcpara%miss_val
          endwhere
	 ! Composit innovs
        call meaneqsamp2(hilf1d,null,rcpara%nmax,cbi-left_maxlen,left_maxlen,right_maxlen,rcpara%snht_increment,rcpara%miss_val,new_max_miss,critical_dates,0, &
        tsa(cbi,ip,ipar),splus(cbi,ip,ipar),sminus(cbi,ip,ipar),sprms(cbi,ip,ipar),smrms(cbi,ip,ipar),spcount(cbi,ip,ipar),smcount(cbi,ip,ipar),rcpara%month) !in file homtests.f90 line 1050
        compplus(cbi,ip,ipar)=splus(cbi,ip,ipar)
        compminus(cbi,ip,ipar)=sminus(cbi,ip,ipar)

        if(cbi .eq. toindex(19880202,rcpara)) then
          write(*,*) use_composite
        endif
        if(use_composite) then
          where(stfgm(:,ip,ipar) .ne. rcpara%miss_val)          
            hilf1d=stfgm(:,ip,ipar)
          elsewhere
            hilf1d=rcpara%miss_val
          endwhere
        !!$ call omp_set_lock(omp_lp)
           if(ipar .eq. 1 .and. ldebug) write(*,*) cname,todate(cbi,rcpara),ip,abs(bgrplus(cbi,ip,ipar)-bgrminus(cbi,ip,ipar)),abs(sbgrplus(cbi,ip,ipar)-sbgrminus(cbi,ip,ipar)),' used composite'
        !!$ call omp_unset_lock(omp_lp)
        else
          where(tfgm(:,ip,ipar) .ne. rcpara%miss_val)          
            hilf1d=tfgm(:,ip,ipar)
          elsewhere
            hilf1d=rcpara%miss_val
          endwhere
        !!$ call omp_set_lock(omp_lp)
           if(ipar .eq. 1 .and. ldebug) write(*,*) cname,todate(cbi,rcpara),ip,abs(bgrplus(cbi,ip,ipar)-bgrminus(cbi,ip,ipar)),abs(sbgrplus(cbi,ip,ipar)-sbgrminus(cbi,ip,ipar)),' abandoned composite'
        !!$ call omp_unset_lock(omp_lp)
	!Innovs mit bg oder composite
          call meaneqsamp2(hilf1d,null,rcpara%nmax,cbi-left_maxlen,left_maxlen,right_maxlen,rcpara%snht_increment,rcpara%miss_val,new_max_miss,critical_dates,ncritical, &
        tsa(cbi,ip,ipar),splus(cbi,ip,ipar),sminus(cbi,ip,ipar),sprms(cbi,ip,ipar),smrms(cbi,ip,ipar),spcount(cbi,ip,ipar),smcount(cbi,ip,ipar),rcpara%month) !in file homtests.f90 line 1050
        endif


        !!$ call omp_set_lock(omp_lp)
if(ldebug) write(*,*) cname,todate(cbi,rcpara),ip,ipar,splus(cbi,ip,ipar),sminus(cbi,ip,ipar),' splus',spcount(cbi,ip,ipar),smcount(cbi,ip,ipar),min(left_maxlen,right_maxlen)-new_max_miss
        !!$ call omp_unset_lock(omp_lp)

        if(use_composite .and. (splus(cbi,ip,ipar) .eq. rcpara%miss_val .or. sminus(cbi,ip,ipar) .eq. rcpara%miss_val)) then

          where(tfgm(:,ip,ipar) .ne. rcpara%miss_val)          
            hilf1d=tfgm(:,ip,ipar)
          elsewhere
            hilf1d=rcpara%miss_val
          endwhere
    !Innovs wenn keine Mittelwerte für Composit
        call meaneqsamp2(hilf1d,null,rcpara%nmax,cbi-left_maxlen,left_maxlen,right_maxlen,rcpara%snht_increment,rcpara%miss_val,new_max_miss,critical_dates,0, &
        tsa(cbi,ip,ipar),splus(cbi,ip,ipar),sminus(cbi,ip,ipar),sprms(cbi,ip,ipar),smrms(cbi,ip,ipar),spcount(cbi,ip,ipar),smcount(cbi,ip,ipar),rcpara%month) !in file homtests.f90 line 1050
        !!$ call omp_set_lock(omp_lp)
if(ldebug) write(*,*) cname,todate(cbi,rcpara),ip,ipar,splus(cbi,ip,ipar),sminus(cbi,ip,ipar),'abandoned composite, retry with bg: splus',spcount(cbi,ip,ipar),smcount(cbi,ip,ipar),min(left_maxlen,right_maxlen)-new_max_miss
        !!$ call omp_unset_lock(omp_lp)
        endif

        plus(cbi,ip,ipar)=splus(cbi,ip,ipar)
        minus(cbi,ip,ipar)=sminus(cbi,ip,ipar)
        if(splus(cbi,ip,ipar) .ne. rcpara%miss_val .and. sminus(cbi,ip,ipar) .ne. rcpara%miss_val) then 
          rplus(cbi,ip,ipar)=1.96*sqrt((sprms(cbi,ip,ipar)-splus(cbi,ip,ipar)*splus(cbi,ip,ipar))/spcount(cbi,ip,ipar))
          rminus(cbi,ip,ipar)=1.96*sqrt((smrms(cbi,ip,ipar)-sminus(cbi,ip,ipar)*sminus(cbi,ip,ipar))/spcount(cbi,ip,ipar))
        else
         rplus(cbi,ip,ipar)=0. 
         rminus(cbi,ip,ipar)=0. 
        endif

        if(ipar .eq. 1) then
         if(use_composite .and. sbgrplus(cbi,ip,ipar) .ne. rcpara%miss_val .and. sbgrminus(cbi,ip,ipar) .ne. rcpara%miss_val) then
           bgrplus(cbi,ip,ipar)=sbgrplus(cbi,ip,ipar)
           bgrminus(cbi,ip,ipar)=sbgrminus(cbi,ip,ipar)
         endif
        endif

      enddo
    enddo

    do ip=1,rcpara%pmax
      risky(ip)=.false.
      do ipar=1,rcpara%parmax
        if(smcount(cbi,ip,ipar) .gt. min(left_maxlen,right_maxlen)-new_max_miss .and. smcount(cbi,ip,ipar) .lt. rcpara%snht_maxlen/2. &
.and. smcount(cbi,ip,3-ipar) .lt. min(left_maxlen,right_maxlen)-new_max_miss .and. smcount(cbi,ip,3-ipar) .gt. 10) risky(ip)=.true.
!!        if(spcount(cbi,ip,ipar) .gt. min(left_maxlen,right_maxlen)-new_max_miss .and. spcount(cbi,ip,ipar) .lt. rcpara%snht_maxlen/2. .and. spcount(cbi,ip,3-ipar) .lt. min(left_maxlen,right_maxlen)-new_max_miss .and. spcount(cbi,ip,3-ipar) .gt. 10) risky=.true.
      enddo


   enddo
      if(count(risky) .lt. -3) then 
        do ip=1,rcpara%pmax
          if(risky(ip)) then
        !!$ call omp_set_lock(omp_lp)
        if(ldebug) write(*,*) cname,todate(cbi,rcpara),ip,ipar,smcount(cbi,ip,:),'adjustment too risky'
        splus(cbi,ip,:)=rcpara%miss_val
        plus(cbi,ip,:)=rcpara%miss_val
        sminus(cbi,ip,:)=rcpara%miss_val
        sminus(cbi,ip,:)=rcpara%miss_val
      !!$ call omp_unset_lock(omp_lp)
          endif
        enddo
      endif
     


!! check if unadjusted anomalies in 850 hPa agree with HadAT anomalies  
!! if yes or if bg anomalies disagree with HadAT anomalies, adjustment of lower levels is much more conservative
!! in adjust_series
!!
!! Leo Haimberger, 11 05 2007

    splus(:,15:16,:)=-999.
    sminus(:,15:16,:)=-999.
    new_max_miss=min(left_maxlen-3*(rcpara%snht_maxlen/2-rcpara%max_miss),right_maxlen-3*(rcpara%snht_maxlen/2-rcpara%max_miss))
    if(new_max_miss .lt. 1) new_max_miss=1
    where(tfgm(:,14,:).eq. rcpara%miss_val)
         dailycrut2=rcpara%miss_val
         dailybgdiff=rcpara%miss_val
    endwhere
    do ipar=1,rcpara%parmax
       call meaneqsamp2(null,dailycrut2(:,ipar),rcpara%nmax,cbi-left_maxlen,left_maxlen,right_maxlen,rcpara%snht_increment,rcpara%miss_val,new_max_miss,critical_dates,0, &
        tsa(cbi,15,ipar),splus(cbi,15,ipar),sminus(cbi,15,ipar),sprms(cbi,15,ipar),smrms(cbi,15,ipar),spcount(cbi,15,ipar),smcount(cbi,15,ipar),rcpara%month) !in file homtests.f90 line 1050
       call meaneqsamp2(null,dailybgdiff(:,ipar),rcpara%nmax,cbi-left_maxlen,left_maxlen,right_maxlen,rcpara%snht_increment,rcpara%miss_val,new_max_miss,critical_dates,0, &
        tsa(cbi,16,ipar),splus(cbi,16,ipar),sminus(cbi,16,ipar),sprms(cbi,16,ipar),smrms(cbi,16,ipar),spcount(cbi,16,ipar),smcount(cbi,16,ipar),rcpara%month) !in file homtests.f90 line 1050

        if(splus(cbi,15,ipar) .ne. rcpara%miss_val .and. sminus(cbi,15,ipar) .ne. rcpara%miss_val &
     .and. spcount(cbi,15,ipar) .ne. 0 .and. sprms(cbi,15,ipar)-splus(cbi,15,ipar)*splus(cbi,15,ipar) .gt. 0.) then 

          rplus(cbi,15,ipar)=sqrt(30.)*1.96*sqrt((sprms(cbi,15,ipar)-splus(cbi,15,ipar)*splus(cbi,15,ipar))/spcount(cbi,15,ipar)) 
        else
         rplus(cbi,15,ipar)=0. 
        endif

        if(spcount(cbi,15,ipar) .ne. 0 .and. smrms(cbi,15,ipar)-sminus(cbi,15,ipar)*sminus(cbi,15,ipar) .gt. 0.) then 
          rminus(cbi,15,ipar)=sqrt(30.)*1.96*sqrt((smrms(cbi,15,ipar)-sminus(cbi,15,ipar)*sminus(cbi,15,ipar))/spcount(cbi,15,ipar))
        else
         rminus(cbi,15,ipar)=0. 
        endif
        !!$ call omp_set_lock(omp_lp)
if(ldebug)       write(*,*) 'surface break uncertainty',rplus(cbi,15,ipar),rminus(cbi,15,ipar)
      !!$ call omp_unset_lock(omp_lp)
    enddo

  !! execute the adjustment

!    if(cbi .lt. indext799) then
     
     if (.false. .and. lat .lt. -55. .and. use_composite) then
      call adjust_series_comp(rcpara,iteration,cname,lon,  chosenbreaks,i,attributionfactors,tm,tfgm,stm,stfgm,tgps,cardsmeta,era40meta,eracorrs,rasocorrs,&
rasobreaks,rasobreakuncertainties,delemask,plus,minus,splus,sminus,rplus,rminus,compplus,compminus,radplus,radminus,bgrplus,bgrminus,plusmean,minusmean,dailycrut2,midx,lasts,gcount,protunit,alt) 
!bei use_composit sind bgrplus und bgrminus gleich sbgrplus und sbgrminus und plus =splus, minus=sminus
     else
      call adjust_series(rcpara,iteration,cname,lon, chosenbreaks,i,attributionfactors,tm,tfgm,stm,stfgm,tgps,cardsmeta,era40meta,eracorrs,rasocorrs,&
rasobreaks,rasobreakuncertainties,delemask,plus,minus,splus,sminus,rplus,rminus,compplus,compminus,radplus,radminus,bgrplus,bgrminus,plusmean,minusmean,dailycrut2,midx,lasts,gcount,protunit,alt) 
     endif
!    endif

  endif
enddo


if(iter .eq. 0 .and. (iteration .gt. 1 .or. iteration .eq. rcpara%maxiter)) then
    istart=floor((rcpara%mmax/12-2)*365.25) 
    write(*,*) 'late break after',istart
    do ipar=1,rcpara%parmax
     do i=istart,rcpara%nmax
      if(rasocorrs(i,8,ipar) .ne. rasocorrs(i-1,8,ipar) .and. rasocorrs(i-1,8,ipar) .ne. 0) then
        ini_correct(ipar)=.true.
        write(*,*) cname,'late break'
        exit
      endif
     enddo
    enddo
!! correct most recent part of series
!! do this after breakpoint correction since then a longer interval for more accurate
!! estimation can be used.
!      call correct_mostrecent(rcpara,iname,midx,lasts,gcount,tfgm,stfgm,needs_composite,ini_correct,rasocorrs,rasobreaks,rasobreakuncertainties) !in this file line681
      if (rcpara%initial_adjust(1:3).eq. 'all') then
         ini_correct2=.true.
      else
         ini_correct2=ini_correct
      endif
      call correct_mostrecent(rcpara,iname,midx,lasts,gcount,tfgm,stfgm,tgps,needs_composite,ini_correct2,rasocorrs,rasobreaks,rasobreakuncertainties) !in this file line681
endif

where(tfgm .ne. rcpara%miss_val)
  tfgm=tfgmorig+rasocorrs
  tm=tmorig-rasocorrs
elsewhere
!!  tfgmorig=rcpara%miss_val
!!  tm=rcpara%miss_val
endwhere
where(tfgm .ne. rcpara%miss_val .and. stfgm .ne. rcpara%miss_val)
  stfgm=stfgm+rasocorrs
elsewhere
!!  tfgmorig=rcpara%miss_val
!!  tm=rcpara%miss_val
endwhere
where(tgps .ne. rcpara%miss_val .and. tgps .ne. rcpara%miss_val)
  tgps=tgps+rasocorrs
elsewhere
!!  tfgmorig=rcpara%miss_val
!!  tm=rcpara%miss_val
endwhere

!! write surface temperature anomaly into dailycrut2
do ipar=1,rcpara%parmax
  call expand(rcpara,crut2(floor((lon+180.)/5)+1,floor((lat+89.999)/5)+1,:) ,dailycrut2(:,ipar)) !in file rfcor.f90 line 88
enddo

if(any(rasocorrs .ne. rcpara%miss_val .and. abs(rasocorrs) .gt. 25.)) then
  !!$ call omp_set_lock(omp_lp)
  write(*,*) cname,': rasocorrs invalid',maxval(rasocorrs),minval(rasocorrs)
  !!$ call omp_unset_lock(omp_lp)
!  call exit(1)
endif

if(rcpara%extended_output .eq. 'Y') then
  filename=TRIM(rcpara%prefix)//cname//'/'//cname//'.dump' 
  if(iname .eq. 1001) then
    write(*,*) 'test',count(tm .ne. rcpara%miss_val),count(tsa_of_mean .ne. rcpara%miss_val)
  endif
  call save_state(iunit,filename,rcpara,tm,mean_of_tsa,tsa_of_mean,apriori_probs,break_probs,breakmean_probs,chosenbreaks) !in file rfcor.f90 line 690
endif 

!!if(iter .eq. 0) then
!!  apriori_probs=rcpara%ap_prob_default
!!else
!!endif

if(iter .eq. 0 .and. chosenbreaks(1) .eq. 0) chosenbreaks(1)=rcpara%nmax
iter=iter+1
enddo

tsa_of_mean=rcpara%miss_val
call analyze_ts(iname,tfgm,tm,stm,tfg,stfgm,dailycrut2,lim,rcpara,protunit,apriori_probs,apriori_probs_rad,tsa_of_mean,breakmean_probs)

tsamax=0
tsamaxloc=0
do i=1,12
 tsamax(i)=maxval(tsa_of_mean(:,i))
 if(tsamax(i) .ne. rcpara%miss_val) tsamaxloc(i)=todate(maxloc(tsa_of_mean(:,i),1),rcpara)
enddo
if(any(tsamax((/1,2,3,4,5/)) .gt. 300)) then
  write(*,'(A5,A9,12F9.1)') cname,': tsamax:',tsamax
  write(*,'(A5,A9,12I9)') cname,': maxloc:',tsamaxloc
endif
!filename=trim(rcpara%prefix)//cname//'/feedbackglobbincorrsave'//cname//'.nc'
!call write_sonde_corr_daily(iunit,filename,rcpara,err,rasobreaks,rasobreakuncertainties)

!call write_sonde_corr_daily_nc(trim(filename),rcpara,istat,err,rasocorrs,lon,lat,rasobreaks,rasobreakuncertainties) !in file read_txt_write_nc line 468
call corr_stats(rcpara,tm,tfgm,stm,stfgm,protunit) !in file rfcor.f90 line 746

return
end subroutine correct_break


!! calc_mean performs pressure averaging or log(pressure) averaging
!! weights can be used to exclude certain levels from averaging
!! useful for levels that are often missing, such as 10,20,250,925,1000 hPa
!!
!! L. Haimberger 08072004
!! 
!! logmean=1 -> log(pressure) average
subroutine calc_mean(rcpara,t,tmean,weights,dim3,logmean,itol) !

implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer i,ipar,min_av
integer,intent(in) :: dim3,itol,logmean
real(kind=JPRM),intent(in) :: weights(:)
real(kind=JPRM) :: player(rcpara%pmax)
real(kind=JPRM),intent(in) :: t(:,:,:)
real(kind=JPRM),intent(out) :: tmean(:,:)
real(kind=JPRM) :: div(rcpara%nmax)
integer :: tcount(rcpara%nmax)

tmean=0.

if(logmean .eq. 0) then
  do i=2,rcpara%pmax-1
    player(i)=(rcpara%plevs(i+1)+rcpara%plevs(i))/2.-(rcpara%plevs(i-1)+rcpara%plevs(i))/2.
  enddo
  player(1)=(rcpara%plevs(2)+rcpara%plevs(1))/2.-rcpara%plevs(1)/2.
  player(rcpara%pmax)=rcpara%plevs(rcpara%pmax)-(rcpara%plevs(rcpara%pmax-1)+rcpara%plevs(rcpara%pmax))/2.
else
  do i=2,rcpara%pmax-1
    player(i)=(log(rcpara%plevs(i+1))+log(rcpara%plevs(i)))/2.-(log(rcpara%plevs(i-1))+log(rcpara%plevs(i)))/2.
  enddo
  player(1)=(log(rcpara%plevs(2))+log(rcpara%plevs(1)))/2.-log(rcpara%plevs(1))/2.
  player(rcpara%pmax)=log(rcpara%plevs(rcpara%pmax))-(log(rcpara%plevs(rcpara%pmax-1))+log(rcpara%plevs(rcpara%pmax)))/2.
endif

min_av=count(weights .gt. 0.)-itol
!!write(*,*) 'min_av',min_av

do ipar=1,dim3
 tcount=0
 div=0.
 do i=1,rcpara%pmax
  if(weights(i) .gt. 0.) then 
    where(t(:,i,ipar) .ne. rcpara%miss_val)
      tmean(:,ipar)=tmean(:,ipar)+player(i)*weights(i)*t(:,i,ipar)
      div=div+player(i)*weights(i) 
      tcount=tcount+1
    endwhere
  endif
 enddo
 where(tcount .ge. maxval(tcount) .and. tcount .gt. 0)
   tmean(:,ipar)=tmean(:,ipar)/div
 elsewhere
   tmean(:,ipar)=rcpara%miss_val
 endwhere
enddo


return
end subroutine calc_mean

subroutine analyze_ts(iname,tfgm,tm,stm,tfg,stfgm,dailycrut2,lim,rcpara,protunit,apriori_probs,apriori_probs_rad,tsa_of_mean,breakmean_probs)

implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer protunit,ic,istart,istop,its,iunit,i,ib,ilayer,l,idiv,k
integer ipar,ip,iter,isurf
integer,intent(in) :: iname,lim(2,2)

real(kind=JPRM),intent(in) :: tm(rcpara%nmax,rcpara%pmax,rcpara%parmax),stm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tfgm(rcpara%nmax, rcpara%pmax,rcpara%parmax),stfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tfg(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: weights(rcpara%pmax)

real(kind=JPRM) :: plus(rcpara%nmax,rcpara%pmax,rcpara%parmax), minus(rcpara%nmax,rcpara%pmax,rcpara%parmax),prms(rcpara%nmax,rcpara%pmax,rcpara%parmax),mrms(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: splus(rcpara%nmax,rcpara%pmax,rcpara%parmax), sminus(rcpara%nmax,rcpara%pmax,rcpara%parmax),sprms(rcpara%nmax,rcpara%pmax,rcpara%parmax),smrms(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: rplus(rcpara%nmax,rcpara%pmax,2), rminus(rcpara%nmax,rcpara%pmax,2)
real(kind=JPRM) :: bgrplus(rcpara%nmax,rcpara%pmax,1), bgrminus(rcpara%nmax,rcpara%pmax,1), bgrprms(rcpara%nmax,rcpara%pmax,1), bgrmrms(rcpara%nmax,rcpara%pmax,1)
real(kind=JPRM) :: sbgrplus(rcpara%nmax,rcpara%pmax,1), sbgrminus(rcpara%nmax,rcpara%pmax,1),sbgrprms(rcpara%nmax,rcpara%pmax,1),sbgrmrms(rcpara%nmax,rcpara%pmax,1)
real(kind=JPRM) :: plusmean(rcpara%nmax,rcpara%parmax), minusmean(rcpara%nmax,rcpara%parmax),prmsmean(rcpara%nmax,rcpara%parmax),mrmsmean(rcpara%nmax,rcpara%parmax)
real(kind=JPRM) :: splusmean(rcpara%nmax,rcpara%parmax), sminusmean(rcpara%nmax,rcpara%parmax),sprmsmean(rcpara%nmax,rcpara%parmax),smrmsmean(rcpara%nmax,rcpara%parmax) 
real(kind=JPRM) :: rplusmean(rcpara%nmax,1),rminusmean(rcpara%nmax,1), rprmsmean(rcpara%nmax,1),rmrmsmean(rcpara%nmax,1),tsahilf(rcpara%nmax),tsahilfx(rcpara%nmax),tsahilf2(rcpara%nmax)
real(kind=JPRM) :: bgrplusmean(rcpara%nmax,1),bgrminusmean(rcpara%nmax,1), bgrprmsmean(rcpara%nmax,1),bgrmrmsmean(rcpara%nmax,1)
real(kind=JPRM) :: tmean(rcpara%nmax,rcpara%parmax,2), tfgmean(rcpara%nmax,rcpara%parmax,2),stfgmean(rcpara%nmax,rcpara%parmax,2), fgmean(rcpara%nmax,rcpara%parmax,2),stmean(rcpara%nmax,rcpara%parmax,2), tmhilf(rcpara%nmax,rcpara%parmax)
real(kind=JPRM) :: compplus(rcpara%nmax,rcpara%pmax,rcpara%parmax), compminus(rcpara%nmax,rcpara%pmax,rcpara%parmax)

integer(kind=JPRM) :: pcount(rcpara%nmax,rcpara%pmax,rcpara%parmax), mcount(rcpara%nmax,rcpara%pmax,rcpara%parmax),spcount(rcpara%nmax,rcpara%pmax, rcpara%parmax),smcount(rcpara%nmax,rcpara%pmax,rcpara%parmax)
integer(kind=JPRM) :: rpcount(rcpara%nmax,rcpara%pmax,1), rmcount(rcpara%nmax,rcpara%pmax,1),bgrpcount(rcpara%nmax,rcpara%pmax,1), bgrmcount(rcpara%nmax,rcpara%pmax,1),sbgrpcount(rcpara%nmax,rcpara%pmax,1), sbgrmcount(rcpara%nmax,rcpara%pmax,1)
integer(kind=JPRM) :: pcountmean(rcpara%nmax,rcpara%parmax,rcpara%parmax), mcountmean(rcpara%nmax,rcpara%parmax,rcpara%parmax)

real(kind=JPRM) :: tsa(rcpara%nmax,rcpara%pmax,rcpara%probmax), tsa_of_meanx(rcpara%nmax,rcpara%probmax),mean_of_tsa(rcpara%nmax,rcpara%probmax)
real(kind=JPRM), intent(out):: tsa_of_mean(rcpara%nmax,rcpara%probmax),breakmean_probs(rcpara%nmax,rcpara%probmax)
integer(kind=JPRM) :: meanprobbreaks(rcpara%nmax,rcpara%probmax), probmeanbreaks(rcpara%nmax,rcpara%probmax),probbreaks(rcpara%nmax,rcpara%probmax),chosenbreaks(rcpara%nmax)
real(kind=JPRM) ::  meanbreak_probs(rcpara%nmax,rcpara%probmax),break_probs(rcpara%nmax,rcpara%probmax)

real(kind=JPRM),intent(in) :: apriori_probs(rcpara%nmax),apriori_probs_rad(rcpara%nmax)

real(kind=JPRM) :: hilfdn(rcpara%nmax,rcpara%parmax),hilf1d(rcpara%nmax), hilf2d(rcpara%nmax,rcpara%pmax),null(rcpara%nmax)
real(kind=JPRM) :: tsa_break_threshs(rcpara%probmax), prob_break_threshs(rcpara%probmax), bprof(1,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: breakprofile(rcpara%pmax,rcpara%parmax),break_thresh,f_val(rcpara%parmax)

real(kind=JPRM),intent(in) :: dailycrut2(rcpara%nmax,rcpara%parmax)

real(kind=JPRM) :: pfill,mfill,swap,mcsave,climate(365,rcpara%parmax), anom(rcpara%nmax,rcpara%parmax),bgclimate(365,rcpara%parmax),bganom(rcpara%nmax,rcpara%parmax)

integer :: err,mindex,xlev,new_maxlen,new_max_miss,j,left_maxlen, right_maxlen, left, right,philf(rcpara%nmax),mhilf(rcpara%nmax),philfx(rcpara%nmax),mhilfx(rcpara%nmax),mgood(rcpara%pmax,rcpara%parmax,3)
integer :: critical_dates(20),ncritical,cbi




ldebug=.false.
tsa_of_mean=rcpara%miss_val
breakmean_probs=rcpara%miss_val
stmean=rcpara%miss_val
stfgmean=rcpara%miss_val
tsahilf=rcpara%miss_val
null=0
critical_dates=0
istart=1
istop=rcpara%nmax

do ilayer=1,2 !259

   do ic=lim(1,ilayer),lim(2,ilayer) !261
     weights=rcpara%mweights
     weights(1:maxval((/1,ic-1/)))=0
     weights(ic+2:rcpara%pmax)=0
     call calc_mean(rcpara,tfgm,tfgmean(:,:,ilayer),weights,rcpara%parmax,1,0) !in this file line 929
     call calc_mean(rcpara,tm,tmean(:,:,ilayer),weights,rcpara%parmax,1,0) !in this file line 929
     call calc_mean(rcpara,stm,stmean(:,:,ilayer),weights,rcpara%parmax,1,0) !in this file line 929
     call calc_mean(rcpara,tfg,fgmean(:,:,ilayer),weights,rcpara%parmax,1,0) !in this file line 929
     call calc_mean(rcpara,stfgm,stfgmean(:,:,ilayer),weights,rcpara%parmax,1,0) !in this file line 929

if(iname .eq. 01001) then
   write(*,*) 'means',count(tfgmean .ne. rcpara%miss_val),count(tmean .ne. rcpara%miss_val),count(stmean .ne. rcpara%miss_val),count(fgmean .ne. rcpara%miss_val),count(stfgmean .ne. rcpara%miss_val)
   write(*,'(A8,16F4.1,3I8)') 'weights:',weights,ic,ilayer,count(tfgm.ne. rcpara%miss_val)
!   do i=1,rcpara%nmax
!      if(any(tfgmean(i,:,ilayer) .ne. rcpara%miss_val)) write(*,*) 'tfgmean ',i,tfgmean(i,:,ilayer)
!   enddo
!   stop
endif

     do ipar=1,rcpara%parmax !273
     !Schichtmittel bg
       call snhteqsamp2(tfgmean(:,ipar,ilayer),null,rcpara%nmax,istart,istop,rcpara%snht_maxlen,rcpara%snht_increment,rcpara%miss_val,rcpara%max_miss,critical_dates,0, & 
         tsahilf,plusmean(:,ipar),minusmean(:,ipar),prmsmean(:,ipar),mrmsmean(:,ipar),philf,mhilf,rcpara%month) !in file homtests.f90 line 654
!       call snhteqsamp2x(tfgmean(:,ipar,ilayer),null,rcpara%nmax,istart,istop,rcpara%snht_maxlen,rcpara%snht_increment,rcpara%miss_val,rcpara%max_miss,critical_dates,0, & 
!         tsahilf,plusmean(:,ipar),minusmean(:,ipar),prmsmean(:,ipar),mrmsmean(:,ipar),philf,mhilf,rcpara%month) !in file homtests.f90 line 654

         do i=1,rcpara%nmax
           if(tsahilf(i) .gt. tsa_of_mean(i,2*(ilayer-1)+ipar)) tsa_of_mean(i,2*(ilayer-1)+ipar)=tsahilf(i)
!           if(tsahilfx(i) .gt. tsa_of_meanx(i,2*(ilayer-1)+ipar)) tsa_of_meanx(i,2*(ilayer-1)+ipar)=tsahilfx(i)
           if(mhilf(i) .gt. mcountmean(i,ipar,ilayer)) mcountmean(i,ipar,ilayer)=mhilf(i)
           if(philf(i) .gt. pcountmean(i,ipar,ilayer)) pcountmean(i,ipar,ilayer)=philf(i)
        enddo

       if(any(stfgmean(:,ipar,ilayer) .ne. rcpara%miss_val)) then 
       !Schichtmittel composit
       call snhteqsamp2(stfgmean(:,ipar,ilayer),null,rcpara%nmax,istart,istop,rcpara%snht_maxlen,rcpara%snht_increment,rcpara%miss_val,rcpara%max_miss,critical_dates,0, & 
         tsahilf,plusmean(:,ipar),minusmean(:,ipar),prmsmean(:,ipar),mrmsmean(:,ipar),philf,mhilf,rcpara%month) !in file homtests.f90 line 654
         do i=1,rcpara%nmax
           if(tsahilf(i) .gt. tsa_of_mean(i,2*(ilayer-1)+ipar+6)) tsa_of_mean(i,2*(ilayer-1)+ipar+6)=tsahilf(i)
!!           if(mhilf(i) .gt. sfgmcountmean(i,ipar,ilayer)) sfgmcountmean(i,ipar,ilayer)=mhilf(i)
!!           if(philf(i) .gt. sfgpcountmean(i,ipar,ilayer)) sfgpcountmean(i,ipar,ilayer)=philf(i)
        enddo
        endif
   
        if(ipar .eq. 1) then
          tsahilf=rcpara%miss_val
          call average_rad(rcpara%nmax,tmean(:,1,ilayer),hilf1d,0,rcpara%miss_val) !in file rfcor.f90 line 1583
	!Strahlung Schicktmittel Radiosonde
            if(rcpara%prob_method .eq. 1) call snhteqsamp2(tmean(:,2,ilayer),hilf1d,rcpara%nmax,istart,istop,rcpara%snht_maxlen,rcpara%snht_increment,rcpara%miss_val,rcpara%max_miss,critical_dates,0,&
    tsahilf,rplusmean(:,ipar),rminusmean(:,ipar),rprmsmean(:,ipar),rmrmsmean(:,ipar),philf,mhilf,rcpara%month) !in file homtests.f90 line 654
         do i=1,rcpara%nmax
           if(tsahilf(i) .gt. tsa_of_mean(i,5)) tsa_of_mean(i,5)=tsahilf(i)
!!           if(mhilf(i) .gt. rmcountmean(i,ipar)) rmcountmean(i,ipar)=mhilf(i)
!!           if(philf(i) .gt. rpcountmean(i,ipar)) rpcountmean(i,ipar)=philf(i)
        enddo
          call average_rad(rcpara%nmax,fgmean(:,1,ilayer),hilf1d,0,rcpara%miss_val) !in file rfcor.f90 line 1583
	 !Strahlung Schichtmittel bg
          if(rcpara%prob_method .eq. 1) call snhteqsamp2(fgmean(:,2,ilayer),hilf1d,rcpara%nmax,istart,istop,rcpara%snht_maxlen,rcpara%snht_increment,rcpara%miss_val,rcpara%max_miss,critical_dates,0, &
    tsahilf,rplusmean(:,ipar),rminusmean(:,ipar),rprmsmean(:,ipar),rmrmsmean(:,ipar),philf,mhilf,rcpara%month) !in file homtests.f90 line 654
         do i=1,rcpara%nmax
           if(tsahilf(i) .gt. tsa_of_mean(i,6)) tsa_of_mean(i,6)=tsahilf(i)
!!           if(mhilf(i) .gt. rmcountmean(i,ipar)) rmcountmean(i,ipar)=mhilf(i)
!!           if(philf(i) .gt. rpcountmean(i,ipar)) rpcountmean(i,ipar)=philf(i)
         enddo
          call average_rad(rcpara%nmax,stmean(:,1,ilayer),hilf1d,0,rcpara%miss_val) !in file rfcor.f90 line 1583
	 !Strahlung Schichtmittel Composit
          if(rcpara%prob_method .eq. 1) call snhteqsamp2(stmean(:,2,ilayer),hilf1d,rcpara%nmax,istart,istop,rcpara%snht_maxlen,rcpara%snht_increment,rcpara%miss_val,rcpara%max_miss,critical_dates,0, &
    tsahilf,rplusmean(:,ipar),rminusmean(:,ipar),rprmsmean(:,ipar),rmrmsmean(:,ipar),philf,mhilf,rcpara%month) !in file homtests.f90 line 654
         do i=1,rcpara%nmax
           if(tsahilf(i) .gt. tsa_of_mean(i,12)) tsa_of_mean(i,12)=tsahilf(i)
!!           if(mhilf(i) .gt. rmcountmean(i,ipar)) rmcountmean(i,ipar)=mhilf(i)
!!           if(philf(i) .gt. rpcountmean(i,ipar)) rpcountmean(i,ipar)=philf(i)
         enddo
        endif
      enddo !273
   enddo !261

!! compare with HadCRUT3 surface data
   if(ilayer .eq. 2) then !327
    do ipar=1,rcpara%parmax
     call snhteqsamp2(null,dailycrut2,rcpara%nmax,istart,istop,rcpara%snht_maxlen,rcpara%snht_increment,rcpara%miss_val,rcpara%max_miss,critical_dates,0, & 
         tsahilf,plusmean(:,ipar),minusmean(:,ipar),prmsmean(:,ipar),mrmsmean(:,ipar),philf,mhilf,rcpara%month) !in file homtests.f90 line 654
     isurf=0
     do i=1,rcpara%nmax
       if(tsahilf(i) .gt. tsa_of_mean(i,2*(ilayer-1)+ipar) .and. tsahilf(i) .gt. rcpara%break_thresh) then
!!         tsa_of_mean(i,2*(ilayer-1)+ipar)=tsahilf(i)
         isurf=isurf+1
       endif
     enddo
      !!$ call omp_set_lock(omp_lp)
       if(isurf .gt. 0) write(*,*) iname,ipar,' surface break ',isurf
      !!$ call omp_unset_lock(omp_lp)
    enddo   
  endif !327

do ipar=1,rcpara%parmax !344

  if(ilayer .eq. 2) then
!!    where(tsa_of_mean(:,5) .eq. rcpara%miss_val .and. apriori_probs .lt. 0.01) 
!!      apriori_probs=apriori_probs+0.01
!!    endwhere
  endif
   call bayes_break_simple(rcpara,apriori_probs,tsa_of_mean(:,2*(ilayer-1)+ipar),breakmean_probs(:,2*(ilayer-1)+ipar),rcpara%break_thresh,1,1,protunit) !in file rfcor.f90 line 865

   call bayes_break_simple(rcpara,apriori_probs,tsa_of_mean(:,2*(ilayer-1)+ipar+6),breakmean_probs(:,2*(ilayer-1)+ipar+6),rcpara%break_thresh,1,1,protunit) !in file rfcor.f90 line 865

 enddo !344

enddo !!ilayer  (259)

!! radiation error
call bayes_break_simple(rcpara,apriori_probs_rad,tsa_of_mean,breakmean_probs(:,5),rcpara%break_thresh_rad,5,5,protunit) !in file rfcor.f90 line 865
call bayes_break_simple(rcpara,null+rcpara%ap_prob_default,tsa_of_mean,breakmean_probs(:,6),rcpara%break_thresh_rad,6,6,protunit) !in file rfcor.f90 line 865

  if(iname .eq. 1001) then
    write(*,'(A20,7I8)') 'test in analyze',count(tm .ne. rcpara%miss_val),count(tfgm .ne. rcpara%miss_val),count(stm .ne. rcpara%miss_val),count(tfg .ne. rcpara%miss_val),count(stfgm .ne. rcpara%miss_val),count(dailycrut2 .ne. rcpara%miss_val),count(tsa_of_mean .ne. rcpara%miss_val)

  endif

return
end subroutine analyze_ts

end module correct_breaks2
