subroutine calc_profile(rcpara,iname,breakloc,plus,minus,rplus,rminus,breakprofile,f_val,kplus1,dim3)

use rfmod

implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer ip,ipar,mc,ip2,i
integer(kind=JPRM)        :: kplus1,dim3
integer(kind=JPRM) :: breakloc,pindex,mindex,iname
!!real(kind=JPRM) :: tm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax),stm(rcpara%nmax,rcpara%pmax,rcpara%parmax),stfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM),intent(in) :: plus(rcpara%nmax,rcpara%pmax,dim3),minus(rcpara%nmax,rcpara%pmax,dim3)
!!real(kind=JPRM) :: splus(rcpara%nmax,rcpara%pmax,rcpara%parmax),sminus(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM),intent(in) :: rplus(rcpara%nmax,rcpara%pmax,dim3),rminus(rcpara%nmax,rcpara%pmax,dim3)
!!real(kind=JPRM) :: plusmean(rcpara%nmax,rcpara%parmax),minusmean(rcpara%nmax,rcpara%parmax)
real(kind=JPRM) :: breakprofile(rcpara%pmax,dim3),f_val(dim3)

real            :: philf(rcpara%pmax),lplevs(rcpara%pmax),bhilf(rcpara%pmax),bhsave(rcpara%pmax,2),dmiss_val,res(21),bex(rcpara%pmax),thresh
integer(kind=JPRM)         :: ifail,goodlev,gindex(rcpara%pmax)
real            :: x(rcpara%pmax),y(rcpara%pmax),w(rcpara%pmax),work1(3*rcpara%pmax),work2(2*6)
real            :: s(6),XCAP,fit(rcpara%pmax),lfit(rcpara%pmax),nfit(rcpara%pmax,dim3),polyrms(dim3),polyrmsfit(dim3)
real, allocatable :: a(:,:)
logical         :: mask(rcpara%pmax),mask2(rcpara%pmax),nmask(rcpara%pmax,dim3)

logical svd,firstgood,log1,log2,log3
integer isx(1),idf,ldq,irank,polyorder,glevs(rcpara%parmax),polydev(dim3),devthresh
real            :: AHILF(rcpara%pmax) !!rss,b(2),se(2),cov(2*(2+1)/2),h(rcpara%pmax),qr(rcpara%pmax,3),dqr(2*2+2*2),wk(5*(2-1)+2*2),


work1=0.
isx=1
ldq=rcpara%pmax
lplevs=log(rcpara%plevs)
dmiss_val=rcpara%miss_val

thresh=rcpara%smooth_fit_thresh

polyrms=thresh+1
polydev=rcpara%pmax
devthresh=1
polyorder=kplus1
if(iname .gt. 50000 .and. iname .lt. 60000) polyorder=6
goodlev=rcpara%pmax

if(breakloc .lt. 5) then
  polyorder=goodlev-1 !! high order interpolation for adjustment at end or before gaps
  if(polyorder .gt. 8) polyorder=8
endif

!!do while(any(polyrms .gt. thresh) .and. polyorder .le. goodlev-1 .and. polyorder .le. 8)
do while(any(polydev .gt. devthresh) .and. polyorder .le. goodlev-1 .and. polyorder .le. 8)

breakprofile=rcpara%miss_val
nfit=rcpara%miss_val
polyrms=0.
polyrmsfit=0.
allocate(a(polyorder,polyorder))

do ipar=1,dim3
  goodlev=0
  f_val=0.
  bhilf=rcpara%miss_val
  bex=rcpara%miss_val
  do ip=1,rcpara%pmax
    log1=rcpara%smoothweights(ip) .gt. 0.
!!    log2= count(plus(1:breakloc,ip,ipar) .ne. rcpara%miss_val) .gt. 0.
!!    write(*,*) breakloc,rcpara%nmax,ip,ipar
!!    write(*,*) minus(breakloc:rcpara%nmax,ip,ipar)
!!    nfit=count(minus(breakloc:rcpara%nmax,ip,ipar) .ne. rcpara%miss_val)
!!    log3=count(minus(breakloc:rcpara%nmax,ip,ipar) .ne. rcpara%miss_val) .gt. 0.
    if(log1) then
!!    if(log1 .and. log2 .and. log3) then
!! with this line large data gaps are not permitted
      if(plus(breakloc,ip,ipar) .ne. rcpara%miss_val .and. minus(breakloc,ip,ipar) .ne. rcpara%miss_val ) then 
!! with this line large data gaps are permitted
!!      if(plus(breakloc,ip,ipar) .ne. 0. .or. minus(breakloc,ip,ipar) .ne. 0.) then 
        pindex=breakloc
        do while(plus(pindex,ip,ipar) .eq. rcpara%miss_val .and. pindex .gt. 1)
          pindex=pindex-1
        enddo
        mindex=breakloc
        do while(minus(mindex,ip,ipar) .eq. rcpara%miss_val .and. mindex .lt. rcpara%nmax)
          mindex=mindex+1
        enddo
        if(mindex .ne. breakloc .or. pindex .ne. breakloc) then
!!           write(*,*) 'no profile at',breakloc
           call exit(1)
        endif
        if(abs(plus(pindex,ip,ipar)-minus(mindex,ip,ipar)) .lt. 15.) then
          goodlev=goodlev+1
          bhilf(goodlev)=plus(pindex,ip,ipar)-minus(mindex,ip,ipar)
          nfit(ip,ipar)=bhilf(goodlev)
!! get plus,minus from further away only once
!!        if(mindex .ne. breakloc .or. pindex .ne. breakloc) then
!!          plus(pindex:mindex,ip,ipar)=0.
!!          minus(pindex:mindex,ip,ipar)=0.
!!        endif
          philf(goodlev)=lplevs(ip)
          w(goodlev)=rcpara%smoothweights(ip)
          gindex(goodlev)=ip
        endif
      endif
    endif
  enddo
  i=1
  if(goodlev .gt. 0) then 
!!    bhilf(goodlev)=0. ! no correction near surface
    do while (exp(philf(i)) .lt. 400. .and. i .lt. goodlev)
      i=i+1
    enddo
  endif
  if(i .gt. 3 .and. count(abs(bhilf(1:i)) .gt. 1.E-30) .gt. 4) then
  
!!    call g02ccf(goodlev,philf,bhilf,dmiss_val,dmiss_val,res,ifail)
!!    write(*,*) i,philf(1:i),'x',bhilf(1:i)
    call g02ccf(i,philf(1:i),bhilf(1:i),dmiss_val,dmiss_val,res,ifail)
!!    print*,'x',res(6),'+',res(7),' F= ',res(15)
    f_val(ipar)=res(15)
    lfit=res(7)+log(rcpara%plevs)*res(6)
!!    call g02daf('M','W',goodlev,philf,goodlev,1,isx,2,bhilf,W,rss,idf,b,se,cov,resid,h,qr,goodlev,svd,irank,dqr,0.d0,wk,ifail)
!!    lfit=b(1)+log(rcpara%plevs)*b(2)
    
    ip=1
    do while(lplevs(ip) .lt. philf(1))
      if(rcpara%smoothweights(ip) .gt. 0. .and. abs(lfit(ip)) .lt. 15.) then 
        bex(ip)=lfit(ip)
!!        write(*,*) iname,breakloc,rcpara%plevs(ip),bex(ip),' extrapolated'
      else
        bex(ip)=0.
      endif
      ip=ip+1
    enddo
    ip=ip-2
    do while(ip .gt. 1)
      bex(ip)=bex(ip)/(1.+ip)
      ip=ip-1
    enddo
    do ip=1,goodlev
      where(lplevs .eq. philf(ip)) 
        bex=bhilf(ip)
      endwhere
    enddo
    if(goodlev .gt. 8) then
      goodlev=0
      firstgood=.true.
      do ip=1,rcpara%pmax
        if(bex(ip) .ne. rcpara%miss_val .and. rcpara%smoothweights(ip) .ne. 0. ) then
          goodlev=goodlev+1
          philf(goodlev)=lplevs(ip)
          if(bex(ip) .ne. bhilf(goodlev)) then 
            w(goodlev)=rcpara%smoothweights(ip)  !!/3. ! extrapolated values have less weight
          else
            w(goodlev)=rcpara%smoothweights(ip)
            if(firstgood) then
!!              w(goodlev)=w(goodlev)*3
              firstgood=.false.
            endif
          endif
          bhilf(goodlev)=bex(ip)
          nfit(ip,ipar)=bex(ip)
        endif
      enddo
    endif

    if(goodlev .gt. 0) then
!!       bhilf(goodlev)=0. ! no correction near surface
!!       where(philf(goodlev) .eq. lplevs) nfit(:,ipar)=0.
    endif

    if(kplus1 .gt. 0 .and. polyorder .lt. goodlev) then
!! break profile with smoothed polynomial
      fit=rcpara%miss_val
      ifail=1
      call E02ADF(goodlev, polyorder, polyorder,philf(1:goodlev),bhilf(1:goodlev), W(1:goodlev), WORK1, WORK2, A, S, IFAIL)
      if(ifail .eq. 1) then
    !!$ call omp_set_lock(omp_lp)
        write(*,*) 'E02ADF failed: ',goodlev,polyorder,bhilf
    !!$ call omp_unset_lock(omp_lp)
      endif
      
      AHILF(1:polyorder)=A(polyorder,:)
      do ip=1,rcpara%pmax
        if(lplevs(ip) .gt. philf(goodlev)) then 
          breakprofile(ip,ipar)=rcpara%miss_val
        endif
!!        if(lplevs(ip) .lt. philf(1)) then
!!          if(goodlev .gt. 8) then !extrapolate only if most levels are available
!!            breakprofile(ip,ipar)=lfit(ip)
!!            breakprofile(ip,ipar)=rcpara%miss_val
!!          else
!!            breakprofile(ip,ipar)=rcpara%miss_val
!!          endif
!!        endif
        if(lplevs(ip) .ge. philf(1) .and. lplevs(ip) .le. philf(goodlev)) then 
          XCAP=((lplevs(ip)-philf(1))-(philf(goodlev)-lplevs(ip)))/(philf(goodlev)-philf(1))
          ifail=1
          call E02AEF(polyorder, AHILF, XCAP, fit(ip), IFAIL)
          if(ifail .eq. 1) then
    !!$ call omp_set_lock(omp_lp)
            write(*,*) XCAP,ip,lplevs(ip),philf(goodlev)
    !!$ call omp_unset_lock(omp_lp)
          endif
          breakprofile(ip,ipar)=fit(ip)
        endif
      enddo
      mask=fit .ne. rcpara%miss_val .and. nfit(:,ipar) .ne. rcpara%miss_val .and. rcpara%plevs .ge. 300.0
      mask2=rcpara%plevs .lt. 850 .and. fit .ne. rcpara%miss_val .and. nfit(:,ipar) .ne. rcpara%miss_val .and. rplus(pindex,:,ipar) .ne. 0 .and. &
            rminus(mindex,:,ipar) .ne. 0. .and. abs(nfit(:,ipar)-fit) .gt. 1.96*(rplus(pindex,:,ipar)+rminus(mindex,:,ipar))/2.
      mc=count(mask)
      if(mc .gt. 0) then
        polyrms(ipar)=sqrt(sum((fit-nfit(:,ipar))*(fit-nfit(:,ipar)),mask)/mc)
        polyrmsfit(ipar)=sqrt(sum((fit*fit),mask)/mc)

      endif
      polydev(ipar)=count(mask2)
   else
      if(kplus1 .eq. -1) then
!! linear break profile
        breakprofile(1:gindex(goodlev),ipar)=lfit
      else
!! unsmoothed break profile
        do ip=1,rcpara%pmax
          do ip2=1,goodlev
            if(lplevs(ip) .eq. philf(ip2)) then 
              breakprofile(ip,ipar)=bhilf(ip2)
            endif
          enddo
        enddo
      endif
    endif
  else
    breakprofile(:,ipar)=rcpara%miss_val
  endif

  if(goodlev .gt. 0) then
    glevs(ipar)=goodlev
  else
    glevs(ipar)=rcpara%pmax
  endif
  bhsave(:,ipar)=bhilf
enddo

goodlev=minval(glevs)

!!if(any(polyrms .gt. thresh) .and. polyorder .lt. goodlev) then
if(any(polydev .gt. devthresh) .and. polyorder .lt. goodlev) then
!!  write(*,'(A11,I5,I6,I2,2I2,6F7.3)') 'Bad fit at ',iname,breakloc,polyorder,polydev,polyrms,polyrmsfit,f_val
  nmask=nfit .ne. rcpara%miss_val 
  where(nmask) breakprofile=nfit

else

  mc=count(nfit .eq. rcpara%miss_val .and. breakprofile .ne. rcpara%miss_val)
  
  nmask=nfit .ne. rcpara%miss_val
  where(nmask) breakprofile=nfit

!!  write(*,'(I2, A23,I5,I6,I2,2I2,6F7.3)') mc,' fitted values used at ',iname,breakloc,polyorder,polydev,polyrms,polyrmsfit,f_val

endif


deallocate(a)
polyorder=polyorder+1

enddo !! polyrms .gt. thresh

if(any(breakprofile .ne. rcpara%miss_val .and. abs(breakprofile) .gt. 15.)) then
    !!$ call omp_set_lock(omp_lp)
  write(*,*) iname,': break too large',breakloc
  write(*,*) 'breakprofile',breakprofile
  write(*,*) 'bhilf',bhsave
  write(*,*) 'abandoning break profile',bhsave
    !!$ call omp_unset_lock(omp_lp)
  breakprofile=rcpara%miss_val
!!  call exit(1)
endif

return
end subroutine calc_profile

