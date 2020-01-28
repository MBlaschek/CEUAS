module rfcorio

use rfmod
!!use homtests
!!use homtestsamp

!! switch on/off additional debugging information
logical ldebug
!! verbose is used to control output into unit protunit
!! verbose=0: only error messages
!! verbose=1: output for later statistics turned on
!! verbose=2: all output into protunit turned on
!! default is 1 
integer verbose

contains

subroutine read_firstiter(rcpara,iunit,filename,bindex,bi,ttmc,rios) !


implicit none

type(rasocor_namelist),intent(in) :: rcpara

character*(*) filename
integer iunit,err,i,rios,bi,ipmax,iparmax,ip,ipar,bindex(rcpara%brmax)
real(kind=JPRM) :: hilfcorr(rcpara%brmax,rcpara%pmax,rcpara%parmax),ttmc(rcpara%nmax,rcpara%pmax,rcpara%parmax)

    rios=0

     open(iunit,file=filename,form='unformatted',status='old',action='read',iostat=rios)

!    if(ios .eq. 0 .or. first) then 
     if(rios .eq. 0) then 
      err=0

      read(iunit) bi,ipmax,iparmax
if (bi .gt. 50 .or. bi .lt. 1) then  
        !!$ call omp_set_lock(omp_lp)
        print*,'error: ',filename,bi,ipmax,iparmax
        !!$ call omp_unset_lock(omp_lp)
endif

      if(bi .gt. rcpara%brmax .or. ipmax .ne. rcpara%pmax .or. iparmax .ne. rcpara%parmax) then
        !!$ call omp_set_lock(omp_lp)
        print*,'error: ',filename,bi,ipmax,iparmax
        !!$ call omp_unset_lock(omp_lp)
      endif
      read(iunit) bindex(1:bi)
      read(iunit) hilfcorr(1:bi,:,:)
      
      bindex(bi+1)=rcpara%nmax

      close(iunit)


      ttmc=0.
      do ipar=1,rcpara%parmax
        do ip=1,rcpara%pmax
          do i=1,bi-1
            ttmc(bindex(i):bindex(i+1),ip,ipar)=hilfcorr(i,ip,ipar)
            
         enddo
        enddo
      enddo
     else
       ttmc=0.
     endif
end subroutine read_firstiter

subroutine read_oldrich_firstiter(rcpara,iunit,filename,bindex,bi,tbindex,tbi,ttmc,rios) !


implicit none

type(rasocor_namelist),intent(in) :: rcpara

character*(*) filename
integer iunit,err,i,rios,bi,tbi,ipmax,iparmax,ip,ipar,bindex(rcpara%brmax),tbindex(rcpara%brmax)
real(kind=JPRM) :: hilfcorr(rcpara%brmax,rcpara%pmax,rcpara%parmax),ttmc(rcpara%nmax,rcpara%pmax,rcpara%parmax)

    rios=0
    hilfcorr=0
    bindex=0

     open(iunit,file=filename,form='unformatted',status='old',action='read',iostat=rios)

!    if(ios .eq. 0 .or. first) then 
     if(rios .eq. 0) then 
      err=0

      read(iunit) bi,ipmax,iparmax
if (bi .gt. 50 .or. bi .lt. 1) then  
        !!$ call omp_set_lock(omp_lp)
        print*,'error: ',filename,bi,ipmax,iparmax
        !!$ call omp_unset_lock(omp_lp)
endif

      if(bi .gt. rcpara%brmax .or. ipmax .ne. rcpara%pmax .or. iparmax .ne. rcpara%parmax) then
        !!$ call omp_set_lock(omp_lp)
        print*,'error: ',filename,bi,ipmax,iparmax
        !!$ call omp_unset_lock(omp_lp)
      endif
      read(iunit) bindex(1:bi)
      read(iunit) hilfcorr(1:bi,:,:)
      
      if(tbindex(tbi) .ne. rcpara%nmax-2) then
         bi=bi-1
      do ipar=1,rcpara%parmax
        do ip=1,rcpara%pmax
          do i=1,bi
            hilfcorr(i,ip,ipar)=hilfcorr(i,ip,ipar)-hilfcorr(bi,ip,ipar)           
         enddo
        enddo
      enddo
      endif
        
      bindex(bi+1)=rcpara%nmax

      close(iunit)


      ttmc=0.
      do ipar=1,rcpara%parmax
        do ip=1,rcpara%pmax
          do i=1,bi-1
            ttmc(bindex(i):bindex(i+1),ip,ipar)=hilfcorr(i,ip,ipar)
            
         enddo
        enddo
      enddo
     else
       ttmc=0.
     endif
end subroutine read_oldrich_firstiter

subroutine read_IRI(iunit,filename,rcpara,type,plevs,profile,err) !

implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer iunit,err,i,j,k,l
integer statanf,statend,ianf,iend,tanf,tend,panf,pend
integer :: statnr,wmonrs(rcpara%statmax),plevs(rcpara%pmax)

real(kind=JPRM) :: profile(rcpara%pmax,rcpara%parmax),prof(rcpara%pmax)

character*(*) filename
character*106 zeile
character*(*) type
character*2 cpmax


  open(iunit,file=trim(filename),form='formatted',status='old',action='read',err=120)
  err=0
  do while(.true.)
    read(iunit,'(A106)',err=120,end=100) zeile
   if(zeile(1:1) .ne. '#') then
!!$ call omp_set_lock(omp_lp)
      if(index(zeile,trim(type)) .ne. 0) then 
        read(iunit,'(A106)',err=120,end=100) zeile
        write(cpmax,'(I2.2)') rcpara%pmax
        read(zeile,'('//cpmax//'I5)') plevs 
        read(iunit,'(A106)',err=120,end=100) zeile
        read(zeile,'('//cpmax//'F5.1)') prof
        profile(:,1)=prof 
        read(iunit,'(A106)',err=120,end=100) zeile
        read(zeile,'('//cpmax//'F5.1)') prof
        profile(:,2)=prof 
      endif
    endif
  enddo

100  close(iunit)
return

120 continue
!!$ call omp_set_lock(omp_lp)
  print*,'could not read ',filename,err
!!$ call omp_unset_lock(omp_lp)
call abort

return
end subroutine read_IRI

subroutine savespagarr(rcpara,iunit,filename,spagarr,tbindex,tbi,xdrarr,wmonrs,dists,index,istatmax) !

implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer iunit,tbi,ios,i,l,k,j,m,n,istatmax
logical lfound
character*80 filename

real(kind=JPRM) :: spagarr(rcpara%brmax,rcpara%pmax,rcpara%parmax,rcpara%cachemax),xdrarr(rcpara%brmax,rcpara%pmax,rcpara%parmax),disthilf(rcpara%statmax)
real            :: dists(rcpara%statmax) 
real(kind=JPRM) :: spaghilf(tbi,rcpara%pmax,rcpara%parmax,rcpara%cachemax),xdrarrhilf(tbi,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: shilf(tbi*rcpara%pmax*rcpara%parmax*200)
integer         :: idx(tbi*rcpara%pmax*rcpara%parmax*200)
integer         :: wmonrs(rcpara%statmax),wmonrhilf(tbi,rcpara%statmax),index(rcpara%cachemax),tbindex(rcpara%brmax)

open(iunit,file=filename,form='unformatted',iostat=ios)
if(ios .ne. 0) then
!!!$ call omp_set_lock(omp_lp)
  write(*,*) 'versager',ios
!!!$ call omp_unset_lock(omp_lp)
else

l=0
do i=1,istatmax
  lfound=.false.
tb: do k=1,rcpara%parmax
      do m=1,rcpara%pmax
       do j=1,tbi
        if(spagarr(j,m,k,i) .ne. rcpara%miss_val) then
          lfound=.true.
          exit tb
        endif
       enddo
      enddo
    enddo tb
!if(any(spagarr(1:tbi,:,:,i) .ne. rcpara%miss_val)) then
 if(lfound) then
    l=l+1
    spaghilf(:,:,:,l)=spagarr(1:tbi,:,:,i)
    wmonrhilf(:,l)=wmonrs(index(i))
    do k=1,tbi
      if(.not. any(spaghilf(k,:,:,l) .ne. rcpara%miss_val)) wmonrhilf(k,l)=-1
    enddo
  endif
enddo

xdrarrhilf=xdrarr(1:tbi,:,:)
    
 write(iunit) tbi,l,rcpara%pmax,rcpara%parmax

 if(tbi .gt. 0 .and. l .gt. 0) then
   write(iunit) tbindex(1:tbi)
   write(iunit) wmonrhilf(:,1:l)
   n=0
   do i=1,tbi
     do j=1,rcpara%pmax
       do k=1,rcpara%parmax
         do m=1,l
           if(spaghilf(i,j,k,m) .ne. rcpara%miss_val) then
              n=n+1
              shilf(n)=spaghilf(i,j,k,m)
              idx(n)=(m-1)*rcpara%parmax*rcpara%pmax*tbi+(k-1)*rcpara%pmax*tbi+(j-1)*tbi+i-1
           endif
         enddo
       enddo
     enddo
   enddo
!   write(iunit) spaghilf(:,:,:,1:l)
   write(iunit) xdrarrhilf   
   write(iunit) n
   write(iunit) idx(1:n)
   write(iunit) shilf(1:n)
 endif
 

 close(iunit)
endif

return
end subroutine savespagarr

! reads sorted list of intervals with bad data
subroutine read_bad_intervals(filename,rcpara,bad_intervals)
implicit none

type(rasocor_namelist),intent(in) :: rcpara
integer :: iunit,l,dstart,dstop,iplev,ic,statnr,ipar
character*(*) filename
character*80 zeile
integer,allocatable:: bad_intervals(:,:)
real(kind=JPRM) :: ri,raw,ra

iunit=55
open(iunit,file=filename,err=100)

l=0
do while(.true.)
   read(iunit,'(A6)',end=97,err=97) zeile
   ic=iachar(zeile(1:1))
   if (ic>47 .and. ic<58) l=l+1 
enddo
97 continue
if (l==0) then
   allocate(bad_intervals(5,1))
   return
else
   allocate(bad_intervals(5,l))
   rewind(iunit)
   l=0
   do while(.true.)
     read(iunit,'(A80)',end=98,err=98) zeile
     ic=iachar(zeile(1:1))
     if (ic>47 .and. ic<58) then
       l=l+1 
       read(zeile,'(I6,2I9,I3,I3)',end=98,err=98) statnr,dstart,dstop,iplev,ipar
       bad_intervals(1,l)=statnr
       bad_intervals(2,l)=toindex(dstart,rcpara)
       bad_intervals(3,l)=toindex(dstop,rcpara)
       bad_intervals(4,l)=iplev
       bad_intervals(5,l)=ipar
     endif
     
   enddo

endif

98 continue
write(*,*) 'read ',l,' bad intervals'
close(iunit)
return
100 continue
allocate(bad_intervals(4,1))
write(*,*) 'bad_intervals not found'
return
end subroutine read_bad_intervals

subroutine read_alarm(iunit,filename,alarm)

implicit none

integer :: iunit,l,loc(3),statnr
integer :: alarm(10000,3)
character*(*) filename
real(kind=JPRM) :: ri,raw,ra

open(iunit,file=filename,err=100)

l=0
do while(.true.)
   l=l+1 
 read(iunit,'(I5,21x,3F7.1,3I6)',end=98,err=98) statnr,ri,raw,ra,loc
 alarm(l,1)=statnr
 alarm(l,2)=loc(1)
! if(ri .gt. 2*ra .and. ri .gt. 1.1*raw) then 
!   alarm(l,3)=0
! else
   alarm(l,3)=1
! endif
enddo

98 continue
write(*,*) 'read ',l,' alarm lines'
close(iunit)
100 continue
return
end subroutine read_alarm


subroutine read_hadCRUT3(iunit,filename,crut2,rcpara,absfile) !


implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer iunit,err,jahr,monat,i,j,k

real(kind=JPRM) :: crut2(72,36,rcpara%mmax),tmonat(72)

character*(*) filename
character*(*)  absfile
character*200 zeile


crut2=0.
  open(iunit,file=trim(filename),form='formatted',status='old',action='read',err=120)
  err=0
  do while(.true.)
    read(iunit,'(A200)',err=120,end=100) zeile
    read(zeile(1:12),'(2I6)') jahr,monat
    do i=1,36
      read(iunit,'(72(E10.3,1x))') tmonat
      if(jahr .ge. rcpara%startdate/10000 .and. (jahr-rcpara%startdate/10000)*12+monat .lt. rcpara%mmax) crut2(:,i,(jahr-rcpara%startdate/10000)*12+monat)=tmonat
    enddo
  enddo

100 close(iunit)

  open(iunit,file=trim(absfile),form='formatted',status='old',action='read',err=120)
    read(iunit,'(A200)',err=120,end=200) zeile
    read(zeile(1:12),'(2I6)') jahr,monat
    do j=1,36
      read(iunit,'(72(E10.3,1x))') tmonat
      do k=1,rcpara%mmax
         do i=1,72
           if(crut2(i,j,k) .gt. -1.e29) then
             crut2(i,j,k)=crut2(i,j,k)+273.15+tmonat(i)
           endif
         enddo        
      enddo
    enddo


200 close(iunit)
 return

120 continue
!!$ call omp_set_lock(omp_lp)
  print*,'could not read ',filename,err
!!$ call omp_unset_lock(omp_lp)
call abort

end subroutine read_hadCRUT3


subroutine read_era40_meta(iunit,filename,ni,era40meta,year,month,day,err) !

implicit none

integer ni,iunit,err,iyear,imon,iday,jumpsize
integer :: era40meta(ni),year(ni),month(ni),day(ni)

character*(*) filename
character*77 zeile

  iunit=iunit

  open(iunit,file=filename,form='formatted',status='old',action='read',err=120)
  err=0
  era40meta=0
  do while(.true.)
    read(iunit,'(A77)',err=120,end=100) zeile
    if(zeile(1:1) .eq. '1' .or. zeile(1:1) .eq. '2') then

      read(zeile,'(I4,2I2.2,I3)') iyear,imon,iday,jumpsize
      where(iyear .eq. year .and. imon .eq. month .and. iday .eq. day) era40meta=jumpsize

    endif
  enddo

100  close(iunit)
return

120 print*,'could not read ',filename
call abort

return
end subroutine read_era40_meta

subroutine read_raso_blacklist(iunit,filename,rcpara,statnr,tm,tfgm,err) !

implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer iunit,err,i,j,k,l
integer statanf,statend,ianf,iend,tanf,tend,panf,pend
integer :: statnr,wmonrs(rcpara%statmax)

real(kind=JPRM) :: tm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax)

character*(*) filename
character*106 zeile
character*60 comment


  open(iunit,file=trim(filename),form='formatted',status='old',action='read',err=120)
  err=0
  do while(.true.)
    read(iunit,'(A106)',err=120,end=100) zeile
   if(zeile(1:1) .ne. '#') then
!!$ call omp_set_lock(omp_lp)
      read(zeile,'(2(I5,1X),2(I8,1X),2(I2,1X),2(I4,1X),A60)') statanf,statend,ianf,iend,tanf,tend,panf,pend,comment
!!$ call omp_unset_lock(omp_lp)
!!      do i=1,rcpara%statmax
        if(statanf .ne. 0 .and. statnr .ge. statanf .and. statnr .le. statend) then
            
          do j=1,rcpara%nmax
            if(rcpara%year(j)*10000+rcpara%month(j)*100+rcpara%day(j) .ge. ianf .and. rcpara%year(j)*10000+rcpara%month(j)*100+rcpara%day(j) .le. iend) then
              do k=1,rcpara%parmax
                if((k-1)*12 .ge. tanf .and. (k-1)*12 .le. tend) then
                   do l=1,rcpara%pmax 
                      if(rcpara%plevs(l) .ge. panf .and. rcpara%plevs(l) .lt. pend) then
                         tm(j,l,k)=rcpara%miss_val
                         tfgm(j,l,k)=rcpara%miss_val
                      endif
                   enddo
                endif
              enddo
            endif
            
          enddo
        endif
!!      enddo
    endif
  enddo

100  close(iunit)
return

120 continue
!!$ call omp_set_lock(omp_lp)
  print*,'could not read ',filename,err
!!$ call omp_unset_lock(omp_lp)
call abort

return
end subroutine read_raso_blacklist

subroutine erastations(wmonrs,wmolats,wmolons,statmax,wmostats,filename) !

implicit none

character filename*80
character clat*9,clon*8,cstat*35,cdum*80,cdum2*80
integer iunit,istat,obs,statmax,i,l,wmostats,height,icountry,idum
integer wmonrs(statmax)
real(kind=JPRM) wmolats(statmax),wmolons(statmax),lat,lon
!!character*20 wmonames(statmax)
logical lmissing
character*200 zeile

iunit=20
open(iunit,file=filename,status='old',err=30)

do while (.true.)
  READ(iunit,'(A200)', END=20) zeile
  print*,zeile
  if(zeile(10:10) == '#') then
      READ(zeile,'(I5,I4,1x,I5,F8.2,F8.2,I8)') l,Idum,istat,lat,lon,height
  else
      READ(zeile,'(I5,I4,I7,F7.2,F8.2)') l,Idum,istat,lat,lon
  endif

  if(lat .ne. 0 .or. lon .ne. 0) then 
  wmonrs(l)=istat
  wmolats(l)=lat
  wmolons(l)=lon

  endif
  
enddo

20 close(iunit)
wmostats=l

return

30 continue
             !!$ call omp_set_lock(omp_lp)
write(*,*) 'error opening '//filename
             !!$ call omp_unset_lock(omp_lp)
call abort

end subroutine erastations

subroutine read_country(iunit,filename,nmem,ngroup,stgroups,stgroupnames,err)

implicit none

integer nmem,ngroup,imem,igroup,err,iunit
character*(*) filename
integer(kind=JPRM) :: stgroups(nmem,ngroup)
character*(*)  ::  stgroupnames(ngroup)

  iunit=iunit

  open(iunit,file=filename,form='unformatted',status='old',action='read',err=120)
  err=0
  read(iunit) imem,igroup
  read(iunit) stgroups
  read(iunit) stgroupnames

  close(iunit)

return

120 print*,'could not open file ',filename
  err=1

  return

end subroutine read_country

subroutine read_sonde(iunit,filename,rcpara,tm,tanm,tfgm,tbcm,solarangles,rtype,err,notbc) !

implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer ini,il,ipmax,iparmax,i,iunit,err,imax,ios,qc_count,bl_count,iflag,ifak,ip,ipar
integer,allocatable :: goodindex(:)
character*(*) filename
real(kind=JPRM) :: tm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tanm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tbcm(rcpara%nmax,rcpara%pmax,rcpara%parmax),solarangles(rcpara%nmax,rcpara%parmax)
real(kind=JPRM) :: tflags(rcpara%nmax,rcpara%pmax,rcpara%parmax)
!!ttypes(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: hilf2(rcpara%pmax,rcpara%parmax)
real(kind=JPRM),allocatable :: hilf(:,:,:),hilf2d(:,:)
logical :: mask(rcpara%pmax,rcpara%parmax)
integer :: rtype(rcpara%nmax),rrcor(rcpara%nmax),val,mr,nancount
logical, optional :: notbc


  tm=rcpara%miss_val
  tanm=rcpara%miss_val
  tfgm=rcpara%miss_val
  tbcm=rcpara%miss_val
!!  solarangles=rcpara%miss_val
  tflags=rcpara%miss_val
!!  ttypes=rcpara%miss_val

  iflag=0

  iunit=iunit
  open(iunit,file=filename,form='unformatted',status='old',action='read',err=120)
  err=0
  read(iunit) ini,il,ipmax,iparmax

  allocate(hilf(il,ipmax,iparmax))
  allocate(hilf2d(il,iparmax))
  allocate(goodindex(il))

  read(iunit) goodindex

  imax=0
  do i=1,il
    if(goodindex(i) .le. rcpara%nmax) imax=imax+1
  enddo

  read(iunit,end=200) hilf
  nancount=count(hilf .ne. hilf)
  if(nancount .gt. 0) then
              !!$ call omp_set_lock(omp_lp)
   write(*,*) nancount,' NaNs found in ',filename
             !!$ call omp_unset_lock(omp_lp)
    where(hilf .ne. hilf) hilf=rcpara%miss_val
  endif
  do ipar=1,rcpara%parmax
  do ip=1,rcpara%pmax
  do i=1,imax
    tm(goodindex(i),ip,ipar)=hilf(i,ip,ipar)
  enddo
  enddo
  enddo

  read(iunit) hilf
  where(hilf .eq. -rcpara%miss_val)
   hilf=rcpara%miss_val
  endwhere
  do ipar=1,rcpara%parmax
  do ip=1,rcpara%pmax
  do i=1,imax
    tanm(goodindex(i),ip,ipar)=hilf(i,ip,ipar)
  enddo
  enddo
  enddo

  read(iunit) hilf
  where(hilf .eq. -rcpara%miss_val)
   hilf=rcpara%miss_val
  endwhere
  do ipar=1,rcpara%parmax
  do ip=1,rcpara%pmax
  do i=1,imax
    tfgm(goodindex(i),ip,ipar)=hilf(i,ip,ipar)
  enddo
  enddo
  enddo

if(solarangles(1,1) .ne. 0) then
  read(iunit) hilf2d
!  do i=1,imax
!!    solarangles(goodindex(i),:)=hilf(i,1,:)
!  enddo
endif

  read(iunit) hilf
  where(hilf .eq. -rcpara%miss_val)
   hilf=rcpara%miss_val
  endwhere
  do ipar=1,rcpara%parmax
  do ip=1,rcpara%pmax
  do i=1,imax
    tbcm(goodindex(i),ip,ipar)=hilf(i,ip,ipar)
  enddo
  enddo
  enddo

  if(any(tfgm .eq. -rcpara%miss_val)) then
    !!$ call omp_set_lock(omp_lp)
    write(*,*) 'Achtung'
    !!$ call omp_unset_lock(omp_lp)
  endif

  if(.not. present(notbc)) then 
    do i=1,imax
      tbcm(goodindex(i),:,:)=hilf(i,:,:)
      mask=tbcm(goodindex(i),:,:) .ne. rcpara%miss_val .and. tfgm(goodindex(i),:,:) .ne. rcpara%miss_val

!! reverse effect of operational bias correction

      where(mask)
        tfgm(goodindex(i),:,:)=tfgm(goodindex(i),:,:)+tbcm(goodindex(i),:,:)
        tanm(goodindex(i),:,:)=tanm(goodindex(i),:,:)+tbcm(goodindex(i),:,:)
      endwhere
    enddo
  else
    if(.not. notbc) then
    do i=1,imax
      tbcm(goodindex(i),:,:)=hilf(i,:,:)
      mask=tbcm(goodindex(i),:,:) .ne. rcpara%miss_val .and. tfgm(goodindex(i),:,:) .ne. rcpara%miss_val

!! reverse effect of operational bias correction

      where(mask)
        tfgm(goodindex(i),:,:)=tfgm(goodindex(i),:,:)+tbcm(goodindex(i),:,:)
        tanm(goodindex(i),:,:)=tanm(goodindex(i),:,:)+tbcm(goodindex(i),:,:)
      endwhere
    enddo
   endif
  endif

!! read qality control flags
!! 100000 = final flag
!! 10000  = first guess flag
!! 1000   = departure flag
!! 100    = analysis quality control flag
!! 10     = blacklist flag
!! 1      = analysis TEMP datum event
  read(iunit,iostat=ios) hilf
  if( ios .eq. 0) then
    qc_count=0
    bl_count=0
    do i=1,imax
      mask=.false.
      if(rcpara%qc .eq. 'T') then
        mask=mask .or. (tfgm(goodindex(i),:,:) .ne. rcpara%miss_val .and. abs(tfgm(goodindex(i),:,:)) .gt. rcpara%qcthresh)
      else
        tflags(goodindex(i),:,:)=hilf(i,:,:)
        ifak=100000
        do iflag=1,5
          if(iflag .eq. 5 .and. rcpara%qc .eq. 'B') then
            mask=mask .or. (tfgm(goodindex(i),:,:) .ne. rcpara%miss_val .and. abs(tfgm(goodindex(i),:,:)) .gt. rcpara%qcthresh)
          else
            hilf2=nint(hilf(i,:,:))/ifak
            mask=mask .or. hilf2 .gt. 2
            hilf(i,:,:)=hilf(i,:,:)-hilf2*ifak
            ifak=ifak/10
          endif
        enddo
      endif
      qc_count=qc_count+count(mask)
      where(mask)
       tm(goodindex(i),:,:)=rcpara%miss_val
       tbcm(goodindex(i),:,:)=rcpara%miss_val
       tfgm(goodindex(i),:,:)=rcpara%miss_val
       tanm(goodindex(i),:,:)=rcpara%miss_val
      endwhere
      if(any(abs(tfgm(goodindex(i),:,:)) .gt. 2*rcpara%qcthresh .and. tfgm(goodindex(i),:,:) .ne. rcpara%miss_val)) then
    !!$ call omp_set_lock(omp_lp)
        write(*,*) 'suspect values found',tfgm(goodindex(i),:,:)
    !!$ call omp_unset_lock(omp_lp)
      endif
      if(iflag .eq. 5 .and. rcpara%qc .ne. 'T') then
        hilf2=nint(hilf(i,:,:))/ifak
        mask=hilf2 .gt. 2
        bl_count=bl_count+count(mask)
      endif
    enddo
  else
    return
  endif
!!  write(*,'(A17,I6,A17,I6,A17)') 'ERA40 data: ', qc_count,' rejected by QC; ',bl_count,' blacklisted'
 

!! read radiosonde type information

  rtype=nint(rcpara%miss_val)
  rrcor=nint(rcpara%miss_val)
  read(iunit,iostat=ios) hilf
  if(ios .eq. 0) then
   do ipar=1,rcpara%parmax
   do ip=1,rcpara%pmax
    do i=1,imax
!!      ttypes(goodindex(i),ip,ipar)=hilf(i,ip,ipar)
      if(hilf(i,ip,ipar) .gt. rrcor(goodindex(i))) then
         rrcor(goodindex(i))=hilf(i,ip,ipar)
      endif
    enddo
   enddo
   enddo
   do i=1,imax
     rtype(goodindex(i))=mod(rrcor(goodindex(i))+1000,1000)
     rrcor(goodindex(i))=rrcor(goodindex(i))/1000
   enddo    
  endif

!!  call rtypeinfo(rtype,rcpara,mr)

  close(iunit)

  deallocate(goodindex,hilf2d,hilf)
  return
200 continue
 print*,'could open but could not read ',filename
 return
120 continue
   !!$ call omp_set_lock(omp_lp)
   print*,'could not open file ',filename
   !!$ call omp_unset_lock(omp_lp)
  err=1

  return

end subroutine read_sonde

! calculates adjustments between ERA-Interim and ERA-40 background departures
subroutine eiminuse40(rcpara,lat,olon,ominuse40,adjust)

implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer ini,il,ipmax,iparmax,i,iunit,err,imax,ip,ipar,ios,iflag,bl_count,qc_count,ifak,j,isave,jsave,izy,ic,min,max,imon
real(kind=JPRM) :: lat,lon,olon

real(kind=JPRM),intent(in) :: ominuse40(rcpara%ni,rcpara%nj,rcpara%pmax,rcpara%parmax,12)
real(kind=JPRM) :: glon(rcpara%ni+1),glat(rcpara%nj),w(2,2),glj
real(kind=JPRM) :: adjust(rcpara%pmax,rcpara%parmax,12),wsum


if (rcpara%fgdepname .ne. 'fg_dep') then
   adjust=0.
   return
endif

do i=1,rcpara%nj
 glat(i)=-90.+(i-1)*1.
enddo
do i=1,rcpara%ni+1
 glon(i)=(i-1)*1.
enddo

j=0
do while(glat(j+1) .le. lat)
  j=j+1
enddo
if(j .eq. 0) then
  lat=0
  glj=0
else
  glj=glat(j)
endif
jsave=j

i=1
lon=olon
if(lon .lt. 0.) lon=lon+360.
do while(glon(i) .le. lon)
  i=i+1
enddo
i=i-1
isave=i
izy=i+1
if(i .eq. rcpara%ni) izy=1
!!write(*,*) lon,i,izy,glon(i)
 
w(1,1)=(glon(izy)-lon)/(glon(izy)-glon(i))*(glat(j+1)-lat)/(glat(j+1)-glj)
w(1,2)=(glon(izy)-lon)/(glon(izy)-glon(i))*(lat-glj)/(glat(j+1)-glj)
w(2,1)=(lon-glon(i))/(glon(izy)-glon(i))*(glat(j+1)-lat)/(glat(j+1)-glj)
w(2,2)=(lon-glon(i))/(glon(izy)-glon(i))*(lat-glj)/(glat(j+1)-glj)

!write(*,*) 'Sum:',sum(w)

do imon=1,12
  do ipar=1,rcpara%parmax
    do ip=1,rcpara%pmax
      wsum=0
      adjust(ip,ipar,imon)=0
      if(ominuse40(isave,jsave,ip,ipar,imon) /= rcpara%miss_val) then
        adjust(ip,ipar,imon)=adjust(ip,ipar,imon)+w(1,1)*ominuse40(isave,jsave,ip,ipar,imon)
        wsum=wsum+w(1,1)
      endif
      if(   ominuse40(isave,jsave+1,ip,ipar,imon) /= rcpara%miss_val) then
adjust(ip,ipar,imon)=adjust(ip,ipar,imon)+w(1,2)*ominuse40(isave,jsave+1,ip,ipar,imon)
        wsum=wsum+w(1,2)
      endif
      if(ominuse40(izy,jsave,ip,ipar,imon) /= rcpara%miss_val) then
adjust(ip,ipar,imon)=adjust(ip,ipar,imon)+w(2,1)*ominuse40(izy,jsave,ip,ipar,imon)
        wsum=wsum+w(2,1)
      endif
      if(ominuse40(izy,jsave+1,ip,ipar,imon) /= rcpara%miss_val) then
adjust(ip,ipar,imon)=adjust(ip,ipar,imon)+w(2,2)*ominuse40(izy,jsave+1,ip,ipar,imon)
        wsum=wsum+w(2,2)
      endif
      if(wsum >0) then
         adjust(ip,ipar,imon)=adjust(ip,ipar,imon)/wsum
      else
        adjust(ip,ipar,imon)=rcpara%miss_val   
      endif    
    enddo
  enddo
enddo

return
end subroutine eiminuse40
  
!! read_sonde_oper reads operational feedback data and
!! performs the "merge" with ERA-40 feedbackdata using
!! a field of operational bg minus ERA40 bg for
!! adjustment of the temperatures between the two.
!! The merge is performed such that the ERA-40 feedback
!! is adjusted to the operational feedback.
!! 
!! Leo Haimberger
!! 20050201
subroutine read_sonde_oper(iunit,filename,rcpara,lat,olon,tm,tanm,tfgm,tbcm,solarangles,ominuse40,adjust,rtype,err,notbc) !

implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer ini,il,ipmax,iparmax,i,iunit,err,imax,ip,ipar,ios,iflag,bl_count,qc_count,ifak,j,isave,jsave,izy,ic,min,max
integer counts(rcpara%pmax,rcpara%parmax)
integer,allocatable :: goodindex(:)
character*(*) filename
real(kind=JPRM) :: tm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tanm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tbcm(rcpara%nmax,rcpara%pmax,rcpara%parmax),solarangles(rcpara%nmax,rcpara%parmax)
real(kind=JPRM) :: lat,lon,olon
real(kind=JPRM), allocatable :: overlap(:,:),hilf(:,:,:),hilf2d(:,:)
real(kind=JPRM) :: rms(rcpara%nmax),diff(rcpara%pmax,rcpara%parmax)
logical :: mask(rcpara%pmax,rcpara%parmax),omask(365)
real(kind=JPRM) :: tflags(rcpara%nmax,rcpara%pmax,rcpara%parmax),ttypes(rcpara%nmax,rcpara%pmax,rcpara%parmax),hilf2(rcpara%pmax,rcpara%parmax)

real(kind=JPRM),intent(in) :: ominuse40(rcpara%ni,rcpara%nj,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: glon(rcpara%ni+1),glat(rcpara%nj),w(2,2),glj
real(kind=JPRM) :: adjust(rcpara%pmax,rcpara%parmax)
integer :: rtype(rcpara%nmax),rrcor(rcpara%nmax),val,mr,nancount
logical,optional :: notbc

  iflag=0
  iunit=iunit
  open(iunit,file=filename,form='unformatted',status='old',action='read',iostat=err)

do i=1,rcpara%nj
 glat(i)=-90.+(i-1)*1.
enddo
do i=1,rcpara%ni+1
 glon(i)=(i-1)*1.
enddo

j=0
do while(glat(j+1) .le. lat)
  j=j+1
enddo
if(j .eq. 0) then
  lat=0
  glj=0
else
  glj=glat(j)
endif
jsave=j

i=1
lon=olon
if(lon .lt. 0.) lon=lon+360.
do while(glon(i) .le. lon)
  i=i+1
enddo
i=i-1
isave=i
izy=i+1
if(i .eq. rcpara%ni) izy=1
!!write(*,*) lon,i,izy,glon(i)
 
w(1,1)=(glon(izy)-lon)/(glon(izy)-glon(i))*(glat(j+1)-lat)/(glat(j+1)-glj)
w(1,2)=(glon(izy)-lon)/(glon(izy)-glon(i))*(lat-glj)/(glat(j+1)-glj)
w(2,1)=(lon-glon(i))/(glon(izy)-glon(i))*(glat(j+1)-lat)/(glat(j+1)-glj)
w(2,2)=(lon-glon(i))/(glon(izy)-glon(i))*(lat-glj)/(glat(j+1)-glj)

  do ipar=1,rcpara%parmax
    do ip=1,rcpara%pmax
      adjust(ip,ipar)=w(1,1)*ominuse40(isave,jsave,ip,ipar)+w(1,2)*ominuse40(isave,jsave+1,ip,ipar)+w(2,1)*ominuse40(izy,jsave,ip,ipar)+w(2,2)*ominuse40(izy,jsave+1,ip,ipar)
!!      where(tfgm(1:goodindex(1)-1,ip,ipar) .ne. rcpara%miss_val)
!!        tfgm(1:goodindex(1)-1,ip,ipar)=tfgm(1:goodindex(1)-1,ip,ipar)+adjust(ip,ipar)
!!      endwhere
    enddo
  enddo

  if(err .ne. 0) then
!!    write(*,*) filename,' not found, calculated bg adjustment only'
    return
  endif
  
  read(iunit) ini,il,ipmax,iparmax
  allocate(goodindex(il),overlap(rcpara%pmax,rcpara%parmax),hilf2d(il,rcpara%parmax),hilf(il,rcpara%pmax,rcpara%parmax))

  read(iunit) goodindex

  overlap=rcpara%miss_val

  imax=0
  do i=1,il
    if(goodindex(i) .le. rcpara%nmax) imax=imax+1
  enddo

  read(iunit) hilf
  nancount=count(hilf .ne. hilf)
  if(nancount .gt. 0) then
    !!$ call omp_set_lock(omp_lp)
    write(*,*) nancount,' NaNs found in ',filename
    !!$ call omp_unset_lock(omp_lp)
    where(hilf .ne. hilf) hilf=rcpara%miss_val
  endif
  do i=1,imax
    tm(goodindex(i),:,:)=hilf(i,:,:)
  enddo

  read(iunit) hilf
  do i=1,imax
    tanm(goodindex(i),:,:)=hilf(i,:,:)
  enddo


!! no ERA data after 2000
  tfgm(16072:rcpara%nmax,:,:)=rcpara%miss_val
  read(iunit) hilf

!! concatenate operational data   
  do i=1,imax
   tfgm(goodindex(i),:,:)=hilf(i,:,:)
  enddo

if(solarangles(1,1) .ne. 0) then
  read(iunit) hilf2d
  do i=1,imax
!!    solarangles(goodindex(i),:)=hilf(i,1,:)
  enddo
endif

  tbcm(16072:rcpara%nmax,:,:)=rcpara%miss_val
  read(iunit) hilf
  if(.not. present(notbc)) then 
    do i=1,imax
      tbcm(goodindex(i),:,:)=hilf(i,:,:)
      mask=tbcm(goodindex(i),:,:) .ne. rcpara%miss_val .and. tfgm(goodindex(i),:,:) .ne. rcpara%miss_val

!! reverse effect of operational bias correction

      where(mask)
        tfgm(goodindex(i),:,:)=tfgm(goodindex(i),:,:)+tbcm(goodindex(i),:,:)
        tanm(goodindex(i),:,:)=tanm(goodindex(i),:,:)+tbcm(goodindex(i),:,:)
      endwhere
    enddo
  else
    if(.not. notbc) then
    do i=1,imax
      tbcm(goodindex(i),:,:)=hilf(i,:,:)
      mask=tbcm(goodindex(i),:,:) .ne. rcpara%miss_val .and. tfgm(goodindex(i),:,:) .ne. rcpara%miss_val

!! reverse effect of operational bias correction

      where(mask)
        tfgm(goodindex(i),:,:)=tfgm(goodindex(i),:,:)+tbcm(goodindex(i),:,:)
        tanm(goodindex(i),:,:)=tanm(goodindex(i),:,:)+tbcm(goodindex(i),:,:)
      endwhere
    enddo
   endif
  endif

 if(solarangles(1,1) .ne. 0) then
!! read qality control flags
!! 100000 = final flag
!! 10000  = first guess flag
!! 1000   = departure flag
!! 100    = analysis quality control flag
!! 10     = blacklist flag
!! 1      = analysis TEMP datum event
  tflags=rcpara%miss_val
  read(iunit,iostat=ios) hilf
  if(ios .eq. 0) then
    qc_count=0
    bl_count=0
    do i=1,imax
      mask=.false.
      if(rcpara%qc .eq. 'T') then
        mask=mask .or. (tfgm(goodindex(i),:,:) .ne. rcpara%miss_val .and. abs(tfgm(goodindex(i),:,:)) .gt. rcpara%qcthresh)
      else
        tflags(goodindex(i),:,:)=hilf(i,:,:)
        ifak=100000
        do iflag=1,5
          if(iflag .eq. 5 .and. rcpara%qc .eq. 'B') then
            mask=mask .or. (tfgm(goodindex(i),:,:) .ne. rcpara%miss_val .and. abs(tfgm(goodindex(i),:,:)) .gt. rcpara%qcthresh)
          else
            hilf2=nint(hilf(i,:,:))/ifak
            mask=mask .or. hilf2 .gt. 2
            hilf(i,:,:)=hilf(i,:,:)-hilf2*ifak
            ifak=ifak/10
          endif
        enddo
      endif
      qc_count=qc_count+count(mask)
      where(mask)
       tm(goodindex(i),:,:)=rcpara%miss_val
       tbcm(goodindex(i),:,:)=rcpara%miss_val
       tfgm(goodindex(i),:,:)=rcpara%miss_val
       tanm(goodindex(i),:,:)=rcpara%miss_val
      endwhere
      if(iflag .eq. 5 .and. rcpara%qc .ne. 'T') then
        hilf2=nint(hilf(i,:,:))/ifak
        mask=hilf2 .gt. 2
        bl_count=bl_count+count(mask)
      endif
    enddo
  endif
!!  write(*,'(A17,I6,A17,I6,A17)') 'Operational data: ', qc_count,' rejected by QC; ',bl_count,' blacklisted'


!! read radiosonde type information
  rtype(16072:rcpara%nmax)=nint(rcpara%miss_val)
  rrcor(16072:rcpara%nmax)=nint(rcpara%miss_val)
  read(iunit,iostat=ios) hilf
  if(ios .eq. 0) then
    do i=1,imax
      ttypes(goodindex(i),:,:)=hilf(i,:,:)
      if(any(hilf(i,:,:) .ne. rcpara%miss_val)) then
        val=nint(maxval(hilf(i,:,:)))
        rtype(goodindex(i))=mod(val+1000,1000)
        rrcor(goodindex(i))=val/1000
      endif
    enddo
  endif
  
!!  call rtypeinfo(rtype,rcpara,mr)

  endif

  close(iunit)

  deallocate(hilf,hilf2d,goodindex,overlap)

  return

120 continue
  !!$ call omp_set_lock(omp_lp)
  print*,'could not open file ',filename
  !!$ call omp_unset_lock(omp_lp)
  err=1

  return

end subroutine read_sonde_oper

subroutine rtypeinfo(rtype,rcpara,mr,istat,imax) !

implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer,save :: maxe(5000)=0
integer i,j,ic,max,min,rtype(rcpara%nmax),mr,imax,istat,itypemin(1000),itypemax(1000)

  if(maxe(istat) .eq. 0) then
  mr=1
  imax=1
!!  write(*,*) 'Radiosonde types used:'
itypemin=rcpara%nmax
itypemax=1
      do j=1,rcpara%nmax
        if(rtype(j) .gt. 10 .and. rtype(j) .lt. 90) then
          if(itypemin(rtype(j)) .gt. j) itypemin(rtype(j))=j
          if(itypemax(rtype(j)) .lt. j) itypemax(rtype(j))=j
        endif
      enddo
!!      write(*,*) i,':',ic,' min:',min,' max:',max
      do i=11,89
        if(itypemax(i) .gt. mr .and. itypemax(i)-itypemin(i) .gt. 100) then
          imax=i
          mr=itypemax(i) 
        endif
      enddo

  maxe(istat)=imax

  else
    imax=maxe(istat)
  endif

!!    !$ call omp_set_lock(omp_lp)
!!  write(*,*) 'Latest radiosonde type used:',imax,' at',mr
!!    !$ call omp_unset_lock(omp_lp)
return
end subroutine rtypeinfo

subroutine create_meta(cardsmeta,rtype,iname,rcpara,iter) !

implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer i,j,l,ic,imax(100),imin(100),good(100),rtype(rcpara%nmax),mr,minstat,iname,swap,ilatest
integer cardsmeta(rcpara%nmax,1,rcpara%nmeta),iter

!! initialize place where GTS information is stored.
  cardsmeta(:,1,7)=0

  mr=toindex(19810101,rcpara)
  imax=1
  minstat=rcpara%snht_maxlen/2-rcpara%max_miss
  imin=rcpara%nmax
  l=0
!!  write(*,*) 'Radiosonde types used:'
!! Types used based on WMO BUFR table C-2
  do i=11,89
    ic=count(rtype .eq. i)
    if(ic .gt. minstat) then
      imin(i)=rcpara%nmax
      imax(i)=1
      do j=mr,rcpara%nmax
        if(rtype(j) .eq. i .and. rtype(j) .gt. 10 .and. rtype(j) .lt. 90) then
          if(j .lt. imin(i) .and. j .lt. rcpara%nmax-minstat) then
            if(count(rtype(j:j+minstat) .eq. i) .gt. minstat/2) then
!!               write(*,*) 'imin ', i,j,imin(i),count(rtype(j:j+minstat) .eq. i)
              imin(i)=j
            endif
          endif
          if(imax(i) .lt. j .and. j .gt. minstat) then
            if(count(rtype(j-minstat:j) .eq. i) .gt. minstat/2)then
!!              write(*,*) 'imax ',i,j,imax(i),count(rtype(j-minstat:j) .eq. i)
              imax(i)=j
            endif
          endif
        endif
      enddo
      if(imin(i) .lt. rcpara%nmax) then
    !!$ call omp_set_lock(omp_lp)
        write(*,*) ' GTS: Radiosonde type used at station ',iname,' : ',i,' from ',todate(imin(i),rcpara), ' to ',todate(imax(i),rcpara)
      l=l+1
      good(l)=i
    !!$ call omp_unset_lock(omp_lp) 
      endif
    endif
  enddo

!! changes specific to VIZ sondes. If sonde type is reported, this change is overruled
meta: do i=10000,rcpara%nmax
  if(cardsmeta(i,1,2) .eq. 912) exit meta
enddo meta
ilatest=toindex(19970601,rcpara)
if(i .lt. ilatest) then
   cardsmeta(ilatest,1,1:2)=952 !! change to VIZ B2 assumed
endif

  if(l .lt. 2) then
!!    !$ call omp_set_lock(omp_lp)
!!      write(*,*) ' No GTS documented breaks for station ',iname,' found '
!!    !$ call omp_unset_lock(omp_lp) 
    return
  endif

    
!! bubble sort good to remove first radiosonde type below
  do i=1,l
    do j=i,l
      if(imin(good(j)) .lt. imin(good(i))) then
        swap=good(j)
        good(j)=good(i)
        good(i)=swap
      endif
    enddo
  enddo


!! discard first radiosonde type as metadata since it only indicates first valid type reporting but not necessarily
!! a sonde change
  if(count(cardsmeta(imin(good(1)):rcpara%nmax,1,2) .gt. 0) .gt. 0) then
    !!$ call omp_set_lock(omp_lp)
      write(*,*) ' CARDS sonde change metadata removed for station ',iname,' from ',todate(imin(good(1)),rcpara), ' onwards '
    !!$ call omp_unset_lock(omp_lp) 
    do i=1,l
      cardsmeta(imin(good(i)):imax(good(i)),1,2)=0
    enddo
  endif
  do i=2,l
!! do not assume change if gap in GTS documentation is too large (>2 years)
    if(imax(good(i-1)) .gt. imin(good(i))-rcpara%snht_maxlen/2) then
    !!$ call omp_set_lock(omp_lp)
      if(good(i-1) .eq. 28 .and. good(i) .eq. 11 &
    .or. good(i-1) .eq. 29 .and. good(i) .eq. 12 &
    .or. good(i-1) .eq. 27 .and. good(i) .eq. 45 &
    .or. good(i-1) .eq. 18 .and. (good(i) .eq. 31 .or. good(i) .eq. 37) &
    .or. good(i-1) .eq. 16 .and. good(i) .eq. 37 &
    .or. good(i-1) .eq. 19 .and. good(i) .eq. 32 &
    .or. good(i-1) .eq. 20 .and. good(i) .eq. 28 &
    .or. good(i-1) .eq. 82 .and. (good(i) .ge. 60 .and. good(i) .le. 63)&
    .or. good(i-1) .eq. 81 .and. (good(i) .ge. 60 .and. good(i) .le. 63)&
    ) then
        write(*,*) ' spurious GTS break due to change to FM94 BUFR ',iname,todate(imin(good(i)),rcpara),' from',good(i-1),' to',good(i)
      else  
        if(iter .eq. 3) write(*,'(a,I5,i8,a,i3,a,i3)') ' GTS break added to metadata ',iname,todate(imin(good(i)),rcpara),' from',good(i-1),' to',good(i)
      cardsmeta(imin(good(i)),1,7)=good(i)
      cardsmeta(imin(good(i)),1,2)=good(i)
      endif
    !!$ call omp_unset_lock(omp_lp) 
     endif
  enddo
    !!$ call omp_set_lock(omp_lp)
!!  write(*,*) 'Latest radiosonde type used:',imax,' at',mr
    !!$ call omp_unset_lock(omp_lp)
return
end subroutine create_meta

subroutine read_sonde_daily(iunit,filename,rcpara,err,feld1,feld2,feld3) !

implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer ini,il,ipmax,iparmax,i,iunit,err,imax,ip,ipar
integer goodindex(rcpara%nmax)
character*80 filename
real(kind=JPRM) :: feld1(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM),optional :: feld2(rcpara%nmax,rcpara%pmax,rcpara%parmax),feld3(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM),allocatable :: hilf(:,:,:)
logical l_open

    !!$ call omp_set_lock(omp_lp)
!write(*,*) rcpara%nmax,rcpara%pmax,rcpara%parmax
    !!$ call omp_unset_lock(omp_lp)

  feld1=rcpara%miss_val
  if(present(feld2)) feld2=rcpara%miss_val
  if(present(feld3)) feld3=rcpara%miss_val

  iunit=iunit 
  i=0
!11  inquire(file=filename,opened=l_open)
!  if(l_open) then
!     write(*,*) filename,' already open, retrying'
!     i=i+1
!     if(i .gt. 1000) stop '1000 retries'
!     goto 11
!  endif
  open(iunit,file=filename,form='unformatted',status='old',action='read',err=120)
  err=0
  read(iunit) ini,il,ipmax,iparmax
  read(iunit) goodindex(1:il)

  allocate(hilf(il,ipmax,iparmax))
  hilf=rcpara%miss_val

  imax=0
  do i=1,il
    if(goodindex(i) .le. rcpara%nmax) imax=imax+1
  enddo

  read(iunit) hilf
  do ipar=1,rcpara%parmax
  do ip=1,rcpara%pmax
  do i=1,imax
    feld1(goodindex(i),ip,ipar)=hilf(i,ip,ipar)
  enddo
  enddo
  enddo

  if(present(feld2)) then 
    read(iunit) hilf
   do ipar=1,rcpara%parmax
  do ip=1,rcpara%pmax
  do i=1,imax
      feld2(goodindex(i),ip,ipar)=hilf(i,ip,ipar)
    enddo
   enddo
  enddo
 endif

  if(present(feld3)) then 
    read(iunit) hilf
   do ipar=1,rcpara%parmax
  do ip=1,rcpara%pmax
    do i=1,imax
      feld3(goodindex(i),ip,ipar)=hilf(i,ip,ipar)
    enddo
   enddo
  enddo
 endif

  deallocate(hilf)

  close(iunit)

  return

120 print*,'could not open file ',filename
  err=1

  return

end subroutine read_sonde_daily

subroutine read_sonde_new(iunit,filename,ni,pmax,parmax,tm,tanm,tfgm,tbcm,solarangles,miss_val,err)

implicit none

integer ni,pmax,parmax,ini,il,ipmax,iparmax,i,iunit,err,imax,ipar
integer goodindex(ni),bi(parmax)
character*(*) filename
real(kind=JPRM) :: tm(ni,pmax,parmax),tanm(ni,pmax,parmax),tfgm(ni,pmax,parmax),tbcm(ni,pmax,parmax),solarangles(ni,parmax)
real(kind=JPRM) :: miss_val
real(kind=JPRM),allocatable :: hilf(:,:)


  tm=miss_val
  tanm=miss_val
  tfgm=miss_val
  tbcm=miss_val
  solarangles=miss_val

  iunit=iunit
  open(iunit,file=filename,form='unformatted',status='old',action='read',err=120)
  err=0
  read(iunit) ini,ipmax,iparmax
  read(iunit) bi(1:iparmax)

  do ipar=1,iparmax

  allocate(hilf(bi(ipar),ipmax))

  read(iunit) goodindex(1:bi(ipar))

  imax=0
  do i=1,bi(ipar)
    if(goodindex(i) .le. ni) imax=imax+1
  enddo

  read(iunit) hilf
  do i=1,imax
    tm(goodindex(i),:,ipar)=hilf(i,:)
  enddo

  read(iunit) hilf
  do i=1,imax
    tanm(goodindex(i),:,ipar)=hilf(i,:)
  enddo

  read(iunit) hilf
  do i=1,imax
    tfgm(goodindex(i),:,ipar)=hilf(i,:)
  enddo

  read(iunit) hilf(:,1)
  do i=1,imax
    solarangles(goodindex(i),ipar)=hilf(i,1)
  enddo

  read(iunit) 
  do i=1,imax
    tbcm(goodindex(i),:,ipar)=hilf(i,:)
  enddo

  deallocate(hilf)

  enddo

!! reverse effect of bias correction
  where(tfgm .ne. -999.) tfgm=tfgm+tbcm

  close(iunit)

  return

120 print*,'could not open file ',filename
  err=1

  return

end subroutine read_sonde_new

subroutine write_sonde_corr_daily(iunit,filename,rcpara,err,rasocorr,eracorr,reserve) !

implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer ini,ipmax,iparmax,i,iunit,err,imax,l,bi,ip,ipar
integer index(rcpara%mmax)
character*(*) filename
real(kind=JPRM) :: rasocorr(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM),optional :: eracorr(rcpara%nmax,rcpara%pmax,rcpara%parmax),reserve(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: hilfcorr(rcpara%mmax,rcpara%pmax,rcpara%parmax)
logical jump

  iunit=iunit

!!  print*,filename
  open(iunit,file=filename,form='unformatted',err=120)
  err=0

  bi=1
  hilfcorr(1,:,:)=rasocorr(1,:,:)
  index(1)=1
  do i=2,rcpara%nmax-1
    jump=.false.
    do ipar=1,rcpara%parmax
      ip=1
      do while(ip .le. rcpara%pmax .and. .not. jump)
        jump=abs(rasocorr(i,ip,ipar)-rasocorr(i+1,ip,ipar)).gt. 0.001
        ip=ip+1
      enddo
    enddo
    if(jump .and. bi .lt. rcpara%mmax) then
      bi=bi+1
      hilfcorr(bi,:,:)=rasocorr(i+1,:,:)
      index(bi)=i+1
    endif
  enddo

  if(bi .eq. rcpara%mmax) write(*,*) 'Warning: Maximum allowed break number ',rcpara%mmax,' reached'

  write(iunit) bi,rcpara%pmax,rcpara%parmax
  write(iunit) index(1:bi)
  write(iunit) hilfcorr(1:bi,:,:)

  if(present(eracorr)) then
    bi=1
    hilfcorr(1,:,:)=eracorr(1,:,:)
    index(1)=1
    do i=2,rcpara%nmax-1
      if(any(eracorr(i,:,:) .ne. eracorr(i+1,:,:)) .and. bi .lt. rcpara%mmax) then
        bi=bi+1
        hilfcorr(bi,:,:)=eracorr(i+1,:,:)
        index(bi)=i+1
      endif
    enddo

    write(iunit) bi,rcpara%pmax,rcpara%parmax
    write(iunit) index(1:bi)
    write(iunit) hilfcorr(1:bi,:,:)
  endif

  if(present(reserve)) then
    bi=1
    hilfcorr(1,:,:)=eracorr(1,:,:)
    index(1)=1
    do i=2,rcpara%nmax-1
      if(any(reserve(i,:,:) .ne. reserve(i+1,:,:)) .and. bi .lt. rcpara%mmax) then
        bi=bi+1
        hilfcorr(bi,:,:)=reserve(i+1,:,:)
        index(bi)=i+1
      endif
    enddo

    write(iunit) bi,rcpara%pmax,rcpara%parmax
    write(iunit) index(1:bi)
    write(iunit) hilfcorr(1:bi,:,:)
  endif

  close(iunit)

  return

120 continue
    !!$ call omp_set_lock(omp_lp)
print*,'could not open file ',filename
    !!$ call omp_unset_lock(omp_lp)
  err=1

  return

end subroutine write_sonde_corr_daily

subroutine write_sonde_corr_daily8(iunit,filename,rcpara,err,rasocorr) !

implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer ini,ipmax,iparmax,i,iunit,err,imax,l,bi,ip,ipar
integer index(rcpara%mmax)
character*(*) filename
real(kind=JPRM) :: rasocorr(rcpara%nmax,rcpara%pmax,rcpara%parmax,8)
real(kind=JPRM) :: hilfcorr(rcpara%mmax,rcpara%pmax,rcpara%parmax,8)
logical jump

  iunit=iunit

!!  print*,filename
  open(iunit,file=filename,form='unformatted',err=120)
  err=0

  bi=1
  hilfcorr(1,:,:,:)=rasocorr(1,:,:,:)
  index(1)=1
  do i=2,rcpara%nmax-1
    jump=.false.
    ipar=1
    do while(ipar .le. rcpara%parmax .and. .not. jump)
      l=1
      do while(l .le. 8 .and. .not. jump)
        ip=1
        do while(ip .le. rcpara%pmax .and. .not. jump)
          jump=abs(rasocorr(i,ip,ipar,l)-rasocorr(i+1,ip,ipar,l)).gt. 0.001
          ip=ip+1
        enddo
        l=l+1
      enddo
      ipar=ipar+1
    enddo
    if(jump .and. bi .lt. rcpara%mmax) then
      bi=bi+1
      hilfcorr(bi,:,:,:)=rasocorr(i+1,:,:,:)
      index(bi)=i+1
    endif
  enddo

  if(bi .eq. rcpara%mmax) write(*,*) 'Warning: Maximum allowed break number ',rcpara%mmax,' reached'

  write(iunit) bi,rcpara%pmax,rcpara%parmax
  write(iunit) index(1:bi)
  write(iunit) hilfcorr(1:bi,:,:,:)


  close(iunit)

  return

120 continue
    !!!$ call omp_set_lock(omp_lp)
print*,'could not open file ',filename
    !!!$ call omp_unset_lock(omp_lp)
  err=1

  return

end subroutine write_sonde_corr_daily8

subroutine read_sonde_corr_daily_IFS(iunit,filename,rcpara,estatindex,err,ifs_rasocorrs,ifs_index,ifs_eracorrs,reserve)

implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer ini,ipmax,iparmax,i,iunit,err,imax,l,bi,estatindex
integer index(rcpara%mmax)
character*(*) filename
real(kind=JPRM) :: hilfcorr(rcpara%mmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: ifs_rasocorrs(rcpara%mmax,rcpara%pmax,rcpara%parmax,rcpara%statmax)
real(kind=JPRM),optional :: ifs_eracorrs(rcpara%mmax,rcpara%pmax,rcpara%parmax,rcpara%statmax),reserve(rcpara%mmax,rcpara%pmax,rcpara%parmax,rcpara%statmax)
integer(kind=JPRM) :: ifs_index(rcpara%mmax,rcpara%statmax)

  open(iunit,file=filename,form='unformatted',status='old',action='read',err=120)
  err=0

  read(iunit) bi,ipmax,iparmax
  if(bi .gt. rcpara%brmax .or. ipmax .ne. rcpara%pmax .or. iparmax .ne. rcpara%parmax) print*,'error: ',filename,bi,ipmax,iparmax
  read(iunit) index(1:bi)
  read(iunit) hilfcorr(1:bi,:,:)

  do i=1,bi
    ifs_rasocorrs(i,:,:,estatindex)=hilfcorr(i,:,:)
    ifs_index(i,estatindex)=index(i)
  enddo

  write(*,*),'readIFS:',filename,bi,estatindex,ifs_index(1:bi,estatindex)
  if(present(ifs_eracorrs)) then
    read(iunit) bi,ipmax,iparmax
  if(bi .gt. rcpara%brmax .or. ipmax .ne. rcpara%pmax .or. iparmax .ne. rcpara%parmax) print*,'error: ',filename,bi,ipmax,iparmax
    read(iunit) index(1:bi)
    read(iunit) hilfcorr(1:bi,:,:)

    do i=1,bi
      ifs_eracorrs(i,:,:,estatindex)=hilfcorr(i,:,:)
!!    ifs_index(i,estatindex)=index(i)
    enddo
  endif

  if(present(reserve)) then
    read(iunit) bi,ipmax,iparmax
  if(bi .gt. rcpara%brmax .or. ipmax .ne. rcpara%pmax .or. iparmax .ne. rcpara%parmax) print*,'error: ',filename,bi,ipmax,iparmax
    read(iunit) index(1:bi)
    read(iunit) hilfcorr(1:bi,:,:)

    do i=1,bi
      reserve(i,:,:,estatindex)=hilfcorr(i,:,:)
!!    ifs_index(i,estatindex)=index(i)
    enddo
  endif

return
120 if(err .ne. 1) print*,'could not open file ',filename
  err=1
return
end subroutine read_sonde_corr_daily_IFS

subroutine write_leobiascor_IFS_table(iunit,filename,rcpara,estatmax,ifs_rasocorrs,ifs_index,ewmonrs,ewmolats,ewmolons,err) !

implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer ini,ipmax,iparmax,i,iunit,err,imax,l,ib,estatmax,istat,idum
integer index(50)
character*80 filename
character :: ccorr*1,cpmax*2
real(kind=JPRM) :: ifs_rasocorrs(rcpara%mmax,rcpara%pmax,rcpara%parmax,estatmax)
integer(kind=JPRM) :: ifs_index(rcpara%mmax,estatmax)
integer :: ewmonrs(estatmax)
real(kind=JPRM) :: ewmolats(estatmax),ewmolons(estatmax)

!!  filename=trim(rcpara%prefix)//'leobiascor.t'
  open(iunit,file=filename,form='formatted',err=120)
  err=0

    !!$ call omp_set_lock(omp_lp)
  write(iunit,'(A5,I5,I3,16F6.0)') rcpara%version,estatmax,rcpara%pmax,rcpara%plevs
    !!$ call omp_unset_lock(omp_lp)
  do istat=1,estatmax
    ccorr='N'
    ib=1
    if(ifs_index(ib,istat) .gt. 0 .and. .not. any(ifs_rasocorrs(ib,:,:,istat) .eq. rcpara%miss_val)) then
      ccorr='Y'
      do while(ifs_index(ib,istat) .gt. 0 .and. .not. any(ifs_rasocorrs(ib,:,:,istat) .eq. rcpara%miss_val)) 
        ib=ib+1
      enddo
      ib=ib-1
    write(iunit,'(I5,I6,2F9.2,1X,A1,I3)') istat,ewmonrs(istat),ewmolats(istat),ewmolons(istat),ccorr,ib
    do i=1,ib
      write(cpmax,'(I2.2)') rcpara%pmax
      write(iunit,'(I4,3I2.2,'//cpmax//'F6.2)') rcpara%year(ifs_index(i,istat)),rcpara%month(ifs_index(i,istat)),rcpara%day(ifs_index(i,istat)),0,ifs_rasocorrs(i,:,1,istat)
      write(iunit,'(I4,3I2.2,'//cpmax//'F6.2)') rcpara%year(ifs_index(i,istat)),rcpara%month(ifs_index(i,istat)),rcpara%day(ifs_index(i,istat)),12,ifs_rasocorrs(i,:,2,istat)

    enddo

    else
      write(iunit,'(I5,I6,2F9.2,1X,A1,I3)') istat,ewmonrs(istat),ewmolats(istat),ewmolons(istat),ccorr,0
    endif
  enddo

  close(iunit)
return
120 print*,'could not open file ',filename
  err=1
return
end subroutine write_leobiascor_IFS_table

!! read_leobiascor_IFS_table reads table leobiascor.t
!! and yields estimated radiosonde temperature biases for
!! a given radiosonde and a given date
!!
!! input:
!! iunit    = filenr
!! filename = path to leo's bias correction table. Table is read only  once  
!!            (at first call of the subroutine)
!! wmonr    = WMO station identifier
!! year,month,day,time specify time of radiosonde launch
!! plevs [hPa]    = pressure levels. Routine assumes 16 levels
!! pmax [hPa]    =  array size of plev
!!
!! output:
!! lcorr    = correction flag; .t. if station wmonr has been corrected
!! bias     = contains temperature biases at pressure levels given in plevs
!!            contains zeros of err .gt. 0
!! err      = 0 no error
!!          = 1 could not read file filename
!!          = 2 station wmonr not found in table. Only those wmonrs as defined
!!              in table stgroups.t are defined.
!!          = 3 at least one pressure level requested is not available
!!              a zero bias is returned at these levels
!!          = 4 invalid time, format=YYYYMMDDHH, only 00,12 is allowed
!!
!! Leo Haimberger, 21 September 2004
!!
subroutine read_leobiascor_IFS_table(iunit,filename,wmonr,year,month,day,time,plevs,pmax,lcorr,bias,version,err,startdate)

implicit none

save

integer ni,pmax,ini,ipmax,iparmax,i,iunit,err,imax,l,istat,istatmax,wmonr
integer, parameter :: brmax=10
integer, parameter :: parmax=2
integer, parameter :: estatmax=3070
integer index(50)
integer  iyear,imonth,itime,iday,datum,idum
integer,optional :: startdate

character*(*) filename
character :: ccorr*1,cpmax*2,version*5
real(kind=JPRM) :: ifs_rasocorrs(brmax,pmax,parmax,estatmax),iplevs(pmax),plevs(pmax),bias(pmax)
integer(kind=JPRM) :: ifs_index(brmax,estatmax),ib(estatmax)
real(kind=JPRM) ::  miss_val
integer year,month,day,time
integer :: ewmonrs(estatmax)
real(kind=JPRM) :: ewmolats(estatmax),ewmolons(estatmax)
logical lcorr,lread


  if(.not. lread) then
  open(iunit,file=filename,form='formatted',status='old',action='read',err=120)
  err=0

  read(iunit,'(A5,I5,I3,16F6.0)') version,istatmax,ipmax,iplevs
  write(cpmax,'(I2.2)') ipmax

  do istat=1,estatmax
    read(iunit,'(I5,I6,2F9.2,1X,A1,I3)') idum,ewmonrs(istat),ewmolats(istat),ewmolons(istat),ccorr,ib(istat)
!!    write(*,*) idum,ewmonrs(istat),ewmolats(istat),ewmolons(istat),ccorr,ib(istat)
    if(ccorr .eq. "Y") then
      do i=1,ib(istat)
        read(iunit,'(I4,3I2.2,'//cpmax//'F6.2)') iyear,imonth,iday,itime,ifs_rasocorrs(i,:,1,istat)
        ifs_index(i,istat)=iyear*1000000+imonth*10000+iday*100
        read(iunit,'(I4,3I2.2,'//cpmax//'F6.2)') iyear,imonth,iday,itime,ifs_rasocorrs(i,:,2,istat)
      enddo
    endif
  enddo

  close(iunit)
  endif

  bias=0.
  lcorr=.false.

  if(pmax .ne. 16 .or. any(plevs .ne. iplevs)) then
     err=3
     return
  endif

  datum=year*1000000+month*10000+day*100
  if(present(startdate)) then 
  if(datum .lt. startdate .or. datum .gt. 2002083112 .or. time .ne. 0 .and. time .ne. 12) then
    err=4
    return
  endif
  else
  if(datum .lt. 19570101 .or. datum .gt. 2002083112 .or. time .ne. 0 .and. time .ne. 12) then
    err=4
    return
  endif
  endif

  do istat=1,estatmax
    if(wmonr .eq. ewmonrs(istat)) then
      if(ib(istat) .eq. 0) then
        lcorr=.false.
      else
       do i=1,ib(istat)
        if(ifs_index(i,istat) .lt. datum) then
          if(time .eq. 0) then
            bias=ifs_rasocorrs(i,:,1,istat)
          endif
          if(time .eq. 12) then
            bias=ifs_rasocorrs(i,:,2,istat)
          endif
          lcorr=.true.
        endif
       enddo
       exit
      endif
    endif
  enddo

  if(istat .eq. estatmax+1) then
    err=2
    return
  endif

return
120 print*,'could not open file ',filename
  err=1
return
end subroutine read_leobiascor_IFS_table

!! read_leobiascor_IFS_table reads table leobiascor.t
!! and yields estimated radiosonde temperature biases for
!! a given radiosonde and a given date
!!
!! input:
!! iunit    = filenr
!! filename = path to leo's bias correction table. Table is read only  once  
!!            (at first call of the subroutine)
!! wmonr    = WMO station identifier
!! year,month,day,time specify time of radiosonde launch
!! plevs [hPa]    = pressure levels. Routine assumes 16 levels
!! pmax [hPa]    =  array size of plev
!!
!! output:
!! lcorr    = correction flag; .t. if station wmonr has been corrected
!! bias     = contains temperature biases at pressure levels given in plevs
!!            contains zeros of err .gt. 0
!! err      = 0 no error
!!          = 1 could not read file filename
!!          = 2 station wmonr not found in table. Only those wmonrs as defined
!!              in table stgroups.t are defined.
!!          = 3 at least one pressure level requested is not available
!!              a zero bias is returned at these levels
!!          = 4 invalid time, format=YYYYMMDDHH, only 00,12 is allowed
!!
!! Leo Haimberger, 21 September 2004
!!
subroutine read_biascor(iunit,filename,corr_profiles,corrtimes,ewmonrs,plevs,pmax,parmax,brmax,estatmax,version,err)

implicit none

integer ni,pmax,ini,ipmax,iparmax,i,iunit,err,imax,l,istat,istatmax,brmax,estatmax,parmax
integer index(brmax)
integer  iyear,imonth,itime,iday,datum,idum

character*(*) filename
character :: ccorr*1,cpmax*2,version*5
real(kind=JPRM) :: corr_profiles(pmax,parmax,brmax,estatmax),iplevs(pmax),plevs(pmax),bias(pmax)
integer(kind=JPRM) :: corrtimes(brmax,estatmax),ib(estatmax)
real(kind=JPRM) ::  miss_val
integer year,month,day,time
integer :: ewmonrs(estatmax)
real(kind=JPRM) :: ewmolats(estatmax),ewmolons(estatmax)
logical lcorr,lread

  if(.not. lread) then
  open(iunit,file=filename,form='formatted',status='old',action='read',err=120)
  err=0

  read(iunit,'(A5,I5,I3,16F6.0)') version,istatmax,ipmax,iplevs
  write(cpmax,'(I2.2)') ipmax

  do istat=1,estatmax
    read(iunit,'(I5,I6,2F9.2,1X,A1,I3)') idum,ewmonrs(istat),ewmolats(istat),ewmolons(istat),ccorr,ib(istat)
    if(ccorr .eq. "Y") then
      do i=1,ib(istat)
        read(iunit,'(I4,3I2.2,'//cpmax//'F6.2)') iyear,imonth,iday,itime,corr_profiles(:,1,i,istat)
        corrtimes(i,istat)=iyear*1000000+imonth*10000+iday*100
        read(iunit,'(I4,3I2.2,'//cpmax//'F6.2)') iyear,imonth,iday,itime,corr_profiles(:,2,i,istat)
      enddo
    endif
  enddo

  close(iunit)
  endif

return
120 print*,'could not open file ',filename
  err=1
return
end subroutine read_biascor

subroutine expand_biascor(corr_profiles,corrtimes,rasocorr,statnr,wmonrs,plevs,nmax,pmax,parmax,brmax,estatmax)

implicit none

integer nmax,pmax,parmax,brmax,estatmax,i,istat,istatmax,wmonrs,statnr
integer index

real(kind=JPRM) :: corr_profiles(pmax,parmax,brmax,estatmax),iplevs(pmax),plevs(pmax),bias(pmax)
real(kind=JPRM) :: rasocorr(nmax,pmax,parmax)
integer(kind=JPRM) :: corrtimes(brmax,estatmax),ib
real(kind=JPRM) ::  miss_val
integer :: ewmonrs(estatmax)
integer(kind=4) :: year(20000),month(20000),day(20000),time(20000)
integer ::         year1,month1,day1,time1,tdiff
real(kind=JPRM) :: ewmolats(estatmax),ewmolons(estatmax)

do istat=1,estatmax
  if(statnr .eq. wmonrs(istat)) exit
enddo
if(statnr .ne. wmonrs(istat)) write(*,*) statnr,' not found in station table'

year(1)=1957
month(1)=1
day(1)=1
time(1)=00
tdiff=24
year1=1957
index=1
do while(index .lt. nmax)
  CALL DATES(Year(index),Month(index),Day(index),Time(index),Tdiff, &
             Year1,Month1,Day1,Time1)
  index=index+1
  year(index)=year1
  month(index)=month1
  day(index)=day1
  time(index)=time1
!!  write(*,*) year1,month1,day1,time1
enddo

ib=1
index=1
do while(corrtimes(ib+1,istat) .gt. 0 .and. ib .lt. brmax)
  do while(year(index)*1000000+month(index)*10000+day(index)*100 .lt. corrtimes(ib+1,istat))
    rasocorr(index,:,:)=corr_profiles(:,:,ib,istat)
    index=index+1
  enddo
  ib=ib+1
enddo

return

end subroutine expand_biascor

subroutine write_igrasonde_daily(iunit,filename,ni,pmax,parmax,tm,tanm,tfgm,tfg12m,miss_val,err)

implicit none

integer ni,pmax,parmax,ini,ipmax,iparmax,i,iunit,err,imax,l,bi
integer goodindex(ni)
character*(*) filename
real(kind=JPRM) :: tm(ni,pmax,parmax),tfgm(ni,pmax,parmax),tanm(ni,pmax,parmax),tfg12m(ni,pmax,parmax)
real(kind=JPRM) ::  miss_val

  iunit=iunit

  open(iunit,file=filename,form='unformatted',err=120)
  err=0

  bi=0
  do i=1,ni
    if(any(tm(i,:,:) .ne. miss_val)) then
      bi=bi+1
      goodindex(bi)=i
      tm(bi,:,:)=tm(i,:,:)
      tanm(bi,:,:)=tanm(i,:,:)
      tfgm(bi,:,:)=tfgm(i,:,:)
      tfg12m(bi,:,:)=tfg12m(i,:,:)
    endif
  enddo

  write(iunit) ni,bi,pmax,parmax
  write(iunit) goodindex(1:bi)
  write(iunit) tm(1:bi,:,:)
  write(iunit) tanm(1:bi,:,:)
  write(iunit) tfgm(1:bi,:,:)
  write(iunit) tfg12m(1:bi,:,:)

  close(iunit)

  return

120 print*,'could not open file ',filename
  err=1

  return

end subroutine write_igrasonde_daily

subroutine write_igrasonde_daily_new(iunit,filename,ni,pmax,parmax,tm,tanm,tfgm,tfg12m,miss_val,err) !

implicit none

integer ni,pmax,parmax,ini,ipmax,iparmax,i,iunit,err,imax,l,ipar,bi(parmax)
integer goodindex(ni,parmax),ihilf(ni)
character*(*) filename
real(kind=JPRM) :: tm(ni,pmax,parmax)
real(kind=JPRM) :: tfgm(ni,pmax,parmax),tanm(ni,pmax,parmax),tfg12m(ni,pmax,parmax)
real(kind=JPRM) ::  miss_val,rhilf(ni)

  iunit=iunit

  open(iunit,file=filename,form='unformatted',err=120)
  err=0

  do ipar=1,parmax
  bi(ipar)=0
  do i=1,ni
    if(any(tm(i,:,ipar) .ne. miss_val)) then
      bi(ipar)=bi(ipar)+1
      goodindex(bi(ipar),ipar)=i
      tm(bi(ipar),:,ipar)=tm(i,:,ipar)
      tanm(bi(ipar),:,ipar)=tanm(i,:,ipar)
      tfgm(bi(ipar),:,ipar)=tfgm(i,:,ipar)
      tfg12m(bi(ipar),:,ipar)=tfg12m(i,:,ipar)
    endif
  enddo
  enddo

  write(iunit) ni,pmax,parmax
  write(iunit) bi
!!  write(*,*) filename,bi
  do ipar=1,parmax
    write(iunit) goodindex(1:bi(ipar),ipar) 
    write(iunit) tm(1:bi(ipar),:,ipar)
    write(iunit) tanm(1:bi(ipar),:,ipar)
    write(iunit) tfgm(1:bi(ipar),:,ipar)
    write(iunit) tfg12m(1:bi(ipar),:,ipar)
  enddo


  close(iunit)

  return

120 print*,'could not open file ',filename
  err=1

  return

end subroutine write_igrasonde_daily_new

subroutine makemonthly(rcpara,tm,tmmon,thresh) !

implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer i,iunit,err,l,iy,index,imon,thresh,ipar,ip,imod,tgm

real(kind=JPRM),intent(in) :: tm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: tmmon(rcpara%mmax,rcpara%pmax,rcpara%parmax)

  do ipar=1,rcpara%parmax
   do ip=1,rcpara%pmax
    index=1
    do imon=1,rcpara%mmax
      iy=rcpara%startdate/10000+(imon-1)/12
      imod=mod(imon-1,12)+1
      tmmon(imon,ip,ipar)=0.
      tgm=0
      do while(rcpara%month(index) .eq. imod .and. index .lt. rcpara%nmax)
        if(tm(index,ip,ipar) .ne. rcpara%miss_val) then
          tmmon(imon,ip,ipar)=tmmon(imon,ip,ipar)+tm(index,ip,ipar)
          tgm=tgm+1
        endif
        index=index+1
      enddo
        
      if(tgm .ge. thresh) then
         tmmon(imon,ip,ipar)=tmmon(imon,ip,ipar)/tgm
      else
         tmmon(imon,ip,ipar)=rcpara%miss_val
      endif   
    enddo
   enddo
  enddo

  return

end subroutine makemonthly

subroutine write_sonde_monthly(iunit,filename,rcpara,tm,err,rasocorr,eracorr) !

implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer i,iunit,err,l,iy,index,imon
character*(*) filename
real(kind=JPRM),intent(in) :: tm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM),intent(in),optional ::rasocorr(rcpara%nmax,rcpara%pmax,rcpara%parmax),eracorr(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: itm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: tmmon(rcpara%mmax,rcpara%pmax,rcpara%parmax),rasocorrmon(rcpara%mmax,rcpara%pmax,rcpara%parmax),eracorrmon(rcpara%mmax,rcpara%pmax,rcpara%parmax)
integer         :: goodmon(rcpara%mmax,rcpara%pmax,rcpara%parmax),imod,tgm,rgm,egm,ip,ipar
logical         :: lrc,ler

  iunit=iunit

  open(iunit,file=filename,form='unformatted',err=120)
  err=0

  lrc=present(rasocorr)
  ler=present(eracorr)
  do ipar=1,rcpara%parmax
   do ip=1,rcpara%pmax
    index=1
    do imon=1,rcpara%mmax
      iy=rcpara%startdate+(imon-1)/12
      imod=mod(imon-1,12)+1
      tmmon(imon,ip,ipar)=0.
      rasocorrmon(imon,ip,ipar)=0.
      eracorrmon(imon,ip,ipar)=0.
      tgm=0
      rgm=0
      egm=0
      do while(rcpara%month(index) .eq. imod .and. index .lt. rcpara%nmax)
        if(tm(index,ip,ipar) .ne. rcpara%miss_val) then
          tmmon(imon,ip,ipar)=tmmon(imon,ip,ipar)+tm(index,ip,ipar)
          tgm=tgm+1
        endif
        if(lrc) then
        if(rasocorr(index,ip,ipar) .ne. rcpara%miss_val) then
          rasocorrmon(imon,ip,ipar)=rasocorrmon(imon,ip,ipar)+rasocorr(index,ip,ipar)
          rgm=rgm+1
        endif
        endif
        if(ler) then 
        if(eracorr(index,ip,ipar) .ne. rcpara%miss_val) then
          eracorrmon(imon,ip,ipar)=eracorrmon(imon,ip,ipar)+eracorr(index,ip,ipar)
          egm=egm+1
        endif
        endif
        index=index+1
!!        write(*,*) index,imon,iy, mod(imon-1,12)+1,sum(goodmon(imon,ip,ipar))
      enddo
        
      goodmon(imon,ip,ipar)=tgm
      if(tgm .gt. 0) then
         tmmon(imon,ip,ipar)=tmmon(imon,ip,ipar)/tgm
      else
         tmmon(imon,ip,ipar)=rcpara%miss_val
      endif   
      if(rgm .gt. 0) then
        rasocorrmon(imon,ip,ipar)=rasocorrmon(imon,ip,ipar)/rgm
      else
        rasocorrmon(imon,ip,ipar)=rcpara%miss_val
      endif
      if(egm .gt. 0) then
        eracorrmon(imon,ip,ipar)=eracorrmon(imon,ip,ipar)/egm
      else
        eracorrmon(imon,ip,ipar)=rcpara%miss_val
      endif   
    enddo
   enddo
  enddo


    write(iunit) rcpara%mmax,rcpara%pmax,rcpara%parmax
    write(iunit) tmmon
    if(lrc) write(iunit) rasocorrmon
    if(ler)write(iunit) eracorrmon
!!    write(*,*) trim(filename)//':',count(tm .ne. rcpara%miss_val),' good values'
!!    write(*,*) trim(filename)//':',sum(goodmon),' good values'
    write(iunit) goodmon
    close(iunit)

  return

120 print*,'could not open file ',filename
  err=1

return
end subroutine write_sonde_monthly

!! read_sonde_monthly for reading monthly mean  RAOBCORE-raw input 
!! and adjustments
!! first month is 195701
!! mmax = maximum number of months (588)
!! pmax = maximum number of pressure levels. Pressure levels are:
!! 10.,20.,30.,50.,70.,100.,150.,200.,250.,300.,400.,500.,700.,850.,925.,1000.
!! parmax is maximum number of ascents per day (2 for 00/12GMT)
!! real*4 arrays are:
!! tmmon(mmax,pmax,parmax) = adjusted RAOBCORE temperatures from merged dataset
!! rasocorrmon(mmax,pmax,parmax) = adjustments applied. Subtracting rasocorrmon from tmmon yields unadjusted temperatures from merged RAOBCORE input dataset
!! eracorrmon(mmax,pmax,parmax) = adjustments applied to ERA-40 bg before using it as a reference. 
!! missing value = -999.
!! RAOBCORE data are big_endian, like the IBM machines at ECMWF. 
!! Therefore the routine has to be compiled with -byteswapio (pgf90) or
!! -convert big_endian (ifort)
subroutine read_sonde_monthly(iunit,filename,mmax,pmax,parmax,tmmon,rasocorrmon,eracorrmon,goodmon,err)

implicit none

integer iunit,err,mmax,pmax,parmax,lmmax,lpmax,lparmax
character*(*) filename
real(kind=4) :: tmmon(mmax,pmax,parmax),rasocorrmon(mmax,pmax,parmax),eracorrmon(mmax,pmax,parmax)
integer         :: goodmon(mmax,pmax,parmax)

  iunit=iunit

  open(iunit,file=filename,form='unformatted',err=120)
  err=0


  read(iunit) lmmax,lpmax,lparmax
  if(lmmax .ne. mmax .or. lpmax .ne. pmax .or. lparmax .ne. parmax) then
    write(*,*) 'dimensions of file        ' ,lmmax,lpmax,lparmax
    write(*,*) 'do not agree with expected' ,mmax,pmax,parmax
  endif
  read(iunit) tmmon
  read(iunit) rasocorrmon
  read(iunit) eracorrmon
  read(iunit) goodmon
  close(iunit)

  return

120 print*,'could not open file ',filename
  err=1

return
end subroutine read_sonde_monthly


subroutine write_sonde_daily(iunit,filename,rcpara,err,feld1,feld2,feld3) !

implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer ini,ipmax,iparmax,i,iunit,err,imax,l,il,ios
integer goodindex(rcpara%nmax)
character*(*) filename
real(kind=JPRM) :: feld1(rcpara%nmax,rcpara%pmax,rcpara%parmax),hilf(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM),optional :: feld2(rcpara%nmax,rcpara%pmax,rcpara%parmax),feld3(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM),allocatable :: hilf2(:,:,:)
logical l_open

i=0
!11  inquire(file=filename,opened=l_open)
!  if(l_open) then
!    i=i+1
!    write(*,*) filename, ' already open, retrying'
!    if(i .gt. 1000) stop 'giving up after 1000 retries'
!    goto 11
!  endif

  open(iunit,file=filename,form='unformatted',action='write',iostat=ios,err=120)
  err=0

  il=1
  do i=1,rcpara%nmax
    if(any(feld1(i,:,:) .gt. -999.)) then
      goodindex(il)=i
      hilf(il,:,:)=feld1(i,:,:)
      if(present(feld2)) then 
        feld2(il,:,:)=feld2(i,:,:)
      endif
      if(present(feld3)) then 
        feld3(il,:,:)=feld3(i,:,:)
      endif
      il=il+1
    endif
  enddo
  il=il-1
  write(iunit) rcpara%nmax,il,rcpara%pmax,rcpara%parmax
  write(iunit) goodindex(1:il)
  allocate(hilf2(il,rcpara%pmax,rcpara%parmax))
  hilf2=hilf(1:il,:,:)
  write(iunit) hilf2
  if(present(feld2)) then
    hilf2=feld2(1:il,:,:)
    write(iunit) hilf2
  endif
  if(present(feld3)) write(iunit) feld3(1:il,:,:)

!! write(iunit) newtbcs(1:indexmax,:,1,:),newradbcs(1:indexmax,:,1,:)
  deallocate(hilf2)

  close(iunit)
return

120 print*,'could not open file ',filename , err, ios
  err=1
return
end subroutine write_sonde_daily

subroutine read_supersonde(iunit,filename,ni,pmax,parmax,tm,tanm,tfgm,tbcm,tnumbers,miss_val,err)

implicit none

integer ni,pmax,parmax,iunit,err
character*(*) filename
real(kind=JPRM) :: tm(ni,pmax,parmax),tanm(ni,pmax,parmax),tfgm(ni,pmax,parmax),tbcm(ni,pmax,parmax)
integer(kind=JPRM) :: tnumbers(ni,pmax,parmax)
real(kind=JPRM) :: miss_val
  iunit=iunit

  open(iunit,file=filename,form='unformatted',status='old',action='read',err=120)
  err=0
  read(iunit) tm,tanm,tfgm,tbcm,tnumbers
  close(iunit)

  return

120 print*,'could not open file ',filename
  err=1

  return

end subroutine read_supersonde

subroutine read_sonde_cards(iunit,filename,ni,pmax,parmax,tm,tanm,tfgm,tfg12,err)

implicit none

integer ni,pmax,parmax,iunit,err,ini,ipmax,iparmax
character*(*) filename
real(kind=JPRM) :: tm(ni,pmax,parmax),tanm(ni,pmax,parmax),tfgm(ni,pmax,parmax),tfg12(ni,pmax,parmax)
integer(kind=JPRM) :: tnumbers(ni,pmax,parmax)

  iunit=iunit

  open(iunit,file=filename,form='unformatted',status='old',action='read',err=120)
  err=0
  read(iunit) ini,ipmax,iparmax
  read(iunit) tm,tanm,tfgm,tfg12
  close(iunit)

  return

120 print*,'could not open file ',filename
  err=1

  return

end subroutine read_sonde_cards


subroutine read_cards_meta(iunit,rcpara,wmostats,statnr,wmonrs,cardsmeta,index,hilf,il,err) !

implicit none

type(rasocor_namelist),intent(in) :: rcpara
integer iunit,err,imeta,wmostats,i,istat,it,statnr,statindex,ios,iter,il,istatsave,imetasave,upperbound
character*80 filename
integer(kind=JPRM) :: cardsmeta(23010,1,rcpara%nmeta),wmonrs(wmostats)
integer(kind=JPRM) :: index(3000000),hilf(3000000),ihilf
integer,save :: ini,iwmostats
integer,allocatable,save :: stats(:,:),metas(:)
  cardsmeta=0

  statindex=0
  do i=1,wmostats
    if(statnr .eq. wmonrs(i)) statindex=i
  enddo
  
  if(statindex .eq. 0) then 
    print*,'no metadata for ',statnr
    return
  else
!!    print*,'reading metadata for station',statnr,' index ',statindex
  endif

  if(iunit .ne. 0) then 
    if(wmostats .eq. 963) then
      open(iunit,file='./cardsmeta.0963',form='unformatted',status='old',action='read',err=120)
    else
      open(iunit,file='./cardsmeta',form='unformatted',status='old',action='read',err=120)
    endif
  err=0
    read(iunit,err=121,iostat=ios) ini,iwmostats,imeta,il
    read(iunit,err=121,iostat=ios) index(1:il)
    read(iunit,err=121,iostat=ios) hilf(1:il)
    close(iunit)

    allocate(metas(9),stats(wmostats+1,9))
    imetasave=0
    istatsave=0
    metas=0
    stats=0
    do i=1,il
      imeta=index(i)/ini/(wmostats-0)+1
      if(imeta .gt. imetasave) then
        imetasave=imeta
        metas(imeta)=i
      endif
      istat=(index(i)-(imeta-1)*ini*(wmostats-0))/ini+1
      if(istat .ne. istatsave) then
        istatsave=istat
        stats(istat,imeta)=i
      endif
    enddo
  endif

  do imeta=1,rcpara%nmeta
    ihilf=statindex+1
    do while (ihilf .lt. wmostats .and. stats(ihilf,imeta) .eq. 0)
      ihilf=ihilf+1
    enddo    
    if(ihilf .lt. wmostats) then
      upperbound=stats(ihilf,imeta)-1
    else
      upperbound=il
    endif
!!    write(*,*) imeta,statindex,stats(statindex,imeta),upperbound
    if(stats(statindex,imeta) .gt. 0) then 
      do i=stats(statindex,imeta),upperbound
        istat=(index(i)-(imeta-1)*ini*(wmostats-0))/ini+1
        it=index(i)-(imeta-1)*ini*(wmostats-0)-ini*(istat-1)+1
        if(istat .ne. statindex) then
!!     !$ call omp_set_lock(omp_lp)
!!         write(*,*) 'wrong station',istat,statindex,i,index(i)
!!     !$ call omp_unset_lock(omp_lp)
        else
          cardsmeta(it,1,imeta)=hilf(i)
        endif
      enddo
    endif
  enddo

  return

120 continue
    !!$ call omp_set_lock(omp_lp)
print*,'could not open file ' !!,filename
    !!$ call omp_unset_lock(omp_lp)
  err=1

  return
121 continue
    !!$ call omp_set_lock(omp_lp)
print*,'could not read file ' !!,filename
    !!$ call omp_unset_lock(omp_lp)
  err=1

  return

end subroutine read_cards_meta

subroutine read_cards_meta_schroeder(iunit,rcpara,statnr,wmonrs,wmolons,wmolats,wmostats,meta,err) !

use rfmod
use sort
implicit none


type(rasocor_namelist),intent(in) :: rcpara
type(metadata) :: meta
integer iunit,err,wmostats,i,istat,it,statnr,statindex,ios,statid,statidsave,bsearch,idatesave,dindex,l,lmax,mstat,merged
character*80 filename,zeile
integer(kind=JPRM) :: wmonrs(:),id,zindex(8),idate
real :: lwmonrs(wmostats),lat,lon
integer :: lindex(wmostats),llwmonrs(wmostats)
character*3 :: rscode
character*2 chmod
character*9 key
 real(kind=JPRM) :: wmolats(:),wmolons(:)
integer :: tsamaxloc(rcpara%probmax),fgbreaks(3),index
real(kind=JPRM) :: tsamax(rcpara%probmax)
logical :: trusted(5000),ldebug

  ldebug=.false.
  lwmonrs=wmonrs(1:wmostats)
  call qsort(lwmonrs,lindex)
  llwmonrs=floor(lwmonrs)
  statidsave=0
  idatesave=0
  meta.rscodes='   '
  meta.rsicodes=0
  merged=0

  meta.trusted=.false.
  if(iunit .ne. 0) then 
    filename='./vapor.instruments.1'
    open(iunit,file=trim(filename),status='old',action='read',err=120)
  err=0
    lmax=0
    do while (.not. eof(iunit))
       lmax=lmax+1
       read(iunit,'(I5,13X,A3)',err=120,iostat=ios) meta.rsicodes(lmax),meta.rscodes(lmax)
       if((index(meta.rscodes(lmax),'V') .eq. 1 .or. index(meta.rscodes(lmax),'J').eq. 1) .and. index('12',meta.rscodes(lmax)(2:2))<1 .and. meta.rscodes(lmax).ne.'V8A') then
         meta.trusted(lmax)=.true.
       endif
    enddo

    filename='./vapor.obsolete.1'
    open(iunit,file=trim(filename),status='old',action='read',err=120)
    do while (.not. eof(iunit))
       lmax=lmax+1
       read(iunit,'(I5,7X,A2)',err=120,iostat=ios) meta.rsicodes(lmax),meta.rscodes(lmax)
       if(index(meta.rscodes(lmax),'V') .eq. 1) meta.trusted(lmax)=.true.
    enddo
    
    id=0
    do i=1,lmax
      if (meta.trusted(i)) id=id+1
    enddo
    write(*,*) 'Trusted stations:',id
    filename='./vapor.library.1'
    open(iunit,file=trim(filename),status='old',action='read',err=120)
  err=0
    zeile=''
    do while (.not. eof(iunit))

       read(iunit,'(A80)',err=120,iostat=ios) zeile
       id=1
       zindex(1)=1
       do i=1,80
          if(zeile(i:i) .eq. ':') then
             id=id+1
             zindex(id)=i
             if(id .eq. 8) exit
          endif
       enddo
       if(index(zeile,'91165') .ne. 0) then
         print*,zeile
       endif
       read(zeile(zindex(6)+1:zindex(7)-5),'(I8)') idate
       if(idate .eq. 88888888 .or. idate .eq. 19000000) cycle
       if(zindex(8)-zindex(7) .ne. 3) cycle
       read(zeile(zindex(7)+1:zindex(8)-1),'(A2)') chmod
       if(chmod(1:1) .eq. '#' .or. chmod(1:1) .eq. '/') cycle
       if(chmod(2:2) .eq. '-' .or. chmod(2:2) .eq. '/') cycle
       if(index(zeile(zindex(2)+1:zindex(4)-1),'?') .eq. 0) then
       read(zeile(zindex(2)+1:zindex(3)-1),'(F8.2)'),lat
       read(zeile(zindex(3)+1:zindex(4)-1),'(F8.2)'),lon
       else
         lat=rcpara%miss_val
         lon=rcpara%miss_val
       endif
       read(zeile(zindex(5)+1:zindex(6)-1),'(A3)'),rscode
       read(zeile(1:5),'(I5)'), statid
       if(statid .ne. statidsave) then
          istat=bsearchi(statid,llwmonrs)
          if(istat .eq. -1) then
! station id not found, but may be merged with other station ID
! check if nearby station exists
! slow - linear search
            do mstat=1,rcpara%statmax
              if(abs(wmolats(mstat)-lat) .lt. 0.3 .and. abs(wmolons(mstat)-lon) .lt. 0.3) then
                write(*,*) 'metadata of ',statid,'merged into', wmonrs(mstat)
                exit
              endif 
            enddo  
            if (mstat .lt. rcpara%statmax) then
               istat=mstat    
               merged=merged+1      
            endif
          else
            istat=lindex(istat)
          endif
          statidsave=statid
       endif
       if(idate/10000 .ge. rcpara%startdate/10000 .and.  idate/10000 .lt. 2100 .and. istat .ne. -1) then
          if(mod(idate,10000) .eq. 0) idate=idate+101
          if(mod(idate,10000) .gt. 1231) idate=(idate/10000)*10000+101
          if(mod(idate,100) .eq. 0) idate=idate+1
          if(mod(idate,100) .gt. 31) idate=idate-mod(idate,100)+1

          dindex=toindex(idate,rcpara,.false.)
          if(dindex .ne. -1) then
            l=1
            do while(l .le. lmax)
             if(rscode .eq. meta.rscodes(l)) exit
             l=l+1
            enddo
            if(l .gt. lmax .or. l .eq. 1) then
              write(*,*) 'rs type ',rscode,' not found'
            else
              meta.cardsmeta_s(dindex,istat)=l
            endif
          endif
          if(ldebug) print*,idate,wmonrs(istat),dindex
          if(idate .eq. idatesave) then 
             if(ldebug)  print*,'date twice'
          else
            idatesave=idate
          endif
          if(dindex .eq. -1) then
             print*,'wrong date'
          endif
       endif
   
    enddo
    print *,count(meta.cardsmeta_s .ne. 0),' metadata events read'
    print *,merged,' stations merged'
    close(iunit)

    filename='./jumps'
    open(iunit,file=trim(filename),status='old',action='read',err=120)
    err=0
    zeile=''
    fgbreaks=(/toindex(19741231,rcpara),toindex(19760901,rcpara),toindex(19790101,rcpara)/)
    do while (.not. eof(iunit))
      read(iunit,'(I5,A9,12F9.1)') statid,key,tsamax
      read(iunit,'(I5,A9,12I9)') statid,key,tsamaxloc
          istat=bsearchi(statid,llwmonrs)
          do i=1,12
            if(tsamaxloc(i) .gt. 0) tsamaxloc(i)=toindex(tsamaxloc(i),rcpara)
          enddo
          if(any(tsamax(7:12) .ne. rcpara%miss_val)) then
          do i=1,4
            if(tsamax(i) .gt. 300 .and. .not. any(abs(tsamaxloc(i)-fgbreaks).lt. 730)) then
               meta.cardsmeta_s(tsamaxloc(i),lindex(istat))=1000
write(*,*) 'jumpmeta',statid,todate(tsamaxloc(i),rcpara)
            endif
            if(tsamax(i) .gt. 400 ) then
               meta.cardsmeta_s(tsamaxloc(i),lindex(istat))=1000
write(*,*) 'jumpmeta',statid,todate(tsamaxloc(i),rcpara)
            endif
          enddo
            if(tsamax(5) .gt. 100 .and. .not. any(abs(tsamaxloc(5)-fgbreaks).lt. 730)) then
               meta.cardsmeta_s(tsamaxloc(5),lindex(istat))=1000
write(*,*) 'jumpmeta5',statid,todate(tsamaxloc(5),rcpara)
            endif
          endif
    enddo
  endif
  
!!$  open(iunit,file='schroed.dump',form='unformatted')
!!$  write(iunit) rcpara%nmax,wmostats
!!$  write(iunit) lwmonrs
!!$  write(iunit) meta.cardsmeta
!!$  close(iunit)
 
  return

120 continue
    !!$ call omp_set_lock(omp_lp)
print*,'read_cards_meta_schroeder: could not open file ' ,filename
    !!$ call omp_unset_lock(omp_lp)
return
end subroutine read_cards_meta_schroeder


subroutine read_corr(iunit,filename,tm,tfgm,cstart,cstop,plevs,year,month,day,ni,pmax,parmax,miss_val)

integer ni,pmax,parmax,istart,istop,i,ip,ipar,iunit,ips,ios,cstart,cstop
integer year(20000),month(20000),day(20000)
character*(*) filename

real(kind=JPRM) :: tfgm(ni,pmax,parmax),tm(ni,pmax,parmax),corr(parmax),plevs(pmax),miss_val

1 open(iunit,file=filename,form='formatted',status='old',action='read',iostat=ios,err=120)

istart=1
do i=1,ni
  if(year(i)*10000+month(i)*100+day(i) .eq. cstart) istart=i
  if(year(i)*10000+month(i)*100+day(i) .eq. cstop) istop=i
enddo

l=1
do while(.true.)
  read(iunit,'(I6,2F6.2)',end=100) ips,corr
  do while(plevs(l) .ne. ips)
    l=l+1
  enddo
  
  where(tfgm(istart:istop,l,1) .ne. miss_val) tfgm(istart:istop,l,1)=tfgm(istart:istop,l,1)+corr(1)  
  where(tfgm(istart:istop,l,2) .ne. miss_val) tfgm(istart:istop,l,2)=tfgm(istart:istop,l,2)+corr(2)  
  
enddo

100 continue
!!print*,l
close(iunit)

return

120 print*,'could not read ',filename, ' trying again'
if(ios .eq. 23) goto 1

return
end subroutine read_corr

subroutine read_initial_correct(iunit,filename,rcpara,wmostats,wmonrs,ini_correct,err)

implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer wmostats,wmonrs(wmostats),err,iunit,ios,ianf,iend
logical ini_correct(rcpara%parmax:wmostats)
character*(*) filename
character*100 zeile

open(iunit,file=filename,form='formatted',status='old',action='read',iostat=ios,err=120)

ini_correct=.false.
do while(.true.)
1  read(iunit,'(A100)',end=100) zeile
  if(zeile(1:1) .eq. '#') goto 1
!!$ call omp_set_lock(omp_lp)
  read(zeile,'(I5,1X,I5)') ianf,iend
!!$ call omp_unset_lock(omp_lp)
  where(wmonrs .ge. ianf .and. wmonrs .le. iend)
    ini_correct=.true.
  endwhere
enddo

100 continue
!!print*,l
close(iunit)

return

120 print*,'read_initial_correct: could not read ',filename

return
end subroutine read_initial_correct

subroutine read_corr_bin(iunit,filename,brcorrs,brcorrdates,brmax,cstart,cstop,rcpara,statmax)

type(rasocor_namelist),intent(in):: rcpara

integer istart,istop,i,ip,ipar,iunit,ips,ios,cstart,cstop,statmax,istatmax,ipmax,iparmax,brmax
character*(*) filename

real(kind=JPRM) :: brcorrs(statmax,rcpara%pmax,rcpara%parmax,brmax),corr(statmax,rcpara%pmax,rcpara%parmax)
integer(kind=JPRM) :: brcorrdates(2,brmax)

1 open(iunit,file=filename,form='unformatted',status='old',action='read',iostat=ios,err=120)

read(iunit) istatmax,ipmax,iparmax
read(iunit) corr

i=1
do while(brcorrs(1,1,1,i) .ne. rcpara%miss_val)
  i=i+1
enddo
brcorrs(:,:,:,i)=corr
brcorrdates(1,i)=cstart
brcorrdates(2,i)=cstop

close(iunit)

return

120 print*,'could not read '//filename
call abort

return
end subroutine read_corr_bin

end module rfcorio
