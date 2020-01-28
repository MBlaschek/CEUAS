subroutine raso_correct_igra_ei(rcpara,wmonrs,wmolats,wmolons,wmonames,wmostats,iwmonrs,iwmolats,iwmolons,iwmonames,iwmostats,istatmax,ominuse40,used,needed,statnr,istat)!

use rfmod
use rfcor
use txtnc

implicit none

type(rasocor_namelist),intent(in) :: rcpara
 
integer,intent(in) :: istat

integer         :: maxlen,max_miss,i,select_mode,iunit,ilmeta,istatmax,eerr,eerr2,ierr,ierr2,ip,ipar,l,m,k
real(kind=JPRM) :: locsig,break_fak

real(kind=JPRM) :: diff(rcpara%nmax),plus(rcpara%nmax),minus(rcpara%nmax),prms(rcpara%nmax),mrms(rcpara%nmax),tsa(rcpara%nmax)
integer         :: pcount(rcpara%nmax),mcount(rcpara%nmax)


real(kind=JPRM) :: tfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tanm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tbcm(rcpara%nmax,rcpara%pmax,rcpara%parmax),stfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: stype(rcpara%nmax),stype_tmp(rcpara%nmax)
real(kind=JPRM) :: tfgmei(rcpara%nmax,rcpara%pmax,rcpara%parmax),tanmei(rcpara%nmax,rcpara%pmax,rcpara%parmax),tmei(rcpara%nmax,rcpara%pmax,rcpara%parmax),tbcmei(rcpara%nmax,rcpara%pmax,rcpara%parmax),fflags(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM),allocatable,dimension(:,:,:,:) :: arr
real(kind=JPRM),allocatable,dimension(:) :: sstype
INTEGER,allocatable,dimension(:,:) :: hdatum,hhours
real(kind=JPRM) :: itfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax),itanm(rcpara%nmax,rcpara%pmax,rcpara%parmax),itm(rcpara%nmax,rcpara%pmax,rcpara%parmax),itfg12m(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: mtfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax),mtanm(rcpara%nmax,rcpara%pmax,rcpara%parmax),mtm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM)  :: solarangles(rcpara%nmax,rcpara%parmax)
integer(kind=JPRM)  :: tnumbers(rcpara%nmax,rcpara%pmax,rcpara%parmax)

real(kind=JPRM),intent(in) :: ominuse40(rcpara%ni,rcpara%nj,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: ominuse40_an(rcpara%ni,rcpara%nj,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: adjust(rcpara%pmax,rcpara%parmax),adjust_an(rcpara%pmax,rcpara%parmax),eps
real :: dists(rcpara%statmax),idists(istatmax)
integer rtype(rcpara%nmax)

integer,intent(in) :: wmonrs(rcpara%statmax),iwmonrs(istatmax)
real(kind=JPRM),intent(in) :: wmolats(rcpara%statmax),wmolons(rcpara%statmax),iwmolats(istatmax),iwmolons(istatmax)
character*20,intent(in) :: wmonames(rcpara%statmax),iwmonames(istatmax)
!!character*64 :: stgroupnames(:)
character prefix*20,binprefix*20,startdate*8,cdatum*8,cstatnr*6,igcstatnr*6,filename*80,cgstat*5
integer tdiff,eramax
integer,intent(in) :: statnr,iwmostats,wmostats
integer nmeta,nmem,ngroup,err,snht_increment,mmax,first,last,snht_maxlen
integer smooth_method,prob_method,protunit,j,iswap
logical ex
integer index(50)
integer needed(rcpara%statmax),is,mini,emini

logical*1 used(rcpara%nmax,rcpara%pmax,rcpara%parmax,2*rcpara%statmax/3),merged(rcpara%nmax,rcpara%pmax,rcpara%parmax)

real(kind=JPRM) :: plevs(rcpara%pmax),mindists,emindists
INTEGER :: datum(rcpara%nmax,1), hours(rcpara%nmax,rcpara%parmax)


iunit=statnr+100000 !! this should avoid clashes in parallel runs

write(igcstatnr,'(I6.6)') statnr

stype=rcpara%miss_val
call sphdist(iwmolats(istat),iwmolons(istat),wmolats,wmolons,dists,rcpara%statmax) !in file rfmod.f90 line 2257
call sphdist(iwmolats(istat),iwmolons(istat),iwmolats,iwmolons,idists,istatmax) !in file rfmod.f90 line 2257
eps=0.01
mindists=180.
mini=2000
do i=1,istatmax
  if(idists(i) .gt. eps .and. idists(i) .lt. mindists) then
     mindists=idists(i)
     mini=i
  endif
enddo
emindists=180.
emini=2000
do i=1,rcpara%statmax
  if(dists(i) .gt. eps .and. dists(i) .lt. emindists) then
     emindists=dists(i)
     emini=i
  endif
enddo
solarangles(1,1)=0.
write(*,*) igcstatnr//': Nearest IGRA neighbour',mindists,' away'
!!call read_igrasonde_daily_new(iunit,'igrafgbin'//igcstatnr,rcpara%nmax,rcpara%pmax,rcpara%parmax,itm,itanm,itfgm,itfg12m,rcpara%miss_val,ierr)
call read_sonde(iunit,trim(rcpara%prefix)//'/igrafgbinera'//igcstatnr,rcpara,itm,itanm,itfgm,itfg12m,solarangles,rtype,ierr2,.true.) !in file rfcorio.f90 line 373

l=toindex(rcpara%switchdate,rcpara)
!l=toindex(19890100,rcpara)
!itm(l:rcpara%nmax,:,:)=rcpara%miss_val
!itfgm(l:rcpara%nmax,:,:)=rcpara%miss_val

i=1
do while(rcpara%year(i) .lt. 2013) 
  i=i+1
enddo

write(*,*) count(itm(366:i-20,:,:) .ne. rcpara%miss_val),count(itfgm(366:i-20,:,:) .ne. rcpara%miss_val),count(itanm(366:i-20,:,:) .ne. rcpara%miss_val),count(itfg12m(366:i-20,:,:) .ne. rcpara%miss_val)
!! calculate analysis adjustment oper/era40
!!call read_sonde_oper(iunit,trim(rcpara%prefix)//'/igrax'//igcstatnr,rcpara,iwmolats(istat),iwmolons(istat),itm,itanm,itfgm,itfg12m,solarangles,ominuse40_an,adjust_an,rtype,ierr2,.true.) !in file rfcorio.f90 line 611

call read_sonde(iunit,trim(rcpara%prefix)//'../igraei/igrafgbinoper'//igcstatnr,rcpara,mtm,mtanm,mtfgm,itfg12m,solarangles,rtype,ierr2,.true.) !in file rfcorio.f90 line 373
mtm(1:l-1,:,:)=itm(1:l-1,:,:)
mtfgm(1:l-1,:,:)=itfgm(1:l-1,:,:)
mtanm(1:l-1,:,:)=itanm(1:l-1,:,:)

k=toindex(19881231,rcpara)
call read_sonde(iunit,trim(rcpara%prefix)//'../igraei79/igrafgbinoper'//igcstatnr,rcpara,itm,itanm,itfgm,itfg12m,solarangles,rtype,ierr2,.true.) !in file rfcorio.f90 line 373
mtm(l:k,:,:)=itm(l:k,:,:)
mtfgm(l:k,:,:)=itfgm(l:k,:,:)
mtanm(l:k,:,:)=itanm(l:k,:,:)

mtm(20455+365:20455+366,:,:)=rcpara%miss_val ! 31.12.2013 fehlt in ERA-Interim
mtfgm(20455+365:20455+366,:,:)=rcpara%miss_val ! 31.12.2013 fehlt in ERA-Interim
mtanm(20455+365:20455+366,:,:)=rcpara%miss_val ! 31.12.2013 fehlt in ERA-Interim

l=toindex(rcpara%switchdate,rcpara)
!call read_sonde_oper(iunit,trim(rcpara%prefix)//'/igrafgbinoper2'//igcstatnr,rcpara,iwmolats(istat),iwmolons(istat),mtm,mtanm,mtfgm,itfg12m,solarangles,ominuse40,adjust,rtype,ierr2,.true.) !in file rfcorio.f90 line 611

write(*,*) count(itm(366:i-20,:,:) .ne. rcpara%miss_val),count(itfgm(366:i-20,:,:) .ne. rcpara%miss_val),count(itanm(366:i-20,:,:) .ne. rcpara%miss_val),count(itfg12m(366:i-20,:,:) .ne. rcpara%miss_val)

eramax=0
do i=1,rcpara%statmax
  if(dists(i) .lt. mindists .and. dists(i) .lt. 0.5) then
!! this is to avoid merging with station with a stationID that is included
!! in IGRA but has not yet been merged. Example:
!! Funchal has ID 60018 since 2002 and had 60020 (and a slightly different location) before that date
!! Station should be merged with 60020 and not with 60018
    if(iwmonrs(istat) .eq. wmonrs(i) .or. .not. any(iwmonrs(istat+1:istatmax) .eq. wmonrs(i))) then
      eramax=eramax+1
      index(eramax)=i
      write(*,'(2I6,5F8.2)') iwmonrs(istat),wmonrs(i),iwmolats(istat),iwmolons(istat),dists(i),wmolats(i),wmolons(i)
    else
      write(*,*) wmonrs(i),' will be merged later'
    endif
  endif
enddo

!! do a bubblesort, is efficient enough in this case
if(eramax .gt. 1) then
  do j=1,eramax-1
   do i=1,eramax-1
    if(dists(index(i)) .gt. dists(index(i+1))) then
      iswap=index(i) 
      index(i)=index(i+1)
      index(i+1)=iswap
    endif
   enddo
  enddo
endif

if(eramax .eq. 0) then
  write(*,*) igcstatnr, ' NO MATCHES: ',mindists,iwmonrs(mini),emindists,wmonrs(emini)
  if(emindists .lt. 0.5) then
    write(*,*) igcstatnr, 'NO MERGE since probably included in other station ID'
    return
  else
    write(*,*) igcstatnr, 'IGRA data treated as merged'
  endif 
else
  write(*,*) igcstatnr, ' MATCHES: ',eramax,'INDEX:',index(1:eramax)
endif

merged=.false. !!  true if ERA-40 value from nearer station has already been substituted

!mtm=itm
!mtfgm=itfgm

solarangles=rcpara%miss_val
do i=1,eramax
  write(cstatnr,'(I6.6)') wmonrs(index(i))

  filename=trim(rcpara%prefix)//'/feedbackglobbinsaveera'//cstatnr
  write(*,*) 'ERA: ',filename
  call read_sonde(iunit,filename,rcpara,tm,tanm,tfgm,tbcm,solarangles,rtype,eerr) !in file rfcorio.f90 line 373
  write(*,*) igcstatnr,'ERA-40',count(tm .ne. rcpara%miss_val),count(tfgm .ne. rcpara%miss_val) 
tm(l:rcpara%nmax,:,:)=rcpara%miss_val
tfgm(l:rcpara%nmax,:,:)=rcpara%miss_val
tanm(l:rcpara%nmax,:,:)=rcpara%miss_val

! ACHTUNG 1 statt l
hours(1:rcpara%nmax,:)=rcpara%miss_val

!tbcm(l:rcpara%nmax,:,:)=rcpara%miss_val

!! calculate analysis adjustment oper/era40
!  filename=trim(rcpara%prefix)//cstatnr//'/feedbackx'//cstatnr
!  call read_sonde_oper(iunit,filename,rcpara,wmolats(index(i)),wmolons(index(i)),tm,tanm,tfgm,tbcm,solarangles,ominuse40_an,adjust_an,rtype,eerr2) !in file rfcorio.f90 line 611
!  filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbinsaveoper'//cstatnr
!  call read_sonde_oper(iunit,filename,rcpara,wmolats(index(i)),wmolons(index(i)),tm,tanm,tfgm,tbcm,solarangles,ominuse40,adjust,rtype,eerr2) !in file rfcorio.f90 line 611

  filename=trim(rcpara%prefix)//'../ei/'//cstatnr//'/'//cstatnr//'_t.nc'
  CALL read_odb_nc(filename,rcpara,0,eerr2,tmei,tfgmei,tbcm=tbcmei,fflags=fflags,tanm=tanmei,stype=stype_tmp,hours=hours) !subroutine in read_txt_write_nc.f90
  where(stype .eq. rcpara%miss_val)
    stype=stype_tmp
  endwhere
!  print*,stype
   do ipar=1,rcpara%parmax
     do ip=1,rcpara%pmax
       do j=l,rcpara%nmax
          if(tfgmei(j,ip,ipar) .ne. rcpara%miss_val) then
            tfgm(j,ip,ipar)=tfgmei(j,ip,ipar)
            tanm(j,ip,ipar)=tanmei(j,ip,ipar)
            tm(j,ip,ipar)=tmei(j,ip,ipar)
!            tbcm(j,ip,ipar)=tbcmei(j,ip,ipar)
          endif
!          tfgm(i,ip,ipar)=flags(i,ip,ipar)
       enddo
     enddo
   enddo   



tbcm=0.

  write(*,*) igcstatnr,'ERA-40',count(tm .ne. rcpara%miss_val),count(tfgm .ne. rcpara%miss_val)

  eerr=eerr*eerr2

  if(eerr .eq. 0) then

!! merge IGRA and ERA-40
!! rule: ERA-40 has precedence over IGRA since it has more stringent QC
!!       If 2 or more ERA-40 data available use nearest one.

  if(needed(index(i)) .eq. 0) then
    needed(index(i))=maxval(needed)+1
    write(*,*) 'needed increased to',needed(index(i))
    where(tm .ne. rcpara%miss_val .and. .not. merged)
      mtm=tm
      mtfgm=tfgm
      mtanm=tanm
!!      used(:,:,:,needed(index(i)))=.true.
      merged=.true. 
    endwhere
  else
    write(*,*) 'ERA-40 station',wmonrs(index(i)),' matches at least twice - used only once'
 
  endif

  endif
  
  write(*,*) igcstatnr, count(tm .ne. rcpara%miss_val),' values in ',cstatnr

enddo

  write(*,*) igcstatnr, count(merged),' values merged in from ERA data'
  write(*,*) igcstatnr, ' original:',count(itm .ne. rcpara%miss_val),' merged:',count(mtm .ne. rcpara%miss_val)


!! do some final quality control
write(*,*) igcstatnr,' obs-bg values above 20:',count(abs(mtfgm) .gt. 20. .and. mtfgm .ne. rcpara%miss_val)
write(*,*) igcstatnr,' obs-bg values above 30:',count(abs(mtfgm) .gt. 30. .and. mtfgm .ne. rcpara%miss_val)

where(abs(mtfgm) .gt. 20. .and. mtfgm .ne. rcpara%miss_val)
  mtm=rcpara%miss_val
  tbcm=rcpara%miss_val
  mtfgm=rcpara%miss_val
  mtanm=rcpara%miss_val
endwhere

!call write_sonde_monthly(iunit,filename,rcpara,mtm,err,tanm,tbcm) !in file rfcorio.f90 line 1823

tm=mtm !! mtm destroyed in call to write_sonde_daily but needed in write_sonde_monthly
tfgm=mtfgm
tanm=mtanm
filename=trim(rcpara%prefix)//'../ei/'//igcstatnr//'/feedbackmerged'//igcstatnr//'.nc'
!call write_igrasonde_daily_new(iunit,filename,rcpara%nmax,rcpara%pmax,rcpara%parmax,mtm,mtanm,mtfgm,itfg12m,rcpara%miss_val,err) !in file rfcorio.f90 line 1730

l=0
do i=1,rcpara%nmax
  if(any(tm(i,:,:) .ne. rcpara%miss_val)) then
    l=l+1
    datum(l,1)=i
  endif
enddo

ALLOCATE(arr(l,rcpara%pmax,rcpara%parmax,5),hdatum(l,1), hhours(l,rcpara%parmax),sstype(l))
! es wird obs-bg abgespeichert, nicht bg-obs!
where(tfgm .ne. -999.) tfgm=-tfgm
where(tanm .ne. -999.) tanm=-tanm
do i=1,l
  arr(i,:,:,2)=tfgm(datum(i,1),:,:)
  arr(i,:,:,1)=tm(datum(i,1),:,:)
  arr(i,:,:,3)=tbcm(datum(i,1),:,:)
  arr(i,:,:,4)=fflags(datum(i,1),:,:)
  arr(i,:,:,5)=tanm(datum(i,1),:,:)
  hhours(i,:)=hours(datum(i,1),:)
  sstype(i)=stype(datum(i,1))
enddo
hdatum=datum(1:l,:)

!print*,stype(datum(1:l,1))
call txttonc(rcpara,0,filename,arr,(/'temperatures','fg_dep','bias','flags','an_dep','s_type'/),l,5, hdatum, hhours, iwmolats(istat),iwmolons(istat), 0,&
(/'merged temperature','merged fg-departures','bias correction','flags-final,fg,depar, blacklist, andepar2big, fg2big','merged an-departures','sonde types'/),sstype)

  CALL read_odb_nc(filename,rcpara,0,eerr2,tmei,tfgmei,tbcmei,fflags,tanmei,stype_tmp,hours=hours) !subroutine in read_txt_write_nc.f90

filename=trim(rcpara%prefix)//'../ei/'//igcstatnr//'/feedbackbc'//igcstatnr//'.nc'
call write_sonde_monthly_nc(filename,rcpara,tbcmei,istat,err, iwmolats(istat),iwmolons(istat))

filename=trim(rcpara%prefix)//'../ei/'//igcstatnr//'/feedbackan'//igcstatnr//'.nc'
call write_sonde_monthly_nc(filename,rcpara,tanm,istat,err, iwmolats(istat),iwmolons(istat))

deallocate(arr,hdatum,hhours)

return
end


