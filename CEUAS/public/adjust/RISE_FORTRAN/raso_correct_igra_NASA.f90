subroutine raso_correct_igra_NASA(rcpara,wmonrs,wmolats,wmolons,wmonames,wmostats,iwmonrs,iwmolats,iwmolons,iwmonames,iwmostats,istatmax,ominuse40,used,needed,statnr,istat,priority) !

use rfmod
use rfcor

implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer,intent(in) :: istat
character*(*),intent(in)  :: priority

integer         :: maxlen,max_miss,i,select_mode,iunit,ilmeta,istatmax,eerr,eerr2,ierr,ierr2,ip,ipar
real(kind=JPRM) :: locsig,break_fak

real(kind=JPRM) :: diff(rcpara%nmax),plus(rcpara%nmax),minus(rcpara%nmax),prms(rcpara%nmax),mrms(rcpara%nmax),tsa(rcpara%nmax)
integer         :: pcount(rcpara%nmax),mcount(rcpara%nmax)


real(kind=JPRM) :: tfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tanm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tm(rcpara%nmax,rcpara%pmax,rcpara%parmax),tbcm(rcpara%nmax,rcpara%pmax,rcpara%parmax),stfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: stanm(rcpara%nmax,rcpara%pmax,rcpara%parmax),stm(rcpara%nmax,rcpara%pmax,rcpara%parmax),stbcm(rcpara%nmax,rcpara%pmax,rcpara%parmax),eracorrs(rcpara%nmax,rcpara%pmax,rcpara%parmax),rasocorrs(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: itfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax),itanm(rcpara%nmax,rcpara%pmax,rcpara%parmax),itm(rcpara%nmax,rcpara%pmax,rcpara%parmax),itfg12m(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: itfgmNASA(rcpara%nmax,rcpara%pmax,rcpara%parmax),itanmNASA(rcpara%nmax,rcpara%pmax,rcpara%parmax),itmNASA(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: mtfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax),mtanm(rcpara%nmax,rcpara%pmax,rcpara%parmax),mtm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM)  :: solarangles(rcpara%nmax,rcpara%parmax)
integer(kind=JPRM)  :: tnumbers(rcpara%nmax,rcpara%pmax,rcpara%parmax)

real(kind=JPRM),intent(in) :: ominuse40(rcpara%ni,rcpara%nj,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) ::ominuse40_an(rcpara%ni,rcpara%nj,rcpara%pmax,rcpara%parmax)
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


iunit=statnr+100000 !! this should avoid clashes in parallel runs

write(igcstatnr,'(I6.6)') statnr

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
call read_sonde(iunit,trim(rcpara%prefix)//igcstatnr//'/igrafgbinera'//igcstatnr,rcpara,itm,itanm,itfgm,itfg12m,solarangles,rtype,ierr2,.true.) !in file rfcorio.f90 line 373

i=1
do while(rcpara%year(i) .lt. 2007) 
  i=i+1
enddo

write(*,*) count(itm(366:i-20,:,:) .ne. rcpara%miss_val),count(itfgm(366:i-20,:,:) .ne. rcpara%miss_val),count(itanm(366:i-20,:,:) .ne. rcpara%miss_val),count(itfg12m(366:i-20,:,:) .ne. rcpara%miss_val)
!! calculate analysis adjustment oper/era40
call read_sonde_oper(iunit,trim(rcpara%prefix)//igcstatnr//'/igrax'//igcstatnr,rcpara,iwmolats(istat),iwmolons(istat),itm,itanm,itfgm,itfg12m,solarangles,ominuse40_an,adjust_an,rtype,ierr2,.true.) !in file rfcorio.f90 line 611

call read_sonde_oper(iunit,trim(rcpara%prefix)//igcstatnr//'/igrafgbinoper'//igcstatnr,rcpara,iwmolats(istat),iwmolons(istat),itm,itanm,itfgm,itfg12m,solarangles,ominuse40,adjust,rtype,ierr2,.true.) !in file rfcorio.f90 line 611

write(*,*) count(itm(366:i-20,:,:) .ne. rcpara%miss_val),count(itfgm(366:i-20,:,:) .ne. rcpara%miss_val),count(itanm(366:i-20,:,:) .ne. rcpara%miss_val),count(itfg12m(366:i-20,:,:) .ne. rcpara%miss_val)

call read_igra_NASA(iunit,igcstatnr,rcpara,itmNASA) !in this file line 299

write(*,*) igcstatnr,': Non-NASA values',count(itm .ne. rcpara%miss_val .and. itmNASA .eq. rcpara%miss_val)
write(*,*) igcstatnr,': New-NASA values',count(itm .ne. rcpara%miss_val .and. itmNASA .ne. rcpara%miss_val)
write(*,*) igcstatnr,': Orphaned-NASA values',count(itm .eq. rcpara%miss_val .and. itmNASA .ne. rcpara%miss_val)

do ipar=1,rcpara%parmax
  do ip=1,rcpara%pmax
    do i=1,rcpara%nmax
      if(itm(i,ip,ipar) .ne. rcpara%miss_val .and. itmNASA(i,ip,ipar) .ne. rcpara%miss_val) then
        itfg12m(i,ip,ipar)=itmNASA(i,ip,ipar)-itm(i,ip,ipar)
        itm(i,ip,ipar) = itmNASA(i,ip,ipar)
        itfgm(i,ip,ipar)=itfgm(i,ip,ipar)-itfg12m(i,ip,ipar)
        itanm(i,ip,ipar)=itanm(i,ip,ipar)-itfg12m(i,ip,ipar)
      endif
    enddo
  enddo
enddo

!!where(itm .ne. rcpara%miss_val .and. itmNASA .ne. rcpara%miss_val) 
!!  itfg12m=itmNASA-itm
!!  itfgmNASA=itfgm+itfg12m
!!  itanmNASA=itanm+itfg12m
!!elsewhere
!!  itfg12m=-5.
!!  itfgmNASA=rcpara%miss_val
!!  itanmNASA=rcpara%miss_val
!!  itfgm=-5.
!!  itanm=-5.
!!endwhere

!!itfgm=itfgmNASA
!!itanm=itanmNASA
!!itm=itmNASA

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

mtm=itm
mtanm=itanm
mtfgm=itfgm

solarangles=rcpara%miss_val
do i=1,eramax
  write(cstatnr,'(I6.6)') wmonrs(index(i))

  filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbinsaveera'//cstatnr
  write(*,*) 'ERA: ',filename
  call read_sonde(iunit,filename,rcpara,tm,tanm,tfgm,tbcm,solarangles,rtype,eerr) !in file rfcorio.f90 line 373
  write(*,*) igcstatnr,'ERA-40',count(tm .ne. rcpara%miss_val),count(tfgm .ne. rcpara%miss_val)

!! calculate analysis adjustment oper/era40
  filename=trim(rcpara%prefix)//cstatnr//'/feedbackx'//cstatnr
  call read_sonde_oper(iunit,filename,rcpara,wmolats(index(i)),wmolons(index(i)),tm,tanm,tfgm,tbcm,solarangles,ominuse40_an,adjust_an,rtype,eerr2) !in file rfcorio.f90 line 611
  filename=trim(rcpara%prefix)//cstatnr//'/feedbackglobbinsaveoper'//cstatnr
  call read_sonde_oper(iunit,filename,rcpara,wmolats(index(i)),wmolons(index(i)),tm,tanm,tfgm,tbcm,solarangles,ominuse40,adjust,rtype,eerr2) !in file rfcorio.f90 line 611
  write(*,*) igcstatnr,'ERA-40',count(tm .ne. rcpara%miss_val),count(tfgm .ne. rcpara%miss_val)

  eerr=eerr*eerr2

  if(eerr .eq. 0) then


    if(priority .eq. 'IGRA') then
      !! merge IGRA and ERA-40
      !! rule: IGRA has precedence over IGRA 
  
      if(needed(index(i)) .eq. 0) then
        needed(index(i))=maxval(needed)+1
        write(*,*) 'needed increased to',needed(index(i))
        where(mtm .eq. rcpara%miss_val .and. tm .ne. rcpara%miss_val .and. .not. merged)
          mtm=tm
          mtanm=tanm
          mtfgm=tfgm
!!      used(:,:,:,needed(index(i)))=.true.
          merged=.true. 
        endwhere
      else
        write(*,*) 'ERA-40 station',wmonrs(index(i)),' matches at least twice - used only once'
 
       endif
     else
       !! merge IGRA and ERA-40
      !! rule: ERA-40 has precedence over IGRA
      !!       If 2 or more ERA-40 data available use nearest one.

      if(needed(index(i)) .eq. 0) then
        needed(index(i))=maxval(needed)+1
        write(*,*) 'needed increased to',needed(index(i))
        where(tm .ne. rcpara%miss_val .and. .not. merged)
          mtm=tm
          mtanm=tanm
          mtfgm=tfgm
!!      used(:,:,:,needed(index(i)))=.true.
          merged=.true. 
        endwhere
      else
        write(*,*) 'ERA-40 station',wmonrs(index(i)),' matches at least twice - used only once'
     endif
   endif
  endif
  
  write(*,*) igcstatnr, count(tm .ne. rcpara%miss_val),' values in ',cstatnr

 enddo

  write(*,*) igcstatnr, count(merged),' values merged in from ERA data'
  write(*,*) igcstatnr, ' original:',count(itm .ne. rcpara%miss_val),' merged:',count(mtm .ne. rcpara%miss_val)

filename=trim(rcpara%prefix)//igcstatnr//'/feedbackmergedmonNASA'//igcstatnr
  tbcm=rcpara%miss_val
  eracorrs=0.
  do ipar=1,rcpara%parmax
    do ip=1,rcpara%pmax
       eracorrs(1:16072,ip,ipar)=adjust(ip,ipar)
    enddo
  enddo
  where(mtfgm .ne. rcpara%miss_val .and. mtm .ne. rcpara%miss_val)
    tbcm=mtm+mtfgm+eracorrs
  endwhere
  do ipar=1,rcpara%parmax
    do ip=1,rcpara%pmax
       eracorrs(1:16072,ip,ipar)=adjust_an(ip,ipar)
    enddo
  enddo
  where(mtanm .ne. rcpara%miss_val .and. mtm .ne. rcpara%miss_val)
    tanm=mtm+mtanm+eracorrs
  endwhere

!! do some final quality control
write(*,*) igcstatnr,' obs-bg values above 20:',count(abs(mtfgm) .gt. 20. .and. mtfgm .ne. rcpara%miss_val)
write(*,*) igcstatnr,' obs-bg values above 30:',count(abs(mtfgm) .gt. 30. .and. mtfgm .ne. rcpara%miss_val)

where(abs(mtfgm) .gt. 20. .and. mtfgm .ne. rcpara%miss_val)
  mtm=rcpara%miss_val
  tanm=rcpara%miss_val
  tbcm=rcpara%miss_val
  mtanm=rcpara%miss_val
  mtfgm=rcpara%miss_val
  itfg12m=rcpara%miss_val
endwhere

call write_sonde_monthly(iunit,filename,rcpara,mtm,err,tanm,tbcm) !in file rfcorio.f90 line 1823

!!tm=mtm ! mtm destroyed in call to write_sonde_daily but needed in write_sonde_monthly
!!tanm=mtanm
!!tfgm=mtfgm
tbcm=itfg12m
filename=trim(rcpara%prefix)//igcstatnr//'/feedbackmergedNASA'//igcstatnr
where(mtm .ne. rcpara%miss_val) 
  itfg12m=itmNASA-mtm
elsewhere
  itfg12m=-5.
endwhere

call write_igrasonde_daily_new(iunit,filename,rcpara%nmax,rcpara%pmax,rcpara%parmax,mtm,mtanm,mtfgm,itfg12m,rcpara%miss_val,err) !in file rfcorio.f90 line 1730

if(wmostats .ne. iwmostats) then
!! IGRA
  filename=trim(rcpara%prefix)//igcstatnr//'/igracorrmonNASA'//igcstatnr
  call write_sonde_monthly(iunit,filename,rcpara,itm,err) !in file rfcorio.f90 line 1823
endif

return
end


subroutine read_igra_NASA(iunit,cstatnr,rcpara,ts) !

use rfmod

implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer imax,first,last,flagmax,iret,statmax,i,l,k,cl,it
integer curmon,curyear,curday,curhour,numrec,relset,index,iunit,igstatnr
integer lgood,lbad,ini,mlvls,lstats,iex,indexmax,pm,tm,istart,wmostats,il
integer iwmo,irem,mdays,ipres,iheight,itemp
integer NGI,NGJ,ios,al
integer pres(rcpara%pmax),height(rcpara%pmax),temp(rcpara%pmax)
real(kind=JPRM) ts(rcpara%nmax,rcpara%pmax,1,rcpara%parmax)
real(kind=JPRM) wmolons(rcpara%statmax),wmolats(rcpara%statmax)
character*20 :: wmonames(rcpara%statmax)
character*36 zeile
CHARACTER*30 filename
character*6 cstatnr
character*4 oper

integer err
integer wmonrs(rcpara%statmax)

REAL(KIND=JPRM) MPLEVS(rcpara%pmax)
logical ex

first=1
index=1
indexmax=rcpara%nmax

filename='../igramerged/'//cstatnr//'.dat'
  open(iunit,file=trim(filename),form='formatted',err=102)

ts=-999.

MPLEVS=rcpara%plevs*100.

l=0
lgood=0
lbad=0
ini=1
do while(.true.)
10 continue

    read(iunit,'(A24)',end=20) zeile

    if(zeile(1:1) .eq. '#') then
!!$ call omp_set_lock(omp_lp)
      read(zeile,'(1X,I5,I4,3I2,I4,I4)',iostat=ios) igstatnr,curyear,curmon,curday,curhour,relset,numrec
      if(ios .ne. 0) then
         print*,ios,zeile
      endif
!!$ call omp_unset_lock(omp_lp)
      mlvls=0
      do i=1,numrec
        read(iunit,'(A36)') zeile
        if(zeile(1:1) .eq. '1' .and. curyear .ge. 1957) then
!!$ call omp_set_lock(omp_lp)
          read(zeile,'(2X,I6,3(1X,I5))',iostat=ios) ipres,iheight,itemp
!!$ call omp_unset_lock(omp_lp)
      if(ios .ne. 0) then
!!$ call omp_set_lock(omp_lp)
         print*,ios,zeile
!!$ call omp_unset_lock(omp_lp)
      endif
          if(any(mplevs .eq. ipres)) then
            mlvls=mlvls+1
            pres(mlvls)=ipres
            height(mlvls)=iheight
            temp(mlvls)=itemp
          endif
        else
!!$ call omp_set_lock(omp_lp)
          if(zeile(1:1) .eq. '#') write(*,*) 'error reading ',filename
!!$ call omp_unset_lock(omp_lp)
        endif
      enddo
    else
!!$ call omp_set_lock(omp_lp)
      write(*,*) 'error reading ',filename
!!$ call omp_unset_lock(omp_lp)
    endif

    if(curyear .lt. 1957 .or. curyear .gt. 2008) goto 10
  
    do while(index .lt. rcpara%nmax .and. rcpara%year(index) .lt. curyear)
      index=index+1
    enddo
    do while(index .lt. rcpara%nmax .and. rcpara%month(index) .lt. curmon)
      index=index+1
    enddo
    do while(index .lt. rcpara%nmax .and. rcpara%day(index) .lt. curday)
      index=index+1
    enddo

  if(curhour .lt. 3 .or. curhour .ge. 21) then
     it=1
     if(curhour .ge. 21) then
       index=index+1 !! launch valid next day
     endif 
    else if   (curhour .ge. 3 .and. curhour .lt. 9) then
     it=2
    else if     (curhour .ge. 9 .and. curhour .lt. 15) then
     it=3
    else if    (curhour .ge. 15 .and. curhour .lt. 21) then
     it=4
    else
       goto 10
    endif

    if(index .gt. rcpara%nmax) goto 20

    if( (it .eq. 2 .or. it .eq. 4).and. rcpara%parmax.eq. 2) goto 10
    if(it.eq.3 .and. rcpara%parmax .eq. 2) it=2
    
!!    write(*,*) index,it

    irem=1
    do cl=mlvls,1,-1
      do i=irem,rcpara%pmax
        if(pres(cl).eq.mplevs(i)) then
          irem=i+1
          if(temp(cl) .ne. -9999 .and. temp(cl).ne.-8888) then
!!            if(qt(cl) .eq. 1 .or. qt(cl) .eq. 4 .or. qt(cl) .eq. 5 ) then 
              ts(index,i,1,it)=temp(cl)/10.+273.2
              lgood=lgood+1
          else
              ts(index,i,1,it)=-999.
              lbad =lbad +1
          endif
        endif
      enddo       
    enddo
!!    write(*,*) lgood,lbad  
  l=l+1
  index=index-1 !! to be safe for launches after 21h where we have increased index by 1
!!$ call omp_set_lock(omp_lp)
  if(mod(l,1000) .eq. 0) print*,l,'records read'
!!$ call omp_unset_lock(omp_lp)
enddo

20 continue
!!$ call omp_set_lock(omp_lp)
print*,l,'records read'
write(*,*) 'statnr: ',cstatnr,',',lgood,' good values',lbad,'bad or unchecked values'
!!$ call omp_unset_lock(omp_lp)
close(iunit)



     return
102  write(*,*) ' Could not open NASA file' , filename
     
     return
end subroutine read_igra_NASA
