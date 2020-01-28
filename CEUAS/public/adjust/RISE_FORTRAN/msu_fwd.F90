module tosat
  !
  ! This module is based on RTTOVS example_fwd.F90
  ! It is intended for converting RICH/RAOBCORE
  ! radiosonde pressure level records
  ! into MSU or AMSU brightness temperatures.
  ! 
  ! Leo Haimberger, 30 June 2011
  !
  !
  ! Copyright:
  !    This software was developed within the context of
  !    the EUMETSAT Satellite Application Facility on
  !    Numerical Weather Prediction (NWP SAF), under the
  !    Cooperation Agreement dated 25 November 1998, between
  !    EUMETSAT and the Met Office, UK, by one or more partners
  !    within the NWP SAF. The partners in the NWP SAF are
  !    the Met Office, ECMWF, KNMI and MeteoFrance.
  !
  !    Copyright 2010, EUMETSAT, All Rights Reserved.
  !
  !     *************************************************************
  !
  !     TEST PROGRAM FOR RTTOV SUITE FORWARD MODEL ONLY
  !          RTTOV VERSION 10
  ! To run this program you must have the following files
  ! either resident in the same directory or set up as a
  ! symbolic link:
  !   prof.dat                      --  input profile
  !   rtcoef_platform_id_sensor.dat --  coefficient file to match
  !   the sensor you request in the input dialogue
  ! The script run_example_fwd.sh may be used to run this program.
  ! The output is generated in a file called example_fwd_output.dat.
  !
  !
  ! If the user wants to use this example to create his own
  ! program he will have to modify the code between
  ! comment lines of that kind:
  !     !================================
  !     !======Read =====start===========
  !          code to be modified
  !     !======Read ===== end ===========
  !     !================================
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date        Comment
  ! -------   ----        -------
  !  1.0    27/04/2004   orginal (based on tstrad) P. Brunel
  !  1.1    09/08/2004   modified to allow for variable no. channels/per profile
  !                       R. Saunders
  !  1.2    13/04/2007   Modified for RTTOV-90
  !  1.3    31/07/2007   Modified for RTTOV-91 R Saunders
  !  1.4    11/10/2007   Parallel version P.Marguinaud
  !  2.0    25/06/2010   Modified for RTTOV-10 J Hocking
  !
  ! Code Description:
  !   Language:          Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !     Documenting Exchangeable Fortran 90 Code".
  !
  Use rttov_const, Only :  &
       & errorstatus_success,&
       & errorstatus_warning,&
       & errorstatus_fatal  ,& 
       & q_mixratio_to_ppmv

  Use rttov_types, Only :   &
       & rttov_options,     &
       & rttov_coefs,       &
       & profile_Type,      &
       & transmission_Type, &
       & radiance_Type,     &
       & rttov_chanprof

  Use parkind1, Only : jpim, jprb, jplm,jprm
  !
  Implicit None
  !
!#ifdef _RTTOV_EXAMPLE_FWD_PARALLEL
#include "rttov_parallel_direct.interface"
#define rttov_direct rttov_parallel_direct
!#else
!#include "rttov_direct.interface"
!#endif
#include "rttov_setup.interface"
#include "rttov_copy_prof.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_alloc_rad.interface"
#include "rttov_alloc_transmission.interface"
#include "rttov_alloc_prof.interface"
#include "rttov_errorreport.interface"


!!$  Integer (Kind=jpim) :: nprof            ! Number of profiles per call
  Integer(Kind=jpim)  :: nrttovid          ! Number of sensor coeff files to read (we use only one sensor)
  Integer(Kind=jpim)  :: nchannels         ! Number channels (we use MSU 2-4)

  ! RTTOV_errorhandling interface
  !====================
  Integer(Kind=jpim) :: Err_Unit        ! Logical error unit (<0 for default)
  Integer(Kind=jpim) :: verbosity_level ! (<0 for default)

  ! RTTOV_setup interface
  !====================
  Integer(Kind=jpim)              :: setup_errorstatus  ! setup return code
  Integer(Kind=jpim), Allocatable :: instrument(:,:)    ! platform id, sat id and sensor id
  Integer(Kind=jpim), Allocatable :: msumask(:,:)    ! pressure levels used for msu bt calculation
  Type(rttov_options)             :: opts               ! Options structure
  Type(rttov_coefs), Allocatable  :: coefs(:)           ! Coefficients structure
! ecskin contains surface input parameters for surface emission at RS stations
type::ecskin
  real(kind=JPRM),allocatable :: stnr(:),lats(:),lons(:),months(:),hours(:)
  real(kind=JPRM),allocatable :: tsk(:,:,:)
  real(kind=JPRM),allocatable :: elev(:)
  real(kind=JPRM),allocatable :: lsm(:)
  real(kind=JPRM),allocatable :: ci(:,:,:)
  real(kind=JPRM),allocatable :: t2(:,:,:) 
  real(kind=JPRM),allocatable :: p(:,:,:) 
end type ecskin 
    
  ! RTTOV interface
  !====================
!!  Integer(Kind=jpim)                :: nchannels, nchanprof
!!$  Logical(Kind=jplm), Allocatable   :: calcemis(:)
!!$  Real(Kind=jprb), Allocatable      :: emissivity_in (:)
!!$  Real(Kind=jprb), Allocatable      :: emissivity_out (:)
!!$  Type(transmission_Type)           :: transmission ! transmittances and layer optical depths
!!$  Type(radiance_Type)               :: radiance

!  Integer(Kind=jpim) :: alloc_status(20)
!  Integer(Kind=jpim) :: errorstatus
!  Real(Kind=jprb), Allocatable :: emissivity(:)
!  Character (len=80) :: errMessage

  ! variables for input
  !====================
  Integer(Kind=jpim), Parameter :: mxchn = 90 ! max number of channels
  Integer(Kind=jpim) :: input_chan(mxchn)
  Real(Kind=jprb)    :: input_ems(mxchn)
!!$  Real(Kind=jprb)    :: zenith
!!$  Real(Kind=jprb)    :: azimut
!!$  Real(Kind=jprb)    :: zerht
!!$  Real(Kind=jprb)    :: sunzang
!!$  Real(Kind=jprb)    :: sunazang
!!$!  Integer(Kind=jpim) :: nlevels
!!$  Integer(Kind=jpim) :: ivch, ich
!!$  Integer(Kind=jpim) :: asw          ! allocate or deallocate switch
!!$  Real(Kind=jprb)    :: ems_val
!!$  Integer(Kind=jpim) :: isurf
!!$  Integer(Kind=jpim) :: nwater
!!$  ! loop variables

!!  type(ecskin) :: skin

  !- End of header --------------------------------------------------------
  
contains

!ecskin_init 
subroutine ecskin_init(skin,starty,endy,startdate,mmax,crut2, path, &
                       suff,miss_val,teclimgs)
use netcdf

type(ecskin) :: skin

INTEGER :: ncid,  err, status, ndims !Anzahl der Dimensionen
INTEGER :: nvars, natts !Anzahl der Variablen, Attribute
INTEGER :: ii, jj, kk, i,j,numdat, presslen, hourlen, vals,l,bi,ip,ipar !Anzahl der Werte (zeitlich)
INTEGER :: starty,endy,mmax,raobstarty,startdate
REAL(KIND=JPRM), ALLOCATABLE  :: mhilf(:),thilf(:,:,:)
INTEGER, ALLOCATABLE ::  dimlen(:) !laenge der Dimensionen
INTEGER, ALLOCATABLE :: dimids(:,:) !Vektor der Dimension ids fuer jede Variable
INTEGER, ALLOCATABLE ::  dimnum(:) !Anzahl der Dimensionen einer Variablen
INTEGER, ALLOCATABLE :: datum(:,:)
CHARACTER*50, ALLOCATABLE :: dimnames(:) !Namen der Dimensionen
CHARACTER*50, ALLOCATABLE :: varnames(:) !Variablennamen
CHARACTER*4 :: cstarty,cendy
CHARACTER*(*) :: path  !,cstartdate
INTEGER :: vartype,k !Variablentyp
REAL(KIND=JPRM) :: miss_val,fak,sum
REAL(KIND=JPRM),intent(in) ::crut2(72,36,mmax)
REAL(KIND=JPRM),intent(in),optional ::teclimgs(:,:,:)
character*(*) :: suff

write(cstarty,'(I4)') starty
write(cendy,'(I4)') endy
!read(cstartdate,'(I8)') startdate
100 status=NF90_OPEN(TRIM(path//'TSK'//cstarty//cendy//suff//'.nc'),NF90_NOWRITE,ncid) !nc-file oeffnen
IF (status /= NF90_NOERR) THEN
	WRITE(*,*) 'Error opening netcdf!', TRIM(path//'TSK'//cstarty//cendy//suff//'.nc')
        STOP
	RETURN
END IF
status=NF90_Inquire(ncid,ndims,nvars,natts) !Anzahl der Dimensionen, Variablen und Attributen abfragen
IF (status /= NF90_NOERR) THEN
	WRITE(*,*) 'Error Inquire File ', TRIM(path//'TSK'//cstarty//cendy//suff//'.nc')
	status=NF90_CLOSE(ncid)
        STOP
	RETURN
END IF
ALLOCATE (dimnames(ndims), varnames(nvars), dimlen(ndims), dimnum(nvars), dimids(nvars,NF90_MAX_VAR_DIMS))
numdat=1
DO ii=1,ndims
	status=NF90_Inquire_Dimension(ncid,ii, dimnames(ii), dimlen(ii))
	IF (status /= NF90_NOERR) WRITE(*,*) 'Error Inquire Dimension', ii, 'in File', TRIM(path//'TSK'//cstarty//cendy//suff//'.nc')
END DO


! buffers for netcdf content
ALLOCATE(mhilf(dimlen(1)),thilf(dimlen(1),dimlen(2),dimlen(3)))
! skin record consistent with rcpara record
raobstarty=startdate/10000
if(mmax .ne. (starty-raobstarty)*12+dimlen(1)) then
   write(*,*) 'Inconsistent skin temperature files'
   write(*,*) starty,endy,raobstarty,mmax,dimlen
   stop
endif
ALLOCATE(skin%lats(dimlen(3)),skin%lons(dimlen(3)),&
         skin%stnr(dimlen(3)),skin%months(mmax), &
         skin%hours(dimlen(2)),&
         skin%tsk(mmax,dimlen(2),dimlen(3)),&
         skin%ci(mmax,dimlen(2),dimlen(3)), &
         skin%t2(mmax,dimlen(2),dimlen(3)), &
         skin%p(mmax,dimlen(2),dimlen(3)), &
         skin%elev(dimlen(3)),skin%lsm(dimlen(3)))

DO jj=1,nvars
	status=NF90_Inquire_Variable(ncid,jj, varnames(jj), vartype, dimnum(jj), dimids(jj,:))
	IF (status /= NF90_NOERR) WRITE(*,*) 'Error Inquire Variable', jj, 'in File', TRIM(path//'TSK'//cstarty//cendy//suff//'.nc')
        if(jj .eq. 1) status=NF90_GET_VAR(ncid,jj,mhilf)
        if(jj .eq. 2) status=NF90_GET_VAR(ncid,jj,skin%hours)
        if(jj .eq. 3) status=NF90_GET_VAR(ncid,jj,skin%stnr)
        if(jj .eq. 4) status=NF90_GET_VAR(ncid,jj,skin%lons)
        if(jj .eq. 5) status=NF90_GET_VAR(ncid,jj,skin%lats)
        if(jj .eq. 6) status=NF90_GET_VAR(ncid,jj,thilf)
	
        IF (status /= NF90_NOERR) WRITE(*,*) &
 'Error getting Variable rasocorr from File', TRIM(path//'TSK'//cstarty//cendy//suff//'.nc')

END DO
status=NF90_CLOSE(ncid)
skin%months((starty-raobstarty)*12+1:mmax)=mhilf
do i=1,(starty-raobstarty)*12
  skin%months(i)=(raobstarty+(i-1)/12)*100+mod(i-1,12)+1
enddo
skin%tsk((starty-raobstarty)*12+1:mmax,:,:)=thilf
skin%tsk(1:(starty-raobstarty)*12,:,:)=miss_val

status=NF90_OPEN(TRIM(path//'CI__'//cstarty//cendy//suff//'.nc'),NF90_NOWRITE,ncid) !nc-file oeffnen
 status=NF90_GET_VAR(ncid,6,thilf)
status=NF90_CLOSE(ncid)
skin%ci((starty-raobstarty)*12+1:mmax,:,:)=thilf
skin%ci(1:(starty-raobstarty)*12,:,:)=miss_val
!status=NF90_OPEN(TRIM(path//'T2__'//cstarty//cendy//suff//'.nc'),NF90_NOWRITE,n!cid) !nc-file oeffnen
! status=NF90_GET_VAR(ncid,6,thilf)
!status=NF90_CLOSE(ncid)
!skin%t2((starty-raobstarty)*12+1:mmax,:,:)=thilf
!skin%t2(1:(starty-raobstarty)*12,:,:)=miss_val
status=NF90_OPEN(TRIM(path//'LSP_'//cstarty//cendy//suff//'.nc'),NF90_NOWRITE,ncid) !nc-file oeffnen
 status=NF90_GET_VAR(ncid,6,thilf)
status=NF90_CLOSE(ncid)
skin%p((starty-raobstarty)*12+1:mmax,:,:)=thilf
skin%p(1:(starty-raobstarty)*12,:,:)=miss_val
status=NF90_OPEN(TRIM(path//'ELEV2000010100'//suff//'.nc'),NF90_NOWRITE,ncid) !nc-file oeffnen
 status=NF90_GET_VAR(ncid,2,skin%elev)
status=NF90_CLOSE(ncid)
status=NF90_OPEN(TRIM(path//'LSM_2000010100'//suff//'.nc'),NF90_NOWRITE,ncid) !nc-file oeffnen
 status=NF90_GET_VAR(ncid,2,skin%lsm)
status=NF90_CLOSE(ncid)

! replace tskin, t2 with CRU temperatures
fak=dimlen(2)/72.
do j=1,dimlen(3)
  sum=0.
  k=0
  do i=1,mmax
    if(skin%tsk(i,1,j).ne. miss_val) then
      sum=sum+skin%tsk(i,1,j)
      k=k+1
    endif
    if(skin%tsk(i,2,j).ne. miss_val) then
      sum=sum+skin%tsk(i,2,j)
      k=k+1
    endif
  enddo
  if(k>0) then
    sum=sum/k
  endif

  do i=(starty-raobstarty)*12+1,mmax
    if(crut2(floor((skin%lons(j)+180.)/5)+1,floor((90.+skin%lats(j))/5)+1,i) .ne. -999.) then
      skin%tsk(i,:,j)=crut2(floor((skin%lons(j)+180.)/5)+1,floor((90.+skin%lats(j))/5)+1,i)+sum
      skin%t2(i,:,j)=skin%tsk(i,:,j)
    else
      skin%tsk(i,:,j)=miss_val
      skin%t2(i,:,j)=miss_val
    endif
enddo
enddo

DEALLOCATE(mhilf,thilf)
DEALLOCATE(dimids,dimnum,dimlen,varnames,dimnames)

return
end subroutine ecskin_init

subroutine msu_fwd_init(instrument,isurf,nwater,zenith,input_chan,input_ems,nplev)

  use parkind1

  Integer(Kind=jpim) :: instrument(:,:)    ! platform id, sat id and sensor id
  Integer(Kind=jpim) :: input_chan(mxchn)
  Real(Kind=jprb)    :: input_ems(mxchn)
  Integer(Kind=jpim) :: isurf,nwater 
  Real(Kind=jprm)    :: zenith
  Integer(Kind=jpim) :: j

  Integer(Kind=jpim) :: alloc_status(20)
  Integer(Kind=jpim) :: asw
  Integer(Kind=jpim) :: errorstatus
  Character (len=80) :: errMessage
  Integer(Kind=jpim) :: nplev

  errorstatus     = 0_jpim
!  alloc_status(:) = 0_jpim
!  allocate (instrument(3,nrttovid),stat= alloc_status(1))

  !=====================================================
  !========== Interactive inputs == start ==============
!!$  Write(0,*) 'enter platform number'
!!$  Read(*,*) instrument(1,nrttovid)
!!$  Write(0,*) 'enter satellite number '
!!$  Read(*,*) instrument(2,nrttovid)
!!$  Write(0,*) 'enter instrument number'
!!$  Read(*,*) instrument(3,nrttovid)
!!$  Write(0,*) 'enter surface type (0=land, 1=sea, 2=ice/snow)'
!!$  Read(*,*) isurf
!!$  Write(0,*) 'enter water type (0=fresh water, 1=ocean water)'
!!$  Read(*,*) nwater
!!$  Write(0,*) 'enter number of profile levels'
!!$  Read(*,*) nlevels
!!$  Write(0,*) 'enter zenith angle in degrees'
!!$  Read(*,*) zenith
  !

! pressure levels used for MSU BT calculations
  allocate(msumask(nplev,4))
  if(nplev .eq. 14) then 
  msumask(1:nplev,1)=(/2,3,4,5,6,7,8,9,10,11,12,13,14,16/)
  msumask(1:nplev,2)=(/2,3,4,5,6,7,8,9,10,11,12,13,14,16/)
  msumask(1:nplev,3)=(/1,2,3,4,5,6,7,8,9,10,11,12,13,14/)
  msumask(1:nplev,4)=(/1,2,3,4,5,6,7,8,9,10,11,12,13,14/)
  msumask(1:nplev,4)=msumask(1:nplev,4)-0
  endif
  if(nplev .eq. 13) then 
  msumask(1:nplev,1)=(/3,4,5,6,7,8,9,10,11,12,13,14,16/)
  msumask(1:nplev,2)=(/3,4,5,6,7,8,9,10,11,12,13,14,16/)
  msumask(1:nplev,3)=(/2,3,4,5,6,7,8,9,10,11,12,13,14/)
  msumask(1:nplev,4)=(/2,3,4,5,6,7,8,9,10,11,12,13,14/)
  msumask(1:nplev,4)=msumask(1:nplev,4)-0
  endif
  if(nplev .eq. 12) then 
  msumask(1:nplev,1)=(/4,5,6,7,8,9,10,11,12,13,14,16/)
  msumask(1:nplev,2)=(/4,5,6,7,8,9,10,11,12,13,14,16/)
  msumask(1:nplev,3)=(/3,4,5,6,7,8,9,10,11,12,13,14/)
  msumask(1:nplev,4)=(/3,4,5,6,7,8,9,10,11,12,13,14/)
  msumask(1:nplev,4)=msumask(1:nplev,4)-0
  endif
  if(nplev .eq. 11) then 
  msumask(1:nplev,1)=(/4,5,6,7,8,9,10,11,12,13,14/)
  msumask(1:nplev,2)=(/4,5,6,7,8,9,10,11,12,13,14/)
  msumask(1:nplev,3)=(/3,4,5,6,7,8,9,10,11,12,13/)
  msumask(1:nplev,4)=(/3,4,5,6,7,8,9,10,11,12,13/)
!  msumask(1:nplev,4)=msumask(1:nplev,2)
  endif

!  nmsumask=(/12,12,12,12/)

  ! Prescribe other inputs
!!$  azimut = 0._jprb   ! Satellite azimuth angle
!!$  sunzang = 0._jprb  ! solar zenith angle
!!$  sunazang = 0._jprb ! solar azimuth angle
!!  lat = 0._jprb      ! profile latitude 
!!$  zerht = 0._jprb    ! elevation of surface
  !
  !
  ! Initialise options structure
  opts % addrefrac  = .true.
  opts % addinterp  = .true.
  opts % addsolar   = .false.
  opts % addclouds  = .false.
  opts % addaerosl  = .false.
  opts % ozone_data = .False.       ! we have an ozone profile  
  opts % co2_data   = .False.      ! we do not have profiles
  opts % n2o_data   = .False.      ! for any other constituents
  opts % ch4_data   = .False.      !
  opts % co_data    = .False.      !
  opts % clw_data   = .False.      !

  opts%addinterp =  .true.
  opts%do_checkinput = .false. !.true.
 
  !========== Interactive inputs == end ==============
  !===================================================
  
  !Initialise error management with default value for
  !the error unit number and all error message output
  Err_unit = -1
  verbosity_level = 3_jpim
  Call rttov_errorhandling(Err_unit, verbosity_level)

  allocate (coefs(nrttovid),stat= alloc_status(1))
  If( (alloc_status(1) /= 0) ) then
     errorstatus = errorstatus_fatal
     Write( errMessage, '( "mem allocation error for coefs")' )
     Call Rttov_ErrorReport (errorstatus, errMessage, "msu_fwd_init")
     Stop
  End If

  !Read and initialise coefficients
  Call rttov_setup (&
      & setup_errorstatus,  &! out
      & Err_unit,           &! in
      & verbosity_level,    &! in
      & opts,               &! in
      & coefs(nrttovid),    &! out
      & instrument)          ! in
  if( setup_errorstatus /= errorstatus_success ) then
     write ( *,* ) 'rttov_setup fatal error'
     stop
  endif

  ! security if input number of channels is higher than number
  ! stored in coeffs
!  If( nchannels > coefs(nrttovid) % coef % fmv_chn ) Then
!      nchannels = coefs(nrttovid) % coef % fmv_chn
!      nchan(nprof) = coefs(nrttovid) % coef % fmv_chn
!  Endif

  end subroutine msu_fwd_init

  subroutine msu_fwd_calc(profiles,bt,nprof,nlevels)

  use parkind1
  use rfmod

  ! RTTOV interface
  !====================
  Integer(Kind=jpim)                :: rttov_errorstatus  ! rttov error return code
  Logical(Kind=jplm), Allocatable   :: calcemis(:)
  Real(Kind=jprb), Allocatable      :: emissivity_in (:)
  Real(Kind=jprb), Allocatable      :: emissivity_out (:)
  Type(transmission_Type)           :: transmission ! transmittances and layer optical depths
  Type(radiance_Type)               :: radiance
  Type(rttov_chanprof), Allocatable :: chanprof(:)

  Integer(Kind=jpim) :: alloc_status(20)
  Integer(Kind=jpim) :: errorstatus
  Real(Kind=jprb), Allocatable :: emissivity(:)
  Integer(Kind=jpim), Allocatable :: nchan(:) ! number of channels per profile
  Character (len=80) :: errMessage

  Type(profile_Type)   :: profiles(:)
  real(kind=jprm)  :: bt(nprof,nchannels)  
  Integer(Kind=jpim) :: nprof,nlevels,nchanprof
  Integer(Kind=jpim) :: j
  Integer(Kind=jpim) :: asw

  character*20       :: NameOfRoutine="msu_fwd_calc"

  Integer(Kind=jpim) :: jch
  Integer(Kind=jpim) :: np, nch
  Integer(Kind=jpim) :: ilev, nprint
  Integer(Kind=jpim) :: iprof, joff

  nchanprof=nchannels*nprof

  Allocate (nchan(nprof))
  
  nchan(:) = nchannels
  nchanprof=SUM(nchan(:))  ! here it is also nchannels * nprof
  
  !Pack channels and emmissivity arrays
  Allocate(chanprof(nchanprof),stat= alloc_status(1))   ! Note these array sizes nchan can vary per profile
  Allocate( emissivity ( nchanprof ) ,stat= alloc_status(2))
!!$  Allocate(emissivity_out(nchanprof))
  
  ! Build the list of profile indices
  nch = 0_jpim
  Do j = 1 , nprof
    DO  jch = 1,nchan(j)
      nch = nch + 1_jpim
      chanprof(nch)%prof = j
      chanprof(nch)%chan = input_chan(jch)
      emissivity(nch)    = input_ems(jch)
    End Do
  End Do
  

  asw = 1 ! allocate
  call rttov_alloc_rad( &
      & errorstatus,    &
      & nchanprof,      &
      & radiance,       &
      & nlevels-1_jpim,      &
      & asw)
  If( errorstatus /= errorstatus_success) Then
     errorstatus = errorstatus_fatal
     Write( errMessage, '( "mem allocation error for radiance arrays")' )
     Call Rttov_ErrorReport (errorstatus, errMessage, NameOfRoutine)
     Stop
  Endif
  
  Allocate( calcemis ( nchanprof ) ,stat= alloc_status(3))
  Allocate( emissivity_in ( nchanprof ) ,stat= alloc_status(4))
  Allocate( emissivity_out ( nchanprof ) ,stat= alloc_status(5))
  If( Any(alloc_status(1:5) /= 0) ) Then
     errorstatus = errorstatus_fatal
     Write( errMessage, '( "mem allocation error in emissivity arrays")' )
     Call Rttov_ErrorReport (errorstatus, errMessage,  NameOfRoutine)
     Stop
  End If

  ! allocate transmittance structure
  call rttov_alloc_transmission( &
      & errorstatus,             &
      & transmission,            &
      & nlevels-1_jpim,          &
      & nchanprof,               &
      & asw,                     &
      & init = .true._jplm)  
  If( errorstatus /= errorstatus_success) Then
     errorstatus = errorstatus_fatal
     Write( errMessage, '( "mem allocation error for transmission arrays")' )
       Call Rttov_ErrorReport (errorstatus, errMessage, NameOfRoutine)
     Stop
  Endif

  ! save input values of emissivities for all calculations
  ! calculate emissivity where the input emissivity value is less than 0.01
  emissivity_in(:) = emissivity(:)
  calcemis(:) = emissivity(:) < 0.01_JPRB

  ! Call RTTOV forward model
!!$omp critical
!!$ call omp_set_lock(omp_lp(2885))
!  call rttov_direct(         &

  call rttov_parallel_direct(         &
        & rttov_errorstatus, &! out   error flag
        & chanprof,          &! in    channel and profile index structure
        & opts,              &! in    options structure
        & profiles,          &! in    profile array
        & coefs(nrttovid),   &! in    coefficients strucutre
        & calcemis,          &! in    flag for intermal emissivity calc
        & emissivity_in,     &! in    input emissivities per channel
        & emissivity_out,    &! out   emissivities used by RTTOV per channel
        & transmission,      &! inout computed transmittances
        & radiance)           ! inout computed radiances
   
  If ( rttov_errorstatus /= errorstatus_success ) Then
     Write ( 0, * ) 'rttov_direct error'
     Stop
  End If
!!$ call omp_unset_lock(omp_lp(2885))
!!$omp end critical

!!$  ! transfer data to printing arrays
!!$  Allocate(pr_radcld(nchannels), stat= alloc_status(1))
!!$  Allocate(pr_trans(nchannels), stat= alloc_status(2))
!!$  Allocate(pr_emis(nchannels), stat= alloc_status(3))
!!$  Allocate(pr_trans_lev(nlevels,nchannels), stat= alloc_status(4))
!!$  If( Any(alloc_status /= 0) ) Then
!!$     errorstatus = errorstatus_fatal
!!$     Write( errMessage, '( "mem allocation error for printing arrays")' )
!!$     Call Rttov_ErrorReport (errorstatus, errMessage, NameOfRoutine)
!!$     Stop
!!$  End If

  Do iprof = 1, nprof
!!$    pr_radcld(:) = 0.0_JPRB
!!$    pr_trans(:) = 0.0_JPRB
!!$    pr_emis(:) = 0.0_JPRB
!!$    pr_trans_lev(:,:) = 0.0_JPRB
    !
    joff = (iprof-1_jpim) * nchannels
!!$    Do j = 1 , nchannels
!!$      pr_radcld(j) = radiance % cloudy(j+joff)
!!$      pr_trans(j)  = Transmission % tau_total(j+joff)
!!$      pr_emis(j)   = emissivity_out(j+joff )
!!$      Do ilev = 1 ,  nlevels
!!$          pr_trans_lev(ilev,j) = Transmission % tau_levels(ilev,J+joff)
!!$      Enddo
!!$    Enddo
!!$    !
!!$    !     OUTPUT RESULTS
!!$    !
!!$    NPRINT = 1+ Int((nchannels-1)/10)
!!$    Write(*,*)' -----------------'
!!$    Write(*,*)' Instrument ', instrument(3,nrttovid)
!!$    Write(*,*)' -----------------'
!!$    Write(*,*)' '
!!$    Write(*,*)' Profile ',iprof
!!$    Write(*,*)' '
!!$    WRITE(*,777) profiles(iprof)%zenangle,profiles(iprof)%azangle, &
!!$        & csun,profiles(iprof)%sunzenangle,profiles(iprof)%sunazangle,profiles(iprof)%skin%surftype,&
!!$        & profiles(iprof)%skin%watertype,profiles(iprof)%latitude,profiles(iprof)%elevation,cref,caer,ccld,&
!!$        & instrument(2,nrttovid)
!!$  
!       WRITE(*,*)'CHANNELS PROCESSED:'
!       WRITE(*,111) (chanprof(j) % chan, j = 1+joff,nchannels+joff)
!       WRITE (*,*)' '
!       Write(*,222) (radiance % bt(j), j = 1+joff,nchannels+joff)
!       Write(*,*)' '
       bt(iprof,1:nchannels)=radiance % bt(1+joff:nchannels+joff)
!!$    Write(IOOUT,*)'CALCULATED RADIANCES: SAT =', instrument(2,nrttovid)
!!$    Write(IOOUT,222) (radiance % total(j), j = 1+joff,nchannels+joff)
!!$    Write(IOOUT,*)' '
!!$    Write(IOOUT,*)'CALCULATED OVERCAST RADIANCES: SAT =', instrument(2,nrttovid)
!!$    Write(IOOUT,222) (pr_radcld(j), j = 1,nchannels)
!!$    Write (IOOUT,*)' '
!!$    Write(IOOUT,*)'CALCULATED SURFACE TO SPACE TRANSMITTANCE: S'&
!!$              & ,'AT =',instrument(2,nrttovid) 
!!$    Write(IOOUT,4444) (pr_trans(j), j = 1,nchannels)
!!$    Write (IOOUT,*)' '
!!$    Write(IOOUT,*)'CALCULATED SURFACE EMISSIVITIES '&
!!$              & ,'SAT =',instrument(2,nrttovid) 
!!$    Write(IOOUT,444) (pr_emis(j), j = 1,nchannels)
!!$    !
!!$    !
!!$    If(nchan(nprof) <= 20)Then
!!$      Do  NP = 1 , NPRINT
!!$          Write (IOOUT,*)' '
!!$          Write (IOOUT,*)'Level to space transmittances for channels'
!!$          Write(IOOUT,1115) (chanprof(j+joff) % chan,&
!!$                  & J = 1+joff+(NP-1)*10,Min(Int(10+joff+(NP-1)*10,jpim),nchannels)) 
!!$          Do  ILEV = 1 , NLEVELS
!!$            Write(IOOUT,4445)ILEV,(pr_trans_lev(ilev,J),&
!!$                    & J = 1+(NP-1)*10,Min(Int(10+(NP-1)*10,jpim),nchannels)) 
!!$          End Do
!!$          Write(IOOUT,1115) (chanprof(j+ joff) % chan,&
!!$                  & J = 1+joff+(NP-1)*10,Min(Int(10+joff+(NP-1)*10,jpim),nchannels)) 
!!$      End Do
!!$    Endif
  End Do
  !
  ! Deallocate arrays
  deallocate( chanprof,       stat=alloc_status(1))
  deallocate( emissivity,     stat=alloc_status(2))
  deallocate( emissivity_in,  stat=alloc_status(3))
  deallocate( emissivity_out, stat=alloc_status(4))
  deallocate( calcemis,       stat=alloc_status(5))
  deallocate( nchan,       stat=alloc_status(6))

!!$  ! dealloc printing arrays
!!$  deallocate( pr_radcld,      stat=alloc_status(6))
!!$  deallocate( pr_trans,       stat=alloc_status(7))
!!$  deallocate( pr_emis,        stat=alloc_status(8))
!!$  deallocate( pr_trans_lev,   stat=alloc_status(9))
  If( any(alloc_status(1:6) /= 0) ) then
    errorstatus = errorstatus_fatal
    Write( errMessage, '( "mem deallocation error")' )
    Call Rttov_ErrorReport (errorstatus, errMessage, NameOfRoutine)
    Stop
  End If

  asw = 0 ! deallocate radiance arrays
  call rttov_alloc_rad (errorstatus,nchannels,radiance,nlevels-1_jpim,asw)
  If(errorstatus /= errorstatus_success) Then
    Write( errMessage, '( "radiance deallocation error")' )
    Call Rttov_ErrorReport (errorstatus, errMessage, NameOfRoutine)
  Endif

  asw = 0 ! deallocate transmission arrays
  call rttov_alloc_transmission (errorstatus,transmission,nlevels-1_jpim,nchannels,asw)
  If(errorstatus /= errorstatus_success) Then
    Write( errMessage, '( "radiance deallocation error")' )
    Call Rttov_ErrorReport (errorstatus, errMessage, NameOfRoutine)
  Endif
  
  asw = 0 ! deallocate profile arrays
  call rttov_alloc_prof (errorstatus,nprof,profiles,nlevels,opts,asw)
!  deallocate(profiles,stat=alloc_status(1))
  If(errorstatus /= errorstatus_success .or. alloc_status(1) /= 0) Then
    Write( errMessage, '( "profile deallocation error")' )
    Call Rttov_ErrorReport (errorstatus, errMessage, NameOfRoutine)
  Endif


111  FORMAT(1X,10I8)
1115 Format(3X,10I8)
222  Format(1X,10F8.2)
777 FORMAT( &
       & ' ZENITH ANGLE       =',F7.2,/ &
       & ' AZIMUTH ANGLE      =',F7.2,/&
       & ' SOLAR RADIATION    =',A7,/&
       & ' SOLAR ZENITH ANGLE =',F7.2,/&
       & ' SOLAR AZIMUTH ANGLE=',F7.2,/ &
       & ' SURFACE TYPE       =',I7,/&
       & ' WATER TYPE         =',I7,/ &
       & ' LATITUDE           =',F7.2,/&
       & ' ELEVATION          =',F7.2/&
       & ' REFRACTION         =',A7,/&
       & ' AEROSOLS           =',A7,/&
       & ' CLOUDS             =',A7,//,&
       &'CALCULATED BRIGHTNESS TEMPERATURES: SAT =',I2 )

  end subroutine msu_fwd_calc

  subroutine msu_fwd_deallocate

  Integer(Kind=jpim) :: alloc_status(20)
  character*20  :: NameOfRoutine="msu_fwd_deallocate"
  Integer(Kind=jpim) :: errorstatus
  Character (len=80) :: errMessage

  Call rttov_dealloc_coefs (errorstatus, coefs(nrttovid))
  If(errorstatus /= errorstatus_success) Then
    Write( errMessage, '( "coefs deallocation error")' )
    Call Rttov_ErrorReport (errorstatus, errMessage, NameOfRoutine)
  Endif
  Deallocate(coefs,stat=alloc_status(1))
  If(errorstatus /= errorstatus_success .or. alloc_status(1) /= 0) Then
    Write( errMessage, '( "coefficients array deallocation error")' )
    Call Rttov_ErrorReport (errorstatus, errMessage, NameOfRoutine)
  Endif
   

111  FORMAT(1X,10I8)
1115 Format(3X,10I8)
222  Format(1X,10F8.2)
444  Format(1X,10F8.3)
4444 Format(1X,10F8.4)
4445 Format(1X,I2,10F8.4)

End subroutine msu_fwd_deallocate

end module tosat
