      MODULE RWGRIB2

      CONTAINS

      SUBROUTINE READLATLON(FILENAME,FELD,MAXL,MAXB,MLEVEL,MSTRIDE,MPAR,LOCK)

      USE GRIB_API

      IMPLICIT NONE

        integer                            ::  ifile
        integer                            ::  iret
   integer                            ::  n,mk
   integer                            ::  i
   integer,dimension(:),allocatable   ::  igrib
   integer                            ::  numberOfPointsAlongAParallel
   integer                            ::  numberOfPointsAlongAMeridian
   real, dimension(:), allocatable    ::  values
   integer                            ::  numberOfValues
   real,dimension(maxl,maxb,mlevel)   ::  feld  
   integer,intent(in),optional        :: lock
integer:: maxl,maxb,mlevel,mstride,mpar(mstride)  ,ioffset,irest  
character*(*):: filename                             

!!$OMP CRITICAL
if(present(lock)) then
!$ call omp_set_lock(lock)
endif
   write(*,*) filename
   call grib_open_file(ifile, TRIM(FILENAME),'r')
 
   ! count the messages in the file
   call grib_count_in_file(ifile,n)
   allocate(igrib(n))
   igrib=-1
 
   ! Load the messages from the file.
   DO i=1,n
      call grib_new_from_file(ifile,igrib(i), iret)
   END DO
   ! we can close the file
   call grib_close_file(ifile)
if(present(lock)) then
!$ call omp_unset_lock(lock)
endif
!!$OMP END CRITICAL

   ! Loop on all the messages in memory
   DO i=1,n
!      write(*,*) 'processing message number ',i
      !     get as a integer
      call grib_get(igrib(i),'numberOfPointsAlongAParallel', &
           numberOfPointsAlongAParallel)
 
      !     get as a integer
      call grib_get(igrib(i),'numberOfPointsAlongAMeridian', &
           numberOfPointsAlongAMeridian)

      call grib_get(igrib(i),'numberOfVerticalCoordinateValues',mk)

      call grib_get_size(igrib(i),'values',numberOfValues)
!      write(*,*) 'numberOfValues=',numberOfValues
 
      allocate(values(numberOfValues), stat=iret)
      !     get data values
      call grib_get(igrib(i),'values',values)

      IOFFSET=mod(i-1,MSTRIDE)*(mk/2-1)
      feld(:,:,IOFFSET+(i-1)/MSTRIDE+1)=reshape(values,(/maxl,maxb/))

   END DO
 
   DO i=1,n
     call grib_release(igrib(i))
   END DO
 
   deallocate(values)
   deallocate(igrib)

      END SUBROUTINE READLATLON

      SUBROUTINE WRITELATLON(iunit,igrib,ogrib,FELD,MAXL,MAXB,MLEVEL,&
      LEVMIN,LEVMAX,MSTRIDE,MPAR)

      USE GRIB_API

      IMPLICIT NONE

      INTEGER IFIELD,MLEVEL,MNAUF,I,J,K,L,MSTRIDE,IERR,JOUT
      INTEGER MPAR(MSTRIDE),MAXL,MAXB,LEVMIN,LEVMAX
      INTEGER IUNIT,igrib,ogrib
      REAL ZSEC4(MAXL*MAXB)
      REAL    FELD(MAXL,MAXB,MLEVEL)


      DO k=LEVMIN,LEVMAX
        call grib_set(igrib,"level",k)
        DO j=1,MSTRIDE
         call grib_set(igrib,"paramId",MPAR(j))
         zsec4(1:maxl*maxb)=RESHAPE(FELD(:,:,k),(/maxl*maxb/))
         call grib_set(igrib,"values",zsec4)

         call grib_write(igrib,iunit)

        ENDDO
      ENDDO



      END SUBROUTINE WRITELATLON

      SUBROUTINE READSPECTRAL(FILENAME,CXMN,MNAUF,MLEVEL,&
        MAXLEV,MSTRIDE,MPAR,A,B)

      USE GRIB_API

      IMPLICIT NONE


        integer                            ::  ifile
   integer                            ::  iret
   integer                            ::  n,mk
   integer                            ::  i,j
   integer,dimension(:),allocatable   ::  igrib
   real, dimension(:), allocatable    ::  values
   integer                            ::  numberOfValues,maxlev
   REAL :: A(MAXLEV+1),B(MAXLEV+1),pv(2*MAXLEV+2)
   REAL:: CXMN(0:(MNAUF+1)*(MNAUF+2)-1,MLEVEL)
integer:: maxl,maxb,mlevel,mstride,mpar(mstride),mnauf,ioffset
character*(*):: filename                             
 
   call grib_open_file(ifile, TRIM(FILENAME),'r')
 
   ! count the messages in the file
   call grib_count_in_file(ifile,n)
   allocate(igrib(n))
   igrib=-1
 
   ! Load the messages from the file.
   DO i=1,n
      call grib_new_from_file(ifile,igrib(i), iret)
   END DO
 
   ! we can close the file
   call grib_close_file(ifile)
 
   ! Loop on all the messages in memory
   DO i=1,n
      write(*,*) 'processing message number ',i
      !     get as a integer
      call grib_get(igrib(i),'pentagonalResolutionParameterJ', j)

      call grib_get_size(igrib(i),'values',numberOfValues)
      write(*,*) 'numberOfValues=',numberOfValues
 
      call grib_get(igrib(i),'numberOfVerticalCoordinateValues',mk)

      call grib_get(igrib(i),'pv',pv)

      allocate(values(numberOfValues), stat=iret)
      !     get data values
      call grib_get(igrib(i),'values',values)

      IOFFSET=mod(i-1,MSTRIDE)*(mk/2-1)
           CXMN(:,IOFFSET+(i-1)/MSTRIDE+1)=values(1:(MNAUF+1)*(MNAUF+2))

   END DO
 
   DO i=1,n
     call grib_release(igrib(i))
   END DO
 
   deallocate(values)
   deallocate(igrib)



        A=pv(1:1+MAXLEV)
        B=pv(2+MAXLEV:2*MAXLEV+1)

      END SUBROUTINE READSPECTRAL

      END MODULE RWGRIB2
