! 7.2.2008 reads station .txt files, writes netcdf
! needs inputfile 'filelist.txt' with the stationfilenames
! compile with netcdf libraries
! flags gespeichert in nc als final*1000000+fg*100000+depar*10000+blacklist*1000+reject*100+depar2big*10+fg2big
MODULE var_def
USE RFMOD
USE RFCORIO
IMPLICIT NONE

!include 'mpif.h'

TYPE, PUBLIC :: atts
	CHARACTER*120 :: title
	CHARACTER*16, ALLOCATABLE :: attname(:)
	CHARACTER*60, ALLOCATABLE :: long_name(:)
	CHARACTER*7, ALLOCATABLE :: units(:)
	CHARACTER*120, ALLOCATABLE :: cell_methods(:)
	REAL(KIND=JPRM), ALLOCATABLE :: missing_val(:)
	REAL(KIND=JPRM), ALLOCATABLE :: valid_range(:,:)
	
!	REAL, ALLOCATABLE :: valid_max(:)
END TYPE atts
!REAL :: lat, lon 
!INTEGER, ALLOCATABLE :: datum(:), hours(:,:)

!this is the 20CR station
TYPE  :: twenty_CR
   INTEGER :: varno !type of observation e
                    ! 2 = Upper Air Temperature [k]
                    ! 3 = U component of wind   [m/s]
                    ! 4 = V component of Wind   [m/s]
                    ! 7 = q Specific humidity   [kg/kg]  
                    !29 = RH Relative humidity  [Numeric]
   REAL, ALLOCATABLE , DIMENSION(:):: lat      !latitude as 180 point each 2 degrees (0 - 358)
   REAL, ALLOCATABLE , DIMENSION(:):: lon      !latitude as 91 point each 2 degrees  (90- 0 - -90)
   INTEGER, ALLOCATABLE , DIMENSION(:):: date  !date yyyy as time_index function
   INTEGER, ALLOCATABLE , DIMENSION(:):: hour  !0 6 12 18 as time_index function
   INTEGER, ALLOCATABLE , DIMENSION(:):: pressure_layers !in the 20 CR there are 24 press_lvl
   INTEGER*2, ALLOCATABLE,DIMENSION(:,:,:,:) :: short_data_PL!here I store the values as (lon,lat,level,time_index)function for the pressure level data
   INTEGER*2, ALLOCATABLE,DIMENSION(:,:,:) :: short_data_S!here I store the values as (lon,lat,time_index)function for the surface data
                                                    !they are as short integer that I will convert in REAL
END TYPE twenty_CR

END MODULE var_def

MODULE txtnc
USE var_def

logical :: lserial=.true.

	CONTAINS

SUBROUTINE txt_to_nc(rcpara)
USE netcdf
USE var_def
USE rfmod
IMPLICIT NONE

type(rasocor_namelist) rcpara
REAL (kind=JPRM):: lat, lon 
INTEGER :: obstype, codetype, station
INTEGER :: itemp, ios, ios2, ss
INTEGER :: vals,l,i,filestat,istat
INTEGER, ALLOCATABLE :: datum(:), hours(:,:)
REAL (kind=JPRM) :: temp, biascorr, fg_depar,an_depar
REAL (kind=JPRM),ALLOCATABLE :: arr(:,:,:,:)
REAL (kind=JPRM), DIMENSION(rcpara%nmax,rcpara%pmax,rcpara%parmax) :: temperaturesall, fg_all, bias_all, flags_all,an_all
INTEGER :: datumall(rcpara%nmax)
INTEGER :: hoursall(rcpara%nmax,rcpara%parmax)
!filename variables
CHARACTER :: infile*60,inf*60,user*60,cstatnr*5
CHARACTER :: outfile*80

INTEGER :: alt, datumneu, hour, press  !Pressure in pa
INTEGER :: datumalt, time !zaehlt Datumsindex
INTEGER :: hr !hourindex, kann 1 oder 2 werden
INTEGER :: p !druckniveauindex
INTEGER :: indexcor ! wird -1 bei hour >=21, sonst 0
!INTEGER :: vals
!flag variables from inputfile
INTEGER :: final, fg, depar, reject, fg2big
INTEGER ::  blacklist,depar2big ! analysis departure too big
!gesamtflag
INTEGER :: flagges 
CHARACTER :: zeile*310

call getenv('HOME',user)
!OPEN (1,FILE=trim(user)//'/tables/filelist.txt')
OPEN (1,FILE='filelist.txt')
!read filelist
!start loop over files

istat=1
ios2=0
ss=0
DO WHILE (ios2==0)
ss=ss+1
!	IF (ss .GT. 1438) THEN
	WRITE(*,*) ss
!	END IF
	READ(1,'(A60)',IOSTAT=ios2) infile
	IF (ios2 .NE. 0) EXIT
!create outfilename 
        if(len(trim(infile)) .eq. 10) then 
          inf='0'//infile
  	  outfile= inf(1:LEN(TRIM(inf))-3)//'nc'
        else
          inf=infile
	  outfile= infile(1:LEN(TRIM(infile))-3)//'nc'
        endif
        cstatnr=outfile(LEN(TRIM(inf))-10:LEN(TRIM(inf))-6)
        outfile='/scratch/srvx7/leo/scratch/ei/'//cstatnr//'/'//trim(outfile)
!        outfile='./'//trim(outfile)
        write(*,*) infile,cstatnr,outfile
        read(cstatnr,*) filestat
!Initialisieren
	CALL init(rcpara, datumalt, datumneu, datumall, hoursall, time, hr, p, vals, indexcor, temperaturesall, fg_all, an_all, bias_all, flags_all, lat, lon, alt)
!open file
	OPEN (2,FILE=TRIM(infile))
	ios=0
!read values
        l=0
        datumall=rcpara%miss_val
	DO WHILE (ios==0)
		READ (2,'(A310)', IOSTAT=ios) zeile
		IF (ios .NE. 0) EXIT
                IF(index(zeile,'NULL') .eq. 0) THEN
			IF(index(zeile(149:166),'.') .EQ. 0) THEN
				READ (zeile,FMT=978,err=500) datumneu, hour, obstype, codetype, station, lat, lon, alt, press, itemp, biascorr, fg_depar, an_depar, blacklist, reject, depar, fg, final,  depar2big, fg2big
				temp=REAL(itemp)
			ELSE	
				READ (zeile,FMT=977,err=500) datumneu, hour, obstype, codetype, station, lat, lon, alt, press, temp, biascorr, fg_depar, an_depar, blacklist, reject, depar, fg, final,  depar2big, fg2big
                           
			END IF
		ELSE
		  CYCLE
		ENDIF
		
!		WRITE(*,*)datum, hour, obstype, codetype, station, lat, lon, stalt, press, temp, biascorr, fg_depar, blacklist, varqc, depar, fg, final, varqcper, depar2big
		IF (filestat .eq. station .and. (hour .LE. 3 .OR. (hour .GE. 9 .AND. hour .LE. 15) .OR. hour .GE. 21)) THEN
			CALL calc_ind_flags(rcpara,datumneu, hour, press, time, hr, p, indexcor, blacklist, reject, depar, fg, final, depar2big, fg2big, flagges) 
!Put values into variables
		
			temperaturesall(time,p,hr)=temp
			fg_all(time,p,hr)=fg_depar
			an_all(time,p,hr)=an_depar
			bias_all(time,p,hr)=biascorr
			flags_all(time,p,hr)=flagges
			hoursall(time,hr)=hour
			datumall(time)=time
!convert to small, one dimension array
!CALL dim3todim1(rcpara,tm, temperatures, vals, ind)
		END IF
                goto 501
500             write(*,*) cstatnr,' read error',zeile
501	END DO ! end reading loop
	CLOSE(2)

        l=0
        do i=1,rcpara%nmax
          if(datumall(i) .eq. i) then
            l=l+1
          endif
        enddo
	vals=l
	ALLOCATE(arr(vals,rcpara%pmax,rcpara%parmax,5),datum(vals), hours(vals,rcpara%parmax))
        l=0
        do i=1,rcpara%nmax
          if(datumall(i) .eq. i) then
            l=l+1
	    arr(l,:,:,1)=temperaturesall(i,:,:)
       	    arr(l,:,:,2)=fg_all(i,:,:)
	    arr(l,:,:,3)=bias_all(i,:,:)
	    arr(l,:,:,4)=flags_all(i,:,:)
       	    arr(l,:,:,5)=an_all(i,:,:)
	    datum(l)=datumall(i)
	    hours(l,:)=hoursall(i,:)
          endif
        enddo
        
!write netcdf
	IF (vals .GT. 0)THEN
		CALL txttonc(rcpara,istat,outfile,arr,(/'temperatures','fg_dep','bias','flags','an_dep'/),vals, 5,datum, hours, lat, lon, alt,(/'temperature','fg-departures','bias correction','flags-final,fg,depar, blacklist, andepar2big, fg2big','an-departure'/))
	END IF
	DEALLOCATE (arr, hours, datum)
END DO !End loop over files
	CLOSE(1)
977	FORMAT(2x,I8,1x,I2,16x,I1,11x,I2,8x,I5,2(15x,F7.2),17x,I7,14x,I6,6x,E16.12,3(1x,E21.17),29x,2I1,23x,3I4,22x,2I1)
978	FORMAT(2x,I8,1x,I2,16x,I1,11x,I2,8x,I5,2(15x,F7.2),17x,I7,14x,I6,6x,I16,3(1x,E21.17),29x,2I1,23x,3I4,22x,2I1)
END SUBROUTINE txt_to_nc

! Read files with all vars (T,q,rh,u,v) 
SUBROUTINE txt_to_nc_MERRA(rcpara)
USE netcdf
USE var_def
USE rfmod
IMPLICIT NONE

type(rasocor_namelist) rcpara
REAL (kind=JPRM):: lat, lon 
INTEGER :: obstype, codetype, station
INTEGER :: itemp, ios, ios2, ss
INTEGER :: vals,l,i,filestat,istat
INTEGER, ALLOCATABLE :: datum(:), hours(:,:)
REAL (kind=JPRM) :: temp, biascorr, fg_depar,an_depar,fg, press,final
REAL (kind=JPRM),ALLOCATABLE :: arr(:,:,:,:),stypes(:)
REAL (kind=JPRM), DIMENSION(rcpara%nmax,rcpara%pmax,rcpara%parmax) :: temperaturesall, fg_all, bias_all, flags_all,an_all
REAL (kind=JPRM), DIMENSION(rcpara%nmax) :: sonde_type_all
INTEGER :: datumall(rcpara%nmax)
INTEGER :: hoursall(rcpara%nmax,rcpara%parmax)
!filename variables
CHARACTER :: infile*60,inf*60,user*60,cstatnr*5
CHARACTER :: outfile*80

INTEGER :: alt, datumneu, hour  !Pressure in pa
INTEGER :: datumalt, time !zaehlt Datumsindex
INTEGER :: hr !hourindex, kann 1 oder 2 werden
INTEGER :: p !druckniveauindex
INTEGER :: indexcor ! wird -1 bei hour >=21, sonst 0
!INTEGER :: vals
!flag variables from inputfile
INTEGER :: sonde_type,varno
INTEGER ::  blacklist,depar2big ! analysis departure too big
!gesamtflag
INTEGER :: flagges 
CHARACTER :: zeile*360,csonde_type*4,chilf*32

call getenv('HOME',user)
!OPEN (1,FILE=trim(user)//'/tables/filelist.txt')
!OPEN (1,FILE='filelist.txt')
!read filelist
!start loop over files

call getarg(2,infile)
istat=1
ios2=0
ss=0
!DO WHILE (ios2==0)
ss=ss+1
!	IF (ss .GT. 1438) THEN
	WRITE(*,*) ss
!	END IF
!	READ(1,'(A60)',IOSTAT=ios2) infile
!	IF (ios2 .NE. 0) EXIT
!create outfilename 
        if(len(trim(infile)) .eq. 10) then 
          inf='0'//infile
  	  outfile= inf(1:LEN(TRIM(inf))-3)//'nc'
        else
          inf=infile
	  outfile= infile(1:LEN(TRIM(infile))-3)//'nc'
        endif
        cstatnr=outfile(LEN(TRIM(inf))-10:LEN(TRIM(inf))-6)
        outfile='/vgc/srvx7/leo/scratch/MERRA/'//cstatnr//'/'//trim(outfile)
!        outfile='./'//trim(outfile)
        write(*,*) infile,cstatnr,outfile
        read(cstatnr,*) filestat
!Initialisieren
!	CALL init(rcpara, datumalt, datumneu, datumall, hoursall, time, hr, p, vals, indexcor, temperaturesall, fg_all, an_all, bias_all, flags_all, lat, lon, alt, sonde_type_all)
!open file
        temperaturesall=rcpara%miss_val
        fg_all=rcpara%miss_val
        an_all=rcpara%miss_val
        bias_all=rcpara%miss_val
        sonde_type_all=rcpara%miss_val
        flags_all=0

	OPEN (2,FILE=TRIM(infile),iostat=ios)
	ios=0
!read values
        l=0
        datumall=rcpara%miss_val
	DO WHILE (ios==0)
		READ (2,'(A310)', IOSTAT=ios) zeile
		IF (ios .NE. 0) EXIT
		  READ (zeile,FMT=980) press,datumneu , hour, temp,fg,fg_depar,an_depar, final,biascorr
                     if(hour .gt. 21 .or. hour .le. 3 .or. hour .gt. 9 .and. hour .le. 15) then
                        sonde_type=0
                        time=toindex(datumneu,rcpara)

                        if(hour .gt. 21) then
                          time=time+1
                        endif
		        if(hour .gt. 21 .or. hour .le. 3) then
                          hr=1
                        else
                          hr=2
                        endif
                        do p=1,rcpara%pmax
                          if(press .eq. rcpara%plevs(p)) exit
                        enddo
			temperaturesall(time,p,hr)=temp
			fg_all(time,p,hr)=fg_depar
			an_all(time,p,hr)=an_depar
			sonde_type_all(time)=sonde_type
			bias_all(time,p,hr)=biascorr
			flags_all(time,p,hr)=1-floor(final)
			hoursall(time,hr)=hour
			datumall(time)=time
                     endif
                goto 501
500             write(*,'(a)') cstatnr,' read error',zeile
501	END DO ! end reading loop
	CLOSE(2)
        l=0
        do i=1,rcpara%nmax
          if(datumall(i) .eq. i) then
            l=l+1
          endif
        enddo
	vals=l
	ALLOCATE(arr(vals,rcpara%pmax,rcpara%parmax,5),datum(vals), hours(vals,rcpara%parmax),stypes(vals))
        l=0
        do i=1,rcpara%nmax
          if(datumall(i) .eq. i) then
            l=l+1
	    arr(l,:,:,1)=temperaturesall(i,:,:)
       	    arr(l,:,:,2)=fg_all(i,:,:)
	    arr(l,:,:,3)=bias_all(i,:,:)
	    arr(l,:,:,4)=flags_all(i,:,:)
       	    arr(l,:,:,5)=an_all(i,:,:)
	    datum(l)=datumall(i)
	    hours(l,:)=hoursall(i,:)
            stypes(l)=sonde_type_all(i)
          endif
        enddo
        print*,vals
!write netcdf
	IF (vals .GT. 0)THEN
		CALL txttonc(rcpara,istat,outfile,arr,(/'temperatures','fg_dep','bias','flags','an_dep','s_type'/),vals, 5,datum, hours, lat, lon, alt,(/'temperature','fg-departures','bias correction','flags-final,fg,depar, blacklist, andepar2big, fg2big','an-departure','sonde type'/),stypes)
	END IF
	DEALLOCATE (arr, hours, datum)
!END DO !End loop over files
	CLOSE(1)
!  01001     70.940   351.340   925.000   199912312303  272.450  273.124   -0.674   -0.043     1.000
980	FORMAT(30x,F9.3,2x,I8,I2,2x,4F9.3,F11.3,F7.1)
!980	FORMAT(31x,F9.3,4x,I8,I2,5I9.3)
END SUBROUTINE txt_to_nc_MERRA

! Read files with all vars (T,q,rh,u,v) 
SUBROUTINE txt_to_nc_wtype(rcpara)
USE netcdf
USE var_def
USE rfmod
IMPLICIT NONE

type(rasocor_namelist) rcpara
REAL (kind=JPRM):: lat, lon 
INTEGER :: obstype, codetype, station,fstation
INTEGER :: itemp, ios, ios2, ss
INTEGER :: vals,l,i,filestat,istat
INTEGER, ALLOCATABLE :: datum(:), hours(:,:)
REAL (kind=JPRM) :: temp, biascorr, fg_depar,an_depar
REAL (kind=JPRM),ALLOCATABLE :: arr(:,:,:,:),stypes(:)
REAL (kind=JPRM), DIMENSION(rcpara%nmax,rcpara%pmax,rcpara%parmax) :: temperaturesall, fg_all, bias_all, flags_all,an_all
REAL (kind=JPRM), DIMENSION(rcpara%nmax) :: sonde_type_all
INTEGER :: datumall(rcpara%nmax)
INTEGER :: hoursall(rcpara%nmax,rcpara%parmax)
!filename variables
CHARACTER :: infile*60,inf*60,user*60,cstatnr*6
CHARACTER :: outfile*80

INTEGER :: alt, datumneu, hour, press  !Pressure in pa
INTEGER :: datumalt, time !zaehlt Datumsindex
INTEGER :: hr !hourindex, kann 1 oder 2 werden
INTEGER :: p !druckniveauindex
INTEGER :: indexcor ! wird -1 bei hour >=21, sonst 0
!INTEGER :: vals
!flag variables from inputfile
INTEGER :: final, fg, depar, reject, fg2big,sonde_type,varno,status,event
INTEGER ::  blacklist,depar2big ! analysis departure too big
!gesamtflag
INTEGER :: flagges 
CHARACTER :: zeile*360,csonde_type*4,chilf*32,cstatid*8

call getenv('HOME',user)
!OPEN (1,FILE=trim(user)//'/tables/filelist.txt')
!OPEN (1,FILE='filelist.txt')
!read filelist
!start loop over files

call getarg(2,infile)
istat=1
ios2=0
ss=0
!DO WHILE (ios2==0)
ss=ss+1
!	IF (ss .GT. 1438) THEN
	WRITE(*,*) ss
!	END IF
!	READ(1,'(A60)',IOSTAT=ios2) infile
!	IF (ios2 .NE. 0) EXIT
!create outfilename 
        print*, len(trim(infile))
        if(len(trim(infile)) .eq. 11) then 
          inf='0'//infile
  	  outfile= inf(1:LEN(TRIM(inf))-3)//'nc'
        else
          inf=infile
	  outfile= infile(1:LEN(TRIM(infile))-3)//'nc'
        endif
        
        cstatnr=outfile(1:LEN(TRIM(inf))-6)
        outfile='/vgc/srvx7/leo/scratch/ei6/'//trim(cstatnr)//'/'//trim(outfile)
!        outfile='./'//trim(outfile)
        write(*,*) infile,cstatnr,outfile
        read(cstatnr,*) filestat
!Initialisieren
	CALL init(rcpara, datumalt, datumneu, datumall, hoursall, time, hr, p, vals, indexcor, temperaturesall, fg_all, an_all, bias_all, flags_all, lat, lon, alt, sonde_type_all)
!open file
	OPEN (2,FILE=TRIM(infile),iostat=ios)
        infile='0'//trim(infile)
	ios=0
!read values
        l=0
        datumall=rcpara%miss_val
	DO WHILE (ios==0)
		READ (2,'(A310)', IOSTAT=ios) zeile
		IF (ios .NE. 0) EXIT
                IF(zeile(22:25) .eq. 'NULL') zeile=zeile(1:21)//'0'//trim(zeile(23:300))
                ios=index(zeile,'NULL')
                do while(ios .ne. 0)
                  zeile(ios:ios+3)='-999'
                  ios=index(zeile,'NULL')
                enddo
              
				READ (zeile,*,iostat=ios) datumneu, hour, obstype, codetype, sonde_type, cstatid, lat, lon, alt, press, varno, temp, biascorr, fg_depar, an_depar!, final, reject, event !, fg, final
                           
                if(ios > 0) then
write(*,*) 'could not read',zeile
ios=0
cycle
endif
                if (len(cstatid)<8) then
                   write(*,*) cstatid
                   cycle
                endif
                if (cstatid(2:2).eq.':') then
                  read(cstatid(3:8),*) station
                else
                  read(cstatid(3:8),*,iostat=ios) station
                  if (ios .ne. 0) then
                      write(*,*) cstatid, ios
                  endif
                endif
                read(infile(1:6),*) fstation
                if(fstation-station==100000) station=fstation

                hour=int((hour/10000))
                depar2big=0
                fg2big=0
                depar=0
                blacklist=0
                fg=0
                final=0
                reject=0
                event=0
		
!		WRITE(*,*)datum, hour, obstype, codetype, station, lat, lon, stalt, press, temp, biascorr, fg_depar, blacklist, varqc, depar, fg, final, varqcper, depar2big
		IF (varno.eq.2 .and. (hour .LE. 3 .OR. (hour .GE. 9 .AND. hour .LE. 15) .OR. hour .GE. 21)) THEN
			CALL calc_ind_flags(rcpara,datumneu, hour, press, time, hr, p, indexcor, blacklist, reject, depar, fg, final, depar2big, fg2big, flagges) 
!Put values into variables
!!$if (p .eq. 10) then
!!$print*,zeile
!!$print*,temp
!!$endif		
if (p.gt.0) then
			temperaturesall(time,p,hr)=temp
			fg_all(time,p,hr)=fg_depar
			an_all(time,p,hr)=an_depar
			sonde_type_all(time)=sonde_type
			bias_all(time,p,hr)=biascorr
			flags_all(time,p,hr)=flagges
			hoursall(time,hr)=hour
			datumall(time)=time
endif
!convert to small, one dimension array
!CALL dim3todim1(rcpara,tm, temperatures, vals, ind)
		END IF
                goto 501
500             write(*,'(a)') cstatnr,' read error',zeile
501	END DO ! end reading loop
	CLOSE(2)
        l=0
        do i=1,rcpara%nmax
          if(datumall(i) .eq. i) then
            l=l+1
          endif
        enddo
	vals=l
	ALLOCATE(arr(vals,rcpara%pmax,rcpara%parmax,5),datum(vals), hours(vals,rcpara%parmax),stypes(vals))
        l=0
        do i=1,rcpara%nmax
          if(datumall(i) .eq. i) then
            l=l+1
	    arr(l,:,:,1)=temperaturesall(i,:,:)
       	    arr(l,:,:,2)=fg_all(i,:,:)
	    arr(l,:,:,3)=bias_all(i,:,:)
	    arr(l,:,:,4)=flags_all(i,:,:)
       	    arr(l,:,:,5)=an_all(i,:,:)
	    datum(l)=datumall(i)
	    hours(l,:)=hoursall(i,:)
            stypes(l)=sonde_type_all(i)
          endif
        enddo
        print*,vals
!write netcdf 

	IF (vals .GT. 0)THEN
		CALL txttonc(rcpara,istat,outfile,arr,(/'temperatures','fg_dep','bias','flags','an_dep','s_type'/),vals, 5,datum, hours, lat, lon, alt,(/'temperature','fg-departures','bias correction','flags-final,fg,depar, blacklist, andepar2big, fg2big','an-departure','sonde type'/),stypes)
	END IF
	DEALLOCATE (arr, hours, datum)
!END DO !End loop over files
	CLOSE(1)
977	FORMAT(2x,I8,1x,I2,16x,I1,10x,I2,10x,A4,8x,I5,2(15x,F7.2),17x,I5,16x,I6,11x,I2,E23.15,3(1x,E21.17),29x,2I1,23x,3I4,22x,2I1)
978	FORMAT(2x,I8,1x,I2,16x,I1,10x,I2,10x,A4,8x,I5,2(15x,F7.2),17x,I5,16x,I6,11x,I2,7x,I16,3(1x,E21.17),29x,2I1,23x,3I4,22x,2I1)
979	FORMAT(2x,I8,1x,I2,16x,I1,10x,I2,10x,A4,8x,I5,2(15x,F7.2),17x,I5,16x,I6,11x,I2,4(1x,E21.17),13x,5I4)
980	FORMAT(2x,I8,1x,I2,16x,I1,10x,I2,10x,A4,8x,I5,2(15x,F7.2),17x,I5,16x,I6,11x,I2,88x,13x,5I4)
END SUBROUTINE txt_to_nc_wtype

! Read files with all vars (T,q,rh,u,v) 
SUBROUTINE txt_to_nc_presat(rcpara)
USE netcdf
USE var_def
USE rfmod
IMPLICIT NONE

type(rasocor_namelist) rcpara
REAL (kind=JPRM):: lat, lon 
INTEGER :: obstype, codetype, station
INTEGER :: itemp, ios, ios2, ss
INTEGER :: vals,l,i,filestat,istat
INTEGER, ALLOCATABLE :: datum(:), hours(:,:)
REAL (kind=JPRM) :: temp, biascorr, fg_depar,an_depar
REAL (kind=JPRM),ALLOCATABLE :: arr(:,:,:,:),stypes(:)
REAL (kind=JPRM), DIMENSION(rcpara%nmax,rcpara%pmax,rcpara%parmax) :: temperaturesall, fg_all, bias_all, flags_all,an_all
REAL (kind=JPRM), DIMENSION(rcpara%nmax) :: sonde_type_all
INTEGER :: datumall(rcpara%nmax)
INTEGER :: hoursall(rcpara%nmax,rcpara%parmax)
!filename variables
CHARACTER :: infile*60,inf*60,user*60,cstatnr*6
CHARACTER :: outfile*80

INTEGER :: alt, datumneu, hour, press  !Pressure in pa
INTEGER :: datumalt, time !zaehlt Datumsindex
INTEGER :: hr !hourindex, kann 1 oder 2 werden
INTEGER :: p !druckniveauindex
INTEGER :: indexcor ! wird -1 bei hour >=21, sonst 0
!INTEGER :: vals
!flag variables from inputfile
INTEGER :: final, fg, depar, reject, fg2big,sonde_type,varno,status,event
INTEGER ::  blacklist,depar2big ! analysis departure too big
!gesamtflag
INTEGER :: flagges
CHARACTER :: zeile*360,csonde_type*10,chilf*32,cstation*8

call getenv('HOME',user)
!OPEN (1,FILE=trim(user)//'/tables/filelist.txt')
!OPEN (1,FILE='filelist.txt')
!read filelist
!start loop over files

call getarg(2,infile)
istat=1
ios2=0
ss=0
!DO WHILE (ios2==0)
ss=ss+1
!	IF (ss .GT. 1438) THEN
	WRITE(*,*) ss
!	END IF
!	READ(1,'(A60)',IOSTAT=ios2) infile
!	IF (ios2 .NE. 0) EXIT
!create outfilename 
        if(len(trim(infile)) .eq. 10) then 
          inf='0'//infile
  	  outfile= inf(1:LEN(TRIM(inf))-3)//'nc'
        else
          inf=infile
	  outfile= infile(1:LEN(TRIM(infile))-3)//'nc'
        endif
        cstatnr=outfile(LEN(TRIM(inf))-11:LEN(TRIM(inf))-6)
        outfile='/vgc/srvx7/leo/scratch/presat/'//cstatnr//'/'//trim(outfile)
!        outfile='./'//trim(outfile)
        write(*,*) infile,cstatnr,outfile
        read(cstatnr,*) filestat
!Initialisieren
	CALL init(rcpara, datumalt, datumneu, datumall, hoursall, time, hr, p, vals, indexcor, temperaturesall, fg_all, an_all, bias_all, flags_all, lat, lon, alt, sonde_type_all)
!open file
	OPEN (2,FILE=TRIM(infile),iostat=ios)
	ios=0
!read values
        l=0
        datumall=rcpara%miss_val
	DO WHILE (ios==0)
		READ (2,'(A310)', IOSTAT=ios) zeile
		IF (ios .NE. 0) EXIT
                IF(zeile(22:25) .eq. 'NULL') zeile=zeile(1:21)//'0'//trim(zeile(23:300))
                ios=index(zeile,'NULL')
                do while(ios .ne. 0)
                  zeile(ios:ios+3)='-999'
                  ios=index(zeile,'NULL')
                enddo
              
				READ (zeile,*,iostat=ios) datumneu, hour, obstype, codetype, csonde_type, cstation, lat, lon, alt, press, varno, temp, biascorr, fg_depar, an_depar !, final, reject, event !, fg, final
                           
                if(ios .ne. 0) then
write(*,*) 'could not read',zeile
endif
                
                IF (cstation(2:2) == ':') THEN
                  READ(cstation(3:8),*) station
                ELSE
                  READ(cstation(1:6),*) station
                ENDIF
                IF (index(csonde_type,'NULL')>-1) THEN
                  sonde_type=-99999999
                ELSE
                  read(csonde_type,*) sonde_type
                ENDIF
                hour=int((hour/10000))
                depar2big=0
                fg2big=0
                depar=0
                blacklist=0
                fg=0
                final=0
                reject=0
                event=0
		
!		WRITE(*,*)datum, hour, obstype, codetype, station, lat, lon, stalt, press, temp, biascorr, fg_depar, blacklist, varqc, depar, fg, final, varqcper, depar2big
		IF ((filestat .eq. station .or. filestat .eq. (station+100000)) .and. varno.eq.2 .and. (hour .LE. 3 .OR. (hour .GE. 9 .AND. hour .LE. 15) .OR. hour .GE. 21)) THEN
			CALL calc_ind_flags(rcpara,datumneu, hour, press, time, hr, p, indexcor, blacklist, reject, depar, fg, final, depar2big, fg2big, flagges) 
!Put values into variables
		
			temperaturesall(time,p,hr)=temp
			fg_all(time,p,hr)=fg_depar
			an_all(time,p,hr)=an_depar
			sonde_type_all(time)=sonde_type
			bias_all(time,p,hr)=biascorr
			flags_all(time,p,hr)=flagges
			hoursall(time,hr)=hour
			datumall(time)=time
!convert to small, one dimension array
!CALL dim3todim1(rcpara,tm, temperatures, vals, ind)
		END IF
                goto 501
500             write(*,'(a)') cstatnr,' read error',zeile
501	END DO ! end reading loop
	CLOSE(2)
        l=0
        do i=1,rcpara%nmax
          if(datumall(i) .eq. i) then
            l=l+1
          endif
        enddo
	vals=l
	ALLOCATE(arr(vals,rcpara%pmax,rcpara%parmax,5),datum(vals), hours(vals,rcpara%parmax),stypes(vals))
        l=0
        do i=1,rcpara%nmax
          if(datumall(i) .eq. i) then
            l=l+1
	    arr(l,:,:,1)=temperaturesall(i,:,:)
       	    arr(l,:,:,2)=fg_all(i,:,:)
	    arr(l,:,:,3)=bias_all(i,:,:)
	    arr(l,:,:,4)=flags_all(i,:,:)
       	    arr(l,:,:,5)=an_all(i,:,:)
	    datum(l)=datumall(i)
	    hours(l,:)=hoursall(i,:)
            stypes(l)=sonde_type_all(i)
          endif
        enddo
        print*,vals
!write netcdf 

	IF (vals .GT. 0)THEN
		CALL txttonc(rcpara,istat,outfile,arr,(/'temperatures','fg_dep','bias','flags','an_dep','s_type'/),vals, 5,datum, hours, lat, lon, alt,(/'temperature','fg-departures','bias correction','flags-final,fg,depar, blacklist, andepar2big, fg2big','an-departure','sonde type'/),stypes)
	END IF
	DEALLOCATE (arr, hours, datum)
!END DO !End loop over files
	CLOSE(1)
977	FORMAT(2x,I8,1x,I2,16x,I1,10x,I2,10x,A4,8x,I5,2(15x,F7.2),17x,I5,16x,I6,11x,I2,E23.15,3(1x,E21.17),29x,2I1,23x,3I4,22x,2I1)
978	FORMAT(2x,I8,1x,I2,16x,I1,10x,I2,10x,A4,8x,I5,2(15x,F7.2),17x,I5,16x,I6,11x,I2,7x,I16,3(1x,E21.17),29x,2I1,23x,3I4,22x,2I1)
979	FORMAT(2x,I8,1x,I2,16x,I1,10x,I2,10x,A4,8x,I5,2(15x,F7.2),17x,I5,16x,I6,11x,I2,4(1x,E21.17),13x,5I4)
980	FORMAT(2x,I8,1x,I2,16x,I1,10x,I2,10x,A4,8x,I5,2(15x,F7.2),17x,I5,16x,I6,11x,I2,88x,13x,5I4)
END SUBROUTINE txt_to_nc_presat

! Read files with all vars (T,q,rh,u,v) 
SUBROUTINE txt_to_nc_oper(rcpara)
USE netcdf
USE var_def
USE rfmod
IMPLICIT NONE

type(rasocor_namelist) rcpara
REAL (kind=JPRM):: lat, lon 
INTEGER :: obstype, codetype, station
INTEGER :: itemp, ios, ios2, ss
INTEGER :: vals,l,i,filestat,istat
INTEGER, ALLOCATABLE :: datum(:), hours(:,:)
REAL (kind=JPRM) :: temp, biascorr, fg_depar,an_depar
REAL (kind=JPRM),ALLOCATABLE :: arr(:,:,:,:),stypes(:)
REAL (kind=JPRM), DIMENSION(rcpara%nmax,rcpara%pmax,rcpara%parmax) :: temperaturesall, fg_all, bias_all, flags_all,an_all
REAL (kind=JPRM), DIMENSION(rcpara%nmax) :: sonde_type_all
INTEGER :: datumall(rcpara%nmax)
INTEGER :: hoursall(rcpara%nmax,rcpara%parmax)
!filename variables
CHARACTER :: infile*60,inf*60,user*60,cstatnr*5
CHARACTER :: outfile*90

INTEGER :: alt, datumneu, hour, press  !Pressure in pa
INTEGER :: datumalt, time !zaehlt Datumsindex
INTEGER :: hr !hourindex, kann 1 oder 2 werden
INTEGER :: p !druckniveauindex
INTEGER :: indexcor ! wird -1 bei hour >=21, sonst 0
!INTEGER :: vals
!flag variables from inputfile
INTEGER :: final, fg, depar, reject, fg2big,sonde_type,varno,status,event
INTEGER ::  blacklist,depar2big ! analysis departure too big
!gesamtflag
INTEGER :: flagges 
CHARACTER :: zeile*360,csonde_type*4,chilf*32

call getenv('HOME',user)
!OPEN (1,FILE=trim(user)//'/tables/filelist.txt')
!OPEN (1,FILE='filelist.txt')
!read filelist
!start loop over files

call getarg(2,infile)
istat=1
ios2=0
ss=0
!DO WHILE (ios2==0)
ss=ss+1
!	IF (ss .GT. 1438) THEN
	WRITE(*,*) ss
!	END IF
!	READ(1,'(A60)',IOSTAT=ios2) infile
!	IF (ios2 .NE. 0) EXIT
!create outfilename 
        if(len(trim(infile)) .eq. 10) then 
          inf='0'//infile
  	  outfile= inf(1:LEN(TRIM(inf))-3)//'nc'
        else
          inf=infile
	  outfile= infile(1:LEN(TRIM(infile))-3)//'nc'
        endif
        cstatnr=outfile(LEN(TRIM(inf))-10:LEN(TRIM(inf))-6)
        outfile='/vgc/srvx7/leo/scratch/oper/'//cstatnr//'/'//trim(outfile)
!        outfile='./'//trim(outfile)
        write(*,*) infile,cstatnr,outfile
        read(cstatnr,*) filestat
!Initialisieren
	CALL init(rcpara, datumalt, datumneu, datumall, hoursall, time, hr, p, vals, indexcor, temperaturesall, fg_all, an_all, bias_all, flags_all, lat, lon, alt, sonde_type_all)
!open file
	OPEN (2,FILE=TRIM(infile),iostat=ios)
	ios=0
!read values
        l=0
        datumall=rcpara%miss_val
	DO WHILE (ios==0)
		READ (2,'(A310)', IOSTAT=ios) zeile
		IF (ios .NE. 0) EXIT
                IF(zeile(22:25) .eq. 'NULL') zeile=zeile(1:21)//'0'//trim(zeile(23:300))
                ios=index(zeile,'NULL')
                do while(ios .ne. 0)
                  zeile(ios:ios+3)='-999'
                  ios=index(zeile,'NULL')
                enddo
              
				READ (zeile,*,iostat=ios) datumneu, hour, obstype, codetype, sonde_type, station, lat, lon, press, varno, temp, biascorr, fg_depar, an_depar !, fg, final
                           
                if(ios .ne. 0) then
write(*,*) 'could not read',zeile
endif
                hour=int((hour/10000))
                depar2big=0
                fg2big=0
                depar=0
                blacklist=0
                fg=0
                final=0
                reject=0
                event=0
		
!		WRITE(*,*)datum, hour, obstype, codetype, station, lat, lon, stalt, press, temp, biascorr, fg_depar, blacklist, varqc, depar, fg, final, varqcper, depar2big
		IF (filestat .eq. station .and. varno.eq.2 .and. (hour .LE. 3 .OR. (hour .GE. 9 .AND. hour .LE. 15) .OR. hour .GE. 21)) THEN
			CALL calc_ind_flags(rcpara,datumneu, hour, press, time, hr, p, indexcor, blacklist, reject, depar, fg, final, depar2big, fg2big, flagges) 
!Put values into variables
		
			temperaturesall(time,p,hr)=temp
			fg_all(time,p,hr)=fg_depar
			an_all(time,p,hr)=an_depar
			sonde_type_all(time)=sonde_type
			bias_all(time,p,hr)=biascorr
			flags_all(time,p,hr)=flagges
			hoursall(time,hr)=hour
			datumall(time)=time
!convert to small, one dimension array
!CALL dim3todim1(rcpara,tm, temperatures, vals, ind)
		END IF
                goto 501
500             write(*,'(a)') cstatnr,' read error',zeile
501	END DO ! end reading loop
	CLOSE(2)
        l=0
        do i=1,rcpara%nmax
          if(datumall(i) .eq. i) then
            l=l+1
          endif
        enddo
	vals=l
	ALLOCATE(arr(vals,rcpara%pmax,rcpara%parmax,5),datum(vals), hours(vals,rcpara%parmax),stypes(vals))
        l=0
        do i=1,rcpara%nmax
          if(datumall(i) .eq. i) then
            l=l+1
	    arr(l,:,:,1)=temperaturesall(i,:,:)
       	    arr(l,:,:,2)=fg_all(i,:,:)
	    arr(l,:,:,3)=bias_all(i,:,:)
	    arr(l,:,:,4)=flags_all(i,:,:)
       	    arr(l,:,:,5)=an_all(i,:,:)
	    datum(l)=datumall(i)
	    hours(l,:)=hoursall(i,:)
            stypes(l)=sonde_type_all(i)
          endif
        enddo
        print*,vals
!write netcdf
	IF (vals .GT. 0)THEN
		CALL txttonc(rcpara,istat,outfile,arr,(/'temperatures','fg_dep','bias','flags','an_dep','s_type'/),vals, 5,datum, hours, lat, lon, alt,(/'temperature','fg-departures','bias correction','flags-final,fg,depar, blacklist, andepar2big, fg2big','an-departure','sonde type'/),stypes)
	END IF
        outfile='/vgc/srvx7/leo/scratch/oper/'//cstatnr//'/opermon'//cstatnr//'.nc'
write(*,*) outfile
        call write_sonde_monthly_nc(outfile,rcpara,temperaturesall,istat,ios,lon,lat,bias_all,fg_all) ! schreibt Files mit Monatsmittel
	DEALLOCATE (arr, hours, datum)
!END DO !End loop over files
	CLOSE(1)
977	FORMAT(2x,I8,1x,I2,16x,I1,10x,I2,10x,A4,8x,I5,2(15x,F7.2),17x,I5,16x,I6,11x,I2,E23.15,3(1x,E21.17),29x,2I1,23x,3I4,22x,2I1)
978	FORMAT(2x,I8,1x,I2,16x,I1,10x,I2,10x,A4,8x,I5,2(15x,F7.2),17x,I5,16x,I6,11x,I2,7x,I16,3(1x,E21.17),29x,2I1,23x,3I4,22x,2I1)
979	FORMAT(2x,I8,1x,I2,16x,I1,10x,I2,10x,A4,8x,I5,2(15x,F7.2),17x,I5,16x,I6,11x,I2,4(1x,E21.17),13x,5I4)
980	FORMAT(2x,I8,1x,I2,16x,I1,10x,I2,10x,A4,8x,I5,2(15x,F7.2),17x,I5,16x,I6,11x,I2,88x,13x,5I4)
END SUBROUTINE txt_to_nc_oper

SUBROUTINE init(rcpara, datumalt, datumneu, datumall, hoursall, time, hr, p, vals, indexcor, temperaturesall, fg_all, an_all, bias_all, flags_all, lat, lon, alt , sonde_type_all)
USE rfmod
IMPLICIT NONE
type(rasocor_namelist) rcpara
REAL (kind=JPRM):: lat, lon 
REAL (kind=JPRM), DIMENSION(rcpara%nmax,rcpara%pmax,rcpara%parmax) :: temperaturesall, fg_all, bias_all, an_all,flags_all
REAL (kind=JPRM),OPTIONAL, DIMENSION(rcpara%nmax) ::sonde_type_all
INTEGER :: datumall(rcpara%nmax)
INTEGER :: hoursall(rcpara%nmax,rcpara%parmax)
INTEGER :: alt, datumneu
INTEGER :: datumalt, time !zaehlt Datumsindex
INTEGER :: hr !hourindex, kann 1 oder 2 werden
INTEGER :: p !druckniveauindex
INTEGER :: indexcor ! wird -1 bei hour >=21, sonst 0
INTEGER :: vals
INTEGER :: flagges 

	datumalt=0
	datumneu=0
	datumall=rcpara%miss_val
	time=0
	hr=0
	p=0
	vals=0
	indexcor=0
	alt=-999
	hoursall=-999
	flagges=-999
	lat=rcpara%miss_val
	lon=rcpara%miss_val
	temperaturesall=rcpara%miss_val
	fg_all=rcpara%miss_val
	an_all=rcpara%miss_val
	bias_all=rcpara%miss_val
	flags_all=rcpara%miss_val
	if(present(sonde_type_all)) sonde_type_all=rcpara%miss_val
END SUBROUTINE

SUBROUTINE calc_ind_flags(rcpara,datumneu, hour, press, time, hr, p, indexcor, blacklist, reject, depar, fg, final, depar2big, fg2big, flagges) 

IMPLICIT NONE
type(rasocor_namelist) rcpara
INTEGER :: datumneu, hour, press  !Pressure in pa
INTEGER :: datumalt, time !zaehlt Datumsindex
INTEGER :: hr !hourindex, kann 1 oder 2 werden
INTEGER :: p !druckniveauindex
INTEGER :: indexcor ! wird -1 bei hour >=21, sonst 0
!flag variables from inputfile
INTEGER :: final, fg, depar, reject, fg2big
INTEGER ::  blacklist,depar2big ! analysis departure too big
!gesamtflag
INTEGER :: flagges 
! calculate indices

        time=toindex(datumneu,rcpara)
        if(time .lt. 1 .or. time .gt. rcpara%nmax) then
          print*,"spurious date",datumneu
        endif

	IF (hour .LE. 3) THEN
		hr=1
		indexcor=0
	ELSEIF (hour .GE. 9 .AND. hour .LE. 15) THEN
		hr=2
	ELSEIF (hour .GE. 21) THEN
		hr=1
		time=time+1
	END IF
        p=0
	IF (press==1000) p=1
	IF (press==2000) p=2
	IF (press==3000) p=3
	IF (press==5000) p=4
	IF (press==7000) p=5
	IF (press==10000) p=6
	IF (press==15000) p=7
	IF (press==20000) p=8
	IF (press==25000) p=9
	IF (press==30000) p=10
	IF (press==40000) p=11
	IF (press==50000) p=12
	IF (press==70000) p=13
	IF (press==85000) p=14
	IF (press==92500) p=15
	IF (press==100000) p=16
        if(p .lt. 0 .or. p .gt. rcpara%pmax) print*,"spurious pindex",p,press
        if(hr .lt. 1 .or. hr .gt. 2) print*,"spurious parindex",hr
!calculate flags
	IF (depar==10) depar=2
	IF (depar==11) depar=3
	IF (depar>11 .OR. depar<0) WRITE(*,*) 'unknown depar flag', depar
	IF (fg==10) fg=2
	IF (fg==11) fg=3
	IF (fg>11 .OR. fg<0) WRITE(*,*) 'unknown fg flag', fg
	IF (final==10) final=2
	IF (final==11) final=3
	IF (final>11 .OR. final<0) WRITE(*,*) 'unknown final flag', final
	depar2big=depar2big*10
	reject=reject*100
	blacklist=blacklist*1000
	depar=depar*10000
	fg=fg*100000
	final=final*1000000
	flagges=final+fg+depar+blacklist+reject+depar2big+fg2big

END SUBROUTINE

SUBROUTINE txttonc(rcpara,mstat,outfile,arr,names,vals,numpar,datum, hours, lat, lon, alt,longnames,stype,stname)
USE netcdf
USE var_def
USE rfmod
IMPLICIT NONE

type(rasocor_namelist) rcpara
type(atts) attribs
INTEGER :: vals, numdat,istat
INTEGER,intent(in) :: mstat
REAL(kind=JPRM) :: temperatures(vals,rcpara%pmax,rcpara%parmax), fg_dep(vals,rcpara%pmax,rcpara%parmax), biasc(vals,rcpara%pmax,rcpara%parmax), flags(vals,rcpara%pmax,rcpara%parmax)
REAL(kind=JPRM),optional :: stype(vals,1)
INTEGER :: datum(vals,1), hours(vals,rcpara%parmax)
INTEGER :: alt
REAL(kind=JPRM) :: lat, lon
INTEGER :: numpar, maxatts,pmax
CHARACTER*80 :: outfile
CHARACTER*(*),optional:: longnames(numpar)
REAL(kind=JPRM) :: arr(vals,rcpara%pmax,rcpara%parmax,numpar)
CHARACTER*12 :: names(numpar)
CHARACTER*50,OPTIONAL:: stname
numdat=1
maxatts=4
if(lserial) then
  istat=rcpara%statmax+1
else
  istat=mstat
endif
ALLOCATE (attribs%attname(maxatts), attribs%long_name(15), attribs%units(15), attribs%missing_val(15), attribs%valid_range(15,2))
!ALLOCATE (arr(vals,rcpara%pmax,rcpara%parmax,numpar), names(numpar))
!names = (/'temperatures','fg_dep','bias','flags'/)
!test differenz 12-00

attribs%attname=(/'long_name','units','missing_value','valid_range'/)
if(present(longnames)) then
attribs%long_name=longnames
else
attribs%long_name=names
endif
attribs%units=(/'K','K','K','','K','K','K','K','K','K','K','K','K','K','K'/)
attribs%missing_val=rcpara%miss_val
attribs%valid_range(:,1)=(/0.,-30.,-20.,0.,-30.,-30.,-30.,-30.,-30.,-30.,-30.,-30.,-30.,-30.,-30./)
attribs%valid_range(:,2)=(/400.,30.,20.,300000.,30.,30.,30.,30.,30.,30.,30.,30.,30.,30.,30./)
attribs%title='station daily temperature series'
! create netcdf file
!!$omp critical
if(vals .gt. 0) then
	CALL gen_nc(rcpara,istat,outfile,numpar,maxatts,vals, arr, names, attribs, lat, lon,datum, numdat, alt, hours,stype=stype,stname=stname)  
endif
!!$omp end critical

deallocate(attribs%attname,attribs%long_name,attribs%units,attribs%missing_val,attribs%valid_range)
END SUBROUTINE txttonc

SUBROUTINE gen_nc(rcpara,mstat,outfile,numpar,maxatts,vals, arr, names, attribs, lat, lon, datum, numdat, alt, hours,stype,stname) !
  !creates netcdffile, with dimensions: station=1, time=vals, pressure=rcpara%pmax, hour=2
  USE var_def
  USE rfmod
  USE typeSizes
  USE netcdf
  IMPLICIT NONE

  type(rasocor_namelist) rcpara
  type(atts) attribs
  INTEGER :: ncid, statid, timid, pressid, hourid, latvarid, lonvarid, altvarid, hoursvarid, pressvarid, datumvarid, numdatid, climboundvarid,stypevarid
  INTEGER :: vals !number of available days for a station
  INTEGER :: numpar !number of arrays with (vals,pmax, parmax) dimension to write to netcdf
  INTEGER :: maxatts !maximum number of attributes per variable
  INTEGER :: arrvarid(numpar)
  INTEGER :: status, ii, jj, numdat,istat,form(4)
  INTEGER,intent(in) :: mstat
  INTEGER :: datum(vals, numdat) 
  INTEGER, ALLOCATABLE :: dhilf(:),xhours(:,:)
  INTEGER, OPTIONAL :: alt
  INTEGER, OPTIONAL :: hours(vals,rcpara%parmax)
  REAL(kind=JPRM),optional :: stype(vals)
  REAL(kind=JPRM) :: arr(vals,rcpara%pmax,rcpara%parmax,numpar) ! vals,pmax,parmax,numpar
  REAL(kind=JPRM) :: lat, lon
  CHARACTER*(*) :: outfile
  CHARACTER*(*) :: names(numpar) !variable names
  CHARACTER*50,optional :: stname
  CHARACTER*8 :: today, date !Funktion, die das aktuelle Datum  ausgibt
  CHARACTER*8 :: today2
  !REAL, PARAMETER :: filval=-999.0
  INTEGER :: pmax !number of vertical levels
  character*4::cstartyear

  if(lserial) then
     istat=rcpara%statmax+1
  else
     istat=mstat+4000
  endif
  !istat=rcpara%statmax+1
  !istat=mstat+4000

  !$ call omp_set_lock(omp_lp(istat))
  !print*,istat,trim(outfile)
  status=NF90_CREATE(TRIM(outfile),NF90_NETCDF4,ncid)
  !write(*,*) 'modulo',mod(ncid,65536),ncid/65536

  IF (status /= NF90_NOERR) THEN
     PRINT *,'Error generating netcdf ', outfile
     status=NF90_CLOSE(ncid)
     !$ call omp_unset_lock(omp_lp(istat))
     RETURN
  END IF
!!$ call omp_unset_lock(omp_lp(istat))

  ! define dimensions
  status=NF90_DEF_DIM(ncid,'station',1,statid) !nur eine Station pro File
  IF (status /= NF90_NOERR) PRINT *,'Error generating station dimension! in ', outfile
  status=NF90_DEF_DIM(ncid,'numdat',numdat,numdatid)
  IF (status /= NF90_NOERR) PRINT *,'Error generating numdat dimension!'
  status=NF90_DEF_DIM(ncid,'time',vals,timid) !vals ist die Azahl der existierenden Werte fuer die Station
  IF (status /= NF90_NOERR) PRINT *,'Error generating time dimension!'
  status=NF90_DEF_DIM(ncid,'pressure',rcpara%pmax,pressid) !pmax =maximale Anzahl der Druckniveaus, ist normalerweise 16
  IF (status /= NF90_NOERR) PRINT *,'Error generating pressure dimension!'
  status=NF90_DEF_DIM(ncid,'hour',rcpara%parmax,hourid)
  IF (status /= NF90_NOERR) PRINT *,'Error generating hour dimension!'

  ! define variables
  status=NF90_DEF_VAR(ncid,'lat',NF90_FLOAT,statid,latvarid) !latitude (real)
  IF (status /= NF90_NOERR) PRINT *,'Error generating latitude variable!'
  status=NF90_DEF_VAR(ncid,'lon',NF90_FLOAT,statid,lonvarid)
  ! longitude (real)
  IF (status /= NF90_NOERR) PRINT *,'Error generating longitude variable!'
  IF (present(alt)) THEN
     status=NF90_DEF_VAR(ncid,'alt',NF90_FLOAT,statid,altvarid) !altitude (real)
     IF (status /= NF90_NOERR) PRINT *,'Error generating altitude variable!'
  END IF
  status=NF90_DEF_VAR(ncid,'press',NF90_FLOAT,pressid, pressvarid) !pressure
  IF (status /= NF90_NOERR) PRINT *,'Error generating pressure variable!'
  status=NF90_DEF_VAR(ncid,'datum',NF90_INT,(/timid,numdatid/), datumvarid) !datum
  IF (present(hours)) THEN
     status=NF90_DEF_VAR(ncid,'hours',NF90_INT,(/timid,hourid/),hoursvarid)! launch time
     IF (status /= NF90_NOERR) PRINT *,'Error generating hours variable!'
  END IF
  DO ii=1,numpar
     status=NF90_DEF_VAR(ncid,TRIM(names(ii)),NF90_FLOAT,(/timid,pressid,hourid/),arrvarid(ii)) 
     !      		status=NF90_DEF_VAR(ncid,TRIM(names(ii)),NF90_FLOAT,(/timid,pressid,hourid/),arrvarid(ii),deflate_level=5) 
     IF (status /= NF90_NOERR) PRINT *,'Error generating ', TRIM(names(ii)), ' variable!'
  END DO

  IF (present(stype)) status=NF90_DEF_VAR(ncid,'s_type',NF90_FLOAT,(/timid/), stypevarid) !datum
  IF (status /= NF90_NOERR) PRINT *,'Error generating  variable!'

  ! define attributes
  status=NF90_PUT_ATT(ncid,latvarid,'long_name','station latitude')  
  IF (status /= NF90_NOERR) PRINT *,'Error giveing long_name attribut to station latitude'
  status=NF90_PUT_ATT(ncid,latvarid,'units','degrees_north') 
  IF (status /= NF90_NOERR) PRINT *,'Error giveing units attribut to station latitude' 
  status=NF90_PUT_ATT(ncid,latvarid,'axis','Y')      
  status=NF90_PUT_ATT(ncid,latvarid,'valid_range',(/-90.0, 90.0/))      
  status=NF90_PUT_ATT(ncid,latvarid,'missing_value',rcpara%miss_val)    
  !      PRINT *, NF90_FILL_REAL  

  status=NF90_PUT_ATT(ncid,lonvarid,'long_name','station longitude')      
  status=NF90_PUT_ATT(ncid,lonvarid,'units','degrees_east')      
  status=NF90_PUT_ATT(ncid,lonvarid,'axis','X')      
  status=NF90_PUT_ATT(ncid,lonvarid,'valid_range',(/-180.0, 180.0/))      
  status=NF90_PUT_ATT(ncid,lonvarid,'missing_value',rcpara%miss_val) 
  IF (present(alt)) THEN     
     status=NF90_PUT_ATT(ncid,altvarid,'long_name','station altitude')
     status=NF90_PUT_ATT(ncid,altvarid,'units','m')
     status=NF90_PUT_ATT(ncid,altvarid,'axis','Z')
     status=NF90_PUT_ATT(ncid,altvarid,'valid_range',(/-100, 8000/))
     status=NF90_PUT_ATT(ncid,altvarid,'comment','station altitude in meter msl')
     status=NF90_PUT_ATT(ncid,altvarid,'missing_value',rcpara%miss_val) 
  END IF
  status=NF90_PUT_ATT(ncid,datumvarid,'long_name','datum') 
  write(cstartyear,'(I4)') rcpara%startdate/10000     
  status=NF90_PUT_ATT(ncid,datumvarid,'units','days since '//cstartyear//'-01-01 0:0:0')      
  status=NF90_PUT_ATT(ncid,datumvarid,'axis','T')      
  status=NF90_PUT_ATT(ncid,datumvarid,'calendar','gregorian') 
  status=NF90_PUT_ATT(ncid,datumvarid,'missing_value',rcpara%miss_val)
  if (rcpara%plevs(1)<10.) then   ! assume MSU layers
     status=NF90_PUT_ATT(ncid,pressvarid,'long_name','MSU layers')      
     status=NF90_PUT_ATT(ncid,pressvarid,'units','')      
     status=NF90_PUT_ATT(ncid,pressvarid,'axis','Z')      
     status=NF90_PUT_ATT(ncid,pressvarid,'valid_range', (/1, 4/))      
  else
     status=NF90_PUT_ATT(ncid,pressvarid,'long_name','pressure levels')      
     status=NF90_PUT_ATT(ncid,pressvarid,'units','hPa')      
     status=NF90_PUT_ATT(ncid,pressvarid,'axis','Z')      
     status=NF90_PUT_ATT(ncid,pressvarid,'valid_range', (/0.00, 1100.0/))      
  endif
  status=NF90_PUT_ATT(ncid,pressvarid,'missing_value',rcpara%miss_val)
  IF (present(hours)) THEN
     status=NF90_PUT_ATT(ncid,hoursvarid,'long_name','launch time')
     status=NF90_PUT_ATT(ncid,hoursvarid,'units','hr')
     status=NF90_PUT_ATT(ncid,hoursvarid,'valid_range',(/00, 23/))
     status=NF90_PUT_ATT(ncid,hoursvarid,'missing_value',floor(rcpara%miss_val))
  END IF
  IF (present(stype)) THEN
     status=NF90_PUT_ATT(ncid,stypevarid,'long_name','Radiosonde Type')
     status=NF90_PUT_ATT(ncid,stypevarid,'units','')
     status=NF90_PUT_ATT(ncid,stypevarid,'valid_range',(/00., 1000./))
     status=NF90_PUT_ATT(ncid,stypevarid,'missing_value',rcpara%miss_val)
  END IF
  DO ii=1,numpar
     IF (TRIM(names(ii)) .EQ. 'flags') THEN
        status=NF90_PUT_ATT(ncid,arrvarid(ii),'flag_values','0000000 3331111')
        status=NF90_PUT_ATT(ncid,arrvarid(ii),'flag_meanings','first_digit_final_flag second_digit_first_guess_flag third_digit_departure_flag fourth_digit_blacklist_flag fifth_digit_reject_flag sixth_digit_departure_too_big seventh_digit_first_guess_too_big' )
     END IF
     DO jj=1,maxatts
        IF (attribs%attname(jj)(1:9) .EQ. 'long_name' .AND. attribs%long_name(ii) .NE. '') THEN
           status=NF90_PUT_ATT(ncid,arrvarid(ii),TRIM(attribs%attname(jj)),TRIM(attribs%long_name(ii)))
        ELSEIF (attribs%attname(jj)(1:5) .EQ. 'units' .AND. attribs%units(ii) .NE. '') THEN
           status=NF90_PUT_ATT(ncid,arrvarid(ii),TRIM(attribs%attname(jj)),TRIM(attribs%units(ii)))
        ELSEIF (attribs%attname(jj)(1:13) .EQ. 'missing_value') THEN
           status=NF90_PUT_ATT(ncid,arrvarid(ii),TRIM(attribs%attname(jj)),attribs%missing_val(ii))
        ELSEIF (attribs%attname(jj)(1:11) .EQ. 'valid_range' .AND. attribs%valid_range(ii,1) .NE. rcpara%miss_val .AND. attribs%valid_range(ii,2) .NE. rcpara%miss_val) THEN
           status=NF90_PUT_ATT(ncid,arrvarid(ii),TRIM(attribs%attname(jj)),(/attribs%valid_range(ii,1),attribs%valid_range(ii,2)/))
        ELSEIF (attribs%attname(jj)(1:12) .EQ. 'cell_methods') THEN
           IF (attribs%cell_methods(ii) .NE. '') THEN
              status=NF90_PUT_ATT(ncid,arrvarid(ii),TRIM(attribs%attname(jj)),TRIM(attribs%cell_methods(ii)))
           END IF
        ELSEIF (TRIM(attribs%attname(jj)) .NE. 'long_name' .AND. TRIM(attribs%attname(jj)) .NE. 'units' .AND. TRIM(attribs%attname(jj)) .NE. 'missing_value' .AND. TRIM(attribs%attname(jj)) .NE. 'valid_range' .AND. TRIM(attribs%attname(jj)) .NE. 'cell_methods') THEN
           WRITE(*,*) 'error writing attributes to nc-file, unknown attribute, (in Subroutine gen_nc in File read_txt_write_nc.f90)'

        END IF
     END DO
  END DO
  status=NF90_PUT_ATT(ncid,NF90_GLOBAL,'Conventions','CF-1.1')
  status=NF90_PUT_ATT(ncid,NF90_GLOBAL,'title',TRIM(attribs%title))
  status=NF90_PUT_ATT(ncid,NF90_GLOBAL,'institution','University of Vienna')
  today=date() !generating date string
  today2=today(7:8)//'/'//today(1:2)//'/'//today(4:5)
  status=NF90_PUT_ATT(ncid,NF90_GLOBAL,'history',trim(today2))
  status=NF90_PUT_ATT(ncid,NF90_GLOBAL,'source','radiosonde, ERA-Interim, ERA-40, RAOBCORE')
  status=NF90_PUT_ATT(ncid,NF90_GLOBAL,'references','www.univie.ac.at/theoret-met/research/raobcore')
  if(present(stname)) then
     status=NF90_PUT_ATT(ncid,NF90_GLOBAL,'Stationname',trim(stname))
  endif
  ! end definitions
  status=NF90_ENDDEF(ncid)
  IF (status /= NF90_NOERR) PRINT *,'Error ending definition phase!',outfile, NF90_STRERROR(status)

  ! assign variables
  status=NF90_PUT_VAR(ncid,latvarid,lat)
  IF (status /= NF90_NOERR) PRINT *,'could not write lat', NF90_STRERROR(status)
  status=NF90_PUT_VAR(ncid,lonvarid,lon)
  IF (present(alt)) THEN
     status=NF90_PUT_VAR(ncid,altvarid,alt)
     IF (status /= NF90_NOERR) PRINT *,'could not write alt', NF90_STRERROR(status)
  END IF
  status=NF90_PUT_VAR(ncid,pressvarid,rcpara%plevs(1:rcpara%pmax))
  IF (status /= NF90_NOERR) THEN
     PRINT *,'could not write press', NF90_STRERROR(status)
  ENDIF
  !print*,'datum',vals,numdat,datumvarid,shape(datum),datum(vals,:)

  !      call dput(ncid,datumvarid,datum(1:vals,:),status,vals,numpar)
  allocate(dhilf(vals))
  dhilf=datum(1:vals,1)
  status=NF90_PUT_VAR(ncid,datumvarid,dhilf)
  deallocate(dhilf)
  IF (status /= NF90_NOERR) PRINT *,'could not write datum', NF90_STRERROR(status)
  !print*,datum(vals,:)
  !write(*,*) hours(vals,:),hoursvarid,timid,hourid,size(hours)

  IF (present(hours)) THEN
     status=NF90_PUT_VAR(ncid,hoursvarid,hours(1:vals,:))
     print*,nf90_strerror(status)
     !                allocate(xhours(vals,2))
     !                xhours=hours(1:vals,:)
     !      		status=NF90_PUT_VAR(ncid,hoursvarid,xhours)
     !                deallocate(xhours)
     !                print*,nf90_strerror(status)
  END IF
  DO ii=1,numpar
     status=NF90_PUT_VAR(ncid,arrvarid(ii),arr(:,:,:,ii)) 
     IF (status /= NF90_NOERR) PRINT *,'Error writing ', TRIM(names(ii)), ' variable!'
  END DO

  if(present(stype)) then 
     !        print*,'stype',stype(vals),stypevarid
     status=NF90_PUT_VAR(ncid,stypevarid,stype(1:vals))
     IF (status /= NF90_NOERR) PRINT *,'could not write stype', NF90_STRERROR(status)
  endif
  status=NF90_CLOSE(ncid)
  IF (status /= NF90_NOERR) PRINT *,'Error closing netcdf!'
  !$ call omp_unset_lock(omp_lp(istat))

  RETURN
END SUBROUTINE gen_nc

SUBROUTINE gen_nc_3D(outfile,dimnames,varnames, attribs, lon, lat, datum, arr, plevs,hours) !
!creates netcdffile, with dimensions: 
USE var_def
USE rfmod
USE netcdf
IMPLICIT NONE

type(rasocor_namelist) rcpara
type(atts) attribs
INTEGER :: ncid
INTEGER :: numpar !number of arrays with (vals,pmax, parmax) dimension to write to netcdf
INTEGER :: maxatts !maximum number of attributes per variable
INTEGER :: status, ii, jj, numdat,istat,form(4),i,varc
INTEGER :: datum(:) 
INTEGER, OPTIONAL :: hours(:)
REAL(kind=JPRM),OPTIONAL :: plevs(:)
INTEGER,allocatable :: dimid(:),varid(:)
INTEGER :: s(5),sd,vd
REAL(kind=JPRM) :: arr(:,:,:,:,:)
REAL(kind=JPRM) :: lat(:), lon(:)
CHARACTER*(*) :: outfile
CHARACTER*(*) :: dimnames(:),varnames(:) !variable names
CHARACTER*8 :: today, date !Funktion, die das aktuelle Datum ausgibt
CHARACTER*8 :: today2

!REAL, PARAMETER :: filval=-999.0
INTEGER :: pmax !number of vertical levels

status=NF90_CREATE(TRIM(outfile),NF90_CLOBBER,ncid)

      	IF (status /= NF90_NOERR) THEN
      		PRINT *,'Error generating netcdf ', outfile
		RETURN
	END IF	

! define dimensions
      s=shape(arr)
      sd=size(dimnames,1)
      vd=size(varnames,1)
      allocate(dimid(sd),varid(vd))
      
      do i=1,sd
        status=NF90_DEF_DIM(ncid,trim(dimnames(i)),s(i),dimid(i)) !nur eine Station pro File
        IF (status /= NF90_NOERR) PRINT *,'Error generating station dimension! in ', outfile
      enddo
     
      ! define variables
      do i=1,vd-1
        status=NF90_DEF_VAR(ncid,varnames(i),NF90_FLOAT,dimid(i),varid(i)) !latitude (real)
        IF (status /= NF90_NOERR) PRINT *,'Error generating '//trim(varnames(i))//' variable!'
      enddo
        status=NF90_DEF_VAR(ncid,trim(varnames(i)),NF90_FLOAT,dimid(1:sd),varid(i)) !latitude (real)
        IF (status /= NF90_NOERR) PRINT *,'Error generating '//trim(varnames(i))//' variable!'
      
      
! define attributes
!      do i=1,vd
!        status=NF90_PUT_ATT(ncid,varid(i),attnames(1),varnames(i))  
!        status=NF90_PUT_ATT(ncid,varid(i),attnames(2),varnames(i))  
!      enddo
      status=NF90_PUT_ATT(ncid,NF90_GLOBAL,'Conventions','CF-1.1')
      status=NF90_PUT_ATT(ncid,NF90_GLOBAL,'title',TRIM(attribs%title))
      status=NF90_PUT_ATT(ncid,NF90_GLOBAL,'institution','University of Vienna')
      today=date() !generating date string
      today2=today(7:8)//'/'//today(1:2)//'/'//today(4:5)
      status=NF90_PUT_ATT(ncid,NF90_GLOBAL,'history',today2)
      status=NF90_PUT_ATT(ncid,NF90_GLOBAL,'source','radiosonde, ERA-Interim, ERA-40, RAOBCORE')
      status=NF90_PUT_ATT(ncid,NF90_GLOBAL,'references','www.univie.ac.at/theoret-met/research/raobcore')
! end definitions
      status=NF90_ENDDEF(ncid)
      IF (status /= NF90_NOERR) PRINT *,'Error ending definition phase!', NF90_STRERROR(status)

! assign variables
      IF (status /= NF90_NOERR) PRINT *,'could not write lat', NF90_STRERROR(status)
        status=NF90_PUT_VAR(ncid,varid(1),lon)
        status=NF90_PUT_VAR(ncid,varid(2),lat)
        varc=3
        if(present(plevs)) then
          status=NF90_PUT_VAR(ncid,varid(varc),plevs)
          varc=varc+1
        endif
        if(present(hours)) then
          status=NF90_PUT_VAR(ncid,varid(varc),hours)
          varc=varc+1
        endif
        status=NF90_PUT_VAR(ncid,varid(varc),datum)
        varc=varc+1
     	status=NF90_PUT_VAR(ncid,varid(varc),arr) 
      	IF (status /= NF90_NOERR) PRINT *,'Error writing ', TRIM(varnames(4)), ' variable!'
      
      status=NF90_CLOSE(ncid)
      IF (status /= NF90_NOERR) PRINT *,'Error closing netcdf!'

     RETURN
END SUBROUTINE gen_nc_3D

!vals, ind und temperatures sind dann im netcdf zu speichern

subroutine write_sonde_corr_daily_nc(filename,rcpara,mstat,err,rasocorr,lon,lat,rasobreak,breakuncertainty,stname) !
!schreibt Files mit Korrekturen
USE var_def
USE rfmod
implicit none

type(rasocor_namelist),intent(in) :: rcpara
type(atts) attribs
integer ini, ipmax, iparmax, i, err, imax, l, bi, di, gi, vals, ip, ipar, jj,istat
INTEGER,intent(in) :: mstat
integer indexb(rcpara%nmax), indexd(rcpara%nmax), indexg(rcpara%nmax)
character*(*) filename
real(kind=JPRM) :: rasocorr(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM),optional :: rasobreak(rcpara%nmax,rcpara%pmax,rcpara%parmax),breakuncertainty(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: hilfcorr(rcpara%mmax,rcpara%pmax,rcpara%parmax,3)
logical jump
REAL (kind=JPRM):: lon, lat
INTEGER :: numpar, maxatts
INTEGER, ALLOCATABLE :: datum(:,:)
REAL (kind=JPRM), ALLOCATABLE :: arr(:,:,:,:)
CHARACTER*12, ALLOCATABLE :: names(:)
CHARACTER*50, OPTIONAL:: stname
!Berechne Parameteranzahl
lserial=.true.
if(lserial) then
  istat=rcpara%statmax+1
else
  istat=mstat
endif
!istat=rcpara%statmax+1
!istat=mstat  !not +4000 since it calls gen_nc
numpar=1
IF (present(rasobreak)) numpar=numpar+1
IF (present(breakuncertainty)) numpar=numpar+1

!!  print*,filename
!  open(iunit,file=filename,form='unformatted',err=120)
!Initialisieren

write(*,*) 'writing x'//trim(filename)//'x'
  err=0
hilfcorr=rcpara%miss_val
bi=0
di=0
gi=0
!Korrekturen in hilfcorr speichern
  bi=1
  hilfcorr(1,:,:,1)=rasocorr(1,:,:)
  indexb=rcpara%miss_val
  indexb(1)=1
  do i=2,rcpara%nmax-1
    jump=.false.
    do ipar=1,rcpara%parmax
      ip=1
      do while(ip .le. rcpara%pmax-2 .and. .not. jump)
        jump=abs(rasocorr(i,ip,ipar)-rasocorr(i+1,ip,ipar)).gt. 0.001
        ip=ip+1
      enddo
    enddo
    if(jump) then
      bi=bi+1
      hilfcorr(bi,:,:,1)=rasocorr(i+1,:,:)
      if(present(rasobreak)) then
        hilfcorr(bi,:,:,2)=rasobreak(i+1,:,:)
      endif
      if(present(breakuncertainty)) then
        hilfcorr(bi,:,:,3)=breakuncertainty(i+1,:,:)
      endif
      indexb(bi)=i+1
    endif
    if(bi .eq. rcpara%mmax) write(*,*) 'Warning: Maximum allowed break number ',rcpara%mmax,' reached'
  enddo


maxatts=4
ALLOCATE(arr(bi,rcpara%pmax,rcpara%parmax,numpar), names(numpar), datum(bi,numpar))
ALLOCATE (attribs%attname(maxatts), attribs%long_name(numpar), attribs%units(numpar), attribs%missing_val(numpar), attribs%valid_range(numpar,2))
datum=rcpara%miss_val

attribs%title='radiosonde corrections'
attribs%attname=(/'long_name','units','missing_value','valid_range'/)
attribs%units='K'
attribs%missing_val=rcpara%miss_val
attribs%valid_range=rcpara%miss_val

!hilfcorr in arr speichern
arr(:,:,:,1)=hilfcorr(1:bi,:,:,1)
names(1)='rasocorr'
attribs%long_name(1)='raso_correct'
datum(:,1)= indexb(1:bi)
IF (numpar .ge. 2 ) THEN
	arr(:,:,:,2)=hilfcorr(1:vals,:,:,2)
	names(2)='rasobreak'
	attribs%long_name(2)='rasobreak'
	datum(:,2)=datum(:,1)
END IF
IF (numpar==3) THEN
	arr(:,:,:,3)=hilfcorr(1:vals,:,:,3)
	names(3)='breakuncert'
	attribs%long_name(3)='breakuncertainty'
	datum(:,3)=datum(:,1)
END IF
!!$omp critical
	CALL gen_nc(rcpara,istat,filename,numpar,maxatts,bi, arr, names, attribs, lat, lon, datum, numpar,stname=stname)  !in this file line 278
!!$omp end critical
  
  return

DEALLOCATE(arr,datum,names,attribs%attname,attribs%long_name,attribs%units,attribs%missing_val,attribs%valid_range)
120 continue
    !!$ call omp_set_lock(omp_lp)
print*,'could not open file ',filename    !!$ call omp_unset_lock(omp_lp)
  err=1

  return

end subroutine write_sonde_corr_daily_nc

SUBROUTINE read_odb_nc(filename,rcpara,mstat,err3,tm,tfgm,tbcm,fflags,tanm,e20c0,stype,hours,stname,bad_intervals,tgps,alt)
  !leseroutine zum lesen der odb-netcdf files
  USE netcdf
  USE rfmod
  IMPLICIT NONE
  type(rasocor_namelist) rcpara
  INTEGER :: status, e20status,gpsstatus,ncid,mask(4) !ID des nc - Files, wird beim oeffnen zugewiesen
  INTEGER :: datumvarid, fg_depvarid, biascorr_varid, tempvarid, flagsvarid, timid,an_depvarid,s_typevarid,e20c0_depvarid,tgps_depvarid,hoursvarid
  CHARACTER :: filename*80,cunits*80 !name+Pfad des nc-files
  CHARACTER*30 :: timedimname
  REAL (kind=JPRM), ALLOCATABLE :: temp(:,:,:), fg_depar(:,:,:), bias_cor(:,:,:), flags(:,:,:), flags2(:,:,:,:),an_dep(:,:,:),e20c0_dep(:,:,:),tgps_dep(:,:,:),sonde_type(:)
  INTEGER, ALLOCATABLE :: datum(:),lhours(:,:)
  REAL(kind=JPRM) :: tm(rcpara%nmax,rcpara%pmax,rcpara%parmax), tfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
  REAL(kind=JPRM) :: tmmon(rcpara%mmax,rcpara%pmax,rcpara%parmax),anom(rcpara%mmax),climate(12)
  INTEGER,OPTIONAL :: hours(rcpara%nmax,rcpara%parmax)
  REAL(kind=JPRM),OPTIONAL :: tanm(rcpara%nmax,rcpara%pmax,rcpara%parmax),stype(rcpara%nmax)
  REAL(kind=JPRM),OPTIONAL :: e20c0(rcpara%nmax,rcpara%pmax,rcpara%parmax)
  REAL(kind=JPRM),OPTIONAL :: tgps(rcpara%nmax,rcpara%pmax,rcpara%parmax)
  REAL(kind=JPRM),optional :: tbcm(rcpara%nmax,rcpara%pmax,rcpara%parmax), fflags(rcpara%nmax,rcpara%pmax,rcpara%parmax)
  REAL(kind=JPRM),optional :: alt
  INTEGER :: err3,jj,istat,i,j,k,id,im,iy,ioffset,ip,ipar,l,statnr,ios,ystart,ystop
  INTEGER,intent(in) :: mstat
  INTEGER :: ii(1) !zaehlvariable
  INTEGER :: datecount !zaehlvariable
  INTEGER :: datum2 !Datum aus den rcpara-arrays year, month und day
  INTEGER :: maxtime, datumyear,altvarid,ofg_depvarid
  INTEGER,OPTIONAL :: bad_intervals(:,:)
  CHARACTER*5 :: uname
  CHARACTER*50,optional :: stname
  logical:: lanp
  !ALLOCATE(tm(rcpara%nmax,rcpara%pmax,rcpara%parmax), tanm(rcpara%nmax,rcpara%pmax,rcpara%parmax), tfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax), itfg12m(rcpara%nmax,rcpara%pmax,rcpara%parmax))
  !initialize
  lserial=.true.
  if(lserial) then
     istat=rcpara%statmax+1
  else
     istat=mstat+4000
  endif
  !istat=rcpara%statmax+1
  !istat=mstat+4000
  if(istat .eq. 0) then
     write(*,*) 'istat zero!'
  endif
  tm=rcpara%miss_val
  tfgm=rcpara%miss_val
  if(present(tbcm)) tbcm=rcpara%miss_val
  if(present(fflags)) fflags=rcpara%miss_val
  if(present(tanm)) tanm=rcpara%miss_val
  if(present(e20c0)) e20c0=rcpara%miss_val
  if(present(tgps)) tgps=rcpara%miss_val
  if(present(stype)) stype=rcpara%miss_val
  if(present(alt)) alt=rcpara%miss_val
  !write(*,*) istat
  !$ call omp_set_lock(omp_lp(istat))
!!$omp critical
  status=NF90_OPEN(TRIM(filename),NF90_NOWRITE,ncid) !nc-file oeffnen
!!$omp end critical
  !write(*,*) 'rmodulo',mod(ncid,65536),ncid/65536
  IF (status /= NF90_NOERR) THEN
     WRITE(*,*) 'Error opening netcdf!', TRIM(filename)
     status=NF90_CLOSE(ncid)
     write(*,*) nf90_strerror(status)
     err3=1
     !$ call omp_unset_lock(omp_lp(istat))
     RETURN
  END IF
!!$ call omp_unset_lock(omp_lp(istat))

  !get time dimension ID
!!!$omp critical
  status=NF90_INQ_DIMID(ncid,'time',timid)
!!!$omp end critical
  IF (status /= NF90_NOERR) THEN
     WRITE(*,*) 'Error inquire dimension ID', TRIM(filename), ' maybe file empty'
     status=NF90_CLOSE(ncid)
     IF (status /= NF90_NOERR) WRITE(*,*) 'Error closeing netcdf!', TRIM(filename)
     !$ call omp_unset_lock(omp_lp(istat))
     RETURN
  END IF
  !get time dimenstion length
!!!$omp critical
  status=NF90_INQUIRE_DIMENSION(ncid,timid,timedimname,maxtime)
!!!$omp end critical
  IF (status /= NF90_NOERR) THEN
     WRITE(*,*) 'Error inquire time dimension length', TRIM(filename), 'maybe file empty'
     status=NF90_CLOSE(ncid)
     !$ call omp_unset_lock(omp_lp(istat))
     RETURN
  END IF
  !get variable IDs
!!!$omp critical
  status=NF90_INQ_VARID(ncid,'datum',datumvarid)
!!!$omp end critical
  IF (status /= NF90_NOERR) THEN
     WRITE(*,*) 'Error inquire variable ID', TRIM(filename), 'maybe file empty'
     status=NF90_CLOSE(ncid)
     !$ call omp_unset_lock(omp_lp(istat))
     RETURN
  END IF
  status=NF90_INQ_VARID(ncid,'hours',hoursvarid)
!!!$omp end critical
  IF (status /= NF90_NOERR) THEN
     WRITE(*,*) 'Error inquire variable ID', TRIM(filename), 'maybe file empty'
     status=NF90_CLOSE(ncid)
     !$ call omp_unset_lock(omp_lp(istat))
     RETURN
  END IF
!!!$omp critical
  status=NF90_INQ_VARID(ncid,'temperatures', tempvarid)
  status=NF90_INQ_VARID(ncid,trim(rcpara%fgdepname), fg_depvarid)
  status=NF90_INQ_VARID(ncid,'fg_dep', ofg_depvarid)
  e20status=0
  gpsstatus=0
  if(present(tbcm)) status=NF90_INQ_VARID(ncid,'bias', biascorr_varid)
  if(status /=0) WRITE(*,*) trim(filename)//'- bias correction missing'
  if(present(fflags)) status=NF90_INQ_VARID(ncid,'flags', flagsvarid)
  if(status /=0) WRITE(*,*) trim(filename)//'- quality control flags missing'
  if(present(tanm)) status=NF90_INQ_VARID(ncid,'an_dep', an_depvarid)
  if(status /=0) WRITE(*,*) trim(filename)//'- analysis departures missing'
  if(present(e20c0)) e20status=NF90_INQ_VARID(ncid,'e20c_0', e20c0_depvarid)
  if(e20status /=0) THEN
     WRITE(*,*) trim(filename)//'- e20c_0 departures missing'
  endif
  if(index(rcpara%initial_adjust,'wet')>0) then
     gpsstatus=NF90_INQ_VARID(ncid,'erai_fggpswetdep', tgps_depvarid)
  else
     gpsstatus=NF90_INQ_VARID(ncid,'erai_fggpsdep', tgps_depvarid)
  endif
  if(present(alt)) status=NF90_INQ_VARID(ncid,'alt', altvarid)
  if(gpsstatus /=0) THEN
     WRITE(*,*) trim(filename)//'- gps departures missing'
  endif
   if(present(stype)) status=NF90_INQ_VARID(ncid,'s_type', s_typevarid)
  if(status /=0) WRITE(*,*) trim(filename)//'- sonde types missing'
  !allocate variables with time dimension
  ALLOCATE(datum(maxtime))
  status=NF90_GET_VAR(ncid,datumvarid, datum)
!!!$omp end critical
  IF (status /= NF90_NOERR) then
     WRITE(*,*) 'Error reading date from nc'
     status=NF90_CLOSE(ncid)
     IF (status /= NF90_NOERR) WRITE(*,*) 'Error closeing netcdf!', TRIM(filename)
     !$ call omp_unset_lock(omp_lp(istat))
     DEALLOCATE(datum)
     RETURN
  END IF
  uname="units"
  status=NF90_GET_ATT(ncid,datumvarid,"units",cunits) !'days since '//cstartyear//'-01-01 0:0:0'

  if(present(stname)) then
     status=NF90_GET_ATT(ncid,NF90_GLOBAL,"Stationnname",stname)
     if (status .ne. 0) then
        status=NF90_GET_ATT(ncid,NF90_GLOBAL,'Stationname',stname) 
     endif
  endif

  i=index(cunits,'-')-4
  cunits=cunits(i:i+9)
  read(cunits,'(I4,1x,I2,1x,I2)') iy,im,id
  ioffset=toindex(iy*10000+im*100+id,rcpara)     


  if(present(hours)) then
     ALLOCATE(lhours(maxtime,rcpara%parmax))
     status=NF90_GET_VAR(ncid,hoursvarid, lhours)
!!!$omp end critical
     IF (status /= NF90_NOERR) then
        WRITE(*,*) 'Error reading date from nc'
	status=NF90_CLOSE(ncid)
	IF (status /= NF90_NOERR) WRITE(*,*) 'Error closeing netcdf!', TRIM(filename)
 !$ call omp_unset_lock(omp_lp(istat))
        DEALLOCATE(lhours)
	RETURN
     ELSE

        hours(datum,:)=lhours

        DEALLOCATE(lhours)
     END IF
  endif

  !allocate variables with time dimension
  ALLOCATE(temp(maxtime,rcpara%pmax,rcpara%parmax))
  !initialize
  !temp=rcpara%miss_val
  !fg_depar=rcpara%miss_val
  !get variables
!!!$omp critical
  status=NF90_GET_VAR(ncid,tempvarid,temp)
  where (isnan(temp))
     temp=rcpara%miss_val
  endwhere
!!!$omp end critical
  !WRITE(*,*) shape(temp), size(shape(temp))
  IF (status /= NF90_NOERR) THEN
     WRITE(*,*) TRIM(NF90_STRERROR(status)), '   Error reading temperature from nc'
     status=NF90_CLOSE(ncid)
     IF (status /= NF90_NOERR) WRITE(*,*) 'Error closeing netcdf!', TRIM(filename)
     !$ call omp_unset_lock(omp_lp(istat))
     DEALLOCATE(temp,datum)
     RETURN
  END IF

  !initialize
  !allocate variables with time dimension
  ALLOCATE(fg_depar(maxtime,rcpara%pmax,rcpara%parmax))
  status=NF90_GET_VAR(ncid,fg_depvarid,fg_depar)
  i=toindex(20170101,rcpara)
  if(rcpara%fgdepname .ne. 'fg_dep' .and. datum(maxtime)>=i) then
    ALLOCATE(an_dep(maxtime,rcpara%pmax,rcpara%parmax))
    status=NF90_GET_VAR(ncid,ofg_depvarid,an_dep)
    j=maxtime
    do while(datum(j)>=i)
      j=j-1
    enddo
      
    fg_depar(j:,:,:)=-an_dep(j:,:,:)
    DEALLOCATE(an_dep)
  endif
  where (isnan(fg_depar))
     fg_depar=rcpara%miss_val
  endwhere
  if (rcpara%fgdepname .ne. 'fg_dep') then
     WHERE (fg_depar .ne. rcpara%miss_val)
        fg_depar=-fg_depar
     ELSEWHERE
        fg_depar= rcpara%miss_val
     ENDWHERE
  endif

  IF (status /= NF90_NOERR) then
     WRITE(*,*) 'Error reading fg_departures from nc'
     status=NF90_CLOSE(ncid)
     IF (status /= NF90_NOERR) WRITE(*,*) 'Error closeing netcdf!', TRIM(filename)
     !$ call omp_unset_lock(omp_lp(istat))
     DEALLOCATE(temp,datum,fg_depar)
     RETURN
  END IF

  if(present(tbcm)) then
     ALLOCATE(bias_cor(maxtime,rcpara%pmax,rcpara%parmax))
     !  bias_cor=rcpara%miss_val
     status=NF90_GET_VAR(ncid,biascorr_varid,bias_cor)
     where (isnan(bias_cor))
        bias_cor=rcpara%miss_val
     endwhere
     IF (status /= NF90_NOERR) WRITE(*,*) 'Error reading bias correction from nc'
  endif
  if(present(fflags) .or. rcpara%qc .EQ. 'B') then
     ALLOCATE(flags(maxtime,rcpara%pmax,rcpara%parmax), flags2(maxtime,rcpara%pmax,rcpara%parmax,6))
     !  flags=rcpara%miss_val
     flags2=rcpara%miss_val
     status=NF90_GET_VAR(ncid,flagsvarid,flags)
     IF (status /= NF90_NOERR) WRITE(*,*) 'Error reading flags from nc'
  endif
  if(present(tanm)) then
     ALLOCATE(an_dep(maxtime,rcpara%pmax,rcpara%parmax))
     !  an_dep=rcpara%miss_val
     status=NF90_GET_VAR(ncid,an_depvarid,an_dep)
     where (isnan(an_dep))
        an_dep=rcpara%miss_val
     endwhere
     IF (status /= NF90_NOERR) WRITE(*,*) 'Error reading analysis departures from nc'
  endif
  if(present(e20c0) .and. e20status==0) then
     ALLOCATE(e20c0_dep(maxtime,rcpara%pmax,rcpara%parmax))
     !  an_dep=rcpara%miss_val
     status=NF90_GET_VAR(ncid,e20c0_depvarid,e20c0_dep)
     IF (status /= NF90_NOERR) WRITE(*,*) 'Error reading analysis departures from nc'
  endif
  if(present(alt)) then
     status=NF90_GET_VAR(ncid,altvarid,alt)
     if (status/=0) write (*,*) 'No station height'
  endif
  if(present(tgps) .and. gpsstatus==0) then
     ALLOCATE(tgps_dep(maxtime,rcpara%pmax,rcpara%parmax))
     !  an_dep=rcpara%miss_val
     status=NF90_GET_VAR(ncid,tgps_depvarid,tgps_dep)
     where (isnan(tgps_dep))
        tgps_dep=rcpara%miss_val
     endwhere
     IF (status /= NF90_NOERR) then
       WRITE(*,*) 'Error reading gps departures from nc'
     ENDIF
  endif
  if(present(stype)) then
     ALLOCATE(sonde_type(maxtime))
     !  an_dep=rcpara%miss_val
     status=NF90_GET_VAR(ncid,s_typevarid,sonde_type)
     IF (status /= NF90_NOERR) WRITE(*,*) 'Error reading analysis departures from nc'
  endif
  status=NF90_CLOSE(ncid)
  !$ call omp_unset_lock(omp_lp(istat))

  if(present(fflags) .or. rcpara%qc .EQ. 'B') then
     flags2(:,:,:,1)=FLOOR(flags/1000000.)
     flags2(:,:,:,2)=FLOOR((flags-flags2(:,:,:,1)*1000000)/100000.)
     flags2(:,:,:,3)=FLOOR((flags-flags2(:,:,:,1)*1000000-flags2(:,:,:,2)*100000)/10000.)
     flags2(:,:,:,4)=FLOOR((flags-flags2(:,:,:,1)*1000000-flags2(:,:,:,2)*100000-flags2(:,:,:,3)*10000)/1000.)
     flags2(:,:,:,5)=FLOOR((flags-flags2(:,:,:,1)*1000000-flags2(:,:,:,2)*100000-flags2(:,:,:,3)*10000-flags2(:,:,:,4)*1000)/100.)
  endif

  IF (status /= NF90_NOERR) WRITE(*,*) 'Error closeing netcdf!', TRIM(filename)
  ! temp ist (obs+biascor), wir wollen obs 
  if(present(tbcm)) then
     WHERE (temp .NE. rcpara%miss_val .AND. bias_cor .NE. rcpara%miss_val)
	temp=temp-bias_cor
     END WHERE
  endif

  ! fg_depar ist (obs+biascor)-bg, wir wollen bg-obs 
  if(present(tbcm)) then
     WHERE (fg_depar .NE.rcpara%miss_val .AND. bias_cor .NE. rcpara%miss_val)
	fg_depar=-(fg_depar-bias_cor)
     END WHERE
  endif
  WHERE (ABS(fg_depar) .GT. 20.)
     temp=rcpara%miss_val
     fg_depar=rcpara%miss_val
  END WHERE
  IF (rcpara%qc .EQ. 'B') THEN
     WHERE (flags2(:,:,:,5) .GE. 1 )
        temp=rcpara%miss_val
        fg_depar=rcpara%miss_val
     END WHERE
  END IF


  IF(PRESENT(TANM) .AND. PRESENT(TBCM)) THEN
     ! an_depar ist (obs+biascor)-an, wir wollen an-obs 
     WHERE (an_dep .NE.rcpara%miss_val .AND. bias_cor .NE. rcpara%miss_val)
	an_dep=-(an_dep-bias_cor)
     END WHERE
     WHERE (ABS(an_dep) .GT. 20.)
	an_dep=rcpara%miss_val
     END WHERE
     IF (rcpara%qc .EQ. 'B') THEN
	WHERE (flags2(:,:,:,5) .GE. 1 )
           an_dep=rcpara%miss_val
	END WHERE
     END IF
  ENDIF
  IF(PRESENT(e20c0) .and. e20status==0) THEN
     ! an_depar ist (obs+biascor)-an, wir wollen an-obs 
     WHERE (e20c0_dep .NE.rcpara%miss_val)
	e20c0_dep=-e20c0_dep
     END WHERE
     !WHERE (ABS(e20c0_dep) .GT. 20.)
     !	e20c0_dep=rcpara%miss_val
     !END WHERE
  ENDIF
  IF(PRESENT(tgps) .and. gpsstatus==0) THEN
     ! an_depar ist (obs+biascor)-an, wir wollen an-obs 
     WHERE (tgps_dep .NE.rcpara%miss_val)
	tgps_dep=-tgps_dep
     END WHERE
     !WHERE (ABS(e20c0_dep) .GT. 20.)
     !	e20c0_dep=rcpara%miss_val
     !END WHERE
  ENDIF
  !tempdiff=temp(:,:,2)-temp(:,:,1)

  datecount=1 !zaehlt die Tage mit Daten
  !datumyear=FLOOR(datum(1)/10000.)
  !ii=MINLOC(rcpara%year,MASK= rcpara%year .GE. datumyear) !zaehlt den Zeitindex, beginnt im ersten Jahr

  jj=datum(1)
  jj=count(datum(1:maxtime) .le. 0 .or. datum(1:maxtime) .gt. rcpara%nmax)
  if(jj .gt. 0) then
     write(*,*) 'datum',maxtime,datum(1),datum(maxtime),jj
     write(*,*) trim(filename),' date read failure'
     !  deallocate(temp, datum, fg_depar)
     !  if(present(tbcm)) DEALLOCATE(bias_cor)
     !if(present(fflags))  DEALLOCATE( flags, flags2)
     !if(present(tanm))   deallocate(an_dep)
     !  return
  endif

  do k=1,rcpara%parmax
     do j=1,rcpara%pmax
        do i=1,maxtime
           if(datum(i)+ioffset.gt.0 .and. datum(i)+ioffset.le.rcpara%nmax) then
              tm(datum(i)+ioffset,j,k)=temp(i,j,k)
              tfgm(datum(i)+ioffset,j,k)=fg_depar(i,j,k)
           endif
        enddo
     enddo
  enddo

!!$  call makemonthly(rcpara,tm,tmmon,5) !
!!$  !      anom=tmmon
!!$  i=datum(1)+ioffset
!!$  ystart=(rcpara%year(i)-1900)*12+1
!!$  i=datum(maxtime)+ioffset-1
!!$  ystop=(rcpara%year(i)-1900+1)*12
!!$  print*,maxtime,datum(maxtime),i,ystop+1-ystart
!!$  
!!$
!!$  do ipar=1,rcpara%parmax
!!$     do ip=1,rcpara%pmax
!!$        i=ystop+1-ystart
!!$        jj=12
!!$        if (ystop>rcpara%mmax) then
!!$           write(*,*) 'ALARM',filename(1:20),i,datum(maxtime),ystop
!!$           exit
!!$        endif
!!$        call anomaly(tmmon(ystart:ystop,ip,ipar),anom(ystart:ystop),climate,i,jj,rcpara%miss_val,10) !in file rfmod.f90 line 2311
!!$        if(count(climate .ne. rcpara%miss_val)>5) then
!!$           do i=1,maxtime
!!$              if( climate(rcpara%month(datum(i)+ioffset)) .ne. rcpara%miss_val .and. tm(datum(i)+ioffset,ip,ipar).ne. rcpara%miss_val) then
!!$                 tm(datum(i)+ioffset,ip,ipar)=tm(datum(i)+ioffset,ip,ipar)-climate(rcpara%month(datum(i)+ioffset))
!!$              else
!!$                 tm(datum(i)+ioffset,ip,ipar)=rcpara%miss_val
!!$              endif
!!$           enddo
!!$        endif
!!$     enddo
!!$  enddo


  if(present(tbcm)) then 
     do k=1,rcpara%parmax
        do j=1,rcpara%pmax
           do i=1,maxtime
              if(datum(i)+ioffset.gt.0 .and. datum(i)+ioffset.le.rcpara%nmax) then
                 tbcm(datum(i)+ioffset,j,k)=bias_cor(i,j,k)
              endif
           enddo
        enddo
     enddo
     DEALLOCATE(bias_cor)
  endif
  if(present(fflags)) then 
     do k=1,rcpara%parmax
        do j=1,rcpara%pmax
           do i=1,maxtime
              if(datum(i)+ioffset.gt.0 .and. datum(i)+ioffset.le.rcpara%nmax) then
                 fflags(datum(i)+ioffset,j,k)=flags(i,j,k)
              endif
           enddo
        enddo
     enddo
     DEALLOCATE( flags, flags2)
  endif

  if(present(tanm)) then
     do k=1,rcpara%parmax
        do j=1,rcpara%pmax
           do i=1,maxtime
              if(datum(i)+ioffset.gt.0 .and. datum(i)+ioffset.le.rcpara%nmax) then
                 tanm(datum(i)+ioffset,j,k)=an_dep(i,j,k)
              endif
           enddo
        enddo
     enddo
     deallocate(an_dep)
  endif
  if(present(e20c0) .and. e20status==0) then
     do k=1,rcpara%parmax
        do j=1,rcpara%pmax
           do i=1,maxtime
              if(datum(i)+ioffset.gt.0 .and. datum(i)+ioffset.le.rcpara%nmax) then
                 e20c0(datum(i)+ioffset,j,k)=e20c0_dep(i,j,k)
              endif
           enddo
        enddo
     enddo
     deallocate(e20c0_dep)
  endif
  if(present(tgps) .and. gpsstatus==0) then
     do k=1,rcpara%parmax
        do j=1,rcpara%pmax
           do i=1,maxtime
              if(datum(i)+ioffset.gt.0 .and. datum(i)+ioffset.le.rcpara%nmax) then
                 tgps(datum(i)+ioffset,j,k)=tgps_dep(i,j,k)
              endif
           enddo
        enddo
     enddo
     deallocate(tgps_dep)
  endif
  if(present(stype)) then
     !print*,'sondetype',sonde_type
     do i=1,maxtime
        if(datum(i)+ioffset.gt.0 .and. datum(i)+ioffset.le.rcpara%nmax) then
           stype(datum(i)+ioffset)=sonde_type(i)
        endif
     enddo
     deallocate(sonde_type)
  endif
  DEALLOCATE(temp, datum, fg_depar)
  if(present(bad_intervals)) then
     i=len(trim(filename))
     ios=1
     do while (i>6 .and. ios .ne. 0)
       read(filename(i-6:i),'(I6)',iostat=ios) statnr
       i=i-1
     enddo
     lanp=present(tanm)
     do l=1,size(bad_intervals,2)
        if (bad_intervals(1,l).eq. statnr ) then
           do ipar=1,rcpara%parmax
              if (bad_intervals(5,l).eq. ipar ) then
                 do ip=1,bad_intervals(4,l)
                    do i=bad_intervals(2,l),bad_intervals(3,l)
                       tm(i,ip,ipar)=rcpara%miss_val
                       tfgm(i,ip,ipar)=rcpara%miss_val
                       if(lanp) tanm(i,ip,ipar)=rcpara%miss_val
                    enddo
                 enddo
              endif
           enddo
        endif
     enddo
  endif
  err3=0

  RETURN
END SUBROUTINE read_odb_nc

SUBROUTINE read_odb_nc_json(filename,rcpara,cstat,ipl,itime,version,selector,exper,nc_files,inc,adj)
!leseroutine zum lesen der odb-netcdf files
USE netcdf
USE rfmod
IMPLICIT NONE
type(rasocor_namelist) rcpara
INTEGER :: status, ncid,mask(4) !ID des nc - Files, wird beim oeffnen zugewiesen
INTEGER :: datumvarid, fg_depvarid, biascorr_varid, tempvarid, flagsvarid, timid,an_depvarid
CHARACTER :: filename*80 !name+Pfad des nc-files
CHARACTER*30 :: timedimname
REAL (kind=JPRM), ALLOCATABLE :: temp(:,:,:)
INTEGER, ALLOCATABLE :: datum(:)
INTEGER :: err3,jj,istat,i,n,j,lpos
INTEGER,intent(in) :: ipl,itime(:),inc
character*5,intent(in) :: cstat
INTEGER :: ii(1) !zaehlvariable
INTEGER :: datecount !zaehlvariable
INTEGER :: datum2 !Datum aus den rcpara-arrays year, month und day
INTEGER :: maxtime, datumyear,iunit
CHARACTER(*) :: selector(4)
CHARACTER(*),intent(in) :: version,exper,nc_files(:),adj
character*2 :: suff(2)
logical :: lfound

istat=rcpara%statmax+1
if(istat .eq. 0) then
  write(*,*) 'istat zero!'
endif

!$ call omp_set_lock(omp_lp(istat))
!!!$omp critical
status=NF90_OPEN(TRIM(filename),NF90_NOWRITE,ncid) !nc-file oeffnen
!!!$omp end critical
IF (status /= NF90_NOERR) THEN
	WRITE(*,*) 'Error opening netcdf!', TRIM(filename)
	status=NF90_CLOSE(ncid)
        err3=1
!$ call omp_unset_lock(omp_lp(istat))
	RETURN
END IF
!!$ call omp_unset_lock(omp_lp(istat))

!get time dimension ID
status=NF90_INQ_DIMID(ncid,'time',timid)
IF (status /= NF90_NOERR) THEN
	WRITE(*,*) 'Error inquire dimension ID', TRIM(filename), ' maybe file empty'
	status=NF90_CLOSE(ncid)
	IF (status /= NF90_NOERR) WRITE(*,*) 'Error closeing netcdf!', TRIM(filename)
!$ call omp_unset_lock(omp_lp(istat))
	RETURN
END IF
!get time dimenstion length
status=NF90_INQUIRE_DIMENSION(ncid,timid,timedimname,maxtime)
!maxtime=20
IF (status /= NF90_NOERR) THEN
	WRITE(*,*) 'Error inquire time dimension length', TRIM(filename), 'maybe file empty'
	status=NF90_CLOSE(ncid)
!$ call omp_unset_lock(omp_lp(istat))
	RETURN
END IF

! NOTE: Before we write anything to JSON file, we check if a varID matches an
! element of the selector array
lfound=.false.
n=size(selector)
do i=1,n
  if(selector(i) .eq. "") exit
enddo
n=i-1
do i=1,n
 status=NF90_INQ_VARID(ncid,trim(selector(i)), tempvarid)
 if(status ==0) then 
   lfound=.true.
   write(*,*) trim(nc_files(inc)),' ',trim(selector(i)),' matches'
 endif
enddo

iunit=1
if(inc == 1) then
  open(iunit,file='/home/srvx2/leo/public_html/jquery/flot/examples/out.json')
  write(iunit,*) '{'
endif

if (.not. lfound) then
  write(*,*) trim(nc_files(inc)),' no match found'
  if(inc == size(nc_files)) then
    write(iunit,*) '}'
    close(iunit)
  endif
!$ call omp_unset_lock(omp_lp(istat))
  return
endif

if(inc > 1) write(iunit,*) ','
write(iunit,*) '"'//trim(nc_files(inc))//'": {'
write(iunit,*) '"pressure_level":',rcpara%plevs(ipl),','
write(iunit,*) '"hour":',(itime(1)-1)*12,','
write(iunit,*) '"statid":"',cstat,'",'
write(iunit,*) '"startdate":',rcpara%startdate,','
write(iunit,*) '"Adjustment_Method": "'//trim(adj)//'",'
write(iunit,*) '"Version":','"'//version(2:len(trim(version)))//'"',','
write(iunit,*) '"Experiment":','"'//trim(exper)//'"',','

!get variable IDs
status=NF90_INQ_VARID(ncid,'datum',datumvarid)
IF (status /= NF90_NOERR) THEN
	WRITE(*,*) 'Error inquire variable ID', TRIM(filename), 'maybe file empty'
	status=NF90_CLOSE(ncid)
!$ call omp_unset_lock(omp_lp(istat))
	RETURN
END IF
ALLOCATE(datum(maxtime))
status=NF90_GET_VAR(ncid,datumvarid, datum)
IF (status /= NF90_NOERR) WRITE(*,*) 'Error reading date from nc'
call jsonwrite_i(iunit,"time","datum",maxtime,datum)

ALLOCATE(temp(maxtime,rcpara%pmax,rcpara%parmax))
n=size(selector)
do i=1,n
  if(selector(i) .eq. "") exit
enddo
n=i-1
do i=1,n
 status=NF90_INQ_VARID(ncid,trim(selector(i)), tempvarid)
 if(status ==0) then 
  temp=rcpara%miss_val
!get variables
  status=NF90_GET_VAR(ncid,tempvarid,temp)
!WRITE(*,*) shape(temp), size(shape(temp))
  IF (status /= NF90_NOERR) THEN
	WRITE(*,*) TRIM(NF90_STRERROR(status)), '   Error reading temperature from nc'
	status=NF90_CLOSE(ncid)
	IF (status /= NF90_NOERR) WRITE(*,*) 'Error closeing netcdf!', TRIM(filename)
!$ call omp_unset_lock(omp_lp(istat))
        DEALLOCATE(temp)
	RETURN
  END IF
  suff=(/"  ","_2"/)
  j=1
  do while(itime(j) .ne. -1 .and. j<3)
  lpos=index(trim(selector(i)),'-')
  if(lpos>=0) selector(i)(lpos:lpos)='_'
  lpos=index(trim(selector(i)),' ')
  if(lpos>=0) selector(i)(lpos:lpos)='_'

  WRITE(iunit,*) ','
    call jsonwrite_f(iunit,"",trim(selector(i))//trim(suff(j)),maxtime,temp(:,ipl,itime(j)))
    j=j+1
  enddo
 endif
enddo
status=NF90_CLOSE(ncid)
write(iunit,*) '}'
if(inc == size(nc_files)) then
  write(iunit,*) '}'
  close(iunit)
endif
!$ call omp_unset_lock(omp_lp(istat))

IF (status /= NF90_NOERR) WRITE(*,*) 'Error closeing netcdf!', TRIM(filename)


RETURN
END SUBROUTINE read_odb_nc_json

subroutine jsonwrite_i(iunit,d_label,var_label,maxtime,datum)
implicit none
character*(*) d_label,var_label
integer maxtime,iunit,i
integer datum(maxtime)

if(d_label .ne. '') write(iunit,*) '"'//d_label//'":',maxtime,','
write(iunit,*) '"'//var_label//'": ['
do i=1,maxtime-1
  write(iunit,'(I5,A1)') datum(i),','
enddo
  write(iunit,'(I5)') datum(i)
write(iunit,*) ']'
return
end subroutine jsonwrite_i

subroutine jsonwrite_f(iunit,d_label,var_label,maxtime,datum)
implicit none
character*(*) d_label,var_label
integer maxtime,iunit,i
real(kind=4) datum(maxtime)

if(d_label .ne. '') write(iunit,*) '"'//d_label//'":',maxtime,','
write(iunit,*) '"'//var_label//'": ['
do i=1,maxtime-1
  write(iunit,'(F7.2,A1)') datum(i),','
enddo
  write(iunit,'(F7.2)') datum(i)
write(iunit,*) ']'
return
end subroutine jsonwrite_f

subroutine jsonwrite_s(iunit,d_label,var_label,maxtime,datum)
implicit none
character*(*) d_label,var_label
integer maxtime,iunit,i
character*(*) datum(maxtime)

if(d_label .ne. '') write(iunit,*) '"'//d_label//'"',maxtime,','
write(iunit,*) '"'//var_label//'": ['
do i=1,maxtime-1
  write(iunit,*) datum(i),','
enddo
  write(iunit,*) datum(i)
write(iunit,*) ']'
return
end subroutine jsonwrite_s


SUBROUTINE read_bgcorr_nc(filename,rcpara,mstat,err3,tm,names,tfgm)
!leseroutine zum lesen der odb-netcdf files
USE netcdf
USE rfmod
IMPLICIT NONE
type(rasocor_namelist) rcpara
INTEGER :: status, ncid,mask(4) !ID des nc - Files, wird beim oeffnen zugewiesen
INTEGER :: datumvarid, fg_depvarid, biascorr_varid, tempvarid, flagsvarid, timid
INTEGER,INTENT(IN) :: mstat
CHARACTER :: filename*80 !name+Pfad des nc-files
CHARACTER*12 :: names(2)
CHARACTER*30 :: timedimname
REAL (kind=JPRM), ALLOCATABLE :: temp(:,:,:), fg_depar(:,:,:)
INTEGER, ALLOCATABLE :: datum(:)
REAL(kind=JPRM) :: tm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
REAL(kind=JPRM),optional :: tfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
INTEGER :: err3,jj,istat
INTEGER :: ii(1) !zaehlvariable
INTEGER :: datecount !zaehlvariable
INTEGER :: datum2 !Datum aus den rcpara-arrays year, month und day
INTEGER :: maxtime, datumyear
!ALLOCATE(tm(rcpara%nmax,rcpara%pmax,rcpara%parmax), tanm(rcpara%nmax,rcpara%pmax,rcpara%parmax), tfgm(rcpara%nmax,rcpara%pmax,rcpara%parmax), itfg12m(rcpara%nmax,rcpara%pmax,rcpara%parmax))
!initialize
lserial=.true.
if(lserial) then
  istat=rcpara%statmax+1
else
  istat=mstat+4000
endif
!istat=rcpara%statmax+1
!istat=mstat+4000
if(istat .eq. 0) then
  write(*,*) 'istat zero!'
endif
tm=rcpara%miss_val
!!!!!!!!!!!!!!
tm=0
return
!!!!!!!!!!!!!!
if(present(tfgm)) tfgm=rcpara%miss_val

!$ call omp_set_lock(omp_lp(istat))
!!!$omp critical
status=NF90_OPEN(TRIM(filename),NF90_NOWRITE,ncid) !nc-file oeffnen
!!!$omp end critical
IF (status /= NF90_NOERR) THEN
	WRITE(*,*) 'Error opening netcdf!', TRIM(filename)
	status=NF90_CLOSE(ncid)
        err3=1
!$ call omp_unset_lock(omp_lp(istat))
	RETURN
END IF
!!$ call omp_unset_lock(omp_lp(istat))
!get time dimension ID
!!!$omp critical
status=NF90_INQ_DIMID(ncid,'time',timid)
!!!$omp end critical
IF (status /= NF90_NOERR) THEN
	WRITE(*,*) 'Error inquire dimension ID', TRIM(filename), ' maybe file empty'
	status=NF90_CLOSE(ncid)
	IF (status /= NF90_NOERR) WRITE(*,*) 'Error closeing netcdf!', TRIM(filename)
!$ call omp_unset_lock(omp_lp(istat))
	RETURN
END IF
!get time dimenstion length
!!!$omp critical
status=NF90_INQUIRE_DIMENSION(ncid,timid,timedimname,maxtime)
!!!$omp end critical
IF (status /= NF90_NOERR) THEN
	WRITE(*,*) 'Error inquire maxtime ', TRIM(filename), 'maybe file empty'
	status=NF90_CLOSE(ncid)
!$ call omp_unset_lock(omp_lp(istat))
	RETURN
END IF
!get variable IDs
!!!$omp critical
status=NF90_INQ_VARID(ncid,'datum',datumvarid)
!!!$omp end critical
IF (status /= NF90_NOERR) THEN
	WRITE(*,*) 'Error inquire variable ID', TRIM(filename), 'maybe file empty'
	status=NF90_CLOSE(ncid)
!$ call omp_unset_lock(omp_lp(istat))
	RETURN
END IF
!!!$omp critical
status=NF90_INQ_VARID(ncid,names(1), tempvarid)
!!!$omp end critical
if(present(tfgm)) status=NF90_INQ_VARID(ncid,names(2), fg_depvarid)
!allocate variables with time dimension
ALLOCATE(temp(maxtime,rcpara%pmax,rcpara%parmax), datum(maxtime), fg_depar(maxtime,rcpara%pmax,rcpara%parmax))
!initialize
!temp=rcpara%miss_val
!fg_depar=rcpara%miss_val
!get variables
!!!$omp critical
status=NF90_GET_VAR(ncid,tempvarid,temp)
!!!$omp end critical
!WRITE(*,*) shape(temp), size(shape(temp))
IF (status /= NF90_NOERR) THEN
	WRITE(*,*) TRIM(NF90_STRERROR(status)), '   Error reading temperature from nc'
	status=NF90_CLOSE(ncid)
	IF (status /= NF90_NOERR) WRITE(*,*) 'Error closeing netcdf!', TRIM(filename)
!$ call omp_unset_lock(omp_lp(istat))
        DEALLOCATE(temp, datum, fg_depar)
	RETURN
END IF
if(present(tfgm)) then 
  status=NF90_GET_VAR(ncid,fg_depvarid,fg_depar)
  IF (status /= NF90_NOERR) WRITE(*,*) 'Error reading fg_departures from nc'
endif
status=NF90_GET_VAR(ncid,datumvarid, datum)
IF (status /= NF90_NOERR) WRITE(*,*) 'Error reading date from nc'
status=NF90_CLOSE(ncid)
!$ call omp_unset_lock(omp_lp(istat))

IF (status /= NF90_NOERR) WRITE(*,*) 'Error closeing netcdf!', TRIM(filename)
! temp ist (obs+biascor), wir wollen obs 

!datecount=1 !zaehlt die Tage mit Daten
!jj=datum(1)
!DO WHILE (datum(datecount) .LE. rcpara%nmax .AND. datecount .LE. maxtime)
		tm(datum(1:maxtime),:,:)=temp(1:maxtime,:,:)
!		datecount=datecount+1
!END DO

if(present(tfgm)) then
!  datecount=1 !zaehlt die Tage mit Daten
!  jj=datum(1)
!  DO WHILE (datum(datecount) .LE. rcpara%nmax .AND. datecount .LE. maxtime)
		tfgm(datum(1:maxtime),:,:)=fg_depar(1:maxtime,:,:)
!		datecount=datecount+1
!!  END DO
endif

DEALLOCATE(temp, datum, fg_depar)
err3=0

RETURN
END SUBROUTINE read_bgcorr_nc

SUBROUTINE read_sonde_corr_daily_IFS_nc(filename, rcpara,mstat,  err, ifs_rasocorrs, ifs_index, bi,ifs_rasobreaks, ifs_breakuncertainties) !
USE rfmod
USE netcdf
IMPLICIT NONE
type(rasocor_namelist), intent(in) :: rcpara
INTEGER :: ncid,  err, status, ndims !Anzahl der Dimensionen
INTEGER :: nvars, natts !Anzahl der Variablen, Attribute
INTEGER :: ii, jj, kk, i,numdat, presslen, hourlen, vals,l,bi !Anzahl der Werte (zeitlich)
INTEGER :: ind1(1), ind2(1), ind3(1)
INTEGER, ALLOCATABLE ::  dimlen(:) !laenge der Dimensionen
INTEGER, ALLOCATABLE :: dimids(:,:) !Vektor der Dimension ids fuer jede Variable
INTEGER, ALLOCATABLE ::  dimnum(:) !Anzahl der Dimensionen einer Variablen
INTEGER, ALLOCATABLE :: datum(:,:)
CHARACTER*50, ALLOCATABLE :: dimnames(:) !Namen der Dimensionen
CHARACTER*50, ALLOCATABLE :: varnames(:) !Variablennamen
INTEGER :: vartype !Variablentyp
CHARACTER* (*) filename
REAL(kind=JPRM) :: ifs_rasocorrs(rcpara%mmax,rcpara%pmax,rcpara%parmax)
REAL(kind=JPRM),ALLOCATABLE :: arr(:,:,:)
!REAL :: hilfcorr(rcpara%mmax,rcpara%pmax,rcpara%parmax)
REAL(kind=JPRM),optional :: ifs_rasobreaks(rcpara%mmax,rcpara%pmax,rcpara%parmax), ifs_breakuncertainties(rcpara%mmax,rcpara%pmax,rcpara%parmax)
INTEGER,intent(in) :: mstat
INTEGER(kind=JPRM) :: ifs_index(rcpara%mmax)
INTEGER :: datumyear, datummonth, datumday,istat
lserial=.true.
if(lserial) then
  istat=rcpara%statmax+1
else
  istat=mstat+4000
endif
!istat=rcpara%statmax+1
!istat=mstat+4000
ifs_index=0
!$ call omp_set_lock(omp_lp(istat))
!!!$omp critical
status=NF90_OPEN(TRIM(filename),NF90_NOWRITE,ncid) !nc-file oeffnen
!!!$omp end critical
IF (status /= NF90_NOERR) THEN
	WRITE(*,*) 'Error opening netcdf!', TRIM(filename)
	err=1
!$ call omp_unset_lock(omp_lp(istat))	
	RETURN
END IF
!!$ call omp_unset_lock(omp_lp(istat))	
status=NF90_Inquire(ncid,ndims,nvars,natts) !Anzahl der Dimensionen, Variablen und Attributen abfragen
IF (status /= NF90_NOERR) THEN
	WRITE(*,*) 'Error Inquire File ', TRIM(filename)
	status=NF90_CLOSE(ncid)
	err=1
!$ call omp_unset_lock(omp_lp(istat))
	RETURN
END IF
ALLOCATE (dimnames(ndims), varnames(nvars), dimlen(ndims), dimnum(nvars), dimids(nvars,NF90_MAX_VAR_DIMS))
numdat=1
DO ii=1,ndims
	status=NF90_Inquire_Dimension(ncid,ii, dimnames(ii), dimlen(ii))
	IF (status /= NF90_NOERR) WRITE(*,*) 'Error Inquire Dimension', ii, 'in File', TRIM(filename)
	IF (dimnames(ii) .EQ. 'time' ) THEN
		vals=dimlen(ii)
	ELSEIF (dimnames(ii) .EQ. 'pressure') THEN
		presslen=dimlen(ii)
	ELSEIF (dimnames(ii) .EQ. 'hour') THEN
		hourlen=dimlen(ii)
	ELSEIF (dimnames(ii) .EQ. 'numdat') THEN
		numdat=dimlen(ii)
	END IF
END DO

IF (vals .GT. rcpara%mmax .OR. presslen .NE. rcpara%pmax .OR. hourlen .NE. rcpara%parmax) THEN
	WRITE(*,*) 'Error: file ', TRIM(filename), 'has dimensions time=', vals, ' pressure=', presslen, ' hour=', hourlen, ' and needs dimensions ', rcpara%mmax, rcpara%pmax, rcpara%parmax, 'in subroutine read_sonde_corr_daily_IFS_nc in file read_txt_write_nc.f90'
END IF

ALLOCATE(datum(vals,numdat),arr(vals,rcpara%pmax,rcpara%parmax))
DO jj=1,nvars
	status=NF90_Inquire_Variable(ncid,jj, varnames(jj), vartype, dimnum(jj), dimids(jj,:))
	IF (status /= NF90_NOERR) WRITE(*,*) 'Error Inquire Variable', jj, 'in File', TRIM(filename)
	IF (varnames(jj) .EQ. 'rasocorr') THEN
		status=NF90_GET_VAR(ncid,jj,arr)
                ifs_rasocorrs=0.
                do i=1,vals 
                  ifs_rasocorrs(i,:,:)=arr(i,:,:)
                enddo
		IF (status /= NF90_NOERR) WRITE(*,*) 'Error getting Variable rasocorr from File', TRIM(filename)
	ELSEIF (varnames(jj) .EQ. 'datum') THEN
		status=NF90_GET_VAR(ncid,jj,datum(1:vals,1))
		IF (status /= NF90_NOERR) WRITE(*,*) 'Error getting Variable datum from File',TRIM(filename)
	ELSEIF (varnames(jj) .EQ. 'rasobreak' .and. present(ifs_rasobreaks)) THEN
		status=NF90_GET_VAR(ncid,jj,arr)
                ifs_rasobreaks=0.
                do i=1,vals 
                  ifs_rasobreaks(i,:,:)=arr(i,:,:)
                enddo
		IF (status /= NF90_NOERR) WRITE(*,*) 'Error getting Variable rasobreak from File', TRIM(filename)
	ELSEIF (varnames(jj) .EQ. 'breakuncertainty' .and. present(ifs_breakuncertainties)) THEN
		status=NF90_GET_VAR(ncid,jj,arr)
                ifs_breakuncertainties=0.
                do i=1,vals 
                  ifs_breakuncertainties(i,:,:)=arr(i,:,:)
                enddo
		IF (status /= NF90_NOERR) WRITE(*,*) 'Error getting Variable breaksuncertainty from File', TRIM(filename)
	END IF
END DO

ifs_index(1:vals)=datum(1:vals,1)
bi=vals

status=NF90_CLOSE(ncid)
!$ call omp_unset_lock(omp_lp(istat))

DEALLOCATE(arr,datum,dimids,dimnum,dimlen,varnames,dimnames)
err=0

END SUBROUTINE read_sonde_corr_daily_IFS_nc

SUBROUTINE read_sonde_corr_daily_nc(filename, rcpara,mstat,  err, ifs_rasocorrs, ifs_index, bi,ifs_rasobreaks, ifs_breakuncertainties) !
USE rfmod
USE netcdf
IMPLICIT NONE
type(rasocor_namelist), intent(in) :: rcpara
INTEGER :: ncid,  err, status, ndims !Anzahl der Dimensionen
INTEGER :: nvars, natts !Anzahl der Variablen, Attribute
INTEGER :: ii, jj, kk, i,numdat, presslen, hourlen, vals,l,bi,ip,ipar !Anzahl der Werte (zeitlich)
INTEGER :: ind1(1), ind2(1), ind3(1),istat
INTEGER,intent(in) :: mstat
INTEGER, ALLOCATABLE ::  dimlen(:) !laenge der Dimensionen
INTEGER, ALLOCATABLE :: dimids(:,:) !Vektor der Dimension ids fuer jede Variable
INTEGER, ALLOCATABLE ::  dimnum(:) !Anzahl der Dimensionen einer Variablen
INTEGER, ALLOCATABLE :: datum(:,:)
CHARACTER*50, ALLOCATABLE :: dimnames(:) !Namen der Dimensionen
CHARACTER*50, ALLOCATABLE :: varnames(:) !Variablennamen
INTEGER :: vartype !Variablentyp
CHARACTER* (*) filename
REAL(kind=JPRM) :: ifs_rasocorrs(rcpara%mmax,rcpara%pmax,rcpara%parmax)
REAL(kind=JPRM),ALLOCATABLE :: arr(:,:,:)
!REAL :: hilfcorr(rcpara%nmax,rcpara%pmax,rcpara%parmax)
REAL(kind=JPRM),optional :: ifs_rasobreaks(rcpara%mmax,rcpara%pmax,rcpara%parmax), ifs_breakuncertainties(rcpara%mmax,rcpara%pmax,rcpara%parmax)
INTEGER(kind=JPRM) :: ifs_index(rcpara%mmax)
INTEGER :: datumyear, datummonth, datumday
err=1
if(lserial) then
  istat=rcpara%statmax+1
else
  istat=mstat+4000
endif
!istat=rcpara%statmax+1
!istat=mstat+4000
ifs_index=0
!$ call omp_set_lock(omp_lp(istat))
!!!$omp critical
100 status=NF90_OPEN(TRIM(filename),NF90_NOWRITE,ncid) !nc-file oeffnen
!!!$omp end critical
IF (status /= NF90_NOERR) THEN
	WRITE(*,*) 'Error opening netcdf!', TRIM(filename)
!$ call omp_unset_lock(omp_lp(istat))
	RETURN
END IF
!!$ call omp_unset_lock(omp_lp(istat))
status=NF90_Inquire(ncid,ndims,nvars,natts) !Anzahl der Dimensionen, Variablen und Attributen abfragen
IF (status /= NF90_NOERR) THEN
	WRITE(*,*) 'Error Inquire File ', TRIM(filename)
	status=NF90_CLOSE(ncid)
!$ call omp_unset_lock(omp_lp(istat))
	RETURN
END IF
ALLOCATE (dimnames(ndims), varnames(nvars), dimlen(ndims), dimnum(nvars), dimids(nvars,NF90_MAX_VAR_DIMS))
numdat=1
DO ii=1,ndims
	status=NF90_Inquire_Dimension(ncid,ii, dimnames(ii), dimlen(ii))
	IF (status /= NF90_NOERR) WRITE(*,*) 'Error Inquire Dimension', ii, 'in File', TRIM(filename)
	IF (dimnames(ii) .EQ. 'time' ) THEN
		vals=dimlen(ii)
	ELSEIF (dimnames(ii) .EQ. 'pressure') THEN
		presslen=dimlen(ii)
	ELSEIF (dimnames(ii) .EQ. 'hour') THEN
		hourlen=dimlen(ii)
	ELSEIF (dimnames(ii) .EQ. 'numdat') THEN
		numdat=dimlen(ii)
	END IF
END DO

IF (vals .GT. rcpara%nmax .OR. presslen .NE. rcpara%pmax .OR. hourlen .NE. rcpara%parmax) THEN
	WRITE(*,*) 'Error: file ', TRIM(filename), 'has dimensions time=', vals, ' pressure=', presslen, ' hour=', hourlen, ' and needs dimensions ', rcpara%nmax, rcpara%pmax, rcpara%parmax, 'in subroutine read_sonde_corr_daily_IFS_nc in file read_txt_write_nc.f90'
END IF

ALLOCATE(datum(vals,numdat),arr(vals,rcpara%pmax,rcpara%parmax))
DO jj=1,nvars
	status=NF90_Inquire_Variable(ncid,jj, varnames(jj), vartype, dimnum(jj), dimids(jj,:))
	IF (status /= NF90_NOERR) WRITE(*,*) 'Error Inquire Variable', jj, 'in File', TRIM(filename)
	IF (varnames(jj) .EQ. 'rasocorr') THEN
             status=NF90_GET_VAR(ncid,jj,arr)
             ifs_rasocorrs=0.
             do i=1,vals
               ifs_rasocorrs(i,:,:)=arr(i,:,:)
             enddo
		IF (status /= NF90_NOERR) WRITE(*,*) 'Error getting Variable rasocorr from File', TRIM(filename)
	ELSEIF (varnames(jj) .EQ. 'datum') THEN
		status=NF90_GET_VAR(ncid,jj,datum(1:vals,1))
		IF (status /= NF90_NOERR) WRITE(*,*) 'Error getting Variable datum from File',TRIM(filename)
	ELSEIF (varnames(jj) .EQ. 'rasobreak' .and. present(ifs_rasobreaks)) THEN
		status=NF90_GET_VAR(ncid,jj,arr)
                do i=1,vals 
                  ifs_rasobreaks(i,:,:)=arr(i,:,:)
                enddo
		IF (status /= NF90_NOERR) WRITE(*,*) 'Error getting Variable rasobreak from File', TRIM(filename)
	ELSEIF (varnames(jj) .EQ. 'breakuncertainty' .and. present(ifs_breakuncertainties)) THEN
		status=NF90_GET_VAR(ncid,jj,arr)
                do i=1,vals 
                  ifs_breakuncertainties(i,:,:)=arr(i,:,:)
                enddo
		IF (status /= NF90_NOERR) WRITE(*,*) 'Error getting Variable breaksuncertainty from File', TRIM(filename)
	END IF
END DO

ifs_index(1:vals)=datum(1:vals,1)
ifs_index(vals+1)=rcpara%nmax
bi=vals

status=NF90_CLOSE(ncid)
!$ call omp_unset_lock(omp_lp(istat))

DEALLOCATE(arr,datum,dimids,dimnum,dimlen,varnames,dimnames)
err=0

RETURN
END SUBROUTINE read_sonde_corr_daily_nc

! Calculate monthly Brightness temperatures with RTTOV and write to NetCDF
! Leo Haimberger, 4 July 2011
!
subroutine write_sonde_monthly_bt_nc(filename,rcpara,tm,istat,err,lon,lat,ichan,skin,rasocorr,eracorr,ancorr,stname) ! schreibt Files mit Monatsmittel
USE var_def
USE rfmod
USE tosat
implicit none

type::itype
integer(kind=jpim) :: i,k
end type itype 

type(rasocor_namelist),intent(in) :: rcpara
type(rasocor_namelist) :: satpara
type(atts) attribs
type(ecskin) :: skin
integer i,k,iunit,err,l,iy,index,imon,j,m
character*(*) filename
real(kind=JPRM),intent(in) :: tm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM),intent(in),optional ::rasocorr(rcpara%nmax,rcpara%pmax,rcpara%parmax),eracorr(rcpara%nmax,rcpara%pmax,rcpara%parmax),ancorr(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: itm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: tmmon(rcpara%mmax,rcpara%pmax,rcpara%parmax),rasocorrmon(rcpara%mmax,rcpara%pmax,rcpara%parmax),eracorrmon(rcpara%mmax,rcpara%pmax,rcpara%parmax),ancorrmon(rcpara%mmax,rcpara%pmax,rcpara%parmax)
integer         :: goodmon(rcpara%mmax,rcpara%pmax,rcpara%parmax)
integer         :: imod,tgm,rgm,egm,agm,ip,ipar,iprof,thresh,sc,mkj
logical         :: lrc,ler,lan
REAL (kind=JPRM):: lon, lat,cmax,ssum
REAL(kind=JPRM), ALLOCATABLE :: hilf(:,:,:),c(:)
INTEGER :: numpar, maxatts,ind(rcpara%mmax)
TYPE(itype) hindex(rcpara%mmax*rcpara%parmax)
INTEGER, ALLOCATABLE :: datum(:,:)
REAL(kind=JPRM), ALLOCATABLE :: arr(:,:,:,:)
CHARACTER*12, ALLOCATABLE :: names(:)
CHARACTER*50, OPTIONAL :: stname
INTEGER,intent(in) :: istat

!msu
  Type(profile_Type), Allocatable   :: profiles(:)  
  character*20 :: NameOfRoutine='write_sonde_monthly_bt_nc'  
  Integer(Kind=jpim)                :: rttov_errorstatus  ! rttov error return code
  Integer(Kind=jpim)                :: nprof,nlevels
  Integer(Kind=jpim) :: alloc_status(20)
  Integer(Kind=jpim) :: asw
  Integer(Kind=jpim) :: ichan
  Integer(Kind=jpim) :: errorstatus
  Character (len=80) :: errMessage
 REAL (kind=JPRM),allocatable:: bt(:,:)

  err=0

  lrc=present(rasocorr)
  ler=present(eracorr)
  lan=present(ancorr)

  satpara=rcpara
  satpara%plevs(1:nchannels)=input_chan(1:nchannels)
  satpara%pmax=nchannels


!Parameteranzahl
numpar=1
IF(lrc) numpar=numpar+1
IF(ler) numpar=numpar+1
IF(lan) numpar=numpar+1

  thresh=2

  allocate(names(numpar), datum(rcpara%mmax,numpar))

  goodmon=0

  do ipar=1,rcpara%parmax
   do ip=1,rcpara%pmax
    index=1
    do imon=1,rcpara%mmax
      iy=rcpara%startdate/10000+(imon-1)/12
      imod=mod(imon-1,12)+1
      datum(imon,:)=toindex(100*imod+10000*iy+1,rcpara)
      tmmon(imon,ip,ipar)=0.
      rasocorrmon(imon,ip,ipar)=0.
      eracorrmon(imon,ip,ipar)=0.
      ancorrmon(imon,ip,ipar)=0.
      tgm=0
      rgm=0
      egm=0
      agm=0
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
        if(lan) then 
        if(ancorr(index,ip,ipar) .ne. rcpara%miss_val) then
          ancorrmon(imon,ip,ipar)=ancorrmon(imon,ip,ipar)+ancorr(index,ip,ipar)
          agm=agm+1
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
      if(agm .gt. 0) then
        ancorrmon(imon,ip,ipar)=ancorrmon(imon,ip,ipar)/agm
      else
        ancorrmon(imon,ip,ipar)=rcpara%miss_val
      endif   
    enddo
   enddo
  enddo

    nlevels=size(msumask,1)

    cmax=0
    allocate(c(nlevels))
    j=0
    do i=1,nlevels
      c(i)=count(goodmon(:,msumask(i,ichan),:) .gt. thresh)
      if(c(i) .gt. cmax) cmax=c(i)
      if(c(i) .lt. 0.5*cmax .and. msumask(i,ichan) .gt. 10 .and. j .eq. 0) j=i-1
    enddo
    if(j .gt. 0) nlevels=j  ! change nlevels for high terrain stations
    deallocate(c)

    j=0
    do k=1,rcpara%parmax
    do i=1,rcpara%mmax
!     if(goodmon(i,3,k) .gt. 0 .and. goodmon(i,ichan,k) .gt. 0 .and. goodmon(i,5,k) .gt. 0) then 
! upmost level should be there
     if(goodmon(i,msumask(1,ichan),k) .gt. thresh .and. skin%tsk(i,k,istat) .ne. rcpara%miss_val .and. skin%p(i,k,istat) .gt. 0) then
       m=0
       do l=1,nlevels
         mkj=msumask(l,ichan)
        if(tmmon(i,mkj,k) .ne. rcpara%miss_val) then
          m=m+1
        else
! These are the 70 and 250 hPa levels
          if(mkj .eq. 5 .or. mkj .eq. 9) then
            if(tmmon(i,mkj-1,k) .ne. rcpara%miss_val .and. &
               tmmon(i,mkj+1,k) .ne. rcpara%miss_val) then

!               select case(l)
!                 case(1) 
                   tmmon(i,mkj,k)=0.5*(tmmon(i,mkj-1,k)+tmmon(i,mkj+1,k))
!                 case(2)
!                   tmmon(i,mkj,k)=0.5*(tmmon(i,mkj-1,k)+tmmon(i,mkj+1,k))
                   if(lrc) rasocorrmon(i,mkj,k)=0.5*(rasocorrmon(i,mkj-1,k)+rasocorrmon(i,mkj+1,k))
!                 case(3)
                   if(ler) eracorrmon(i,mkj,k)=0.5*(eracorrmon(i,mkj-1,k)+eracorrmon(i,mkj+1,k))
!                 case(4)
                   if(lan) ancorrmon(i,mkj,k)=0.5*(ancorrmon(i,mkj-1,k)+ancorrmon(i,mkj+1,k))
!               end select
                       
               m=m+1
            endif
          endif
         endif
       enddo
       if(m .ge. nlevels) then
       j=j+1
           hindex(j)%i=i
         hindex(j)%k=k
!         write(*,*) 'accepted',m,nlevels,i,k
       else
!         write(*,*) 'rejected',m,nlevels,i,k
       endif
     endif
    enddo
    enddo
    nprof=j*numpar
    
    if(nprof .gt. 0) then

    allocate(bt(nprof,nchannels),stat= alloc_status(1))
    allocate(hilf(rcpara%mmax,rcpara%parmax,nchannels),stat= alloc_status(3))
    allocate( profiles(nprof),stat= alloc_status(2))

    if(any(alloc_status(1:3) .ne. 0)) THEN
      WRITE(*,*) 'write_monthly_bt allocation error'
      STOP 
    endif

    asw=1
      call rttov_alloc_prof(         &
      & errorstatus,             &
      & nprof,                   &
      & profiles,                &
      & nlevels,                 &
      & opts,                    &
      & asw,                     &
      & coefs = coefs(nrttovid), &
      & init = .true._jplm       )
      If( errorstatus /= errorstatus_success ) Then
        errorstatus = errorstatus_fatal
        Write( errMessage, '( "mem allocation error for profile arrays")' )
        Call Rttov_ErrorReport (errorstatus, errMessage, NameOfRoutine)
        Stop
      Endif

    iprof=0
    do l=1,numpar
    do j=1,rcpara%parmax
      ssum=0.
      sc=0.
      do i=1,rcpara%mmax
        if(skin%ci(i,j,istat) .ne. -999.) then
          ssum=ssum+skin%ci(i,j,istat)/10000
          sc=sc+1
        endif
      enddo
    do i=1,rcpara%mmax
!!$     if (ichan .eq. 2) write(*,*) i,j,istat,ichan,goodmon(i,msumask(1,ichan),j),tmmon(i,msumask(1:nlevels,ichan),j),skin%tsk(i,j,istat)
     if(goodmon(i,msumask(1,ichan),j) .gt. thresh .and. skin%tsk(i,j,istat) .ne. rcpara%miss_val .and. skin%p(i,j,istat) .gt. 0) then
      m=0
      do k=1,nlevels
        mkj=msumask(k,ichan)
        if(tmmon(i,mkj,j) .ne. rcpara%miss_val) then
          m=m+1
       endif

      enddo
      if(m .ge. nlevels) then 
      iprof=iprof+1
!      iprof=(l-1)*rcpara%mmax*rcpara%parmax+(j-1)*rcpara%mmax+i
      m=0
      do k=1,nlevels
        if(tmmon(i,msumask(k,ichan),j) .ne. rcpara%miss_val) then
          m=m+1
          profiles(iprof)%p(m)=rcpara%plevs(msumask(k,ichan))
          if(l .eq. 1) then 
            profiles(iprof)%t(m)=tmmon(i,msumask(k,ichan),j)
!            if(m.eq.1) then
!              write(*,*) 'iprof',l,i,iprof,profiles(iprof)%t(m)
!            endif
          endif
          if(l .eq. 2) then
            profiles(iprof)%t(m)=tmmon(i,msumask(k,ichan),j)-rasocorrmon(i,msumask(k,ichan),j)
!            if(m.eq.1) then
!              write(*,*) 'iprof',l,i,iprof,profiles(iprof)%t(m)
!            endif
          endif
          if(l .eq. 3) then
            profiles(iprof)%t(m)=eracorrmon(i,msumask(k,ichan),j)
!            if(m.eq.1) then
!              write(*,*) 'iprof',l,i,iprof,profiles(iprof)%t(m)
!            endif
          endif
          if(l .eq. 4) then
            profiles(iprof)%t(m)=ancorrmon(i,msumask(k,ichan),j)
!            if(m.eq.1) then
!              write(*,*) 'iprof',l,i,iprof,profiles(iprof)%t(m)
!            endif
          endif
        endif
!      if(iprof .eq. 1) then 
!      else
!        profiles(iprof)%q=profiles(iprof-1)%q
!      endif
      enddo

        profiles(iprof)%q(1:m)=1.!6.1121*exp((18.678-(profiles(iprof)%t(1:nlevels)-273.15)/234.5)*&
                             !      (profiles(iprof)%t(1:nlevels)-273.15)/(257.14+(profiles(iprof)%t(1:nlevels)-273.15)))/profiles(iprof)%p(1:nlevels)*1.e4

      profiles(iprof) % s2m % t= skin%t2(i,j,istat)
      profiles(iprof) % skin % t =skin%tsk(i,j,istat)
      if(profiles(iprof)%q(m) .ne. rcpara%miss_val) then
        profiles(iprof) % s2m % q= profiles(iprof)%q(m)
      else
        profiles(iprof) % s2m % q= 10000.
      endif
!      if(skin%p(i,j,istat) .lt. 0) then
!      write(*,*) i,j,istat,skin%p(i,j,istat),'false p'
!      endif
!      if(istat .eq. 2304) then
!      write(*,*) i,j,istat,skin%p(i,j,istat)
!      endif
      profiles(iprof) % s2m % p= skin%p(i,j,istat)/100. ! hPa
      profiles(iprof) % s2m % u= 0. 
      profiles(iprof) % s2m % v= 0. 
      profiles(iprof) % skin % fastem =(/3.0, 5.0, 15.0, 0.1, 0.3/)
      profiles(iprof) % skin % surftype=1 ! sea
      if(skin%lsm(istat) .gt. 0.5) then
        profiles(iprof) % skin % surftype=0 !land
      endif

! shift from E40 to EI caused a hard to find change here -- do not
! use time varying surface type
!      if(skin%ci(i,j,istat) .ne. -999. .and. skin%ci(i,j,istat)/10000 .gt. 0.5) then
!        profiles(iprof) % skin % surftype=2 !sea ice
!      endif   
! instead use a mean value averaged over time period
      if(sc .gt. 0 .and. ssum/sc .gt. 0.5) then
         profiles(iprof) % skin % surftype=2 !sea ice
      endif

      profiles(iprof) % skin % watertype=1 ! ocean water 
      profiles(iprof) % skin % fastem =(/3.0, 5.0, 15.0, 0.1, 0.3/) ! default values
      profiles(iprof) % elevation = skin%elev(istat)/1000 !km
      profiles(iprof) % latitude = lat
!!$    if(ichan .eq. 2) then
!!$    write(*,*) iprof,i,j,istat
!!$    write(*,*) profiles(iprof)%t,profiles(iprof)%s2m,profiles(iprof)%skin,profiles(iprof)%elevation,profiles(iprof)%latitude
!!$    endif
      endif
      endif
    enddo
    enddo
    enddo

    do i=1,nprof
      if(profiles(i)%skin%fastem(1) .eq. 0) then
        write(*,*) i,profiles(i)%skin%fastem
      endif
    enddo
!!$ call omp_set_lock(omp_lp(rcpara%statmax+1))
    call msu_fwd_calc(profiles,bt,nprof,nlevels)
!!$ call omp_unset_lock(omp_lp(rcpara%statmax+1))

!!$if(ichan .eq. 2) then
!!$  do i=1,nprof
!!$    write(*,*) i,bt(i,:)
!!$  enddo
!!$endif
    deallocate(profiles)
  where(isnan(bt)) 
    bt=rcpara%miss_val
  endwhere
  where(bt .gt. 1.e30) 
    bt=rcpara%miss_val
  endwhere


  hilf=rcpara%miss_val
  do i=1,nprof/numpar
    hilf(hindex(i)%i,hindex(i)%k,:)=bt(i,:)
  enddo
  

j=0 
 do i=1,rcpara%mmax
   if(any(hilf(i,:,:) .ne. rcpara%miss_val)) then
     j=j+1
     datum(j,:)=datum(i,:)
     ind(j)=i
   endif
 enddo

 if(j .eq. 0) then

  err=1

 else

  maxatts=5
  
  ALLOCATE (attribs%attname(maxatts), attribs%long_name(numpar), attribs%units(numpar), attribs%missing_val(numpar), attribs%valid_range(numpar,2),attribs%cell_methods(numpar))
  
  attribs%attname=(/'long_name','units','missing_value','valid_range','cell_methods'/)
  attribs%units='K'
  attribs%missing_val=rcpara%miss_val
  attribs%valid_range=rcpara%miss_val
!  datum=rcpara%miss_val
  attribs%valid_range(1,:)=(/0, 400/)
  IF (numpar .GE. 2) THEN
  	attribs%valid_range(2,:)=(/-20, 20/)
  END IF
  IF (numpar .GE. 2) THEN
  	attribs%valid_range(2,:)=(/-20, 20/)
        attribs%cell_methods(2)=''
  END IF
  attribs%cell_methods='time: mean over months'
 
  
 ALLOCATE(arr(j,satpara%pmax,satpara%parmax,numpar))
  do i=1,satpara%pmax
  do k=1,satpara%parmax
     arr(:,i,k,1)=hilf(ind(1:j),k,i)
   enddo
   enddo
  names(1)='montemp'
  attribs%long_name(1)='monthly_uncorrected_bt'
  
  IF (numpar .ge. 2 ) THEN
        
  do i=1,nprof/numpar
    hilf(hindex(i)%i,hindex(i)%k,:)=bt(i+nprof/numpar,:)
  enddo
  do k=1,satpara%parmax
   do i=1,satpara%pmax
     arr(:,i,k,2)=hilf(ind(1:j),k,i)
   enddo
  enddo
	names(2)='rasocorrmon'
	attribs%long_name(2)='monthly_corrected_raso_bt'
  END IF
  IF (numpar.ge.3) THEN
  do i=1,nprof/numpar
    hilf(hindex(i)%i,hindex(i)%k,:)=bt(i+2*nprof/numpar,:)
  enddo
  do k=1,satpara%parmax
   do i=1,satpara%pmax
     arr(:,i,k,3)=hilf(ind(1:j),k,i)
   enddo
  enddo
	names(3)='eracorrmon'
	attribs%long_name(3)='monthly_corrected_era_bt'
  END IF
  IF (numpar.eq.4) THEN
  do i=1,nprof/numpar
    hilf(hindex(i)%i,hindex(i)%k,:)=bt(i+3*nprof/numpar,:)
  enddo
  do k=1,satpara%parmax
   do i=1,satpara%pmax
     arr(:,i,k,4)=hilf(ind(1:j),k,i)
   enddo
  enddo
	names(4)='ancorrmon'
	attribs%long_name(3)='monthly_corrected_eraan_bt'
  END IF
  IF (filename(9:30) .EQ. 'feedbackglobbincorrmon' ) THEN
	  attribs%title='monthly radiosonde brightness temperatures, corrections, backgroundcorrection'
	  attribs%long_name(1)='monthly_temperatures'
  ELSEIF (filename(9:25) .EQ. 'feedbackglobbgmon') THEN
  	attribs%title='corrected background monthly'
	attribs%long_name(1)='monthly_corrected_background_temperatures'
  ELSEIF (filename(9:30) .EQ. 'feedbackglobbincompmon') THEN
  	attribs%title='monthly radiosonde temperatures, corrections, backgroundcorrection with composits'
  END IF


!!$omp critical
  CALL gen_nc(satpara,istat,filename,numpar,maxatts,j,arr, names, attribs, lat, lon, datum(1:j,:), numpar,stname=stname) !in this file line 276
!!$omp end critical
  DEALLOCATE(arr)
  DEALLOCATE (attribs%attname, attribs%long_name, attribs%units, attribs%missing_val, attribs%valid_range,attribs%cell_methods)
 endif
    
  DEALLOCATE (bt,hilf)
  else
    err=1
  endif

  DEALLOCATE(names,datum)

  return

end subroutine write_sonde_monthly_bt_nc

subroutine write_eclim_monthly_nc(filename,rcpara,teclimgens,err,talon,talat,ichan,gskin) ! schreibt Files mit Monatsmittel
USE var_def
USE rfmod
use tosat
implicit none

type(rasocor_namelist),intent(in) :: rcpara
type(ecskin) :: gskin
type(atts) attribs
integer i,iunit,err,l,iy,index,imon,j,k,m,istat,ni,nj,ii,jj
character*(*) filename
real(kind=JPRM),intent(in) :: teclimgens(:,:,:,:,:)
integer         :: ip,ipar,iprof,mkj,varc
REAL (kind=JPRM):: talon(:), talat(:),ssum,sc
INTEGER :: numpar, maxatts,ind(rcpara%mmax),tsize(5)
INTEGER, ALLOCATABLE :: datum(:),idxarr(:,:)
REAL(kind=JPRM), ALLOCATABLE :: arr(:,:,:,:,:),hilf(:,:,:)
CHARACTER*12, ALLOCATABLE :: names(:)

!msu
  Type(profile_Type), Allocatable   :: profiles(:)  
  character*20 :: NameOfRoutine='write_sonde_monthly_bt_nc'  
  Integer(Kind=jpim)                :: rttov_errorstatus  ! rttov error return code
  Integer(Kind=jpim)                :: nprof,nlevels
  Integer(Kind=jpim) :: alloc_status(20)
  Integer(Kind=jpim) :: asw
  Integer(Kind=jpim) :: ichan
  Integer(Kind=jpim) :: errorstatus
  Character (len=80) :: errMessage
 REAL (kind=JPRM),allocatable:: bt(:,:)


  err=0


tsize=shape(teclimgens)
numpar=1
ni=tsize(1)
nj=tsize(2)

allocate(names(numpar), datum(rcpara%mmax))
  
    nlevels=size(msumask,1)

    nprof=rcpara%mmax*ni*numpar
    


 ALLOCATE(arr(ni,nj,rcpara%mmax,numpar,1))
 arr=rcpara%miss_val
    allocate(idxarr(4,nprof))

if(.true.) then
    do l=1,numpar
    do jj=1,nj


    allocate(profiles(nprof),stat= alloc_status(1))

    if((alloc_status(1) .ne. 0)) THEN
      WRITE(*,*) 'write_monthly_bt allocation error'
      STOP 
    endif

    asw=1
      call rttov_alloc_prof(         &
      & errorstatus,             &
      & nprof,                   &
      & profiles,                &
      & nlevels,                 &
      & opts,                    &
      & asw,                     &
      & coefs = coefs(nrttovid), &
      & init = .true._jplm       )
      If( errorstatus /= errorstatus_success ) Then
        errorstatus = errorstatus_fatal
        Write( errMessage, '( "mem allocation error for profile arrays")' )
        Call Rttov_ErrorReport (errorstatus, errMessage, NameOfRoutine)
        Stop
      Endif

    iprof=0
    do ii=1,ni
    istat=(jj-1)*ni+ii
    ssum=0.
    sc=0
    do k=1,rcpara%mmax
    do j=1,2
      if(gskin%ci(k,j,istat) .ne. rcpara%miss_val .and. gskin%ci(k,j,istat) .ne. 9999) then
        ssum=ssum+gskin%ci(k,j,istat)
        sc=sc+1
      endif
    enddo
    enddo
    do i=1,rcpara%mmax
!      if(gskin%p(i,1,istat) .ne. rcpara%miss_val .and. teclimgens(i,ii,jj,1,l) .ne. rcpara%miss_val) then
      if(gskin%p(i,1,istat) .ne. rcpara%miss_val .and. teclimgens(ii,jj,1,i,l) .ne. rcpara%miss_val) then
      iprof=iprof+1
      idxarr(:,iprof)=(/ii,jj,i,l/)
      do k=1,nlevels
        profiles(iprof)%p(k)=rcpara%plevs(msumask(k,ichan))
!        profiles(iprof)%t(k)=teclimgens(i,ii,jj,msumask(k,ichan),l)
        profiles(iprof)%t(k)=teclimgens(ii,jj,msumask(k,ichan),i,l)
      enddo

        profiles(iprof)%q(1:nlevels)=1.!6.1121*exp((18.678-(profiles(iprof)%t(1:nlevels)-273.15)/234.5)*&
                             !      (profiles(iprof)%t(1:nlevels)-273.15)/(257.14+(profiles(iprof)%t(1:nlevels)-273.15)))/profiles(iprof)%p(1:nlevels)*1.e4

      profiles(iprof) % s2m % t= sum(gskin%t2(i,1:2,istat))/2
      profiles(iprof) % skin % t =sum(gskin%tsk(i,1:2,istat))/2
      if(profiles(iprof)%q(nlevels) .ne. rcpara%miss_val) then
        profiles(iprof) % s2m % q= profiles(iprof)%q(nlevels)
      else
        profiles(iprof) % s2m % q= 10000.
      endif
      profiles(iprof) % s2m % p= sum(gskin%p(i,1:2,istat))/100./2 ! hPa
      profiles(iprof) % s2m % u= 0. 
      profiles(iprof) % s2m % v= 0. 
      profiles(iprof) % skin % fastem =(/3.0, 5.0, 15.0, 0.1, 0.3/)
      profiles(iprof) % skin % surftype=1 ! sea
      if(gskin%lsm(istat) .gt. 0.5) then
        profiles(iprof) % skin % surftype=0 !land
      endif

! shift from E40 to EI caused a hard to find change here -- do not
! use time varying surface type
!      if(skin%ci(i,j,istat) .ne. -999. .and. skin%ci(i,j,istat)/10000 .gt. 0.5) then
!        profiles(iprof) % skin % surftype=2 !sea ice
!      endif   
! instead use a mean value averaged over time period
      if(sc .gt. 0 .and. ssum/sc .gt. 0.5) then
         profiles(iprof) % skin % surftype=2 !sea ice
      endif

      profiles(iprof) % skin % watertype=1 ! ocean water 
      profiles(iprof) % skin % fastem =(/3.0, 5.0, 15.0, 0.1, 0.3/) ! default values
      profiles(iprof) % elevation = gskin%elev(istat)/1000 !km
      profiles(iprof) % latitude = talat(jj)
      profiles(iprof) % longitude = talon(ii)
!    if(ichan .eq. 4) then
!    write(*,*) iprof,i,j,istat
!    write(*,*) 't',profiles(iprof)%t,'s2m',profiles(iprof)%s2m,'skin',profiles(iprof)%skin,'elevation',profiles(iprof)%elevation,'lat',profiles(iprof)%latitude

     endif
    enddo
    enddo


    allocate(bt(iprof,nchannels),stat= alloc_status(1))

    call msu_fwd_calc(profiles(1:iprof),bt,iprof,nlevels)
    deallocate(profiles)
  where(isnan(bt)) 
    bt=rcpara%miss_val
  endwhere
  where(bt .gt. 1.e30) 
    bt=rcpara%miss_val
  endwhere

do k=1,iprof
  arr(idxarr(1,k),idxarr(2,k),idxarr(3,k),idxarr(4,k),1)=bt(k,ichan-1)
enddo
    enddo
    enddo
endif


  maxatts=5
  
  ALLOCATE (attribs%attname(maxatts), attribs%long_name(numpar), attribs%units(numpar), attribs%missing_val(numpar), attribs%valid_range(numpar,2),attribs%cell_methods(numpar))
  
  attribs%attname=(/'long_name','units','missing_value','valid_range','cell_methods'/)
  attribs%units='K'
  attribs%missing_val=rcpara%miss_val
  attribs%valid_range=rcpara%miss_val
!  datum=rcpara%miss_val
  attribs%valid_range(1,:)=(/0, 400/)
  attribs%cell_methods='time: mean over months'
 
!!$  l=0
!!$  do ipar=1,numpar
!!$  do ij=1,nj
!!$  do ii=1,ni
!!$  do im=1,rcpara%mmax
!!$    l=l+1
!!$    arr(idxarr(,ii,ij,ipar)=bt(l,ichan)
!!$  enddo
!!$  enddo
!!$  enddo
!!$  enddo

  names(1)='montemp'
  attribs%title='monthly ERA-20CM temperatures'
  attribs%long_name(1)='monthly_temperatures'
  if(any(abs(arr) .gt. 1.e10)) then
    write(*,*) 'invalid value'
    call abort
  endif


  call gen_nc_3D(filename,(/'months','longitude','latitude'/),(/'months      ','longitude','latitude','brightness_temperature'/), attribs, talon, talat, datum, arr) !
    
  DEALLOCATE (attribs%attname, attribs%long_name, attribs%units, attribs%missing_val, attribs%valid_range,attribs%cell_methods)
  DEALLOCATE(arr,names,datum)

  return

end subroutine write_eclim_monthly_nc

subroutine write_sonde_monthly_nc(filename,rcpara,tm,istat,err,lon,lat,rasocorr,eracorr,stname) ! schreibt Files mit Monatsmittel
USE var_def
USE rfmod
implicit none

type(rasocor_namelist),intent(in) :: rcpara
type(atts) attribs
integer i,iunit,err,l,iy,index,imon,j
character*(*) filename
real(kind=JPRM),intent(in) :: tm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM),intent(in),optional ::rasocorr(rcpara%nmax,rcpara%pmax,rcpara%parmax),eracorr(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: itm(rcpara%nmax,rcpara%pmax,rcpara%parmax)
real(kind=JPRM) :: tmmon(rcpara%mmax,rcpara%pmax,rcpara%parmax),rasocorrmon(rcpara%mmax,rcpara%pmax,rcpara%parmax),eracorrmon(rcpara%mmax,rcpara%pmax,rcpara%parmax)
integer         :: goodmon(rcpara%mmax,rcpara%pmax,rcpara%parmax),imod,tgm,rgm,egm,ip,ipar
logical         :: lrc,ler
REAL (kind=JPRM):: lon, lat
INTEGER :: numpar, maxatts,ind(rcpara%mmax)
INTEGER, ALLOCATABLE :: datum(:,:)
REAL(kind=JPRM), ALLOCATABLE :: arr(:,:,:,:)
CHARACTER*12, ALLOCATABLE :: names(:)
CHARACTER*50,OPTIONAL :: stname
INTEGER,intent(in) :: istat
  err=0

  lrc=present(rasocorr)
  ler=present(eracorr)
!Parameteranzahl
numpar=2
IF(lrc) numpar=numpar+1
IF(ler) numpar=numpar+1

allocate(names(numpar), datum(rcpara%mmax,numpar))
  goodmon=0

  do ipar=1,rcpara%parmax
   do ip=1,rcpara%pmax
    index=1
    do imon=1,rcpara%mmax
      iy=rcpara%startdate/10000+(imon-1)/12
      imod=mod(imon-1,12)+1
      datum(imon,:)=toindex(100*imod+10000*iy+1,rcpara)
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

  
  maxatts=5
  
  ALLOCATE (attribs%attname(maxatts), attribs%long_name(numpar), attribs%units(numpar), attribs%missing_val(numpar), attribs%valid_range(numpar,2),attribs%cell_methods(numpar))
  
  attribs%attname=(/'long_name','units','missing_value','valid_range','cell_methods'/)
  attribs%units='K'
  attribs%units(2)=''
  attribs%missing_val=rcpara%miss_val
  attribs%valid_range=rcpara%miss_val
!  datum=rcpara%miss_val
  attribs%valid_range(1,:)=(/0, 400/)
  attribs%valid_range(2,:)=(/0, 31/)
  IF (numpar .GE. 3) THEN
  	attribs%valid_range(3,:)=(/-20, 20/)
  END IF
  IF (numpar .GE. 4) THEN
  	attribs%valid_range(4,:)=(/-20, 20/)
  END IF
  attribs%cell_methods='time: mean over months'
  attribs%cell_methods(2)=''
 
j=0 
 do i=1,rcpara%mmax
   if(any(goodmon(i,:,:) .ne. 0)) then
     j=j+1
     datum(j,:)=datum(i,:)
     ind(j)=i
   endif
 enddo
 if(j .eq. 0) then
  err=1
  return
 endif
 ALLOCATE(arr(j,rcpara%pmax,rcpara%parmax,numpar))
 arr(:,:,:,1)=tmmon(ind(1:j),:,:)
  names(1)='montemp'
  arr(:,:,:,2)=goodmon(ind(1:j),:,:)
  names(2)='goodmon'
  attribs%long_name(2)='number_of_values_in_month'
  
  IF (numpar==3 .AND. ler) THEN
  	arr(:,:,:,3)=eracorrmon(ind(1:j),:,:)
	names(3)='eracorrmon'
	attribs%long_name(3)='monthly_era_correction'
  ELSEIF (numpar==3 .AND. lrc .OR. numpar==4) THEN
  	arr(:,:,:,3)=rasocorrmon(ind(1:j),:,:)
	names(3)='rasocorrmon'
	attribs%long_name(3)='monthly_raso_correction'
  END IF
  IF (numpar==4) THEN
  	arr(:,:,:,4)=eracorrmon(ind(1:j),:,:)
	names(4)='eracorrmon'
	attribs%long_name(4)='monthly_era_correction'
  END IF
  IF (filename(9:30) .EQ. 'feedbackglobbincorrmon' ) THEN
	  attribs%title='monthly radiosonde temperatures, corrections, backgroundcorrection'
	  attribs%long_name(1)='monthly_temperatures'
  ELSEIF (filename(9:25) .EQ. 'feedbackglobbgmon') THEN
  	attribs%title='corrected background monthly'
	attribs%long_name(1)='monthly_corrected_background_temperatures'
  ELSEIF (filename(9:30) .EQ. 'feedbackglobbincompmon') THEN
  	attribs%title='monthly radiosonde temperatures, corrections, backgroundcorrection with composits'
  END IF
  if(any(abs(arr) .gt. 1.e10)) then
    write(*,*) 'invalid value'
    return
    call abort
  endif
!!$omp critical
  CALL gen_nc(rcpara,istat,filename,numpar,maxatts,j,arr, names, attribs, lat, lon, datum(1:j,:), numpar,stname=stname) !in this file line 276
!!$omp end critical
    
  DEALLOCATE (attribs%attname, attribs%long_name, attribs%units, attribs%missing_val, attribs%valid_range,attribs%cell_methods)
  DEALLOCATE(arr,names,datum)

  return

end subroutine write_sonde_monthly_nc


INTEGER FUNCTION read_20CR_PL_LH(current_year,current_varno,cr20,rcpara)

  use netcdf
  use rfmod

  !EXTERNAL
  type(rasocor_namelist),intent(in) :: rcpara
  INTEGER :: current_year
  INTEGER :: current_varno
  INTEGER :: time_dim !time dimension = numbers of items
  REAL(KIND=JPRM) :: cr20(rcpara%ni,rcpara%nj,rcpara%pmax,rcpara%parmax,1)
  !INTERNAL
  TYPE(twenty_CR):: tw_CR_station
  INTEGER :: status,status1
  INTEGER :: ind
  CHARACTER (len = 100):: path_name
  CHARACTER (len = 50 ):: file_name
  CHARACTER (len = 10 ):: string_var
  CHARACTER (len = 10 ):: string_year
  

  INTEGER :: ncid !!NetCDF file Id number
                !!number of dimensions
                         !!number of variables
                                   !!Number of global attributes
    INTEGER :: N_dim, N_var, N_gatt
    CHARACTER (len = 100):: title
    INTEGER :: lat_index_varid, lon_index_varid
    CHARACTER (len = 50) :: lon_name,lat_name,lev_name,pres_l_name,time_name,air_name
    INTEGER :: lon_varid,lat_varid,lev_varid,time_varid,air_varid
    INTEGER :: lon_dim,lat_dim,lev_dim
    INTEGER :: lat,lon,lev,pres_l_varid,pres_l_dim,air_dim
    REAL :: longitude(180)               !longitude
    REAL :: latitude(91)                 !atitude
    REAL :: p_levels(24)                  !pressure levels
    REAL*8 , ALLOCATABLE,DIMENSION(:) :: time_index!!real double precision
    !REAL , ALLOCATABLE :: temp_arr(:,:,:,:)
    INTEGER*2, ALLOCATABLE, DIMENSION(:):: short_data_PL !data from the 20CR as short integer
    INTEGER :: i,j,k,l,z

    
    !now I have to create the file name that I have to use for the path_name
    if(current_varno .eq. 2)then
       string_var= 'air'
    elseif(current_varno .eq. 3)then
       string_var= 'uwnd'
    elseif(current_varno .eq. 4)then
       string_var= 'vwnd'
    elseif(current_varno .eq. 0)then
       string_var= 'hgt'
    else
           write(*,9000)trim(file_name)
9000 FORMAT(/, "      ---->INCORRECT FILE CAPTION<---- ", A15)
       STOP
    end if
    write(string_year,'(I4)')current_year
    !now I have the file name
    file_name = trim(string_var)//"."//trim(string_year)//".nc"
    path_name = "/scratch/srvx7/20CR/"//trim(file_name)
    !new path /scratch/srvx7/20CR
   !OLD  path_name = "/home/srvx9/scratch/20CR/"//trim(file_name)
    !now I need to extract the file name
    !IND=1
    !check: do while(IND > 0)
    !   IND=index(path_name,'.nc')
    !   IF(IND >0) THEN
    !      file_name=path_name(IND-8:IND+2)
    !      IND = -1
    !   ELSE
    !!      WRITE(*,*)'THIS FILES IS NOT CORRECT'
    !      WRITE(*,*)'PROGRAM TERMINATED'
    !      STOP
    !   ENDIF
    !end do check
    write(*,900)trim(file_name)
900 FORMAT(/, "      ----> I m working on File-> ", A15)
    
   !Here I set the VARNO for the 20CR station
    tw_CR_station%varno = current_varno

    !OPEN NETCDF files
    status = nf90_open(path_name, nf90_nowrite,ncid)
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error opening file", path_name
       STOP
    end if
    
    !inquire file:->
    !dimension Number, dimension Variables, dimension global setting
    status = nf90_inquire(ncid,N_dim, N_var, N_gatt )
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error Inquire file: ",  file_name
       STOP
    end if
    
    !write(*,*)N_dimsin, N_varsin, N_gatts_in
    
    !!NOW I READ THE GLOBAL ATTRIBUTE 
    status =  nf90_get_att(ncid,nf90_global,"title",title)
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error writing Global Attribute ->","title"
       STOP
    end if
    WRITE(*,*)"title-> ",trim(title)
    
    
    
    !LONGITUDE
    status = nf90_inq_varid(ncid,"lon",lon_varid)
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error inquire Longitude variable id"
       STOP
    end if
    
    !inquire dimension
    !status=NF90_INQUIRE_DIMENSION(ncid,1,lon_name,lon_dim)
    lon_varid=1
    status=NF90_INQUIRE_DIMENSION(ncid,lon_varid,lon_name,lon_dim)
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error inquire Longitude  dimension"
       STOP
    end if
    ALLOCATE (tw_CR_station%lon(lon_dim), STAT = status)
    if(status .eq. 0)then
       tw_CR_station%lon=0
    else
       WRITE(*,*)"Allocate memory ERROR, variable ""longitude"""
       STOP
    end if
    !longitude=0
    !get the variable
    !status = nf90_get_var(ncid,3,tw_CR_station%lon)
    lon_varid=3
    status = nf90_get_var(ncid,lon_varid,tw_CR_station%lon)
    !status = nf90_get_var(ncid,3,longitude)
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error get Longitude variable"
       STOP
    end if
    
    !LATITUDE
    status = nf90_inq_varid(ncid,"lat",lat_varid)
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error inquire Latitude variable id"
       STOP
    end if
    !inquire dimension
    !status=NF90_INQUIRE_DIMENSION(ncid,2,lat_name,lat_dim)
    lat_varid=2
    status=NF90_INQUIRE_DIMENSION(ncid,2,lat_name,lat_dim)
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error inquire Latutude dimension"
       STOP
    end if
    
    ALLOCATE (tw_CR_station%lat(lat_dim), STAT = status)
    if(status .eq. 0)then
       tw_CR_station%lat=0
    else
       WRITE(*,*)"Allocate memory ERROR, variable ""latitude"""
       STOP
    end if
    
    !get the variable
    !status = nf90_get_var(ncid,2,tw_CR_station%lat)
    status = nf90_get_var(ncid,lat_varid,tw_CR_station%lat)
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error get Latitude variable"
       STOP
    end if
    
    !PRESSURE LAYERS
    status = nf90_inq_varid(ncid,"level",pres_l_varid)
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error inquire Pressure Layers variable id"
       STOP
    end if
    !inquire dimension
    !status=NF90_INQUIRE_DIMENSION(ncid,3,pres_l_name,pres_l_dim)
    pres_l_varid = 3!there is a NetCDF file incongruence
    status=NF90_INQUIRE_DIMENSION(ncid,pres_l_varid,pres_l_name,pres_l_dim)
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error inquire Pressure Layers  dimension"
       STOP
    end if

 ALLOCATE (tw_CR_station%pressure_layers(pres_l_dim), STAT = status)
    if(status .eq. 0)then
       tw_CR_station%pressure_layers=0
    else
       WRITE(*,*)"Allocate memory ERROR, variable ""pressure_layers"""
       STOP
    end if
  !get the variable
    !status = nf90_get_var(ncid,1,tw_CR_station%pressure_layers)
     pres_l_varid = 1
    status = nf90_get_var(ncid,pres_l_varid,tw_CR_station%pressure_layers)
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error get Pressure Layers variable"
       STOP
    end if
    !TIME
  status = nf90_inq_varid(ncid,"time",time_varid)
  if(status .NE. NF90_NOERR) then
     PRINT*, "Error inquire Time variable id"
     STOP
  end if
  !inquire dimension
  !status=NF90_INQUIRE_DIMENSION(ncid,4,time_name,time_dim)
  time_varid =4
  status=NF90_INQUIRE_DIMENSION(ncid,time_varid,time_name,time_dim)
  if(status .NE. NF90_NOERR) then
     PRINT*, "Error inquire Time  dimension"
     STOP
  end if

  !allocate the memory
  ALLOCATE (time_index(time_dim), STAT = status)
  if(status .NE. NF90_NOERR) then
     PRINT*, "Allocation error ->",TRIM(time_name)
     STOP
  end if
  time_index=0
  !get the variable
  !status = nf90_get_var(ncid,4,time_index)
  status = nf90_get_var(ncid,time_varid,time_index)
  if(status .NE. NF90_NOERR) then
     PRINT*, "Error get Time variable"
     STOP
  end if
  !Now I need to convert this index in data and time
  !the correct data and time is stored in the tw_CR_station type
 ! status = index_generate(file_name,time_dim,time_index(1),time_index(time_dim),tw_CR_station)

  !VARIABLE string_var= 'air','uwnd','vwnd'
  status = nf90_inq_varid(ncid,string_var,air_varid)
  if(status .NE. NF90_NOERR) then
     PRINT*, "Error inquire", string_var ," variable id"
     STOP
  end if

  !allocate the memory
  ALLOCATE (tw_CR_station%short_data_PL(lon_dim,lat_dim,pres_l_dim,time_dim), STAT = status)
   if(status .NE. NF90_NOERR) then
     PRINT*, "Allocation error ->",TRIM(air_name)
     STOP
  end if
  tw_CR_station%short_data_PL=0
  
  !get the variable
  status = nf90_get_var(ncid,5,tw_CR_station%short_data_PL)
  if(status .NE. NF90_NOERR) then
     PRINT*, NF90_STRERROR(status)
     STOP
  end if  
  
  !Close file
      status=NF90_CLOSE(ncid)
      if(status .NE. NF90_NOERR) then
         PRINT*, "Closing Error, with file: ",file_name
         STOP
      end if
      
      write(*,800)trim(file_name)
800   FORMAT(/, "       --> Reading NetCDF Files: ",A15," ----> OK")

  do l=1,365
  do m=1,2
!$OMP parallel do private(k,z,i,j)
  do k=1,rcpara%pmax
    do z=1,pres_l_dim
      if(tw_CR_station%pressure_layers(z) .eq. rcpara%plevs(k)) exit
    enddo
    if(z .ne. pres_l_dim+1) then
    do j=1,lat_dim
      do i=1,lon_dim
        cr20(2*i-1,2*j-1,k,m,1)=cr20(2*i-1,2*j-1,k,m,1)+0.01*tw_CR_station%short_data_PL(i,lat_dim+1-j,z,(l-1)*4+2*(m-1)+1)+427.66
        cr20(2*i,2*j-1,k,m,1)=cr20(2*i,2*j-1,k,m,1)+0.01*tw_CR_station%short_data_PL(i,lat_dim+1-j,z,(l-1)*4+2*(m-1)+1)+427.66
        if(j .lt. lat_dim) then
        cr20(2*i-1,2*j,k,m,1)=cr20(2*i-1,2*j,k,m,1)+0.01*tw_CR_station%short_data_PL(i,lat_dim+1-j,z,(l-1)*4+2*(m-1)+1)+427.66
        cr20(2*i,2*j,k,m,1)=cr20(2*i,2*j,k,m,1)+0.01*tw_CR_station%short_data_PL(i,lat_dim+1-j,z,(l-1)*4+2*(m-1)+1)+427.66
        endif
      enddo
    enddo
    endif
  enddo
!$OMP end parallel do
  enddo
  enddo
  cr20=cr20/365

  deallocate(tw_CR_station%short_data_PL,tw_CR_station%pressure_layers,tw_CR_station%lon,tw_CR_station%lat)

read_20CR_PL_LH = 0
END FUNCTION read_20CR_PL_LH

subroutine read_hadCRUT4(filename,crut4,rcpara)

use netcdf
use rfmod

implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer iunit,err,jahr,monat,i,j,k,status,ncid,timid,maxtime,timevarid,tempvarid
integer crustart,cruoffset,crut4stop,crustop,dimids(3)
character*10 :: timedimname,syear
character*100:: filename

real(kind=JPRM) :: crut4(72,36,rcpara%mmax)
real(kind=JPRM), ALLOCATABLE :: tanomaly(:,:,:)

status=NF90_OPEN(TRIM(filename),NF90_NOWRITE,ncid) !nc-file oeffnen
status=NF90_INQ_DIMID(ncid,'time',timid)
timedimname='time'
status=NF90_INQUIRE_DIMENSION(ncid,timid,timedimname,maxtime)
status=NF90_INQ_VARID(ncid,'time', timevarid)
if(status .ne. 0) then
  write(*,*) nf90_strerror(status)
endif
status=NF90_INQ_VARID(ncid,'temperature_anomaly', tempvarid)
status=NF90_INQUIRE_VARIABLE(ncid,tempvarid,dimids=dimids)
status=NF90_GET_ATT(ncid,timevarid,'start_year',crustart)
!read(syear,*) crustart
cruoffset=(rcpara%startdate/10000-crustart)*12
if (maxtime<cruoffset+rcpara%mmax) then
  crustop=maxtime
  crut4stop=maxtime-cruoffset
else
  crustop=cruoffset+rcpara%mmax
  crut4stop=rcpara%mmax
endif
allocate(tanomaly(72,36,maxtime))
status=NF90_GET_VAR(ncid,tempvarid, tanomaly)
where (tanomaly.eq. -1.e30)
  tanomaly=-999.
endwhere
crut4(:,:,1:crut4stop)=tanomaly(:,:,cruoffset:crustop)
deallocate(tanomaly)

end subroutine read_hadCRUT4


END MODULE txtnc
