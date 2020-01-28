
MODULE DATA_setting
!this module develope a derived type 'station' where it is possible to store teh data and then it manage the data and create time series and the short time series. 

USE netcdf

USE error_management

USE JSON_file

IMPLICIT NONE


!Global variables for the program 
TYPE,PUBLIC :: global_settings
     CHARACTER(len =100) ::title(4)= (/"RAOBCORE","RADIOSONDE INPUT DATA","ERA-INTERIM","www.univie.ac.at/theoret-met/research/raobcore"/)
     INTEGER :: n_days= 45000         !!max number of days = 45000
     INTEGER :: n_layers= 16          !!maximun number of pressure layer +1 for the error
     INTEGER :: n_obs=2               !!maximum number of observation for day
     INTEGER :: start_day = 19000101  !!fist day of observation
     INTEGER :: end_day   = 20193112  !!last day of observation
     INTEGER, ALLOCATABLE :: date(:)  !!useful to create the data index
     INTEGER :: ts_initialize = -999  !!value to initialize all the time series 
                                      !!or for missing value
     INTEGER, DIMENSION(5) :: varno_in = (/2,3,4,7,29/) 
                                   !!varno_in(1)= 2 Temperture           [k]
                                   !!varno_in(2)= 3 u component of wind  [m/s]
                                   !!varno_in(3)= 4 v component of wind  [m/s]
                                   !!varno_in(4)= 7 Specific umidity    [kg/kg]
                                   !!varno_in(5)=29 Realtive Humidity   [Numeric]

     INTEGER, DIMENSION(2):: time_obs = (/000000,120000/)!!observation at 0 and 12 hour
     
  INTEGER, DIMENSION(16):: pressure_layers = (/1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000,92500,100000/)

   !!100000 Pa -> press_layer = 1
   !! 92500 Pa -> press_layer = 2
   !! 85000 Pa -> press_layer = 3
   !! 70000 Pa -> press_layer = 4
   !! 50000 Pa -> press_layer = 5
   !! 40000 Pa -> press_layer = 6
   !! 30000 Pa -> press_layer = 7
   !! 25000 Pa -> press_layer = 8
   !! 20000 Pa -> press_layer = 9
   !! 15000 Pa -> press_layer =10
   !! 10000 Pa -> press_layer =11
   !!  7000 Pa -> press_layer =12
   !!  5000 Pa -> press_layer =13
   !!  3000 Pa -> press_layer =14
   !!  2000 Pa -> press_layer =15
   !!  1000 Pa -> press_layer =16
  
  
END TYPE global_settings


!make the station
TYPE, PUBLIC ::  station
   INTEGER :: date  ! date
   INTEGER :: time  ! time observation 00 or 12
   INTEGER :: obstype ! 5 = radiosonde
   INTEGER :: codetype !tipe of code = 35
   INTEGER :: varno !type of observation e
                    ! 2 = Upper Air Temperature [k]
                    ! 3 = U component of wind   [m/s]
                    ! 4 = V component of Wind   [m/s]
                    ! 7 = q Specific humidity   [kg/kg]  
                    !29 = RH Relative humidity  [Numeric]
   INTEGER :: ident !identification number
   INTEGER :: ident_NOwmo
   REAL    :: lat, long !latitude and longitude
   INTEGER :: stalt!statution altitude
   INTEGER :: press !pressure level
   REAL(KIND = 4) :: obsvalue !observated value
   REAL :: biascorrect !correction value
   REAL :: fg_depar !
   REAL :: an_depar !
   INTEGER :: status !
   INTEGER :: anflag!
   INTEGER :: event1!
   INTEGER :: sonde_type
END TYPE station

!this is the 20CR station
TYPE  :: twenty_CR
   INTEGER :: varno !type of observation e
                    ! 2 = Upper Air Temperature [k]
                    ! 3 = U component of wind   [m/s]
                    ! 4 = V component of Wind   [m/s]
                    ! 7 = q Specific humidity   [kg/kg]  
                    !29 = RH Relative humidity  [Numeric]
   REAL:: scale_factor
   REAL :: offset
   REAL, ALLOCATABLE , DIMENSION(:):: lat      !latitude as 180 point each 2 degrees (0 - 358)
   REAL, ALLOCATABLE , DIMENSION(:):: lon      !latitude as 91 point each 2 degrees  (90- 0 - -90)
   INTEGER, ALLOCATABLE , DIMENSION(:):: date  !date yyyy as time_index function
   INTEGER, ALLOCATABLE , DIMENSION(:):: hour  !0 6 12 18 as time_index function
   INTEGER, ALLOCATABLE , DIMENSION(:):: pressure_layers !in the 20 CR there are 24 press_lvl
   INTEGER*2, ALLOCATABLE,DIMENSION(:,:,:,:) :: short_data_PL!here I store the values as (lon,lat,level,time_index)function for the pressure level data
   INTEGER*2, ALLOCATABLE,DIMENSION(:,:,:) :: short_data_S!here I store the values as (lon,lat,time_index)function for the surface data
                                                    !they are as short integer that I will convert in REAL
END TYPE twenty_CR
!this is the 20CR station elaborated 
TYPE  :: twenty_CR2
   INTEGER :: varno !type of observation e
                    ! 2 = Upper Air Temperature [k]
                    ! 3 = U component of wind   [m/s]
                    ! 4 = V component of Wind   [m/s]
                    ! 7 = q Specific humidity   [kg/kg]  
                    !29 = RH Relative humidity  [Numeric]
   REAL:: scale_factor
   REAL :: offset
   REAL, ALLOCATABLE , DIMENSION(:):: lat      !latitude as 180 point each 2 degrees (0 - 358)
   REAL, ALLOCATABLE , DIMENSION(:):: lon      !latitude as 91 point each 2 degrees  (90- 0 - -90)
   INTEGER, ALLOCATABLE , DIMENSION(:):: date  !date yyyy as time_index function
   INTEGER, ALLOCATABLE , DIMENSION(:):: hour  !0 6 12 18 as time_index function
   INTEGER, ALLOCATABLE , DIMENSION(:):: pressure_layers !in the 20 CR there are 16 press_lvl
   INTEGER, ALLOCATABLE,DIMENSION(:,:,:,:) :: short_data!here I store the values as (lon,lat,level,time_index)function
                                                    !they are as short integer that I will convert in REAL
END TYPE twenty_CR2

TYPE  :: grid4_20CR
   INTEGER,dimension(5):: varno !type of observation 
                    ! 2  = Upper Air Temperature [k]
                    ! 0  = geopotential height
                    ! 22 = surface Termperature  [K]    
                    ! 3  = U component of wind   [m/s]
                    ! 4  = V component of Wind   [m/s]
                    ! 7  = q Specific humidity   [kg/kg]  
                    !29  = RH Relative humidity  [Numeric]
     REAL, DIMENSION (4) :: lat      !latitude as 4 points 
     REAL, DIMENSION (4) :: lon      !latitude as 4 points
     REAL:: scale_factor
     REAL :: offset
     !!from lat and lon combination 11-> base left, 22 -> top left, 33 -> base right, 44 -> top right
     INTEGER, ALLOCATABLE , DIMENSION(:):: date  !date yyyy as time_index function
     INTEGER, ALLOCATABLE , DIMENSION(:):: hour  !0 6 12 18 as time_index function
     INTEGER, ALLOCATABLE , DIMENSION(:):: pressure_layers !the standard 16 pressure layers
     REAL, ALLOCATABLE,DIMENSION(:,:,:,:) :: grids4_data!here I store the values as (varno_index,grid_point,level,time_index)function temperature(press_level)
                                                      !hgt(press_level)
     REAL, ALLOCATABLE,DIMENSION(:,:) :: grids4_data_sT!here I store the values as (grid_point,time_index)function-> surface Temperature
     REAL, ALLOCATABLE,DIMENSION(:,:,:) :: interpolated_data!here I store the values as (varno_index,level,time_index)function 
     REAL, ALLOCATABLE,DIMENSION(:) :: interpolated_data_ST!here I store the values as (time_index)function
END TYPE grid4_20CR

!make the timeseries (complete timeserie with the missing value)
TYPE, PUBLIC :: timeseries
     REAL, ALLOCATABLE , DIMENSION(:,:,:)::obs !!observation tipe
     REAL, ALLOCATABLE , DIMENSION(:,:,:):: biascorrect !correction valu
     REAL, ALLOCATABLE , DIMENSION(:,:,:):: fg_depar
     REAL, ALLOCATABLE , DIMENSION(:,:,:):: an_depar
     INTEGER, ALLOCATABLE , DIMENSION(:) :: sonde_type
     INTEGER, ALLOCATABLE , DIMENSION(:,:,:):: status
     INTEGER, ALLOCATABLE , DIMENSION(:,:,:):: anflag
     INTEGER, ALLOCATABLE , DIMENSION(:,:,:):: event1

END TYPE timeseries


!!make a compact timeserie with only the available and good data.
TYPE, PUBLIC :: compact_timeseries
     INTEGER :: g_days !! number of days with data
     INTEGER,ALLOCATABLE , DIMENSION(:)::     index !!observation tipe  
     REAL, ALLOCATABLE , DIMENSION(:,:,:)::   obs !!observation tipe 
     REAL, ALLOCATABLE , DIMENSION(:,:,:)::   biascorrect !correction valu
     REAL, ALLOCATABLE , DIMENSION(:,:,:)::   fg_depar
     REAL, ALLOCATABLE , DIMENSION(:,:,:)::   an_depar
     INTEGER , ALLOCATABLE , DIMENSION(:):: sonde_type
     INTEGER, ALLOCATABLE , DIMENSION(:,:,:)::status
     INTEGER, ALLOCATABLE , DIMENSION(:,:,:)::anflag
     INTEGER, ALLOCATABLE , DIMENSION(:,:,:)::event1

END TYPE compact_timeseries

!!make the complete_station
TYPE, PUBLIC :: complete_station   
   INTEGER :: ident !identification number
   REAL    :: lat, long !latitude and longitude
   INTEGER :: stalt!
   type(global_settings) :: g_s ! global settings used
   type(compact_timeseries), DIMENSION(5) :: all_ts !!here I can stored all the time series
   
END TYPE complete_station

!!to check if I have already this station complete
TYPE, PUBLIC :: existing_station   
   INTEGER :: date  !date
   INTEGER :: time=0  ! time observation 00 or 12
   INTEGER :: varno !type of observation e
                    ! 2 = Upper Air Temperature [k]
                    ! 3 = U component of wind   [m/s]
                    ! 4 = V component of Wind   [m/s]
                    ! 7 = q Specific humidity   [kg/kg]  
                    !29 = RH Relative humidity  [Numeric]
   INTEGER :: press =0  !pressure level
   
END TYPE existing_station


!!Here you can find:
!! a function to read the data file
!! a function to create the corret time-index
!! a function to convert index to date
!! a function to convert date to index
!! a function to initialize the timeserie
!! a function to build up the timeserie and the compact_timeserie
!! a function to build up the complete station type. Here are stored: 
              !! the general information(, stat_index,lat, long, statalt)
              !! and the compact_time_series in an array named ->all_ts


CONTAINS

INTEGER FUNCTION read_20CR_PL(current_year,current_varno,tw_CR_station,time_dim)
  !EXTERNAL
  INTEGER :: current_year
  INTEGER :: current_varno
  TYPE(twenty_CR):: tw_CR_station
  INTEGER :: time_dim !time dimension = numbers of items
  !INTERNAL
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
    INTEGER :: lon_dimid,lat_dimid,pres_l_dimid, time_dimid
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
    INTEGER :: i,j,k,z
    
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
    !path_name = "/home/srvx7/lorenzo/scratch/test_20CR/"//trim(file_name)
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
    status =nf90_inq_dimid(ncid, "lon", lon_dimid)
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error inquire Longitude dimension id"
       STOP
    end if
    !inquire dimension
    status=NF90_INQUIRE_DIMENSION(ncid,lon_dimid,lon_name,lon_dim)
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
    !inquire the variable
    status = nf90_inq_varid(ncid,"lon",lon_varid)
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error inquire Longitude variable id"
       STOP
    end if
    !get the variable
    status = nf90_get_var(ncid,lon_varid,tw_CR_station%lon)
    !status = nf90_get_var(ncid,3,longitude)
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error get Longitude variable"
       STOP
    end if
    
    !LATITUDE
    status =nf90_inq_dimid(ncid, "lat", lat_dimid)
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error inquire Latitude dimension id"
       STOP
    end if
    !inquire dimension
    status=NF90_INQUIRE_DIMENSION(ncid,lat_dimid,lon_name,lat_dim)
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error inquire Latitude  dimension"
       STOP
    end if
    ALLOCATE (tw_CR_station%lat(lat_dim), STAT = status)
    if(status .eq. 0)then
       tw_CR_station%lat=0
    else
       WRITE(*,*)"Allocate memory ERROR, variable ""latitude"""
       STOP
    end if
    !inquire the variable
    status = nf90_inq_varid(ncid,"lat",lat_varid)
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error inquire Latitude variable id"
       STOP
    end if   
    !get the variable
    status = nf90_get_var(ncid,lat_varid,tw_CR_station%lat)
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error get Latitude variable"
       STOP
    end if
    
    !PRESSURE LAYERS
    status =nf90_inq_dimid(ncid, "level", pres_l_dimid)
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error inquire Pressure dimension id"
       STOP
    end if
    !inquire dimension
    status=NF90_INQUIRE_DIMENSION(ncid,pres_l_dimid,pres_l_name,pres_l_dim)
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
     !inquire the variable
    status = nf90_inq_varid(ncid,"level",pres_l_varid)
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error inquire Pressure Layers variable id"
       STOP
    end if
    !get the variable
    status = nf90_get_var(ncid,pres_l_varid,tw_CR_station%pressure_layers)
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error get Pressure Layers variable"
       STOP
    end if

    !TIME
    status =nf90_inq_dimid(ncid, "time", time_dimid)
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error inquire Time dimension id"
       STOP
    end if
    !inquire dimension
    status=NF90_INQUIRE_DIMENSION(ncid,time_dimid,time_name,time_dim)
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
  !inquire the variable
  status = nf90_inq_varid(ncid,"time",time_varid)
  if(status .NE. NF90_NOERR) then
     PRINT*, "Error inquire Time variable id"
     STOP
  end if
  !get the variable
  status = nf90_get_var(ncid,time_varid,time_index)
  if(status .NE. NF90_NOERR) then
     PRINT*, "Error get Time variable"
     STOP
  end if
  !Now I need to convert this index in data and time
  !the correct data and time is stored in the tw_CR_station type
!  status = index_generate(file_name,time_dim,time_index(1),time_index(time_dim),tw_CR_station)

  !VARIABLE string_var= 'air','uwnd','vwnd', 'hgt'
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
  status = nf90_get_var(ncid,air_varid,tw_CR_station%short_data_PL)
  if(status .NE. NF90_NOERR) then
     PRINT*, NF90_STRERROR(status)
     STOP
  end if  
  
  status=nf90_get_att(ncid, air_varid, 'add_offset', tw_CR_station%offset)
  if(status .NE. NF90_NOERR) then
     PRINT*, NF90_STRERROR(status)
     STOP
  end if
 status=nf90_get_att(ncid, air_varid, 'scale_factor', tw_CR_station%scale_factor)
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
800   FORMAT(/, "      --> Reading NetCDF Files: ",A15," ----> OK")

read_20CR_PL = 0
END FUNCTION read_20CR_PL


INTEGER FUNCTION read_20CR_S(current_year,current_varno,tw_CR_station,time_dim)
  !EXTERNAL
  INTEGER :: current_year
  INTEGER :: current_varno
  TYPE(twenty_CR):: tw_CR_station
  INTEGER :: time_dim !time dimension = numbers of items
  !INTERNAL
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
    INTEGER :: lon_dimid,lat_dimid,pres_l_dimid, time_dimid
    CHARACTER (len = 50) :: lon_name,lat_name,lev_name,pres_l_name,time_name,air_name
    INTEGER :: lon_varid,lat_varid,lev_varid,time_varid,air_varid
    INTEGER :: lon_dim,lat_dim,lev_dim
    INTEGER :: lat,lon,lev,pres_l_varid,pres_l_dim,air_dim
    REAL :: longitude(180)               !longitude
    REAL :: latitude(91)                 !atitude
    REAL :: p_levels(24)                  !pressure levels
    REAL*8 , ALLOCATABLE,DIMENSION(:) :: time_index!!real double precision
    !REAL , ALLOCATABLE :: temp_arr(:,:,:,:)
    INTEGER*2, ALLOCATABLE, DIMENSION(:,:,:):: short_data_S !data from the 20CR as short integer
    INTEGER :: i,j,k,z
    !now I have to create the file name that I have to use for the path_name
    if(current_varno .eq. 22)then      
       string_var= 'air.sig995'
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
     status =nf90_inq_dimid(ncid, "lon", lon_dimid)
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error inquire Longitude dimension id"
       STOP
    end if
    !inquire dimension
    status=NF90_INQUIRE_DIMENSION(ncid,lon_dimid,lon_name,lon_dim)
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
    
    status = nf90_inq_varid(ncid,"lon",lon_varid)
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error inquire Longitude variable id"
       STOP
    end if
     status = nf90_get_var(ncid,lon_varid,tw_CR_station%lon)
    !status = nf90_get_var(ncid,3,longitude)
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error get Longitude variable"
       STOP
    end if
    
    !LATITUDE
     status =nf90_inq_dimid(ncid, "lat", lat_dimid)
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error inquire Latitude dimension id"
       STOP
    end if
    !inquire dimension
    status=NF90_INQUIRE_DIMENSION(ncid,lat_dimid,lat_name,lat_dim)
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
    status = nf90_inq_varid(ncid,"lat",lat_varid)
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error inquire Latitude variable id"
       STOP
    end if
       
    !get the variable 
    status = nf90_get_var(ncid,lat_varid,tw_CR_station%lat)
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error get Latitude variable"
       STOP
    end if
    
    !PRESSURE LAYERS
    !if Im working with surface temperature data
    !I haven t pressure levels
    pres_l_dim =1
 ALLOCATE (tw_CR_station%pressure_layers(pres_l_dim), STAT = status)
    if(status .eq. 0)then
       tw_CR_station%pressure_layers=0
    else
       WRITE(*,*)"Allocate memory ERROR, variable ""pressure_layers"""
       STOP
    end if
    tw_CR_station%pressure_layers = 0!pressure surface layer

    !TIME
    status =nf90_inq_dimid(ncid, "time", time_dimid)
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error inquire Time dimension id"
       STOP
    end if
     !inquire dimension

    status=NF90_INQUIRE_DIMENSION(ncid,time_dimid,time_name,time_dim)
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
    status = nf90_inq_varid(ncid,"time",time_varid)
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error inquire Time variable id"
       STOP
    end if
    
    !get the variable
    status = nf90_get_var(ncid,time_varid,time_index)
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error get Time variable"
       STOP
    end if
    !Now I need to convert this index in data and time
  !the correct data and time is stored in the tw_CR_station type
!    status = index_generate(file_name,time_dim,time_index(1),time_index(time_dim),tw_CR_station)
    
    !VARIABLE string_var= 'air'
    status = nf90_inq_varid(ncid,string_var(1:3),air_varid)
    if(status .NE. NF90_NOERR) then
     PRINT*, "Error inquire", string_var ," variable id"
     STOP
  end if

  !allocate the memory
 ! ALLOCATE (tw_CR_station%short_data_S(time_dim,lat_dim,lon_dim),STAT = status)
ALLOCATE (tw_CR_station%short_data_S(lon_dim,lat_dim,time_dim),STAT = status)
   if(status .NE. NF90_NOERR) then
     PRINT*, "Allocation error ->",TRIM(air_name)
     STOP
  end if
  tw_CR_station%short_data_S=0

  !get the variable
  status = nf90_inq_varid(ncid,"air",air_varid)
    if(status .NE. NF90_NOERR) then
       PRINT*, "Error inquire Time variable id"
       STOP
    end if
  status = nf90_get_var(ncid,air_varid,tw_CR_station%short_data_S)
  if(status .NE. NF90_NOERR) then
     PRINT*, NF90_STRERROR(status)
     STOP
  end if
  status=nf90_get_att(ncid, air_varid, 'add_offset', tw_CR_station%offset)
  if(status .NE. NF90_NOERR) then
     PRINT*, NF90_STRERROR(status)
     STOP
  end if
  status=nf90_get_att(ncid, air_varid, 'scale_factor', tw_CR_station%scale_factor)
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
800   FORMAT(/, "      --> Reading NetCDF Files: ",A15," ----> OK")

read_20CR_S = 0
END FUNCTION read_20CR_S


 !!~~~~~~~~~ FUNCTION TO READ NetCTF FILES~~~~~~~~~~~~!!

  !!NB this function needs as input an array file_name_in(5) where inside there are the all the file that will be read
  !!(of course the NetCDF files should be of the same station..they are maximum 5 because 5 are the varno we are interested in)
  !! net_complete_stat is the Fortran Type object reconstructed by the function from the Netcdf file/s

  INTEGER FUNCTION read_ncdf(file_name_in,net_complete_stat,routine_name, g_s, g_s_in)  
    IMPLICIT NONE
    !! DUMMY OBJECT declaration  
    ! This is the name of the data file we will read
    !!dummy object for the NetCDF files usefull to read the files
    CHARACTER(len= *) :: file_name_in(1)

    !!dummy object for the numer of Varno I m interested in
   ! INTEGER,dimension(5) :: Num_varno


    !!dummy objec for the complete_stat
    TYPE(complete_station)::net_complete_stat

    TYPE(global_settings_raobcore),OPTIONAL :: g_s,g_s_in !type used to store all the global setting

    INTEGER :: ncid !!NetCDF file Id number
               !!number of dimensions
                         !!number of variables
                                   !!Number of global attributes
    INTEGER :: N_dimsin, N_varsin, N_gatts_in

    INTEGER :: stat_index_varid    !!station index var id
    INTEGER :: lat_varid           !!latitude var id
    INTEGER :: lon_varid           !!longitude var id
    INTEGER :: alt_varid           !!altitude var id
    INTEGER :: pres_l_varid        !!pressure layers var id
    INTEGER :: pres_l_dim          !!pressure layers dimension
    CHARACTER(len = 50)::pres_l_name!!pressure layers name
    INTEGER :: obs_time_varid      !!observation time var id
    INTEGER :: obs_time_dim        !!observation time dimension 
    CHARACTER(len = 50) :: obs_time_name!!observation time name
    INTEGER :: date_varid          !!number of total day var id
    INTEGER :: date_dim            !!dimension of data variable
    CHARACTER(len = 50)::date_name !!data name

    INTEGER, DIMENSION(5):: varno_varid !!varno number var id
    CHARACTER(len = 50) :: name_varno(5)=(/"Varno_Temp","Varno_U_Wind","Varno_V_Wind","Varno_S_Humidity","Varno_R_Humidity"/)!!name varno variable 
    INTEGER :: index_days_varid       !!number of good and available days
    INTEGER :: index_dim             !!dimension of index variable
    CHARACTER(len = 50)::index_name  !!index name

    CHARACTER(len = 50)::string       !!string_global attribute

  !!Time series
    !!Observation
    INTEGER:: temp_varid
    INTEGER:: biasc_varid
    INTEGER:: fg_depar_varid
    INTEGER:: an_depar_varid
    !!Fags
    INTEGER:: status_varid
    INTEGER:: anflag_varid
    INTEGER:: event1_varid

    !Global attribute
    INTEGER :: conv_varid

    !!Name_variables
    CHARACTER(len = 20):: variable_names(5)= (/"Temperature","U_Wind","V_Wind","RH","SH"/)

    !!status variable declaration
    INTEGER :: status
    !!counters and dummy
    INTEGER :: i,j
    INTEGER ::count !It is necessary to initialize to -1
                       ! because insede the program I use the condition cout >0 as loop  exit 
    INTEGER, dimension(2) ::start_end=(/0,0/) !for the start and end day
   ! INTEGER ::dim
    !!number of file I will read
    INTEGER :: dim
    !!variables for the time initializing
    character (len = 20):: date, time, zone
    character(len = 50) ::today
    character(len = 200) ::source
    
    !internal variables
    CHARACTER(len = 200) ::routine_name,error_mex,temp_str
    !logical
    LOGICAL :: debugging_mode
    
    !program name
    routine_name="read_ncdf"
    debugging_mode=.false.
!!$    debug_check:do i=1,g_s_in%debug_mode%dim
!!$       if(g_s_in%debug_mode%routine_name(i) == routine_name)then
!!$          debugging_mode=.true.
!!$          if(debugging_mode == .true.)then
!!$             exit debug_check
!!$          endif
!!$       endif
!!$    enddo debug_check
    !It is necessary to initialize to -1
    ! because insede the program I use the condition cout >0 as loop  exit 
    count= -1
    !!big loop to write in the Net_Complete_station type
    j = 0
    put_Net_Stat:  do j = 1, size(file_name_in)
       
       !WRITE(*,*)"Filename(",j,") -->",trim(file_name_in(j)),"<--"
       !WRITE(*,*)"len_trim(Filename)",len_trim(file_name_in(j))
       
       if (trim(file_name_in(j)) .NE. "") then
          
          !open the file
!          !$OMP CRITICAL
          status = nf90_open(file_name_in(j), nf90_nowrite,ncid)
          
          error_mex= "Error opening file->"//file_name_in(j)
          call error(routine_name,status,2,error_mex)
          
          
          !inquire file:->
          !dimension Number, dimension Variables, dimension global setting
          status = nf90_inquire(ncid,N_dimsin, N_varsin, N_gatts_in )
          error_mex= "Error Inquire file->"//file_name_in(j)
          call error(routine_name,status,2,error_mex)
          
          
          !write(*,*)N_dimsin, N_varsin, N_gatts_in
          
          !!NOW IPUT THE GLOBAL ATTRIBUTE IN THE NET_STATION
          status =  nf90_get_att(ncid,nf90_global,"title",g_s%title(1))
          error_mex= "Error writing Global Attribute ->"//"title"
          call error(routine_name,status,2,error_mex)

          status =  nf90_get_att(ncid,nf90_global,"source",g_s%title(2))
          error_mex= "Error writing Global Attribute ->"//"source"
          call error(routine_name,status,2,error_mex)
          
          status =  nf90_get_att(ncid,nf90_global,"datatype",g_s%title(3))
          error_mex= "Error writing Global Attribute ->"//"datatype"
          call error(routine_name,status,2,error_mex)

          status =  nf90_get_att(ncid,nf90_global,"references",g_s%title(4))
          error_mex= "Error writing Global Attribute ->"//"references"
          call error(routine_name,status,2,error_mex)


      !inquire the variable identification number
      !and inquire the variable to put the value in to the station Type
      !Station Index
      global_set:do while(count < 0) ! if (j .EQ.1) then
         count = 1
         status = nf90_inq_varid(ncid,"station",stat_index_varid)
         error_mex= "Error inquire station variable index"
         call error(routine_name,status,2,error_mex)
         
         status = nf90_get_var(ncid,stat_index_varid,net_complete_stat%ident)
         error_mex= "Error get Station index variable in Complete Station Type"
         call error(routine_name,status,2,error_mex)
         
          !Latitude
         status = nf90_inq_varid(ncid,"lat",lat_varid)
         error_mex= "Error inquire Latitude variable id"
         call error(routine_name,status,2,error_mex)

         status = nf90_get_var(ncid,lat_varid,net_complete_stat%lat)
         error_mex= "Error get Latitude variable in Complete Station Type"
         call error(routine_name,status,2,error_mex)
         
         !Longitude
         status = nf90_inq_varid(ncid,"lon",lon_varid)
         error_mex= "Error inquire Longitude variable id"
         call error(routine_name,status,2,error_mex)

         status = nf90_get_var(ncid,lon_varid,net_complete_stat%long)
         error_mex= "Error get Longitude variable in Complete Station Type"
         call error(routine_name,status,2,error_mex)

         
         !Altitude
         status = nf90_inq_varid(ncid,"alt",alt_varid)
         error_mex= "Error inquire Altitude variable id"
         call error(routine_name,status,2,error_mex)

         status = nf90_get_var(ncid,alt_varid,net_complete_stat%stalt)
         error_mex= "Error get Altitude variable in Complete Station Type"
         call error(routine_name,status,2,error_mex)
         
         !Pressure Layers
         status = nf90_inq_varid(ncid,"pressure_layers",pres_l_varid)
         error_mex= "Error inquire Pressure Layers variable id"
         call error(routine_name,status,2,error_mex)
         
         !inquire dimension
         status=NF90_INQUIRE_DIMENSION(ncid,pres_l_varid,pres_l_name,pres_l_dim)
         error_mex= "Error inquire Pressure Layers  dimension"
         call error(routine_name,status,2,error_mex)

         status = nf90_get_var(ncid,pres_l_varid,g_s%pressure_layers)
         error_mex= "Error get Pressure Layers variable in Complete Station Type"
         call error(routine_name,status,2,error_mex)
        
         
         !Observation Time
         status = nf90_inq_varid(ncid,"obs_time",obs_time_varid)
         error_mex= "Error inquire Observation Time variable id"
         call error(routine_name,status,2,error_mex)

         !inquire dimension
         status=NF90_INQUIRE_DIMENSION(ncid,obs_time_varid,obs_time_name,obs_time_dim)
         error_mex= "Error inquire Observation Time dimension"
         call error(routine_name,status,2,error_mex)

         status = nf90_get_var(ncid,obs_time_varid,g_s%time_obs)
         error_mex= "Error get Observation Time variable in Complete Station Type"
         call error(routine_name,status,2,error_mex)

         !Date
         status = nf90_inq_varid(ncid,"date",date_varid)
         error_mex= "Error inquire Date variable id"
         call error(routine_name,status,2,error_mex)

         !inquire dimension
         status=NF90_INQUIRE_DIMENSION(ncid,date_varid,date_name,date_dim)
         error_mex= "Error inquire Date dimension"
         call error(routine_name,status,2,error_mex)

         !allocate the memory in the Station Type
         ALLOCATE (g_s%time_window%date(date_dim), STAT = status)
         WRITE(temp_str,'(I1)')j
         error_mex= "Allocation error in Station Type->"//TRIM(date_name)// "Loop:"//temp_str(1:1)
         call error(routine_name,status,2,error_mex)
         
         !get the variable in the Statyion Type
         status = nf90_get_var(ncid,date_varid,g_s%time_window%date)
         error_mex= "Error get Date variable in Complete Station Type"
         call error(routine_name,status,2,error_mex)

         !star day
         status = nf90_get_att(ncid,date_varid,"valid_range",start_end)
         error_mex= "Error get Satrt Day variable in Complete Station Type"
         call error(routine_name,status,2,error_mex)
         g_s%time_window%start_day =  start_end(1)
         g_s%time_window%end_day  = start_end(2)

         !Varno
         i = 0 
         do i = 1,5
            status = nf90_inq_varid(ncid,name_varno(i),varno_varid(i))
            WRITE(temp_str,'(I1)')
            error_mex="Error inqure Varno variable id : "//temp_str
            call error(routine_name,status,2,error_mex)
            status = nf90_get_var(ncid,varno_varid(i),g_s%varno_in(i))
            WRITE(temp_str,'(I1)')
            error_mex="Error get arno variable in Complete Station Type : "//temp_str
            call error(routine_name,status,2,error_mex)
           
         end do
      end do global_set
      !ed if global_set
     


      !!INTERNAL STRUCTURE : DIFFERENCE FOR EACH TIME SERIE


      !Index Days
      status = nf90_inq_varid(ncid,"index_days",index_days_varid)
      WRITE(temp_str,'(I1)')j
      error_mex="Error inquire Index Days variable id"// "Varno_index = "//temp_str
      call error(routine_name,status,2,error_mex)

      !inquire dimension
      status=NF90_INQUIRE_DIMENSION(ncid,index_days_varid,index_name,index_dim)
      WRITE(temp_str,'(I1)')j
            error_mex="Error inquire Index dimension"// "Varno_index = "//temp_str
            call error(routine_name,status,2,error_mex)
      
      !!the index dimension is the good days number
      net_complete_stat%all_ts(j)%g_days = index_dim
      
!allocate the memory in the Station Type
      ALLOCATE (net_complete_stat%all_ts(j)%index(index_dim), STAT = status)
      WRITE(temp_str,'(I1)')j
      error_mex= "Allocation error in Station Type->"//index_name//"Varno_index = "//temp_str
      call error(routine_name,status,2,error_mex)
      status = nf90_get_var(ncid,index_days_varid,net_complete_stat%all_ts(j)%index)
      WRITE(temp_str,'(I1)')j
      error_mex="Error get Index Days variable in Complete Station Type"// "Varno_index = "//temp_str
      call error(routine_name,status,2,error_mex)
      

      !!TIME SERIES!!
      !!OBSERVATIONS
      
      !OBSERVATION TYPE -> TEMPERATURE; U_WIND; V_WIND; RH; SH
      !inquire temperature variable
      status = nf90_inq_varid(ncid,variable_names(j),temp_varid)
      error_mex="Error inquire "//variable_names(j)//" variable id"
      call error(routine_name,status,2,error_mex)
      
      !allocate the memory in the Station Type
      ALLOCATE (net_complete_stat%all_ts(j)%obs(index_dim,pres_l_dim,obs_time_dim ), STAT = status)
      error_mex="Allocation error in Station Type->"//variable_names(j)
      call error(routine_name,status,2,error_mex)
      
      !get variable Temperature
      status = nf90_get_var(ncid,temp_varid,net_complete_stat%all_ts(j)%obs)
      error_mex= "Error get "//variable_names(j)//" variable in Complete Station Type"
      call error(routine_name,status,2,error_mex)
      
      !Biascorrect
      !inquire biascorrect variable
      status = nf90_inq_varid(ncid,"biascorrect",biasc_varid)
      error_mex="Error inquire Biascorrect variable id"
      call error(routine_name,status,2,error_mex)
      
      !allocate the memory in the Station Type
      ALLOCATE (net_complete_stat%all_ts(j)%biascorrect(index_dim,pres_l_dim,obs_time_dim ), STAT = status)
      error_mex="Allocation error in Station Type->"//"Biascorrect"
      call error(routine_name,status,2,error_mex)
      
      !get variable Biascorrect
      status = nf90_get_var(ncid,biasc_varid,net_complete_stat%all_ts(j)%biascorrect)
      error_mex= "Error get Temperature variable in Complete Station Type"
      call error(routine_name,status,2,error_mex)
      

      !Fg_depar
      !inquire fg_depar variable
      status = nf90_inq_varid(ncid,"fg_depar",fg_depar_varid)
      error_mex= "Error inquire Fg_depar variable id"
      call error(routine_name,status,2,error_mex)
      
      
      !allocate the memory in the Station Type
      ALLOCATE (net_complete_stat%all_ts(j)%fg_depar(index_dim,pres_l_dim,obs_time_dim ), STAT = status)
      error_mex= "Allocation error in Station Type->"//"Fg_depar"
      call error(routine_name,status,2,error_mex)
      
      !get variable Fg_depar
      status = nf90_get_var(ncid,fg_depar_varid,net_complete_stat%all_ts(j)%fg_depar)
      error_mex= "Error get Fg_depar variable in Complete Station Type"
      call error(routine_name,status,2,error_mex)
      
      !An_depar
      !inquire an_depar variable
      status = nf90_inq_varid(ncid,"an_depar",an_depar_varid)
      error_mex= "Error inquire An_depar variable id"
      call error(routine_name,status,2,error_mex)
      
      !allocate the memory in the Station Type
      ALLOCATE (net_complete_stat%all_ts(j)%an_depar(index_dim,pres_l_dim,obs_time_dim ), STAT = status)
      error_mex= "Allocation error in Station Type->"//"An_depar"
      call error(routine_name,status,2,error_mex)
      
      !get variable Fg_depar
      status = nf90_get_var(ncid,an_depar_varid,net_complete_stat%all_ts(j)%an_depar)
      error_mex="Error get An_depar variable in Complete Station Type"
      call error(routine_name,status,2,error_mex)      
      
      !!FLAGS
      
      !Status
      !inquire status variable
      status = nf90_inq_varid(ncid,"status",status_varid)
      error_mex= "Error inquire Status variable id"
      call error(routine_name,status,2,error_mex)
      
      !allocate the memory in the Station Type
      ALLOCATE (net_complete_stat%all_ts(j)%status(index_dim,pres_l_dim,obs_time_dim ), STAT = status)
      error_mex= "Allocation error in Station Type->"//"Status"
      call error(routine_name,status,2,error_mex)
      
      !get variable status
      status = nf90_get_var(ncid,status_varid,net_complete_stat%all_ts(j)%status)
      error_mex="Error get Status variable in Complete Station Type"
      call error(routine_name,status,2,error_mex)
      
      !Anflag
      !inquire  Anflag variable
      status = nf90_inq_varid(ncid,"anflag",anflag_varid)
      error_mex="Error inquire Anflag variable id"
      call error(routine_name,status,2,error_mex)

      !allocate the memory in the Station Type
      ALLOCATE (net_complete_stat%all_ts(j)%anflag(index_dim,pres_l_dim,obs_time_dim ), STAT = status)
      error_mex="Allocation error in Station Type->"//"Anflag"
      call error(routine_name,status,2,error_mex)
      
      !get variable anflag
      status = nf90_get_var(ncid,anflag_varid,net_complete_stat%all_ts(j)%anflag)
      error_mex= "Error get Anflag_varid variable in Complete Station Type"
      call error(routine_name,status,2,error_mex)
      
      !Event1
      !inquire event1  variable
      status = nf90_inq_varid(ncid,"event1",event1_varid)
      error_mex="Error inquire Anflag variable id"
      call error(routine_name,status,2,error_mex)

      !allocate the memory in the Station Type
      ALLOCATE (net_complete_stat%all_ts(j)%event1(index_dim,pres_l_dim,obs_time_dim ), STAT = status)
      error_mex="Allocation error in Station Type->"//"event1"
      call error(routine_name,status,2,error_mex)

      !get variable anflag
      status = nf90_get_var(ncid,event1_varid,net_complete_stat%all_ts(j)%event1)
      error_mex="Error get event1 variable in Complete Station Type"
      call error(routine_name,status,2,error_mex) 
     
      
      !Close file
      status=NF90_CLOSE(ncid)
      error_mex="Closing Error, with file: "//file_name_in(j)
      call error(routine_name,status,2,error_mex)
!      !$OMP END CRITICAL
      if(debugging_mode)then
900      FORMAT(/, "      --> Reading NetCDF Files: ...",A25,"----> OK")
         
         if (len_trim(file_name_in(j)).gt. 25)then
            write(*,900)file_name_in(j)(len_trim(file_name_in(j))-25:len_trim(file_name_in(j)))
            !900   FORMAT(/, "      --> Reading NetCDF Files: ",A25,"----> OK") 
         else
            write(*,900)file_name_in(j)
         end if
         
      end if
   endif
   

end do put_Net_Stat
   
   !!I return status function value
   read_ncdf = status

 END FUNCTION read_ncdf
 


END MODULE DATA_setting
