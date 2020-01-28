!This modure write and read a JSON file
! and organize a derived type where store the Json file infortation. 

MODULE JSON_file

 USE error_management ! file: error_management.f90
                      !In this module are defined the error messages used
                      !in the program.

!Description
!The module is built for:
!                         read from a JSON file the input parameters for the program and store them in the so called global_setting structure;
!                         write the so called global_setting structure in a file respecting the JSON format



!INPUT FILES:
! It needs only the global setting file, as json file

!OUTPUT FILES:

!Current Code Owner: Haimberger Leopold and Ramella Pralungo Lorenzo
!                    

!History:
! Version     DATE          Comment
! ------    ----------      -------
! 1.0       27.08.2011      Haimberger Leopold and Ramella Pralungo Lorenzo



!Code Description
! Language : Fortran 90.
!Software Standards: "European Standards for Writing and
!                     Documenting Exchangeable Fortran 90 Code".

!Global variables for the program 

  TYPE, PUBLIC :: stat_info
     INTEGER :: start=0
     INTEGER :: end=099999
     INTEGER :: number_of_station=3070
     INTEGER :: available_stations=5423
  END TYPE stat_info
  
  TYPE,PUBLIC ::swich_archive
     INTEGER :: ERA_INT = 19890101
     INTEGER :: ERA_FT = 19580101
     INTEGER :: tw_CR = 19000101
  END TYPE swich_archive

  TYPE,PUBLIC :: sat_time_window
     INTEGER :: start_day = 19710101  !!fist day of satellite adjustment
     INTEGER :: end_day   = 19781231  !!last day of satellite adjustment
  END TYPE sat_time_window

  TYPE,PUBLIC :: time_window
     INTEGER :: start_day = 19000101  !!fist day of observation
     INTEGER :: end_day   = 20191231  !!last day of observation
     INTEGER, ALLOCATABLE :: date(:)  !!useful to create the data index
     INTEGER, ALLOCATABLE :: year(:)  !!useful to create the data index
     INTEGER, ALLOCATABLE :: month(:)  !!useful to create the data index
     INTEGER, ALLOCATABLE :: day(:)  !!useful to create the data index
     INTEGER, ALLOCATABLE :: n_day_x_month(:)!useful to create the data index
     INTEGER :: first_month = 1
     INTEGER :: last_month = 1350
     TYPE(swich_archive) :: Swich_archive
     TYPE(sat_time_window):: Sat_time_window
     TYPE(sat_time_window):: NOAA_4
     INTEGER :: archive = 1 !1=20CR, 2=ERA40, 3=ERAINT, 12=20CR+ERA40, 13=20CR+ERAINT,123=20CR+ERA40+ERAINT
  END TYPE time_window
  

  
  TYPE,PUBLIC :: path
     CHARACTER(len =200) :: MERGED = "/vgc/srvx7/lorenzo/scratch/RAOBCORE_ARCHIVE/MERGED_ARCHIVE"
     CHARACTER(len =200) :: ERA_interim = "/vgc/srvx7/lorenzo/scratch/RAOBCORE_ARCHIVE/ERA_INTERIM"
     CHARACTER(len =200) :: ERA_forty = "/vgc/srvx7/lorenzo/scratch/RAOBCORE_ARCHIVE/ERA_40"
     CHARACTER(len =200) :: ERA_forty_analysis= "/home/srvx7/lorenzo/scratch/odbs/Stationfiles/ERA_40_analysis"
     CHARACTER(len =200) :: CHUAN = "/vgc/srvx7/lorenzo/scratch/RAOBCORE_ARCHIVE/CHUAN"
     CHARACTER(len =200) :: tw_CR = "/vgc/srvx7/lorenzo/scratch/RAOBCORE_ARCHIVE/LONG_20CR/JOINED_LONG_twCR/20CR"
     CHARACTER(len =200) :: output = "/vgc/srvx7/lorenzo/scratch/experiments/test"
END TYPE path
  
  TYPE, PUBLIC :: dimensions
     INTEGER :: n_days = 45000 !number of days
     INTEGER :: n_layers=16 !pressure levels
     INTEGER :: n_obs=2 !obstervation time
     INTEGER :: prob_max=16 ! probability time series
     INTEGER :: n_metadata=9 !number of metadata used
     INTEGER :: montly_threshold=9 ! minimum number of days to be considered a good month
     INTEGER :: max_expected_breaks=10!maximum number of expected breaks for each timeserie
     INTEGER :: min_n_days_montly_mean=10!minimum number of day to calculate montly_mean
     END TYPE dimensions

     TYPE, PUBLIC ::rem_clima
        INTEGER :: swich=0
        INTEGER :: minimum_n_year=5
     END TYPE rem_clima
     
TYPE, PUBLIC :: SNHT_param
   INTEGER :: increment=30 
   INTEGER :: max_len=1460 !max len of analysis windows 0> 2 years is the standard
   INTEGER :: mean_maxlen=2920
   INTEGER :: RAOB_max_miss = 650!max number of missing values => RAOBCORE
   !INTEGER :: RICH_max_miss = 450!max number of missing values => RICH
END TYPE  SNHT_param

TYPE, PUBLIC:: break_param
   !INTEGER:: threshold=50 !break threshold value
   INTEGER:: threshold_T=35 !break threshold value for Temperature
   INTEGER:: threshold_WD=18 !break threshold value for WindDirection
   INTEGER:: threshold_WS=50 !break threshold value for WindDirection
   REAL :: threshold_probability=0.5 !reak threshold probability
   INTEGER :: threshold_rad=23 !day-night differece => radiation
   REAL :: fak=0.7
   !REAL :: ap_prob_default=0!!!!!!!!OLD TO BE DELETED
   REAL :: ap_prob_default_T = 0.02 ! a priori probability => TEMPERATURE
   REAL :: ap_prob_default_WDWS = 0.075 ! a priori probability => WIND DIRECTION & SPEED
   INTEGER :: locsig=600 ! minimum distance between 2 break points (in samples)
   INTEGER :: br_max = 50 !max numbers of break points
   REAL :: smooth_fit_threshold = 0.2
   INTEGER :: bg_correction_factor = 1 !or 0 or 1
   INTEGER :: qc = 1 !quality control; 0=> NO; 1=> YES
   INTEGER :: old = 20000101 !initial adjustment ends here
   !TYPE(sat_time_window):: sat_time_wind
   INTEGER :: use_meta = 1 !use of metadata could be 0 1 2
   REAL :: sig_thresh_T = 0.5
   REAL :: sig_thresh_WS = 2.0
   REAL :: sig_thresh_mean_WD = 3.0
   REAL :: sig_thresh_rad = 0.3
   INTEGER :: bg_hom= 2 ! it tells the box in which we calculate the global departures
   INTEGER :: wide_break_detection_time_fac=3
   INTEGER :: fine_break_detection_time_fac=8
END TYPE break_param

TYPE,PUBLIC :: RAOBCORE_param
   INTEGER :: smooth_method=2
   INTEGER :: prob_method = 1 
END TYPE RAOBCORE_param

TYPE,PUBLIC :: RICH_param
      INTEGER :: RICH_max_miss = 450!max number of missing values => RICH
      INTEGER :: weight_distance = 3000 ! 
      INTEGER :: minst =2
      INTEGER ::ens = 0 
END TYPE RICH_param

TYPE, PUBLIC ::g_sonde_t
   INTEGER :: dim=27
   INTEGER, ALLOCATABLE::good_sondes(:)
END TYPE g_sonde_t

TYPE, PUBLIC :: priority_time_series
   CHARACTER(len=50) :: ts_name(12)=(/"fg_dep", "an_dep", "tropo", "strato","obs_day_night_diff","fg_day_night_diff","an_day_night_diff","composite","","","",""/)
   INTEGER, DIMENSION(12) :: ts_priority=(/5,1,2,3,4,0,0,0,0,0,0,0/)
END TYPE priority_time_series

TYPE, PUBLIC :: wind_param
   INTEGER:: min_WS_for_WD_dep=1
   INTEGER::break_det_method_mean_vs_3max=1
   INTEGER, DIMENSION(2) :: wind_press_layers=(/6,16/)
   INTEGER:: min_num_level_00_12=5
   REAL :: min_mean_WD_profile=3
   REAL :: st_dev_vertical_profile=1.5
END type wind_param

TYPE, PUBLIC :: d_mode
   INTEGER :: dim=1
   CHARACTER (len=200) , ALLOCATABLE::routine_name(:)
END TYPE d_mode


TYPE,PUBLIC :: global_settings_raobcore
     CHARACTER(len =100) ::title(4)= (/"RAOBCORE","RADIOSONDE DATA","homogenized using 20CR","www.univie.ac.at/theoret-met/research/raobcore"/)
     !INTEGER :: n_days= 40000         !!max number of days = 40000
     !INTEGER :: n_layers= 16          !!maximun number of pressure layer +1 for the error
     !INTEGER :: n_obs=2               !!maximum number of observation for day

     !INTEGER, ALLOCATABLE :: date(:)  !!useful to create the data index
     INTEGER :: ts_initialize = -999  !!value to initialize all the time series 
                                      !!or for missing value
     INTEGER, DIMENSION(7) :: varno_in = (/2,3,4,7,29,111,112/) 
                                   !!varno_in(1)= 2  Temperture           [k]
                                   !!varno_in(2)= 3  u component of wind  [m/s]
                                   !!varno_in(3)= 4  v component of wind  [m/s]
                                   !!varno_in(4)= 7  Specific umidity    [kg/kg]
                                   !!varno_in(5)=29  Realtive Humidity   [Numeric]
                                   !!varno_in(6)=111 WindDirection       [°]
                                   !!varno_in(7)=112 WindSpeed       [m/s]
     INTEGER, DIMENSION(7) :: var_check=(/1,1,1,0,0,1,1/)
     INTEGER, DIMENSION(2) :: time_obs = (/000000,120000/)!!observation at 0 and 12 hour
     
  INTEGER, DIMENSION(16) :: pressure_layers = (/1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,70000,85000,92500,100000/)

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
  TYPE(stat_info) ::Stat_info
  TYPE(dimensions)::Dimensions
  TYPE(time_window)::Time_window
  TYPE(rem_clima):: Remove_clima
  TYPE(SNHT_param)::SNHT_Param 
  TYPE(break_param) :: Break_Param
  TYPE(RICH_param) :: RICH_param
  TYPE(RAOBCORE_param) :: RAOBCORE_param
  TYPE(path) :: Path
  TYPE(g_sonde_t) :: Good_sondes_type
  INTEGER, DIMENSION(16) :: Pressure_layers_weights=(/1,1,1,1,1,1,1,1,0,1,1,1,1,1,0,0/)
  INTEGER, DIMENSION(16) :: Pressure_layers_smooth_weights=(/1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0/)
  REAL, DIMENSION(7) :: Max_break_size =(/15.0,0.0,0.0,15.0,15.0,60.0,30.0/)
  REAL, DIMENSION(16) ::Reduction_factor=(/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.5,0.0,0.0,0.0/)
  TYPE(wind_param)::Wind_parameters
  INTEGER :: press_lats_adj_most_recent=850
  TYPE(priority_time_series):: priority_time_series
  CHARACTER(len=2)::innovation="NN"
  CHARACTER(len=1)::extended_output="n"
  INTEGER :: max_iter =2
  TYPE(d_mode):: debug_mode
  
  
END TYPE global_settings_raobcore



contains
  

!!~~~~~~~~~ FUNCTION TO WRITE JSON  FILES~~~~~~~~~~~~!!

  INTEGER FUNCTION write_json(filename, routine_name,g_s)
    
    IMPLICIT NONE
    
    TYPE(global_settings_raobcore) :: g_s !! global settings variables
    
    CHARACTER(len = * ) :: filename !file name
    !CHARACTER(len = 350) :: line    
    INTEGER :: status               !opening status
    INTEGER :: n_lines              !number of lines


    INTEGER :: i,ii, lng                    !counters
    CHARACTER(len = 50) :: string        !only for internal operation
    CHARACTER(len = 50) :: string_s        !only for internal operation
    CHARACTER(len = 150):: temp          !only for internal operation
    !CHARACTER(len = 150)::aus,aus1,aus2 !only for internal operation
                                         !not used
    CHARACTER (len = 50) ::cform         !format name
    !filename = "jsonwrite.json"
    INTEGER:: long_s,long_ss
    CHARACTER (len = 200) ::  routine_name,mex
    CHARACTER(len=2) :: F1, F2
    !The routine name
    routine_name ="write_json"
    

    WRITE(*,*)""
    
    WRITE(*,*)"       -->I am working to write the Json file: ", trim(filename)
    
    OPEN(UNIT = 1, FILE = trim(filename), STATUS = "REPLACE", ACTION = "WRITE", IOSTAT= status)
    file: if(status .EQ. 0) then 
       !a json file start with "{"
       WRITE(1,*)"{"
       !PARAMETER!
       !!TITLE
       WRITE(1,*)"""title""",":"
       do i= 1,4
          if (i .EQ. 1) then
             WRITE(1,*)"[" 
          end if
          if(i .EQ. 4) then
             WRITE(1,*)"""", trim(g_s%title(i)),"""","]" 
          else  
             WRITE(1,*)"""", trim(g_s%title(i)),"""",","
          end if
       end do
       WRITE(1,*)","
       
 
       

       
       ! FUNCTION write_int_to_json_string(unit,title_str,body_str,body_end)
       WRITE(1,*)'"stat_info"',':'
       WRITE(1,*)"[" 
       status = write_int_to_json_string(1,trim('start'),g_s%Stat_info%start,0)
       status = write_int_to_json_string(1,trim('end'),g_s%Stat_info%end,0)
       status = write_int_to_json_string(1,trim('number_of_station'),g_s%Stat_info%number_of_station,0)
       status = write_int_to_json_string(1,trim('available_stations'),g_s%Stat_info%available_stations,1)
       WRITE(1,*) "]" 
       WRITE(1,*)","

       !DIMENSIONS
        WRITE(1,*)'"Dimensions"',':'
        WRITE(1,*)"[" 
        status = write_int_to_json_string(1,trim('n_days'),g_s%Dimensions%n_days,0)
        status = write_int_to_json_string(1,trim('n_layers'),g_s%Dimensions%n_layers,0)
        status = write_int_to_json_string(1,trim('n_obs'),g_s%Dimensions%n_obs,0)
        status = write_int_to_json_string(1,trim('prob_max'),g_s%Dimensions%prob_max,0)
        status = write_int_to_json_string(1,trim('n_metadata'),g_s%Dimensions%n_metadata,0)
        status = write_int_to_json_string(1,trim('montly_threshold'),g_s%Dimensions%montly_threshold,0)
        status = write_int_to_json_string(1,trim('max_expected_breaks'),g_s%Dimensions%max_expected_breaks,0)
        status = write_int_to_json_string(1,trim('min_n_days_montly_mean'),g_s%Dimensions%min_n_days_montly_mean,1)

        WRITE(1,*) "]" 
        WRITE(1,*)","

        !TIME WINDOW
        WRITE(1,*)'"Time_Window"',':'
        WRITE(1,*)"["
        status = write_int_to_json_string(1,trim('start_day'),g_s%Time_window%start_day,0)
        status = write_int_to_json_string(1,trim('end_day'),g_s%Time_window%end_day,0)
        status = write_int_to_json_string(1,trim('first_month'),g_s%Time_window%first_month,0)
        status = write_int_to_json_string(1,trim('last_month'),g_s%Time_window%last_month,0)
        status = write_int_to_json_string(1,trim('swich_CR'),g_s%Time_window%Swich_archive%tw_CR,0)
        status = write_int_to_json_string(1,trim('swich_ERA_FT'),g_s%Time_window%Swich_archive%ERA_FT,0)
        status = write_int_to_json_string(1,trim('swich_ERA_INT'),g_s%Time_window%Swich_archive%ERA_INT,0)
        status = write_int_to_json_string(1,trim('Sat_time_window_start_day'),g_s%Time_window%Sat_time_window%start_day,0)
        status = write_int_to_json_string(1,trim('Sat_time_window_end_day'),g_s%Time_window%Sat_time_window%end_day,0)
        g_s%Time_window%NOAA_4%start_day=19750101
        status = write_int_to_json_string(1,trim('NOAA_4_start_day'),g_s%Time_window%NOAA_4%start_day,0)
        g_s%Time_window%NOAA_4%end_day=19760901
        status = write_int_to_json_string(1,trim('NOAA_4_end_day'),g_s%Time_window%NOAA_4%end_day,0)
        status = write_int_to_json_string(1,trim('archive'),g_s%Time_window%archive,1)
        
        WRITE(1,*) "]" 
        WRITE(1,*)","

        !GENERAL INFO
        
        !!TIME SERIE INITIALIZE
        WRITE(1,*)'"time_serie_initialize"',':'
        WRITE(1,*)"["
        status = write_int_arr_line_to_json_string(1,g_s%ts_initialize,1)
        WRITE(1,*) "]" 
        WRITE(1,*)","
        
        !!VARNO IN
        WRITE(1,*)'"Varno_in"',':'
        WRITE(1,*)"["
        do i= 1,size(g_s%varno_in)
           if(i .ne. size(g_s%varno_in) ) then
              status = write_int_arr_line_to_json_string(1,g_s%varno_in(i),0)   
           else
              status = write_int_arr_line_to_json_string(1,g_s%varno_in(i),1)
           end if
        end do
        WRITE(1,*) "]" 
        WRITE(1,*)","
        !!VARNO CHECK 
        WRITE(1,*)'"Var_check"',':'
        WRITE(1,*)"["
        do i= 1,size(g_s%var_check)
           if(i .ne. size(g_s%var_check)) then
              status = write_int_arr_line_to_json_string(1,g_s%var_check(i),0)   
           else
              status = write_int_arr_line_to_json_string(1,g_s%var_check(i),1)
           end if
        end do
        WRITE(1,*) "]" 
        WRITE(1,*)","
        !!TIME OBS
        WRITE(1,*)'"Observation_time"', ":"
        WRITE(1,*)"[" 
       do i= 1,2
   

       if(i .eq. 2) then
           status = write_int_arr_line_to_json_string(1,g_s%time_obs(i),1)
       else
             status = write_int_arr_line_to_json_string(1,g_s%time_obs(i),0)
          end if
       end do
       WRITE(1,*) "]"
       WRITE(1,*)","


       !!TIME Pressure Layers
       WRITE(1,*)'"Pressure_layers"', ":"
       WRITE(1,*)"["
       do i= 1,g_s%Dimensions%n_layers
          if(i .EQ. 16) then
             status = write_int_arr_line_to_json_string(1,g_s%pressure_layers(i),1)
          else
        status = write_int_arr_line_to_json_string(1,g_s%pressure_layers(i),0)
     end if
  end do
  WRITE(1,*) "]"
    WRITE(1,*)","
    !Path
  WRITE(1,*)"""paths""",":"
  WRITE(1,*)"["
  ! FUNCTION write_string_to_json_string(unit,title_str,body_str,body_end)
  status = write_string_to_json_string(1,trim('MERGED'),trim(g_s%Path%MERGED),0)
  status = write_string_to_json_string(1,trim('ERA_interim'),trim(g_s%Path%ERA_interim),0)
  status = write_string_to_json_string(1,trim('ERA_forty'),trim(g_s%path%era_forty),0)
  status = write_string_to_json_string(1,trim('ERA_forty_analysis'),trim(g_s%path%era_forty_analysis),0)
  status = write_string_to_json_string(1,trim('CHUAN'),trim(g_s%path%chuan),0)
  status = write_string_to_json_string(1,trim('20CR'),trim(g_s%path%tw_cr),0)
  status = write_string_to_json_string(1,trim('OUTPUT'),trim(g_s%path%output),1)
  WRITE(1,*)"]"
  WRITE(1,*)","
  
  !REMOVE CLIMA
  WRITE(1,*)"""Remove_clima""", ":"
  WRITE(1,*)"[" 
  status = write_int_to_json_string(1,trim('Swich'),g_s%remove_clima%swich,0)
  status = write_int_to_json_string(1,trim('Minimum_n_year'),g_s%remove_clima%minimum_n_year,1)
  WRITE(1,*)"]"
  WRITE(1,*)","
  
  !!SNHT_parameter
  WRITE(1,*)"""SNHT_Param""", ":"
  WRITE(1,*)"[" 
  status = write_int_to_json_string(1,trim('Increment'),g_s%SNHT_Param%increment,0)
  status = write_int_to_json_string(1,trim('Max_len'),g_s%SNHT_Param%max_len,0)
  status = write_int_to_json_string(1,trim('Mean_maxlen'),g_s%SNHT_Param%mean_maxlen,0)
  status = write_int_to_json_string(1,trim('RAOB_max_miss'),g_s%SNHT_Param%RAOB_max_miss,1)
  WRITE(1,*)"]"
  WRITE(1,*)","
  
  !!BREAK_parameters
  WRITE(1,*)"""Break_Param""", ":"
  WRITE(1,*)"[" 
  !status = write_int_to_json_string(1,trim('Threshold'),g_s%Break_Param%threshold,0)
  status = write_int_to_json_string(1,trim('Threshold_T'),g_s%Break_Param%threshold_T,0)
  status = write_int_to_json_string(1,trim('Threshold_WD'),g_s%Break_Param%threshold_WD,0)
    status = write_int_to_json_string(1,trim('Threshold_WS'),g_s%Break_Param%threshold_WS,0)
  status = write_real_to_json_string(1,trim('Threshold_probability'),g_s%Break_Param%threshold_probability,0)
  status = write_int_to_json_string(1,trim('Threshold_rad'),g_s%Break_Param%threshold_rad,0)
  status = write_real_to_json_string(1,trim('Fak'),g_s%Break_Param%fak,0)
  status = write_real_to_json_string(1,trim('Ap_prob_default_T'),g_s%Break_Param%ap_prob_default_T,0)
  status = write_real_to_json_string(1,trim('Ap_prob_default_WDWS'),g_s%Break_Param%ap_prob_default_WDWS,0)
  !status = write_int_to_json_string(1,trim('Threshold_rad'),g_s%Break_Param%threshold_rad,0)
  status = write_int_to_json_string(1,trim('Locsig'),g_s%Break_Param%locsig,0)
  status = write_int_to_json_string(1,trim('Br_max'),g_s%Break_Param%br_max,0)
  status = write_real_to_json_string(1,trim('Smooth_fit_threshold'),g_s%Break_Param%smooth_fit_threshold,0)
  status = write_int_to_json_string(1,trim('Bg_correction_factor'),g_s%Break_Param%bg_correction_factor,0)
  status = write_int_to_json_string(1,trim('Qc'),g_s%Break_Param%qc,0)
  status = write_int_to_json_string(1,trim('Old'),g_s%Break_Param%old,0)
  status = write_int_to_json_string(1,trim('Use_meta'),g_s%Break_Param%use_meta,0)
  status = write_real_to_json_string(1,trim('Sig_thresh_T'),g_s%Break_Param%sig_thresh_T,0)
  status = write_real_to_json_string(1,trim('Sig_thresh_WS'),g_s%Break_Param%sig_thresh_WS,0)
  status = write_real_to_json_string(1,trim('Sig_thresh_mean_WD'),g_s%Break_Param%sig_thresh_mean_WD,0)
  status = write_real_to_json_string(1,trim('Sig_thresh_rad'),g_s%Break_Param%sig_thresh_rad,0)
  status = write_int_to_json_string(1,trim('Bg_hom'),g_s%Break_Param%bg_hom,0)
  status = write_int_to_json_string(1,trim('Wide_break_detection_time_fac'),g_s%Break_Param%wide_break_detection_time_fac,0)
  status = write_int_to_json_string(1,trim('Fine_break_detection_time_fac'),g_s%Break_Param%fine_break_detection_time_fac,1)
  WRITE(1,*) "]"
  WRITE(1,*)","
  
  !RAOBCORE PARAMETERS
  WRITE(1,*)'"RAOBCORE_param"', ":"
  WRITE(1,*)"[" 
  status = write_int_to_json_string(1,trim('smooth_method'),g_s%RAOBCORE_param%smooth_method,0)
  status = write_int_to_json_string(1,trim('prob_method'),g_s%RAOBCORE_param%prob_method,1)
  WRITE(1,*) "]"
  WRITE(1,*)","
  
       !RICH PARAMETERS
       WRITE(1,*)'"RICH_param"', ":"
       WRITE(1,*)"[" 
       status = write_int_to_json_string(1,trim('RICH_max_miss'),g_s%RICH_Param%RICH_max_miss,0)
       status = write_int_to_json_string(1,trim('Weight_distance'),g_s%RICH_Param%weight_distance,0)
       status = write_int_to_json_string(1,trim('Minst'),g_s%RICH_Param%minst,0)
       status = write_int_to_json_string(1,trim('Ens'),g_s%RICH_Param%ens,1)
       WRITE(1,*) "]"
       WRITE(1,*)","

       !GOOD SONDES
       !WRITE(1,*)'"Good_sondes_dim"',":"
       !WRITE(1,*)"[" 
       !status = write_int_arr_line_to_json_string(1,g_s%good_sondes_type%dim,1)
       !WRITE(1,*) "]"
       !WRITE(1,*)","
       WRITE(1,*)'"Good_sondes"',":"
       WRITE(1,*)"["
       status = write_int_to_json_string(1,trim('dim'),g_s%good_sondes_type%dim,0)
       ALLOCATE(g_s%good_sondes_type%good_sondes(g_s%good_sondes_type%dim), STAT = status)
       mex="ALLOCATE(g_s%good_sondes_type%good_sondes(g_s%good_sondes_type%dim), STAT = status)"
       CALL  error(routine_name,status,2,mex)
       g_s%good_sondes_type%good_sondes=(/37,52,53,60,61,62,63,66,67,71,72,73,74,78,79,80,81,82,83,47,55,56,26,76,85,86,87/)
       do i= 1,g_s%good_sondes_type%dim
          if(i .EQ. g_s%good_sondes_type%dim) then
             status = write_int_arr_line_to_json_string(1,g_s%good_sondes_type%good_sondes(i),1)
          else
             status = write_int_arr_line_to_json_string(1,g_s%good_sondes_type%good_sondes(i),0)
          end if
       end do
        DEALLOCATE(g_s%good_sondes_type%good_sondes, STAT = status)
       mex="DEALLOCATE(g_s%good_sondes_type%good_sondes, STAT = status)"
       CALL  error(routine_name,status,2,mex)
       WRITE(1,*) "]"
       WRITE(1,*)","
       !!PRESSURE_LAYERS_WEIGHTS
       WRITE(1,*)'"Pressure_layers_weights"',":"
       WRITE(1,*)"["
       do i=1,g_s%Dimensions%n_layers
          if(i .EQ. g_s%Dimensions%n_layers) then
             status = write_int_arr_line_to_json_string(1,g_s%pressure_layers_weights(i),1)
          else
             status = write_int_arr_line_to_json_string(1,g_s%pressure_layers_weights(i),0)
          end if
       enddo
       WRITE(1,*) "]"
       WRITE(1,*)","
       !!Pressure_layers_smooth_weights
       WRITE(1,*)'"Pressure_layers_smooth_weights"',":"
       WRITE(1,*)"["
       do i=1,g_s%Dimensions%n_layers
          if(i .EQ. g_s%Dimensions%n_layers) then
             status = write_int_arr_line_to_json_string(1,g_s%pressure_layers_smooth_weights(i),1)
          else
             status = write_int_arr_line_to_json_string(1,g_s%pressure_layers_smooth_weights(i),0)
          end if
       enddo
       WRITE(1,*) "]"
       WRITE(1,*)","
       !Reduction_factor
       WRITE(1,*)'"Reduction_factor"',":"
       WRITE(1,*)"["
       F1='3'
       F2='1'
       do i=1,g_s%Dimensions%n_layers
          if(i .EQ. g_s%Dimensions%n_layers) then
             status = write_real_arr_line_to_json_string(1,g_s%reduction_factor(i),F1,F2,1)
          else
             status = write_real_arr_line_to_json_string(1,g_s%reduction_factor(i),F1,F2,0)
          end if
       enddo
       WRITE(1,*) "]"
       WRITE(1,*)","
       !!Max_break_size
       WRITE(1,*)'"Max_break_size"',":"
       WRITE(1,*)"["
       F1='5'
       F2='1'
       do i=1,size(g_s%varno_in)
          if(i .EQ. size(g_s%varno_in)) then
             status = write_real_arr_line_to_json_string(1,g_s%max_break_size(i),F1,F2,1)
          else
             status = write_real_arr_line_to_json_string(1,g_s%max_break_size(i),F1,F2,0)
          end if
       enddo
       WRITE(1,*) "]"
       WRITE(1,*)","

       !!LAST PRESSURE I AGGIUST to the MOST RECENT
       WRITE(1,*)"""Press_lats_adj_most_recent""", ":"
       WRITE(1,*)"[" 
       status = write_int_arr_line_to_json_string(1, g_s%press_lats_adj_most_recent,1)   
       WRITE(1,*)"]"
       WRITE(1,*)","
       !!TIME SERIES NAME
       WRITE(1,*)'"Time_series_name"',":"
       WRITE(1,*)"["
       do i=1,size(g_s%priority_time_series%ts_name)
          if(i .EQ. size(g_s%priority_time_series%ts_name)) then
             status = write_string_arr_line_to_json_string(1,g_s%priority_time_series%ts_name(i),1)
          else
             status = write_string_arr_line_to_json_string(1,g_s%priority_time_series%ts_name(i),0)
          end if
       enddo
       WRITE(1,*) "]"
       WRITE(1,*)","
       !!TIME SERIES PRIORITY
       WRITE(1,*)'"Time_series_priority"',":"
       WRITE(1,*)"["
       do i=1,size(g_s%priority_time_series%ts_name)
          if(i .EQ. size(g_s%priority_time_series%ts_name)) then
             status = write_int_arr_line_to_json_string(1,g_s%priority_time_series%ts_priority(i),1)
          else
             status = write_int_arr_line_to_json_string(1,g_s%priority_time_series%ts_priority(i),0)
          end if
       enddo
       WRITE(1,*) "]"
       WRITE(1,*)","

       !WIND PARAMETER
       WRITE(1,*)'"Wind_parameters"',':'
       WRITE(1,*)"["
       status = write_int_to_json_string(1,trim('min_WS_for_WD_dep'),g_s%Wind_parameters%min_WS_for_WD_dep,0)
       status = write_int_to_json_string(1,trim('break_det_method_mean_vs_3max'),g_s%Wind_parameters%break_det_method_mean_vs_3max,0)
       status = write_int_to_json_string(1,trim('first_press_lvl'),g_s%Wind_parameters%wind_press_layers(1),0)
       status = write_int_to_json_string(1,trim('last_press_lvl'),g_s%Wind_parameters%wind_press_layers(2),0)
       status = write_int_to_json_string(1,trim('min_num_level_00_12'),g_s%Wind_parameters%min_num_level_00_12,0)
       status = write_real_to_json_string(1,trim('min_mean_WD_profile'),g_s%Wind_parameters%min_mean_WD_profile,0)
       status = write_real_to_json_string(1,trim('st_dev_vertical_profile'),g_s%Wind_parameters%st_dev_vertical_profile,1)
       WRITE(1,*) "]" 
       WRITE(1,*)","
       
       !GENERAL INFO
       !!INNOVATION
       WRITE(1,*)'"Innovation"', ":"
       WRITE(1,*)"[" 
       status = write_string_arr_line_to_json_string(1,g_s%innovation,1)
       WRITE(1,*) "]"
       WRITE(1,*)","
       !!EXTENDED OUTPUT
       WRITE(1,*)'"Extended_output"', ":"
       WRITE(1,*)"[" 
       status = write_string_arr_line_to_json_string(1,g_s%extended_output,1)
       WRITE(1,*)"]"
       WRITE(1,*)"," 
       
       !!MAX_ITER
       WRITE(1,*)"""Max_iter""", ":"
       WRITE(1,*)"[" 
       status = write_int_arr_line_to_json_string(1, g_s%max_iter,1)   
       WRITE(1,*)"]"
       WRITE(1,*)"," 

       !DEBUG MODED
       WRITE(1,*)"""Debug_mode""", ":"
       WRITE(1,*)"[" 
       g_s%debug_mode%dim=45
       status = write_int_to_json_string(1,trim('dim'),g_s%debug_mode%dim,0)
       ALLOCATE(g_s%debug_mode%routine_name(g_s%debug_mode%dim), STAT = status)
       mex="ALLOCATE(g_s%debug_mode%routine_name(g_s%debug_mode%dim), STAT = status)"
       CALL  error(routine_name,status,2,mex)
       i=1
       g_s%debug_mode%routine_name(i) ="rasocorrect_main_new"
       i=i+1
       g_s%debug_mode%routine_name(i) ="raso_correct_new"
       i=i+1
       g_s%debug_mode%routine_name(i) ="XX_apriori_prob_new"
       i=i+1
       g_s%debug_mode%routine_name(i) ="XX_adjust_bg_ei_new"
       i=i+1
       g_s%debug_mode%routine_name(i) ="XX_write_ncdf"
       i=i+1
       g_s%debug_mode%routine_name(i) ="XX_read_ncdf"
       i=i+1
       g_s%debug_mode%routine_name(i) ="XX_write_sonde_montly_netcdf"
       i=i+1
       g_s%debug_mode%routine_name(i) ="XX_write_WDWS_sonde_montly_netcdf"
       i=i+1
       g_s%debug_mode%routine_name(i) ="XX_fill_attributes"
       i=i+1
       g_s%debug_mode%routine_name(i) ="XX_gen_nc"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_apriori_prob_new"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_write_apriori_prob_netcdf"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_gen_array_nc"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_correct_break_new"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_make_montly_new"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_make_cru_an_Tmontly_new"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_anomaly_new"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_calc_mean_new"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_average_rad_new"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_correct_break_T_new"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_correct_break_WD_new"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_correct_break_WS_new"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_bayes_break_simple_new"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_locate_combined_breaks_T_new"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_locate_combined_breaks_WD_new"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_locate_combined_breaks_WS_new"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_select_breaks_T_new"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_select_breaks_WD_new"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_select_breaks_WS_new"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_adjust_series_T_new"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_adjust_series_WDS_new"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_calc_profile_new"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_test_significance_T_new"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_test_significance_WS_new"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_test_significance_WD_new"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_correct_mostrecent_new"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_save_state_netcdf_new"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_corr_stats_new"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_calc_gapbreak_new"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_gen_n_array_nc"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_write_output_nc"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_anomaly_obs_bg_climatology"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_WRITING_STAT_NETCDF"
       i=i+1
       g_s%debug_mode%routine_name(i)="XX_WIND_extract2_ts_netcdf"
        i=i+1
       g_s%debug_mode%routine_name(i)="XX_read_montly_nc"

       do i= 1,g_s%debug_mode%dim
          if(i .EQ. g_s%debug_mode%dim) then
             status = write_string_to_json_string(1,trim(''),trim(g_s%debug_mode%routine_name(i)),1)
             
          else
             status = write_string_to_json_string(1,trim(''),trim(g_s%debug_mode%routine_name(i)),0)
          end if
       end do
       DEALLOCATE(g_s%debug_mode%routine_name, STAT = status)
       mex="DEALLOCATE(g_s%debug_mode%routine_name, STAT = status)"
       CALL  error(routine_name,status,2,mex)
       WRITE(1,*)"]"
       
    

       !a json file end with "}"
       WRITE(1,*)"}"
    else file
       WRITE(*,*) "Error opening file ", filename
    end if file
    
    !closing file
 CLOSE (UNIT = 1,  IOSTAT = status)
 
 
 !!I return status function value
 write_json = status
 
END FUNCTION write_json



!!~~~~~~~~~ FUNCTION TO READ JSON FILES~~~~~~~~~~~~!!
INTEGER FUNCTION read_json(filename,routine_name,g_s_read)
  
  IMPLICIT NONE
  
  TYPE(global_settings_raobcore) :: g_s_read !! global settings variables
    
  CHARACTER(len = *) :: filename !file name
  CHARACTER(len = 100000) :: long_string
  CHARACTER(len = 200) :: temp
  CHARACTER(len = 100) :: cform
  CHARACTER(len = 100) :: aus,case
  
  INTEGER :: status               !opening status
  INTEGER :: status11             !reading file
  INTEGER :: status1              !internal operations
  INTEGER :: n_lines              !number of lines
  INTEGER :: value

  !COUNTERS
  INTEGER :: count,i,j,k,z,l, temp_long, loop_count, loop_incr
  INTEGER :: pos, pos_start, pos_end,pos_end_int   !!these are string position
  INTEGER :: pos1, pos2, pos_1, pos_2
  INTEGER :: pos3, pos4
  INTEGER :: pos5, pos6
  CHARACTER (len = 200) ::  routine_name, error_mex
  !The routine name
  routine_name ="read_json"
  
  WRITE(*,*)""
  error_mex=" working to read the Json file: "//filename(1:len_trim(filename))
  call error(routine_name,0,0,error_mex)
    
  OPEN(UNIT = 1, FILE = filename, STATUS = "OLD", ACTION = "READ", IOSTAT= status)
  !if the open file is well done, I work on
  file:  if (status .EQ. 0) then
     
     !1) Read line by line and put al in a huge string
     count = 1 !string position
     do       
        !here I suppose the lines have less than 100 character
        READ(1,'(A200)', IOSTAT = status11) temp
        i = 0
        long_string(count:count+len_trim(temp)) = trim(temp)
        count = count + len_trim(temp)
        if (status11 .NE. 0)then
           status11 = 0
           EXIT
        end if
     end do
     
     
     i = 1
     count = 0
     status1 = 0
     pos_start = 0 !{ position in the file
     pos_end = 0   !} position in the file
     
     !Now I have the huge string..I analyse it.
     !Here I check if I have a JSON file
     !the JSON file start with "{" and end with "}"
     !I search them
     !write(*,*)len_trim(long_string)
     do while (i < len_trim(long_string) )
        if(long_string(i:i) .EQ. "{") then
           pos_start = i
           EXIT
        end if
        i = i +1
     end do
     if(pos_start .eq.0)then
        error_mex="NO JSON FORMAT -> TERMINATED PROGRAM"
        call error(routine_name,-999,2,error_mex)
          read_json=1
          RETURN
        end if
     i = pos_start
     do while(i < len_trim(long_string) )
        if(long_string(i:i) .EQ. "}") then
           pos_end = i
           EXIT
        end if
        i = i+1
     end do
     if(pos_end .eq.0)then
        error_mex="NO JSON FORMAT -> TERMINATED PROGRAM"
        call error(routine_name,-999,2,error_mex)
        read_json=1
        RETURN
     end if
     
     !check if the file is a Json file
     if(pos_end > pos_start ) then
        error_mex="----> Valid JSON file"  
        call error(routine_name,0,2,error_mex)
     else 

        error_mex= filename//' is NOT a  JSON fil!'
        call error(routine_name,-999,2,error_mex)
        status1 = 1
        STOP
     end if
     
     !now I start with the string analysis
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !if there are some exceptions, loop_count check them
     loop_count = 0
     !pos_end -2 : for sure at the end I have ]} and it does not matter
     do while(pos_start < (pos_end-2))
        !Now I check the """ so I know the argument
        i = 1
        do i = pos_start, pos_end
           if (long_string(i:i) .EQ. '"') then
              pos1 = i
              EXIT
           end if
        end do
        i = 0
        do i = pos1+1, pos_end
           if (long_string(i:i) .EQ. '"') then
              pos2 = i
              EXIT
           end if
        end do
        !!what there is beetwen "  " is my argument
        !!and now I inquire it
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
        !write(*,*)"Item_title ->",long_string(pos1+1:pos2-1)
        SELECT CASE (long_string(pos1+1:pos2-1))
        CASE("title")
           case = long_string(pos1+1:pos2-1)
           !INPUT
           !long_strin: all the json file as a string
           !int_case: the actual case : "title","n_days", "n_layers","n_obs".....
           !pos2: from where the function starts to read ->after the " of the int_case 
           !OUTPUT:
           !pos: return the actual position in the long_string
           !INTEGER FUNCTION inquire_json_string(long_string, int_case, pos2,pos)
           status =  inquire_json_string(long_string, case , pos2 ,pos_end_int,g_s_read)
           status11 = status11 + status
        CASE("stat_info")
           case = long_string(pos1+1:pos2-1)
           status =  inquire_json_string(long_string, case , pos2,pos_end_int,g_s_read)
           status11 = status11 + status
        CASE("Dimensions")
           case = long_string(pos1+1:pos2-1)
           status =  inquire_json_string(long_string, case , pos2,pos_end_int,g_s_read)
           status11 = status11 + status 
        CASE("Time_Window")
           case = long_string(pos1+1:pos2-1)
           status =  inquire_json_string(long_string, case , pos2,pos_end_int,g_s_read)
           status11 = status11 + status       
        CASE("time_serie_initialize")
           Case = Long_string(Pos1+1 : pos2-1)
           status =  inquire_json_string(long_string, case , pos2,pos_end_int,g_s_read)
           status11 = status11 + status
        CASE("Varno_in")
           Case = Long_string(Pos1+1 : pos2-1)
           status =  inquire_json_string(long_string, case , pos2,pos_end_int,g_s_read)
           status11 = status11 + status
        CASE("Var_check")
           Case = Long_string(Pos1+1 : pos2-1)
           status =  inquire_json_string(long_string, case , pos2,pos_end_int,g_s_read)
           status11 = status11 + status 
        CASE("Observation_time")
           Case = Long_string(Pos1+1 : pos2-1)
           status =  inquire_json_string(long_string, case , pos2,pos_end_int,g_s_read)
           status11 = status11 + status  
        CASE("Pressure_layers")
           Case = Long_string(Pos1+1 : pos2-1)
           status =  inquire_json_string(long_string, case , pos2,pos_end_int,g_s_read)
           status11 = status11 + status
      CASE("paths")
           Case = Long_string(Pos1+1 : pos2-1)
           status =  inquire_json_string(long_string, case , pos2,pos_end_int,g_s_read)
           status11 = status11 + status
        CASE("Remove_clima")
           Case = Long_string(Pos1+1 : pos2-1)
           status =  inquire_json_string(long_string, case , pos2,pos_end_int,g_s_read)
           status11 = status11 + status   
       CASE("SNHT_Param")
           Case = Long_string(Pos1+1 : pos2-1)
           status =  inquire_json_string(long_string, case , pos2,pos_end_int,g_s_read)
           status11 = status11 + status     
       CASE("Break_Param")
           Case = Long_string(Pos1+1 : pos2-1)
           status =  inquire_json_string(long_string, case , pos2,pos_end_int,g_s_read)
           status11 = status11 + status
        CASE("RAOBCORE_param")
           Case = Long_string(Pos1+1 : pos2-1)
           status =  inquire_json_string(long_string, case , pos2,pos_end_int,g_s_read)
           status11 = status11 + status
       CASE("RICH_param")
           Case = Long_string(Pos1+1 : pos2-1)
           status =  inquire_json_string(long_string, case , pos2,pos_end_int,g_s_read)
           status11 = status11 + status 
        CASE("Good_sondes_dim")
           Case = Long_string(Pos1+1 : pos2-1)
           status =  inquire_json_string(long_string, case , pos2,pos_end_int,g_s_read)
           status11 = status11 + status
        CASE("Good_sondes")
           Case = Long_string(Pos1+1 : pos2-1)
           status =  inquire_json_string(long_string, case , pos2,pos_end_int,g_s_read)
           status11 = status11 + status
        CASE("Pressure_layers_weights")
           Case = Long_string(Pos1+1 : pos2-1)
           status =  inquire_json_string(long_string, case , pos2,pos_end_int,g_s_read)
           status11 = status11 + status
        CASE("Pressure_layers_smooth_weights")
           Case = Long_string(Pos1+1 : pos2-1)
           status =  inquire_json_string(long_string, case , pos2,pos_end_int,g_s_read)
           status11 = status11 + status
        CASE("Max_break_size")
           Case = Long_string(Pos1+1 : pos2-1)
           status =  inquire_json_string(long_string, case , pos2,pos_end_int,g_s_read)  
           status11 = status11 + status
        CASE("Reduction_factor")
           Case = Long_string(Pos1+1 : pos2-1)
           status =  inquire_json_string(long_string, case , pos2,pos_end_int,g_s_read)  
           status11 = status11 + status 
        CASE("Press_lats_adj_most_recent")
           Case = Long_string(Pos1+1 : pos2-1)
           status =  inquire_json_string(long_string, case , pos2,pos_end_int,g_s_read)  
           status11 = status11 + status   
        CASE("Time_series_name")
           Case = Long_string(Pos1+1 : pos2-1)
           status =  inquire_json_string(long_string, case , pos2,pos_end_int,g_s_read)
           status11 = status11 + status
        CASE("Time_series_priority")
           Case = Long_string(Pos1+1 : pos2-1)
           status =  inquire_json_string(long_string, case , pos2,pos_end_int,g_s_read)
           status11 = status11 + status    
        CASE("Innovation")
           Case = Long_string(Pos1+1 : pos2-1)
           status =  inquire_json_string(long_string, case , pos2,pos_end_int,g_s_read)
           status11 = status11 + status
        CASE("Extended_output")
           Case = Long_string(Pos1+1 : pos2-1)
           status =  inquire_json_string(long_string, case , pos2,pos_end_int,g_s_read)
           status11 = status11 + status
        CASE("Max_iter")
           Case = Long_string(Pos1+1 : pos2-1)
           status =  inquire_json_string(long_string, case , pos2,pos_end_int,g_s_read)
           status11 = status11 + status
        CASE("Wind_parameters")
           Case = Long_string(Pos1+1 : pos2-1)
           status =  inquire_json_string(long_string, case , pos2,pos_end_int,g_s_read)
           status11 = status11 + status
        CASE("Debug_mode")
           Case = Long_string(Pos1+1 : pos2-1)
           status =  inquire_json_string(long_string, case , pos2,pos_end_int,g_s_read)
           status11 = status11 + status
        CASE DEFAULT
           !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!here there is an infinite loop since it doesn't update the pos!!!
!!!!!!!!STOP the program!!!!!!!!
           error_mex="WARINING -> No type < "//long_string(pos1+1:pos2-1)//" > in the global setting"
           call error(routine_name,-999,2,error_mex)
           !STOP
        END SELECT
        
        pos_start = pos_end_int
        loop_incr = pos_start
        if (loop_incr .EQ. pos_end_int )then
           loop_count = loop_count +1
           !WRITE(*,*)"loop_count------->",loop_count
           if (loop_count .EQ. 1)then
              !EXIT
              cycle
           end if
        end if
       ! loop_incr = pos_start
     end do 

     !Here I fill in the date ALLOCATABLE allary
     !The allocation is done in the function
     status =  Date_filling(g_s_read)
  else file
     error_mex= "Error opening file--> "// filename
     call error(routine_name,-999,2,error_mex)
     status1 = 1
  end if file
  
  
  !closing file
  CLOSE (UNIT = 1)
  
  if(status .EQ. -1) then
     status =0 !succesful
  end if
  
  !!I return status function value
  read_json = status + status11+ status1
  
END FUNCTION read_json



!!~~~~~~~~~ FUNCTION INQUIRE JSON STRING~~~~~~~~~~~~!!
!INPUT
!long_strin: all the json file as a string
!int_case: the actual case : "title","n_days", "n_layers","n_obs".....
!pos2: from where the function starts to read ->after the " of the int_case 
!OUTPUT:
!pos: return the actual position in the long_string
INTEGER FUNCTION inquire_json_string(long_string, int_case, pos2,pos,g_s_read) 
  IMPLICIT NONE
  
  CHARACTER(len = 10000) :: long_string
  CHARACTER(len = 30) :: int_case
  INTEGER :: pos2
  
  
  !internal
  TYPE(global_settings_raobcore) :: g_s_read !! global settings variables
  
  CHARACTER(len = 200) :: temp, name,temp_string
  CHARACTER(len = 100) :: cform
  CHARACTER(len = 100) :: aus,aus1
  
  INTEGER :: value, status
  REAL::r_value
  !COUNTERS
  INTEGER :: count,i,j,k,kk,kkk,z,l, temp_long, check_string,ii
  INTEGER :: pos, pos_start=0, pos_end=0,flag_real
  INTEGER :: pos_1=0, pos_2=0
  INTEGER :: pos3=0, pos4=0
  INTEGER :: pos5=0, pos6=0, pos_55=0
  INTEGER :: pos_2points=0,pos_2points_bis=0, first_pos=0,final_pos=0,comma_pos=0
  INTEGER :: dim =0, dim1=0
  INTEGER :: item
  !TEMPORARY
  INTEGER :: temp_int
  REAL:: temp_real
  CHARACTER(len=200)::error_mex,routine_name
  
  routine_name="inquire_json_string"
  status =0
  pos3 = pos2 !because the I need two differnt case 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! the internal extension beetwen []
  j = 0
  do j = pos2,len_trim(long_string)
     if (long_string(j:j) .EQ. '[') then
        pos3 = j
        EXIT
     end if
  end do
  j = 0

  do j = pos3+1,len_trim(long_string)
     if (long_string(j:j) .EQ. ']') then
        pos4 = j
        EXIT
     end if
  end do
  pos = pos3+1
  z = 0
  check_string = 0

  !now I want to count the number of items in the array
  k=0
  !if item =1 => scalar else I have an array
  item=1
  do k = pos3,pos4
     if(long_string(k:k) .EQ. ',') then
        item=item+1
     endif
  end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!pos = start
  !here I check if I have a string array or a integer array
     do while(pos < pos4)
       if(long_string(pos:pos) .EQ. '"')then
          check_string = 1
          EXIT
       end if
       pos = pos +1
    end do
    
    !now I start with the analysis
    pos = pos3+1
    if(check_string .EQ. 1 ) then
    !string analysis
   item_analysis: do while (pos<pos4)
       !In z I have the item counter
       z = z +1
       if(z .gt. item)then
          write(*,*)"WARNING -> probably ITEMS mistake"
       endif
       k =0
       do k = pos3+1,pos4-1
          
          !WRITE(*,*)"k",k
          !WRITE(*,*)"long_string(k:k)",long_string(k:k)
          if (long_string(k:k) .EQ. '"') then
             pos5 = k
             EXIT
          end if
       end do
       k = 0
       do k = pos5+1,pos4-1
          if (long_string(k:k) .EQ. '"') then
             pos6 = k
             EXIT
          end if
       end do
       pos_2points=0
!now I have to find the ":" in the lines
       do k = pos5+1,pos6-1
          if (long_string(k:k) .EQ. ':') then
             pos_2points = k  
             exit
          end if
       end do
       
!If pos_2points=0 -> I have simply a string
       name=""
       SELECT CASE (int_case)
          !this is the only string I have
       CASE("title")
          !here I should have just one :
          status =json_Line_to_INT_REAL_STRING(long_string,pos5,pos6,pos_2points,name,temp_int,temp_real,temp_string)
          !write(*,*)"TITLE->temp_string->",temp_string
          g_s_read%title(z)=temp_string! long_string(pos5+1 : pos6-1)
       !WRITE(*,*)long_string(pos5+1:pos6-1)
       CASE("stat_info")         
          status =json_Line_to_INT_REAL_STRING(long_string,pos5,pos6,pos_2points,name,temp_int,temp_real,temp_string)
          
          !WRITE(*,*)"name -",trim(name),"-"
          if(trim(name) .eq. 'start')then
             !WRITE(*,*)"name -",trim(name),"-"
             !WRITE(*,*)"",temp_int
             g_s_read%Stat_info%start=temp_int
          elseif(trim(name) .eq. 'end')then
             g_s_read%Stat_info%end=temp_int
          elseif(trim(name) .eq. 'number_of_station')then
             g_s_read%Stat_info%number_of_station=temp_int
          elseif(trim(name) .eq. 'available_stations')then
             g_s_read%Stat_info%available_stations=temp_int
          else
             write(*,*)"WARNING -> UNKNOWN NAME->",name
          end if
          
       CASE("Dimensions")
          
          status =json_Line_to_INT_REAL_STRING(long_string,pos5,pos6,pos_2points,name,temp_int,temp_real,temp_string)
          
          !WRITE(*,*)"name -",trim(name),"-"
          if(trim(name) .eq. 'n_days')then
             !WRITE(*,*)"name -",trim(name),"-"
             !WRITE(*,*)"",temp_int
             g_s_read%dimensions%n_days=temp_int
          elseif(trim(name) .eq. 'n_layers')then
             g_s_read%dimensions%n_layers=temp_int
          elseif(trim(name) .eq. 'n_obs')then
             g_s_read%dimensions%n_obs=temp_int
          elseif(trim(name) .eq. 'prob_max')then
             g_s_read%dimensions%prob_max=temp_int
          elseif(trim(name) .eq. 'n_metadata')then
             g_s_read%dimensions%n_metadata=temp_int
          elseif(trim(name) .eq. 'montly_threshold')then
             g_s_read%dimensions%montly_threshold=temp_int  
          elseif(trim(name) .eq. 'max_expected_breaks')then
             g_s_read%dimensions%max_expected_breaks=temp_int
          elseif(trim(name) .eq. 'min_n_days_montly_mean')then
             g_s_read%dimensions%min_n_days_montly_mean=temp_int
          else
             write(*,*)"WARNING -> UNKNOWN NAME->",name
          end if
!!$          
       CASE("Time_Window")
          status =json_Line_to_INT_REAL_STRING(long_string,pos5,pos6,pos_2points,name,temp_int,temp_real,temp_string)
          !WRITE(*,*)"name -",trim(name),"-"
          if(trim(name) .eq. 'start_day')then
             !WRITE(*,*)"name -",trim(name),"-"
             !WRITE(*,*)"",temp_int
             g_s_read%time_window%start_day=temp_int
          elseif(trim(name) .eq. 'end_day')then
             g_s_read%time_window%end_day=temp_int
          elseif(trim(name) .eq. 'first_month')then
             g_s_read%time_window%first_month=temp_int
          elseif(trim(name) .eq. 'last_month')then
             g_s_read%time_window%last_month=temp_int
          elseif(trim(name) .eq. 'swich_CR')then
             g_s_read%time_window%swich_archive%tw_cr=temp_int
          elseif(trim(name) .eq. 'swich_ERA_FT')then
             g_s_read%time_window%swich_archive%era_ft=temp_int
          elseif(trim(name) .eq. 'swich_ERA_INT')then
             g_s_read%time_window%swich_archive%era_int=temp_int
          elseif(trim(name) .eq. 'Sat_time_window_start_day')then
             g_s_read%time_window%sat_time_window%start_day=temp_int
          elseif(trim(name) .eq. 'Sat_time_window_end_day')then
             g_s_read%time_window%sat_time_window%end_day=temp_int
          elseif(trim(name) .eq. 'NOAA_4_start_day')then
             g_s_read%Time_window%NOAA_4%start_day=temp_int
          elseif(trim(name) .eq. 'NOAA_4_end_day')then
             g_s_read%Time_window%NOAA_4%end_day=temp_int
          elseif(trim(name) .eq. 'archive')then
             g_s_read%time_window%archive=temp_int
          else
             write(*,*)"WARNING -> UNKNOWN NAME->",name
          end if
          
       CASE("time_serie_initialize")
          !write(*,*)"name -",long_string(pos5+1 : pos6-1),"-"
          status =json_Line_to_INT_REAL_STRING(long_string,pos5,pos6,pos_2points,name,temp_int,temp_real,temp_string)
          g_s_read%ts_initialize=temp_int
          !write(*,*)"temp_int",temp_int
       CASE("Varno_in")
          !write(*,*)"name -",long_string(pos5+1 : pos6-1),"-"
          status =json_Line_to_INT_REAL_STRING(long_string,pos5,pos6,pos_2points,name,temp_int,temp_real,temp_string)
          g_s_read%varno_in(z)=temp_int
          !write(*,*)"temp_int",temp_int 
       CASE("Var_check")
          !write(*,*)"name -",long_string(pos5+1 : pos6-1),"-"
          status =json_Line_to_INT_REAL_STRING(long_string,pos5,pos6,pos_2points,name,temp_int,temp_real,temp_string)
          g_s_read%var_check(z)=temp_int
          !write(*,*)"temp_int",temp_int
       CASE("Observation_time")
          status =json_Line_to_INT_REAL_STRING(long_string,pos5,pos6,pos_2points,name,temp_int,temp_real,temp_string)
          g_s_read%time_obs(z)=temp_int
       CASE("Pressure_layers")
          status =json_Line_to_INT_REAL_STRING(long_string,pos5,pos6,pos_2points,name,temp_int,temp_real,temp_string)
          g_s_read%pressure_layers(z)=temp_int 
       CASE("paths")
          status =json_Line_to_INT_REAL_STRING(long_string,pos5,pos6,pos_2points,name,temp_int,temp_real,temp_string)
          !write(*,*)"name -",trim(temp_string),"-"
          if(trim(name) .eq. "MERGED")then
             g_s_read%path%merged=temp_string
          elseif(name .eq. "ERA_interim")then
             g_s_read%path%era_interim=temp_string
          elseif(trim(name) .eq. "ERA_forty")then
             g_s_read%path%era_forty=temp_string  
          elseif(trim(name) .eq. "ERA_forty_analysis")then
             g_s_read%path%era_forty_analysis=temp_string  
          elseif(trim(name) .eq. "CHUAN")then
             g_s_read%path%chuan=temp_string  
          elseif(trim(name) .eq. "20CR")then
             g_s_read%path%tw_cr=temp_string 
          elseif(trim(name) .eq. "OUTPUT")then
             g_s_read%path%output=temp_string
          else
             write(*,*)"WARNING -> UNKNOWN NAME->",name
          end if
          
       CASE("Remove_clima")
          status =json_Line_to_INT_REAL_STRING(long_string,pos5,pos6,pos_2points,name,temp_int,temp_real,temp_string)
          if(trim(name) .eq. 'Swich')then
             g_s_read%remove_clima%swich=temp_int
          elseif(trim(name) .eq. 'Minimum_n_year')then
             g_s_read%remove_clima%minimum_n_year=temp_int
          else
             write(*,*)"WARNING -> UNKNOWN NAME->",name
          endif
          
       CASE("SNHT_Param")
          status =json_Line_to_INT_REAL_STRING(long_string,pos5,pos6,pos_2points,name,temp_int,temp_real,temp_string)
          
          !WRITE(*,*)"name -",trim(name),"-"
          if(trim(name) .eq. 'Increment')then
             g_s_read%snht_param%increment=temp_int
          elseif(trim(name) .eq. 'Max_len')then
             g_s_read%snht_param%max_len=temp_int
          elseif(trim(name) .eq. 'Mean_maxlen')then
             g_s_read%snht_param%mean_maxlen=temp_int
          elseif(trim(name) .eq. 'RAOB_max_miss')then
             g_s_read%snht_param%raob_max_miss=temp_int
          else
             write(*,*)"WARNING -> UNKNOWN NAME->",name
          endif
          
       CASE("Break_Param")
          status =json_Line_to_INT_REAL_STRING(long_string,pos5,pos6,pos_2points,name,temp_int,temp_real,temp_string)
          !WRITE(*,*)"name -",trim(name),"-"
          !if(trim(name) .eq.'Threshold')then
          !   !write(*,*)"g_s_read%break_param%threshold_probability",temp_int
          !   g_s_read%break_param%threshold=temp_int
          if(trim(name) .eq.'Threshold_T')then
             !write(*,*)"g_s_read%break_param%threshold_probability",temp_int
             g_s_read%break_param%threshold_T=temp_int
          elseif(trim(name) .eq.'Threshold_WD')then
             !write(*,*)"g_s_read%break_param%threshold_probability",temp_int
             g_s_read%break_param%threshold_WD=temp_int
          elseif(trim(name) .eq.'Threshold_WS')then
             !write(*,*)"g_s_read%break_param%threshold_probability",temp_int
             g_s_read%break_param%threshold_WS=temp_int
          elseif(trim(name).eq. 'Threshold_probability')then
             !write(*,*)"g_s_read%break_param%threshold_probability",temp_real
             g_s_read%break_param%threshold_probability=temp_real 
          elseif(trim(name) .eq. 'Threshold_rad')then
             !write(*,*)"g_s_read%break_param%threshold_rad",temp_int
             g_s_read%break_param%threshold_rad=temp_int
          elseif(trim(name) .eq. 'Fak')then
             !write(*,*)"g_s_read%break_param%fak",temp_real
             g_s_read%break_param%fak=temp_real
          elseif(trim(name) .eq. 'Ap_prob_default_T')then
             !write(*,*)"g_s_read%break_param%ap_prob_default_T",temp_real
             g_s_read%break_param%ap_prob_default_T=temp_real
          elseif(trim(name) .eq. 'Ap_prob_default_WDWS')then
             !write(*,*)"g_s_read%break_param%ap_prob_default_WDWS",temp_real
             g_s_read%break_param%ap_prob_default_WDWS=temp_real
          elseif(trim(name) .eq. 'Locsig')then
             !write(*,*)"g_s_read%break_param%locsig",temp_int
             g_s_read%break_param%locsig=temp_int
          elseif(trim(name) .eq. 'Br_max')then
             !write(*,*)"g_s_read%break_param%br_max",temp_int
             g_s_read%break_param%br_max=temp_int
          elseif(trim(name) .eq. 'Smooth_fit_threshold')then
             !write(*,*)"g_s_read%break_param%smooth_fit_threshold",temp_real
             g_s_read%break_param%smooth_fit_threshold=temp_real
          elseif(trim(name) .eq. 'Bg_correction_factor')then
             !write(*,*)"Bg_correction_factor",temp_int
             g_s_read%break_param%bg_correction_factor=temp_int 
          elseif(trim(name) .eq. 'Qc')then
             !write(*,*)"Qc",temp_int
             g_s_read%break_param%qc=temp_int
          elseif(trim(name) .eq. 'Old')then
             !write(*,*)"Old",temp_int
             g_s_read%break_param%old=temp_int
          elseif(trim(name) .eq. 'Use_meta')then
             !write(*,*)"Use_meta",temp_int
             g_s_read%break_param%use_meta=temp_int
          elseif(trim(name) .eq. 'Sig_thresh_T')then
             !write(*,*)"Sig_tresh_T",temp_real
             g_s_read%break_param%sig_thresh_T=temp_real
          elseif(trim(name) .eq. 'Sig_thresh_WS')then
             !write(*,*)"Sig_tresh_WS",temp_real
          elseif(trim(name) .eq. 'Sig_thresh_mean_WD')then
             !write(*,*)"Sig_tresh_mean_WD",temp_real  
             g_s_read%break_param%sig_thresh_mean_WD=temp_real
          elseif(trim(name) .eq. 'Sig_thresh_rad')then
             !write(*,*)"Sig_thresh_rad",temp_real
             g_s_read%break_param%sig_thresh_rad=temp_real
          elseif(trim(name) .eq. 'Bg_hom')then
             !write(*,*)"Bg_hom",temp_int
             g_s_read%break_param%bg_hom=temp_int
          elseif(trim(name) .eq. 'Wide_break_detection_time_fac')then
             !write(*,*)"Wide_break_detection_time_fac",temp_int
             g_s_read%break_param%wide_break_detection_time_fac=temp_int
          elseif(trim(name) .eq. 'Fine_break_detection_time_fac')then
             !write(*,*)"Fine_break_detection_time_fac",temp_int
             g_s_read%break_param%fine_break_detection_time_fac=temp_int 
          else
             write(*,*)"WARNING -> UNKNOWN NAME->",name
          endif
       CASE("RAOBCORE_param")
          status =json_Line_to_INT_REAL_STRING(long_string,pos5,pos6,pos_2points,name,temp_int,temp_real,temp_string)
          !WRITE(*,*)"name -",trim(name),"-"
          if(trim(name) .eq. 'smooth_method')then 
             !write(*,*)"g_s_read%raobcore_param%smooth_method",temp_int
             g_s_read%raobcore_param%smooth_method=temp_int
          elseif(trim(name) .eq. 'prob_method')then
             !write(*,*)"g_s_read%raobcore_param%prob_method",temp_int
             g_s_read%raobcore_param%prob_method=temp_int
          else
             write(*,*)"WARNING -> UNKNOWN NAME->",name
          endif
       CASE("RICH_param")
          status =json_Line_to_INT_REAL_STRING(long_string,pos5,pos6,pos_2points,name,temp_int,temp_real,temp_string)
          !WRITE(*,*)"name -",trim(name),"-"
          if(trim(name) .eq. 'RICH_max_miss')then
             !write(*,*)"g_s_read%rich_param%rich_max_miss",temp_int
             g_s_read%rich_param%rich_max_miss=temp_int
          elseif(trim(name) .eq. 'Weight_distance')then
             !write(*,*)"g_s_read%rich_param%Weight_distance",temp_int
             g_s_read%rich_param%Weight_distance=temp_int 
          elseif(trim(name) .eq. 'Minst')then
             !write(*,*)"g_s_read%rich_param%mints",temp_int
             g_s_read%rich_param%minst=temp_int
          elseif(trim(name) .eq. 'Ens')then
             !write(*,*)"g_s_read%rich_param%ens",temp_int
             g_s_read%rich_param%ens=temp_int
          else
             write(*,*)"WARNING -> UNKNOWN NAME->",name
          endif
       CASE("Good_sondes_dim")
          status =json_Line_to_INT_REAL_STRING(long_string,pos5,pos6,pos_2points,name,temp_int,temp_real,temp_string)
          !WRITE(*,*)"Innovation -",trim(temp_string),"-"
          g_s_read%good_sondes_type%dim=temp_int

       CASE("Good_sondes")
          status =json_Line_to_INT_REAL_STRING(long_string,pos5,pos6,pos_2points,name,temp_int,temp_real,temp_string)
          if(trim(name) .eq. 'dim')then
             g_s_read%good_sondes_type%dim=temp_int
          else
             !if I'm at the first iteration->I Have to allocate the array
             if(.not. allocated(g_s_read%good_sondes_type%good_sondes))then
                ALLOCATE(g_s_read%good_sondes_type%good_sondes(g_s_read%good_sondes_type%dim), STAT= status)
                error_mex = "ALLOCATE(g_s_read%good_sondes_type%good_sondes(g_s_read%good_sondes_type%dim), STAT= status)"
                call error(routine_name,status,2,error_mex)
                g_s_read%good_sondes_type%good_sondes=g_s_read%ts_initialize
             endif
             num_wr: do i=1,g_s_read%good_sondes_type%dim
                if( g_s_read%good_sondes_type%good_sondes(i) == g_s_read%ts_initialize)then
                   g_s_read%good_sondes_type%good_sondes(i)=temp_int
                   exit num_wr
                endif
             enddo num_wr
          endif
       CASE("Pressure_layers_weights")
          status =json_Line_to_INT_REAL_STRING(long_string,pos5,pos6,pos_2points,name,temp_int,temp_real,temp_string)
          g_s_read%pressure_layers_weights(z)=temp_int

       CASE("Pressure_layers_smooth_weights")
          status =json_Line_to_INT_REAL_STRING(long_string,pos5,pos6,pos_2points,name,temp_int,temp_real,temp_string)
          g_s_read%pressure_layers_smooth_weights(z)=temp_int 
       CASE("Max_break_size")
          status =json_Line_to_INT_REAL_STRING(long_string,pos5,pos6,pos_2points,name,temp_int,temp_real,temp_string)
          g_s_read%max_break_size(z)=temp_real
       CASE("Reduction_factor") 
          status =json_Line_to_INT_REAL_STRING(long_string,pos5,pos6,pos_2points,name,temp_int,temp_real,temp_string)
          g_s_read%reduction_factor(z)=temp_real
       CASE("Press_lats_adj_most_recent")
            status =json_Line_to_INT_REAL_STRING(long_string,pos5,pos6,pos_2points,name,temp_int,temp_real,temp_string)
          g_s_read%press_lats_adj_most_recent=temp_int
       CASE("Time_series_name")
          status =json_Line_to_INT_REAL_STRING(long_string,pos5,pos6,pos_2points,name,temp_int,temp_real,temp_string)
          g_s_read%priority_time_series%ts_name(z)=temp_string
       CASE("Time_series_priority")
          status =json_Line_to_INT_REAL_STRING(long_string,pos5,pos6,pos_2points,name,temp_int,temp_real,temp_string)
          g_s_read%priority_time_series%ts_priority(z)=temp_int
       CASE("Innovation")
          status =json_Line_to_INT_REAL_STRING(long_string,pos5,pos6,pos_2points,name,temp_int,temp_real,temp_string)
          !WRITE(*,*)"Innovation -",trim(temp_string),"-"
          g_s_read%innovation=temp_string
          !g_s_read%innovation= long_string(pos5+1 : pos6-1)       
       CASE("Extended_output")
          status =json_Line_to_INT_REAL_STRING(long_string,pos5,pos6,pos_2points,name,temp_int,temp_real,temp_string)
          !WRITE(*,*)"Extended_output -",trim(temp_string),"-"
          g_s_read%extended_output=temp_string
       CASE("Max_iter")
          status =json_Line_to_INT_REAL_STRING(long_string,pos5,pos6,pos_2points,name,temp_int,temp_real,temp_string)
          !WRITE(*,*)"Max_iter -",temp_int,"-"
          g_s_read%max_iter=temp_int
       CASE("Wind_parameters")
          
          status =json_Line_to_INT_REAL_STRING(long_string,pos5,pos6,pos_2points,name,temp_int,temp_real,temp_string)

          !WRITE(*,*)"name -",trim(name),"-"
          if(trim(name) .eq. 'min_WS_for_WD_dep')then
             g_s_read%Wind_parameters%min_WS_for_WD_dep=temp_int
          elseif(trim(name) .eq. 'break_det_method_mean_vs_3max')then
             g_s_read%Wind_parameters%break_det_method_mean_vs_3max=temp_int
          elseif(trim(name) .eq. 'first_press_lvl')then
             g_s_read%Wind_parameters%wind_press_layers(1)=temp_int
          elseif(trim(name) .eq. 'last_press_lvl')then
             g_s_read%Wind_parameters%wind_press_layers(2)=temp_int
          elseif(trim(name) .eq. 'min_num_level_00_12')then
             g_s_read%Wind_parameters%min_num_level_00_12=temp_int
             elseif(trim(name) .eq. 'min_mean_WD_profile')then
                g_s_read%Wind_parameters%min_mean_WD_profile=temp_real
          elseif(trim(name) .eq. 'st_dev_vertical_profile')then
             g_s_read%Wind_parameters%st_dev_vertical_profile=temp_real
          else
             write(*,*)"WARNING -> UNKNOWN NAME->",name
          end if
       CASE("Debug_mode")
          !status =json_Line_to_INT_REAL_STRING(long_string,pos5,pos6,pos_2points,name,temp_int,temp_real,temp_string)
          status =json_Line_to_INT_REAL_STRING(long_string,pos5,pos6,pos_2points,name,temp_int,temp_real,temp_string)
          if(trim(name) .eq. 'dim')then
             g_s_read%debug_mode%dim=temp_int
          else
             if(.not. allocated(g_s_read%debug_mode%routine_name))then
                ALLOCATE(g_s_read%debug_mode%routine_name(g_s_read%debug_mode%dim))
                error_mex="ALLOCATE(g_s_read%debug_mode%routine_name(g_s_read%debug_mode%dim))"
                call error(routine_name,status,2,error_mex)
                g_s_read%debug_mode%routine_name=""
             endif
             string_wr: do i=1,g_s_read%debug_mode%dim
                if(g_s_read%debug_mode%routine_name(i) == "")then
                   g_s_read%debug_mode%routine_name(i)=temp_string
                   exit string_wr
                endif
             enddo string_wr
          endif
       CASE DEFAULT
          WRITE(*,*) "Unknown case -->",int_case 
          status = 1
       END SELECT
       
       if(long_string(pos6+1:pos6+1).ne. ',' )then
          !WRITE(*,*)"END of ARRAY"
          exit
       else
          !WRITE(*,*)"I'm in the array body"
          pos3 = pos6+1
          pos = pos6+1
       endif
       
    end do item_analysis
    
 end if
 pos_start = pos4+1
 pos=pos4+1

!!I return status function value
    inquire_json_string  = status 

END FUNCTION inquire_json_string


INTEGER FUNCTION write_string_to_json_string(unit,title_str,body_str,end)
  
  INTEGER :: unit,end !end = 0 no end, else I am at the end
  CHARACTER(len = *)::title_str,body_str
  INTEGER :: len1,len2
  CHARACTER(len = 100) :: str1,str2, cform


  len1= 1 + len_trim(title_str)
  len2= 1 + len_trim(body_str)
  WRITE(str1, '(I2)') len1 
  WRITE(str2, '(I2)') len2
  
  if(end .eq. 0)then
     if(len1 == 1 ) then
        cform = "(A10,A"//trim(str2)//",A1,A1)"
        WRITE(unit,trim(cform))'"',trim(body_str),'"',','
     else
        cform = "(A10,A"//trim(str1)//",A2,A"//trim(str2)//",A1,A1)"
        WRITE(unit,trim(cform))'"',trim(title_str),': ', trim(body_str),'"',','
     endif
  else
     if(len1==1)then
        cform = "(A10,A"//trim(str2)//",A1)"
        WRITE(unit,trim(cform))"""",trim(body_str),""""
     else
        cform = "(A10,A"//trim(str1)//",A2,A"//trim(str2)//",A1)"
        WRITE(unit,trim(cform))"""",trim(title_str),": ", trim(body_str),""""
     endif
  endif
write_string_to_json_string=0

END FUNCTION write_string_to_json_string

INTEGER FUNCTION write_int_to_json_string(unit,title_str,body_int,end)
  INTEGER :: unit
  CHARACTER (len = *)::title_str
  INTEGER :: len1,len2,body_int
 CHARACTER(len =200) :: str1,str2, cform,body_str
  INTEGER ::end !end = 0 no end, else I am at the end
 
  len1= 1+len_trim(title_str)
  write(body_str,'(I10)')body_int
  len2=1+len_trim(body_str)
  WRITE(str1, '(I2)') len1
  WRITE(str2, '(I2)') len2
  
  if(end .eq. 0)then
  !cform = "(A10,A"//trim(str1)//",A2,A"//trim(str2)//",A1,A1)"
  cform = "(A1,A30,A2,A"//trim(str2)//",A1,A1)"
  WRITE(unit,trim(cform))'"',trim(title_str),': ', trim(body_str),'"',','
else
   !cform = "(A10,A"//trim(str1)//",A2,A"//trim(str2)//",A1)"
   cform = "(A1,A30,A2,A"//trim(str2)//",A1)"
   WRITE(unit,trim(cform))'"',trim(title_str),': ', trim(body_str),'"'
endif
write_int_to_json_string=0

END FUNCTION write_int_to_json_string

INTEGER FUNCTION write_int_arr_line_to_json_string(unit,body_int,end)
  INTEGER :: unit
  INTEGER :: len1,len2,body_int
  CHARACTER(len =200) :: str1,str2, cform,body_str
  INTEGER ::end !end = 0 no end, else I am at the end
 
  write(body_str,'(I10)')body_int
  len2=1+len_trim(body_str)
 
  WRITE(str2, '(I2)') len2
  
  if(end .eq. 0)then
  !cform = "(A10,A"//trim(str1)//",A2,A"//trim(str2)//",A1,A1)"
  cform = "(A1,A5,A"//trim(str2)//",A1,A1)"
  WRITE(unit,trim(cform))'"',trim(''), trim(body_str),'"',','
else
   !cform = "(A10,A"//trim(str1)//",A2,A"//trim(str2)//",A1)"
   cform = "(A1,A5,A"//trim(str2)//",A1)"
   WRITE(unit,trim(cform))'"',trim(''), trim(body_str),'"'
endif
write_int_arr_line_to_json_string=0

END FUNCTION write_int_arr_line_to_json_string

INTEGER FUNCTION write_real_arr_line_to_json_string(unit,body_real,F1,F2,end)
  INTEGER :: unit
  CHARACTER(len=2) :: F1, F2
  INTEGER :: len1,len2
  REAL :: body_real
  CHARACTER(len =200) :: str1,str2, cform,body_str
  INTEGER ::end !end = 0 no end, else I am at the end
 
  write(body_str,'(F'//trim(F1)//'.'//trim(F2)//')')body_real
  len2=1+len_trim(body_str)
 
  WRITE(str2, '(I2)') len2
  
  if(end .eq. 0)then
  cform = "(A1,A5,A"//trim(str2)//",A1,A1)"
  WRITE(unit,trim(cform))'"',trim(''), trim(body_str),'"',','
else     
   cform = "(A1,A5,A"//trim(str2)//",A1)"
   WRITE(unit,trim(cform))'"',trim(''), trim(body_str),'"'
endif
write_real_arr_line_to_json_string=0

END FUNCTION write_real_arr_line_to_json_string


INTEGER FUNCTION write_string_arr_line_to_json_string(unit,body_str,end)
  INTEGER :: unit 
  INTEGER :: len1,len2
  CHARACTER(len =*) :: body_str
  CHARACTER(len =200) :: str1,str2, cform
  INTEGER ::end !end = 0 no end, else I am at the end


  len2= 1 + len_trim(body_str)
  WRITE(str2, '(I2)') len2
  
  if(end .eq. 0)then
  !cform = "(A10"//trim(str2)//",A1,A1)"
     cform = "(A1,A"//trim(str2)//",A1,A1)"
  WRITE(unit,trim(cform))'"', trim(body_str),'"',','
else
   cform = "(A1,A"//trim(str2)//",A1)"
   !cform = "(A3,A1)"!"(A"//trim(str2)//",A1)"

   !write(*,"(A1,A3,A1)")"""", trim(body_str),""""
   WRITE(unit,trim(cform))"""", trim(body_str),""""
   
endif

write_string_arr_line_to_json_string=0

END FUNCTION write_string_arr_line_to_json_string


INTEGER FUNCTION write_real_to_json_string(unit,title_str,body_real,end)
  INTEGER :: unit
  CHARACTER (len = *)::title_str
  INTEGER :: len1,len2
  CHARACTER(len =200) :: str1,str2, cform,body_str
  INTEGER ::end !end = 0 no end, else I am at the end
  REAL :: body_real

  len1= 5+len_trim(title_str)
  write(body_str,'(F10.4)')body_real
  len2=5+len_trim(body_str)
  WRITE(str1, '(I2)') len1
  WRITE(str2, '(I2)') len2
  
  if(end .eq. 0)then
  !cform = "(A10,A"//trim(str1)//",A2,A"//trim(str2)//",A1,A1)"
   cform = "(A1,A30,A2,A"//trim(str2)//",A1,A1)"
  WRITE(unit,trim(cform))'"',trim(title_str),': ', trim(body_str),'"',','
else
   !cform = "(A10,A"//trim(str1)//",A2,A"//trim(str2)//",A1)"
   cform = "(A1,A30,A2,A"//trim(str2)//",A1)"
   WRITE(unit,trim(cform))'"',trim(title_str),': ', trim(body_str),'"'
endif
write_real_to_json_string=0

END FUNCTION write_real_to_json_string

INTEGER FUNCTION json_Line_to_INT_REAL_STRING(long_string,pos5,pos6,pos_2points,name,temp_int,temp_real, temp_string)

CHARACTER (len = *) :: long_string
INTEGER :: pos5
INTEGER :: pos6
INTEGER :: pos_2points
INTEGER :: temp_int
REAL :: temp_real
CHARACTER (len =*) ::name, temp_string

INTEGER :: first_pos=0,final_pos=0,comma_pos=0
INTEGER ::dim,dim1
CHARACTER (len =200)::aus,aus1,cform

INTEGER ::k,kk
  

temp_int =-99
temp_real =-99
temp_string=''
!before of all I Have to understand if I'm working with one item without name or with name
points: if(pos_2points .eq.0)then
   !If I'm here I have only to read the items...
   !but It could be a string or a number
   k=0
            !extract_json_line_content(long_string,init_pos,end_pos,temp_int,temp_real,temp_string)
   status = extract_json_line_content(long_string,pos5,pos6,temp_int,temp_real,temp_string)

   !do k=pos5+1,pos6-1
      
   !enddo
else
 !I have the items with the name separated from the value with a :
!Now I have to understand the name
!pos5+1 oin order to avoid the "
k=0
do k = pos5+1,pos_2points-1
   if(long_string(k:k).ne. ' ')then
      pos5=k
      EXIT
   endif
enddo
name = trim(long_string(pos5:pos_2points-1))


!!!In order to estract the vaalues from the json long_string between the init_pos and the end_pos in the correct data type: INTEGER, REAL and STRING
            !extract_json_line_content(long_string,init_pos,end_pos,temp_int,temp_real,temp_string)
   status = extract_json_line_content(long_string,pos_2points,pos6,temp_int,temp_real,temp_string)



endif points
json_Line_to_INT_REAL_STRING=0
END FUNCTION json_Line_to_INT_REAL_STRING


INTEGER FUNCTION extract_json_line_content(long_string,init_pos,end_pos,temp_int,temp_real,temp_string)

CHARACTER (len = *) :: long_string
INTEGER :: pos5
INTEGER :: pos6
INTEGER :: pos_2points
INTEGER :: temp_int
REAL :: temp_real
CHARACTER (len =*) :: temp_string

INTEGER :: first_pos=0,final_pos=0,comma_pos=0
INTEGER ::dim,dim1
CHARACTER (len =200)::aus,aus1,cform

INTEGER ::k,kk

INTEGER :: init_pos
INTEGER :: end_pos
INTEGER :: sign_check=0

!!$!now I have to read the really info stored as number or string
!!$!I initialize the position for the text message

first_pos = 0
final_pos =0
 sign_check=0
k=0
do_char_num:  do k =init_pos+1,end_pos-1
   !write(*,*)"long_string(k:k)",long_string(k:k)
   !write(*,*)"ACHAR(long_string(k:k))",IACHAR(long_string(k:k))
   !now I have to understand if I have a positive or negative number
   !the sign use 1 digits
  
   sign:if(long_string(k:k) .eq. '-')then
      sign_check=1
   endif sign
   if_char_num: if (long_string(k:k) .NE. ' '  .and. (IACHAR(long_string(k:k)).le. 57 .and. IACHAR(long_string(k:k)).ge. 48)) then
      !I'm sure I have a number
      first_pos=k
      do kk =first_pos,end_pos-1
         if (long_string(kk:kk) .NE. ' ')then
            final_pos = kk
         end if
      end do
      
    !  write(*,*)"Number"
    !  write(*,*)"ACHAR(long_string(k:k))",IACHAR(long_string(k:k))
     ! write(*,*)"value =", long_string(first_pos:final_pos)
     
      !now If between first_pos and final_pos
      !I have to search for the . or , => REAL NUMBER else INTEGER NUMBER
      
      comma_pos=0
      comma: do kk= first_pos,final_pos
         !write(*,*)"long_string(kk:kk)",long_string(kk:kk)
         if (long_string(kk:kk) .eq. '.' .or. long_string(kk:kk) .eq. ',')then
            comma_pos = kk
            EXIT comma
         end if
      enddo comma
      !finally the format
      if(comma_pos .eq. 0) then
         !INTEGER
         dim = final_pos-first_pos+1+sign_check
         write(aus,'(I1)')dim
         cform = '(I'//trim(aus)//')'
         !write(*,*)"cform",trim(cform)
         !the conversion string to integer
         !here my INTEGER number
         READ(long_string(first_pos-sign_check:final_pos),'(I3)')temp_int
         READ(long_string(first_pos-sign_check:final_pos),trim(cform))temp_int
         !write(*,*)"temp_int",temp_int
         
      else
         !REAL
         dim = final_pos-first_pos+1+sign_check+1! the +1 is for the comma position
         write(aus,'(I1)')dim
         dim1= final_pos-comma_pos
         !WRITE( aus1,'(I10)') final_pos-comma_pos
         write(aus1,'(I1)')dim1
         cform = '(F'//trim(aus)//"."//trim(aus1)//")'"
         !write(*,*)"cform-",trim(cform),'-'
         READ(long_string(first_pos:final_pos),trim(cform))temp_real
         !write(*,trim(cform))temp_real
         
      endif
      
      EXIT do_char_num
      
      
   elseif (long_string(k:k) .NE. ' '  .and. ( ( IACHAR(long_string(k:k)).ge. 65 .and. IACHAR(long_string(k:k)).le. 90).or.(IACHAR(long_string(k:k)).ge. 97 .and. IACHAR(long_string(k:k)).le. 122)  .or. (IACHAR(long_string(k:k)).eq.47) )) then
      !I'm sure I have a character 
      !write(*,*)"string"
      !write(*,*)"ACHAR(long_string(k:k))",IACHAR(long_string(k:k))
      first_pos=k
      do kk =first_pos+1,end_pos!-1
         if (long_string(kk:kk) .NE. ' ')then
            final_pos = kk
         end if
      end do
      !now all what I have between first and final position
      !I have to store in the string variable
      
      temp_string= long_string(first_pos:final_pos-1)
      !write(*,*)temp_string
      EXIT do_char_num
   endif if_char_num
end do do_char_num




  extract_json_line_content=0
END FUNCTION extract_json_line_content

INTEGER FUNCTION Date_filling(g_s)
  IMPLICIT NONE
  TYPE(global_settings_raobcore) :: g_s !! global settings variables
  INTEGER :: status          !!status of the function
  INTEGER, ALLOCATABLE :: year(:),month(:),day(:),time(:)
  INTEGER :: ini_day,ini_month,ini_year,index,tdiff
  INTEGER :: year1,month1,day1,time1,mod,month_save,month_index
  CHARACTER (len = 100):: routine_name
  CHARACTER (len = 200):: error_mex,toindex_name

  !the routine_name
  routine_name = "Date_filling"
  !!I set the status= 0
  status = 0
  
  ini_year=g_s%time_window%start_day/10000
  ini_month=(g_s%time_window%start_day-ini_year*10000)/100
  ini_day=mod(g_s%time_window%start_day,100)
  
  AllOCATE(year(g_s%dimensions%n_days),month(g_s%dimensions%n_days),day(g_s%dimensions%n_days),time(g_s%dimensions%n_days), STAT = status)
  if(status ==0)then
     error_mex=""
  else
     error_mex = "  AllOCATE(year(g_s%dimensions%n_days),month(g_s%dimensions%n_days),day(g_s%dimensions%n_days),time(g_s%dimensions%n_days), STAT = status)"
  endif
  call error(routine_name,status,2,error_mex)

  allocate(g_s%time_window%date(g_s%dimensions%n_days),g_s%time_window%year(g_s%dimensions%n_days),g_s%time_window%month(g_s%dimensions%n_days),g_s%time_window%day(g_s%dimensions%n_days),g_s%time_window%n_day_x_month(g_s%time_window%last_month) ,STAT = status)
  error_mex = " allocate(g_s%time_window%date(g_s%dimensions%n_days),g_s%time_window%year(g_s%dimensions%n_days),g_s%time_window%month(g_s%dimensions%n_days),g_s%time_window%day(g_s%dimensions%n_days),g_s%time_window%n_day_x_month(g_s%time_window%last_month),STAT = status)"
  call error(routine_name,status,2,error_mex)
  
  !INITIALIZE
  g_s%time_window%date=g_s%ts_initialize
  g_s%time_window%year=g_s%ts_initialize
  g_s%time_window%month=g_s%ts_initialize
  g_s%time_window%day=g_s%ts_initialize
  g_s%time_window%n_day_x_month=g_s%ts_initialize  



  year(1)=ini_year
  month(1)=ini_month
  day(1)=ini_day
  time(1)=00
  tdiff=24
  index=1  
  g_s%time_window%date(index)=year(index)*10000+month(index)*100+day(index)
  g_s%time_window%year(index)=year(index)
  g_s%time_window%month(index)=month(index)
  g_s%time_window%day(index)=day(index)
  month_save=month(1)
  month_index=1
 g_s%time_window%n_day_x_month(month_index)=1
  do while(index <= g_s%dimensions%n_days)
     !write(*,*)"index",index
     !the dates subroutine is in the file dates.f
     CALL DATES(Year(index),Month(index),Day(index),Time(index),Tdiff, &
          Year1,Month1,Day1,Time1)
     index=index+1
     year(index)=year1
     month(index)=month1
     day(index)=day1
     time(index)=time1
     !write(*,*) "year1,month1,day1,time1,index",year1,month1,day1,time1,index
     !WRITE(*,*)"year1*10000+month1*100+day1",year1*10000+month1*100+day1
     g_s%time_window%date(index)=year(index)*10000+month(index)*100+day(index)
     g_s%time_window%year(index)=year(index)
     g_s%time_window%month(index)=month(index)
     g_s%time_window%day(index)=day(index)
     if(month(index) /= month_save .and. (month_index+1) <= g_s%time_window%last_month)then
        month_index = month_index+1
        !write(*,*)"g_s%time_window%date(index)", g_s%time_window%date(index)
        !write(*,*)"month_index",month_index
        g_s%time_window%n_day_x_month(month_index)=toindex_lo(100*month(index)+10000*year(index)+1,toindex_name,g_s)
        month_save=month(index)
     endif

     if(year1*10000+month1*100+day1 .eq. g_s%time_window%end_day) then
        exit
     end if
  enddo
  index=index-1

  !FILLING
  !g_s%time_window%date=year*10000+month*100+day
  !g_s%time_window%year=year
  !g_s%time_window%month=month
  !g_s%time_window%day=day
  

  !DEALLOCATE year, month and day
  DEALLOCATE(year, month, day, STAT= status)
  if (status == 0)then
     error_mex=""
  else
     error_mex = "DEALLOCATE(year, month, day, STAT= status)"
  endif
  call error(routine_name,status,2,error_mex) 
  
  
  Date_filling =status
END FUNCTION Date_filling


INTEGER FUNCTION create_date_index(start_date,n_days,date_index,g_s)
  IMPLICIT NONE
  INTEGER :: start_date,n_days
  INTEGER, ALLOCATABLE :: date_index(:)
  TYPE(global_settings_raobcore) :: g_s !! global settings variables
  INTEGER :: status          !!status of the function
  INTEGER, ALLOCATABLE :: year(:),month(:),day(:),time(:)
  INTEGER :: ini_day,ini_month,ini_year,index,tdiff
  INTEGER :: year1,month1,day1,time1,mod
  CHARACTER (len = 100):: routine_name
  CHARACTER (len = 200):: error_mex

  !the routine_name
  routine_name = "create_date_index"
  !!I set the status= 0
  status = 0
  
  ini_year=start_date/10000
  ini_month=(start_date-ini_year*10000)/100
  ini_day=mod(start_date,100)
  
  AllOCATE(year(n_days),month(n_days),day(n_days),time(n_days), STAT = status)
  error_mex = "  AllOCATE(year(n_days),month(n_days),day(n_days),time(n_days), STAT = status)"
 
  call error(routine_name,status,2,error_mex)
     
  year(1)=ini_year
  month(1)=ini_month
  day(1)=ini_day
  time(1)=00
  tdiff=24
  index=1
  do while(index <= n_days)
     !the dates subroutine is in the file dates.f
     CALL DATES(Year(index),Month(index),Day(index),Time(index),Tdiff, &
          Year1,Month1,Day1,Time1)
     index=index+1
     year(index)=year1
     month(index)=month1
     day(index)=day1
     time(index)=time1
     !write(*,*) "year1,month1,day1,time1,index",year1,month1,day1,time1,index
     !WRITE(*,*)"year1*10000+month1*100+day1",year1*10000+month1*100+day1
     if(year1*10000+month1*100+day1 .eq. g_s%time_window%end_day) then
        exit
     end if
  enddo
  index=index-1
  allocate(date_index(n_days),STAT = status)
  
  error_mex = "allocate(date_index(n_days),STAT = status)"
  
  call error(routine_name,status,2,error_mex)
  !INITIALIZING
  date_index=g_s%ts_initialize
  !FILLING
  date_index=year*10000+month*100+day
  !DEALLOCATE year, month and day
  DEALLOCATE(year, month, day, STAT= status)
  
  error_mex = "DEALLOCATE(year, month, day, STAT= status)"
  call error(routine_name,status,2,error_mex) 
  
  create_date_index=status
END FUNCTION create_date_index


INTEGER function todate_lo(index,routine_name,g_s)
IMPLICIT NONE
INTEGER:: index
CHARACTER (len = 100) :: routine_name
TYPE(global_settings_raobcore) :: g_s !type used to store all the global setting

!INTERNAL
INTEGER :: status          !!status of the function
CHARACTER(len=5)::index_str
INTEGER, ALLOCATABLE :: year(:),month(:),day(:),time(:)
INTEGER :: ini_day,ini_month,ini_year,tdiff,ind
INTEGER :: year1,month1,day1,time1,mod
CHARACTER(len=200)::error_mex

!the routine_name
routine_name = "todate_lo"
!!I set the status= 0
status = 0


index_check:if((index >  0) .and. (index <= g_s%dimensions%n_days)) then 
   todate_lo=g_s%time_window%year(index)*10000+g_s%time_window%month(index)*100+g_s%time_window%day(index)
else
   WRITE(index_str,'(I5)')index
   error_mex = "UNCORRECT input index->"//index_str
   status=-999
   call error(routine_name,status,2,error_mex)
endif index_check

return
end function todate_lo

integer function toindex_lo(date,routine_name,g_s)
IMPLICIT NONE
!EXTERNAL
INTEGER:: date
CHARACTER (len = 100) :: routine_name
TYPE(global_settings_raobcore) :: g_s !type used to store all the global setting

!INTERNAL
INTEGER :: status          !!status of the function
INTEGER, ALLOCATABLE :: year(:),month(:),day(:),time(:)
INTEGER :: ini_day,ini_month,ini_year,index,tdiff
INTEGER :: year1,month1,day1,time1,mod
CHARACTER(len=8)::date_str
CHARACTER(len=200)::error_mex
!the routine_name
routine_name = "toindex"
!!I set the status= 0
status = 0


!I check the date is inside my time window
date_check:if((date >= g_s%time_window%start_day) .and. (date <= g_s%time_window%end_day)) then
   index=1
   do while (date > g_s%time_window%date(index))
      index=index+1
   enddo
   toindex_lo=index
   
else 
   WRITE(date_str,'(I8.8)')date
   error_mex = "UNCORRECT input date->"//date_str
   status=-999
   call error(routine_name,status,2,error_mex) 
endif date_check
return
end function toindex_lo


END MODULE JSON_file
