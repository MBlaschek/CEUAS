SUBROUTINE BIASCOR_ERA40 (LD_LBC,K_NST,P_SPP,P_ST,P_STNBC, &
 & CD_CIDENT,K_IY,K_IM,K_ID,K_IH,P_RLAT,P_RLON,K_IRSTYP,  &
 & K_IMISS,P_RMISS,CL_SOLAR,CL_HOMOGEN)  

!!leo USE PARKIND1  ,ONLY : JPIM     ,JPRB     ,JPRM
!!leo USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
                     
!!**** *BIASCOR_ERA40* 

!!       PURPOSE.
!!      ------------

!!       BIAS CORRECTION OF RADIOSONDE TEMPERATURES FOR 

!!          EEEEE  RRRR      A         4    000
!!          E      R   R    A A       44   0   0
!!          EEEE   RRRR    A   A     4 4   0   0
!!          E      R R     AAAAA    44444  0   0
!!          EEEEE  R   R   A   A       4    000

!!       INTERFACE.
!!      ------------

!!          CALL BIASCOR_ERA40 (LBC,NST,SPP,ST,STNBC,
!!     X                  CIDENT,IM,ID,IH,RLAT,RLON,IRSTYP,
!!     X                  IMISS,RMISS)

!!        INPUT
!!         NST       NUMBER OF LEVELS (4 BYTE INTEGER)
!!         SPP       ARRAY WITH PRESSURE VALUES (Pa) (8 BYTE REAL)
!!         ST        ARRAY WITH T VALUES (8 BYTE REAL)
!!         CIDENT    STATION IDENTIFIER  (CHARACTER)
!!         IM        MONTH (4 BYTE INTEGER)
!!         ID        DAY   (4 BYTE INTEGER)
!!         IH        HOUR  (4 BYTE INTEGER)
!!         RLAT      LATITUDE (8 BYTE REAL)
!!         RLON      LONGITUDE (8 BYTE REAL)
!!         IRSTYP    RADIOSONDE TYPE (BUFR CODE TABLE 002011) 
!!                   (Not need for ERA)
!!         IMISS     MISSING VALUE FOR INTEGERS (4BYTE INTEGER)
!!         RMISS     MISSING VALUE FOR REALS (8 BYTE REAL) 

!!        OUTPUT
!!         LBC       LOGICAL TO STATE IF BIAS COR. WAS SUCCESSFUL
!!         STNBC     ARRAY WITH BIAS CORRECTED T VALUES

!!       METHOD.
!!      ---------

!!        READ BIAS CORRECTION TABLES:
!!            STGROUP.T : STATION GROUP TABLE
!!            CORRECT.T : ARRAY OF CORRECTIONS DEPENDING ON
!!                        SOLAR ELEVATION AND PRESSURE LEVEL
!!            COUNTRY.T : DEFINITION OF CATEGORIES

!!        1ST A STATION GROUP 0F THE DATA IS DETECTED.
!!        2ND A CATEGORY INCLUDING THE STATION GROUP IS DETECTED.
!!        3RD IF THE CATEGORY IS ONE FOR CORRECTION, APPLY CORRECTION.

!!        FIRST 32 CHARACTERS ARE COMPARED TO DETECT CATEGORY.
!!        'CATG' FROM COUNTRY.T AND 'YMNAMBC' FROM CORRECT.T

!!       EXTERNALS.
!!      ------------

!!        DIURNAL          CALCULATE SOLAR ELEVATION
!!        PNTERP           INTERPOLATE TO PRESSURE LEVEL

!!       REFERENCE.
!!      ------------

!!        RADIOSONDE BIAS CORRECTION
!!        OD MEMORANDUM BY B. STRAUSS  22.12.92

!!       AUTHOR.
!!      ---------

!!        B. NORRIS       JULY  1991   APPLY CORRECTIONS TO AOF

!!       MODIFICATIONS.
!!      ----------------

!!        M. DRAGOSAVAC   APRIL 1993   IMPLEMENT IN PRE-PROCESSING
!!        B. NORRIS       JULY  1998   INTERFACE TO DATA ASSIMILATION
!!        K. ONOGI        MARCH 2000   MODIFIED FOR ERA-40
!!        K. ONOGI      OCTOBER 2000   MAKE NUMBER OF B.COR. CATEGORIES 8 TO 4
!!                                     AS SAME AS THE CORRECTION TABLE
!!        S. SAARINEN  NOVEMBER 2000   CLEANING UP PLUS IMPLICIT NONE
!!        M. Hamrud      01-Oct-2003   CY28 Cleaning
!!      L. Haimberger  NOVEMBER 2004   Allow use of tables from  radiosonde homogenization 
!!                                    (see ERA-40 Project report series 30)
!!-------------------------------------------------------------------

IMPLICIT NONE

!!     Table values are expected to be saved
!!     after they are read

 INTEGER,PARAMETER :: JPIM=4
 INTEGER,PARAMETER :: JPRM=4
 INTEGER,PARAMETER :: JPRB=4

SAVE

LOGICAL           ,INTENT(OUT)   :: LD_LBC 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_NST 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_SPP(K_NST) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_ST(K_NST) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_STNBC(K_NST) 
CHARACTER(LEN=*)  ,INTENT(IN)    :: CD_CIDENT 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_IM 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_ID 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_IH 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_IY 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_RLAT 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_RLON 
INTEGER(KIND=JPIM)               :: K_IRSTYP !! Argument NOT used
INTEGER(KIND=JPIM)               :: K_IMISS !! Argument NOT used
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_RMISS 
CHARACTER(LEN=5)   ,INTENT(IN)    :: CL_SOLAR,CL_HOMOGEN
INTEGER(KIND=JPIM)::I_NLEVBC   ,I_NCORBC  ,I_NCORTB   
PARAMETER (I_NLEVBC=16,I_NCORBC=4,I_NCORTB=4)
INTEGER(KIND=JPIM)::I_NSTNBC     ,I_NSONBC    ,I_NSSNBC     
PARAMETER (I_NSTNBC=3100,I_NSONBC=200,I_NSSNBC=300)

INTEGER(KIND=JPIM)::I_NXGC      , I_NXGP      
PARAMETER (I_NXGC=20000, I_NXGP=4000)
INTEGER(KIND=JPIM)::I_NXDG    , I_NXCT      
PARAMETER (I_NXDG=100, I_NXCT=1000)
INTEGER(KIND=JPIM)::I_MTBL   , I_MCOR   , I_MSGT    , I_MLEO   
PARAMETER (I_MTBL=65, I_MCOR=66, I_MSGT=67, I_MLEO=68)

REAL(KIND=JPRM)     , DIMENSION (I_NCORBC,I_NLEVBC,I_NSONBC) :: Z_RCORBC
INTEGER(KIND=JPIM)  , DIMENSION (I_NLEVBC,I_NSONBC) :: ILEVBC
REAL(KIND=JPRM)     , DIMENSION (I_NCORTB) :: Z_WKCOR
INTEGER(KIND=JPIM) :: IPLAT,IPLON,IRLAT,IRLON
CHARACTER(LEN= 5) :: CL_YDUMMYBC
CHARACTER(LEN=32), DIMENSION (I_NSONBC) :: CL_YSNAMBC
CHARACTER(LEN=58) :: CL_CDATE

LOGICAL :: LL_FILES_NOT_READ = .TRUE.
INTEGER(KIND=JPIM) :: I_NSOBC
      
!! --- ARRAYS FOR STATION GROUP TABLE
REAL(KIND=JPRB)     , DIMENSION (I_NXGC) :: Z_PLAT,Z_PLON
CHARACTER(LEN= 6), DIMENSION (I_NXGC) :: CL_CSNG,CL_CSNO
CHARACTER(LEN= 6) :: CL_C1,CL_C4
CHARACTER(LEN= 1) :: CL_C0
CHARACTER(LEN= 5) :: CL_C2
CHARACTER(LEN=19) :: CL_C3
INTEGER(KIND=JPIM)  , DIMENSION (I_NXGP) :: JSGC,JEGC,I_MREP
REAL(KIND=JPRM)     , DIMENSION (I_NXGP) :: Z_QLAT,Z_QLON
INTEGER(KIND=JPIM)  , DIMENSION (I_NXCT) :: I_NTGT
INTEGER(KIND=JPIM)  , DIMENSION (I_NXGP,I_NXCT) :: I_MTGT

!! --- ARRAYS FOR INDEX TABLE
INTEGER(KIND=JPIM)  , DIMENSION (I_NXDG,I_NXCT) :: IDSTA,IDEND, &
 & I_LATS,I_LATE,I_LONS,I_LONE  
INTEGER(KIND=JPIM)  , DIMENSION (I_NXDG) :: I_NCMB,I_NCDU
INTEGER(KIND=JPIM)  , DIMENSION (I_NXCT) :: I_NCCC,I_MCAT
CHARACTER(LEN=64), DIMENSION (I_NXDG,I_NXCT) :: CL_CNTRY
CHARACTER(LEN=64), DIMENSION (I_NXCT) :: CL_CDG
CHARACTER(LEN=64) ::  CL_C64
CHARACTER(LEN=18) ::  CLTLN
CHARACTER(LEN=32), DIMENSION (I_NXCT) :: CL_CATG

LOGICAL, PARAMETER :: LGRP = .TRUE.

!! --- LEO VARIABLES
INTEGER(KIND=JPIM) :: CLEOUNIT   = 11,ierr,icsn
LOGICAL :: CLEO_FOUND = .FALSE.
CHARACTER(LEN=50) :: CLEONAME = 'table4'
CHARACTER(LEN=120) :: zeile
CHARACTER(LEN=1)  :: CCORR
CHARACTER(LEN=2)  :: CPMAX
CHARACTER(LEN=4)  :: CDSTRING
INTEGER(KIND=JPIM), PARAMETER :: pmax    = 16
INTEGER(KIND=JPIM), PARAMETER :: parmax   = 2
INTEGER(KIND=JPIM), PARAMETER :: estatmax = 3070
REAL(KIND=JPRB) :: ifs_rasocorrs(pmax,parmax,estatmax)
INTEGER(KIND=JPIM) :: ewmonrs(estatmax)
INTEGER(KIND=JPIM) :: ISTAT,IB
REAL(KIND=JPRB)    :: LEO_PLEV(pmax),iplevs(pmax),cleo_bias(pmax)
real(kind=JPRB)    :: ewmolats(estatmax),ewmolons(estatmax)
real(kind=JPRB)    :: hilf(pmax,parmax)
integer  iyear,imonth,itime,iday,idatum,datum,idum,ip,ipmax,istatmax

      
!! --- MISCELLANEOUS DECLARATIONS
INTEGER(KIND=JPIM) :: IOS, I, IR,JDP,IDS,IDE,I_LAS,I_LAE,I_LOS,I_LOE
INTEGER(KIND=JPIM) :: I_M, JDG, I_N, I_NGC, I_NGP, I1, I2,I6
REAL(KIND=JPRM) :: Z_R1,Z_R2
INTEGER(KIND=JPIM) :: ISY,ISM,ISD,ISH,IEY,IEM,IED,IEH,I5
INTEGER(KIND=JPIM) :: I_NCNT, J, I_KG, I_KM, IXT, ISON, I_K, ISNG, ICOL
REAL(KIND=JPRB) :: Z_ANGLEBC, PP, Z_POS, Z_RADD, Z_LEO
INTEGER(KIND=JPIM) :: ILOW, IHGH,KIHSAVE
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!!leo #include "abor1.intfb.h"
!!leo #include "diurnal.intfb.h"
!!leo #include "pnterp.intfb.h"

!!  -------------------------------------------------------------------

!!     1.  OPEN AND READ TABLE FILES
!!     -------------------------------

!!leo IF (LHOOK) CALL DR_HOOK('BIASCOR_ERA40',0,ZHOOK_HANDLE)

100 CONTINUE

LD_LBC=.FALSE.
IF(KIHSAVE .ne. K_IY) THEN
   LL_FILES_NOT_READ=.true.
   KIHSAVE=K_IY
ENDIF
IF (LL_FILES_NOT_READ) THEN

WRITE(CDSTRING,'(I4)') K_IY

!!        WRITE(*,*)
!!        WRITE(*,*) '------------------------------'
!!        WRITE(*,*) '  RADIOSONDE BIAS CORRECTION  '
!!        WRITE(*,*) '------------------------------'
!!        WRITE(*,*)

  IF(CL_SOLAR(1:2) .eq. 'ON' ) THEN

!!    OPEN(UNIT=I_MTBL, FILE='table1', &
         OPEN(UNIT=I_MTBL, FILE='/home/srvx2/leo/tables/rsbias/country.t', &
   & IOSTAT=IOS,ERR=965,STATUS='OLD')  

!!    OPEN(UNIT=I_MCOR, FILE='table2', &
!!         OPEN(UNIT=I_MCOR, FILE='/home/srvx2/leo/tables/rsbias/T_correct'//CDSTRING//'010100', &
         OPEN(UNIT=I_MCOR, FILE='/home/srvx2/leo/tables/rsbias/wrk_biastables_0001_new/corcand_'//CDSTRING//'.t', &
   & IOSTAT=IOS,ERR=966,STATUS='OLD')  

!!    OPEN(UNIT=I_MSGT, FILE='table3', &
         OPEN(UNIT=I_MSGT, FILE='/home/srvx2/leo/tables/rsbias/stgroup.t', &
   & IOSTAT=IOS,ERR=967,STATUS='OLD')  



!!       1.1 READ CATEGORY DEFINITION TABLE

  DO I=1,2000
    READ(I_MTBL,'(I4)') IR
    IF (IR == 10) GOTO 10
  ENDDO
  10   CONTINUE
  I_NCNT = 0
  DO I=1,I_NXDG
    I_NCMB(I) = 0
  ENDDO
  DO I=1,I_NXCT
    READ(I_MTBL,'(I4,I3,2I6,1X,2I3,2I4,1X,A64)') &
     & IR,JDP,IDS,IDE,I_LAS,I_LAE,I_LOS,I_LOE,CL_C64  
    IF (IR == 999) GOTO 20
    IF (IR == 1) THEN
      IF (I_LAS == 0.AND.I_LAE == 0.AND.I_LOS == 0.AND.I_LOE == 0) THEN
        I_LAS =  -90
        I_LAE =   90
        I_LOS = -180
        I_LOE =  180
      ENDIF
      IF (.NOT.LGRP) JDP = 0
      IF (JDP == 0) THEN
        I_NCNT = I_NCNT+1
        I_NCCC(I_NCNT) = 1
        I_MCAT(I_NCNT) = 0
        IDSTA(1,I_NCNT) = IDS
        IDEND(1,I_NCNT) = IDE
        I_LATS(1,I_NCNT)  = I_LAS
        I_LATE(1,I_NCNT)  = I_LAE
        I_LONS(1,I_NCNT)  = I_LOS
        I_LONE(1,I_NCNT)  = I_LOE
        CL_CNTRY(1,I_NCNT) = CL_C64
      ELSE
        IF (I_NCMB(JDP) == 0) THEN
          I_NCNT = I_NCNT+1
          I_NCDU(JDP) = I_NCNT
          I_NCMB(JDP) = 1
          I_NCCC(I_NCDU(JDP)) = 1
          I_MCAT(I_NCNT) = JDP
          IDSTA(I_NCMB(JDP),I_NCDU(JDP)) = IDS
          IDEND(I_NCMB(JDP),I_NCDU(JDP)) = IDE
          I_LATS(I_NCMB(JDP),I_NCDU(JDP))  = I_LAS
          I_LATE(I_NCMB(JDP),I_NCDU(JDP))  = I_LAE
          I_LONS(I_NCMB(JDP),I_NCDU(JDP))  = I_LOS
          I_LONE(I_NCMB(JDP),I_NCDU(JDP))  = I_LOE
          CL_CNTRY(I_NCMB(JDP),I_NCDU(JDP)) = CL_C64
        ELSE
          I_NCMB(JDP) = I_NCMB(JDP)+1
          I_NCCC(I_NCDU(JDP)) = I_NCCC(I_NCDU(JDP))+1
          IDSTA(I_NCMB(JDP),I_NCDU(JDP)) = IDS
          IDEND(I_NCMB(JDP),I_NCDU(JDP)) = IDE
          I_LATS(I_NCMB(JDP),I_NCDU(JDP))  = I_LAS
          I_LATE(I_NCMB(JDP),I_NCDU(JDP))  = I_LAE
          I_LONS(I_NCMB(JDP),I_NCDU(JDP))  = I_LOS
          I_LONE(I_NCMB(JDP),I_NCDU(JDP))  = I_LOE
          CL_CNTRY(I_NCMB(JDP),I_NCDU(JDP)) = CL_C64
        ENDIF
      ENDIF
    ENDIF
  ENDDO
  20   CONTINUE

!!       1.1.1 PUT LAT LON LIMIT ON CNTRY

  DO I_M=1,I_NCNT
    DO I=1,I_NCCC(I_M)
      IF (I_LATS(I,I_M) /= -90.OR.I_LATE(I,I_M) /= 90.OR. &
         & I_LONS(I,I_M) /= -180.OR.I_LONE(I,I_M) /= 180) THEN  
        CLTLN = 'xxS-xxN,xxxW-xxxE '
        WRITE(CLTLN( 1: 2),'(I2)') ABS(I_LATS(I,I_M))
        IF (I_LATS(I,I_M) < 0) CLTLN( 3: 3) = 'S'
        IF (I_LATS(I,I_M) == 0) CLTLN( 3: 3) = ' '
        IF (I_LATS(I,I_M) > 0) CLTLN( 3: 3) = 'N'
        WRITE(CLTLN( 5: 6),'(I2)') ABS(I_LATE(I,I_M))
        IF (I_LATE(I,I_M) < 0) CLTLN( 7: 7) = 'S'
        IF (I_LATE(I,I_M) == 0) CLTLN( 7: 7) = ' '
        IF (I_LATE(I,I_M) > 0) CLTLN( 7: 7) = 'N'
        WRITE(CLTLN( 9:11),'(I3)') ABS(I_LONS(I,I_M))
        IF (I_LONS(I,I_M) < 0) CLTLN(12:12) = 'W'
        IF (I_LONS(I,I_M) == 0) CLTLN(12:12) = ' '
        IF (I_LONS(I,I_M) > 0) CLTLN(12:12) = 'E'
        WRITE(CLTLN(14:16),'(I3)') ABS(I_LONE(I,I_M))
        IF (I_LONE(I,I_M) < 0) CLTLN(17:17) = 'W'
        IF (I_LONE(I,I_M) == 0) CLTLN(17:17) = ' '
        IF (I_LONE(I,I_M) > 0) CLTLN(17:17) = 'E'
!!               CNTRY(I,M) = CNTRY(I,M)(1:20)//CLTLN
        CL_CNTRY(I,I_M) = CLTLN//CL_CNTRY(I,I_M)(1:20)

      ENDIF
    ENDDO
  ENDDO

!!       1.1.2 GATHER SOME GROUPS INTO ONE GROUP

  DO I=1,I_NXCT
    CL_CDG(I) = ' '
  ENDDO
  IF (LGRP) THEN
    DO I=1,2000
      READ(I_MTBL,'(I4)') IR
      IF (IR == 20) GOTO 30
    ENDDO
    30     CONTINUE
    DO I=1,I_NXDG
      READ(I_MTBL,'(I4,I3,1X,A64)') IR,JDG,CL_C64
      IF (IR == 999) GOTO 40
      IF (IR == 1) THEN
        CL_CDG(I_NCDU(JDG)) = CL_C64
      ENDIF
    ENDDO
    40     CONTINUE
  ELSE
    DO I=1,I_NXDG
      I_NCDU(I)=0
    ENDDO
  ENDIF
  DO I_N=1,I_NCNT
    IF (I_MCAT(I_N) == 0) THEN
      CL_CATG(I_N) = CL_CNTRY(1,I_N)(1:32)
    ELSE
      CL_CATG(I_N) = CL_CDG(I_NCDU(I_MCAT(I_N)))(1:32)
    ENDIF
  ENDDO

!!       1.1.3 MONITOR THE CATEGORIES
!!       DO N=1,NCNT
!!         DO I=1,NXDG
!!         IF (IDSTA(I,N).NE.0) &
!!       & WRITE(*,'(I5,2I7,4I5,1X,A32,2X,A32)') &
!!             & I,IDSTA(I,N),IDEND(I,N),LATS(I,N),LATE(I,N), &
!!             & LONS(I,N),LONE(I,N),CNTRY(I,N)(1:32),CATG(N)
!!         END DO
!!       END DO

!!       1.2 READ STATION GROUP TABLE
  I_NGC = 0  !! NUMBER OF RECORD
  I_NGP = 0  !! NUMBER OF GROUP
  350   CONTINUE
  READ(I_MSGT, &
   & '(I5,I4,A1,A6,F7.2,F8.2,I8,2(1X,I4,3I2),I7,2X,A6)', &
   & END=190) &
   & I1,I2,CL_C0,CL_C1,Z_R1,Z_R2,I6, &
   & ISY,ISM,ISD,ISH,IEY,IEM,IED,IEH,I5,CL_C4  
!!           WRITE(*, &
!!           & '(I5,I4,A1,A6,F7.2,F8.2,A5,2(1X,4I2),I7,2X,A19,2X,A6)') &
!!           &  I1,I2,C0,C1,R1,R2,C2, &
!!           &  ISY,ISM,ISD,ISH,IEY,IEM,IED,IEH,I5,C3,C4
  I_NGC = I_NGC+1
  Z_PLAT(I_NGC) = Z_R1
  Z_PLON(I_NGC) = Z_R2
  CL_CSNO(I_NGC) = CL_C1
  IF (I2 == 1) THEN
    I_NGP = I_NGP+1
    JSGC(I_NGP) = I_NGC
  ENDIF
  JEGC(I_NGP) = I_NGC
  IF (CL_C4 == '      ') THEN
    CL_CSNG(I_NGC) = CL_C1
  ELSE
    CL_CSNG(I_NGC) = CL_C4
  ENDIF
  IF (CL_C0 == '#') THEN
    I_MREP(I_NGP) = I_NGC
    Z_QLAT(I_NGP) = Z_R1
    Z_QLON(I_NGP) = Z_R2
  ENDIF
  GOTO 350
  190   CONTINUE

!!       1.3  READ CORRECTION TABLE

  I_NSOBC=0
  READ (I_MCOR,'(A)') CL_CDATE

  IF(CL_HOMOGEN(1:2) .EQ. 'ON') THEN
!! If L. Haimbergers correction table is used, it has to be used with a special version of
!! solar angle bias correction tables. See ERA-40 project report series Nr. 30 by L. Haimberger
    IF(CL_CDATE(40:58) .NE. 'with homogenization') THEN
      CALL ABOR1('BIAS CORRECTION ERROR: wrong table2 for solar angle correction + homogenization')
    ENDIF
  ELSE
    IF(CL_CDATE(40:58) .EQ. 'with homogenization') THEN
      CALL ABOR1('BIAS CORRECTION ERROR: wrong table2 for solar correction only')
    ENDIF
  ENDIF

  READ (I_MCOR,*)

  132   CONTINUE
  IF (I_NSOBC < I_NSONBC) THEN
    I_NSOBC=I_NSOBC+1
    READ (I_MCOR,'(A)',END=138) CL_YSNAMBC(I_NSOBC)
!!           WRITE(*,'(2A)') YSNAMBC(NSOBC)

    DO J=1,I_NLEVBC
      READ (I_MCOR,'(I5,8F8.2)',END=138) &
       & ILEVBC(J,I_NSOBC),(Z_WKCOR(I_K),I_K=1,I_NCORTB)  

!!         WRITE(*,'(I5,8F8.2)') &
!!               & ILEVBC(J,NSOBC),(WKCOR(K),K=1,NCORTB)

!!     CORRECTION TABLE
!!     ~-7.5 : -7.5~7.5 : 7.5~22.5 : 22.5~

      DO I=1,I_NCORBC
        Z_RCORBC(I,J,I_NSOBC) = Z_WKCOR(I)
      ENDDO

    ENDDO
    READ (I_MCOR,'(A)',END=138) CL_YDUMMYBC
    GOTO 132
  ELSE
    WRITE (0,'(1X,A,I5)') &
     & ' BIAS CORRECTION ERROR: DIMENSION NSON TOO SMALL :', &
     & I_NSONBC  
    CALL ABOR1('BIAS CORRECTION ERROR: DIMENSION NSON TOO SMALL')
    CALL ABORT
  ENDIF
  138   CONTINUE
  I_NSOBC = I_NSOBC-1

!!         WRITE (0,*) &
!!         & ' BIAS CORRECTION: NUMBER OF CATEGORIES IN CORRECTION TABLE ', &
!!         &  NSOBC


        !!!WRITE(6,*)'SUCCESSFULLY READ TABLE 1-3'

  CLOSE (I_MTBL)
  CLOSE (I_MCOR)
  CLOSE (I_MSGT)

  ENDIF !! CL_SOLAR

  IF(CL_HOMOGEN(1:2) .EQ. 'ON') THEN

!!         OPEN(UNIT=MSGT, FILE='leobiascor.t', &
  OPEN(UNIT=I_MLEO, FILE='table4', &
   & IOSTAT=IOS,ERR=968,STATUS='OLD')  

!! -------------------------------------
!! 
!!     Read L. Haimbergers's homogenized temp corrections     
!! -------------------------------------
   
 LEO_PLEV=(/10.,20.,30.,50.,70.,100.,150.,200.,250.,300.,400.,500.,700.,850.,925.,1000./)
 do ip=1,pmax
    ilevbc(ip,1)=leo_plev(pmax-ip+1)
 enddo

      cleo_bias=0.

!! read and print comment line(s) in header, if available
  zeile='#'
  do while(zeile(1:1) .eq. '#')
    read(I_MLEO,'(A120)') zeile
    if(zeile(1:1) .eq. '#') write(*,*) zeile
  enddo

  read(zeile,'(I5,I3,16F6.0)') istatmax,ipmax,iplevs

  if(pmax .ne. ipmax .or. any(LEO_PLEV .ne. iplevs) .or. istatmax .ne. estatmax) then
     write(6,*)'istatmax,estatmax',istatmax,estatmax
     write(6,*)'pmax,ipmax',pmax,ipmax
     write(6,*)'plevs',LEO_PLEV
     write(6,*)'iplevs',iplevs
 CALL ABOR1('WRONG DIMENSIONS OF LEOBIASCOR.T')
  endif

  DATUM=K_IY*1000000+K_IM*10000+K_ID*100

  write(*,*) 'HOMOGENIZING BIAS CORRECTION:',datum
  write(cpmax,'(I2.2)') ipmax

  do istat=1,estatmax
    read(I_MLEO,'(I5,I6,2F9.2,1X,A1,I3)',ERR=968) idum,ewmonrs(istat),ewmolats(istat),ewmolons(istat),ccorr,ib

    if(ccorr .eq. "Y") then
      do i=1,ib
        read(I_MLEO,'(I4,3I2.2,'//cpmax//'F6.2)',ERR=968) iyear,imonth,iday,itime,hilf(:,1)
        read(I_MLEO,'(I4,3I2.2,'//cpmax//'F6.2)',ERR=968) iyear,imonth,iday,itime,hilf(:,2)
        idatum=iyear*1000000+imonth*10000+iday*100
        if(idatum .lt. datum) then
          do ip=1,pmax
            ifs_rasocorrs(ip,:,istat)=hilf(pmax-ip+1,:)
          enddo
        endif
      enddo
    write(*,*) 'read ',ewmonrs(istat),ewmolats(istat),ewmolons(istat),ccorr,ib,ifs_rasocorrs(5,:,istat)
    endif

  enddo

  CLOSE (I_MLEO)

  ENDIF !! CL_HOMOGEN

  LL_FILES_NOT_READ = .FALSE.
ENDIF

!!-----------------------------------------------------------------------------


!!     1.4  INITIALLY SET CORRECTED = ORIGINAL

DO I=1,K_NST
  P_STNBC(I)=P_ST(I)
ENDDO

IF(CL_SOLAR(1:2) .eq. 'ON') THEN

!!  ------------------------------------------------------------------ 

!!  Solar angle bias correction
!!
!!     2.   DETERMINE TABLES AND TARGET

!!     2.1  DETECT THE STATION GROUP
I_KG = 0 ; I_KM = 0
DO I=1,I_NGP
  DO J=JSGC(I),JEGC(I)
    IPLAT = NINT(Z_PLAT(J)*100)
    IPLON = NINT(Z_PLON(J)*100)
    IRLAT = NINT(P_RLAT*100)
    IRLON = NINT(P_RLON*100)
    IF (IRLAT == IPLAT.AND.IRLON == IPLON) THEN
      IF (CD_CIDENT(1:5) == CL_CSNO(J)(1:5)) THEN
        I_KG = I               !! STATION GROUP
        I_KM = J-JSGC(I)+1     !! MEMBER 
        GOTO 440
      ENDIF
    ENDIF
  ENDDO
ENDDO
440 CONTINUE

IF (I_KG*I_KM == 0) THEN
  WRITE(*,'(A,A6,F7.2,F8.2,I6,I3,2A)') &
   & 'BIASCOR: K_ID,LAT,LON,GRP,MEM= ',CD_CIDENT,P_RLAT,P_RLON,I_KG,I_KM, &
   & ' | THE SPECIFIED DATA WAS NOT FOUND ', &
   & 'IN THE STATION GROUP DEFINITION TABLE.'  
  GOTO 500
      !! THE SPECIFIED DATA WAS NOT FOUND IN THE STATION GROUP DEFINITION TABLE.
ENDIF

!!     2.2  PICK UP THE TARGET CATEGORY

DO I_M=1,I_NCNT
  I_NTGT(I_M) = 0
  READ(CL_CSNG(I_MREP(I_KG))(1:5),'(I5)') ISNG
  IXT = 0
  DO I=1,I_NCCC(I_M)
    IF (ISNG >= IDSTA(I,I_M).AND.ISNG <= IDEND(I,I_M)) THEN
      IF ( (I_LONS(I,I_M) <= I_LONE(I,I_M).AND. &
         & (Z_QLAT(I_KG) >= REAL(I_LATS(I,I_M)).AND. &
         & Z_QLAT(I_KG) <= REAL(I_LATE(I,I_M)).AND. &
         & Z_QLON(I_KG) >= REAL(I_LONS(I,I_M)).AND. &
         & Z_QLON(I_KG) <= REAL(I_LONE(I,I_M))))    &
         & .OR.                               &
         & (I_LONS(I,I_M) > I_LONE(I,I_M).AND.       &
         & (Z_QLAT(I_KG) >= REAL(I_LATS(I,I_M)).AND. &
         & Z_QLAT(I_KG) <= REAL(I_LATE(I,I_M)).AND. &
         & Z_QLON(I_KG) >= REAL(I_LONS(I,I_M)).OR.  &
         & Z_QLON(I_KG) <= REAL(I_LONE(I,I_M)))) ) THEN  
        IXT = IXT+1
      ENDIF
    ENDIF
  ENDDO
  IF (IXT > 0) THEN
    I_NTGT(I_M) = I_NTGT(I_M)+1
    I_MTGT(I_NTGT(I_M),I_M) = I_M
  ENDIF
ENDDO

!!     DO M=1,NCNT
!!     DO K=1,NTGT(M)
!!       WRITE(*,*) M,K,' ',CATG(MTGT(NTGT(M),M))
!!     END DO
!!     END DO

200 CONTINUE

!!     2.3  DETERMINE SECTION OF CORRECT.T

DO I_M=1,I_NCNT
  DO I_K=1,I_NTGT(I_M)
    DO I=1,I_NSOBC
      IF (CL_CATG(I_M) == CL_YSNAMBC(I)) THEN
        ISON=I
        GO TO 234
      ENDIF
    ENDDO

    WRITE(*,'(A,A6,F7.2,F8.2,I6,I3,A,I3,2X,A)') &
     & 'BIASCOR: K_ID,LAT,LON,P_ST-GP&MB=', &
     & CD_CIDENT,P_RLAT,P_RLON,I_KG,I_KM, &
     & ' | NO COR. TABLE FOR CATEGORY NO.',I_M,CL_CATG(I_M)    

    GOTO 400
    234   CONTINUE
                                       
!!     2.4  CALCULATE SOLAR ELEVATION

    CALL DIURNAL ( K_IM, K_ID, K_IH, P_RLAT,P_RLON, Z_ANGLEBC )
                                   
!!     2.5  CALCULATE CORRECT COLUMN OF CORRECTION TABLE

    IF (Z_ANGLEBC < -7.5_JPRB)                          ICOL=1
    IF (Z_ANGLEBC >= -7.5_JPRB.AND.Z_ANGLEBC < 7.5_JPRB) ICOL=2
    IF (Z_ANGLEBC >= 7.5_JPRB.AND.Z_ANGLEBC < 22.5_JPRB) ICOL=3
    IF (Z_ANGLEBC >= 22.5_JPRB)                          ICOL=4

    WRITE(*,'(A,A6,F7.2,F8.2,I5,I3,I6,I3,A,I3,2X,A,A,F6.2,A,I2)') &
     & 'BIASCOR: K_ID,LAT,LON,P_ST-GP&MB=', &
     & CD_CIDENT,P_RLAT,P_RLON,K_IY,K_IM,I_KG,I_KM, &
     & ' | COR. TABLE NO.',ISON,CL_YSNAMBC(ISON), &
     & ' | SOL.EV.=',Z_ANGLEBC, &
     & ' | CLASS',ICOL   

!!     2.6  LOOP OVER LEVELS

    DO I=1,K_NST
      IF (P_SPP(I) /= P_RMISS.AND.P_ST(I) /= P_RMISS) THEN

!!     2.7  CALCULATE CORRECT ROW OF CORRECTION TABLE

        PP=P_SPP(I)/100.0
        CALL PNTERP (ILEVBC(1,ISON),I_NLEVBC,PP,Z_POS)
        ILOW=Z_POS
        IHGH=ILOW+1


!!  -------------------------------------------------------------
                                         
!!                3.  Solar angle BIAS CORRECTION
!!               ---------------------

        Z_RADD=Z_RCORBC(ICOL,ILOW,ISON)*(REAL(IHGH)-Z_POS) &
         & +Z_RCORBC(ICOL,IHGH,ISON)*(Z_POS-REAL(ILOW))  


        P_STNBC(I)=P_ST(I) - Z_RADD
!!                         !!!     

!!   !!! NOTICE !!!

!!   Please note that the values in the bias correction table have
!!   opposite sign of correction (the same sign of biases)
!!   Values taken from the tables must be subtracted.
!!   This is opposite sign from operational bias correction.
!!   In operation, values are added to observations.

!!   We decided to use oppsite from operational one because we need 
!!   manual adjustment of correction tables.
!!   It was decided to avoid manual errors (confusion of signs).

!!   Please allow us to use opposite sign from operational 'biascor'.

!!                                                 2000.3.22 K.ONOGI

!!              WRITE(*,'(2A,2X,3F10.2)') &
!!              & 'BIASCOR:',CD_CIDENT(1:5),P_SPP(I),P_ST(I),P_STNBC(I)

      ENDIF
    ENDDO

!!       WRITE(*,'(A,10F7.2,;,/,(10(5X,10F7.2/)))') 'OBS.:',(ST(I),I=1,NST)
!!       WRITE(*,'(A,10F7.2,;,/,(10(5X,10F7.2/)))') 'COR.:',(STNBC(I),I=1,NST)

    LD_LBC=.TRUE.

!!--------------------------------------------------------------------
                                              
!!                   4.  EXIT
!!                  -----------

    400   CONTINUE

  ENDDO
ENDDO

ENDIF !! CL_SOLAR

IF(CL_HOMOGEN(1:2) .EQ. 'ON') THEN


!!     5. Use corrections determined by homogenization of radiosondes
!!

          READ(CD_CIDENT(1:5),'(I5)',IOSTAT=IERR) ICSN
          IF (IERR /= 0 ) THEN
            WRITE(6,*)'Could not make integer of ',CD_CIDENT(1:5)
          ENDIF
     
         cleo_found=.false.
         cleo_bias=0.
         ISTAT_LOOP : do istat=1,estatmax
           if(ICSN .eq. ewmonrs(istat)) then
     WRITE(*,*)'Station found:',ICSN
           if(K_IH .eq. 0) then
             cleo_bias=ifs_rasocorrs(:,1,istat)
           endif
           if(K_IH .eq. 12) then
             cleo_bias=ifs_rasocorrs(:,2,istat)
           endif
           cleo_found=.true.
           exit istat_loop
        endif
     enddo ISTAT_LOOP

         WRITE(*,'(A27,I5,A2,16F6.2)') 'HOMOGENIZING CORRECTION FOR',ICSN,': ',cleo_bias

!!     2.6  LOOP OVER LEVELS

    DO I=1,K_NST
      IF (P_SPP(I) /= P_RMISS.AND.P_ST(I) /= P_RMISS) THEN

!!     2.7  CALCULATE CORRECT ROW OF CORRECTION TABLE

        PP=P_SPP(I)/100.0
        CALL PNTERP (ILEVBC(1,1),I_NLEVBC,PP,Z_POS)
        ILOW=Z_POS
        IHGH=ILOW+1


!!        write(*,*) 'leo only',PP,ILOW,LEO_PLEV(pmax-ILOW+1),cleo_bias(ILOW),IHGH,LEO_PLEV(pmax-IHGH+1),cleo_bias(IHGH)

!!  -------------------------------------------------------------
                                         
!!                6.  BIAS CORRECTION using tables from homogenization
!!               ---------------------

        Z_LEO=cleo_bias(ILOW)*(REAL(IHGH)-Z_POS) &
         & +cleo_bias(IHGH)*(Z_POS-REAL(ILOW))  

        P_STNBC(I)=P_STNBC(I) - Z_LEO
!!                            !!!  

!!   !!! NOTICE !!!

!!   Please note that the values in the bias correction table have
!!   opposite sign of correction (the same sign of biases)
!!   Values taken from the tables must be subtracted.
!!   This is opposite sign from operational bias correction.
!!   In operation, values are added to observations.

!!   We decided to use oppsite from operational one because we need 
!!   manual adjustment of correction tables.
!!   It was decided to avoid manual errors (confusion of signs).

!!   Please allow us to use opposite sign from operational 'biascor'.

!!                                                 2000.3.22 K.ONOGI

!!              WRITE(*,'(2A,2X,3F10.2)') &
!!              & 'BIASCOR:',CIDENT(1:5),SPP(I),ST(I),STNBC(I)

      ENDIF
    ENDDO

!!       WRITE(*,'(A,10F7.2,;,/,(10(5X,10F7.2/)))') 'OBS.:',(ST(I),I=1,NST)
!!       WRITE(*,'(A,10F7.2,;,/,(10(5X,10F7.2/)))') 'COR.:',(STNBC(I),I=1,NST)

    LD_LBC=.TRUE.

ENDIF

500 CONTINUE

!!          WRITE(*,*) CIDENT,LBC

!!leo IF (LHOOK) CALL DR_HOOK('BIASCOR_ERA40',1,ZHOOK_HANDLE)
RETURN

965 CONTINUE
 CALL ABOR1('ERROR ON I/O OF STGROUP.T')
966 CONTINUE
 CALL ABOR1('ERROR ON I/O OF CORRECT.T')
967 CONTINUE
 CALL ABOR1('ERROR ON I/O OF COUNTRY.T')
968  CONTINUE
 CALL ABOR1('ERROR ON I/O OF LEOBIASCOR.T')
!!leo IF (LHOOK) CALL DR_HOOK('BIASCOR_ERA40',1,ZHOOK_HANDLE)

END SUBROUTINE BIASCOR_ERA40

