SUBROUTINE PNTERP (ILEV,K_NLEV,PP,POS)

!!**** *PNTERP*

!!       PURPOSE.
!!      ----------

!!           INTERPOLATE POSITION OF PRESSURE OR HEIGHT
!!               WITHIN TABLE OF UPPER AIR LEVELS

!!       INTERFACE.
!!      ------------

!!         CALL PNTERP (ILEV,NLEV,PP,POS)

!!        INPUT
!!         ILEV   -  TABLE OF STANDARD PRESSURE OR HEIGHT LEVELS
!!         NLEV   -  DIMENSION OF ILEV
!!         PP     -  PRESSURE LEVEL

!!        OUTPUT
!!         IERR   -  ERROR CODE
!!         POS    -  INTERPOLATED POSITION WITHIN TABLE

!!       METHOD.
!!      ---------

!!           IF VALUES OF ILEV INCREASE, HEIGHT ASSUMED, LINEAR IN Z
!!           IF VALUES OF ILEV DECREASE, PRESSURE ASSUMED, LINEAR IN LOGP

!!       EXTERNALS.
!!      -----------

!!         NONE

!!       REFERENCE.
!!      ------------

!!         NONE

!!       AUTHOR.
!!      ---------

!!         B. NORRIS,  ECMWF,  JANUARY 1989.

!!       MODIFICATIONS.
!!      ----------------
!!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!!         NONE

IMPLICIT NONE
INTEGER,PARAMETER :: JPIM = 4
INTEGER ,PARAMETER:: JPRM = 4

INTEGER(KIND=JPIM),INTENT(IN)    :: ILEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_NLEV 
REAL(KIND=JPRM)   ,INTENT(IN)    :: PP 
REAL(KIND=JPRM)   ,INTENT(OUT)   :: POS 
REAL(KIND=JPRM) :: Z_EPS 

DIMENSION ILEV(K_NLEV)
LOGICAL :: LL_LDOWN
DATA Z_EPS/0.000001_JPRM/

INTEGER(KIND=JPIM) :: IERR, JLEV, JPOS
REAL(KIND=JPRM) :: ZHOOK_HANDLE

!!    -------------------------------------------------------------------

!!                       1.  SETUP
!!                     ------------

100 CONTINUE
IERR=0
LL_LDOWN=.FALSE.
IF(ILEV(1) > ILEV(K_NLEV)) LL_LDOWN=.TRUE.

!!    -------------------------------------------------------------------

!!           2.  CHECK WHETHER VALUE IS OUTSIDE TABLE LIMITS
!!         --------------------------------------------------

200 CONTINUE
IF(LL_LDOWN) THEN
  IF(PP >= REAL(ILEV(1))) THEN
    POS=1.0_JPRM+Z_EPS
    GO TO 400
  ENDIF
  IF(PP <= REAL(ILEV(K_NLEV))) THEN
    POS=REAL(K_NLEV)-Z_EPS
    GO TO 400
  ENDIF
ELSE
  IF(PP <= REAL(ILEV(1))) THEN
    POS=1.0_JPRM+Z_EPS
    GO TO 400
  ENDIF
  IF(PP >= REAL(ILEV(K_NLEV))) THEN
    POS=REAL(K_NLEV)-Z_EPS
    GO TO 400
  ENDIF
ENDIF

!!    ------------------------------------------------------------------

!!          3.   VALUE IS BETWEEN MAX AND MIN TABLE VALUES
!!         ------------------------------------------------

300 CONTINUE

!!         3.2   VALUE IS BETWEEN 2 TABLE VALUES

DO 322 JLEV=1,K_NLEV-1
IF(LL_LDOWN) THEN
  IF(PP <= REAL(ILEV(JLEV)).AND.PP >= REAL(ILEV(JLEV+1)))THEN
    JPOS=JLEV
    GO TO 324
  ENDIF
ELSE
  IF(PP >= REAL(ILEV(JLEV)).AND.PP <= REAL(ILEV(JLEV+1)))THEN
    JPOS=JLEV
    GO TO 324
  ENDIF
ENDIF
322 CONTINUE
324 CONTINUE

!!           3.3   INTERPOLATE

IF(LL_LDOWN) THEN
  POS=REAL(JPOS)+&
   & (LOG(REAL(ILEV(JPOS)))-LOG(PP))/&
   & (LOG(REAL(ILEV(JPOS)))-LOG(REAL(ILEV(JPOS+1))))  
ELSE
  POS=REAL(JPOS)+&
   & (REAL(ILEV(JPOS))-PP)/&
   & (REAL(ILEV(JPOS))-REAL(ILEV(JPOS+1)))  
ENDIF

!!    ------------------------------------------------------------------

!!                   4.   EXIT
!!                  -----------

400 CONTINUE
END SUBROUTINE PNTERP
