      module lgribread
      contains

      SubROUTINE GRREAD(KEY,YOPER,FARRAY,NI,NJ,NK)

      REAL*8 FARRAY(*)
      CHARACTER*(*) KEY,YOPER
      INTEGER NI,NJ,NK,LLLA,LLLO,URLA,URLO
      character*10 termin
      
   
      CALL GRRE20(KEY,YOPER,FARRAY,NI,NJ,NK,LLLA,LLLO,URLA,URLO,TERMIN)
!!         write(*,*) 'nach grre20'
      RETURN
      END SubROUTINE GRREAD

      SubROUTINE GRRE20(KEY,YOPER,FARRAY,NI,NJ,NK,LLLA,LLLO, URLA,URLO,TERMIN) !

	IMPLICIT NONE
			
      INTEGER LNBLNK
      INTEGER JPACK,JPUNP,J,L,JJ
      INTEGER IWORD,NUMERR,I,IUNIT,IFIELD,ILREC 
      INTEGER NNI,NNJ,NNK,NI,NJ,NK,LLLA,LLLO,URLA,URLO
!!      INTEGER PBCLOSE
 
      PARAMETER (JPACK=1560000)
      PARAMETER (JPUNP=1560000*4)
      INTEGER ISEC0(20)
      INTEGER ISEC1(15120)
      INTEGER ISEC2(15120)
      INTEGER ISEC3(20)
      INTEGER ISEC4(16400),i11
!!      
      REAL*8 ZSEC2(15120)
      REAL*8 ZSEC3(15120)
      REAL   FARRAY(NI*NJ*NK),mean,rms,sig,a
!!      
      CHARACTER*1 YOPER,YOP,swap1,swap2

      character*1 :: INBUF(JPUNP)
	
      integer leni,ierr
      
      CHARACTER*80 KEY
      CHARACTER*10 TERMIN

!      REAL*8 ::  DARRAY(JPACK)
 
   
write(*,*) key
      IUNIT=1
      CALL PBOPEN(IUNIT,KEY(1:LEN(TRIM(KEY))),'R',IERR)
     
      IF (IERR.LT.0) THEN
        WRITE(*,*) ' COULD NOT OPEN GRIB FILE' , KEY
        WRITE(*,*) 'ERRNO: ', IUNIT
        GOTO 9999
      ENDIF
      IFIELD=0
      
50		IFIELD=IFIELD+1
           
               CALL PBGRIB(IUNIT,INBUF,JPACK,LENI,IERR)
        IF (IERR.EQ.0 .and. LENI .gt. 0) THEN
!!          WRITE(*,*) 'GOT FIELD: ',IFIELD, ' LENGTH IS ',LENI
        ELSE IF (IERR.EQ.-1 .OR. LENI .EQ. 0) THEN
!!          WRITE(*,*) IFIELD,' EOF ...'
        	  GOTO 60
        ELSE
          WRITE(*,*) ' ERROR:',IERR
        ENDIF

!!                                                                               
        IF (IERR.NE.0) GOTO 9999
        ILREC = LENI

      IERR = 1
     
!!      WRITE(*,*) 'Vor gribex',ifield
      YOP=YOPER
100   CALL GRIBEX (ISEC0,ISEC1,ISEC2,ZSEC2,ISEC3,ZSEC3,ISEC4, FARRAY((ifield-1)*NI*NJ+1), &
     JPUNP,INBUF,JPACK,IWORD,YOP,IERR)

!!       write(*,*) DARRAY(1:10)
!!       write(*,*) ifield

       if(yoper .eq. 'D') then
!!         print*,ifield,isec4(1)
        l=1
!!        do jj=1,nj
!!         do i=1,ni
!!           FARRAY(i,jj,ifield)=DARRAY(l)
           l=l+1
!!         enddo
!!         write(*,*) jj,farray(1,1,1)
!!        enddo
       endif
      
!!      WRITE(*,*) 'Nach gribex',ifield

      IF (IERR.NE.0) THEN
        NUMERR = NUMERR + 1
      ENDIF

 
      WRITE(TERMIN(1:2),'(I2)') ISEC1(21)-1
      WRITE(TERMIN(3:4),'(I2)') ISEC1(10)
      WRITE(TERMIN(5:6),'(I2)') ISEC1(11)
      WRITE(TERMIN(7:8),'(I2)') ISEC1(12)
      WRITE(TERMIN(9:10),'(I2)') ISEC1(13)
      if(termin(5:5).eq.' ') termin(5:5)='0'
         if(termin(7:7).eq.' ') termin(7:7)='0'
       if(termin(9:9).eq.' ') termin(9:9)='0'

      LLLA=ISEC2(4)
    
      LLLO=ISEC2(5)
      URLA=ISEC2(7)
      URLO=ISEC2(8)

!!**** Print Var-Type of Field-Components: [integer/real]
!!      write(*,*) 'Var-Type = Data-Repr. = ', ISEC4(5)

!!     Print section 0 . Indicator section.

!!		  CALL GRPRS0 (ISEC0)

!!     Print section 1 , product definition section.

!!      CALL GRPRS1 (ISEC0,ISEC1)

!!     Print section 2 , grid definition section.

!!      CALL GRPRS2 (ISEC0,ISEC2,ZSEC2)

!!      write(*,*) farray(FARRAY((IFIELD-1)*ISEC4(1)+1:(IFIELD-1)*ISEC4(1)+10))

      GOTO 50
     	
60		CONTINUE

!!	Mit dieser Abfrage erkennt man, ob es sich um Gitterpunktsdaten
!!	oder Spherical Harmonics handelt.

			IF(ISEC2(1).NE.50) THEN
			  IF(NI.EQ.0.OR.NJ.EQ.0 .OR. YOPER .EQ. "I") THEN
					NI=ISEC2(2)
					NJ=ISEC2(3)
				ELSE
					NNI=ISEC2(2)
					NNJ=ISEC2(3)
					IF((NNI.NE.NI.OR.NNJ.NE.NJ) .AND. YOPER .EQ. "D") THEN
						WRITE(*,*) " FALSE DIMENSIONS: ",NNI,NNJ
!!						STOP
					ENDIF
				ENDIF
				NNK=IFIELD-1
!!	Reduziertes Gitter?
				IF(NI.LT.0) THEN
					NI=0
					DO 222 I=23,23+NJ-1
						NI=NI+ISEC2(I)
222				CONTINUE
          NJ=1		
				ENDIF
			ELSE
				NI=(ISEC2(2)+1)*(ISEC2(2)+2)
				NJ=IFIELD-1
			NNK=1
			ENDIF
			
			IF(NNK.NE.NK.AND.NK.NE.0 .AND. YOPER .EQ. "D") THEN
				WRITE(*,*) " FALSE DIMENSION: ",NNK
				CALL EXIT(1)
			ELSE
			  NK=NNK
			ENDIF
			
			call PBCLOSE(IUNIT,IERR)
!!      write(*,9002) NUMERR
       if(yoper .eq. 'D') then
         write(*,*) ni,nj,nk
         CALL STATIS(NI,NJ,NK,FARRAY,RMS,mean,SIG)
         write(*,*) key,mean,rms,sig
       endif
       RETURN
      
 9002 FORMAT (1H ,'Number of decoding errors = ',I9)
 9004 FORMAT (1H ,'Error return code = ',I4)
 9999 CALL EXIT(1)
      END SubROUTINE GRRE20


      SUBROUTINE GRWRIT(KEY,FARRAY,NI,NJ,NK,NUMERR)
      CHARACTER*(*) KEY
      CHARACTER*10 TERMIN
      INTEGER NUMERR,NI,NJ,NK,I,LLLO,LLLA,URLO,URLA,IPAR
      REAL FARRAY(NI,NJ,NK)
      
      CALL GRWR20(KEY,FARRAY,NI,NJ,NK,NUMERR,LLLA,LLLO,URLA,URLO,IPAR)

      RETURN
      END SUBROUTINE GRWRIT

      SUBROUTINE GRWRSH(KEY,FARRAY,NI,NJ,NK,NUMERR,MM,IPAR)
      
      IMPLICIT NONE
      
      INTEGER JPACK
      INTEGER IWORD,NUMERR,NI,NJ,NK,I,K,IPAR,MM
      CHARACTER*(*) KEY
      CHARACTER*10 TERMIN
      CHARACTER*80 GSFILEN
      REAL*8 FARRAY(NI,NJ,NK)
!!
      LOGICAL LDEBU
      PARAMETER ( LDEBU=.FALSE. )
!!      PARAMETER ( LDEBU=.TRUE. )
!!
!!
!!     Arrays are dimensioned to accommodate T319/N160 data volumes.
!!      INTEGER JPACK
      PARAMETER ( JPACK = 70000 )
!!
!!     This parameter holds the 'number of bytes per word'
      INTEGER JPNBYTE

      PARAMETER( JPNBYTE = 8 )
!!
!!
      REAL*8 ZSEC2(512), ZSEC3(512), ZSEC4(JPACK*JPNBYTE)
      INTEGER ISEC0(2), ISEC1(512), ISEC2(512), ISEC3(2), ISEC4(512)
      INTEGER NMESSAGE, ISTATUS, ILENB, IPUNP
      INTEGER KUNIT, LENGTH,KL,IT,LNBLNK,J,LEN_TRIM
!! 
!!
      CHARACTER*1 INBUFF(JPACK*JPNBYTE)
!!
      CHARACTER*1 YOPER,swap1,swap2
!!
      integer itref,irec
!!     
      integer nstep,ncount
!!         
      write(*,*) key
      KL=LEN_TRIM(KEY)
      write(*,*) kl
      IF(KL.LT.14) WRITE(*,*) 'Invalid File Name ',KEY
      IF(KEY(KL-12:KL-12).EQ.'/') THEN
        IT=KL-7
      ELSE IF(KEY(KL-14:KL-14).EQ.'/' .or. kl .eq. 14) THEN
        IT=KL-9
      ENDIF
      
      IF(IT.EQ.KL-9 ) THEN
        TERMIN=KEY(IT:KL)
      ELSE
        IF(KEY(IT:IT+1).LT.'70') THEN
        TERMIN(1:2)='20'
        TERMIN(3:10)=KEY(IT:KL)
        ELSE
        TERMIN(1:2)='19'
        TERMIN(3:10)=KEY(IT:KL)
        ENDIF
      ENDIF

!!     Clear error counter.
      NUMERR = 0
!!
!!     Set message counter
      NMESSAGE = 0
!!
!!
!!     Read dummy GRIB file to get (default) GRIB headers
!!
      ITREF  = 0
      NSTEP  = 0
      NCOUNT = 0
!!
      IF(LDEBU) THEN
        WRITE(*,*) ' '
        WRITE(*,*) 'QG2grib: Open dummy GRIB file    = ', ITREF
        WRITE(*,*) ' '
      ENDIF
!!
!!     Lengths of INBUFF and PSEC4
      ILENB = JPACK
      IPUNP = JPACK * JPNBYTE
!!
!!     Open data file for reading
      CALL GETENV('GRIB_SH_SAMPLE',GSFILEN)
      CALL PBOPEN( KUNIT,GSFILEN, 'R', ISTATUS)
      IF(LDEBU) WRITE(*,*) 'QG2grib: After PBOPEN, status      = ', ISTATUS
      IF (ISTATUS .NE. 0) STOP 'QG2grib: PBOPEN failed'
      CALL PBGRIB( KUNIT, INBUFF, JPACK*JPNBYTE, LENGTH, ISTATUS)
      IF(LDEBU) THEN
         WRITE(*,*) ' '
!!
        WRITE(*,*) '***********************************************'
        WRITE(*,*) ' '
        WRITE(*,*) 'QG2grib: PBGRIB return status      = ',ISTATUS
      ENDIF
!!
      IF( ISTATUS.EQ.-1) THEN
!!!         
!!$$$         GO TO 222
      END IF

      IF( ISTATUS.LT.-1) STOP 'QG2grib: PBGRIB error reading file'
!!
      NMESSAGE = NMESSAGE+1
      NCOUNT = NCOUNT + 1
!!
      IF(LDEBU) THEN 
        WRITE(*,*) 'QG2grib: Message number            = ', NMESSAGE
        WRITE(*,*) 'QG2grib: Length of message         = ', LENGTH
        WRITE(*,*) 'QG2grib: Time step count           = ', NCOUNT
      ENDIF
!!
!!     'D' function to unpack entire GRIB message.
      YOPER = 'D'

!!
!!     Force a return from GRIBEX even if an error occurs.
      ISTATUS = 1
      CALL GRIBEX( ISEC0, ISEC1, ISEC2, ZSEC2, ISEC3, ZSEC3, &
     ISEC4,ZSEC4, IPUNP, INBUFF, ILENB, IWORD, YOPER, ISTATUS)
!!
!!     Check return code.
      IF(LDEBU) THEN
        WRITE(*,*) 'QG2grib: GRIBEX return code        = ', ISTATUS
        WRITE(*,*) ' '
        WRITE(*,*) '***********************************************'
        IF (ISTATUS.EQ.-6) WRITE(*,*) 'QG2grib: Pseudo-grib data found.'
      ENDIF
!!
!!
      IF (ISTATUS.GT.0) THEN
        NUMERR = NUMERR + 1
        GO TO 300
      ENDIF
!!
!!     Print section 0 and 1.
      IF (LDEBU) CALL GRPRS0( ISEC0)
      IF (LDEBU) CALL GRPRS1( ISEC0, ISEC1)
!!
!!     Print section 2 if present.
      IF(LDEBU) THEN
        IF (ISEC1(5).EQ.0.OR.ISEC1(5).EQ.64) THEN
          WRITE(*,*) ' '
          WRITE(*,*) 'QG2grib: No section 2 in GRIB message.'
        ELSE
          CALL GRPRS2( ISEC0, ISEC2, ZSEC2)
        ENDIF
!!
!!     Print section 3 if present.
        IF (ISEC1(5).EQ.0.OR.ISEC1(5).EQ.128) THEN
          WRITE(*,*) ' '
          WRITE(*,*) 'QG2grib: No section 3 in GRIB message.'
        ELSE
          CALL GRPRS3( ISEC0, ISEC3, ZSEC3)
        ENDIF
      ENDIF
!!
!!     Print section 4.
      IF (LDEBU) CALL GRPRS4( ISEC0, ISEC4, ZSEC4)
!!
!!
      call pbclose(KUNIT,ISTATUS)
!!
      IF(LDEBU) WRITE(*,*)'QG2grib: PBCLOSE return code        = ', ISTATUS

      call pbopen(kunit, KEY, 'w', istatus)
!!
      IF(LDEBU) WRITE(*,*)'QG2grib: PBOPEN return code        = ', ISTATUS
!!
!!
!!     Modify GRIB header for LL files
!!     

      READ(TERMIN(1:2),'(I2)') ISEC1(21)
      READ(TERMIN(3:4),'(I2)') ISEC1(10)
      READ(TERMIN(5:6),'(I2)') ISEC1(11)
      READ(TERMIN(7:8),'(I2)') ISEC1(12)
      READ(TERMIN(9:10),'(I2)') ISEC1(13)
      ISEC1(21)=ISEC1(21)+1


      IF(IPAR.EQ.0) IPAR=135
      ISEC1(6)=IPAR

      ISEC4(1)=(MM+1)*(MM+2)
!!      ISEC4(2)=16 !  24 bit

!!      ISEC2 (12) = 0
!!      ISEC2 (17) = 0
!!
!!
!!     Print section 0 and 1.
      IF (LDEBU) CALL GRPRS0( ISEC0)
      IF (LDEBU) CALL GRPRS1( ISEC0, ISEC1)
!!
!!     Print section 2 if present.
      IF(LDEBU) THEN    
        IF (ISEC1(5).EQ.0.OR.ISEC1(5).EQ.64) THEN
          WRITE(*,*) ' '
          WRITE(*,*) 'QG2grib: No section 2 in GRIB message.'
        ELSE
          CALL GRPRS2( ISEC0, ISEC2, ZSEC2)
        ENDIF
!!
!!     Print section 3 if present.
        IF (ISEC1(5).EQ.0.OR.ISEC1(5).EQ.128) THEN
          WRITE(*,*) ' '
          WRITE(*,*) 'QG2grib: No section 3 in GRIB message.'
        ELSE
          CALL GRPRS3( ISEC0, ISEC3, ZSEC3)
        ENDIF
      ENDIF

      do irec = 1, NK
!!

      ISEC1(8)=irec

      ILENB=ISEC4(1)*ISEC4(2)/8+1000
      YOPER = 'C'
      K=0
      DO 22 J=0,MM
      DO 22 I=I,MM
        ZSEC4(K+1)=FARRAY(K,1,IREC)
        K=K+1
 22   CONTINUE

!!
!!     Print section 4.
      IF (LDEBU) CALL GRPRS4( ISEC0, ISEC4, ZSEC4)

      CALL GRIBEX (ISEC0,ISEC1,ISEC2,ZSEC2,ISEC3,ZSEC3,ISEC4,ZSEC4, &
      NI*NJ*JPNBYTE,INBUFF,ILENB,IWORD,YOPER,ISTATUS)

!!
      call pbwrite (KUNIT, INBUFF, ILENB, ISTATUS)
!!
      IF(LDEBU) THEN
        WRITE(*,*) 'QG2grib: PBWRITE return code        = ', ISTATUS
        WRITE(*,*) 'QG2grib: Record added               = ', irec
      ENDIF
!!
      end do

      call pbclose(KUNIT,ISTATUS)
!!
      IF(LDEBU) WRITE(*,*)'QG2grib: PBCLOSE return code        = ', ISTATUS
!!
 999   continue    
      RETURN

!!     Error reading input file.
  300 STOP 'GRDEMO: Read error'
!!

      END SUBROUTINE GRWRSH

      SUBROUTINE GRWR20(KEY,FARRAY,NI,NJ,NK,NUMERR,LLLA,LLLO,URLA,URLO,IPAR)
      
      IMPLICIT NONE
      
      INTEGER JPACK
      INTEGER IWORD,NUMERR,NI,NJ,NK,I,IPAR
      INTEGER LLLO,LLLA,URLO,URLA
      CHARACTER*(*) KEY
      CHARACTER*10 TERMIN
      CHARACTER*80 GSFILEN
      REAL*8 FARRAY(NI,NJ,NK),MW,RMS,SIG
!!
      LOGICAL LDEBU
      PARAMETER ( LDEBU=.FALSE. )
!!      PARAMETER ( LDEBU=.TRUE. )
!!
!!
!!     Arrays are dimensioned to accommodate T319/N160 data volumes.
!!      INTEGER JPACK
      PARAMETER ( JPACK = 70000 )
!!
!!     This parameter holds the 'number of bytes per word'
      INTEGER JPNBYTE

      PARAMETER( JPNBYTE = 8 )
!!
!!
      REAL*8 ZSEC2(512), ZSEC3(512), ZSEC4(JPACK*JPNBYTE)
      INTEGER ISEC0(2), ISEC1(512), ISEC2(512), ISEC3(2), ISEC4(512)
      INTEGER NMESSAGE, ISTATUS, ILENB, IPUNP
      INTEGER KUNIT, LENGTH,KL,IT,LNBLNK,J,LEN_TRIM
!!
!!
      CHARACTER*1 INBUFF(JPACK*JPNBYTE)
!!
      CHARACTER*1 YOPER,swap1,swap2
!!
      integer itref,irec
!!     
      integer nstep,ncount
!!         
      write(*,*) key
      KL=LEN_TRIM(KEY)
      write(*,*) kl
      IF(KL.LT.14) WRITE(*,*) 'Invalid File Name ',KEY
      IF(KEY(KL-12:KL-12).EQ.'/') THEN
        IT=KL-7
      ELSE IF(KEY(KL-14:KL-14).EQ.'/' .or. kl .eq. 14) THEN
        IT=KL-9
      ENDIF
      
      IF(IT.EQ.KL-9) THEN
        TERMIN=KEY(IT:KL)
      ELSE
        IF(KEY(IT:IT+1).LT.'70') THEN
        TERMIN(1:2)='20'
        TERMIN(3:10)=KEY(IT:KL)
        ELSE
        TERMIN(1:2)='19'
        TERMIN(3:10)=KEY(IT:KL)
        ENDIF
      ENDIF

!!     Clear error counter.
      NUMERR = 0
!!
!!     Set message counter
      NMESSAGE = 0
!!
!!
!!     Read dummy GRIB file to get (default) GRIB headers
!!
      ITREF  = 0
      NSTEP  = 0
      NCOUNT = 0
!!
      IF(LDEBU) THEN
        WRITE(*,*) ' '
        WRITE(*,*) 'QG2grib: Open dummy GRIB file    = ', ITREF
        WRITE(*,*) ' '
      ENDIF
!!
!!     Lengths of INBUFF and PSEC4
      ILENB = JPACK
      IPUNP = JPACK * 8
!!
!!     Open data file for reading
      CALL GETENV('GRIB_SAMPLE',GSFILEN)
      CALL PBOPEN( KUNIT,GSFILEN, 'R', ISTATUS)
      IF(LDEBU) WRITE(*,*) 'QG2grib: After PBOPEN, status      = ', ISTATUS
      IF (ISTATUS .NE. 0) STOP 'QG2grib: PBOPEN failed'
      CALL PBGRIB( KUNIT, INBUFF, JPACK*JPNBYTE, LENGTH, ISTATUS)
      IF(LDEBU) THEN
         WRITE(*,*) ' '
!!
        WRITE(*,*) '***********************************************'
        WRITE(*,*) ' '
        WRITE(*,*) 'QG2grib: PBGRIB return status      = ',ISTATUS
      ENDIF
!!
      IF( ISTATUS.EQ.-1) THEN
!!!         
!!$$$         GO TO 222
      END IF

      IF( ISTATUS.LT.-1) STOP 'QG2grib: PBGRIB error reading file'
!!
      NMESSAGE = NMESSAGE+1
      NCOUNT = NCOUNT + 1
!!
      IF(LDEBU) THEN 
        WRITE(*,*) 'QG2grib: Message number            = ', NMESSAGE
        WRITE(*,*) 'QG2grib: Length of message         = ', LENGTH
        WRITE(*,*) 'QG2grib: Time step count           = ', NCOUNT
      ENDIF
!!
!!     'D' function to unpack entire GRIB message.
      YOPER = 'D'
!!
!!     Force a return from GRIBEX even if an error occurs.
      ISTATUS = 1
      CALL GRIBEX( ISEC0, ISEC1, ISEC2, ZSEC2, ISEC3, ZSEC3, &
      ISEC4,ZSEC4, IPUNP, INBUFF, ILENB, IWORD, YOPER, ISTATUS)
!!
!!     Check return code.
      IF(LDEBU) THEN
        WRITE(*,*) 'QG2grib: GRIBEX return code        = ', ISTATUS
        WRITE(*,*) ' '
        WRITE(*,*) '***********************************************'
        IF (ISTATUS.EQ.-6) WRITE(*,*) 'QG2grib: Pseudo-grib data found.'
      ENDIF
!!
!!
      IF (ISTATUS.GT.0) THEN
        NUMERR = NUMERR + 1
        GO TO 300
      ENDIF
!!
!!     Print section 0 and 1.
      IF (LDEBU) CALL GRPRS0( ISEC0)
      IF (LDEBU) CALL GRPRS1( ISEC0, ISEC1)
!!
!!     Print section 2 if present.
      IF(LDEBU) THEN
        IF (ISEC1(5).EQ.0.OR.ISEC1(5).EQ.64) THEN
          WRITE(*,*) ' '
          WRITE(*,*) 'QG2grib: No section 2 in GRIB message.'
        ELSE
          CALL GRPRS2( ISEC0, ISEC2, ZSEC2)
        ENDIF
!!
!!     Print section 3 if present.
        IF (ISEC1(5).EQ.0.OR.ISEC1(5).EQ.128) THEN
          WRITE(*,*) ' '
          WRITE(*,*) 'QG2grib: No section 3 in GRIB message.'
        ELSE
          CALL GRPRS3( ISEC0, ISEC3, ZSEC3)
        ENDIF
      ENDIF
!!
!!     Print section 4.
      IF (LDEBU) CALL GRPRS4( ISEC0, ISEC4, ZSEC4)
!!
!!
      call pbclose(KUNIT,ISTATUS)
!!
      IF(LDEBU) WRITE(*,*)'QG2grib: PBCLOSE return code        = ', ISTATUS

      call pbopen(kunit, KEY, 'w', istatus)
!!
      IF(LDEBU) WRITE(*,*)'QG2grib: PBOPEN return code        = ', ISTATUS
!!
!!
!!     Modify GRIB header for LL files
!!     
      ISEC2 (1)  = 0
      ISEC2 (2)  = NI  
      ISEC2 (3)  = NJ

      READ(TERMIN(1:2),'(I2)') ISEC1(21)
      READ(TERMIN(3:4),'(I2)') ISEC1(10)
      READ(TERMIN(5:6),'(I2)') ISEC1(11)
      READ(TERMIN(7:8),'(I2)') ISEC1(12)
      READ(TERMIN(9:10),'(I2)') ISEC1(13)
      ISEC1(21)=ISEC1(21)+1
      IF(ISEC1(12) .EQ. 0) ISEC1(12)=1
      IF(ISEC1(11) .EQ. 0) ISEC1(11)=1
      WRITE(*,*) 'ISEC1',ISEC1(1:30)

      if(ni .eq. 35718) then 
      else
        ISEC2(4)=LLLA
        ISEC2(5)=LLLO
        ISEC2(7)=URLA
        ISEC2(8)=URLO
        ISEC2(9)=(URLO-LLLO)/NI
        ISEC2(10)=(URLA-LLLA)/NJ
      endif
      WRITE(*,*) 'ISEC2:',ISEC2(4:10)
      ISEC2 (11) = 64

      IF(IPAR.EQ.0) IPAR=135
      ISEC1(6)=IPAR

      IF(LLLA.EQ.0.AND.URLA.EQ.0) THEN
        IF(NI.EQ.32.AND.NJ.EQ.32) THEN
        
          ISEC2 (4)  = 62890
          ISEC2 (5)  = -9865
!!!!      ISEC2 (6)  = 10000000
          ISEC2 (7)  = 35436
          ISEC2 (8)  = 33771
          ISEC2 (9) = 1406
          ISEC2 (10)  = 871

      ELSE

        
          ISEC2 (4)  = -90000
          ISEC2 (5)  = 1500
!!!!      ISEC2 (6)  = 10000000
          ISEC2 (7)  = 90000
          ISEC2 (8)  = 358500
          ISEC2 (9)  = 3000
          ISEC2 (10) = 3000
 
        ENDIF
      ENDIF

      ISEC4(1)=NI*NJ
!!      ISEC4(2)=24 !  24 bit

!!      ISEC2 (12) = 0
!!      ISEC2 (17) = 0
!!
!!
!!     Print section 0 and 1.
      IF (LDEBU) CALL GRPRS0( ISEC0)
      IF (LDEBU) CALL GRPRS1( ISEC0, ISEC1)
!!
!!     Print section 2 if present.
      IF(LDEBU) THEN    
        IF (ISEC1(5).EQ.0.OR.ISEC1(5).EQ.64) THEN
          WRITE(*,*) ' '
          WRITE(*,*) 'QG2grib: No section 2 in GRIB message.'
        ELSE
          CALL GRPRS2( ISEC0, ISEC2, ZSEC2)
        ENDIF
!!
!!     Print section 3 if present.
        IF (ISEC1(5).EQ.0.OR.ISEC1(5).EQ.128) THEN
          WRITE(*,*) ' '
          WRITE(*,*) 'QG2grib: No section 3 in GRIB message.'
        ELSE
          CALL GRPRS3( ISEC0, ISEC3, ZSEC3)
        ENDIF
      ENDIF
!!
!!     Print section 4.
      IF (LDEBU) CALL GRPRS4( ISEC0, ISEC4, FARRAY)

      do irec = 1, NK
!!

      ISEC1(8)=irec

      ILENB=NI*NJ*ISEC4(2)/8+700
      YOPER = 'C'
      DO 22 J=1,NJ
      DO 22 I=1,NI
        ZSEC4((J-1)*NI+I)=FARRAY(I,J,IREC)
 22   CONTINUE

!!     CALL STATIS(NI,NJ,1,FARRAY(:,:,IREC),RMS,MW,SIG)
!!     WRITE(*,*) 'In grwr20',NI,NJ,NK,mw,rms,sig
!!     CALL STATIS(NI,NJ,1,ZSEC4,RMS,MW,SIG)
!!     WRITE(*,*) 'In grwr20',NI,NJ,NK,mw,rms,sig

      CALL GRIBEX (ISEC0,ISEC1,ISEC2,ZSEC2,ISEC3,ZSEC3,ISEC4, &
      ZSEC4,NI*NJ*8,INBUFF,ILENB,IWORD,YOPER,ISTATUS)

!!     WRITE(*,*) 'ISTATUS',ISTATUS,ILENB

       do 112 i=1,ilenb,4
         swap1=inbuff(i)
         swap2=inbuff(i+1)
 112  continue

!!
      call pbwrite (KUNIT, INBUFF, ILENB, ISTATUS)
!!
      IF(LDEBU) THEN
        WRITE(*,*) 'QG2grib: PBWRITE return code        = ', ISTATUS
        WRITE(*,*) 'QG2grib: Record added               = ', irec
      ENDIF
!!
      end do

      call pbclose(KUNIT,ISTATUS)
!!
      IF(LDEBU) WRITE(*,*)'QG2grib: PBCLOSE return code        = ', ISTATUS
!!
 999   continue     
      RETURN

!!     Error reading input file.
  300 STOP 'GRDEMO: Read error'
!!

      END SUBROUTINE GRWR20

      SUBROUTINE STATIS (NI,NJ,NK,PHI,RMS,MW,SIG)
	

      REAL*8 PHI(NI*NJ*NK)
      REAL*8 SIG,MW,RMS,P
      REAL*8 R,M
 
      N=NI*NJ*NK
 
      R=0.D0
			M=0.D0
 
      DO 10 K=1,N
	      P=PHI(K)
        R=R+P*P
 				M=M+P
10    CONTINUE

      RMS=SQRT(R/N)
	MW=M/N
	
	IF(RMS*RMS-MW*MW.LT.0.) THEN	  
	  SIG=0.0
	ELSE
	  SIG=SQRT(RMS*RMS-MW*MW)
	ENDIF
 
      RETURN
      END subroutine statis

      END module lgribread


