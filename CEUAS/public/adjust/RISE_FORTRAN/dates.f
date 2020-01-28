      SUBROUTINE DATES(Year,Month,Day,Time,Tdiff,Year1,Month1,Day1,
     +                       Time1)
C
C**********************************************************************
C
C  Name        :  Dates
C
C  Purpose     :  Computes a new date(Year,Month,Day and Hour)
C                 from a given date and time-difference in
C                 hours ( Positive or negative>

C  Description :
C
C  Restrictions:
C
C
C  Environment : Type          Name         Comment
C
C                Parameter     Year         In Integer
C                Parameter     Month        In Integer
C                Parameter     Day          In Integer
C                Parameter     Time         In Integer
C                Parameter     Tdiff        In Integer
C                Parameter     Year1        Out Integer
C                Parameter     Month1       Out Integer
C                Parameter     Day1         Out Integer
C                Parameter     Time1        Out Integer
C                Intr. Func.   MOD          Fortran
C                Intr. Func.   IABS         Fortran
C
C  History     : Created   :  Bo Lindgren
C
C                Docum.    :  Karin Kirkland    811102
C                             Bo Strandberg     910830
C
C                Rev.      :
C
C**********************************************************************
C
      IMPLICIT NONE

      INTEGER Year,Month,Day,Time,Tdiff,Year1,Month1,Day1,Time1
      INTEGER Ileap,Idelt,Itil,Iad,Ix,Idmax
      DIMENSION Idmax(12)
      DATA Idmax/31,28,31,30,31,30,31,31,30,31,30,31/

      IF(MOD(Year,4).EQ.0.AND..NOT.MOD(Year,100).EQ.0.OR.MOD(Year,400)
     +.EQ.0) GOTO 100
      GOTO 600

  100 IF(Month .NE. 2) GOTO 610
      IF(Day .EQ. 29) GOTO 120
  110 IF(Tdiff) 600,120,120
  120 Ileap = 1
      GOTO 299

  610 IF(Month .NE. 3) GOTO 600
  300 IF(Tdiff) 120,600,600
  600 Ileap= 0
  299 Idelt = Tdiff + Time
      IF(Idelt) 620,130,310
  620 IF(MOD(IABS(Idelt),24) .EQ. 0) GOTO 310
      Itil = 1
      GOTO 111

  310 Itil = 0
  111 Iad = Idelt/24-Itil
      GOTO 320

  130 Time1 = 0
      Day1 = Day
  132 Month1 = Month
      GOTO 155

  320 Time1 = MOD(Idelt,24)
      IF(Time1) 321,322,322
  321 Time1 = Time1+24
  322 Day1 = Day+Iad
  330 Ix= Idmax(Month)+Ileap-Day1
      IF(Ix) 630,140,140
  630 IF(Month-12) 635,640,635
  635 Month1 = Month+1
      Year1 = Year
      GOTO 650

  640 Month1 = 1
      Year1 = Year+1
  650 Day1 = -Ix
      GOTO 160

  140 IF(Day1-1) 141,132,132
  141 IF(Month-1) 145,150,145
  145 Day1 = Day1+Idmax(Month-1)+Ileap
      Month1 = Month-1
  155 Year1 = Year
      GOTO 160

  150 Month1 = 12
      Year1 = Year-1
      Day1 = IDMAX(12)+Day1
  160 RETURN
c
      END
