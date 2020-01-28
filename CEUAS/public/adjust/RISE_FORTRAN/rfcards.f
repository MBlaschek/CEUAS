      module rfcards

      integer, parameter :: JPRB=4

      contains


      SUBROUTINE readsounding(infileunit,EOFflag,msgfile,addrec, 
     +  bar,barn,cardsnum,sonde,clouds,corp,corrh,cort,cortd,corw,corz,
     +  curday,curhour,curmon,curyear,datas,dpdp,elev,etime,hgt,
     +  lat,lon,lqual,mlvls,obt,press,qc,qd,qe,qp,qrh,qt,qw,qz,
     +  rh,temp,tlvl,wdir,wspd)

!! Purpose: To read metadata and data for one sounding and store them in 
!!          appropriate variables.
!!          When a header record is found, parse it if
!!          possible and read and parse the subsequent data lines, 
!!          if possible. 

!!          The program is meant to be used as a subroutine when 
!!          reading multiple soundings. (As written, the program 
!!          reads the first sounding in a data file.)
!!
!! Version 1.0 - Imke Durre - 27 August 2001



!! Sounding variables
      INTEGER addrec      !! Number of additional sounding records 
                          !! for the same sounding; introduced when 
                          !! the number of levels per sounding was 
                          !! limited to 175 (from header)
      INTEGER bar         !! Code specifying the type of number 
                          !! identifying the sonde (from header)
      CHARACTER*20 barn   !! Sonde serial number (from header)
      INTEGER*8 cardsnum    !! WMO station code plus CARDS digit (header)
      CHARACTER*9 clouds  !! Code giving weather and cloud
                          !! observations (from header)
      INTEGER corrh       !! Code for type of correction applied to 
                          !! relative humidity values (from header)
      INTEGER corp        !! Code for type of correction applied to 
                          !! pressure values (from header)
      INTEGER cort        !! Code for type of correction applied to
                          !! temperature values (from header)
      INTEGER cortd       !! Code for type of correction applied to 
                          !! dewpoint values (from header)
      INTEGER corw        !! Code for type of correction applied to 
                          !! direction and speed values (from header)
      INTEGER corz        !! Code for type of correction applied to 
                          !! height values (from header)
      INTEGER curyear     !! year of current record (from header)
      INTEGER curmon      !! month of current record (from header)
      INTEGER curday      !! day of current record (from header)
      INTEGER curhour     !! hour of current record (from header)
      INTEGER datas       !! Data source (from header)
      INTEGER dpdp(500)   !! Dewpoint depression (deg C * 10)
                          !! at each level (from data records)
      INTEGER*8 elev      !! Site elevation (in m*10; from header)
      INTEGER etime(500)  !! Time elapsed since the release of the
                          !! sonde in whole minutes and seconds (mmmss)
                          !! at each level (from data records)
      INTEGER*8 hgt(500)    !! Geopotential height (m) at each level 
                          !! (from data records)
      INTEGER*8 lat       !! Latitude (in degrees*10**5; from header)
      INTEGER*8 lon       !! Longitude (in degrees*10**5; from header)
      INTEGER lqual(500)  !! Level quality indicator at each level
                          !! (from data records)
      INTEGER mlvls       !! # of levels in sounding (from header)
      INTEGER obt         !! Code for the type of observation (header)
      INTEGER*8 press(500)  !! Pressure (Pa) at each level (data records)
      INTEGER qc          !! Code referencing mode of decoding 
                          !! quality flags (from header)
      INTEGER qdum(500)   !! "Internal NCDC use" field at each level
                          !! (from data records)
      INTEGER qd(500)     !! Quality code for dewpoint depression
                          !! at each level (from data records)
      INTEGER qe(500)     !! Quality code for elapsed time at
                          !! each level (from data records)
      INTEGER qp(500)     !! Quality code for pressure at each level
                          !! (from data records)
      INTEGER qrh(500)     !! Quality code for humidity at each level
                          !! (from data records)
      INTEGER qt(500)     !! Quality code for temperature at each level
                          !! (from data records)
      INTEGER qw(500)     !! Quality code for wind values at each level
                          !! (from data records)
      INTEGER qz(500)     !! Quality code for height at each level
                          !! (from data records)
      INTEGER rh(500)     !! Relative humidity (%*10) at each level 
                          !! (from data records)
      INTEGER rtime       !! Actual time of release of sonde (hour and
                          !! minutes UTC) (from header)
      INTEGER sonde       !! Type of sonde used for the current flight
                          !! (from header)
      INTEGER stni        !! Station number indicator (from header)
      CHARACTER stnnum*8  !! Station number of the type indicated 
                          !! by stni (from header)
      INTEGER temp(500)   !! Temperature (deg C * 10) at each level
                          !! (from data records)
      INTEGER tlvl(500)   !! level type indicator at each level
                          !! (from data records)
      INTEGER wdir(500)   !! Wind direction (deg) at each level 
                          !! (from data records)
      INTEGER wspd(500)   !! Wind speed (m/s * 10) at each level
                          !! (from data records)

!! Other variables
      CHARACTER datalvl*56  !! Current data record 
      CHARACTER record*120  !! Current header record
      INTEGER err           !! Error status code for header
      INTEGER infileunit   !! unit number of input data file 
      INTEGER irec          !! No. characters in current record
      INTEGER istat         !! I/O status code
      INTEGER letrim       !! INTEGER FUNCTION letrim
      INTEGER lvl           !! Index of the current data record in the sounding
      INTEGER msgfile       !! unit number of messages file
      LOGICAL EOFflag       !! .TRUE. if end of file reached
      LOGICAL header       !! .TRUE. if current record is a header 

      EOFflag = .FALSE.
!! Set header and data variables to missing
      CALL set_header_missing(addrec,bar,barn,cardsnum,clouds,
     +  corp,corrh,cort,cortd,corw,corz,curday,curhour,curmon,curyear,
     +  datas,elev,lat,lon,mlvls,obt,qc,rtime,sonde,stni,stnnum)
      DO lvl = 1,500
        CALL set_data_missing(dpdp(lvl),etime(lvl),hgt(lvl),
     +    lqual(lvl),press(lvl),qd(lvl),qdum(lvl),qe(lvl),qp(lvl),
     +    qrh(lvl),qt(lvl),qw(lvl),qz(lvl),rh(lvl),temp(lvl),tlvl(lvl),
     +    wdir(lvl),wspd(lvl))
      ENDDO

!! Read one record
      READ(infileunit,100,IOSTAT=istat) record

      IF (istat .gt. 0) THEN !! Record is unreadable.
!! Write appropriate error message to messages file.
         WRITE(msgfile,200)
      ELSE IF (istat .lt. 0) THEN !! End of file has been reached.
         EOFflag = .TRUE.
      ELSE   !! Header record is readable
         IF (INDEX(record,ACHAR(13)) .ne. 0) THEN
!! If record contains a carriage return character (ASCII code 13), 
!! write the appropriate error message to the messages file and stop.
            WRITE(msgfile,300)
            WRITE(msgfile,'(a)') ACHAR(7)
            STOP
         ENDIF

         irec = letrim(record)  !! length of record minus trailing blanks
         IF (header(irec,record(1:1))) THEN
!! If the record is a valid header record, parse it.
            CALL parse_header(record,err,addrec,bar,barn,cardsnum,
     +        clouds,corp,corrh,cort,cortd,corw,corz,curday,curhour,
     +        curmon,curyear,datas,elev,lat,lon,mlvls,obt,qc,rtime,
     +        sonde,stni,stnnum)

            IF (err .eq. 1) then   !! header only partially parsable
              WRITE(msgfile,400) curyear, curmon, curday, curhour,record
              CALL set_header_missing(addrec,bar,barn,cardsnum,clouds,
     +          corp,corrh,cort,cortd,corw,corz,curday,curhour,curmon,
     +          curyear,datas,elev,lat,lon,mlvls,obt,qc,rtime,
     +          sonde,stni,stnnum)
            elseIF (err .eq. 2) then   !! header unparsable
              WRITE(msgfile,500) record
              CALL set_header_missing(addrec,bar,barn,cardsnum,clouds,
     +          corp,corrh,cort,cortd,corw,corz,curday,curhour,curmon,
     +          curyear,datas,elev,lat,lon,mlvls,obt,qc,rtime,
     +          sonde,stni,stnnum)
            ELSEIF (err .eq. 0) then  !! if header parsable

              IF (mlvls .gt. 500) WRITE(msgfile,600) curyear,curmon,
     +          curday,curhour,record
              DO lvl=1,mlvls  !! For each level in the sounding
!! Read a data record from the file
                READ (infileunit,700,IOSTAT=istat) datalvl
!! Try to parse the data record
                CALL parse_data(record,datalvl,lvl,msgfile,istat,
     +            curday,curhour,curmon,curyear,dpdp(lvl),etime(lvl),
     +            hgt(lvl),lqual(lvl),press(lvl),qd(lvl),qdum(lvl),
     +            qe(lvl),qp(lvl),qrh(lvl),qt(lvl),qw(lvl),qz(lvl),
     +            rh(lvl),temp(lvl),tlvl(lvl),wdir(lvl),wspd(lvl))
              ENDDO
           ENDIF
        ELSE 
          WRITE(msgfile,800) record
          CALL set_header_missing(addrec,bar,barn,cardsnum,clouds,
     +      corp,corrh,cort,cortd,corw,corz,curday,curhour,curmon,
     +      curyear,datas,elev,lat,lon,mlvls,obt,qc,rtime,sonde,
     +      stni,stnnum)
        ENDIF
      ENDIF

 100  FORMAT(A120)
 200  FORMAT('Header record cannot be read from file')
 300  FORMAT(9x,'Error: Sounding file contains carriage returns')
 400   FORMAT('Error in header record, sounding skipped: ',
     +   i5,3i3,/,a108)
 500  FORMAT('Unparsable header record, sounding skipped:',/,a108)
 600  FORMAT('Number of levels exceeds bounds of sounding arrays',
     +  i5,3i3,/,a108)
 700  FORMAT(A56)
 800  FORMAT('Invalid header record:',/,a108)

      END SUBROUTINE readsounding

!!--------------------------------------------------------------------

      SUBROUTINE set_header_missing(addrec,bar,barn,cardsnum,clouds,
     +  corp,corrh,cort,cortd,corw,corz,curday,curhour,curmon,curyear,
     +  datas,elev,lat,lon,mlvls,obt,qc,rtime,sonde,stni,stnnum)

!! Purpose: To set all header variables to missing.

!! Arguments
      INTEGER addrec       !! Number of additional sounding records 
                           !! for the same sounding; introduced when 
                           !! the number of levels per sounding was 
                           !! limited to 175
      INTEGER bar          !! Code specifying the type of number
                           !! identifying the sonde 
      CHARACTER*20 barn    !! Sonde serial number
      INTEGER*8 cardsnum    !! WMO station code plus CARDS digit
      CHARACTER*9 clouds   !! Code giving weather and cloud
                           !! observations
      INTEGER corp         !! Code for type of correction applied to 
                           !! pressure values
      INTEGER corrh         !! Code for type of correction applied to 
                           !! relative humidity values
      INTEGER cort         !! Code for type of correction applied to
                           !! temperature values
      INTEGER cortd        !! Code for type of correction applied to 
                           !! dewpoint values
      INTEGER corw         !! Code for type of correction applied to 
                           !! wind direction and speed values
      INTEGER corz         !! Code for type of correction applied to
                           !! height values
      INTEGER curyear      !! year of current record
      INTEGER curmon       !! month of current record
      INTEGER curday       !! day of current record
      INTEGER curhour      !! hour of current record
      INTEGER datas        !! Data source
      INTEGER*8 elev       !! Site elevation (in m*10)
      CHARACTER*1 id       !! First character in sounding header
      INTEGER*8 lat        !! Latitude (in degrees*10**5)
      INTEGER*8 lon        !! Longitude (in degrees*10**5)
      INTEGER mlvls        !! # of levels in sounding
      INTEGER obt          !! Code for the type of observation
      INTEGER qc           !! Code referencing mode of decoding quality flags
      INTEGER rtime        !! Actual time of release of sonde (hour and
                           !! minutes UTC)
      INTEGER sonde        !! Type of sonde used for the current flight
      INTEGER stni         !! Number of the station type
      CHARACTER stnnum*8   !! Station number

      id = ' '
      cardsnum = 999999
      stni = 9
      stnnum = '99999999'
      lat = 9999999
      lon = 99999999
      elev = 99999
      curyear = 9999
      curmon = 99
      curday = 99
      curhour = 99
      rtime = 9999
      clouds = '999999999'
      obt = 99
      bar = 9
      barn = '99999999999999999999'
      sonde = 999
      qc = 9
      datas = 99
      corp = 99
      corrh = 99
      cort = 99
      cortd = 99
      corw = 99
      corz = 99

!! For convenience, addrec and mlvls are set to zero.
      addrec = 0
      mlvls = 0

      END SUBROUTINE

!!--------------------------------------------------------------------

      SUBROUTINE set_data_missing(curdpdp,curetime,curhgt,curlqual,
     +  curpress,curqd,curqdum,curqe,curqp,curqrh,curqt,curqw,curqz,
     +  currh,curtemp,curtlvl,curwdir,curwspd)

!! Purpose: To set values in the data variable arrays at the current
!!          level in the sounding to missing.

!! Arguments
      INTEGER curdpdp   !! Dewpoint depression (deg C * 10)
                        !! at current level
      INTEGER curetime  !! Time elapsed since the release of the
                        !! sonde in whole minutes and seconds (mmmss)
                        !! at current level
      INTEGER*8 curhgt    !! Geopotential height of current level (m)
      INTEGER curlqual  !! Level quality indicator at current level
      INTEGER*8 curpress  !! Pressure at current level (Pa)
      INTEGER curqd     !! Quality code for dewpoint-depression 
                        !! at current level
      INTEGER curqdum   !! "Internal NCDC use" field at current level
      INTEGER curqe     !! Quality code for elapsed time value
                        !! at current level
      INTEGER curqp     !! Quality code for pressure at current level
      INTEGER curqrh     !! Quality code for humidity at current level
      INTEGER curqt     !! Quality code for temperature at current level
      INTEGER curqw     !! Quality code for wind speed and direction
                        !! at current level
      INTEGER curqz     !! Quality code for height at current level
      INTEGER currh     !! Relative humidity at current level (%*10)
      INTEGER curtemp   !! Temperature at current level (deg C * 10)
      INTEGER curtlvl   !! Level type indicator at current level
      INTEGER curwdir   !! Wind direction at current level (deg)
      INTEGER curwspd   !! Wind speed at current level (m/s * 10)

      curlqual = 9
      curetime = 99999
      curpress = 999999
      curhgt = -999999
      curtemp = 9999
      currh = 9999
      curdpdp = 999
      curwdir = 999
      curwspd = 9999
      curtlvl = 99
      curqe = 9 
      curqp = 9
      curqz = 9
      curqt = 9
      curqrh = 9
      curqd = 9
      curqw = 9
      curqdum = 99

      end SUBROUTINE

!!--------------------------------------------------------------------

      SUBROUTINE parse_header(record,parse_err,addrec,bar,barn,
     +  cardsnum,clouds,corp,corrh,cort,cortd,corw,corz,curday,curhour,
     +  curmon,curyear,datas,elev,lat,lon,mlvls,obt,qc,rtime,
     +  sonde,stni,stnnum)

!! Purpose: To interpret a sounding header record.  Two attempts are
!!          made to interpret the header.  The first tries to parse
!!          the entire header (108 characters).  This parsing must be
!!          successful if the sounding is to be used.  If it is not,
!!          then the second attempt tries to parse through character
!!          48, which gives the date and time of the sounding.  If
!!          successful, the date and time of the sounding are written
!!          to the error file.  If it is not, a message is written to
!!          the error file that an unidentifiable sounding has been
!!          encountered.  (Consult NCDC's TDF63 format documentation
!!          for further details on interpreting a sounding header).

!! Arguments required as input:
      CHARACTER*(*) record  !! Header record of current sounding

!! Arguments parsed from header record and returned to main program:
      INTEGER addrec       !! Number of additional sounding records 
                           !! for the same sounding; introduced when 
                           !! the number of levels per sounding was 
                           !! limited to 175
      INTEGER bar          !! Code specifying the type of number
                           !! identifying the sonde 
      CHARACTER*20 barn    !! Sonde serial number
      INTEGER*8 cardsnum    !! WMO station code plus CARDS digit
      CHARACTER*9 clouds   !! Code giving weather and cloud
                           !! observations
      INTEGER corp         !! Code for type of correction applied to 
                           !! pressure values
      INTEGER corrh         !! Code for type of correction applied to 
                           !! relative humidity values
      INTEGER cort         !! Code for type of correction applied to
                           !! temperature values
      INTEGER cortd        !! Code for type of correction applied to 
                           !! dewpoint values
      INTEGER corw         !! Code for type of correction applied to 
                           !! wind direction and speed values
      INTEGER corz         !! Code for type of correction applied to
                           !! height values
      INTEGER curyear      !! year of current record
      INTEGER curmon       !! month of current record
      INTEGER curday       !! day of current record
      INTEGER curhour      !! hour of current record
      INTEGER datas        !! Data source
      INTEGER*8 elev       !! Site elevation (tenths of a meter, MSL)
      CHARACTER*1 id       !! First character in sounding header
      INTEGER*8 lat        !! Latitude (in degrees*10**5, signed)
      INTEGER*8 lon        !! Longitude (in degrees*10**5, signed)
      INTEGER mlvls        !! # of levels in sounding
      INTEGER obt          !! Code for the type of observation
      INTEGER qc           !! Code referencing mode of decoding quality flags
      INTEGER rtime        !! Actual time of release of sonde (hour and
                           !! minutes UTC)
      INTEGER sonde        !! Type of sonde used for the current flight
      INTEGER stni         !! Number of the station type
      CHARACTER stnnum*8   !! Station number

!! Additional arguments
      INTEGER parse_err    !! Error status of parsing header

!! Local Variables
      INTEGER hemi_sign    !! Algebraic sign associated with a hemisphere
                           !! (N and E are +; S and W are -)
      INTEGER istat        !! Read error status code
      INTEGER*8 latu       !! Site latitude, (100000ths of a deg, unsigned)
      CHARACTER*1 latns    !! Hemisphere of sounding station
                           !! (north, south)
      INTEGER*8 lonu       !! Site longitude (100000ths of a deg, unsigned)
      CHARACTER*1 longew   !! Hemisphere of sounding station
                           !! (east, west)


!! The incoming header is contained in a 108-character record.
!! Treating "record" as an internal file, the first 108 characters are
!! parsed into 29 values.
      READ(record,100,IOSTAT=istat) id, cardsnum, stni, stnnum, latu, 
     +  latns, lonu, longew, elev, curyear, curmon, curday, curhour, 
     +  rtime, clouds, obt, bar, barn, sonde, qc, datas, corp, corz, 
     +  cort, corrh, cortd, corw, addrec, mlvls

      IF (istat .eq. 0) THEN  !! parsing successful
        parse_err = 0

!! if latns = 'S' (latitude in Southern Hemisphere), 
!! change the latitude to a negative number.
        IF (latns .eq. 'N') THEN
          hemi_sign = 1.0
        ELSE
          hemi_sign = -1.0
        ENDIF
        lat = hemi_sign*latu

!! if lonew = 'W' (longitude in Western Hemisphere), 
!! change the longitude to a negative number.
        IF (longew .eq. 'E') THEN
          hemi_sign = 1.0
        ELSE
          hemi_sign = -1.0
        ENDIF
        lon = hemi_sign*lonu

      ELSE  !! if first parsing attempt unsuccessful

!! If the full header cannot be parsed, a second attempt is made to
!! parse it up through the date variables.  If this is successful,
!! then the date of the erroneous sounding can be written to the error
!! file.  Otherwise, a message can be written to the error file that an
!! unidentifiable sounding is erroneous (see main program).

         READ(record,200,IOSTAT=istat) id, cardsnum, stni, stnnum, 
     +     latu, latns, lonu, longew, elev, curyear, curmon, curday, 
     +     curhour
         IF (istat .eq. 0) THEN  !! partial parsing successful
            parse_err = 1
         ELSE  !! partial parsing unsuccessful
            parse_err = 2
         ENDIF
      ENDIF

 100  FORMAT(a1,i6,i1,a8,i7,a1,i8,a1,i5,i4,3i2,i4,a9,i2,i1,a20,i3,i1, 
     +    7i2,2i3)
 200  FORMAT(a1,i6,i1,a8,i7,a1,i8,a1,i5,i4,3i2)

      END SUBROUTINE

!!--------------------------------------------------------------------

      SUBROUTINE parse_data(header,datalvl,lvl,msgfile,istat,curday,
     +  curhour,curmon,curyear,curdpdp,curetime,curhgt,curlqual,
     +  curpress,curqd,curqdum,curqe,curqp,curqrh,curqt,curqw,curqz,
     +  currh,curtemp,curtlvl,curwdir,curwspd)

!! Purpose: To parse the data at one level of one sounding and 
!!          store them in the appropriate arrays

!! Arguments required as input from main program:
      CHARACTER*(*) datalvl  !! Line of data in current sounding
                             !! (passed as input from main program)
      INTEGER curyear        !! year of current record (from header)
      INTEGER curmon         !! month of current record (from header)
      INTEGER curday         !! day of current record (from header)
      INTEGER curhour        !! hour of current record (from header)
      CHARACTER*(*) header   !! header record of sounding
      INTEGER istat          !! Error status code of reading data
                             !! record from file 
      INTEGER lvl            !! Index of the current data record 
                             !! in the sounding
      INTEGER msgfile        !! Unit number of messages file

!! Arguments read from data record
      INTEGER curdpdp   !! Dewpoint depression (deg C * 10)
                        !! at current level
      INTEGER curetime  !! Time elapsed since the release of the
                        !! sonde in whole minutes and seconds (mmmss)
                        !! at current level
      INTEGER*8 curhgt    !! Geopotential height of current level (m)
      INTEGER curlqual  !! Level quality indicator at current level
      INTEGER*8 curpress  !! Pressure at current level (Pa)
      INTEGER curqd     !! Quality code for dewpoint-depression 
                        !! at current level
      INTEGER curqdum   !! "Internal NCDC use" field at current level
      INTEGER curqe     !! Quality code for elapsed time value
                        !! at current level
      INTEGER curqp     !! Quality code for pressure at current level
      INTEGER curqrh    !! Quality code for humidity at current level
      INTEGER curqt     !! Quality code for temperature at current level
      INTEGER curqw     !! Quality code for wind speed and direction
                        !! at current level
      INTEGER curqz     !! Quality code for height at current level
      INTEGER currh     !! Relative humidity at current level (%*10)
      INTEGER curtemp   !! Temperature at current level (deg C * 10)
      INTEGER curtlvl   !! Level type indicator at current level
      INTEGER curwdir   !! Wind direction at current level (deg)
      INTEGER curwspd   !! Wind speed at current level (m/s * 10)

!! Local variables
      INTEGER irec        !! Length of string datalvl
      INTEGER jstat       !! Read error status code from reading data 
                          !! record as internal file
      INTEGER letrim     !! Function letrim
      LOGICAL datarecord  !! .TRUE. if record is a data record


!! Set all data variables to missing
      irec = letrim(datalvl)  !! length of datalvl
      IF (istat .gt. 0) THEN !! Data record unreadable.
!! Write appropriate error message to messages file.
        WRITE(msgfile,200) curyear,curmon,curday,curhour,lvl,
     +    header, datalvl
      ELSE IF (istat .lt. 0) THEN !! End of file reached.
!! Write appropriate error message to messages file.
        WRITE(msgfile,300) curyear,curmon,curday,curhour,lvl,header
      ELSE !! Data record readable.
        IF (datarecord(irec,datalvl(1:1))) then
!! If the data record is valid, try to parse it by reading from 
!! it as an internal file.
          READ (datalvl,400,IOSTAT=jstat) curlqual, curetime, 
     +      curpress, curhgt, curtemp, currh, curdpdp, curwdir,
     +      curwspd, curtlvl, curqe, curqp, curqz, curqt, curqrh,
     +      curqd, curqw, curqdum
          IF (jstat .gt. 0) THEN !! Data record unparsable.
            WRITE(msgfile,500) curyear, curmon, curday, curhour, lvl,
     +        header, datalvl
            CALL set_data_missing(curdpdp,curetime,curhgt,curlqual,
     +        curpress,curqd,curqdum,curqe,curqp,curqrh,curqt,
     +        curqw,curqz,currh,curtemp,curtlvl,curwdir,curwspd)
          ENDIF
        ELSE  !! Data record not valid
!! Write appropriate error message to messages file.
          IF (datalvl(1:1) .eq. '#') THEN   
            WRITE(msgfile,600) curyear, curmon, curday, curhour,
     +        lvl, header, datalvl
          ELSE
            WRITE(msgfile,700) curyear, curmon, curday, curhour,
     +        lvl, header, datalvl
          ENDIF
        ENDIF
      ENDIF

 200  FORMAT('Data record cannot be read from file: ',i5,3i3,i4,
     +  /,a108,/,a56)
 300  FORMAT('End of file reached while reading data records',
     +  i5,3i3,i4,/,a108)
 400  FORMAT(i1,i5,i6,i7,i5,i4,i3,i3,i4,i2,8i2)
 500  FORMAT('Data record unparsable: ',i5,3i3,i4,/,a108,/,a56)
 600  FORMAT('Header record reached prematurely:',i5,3i3,i4,
     +  /,a108,/,a56)
 700  FORMAT('Data record invalid: ',i5,3i3,i4,/,a108,/,a56)

      END SUBROUTINE

!!--------------------------------------------------------------------

      end module rfcards

      INTEGER FUNCTION letrim(cstring)

!! Purpose: To find the length of a character string minus any 
!!          trailing blanks.

      character*(*) cstring

      letrim = len(cstring)
      DO WHILE (cstring(letrim:letrim) .eq. ' ')
        letrim = letrim-1
      ENDDO

      end function

!!--------------------------------------------------------------------


      LOGICAL FUNCTION datarecord(irec,marker)

!! Purpose: To determine whether the current record is a valid data
!!          record.  To be valid it must begin with a 0, 1, 2, 8, or 9
!!          and contain 56 characters.

!! Arguments
      INTEGER irec        !! Number of characters in the current record
      CHARACTER*1 marker  !! The first character of the record

      IF (irec .eq. 56 .AND. marker .ne. '#') then
        datarecord =  .TRUE.
      ELSE
        datarecord = .FALSE.
      ENDIF

      END function datarecord

!!--------------------------------------------------------------------

 
      LOGICAL  FUNCTION header(irec,marker)

!! Purpose: To determine whether the current record is a valid sounding
!!          header.  To be valid it must begin with the # character and
!!          contain 108 characters.

!! Arguments
      CHARACTER*1 marker   !! The '#' character
      INTEGER irec   !! Number of characters in the current record

      IF (irec .eq. 108 .AND. marker .eq. '#') THEN
        header =  .TRUE.
      ELSE
        header = .FALSE.
      ENDIF

      END  function

!!-----------------------------------------------------------------------
