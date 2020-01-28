MODULE ERROR_MANAGEMENT

contains


subroutine error(routine_name,status,severity,file_name)
!severity 0 => simple possible mistake
!severity 1 => WARING, but the program can still run
!severity 2 => ERROR , the program must be stopped
!file_name or error_message
use ifcore
IMPLICIT NONE
CHARACTER (len = 100) :: routine_name
CHARACTER (len = 200),OPTIONAL :: file_name
INTEGER :: status,severity
CHARACTER (len =200)::cform_1,cform_2,cform, message
INTEGER :: dim_cform_1, dim_cform_2


!If status == 0 I don't have to check anything
!but I could print a good message
if ((status == 0) .and. present(file_name))then
   if( len_trim(file_name) > 1 .and. severity == 0 )then
      
      call format(routine_name,cform_1,dim_cform_1)
      call format(file_name,cform_2, dim_cform_2)
      !write(*,*)"cform_1",trim(cform_1),"!"
      !write(*,*)"cform_2",trim(cform_2),"!"
      !write(*,*)"len_trim(cform_1)",len_trim(cform_1)
      !write(*,*)'(' //trim(cform_1(1:dim_cform_1)) //',A2,'// trim(cform_2(1:dim_cform_2))// ')'
      cform='(' //trim(cform_1(1:dim_cform_1)) //',A2,'// trim(cform_2(1:dim_cform_2))// ')'
      write(*,cform)trim(routine_name),'->',file_name
   endif
elseif(status /= 0) then
   sever: if(severity == 0) then
      message = "Possible incongruence in the program, check ->"
      call format(message,cform_1, dim_cform_1)
      call format(routine_name,cform_2,dim_cform_2)
      cform='('//trim(cform_1(1:dim_cform_1))//','//trim(cform_2(1:dim_cform_2))//')'
      WRITE(*,cform)trim(message),trim(routine_name)
      
      message="Error message :"
      call format(message,cform_1, dim_cform_1)
      call format(file_name,cform_2, dim_cform_2)
      cform='('//trim(cform_1(1:dim_cform_1))//','//trim(cform_2(1:dim_cform_2))//')'
      WRITE(*,cform)message,file_name
      
   else if(severity == 1)then
      message= "WARINING, check the routine ->"
      call format(message,cform_1, dim_cform_1)
      call format(routine_name,cform_2, dim_cform_2)
      cform='('//trim(cform_1(1:dim_cform_1))//','//trim(cform_2(1:dim_cform_2))//')'
      WRITE(*,cform)message,routine_name
      message="Error message :"
      call format(message,cform_1, dim_cform_1)
      call format(file_name,cform_2, dim_cform_2)
      cform='('//trim(cform_1(1:dim_cform_1))//','//trim(cform_2(1:dim_cform_2))//')'
      WRITE(*,cform)message,file_name
   elseif(severity ==2)then
      err: if(status /= 0 .and.  status /= -1)then
         message="ERROR, status ="
         call format(message,cform_1, dim_cform_1)
         call format(routine_name,cform_2, dim_cform_2)
         cform='('//trim(cform_1(1:dim_cform_1))//','//'I4,A4,'//trim(cform_2(1:dim_cform_2))//')'
         WRITE(*,cform)message,status," in ",routine_name
         if(status == 29)then
            message=" could not be found during an open operation."
            call format(message,cform_1, dim_cform_1)
            call format(file_name,cform_2, dim_cform_2)
            cform='(A4'//trim(cform_1(1:dim_cform_1))//','//trim(cform_2(1:dim_cform_2))//')'
            WRITE(*,'(A5,A200,A40)')"The", trim(file_name),message
            call tracebackqq()
         elseif(status == 30)then
            message=" Open failure file ->"
            call format(message,cform_1, dim_cform_1)
            call format(file_name,cform_2, dim_cform_2)
            cform='('//trim(cform_1(1:dim_cform_1))//trim(cform_2(1:dim_cform_2))//')'
            WRITE(*,cform)message,trim(file_name)
            call tracebackqq()
         elseif(status == 38)then
            message="REWIND error ->"
            call format(message,cform_1, dim_cform_1)
            call format(file_name,cform_2, dim_cform_2)
            cform='('//trim(cform_1(1:dim_cform_1))//','//trim(cform_2(1:dim_cform_2))//')'
            WRITE(*,cform)message,trim(file_name)
            call tracebackqq()
         elseif(status == 39)then
            message= "ENDFILE error ->"
            call format(message,cform_1, dim_cform_1)
            call format(file_name,cform_2, dim_cform_2)
            cform='('//trim(cform_1(1:dim_cform_1))//','//trim(cform_2(1:dim_cform_2))//')'
            WRITE(*,cform)message,trim(file_name)
            call tracebackqq()
            !endif
         elseif(status == 67)then
            message= "Input statement requires too much data ->"
            call format(message,cform_1, dim_cform_1)
            call format(file_name,cform_2, dim_cform_2)
            cform='('//trim(cform_1(1:dim_cform_1))//','//trim(cform_2(1:dim_cform_2))//')'
            WRITE(*,cform)message,trim(file_name)
            call tracebackqq()
         elseif(status == 151)then
            message="Allocatable array or pointer is already allocated"
            call format(message,cform_1, dim_cform_1)
            call format(file_name,cform_2, dim_cform_2)
            !cform ='('//cform_1(1:dim_cform_1)//')'
            cform='('//trim(cform_1(1:dim_cform_1))//','//trim(cform_2(1:dim_cform_2))//')'
            WRITE(*,cform)message,trim(file_name)
            call tracebackqq()
         elseif(status == 153)then
            message="Allocatable array or pointer is not allocated"
            call format(message,cform_1, dim_cform_1)
            call format(file_name,cform_2, dim_cform_2)
            !cform ='('//cform_1(1:dim_cform_1)//')'
            cform='('//trim(cform_1(1:dim_cform_1))//','//trim(cform_2(1:dim_cform_2))//')'
            WRITE(*,cform)message,trim(file_name)
            call tracebackqq()
         elseif(status >=985 .and. status <=1000)then
            message="ALLOCATION- DEALLOCATION ERROR"
            call format(message,cform_1, dim_cform_1)
            call format(file_name,cform_2, dim_cform_2)
            !cform ='('//cform_1(1:dim_cform_1)//')'
            cform='('//trim(cform_1(1:dim_cform_1))//','//trim(cform_2(1:dim_cform_2))//')'
            WRITE(*,cform)message,trim(file_name)
            call tracebackqq()
         else
            message= "UNKNOWN error ->"
            call format(message,cform_1, dim_cform_1)
            call format(file_name,cform_2, dim_cform_2)
            cform='('//trim(cform_1(1:dim_cform_1))//','//trim(cform_2(1:dim_cform_2))//')'
            WRITE(*,cform)message,trim(file_name)
            call tracebackqq()
         endif
      elseif(status == -1)then
         !WRITE(*,*)"END of FILE  -->",file_name
      endif err
      
   endif sever
end if
return
end subroutine error

subroutine format(string,cform,dim_cform)
  IMPLICIT NONE
  CHARACTER (len =*)::string,cform
  CHARACTER (len =200)::aus,error_mex, routine_name
  INTEGER :: dim,dim_cform
  
  routine_name="format"
  dim =len_trim(string)
  if(dim >= 0)then
     if(dim <= 9) then
        write(aus,'(I1)')dim
        dim_cform=1+1
     elseif(dim>=10 .and. dim <= 99)then
        write(aus,'(I2)')dim
        dim_cform=1+2
     else
        write(aus,'(I3)')dim
        dim_cform=1+3
     endif
     cform(1:dim_cform) = 'A'//trim(aus)
  else
     error_mex="string dimesnion must be greater than 0"
     call error(routine_name,-999,2,error_mex)
  endif
end subroutine format


END MODULE ERROR_MANAGEMENT
