subroutine abor1(message)

implicit none

character*(*) message

  write(*,*) message
  call abort

end subroutine abor1
