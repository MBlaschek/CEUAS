program rcmstart

use rfmod, only: rasocor_namelist,JPRM
use  rfcor 

implicit none

type(rasocor_namelist) rcpara
character*80 filename

  call getarg(1,filename)
  call rcpara_ini(rcpara) !in file rfmod.f90
  call rasocorrect_main(rcpara) !in file rasocorrect_main.f90


  stop

end program rcmstart


