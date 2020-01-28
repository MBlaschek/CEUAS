subroutine dput(ncid,datumvarid,datum,status,ni,nj)

use netcdf
integer ncid,datumvarid,status,ni,nj
integer datum(ni,nj)
print*,ncid,datumvarid,shape(datum)
      status=NF90_PUT_VAR(ncid,datumvarid,datum(:,1))
      IF (status /= NF90_NOERR) PRINT *,'could not write datum', NF90_STRERROR(status)
return
end subroutine dput
