module rfmod

use parkind1

type,public:: rasocor_namelist
  integer(kind=4) :: statmax,pmax,parmax,nmax,mmax,first,last,snht_maxlen,snht_increment,max_miss,rimax_miss,nmem,ngroup,nmeta,probmax
  integer(kind=4) :: bg_final(2),bg_initial(2),use_meta,bghom,minst
  integer(kind=4) :: anfsonde,endsonde,maxiter,smooth_method,brmax,select_mode,prob_method,bgbrmax,ni,nj,mean_maxlen,old,cachemax
  real(kind=JPRM) :: miss_val,break_thresh,break_thresh_rad,break_thresh_prob,locsig,break_fak,ap_prob_default,bg_correction_factor
  real(kind=JPRM) :: smooth_fit_thresh,sig_thresh,sig_thresh_rad
  real(kind=JPRM),allocatable :: mweights(:),tsaweights(:),smoothweights(:),plevs(:)
  real(kind=JPRM) :: qcthresh,weight_distance
  integer(kind=4),allocatable :: year(:),month(:),day(:),time(:)
  character*60 prefix
!  character*8  startdate
  character*2  innov,ens
  character*1  extended_output
  character*1  qc
  character*6  version
  logical downwards
  integer(kind=4)      switchdate,startdate
  character*20 fgdepname
  character*9  initial_adjust
end type rasocor_namelist

type,public:: cacherecord
  integer ::                     vals
  integer(kind=4),allocatable :: index(:)
  real(kind=JPRM),allocatable :: feld(:,:,:)
end type cacherecord

type,public:: metadata
  integer(kind=JPRM),allocatable :: cardsmeta_s(:,:)
  integer(kind=JPRM),allocatable :: rsicodes(:)
  logical(kind=JPRM),allocatable :: trusted(:)
  character*3,allocatable ::  rscodes(:)
end type metadata

integer*8, public :: omp_lp(8000)

contains

subroutine rcpara_ini(rcpara)

implicit none


character*80 parafile

type(rasocor_namelist) rcpara

integer statmax,pmax,parmax,nmax,mmax,first,last,snht_maxlen,&
snht_increment,max_miss,rimax_miss,mean_maxlen,nmem,ngroup,nmeta,probmax,&
prob_method,anfsonde,endsonde,maxiter,smooth_method,brmax,old,bghom,i
integer anfmon,anfyear,anfday,tdiff,index,year1,month1,day1,time1,ni,nj,cachemax,bg_initial(2),bg_final(2),use_meta,minst,switchdate
logical downwards

real(kind=JPRM) :: miss_val, &
ap_prob_default,break_thresh,break_thresh_prob,break_thresh_rad,break_fak,&
locsig,bg_correction_factor,smooth_fit_thresh,qcthresh,weight_distance,sig_thresh,sig_thresh_rad

  character*60 prefix
!  character*8  startdate
  integer       startdate
  character*2  innov,ens
  character*1  extended_output
  character*1  qc
  character*9  version,initial_adjust
  character*20 fgdepname

  integer(kind=4),allocatable :: year(:),month(:),day(:),time(:)

namelist /rfpar/statmax,pmax,parmax,nmax,mmax,first,last,startdate,prefix,snht_maxlen,snht_increment,miss_val,max_miss,rimax_miss,mean_maxlen,nmem,ngroup,nmeta,probmax,anfsonde,endsonde,maxiter,smooth_method,prob_method,brmax, &
break_thresh,break_thresh_prob,break_thresh_rad,break_fak,locsig,ap_prob_default,innov,extended_output,bg_correction_factor,smooth_fit_thresh,qc,&
old,weight_distance,version,bg_initial,bg_final,use_meta,sig_thresh,sig_thresh_rad,bghom,minst,ens,downwards,switchdate,fgdepname,initial_adjust

statmax=-1
pmax=-1
parmax=-1
nmax=-1
mmax=-1
first=-1
last=-1
snht_increment=-1
snht_maxlen=-1
mean_maxlen=-1
max_miss=-1
rimax_miss=-1
nmem=-1
ngroup=-1
nmeta=-1
probmax=-1
anfsonde=0
endsonde=99999
maxiter=-1
smooth_method=-1
prob_method=-1
innov='N'
ap_prob_default=-1
break_thresh=-1
break_thresh_prob=-1
break_thresh_rad=-1
break_fak=-1
locsig=-1
brmax=-1
extended_output='N'
bg_correction_factor=-1
smooth_fit_thresh=-1
old=-1
qc='X'
qcthresh=10.
weight_distance=-1.
cachemax=3000
version=''
bg_initial=0
bg_final=0
use_meta=-1
sig_thresh=-1.0
sig_thresh_rad=-1.0
bghom=-1
minst=-1
ens='-1'
downwards=.false.
switchdate=0
fgdepname='fg_dep'
initial_adjust='missin'

call getarg(1,parafile)
open(11,file=trim(parafile),form='formatted',status='old',action='read')
read(11,rfpar)
close(11)

allocate(rcpara%year(nmax),rcpara%month(nmax),rcpara%day(nmax),rcpara%time(nmax))
allocate(year(nmax),month(nmax),day(nmax),time(nmax))

allocate(rcpara%mweights(pmax),rcpara%tsaweights(pmax),rcpara%smoothweights(pmax),rcpara%plevs(pmax))

rcpara%plevs=(/10.,20.,30.,50.,70.,100.,150.,200.,250.,300.,400.,500.,700.,850.,925.,1000./)
rcpara%mweights=(/1.,1.,1.,1.,1.,1.,1.,1.,0.,1.,1.,1.,1.,1.0,0.,0./) ! 10,20,250,925,1000 nicht
rcpara%smoothweights=(/1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0./) ! 250,925,1000 nicht, 850 almost fixed (to zero) 

rcpara%statmax=statmax
rcpara%pmax=pmax
rcpara%parmax=parmax
rcpara%nmax=nmax
rcpara%mmax=mmax
rcpara%cachemax=cachemax
rcpara%first=first
rcpara%last=last
rcpara%snht_maxlen=snht_maxlen
rcpara%snht_increment=snht_increment
rcpara%miss_val=miss_val
rcpara%max_miss=max_miss
rcpara%rimax_miss=rimax_miss
rcpara%mean_maxlen=mean_maxlen
rcpara%nmem=nmem
rcpara%ngroup=ngroup
rcpara%nmeta=nmeta 
rcpara%probmax=probmax
rcpara%prob_method=prob_method
rcpara%anfsonde=anfsonde
rcpara%endsonde=endsonde
rcpara%maxiter=maxiter
rcpara%smooth_method=smooth_method
rcpara%brmax=brmax
rcpara%prefix=prefix
rcpara%startdate=startdate
rcpara%ap_prob_default=ap_prob_default
rcpara%break_thresh=break_thresh
rcpara%break_thresh_prob=break_thresh_prob
rcpara%break_thresh_rad=break_thresh_rad
rcpara%break_fak=break_fak
rcpara%locsig=locsig
rcpara%innov=innov
rcpara%extended_output=extended_output
rcpara%bg_correction_factor=bg_correction_factor
rcpara%smooth_fit_thresh=smooth_fit_thresh
rcpara%qc=qc
rcpara%old=old
rcpara%qcthresh=qcthresh
rcpara%weight_distance=weight_distance
rcpara%version=version
rcpara%use_meta=use_meta
rcpara%sig_thresh=sig_thresh
rcpara%sig_thresh_rad=sig_thresh_rad
rcpara%bghom=bghom
rcpara%minst=minst
rcpara%ens=ens
rcpara%downwards=downwards
rcpara%switchdate=switchdate
if(bg_final(1) .gt. rcpara%switchdate-1) then
  bg_final(1)=rcpara%switchdate-1
endif
rcpara%fgdepname=fgdepname
rcpara%initial_adjust=initial_adjust
!print *, rcpara

rcpara%ni=360
rcpara%nj=181
ni=360
nj=181

if(statmax .eq. -1 .or. pmax .eq. -1 .or. parmax .eq. -1 .or. nmax .eq. -1 .or. mmax .eq. -1 .or. first .eq. -1 .or. last .eq. -1 &
.or. snht_increment .eq. -1 .or. snht_maxlen .eq. -1 .or. mean_maxlen .eq. -1 .or. max_miss .eq. -1 .or. rimax_miss .eq. -1 .or. nmem .eq. -1 .or. ngroup .eq. -1 .or. nmeta .eq. -1 .or. &
probmax .eq. -1 .or. ap_prob_default .eq. -1 .or. break_thresh .eq. -1 .or. break_thresh_prob .eq. -1 .or. break_thresh_rad .eq. -1 &
.or. locsig .eq. -1 .or. brmax .eq. -1 .or. innov .eq. 'N' .or. smooth_fit_thresh .eq. -1 .or. bg_correction_factor .eq. -1 .or. qc .eq. 'X' .or. old .eq. -1 .or. weight_distance .eq. -1. &
.or. version .eq. '' .or. bg_initial(1) .eq. 0 .or. bg_final(1) .eq. 0 .or. use_meta .eq. -1 .or. sig_thresh .eq. -1.0 .or. sig_thresh_rad .eq. -1.0 .or. bghom .eq. -1 .or. minst .eq. -1 .or. ens .eq. '-1' .or. switchdate .eq. 0 .or. initial_adjust .eq. 'missin') then
  write(*,*) 'At least one essential parameter is not set in namelist'
  write(*,*) statmax,pmax,parmax,nmax,mmax,first,last,snht_increment,snht_maxlen,max_miss,rimax_miss,nmem,ngroup,nmeta,probmax, &
ap_prob_default,break_thresh,break_thresh_prob,break_thresh_rad,locsig,brmax,innov,old,'version ',version,bg_initial,bg_final,use_meta,sig_thresh,sig_thresh_rad,bghom,minst
  call exit(1)
else
  write(*,*) 'Namelist parameters:'
  write(*,*) statmax,pmax,parmax,nmax,mmax,first,last,snht_increment,snht_maxlen,max_miss,nmem,ngroup,nmeta,probmax, ap_prob_default,break_thresh,&
break_thresh_prob,break_thresh_rad,locsig,brmax,innov,old,'version ',version,bg_initial,bg_final,use_meta,sig_thresh,sig_thresh_rad,bghom,minst,ens,switchdate,fgdepname,initial_adjust

endif

!read(startdate,'(I4,2I2.2)') anfyear,anfmon,anfday
anfyear=startdate/10000
anfmon=(startdate-anfyear*10000)/100
anfday=mod(startdate,100)
year(1)=anfyear
month(1)=anfmon
day(1)=anfday
time(1)=00
tdiff=24
year1=anfyear
index=1
do while(index .lt. nmax)
  CALL DATES(Year(index),Month(index),Day(index),Time(index),Tdiff, &
             Year1,Month1,Day1,Time1)
  index=index+1
  year(index)=year1
  month(index)=month1
  day(index)=day1
  time(index)=time1
!  write(*,*) year1,month1,day1,time1,index
  if(year1 .eq. 2030) stop 'nmax or startdate are invalid'
enddo
index=index-1
rcpara%year=year
rcpara%month=month
rcpara%day=day
rcpara%time=time

do i=1,2
  if(bg_initial(i) .gt. 0 .and. bg_final(i) .gt. 0) then
    rcpara%bg_initial(i)=toindex(bg_initial(i),rcpara)
    rcpara%bg_final(i)=toindex(bg_final(i),rcpara)
  else
    rcpara%bg_initial(i)=1
    rcpara%bg_final(i)=1
  endif

enddo


return
end subroutine rcpara_ini


subroutine rcparapmaxparmax(rcpara,pmax,parmax) !

implicit none

type(rasocor_namelist) rcpara

integer pmax,parmax

rcpara%pmax=pmax
rcpara%parmax=parmax

return
end subroutine rcparapmaxparmax

subroutine ulfcorr(rcpara,istat,wmonrs,wmolons,wmolats,wmostats,newbc)

implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer,intent(in) :: wmostats,istat
integer,intent(in) :: wmonrs(rcpara%statmax)
real(kind=JPRM),intent(in) :: wmolats(rcpara%statmax),wmolons(rcpara%statmax)

integer i,ip,ipar,idat,idelt
real(kind=JPRM) :: newbc(rcpara%nmax,rcpara%pmax,rcpara%parmax)

LOGICAL              :: LD_LBC = .FALSE.
REAL(KIND=JPRM)      :: P_SPP(rcpara%pmax),P_ST(rcpara%pmax)		! values
REAL(KIND=JPRM)      :: P_STNBC(rcpara%pmax)	! Corrected values
CHARACTER(LEN=5)     :: CD_CIDENT		! Station ID
CHARACTER(LEN=5)     :: CL_solar ='ON'	 	! solar correction switch
CHARACTER(LEN=5)     :: CL_homogen ='OFF'       ! homogenization switch
INTEGER(KIND=JPIM)   :: K_IY = 1986
INTEGER(KIND=JPIM)   :: K_IM = 06
INTEGER(KIND=JPIM)   :: K_ID = 01
INTEGER(KIND=JPIM)   :: K_IH = 12
INTEGER(KIND=JPIM)   :: K_IRSTYP = 0
INTEGER(KIND=JPIM)   :: K_IMISS = 0
REAL(KIND=JPRM)      :: P_RLAT
REAL(KIND=JPRM)      :: P_RLON
REAL(KIND=JPRM)      :: P_RMISS = 99999.0

write(CD_CIDENT,'(I6.6)') wmonrs(istat)
do ip=1,rcpara%pmax
  P_SPP(rcpara%pmax-ip+1)=rcpara%plevs(ip)*100.
enddo
P_ST=273.
newbc=rcpara%miss_val
idelt=20
do i=366+idelt/2,rcpara%nmax,idelt
  idat=todate(i,rcpara)
  K_IY=idat/10000
  K_IM=(idat-K_IY*10000)/100
  K_ID=mod(idat,100)
  do ipar=1,rcpara%parmax
    K_IH=(ipar-1)*12
    if(K_IY .lt. 2005) CALL BIASCOR_ERA40 (LD_LBC,rcpara%pmax,P_SPP,P_ST,P_STNBC, &
    CD_CIDENT,K_IY,K_IM,K_ID,K_IH,wmolats(istat),wmolons(istat),K_IRSTYP,K_IMISS,P_RMISS,CL_solar,CL_homogen)

   do ip=1,rcpara%pmax
     newbc(i-idelt/2:i+idelt/2-1,ip,ipar)=P_STNBC(rcpara%pmax-ip+1)-P_ST(rcpara%pmax-ip+1)
   enddo
  enddo
enddo

end subroutine ulfcorr


integer function todate(index,rcpara)
implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer index

if(index .gt. 0 .and. index .le. rcpara%nmax) then 
  todate=rcpara%year(index)*10000+rcpara%month(index)*100+rcpara%day(index)
else
  todate=0
endif

return
end function todate

integer function toindex(date,rcpara,abo)
implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer index,date,istart,istop,idate
logical found,ab
logical,optional :: abo

if(.not. present(abo)) then
  ab=.true.
else
  ab=abo
endif

istart=1
istop=rcpara%nmax
index=1
found=.false.
if(date .gt. rcpara%startdate .and. date .lt. 20281231) then 
do while(istop .gt. istart+1)
  index=(istart+istop)/2
  idate=rcpara%year(index)*10000+rcpara%month(index)*100+rcpara%day(index)
  if(date .eq. idate) then
    found=.true.
    exit
  endif
  if(date .gt. idate) then
    istart=index
  else
    istop=index
  endif
enddo
  if(found) then
    toindex=index
  else
    write(*,*) 'invalid date ',date,index,istart,istop
    if(ab) then
      call abort
    endif
    toindex=-1
  endif
else
  if(date .eq. rcpara%startdate) then
    toindex=1
  else
    write(*,*) 'invalid date ',date
    if(ab) then
      call abort
    endif
    toindex=-1
  endif
endif

return
end function toindex

integer function toindexold(date,rcpara)
implicit none

type(rasocor_namelist),intent(in) :: rcpara

integer index,date

index=1
if(date .gt. rcpara%startdate .and. date .lt. 20281231) then 
  do while (rcpara%year(index)*10000 .lt. date-20000)
    index=index+365
  enddo
  do while (rcpara%year(index)*10000+rcpara%month(index)*100+rcpara%day(index) .lt. date) 
    index=index+1
  enddo
  toindexold=index
else
  toindexold=0
endif

return
end function toindexold

integer function bsearchi(xin,gx)

implicit none

integer :: i,istart,istop,xin,gx(:)

istop=size(gx)    
        istart=1
        do while (istart .lt. istop)
           i=(istart+istop)/2
           if(xin .eq. gx(i)) then  !found what we want
              bsearchi=i
              return
           endif
           if(istart+1.eq. istop) exit
           if(xin .gt. gx(i)) then !xin is in second half of array
              istart=i
           else !xin is in first half of array
              istop=i
           endif
        enddo
        bsearchi=-1
        return

end function bsearchi

subroutine monanfg(mon,tan,tanfg,tanfgall,wmonrs,wmolons,wmolats, &
               mmax,pmax,wmostats,tmax,NGI,NGJ,prefix)

implicit none

integer pmax,tmax,mmax,wmostats
integer mon,i,iunit,it,is,NGI,NGJ,loni,lati,loniplus1,ip

integer wmonrs(:)
real(kind=JPRM) :: wmolats(:),wmolons(:)

logical ex

real(kind=JPRM) :: tan(NGI,NGJ,pmax,tmax)
real(kind=JPRM) :: tanfg(pmax,wmostats,tmax)
real(kind=JPRM) :: tanfgall(wmostats,pmax,mmax,tmax)
real(kind=JPRM) :: w1(pmax),w2(pmax),w3(pmax)
real(kind=JPRM) :: lonfrac,latfrac

character*(*) prefix

  print *,prefix
!  print*,wmonrs

  iunit=20
  inquire(file=prefix,exist=ex)
  if(ex) then
  open(iunit,file=prefix,form='unformatted')
  read(iunit) wmostats,pmax,mmax,tmax
  read(iunit) tanfgall
  close(iunit)
  else
    tanfgall=-999.
  endif

do it=1,tmax
  do is=1,wmostats

    latfrac=wmolats(is)-floor(wmolats(is))
    lonfrac=wmolons(is)-floor(wmolons(is))

    lati=floor(wmolats(is))+91
    if(wmolons(is) .lt. 0.) then
      loni=361+floor(wmolons(is))
    else
      loni=1+floor(wmolons(is))
    endif
    loniplus1=loni+1
    if(loniplus1 .eq. 361) loniplus1=1

! simple bilinear interpolation to station location

    if(latfrac .eq. 0 .and. lonfrac .eq. 0) then
      tanfg(:,is,it)=tan(loni,lati,:,it)
    else

      w1=tan(loni,lati,:,it)*(1.-lonfrac)+tan(loniplus1,lati,:,it)*lonfrac
      w2=tan(loni,lati+1,:,it)*(1.-lonfrac)+tan(loniplus1,lati+1,:,it)*lonfrac
      w3=w1*(1.-latfrac)+w2*latfrac

      tanfg(:,is,it)=w3

    endif

  enddo
enddo

do it=1,tmax
  do ip=1,pmax
    tanfgall(:,ip,mon,it)=tanfg(ip,:,it)
  enddo
enddo

  open(iunit,file=prefix,form='unformatted')
  write(iunit) wmostats,pmax,mmax,tmax
  write(iunit) tanfgall
  close(iunit)

return
end subroutine monanfg


subroutine momean(mon,values,times,statnrs,wmonrs,ts,tbcs,newtbcs,newradbcs,solarangles,tanfeeds,tfgfeeds,tflags, &
               goodvals,tm,tbcm,newtbcm,newradbcm,solaranglem,tanfeedm,tfgfeedm,tsig,tanfeedsig,tfgfeedsig, &
               statmax,pmax,parmax,nmax,mmax,flagmax,tmax,wmostats,prefix)

implicit none

integer statmax,pmax,parmax,nmax,mmax,tmax,first,last,flagmax,wmostats
integer mon,i,iunit,l,height,time,statnr,it,tim100,ip,ik,is,nv,iret

integer values(pmax,statmax,tmax)
integer times(nmax,pmax,statmax,tmax)
integer statnrs(statmax,tmax)

integer wmonrs(statmax)

logical ex

real(kind=JPRM) ts(nmax,pmax,statmax,tmax)
real(kind=JPRM) tbcs(nmax,pmax,statmax,tmax)
real(kind=JPRM) newtbcs(nmax,pmax,statmax,tmax)
real(kind=JPRM) newradbcs(nmax,pmax,statmax,tmax)
real(kind=JPRM) solarangles(nmax,statmax,tmax)
real(kind=JPRM) tanfeeds(nmax,pmax,statmax,tmax)
real(kind=JPRM) tfgfeeds(nmax,pmax,statmax,tmax)
integer*1 tflags !(flagmax,nmax,pmax,statmax,tmax)

integer goodvals(pmax,statmax,tmax)
real(kind=JPRM) tm(pmax,statmax,tmax)
real(kind=JPRM) tbcm(pmax,statmax,tmax)
real(kind=JPRM) newtbcm(pmax,statmax,tmax)
real(kind=JPRM) newradbcm(pmax,statmax,tmax)
real(kind=JPRM) solaranglem(statmax,tmax)
real(kind=JPRM) tanfeedm(pmax,statmax,tmax)
real(kind=JPRM) tfgfeedm(pmax,statmax,tmax)

real(kind=JPRM) tsig(pmax,statmax,tmax)
real(kind=JPRM) tanfeedsig(pmax,statmax,tmax)
real(kind=JPRM) tfgfeedsig(pmax,statmax,tmax)

real(kind=JPRM) tzwisch(nmax),mean

character*8 cdatum
character*2 czeit
character*(*) prefix

!write(cdatum,'(I4,I2.2,I2.2)') 1957+mon/12,mod(mon,12)+1,0

do it=1,tmax
  do is=1,wmostats
    do ip=1,pmax
      nv=values(ip,is,it)
      if(nv .gt. 1) then 
        tzwisch(1:nv)=ts(1:nv,ip,is,it)
        tm(ip,is,it)=sum(tzwisch(1:nv))/nv
        mean=tm(ip,is,it)
        tsig(ip,is,it)=sqrt(dot_product(tzwisch(1:nv),tzwisch(1:nv))/(nv-1)-mean*mean)
      
        tzwisch(1:nv)=tbcs(1:nv,ip,is,it)
        tbcm(ip,is,it)=sum(tzwisch(1:nv))/nv
        mean=tbcm(ip,is,it)
!        tbcsig(ip,is,it)=sqrt(dot_product(tzwisch(1:nv),tzwisch(1:nv))/(nv-1)-mean*mean)

        tzwisch(1:nv)=newtbcs(1:nv,ip,is,it)
        newtbcm(ip,is,it)=sum(tzwisch(1:nv))/nv
        mean=newtbcm(ip,is,it)
!        tbcsig(ip,is,it)=sqrt(dot_product(tzwisch(1:nv),tzwisch(1:nv))/(nv-1)-mean*mean)

        tzwisch(1:nv)=newradbcs(1:nv,ip,is,it)
        newradbcm(ip,is,it)=sum(tzwisch(1:nv))/nv
        mean=newradbcm(ip,is,it)
!        tbcsig(ip,is,it)=sqrt(dot_product(tzwisch(1:nv),tzwisch(1:nv))/(nv-1)-mean*mean)

        if(ip .eq. 1) then 
          solaranglem(is,it)=sum(solarangles(1:nv,is,it))/nv
        endif

        tzwisch(1:nv)=tanfeeds(1:nv,ip,is,it)
        tanfeedm(ip,is,it)=sum(tzwisch(1:nv))/nv
        mean=tanfeedm(ip,is,it)
        tanfeedsig(ip,is,it)=sqrt(dot_product(tzwisch(1:nv),tzwisch(1:nv))/(nv-1)-mean*mean)

        tzwisch(1:nv)=tfgfeeds(1:nv,ip,is,it)
        tfgfeedm(ip,is,it)=sum(tzwisch(1:nv))/nv
        mean=tfgfeedm(ip,is,it)
        tfgfeedsig(ip,is,it)=sqrt(dot_product(tzwisch(1:nv),tzwisch(1:nv))/(nv-1)-mean*mean)
  
      endif 
      if(nv.eq.1) then
        tm(ip,is,it)=ts(1,ip,is,it)
        tsig(ip,is,it)=-999.
        tbcm(ip,is,it)=tbcs(1,ip,is,it)
        newtbcm(ip,is,it)=newtbcs(1,ip,is,it)
        newradbcm(ip,is,it)=newradbcs(1,ip,is,it)
        solaranglem(is,it)=solarangles(1,is,it)
        tanfeedm(ip,is,it)=tanfeeds(1,ip,is,it)
        tanfeedsig(ip,is,it)=-999.
        tfgfeedm(ip,is,it)=tfgfeeds(1,ip,is,it)
        tfgfeedsig(ip,is,it)=-999.
      endif     
    enddo
  enddo
enddo
  print *,prefix
!  print*,wmonrs
  iunit=20
  open(iunit,file=prefix,form='unformatted')
  write(iunit) pmax,wmostats,tmax
  write(iunit) values(:,1:wmostats,:),tm(:,1:wmostats,:),tbcm(:,1:wmostats,:),newtbcm(:,1:wmostats,:),newradbcm(:,1:wmostats,:),tanfeedm(:,1:wmostats,:),tfgfeedm(:,1:wmostats,:)
  write(iunit) solaranglem(1:wmostats,:)
  write(iunit) tsig(:,1:wmostats,:),tanfeedsig(:,1:wmostats,:),tfgfeedsig(:,1:wmostats,:)
  close(iunit)

return
end subroutine momean

subroutine momeancards(first,last,values,times,statnr,ts,tbcs,tflags, &
               tm,tbcm,tsig,tbcsig,corts,sondetypes, &
               pmax,nmax,mmax,flagmax,tmax)

implicit none

integer pmax,nmax,mmax,tmax,first,last,flagmax,statnr
integer mon,i,iunit,l,height,time,it,tim100,ip,ik,nv,iret,im

integer values(pmax,mmax,tmax)
integer times(nmax,pmax,mmax,tmax)
integer corts(mmax),sondetypes(mmax)

logical ex

real(kind=JPRM) ts(nmax,pmax,mmax,tmax)
real(kind=JPRM) tbcs(nmax,pmax,mmax,tmax)
integer*1 tflags !(flagmax,nmax,pmax,mmax,tmax)

real(kind=JPRM) tm(pmax,mmax,tmax)
real(kind=JPRM) tbcm(pmax,mmax,tmax)

real(kind=JPRM) tsig(pmax,mmax,tmax)
real(kind=JPRM) tbcsig(pmax,mmax,tmax)

real(kind=JPRM) tzwisch(nmax),mean

character*8 cdatum
character*2 czeit
character*6 cstatnr

!!$ call omp_set_lock(omp_lp)
write(cstatnr,'(I6.6)') statnr
!!$ call omp_unset_lock(omp_lp)
do it=1,tmax
  do im=first,last
    do ip=1,pmax
      nv=values(ip,im,it)
      if(nv .gt. 1) then 
        tzwisch(1:nv)=ts(1:nv,ip,im,it)
        tm(ip,im,it)=sum(tzwisch(1:nv))/nv
        mean=tm(ip,im,it)
        tsig(ip,im,it)=sqrt(dot_product(tzwisch(1:nv),tzwisch(1:nv))/(nv-1)-mean*mean)
      
        tzwisch(1:nv)=tbcs(1:nv,ip,im,it)
        tbcm(ip,im,it)=sum(tzwisch(1:nv))/nv
        mean=tbcm(ip,im,it)
        tbcsig(ip,im,it)=sqrt(dot_product(tzwisch(1:nv),tzwisch(1:nv))/(nv-1)-mean*mean)
 
      else
!       if(values(ip,im,it).eq.1) print*,'values(',ip,im,it,') set to zero'
       values(ip,im,it)=0
       tm(ip,im,it)=-999. 
       tsig(ip,im,it)=-999. 
       tbcm(ip,im,it)=-999. 
       tbcsig(ip,im,it)=-999. 
      endif     
    enddo
  enddo
!  write(*,'(540F8.2,540I8)') tm(12,:,1),values(12,:,1)
  write(czeit,'(I2.2)') it
  print *,cdatum
!  print*,wmonrs
enddo
  iunit=20
  open(iunit,file='cards'//cstatnr//'bin',form='unformatted')
  write(iunit) pmax,mmax,first,last,statnr
  write(iunit) values,tm,tbcm
!  write(iunit) tsig,tbcsig
  write(iunit) corts,sondetypes
  close(iunit)

return
end subroutine momeancards

subroutine readmomean(first,last,it,values,tm,tbcm,tanfeedm,tfgfeedm,tsig,tbcsig,tanfeedsig,tfgfeedsig, &
               statmax,pmax,parmax,nmax,mmax,flagmax,tmax,wmostats)

implicit none

integer statmax,pmax,parmax,nmax,mmax,tmax,first,last,flagmax,wmostats
integer mon,i,iunit,l,height,time,statnr,it,tim100,ip,ik,is,nv

integer values(pmax,statmax,tmax)
integer statnrs(statmax,tmax)

integer wmonrs(statmax)

real(kind=JPRM) tm(pmax,statmax,tmax)
real(kind=JPRM) tbcm(pmax,statmax,tmax)
real(kind=JPRM) tanfeedm(pmax,statmax,tmax)
real(kind=JPRM) tfgfeedm(pmax,statmax,tmax)

real(kind=JPRM) tsig(pmax,statmax,tmax)
real(kind=JPRM) tbcsig(pmax,statmax,tmax)
real(kind=JPRM) tanfeedsig(pmax,statmax,tmax)
real(kind=JPRM) tfgfeedsig(pmax,statmax,tmax)

character*8 cdatum
character*2 czeit


do mon=first,last
  write(cdatum,'(I4,I2.2,I2.2)') 1957+mon/12,mod(mon,12)+1,0
  write(czeit,'(I2.2)') it
  print *,cdatum

  open(iunit,file='feedbackglobbin'//cdatum//czeit,form='unformatted',status='old',err=20)
  read(iunit,err=20,end=20) pmax,wmostats
  read(iunit,err=20,end=20) values(:,:,mon),tm(:,:,mon),tbcm(:,:,mon),tanfeedm(:,:,mon),tfgfeedm(:,:,mon)
  read(iunit,err=20,end=19) tsig(:,:,mon),tanfeedsig(:,:,mon),tfgfeedsig(:,:,mon)
  close(iunit)
  19 goto 21
  20 print*, 'could not open/read feedbackglobbin'//cdatum//czeit
  21 continue
enddo

return

end subroutine readmomean

subroutine readmomeanalltimes(first,last,it,values,tm,tbcm,tanfeedm,tfgfeedm,values12,tm12,tbcm12,tanfeedm12,tfgfeedm12, &
               statmax,pmax,parmax,nmax,mmax,flagmax,tmax,wmostats)

implicit none

integer statmax,pmax,parmax,nmax,mmax,tmax,first,last,flagmax,wmostats
integer mon,i,iunit,l,height,time,statnr,it,tim100,ip,ik,is,nv,ios

integer values(pmax,statmax,tmax),values12(pmax,statmax,tmax)
integer statnrs(statmax,tmax)

integer wmonrs(statmax)
logical ex

real(kind=JPRM) tm(pmax,statmax,tmax)
real(kind=JPRM) tbcm(pmax,statmax,tmax)
real(kind=JPRM) tanfeedm(pmax,statmax,tmax)
real(kind=JPRM) tfgfeedm(pmax,statmax,tmax)

real(kind=JPRM) tm12(pmax,statmax,tmax)
real(kind=JPRM) tbcm12(pmax,statmax,tmax)
real(kind=JPRM) tanfeedm12(pmax,statmax,tmax)
real(kind=JPRM) tfgfeedm12(pmax,statmax,tmax)

character*8 cdatum
character*2 czeit

iunit=20

do mon=first,last
  write(cdatum,'(I4,I2.2,I2.2)') 1957+mon/12,mod(mon,12)+1,0
  write(czeit,'(I2.2)') 1

  inquire(file='feedbackglobbin'//cdatum//czeit,exist=ex)
  print *,cdatum//czeit,ex
  open(iunit,file='feedbackglobbin'//cdatum//czeit,form='unformatted',status='old',iostat=ios,err=20)
  read(iunit,err=20,end=20) pmax,wmostats
  read(iunit,err=20,end=20) values(:,:,mon),tm(:,:,mon),tbcm(:,:,mon),tanfeedm(:,:,mon),tfgfeedm(:,:,mon)
  close(iunit)
  if(it.lt.2) goto 21
  write(czeit,'(I2.2)') 2
!  print *,cdatum//czeit
  open(iunit,file='feedbackglobbin'//cdatum//czeit,form='unformatted',status='old',err=20)
  read(iunit,err=20,end=20) pmax,wmostats
  read(iunit,err=20,end=20) values12(:,:,mon),tm12(:,:,mon),tbcm12(:,:,mon),tanfeedm12(:,:,mon),tfgfeedm12(:,:,mon)
  close(iunit)
  19 goto 21
  20 print*, 'could not open/read feedbackglobbin'//cdatum//czeit,' iostat=',ios
  21 continue
enddo

return

end subroutine readmomeanalltimes

subroutine  makecorr(psvals,values,times,statnrs,wmonrs,lats,lons,ts,newtbcs,newradbcs,solarangles, &
               statmax,pmax,parmax,nmax,mmax,flagmax,tmax,prefix,cdatum)

implicit none

INTEGER, PARAMETER :: JPRM= 4

integer statmax,pmax,parmax,nmax,mmax,tmax,first,last,flagmax
integer mon,i,iunit,l,height,time,statnr,it,ip,ik,iret,ios,iv,is
integer tim100

integer values(pmax,statmax,tmax)
integer times(nmax,pmax,statmax,tmax)
integer statnrs(statmax,tmax)

integer flags(flagmax),stattype

real(kind=JPRM) lats(nmax,statmax,tmax)
real(kind=JPRM) lons(nmax,statmax,tmax)
real(kind=JPRM) ts(nmax,pmax,statmax,tmax)
real(kind=JPRM) newtbcs(nmax,pmax,statmax,tmax)
real(kind=JPRM) newradbcs(nmax,pmax,statmax,tmax)
real(kind=JPRM) solarangles(nmax,statmax,tmax)

real(kind=JPRM) lat,lon,pres

integer wmonrs(statmax)
real(kind=JPRM) wmolats(statmax),wmolons(statmax)

real(kind=JPRM) ps,psvals(pmax)

character*(*) cdatum
character*(*) prefix
character*6 cstatnr

logical ex

INTEGER, PARAMETER :: JPIM = 4
INTEGER, PARAMETER :: MAXDAYS=31
INTEGER(KIND=JPIM) :: K_IM_HILF(MAXDAYS)
INTEGER(KIND=JPIM) :: K_ID_HILF(MAXDAYS)
INTEGER(KIND=JPIM) :: K_IH_HILF(MAXDAYS)
REAL(KIND=JPRM) :: THILF(pmax,MAXDAYS),LATHILF(MAXDAYS),LONHILF(MAXDAYS)
REAL(KIND=JPRM) :: NEWTBCHILF(pmax,MAXDAYS),NEWRADBCHILF(pmax,MAXDAYS)
REAL(KIND=JPRM) :: SOLARANGLEHILF(MAXDAYS),NEWTBCCOL(pmax)
REAL(KIND=JPRM) :: psvals_pa(pmax)
LOGICAL              :: LD_LBC = .FALSE.,LL_FILES_NOT_READ,DAYHILF(MAXDAYS)
INTEGER(KIND=JPIM)   :: K_IM
INTEGER(KIND=JPIM)   :: K_ID
INTEGER(KIND=JPIM)   :: K_IH
REAL(KIND=JPRM)      :: P_RMISS = -999.
INTEGER(KIND=JPIM)   :: I_MTBL,I_MCOR,I_MCORRAD,I_MSGT,INDEX,INDEXMAX

!write(cdatum,'(I4,I2.2)') 1957+mon/12,mod(mon,12)+1

psvals_pa=psvals*100.

LL_FILES_NOT_READ=.TRUE.
print *,prefix
!print*,wmonrs
I_MTBL=20
I_MCOR=21
I_MSGT=22
I_MCORRAD=23
i=0

!  OPEN(UNIT=I_MTBL, FILE='/era/work/erl/recalc_biases_era40/radiation_only/scr/country.t', &
  OPEN(UNIT=I_MTBL, FILE='country.t', &
   & IOSTAT=IOS,ERR=967,STATUS='OLD')  

!  OPEN(UNIT=I_MCOR, FILE='/era/work/erl/recalc_biases_era40/rad_and_mean/tab/corcand_'//cdatum(1:4)//'.t', &
  OPEN(UNIT=I_MCOR, FILE='rad_and_mean/corcand_'//cdatum(1:4)//'.t', &
   & IOSTAT=IOS,ERR=967,STATUS='OLD')  

!  OPEN(UNIT=I_MCORRAD, FILE='/era/work/erl/recalc_biases_era40/radiation_only/tab/corcand_'//cdatum(1:4)//'.t', &
  OPEN(UNIT=I_MCORRAD, FILE='radiation_only/corcand_'//cdatum(1:4)//'.t', &
   & IOSTAT=IOS,ERR=967,STATUS='OLD')  

!  OPEN(UNIT=I_MSGT, FILE='/era/work/erl/recalc_biases_era40/radiation_only/scr/stgroup.t', &
    OPEN(UNIT=I_MSGT, FILE='stgroup.t', &
   & IOSTAT=IOS,ERR=967,STATUS='OLD')  

DO IT=1,TMAX
  DO IS=1,STATMAX
    WRITE(CSTATNR,'(I6.6)') WMONRS(IS)
    THILF=-999.
    DAYHILF=.FALSE.
    INDEXMAX=0
    DO IP=1,PMAX
      DO IV=1,VALUES(IP,IS,IT)
        K_IM=mod(TIMES(IV,IP,IS,IT)/10000,100)
        K_ID=mod(TIMES(IV,IP,IS,IT)/100,100)
        K_IH=mod(TIMES(IV,IP,IS,IT),100)
        INDEX=K_ID
        DAYHILF(INDEX)=.TRUE.
        IF(K_ID.GT.INDEXMAX) INDEXMAX=K_ID
        THILF(IP,INDEX)=TS(IV,IP,IS,IT)
        LATHILF(INDEX)=LATS(IV,IS,IT)
        LONHILF(INDEX)=LONS(IV,IS,IT)
        K_IM_HILF(INDEX)=K_IM
        K_ID_HILF(INDEX)=K_ID
        K_IH_HILF(INDEX)=K_IH
      ENDDO
    ENDDO
    DO INDEX=1,INDEXMAX
     IF(DAYHILF(INDEX)) THEN
      CALL BIASCOR_ERA40 (LD_LBC,pmax,PSVALS_PA,THILF(:,INDEX),NEWTBCHILF(:,index),NEWRADBCHILF(:,index),SOLARANGLEHILF(index), &
   CSTATNR,K_IM_HILF(INDEX),K_ID_HILF(INDEX),K_IH_HILF(INDEX),LATHILF(INDEX),LONHILF(INDEX),&
   P_RMISS,I_MTBL,I_MCOR,I_MCORRAD,I_MSGT,LL_FILES_NOT_READ)
      IF (LD_LBC .and. .false. ) THEN
        WRITE(6,'('' RS TEMPERATURE BIAS CORRECTION DONE FOR:'')')
        WRITE(6,'(''    LAT, LON, STID: '',&
       & 2(F10.2,1X),A)') LATHILF(INDEX),LONHILF(INDEX),CSTATNR
        WRITE(6,'(16F7.1)') THILF(:,index)
        WRITE(6,'(16F7.1)') NEWTBCHILF(:,index)
        WRITE(6,'(16F7.1)') NEWRADBCHILF(:,index)
      ENDIF
     ENDIF
    ENDDO
    DO IP=1,PMAX
      DO IV=1,VALUES(IP,IS,IT)
        K_ID=mod(TIMES(IV,IP,IS,IT)/100,100)
        INDEX=K_ID
        NEWTBCS(IV,IP,IS,IT)=NEWTBCHILF(IP,INDEX)-THILF(IP,INDEX)
        NEWRADBCS(IV,IP,IS,IT)=NEWRADBCHILF(IP,INDEX)-THILF(IP,INDEX)
        SOLARANGLES(IV,IS,IT)=SOLARANGLEHILF(INDEX)
      ENDDO
    ENDDO
  ENDDO
ENDDO

return
967 continue
  write(*,*) 'bias correction table(s) could not be opened'
stop
return
end subroutine makecorr


subroutine makebin(mon,psvals,values,times,statnrs,wmonrs,wmolats,wmolons,lats,lons,ts,tbcs,tanfeeds,tfgfeeds,tflags, &
               statmax,pmax,parmax,nmax,mmax,flagmax,tmax,prefix,cdatum)

! externals: dates

implicit none

integer statmax,pmax,parmax,nmax,mmax,tmax,first,last,flagmax
integer mon,i,iunit,l,height,time,statnr,oldstatnr,it,ip,ik,iret
integer tim100

integer values(pmax,statmax,tmax)
integer times(nmax,pmax,statmax,tmax)
integer statnrs(statmax,tmax)

integer flags(flagmax),stattype

real(kind=JPRM) lats(nmax,statmax,tmax)
real(kind=JPRM) lons(nmax,statmax,tmax)
real(kind=JPRM) ts(nmax,pmax,statmax,tmax)
real(kind=JPRM) tbcs(nmax,pmax,statmax,tmax)
real(kind=JPRM) tanfeeds(nmax,pmax,statmax,tmax)
real(kind=JPRM) tfgfeeds(nmax,pmax,statmax,tmax)
integer*1 tflags !(flagmax,nmax,pmax,statmax,tmax)

real(kind=JPRM) lat,lon,pres,oldlat,oldlon

integer wmonrs(statmax)
real(kind=JPRM) wmolats(statmax),wmolons(statmax)

real(kind=JPRM) ps,psvals(pmax),measurements(4)

character*(*) cdatum
character*(*) prefix
character*249 zeile

logical ex

! --- ARRAYS FOR STATION GROUP TABLE
INTEGER, PARAMETER :: JPRM=4
INTEGER, PARAMETER :: JPIM=4
INTEGER(KIND=JPIM)::I_NXGC      , I_NXGP      
PARAMETER (I_NXGC=15000, I_NXGP=3100)
REAL(KIND=JPRM)     , DIMENSION (I_NXGC) :: Z_PLAT,Z_PLON
REAL(KIND=JPRM)     , DIMENSION (I_NXGP) :: Z_QLAT,Z_QLON
CHARACTER(LEN= 6), DIMENSION (I_NXGC) :: CL_CSNG,CL_CSNO
CHARACTER(LEN= 5), DIMENSION (I_NXGC) :: CL_CWMO
CHARACTER(LEN= 6) :: CL_C1,CL_C4,CSTATNR
CHARACTER(LEN= 1) :: CL_C0
CHARACTER(LEN= 5) :: CL_C2
CHARACTER(LEN=19) :: CL_C3
CHARACTER(LEN=21) :: CDATSOURCE
INTEGER(KIND=JPIM)  , DIMENSION (I_NXGP) :: JSGC,JEGC,I_MREP
INTEGER(KIND=JPIM)  , DIMENSION (I_NXGC) :: I_EC,I_NC,I_JM
REAL(KIND=JPRM)    :: LAT4,LON4,Z_R1,Z_R2
INTEGER(KIND=JPIM) :: IOS,I_MSGT
INTEGER(KIND=JPIM) :: I_NGC,I_NGP,I_KM,I_KG,J,I5,I1,I2,ISY,ISM,ISD,ISH,IEY,IEM,IED,IEH,J_EC,J_NC,J_JM,I_KJ,IGOOD,IBAD,IRELOC,IUSED

!       1.2 READ STATION GROUP TABLE

  I_MSGT=21
  OPEN(UNIT=I_MSGT, FILE='stgroup.t', &
   & IOSTAT=IOS,ERR=967,STATUS='OLD')
  I_NGC = 0  ! NUMBER OF RECORD
  I_NGP = 0  ! NUMBER OF GROUP
  350   CONTINUE
  READ(I_MSGT, &
   & '(I5,I4,A1,A6,F7.2,F8.2,A8,2(1X,I4,3I2),I7,23X,A6)', &
   & END=190) &
   & I1,I2,CL_C0,CL_C1,Z_R1,Z_R2,CL_C2, &
   & ISY,ISM,ISD,ISH,IEY,IEM,IED,IEH,I5,CL_C4  

  I_NGC = I_NGC+1
  Z_PLAT(I_NGC) = Z_R1
  Z_PLON(I_NGC) = Z_R2
  CL_CSNO(I_NGC) = CL_C1
!  READ(CDATSOURCE,'(3I1)') J_EC,J_NC,J_JM
!  I_EC(I_NGC)=J_EC
!  I_NC(I_NGC)=J_NC
!  I_JM(I_NGC)=J_JM
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
    CL_C4=CL_CSNG(I_NGC)
    CL_CWMO(I_NGP) = CL_C4(1:5)
    I_MREP(I_NGP) = I_NGC
    Z_QLAT(I_NGP) = Z_R1
    Z_QLON(I_NGP) = Z_R2
  ENDIF
  GOTO 350
  190   CONTINUE

!write(cdatum,'(I4,I2.2)') 1957+mon/12,mod(mon,12)+1

print *,prefix
!print*,wmonrs
iunit=20
inquire(file=prefix//cdatum(1:6),exist=ex)
if(.not.ex) then
  inquire(file=prefix//cdatum(1:6)//'.gz',exist=ex)
  if(.not.ex) then
    call system('ecp ec:feedbacknew/'//prefix//cdatum(1:6)//'.gz .',iret)
    if(iret .ne. 0) then 
      print*,prefix//cdatum//'.gz not found'
      return
    endif
  endif
  call system('gzip -d '//prefix//cdatum(1:6)//'.gz',iret)
  if(iret .ne. 0) then 
    print*,'gzip -d -f '//prefix//cdatum(1:6)//'.gz failed'
    return
  endif
endif


open(iunit,file=prefix//cdatum(1:6),status='old',err=20)
IGOOD=0
IRELOC=0
IBAD=0
IUSED=0
i=0
do while (.true.)
10  continue 
    read(iunit,'(A93)',end=20) zeile

    i=i+1
    if(zeile(12:14) .ne. '101') goto 10
    if(zeile(15:20) .eq. '******') goto 10
    read(zeile(42:49),'(F8.2)') ps
    ip=1
    do while(ip.lt.15 .and. ps.ne.psvals(ip))
      ip=ip+1
    enddo
    if(ip.eq.15) goto 10
!    if(index(zeile,'65.00') .ne. 0 .or.  index(zeile,'67.00') .ne. 0) goto 10
!    print*,ps
!    if(zeile(37:41) .eq. '*****') then 
!       zeile(37:41)=' -999'
!       print *,zeile
!    endif

    read(zeile(1:49),'(I10,I4,I6,2F8.2,I5,F8.2)') time,stattype,statnr,lat,lon,height,pres

    if(statnr .ne. oldstatnr .or. lat .ne. oldlat .or. lon .ne. oldlon) then 
      lat4=lat
      lon4=lon
      WRITE(cstatnr,'(I6.6)') statnr
      call identify_station(cstatnr,CL_CSNO,lat4,lon4,Z_PLAT,Z_PLON,JSGC,JEGC,I_NGP,I_NGC,I_KG,I_KM,I_KJ)
      oldstatnr=statnr
      oldlat=lat
      oldlon=lon
    endif

    IF(I_KG*I_KM .NE. 0) THEN
     if(cstatnr(1:6) .ne. CL_CWMO(I_KG)) THEN
!      WRITE(*,*) cstatnr(1:6),CL_CSNO(I_KJ),CL_CWMO(I_KG)
      READ(CL_CWMO(I_KG),'(I5)') statnr
      lat=Z_QLAT(I_KG)
      lon=Z_QLON(I_KG)
!      WRITE(*,*) 'statnr changed to ',statnr,lat,lon
      IRELOC=IRELOC+1
     ELSE
      IGOOD=IGOOD+1
     ENDIF
    ELSE
!      WRITE(*,*) cstatnr(1:6),' ',CL_CSNO(I_KJ),' ',CL_CWMO(I_KG)
      IBAD=IBAD+1
    ENDIF

    l=1
    do while(wmonrs(l).ne.0 .and. statnr.ne.wmonrs(l))
      l=l+1
    enddo
    if(wmonrs(l).eq.0) goto 10
    
!    if(abs(lat-wmolats(l)) .gt. 2. .or. abs(lon-wmolons(l)) .gt. 2.) then
!      print*,wmonrs(l),lat,lon,wmolats(l),wmolons(l),' discarded'
!      goto 10
!    endif

!    print*,i,statnr,wmonrs(i)

! old ASCII feedback file
!    nz=zeile(46:53)//zeile(46+6*8:53+6*8)//zeile(46+12*8:53+12*8)//&
!       zeile(190:191)//zeile(202:203)//zeile(214:215)//zeile(226:227)//zeile(238:239)
!    read(nz,'(3F8.2,5I2)') measurements(1:3),flags(1:5)

! new ASCII feedback file
    read(zeile(50:93),'(4F8.2,6I2)') measurements,flags

    if(any(flags(1:4) .eq.3) .or. any(measurements .eq. -999.)) then  ! all observations rejected by QC or with no feedback information are not used - blacklisted values are used
      goto 10
    endif

    
    tim100=mod(time,100)

    if(tim100 .le. 3 .or. tim100 .gt. 21) then
     it=1
    else if   (tim100 .gt. 3 .and. tim100 .le. 9) then
     it=2
    else if     (tim100 .gt. 9 .and. tim100 .le. 15) then
     it=3
    else if    (tim100 .gt. 15 .and. tim100 .lt. 21) then
     it=4
    else
       goto 10
    endif
    
    if( (it .eq. 2 .or. it .eq. 4).and. tmax.eq. 2) goto 10

    if(it.eq.3 .and. tmax .eq. 2) it=2
    

      do ik=1,values(ip,l,it)
        if(times(ik,ip,l,it).eq.time .or. times(ik,ip,l,it).eq.time-1 .or.&
           times(ik,ip,l,it).eq.time+1 .or. times(ik,ip,l,it).eq.time-2 .or.&
           times(ik,ip,l,it).eq.time+2 .or. times(ik,ip,l,it).eq.time-77 .or. &
           times(ik,ip,l,it).eq.time+77.or. times(ik,ip,l,it).eq.time-78 .or. &
           times(ik,ip,l,it).eq.time+78) then
!          print*,'duplicate entry',statnr,times(ik,ip,l,it),time,measurements(1)
          goto 10
        endif
      enddo

      if(values(ip,l,it) .lt. 31) then
!        print*,ip,l,it, values(ip,l,it)
        values(ip,l,it)=values(ip,l,it)+1
!        print*,ip,l,it, values(ip,l,it)
      else
        print*,pres,statnr,it,time,'too many values'
        goto 10
      endif
    
!      print*,zeile,ip,l,it,values(ip,l,it)
      times(values(ip,l,it),ip,l,it)=time
      stop 'MAKEBIN: probably wrong date for measurements >22h - please correct me'

      lats(values(ip,l,it),l,it)=lat
      lons(values(ip,l,it),l,it)=lon

      ts(values(ip,l,it),ip,l,it)=measurements(1)

! bias correction is observation_presented_to_analysis-original_observation
      tbcs(values(ip,l,it),ip,l,it)=measurements(2)-measurements(1)
        tfgfeeds(values(ip,l,it),ip,l,it)=-measurements(3) ! fg-obs instead of obs-fg
      tanfeeds(values(ip,l,it),ip,l,it)=-measurements(4) ! an-obs instead of ob

!    tflags(:,values(ip,l,it),ip,l,it)=flags

!    read(zeile(46:249),'(18F8.2,30I2)') measurements,flags

!    if(i .gt. 100) stop
enddo
20 print *,i,' lines read ',IGOOD,' obs in place',IRELOC,' obs relocated',IBAD,' obs not identified'
close(iunit)
call system('/bin/rm  '//prefix//cdatum(1:6),iret)
if(iret .ne. 0) print*,prefix//cdatum(1:6)//' not removed'

return
967 STOP 'COULD NOT OPEN STGROUPS.T'

end subroutine makebin

subroutine makebin_mike(mon,psvals,values,times,statnrs,wmonrs,wmolats,wmolons,lats,lons,ts,tbcs,tanfeeds,tfgfeeds,tflags, &
               statmax,pmax,parmax,nmax,mmax,flagmax,tmax,prefix,cdatum)

! externals: dates

implicit none

integer statmax,pmax,parmax,nmax,mmax,tmax,first,last,flagmax
integer mon,i,iunit,l,height,time,statnr,oldstatnr,it,ip,ik,iret
integer tim100

integer values(pmax,statmax,tmax)
integer times(nmax,pmax,statmax,tmax)
integer statnrs(statmax,tmax)

integer flags(flagmax),stattype

real(kind=JPRM) lats(nmax,statmax,tmax)
real(kind=JPRM) lons(nmax,statmax,tmax)
real(kind=JPRM) ts(nmax,pmax,statmax,tmax)
real(kind=JPRM) tbcs(nmax,pmax,statmax,tmax)
real(kind=JPRM) tanfeeds(nmax,pmax,statmax,tmax)
real(kind=JPRM) tfgfeeds(nmax,pmax,statmax,tmax)
integer*1 tflags !(flagmax,nmax,pmax,statmax,tmax)

real(kind=JPRM) lat,lon,pres,oldlat,oldlon

integer wmonrs(statmax)
real(kind=JPRM) wmolats(statmax),wmolons(statmax)

real(kind=JPRM) ps,psvals(pmax),measurements(4)

character*(*) cdatum
character*(*) prefix
character*249 zeile

logical ex

! --- ARRAYS FOR STATION GROUP TABLE
INTEGER, PARAMETER :: JPRM=4
INTEGER, PARAMETER :: JPIM=4
INTEGER(KIND=JPIM)::I_NXGC      , I_NXGP      
PARAMETER (I_NXGC=15000, I_NXGP=3100)
REAL(KIND=JPRM)     , DIMENSION (I_NXGC) :: Z_PLAT,Z_PLON
REAL(KIND=JPRM)     , DIMENSION (I_NXGP) :: Z_QLAT,Z_QLON
CHARACTER(LEN= 6), DIMENSION (I_NXGC) :: CL_CSNG,CL_CSNO
CHARACTER(LEN= 5), DIMENSION (I_NXGC) :: CL_CWMO
CHARACTER(LEN= 6) :: CL_C1,CL_C4,CSTATNR
CHARACTER(LEN= 1) :: CL_C0
CHARACTER(LEN= 5) :: CL_C2
CHARACTER(LEN=19) :: CL_C3
CHARACTER(LEN=21) :: CDATSOURCE
INTEGER(KIND=JPIM)  , DIMENSION (I_NXGP) :: JSGC,JEGC,I_MREP
INTEGER(KIND=JPIM)  , DIMENSION (I_NXGC) :: I_EC,I_NC,I_JM
REAL(KIND=JPRM)    :: LAT4,LON4,Z_R1,Z_R2
INTEGER(KIND=JPIM) :: IOS,I_MSGT
INTEGER(KIND=JPIM) :: I_NGC,I_NGP,I_KM,I_KG,J,I5,I1,I2,ISY,ISM,ISD,ISH,IEY,IEM,IED,IEH,J_EC,J_NC,J_JM,I_KJ,IGOOD,IBAD,IRELOC,IUSED

!       1.2 READ STATION GROUP TABLE

  I_MSGT=21
  OPEN(UNIT=I_MSGT, FILE='stgroup.t', &
   & IOSTAT=IOS,ERR=967,STATUS='OLD')
  I_NGC = 0  ! NUMBER OF RECORD
  I_NGP = 0  ! NUMBER OF GROUP
  350   CONTINUE
  READ(I_MSGT, &
   & '(I5,I4,A1,A6,F7.2,F8.2,A8,2(1X,I4,3I2),I7,23X,A6)', &
   & END=190) &
   & I1,I2,CL_C0,CL_C1,Z_R1,Z_R2,CL_C2, &
   & ISY,ISM,ISD,ISH,IEY,IEM,IED,IEH,I5,CL_C4  

  I_NGC = I_NGC+1
  Z_PLAT(I_NGC) = Z_R1
  Z_PLON(I_NGC) = Z_R2
  CL_CSNO(I_NGC) = CL_C1
!  READ(CDATSOURCE,'(3I1)') J_EC,J_NC,J_JM
!  I_EC(I_NGC)=J_EC
!  I_NC(I_NGC)=J_NC
!  I_JM(I_NGC)=J_JM
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
    CL_C4=CL_CSNG(I_NGC)
    CL_CWMO(I_NGP) = CL_C4(1:5)
    I_MREP(I_NGP) = I_NGC
    Z_QLAT(I_NGP) = Z_R1
    Z_QLON(I_NGP) = Z_R2
  ENDIF
  GOTO 350
  190   CONTINUE

!write(cdatum,'(I4,I2.2)') 1957+mon/12,mod(mon,12)+1

print *,prefix
!print*,wmonrs
iunit=20
inquire(file=prefix//cdatum(1:6),exist=ex)
if(.not.ex) then
  inquire(file=prefix//cdatum(1:6)//'.gz',exist=ex)
  if(.not.ex) then
    call system('ecp ec:feedbacknew/'//prefix//cdatum(1:6)//'.gz .',iret)
    if(iret .ne. 0) then 
      print*,prefix//cdatum//'.gz not found'
      return
    endif
  endif
  call system('gzip -d '//prefix//cdatum(1:6)//'.gz',iret)
  if(iret .ne. 0) then 
    print*,'gzip -d -f '//prefix//cdatum(1:6)//'.gz failed'
    return
  endif
endif


open(iunit,file=prefix//cdatum(1:6),status='old',err=20)
IGOOD=0
IRELOC=0
IBAD=0
IUSED=0
i=0
do while (.true.)
10  continue 
    read(iunit,'(A93)',end=20) zeile

    i=i+1
    if(zeile(12:14) .ne. '101') goto 10
    if(zeile(15:20) .eq. '******') goto 10
    read(zeile(42:49),'(F8.2)') ps
    ip=1
    do while(ip.lt.15 .and. ps.ne.psvals(ip))
      ip=ip+1
    enddo
    if(ip.eq.15) goto 10
!    if(index(zeile,'65.00') .ne. 0 .or.  index(zeile,'67.00') .ne. 0) goto 10
!    print*,ps
!    if(zeile(37:41) .eq. '*****') then 
!       zeile(37:41)=' -999'
!       print *,zeile
!    endif

    read(zeile(1:49),'(I10,I4,I6,2F8.2,I5,F8.2)') time,stattype,statnr,lat,lon,height,pres

    if(statnr .ne. oldstatnr .or. lat .ne. oldlat .or. lon .ne. oldlon) then 
      lat4=lat
      lon4=lon
      WRITE(cstatnr,'(I6.6)') statnr
      call identify_station(cstatnr,CL_CSNO,lat4,lon4,Z_PLAT,Z_PLON,JSGC,JEGC,I_NGP,I_NGC,I_KG,I_KM,I_KJ)
      oldstatnr=statnr
      oldlat=lat
      oldlon=lon
    endif

    IF(I_KG*I_KM .NE. 0) THEN
     if(cstatnr(1:6) .ne. CL_CWMO(I_KG)) THEN
!      WRITE(*,*) cstatnr(1:6),CL_CSNO(I_KJ),CL_CWMO(I_KG)
      READ(CL_CWMO(I_KG),'(I5)') statnr
      lat=Z_QLAT(I_KG)
      lon=Z_QLON(I_KG)
!      WRITE(*,*) 'statnr changed to ',statnr,lat,lon
      IRELOC=IRELOC+1
     ELSE
      IGOOD=IGOOD+1
     ENDIF
    ELSE
!      WRITE(*,*) cstatnr(1:6),' ',CL_CSNO(I_KJ),' ',CL_CWMO(I_KG)
      IBAD=IBAD+1
    ENDIF

    l=1
    do while(wmonrs(l).ne.0 .and. statnr.ne.wmonrs(l))
      l=l+1
    enddo
    if(wmonrs(l).eq.0) goto 10
    
!    if(abs(lat-wmolats(l)) .gt. 2. .or. abs(lon-wmolons(l)) .gt. 2.) then
!      print*,wmonrs(l),lat,lon,wmolats(l),wmolons(l),' discarded'
!      goto 10
!    endif

!    print*,i,statnr,wmonrs(i)

! old ASCII feedback file
!    nz=zeile(46:53)//zeile(46+6*8:53+6*8)//zeile(46+12*8:53+12*8)//&
!       zeile(190:191)//zeile(202:203)//zeile(214:215)//zeile(226:227)//zeile(238:239)
!    read(nz,'(3F8.2,5I2)') measurements(1:3),flags(1:5)

! new ASCII feedback file
    read(zeile(50:93),'(4F8.2,6I2)') measurements,flags

    if(any(flags(1:4) .eq.3) .or. any(measurements .eq. -999.)) then  ! all observations rejected by QC or with no feedback information are not used - blacklisted values are used
      goto 10
    endif

    
    tim100=mod(time,100)

    if(tim100 .le. 3 .or. tim100 .gt. 21) then
     it=1
    else if   (tim100 .gt. 3 .and. tim100 .le. 9) then
     it=2
    else if     (tim100 .gt. 9 .and. tim100 .le. 15) then
     it=3
    else if    (tim100 .gt. 15 .and. tim100 .lt. 21) then
     it=4
    else
       goto 10
    endif
    
    if( (it .eq. 2 .or. it .eq. 4).and. tmax.eq. 2) goto 10

    if(it.eq.3 .and. tmax .eq. 2) it=2
    

      do ik=1,values(ip,l,it)
        if(times(ik,ip,l,it).eq.time .or. times(ik,ip,l,it).eq.time-1 .or.&
           times(ik,ip,l,it).eq.time+1 .or. times(ik,ip,l,it).eq.time-2 .or.&
           times(ik,ip,l,it).eq.time+2 .or. times(ik,ip,l,it).eq.time-77 .or. &
           times(ik,ip,l,it).eq.time+77.or. times(ik,ip,l,it).eq.time-78 .or. &
           times(ik,ip,l,it).eq.time+78) then
!          print*,'duplicate entry',statnr,times(ik,ip,l,it),time,measurements(1)
          goto 10
        endif
      enddo

      if(values(ip,l,it) .lt. 31) then
!        print*,ip,l,it, values(ip,l,it)
        values(ip,l,it)=values(ip,l,it)+1
!        print*,ip,l,it, values(ip,l,it)
      else
        print*,pres,statnr,it,time,'too many values'
        goto 10
      endif
    
!      print*,zeile,ip,l,it,values(ip,l,it)
      times(values(ip,l,it),ip,l,it)=time
      stop 'MAKEBIN_mike: probably wrong date for measurements >22h - please correct me'

      lats(values(ip,l,it),l,it)=lat
      lons(values(ip,l,it),l,it)=lon

      ts(values(ip,l,it),ip,l,it)=measurements(1)

! bias correction is observation_presented_to_analysis-original_observation
      tbcs(values(ip,l,it),ip,l,it)=measurements(2)-measurements(1)
        tfgfeeds(values(ip,l,it),ip,l,it)=-measurements(3) ! fg-obs instead of obs-fg
      tanfeeds(values(ip,l,it),ip,l,it)=-measurements(4) ! an-obs instead of ob

!    tflags(:,values(ip,l,it),ip,l,it)=flags

!    read(zeile(46:249),'(18F8.2,30I2)') measurements,flags

!    if(i .gt. 100) stop
enddo
20 print *,i,' lines read ',IGOOD,' obs in place',IRELOC,' obs relocated',IBAD,' obs not identified'
close(iunit)
call system('/bin/rm  '//prefix//cdatum(1:6),iret)
if(iret .ne. 0) print*,prefix//cdatum(1:6)//' not removed'

return
967 STOP 'COULD NOT OPEN STGROUPS.T'

end subroutine makebin_mike

subroutine makebinall(first,last,year,month,day,time,psvals,statnr,wmonrs,wmolats,wmolons,ts,tbcs,newtbcs,newradbcs,solarangles,tanfeeds,tfgfeeds,tflags, &
               statmax,pmax,parmax,nmax,indexmax,flagmax,tmax,prefix,lsave)

! externals: dates

implicit none

INTEGER, PARAMETER :: JPIM = 4
INTEGER, PARAMETER :: JPRM= 4

integer statmax,pmax,parmax,nmax,mmax,tmax,first,last,flagmax,indexmax
integer i,iunit,l,height,it,ip,ik,iret,time2,il
integer tim100,day100,year100,mon100,tim101,day101,year101,mon101
integer month(nmax),year(nmax),day(nmax),time(nmax),goodindex(indexmax)

integer statnr,statnr2

integer flags(flagmax),stattype,index

integer times(nmax,1,tmax)

logical lsave ! efficient storage ??

real(kind=JPRM) ts(nmax,pmax,1,tmax)
real(kind=JPRM) tbcs(nmax,pmax,1,tmax)
real(kind=JPRM) newtbcs(nmax,pmax,1,tmax)
real(kind=JPRM) newradbcs(nmax,pmax,1,tmax)
real(kind=JPRM) solarangles(nmax,1,tmax)
real(kind=JPRM) tanfeeds(nmax,pmax,1,tmax)
real(kind=JPRM) tfgfeeds(nmax,pmax,1,tmax)
integer*1 tflags !(flagmax,nmax,pmax,statmax,tmax)

real(kind=JPRM) lat,lon,pres

integer wmonrs(statmax)
real(kind=JPRM) wmolats(statmax),wmolons(statmax)

real(kind=JPRM) ps,psvals(pmax),measurements(4)

character*6   cstatnr
character*(*) prefix
character*249 zeile
character*1 czeit
character*10 cdatum

logical ex

REAL(KIND=JPRM) :: THILF(pmax)
REAL(KIND=JPRM) :: NEWTBCHILF(pmax),NEWRADBCHILF(pmax)
real(KIND=JPRM) :: lats(nmax,1,tmax)
real(KIND=JPRM) :: lons(nmax,1,tmax)

LOGICAL              :: LD_LBC = .FALSE.,LL_FILES_NOT_READ
INTEGER(KIND=JPIM)   :: K_IY,K_IY_OLD
INTEGER(KIND=JPIM)   :: K_IM
INTEGER(KIND=JPIM)   :: K_ID
INTEGER(KIND=JPIM)   :: K_IH
REAL(KIND=JPRM)      :: P_RMISS = -999.,psvals_pa(pmax)
INTEGER(KIND=JPIM)   :: I_MTBL,I_MCOR,I_MCORRAD,I_MSGT,IOS


times=0

psvals_pa=psvals*100.
write(cstatnr,'(I6.6)') statnr

print *,prefix
print*,'x'//prefix//cstatnr//'x'
!print*,wmonrs
iunit=20
i=0
inquire(file=prefix//cstatnr,exist=ex)
if(.not.ex) then
  inquire(file=prefix//cstatnr//'.gz',exist=ex)
  if(.not.ex) then
    call system('ecp ec:bufrseriesguan/'//prefix//cstatnr//'.gz .',iret)
    if(iret .ne. 0) then 
      print*,prefix//cstatnr//'.gz not found'
      return
    endif
  endif
  call system('gzip -d '//prefix//cstatnr//'.gz',iret)
  if(iret .ne. 0) then 
    print*,'gzip -d -f '//prefix//cstatnr//'.gz failed'
    return
  endif
endif

open(iunit,file=prefix//cstatnr,status='old',err=20)

index=1

    l=1
do while (.true.)
10  continue 
    read(iunit,'(A93)',end=20) zeile

    if(zeile(12:14) .ne. '101') goto 10
    if(zeile(15:20) .eq. '******') goto 10
    read(zeile(42:49),'(F8.2)') ps
    ip=1
    do while(ip.le.pmax .and. ps.ne.psvals(ip))
      ip=ip+1
    enddo
    if(ip.eq.pmax+1) goto 10

    read(zeile(1:49),'(I10,I4,I6,2F8.2,I5,F8.2)') time2,stattype,statnr2,lat,lon,height,pres

!    if(statnr2 .ne. statnr) stop 'wrong station number'

!    do while(wmonrs(l).ne.0 .and. statnr.ne.wmonrs(l))
!      l=l+1
!    enddo
!    if(wmonrs(l).eq.0) goto 10
!    if(abs(lat-wmolats(l)) .gt. 2. .or. abs(lon-wmolons(l)) .gt. 2.) then
!      print*,wmonrs(l),lat,lon,wmolats(l),wmolons(l),' discarded'
!      goto 10
!    endif

!    print*,i,statnr,wmonrs(i)

! old ASCII feedback file
!    nz=zeile(46:53)//zeile(46+6*8:53+6*8)//zeile(46+12*8:53+12*8)//&
!       zeile(190:191)//zeile(202:203)//zeile(214:215)//zeile(226:227)//zeile(238:239)
!    read(nz,'(3F8.2,5I2)') measurements(1:3),flags(1:5)

! new ASCII feedback file
    read(zeile(50:93),'(4F8.2,6I2)') measurements,flags

    if(any(flags(1:4) .eq.3) .or. any(measurements .eq. -999.)) then  ! all observations rejected by QC or with no feedback information are not used - blacklisted values are used
      goto 10
    endif

    year100=time2/1000000
    mon100=mod(time2/10000,100)
    day100=mod(time2/100,100)
    tim100=mod(time2,100)

    if(tim100 .le. 3 .or. tim100 .gt. 21) then
     it=1
     if(tim100 .gt. 21) then
       call dates(year100,mon100,day100,tim100,24,year101,mon101,day101,tim101)
       year100=year101
       mon100=mon101
       day100=day101
     endif
    else if   (tim100 .gt. 3 .and. tim100 .le. 9) then
     it=2
    else if     (tim100 .gt. 9 .and. tim100 .le. 15) then
     it=3
    else if    (tim100 .gt. 15 .and. tim100 .lt. 21) then
     it=4
    else
       goto 10
    endif
    
    if( (it .eq. 2 .or. it .eq. 4).and. tmax.eq. 2) goto 10

    if(it.eq.3 .and. tmax .eq. 2) it=2
    
    do while(year(index) .lt. year100)
      index=index+1
    enddo
    do while(month(index) .lt. mon100)
      index=index+1
    enddo
    do while(day(index) .lt. day100)
      index=index+1
    enddo

    times(index,l,it)=time2
    lats(index,l,it)=lat
    lons(index,l,it)=lon

    ts(index,ip,l,it)=measurements(1)
      stop 'MAKEBINALL: probably wrong date for measurements >22h - please correct me'

! bias correction is observation_presented_to_analysis-original_observation
      tbcs(index,ip,l,it)=measurements(2)-measurements(1)
      tfgfeeds(index,ip,l,it)=-measurements(3) ! fg-obs instead of obs-fg
      tanfeeds(index,ip,l,it)=-measurements(4) ! an-obs instead of ob

  i=i+1
enddo
20 print *,i,' lines read'
close(iunit)

I_MTBL=20
I_MCOR=21
I_MCORRAD=23
I_MSGT=22
i=0


do it=1,tmax
  K_IY_OLD=1957
  do index=1,indexmax
   IF(TIMES(index,1,IT).GT.0) THEN
    K_IY=TIMES(index,1,IT)/1000000
    K_IM=mod(TIMES(index,1,IT)/10000,100)
    K_ID=mod(TIMES(index,1,IT)/100,100)
    K_IH=mod(TIMES(index,1,IT),100)
    IF(K_IY_OLD .lt. K_IY) THEN
      write(cdatum,'(I10)') TIMES(index,1,IT)

      write(*,*) it,' ',cdatum
!      OPEN(UNIT=I_MTBL, FILE='/era/work/eru/recalc_biases_era40/radiation_only/scr/country.t', &
      OPEN(UNIT=I_MTBL, FILE='country.t', &
   & IOSTAT=IOS,ERR=967,STATUS='OLD')  
!      write(*,*) 'country.t opened'

!      OPEN(UNIT=I_MCOR, FILE='/era/work/eru/recalc_biases_era40/radiation_only/tab/corcand_'//cdatum(1:4)//'.t', &
      OPEN(UNIT=I_MCOR, FILE='rad_and_mean/corcand_'//cdatum(1:4)//'.t', &
   & IOSTAT=IOS,ERR=967,STATUS='OLD')  
!      write(*,*) 'rad_and_mean/corcand_'//cdatum(1:4)//'.t opened'

!  OPEN(UNIT=I_MCORRAD, FILE='/era/work/erl/recalc_biases_era40/radiation_only/tab/corcand_'//cdatum(1:4)//'.t', &
      OPEN(UNIT=I_MCORRAD, FILE='radiation_only/corcand_'//cdatum(1:4)//'.t', &
   & IOSTAT=IOS,ERR=967,STATUS='OLD')  
!      write(*,*) 'radiation_only/corcand_'//cdatum(1:4)//'.t opened'

!      OPEN(UNIT=I_MSGT, FILE='/era/work/eru/recalc_biases_era40/radiation_only/scr/stgroup.t', &
      OPEN(UNIT=I_MSGT, FILE='stgroup.t', &
   & IOSTAT=IOS,ERR=967,STATUS='OLD')
!      write(*,*) 'stgroup.t opened'

      LL_FILES_NOT_READ=.TRUE.
!      WRITE(*,*) cdatum
    ENDIF

    K_IY_OLD=K_IY
    THILF=TS(index,:,1,it)
    CALL BIASCOR_ERA40 (LD_LBC,pmax,PSVALS_PA,THILF,NEWTBCHILF,NEWRADBCHILF,SOLARANGLES(index,1,it), &
   CSTATNR,K_IM,K_ID,K_IH,LATS(INDEX,1,IT),LONS(INDEX,1,IT),&
   P_RMISS,I_MTBL,I_MCOR,I_MCORRAD,I_MSGT,LL_FILES_NOT_READ)
    NEWTBCS(index,:,1,it)=NEWTBCHILF-THILF
    NEWRADBCS(index,:,1,it)=NEWRADBCHILF-THILF
    IF (LD_LBC .and. .false.) THEN
        WRITE(6,'('' RS TEMPERATURE BIAS CORRECTION DONE FOR:'')')
        WRITE(6,'(''    LAT, LON, STID: '',&
       & 2(F10.2,1X),A)') LATS(INDEX,1,IT),LONS(INDEX,1,IT),CSTATNR,K_IM,K_ID,K_IH
        WRITE(6,'(16F7.1)')TS(INDEX,:,1,IT)
        WRITE(6,'(16F7.1)')NEWTBCS(INDEX,:,1,IT)
    ENDIF
   ENDIF
  ENDDO
ENDDO

if(lsave) then
  il=1
  do i=1,indexmax
    if(any(ts(i,:,:,:) .gt. -999.)) then
      goodindex(il)=i
      ts(il,:,:,:)=ts(i,:,:,:)
      tanfeeds(il,:,:,:)=tanfeeds(i,:,:,:)
      tfgfeeds(il,:,:,:)=tfgfeeds(i,:,:,:)
      tbcs(il,:,:,:)=tbcs(i,:,:,:)
      solarangles(il,:,:)=solarangles(i,:,:)
      il=il+1
    endif
  enddo
  il=il-1
  open(iunit,file=prefix//'binsave'//cstatnr,form='unformatted')
  write(iunit) indexmax,il,pmax,tmax
  write(iunit) goodindex(1:il)
  write(iunit) ts(1:il,:,1,:)
  write(iunit) tanfeeds(1:il,:,1,:)
  write(iunit) tfgfeeds(1:il,:,1,:)
  write(iunit) solarangles(1:il,1,:)
  write(iunit) tbcs(1:il,:,1,:)

! write(iunit) newtbcs(1:indexmax,:,1,:),newradbcs(1:indexmax,:,1,:)
  close(iunit)
else
  open(iunit,file=prefix//'bin'//cstatnr,form='unformatted')
  write(iunit) indexmax,pmax,tmax
  write(iunit) ts(1:indexmax,:,1,:),tanfeeds(1:indexmax,:,1,:),tfgfeeds(1:indexmax,:,1,:),solarangles(1:indexmax,1,:),tbcs(1:indexmax,:,1,:)
! write(iunit) newtbcs(1:indexmax,:,1,:),newradbcs(1:indexmax,:,1,:)
  close(iunit)
endif

call system('/bin/rm  '//prefix//cstatnr,iret)
if(iret .ne. 0) print*,prefix//cstatnr//' not removed'

return
967 write(*,*) 'error opening bias correction tables'
call exit(1)
return
end subroutine makebinall

subroutine makebinallgeopot(first,last,year,month,day,time,psvals,statnr,wmonrs,wmolats,wmolons,ts,solarangles,tanfeeds,tfgfeeds,tflags, &
               statmax,pmax,parmax,nmax,indexmax,flagmax,tmax,prefix)

! externals: dates

implicit none

INTEGER, PARAMETER :: JPIM = 4
INTEGER, PARAMETER :: JPRM= 4

integer statmax,pmax,parmax,nmax,mmax,tmax,first,last,flagmax,indexmax
integer i,iunit,l,height,it,ip,ik,iret,time2
integer tim100,day100,year100,mon100,tim101,day101,year101,mon101
integer month(nmax),year(nmax),day(nmax),time(nmax)

 integer statnr,statnr2

integer flags(flagmax),stattype,index

integer times(nmax,1,tmax)

real(kind=JPRM) ts(nmax,pmax,1,tmax)
!real(kind=JPRM) tbcs(nmax,pmax,1,tmax)
!real(kind=JPRM) newtbcs(nmax,pmax,1,tmax)
!real(kind=JPRM) newradbcs(nmax,pmax,1,tmax)
real(kind=JPRM) solarangles(nmax,1,tmax)
real(kind=JPRM) tanfeeds(nmax,pmax,1,tmax)
real(kind=JPRM) tfgfeeds(nmax,pmax,1,tmax)
integer*1 tflags !(flagmax,nmax,pmax,statmax,tmax)

real(kind=JPRM) lat,lon,pres

integer wmonrs(statmax)
real(kind=JPRM) wmolats(statmax),wmolons(statmax)

real(kind=JPRM) ps,psvals(pmax),measurements(4),obssource

character*6   cstatnr
character*(*) prefix
character*280 zeile
character*1 czeit
character*10 cdatum
character*40 hilf

logical ex

REAL(KIND=JPRM) :: THILF(pmax)
REAL(KIND=JPRM) :: NEWTBCHILF(pmax),NEWRADBCHILF(pmax)
real(KIND=JPRM) :: lats(nmax,1,tmax)
real(KIND=JPRM) :: lons(nmax,1,tmax)

LOGICAL              :: LD_LBC = .FALSE.,LL_FILES_NOT_READ
INTEGER(KIND=JPIM)   :: K_IY,K_IY_OLD
INTEGER(KIND=JPIM)   :: K_IM
INTEGER(KIND=JPIM)   :: K_ID
INTEGER(KIND=JPIM)   :: K_IH
REAL(KIND=JPRM)      :: P_RMISS = -999.,psvals_pa(pmax)
INTEGER(KIND=JPIM)   :: I_MTBL,I_MCOR,I_MCORRAD,I_MSGT,IOS


times=0

psvals_pa=psvals*100.
write(cstatnr,'(I6.6)') statnr

print *,prefix
print*,'x'//prefix//cstatnr//'x'
!print*,wmonrs
iunit=20
i=0
inquire(file=prefix//cstatnr,exist=ex)
if(.not.ex) then
  inquire(file=prefix//cstatnr//'.gz',exist=ex)
  if(.not.ex) then
    call system('ecp ec:bufrseriesguan/'//prefix//cstatnr//'.gz .',iret)
    if(iret .ne. 0) then 
      print*,prefix//cstatnr//'.gz not found'
      return
    endif
  endif
  call system('gzip -d '//prefix//cstatnr//'.gz',iret)
  if(iret .ne. 0) then 
    print*,'gzip -d -f '//prefix//cstatnr//'.gz failed'
    return
  endif
endif

open(iunit,file=prefix//cstatnr,status='old',err=20)

index=1

    l=1
do while (.true.)
10  continue 
    read(iunit,'(A280)',end=20) zeile

    if(zeile(12:14) .ne. '101') goto 10
    if(zeile(18:23) .eq. '******') goto 10
    read(zeile(45:52),'(F8.2)') ps
    ip=1
    do while(ip.le.pmax .and. ps.ne.psvals(ip))
      ip=ip+1
    enddo
    if(ip.eq.pmax+1) goto 10

    read(zeile(1:52),'(I10,I4,I3,I6,2F8.2,I5,F8.2)') time2,stattype,obssource,statnr2,lat,lon,height,pres

! write statement in fbdecodemike
!
!                 WRITE(10,
!     *'(I4,3I2.2,I4,I3,I6,2F8.2,I6.6F8.2,6I2,2(3F8.2,5I2),
!     *7F8.2,5I2,4F10.2,5I2)')
!     *nint(values(5)),nint(values(6)),nint(values(7)),nint(values(8)),
!     *ksec1(7),ksec1(21),nint(values(1)*1000+values(2)),values(10),
!     *values(11),nint(values(12)),pres(2,il)/100.,t(1:4,il),
!     *nint(t(6:11,il)),sh(2:4,il)*10000,nint(sh(6:10,il)),
!     *rh(2:4,il),nint(rh(6:10,il)),ws(1,il),
!     *u(2:4,il),v(2:4,il),
!     *nint(u(6:10,il)),geopot(1:4,il),
!     *nint(geopot(6:10,il))

! Read temperature information, we use only flags here to find spot bad temperature
! information which indicates also erroneous geopotential information. 
! Geopotential flags useless since geopotential is blacklisted
    read(zeile(53:96),'(4F8.2,6I2)') measurements,flags

    if(any(flags(1:4) .eq.3) .or. any(measurements .eq. -999.)) then  ! all observations rejected by QC or with no feedback information are not used - blacklisted values are used
      goto 10
    endif
! Read geopotential information
    hilf=zeile(231:270)
    read(hilf,'(4F10.2)') measurements
    if(any(flags(1:4) .eq.3) .or. any(measurements .eq. -999.)) then  ! all observations rejected by QC or with no feedback information are not used - blacklisted values are used
      goto 10
    endif

! measurenments now contains geopotential or height information

    year100=time2/1000000
    mon100=mod(time2/10000,100)
    day100=mod(time2/100,100)
    tim100=mod(time2,100)

    if(tim100 .le. 3 .or. tim100 .gt. 21) then
     it=1
     if(tim100 .gt. 21) then
       call dates(year100,mon100,day100,tim100,24,year101,mon101,day101,tim101)
       year100=year101
       mon100=mon101
       day100=day101
     endif
    else if   (tim100 .gt. 3 .and. tim100 .le. 9) then
     it=2
    else if     (tim100 .gt. 9 .and. tim100 .le. 15) then
     it=3
    else if    (tim100 .gt. 15 .and. tim100 .lt. 21) then
     it=4
    else
       goto 10
    endif
    
    if( (it .eq. 2 .or. it .eq. 4).and. tmax.eq. 2) goto 10

    if(it.eq.3 .and. tmax .eq. 2) it=2
    
    do while(year(index) .lt. year100)
      index=index+1
    enddo
    do while(month(index) .lt. mon100)
      index=index+1
    enddo
    do while(day(index) .lt. day100)
      index=index+1
    enddo

    times(index,l,it)=time2
    lats(index,l,it)=lat
    lons(index,l,it)=lon

    ts(index,ip,l,it)=measurements(2) ! we want height, not geopotential
      stop 'MAKEBIN: probably wrong date for measurements >22h - please correct me'

      tfgfeeds(index,ip,l,it)=-measurements(3) ! fg-obs instead of obs-fg
      tanfeeds(index,ip,l,it)=-measurements(4) ! an-obs instead of ob

  i=i+1
enddo
20 print *,i,' lines read'
close(iunit)

I_MTBL=20
I_MCOR=21
I_MCORRAD=23
I_MSGT=22
i=0


do it=1,tmax
  K_IY_OLD=1957
  do index=1,indexmax
   IF(TIMES(index,1,IT).GT.0) THEN
    K_IY=TIMES(index,1,IT)/1000000
    K_IM=mod(TIMES(index,1,IT)/10000,100)
    K_ID=mod(TIMES(index,1,IT)/100,100)
    K_IH=mod(TIMES(index,1,IT),100)
    IF(K_IY_OLD .lt. K_IY) THEN
      write(cdatum,'(I10)') TIMES(index,1,IT)

      write(*,*) it,' ',cdatum
!      OPEN(UNIT=I_MTBL, FILE='/era/work/eru/recalc_biases_era40/radiation_only/scr/country.t', &
      OPEN(UNIT=I_MTBL, FILE='country.t', &
   & IOSTAT=IOS,ERR=967,STATUS='OLD')  
!      write(*,*) 'country.t opened'

!      OPEN(UNIT=I_MCOR, FILE='/era/work/eru/recalc_biases_era40/radiation_only/tab/corcand_'//cdatum(1:4)//'.t', &
      OPEN(UNIT=I_MCOR, FILE='rad_and_mean/corcand_'//cdatum(1:4)//'.t', &
   & IOSTAT=IOS,ERR=967,STATUS='OLD')  
!      write(*,*) 'rad_and_mean/corcand_'//cdatum(1:4)//'.t opened'

!  OPEN(UNIT=I_MCORRAD, FILE='/era/work/erl/recalc_biases_era40/radiation_only/tab/corcand_'//cdatum(1:4)//'.t', &
      OPEN(UNIT=I_MCORRAD, FILE='radiation_only/corcand_'//cdatum(1:4)//'.t', &
   & IOSTAT=IOS,ERR=967,STATUS='OLD')  
!      write(*,*) 'radiation_only/corcand_'//cdatum(1:4)//'.t opened'

!      OPEN(UNIT=I_MSGT, FILE='/era/work/eru/recalc_biases_era40/radiation_only/scr/stgroup.t', &
      OPEN(UNIT=I_MSGT, FILE='stgroup.t', &
   & IOSTAT=IOS,ERR=967,STATUS='OLD')
!      write(*,*) 'stgroup.t opened'

      LL_FILES_NOT_READ=.TRUE.
!      WRITE(*,*) cdatum
    ENDIF

    K_IY_OLD=K_IY
    THILF=TS(index,:,1,it)
    CALL BIASCOR_ERA40 (LD_LBC,pmax,PSVALS_PA,THILF,NEWTBCHILF,NEWRADBCHILF,SOLARANGLES(index,1,it), &
   CSTATNR,K_IM,K_ID,K_IH,LATS(INDEX,1,IT),LONS(INDEX,1,IT),&
   P_RMISS,I_MTBL,I_MCOR,I_MCORRAD,I_MSGT,LL_FILES_NOT_READ)
!    NEWTBCS(index,:,1,it)=NEWTBCHILF-THILF
!    NEWRADBCS(index,:,1,it)=NEWRADBCHILF-THILF
    IF (LD_LBC .and. .false.) THEN
        WRITE(6,'('' RS TEMPERATURE BIAS CORRECTION DONE FOR:'')')
        WRITE(6,'(''    LAT, LON, STID: '',&
       & 2(F10.2,1X),A)') LATS(INDEX,1,IT),LONS(INDEX,1,IT),CSTATNR,K_IM,K_ID,K_IH
        WRITE(6,'(16F7.1)')TS(INDEX,:,1,IT)
!        WRITE(6,'(16F7.1)')NEWTBCS(INDEX,:,1,IT)
    ENDIF
   ENDIF
  ENDDO
ENDDO

open(iunit,file=prefix//'bin'//cstatnr,form='unformatted')
write(iunit) indexmax,pmax,tmax
write(iunit) ts(1:indexmax,:,1,:),tanfeeds(1:indexmax,:,1,:),tfgfeeds(1:indexmax,:,1,:),solarangles(1:indexmax,1,:)!,tbcs(1:indexmax,:,1,:)
! write(iunit) newtbcs(1:indexmax,:,1,:),newradbcs(1:indexmax,:,1,:)
close(iunit)
call system('/bin/rm  '//prefix//cstatnr,iret)
if(iret .ne. 0) print*,prefix//cstatnr//' not removed'

return
967 write(*,*) 'error opening bias correction tables'
call exit(1)
return
end subroutine makebinallgeopot

subroutine read_bg_daily(iunit,filename,ni,pmax,parmax,tfgm,miss_val,err) !

implicit none

integer ni,pmax,parmax,ini,il,ipmax,iparmax,i,iunit,err,imax, ios 
integer goodindex(ni)
character*(*) filename
real(kind=JPRM) :: tfgm(ni,pmax,parmax),miss_val
real(kind=JPRM),allocatable :: hilf(:,:,:)


  tfgm=miss_val

  iunit=iunit
  open(iunit,file=filename,form='unformatted',status='old',action='read',iostat=ios, err=120)
  err=0
  read(iunit,err=120,end=120) ini,il,ipmax,iparmax


  if(il.eq.0) then 
    close(iunit)
    return
  endif
  read(iunit,err=120,end=120) goodindex(1:il)

  imax=0
  do i=1,il
    if(goodindex(i) .le. ni) imax=imax+1
  enddo

  allocate(hilf(il,ipmax,iparmax))

  read(iunit,err=119,end=119) hilf
  do i=1,imax
    tfgm(goodindex(i),:,:)=hilf(i,:,:)
  enddo


  close(iunit)

  deallocate(hilf)
  return

119 print*,'could not read file ',filename, err
  deallocate(hilf)
120 print*,'could not open/read file ',filename, err, ios
  err=1

  return

end subroutine read_bg_daily


subroutine read_igrasonde_daily(iunit,filename,ni,pmax,parmax,tm,tanm,tfgm,t12m,miss_val,ilmin,err)

implicit none

integer ni,pmax,parmax,ini,il,ipmax,iparmax,i,iunit,err,imax,ilmin
integer goodindex(ni)
character*(*) filename
real(kind=JPRM) :: tm(ni,pmax,parmax),tanm(ni,pmax,parmax),tfgm(ni,pmax,parmax),t12m(ni,pmax,parmax),miss_val
real(kind=JPRM),allocatable :: hilf(:,:,:)


  iunit=iunit
  open(iunit,file=filename,form='unformatted',status='old',action='read',err=120)
  err=0
  read(iunit,err=120,end=120) ini,il,ipmax,iparmax


  if(il.le.ilmin) then 
    close(iunit)
    return
  endif

  tm=miss_val
  tanm=miss_val
  tfgm=miss_val
  t12m=miss_val

  read(iunit,err=120,end=120) goodindex(1:il)

  imax=0
  do i=1,il
    if(goodindex(i) .le. ni) imax=imax+1
  enddo

  allocate(hilf(il,ipmax,iparmax))

  read(iunit,err=119,end=119) hilf
  do i=1,imax
    tm(goodindex(i),:,:)=hilf(i,:,:)
  enddo

  read(iunit,err=119,end=119) hilf 
  do i=1,imax
    tanm(goodindex(i),:,:)=hilf(i,:,:)
  enddo

  read(iunit,err=119,end=119) hilf
  do i=1,imax
    tfgm(goodindex(i),:,:)=hilf(i,:,:)
  enddo

  read(iunit,err=119,end=119) hilf
  do i=1,imax
    t12m(goodindex(i),:,:)=hilf(i,:,:)
  enddo


  close(iunit)

  deallocate(hilf)
  return

119 print*,'could not read file ',filename
  deallocate(hilf)
120 print*,'could not open/read file ',filename
  err=1

  return

end subroutine read_igrasonde_daily

subroutine read_igrasonde_daily_new(iunit,filename,ni,pmax,parmax,tm,tanm,tfgm,t12m,miss_val,ilmin,err) !

implicit none

integer ni,pmax,parmax,ini,ipmax,iparmax,i,iunit,err,imax,ipar,ip,ios,ilmin
integer goodindex(ni),bi(parmax),count
character*80 filename
real(kind=JPRM) :: tm(ni,pmax,parmax),tanm(ni,pmax,parmax),tfgm(ni,pmax,parmax),t12m(ni,pmax,parmax),miss_val
real(kind=JPRM),allocatable :: hilf(:,:)



  count=0
11  open(iunit,file=filename,form='unformatted',status='old',action='read',iostat=ios)
  if(ios .ne. 0) then
    if(ios .eq. 4052) then
      count=count+1
      goto 11
    else
      goto 120
    endif
  endif
  err=0

  read(iunit,iostat=ios) ini,ipmax,iparmax
  read(iunit,iostat=ios) bi(1:iparmax)

  if(bi(1) .le. ilmin .and. bi(2) .le. ilmin) then

   write(*,'(A11,I4,A14,A80,I2)'),'fewer than ',ilmin,'days found in ',filename, err
   err=2

   return

  endif

 
  tm=miss_val
  tanm=miss_val
  tfgm=miss_val
  t12m=miss_val

  do ipar=1,iparmax

  allocate(hilf(bi(ipar),ipmax))

  read(iunit,iostat=ios) goodindex(1:bi(ipar))

  imax=0
  do i=1,bi(ipar)
    if(goodindex(i) .le. ni) imax=imax+1
  enddo

  read(iunit,iostat=ios) hilf
  do ip=1,pmax
  do i=1,imax
    tm(goodindex(i),ip,ipar)=hilf(i,ip)
  enddo
  enddo

  read(iunit,iostat=ios) hilf
  do ip=1,pmax
  do i=1,imax
    tanm(goodindex(i),ip,ipar)=hilf(i,ip)
  enddo
  enddo

  read(iunit,iostat=ios) hilf
  do ip=1,pmax
  do i=1,imax
    tfgm(goodindex(i),ip,ipar)=hilf(i,ip)
  enddo
  enddo
  if(trim(filename) .eq. './feedbackmerged10307') then
    print*,'found'
  endif

  read(iunit,iostat=ios) hilf
  do ip=1,pmax
  do i=1,imax
    t12m(goodindex(i),ip,ipar)=hilf(i,ip)
  enddo
  enddo

  if(ios .ne. 0) then
  !!$ call omp_set_lock(omp_lp)
  write(*,*) filename,' error in read',ios
  call exit(1)
  !!$ call omp_unset_lock(omp_lp)
  endif
  deallocate(hilf)
 
  enddo

  close(iunit)

    do ipar=1,parmax
    do ip=1,pmax
  do i=1,ni
    if(tfgm(i,ip,ipar) .ne. miss_val .and. abs(tfgm(i,ip,ipar)) .gt. 30.) then
     !!$ call omp_set_lock(omp_lp)
       write(*,*) filename,i,ip,ipar,' erroneous ',tfgm(i,ip,ipar),tm(i,ip,ipar)
     !!$ call omp_unset_lock(omp_lp)
     tfgm(i,ip,ipar)=miss_val
     tm(i,ip,ipar)=miss_val
    endif 
  enddo
  enddo
  enddo
  

  return

120 continue
  !!$ call omp_set_lock(omp_lp)
  print*,'read_igrasonde_daily_new: could not open file ',filename,ios
  !!$ call omp_unset_lock(omp_lp)
  err=1

  return

end subroutine read_igrasonde_daily_new

subroutine cards(wmonrs,wmolats,wmolons,wmonames,statmax,wmostats) !

implicit none

character clat*9,clon*8,cstat*20,cdum*80,cdum2*80
integer iunit,istat,obs,statmax,i,l,wmostats
integer wmonrs(statmax)
real(kind=JPRM) wmolats(statmax),wmolons(statmax),lat,lon
character*20 wmonames(statmax)
logical lmissing
integer missing(51)

data missing /1661,3693,8301,8508,8522,14240,15730,16320,29839,33317,40948,41923,43353,47881,57290,59553,60760,61415,62721,67341,80413,82599,82900,89001, &
61998,63450,67964,71913,72293,85469,89022,89050,89564,61901,76654,78583,78762,82397,84008,85934,89002,89055,89512,89592,89642,91557,91801,92035,96315,96935,98223/

iunit=20
open(iunit,file='cards_stations.txt',status='old',err=30)
do i=1,10
  read(iunit,'(A1)') cdum2
enddo
l=1
do while (.true.)
  read(iunit,'(I5,7X,A20,18X,A8,2X,A9,24X,I5)',end=20) istat,cstat,clat,clon,obs

  read(clat,'(F8.4)') lat
  read(clon,'(F9.4)') lon

  lmissing=.false.
  do i=1,51
    if(istat .eq. missing(i)) lmissing=.true.
  enddo

  if(lat .ne. 0 .and. lon .ne. 0 .and. obs .gt. 10000 .or. lmissing) then 
  wmonrs(l)=istat
  wmolats(l)=lat
  wmolons(l)=lon
  wmonames(l)=cstat
!  write(*,'(I3,I6,F8.2,F8.2,A20)') l,wmonrs(l),wmolats(l),wmolons(l),wmonames(l)

  l=l+1
  endif
  
enddo

20 close(iunit)
wmostats=l

return

30 write(*,*) 'error opening cards_stations.txt'
call exit(1)

end subroutine cards

subroutine igra(wmonrs,wmolats,wmolons,wmonames,statmax,wmostats) !

implicit none

character clat*9,clon*8,cstat*35,cdum*80,cdum2*80
integer iunit,istat,obs,statmax,i,l,wmostats,height,icountry
integer wmonrs(statmax)
real(kind=JPRM) wmolats(statmax),wmolons(statmax),lat,lon
character*20 wmonames(statmax)
!character*(*) prefix
logical lmissing
!integer missing(51)

!data missing /1661,3693,8301,8508,8522,14240,15730,16320,29839,33317,40948,41923,43353,47881,57290,59553,60760,61415,62721,67341,80413,82599,82900,89001, &
!61998,63450,67964,71913,72293,85469,89022,89050,89564,61901,76654,78583,78762,82397,84008,85934,89002,89055,89512,89592,89642,91557,91801,92035,96315,96935,98223/

iunit=20
open(iunit,file='igra.station.list',status='old',err=30)

l=0
do while (.true.)
  read(iunit,'(I3,1X,I5,2X,A35,F7.2,F8.2,I5)',end=20) icountry,istat,cstat,lat,lon,height


  if(lat .ne. 0 .or. lon .ne. 0) then 
  l=l+1
  wmonrs(l)=istat
  wmolats(l)=lat
  wmolons(l)=lon
  wmonames(l)=cstat
!  write(*,'(I3,I6,F8.2,F8.2,A20)') l,wmonrs(l),wmolats(l),wmolons(l),wmonames(l)

  endif
  
enddo

20 close(iunit)
wmostats=l

return

30 write(*,*) 'error opening '//'igra.station.list'
call exit(1)

end subroutine igra

subroutine sphdist(lat,lon,tlats,tlons,dists,statmax) !

implicit none

integer i,j,statmax
real(kind=JPRM) :: lat,lon,tlats(statmax),tlons(statmax)

real :: sinlats(statmax),coslats(statmax),sinlons(statmax),coslons(statmax),vecs(3,statmax),dists(statmax),vec(3),pi,dp

  pi=acos(-1.d0)
  coslats=cos(tlats*pi/180)
  sinlats=sin(tlats*pi/180)
  coslons=cos(tlons*pi/180)
  sinlons=sin(tlons*pi/180)

  vecs(1,:)=coslats*coslons
  vecs(2,:)=coslats*sinlons
  vecs(3,:)=sinlats


vec(1)=cos(lat*pi/180)*cos(lon*pi/180)
vec(2)=cos(lat*pi/180)*sin(lon*pi/180)
vec(3)=sin(lat*pi/180)

do i=1,statmax
  dp=dot_product(vec,vecs(:,i))
  if(dp .gt. 1.) dp=1.
  if(dp .lt. -1.) dp=-1.
  dists(i)=acos(dp)/pi*180.
enddo

return
end subroutine sphdist

subroutine sphdistlat(lat,lon,tlats,tlons,dists,statmax)

implicit none

integer i,j,statmax
real(kind=JPRM) :: lat,lon,tlats(statmax),tlons(statmax)

real :: londiff,latdiff,dists(statmax),vec(3),pi,dp

pi=acos(-1.d0)
do i=1,statmax
  londiff=abs(lon-tlons(i))
  if(londiff .gt. 180.) londiff=360.-londiff
  latdiff=abs(lat-tlats(i))
  dists(i)=latdiff+0.1*londiff*cos((lat+tlats(i))/2.*pi/180.)
enddo

return
end subroutine sphdistlat

subroutine anomaly(tin,tout,tmonmean,mmax,mdiff,miss_val,thresh) !
!
! this subroutine calculates anomalies from the long-term monthly means
!

real(kind=JPRM) tin(mmax),tout(mmax)
real(kind=JPRM) tmonmean(mdiff),miss_val
integer im,iy,il,tmonval(mdiff),thresh

if(mdiff.ne.1.and.mdiff.ne.12) print*,'mdiff should be for anomalies from'// &
' the annual mean or 12 for anomalies from the long term monthly mean'

tmonmean=0.
tmonval=0

il=0
do iy=1,mmax/mdiff
  do im=1,mdiff    
    il=il+1
    if(tin(il).gt.miss_val) then
      tmonmean(im)=tmonmean(im)+tin(il)
      tmonval(im)=tmonval(im)+1
    endif
  enddo
enddo

do im=1,mdiff
    if(tmonval(im) .ge. thresh) then 
      tmonmean(im)=tmonmean(im)/tmonval(im)
    else
      tmonmean(im)=miss_val
    endif
enddo

il=0
do iy=1,mmax/mdiff
  do im=1,mdiff    
    il=il+1
    if(tin(il).gt.miss_val .and. tmonmean(im) .gt. miss_val) then
      tout(il)=tin(il)-tmonmean(im)
    else
      tout(il)=miss_val
    endif
  enddo
enddo

end subroutine anomaly
subroutine makebinallwind(first,last,year,month,day,time,psvals,statnr,wmonrs,wmolats,wmolons,us,uanfeeds,ufgfeeds,vs,vanfeeds,vfgfeeds,uflags, &
               statmax,pmax,parmax,nmax,mmax,indexmax,flagmax,tmax,prefix,lsave,L_ERA40)

! externals: dates

implicit none

INTEGER, PARAMETER :: JPIM = 4
INTEGER, PARAMETER :: JPRM= 4

integer statmax,pmax,parmax,nmax,mmax,tmax,first,last,flagmax,indexmax
integer i,iunit,l,height,it,ip,ik,iret,time2,il,err,iunit2
integer tim100,day100,year100,mon100,tim101,day101,year101,mon101
integer month(nmax),year(nmax),day(nmax),time(nmax),goodindex(indexmax)

integer statnr,statnr2

integer flags(flagmax),stattype,index,indexmax2,il2,pmax2,tmax2

integer times(nmax,1,tmax)

logical lsave ! efficient storage ??
logical L_ERA40

real(kind=JPRM) us(nmax,pmax,1,tmax)
real(kind=JPRM) vs(nmax,pmax,1,tmax)
real(kind=JPRM) uanfeeds(nmax,pmax,1,tmax)
real(kind=JPRM) ufgfeeds(nmax,pmax,1,tmax)
real(kind=JPRM) vanfeeds(nmax,pmax,1,tmax)
real(kind=JPRM) vfgfeeds(nmax,pmax,1,tmax)
real(kind=JPRM) uflags(nmax,pmax,1,tmax)
!real(kind=JPRM) ttypes(nmax,pmax,1,tmax)

real(kind=JPRM) usm(nmax,pmax,1,tmax)
real(kind=JPRM) vsm(nmax,pmax,1,tmax)
real(kind=JPRM) uansm(nmax,pmax,1,tmax)
real(kind=JPRM) ufgsm(nmax,pmax,1,tmax)
real(kind=JPRM) vansm(nmax,pmax,1,tmax)
real(kind=JPRM) vfgsm(nmax,pmax,1,tmax)

real(kind=JPRM),allocatable :: us2(:,:,:,:)
real(kind=JPRM),allocatable :: uans2(:,:,:,:)
real(kind=JPRM),allocatable :: ufgs2(:,:,:,:)
real(kind=JPRM),allocatable :: vs2(:,:,:,:)
real(kind=JPRM),allocatable :: vans2(:,:,:,:)
real(kind=JPRM),allocatable :: vfgs2(:,:,:,:)
integer,allocatable :: uflags2(:,:,:,:)
integer,allocatable :: goodindex2(:)

real(kind=JPRM) lat,lon,pres

integer wmonrs(statmax)
real(kind=JPRM) wmolats(statmax),wmolons(statmax)

real(kind=JPRM) ps,psvals(pmax),measurements(6)

character*6   cstatnr
character*(*) prefix
character*249 zeile
character*1 czeit
character*10 cdatum

logical ex

REAL(KIND=JPRM) :: THILF(pmax)
REAL(KIND=JPRM) :: NEWTBCHILF(pmax),NEWRADBCHILF(pmax)
real(KIND=JPRM) :: lats(nmax,1,tmax)
real(KIND=JPRM) :: lons(nmax,1,tmax)

LOGICAL              :: LD_LBC = .FALSE.,LL_FILES_NOT_READ
INTEGER(KIND=JPIM)   :: K_IY,K_IY_OLD
INTEGER(KIND=JPIM)   :: K_IM
INTEGER(KIND=JPIM)   :: K_ID
INTEGER(KIND=JPIM)   :: K_IH
REAL(KIND=JPRM)      :: P_RMISS = -999.,psvals_pa(pmax)
INTEGER(KIND=JPIM)   :: I_MTBL,I_MCOR,I_MCORRAD,I_MSGT,IOS
INTEGER(KIND=JPIM) :: RTYPE,RCOMP


times=0

psvals_pa=psvals*100.
write(cstatnr,'(I6.6)') statnr

print *,prefix
print*,'x'//prefix//cstatnr//'x'
!print*,wmonrs
iunit=20
i=0
inquire(file=prefix//cstatnr,exist=ex)
if(.not.ex) then
  inquire(file=prefix//cstatnr//'.gz',exist=ex)
  call system('gzip -d '//prefix//cstatnr//'.gz',iret)
  if(iret .ne. 0) then 
    print*,'gzip -d -f '//prefix//cstatnr//'.gz failed'
    return
  endif
endif

open(iunit,file=prefix//cstatnr,status='old',err=20)

index=1

    l=1
do while (.true.)
10  continue 
    read(iunit,'(A249)',end=20) zeile

    if(zeile(12:14) .ne. '101' .and.  zeile(12:14) .ne.' 91' .and. zeile(12:14) .ne. '102' .and. zeile(12:14) .ne.' 92') goto 10
    if(zeile(15:20) .eq. '******') goto 10
    read(zeile(45:52),'(F8.2)') ps
    ip=1
    do while(ip.le.pmax .and. ps.ne.psvals(ip))
      ip=ip+1
    enddo
    if(ip.eq.pmax+1) goto 10

    read(zeile(1:52),'(I10,I4,3X,I6,2F8.2,I5,F8.2)') time2,stattype,statnr2,lat,lon,height,pres
   
    read(zeile(173:230),'(6F8.2,5I2)') measurements,flags

    year100=time2/1000000
    mon100=mod(time2/10000,100)
    day100=mod(time2/100,100)
    tim100=mod(time2,100)

    if(tim100 .le. 3 .or. tim100 .gt. 21) then
     it=1
     if(tim100 .gt. 21) then
       call dates(year100,mon100,day100,tim100,24,year101,mon101,day101,tim101)
       year100=year101
       mon100=mon101
       day100=day101
     endif
    else if   (tim100 .gt. 3 .and. tim100 .le. 9) then
     it=2
    else if     (tim100 .gt. 9 .and. tim100 .le. 15) then
     it=3
    else if    (tim100 .gt. 15 .and. tim100 .lt. 21) then
     it=4
    else
       goto 10
    endif
    
    if( (it .eq. 2 .or. it .eq. 4).and. tmax.eq. 2) goto 10

    if(it.eq.3 .and. tmax .eq. 2) it=2
    
    do while(year(index) .lt. year100)
      index=index+1
    enddo
    do while(month(index) .lt. mon100)
      index=index+1
    enddo
    do while(day(index) .lt. day100)
      index=index+1
    enddo

    times(index,l,it)=time2
    lats(index,l,it)=lat
    lons(index,l,it)=lon
      stop 'MAKEBINALLWIND: probably wrong date for measurements >22h - please correct me'

    us(index,ip,l,it)=-999.
    ufgfeeds(index,ip,l,it)=-999.
    uanfeeds(index,ip,l,it)=-999.
    vs(index,ip,l,it)=-999.
    vfgfeeds(index,ip,l,it)=-999.
    vanfeeds(index,ip,l,it)=-999.

    if(.not. any(flags .gt. 1)) then
 
    if(measurements(1) .ne. -999.)  us(index,ip,l,it)=measurements(1)
    if(measurements(1) .ne. -999.)  vs(index,ip,l,it)=measurements(4)

    if(measurements(2) .ne. -999.) ufgfeeds(index,ip,l,it)=-measurements(2) ! fg-obs instead of obs-fg
    if(measurements(3) .ne. -999.) uanfeeds(index,ip,l,it)=-measurements(3) ! an-obs instead of ob
    if(measurements(5) .ne. -999.) vfgfeeds(index,ip,l,it)=-measurements(5) ! fg-obs instead of obs-fg
    if(measurements(6) .ne. -999.) vanfeeds(index,ip,l,it)=-measurements(6) ! an-obs instead of ob
!    ttypes(index,ip,l,it)=rcomp*1000+rtype ! 
    endif 
    uflags(index,ip,l,it)=flags(1)*10000+flags(2)*1000+flags(3)*100+flags(4)*10+flags(5)

  i=i+1
enddo
20 print *,i,' lines read'
close(iunit)

  usm=us
  uansm=uanfeeds
  ufgsm=ufgfeeds
  vsm=vs
  vansm=vanfeeds
  vfgsm=vfgfeeds

  if(.not. L_ERA40) then
    iunit2=22
    inquire(file='../../bufrseriesall/all/'//prefix//'binsave'//cstatnr,exist=ex)
    if(ex) then
    open(iunit2,file='../../bufrseriesall/all/'//prefix//'binsave'//cstatnr,action='read',status='old',form='unformatted')      
      read(iunit2) indexmax2,il2,pmax2,tmax2
allocate(goodindex2(il2))
allocate(us2(il2,pmax2,1,tmax2),uans2(il2,pmax2,1,tmax2),ufgs2(il2,pmax2,1,tmax2))
allocate(vs2(il2,pmax2,1,tmax2),vans2(il2,pmax2,1,tmax2),vfgs2(il2,pmax2,1,tmax2))
allocate(uflags2(il2,pmax2,1,tmax2))
      read(iunit2) goodindex2
      read(iunit2) us2
      read(iunit2) uans2
      read(iunit2) ufgs2
      read(iunit2) vs2
      read(iunit2) vans2
      read(iunit2) vfgs2
      read(iunit2) uflags2
   
  close(iunit2)

  il=1
  do while(goodindex2(il) .lt. 16072 .and. il .le. il2) 
     il=il+1
  enddo

  usm(goodindex2(1:il-1),:,:,:)=us2(1:il-1,:,:,:)
  uansm(goodindex2(1:il-1),:,:,:)=uans2(1:il-1,:,:,:)
  ufgsm(goodindex2(1:il-1),:,:,:)=ufgs2(1:il-1,:,:,:)
  vsm(goodindex2(1:il-1),:,:,:)=vs2(1:il-1,:,:,:)
  vansm(goodindex2(1:il-1),:,:,:)=vans2(1:il-1,:,:,:)
  vfgsm(goodindex2(1:il-1),:,:,:)=vfgs2(1:il-1,:,:,:)

    else
     write(*,*) cstatnr//': no ERA-40 daily data found'
    endif
  endif

  call write_sonde_wind_monthly(iunit,prefix//'binumon'//cstatnr,usm,uansm,ufgsm,indexmax,mmax,pmax,tmax,month,P_RMISS,err)
  call write_sonde_wind_monthly(iunit,prefix//'binvmon'//cstatnr,vsm,vansm,vfgsm,indexmax,mmax,pmax,tmax,month,P_RMISS,err)

  il=1
  do i=1,indexmax
    if(any(us(i,:,:,:) .gt. -999.)) then
      goodindex(il)=i
      us(il,:,:,:)=us(i,:,:,:)
      vs(il,:,:,:)=vs(i,:,:,:)
      uanfeeds(il,:,:,:)=uanfeeds(i,:,:,:)
      ufgfeeds(il,:,:,:)=ufgfeeds(i,:,:,:)
      vanfeeds(il,:,:,:)=vanfeeds(i,:,:,:)
      vfgfeeds(il,:,:,:)=vfgfeeds(i,:,:,:)
      uflags(il,:,:,:)=uflags(i,:,:,:)
       il=il+1
    endif
  enddo
  il=il-1
  if(L_ERA40) then 
    open(iunit,file=prefix//'binsave'//cstatnr,form='unformatted')
  else
    open(iunit,file=prefix//'binsaveoper'//cstatnr,form='unformatted')
  endif
  write(iunit) indexmax,il,pmax,tmax
  write(iunit) goodindex(1:il)
  write(iunit) us(1:il,:,1,:)
  write(iunit) uanfeeds(1:il,:,1,:)
  write(iunit) ufgfeeds(1:il,:,1,:)
  write(iunit) vs(1:il,:,1,:)
  write(iunit) vanfeeds(1:il,:,1,:)
  write(iunit) vfgfeeds(1:il,:,1,:)
  write(iunit) uflags(1:il,:,1,:)

! write(iunit) newtbcs(1:indexmax,:,1,:),newradbcs(1:indexmax,:,1,:)
  close(iunit)

!call system('/bin/rm  '//prefix//cstatnr,iret)
!if(iret .ne. 0) print*,prefix//cstatnr//' not removed'

return
967 write(*,*) 'error opening bias correction tables'
call exit(1)
return
end subroutine makebinallwind

subroutine makebinwind(mon,psvals,values,times,statnrs,wmonrs,wmolats,wmolons,lats,lons,ts,tbcs,tanfeeds,tfgfeeds,tflags, &
               statmax,pmax,parmax,nmax,mmax,flagmax,tmax,L_ERA40,prefix,cdatum)

! externals: dates

implicit none

integer statmax,pmax,parmax,nmax,mmax,tmax,first,last,flagmax
integer mon,i,iunit,l,height,time,statnr,oldstatnr,it,ip,ik,iret
integer tim100

integer values(pmax,statmax,tmax)
integer times(nmax,pmax,statmax,tmax)
integer statnrs(statmax,tmax)

integer flags(flagmax),stattype

real(kind=JPRM) lats(nmax,statmax,tmax)
real(kind=JPRM) lons(nmax,statmax,tmax)
real(kind=JPRM) ts(nmax,pmax,statmax,tmax)
real(kind=JPRM) tbcs(nmax,pmax,statmax,tmax)
real(kind=JPRM) tanfeeds(nmax,pmax,statmax,tmax)
real(kind=JPRM) tfgfeeds(nmax,pmax,statmax,tmax)
integer*1 tflags !(flagmax,nmax,pmax,statmax,tmax)

real(kind=JPRM) lat,lon,pres,oldlat,oldlon

integer wmonrs(statmax)
real(kind=JPRM) wmolats(statmax),wmolons(statmax)

real(kind=JPRM) ps,psvals(pmax),measurements(4)

character*(*) cdatum
character*(*) prefix
character*249 zeile

logical ex,L_ERA40

! --- ARRAYS FOR STATION GROUP TABLE
INTEGER, PARAMETER :: JPRM=4
INTEGER, PARAMETER :: JPIM=4
INTEGER(KIND=JPIM)::I_NXGC      , I_NXGP      
PARAMETER (I_NXGC=15000, I_NXGP=3100)
REAL(KIND=JPRM)     , DIMENSION (I_NXGC) :: Z_PLAT,Z_PLON
REAL(KIND=JPRM)     , DIMENSION (I_NXGP) :: Z_QLAT,Z_QLON
CHARACTER(LEN= 6), DIMENSION (I_NXGC) :: CL_CSNG,CL_CSNO
CHARACTER(LEN= 5), DIMENSION (I_NXGC) :: CL_CWMO
CHARACTER(LEN= 6) :: CL_C1,CL_C4,CSTATNR
CHARACTER(LEN= 1) :: CL_C0
CHARACTER(LEN= 5) :: CL_C2
CHARACTER(LEN=19) :: CL_C3
CHARACTER(LEN=21) :: CDATSOURCE
INTEGER(KIND=JPIM)  , DIMENSION (I_NXGP) :: JSGC,JEGC,I_MREP
INTEGER(KIND=JPIM)  , DIMENSION (I_NXGC) :: I_EC,I_NC,I_JM
REAL(KIND=JPRM)    :: LAT4,LON4,Z_R1,Z_R2
INTEGER(KIND=JPIM) :: IOS,I_MSGT
INTEGER(KIND=JPIM) :: I_NGC,I_NGP,I_KM,I_KG,J,I5,I1,I2,ISY,ISM,ISD,ISH,IEY,IEM,IED,IEH,J_EC,J_NC,J_JM,I_KJ,IGOOD,IBAD,IRELOC,IUSED

!       1.2 READ STATION GROUP TABLE

  I_MSGT=21
  OPEN(UNIT=I_MSGT, FILE='stgroup.t', &
   & IOSTAT=IOS,ERR=967,STATUS='OLD')
  I_NGC = 0  ! NUMBER OF RECORD
  I_NGP = 0  ! NUMBER OF GROUP
  350   CONTINUE
  READ(I_MSGT, &
   & '(I5,I4,A1,A6,F7.2,F8.2,A8,2(1X,I4,3I2),I7,23X,A6)', &
   & END=190) &
   & I1,I2,CL_C0,CL_C1,Z_R1,Z_R2,CL_C2, &
   & ISY,ISM,ISD,ISH,IEY,IEM,IED,IEH,I5,CL_C4  

  I_NGC = I_NGC+1
  Z_PLAT(I_NGC) = Z_R1
  Z_PLON(I_NGC) = Z_R2
  CL_CSNO(I_NGC) = CL_C1
!  READ(CDATSOURCE,'(3I1)') J_EC,J_NC,J_JM
!  I_EC(I_NGC)=J_EC
!  I_NC(I_NGC)=J_NC
!  I_JM(I_NGC)=J_JM
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
    CL_C4=CL_CSNG(I_NGC)
    CL_CWMO(I_NGP) = CL_C4(1:5)
    I_MREP(I_NGP) = I_NGC
    Z_QLAT(I_NGP) = Z_R1
    Z_QLON(I_NGP) = Z_R2
  ENDIF
  GOTO 350
  190   CONTINUE

!write(cdatum,'(I4,I2.2)') 1957+mon/12,mod(mon,12)+1

print *,prefix
!print*,wmonrs
iunit=20
inquire(file=prefix//cdatum(1:6),exist=ex)
if(.not.ex) then
  inquire(file=prefix//cdatum(1:6)//'.gz',exist=ex)
!  if(.not.ex) then
!    call system('ecp ec:feedbacknew/'//prefix//cdatum(1:6)//'.gz .',iret)
!    if(iret .ne. 0) then 
!      print*,prefix//cdatum//'.gz not found'
!      return
!    endif
!  endif
  call system('gzip -d '//prefix//cdatum(1:6)//'.gz',iret)
  if(iret .ne. 0) then 
    print*,'gzip -d -f '//prefix//cdatum(1:6)//'.gz failed'
    return
  endif
endif


open(iunit,file=prefix//cdatum(1:6),status='old',err=20)
IGOOD=0
IRELOC=0
IBAD=0
IUSED=0
i=0
do while (.true.)
10  continue 
    read(iunit,'(A93)',end=20) zeile

    i=i+1
    if(zeile(12:14) .ne. '101') goto 10
    if(zeile(15:20) .eq. '******') goto 10
    read(zeile(42:49),'(F8.2)') ps
    ip=1
    do while(ip.lt.15 .and. ps.ne.psvals(ip))
      ip=ip+1
    enddo
    if(ip.eq.15) goto 10
!    if(index(zeile,'65.00') .ne. 0 .or.  index(zeile,'67.00') .ne. 0) goto 10
!    print*,ps
!    if(zeile(37:41) .eq. '*****') then 
!       zeile(37:41)=' -999'
!       print *,zeile
!    endif

    read(zeile(1:49),'(I10,I4,I6,2F8.2,I5,F8.2)') time,stattype,statnr,lat,lon,height,pres

    if(statnr .ne. oldstatnr .or. lat .ne. oldlat .or. lon .ne. oldlon) then 
      lat4=lat
      lon4=lon
      WRITE(cstatnr,'(I6.6)') statnr
      call identify_station(cstatnr,CL_CSNO,lat4,lon4,Z_PLAT,Z_PLON,JSGC,JEGC,I_NGP,I_NGC,I_KG,I_KM,I_KJ)
      oldstatnr=statnr
      oldlat=lat
      oldlon=lon
    endif

    IF(I_KG*I_KM .NE. 0) THEN
     if(cstatnr(1:6) .ne. CL_CWMO(I_KG)) THEN
!      WRITE(*,*) cstatnr(1:6),CL_CSNO(I_KJ),CL_CWMO(I_KG)
      READ(CL_CWMO(I_KG),'(I5)') statnr
      lat=Z_QLAT(I_KG)
      lon=Z_QLON(I_KG)
!      WRITE(*,*) 'statnr changed to ',statnr,lat,lon
      IRELOC=IRELOC+1
     ELSE
      IGOOD=IGOOD+1
     ENDIF
    ELSE
!      WRITE(*,*) cstatnr(1:6),' ',CL_CSNO(I_KJ),' ',CL_CWMO(I_KG)
      IBAD=IBAD+1
    ENDIF

    l=1
    if(statmax .gt. 3000) then ! all stations
      do while(statnr.ne.wmonrs(l) .and. l .le. statmax)
        l=l+1
      enddo
      if(l .gt. statmax) goto 10
      if(wmonrs(l).eq.0) goto 10
    else ! only CARDS stations
      do while(wmonrs(l).ne.0 .and. statnr.ne.wmonrs(l))
        l=l+1
      enddo
      if(wmonrs(l).eq.0) goto 10
    endif
    
!    if(abs(lat-wmolats(l)) .gt. 2. .or. abs(lon-wmolons(l)) .gt. 2.) then
!      print*,wmonrs(l),lat,lon,wmolats(l),wmolons(l),' discarded'
!      goto 10
!    endif

!    print*,i,statnr,wmonrs(i)

! old ASCII feedback file
!    nz=zeile(46:53)//zeile(46+6*8:53+6*8)//zeile(46+12*8:53+12*8)//&
!       zeile(190:191)//zeile(202:203)//zeile(214:215)//zeile(226:227)//zeile(238:239)
!    read(nz,'(3F8.2,5I2)') measurements(1:3),flags(1:5)

! new ASCII feedback file
    read(zeile(50:93),'(4F8.2,6I2)') measurements,flags

    if(any(flags(1:4) .eq.3) .or. any(measurements .eq. -999.)) then  ! all observations rejected by QC or with no feedback information are not used - blacklisted values are used
      goto 10
    endif

    
    tim100=mod(time,100)

    if(tim100 .le. 3 .or. tim100 .gt. 21) then
     it=1
    else if   (tim100 .gt. 3 .and. tim100 .le. 9) then
     it=2
    else if     (tim100 .gt. 9 .and. tim100 .le. 15) then
     it=3
    else if    (tim100 .gt. 15 .and. tim100 .lt. 21) then
     it=4
    else
       goto 10
    endif
    
    if( (it .eq. 2 .or. it .eq. 4).and. tmax.eq. 2) goto 10

    if(it.eq.3 .and. tmax .eq. 2) it=2
    

      do ik=1,values(ip,l,it)
        if(times(ik,ip,l,it).eq.time .or. times(ik,ip,l,it).eq.time-1 .or.&
           times(ik,ip,l,it).eq.time+1 .or. times(ik,ip,l,it).eq.time-2 .or.&
           times(ik,ip,l,it).eq.time+2 .or. times(ik,ip,l,it).eq.time-77 .or. &
           times(ik,ip,l,it).eq.time+77.or. times(ik,ip,l,it).eq.time-78 .or. &
           times(ik,ip,l,it).eq.time+78) then
!          print*,'duplicate entry',statnr,times(ik,ip,l,it),time,measurements(1)
          goto 10
        endif
      enddo

      if(values(ip,l,it) .lt. 31) then
!        print*,ip,l,it, values(ip,l,it)
        values(ip,l,it)=values(ip,l,it)+1
!        print*,ip,l,it, values(ip,l,it)
      else
        print*,pres,statnr,it,time,'too many values'
        goto 10
      endif
    
!      print*,zeile,ip,l,it,values(ip,l,it)
      times(values(ip,l,it),ip,l,it)=time

      lats(values(ip,l,it),l,it)=lat
      lons(values(ip,l,it),l,it)=lon

      ts(values(ip,l,it),ip,l,it)=measurements(1)
      stop 'MAKEBINWIND: probably wrong date for measurements >22h - please correct me'

! bias correction is observation_presented_to_analysis-original_observation
      tbcs(values(ip,l,it),ip,l,it)=measurements(2)-measurements(1)
        tfgfeeds(values(ip,l,it),ip,l,it)=-measurements(3) ! fg-obs instead of obs-fg
      tanfeeds(values(ip,l,it),ip,l,it)=-measurements(4) ! an-obs instead of ob

!    tflags(:,values(ip,l,it),ip,l,it)=flags

!    read(zeile(46:249),'(18F8.2,30I2)') measurements,flags

!    if(i .gt. 100) stop
enddo
20 print *,i,' lines read ',IGOOD,' obs in place',IRELOC,' obs relocated',IBAD,' obs not identified'
close(iunit)
call system('/bin/rm  '//prefix//cdatum(1:6),iret)
if(iret .ne. 0) print*,prefix//cdatum(1:6)//' not removed'

return
967 STOP 'COULD NOT OPEN STGROUPS.T'

end subroutine makebinwind

subroutine write_sonde_wind_monthly(iunit,filename,us,uans,ufgs,nmax,mmax,pmax,parmax,month,miss_val,err)

implicit none

integer i,iunit,err,l,iy,index,imon
integer,intent(in) :: nmax,mmax,pmax,parmax,month(20000)
character*(*) filename
real(kind=JPRM),intent(in) :: us(nmax,pmax,parmax),uans(nmax,pmax,parmax),ufgs(nmax,pmax,parmax),miss_val
real(kind=JPRM) :: itm(nmax,pmax,parmax)
real(kind=JPRM) :: tmmon(mmax,pmax,parmax),rasocorrmon(mmax,pmax,parmax),eracorrmon(mmax,pmax,parmax)
integer         :: goodmon(mmax,pmax,parmax),imod,tgm,rgm,egm,ip,ipar

  iunit=iunit

  open(iunit,file=filename,form='unformatted',err=120)
  err=0

    write(iunit) mmax,pmax,parmax

do i=1,3
  if(i .eq. 1) itm=us
  if(i .eq. 2) itm=uans
  if(i .eq. 3) itm=ufgs
 
  do ipar=1,parmax
   do ip=1,pmax
    index=1
    do imon=1,mmax
      iy=1957+(imon-1)/12
      imod=mod(imon-1,12)+1
      tmmon(imon,ip,ipar)=0.
      tgm=0
      do while(month(index) .eq. imod .and. index .lt. nmax)
        if(itm(index,ip,ipar) .ne. miss_val) then
          tmmon(imon,ip,ipar)=tmmon(imon,ip,ipar)+itm(index,ip,ipar)
          tgm=tgm+1
        endif
        index=index+1
!        write(*,*) index,imon,iy, mod(imon-1,12)+1,sum(goodmon(imon,ip,ipar))
      enddo
        
      goodmon(imon,ip,ipar)=tgm
      if(tgm .gt. 0) then
         tmmon(imon,ip,ipar)=tmmon(imon,ip,ipar)/tgm
      else
         tmmon(imon,ip,ipar)=miss_val
      endif   
    enddo
   enddo
  enddo
    write(iunit) tmmon
enddo

    write(iunit) goodmon
    close(iunit)

  return

120 print*,'could not open file ',filename
  err=1

return
end subroutine write_sonde_wind_monthly

SUBROUTINE BIASCOR_ERA40_changed (LD_LBC,K_NST,P_SPP,P_ST,P_STNBC, &
 & P_SOLBC,KSON,P_ANGLE,CD_CIDENT,K_IY,K_IM,K_ID,K_IH,P_RLAT,P_RLON,K_IRSTYP,  &
 & K_IMISS,P_RMISS,LDSOLAR,LDHOMOGEN,Z_RASOCORRS,rcpara)  

type(rasocor_namelist),intent(in) :: rcpara

!leoUSE PARKIND1  ,ONLY : JPIM     ,JPRM
!leoUSE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
!leoUSE YOMOML
                   
!**** *BIASCOR_ERA40* 

!       PURPOSE.
!      ------------%

!       BIAS CORRECTION OF RADIOSONDE TEMPERATURES FOR 

!          EEEEE  RRRR      A         4    000
!          E      R   R    A A       44   0   0
!          EEEE   RRRR    A   A     4 4   0   0
!          E      R R     AAAAA    44444  0   0
!          EEEEE  R   R   A   A       4    000

!       INTERFACE.
!      ------------

!          CALL BIASCOR_ERA40 (LBC,NST,SPP,ST,STNBC,
!     X                  CIDENT,IM,ID,IH,RLAT,RLON,IRSTYP,
!     X                  IMISS,RMISS)

!        INPUT
!         NST       NUMBER OF LEVELS (4 BYTE INTEGER)
!         SPP       ARRAY WITH PRESSURE VALUES (Pa) (8 BYTE REAL)
!         ST        ARRAY WITH T VALUES (8 BYTE REAL)
!         CIDENT    STATION IDENTIFIER  (CHARACTER)
!         IM        MONTH (4 BYTE INTEGER)
!         ID        DAY   (4 BYTE INTEGER)
!         IH        HOUR  (4 BYTE INTEGER)
!         RLAT      LATITUDE (8 BYTE REAL)
!         RLON      LONGITUDE (8 BYTE REAL)
!         IRSTYP    RADIOSONDE TYPE (BUFR CODE TABLE 002011) 
!                   (Not need for ERA)
!         IMISS     MISSING VALUE FOR INTEGERS (4BYTE INTEGER)
!         RMISS     MISSING VALUE FOR REALS (8 BYTE REAL) 

!        OUTPUT
!         LBC       LOGICAL TO STATE IF BIAS COR. WAS SUCCESSFUL
!         STNBC     ARRAY WITH BIAS CORRECTED T VALUES

!       METHOD.
!      ---------

!        READ BIAS CORRECTION TABLES:
!            STGROUP.T : STATION GROUP TABLE
!            CORRECT.T : ARRAY OF CORRECTIONS DEPENDING ON
!                        SOLAR ELEVATION AND PRESSURE LEVEL
!            COUNTRY.T : DEFINITION OF CATEGORIES

!        1ST A STATION GROUP 0F THE DATA IS DETECTED.
!        2ND A CATEGORY INCLUDING THE STATION GROUP IS DETECTED.
!        3RD IF THE CATEGORY IS ONE FOR CORRECTION, APPLY CORRECTION.

!        FIRST 32 CHARACTERS ARE COMPARED TO DETECT CATEGORY.
!        'CATG' FROM COUNTRY.T AND 'YMNAMBC' FROM CORRECT.T

!       EXTERNALS.
!      ------------

!        DIURNAL          CALCULATE SOLAR ELEVATION
!        PNTERP           INTERPOLATE TO PRESSURE LEVEL

!       REFERENCE.
!      ------------

!        RADIOSONDE BIAS CORRECTION
!        OD MEMORANDUM BY B. STRAUSS  22.12.92

!       AUTHOR.
!      ---------

!        B. NORRIS       JULY  1991   APPLY CORRECTIONS TO AOF

!       MODIFICATIONS.
!      ----------------

!        M. DRAGOSAVAC   APRIL 1993   IMPLEMENT IN PRE-PROCESSING
!        B. NORRIS       JULY  1998   INTERFACE TO DATA ASSIMILATION
!        K. ONOGI        MARCH 2000   MODIFIED FOR ERA-40
!        K. ONOGI      OCTOBER 2000   MAKE NUMBER OF B.COR. CATEGORIES 8 TO 4
!                                     AS SAME AS THE CORRECTION TABLE
!        S. SAARINEN  NOVEMBER 2000   CLEANING UP PLUS IMPLICIT NONE
!        M. Hamrud      01-Oct-2003   CY28 Cleaning
!      L. Haimberger  NOVEMBER 2004   Allow use of tables from  radiosonde homogenization 
!                                    (see ERA-40 Project report series 23)
!        S. SAARINEN       MAY 2005   Removed standalone SAVE-stmt to allow Dr.Hook
!-------------------------------------------------------------------

!IMPLICIT NONE

!     Table values are expected to be saved
!     after they are read

! INTEGER,PARAMETER :: JPIM=4
! INTEGER,PARAMETER :: JPRM=4
! INTEGER,PARAMETER :: JPRM=4

LOGICAL           ,INTENT(OUT)   :: LD_LBC 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_NST 
REAL(KIND=JPRM)   ,INTENT(IN)    :: P_SPP(K_NST) 
REAL(KIND=JPRM)   ,INTENT(IN)    :: P_ST(K_NST) 
REAL(KIND=JPRM)   ,INTENT(OUT)   :: P_STNBC(K_NST) 
REAL(KIND=JPRM)   ,INTENT(OUT)   :: P_SOLBC(K_NST) 
CHARACTER(LEN=5)  ,INTENT(IN)    :: CD_CIDENT 
INTEGER(KIND=JPIM),INTENT(OUT)    :: KSON 
REAL(KIND=JPRM),INTENT(OUT)    :: P_ANGLE 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_IM 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_ID 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_IH 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_IY 
REAL(KIND=JPRM)   ,INTENT(IN)    :: P_RLAT 
REAL(KIND=JPRM)   ,INTENT(IN)    :: P_RLON 
INTEGER(KIND=JPIM)               :: K_IRSTYP ! Argument NOT used
INTEGER(KIND=JPIM)               :: K_IMISS ! Argument NOT used
REAL(KIND=JPRM)   ,INTENT(IN)    :: P_RMISS 
LOGICAL   ,INTENT(IN)            :: LDSOLAR,LDHOMOGEN

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

REAL(KIND=JPRM)     , DIMENSION (I_NCORBC,I_NLEVBC,I_NSONBC) ,SAVE:: Z_RCORBC
INTEGER(KIND=JPIM)  , DIMENSION (I_NLEVBC,I_NSONBC) ,SAVE:: ILEVBC
REAL(KIND=JPRM)     , DIMENSION (I_NCORTB) ,SAVE:: Z_WKCOR
INTEGER(KIND=JPIM) ,SAVE:: IPLAT,IPLON,IRLAT,IRLON
CHARACTER(LEN= 5) ,SAVE:: CL_YDUMMYBC,CL_VERSION,CL_ISTAT
CHARACTER(LEN=32) , DIMENSION (I_NSONBC) ,SAVE:: CL_YSNAMBC
CHARACTER(LEN=58) ,SAVE:: CL_CDATE
CHARACTER(LEN= 2) ,SAVE:: CL_CPMAX
CHARACTER(LEN= 3) ,SAVE:: CL_F
CHARACTER(LEN= 8)      :: CLPUINT
CHARACTER(LEN=13)      :: CLR='BIASCOR_ERA40'
CHARACTER(LEN=4) :: CLYEAR

LOGICAL ,SAVE:: LL_FILES_NOT_READ,LL_DEBUG=.false.
INTEGER(KIND=JPIM) ,SAVE:: I_NSOBC
      
! --- ARRAYS FOR STATION GROUP TABLE
REAL(KIND=JPRM)     , DIMENSION (I_NXGC) ,SAVE:: Z_PLAT,Z_PLON
CHARACTER(LEN= 6), DIMENSION (I_NXGC) ,SAVE:: CL_CSNG,CL_CSNG2,CL_CSNO
CHARACTER(LEN= 6) ,SAVE:: CL_C1,CL_C4
CHARACTER(LEN= 1) ,SAVE:: CL_C0
INTEGER(KIND=JPIM)  , DIMENSION (I_NXGP) ,SAVE:: JSGC,JEGC,I_MREP
REAL(KIND=JPRM)     , DIMENSION (I_NXGP) ,SAVE:: Z_QLAT,Z_QLON
INTEGER(KIND=JPIM)  , DIMENSION (I_NXCT) ,SAVE:: I_NTGT
INTEGER(KIND=JPIM)  , DIMENSION (I_NXGP,I_NXCT) ,SAVE:: I_MTGT

! --- ARRAYS FOR INDEX TABLE
INTEGER(KIND=JPIM)  , DIMENSION (I_NXDG,I_NXCT) ,SAVE:: IDSTA,IDEND, &
 & I_LATS,I_LATE,I_LONS,I_LONE  
INTEGER(KIND=JPIM)  , DIMENSION (I_NXDG) ,SAVE:: I_NCMB,I_NCDU
INTEGER(KIND=JPIM)  , DIMENSION (I_NXCT) ,SAVE:: I_NCCC,I_MCAT
CHARACTER(LEN=64), DIMENSION (I_NXDG,I_NXCT) ,SAVE:: CL_CNTRY
CHARACTER(LEN=64), DIMENSION (I_NXCT) ,SAVE:: CL_CDG
CHARACTER(LEN=64) ,SAVE::  CL_C64
CHARACTER(LEN=18) ,SAVE::  CLTLN
CHARACTER(LEN=32), DIMENSION (I_NXCT) ,SAVE:: CL_CATG

LOGICAL, PARAMETER :: LGRP = .TRUE.

! --- LEO VARIABLES
INTEGER(KIND=JPIM) :: ierr,icsn
LOGICAL :: CLEO_FOUND = .FALSE.
CHARACTER(LEN=50) :: CLEONAME = 'table4'
CHARACTER(LEN=120) :: CL_ZEILE
INTEGER(KIND=JPIM), PARAMETER :: IPMAX    = 16
INTEGER(KIND=JPIM), PARAMETER :: IPARMAX   = 2
INTEGER(KIND=JPIM), PARAMETER :: IESTATMAX = 3070
INTEGER(KIND=JPIM), PARAMETER :: NMAX = 18993
CHARACTER(LEN=1),SAVE  :: CL_CCORR(IESTATMAX)
REAL(KIND=JPRM) :: Z_RASOCORRS(NMAX,IPMAX,IPARMAX)
INTEGER(KIND=JPIM),SAVE :: IEWMONRS(IESTATMAX)
INTEGER(KIND=JPIM),SAVE :: ISTAT,IB
REAL(KIND=JPRM),SAVE    :: ZLEO_PLEV(IPMAX),ZIPLEVS(IPMAX),ZLEO_BIAS(IPMAX)
real(kind=JPRM),SAVE    :: ZWMOLATS(IESTATMAX),ZWMOLONS(IESTATMAX)
real(kind=JPRM)    :: ZHILF(IPMAX,IPARMAX),ZHILFOLD(IPMAX,IPARMAX)
integer(KIND=JPIM),SAVE :: iyear,imonth,itime,iday,idatum,idatumold,i_datum,idum,ip,IIPMAX,istatmax

      
! --- MISCELLANEOUS DECLARATIONS
INTEGER(KIND=JPIM) ,SAVE:: IOS, I, IR,JDP,IDS,IDE,I_LAS,I_LAE,I_LOS,I_LOE
INTEGER(KIND=JPIM) ,SAVE:: I_M, JDG, I_N, I_NGC, I_NGP, I1, I2,I6
REAL(KIND=JPRM) ,SAVE:: Z_R1,Z_R2
INTEGER(KIND=JPIM) ,SAVE:: ISY,ISM,ISD,ISH,IEY,IEM,IED,IEH,I5
INTEGER(KIND=JPIM) ,SAVE:: I_NCNT, J, I_KG, I_KM, IXT, ISON, I_K, ISNG, ICOL
REAL(KIND=JPRM) ,SAVE:: Z_ANGLEBC, PP, Z_POS, Z_RADD
INTEGER(KIND=JPIM) ,SAVE:: ILOW, IHGH,ianf,iend,ipar
REAL(KIND=JPRM) :: ZHOOK_HANDLE
REAL(KIND=JPRM) :: Z_LEO(K_NST)

LOGICAL EX

!leo#include "abor1.intfb.h"
!leo#include "diurnal.intfb.h"
!leo#include "pnterp.intfb.h"

!  -------------------------------------------------------------------

!     1.  OPEN AND READ TABLE FILES
!     -------------------------------

!leo IF (LHOOK) CALL DR_HOOK(CLR,0,ZHOOK_HANDLE)

! CALL OML_SET_LOCK()
100 CONTINUE

WRITE(CL_F,'(I3.3)') K_NST
LD_LBC=.FALSE.
IF(K_IMISS == 1) THEN
  LL_FILES_NOT_READ=.true.
  K_IMISS=0
  I_KM=0
  I_KG=0
ENDIF
if(ll_files_not_read) then
!         OPEN(UNIT=MCOR, FILE='/era/data/erk/comtable/correct.t', &

    WRITE(CLYEAR,'(I4)') K_IY
    INQUIRE(FILE='/home/srvx9/rs_common/tables/ulfV1_1/T_correct'//CLYEAR//'010100',EXIST=EX)
    if(ex) then 
    OPEN(UNIT=I_MCOR, FILE='/home/srvx9/rs_common/tables/ulfV1_1/T_correct'//CLYEAR//'010100', &
   & IOSTAT=IOS,ERR=966,STATUS='OLD')  

!       1.3  READ CORRECTION TABLE

  I_NSOBC=0
  READ (I_MCOR,'(A)') CL_CDATE

  IF(LDSOLAR .and. LDHOMOGEN) THEN

!  write(*,*) 'reading biascor.t'
! If L. Haimbergers correction table is used, it has to be used with a special version of
! solar angle bias correction tables. See ERA-40 project report series Nr. 23 by L. Haimberger
    IF(CL_CDATE(52:52) /= '1') THEN
      CALL ABOR1(CLR//': wrong table2 for solar angle correction + homogenization')
    ENDIF
  ELSE
    IF(LDSOLAR .and. CL_CDATE(52:52) == '1') THEN
      CALL ABOR1(CLR//': wrong table2 for solar correction only')
    ENDIF
  ENDIF

  READ (I_MCOR,*)

  DO WHILE(.TRUE.)
    IF (I_NSOBC < I_NSONBC) THEN
      I_NSOBC=I_NSOBC+1
      READ (I_MCOR,'(A)',END=138) CL_YSNAMBC(I_NSOBC)
!           WRITE(*,'(2A)') YSNAMBC(NSOBC)

      DO J=1,I_NLEVBC
        READ (I_MCOR,'(I5,8F8.2)',END=138) &
       & ILEVBC(J,I_NSOBC),(Z_WKCOR(I_K),I_K=1,I_NCORTB)  

!         WRITE(*,'(I5,8F8.2)') &
!               & ILEVBC(J,NSOBC),(WKCOR(K),K=1,NCORTB)

!     CORRECTION TABLE
!     ~-7.5 : -7.5~7.5 : 7.5~22.5 : 22.5~

        DO I=1,I_NCORBC
          Z_RCORBC(I,J,I_NSOBC) = Z_WKCOR(I)
        ENDDO

      ENDDO
      READ (I_MCOR,'(A)',END=138) CL_YDUMMYBC
    ELSE
      WRITE (0,'(1X,A,I5)') &
     & CLR//': DIMENSION NSON TOO SMALL :', &
     & I_NSONBC  
      CALL ABOR1(CLR//': DIMENSION NSON TOO SMALL')
      CALL ABORT
    ENDIF
  ENDDO
  138   CONTINUE
  I_NSOBC = I_NSOBC-1

!         WRITE (0,*) &
!         & CLR//': NUMBER OF CATEGORIES IN CORRECTION TABLE ', &
!         &  NSOBC


        !!WRITE(6,*)'SUCCESSFULLY READ TABLE 1-3'
  CLOSE (I_MCOR)
  else
    Z_RCORBC=0.
    CL_CDATE =''
  endif

  LL_FILES_NOT_READ = .FALSE.

IF (CL_ISTAT .ne. CD_CIDENT ) THEN

!        WRITE(*,*)
!        WRITE(*,*) '------------------------------'
!        WRITE(*,*) '  RADIOSONDE BIAS CORRECTION  '
!        WRITE(*,*) '------------------------------'
!        WRITE(*,*)


!         OPEN(UNIT=MTBL, FILE='/era/data/erk/comtable/country.t', &
    OPEN(UNIT=I_MTBL, FILE='/home/srvx9/rs_common/tables/country.t', &
   & IOSTAT=IOS,ERR=965,STATUS='OLD')  


!         OPEN(UNIT=MSGT, FILE='/era/data/erk/comtable/stgroup.t', &
    OPEN(UNIT=I_MSGT, FILE='/home/srvx9/rs_common/tables/stgroup.t', &
   & IOSTAT=IOS,ERR=967,STATUS='OLD')  



!       1.1 READ CATEGORY DEFINITION TABLE

  I=1
  IR=0
  DO WHILE(I <= 2000 .and. IR /=10)
    READ(I_MTBL,'(I4)') IR
    I=I+1
  ENDDO
  I_NCNT = 0
  DO I=1,I_NXDG
    I_NCMB(I) = 0
  ENDDO
  NXCT: DO I=1,I_NXCT
    READ(I_MTBL,'(I4,I3,2I6,1X,2I3,2I4,1X,A64)') &
     & IR,JDP,IDS,IDE,I_LAS,I_LAE,I_LOS,I_LOE,CL_C64  
    IF (IR == 999) EXIT NXCT
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
  ENDDO NXCT

!       1.1.1 PUT LAT LON LIMIT ON CNTRY

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
!               CNTRY(I,M) = CNTRY(I,M)(1:20)//CLTLN
        CL_CNTRY(I,I_M) = CLTLN//CL_CNTRY(I,I_M)(1:20)

      ENDIF
    ENDDO
  ENDDO

!       1.1.2 GATHER SOME GROUPS INTO ONE GROUP

  DO I=1,I_NXCT
    CL_CDG(I) = ' '
  ENDDO
  IF (LGRP) THEN
    IR=0
    I=1
    DO WHILE(I <=2000 .and. IR /=20)
      READ(I_MTBL,'(I4)') IR
      I=I+1
    ENDDO
NXDG:    DO I=1,I_NXDG
      READ(I_MTBL,'(I4,I3,1X,A64)') IR,JDG,CL_C64
      IF (IR == 999) EXIT NXDG
      IF (IR == 1) THEN
        CL_CDG(I_NCDU(JDG)) = CL_C64
      ENDIF
    ENDDO NXDG
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

!       1.1.3 MONITOR THE CATEGORIES
!       DO N=1,NCNT
!         DO I=1,NXDG
!         IF (IDSTA(I,N)/=0) &
!       & WRITE(*,'(I6.6,2I7,4I5,1X,A32,2X,A32)') &
!             & I,IDSTA(I,N),IDEND(I,N),LATS(I,N),LATE(I,N), &
!             & LONS(I,N),LONE(I,N),CNTRY(I,N)(1:32),CATG(N)
!         ENDDO
!       ENDDO

!       1.2 READ STATION GROUP TABLE
  I_NGC = 0  ! NUMBER OF RECORD
  I_NGP = 0  ! NUMBER OF GROUP
  DO WHILE(.TRUE.)
  READ(I_MSGT, &
   & '(I6.6,I4,A1,A6,F7.2,F8.2,I8,2(1X,I4,3I2),I7,2X,A6)', &
   & END=190) &
   & I1,I2,CL_C0,CL_C1,Z_R1,Z_R2,I6, &
   & ISY,ISM,ISD,ISH,IEY,IEM,IED,IEH,I5,CL_C4  
!           WRITE(*, &
!           & '(I6.6,I4,A1,A6,F7.2,F8.2,A5,2(1X,4I2),I7,2X,A19,2X,A6)') &
!           &  I1,I2,C0,C1,R1,R2,C2, &
!           &  ISY,ISM,ISD,ISH,IEY,IEM,IED,IEH,I5,C3,C4
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
  CL_CSNG2(I_NGC) = CL_C1
  IF (CL_C0 == '#') THEN
    I_MREP(I_NGP) = I_NGC
    Z_QLAT(I_NGP) = Z_R1
    Z_QLON(I_NGP) = Z_R2
  ENDIF
  ENDDO
  190   CONTINUE



  CLOSE (I_MTBL)
  CLOSE (I_MSGT)


  IF(LDHOMOGEN) THEN

!         OPEN(UNIT=MSGT, FILE='biascor.t', &
  OPEN(UNIT=I_MLEO, FILE='leobiascor.t', &
   & IOSTAT=IOS,ERR=968,STATUS='OLD')  

! -------------------------------------
! 
!     Read L. Haimbergers's homogenizing temp corrections
!     http://www.univie.ac.at/theoret-met/research/RAOBCORE/     
! -------------------------------------
   
 ZLEO_PLEV=(/10.,20.,30.,50.,70.,100.,150.,200.,250.,300.,400.,500.,700.,850.,925.,1000./)
 do ip=1,IPMAX
    ilevbc(ip,1)=ZLEO_PLEV(IPMAX-ip+1)
 enddo

      ZLEO_BIAS=0.

! read and print comment line(s) in header, if available
  CL_ZEILE='#'
  do while(CL_ZEILE(1:1) == '#')
    read(I_MLEO,'(A120)') CL_ZEILE
    if(CL_ZEILE(1:1) == '#') write(*,*) CL_ZEILE
  enddo

  read(CL_ZEILE,'(A5,I5,I3)')  CL_VERSION,istatmax,IIPMAX
  write(CL_CPMAX,'(I2.2)') IIPMAX
  read(CL_ZEILE,'(A5,I5,I3,'//CL_CPMAX//'F6.0)') CL_VERSION,istatmax,IIPMAX,ZIPLEVS

  if(IPMAX /= IIPMAX .or. any(ZLEO_PLEV /= ZIPLEVS) .or. istatmax /= IESTATMAX) then
     write(6,*) CLR//': istatmax,IESTATMAX',istatmax,IESTATMAX
     write(6,*) CLR//': IPMAX,IIPMAX',IPMAX,IIPMAX
     write(6,*) CLR//': plevs',ZLEO_PLEV
     write(6,*) CLR//': ZIPLEVS',ZIPLEVS
     CALL ABOR1(CLR//': WRONG DIMENSIONS OF BIASCOR.T')
  endif
  if(LDSOLAR .and. CL_VERSION /= CL_CDATE(54:58) .and. CL_CDATE /='' ) then
    CALL ABOR1(CLR//': table2 and table4 version mismatch:'//CL_VERSION//':'//CL_CDATE(54:58)//':')
  endif

  I_DATUM=K_IY*1000000+K_IM*10000+K_ID*100

  z_rasocorrs=0.
  do istat=1,IESTATMAX
    read(I_MLEO,'(I5,I6,2F9.2,1X,A1,I3)',ERR=968) idum,IEWMONRS(istat),ZWMOLATS(istat),ZWMOLONS(istat),CL_CCORR(istat),ib

    if(CL_CCORR(istat) == "Y") then
      ZHILF=0.
      idatum=1957010100
      do i=1,ib
        ZHILFOLD=ZHILF
        idatumold=idatum
        read(I_MLEO,'(I4,3I2.2,'//CL_CPMAX//'F6.2)',ERR=968) iyear,imonth,iday,itime,ZHILF(:,1)
        read(I_MLEO,'(I4,3I2.2,'//CL_CPMAX//'F6.2)',ERR=968) iyear,imonth,iday,itime,ZHILF(:,2)
        idatum=iyear*1000000+imonth*10000+iday*100
        ianf=toindex(idatumold/100,rcpara)
        iend=toindex(idatum/100,rcpara)
        write(CL_ISTAT,'(I6.6)') iewmonrs(istat)
        IF(CD_CIDENT .eq. CL_ISTAT) then
          do ip=1,rcpara%PMAX
            do ipar=1,rcpara%PARMAX
            Z_RASOCORRS(ianf:iend,ip,ipar)=ZHILFOLD(IPMAX-ip+1,ipar)
            enddo
          enddo
        ENDIF
      enddo
!     write(*,*) 'read ',IEWMONRS(istat),ZWMOLATS(istat),ZWMOLONS(istat),CL_CCORR(istat),ib,Z_RASOCORRS(5,:,istat)
    endif

  enddo
  
  CL_ISTAT=CD_CIDENT

  CLOSE (I_MLEO)

  ENDIF ! LDHOMOGEN

ENDIF

ENDIF

!-----------------------------------------------------------------------------


!     1.4  INITIALLY SET CORRECTED = ORIGINAL

DO I=1,K_NST
  P_STNBC(I)=P_ST(I)
  P_SOLBC(I)=P_ST(I)
ENDDO

!  ------------------------------------------------------------------ 

!  Solar angle bias correction
!
!     2.   DETERMINE TABLES AND TARGET

!     2.1  DETECT THE STATION GROUP

IF(I_KG == 0 .and. I_KM == 0) THEN
I_KG = 0 ; I_KM = 0
NGP: DO I=1,I_NGP
  DO J=JSGC(I),JEGC(I)
    IPLAT = NINT(Z_PLAT(J)*100)
    IPLON = NINT(Z_PLON(J)*100)
    IRLAT = NINT(P_RLAT*100)
    IRLON = NINT(P_RLON*100)
    IF (IRLAT == IPLAT.AND.IRLON == IPLON) THEN
      IF (CD_CIDENT(1:5) == CL_CSNO(J)(1:5)) THEN
        I_KG = I               ! STATION GROUP
        I_KM = J-JSGC(I)+1     ! MEMBER 
        EXIT NGP
      ENDIF
    ENDIF
  ENDDO
ENDDO NGP

IF (I_KG*I_KM == 0) THEN
  WRITE(*,'(A,A6,F7.2,F8.2,I6,I3,2A)') &
   & CLR//': K_ID,LAT,LON,GRP,MEM= ',CD_CIDENT,P_RLAT,P_RLON,I_KG,I_KM, &
   & ' | THE SPECIFIED DATA WAS NOT FOUND ', &
   & 'IN THE STATION GROUP DEFINITION TABLE.'  
      ! THE SPECIFIED DATA WAS NOT FOUND IN THE STATION GROUP DEFINITION TABLE.
! CALL OML_UNSET_LOCK()
!leo  IF (LHOOK) CALL DR_HOOK('BIASCOR_ERA40',1,ZHOOK_HANDLE)
!  RETURN

ENDIF

ENDIF

IF(LDSOLAR .AND. I_KG*I_KM > 0 ) THEN

!     2.2  PICK UP THE TARGET CATEGORY

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

!     DO M=1,NCNT
!     DO K=1,NTGT(M)
!       WRITE(*,*) M,K,' ',CATG(MTGT(NTGT(M),M))
!     ENDDO
!     ENDDO

200 CONTINUE

!     2.3  DETERMINE SECTION OF CORRECT.T

  DO I_M=1,I_NCNT
    DO I_K=1,I_NTGT(I_M)
      ISON=0
NSOBC:DO I=1,I_NSOBC
        IF (CL_CATG(I_M) == CL_YSNAMBC(I)) THEN
          ISON=I
          EXIT NSOBC
        ENDIF
      ENDDO NSOBC

    IF(ISON == 0) THEN
      WRITE(*,'(A,A6,F7.2,F8.2,I6,I3,A,I3,2X,A)') &
     & CLR//': K_ID,LAT,LON,P_ST-GP&MB=', &
     & CD_CIDENT,P_RLAT,P_RLON,I_KG,I_KM, &
     & ' | NO COR. TABLE FOR CATEGORY NO.',I_M,CL_CATG(I_M)    
    ELSE
                                       
!     2.4  CALCULATE SOLAR ELEVATION

    CALL DIURNAL ( K_IM, K_ID, K_IH, P_RLAT,P_RLON, Z_ANGLEBC )
                                   
!     2.5  CALCULATE CORRECT COLUMN OF CORRECTION TABLE

    IF (Z_ANGLEBC < -7.5_JPRM)                          ICOL=1
    IF (Z_ANGLEBC >= -7.5_JPRM.AND.Z_ANGLEBC < 7.5_JPRM) ICOL=2
    IF (Z_ANGLEBC >= 7.5_JPRM.AND.Z_ANGLEBC < 22.5_JPRM) ICOL=3
    IF (Z_ANGLEBC >= 22.5_JPRM)                          ICOL=4

    if(K_IY .eq. 1989 .and. K_ID .eq. 1) THEN
     WRITE(*,'(A,A6,F7.2,F8.2,I6,I3,A,I3,2X,A,A,F6.2,A,I2,F8.2)') &
     & CLR//': K_ID,LAT,LON,P_ST-GP&MB=', &
     & CD_CIDENT,P_RLAT,P_RLON,I_KG,I_KM, &
     & ' | COR. TABLE NO.',ISON,CL_YSNAMBC(ISON), &
     & ' | SOL.EV.=',Z_ANGLEBC, &
     & ' | CLASS',ICOL, Z_RCORBC(ICOL,13,ISON)  

   ENDIF

!   IF(K_IY .eq. 1989) THEN
     KSON=ICOL
     P_ANGLE=Z_ANGLEBC

!   ENDIF
    
    if(LL_DEBUG) WRITE(*,'(A34,I6.6,A2,I2,'//CL_CPMAX//'F6.2)') CLR//':   RADIATION TABLE: ',ICSN,'  ',K_IH,Z_RCORBC(ICOL,:,ISON)


!     2.6  LOOP OVER LEVELS

    DO I=1,K_NST
      IF (P_SPP(I) /= P_RMISS.AND.P_ST(I) /= P_RMISS) THEN

!     2.7  CALCULATE CORRECT ROW OF CORRECTION TABLE

        PP=P_SPP(I)/100.0
        CLPUINT = 'hPa'
        CALL PNTERP (ILEVBC(1,ISON),I_NLEVBC,PP,Z_POS,CLPUINT)
        ILOW=Z_POS
        IHGH=ILOW+1


!  -------------------------------------------------------------
                                         
!                3.  Solar angle BIAS CORRECTION
!               ---------------------

        Z_RADD=Z_RCORBC(ICOL,ILOW,ISON)*(REAL(IHGH)-Z_POS) &
         & +Z_RCORBC(ICOL,IHGH,ISON)*(Z_POS-REAL(ILOW))  


        P_STNBC(I)=P_ST(I) - Z_RADD
        P_SOLBC(I)=P_ST(I) - Z_RADD
!                         !!!     

!   !!! NOTICE !!!

!   Please note that the values in the bias correction table have
!   opposite sign of correction (the same sign of biases)
!   Values taken from the tables must be subtracted.
!   This is opposite sign from operational bias correction.
!   In operation, values are added to observations.

!   We decided to use oppsite from operational one because we need q
!   manual adjustment of correction tables.
!   It was decided to avoid manual errors (confusion of signs).

!   Please allow us to use opposite sign from operational 'biascor'.

!                                                 2000.3.22 K.ONOGI

      ENDIF
    ENDDO

    IF(LL_DEBUG) WRITE(*,'(A31,I6.6,A2,I2,'//CL_F//'F6.2)') CLR//':   RADIATION ADJ ',ICSN,': ',K_IH,P_STNBC(1:K_NST)-P_ST(1:K_NST)

    LD_LBC=.TRUE.

!--------------------------------------------------------------------
                                              
!                   4.  EXIT
!                  -----------

    ENDIF

  ENDDO
ENDDO

ENDIF ! LDSOLAR

IF(LDHOMOGEN) THEN


!     5. Use corrections determined by homogenization of radiosondes
!

          READ(CD_CIDENT(1:5),'(I5)',IOSTAT=IERR) ICSN
          IF (IERR /= 0 ) THEN
            WRITE(6,*) CLR//': Could not make integer of ',CD_CIDENT(1:5)
          ENDIF
     
!         READ(CL_CSNG2(I_MREP(I_KG))(1:5),'(I5)') ISNG
!         WRITE(*,'(A31,I6.6,A2,3I6)') CLR//':    HOMOGENEITY ADJ ',ICSN,':',I_KG,I_KM,ISNG
!         ICSN=ISNG
         cleo_found=.false.
         ZLEO_BIAS=0.
         index=toindex(K_IY*10000+K_IM*100+K_ID,rcpara)
         ISTAT_LOOP : do istat=1,IESTATMAX
           if(ICSN == IEWMONRS(istat) .and. CL_CCORR(istat) == 'Y') then
!             WRITE(*,*)'Station found:',ICSN
           if(K_IH < 3 .or. K_IH >= 21) then
             ZLEO_BIAS=Z_RASOCORRS(index,:,1)
           endif
           if(K_IH >= 9 .and. K_IH < 15) then
             ZLEO_BIAS=Z_RASOCORRS(index,:,2)
           endif
           if(K_IH >= 3 .and. K_IH < 9 .or. K_IH >= 15 .and. K_IH < 21) then
             ZLEO_BIAS=0.5*(Z_RASOCORRS(index,:,1)+Z_RASOCORRS(index,:,2))
           endif
           cleo_found=.true.
           exit istat_loop
        endif
     enddo ISTAT_LOOP

     IF(cleo_found) then
       if(LL_DEBUG) WRITE(*,'(A34,I6.6,A2,I2,'//CL_CPMAX//'F6.2)') CLR//': HOMOGENEITY TABLE: ',ICSN,'  ',K_IH,ZLEO_BIAS
     ELSE
       if(LL_DEBUG) WRITE(*,'(A34,I6.6,A2,I2,A14)') CLR//': HOMOGENEITY TABLE: ',ICSN,'  ',K_IH,' NOT AVAILABLE'
     ENDIF

!     2.6  LOOP OVER LEVELS

    DO I=1,K_NST
      IF (P_SPP(I) /= P_RMISS.AND.P_ST(I) /= P_RMISS) THEN

!     2.7  CALCULATE CORRECT ROW OF CORRECTION TABLE

        PP=P_SPP(I)/100.0
        CLPUINT = 'hPa'
        CALL PNTERP (ILEVBC(1,1),I_NLEVBC,PP,Z_POS,CLPUINT)
        ILOW=Z_POS
        IHGH=ILOW+1


!  -------------------------------------------------------------
                                         
!                6.  BIAS CORRECTION using tables from homogenization
!               ---------------------

        Z_LEO(I)=ZLEO_BIAS(ILOW)*(REAL(IHGH)-Z_POS) &
         & +ZLEO_BIAS(IHGH)*(Z_POS-REAL(ILOW))  

        P_STNBC(I)=P_STNBC(I) - Z_LEO(I)
!                            !!!  

!   !!! NOTICE !!!

!   Please note that the values in the bias correction table have
!   opposite sign of correction (the same sign of biases)
!   Values taken from the tables must be subtracted.
!   This is opposite sign from operational bias correction.
!   In operation, values are added to observations.

!   We decided to use oppsite from operational one because we need 
!   manual adjustment of correction tables.
!   It was decided to avoid manual errors (confusion of signs).

!   Please allow us to use opposite sign from operational 'biascor'.

!                                                 2000.3.22 K.ONOGI

      ENDIF
    ENDDO

    IF(LL_DEBUG) WRITE(*,'(A31,I6.6,A2,I2,'//CL_F//'F6.2)') CLR//': HOMOGENEITY ADJ ',ICSN,': ',K_IH,-Z_LEO
    IF(LL_DEBUG) WRITE(*,'(A31,I6.6,A2,I2,'//CL_F//'F6.2)') CLR//':       TOTAL ADJ ',ICSN,': ',K_IH,P_STNBC-P_ST

    LD_LBC=.TRUE.

ENDIF

!          WRITE(*,*) CIDENT,LBC

! CALL OML_UNSET_LOCK()
!leoIF (LHOOK) CALL DR_HOOK('BIASCOR_ERA40',1,ZHOOK_HANDLE)
RETURN

965 CONTINUE
 CALL ABOR1('ERROR ON I/O OF STGROUP.T')
966 CONTINUE
 CALL ABOR1('ERROR ON I/O OF CORRECT.T')
967 CONTINUE
 CALL ABOR1('ERROR ON I/O OF COUNTRY.T')
968  CONTINUE
 CALL ABOR1('ERROR ON I/O OF BIASCOR.T')
!leo IF (LHOOK) CALL DR_HOOK('BIASCOR_ERA40',1,ZHOOK_HANDLE)

END SUBROUTINE BIASCOR_ERA40_changed

SUBROUTINE BIASCOR_ERA40_split (LD_LBC,K_NST,P_SPP,P_ST,P_STNBC, &
 & P_SOLBC,KSON,P_ANGLE,CD_CIDENT,K_IY,K_IM,K_ID,K_IH,P_RLAT,P_RLON,K_IRSTYP,  &
 & K_IMISS,P_RMISS,LDSOLAR,LDHOMOGEN,Z_RASOCORRS,rcpara)  

type(rasocor_namelist),intent(in) :: rcpara

!leoUSE PARKIND1  ,ONLY : JPIM     ,JPRM
!leoUSE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
!leoUSE YOMOML
                   
!**** *BIASCOR_ERA40* 

!       PURPOSE.
!      ------------%

!       BIAS CORRECTION OF RADIOSONDE TEMPERATURES FOR 

!          EEEEE  RRRR      A         4    000
!          E      R   R    A A       44   0   0
!          EEEE   RRRR    A   A     4 4   0   0
!          E      R R     AAAAA    44444  0   0
!          EEEEE  R   R   A   A       4    000

!       INTERFACE.
!      ------------

!          CALL BIASCOR_ERA40 (LBC,NST,SPP,ST,STNBC,
!     X                  CIDENT,IM,ID,IH,RLAT,RLON,IRSTYP,
!     X                  IMISS,RMISS)

!        INPUT
!         NST       NUMBER OF LEVELS (4 BYTE INTEGER)
!         SPP       ARRAY WITH PRESSURE VALUES (Pa) (8 BYTE REAL)
!         ST        ARRAY WITH T VALUES (8 BYTE REAL)
!         CIDENT    STATION IDENTIFIER  (CHARACTER)
!         IM        MONTH (4 BYTE INTEGER)
!         ID        DAY   (4 BYTE INTEGER)
!         IH        HOUR  (4 BYTE INTEGER)
!         RLAT      LATITUDE (8 BYTE REAL)
!         RLON      LONGITUDE (8 BYTE REAL)
!         IRSTYP    RADIOSONDE TYPE (BUFR CODE TABLE 002011) 
!                   (Not need for ERA)
!         IMISS     MISSING VALUE FOR INTEGERS (4BYTE INTEGER)
!         RMISS     MISSING VALUE FOR REALS (8 BYTE REAL) 

!        OUTPUT
!         LBC       LOGICAL TO STATE IF BIAS COR. WAS SUCCESSFUL
!         STNBC     ARRAY WITH BIAS CORRECTED T VALUES

!       METHOD.
!      ---------

!        READ BIAS CORRECTION TABLES:
!            STGROUP.T : STATION GROUP TABLE
!            CORRECT.T : ARRAY OF CORRECTIONS DEPENDING ON
!                        SOLAR ELEVATION AND PRESSURE LEVEL
!            COUNTRY.T : DEFINITION OF CATEGORIES

!        1ST A STATION GROUP 0F THE DATA IS DETECTED.
!        2ND A CATEGORY INCLUDING THE STATION GROUP IS DETECTED.
!        3RD IF THE CATEGORY IS ONE FOR CORRECTION, APPLY CORRECTION.

!        FIRST 32 CHARACTERS ARE COMPARED TO DETECT CATEGORY.
!        'CATG' FROM COUNTRY.T AND 'YMNAMBC' FROM CORRECT.T

!       EXTERNALS.
!      ------------

!        DIURNAL          CALCULATE SOLAR ELEVATION
!        PNTERP           INTERPOLATE TO PRESSURE LEVEL

!       REFERENCE.
!      ------------

!        RADIOSONDE BIAS CORRECTION
!        OD MEMORANDUM BY B. STRAUSS  22.12.92

!       AUTHOR.
!      ---------

!        B. NORRIS       JULY  1991   APPLY CORRECTIONS TO AOF

!       MODIFICATIONS.
!      ----------------

!        M. DRAGOSAVAC   APRIL 1993   IMPLEMENT IN PRE-PROCESSING
!        B. NORRIS       JULY  1998   INTERFACE TO DATA ASSIMILATION
!        K. ONOGI        MARCH 2000   MODIFIED FOR ERA-40
!        K. ONOGI      OCTOBER 2000   MAKE NUMBER OF B.COR. CATEGORIES 8 TO 4
!                                     AS SAME AS THE CORRECTION TABLE
!        S. SAARINEN  NOVEMBER 2000   CLEANING UP PLUS IMPLICIT NONE
!        M. Hamrud      01-Oct-2003   CY28 Cleaning
!      L. Haimberger  NOVEMBER 2004   Allow use of tables from  radiosonde homogenization 
!                                    (see ERA-40 Project report series 23)
!        S. SAARINEN       MAY 2005   Removed standalone SAVE-stmt to allow Dr.Hook
!-------------------------------------------------------------------

!IMPLICIT NONE

!     Table values are expected to be saved
!     after they are read

! INTEGER,PARAMETER :: JPIM=4
! INTEGER,PARAMETER :: JPRM=4
! INTEGER,PARAMETER :: JPRM=4

LOGICAL           ,INTENT(OUT)   :: LD_LBC 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_NST 
REAL(KIND=JPRM)   ,INTENT(IN)    :: P_SPP(K_NST) 
REAL(KIND=JPRM)   ,INTENT(IN)    :: P_ST(K_NST) 
REAL(KIND=JPRM)   ,INTENT(OUT)   :: P_STNBC(K_NST) 
REAL(KIND=JPRM)   ,INTENT(OUT)   :: P_SOLBC(K_NST) 
CHARACTER(LEN=5)  ,INTENT(IN)    :: CD_CIDENT 
INTEGER(KIND=JPIM),INTENT(OUT)    :: KSON 
REAL(KIND=JPRM),INTENT(OUT)    :: P_ANGLE 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_IM 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_ID 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_IH 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_IY 
REAL(KIND=JPRM)   ,INTENT(IN)    :: P_RLAT 
REAL(KIND=JPRM)   ,INTENT(IN)    :: P_RLON 
INTEGER(KIND=JPIM)               :: K_IRSTYP ! Argument NOT used
INTEGER(KIND=JPIM)               :: K_IMISS ! Argument NOT used
REAL(KIND=JPRM)   ,INTENT(IN)    :: P_RMISS 
LOGICAL   ,INTENT(IN)            :: LDSOLAR,LDHOMOGEN

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

REAL(KIND=JPRM)     , DIMENSION (I_NCORBC,I_NLEVBC,I_NSONBC) ,SAVE:: Z_RCORBC
INTEGER(KIND=JPIM)  , DIMENSION (I_NLEVBC,I_NSONBC) ,SAVE:: ILEVBC
REAL(KIND=JPRM)     , DIMENSION (I_NCORTB) ,SAVE:: Z_WKCOR
INTEGER(KIND=JPIM) ,SAVE:: IPLAT,IPLON,IRLAT,IRLON
CHARACTER(LEN= 5) ,SAVE:: CL_YDUMMYBC,CL_VERSION,CL_ISTAT
CHARACTER(LEN=32) , DIMENSION (I_NSONBC) ,SAVE:: CL_YSNAMBC
CHARACTER(LEN=58) ,SAVE:: CL_CDATE
CHARACTER(LEN= 2) ,SAVE:: CL_CPMAX,CLHOUR
CHARACTER(LEN= 3) ,SAVE:: CL_F
CHARACTER(LEN= 8)      :: CLPUINT
CHARACTER(LEN=13)      :: CLR='BIASCOR_ERA40'
CHARACTER(LEN=4) :: CLYEAR

LOGICAL ,SAVE:: LL_FILES_NOT_READ,LL_DEBUG=.false.
INTEGER(KIND=JPIM) ,SAVE:: I_NSOBC
      
! --- ARRAYS FOR STATION GROUP TABLE
REAL(KIND=JPRM)     , DIMENSION (I_NXGC) ,SAVE:: Z_PLAT,Z_PLON
CHARACTER(LEN= 6), DIMENSION (I_NXGC) ,SAVE:: CL_CSNG,CL_CSNG2,CL_CSNO
CHARACTER(LEN= 6) ,SAVE:: CL_C1,CL_C4
CHARACTER(LEN= 1) ,SAVE:: CL_C0
INTEGER(KIND=JPIM)  , DIMENSION (I_NXGP) ,SAVE:: JSGC,JEGC,I_MREP
REAL(KIND=JPRM)     , DIMENSION (I_NXGP) ,SAVE:: Z_QLAT,Z_QLON
INTEGER(KIND=JPIM)  , DIMENSION (I_NXCT) ,SAVE:: I_NTGT
INTEGER(KIND=JPIM)  , DIMENSION (I_NXGP,I_NXCT) ,SAVE:: I_MTGT

! --- ARRAYS FOR INDEX TABLE
INTEGER(KIND=JPIM)  , DIMENSION (I_NXDG,I_NXCT) ,SAVE:: IDSTA,IDEND, &
 & I_LATS,I_LATE,I_LONS,I_LONE  
INTEGER(KIND=JPIM)  , DIMENSION (I_NXDG) ,SAVE:: I_NCMB,I_NCDU
INTEGER(KIND=JPIM)  , DIMENSION (I_NXCT) ,SAVE:: I_NCCC,I_MCAT
CHARACTER(LEN=64), DIMENSION (I_NXDG,I_NXCT) ,SAVE:: CL_CNTRY
CHARACTER(LEN=64), DIMENSION (I_NXCT) ,SAVE:: CL_CDG
CHARACTER(LEN=64) ,SAVE::  CL_C64
CHARACTER(LEN=18) ,SAVE::  CLTLN
CHARACTER(LEN=32), DIMENSION (I_NXCT) ,SAVE:: CL_CATG

LOGICAL, PARAMETER :: LGRP = .TRUE.

! --- LEO VARIABLES
INTEGER(KIND=JPIM) :: ierr,icsn
LOGICAL :: CLEO_FOUND = .FALSE.
CHARACTER(LEN=50) :: CLEONAME = 'table4'
CHARACTER(LEN=120) :: CL_ZEILE
INTEGER(KIND=JPIM), PARAMETER :: IPMAX    = 16
INTEGER(KIND=JPIM), PARAMETER :: IPARMAX   = 2
INTEGER(KIND=JPIM), PARAMETER :: IESTATMAX = 3070
INTEGER(KIND=JPIM), PARAMETER :: NMAX = 18993
CHARACTER(LEN=1),SAVE  :: CL_CCORR(IESTATMAX)
REAL(KIND=JPRM) :: Z_RASOCORRS(NMAX,IPMAX,IPARMAX)
INTEGER(KIND=JPIM),SAVE :: IEWMONRS(IESTATMAX)
INTEGER(KIND=JPIM),SAVE :: ISTAT,IB
REAL(KIND=JPRM),SAVE    :: ZLEO_PLEV(IPMAX),ZIPLEVS(IPMAX),ZLEO_BIAS(IPMAX)
real(kind=JPRM),SAVE    :: ZWMOLATS(IESTATMAX),ZWMOLONS(IESTATMAX)
real(kind=JPRM)    :: ZHILF(IPMAX,IPARMAX),ZHILFOLD(IPMAX,IPARMAX)
integer(KIND=JPIM),SAVE :: iyear,imonth,itime,iday,idatum,idatumold,i_datum,idum,ip,IIPMAX,istatmax

      
! --- MISCELLANEOUS DECLARATIONS
INTEGER(KIND=JPIM) ,SAVE:: IOS, I, IR,JDP,IDS,IDE,I_LAS,I_LAE,I_LOS,I_LOE
INTEGER(KIND=JPIM) ,SAVE:: I_M, JDG, I_N, I_NGC, I_NGP, I1, I2,I6
REAL(KIND=JPRM) ,SAVE:: Z_R1,Z_R2
INTEGER(KIND=JPIM) ,SAVE:: ISY,ISM,ISD,ISH,IEY,IEM,IED,IEH,I5
INTEGER(KIND=JPIM) ,SAVE:: I_NCNT, J, I_KG, I_KM, IXT, ISON, I_K, ISNG, ICOL
REAL(KIND=JPRM) ,SAVE:: Z_ANGLEBC, PP, Z_POS, Z_RADD
INTEGER(KIND=JPIM) ,SAVE:: ILOW, IHGH,ianf,iend,ipar
REAL(KIND=JPRM) :: ZHOOK_HANDLE
REAL(KIND=JPRM) :: Z_LEO(K_NST)

LOGICAL EX

!leo#include "abor1.intfb.h"
!leo#include "diurnal.intfb.h"
!leo#include "pnterp.intfb.h"

!  -------------------------------------------------------------------

!     1.  OPEN AND READ TABLE FILES
!     -------------------------------

!leo IF (LHOOK) CALL DR_HOOK(CLR,0,ZHOOK_HANDLE)

! CALL OML_SET_LOCK()
100 CONTINUE

WRITE(CL_F,'(I3.3)') K_NST
LD_LBC=.FALSE.
IF(K_IMISS == 1) THEN
  LL_FILES_NOT_READ=.true.
  K_IMISS=0
  I_KM=0
  I_KG=0
ENDIF
if(ll_files_not_read) then
!         OPEN(UNIT=MCOR, FILE='/era/data/erk/comtable/correct.t', &

    WRITE(CLYEAR,'(I4)') K_IY
    WRITE(CLHOUR,'(I2.2)') K_IH
    INQUIRE(FILE='/home/srvx2/leo/tables/ulf_split/T_correct_'//CLHOUR//'_'//CLYEAR//'010100',EXIST=EX)
    if(ex) then 
    OPEN(UNIT=I_MCOR,FILE='/home/srvx2/leo/tables/ulf_split/T_correct_'//CLHOUR//'_'//CLYEAR//'010100' , &
   & IOSTAT=IOS,ERR=966,STATUS='OLD')  

!       1.3  READ CORRECTION TABLE

  I_NSOBC=0
  READ (I_MCOR,'(A)') CL_CDATE

  IF(LDSOLAR .and. LDHOMOGEN) THEN

!  write(*,*) 'reading biascor.t'
! If L. Haimbergers correction table is used, it has to be used with a special version of
! solar angle bias correction tables. See ERA-40 project report series Nr. 23 by L. Haimberger
    IF(CL_CDATE(52:52) /= '1') THEN
      CALL ABOR1(CLR//': wrong table2 for solar angle correction + homogenization')
    ENDIF
  ELSE
    IF(LDSOLAR .and. CL_CDATE(52:52) == '1') THEN
      CALL ABOR1(CLR//': wrong table2 for solar correction only')
    ENDIF
  ENDIF

  READ (I_MCOR,*)

  DO WHILE(.TRUE.)
    IF (I_NSOBC < I_NSONBC) THEN
      I_NSOBC=I_NSOBC+1
      READ (I_MCOR,'(A)',END=138) CL_YSNAMBC(I_NSOBC)
!           WRITE(*,'(2A)') YSNAMBC(NSOBC)

      DO J=1,I_NLEVBC
        READ (I_MCOR,'(I5,8F8.2)',END=138) &
       & ILEVBC(J,I_NSOBC),(Z_WKCOR(I_K),I_K=1,I_NCORTB)  

!         WRITE(*,'(I5,8F8.2)') &
!               & ILEVBC(J,NSOBC),(WKCOR(K),K=1,NCORTB)

!     CORRECTION TABLE
!     ~-7.5 : -7.5~7.5 : 7.5~22.5 : 22.5~

        DO I=1,I_NCORBC
          Z_RCORBC(I,J,I_NSOBC) = Z_WKCOR(I)
        ENDDO

      ENDDO
      READ (I_MCOR,'(A)',END=138) CL_YDUMMYBC
    ELSE
      WRITE (0,'(1X,A,I5)') &
     & CLR//': DIMENSION NSON TOO SMALL :', &
     & I_NSONBC  
      CALL ABOR1(CLR//': DIMENSION NSON TOO SMALL')
      CALL ABORT
    ENDIF
  ENDDO
  138   CONTINUE
  I_NSOBC = I_NSOBC-1

!         WRITE (0,*) &
!         & CLR//': NUMBER OF CATEGORIES IN CORRECTION TABLE ', &
!         &  NSOBC


        !!WRITE(6,*)'SUCCESSFULLY READ TABLE 1-3'
  CLOSE (I_MCOR)
  else
    Z_RCORBC=0.
    CL_CDATE =''
  endif

  LL_FILES_NOT_READ = .FALSE.

IF (CL_ISTAT .ne. CD_CIDENT ) THEN

!        WRITE(*,*)
!        WRITE(*,*) '------------------------------'
!        WRITE(*,*) '  RADIOSONDE BIAS CORRECTION  '
!        WRITE(*,*) '------------------------------'
!        WRITE(*,*)


!         OPEN(UNIT=MTBL, FILE='/era/data/erk/comtable/country.t', &
    OPEN(UNIT=I_MTBL, FILE='/home/srvx9/rs_common/tables/ulf_split/country.t', &
   & IOSTAT=IOS,ERR=965,STATUS='OLD')  


!         OPEN(UNIT=MSGT, FILE='/era/data/erk/comtable/stgroup.t', &
    OPEN(UNIT=I_MSGT, FILE='/home/srvx9/rs_common/tables/stgroup.t', &
   & IOSTAT=IOS,ERR=967,STATUS='OLD')  



!       1.1 READ CATEGORY DEFINITION TABLE

  I=1
  IR=0
  DO WHILE(I <= 2000 .and. IR /=10)
    READ(I_MTBL,'(I4)') IR
    I=I+1
  ENDDO
  I_NCNT = 0
  DO I=1,I_NXDG
    I_NCMB(I) = 0
  ENDDO
  NXCT: DO I=1,I_NXCT
    READ(I_MTBL,'(I4,I3,2I6,1X,2I3,2I4,1X,A64)') &
     & IR,JDP,IDS,IDE,I_LAS,I_LAE,I_LOS,I_LOE,CL_C64  
    IF (IR == 999) EXIT NXCT
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
  ENDDO NXCT

!       1.1.1 PUT LAT LON LIMIT ON CNTRY

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
!               CNTRY(I,M) = CNTRY(I,M)(1:20)//CLTLN
        CL_CNTRY(I,I_M) = CLTLN//CL_CNTRY(I,I_M)(1:20)

      ENDIF
    ENDDO
  ENDDO

!       1.1.2 GATHER SOME GROUPS INTO ONE GROUP

  DO I=1,I_NXCT
    CL_CDG(I) = ' '
  ENDDO
  IF (LGRP) THEN
    IR=0
    I=1
    DO WHILE(I <=2000 .and. IR /=20)
      READ(I_MTBL,'(I4)') IR
      I=I+1
    ENDDO
NXDG:    DO I=1,I_NXDG
      READ(I_MTBL,'(I4,I3,1X,A64)') IR,JDG,CL_C64
      IF (IR == 999) EXIT NXDG
      IF (IR == 1) THEN
        CL_CDG(I_NCDU(JDG)) = CL_C64
      ENDIF
    ENDDO NXDG
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

!       1.1.3 MONITOR THE CATEGORIES
!       DO N=1,NCNT
!         DO I=1,NXDG
!         IF (IDSTA(I,N)/=0) &
!       & WRITE(*,'(I6.6,2I7,4I5,1X,A32,2X,A32)') &
!             & I,IDSTA(I,N),IDEND(I,N),LATS(I,N),LATE(I,N), &
!             & LONS(I,N),LONE(I,N),CNTRY(I,N)(1:32),CATG(N)
!         ENDDO
!       ENDDO

!       1.2 READ STATION GROUP TABLE
  I_NGC = 0  ! NUMBER OF RECORD
  I_NGP = 0  ! NUMBER OF GROUP
  DO WHILE(.TRUE.)
  READ(I_MSGT, &
   & '(I6.6,I4,A1,A6,F7.2,F8.2,I8,2(1X,I4,3I2),I7,2X,A6)', &
   & END=190) &
   & I1,I2,CL_C0,CL_C1,Z_R1,Z_R2,I6, &
   & ISY,ISM,ISD,ISH,IEY,IEM,IED,IEH,I5,CL_C4  
!           WRITE(*, &
!           & '(I6.6,I4,A1,A6,F7.2,F8.2,A5,2(1X,4I2),I7,2X,A19,2X,A6)') &
!           &  I1,I2,C0,C1,R1,R2,C2, &
!           &  ISY,ISM,ISD,ISH,IEY,IEM,IED,IEH,I5,C3,C4
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
  CL_CSNG2(I_NGC) = CL_C1
  IF (CL_C0 == '#') THEN
    I_MREP(I_NGP) = I_NGC
    Z_QLAT(I_NGP) = Z_R1
    Z_QLON(I_NGP) = Z_R2
  ENDIF
  ENDDO
  190   CONTINUE



  CLOSE (I_MTBL)
  CLOSE (I_MSGT)


  IF(LDHOMOGEN) THEN

!         OPEN(UNIT=MSGT, FILE='biascor.t', &
  OPEN(UNIT=I_MLEO, FILE='leobiascor.t', &
   & IOSTAT=IOS,ERR=968,STATUS='OLD')  

! -------------------------------------
! 
!     Read L. Haimbergers's homogenizing temp corrections
!     http://www.univie.ac.at/theoret-met/research/RAOBCORE/     
! -------------------------------------
   
 ZLEO_PLEV=(/10.,20.,30.,50.,70.,100.,150.,200.,250.,300.,400.,500.,700.,850.,925.,1000./)
 do ip=1,IPMAX
    ilevbc(ip,1)=ZLEO_PLEV(IPMAX-ip+1)
 enddo

      ZLEO_BIAS=0.

! read and print comment line(s) in header, if available
  CL_ZEILE='#'
  do while(CL_ZEILE(1:1) == '#')
    read(I_MLEO,'(A120)') CL_ZEILE
    if(CL_ZEILE(1:1) == '#') write(*,*) CL_ZEILE
  enddo

  read(CL_ZEILE,'(A5,I5,I3)')  CL_VERSION,istatmax,IIPMAX
  write(CL_CPMAX,'(I2.2)') IIPMAX
  read(CL_ZEILE,'(A5,I5,I3,'//CL_CPMAX//'F6.0)') CL_VERSION,istatmax,IIPMAX,ZIPLEVS

  if(IPMAX /= IIPMAX .or. any(ZLEO_PLEV /= ZIPLEVS) .or. istatmax /= IESTATMAX) then
     write(6,*) CLR//': istatmax,IESTATMAX',istatmax,IESTATMAX
     write(6,*) CLR//': IPMAX,IIPMAX',IPMAX,IIPMAX
     write(6,*) CLR//': plevs',ZLEO_PLEV
     write(6,*) CLR//': ZIPLEVS',ZIPLEVS
     CALL ABOR1(CLR//': WRONG DIMENSIONS OF BIASCOR.T')
  endif
  if(LDSOLAR .and. CL_VERSION /= CL_CDATE(54:58) .and. CL_CDATE /='' ) then
    CALL ABOR1(CLR//': table2 and table4 version mismatch:'//CL_VERSION//':'//CL_CDATE(54:58)//':')
  endif

  I_DATUM=K_IY*1000000+K_IM*10000+K_ID*100

  z_rasocorrs=0.
  do istat=1,IESTATMAX
    read(I_MLEO,'(I5,I6,2F9.2,1X,A1,I3)',ERR=968) idum,IEWMONRS(istat),ZWMOLATS(istat),ZWMOLONS(istat),CL_CCORR(istat),ib

    if(CL_CCORR(istat) == "Y") then
      ZHILF=0.
      idatum=1957010100
      do i=1,ib
        ZHILFOLD=ZHILF
        idatumold=idatum
        read(I_MLEO,'(I4,3I2.2,'//CL_CPMAX//'F6.2)',ERR=968) iyear,imonth,iday,itime,ZHILF(:,1)
        read(I_MLEO,'(I4,3I2.2,'//CL_CPMAX//'F6.2)',ERR=968) iyear,imonth,iday,itime,ZHILF(:,2)
        idatum=iyear*1000000+imonth*10000+iday*100
        ianf=toindex(idatumold/100,rcpara)
        iend=toindex(idatum/100,rcpara)
        write(CL_ISTAT,'(I6.6)') iewmonrs(istat)
        IF(CD_CIDENT .eq. CL_ISTAT) then
          do ip=1,rcpara%PMAX
            do ipar=1,rcpara%PARMAX
            Z_RASOCORRS(ianf:iend,ip,ipar)=ZHILFOLD(IPMAX-ip+1,ipar)
            enddo
          enddo
        ENDIF
      enddo
!     write(*,*) 'read ',IEWMONRS(istat),ZWMOLATS(istat),ZWMOLONS(istat),CL_CCORR(istat),ib,Z_RASOCORRS(5,:,istat)
    endif

  enddo
  
  CL_ISTAT=CD_CIDENT

  CLOSE (I_MLEO)

  ENDIF ! LDHOMOGEN

ENDIF

ENDIF

!-----------------------------------------------------------------------------


!     1.4  INITIALLY SET CORRECTED = ORIGINAL

DO I=1,K_NST
  P_STNBC(I)=P_ST(I)
  P_SOLBC(I)=P_ST(I)
ENDDO

!  ------------------------------------------------------------------ 

!  Solar angle bias correction
!
!     2.   DETERMINE TABLES AND TARGET

!     2.1  DETECT THE STATION GROUP

IF(I_KG == 0 .and. I_KM == 0) THEN
I_KG = 0 ; I_KM = 0
NGP: DO I=1,I_NGP
  DO J=JSGC(I),JEGC(I)
    IPLAT = NINT(Z_PLAT(J)*100)
    IPLON = NINT(Z_PLON(J)*100)
    IRLAT = NINT(P_RLAT*100)
    IRLON = NINT(P_RLON*100)
    IF (abs(IRLAT-IPLAT) <50 .AND. abs(IRLON - IPLON) <50) THEN
      IF (CD_CIDENT(1:5) == CL_CSNO(J)(1:5)) THEN
        I_KG = I               ! STATION GROUP
        I_KM = J-JSGC(I)+1     ! MEMBER 
        EXIT NGP
      ENDIF
    ENDIF
  ENDDO
ENDDO NGP

IF (I_KG*I_KM == 0) THEN
  WRITE(*,'(A,A6,F7.2,F8.2,I6,I3,2A)') &
   & CLR//': K_ID,LAT,LON,GRP,MEM= ',CD_CIDENT,P_RLAT,P_RLON,I_KG,I_KM, &
   & ' | THE SPECIFIED DATA WAS NOT FOUND ', &
   & 'IN THE STATION GROUP DEFINITION TABLE.'  
      ! THE SPECIFIED DATA WAS NOT FOUND IN THE STATION GROUP DEFINITION TABLE.
! CALL OML_UNSET_LOCK()
!leo  IF (LHOOK) CALL DR_HOOK('BIASCOR_ERA40',1,ZHOOK_HANDLE)
!  RETURN

ENDIF

ENDIF

IF(LDSOLAR .AND. I_KG*I_KM > 0 ) THEN

!     2.2  PICK UP THE TARGET CATEGORY

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

!     DO M=1,NCNT
!     DO K=1,NTGT(M)
!       WRITE(*,*) M,K,' ',CATG(MTGT(NTGT(M),M))
!     ENDDO
!     ENDDO

200 CONTINUE

!     2.3  DETERMINE SECTION OF CORRECT.T

  DO I_M=1,I_NCNT
    DO I_K=1,I_NTGT(I_M)
      ISON=0
NSOBC:DO I=1,I_NSOBC
        IF (CL_CATG(I_M) == CL_YSNAMBC(I)) THEN
          ISON=I
          EXIT NSOBC
        ENDIF
      ENDDO NSOBC

    IF(ISON == 0) THEN
      WRITE(*,'(A,A6,F7.2,F8.2,I6,I3,A,I3,2X,A)') &
     & CLR//': K_ID,LAT,LON,P_ST-GP&MB=', &
     & CD_CIDENT,P_RLAT,P_RLON,I_KG,I_KM, &
     & ' | NO COR. TABLE FOR CATEGORY NO.',I_M,CL_CATG(I_M)    
    ELSE
                                       
!     2.4  CALCULATE SOLAR ELEVATION

    CALL DIURNAL ( K_IM, K_ID, K_IH, P_RLAT,P_RLON, Z_ANGLEBC )
                                   
!     2.5  CALCULATE CORRECT COLUMN OF CORRECTION TABLE

    IF (Z_ANGLEBC < -7.5_JPRM)                          ICOL=1
    IF (Z_ANGLEBC >= -7.5_JPRM.AND.Z_ANGLEBC < 7.5_JPRM) ICOL=2
    IF (Z_ANGLEBC >= 7.5_JPRM.AND.Z_ANGLEBC < 22.5_JPRM) ICOL=3
    IF (Z_ANGLEBC >= 22.5_JPRM)                          ICOL=4

    if(K_IY .eq. 1989 .and. K_ID .eq. 1) THEN
     WRITE(*,'(A,A6,F7.2,F8.2,I6,I3,A,I3,2X,A,A,F6.2,A,I2,F8.2)') &
     & CLR//': K_ID,LAT,LON,P_ST-GP&MB=', &
     & CD_CIDENT,P_RLAT,P_RLON,I_KG,I_KM, &
     & ' | COR. TABLE NO.',ISON,CL_YSNAMBC(ISON), &
     & ' | SOL.EV.=',Z_ANGLEBC, &
     & ' | CLASS',ICOL, Z_RCORBC(ICOL,13,ISON)  

   ENDIF

!   IF(K_IY .eq. 1989) THEN
     KSON=ICOL
     P_ANGLE=Z_ANGLEBC

!   ENDIF
    
    if(LL_DEBUG) WRITE(*,'(A34,I6.6,A2,I2,'//CL_CPMAX//'F6.2)') CLR//':   RADIATION TABLE: ',ICSN,'  ',K_IH,Z_RCORBC(ICOL,:,ISON)


!     2.6  LOOP OVER LEVELS

    DO I=1,K_NST
      IF (P_SPP(I) /= P_RMISS.AND.P_ST(I) /= P_RMISS) THEN

!     2.7  CALCULATE CORRECT ROW OF CORRECTION TABLE

        PP=P_SPP(I)/100.0
        CLPUINT = 'hPa'
        CALL PNTERP (ILEVBC(1,ISON),I_NLEVBC,PP,Z_POS,CLPUINT)
        ILOW=Z_POS
        IHGH=ILOW+1


!  -------------------------------------------------------------
                                         
!                3.  Solar angle BIAS CORRECTION
!               ---------------------

        Z_RADD=Z_RCORBC(ICOL,ILOW,ISON)*(REAL(IHGH)-Z_POS) &
         & +Z_RCORBC(ICOL,IHGH,ISON)*(Z_POS-REAL(ILOW))  


        P_STNBC(I)=P_ST(I) - Z_RADD
        P_SOLBC(I)=P_ST(I) - Z_RADD
!                         !!!     

!   !!! NOTICE !!!

!   Please note that the values in the bias correction table have
!   opposite sign of correction (the same sign of biases)
!   Values taken from the tables must be subtracted.
!   This is opposite sign from operational bias correction.
!   In operation, values are added to observations.

!   We decided to use oppsite from operational one because we need q
!   manual adjustment of correction tables.
!   It was decided to avoid manual errors (confusion of signs).

!   Please allow us to use opposite sign from operational 'biascor'.

!                                                 2000.3.22 K.ONOGI

      ENDIF
    ENDDO

    IF(LL_DEBUG) WRITE(*,'(A31,I6.6,A2,I2,'//CL_F//'F6.2)') CLR//':   RADIATION ADJ ',ICSN,': ',K_IH,P_STNBC(1:K_NST)-P_ST(1:K_NST)

    LD_LBC=.TRUE.

!--------------------------------------------------------------------
                                              
!                   4.  EXIT
!                  -----------

    ENDIF

  ENDDO
ENDDO

ENDIF ! LDSOLAR

IF(LDHOMOGEN) THEN


!     5. Use corrections determined by homogenization of radiosondes
!

          READ(CD_CIDENT(1:5),'(I5)',IOSTAT=IERR) ICSN
          IF (IERR /= 0 ) THEN
            WRITE(6,*) CLR//': Could not make integer of ',CD_CIDENT(1:5)
          ENDIF
     
!         READ(CL_CSNG2(I_MREP(I_KG))(1:5),'(I5)') ISNG
!         WRITE(*,'(A31,I6.6,A2,3I6)') CLR//':    HOMOGENEITY ADJ ',ICSN,':',I_KG,I_KM,ISNG
!         ICSN=ISNG
         cleo_found=.false.
         ZLEO_BIAS=0.
         index=toindex(K_IY*10000+K_IM*100+K_ID,rcpara)
         ISTAT_LOOP : do istat=1,IESTATMAX
           if(ICSN == IEWMONRS(istat) .and. CL_CCORR(istat) == 'Y') then
!             WRITE(*,*)'Station found:',ICSN
           if(K_IH < 3 .or. K_IH >= 21) then
             ZLEO_BIAS=Z_RASOCORRS(index,:,1)
           endif
           if(K_IH >= 9 .and. K_IH < 15) then
             ZLEO_BIAS=Z_RASOCORRS(index,:,2)
           endif
           if(K_IH >= 3 .and. K_IH < 9 .or. K_IH >= 15 .and. K_IH < 21) then
             ZLEO_BIAS=0.5*(Z_RASOCORRS(index,:,1)+Z_RASOCORRS(index,:,2))
           endif
           cleo_found=.true.
           exit istat_loop
        endif
     enddo ISTAT_LOOP

     IF(cleo_found) then
       if(LL_DEBUG) WRITE(*,'(A34,I6.6,A2,I2,'//CL_CPMAX//'F6.2)') CLR//': HOMOGENEITY TABLE: ',ICSN,'  ',K_IH,ZLEO_BIAS
     ELSE
       if(LL_DEBUG) WRITE(*,'(A34,I6.6,A2,I2,A14)') CLR//': HOMOGENEITY TABLE: ',ICSN,'  ',K_IH,' NOT AVAILABLE'
     ENDIF

!     2.6  LOOP OVER LEVELS

    DO I=1,K_NST
      IF (P_SPP(I) /= P_RMISS.AND.P_ST(I) /= P_RMISS) THEN

!     2.7  CALCULATE CORRECT ROW OF CORRECTION TABLE

        PP=P_SPP(I)/100.0
        CLPUINT = 'hPa'
        CALL PNTERP (ILEVBC(1,1),I_NLEVBC,PP,Z_POS,CLPUINT)
        ILOW=Z_POS
        IHGH=ILOW+1


!  -------------------------------------------------------------
                                         
!                6.  BIAS CORRECTION using tables from homogenization
!               ---------------------

        Z_LEO(I)=ZLEO_BIAS(ILOW)*(REAL(IHGH)-Z_POS) &
         & +ZLEO_BIAS(IHGH)*(Z_POS-REAL(ILOW))  

        P_STNBC(I)=P_STNBC(I) - Z_LEO(I)
!                            !!!  

!   !!! NOTICE !!!

!   Please note that the values in the bias correction table have
!   opposite sign of correction (the same sign of biases)
!   Values taken from the tables must be subtracted.
!   This is opposite sign from operational bias correction.
!   In operation, values are added to observations.

!   We decided to use oppsite from operational one because we need 
!   manual adjustment of correction tables.
!   It was decided to avoid manual errors (confusion of signs).

!   Please allow us to use opposite sign from operational 'biascor'.

!                                                 2000.3.22 K.ONOGI

      ENDIF
    ENDDO

    IF(LL_DEBUG) WRITE(*,'(A31,I6.6,A2,I2,'//CL_F//'F6.2)') CLR//': HOMOGENEITY ADJ ',ICSN,': ',K_IH,-Z_LEO
    IF(LL_DEBUG) WRITE(*,'(A31,I6.6,A2,I2,'//CL_F//'F6.2)') CLR//':       TOTAL ADJ ',ICSN,': ',K_IH,P_STNBC-P_ST

    LD_LBC=.TRUE.

ENDIF

!          WRITE(*,*) CIDENT,LBC

! CALL OML_UNSET_LOCK()
!leoIF (LHOOK) CALL DR_HOOK('BIASCOR_ERA40',1,ZHOOK_HANDLE)
RETURN

965 CONTINUE
 CALL ABOR1('ERROR ON I/O OF STGROUP.T')
966 CONTINUE
 CALL ABOR1('ERROR ON I/O OF CORRECT.T')
967 CONTINUE
 CALL ABOR1('ERROR ON I/O OF COUNTRY.T')
968  CONTINUE
 CALL ABOR1('ERROR ON I/O OF BIASCOR.T')
!leo IF (LHOOK) CALL DR_HOOK('BIASCOR_ERA40',1,ZHOOK_HANDLE)

END SUBROUTINE BIASCOR_ERA40_split

subroutine grib_and_fill_noaa(rcpara,te20c,e20path,fieldsm,wmolons,wmolats,wmostats,ni,nj,nk,iens,iy,iref)

 use data_setting

implicit none

  INTEGER :: current_year
  INTEGER :: current_varno
  TYPE(twenty_CR):: tw_CR_station,tw_CR_station_s
  INTEGER :: time_dim !time dimension = numbers of items

type(rasocor_namelist) rcpara
real(kind=jprm),intent(in) :: wmolons(:),wmolats(:)
real(kind=jprm) :: fieldsm(:,:,:,:,:)
type(cacherecord) :: te20c(:,:)
real:: fields(ni,nj,rcpara%pmax*366*4)
integer,intent(in) :: ni,nj,nk,iens,iy,wmostats,iref
integer :: istat,iend,istart,mstride,mpar(1),i,k,il,nnk,status,fcount,idate,im,imold,ip,it,index
integer :: idxp(rcpara%pmax)
character cdate*4,filename*80,cens*1
character*(*) :: e20path
logical:: ex3

           mstride=1
           write(cens,'(I1.1)') iens-1
           write(cdate,'(I4)') iy
             filename=trim(e20path)//'/air.'//cdate//'.nc'
             write(*,*) filename
             inquire(file=filename,exist=ex3)
             if(ex3) then
current_varno = 0 !this is the index for geopotential height
status = read_20CR_PL(iy,2,tw_CR_station,time_dim)!in file new_station.f90

 idxp=[24,23,22,21,20,19,18,17,16,15,13,11,7,4,2,1]

 istart=toindex(iy*10000+100+1,rcpara)
 iend=toindex((iy+1)*10000+100+1,rcpara)-1
 index=istart
 fieldsm(:,:,:,(iy-iref)*12+1:(iy-iref+1)*12,iens)=0._JPRB
 imold=1
 fcount=0
 do it=1,time_dim
   idate=todate(index+(it-1)/4,rcpara)
   im=(idate-iy*10000)/100
   if(imold .ne. im) then
     fieldsm(:,:,:,il,iens)=fieldsm(:,:,:,il,iens)/fcount
     fcount=0
   endif
   fcount=fcount+1
   il=(iy-iref)*12+im
   do ip=1,rcpara%pmax
     fieldsm(:,:,ip,il,iens)=fieldsm(:,:,ip,il,iens)+tw_CR_station%short_data_pl(:,:,idxp(ip),it)*tw_CR_station%scale_factor+tw_CR_station%offset
     fields(:,:,(it-1)*rcpara%pmax+ip)=tw_CR_station%short_data_pl(:,:,idxp(ip),it)*tw_CR_station%scale_factor+tw_CR_station%offset
   enddo
   imold=im
 enddo
!for Surface Temperature I don t have varno stored in gs
!current_varno = 22 !this is the index for surface temperature
!status = read_20CR_S(iy,current_varno,tw_CR_station_s,time_dim)!in file
!             fields

!                il=(iy-1957)*12+im
!                nnk=nk/4/rcpara%pmax
!                fieldsm(:,:,:,il,1)=0_JPRB
!                do k=1,nnk
!                do i=1,rcpara%pmax
!                  fieldsm(:,:,i,il,1)=fieldsm(:,:,i,il,1)+fields(:,:,(k-1)*rcpara%pmax+i)
!                enddo
!                enddo
!                do i=1,rcpara%pmax
!                  fieldsm(:,:,i,il,1)=fieldsm(:,:,i,il,1)/nnk
!                enddo
                do istat=1,wmostats
                  call fill_tec(te20c,fields,wmolons,wmolats,istat,istart,iend,ni,iens)
                enddo 
             endif
end subroutine grib_and_fill_noaa

subroutine grib_and_fill(rcpara,te20c,e20path,fieldsm,wmolons,wmolats,wmostats,ni,nj,nk,iens,iy,im)

 use rwgrib2

implicit none

type(rasocor_namelist) rcpara
real(kind=jprm),intent(in) :: wmolons(:),wmolats(:)
real(kind=jprm) :: fieldsm(:,:,:,:,:)
type(cacherecord) :: te20c(:,:)
real:: fields(ni,nj,nk)
integer,intent(in) :: ni,nj,nk,iens,iy,im,wmostats
integer :: istat,iend,istart,mstride,mpar(1),i,k,il,nnk
character cdate*6,filename*80,cens*2
character*(*) :: e20path
logical:: ex3

           mstride=1
           write(cens,'(A1,I1.1)') '.',iens-1
           if (iens .eq. 1) then
             cens=''
           endif
           write(cdate,'(I4,I2.2)') iy,im
             filename=trim(e20path)//'/era20C.'//cdate//'.130'//trim(cens)
             write(*,*) filename
             inquire(file=filename,exist=ex3)
             if(ex3) then
                call READLATLON(FILENAME,fields,ni,nj,nk,mstride,mpar,omp_lp(1))

                il=(iy-rcpara%startdate/10000)*12+im
                nnk=nk/4/rcpara%pmax
                fieldsm(:,:,:,il,1)=0_JPRB
                do k=1,nnk
                do i=1,rcpara%pmax
                  fieldsm(:,:,i,il,1)=fieldsm(:,:,i,il,1)+fields(:,:,(k-1)*rcpara%pmax+i)
                enddo
                enddo
                do i=1,rcpara%pmax
                  fieldsm(:,:,i,il,1)=fieldsm(:,:,i,il,1)/nnk
                enddo
                 istart=toindex(iy*10000+im*100+1,rcpara)
                 if(im .lt. 12) then
                   iend=toindex(iy*10000+(im+1)*100+1,rcpara)-1
                 else
                   iend=toindex((iy+1)*10000+1*100+1,rcpara)-1
                 endif
                do istat=1,wmostats
                  call fill_tec(te20c,fields,wmolons,wmolats,istat,istart,iend,ni,iens)
                enddo 
             endif
end subroutine grib_and_fill

subroutine fill_tec(te20c,fields,wmolons,wmolats,istat,istart,iend,ni,iens)

implicit none

real(kind=jprm) :: lrd,uod
real(kind=jprm),intent(in) :: wmolons(:),wmolats(:)
type(cacherecord) :: te20c(:,:)
real,intent(in):: fields(:,:,:)
integer :: l,lr,uo,uop,lrp,stride0,stride1
integer,intent(in) :: istat,iend,istart,ni,iens
                  if(te20c(istat,1)%vals .gt. 0) then
                    do l=1,te20c(istat,1)%vals
                      if(te20c(istat,1)%index(l) .ge. istart .and. te20c(istat,1)%index(l) .le. iend) exit
                    enddo
                    if(l .le. te20c(istat,1)%vals ) then
                      if(wmolons(istat) .lt. 0) then
                         lr=floor((wmolons(istat)+360)/2.0)+1
                         lrd=((wmolons(istat)+360)/2.0-floor(wmolons(istat)+360)/2.0)
                         lrp=lr+1
                         if(lr == ni) lrp=1
                      else
                         lr=floor(wmolons(istat)/2.0)+1
                         lrd=((wmolons(istat))/2.0-floor(wmolons(istat)/2.0))
                         lrp=lr+1
                      endif
                      uo=floor((-wmolats(istat)+90)/2.0)+1
                      uod=((-wmolats(istat)+90)/2.0-floor((-wmolats(istat)+90)/2.0))
                      uop=uo+1
                      do while(te20c(istat,1)%index(l) .le. iend )
                      stride0=(te20c(istat,1)%index(l)-istart)*4*16+1
                      stride1=((te20c(istat,1)%index(l)-istart)*4+2)*16+1
                      te20c(istat,iens)%feld(l,:,1)=&
                        fields(lr,uo,stride0:stride0+15)*(1-lrd)*(1-uod)+&
                        fields(lrp,uo,stride0:stride0+15)*(lrd)*(1-uod)+&
                        fields(lr,uop,stride0:stride0+15)*(1-lrd)*(uod)+&
                        fields(lrp,uop,stride0:stride0+15)*(lrd)*(uod)
                      te20c(istat,iens)%feld(l,:,2)=&
                        fields(lr,uo,stride1:stride1+15)*(1-lrd)*(1-uod)+&
                        fields(lrp,uo,stride1:stride1+15)*(lrd)*(1-uod)+&
                        fields(lr,uop,stride1:stride1+15)*(1-lrd)*(uod)+&
                        fields(lrp,uop,stride1:stride1+15)*(lrd)*(uod)
                      l=l+1
                      if(l .gt. te20c(istat,1)%vals) exit
                      enddo
                    endif
                  endif
return
end subroutine fill_tec

end module rfmod
