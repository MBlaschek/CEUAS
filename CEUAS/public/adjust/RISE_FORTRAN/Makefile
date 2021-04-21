
#include ~/aix.mk

F90=$(F90C)
FC=$(F90C)

RTTOV= #RTTOV
RTTOVLIB=  -L$(HOME)/projects/rttov10/lib -lrttov10.1.0_test -lrttov10.1.0_other -lrttov10.1.0_coef_io  -lrttov10.1.0_parallel -lrttov10.1.0_main
RTTOVFLAGS= -I$(HOME)/projects/rttov10/include  -I$(HOME)/projects/rttov10/mod

.SUFFIXES: .o .f .F90 .f90 .F 
.f.o:
	$(F90) $(FFLAGS) -c $<

.f90.o:
	$(F90) $(FFLAGS) -c $<

.F90.o:
	$(F90) $(FFLAGS) -c $<

.f90.mod:
	$(F90) $(FFLAGS) -c $<
.mod.o:
	echo 

OPT	= -g -no-prec-div #xHost

GRIB=  -I/usr/local/ifort/include -L/usr/local/ifort/lib -Bstatic -lgrib_api_f90 -lgrib_api_f77  -lgrib_api -Bdynamic  -ljasper
GRIB=  -I/usr/local/ifort/grib1.12.3/include -L/usr/local/ifort/grib1.12.3/lib -Bstatic -lgrib_api_f90 -lgrib_api_f77  -lgrib_api -Bdynamic  -ljasper

F90C=ifort   $(OPT)  -r8  -132  -convert big_endian -traceback -nowarn -assume buffered_io -I.  -I$(NETCDF)/include $(GRIB)
FFLAGS =  $(OPT) $(IPO)  $(PROC) -qopenmp -cpp  -D$(RTTOV) -I/usr/local/include -I.  -I. $(RTTOVFLAGS) # -pg
F90FLAGS =  $(OPT)  $(IPO) $(PROC) -I/usr/local/include -I. $(RTTOVFLAGS)   #-pg
LDFLAGS =  $(OPT)  $(IPO) -qopenmp $(PROC) #-i-static #-pg
LNETCDF =  -L/usr/local/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lcurl

NAGLIB = #-lnag_r8



raso_correct_nc: dput.o  qsort.o grread.o biascor_era40new.o abor1.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o rwGRIB2.o error_management.o reading_writing_g_setting.o new_station.o rfmod.o calc_profile.o  homtests.o homtestsamp.o rfcorio.o  msu_fwd.o read_txt_write_nc.o rfcor.o  rfcomp_1.o correct_mostrecent.o rfcomp2.o correct_breaks2.o raso_correct.o raso_correct_igra_ei.o raso_correct_igra_NASA.o comp_correct.o rasocorrect_main.o rasocorrect_ini.o
	$(FC)  $(LDFLAGS) -o raso_correct_nc dput.o qsort.o abor1.o biascor_era40new.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o  error_management.o reading_writing_g_setting.o new_station.o rfmod.o rfcorio.o homtests.o homtestsamp.o rfcor.o calc_profile.o  rfcomp_1.o correct_mostrecent.o rfcomp2.o correct_breaks2.o raso_correct.o raso_correct_igra_ei.o raso_correct_igra_NASA.o comp_correct.o rasocorrect_main.o rasocorrect_ini.o read_txt_write_nc.o rwGRIB2.o grread.o  msu_fwd.o  $(LNETCDF) $(GRIB) -L/home/imgw/lib64 -lMvEmos_3.6 $(NAGLIB) $(RTTOVLIB)


