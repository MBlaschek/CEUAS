
include ~/aix.mk

F90=$(F90C)

.F90.o:
	$(F90) $(FFLAGS) -c $<

.f90.mod:
	$(F90) $(FFLAGS) -c $<
.mod.o:
	echo 

#make 'OPT= -O3 -ipo -g -fp
OPT	= -g -O3 -no-prec-div #xHost 
PROC    = 
#LIB     = /home/imgw/lib64
#LIB     = /usr/local/lib/
GRIB=  -I/usr/local/ifort/include -L/usr/local/ifort/lib -Bstatic -lgrib_api_f90 -lgrib_api_f77  -lgrib_api -Bdynamic  -ljasper
GRIB=  -I/usr/local/ifort/grib1.12.3/include -L/usr/local/ifort/grib1.12.3/lib -Bstatic -lgrib_api_f90 -lgrib_api_f77  -lgrib_api -Bdynamic  -ljasper

FC=ifort  $(OPT) -r8  -132 -fopenmp  -convert big_endian -traceback  -nowarn -assume buffered_io -I. $(GRIB)
F90C=ifort   $(OPT)  -r8  -132 -fopenmp -convert big_endian -traceback -nowarn -assume buffered_io -I.  $(GRIB)
#FC=gfortran -openmp -fdefault-real-8  -ffree-line-length-none -ffixed-line-length-none -fconvert=big-endian 
#F90C=gfortran  -openmp  -fdefault-real-8  -ffree-line-length-none -ffixed-line-length-none -fconvert=big-endian
#FC=ifort -132 -convert big_endian -traceback  -nowarn -assume buffered_io
#F90C=ifort -132 -convert big_endian -traceback -nowarn -assume buffered_io
FFLAGS =  $(OPT) $(IPO)  $(PROC)  -I/usr/local/include -I.  -I. -I$(HOME)/projects/rttov10/include  -I$(HOME)/projects/rttov10/mod # -pg
F90FLAGS =  $(OPT)  $(IPO) $(PROC) -I/usr/local/include -I. -fopenmp  -I$(HOME)/projects/rttov10/include  -I$(HOME)/projects/rttov10/mod   #-pg
LDFLAGS =  $(OPT) -fopenmp $(IPO) $(PROC) #-i-static #-pg
FLIBSf  = -L$(LIB) -lfspectralf -L/home/imgw/software/NAG  -lnag_r8 -L/home/imgw/lib64 -L/home/imgw/software/emos_000310 -lemosR64A64_ifort
NAGLIB = -lnag_r8
#NETCDF =  -L$(LIB) -lnetcdff -lnetcdf  -Bdynamic -lhdf5_hl -lhdf5  -lz -lcurl
NETCDF =  -L/usr/local/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lcurl
#NETCDF =  /usr/local/lib/libnetcdff.a /usr/local/lib/libnetcdf.a  /usr/local/lib/libhdf5_hl.a /usr/local/lib/libhdf5.a  -Bdynamic -lz -lcurl
SEGCLUST= -L$(LIB) -Bstatic -lsegclust -Bdynamic  -lgsl -lgslcblas -lstdc++
RTTOV=  -L$(HOME)/projects/rttov10/lib -lrttov10.1.0_test -lrttov10.1.0_other -lrttov10.1.0_parallel -lrttov10.1.0_coef_io -lrttov10.1.0_main
#GRIB= -Bstatic -L$(HOME)/projects/grib_api-1.10.4_ifort/fortran/.libs -lgrib_api_f90 -lgrib_api_f77  -L$(HOME)/projects/grib_api-1.10.4_ifort/src/.libs -lgrib_api -Bdynamic  -ljasper
#GRIB=  -I/usr/local/ifort/include -L/usr/local/ifort -lgrib_api_f90 -lgrib_api_f77  -lgrib_api -lgrib_api_f90 -lgrib_api_f77  -lgrib_api -Bdynamic  -ljasper

$(LIB)/libfspectralf.a: grtophf.o phtogrf.o cfft.o gribinf.o grread.o cfftd.o
	ar rc $(LIB)/libfspectralf.a grtophf.o phtogrf.o cfft.o gribinf.o grread.o cfftd.o


$(BIN)/grtouf: grtouf.o $(LIB)/libfspectralf.a
	$(FC) $(LDFLAGS) -o $(BIN)/grtouf grtouf.o $(FLIBSf)

$(BIN)/fbdecode: fbdecode_1.0.o
	$(FC) $(LDFLAGS) -g -o $(BIN)/fbdecode fbdecode_1.0.o $(FLIBSf)


$(BIN)/fbdecodemike: fbdecode_mike.o
	$(FC) $(LDFLAGS) -g -o $(BIN)/fbdecodemike fbdecode_mike.o $(FLIBSf)

$(BIN)/fbdecode2: fbdecode_1.0.o
	$(FC) $(LDFLAGS) -g -o $(BIN)/fbdecode2 fbdecode_1.0.o $(FLIBSf)

$(BIN)/fbreduce: fbreduce.o
	$(FC) $(LDFLAGS) -g -o $(BIN)/fbreduce fbreduce.o $(FLIBSf)


$(BIN)/fbvanda: fbvanda.o
	$(FC) $(LDFLAGS) -g -o $(BIN)/fbvanda fbvanda.o $(FLIBSf)

$(BIN)/fbreduceguan: fbreduce_guan.o identify_station.o
	$(FC) $(LDFLAGS) -g -o $(BIN)/fbreduceguan identify_station.o fbreduce_guan.o $(FLIBSf)

$(BIN)/fbreduce_timeseries: fbreduce_timeseries.o
	$(FC) $(LDFLAGS) -g -o $(BIN)/fbreduce_timeseries identify_station.o fbreduce_timeseries.o $(FLIBSf)

$(BIN)/readfeed: biascor_era40.o diurnal.o pnterp.o rfmod.o dates.o identify_station.o readfeed.o
	$(FC) $(LDFLAGS) -g -o $(BIN)/readfeed readfeed.o rfmod.o biascor_era40.o identify_station.o diurnal.o pnterp.o dates.o $(FLIBSf)

$(BIN)/readfeeddaily: biascor_era40.o diurnal.o pnterp.o  identify_station.o rfmod.o dates.o readfeeddaily.o 
	$(FC) $(LDFLAGS) -g -o $(BIN)/readfeeddaily readfeeddaily.o identify_station.o rfmod.o dates.o biascor_era40.o diurnal.o pnterp.o $(FLIBSf)

$(BIN)/readfeeddailysave: biascor_era40.o diurnal.o pnterp.o  identify_station.o rfmod.o dates.o readfeeddailysave.o 
	$(FC) $(LDFLAGS) -g -o $(BIN)/readfeeddailysave readfeeddailysave.o identify_station.o rfmod.o dates.o biascor_era40.o diurnal.o pnterp.o $(FLIBSf)

$(BIN)/readfeeddailywind: biascor_era40.o diurnal.o pnterp.o  identify_station.o rfmod.o dates.o readfeeddailywind.o 
	$(FC) $(LDFLAGS) -g -o $(BIN)/readfeeddailywind readfeeddailywind.o identify_station.o rfmod.o dates.o biascor_era40.o diurnal.o pnterp.o $(FLIBSf)

	$(FC) $(LDFLAGS) -g -o $(BIN)/readfeeddailygeopot readfeeddailygeopot.o identify_station.o rfmod.o dates.o biascor_era40.o diurnal.o pnterp.o $(FLIBSf)

$(BIN)/cardsread: dates.o biascor_era40.o diurnal.o pnterp.o rfmod.o rfcards.o identify_station.o cards_read.o
	$(FC) $(LDFLAGS) -g -o $(BIN)/cardsread  rfmod.o rfcards.o dates.o biascor_era40.o diurnal.o pnterp.o identify_station.o cards_read.o $(FLIBSf)

$(BIN)/addfb_to_CARDS: dates.o biascor_era40.o diurnal.o pnterp.o rfmod.o rfcards.o identify_station.o gribinf90.o addfb_to_cards_2.o
	$(FC) $(LDFLAGS) -g -o $(BIN)/addfb_to_CARDS  rfmod.o rfcards.o dates.o biascor_era40.o diurnal.o pnterp.o identify_station.o gribinf90.o addfb_to_cards_2.o $(FLIBSf)

$(BIN)/extract_bg_merged: dates.o biascor_era40.o diurnal.o pnterp.o rfmod.o rfcards.o identify_station.o gribinf90.o extract_bg_merged.o
	$(FC) $(LDFLAGS) -g -o $(BIN)/extract_bg_merged  rfmod.o rfcards.o dates.o biascor_era40.o diurnal.o pnterp.o identify_station.o gribinf90.o extract_bg_merged.o $(FLIBSf)

$(BIN)/addfb_to_IGRA: dates.o biascor_era40.o diurnal.o pnterp.o rfmod.o rfcards.o identify_station.o gribinf90.o addfb_to_igra.o
	$(FC) $(LDFLAGS) -g -o $(BIN)/addfb_to_IGRA  rfmod.o rfcards.o dates.o biascor_era40.o diurnal.o pnterp.o identify_station.o gribinf90.o addfb_to_igra.o $(FLIBSf)

$(BIN)/analyse_fg: biascor_era40.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o rfmod.o gribinf90.o analyse_fg.o
	$(FC) $(LDFLAGS) -g -o $(BIN)/analyse_fg  biascor_era40.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o rfmod.o gribinf90.o analyse_fg.o $(FLIBSf)

$(BIN)/plotfeed:  rfmod.o statmod.o magroutines.o plotfeed.o 
	$(CC) $(LDFLAGS) -g -o $(BIN)/plotfeed plotfeed.o magroutines.o  rfmod.o statmod.o $(MAGLIB) $(EMOSLIB) $(ECLIB) $(NAGLIB)

$(BIN)/plottest: plottest.o
	$(CC) $(LDFLAGS) -g -o $(BIN)/plottest plottest.o $(MAGLIB) $(EMOSLIB) $(ECLIB)

$(BIN)/test_kalman: test_kalman.o
	$(FC) $(LDFLAGS) -g -o $(BIN)/test_kalman test_kalman.o $(NAGLIB)

$(BIN)/test_kalman_step: test_kalman_step.o
	$(FC) $(LDFLAGS) -g -o $(BIN)/test_kalman_step test_kalman_step.o $(NAGLIB)

$(BIN)/test_kalman_abc: test_kalman_abc.o
	$(FC) $(LDFLAGS) -g -o $(BIN)/test_kalman_abc test_kalman_abc.o $(NAGLIB)

$(BIN)/daily_kalman: daily_kalman.o
	$(FC) $(LDFLAGS) -g -o $(BIN)/daily_kalman daily_kalman.o $(NAGLIB)

$(BIN)/dailykalmancoupled: dailykalmancoupled.o
	$(FC) $(LDFLAGS) -g -o $(BIN)/dailykalmancoupled dailykalmancoupled.o $(NAGLIB)

grread.o: grread.f90
	$(FC) -c -g grread.f90

rfcards.o: rfcards.f
	$(FC) -c $(OPT) rfcards.f

biascor_era40new.o: biascor_era40new.f
	$(FC) -c $(OPT) -free biascor_era40new.f

abor1.o: abor1.f
	$(FC) -c $(OPT) -funderscoring -free abor1.f

homtests.o: homtests.f90 
	$(FC) -c -g  $(FFLAGS) $(OPT)  $(PROC) homtests.f90 

#rfcor.o: rfcor.f90
#	$(FC) -c -g -rcd -fp-model fast rfcor.f90


$(BIN)/comp_correctleo: biascor_era40new.o abor1.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o rfmod.o homtests.o homtestsamp.o rfcorio.o rfcor.o  calc_profile.o correct_breaks2.o rfcomp_1.o rfcomp2.o raso_correct.o raso_correct_igra.o  grread.o rasocorrect_main.o rasocorrect_ini.o
	$(FC) $(LDFLAGS) -o $(BIN)/comp_correctleo abor1.o biascor_era40new.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o rfmod.o rfcorio.o homtests.o homtestsamp.o rfcor.o calc_profile.o correct_breaks2.o rfcomp_1.o rfcomp2.o raso_correct.o raso_correct_igra.o rasocorrect_main.o rasocorrect_ini.o grread.o -L/home/imgw/lib64 -lMvEmos_3.6 $(NAGLIB)


$(BIN)/ulfcor: abor1.o biascor_era40new.o  rfcards.o identify_station.o dates.o diurnal.o pnterp.o rfmod.o  rfcorio.o grread.o read_ulfs.o ulfcor.o rasocorrect_ini.o
	$(FC) $(LDFLAGS) -o $(BIN)/ulfcor rfmod.o biascor_era40new.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o rfcorio.o  abor1.o grread.o read_ulfs.o ulfcor.o  rasocorrect_ini.o -L/home/imgw/lib64/pathscale -lMvEmos $(NAGLIB) # $(INTEL)

$(BIN)/ulfcorsplit: abor1.o biascor_era40new.o  rfcards.o identify_station.o dates.o diurnal.o pnterp.o rfmod.o  rfcorio.o grread.o read_ulfs_split.o ulfcor.o rasocorrect_ini.o
	$(FC) $(LDFLAGS) -o $(BIN)/ulfcorsplit rfmod.o biascor_era40new.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o rfcorio.o  abor1.o grread.o read_ulfs_split.o ulfcor.o  rasocorrect_ini.o -L/home/imgw/lib64/pathscale -lMvEmos $(NAGLIB) # $(INTEL)

#$(BIN)/raso_correct: biascor_era40new.o  grread.o abor1.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o rfmod.o homtests.o homtestsamp.o rfcorio.o rfcor.o  calc_profile.o correct_breaks2.o raso_correct.o raso_correct_igra.o raso_correct_igra_NASA.o rasocorrect_main.o rasocorrect_ini.o
#	$(FC) $(LDFLAGS) -o $(BIN)/raso_correct abor1.o biascor_era40new.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o rfmod.o rfcorio.o homtests.o homtestsamp.o rfcor.o calc_profile.o correct_breaks2.o raso_correct.o raso_correct_igra.o raso_correct_igra_NASA.o rasocorrect_main.o rasocorrect_ini.o grread.o -L/home/imgw/lib64/pathscale -lMvEmos $(NAGLIB) # $(INTEL)

$(BIN)/bg_correct: grread.o biascor_era40new.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o rfmod.o homtests.o homtestsamp.o rfcor.o rfcorio.o correct_breaks2.o abor1.o  raso_correct.o bgcorrect_main.o
	$(FC) $(LDFLAGS) -o $(BIN)/bg_correct biascor_era40new.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o rfmod.o homtests.o homtestsamp.o rfcor.o rfcorio.o correct_breaks2.o abor1.o  raso_correct.o bgcorrect_main.o grread.o -L$(LIB) $(NAGLIB) -L/home/imgw/metview/lib  -lMvEmos

$(BIN)/bg_correct_eqarea: grread.o biascor_era40new.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o rfmod.o homtests.o homtestsamp.o rfcor.o rfcorio.o correct_breaks2.o abor1.o calc_profile.o raso_correct.o bgcorrect_eqarea.o
	$(FC) $(LDFLAGS) -o $(BIN)/bg_correct_eqarea biascor_era40new.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o rfmod.o homtests.o homtestsamp.o rfcor.o rfcorio.o correct_breaks2.o abor1.o calc_profile.o  raso_correct.o bgcorrect_eqarea.o grread.o -L$(LIB) $(NAGLIB) -L/home/imgw/metview/lib  -lMvEmos

$(BIN)/bg_correct_hom: grread.o biascor_era40new.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o rfmod.o homtests.o homtestsamp.o rfcor.o rfcorio.o correct_breaks2.o abor1.o  raso_correct.o bgcorrect_main.o calc_profile.o
	$(FC) $(LDFLAGS) -g -o $(BIN)/bg_correct_hom biascor_era40new.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o rfmod.o homtests.o homtestsamp.o rfcor.o rfcorio.o correct_breaks2.o abor1.o  raso_correct.o bgcorrect_main.o grread.o calc_profile.o -L$(LIB) $(NAGLIB) -L/home/imgw/metview/lib  -lMvEmos

$(BIN)/test_bbs: rfcards.o identify_station.o dates.o diurnal.o pnterp.o rfmod.o homtests.o homtestsamp.o rfcor.o test_bbs.o 
	$(FC) $(LDFLAGS) -g -o $(BIN)/test_bbs biascor_era40.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o rfmod.o homtests.o homtestsamp.o rfcor.o test_bbs.o $(NAGLIB) 

$(BIN)/typetest: biascor_era40.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o rfmod.o typetest.o
	$(FC) $(LDFLAGS) -g -o $(BIN)/typetest typetest.o rfmod.o  biascor_era40.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o $(NAGLIB) 

$(BIN)/read_biascor: read_biascor.o
	$(FC) $(LDFLAGS) -g -o $(BIN)/read_biascor read_biascor.o

smoothgausslb.o:smoothgausslb.F
	$(FC) $(FFLAGS) -c -fixedform smoothgausslb.F
$(BIN)/smoothgausslb: smoothgausslb.o $(LIB)/libfspectralf.a utils.o IO2000.o
	$(FC) $(FFLAGS)  -o $(BIN)/smoothgausslb smoothgausslb.o  utils.o IO2000.o $(FLIBSf)

#$(BIN)/raso_correcti:  grread.o biascor_era40new.o abor1.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o rfmod.o homtests.o homtestsamp.o rfcorio.o rfcor.o calc_profile.o  rfcomp_1.o rfcomp2.o correct_breaks2.o raso_correct.o raso_correct_igra.o raso_correct_igra_NASA.o comp_correct.o rasocorrect_main.o rasocorrect_ini.o
#	$(FC) $(LDFLAGS) -o $(BIN)/raso_correcti abor1.o biascor_era40new.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o rfmod.o rfcorio.o homtests.o homtestsamp.o rfcor.o calc_profile.o  rfcomp_1.o rfcomp2.o correct_breaks2.o raso_correct.o raso_correct_igra.o raso_correct_igra_NASA.o comp_correct.o rasocorrect_main.o rasocorrect_ini.o grread.o -L/home/imgw/lib64 -lMvEmos_3.6 $(NAGLIB) 

raso_correct_nc: dput.o  qsort.o grread.o biascor_era40new.o abor1.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o rwGRIB2.o error_management.o reading_writing_g_setting.o new_station.o rfmod.o  homtests.o homtestsamp.o rfcorio.o  msu_fwd.o read_txt_write_nc.o rfcor.o calc_profile.o  rfcomp_1.o correct_mostrecent.o rfcomp2.o correct_breaks2.o raso_correct.o raso_correct_igra_ei.o raso_correct_igra_NASA.o comp_correct.o rasocorrect_main.o rasocorrect_ini.o
	$(FC) -openmp $(LDFLAGS) -o raso_correct_nc dput.o qsort.o abor1.o biascor_era40new.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o  error_management.o reading_writing_g_setting.o new_station.o rfmod.o rfcorio.o homtests.o homtestsamp.o rfcor.o calc_profile.o  rfcomp_1.o correct_mostrecent.o rfcomp2.o correct_breaks2.o raso_correct.o raso_correct_igra_ei.o raso_correct_igra_NASA.o comp_correct.o rasocorrect_main.o rasocorrect_ini.o read_txt_write_nc.o rwGRIB2.o grread.o  msu_fwd.o  $(NETCDF) $(GRIB) -L/home/imgw/lib64 -lMvEmos_3.6 $(NAGLIB) $(RTTOV)

read_era5_nc.o:   read_era5_nc.f90
	$(FC) $(FFLAGS) -c  read_era5_nc.f90

$(BIN)/read_era5_nc:  read_era5_nc.o
	$(FC) $(LDFLAGS) -o $(BIN)/read_era5_nc read_era5_nc.o -Bstatic -L/usr/local/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5  -Bdynamic -lz -lcurl


raso_correct_nc_test: dput.o  qsort.o grread.o biascor_era40new.o abor1.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o rwGRIB2.o error_management.o reading_writing_g_setting.o new_station.o rfmod.o  homtestsx.o homtestsamp.o rfcorio.o  msu_fwd.o read_txt_write_nc.o rfcor.o calc_profile.o  rfcomp_1.o correct_mostrecent.o rfcomp2.o correct_breaks2.o raso_correct.o raso_correct_lo.o raso_correct_igra_ei.o raso_correct_igra_NASA.o comp_correct.o rasocorrect_main.o rasocorrect_ini.o
	$(FC) $(LDFLAGS) -o raso_correct_nc_test dput.o qsort.o abor1.o biascor_era40new.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o  error_management.o reading_writing_g_setting.o new_station.o rfmod.o rfcorio.o homtestsx.o homtestsamp.o rfcor.o calc_profile.o  rfcomp_1.o correct_mostrecent.o rfcomp2.o correct_breaks2.o raso_correct.o raso_correct_lo.o raso_correct_igra_ei.o raso_correct_igra_NASA.o comp_correct.o rasocorrect_main.o rasocorrect_ini.o read_txt_write_nc.o rwGRIB2.o grread.o  msu_fwd.o  $(GRIB) -L/home/imgw/lib64 -lMvEmos_3.6 $(NAGLIB) $(RTTOV) $(NETCDF)

make_e20cfb: dput.o  qsort.o rwGRIB2.o biascor_era40new.o abor1.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o  new_station.o  error_management.o rfmod.o homtests.o homtestsamp.o rfcorio.o msu_fwd.o read_txt_write_nc.o rfcor.o  make_e20cfb_presat.o make_e20cfb.o make_noaa20cfb.o make_noaa20fb_presat.o make_e20cfb_ini.o
	$(FC) $(LDFLAGS) -o make_e20cfb dput.o qsort.o abor1.o biascor_era40new.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o   new_station.o error_management.o make_noaa20cfb.o make_e20cfb.o make_e20cfb_presat.o make_noaa20fb_presat.o make_e20cfb_ini.o  rfmod.o rfcorio.o read_txt_write_nc.o rwGRIB2.o  msu_fwd.o $(GRIB) -L/home/imgw/lib64 -lMvEmos_3.6 $(NAGLIB) $(RTTOV) $(NETCDF)

average_richens:  dates.o  qsort.o  biascor_era40new.o abor1.o diurnal.o pnterp.o rfmod.o rfcorio.o   msu_fwd.o  read_txt_write_nc.o  calc_profile.o  identify_station.o homtests.o homtestsamp.o rfcor.o average_richens.o
	$(FC) $(LDFLAGS) -o average_richens dates.o  qsort.o  biascor_era40new.o abor1.o diurnal.o pnterp.o rfmod.o rfcorio.o   msu_fwd.o read_txt_write_nc.o calc_profile.o  identify_station.o homtests.o homtestsamp.o rfcor.o  average_richens.o  -L/home/imgw/lib64 -lMvEmos_3.6 $(NAGLIB) $(RTTOV) $(NETCDF)


$(BIN)/odb_to_nc:  read_txt_write_nc.o msu_fwd.o odb_to_nc_ini.o abor1.o biascor_era40new.o rfcards.o  identify_station.o error_management.o reading_writing_g_setting.o new_station.o dates.o diurnal.o pnterp.o rfmod.o rwGRIB2.o
	$(FC) $(LDFLAGS) -o $(BIN)/odb_to_nc odb_to_nc_ini.o read_txt_write_nc.o msu_fwd.o rfmod.o abor1.o biascor_era40new.o rfcards.o identify_station.o error_management.o reading_writing_g_setting.o new_station.o dates.o diurnal.o pnterp.o rwGRIB2.o  $(GRIB) $(RTTOV) $(NETCDF)

$(BIN)/odb_to_nc_presat:  read_txt_write_nc.o msu_fwd.o odb_to_nc_ini_presat.o abor1.o biascor_era40new.o rfcards.o  identify_station.o error_management.o reading_writing_g_setting.o new_station.o dates.o diurnal.o pnterp.o rfmod.o rwGRIB2.o
	$(FC) $(LDFLAGS) -o $(BIN)/odb_to_nc_presat odb_to_nc_ini_presat.o read_txt_write_nc.o msu_fwd.o rfmod.o abor1.o biascor_era40new.o rfcards.o identify_station.o error_management.o reading_writing_g_setting.o new_station.o dates.o diurnal.o pnterp.o rwGRIB2.o  $(GRIB) $(RTTOV) $(NETCDF)

$(BIN)/odb_to_nc_oper:  read_txt_write_nc.o msu_fwd.o odb_to_nc_ini_oper.o abor1.o biascor_era40new.o rfcards.o identify_station.o error_management.o reading_writing_g_setting.o new_station.o dates.o diurnal.o pnterp.o rfmod.o rwGRIB2.o
	$(FC) $(LDFLAGS) -o $(BIN)/odb_to_nc_oper odb_to_nc_ini_oper.o read_txt_write_nc.o msu_fwd.o rfmod.o abor1.o biascor_era40new.o rfcards.o identify_station.o error_management.o reading_writing_g_setting.o new_station.o dates.o diurnal.o pnterp.o rwGRIB2.o  $(GRIB) $(RTTOV) $(NETCDF)

$(BIN)/odb_to_nc_MERRA:  read_txt_write_nc.o msu_fwd.o odb_to_nc_ini_MERRA.o abor1.o biascor_era40new.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o rfmod.o
	$(FC) $(LDFLAGS) -o $(BIN)/odb_to_nc_MERRA odb_to_nc_ini_MERRA.o read_txt_write_nc.o msu_fwd.o rfmod.o abor1.o biascor_era40new.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o  rwGRIB2.o   new_station.o  error_management.o reading_writing_g_setting.o $(GRIB) $(RTTOV) $(NETCDF)

$(BIN)/migrate_RICH:  read_txt_write_nc.o  rfcorio.o rfmod.o migrate_rich.o migrate_station.o migrate_RICH_ini.o abor1.o biascor_era40new.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o
	$(FC) $(LDFLAGS) -o $(BIN)/migrate_RICH migrate_rich.o migrate_station.o migrate_RICH_ini.o read_txt_write_nc.o rfcorio.o rfmod.o abor1.o biascor_era40new.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o $(NETCDF)

$(BIN)/richvis_old:  read_txt_write_nc.o  rfcorio.o rfmod.o ncto_richvis_json.o  richvis_ini.o abor1.o biascor_era40new.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o
	$(FC) $(LDFLAGS) -o $(BIN)/richvis_old ncto_richvis_json.o richvis_ini.o read_txt_write_nc.o rfcorio.o rfmod.o abor1.o biascor_era40new.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o $(NETCDF)

$(BIN)/richvis_leo:  read_txt_write_nc.o  rfcorio.o rfmod.o ncto_richvis_json.o  richvis_ini.o abor1.o biascor_era40new.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o
	$(FC) $(LDFLAGS) -o $(BIN)/richvis_old ncto_richvis_json.o richvis_ini.o read_txt_write_nc.o rfcorio.o rfmod.o abor1.o biascor_era40new.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o $(NETCDF)

$(BIN)/corrsave_to_nc:  read_txt_write_nc.o  rfcorio.o rfmod.o corrsave_to_nc.o abor1.o biascor_era40new.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o
	$(FC) $(LDFLAGS) -o $(BIN)/corrsave_to_nc  corrsave_to_nc.o read_txt_write_nc.o rfcorio.o rfmod.o abor1.o biascor_era40new.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o $(NETCDF)

$(BIN)/daily_to_monthly: error_management.o reading_writing_g_setting.o  new_station.o   rwGRIB2.o rfmod.o  qsort.o  msu_fwd.o rfcorio.o read_txt_write_nc.o abor1.o biascor_era40new.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o homtests.o homtestsamp.o rfcor.o  daily_to_monthly.o
	ifort $(LDFLAGS) -o $(BIN)/daily_to_monthly  daily_to_monthly.o error_management.o read_txt_write_nc.o rfcorio.o abor1.o biascor_era40new.o rfcards.o identify_station.o dates.o diurnal.o pnterp.o msu_fwd.o new_station.o qsort.o rwGRIB2.o rfmod.o $(GRIB) -Bstatic -L/usr/local/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5  -Bdynamic -lz -lcurl  $(RTTOV) -lopenjpeg






