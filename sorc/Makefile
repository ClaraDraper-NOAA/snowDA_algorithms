SHELL=	/bin/ksh
FC      =$(FCMP)
LIBS    =$(LIBSM) -L/scratch2/BMC/gsienkf/Tseganeh.Gichamo/libs -llapack
OBJS	=num_parthds.o snow_da.o ufs_land_update_snow.o scf_density.o snow_analysis.o driver_snowda.o

####

CMD     =driver_snowda
$(CMD):	$(OBJS)
	${FC} ${DEBUG} $(FFLAGS) $(LDFLG) $(OBJS) $(LIBS) -o $(CMD)

num_parthds.o:	num_parthds.f90
	$(FC) $(FFLAGS) $(OMPFLAG) -c num_parthds.f90

snow_da.o: snow_da.f90
	$(FC)  $(FFLAGS) -c snow_da.f90
	
ufs_land_update_snow.o: ufs_land_update_snow.f90
	 $(FC)  $(FFLAGS) -c ufs_land_update_snow.f90

scf_density.o: scf_density.f90
	 $(FC)  $(FFLAGS) -c scf_density.f90

snow_analysis.o: snow_da.o ufs_land_update_snow.o scf_density.o snow_analysis.f90
	$(FC)  $(FFLAGS) -c snow_analysis.f90

driver_snowda.o: snow_analysis.o driver_snowda.f90
	$(FC)  $(FFLAGS) -c driver_snowda.f90

install:
	-cp -p $(CMD) ../../exec/.
clean:
	-rm -f $(OBJS) *.mod $(CMD)
