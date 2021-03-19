SHELL=	/bin/ksh
FC      =$(FCMP)
LIBS    =$(LIBSM) -L/scratch2/BMC/gsienkf/Tseganeh.Gichamo/libs -llapack
OBJS	=num_parthds.o driver_snowOI.o snow_da.o snow_analysis.o

CMD     =snowOI
$(CMD):	$(OBJS)
	${FC} ${DEBUG} $(FFLAGS) $(LDFLG) $(OBJS) $(LIBS) -o $(CMD)

num_parthds.o:	num_parthds.f90
	$(FC) $(FFLAGS) $(OMPFLAG) -c num_parthds.f90

snow_da.o: snow_da.f90
	$(FC)  $(FFLAGS) $(OMPFLAG) -c snow_da.f90

snow_analysis.o: snow_da.o snow_analysis.f90
	$(FC)  $(FFLAGS) $(OMPFLAG) -c snow_analysis.f90

driver_snowOI.o: snow_analysis.o driver_snowOI.f90
	$(FC)  $(FFLAGS) -c driver_snowOI.f90

install:
	-cp -p $(CMD) ../exec/.
clean:
	-rm -f $(OBJS) *.mod $(CMD)
