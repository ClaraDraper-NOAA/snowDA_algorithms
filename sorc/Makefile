SHELL=	/bin/ksh
FC      =$(FCMP)
LIBS    =$(LIBSM) -L/scratch2/BMC/gsienkf/Tseganeh.Gichamo/libs -llapack
OBJS	=num_parthds.o driver_fIMS.o snow_da.o IMSaggregate_mod.o

CMD     =calcfIMS
$(CMD):	$(OBJS)
	${FC} ${DEBUG} $(FFLAGS) $(LDFLG) $(OBJS) $(LIBS) -o $(CMD)

num_parthds.o:	num_parthds.f90
	$(FC) $(FFLAGS) $(OMPFLAG) -c num_parthds.f90

snow_da.o: snow_da.f90
	$(FC)  $(FFLAGS) $(OMPFLAG) -c snow_da.f90

IMSaggregate_mod.o: snow_da.o IMSaggregate_mod.f90
	$(FC)  $(FFLAGS) $(OMPFLAG) -c IMSaggregate_mod.f90

driver_fIMS.o: IMSaggregate_mod.o driver_fIMS.f90
	$(FC)  $(FFLAGS) -c driver_fIMS.f90

install:
	-cp -p $(CMD) ../exec/.
clean:
	-rm -f $(OBJS) *.mod $(CMD)
