F95=ftn
FFLAGS=-cpp -O3 -I/N/u/nbcallah/BigRed2/libraries/include/
VPATH=modules
MODOBJ=constants.o forcesAndPotential.o symplecticInt.o testSubroutines.o #lyapunov.o trackGeometry.o
COBJ=modules/fields_fortran.o
LFLAGS=-L/N/u/nbcallah/BigRed2/libraries/lib -lnlopt

#all: symplecticInt trajGen henonHeiles symplecticIntNoPert symplecticIntNoPertEscape #symplecticIntPertEscapeFieldStrength biasedVsqDescent

#symplecticInt: $(MODOBJ) symplecticStep-currentRings-MPI.o
#	$(F95) $(MODOBJ) symplecticStep-currentRings-MPI.o -o symplecticInt $(LFLAGS)

#symplecticIntNoPert: $(MODOBJ) symplecticStep-NoPert.o
#	$(F95) $(MODOBJ) symplecticStep-NoPert.o -o symplecticIntNoPert $(LFLAGS)
	
#symplecticIntNoPertEscape: $(MODOBJ) symplecticStep-NoPert-Escape.o
#	$(F95) $(MODOBJ) symplecticStep-NoPert-Escape.o -o symplecticIntNoPertEscape $(LFLAGS)
	
#symplecticIntPertEscapeFieldStrength: $(MODOBJ) symplecticStep-WindowPert-Escape-FieldStrengthScan.o
#	$(F95) $(MODOBJ) symplecticStep-WindowPert-Escape-FieldStrengthScan.o -o symplecticIntPertEscapeFieldStrength $(LFLAGS)

#trajGen: $(MODOBJ) symplecticTrajectory.o
#	$(F95) $(MODOBJ) symplecticTrajectory.o -o trajGen $(LFLAGS)

all: test_C_eval

test_C_eval: $(MODOBJ) test_C_eval.o
	$(F95) $(MODOBJ) $(COBJ) test_C_eval.o -o test_C_eval $(LFLAGS)

%.o:%.f95 modules/constants.h
	$(F95) $(FFLAGS) -c $< -o $@
	
clean:
	find . -type f -name '*.o' -delete
	find . -type f -name '*.mod' -delete
	rm symplecticInt