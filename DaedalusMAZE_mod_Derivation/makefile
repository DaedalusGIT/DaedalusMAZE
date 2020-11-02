FCOMP = mpif90
FCFLAGS = -fopenmp
OPTFLAGS = -O2
PROGRAM =  JouleHeatingDerivation
SRCS = DeriveJH.F03 apex.f90
OBJECTS = $(SRCS:.F03 .F =.o)
EXTRAS =   -I/usr/include -L/usr/lib -lnetcdff -Wl,-Bsymbolic-functions -Wl,-z,relro -Wl,-z,now -lnetcdf -lnetcdf
LIBS =
DEBUG =


all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(FCOMP) $(FCFLAGS) $(DEBUG) -o $@ $^ $(LDFLAGS) 	$(LIBS)	$(EXTRAS) $(OPTFLAGS)

%.o: %.F02 %.f
	$(FCOMP) $(FCFLAGS) -c $<


.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAM) logfile.txt
