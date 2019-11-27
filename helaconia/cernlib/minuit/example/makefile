FFLAGS= -O
include ../../../input/fortran_compiler
#FC=gfortran
LIBDIR=../../../lib/

example1: example1.o
	$(FC) $(FFLAGS) -o example1 example1.o -L$(LIBDIR) -lminuit

example2: example2.o
	$(FC) $(FFLAGS) -o example2 example2.o -L$(LIBDIR) -lminuit

example3: example3.o
	$(FC) $(FFLAGS) -o example3 example3.o -L$(LIBDIR) -lminuit

.f.o:
	$(FC) $(FFLAGS) -c -o $@ $<
 
clean:
	rm -f *.o *~ *.out
	rm example1 example2 example3