########################
# for HP
include ../../input/fortran_compiler
F77       = $(FC)
FOPT      = 
AR = ar
########################
# for Hitachi Compiler on HP
#F77       = /usr/lib/hf_FORT/f77
#FOPT      = 
#AR = ar
########################
# for SunOS 5.5.1 Generic_103640-03 sun4u sparc SUNW,Ultra-Enterprise
#F77       = f77
#FOPT      = 
#AR = /usr/ccs/bin/ar
########################

MOBJ1      = main1.o
MOBJ2      = main2.o
QOBJ      = QEDPS.o

MOBJ1      = main1.o
MOBJ2      = main2.o
QOBJ      = QEDPS.o

QEDPS1     = qedps1
QEDPS2     = qedps2
LIBQEDPS  = libqedps.a

all:${QEDPS1} ${QEDPS2}

$(QEDPS1): $(MOBJ1) $(QEDPSLIB)
	$(F77) $(FOPT) -o $(QEDPS1) $(MOBJ1) -L./ -lqedps

$(QEDPS2): $(MOBJ2) $(QEDPSLIB)
	$(F77) $(FOPT) -o $(QEDPS2) $(MOBJ2) -L./ -lqedps
install:${LIBQEDPS}

$(LIBQEDPS): $(QOBJ)
	$(AR) rc $(LIBQEDPS) $(QOBJ)

.SUFFIXES: .f .o

.f.o:
	$(F77) $(FOPT) -c $*.f

clean:
	rm -f ${QEDPS1} $(QEDPS2)  ${LIBQEDPS} *.o junk* tmp* core 

cleaninstall:
	rm -f ${LIBQEDPS}

