############################################################################
#
#  Program:         ScaLAPACK
#
#  Module:          Makefile
#
#  Purpose:         Tools - LAPACK Makefile
#
#  Creation date:   March 20, 1995
#
#  Modified:        February 16, 2000
#
#  Send bug reports, comments or suggestions to scalapack@cs.utk.edu
#
############################################################################

include ../../SLmake.inc

ALLAUX = icopy.o 

SCATGEN = slatm1.o slaran.o slarnd.o

SMATGEN = slatms.o slagge.o slagsy.o slarot.o

CMATGEN = clarnv.o clatm1.o clatms.o clagge.o claghe.o clagsy.o clarot.o \
          clarnd.o

DZATGEN = dlatm1.o dlaran.o dlarnd.o

DMATGEN = dlatms.o dlagge.o dlagsy.o dlarot.o

ZMATGEN = zlarnv.o zlatm1.o zlatms.o zlagge.o zlaghe.o zlagsy.o zlarot.o \
          zlarnd.o

all : single complex double complex16

single: $(ALLAUX) $(SMATGEN) $(SCATGEN)
	$(ARCH) $(ARCHFLAGS) ../../$(SCALAPACKLIB) $(SLASRC) $(ALLAUX) $(SCLAUX) \
	$(SMATGEN) $(SCATGEN)
	$(RANLIB) ../../$(SCALAPACKLIB)

complex: $(ALLAUX) $(CMATGEN) $(SCATGEN)
	$(ARCH) $(ARCHFLAGS) ../../$(SCALAPACKLIB) $(CLASRC) $(ALLAUX) $(SCLAUX) \
	$(CMATGEN) $(SCATGEN)
	$(RANLIB) ../../$(SCALAPACKLIB)

double: $(ALLAUX) $(DMATGEN) $(DZATGEN)
	$(ARCH) $(ARCHFLAGS) ../../$(SCALAPACKLIB) $(DLASRC) $(ALLAUX) $(DZLAUX) \
	$(DMATGEN) $(DZATGEN)
	$(RANLIB) ../../$(SCALAPACKLIB)

complex16: $(ALLAUX) $(ZMATGEN) $(DZATGEN)
	$(ARCH) $(ARCHFLAGS) ../../$(SCALAPACKLIB) $(ZLASRC) $(ALLAUX) $(DZLAUX) \
	$(ZMATGEN) $(DZATGEN)
	$(RANLIB) ../../$(SCALAPACKLIB)

clean :
	rm -f *.o

slamch.o:
	$(FC) -c $(NOOPT) slamch.f

dlamch.o:
	$(FC) -c $(NOOPT) dlamch.f

.f.o : ; $(FC) -c $(FCFLAGS) $*.f

