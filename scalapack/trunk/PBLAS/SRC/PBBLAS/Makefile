############################################################################
#
#  Program:         PBLAS -- (version 2.0 beta)
#
#  Module:          Makefile
#
#  Purpose:         PB-BLAS remaining source Makefile
#
#  Creation date:   August 24, 1997
#
#  Modified:        February 15, 2000
#
#  Send bug reports, comments or suggestions to scalapack@cs.utk.edu
#
############################################################################

include ../../../SLmake.inc

############################################################################
#
#  The library can be set up to include routines for any combination
#  of the four precisions.  First, modify the ../../../SLmake.inc file
#  definitions to match your compiler and the options to be used.
#  Then to create or add to the library, enter make followed by one or
#  more of the precisions desired.  Some examples:
#       make single
#       make single complex
#       make single double complex complex16
#  Alternatively, the command
#       make
#  without any arguments creates a library of all four precisions.
#
#  To remove the object files after the library is created, enter
#       make clean
#
############################################################################

all: single double complex complex16

PBSBLASAUX = pbstran.o pbsmatadd.o pbstrsrt.o pbstrget.o \
             pbstrnv.o pbsvecadd.o pbstrst1.o

PBCBLASAUX = pbctran.o pbcmatadd.o pbctrsrt.o pbctrget.o \
             pbctrnv.o pbcvecadd.o pbctrst1.o

PBDBLASAUX = pbdtran.o pbdmatadd.o pbdtrsrt.o pbdtrget.o \
             pbdtrnv.o pbdvecadd.o pbdtrst1.o

PBZBLASAUX = pbztran.o pbzmatadd.o pbztrsrt.o pbztrget.o \
             pbztrnv.o pbzvecadd.o pbztrst1.o

#---------------------------------------------------------------------------

single: $(PBSBLASAUX)
	$(ARCH) $(ARCHFLAGS) ../../../$(SCALAPACKLIB) $(PBSBLASAUX)
	$(RANLIB) ../../../$(SCALAPACKLIB)

double: $(PBDBLASAUX)
	$(ARCH) $(ARCHFLAGS) ../../../$(SCALAPACKLIB) $(PBDBLASAUX)
	$(RANLIB) ../../../$(SCALAPACKLIB)

complex: $(PBCBLASAUX)
	$(ARCH) $(ARCHFLAGS) ../../../$(SCALAPACKLIB) $(PBCBLASAUX)
	$(RANLIB) ../../../$(SCALAPACKLIB)

complex16: $(PBZBLASAUX)
	$(ARCH) $(ARCHFLAGS) ../../../$(SCALAPACKLIB) $(PBZBLASAUX)
	$(RANLIB) ../../../$(SCALAPACKLIB)

#---------------------------------------------------------------------------

clean :
	rm -f *.o

.f.o : ; $(FC) -c $(FCFLAGS) $*.f
