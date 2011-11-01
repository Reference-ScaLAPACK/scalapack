############################################################################
#
#  Program:         ScaLAPACK Redistrib
#
#  Module:          Makefile
#
#  Purpose:         Redistribution Sources Makefile
#
#  Creation date:   March 20, 1995
#
#  Modified:        February 15, 2000
#
#  Send bug reports, comments or suggestions to scalapack@cs.utk.edu
#
############################################################################

include ../../SLmake.inc

#######################################################################
#  This is the makefile to create a library for redistribution.
#  The files are organized as follows:
#       ALLAUX -- Auxiliary routines called from all precisions
#       IMRSRC -- Integer REDIST routines
#       SMRSRC -- Single precision real REDIST routines
#       CMRSRC -- Single precision complex REDIST routines
#       DMRSRC -- Double precision real REDIST routines
#       ZMRSRC -- Double precision complex REDIST routines
#
#  The library can be set up to include routines for any combination
#  of the four precisions.  First, modify the ARCH, ARCHFLAGS, RANLIB,
#  CC and CFLAGS definitions in ../../SLmake.inc to match your library
#  archiver, compiler and the options to be used.
#  Then to create or add to the library, enter make followed by one or
#  more of the precisions desired.  Some examples:
#       make single
#       make single complex
#       make single double complex complex16
#  Alternatively, the command
#       make
#  without any arguments creates a library of all four precisions.
#  The name of the library is defined by ../../$(SCALAPACKLIB) in
#  ../../SLmake.inc.
#
#  To remove the object files after the library is created, enter
#       make clean
#
#######################################################################

ALLAUX = pgemraux.o

IMRSRC = pigemr.o pigemr2.o pitrmr.o pitrmr2.o 

SMRSRC = psgemr.o psgemr2.o pstrmr.o pstrmr2.o 

CMRSRC = pcgemr.o pcgemr2.o pctrmr.o pctrmr2.o 

DMRSRC = pdgemr.o pdgemr2.o pdtrmr.o pdtrmr2.o 

ZMRSRC = pzgemr.o pzgemr2.o pztrmr.o pztrmr2.o 

all: integer single complex double complex16

lib: all

integer: $(IMRSRC) $(ALLAUX)
	$(ARCH) $(ARCHFLAGS) ../../$(SCALAPACKLIB) $(IMRSRC) $(ALLAUX)
	$(RANLIB) ../../$(SCALAPACKLIB)

single: $(SMRSRC) $(ALLAUX)
	$(ARCH) $(ARCHFLAGS) ../../$(SCALAPACKLIB) $(SMRSRC) $(ALLAUX)
	$(RANLIB) ../../$(SCALAPACKLIB)

complex: $(CMRSRC) $(ALLAUX)
	$(ARCH) $(ARCHFLAGS) ../../$(SCALAPACKLIB) $(CMRSRC) $(ALLAUX)
	$(RANLIB) ../../$(SCALAPACKLIB)

double: $(DMRSRC) $(ALLAUX)
	$(ARCH) $(ARCHFLAGS) ../../$(SCALAPACKLIB) $(DMRSRC) $(ALLAUX)
	$(RANLIB) ../../$(SCALAPACKLIB)

complex16: $(ZMRSRC) $(ALLAUX)
	$(ARCH) $(ARCHFLAGS) ../../$(SCALAPACKLIB) $(ZMRSRC) $(ALLAUX)
	$(RANLIB) ../../$(SCALAPACKLIB)

clean :
	rm -f *.o

.c.o : ; $(CC) -c $(CDEFS) $(CCFLAGS) $*.c
