############################################################################
#
#  Program:         PBLAS -- (version 2.0)
#
#  Module:          Makefile
#
#  Purpose:         SRC Makefile
#
#  Creation date:   April 1, 1998
#
#  Modified:        February 15, 2000
#
#  Send bug reports, comments or suggestions to scalapack@cs.utk.edu
#
############################################################################

include ../../SLmake.inc

############################################################################
#
#  This is the makefile to create a library for the PBLAS.
#  The files are grouped as follows:
#
#      PSBLAS1 -- Single precision real    PBLAS1 routines
#      PCBLAS1 -- Single precision complex PBLAS1 routines
#      PDBLAS1 -- Double precision real    PBLAS1 routines
#      PZBLAS1 -- Double precision complex PBLAS1 routines
#
#      PSBLAS2 -- Single precision real    PBLAS2 routines
#      PCBLAS2 -- Single precision complex PBLAS2 routines
#      PDBLAS2 -- Double precision real    PBLAS2 routines
#      PZBLAS2 -- Double precision complex PBLAS2 routines
#
#      PSBLAS3 -- Single precision real    PBLAS3 routines
#      PCBLAS3 -- Single precision complex PBLAS3 routines
#      PDBLAS3 -- Double precision real    PBLAS3 routines
#      PZBLAS3 -- Double precision complex PBLAS3 routines
#
#  The library can be set up to include routines for any combination
#  of the four precisions.  First, modify the ../../SLmake.inc file
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

APPBLAS = pilaenv.o

#---------------------------------------------------------------------------
#  Comment out the next 5 definitions if you already have the Level 1 PBLAS.
#---------------------------------------------------------------------------

PIBLAS1 = picopy_.o

PSBLAS1 = psswap_.o psscal_.o  pscopy_.o  psaxpy_.o psdot_.o  psnrm2_.o \
          psasum_.o psamax_.o

PCBLAS1 = pcswap_.o pcscal_.o  pcsscal_.o pccopy_.o pcaxpy_.o pcdotu_.o \
          pcdotc_.o pscnrm2_.o pscasum_.o pcamax_.o

PDBLAS1 = pdswap_.o pdscal_.o  pdcopy_.o  pdaxpy_.o pddot_.o  pdnrm2_.o \
          pdasum_.o pdamax_.o

PZBLAS1 = pzswap_.o pzscal_.o  pzdscal_.o pzcopy_.o pzaxpy_.o pzdotu_.o \
          pzdotc_.o pdznrm2_.o pdzasum_.o pzamax_.o

#---------------------------------------------------------------------------
#  Comment out the next 4 definitions if you already have the Level 2 PBLAS.
#---------------------------------------------------------------------------

PSBLAS2 = psgemv_.o  psger_.o   pssymv_.o  pssyr_.o   pssyr2_.o  pstrmv_.o \
          pstrsv_.o  psagemv_.o psasymv_.o psatrmv_.o

PCBLAS2 = pcgemv_.o  pcgerc_.o  pcgeru_.o  pchemv_.o  pcher_.o   pcher2_.o \
          pctrmv_.o  pctrsv_.o  pcagemv_.o pcahemv_.o pcatrmv_.o

PDBLAS2 = pdgemv_.o  pdger_.o   pdsymv_.o  pdsyr_.o   pdsyr2_.o  pdtrmv_.o \
          pdtrsv_.o  pdagemv_.o pdasymv_.o pdatrmv_.o

PZBLAS2 = pzgemv_.o  pzgerc_.o  pzgeru_.o  pzhemv_.o  pzher_.o   pzher2_.o \
          pztrmv_.o  pztrsv_.o  pzagemv_.o pzahemv_.o pzatrmv_.o

#---------------------------------------------------------------------------
#  Comment out the next 4 definitions if you already have the Level 3 PBLAS.
#---------------------------------------------------------------------------

PSBLAS3 = psgeadd_.o psgemm_.o pssymm_.o  pssyr2k_.o pssyrk_.o  pstradd_.o \
          pstran_.o  pstrmm_.o pstrsm_.o

PCBLAS3 = pcgeadd_.o pcgemm_.o pchemm_.o  pcher2k_.o pcherk_.o  pcsymm_.o  \
          pcsyr2k_.o pcsyrk_.o pctradd_.o pctranc_.o pctranu_.o pctrmm_.o  \
          pctrsm_.o

PDBLAS3 = pdgeadd_.o pdgemm_.o pdsymm_.o  pdsyr2k_.o pdsyrk_.o  pdtradd_.o \
          pdtran_.o  pdtrmm_.o pdtrsm_.o

PZBLAS3 = pzgeadd_.o pzgemm_.o pzhemm_.o  pzher2k_.o pzherk_.o  pzsymm_.o  \
          pzsyr2k_.o pzsyrk_.o pztradd_.o pztranc_.o pztranu_.o pztrmm_.o  \
          pztrsm_.o

#---------------------------------------------------------------------------

PIBLAS = $(PIBLAS1) $(APPBLAS)
PSBLAS = $(PSBLAS1) $(PSBLAS2) $(PSBLAS3) $(APPBLAS)
PCBLAS = $(PCBLAS1) $(PCBLAS2) $(PCBLAS3) $(APPBLAS)
PDBLAS = $(PDBLAS1) $(PDBLAS2) $(PDBLAS3) $(APPBLAS)
PZBLAS = $(PZBLAS1) $(PZBLAS2) $(PZBLAS3) $(APPBLAS)

integer: $(PIBLAS)
	$(ARCH) $(ARCHFLAGS) ../../$(SCALAPACKLIB) $(PIBLAS)
	$(RANLIB) ../../$(SCALAPACKLIB)

single: integer $(PSBLAS)
	( cd PBBLAS;  $(MAKE) single )
	( cd PTZBLAS; $(MAKE) single )
	( cd PTOOLS;  $(MAKE) single )
	$(ARCH) $(ARCHFLAGS) ../../$(SCALAPACKLIB) $(PSBLAS)
	$(RANLIB) ../../$(SCALAPACKLIB)

double: integer $(PDBLAS)
	( cd PBBLAS;  $(MAKE) double )
	( cd PTZBLAS; $(MAKE) double )
	( cd PTOOLS;  $(MAKE) double )
	$(ARCH) $(ARCHFLAGS) ../../$(SCALAPACKLIB) $(PDBLAS)
	$(RANLIB) ../../$(SCALAPACKLIB)

complex: integer $(PCBLAS)
	( cd PBBLAS;  $(MAKE) complex )
	( cd PTZBLAS; $(MAKE) complex )
	( cd PTOOLS;  $(MAKE) complex )
	$(ARCH) $(ARCHFLAGS) ../../$(SCALAPACKLIB) $(PCBLAS)
	$(RANLIB) ../../$(SCALAPACKLIB)

complex16: integer $(PZBLAS)
	( cd PBBLAS;  $(MAKE) complex16 )
	( cd PTZBLAS; $(MAKE) complex16 )
	( cd PTOOLS;  $(MAKE) complex16 )
	$(ARCH) $(ARCHFLAGS) ../../$(SCALAPACKLIB) $(PZBLAS)
	$(RANLIB) ../../$(SCALAPACKLIB)

clean :
	rm -f *.o
	( cd PBBLAS;  $(MAKE) clean )
	( cd PTZBLAS; $(MAKE) clean )
	( cd PTOOLS;  $(MAKE) clean )

.f.o : ; $(FC) -c $(FCFLAGS) $*.f

.c.o : ; $(CC) -c $(CDEFS) $(CCFLAGS) $*.c
#
# To compile  without  input  argument checking replace previous line by
# .c.o : ; $(CC) -c $(CCFLAGS)  -DNO_ARGCHK $*.c
#
# Note that when the PBLAS routines have been compiled without input ar-
# gument checking, it is impossible to test PBLAS error exits. Trying to
# test PBLAS error exits will then cause the tester executable to hang.
#
