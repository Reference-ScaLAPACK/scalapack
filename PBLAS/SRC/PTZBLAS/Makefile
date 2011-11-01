############################################################################
#
#  Program:         PBLAS -- (version 2.0)
#
#  Module:          Makefile
#
#  Purpose:         PTZBLAS source Makefile
#
#  Creation date:   April 1, 1998
#
#  Modified:
#
#  Send bug reports, comments or suggestions to scalapack@cs.utk.edu
#
############################################################################

include ../../../SLmake.inc

############################################################################
#
#  This is the makefile to create a library for the PBTZBLAS tools.
#  The files are grouped as follows:
#
#     ALLPBTZ -- Auxiliary routines for Level 1, 2 and 3 PBLAS
#
#       SPBTZ -- Single precision real    PBLAS Level F77 tools routines
#       CPBTZ -- Single precision complex PBLAS Level F77 tools routines
#       DPBTZ -- Double precision real    PBLAS Level F77 tools routines
#       ZPBTZ -- Double precision complex PBLAS Level F77 tools routines
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
#  The library is called
#       pblas_$(PLAT).a
#  (see ../../../SLmake.inc for the definition of $(PLAT)).
#
#  To remove the object files after the library is created, enter
#       make clean
#
############################################################################

all: single double complex complex16

APBTZ = \
    pxerbla.o

IPBTZ = \
    immadd.o       immdda.o       immtadd.o      immddat.o

SPBTZ = \
    svasum.o       sset.o         scshft.o       srshft.o     \
    svvdot.o       smmadd.o       smmcadd.o      smmtadd.o    \
    smmtcadd.o     smmdda.o       smmddac.o      smmddat.o    \
    smmddact.o     stzpad.o       stzpadcpy.o    stzscal.o    \
    sagemv.o       sasymv.o       satrmv.o

SCPBTZ = \
    sasqrtb.o      sascal.o

CPBTZ = \
    scvasum.o      cset.o         ccshft.o       crshft.o     \
    cvvdotu.o      cvvdotc.o      cmmadd.o       cmmcadd.o    \
    cmmtadd.o      cmmtcadd.o     cmmdda.o       cmmddac.o    \
    cmmddat.o      cmmddact.o     ctzpad.o       ctzpadcpy.o  \
    chescal.o      ctzscal.o      ctzcnjg.o      cagemv.o     \
    cahemv.o       catrmv.o       casymv.o       csymv.o      \
    csyr.o         csyr2.o

DPBTZ = \
    dvasum.o       dset.o         dcshft.o       drshft.o     \
    dvvdot.o       dmmadd.o       dmmcadd.o      dmmtadd.o    \
    dmmtcadd.o     dmmdda.o       dmmddac.o      dmmddat.o    \
    dmmddact.o     dtzpad.o       dtzpadcpy.o    dtzscal.o    \
    dagemv.o       dasymv.o       datrmv.o

DZPBTZ = \
    dasqrtb.o      dascal.o

ZPBTZ = \
    dzvasum.o      zset.o         zcshft.o       zrshft.o     \
    zvvdotu.o      zvvdotc.o      zmmadd.o       zmmcadd.o    \
    zmmtadd.o      zmmtcadd.o     zmmdda.o       zmmddac.o    \
    zmmddat.o      zmmddact.o     ztzpad.o       ztzpadcpy.o  \
    zhescal.o      ztzscal.o      ztzcnjg.o      zagemv.o     \
    zahemv.o       zatrmv.o       zasymv.o       zsymv.o      \
    zsyr.o         zsyr2.o

#---------------------------------------------------------------------------

single:    $(APBTZ) $(IPBTZ) $(SPBTZ) $(SCPBTZ)
	$(ARCH) $(ARCHFLAGS) ../../../$(SCALAPACKLIB) $(APBTZ) $(IPBTZ) $(SPBTZ) $(SCPBTZ)
	$(RANLIB) ../../../$(SCALAPACKLIB)

complex:   $(APBTZ) $(IPBTZ) $(CPBTZ) $(SCPBTZ)
	$(ARCH) $(ARCHFLAGS) ../../../$(SCALAPACKLIB) $(APBTZ) $(IPBTZ) $(CPBTZ) $(SCPBTZ)
	$(RANLIB) ../../../$(SCALAPACKLIB)

double:    $(APBTZ) $(IPBTZ) $(DPBTZ) $(DZPBTZ)
	$(ARCH) $(ARCHFLAGS) ../../../$(SCALAPACKLIB) $(APBTZ) $(IPBTZ) $(DPBTZ) $(DZPBTZ)
	$(RANLIB) ../../../$(SCALAPACKLIB)

complex16: $(APBTZ) $(IPBTZ) $(ZPBTZ) $(DZPBTZ)
	$(ARCH) $(ARCHFLAGS) ../../../$(SCALAPACKLIB) $(APBTZ) $(IPBTZ) $(ZPBTZ) $(DZPBTZ)
	$(RANLIB) ../../../$(SCALAPACKLIB)

#---------------------------------------------------------------------------

clean:
	rm -f *.o

.f.o : ; $(FC) -c $(FCFLAGS) $*.f
.c.o : ; $(CC) -c $(CDEFS) $(CCFLAGS) $*.c
