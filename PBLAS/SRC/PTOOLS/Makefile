############################################################################
#
#  Program:         PBLAS -- (version 2.0)
#
#  Module:          Makefile
#
#  Purpose:         PTOOLS makefile
#
#  Creation date:   April 1, 1998
#
#  Modified:        February 15, 2000
#
#  Send bug reports, comments or suggestions to scalapack@cs.utk.edu
#
############################################################################

include ../../../SLmake.inc

############################################################################
#
#  This is the makefile to create a library for the PBLAS F77 tools.
#  The files are grouped as follows:
#
#     ALLTOOLS -- Auxiliary routines for Level 1, 2 and 3 PBLAS
#
#       STOOLS -- Single precision real PBLAS Level F77 tools routines
#       CTOOLS -- Single precision complex PBLAS Level F77 tools routines
#       DTOOLS -- Double precision real PBLAS Level F77 tools routines
#       ZTOOLS -- Double precision complex PBLAS Level F77 tools routines
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

#---------------------------------------------------------------------------
#  The following line defines auxiliary routines needed by the PBLAS tools
#---------------------------------------------------------------------------

ALLCTOOLS = \
    PB_CGatherV.o    PB_CInV.o        PB_CInV2.o       PB_CInOutV.o     \
    PB_CInOutV2.o    PB_COutV.o       PB_CScatterV.o   PB_CVMinit.o     \
    PB_CVMloc.o      PB_CVMnpq.o      PB_CVMpack.o     PB_CVMswp.o      \
    PB_CVMupdate.o   PB_CVMcontig.o   PB_Cabort.o      PB_Cainfog2l.o   \
    PB_Cbinfo.o      PB_Cchkmat.o     PB_Cchkvec.o     PB_Cconjg.o      \
    PB_Cgetbuf.o     PB_Cinfog2l.o    PB_Citypeset.o   PB_Cgcd.o        \
    PB_Clcm.o        PB_Cmalloc.o     PB_Cnumroc.o     PB_Cg2lrem.o     \
    PB_Cindxg2p.o    PB_Cnnxtroc.o    PB_Cnpreroc.o    PB_CpswapNN.o    \
    PB_CpswapND.o    PB_Cpdot11.o     PB_CpdotNN.o     PB_CpdotND.o     \
    PB_CpaxpbyNN.o   PB_CpaxpbyND.o   PB_CpaxpbyDN.o   PB_Cpaxpby.o     \
    PB_CpgemmBC.o    PB_CpgemmAC.o    PB_CpgemmAB.o    PB_Cplaprnt.o    \
    PB_Cplapad.o     PB_Cplapd2.o     PB_Cplascal.o    PB_Cplasca2.o    \
    PB_Cplacnjg.o    PB_Cpsym.o       PB_CpsymmAB.o    PB_CpsymmBC.o    \
    PB_Cpsyr.o       PB_CpsyrkA.o     PB_CpsyrkAC.o    PB_Cpsyr2.o      \
    PB_Cpsyr2kA.o    PB_Cpsyr2kAC.o   PB_Cptrm.o       PB_Cpgeadd.o     \
    PB_Cptradd.o     PB_Cptran.o      PB_CptrmmAB.o    PB_CptrmmB.o     \
    PB_Cptrsm.o      PB_CptrsmAB.o    PB_CptrsmAB0.o   PB_CptrsmAB1.o   \
    PB_CptrsmB.o     PB_Cptrsv.o      PB_Ctop.o        PB_Ctzahemv.o    \
    PB_Ctzasymv.o    PB_Ctzatrmv.o    PB_Ctzhemm.o     PB_Ctzhemv.o     \
    PB_Ctzher.o      PB_Ctzherk.o     PB_Ctzher2.o     PB_Ctzher2k.o    \
    PB_Ctzsymm.o     PB_Ctzsymv.o     PB_Ctzsyr.o      PB_Ctzsyrk.o     \
    PB_Ctzsyr2.o     PB_Ctzsyr2k.o    PB_Ctztrmm.o     PB_Ctztrmv.o     \
    PB_Cwarn.o       PB_freebuf_.o    PB_topget_.o     PB_topset_.o     \
    PB_Cdescset.o    PB_Cdescribe.o   PB_CargFtoC.o    PB_Cfirstnb.o    \
    PB_Clastnb.o     PB_Cspan.o

#---------------------------------------------------------------------------

ALLTOOLS = $(ALLCTOOLS)

#---------------------------------------------------------------------------
#  The C PBLAS tools
#---------------------------------------------------------------------------

SCTOOLS = PB_Cstypeset.o

DCTOOLS = PB_Cdtypeset.o

CCTOOLS = PB_Cctypeset.o

ZCTOOLS = PB_Cztypeset.o

#---------------------------------------------------------------------------

STOOLS = $(SCTOOLS)
CTOOLS = $(CCTOOLS) $(SCTOOLS)
DTOOLS = $(DCTOOLS)
ZTOOLS = $(ZCTOOLS) $(DCTOOLS)

#---------------------------------------------------------------------------

SPTOOLS = $(ALLTOOLS) $(STOOLS)
CPTOOLS = $(ALLTOOLS) $(CTOOLS)
DPTOOLS = $(ALLTOOLS) $(DTOOLS)
ZPTOOLS = $(ALLTOOLS) $(ZTOOLS)

single: $(SPTOOLS)
	$(ARCH) $(ARCHFLAGS) ../../../$(SCALAPACKLIB) $(SPTOOLS)
	$(RANLIB) ../../../$(SCALAPACKLIB)

complex: $(CPTOOLS)
	$(ARCH) $(ARCHFLAGS) ../../../$(SCALAPACKLIB) $(CPTOOLS)
	$(RANLIB) ../../../$(SCALAPACKLIB)

double: $(DPTOOLS)
	$(ARCH) $(ARCHFLAGS) ../../../$(SCALAPACKLIB) $(DPTOOLS)
	$(RANLIB) ../../../$(SCALAPACKLIB)

complex16: $(ZPTOOLS)
	$(ARCH) $(ARCHFLAGS) ../../../$(SCALAPACKLIB) $(ZPTOOLS)
	$(RANLIB) ../../../$(SCALAPACKLIB)

clean:
	rm -f *.o

.f.o : ; $(FC) -c $(FCFLAGS) $*.f
.c.o : ; $(CC) -c $(CDEFS) $(CCFLAGS) $*.c
