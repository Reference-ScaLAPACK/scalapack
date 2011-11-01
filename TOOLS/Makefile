############################################################################
#
#  Program:         ScaLAPACK
#
#  Module:          Makefile
#
#  Purpose:         Tools Makefile
#
#  Creation date:   March 20, 1995
#
#  Modified:        February 15, 2000
#
#  Send bug reports, comments or suggestions to scalapack@cs.utk.edu
#
############################################################################

include ../SLmake.inc

ATOOLS = iceil.o      ilacpy.o     ilcm.o       indxg2p.o    indxg2l.o    \
         indxl2g.o    infog1l.o    infog2l.o    npreroc.o    numroc.o     \
         chk1mat.o    pchkxmat.o   sltimer.o    desc_convert.o \
         descinit.o   descset.o    reshape.o    SL_gridreshape.o \
         SL_init.o

ITOOLS = picol2row.o  pirow2col.o  pilaprnt.o   pitreecomb.o  pifillpad.o  \
         pichekpad.o  pielset.o    pielset2.o   pielget.o

STOOLS = dsnrm2.o     dsasum.o     slatcpy.o    ssdot.o       smatadd.o   \
         psmatadd.o   pscol2row.o  psrow2col.o  pslaprnt.o    pstreecomb.o \
         psfillpad.o  pschekpad.o  pselset.o    pselset2.o    pselget.o \
         pslaread.o   pslawrite.o

DTOOLS = dddot.o      dlatcpy.o    dmatadd.o    pdmatadd.o    pdcol2row.o \
         pdrow2col.o  pdlaprnt.o   pdtreecomb.o pdfillpad.o   pdchekpad.o \
         pdelset.o    pdelset2.o   pdelget.o \
         pdlaread.o   pdlawrite.o

CTOOLS = dscnrm2.o    dscasum.o    ccdotu.o     ccdotc.o      clatcpy.o \
         cmatadd.o    pcmatadd.o   pccol2row.o  pcrow2col.o   pclaprnt.o \
         pctreecomb.o pcfillpad.o  pcchekpad.o  pcelset.o     pcelset2.o \
         pcelget.o \
         pclaread.o   pclawrite.o

ZTOOLS = zzdotu.o     zzdotc.o     zlatcpy.o    zmatadd.o     pzmatadd.o \
         pzcol2row.o  pzrow2col.o  pzlaprnt.o   pztreecomb.o  pzfillpad.o \
         pzchekpad.o  pzelset.o    pzelset2.o   pzelget.o \
         pzlaread.o   pzlawrite.o

all: single double complex complex16

slapackaux:
	( cd LAPACK; $(MAKE) single )

dlapackaux:
	( cd LAPACK; $(MAKE) double )

clapackaux:
	( cd LAPACK; $(MAKE) complex )

zlapackaux:
	( cd LAPACK; $(MAKE) complex16 )

integer: $(ATOOLS) $(ITOOLS)
	$(ARCH) $(ARCHFLAGS) ../$(SCALAPACKLIB) $(ATOOLS) $(ITOOLS)
	$(RANLIB) ../$(SCALAPACKLIB)

single: slapackaux integer $(STOOLS)
	$(ARCH) $(ARCHFLAGS) ../$(SCALAPACKLIB) $(STOOLS)
	$(RANLIB) ../$(SCALAPACKLIB)

double: dlapackaux integer $(DTOOLS)
	$(ARCH) $(ARCHFLAGS) ../$(SCALAPACKLIB) $(DTOOLS)
	$(RANLIB) ../$(SCALAPACKLIB)

complex: clapackaux integer $(CTOOLS)
	$(ARCH) $(ARCHFLAGS) ../$(SCALAPACKLIB) $(CTOOLS)
	$(RANLIB) ../$(SCALAPACKLIB)

complex16: zlapackaux integer $(ZTOOLS)
	$(ARCH) $(ARCHFLAGS) ../$(SCALAPACKLIB) $(ZTOOLS)
	$(RANLIB) ../$(SCALAPACKLIB)

clean :
	rm -f *.o
	( cd LAPACK; $(MAKE) clean )

.f.o : ; $(FC) -c $(FCFLAGS) $*.f

.c.o : ; $(CC) -c $(CDEFS) $(CCFLAGS) $*.c
