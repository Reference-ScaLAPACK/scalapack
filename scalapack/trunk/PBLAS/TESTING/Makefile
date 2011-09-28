############################################################################
#
#  Program:         PBLAS -- (version 2.0)
#
#  Module:          Makefile
#
#  Purpose:         Testing Makefile
#
#  Creation date:   April 1, 1998
#
#  Modified:        February 15, 2000
#
#  Send bug reports, comments or suggestions to scalapack@cs.utk.edu
#
############################################################################

include ../../SLmake.inc

sPBLAS1exe    = xspblas1tst
dPBLAS1exe    = xdpblas1tst
cPBLAS1exe    = xcpblas1tst
zPBLAS1exe    = xzpblas1tst

sPBLAS2exe    = xspblas2tst
dPBLAS2exe    = xdpblas2tst
cPBLAS2exe    = xcpblas2tst
zPBLAS2exe    = xzpblas2tst

sPBLAS3exe    = xspblas3tst
dPBLAS3exe    = xdpblas3tst
cPBLAS3exe    = xcpblas3tst
zPBLAS3exe    = xzpblas3tst

spb1tst       = $(PBLASTSTdir)/$(sPBLAS1exe)
dpb1tst       = $(PBLASTSTdir)/$(dPBLAS1exe)
cpb1tst       = $(PBLASTSTdir)/$(cPBLAS1exe)
zpb1tst       = $(PBLASTSTdir)/$(zPBLAS1exe)

spb2tst       = $(PBLASTSTdir)/$(sPBLAS2exe)
dpb2tst       = $(PBLASTSTdir)/$(dPBLAS2exe)
cpb2tst       = $(PBLASTSTdir)/$(cPBLAS2exe)
zpb2tst       = $(PBLASTSTdir)/$(zPBLAS2exe)

spb3tst       = $(PBLASTSTdir)/$(sPBLAS3exe)
dpb3tst       = $(PBLASTSTdir)/$(dPBLAS3exe)
cpb3tst       = $(PBLASTSTdir)/$(cPBLAS3exe)
zpb3tst       = $(PBLASTSTdir)/$(zPBLAS3exe)

pbtcom        = pblastst.o PB_Cwarn.o PB_Cabort.o

spbtcom       = psblastst.o slamch.o $(pbtcom)
dpbtcom       = pdblastst.o dlamch.o $(pbtcom)
cpbtcom       = pcblastst.o slamch.o $(pbtcom)
zpbtcom       = pzblastst.o dlamch.o $(pbtcom)

spb1t         = psblas1tst.o $(spbtcom)
dpb1t         = pdblas1tst.o $(dpbtcom)
cpb1t         = pcblas1tst.o $(cpbtcom)
zpb1t         = pzblas1tst.o $(zpbtcom)

spb2t         = psblas2tst.o $(spbtcom)
dpb2t         = pdblas2tst.o $(dpbtcom)
cpb2t         = pcblas2tst.o $(cpbtcom)
zpb2t         = pzblas2tst.o $(zpbtcom)

spb3t         = psblas3tst.o $(spbtcom)
dpb3t         = pdblas3tst.o $(dpbtcom)
cpb3t         = pcblas3tst.o $(cpbtcom)
zpb3t         = pzblas3tst.o $(zpbtcom)

all : single double complex complex16

single:    PblasErrorHandler $(spb1tst) $(spb2tst) $(spb3tst)

double:    PblasErrorHandler $(dpb1tst) $(dpb2tst) $(dpb3tst)

complex:   PblasErrorHandler $(cpb1tst) $(cpb2tst) $(cpb3tst)

complex16: PblasErrorHandler $(zpb1tst) $(zpb2tst) $(zpb3tst)

PblasErrorHandler:
	rm -f PB_Cwarn.c  PB_Cwarn.o
	ln -s $(PBLASdir)/SRC/PTOOLS/PB_Cwarn.c  PB_Cwarn.c
	$(CC) -c $(CCFLAGS)  -DTestingPblas PB_Cwarn.c
	rm -f PB_Cabort.c PB_Cabort.o
	ln -s $(PBLASdir)/SRC/PTOOLS/PB_Cabort.c PB_Cabort.c
	$(CC) -c $(CCFLAGS)  -DTestingPblas PB_Cabort.c

$(PBLASTSTdir)/PSBLAS1TST.dat: PSBLAS1TST.dat
	cp PSBLAS1TST.dat $(PBLASTSTdir)
$(PBLASTSTdir)/PDBLAS1TST.dat: PDBLAS1TST.dat
	cp PDBLAS1TST.dat $(PBLASTSTdir)
$(PBLASTSTdir)/PCBLAS1TST.dat: PCBLAS1TST.dat
	cp PCBLAS1TST.dat $(PBLASTSTdir)
$(PBLASTSTdir)/PZBLAS1TST.dat: PZBLAS1TST.dat
	cp PZBLAS1TST.dat $(PBLASTSTdir)

$(spb1tst) : $(SCALAPACKLIB) $(spb1t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(spb1tst) $(spb1t) $(LIBS)
	$(MAKE) $(PBLASTSTdir)/PSBLAS1TST.dat
$(dpb1tst) : $(SCALAPACKLIB) $(dpb1t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(dpb1tst) $(dpb1t) $(LIBS)
	$(MAKE) $(PBLASTSTdir)/PDBLAS1TST.dat
$(cpb1tst) : $(SCALAPACKLIB) $(cpb1t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(cpb1tst) $(cpb1t) $(LIBS)
	$(MAKE) $(PBLASTSTdir)/PCBLAS1TST.dat
$(zpb1tst) : $(SCALAPACKLIB) $(zpb1t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(zpb1tst) $(zpb1t) $(LIBS)
	$(MAKE) $(PBLASTSTdir)/PZBLAS1TST.dat

$(PBLASTSTdir)/PSBLAS2TST.dat: PSBLAS2TST.dat
	cp PSBLAS2TST.dat $(PBLASTSTdir)
$(PBLASTSTdir)/PDBLAS2TST.dat: PDBLAS2TST.dat
	cp PDBLAS2TST.dat $(PBLASTSTdir)
$(PBLASTSTdir)/PCBLAS2TST.dat: PCBLAS2TST.dat
	cp PCBLAS2TST.dat $(PBLASTSTdir)
$(PBLASTSTdir)/PZBLAS2TST.dat: PZBLAS2TST.dat
	cp PZBLAS2TST.dat $(PBLASTSTdir)

$(spb2tst) : $(SCALAPACKLIB) $(spb2t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(spb2tst) $(spb2t) $(LIBS)
	$(MAKE) $(PBLASTSTdir)/PSBLAS2TST.dat
$(dpb2tst) : $(SCALAPACKLIB) $(dpb2t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(dpb2tst) $(dpb2t) $(LIBS)
	$(MAKE) $(PBLASTSTdir)/PDBLAS2TST.dat
$(cpb2tst) : $(SCALAPACKLIB) $(cpb2t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(cpb2tst) $(cpb2t) $(LIBS)
	$(MAKE) $(PBLASTSTdir)/PCBLAS2TST.dat
$(zpb2tst) : $(SCALAPACKLIB) $(zpb2t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(zpb2tst) $(zpb2t) $(LIBS)
	$(MAKE) $(PBLASTSTdir)/PZBLAS2TST.dat

$(PBLASTSTdir)/PSBLAS3TST.dat: PSBLAS3TST.dat
	cp PSBLAS3TST.dat $(PBLASTSTdir)
$(PBLASTSTdir)/PDBLAS3TST.dat: PDBLAS3TST.dat
	cp PDBLAS3TST.dat $(PBLASTSTdir)
$(PBLASTSTdir)/PCBLAS3TST.dat: PCBLAS3TST.dat
	cp PCBLAS3TST.dat $(PBLASTSTdir)
$(PBLASTSTdir)/PZBLAS3TST.dat: PZBLAS3TST.dat
	cp PZBLAS3TST.dat $(PBLASTSTdir)

$(spb3tst) : $(SCALAPACKLIB) $(spb3t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(spb3tst) $(spb3t) $(LIBS)
	$(MAKE) $(PBLASTSTdir)/PSBLAS3TST.dat
$(dpb3tst) : $(SCALAPACKLIB) $(dpb3t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(dpb3tst) $(dpb3t) $(LIBS)
	$(MAKE) $(PBLASTSTdir)/PDBLAS3TST.dat
$(cpb3tst) : $(SCALAPACKLIB) $(cpb3t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(cpb3tst) $(cpb3t) $(LIBS)
	$(MAKE) $(PBLASTSTdir)/PCBLAS3TST.dat
$(zpb3tst) : $(SCALAPACKLIB) $(zpb3t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(zpb3tst) $(zpb3t) $(LIBS)
	$(MAKE) $(PBLASTSTdir)/PZBLAS3TST.dat

$(spb1t): $(FRC)
$(dpb1t): $(FRC)
$(cpb1t): $(FRC)
$(zpb1t): $(FRC)

$(spb2t): $(FRC)
$(dpb2t): $(FRC)
$(cpb2t): $(FRC)
$(zpb2t): $(FRC)

$(spb3t): $(FRC)
$(dpb3t): $(FRC)
$(cpb3t): $(FRC)
$(zpb3t): $(FRC)

FRC:
	@FRC=$(FRC)

clean :
	rm -f *.o

slamch.o:
	$(F77) -c $(NOOPT) slamch.f

dlamch.o:
	$(F77) -c $(NOOPT) dlamch.f

.f.o : ; $(F77) -c $(F77FLAGS) $*.f

# .c.o : ; $(CC) -c $(CCFLAGS) -DTestingPblas $*.c
