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

spb1tst       = $(sPBLAS1exe)
dpb1tst       = $(dPBLAS1exe)
cpb1tst       = $(cPBLAS1exe)
zpb1tst       = $(zPBLAS1exe)

spb2tst       = $(sPBLAS2exe)
dpb2tst       = $(dPBLAS2exe)
cpb2tst       = $(cPBLAS2exe)
zpb2tst       = $(zPBLAS2exe)

spb3tst       = $(sPBLAS3exe)
dpb3tst       = $(dPBLAS3exe)
cpb3tst       = $(cPBLAS3exe)
zpb3tst       = $(zPBLAS3exe)

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

PB_Cwarn.o:
	$(CC) -c $(CDEFS) $(CCFLAGS) -I../SRC -o PB_Cwarn.o -DTestingPblas ../SRC/PTOOLS/PB_Cwarn.c

PB_Cabort.o:
	$(CC) -c $(CDEFS) $(CCFLAGS) -I../SRC -o PB_Cabort.o -DTestingPblas ../SRC/PTOOLS/PB_Cabort.c

PblasErrorHandler:  PB_Cwarn.o  PB_Cwarn.o

$(spb1tst) : ../../$(SCALAPACKLIB) $(spb1t)
	$(FCLOADER) $(FCLOADFLAGS) -o $(spb1tst) $(spb1t) ../../$(SCALAPACKLIB) $(LIBS)

$(dpb1tst) : ../../$(SCALAPACKLIB) $(dpb1t)
	$(FCLOADER) $(FCLOADFLAGS) -o $(dpb1tst) $(dpb1t) ../../$(SCALAPACKLIB) $(LIBS)

$(cpb1tst) : ../../$(SCALAPACKLIB) $(cpb1t)
	$(FCLOADER) $(FCLOADFLAGS) -o $(cpb1tst) $(cpb1t) ../../$(SCALAPACKLIB) $(LIBS)

$(zpb1tst) : ../../$(SCALAPACKLIB) $(zpb1t)
	$(FCLOADER) $(FCLOADFLAGS) -o $(zpb1tst) $(zpb1t) ../../$(SCALAPACKLIB) $(LIBS)

$(spb2tst) : ../../$(SCALAPACKLIB) $(spb2t)
	$(FCLOADER) $(FCLOADFLAGS) -o $(spb2tst) $(spb2t) ../../$(SCALAPACKLIB) $(LIBS)

$(dpb2tst) : ../../$(SCALAPACKLIB) $(dpb2t)
	$(FCLOADER) $(FCLOADFLAGS) -o $(dpb2tst) $(dpb2t) ../../$(SCALAPACKLIB) $(LIBS)

$(cpb2tst) : ../../$(SCALAPACKLIB) $(cpb2t)
	$(FCLOADER) $(FCLOADFLAGS) -o $(cpb2tst) $(cpb2t) ../../$(SCALAPACKLIB) $(LIBS)

$(zpb2tst) : ../../$(SCALAPACKLIB) $(zpb2t)
	$(FCLOADER) $(FCLOADFLAGS) -o $(zpb2tst) $(zpb2t) ../../$(SCALAPACKLIB) $(LIBS)

$(spb3tst) : ../../$(SCALAPACKLIB) $(spb3t)
	$(FCLOADER) $(FCLOADFLAGS) -o $(spb3tst) $(spb3t) ../../$(SCALAPACKLIB) $(LIBS)

$(dpb3tst) : ../../$(SCALAPACKLIB) $(dpb3t)
	$(FCLOADER) $(FCLOADFLAGS) -o $(dpb3tst) $(dpb3t) ../../$(SCALAPACKLIB) $(LIBS)

$(cpb3tst) : ../../$(SCALAPACKLIB) $(cpb3t)
	$(FCLOADER) $(FCLOADFLAGS) -o $(cpb3tst) $(cpb3t) ../../$(SCALAPACKLIB) $(LIBS)

$(zpb3tst) : ../../$(SCALAPACKLIB) $(zpb3t)
	$(FCLOADER) $(FCLOADFLAGS) -o $(zpb3tst) $(zpb3t) ../../$(SCALAPACKLIB) $(LIBS)

clean :
	rm -f *.o x*

slamch.o:
	$(FC) -c $(NOOPT) slamch.f

dlamch.o:
	$(FC) -c $(NOOPT) dlamch.f

.f.o : ; $(FC) -c $(FCFLAGS) $*.f

.c.o : ; $(CC) -c $(CDEFS) $(CCFLAGS) -DTestingPblas $*.c
