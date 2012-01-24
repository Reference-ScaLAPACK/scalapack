############################################################################
#
#  Program:         PBLAS -- (version 2.0)
#
#  Module:          Makefile
#
#  Purpose:         Timing Makefile
#
#  Creation date:   April 1, 1998
#
#  Modified:        February 15, 2000
#
#  Send bug reports, comments or suggestions to scalapack@cs.utk.edu
#
############################################################################

include ../../SLmake.inc

sPBLAS1exe    = xspblas1tim
dPBLAS1exe    = xdpblas1tim
cPBLAS1exe    = xcpblas1tim
zPBLAS1exe    = xzpblas1tim

sPBLAS2exe    = xspblas2tim
dPBLAS2exe    = xdpblas2tim
cPBLAS2exe    = xcpblas2tim
zPBLAS2exe    = xzpblas2tim

sPBLAS3exe    = xspblas3tim
dPBLAS3exe    = xdpblas3tim
cPBLAS3exe    = xcpblas3tim
zPBLAS3exe    = xzpblas3tim

spb1tim       = $(sPBLAS1exe)
dpb1tim       = $(dPBLAS1exe)
cpb1tim       = $(cPBLAS1exe)
zpb1tim       = $(zPBLAS1exe)

spb2tim       = $(sPBLAS2exe)
dpb2tim       = $(dPBLAS2exe)
cpb2tim       = $(cPBLAS2exe)
zpb2tim       = $(zPBLAS2exe)

spb3tim       = $(sPBLAS3exe)
dpb3tim       = $(dPBLAS3exe)
cpb3tim       = $(cPBLAS3exe)
zpb3tim       = $(zPBLAS3exe)

pbtcom        = pblastim.o PB_Cwarn.o PB_Cabort.o

spb1t         = psblas1tim.o psblastim.o $(pbtcom)
dpb1t         = pdblas1tim.o pdblastim.o $(pbtcom)
cpb1t         = pcblas1tim.o pcblastim.o $(pbtcom)
zpb1t         = pzblas1tim.o pzblastim.o $(pbtcom)

spb2t         = psblas2tim.o psblastim.o $(pbtcom)
dpb2t         = pdblas2tim.o pdblastim.o $(pbtcom)
cpb2t         = pcblas2tim.o pcblastim.o $(pbtcom)
zpb2t         = pzblas2tim.o pzblastim.o $(pbtcom)

spb3t         = psblas3tim.o psblastim.o $(pbtcom)
dpb3t         = pdblas3tim.o pdblastim.o $(pbtcom)
cpb3t         = pcblas3tim.o pcblastim.o $(pbtcom)
zpb3t         = pzblas3tim.o pzblastim.o $(pbtcom)

all : single double complex complex16

single:    PblasErrorHandler $(spb1tim) $(spb2tim) $(spb3tim)

double:    PblasErrorHandler $(dpb1tim) $(dpb2tim) $(dpb3tim)

complex:   PblasErrorHandler $(cpb1tim) $(cpb2tim) $(cpb3tim)

complex16: PblasErrorHandler $(zpb1tim) $(zpb2tim) $(zpb3tim)

PB_Cwarn.o:
	$(CC) -c $(CDEFS) $(CCFLAGS) -DTestingPblas -I../SRC -o PB_Cwarn.o ../SRC/PTOOLS/PB_Cwarn.c

PB_Cabort.o:
	$(CC) -c $(CDEFS) $(CCFLAGS) -DTestingPblas -I../SRC -o PB_Cabort.o ../SRC/PTOOLS/PB_Cabort.c

PblasErrorHandler: PB_Cwarn.o PB_Cabort.o

$(spb1tim) : ../../$(SCALAPACKLIB) $(spb1t)
	$(FCLOADER) $(FCLOADFLAGS) -o $(spb1tim) $(spb1t) ../../$(SCALAPACKLIB) $(LIBS)

$(dpb1tim) : ../../$(SCALAPACKLIB) $(dpb1t)
	$(FCLOADER) $(FCLOADFLAGS) -o $(dpb1tim) $(dpb1t) ../../$(SCALAPACKLIB) $(LIBS)

$(cpb1tim) : ../../$(SCALAPACKLIB) $(cpb1t)
	$(FCLOADER) $(FCLOADFLAGS) -o $(cpb1tim) $(cpb1t) ../../$(SCALAPACKLIB) $(LIBS)

$(zpb1tim) : ../../$(SCALAPACKLIB) $(zpb1t)
	$(FCLOADER) $(FCLOADFLAGS) -o $(zpb1tim) $(zpb1t) ../../$(SCALAPACKLIB) $(LIBS)

$(spb2tim) : ../../$(SCALAPACKLIB) $(spb2t)
	$(FCLOADER) $(FCLOADFLAGS) -o $(spb2tim) $(spb2t) ../../$(SCALAPACKLIB) $(LIBS)

$(dpb2tim) : ../../$(SCALAPACKLIB) $(dpb2t)
	$(FCLOADER) $(FCLOADFLAGS) -o $(dpb2tim) $(dpb2t) ../../$(SCALAPACKLIB) $(LIBS)

$(cpb2tim) : ../../$(SCALAPACKLIB) $(cpb2t)
	$(FCLOADER) $(FCLOADFLAGS) -o $(cpb2tim) $(cpb2t) ../../$(SCALAPACKLIB) $(LIBS)

$(zpb2tim) : ../../$(SCALAPACKLIB) $(zpb2t)
	$(FCLOADER) $(FCLOADFLAGS) -o $(zpb2tim) $(zpb2t) ../../$(SCALAPACKLIB) $(LIBS)

$(spb3tim) : ../../$(SCALAPACKLIB) $(spb3t)
	$(FCLOADER) $(FCLOADFLAGS) -o $(spb3tim) $(spb3t) ../../$(SCALAPACKLIB) $(LIBS)

$(dpb3tim) : ../../$(SCALAPACKLIB) $(dpb3t)
	$(FCLOADER) $(FCLOADFLAGS) -o $(dpb3tim) $(dpb3t) ../../$(SCALAPACKLIB) $(LIBS)

$(cpb3tim) : ../../$(SCALAPACKLIB) $(cpb3t)
	$(FCLOADER) $(FCLOADFLAGS) -o $(cpb3tim) $(cpb3t) ../../$(SCALAPACKLIB) $(LIBS)

$(zpb3tim) : ../../$(SCALAPACKLIB) $(zpb3t)
	$(FCLOADER) $(FCLOADFLAGS) -o $(zpb3tim) $(zpb3t) ../../$(SCALAPACKLIB) $(LIBS)

clean :
	rm -f *.o x*

.f.o : ; $(FC) -c $(FCFLAGS) $*.f

 .c.o : ; $(CC) -c $(CDEFS) $(CCFLAGS) -DTestingPblas $*.c
