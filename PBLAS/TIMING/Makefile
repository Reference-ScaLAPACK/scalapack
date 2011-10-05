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
	$(CC) -c $(CCFLAGS) -DTestingPblas -I../SRC -o PB_Cwarn.o ../SRC/PTOOLS/PB_Cwarn.c

PB_Cabort.o:
	$(CC) -c $(CCFLAGS) -DTestingPblas -I../SRC -o PB_Cabort.o ../SRC/PTOOLS/PB_Cabort.c

PblasErrorHandler: PB_Cwarn.o PB_Cabort.o

$(spb1tim) : ../../$(SCALAPACKLIB) $(spb1t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(spb1tim) $(spb1t) ../../$(SCALAPACKLIB) $(LIBS)

$(dpb1tim) : ../../$(SCALAPACKLIB) $(dpb1t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(dpb1tim) $(dpb1t) ../../$(SCALAPACKLIB) $(LIBS)

$(cpb1tim) : ../../$(SCALAPACKLIB) $(cpb1t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(cpb1tim) $(cpb1t) ../../$(SCALAPACKLIB) $(LIBS)

$(zpb1tim) : ../../$(SCALAPACKLIB) $(zpb1t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(zpb1tim) $(zpb1t) ../../$(SCALAPACKLIB) $(LIBS)

$(spb2tim) : ../../$(SCALAPACKLIB) $(spb2t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(spb2tim) $(spb2t) ../../$(SCALAPACKLIB) $(LIBS)

$(dpb2tim) : ../../$(SCALAPACKLIB) $(dpb2t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(dpb2tim) $(dpb2t) ../../$(SCALAPACKLIB) $(LIBS)

$(cpb2tim) : ../../$(SCALAPACKLIB) $(cpb2t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(cpb2tim) $(cpb2t) ../../$(SCALAPACKLIB) $(LIBS)

$(zpb2tim) : ../../$(SCALAPACKLIB) $(zpb2t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(zpb2tim) $(zpb2t) ../../$(SCALAPACKLIB) $(LIBS)

$(spb3tim) : ../../$(SCALAPACKLIB) $(spb3t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(spb3tim) $(spb3t) ../../$(SCALAPACKLIB) $(LIBS)

$(dpb3tim) : ../../$(SCALAPACKLIB) $(dpb3t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(dpb3tim) $(dpb3t) ../../$(SCALAPACKLIB) $(LIBS)

$(cpb3tim) : ../../$(SCALAPACKLIB) $(cpb3t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(cpb3tim) $(cpb3t) ../../$(SCALAPACKLIB) $(LIBS)

$(zpb3tim) : ../../$(SCALAPACKLIB) $(zpb3t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(zpb3tim) $(zpb3t) ../../$(SCALAPACKLIB) $(LIBS)


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
	rm -f *.o x*

.f.o : ; $(F77) -c $(F77FLAGS) $*.f

# .c.o : ; $(CC) -c $(CCFLAGS) -DTestingPblas $*.c
