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

spb1tim       = $(PBLASTSTdir)/$(sPBLAS1exe)
dpb1tim       = $(PBLASTSTdir)/$(dPBLAS1exe)
cpb1tim       = $(PBLASTSTdir)/$(cPBLAS1exe)
zpb1tim       = $(PBLASTSTdir)/$(zPBLAS1exe)

spb2tim       = $(PBLASTSTdir)/$(sPBLAS2exe)
dpb2tim       = $(PBLASTSTdir)/$(dPBLAS2exe)
cpb2tim       = $(PBLASTSTdir)/$(cPBLAS2exe)
zpb2tim       = $(PBLASTSTdir)/$(zPBLAS2exe)

spb3tim       = $(PBLASTSTdir)/$(sPBLAS3exe)
dpb3tim       = $(PBLASTSTdir)/$(dPBLAS3exe)
cpb3tim       = $(PBLASTSTdir)/$(cPBLAS3exe)
zpb3tim       = $(PBLASTSTdir)/$(zPBLAS3exe)

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

PblasErrorHandler:
	rm -f PB_Cwarn.c  PB_Cwarn.o
	ln -s $(PBLASdir)/SRC/PTOOLS/PB_Cwarn.c  PB_Cwarn.c
	$(CC) -c $(CCFLAGS) -DTestingPblas PB_Cwarn.c
	rm -f PB_Cabort.c PB_Cabort.o
	ln -s $(PBLASdir)/SRC/PTOOLS/PB_Cabort.c PB_Cabort.c
	$(CC) -c $(CCFLAGS) -DTestingPblas PB_Cabort.c

$(PBLASTSTdir)/PSBLAS1TIM.dat: PSBLAS1TIM.dat
	cp PSBLAS1TIM.dat $(PBLASTSTdir)
$(PBLASTSTdir)/PDBLAS1TIM.dat: PDBLAS1TIM.dat
	cp PDBLAS1TIM.dat $(PBLASTSTdir)
$(PBLASTSTdir)/PCBLAS1TIM.dat: PCBLAS1TIM.dat
	cp PCBLAS1TIM.dat $(PBLASTSTdir)
$(PBLASTSTdir)/PZBLAS1TIM.dat: PZBLAS1TIM.dat
	cp PZBLAS1TIM.dat $(PBLASTSTdir)

$(spb1tim) : $(SCALAPACKLIB) $(spb1t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(spb1tim) $(spb1t) $(LIBS)
	$(MAKE) $(PBLASTSTdir)/PSBLAS1TIM.dat
$(dpb1tim) : $(SCALAPACKLIB) $(dpb1t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(dpb1tim) $(dpb1t) $(LIBS)
	$(MAKE) $(PBLASTSTdir)/PDBLAS1TIM.dat
$(cpb1tim) : $(SCALAPACKLIB) $(cpb1t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(cpb1tim) $(cpb1t) $(LIBS)
	$(MAKE) $(PBLASTSTdir)/PCBLAS1TIM.dat
$(zpb1tim) : $(SCALAPACKLIB) $(zpb1t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(zpb1tim) $(zpb1t) $(LIBS)
	$(MAKE) $(PBLASTSTdir)/PZBLAS1TIM.dat

$(PBLASTSTdir)/PSBLAS2TIM.dat: PSBLAS2TIM.dat
	cp PSBLAS2TIM.dat $(PBLASTSTdir)
$(PBLASTSTdir)/PDBLAS2TIM.dat: PDBLAS2TIM.dat
	cp PDBLAS2TIM.dat $(PBLASTSTdir)
$(PBLASTSTdir)/PCBLAS2TIM.dat: PCBLAS2TIM.dat
	cp PCBLAS2TIM.dat $(PBLASTSTdir)
$(PBLASTSTdir)/PZBLAS2TIM.dat: PZBLAS2TIM.dat
	cp PZBLAS2TIM.dat $(PBLASTSTdir)

$(spb2tim) : $(SCALAPACKLIB) $(spb2t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(spb2tim) $(spb2t) $(LIBS)
	$(MAKE) $(PBLASTSTdir)/PSBLAS2TIM.dat
$(dpb2tim) : $(SCALAPACKLIB) $(dpb2t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(dpb2tim) $(dpb2t) $(LIBS)
	$(MAKE) $(PBLASTSTdir)/PDBLAS2TIM.dat
$(cpb2tim) : $(SCALAPACKLIB) $(cpb2t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(cpb2tim) $(cpb2t) $(LIBS)
	$(MAKE) $(PBLASTSTdir)/PCBLAS2TIM.dat
$(zpb2tim) : $(SCALAPACKLIB) $(zpb2t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(zpb2tim) $(zpb2t) $(LIBS)
	$(MAKE) $(PBLASTSTdir)/PZBLAS2TIM.dat

$(PBLASTSTdir)/PSBLAS3TIM.dat: PSBLAS3TIM.dat
	cp PSBLAS3TIM.dat $(PBLASTSTdir)
$(PBLASTSTdir)/PDBLAS3TIM.dat: PDBLAS3TIM.dat
	cp PDBLAS3TIM.dat $(PBLASTSTdir)
$(PBLASTSTdir)/PCBLAS3TIM.dat: PCBLAS3TIM.dat
	cp PCBLAS3TIM.dat $(PBLASTSTdir)
$(PBLASTSTdir)/PZBLAS3TIM.dat: PZBLAS3TIM.dat
	cp PZBLAS3TIM.dat $(PBLASTSTdir)

$(spb3tim) : $(SCALAPACKLIB) $(spb3t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(spb3tim) $(spb3t) $(LIBS)
	$(MAKE) $(PBLASTSTdir)/PSBLAS3TIM.dat
$(dpb3tim) : $(SCALAPACKLIB) $(dpb3t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(dpb3tim) $(dpb3t) $(LIBS)
	$(MAKE) $(PBLASTSTdir)/PDBLAS3TIM.dat
$(cpb3tim) : $(SCALAPACKLIB) $(cpb3t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(cpb3tim) $(cpb3t) $(LIBS)
	$(MAKE) $(PBLASTSTdir)/PCBLAS3TIM.dat
$(zpb3tim) : $(SCALAPACKLIB) $(zpb3t)
	$(F77LOADER) $(F77LOADFLAGS) -o $(zpb3tim) $(zpb3t) $(LIBS)
	$(MAKE) $(PBLASTSTdir)/PZBLAS3TIM.dat

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

.f.o : ; $(F77) -c $(F77FLAGS) $*.f

# .c.o : ; $(CC) -c $(CCFLAGS) -DTestingPblas $*.c
