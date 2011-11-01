############################################################################
#
#  Program:         ScaLAPACK
#
#  Module:          Makefile
#
#  Purpose:         Redistribution Testing Makefile
#
#  Creation date:   March 20, 1995
#
#  Modified:        February 15, 2000
#
#  Send bug reports, comments or suggestions to scalapack@cs.utk.edu
#
############################################################################

include ../../SLmake.inc

igemrexe = xigemr
sgemrexe = xsgemr
dgemrexe = xdgemr
cgemrexe = xcgemr
zgemrexe = xzgemr
itrmrexe = xitrmr
strmrexe = xstrmr
dtrmrexe = xdtrmr
ctrmrexe = xctrmr
ztrmrexe = xztrmr

igemr = pigemrdrv.o
sgemr = psgemrdrv.o
dgemr = pdgemrdrv.o
cgemr = pcgemrdrv.o
zgemr = pzgemrdrv.o

itrmr = pitrmrdrv.o
strmr = pstrmrdrv.o
dtrmr = pdtrmrdrv.o
ctrmr = pctrmrdrv.o
ztrmr = pztrmrdrv.o

all: exe
exe: integer single double complex complex16

integer: $(igemrexe) $(itrmrexe)

single: $(sgemrexe) $(strmrexe) 

double: $(dgemrexe) $(dtrmrexe)

complex: $(cgemrexe) $(ctrmrexe)

complex16: $(zgemrexe) $(ztrmrexe)

$(igemrexe): ../../$(SCALAPACKLIB) $(igemr)
	$(CCLOADER) $(CCLOADFLAGS) -o $(igemrexe) $(igemr) ../../$(SCALAPACKLIB) $(LIBS)

$(sgemrexe): ../../$(SCALAPACKLIB) $(sgemr)
	$(CCLOADER) $(CCLOADFLAGS) -o $(sgemrexe) $(sgemr) ../../$(SCALAPACKLIB) $(LIBS)

$(dgemrexe): ../../$(SCALAPACKLIB) $(dgemr)
	$(CCLOADER) $(CCLOADFLAGS) -o $(dgemrexe) $(dgemr) ../../$(SCALAPACKLIB) $(LIBS)

$(cgemrexe): ../../$(SCALAPACKLIB) $(cgemr)
	$(CCLOADER) $(CCLOADFLAGS) -o $(cgemrexe) $(cgemr) ../../$(SCALAPACKLIB) $(LIBS)

$(zgemrexe): ../../$(SCALAPACKLIB) $(zgemr)
	$(CCLOADER) $(CCLOADFLAGS) -o $(zgemrexe) $(zgemr) ../../$(SCALAPACKLIB) $(LIBS)

$(itrmrexe): ../../$(SCALAPACKLIB) $(itrmr)
	$(CCLOADER) $(CCLOADFLAGS) -o $(itrmrexe) $(itrmr) ../../$(SCALAPACKLIB) $(LIBS)

$(strmrexe): ../../$(SCALAPACKLIB) $(strmr)
	$(CCLOADER) $(CCLOADFLAGS) -o $(strmrexe) $(strmr) ../../$(SCALAPACKLIB) $(LIBS)

$(dtrmrexe): ../../$(SCALAPACKLIB) $(dtrmr)
	$(CCLOADER) $(CCLOADFLAGS) -o $(dtrmrexe) $(dtrmr) ../../$(SCALAPACKLIB) $(LIBS)

$(ctrmrexe): ../../$(SCALAPACKLIB) $(ctrmr)
	$(CCLOADER) $(CCLOADFLAGS) -o $(ctrmrexe) $(ctrmr) ../../$(SCALAPACKLIB) $(LIBS)

$(ztrmrexe): ../../$(SCALAPACKLIB) $(ztrmr)
	$(CCLOADER) $(CCLOADFLAGS) -o $(ztrmrexe) $(ztrmr) ../../$(SCALAPACKLIB) $(LIBS)

clean :
	rm -f *.o x*

.c.o : ; $(CC) -c $(CDEFS) $(CCFLAGS) $*.c
