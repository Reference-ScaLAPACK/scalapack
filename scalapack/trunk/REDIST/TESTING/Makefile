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

igemrexe = $(REDISTTSTdir)/xigemr
sgemrexe = $(REDISTTSTdir)/xsgemr
dgemrexe = $(REDISTTSTdir)/xdgemr
cgemrexe = $(REDISTTSTdir)/xcgemr
zgemrexe = $(REDISTTSTdir)/xzgemr
itrmrexe = $(REDISTTSTdir)/xitrmr
strmrexe = $(REDISTTSTdir)/xstrmr
dtrmrexe = $(REDISTTSTdir)/xdtrmr
ctrmrexe = $(REDISTTSTdir)/xctrmr
ztrmrexe = $(REDISTTSTdir)/xztrmr

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

$(TESTINGdir)/GEMR2D.dat: GEMR2D.dat
	cp GEMR2D.dat $(TESTINGdir)

$(igemrexe): $(SCALAPACKLIB) $(igemr)
	$(CCLOADER) $(CCLOADFLAGS) -o $(igemrexe) $(igemr) $(PRLIBS)
	$(MAKE) $(TESTINGdir)/GEMR2D.dat

$(sgemrexe): $(SCALAPACKLIB) $(sgemr)
	$(CCLOADER) $(CCLOADFLAGS) -o $(sgemrexe) $(sgemr) $(PRLIBS)
	$(MAKE) $(TESTINGdir)/GEMR2D.dat

$(dgemrexe): $(SCALAPACKLIB) $(dgemr)
	$(CCLOADER) $(CCLOADFLAGS) -o $(dgemrexe) $(dgemr) $(PRLIBS)
	$(MAKE) $(TESTINGdir)/GEMR2D.dat

$(cgemrexe): $(SCALAPACKLIB) $(cgemr)
	$(CCLOADER) $(CCLOADFLAGS) -o $(cgemrexe) $(cgemr) $(PRLIBS)
	$(MAKE) $(TESTINGdir)/GEMR2D.dat

$(zgemrexe): $(SCALAPACKLIB) $(zgemr)
	$(CCLOADER) $(CCLOADFLAGS) -o $(zgemrexe) $(zgemr) $(PRLIBS)
	$(MAKE) $(TESTINGdir)/GEMR2D.dat

$(TESTINGdir)/TRMR2D.dat: TRMR2D.dat
	cp TRMR2D.dat $(TESTINGdir)

$(itrmrexe): $(SCALAPACKLIB) $(itrmr)
	$(CCLOADER) $(CCLOADFLAGS) -o $(itrmrexe) $(itrmr) $(PRLIBS)
	$(MAKE) $(TESTINGdir)/TRMR2D.dat

$(strmrexe): $(SCALAPACKLIB) $(strmr)
	$(CCLOADER) $(CCLOADFLAGS) -o $(strmrexe) $(strmr) $(PRLIBS)
	$(MAKE) $(TESTINGdir)/TRMR2D.dat

$(dtrmrexe): $(SCALAPACKLIB) $(dtrmr)
	$(CCLOADER) $(CCLOADFLAGS) -o $(dtrmrexe) $(dtrmr) $(PRLIBS)
	$(MAKE) $(TESTINGdir)/TRMR2D.dat

$(ctrmrexe): $(SCALAPACKLIB) $(ctrmr)
	$(CCLOADER) $(CCLOADFLAGS) -o $(ctrmrexe) $(ctrmr) $(PRLIBS)
	$(MAKE) $(TESTINGdir)/TRMR2D.dat

$(ztrmrexe): $(SCALAPACKLIB) $(ztrmr)
	$(CCLOADER) $(CCLOADFLAGS) -o $(ztrmrexe) $(ztrmr) $(PRLIBS)
	$(MAKE) $(TESTINGdir)/TRMR2D.dat

$(igemr): $(FRC)
$(sgemr): $(FRC)
$(dgemr): $(FRC)
$(cgemr): $(FRC)
$(zgemr): $(FRC)

$(itrmr): $(FRC)
$(strmr): $(FRC)
$(dtrmr): $(FRC)
$(ctrmr): $(FRC)
$(ztrmr): $(FRC)

FRC:
	@FRC=$(FRC)

clean :
	rm -f *.o

.c.o : ; $(CC) -c $(CCFLAGS) $(CDEFS) $*.c
