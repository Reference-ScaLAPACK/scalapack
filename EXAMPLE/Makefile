include ../SLmake.inc

TESTOBJS = psscaex.o pdscaexinfo.o
TESTOBJD = pdscaex.o pdscaexinfo.o
TESTOBJC = pcscaex.o pdscaexinfo.o
TESTOBJZ = pzscaex.o pdscaexinfo.o

all: pdscaex pcscaex pzscaex psscaex

single: psscaex

double: pdscaex

complex: pcscaex

complex16: pzscaex

psscaex: $(TESTOBJS)
	$(FCLOADER) $(FCLOADFLAGS) -o xsscaex $(TESTOBJS) ../$(SCALAPACKLIB) $(LIBS)

pzscaex: $(TESTOBJZ)
	$(FCLOADER) $(FCLOADFLAGS) -o xzscaex $(TESTOBJZ) ../$(SCALAPACKLIB) $(LIBS)

pcscaex: $(TESTOBJC)
	$(FCLOADER) $(FCLOADFLAGS) -o xcscaex $(TESTOBJC) ../$(SCALAPACKLIB) $(LIBS)

pdscaex: $(TESTOBJD)
	$(FCLOADER) $(FCLOADFLAGS) -o xdscaex $(TESTOBJD) ../$(SCALAPACKLIB) $(LIBS)

clean :
	rm -f $(TESTOBJS) $(TESTOBJD) $(TESTOBJZ) $(TESTOBJC) xsscaex xzscaex xcscaex xdscaex

.f.o : ; $(FC) -c $(FCFLAGS) $*.f

.c.o : ; $(CC) -c $(CDEFS) $(CCFLAGS) $*.c

