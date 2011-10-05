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
	$(F77LOADER) $(F77LOADFLAGS) -o xsscaex $(TESTOBJS) ../$(SCALAPACKLIB) $(LIBS)

pzscaex: $(TESTOBJZ)
	$(F77LOADER) $(F77LOADFLAGS) -o xzscaex $(TESTOBJZ) ../$(SCALAPACKLIB) $(LIBS)

pcscaex: $(TESTOBJC)
	$(F77LOADER) $(F77LOADFLAGS) -o xcscaex $(TESTOBJC) ../$(SCALAPACKLIB) $(LIBS)

pdscaex: $(TESTOBJD)
	$(F77LOADER) $(F77LOADFLAGS) -o xdscaex $(TESTOBJD) ../$(SCALAPACKLIB) $(LIBS)

clean :
	rm -f $(TESTOBJS) $(TESTOBJD) $(TESTOBJZ) $(TESTOBJC) xsscaex xzscaex xcscaex xdscaex

.f.o : ; $(F77) -c $(F77FLAGS) $*.f

.c.o : ; $(CC) -c $(CCFLAGS) $*.c

