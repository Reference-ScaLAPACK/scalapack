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

$(TESTINGdir)/SCAEX.dat: SCAEX.dat
	cp SCAEX.dat $(TESTINGdir)

$(TESTINGdir)/SSCAEXMAT.dat: SSCAEXMAT.dat
	cp SSCAEXMAT.dat $(TESTINGdir)

$(TESTINGdir)/SSCAEXRHS.dat: SSCAEXRHS.dat
	cp SSCAEXRHS.dat $(TESTINGdir)

psscaex: $(TESTOBJS)
	$(F77LOADER) $(F77LOADFLAGS) -o $(TESTINGdir)/xsscaex $(TESTOBJS) $(LIBS)
	$(MAKE) $(TESTINGdir)/SCAEX.dat
	$(MAKE) $(TESTINGdir)/SSCAEXMAT.dat
	$(MAKE) $(TESTINGdir)/SSCAEXRHS.dat

$(TESTINGdir)/ZSCAEXMAT.dat: ZSCAEXMAT.dat
	cp ZSCAEXMAT.dat $(TESTINGdir)

$(TESTINGdir)/ZSCAEXRHS.dat: ZSCAEXRHS.dat
	cp ZSCAEXRHS.dat $(TESTINGdir)

pzscaex: $(TESTOBJZ)
	$(F77LOADER) $(F77LOADFLAGS) -o $(TESTINGdir)/xzscaex $(TESTOBJZ) $(LIBS)
	$(MAKE) $(TESTINGdir)/SCAEX.dat
	$(MAKE) $(TESTINGdir)/ZSCAEXMAT.dat
	$(MAKE) $(TESTINGdir)/ZSCAEXRHS.dat

$(TESTINGdir)/CSCAEXMAT.dat: CSCAEXMAT.dat
	cp CSCAEXMAT.dat $(TESTINGdir)

$(TESTINGdir)/CSCAEXRHS.dat: CSCAEXRHS.dat
	cp CSCAEXRHS.dat $(TESTINGdir)

pcscaex: $(TESTOBJC)
	$(F77LOADER) $(F77LOADFLAGS) -o $(TESTINGdir)/xcscaex $(TESTOBJC) $(LIBS)
	$(MAKE) $(TESTINGdir)/SCAEX.dat
	$(MAKE) $(TESTINGdir)/CSCAEXMAT.dat
	$(MAKE) $(TESTINGdir)/CSCAEXRHS.dat

$(TESTINGdir)/DSCAEXMAT.dat: DSCAEXMAT.dat
	cp DSCAEXMAT.dat $(TESTINGdir)

$(TESTINGdir)/DSCAEXRHS.dat: DSCAEXRHS.dat
	cp DSCAEXRHS.dat $(TESTINGdir)

pdscaex: $(TESTOBJD)
	$(F77LOADER) $(F77LOADFLAGS) -o $(TESTINGdir)/xdscaex $(TESTOBJD) $(LIBS)
	$(MAKE) $(TESTINGdir)/SCAEX.dat
	$(MAKE) $(TESTINGdir)/DSCAEXMAT.dat
	$(MAKE) $(TESTINGdir)/DSCAEXRHS.dat

clean :
	rm -f *.o

.f.o : ; $(F77) -c $(F77FLAGS) $*.f

.c.o : ; $(CC) -c $(CCFLAGS) $*.c

