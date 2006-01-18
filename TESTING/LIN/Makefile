############################################################################
#
#  Program:         ScaLAPACK
#
#  Module:          Makefile
#
#  Purpose:         Linear Equations Testing Makefile
#
#  Creation date:   March 20, 1995
#
#  Modified:
#
#  Send bug reports, comments or suggestions to scalapack@cs.utk.edu
#
############################################################################

include ../../SLmake.inc

sluexe = $(TESTINGdir)/xslu
dluexe = $(TESTINGdir)/xdlu
cluexe = $(TESTINGdir)/xclu
zluexe = $(TESTINGdir)/xzlu

sdbluexe = $(TESTINGdir)/xsdblu
ddbluexe = $(TESTINGdir)/xddblu
cdbluexe = $(TESTINGdir)/xcdblu
zdbluexe = $(TESTINGdir)/xzdblu

sdtluexe = $(TESTINGdir)/xsdtlu
ddtluexe = $(TESTINGdir)/xddtlu
cdtluexe = $(TESTINGdir)/xcdtlu
zdtluexe = $(TESTINGdir)/xzdtlu

sgbluexe = $(TESTINGdir)/xsgblu
dgbluexe = $(TESTINGdir)/xdgblu
cgbluexe = $(TESTINGdir)/xcgblu
zgbluexe = $(TESTINGdir)/xzgblu

slltexe = $(TESTINGdir)/xsllt
dlltexe = $(TESTINGdir)/xdllt
clltexe = $(TESTINGdir)/xcllt
zlltexe = $(TESTINGdir)/xzllt

spblltexe = $(TESTINGdir)/xspbllt
dpblltexe = $(TESTINGdir)/xdpbllt
cpblltexe = $(TESTINGdir)/xcpbllt
zpblltexe = $(TESTINGdir)/xzpbllt

sptlltexe = $(TESTINGdir)/xsptllt
dptlltexe = $(TESTINGdir)/xdptllt
cptlltexe = $(TESTINGdir)/xcptllt
zptlltexe = $(TESTINGdir)/xzptllt

sinvexe = $(TESTINGdir)/xsinv
dinvexe = $(TESTINGdir)/xdinv
cinvexe = $(TESTINGdir)/xcinv
zinvexe = $(TESTINGdir)/xzinv

sqrexe = $(TESTINGdir)/xsqr
dqrexe = $(TESTINGdir)/xdqr
cqrexe = $(TESTINGdir)/xcqr
zqrexe = $(TESTINGdir)/xzqr

slsexe = $(TESTINGdir)/xsls
dlsexe = $(TESTINGdir)/xdls
clsexe = $(TESTINGdir)/xcls
zlsexe = $(TESTINGdir)/xzls

smatgen = psmatgen.o pmatgeninc.o
dmatgen = pdmatgen.o pmatgeninc.o
cmatgen = pcmatgen.o pmatgeninc.o
zmatgen = pzmatgen.o pmatgeninc.o

slinchk = pslaschk.o pslafchk.o
dlinchk = pdlaschk.o pdlafchk.o
clinchk = pclaschk.o pclafchk.o
zlinchk = pzlaschk.o pzlafchk.o

slu = psludriver.o psluinfo.o psgetrrv.o $(smatgen) $(slinchk)
dlu = pdludriver.o pdluinfo.o pdgetrrv.o $(dmatgen) $(dlinchk)
clu = pcludriver.o pcluinfo.o pcgetrrv.o $(cmatgen) $(clinchk)
zlu = pzludriver.o pzluinfo.o pzgetrrv.o $(zmatgen) $(zlinchk)

sdblu = psdbdriver.o psdbinfo.o psdblaschk.o psdbmv1.o psbmatgen.o $(smatgen)
ddblu = pddbdriver.o pddbinfo.o pddblaschk.o pddbmv1.o pdbmatgen.o $(dmatgen)
cdblu = pcdbdriver.o pcdbinfo.o pcdblaschk.o pcdbmv1.o pcbmatgen.o $(cmatgen)
zdblu = pzdbdriver.o pzdbinfo.o pzdblaschk.o pzdbmv1.o pzbmatgen.o $(zmatgen)

sdtlu = psdtdriver.o psdtinfo.o psdtlaschk.o psdbmv1.o psbmatgen.o $(smatgen)
ddtlu = pddtdriver.o pddtinfo.o pddtlaschk.o pddbmv1.o pdbmatgen.o $(dmatgen)
cdtlu = pcdtdriver.o pcdtinfo.o pcdtlaschk.o pcdbmv1.o pcbmatgen.o $(cmatgen)
zdtlu = pzdtdriver.o pzdtinfo.o pzdtlaschk.o pzdbmv1.o pzbmatgen.o $(zmatgen)

sgblu = psgbdriver.o psgbinfo.o psdblaschk.o psgbmv1.o psbmatgen.o $(smatgen)
dgblu = pdgbdriver.o pdgbinfo.o pddblaschk.o pdgbmv1.o pdbmatgen.o $(dmatgen)
cgblu = pcgbdriver.o pcgbinfo.o pcdblaschk.o pcgbmv1.o pcbmatgen.o $(cmatgen)
zgblu = pzgbdriver.o pzgbinfo.o pzdblaschk.o pzgbmv1.o pzbmatgen.o $(zmatgen)


sllt = pslltdriver.o pslltinfo.o pspotrrv.o $(smatgen) $(slinchk)
dllt = pdlltdriver.o pdlltinfo.o pdpotrrv.o $(dmatgen) $(dlinchk)
cllt = pclltdriver.o pclltinfo.o pcpotrrv.o $(cmatgen) $(clinchk)
zllt = pzlltdriver.o pzlltinfo.o pzpotrrv.o $(zmatgen) $(zlinchk)

spbllt = pspbdriver.o pspbinfo.o pspblaschk.o pspbmv1.o psbmatgen.o $(smatgen)
dpbllt = pdpbdriver.o pdpbinfo.o pdpblaschk.o pdpbmv1.o pdbmatgen.o $(dmatgen)
cpbllt = pcpbdriver.o pcpbinfo.o pcpblaschk.o pcpbmv1.o pcbmatgen.o $(cmatgen)
zpbllt = pzpbdriver.o pzpbinfo.o pzpblaschk.o pzpbmv1.o pzbmatgen.o $(zmatgen)

sptllt = psptdriver.o psptinfo.o psptlaschk.o pspbmv1.o psbmatgen.o $(smatgen)
dptllt = pdptdriver.o pdptinfo.o pdptlaschk.o pdpbmv1.o pdbmatgen.o $(dmatgen)
cptllt = pcptdriver.o pcptinfo.o pcptlaschk.o pcpbmv1.o pcbmatgen.o $(cmatgen)
zptllt = pzptdriver.o pzptinfo.o pzptlaschk.o pzpbmv1.o pzbmatgen.o $(zmatgen)

sinv = psinvdriver.o psinvinfo.o psinvchk.o $(smatgen)
dinv = pdinvdriver.o pdinvinfo.o pdinvchk.o $(dmatgen)
cinv = pcinvdriver.o pcinvinfo.o pcinvchk.o $(cmatgen)
zinv = pzinvdriver.o pzinvinfo.o pzinvchk.o $(zmatgen)

sqr = psqrdriver.o psqrinfo.o psgeqrrv.o psgeqlrv.o psgelqrv.o \
      psgerqrv.o   pstzrzrv.o pslafchk.o $(smatgen)
dqr = pdqrdriver.o pdqrinfo.o pdgeqrrv.o pdgeqlrv.o pdgelqrv.o \
      pdgerqrv.o   pdtzrzrv.o pdlafchk.o $(dmatgen)
cqr = pcqrdriver.o pcqrinfo.o pcgeqrrv.o pcgeqlrv.o pcgelqrv.o \
      pcgerqrv.o   pctzrzrv.o pclafchk.o $(cmatgen)
zqr = pzqrdriver.o pzqrinfo.o pzgeqrrv.o pzgeqlrv.o pzgelqrv.o \
      pzgerqrv.o   pztzrzrv.o pzlafchk.o $(zmatgen)

sls = pslsdriver.o pslsinfo.o psqrt13.o psqrt14.o psqrt16.o \
      psqrt17.o $(smatgen)
dls = pdlsdriver.o pdlsinfo.o pdqrt13.o pdqrt14.o pdqrt16.o \
      pdqrt17.o $(dmatgen)
cls = pclsdriver.o pclsinfo.o pcqrt13.o pcqrt14.o pcqrt16.o \
      pcqrt17.o $(cmatgen)
zls = pzlsdriver.o pzlsinfo.o pzqrt13.o pzqrt14.o pzqrt16.o \
      pzqrt17.o $(zmatgen)

all : single double complex complex16

single: $(sluexe) $(sdbluexe) $(sdtluexe) $(sgbluexe) $(slltexe) \
        $(spblltexe) $(sptlltexe) $(sinvexe) $(sqrexe) $(slsexe)

double: $(dluexe) $(ddbluexe) $(ddtluexe) $(dgbluexe) $(dlltexe) \
        $(dpblltexe) $(dptlltexe) $(dinvexe) $(dqrexe) $(dlsexe)

complex: $(cluexe) $(cdbluexe) $(cdtluexe) $(cgbluexe) $(clltexe) \
         $(cpblltexe) $(cptlltexe) $(cinvexe) $(cqrexe) $(clsexe)

complex16: $(zluexe) $(zdbluexe) $(zdtluexe) $(zgbluexe) $(zlltexe) \
           $(zpblltexe) $(zptlltexe) $(zinvexe) $(zqrexe) $(zlsexe)

$(TESTINGdir)/LU.dat: ../LU.dat
	cp ../LU.dat $(TESTINGdir)

$(sluexe): $(SCALAPACKLIB) $(slu)
	$(F77LOADER) $(F77LOADFLAGS) -o $(sluexe) $(slu) $(LIBS)
	$(MAKE) $(TESTINGdir)/LU.dat
$(dluexe): $(SCALAPACKLIB) $(dlu)
	$(F77LOADER) $(F77LOADFLAGS) -o $(dluexe) $(dlu) $(LIBS)
	$(MAKE) $(TESTINGdir)/LU.dat
$(cluexe): $(SCALAPACKLIB) $(clu)
	$(F77LOADER) $(F77LOADFLAGS) -o $(cluexe) $(clu) $(LIBS)
	$(MAKE) $(TESTINGdir)/LU.dat
$(zluexe): $(SCALAPACKLIB) $(zlu)
	$(F77LOADER) $(F77LOADFLAGS) -o $(zluexe) $(zlu) $(LIBS)
	$(MAKE) $(TESTINGdir)/LU.dat

$(TESTINGdir)/BLU.dat: ../BLU.dat
	cp ../BLU.dat $(TESTINGdir)

$(sdbluexe): $(SCALAPACKLIB) $(sdblu)
	$(F77LOADER) $(F77LOADFLAGS) -o $(sdbluexe) $(sdblu) $(LIBS)
	$(MAKE) $(TESTINGdir)/BLU.dat
$(ddbluexe): $(SCALAPACKLIB) $(ddblu)
	$(F77LOADER) $(F77LOADFLAGS) -o $(ddbluexe) $(ddblu) $(LIBS)
	$(MAKE) $(TESTINGdir)/BLU.dat
$(cdbluexe): $(SCALAPACKLIB) $(cdblu)
	$(F77LOADER) $(F77LOADFLAGS) -o $(cdbluexe) $(cdblu) $(LIBS)
	$(MAKE) $(TESTINGdir)/BLU.dat
$(zdbluexe): $(SCALAPACKLIB) $(zdblu)
	$(F77LOADER) $(F77LOADFLAGS) -o $(zdbluexe) $(zdblu) $(LIBS)
	$(MAKE) $(TESTINGdir)/BLU.dat

$(sdtluexe): $(SCALAPACKLIB) $(sdtlu)
	$(F77LOADER) $(F77LOADFLAGS) -o $(sdtluexe) $(sdtlu) $(LIBS)
	$(MAKE) $(TESTINGdir)/BLU.dat
$(ddtluexe): $(SCALAPACKLIB) $(ddtlu)
	$(F77LOADER) $(F77LOADFLAGS) -o $(ddtluexe) $(ddtlu) $(LIBS)
	$(MAKE) $(TESTINGdir)/BLU.dat
$(cdtluexe): $(SCALAPACKLIB) $(cdtlu)
	$(F77LOADER) $(F77LOADFLAGS) -o $(cdtluexe) $(cdtlu) $(LIBS)
	$(MAKE) $(TESTINGdir)/BLU.dat
$(zdtluexe): $(SCALAPACKLIB) $(zdtlu)
	$(F77LOADER) $(F77LOADFLAGS) -o $(zdtluexe) $(zdtlu) $(LIBS)
	$(MAKE) $(TESTINGdir)/BLU.dat

$(sgbluexe): $(SCALAPACKLIB) $(sgblu)
	$(F77LOADER) $(F77LOADFLAGS) -o $(sgbluexe) $(sgblu) $(LIBS)
	$(MAKE) $(TESTINGdir)/BLU.dat
$(dgbluexe): $(SCALAPACKLIB) $(dgblu)
	$(F77LOADER) $(F77LOADFLAGS) -o $(dgbluexe) $(dgblu) $(LIBS)
	$(MAKE) $(TESTINGdir)/BLU.dat
$(cgbluexe): $(SCALAPACKLIB) $(cgblu)
	$(F77LOADER) $(F77LOADFLAGS) -o $(cgbluexe) $(cgblu) $(LIBS)
	$(MAKE) $(TESTINGdir)/BLU.dat
$(zgbluexe): $(SCALAPACKLIB) $(zgblu)
	$(F77LOADER) $(F77LOADFLAGS) -o $(zgbluexe) $(zgblu) $(LIBS)
	$(MAKE) $(TESTINGdir)/BLU.dat


$(TESTINGdir)/LLT.dat: ../LLT.dat
	cp ../LLT.dat $(TESTINGdir)

$(slltexe): $(SCALAPACKLIB) $(sllt)
	$(F77LOADER) $(F77LOADFLAGS) -o $(slltexe) $(sllt) $(LIBS)
	$(MAKE) $(TESTINGdir)/LLT.dat
$(dlltexe): $(SCALAPACKLIB) $(dllt)
	$(F77LOADER) $(F77LOADFLAGS) -o $(dlltexe) $(dllt) $(LIBS)
	$(MAKE) $(TESTINGdir)/LLT.dat
$(clltexe): $(SCALAPACKLIB) $(cllt)
	$(F77LOADER) $(F77LOADFLAGS) -o $(clltexe) $(cllt) $(LIBS)
	$(MAKE) $(TESTINGdir)/LLT.dat
$(zlltexe): $(SCALAPACKLIB) $(zllt)
	$(F77LOADER) $(F77LOADFLAGS) -o $(zlltexe) $(zllt) $(LIBS)
	$(MAKE) $(TESTINGdir)/LLT.dat

$(TESTINGdir)/BLLT.dat: ../BLLT.dat
	cp ../BLLT.dat $(TESTINGdir)

$(spblltexe): $(SCALAPACKLIB) $(spbllt)
	$(F77LOADER) $(F77LOADFLAGS) -o $(spblltexe) $(spbllt) $(LIBS)
	$(MAKE) $(TESTINGdir)/BLLT.dat
$(dpblltexe): $(SCALAPACKLIB) $(dpbllt)
	$(F77LOADER) $(F77LOADFLAGS) -o $(dpblltexe) $(dpbllt) $(LIBS)
	$(MAKE) $(TESTINGdir)/BLLT.dat
$(cpblltexe): $(SCALAPACKLIB) $(cpbllt)
	$(F77LOADER) $(F77LOADFLAGS) -o $(cpblltexe) $(cpbllt) $(LIBS)
	$(MAKE) $(TESTINGdir)/BLLT.dat
$(zpblltexe): $(SCALAPACKLIB) $(zpbllt)
	$(F77LOADER) $(F77LOADFLAGS) -o $(zpblltexe) $(zpbllt) $(LIBS)
	$(MAKE) $(TESTINGdir)/BLLT.dat

$(sptlltexe): $(SCALAPACKLIB) $(sptllt)
	$(F77LOADER) $(F77LOADFLAGS) -o $(sptlltexe) $(sptllt) $(LIBS)
	$(MAKE) $(TESTINGdir)/BLLT.dat
$(dptlltexe): $(SCALAPACKLIB) $(dptllt)
	$(F77LOADER) $(F77LOADFLAGS) -o $(dptlltexe) $(dptllt) $(LIBS)
	$(MAKE) $(TESTINGdir)/BLLT.dat
$(cptlltexe): $(SCALAPACKLIB) $(cptllt)
	$(F77LOADER) $(F77LOADFLAGS) -o $(cptlltexe) $(cptllt) $(LIBS)
	$(MAKE) $(TESTINGdir)/BLLT.dat
$(zptlltexe): $(SCALAPACKLIB) $(zptllt)
	$(F77LOADER) $(F77LOADFLAGS) -o $(zptlltexe) $(zptllt) $(LIBS)
	$(MAKE) $(TESTINGdir)/BLLT.dat

$(TESTINGdir)/INV.dat: ../INV.dat
	cp ../INV.dat $(TESTINGdir)

$(sinvexe): $(SCALAPACKLIB) $(sinv)
	$(F77LOADER) $(F77LOADFLAGS) -o $(sinvexe) $(sinv) $(LIBS)
	$(MAKE) $(TESTINGdir)/INV.dat
$(dinvexe): $(SCALAPACKLIB) $(dinv)
	$(F77LOADER) $(F77LOADFLAGS) -o $(dinvexe) $(dinv) $(LIBS)
	$(MAKE) $(TESTINGdir)/INV.dat
$(cinvexe): $(SCALAPACKLIB) $(cinv)
	$(F77LOADER) $(F77LOADFLAGS) -o $(cinvexe) $(cinv) $(LIBS)
	$(MAKE) $(TESTINGdir)/INV.dat
$(zinvexe): $(SCALAPACKLIB) $(zinv)
	$(F77LOADER) $(F77LOADFLAGS) -o $(zinvexe) $(zinv) $(LIBS)
	$(MAKE) $(TESTINGdir)/INV.dat

$(TESTINGdir)/QR.dat: ../QR.dat
	cp ../QR.dat $(TESTINGdir)

$(sqrexe): $(SCALAPACKLIB) $(sqr)
	$(F77LOADER) $(F77LOADFLAGS) -o $(sqrexe) $(sqr) $(LIBS)
	$(MAKE) $(TESTINGdir)/QR.dat
$(dqrexe): $(SCALAPACKLIB) $(dqr)
	$(F77LOADER) $(F77LOADFLAGS) -o $(dqrexe) $(dqr) $(LIBS)
	$(MAKE) $(TESTINGdir)/QR.dat
$(cqrexe): $(SCALAPACKLIB) $(cqr)
	$(F77LOADER) $(F77LOADFLAGS) -o $(cqrexe) $(cqr) $(LIBS)
	$(MAKE) $(TESTINGdir)/QR.dat
$(zqrexe): $(SCALAPACKLIB) $(zqr)
	$(F77LOADER) $(F77LOADFLAGS) -o $(zqrexe) $(zqr) $(LIBS)
	$(MAKE) $(TESTINGdir)/QR.dat

$(TESTINGdir)/LS.dat: ../LS.dat
	cp ../LS.dat $(TESTINGdir)

$(slsexe): $(SCALAPACKLIB) $(sls)
	$(F77LOADER) $(F77LOADFLAGS) -o $(slsexe) $(sls) $(LIBS)
	$(MAKE) $(TESTINGdir)/LS.dat
$(dlsexe): $(SCALAPACKLIB) $(dls)
	$(F77LOADER) $(F77LOADFLAGS) -o $(dlsexe) $(dls) $(LIBS)
	$(MAKE) $(TESTINGdir)/LS.dat
$(clsexe): $(SCALAPACKLIB) $(cls)
	$(F77LOADER) $(F77LOADFLAGS) -o $(clsexe) $(cls) $(LIBS)
	$(MAKE) $(TESTINGdir)/LS.dat
$(zlsexe): $(SCALAPACKLIB) $(zls)
	$(F77LOADER) $(F77LOADFLAGS) -o $(zlsexe) $(zls) $(LIBS)
	$(MAKE) $(TESTINGdir)/LS.dat

$(slu): $(FRC)
$(dlu): $(FRC)
$(clu): $(FRC)
$(zlu): $(FRC)

$(sdblu): $(FRC)
$(ddblu): $(FRC)
$(cdblu): $(FRC)
$(zdblu): $(FRC)

$(sdtlu): $(FRC)
$(ddtlu): $(FRC)
$(cdtlu): $(FRC)
$(zdtlu): $(FRC)


$(sllt): $(FRC)
$(dllt): $(FRC)
$(cllt): $(FRC)
$(zllt): $(FRC)

$(spbllt): $(FRC)
$(dpbllt): $(FRC)
$(cpbllt): $(FRC)
$(zpbllt): $(FRC)

$(sptllt): $(FRC)
$(dptllt): $(FRC)
$(cptllt): $(FRC)
$(zptllt): $(FRC)

$(sinv): $(FRC)
$(dinv): $(FRC)
$(cinv): $(FRC)
$(zinv): $(FRC)

$(sqr): $(FRC)
$(dqr): $(FRC)
$(cqr): $(FRC)
$(zqr): $(FRC)

$(sls): $(FRC)
$(dls): $(FRC)
$(cls): $(FRC)
$(zls): $(FRC)

FRC:
	@FRC=$(FRC)

clean :
	rm -f *.o

psludriver.o: psludriver.f
	$(F77) $(DRVOPTS) -c $<
pdludriver.o: pdludriver.f
	$(F77) $(DRVOPTS) -c $<
pcludriver.o: pcludriver.f
	$(F77) $(DRVOPTS) -c $<
pzludriver.o: pzludriver.f
	$(F77) $(DRVOPTS) -c $<

psdbdriver.o: psdbdriver.f
	$(F77) $(DRVOPTS) -c $<
pddbdriver.o: pddbdriver.f
	$(F77) $(DRVOPTS) -c $<
pcdbdriver.o: pcdbdriver.f
	$(F77) $(DRVOPTS) -c $<
pzdbdriver.o: pzdbdriver.f
	$(F77) $(DRVOPTS) -c $<

psdtdriver.o: psdtdriver.f
	$(F77) $(DRVOPTS) -c $<
pddtdriver.o: pddtdriver.f
	$(F77) $(DRVOPTS) -c $<
pcdtdriver.o: pcdtdriver.f
	$(F77) $(DRVOPTS) -c $<
pzdtdriver.o: pzdtdriver.f
	$(F77) $(DRVOPTS) -c $<

psgbdriver.o: psgbdriver.f
	$(F77) $(DRVOPTS) -c $<
pdgbdriver.o: pdgbdriver.f
	$(F77) $(DRVOPTS) -c $<
pcgbdriver.o: pcgbdriver.f
	$(F77) $(DRVOPTS) -c $<
pzgbdriver.o: pzgbdriver.f
	$(F77) $(DRVOPTS) -c $<

pslltdriver.o: pslltdriver.f
	$(F77) $(DRVOPTS) -c $<
pdlltdriver.o: pdlltdriver.f
	$(F77) $(DRVOPTS) -c $<
pclltdriver.o: pclltdriver.f
	$(F77) $(DRVOPTS) -c $<
pzlltdriver.o: pzlltdriver.f
	$(F77) $(DRVOPTS) -c $<

pspbdriver.o: pspbdriver.f
	$(F77) $(DRVOPTS) -c $<
pdpbdriver.o: pdpbdriver.f
	$(F77) $(DRVOPTS) -c $<
pcpbdriver.o: pcpbdriver.f
	$(F77) $(DRVOPTS) -c $<
pzpbdriver.o: pzpbdriver.f
	$(F77) $(DRVOPTS) -c $<

psptdriver.o: psptdriver.f
	$(F77) $(DRVOPTS) -c $<
pdptdriver.o: pdptdriver.f
	$(F77) $(DRVOPTS) -c $<
pcptdriver.o: pcptdriver.f
	$(F77) $(DRVOPTS) -c $<
pzptdriver.o: pzptdriver.f
	$(F77) $(DRVOPTS) -c $<

psinvdriver.o: psinvdriver.f
	$(F77) $(DRVOPTS) -c $<
pdinvdriver.o: pdinvdriver.f
	$(F77) $(DRVOPTS) -c $<
pcinvdriver.o: pcinvdriver.f
	$(F77) $(DRVOPTS) -c $<
pzinvdriver.o: pzinvdriver.f
	$(F77) $(DRVOPTS) -c $<

psqrdriver.o: psqrdriver.f
	$(F77) $(DRVOPTS) -c $<
pdqrdriver.o: pdqrdriver.f
	$(F77) $(DRVOPTS) -c $<
pcqrdriver.o: pcqrdriver.f
	$(F77) $(DRVOPTS) -c $<
pzqrdriver.o: pzqrdriver.f
	$(F77) $(DRVOPTS) -c $<

pslsdriver.o: pslsdriver.f
	$(F77) $(DRVOPTS) -c $<
pdlsdriver.o: pdlsdriver.f
	$(F77) $(DRVOPTS) -c $<
pclsdriver.o: pclsdriver.f
	$(F77) $(DRVOPTS) -c $<
pzlsdriver.o: pzlsdriver.f
	$(F77) $(DRVOPTS) -c $<

.f.o : ; $(F77) -c $(F77FLAGS) $*.f
