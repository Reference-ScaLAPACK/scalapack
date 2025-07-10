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

sluexe = ../xslu
dluexe = ../xdlu
cluexe = ../xclu
zluexe = ../xzlu

sdbluexe = ../xsdblu
ddbluexe = ../xddblu
cdbluexe = ../xcdblu
zdbluexe = ../xzdblu

sdtluexe = ../xsdtlu
ddtluexe = ../xddtlu
cdtluexe = ../xcdtlu
zdtluexe = ../xzdtlu

sgbluexe = ../xsgblu
dgbluexe = ../xdgblu
cgbluexe = ../xcgblu
zgbluexe = ../xzgblu

slltexe = ../xsllt
dlltexe = ../xdllt
clltexe = ../xcllt
zlltexe = ../xzllt

spblltexe = ../xspbllt
dpblltexe = ../xdpbllt
cpblltexe = ../xcpbllt
zpblltexe = ../xzpbllt

sptlltexe = ../xsptllt
dptlltexe = ../xdptllt
cptlltexe = ../xcptllt
zptlltexe = ../xzptllt

sinvexe = ../xsinv
dinvexe = ../xdinv
cinvexe = ../xcinv
zinvexe = ../xzinv

sqrexe = ../xsqr
dqrexe = ../xdqr
cqrexe = ../xcqr
zqrexe = ../xzqr

slsexe = ../xsls
dlsexe = ../xdls
clsexe = ../xcls
zlsexe = ../xzls

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

$(sluexe): ../../$(SCALAPACKLIB) $(slu)
	$(FCLOADER) $(FCLOADFLAGS) -o $(sluexe) $(slu) ../../$(SCALAPACKLIB) $(LIBS)

$(dluexe): ../../$(SCALAPACKLIB) $(dlu)
	$(FCLOADER) $(FCLOADFLAGS) -o $(dluexe) $(dlu) ../../$(SCALAPACKLIB) $(LIBS)

$(cluexe): ../../$(SCALAPACKLIB) $(clu)
	$(FCLOADER) $(FCLOADFLAGS) -o $(cluexe) $(clu) ../../$(SCALAPACKLIB) $(LIBS)

$(zluexe): ../../$(SCALAPACKLIB) $(zlu)
	$(FCLOADER) $(FCLOADFLAGS) -o $(zluexe) $(zlu) ../../$(SCALAPACKLIB) $(LIBS)

$(sdbluexe): ../../$(SCALAPACKLIB) $(sdblu)
	$(FCLOADER) $(FCLOADFLAGS) -o $(sdbluexe) $(sdblu) ../../$(SCALAPACKLIB) $(LIBS)

$(ddbluexe): ../../$(SCALAPACKLIB) $(ddblu)
	$(FCLOADER) $(FCLOADFLAGS) -o $(ddbluexe) $(ddblu) ../../$(SCALAPACKLIB) $(LIBS)

$(cdbluexe): ../../$(SCALAPACKLIB) $(cdblu)
	$(FCLOADER) $(FCLOADFLAGS) -o $(cdbluexe) $(cdblu) ../../$(SCALAPACKLIB) $(LIBS)

$(zdbluexe): ../../$(SCALAPACKLIB) $(zdblu)
	$(FCLOADER) $(FCLOADFLAGS) -o $(zdbluexe) $(zdblu) ../../$(SCALAPACKLIB) $(LIBS)

$(sdtluexe): ../../$(SCALAPACKLIB) $(sdtlu)
	$(FCLOADER) $(FCLOADFLAGS) -o $(sdtluexe) $(sdtlu) ../../$(SCALAPACKLIB) $(LIBS)

$(ddtluexe): ../../$(SCALAPACKLIB) $(ddtlu)
	$(FCLOADER) $(FCLOADFLAGS) -o $(ddtluexe) $(ddtlu) ../../$(SCALAPACKLIB) $(LIBS)

$(cdtluexe): ../../$(SCALAPACKLIB) $(cdtlu)
	$(FCLOADER) $(FCLOADFLAGS) -o $(cdtluexe) $(cdtlu) ../../$(SCALAPACKLIB) $(LIBS)

$(zdtluexe): ../../$(SCALAPACKLIB) $(zdtlu)
	$(FCLOADER) $(FCLOADFLAGS) -o $(zdtluexe) $(zdtlu) ../../$(SCALAPACKLIB) $(LIBS)

$(sgbluexe): ../../$(SCALAPACKLIB) $(sgblu)
	$(FCLOADER) $(FCLOADFLAGS) -o $(sgbluexe) $(sgblu) ../../$(SCALAPACKLIB) $(LIBS)

$(dgbluexe): ../../$(SCALAPACKLIB) $(dgblu)
	$(FCLOADER) $(FCLOADFLAGS) -o $(dgbluexe) $(dgblu) ../../$(SCALAPACKLIB) $(LIBS)

$(cgbluexe): ../../$(SCALAPACKLIB) $(cgblu)
	$(FCLOADER) $(FCLOADFLAGS) -o $(cgbluexe) $(cgblu) ../../$(SCALAPACKLIB) $(LIBS)

$(zgbluexe): ../../$(SCALAPACKLIB) $(zgblu)
	$(FCLOADER) $(FCLOADFLAGS) -o $(zgbluexe) $(zgblu) ../../$(SCALAPACKLIB) $(LIBS)

$(slltexe): ../../$(SCALAPACKLIB) $(sllt)
	$(FCLOADER) $(FCLOADFLAGS) -o $(slltexe) $(sllt) ../../$(SCALAPACKLIB) $(LIBS)

$(dlltexe): ../../$(SCALAPACKLIB) $(dllt)
	$(FCLOADER) $(FCLOADFLAGS) -o $(dlltexe) $(dllt) ../../$(SCALAPACKLIB) $(LIBS)

$(clltexe): ../../$(SCALAPACKLIB) $(cllt)
	$(FCLOADER) $(FCLOADFLAGS) -o $(clltexe) $(cllt) ../../$(SCALAPACKLIB) $(LIBS)

$(zlltexe): ../../$(SCALAPACKLIB) $(zllt)
	$(FCLOADER) $(FCLOADFLAGS) -o $(zlltexe) $(zllt) ../../$(SCALAPACKLIB) $(LIBS)

$(spblltexe): ../../$(SCALAPACKLIB) $(spbllt)
	$(FCLOADER) $(FCLOADFLAGS) -o $(spblltexe) $(spbllt) ../../$(SCALAPACKLIB) $(LIBS)

$(dpblltexe): ../../$(SCALAPACKLIB) $(dpbllt)
	$(FCLOADER) $(FCLOADFLAGS) -o $(dpblltexe) $(dpbllt) ../../$(SCALAPACKLIB) $(LIBS)

$(cpblltexe): ../../$(SCALAPACKLIB) $(cpbllt)
	$(FCLOADER) $(FCLOADFLAGS) -o $(cpblltexe) $(cpbllt) ../../$(SCALAPACKLIB) $(LIBS)

$(zpblltexe): ../../$(SCALAPACKLIB) $(zpbllt)
	$(FCLOADER) $(FCLOADFLAGS) -o $(zpblltexe) $(zpbllt) ../../$(SCALAPACKLIB) $(LIBS)

$(sptlltexe): ../../$(SCALAPACKLIB) $(sptllt)
	$(FCLOADER) $(FCLOADFLAGS) -o $(sptlltexe) $(sptllt) ../../$(SCALAPACKLIB) $(LIBS)

$(dptlltexe): ../../$(SCALAPACKLIB) $(dptllt)
	$(FCLOADER) $(FCLOADFLAGS) -o $(dptlltexe) $(dptllt) ../../$(SCALAPACKLIB) $(LIBS)

$(cptlltexe): ../../$(SCALAPACKLIB) $(cptllt)
	$(FCLOADER) $(FCLOADFLAGS) -o $(cptlltexe) $(cptllt) ../../$(SCALAPACKLIB) $(LIBS)

$(zptlltexe): ../../$(SCALAPACKLIB) $(zptllt)
	$(FCLOADER) $(FCLOADFLAGS) -o $(zptlltexe) $(zptllt) ../../$(SCALAPACKLIB) $(LIBS)

$(sinvexe): ../../$(SCALAPACKLIB) $(sinv)
	$(FCLOADER) $(FCLOADFLAGS) -o $(sinvexe) $(sinv) ../../$(SCALAPACKLIB) $(LIBS)

$(dinvexe): ../../$(SCALAPACKLIB) $(dinv)
	$(FCLOADER) $(FCLOADFLAGS) -o $(dinvexe) $(dinv) ../../$(SCALAPACKLIB) $(LIBS)

$(cinvexe): ../../$(SCALAPACKLIB) $(cinv)
	$(FCLOADER) $(FCLOADFLAGS) -o $(cinvexe) $(cinv) ../../$(SCALAPACKLIB) $(LIBS)

$(zinvexe): ../../$(SCALAPACKLIB) $(zinv)
	$(FCLOADER) $(FCLOADFLAGS) -o $(zinvexe) $(zinv) ../../$(SCALAPACKLIB) $(LIBS)

$(sqrexe): ../../$(SCALAPACKLIB) $(sqr)
	$(FCLOADER) $(FCLOADFLAGS) -o $(sqrexe) $(sqr) ../../$(SCALAPACKLIB) $(LIBS)

$(dqrexe): ../../$(SCALAPACKLIB) $(dqr)
	$(FCLOADER) $(FCLOADFLAGS) -o $(dqrexe) $(dqr) ../../$(SCALAPACKLIB) $(LIBS)

$(cqrexe): ../../$(SCALAPACKLIB) $(cqr)
	$(FCLOADER) $(FCLOADFLAGS) -o $(cqrexe) $(cqr) ../../$(SCALAPACKLIB) $(LIBS)

$(zqrexe): ../../$(SCALAPACKLIB) $(zqr)
	$(FCLOADER) $(FCLOADFLAGS) -o $(zqrexe) $(zqr) ../../$(SCALAPACKLIB) $(LIBS)

$(slsexe): ../../$(SCALAPACKLIB) $(sls)
	$(FCLOADER) $(FCLOADFLAGS) -o $(slsexe) $(sls) ../../$(SCALAPACKLIB) $(LIBS)

$(dlsexe): ../../$(SCALAPACKLIB) $(dls)
	$(FCLOADER) $(FCLOADFLAGS) -o $(dlsexe) $(dls) ../../$(SCALAPACKLIB) $(LIBS)

$(clsexe): ../../$(SCALAPACKLIB) $(cls)
	$(FCLOADER) $(FCLOADFLAGS) -o $(clsexe) $(cls) ../../$(SCALAPACKLIB) $(LIBS)

$(zlsexe): ../../$(SCALAPACKLIB) $(zls)
	$(FCLOADER) $(FCLOADFLAGS) -o $(zlsexe) $(zls) ../../$(SCALAPACKLIB) $(LIBS)

clean :
	rm -f *.o

psludriver.o: psludriver.f
	$(FC) $(FCFLAGS) -c $<
pdludriver.o: pdludriver.f
	$(FC) $(FCFLAGS) -c $<
pcludriver.o: pcludriver.f
	$(FC) $(FCFLAGS) -c $<
pzludriver.o: pzludriver.f
	$(FC) $(FCFLAGS) -c $<

psdbdriver.o: psdbdriver.f
	$(FC) $(FCFLAGS) -c $<
pddbdriver.o: pddbdriver.f
	$(FC) $(FCFLAGS) -c $<
pcdbdriver.o: pcdbdriver.f
	$(FC) $(FCFLAGS) -c $<
pzdbdriver.o: pzdbdriver.f
	$(FC) $(FCFLAGS) -c $<

psdtdriver.o: psdtdriver.f
	$(FC) $(FCFLAGS) -c $<
pddtdriver.o: pddtdriver.f
	$(FC) $(FCFLAGS) -c $<
pcdtdriver.o: pcdtdriver.f
	$(FC) $(FCFLAGS) -c $<
pzdtdriver.o: pzdtdriver.f
	$(FC) $(FCFLAGS) -c $<

psgbdriver.o: psgbdriver.f
	$(FC) $(FCFLAGS) -c $<
pdgbdriver.o: pdgbdriver.f
	$(FC) $(FCFLAGS) -c $<
pcgbdriver.o: pcgbdriver.f
	$(FC) $(FCFLAGS) -c $<
pzgbdriver.o: pzgbdriver.f
	$(FC) $(FCFLAGS) -c $<

pslltdriver.o: pslltdriver.f
	$(FC) $(FCFLAGS) -c $<
pdlltdriver.o: pdlltdriver.f
	$(FC) $(FCFLAGS) -c $<
pclltdriver.o: pclltdriver.f
	$(FC) $(FCFLAGS) -c $<
pzlltdriver.o: pzlltdriver.f
	$(FC) $(FCFLAGS) -c $<

pspbdriver.o: pspbdriver.f
	$(FC) $(FCFLAGS) -c $<
pdpbdriver.o: pdpbdriver.f
	$(FC) $(FCFLAGS) -c $<
pcpbdriver.o: pcpbdriver.f
	$(FC) $(FCFLAGS) -c $<
pzpbdriver.o: pzpbdriver.f
	$(FC) $(FCFLAGS) -c $<

psptdriver.o: psptdriver.f
	$(FC) $(FCFLAGS) -c $<
pdptdriver.o: pdptdriver.f
	$(FC) $(FCFLAGS) -c $<
pcptdriver.o: pcptdriver.f
	$(FC) $(FCFLAGS) -c $<
pzptdriver.o: pzptdriver.f
	$(FC) $(FCFLAGS) -c $<

psinvdriver.o: psinvdriver.f
	$(FC) $(FCFLAGS) -c $<
pdinvdriver.o: pdinvdriver.f
	$(FC) $(FCFLAGS) -c $<
pcinvdriver.o: pcinvdriver.f
	$(FC) $(FCFLAGS) -c $<
pzinvdriver.o: pzinvdriver.f
	$(FC) $(FCFLAGS) -c $<

psqrdriver.o: psqrdriver.f
	$(FC) $(FCFLAGS) -c $<
pdqrdriver.o: pdqrdriver.f
	$(FC) $(FCFLAGS) -c $<
pcqrdriver.o: pcqrdriver.f
	$(FC) $(FCFLAGS) -c $<
pzqrdriver.o: pzqrdriver.f
	$(FC) $(FCFLAGS) -c $<

pslsdriver.o: pslsdriver.f
	$(FC) $(FCFLAGS) -c $<
pdlsdriver.o: pdlsdriver.f
	$(FC) $(FCFLAGS) -c $<
pclsdriver.o: pclsdriver.f
	$(FC) $(FCFLAGS) -c $<
pzlsdriver.o: pzlsdriver.f
	$(FC) $(FCFLAGS) -c $<

.f.o : ; $(FC) -c $(FCFLAGS) $*.f
