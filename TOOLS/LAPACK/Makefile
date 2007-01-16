############################################################################
#
#  Program:         ScaLAPACK
#
#  Module:          Makefile
#
#  Purpose:         Tools - LAPACK Makefile
#
#  Creation date:   March 20, 1995
#
#  Modified:        February 16, 2000
#
#  Send bug reports, comments or suggestions to scalapack@cs.utk.edu
#
############################################################################

include ../../SLmake.inc

ALLAUX = icopy.o ilaenv.o lsame.o  lsamen.o xerbla.o

SCLAUX = \
   slabad.o slacpy.o sladiv.o slae2.o slaed0.o slaed1.o slaed2.o \
   slaed3.o slaed4.o slaed5.o slaed6.o slaed7.o slaed8.o slaed9.o \
   slaeda.o slagtf.o slagts.o slamch.o slamrg.o \
   slanst.o slapy2.o slapy3.o slarnv.o slartg.o slaruv.o slas2.o \
   slascl.o slaset.o slasq1.o slasq2.o slasq3.o slasq4.o slasq5.o \
   slasq6.o slasrt.o slassq.o slasv2.o ssterf.o ssteqr.o

DZLAUX = \
   dlabad.o dlacpy.o dladiv.o dlae2.o dlaed0.o dlaed1.o dlaed2.o \
   dlaed3.o dlaed4.o dlaed5.o dlaed6.o dlaed7.o dlaed8.o dlaed9.o \
   dlaeda.o dlagtf.o dlagts.o dlamch.o dlamrg.o \
   dlanst.o dlapy2.o dlapy3.o dlarnv.o dlartg.o dlaruv.o dlas2.o \
   dlascl.o dlaset.o dlasq1.o dlasq2.o dlasq3.o dlasq4.o dlasq5.o \
   dlasq6.o dlasrt.o dlassq.o dlasv2.o dsterf.o dsteqr.o

SLASRC = \
   sbdsqr.o sgbtf2.o sgbtrf.o sgetf2.o sgetrf.o sgetrs.o slaev2.o \
   slahqr.o slange.o slanhs.o slanv2.o slarf.o  slarfb.o slarfg.o slarft.o \
   slasr.o  slaswp.o sorm2l.o sorm2r.o sormql.o sormqr.o sormtr.o sstein.o \
   stbtrs.o spbtrf.o spttrf.o strtrs.o spotrf.o spbtf2.o spotf2.o \
   ssytrd.o ssytd2.o sstedc.o slatrd.o 

DLASRC = \
   dbdsqr.o dgbtf2.o dgbtrf.o dgetf2.o dgetrf.o dgetrs.o dlaev2.o \
   dlahqr.o dlange.o dlanhs.o dlanv2.o dlarf.o  dlarfb.o dlarfg.o dlarft.o \
   dlasr.o  dlaswp.o dorm2l.o dorm2r.o dormql.o dormqr.o dormtr.o dstein.o \
   dtbtrs.o dpbtrf.o dpttrf.o dtrtrs.o dpotrf.o dpbtf2.o dpotf2.o  \
   dsytrd.o dsytd2.o dstedc.o dlatrd.o 

CLASRC = \
   cgbtf2.o cgbtrf.o cgetf2.o cgetrf.o cgetrs.o \
   clacgv.o clacpy.o cladiv.o claev2.o clarf.o  clarfb.o clarfg.o clarft.o \
   clartg.o clasr.o claswp.o clascl.o claset.o classq.o cstein.o csymv.o \
   cunm2l.o cunm2r.o cunmql.o cunmqr.o cunmtr.o \
   icmax1.o scsum1.o \
   ctbtrs.o cpbtrf.o cpttrf.o ctrtrs.o cpotrf.o cpbtf2.o cpotf2.o \
   chetrd.o chetd2.o clatrd.o cstedc.o clanhs.o crot.o

ZLASRC = \
   zgbtf2.o zgbtrf.o zgetf2.o zgetrf.o zgetrs.o \
   zlacgv.o zlacpy.o zladiv.o zlaev2.o zlarf.o  zlarfb.o zlarfg.o zlarft.o \
   zlartg.o zlasr.o  zlaswp.o zlascl.o zlaset.o zlassq.o zstein.o zsymv.o \
   zunm2l.o zunm2r.o zunmql.o zunmqr.o zunmtr.o \
   izmax1.o dzsum1.o \
   ztbtrs.o zpbtrf.o zpttrf.o ztrtrs.o zpotrf.o zpbtf2.o zpotf2.o \
   zhetrd.o zhetd2.o zlatrd.o zstedc.o zlanhs.o zrot.o \
   cbdsqr.o zbdsqr.o

SCATGEN = slatm1.o slaran.o slarnd.o

SMATGEN = slatms.o slagge.o slagsy.o slarot.o

CMATGEN = clarnv.o clatm1.o clatms.o clagge.o claghe.o clagsy.o clarot.o \
          clarnd.o

DZATGEN = dlatm1.o dlaran.o dlarnd.o

DMATGEN = dlatms.o dlagge.o dlagsy.o dlarot.o

ZMATGEN = zlarnv.o zlatm1.o zlatms.o zlagge.o zlaghe.o zlagsy.o zlarot.o \
          zlarnd.o

all : single complex double complex16

single: $(SLASRC) $(ALLAUX) $(SCLAUX) $(SMATGEN) $(SCATGEN)
	$(ARCH) $(ARCHFLAGS) $(SCALAPACKLIB) $(SLASRC) $(ALLAUX) $(SCLAUX) \
	$(SMATGEN) $(SCATGEN)
	$(RANLIB) $(SCALAPACKLIB)

complex: $(CLASRC) $(ALLAUX) $(SCLAUX) $(CMATGEN) $(SCATGEN)
	$(ARCH) $(ARCHFLAGS) $(SCALAPACKLIB) $(CLASRC) $(ALLAUX) $(SCLAUX) \
	$(CMATGEN) $(SCATGEN)
	$(RANLIB) $(SCALAPACKLIB)

double: $(DLASRC) $(ALLAUX) $(DZLAUX) $(DMATGEN) $(DZATGEN)
	$(ARCH) $(ARCHFLAGS) $(SCALAPACKLIB) $(DLASRC) $(ALLAUX) $(DZLAUX) \
	$(DMATGEN) $(DZATGEN)
	$(RANLIB) $(SCALAPACKLIB)

complex16: $(ZLASRC) $(ALLAUX) $(DZLAUX) $(ZMATGEN) $(DZATGEN)
	$(ARCH) $(ARCHFLAGS) $(SCALAPACKLIB) $(ZLASRC) $(ALLAUX) $(DZLAUX) \
	$(ZMATGEN) $(DZATGEN)
	$(RANLIB) $(SCALAPACKLIB)

$(ALLAUX): $(FRC)
$(SCLAUX): $(FRC)
$(DZLAUX): $(FRC)
$(SLASRC): $(FRC)
$(CLASRC): $(FRC)
$(DLASRC): $(FRC)
$(ZLASRC): $(FRC)

$(SCATGEN): $(FRC)
$(SMATGEN): $(FRC)
$(CMATGEN): $(FRC)
$(DZATGEN): $(FRC)
$(DMATGEN): $(FRC)
$(ZMATGEN): $(FRC)

FRC:
	@FRC=$(FRC)

clean :
	rm -f *.o

slamch.o:
	$(F77) -c $(NOOPT) slamch.f

dlamch.o:
	$(F77) -c $(NOOPT) dlamch.f

.f.o : ; $(F77) -c $(F77FLAGS) $*.f

