############################################################################
#
#  Program:         ScaLAPACK
#
#  Module:          Makefile
#
#  Purpose:         SRC Makefile
#
#  Creation date:   March 20, 1995
#
#  Modified:
#
#  Send bug reports, comments or suggestions to scalapack@cs.utk.edu
#
############################################################################

include ../SLmake.inc

#######################################################################
#  This is the makefile to create a library for ScaLAPACK.
#  The files are organized as follows:
#       ALLAUX -- Auxiliary routines called from all precisions
#       SCLAUX -- Auxiliary routines called from both REAL and COMPLEX
#       DZLAUX -- Auxiliary routines called from both DOUBLE PRECISION
#                 and COMPLEX*16
#       SLASRC -- Single precision real ScaLAPACK routines
#       CLASRC -- Single precision complex ScaLAPACK routines
#       DLASRC -- Double precision real ScaLAPACK routines
#       ZLASRC -- Double precision complex ScaLAPACK routines
#
#  The library can be set up to include routines for any combination
#  of the four precisions.  First, modify the ARCH, ARCHFLAGS, RANLIB,
#  FC and FCFLAGS definitions in ../SLmake.inc to match your library
#  archiver, compiler and the options to be used.
#  Then to create or add to the library, enter make followed by one or
#  more of the precisions desired.  Some examples:
#       make single
#       make single complex
#       make single double complex complex16
#  Alternatively, the command
#       make
#  without any arguments creates a library of all four precisions.
#  The name of the library is defined by ../$(SCALAPACKLIB) in
#  ../SLmake.inc and is created at the next higher directory level.
#
#  To remove the object files after the library is created, enter
#       make clean
#
#######################################################################

ALLAUX = pjlaenv.o pbchkvect.o getpbbuf.o pilaenvx.o piparmq.o pilaver.o pmpim2.o pmpcol.o

SCLAUX = \
   pslabad.o pslaed0.o pslaed1.o pslaed2.o pslaed3.o pslaedz.o pslaiect.o \
   pslamch.o pslared1d.o pslasrt.o psstebz.o psstedc.o slapst.o slasrt2.o \
   sstein2.o \
   slar1va.o slarrb2.o slarrd2.o slarre2.o slarre2a.o slarrf2.o \
   slarrv2.o sstegr2.o sstegr2a.o sstegr2b.o slamov.o clamov.o

DZLAUX = \
   pdlabad.o pdlaed0.o pdlaed1.o pdlaed2.o pdlaed3.o pdlaedz.o pdlaiect.o \
   pdlamch.o pdlared1d.o pdlasrt.o pdstebz.o pdstedc.o dlapst.o dlasrt2.o \
   dstein2.o \
   dlar1va.o dlarrb2.o dlarrd2.o dlarre2.o dlarre2a.o dlarrf2.o \
   dlarrv2.o dstegr2.o dstegr2a.o dstegr2b.o dlamov.o zlamov.o

SLASRC = \
   psdbsv.o  psdbtrf.o psdbtrs.o psdbtrsv.o \
   psdtsv.o  psdttrf.o psdttrs.o psdttrsv.o \
   psgbsv.o  psgbtrf.o psgbtrs.o \
   psgebd2.o psgebrd.o psgecon.o pstrcon.o psgeequ.o psgehd2.o psgehrd.o \
   psgelq2.o psgelqf.o psgels.o  psgeql2.o psgeqlf.o psgeqpf.o psgeqr2.o \
   psgeqrf.o psgerfs.o psgerq2.o psgerqf.o                               \
   psgesv.o  psgesvd.o psgesvx.o psgetf2.o psgetrf.o psgetri.o psgetrs.o \
   psggqrf.o psggrqf.o pslabrd.o \
   pslacon.o pslacp2.o pslacpy.o pslahrd.o pslange.o pslanhs.o pslansy.o \
   pslantr.o pslapiv.o pslapv2.o pslaqge.o pslaqsy.o pslarf.o  pslarfb.o \
   pslarfg.o pslarft.o pslase2.o pslaset.o pslascl.o pslassq.o pslaswp.o \
   pslatra.o pslatrd.o pslatrs.o pslauu2.o pslauum.o psorg2l.o psorg2r.o \
                       psorgl2.o psorglq.o psorgql.o psorgqr.o psorgr2.o \
   psorgrq.o           psorm2l.o psorm2r.o psormbr.o psormhr.o psorml2.o \
   psormlq.o psormql.o psormqr.o psormr2.o psormrq.o psormtr.o pspocon.o \
   pspbsv.o  pspbtrf.o pspbtrs.o pspbtrsv.o \
   psptsv.o  pspttrf.o pspttrs.o pspttrsv.o \
   pspoequ.o psporfs.o psposv.o  psposvx.o pspotf2.o pspotrf.o pspotri.o \
   pspotrs.o psrscl.o  psstein.o pssyev.o  pssyevd.o pssyevx.o \
   pssygs2.o pssygst.o \
   pssygvx.o pssyngst.o pssyntrd.o  pssyttrd.o pssytd2.o pssytrd.o \
   pstrti2.o pstrtri.o pstrtrs.o \
   pslaevswp.o \
   pslarzb.o pslarzt.o pslarz.o pslatrz.o pstzrzf.o psormr3.o psormrz.o \
   pslahqr.o pslaconsb.o pslacp3.o pslawil.o \
   pslasmsub.o pslared2d.o pslamr1d.o slaref.o slamsh.o slasorte.o ssteqr2.o \
   sdbtf2.o  sdbtrf.o  sdttrf.o sdttrsv.o spttrsv.o strmvt.o pssyevr.o \
   bslaapp.o bslaexc.o bstrexc.o pstrord.o pstrsen.o psgebal.o pshseqr.o \
   pslamve.o pslaqr0.o pslaqr1.o pslaqr2.o pslaqr3.o pslaqr4.o pslaqr5.o \
   psrot.o slaqr6.o

CLASRC = \
   pcdbsv.o  pcdbtrf.o pcdbtrs.o pcdbtrsv.o \
   pcdtsv.o  pcdttrf.o pcdttrs.o pcdttrsv.o \
   pcgbsv.o  pcgbtrf.o pcgbtrs.o \
   pcgebd2.o pcgebrd.o pcgecon.o pctrcon.o pcgeequ.o pcgehd2.o pcgehrd.o \
   pcgelq2.o pcgelqf.o  pcgels.o pcgeql2.o pcgeqlf.o pcgeqpf.o pcgeqr2.o \
   pcgeqrf.o pcgerfs.o pcgerq2.o pcgerqf.o                               \
   pcgesv.o  pcgesvd.o pcgesvx.o pcgetf2.o pcgetrf.o pcgetri.o pcgetrs.o \
   pcggqrf.o \
   pcggrqf.o pcheev.o  pcheevd.o pcheevx.o pchegs2.o pchegst.o pchegvx.o \
   pchengst.o pchentrd.o pchettrd.o pchetd2.o \
   pchetrd.o pclabrd.o pclacon.o pclacgv.o pclacp2.o pclacpy.o pclahrd.o \
   pclahqr.o pclaconsb.o pclasmsub.o pclacp3.o pclawil.o pcrot.o \
   pclange.o pclanhe.o pclanhs.o pclansy.o pclantr.o pclapiv.o pclapv2.o \
   pclaqge.o pclaqsy.o pclarf.o  pclarfb.o pclarfc.o pclarfg.o pclarft.o \
   pclascl.o pclase2.o pclaset.o pclassq.o pclaswp.o pclatra.o pclatrd.o \
   pclatrs.o pclauu2.o pclauum.o pcpocon.o pcpoequ.o pcporfs.o pcposv.o  \
   pcpbsv.o  pcpbtrf.o pcpbtrs.o pcpbtrsv.o \
   pcptsv.o  pcpttrf.o pcpttrs.o pcpttrsv.o \
   pcposvx.o pcpotf2.o pcpotrf.o pcpotri.o pcpotrs.o pcsrscl.o pcstein.o \
   pctrevc.o pctrti2.o pctrtri.o pctrtrs.o pcung2l.o pcung2r.o \
             pcungl2.o pcunglq.o pcungql.o pcungqr.o pcungr2.o pcungrq.o \
             pcunm2l.o pcunm2r.o pcunmbr.o pcunmhr.o pcunml2.o pcunmlq.o \
   pcunmql.o pcunmqr.o pcunmr2.o pcunmrq.o pcunmtr.o \
   pclaevswp.o \
   pclarzb.o pclarzt.o pclarz.o  pclarzc.o pclatrz.o pctzrzf.o \
   pclattrs.o \
   pcunmr3.o pcunmrz.o pcmax1.o pscsum1.o pclamr1d.o \
   cdbtf2.o  cdbtrf.o  cdttrf.o cdttrsv.o  cpttrsv.o csteqr2.o ctrmvt.o \
   clamsh.o  claref.o  clanv2.o clahqr2.o pcheevr.o

DLASRC = \
   pddbsv.o  pddbtrf.o pddbtrs.o pddbtrsv.o \
   pddtsv.o  pddttrf.o pddttrs.o pddttrsv.o \
   pdgbsv.o  pdgbtrf.o pdgbtrs.o \
   pdgebd2.o pdgebrd.o pdgecon.o pdtrcon.o pdgeequ.o pdgehd2.o pdgehrd.o \
   pdgelq2.o pdgelqf.o pdgels.o  pdgeql2.o pdgeqlf.o pdgeqpf.o pdgeqr2.o \
   pdgeqrf.o pdgerfs.o pdgerq2.o pdgerqf.o                               \
   pdgesv.o  pdgesvd.o pdgesvx.o pdgetf2.o pdgetrf.o pdgetri.o pdgetrs.o \
   pdggqrf.o pdggrqf.o pdlabrd.o \
   pdlacon.o pdlacp2.o pdlacpy.o pdlahrd.o pdlange.o pdlanhs.o pdlansy.o \
   pdlantr.o pdlapiv.o pdlapv2.o pdlaqge.o pdlaqsy.o pdlarf.o  pdlarfb.o \
   pdlarfg.o pdlarft.o pdlase2.o pdlaset.o pdlascl.o pdlassq.o pdlaswp.o \
   pdlatra.o pdlatrd.o pdlatrs.o pdlauu2.o pdlauum.o pdorg2l.o pdorg2r.o \
                       pdorgl2.o pdorglq.o pdorgql.o pdorgqr.o pdorgr2.o \
   pdorgrq.o           pdorm2l.o pdorm2r.o pdormbr.o pdormhr.o pdorml2.o \
   pdormlq.o pdormql.o pdormqr.o pdormr2.o pdormrq.o pdormtr.o pdpocon.o \
   pdpbsv.o  pdpbtrf.o pdpbtrs.o pdpbtrsv.o \
   pdptsv.o  pdpttrf.o pdpttrs.o pdpttrsv.o \
   pdpoequ.o pdporfs.o pdposv.o  pdposvx.o pdpotf2.o pdpotrf.o pdpotri.o \
   pdpotrs.o pdrscl.o  pdstein.o pdsyev.o  pdsyevd.o pdsyevx.o \
   pdsygs2.o pdsygst.o \
   pdsygvx.o pdsyngst.o pdsyntrd.o pdsyttrd.o pdsytd2.o pdsytrd.o pdtrti2.o \
   pdtrtri.o pdtrtrs.o \
   pdlaevswp.o \
   pdlarzb.o pdlarzt.o pdlarz.o pdlatrz.o pdtzrzf.o pdormr3.o pdormrz.o \
   pdlahqr.o pdlaconsb.o pdlacp3.o pdlawil.o \
   pdlasmsub.o pdlared2d.o pdlamr1d.o dlaref.o dlamsh.o dlasorte.o dsteqr2.o \
   ddbtf2.o  ddbtrf.o  ddttrf.o ddttrsv.o dpttrsv.o dtrmvt.o pdsyevr.o \
   bdlaapp.o bdlaexc.o bdtrexc.o dlaqr6.o pdtrord.o \
   pdtrsen.o pdgebal.o pdhseqr.o pdlamve.o pdlaqr0.o pdlaqr1.o pdlaqr2.o \
   pdlaqr3.o pdlaqr4.o pdlaqr5.o pdrot.o

ZLASRC = \
   pzdbsv.o  pzdbtrf.o pzdbtrs.o pzdbtrsv.o \
   pzdtsv.o  pzdttrf.o pzdttrs.o pzdttrsv.o \
   pzgbsv.o  pzgbtrf.o pzgbtrs.o \
   pzgebd2.o pzgebrd.o pzgecon.o pztrcon.o pzgeequ.o pzgehd2.o pzgehrd.o \
   pzgelq2.o pzgelqf.o pzgels.o  pzgeql2.o pzgeqlf.o pzgeqpf.o pzgeqr2.o \
   pzgeqrf.o pzgerfs.o pzgerq2.o pzgerqf.o                               \
   pzgesv.o  pzgesvd.o pzgesvx.o pzgetf2.o pzgetrf.o pzgetri.o pzgetrs.o \
   pzggqrf.o \
   pzggrqf.o pzheev.o  pzheevd.o pzheevx.o pzhegs2.o pzhegst.o pzhegvx.o \
   pzhengst.o pzhentrd.o pzhettrd.o pzhetd2.o \
   pzhetrd.o pzlabrd.o pzlacon.o pzlacgv.o pzlacp2.o pzlacpy.o pzlahrd.o \
   pzlahqr.o pzlaconsb.o pzlasmsub.o pzlacp3.o pzlawil.o pzrot.o \
   pzlange.o pzlanhe.o pzlanhs.o pzlansy.o pzlantr.o pzlapiv.o pzlapv2.o \
   pzlaqge.o pzlaqsy.o pzlarf.o  pzlarfb.o pzlarfc.o pzlarfg.o pzlarft.o \
   pzlascl.o pzlase2.o pzlaset.o pzlassq.o pzlaswp.o pzlatra.o pzlatrd.o \
   pzlattrs.o \
   pzlatrs.o pzlauu2.o pzlauum.o pzpocon.o pzpoequ.o pzporfs.o pzposv.o  \
   pzpbsv.o  pzpbtrf.o pzpbtrs.o pzpbtrsv.o \
   pzptsv.o  pzpttrf.o pzpttrs.o pzpttrsv.o \
   pzposvx.o pzpotf2.o pzpotrf.o pzpotri.o pzpotrs.o pzdrscl.o pzstein.o \
   pztrevc.o pztrti2.o pztrtri.o pztrtrs.o pzung2l.o pzung2r.o \
             pzungl2.o pzunglq.o pzungql.o pzungqr.o pzungr2.o pzungrq.o \
             pzunm2l.o pzunm2r.o pzunmbr.o pzunmhr.o pzunml2.o pzunmlq.o \
   pzunmql.o pzunmqr.o pzunmr2.o pzunmrq.o pzunmtr.o \
   pzlaevswp.o \
   pzlarzb.o pzlarzt.o pzlarz.o  pzlarzc.o pzlatrz.o pztzrzf.o \
   pzunmr3.o pzunmrz.o pzmax1.o pdzsum1.o pzlamr1d.o \
   zdbtf2.o  zdbtrf.o  zdttrf.o zdttrsv.o  zpttrsv.o zsteqr2.o ztrmvt.o \
   zlamsh.o  zlaref.o  zlanv2.o zlahqr2.o pzheevr.o

all: single complex double complex16

single: $(SLASRC) $(SCLAUX) $(ALLAUX)
	$(ARCH) $(ARCHFLAGS) ../$(SCALAPACKLIB) $(SLASRC) $(SCLAUX) \
	$(ALLAUX)
	$(RANLIB) ../$(SCALAPACKLIB)

complex: $(CLASRC) $(SCLAUX) $(ALLAUX)
	$(ARCH) $(ARCHFLAGS) ../$(SCALAPACKLIB) $(CLASRC) $(SCLAUX) \
	$(ALLAUX)
	$(RANLIB) ../$(SCALAPACKLIB)

double: $(DLASRC) $(DZLAUX) $(ALLAUX)
	$(ARCH) $(ARCHFLAGS) ../$(SCALAPACKLIB) $(DLASRC) $(DZLAUX) \
	$(ALLAUX)
	$(RANLIB) ../$(SCALAPACKLIB)

complex16: $(ZLASRC) $(DZLAUX) $(ALLAUX)
	$(ARCH) $(ARCHFLAGS) ../$(SCALAPACKLIB) $(ZLASRC) $(DZLAUX) \
	$(ALLAUX)
	$(RANLIB) ../$(SCALAPACKLIB)

slamov.o: slamov.c
dlamov.o: dlamov.c
clamov.o: clamov.c
zlamov.o: zlamov.c
slamov.o dlamov.o clamov.o zlamov.o: lamov.h
	$(CC) -c $(CFLAGS) $(CDEFS) $(@:.o=.c) -o $@

clean :
	rm -f *.o

.f.o : 
	$(FC) -c $(FCFLAGS) $*.f

.c.o : 
	$(CC) -c $(CDEFS) $(CCFLAGS) $*.c
