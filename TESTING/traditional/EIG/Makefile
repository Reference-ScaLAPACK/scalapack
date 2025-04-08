############################################################################
#
#  Program:         ScaLAPACK
#
#  Module:          Makefile
#
#  Purpose:         Eigenroutine Testing Makefile
#
#  Creation date:   February 20, 2000
#
#  Modified:        November, 2011
#
#  Send bug reports, comments or suggestions to scalapack@cs.utk.edu
#
############################################################################

include ../../SLmake.inc

shrdexe = ../xshrd
dhrdexe = ../xdhrd
chrdexe = ../xchrd
zhrdexe = ../xzhrd

strdexe = ../xstrd
dtrdexe = ../xdtrd
ctrdexe = ../xctrd
ztrdexe = ../xztrd

sbrdexe = ../xsbrd
dbrdexe = ../xdbrd
cbrdexe = ../xcbrd
zbrdexe = ../xzbrd


ssepexe = ../xssep
dsepexe = ../xdsep
csepexe = ../xcsep
zsepexe = ../xzsep

sgsepexe = ../xsgsep
dgsepexe = ../xdgsep
cgsepexe = ../xcgsep
zgsepexe = ../xzgsep

ssvdexe = ../xssvd
dsvdexe = ../xdsvd

snepexe = ../xsnep
dnepexe = ../xdnep
cnepexe = ../xcnep
znepexe = ../xznep

cevcexe = ../xcevc
zevcexe = ../xzevc

ssyevrexe = ../xssyevr
dsyevrexe = ../xdsyevr
cheevrexe = ../xcheevr
zheevrexe = ../xzheevr

shseqrexe = ../xshseqr
dhseqrexe = ../xdhseqr

smatgen = psmatgen.o pmatgeninc.o
dmatgen = pdmatgen.o pmatgeninc.o
cmatgen = pcmatgen.o pmatgeninc.o
zmatgen = pzmatgen.o pmatgeninc.o

sbrd = psbrddriver.o psbrdinfo.o psgebdrv.o pslafchk.o $(smatgen)
dbrd = pdbrddriver.o pdbrdinfo.o pdgebdrv.o pdlafchk.o $(dmatgen)
cbrd = pcbrddriver.o pcbrdinfo.o pcgebdrv.o pclafchk.o $(cmatgen)
zbrd = pzbrddriver.o pzbrdinfo.o pzgebdrv.o pzlafchk.o $(zmatgen)

shrd = pshrddriver.o pshrdinfo.o psgehdrv.o pslafchk.o $(smatgen)
dhrd = pdhrddriver.o pdhrdinfo.o pdgehdrv.o pdlafchk.o $(dmatgen)
chrd = pchrddriver.o pchrdinfo.o pcgehdrv.o pclafchk.o $(cmatgen)
zhrd = pzhrddriver.o pzhrdinfo.o pzgehdrv.o pzlafchk.o $(zmatgen)

strd = pstrddriver.o psttrdtester.o pslatran.o pstrdinfo.o pssytdrv.o \
       pslafchk.o xpjlaenv.o $(smatgen)
dtrd = pdtrddriver.o pdttrdtester.o pdlatran.o pdtrdinfo.o pdsytdrv.o \
       pdlafchk.o xpjlaenv.o $(dmatgen)
ctrd = pctrddriver.o pcttrdtester.o pclatran.o pctrdinfo.o pchetdrv.o \
       pclafchk.o xpjlaenv.o $(cmatgen)
ztrd = pztrddriver.o pzttrdtester.o pzlatran.o pztrdinfo.o pzhetdrv.o \
       pzlafchk.o xpjlaenv.o $(zmatgen)

ssvd =  pssvddriver.o pslagge.o pssvdchk.o pssvdcmp.o pssvdtst.o $(smatgen)
dsvd =  pdsvddriver.o pdlagge.o pdsvdchk.o pdsvdcmp.o pdsvdtst.o $(dmatgen)

ssep =  psseptst.o pssepsubtst.o pssepchk.o pssepqtq.o pslagsy.o \
        pslatms.o  pslasizesep.o pslasizesyevx.o pssepdriver.o \
        pssepreq.o pssepinfo.o pslasizesyev.o \
        pssqpsubtst.o pslasizesqp.o pssdpsubtst.o $(smatgen)
dsep =  pdseptst.o pdsepsubtst.o pdsepchk.o pdsepqtq.o pdlagsy.o \
        pdlatms.o  pdlasizesep.o pdlasizesyevx.o pdsepdriver.o \
        pdsepreq.o pdsepinfo.o pdlasizesyev.o \
        pdsqpsubtst.o pdlasizesqp.o pdsdpsubtst.o $(dmatgen)
csep =  pcseptst.o pcsepsubtst.o pcsepchk.o pcsepqtq.o pclagsy.o \
        pclatms.o  pclasizesep.o pclasizeheevx.o pcsepdriver.o \
        pcsepreq.o pssepinfo.o pcsdpsubtst.o $(cmatgen)
zsep =  pzseptst.o pzsepsubtst.o pzsepchk.o pzsepqtq.o pzlagsy.o \
        pzlatms.o  pzlasizesep.o pzlasizeheevx.o pzsepdriver.o \
        pzsepreq.o pdsepinfo.o pzsdpsubtst.o $(zmatgen)

sgsep =  psgseptst.o psgsepsubtst.o psgsepchk.o pslagsy.o \
         pslatms.o pslasizesyevx.o pslasizegsep.o pslasizesep.o \
         psgsepdriver.o psgsepreq.o pssepinfo.o $(smatgen)
dgsep =  pdgseptst.o pdgsepsubtst.o pdgsepchk.o pdlagsy.o \
         pdlatms.o pdlasizesyevx.o pdlasizegsep.o pdlasizesep.o \
         pdgsepdriver.o pdgsepreq.o pdsepinfo.o $(dmatgen)
cgsep =  pcgseptst.o pcgsepsubtst.o pcgsepchk.o pclagsy.o \
         pclatms.o pclasizegsep.o pclasizeheevx.o pclasizesep.o \
         pcgsepdriver.o pcgsepreq.o pssepinfo.o $(cmatgen)
zgsep =  pzgseptst.o pzgsepsubtst.o pzgsepchk.o pzlagsy.o \
         pzlatms.o pzlasizegsep.o pzlasizeheevx.o pzlasizesep.o \
         pzgsepdriver.o pzgsepreq.o pdsepinfo.o $(zmatgen)

snep =  psnepdriver.o psnepinfo.o psnepfchk.o $(smatgen)
dnep =  pdnepdriver.o pdnepinfo.o pdnepfchk.o $(dmatgen)
cnep =  pcnepdriver.o pcnepinfo.o pcnepfchk.o $(cmatgen)
znep =  pznepdriver.o pznepinfo.o pznepfchk.o $(zmatgen)

cevc =  pcevcdriver.o pcevcinfo.o pcget22.o $(cmatgen)
zevc =  pzevcdriver.o pzevcinfo.o pzget22.o $(zmatgen)

ssyevr =  pslasizesepr.o pslasizesyevr.o psseprdriver.o psseprreq.o psseprsubtst.o \
pssepchk.o pssepqtq.o pslatms.o psseprtst.o pssepinfo.o pslagsy.o pslasizesep.o $(smatgen)
dsyevr =  pdlasizesepr.o pdlasizesyevr.o pdseprdriver.o pdseprreq.o pdseprsubtst.o \
pdsepchk.o pdsepqtq.o pdlatms.o pdseprtst.o pdsepinfo.o pdlagsy.o pdlasizesep.o $(dmatgen) 
cheevr =  pclasizesepr.o pclasizeheevr.o pcseprdriver.o pcseprreq.o pcseprsubtst.o \
pcsepchk.o pcsepqtq.o pclatms.o pcseprtst.o pssepinfo.o pclagsy.o pclasizesep.o $(cmatgen)
zheevr =  pzlasizesepr.o pzlasizeheevr.o pzseprdriver.o pzseprreq.o pzseprsubtst.o \
pzsepchk.o pzsepqtq.o pzlatms.o pzseprtst.o pdsepinfo.o pzlagsy.o pzlasizesep.o $(zmatgen)

shseqr = pshseqrdriver.o psmatgen2.o $(smatgen)
dhseqr = pdhseqrdriver.o pdmatgen2.o $(dmatgen)

all : single double complex complex16

single:    $(shrdexe) $(strdexe) $(sbrdexe) $(ssepexe) $(sgsepexe) $(snepexe) $(ssvdexe) $(ssyevrexe) $(shseqrexe)
double:    $(dhrdexe) $(dtrdexe) $(dbrdexe) $(dsepexe) $(dgsepexe) $(dnepexe) $(dsvdexe) $(dsyevrexe) $(dhseqrexe)
complex:   $(chrdexe) $(ctrdexe) $(cbrdexe) $(csepexe) $(cgsepexe) $(cnepexe) $(cevcexe) $(cheevrexe)
complex16: $(zhrdexe) $(ztrdexe) $(zbrdexe) $(zsepexe) $(zgsepexe) $(znepexe) $(zevcexe) $(zheevrexe)

$(sbrdexe) : ../../$(SCALAPACKLIB) $(sbrd)
	$(FCLOADER) $(FCLOADFLAGS) -o $(sbrdexe) $(sbrd) ../../$(SCALAPACKLIB) $(LIBS)

$(dbrdexe) : ../../$(SCALAPACKLIB) $(dbrd)
	$(FCLOADER) $(FCLOADFLAGS) -o $(dbrdexe) $(dbrd) ../../$(SCALAPACKLIB) $(LIBS)

$(cbrdexe) : ../../$(SCALAPACKLIB) $(cbrd)
	$(FCLOADER) $(FCLOADFLAGS) -o $(cbrdexe) $(cbrd) ../../$(SCALAPACKLIB) $(LIBS)

$(zbrdexe) : ../../$(SCALAPACKLIB) $(zbrd)
	$(FCLOADER) $(FCLOADFLAGS) -o $(zbrdexe) $(zbrd) ../../$(SCALAPACKLIB) $(LIBS)

$(shrdexe) : ../../$(SCALAPACKLIB) $(shrd)
	$(FCLOADER) $(FCLOADFLAGS) -o $(shrdexe) $(shrd) ../../$(SCALAPACKLIB) $(LIBS)

$(dhrdexe) : ../../$(SCALAPACKLIB) $(dhrd)
	$(FCLOADER) $(FCLOADFLAGS) -o $(dhrdexe) $(dhrd) ../../$(SCALAPACKLIB) $(LIBS)

$(chrdexe) : ../../$(SCALAPACKLIB) $(chrd)
	$(FCLOADER) $(FCLOADFLAGS) -o $(chrdexe) $(chrd) ../../$(SCALAPACKLIB) $(LIBS)

$(zhrdexe) : ../../$(SCALAPACKLIB) $(zhrd)
	$(FCLOADER) $(FCLOADFLAGS) -o $(zhrdexe) $(zhrd) ../../$(SCALAPACKLIB) $(LIBS)

$(strdexe) : ../../$(SCALAPACKLIB) $(strd)
	$(FCLOADER) $(FCLOADFLAGS) -o $(strdexe) $(strd) ../../$(SCALAPACKLIB) $(LIBS)

$(dtrdexe) : ../../$(SCALAPACKLIB) $(dtrd)
	$(FCLOADER) $(FCLOADFLAGS) -o $(dtrdexe) $(dtrd) ../../$(SCALAPACKLIB) $(LIBS)

$(ctrdexe) : ../../$(SCALAPACKLIB) $(ctrd)
	$(FCLOADER) $(FCLOADFLAGS) -o $(ctrdexe) $(ctrd) ../../$(SCALAPACKLIB) $(LIBS)

$(ztrdexe) : ../../$(SCALAPACKLIB) $(ztrd)
	$(FCLOADER) $(FCLOADFLAGS) -o $(ztrdexe) $(ztrd) ../../$(SCALAPACKLIB) $(LIBS)

$(ssvdexe) : ../../$(SCALAPACKLIB) $(ssvd)
	$(FCLOADER) $(FCLOADFLAGS) -o $(ssvdexe) $(ssvd) ../../$(SCALAPACKLIB) $(LIBS)

$(dsvdexe) : ../../$(SCALAPACKLIB) $(dsvd)
	$(FCLOADER) $(FCLOADFLAGS) -o $(dsvdexe) $(dsvd) ../../$(SCALAPACKLIB) $(LIBS)

$(ssepexe) : ../../$(SCALAPACKLIB) $(ssep)
	$(FCLOADER) $(FCLOADFLAGS) -o $(ssepexe) $(ssep) ../../$(SCALAPACKLIB) $(LIBS)

$(dsepexe) : ../../$(SCALAPACKLIB) $(dsep)
	$(FCLOADER) $(FCLOADFLAGS) -o $(dsepexe) $(dsep) ../../$(SCALAPACKLIB) $(LIBS)

$(csepexe) : ../../$(SCALAPACKLIB) $(csep)
	$(FCLOADER) $(FCLOADFLAGS) -o $(csepexe) $(csep) ../../$(SCALAPACKLIB) $(LIBS)

$(zsepexe) : ../../$(SCALAPACKLIB) $(zsep)
	$(FCLOADER) $(FCLOADFLAGS) -o $(zsepexe) $(zsep) ../../$(SCALAPACKLIB) $(LIBS)

$(sgsepexe) : ../../$(SCALAPACKLIB) $(sgsep)
	$(FCLOADER) $(FCLOADFLAGS) -o $(sgsepexe) $(sgsep) ../../$(SCALAPACKLIB) $(LIBS)

$(dgsepexe) : ../../$(SCALAPACKLIB) $(dgsep)
	$(FCLOADER) $(FCLOADFLAGS) -o $(dgsepexe) $(dgsep) ../../$(SCALAPACKLIB) $(LIBS)

$(cgsepexe) : ../../$(SCALAPACKLIB) $(cgsep)
	$(FCLOADER) $(FCLOADFLAGS) -o $(cgsepexe) $(cgsep) ../../$(SCALAPACKLIB) $(LIBS)

$(zgsepexe) : ../../$(SCALAPACKLIB) $(zgsep)
	$(FCLOADER) $(FCLOADFLAGS) -o $(zgsepexe) $(zgsep) ../../$(SCALAPACKLIB) $(LIBS)

$(snepexe) : ../../$(SCALAPACKLIB) $(snep)
	$(FCLOADER) $(FCLOADFLAGS) -o $(snepexe) $(snep) ../../$(SCALAPACKLIB) $(LIBS)

$(dnepexe) : ../../$(SCALAPACKLIB) $(dnep)
	$(FCLOADER) $(FCLOADFLAGS) -o $(dnepexe) $(dnep) ../../$(SCALAPACKLIB) $(LIBS)

$(cnepexe) : ../../$(SCALAPACKLIB) $(cnep)
	$(FCLOADER) $(FCLOADFLAGS) -o $(cnepexe) $(cnep) ../../$(SCALAPACKLIB) $(LIBS)

$(znepexe) : ../../$(SCALAPACKLIB) $(znep)
	$(FCLOADER) $(FCLOADFLAGS) -o $(znepexe) $(znep) ../../$(SCALAPACKLIB) $(LIBS)

$(cevcexe): ../../$(SCALAPACKLIB) $(cevc)
	$(FCLOADER) $(FCLOADFLAGS) -o $(cevcexe) $(cevc) ../../$(SCALAPACKLIB) $(LIBS)

$(zevcexe): ../../$(SCALAPACKLIB) $(zevc)
	$(FCLOADER) $(FCLOADFLAGS) -o $(zevcexe) $(zevc) ../../$(SCALAPACKLIB) $(LIBS)

$(ssyevrexe): ../../$(SCALAPACKLIB) $(ssyevr)
	$(FCLOADER) $(FCLOADFLAGS) -o $(ssyevrexe) $(ssyevr) ../../$(SCALAPACKLIB) $(LIBS)

$(dsyevrexe): ../../$(SCALAPACKLIB) $(dsyevr)
	$(FCLOADER) $(FCLOADFLAGS) -o $(dsyevrexe) $(dsyevr) ../../$(SCALAPACKLIB) $(LIBS)

$(cheevrexe): ../../$(SCALAPACKLIB) $(cheevr)
	$(FCLOADER) $(FCLOADFLAGS) -o $(cheevrexe) $(cheevr) ../../$(SCALAPACKLIB) $(LIBS)

$(zheevrexe): ../../$(SCALAPACKLIB) $(zheevr)
	$(FCLOADER) $(FCLOADFLAGS) -o $(zheevrexe) $(zheevr) ../../$(SCALAPACKLIB) $(LIBS)

$(shseqrexe): ../../$(SCALAPACKLIB) $(shseqr)
	$(FCLOADER) $(FCLOADFLAGS) -o $(shseqrexe) $(shseqr) ../../$(SCALAPACKLIB) $(LIBS)

$(dhseqrexe): ../../$(SCALAPACKLIB) $(dhseqr)
	$(FCLOADER) $(FCLOADFLAGS) -o $(dhseqrexe) $(dhseqr) ../../$(SCALAPACKLIB) $(LIBS)

clean:
	rm -f *.o

.f.o:
	$(FC) -c $(FCFLAGS) $<
