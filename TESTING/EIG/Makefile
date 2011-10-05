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
#  Modified:
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

all : single double complex complex16

single:    $(shrdexe) $(strdexe) $(sbrdexe) $(ssepexe) $(sgsepexe) $(snepexe) $(ssvdexe) $(ssyevrexe)
double:    $(dhrdexe) $(dtrdexe) $(dbrdexe) $(dsepexe) $(dgsepexe) $(dnepexe) $(dsvdexe) $(dsyevrexe)
complex:   $(chrdexe) $(ctrdexe) $(cbrdexe) $(csepexe) $(cgsepexe) $(cnepexe) $(cevcexe) $(cheevrexe)
complex16: $(zhrdexe) $(ztrdexe) $(zbrdexe) $(zsepexe) $(zgsepexe) $(znepexe) $(zevcexe) $(zheevrexe)

$(sbrdexe) : ../../$(SCALAPACKLIB) $(sbrd)
	$(F77LOADER) $(F77LOADFLAGS) -o $(sbrdexe) $(sbrd) ../../$(SCALAPACKLIB) $(LIBS)

$(dbrdexe) : ../../$(SCALAPACKLIB) $(dbrd)
	$(F77LOADER) $(F77LOADFLAGS) -o $(dbrdexe) $(dbrd) ../../$(SCALAPACKLIB) $(LIBS)

$(cbrdexe) : ../../$(SCALAPACKLIB) $(cbrd)
	$(F77LOADER) $(F77LOADFLAGS) -o $(cbrdexe) $(cbrd) ../../$(SCALAPACKLIB) $(LIBS)

$(zbrdexe) : ../../$(SCALAPACKLIB) $(zbrd)
	$(F77LOADER) $(F77LOADFLAGS) -o $(zbrdexe) $(zbrd) ../../$(SCALAPACKLIB) $(LIBS)

$(shrdexe) : ../../$(SCALAPACKLIB) $(shrd)
	$(F77LOADER) $(F77LOADFLAGS) -o $(shrdexe) $(shrd) ../../$(SCALAPACKLIB) $(LIBS)

$(dhrdexe) : ../../$(SCALAPACKLIB) $(dhrd)
	$(F77LOADER) $(F77LOADFLAGS) -o $(dhrdexe) $(dhrd) ../../$(SCALAPACKLIB) $(LIBS)

$(chrdexe) : ../../$(SCALAPACKLIB) $(chrd)
	$(F77LOADER) $(F77LOADFLAGS) -o $(chrdexe) $(chrd) ../../$(SCALAPACKLIB) $(LIBS)

$(zhrdexe) : ../../$(SCALAPACKLIB) $(zhrd)
	$(F77LOADER) $(F77LOADFLAGS) -o $(zhrdexe) $(zhrd) ../../$(SCALAPACKLIB) $(LIBS)

$(strdexe) : ../../$(SCALAPACKLIB) $(strd)
	$(F77LOADER) $(F77LOADFLAGS) -o $(strdexe) $(strd) ../../$(SCALAPACKLIB) $(LIBS)

$(dtrdexe) : ../../$(SCALAPACKLIB) $(dtrd)
	$(F77LOADER) $(F77LOADFLAGS) -o $(dtrdexe) $(dtrd) ../../$(SCALAPACKLIB) $(LIBS)

$(ctrdexe) : ../../$(SCALAPACKLIB) $(ctrd)
	$(F77LOADER) $(F77LOADFLAGS) -o $(ctrdexe) $(ctrd) ../../$(SCALAPACKLIB) $(LIBS)

$(ztrdexe) : ../../$(SCALAPACKLIB) $(ztrd)
	$(F77LOADER) $(F77LOADFLAGS) -o $(ztrdexe) $(ztrd) ../../$(SCALAPACKLIB) $(LIBS)

$(ssvdexe) : ../../$(SCALAPACKLIB) $(ssvd)
	$(F77LOADER) $(F77LOADFLAGS) -o $(ssvdexe) $(ssvd) ../../$(SCALAPACKLIB) $(LIBS)

$(dsvdexe) : ../../$(SCALAPACKLIB) $(dsvd)
	$(F77LOADER) $(F77LOADFLAGS) -o $(dsvdexe) $(dsvd) ../../$(SCALAPACKLIB) $(LIBS)

$(ssepexe) : ../../$(SCALAPACKLIB) $(ssep)
	$(F77LOADER) $(F77LOADFLAGS) -o $(ssepexe) $(ssep) ../../$(SCALAPACKLIB) $(LIBS)

$(dsepexe) : ../../$(SCALAPACKLIB) $(dsep)
	$(F77LOADER) $(F77LOADFLAGS) -o $(dsepexe) $(dsep) ../../$(SCALAPACKLIB) $(LIBS)

$(csepexe) : ../../$(SCALAPACKLIB) $(csep)
	$(F77LOADER) $(F77LOADFLAGS) -o $(csepexe) $(csep) ../../$(SCALAPACKLIB) $(LIBS)

$(zsepexe) : ../../$(SCALAPACKLIB) $(zsep)
	$(F77LOADER) $(F77LOADFLAGS) -o $(zsepexe) $(zsep) ../../$(SCALAPACKLIB) $(LIBS)

$(sgsepexe) : ../../$(SCALAPACKLIB) $(sgsep)
	$(F77LOADER) $(F77LOADFLAGS) -o $(sgsepexe) $(sgsep) ../../$(SCALAPACKLIB) $(LIBS)

$(dgsepexe) : ../../$(SCALAPACKLIB) $(dgsep)
	$(F77LOADER) $(F77LOADFLAGS) -o $(dgsepexe) $(dgsep) ../../$(SCALAPACKLIB) $(LIBS)

$(cgsepexe) : ../../$(SCALAPACKLIB) $(cgsep)
	$(F77LOADER) $(F77LOADFLAGS) -o $(cgsepexe) $(cgsep) ../../$(SCALAPACKLIB) $(LIBS)

$(zgsepexe) : ../../$(SCALAPACKLIB) $(zgsep)
	$(F77LOADER) $(F77LOADFLAGS) -o $(zgsepexe) $(zgsep) ../../$(SCALAPACKLIB) $(LIBS)

$(snepexe) : ../../$(SCALAPACKLIB) $(snep)
	$(F77LOADER) $(F77LOADFLAGS) -o $(snepexe) $(snep) ../../$(SCALAPACKLIB) $(LIBS)

$(dnepexe) : ../../$(SCALAPACKLIB) $(dnep)
	$(F77LOADER) $(F77LOADFLAGS) -o $(dnepexe) $(dnep) ../../$(SCALAPACKLIB) $(LIBS)

$(cnepexe) : ../../$(SCALAPACKLIB) $(cnep)
	$(F77LOADER) $(F77LOADFLAGS) -o $(cnepexe) $(cnep) ../../$(SCALAPACKLIB) $(LIBS)

$(znepexe) : ../../$(SCALAPACKLIB) $(znep)
	$(F77LOADER) $(F77LOADFLAGS) -o $(znepexe) $(znep) ../../$(SCALAPACKLIB) $(LIBS)

$(cevcexe) : ../../$(SCALAPACKLIB) $(cevc)
	$(F77LOADER) $(F77LOADFLAGS) -o $(cevcexe) $(cevc) ../../$(SCALAPACKLIB) $(LIBS)

$(zevcexe) : ../../$(SCALAPACKLIB) $(zevc)
	$(F77LOADER) $(F77LOADFLAGS) -o $(zevcexe) $(zevc) ../../$(SCALAPACKLIB) $(LIBS)

$(ssyevrexe): ../../$(SCALAPACKLIB) $(ssyevr)
	$(F77LOADER) $(F77LOADFLAGS) -o $(ssyevrexe) $(ssyevr) ../../$(SCALAPACKLIB) $(LIBS)

$(dsyevrexe): ../../$(SCALAPACKLIB) $(dsyevr)
	$(F77LOADER) $(F77LOADFLAGS) -o $(dsyevrexe) $(dsyevr) ../../$(SCALAPACKLIB) $(LIBS)

$(cheevrexe): ../../$(SCALAPACKLIB) $(cheevr)
	$(F77LOADER) $(F77LOADFLAGS) -o $(cheevrexe) $(cheevr) ../../$(SCALAPACKLIB) $(LIBS)

$(zheevrexe): ../../$(SCALAPACKLIB) $(zheevr)
	$(F77LOADER) $(F77LOADFLAGS) -o $(zheevrexe) $(zheevr) ../../$(SCALAPACKLIB) $(LIBS)

$(sbrd): $(FRC)
$(dbrd): $(FRC)
$(cbrd): $(FRC)
$(zbrd): $(FRC)

$(shrd): $(FRC)
$(dhrd): $(FRC)
$(chrd): $(FRC)
$(zhrd): $(FRC)

$(strd): $(FRC)
$(dtrd): $(FRC)
$(ctrd): $(FRC)
$(ztrd): $(FRC)

$(ssvd): $(FRC)
$(dsvd): $(FRC)

$(ssep): $(FRC)
$(dsep): $(FRC)
$(csep): $(FRC)
$(zsep): $(FRC)

$(sgsep): $(FRC)
$(dgsep): $(FRC)
$(cgsep): $(FRC)
$(zgsep): $(FRC)

$(snep): $(FRC)
$(dnep): $(FRC)
$(cnep): $(FRC)
$(znep): $(FRC)

$(cevc): $(FRC)
$(zevc): $(FRC)

FRC:
	@FRC=$(FRC)

clean :
	rm -f *.o

psbrddriver.o: psbrddriver.f
	$(F77) $(DRVOPTS) -c $<
pdbrddriver.o: pdbrddriver.f
	$(F77) $(DRVOPTS) -c $<
pcbrddriver.o: pcbrddriver.f
	$(F77) $(DRVOPTS) -c $<
pzbrddriver.o: pzbrddriver.f
	$(F77) $(DRVOPTS) -c $<

pshrddriver.o: pshrddriver.f
	$(F77) $(DRVOPTS) -c $<
pdhrddriver.o: pdhrddriver.f
	$(F77) $(DRVOPTS) -c $<
pchrddriver.o: pchrddriver.f
	$(F77) $(DRVOPTS) -c $<
pzhrddriver.o: pzhrddriver.f
	$(F77) $(DRVOPTS) -c $<

pstrddriver.o: pstrddriver.f
	$(F77) $(DRVOPTS) -c $<
pdtrddriver.o: pdtrddriver.f
	$(F77) $(DRVOPTS) -c $<
pctrddriver.o: pctrddriver.f
	$(F77) $(DRVOPTS) -c $<
pztrddriver.o: pztrddriver.f
	$(F77) $(DRVOPTS) -c $<

pssvddriver.o: pssvddriver.f
	$(F77) $(DRVOPTS) -c $<
pdsvddriver.o: pdsvddriver.f
	$(F77) $(DRVOPTS) -c $<

pssepdriver.o: pssepdriver.f
	$(F77) $(DRVOPTS) -c $<
pdsepdriver.o: pdsepdriver.f
	$(F77) $(DRVOPTS) -c $<
pcsepdriver.o: pcsepdriver.f
	$(F77) $(DRVOPTS) -c $<
pzsepdriver.o: pzsepdriver.f
	$(F77) $(DRVOPTS) -c $<

psgsepdriver.o: psgsepdriver.f
	$(F77) $(DRVOPTS) -c $<
pdgsepdriver.o: pdgsepdriver.f
	$(F77) $(DRVOPTS) -c $<
pcgsepdriver.o: pcgsepdriver.f
	$(F77) $(DRVOPTS) -c $<
pzgsepdriver.o: pzgsepdriver.f
	$(F77) $(DRVOPTS) -c $<

psnepdriver.o: psnepdriver.f
	$(F77) $(DRVOPTS) -c $<
pdnepdriver.o: pdnepdriver.f
	$(F77) $(DRVOPTS) -c $<
pcnepdriver.o: pcnepdriver.f
	$(F77) $(DRVOPTS) -c $<
pznepdriver.o: pznepdriver.f
	$(F77) $(DRVOPTS) -c $<

pcevcdriver.o: pcevcdriver.f
	$(F77) $(DRVOPTS) -c $<
pzevcdriver.o: pzevcdriver.f
	$(F77) $(DRVOPTS) -c $<

.f.o : ; $(F77) -c $(F77FLAGS) $*.f
