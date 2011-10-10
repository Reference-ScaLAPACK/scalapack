dlvl = ../..
include $(dlvl)/SLmake.inc

#  --------------------------
#  The communication routines
#  --------------------------
comm = igesd2d_.o sgesd2d_.o dgesd2d_.o cgesd2d_.o zgesd2d_.o \
       itrsd2d_.o strsd2d_.o dtrsd2d_.o ctrsd2d_.o ztrsd2d_.o \
       igerv2d_.o sgerv2d_.o dgerv2d_.o cgerv2d_.o zgerv2d_.o \
       itrrv2d_.o strrv2d_.o dtrrv2d_.o ctrrv2d_.o ztrrv2d_.o \
       igebs2d_.o sgebs2d_.o dgebs2d_.o cgebs2d_.o zgebs2d_.o \
       igebr2d_.o sgebr2d_.o dgebr2d_.o cgebr2d_.o zgebr2d_.o \
       itrbs2d_.o strbs2d_.o dtrbs2d_.o ctrbs2d_.o ztrbs2d_.o \
       itrbr2d_.o strbr2d_.o dtrbr2d_.o ctrbr2d_.o ztrbr2d_.o \
       igsum2d_.o sgsum2d_.o dgsum2d_.o cgsum2d_.o zgsum2d_.o \
       igamx2d_.o sgamx2d_.o dgamx2d_.o cgamx2d_.o zgamx2d_.o \
       igamn2d_.o sgamn2d_.o dgamn2d_.o cgamn2d_.o zgamn2d_.o

#  --------------------
#  The support routines
#  --------------------
supp = blacs_setup_.o blacs_set_.o blacs_get_.o \
       blacs_abort_.o blacs_exit_.o blacs_pnum_.o blacs_pcoord_.o \
       ksendid_.o krecvid_.o kbsid_.o kbrid_.o \
       dcputime00_.o dwalltime00_.o blacs_pinfo_.o \
       blacs_init_.o blacs_map_.o blacs_free_.o blacs_grid_.o blacs_info_.o \
       blacs_barr_.o sys2blacs_.o blacs2sys_.o free_handle_.o


#  ----------------------------
#  The fortran and C interfaces
#  ----------------------------
Fintobj   = $(comm) $(supp)
Cintobj   = $(comm:.o=.oo) $(supp:.o=.oo)

#  ---------------------
#  The internal routines
#  ---------------------
internal = BI_HypBS.o BI_HypBR.o BI_IdringBS.o BI_IdringBR.o \
           BI_MpathBS.o BI_MpathBR.o BI_SringBS.o BI_SringBR.o \
           BI_TreeBS.o BI_TreeBR.o \
           BI_Ssend.o BI_Rsend.o BI_Srecv.o BI_Asend.o BI_Arecv.o \
           BI_TreeComb.o BI_BeComb.o BI_MringComb.o \
           BI_ArgCheck.o BI_TransDist.o BI_GetBuff.o BI_UpdateBuffs.o \
           BI_EmergencyBuff.o BI_BlacsErr.o BI_BlacsWarn.o BI_BlacsAbort.o \
           BI_BuffIsFree.o BI_imvcopy.o BI_smvcopy.o BI_dmvcopy.o \
           BI_ivmcopy.o BI_svmcopy.o BI_dvmcopy.o \
           BI_Pack.o BI_Unpack.o BI_GetMpiGeType.o BI_GetMpiTrType.o \
           BI_ivvsum.o BI_svvsum.o BI_dvvsum.o BI_cvvsum.o BI_zvvsum.o \
           BI_ivvamx.o BI_svvamx.o BI_dvvamx.o BI_cvvamx.o BI_zvvamx.o \
           BI_ivvamx2.o BI_svvamx2.o BI_dvvamx2.o BI_cvvamx2.o BI_zvvamx2.o \
           BI_ivvamn.o BI_svvamn.o BI_dvvamn.o BI_cvvamn.o BI_zvvamn.o \
           BI_ivvamn2.o BI_svvamn2.o BI_dvvamn2.o BI_cvvamn2.o BI_zvvamn2.o \
           BI_iMPI_amx.o BI_sMPI_amx.o BI_dMPI_amx.o BI_cMPI_amx.o \
           BI_zMPI_amx.o BI_iMPI_amx2.o BI_sMPI_amx2.o BI_dMPI_amx2.o \
           BI_cMPI_amx2.o BI_zMPI_amx2.o BI_iMPI_amn.o BI_sMPI_amn.o \
           BI_dMPI_amn.o BI_cMPI_amn.o BI_zMPI_amn.o BI_iMPI_amn2.o \
           BI_sMPI_amn2.o BI_dMPI_amn2.o BI_cMPI_amn2.o BI_zMPI_amn2.o \
           BI_cMPI_sum.o BI_zMPI_sum.o BI_ContxtNum.o BI_GlobalVars.o \
           BI_TransUserComm.o

lib : all

$(internal): Bdef.h Bconfig.h

$(Fintobj): Bdef.h Bconfig.h

#  ---------------------------------------
#  Make both C and fortran interface BLACS
#  ---------------------------------------
all : $(internal) $(Fintobj) $(Cintobj)
	$(ARCH) $(ARCHFLAGS) ../../$(SCALAPACKLIB) $(internal) $(Fintobj) $(Cintobj)
	$(RANLIB) ../../$(SCALAPACKLIB)

#  -----------------------
#  Delete the object files
#  -----------------------
clean :
	rm -f $(Cintobj) $(Fintobj) $(internal)

#  ------------------------------------------------------------------------
#  We move C .o files to .oo so that we can use the portable suffix rule for
#  compilation, and still have them coexist with the fortran interface
#  .o files.
#  ------------------------------------------------------------------------

.SUFFIXES: .o .oo
.c.oo:
	$(CC) -o C$*.o -c $(CDEFS) $(CCFLAGS) -DCallFromC $<
	mv C$*.o $*.oo
.c.o:
	$(CC) -c $(CDEFS) $(CCFLAGS) $<

