#  --------------------------
#  The communication routines
#  --------------------------
set(comm
igesd2d_.c itrsd2d_.c igerv2d_.c itrrv2d_.c igebs2d_.c igebr2d_.c itrbs2d_.c itrbr2d_.c igsum2d_.c igamx2d_.c igamn2d_.c
dgsum2d_.c dgamx2d_.c dgamn2d_.c
)
if(SCALAPACK_BUILD_SINGLE)
  list(APPEND comm sgesd2d_.c strsd2d_.c sgerv2d_.c strrv2d_.c sgebs2d_.c sgebr2d_.c strbs2d_.c strbr2d_.c sgsum2d_.c sgamx2d_.c sgamn2d_.c)
endif()
if(SCALAPACK_BUILD_DOUBLE)
  list(APPEND comm dgesd2d_.c dtrsd2d_.c dgerv2d_.c dtrrv2d_.c dgebs2d_.c dgebr2d_.c dtrbs2d_.c dtrbr2d_.c)
endif()
if(SCALAPACK_BUILD_COMPLEX)
  list(APPEND comm cgesd2d_.c ctrsd2d_.c cgerv2d_.c ctrrv2d_.c cgebs2d_.c cgebr2d_.c ctrbs2d_.c ctrbr2d_.c cgsum2d_.c cgamx2d_.c cgamn2d_.c)
endif()
if(SCALAPACK_BUILD_COMPLEX16)
  list(APPEND comm zgesd2d_.c ztrsd2d_.c zgerv2d_.c ztrrv2d_.c zgebs2d_.c zgebr2d_.c ztrbs2d_.c ztrbr2d_.c zgsum2d_.c zgamx2d_.c zgamn2d_.c)
endif()
#  --------------------
#  The support routines
#  --------------------
set (supp
       blacs_setup_.c blacs_set_.c blacs_get_.c
       blacs_abort_.c blacs_exit_.c blacs_pnum_.c blacs_pcoord_.c
       ksendid_.c krecvid_.c kbsid_.c kbrid_.c
       dcputime00_.c dwalltime00_.c blacs_pinfo_.c
       blacs_init_.c blacs_map_.c blacs_free_.c blacs_grid_.c blacs_info_.c
       blacs_barr_.c sys2blacs_.c blacs2sys_.c free_handle_.c)

#  ---------------------
#  The internal routines
#  ---------------------
set (internal
           BI_HypBS.c BI_HypBR.c BI_IdringBS.c BI_IdringBR.c
           BI_MpathBS.c BI_MpathBR.c BI_SringBS.c BI_SringBR.c
           BI_TreeBS.c BI_TreeBR.c
           BI_Ssend.c BI_Rsend.c BI_Srecv.c BI_Asend.c BI_Arecv.c
           BI_TreeComb.c BI_BeComb.c BI_MringComb.c
           BI_ArgCheck.c BI_TransDist.c BI_GetBuff.c BI_UpdateBuffs.c
           BI_EmergencyBuff.c BI_BlacsErr.c BI_BlacsWarn.c BI_BlacsAbort.c
           BI_BuffIsFree.c BI_imvcopy.c BI_smvcopy.c BI_dmvcopy.c
           BI_ivmcopy.c BI_svmcopy.c BI_dvmcopy.c
           BI_Pack.c BI_Unpack.c BI_GetMpiGeType.c BI_GetMpiTrType.c
           BI_ivvsum.c BI_ivvamx.c BI_ivvamx2.c BI_ivvamn.c BI_ivvamn2.c
           BI_iMPI_amx.c BI_iMPI_amx2.c BI_iMPI_amn.c BI_iMPI_amn2.c
           BI_ContxtNum.c BI_GlobalVars.c BI_TransUserComm.c
           BI_dvvsum.c BI_dvvamx.c BI_dvvamx2.c BI_dvvamn.c BI_dvvamn2.c
           BI_dMPI_amx.c BI_dMPI_amx2.c BI_dMPI_amn.c BI_dMPI_amn2.c
           )

if(SCALAPACK_BUILD_SINGLE)
  list(APPEND internal
       BI_svvsum.c BI_svvamx.c BI_svvamx2.c BI_svvamn.c BI_svvamn2.c
       BI_sMPI_amx.c BI_sMPI_amx2.c BI_sMPI_amn.c BI_sMPI_amn2.c
       BI_smvcopy.c BI_svmcopy.c
       )
endif()
if(SCALAPACK_BUILD_COMPLEX)
  list(APPEND internal
       BI_cvvsum.c BI_cvvamx.c BI_cvvamx2.c BI_cvvamn.c BI_cvvamn2.c
       BI_cMPI_amx.c BI_cMPI_amx2.c BI_cMPI_amn.c BI_cMPI_amn2.c BI_cMPI_sum.c
       )
endif()
if(SCALAPACK_BUILD_COMPLEX16)
  list(APPEND internal
       BI_zvvsum.c BI_zvvamx.c BI_zvvamx2.c BI_zvvamn.c BI_zvvamn2.c
       BI_zMPI_amx.c BI_zMPI_amx2.c BI_zMPI_amn.c BI_zMPI_amn2.c BI_zMPI_sum.c
       )
endif()
# Compile the BLACS entry points twice to provide both Fortran and C APIs.
set(_blacs_interface_sources ${comm} ${supp})
add_library(scalapack-blacs-fortran OBJECT ${_blacs_interface_sources})
add_library(scalapack-blacs-c OBJECT ${_blacs_interface_sources})
target_compile_definitions(scalapack-blacs-c PRIVATE CallFromC)

foreach(_blacs_target IN ITEMS scalapack-blacs-fortran scalapack-blacs-c)
  set_target_properties(${_blacs_target} PROPERTIES POSITION_INDEPENDENT_CODE ON)
  target_link_libraries(${_blacs_target} PRIVATE MPI::MPI_C)
  target_sources(scalapack PRIVATE $<TARGET_OBJECTS:${_blacs_target}>)
endforeach()

target_sources(scalapack PRIVATE ${internal})
