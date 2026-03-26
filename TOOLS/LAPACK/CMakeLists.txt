set (ALLAUX icopy.f)

set (SCATGEN slatm1.f slaran.f slarnd.f)

set (SMATGEN slatms.f slagge.f slagsy.f slarot.f)

set (CMATGEN clarnv.f clatm1.f clatms.f clagge.f claghe.f clagsy.f clarot.f clarnd.f)

set (DZATGEN dlatm1.f dlaran.f dlarnd.f)

set (DMATGEN dlatms.f dlagge.f dlagsy.f dlarot.f)

set (ZMATGEN zlarnv.f zlatm1.f zlatms.f zlagge.f zlaghe.f zlagsy.f zlarot.f zlarnd.f)

set(extra_lapack ${ALLAUX})
if(BUILD_SINGLE)
  list(APPEND extra_lapack ${SMATGEN})
endif()
if(BUILD_SINGLE OR BUILD_COMPLEX)
  list(APPEND extra_lapack ${SCATGEN})
endif()
if(BUILD_DOUBLE)
  list(APPEND extra_lapack ${DMATGEN})
endif()
if(BUILD_DOUBLE OR BUILD_COMPLEX16)
  list(APPEND extra_lapack ${DZATGEN})
endif()
if(BUILD_COMPLEX)
  list(APPEND extra_lapack ${CMATGEN})
endif()
if(BUILD_COMPLEX16)
  list(APPEND extra_lapack  ${ZMATGEN})
endif()
