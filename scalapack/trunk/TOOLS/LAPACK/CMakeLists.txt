set (ALLAUX icopy.f)

set (SCATGEN slatm1.f slaran.f slarnd.f)

set (SMATGEN slatms.f slagge.f slagsy.f slarot.f)

set (CMATGEN clarnv.f clatm1.f clatms.f clagge.f claghe.f clagsy.f clarot.f clarnd.f)

set (DZATGEN dlatm1.f dlaran.f dlarnd.f)

set (DMATGEN dlatms.f dlagge.f dlagsy.f dlarot.f)

set (ZMATGEN zlarnv.f zlatm1.f zlatms.f zlagge.f zlaghe.f zlagsy.f zlarot.f zlarnd.f)

set (extra_lapack 
     ${ALLAUX} ${SCATGEN} ${SMATGEN} ${CMATGEN} ${DZATGEN} ${DMATGEN} ${ZMATGEN})
