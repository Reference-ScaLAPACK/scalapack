set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${SCALAPACK_BINARY_DIR}/TESTING)

set (smatgen psmatgen.f pmatgeninc.f)
set (dmatgen pdmatgen.f pmatgeninc.f)
set (cmatgen pcmatgen.f pmatgeninc.f)
set (zmatgen pzmatgen.f pmatgeninc.f)

add_executable(xsbrd psbrddriver.f psbrdinfo.f psgebdrv.f pslafchk.f ${smatgen})
add_executable(xdbrd pdbrddriver.f pdbrdinfo.f pdgebdrv.f pdlafchk.f ${dmatgen})
add_executable(xcbrd pcbrddriver.f pcbrdinfo.f pcgebdrv.f pclafchk.f ${cmatgen})
add_executable(xzbrd pzbrddriver.f pzbrdinfo.f pzgebdrv.f pzlafchk.f ${zmatgen})

add_executable(xshrd pshrddriver.f pshrdinfo.f psgehdrv.f pslafchk.f ${smatgen})
add_executable(xdhrd pdhrddriver.f pdhrdinfo.f pdgehdrv.f pdlafchk.f ${dmatgen})
add_executable(xchrd pchrddriver.f pchrdinfo.f pcgehdrv.f pclafchk.f ${cmatgen})
add_executable(xzhrd pzhrddriver.f pzhrdinfo.f pzgehdrv.f pzlafchk.f ${zmatgen})

add_executable(xstrd pstrddriver.f psttrdtester.f pslatran.f pstrdinfo.f pssytdrv.f pslafchk.f xpjlaenv.f ${smatgen})
add_executable(xdtrd pdtrddriver.f pdttrdtester.f pdlatran.f pdtrdinfo.f pdsytdrv.f pdlafchk.f xpjlaenv.f ${dmatgen})
add_executable(xctrd pctrddriver.f pcttrdtester.f pclatran.f pctrdinfo.f pchetdrv.f pclafchk.f xpjlaenv.f ${cmatgen})
add_executable(xztrd pztrddriver.f pzttrdtester.f pzlatran.f pztrdinfo.f pzhetdrv.f pzlafchk.f xpjlaenv.f ${zmatgen})

add_executable(xssvd pssvddriver.f pslagge.f pssvdchk.f pssvdcmp.f pssvdtst.f ${smatgen})
add_executable(xdsvd pdsvddriver.f pdlagge.f pdsvdchk.f pdsvdcmp.f pdsvdtst.f ${dmatgen})

add_executable(xssep psseptst.f pssepsubtst.f pssepchk.f pssepqtq.f pslagsy.f pslatms.f pslasizesep.f pslasizesyevx.f pssepdriver.f pssepreq.f pssepinfo.f pslasizesyev.f pssqpsubtst.f pslasizesqp.f pssdpsubtst.f ${smatgen})
add_executable(xdsep pdseptst.f pdsepsubtst.f pdsepchk.f pdsepqtq.f pdlagsy.f pdlatms.f pdlasizesep.f pdlasizesyevx.f pdsepdriver.f pdsepreq.f pdsepinfo.f pdlasizesyev.f pdsqpsubtst.f pdlasizesqp.f pdsdpsubtst.f ${dmatgen})
add_executable(xcsep pcseptst.f pcsepsubtst.f pcsepchk.f pcsepqtq.f pclagsy.f pclatms.f pclasizesep.f pclasizeheevx.f pcsepdriver.f pcsepreq.f pssepinfo.f pcsdpsubtst.f ${cmatgen})
add_executable(xzsep pzseptst.f pzsepsubtst.f pzsepchk.f pzsepqtq.f pzlagsy.f pzlatms.f pzlasizesep.f pzlasizeheevx.f pzsepdriver.f pzsepreq.f pdsepinfo.f pzsdpsubtst.f ${zmatgen})

add_executable(xsgsep psgseptst.f psgsepsubtst.f psgsepchk.f pslagsy.f pslatms.f pslasizesyevx.f pslasizegsep.f pslasizesep.f psgsepdriver.f psgsepreq.f pssepinfo.f ${smatgen})
add_executable(xdgsep pdgseptst.f pdgsepsubtst.f pdgsepchk.f pdlagsy.f pdlatms.f pdlasizesyevx.f pdlasizegsep.f pdlasizesep.f pdgsepdriver.f pdgsepreq.f pdsepinfo.f ${dmatgen})
add_executable(xcgsep pcgseptst.f pcgsepsubtst.f pcgsepchk.f pclagsy.f pclatms.f pclasizegsep.f pclasizeheevx.f pclasizesep.f pcgsepdriver.f pcgsepreq.f pssepinfo.f ${cmatgen})
add_executable(xzgsep pzgseptst.f pzgsepsubtst.f pzgsepchk.f pzlagsy.f pzlatms.f pzlasizegsep.f pzlasizeheevx.f pzlasizesep.f pzgsepdriver.f pzgsepreq.f pdsepinfo.f ${zmatgen})

add_executable(xsnep psnepdriver.f psnepinfo.f psnepfchk.f ${smatgen})
add_executable(xdnep pdnepdriver.f pdnepinfo.f pdnepfchk.f ${dmatgen})
add_executable(xcnep pcnepdriver.f pcnepinfo.f pcnepfchk.f ${cmatgen})
add_executable(xznep pznepdriver.f pznepinfo.f pznepfchk.f ${zmatgen})

add_executable(xcevc pcevcdriver.f pcevcinfo.f pcget22.f ${cmatgen})
add_executable(xzevc pzevcdriver.f pzevcinfo.f pzget22.f ${zmatgen})

add_executable(xssyevr pslasizesepr.f pslasizesyevr.f psseprdriver.f psseprreq.f psseprsubtst.f
pssepchk.f pssepqtq.f pslatms.f psseprtst.f pssepinfo.f pslagsy.f pslasizesep.f ${smatgen})
add_executable(xdsyevr  pdlasizesepr.f pdlasizesyevr.f pdseprdriver.f pdseprreq.f pdseprsubtst.f
pdsepchk.f pdsepqtq.f pdlatms.f pdseprtst.f pdsepinfo.f pdlagsy.f pdlasizesep.f ${dmatgen}) 
add_executable(xcheevr pclasizesepr.f pclasizeheevr.f pcseprdriver.f pcseprreq.f pcseprsubtst.f
pcsepchk.f pcsepqtq.f pclatms.f pcseprtst.f pssepinfo.f pclagsy.f pclasizesep.f ${cmatgen})
add_executable(xzheevr  pzlasizesepr.f pzlasizeheevr.f pzseprdriver.f pzseprreq.f pzseprsubtst.f
pzsepchk.f pzsepqtq.f pzlatms.f pzseprtst.f pdsepinfo.f pzlagsy.f pzlasizesep.f ${zmatgen})

add_executable(xshseqr pshseqrdriver.f psmatgen2.f ${cmatgen})
add_executable(xdhseqr pdhseqrdriver.f pdmatgen2.f ${cmatgen})

target_link_libraries(xsbrd scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(xdbrd scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(xcbrd scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(xzbrd scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

target_link_libraries(xshrd scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(xdhrd scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(xchrd scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(xzhrd scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

target_link_libraries(xstrd scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(xdtrd scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(xctrd scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(xztrd scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

target_link_libraries(xssvd scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(xdsvd scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

target_link_libraries(xssep scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(xdsep scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(xcsep scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(xzsep scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

target_link_libraries(xsgsep scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(xdgsep scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(xcgsep scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(xzgsep scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

target_link_libraries(xsnep scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(xdnep scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(xcnep scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(xznep scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

target_link_libraries(xcevc scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(xzevc scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

target_link_libraries(xssyevr scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(xdsyevr scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(xcheevr scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(xzheevr scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

target_link_libraries(xshseqr scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(xdhseqr scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

