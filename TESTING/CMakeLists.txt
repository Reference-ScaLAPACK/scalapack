add_subdirectory(LIN)
add_subdirectory(EIG)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${SCALAPACK_BINARY_DIR}/TESTING)

file(COPY BLLT.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY BLU.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY BRD.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY EVC.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY HRD.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY INV.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY LLT.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY LS.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY LU.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY NEP.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY QR.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY SEP.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY SEPR.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY SVD.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY TRD.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

add_test(xslu ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xslu)
add_test(xdlu ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xdlu)
add_test(xclu ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xclu)
add_test(xzlu ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xzlu)

add_test(xsdblu ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xsdblu)
add_test(xddblu ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xddblu)
add_test(xcdblu ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xcdblu)
add_test(xzdblu ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xzdblu)

add_test(xsdtlu ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xsdtlu)
add_test(xddtlu ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xddtlu)
add_test(xcdtlu ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xcdtlu)
add_test(xzdtlu ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xzdtlu)

add_test(xsgblu ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xsgblu)
add_test(xdgblu ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xdgblu)
add_test(xcgblu ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xcgblu)
add_test(xzgblu ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xzgblu)

add_test(xsllt ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xsllt)
add_test(xdllt ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xdllt)
add_test(xcllt ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xcllt)
add_test(xzllt ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xzllt)

add_test(xspbllt ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xspbllt)
add_test(xdpbllt ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xdpbllt)
add_test(xcpbllt ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xcpbllt)
add_test(xzpbllt ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xzpbllt)

add_test(xsptllt ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xsptllt)
add_test(xdptllt ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xdptllt)
add_test(xcptllt ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xcptllt)
add_test(xzptllt ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xzptllt)

add_test(xsinv ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xsinv)
add_test(xdinv ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xdinv)
add_test(xcinv ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xcinv)
add_test(xzinv ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xzinv)

add_test(xsqr ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xsqr)
add_test(xdqr ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xdqr)
add_test(xcqr ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xcqr)
add_test(xzqr ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xzqr)

add_test(xsbrd ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xsbrd)
add_test(xdbrd ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xdbrd)
add_test(xcbrd ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xcbrd)
add_test(xzbrd ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xzbrd)

add_test(xshrd ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xshrd)
add_test(xdhrd ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xdhrd)
add_test(xchrd ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xchrd)
add_test(xzhrd ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xzhrd)

add_test(xstrd ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xstrd)
add_test(xdtrd ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xdtrd)
add_test(xctrd ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xctrd)
add_test(xztrd ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xztrd)

add_test(xssvd ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xssvd)
add_test(xdsvd ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xdsvd)

add_test(xssep ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xssep)
add_test(xdsep ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xdsep)
add_test(xcsep ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xcsep)
add_test(xzsep ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xzsep)

add_test(xsgsep ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xsgsep)
add_test(xdgsep ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xdgsep)
add_test(xcgsep ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xcgsep)
add_test(xzgsep ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xzgsep)

add_test(xsnep ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xsnep)
add_test(xdnep ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xdnep)
add_test(xcnep ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xcnep)
add_test(xznep ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xznep)

add_test(xcevc ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xcevc)
add_test(xzevc ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xzevc)

add_test(xssyevr ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xssyevr)
add_test(xdsyevr ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xdsyevr)
add_test(xcheevr ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xcheevr)
add_test(xzheevr ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xzheevr)

add_test(xshseqr ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xshseqr)
add_test(xdhseqr ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xdhseqr)

