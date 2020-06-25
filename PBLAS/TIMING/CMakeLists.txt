file(COPY ../SRC/PTOOLS/PB_Cwarn.c DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
file(COPY ../SRC/PTOOLS/PB_Cabort.c DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})

set (PblasErrorHandler  PB_Cwarn.c PB_Cabort.c)
set (pbtcom pblastim.f ${PblasErrorHandler})
set (spbtcom psblastim.f ${pbtcom})
set (dpbtcom pdblastim.f ${pbtcom})
set (cpbtcom pcblastim.f ${pbtcom})
set (zpbtcom pzblastim.f ${pbtcom})
   
set_property(
   SOURCE ${PblasErrorHandler}
   APPEND PROPERTY COMPILE_DEFINITIONS TestingPblas 
   )

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${SCALAPACK_BINARY_DIR}/PBLAS/TIMING)

file(COPY PCBLAS1TIM.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY PCBLAS2TIM.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY PCBLAS3TIM.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY PDBLAS1TIM.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY PDBLAS2TIM.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY PDBLAS3TIM.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY PSBLAS1TIM.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY PSBLAS2TIM.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY PSBLAS3TIM.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY PZBLAS1TIM.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY PZBLAS2TIM.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY PZBLAS3TIM.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})


add_executable(spb1tim psblas1tim.f ${spbtcom})
add_executable(dpb1tim pdblas1tim.f ${dpbtcom})
add_executable(cpb1tim pcblas1tim.f ${cpbtcom})
add_executable(zpb1tim pzblas1tim.f ${zpbtcom})

add_executable(spb2tim psblas2tim.f ${spbtcom})
add_executable(dpb2tim pdblas2tim.f ${dpbtcom})
add_executable(cpb2tim pcblas2tim.f ${cpbtcom})
add_executable(zpb2tim pzblas2tim.f ${zpbtcom})

add_executable(spb3tim psblas3tim.f ${spbtcom})
add_executable(dpb3tim pdblas3tim.f ${dpbtcom})
add_executable(cpb3tim pcblas3tim.f ${cpbtcom})
add_executable(zpb3tim pzblas3tim.f ${zpbtcom})

target_link_libraries(spb1tim scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(dpb1tim scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(cpb1tim scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(zpb1tim scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

target_link_libraries(spb2tim scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(dpb2tim scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(cpb2tim scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(zpb2tim scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

target_link_libraries(spb3tim scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(dpb3tim scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(cpb3tim scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(zpb3tim scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

add_test(spb1tim ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./spb1tim)
add_test(dpb1tim ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./dpb1tim)
add_test(cpb1tim ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./cpb1tim)
add_test(zpb1tim ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./zpb1tim)

add_test(spb2tim ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./spb2tim)
add_test(dpb2tim ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./dpb2tim)
add_test(cpb2tim ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./cpb2tim)
add_test(zpb2tim ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./zpb2tim)

add_test(spb3tim ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./spb3tim)
add_test(dpb3tim ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./dpb3tim)
add_test(cpb3tim ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./cpb3tim)
add_test(zpb3tim ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./zpb3tim)

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=legacy" )  # local to this directory
endif()
