set(FTestObj
   blacstest.f btprim.f tools.f)

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
    set_source_files_properties(blacstest.f PROPERTIES COMPILE_FLAGS "-std=legacy")
endif()

add_executable(xFbtest ${FTestObj})
target_link_libraries(xFbtest scalapack)

set(CTestObj
   Cbt.c)

set_property(
   SOURCE Cbt.c
   APPEND PROPERTY COMPILE_DEFINITIONS BTCINTFACE
   )

add_executable(xCbtest ${CTestObj} ${FTestObj})
target_link_libraries(xCbtest scalapack)

file(COPY bsbr.dat DESTINATION ${SCALAPACK_BINARY_DIR}/BLACS/TESTING)
file(COPY bt.dat DESTINATION ${SCALAPACK_BINARY_DIR}/BLACS/TESTING)
file(COPY comb.dat DESTINATION ${SCALAPACK_BINARY_DIR}/BLACS/TESTING)
file(COPY sdrv.dat DESTINATION ${SCALAPACK_BINARY_DIR}/BLACS/TESTING)

# We could run the BLACS TESTING the following way
# But BLACS TESTING are TESTING anormal exit so even if they pass,
# CTest will determine they fail
#add_test(xFbtest0 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xFbtest)
#add_test(xCbtest0 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./xCbtest)

add_test(xCbtest
     ${CMAKE_COMMAND}
    -DMPIEXEC=${MPIEXEC}
    -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
    -DTEST_PROG=xCbtest
    -DOUTPUTDIR=${SCALAPACK_BINARY_DIR}/BLACS/TESTING
    -DRUNTIMEDIR=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    -DSOURCEDIR=${CMAKE_CURRENT_SOURCE_DIR}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/runtest.cmake
    )

add_test(xFbtest
     ${CMAKE_COMMAND}
    -DMPIEXEC=${MPIEXEC}
    -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
    -DTEST_PROG=xFbtest
    -DOUTPUTDIR=${SCALAPACK_BINARY_DIR}/BLACS/TESTING
    -DRUNTIMEDIR=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    -DSOURCEDIR=${CMAKE_CURRENT_SOURCE_DIR}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/runtest.cmake
    )
